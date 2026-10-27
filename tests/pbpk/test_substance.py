"""The endogenous baseline of a substance (`Substance.baseline`).

A substance with a baseline starts at its physiological concentration: every
plasma species at the baseline, every other species of the blood and the
tissue at the equilibrium of its exchange with plasma, the lumen, the tubule
and the sinks at 0. Without production and elimination the body is then
stationary, in every form of the organ the library builds.
"""

import dataclasses
from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING, Any

import libsbml
import numpy as np
import pytest

from pkpdsim.pbpk import (
    MM,
    Binding,
    Erythrocytes,
    Formula,
    Inhibition,
    Kp,
    KpMethod,
    Metabolism,
    Mode,
    Observer,
    Organ,
    Part,
    PBPKDefinition,
    Pool,
    Route,
    Substance,
    U,
    Value,
    Variants,
    build,
    compose,
    ids,
)
from pkpdsim.pbpk.variants.intestine import SingleLumen
from pkpdsim.pbpk.variants.kidney_nephron import NephronKidney
from pkpdsim.pbpk.variants.liver_zonated import ZonatedLiver

from .conftest import (
    assert_reference_sbml,
    cardiac_pd,
    individual,
    load,
    reference_definition,
)

if TYPE_CHECKING:
    import roadrunner

S = "end"
B = 0.05
DAY = 24 * 60.0


def _endogenous(**fields: Any) -> Substance:  # noqa: ANN401
    """Get an endogenous substance: baseline 0.05 mM, Kp 2, no elimination."""
    values: dict[str, Any] = {
        "sid": S,
        "name": "endogenous substance",
        "mr": 150.0,
        "binding": Binding(fu=0.5, protein="albumin"),
        "distribution": Kp(2.0, ftissue=1.0),
        "baseline": B,
    }
    values.update(fields)
    return Substance(**values)


def _definition(sid: str, s: Substance, **fields: Any) -> PBPKDefinition:  # noqa: ANN401
    """Get the definition of one substance."""
    return PBPKDefinition(sid=sid, substances=[s], **fields)


def _species(r: "roadrunner.RoadRunner") -> list[str]:
    """The floating species of a loaded model."""
    return list(r.model.getFloatingSpeciesIds())


def _drift(r: "roadrunner.RoadRunner") -> dict[str, float]:
    """Simulate 24 hr, the largest relative change of every species (absolute at 0)."""
    sids = _species(r)
    start = np.array([r[f"[{sid}]"] for sid in sids])
    result = r.simulate(0, DAY, 49, selections=[f"[{sid}]" for sid in sids])
    values = np.asarray(result)
    change = np.max(np.abs(values - start), axis=0)
    scale = np.where(start > 0, start, 1.0)
    return dict(zip(sids, change / scale, strict=True))


def _assert_stationary(r: "roadrunner.RoadRunner") -> None:
    """Every species stays at its start over 24 hr, rtol 1e-6 (atol 1e-12 at 0)."""
    sids = _species(r)
    zero = {sid for sid in sids if r[f"[{sid}]"] == 0.0}
    drift = _drift(r)
    moved = {sid: d for sid, d in drift.items() if d > (1e-12 if sid in zero else 1e-6)}
    assert not moved, moved


def _local(sid: str) -> str:
    """The local id of a flat id, `MU__Cti_end` -> `Cti_end`."""
    return sid.rsplit("__", 1)[-1]


_ORGAN = {"GU": Organ.SMALL_INTESTINE.value}
"""The organ of the Kp of a submodel whose id is no organ token."""


def _expected_starts(r: "roadrunner.RoadRunner") -> dict[str, float]:
    """The start of every species from the kind of its part, by the coefficients of the model.

    Plasma (organs and pools) the baseline, blood cells `Kbc*B`, merged
    tissue `Kp*B` (of the organ), interstitium `K_int_pls*B`, cells
    `K_cell_pls*B`, the flow-limited organ `K_or*B`, the mucosa `K_mu*B`,
    everything else (lumen, tubule, sinks) 0.
    """
    b = r[ids.baseline(S)]
    kbc = ids.kbc(S)
    top = set(r.model.getGlobalParameterIds())
    expected: dict[str, float] = {}
    for sid in _species(r):
        organ, _, local = sid.rpartition("__")
        prefix = f"{organ}__" if organ else ""
        part = local.split("_", 1)[0]
        if not sid.endswith(f"_{S}"):
            expected[sid] = 0.0  # another substance, without baseline
        elif sid in {ids.pool_species(p, S) for p in Pool}:
            expected[sid] = b
        elif sid in {ids.pool_species(p, S, Part.BLOOD_CELLS) for p in Pool}:
            expected[sid] = r[kbc] * b
        elif part == "Cpl":
            expected[sid] = b
        elif part == "Cbc":
            expected[sid] = r[kbc] * b
        elif part == "Cti":
            # the Kp of the organ (override, method) or of the substance; the
            # wall of the intestine (GU) has the Kp of the small intestine
            own = ids.kp_organ(Organ(_ORGAN.get(organ, organ.lower())), S)
            expected[sid] = r[own if own in top else ids.kp(S)] * b
        elif part == "Cis":
            expected[sid] = r[f"{prefix}{ids.k_int_pls(S)}"] * b
        elif part == "Cic":
            expected[sid] = r[f"{prefix}{ids.k_cell_pls(S)}"] * b
        elif part == "Cor":
            expected[sid] = r[f"{prefix}{ids.k_organ(S)}"] * b
        elif part == "Cmu":
            expected[sid] = r[f"{prefix}{ids.mucosa_coefficient(S)}"] * b
        else:
            expected[sid] = 0.0
    return expected


def _assert_starts(r: "roadrunner.RoadRunner") -> None:
    """Every species starts at the equilibrium of its part with the baseline."""
    for sid, value in _expected_starts(r).items():
        assert r[f"[{sid}]"] == pytest.approx(value, rel=1e-12, abs=1e-15), sid


def test_a_baseline_without_turnover_is_stationary(tmp_path: Path) -> None:
    """Baseline 0.05 mM, Kp 2, no elimination, no dose: every species stays at its start.

    The plasma species (pools, organs) start at 0.05 mM, the tissue species
    at `Kp*B` = 0.1 mM, the arm vein at the baseline.
    """
    r = load(build(_definition("baseline", _endogenous()), tmp_path))
    assert r[ids.baseline(S)] == B
    sids = _species(r)
    plasma = [sid for sid in sids if _local(sid).startswith(("Cpl_", "Cve_", "Car_"))]
    tissue = [sid for sid in sids if _local(sid).startswith("Cti_")]
    assert len(plasma) > 10
    assert len(tissue) > 10
    for sid in plasma:
        assert r[f"[{sid}]"] == pytest.approx(B, rel=1e-12), sid
    for sid in tissue:
        assert r[f"[{sid}]"] == pytest.approx(2 * B, rel=1e-12), sid
    assert r[f"[{ids.pool_species(Pool.VENOUS, S)}]"] == pytest.approx(B, rel=1e-12)
    assert r[ids.arm(S)] == pytest.approx(B, rel=1e-12)
    _assert_starts(r)
    _assert_stationary(r)
    assert r[ids.arm(S)] == pytest.approx(B, rel=1e-6)


_INHIBITION = Inhibition("CYP3A4", "competitive", ki=0.01)
"""An inhibition of the intestinal CYP3A4, which gives the substance mucosa species."""


def _with_mucosa(sid: str, s: Substance) -> PBPKDefinition:
    """A definition in which the endogenous substance has mucosa species, without elimination.

    A drug (no baseline, no dose) is metabolized by CYP3A4 of the small
    intestine, which the endogenous substance inhibits: the inhibitor of an
    intestinal protein has mucosa species; the drug stays at 0.
    """
    drug = Substance("drg", name="drug", mr=300.0, distribution=Kp(1.0, ftissue=1.0))
    gut = Metabolism(
        sid="DRGGU",
        organ=Organ.SMALL_INTESTINE,
        substrate="drg",
        products={},
        kinetics=MM(vmax=1e-3, km=0.01),
        protein="CYP3A4",
    )
    return PBPKDefinition(
        sid=sid,
        substances=[s, drug],
        chemistry=[gut],
        routes=[Route.oral(S), Route.iv(S)],
    )


_MODES: dict[str, Callable[[], PBPKDefinition]] = {
    "perfusion, erythrocytes bp": lambda: _definition(
        "bl_perfusion",
        _endogenous(
            distribution=Kp(2.0, mode=Mode.PERFUSION),
            erythrocytes=Erythrocytes("equilibrium", bp=1.3),
        ),
    ),
    "tissue, erythrocytes kbc, overrides, organs, ps": lambda: _definition(
        "bl_tissue",
        _endogenous(
            distribution=Kp(
                2.0,
                overrides={Organ.ADIPOSE: 5.0},
                organs=frozenset(
                    {Organ.ADIPOSE, Organ.LIVER, Organ.MUSCLE, Organ.KIDNEY}
                ),
                ftissue=1.0,
                ps={Organ.LIVER: 40.0},
            ),
            erythrocytes=Erythrocytes("equilibrium", kbc=1.5),
        ),
    ),
    "permeability, kinetic erythrocytes": lambda: _definition(
        "bl_permeability",
        _endogenous(
            distribution=Kp(3.0, mode=Mode.PERMEABILITY, ftissue=1.0),
            erythrocytes=Erythrocytes("kinetic", kbc=2.0, k_exchange=1.0),
        ),
    ),
    "kp method": lambda: _definition(
        "bl_method",
        _endogenous(
            distribution=KpMethod(
                "rodgers_rowland",
                logp=1.0,
                pka=(),
                compound="neutral",
                ftissue=1.0,
            ),
        ),
    ),
    "vascular only": lambda: _definition("bl_vascular", _endogenous(distribution=None)),
    "zonated liver, nephron kidney, single lumen, oral route": lambda: _definition(
        "bl_variants",
        _endogenous(distribution=Kp(2.0, mode=Mode.PERMEABILITY, ftissue=1.0)),
        routes=[Route.oral(S)],
        variants=Variants(
            organs={
                Organ.LIVER: ZonatedLiver(n_zones=3),
                Organ.KIDNEY: NephronKidney(),
                Organ.SMALL_INTESTINE: SingleLumen(),
            }
        ),
    ),
    "segmented intestine with mucosa, flow-limited": lambda: _with_mucosa(
        "bl_mucosa",
        _endogenous(
            distribution=Kp(2.0, mode=Mode.PERFUSION),
            erythrocytes=Erythrocytes("equilibrium", kbc=1.5),
            inhibits=(_INHIBITION,),
        ),
    ),
    "segmented intestine with mucosa, kinetic erythrocytes": lambda: _with_mucosa(
        "bl_mucosa_bc",
        _endogenous(
            erythrocytes=Erythrocytes("kinetic", kbc=2.0, k_exchange=1.0),
            inhibits=(_INHIBITION,),
        ),
    ),
}


_PARTS = {
    "perfusion, erythrocytes bp": {"Cor"},
    "tissue, erythrocytes kbc, overrides, organs, ps": {"Cpl", "Cti"},
    "permeability, kinetic erythrocytes": {"Cpl", "Cbc", "Cis", "Cic"},
    "kp method": {"Cpl", "Cti"},
    "vascular only": {"Cpl"},
    "zonated liver, nephron kidney, single lumen, oral route": {"Cpl", "Cis", "Cic"},
    "segmented intestine with mucosa, flow-limited": {"Cor", "Cmu"},
    "segmented intestine with mucosa, kinetic erythrocytes": {
        "Cpl",
        "Cbc",
        "Cti",
        "Cmu",
    },
}
"""The parts of the organs whose species each case builds for the substance."""


@pytest.mark.parametrize("case", list(_MODES))
def test_a_baseline_is_stationary_in_every_mode(tmp_path: Path, case: str) -> None:
    """Every form of the organ, erythrocyte mode and variant starts at the equilibrium.

    The case builds the species of its parts (`_PARTS`), and every species
    starts at the equilibrium of its part with the baseline.
    """
    r = load(build(_MODES[case](), tmp_path))
    parts = {
        _local(sid).split("_", 1)[0] for sid in _species(r) if sid.endswith(f"_{S}")
    }
    assert _PARTS[case] <= parts
    _assert_starts(r)
    assert r[ids.arm(S)] == pytest.approx(B, rel=1e-12)
    _assert_stationary(r)


def test_the_cases_cover_every_part() -> None:
    """Together the cases build plasma, blood cells, tissue, interstitium, cells, organ, mucosa."""
    assert set(_PARTS) == set(_MODES)
    covered = set().union(*_PARTS.values())
    assert covered == {"Cpl", "Cbc", "Cti", "Cis", "Cic", "Cor", "Cmu"}


def test_a_baseline_follows_a_preinit_impairment(tmp_path: Path) -> None:
    """Cirrhosis, CKD and heart failure before the simulation: the start follows them.

    `fu`, the coefficients and the functional liver change with the
    impairment; the species read the coefficients, so the body stays
    stationary.
    """
    d = _definition(
        "bl_impaired",
        _endogenous(distribution=Kp(3.0, mode=Mode.PERMEABILITY, ftissue=1.0)),
    )
    path = build(d, tmp_path).flat_path
    r = individual(
        path,
        **{ids.CP_SEVERITY: 2.0, ids.F_RENAL_FUNCTION: 0.3, ids.HF_SEVERITY: 2.0},
    )
    healthy = individual(path)
    assert r[ids.fu(S)] > healthy[ids.fu(S)]
    assert r["[LI__Cis_end]"] != pytest.approx(healthy["[LI__Cis_end]"])
    r.integrator.absolute_tolerance = 1e-12
    r.integrator.relative_tolerance = 1e-10
    _assert_starts(r)
    _assert_stationary(r)


def test_a_formula_baseline_follows_a_preinit_change(tmp_path: Path) -> None:
    """A baseline of production over clearance: a change of the clearance moves every start.

    `baseline_<s> = 0.006 mmole_per_min / GFR`; the renal function 0.5 halves
    the GFR before the simulation, so every species starts at twice its
    healthy concentration.
    """
    d = _definition(
        "bl_formula",
        _endogenous(baseline=Formula("0.006 mmole_per_min / GFR", U.mM)),
    )
    path = build(d, tmp_path).flat_path
    healthy = individual(path)
    impaired = individual(path, **{ids.F_RENAL_FUNCTION: 0.5})
    b = healthy[ids.baseline(S)]
    assert b == pytest.approx(0.006 / healthy[ids.GFR], rel=1e-12)
    assert impaired[ids.baseline(S)] == pytest.approx(2 * b, rel=1e-9)
    for sid in _species(healthy):
        start = healthy[f"[{sid}]"]
        assert impaired[f"[{sid}]"] == pytest.approx(2 * start, rel=1e-9), sid
    assert healthy[f"[{ids.pool_species(Pool.VENOUS, S)}]"] == pytest.approx(b)
    assert healthy["[MU__Cti_end]"] == pytest.approx(2 * b)


def test_a_formula_baseline_may_read_an_observer(tmp_path: Path) -> None:
    """Production over a clearance observer of the definition; the renal function is dynamic.

    With `f_renal_function` dynamic the GFR is a rule, the observer of the
    clearance reads it; the baseline reads the observer at the start, the
    species keep their initial assignments.
    """
    clearance = Observer(
        "CL_end", f"{ids.fu(S)} * {ids.GFR}", U.l_per_min, "clearance of end"
    )
    d = _definition(
        "bl_observer",
        _endogenous(baseline=Formula("0.003 mmole_per_min / CL_end", U.mM)),
        observers=[clearance],
        dynamic=(ids.F_RENAL_FUNCTION,),
    )
    result = build(d, tmp_path)
    r = load(result)
    b = 0.003 / (r[ids.fu(S)] * r[ids.GFR])
    assert r[ids.baseline(S)] == pytest.approx(b, rel=1e-12)
    assert r[f"[{ids.pool_species(Pool.VENOUS, S)}]"] == pytest.approx(b)
    model = libsbml.readSBMLFromFile(str(result.flat_path)).getModel()
    species = {sp.getId() for sp in model.getListOfSpecies()}
    assert not [rule for rule in model.getListOfRules() if rule.getId() in species]
    _assert_stationary(r)
    # the start follows the renal function before the simulation
    impaired = individual(result.flat_path, **{ids.F_RENAL_FUNCTION: 0.5})
    assert impaired[ids.baseline(S)] == pytest.approx(
        0.003 / (impaired[ids.fu(S)] * impaired[ids.GFR]), rel=1e-9
    )
    assert impaired[ids.GFR] == pytest.approx(0.5 * r[ids.GFR], rel=1e-9)


def test_a_dynamic_coefficient_keeps_the_species_initial_assignments(
    tmp_path: Path,
) -> None:
    """With `CP_severity` dynamic, fu and the coefficients are rules, the species start once."""
    d = _definition(
        "bl_dynamic",
        _endogenous(distribution=Kp(3.0, mode=Mode.PERMEABILITY, ftissue=1.0)),
        dynamic=(ids.CP_SEVERITY,),
    )
    result = build(d, tmp_path)
    model = libsbml.readSBMLFromFile(str(result.flat_path)).getModel()
    species = {sp.getId() for sp in model.getListOfSpecies()}
    assert model.getRule(f"LI__{ids.k_int_pls(S)}") is not None
    assert not [rule for rule in model.getListOfRules() if rule.getId() in species]
    assert model.getInitialAssignment("LI__Cis_end") is not None
    r = load(result)
    _assert_starts(r)
    _assert_stationary(r)


def test_a_formula_baseline_refuses_what_changes_during_the_simulation() -> None:
    """A baseline which reads a species, directly or through an observer, is refused."""
    pool = ids.pool_species(Pool.VENOUS, S)
    with pytest.raises(ValueError, match=rf"baseline_{S}.*{pool}"):
        compose(_definition("bl_species", _endogenous(baseline=Formula(pool, U.mM))))
    observer = Observer("C_obs", f"2 dimensionless * {pool}", U.mM, "twice the pool")
    with pytest.raises(ValueError, match=rf"baseline_{S}.*C_obs.*{pool}"):
        compose(
            _definition(
                "bl_cycle",
                _endogenous(baseline=Formula("C_obs", U.mM)),
                observers=[observer],
            )
        )
    with pytest.raises(ValueError, match=rf"baseline_{S}.*unknown.*k_prod"):
        compose(
            _definition("bl_unknown", _endogenous(baseline=Formula("k_prod", U.mM)))
        )


def test_a_formula_baseline_refuses_a_pd_output(tmp_path: Path) -> None:
    """A PD output reads the species: a baseline must not read it nor what follows it."""
    pd = cardiac_pd(reads=S)
    for read in (ids.F_CARDIAC_FUNCTION, ids.CO):
        d = _definition(
            "bl_pd",
            _endogenous(baseline=Formula(f"0.006 mmole_per_min / {read}", U.mM)),
            pd=[pd],
        )
        with pytest.raises(ValueError, match=rf"baseline_{S}.*{read} \(a PD output"):
            compose(d)
    observer = Observer(
        "CO_obs",
        f"{ids.CO} * (1 dimensionless + 0 per_min * time)",
        U.l_per_min,
        "cardiac output",
    )
    d = _definition(
        "bl_pd_observer",
        _endogenous(baseline=Formula("0.006 mmole_per_min / CO_obs", U.mM)),
        pd=[pd],
        observers=[observer],
    )
    with pytest.raises(ValueError, match=r"CO \(a PD output.*through observer CO_obs"):
        compose(d)
    # without the PD component the same formula is valid
    build(
        _definition(
            "bl_no_pd",
            _endogenous(baseline=Formula(f"0.006 mmole_per_min / {ids.CO}", U.mM)),
        ),
        tmp_path,
    )


def test_the_baseline_is_validated() -> None:
    """A baseline which is not positive or a `Value` in another unit is refused."""
    with pytest.raises(ValueError, match="baseline must be positive"):
        _definition("bl_negative", _endogenous(baseline=-1.0)).validate()
    with pytest.raises(ValueError, match="baseline must be in"):
        _definition(
            "bl_unit", _endogenous(baseline=Value(0.05, U.mmole_per_min))
        ).validate()
    _definition("bl_value", _endogenous(baseline=Value(0.05, U.mM))).validate()


def test_the_baseline_is_annotated(tmp_path: Path) -> None:
    """`baseline_<s>` is a parameter in mM with an SBO term and a PBPKO term."""
    result = build(_definition("bl_annotated", _endogenous()), tmp_path)
    model = libsbml.readSBMLFromFile(str(result.comp_path)).getModel()
    p = model.getParameter(ids.baseline(S))
    assert p is not None
    assert p.getValue() == B
    assert p.getUnits() == U.mM.sid
    assert p.isSetSBOTerm()
    assert "PBPKO" in libsbml.XMLNode.convertXMLNodeToString(p.getAnnotation())


# The flat and comp model of `_unchanged()` built at pkpdsim f008d56 (before
# the baseline): the numbers of species, parameters, initial assignments,
# rules, reactions and compartments. Pinned as counts and values, not as a
# hash of the file, which every unrelated change of the library would move.
_PINNED = {
    "comp": (9, 236, 145, 12, 3, 9),
    "flat": (65, 404, 261, 68, 108, 77),
}


def _unchanged(baseline: float | None = None) -> PBPKDefinition:
    """The fixed small definition of the pinned counts: permeability, kinetic erythrocytes."""
    return _definition(
        "pin",
        _endogenous(
            distribution=Kp(2.0, mode=Mode.PERMEABILITY, ftissue=1.0),
            erythrocytes=Erythrocytes("kinetic", kbc=2.0, k_exchange=1.0),
            baseline=baseline,
        ),
    )


def _counts(model: libsbml.Model) -> tuple[int, ...]:
    """The numbers of species, parameters, initial assignments, rules, reactions, compartments."""
    return (
        model.getNumSpecies(),
        model.getNumParameters(),
        model.getNumInitialAssignments(),
        model.getNumRules(),
        model.getNumReactions(),
        model.getNumCompartments(),
    )


def test_a_substance_without_baseline_is_unchanged(tmp_path: Path) -> None:
    """Without a baseline: the pinned counts, every species at 0, no `baseline_` id.

    With a baseline the model differs only by the parameter `baseline_<s>`
    and the initial assignments of the species.
    """
    result = build(_unchanged(), tmp_path / "plain")
    comp = libsbml.readSBMLFromFile(str(result.comp_path))
    flat = libsbml.readSBMLFromFile(str(result.flat_path))
    assert _counts(comp.getModel()) == _PINNED["comp"]
    assert _counts(flat.getModel()) == _PINNED["flat"]
    model = flat.getModel()
    species = {sp.getId(): sp for sp in model.getListOfSpecies()}
    for sid, sp in species.items():
        value = (
            sp.getInitialConcentration()
            if sp.isSetInitialConcentration()
            else sp.getInitialAmount()
        )
        assert value == 0.0, sid
    assert not [
        ia for ia in model.getListOfInitialAssignments() if ia.getSymbol() in species
    ]
    assert not [
        el.getId()
        for el in [*model.getListOfParameters(), *comp.getModel().getListOfParameters()]
        if "baseline" in el.getId()
    ]

    with_baseline = build(_unchanged(baseline=B), tmp_path / "baseline")
    for path, plain in (
        (with_baseline.comp_path, result.comp_path),
        (with_baseline.flat_path, result.flat_path),
    ):
        assert _without_baseline(path) == _without_baseline(plain)


def _without_baseline(path: Path) -> str:
    """The SBML of a model without the baseline: its parameters and the starts of the species.

    The parameter `baseline_<s>` (top model and inputs of the organs), its
    port and replacements, and the initial assignments of the species are
    removed; the initial values of the species set to 0.
    """
    doc = libsbml.readSBMLFromFile(str(path))
    models = [doc.getModel()]
    comp = doc.getPlugin("comp")
    if comp is not None:
        models += [
            comp.getModelDefinition(k) for k in range(comp.getNumModelDefinitions())
        ]
    for model in models:
        species = {sp.getId() for sp in model.getListOfSpecies()}
        for ia in list(model.getListOfInitialAssignments()):
            if ia.getSymbol() in species or "baseline" in ia.getSymbol():
                model.removeInitialAssignment(ia.getSymbol())
        for p in list(model.getListOfParameters()):
            if p.getId().startswith("baseline_"):
                model.removeParameter(p.getId())
        plugin = model.getPlugin("comp")
        if plugin is not None:
            for k in reversed(range(plugin.getNumPorts())):
                if "baseline" in plugin.getPort(k).getId():
                    plugin.removePort(k)
        for sp in model.getListOfSpecies():
            if sp.isSetInitialConcentration() or not sp.isSetInitialAmount():
                sp.setInitialConcentration(0.0)
    return libsbml.writeSBMLToString(doc)


def test_the_sbml_without_baseline_is_that_of_f008d56(tmp_path: Path) -> None:
    """Without a baseline the comp and flat SBML are those of f008d56, line by line.

    The order of the elements counts: an input of an organ created earlier
    (a coefficient a start would read) moves lines of the comp model
    (`conftest.assert_reference_sbml`).
    """
    assert_reference_sbml(build(reference_definition(), tmp_path))


def test_dataclass_default() -> None:
    """A substance has no baseline by default."""
    s = dataclasses.replace(_endogenous(), baseline=None)
    assert s.baseline is None
    assert Substance("x", name="x", mr=1.0).baseline is None
