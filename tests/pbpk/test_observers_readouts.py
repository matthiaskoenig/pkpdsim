"""Tests of the readouts: organ elements, intestine totals and observers of constants.

An observer of a definition reads a species or an unbound concentration of
an organ by its flat id (`MU__Cti_drg`), exposed by the organ as a PD input
is (one rule); the intestine adds the lumen amount and the summed fluxes of
absorption, efflux and loss to the feces on request
(`PBPKDefinition.readouts`); an observer of constants only is a constant with
an initial assignment, in a unit of its own (per body surface area).
"""

import copy
import dataclasses
from pathlib import Path

import libsbml
import numpy as np
import pytest
from sbmlsim.result.scan import end_values, start_values
from sbmlsim.simulation import Simulation
from sbmlsim.simulator import Simulator
from sbmlutils.factory import ModelDefinition, Parameter, Reaction, Species

from pkpdsim.pbpk import (
    CustomReaction,
    Direction,
    Excretion,
    FirstOrder,
    Membrane,
    ModelParameter,
    Observer,
    Organ,
    OrganComponent,
    Part,
    PBPKDefinition,
    Route,
    SegmentedIntestine,
    SingleLumen,
    Transport,
    U,
    Variants,
    compose,
    ids,
)
from pkpdsim.pbpk.builder import build
from pkpdsim.pbpk.equilibration import accumulating_states, equilibrated_simulation
from pkpdsim.pbpk.observers import _flux_sum, intestine_readouts
from pkpdsim.pbpk.units import (
    EXTRA_UNITS,
    ML_PER_MIN_1P73M2,
    UNIT_LABELS,
    axis,
    unit_label,
)
from pkpdsim.pbpk.variants.segments import segments

from .conftest import (
    assert_reference_sbml,
    cardiac_pd,
    compiled,
    definition,
    individual,
    individuals,
    load,
    reference_definition,
    simple_drug,
)
from .test_pd import _organ_definitions, _organ_reader

S, M = "drg", "drg_m"
TIMES = np.linspace(0.0, 600.0, 61)
"""The time points of the curves (min)."""


def _model(path: Path) -> libsbml.Model:
    """Read the model of an SBML file."""
    return libsbml.readSBMLFromFile(str(path)).getModel()


def _iv(tmp_path: Path, **overrides: object) -> Path:
    """Build the simple definition with the IV route only, the flat model."""
    return build(definition(routes=[Route.iv(S)], **overrides), tmp_path).flat_path


def _simulate(path: Path, cols: list[str], dose: str | None = None) -> np.ndarray:
    """Simulate a dose of 10 mg (IV unless named), the columns without time."""
    r = individual(path)
    r.integrator.absolute_tolerance = 1e-12
    r.integrator.relative_tolerance = 1e-10
    r[dose or ids.iv_dose(S)] = 10.0
    return r.simulate(times=TIMES, selections=["time", *cols])[:, 1:]


# --- species and unbound concentrations of the organs -------------------------


def test_an_observer_reads_a_species_of_an_organ(tmp_path: Path) -> None:
    """An observer names a species of an organ by its flat id and equals it.

    The organ exposes the species as the parameter `C<o>_<part>_<s>` of the
    top model (`Cmu_ti_drg`), which the rule of the observer reads; the
    species keeps its id.
    """
    local = ids.species(Part.TISSUE, S)
    flat, top = ids.flat(Organ.MUSCLE, local), ids.organ_species(Organ.MUSCLE, local)
    observer = Observer("Cmu_twice", f"2 dimensionless * {flat}", U.mM, "twice")
    path = _iv(tmp_path, observers=[observer])
    model = _model(path)
    assert model.getSpecies(flat) is not None
    rule = model.getAssignmentRuleByVariable("Cmu_twice")
    assert libsbml.formulaToL3String(rule.getMath()) == f"2 dimensionless * {top}"
    s = _simulate(path, ["Cmu_twice", top, f"[{flat}]"])
    assert s[:, 2].max() > 0
    np.testing.assert_allclose(s[:, 1], s[:, 2], rtol=1e-12, atol=0.0)
    np.testing.assert_allclose(s[:, 0], 2 * s[:, 2], rtol=1e-12, atol=0.0)


def test_an_observer_and_a_pd_component_read_one_exposure(tmp_path: Path) -> None:
    """A PD input and an observer of the same species read one output of the organ."""
    flat = ids.flat(Organ.MUSCLE, ids.species(Part.TISSUE, S))
    top = ids.organ_species(Organ.MUSCLE, ids.species(Part.TISSUE, S))
    observer = Observer("Cmu_copy", flat, U.mM, "copy")
    path = _iv(tmp_path, observers=[observer], pd=[_organ_reader(flat)])
    model = _model(path)
    rules = [r.getVariable() for r in model.getListOfRules()]
    assert rules.count(top) == 1
    s = _simulate(path, ["Cmu_copy", "org__E", top])
    np.testing.assert_allclose(s[:, 0], s[:, 2], rtol=1e-12, atol=0.0)
    np.testing.assert_allclose(s[:, 1], 2 * s[:, 2], rtol=1e-12, atol=0.0)


def test_the_unbound_concentration_of_an_organ_part(tmp_path: Path) -> None:
    """The unbound concentration of an organ part is exposed as `Cu_<o>_<part>_<s>`.

    It is the fu-weighted total, `fu/Kp*C` of the merged tissue, for an
    observer and for a PD component alike.
    """
    unbound = ids.flat(Organ.LIVER, ids.unbound(Part.TISSUE, S))
    total = ids.flat(Organ.LIVER, ids.species(Part.TISSUE, S))
    top = ids.organ_species(Organ.LIVER, ids.unbound(Part.TISSUE, S))
    assert top == "Cu_li_ti_drg"
    observer = Observer("Cli_unbound", unbound, U.mM, "unbound liver")
    path = _iv(tmp_path, observers=[observer], pd=[_organ_reader(unbound)])
    r = individual(path)
    fu, kp = r[ids.fu(S)], r[ids.kp(S)]
    s = _simulate(path, ["Cli_unbound", top, "org__E", f"[{total}]"])
    assert s[:, 3].max() > 0
    np.testing.assert_allclose(s[:, 0], fu / kp * s[:, 3], rtol=1e-10, atol=0.0)
    np.testing.assert_allclose(s[:, 1], s[:, 0], rtol=1e-12, atol=0.0)
    np.testing.assert_allclose(s[:, 2], 2 * s[:, 0], rtol=1e-12, atol=0.0)


def _pool(local: str, amount: bool = True) -> CustomReaction:
    """A custom reaction of the muscle into a species `local` (an amount by default)."""
    pool = Species(
        local,
        "Vti",
        initialConcentration=None if amount else 0.0,
        initialAmount=0.0 if amount else None,
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=amount,
    )
    return CustomReaction(
        "POOL",
        Organ.MUSCLE,
        Part.TISSUE,
        f"Cti_drg -> {local}",
        "0.001 l_per_min * Cti_drg",
        species=[pool],
    )


def test_an_observer_reads_an_amount_of_an_organ(tmp_path: Path) -> None:
    """An amount species `A<x>` of an organ is exposed as `A<o>_<x>` in its unit (mmole)."""
    top = ids.organ_species(Organ.MUSCLE, "Apool_drg")
    assert top == "Amu_pool_drg"
    d = definition(
        routes=[Route.iv(S)],
        observers=[Observer("Apool_mg", f"{ids.mr(S)} * MU__Apool_drg", U.mg, "mg")],
    )
    d = dataclasses.replace(d, chemistry=[*d.chemistry, _pool("Apool_drg")])
    path = build(d, tmp_path).flat_path
    model = _model(path)
    # the flat model defines every unit once, so the parameter of an organ names
    # the unit of the model (`mmole`, not the copy `MU__mmole` of the flattening)
    assert model.getParameter(top).getUnits() == "mmole"
    s = _simulate(path, ["Apool_mg", top, "MU__Apool_drg"])
    assert s[-1, 2] > 0
    np.testing.assert_allclose(s[:, 1], s[:, 2], rtol=1e-12, atol=0.0)
    np.testing.assert_allclose(s[:, 0], 300.0 * s[:, 2], rtol=1e-12, atol=0.0)


@pytest.mark.parametrize(
    ("read", "match"),
    [
        ("MU__Vti", r"observer 'x' reads 'MU__Vti', which is no species of the organ"),
        ("MU__Car_drg", r"'MU__Car_drg', an element of the body .* read 'Car_drg'"),
        ("XX__Cti_drg", r"'x' reads \['XX__Cti_drg'\], which are no values"),
    ],
    ids=["compartment", "input of the organ", "unknown submodel"],
)
def test_an_observer_of_no_readable_element_of_an_organ_raises(
    read: str, match: str
) -> None:
    """An observer which reads no species or unbound concentration of an organ raises."""
    with pytest.raises(ValueError, match=match):
        compose(definition(observers=[Observer("x", read, U.mM, "x")]))


@pytest.mark.parametrize(
    ("local", "amount", "match"),
    [
        ("Xpool_drg", True, r"'MU__Xpool_drg', an amount .* no amount of the id"),
        ("Apool_drg", False, r"'MU__Apool_drg', a species .* no concentration"),
    ],
    ids=["amount", "concentration"],
)
def test_an_observer_of_a_species_outside_the_id_scheme_raises(
    local: str, amount: bool, match: str
) -> None:
    """An amount whose id is no `A...` or a concentration whose id is no `C...` is named."""
    d = definition(observers=[Observer("x", f"MU__{local}", U.mM, "x")])
    d = dataclasses.replace(d, chemistry=[*d.chemistry, _pool(local, amount)])
    with pytest.raises(ValueError, match=match):
        compose(d)


def test_an_observer_of_the_top_model_changes_no_organ() -> None:
    """Observers which read the top model only add nothing to the organs."""
    plain = _organ_definitions(compose(definition()))
    observer = Observer("Carm_twice", f"2 dimensionless * {ids.arm(S)}", U.mM, "x")
    assert _organ_definitions(compose(definition(observers=[observer]))) == plain


# --- observers of constants ---------------------------------------------------


def test_an_observer_of_constants_is_an_initial_assignment(tmp_path: Path) -> None:
    """An observer of constants is a constant with an initial assignment.

    The GFR per 1.73 m2 body surface area in ml/min/1.73 m2, a unit outside
    `U` which the model defines because the observer uses it; it equals its
    formula at the start after a change of the individual before the
    simulation. An observer of a constant observer is one too; an observer
    of a state stays a rule.
    """
    observers = [
        Observer(
            "GFRbsa",
            "GFR / BSA * 1730 ml_per_min_1p73m2 * 1 m2 / 1 l_per_min",
            ML_PER_MIN_1P73M2,
            "GFR per 1.73 m2",
        ),
        Observer("GFRbsa2", "2 dimensionless * GFRbsa", ML_PER_MIN_1P73M2, "x2"),
        Observer(
            "Carm_gfr", f"{ids.arm(S)} * GFRbsa / 1 ml_per_min_1p73m2", U.mM, "state"
        ),
    ]
    path = _iv(tmp_path, observers=observers)
    model = _model(path)
    unit = model.getUnitDefinition(ML_PER_MIN_1P73M2.sid)
    assert unit is not None
    for sid in ("GFRbsa", "GFRbsa2"):
        assert model.getParameter(sid).getConstant()
        assert model.getInitialAssignmentBySymbol(sid) is not None
        assert model.getAssignmentRuleByVariable(sid) is None
        assert model.getParameter(sid).getUnits() == ML_PER_MIN_1P73M2.sid
    assert model.getAssignmentRuleByVariable("Carm_gfr") is not None
    person = individuals(path)
    for changes in ({}, {"BW": 100.0, "HEIGHT": 190.0, "f_renal_function": 0.5}):
        r = person(**changes)
        expected = r[ids.GFR] / r[ids.BSA] * 1730.0
        if not changes:
            assert 60.0 < expected < 150.0  # a healthy GFR in ml/min/1.73 m2
        assert r["GFRbsa"] == pytest.approx(expected, rel=1e-12)
        assert r["GFRbsa2"] == pytest.approx(2 * expected, rel=1e-12)


def test_the_unit_per_body_surface_area() -> None:
    """`ml_per_min_1p73m2` is 1e-3/1.73 l/min/m2 and no unit of `U`."""
    doc = libsbml.SBMLDocument(3, 2)
    ML_PER_MIN_1P73M2.create_sbml(model=doc.createModel())
    unit = doc.getModel().getUnitDefinition("ml_per_min_1p73m2")
    kinds = {
        libsbml.UnitKind_toString(u.getKind()): (u.getMultiplier(), u.getExponent())
        for u in unit.getListOfUnits()
    }
    assert kinds["litre"][0] == pytest.approx(1e-3 / 1.73, rel=1e-12)
    assert kinds["second"] == (60.0, -1.0)
    assert kinds["metre"] == (1.0, -2.0)
    assert ML_PER_MIN_1P73M2.sid not in {
        getattr(u, "sid", None) for _, u in U.attributes()
    }


def test_the_unit_per_body_surface_area_is_defined_wherever_it_is_used(
    tmp_path: Path,
) -> None:
    """A unit outside `U` is a unit definition of every model which uses it.

    As the unit of a parameter of the definition, of a parameter of a custom
    reaction of an organ, or only as the unit of a number of a formula.
    """
    unit = ML_PER_MIN_1P73M2.sid
    reaction = dataclasses.replace(
        _pool("Apool_drg"),
        parameters=[Parameter("k_bsa", 1.0, ML_PER_MIN_1P73M2, constant=True)],
    )
    d = definition(
        routes=[Route.iv(S)],
        parameters=[
            ModelParameter("eGFR_ref", 90.0, ML_PER_MIN_1P73M2, "reference eGFR")
        ],
        observers=[
            Observer(
                "Carm_one",
                f"{ids.arm(S)} * 1 {unit} / 1 {unit}",
                U.mM,
                "arm vein",
            )
        ],
    )
    d = dataclasses.replace(d, chemistry=[*d.chemistry, reaction])
    result = build(d, tmp_path)
    comp = libsbml.readSBMLFromFile(str(result.comp_path))
    assert comp.getModel().getUnitDefinition(unit) is not None
    muscle = comp.getPlugin("comp").getModelDefinition("standard_mu")
    assert muscle.getUnitDefinition(unit) is not None
    flat = _model(result.flat_path)
    assert flat.getParameter("eGFR_ref").getUnits() == unit
    assert individual(result.flat_path)["eGFR_ref"] == 90.0
    # only in a formula
    only = definition(routes=[Route.iv(S)], observers=d.observers)
    model = _model(build(only, tmp_path / "formula").flat_path)
    assert model.getUnitDefinition(unit) is not None


def test_a_pd_component_gets_the_unit_per_area_it_uses(tmp_path: Path) -> None:
    """A PD component in the units of `U` which uses a unit outside `U` defines it.

    As the unit of a parameter of its model, as an organ does; a component
    which uses none keeps its units.
    """
    unit = ML_PER_MIN_1P73M2.sid
    plain = cardiac_pd()
    model = copy.copy(plain.model)
    model.parameters = [
        *(plain.model.parameters or []),
        Parameter("eGFR_pd", 90.0, ML_PER_MIN_1P73M2, constant=True),
    ]
    pd = dataclasses.replace(plain, model=model)
    result = build(definition(routes=[Route.iv(S)], pd=[pd]), tmp_path)
    comp = libsbml.readSBMLFromFile(str(result.comp_path))
    md = comp.getPlugin("comp").getModelDefinition(str(model.sid))
    assert md.getUnitDefinition(unit) is not None
    assert _model(result.flat_path).getParameter(f"{model.sid}__eGFR_pd") is not None
    plain_result = build(
        definition(routes=[Route.iv(S)], pd=[plain]), tmp_path / "plain"
    )
    comp = libsbml.readSBMLFromFile(str(plain_result.comp_path))
    md = comp.getPlugin("comp").getModelDefinition(str(plain.model.sid))
    assert md.getUnitDefinition(unit) is None


def test_an_observer_in_the_units_of_u_adds_no_unit(tmp_path: Path) -> None:
    """The unit definitions of a model with observers in `U` are those of `U`."""
    plain = _model(_iv(tmp_path / "plain"))
    observer = Observer("Carm_twice", f"2 dimensionless * {ids.arm(S)}", U.mM, "x")
    observed = _model(_iv(tmp_path / "observed", observers=[observer]))
    units = [u.getId() for u in plain.getListOfUnitDefinitions()]
    assert [u.getId() for u in observed.getListOfUnitDefinitions()] == units


# --- the readouts of the intestine --------------------------------------------


def _efflux() -> Transport:
    """The apical P-gp efflux of the small intestine, first order (l/min)."""
    return Transport(
        sid="PGP_DRG",
        organ=Organ.SMALL_INTESTINE,
        substance=S,
        direction=Direction.EFFLUX,
        membrane=Membrane.APICAL,
        kinetics=FirstOrder(k=0.5, per_volume=False),
        protein="ABCB1",
    )


READOUTS = tuple(
    f(S)
    for f in (ids.lumen_amount, ids.absorption_rate, ids.efflux_rate, ids.feces_rate)
)
"""The readouts of the intestine of the drug."""


def _oral(intestine: OrganComponent, **overrides: object) -> PBPKDefinition:
    """The drug with an oral solution, a P-gp efflux and fabs 0.8, without renal excretion."""
    drug = dataclasses.replace(simple_drug(), excretion=Excretion())
    fields: dict[str, object] = {
        "substances": [drug, *definition().substances[1:]],
        "routes": [Route.oral_solution(S)],
        "chemistry": [*definition().chemistry, _efflux()],
        "variants": Variants(organs={Organ.SMALL_INTESTINE: intestine}),
        "readouts": READOUTS,
    }
    fields.update(overrides)
    return definition(**fields)


def _reactions(model: libsbml.Model, reactant: set[str], product: set[str]) -> list:
    """The reactions with a reactant in `reactant` and a product in `product`."""
    return [
        r.getId()
        for r in model.getListOfReactions()
        if {x.getSpecies() for x in r.getListOfReactants()} & reactant
        and {x.getSpecies() for x in r.getListOfProducts()} & product
    ]


@pytest.mark.parametrize(
    ("intestine", "tokens"),
    [
        (SegmentedIntestine(fabs=0.8), ids.GUT_SEGMENTS),
        (SingleLumen(fabs=0.8), (None,)),
    ],
    ids=["segmented", "single_lumen"],
)
def test_the_intestine_readouts_are_the_sums_of_the_segments(
    tmp_path: Path, intestine: OrganComponent, tokens: tuple[str | None, ...]
) -> None:
    """The lumen amount is the sum of the segments, the rates the sums of the fluxes.

    `Alumen_<s> = sum V*C` of the lumen segments; `Rabs_<s>` the fluxes from
    the lumen into the wall, `Reff_<s>` those from the wall into the lumen
    (the apical efflux), `Rfeces_<s>` the fluxes into the feces, which
    integrate to the feces.
    """
    result = build(_oral(intestine), tmp_path)
    model = _model(result.flat_path)
    volumes = {
        ids.flat(Organ.SMALL_INTESTINE, ids.lumen(S, seg)): ids.flat(
            Organ.SMALL_INTESTINE, ids.lumen_volume(seg)
        )
        for seg in tokens
    }
    wall = {
        ids.flat(Organ.SMALL_INTESTINE, local)
        for local in (
            *(ids.mucosa(S, seg) for seg in tokens),
            *(ids.species(part, S) for part in Part),
        )
    }
    absorbed = _reactions(model, set(volumes), wall)
    effluxed = _reactions(model, wall, set(volumes))
    feces = _reactions(model, set(volumes), {ids.feces(S)})
    assert len(absorbed) == len([t for t in tokens if t != "st"])
    assert effluxed
    assert all(r.startswith("GU__PGP_DRG") for r in effluxed)
    assert feces
    cols = [*READOUTS, *(f"[{c}]" for c in volumes), *volumes.values()]
    r = load(result)
    r[ids.sol_dose(S)] = 10.0
    s = r.simulate(
        times=TIMES, selections=["time", *cols, *absorbed, *effluxed, *feces]
    )
    lumen = sum(s[f"[{c}]"] * s[v] for c, v in volumes.items())
    assert s[READOUTS[0]].max() > 0
    np.testing.assert_allclose(s[READOUTS[0]], lumen, rtol=1e-10, atol=1e-15)
    for readout, fluxes in zip(READOUTS[1:], (absorbed, effluxed, feces), strict=True):
        assert s[readout].max() > 0
        expected = sum(s[f] for f in fluxes)
        np.testing.assert_allclose(s[readout], expected, rtol=1e-10, atol=1e-15)
    r.reset()
    r[ids.sol_dose(S)] = 10.0
    fine = r.simulate(0, 600, 6001, selections=["time", READOUTS[3], ids.feces(S)])
    integral = np.trapezoid(fine[READOUTS[3]], fine["time"])
    assert integral == pytest.approx(fine[ids.feces(S)][-1], rel=1e-3)


def test_the_intestine_readouts_close_the_mass_balance_of_the_lumen(
    tmp_path: Path,
) -> None:
    """The dose is the lumen, what is absorbed net of the efflux and the feces.

    `D = SOLDOSE/Mr + Alumen + int(Rabs - Reff) + Afeces` at every time, a
    check of the readouts which is independent of how they select the fluxes.
    """
    r = load(build(_oral(SegmentedIntestine(fabs=0.8)), tmp_path))
    r[ids.sol_dose(S)] = 10.0
    mr = r[ids.mr(S)]
    s = r.simulate(
        0, 600, 6001, selections=["time", *READOUTS, ids.feces(S), ids.sol_dose(S)]
    )
    net = s[READOUTS[1]] - s[READOUTS[2]]
    absorbed = np.concatenate(
        [[0.0], np.cumsum(np.diff(s["time"]) * (net[1:] + net[:-1]) / 2)]
    )
    total = s[ids.sol_dose(S)] / mr + s[READOUTS[0]] + absorbed + s[ids.feces(S)]
    assert s[READOUTS[2]].max() > 0
    np.testing.assert_allclose(total, 10.0 / mr, rtol=1e-4)


def test_the_flux_sums_weight_the_stoichiometry() -> None:
    """A reaction counts with the net stoichiometry of the species it moves into."""
    md = ModelDefinition(
        "gu",
        reactions=[
            Reaction("R1", equation="Clumen_du_drg -> 2 Cmu_du_drg"),
            Reaction("R2", equation="Cmu_du_drg -> Clumen_du_drg"),
            Reaction("R3", equation="Clumen_rec_drg -> 3 Afeces_drg"),
            Reaction("R4", equation="Afeces_drg -> Clumen_rec_drg"),
        ],
    )
    lumen = {"Clumen_du_drg", "Clumen_rec_drg"}
    assert _flux_sum(md, "absorption", S, lumen) == "2 dimensionless * R1"
    assert _flux_sum(md, "efflux", S, lumen) == "R2"
    assert _flux_sum(md, "feces", S, lumen) == "3 dimensionless * R3 - R4"


def test_the_readouts_change_only_the_intestine() -> None:
    """The readouts add outputs to the intestine only; every other organ is unchanged."""
    plain = _organ_definitions(compose(_oral(SegmentedIntestine(), readouts=())))
    read = _organ_definitions(compose(_oral(SegmentedIntestine())))
    changed = {sid for sid in plain if read[sid] != plain[sid]}
    assert len(changed) == 1
    assert next(iter(changed)).endswith("_si")


@pytest.mark.parametrize(
    ("overrides", "match"),
    [
        ({"readouts": ("Xlumen_drg",)}, r"readout 'Xlumen_drg' is no readout"),
        ({"readouts": ("Alumen_drg", "Alumen_drg")}, r"duplicate readouts"),
        ({"readouts": ("Alumen_nope",)}, r"readout 'Alumen_nope' is no readout"),
    ],
    ids=["unknown", "twice", "unknown substance"],
)
def test_an_invalid_readout_is_refused_by_validate(
    overrides: dict[str, object], match: str
) -> None:
    """A readout which is no readout of a substance, or one given twice, is refused."""
    with pytest.raises(ValueError, match=match):
        _oral(SegmentedIntestine(), **overrides).validate()


def test_a_readout_the_intestine_cannot_compute_is_refused() -> None:
    """A readout without a lumen, of a substance not in the lumen or without efflux raises."""
    standard = Variants(organs={Organ.SMALL_INTESTINE: "standard"})
    with pytest.raises(ValueError, match=r"\['Alumen_drg'\] need an intestine with"):
        compose(
            definition(
                routes=[Route.iv(S)], variants=standard, readouts=("Alumen_drg",)
            )
        )
    with pytest.raises(ValueError, match=r"'Alumen_drg_m': .* no lumen species"):
        compose(_oral(SegmentedIntestine(), readouts=(ids.lumen_amount(M),)))
    no_efflux = _oral(
        SegmentedIntestine(), chemistry=definition().chemistry, readouts=("Reff_drg",)
    )
    with pytest.raises(ValueError, match=r"'Reff_drg': .* no flux"):
        compose(no_efflux)


def test_an_observer_of_an_unrequested_readout_names_the_request() -> None:
    """An observer which reads a readout the definition does not request is told so."""
    d = _oral(
        SegmentedIntestine(),
        readouts=(),
        observers=[Observer("x", "2 dimensionless * Alumen_drg", U.mmole, "x")],
    )
    with pytest.raises(ValueError, match=r"\['Alumen_drg'\].*PBPKDefinition.readouts"):
        compose(d)


def test_a_delta_of_an_organ_output(tmp_path: Path) -> None:
    """The urinary excretion rate and a readout, outputs of organs, take deltas."""
    d = definition(routes=[Route.iv(S)])
    drug = dataclasses.replace(d.substances[0], baseline=0.01)
    d = dataclasses.replace(
        d, substances=[drug, *d.substances[1:]], deltas=(ids.urine_rate(S),)
    )
    r = individual(build(d, tmp_path / "urine").flat_path)
    start = r[ids.urine_rate(S)]
    assert start > 0
    assert r[ids.delta_base(ids.urine_rate(S))] == pytest.approx(start, rel=1e-12)
    s = r.simulate(
        0, 600, 61, selections=[ids.urine_rate(S), ids.delta(ids.urine_rate(S))]
    )
    np.testing.assert_allclose(s[:, 1], s[:, 0] - start, rtol=1e-9, atol=1e-15)
    assert s[-1, 1] < 0
    oral = _oral(SegmentedIntestine(), deltas=(READOUTS[0],))
    r = individual(build(oral, tmp_path / "lumen").flat_path)
    r[ids.sol_dose(S)] = 10.0
    s = r.simulate(0, 600, 61, selections=[READOUTS[0], ids.delta(READOUTS[0])])
    assert s[:, 0].max() > 0
    np.testing.assert_allclose(s[:, 1], s[:, 0], rtol=1e-12, atol=1e-15)


def test_the_readouts_of_combined_definitions_are_their_union() -> None:
    """`combine` takes the readouts of every definition once."""
    a = _oral(SegmentedIntestine(), readouts=("Alumen_drg",))
    b = _oral(SegmentedIntestine(), readouts=("Alumen_drg", "Rabs_drg"))
    assert PBPKDefinition.combine(a, b, sid="both").readouts == (
        "Alumen_drg",
        "Rabs_drg",
    )


def test_gut_segments_are_the_segments_of_the_segmented_intestine() -> None:
    """`ids.GUT_SEGMENTS` are the segments of the segmented intestine in transit order."""
    assert tuple(seg.token for seg in segments()) == ids.GUT_SEGMENTS
    assert ids.GUT_SEGMENTS[0] == "st"
    assert ids.GUT_SEGMENTS[-1] == "rec"


# --- byte identity -------------------------------------------------------------


def test_without_readouts_the_model_is_unchanged(tmp_path: Path) -> None:
    """Empty readouts and no observers give the SBML of the stored reference."""
    d = dataclasses.replace(reference_definition(), readouts=(), observers=())
    assert_reference_sbml(build(d, tmp_path))


ABSORBED = (ids.absorption_rate(S), ids.absorbed_amount(S), ids.absorbed_fraction(S))
"""The absorption rate, the amount absorbed and the gross absorbed fraction."""


@pytest.mark.parametrize(
    ("intestine", "route", "above"),
    [
        (SegmentedIntestine(fabs=0.8), Route.oral_solution(S), 1.0),
        (SegmentedIntestine(fabs=0.8), Route.oral(S), 1.0),
        (SingleLumen(fabs=0.8), Route.oral_solution(S), 0.5),
    ],
    ids=["segmented solution", "segmented tablet", "single_lumen solution"],
)
def test_the_absorbed_amount_and_the_gross_absorbed_fraction(
    tmp_path: Path, intestine: OrganComponent, route: Route, above: float
) -> None:
    """`Aabs_<s>` integrates `Rabs_<s>`; `Fa_<s>` is it over the oral dose given.

    Two doses, 10 mg at the start and 5 mg after 5 hours: the fraction is
    over 10 mg until the second dose and over 15 mg after it, 0 without a
    dose; gross, so the drug which the efflux of the wall returns into the
    lumen and which is absorbed again counts again.
    """
    d = _oral(intestine, routes=[route], readouts=ABSORBED)
    r = load(build(d, tmp_path))
    dose = ids.po_dose(S) if route.kind == "oral" else ids.sol_dose(S)
    mr = r[ids.mr(S)]
    s = r.simulate(0, 600, 6001, selections=["time", *ABSORBED])
    assert np.all(s[ABSORBED[1]] == 0.0)
    assert np.all(s[ABSORBED[2]] == 0.0)
    r.reset()
    r[dose] = 10.0
    first = r.simulate(0, 300, 3001, selections=["time", *ABSORBED])
    r[dose] = r[dose] + 5.0
    second = r.simulate(300, 1440, 11401, selections=["time", *ABSORBED])
    for part, given in ((first, 10.0), (second, 15.0)):
        np.testing.assert_allclose(
            part[ABSORBED[2]][1:],
            part[ABSORBED[1]][1:] / (given / mr),
            rtol=1e-9,
            atol=1e-14,
        )
    rate = np.concatenate([first[ABSORBED[0]], second[ABSORBED[0]][1:]])
    time = np.concatenate([first["time"], second["time"][1:]])
    integral = np.trapezoid(rate, time)
    assert second[ABSORBED[1]][-1] == pytest.approx(integral, rel=1e-4)
    # the efflux of the wall returns drug into the lumen, which is absorbed
    # again: the gross fraction counts it every time, above 1 for the segments
    # (the efflux of every mucosa segment), above 0.5 for the single lumen
    assert second[ABSORBED[2]][-1] > above
    assert np.all(np.diff(first[ABSORBED[1]]) >= -1e-15)


def test_the_absorbed_fraction_alone_keeps_its_amount_in_the_intestine(
    tmp_path: Path,
) -> None:
    """`Fa_<s>` alone: the amount and the applied dose are states of the intestine only.

    The equilibration sets them back as the sinks, they accumulate.
    """
    result = build(
        _oral(SegmentedIntestine(), readouts=(ids.absorbed_fraction(S),)), tmp_path
    )
    doc = libsbml.readSBMLFromFile(str(result.flat_path))
    found = {p.getId() for p in doc.getModel().getListOfParameters()}
    gu = {
        ids.flat(Organ.SMALL_INTESTINE, ids.absorbed_amount(S)),
        ids.flat(Organ.SMALL_INTESTINE, ids.oral_dose_applied(S)),
    }
    assert gu <= found
    assert ids.absorbed_fraction(S) in found
    assert ids.absorbed_amount(S) not in found
    assert gu <= set(accumulating_states(result.flat_path.read_text(encoding="utf-8")))


@pytest.mark.parametrize(
    "readouts",
    [(ids.absorbed_fraction(S),), (ids.absorbed_amount(S), ids.absorbed_fraction(S))],
    ids=["fraction alone", "amount and fraction"],
)
def test_the_equilibration_sets_back_the_absorbed_amount_of_an_endogenous_substance(
    tmp_path: Path, readouts: tuple[str, ...]
) -> None:
    """An endogenous substance in the lumen at baseline: its readouts start at 0.

    The apical efflux of the wall moves the endogenous substance into the
    lumen at its baseline, and it is absorbed again, so its absorbed amount
    grows during the equilibration without a steady state. The equilibration
    sets it and the oral dose applied back at the start, as the sinks, so the
    simulation starts at 0 and needs no reset.

    After the equilibration the lumen is at its steady state with the wall,
    so the amount which is absorbed again over the simulation is the
    difference of two nearly equal fluxes: 4.6e-18 mmole over 60 min and
    3.7e-16 over 6000 min, the noise of the integrator, whose sign is that of
    the platform (-1.4e-17 on a macos runner). The test therefore asserts
    that the readout starts at exactly 0 and stays at 0 within that noise,
    which is what the setting back has to achieve; that it grows at all is
    the assertion `before > 1e-6` of the simulation without the
    equilibration, over 1e5 min.
    """
    drug = dataclasses.replace(simple_drug(), excretion=Excretion(), baseline=0.001)
    d = _oral(
        SegmentedIntestine(),
        substances=[drug, *definition().substances[1:]],
        readouts=readouts,
    )
    path = build(d, tmp_path).flat_path
    states = set(accumulating_states(path.read_text(encoding="utf-8")))
    amount, applied = (
        sid
        for local in (ids.absorbed_amount(S), ids.oral_dose_applied(S))
        for sid in (local, ids.flat(Organ.SMALL_INTESTINE, local))
        if sid in states
    )
    simulation = equilibrated_simulation(
        Simulation(end=60, preinit_changes={ids.sol_dose(S): 0.0}), path
    )
    assert simulation.changes[0].values[amount] == 0.0
    assert simulation.changes[0].values[applied] == 0.0
    r = compiled(path)
    r.integrator.absolute_tolerance = 1e-14
    r.integrator.relative_tolerance = 1e-10
    before = r.simulate(0, 1e5, 2, selections=["time", amount])[amount][-1]
    assert before > 1e-6
    result = Simulator(
        n_workers=1, absolute_tolerance=1e-14, relative_tolerance=1e-10
    ).run(path, simulation)
    assert float(start_values(result, amount).item()) == 0.0
    assert float(start_values(result, applied).item()) == 0.0
    assert abs(float(end_values(result, amount).item())) < 1e-12


def test_the_absorbed_fraction_needs_an_oral_dose() -> None:
    """A substance in the lumen without an oral dose state has no absorbed fraction of it."""
    md = ModelDefinition(
        "gu",
        species=[
            Species("Clumen_du_drg", compartment="Vlumen_du", initialConcentration=0.0)
        ],
        reactions=[Reaction("J_abs_du_drg", equation="Clumen_du_drg -> Cmu_du_drg")],
    )
    with pytest.raises(ValueError, match=r"'Fa_drg': the substance 'drg' has no oral"):
        intestine_readouts(md, (ids.absorbed_fraction(S),), (S,))


def test_the_unit_per_area_has_a_text_on_the_figures() -> None:
    """An axis in `ml_per_min_1p73m2` reads `[ml/min/1.73 m²]`, its unit stays the id.

    sbmlsim writes a unit id of the model as it is; the axis helper of the
    library writes the text of the unit and keeps the unit of the conversion.
    """
    from sbmlsim.plot import Axis

    assert Axis(label="GFR", unit=ML_PER_MIN_1P73M2.sid).name == (
        "GFR [ml_per_min_1p73m2]"
    )
    gfr = axis("GFR", ML_PER_MIN_1P73M2.sid, scale="log")
    assert gfr.name == "GFR [ml/min/1.73\u00a0m²]"
    assert gfr.unit == ML_PER_MIN_1P73M2.sid
    assert gfr.label == "GFR"
    assert "log" in str(gfr.scale).lower()
    assert unit_label(ML_PER_MIN_1P73M2.sid) == UNIT_LABELS["ml_per_min_1p73m2"]
    assert axis("time", "hr").name == Axis(label="time", unit="hr").name
    assert axis("fraction", "dimensionless").name == "fraction"
    assert set(UNIT_LABELS) == {u.sid for u in EXTRA_UNITS}


def test_the_text_of_a_unit_does_not_wrap() -> None:
    """The words of a unit text are joined by a non-breaking space.

    sbmlsim wraps a label longer than its axis at its spaces (`_wrap_label`
    of its matplotlib serializer), which broke `ml/min/1.73 m²` between
    `1.73` and `m²`. The label of an axis is wrapped by sbmlsim's own
    wrapper on axes from one character to the full label long: the unit
    text stays on one line.
    """
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from sbmlsim.plot.serialization_matplotlib import MatplotlibFigureSerializer

    fig = Figure()
    canvas = FigureCanvasAgg(fig)
    renderer = canvas.get_renderer()
    ax = fig.add_subplot()
    for unit, text in UNIT_LABELS.items():
        assert not any(c.isspace() and c != "\u00a0" for c in text), unit
        name = axis("creatinine clearance", unit).name
        assert name is not None
        ax.set_xlabel(name)
        full = ax.xaxis.label.get_window_extent(renderer).width
        wrapped = set()
        for k in range(1, 41):
            ax.set_xlabel(name)
            MatplotlibFigureSerializer._wrap_label(
                ax.xaxis.label, full * k / 40, False, renderer
            )
            lines = ax.xaxis.label.get_text().split("\n")
            assert any(f"[{text}]" in line for line in lines), (unit, k, lines)
            wrapped.add(len(lines))
        # the label is wrapped into several lines on the short axes
        assert max(wrapped) > 1, unit
    assert UNIT_LABELS["ml_per_min_1p73m2"] == "ml/min/1.73\u00a0m²"
