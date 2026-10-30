"""Tests of the observers: assignment rules of the top model, no states.

The observers read the pool species and the outflow concentrations of the
organs: the unbound and the blood concentration of every pool, the hepatic
vein and the arm vein, the flow-weighted mix of the venous outflows of the
organs of a preset of `data/arm_vein.toml`. The kidney exposes the urinary
excretion rate `Rurine_<s>`, and a definition adds observers of its own.
"""

import dataclasses
from pathlib import Path
from typing import TYPE_CHECKING

import libsbml
import numpy as np
import pytest

from pkpdsim.pbpk import (
    CustomReaction,
    Direction,
    Erythrocytes,
    FirstOrder,
    Kp,
    Membrane,
    Mode,
    Observer,
    Organ,
    Part,
    PBPKDefinition,
    Pool,
    Route,
    Transport,
    U,
    Variants,
    compose,
    ids,
    load,
)
from pkpdsim.pbpk.builder import build
from pkpdsim.pbpk.observers import arm_vein_weights

from .conftest import definition, individual, simple_drug

if TYPE_CHECKING:
    import roadrunner

S = "drg"
"""The substance of the tests."""

TIMES = np.linspace(0.0, 240.0, 49)
"""The time points of the curves (min)."""

DEFAULT = arm_vein_weights("pksim_default")
"""The weights of the default preset of the arm vein, skin and muscle."""


def _run(path: Path, cols: list[str], **changes: float) -> np.ndarray:
    """Simulate an IV bolus of 10 mg, the columns without time."""
    r = individual(path, **changes)
    r.integrator.absolute_tolerance = 1e-12
    r.integrator.relative_tolerance = 1e-10
    r[ids.iv_dose(S)] = 10.0
    return r.simulate(times=TIMES, selections=["time", *cols])[:, 1:]


def _build(tmp_path: Path, **overrides: object) -> Path:
    """Build the simple definition with IV route only."""
    d = definition(routes=[Route.iv(S)], **overrides)
    return build(d, tmp_path).flat_path


def _rr(path: Path) -> "roadrunner.RoadRunner":
    """Load a model without changes."""
    return individual(path)


def test_arm_vein_weighted_outflow(tmp_path: Path) -> None:
    """`Carm` is the default mix `w_sk*Cout_sk + w_mu*Cout_mu` of skin and muscle.

    The PK-Sim documentation gives about 70 % skin and 30 % muscle.
    """
    assert set(DEFAULT) == {Organ.SKIN, Organ.MUSCLE}
    w_sk, w_mu = DEFAULT[Organ.SKIN].value, DEFAULT[Organ.MUSCLE].value
    assert (w_sk, w_mu) == (0.7, 0.3)
    path = _build(tmp_path)
    sk, mu = ids.organ_outflow(Organ.SKIN, S), ids.organ_outflow(Organ.MUSCLE, S)
    carm, csk, cmu = _run(path, [ids.arm(S), sk, mu]).T
    np.testing.assert_allclose(carm, w_sk * csk + w_mu * cmu, rtol=1e-12, atol=0)
    assert np.max(carm) > 0


@pytest.mark.parametrize("preset", ["pksim_default", "levitt2004", "musther2015"])
def test_arm_vein_presets(tmp_path: Path, preset: str) -> None:
    """Every preset is a mix of outflows (weights sum to 1) which `Carm` follows."""
    weights = arm_vein_weights(preset)
    assert sum(w.value for w in weights.values()) == pytest.approx(1.0, abs=1e-12)
    path = _build(tmp_path, arm_vein=preset)
    cols = [
        f"[{ids.pool_species(site, S)}]"
        if isinstance(site, Pool)
        else ids.organ_outflow(site, S)
        for site in weights
    ]
    result = _run(path, [ids.arm(S), *cols])
    expected = sum(w.value * result[:, i + 1] for i, w in enumerate(weights.values()))
    np.testing.assert_allclose(result[:, 0], expected, rtol=1e-12, atol=0)


def test_arterialised_is_arterial(tmp_path: Path) -> None:
    """The arterialised hand vein is the arterial plasma."""
    path = _build(tmp_path, arm_vein="arterialised")
    carm, car = _run(path, [ids.arm(S), f"[{ids.pool_species(Pool.ARTERIAL, S)}]"]).T
    np.testing.assert_allclose(carm, car, rtol=1e-12, atol=0)


def test_arm_vein_weights_are_parameters(tmp_path: Path) -> None:
    """The weights are parameters with their source, which a model may change."""
    path = _build(tmp_path)
    r = _rr(path)
    for site, w in DEFAULT.items():
        assert r[ids.arm_weight(site)] == w.value
    doc = libsbml.readSBMLFromFile(str(path))
    notes = doc.getModel().getParameter(ids.arm_weight(Organ.SKIN)).getNotesString()
    assert "PK-Sim" in notes


def test_unknown_arm_vein_preset_raises() -> None:
    """The preset is validated with the definition."""
    with pytest.raises(ValueError, match="arm_vein 'nope'"):
        definition(arm_vein="nope").validate()


def test_arm_vein_is_observer(tmp_path: Path) -> None:
    """The observer changes no state: the curves are the same with every preset.

    No species, reaction or rate rule reads an observer; it is an assignment
    rule which only other observers read.
    """
    species = [
        f"[{ids.pool_species(Pool.VENOUS, S)}]",
        "[LI__Cti_drg]",
        "[MU__Cti_drg]",
    ]
    curves = [
        _run(_build(tmp_path / preset, arm_vein=preset), species)
        for preset in ("pksim_default", "arterialised")
    ]
    np.testing.assert_array_equal(curves[0], curves[1])

    path = _build(tmp_path / "check")
    model = libsbml.readSBMLFromFile(str(path)).getModel()
    observers = {
        f(S)
        for f in (
            ids.arm,
            ids.arm_free,
            ids.arm_blood,
            ids.hepatic_vein,
            ids.hepatic_vein_free,
            ids.hepatic_vein_blood,
        )
    }
    observers |= {
        ids.flat(o, ids.unbound(p, S))
        for o in Organ
        for p in (Part.TISSUE, Part.INTERSTITIAL, Part.INTRACELLULAR, Part.ORGAN)
    }
    observers |= {f(p, S) for p in Pool for f in (ids.free, ids.blood)}
    assert all(model.getSpecies(o) is None for o in observers)
    readers = [
        libsbml.formulaToL3String(r.getKineticLaw().getMath())
        for r in model.getListOfReactions()
    ] + [
        libsbml.formulaToL3String(rule.getMath())
        for rule in model.getListOfRules()
        if rule.getVariable() not in observers
    ]
    for formula in readers:
        names = set(formula.replace("(", " ").replace(")", " ").split())
        assert not names & observers, formula


@pytest.mark.parametrize(
    "ery",
    [
        None,
        Erythrocytes("equilibrium", bp=1.3),
        Erythrocytes("kinetic", kbc=2.0, k_exchange=1.0),
    ],
    ids=["plasma", "equilibrium", "kinetic"],
)
def test_free_and_blood_observers(tmp_path: Path, ery: Erythrocytes | None) -> None:
    """Unbound `fu*C` and blood `(1-HCT)*C + HCT*Cbc` of the pools and the arm vein.

    Without kinetic blood cells the blood cells hold `Kbc*C`, the blood
    `BP*C` (`BP = 1 - HCT` without erythrocyte distribution); kinetic blood
    cells are their species, at the arm vein the flow-weighted mix of the
    blood cell outflows of skin and muscle.
    """
    drug = dataclasses.replace(simple_drug(), erythrocytes=ery)
    path = _build(tmp_path, substances=[drug, *definition().substances[1:]])
    r = _rr(path)
    fu, bp, hct = r[ids.fu(S)], r[ids.bp(S)], r[ids.HCT]
    for pool in Pool:
        c = f"[{ids.pool_species(pool, S)}]"
        cols = [c, ids.free(pool, S), ids.blood(pool, S)]
        if ery is not None and ery.mode == "kinetic":
            cols.append(f"[{ids.pool_species(pool, S, part=Part.BLOOD_CELLS)}]")
        result = _run(path, cols)
        np.testing.assert_allclose(result[:, 1], fu * result[:, 0], rtol=1e-12)
        if ery is not None and ery.mode == "kinetic":
            blood = (1 - hct) * result[:, 0] + hct * result[:, 3]
        else:
            blood = bp * result[:, 0]
        np.testing.assert_allclose(result[:, 2], blood, rtol=1e-12)

    cols = [ids.arm(S), ids.arm_free(S), ids.arm_blood(S)]
    if ery is not None and ery.mode == "kinetic":
        cols += [
            ids.organ_outflow_bc(Organ.SKIN, S),
            ids.organ_outflow_bc(Organ.MUSCLE, S),
        ]
    result = _run(path, cols)
    np.testing.assert_allclose(result[:, 1], fu * result[:, 0], rtol=1e-12)
    if ery is not None and ery.mode == "kinetic":
        cbc = (
            DEFAULT[Organ.SKIN].value * result[:, 3]
            + DEFAULT[Organ.MUSCLE].value * result[:, 4]
        )
        blood = (1 - hct) * result[:, 0] + hct * cbc
        assert np.max(cbc) > 0
    else:
        blood = bp * result[:, 0]
    np.testing.assert_allclose(result[:, 2], blood, rtol=1e-12)

    f_shunt = 0.3
    cols = [ids.hepatic_vein(S), ids.hepatic_vein_free(S), ids.hepatic_vein_blood(S)]
    if ery is not None and ery.mode == "kinetic":
        cols += [
            ids.organ_outflow_bc(Organ.LIVER, S),
            f"[{ids.pool_species(Pool.ARTERIAL, S, Part.BLOOD_CELLS)}]",
            f"[{ids.pool_species(Pool.PORTAL, S, Part.BLOOD_CELLS)}]",
        ]
    result = _run(path, cols, f_shunt=f_shunt)
    np.testing.assert_allclose(result[:, 1], fu * result[:, 0], rtol=1e-12)
    if ery is not None and ery.mode == "kinetic":
        qh, qha, qpo = r[ids.QH], r[ids.QHA], r[ids.QPO]
        cbc = (
            (1 - f_shunt) * qh * result[:, 3]
            + f_shunt * (qha * result[:, 4] + qpo * result[:, 5])
        ) / qh
        blood = (1 - hct) * result[:, 0] + hct * cbc
        assert np.max(cbc) > 0
    else:
        blood = bp * result[:, 0]
    np.testing.assert_allclose(result[:, 2], blood, rtol=1e-12)


def test_unbound_tissue_observers(tmp_path: Path) -> None:
    """Every organ has the unbound concentration of each tissue part of a substance.

    `Cu_<part>_<s> = fu/K*C` with the coefficient of the part to plasma: the
    merged tissue `Kp`, interstitium and cells `K_int_pls`, `K_cell_pls`, the
    flow-limited organ `K_or`.
    """
    full = dataclasses.replace(
        simple_drug(), distribution=Kp(2.0, mode=Mode.PERMEABILITY, ftissue=1.0)
    )
    perfusion = dataclasses.replace(
        simple_drug(), distribution=Kp(2.0, mode=Mode.PERFUSION)
    )
    cases = [
        (simple_drug(), {Part.TISSUE: ids.kp(S)}),
        (
            full,
            {
                Part.INTERSTITIAL: ids.k_int_pls(S),
                Part.INTRACELLULAR: ids.k_cell_pls(S),
            },
        ),
        (perfusion, {Part.ORGAN: ids.k_organ(S)}),
    ]
    for i, (drug, coefficients) in enumerate(cases):
        path = _build(tmp_path / str(i), substances=[drug], chemistry=[])
        r = _rr(path)
        fu = r[ids.fu(S)]
        parameters = set(r.model.getGlobalParameterIds())
        for organ in (Organ.LIVER, Organ.MUSCLE):
            for part, k in coefficients.items():
                c = f"[{ids.flat(organ, ids.species(part, S))}]"
                cu = ids.flat(organ, ids.unbound(part, S))
                coefficient = r[k] if k == ids.kp(S) else r[ids.flat(organ, k)]
                result = _run(path, [c, cu])
                np.testing.assert_allclose(
                    result[:, 1], fu / coefficient * result[:, 0], rtol=1e-12
                )
                assert np.max(result[:, 1]) > 0
            absent = {
                Part.TISSUE,
                Part.INTERSTITIAL,
                Part.INTRACELLULAR,
                Part.ORGAN,
            } - set(coefficients)
            for part in absent:
                assert ids.flat(organ, ids.unbound(part, S)) not in parameters


def test_arm_vein_data_sources() -> None:
    """Every preset of the data file names its sources; the weights are fractions."""
    data = load("arm_vein")
    assert {"pksim_default", "levitt2004", "musther2015", "arterialised"} <= set(data)
    for preset in data:
        for w in arm_vein_weights(preset).values():
            assert 0 < w.value <= 1
            assert w.source


def test_definition_arm_vein_default() -> None:
    """The PK-Sim default is the preset of a definition which names none."""
    assert PBPKDefinition(sid="x", substances=[]).arm_vein == "pksim_default"


# --- the urinary excretion rate ----------------------------------------------


_SECRETION = Transport(
    "SEC", Organ.KIDNEY, S, Direction.EFFLUX, Membrane.APICAL, FirstOrder(k=0.05)
)
"""The tubular secretion of the drug, an apical efflux of the kidney."""


def _urine_run(path: Path, cols: list[str]) -> np.ndarray:
    """Simulate an IV bolus of 10 mg for 10 hours, with time."""
    r = individual(path)
    r.integrator.absolute_tolerance = 1e-12
    r.integrator.relative_tolerance = 1e-10
    r[ids.iv_dose(S)] = 10.0
    return r.simulate(0, 600, 6001, selections=["time", *cols])


def test_urine_rate_integrates_to_the_urine(tmp_path: Path) -> None:
    """The urinary excretion rate integrates to the amount in the urine."""
    s = _urine_run(_build(tmp_path), [ids.urine_rate(S), ids.urine(S)])
    assert s[-1, 2] > 0
    assert np.trapezoid(s[:, 1], s[:, 0]) == pytest.approx(s[-1, 2], rel=1e-3)


def test_urine_rate_sums_filtration_and_secretion(tmp_path: Path) -> None:
    """`Rurine_<s>` is the filtration plus the apical efflux of the kidney."""
    path = _build(tmp_path, chemistry=[_SECRETION])
    gfr = ids.flat(Organ.KIDNEY, ids.filtration_reaction(S))
    sec = ids.flat(Organ.KIDNEY, "SEC")
    s = _urine_run(path, [ids.urine_rate(S), gfr, sec])
    assert np.all(s[1:, 3] > 0)
    np.testing.assert_allclose(s[:, 1], s[:, 2] + s[:, 3], rtol=1e-12, atol=0)


_REABSORPTION = CustomReaction(
    "REABS",
    Organ.KIDNEY,
    Part.PLASMA,
    f"Aurine_{S} -> Cpl_{S} [Car_{S}]",
    f"0.5 dimensionless * fu_{S} * GFR * Car_{S}",
    [],
)
"""A reabsorption of half the filtrate, a custom reaction from the urine of the kidney."""


def test_urine_rate_is_net_of_a_reaction_from_the_urine(tmp_path: Path) -> None:
    """A custom reaction which takes from the urine is subtracted from `Rurine_<s>`.

    The rate is the net flux into the urine, which integrates to the amount
    in the urine: the filtration less the reabsorption, half of it.
    """
    path = _build(tmp_path, chemistry=[_REABSORPTION])
    gfr = ids.flat(Organ.KIDNEY, ids.filtration_reaction(S))
    reabs = ids.flat(Organ.KIDNEY, "REABS")
    s = _urine_run(path, [ids.urine_rate(S), gfr, reabs, ids.urine(S)])
    assert np.all(s[1:, 3] > 0)
    np.testing.assert_allclose(s[:, 1], s[:, 2] - s[:, 3], rtol=1e-12, atol=1e-15)
    np.testing.assert_allclose(s[:, 1], 0.5 * s[:, 2], rtol=1e-9, atol=1e-15)
    assert np.trapezoid(s[:, 1], s[:, 0]) == pytest.approx(s[-1, 4], rel=1e-3)


def test_urine_rate_is_unchanged_without_a_reaction_from_the_urine(
    tmp_path: Path,
) -> None:
    """Without a reaction from the urine the rule of `Rurine_<s>` is that of before.

    The formula sums the fluxes into the urine and the notes end as before
    the net rate (the SBML of the existing models is byte-identical).
    """
    model = libsbml.readSBMLFromFile(
        str(_build(tmp_path, chemistry=[_SECRETION]))
    ).getModel()
    rule = model.getAssignmentRuleByVariable(ids.urine_rate(S))
    assert libsbml.formulaToL3String(rule.getMath()) == "KI__J_gfr_drg + KI__SEC"
    notes = model.getParameter(ids.urine_rate(S)).getNotesString()
    assert "the sum of the fluxes of the kidney into the urine." in notes
    assert "less those out of it" not in notes


def test_urine_rate_of_the_nephron_kidney(tmp_path: Path) -> None:
    """In the nephron kidney `Rurine_<s>` is the outflow of the collecting duct."""
    path = _build(
        tmp_path,
        chemistry=[_SECRETION],
        variants=Variants(organs={Organ.KIDNEY: "nephron"}),
    )
    cd = ids.flat(Organ.KIDNEY, ids.tubule_reaction(S, "cd"))
    s = _urine_run(path, [ids.urine_rate(S), cd, ids.urine(S)])
    np.testing.assert_allclose(s[:, 1], s[:, 2], rtol=1e-12, atol=0)
    assert np.trapezoid(s[:, 1], s[:, 0]) == pytest.approx(s[-1, 3], rel=1e-3)


def test_urine_rate_only_for_the_urine(tmp_path: Path) -> None:
    """A substance which reaches no urine has no urinary excretion rate."""
    model = libsbml.readSBMLFromFile(str(_build(tmp_path))).getModel()
    assert model.getParameter(ids.urine_rate(S)) is not None
    assert model.getParameter(ids.urine_rate("drg_m")) is None


# --- the changes of observers over their start --------------------------------


def _endogenous_drug() -> PBPKDefinition:
    """The simple definition with a baseline of the drug (it decays from it)."""
    d = definition(routes=[Route.iv(S)])
    drug = dataclasses.replace(d.substances[0], baseline=0.01)
    return dataclasses.replace(d, substances=[drug, *d.substances[1:]])


def test_a_delta_is_the_change_of_an_observer_over_its_start(tmp_path: Path) -> None:
    """`<id>_delta` is `<id> - <id>_0` with the value of the observer at the start.

    The arm vein of an endogenous substance and an observer of the definition:
    the change is 0 at the start, also for an individual changed before the
    simulation (the start value is an initial assignment of the observer),
    and follows the observer after it.
    """
    d = dataclasses.replace(
        _endogenous_drug(),
        observers=[
            Observer("Ctwice_drg", f"2 dimensionless * {ids.arm(S)}", U.mM, "x")
        ],
        deltas=(ids.arm(S), "Ctwice_drg"),
    )
    path = build(d, tmp_path).flat_path
    for changes in ({}, {"BW": 110.0, "f_renal_function": 0.3}):
        r = individual(path, **changes)
        r.integrator.absolute_tolerance = 1e-14
        r.integrator.relative_tolerance = 1e-10
        start = r[ids.arm(S)]
        assert start > 0
        assert r[ids.delta_base(ids.arm(S))] == pytest.approx(start, rel=1e-12)
        # 0 to the rounding of the observer, relative to its magnitude (the
        # rounding differs with the JIT level)
        atol = 1e-12 * start
        assert r[ids.delta(ids.arm(S))] == pytest.approx(0.0, abs=atol)
        cols = [ids.arm(S), ids.delta(ids.arm(S)), ids.delta("Ctwice_drg")]
        s = r.simulate(0, 600, 61, selections=cols)
        np.testing.assert_allclose(s[:, 1], s[:, 0] - start, rtol=1e-9, atol=atol)
        np.testing.assert_allclose(s[:, 2], 2 * s[:, 1], rtol=1e-9, atol=2 * atol)
        assert s[-1, 1] < 0  # the drug decays from its baseline


def test_a_delta_of_no_observer_is_refused() -> None:
    """A delta which is no observer of the top model or whose ids exist raises, naming them."""
    with pytest.raises(ValueError, match=r"\['BW', 'Cnone_drg'\]"):
        compose(dataclasses.replace(_endogenous_drug(), deltas=("BW", "Cnone_drg")))
    d = dataclasses.replace(
        _endogenous_drug(),
        observers=[Observer("Carm_drg_0", f"{ids.arm(S)}", U.mM, "x")],
        deltas=(ids.arm(S),),
    )
    with pytest.raises(ValueError, match="Carm_drg_0"):
        compose(d)
    with pytest.raises(ValueError, match="duplicate deltas"):
        dataclasses.replace(_endogenous_drug(), deltas=("A", "A")).validate()


def test_without_deltas_the_model_is_unchanged(tmp_path: Path) -> None:
    """Empty deltas give the SBML of before the field (the reference of f008d56).

    The reference definition with `deltas=()` builds the stored comp and flat
    SBML of f008d56, which predates the deltas, line by line.
    """
    from .conftest import assert_reference_sbml, reference_definition

    d = dataclasses.replace(reference_definition(), deltas=())
    assert_reference_sbml(build(d, tmp_path))


# --- the observers of a definition -------------------------------------------


def test_observer_sums_two_sinks(tmp_path: Path) -> None:
    """An observer of the definition is an assignment rule of the top model."""
    observer = Observer(
        "Aurine_sum", "Aurine_drg + Aurine_drg_m", U.mmole, "drug and metabolite"
    )
    path = _build(tmp_path, observers=[observer])
    cols = ["Aurine_sum", ids.urine(S), ids.urine("drg_m")]
    s = _urine_run(path, cols)
    assert s[-1, 1] > 0
    np.testing.assert_allclose(s[:, 1], s[:, 2] + s[:, 3], rtol=1e-12, atol=0)
    model = libsbml.readSBMLFromFile(str(path)).getModel()
    parameter = model.getParameter("Aurine_sum")
    assert not parameter.getConstant()
    assert parameter.getName() == "drug and metabolite"


def test_observer_reads_another_observer(tmp_path: Path) -> None:
    """An observer may read an observer of the library and one of the definition."""
    observers = [
        Observer("Carm_uM", "1000 dimensionless * Carm_drg", U.mM, "arm vein x 1000"),
        Observer("Carm_ratio", "Carm_uM / Cve_drg", U.dimensionless, "ratio"),
    ]
    path = _build(tmp_path, observers=observers)
    s = _urine_run(
        path, ["Carm_ratio", ids.arm(S), f"[{ids.pool_species(Pool.VENOUS, S)}]"]
    )
    np.testing.assert_allclose(s[1:, 1], 1000 * s[1:, 2] / s[1:, 3], rtol=1e-10, atol=0)


def test_observer_of_constants_is_a_constant() -> None:
    """An observer which reads no state is a constant with an initial assignment (R53).

    So is an observer which reads only observers of constants; one which reads a
    state as well is an assignment rule (`test_observers_readouts`).
    """
    observers = [
        Observer("x1", "2 dimensionless * HEIGHT", U.cm, "x1"),
        Observer("x2", "x1 + Carm_drg * 0 cm / 1 mM", U.cm, "x2"),
        Observer("x3", "x1 / 2 dimensionless", U.cm, "x3"),
    ]
    model = compose(definition(observers=observers))
    rules = {str(r.variable) for r in model.rules}
    assignments = {str(a.symbol) for a in model.assignments}
    assert {"x1", "x3"} <= assignments
    assert "x2" in rules
    assert not {"x1", "x3"} & rules


def test_observer_cycle_is_rejected() -> None:
    """Observers which read each other in a cycle are rejected, naming them."""
    observers = [
        Observer("a1", "a2 + Carm_drg", U.mM, "a1"),
        Observer("a2", "a1 + Carm_drg", U.mM, "a2"),
    ]
    with pytest.raises(ValueError, match=r"\['a1', 'a2'\] read each other"):
        definition(observers=observers).validate()


def test_compose_checks_the_observers() -> None:
    """`compose` makes the checks of the observers on the model it composes."""
    from pkpdsim.pbpk import compose

    observer = Observer("x", "Cnope + 1 dimensionless", U.mM, "x")
    with pytest.raises(ValueError, match="Cnope"):
        compose(definition(observers=[observer]))
