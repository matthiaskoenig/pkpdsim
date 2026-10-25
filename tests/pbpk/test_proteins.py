"""Tests of the activity of the proteins, their inhibition and induction.

A process of a protein runs with the base activity `f_<protein>` of its organ;
a reversible inhibitor enters the rate with its unbound concentration in the
compartment of the process (competitive on Km, uncompetitive on Vmax and Km,
noncompetitive on Vmax, mixed on both), a mechanism-based inhibitor or an
inducer changes the turnover state `E_<protein>`, which multiplies the rates.
"""

import dataclasses
import math
from pathlib import Path
from typing import TYPE_CHECKING

import libsbml
import numpy as np
import pytest

from pkpdsim.pbpk import (
    MM,
    Binding,
    Direction,
    FirstOrder,
    Hill,
    Induction,
    Inhibition,
    Kp,
    Membrane,
    Metabolism,
    Organ,
    Part,
    Route,
    Substance,
    Transport,
    ids,
)
from pkpdsim.pbpk.builder import BuildResult, build
from pkpdsim.pbpk.proteins import turnover

from .conftest import definition, individual, load, simple_drug

if TYPE_CHECKING:
    import roadrunner

S, INH = "drg", "inh"
"""The substrate and the inhibitor (or inducer) of the tests."""

VMAX, KM, KI = 1e-3, 0.01, 0.005
"""Vmax (mmole/min/l), Km and Ki (mM) of the tests."""

DAY = 24 * 60.0
"""One day (min)."""

PROCESS = ids.flat(Organ.LIVER, "DRG2M")
"""The flat id of the liver metabolism of the substrate."""


def _inhibitor(
    inhibits: tuple[Inhibition, ...] = (),
    induces: tuple[Induction, ...] = (),
    distributes: bool = True,
) -> Substance:
    """The inhibitor: fu 0.5, Kp 1.5 (merged tissue), not eliminated."""
    return Substance(
        sid=INH,
        name="inhibitor",
        mr=400.0,
        binding=Binding(fu=0.5),
        distribution=Kp(1.5, ftissue=1.0) if distributes else None,
        inhibits=inhibits,
        induces=induces,
    )


def _build(
    tmp_path: Path,
    inhibitor: Substance,
    kinetics: FirstOrder | MM | Hill | None = None,
    chemistry: list[Metabolism | Transport] | None = None,
) -> BuildResult:
    """Build the substrate with its CYP3A4 liver metabolism and the inhibitor."""
    if chemistry is None:
        chemistry = [
            Metabolism(
                sid="DRG2M",
                organ=Organ.LIVER,
                substrate=S,
                products={},
                kinetics=kinetics or MM(vmax=VMAX, km=KM),
                protein="CYP3A4",
            )
        ]
    d = definition(
        sid="ddi",
        substances=[simple_drug(), inhibitor],
        chemistry=chemistry,
        routes=[Route.iv(S), Route.iv(INH)],
    )
    return build(d, tmp_path)


def _liver(part: Part, s: str) -> str:
    """The flat id of the concentration of a species of the liver."""
    return f"[{ids.flat(Organ.LIVER, ids.species(part, s))}]"


def _state(r: "roadrunner.RoadRunner", c: float, i: float) -> tuple[float, float]:
    """Set the liver tissue of substrate and inhibitor, get their unbound concentrations."""
    r[_liver(Part.TISSUE, S)] = c
    r[_liver(Part.TISSUE, INH)] = i
    return r[ids.fu(S)] / r[ids.kp(S)] * c, r[ids.fu(INH)] / r[ids.kp(INH)] * i


def _rate_at(
    tmp_path: Path, inhibition: Inhibition, kinetics: FirstOrder | MM | Hill
) -> tuple[float, float, float, float]:
    """The rate of the metabolism at a fixed state.

    Returns:
        The rate (mmole/min), the liver volume, the unbound substrate and
        inhibitor concentrations.
    """
    r = load(_build(tmp_path, _inhibitor((inhibition,)), kinetics))
    cu, iu = _state(r, 0.02, 0.01)
    return r[PROCESS], r[ids.organ_volume(Organ.LIVER)], cu, iu


def test_competitive_inhibition_rate(tmp_path: Path) -> None:
    """Competitive: `Vmax*V*C/(Km*(1 + I/Ki) + C)` on the unbound concentrations."""
    rate, v, cu, iu = _rate_at(
        tmp_path, Inhibition("CYP3A4", "competitive", KI), MM(VMAX, KM)
    )
    assert iu > 0
    assert rate == pytest.approx(VMAX * v * cu / (KM * (1 + iu / KI) + cu), rel=1e-12)


def test_noncompetitive_inhibition_rate(tmp_path: Path) -> None:
    """Noncompetitive: `Vmax/(1 + I/Ki)`, Km unchanged."""
    rate, v, cu, iu = _rate_at(
        tmp_path, Inhibition("CYP3A4", "noncompetitive", KI), MM(VMAX, KM)
    )
    assert rate == pytest.approx(VMAX / (1 + iu / KI) * v * cu / (KM + cu), rel=1e-12)


def test_uncompetitive_inhibition_rate(tmp_path: Path) -> None:
    """Uncompetitive: `Vmax_app = Vmax/(1 + I/Ki)`, `Km_app = Km/(1 + I/Ki)`."""
    rate, v, cu, iu = _rate_at(
        tmp_path, Inhibition("CYP3A4", "uncompetitive", KI), MM(VMAX, KM)
    )
    f = 1 + iu / KI
    assert rate == pytest.approx(VMAX / f * v * cu / (KM / f + cu), rel=1e-12)


def test_mixed_inhibition_rate(tmp_path: Path) -> None:
    """Mixed: `Vmax_app = Vmax/(1 + I/Kiu)`, `Km_app = Km*(1 + I/Kic)/(1 + I/Kiu)`."""
    kiu = 0.02
    rate, v, cu, iu = _rate_at(
        tmp_path, Inhibition("CYP3A4", "mixed", KI, kiu=kiu), MM(VMAX, KM)
    )
    fc, fu = 1 + iu / KI, 1 + iu / kiu
    assert rate == pytest.approx(VMAX / fu * v * cu / (KM * fc / fu + cu), rel=1e-12)


def test_mixed_inhibition_hill_rate(tmp_path: Path) -> None:
    """A Hill process takes the apparent Vmax and Km of the inhibition."""
    kiu, n = 0.02, 1.7
    rate, v, cu, iu = _rate_at(
        tmp_path, Inhibition("CYP3A4", "mixed", KI, kiu=kiu), Hill(VMAX, KM, n)
    )
    fc, fu = 1 + iu / KI, 1 + iu / kiu
    km_app = KM * fc / fu
    assert rate == pytest.approx(VMAX / fu * v * cu**n / (km_app**n + cu**n), rel=1e-12)


@pytest.mark.parametrize(
    ("inhibition", "inhibited"),
    [
        (Inhibition("CYP3A4", "competitive", KI), True),
        (Inhibition("CYP3A4", "noncompetitive", KI), True),
        (Inhibition("CYP3A4", "uncompetitive", KI), False),
    ],
)
def test_first_order_inhibition_rate(
    tmp_path: Path, inhibition: Inhibition, inhibited: bool
) -> None:
    """A first order process (`C << Km`) is `k*C/(1 + I/Ki)`; uncompetitive leaves it."""
    k = 0.5
    rate, _, cu, iu = _rate_at(tmp_path, inhibition, FirstOrder(k))
    expected = k * cu / (1 + iu / KI) if inhibited else k * cu
    assert rate == pytest.approx(expected, rel=1e-12)


def test_inhibitors_add(tmp_path: Path) -> None:
    """Two competitive inhibitors add their terms, `Km*(1 + I1/Ki1 + I2/Ki2)`."""
    second = Substance(
        sid="inh2",
        name="second inhibitor",
        mr=200.0,
        inhibits=(Inhibition("CYP3A4", "competitive", 0.05),),
    )
    d = definition(
        sid="ddi2",
        substances=[
            simple_drug(),
            _inhibitor((Inhibition("CYP3A4", "competitive", KI),)),
            second,
        ],
        chemistry=[
            Metabolism(
                sid="DRG2M",
                organ=Organ.LIVER,
                substrate=S,
                products={},
                kinetics=MM(VMAX, KM),
                protein="CYP3A4",
            )
        ],
        routes=[Route.iv(S)],
    )
    r = load(build(d, tmp_path))
    cu, iu = _state(r, 0.02, 0.01)
    r[_liver(Part.PLASMA, "inh2")] = 0.03
    i2 = 0.03  # vascular only, fu 1: the plasma concentration
    v = r[ids.organ_volume(Organ.LIVER)]
    assert r[PROCESS] == pytest.approx(
        VMAX * v * cu / (KM * (1 + iu / KI + i2 / 0.05) + cu), rel=1e-12
    )


def _model(result: BuildResult) -> libsbml.Model:
    """The flat model of a build."""
    return libsbml.readSBMLFromFile(str(result.flat_path)).getModel()


def _kinetic_law(result: BuildResult, organ: Organ, rid: str) -> str:
    """The formula of the kinetic law of a reaction of the flat model."""
    model = _model(result)
    reaction = model.getReaction(ids.flat(organ, rid))
    return libsbml.formulaToL3String(reaction.getKineticLaw().getMath())


def test_inhibitor_without_tissue_acts_from_plasma(tmp_path: Path) -> None:
    """An inhibitor without a species in the part of the process acts from plasma."""
    inhibitor = _inhibitor(
        (Inhibition("CYP3A4", "competitive", KI),), distributes=False
    )
    law = _kinetic_law(_build(tmp_path, inhibitor), Organ.LIVER, "DRG2M")
    assert ids.flat(Organ.LIVER, ids.species(Part.PLASMA, INH)) in law
    assert ids.flat(Organ.LIVER, ids.species(Part.TISSUE, INH)) not in law


def test_transporter_inhibitor_acts_where_the_substrate_is(tmp_path: Path) -> None:
    """Basolateral influx: the inhibitor of the interstitium (its merged tissue);
    apical influx: the inhibitor of the lumen of the segment."""
    chemistry: list[Metabolism | Transport] = [
        Transport(
            sid="OATP_DRG",
            organ=Organ.LIVER,
            substance=S,
            protein="SLCO1B1",
            direction=Direction.INFLUX,
            membrane=Membrane.BASOLATERAL,
            kinetics=MM(VMAX, KM),
        ),
        Transport(
            sid="OCT1_DRG",
            organ=Organ.SMALL_INTESTINE,
            substance=S,
            protein="SLC22A1",
            direction=Direction.INFLUX,
            membrane=Membrane.APICAL,
            kinetics=MM(VMAX, KM),
        ),
    ]
    inhibitor = _inhibitor(
        (
            Inhibition("SLCO1B1", "competitive", KI),
            Inhibition("SLC22A1", "competitive", KI),
        )
    )
    path = _build(tmp_path, inhibitor, chemistry=chemistry)
    liver = _kinetic_law(path, Organ.LIVER, "OATP_DRG")
    assert ids.flat(Organ.LIVER, ids.species(Part.TISSUE, INH)) in liver
    gut = _kinetic_law(
        path, Organ.SMALL_INTESTINE, ids.segment_process("OCT1_DRG", "du")
    )
    assert ids.flat(Organ.SMALL_INTESTINE, ids.lumen(INH, "du")) in gut


def test_activity_base_factor(tmp_path: Path) -> None:
    """`LI__f_cyp3a4 = 0.5` halves the clearance of a linear process.

    The base activity is the constant of the organ which genotype,
    impairment and the experiments set (`f_<protein>`, flat `LI__f_cyp3a4`).
    """
    k = 0.5
    path = _build(tmp_path, _inhibitor(), FirstOrder(k)).flat_path
    f = ids.flat(Organ.LIVER, ids.activity("CYP3A4"))
    clearance = {}
    for activity in (1.0, 0.5):
        r = individual(path, **{f: activity})
        r[ids.iv_dose(S)] = 10.0
        r.simulate(0, 60, 2)
        cu = r[ids.flat(Organ.LIVER, ids.unbound(Part.TISSUE, S))]
        clearance[activity] = r[PROCESS] / cu
    assert clearance[1.0] == pytest.approx(k, rel=1e-10)
    assert clearance[0.5] == pytest.approx(0.5 * k, rel=1e-10)


def _turnover_run(
    tmp_path: Path, inhibitor: Substance, days: float = 30.0
) -> tuple["roadrunner.RoadRunner", np.ndarray]:
    """Dose the (not eliminated) inhibitor and follow the liver CYP3A4 state."""
    r = load(_build(tmp_path, inhibitor))
    r[ids.iv_dose(INH)] = 100.0
    e = ids.flat(Organ.LIVER, ids.enzyme("CYP3A4"))
    r.timeCourseSelections = ["time", e]
    result = r.simulate(0, days * DAY, 301)
    return r, np.asarray(result[e])


def test_mechanism_based_inactivation_reduces_enzyme(tmp_path: Path) -> None:
    """MBI lowers `E` to `kdeg/(kdeg + kinact*INH/(KI + INH))` and inhibits competitively."""
    kinact = 1e-3
    inhibitor = _inhibitor(
        (Inhibition("CYP3A4", "mechanism", KI, kinact=kinact),),
    )
    r, e = _turnover_run(tmp_path, inhibitor)
    assert e[0] == 1.0
    assert np.all(np.diff(e) <= 1e-12)
    iu = r[ids.flat(Organ.LIVER, ids.unbound(Part.TISSUE, INH))]
    kdeg = r[ids.flat(Organ.LIVER, ids.kdeg("CYP3A4"))]
    expected = kdeg / (kdeg + kinact * iu / (KI + iu))
    assert 0.0 < expected < 0.9
    assert e[-1] == pytest.approx(expected, rel=1e-6)
    cu, iu = _state(r, 0.02, r[_liver(Part.TISSUE, INH)])
    v = r[ids.organ_volume(Organ.LIVER)]
    assert r[PROCESS] == pytest.approx(
        e[-1] * VMAX * v * cu / (KM * (1 + iu / KI) + cu), rel=1e-9
    )


def test_induction_increases_enzyme_to_steady_state(tmp_path: Path) -> None:
    """Induction raises `E` to `1 + Emax*INH/(EC50 + INH)`, which multiplies the rate."""
    emax, ec50 = 4.0, 0.002
    inhibitor = _inhibitor(induces=(Induction("CYP3A4", emax, ec50),))
    r, e = _turnover_run(tmp_path, inhibitor)
    assert np.all(np.diff(e) >= -1e-12)
    iu = r[ids.flat(Organ.LIVER, ids.unbound(Part.TISSUE, INH))]
    assert e[-1] == pytest.approx(1 + emax * iu / (ec50 + iu), rel=1e-6)
    cu, _ = _state(r, 0.02, r[_liver(Part.TISSUE, INH)])
    v = r[ids.organ_volume(Organ.LIVER)]
    assert r[PROCESS] == pytest.approx(e[-1] * VMAX * v * cu / (KM + cu), rel=1e-9)


def test_kdeg_of_liver_and_intestine(tmp_path: Path) -> None:
    """The turnover of the liver (36 h) and of the small intestine (23 h), cited."""
    chemistry: list[Metabolism | Transport] = [
        Metabolism(
            sid=f"DRG2M_{o.value}",
            organ=o,
            substrate=S,
            products={},
            kinetics=MM(VMAX, KM),
            protein="CYP3A4",
        )
        for o in (Organ.LIVER, Organ.SMALL_INTESTINE, Organ.KIDNEY)
    ]
    inhibitor = _inhibitor((Inhibition("CYP3A4", "mechanism", KI, kinact=1e-3),))
    model = _model(_build(tmp_path, inhibitor, chemistry=chemistry))
    half_life = {}
    for organ in (Organ.LIVER, Organ.SMALL_INTESTINE, Organ.KIDNEY):
        p = model.getParameter(ids.flat(organ, ids.kdeg("CYP3A4")))
        half_life[organ] = math.log(2) / p.getValue() / 60.0
        assert p.getValue() == pytest.approx(turnover(organ).value)
    assert half_life[Organ.LIVER] == pytest.approx(36.0, abs=0.5)
    assert half_life[Organ.SMALL_INTESTINE] == pytest.approx(23.0, abs=0.1)
    assert half_life[Organ.KIDNEY] == half_life[Organ.LIVER]
    notes = model.getParameter(
        ids.flat(Organ.LIVER, ids.kdeg("CYP3A4"))
    ).getNotesString()
    assert "10.1016/j.ejps.2011.04.008" in notes
    notes = model.getParameter(
        ids.flat(Organ.SMALL_INTESTINE, ids.kdeg("CYP3A4"))
    ).getNotesString()
    assert "10.1016/S0009-9236(03)00118-8" in notes


def test_interactions_only_where_a_process_uses_the_protein(tmp_path: Path) -> None:
    """Nothing is created for a protein in an organ without a process of it.

    The constants of the interactions are parameters of the body, of a
    protein which some process uses; the organs without a process of the
    protein have no activity, no turnover state and no inputs of it.
    """
    inhibitor = _inhibitor(
        (
            Inhibition("CYP3A4", "mechanism", KI, kinact=1e-3),
            Inhibition("CYP2C9", "competitive", KI),
        ),
        (Induction("CYP3A4", 2.0, 0.01),),
    )
    model = _model(_build(tmp_path, inhibitor))
    present = {p.getId() for p in model.getListOfParameters()}
    assert {
        ids.ki("CYP3A4", INH),
        ids.kinact("CYP3A4", INH),
        ids.emax("CYP3A4", INH),
        ids.ec50("CYP3A4", INH),
    } <= present
    assert ids.ki("CYP2C9", INH) not in present
    liver = {
        ids.flat(Organ.LIVER, x)
        for x in (ids.activity("CYP3A4"), ids.enzyme("CYP3A4"), ids.kdeg("CYP3A4"))
    }
    assert liver <= present
    gut = (ids.activity("CYP3A4"), ids.enzyme("CYP3A4"), ids.activity("CYP2C9"))
    assert not {ids.flat(Organ.SMALL_INTESTINE, x) for x in gut} & present
    assert model.getRateRule(ids.flat(Organ.LIVER, ids.enzyme("CYP3A4")))


def test_no_turnover_state_for_reversible_inhibition(tmp_path: Path) -> None:
    """A reversible inhibition adds no state: the rates carry its terms."""
    inhibitor = _inhibitor((Inhibition("CYP3A4", "noncompetitive", KI),))
    model = _model(_build(tmp_path, inhibitor))
    assert model.getParameter(ids.flat(Organ.LIVER, ids.enzyme("CYP3A4"))) is None
    assert not [
        r for r in model.getListOfRules() if r.isRate() and "cyp3a4" in r.getVariable()
    ]


@pytest.mark.parametrize(
    ("inhibition", "message"),
    [
        (Inhibition("CYP3A4", "mixed", KI), "kiu is required"),
        (Inhibition("CYP3A4", "competitive", KI, kiu=0.1), "kiu is required"),
        (Inhibition("CYP3A4", "mechanism", KI), "kinact is required"),
        (Inhibition("CYP3A4", "competitive", -1.0), "ki must be positive"),
    ],
)
def test_invalid_inhibition(inhibition: Inhibition, message: str) -> None:
    """An inhibition needs exactly the constants of its kind, positive."""
    d = definition(substances=[simple_drug(), _inhibitor((inhibition,))])
    with pytest.raises(ValueError, match=message):
        d.validate()


def test_one_inhibition_per_protein() -> None:
    """A substance inhibits a protein once (MBI includes the competitive part)."""
    twice = (
        Inhibition("CYP3A4", "competitive", KI),
        Inhibition("CYP3A4", "mechanism", KI, kinact=1e-3),
    )
    d = definition(substances=[simple_drug(), _inhibitor(twice)])
    with pytest.raises(ValueError, match="several inhibitions"):
        d.validate()


def test_autoinhibition(tmp_path: Path) -> None:
    """A substrate which inactivates its own enzyme lowers its clearance over time.

    The mechanism-based inhibitor of CYP3A4 is the substrate of the CYP3A4
    metabolism (time-dependent autoinhibition); its species is the reactant
    of the reaction, not repeated as modifier.
    """
    drug = dataclasses.replace(
        simple_drug(),
        inhibits=(Inhibition("CYP3A4", "mechanism", KI, kinact=1e-2),),
    )
    d = definition(
        sid="auto",
        substances=[drug, _inhibitor()],
        chemistry=[
            Metabolism(
                sid="DRG2M",
                organ=Organ.LIVER,
                substrate=S,
                products={},
                kinetics=MM(VMAX, KM),
                protein="CYP3A4",
            )
        ],
        routes=[Route.iv(S)],
    )
    result = build(d, tmp_path)
    reaction = _model(result).getReaction(PROCESS)
    assert reaction.getNumModifiers() == 0
    r = load(result)
    r[ids.iv_dose(S)] = 100.0
    e = ids.flat(Organ.LIVER, ids.enzyme("CYP3A4"))
    r.timeCourseSelections = ["time", e]
    values = np.asarray(r.simulate(0, DAY, 50)[e])
    assert values[0] == 1.0
    assert values[-1] < 0.99
