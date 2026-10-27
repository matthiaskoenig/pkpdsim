"""Tests of the specifics of the standard kidney and liver and of the bile.

The kidney filters the unbound plasma concentration into the urine,
`fu*GFR*Cpl*(1 - f_reabsorbed)`, and the plasma leaving it carries the plasma
flow less the filtration, so the filtration clearance of the arterial plasma
is `fu*GFR*(1 - f_reabsorbed)`; the liver has the portosystemic shunt
`f_shunt`, which bypasses the sinusoids, and the functional fraction
`f_functional`, which scales its tissue and the capacity of its processes;
the apical efflux of the liver writes into the bile, which drains into the
intestinal lumen.
"""

import dataclasses
from pathlib import Path
from typing import TYPE_CHECKING

import libsbml
import numpy as np
import pytest

from pkpdsim.pbpk import (
    Binding,
    Direction,
    Erythrocytes,
    Excretion,
    FirstOrder,
    Kp,
    Membrane,
    Metabolism,
    Organ,
    Part,
    PBPKDefinition,
    Pool,
    Route,
    Substance,
    Transport,
    ids,
)
from pkpdsim.pbpk.builder import build

from .conftest import definition, dose_amount, individual, simple_drug, total_amount

if TYPE_CHECKING:
    import roadrunner

S = "drg"
"""The substance of the tests."""

DAY = 24 * 60.0
"""One day (min)."""

INFUSION = 1.0
"""The infusion rate of the steady state tests (mg/min)."""


def _renal_drug(fraction_reabsorbed: float = 0.0, fu: float = 0.5) -> Substance:
    """A vascular drug (fu 0.5) eliminated by filtration only.

    Without tissue distribution the volume is the plasma, so the steady state
    of an infusion is reached within hours.
    """
    return Substance(
        sid=S,
        name="drug",
        mr=300.0,
        binding=Binding(fu=fu, protein="albumin"),
        excretion=Excretion(urine=True, fraction_reabsorbed=fraction_reabsorbed),
    )


def _renal(tmp_path: Path, fraction_reabsorbed: float = 0.0, fu: float = 0.5) -> Path:
    """Build the renal drug, given intravenously, without chemistry."""
    d = definition(
        sid="renal",
        substances=[_renal_drug(fraction_reabsorbed, fu)],
        chemistry=[],
        routes=[Route.iv(S)],
    )
    return build(d, tmp_path).flat_path


def _load(path: Path, **changes: float) -> "roadrunner.RoadRunner":
    """Load a model with changes before its initialization, tight tolerances."""
    r = individual(path, **changes)
    r.integrator.absolute_tolerance = 1e-12
    r.integrator.relative_tolerance = 1e-10
    return r


def _kidney_plasma(r: "roadrunner.RoadRunner") -> float:
    """The plasma concentration of the kidney (mM), `KI__Cpl_<s>` (R1)."""
    return r[f"[{ids.flat(Organ.KIDNEY, ids.species(Part.PLASMA, S))}]"]


def _filtration(r: "roadrunner.RoadRunner") -> float:
    """The rate of the filtration into the urine (mmole/min)."""
    return r[ids.flat(Organ.KIDNEY, ids.filtration_reaction(S))]


def test_renal_clearance_fu_gfr(tmp_path: Path) -> None:
    """At the steady state of an infusion the urine rate over `KI__Cpl` is `fu*GFR`.

    The filtration is the only elimination: the urine rate is the infusion.
    """
    r = _load(_renal(tmp_path))
    r[ids.infusion_rate(S)] = INFUSION
    r.simulate(0, 2 * DAY, 2)
    rate = _filtration(r)
    assert rate == pytest.approx(INFUSION / r[ids.mr(S)], rel=1e-8)
    assert rate / _kidney_plasma(r) == pytest.approx(
        r[ids.fu(S)] * r[ids.GFR], rel=1e-10
    )
    assert r[ids.urine(S)] > 0


@pytest.mark.parametrize(
    ("fu", "fraction_reabsorbed"), [(1.0, 0.0), (0.2, 0.3)], ids=["free", "bound"]
)
def test_filtration_clearance_of_the_arterial_plasma(
    tmp_path: Path, fu: float, fraction_reabsorbed: float
) -> None:
    """The systemic filtration clearance is `fu*GFR*(1 - f_reabsorbed)`.

    The filtered amount leaves the plasma of the kidney with the filtrate, the
    plasma leaving the kidney carries the rest: at the steady state of an
    infusion the plasma of the kidney is the arterial plasma, and the urine
    rate over the arterial (and the arm vein) plasma is the clearance, not
    `fu*GFR*Qp/(Qp + fu*GFR)` of a kidney which returns the whole plasma flow.
    """
    r = _load(_renal(tmp_path, fraction_reabsorbed, fu))
    r[ids.infusion_rate(S)] = INFUSION
    r.simulate(0, 4 * DAY, 2)
    rate = r[ids.urine_rate(S)]
    assert rate == pytest.approx(INFUSION / r[ids.mr(S)], rel=1e-8)
    clearance = fu * r[ids.GFR] * (1 - fraction_reabsorbed)
    arterial = r[f"[{ids.pool_species(Pool.ARTERIAL, S)}]"]
    assert rate / arterial == pytest.approx(clearance, rel=1e-6)
    assert rate / r[ids.arm(S)] == pytest.approx(clearance, rel=1e-6)
    assert _kidney_plasma(r) == pytest.approx(arterial, rel=1e-6)


def test_reabsorption_fraction(tmp_path: Path) -> None:
    """The reabsorbed fraction of the filtrate stays in the plasma.

    `f_reabsorbed_<s>` is a parameter of the top model, which a model fits.
    """
    path = _renal(tmp_path, fraction_reabsorbed=0.4)
    r = _load(path)
    assert r[ids.fraction_reabsorbed(S)] == 0.4
    r[ids.iv_dose(S)] = 10.0
    r.simulate(0, 60, 2)
    clearance = r[ids.fu(S)] * r[ids.GFR] * 0.6
    assert _filtration(r) / _kidney_plasma(r) == pytest.approx(clearance, rel=1e-10)
    r = _load(path, **{ids.fraction_reabsorbed(S): 1.0})
    r[ids.iv_dose(S)] = 10.0
    r.simulate(0, 60, 2)
    assert r[ids.urine(S)] == 0.0


def test_no_filtration_without_urine_excretion(tmp_path: Path) -> None:
    """A substance without `Excretion(urine=True)` is not filtered."""
    d = definition(sid="no_urine")
    result = build(d, tmp_path)
    model = libsbml.readSBMLFromFile(str(result.flat_path)).getModel()
    assert model.getReaction(ids.flat(Organ.KIDNEY, ids.filtration_reaction(S)))
    assert (
        model.getReaction(ids.flat(Organ.KIDNEY, ids.filtration_reaction("drg_m")))
        is None
    )
    assert model.getParameter(ids.fraction_reabsorbed(S)) is not None
    assert model.getParameter(ids.fraction_reabsorbed("drg_m")) is None


def _metabolized(tmp_path: Path, **changes: object) -> Path:
    """The simple drug with the liver metabolism, without excretion, IV.

    Args:
        tmp_path: the folder of the build.
        **changes: fields of the drug which replace those of `simple_drug`.
    """
    drug = dataclasses.replace(simple_drug(), excretion=Excretion(), **changes)
    metabolite = Substance(sid="drg_m", name="metabolite", mr=316.0)
    d = definition(sid="liver", substances=[drug, metabolite], routes=[Route.iv(S)])
    return build(d, tmp_path).flat_path


def test_shunt_bypasses_liver(tmp_path: Path) -> None:
    """With `f_shunt` 1 no blood reaches the liver: no hepatic metabolism.

    The shunt flows from the arterial and the portal pool into the venous
    pool, so the amount is conserved at every shunt fraction.
    """
    path = _metabolized(tmp_path)
    for f_shunt, metabolized in ((0.0, True), (0.5, True), (1.0, False)):
        r = _load(path, f_shunt=f_shunt)
        r[ids.iv_dose(S)] = 10.0
        r.simulate(0, DAY, 50)
        amount = total_amount(r, S) + total_amount(r, "drg_m") + dose_amount(r, S)
        assert amount == pytest.approx(10.0 / r[ids.mr(S)], rel=1e-7)
        if metabolized:
            assert total_amount(r, "drg_m") > 1e-3 * amount
        else:
            assert total_amount(r, "drg_m") == pytest.approx(0.0, abs=1e-15)
            assert r["[LI__Cti_drg]"] == pytest.approx(0.0, abs=1e-15)


def test_shunt_mass_balance_kinetic_erythrocytes(tmp_path: Path) -> None:
    """The shunt of the kinetic blood cells conserves the amount as well."""
    ery = Erythrocytes("kinetic", kbc=2.0, k_exchange=1.0)
    r = _load(_metabolized(tmp_path, erythrocytes=ery), f_shunt=0.5)
    r[ids.iv_dose(S)] = 10.0
    r.simulate(0, DAY, 50)
    amount = total_amount(r, S) + total_amount(r, "drg_m") + dose_amount(r, S)
    assert amount == pytest.approx(10.0 / r[ids.mr(S)], rel=1e-7)
    assert r[f"[{ids.pool_species(Pool.VENOUS, S, Part.BLOOD_CELLS)}]"] > 0
    assert (
        r[ids.flat(Organ.LIVER, ids.shunt_reaction(Pool.PORTAL, Part.BLOOD_CELLS, S))]
        > 0
    )


def test_functional_fraction_scales_capacity(tmp_path: Path) -> None:
    """`f_functional` scales the tissue volume and the Vmax share of the liver.

    The rate of the metabolism at the same unbound concentration halves with
    half of the liver functional, and so does the tissue volume; the
    vascular volume of the liver stays.
    """
    path = _metabolized(tmp_path)
    values = {}
    for f in (1.0, 0.5):
        r = _load(path, f_functional=f)
        r["[LI__Cti_drg]"] = 0.001
        values[f] = (r["LI__DRG2M"], r["LI__Vti"], r["LI__Vpl"], r["LI__PS_drg"])
    rate, vti, vpl, ps = values[0.5]
    assert rate == pytest.approx(0.5 * values[1.0][0], rel=1e-12)
    assert vti == pytest.approx(0.5 * values[1.0][1], rel=1e-12)
    assert vpl == values[1.0][2]
    assert ps == pytest.approx(0.5 * values[1.0][3], rel=1e-12)


def test_bile_reaches_lumen(tmp_path: Path) -> None:
    """The apical efflux of the liver writes into the bile, which drains into the duodenum.

    Without an oral route the substance is not reabsorbed (`ka_<s>` 0), so it
    goes on with the transit into the feces; the amount is conserved.
    """
    drug = dataclasses.replace(simple_drug(), excretion=Excretion())
    efflux = Transport(
        sid="DRGBILE",
        organ=Organ.LIVER,
        substance=S,
        direction=Direction.EFFLUX,
        membrane=Membrane.APICAL,
        kinetics=FirstOrder(k=0.5),
    )
    d = PBPKDefinition(
        sid="bile", substances=[drug], chemistry=[efflux], routes=[Route.iv(S)]
    )
    r = _load(build(d, tmp_path).flat_path)
    r[ids.iv_dose(S)] = 10.0
    result = r.simulate(
        0,
        DAY,
        49,
        selections=[
            "time",
            f"[{ids.bile(S)}]",
            f"[{ids.flat(Organ.SMALL_INTESTINE, ids.lumen(S, 'du'))}]",
        ],
    )
    assert np.max(result[:, 1]) > 0
    assert np.max(result[:, 2]) > 0
    assert r[ids.feces(S)] > 0
    assert total_amount(r, S) + dose_amount(r, S) == pytest.approx(
        10.0 / r[ids.mr(S)], rel=1e-7
    )


def test_hepatic_vein_observer(tmp_path: Path) -> None:
    """`Chv_<s>` is the flow-weighted mix of the liver outflow and the shunt."""
    path = _metabolized(tmp_path)
    for f_shunt in (0.0, 0.3):
        r = _load(path, f_shunt=f_shunt)
        r[ids.iv_dose(S)] = 10.0
        cols = [
            ids.hepatic_vein(S),
            ids.organ_outflow(Organ.LIVER, S),
            f"[{ids.pool_species(Pool.ARTERIAL, S)}]",
            f"[{ids.pool_species(Pool.PORTAL, S)}]",
        ]
        result = r.simulate(0, 120, 25, selections=["time", *cols])
        chv, cli, car, cpo = (result[:, i] for i in range(1, 5))
        qh, qha, qpo = r[ids.QH], r[ids.QHA], r[ids.QPO]
        expected = ((1 - f_shunt) * qh * cli + f_shunt * (qha * car + qpo * cpo)) / qh
        np.testing.assert_allclose(chv, expected, rtol=1e-12, atol=0)
        if f_shunt == 0.0:
            np.testing.assert_allclose(chv, cli, rtol=1e-12, atol=0)


# --- per-volume first order kinetics and the tubular capacity ---------------


def _kidney_transport(
    direction: Direction, membrane: Membrane, kinetics: FirstOrder
) -> Transport:
    """A transport `KITR` of the drug in the kidney."""
    return Transport("KITR", Organ.KIDNEY, S, direction, membrane, kinetics)


def _process_model(tmp_path: Path, process: Transport | Metabolism) -> Path:
    """Build the simple drug with one process, given intravenously."""
    d = definition(sid="process", chemistry=[process], routes=[Route.iv(S)])
    return build(d, tmp_path).flat_path


def _clearance(
    path: Path, process: str, species: str, coefficient: str, **changes: float
) -> float:
    """The clearance (l/min) of a process: its rate over its unbound substrate.

    The model is initialized with the changes (sbmlsim evaluates the initial
    assignments which read them again, R21), the substrate species is set
    and the rate read at that state; the unbound concentration is
    `fu*C/coefficient` with the values of the model (a CKD stage changes fu).

    Args:
        path: the flat model.
        process: the flat id of the reaction.
        species: the flat id of the substrate species.
        coefficient: the id of its coefficient to plasma (`K_int_pls_<s>`).
        **changes: the covariates and the impairment.
    """
    r = _load(path, **changes)
    r[f"[{species}]"] = 1e-3
    unbound = r[ids.fu(S)] * 1e-3 / r[coefficient]
    return r[process] / unbound


def _kidney_influx(path: Path, **changes: float) -> float:
    """The clearance of the basolateral influx of the kidney."""
    return _clearance(
        path,
        ids.flat(Organ.KIDNEY, "KITR"),
        ids.flat(Organ.KIDNEY, ids.species(Part.INTERSTITIAL, S)),
        ids.flat(Organ.KIDNEY, ids.k_int_pls(S)),
        **changes,
    )


def _kidney_efflux(path: Path, **changes: float) -> float:
    """The clearance of the apical efflux of the kidney (secretion)."""
    return _clearance(
        path,
        ids.flat(Organ.KIDNEY, "KITR"),
        ids.flat(Organ.KIDNEY, ids.species(Part.INTRACELLULAR, S)),
        ids.flat(Organ.KIDNEY, ids.k_cell_pls(S)),
        **changes,
    )


def _volume(path: Path, organ: Organ, **changes: float) -> float:
    """The volume of an organ of an individual (l)."""
    return _load(path, **changes)[ids.organ_volume(organ)]


def test_first_order_per_volume_scales_with_the_organ(tmp_path: Path) -> None:
    """A first order clearance per volume is `k*V<o>`: it follows the organ of the individual."""
    path = _process_model(
        tmp_path,
        _kidney_transport(
            Direction.INFLUX,
            Membrane.BASOLATERAL,
            FirstOrder(k=0.5),
        ),
    )
    r = _load(path)
    assert r[ids.flat(Organ.KIDNEY, "KITR_k")] == 0.5
    doc = libsbml.readSBMLFromFile(str(path))
    k = doc.getModel().getParameter(ids.flat(Organ.KIDNEY, "KITR_k"))
    assert k.getUnits().rsplit("__", 1)[-1] == "per_min"  # flat: `KI__per_min`
    assert _kidney_influx(path) == pytest.approx(
        0.5 * _volume(path, Organ.KIDNEY), rel=1e-9
    )
    ratio = _kidney_influx(path, BW=150.0) / _kidney_influx(path, BW=75.0)
    expected = _volume(path, Organ.KIDNEY, BW=150.0) / _volume(
        path, Organ.KIDNEY, BW=75.0
    )
    assert expected > 1.1
    assert ratio == pytest.approx(expected, rel=1e-9)


def test_first_order_whole_organ_does_not_scale(tmp_path: Path) -> None:
    """A first order clearance of the whole organ is `k`, whatever the organ volume."""
    path = _process_model(
        tmp_path,
        _kidney_transport(
            Direction.INFLUX, Membrane.BASOLATERAL, FirstOrder(k=0.5, per_volume=False)
        ),
    )
    assert _kidney_influx(path) == pytest.approx(0.5, rel=1e-9)
    ratio = _kidney_influx(path, BW=150.0) / _kidney_influx(path, BW=75.0)
    assert ratio == pytest.approx(1.0, rel=1e-9)


def _ps_model(tmp_path: Path, distribution: Kp) -> Path:
    """Build the simple drug with a distribution of its own, given intravenously."""
    drug = dataclasses.replace(simple_drug(), distribution=distribution)
    d = definition(sid="ps", substances=[drug], chemistry=[], routes=[Route.iv(S)])
    return build(d, tmp_path).flat_path


def test_ps_per_volume_scales_with_the_organ(tmp_path: Path) -> None:
    """The exchange of an organ given per litre of organ is `kps*V<o>`.

    It follows the organ of the individual when the body weight changes, as
    the share of `ftissue` and a first order clearance per volume do.
    """
    path = _ps_model(tmp_path, Kp(2.0, ftissue=1.0, ps={Organ.MUSCLE: 0.05}))
    ps, kps = ids.flat(Organ.MUSCLE, ids.ps(S)), ids.flat(Organ.MUSCLE, ids.kps(S))
    r = _load(path)
    assert r[kps] == 0.05
    doc = libsbml.readSBMLFromFile(str(path))
    unit = doc.getModel().getParameter(kps).getUnits().rsplit("__", 1)[-1]
    assert unit == "per_min"
    assert r[ps] == pytest.approx(0.05 * _volume(path, Organ.MUSCLE), rel=1e-9)
    ratio = _load(path, BW=150.0)[ps] / _load(path, BW=75.0)[ps]
    expected = _volume(path, Organ.MUSCLE, BW=150.0) / _volume(
        path, Organ.MUSCLE, BW=75.0
    )
    assert expected > 1.1
    assert ratio == pytest.approx(expected, rel=1e-9)


def test_ps_of_the_whole_organ_does_not_scale(tmp_path: Path) -> None:
    """With `ps_per_volume=False` the exchange of an organ is its value in l/min."""
    path = _ps_model(
        tmp_path,
        Kp(2.0, ftissue=1.0, ps={Organ.MUSCLE: 0.05}, ps_per_volume=False),
    )
    ps = ids.flat(Organ.MUSCLE, ids.ps(S))
    assert _load(path)[ps] == pytest.approx(0.05, rel=1e-9)
    assert _load(path, BW=150.0)[ps] == pytest.approx(0.05, rel=1e-9)
    assert (
        ids.flat(Organ.MUSCLE, ids.kps(S))
        not in _load(path).model.getGlobalParameterIds()
    )


def test_renal_function_scales_tubular_transport(tmp_path: Path) -> None:
    """The secretion of the kidney follows the renal function (intact nephron hypothesis).

    With half the renal function (a CKD stage, which changes fu) the
    clearance of the apical efflux on its unbound substrate halves.
    """
    path = _process_model(
        tmp_path,
        _kidney_transport(Direction.EFFLUX, Membrane.APICAL, FirstOrder(k=0.5)),
    )
    r = _load(path)
    assert r[ids.flat(Organ.KIDNEY, ids.F_TUBULAR)] == 1.0
    ratio = _kidney_efflux(path, f_renal_function=0.5) / _kidney_efflux(path)
    assert ratio == pytest.approx(0.5, rel=1e-9)


def test_renal_function_leaves_the_liver(tmp_path: Path) -> None:
    """The renal function does not change the capacity of the liver."""
    path = _process_model(
        tmp_path,
        Metabolism("DRG2M", Organ.LIVER, S, {}, FirstOrder(k=0.5, per_volume=True)),
    )

    def liver(**changes: float) -> float:
        return _clearance(
            path,
            ids.flat(Organ.LIVER, "DRG2M"),
            ids.flat(Organ.LIVER, ids.species(Part.TISSUE, S)),
            ids.kp(S),
            **changes,
        )

    assert liver(f_renal_function=0.5) / liver() == pytest.approx(1.0, rel=1e-9)
