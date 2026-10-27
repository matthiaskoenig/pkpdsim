"""Equivalence and analytic tests of the standard organ.

The forms of the organ are limits of each other: the merged tissue with a
large exchange is the flow-limited organ, the full form with large endothelial
and cell exchanges is the merged tissue, kinetic erythrocytes with a fast
exchange are erythrocytes in equilibrium. The equivalence tests compare the
venous plasma curves of an IV bolus (with a linear liver metabolism, so that
the curves decay) after the first 5 min, in which the bolus mixes. A limit is
approached, not reached: the deviation is first order in the reciprocal
exchange, and each test states the exchange which brings it below its
tolerance and why.

The analytic tests run a closed system into its equilibrium or a constant IV
infusion into its steady state and compare with the closed formulas of the
volume of distribution, the well-stirred liver and the extended clearance
model. The expected values are computed from the parameters of the model read
from roadrunner (volumes, fractions, flows, `HCT`, `fu`, `BP`, `PS_vas`); the
formulas are exact for the model, so the tolerances are those of the
integrator.
"""

import dataclasses
import itertools
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pytest

from pkpdsim.pbpk import (
    Direction,
    Erythrocytes,
    Excretion,
    FirstOrder,
    Kp,
    Membrane,
    Metabolism,
    Mode,
    Organ,
    Part,
    Pool,
    Route,
    Transport,
    ids,
)
from pkpdsim.pbpk.builder import build
from pkpdsim.pbpk.distribution import K_VAS

from .conftest import definition, individual, simple_drug, total_amount

if TYPE_CHECKING:
    import roadrunner

S = "drg"
"""The substance of the tests, the simple drug of the test helpers (fu 0.1)."""

CLINT = 5.0
"""The unbound intrinsic clearance of the liver metabolism (l/min)."""

DAY = 24 * 60.0
"""One day (min)."""

STEADY = 10 * DAY
"""The end of the equilibrations (min): the slowest time constant of the tests,
the elimination from the volume of distribution, is about 5 h, so the
equilibrium or steady state is reached far below the tolerances."""

TIMES = np.linspace(0.0, DAY, 2 * 1440 + 1)
"""The time points of the curves (min)."""

MIXED = TIMES >= 5.0
"""The time points of the comparison, after the mixing of the bolus."""

INFUSION = 1.0
"""The infusion rate of the steady state tests (mg/min)."""

KP_ORGANS = {Organ.LIVER: 5.0, Organ.ADIPOSE: 0.5}
"""The coefficients of single organs in the volume of distribution tests."""

PS_ACT = 20.0
"""The active uptake (OATP) on the unbound interstitial concentration (l/min)."""

PS_DIF = 1.0
"""The passive exchange of the liver cells (l/min), the `ps` of the liver."""


def _metabolism() -> Metabolism:
    """The linear liver metabolism `CLINT` on the unbound concentration."""
    return Metabolism(
        sid="DRGX",
        organ=Organ.LIVER,
        substrate=S,
        products={},
        kinetics=FirstOrder(k=CLINT, per_volume=False),
    )


def _model(
    tmp_path: Path,
    sid: str,
    chemistry: list[Metabolism | Transport],
    **changes: object,
) -> Path:
    """Build the drug without excretion, given intravenously.

    Args:
        tmp_path: the folder of the build.
        sid: the id of the model.
        chemistry: the processes.
        **changes: fields of the drug which replace those of `simple_drug`.

    Returns:
        The path of the flat model.
    """
    drug = dataclasses.replace(simple_drug(), excretion=Excretion(), **changes)
    d = definition(
        sid=sid, substances=[drug], chemistry=chemistry, routes=[Route.iv(S)]
    )
    return build(d, tmp_path).flat_path


def _load(path: Path, **changes: float) -> "roadrunner.RoadRunner":
    """Load a model with changes before its initialization, tight tolerances."""
    r = individual(path, **changes)
    r.integrator.absolute_tolerance = 1e-12
    r.integrator.relative_tolerance = 1e-10
    return r


def _curve(r: "roadrunner.RoadRunner") -> np.ndarray:
    """The venous plasma concentration (mM) after an IV bolus of 10 mg."""
    r[ids.iv_dose(S)] = 10.0
    cve = f"[{ids.pool_species(Pool.VENOUS, S)}]"
    return r.simulate(times=TIMES, selections=["time", cve])[:, 1]


def _deviation(a: np.ndarray, b: np.ndarray) -> float:
    """The largest relative deviation of curve `b` from curve `a` after mixing."""
    return float(np.max(np.abs(b[MIXED] - a[MIXED]) / a[MIXED]))


def _liver(local: str) -> str:
    """The flat id of a local id of the liver."""
    return ids.flat(Organ.LIVER, local)


def test_merged_large_ps_equals_perfusion(tmp_path: Path) -> None:
    """The merged tissue converges to the flow-limited organ with `ftissue`.

    The merged tissue lags behind the plasma by its exchange
    `PS_o*(Cu_pl - Cu_ti)`, the flow-limited organ is the limit of an infinite
    exchange. The liver metabolism acts on the unbound tissue concentration,
    behind the exchange of the liver, so the unbound clearance from plasma is
    `1/(1/PS_li + 1/CLint)` instead of `CLint`; `PS_li = ftissue *
    Vti_li/Vti_total_ref` is a few hundred l/min at `ftissue` 1e4 l/min, and
    the terminal phase deviates by about 3 % after a day. The deviation is
    first order in `1/ftissue` (a factor 10 per decade), `ftissue` 1e6 l/min
    brings it below 1e-3.
    """
    perfusion = _curve(
        _load(
            _model(
                tmp_path,
                "perfusion",
                [_metabolism()],
                distribution=Kp(2.0, mode=Mode.PERFUSION),
            )
        )
    )
    merged = _model(
        tmp_path, "merged", [_metabolism()], distribution=Kp(2.0, ftissue=1.0)
    )
    deviations = [
        _deviation(perfusion, _curve(_load(merged, **{ids.ftissue(S): ftissue})))
        for ftissue in (1e4, 1e5, 1e6)
    ]
    assert deviations[-1] < 1e-3
    for coarse, fine in itertools.pairwise(deviations):
        assert coarse / fine == pytest.approx(10.0, rel=0.2)


def test_full_large_ps_equals_merged(tmp_path: Path) -> None:
    """The full form with large exchanges equals the merged tissue.

    The full form has two exchanges where the merged tissue has one: the
    endothelium `PS_vas = k_vas*V` between plasma and interstitium and the
    cell membrane `PS_o` between interstitium and cells. The merged tissue is
    its limit only when both are large:

    - the endothelium is a resistance the merged tissue does not have; at the
      default `k_vas` 1000 1/min it shifts the curve by about 4e-3, in
      proportion to `1/k_vas`, so the test raises `k_vas` to 1e4 1/min in
      every organ and checks that the deviation falls with it;
    - at a finite cell exchange the interstitium is in equilibrium with the
      plasma in the full form (the endothelium is fast) but behind the
      exchange in the merged tissue: the cells have the time constant
      `Vic*K_cell/(PS*fu)`, the merged tissue `Vti*Kp/(PS*fu)`. Both forms
      therefore take the large exchange `ftissue` 1e4 l/min; what remains of
      its finite value (the deviation from the flow-limited organ, above) is
      the same in both forms and cancels.
    """
    merged = _curve(
        _load(
            _model(
                tmp_path, "merged", [_metabolism()], distribution=Kp(2.0, ftissue=1e4)
            )
        )
    )
    full = _model(
        tmp_path,
        "full",
        [_metabolism()],
        distribution=Kp(2.0, mode=Mode.PERMEABILITY, ftissue=1e4),
    )
    # the intestine (submodel GU) replaces the large intestine
    fast = {ids.flat(o, K_VAS): 1e4 for o in Organ if o is not Organ.LARGE_INTESTINE}
    deviation = _deviation(merged, _curve(_load(full, **fast)))
    assert deviation < 1e-3
    assert _deviation(merged, _curve(_load(full))) > 5 * deviation


def test_erythrocytes_fast_kinetic_equals_equilibrium(tmp_path: Path) -> None:
    """Kinetic erythrocytes with a fast exchange equal erythrocytes in equilibrium.

    The blood cells approach `Kbc*Cpl` with the rate `k_exchange/Kbc`
    (5e3 1/min), fast against every flow; the flows carry
    `Q*((1-HCT)*Cpl + HCT*Cbc)`, which is `Q*BP*Cpl` in equilibrium with
    `BP = 1 - HCT + HCT*Kbc`.
    """
    equilibrium, kinetic = (
        _curve(
            _load(
                _model(
                    tmp_path,
                    sid,
                    [_metabolism()],
                    distribution=Kp(2.0, ftissue=1.0),
                    erythrocytes=ery,
                )
            )
        )
        for sid, ery in (
            ("equilibrium", Erythrocytes("equilibrium", kbc=2.0)),
            ("kinetic", Erythrocytes("kinetic", kbc=2.0, k_exchange=1e4)),
        )
    )
    assert _deviation(equilibrium, kinetic) < 1e-3


def _vss(r: "roadrunner.RoadRunner", kbc: float) -> float:
    """The volume of distribution (l) from the parameters of the model.

    `sum(Vpl + Vbc*Kbc + Vti*Kp_o)` over the organs, with
    `Vpl = V*f_vas*(1-HCT)`, `Vbc = V*f_vas*HCT` and `Vti = V*(1-f_vas)`,
    and `Vpl + Vbc*Kbc` over the pools.
    """
    hct = r[ids.HCT]
    blood = 1 - hct + hct * kbc
    vss = sum(r[ids.pool_volume(p)] * blood for p in Pool)
    for o in Organ:
        v, f_vas = r[ids.organ_volume(o)], r[ids.vascular_fraction(o)]
        kp = r[ids.kp_organ(o, S)] if o in KP_ORGANS else r[ids.kp(S)]
        vss += v * f_vas * blood + v * (1 - f_vas) * kp
    return vss


def _closed_vss(r: "roadrunner.RoadRunner") -> float:
    """Amount / venous plasma concentration (l) at the equilibrium of a closed system.

    An IV infusion of 10 mg over 10 min without elimination, then the
    distribution into every organ to the equilibrium.
    """
    r[ids.infusion_rate(S)] = 1.0
    r.simulate(0, 10, 2)
    r[ids.infusion_rate(S)] = 0.0
    r.simulate(10, STEADY, 2)
    return total_amount(r, S) / r[f"[{ids.pool_species(Pool.VENOUS, S)}]"]


@pytest.mark.parametrize(
    ("kbc", "ery"),
    [(0.0, None), (3.0, Erythrocytes("equilibrium", kbc=3.0))],
    ids=["no erythrocytes", "equilibrium"],
)
def test_vss_perfusion(tmp_path: Path, kbc: float, ery: Erythrocytes | None) -> None:
    """At the equilibrium of a closed system amount/Cve is the sum of the capacities.

    Every flow-limited organ holds `V*K_or*Cpl` with
    `K_or = (Vpl + Vbc*Kbc + Vti*Kp)/V`, every pool `V*BP*Cpl`; without
    erythrocyte distribution `Kbc` is 0 (`BP = 1 - HCT`).
    """
    path = _model(
        tmp_path,
        "vss",
        [],
        distribution=Kp(2.0, mode=Mode.PERFUSION, overrides=KP_ORGANS),
        erythrocytes=ery,
    )
    r = _load(path)
    assert _closed_vss(r) == pytest.approx(_vss(r, kbc), rel=1e-9)


@pytest.mark.parametrize("mode", [Mode.TISSUE, Mode.PERMEABILITY])
def test_vss_exchange(tmp_path: Path, mode: Mode) -> None:
    """The merged and the full form have the volume of distribution of the Kp.

    The full form splits the Kp into `K_int_pls` and `K_cell_pls`, kinetic
    erythrocytes hold `Vbc*Kbc*Cpl` in equilibrium. `ftissue` 100 l/min
    brings the slowest organ (adipose) to the equilibrium within hours.
    """
    path = _model(
        tmp_path,
        "vss",
        [],
        distribution=Kp(2.0, mode=mode, ftissue=100.0, overrides=KP_ORGANS),
        erythrocytes=Erythrocytes("kinetic", kbc=3.0, k_exchange=1.0),
    )
    r = _load(path)
    assert _closed_vss(r) == pytest.approx(_vss(r, 3.0), rel=1e-9)


def _steady_state(r: "roadrunner.RoadRunner") -> tuple[float, float]:
    """The elimination rate (mmole/min) and `Car` (mM) at the steady state of the infusion.

    The steady state is checked: the liver metabolism eliminates what the
    infusion gives.
    """
    r[ids.infusion_rate(S)] = INFUSION
    r.simulate(0, STEADY, 2)
    rate = INFUSION / r[ids.mr(S)]
    assert r[_liver("DRGX")] == pytest.approx(rate, rel=1e-9)
    return rate, r[f"[{ids.pool_species(Pool.ARTERIAL, S)}]"]


def _hepatic(r: "roadrunner.RoadRunner", clu: float) -> float:
    """The plasma clearance (l/min) of the liver with the unbound clearance `clu`.

    `Qh*fu*CLu/(Qh + fu*CLu/BP)`, the well-stirred liver.
    """
    qh, fu, bp = r[ids.QH], r[ids.fu(S)], r[ids.bp(S)]
    return qh * fu * clu / (qh + fu * clu / bp)


@pytest.mark.parametrize(
    "ery", [None, Erythrocytes("equilibrium", bp=1.3)], ids=["plasma", "bp"]
)
def test_well_stirred_hepatic_clearance(
    tmp_path: Path, ery: Erythrocytes | None
) -> None:
    """The clearance of a flow-limited liver is the well-stirred clearance.

    The flows carry `Q*BP*C` of the plasma concentration `C`: without
    erythrocyte distribution `BP = 1 - HCT`, the blood cells hold no drug and
    the plasma flow `Q*(1-HCT)` carries it; with erythrocytes in equilibrium
    the blood `Q` carries `BP*C`. The metabolism acts on the unbound
    concentration of the blood leaving the liver, `fu*Cout`. At the steady
    state the portal organs eliminate nothing, the portal vein has the
    arterial concentration `Car`, and the balance of the liver

    `Qh*BP*Car = Qh*BP*Cout + fu*CLint*Cout`

    gives the plasma clearance

    `rate/Car = fu*CLint*Cout/Car = Qh*fu*CLint/(Qh + fu*CLint/BP)`

    and the blood clearance `rate/(BP*Car) = Qh*fub*CLint/(Qh + fub*CLint)`
    with the unbound fraction in blood `fub = fu/BP`. The extraction is
    0.33 (plasma) and 0.18 (`BP` 1.3), so both terms of the denominator
    count.
    """
    path = _model(
        tmp_path,
        "well_stirred",
        [_metabolism()],
        distribution=Kp(2.0, mode=Mode.PERFUSION),
        erythrocytes=ery,
    )
    r = _load(path)
    rate, car = _steady_state(r)
    assert rate / car == pytest.approx(_hepatic(r, CLINT), rel=1e-8)
    qh, fub, bp = r[ids.QH], r[ids.fu(S)] / r[ids.bp(S)], r[ids.bp(S)]
    assert rate / (bp * car) == pytest.approx(
        qh * fub * CLINT / (qh + fub * CLINT), rel=1e-8
    )


def _unbound(r: "roadrunner.RoadRunner", part: Part, k: str) -> float:
    """The unbound concentration `fu/K*C` (mM) of a part of the liver."""
    c = r[f"[{_liver(ids.species(part, S))}]"]
    return r[ids.fu(S)] / r[_liver(k)] * c


def test_extended_clearance(tmp_path: Path) -> None:
    """The clearance of a liver with uptake is the extended clearance.

    The uptake makes the liver the full form: active influx `PS_act*Cu_is`,
    passive exchange `PS_dif*(Cu_is - Cu_ic)` (the `ps` of the liver) and the
    metabolism `CLint*Cu_ic`, all on unbound concentrations. The steady state
    of the cells gives the liver Kp,uu

    `Cu_ic/Cu_is = (PS_act + PS_dif)/(PS_dif + CLint)`

    (the passive exchange is an influx as well, the total influx is
    `PS_act + PS_dif`) and the unbound clearance from the interstitium
    `CLall = (PS_act + PS_dif)*CLint/(PS_dif + CLint)`. The endothelium
    `PS_vas*(Cu_pl - Cu_is)` is in series, so the unbound clearance from
    plasma is `CLu = 1/(1/PS_vas + 1/CLall)`, and the balance of the liver
    plasma, as in the well-stirred liver (`BP = 1 - HCT`), gives

    `rate/Car = Qh*fu*CLu/(Qh + fu*CLu/BP)`.

    The extended clearance model of the literature is the limit of a fast
    endothelium, `CLu = CLall` and `Cu_is = fu*Cpl`: at the default `k_vas`
    it is off by `CLall/PS_vas` (about 1 %, 0.3 % after the extraction), so
    it is checked with `k_vas` raised to 1e6 1/min. The passive influx
    matters: without it the clearance is 1.5 % lower.
    """
    uptake = Transport(
        sid="OATP",
        organ=Organ.LIVER,
        substance=S,
        direction=Direction.INFLUX,
        membrane=Membrane.BASOLATERAL,
        kinetics=FirstOrder(k=PS_ACT, per_volume=False),
    )
    path = _model(
        tmp_path,
        "extended",
        [uptake, _metabolism()],
        distribution=Kp(
            2.0, mode=Mode.PERFUSION, ps={Organ.LIVER: PS_DIF}, ps_per_volume=False
        ),
    )
    kp_uu = (PS_ACT + PS_DIF) / (PS_DIF + CLINT)
    clall = kp_uu * CLINT

    r = _load(path)
    rate, car = _steady_state(r)
    ps_vas = r[_liver(ids.PS_VAS)]
    assert rate / car == pytest.approx(
        _hepatic(r, 1 / (1 / ps_vas + 1 / clall)), rel=1e-8
    )
    cu_is = _unbound(r, Part.INTERSTITIAL, ids.k_int_pls(S))
    cu_ic = _unbound(r, Part.INTRACELLULAR, ids.k_cell_pls(S))
    assert cu_ic / cu_is == pytest.approx(kp_uu, rel=1e-8)

    r = _load(path, **{_liver(K_VAS): 1e6})
    rate, car = _steady_state(r)
    assert rate / car == pytest.approx(_hepatic(r, clall), rel=1e-3)
    cu_pl = r[ids.fu(S)] * r[f"[{_liver(ids.species(Part.PLASMA, S))}]"]
    cu_ic = _unbound(r, Part.INTRACELLULAR, ids.k_cell_pls(S))
    assert cu_ic / cu_pl == pytest.approx(kp_uu, rel=1e-3)


@pytest.mark.parametrize("mode", [Mode.PERFUSION, Mode.TISSUE, Mode.PERMEABILITY])
def test_vss_with_gut_mucosa(tmp_path: Path, mode: Mode) -> None:
    """Gut chemistry adds mucosa species but not capacity: the mucosa is carved out of the wall.

    A CYP3A4 metabolism of the small intestine with its activity 0 (a closed
    system) gives the drug the mucosa of every intestinal segment; wall and
    mucosa together hold what the wall holds without gut chemistry, so the
    volume of distribution is the sum of the capacities of the organs, as
    without (and the same value).
    """
    gut = Metabolism(
        sid="DRG_GU",
        organ=Organ.SMALL_INTESTINE,
        substrate=S,
        products={},
        kinetics=FirstOrder(k=1.0),
        protein="CYP3A4",
    )
    distribution = Kp(2.0, mode=mode, ftissue=100.0, overrides=KP_ORGANS)
    vss = {}
    cases: list[list[Metabolism | Transport]] = [[], [gut]]
    for chemistry in cases:
        path = _model(
            tmp_path, f"vss_gut{len(chemistry)}", chemistry, distribution=distribution
        )
        changes = {ids.flat(Organ.SMALL_INTESTINE, ids.activity("CYP3A4")): 0.0}
        r = _load(path, **(changes if chemistry else {}))
        mucosa = ids.flat(Organ.SMALL_INTESTINE, ids.mucosa(S, "du"))
        assert (mucosa in r.model.getFloatingSpeciesIds()) is bool(chemistry)
        vss[len(chemistry)] = _closed_vss(r)
        assert vss[len(chemistry)] == pytest.approx(_vss(r, 0.0), rel=1e-9)
    assert vss[1] == pytest.approx(vss[0], rel=1e-9)
