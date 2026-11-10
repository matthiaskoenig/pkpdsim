"""The erythrocytes of the zones: their exchange with the plasma and their flow.

The haematocrit splits the lumen into two conduits. A substance which states
how it distributes into the erythrocytes
(`pkpdsim.pbpk.substance.Erythrocytes`, read as
`pkpdsim.zonation.transport.ErythrocyteDistribution`) exchanges between them,
in equilibrium (a coefficient, realized as an exchange which is fast against
the transit of a zone) or kinetically (its own rate); a substance which states
none keeps an erythrocyte species which exchanges with nothing.

The tests use a tracer with a tiny molecular radius, which no fenestra
hinders, so the convection of the two conduits and the exchange between them
set the profile and the axial diffusion is switched off.
"""

from pathlib import Path

import pytest
import roadrunner
from pbpk.conftest import individual
from sbmlutils.factory import create_model

from pkpdsim.pbpk.kp import blood_to_plasma
from pkpdsim.pbpk.substance import Erythrocytes, Substance
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Value
from pkpdsim.zonation import ZonationBuilder, ids, observers
from pkpdsim.zonation.builder import EQUILIBRATION
from pkpdsim.zonation.flow import ImposedFlow, PressureFlow
from pkpdsim.zonation.geometry import Layer
from pkpdsim.zonation.observers import add_observers, observer_objects
from pkpdsim.zonation.settings import Geometry, ZonationSettings
from pkpdsim.zonation.transport import Solute, blood_cell_coefficient

from .conftest import liver_flow, mass_g

HCT = 0.45
"""The haematocrit of the tests, the one of the rat."""

N = 10
"""The zones of the chain."""

EQUILIBRIUM = Erythrocytes(mode="equilibrium", kbc=5.0)
"""Erythrocytes which hold five times the concentration of the plasma."""

KINETIC = Erythrocytes(mode="kinetic", kbc=5.0, k_exchange=0.2)
"""The same coefficient, approached with a rate of the order of the transit of
the whole chain (0.2/min), so the erythrocytes fill along the zones."""


def _tracer(ery: Erythrocytes | None) -> Substance:
    """The tracer with an erythrocyte distribution, or without one."""
    return Substance("tr", "tracer", Value(180.0, U.g_per_mole), erythrocytes=ery)


class _Sized(ZonationBuilder):
    """A builder whose tracer is a point, which every fenestra lets through."""

    def _radius(self, substance: Solute) -> Value:
        """Get the radius of the tracer."""
        return Value(1e-12, U.meter)


def _settings(**kwargs: object) -> ZonationSettings:
    """The settings of a pressure-driven rat chain at the haematocrit of the tests."""
    kwargs.setdefault("flow", PressureFlow())
    kwargs.setdefault("haematocrit", HCT)
    kwargs.setdefault("geometry", Geometry.SINUSOID)
    return ZonationSettings(
        n_zones=N,
        bile=None,
        liver_mass=mass_g(7.229),
        reference_flow=liver_flow(),
        **kwargs,  # ty: ignore[invalid-argument-type]
    )


def _builder(ery: Erythrocytes | None, **kwargs: object) -> _Sized:
    """The builder of a chain of the tracer with an erythrocyte distribution."""
    return _Sized(_settings(**kwargs), [_tracer(ery)], None)


def _reactions(ery: Erythrocytes | None, **kwargs: object) -> set[str]:
    """The reaction ids of the built model."""
    return {str(r.sid) for r in _builder(ery, **kwargs).build().reactions}


def _parameters(ery: Erythrocytes | None, **kwargs: object) -> dict[str, float]:
    """The parameter values of the built model by id."""
    return {
        str(p.sid): float(p.value)  # ty: ignore[invalid-argument-type]
        for p in _builder(ery, **kwargs).build().parameters
    }


def _runner(
    tmp_path: Path,
    ery: Erythrocytes | None,
    *,
    observed: bool = False,
    **kwargs: object,
) -> roadrunner.RoadRunner:
    """A loaded chain of the tracer, the axial diffusion switched off."""
    builder = _builder(ery, **kwargs)
    model = builder.build()
    if observed:
        add_observers(model, observer_objects(builder.settings, [_tracer(ery)]))
    path = tmp_path / "erythrocytes.xml"
    create_model(model, filepath=path, sbml_level=3, sbml_version=2, validate=False)
    r = individual(path)
    r.integrator.absolute_tolerance = 1e-14
    r.integrator.relative_tolerance = 1e-12
    r[ids.diffusion("tr")] = 0.0
    return r


def _steady(
    tmp_path: Path, ery: Erythrocytes | None, **kwargs: object
) -> roadrunner.RoadRunner:
    """The steady state of a chain whose inflowing plasma holds 1 mM."""
    r = _runner(tmp_path, ery, **kwargs)  # ty: ignore[invalid-argument-type]
    r[ids.species("tr", Layer.PLASMA, ids.INFLOW)] = 1.0
    r.simulate(0.0, 400.0, 2)
    return r


def _profile(r: roadrunner.RoadRunner, layer: Layer | str) -> list[float]:
    """The concentrations of a layer over the zones."""
    return [r[f"[{ids.species('tr', layer, ids.zone(k))}]"] for k in range(1, N + 1)]


def _flows(r: roadrunner.RoadRunner, layer: Layer | None) -> list[float]:
    """The flow (l/min) of a layer across every face, the blood flow for `None`.

    The chain evaluates its flows at the faces of its zones, so the first
    entry is the flow the chain receives at its periportal face (`in`) and the
    others the flow which leaves every zone at its pericentral face: `N + 1`
    values.
    """
    sites = (ids.INFLOW, *(ids.zone(k) for k in range(1, N + 1)))
    return [r[ids.zone_flow(layer, site)] for site in sites]


# --- the flow of the two conduits -------------------------------------------


def test_the_erythrocytes_keep_their_flow_along_a_pressure_driven_sinusoid(
    tmp_path: Path,
) -> None:
    """An erythrocyte does not pass a fenestra, so its conduit loses no water.

    The flow of the erythrocytes is the share `Hct` of the blood entering the
    chain across every face, the plasma carries the rest of the flow of the
    leaky tube there, and the two plus the filtrate of the space of Disse add
    up to the inflow. The engine split the local flow of every zone into `Hct`
    and `1 - Hct` before, which let the erythrocyte conduit lose water to
    nowhere.
    """
    r = _steady(tmp_path, None)
    blood = _flows(r, None)
    rbc = _flows(r, Layer.RBC)
    plasma = _flows(r, Layer.PLASMA)
    assert blood[5] < blood[0]  # the tube leaks, so there is something to lose
    assert rbc == pytest.approx([HCT * blood[0]] * (N + 1), rel=1e-12)
    assert plasma == pytest.approx([q - HCT * blood[0] for q in blood], rel=1e-12)
    # the filtrate across the periportal face of the chain is 0, nothing having
    # filtered yet, and every zone filters what its lumen loses over its length
    filtrate = [0.0] + [r[ids.filtrate_flow(ids.zone(k))] for k in range(1, N + 1)]
    total = [p + b + f for p, b, f in zip(plasma, rbc, filtrate, strict=True)]
    assert total == pytest.approx([blood[0]] * (N + 1), rel=1e-12)


def test_a_solute_in_the_erythrocytes_is_not_enriched_by_the_filtration(
    tmp_path: Path,
) -> None:
    """The erythrocyte conduit of a pressure flow keeps the concentration of its inflow.

    The regression test of the conduit: a substance which states no
    erythrocyte distribution rides the erythrocytes and exchanges with
    nothing, so its concentration is that of the inflow in every zone. With
    the share `Hct` of the local flow the conduit lost water where the plasma
    filtered and the solute was enriched by `Q(z1)/Q(z<k>)`, 27 % at the
    middle of a ten zone rat liver.
    """
    r = _runner(tmp_path, None)
    for layer in (Layer.PLASMA, Layer.RBC):
        r[ids.species("tr", layer, ids.INFLOW)] = 1.0
    r.simulate(0.0, 400.0, 2)
    assert _profile(r, Layer.RBC) == pytest.approx([1.0] * N, rel=1e-6)
    assert _profile(r, Layer.PLASMA) == pytest.approx([1.0] * N, rel=1e-6)


# --- the exchange -----------------------------------------------------------


def test_a_substance_without_erythrocytes_has_no_exchange() -> None:
    """The behaviour before the exchange existed: no reaction and no parameter."""
    assert not [r for r in _reactions(None) if r.startswith("RBC_")]
    assert "Kbc_tr" not in _parameters(None)
    assert "kex_bc_tr" not in _parameters(None)


def test_a_substance_with_erythrocytes_exchanges_in_every_zone() -> None:
    """One reaction per zone, with the rate law of the PBPK library."""
    reactions = _reactions(EQUILIBRIUM)
    assert {f"RBC_tr_{ids.zone(k)}" for k in range(1, N + 1)} <= reactions
    (rbc,) = [
        r for r in _builder(EQUILIBRIUM).build().reactions if str(r.sid) == "RBC_tr_z3"
    ]
    assert rbc.formula is not None
    assert rbc.formula.math == (
        "kex_bc_tr * V_rbc_z3 * (tr_plasma_z3 - tr_rbc_z3 / Kbc_tr)"
    )


@pytest.mark.parametrize(
    ("geometry", "haematocrit"),
    [(Geometry.WELL_STIRRED, HCT), (Geometry.SINUSOID, 0.0)],
    ids=["well_stirred", "hct0"],
)
def test_a_lumen_which_is_one_compartment_has_no_exchange(
    geometry: Geometry, haematocrit: float
) -> None:
    """Nothing to exchange without an erythrocyte layer of its own.

    The well-stirred geometry holds the erythrocytes and the plasma in one
    compartment, so they are in equilibrium by construction, and a
    haematocrit of 0 leaves out the erythrocyte layer altogether.
    """
    reactions = _reactions(
        EQUILIBRIUM,
        geometry=geometry,
        haematocrit=haematocrit,
        flow=ImposedFlow(0.015),
    )
    assert not [r for r in reactions if r.startswith("RBC_")]


def test_the_erythrocytes_of_a_substance_in_equilibrium_follow_the_plasma(
    tmp_path: Path,
) -> None:
    """The profile is flat and the erythrocytes hold `Kbc` times the plasma.

    A freely sieved substance in equilibrium with its erythrocytes keeps the
    concentration of the inflow in both conduits of every zone: the plasma
    which filters takes it with it and brings it back, and the erythrocytes
    hold their coefficient everywhere.
    """
    r = _steady(tmp_path, EQUILIBRIUM)
    kbc = r["Kbc_tr"]
    assert kbc == pytest.approx(5.0, rel=1e-12)
    assert _profile(r, Layer.PLASMA) == pytest.approx([1.0] * N, rel=1e-3)
    assert _profile(r, Layer.RBC) == pytest.approx([kbc] * N, rel=1e-3)
    assert _profile(r, Layer.DISSE) == pytest.approx([1.0] * N, rel=1e-3)


def test_the_two_conduits_carry_the_blood_flux_of_the_pbpk_library(
    tmp_path: Path,
) -> None:
    """The chain carries `Q * BP * c_plasma`, the blood flux of the library.

    `BP = 1 - Hct + Hct*Kbc` (`pkpdsim.pbpk.kp.blood_to_plasma`): the sum of
    the two conduits of the engine is the one flux the PBPK library writes for
    erythrocytes in equilibrium, which is what makes the hepatic extraction of
    the two agree.
    """
    r = _steady(tmp_path, EQUILIBRIUM)
    bp = blood_to_plasma(_tracer(EQUILIBRIUM), HCT)
    assert bp == pytest.approx(1.0 - HCT + HCT * 5.0, rel=1e-12)
    # across the periportal face of the chain, where the blood enters at the
    # haematocrit of the settings; from there on the plasma carries the leak
    # of the leaky tube alone, so the split of a face inside the chain is not
    # that of the haematocrit
    face = ids.INFLOW
    flux = (
        r[ids.zone_flow(Layer.PLASMA, face)]
        * r[ids.species("tr", Layer.PLASMA, ids.INFLOW)]
        + r[ids.zone_flow(Layer.RBC, face)]
        * r[ids.species("tr", Layer.RBC, ids.INFLOW)]
    )
    expected = (
        r[ids.zone_flow(None, face)]
        * bp
        * r[ids.species("tr", Layer.PLASMA, ids.INFLOW)]
    )
    assert flux == pytest.approx(expected, rel=1e-12)


def test_the_inflow_of_the_erythrocytes_in_equilibrium_follows_the_plasma(
    tmp_path: Path,
) -> None:
    """`<s>_rbc_in` is no independent input of a substance in equilibrium.

    The mode states that the erythrocytes hold no state of their own, so the
    blood entering the chain carries them in equilibrium with its plasma: an
    assignment rule sets the inflow, at any value of the plasma.
    """
    r = _runner(tmp_path, EQUILIBRIUM)
    for c in (1.0, 0.25):
        r[ids.species("tr", Layer.PLASMA, ids.INFLOW)] = c
        assert r[ids.species("tr", Layer.RBC, ids.INFLOW)] == pytest.approx(
            r["Kbc_tr"] * c, rel=1e-12
        )


def test_the_rate_of_the_equilibrium_is_fast_against_the_transit_of_a_zone() -> None:
    """The derived rate relaxes the difference `EQUILIBRATION` times per transit.

    The relaxation of `c_plasma - c_rbc/Kbc` has the rate
    `kex_bc * (V_rbc/V_plasma + 1/Kbc)`; it is `EQUILIBRATION` reciprocal
    residence times of the fastest zone, which is what the engine means by an
    exchange in equilibrium.
    """
    builder = _builder(EQUILIBRIUM)
    values = _parameters(EQUILIBRIUM)
    transit = builder._transit()
    relaxation = values["kex_bc_tr"] * (HCT / (1.0 - HCT) + 1.0 / values["Kbc_tr"])
    assert transit > 0.0
    assert relaxation * transit == pytest.approx(EQUILIBRATION, rel=1e-12)


def test_the_coefficient_is_inverted_from_the_blood_to_plasma_ratio() -> None:
    """A substance which gives `bp` has `Kbc = (BP - 1 + Hct)/Hct`.

    The relation `BP = 1 - Hct + Hct*Kbc` of the PBPK library, read the other
    way: the engine needs the concentration of the erythrocytes, a study
    measures the ratio of the whole blood. The round trip through
    `blood_to_plasma` of the library returns the ratio.
    """
    ery = Erythrocytes(mode="equilibrium", bp=1.3)
    kbc = blood_cell_coefficient(ery, HCT)
    assert kbc == pytest.approx((1.3 - 1.0 + HCT) / HCT, rel=1e-12)
    assert blood_to_plasma(_tracer(ery), HCT) == pytest.approx(1.3, rel=1e-12)
    assert _parameters(ery)["Kbc_tr"] == pytest.approx(kbc, rel=1e-12)


def test_a_kinetic_exchange_fills_the_erythrocytes_along_the_zones(
    tmp_path: Path,
) -> None:
    """A rate of its own, and the erythrocytes lag behind their equilibrium.

    The point of the kinetic mode: with a rate of the order of the transit of
    the chain the erythrocytes are still filling when the blood leaves, so
    they rise from the periportal to the pericentral end and stay below
    `Kbc * c_plasma`, while the plasma falls because the erythrocytes take
    from it.
    """
    assert _parameters(KINETIC)["kex_bc_tr"] == pytest.approx(0.2, rel=1e-12)
    r = _steady(tmp_path, KINETIC)
    rbc = _profile(r, Layer.RBC)
    plasma = _profile(r, Layer.PLASMA)
    assert rbc == sorted(rbc)
    assert plasma == sorted(plasma, reverse=True)
    assert 0.0 < rbc[0] < rbc[-1] < r["Kbc_tr"] * plasma[-1]
    # the solute is conserved: what the two conduits carry out is what the
    # inflow brought in, whatever the erythrocytes hold on the way
    last = ids.zone(N)
    out = (
        r[ids.zone_flow(Layer.PLASMA, last)] * plasma[-1]
        + r[ids.zone_flow(Layer.RBC, last)] * rbc[-1]
    )
    into = r[ids.zone_flow(Layer.PLASMA, ids.INFLOW)] * 1.0
    assert out == pytest.approx(into, rel=1e-6)


def test_the_balance_of_the_observers_stays_at_zero_with_the_exchange(
    tmp_path: Path,
) -> None:
    """The exchange moves solute between two layers of the organ and creates none."""
    r = _runner(tmp_path, EQUILIBRIUM, observed=True)
    for layer in (Layer.PLASMA, Layer.DISSE, Layer.ENDOTHELIUM, Layer.CELL):
        for k in range(1, N + 1):
            r[f"[{ids.species('tr', layer, ids.zone(k))}]"] = 0.5
    observers.sync_dose(r, [_tracer(EQUILIBRIUM)])
    dose = r[observers.ledger("dose", "tr")]
    r[ids.species("tr", Layer.PLASMA, ids.INFLOW)] = 1.0
    r.simulate(0.0, 20.0, 200)
    scale = dose + r[observers.ledger("in", "tr")]
    assert dose > 0.0
    assert r[observers.ledger("in", "tr")] > 0.0
    assert abs(r["balance_tr"]) <= 1e-10 * scale


# --- what the engine refuses ------------------------------------------------


def test_a_blood_to_plasma_ratio_the_erythrocytes_cannot_hold_is_refused() -> None:
    """A ratio at or below `1 - Hct` needs an erythrocyte concentration of 0 or less."""
    with pytest.raises(ValueError, match="at or below 1 - Hct"):
        blood_cell_coefficient(Erythrocytes(mode="equilibrium", bp=0.5), HCT)
    with pytest.raises(ValueError, match="cannot be inverted at a haematocrit of 0"):
        blood_cell_coefficient(Erythrocytes(mode="equilibrium", bp=1.3), 0.0)


def test_a_chain_which_filters_more_plasma_than_it_carries_is_refused() -> None:
    """A tube without an oncotic pressure cannot carry erythrocytes.

    The leaky tube knows the flow of the lumen and not what is in it: without
    an oncotic counter-pressure it filters three quarters of its inflow away,
    more than the plasma of a haematocrit of 0.45, and the erythrocytes would
    have to pass the fenestrae. The engine refuses the pair instead of letting
    the plasma conduit run negative; an acellular perfusate (a haematocrit of
    0) is the case that pressure belongs to.
    """
    with pytest.raises(ValueError, match="filter more plasma than enters it"):
        _builder(None, flow=PressureFlow(p_onc=Value(0.0, "Pa"))).build()
    _builder(None, flow=PressureFlow(p_onc=Value(0.0, "Pa")), haematocrit=0.0).build()


def test_erythrocytes_in_equilibrium_need_a_flow() -> None:
    """Without a flow the engine has no time scale the exchange could be fast against."""
    with pytest.raises(ValueError, match="Erythrocytes in equilibrium need a flow"):
        _builder(EQUILIBRIUM, flow=ImposedFlow(0.0)).build()


def test_a_kinetic_rate_in_another_unit_is_refused() -> None:
    """The rate of the engine is 1/min, and a value in 1/hr is not read as one."""
    ery = Erythrocytes(mode="kinetic", kbc=5.0, k_exchange=Value(0.2, U.per_hr))
    with pytest.raises(ValueError, match="k_exchange must be in 1/min"):
        _builder(ery).build()


def test_a_kinetic_exchange_without_a_rate_is_refused() -> None:
    """The kinetic mode is the rate: there is nothing to derive it from."""
    ery = Erythrocytes(mode="kinetic", kbc=5.0)
    with pytest.raises(ValueError, match="need a positive exchange rate"):
        _builder(ery).build()
