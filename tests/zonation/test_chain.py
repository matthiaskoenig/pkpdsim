"""Convection, axial diffusion and the boundary of the chain of zones.

The tests use a tracer, a substance which does not react, so no cell model is
needed. The chain is simulated with roadrunner as the experiments are.
"""

import dataclasses
import warnings
from pathlib import Path

import numpy as np
import pytest
import roadrunner
from pbpk.conftest import individual
from sbmlutils.factory import Model, Parameter, create_model
from sbmlutils.metadata import SBO

from pkpdsim.pbpk import MODEL_UNITS
from pkpdsim.pbpk import ids as pbpk_ids
from pkpdsim.pbpk.physiology import Sex, physiology_objects, reference_flow
from pkpdsim.pbpk.placement import value_parameter
from pkpdsim.pbpk.substance import Substance
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Formula, Value
from pkpdsim.species import Species
from pkpdsim.zonation import ZonationBuilder, effective_peclet, ids
from pkpdsim.zonation.flow import ImposedFlow
from pkpdsim.zonation.geometry import Layer, axial_area, fractions, geometry
from pkpdsim.zonation.ids import flow
from pkpdsim.zonation.profiles import Constant
from pkpdsim.zonation.settings import Geometry, ZonationSettings
from pkpdsim.zonation.transport import M2_PER_S, Solute

from .conftest import liver_flow, mass_g

TRACER = Substance("tr", "tracer", Value(100.0, U.g_per_mole))
"""A substance which does not react; the builder of these tests gives it the radius
of an intravascular marker (`_Diffusing._radius`), which cannot cross the fenestrae."""

_MASS = mass_g(1000.0)


class _Diffusing(ZonationBuilder):
    """A builder whose tracer diffuses with a given coefficient (0: it does not)."""

    diffusion_value = 0.0

    def _radius(self, substance: Solute) -> Value:
        """Get a radius above that of a fenestra: the tracer stays in the lumen."""
        return Value(1e-6, U.meter)

    def _diffusion(self, substance: Solute) -> Parameter | None:
        """Get the diffusion of the tracer, `None` for none."""
        if not self.diffusion_value:
            return None
        return value_parameter(
            f"D_{substance.sid}",
            Value(self.diffusion_value, M2_PER_S),
            M2_PER_S,
            "diffusion tracer",
            SBO.DIFFUSION_COEFFICIENT,
        )


class _Permeable(_Diffusing):
    """A builder whose tracer is small enough to cross the fenestrae."""

    def _radius(self, substance: Solute) -> Value:
        """Get the radius of a small solute, 0.365 nm against 53.5 nm of a fenestra."""
        return Value(0.365e-9, U.meter)


def _settings(
    n: int = 10, flow: float | None = None, **kwargs: object
) -> ZonationSettings:
    """The settings of a chain of `n` zones, the flow in l/min."""
    kwargs.setdefault("reference_flow", liver_flow())
    return ZonationSettings(
        n_zones=n,
        flow=ImposedFlow(flow),
        liver_mass=_MASS,
        **kwargs,  # ty: ignore[invalid-argument-type]
    )


def _builder(
    n: int, flow: float, diffusion: float, **kwargs: object
) -> ZonationBuilder:
    """A builder of a tracer chain."""
    builder = _Diffusing(_settings(n, flow, **kwargs), [TRACER], None)
    builder.diffusion_value = diffusion
    return builder


def _tracer_chain(
    n: int, flow: float, diffusion: float, tmp_path: Path, **kwargs: object
) -> roadrunner.RoadRunner:
    """A bolus of the tracer (an intravascular marker) in the plasma of the first zone."""
    model = _builder(n, flow, diffusion, **kwargs).build()
    path = tmp_path / "chain.xml"
    create_model(model, filepath=path, sbml_level=3, sbml_version=2, validate=False)
    r = individual(path)
    r.integrator.absolute_tolerance = 1e-14
    r.integrator.relative_tolerance = 1e-10
    r[f"[{ids.species('tr', Layer.PLASMA, ids.zone(1))}]"] = 1.0
    return r


def _plasma_ids(r: roadrunner.RoadRunner) -> list[str]:
    """The plasma tracer species in the order of the zones."""
    n = sum(
        1 for sid in r.model.getFloatingSpeciesIds() if sid.startswith("tr_plasma_")
    )
    return [ids.species("tr", Layer.PLASMA, ids.zone(k)) for k in range(1, n + 1)]


def _volume(r: roadrunner.RoadRunner) -> float:
    """The plasma volume (l) of all zones."""
    return sum(
        r[ids.volume(Layer.PLASMA, f"z{k}")] for k in range(1, len(_plasma_ids(r)) + 1)
    )


def _total(r: roadrunner.RoadRunner) -> float:
    """The amount (mmole) of the tracer in all layers and zones."""
    total = 0.0
    for sid in r.model.getFloatingSpeciesIds():
        layer_zone = sid.removeprefix("tr_")
        total += r[f"[{sid}]"] * r[f"V_{layer_zone}"]
    return total


def _outflow(r: roadrunner.RoadRunner, horizon: float) -> tuple[np.ndarray, np.ndarray]:
    """The time and the concentration at the pericentral end of the plasma."""
    last = _plasma_ids(r)[-1]
    out = r.simulate(0, horizon, 40001, selections=["time", f"[{last}]"])
    return out[:, 0], out[:, 1]


def _moments(r: roadrunner.RoadRunner) -> tuple[float, float]:
    """The mean and the normalized variance of the outflow of the bolus."""
    tau = _volume(r) / r["Q_plasma_li"]
    t, c = _outflow(r, 12.0 * tau)
    area = np.trapezoid(c, t)
    mean = np.trapezoid(t * c, t) / area
    variance = np.trapezoid((t - mean) ** 2 * c, t) / area
    return float(mean), float(variance / mean**2)


def _mrt(r: roadrunner.RoadRunner) -> float:
    """The mean residence time (min) of the bolus."""
    return _moments(r)[0]


def _outflow_variance(r: roadrunner.RoadRunner) -> float:
    """The normalized variance (variance over the squared mean) of the outflow."""
    return _moments(r)[1]


def _relative_outflow_auc(n: int, tmp_path: Path) -> float:
    """The outflow AUC of a bolus through `n` zones over the AUC of plug flow.

    The dose is the amount of the bolus in the plasma of the first zone, which
    the zone count divides, so the bolus is normalized to the same dose at
    every zone count: the concentration `dose / V_plasma_z1`, a bolus which
    enters the chain at the periportal end and becomes a point source there as
    the zones grow. Its outflow AUC over the whole elution is `dose /
    Q_plasma` for a chain which carries the dose out and holds nothing back,
    the AUC of plug flow and of any other discretization of it, so the
    quotient is 1 and the convergence of the AUC in the zone count is the
    convergence of that quotient.

    Args:
        n: the number of zones.
        tmp_path: the directory of the model.

    Returns:
        The AUC of the outflow concentration over `dose / Q_plasma`.
    """
    dose = 1e-3
    r = _tracer_chain(n, 1.0, 0.0, tmp_path, haematocrit=0.0)
    r[f"[{ids.species('tr', Layer.PLASMA, ids.zone(1))}]"] = (
        dose / r[ids.volume(Layer.PLASMA, ids.zone(1))]
    )
    tau = _volume(r) / r["Q_plasma_li"]
    t, c = _outflow(r, 20.0 * tau)
    return float(np.trapezoid(c, t)) * r["Q_plasma_li"] / dose


def test_a_closed_chain_conserves_the_substance(tmp_path: Path) -> None:
    """Without a flow only the diffusion moves the tracer, and none leaves."""
    r = _tracer_chain(10, 0.0, 1e-9, tmp_path)
    before = _total(r)
    r.simulate(0, 60, 61)
    assert _total(r) == pytest.approx(before, rel=1e-6)
    assert r["[tr_plasma_z10]"] > 0.0, "the diffusion moved nothing"


def test_mean_residence_time_is_volume_over_flow(tmp_path: Path) -> None:
    """The mean of a bolus through the chain is the plasma volume over the flow.

    The tracer is an intravascular marker (1 um, it cannot cross the fenestrae),
    which is why it stays in the plasma: this is no claim about every solute.
    A small solute distributes into the space of Disse and stays longer
    (`test_barrier.py`).

    The first moment of a chain of stirred tanks is `V / Q` whatever the zone
    count, so the only error here is the integration of the curve: measured,
    the quotient is 1.0000000000070 (7e-12), the same number in three runs.
    The tolerance is 1e-8, a margin of a thousand over that and a tenth of
    what the sibling moment tests use, where a zone ten percent too large
    passed at the previous 1/10.
    """
    r = _tracer_chain(10, 1.0, 0.0, tmp_path, haematocrit=0.0)
    assert _mrt(r) == pytest.approx(_volume(r) / 1.0, rel=1e-8)


def test_the_outflow_auc_converges_in_the_zone_count(tmp_path: Path) -> None:
    """The outflow AUC of a tracer at ten zones is within 2 % of forty zones.

    The acceptance criterion of the design for the discretization of the
    chain, on which every model of the engine stands at ten zones. Measured,
    the two agree to 4.3e-11 of each other and each is the AUC of plug flow
    (`dose / Q_plasma`, an oracle outside the chain) to 1e-10, so ten zones
    meet the criterion with nine orders of magnitude to spare; the assertions
    are at 1e-8, a hundredfold margin over the integration of the curve.

    The criterion is met that well because the AUC of a conservative tracer is
    fixed by the mass balance and not by the discretization: whatever the zone
    count smears, the whole dose leaves through the pericentral end at the
    plasma flow. What the test pins is therefore that the chain carries the
    dose out and holds none of it back at any zone count, and the chain's
    discretization error is in the shape of the outflow and not in its area:
    it is the variance of `test_numerical_dispersion_matches_two_n`, which is
    `1/n` and falls with the zone count (0.1 at ten zones, 0.025 at forty).
    """
    coarse = _relative_outflow_auc(10, tmp_path)
    fine = _relative_outflow_auc(40, tmp_path)
    assert coarse == pytest.approx(fine, rel=0.02)
    assert coarse == pytest.approx(fine, rel=1e-8)
    assert coarse == pytest.approx(1.0, rel=1e-8)
    assert fine == pytest.approx(1.0, rel=1e-8)


@pytest.mark.parametrize("n", [4, 10])
def test_numerical_dispersion_matches_two_n(n: int, tmp_path: Path) -> None:
    """The upwind chain of n zones has the Peclet number 2n.

    The tracer is an intravascular marker which cannot cross the fenestrae, so
    only the convection spreads the bolus. The normalized variance of `n`
    stirred tanks in series is `2 / Pe` = `1 / n`, which the chain realizes
    to the integration of the curve: measured 0.10000000006 at ten zones and
    0.2500000000013 at four, the same numbers in three runs. The tolerance is
    1e-6, that of the sibling moment tests and four orders tighter than the
    1/4 which let a Peclet number a quarter off pass; two zone counts, so the
    `1 / n` and not one point of it is pinned.
    """
    variance = _outflow_variance(_tracer_chain(n, 1.0, 0.0, tmp_path, haematocrit=0.0))
    assert variance == pytest.approx(2.0 / (2 * n), rel=1e-6)  # Pe = 2n


@pytest.mark.parametrize("flow,diffusion", [(0.0, 0.0), (1.0, 0.0), (0.0, 1e-9)])
def test_zero_flow_or_zero_diffusion_stays_finite(
    flow: float, diffusion: float, tmp_path: Path
) -> None:
    """A rate law never divides by a zero flow or diffusion."""
    r = _tracer_chain(10, flow, diffusion, tmp_path)
    assert np.all(np.isfinite(r.simulate(0, 60, 61)))


def test_a_diffusion_the_numerical_dispersion_dominates_warns() -> None:
    """The zone count and not the physics sets the gradient then."""
    with pytest.warns(UserWarning, match="numerical dispersion"):
        _builder(10, 1.0, 1e-14).build()


def test_a_diffusion_above_the_numerical_dispersion_does_not_warn() -> None:
    """Nothing to say when the physics sets the gradient."""
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        _builder(10, 1.0, 1e-3).build()


def test_the_warning_fires_at_the_factor_of_ten_and_not_before() -> None:
    """The threshold is pinned on both sides of `DISPERSION_FACTOR`.

    The chain of these tests disperses like 7.429245e-10 m2/s (the arithmetic
    of `test_a_physical_diffusion_lowers_the_peclet_number`: 1 l/min through
    the 0.1009524 l of lumen of 1 kg of liver, ten zones of 300 um), so the
    limit of the warning is a tenth of it, 7.429245e-11 m2/s. A diffusion 1 %
    below the limit warns and one 1 % above it does not, and the message names
    the factor and the dispersion: the two were tested at 1e-14 and 1e-3
    before, five orders either side of the limit, so the factor could move
    without failing.
    """
    limit = 7.429245e-11
    with pytest.warns(
        UserWarning,
        match=r"more than 10 times below the numerical dispersion 7\.43e-10",
    ):
        _builder(10, 1.0, 0.99 * limit).build()
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        _builder(10, 1.0, 1.01 * limit).build()


@pytest.mark.parametrize(
    ("n", "diffusion", "geometry"),
    [
        (10, 0.9e-9, Geometry.SINUSOID),
        (1, 1e-14, Geometry.SINUSOID),
        (10, 1e-14, Geometry.WELL_STIRRED),
    ],
    ids=["same_order", "one_zone", "well_stirred"],
)
def test_the_builds_a_caller_cannot_avoid_do_not_warn(
    n: int, diffusion: float, geometry: Geometry
) -> None:
    """The warning has to be one a correct caller can avoid.

    The numerical dispersion of the engine's own defaults is of the order of
    the aqueous diffusion of a small solute (1.63e-9 m**2/s against 0.9e-9 of
    glucose at ten zones of a rat liver), so a comparison with the bare
    dispersion fired for every physiological build and five test modules muted
    it; it is `DISPERSION_FACTOR` now. A chain without an axial diffusion (one
    zone, the well-stirred geometry) is not compared at all: it has no axial
    gradient to get wrong.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        _builder(n, 1.0, diffusion, geometry=geometry).build()


def _equivalence_chain(
    geo: Geometry, tmp_path: Path, porosity: float, p_endo: float
) -> roadrunner.RoadRunner:
    """A bolus in the carrier of the first zone of a chain of the given geometry.

    One kg of rat liver, ten zones, 1 l/min imposed, no erythrocytes and no
    axial diffusion, and the same geometry constants in both geometries, so
    the geometry is the only difference. `porosity` is the fenestral porosity
    of every zone and `p_endo` the transcellular permeability of the
    endothelium (m/s), the two paths of the barrier; the well-stirred
    geometry has neither, having no barrier. The dose is an amount, so the
    two carriers receive the same dose although their volumes differ.
    """
    constants = dataclasses.replace(
        geometry(Species.RAT), porosity_pp=porosity, porosity_pc=porosity
    )
    builder = _Permeable(
        ZonationSettings(
            n_zones=10,
            geometry=geo,
            flow=ImposedFlow(1.0),
            liver_mass=_MASS,
            haematocrit=0.0,
            reference_flow=liver_flow(),
            geometry_constants=constants,
        ),
        [TRACER],
        None,
    )
    path = tmp_path / f"{geo}_{porosity}_{p_endo}.xml"
    create_model(
        builder.build(), filepath=path, sbml_level=3, sbml_version=2, validate=False
    )
    r = (
        individual(path)
        if geo is Geometry.WELL_STIRRED
        else individual(path, P_endo_in_tr=p_endo, P_endo_out_tr=p_endo)
    )
    r.integrator.absolute_tolerance = 1e-14
    r.integrator.relative_tolerance = 1e-10
    r[ids.species("tr", _carrier(geo), ids.zone(1))] = 1.0
    return r


def _carrier(geo: Geometry) -> Layer | str:
    """The layer the blood carries the tracer in, by the geometry."""
    return ids.VASCULAR if geo is Geometry.WELL_STIRRED else Layer.PLASMA


def _carrier_volume(r: roadrunner.RoadRunner, geo: Geometry) -> float:
    """The volume (l) of every compartment of the chain which holds the tracer."""
    layers = (
        [ids.VASCULAR]
        if geo is Geometry.WELL_STIRRED
        else [Layer.PLASMA, Layer.ENDOTHELIUM, Layer.DISSE]
    )
    return sum(
        r[ids.volume(layer, ids.zone(k))] for layer in layers for k in range(1, 11)
    )


def _pericentral_curve(r: roadrunner.RoadRunner, geo: Geometry) -> np.ndarray:
    """The concentration at the pericentral end over twelve minutes."""
    last = ids.species("tr", _carrier(geo), ids.zone(10))
    out = r.simulate(0, 12.0, 2001, selections=["time", f"[{last}]"])
    return np.asarray(out)[:, 1]


def test_the_sinusoid_with_an_open_barrier_is_the_well_stirred_chain(
    tmp_path: Path,
) -> None:
    """The limit the design states, on the outflow curve of the engine itself.

    "The sinusoid with the diffusion off and the barrier open equals well
    stirred": the well-stirred geometry holds the four layers of the lumen in
    one compartment, so a sinusoid whose barrier does not limit the transfer
    has to give the same outflow. The two builds differ in the geometry alone
    (ten zones, 1 kg of rat liver, 1 l/min imposed, no erythrocytes, no axial
    diffusion) and hold the same volume, 0.1742857 l in both: the 0.106 + 0.028
    + 0.049 of the parenchyma as three compartments per zone and as one.

    The barrier is opened with a fenestral porosity of 1 (a fully fenestrated
    endothelium, the table has 6 to 8 %) and a transcellular permeability of
    1 m/s, which makes the exchange of a zone about 7e4 times its flow.
    Measured, the two curves agree to 9.0e-7 of the peak, and to 9.0e-8 at ten
    times the two permeabilities, so the residual is the finite barrier and
    not a difference of the chains; the assertion is at 1e-5, a margin of ten
    over the measurement.

    The test is not vacuous. Opening only the fenestrae leaves the
    endothelium, 15 % of the volume of the lumen, out of the exchange and the
    curves differ by 0.41 of the peak; a barrier which is nearly closed keeps
    the tracer in the plasma, a third of the volume, and differs by 1.25.
    """
    stirred = _equivalence_chain(Geometry.WELL_STIRRED, tmp_path, 1.0, 1.0)
    opened = _equivalence_chain(Geometry.SINUSOID, tmp_path, 1.0, 1.0)
    volume = _carrier_volume(stirred, Geometry.WELL_STIRRED)
    assert volume == pytest.approx(0.1742857, rel=1e-6)
    assert _carrier_volume(opened, Geometry.SINUSOID) == pytest.approx(
        volume, rel=1e-12
    )

    reference = _pericentral_curve(stirred, Geometry.WELL_STIRRED)
    peak = float(np.max(reference))
    assert peak > 0.0

    def deviation(r: roadrunner.RoadRunner) -> float:
        """The largest difference from the well-stirred curve, over its peak."""
        curve = _pericentral_curve(r, Geometry.SINUSOID)
        return float(np.max(np.abs(curve - reference))) / peak

    assert deviation(opened) < 1e-5
    assert deviation(_equivalence_chain(Geometry.SINUSOID, tmp_path, 1.0, 0.0)) > 0.4
    assert deviation(_equivalence_chain(Geometry.SINUSOID, tmp_path, 1e-9, 0.0)) > 1.0


def test_the_reported_peclet_number_is_the_physical_one() -> None:
    """Without physical diffusion the chain realizes 2n."""
    assert effective_peclet(_settings(n=10), d_eff=0.0) == pytest.approx(20.0)


def test_a_physical_diffusion_lowers_the_peclet_number() -> None:
    """The numerical and the physical dispersion add, by the arithmetic below.

    1 kg of liver at the rat density 1.05 kg/l is 0.9523810 l of parenchyma,
    the lumen 10.6 % of it: 0.1009524 l. A flow of 1 l/min needs 0.1009524 min
    = 6.057143 s to replace it, so the blood crosses the 300 um of the sinusoid
    at 300e-6 / 6.057143 = 4.952830e-5 m/s (49.53 um/s) and `v * L` is
    1.485849e-8 m2/s.

    Ten zones disperse like `v * L / (2 * 10)` = 7.429245e-10 m2/s, so the
    Peclet number `v * L / (d_eff + v * L / 20)` is 1.485849e-8 / 7.429245e-10
    = 20 without a physical diffusion and 1.485849e-8 / 1.742925e-9 = 8.525034
    with 1e-9 m2/s; four zones disperse like 1.857311e-9 m2/s, and with the
    0.9e-9 m2/s of glucose the number is 1.485849e-8 / 2.757311e-9 = 5.388761.

    The three numbers are the division of the two quantities written out and
    not the expression of `peclet`, so a `d_eff` which is scaled before it is
    added (a tenth of it leaves 15.46 at ten zones) or a dispersion without
    its factor of two fails here.
    """
    free = effective_peclet(_settings(n=10, flow=1.0), d_eff=0.0)
    assert free == pytest.approx(20.0, rel=1e-9)
    assert effective_peclet(_settings(n=10, flow=1.0), d_eff=1e-9) == pytest.approx(
        8.525034, rel=1e-6
    )
    assert effective_peclet(_settings(n=4, flow=1.0), d_eff=0.9e-9) == pytest.approx(
        5.388761, rel=1e-6
    )
    assert effective_peclet(_settings(n=10, flow=0.0), d_eff=1e-9) == 0.0


@pytest.mark.parametrize("haematocrit", [None, 0.0, 0.45])
def test_the_chain_runs_with_and_without_erythrocytes(
    haematocrit: float | None, tmp_path: Path
) -> None:
    """The erythrocyte chain exists with a haematocrit above 0 only."""
    builder = _builder(4, 1.0, 1e-9, haematocrit=haematocrit)
    reactions = {r.sid for r in builder.convection()}
    rbc = flow("tr", "rbc_z1", "rbc_z2")
    assert (rbc in reactions) is (haematocrit != 0.0)
    r = _tracer_chain(4, 1.0, 1e-9, tmp_path, haematocrit=haematocrit)
    assert np.all(np.isfinite(r.simulate(0, 5, 51)))


def test_the_plasma_flow_is_the_blood_flow_without_the_erythrocytes(
    tmp_path: Path,
) -> None:
    """The plasma carries the share of the flow which is not erythrocytes."""
    r = _tracer_chain(10, 1.0, 0.0, tmp_path, haematocrit=0.45)
    assert r["Q_plasma_li"] == pytest.approx(0.55)
    assert r["Q_rbc_li"] == pytest.approx(0.45)
    # measured 1.000000000052 of V / Q_plasma, as in
    # `test_mean_residence_time_is_volume_over_flow`
    assert _mrt(r) == pytest.approx(_volume(r) / 0.55, rel=1e-8)


def test_without_erythrocytes_no_reaction_names_them() -> None:
    """Every reaction of the chain refers to a compartment which exists."""
    model = _builder(5, 1.0, 1e-9, haematocrit=0.0).build()
    assert not any("rbc" in str(r.sid) for r in model.reactions or [])
    compartments = {str(c.sid) for c in model.compartments or []}
    assert all(str(r.compartment) in compartments for r in model.reactions or [])


def test_a_single_zone_has_only_the_inflow_and_the_outflow(tmp_path: Path) -> None:
    """No convection and no diffusion between zones, and it simulates."""
    builder = _builder(1, 1.0, 1e-9, haematocrit=0.0)
    assert [str(r.sid) for r in builder.convection()] == [
        "FLOW_tr_plasma_in_plasma_z1",
        "FLOW_tr_plasma_z1_plasma_out",
    ]
    assert builder.axial_diffusion() == []
    r = _tracer_chain(1, 1.0, 1e-9, tmp_path)
    assert np.all(np.isfinite(r.simulate(0, 5, 51)))
    assert r["[tr_plasma_z1]"] < 1.0, "the outflow carries the tracer away"


def test_the_reaction_ids() -> None:
    """The convection and the diffusion follow the id scheme of the engine."""
    builder = _builder(3, 1.0, 1e-9, haematocrit=0.0)
    convection = {str(r.sid) for r in builder.convection()}
    assert convection == {
        "FLOW_tr_plasma_in_plasma_z1",
        "FLOW_tr_plasma_z1_plasma_z2",
        "FLOW_tr_plasma_z2_plasma_z3",
        "FLOW_tr_plasma_z3_plasma_out",
    }
    diffusion = {str(r.sid) for r in builder.axial_diffusion()}
    assert diffusion == {
        f"D_tr_{layer}_z{k}_{layer}_z{k + 1}"
        for layer in ("plasma", "disse")
        for k in (1, 2)
    }


def test_a_zero_flow_or_diffusion_emits_no_reaction() -> None:
    """A rate law which would be zero is not emitted."""
    assert _builder(5, 0.0, 1e-9).convection() == []
    assert _builder(5, 1.0, 0.0).axial_diffusion() == []


def test_the_diffusion_conductance_is_the_axial_area_over_the_zone_length(
    tmp_path: Path,
) -> None:
    """`D * A / dx` of `axial_area`, in l/min, is what the reaction is built on."""
    r = _tracer_chain(10, 0.0, 1e-9, tmp_path)
    geo = geometry(_settings().species)
    fr = dict(fractions(_settings().species))
    area = axial_area(Layer.PLASMA, _MASS, geo, fr).value
    expected = 1e-9 * area / (geo.L_sin.value / 10) * 1e3 * 60.0
    # the rate of the first interface at a unit gradient
    r["[tr_plasma_z2]"] = 0.0
    rate = r["D_tr_plasma_z1_plasma_z2"]
    assert rate == pytest.approx(expected * (1.0 - 0.0), rel=1e-6)


def test_a_haematocrit_of_zero_has_no_erythrocyte_area() -> None:
    """The axial area of a layer without volume is an error, not a KeyError."""
    geo = geometry(_settings().species)
    fr = dict(fractions(_settings().species, 0.0))
    with pytest.raises(ValueError, match="no volume at this haematocrit"):
        axial_area(Layer.RBC, _MASS, geo, fr)


def test_the_well_stirred_geometry_flows_through_its_vascular_zones() -> None:
    """The vascular compartments of the zones are the chain."""
    s = _settings(3, 1.0, geometry=Geometry.WELL_STIRRED)
    reactions = {str(r.sid) for r in ZonationBuilder(s, [TRACER], None).convection()}
    assert "FLOW_tr_vas_z1_vas_z2" in reactions


def test_the_axial_diffusion_profile_scales_the_conductance() -> None:
    """A registered process scales its rate by the zonal factor."""
    s = _settings(3, 1.0, profiles={"axial_diffusion": Constant()})
    model = ZonationBuilder(s, [TRACER], None).build()
    assert ids.profile("axial_diffusion", "z2") in {
        str(p.sid) for p in model.parameters or []
    }


@pytest.mark.parametrize("species", [Species.RAT, Species.HUMAN])
def test_the_reference_flow_is_the_total_liver_flow_of_the_pbpk_library(
    species: Species, tmp_path: Path
) -> None:
    """`reference_flow` equals `QH = QHA + QPO` of the library's physiology model.

    This is the invariant which lets the zonated liver stand in for the PBPK
    liver: both carry the same blood flow. The library's physiology is built
    as SBML and `QH` is read from roadrunner, so the renormalisation of the
    flow fractions over the organs is the library's own. The value is the
    library's (`pkpdsim.pbpk.physiology.reference_flow`), which a caller of
    the engine passes as `ZonationSettings.reference_flow`: the engine reads
    no physiology.
    """
    path = tmp_path / "physiology.xml"
    model = Model(
        sid="physiology",
        units=U,
        model_units=MODEL_UNITS,
        objects=physiology_objects(Sex.MALE, species=species),
    )
    create_model(model, filepath=path, validate=False)
    qh = individual(path)[pbpk_ids.QH]
    assert reference_flow(species).value == pytest.approx(qh, rel=1e-9)


def test_the_dispersion_is_in_the_notes_of_the_model() -> None:
    """A silenced warning does not take the numbers with it."""
    model = _builder(10, 1.0, 1e-9).build()
    assert "Peclet number of 20" in str(model.notes)
    assert "m**2/s" in str(model.notes)


def test_a_formula_flow_names_the_flow_to_use_instead() -> None:
    """A computed flow is refused with the reason."""
    flow = ImposedFlow(Formula("1 + 1", U.l_per_min))
    with pytest.raises(ValueError, match="PressureFlow"):
        ZonationBuilder(ZonationSettings(flow=flow), [TRACER], None).build()


@pytest.mark.parametrize("haematocrit", [1.0, 1.5, -0.1])
def test_a_haematocrit_without_plasma_is_refused(haematocrit: float) -> None:
    """A haematocrit of 1 leaves no plasma, which the solver cannot integrate."""
    with pytest.raises(ValueError, match=r"outside \[0, 1\)"):
        ZonationSettings(haematocrit=haematocrit)


def test_a_negative_flow_is_refused() -> None:
    """A negative flow would reverse the blood chain."""
    with pytest.raises(ValueError, match="negative"):
        _builder(10, -1.0, 0.0).convection()
    with pytest.raises(ValueError, match="-1"):
        effective_peclet(_settings(n=10, flow=-1.0), d_eff=0.0)


def test_the_notes_the_warning_and_the_function_report_one_peclet_number() -> None:
    """The reported number is the one the builder computes, and it is the right one.

    The identity of the two paths says nothing on its own, so the shared value
    is pinned against the 8.525034 of
    `test_a_physical_diffusion_lowers_the_peclet_number`, where the arithmetic
    is written out.
    """
    builder = _builder(10, 1.0, 0.0)
    assert builder.peclet(1e-9) == effective_peclet(builder.settings, 1e-9)
    assert builder.peclet(1e-9) == pytest.approx(8.525034, rel=1e-6)
