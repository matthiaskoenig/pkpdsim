"""Convection, axial diffusion and the boundary of the chain of zones.

The tests use a tracer, a substance which does not react, so no cell model is
needed. The chain is simulated with roadrunner as the experiments are.
"""

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
    """
    r = _tracer_chain(10, 1.0, 0.0, tmp_path, haematocrit=0.0)
    assert _mrt(r) == pytest.approx(_volume(r) / 1.0, rel=1 / 10)


def test_numerical_dispersion_matches_two_n(tmp_path: Path) -> None:
    """The upwind chain of n zones has the Peclet number 2n.

    The tracer is an intravascular marker which cannot cross the fenestrae, so
    only the convection spreads the bolus.
    """
    variance = _outflow_variance(_tracer_chain(10, 1.0, 0.0, tmp_path, haematocrit=0.0))
    assert variance == pytest.approx(2.0 / 20.0, rel=0.25)  # Pe = 2n


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
    the aqueous diffusion of a small solute (1.66e-9 m**2/s against 0.9e-9 of
    glucose at ten zones of a rat liver), so a comparison with the bare
    dispersion fired for every physiological build and five test modules muted
    it; it is `DISPERSION_FACTOR` now. A chain without an axial diffusion (one
    zone, the well-stirred geometry) is not compared at all: it has no axial
    gradient to get wrong.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        _builder(n, 1.0, diffusion, geometry=geometry).build()


def test_the_reported_peclet_number_is_the_physical_one() -> None:
    """Without physical diffusion the chain realizes 2n."""
    assert effective_peclet(_settings(n=10), d_eff=0.0) == pytest.approx(20.0)


def test_a_physical_diffusion_lowers_the_peclet_number() -> None:
    """The numerical and the physical dispersion add."""
    free = effective_peclet(_settings(n=10, flow=1.0), d_eff=0.0)
    assert effective_peclet(_settings(n=10, flow=1.0), d_eff=1e-9) < free
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
    assert _mrt(r) == pytest.approx(_volume(r) / 0.55, rel=1 / 10)


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
    expected = 1e-9 * area / (geo.L_sin / 10) * 1e3 * 60.0
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
    """The reported number is the one the builder computes."""
    builder = _builder(10, 1.0, 0.0)
    assert builder.peclet(1e-9) == effective_peclet(builder.settings, 1e-9)
