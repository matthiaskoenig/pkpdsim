"""The canalicular chain: pericentral to periportal and out into the duct.

The tests use a tracer, a substance which does not react, so no cell model is
needed. The blood does not flow (`flow=0`), so only the bile moves the tracer,
and the tracer is released into the canaliculus of one zone.
"""

import math
from pathlib import Path

import numpy as np
import pytest
import roadrunner
from pbpk.conftest import individual
from sbmlutils.factory import create_model

from pkpdsim.pbpk.substance import Substance
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Value
from pkpdsim.species import Species
from pkpdsim.zonation import ZonationBuilder, ids
from pkpdsim.zonation.flow import ImposedFlow
from pkpdsim.zonation.geometry import Layer, bile_flow
from pkpdsim.zonation.profiles import Weights
from pkpdsim.zonation.settings import BileTransport, Geometry, ZonationSettings

from .conftest import mass_g

pytestmark = pytest.mark.filterwarnings("ignore:.*numerical dispersion:UserWarning")

TRACER = Substance("tr", "tracer", Value(180.0, U.g_per_mole))
DUCT = "tr_bile_duct"
"""The amount (mmole) of the tracer in the duct."""


def _settings(
    n: int = 10,
    mode: BileTransport | None = BileTransport.ADVECTIVE,
    mass: float = 1000.0,
    **kwargs: object,
) -> ZonationSettings:
    """The settings of a chain of `n` zones with a bile and no blood flow."""
    return ZonationSettings(
        n_zones=n,
        flow=ImposedFlow(0.0),
        bile=mode,
        liver_mass=mass_g(mass),
        **kwargs,  # ty: ignore[invalid-argument-type]
    )


def _bile_chain(
    tmp_path: Path,
    n: int = 10,
    released_in_zone: int = 10,
    mode: BileTransport | None = BileTransport.ADVECTIVE,
    mass: float = 1000.0,
    **kwargs: object,
) -> roadrunner.RoadRunner:
    """A model with 1 mM of the tracer in the canaliculus of one zone."""
    builder = ZonationBuilder(_settings(n, mode, mass, **kwargs), [TRACER], None)
    path = tmp_path / "bile.xml"
    create_model(
        builder.build(), filepath=path, sbml_level=3, sbml_version=2, validate=False
    )
    r = individual(path)
    r.integrator.absolute_tolerance = 1e-14
    r.integrator.relative_tolerance = 1e-10
    r[f"[{ids.species('tr', Layer.BILE, ids.zone(released_in_zone))}]"] = 1.0
    return r


def _released(r: roadrunner.RoadRunner, zone: int) -> float:
    """The amount (mmole) released: 1 mM in the canaliculus of the zone."""
    return 1.0 * r[ids.volume(Layer.BILE, ids.zone(zone))]


def _bile_ids(r: roadrunner.RoadRunner) -> list[str]:
    """The canalicular species of the tracer in the order of the zones."""
    n = sum(1 for s in r.model.getFloatingSpeciesIds() if s.startswith("tr_bile_z"))
    return [ids.species("tr", Layer.BILE, ids.zone(k)) for k in range(1, n + 1)]


def _in_bile(r: roadrunner.RoadRunner) -> float:
    """The amount (mmole) of the tracer in the canaliculi of all zones."""
    return sum(
        r[f"[{sid}]"] * r[f"V_{sid.removeprefix('tr_')}"] for sid in _bile_ids(r)
    )


def _arrival_time(out: np.ndarray, column: int, threshold: float = 1e-3) -> float:
    """The first time a concentration exceeds the threshold, infinity if never."""
    above = np.nonzero(out[:, column] > threshold)[0]
    return float(out[above[0], 0]) if above.size else math.inf


def test_bile_runs_pericentral_to_periportal(tmp_path: Path) -> None:
    """A release in z10 reaches z9 before z1: the opposite direction arrives never."""
    r = _bile_chain(tmp_path, n=10, released_in_zone=10)
    out = r.simulate(0, 60, 6001, selections=["time", "[tr_bile_z9]", "[tr_bile_z1]"])
    z9, z1 = _arrival_time(out, 1), _arrival_time(out, 2)
    assert math.isfinite(z9)
    assert math.isfinite(z1)
    assert z1 > z9


def test_nothing_runs_towards_the_pericentral_end(tmp_path: Path) -> None:
    """A release in z5 never reaches z6 to z10 (it is the other direction)."""
    r = _bile_chain(tmp_path, n=10, released_in_zone=5)
    names = [f"[tr_bile_z{k}]" for k in range(6, 11)]
    out = r.simulate(0, 60, 601, selections=["time", *names, "[tr_bile_z1]"])
    assert out[:, 1:6].max() == 0.0
    assert out[:, 6].max() > 1e-3  # while z1, downstream, does receive it


@pytest.mark.parametrize("zone", [10, 5, 1])
@pytest.mark.parametrize("mode", [BileTransport.ADVECTIVE, BileTransport.DIFFUSIVE])
def test_the_duct_collects_everything_released(
    tmp_path: Path, zone: int, mode: BileTransport
) -> None:
    """Every release arrives in the duct, from the pericentral end, the middle and z1."""
    r = _bile_chain(tmp_path, n=10, released_in_zone=zone, mode=mode)
    released = _released(r, zone)
    assert released > 0.0
    # the diffusive chain empties in minutes (L**2 / D of 300 um is 2.5 min)
    r.simulate(0, 600, 101)
    assert r[DUCT] == pytest.approx(released, rel=1e-3)
    assert _in_bile(r) < 1e-3 * released


def test_the_amount_is_conserved_while_it_runs(tmp_path: Path) -> None:
    """The bile and the duct hold the release at every time, the duct only gains."""
    r = _bile_chain(tmp_path, n=10, released_in_zone=7)
    released = _released(r, 7)
    out = r.simulate(0, 3, 31, selections=["time", DUCT])
    assert 0.0 < out[-1, 1] < released  # still on its way at 3 min
    assert np.all(np.diff(out[:, 1]) >= 0.0)
    assert out[-1, 1] + _in_bile(r) == pytest.approx(released, rel=1e-9)


@pytest.mark.parametrize("mass", [1000.0, 2000.0, 7.229])
def test_the_transit_time_is_the_bile_volume_over_the_flow(
    tmp_path: Path, mass: float
) -> None:
    """The mean transit from z10 to the duct is the bile volume over `Q_bile`, at any mass.

    10 tanks of a tenth of the bile volume each in series: per kg of liver
    3.738 ml (0.4 % of the 0.9346 l at 1.07 kg/l) over 2.1607 ml/min (0.01562
    ml/min of 7.229 g of rat liver, times 1000 g / 7.229 g) is 1.730 min. Both
    the volume and the flow are proportional to the mass, so the time is the
    same for 1000 g, 2000 g and the 7.229 g of the reference liver. The mean
    time is the area above the curve of the duct amount.
    """
    r = _bile_chain(tmp_path, n=10, released_in_zone=10, mass=mass)
    released = _released(r, 10)
    # the equality: the bile volume over the flow, read from the model
    volume = sum(r[ids.volume(Layer.BILE, ids.zone(k))] for k in range(1, 11))
    assert volume / r["Q_bile"] == pytest.approx(1.730, rel=2e-3)
    out = r.simulate(0, 40, 40001, selections=["time", DUCT])
    mean = np.trapezoid(1.0 - out[:, 1] / released, out[:, 0])
    assert mean == pytest.approx(1.730, rel=2e-3)


def test_the_bile_flow_follows_the_liver_mass(tmp_path: Path) -> None:
    """`Q_bile` is the 0.01562 ml/min of the 7.229 g rat liver per gram of liver."""
    small = _bile_chain(tmp_path, mass=1000.0)
    big = _bile_chain(tmp_path, mass=2000.0)
    # 0.01562 ml/min / 7.229 g = 2.1607e-3 ml/min per g
    assert small["Q_bile"] == pytest.approx(2.1607e-3 * 1000.0 / 1000.0, rel=1e-4)
    assert big["Q_bile"] == pytest.approx(2.0 * small["Q_bile"], rel=1e-12)


def test_the_bile_flow_of_the_species_are_not_scaled_from_each_other() -> None:
    """The human canalicular flow is 0.27 ml/min of 1.8 kg, the rat 0.01562 ml/min of 7.229 g."""
    human = bile_flow(Species.HUMAN).value  # l/min of the reference liver
    rat = bile_flow(Species.RAT).value
    assert human == pytest.approx(0.27e-3, rel=1e-9)
    assert rat == pytest.approx(0.01562e-3, rel=1e-9)
    # per gram of liver the rat secretes 14.4 times more (2.16 against 0.15 ul/min):
    # an upper bound on the species difference, a rat total (with a ductular part)
    # against a human canalicular value; against the human total of 0.43 ml/min
    # (0.24 ul/min/g) it is about 9
    assert (rat / 7.229) / (human / 1800.0) == pytest.approx(14.4, rel=1e-2)


def test_the_diffusive_mode_needs_no_bile_flow(tmp_path: Path) -> None:
    """The diffusive chain has no `Q_bile`, and still moves the tracer into the duct."""
    r = _bile_chain(tmp_path, n=10, released_in_zone=10, mode=BileTransport.DIFFUSIVE)
    assert "Q_bile" not in r.model.getGlobalParameterIds()
    r.simulate(0, 600, 601)
    assert r[DUCT] > 0.0


def test_diffusion_spreads_both_ways_and_flow_only_one(tmp_path: Path) -> None:
    """A release in z5 reaches z10 by diffusion and never by flow."""
    names = ["time", "[tr_bile_z10]"]
    diffusive = _bile_chain(
        tmp_path, released_in_zone=5, mode=BileTransport.DIFFUSIVE
    ).simulate(0, 600, 601, selections=names)
    assert diffusive[:, 1].max() > 1e-4
    advective = _bile_chain(tmp_path, released_in_zone=5).simulate(
        0, 600, 601, selections=names
    )
    assert advective[:, 1].max() == 0.0


def test_the_first_conductance_of_the_diffusive_chain(tmp_path: Path) -> None:
    """From z1 the duct receives `2 D A / dx * c` (the duct is half a zone away).

    Independent of the builder: `A = f_bile V_par / L` of 1 kg of liver is 0.004
    * 0.9346e-3 m3 / 300e-6 m = 0.01246 m2 (the cross-section of all canaliculi),
    `dx = 30e-6 m`, and the aqueous diffusion of a 180 g/mole solute is read
    from the model (`D_tr`, m2/s). The rate constant `2 D A / dx / V_z1` is
    about 90 per min, so the test reads the first 6 ms.
    """
    r = _bile_chain(tmp_path, n=10, released_in_zone=1, mode=BileTransport.DIFFUSIVE)
    d = r["D_tr"]
    assert 1e-10 < d < 1e-9  # the order of a small solute in water
    area = 0.004 * (1.0 / 1.07 * 1e-3) / 300e-6
    k = 2.0 * d * area / 30e-6 * 1e3 * 60.0  # l/min
    # the amount which left in a short time, the neighbour has not given back yet
    t = 1e-4
    out = r.simulate(0, t, 11, selections=["time", DUCT])
    assert out[-1, 1] == pytest.approx(1.0 * k * t, rel=0.02)


def test_a_model_without_bile_emits_no_bile_and_simulates(tmp_path: Path) -> None:
    """`bile=None`: no bile compartment, reaction or duct, the model runs."""
    builder = ZonationBuilder(_settings(4, None), [TRACER], None)
    assert builder.bile() == []
    path = tmp_path / "nobile.xml"
    create_model(
        builder.build(), filepath=path, sbml_level=3, sbml_version=2, validate=False
    )
    r = individual(path)
    r[f"[{ids.species('tr', Layer.PLASMA, 'z1')}]"] = 1.0
    r.simulate(0, 10, 11)
    names = r.model.getFloatingSpeciesIds() + r.model.getCompartmentIds()
    assert not [s for s in names if "bile" in s]
    assert not [s for s in r.model.getReactionIds() if "bile" in s]


@pytest.mark.parametrize(
    "geometry", [Geometry.SINUSOID, Geometry.WELL_STIRRED], ids=str
)
def test_both_geometries_have_the_chain(tmp_path: Path, geometry: Geometry) -> None:
    """The well-stirred zones keep the canalicular chain."""
    r = _bile_chain(tmp_path, n=4, released_in_zone=4, geometry=geometry)
    released = _released(r, 4)
    r.simulate(0, 600, 101)
    assert r[DUCT] == pytest.approx(released, rel=1e-3)


def test_a_profile_of_the_canalicular_flow_scales_the_zones(tmp_path: Path) -> None:
    """The canalicular flow can be zonated: a slower zone holds the tracer longer."""
    flat = _bile_chain(tmp_path, n=4, released_in_zone=4)
    slow = _bile_chain(
        tmp_path,
        n=4,
        released_in_zone=4,
        profiles={"canalicular_flow": Weights((1.0, 1.0, 1.0, 0.2))},
    )
    # the levels are normalized to the mean 0.8: 1.25, 1.25, 1.25, 0.25
    assert slow["p_canalicular_flow_z4"] == pytest.approx(0.25)
    assert slow["p_canalicular_flow_z3"] == pytest.approx(1.25)
    # at the start only z4 holds the tracer: its flow is the quarter of the flat one
    assert flat["FLOW_tr_bile_z4_bile_z3"] > 0.0
    assert slow["FLOW_tr_bile_z4_bile_z3"] == pytest.approx(
        0.25 * flat["FLOW_tr_bile_z4_bile_z3"], rel=1e-9
    )
    slow.simulate(0, 600, 101)
    assert slow[DUCT] == pytest.approx(_released(slow, 4), rel=1e-3)


@pytest.mark.parametrize(
    ("mode", "process"),
    [
        (BileTransport.ADVECTIVE, "canalicular_diffusion"),
        (BileTransport.DIFFUSIVE, "canalicular_flow"),
        (None, "canalicular_flow"),
    ],
)
def test_a_profile_of_a_process_the_bile_does_not_have_is_refused(
    mode: BileTransport | None, process: str
) -> None:
    """A profile of the diffusion of a bile with a flow (and the reverse) names no process."""
    settings = _settings(4, mode, profiles={process: Weights((1.0, 1.0, 1.0, 2.0))})
    with pytest.raises(ValueError, match="does not have"):
        ZonationBuilder(settings, [TRACER], None)
