"""The flow through the sinusoids: imposed, or driven by a pressure drop.

The pressure-driven flow is the leaky tube of the spt model: the blood flows
through a tube with the hydraulic resistance `W` of the lumen which leaks
through the fenestrae with the resistance `w` of all pores in parallel. Its
flow changes with the position along the tube, so every zone has the flow of
its own position. The tests compute the two limits independently of the
implementation: the Poiseuille flow of a tube without leak, and the closed
form of the spt model (`spt_liver_zonated.py` of pkdb_models, in the
exponential form, where the engine uses hyperbolic functions).
"""

import dataclasses
import itertools
import math
from collections.abc import Sequence
from pathlib import Path

import pytest
import roadrunner
from pbpk.conftest import individual
from sbmlutils.factory import Parameter, create_model

from pkpdsim.pbpk.substance import Substance
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Value
from pkpdsim.species import Species
from pkpdsim.zonation import ZonationBuilder, ids
from pkpdsim.zonation.flow import (
    ImposedFlow,
    PressureFlow,
    calibrated_resistance_factor,
    flow_objects,
    resolved_resistance_factor,
    zone_flows,
)
from pkpdsim.zonation.geometry import SinusoidGeometry, geometry, reference_flow
from pkpdsim.zonation.profiles import Linear, zone_positions
from pkpdsim.zonation.settings import ZonationSettings

from .conftest import mass_g

pytestmark = pytest.mark.filterwarnings("ignore:.*numerical dispersion:UserWarning")

TRACER = Substance("tr", "tracer", Value(100.0, U.g_per_mole))
_MASS = mass_g(1000.0)
_N = 10
_SPT_FLOW = PressureFlow(resistance_factor=10.0)
"""The uncalibrated factor of the spt model, which its formula is checked at."""


def _settings(flow: ImposedFlow | PressureFlow, **kwargs: object) -> ZonationSettings:
    """The settings of a chain of ten zones of a 1 kg liver."""
    return ZonationSettings(
        n_zones=_N,
        flow=flow,
        liver_mass=_MASS,
        **kwargs,  # ty: ignore[invalid-argument-type]
    )


def _build(settings: ZonationSettings, directory: Path) -> roadrunner.RoadRunner:
    """The roadrunner model of the settings, loaded as written."""
    directory.mkdir(parents=True, exist_ok=True)
    model = ZonationBuilder(settings, [TRACER], None).build()
    path = directory / "flow.xml"
    create_model(model, filepath=path, sbml_level=3, sbml_version=2, validate=False)
    return roadrunner.RoadRunner(str(path))


def _simulated(
    flow: ImposedFlow | PressureFlow,
    tmp_path: Path,
    geo: SinusoidGeometry | None = None,
    **changes: float,
) -> roadrunner.RoadRunner:
    """The initialized roadrunner model of a chain with the given flow."""
    settings = _settings(flow, geometry_constants=geo)
    model = ZonationBuilder(settings, [TRACER], None).build()
    path = tmp_path / "flow.xml"
    create_model(model, filepath=path, sbml_level=3, sbml_version=2, validate=False)
    if changes:
        return individual(path, **changes)
    # roadrunner evaluates the initial assignments of the file as written; the
    # values of the file are what `individual` keeps for an unchanged parameter,
    # so the formulas are only tested by loading the file itself
    return roadrunner.RoadRunner(str(path))


def _q_sin(r: roadrunner.RoadRunner) -> list[float]:
    """The flow (m**3/s) of one sinusoid in every zone."""
    return [r[ids.sinusoid_flow(ids.zone(k))] for k in range(1, _N + 1)]


def _without_leak(geo: SinusoidGeometry) -> SinusoidGeometry:
    """The geometry with closed fenestrae."""
    return dataclasses.replace(geo, porosity_pp=0.0, porosity_pc=0.0)


def _poiseuille(geo: SinusoidGeometry, dp: float, nu: float) -> float:
    """The flow (m**3/s) of a tube of radius `r_sin` and length `L_sin`: `dp * pi * r**4 / (8 nu L)`."""
    return dp * math.pi * geo.r_sin**4 / (8.0 * nu * geo.L_sin)


def _spt_q(position: float) -> float:
    """The flow (m**3/s) of the spt model at a position, from its own formula.

    The parameters, the expressions `W`, `w`, `lambda` and `PP_Q` are those of
    the spt model (`spt_liver_zonated.py` lines 346 to 382), written in the
    exponential form of that file.
    """
    pa, pb = 1333.22, 266.64
    nu = 10.0 * 0.0018
    y_sin, y_end, r_fen, n_fen, length = 4.4e-6, 0.165e-6, 53.5e-9, 10e12, 500e-6
    p0 = 0.5 * (pa + pb)
    w_cap = 8.0 * nu / (math.pi * y_sin**4)
    w_pore = 4.0 * nu * y_end / (math.pi**2 * r_fen**4 * y_sin * n_fen)
    lam = math.sqrt(w_pore / w_cap)
    x = position * length
    e_neg, e_pos = math.exp(-length / lam), math.exp(length / lam)
    return (
        -1.0
        / math.sqrt(w_cap * w_pore)
        * (
            (-(pb - p0) + (pa - p0) * e_neg) / (e_neg - e_pos) * math.exp(x / lam)
            - (pb - p0 - (pa - p0) * e_pos) / (e_neg - e_pos) * math.exp(-x / lam)
        )
    )


def _spt_geometry() -> SinusoidGeometry:
    """The rat geometry with the parameters of the spt model.

    The count of the fenestrae is derived from the porosity, so the porosity
    which gives the 1e13 /m**2 of spt is `n_fen * pi * r_fen**2`.
    """
    porosity = 10e12 * math.pi * (53.5e-9) ** 2
    return dataclasses.replace(
        geometry(ZonationSettings().species),
        r_sin=4.4e-6,
        L_sin=500e-6,
        t_endo=0.165e-6,
        r_fen=53.5e-9,
        porosity_pp=porosity,
        porosity_pc=porosity,
    )


def _eq(expected: float | list[float], rel: float = 1e-9) -> object:
    """A relative comparison without an absolute tolerance.

    The flows of one sinusoid are of the order of 1e-14 m**3/s, below the
    default absolute tolerance of `pytest.approx` (1e-12), which would accept
    any flow: only the relative tolerance may decide.
    """
    return pytest.approx(expected, rel=rel, abs=0.0)


def _value_of(objects: Sequence[object], sid: str) -> float:
    """The value of the parameter `sid` among the emitted objects."""
    for o in objects:
        if isinstance(o, Parameter) and o.sid == sid:
            assert isinstance(o.value, float)
            return o.value
    raise AssertionError(
        f"no parameter {sid} in {[getattr(o, 'sid', o) for o in objects]}"
    )


def test_imposed_flow_emits_q_directly() -> None:
    """The flow of the liver is the number given, in l/min."""
    settings = _settings(ImposedFlow(Value(1.0, "ml/min")))
    objects = flow_objects(settings, geometry(settings.species), _MASS)
    assert _value_of(objects, ids.FLOW) == _eq(1e-3)


def test_an_imposed_flow_per_mass_follows_the_liver_mass() -> None:
    """The perfusion states its flow per g of liver (ml/min/g)."""
    flow = ImposedFlow(per_mass=Value(2.0, "ml/min/g"))
    settings = _settings(flow)
    objects = flow_objects(settings, geometry(settings.species), mass_g(10.0))
    assert _value_of(objects, ids.FLOW) == _eq(0.02)


def test_a_flow_and_a_flow_per_mass_exclude_each_other() -> None:
    """Which of the two would win is not for the model to guess."""
    with pytest.raises(ValueError, match="per_mass"):
        ImposedFlow(Value(1.0, "ml/min"), Value(1.0, "ml/min/g"))


def test_the_leaky_tube_reduces_to_poiseuille_without_fenestrae(tmp_path: Path) -> None:
    """Closed fenestrae (porosity 0, so a count of 0) leave the flow of a tube.

    The formula is singular at a count of 0 (`w` diverges, so does `lambda`)
    and its limit is the linear pressure drop of Poiseuille, a flow without
    position. The engine takes that limit exactly at a porosity of 0, and the
    next test approaches it with a vanishing porosity through the formula.
    """
    geo = _without_leak(geometry(ZonationSettings().species))
    flow = PressureFlow(resistance_factor=10.0)
    expected = _poiseuille(geo, 1333.22 - 266.64, 10.0 * 0.0018)
    for q in _q_sin(_simulated(flow, tmp_path, geo)):
        assert q == _eq(expected, rel=1e-9)


def test_the_flow_approaches_poiseuille_as_the_leak_vanishes(tmp_path: Path) -> None:
    """The exponential solution tends to Poiseuille as `lambda` outgrows the tube."""
    base = geometry(ZonationSettings().species)
    expected = _poiseuille(base, 1333.22 - 266.64, 10.0 * 0.0018)
    deviation = []
    for porosity in (1e-4, 1e-6, 1e-8):
        geo = dataclasses.replace(base, porosity_pp=porosity, porosity_pc=porosity)
        flows = _q_sin(_simulated(_SPT_FLOW, tmp_path, geo))
        deviation.append(max(abs(q / expected - 1.0) for q in flows))
    assert deviation[-1] < 1e-6
    assert deviation[0] > 10.0 * deviation[1] > 100.0 * deviation[2] > 0.0
    # and the real porosity is far from it: the leak matters
    leaky = _q_sin(_simulated(_SPT_FLOW, tmp_path))
    assert max(abs(q / expected - 1.0) for q in leaky) > 0.05


def test_the_leaky_tube_matches_the_spt_formula_at_the_spt_parameters(
    tmp_path: Path,
) -> None:
    """Every zone has the spt flow at its own position, `PP_x = position * L`."""
    flows = _q_sin(_simulated(_SPT_FLOW, tmp_path, _spt_geometry()))
    for q, position in zip(flows, zone_positions(_N), strict=True):
        assert q == _eq(_spt_q(position), rel=1e-6)


def test_the_values_of_the_file_are_the_values_of_its_formulas(tmp_path: Path) -> None:
    """What the build writes as a value is what the initial assignment evaluates to.

    A loader which keeps the values of the file (as `individual` does for an
    unchanged parameter) and one which evaluates the formulas give one flow.
    """
    settings = _settings(PressureFlow())
    model = ZonationBuilder(settings, [TRACER], None).build()
    path = tmp_path / "flow.xml"
    create_model(model, filepath=path, sbml_level=3, sbml_version=2, validate=False)
    formulas = roadrunner.RoadRunner(str(path))
    values = individual(path)
    for k in range(1, _N + 1):
        for sid in (ids.sinusoid_flow(ids.zone(k)), ids.zone_flow(None, ids.zone(k))):
            assert values[sid] == _eq(formulas[sid], rel=1e-12)
    assert values[ids.FLOW] == _eq(formulas[ids.FLOW], rel=1e-12)


def test_the_zones_have_different_flows(tmp_path: Path) -> None:
    """The flow depends on the position (the FIXME of spt): the zones differ."""
    flows = _q_sin(_simulated(PressureFlow(), tmp_path))
    assert max(flows) / min(flows) > 1.01


def test_the_default_oncotic_pressure_is_the_mean_of_the_pressures(
    tmp_path: Path,
) -> None:
    """Without a net filtration the flow is symmetric along the tube.

    With `P0 = (p_pp + p_pc) / 2` what leaks out of the first half returns in
    the second half, so the first and the last zone carry the same flow. A
    lower oncotic pressure filters out and the pericentral zone carries less.
    """
    flows = _q_sin(_simulated(PressureFlow(), tmp_path))
    assert flows[0] == _eq(flows[-1], rel=1e-9)
    low = _q_sin(_simulated(PressureFlow(p_onc=Value(0.0, "Pa")), tmp_path))
    assert low[-1] < 0.9 * low[0]


def test_a_higher_pressure_drop_gives_a_higher_flow(tmp_path: Path) -> None:
    """The flow follows the portal pressure, at build and before a simulation."""
    base = _q_sin(_simulated(PressureFlow(), tmp_path))
    built = _q_sin(_simulated(PressureFlow(p_pp=Value(2000.0, "Pa")), tmp_path))
    changed = _q_sin(_simulated(PressureFlow(), tmp_path, P_pp=2000.0))
    assert all(b > a for a, b in zip(base, built, strict=True))
    assert changed == _eq(built, rel=1e-9)


def test_the_flow_of_the_liver_follows_the_mass(tmp_path: Path) -> None:
    """A liver of twice the mass has twice the sinusoids and twice the flow."""
    r = _simulated(PressureFlow(), tmp_path)
    q = r[ids.zone_flow(None, "z4")]
    twice = _simulated(PressureFlow(), tmp_path, M_li=2.0 * r[ids.MASS])
    assert twice[ids.zone_flow(None, "z4")] == _eq(2.0 * q, rel=1e-9)
    assert twice[ids.sinusoid_flow("z4")] == _eq(r[ids.sinusoid_flow("z4")])


def test_the_convection_reads_the_flow_of_its_zone(tmp_path: Path) -> None:
    """The reaction out of a zone is its own flow times its concentration."""
    r = _simulated(PressureFlow(), tmp_path)
    r["[tr_plasma_z3]"] = 2.0
    r["[tr_plasma_z5]"] = 2.0
    for k in (3, 5):
        zone = ids.zone(k)
        rate = r[f"FLOW_tr_plasma_{zone}_plasma_z{k + 1}"]
        assert rate == _eq(2.0 * r[ids.zone_flow("plasma", zone)], rel=1e-9)
    assert r[ids.zone_flow("plasma", "z3")] != _eq(
        r[ids.zone_flow("plasma", "z5")], rel=1e-3
    )


def test_the_zone_flows_split_into_the_shares_of_the_carriers(tmp_path: Path) -> None:
    """Plasma and erythrocytes carry `1 - Hct` and `Hct` of the flow of the zone."""
    r = _simulated(PressureFlow(), tmp_path)
    q = r[ids.zone_flow(None, "z5")]
    h = r[ids.HAEMATOCRIT]
    assert r[ids.zone_flow("plasma", "z5")] == _eq(q * (1.0 - h))
    assert r[ids.zone_flow("rbc", "z5")] == _eq(q * h)


def test_the_flow_of_the_liver_is_the_mean_of_the_zones(tmp_path: Path) -> None:
    """`Q_li` of the pressure flow is the mean flow of the zones, in l/min."""
    r = _simulated(PressureFlow(), tmp_path)
    zones = [r[ids.zone_flow(None, ids.zone(k))] for k in range(1, _N + 1)]
    assert r[ids.FLOW] == _eq(sum(zones) / _N)
    assert r[ids.FLOW] > 0.0


def test_the_flow_of_the_liver_is_in_litre_per_minute(tmp_path: Path) -> None:
    """The flow of the whole liver is the sinusoid flow times the sinusoids."""
    geo = geometry(ZonationSettings().species)
    r = _simulated(PressureFlow(), tmp_path)
    lumen = (r[ids.volume("plasma", "z1")] + r[ids.volume("rbc", "z1")]) * _N * 1e-3
    n_sin = lumen / (math.pi * geo.r_sin**2 * geo.L_sin)
    assert r[ids.zone_flow(None, "z1")] == _eq(
        r[ids.sinusoid_flow("z1")] * n_sin * 1e3 * 60.0, rel=1e-9
    )


def test_a_negative_pressure_drop_is_refused() -> None:
    """The chain does not run backwards; the message names both pressures."""
    flow = PressureFlow(p_pp=Value(100.0, "Pa"), p_pc=Value(900.0, "Pa"))
    with pytest.raises(ValueError, match=r"100.*900|900.*100"):
        ZonationBuilder(_settings(flow), [TRACER], None).build()


def test_the_pressures_take_the_units_they_are_given_in(tmp_path: Path) -> None:
    """10 mmHg and 2 mmHg are 1333.22 Pa and 266.64 Pa."""
    mmhg = PressureFlow(p_pp=Value(10.0, "mmHg"), p_pc=Value(2.0, "mmHg"))
    a = _q_sin(_simulated(mmhg, tmp_path))
    b = _q_sin(_simulated(PressureFlow(), tmp_path))
    assert a == _eq(b, rel=1e-5)


def test_the_table_count_of_fenestrae_agrees_with_the_measured_porosity() -> None:
    """The table's `n_fen` is a cross-check of the count derived from the porosity.

    The count per zone is `porosity / (pi r_fen**2)`, 6.7e12 periportal and
    8.9e12 pericentral at the measured 6 % and 8 % porosity (Braet and Wisse
    2002); the 1e13 of the table must lie within 35 % of that range.
    """
    geo = geometry(ZonationSettings().species)
    derived = [p / (math.pi * geo.r_fen**2) for p in (geo.porosity_pp, geo.porosity_pc)]
    assert derived[0] == pytest.approx(6.67e12, rel=1e-2)
    assert derived[1] == pytest.approx(8.90e12, rel=1e-2)
    assert min(derived) / 1.35 <= geo.n_fen <= max(derived) * 1.35
    assert geo.n_fen == pytest.approx(sum(derived) / 2.0, rel=0.35)


def test_a_formula_flow_is_refused() -> None:
    """A flow which is a formula would have no number to check."""
    from pkpdsim.pbpk.values import Formula

    flow = ImposedFlow(Formula("1 + 1", U.l_per_min))
    with pytest.raises(ValueError, match="PressureFlow"):
        ZonationBuilder(ZonationSettings(flow=flow), [TRACER], None).build()


def _inflow(r: roadrunner.RoadRunner) -> float:
    """The periportal inflow (l/min) of a pressure flow from the values of its model.

    The flow of the leaky tube at `x = 0` from `W_sin`, `lambda_sin`, the
    pressures and the number of sinusoids of the loaded model, in the closed
    form `(a cosh(l) - b) / (W lambda sinh(l))`.
    """
    pop = r[ids.PRESSURE_ONCOTIC]
    a, b = r[ids.PRESSURE_PORTAL] - pop, r[ids.PRESSURE_CENTRAL] - pop
    big = r[ids.LENGTH] / r[ids.DECAY_LENGTH]
    q = (a * math.cosh(big) - b) / (
        r[ids.RESISTANCE_LUMEN] * r[ids.DECAY_LENGTH] * math.sinh(big)
    )
    return q * r[ids.N_SINUSOIDS] * 1e3 * 60.0


@pytest.mark.parametrize("species", [Species.RAT, Species.HUMAN])
def test_both_flow_modes_agree_at_their_defaults(
    species: Species, tmp_path: Path
) -> None:
    """The periportal inflow of a default pressure flow is the default imposed flow."""
    outputs = {}
    for name, flow in (("pressure", PressureFlow()), ("imposed", ImposedFlow())):
        settings = ZonationSettings(n_zones=_N, flow=flow, species=species)
        model = ZonationBuilder(settings, [TRACER], None).build()
        path = tmp_path / f"{name}.xml"
        create_model(model, filepath=path, sbml_level=3, sbml_version=2, validate=False)
        outputs[name] = roadrunner.RoadRunner(str(path))
    expected = reference_flow(species).value
    assert outputs["imposed"][ids.FLOW] == _eq(expected, 1e-12)
    assert _inflow(outputs["pressure"]) == _eq(expected, 1e-9)


def test_the_calibrated_factor_is_per_species_and_replaced_by_a_given_one() -> None:
    """The rat needs about five times the factor of spt; a given factor is kept."""
    rat, human = (
        calibrated_resistance_factor(sp, geometry(sp))
        for sp in (Species.RAT, Species.HUMAN)
    )
    assert 45.0 < rat < 60.0
    assert human != _eq(rat, 0.1)
    assert (
        resolved_resistance_factor(
            PressureFlow(resistance_factor=7.0), geometry(Species.RAT), Species.RAT
        )
        == 7.0
    )


def test_a_porosity_profile_changes_the_sieving_and_not_the_hydraulics(
    tmp_path: Path,
) -> None:
    """The leaky tube is evaluated at the mean porosity of the table.

    The analytic solution assumes a uniform leak, so a zonal `porosity`
    profile scales the fenestral permeability of Task 6 and leaves the flows
    of the zones as they are.
    """
    profile = {"porosity": Linear(pp=2.0, pc=0.5)}
    plain, scaled = (
        _build(_settings(PressureFlow(), **extra), tmp_path / name)
        for name, extra in (("a", {}), ("b", {"profiles": profile}))
    )
    assert _q_sin(scaled) == _eq(_q_sin(plain), 1e-12)
    for r in (plain, scaled):
        r["[tr_plasma_z1]"] = 1.0
    ratio = scaled["FEN_tr_z1"] / plain["FEN_tr_z1"]
    assert ratio == _eq(scaled["p_porosity_z1"], 1e-9)
    assert abs(ratio - 1.0) > 0.5


def test_the_flow_falls_to_the_middle_and_rises_again() -> None:
    """With `P0` at the mean of the pressures the flow is not monotonic.

    Fluid filters out of the first half and is reabsorbed in the second.
    """
    settings = _settings(PressureFlow())
    geo = geometry(settings.species)
    q = zone_flows(settings, geo, _MASS)
    half = _N // 2
    assert all(a > b for a, b in itertools.pairwise(q[:half]))
    assert all(a < b for a, b in itertools.pairwise(q[half:]))
