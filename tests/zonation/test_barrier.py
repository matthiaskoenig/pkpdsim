"""The fenestrated barrier and the radial exchange of the zones.

The tests use a tracer, a substance which does not react, so no cell model is
needed. The tracer starts in the plasma of every zone and the rates of the
radial reactions are read at the start.
"""

import dataclasses
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
from pkpdsim.zonation.flow import ImposedFlow, PressureFlow
from pkpdsim.zonation.geometry import Layer, geometry
from pkpdsim.zonation.profiles import Weights
from pkpdsim.zonation.settings import Geometry, ZonationSettings
from pkpdsim.zonation.transport import M2_PER_S, Solute

from .conftest import liver_flow, mass_g

TRACER = Substance("tr", "tracer", Value(180.0, U.g_per_mole))
"""A small solute, much smaller than a fenestra."""


class _Sized(ZonationBuilder):
    """A builder whose tracer has a given molecular radius (m), 0 for derived."""

    radius = 0.0

    def _radius(self, substance: Solute) -> Value:
        """Get the radius of the tracer."""
        if not self.radius:
            return super()._radius(substance)
        return Value(self.radius, U.meter)


def _builder(
    n: int = 4,
    flow: float = 1.0,
    radius: float = 0.0,
    porosity_pp: float = 0.06,
    porosity_pc: float = 0.08,
    **kwargs: object,
) -> _Sized:
    """A builder of a tracer chain with the given porosity gradient."""
    geo = dataclasses.replace(
        geometry(Species.RAT), porosity_pp=porosity_pp, porosity_pc=porosity_pc
    )
    kwargs.setdefault("reference_flow", liver_flow())
    settings = ZonationSettings(
        n_zones=n,
        flow=ImposedFlow(flow),
        liver_mass=mass_g(1000.0),
        geometry_constants=geo,
        **kwargs,  # ty: ignore[invalid-argument-type]
    )
    builder = _Sized(settings, [TRACER], None)
    builder.radius = radius
    return builder


def _barrier(
    tmp_path: Path,
    endothelial_ps: float = 0.0,
    basolateral_ps: float = 0.0,
    porosity: float | None = None,
    builder: ZonationBuilder | None = None,
    **kwargs: object,
) -> roadrunner.RoadRunner:
    """A tracer in the plasma of every zone, with the transcellular paths set."""
    if porosity is not None:
        kwargs |= {"porosity_pp": porosity, "porosity_pc": porosity}
    builder = builder or _builder(**kwargs)  # ty: ignore[invalid-argument-type]
    path = tmp_path / "barrier.xml"
    create_model(
        builder.build(), filepath=path, sbml_level=3, sbml_version=2, validate=False
    )
    r = individual(
        path,
        P_endo_in_tr=endothelial_ps,
        P_endo_out_tr=endothelial_ps,
        P_hep_tr=basolateral_ps,
    )
    r.integrator.absolute_tolerance = 1e-14
    r.integrator.relative_tolerance = 1e-10
    for k in range(1, builder.settings.n_zones + 1):
        r[f"[{ids.species('tr', Layer.PLASMA, ids.zone(k))}]"] = 1.0
    return r


def _fenestral_fluxes(r: roadrunner.RoadRunner) -> list[float]:
    """The rate of the fenestral path (mmole/min) of every zone."""
    n = sum(1 for sid in r.model.getReactionIds() if sid.startswith("FEN_tr_"))
    return [r[f"FEN_tr_{ids.zone(k)}"] for k in range(1, n + 1)]


def _flux(r: roadrunner.RoadRunner) -> float:
    """The rate of both paths from the plasma into the Disse space, first zone."""
    return r["FEN_tr_z1"] + r["END_tr_z1_in"]


def _arrived(r: roadrunner.RoadRunner, horizon: float) -> float:
    """The amount (mmole) in the space of Disse of the first zone after `horizon` min."""
    r.simulate(0, horizon, 101)
    return r["[tr_disse_z1]"] * r["V_disse_z1"]


def test_the_two_paths_are_parallel_and_additive(tmp_path: Path) -> None:
    """The amount arriving in the Disse space with both paths is about the sum.

    Fenestral path alone, cell path alone (a 60 nm solute has no fenestral
    reaction; the cell path is the two membranes in series) and both, over
    0.01 min with the Disse concentration far below the plasma one. Both is
    strictly greater than either and within 10 % below the sum. It is not
    exact: the paths share the plasma they draw from and the Disse pool they
    fill (which pushes back on both), and the cell path drains the Disse pool
    through its reversible abluminal membrane while its intermediate pool, the
    endothelial cell, fills. The buffering of that pool is pinned by
    `test_the_endothelial_cell_buffers_the_cell_path_early`, the same mechanism
    at a weaker cell path.
    """
    horizon = 0.01
    ps = 1e-7
    fen = _arrived(_barrier(tmp_path, porosity=1e-5), horizon)
    cell_dir = tmp_path / "cell"
    cell_dir.mkdir()
    cell = _arrived(
        _barrier(cell_dir, endothelial_ps=ps, builder=_builder(radius=60e-9)),
        horizon,
    )
    both_dir = tmp_path / "both"
    both_dir.mkdir()
    both = _arrived(_barrier(both_dir, endothelial_ps=ps, porosity=1e-5), horizon)
    assert both > max(fen, cell)
    assert 0.9 * (fen + cell) <= both <= fen + cell


def test_the_endothelial_cell_buffers_the_cell_path_early(tmp_path: Path) -> None:
    """With a weak cell path both paths deliver less than the fenestrae alone, at first.

    The transcellular path is two conductances in series through a compartment
    with volume, the endothelial cell, so early after the tracer arrives the
    cell is a sink which also drains the Disse space through its reversible
    abluminal membrane. Porosity 1e-5, `ps` 1e-8, flow 1 l/min: the amount in
    the Disse space of zone 1 with both paths is 0.9955 of the fenestral-only
    amount at 0.01 min, and 1.0545 of it at 1 min, when the plasma has washed
    out and the filled cell gives its tracer back to the Disse space.
    Additivity is a quasi-steady-state property, not an instantaneous one.
    """
    fen = _barrier(tmp_path, porosity=1e-5)
    both_dir = tmp_path / "both"
    both_dir.mkdir()
    both = _barrier(both_dir, porosity=1e-5, endothelial_ps=1e-8)
    amounts = {}
    for name, r in (("fen", fen), ("both", both)):
        r.simulate(0, 0.01, 11)
        early = r["[tr_disse_z1]"] * r["V_disse_z1"]
        r.simulate(0.01, 1.0, 101)
        amounts[name] = (early, r["[tr_disse_z1]"] * r["V_disse_z1"])
    assert amounts["both"][0] < amounts["fen"][0]
    assert amounts["both"][1] > amounts["fen"][1]


def test_a_solute_larger_than_the_fenestra_uses_only_the_cell_path(
    tmp_path: Path,
) -> None:
    """There is no fenestral reaction at all, by construction."""
    builder = _builder(radius=60e-9)
    assert not [r for r in builder.radial_exchange() if str(r.sid).startswith("FEN_")]
    r = _barrier(tmp_path, endothelial_ps=1e-6, builder=builder)
    assert "FEN_tr_z1" not in r.model.getReactionIds()
    assert r["END_tr_z1_in"] > 0.0


def test_the_porosity_gradient_is_a_profile(tmp_path: Path) -> None:
    """The fenestral flux rises from the periportal to the pericentral zone."""
    fluxes = _fenestral_fluxes(_barrier(tmp_path, porosity_pp=0.06, porosity_pc=0.08))
    assert fluxes[-1] > fluxes[0]
    assert fluxes == sorted(fluxes)


def test_the_default_gradient_is_the_table_one(tmp_path: Path) -> None:
    """Zone permeabilities follow the interpolation of the table."""
    n = 4
    r = _barrier(tmp_path, n=n)
    geo = geometry(Species.RAT)
    for k in range(1, n + 1):
        eps = geo.porosity_pp + (geo.porosity_pc - geo.porosity_pp) * (k - 0.5) / n
        assert r[f"eps_fen_{ids.zone(k)}"] == pytest.approx(eps)


def test_a_profile_scales_the_gradient(tmp_path: Path) -> None:
    """A mean-1 factor multiplies the gradient of the table."""
    plain = _fenestral_fluxes(_barrier(tmp_path))
    scaled = _fenestral_fluxes(
        _barrier(tmp_path, profiles={"porosity": Weights(levels=(2.0, 1.0, 1.0, 1.0))})
    )
    factor = Weights(levels=(2.0, 1.0, 1.0, 1.0)).values(4)
    for a, b, f in zip(plain, scaled, factor, strict=True):
        assert b == pytest.approx(a * f)


def _mismatch(r: roadrunner.RoadRunner, floor: float = 0.0) -> float:
    """The relative difference of plasma and Disse concentration at 0.5 min.

    No absolute tolerance: the tracer washes out of the lumen, so a late sample
    compares two numbers near 0; `floor` asserts there is still tracer.
    """
    r.simulate(0, 0.5, 51)
    plasma = r["[tr_plasma_z1]"]
    assert plasma > floor, "the tracer has washed out, the test would be vacuous"
    return abs(plasma - r["[tr_disse_z1]"]) / plasma


def test_an_open_barrier_equilibrates_plasma_and_disse(tmp_path: Path) -> None:
    """Open paths equilibrate the layers, a nearly closed barrier does not.

    At 0.5 min with the tracer still in the lumen: the open barrier (porosity 1,
    open endothelial path) is within 1e-3, a barrier of porosity 1e-6 (and
    one a 60 nm solute cannot cross) is far from it. The barrier at the
    measured porosity equilibrates a small solute within 1 %: the fenestrae
    exchange it in seconds, which is why the porosity matters for larger solutes.
    """
    openb = _mismatch(_barrier(tmp_path, endothelial_ps=1e3, porosity=1.0), floor=1e-4)
    assert openb < 1e-3
    for name, kwargs in (
        ("tight", {"porosity": 1e-6}),
        ("large", {"builder": _builder(radius=60e-9)}),
    ):
        d = tmp_path / name
        d.mkdir()
        assert _mismatch(_barrier(d, **kwargs)) > 0.5  # ty: ignore[invalid-argument-type]
    measured = tmp_path / "measured"
    measured.mkdir()
    assert _mismatch(_barrier(measured), floor=1e-4) < 1e-2


def test_the_basolateral_path_moves_the_solute_into_the_hepatocyte(
    tmp_path: Path,
) -> None:
    """The hepatocyte fills when the basolateral path is set, and not otherwise."""
    r = _barrier(tmp_path, basolateral_ps=1e-5, porosity=1.0)
    r.simulate(0, 60, 61)
    assert r["[tr_cell_z1]"] > 0.0
    closed = _barrier(tmp_path, porosity=1.0)
    closed.simulate(0, 60, 61)
    assert closed["[tr_cell_z1]"] == 0.0


def test_the_transcellular_paths_default_to_off(tmp_path: Path) -> None:
    """The model input of the cell paths is 0 until a model sets it."""
    path = tmp_path / "default.xml"
    create_model(
        _builder().build(), filepath=path, sbml_level=3, sbml_version=2, validate=False
    )
    r = individual(path)
    for pid in ("P_endo_in_tr", "P_endo_out_tr", "P_hep_tr"):
        assert r[pid] == 0.0
    assert r["P_fen_tr"] > 0.0


def test_the_fenestral_permeability_by_hand(tmp_path: Path) -> None:
    """A radius of 5.35 nm in a rat fenestra gives 3.886e-3 m/s for D = 1e-9.

    rat: t_endo = 0.165e-6 m, r_fen = 53.5e-9 m, so the ratio is 0.1 and
    renkin(0.1) = 0.81 * 0.7916805 = 0.6412612; the permeability of an open
    pore is 1e-9 / 0.165e-6 * 0.6412612 = 6.060606e-3 * 0.6412612 = 3.8864e-3
    m/s.
    """

    class _Diffusing(_Sized):
        """A tracer of a given diffusion coefficient."""

    substance = Substance(
        "tr", "tracer", Value(180.0, U.g_per_mole), diffusion=Value(1e-9, M2_PER_S)
    )
    builder = _Diffusing(_builder().settings, [substance], None)
    builder.radius = 5.35e-9
    assert _barrier(tmp_path, builder=builder)["P_fen_tr"] == pytest.approx(
        3.8864e-3, rel=1e-4
    )


def test_the_lateral_areas_by_hand(tmp_path: Path) -> None:
    """The areas of zone 1 of 4 in a rat liver of 1 kg without erythrocytes.

    The liver is 1 / 1.07 = 0.934579 l, the plasma 10.6 % of it, a quarter of
    that per zone: 0.0247664 l = 2.47664e-5 m3. With the radius 4.4 um of the
    sinusoid the area of the lumen is 2 * 2.47664e-5 / 4.4e-6 = 11.2574 m2.
    The abluminal face of the endothelial cell lies at 4.4 + 0.165 = 4.565 um,
    the interface of the hepatocyte at 4.565 + 1.0 = 5.565 um: 11.2574 *
    4.565 / 4.4 = 11.6796 m2 and 11.2574 * 5.565 / 4.4 = 14.2381 m2.
    """
    r = _barrier(tmp_path, haematocrit=0.0)
    assert r["A_pd_z1"] == pytest.approx(11.2574, rel=1e-4)
    assert r["A_pe_z1"] == pytest.approx(11.6796, rel=1e-4)
    assert r["A_dh_z1"] == pytest.approx(14.2381, rel=1e-4)


def test_the_lateral_area_follows_the_mass(tmp_path: Path) -> None:
    """Doubling the liver mass doubles every lateral area."""
    r = _barrier(tmp_path)
    big = tmp_path / "big"
    big.mkdir()
    doubled = _builder(n=4)
    doubled.settings = dataclasses.replace(doubled.settings, liver_mass=mass_g(2000.0))
    r2 = _barrier(big, builder=doubled)
    for pid in ("A_pd_z1", "A_pe_z1", "A_dh_z1"):
        assert r2[pid] == pytest.approx(2.0 * r[pid])


def test_a_solute_at_the_radius_of_the_fenestra_has_no_fenestral_path(
    tmp_path: Path,
) -> None:
    """A radius equal to `r_fen` closes the sieve exactly: no reaction, finite run."""
    builder = _builder(radius=geometry(Species.RAT).r_fen)
    assert not [r for r in builder.radial_exchange() if str(r.sid).startswith("FEN_")]
    r = _barrier(tmp_path, endothelial_ps=1e-6, builder=builder)
    assert np.all(np.isfinite(r.simulate(0, 10, 11)))


def test_all_paths_at_zero_flow_conserve_the_substance(tmp_path: Path) -> None:
    """Fenestral, endothelial and basolateral exchange only move the tracer."""
    r = _barrier(tmp_path, endothelial_ps=1e-6, basolateral_ps=1e-6, flow=0.0)

    def total() -> float:
        return sum(
            r[f"[{sid}]"] * r[f"V_{sid.removeprefix('tr_')}"]
            for sid in r.model.getFloatingSpeciesIds()
        )

    before = total()
    r.simulate(0, 60, 61)
    assert total() == pytest.approx(before, rel=1e-6)
    for layer in ("endothelium", "disse", "cell"):
        assert r[f"[tr_{layer}_z1]"] > 0.0, f"nothing reached the {layer}"


def test_without_erythrocytes_the_area_is_that_of_the_plasma(tmp_path: Path) -> None:
    """A haematocrit of 0 has no erythrocyte compartment and still exchanges."""
    r = _barrier(tmp_path, haematocrit=0.0)
    geo = geometry(Species.RAT)
    assert r["A_pd_z1"] == pytest.approx(2.0 * r["V_plasma_z1"] * 1e-3 / geo.r_sin)
    assert r["FEN_tr_z1"] > 0.0


def test_a_single_zone_has_the_exchange(tmp_path: Path) -> None:
    """One zone of the sinusoid has the barrier at the mid porosity."""
    r = _barrier(tmp_path, n=1)
    geo = geometry(Species.RAT)
    assert r["eps_fen_z1"] == pytest.approx((geo.porosity_pp + geo.porosity_pc) / 2)
    assert r["FEN_tr_z1"] > 0.0


def test_the_well_stirred_geometry_has_no_radial_exchange() -> None:
    """Its vascular compartment contains the Disse space already."""
    builder = _builder(geometry=Geometry.WELL_STIRRED)
    assert builder.radial_exchange() == []
    model = builder.build()
    assert not [p for p in model.parameters if p.sid.startswith(("P_fen", "A_pd"))]


def test_the_processes_are_registered() -> None:
    """A profile can key the three processes of the barrier."""
    processes = _builder().processes()
    assert {
        "porosity",
        "endothelial_permeability",
        "basolateral_permeability",
        "axial_diffusion",
    } <= processes


def _mean_residence_time(r: roadrunner.RoadRunner) -> float:
    """The mean time (min) a bolus in the first plasma takes to leave the chain."""
    n = sum(1 for sid in r.model.getFloatingSpeciesIds() if "_plasma_" in sid)
    for k in range(1, n + 1):
        r[f"[{ids.species('tr', Layer.PLASMA, ids.zone(k))}]"] = 0.0
    r[f"[{ids.species('tr', Layer.PLASMA, ids.zone(1))}]"] = 1.0
    volume = sum(
        r[f"V_{layer}_z{k}"] for k in range(1, n + 1) for layer in ("plasma", "disse")
    )
    horizon = 15.0 * volume / r["Q_plasma_li"]
    last = ids.species("tr", Layer.PLASMA, ids.zone(n))
    out = r.simulate(0, horizon, 60001, selections=["time", f"[{last}]"])
    t, c = out[:, 0], out[:, 1]
    return float(np.trapezoid(t * c, t) / np.trapezoid(c, t))


def test_a_small_solute_stays_longer_than_the_intravascular_marker(
    tmp_path: Path,
) -> None:
    """A solute which crosses the fenestrae distributes into the space of Disse.

    The marker of 1 um cannot cross a fenestra, so its mean residence time is
    the plasma volume over the flow. A small solute leaves the lumen and
    returns, so it stays longer. As the barrier opens (porosity 1, no
    erythrocytes) the mean residence time approaches the plasma plus the space
    of Disse over the flow, the volume of distribution of the fenestral path
    alone; an open endothelial path would add the endothelium.
    """
    marker_dir = tmp_path / "marker"
    marker_dir.mkdir()
    marker = _mean_residence_time(
        _barrier(
            marker_dir,
            builder=_builder(
                radius=1e-6, haematocrit=0.0, porosity_pp=1.0, porosity_pc=1.0
            ),
        )
    )
    small = _mean_residence_time(_barrier(tmp_path, porosity=1.0, haematocrit=0.0))
    r = _barrier(tmp_path, porosity=1.0, haematocrit=0.0)
    q = r["Q_plasma_li"]
    plasma = sum(r[f"V_plasma_z{k}"] for k in range(1, 5))
    disse = sum(r[f"V_disse_z{k}"] for k in range(1, 5))
    assert marker == pytest.approx(plasma / q, rel=1e-6)
    assert small > marker
    assert small == pytest.approx((plasma + disse) / q, rel=1e-6)


def test_the_fenestral_clearance_is_p_times_a_in_litres_per_minute(
    tmp_path: Path,
) -> None:
    """A permeability in m/s on an area in m2 is a clearance in l/min.

    By hand: P = 5.45e-9 m/s on A = 11.26 m2 is 6.137e-8 m3/s = 6.137e-5 l/s
    = 3.682e-3 l/min. The rate of the reaction at 1 mM in the plasma and 0 in
    the space of Disse is that clearance in mmole/min. The barrier is the rate
    limiting step here (no flow, nothing else is set), and a wrong factor
    (60000 for the units, say) would give 220 l/min or 6e-8 l/min.
    """
    builder = ZonationBuilder(
        ZonationSettings(
            n_zones=1, liver_mass=mass_g(1000.0), reference_flow=liver_flow()
        ),
        [TRACER],
        None,
    )
    path = tmp_path / "barrier.xml"
    create_model(
        builder.build(), filepath=path, sbml_level=3, sbml_version=2, validate=False
    )
    r = roadrunner.RoadRunner(str(path))
    r["P_fen_tr"] = 5.45e-9
    r[ids.porosity("z1")] = 1.0
    r[ids.area("pd", "z1")] = 11.26
    r[f"[{ids.species('tr', Layer.PLASMA, 'z1')}]"] = 1.0
    r[f"[{ids.species('tr', Layer.DISSE, 'z1')}]"] = 0.0
    assert r[ids.exchange("FEN", "tr", "z1", "")] == pytest.approx(3.682e-3, rel=1e-3)


# --- the filtrate of a pressure-driven sinusoid ------------------------------
#
# The leaky tube loses plasma through the fenestrae where its flow falls and
# regains it where the flow rises, so a zone receives more solute than it
# passes on. The filtrate carries the solute it took with it (the solvent drag
# and the chain of the space of Disse), and the test of that is the steady
# profile of a tracer at a constant inflow: flat for a solute the fenestrae do
# not hinder, enriched by `Q(z1)/Q(zk)` for one they stop.

_PRESSURE_PROFILE = (
    1.0,
    1.094652,
    1.175702,
    1.235407,
    1.267174,
    1.267174,
    1.235407,
    1.175702,
    1.094652,
    1.0,
)
"""`Q(z1)/Q(z<k>)` of a default pressure flow at ten zones, the enrichment of a
solute which stays in the lumen (the flow falls to the middle of the sinusoid
and rises back, `pkpdsim.zonation.flow`)."""


def _steady(
    tmp_path: Path,
    radius: float,
    *,
    n: int = 10,
    geometry: Geometry = Geometry.SINUSOID,
    p_onc: float | None = None,
) -> roadrunner.RoadRunner:
    """A pressure-driven chain at the steady state of an inflow of 1 mM.

    The axial diffusion is switched off (`D_tr = 0`), so the convection, the
    barrier and the filtrate alone set the profile and it can be compared with
    the flows of the zones.
    """
    settings = ZonationSettings(
        n_zones=n,
        geometry=geometry,
        flow=PressureFlow(p_onc=p_onc),
        bile=None,
        haematocrit=0.0,
        liver_mass=mass_g(7.229),
        reference_flow=liver_flow(),
    )
    builder = _Sized(settings, [TRACER], None)
    builder.radius = radius
    path = tmp_path / "filtrate.xml"
    create_model(
        builder.build(), filepath=path, sbml_level=3, sbml_version=2, validate=False
    )
    r = individual(path)
    r.integrator.absolute_tolerance = 1e-14
    r.integrator.relative_tolerance = 1e-12
    if ids.diffusion("tr") in r.model.getGlobalParameterIds():
        r[ids.diffusion("tr")] = 0.0
    layer = ids.VASCULAR if geometry is Geometry.WELL_STIRRED else Layer.PLASMA
    r[ids.species("tr", layer, ids.INFLOW)] = 1.0
    r.simulate(0.0, 200.0, 2)
    return r


def _profile(r: roadrunner.RoadRunner, layer: Layer | str, n: int = 10) -> list[float]:
    """The concentrations of a layer over the zones."""
    return [r[f"[{ids.species('tr', layer, ids.zone(k))}]"] for k in range(1, n + 1)]


def test_a_freely_sieved_solute_keeps_the_concentration_of_the_inflow(
    tmp_path: Path,
) -> None:
    """The filtered water takes the solute with it and brings it back.

    A solute which no fenestra hinders (`sigma = 0`) is dragged out of the
    lumen with the plasma which filters and back in where the lumen regains
    flow, so neither the lumen nor the space of Disse is enriched: the profile
    is the concentration of the inflow in every zone. Without the solvent drag
    and the chain of the filtrate it is `_PRESSURE_PROFILE`, a hump of 26 %
    out of nothing.
    """
    r = _steady(tmp_path, radius=1e-12)
    assert _profile(r, Layer.PLASMA) == pytest.approx([1.0] * 10, rel=1e-6)
    assert _profile(r, Layer.DISSE) == pytest.approx([1.0] * 10, rel=1e-6)


def test_a_solute_no_fenestra_lets_through_follows_the_flow_of_the_zones(
    tmp_path: Path,
) -> None:
    """A fully reflected solute stays in the lumen and is concentrated there.

    `sigma = 1`: the water leaves it behind, so the flux `Q*c` through every
    zone is the one of the inflow and the concentration is `Q(z1)/Q(zk)`. This
    is the physics of a solute which cannot cross the endothelium, and the
    reason the drag is scaled by the reflection coefficient and not applied to
    every solute.
    """
    r = _steady(tmp_path, radius=1e-6)
    assert _profile(r, Layer.PLASMA) == pytest.approx(_PRESSURE_PROFILE, rel=1e-5)
    assert _profile(r, Layer.DISSE) == pytest.approx([0.0] * 10, abs=1e-12)


def test_a_well_stirred_zone_carries_the_whole_axial_water_flux(
    tmp_path: Path,
) -> None:
    """One compartment holds the lumen and the space of Disse, so nothing filters out.

    The filtrate stays in the vascular compartment of the zone and flows on
    with it, so the chain carries the flow which enters it and the profile is
    flat; with the flow of the zone it would be `_PRESSURE_PROFILE`.
    """
    r = _steady(tmp_path, radius=1e-12, geometry=Geometry.WELL_STIRRED)
    assert _profile(r, ids.VASCULAR) == pytest.approx([1.0] * 10, rel=1e-9)


@pytest.mark.parametrize(
    "p_onc,lymph", [(None, False), (0.0, True), (2000.0, True)], ids=str
)
def test_only_a_net_filtration_leaves_at_the_pericentral_end(
    p_onc: float | None, lymph: bool, tmp_path: Path
) -> None:
    """The default oncotic pressure filters and reabsorbs the same amount.

    It is the mean of the two pressures, so `Q(z<n>)` is `Q(z1)` and the space
    of Disse gives everything back to the lumen: no reaction takes the
    filtrate out of the organ. A lower or a higher oncotic pressure filters or
    reabsorbs on balance, and that remainder leaves (or enters) at the
    pericentral end and is in the outflow ledger, which keeps the mass balance
    exact (`tests/zonation/test_observers.py`).
    """
    r = _steady(tmp_path, radius=1e-12, p_onc=p_onc)
    out = ids.flow(
        "tr", ids.site(Layer.DISSE, ids.zone(10)), ids.site(Layer.DISSE, ids.OUTFLOW)
    )
    assert (out in r.model.getReactionIds()) is lymph
    assert _profile(r, Layer.PLASMA) == pytest.approx([1.0] * 10, rel=1e-4)


def test_an_imposed_flow_filters_nothing(tmp_path: Path) -> None:
    """Every zone has the same flow, so no water crosses the endothelium.

    The filtrate is the leaky tube's and nothing else emits it: an imposed
    flow has neither the drag nor the chain of the space of Disse nor the
    parameters of the filtrate.
    """
    model = _builder(n=10, flow=1.0, radius=1e-12).build()
    reactions = {str(x.sid) for x in model.reactions}
    assert not [x for x in reactions if x.startswith(("SD_", "FLOW_tr_disse_"))]
    assert not [str(p.sid) for p in model.parameters if str(p.sid).startswith("Q_fil_")]
