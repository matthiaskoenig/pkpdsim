"""The fenestrated barrier and the radial exchange of the zones.

The tests use a tracer, a substance which does not react, so no cell model is
needed. The tracer starts in the plasma of every zone and the rates of the
radial reactions are read at the start.
"""

import dataclasses
from pathlib import Path

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
from pkpdsim.zonation.geometry import Layer, geometry
from pkpdsim.zonation.profiles import Weights
from pkpdsim.zonation.settings import Geometry, ZonationSettings
from pkpdsim.zonation.transport import (
    fenestral_permeability,
    molecular_radius,
)

from .conftest import mass_g

pytestmark = pytest.mark.filterwarnings("ignore:.*numerical dispersion:UserWarning")

TRACER = Substance("tr", "tracer", Value(180.0, U.g_per_mole))
"""A small solute, much smaller than a fenestra."""


class _Sized(ZonationBuilder):
    """A builder whose tracer has a given molecular radius (m), 0 for derived."""

    radius = 0.0

    def _radius(self, substance: Substance) -> Value:
        """Get the radius of the tracer."""
        if not self.radius:
            return super()._radius(substance)
        return Value(self.radius, U.meter)


def _builder(
    n: int = 4,
    radius: float = 0.0,
    porosity_pp: float = 0.06,
    porosity_pc: float = 0.08,
    **kwargs: object,
) -> _Sized:
    """A builder of a tracer chain with the given porosity gradient."""
    geo = dataclasses.replace(
        geometry(Species.RAT), porosity_pp=porosity_pp, porosity_pc=porosity_pc
    )
    settings = ZonationSettings(
        n_zones=n,
        flow=ImposedFlow(1.0),
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


def test_the_two_paths_are_parallel_and_additive(tmp_path: Path) -> None:
    """The cell path adds to the fenestral one."""
    only_fen = _barrier(tmp_path, endothelial_ps=0.0, porosity=0.06)
    both = _barrier(tmp_path, endothelial_ps=1e-6, porosity=0.06)
    assert _flux(both) > _flux(only_fen)
    assert both["FEN_tr_z1"] == pytest.approx(only_fen["FEN_tr_z1"])
    assert _flux(both) == pytest.approx(both["FEN_tr_z1"] + both["END_tr_z1_in"])


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


def test_an_open_barrier_equilibrates_plasma_and_disse(tmp_path: Path) -> None:
    """With open paths the two layers reach one concentration."""
    r = _barrier(tmp_path, endothelial_ps=1e3, porosity=1.0)
    r.simulate(0, 60, 61)
    assert r["[tr_plasma_z1]"] == pytest.approx(r["[tr_disse_z1]"], rel=1e-3)


def test_the_basolateral_path_moves_the_solute_into_the_hepatocyte(
    tmp_path: Path,
) -> None:
    """The hepatocyte fills when the basolateral path is set, and not otherwise."""
    r = _barrier(tmp_path, basolateral_ps=1e-5, porosity=1.0)
    r.simulate(0, 60, 61)
    assert r["[tr_hepatocyte_z1]"] > 0.0
    closed = _barrier(tmp_path, porosity=1.0)
    closed.simulate(0, 60, 61)
    assert closed["[tr_hepatocyte_z1]"] == 0.0


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


def test_the_fenestral_permeability_is_the_derived_one(tmp_path: Path) -> None:
    """The open pore permeability of the model is that of `transport`."""
    r = _barrier(tmp_path)
    geo = geometry(Species.RAT)
    d = molecular_radius(TRACER.mr)
    assert r["P_fen_tr"] > 0.0
    from pkpdsim.zonation.transport import diffusion_of

    expected = fenestral_permeability(diffusion_of(TRACER), d, geo, 1.0)
    assert r["P_fen_tr"] == pytest.approx(expected.value)


def test_the_lateral_area_follows_the_mass(tmp_path: Path) -> None:
    """The area of the interface is 2 V / r_sin of the lumen of the zone."""
    r = _barrier(tmp_path)
    geo = geometry(Species.RAT)
    lumen = r["V_plasma_z1"] + r["V_rbc_z1"]
    assert r["A_pd_z1"] == pytest.approx(2.0 * lumen * 1e-3 / geo.r_sin)
    outer = geo.r_sin + geo.t_endo + geo.t_disse
    assert r["A_dh_z1"] == pytest.approx(r["A_pd_z1"] * outer / geo.r_sin)
    big = tmp_path / "big"
    big.mkdir()
    doubled = _builder(n=4)
    doubled.settings = dataclasses.replace(doubled.settings, liver_mass=mass_g(2000.0))
    r2 = _barrier(big, builder=doubled)
    assert r2["A_pd_z1"] == pytest.approx(2.0 * r["A_pd_z1"])


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
