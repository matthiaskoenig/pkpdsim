"""A solute is the minimal protocol of what rides the chain: an id and a name.

The perfusate quantities of the machine perfusion (oxygen, glucose, lactate, the
enzymes) are no PBPK substances, so they have no molar mass; the engine carries
them through the zones, and a molar mass or a diffusion coefficient adds the
transport which needs it.
"""

from dataclasses import dataclass

import pytest

from pkpdsim.pbpk.substance import Substance
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Q, Value
from pkpdsim.zonation import ZonationBuilder, ids
from pkpdsim.zonation.flow import ImposedFlow
from pkpdsim.zonation.geometry import Layer
from pkpdsim.zonation.settings import Geometry, ZonationSettings
from pkpdsim.zonation.transport import M2_PER_S, Solute

from .conftest import mass_g


@dataclass(frozen=True)
class _Bare:
    """A solute with an id and a name and nothing else."""

    sid: str = "o2"
    name: str = "oxygen"
    mr: Q | None = None
    diffusion: Value | None = None


def _builder(solute: Solute, geometry: Geometry = Geometry.SINUSOID) -> ZonationBuilder:
    settings = ZonationSettings(
        n_zones=3,
        geometry=geometry,
        flow=ImposedFlow(0.01),
        liver_mass=mass_g(1000.0),
    )
    return ZonationBuilder(settings, [solute], None)


def _reactions(solute: Solute, geometry: Geometry = Geometry.SINUSOID) -> set[str]:
    model = _builder(solute, geometry).build()
    return {str(r.sid) for r in model.reactions}


def _parameters(solute: Solute) -> set[str]:
    return {str(p.sid) for p in _builder(solute).build().parameters}


def test_a_substance_is_a_solute() -> None:
    assert isinstance(Substance("s", "s", Value(100.0, U.g_per_mole)), Solute)
    assert isinstance(_Bare(), Solute)


def test_a_bare_solute_has_a_species_in_every_layer_of_every_zone() -> None:
    builder = _builder(_Bare())
    sids = {str(s.sid) for s in builder.build().species}
    expected = {
        ids.species("o2", layer, ids.zone(k))
        for layer in builder.layers()
        for k in (1, 2, 3)
    }
    duct = {ids.species("o2", Layer.BILE, ids.DUCT)}
    assert sids == expected | duct


def test_a_bare_solute_rides_the_chain_and_has_no_diffusion() -> None:
    reactions = _reactions(_Bare())
    assert "FLOW_o2_plasma_z1_plasma_z2" in reactions
    assert not [r for r in reactions if r.startswith(("D_", "FEN_"))]
    assert not [p for p in _parameters(_Bare()) if p.startswith(("D_o2", "P_fen_"))]


def test_a_bare_solute_reaches_the_cells_through_the_membrane_paths() -> None:
    """The molar mass gates the fenestral path and no other one.

    The radius of the molecule, which the molar mass gives, is what the pore
    needs; the transcellular and the basolateral permeabilities are model
    input in m/s, so a perfusate quantity without a molar mass has them and
    can be let into the cells from outside.
    """
    assert {"END_o2_z1_in", "END_o2_z1_out", "PS_o2_z1"} <= _reactions(_Bare())
    assert {"P_endo_in_o2", "P_endo_out_o2", "P_hep_o2"} <= _parameters(_Bare())


def test_a_given_diffusion_without_a_molar_mass_diffuses_but_has_no_pore_path() -> None:
    solute = _Bare(diffusion=Value(1e-9, M2_PER_S))
    assert "D_o2_plasma_z1_plasma_z2" in _reactions(solute)
    assert "D_o2" in _parameters(solute)
    assert not [r for r in _reactions(solute) if r.startswith("FEN_")]


def test_a_molar_mass_derives_the_diffusion_and_the_pore_path() -> None:
    solute = _Bare(mr=Value(180.0, U.g_per_mole))
    assert "D_o2_plasma_z1_plasma_z2" in _reactions(solute)
    assert "FEN_o2_z1" in _reactions(solute)


def test_a_solute_without_a_molar_mass_rides_the_well_stirred_chain() -> None:
    reactions = _reactions(_Bare(), Geometry.WELL_STIRRED)
    assert "FLOW_o2_vas_z1_vas_z2" in reactions


def test_a_non_solute_is_refused() -> None:
    class _NoName:
        sid = "x"

    assert not isinstance(_NoName(), Solute)
    with pytest.raises(TypeError, match="solute"):
        ZonationBuilder(
            ZonationSettings(liver_mass=mass_g(1000.0)),
            [_NoName()],  # ty: ignore[invalid-argument-type]
            None,
        )
