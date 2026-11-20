"""Tests of the two spatial axes of the liver: the lobes of zones.

A lobe is a parallel chain, a zone a step along the sinusoid of a chain, and
the two compose: a site (`ids.site`) carries the lobe and the zone, the zones
of a lobe are perfused in series and no reaction crosses from one lobe into
another. A liver of one lobe of zones is the zonated liver, line by line: the
lobe axis is the axis of length one there, as the zone axis is in a lumped
lobed liver.

The invariant of the capacity: a site holds the share of the capacity of its
lobe (the volume share) times the share of its zone (the placement weight of
the profile, which sums to 1 over the zones), so the sites of the liver sum to
the whole capacity whatever the two axes are.
"""

from pathlib import Path

import pytest

from pkpdsim.pbpk import (
    FirstOrder,
    Metabolism,
    Organ,
    PBPKDefinition,
    Route,
    Variants,
    ZonatedLiver,
    ids,
)
from pkpdsim.pbpk.builder import build
from pkpdsim.pbpk.organ import OrganComponent
from pkpdsim.pbpk.organs import Part, Pool
from pkpdsim.pbpk.variants.liver_lobes import LobedLiver
from pkpdsim.species import Species
from pkpdsim.zonation.profiles import Linear
from pkpdsim.zonation.settings import Geometry

from .conftest import definition, individual
from .test_lobes_build import RAT_TOKENS, WHOLE, _dateless

S = "drg"
PROTEIN = "CYP3A4"


def _definition(
    sid: str, liver: OrganComponent | None, species: Species = Species.RAT
) -> PBPKDefinition:
    return definition(
        sid=sid,
        species=species,
        chemistry=[
            Metabolism(
                sid="DRG2M",
                organ=Organ.LIVER,
                substrate=S,
                products={"drg_m": 1.0},
                kinetics=FirstOrder(k=0.5),
                protein=PROTEIN,
            )
        ],
        routes=[Route.iv(S)],
        variants=Variants(organs={Organ.LIVER: liver})
        if liver is not None
        else Variants(),
    )


def _built(
    tmp_path: Path,
    liver: OrganComponent | None,
    species: Species = Species.RAT,
) -> Path:
    return build(_definition("m", liver, species), tmp_path).flat_path


def _reactions(path: Path) -> list[str]:
    # the instance is held while its ids are read: the model of a roadrunner
    # which is freed first is a dangling handle
    r = individual(path)
    return [str(i) for i in r.model.getReactionIds()]


def test_one_lobe_of_zones_is_the_zonated_liver(tmp_path: Path) -> None:
    """A single lobe carries no lobe token, so the model is the zonated liver."""
    zonation = ZonatedLiver(n_zones=4)
    lobed = _built(tmp_path / "a", LobedLiver(lobes=[WHOLE], zonation=zonation))
    zoned = _built(tmp_path / "b", zonation)
    assert _dateless(lobed) == _dateless(zoned)


def test_the_sites_are_the_product_of_the_axes(tmp_path: Path) -> None:
    r = individual(_built(tmp_path, LobedLiver(zonation=ZonatedLiver(n_zones=3))))
    compartments = set(r.model.getCompartmentIds())
    for token in RAT_TOKENS:
        for k in (1, 2, 3):
            site = ids.site(token, ids.zone(k))
            assert site is not None
            assert (
                ids.flat(Organ.LIVER, ids.compartment(Part.PLASMA, site))
                in compartments
            )


def test_each_lobe_has_its_own_inflow_and_outflow(tmp_path: Path) -> None:
    """The lobes are parallel: each is fed at its periportal end and drains at its pericentral one.

    There is one inflow and one outflow per lobe, whatever the zones of its
    chain, so both are keyed by the lobe and not by the site they reach.
    """
    reactions = _reactions(
        _built(tmp_path, LobedLiver(zonation=ZonatedLiver(n_zones=3)))
    )
    for token in RAT_TOKENS:
        for pool in (Pool.ARTERIAL, Pool.PORTAL):
            inflow = ids.inflow_reaction(pool, Part.PLASMA, S, token)
            assert ids.flat(Organ.LIVER, inflow) in reactions
        outflow = ids.outflow_reaction(Part.PLASMA, S, token)
        assert ids.flat(Organ.LIVER, outflow) in reactions


def test_the_sinusoid_flow_stays_inside_a_lobe(tmp_path: Path) -> None:
    """No reaction carries blood from a zone of one lobe into a zone of another."""
    reactions = _reactions(
        _built(tmp_path, LobedLiver(zonation=ZonatedLiver(n_zones=3)))
    )
    sinusoid = [r for r in reactions if "J_sin_" in r]
    assert sinusoid
    for r in sinusoid:
        assert len([t for t in RAT_TOKENS if f"{t}_z" in r]) == 1, r
    for token in RAT_TOKENS:
        # three zones in a chain are two steps, of every carried substance
        steps = [r for r in sinusoid if f"{token}_z" in r]
        assert len(steps) >= 2
        assert len(steps) % 2 == 0


def test_a_site_volume_is_the_lobe_share_of_the_zone_share(tmp_path: Path) -> None:
    """The zones of a lobe share its volume, and the lobes share the volume of the liver."""
    n = 4
    r = individual(_built(tmp_path, LobedLiver(zonation=ZonatedLiver(n_zones=n))))

    def volume(token: str, k: int) -> float:
        site = ids.site(token, ids.zone(k))
        assert site is not None
        return float(r[ids.flat(Organ.LIVER, ids.compartment(Part.PLASMA, site))])

    total = sum(volume(t, k) for t in RAT_TOKENS for k in range(1, n + 1))
    for token in RAT_TOKENS:
        of_lobe = [volume(token, k) for k in range(1, n + 1)]
        assert max(of_lobe) == pytest.approx(min(of_lobe), rel=1e-12)
        assert sum(of_lobe) / total == pytest.approx(
            r[ids.lobe_volume_fraction(token)], rel=1e-9
        )


def test_the_capacity_of_the_sites_sums_to_the_capacity_of_the_liver(
    tmp_path: Path,
) -> None:
    """The share of a site is its lobe share times its zone weight, so the sites sum to 1."""
    liver = LobedLiver(
        zonation=ZonatedLiver(n_zones=4, profiles={PROTEIN: Linear(1.0, 4.0)})
    )
    r = individual(_built(tmp_path, liver))
    weights = [
        r[ids.flat(Organ.LIVER, ids.placement_weight(PROTEIN.lower(), ids.zone(k)))]
        for k in range(1, 5)
    ]
    assert sum(weights) == pytest.approx(1.0, rel=1e-12)
    shares = [r[ids.lobe_volume_fraction(t)] for t in RAT_TOKENS]
    assert sum(shares) == pytest.approx(1.0, rel=1e-12)


def test_the_zonal_means_are_per_lobe(tmp_path: Path) -> None:
    """Each lobe reports the periportal and the pericentral mean of its own chain."""
    liver = LobedLiver(zonation=ZonatedLiver(n_zones=4, means=True))
    r = individual(_built(tmp_path, liver))
    for token in RAT_TOKENS:
        for periportal in (True, False):
            mean = ids.zonal_mean(Part.PLASMA, S, periportal, token)
            assert r[ids.flat(Organ.LIVER, mean)] == pytest.approx(0.0)


def test_the_axial_diffusion_stays_inside_a_lobe(tmp_path: Path) -> None:
    """The sinusoid diffuses along a lobe and never between two lobes."""
    liver = LobedLiver(
        zonation=ZonatedLiver(n_zones=3, geometry=Geometry.SINUSOID),
        lobes=None,
    )
    reactions = _reactions(_built(tmp_path, liver))
    axial = [r for r in reactions if "__D_" in r]
    assert axial
    for r in axial:
        assert len([t for t in RAT_TOKENS if f"{t}_z" in r]) == 1, r
    for token in RAT_TOKENS:
        # three zones in a chain are two steps, of every diffusing part
        steps = [r for r in axial if f"{token}_z" in r]
        assert len(steps) >= 2
        assert len(steps) % 2 == 0
