"""Tests of the lobed liver: the liver as lobes perfused in parallel.

A lobe is a parallel chain with its own share of the liver mass, its own
hepatic artery and portal vein inflow and its own drainage; the lobes of a
liver sum to the whole organ. With the shares of the supply proportional to
the shares of the mass every lobe has the same perfusion per volume and the
same extraction, so a lobed liver is then the lumped liver: each lobe
extracts `CL/(Q + CL)` of its inflow, whatever its share, and the mix of the
lobe outflows is the outflow of the whole organ. That equivalence is the
first test of the structure.
"""

from pathlib import Path

import numpy as np
import pytest

from pkpdsim.pbpk import (
    FirstOrder,
    Metabolism,
    Organ,
    PBPKDefinition,
    Route,
    U,
    Variants,
    ZonatedLiver,
    ids,
)
from pkpdsim.pbpk.builder import build
from pkpdsim.pbpk.lobes import Lobe, lobe_table
from pkpdsim.pbpk.organ import OrganComponent
from pkpdsim.pbpk.organs import Part, Pool
from pkpdsim.pbpk.variants.liver_lobes import MAX_SITES, LobedLiver
from pkpdsim.sbml.values import Value
from pkpdsim.species import Species

from .conftest import definition, individual

S = "drg"
RAT_TOKENS = ("cl", "lll", "lml", "rl", "rml")
DAY = 24 * 60.0
TIMES = np.linspace(0.0, DAY, 1441)
WHOLE = Lobe("whole", "whole liver", Value(1.0, U.dimensionless, "ModelChoice"))


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
                protein="CYP3A4",
            )
        ],
        routes=[Route.iv(S)],
        variants=Variants(organs={Organ.LIVER: liver})
        if liver is not None
        else Variants(),
    )


def _built(
    tmp_path: Path,
    sid: str,
    liver: OrganComponent | None,
    species: Species = Species.RAT,
) -> Path:
    result = build(_definition(sid, liver, species), tmp_path)
    return result.flat_path


def test_one_lobe_liver_is_the_standard_liver(tmp_path: Path) -> None:
    """A single lobe carries no lobe token, so the model is the standard liver."""
    lobed = _built(tmp_path / "a", "m", LobedLiver(lobes=[WHOLE]))
    standard = _built(tmp_path / "b", "m", None)
    assert _dateless(lobed) == _dateless(standard)


def _dateless(path: Path) -> list[str]:
    return [
        line
        for line in path.read_text().splitlines()
        if "Date" not in line and "dcterms" not in line and "W3CDTF" not in line
    ]


def test_rat_lobes_have_a_compartment_per_lobe(tmp_path: Path) -> None:
    r = individual(_built(tmp_path, "m", LobedLiver()))
    for token in RAT_TOKENS:
        assert (
            ids.flat(Organ.LIVER, ids.compartment(Part.PLASMA, token))
            in r.model.getCompartmentIds()
        )


def test_lobe_volumes_sum_to_the_liver_volume(tmp_path: Path) -> None:
    r = individual(_built(tmp_path, "m", LobedLiver()))
    total = sum(r[ids.lobe_volume(t)] for t in RAT_TOKENS)
    assert total == pytest.approx(r[ids.organ_volume(Organ.LIVER)], rel=1e-9)


def test_lobe_volume_shares_are_those_of_the_table(tmp_path: Path) -> None:
    r = individual(_built(tmp_path, "m", LobedLiver()))
    table = {lobe.token: lobe.mass_fraction.value for lobe in lobe_table(Species.RAT)}
    for token, share in table.items():
        assert r[ids.lobe_volume_fraction(token)] == pytest.approx(share, rel=1e-9)


@pytest.mark.parametrize("flow", ids.HEPATIC_FLOWS)
def test_lobe_flows_sum_to_the_liver_flows(tmp_path: Path, flow: str) -> None:
    r = individual(_built(tmp_path, "m", LobedLiver()))
    total = sum(r[ids.lobe_flow(flow, t)] for t in RAT_TOKENS)
    assert total == pytest.approx(r[flow], rel=1e-9)


def test_perfusion_per_volume_is_uniform_at_baseline(tmp_path: Path) -> None:
    """Every lobe is supplied in proportion to its mass, so the perfusion per volume is one number."""
    r = individual(_built(tmp_path, "m", LobedLiver()))
    q = [r[ids.lobe_perfusion(t)] for t in RAT_TOKENS]
    assert max(q) == pytest.approx(min(q), rel=1e-9)
    assert min(q) > 0.0


def test_the_conductances_are_open_at_baseline(tmp_path: Path) -> None:
    r = individual(_built(tmp_path, "m", LobedLiver()))
    for token in RAT_TOKENS:
        assert r[ids.arterial_conductance(token)] == 1.0
        assert r[ids.portal_conductance(token)] == 1.0
        assert r[ids.lobe_present(token)] == 1.0


def test_the_remnant_is_the_whole_liver_at_baseline(tmp_path: Path) -> None:
    r = individual(_built(tmp_path, "m", LobedLiver()))
    assert r[ids.F_REMNANT] == pytest.approx(1.0, rel=1e-9)


def test_proportional_lobes_match_the_lumped_liver(tmp_path: Path) -> None:
    """Parallel lobes with the supply proportional to the mass are the lumped liver.

    The structural invariant of the lobes: a lobe of the share `w` holds `w`
    of the volume, receives `w` of the flow and holds `w` of the clearance,
    so it extracts `CL/(Q + CL)` of what reaches it, as every other lobe and
    as the whole liver does, and the mix of the lobes is the outflow of the
    lumped liver. The two models are different ODE systems of different size,
    integrated separately, so the agreement is of the solver and not exact;
    `1e-4` is far below any structural mismatch (a wrong weight of a lobe or
    a flow which does not follow its volume moves the curve by percent).
    """
    lobed = individual(_built(tmp_path / "a", "m", LobedLiver()))
    lumped = individual(_built(tmp_path / "b", "m", None))
    curves = []
    for r in (lobed, lumped):
        r[ids.iv_dose(S)] = 10.0
        cve = f"[{ids.pool_species(Pool.VENOUS, S)}]"
        curves.append(r.simulate(times=TIMES, selections=["time", cve])[:, 1])
    a, b = curves
    assert np.allclose(a[5:], b[5:], rtol=1e-4, atol=1e-12)


def test_lobes_on_the_pig_raise(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="pig"):
        _built(tmp_path, "m", LobedLiver(), species=Species.PIG)


def test_too_many_sites_raise() -> None:
    with pytest.raises(ValueError, match=str(MAX_SITES)):
        LobedLiver(lobes=lobe_table(Species.HUMAN), zonation=ZonatedLiver(n_zones=10))


def test_a_lobe_set_which_is_no_set_is_refused() -> None:
    with pytest.raises(ValueError, match="sum to"):
        LobedLiver(lobes=[Lobe("a", "a", Value(0.5, U.dimensionless, "ModelChoice"))])


def test_the_lobe_objects_of_a_variant_are_public() -> None:
    """A user configures a lobed liver from `pkpdsim.pbpk`, as with every variant.

    `ZonatedLiver` and `PressureFlow` are names of the package, so the lobed
    liver and the objects which configure it are too: a study imports the
    variant, its buffer response and its drivers from one place.
    """
    import pkpdsim.pbpk as pbpk

    names = ("LobedLiver", "HABR", "FlowDriven", "PrescribedVolume", "PrescribedFlow")
    for name in names:
        assert name in pbpk.__all__, name
        assert getattr(pbpk, name) is not None, name
