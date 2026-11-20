"""Tests of what a lobed liver reports: the per lobe readouts of the body.

A lobed liver keeps every observer of the liver without lobes, because the
organ still drains one concentration into the venous pool: the lobes change
how that concentration comes about and not what the body reads. What the lobes
add are readouts of their own, which a study of a regional intervention needs:
the concentration leaving every lobe, its perfusion per volume, its volume and
share, and the remnant fraction of the liver.

The exposure of a lobe is the area under its outflow concentration, which the
library computes from the simulated course as it computes every area
(`pkpdsim.pbpk.regional`); the model holds no accumulating AUC of its own.
"""

from pathlib import Path

import pytest

from pkpdsim.pbpk import (
    FirstOrder,
    Metabolism,
    Organ,
    Part,
    PBPKDefinition,
    Route,
    Variants,
    ids,
)
from pkpdsim.pbpk.builder import build
from pkpdsim.pbpk.interventions import PortalLigation
from pkpdsim.pbpk.organ import OrganComponent
from pkpdsim.pbpk.variants.liver_lobes import LobedLiver
from pkpdsim.species import Species

from .conftest import definition, individual
from .test_lobes_build import RAT_TOKENS

S = "drg"


def _definition(liver: OrganComponent | None) -> PBPKDefinition:
    return definition(
        sid="m",
        species=Species.RAT,
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


def _built(tmp_path: Path, liver: OrganComponent | None) -> Path:
    return build(_definition(liver), tmp_path).flat_path


def _names(tmp_path: Path, liver: OrganComponent | None) -> set[str]:
    r = individual(_built(tmp_path, liver))
    return {str(i) for i in r.model.getGlobalParameterIds()}


def _of_the_body(names: set[str]) -> set[str]:
    """The names of the top model: the ids an organ brought carry its prefix."""
    return {name for name in names if "__" not in name}


def test_a_lobed_liver_keeps_every_readout_of_the_body(tmp_path: Path) -> None:
    """The body reads what it read before: the lobes add names and remove none.

    The elements inside the liver do change, as they do with the zones: a
    lobed liver has no single `Vti` or `Cu_ti_<s>` of the whole organ, it has
    one per lobe, exactly as a zonated liver has one per zone. What may not
    change is what leaves the organ, which is the top model.
    """
    without = _of_the_body(_names(tmp_path / "a", None))
    with_lobes = _of_the_body(_names(tmp_path / "b", LobedLiver()))
    assert without <= with_lobes
    assert ids.organ_outflow(Organ.LIVER, S) in with_lobes


def test_every_lobe_reports_the_concentration_leaving_it(tmp_path: Path) -> None:
    names = _names(tmp_path, LobedLiver())
    for token in RAT_TOKENS:
        assert ids.lobe_outflow(Organ.LIVER, token, S) in names


def test_the_outflow_of_the_liver_is_the_mix_of_the_lobes(tmp_path: Path) -> None:
    """The flow weighted mix of the lobe outflows is what the body receives."""
    r = individual(_built(tmp_path, LobedLiver()))
    r[ids.iv_dose(S)] = 10.0
    lobes = [ids.lobe_outflow(Organ.LIVER, t, S) for t in RAT_TOKENS]
    flows = [ids.lobe_flow(ids.QH, t) for t in RAT_TOKENS]
    whole = ids.organ_outflow(Organ.LIVER, S)
    result = r.simulate(0.0, 240.0, 61, selections=["time", whole, *lobes, *flows])
    n = len(RAT_TOKENS)
    for k in range(1, result.shape[0]):
        mix = sum(result[k, 1 + n + j] * result[k, 1 + j] for j in range(n)) / sum(
            result[k, 1 + n + j] for j in range(n)
        )
        assert mix == pytest.approx(result[k, 0 + 1], rel=1e-8)


def test_every_lobe_reports_its_perfusion_and_its_share(tmp_path: Path) -> None:
    names = _names(tmp_path, LobedLiver())
    for token in RAT_TOKENS:
        assert ids.lobe_perfusion(token) in names
        assert ids.lobe_volume(token) in names
        assert ids.lobe_volume_fraction(token) in names
    assert ids.F_REMNANT in names
    assert ids.VLI_LOBES in names
    assert ids.VLI_REF in names


def test_a_liver_without_lobes_reports_no_lobe(tmp_path: Path) -> None:
    names = _names(tmp_path, None)
    assert not [n for n in names if n.startswith("f_vol_")]
    assert ids.F_REMNANT not in names


def test_the_outflow_of_the_liver_is_the_mean_over_the_lobes_which_flow(
    tmp_path: Path,
) -> None:
    """`Cout_li_<s>` is weighted by the flows of the lobes, not by the flow of
    the organ.

    With the supply fractions the flow of a lobe is a share of the flow of the
    liver, so the two agree while a branch of each kind is open; they do not
    when every branch of one kind is closed, and the mix has to stay the mean
    of the lobes which flow. The portal branches of every lobe closed is the
    case: the liver is then perfused by its artery alone.
    """
    r = individual(
        _built(tmp_path, LobedLiver()),
        **PortalLigation(RAT_TOKENS).changes(RAT_TOKENS),
    )
    r[ids.iv_dose(S)] = 10.0
    out = ids.organ_outflow(Organ.LIVER, S)
    lobes = [ids.lobe_outflow(Organ.LIVER, t, S) for t in RAT_TOKENS]
    flows = [ids.lobe_flow(ids.QH, t) for t in RAT_TOKENS]
    r.simulate(0.0, 120.0, 13)
    total = sum(r[f] for f in flows)
    assert total > 0.0
    assert total < r[ids.QH], "the lobes carry the arterial inflow alone"
    weighted = sum(r[f] * r[c] for f, c in zip(flows, lobes, strict=True)) / total
    assert r[out] == pytest.approx(weighted, rel=1e-9)


def test_the_zonal_means_of_a_lobed_liver_keep_the_mean_of_the_whole_organ(
    tmp_path: Path,
) -> None:
    """The lobes add a mean per lobe and keep the mean of the whole liver.

    `Cpl_pp_<s>` and `Cpl_pc_<s>` of the zonated liver are the two means the
    observers of the zonation engine report and what a PD input or an
    observer of a model reads; the lobes add `Cpl_<lobe>_pp_<s>` beside them
    and leave the unprefixed pair in place, the mean of the lobes weighted by
    the volume each lobe has.
    """
    from pkpdsim.pbpk.variants.liver_zonated import ZonatedLiver

    liver = LobedLiver(zonation=ZonatedLiver(n_zones=3, means=True))
    r = individual(_built(tmp_path, liver))
    r[ids.iv_dose(S)] = 10.0
    r.simulate(0.0, 10.0, 3)
    for periportal in (True, False):
        whole = ids.flat(Organ.LIVER, ids.zonal_mean(Part.PLASMA, S, periportal))
        lobes = [
            (
                ids.lobe_volume_fraction(token),
                ids.flat(
                    Organ.LIVER, ids.zonal_mean(Part.PLASMA, S, periportal, token)
                ),
            )
            for token in RAT_TOKENS
        ]
        total = sum(r[f] * r[c] for f, c in lobes)
        assert r[whole] == pytest.approx(total, rel=1e-9)
        assert r[whole] > 0.0
