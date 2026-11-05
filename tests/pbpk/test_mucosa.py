"""Tests of the mucosa of the intestine: flow-limited and permeability-limited.

The mucosa of a substance with chemistry of the small intestine is
flow-limited by default: the mucosal blood flow carries the arterial blood
into the enterocytes, so a gut extraction of the absorbed drug extracts
systemic drug with about `Qmu*E_G` as well. `mucosa_ps` of the intestine
makes the mucosa of a substance permeability-limited (the Qgut model, Yang
2007): the enterocytes see the hybrid flow `Qgut` of the mucosal blood flow
and the permeability clearance of the unbound drug instead of the mucosal
blood flow, and the rest of the mucosal blood passes them, so the gut
availability of the absorbed drug and the systemic clearance of the gut are
separate.
"""

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING

import pytest

from pkpdsim.pbpk import (
    Binding,
    Erythrocytes,
    FirstOrder,
    Metabolism,
    Organ,
    PBPKDefinition,
    Pool,
    Route,
    SegmentedIntestine,
    SingleLumen,
    Substance,
    Variants,
    build,
    ids,
)

from .conftest import definition, load

if TYPE_CHECKING:
    import roadrunner

S, M, OTHER = "drg", "drg_m", "oth"
"""The drug, its gut metabolite and a substance without chemistry of the gut."""

FU, BP = 0.08, 142.0
"""The unbound fraction and the blood to plasma ratio of the drug (sirolimus)."""

SMALL = ("du", "jej1", "jej2", "ile1", "ile2")
"""The segments of the small intestine."""

COLON = ("cae", "asc", "tra", "des", "sig", "rec")
"""The segments of the large intestine."""

CLPERM = 100.0 / 60.0
"""The permeability-surface product of the whole mucosa (l/min)."""

CLINT_PERMEABILITY = 3.8815
"""The unbound intrinsic clearance of the gut (l/min) which gives the
permeability-limited mucosa the gut availability 0.3."""

CLINT_FLOW = 2036.15
"""The unbound intrinsic clearance of the gut (l/min) which gives the
flow-limited mucosa the gut availability 0.3."""


def _gu(local: str) -> str:
    """The flat id of a local id of the intestine."""
    return ids.flat(Organ.SMALL_INTESTINE, local)


def _drug() -> Substance:
    """The drug: strongly bound, in the erythrocytes, no tissue, no excretion."""
    return Substance(
        sid=S,
        name="drug",
        mr=900.0,
        binding=Binding(fu=FU),
        erythrocytes=Erythrocytes("equilibrium", bp=BP),
    )


def _definition(
    sid: str,
    k: float,
    *,
    mucosa_ps: "Mapping[str, float] | None" = None,
    per_volume: bool = False,
    single: bool = True,
) -> PBPKDefinition:
    """The drug with a gut metabolism, optionally permeability-limited.

    Args:
        sid: the id of the definition.
        k: the unbound intrinsic clearance of the whole small intestine (l/min).
        mucosa_ps: the exchange of the mucosa by substance, `None` flow-limited.
        per_volume: the exchange is per litre of mucosa (1/min).
        single: the single lumen (one mucosa segment), else the segments.

    Returns:
        The definition.
    """
    intestine: SingleLumen | SegmentedIntestine
    if mucosa_ps is None:
        intestine = SingleLumen() if single else SegmentedIntestine()
    elif single:
        intestine = SingleLumen(mucosa_ps=mucosa_ps, mucosa_ps_per_volume=per_volume)
    else:
        intestine = SegmentedIntestine(
            mucosa_ps=mucosa_ps, mucosa_ps_per_volume=per_volume
        )
    return definition(
        sid=sid,
        substances=[
            _drug(),
            Substance(sid=M, name="metabolite", mr=900.0),
            Substance(sid=OTHER, name="other drug", mr=300.0),
        ],
        chemistry=[
            Metabolism(
                sid="DRGGU",
                organ=Organ.SMALL_INTESTINE,
                substrate=S,
                products={M: 1.0},
                kinetics=FirstOrder(k=k, per_volume=False),
                protein="CYP3A4",
            )
        ],
        variants=Variants(organs={Organ.SMALL_INTESTINE: intestine}),
        routes=[Route.iv(S), Route.oral_solution(S), Route.iv(OTHER)],
    )


def _processes(segments: Sequence[str | None]) -> list[str]:
    """The flat ids of the gut metabolism in the mucosa segments."""
    return [_gu(ids.segment_process("DRGGU", seg)) for seg in segments]


def _qgut(
    r: "roadrunner.RoadRunner", clperm: float, segments: Sequence[str | None]
) -> float:
    """The hybrid flow of the Qgut model from the model's physiology.

    `Qgut = Qvilli*CLperm/(Qvilli + CLperm)` (Yang 2007) of the mucosal blood
    flow on the plasma concentration (`Qmu*BP`) and the permeability clearance
    of the unbound drug (`fu*PS`), summed over the mucosa segments (the
    permeability of a segment by its share of the mucosa volume).
    """
    total = sum(r[_gu(ids.mucosa_volume(seg))] for seg in segments)
    out = 0.0
    for seg in segments:
        share = r[_gu(ids.mucosa_volume(seg))] / total
        q = r[_gu(ids.mucosa_flow(seg))] * r[ids.bp(S)]
        ps = r[ids.fu(S)] * clperm * share
        out += q * ps / (q + ps)
    return out


def _availability(
    r: "roadrunner.RoadRunner", segments: Sequence[str | None], c: float = 1e-3
) -> float:
    """The gut availability of the absorbed drug, from the fluxes of the mucosa.

    Of what leaves the enterocytes, the share which reaches the blood: the
    model is loaded (every drug species 0, so the mucosa gets nothing from the
    blood) and the mucosa of every segment set to `c`, so the ratio of the
    outflow into the portal vein and the metabolism is `F_G` of the linear
    gut, flow-limited and permeability-limited alike.
    """
    for seg in segments:
        r[f"[{_gu(ids.mucosa(S, seg))}]"] = c
    to_blood = sum(r[_gu(ids.mucosa_reaction("out", S, seg))] for seg in segments)
    metabolism = sum(r[p] for p in _processes(segments))
    return to_blood / (to_blood + metabolism)


def _systemic_blood_clearance(
    r: "roadrunner.RoadRunner", segments: Sequence[str | None], dose: float = 1.0
) -> float:
    """The clearance of the gut on the arterial blood of an iv dose (l/hr).

    The sum of the gut metabolism over the arterial whole-blood concentration
    an hour after an iv dose, when the mucosa follows the blood (its time
    constant is under a minute).
    """
    r[ids.iv_dose(S)] = dose
    r.simulate(0, 60.0, 20)
    car = r[f"[{ids.pool_species(Pool.ARTERIAL, S)}]"] * r[ids.bp(S)]
    return sum(r[p] for p in _processes(segments)) / car * 60.0


def test_the_mucosa_without_the_option_is_flow_limited(tmp_path: Path) -> None:
    """Without `mucosa_ps` the mucosa keeps its flow and has no element of the option.

    The default creates no exchange (`kps_mu_<s>`, `PS_mu_<s>`, `Qgut_<s>`)
    and no bypass of the mucosal blood (`J_mu_blood_<s>`), and keeps the
    arterial inflow and the outflow into the portal vein of the mucosa.
    """
    r = load(build(_definition("flow", CLINT_FLOW), tmp_path))
    created = {
        ids.mucosa_kps(S),
        ids.mucosa_ps_ref(S),
        ids.mucosa_ps(S),
        ids.mucosa_qgut(S),
        ids.mucosa_reaction("blood", S),
    }
    present = set(r.model.getGlobalParameterIds()) | set(r.model.getReactionIds())
    assert not {_gu(i) for i in created} & present
    assert {_gu(ids.mucosa_reaction(k, S)) for k in ("in", "out")} <= present


def test_a_large_permeability_is_the_flow_limited_mucosa(tmp_path: Path) -> None:
    """The limit of a large exchange is the flow-limited mucosa.

    With `PS` far above the mucosal blood flow the hybrid flow `Qgut` is the
    flow `Qmu*BP` and the permeability-limited mucosa has the gut availability
    and the systemic clearance of the flow-limited one (both within 1e-4,
    the relative error of the limit at `fu*PS/(Qmu*BP)` about 1e5).
    """
    flow = load(build(_definition("limit_flow", CLINT_FLOW), tmp_path))
    limit = load(
        build(_definition("limit_ps", CLINT_FLOW, mucosa_ps={S: 1e8}), tmp_path)
    )
    segments: tuple[str | None, ...] = (None,)
    fg_flow = _availability(flow, segments)
    fg_limit = _availability(limit, segments)
    assert fg_limit == pytest.approx(fg_flow, rel=1e-4)
    cl_flow = _systemic_blood_clearance(flow, segments)
    cl_limit = _systemic_blood_clearance(limit, segments)
    assert cl_limit == pytest.approx(cl_flow, rel=1e-4)


def test_the_gut_availability_is_the_qgut_formula(tmp_path: Path) -> None:
    """The gut availability of the absorbed drug is `Qgut/(Qgut + fu*CLint)`.

    The linear case of the Qgut model (Yang 2007): the hybrid flow of the
    mucosal blood flow and the permeability clearance, computed from the
    physiology of the model and the given `CLperm`, with the unbound intrinsic
    clearance of the gut metabolism.
    """
    r = load(
        build(_definition("qgut", CLINT_PERMEABILITY, mucosa_ps={S: CLPERM}), tmp_path)
    )
    qgut = _qgut(r, CLPERM, (None,))
    expected = qgut / (qgut + r[ids.fu(S)] * CLINT_PERMEABILITY)
    assert _availability(r, (None,)) == pytest.approx(expected, rel=1e-9)
    assert r[_gu(ids.mucosa_qgut(S))] == pytest.approx(qgut, rel=1e-9)


def test_the_permeability_limited_mucosa_frees_the_systemic_clearance(
    tmp_path: Path,
) -> None:
    """The smallest failing case of the rapamycin survey.

    A drug with fu 0.08 and a blood to plasma ratio of 142 whose gut
    metabolism leaves the gut availability 0.3: with the flow-limited mucosa
    the gut clears systemic drug with about 20 l/hr of blood, far above the
    total blood clearance of sirolimus (1.5 to 4 l/hr); with the
    permeability-limited mucosa (`CLperm` 100 l/hr of the whole mucosa) the
    same gut availability costs 0.039 l/hr.
    """
    segments: tuple[str | None, ...] = (None,)
    flow = load(build(_definition("gut_flow", CLINT_FLOW), tmp_path))
    assert _availability(flow, segments) == pytest.approx(0.3, abs=0.01)
    cl_flow = _systemic_blood_clearance(flow, segments)
    assert cl_flow > 10.0
    ps = load(
        build(
            _definition("gut_ps", CLINT_PERMEABILITY, mucosa_ps={S: CLPERM}),
            tmp_path,
        )
    )
    assert _availability(ps, segments) == pytest.approx(0.3, abs=0.01)
    cl_ps = _systemic_blood_clearance(ps, segments)
    assert cl_ps < 0.1
    assert cl_ps == pytest.approx(0.0394, abs=0.002)
    assert cl_flow == pytest.approx(20.6, abs=0.5)


def test_the_segments_share_the_permeability_of_the_whole_mucosa(
    tmp_path: Path,
) -> None:
    """A `CLperm` of the whole mucosa is split over the segments by mucosa volume.

    The segmented intestine has the exchange and the hybrid flow per segment;
    the exchanges sum to the given value, and each hybrid flow is below the
    blood flow of its segment.
    """
    r = load(
        build(
            _definition(
                "segments", CLINT_PERMEABILITY, mucosa_ps={S: CLPERM}, single=False
            ),
            tmp_path,
        )
    )
    segments = SMALL + COLON
    total = sum(r[_gu(ids.mucosa_ps(S, seg))] for seg in segments)
    assert total == pytest.approx(CLPERM, rel=1e-12)
    for seg in segments:
        qgut = r[_gu(ids.mucosa_qgut(S, seg))]
        assert 0.0 < qgut < r[_gu(ids.mucosa_flow(seg))] * r[ids.bp(S)]


def test_an_exchange_per_litre_of_mucosa_follows_the_individual(
    tmp_path: Path,
) -> None:
    """`mucosa_ps_per_volume` is the exchange per litre of mucosa, `kps_mu_<s>*Vmu`."""
    result = build(
        _definition(
            "per_volume", CLINT_PERMEABILITY, mucosa_ps={S: 20.0}, per_volume=True
        ),
        tmp_path,
    )
    r = load(result)
    assert r[_gu(ids.mucosa_kps(S))] == pytest.approx(20.0, rel=1e-12)
    assert r[_gu(ids.mucosa_ps(S))] == pytest.approx(
        20.0 * r[_gu(ids.mucosa_volume(None))], rel=1e-12
    )


@pytest.mark.parametrize("single", [True, False])
def test_a_permeability_of_a_substance_without_a_mucosa_is_refused(
    single: bool, tmp_path: Path
) -> None:
    """A `mucosa_ps` of a substance which has no mucosa species names the reason.

    The other drug has no chemistry of the small intestine, so it has no
    mucosa species: the message names the substances which have one.
    """
    d = _definition(
        f"no_mucosa{int(single)}",
        CLINT_PERMEABILITY,
        mucosa_ps={OTHER: CLPERM},
        single=single,
    )
    with pytest.raises(ValueError, match=f"mucosa_ps: '{OTHER}' has no mucosa"):
        build(d, tmp_path)


def test_a_permeability_of_an_unknown_substance_is_refused(tmp_path: Path) -> None:
    """A `mucosa_ps` of an id which is no substance of the definition is refused."""
    d = _definition("unknown", CLINT_PERMEABILITY, mucosa_ps={"nope": CLPERM})
    with pytest.raises(ValueError, match="mucosa_ps: 'nope' is no substance"):
        build(d, tmp_path)


def test_a_permeability_which_is_no_mapping_is_refused() -> None:
    """The option is per substance: a single value for every substance is refused."""
    with pytest.raises(ValueError, match="mucosa_ps: a mapping from substance id"):
        SegmentedIntestine(mucosa_ps=1.0)  # ty: ignore[invalid-argument-type]


def test_the_intestine_with_the_option_is_hashable() -> None:
    """The option is stored as a read-only mapping, the dataclass stays hashable."""
    for intestine in (
        SegmentedIntestine(mucosa_ps={S: 1.0}),
        SingleLumen(mucosa_ps={S: 1.0}),
    ):
        assert hash(intestine)
        assert dict(intestine.mucosa_ps or {}) == {S: 1.0}
