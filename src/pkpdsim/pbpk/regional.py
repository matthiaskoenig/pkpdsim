"""Derive a regional absorption profile with the model from a release study.

A regional absorption study releases a solution of the drug in a segment of
the intestine (an intubation, a remote-controlled capsule) and reports its
exposure relative to an oral solution, e.g. the AUC of apixaban released in
the distal small bowel 0.40 and in the ascending colon 0.16 of the oral
solution (Byon 2018). The factors of a `RegionalAbsorption` profile which
give these ratios depend on the rest of the absorption of the model (`ka`,
`fabs`, an efflux, the first pass), so they are derived with the model:

- `release_auc_ratios` simulates the oral solution (`SOLDOSE_<s>`) and the
  same dose put dissolved into the lumen of each release segment
  (`GU__Clumen_<seg>_<s>`) and gives the AUC of the arm vein relative to the
  oral solution;
- `derive_regional_absorption` finds the factors for which the ratios are the
  reported ones: one factor per release segment, the factor of the profile
  which scales the absorption there (`regional_key`: the segment itself in
  the small intestine, the group `colon` in the large intestine), each the
  root of its ratio (`brentq`) with the others fixed, in turns until no
  factor moves (a release reaches the other regions little, so a few turns
  converge).

Because the ratios depend on the fitted absorption, a model derives the
factors on its adopted fit, writes them into its profile (a `Value` with the
source of the study) and fits again from the adopted fit, until the factors
derived on the fit are those of the model. The model must have the factors
(its profile names them, e.g. with a first guess) and an oral solution of the
substance.
"""

from collections.abc import Generator, Mapping, Sequence
from contextlib import contextmanager
from typing import TYPE_CHECKING

import numpy as np

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.organs import Organ
from pkpdsim.pbpk.variants.lumen_options import ABSORBING, COLON, LARGE

if TYPE_CHECKING:
    import roadrunner

__all__ = [
    "derive_regional_absorption",
    "regional_key",
    "release_auc_ratios",
]

DAY = 24 * 60.0
"""A day [min]."""


def regional_key(segment: str) -> str:
    """The key of the profile whose factor scales the absorption of a release segment.

    Args:
        segment: the token of an absorbing segment (`du` to `rec`).

    Returns:
        The segment itself in the small intestine, the group `colon` in the
        large intestine (caecum to rectum).

    Raises:
        ValueError: for the stomach or an unknown token.
    """
    if segment not in ABSORBING:
        raise ValueError(
            f"'{segment}' is no absorbing segment of the intestine, the segments "
            f"are {list(ABSORBING)}"
        )
    return COLON if segment in LARGE else segment


def _factor_id(substance: str, key: str) -> str:
    """The id of the factor of a key of the profile, `f_ka_<seg>_<s>` or `f_ka_colon_<s>`."""
    if key == COLON:
        return ids.colon_absorption_factor(substance)
    return ids.regional_absorption_factor(substance, key)


def _check(r: "roadrunner.RoadRunner", names: Sequence[str]) -> None:
    """Refuse a model without an element the derivation sets or reads.

    Raises:
        ValueError: naming the missing ids.
    """
    known = set(r.model.getGlobalParameterIds()) | set(r.model.getFloatingSpeciesIds())
    known |= set(r.model.getCompartmentIds())
    missing = sorted(set(names) - known)
    if missing:
        raise ValueError(
            f"The model has no {missing}: the derivation needs the oral solution "
            "of the substance, the lumen of the release segments and the factors "
            "of the profile (give the profile a first value of each)."
        )


@contextmanager
def _condition(
    r: "roadrunner.RoadRunner", changes: Mapping[str, float] | None
) -> Generator[None]:
    """Simulate in the condition of the study, then take it back.

    The changes are initial values (`init(<id>)`), so a reset keeps them and
    the initial assignments which read them follow (the volumes of a body
    weight, `T_st_eff` of the fed state); at the end the initial values of
    before are set again and every value is reset (`reset` of all values).
    Setting an initial value compiles the model again (seconds per change and
    call).
    """
    import roadrunner

    changes = dict(changes or {})
    _check(r, list(changes))
    before = {sid: r.getValue(f"init({sid})") for sid in changes}
    try:
        for sid, value in changes.items():
            r.setValue(f"init({sid})", float(value))
        r.reset(roadrunner.SelectionRecord.ALL)
        yield
    finally:
        for sid, value in before.items():
            r.setValue(f"init({sid})", value)
        r.reset(roadrunner.SelectionRecord.ALL)


def _ratios(
    r: "roadrunner.RoadRunner",
    substance: str,
    segments: Sequence[str],
    dose: float,
    factors: Mapping[str, float] | None,
    end: float,
    points: int,
    observable: str | None,
) -> dict[str, float]:
    """The ratios of `release_auc_ratios` in the condition which `r` is in."""
    import roadrunner

    observable = observable or ids.arm(substance)
    for seg in segments:
        regional_key(seg)
    lumen = {
        seg: (
            ids.flat(Organ.SMALL_INTESTINE, ids.lumen(substance, seg)),
            ids.flat(Organ.SMALL_INTESTINE, ids.lumen_volume(seg)),
        )
        for seg in segments
    }
    _check(
        r,
        [
            ids.sol_dose(substance),
            ids.mr(substance),
            *(x for pair in lumen.values() for x in pair),
            *(_factor_id(substance, k) for k in (factors or {})),
        ],
    )

    def start() -> None:
        r.reset(roadrunner.SelectionRecord.ALL)
        for key, x in (factors or {}).items():
            r[_factor_id(substance, key)] = x

    def area() -> float:
        s = r.simulate(0, end, points, selections=["time", observable])
        return float(np.trapezoid(s[:, 1], s[:, 0]))

    start()
    r[ids.sol_dose(substance)] = dose
    oral = area()
    ratios = {}
    for seg, (species, volume) in lumen.items():
        start()
        r[f"[{species}]"] = dose / r[ids.mr(substance)] / r[volume]
        ratios[seg] = area() / oral
    return ratios


def release_auc_ratios(
    r: "roadrunner.RoadRunner",
    substance: str,
    segments: Sequence[str],
    dose: float,
    *,
    factors: Mapping[str, float] | None = None,
    changes: Mapping[str, float] | None = None,
    end: float = 14 * DAY,
    points: int = 20161,
    observable: str | None = None,
) -> dict[str, float]:
    """Get the AUC of a solution released in each segment relative to the oral solution.

    Every simulation starts from the initial values of the model in the
    condition of the study (`changes`, set as initial values) with the
    factors set; the oral solution is the dose state `SOLDOSE_<s>`, a release
    the dose dissolved in the lumen of the segment (`dose/Mr/Vlumen`). The
    AUC is the trapezoid of `observable` over `points` equidistant times from
    0 to `end`. Afterwards the condition is taken back and every value of `r`
    is reset to its initial value (`reset` of all values): a value set on
    `r` before the call is kept only if it was set as an initial value.

    Args:
        r: the model, a roadrunner of the flat model with the tolerances of
            the derivation (e.g. absolute 1e-14, relative 1e-10).
        substance: the substance id.
        segments: the release segments (`ile2`, `asc`).
        dose: the dose [mg].
        factors: the factors of the profile by key (`regional_key`), those
            of the model if `None`.
        changes: the condition of the study, initial values by id (`FED`,
            `BW`, an activity); none if `None`.
        end: the end of the simulations [min], long enough for the AUC.
        points: the number of time points.
        observable: the id whose AUC is compared, the arm vein `Carm_<s>` by
            default.

    Returns:
        The ratio by release segment.

    Raises:
        ValueError: for a segment which does not absorb, or a model without
            the oral solution, a lumen, a factor or a changed id.
    """
    with _condition(r, changes):
        return _ratios(r, substance, segments, dose, factors, end, points, observable)


def derive_regional_absorption(
    r: "roadrunner.RoadRunner",
    substance: str,
    targets: Mapping[str, float],
    dose: float,
    *,
    changes: Mapping[str, float] | None = None,
    bounds: tuple[float, float] = (1e-4, 10.0),
    xtol: float = 1e-6,
    max_turns: int = 20,
    end: float = 14 * DAY,
    points: int = 20161,
    observable: str | None = None,
) -> dict[str, float]:
    """Derive the factors of a profile which give the ratios of a release study.

    One factor per release segment (`regional_key`: the segment in the small
    intestine, `colon` in the large one, so two release segments of the
    large intestine are refused), each the root of the ratio of its segment
    within `bounds` with the other factors fixed, in the order of `targets`
    and in turns, starting from the factors of the model, until no factor
    moves by more than `xtol`. Afterwards `r` is reset as by
    `release_auc_ratios`: the factors of the model, not the derived ones.

    Args:
        r: the model, see `release_auc_ratios`.
        substance: the substance id.
        targets: the reported AUC of the release relative to the oral
            solution by release segment (`{"ile2": 0.40, "asc": 0.16}`).
        dose: the dose of the study [mg].
        changes: the condition of the study, see `release_auc_ratios`.
        bounds: the bracket of every factor.
        xtol: the tolerance of a factor.
        max_turns: the most turns over the factors.
        end: see `release_auc_ratios`.
        points: see `release_auc_ratios`.
        observable: see `release_auc_ratios`.

    Returns:
        The factors by key of the profile (`{"ile2": 0.149, "colon": 0.0567}`).

    Raises:
        ValueError: for two release segments of one factor, a ratio which no
            factor within `bounds` gives, or a model without what the
            derivation needs.
        RuntimeError: if the factors do not converge in `max_turns` turns.
    """
    from scipy.optimize import brentq

    keys = {seg: regional_key(seg) for seg in targets}
    twice = sorted({k for k in keys.values() if list(keys.values()).count(k) > 1})
    if twice:
        raise ValueError(
            f"The release segments {sorted(targets)} share the factors {twice}: "
            "one release segment per factor of the profile."
        )
    _check(r, [_factor_id(substance, k) for k in keys.values()])
    with _condition(r, changes):
        factors = {k: float(r[_factor_id(substance, k)]) for k in keys.values()}

        def residual(x: float, seg: str) -> float:
            trial = {**factors, keys[seg]: x}
            ratio = _ratios(r, substance, [seg], dose, trial, end, points, observable)[
                seg
            ]
            return ratio - targets[seg]

        for _ in range(max_turns):
            moved = 0.0
            for seg in targets:
                low, high = residual(bounds[0], seg), residual(bounds[1], seg)
                if low * high > 0:
                    raise ValueError(
                        f"No factor of '{keys[seg]}' within {bounds} gives the "
                        f"ratio {targets[seg]} of '{seg}' (the ratios at the bounds "
                        f"{low + targets[seg]:.4g} and {high + targets[seg]:.4g})."
                    )
                x = float(brentq(residual, *bounds, args=(seg,), xtol=xtol))
                moved = max(moved, abs(x - factors[keys[seg]]))
                factors[keys[seg]] = x
            if moved <= xtol:
                return factors
    raise RuntimeError(
        f"The factors {factors} of '{substance}' did not converge in {max_turns} turns."
    )
