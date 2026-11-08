"""The integration of the studies of the database into a consensus with its uncertainty.

The consensus of a quantity of an entity (e.g. the protein abundance of
CYP3A4 in the liver) in a condition (e.g. cirrhosis of Child-Pugh class C) is
the ratio of the condition to the control, integrated over the studies which
report both.

1. The location of a group on the log scale, the log of its median (equal to
   the geometric mean for log-normal data, so means and medians of different
   studies estimate the same ratio), with its standard error, from the
   statistics its source reports, in this order:
   - the median with its 95 % confidence interval: `ln(median)`, the standard
     error `(ln(upper) - ln(lower))/(2*1.96)`;
   - the median with its quartiles (a box plot): the SD of the log of the
     interquartile range of a normal distribution of the log,
     `sigma = (ln(q3) - ln(q1))/1.349`, the standard error of the log of a
     median `sqrt(pi/2)*sigma/sqrt(n)`;
   - the median with the mean and SD: the log-normal SD of the log
     `sigma = sqrt(ln(1 + CV^2))` (CV = SD/mean), the standard error of the
     log of a median `sqrt(pi/2)*sigma/sqrt(n)`;
   - the mean and SD alone: the geometric mean of the log-normal mean and SD,
     `ln(mean) - sigma^2/2`, its standard error `sigma/sqrt(n)`.
   A group without its spread (SD, CV or interval) takes the `sigma` of its
   study and kind of measurement (the median over its other records), or of
   all studies of the kind, and is marked as imputed. A group whose location
   is not positive (a median or mean of 0, below the limit of quantification
   in most samples) has no log and is excluded with the reason.

   A group with many samples below the limit of quantification (counted as
   0) is no log-normal group: its median and quartiles sit on the limit and
   their ratio says nothing. Where a group of a comparison is such a group (a
   first quartile or a median of 0, or fewer than 3/4 of its samples
   detected), both groups of the comparison take their arithmetic means,
   `ln(mean)` with the standard error `CV/sqrt(n)` (the delta method, the CV
   of the SD over the mean, imputed as above without), the ratio of the
   means. A group with at most a quarter of its samples detected (a third
   quartile of 0) says nothing of a ratio: its comparison is excluded.
2. The estimate of a study: the groups of the study which match the
   condition (every field of the condition equal), the group without an
   etiology (or `mixed`) if the study has one, else every matching group,
   strata of disjoint subjects pooled by their inverse variance; minus the
   control group of the study. Its variance is the sum of both.
3. The consensus of the estimates of the studies: the random-effects
   meta-analysis of DerSimonian and Laird, the between-study variance `tau^2`
   from Cochran's Q, so heterogeneity between the studies (populations,
   etiologies, methods) widens the consensus. A single study has no `tau^2`.

The consensus is a log-normal factor: `value = exp(mu)`, the 95 % confidence
interval of the mean ratio `exp(mu +- 1.96*se)` and the uncertainty of the
ratio of a new population `sqrt(se^2 + tau^2)` (the prediction).
"""

import math
import statistics
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from functools import cache

from pkpdsim.database.studies import KINDS, Record, Study, Summary, records, studies

__all__ = [
    "Consensus",
    "Estimate",
    "Exclusion",
    "Level",
    "LevelEstimate",
    "LogLocation",
    "Z",
    "consensus",
    "entities",
    "factor",
    "level",
    "log_location",
    "mean_location",
    "mostly_undetected",
    "random_effects",
    "sigma_log",
    "study_estimate",
    "zero_inflated",
]

Z = 1.959963984540054
"""The 97.5 % quantile of the standard normal distribution."""

_Z_QUARTILE = 0.6744897501960817
"""The 75 % quantile of the standard normal distribution."""

_MIXED = (None, "mixed")
"""The etiologies of a group which is no stratum of an etiology."""


@dataclass(frozen=True)
class LogLocation:
    """The location of a group on the log scale.

    Attributes:
        y: the log of the median (or of the geometric mean).
        se: the standard error of `y`.
        statistic: the statistics `y` and `se` are taken from.
        imputed: the spread is the spread of the other records of its study
            or kind, the group reports none.
    """

    y: float
    se: float
    statistic: str
    imputed: bool = False


@dataclass(frozen=True)
class Exclusion:
    """A record which no estimate uses, with the reason.

    Attributes:
        study: the key of the study.
        group: the key of the group.
        reason: why it is excluded.
    """

    study: str
    group: str
    reason: str


@dataclass(frozen=True)
class Estimate:
    """The estimate of the log ratio of a condition to the control of a study.

    Attributes:
        study: the key of the study.
        kind: the kind of measurement.
        groups: the keys of the groups of the condition.
        control: the key of the control group.
        log_ratio: the log of the ratio to the control.
        se: its standard error.
        statistic: the statistics of the locations.
        imputed: a spread is imputed.
        n: the samples of the groups of the condition.
        n_control: the samples of the control group.
        assumptions: the assumptions on the groups, e.g. their class.
    """

    study: str
    kind: str
    groups: tuple[str, ...]
    control: str
    log_ratio: float
    se: float
    statistic: str
    imputed: bool
    n: int | None
    n_control: int | None
    assumptions: tuple[str, ...] = ()

    @property
    def ratio(self) -> float:
        """The ratio to the control."""
        return math.exp(self.log_ratio)


@dataclass(frozen=True)
class Consensus:
    """The consensus of the studies on the ratio of a condition to the control.

    Attributes:
        kind: the kind of measurement.
        entity: the gene symbol of the entity.
        tissue: the tissue.
        condition: the condition, `(field, value)` pairs.
        estimates: the estimates of the studies.
        excluded: the records which no estimate uses.
        mu: the log of the consensus ratio.
        se: the standard error of `mu`.
        tau2: the between-study variance, 0 for a single study.
        q: Cochran's Q, `None` for a single study.
        i2: the share of the variance between the studies, `None` for a single
            study.
    """

    kind: str
    entity: str
    tissue: str
    condition: tuple[tuple[str, str], ...]
    estimates: tuple[Estimate, ...]
    excluded: tuple[Exclusion, ...]
    mu: float
    se: float
    tau2: float
    q: float | None
    i2: float | None

    @property
    def value(self) -> float:
        """The consensus ratio."""
        return math.exp(self.mu)

    @property
    def ci(self) -> tuple[float, float]:
        """The 95 % confidence interval of the consensus ratio."""
        return math.exp(self.mu - Z * self.se), math.exp(self.mu + Z * self.se)

    @property
    def prediction_se(self) -> float:
        """The standard error of the log ratio of a new population, `sqrt(se^2 + tau^2)`."""
        return math.sqrt(self.se**2 + self.tau2)

    @property
    def sources(self) -> tuple[str, ...]:
        """The keys of the studies of the estimates."""
        return tuple(e.study for e in self.estimates)

    def text(self) -> str:
        """A short statement of the consensus with its interval and studies."""
        lo, hi = self.ci
        parts = [
            f"{self.value:.3g} (95 % CI {lo:.3g}-{hi:.3g}; "
            f"{len(self.estimates)} {'study' if len(self.estimates) == 1 else 'studies'}: "
            + ", ".join(f"{e.study} {e.ratio:.3g}" for e in self.estimates)
            + ")"
        ]
        if self.kind != "protein":
            parts.append(f"from {KINDS[self.kind]}")
        if any(e.imputed for e in self.estimates):
            parts.append("a spread imputed")
        return ", ".join(parts)


def sigma_log(summary: Summary) -> float | None:
    """The SD of the log of a log-normal group of its mean and SD (or CV), or its quartiles.

    Args:
        summary: the statistics of the group.

    Returns:
        `sqrt(ln(1 + CV^2))` of the mean and SD (the SD of the standard error
        of the mean times the root of n, `se*sqrt(n)`, for a group which
        reports the SE alone) or the CV, else `(ln(q3) - ln(q1))/1.349` of
        positive quartiles; `None` without a positive spread.
    """
    cv = None
    sd = summary.sd
    if not sd and summary.se and summary.n:
        sd = summary.se * math.sqrt(summary.n)
    if summary.mean and summary.mean > 0 and sd:
        cv = sd / summary.mean
    elif summary.cv:
        cv = summary.cv
    if cv:
        return math.sqrt(math.log1p(cv**2))
    q1, q3 = summary.q1, summary.q3
    if q1 is not None and q3 is not None and q1 > 0 and q3 > q1:
        return (math.log(q3) - math.log(q1)) / (2 * _Z_QUARTILE)
    return None


@cache
def _sigmas(study: str | None, measurement: str | None, kind: str) -> float | None:
    """The median SD of the log of the records of a measurement, or of a kind."""
    found = [
        s
        for r in records(kind=kind)
        if (study is None or r.study == study)
        and (measurement is None or r.measurement == measurement)
        and (s := sigma_log(r.summary)) is not None
    ]
    return statistics.median(found) if found else None


def _fallback_sigma(record: Record) -> float | None:
    """The SD of the log for a record without a spread: of its measurement, else of its kind."""
    return _sigmas(record.study, record.measurement, record.kind) or _sigmas(
        None, None, record.kind
    )


def log_location(record: Record) -> LogLocation | str:
    """Get the location of a record on the log scale with its standard error.

    Args:
        record: the record of a group.

    Returns:
        The location, or the reason why the record has none.
    """
    s = record.summary
    if s.median is not None and s.median > 0 and s.ci is not None and s.ci[0] > 0:
        return LogLocation(
            y=math.log(s.median),
            se=(math.log(s.ci[1]) - math.log(s.ci[0])) / (2 * Z),
            statistic="median and its 95 % CI",
        )
    n = s.n
    if not n:
        return "no number of samples"
    if (
        s.median is not None
        and s.median > 0
        and s.q1 is not None
        and s.q3 is not None
        and s.q1 > 0
        and s.q3 > s.q1
    ):
        sigma_iqr = (math.log(s.q3) - math.log(s.q1)) / (2 * _Z_QUARTILE)
        return LogLocation(
            y=math.log(s.median),
            se=math.sqrt(math.pi / 2) * sigma_iqr / math.sqrt(n),
            statistic="median with the log-normal SD of its quartiles",
        )
    median = s.median if s.median is not None and s.median > 0 else None
    mean = s.mean if s.mean is not None and s.mean > 0 else None
    if median is None and mean is None:
        return "no positive median or mean (below the limit of quantification)"
    sigma = sigma_log(s)
    imputed = False
    if sigma is None:
        sigma = _fallback_sigma(record)
        imputed = True
        if sigma is None:
            return "no spread and none to impute"
    if median is not None:
        return LogLocation(
            y=math.log(median),
            se=math.sqrt(math.pi / 2) * sigma / math.sqrt(n),
            statistic="median with the log-normal SD of mean and SD",
            imputed=imputed,
        )
    assert mean is not None
    return LogLocation(
        y=math.log(mean) - sigma**2 / 2,
        se=sigma / math.sqrt(n),
        statistic="geometric mean of the log-normal mean and SD",
        imputed=imputed,
    )


def mostly_undetected(summary: Summary) -> bool:
    """Whether at most a quarter of the samples of a group are detected.

    Args:
        summary: the statistics of the group.

    Returns:
        A third quartile of 0, or fewer than a quarter of the samples detected.
    """
    if summary.q3 is not None and summary.q3 <= 0:
        return True
    return (
        summary.n is not None
        and summary.n_detected is not None
        and summary.n_detected < 0.25 * summary.n
    )


def zero_inflated(summary: Summary) -> bool:
    """Whether a group has many samples below the limit of quantification.

    Args:
        summary: the statistics of the group.

    Returns:
        A first quartile or a median of 0, or fewer than 3/4 of the samples
        detected.
    """
    if summary.q1 is not None and summary.q1 <= 0:
        return True
    if summary.median is not None and summary.median <= 0:
        return True
    return (
        summary.n is not None
        and summary.n_detected is not None
        and summary.n_detected < 0.75 * summary.n
    )


@cache
def _cvs(study: str | None, measurement: str | None, kind: str) -> float | None:
    """The median CV (SD over mean) of the records of a measurement, or of a kind."""
    found = [
        r.summary.sd / r.summary.mean
        for r in records(kind=kind)
        if (study is None or r.study == study)
        and (measurement is None or r.measurement == measurement)
        and r.summary.mean
        and r.summary.mean > 0
        and r.summary.sd
    ]
    return statistics.median(found) if found else None


def mean_location(record: Record) -> LogLocation | str:
    """Get the log of the arithmetic mean of a record with its standard error.

    Args:
        record: the record of a group.

    Returns:
        `ln(mean)` with the standard error `CV/sqrt(n)` (the CV imputed from
        its measurement or kind without an SD), or the reason why the record
        has none.
    """
    s = record.summary
    if not s.n:
        return "no number of samples"
    if s.mean is None or s.mean <= 0:
        return "no positive mean"
    imputed = False
    if s.sd:
        cv = s.sd / s.mean
    elif s.cv:
        cv = s.cv
    else:
        found = _cvs(record.study, record.measurement, record.kind) or _cvs(
            None, None, record.kind
        )
        if found is None:
            return "no spread and none to impute"
        cv = found
        imputed = True
    return LogLocation(
        y=math.log(s.mean),
        se=cv / math.sqrt(s.n),
        statistic="arithmetic mean (samples below the limit of quantification)",
        imputed=imputed,
    )


def _matches(fields: Mapping[str, str], condition: Mapping[str, str]) -> bool:
    """Whether the fields of a group have every field of the condition."""
    return all(fields.get(k) == v for k, v in condition.items())


def _pooled(locations: Sequence[LogLocation]) -> LogLocation:
    """The inverse-variance pooled location of strata of disjoint subjects."""
    if len(locations) == 1:
        return locations[0]
    weights = [1 / loc.se**2 for loc in locations]
    total = sum(weights)
    return LogLocation(
        y=sum(w * loc.y for w, loc in zip(weights, locations, strict=True)) / total,
        se=math.sqrt(1 / total),
        statistic="; ".join(sorted({loc.statistic for loc in locations})),
        imputed=any(loc.imputed for loc in locations),
    )


def study_estimate(
    study: Study, kind: str, entity: str, condition: Mapping[str, str]
) -> tuple[Estimate | None, list[Exclusion]]:
    """Get the estimate of a study of the log ratio of a condition to its control.

    Args:
        study: the study.
        kind: the kind of measurement.
        entity: the gene symbol.
        condition: the fields the groups must have, e.g.
            `{"condition": "cirrhosis", "child_pugh": "C"}`.

    Returns:
        The estimate, `None` if the study has no usable group of the
        condition or no usable control, and the records it excludes.
    """
    rows = [
        r
        for m in study.measurements.values()
        if m.kind == kind
        for r in m.records
        if r.entity == entity
    ]
    control = next((r for r in rows if r.group.id == study.control), None)
    matching = [
        r
        for r in rows
        if r.group.id != study.control and _matches(r.group.fields(), condition)
    ]
    if not matching or control is None:
        return None, []
    aggregate = [r for r in matching if r.group.etiology in _MIXED]
    chosen = aggregate or matching
    undetected = [r for r in [*chosen, control] if mostly_undetected(r.summary)]
    if undetected:
        return None, [
            Exclusion(
                study.key, r.group.id, "at most a quarter of the samples detected"
            )
            for r in undetected
        ]
    by_mean = any(zero_inflated(r.summary) for r in [*chosen, control])
    locate = mean_location if by_mean else log_location
    excluded: list[Exclusion] = []
    locations: list[LogLocation] = []
    used: list[Record] = []
    for r in chosen:
        loc = locate(r)
        if isinstance(loc, str):
            excluded.append(Exclusion(study.key, r.group.id, loc))
            continue
        locations.append(loc)
        used.append(r)
    ref = locate(control)
    if isinstance(ref, str):
        excluded.append(Exclusion(study.key, control.group.id, ref))
        return None, excluded
    if not locations:
        return None, excluded
    group = _pooled(locations)
    ns = [r.summary.n for r in used]
    return (
        Estimate(
            study=study.key,
            kind=kind,
            groups=tuple(r.group.id for r in used),
            control=control.group.id,
            log_ratio=group.y - ref.y,
            se=math.sqrt(group.se**2 + ref.se**2),
            statistic=group.statistic,
            imputed=group.imputed or ref.imputed,
            n=sum(ns) if all(n is not None for n in ns) else None,  # ty: ignore[no-matching-overload]
            n_control=control.summary.n,
            assumptions=tuple(r.group.assumption for r in used if r.group.assumption),
        ),
        excluded,
    )


def random_effects(
    estimates: Sequence[Estimate],
) -> tuple[float, float, float, float | None, float | None]:
    """Pool estimates by the random-effects meta-analysis of DerSimonian and Laird.

    Args:
        estimates: the estimates of the studies, at least one.

    Returns:
        `(mu, se, tau2, q, i2)`: the pooled log ratio, its standard error, the
        between-study variance, Cochran's Q and I^2 (`None` for one study).

    Raises:
        ValueError: without an estimate.
    """
    if not estimates:
        raise ValueError("No estimate to pool.")
    y = [e.log_ratio for e in estimates]
    v = [e.se**2 for e in estimates]
    if len(estimates) == 1:
        return y[0], math.sqrt(v[0]), 0.0, None, None
    w = [1 / vi for vi in v]
    sw = sum(w)
    fixed = sum(wi * yi for wi, yi in zip(w, y, strict=True)) / sw
    q = sum(wi * (yi - fixed) ** 2 for wi, yi in zip(w, y, strict=True))
    df = len(estimates) - 1
    c = sw - sum(wi**2 for wi in w) / sw
    tau2 = max(0.0, (q - df) / c)
    ws = [1 / (vi + tau2) for vi in v]
    mu = sum(wi * yi for wi, yi in zip(ws, y, strict=True)) / sum(ws)
    se = math.sqrt(1 / sum(ws))
    i2 = max(0.0, (q - df) / q) if q > 0 else 0.0
    return mu, se, tau2, q, i2


def _key(condition: Mapping[str, str]) -> tuple[tuple[str, str], ...]:
    """The hashable key of a condition."""
    return tuple(sorted(condition.items()))


@cache
def _consensus(
    kind: str, entity: str, tissue: str, condition: tuple[tuple[str, str], ...]
) -> Consensus | None:
    """The cached consensus of a hashable condition."""
    estimates: list[Estimate] = []
    excluded: list[Exclusion] = []
    for study in studies().values():
        if study.tissue != tissue:
            continue
        estimate, out = study_estimate(study, kind, entity, dict(condition))
        excluded += out
        if estimate is not None:
            estimates.append(estimate)
    if not estimates:
        return None
    mu, se, tau2, q, i2 = random_effects(estimates)
    return Consensus(
        kind=kind,
        entity=entity,
        tissue=tissue,
        condition=condition,
        estimates=tuple(estimates),
        excluded=tuple(excluded),
        mu=mu,
        se=se,
        tau2=tau2,
        q=q,
        i2=i2,
    )


def consensus(
    entity: str,
    condition: Mapping[str, str],
    *,
    kind: str = "protein",
    tissue: str = "liver",
) -> Consensus | None:
    """Get the consensus of the studies on the ratio of a condition to the control.

    Args:
        entity: the gene symbol, e.g. `CYP3A4`.
        condition: the fields of the groups, e.g.
            `{"condition": "cirrhosis", "child_pugh": "C"}`.
        kind: the kind of measurement, a key of `KINDS`.
        tissue: the tissue.

    Returns:
        The consensus, `None` if no study reports the condition and its control.

    Raises:
        ValueError: for an unknown kind.
    """
    if kind not in KINDS:
        raise ValueError(f"Unknown kind '{kind}', the kinds are {sorted(KINDS)}.")
    return _consensus(kind, entity, tissue, _key(condition))


def factor(
    entity: str,
    condition: Mapping[str, str],
    *,
    tissue: str = "liver",
    kinds: Sequence[str] = ("protein", "mrna"),
) -> Consensus | None:
    """Get the consensus factor of a condition, the protein first.

    Args:
        entity: the gene symbol.
        condition: the fields of the groups.
        tissue: the tissue.
        kinds: the kinds in their order of preference: the mRNA expression is
            used only without a protein abundance of the condition.

    Returns:
        The consensus of the first kind with an estimate, `None` without.
    """
    for kind in kinds:
        found = consensus(entity, condition, kind=kind, tissue=tissue)
        if found is not None:
            return found
    return None


def entities(
    kind: str = "protein", tissue: str = "liver", condition: str = "cirrhosis"
) -> list[str]:
    """Get the entities with records of a condition.

    Args:
        kind: the kind of measurement.
        tissue: the tissue.
        condition: the condition of the groups.

    Returns:
        The gene symbols, sorted.
    """
    return sorted(
        {
            r.entity
            for r in records(kind=kind, tissue=tissue)
            if r.group.condition == condition
        }
    )


@dataclass(frozen=True)
class LevelEstimate:
    """The level of a quantity in the groups of a condition of a study.

    Attributes:
        study: the key of the study.
        groups: the keys of the groups.
        value: the level, the mean (or the median) of the groups.
        se: its standard error.
        statistic: the statistics of the level.
        n: the samples of the groups.
        unit: the unit of the level.
    """

    study: str
    groups: tuple[str, ...]
    value: float
    se: float
    statistic: str
    n: int | None
    unit: str


@dataclass(frozen=True)
class Level:
    """The consensus of the studies on the level of a quantity in a condition.

    Attributes:
        kind: the kind of measurement.
        entity: the entity, a gene symbol or a quantity.
        tissue: the tissue.
        condition: the condition, `(field, value)` pairs.
        unit: the unit of the levels.
        estimates: the levels of the studies.
        value: the pooled level.
        se: its standard error.
        tau2: the between-study variance, 0 for a single study.
    """

    kind: str
    entity: str
    tissue: str
    condition: tuple[tuple[str, str], ...]
    unit: str
    estimates: tuple[LevelEstimate, ...]
    value: float
    se: float
    tau2: float

    @property
    def ci(self) -> tuple[float, float]:
        """The 95 % confidence interval of the pooled level."""
        return self.value - Z * self.se, self.value + Z * self.se

    @property
    def sources(self) -> tuple[str, ...]:
        """The keys of the studies of the estimates."""
        return tuple(e.study for e in self.estimates)

    def text(self) -> str:
        """A short statement of the level with its interval and studies."""
        lo, hi = self.ci
        n = len(self.estimates)
        return (
            f"{self.value:.3g} {self.unit} (95 % CI {lo:.3g}-{hi:.3g}; {n} "
            f"{'study' if n == 1 else 'studies'}: "
            + ", ".join(f"{e.study} {e.value:.3g}" for e in self.estimates)
            + ")"
        )


def _level_of(record: Record) -> tuple[float, float, str] | None:
    """The level of a group with its standard error, `None` without a spread or n.

    The mean with its SE (or SD over the root of n); else the median with the
    SE of a median of the SD of its quartiles, `sqrt(pi/2)*(q3 - q1)/1.349/sqrt(n)`.
    """
    s = record.summary
    if s.mean is not None:
        if s.se:
            return s.mean, s.se, "mean and SE"
        if s.sd is not None and s.n:
            return s.mean, max(s.sd, 0.0) / math.sqrt(s.n), "mean and SD"
    if s.median is not None and s.q1 is not None and s.q3 is not None and s.n:
        sd = (s.q3 - s.q1) / (2 * _Z_QUARTILE)
        return (
            s.median,
            math.sqrt(math.pi / 2) * sd / math.sqrt(s.n),
            "median and quartiles",
        )
    return None


@cache
def _level(
    kind: str, entity: str, tissue: str, condition: tuple[tuple[str, str], ...]
) -> Level | None:
    """The cached level of a hashable condition."""
    estimates: list[LevelEstimate] = []
    units: set[str] = set()
    for study in studies().values():
        if study.tissue != tissue:
            continue
        rows = [
            r
            for m in study.measurements.values()
            if m.kind == kind
            for r in m.records
            if r.entity == entity and _matches(r.group.fields(), dict(condition))
        ]
        aggregate = [r for r in rows if r.group.etiology in _MIXED]
        chosen = aggregate or rows
        found = [
            (r, lv) for r in chosen if (lv := _level_of(r)) is not None and lv[1] > 0
        ]
        if not found:
            continue
        weights = [1 / lv[1] ** 2 for _, lv in found]
        total = sum(weights)
        value = (
            sum(w * lv[0] for w, (_, lv) in zip(weights, found, strict=True)) / total
        )
        ns = [r.summary.n for r, _ in found]
        units |= {r.unit for r, _ in found}
        estimates.append(
            LevelEstimate(
                study=study.key,
                groups=tuple(r.group.id for r, _ in found),
                value=value,
                se=math.sqrt(1 / total),
                statistic="; ".join(sorted({lv[2] for _, lv in found})),
                n=sum(n for n in ns if n is not None) if all(ns) else None,
                unit=found[0][0].unit,
            )
        )
    if not estimates:
        return None
    if len(units) > 1:
        raise ValueError(f"The levels of {entity} have several units {sorted(units)}.")
    y = [e.value for e in estimates]
    v = [e.se**2 for e in estimates]
    if len(estimates) == 1:
        mu, se, tau2 = y[0], math.sqrt(v[0]), 0.0
    else:
        w = [1 / vi for vi in v]
        sw = sum(w)
        fixed = sum(wi * yi for wi, yi in zip(w, y, strict=True)) / sw
        q = sum(wi * (yi - fixed) ** 2 for wi, yi in zip(w, y, strict=True))
        c = sw - sum(wi**2 for wi in w) / sw
        tau2 = max(0.0, (q - (len(y) - 1)) / c)
        ws = [1 / (vi + tau2) for vi in v]
        mu = sum(wi * yi for wi, yi in zip(ws, y, strict=True)) / sum(ws)
        se = math.sqrt(1 / sum(ws))
    return Level(
        kind=kind,
        entity=entity,
        tissue=tissue,
        condition=condition,
        unit=units.pop(),
        estimates=tuple(estimates),
        value=mu,
        se=se,
        tau2=tau2,
    )


def level(
    entity: str,
    condition: Mapping[str, str],
    *,
    kind: str = "physiology",
    tissue: str = "body",
) -> Level | None:
    """Get the consensus of the studies on the level of a quantity in a condition.

    The level of a study is the mean (or the median) of its group of the
    condition, its strata pooled by their inverse variance; the levels of the
    studies are pooled by the random-effects meta-analysis of DerSimonian and
    Laird on the linear scale, e.g. the MELD score of the Child-Pugh classes.

    Args:
        entity: the entity, a quantity of the kind `physiology` (e.g. `meld`).
        condition: the fields of the groups, e.g. `{"condition": "control"}`.
        kind: the kind of measurement.
        tissue: the tissue, `body` for the quantities of the subjects.

    Returns:
        The level, `None` if no study reports the quantity in the condition.

    Raises:
        ValueError: for an unknown kind, or levels of several units.
    """
    if kind not in KINDS:
        raise ValueError(f"Unknown kind '{kind}', the kinds are {sorted(KINDS)}.")
    return _level(kind, entity, tissue, _key(condition))
