"""The declaration of the contexts of use of a model.

A `Credibility` lists the `ContextOfUse`s of a model. A context of use asks a
question of the model, names the quantities and the domain it is asked for,
the evidence (a selection of the fit mappings with a role) and, optionally,
the rules by which it is accepted. A context of use does not raise: it reports
`problems()`, and `Credibility` raises one `DeclarationError` with the
problems of every context.
"""

from __future__ import annotations

import math
import re
from collections.abc import Collection, Mapping
from dataclasses import dataclass, field
from enum import StrEnum, unique
from types import MappingProxyType
from typing import Final, Protocol, get_type_hints

import pandas as pd
from sbmlsim.fit.evaluation import Criterion, afe, criteria_metrics, gmfe
from sbmlsim.fit.helpers import MappingFilter

from pkpdsim.metadata import PKPDMappingMetaData

__all__ = [
    "DEFAULT_COVARIATES",
    "DOMAIN_FIELDS",
    "GUEST_DELTA",
    "UNCLAIMED",
    "Acceptance",
    "AcceptanceResult",
    "AcceptanceRule",
    "AfeRange",
    "Consequence",
    "ContextOfUse",
    "CouLevel",
    "Credibility",
    "DeclarationError",
    "Domain",
    "Evidence",
    "EvidenceRole",
    "GuestAcceptance",
    "Influence",
    "Level",
    "MaxGmfe",
    "QuantityOfInterest",
    "Ratio",
    "Risk",
    "Selection",
    "Summary",
    "Unclaimed",
]

#: the covariates of the subjects whose simulated values are checked against
#: the reported ones: the target of the model and the measurement type of PK-DB,
#: see `pkpdsim.credibility.inputs`
DEFAULT_COVARIATES: Final[Mapping[str, str]] = MappingProxyType(
    {"BW": "weight", "HEIGHT": "height"}
)

#: the pattern of the id of a context of use
ID_PATTERN: Final = re.compile(r"^[a-z][a-z0-9_]*$")

_COU_LABELS: Final = {
    "C1": "describe",
    "C2": "interpolate",
    "C3": "extrapolate with data",
    "C4": "extrapolate without data",
    "C5": "clinical endpoint",
}


@unique
class CouLevel(StrEnum):
    """The level of a context of use on the ladder from describing to predicting."""

    C1 = "C1"
    C2 = "C2"
    C3 = "C3"
    C4 = "C4"
    C5 = "C5"

    @property
    def label(self) -> str:
        """Get the name of the level, e.g. `interpolate`."""
        return _COU_LABELS[self.value]

    @property
    def decisive_role(self) -> EvidenceRole:
        """Get the role of the evidence which decides: calibration for C1, else validation."""
        if self is CouLevel.C1:
            return EvidenceRole.CALIBRATION
        return EvidenceRole.VALIDATION


@unique
class Level(StrEnum):
    """A level of influence, consequence or risk."""

    LOW = "low"
    MEDIUM = "medium"
    HIGH = "high"


Influence = Level
Consequence = Level
Risk = Level


@unique
class EvidenceRole(StrEnum):
    """What evidence is used for."""

    CALIBRATION = "calibration"
    VALIDATION = "validation"
    LIMITATION = "limitation"


class DeclarationError(ValueError):
    """The declaration of the credibility of a model is invalid.

    Attributes:
        problems: every problem of the declaration.
    """

    def __init__(self, problems: list[str]) -> None:
        """Create the error.

        Args:
            problems: every problem of the declaration.
        """
        self.problems = list(problems)
        super().__init__(
            "The credibility declaration is invalid:\n"
            + "\n".join(f"- {problem}" for problem in self.problems)
        )


@dataclass(frozen=True)
class Unclaimed:
    """The validation points which no other evidence and no limitation selects.

    Attributes:
        within: restricts the points to the ones this filter selects.
    """

    within: MappingFilter | None = None


UNCLAIMED: Final = Unclaimed()

type Selection = MappingFilter | Unclaimed


@unique
class Summary(StrEnum):
    """How the points of an arm of a ratio are summarized, see `Ratio`."""

    AUC = "auc"
    VALUE = "value"


#: the accepted variability of a ratio of 1 of the limits of Guest et al. (2011)
GUEST_DELTA: Final = 1.25


@dataclass(frozen=True, kw_only=True)
class Ratio:
    """The arms of a ratio quantity: a numerator over its reference in a study.

    See `pkpdsim.credibility.ratios`: per fit and study (the experiment) every
    decisive mapping of the observable which `numerator` selects is paired
    with the one which `denominator` selects, and the ratio of their summaries
    is compared for the data and the predictions.

    Attributes:
        observable: the observable of the arms, e.g. `Carm_hctz`; compared
            without the brackets of a concentration.
        numerator: selects the mappings of the arms, e.g. the renal impairment.
        denominator: selects the mapping of the reference arm, e.g. healthy.
        summary: how the points of an arm are summarized.
    """

    observable: str
    numerator: MappingFilter
    denominator: MappingFilter
    summary: Summary = Summary.AUC


@dataclass(frozen=True)
class QuantityOfInterest:
    """A quantity a context of use is asked for.

    Attributes:
        id: the observable of the fit mappings, e.g. `Cve_hctz` or `pk.cmax`;
            for a ratio its name, it is matched against no observable.
        description: what the quantity is.
        ratio: whether the quantity is a ratio between arms: `True` is
            declared and recorded, not assessed; a `Ratio` pairs the arms and
            is assessed with the limits of Guest et al., its observable
            restricts the evidence like the id of an absolute quantity.
    """

    id: str
    description: str
    ratio: bool | Ratio = field(default=False, kw_only=True)

    @property
    def pairing(self) -> Ratio | None:
        """Get the pairing of an assessed ratio, `None` otherwise."""
        return self.ratio if isinstance(self.ratio, Ratio) else None


#: the fields of the metadata which a domain can claim
DOMAIN_FIELDS: Final = (
    "tissue",
    "route",
    "application_form",
    "dosing",
    "fasting",
    "health",
    "coadministration",
    "genotype",
    "phenotype",
)


def _domain_fields() -> dict[str, type]:
    """Get the fields of the metadata which a domain can claim, with their enum."""
    hints = get_type_hints(PKPDMappingMetaData)
    return {name: hints[name] for name in DOMAIN_FIELDS}


@dataclass(frozen=True)
class Domain:
    """The domain a context of use claims.

    Attributes:
        metadata: the claimed values by field of `PKPDMappingMetaData`.
        statement: what is claimed beyond the metadata, e.g. the dose range.
    """

    metadata: Mapping[str, Collection[StrEnum]] = field(
        default_factory=dict, kw_only=True
    )
    statement: str = field(default="", kw_only=True)

    def __post_init__(self) -> None:
        """Store the values of every field as a tuple."""
        object.__setattr__(
            self,
            "metadata",
            {name: tuple(values) for name, values in self.metadata.items()},
        )

    def problems(self) -> list[str]:
        """Get the problems of the domain: an unknown field or a value of another vocabulary."""
        fields = _domain_fields()
        found: list[str] = []
        for name, values in self.metadata.items():
            enum = fields.get(name)
            if enum is None:
                found.append(
                    f"the domain field '{name}' is not a field of the metadata, "
                    f"use one of {', '.join(fields)}"
                )
                continue
            found.extend(
                f"the domain value {value!r} of '{name}' is not a member of "
                f"{enum.__name__}"
                for value in values
                if not isinstance(value, enum)
            )
        return found


@dataclass(frozen=True)
class Evidence:
    """A selection of the data with a role.

    Attributes:
        name: the name of the evidence.
        role: what the evidence is used for.
        selection: the fit mappings it selects.
        description: what the evidence is.
        fit: the subset of the adopted fits it is taken from, `None` for all.
    """

    name: str
    role: EvidenceRole
    selection: Selection
    description: str = ""
    fit: str | None = field(default=None, kw_only=True)


@dataclass(frozen=True)
class AcceptanceResult:
    """The result of an acceptance rule.

    Attributes:
        name: the name of the rule.
        value: the value of the rule on the points, `NaN` if not evaluable.
        limit: the limit, e.g. `>= 80`.
        met: whether the rule is met, `None` if the value is `NaN`.
    """

    name: str
    value: float
    limit: str
    met: bool | None


class AcceptanceRule(Protocol):
    """A rule on the data points of the evidence of a context of use."""

    @property
    def name(self) -> str:
        """Get the name of the rule."""
        ...

    @property
    def limit(self) -> str:
        """Get the limit of the rule, e.g. `>= 80`."""
        ...

    def evaluate(self, points: pd.DataFrame) -> AcceptanceResult:
        """Evaluate the rule on the data points (`DV`, `IPRED`, `DV_SD`)."""
        ...


def _result(name: str, value: float, limit: str, met: bool) -> AcceptanceResult:
    """Get the result of a rule, not evaluable for `NaN`."""
    return AcceptanceResult(name, value, limit, None if math.isnan(value) else met)


@dataclass(frozen=True)
class Acceptance:
    """The percentage of the points a criterion accepts is at least `minimum`.

    Attributes:
        criterion: the criterion, e.g. `Fold(2.0)`.
        minimum: the minimum percentage (0 to 100).
    """

    criterion: Criterion
    minimum: float

    @property
    def name(self) -> str:
        """Get the name of the criterion."""
        return self.criterion.name

    @property
    def limit(self) -> str:
        """Get the limit, e.g. `>= 80`."""
        return f">= {self.minimum:g}"

    def evaluate(self, points: pd.DataFrame) -> AcceptanceResult:
        """Evaluate the percentage of the points the criterion accepts."""
        value = criteria_metrics(points, (self.criterion,))[self.criterion.name]
        return _result(self.name, value, self.limit, value >= self.minimum)


@dataclass(frozen=True)
class MaxGmfe:
    """The geometric mean fold error is at most `maximum`.

    Attributes:
        maximum: the largest accepted GMFE.
    """

    maximum: float

    @property
    def name(self) -> str:
        """Get the name of the rule."""
        return "GMFE"

    @property
    def limit(self) -> str:
        """Get the limit, e.g. `<= 1.5`."""
        return f"<= {self.maximum:g}"

    def evaluate(self, points: pd.DataFrame) -> AcceptanceResult:
        """Evaluate the GMFE of the points."""
        value = gmfe(points)
        return _result(self.name, value, self.limit, value <= self.maximum)


@dataclass(frozen=True)
class AfeRange:
    """The average fold error, the bias, is inside the range.

    Attributes:
        low: the smallest accepted AFE.
        high: the largest accepted AFE.
    """

    low: float
    high: float

    @property
    def name(self) -> str:
        """Get the name of the rule."""
        return "AFE"

    @property
    def limit(self) -> str:
        """Get the limit, e.g. `0.8 to 1.25`."""
        return f"{self.low:g} to {self.high:g}"

    def evaluate(self, points: pd.DataFrame) -> AcceptanceResult:
        """Evaluate the AFE of the points, the bounds included."""
        value = afe(points)
        return _result(self.name, value, self.limit, self.low <= value <= self.high)


@dataclass(frozen=True)
class GuestAcceptance:
    """The percentage of the ratio pairs within the limits of Guest et al. is at least `minimum`.

    Attributes:
        minimum: the minimum percentage (0 to 100).
        delta: the accepted variability at a ratio of 1, see
            `pkpdsim.credibility.ratios.guest_limit`.
    """

    minimum: float
    delta: float = GUEST_DELTA

    @property
    def name(self) -> str:
        """Get the name of the rule."""
        return "Guest"

    @property
    def limit(self) -> str:
        """Get the limit, e.g. `>= 80`."""
        return f">= {self.minimum:g}"

    def evaluate(self, pairs: Collection[bool | None]) -> AcceptanceResult:
        """Evaluate the percentage of the evaluable pairs within the limits.

        Args:
            pairs: whether every pair is within, `None` for one which is not
                evaluable.

        Returns:
            The result, not evaluable without an evaluable pair.
        """
        evaluable = [within for within in pairs if within is not None]
        value = 100.0 * sum(evaluable) / len(evaluable) if evaluable else math.nan
        return _result(self.name, value, self.limit, value >= self.minimum)

    def problems(self) -> list[str]:
        """Get the problems of the rule: a minimum outside 0 to 100, a delta not above 1."""
        found: list[str] = []
        if not 0 <= self.minimum <= 100:
            found.append(
                f"the minimum of the ratio acceptance must be between 0 and 100, "
                f"not {self.minimum:g}"
            )
        if not self.delta > 1:
            found.append(
                f"the delta of the ratio acceptance must be above 1, not {self.delta:g}"
            )
        return found


@dataclass(frozen=True, kw_only=True)
class ContextOfUse:
    """What a model is used for: the question, the domain, the evidence and the acceptance.

    Attributes:
        id: the slug of the context, unique in the model.
        title: the title.
        question: the question of the context, one sentence.
        level: the level on the ladder.
        influence: the influence of the model on the decision.
        consequence: the consequence of a wrong decision.
        quantities: the quantities asked for, at least one which is no ratio.
        domain: the domain claimed, at least one value of the metadata above
            C1.
        evidence: the evidence.
        acceptance: the rules, empty for the default of the risk.
        limitations: the evidence of what the model does not describe.
        rationale: why this level, influence and consequence.
        extrapolation: the argument for a domain without validation evidence.
        min_points: the evaluable decisive points the verdict needs, at least
            1; `None` for the default of the risk, see `minimum_points`.
        ratio_acceptance: the rule of the pairs of the assessed ratio
            quantities, `None` for the default of the risk, see `ratio_rule`.
    """

    id: str
    title: str
    question: str
    level: CouLevel
    influence: Level
    consequence: Level
    quantities: tuple[QuantityOfInterest, ...]
    domain: Domain
    evidence: tuple[Evidence, ...]
    acceptance: tuple[AcceptanceRule, ...] = ()
    limitations: tuple[Evidence, ...] = ()
    rationale: str = ""
    extrapolation: str = ""
    min_points: int | None = None
    ratio_acceptance: GuestAcceptance | None = None

    @property
    def risk(self) -> Level:
        """Get the model risk of the influence and the consequence."""
        from pkpdsim.credibility.risk import model_risk

        return model_risk(self.influence, self.consequence)

    @property
    def rules(self) -> tuple[AcceptanceRule, ...]:
        """Get the acceptance rules: the declared ones or the default of the risk."""
        from pkpdsim.credibility.risk import default_acceptance

        return self.acceptance or default_acceptance(self.risk)

    @property
    def ratio_rule(self) -> GuestAcceptance:
        """Get the rule of the ratio pairs: the declared or the default of the risk."""
        from pkpdsim.credibility.risk import DEFAULT_RATIO_ACCEPTANCE

        return self.ratio_acceptance or DEFAULT_RATIO_ACCEPTANCE[self.risk]

    @property
    def ratios(self) -> tuple[QuantityOfInterest, ...]:
        """Get the assessed ratio quantities, those with a `Ratio`."""
        return tuple(q for q in self.quantities if q.pairing is not None)

    @property
    def minimum_points(self) -> int:
        """Get the evaluable points the verdict needs: the declared or the default of the risk."""
        from pkpdsim.credibility.risk import DEFAULT_MIN_POINTS

        if self.min_points is not None:
            return self.min_points
        return DEFAULT_MIN_POINTS[self.risk]

    @property
    def decisive_role(self) -> EvidenceRole:
        """Get the role of the evidence which decides, see `CouLevel.decisive_role`."""
        return self.level.decisive_role

    def problems(self) -> list[str]:
        """Get the problems of the context of use, each naming the context."""
        who = f"context of use '{self.id}'"
        found: list[str] = []
        if not ID_PATTERN.match(self.id):
            found.append(f"{who}: the id must match {ID_PATTERN.pattern}")
        if not self.question.rstrip().endswith("?"):
            found.append(f"{who}: the question must be a question, ending with '?'")
        if not any(not q.ratio for q in self.quantities):
            # the absolute quantities decide the verdict, a ratio adds to it
            found.append(
                f"{who}: at least one quantity of interest which is no ratio "
                "is required"
            )
        if self.level is not CouLevel.C1 and not any(self.domain.metadata.values()):
            found.append(
                f"{who}: a {self.level} context claims at least one value of the "
                "metadata in its domain"
            )
        roles = {e.role for e in self.evidence}
        if self.level is CouLevel.C1:
            if EvidenceRole.CALIBRATION not in roles:
                found.append(f"{who}: a C1 context needs calibration evidence")
        elif self.level is CouLevel.C4:
            if not self.extrapolation:
                found.append(f"{who}: a C4 context needs an extrapolation argument")
        elif EvidenceRole.VALIDATION not in roles and not self.extrapolation:
            found.append(
                f"{who}: a {self.level} context needs validation evidence "
                "or an extrapolation argument"
            )
        if self.min_points is not None and self.min_points < 1:
            found.append(
                f"{who}: the minimum of points must be at least 1, not "
                f"{self.min_points}"
            )
        unclaimed = [e for e in self.evidence if isinstance(e.selection, Unclaimed)]
        if len(unclaimed) > 1:
            found.append(
                f"{who}: at most one unclaimed evidence, {len(unclaimed)} given"
            )
        found.extend(
            f"{who}: the unclaimed evidence '{e.name}' must have the role "
            f"validation, not {e.role}"
            for e in unclaimed
            if e.role is not EvidenceRole.VALIDATION
        )
        found.extend(
            f"{who}: the limitation '{e.name}' must have the role limitation, "
            f"not {e.role}"
            for e in self.limitations
            if e.role is not EvidenceRole.LIMITATION
        )
        found.extend(f"{who}: {problem}" for problem in self.domain.problems())
        if self.ratio_acceptance is not None:
            found.extend(f"{who}: {p}" for p in self.ratio_acceptance.problems())
        found.extend(
            f"{who}: the ratio '{q.id}' needs the observable of its arms"
            for q in self.ratios
            if q.pairing is not None and not q.pairing.observable.strip()
        )
        return found


@dataclass(frozen=True, kw_only=True)
class Credibility:
    """The contexts of use of a model.

    Attributes:
        contexts: the contexts of use.
        limitations: the model limitations, evidence of what it does not describe.
        notes: notes on the declaration.
        covariates: the covariates of the subjects whose simulated values the
            equivalency of the inputs checks against the reported ones, the
            target of the model and the measurement type of PK-DB of its
            reported value (`DEFAULT_COVARIATES`, the body weight and the
            height), see `pkpdsim.credibility.inputs`; a model whose targets
            have other ids or which simulates further covariates declares
            them.

    Raises:
        DeclarationError: with every problem of the declaration.
    """

    contexts: tuple[ContextOfUse, ...]
    limitations: tuple[Evidence, ...] = ()
    notes: str = ""
    covariates: Mapping[str, str] = field(default=DEFAULT_COVARIATES, hash=False)

    def __post_init__(self) -> None:
        """Check the declaration.

        Raises:
            DeclarationError: with every problem of the declaration.
        """
        found: list[str] = []
        seen: set[str] = set()
        for cou in self.contexts:
            if cou.id in seen:
                found.append(f"context of use '{cou.id}': the id is not unique")
            seen.add(cou.id)
            found.extend(cou.problems())
        found.extend(
            f"model limitation '{e.name}' must have the role limitation, not {e.role}"
            for e in self.limitations
            if e.role is not EvidenceRole.LIMITATION
        )
        found.extend(
            f"covariate '{target}': the target and the measurement type must be "
            "texts which are not empty"
            for target, measurement in self.covariates.items()
            if not (
                isinstance(target, str)
                and target.strip()
                and isinstance(measurement, str)
                and measurement.strip()
            )
        )
        if found:
            raise DeclarationError(found)
        object.__setattr__(self, "covariates", MappingProxyType(dict(self.covariates)))

    def context(self, id: str) -> ContextOfUse:
        """Get a context of use by its id.

        Args:
            id: the id of the context.

        Returns:
            The context.

        Raises:
            KeyError: if there is no such context.
        """
        for cou in self.contexts:
            if cou.id == id:
                return cou
        raise KeyError(
            f"No context of use '{id}'; known: {[c.id for c in self.contexts]}"
        )
