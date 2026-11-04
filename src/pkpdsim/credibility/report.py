"""The outputs of a credibility assessment: the record, the tables and the report.

`write_outputs` writes an assessment into a directory of its own:

- `credibility.json`: the assessment, see `CredibilityAssessment.to_json`.
- `contexts.tsv`: one row per context of use with its level, influence,
  consequence, risk, verdict, the decisive and the evaluable points, the
  minimum of evaluable points, the uncovered quantities, whether it is
  verified and per acceptance rule its value (`<name>`) and limit
  (`<name>_limit`).
- `evidence.tsv`: one row per context of use and evidence and per limitation
  (the limitations of the model with the context `""`): the role, the points
  (`n`), mappings and studies, GMFE, AFE and the criteria.
- `domain.tsv`: one row per context of use, field and value of the domain:
  whether it is claimed, its coverage and its points (`n`).
- `quantities.tsv`: one row per context of use and quantity of interest: its
  description, whether it is a ratio, its coverage and its points.
- `ratios.tsv`: one row per context of use and pair of the arms of a ratio
  quantity: the fit, the study, the arms, the observed and the predicted
  ratio, the limit of Guest et al. and whether the pair is within.
- `verification.tsv`: one row per context of use and check of the
  verification, empty if the assessment was not verified: the factor, the
  compared points (`n`), the largest relative difference and its point,
  whether the acceptance changed, whether the check passed and its problems;
  the checks of the battery are rows of every context of use.
- `temporal.tsv`: one row per context of use and study the model collected
  after the fit of its decisive points: the fit, the day of the fit, of the
  collection and of the release in PK-DB, the points and the evaluable points.
- `uncertainty.tsv`: one row per context of use and acceptance rule over the
  parameter draws of the adopted fits, empty without the propagation: the
  probability that the rule is met with its anytime-valid sequence, the
  decision of the statement at the threshold, the draws and the failed and
  needed draws.
- `uncertainty_quantities.tsv`: one row per context of use and quantity of
  interest which is no ratio: the points, the median fold width of the
  prediction bands and their coverage of the observations.
- `uncertainty_sensitivity.tsv`: one row per context of use, target (a
  quantity or a rule) and fitted parameter: the first-order index of the
  target over the parameter draws.
- `inputs.tsv`: one row per mapping of the adopted fits and simulated input
  (a covariate of the subjects or a dose): the simulated value and where it
  comes from, the reported value and range it is compared with, by which
  groups, individuals or interventions, the status and the ratio.
- `input_ranges.tsv`: one row per context of use and input of its decisive
  mappings: the mappings, the simulated range and the counts of the statuses.
- `studies.tsv`: one row per study with points in the adopted fits, the
  provenance of its data: the folder, the fits and points, the record of its
  collection (the days of the first and the latest collection and how the
  first is known, the release in PK-DB, the source digest, the commit of
  pkdb_data) and the tables which differ from the record.
- `cross_validation.tsv`: one row per cross-validated adopted fit and held-out
  study: its points, the GMFE and the share within 2-fold of the held-out and
  of the fitted predictions, the cost of the refit.
- `battery.tsv`: one row per check of the battery of the builds, empty if it
  did not run: the factor, the status, the cases, the largest deviation and
  where, the statement and the problems.
- `fates.tsv`: one row per build, dose and part of the fate of the dose of the
  battery: the substance, the dose (mmole), the part (mmole) and its share of
  the dose.
- `index.html`: the report (`render_html`), which extends `report_base.html`
  of `sbmlsim.report.templates` and loads nothing from the network.

`print_summary` prints the table of the contexts of use, their verdicts, the
role of the evidence which decides them and the evaluable points with their
minimum, and how many of them are verified.
"""

from __future__ import annotations

import logging
import math
import shutil
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Final

import jinja2
import pandas as pd
from rich import box
from rich.markup import escape
from rich.table import Table
from sbmlsim.console import console
from sbmlsim.report.templates import template_environment

from pkpdsim.credibility.assessment import (
    ContextResult,
    CredibilityAssessment,
    EvidenceResult,
    Verdict,
)
from pkpdsim.credibility.inputs import INPUT_COLUMNS
from pkpdsim.credibility.ratios import RATIO_COLUMNS
from pkpdsim.credibility.risk import DEFAULT_MIN_POINTS
from pkpdsim.credibility.studies import STUDY_COLUMNS
from pkpdsim.credibility.temporal import TEMPORAL_COLUMNS
from pkpdsim.credibility.verification import NUMERICAL_CHECKS, Check

__all__ = [
    "BATTERY_COLUMNS",
    "CONTEXT_COLUMNS",
    "CROSS_VALIDATION_COLUMNS",
    "DOMAIN_COLUMNS",
    "EVIDENCE_COLUMNS",
    "FATE_COLUMNS",
    "INPUT_RANGE_COLUMNS",
    "QUANTITY_COLUMNS",
    "RECORD",
    "TEMPLATE",
    "TEMPLATE_DIR",
    "VERIFICATION_COLUMNS",
    "battery_table",
    "contexts_table",
    "cross_validation_table",
    "domain_table",
    "evidence_table",
    "fates_table",
    "input_ranges_table",
    "inputs_table",
    "print_summary",
    "quantities_table",
    "ratios_table",
    "render_html",
    "sensitivity_table",
    "studies_table",
    "temporal_table",
    "uncertainty_quantities_table",
    "uncertainty_table",
    "verification_table",
    "write_outputs",
]

logger = logging.getLogger(__name__)

#: the directory of the templates of pkpdsim
TEMPLATE_DIR: Final = Path(__file__).resolve().parents[1] / "resources" / "templates"
#: the template of the report
TEMPLATE: Final = "credibility.html"
#: the record of the assessment, which marks a directory of the outputs
RECORD: Final = "credibility.json"

#: the style of a verdict in the console
_VERDICT_STYLES: Final = {
    Verdict.SUPPORTED: "green",
    Verdict.SUPPORTED_WITH_EXTRAPOLATION: "blue",
    Verdict.NOT_SUPPORTED: "red",
    Verdict.INSUFFICIENT_EVIDENCE: "orange3",
    Verdict.NOT_VERIFIED: "magenta",
}

#: the leading columns of `contexts_table`
CONTEXT_COLUMNS: Final = (
    "id",
    "title",
    "level",
    "influence",
    "consequence",
    "risk",
    "verdict",
    "n_decisive",
    "n_evaluable",
    "minimum_points",
    "uncovered_quantities",
    "verified",
)
#: the leading columns of `evidence_table`, the metrics follow
EVIDENCE_COLUMNS: Final = (
    "context",
    "name",
    "role",
    "selection",
    "n",
    "n_mappings",
    "n_studies",
    "outliers_dropped",
    "limitation_subtracted",
)
#: the columns of `domain_table`
DOMAIN_COLUMNS: Final = ("context", "field", "value", "claimed", "coverage", "n")
#: the columns of `quantities_table`
QUANTITY_COLUMNS: Final = (
    "context",
    "quantity",
    "description",
    "ratio",
    "coverage",
    "n_points",
)
#: the columns of `verification_table`
VERIFICATION_COLUMNS: Final = (
    "context",
    "check",
    "factor",
    "n",
    "max_difference",
    "at",
    "acceptance_changed",
    "passed",
    "problems",
)


#: the columns of `cross_validation_table`
CROSS_VALIDATION_COLUMNS: Final = (
    "fit",
    "subset",
    "study",
    "n_points",
    "gmfe_held_out",
    "gmfe_fitted",
    "within_2fold_held_out",
    "within_2fold_fitted",
    "cost",
)
#: the columns of `battery_table`
BATTERY_COLUMNS: Final = (
    "check",
    "factor",
    "status",
    "cases",
    "largest",
    "at",
    "statement",
    "problems",
)
#: the columns of `fates_table`
FATE_COLUMNS: Final = (
    "build",
    "dose",
    "substance",
    "amount",
    "part",
    "value",
    "share",
)


def contexts_table(assessment: CredibilityAssessment) -> pd.DataFrame:
    """Get one row per context of use with the values and limits of its rules.

    Args:
        assessment: the assessment.

    Returns:
        The columns `id`, `title`, `level`, `influence`, `consequence`,
        `risk`, `verdict`, `n_decisive`, `n_evaluable`, `minimum_points`,
        `uncovered_quantities` (the ids, comma separated, empty without),
        `verified` (empty if the assessment was not verified) and per
        acceptance rule, the rule of the ratio pairs (`Guest`) included,
        `<name>` and `<name>_limit`, empty for a context without the rule.
    """
    rows: list[dict[str, object]] = []
    for c in assessment.contexts:
        row: dict[str, object] = {
            "id": c.id,
            "title": c.title,
            "level": c.level.value,
            "influence": c.influence.value,
            "consequence": c.consequence.value,
            "risk": c.risk.value,
            "verdict": c.verdict.value,
            "n_decisive": c.n_decisive,
            "n_evaluable": c.n_evaluable,
            "minimum_points": c.minimum_points,
            "uncovered_quantities": ",".join(c.uncovered_quantities),
            "verified": c.verified,
        }
        for result in [
            *c.acceptance,
            *([] if c.ratio_acceptance is None else [c.ratio_acceptance]),
        ]:
            row[result["name"]] = result["value"]
            row[f"{result['name']}_limit"] = result["limit"]
        rows.append(row)
    return pd.DataFrame(rows, columns=_columns(rows, CONTEXT_COLUMNS))


def _columns(rows: Iterable[Mapping[str, object]], leading: Sequence[str]) -> list[str]:
    """Get the leading columns and the other keys of rows in their first order."""
    columns = list(leading)
    for row in rows:
        columns += [key for key in row if key not in columns]
    return columns


def evidence_table(assessment: CredibilityAssessment) -> pd.DataFrame:
    """Get one row per evidence of a context of use and per limitation.

    Args:
        assessment: the assessment.

    Returns:
        The columns `context` (`""` for a limitation of the model), `name`,
        `role`, `selection`, `n` (the points), `n_mappings`, `n_studies`,
        `outliers_dropped`, `limitation_subtracted` and the metrics: GMFE, AFE,
        the criteria and the values of the acceptance rules of the context,
        empty if not evaluable.
    """
    rows = [
        {
            "context": e.context,
            "name": e.name,
            "role": e.role.value,
            "selection": e.selection,
            "n": e.n_points,
            "n_mappings": e.n_mappings,
            "n_studies": e.n_studies,
            "outliers_dropped": e.outliers_dropped,
            "limitation_subtracted": e.limitation_subtracted,
            **e.metrics,
        }
        for e in (*assessment.evidence, *assessment.limitations)
    ]
    return pd.DataFrame(rows, columns=_columns(rows, EVIDENCE_COLUMNS))


def domain_table(assessment: CredibilityAssessment) -> pd.DataFrame:
    """Get one row per context of use, field and value of its domain.

    Args:
        assessment: the assessment.

    Returns:
        The columns `context`, `field`, `value`, `claimed`, `coverage` (empty
        for a value outside the claim) and `n`, the points, see
        `CoverageEntry`.
    """
    rows = [
        {
            "context": e.context,
            "field": e.field,
            "value": e.value,
            "claimed": e.claimed,
            "coverage": None if e.coverage is None else e.coverage.value,
            "n": e.n_points,
        }
        for e in assessment.coverage
    ]
    return pd.DataFrame(rows, columns=list(DOMAIN_COLUMNS))


def quantities_table(assessment: CredibilityAssessment) -> pd.DataFrame:
    """Get one row per context of use and quantity of interest.

    Args:
        assessment: the assessment.

    Returns:
        The columns `context`, `quantity`, `description`, `ratio`, `coverage`
        (empty for a ratio) and `n_points`, see `QuantityRecord`.
    """
    rows = [
        {
            "context": c.id,
            "quantity": q["id"],
            "description": q["description"],
            "ratio": q["ratio"],
            "coverage": q["coverage"],
            "n_points": q["n_points"],
        }
        for c in assessment.contexts
        for q in c.quantities
    ]
    return pd.DataFrame(rows, columns=list(QUANTITY_COLUMNS))


def verification_table(assessment: CredibilityAssessment) -> pd.DataFrame:
    """Get one row per context of use and check of the verification.

    Args:
        assessment: the assessment.

    Returns:
        The columns of `VERIFICATION_COLUMNS`, see `CheckRecord`, the problems
        joined by `; `; no rows if the assessment was not verified.
    """
    rows = [
        {
            "context": c.id,
            "check": r["check"],
            "factor": r["factor"],
            "n": r["n_points"],
            "max_difference": r["max_difference"],
            "at": r["at"],
            "acceptance_changed": r["acceptance_changed"],
            "passed": r["passed"],
            "problems": "; ".join(r["problems"]),
        }
        for c in assessment.contexts
        for r in c.verification
    ]
    return pd.DataFrame(rows, columns=list(VERIFICATION_COLUMNS))


def ratios_table(assessment: CredibilityAssessment) -> pd.DataFrame:
    """Get one row per context of use and pair of the arms of a ratio quantity.

    Args:
        assessment: the assessment.

    Returns:
        The column `context` and the columns of `RATIO_COLUMNS`, see
        `RatioRecord`; no rows without an assessed ratio.
    """
    rows = [{"context": c.id, **r} for c in assessment.contexts for r in c.ratios]
    return pd.DataFrame(rows, columns=["context", *RATIO_COLUMNS])


def temporal_table(assessment: CredibilityAssessment) -> pd.DataFrame:
    """Get one row per context of use and study collected after its fit.

    Args:
        assessment: the assessment.

    Returns:
        The column `context` and the columns of `TEMPORAL_COLUMNS`, see
        `TemporalStudy`; no rows without a new study.
    """
    rows = [
        {"context": c.id, **study}
        for c in assessment.contexts
        if c.temporal is not None
        for study in c.temporal["studies"]
    ]
    return pd.DataFrame(rows, columns=["context", *TEMPORAL_COLUMNS])


#: the columns of `uncertainty.tsv`
UNCERTAINTY_COLUMNS = (
    "context",
    "rule",
    "limit",
    "estimate",
    "low",
    "high",
    "threshold",
    "decision",
    "n",
    "n_failed",
    "n_needed",
)
#: the columns of `uncertainty_quantities.tsv`
UNCERTAINTY_QUANTITY_COLUMNS = (
    "context",
    "quantity",
    "n_points",
    "band_fold",
    "coverage",
)


def uncertainty_table(assessment: CredibilityAssessment) -> pd.DataFrame:
    """Get one row per context of use and rule over the parameter draws.

    Args:
        assessment: the assessment.

    Returns:
        The columns of `UNCERTAINTY_COLUMNS`, see `RuleProbability`; no rows
        without the propagation.
    """
    rows = [
        {
            "context": c.id,
            "rule": r["name"],
            "threshold": c.uncertainty["threshold"],
            **{k: v for k, v in r.items() if k != "name"},
        }
        for c in assessment.contexts
        if c.uncertainty is not None
        for r in c.uncertainty["rules"]
    ]
    return pd.DataFrame(rows, columns=list(UNCERTAINTY_COLUMNS))


def uncertainty_quantities_table(assessment: CredibilityAssessment) -> pd.DataFrame:
    """Get one row per context of use and quantity with its prediction band.

    Args:
        assessment: the assessment.

    Returns:
        The columns of `UNCERTAINTY_QUANTITY_COLUMNS`, see `QuantityBand`; no
        rows without the propagation.
    """
    rows = [
        {"context": c.id, **q}
        for c in assessment.contexts
        if c.uncertainty is not None
        for q in c.uncertainty["quantities"]
    ]
    return pd.DataFrame(rows, columns=list(UNCERTAINTY_QUANTITY_COLUMNS))


#: the columns of `uncertainty_sensitivity.tsv`
SENSITIVITY_COLUMNS = ("context", "target", "kind", "fit", "parameter", "index")


def sensitivity_table(assessment: CredibilityAssessment) -> pd.DataFrame:
    """Get one row per context of use, target and parameter of the draws.

    Args:
        assessment: the assessment.

    Returns:
        The columns of `SENSITIVITY_COLUMNS`, see `SensitivityIndex`; no rows
        without the propagation.
    """
    rows = [
        {"context": c.id, **i}
        for c in assessment.contexts
        if c.uncertainty is not None
        for i in c.uncertainty["sensitivity"]
    ]
    return pd.DataFrame(rows, columns=list(SENSITIVITY_COLUMNS))


def inputs_table(assessment: CredibilityAssessment) -> pd.DataFrame:
    """Get one row per mapping of the adopted fits and simulated input.

    Args:
        assessment: the assessment.

    Returns:
        The columns of `INPUT_COLUMNS`, see `InputCheck`; the reporting
        groups, individuals and interventions joined by `;`.
    """
    rows = [
        {**check, "reported_by": ";".join(check["reported_by"])}
        for check in assessment.inputs
    ]
    return pd.DataFrame(rows, columns=list(INPUT_COLUMNS))


#: the columns of `input_ranges.tsv`
INPUT_RANGE_COLUMNS = (
    "context",
    "target",
    "kind",
    "unit",
    "n_mappings",
    "low",
    "high",
    "n_equivalent",
    "n_in_range",
    "n_differs",
    "n_not_reported",
)


def input_ranges_table(assessment: CredibilityAssessment) -> pd.DataFrame:
    """Get one row per context of use and input of its decisive mappings.

    Args:
        assessment: the assessment.

    Returns:
        The columns of `INPUT_RANGE_COLUMNS`, see `InputSummary`; no rows
        without the checks of the inputs.
    """
    rows = [
        {"context": c.id, **summary}
        for c in assessment.contexts
        if c.inputs is not None
        for summary in c.inputs["inputs"]
    ]
    return pd.DataFrame(rows, columns=list(INPUT_RANGE_COLUMNS))


def studies_table(assessment: CredibilityAssessment) -> pd.DataFrame:
    """Get one row per study with points in the adopted fits.

    Args:
        assessment: the assessment.

    Returns:
        The columns of `STUDY_COLUMNS`, see `StudyData`; the fits and the
        changed tables joined by `;`.
    """
    rows = [
        {
            **study,
            "fits": ";".join(study["fits"]),
            "changed": (
                None if study["changed"] is None else ";".join(study["changed"])
            ),
        }
        for study in assessment.studies
    ]
    return pd.DataFrame(rows, columns=list(STUDY_COLUMNS))


def cross_validation_table(assessment: CredibilityAssessment) -> pd.DataFrame:
    """Get one row per cross-validated fit and held-out study.

    Args:
        assessment: the assessment.

    Returns:
        The columns `fit`, `subset` and the fields of `HeldOutStudy`; no rows
        without a cross-validation.
    """
    rows = [
        {"fit": cv.fit, "subset": cv.subset, **study.to_dict()}
        for cv in assessment.cross_validations
        for study in cv.studies
    ]
    return pd.DataFrame(rows, columns=list(CROSS_VALIDATION_COLUMNS))


def battery_table(assessment: CredibilityAssessment) -> pd.DataFrame:
    """Get one row per check of the battery of the builds of the model.

    Args:
        assessment: the assessment.

    Returns:
        The columns of `BATTERY_COLUMNS`, see `BatteryResult`, the problems
        joined by `; `; no rows if the battery did not run.
    """
    battery = (
        None if assessment.verification is None else assessment.verification.battery
    )
    rows = [
        {
            "check": r.check.value,
            "factor": r.check.factor,
            "status": r.status.value,
            "cases": r.cases,
            "largest": r.largest,
            "at": r.at,
            "statement": r.statement,
            "problems": "; ".join(r.problems),
        }
        for r in ([] if battery is None else battery.results)
    ]
    return pd.DataFrame(rows, columns=list(BATTERY_COLUMNS))


def fates_table(assessment: CredibilityAssessment) -> pd.DataFrame:
    """Get one row per build, dose and part of the fates of the doses of the battery.

    Args:
        assessment: the assessment.

    Returns:
        The columns of `FATE_COLUMNS`, see `Fate`: the dose and the part in
        mmole and the share of the part in the dose; no rows without.
    """
    battery = (
        None if assessment.verification is None else assessment.verification.battery
    )
    rows = [
        {
            "build": f.build,
            "dose": f.dose,
            "substance": f.substance,
            "amount": f.amount,
            "part": part,
            "value": value,
            "share": value / f.amount,
        }
        for f in ([] if battery is None else battery.fates)
        for part, value in f.parts.items()
    ]
    return pd.DataFrame(rows, columns=list(FATE_COLUMNS))


def _number(value: float | None) -> str:
    """Show a metric with three significant digits, `n/a` if not evaluable."""
    if value is None or math.isnan(value):
        return "n/a"
    return f"{value:.3g}" if abs(value) < 1000 else f"{value:.0f}"


def _metric_names(results: Sequence[EvidenceResult]) -> list[str]:
    """Get the names of the metrics of evidence in their first order."""
    names: list[str] = []
    for result in results:
        names += [name for name in result.metrics if name not in names]
    return names


def _context_view(
    context: ContextResult, assessment: CredibilityAssessment
) -> dict[str, object]:
    """Get what the report shows of a context of use."""
    evidence = [e for e in assessment.evidence if e.context == context.id]
    return {
        "result": context,
        "evidence": evidence,
        "metrics": _metric_names(evidence),
        "coverage": [c for c in assessment.coverage if c.context == context.id],
    }


def render_html(assessment: CredibilityAssessment) -> str:
    """Render the report of an assessment.

    The template `credibility.html` of pkpdsim extends `report_base.html` of
    sbmlsim: the environment of `sbmlsim.report.templates` with the templates
    of pkpdsim first. The sections: the overview (the contexts of use and the
    limitations of the model), one section per context of use, the
    limitations, the unclaimed validation data, the verification and the
    provenance.

    Args:
        assessment: the assessment.

    Returns:
        The HTML of the report.
    """
    env = template_environment()
    if env.loader is None:
        raise RuntimeError("The environment of sbmlsim has no loader.")
    env.loader = jinja2.ChoiceLoader(
        [jinja2.FileSystemLoader(str(TEMPLATE_DIR)), env.loader]
    )
    env.filters["number"] = _number
    verdicts = [v for v in Verdict if any(c.verdict is v for c in assessment.contexts)]
    return env.get_template(TEMPLATE).render(
        assessment=assessment,
        provenance=assessment.provenance,
        contexts=[_context_view(c, assessment) for c in assessment.contexts],
        verdicts=verdicts,
        verdict_counts={
            v: sum(c.verdict is v for c in assessment.contexts) for v in verdicts
        },
        model_limitations=[e for e in assessment.limitations if not e.context],
        default_min_points=[(risk.value, n) for risk, n in DEFAULT_MIN_POINTS.items()],
        limitation_metrics=_metric_names(assessment.limitations),
        verification=assessment.verification,
        check_labels={check.value: check.label for check in Check},
        prediction_checks=[c.value for c in (*NUMERICAL_CHECKS, Check.SETUP)],
    )


def _tsv(table: pd.DataFrame) -> str:
    """Get a table as TSV, missing values empty."""
    return table.to_csv(sep="\t", index=False, lineterminator="\n")


def write_outputs(assessment: CredibilityAssessment, output_dir: Path) -> Path:
    """Write the record, the tables and the report of an assessment.

    A directory of the outputs is replaced as a whole, so no file of an earlier
    assessment is left; a directory is one when it holds `credibility.json`.
    The outputs are rendered before anything is removed.

    Args:
        assessment: the assessment.
        output_dir: the directory, created if missing.

    Returns:
        The directory.

    Raises:
        FileExistsError: if `output_dir` is a file or a non-empty directory
            without `credibility.json`.
    """
    if output_dir.exists() and not output_dir.is_dir():
        raise FileExistsError(f"The output directory '{output_dir}' is a file.")
    if (
        output_dir.is_dir()
        and any(output_dir.iterdir())
        and not (output_dir / RECORD).is_file()
    ):
        raise FileExistsError(
            f"The output directory '{output_dir}' is not empty and holds no "
            f"{RECORD}; it is not replaced."
        )
    files = {
        "contexts.tsv": _tsv(contexts_table(assessment)),
        "evidence.tsv": _tsv(evidence_table(assessment)),
        "domain.tsv": _tsv(domain_table(assessment)),
        "quantities.tsv": _tsv(quantities_table(assessment)),
        "ratios.tsv": _tsv(ratios_table(assessment)),
        "temporal.tsv": _tsv(temporal_table(assessment)),
        "inputs.tsv": _tsv(inputs_table(assessment)),
        "input_ranges.tsv": _tsv(input_ranges_table(assessment)),
        "studies.tsv": _tsv(studies_table(assessment)),
        "uncertainty.tsv": _tsv(uncertainty_table(assessment)),
        "uncertainty_quantities.tsv": _tsv(uncertainty_quantities_table(assessment)),
        "uncertainty_sensitivity.tsv": _tsv(sensitivity_table(assessment)),
        "cross_validation.tsv": _tsv(cross_validation_table(assessment)),
        "verification.tsv": _tsv(verification_table(assessment)),
        "battery.tsv": _tsv(battery_table(assessment)),
        "fates.tsv": _tsv(fates_table(assessment)),
        "index.html": render_html(assessment),
    }
    if output_dir.is_dir():
        logger.info("Replace the credibility outputs in '%s'", output_dir)
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True)
    assessment.to_json(output_dir / RECORD)
    for name, text in files.items():
        (output_dir / name).write_text(text, encoding="utf-8")
    return output_dir


def _verification_line(assessment: CredibilityAssessment) -> str:
    """Get the line of the summary which tells whether the contexts are verified."""
    if assessment.verification is None:
        return "The predictions are not verified (--no-verification)."
    verified = sum(c.verified is True for c in assessment.contexts)
    n = len(assessment.contexts)
    style = "green" if verified == n else "magenta"
    checks = ", ".join(check.value for check in (*NUMERICAL_CHECKS, Check.SETUP))
    line = f"[{style}]{verified} of {n} contexts of use verified[/] ({checks}"
    battery = assessment.verification.battery
    if battery is None:
        return f"{line}; the battery of the builds did not run)."
    applied = [r for r in battery.results if r.passed is not None]
    passed = sum(r.passed is True for r in applied)
    failed = ", ".join(r.check.value for r in applied if r.passed is False)
    return (
        f"{line}; battery of the builds: {passed} of {len(applied)} applicable "
        f"checks passed{f', failed: {failed}' if failed else ''})."
    )


def print_summary(assessment: CredibilityAssessment) -> None:
    """Print the contexts of use and their verdicts on the console of sbmlsim.

    Args:
        assessment: the assessment.
    """
    table = Table(
        title=f"Credibility of {escape(assessment.provenance.title)}",
        box=box.SIMPLE_HEAD,
        # one space between the columns and no edges, the table fits 79
        # columns, the width of a console of windows
        padding=0,
        pad_edge=False,
        show_edge=False,
    )
    # an id is not cut, a column which does not fit folds its text; a verdict
    # wraps between its words (`supported with extrapolation` beside a long id)
    table.add_column(
        "context",
        overflow="fold",
        min_width=max((len(c.id) for c in assessment.contexts), default=0),
    )
    table.add_column("level", no_wrap=True)
    table.add_column("risk", no_wrap=True)
    table.add_column(
        "verdict",
        min_width=max(len(word) for v in Verdict for word in v.label.split()),
    )
    # a description (C1) is decided on the data it was fitted to
    table.add_column("decided on", no_wrap=True)
    table.add_column("points", justify="right", no_wrap=True)
    # the evaluable points and their minimum
    table.add_column("evaluable", justify="right", no_wrap=True)
    for c in assessment.contexts:
        table.add_row(
            escape(c.id),
            c.level.value,
            c.risk.value,
            f"[{_VERDICT_STYLES[c.verdict]}]{c.verdict.label}[/]",
            c.decisive_role.value,
            str(c.n_decisive),
            # a record of pkpdsim 0.17 has no minimum
            f"{c.n_evaluable}/{c.minimum_points}"
            if c.minimum_points
            else str(c.n_evaluable),
        )
    # at most the 79 columns of a console of windows
    console.print(table, width=min(console.width, 79))
    console.print(_verification_line(assessment))
