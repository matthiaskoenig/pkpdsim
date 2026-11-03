"""The outputs of a credibility assessment: the record, the tables and the report.

`write_outputs` writes an assessment into a directory of its own:

- `credibility.json`: the assessment, see `CredibilityAssessment.to_json`.
- `contexts.tsv`: one row per context of use with its level, influence,
  consequence, risk, verdict, the decisive and the evaluable points and per
  acceptance rule its value (`<name>`) and limit (`<name>_limit`).
- `evidence.tsv`: one row per context of use and evidence and per limitation
  (the limitations of the model with the context `""`): the role, the points
  (`n`), mappings and studies, GMFE, AFE and the criteria.
- `domain.tsv`: one row per context of use, field and value of the domain:
  whether it is claimed, its coverage and its points (`n`).
- `index.html`: the report (`render_html`), which extends `report_base.html`
  of `sbmlsim.report.templates` and loads nothing from the network.

`print_summary` prints the table of the contexts of use and their verdicts.
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

__all__ = [
    "CONTEXT_COLUMNS",
    "DOMAIN_COLUMNS",
    "EVIDENCE_COLUMNS",
    "RECORD",
    "TEMPLATE",
    "TEMPLATE_DIR",
    "contexts_table",
    "domain_table",
    "evidence_table",
    "print_summary",
    "render_html",
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


def contexts_table(assessment: CredibilityAssessment) -> pd.DataFrame:
    """Get one row per context of use with the values and limits of its rules.

    Args:
        assessment: the assessment.

    Returns:
        The columns `id`, `title`, `level`, `influence`, `consequence`,
        `risk`, `verdict`, `n_decisive`, `n_evaluable` and per acceptance rule
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
        }
        for result in c.acceptance:
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


def _number(value: float) -> str:
    """Show a metric with three significant digits, `n/a` if not evaluable."""
    if math.isnan(value):
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
    limitations, the unclaimed validation data and the provenance.

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
        limitation_metrics=_metric_names(assessment.limitations),
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


def print_summary(assessment: CredibilityAssessment) -> None:
    """Print the contexts of use and their verdicts on the console of sbmlsim.

    Args:
        assessment: the assessment.
    """
    table = Table(
        title=f"Credibility of {escape(assessment.provenance.title)}",
        box=box.SIMPLE_HEAD,
    )
    # an id and the longest verdict are not cut, a column which does not fit
    # folds its text
    table.add_column(
        "context",
        overflow="fold",
        min_width=max((len(c.id) for c in assessment.contexts), default=0),
    )
    table.add_column("level", no_wrap=True)
    table.add_column("risk", no_wrap=True)
    table.add_column(
        "verdict",
        no_wrap=True,
        min_width=max(len(v.label) for v in Verdict),
    )
    table.add_column("points", justify="right", no_wrap=True)
    table.add_column("evaluable", justify="right", no_wrap=True)
    for c in assessment.contexts:
        table.add_row(
            c.id,
            c.level.value,
            c.risk.value,
            f"[{_VERDICT_STYLES[c.verdict]}]{c.verdict.label}[/]",
            str(c.n_decisive),
            str(c.n_evaluable),
        )
    console.print(table)
