"""The predictive check of the decisive points (`pkpdsim.credibility.variability`).

The record of a context, the outputs and the statement run on the hand-built
table of the test package with hand-built bands; `predictive_checks` and
`assess` run on its adopted fit with a variability of the fitted rate.
"""

import dataclasses
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from credibility_tables import fit_table
from pkpdsim_credibility.credibility import CREDIBILITY
from sbmlsim.fit.variability import Variability
from sbmlsim.simulation.sampling import LogNormal

from pkpdsim.cli import credibility_command
from pkpdsim.credibility import (
    CredibilityAssessment,
    Provenance,
    VariabilitySettings,
    assess,
    assess_table,
    evaluate_fits,
    predictive_checks,
)
from pkpdsim.credibility.report import (
    variability_quantities_table,
    variability_table,
    vpc_table,
    write_outputs,
)
from pkpdsim.credibility.scorecard import scorecard
from pkpdsim.credibility.variability import (
    PredictiveChecks,
    _wilson,
    variability_statement,
)
from pkpdsim.model import PKPDModel

PROVENANCE = Provenance(
    model="credibility",
    title="Toy model with every kind of evidence",
    fits=["DECAY_20261001_120000__abcde"],
    versions={"pkpdsim": "0.0.0"},
    commit=None,
    dirty=None,
    date="2026-11-04T12:00:00+01:00",
)

VARIABILITY = Variability({"k": LogNormal(cv=0.2)})


def _checks(width: float) -> PredictiveChecks:
    """Get bands of the points of the table, `width` around the observations.

    A width below 1 puts every observation outside its band.
    """
    table = fit_table().table
    dv = table["DV"].to_numpy(dtype=float)
    low, high = dv / (1 + width), dv * (1 + width)
    if width < 1:
        low, high = dv * 2.0, dv * 3.0
    points = pd.DataFrame(
        {
            "N": 4,
            "statistic": "mean",
            "PRED": table["IPRED"].to_numpy(dtype=float),
            "low": low,
            "median": dv,
            "high": high,
            "inside": [
                bool(lo <= v <= hi) for v, lo, hi in zip(dv, low, high, strict=True)
            ],
        },
        index=table.index,
    )
    return PredictiveChecks(
        settings=VariabilitySettings(subjects=16, groups=10),
        variability=VARIABILITY,
        points=points,
        vectors={"DECAY": 1},
        studies={"DECAY": 1},
        failed={"DECAY": 0},
    )


def test_the_settings() -> None:
    assert VariabilitySettings().subjects == 64
    with pytest.raises(ValueError, match="fewer than two subjects; no virtual group"):
        VariabilitySettings(subjects=1, groups=0)
    with pytest.raises(ValueError, match="not increasing"):
        VariabilitySettings(quantiles=(0.9, 0.1))


def test_the_wilson_interval() -> None:
    low, high = _wilson(9, 10, 0.95)
    assert low == pytest.approx(0.5958, abs=1e-4)
    assert high == pytest.approx(0.9821, abs=1e-4)
    assert _wilson(10, 10, 0.95)[1] == pytest.approx(1.0)


def test_the_record_of_a_context(tmp_path: Path) -> None:
    assessment = assess_table(
        CREDIBILITY, fit_table(), PROVENANCE, variability=_checks(1.0)
    )
    contexts = {c.id: c for c in assessment.contexts}
    dose = contexts["dose"]
    record = dose.variability
    assert record is not None
    assert record["n_points"] == dose.n_evaluable
    assert record["coverage"] == 1.0
    # every observation inside: wide unless the interval of few points holds 90 %
    low = record["low"]
    assert low is not None
    assert record["agreement"] == ("wide" if low > record["expected"] else "consistent")
    assert record["targets"] == {"k": "LogNormal(cv=0.2)"}
    assert record["n_groups"] == record["n_points"]
    assert [q["quantity"] for q in record["quantities"]] == ["[A]"]
    factor = next(f for f in dose.factors if f["id"] == "equivalency_of_inputs")
    assert factor["status"] == "assessed"
    assert factor["statement"] == variability_statement(record)
    assert "16 subjects at the adopted parameters" in factor["statement"]
    assert "hold 100 %" in factor["statement"]
    # every observation outside: the bands are narrow
    narrow = assess_table(
        CREDIBILITY, fit_table(), PROVENANCE, variability=_checks(0.5)
    )
    assert {c.variability["agreement"] for c in narrow.contexts if c.variability} == {
        "narrow"
    }
    path = tmp_path / "credibility.json"
    assessment.to_json(path)
    read = CredibilityAssessment.from_json(path)
    assert read.to_dict() == assessment.to_dict()
    assert len(read.vpc) == len(fit_table().table)
    # a record of pkpdsim 0.18 has no predictive check
    data = json.loads(path.read_text())
    del data["vpc"]
    for context in data["contexts"]:
        del context["variability"]
    old = CredibilityAssessment.from_dict(data)
    assert old.vpc == [] and old.contexts[0].variability is None
    table = variability_table(assessment)
    assert list(table["context"]) == [c.id for c in assessment.contexts]
    assert table.loc[0, "targets"] == "k LogNormal(cv=0.2)"
    assert set(variability_quantities_table(assessment)["quantity"]) == {"[A]"}
    assert len(vpc_table(assessment)) == len(fit_table().table)
    card = scorecard({"credibility": assessment}).set_index("context")
    assert card.loc["dose", "vpc_coverage"] == 1.0
    assert card.loc["dose", "vpc_agreement"] == record["agreement"]
    output = write_outputs(assessment, tmp_path / "out")
    for name in ("variability.tsv", "variability_quantities.tsv", "vpc.tsv"):
        assert (output / name).is_file()
    html = (output / "index.html").read_text()
    assert "<h4>Variability</h4>" in html
    assert "k LogNormal(cv=0.2)" in html


def test_the_statement_without_a_point() -> None:
    assessment = assess_table(
        CREDIBILITY, fit_table(), PROVENANCE, variability=_checks(1.0)
    )
    record = dict(assessment.contexts[0].variability or {})
    record.update(coverage=None, vectors=8, failed=2)
    statement = variability_statement(record)  # ty: ignore[invalid-argument-type]
    assert "nested in 8 draws of the parameters (2 failed)" in statement
    assert statement.endswith("No evaluable decisive point has a band.")


def test_the_predictive_checks_of_the_adopted_fit(
    credibility_model: PKPDModel,
) -> None:
    fits = evaluate_fits(credibility_model)
    settings = VariabilitySettings(subjects=8, groups=50, seed=2)
    checks = predictive_checks(credibility_model, fits, VARIABILITY, settings)
    points = checks.points
    assert list(points.index) == list(fits.table.index)
    assert set(points["statistic"]) == {"mean"}
    assert (points["N"] == 4).all()
    np.testing.assert_allclose(
        points["PRED"], fits.table["IPRED"].to_numpy(dtype=float), rtol=1e-6
    )
    banded = points.dropna(subset=["low"])
    assert (banded["low"] <= banded["high"]).all()
    assert checks.vectors == {"DECAY": 1}


def test_the_assessment_with_the_variability(
    credibility_model: PKPDModel, tmp_path: Path
) -> None:
    settings = VariabilitySettings(subjects=8, groups=50, seed=2)
    with pytest.raises(ValueError, match="declares no variability"):
        assess(credibility_model, verification=False, variability=settings)
    model = dataclasses.replace(
        credibility_model,
        credibility="pkpdsim_credibility.credibility:CREDIBILITY_VARIABILITY",
    )
    assessment = assess(model, verification=False, variability=settings)
    for context in assessment.contexts:
        record = context.variability
        assert record is not None
        assert record["subjects"] == 8
        assert record["n_points"] <= context.n_evaluable
    assert len(assessment.vpc) == len(evaluate_fits(model).table)
    for args in (
        ["--subjects", "8"],
        ["--variability", "--subjects", "1"],
        ["--variability", "--from-json", str(tmp_path / "credibility.json")],
    ):
        with pytest.raises(SystemExit) as exc:
            credibility_command(dataclasses.replace(model, results_path=tmp_path), args)
        assert exc.value.code == 2


def test_the_variability_between_the_studies(credibility_model: PKPDModel) -> None:
    fits = evaluate_fits(credibility_model)
    variability = Variability(
        {"k": LogNormal(cv=0.1)}, between=Variability({"k": LogNormal(cv=0.2)})
    )
    settings = VariabilitySettings(subjects=4, studies=3, groups=20, seed=2)
    checks = predictive_checks(credibility_model, fits, variability, settings)
    assert checks.studies == {"DECAY": 3}
    assessment = assess_table(CREDIBILITY, fits, PROVENANCE, variability=checks)
    record = assessment.contexts[0].variability
    assert record is not None
    assert (record["studies"], record["between"]) == (3, {"k": "LogNormal(cv=0.2)"})
    statement = variability_statement(record)
    assert (
        "The variability of the subjects (k LogNormal(cv=0.1)) and between the "
        "studies (k LogNormal(cv=0.2)): 3 studies of 4 subjects at the adopted "
        "parameters."
    ) in statement
    assert variability_table(assessment).loc[0, "between"] == "k LogNormal(cv=0.2)"
