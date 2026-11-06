"""The parameter uncertainty of the decisive predictions (`pkpdsim.credibility.uncertainty`).

The draws are tested on records built by hand, the record of a context on the
hand-built table of the test package with predictions per draw, and the
propagation on the adopted fit of the session fixture `credibility_model`.
"""

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from credibility_tables import fit_table
from pkpdsim_credibility.credibility import CREDIBILITY
from sbmlsim.fit.fisher import FisherInformation
from sbmlsim.fit.objects import ParameterScaleType
from sbmlsim.fit.records import ParameterUncertainty, UncertaintyMethod

from pkpdsim.cli import credibility_command
from pkpdsim.credibility import (
    ContextResult,
    CredibilityAssessment,
    Provenance,
    assess,
    assess_table,
)
from pkpdsim.credibility.report import write_outputs
from pkpdsim.credibility.scorecard import compare, scorecard
from pkpdsim.credibility.uncertainty import (
    FitDraws,
    Propagation,
    UncertaintySettings,
    parameter_draws,
    rule_probability,
    sensitivity_indices,
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
SETTINGS = UncertaintySettings(draws=256, seed=1)


def _parameter(pid: str, value: float) -> ParameterUncertainty:
    """Get a parameter on a log scale with a profile interval of 2-fold."""
    return ParameterUncertainty(
        pid=pid,
        value=value,
        unit="1/min",
        lower_bound=value / 100,
        upper_bound=value * 100,
        scale=ParameterScaleType.LOG10,
        alpha=0.95,
        z=1.96,
        method=UncertaintyMethod.PROFILE,
        ci_lower=value / 2,
        ci_upper=value * 2,
        fit_id="PK_20261001_120000__abcde",
        cost=1.0,
        record="PK_20261001_120000__abcde",
    )


def _fisher(correlation: float, matrix: np.ndarray | None = None) -> FisherInformation:
    covariance = np.array([[1.0, correlation], [correlation, 1.0]]) * 0.01
    return FisherInformation(
        opid="PK",
        sid="best",
        pids=["a", "b"],
        values=np.array([1.0, 2.0]),
        scale=ParameterScaleType.LOG10,
        matrix=np.linalg.inv(covariance) if matrix is None else matrix,
        cost=1.0,
        n=100,
        scales=[ParameterScaleType.LOG10, ParameterScaleType.LOG10],
    )


def test_the_draws_have_the_fisher_correlation() -> None:
    parameters = [_parameter("a", 1.0), _parameter("b", 2.0)]
    settings = UncertaintySettings(draws=1024, seed=1)
    draws = parameter_draws(
        "PK_20261001_120000__abcde", parameters, _fisher(0.8), settings
    )
    assert draws.correlated
    assert draws.pids == ("a", "b")
    assert draws.values.shape == (1024, 2)
    logs = np.log(draws.values)
    assert np.corrcoef(logs.T)[0, 1] == pytest.approx(0.8, abs=0.05)
    # the marginals are those of the profile, a median at the value
    assert np.median(draws.values, axis=0) == pytest.approx([1.0, 2.0], rel=0.05)
    assert np.std(logs[:, 0]) == pytest.approx(math.log(2) / 1.96, rel=0.1)
    # the same seed, the same draws
    again = parameter_draws(
        "PK_20261001_120000__abcde", parameters, _fisher(0.8), settings
    )
    assert np.array_equal(again.values, draws.values)


def test_independent_draws() -> None:
    parameters = [_parameter("a", 1.0), _parameter("b", 2.0)]
    singular = _fisher(0.0, matrix=np.array([[1.0, 1.0], [1.0, 1.0]]))
    draws = parameter_draws("PK_20261001_120000__abcde", parameters, singular, SETTINGS)
    assert not draws.correlated
    assert "rank 1 of 2" in draws.note
    none = parameter_draws("PK_20261001_120000__abcde", parameters, None, SETTINGS)
    assert "no Fisher information" in none.note


def test_the_settings() -> None:
    with pytest.raises(ValueError, match="fewer than two draws"):
        UncertaintySettings(draws=1)
    with pytest.raises(ValueError, match="not increasing"):
        UncertaintySettings(quantiles=(0.9, 0.1))
    with pytest.raises(ValueError, match="threshold"):
        UncertaintySettings(threshold=1.0)


def test_the_sensitivity_indices() -> None:
    parameters = [_parameter("a", 1.0), _parameter("b", 2.0)]
    settings = UncertaintySettings(draws=512, seed=2)
    draws = parameter_draws("PK_20261001_120000__abcde", parameters, None, settings)
    y = np.log(draws.values[:, 0])
    y[:3] = np.nan
    a, b = sensitivity_indices("Cve", "quantity", y, {"PK": draws})
    assert (a["parameter"], b["parameter"]) == ("a", "b")
    assert (a["index"] or 0.0) > 0.8
    assert (b["index"] or 0.0) < 0.05
    (constant, _) = sensitivity_indices("C", "quantity", np.ones(512), {"PK": draws})
    assert constant["index"] is None


def test_a_rule_probability() -> None:
    accepted = rule_probability("GMFE", "<= 1.5", [True] * 256, SETTINGS, 1)
    assert (accepted["estimate"], accepted["decision"]) == (1.0, "accepted")
    rejected = rule_probability("GMFE", "<= 1.5", [True, False] * 128, SETTINGS, 1)
    assert rejected["decision"] == "rejected"
    failed = rule_probability("GMFE", "<= 1.5", [True] * 31 + [None], SETTINGS, 1)
    assert (failed["n"], failed["n_failed"]) == (31, 1)
    nothing = rule_probability("GMFE", "<= 1.5", [None] * 4, SETTINGS, 1)
    assert nothing["decision"] == "not_evaluable"
    assert nothing["estimate"] is None


def _propagation(spread: float, failed: int = 0) -> Propagation:
    """Get predictions per draw, the reference times a factor of the spread."""
    table = fit_table().table
    rng = np.random.default_rng(1)
    factors = np.exp(rng.normal(0.0, spread, SETTINGS.draws))
    values = np.outer(table["IPRED"].to_numpy(dtype=float), factors)
    values[:, :failed] = np.nan
    draws = FitDraws(
        fit_id="DECAY_20261001_120000__abcde",
        pids=("k",),
        values=factors[:, None],
        units={"k": "1/min"},
        correlated=False,
        note="nothing to correlate",
    )
    return Propagation(
        settings=SETTINGS,
        draws={"DECAY": draws},
        predictions=pd.DataFrame(values, index=table.index),
        failed={"DECAY": failed},
    )


def _context(assessment: CredibilityAssessment, id: str) -> ContextResult:
    return next(c for c in assessment.contexts if c.id == id)


def _factor(context: ContextResult, id: str) -> dict[str, object]:
    return dict(next(f for f in context.factors if f["id"] == id))


def test_the_record_of_a_context() -> None:
    assessment = assess_table(
        CREDIBILITY, fit_table(), PROVENANCE, propagation=_propagation(0.01)
    )
    dose = _context(assessment, "dose")
    record = dose.uncertainty
    assert record is not None
    assert (record["n_draws"], record["seed"]) == (256, 1)
    assert [f["fit"] for f in record["fits"]] == ["DECAY"]
    assert [r["name"] for r in record["rules"]] == [a["name"] for a in dose.acceptance]
    assert all(r["decision"] == "accepted" for r in record["rules"])
    assert record["n_points"] == dose.n_evaluable
    (quantity,) = record["quantities"]
    assert quantity["band_fold"] == pytest.approx(math.exp(0.01 * 2 * 1.645), rel=0.02)
    factor = _factor(dose, "model_inputs")
    assert factor["status"] == "assessed"
    statement = str(factor["statement"])
    assert "propagated with 256 draws" in statement
    assert "carried most by" in statement
    # the one parameter carries the variance of the quantity and of the rules
    by_target = {i["target"]: i for i in record["sensitivity"]}
    assert set(by_target) == {
        quantity["quantity"],
        *(r["name"] for r in record["rules"]),
    }
    assert by_target[quantity["quantity"]]["kind"] == "quantity"
    assert (by_target[quantity["quantity"]]["index"] or 0.0) > 0.8
    assert "accepted" in statement
    # the verdicts do not change, the factor is not assessed without
    plain = assess_table(CREDIBILITY, fit_table(), PROVENANCE)
    assert [c.verdict for c in assessment.contexts] == [
        c.verdict for c in plain.contexts
    ]
    assert _factor(_context(plain, "dose"), "model_inputs")["status"] == "not_assessed"
    assert all(c.uncertainty is None for c in plain.contexts)


def test_a_wide_uncertainty_and_failed_draws() -> None:
    assessment = assess_table(
        CREDIBILITY, fit_table(), PROVENANCE, propagation=_propagation(0.5, failed=2)
    )
    record = _context(assessment, "dose").uncertainty
    assert record is not None
    assert any(r["decision"] == "rejected" for r in record["rules"])
    assert all(r["n_failed"] == 2 for r in record["rules"])
    assert record["fits"][0]["n_failed"] == 2
    fold = record["quantities"][0]["band_fold"]
    assert fold is not None
    assert fold > 3


def test_the_record_and_the_outputs(tmp_path: Path) -> None:
    assessment = assess_table(
        CREDIBILITY, fit_table(), PROVENANCE, propagation=_propagation(0.01)
    )
    path = tmp_path / "credibility.json"
    assessment.to_json(path)
    assert CredibilityAssessment.from_json(path).to_dict() == assessment.to_dict()
    # a record of pkpdsim 0.18 has no uncertainty
    data = json.loads(path.read_text(encoding="utf-8"))
    for context in data["contexts"]:
        del context["uncertainty"]
    assert all(
        c.uncertainty is None for c in CredibilityAssessment.from_dict(data).contexts
    )
    output = write_outputs(assessment, tmp_path / "out")
    rules = pd.read_csv(output / "uncertainty.tsv", sep="\t")
    assert set(rules["decision"]) == {"accepted"}
    quantities = pd.read_csv(output / "uncertainty_quantities.tsv", sep="\t")
    assert list(quantities.columns) == [
        "context",
        "quantity",
        "n_points",
        "band_fold",
        "coverage",
    ]
    sensitivity = pd.read_csv(output / "uncertainty_sensitivity.tsv", sep="\t")
    assert set(sensitivity["parameter"]) == {"k"}
    html = (output / "index.html").read_text(encoding="utf-8")
    assert "leading parameters" in html
    assert "Parameter uncertainty" in html
    assert "P(met)" in html


def test_the_scorecard() -> None:
    before = scorecard(
        {
            "toy": assess_table(
                CREDIBILITY, fit_table(), PROVENANCE, propagation=_propagation(0.01)
            )
        }
    )
    assert before.set_index("context").loc["dose", "uncertainty_robust"] == True  # noqa: E712
    after = scorecard(
        {
            "toy": assess_table(
                CREDIBILITY, fit_table(), PROVENANCE, propagation=_propagation(0.5)
            )
        }
    )
    regressions = [r for r in compare(before, after) if r.what == "uncertainty_robust"]
    assert "dose" in {r.context for r in regressions}
    old = before.drop(columns=["uncertainty_robust", "band_fold"])
    assert not [r for r in compare(old, after) if r.what == "uncertainty_robust"]


def test_the_propagation_of_the_adopted_fit(credibility_model: PKPDModel) -> None:
    assessment = assess(
        credibility_model,
        verification=False,
        uncertainty=UncertaintySettings(draws=8, seed=3),
    )
    for context in assessment.contexts:
        record = context.uncertainty
        assert record is not None
        (fit,) = record["fits"]
        assert (fit["n_parameters"], fit["n_failed"]) == (1, 0)
        assert record["n_points"] == context.n_evaluable
        assert all(r["n"] == 8 for r in record["rules"])
        assert all(
            q["band_fold"] is None or q["band_fold"] >= 1.0
            for q in record["quantities"]
        )


def test_the_options_of_the_command(
    credibility_model: PKPDModel, tmp_path: Path
) -> None:
    model = credibility_model.replace(results_path=tmp_path)
    for args in (["--draws", "8"], ["--uncertainty", "--draws", "1"]):
        with pytest.raises(SystemExit) as exc:
            credibility_command(model, args)
        assert exc.value.code == 2
