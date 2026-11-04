"""The scorecard and the qualification across models (`pkpdsim.credibility.scorecard`)."""

import dataclasses
from pathlib import Path

import pandas as pd
import pytest
from credibility_tables import fit_table
from test_credibility_ratios import PROVENANCE, _renal

from pkpdsim.credibility import CredibilityAssessment, Verdict, assess_table
from pkpdsim.credibility.scorecard import (
    SCORECARD_COLUMNS,
    Regression,
    compare,
    load_records,
    qualification,
    scorecard,
)
from pkpdsim.verify import verify


def _assessment(model: str, renal_factor: float = 1.1) -> CredibilityAssessment:
    """Get an assessment of the test package as the one of a model."""
    fits = fit_table(factors={("CredA", "fm_renal"): renal_factor}, times=True)
    provenance = dataclasses.replace(PROVENANCE, model=model)
    return assess_table(_renal(), fits, provenance)


def _write(directory: Path, *assessments: CredibilityAssessment) -> Path:
    for assessment in assessments:
        folder = directory / assessment.provenance.model / "credibility"
        folder.mkdir(parents=True)
        assessment.to_json(folder / "credibility.json")
    return directory


def test_load_records(tmp_path: Path) -> None:
    root = _write(tmp_path, _assessment("b"), _assessment("a"))
    records = load_records([root])
    assert list(records) == ["a", "b"]
    single = load_records([root / "a" / "credibility" / "credibility.json"])
    assert list(single) == ["a"]


def test_load_records_refuses_two_records_of_a_model(tmp_path: Path) -> None:
    _write(tmp_path / "one", _assessment("a"))
    _write(tmp_path / "two", _assessment("a"))
    with pytest.raises(ValueError, match="both of the model 'a'"):
        load_records([tmp_path])
    with pytest.raises(FileNotFoundError, match="No record"):
        load_records([])
    with pytest.raises(FileNotFoundError, match="does not exist"):
        load_records([tmp_path / "nothing"])


def test_the_scorecard() -> None:
    table = scorecard({"a": _assessment("a"), "b": _assessment("b", 2.0)})
    assert list(table.columns[: len(SCORECARD_COLUMNS)]) == list(SCORECARD_COLUMNS)
    assert len(table) == 6
    renal = table[table["context"] == "renal"].set_index("model")
    assert renal.loc["a", "verdict"] == Verdict.SUPPORTED.value
    assert renal.loc["b", "verdict"] == Verdict.NOT_SUPPORTED.value
    assert renal.loc["a", "Guest"] == 100.0
    assert renal.loc["b", "Guest"] == 0.0
    assert "within_2fold" in table.columns


def test_compare_finds_the_regressions() -> None:
    baseline = scorecard({"a": _assessment("a"), "b": _assessment("b")})
    current = scorecard({"a": _assessment("a", 2.0)})
    regressions = compare(baseline, current)
    assert (
        Regression(
            model="a",
            context="renal",
            what="verdict",
            baseline="supported",
            current="not_supported",
        )
        in regressions
    )
    assert (
        Regression(
            model="a", context="renal", what="Guest", baseline="100", current="0"
        )
        in regressions
    )
    missing = [r for r in regressions if r.what == "missing"]
    assert {(r.model, r.context) for r in missing} == {
        ("b", "describe"),
        ("b", "dose"),
        ("b", "renal"),
    }
    assert str(missing[0]) == "b/describe: missing present -> missing"
    assert compare(baseline, baseline) == []


@pytest.mark.parametrize(
    ("rule", "before", "after", "worse"),
    [
        ("within_2fold", 90.0, 86.0, False),
        ("within_2fold", 90.0, 84.0, True),
        ("GMFE", 1.3, 1.36, False),
        ("GMFE", 1.3, 1.4, True),
        ("AFE", 1.0, 1.04, False),
        ("AFE", 1.0, 0.9, True),
        ("Guest", 80.0, float("nan"), True),
        ("Guest", float("nan"), 0.0, False),
    ],
)
def test_the_tolerance_of_a_rule(
    rule: str, before: float, after: float, worse: bool
) -> None:
    row = dict.fromkeys(SCORECARD_COLUMNS, "") | {
        "model": "a",
        "context": "c",
        "verdict": "supported",
        "verified": "",
    }
    baseline = pd.DataFrame([row | {rule: before}])
    current = pd.DataFrame([row | {rule: after}])
    assert bool(compare(baseline, current)) is worse


def test_a_context_which_is_no_longer_verified() -> None:
    row = dict.fromkeys(SCORECARD_COLUMNS, "") | {
        "model": "a",
        "context": "c",
        "verdict": "supported",
    }
    baseline = pd.DataFrame([row | {"verified": "True"}])
    current = pd.DataFrame([row | {"verified": "False"}])
    (regression,) = compare(baseline, current)
    assert regression.what == "verified"


def test_the_qualification() -> None:
    table = qualification({"a": _assessment("a"), "b": _assessment("b", 2.0)})
    assert list(table["model"]) == ["all", "a", "b"]
    assert set(table["value"]) == {"RENAL_IMPAIRMENT"}
    pooled = table.iloc[0]
    assert pooled["n_pairs"] == 2
    assert pooled["n_within"] == 1
    assert pooled["within"] == 50.0
    assert pooled["gmfe"] == pytest.approx((1.1 * 2.0) ** 0.5)
    assert qualification({"a": _assessment("a")}, field="genotype").empty


def test_the_command(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    root = _write(tmp_path / "records", _assessment("a"))
    out = tmp_path / "scorecard.tsv"
    qualified = tmp_path / "qualification.tsv"
    assert (
        verify(
            [
                "scorecard",
                str(root),
                "--out",
                str(out),
                "--qualification",
                str(qualified),
            ]
        )
        == 0
    )
    printed = capsys.readouterr().out
    assert "1 models, 3 contexts of use" in printed
    assert "a/renal: supported" in printed
    assert "RENAL_IMPAIRMENT: 1 of 1 ratio pairs within the limits of Guest" in printed
    assert pd.read_csv(qualified, sep="\t")["n_pairs"].tolist() == [1, 1]
    # a worse renal context of a, and b which the baseline does not know
    worse = _write(tmp_path / "worse", _assessment("a", 2.0), _assessment("b"))
    assert verify(["scorecard", str(worse), "--baseline", str(out)]) == 1
    printed = capsys.readouterr().out
    assert "regressions against" in printed
    assert "a/renal: verdict supported -> not_supported" in printed
    assert verify(["scorecard", str(tmp_path / "nothing")]) == 2
