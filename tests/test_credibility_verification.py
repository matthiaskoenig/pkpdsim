"""Tests of the verification of the predictions of a credibility assessment."""

import dataclasses
import json
import math
from pathlib import Path

import libsbml
import numpy as np
import pandas as pd
import pytest
from credibility_tables import fit_table, verification_of
from pkpdsim_credibility.credibility import CREDIBILITY
from sbmlsim.model.tolerances import AbsoluteTolerance

from pkpdsim.cli import credibility_command
from pkpdsim.credibility import (
    Check,
    CheckRecord,
    ContextResult,
    CredibilityAssessment,
    FactorStatus,
    FitTable,
    Provenance,
    Verdict,
    Verification,
    assess,
    assess_table,
    evaluate_fits,
)
from pkpdsim.credibility.assessment import verdict
from pkpdsim.credibility.declaration import AcceptanceResult, CouLevel
from pkpdsim.credibility.report import (
    VERIFICATION_COLUMNS,
    print_summary,
    render_html,
    verification_table,
)
from pkpdsim.credibility.verification import (
    LIMIT,
    SCALE_FLOOR,
    build_problems,
    relative_differences,
    tightened,
    verify,
)
from pkpdsim.model import PKPDModel
from pkpdsim_toy import MODEL_PATH

PROVENANCE = Provenance(
    model="credibility",
    title="Toy model with every kind of evidence",
    fits=["DECAY_fit"],
    versions={"pkpdsim": "0.0.0"},
    commit=None,
    dirty=None,
    date="2026-11-03T12:00:00+01:00",
)

VERIFICATION_FACTORS = ("discretization_error", "numerical_solver_error", "use_error")


def _assess(fits: FitTable, verification: Verification | None) -> CredibilityAssessment:
    return assess_table(CREDIBILITY, fits, PROVENANCE, verification)


def _context(assessment: CredibilityAssessment, id: str) -> ContextResult:
    (result,) = (c for c in assessment.contexts if c.id == id)
    return result


def _check(context: ContextResult, check: Check) -> CheckRecord:
    (record,) = (r for r in context.verification if r["check"] == check.value)
    return record


def test_tightened() -> None:
    assert tightened(1e-6) == AbsoluteTolerance(
        amount=1e-8, concentration=1e-8, other=1e-8
    )
    tolerance = AbsoluteTolerance(
        amount=1e-10, concentration=1e-6, other=1e-4, ids={"A": 1e-3}
    )
    result = tightened(tolerance, 10.0)
    assert result.amount == pytest.approx(1e-11)
    assert result.concentration == pytest.approx(1e-7)
    assert result.other == pytest.approx(1e-5)
    assert result.overrides == {"A": pytest.approx(1e-4)}


def test_relative_differences() -> None:
    table = pd.DataFrame(
        {
            "fit": "F",
            "experiment": "E",
            "mapping": ["m"] * 5 + ["n"] * 2,
            "IPRED": [1.0, 0.5, 0.0, 1e-9, np.nan, 2.0, 2.0],
        }
    )
    alternative = pd.Series([1.01, 0.5, 0.0, 2e-9, 1.0, np.nan, 2.2])
    differences = relative_differences(table, alternative)
    assert differences[0] == pytest.approx(0.01 / 1.01)
    assert differences[1] == 0.0
    # zeros are equal, a difference far below the peak of the curve is none
    assert differences[2] == 0.0
    assert differences[3] == pytest.approx(1e-9 / SCALE_FLOOR)
    # a reference which is not finite is the setup's, an alternative which is
    # not finite is an infinite difference
    assert math.isnan(differences[4])
    assert differences[5] == math.inf
    assert differences[6] == pytest.approx(0.2 / 2.2)


def test_build_problems(tmp_path: Path) -> None:
    assert build_problems({"male": MODEL_PATH}) == []

    doc = libsbml.SBMLDocument(3, 2)
    model = doc.createModel()
    model.setId("inconsistent")
    for sid, unit in (("p", "mole"), ("q", "litre")):
        parameter = model.createParameter()
        parameter.setId(sid)
        parameter.setValue(1.0)
        parameter.setUnits(unit)
        parameter.setConstant(sid == "q")
    rule = model.createAssignmentRule()
    rule.setVariable("p")
    rule.setMath(libsbml.parseL3Formula("q"))
    path = tmp_path / "inconsistent.xml"
    libsbml.writeSBMLToFile(doc, str(path))
    (problem,) = build_problems({"female": path})
    assert problem.startswith("the female build 'inconsistent.xml' has ")
    assert "[10513]" in problem


def test_a_verified_context() -> None:
    fits = fit_table()
    assessment = _assess(fits, verification_of(fits, 1.0 + 1e-6))
    dose = _context(assessment, "dose")
    assert dose.verdict is Verdict.SUPPORTED
    assert dose.verified is True
    assert [r["check"] for r in dose.verification] == [c.value for c in Check]
    tolerances = _check(dose, Check.TOLERANCES)
    assert tolerances["factor"] == "discretization_error"
    assert tolerances["n_points"] == dose.n_decisive == 10
    assert tolerances["max_difference"] == pytest.approx(1e-6, rel=1e-3)
    assert tolerances["acceptance_changed"] is False
    assert tolerances["passed"] is True
    setup = _check(dose, Check.SETUP)
    assert setup["acceptance_changed"] is None
    assert math.isnan(setup["max_difference"])
    assert setup["problems"] == []
    factors = {f["id"]: f for f in dose.factors}
    for factor in VERIFICATION_FACTORS:
        assert factors[factor]["status"] == FactorStatus.ASSESSED.value
        assert factors[factor]["step"] is None
        assert factors[factor]["statement"].endswith(": passed.")
    assert "acceptance unchanged" in factors["discretization_error"]["statement"]
    assert assessment.verification is not None
    assert assessment.verification.limit == LIMIT


def test_an_assessment_without_verification() -> None:
    assessment = _assess(fit_table(), None)
    assert assessment.verification is None
    for context in assessment.contexts:
        assert context.verification == []
        assert context.verified is None
        factors = {f["id"]: f for f in context.factors}
        for factor in VERIFICATION_FACTORS:
            assert factors[factor]["status"] == FactorStatus.NOT_ASSESSED.value
            assert factors[factor]["step"] == 2
    assert verification_table(assessment).empty


def test_a_difference_above_the_limit_is_not_verified() -> None:
    fits = fit_table()
    verification = verification_of(fits, factors={Check.INTEGRATOR: 1.02})
    dose = _context(_assess(fits, verification), "dose")
    assert dose.verdict is Verdict.NOT_VERIFIED
    assert dose.verdict.label == "not verified"
    assert dose.verified is False
    integrator = _check(dose, Check.INTEGRATOR)
    assert integrator["passed"] is False
    assert integrator["max_difference"] == pytest.approx(0.02 / 1.02)
    # the acceptance holds, the numbers are not accurate enough to say so
    assert integrator["acceptance_changed"] is False
    assert _check(dose, Check.TOLERANCES)["passed"] is True
    factors = {f["id"]: f for f in dose.factors}
    assert factors["numerical_solver_error"]["statement"].endswith(": failed.")
    # the limit is the one of the verification
    relaxed = verification_of(fits, factors={Check.INTEGRATOR: 1.02}, limit=0.05)
    assert _context(_assess(fits, relaxed), "dose").verdict is Verdict.SUPPORTED


def test_a_changed_acceptance_is_not_verified() -> None:
    """Predictions at the edge of a criterion: the alternative leaves the fold."""
    fits = fit_table(1.49)
    assert _context(_assess(fits, None), "dose").verdict is Verdict.SUPPORTED
    # the GMFE of the alternative is above 1.5
    verification = verification_of(fits, 1.02, limit=1.0)
    dose = _context(_assess(fits, verification), "dose")
    assert _check(dose, Check.TOLERANCES)["acceptance_changed"] is True
    assert dose.verdict is Verdict.NOT_VERIFIED


def test_not_verified_comes_before_not_supported() -> None:
    fits = fit_table(3.0)
    dose = _context(_assess(fits, verification_of(fits, 1.5)), "dose")
    assert dose.verdict is Verdict.NOT_VERIFIED
    met = [AcceptanceResult("within_2fold", 0.0, ">= 70", False)]
    assert verdict(
        10, met, True, "", level=CouLevel.C2, n_evaluable=10, minimum=1, verified=False
    ) is (Verdict.NOT_VERIFIED)
    # no decisive point is insufficient evidence, whatever the verification
    assert verdict(
        0, met, True, "", level=CouLevel.C2, n_evaluable=0, minimum=1, verified=False
    ) is (Verdict.INSUFFICIENT_EVIDENCE)
    assert verdict(
        10, met, True, "", level=CouLevel.C2, n_evaluable=10, minimum=1, verified=True
    ) is (Verdict.NOT_SUPPORTED)


def test_a_failed_evaluation_fails_its_check() -> None:
    fits = fit_table()
    message = "the evaluation of the fit 'DECAY' failed: CVODE error"
    verification = verification_of(
        fits, failures={Check.INTEGRATOR: {"DECAY": message}}
    )
    dose = _context(_assess(fits, verification), "dose")
    integrator = _check(dose, Check.INTEGRATOR)
    assert integrator["passed"] is False
    assert integrator["problems"] == [message]
    assert dose.verdict is Verdict.NOT_VERIFIED
    # a failure of another fit is not a problem of the context
    other = verification_of(fits, failures={Check.INTEGRATOR: {"PD": message}})
    assert _context(_assess(fits, other), "dose").verified is True


def test_the_setup() -> None:
    fits = fit_table()
    problem = "the male build 'toy.xml' has 1 issue of the SBML: [10513] units"
    dose = _context(_assess(fits, verification_of(fits, problems=(problem,))), "dose")
    setup = _check(dose, Check.SETUP)
    assert setup["passed"] is False
    assert setup["problems"] == [problem]
    assert dose.verdict is Verdict.NOT_VERIFIED

    table = fits.table.copy()
    rows = table.index[
        (table["mapping"] == "fm_dose") & (table["experiment"] == "CredA")
    ]
    table.loc[rows[:2], "IPRED"] = np.nan
    failed = dataclasses.replace(fits, table=table)
    dose = _context(_assess(failed, verification_of(failed)), "dose")
    setup = _check(dose, Check.SETUP)
    assert setup["problems"] == [
        "2 of 10 decisive predictions are not finite, of CredA.fm_dose"
    ]
    assert dose.verdict is Verdict.NOT_VERIFIED
    # a reference which is not finite is no numerical difference
    assert _check(dose, Check.TOLERANCES)["n_points"] == 8


def test_the_json_round_trip_of_a_verified_assessment(tmp_path: Path) -> None:
    fits = fit_table()
    assessment = _assess(fits, verification_of(fits))
    path = tmp_path / "credibility.json"
    assessment.to_json(path)
    data = json.loads(path.read_text())
    (setup,) = (r for r in data["contexts"][0]["verification"] if r["check"] == "setup")
    assert setup["max_difference"] is None
    read = CredibilityAssessment.from_json(path)
    assert read.to_dict() == assessment.to_dict()
    assert read.verification == assessment.verification

    # a record of pkpdsim 0.17 was not verified
    del data["verification"]
    for context in data["contexts"]:
        del context["verification"]
    path.write_text(json.dumps(data))
    old = CredibilityAssessment.from_json(path)
    assert old.verification is None
    assert all(c.verified is None for c in old.contexts)


def test_the_report_of_a_verification() -> None:
    fits = fit_table()
    verification = verification_of(fits, factors={Check.INTEGRATOR: 1.02})
    assessment = _assess(fits, verification)
    table = verification_table(assessment)
    assert list(table.columns) == list(VERIFICATION_COLUMNS)
    assert len(table) == len(assessment.contexts) * len(Check)
    html = render_html(assessment)
    section = html[html.index('id="verification"') : html.index('id="provenance"')]
    assert "CVODE BDF of order at most 2, rtol 1e-08" in section
    assert "Every build is valid SBML with consistent units." in section
    assert "second integration method" in html
    assert 'class="state failed">failed<' in html
    assert 'class="verdict not_verified">not verified<' in html
    unverified = render_html(_assess(fits, None))
    assert "The assessment was not verified." in unverified


def test_verify_needs_the_reports() -> None:
    with pytest.raises(ValueError, match="no reports to verify"):
        verify(PKPDModel.__new__(PKPDModel), fit_table())


def test_verify_the_test_package(credibility_model: PKPDModel) -> None:
    assessment = assess(credibility_model)
    assert assessment.verification is not None
    (fit,) = assessment.verification.fits
    assert fit["reference"].startswith("CVODE BDF, rtol 1e-06")
    assert fit["tolerances"].startswith("CVODE BDF, rtol 1e-08")
    assert fit["integrator"].startswith("CVODE BDF of order at most 2, rtol 1e-08")
    assert assessment.verification.problems == []
    for context in assessment.contexts:
        assert context.verified is True, context.verification
        for record in context.verification[:2]:
            assert record["n_points"] > 0
            assert record["max_difference"] < LIMIT
    # the verdicts of the assessment without verification, see
    # `test_assess_the_test_package`
    assert {c.id: c.verdict for c in assessment.contexts} == {
        "describe": Verdict.SUPPORTED,
        "dose": Verdict.SUPPORTED,
        "renal": Verdict.INSUFFICIENT_EVIDENCE,
    }
    unverified = assess(credibility_model, verification=False)
    assert unverified.verification is None
    assert [c.verdict for c in unverified.contexts] == [
        c.verdict for c in assessment.contexts
    ]


def test_the_command_without_verification(
    credibility_model: PKPDModel, tmp_path: Path
) -> None:
    output_dir = credibility_command(
        credibility_model, ["-o", str(tmp_path / "out"), "--no-verification"]
    )
    record = CredibilityAssessment.from_json(output_dir / "credibility.json")
    assert record.verification is None
    assert (output_dir / "verification.tsv").read_text().strip() == "\t".join(
        VERIFICATION_COLUMNS
    )
    with pytest.raises(SystemExit) as exc:
        credibility_command(
            credibility_model,
            ["--from-json", str(output_dir / "credibility.json"), "--no-verification"],
        )
    assert exc.value.code == 2


def test_the_summary_tells_whether_the_contexts_are_verified(
    capsys: pytest.CaptureFixture[str],
) -> None:
    fits = fit_table()
    print_summary(_assess(fits, None))
    assert "not verified (--no-verification)" in capsys.readouterr().out
    print_summary(_assess(fits, verification_of(fits)))
    assert "3 of 3 contexts of use verified" in capsys.readouterr().out


def test_a_failed_simulation_of_an_alternative_is_a_failure(
    credibility_model: PKPDModel, monkeypatch: pytest.MonkeyPatch
) -> None:
    """sbmlsim raises a `ValueError` for the data of a failed simulation."""
    import pkpdsim.credibility.verification as module

    def failing(report: object, check: Check, tightening: float) -> None:
        if check is Check.INTEGRATOR:
            raise ValueError("'fm_dose': no simulation results")
        raise RuntimeError("CVODE Error: CV_TOO_MUCH_WORK")

    fits = evaluate_fits(credibility_model)
    monkeypatch.setattr(module, "_alternative_points", failing)
    verification = verify(credibility_model, fits)
    (fit,) = fits.reports
    assert verification.failures[Check.TOLERANCES] == {
        fit: f"the evaluation of the fit '{fit}' failed: CVODE Error: CV_TOO_MUCH_WORK"
    }
    assert verification.failures[Check.INTEGRATOR][fit].endswith(
        "no simulation results"
    )
    assert verification.predictions[Check.INTEGRATOR].isna().all()
    assessment = assess_table(CREDIBILITY, fits, PROVENANCE, verification)
    for context in assessment.contexts:
        if context.n_decisive:
            assert context.verdict is Verdict.NOT_VERIFIED
