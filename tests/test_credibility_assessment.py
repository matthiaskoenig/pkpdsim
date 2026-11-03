"""Tests of the assessment of the contexts of use."""

import dataclasses
import json
import math
from pathlib import Path

import pytest
from credibility_tables import fit_table
from pkpdsim_credibility import MODEL
from pkpdsim_credibility.credibility import CREDIBILITY
from pkpdsim_credibility.fitting import DOSE, ZERO
from release_git import git, init_repo

from pkpdsim.credibility import (
    UNCLAIMED,
    AcceptanceResult,
    ContextOfUse,
    ContextResult,
    CouLevel,
    Coverage,
    CoverageEntry,
    Credibility,
    CredibilityAssessment,
    Domain,
    Evidence,
    EvidenceRole,
    FactorStatus,
    FitTable,
    Provenance,
    QuantityOfInterest,
    Unclaimed,
    Verdict,
    assess,
    assess_table,
    describe_selection,
)
from pkpdsim.credibility.assessment import provenance, verdict
from pkpdsim.credibility.factors import FACTORS
from pkpdsim.fit import NamedMappings, metadata_filter, named_filter
from pkpdsim.metadata import Health, Route
from pkpdsim.model import PKPDModel

PROVENANCE = Provenance(
    model="credibility",
    title="Toy model with every kind of evidence",
    fits=["DECAY_fit"],
    versions={"pkpdsim": "0.0.0"},
    commit=None,
    dirty=None,
    date="2026-11-03T12:00:00+01:00",
)


def _assess(
    fits: FitTable, credibility: Credibility = CREDIBILITY
) -> CredibilityAssessment:
    """Assess a declaration on a table."""
    return assess_table(credibility, fits, PROVENANCE)


def _with(*contexts: ContextOfUse) -> Credibility:
    """Get the declaration of the test package with other contexts."""
    return dataclasses.replace(CREDIBILITY, contexts=contexts)


def _context(assessment: CredibilityAssessment, id: str) -> ContextResult:
    """Get the result of a context of use."""
    (result,) = (c for c in assessment.contexts if c.id == id)
    return result


def _coverage(
    assessment: CredibilityAssessment, context: str, field: str, value: str
) -> CoverageEntry:
    """Get the coverage of a value of a context of use."""
    (entry,) = (
        e
        for e in assessment.coverage
        if (e.context, e.field, e.value) == (context, field, value)
    )
    return entry


def test_a_supported_context() -> None:
    dose = _context(_assess(fit_table()), "dose")
    assert dose.verdict is Verdict.SUPPORTED
    assert dose.verdict.label == "supported"
    assert dose.acceptance
    assert all(result["met"] is True for result in dose.acceptance)
    # fm_dose and fm_subject, 5 points each
    assert dose.n_decisive == 10


def test_a_not_supported_context() -> None:
    dose = _context(_assess(fit_table(3.0)), "dose")
    assert dose.verdict is Verdict.NOT_SUPPORTED
    assert dose.verdict.label == "not supported"
    assert any(result["met"] is False for result in dose.acceptance)


def test_an_uncovered_value_without_extrapolation() -> None:
    fits = fit_table()
    table = fits.table.assign(
        health=fits.table.health.replace("RENAL_IMPAIRMENT", "HEALTHY")
    )
    fits = dataclasses.replace(fits, table=table)
    assessment = _assess(fits)
    renal = _context(assessment, "renal")
    assert renal.verdict is Verdict.INSUFFICIENT_EVIDENCE
    assert all(result["met"] is True for result in renal.acceptance)
    uncovered = _coverage(assessment, "renal", "health", "RENAL_IMPAIRMENT")
    assert uncovered.coverage is Coverage.UNCOVERED
    assert uncovered.claimed
    assert uncovered.n_points == 0
    healthy = _coverage(assessment, "renal", "health", "HEALTHY")
    assert healthy.coverage is Coverage.VALIDATED
    assert healthy.n_points == 5

    cou = dataclasses.replace(
        CREDIBILITY.context("renal"), extrapolation="the kidney is a clearance"
    )
    renal = _context(_assess(fits, _with(cou)), "renal")
    assert renal.verdict is Verdict.SUPPORTED_WITH_EXTRAPOLATION
    assert renal.verdict.label == "supported with extrapolation"


def test_a_c4_context_is_never_plainly_supported() -> None:
    """A C4 context extrapolates by definition, even on validation evidence."""
    cou = dataclasses.replace(
        CREDIBILITY.context("dose"),
        level=CouLevel.C4,
        extrapolation="the qualification of the platform",
    )
    assessment = _assess(fit_table(), _with(cou))
    dose = _context(assessment, "dose")
    assert all(result["met"] is True for result in dose.acceptance)
    assert all(
        e.coverage is Coverage.VALIDATED
        for e in assessment.coverage
        if e.context == "dose"
    )
    assert dose.verdict is Verdict.SUPPORTED_WITH_EXTRAPOLATION

    failed = _context(_assess(fit_table(3.0), _with(cou)), "dose")
    assert failed.verdict is Verdict.NOT_SUPPORTED


def test_the_verdict_of_a_level() -> None:
    met = [AcceptanceResult("within_2fold", 100.0, ">= 70", True)]
    for level in CouLevel:
        expected = (
            Verdict.SUPPORTED_WITH_EXTRAPOLATION
            if level is CouLevel.C4
            else Verdict.SUPPORTED
        )
        assert verdict(10, met, True, "argued", level=level) is expected
    assert verdict(10, met, False, "", level=CouLevel.C2) is (
        Verdict.INSUFFICIENT_EVIDENCE
    )


def test_a_ratio_quantity_is_recorded() -> None:
    ratio = QuantityOfInterest("A_ratio", "impaired over healthy", ratio=True)
    cou = dataclasses.replace(
        CREDIBILITY.context("dose"),
        quantities=(*CREDIBILITY.context("dose").quantities, ratio),
    )
    assessment = _assess(fit_table(), _with(cou))
    dose = _context(assessment, "dose")
    assert dose.quantities[-1] == {
        "id": "A_ratio",
        "description": "impaired over healthy",
        "ratio": True,
    }
    assert dose.verdict is Verdict.SUPPORTED
    assert dose.n_decisive == _context(_assess(fit_table()), "dose").n_decisive


def test_a_point_without_the_field_is_no_value_outside_the_claim() -> None:
    fits = fit_table()
    table = fits.table.copy()
    subject = (table.experiment == "CredA") & (table.mapping == "fm_subject")
    table["health"] = table["health"].astype(object)
    table.loc[subject, "health"] = None
    assessment = _assess(dataclasses.replace(fits, table=table))
    assert not [e for e in assessment.coverage if e.value == ""]
    assert _coverage(assessment, "dose", "health", "HEALTHY").n_points == 5


@pytest.mark.parametrize(
    ("selected", "factors"),
    [
        pytest.param(ZERO, {}, id="zeros"),
        pytest.param(DOSE, {("CredA", "fm_dose"): math.nan}, id="failed simulations"),
    ],
)
def test_a_criterion_without_applicable_points(
    selected: NamedMappings, factors: dict[tuple[str, str], float]
) -> None:
    cou = dataclasses.replace(
        CREDIBILITY.context("dose"),
        evidence=(
            Evidence("no ratios", EvidenceRole.VALIDATION, named_filter(selected)),
        ),
    )
    result = _context(_assess(fit_table(factors=factors), _with(cou)), "dose")
    (within,) = (r for r in result.acceptance if r["name"] == "within_2fold")
    assert within["met"] is None
    assert math.isnan(within["value"])
    assert result.n_decisive == 5
    assert result.n_evaluable == 0
    assert result.verdict is Verdict.INSUFFICIENT_EVIDENCE
    assert result.verdict.label == "insufficient evidence"


def _without_ratios(
    fits: FitTable, mapping: str, column: str, value: float
) -> FitTable:
    """Get the table with a value of a column for the rows of a mapping of CredA."""
    table = fits.table.copy()
    table.loc[(table.experiment == "CredA") & (table.mapping == mapping), column] = (
        value
    )
    return dataclasses.replace(fits, table=table)


@pytest.mark.parametrize(
    ("column", "value"),
    [
        pytest.param("IPRED", math.nan, id="failed simulations"),
        pytest.param("DV", 0.0, id="zeros"),
    ],
)
def test_points_without_a_ratio_cover_nothing(column: str, value: float) -> None:
    cou = dataclasses.replace(
        CREDIBILITY.context("renal"),
        evidence=(
            Evidence(
                "doses and renal",
                EvidenceRole.VALIDATION,
                named_filter({"CredA": {"fm_dose", "fm_renal"}}),
            ),
        ),
    )
    fits = _without_ratios(fit_table(), "fm_renal", column, value)
    assessment = _assess(fits, _with(cou))
    renal = _context(assessment, "renal")
    impaired = _coverage(assessment, "renal", "health", "RENAL_IMPAIRMENT")
    assert impaired.coverage is Coverage.UNCOVERED
    assert impaired.n_points == 0
    assert _coverage(assessment, "renal", "health", "HEALTHY").n_points == 5
    assert (renal.n_decisive, renal.n_evaluable) == (10, 5)
    assert all(result["met"] is True for result in renal.acceptance)
    assert renal.verdict is Verdict.INSUFFICIENT_EVIDENCE

    covered = _assess(fit_table(), _with(cou))
    assert _context(covered, "renal").verdict is Verdict.SUPPORTED


def test_no_decisive_evidence_is_insufficient() -> None:
    cou = dataclasses.replace(
        CREDIBILITY.context("dose"),
        evidence=(
            Evidence(
                "fitted", EvidenceRole.CALIBRATION, metadata_filter(route={Route.IV})
            ),
        ),
        extrapolation="the decay is linear",
    )
    result = _context(_assess(fit_table(), _with(cou)), "dose")
    assert result.n_decisive == 0
    assert result.verdict is Verdict.INSUFFICIENT_EVIDENCE


def test_c1_decides_on_calibration() -> None:
    cou = dataclasses.replace(
        CREDIBILITY.context("describe"),
        domain=Domain(metadata={"route": (Route.IV,)}),
    )
    assessment = _assess(fit_table(), _with(cou))
    describe = _context(assessment, "describe")
    assert describe.verdict is Verdict.SUPPORTED
    # CredA.fm_train and CredB.fm_dose
    assert describe.n_decisive == 10
    iv = _coverage(assessment, "describe", "route", "IV")
    assert iv.coverage is Coverage.CALIBRATED
    assert iv.n_points == 10

    # a C2 context on the same evidence is not decided on calibration
    interpolate = dataclasses.replace(
        cou, level=CouLevel.C2, extrapolation="the decay is linear"
    )
    result = _context(_assess(fit_table(), _with(interpolate)), "describe")
    assert result.n_decisive == 0
    assert result.verdict is Verdict.INSUFFICIENT_EVIDENCE


def test_the_limitations_are_reported_apart() -> None:
    assessment = _assess(fit_table())
    (limitation,) = assessment.limitations
    assert limitation.context == ""
    assert limitation.name == "scaled data"
    assert limitation.role is EvidenceRole.LIMITATION
    assert limitation.n_points == 2
    assert limitation.n_mappings == 1
    assert limitation.metrics["within_2fold"] == 0.0
    for factor in (0.01, 1.0, 3.0, 100.0):
        other = _assess(fit_table(factors={("CredA", "fm_limit"): factor}))
        assert [c.acceptance for c in other.contexts] == [
            c.acceptance for c in assessment.contexts
        ]


def test_the_limitations_of_a_context() -> None:
    cou = dataclasses.replace(
        CREDIBILITY.context("dose"),
        limitations=(
            Evidence(
                "scaled",
                EvidenceRole.LIMITATION,
                named_filter({"CredA": {"fm_limit"}}),
            ),
        ),
    )
    assessment = _assess(fit_table(), _with(cou))
    assert [(e.context, e.name) for e in assessment.limitations] == [
        ("", "scaled data"),
        ("dose", "scaled"),
    ]


def test_the_evidence() -> None:
    assessment = _assess(fit_table())
    (doses,) = (e for e in assessment.evidence if e.name == "doses")
    assert doses.context == "dose"
    assert doses.role is EvidenceRole.VALIDATION
    assert doses.selection == "CredA: fm_dose"
    assert (doses.n_points, doses.n_mappings, doses.n_studies) == (5, 1, 1)
    assert doses.metrics["GMFE"] == 1.0
    assert doses.metrics["within_2fold"] == 100.0
    (other,) = (e for e in assessment.evidence if e.name == "other healthy data")
    # the limitation is healthy and intravenous as well
    assert other.limitation_subtracted == 2
    (fitted,) = (e for e in assessment.evidence if e.context == "describe")
    assert (fitted.n_mappings, fitted.n_studies) == (2, 2)
    assert assessment.unclaimed == [
        {"fit": "DECAY", "experiment": "CredB", "mapping": "fm_zero", "n": 5}
    ]


def test_the_values_outside_the_claim() -> None:
    cou = dataclasses.replace(
        CREDIBILITY.context("dose"),
        evidence=(
            Evidence(
                "doses",
                EvidenceRole.VALIDATION,
                named_filter({"CredA": {"fm_dose", "fm_renal"}}),
            ),
        ),
    )
    assessment = _assess(fit_table(), _with(cou))
    outside = _coverage(assessment, "dose", "health", "RENAL_IMPAIRMENT")
    assert outside == CoverageEntry(
        context="dose",
        field="health",
        value="RENAL_IMPAIRMENT",
        coverage=None,
        claimed=False,
        n_points=5,
    )
    assert _coverage(assessment, "dose", "health", "HEALTHY").claimed


def test_the_factors() -> None:
    assessment = _assess(fit_table())
    assert len(FACTORS) == 13
    for context in assessment.contexts:
        assert [f["id"] for f in context.factors] == [f.id for f in FACTORS]
        assert [f["status"] for f in context.factors] == [
            f.status.value for f in FACTORS
        ]
    statuses = {f.id: (f.status, f.step) for f in FACTORS}
    assert statuses["numerical_code_verification"] == (FactorStatus.PLATFORM, None)
    assert statuses["model_inputs"] == (FactorStatus.NOT_ASSESSED, 4)
    assert statuses["output_comparison"] == (FactorStatus.ASSESSED, None)


def test_the_json_round_trip(tmp_path: Path) -> None:
    cou = dataclasses.replace(
        CREDIBILITY.context("dose"),
        evidence=(Evidence("zeros", EvidenceRole.VALIDATION, named_filter(ZERO)),),
    )
    assessment = _assess(fit_table(), _with(cou, CREDIBILITY.context("describe")))
    (zeros,) = (e for e in assessment.evidence if e.name == "zeros")
    assert math.isnan(zeros.metrics["GMFE"])
    path = tmp_path / "credibility.json"
    assessment.to_json(path)
    data = json.loads(path.read_text())
    (written,) = (e for e in data["evidence"] if e["name"] == "zeros")
    assert written["metrics"]["GMFE"] is None
    read = CredibilityAssessment.from_json(path)
    assert read.to_dict() == assessment.to_dict()
    (zeros,) = (e for e in read.evidence if e.name == "zeros")
    assert math.isnan(zeros.metrics["GMFE"])
    assert read.contexts[0].verdict is Verdict.INSUFFICIENT_EVIDENCE
    assert read.provenance == PROVENANCE


def test_describe_selection() -> None:
    assert (
        describe_selection(
            named_filter({"CredB": {"fm_zero", "fm_dose"}, "CredA": {"fm_dose"}})
        )
        == "CredA: fm_dose; CredB: fm_dose, fm_zero"
    )
    assert (
        describe_selection(
            metadata_filter(
                route={Route.PO, Route.IV}, excluded={"health": {Health.HEALTHY}}
            )
        )
        == "route in [IV, PO] and health not in [HEALTHY]"
    )
    assert describe_selection(UNCLAIMED) == "unclaimed validation data"
    assert (
        describe_selection(Unclaimed(within=named_filter(DOSE)))
        == "unclaimed validation data within CredA: fm_dose"
    )

    def by_key(key: str, fit_mapping: object) -> bool:
        return key == "fm_dose"

    assert describe_selection(by_key) == by_key.__qualname__


def test_assess_the_test_package(credibility_model: PKPDModel) -> None:
    assessment = assess(credibility_model)
    assert [c.id for c in assessment.contexts] == ["describe", "dose", "renal"]
    assert assessment.provenance.fits == list(credibility_model.fits)
    assert assessment.provenance.model == credibility_model.name
    assert set(assessment.provenance.versions) == {
        "pkpdsim",
        "sbmlsim",
        "sbmlutils",
        "libroadrunner",
    }
    # renal claims healthy subjects, its evidence is of renal impairment only
    assert {c.id: c.verdict for c in assessment.contexts} == {
        "describe": Verdict.SUPPORTED,
        "dose": Verdict.SUPPORTED,
        "renal": Verdict.INSUFFICIENT_EVIDENCE,
    }
    renal = _context(assessment, "renal")
    assert all(result["met"] is True for result in renal.acceptance)
    assert _coverage(assessment, "renal", "health", "HEALTHY").coverage is (
        Coverage.UNCOVERED
    )


@pytest.mark.parametrize("verdict", list(Verdict))
def test_the_labels(verdict: Verdict) -> None:
    assert verdict.label == verdict.value.replace("_", " ")


def test_the_provenance_of_a_release() -> None:
    """In a release the source commit is given, whether its tree was clean is unknown."""
    result = provenance(MODEL, source_commit="0123abc")
    assert (result.commit, result.dirty) == ("0123abc", None)
    assert result.model == MODEL.name


@pytest.mark.usefixtures("git_env")
def test_the_provenance_is_of_the_folder_of_the_model(tmp_path: Path) -> None:
    root = tmp_path / "repository"
    init_repo(root)
    package = root / "model"
    package.mkdir()
    (package / "credibility.py").write_text("CREDIBILITY = None\n")
    git(root, "add", "-A")
    git(root, "commit", "-m", "init")
    model = dataclasses.replace(MODEL, path=package)
    (root / "other.txt").write_text("a change outside the model\n")
    result = provenance(model)
    assert result.commit == git(root, "rev-parse", "HEAD").strip()
    assert result.dirty is False
    (package / "notes.txt").write_text("a change of the model\n")
    assert provenance(model).dirty is True
    assert provenance(dataclasses.replace(model, path=tmp_path)).commit is None
