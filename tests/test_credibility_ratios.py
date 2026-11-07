"""The ratio quantities and their limits of Guest et al. (`pkpdsim.credibility.ratios`).

The ratio pairs are the evidence extension `ratios`, assessed by
`sbmlsim.credibility.assess_table` with the extensions of pkpdsim.

The hand-built table of the test package has a renal impairment arm (`fm_renal`)
and a healthy arm (`fm_dose`) in the study `CredA`, five points each.
"""

import dataclasses
import math
from pathlib import Path
from types import MappingProxyType

import pandas as pd
import pytest
from credibility_tables import assess, fit_table, verification_of
from pkpdsim_credibility import MODEL
from pkpdsim_credibility.credibility import CREDIBILITY
from sbmlsim.credibility import (
    AcceptanceRecord,
    ContextOfUse,
    ContextResult,
    Credibility,
    CredibilityAssessment,
    DeclarationError,
    Evidence,
    EvidenceRole,
    FactorStatus,
    QuantityOfInterest,
    SelectionError,
    Verdict,
    select,
)
from sbmlsim.credibility.report import write_outputs
from sbmlsim.fit.filters import named_filter

from pkpdsim.credibility import EXTENSIONS, PK_BATTERY, RATIOS
from pkpdsim.credibility.ratios import (
    DEFAULT_RATIO_ACCEPTANCE,
    GuestAcceptance,
    Ratio,
    Summary,
    guest_limit,
    ratio_pairs,
    ratio_rule,
    section_pairs,
    summarize,
)
from pkpdsim.fit import metadata_filter
from pkpdsim.metadata import Health

RENAL = Ratio(
    observable="[A]",
    numerator=metadata_filter(health={Health.RENAL_IMPAIRMENT}),
    denominator=metadata_filter(health={Health.HEALTHY}),
)


def _renal(
    ratio: Ratio = RENAL,
    acceptance: GuestAcceptance | None = None,
    reference: tuple[str, ...] = ("fm_dose",),
) -> Credibility:
    """Get the declaration with the ratio of the renal context assessed.

    The healthy reference arms of `reference` are validation evidence of the
    context as well.
    """
    contexts = []
    for cou in CREDIBILITY.contexts:
        if cou.id == "renal":
            quantities = tuple(
                dataclasses.replace(q, ratio=ratio) if q.ratio else q
                for q in cou.quantities
            )
            healthy = Evidence(
                "healthy reference",
                EvidenceRole.VALIDATION,
                named_filter({"CredA": set(reference)}),
            )
            cou = dataclasses.replace(
                cou,
                quantities=quantities,
                evidence=(*cou.evidence, healthy),
                extensions={} if acceptance is None else {"ratios": acceptance},
                min_points=5,
            )
        contexts.append(cou)
    return dataclasses.replace(CREDIBILITY, contexts=tuple(contexts))


def _context(assessment: CredibilityAssessment, cou: str) -> ContextResult:
    return next(c for c in assessment.contexts if c.id == cou)


def _rule(context: ContextResult) -> AcceptanceRecord | None:
    """Get the row of the rule of the ratio pairs of a context, `None` without."""
    return next((a for a in context.acceptance if a.get("extension") == "ratios"), None)


def _pairs(context: ContextResult) -> list:
    return section_pairs(context.section("ratios"))


@pytest.mark.parametrize(
    ("observed", "limit"),
    [(1.0, 1.25), (2.0, 1.625), (0.5, 1.625), (10.0, 1.925)],
)
def test_guest_limit(observed: float, limit: float) -> None:
    assert guest_limit(observed) == pytest.approx(limit)


def test_guest_limit_of_no_ratio() -> None:
    assert math.isnan(guest_limit(0.0))
    assert math.isnan(guest_limit(math.nan))


def test_summarize() -> None:
    points = pd.DataFrame(
        {"x": [2.0, 0.0, 1.0], "DV": [2.0, 0.0, 2.0], "IPRED": [1.0, 1.0, 1.0]}
    )
    assert summarize(points, Summary.AUC) == (3.0, 2.0)
    assert all(math.isnan(v) for v in summarize(points, Summary.VALUE))
    assert summarize(points.iloc[[0]], Summary.VALUE) == (2.0, 1.0)
    assert all(math.isnan(v) for v in summarize(points.iloc[[0]], Summary.AUC))
    without_times = points.drop(columns="x")
    assert all(math.isnan(v) for v in summarize(without_times, Summary.AUC))


def test_the_pairs_of_a_study() -> None:
    fits = fit_table(factors={("CredA", "fm_renal"): 2.0}, times=True)
    arms = fits.table[fits.table["mapping"].isin(["fm_renal", "fm_dose"])]
    pairs, problems = ratio_pairs("A_ratio", RENAL, arms, fits.lookup)
    assert problems == []
    assert len(pairs) == 1
    pair = pairs[0]
    assert (pair["experiment"], pair["numerator"], pair["denominator"]) == (
        "CredA",
        "fm_renal",
        "fm_dose",
    )
    assert pair["observed"] == pytest.approx(1.0)
    assert pair["predicted"] == pytest.approx(2.0)
    assert pair["limit"] == pytest.approx(1.25)
    assert pair["within"] is False


def test_two_denominators_are_a_problem() -> None:
    fits = fit_table(times=True)
    both = Ratio(
        observable="[A]",
        numerator=metadata_filter(health={Health.RENAL_IMPAIRMENT}),
        denominator=lambda key, mapping: key != "fm_renal",
    )
    pairs, problems = ratio_pairs("A_ratio", both, fits.table, fits.lookup)
    assert pairs == []
    assert problems[0].startswith(
        "the ratio 'A_ratio' has more than one denominator in CredA (DECAY): "
    )


def test_a_ratio_within_the_limits_supports_the_context() -> None:
    fits = fit_table(factors={("CredA", "fm_renal"): 1.1}, times=True)
    renal = _context(assess(_renal(), fits), "renal")
    assert [p["within"] for p in _pairs(renal)] == [True]
    assert (rule := _rule(renal)) is not None
    assert rule["value"] == 100.0
    assert rule["met"] is True
    quantity = next(q for q in renal.quantities if q["ratio"])
    assert quantity["coverage"] == "covered"
    assert quantity["n_points"] == 1
    factors = {f["id"]: f for f in renal.factors}
    relevance = factors["relevance_of_quantities"]
    assert relevance["status"] == FactorStatus.ASSESSED.value
    assert "1 of 1 ratio pair within the limits of Guest" in relevance["statement"]
    assert factors["output_comparison"]["statement"].startswith(
        "4 of 4 acceptance rules"
    )


def test_a_ratio_outside_the_limits_does_not_support_the_context() -> None:
    fits = fit_table(factors={("CredA", "fm_renal"): 2.0}, times=True)
    renal = _context(assess(_renal(), fits), "renal")
    assert (rule := _rule(renal)) is not None
    assert rule["met"] is False
    assert renal.verdict is Verdict.NOT_SUPPORTED


def test_a_ratio_without_a_pair_is_insufficient_evidence() -> None:
    fits = fit_table(times=True)
    nobody = dataclasses.replace(
        RENAL, numerator=metadata_filter(health={Health.CIRRHOSIS})
    )
    renal = _context(assess(_renal(nobody), fits), "renal")
    assert _pairs(renal) == []
    assert (rule := _rule(renal)) is not None
    assert rule["met"] is None
    assert renal.verdict is Verdict.INSUFFICIENT_EVIDENCE
    quantity = next(q for q in renal.quantities if q["ratio"])
    assert quantity["coverage"] == "uncovered"


def test_a_declared_ratio_without_a_pairing_is_not_assessed() -> None:
    fits = fit_table(times=True)
    renal = _context(assess(CREDIBILITY, fits), "renal")
    assert _pairs(renal) == []
    assert _rule(renal) is None
    quantity = next(q for q in renal.quantities if q["ratio"])
    assert quantity["coverage"] is None


def test_the_default_and_a_declared_acceptance() -> None:
    renal = next(c for c in _renal().contexts if c.id == "renal")
    assert ratio_rule(renal) == DEFAULT_RATIO_ACCEPTANCE[renal.risk]
    declared = next(
        c for c in _renal(acceptance=GuestAcceptance(50.0)).contexts if c.id == "renal"
    )
    assert ratio_rule(declared) == GuestAcceptance(50.0)
    assert GuestAcceptance(50.0).limit == ">= 50"


def test_the_declaration_of_a_ratio() -> None:
    with pytest.raises(DeclarationError, match="between 0 and 100"):
        _renal(acceptance=GuestAcceptance(120.0))
    with pytest.raises(DeclarationError, match="needs the observable of its arms"):
        _renal(dataclasses.replace(RENAL, observable=" "))
    renal = next(c for c in _renal().contexts if c.id == "renal")
    only_ratio = dataclasses.replace(
        renal, quantities=tuple(q for q in renal.quantities if q.ratio)
    )
    with pytest.raises(DeclarationError, match="no ratio is required"):
        Credibility(contexts=(only_ratio,))


def test_a_ratio_of_an_unknown_observable() -> None:
    unknown = _renal(dataclasses.replace(RENAL, observable="B"))
    with pytest.raises(
        SelectionError, match="the ratio 'A_ratio' matches no observable"
    ):
        assess(unknown, fit_table(times=True))


def test_the_selection_refuses_two_denominators() -> None:
    two = _renal(reference=("fm_dose", "fm_subject"))
    with pytest.raises(SelectionError, match="more than one denominator in CredA"):
        select(two, fit_table(times=True), EXTENSIONS)


def test_the_outputs_of_the_ratios(tmp_path: Path) -> None:
    fits = fit_table(factors={("CredA", "fm_renal"): 1.1}, times=True)
    assessment = assess(_renal(), fits)
    table = RATIOS.report_tables(assessment.to_dict())["ratios"]
    assert list(table["context"]) == ["renal"]
    assert table.loc[0, "within"]
    record = tmp_path / "credibility.json"
    assessment.to_json(record)
    read = CredibilityAssessment.from_json(record)
    assert _pairs(_context(read, "renal")) == _pairs(_context(assessment, "renal"))
    output = write_outputs(assessment, tmp_path / "out", EXTENSIONS, PK_BATTERY)
    assert (
        (output / "ratios.tsv")
        .read_text(encoding="utf-8")
        .startswith("context\tquantity\t")
    )
    html = (output / "index.html").read_text(encoding="utf-8")
    assert "Ratio pairs" in html
    assert "fm_renal / fm_dose" in html


def test_the_verification_compares_the_ratio() -> None:
    fits = fit_table(factors={("CredA", "fm_renal"): 1.2}, times=True)
    verification = verification_of(fits)
    renal_rows = fits.table["mapping"] == "fm_renal"
    # the alternatives move the renal arm alone, the ratio by 1.1
    predictions = {
        check: series.where(~renal_rows, series * 1.1)
        for check, series in verification.predictions.items()
    }
    verification = dataclasses.replace(
        verification, predictions=MappingProxyType(predictions)
    )
    renal = _context(assess(_renal(), fits, verification), "renal")
    tolerances = next(r for r in renal.verification if r["check"] == "tolerances")
    # 1.2 is within 1.25 of 1, 1.32 with the alternative is not
    assert tolerances["acceptance_changed"] is True
    assert renal.verdict is Verdict.NOT_VERIFIED


def test_a_question_is_still_required() -> None:
    renal = next(c for c in _renal().contexts if c.id == "renal")
    assert isinstance(renal, ContextOfUse)
    assert any(isinstance(q.ratio, Ratio) for q in renal.quantities)
    assert any(not q.ratio for q in renal.quantities)
    assert isinstance(QuantityOfInterest("x", "y").ratio, bool)


def test_a_model_claims_the_ratio_pairings() -> None:
    """The extension `ratios` of a `PKPDModel` claims a `Ratio` and its acceptance."""
    declaration = _renal(acceptance=GuestAcceptance(50.0))
    assert MODEL.replace(credibility=declaration).load_credibility() is declaration
    assert RATIOS.claims(RENAL)
    assert not RATIOS.claims(QuantityOfInterest("x", "y"))  # ty: ignore[invalid-argument-type]


def test_the_ratio_pairs_pool_by_the_metadata() -> None:
    fits = fit_table(factors={("CredA", "fm_renal"): 1.1}, times=True)
    assessment = assess(_renal(), fits)
    sections = [c.section("ratios") for c in assessment.contexts]
    pooled = RATIOS.pool({"a": sections, "b": sections}, {})
    assert pooled is not None
    assert list(pooled["model"]) == ["all", "a", "b"]
    assert set(pooled["value"]) == {"RENAL_IMPAIRMENT"}
    assert list(pooled["n_pairs"]) == [2, 1, 1]
    by_route = RATIOS.pool({"a": sections}, {"field": "route"})
    assert by_route is not None
    assert by_route["value"].tolist() == ["IV", "IV"]
