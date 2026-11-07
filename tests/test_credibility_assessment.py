"""The assessment of a model by sbmlsim with the extensions and the battery of pkpdsim.

`tests/data/credibility/toy-0.19.json` and `toy-uncertainty-0.19.json` are the
records pkpdsim 0.19 wrote for the toy credibility model (`run_credibility`
before the credibility moved into sbmlsim, the second with the propagation of
the uncertainty and the variability, without the verification); the assessment
of the model by the kind `credibility` gives the same verdicts and evidence.
"""

import math
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest
from sbmlsim.credibility import CredibilityAssessment
from sbmlsim.credibility.report import render_html

from pkpdsim.cli import credibility_command
from pkpdsim.credibility import EXTENSIONS, PK_BATTERY
from pkpdsim.credibility.ratios import section_pairs
from pkpdsim.model import PKPDModel

DATA = Path(__file__).parent / "data" / "credibility"

#: the outputs pkpdsim 0.19 wrote, every one of which the command writes again
FILES_0_19 = {
    "battery.tsv",
    "contexts.tsv",
    "credibility.json",
    "cross_validation.tsv",
    "domain.tsv",
    "evidence.tsv",
    "fates.tsv",
    "index.html",
    "input_ranges.tsv",
    "inputs.tsv",
    "quantities.tsv",
    "ratios.tsv",
    "studies.tsv",
    "temporal.tsv",
    "uncertainty_quantities.tsv",
    "uncertainty_sensitivity.tsv",
    "uncertainty.tsv",
    "variability_quantities.tsv",
    "variability.tsv",
    "verification.tsv",
    "vpc.tsv",
}

#: the records of the run of the kind next to the outputs
RUN_RECORDS = {"run.json", "spec.json", "summary.json", "progress.jsonl"}

#: the factors whose statements name the platform, the verification and the
#: battery, which are sbmlsim's now (the generic non-negativity of the battery
#: assesses the numerical code verification of the toy, which 0.19 did not)
PLATFORM_FACTORS = {
    "software_quality_assurance",
    "numerical_code_verification",
    "discretization_error",
    "numerical_solver_error",
    "use_error",
}


def _same(a: object, b: object) -> bool:
    """Compare two values of a record, numbers within the integration."""
    if isinstance(a, float) and isinstance(b, float):
        return (math.isnan(a) and math.isnan(b)) or math.isclose(
            a, b, rel_tol=1e-6, abs_tol=1e-12
        )
    return a == b


def _rows(a: list[Any], b: list[Any]) -> bool:
    return len(a) == len(b) and all(
        x.keys() == y.keys() and all(_same(x[k], y[k]) for k in x)
        for x, y in zip(a, b, strict=True)
    )


def _compare(old: CredibilityAssessment, new: CredibilityAssessment) -> None:
    """Compare the contexts, the evidence and the extensions of two records."""
    assert [c.id for c in new.contexts] == [c.id for c in old.contexts]
    for before, after in zip(old.contexts, new.contexts, strict=True):
        assert after.verdict == before.verdict, after.id
        assert (after.n_decisive, after.n_evaluable, after.minimum_points) == (
            before.n_decisive,
            before.n_evaluable,
            before.minimum_points,
        )
        assert _rows(list(after.acceptance), list(before.acceptance))
        assert after.quantities == before.quantities
        # the sections of the extensions, the fields of 0.19 read as sections;
        # a context without a ratio pairing has no section of the ratios, 0.19
        # wrote an empty list
        assert section_pairs(after.section("ratios")) == section_pairs(
            before.section("ratios")
        )
        for name in ("temporal", "inputs"):
            assert after.section(name) == before.section(name), (after.id, name)
        # the statement of a factor which is not assessed is the one of sbmlsim
        for a, b in zip(after.factors, before.factors, strict=True):
            assert a["id"] == b["id"]
            if a["id"] in PLATFORM_FACTORS:
                continue
            assert a["status"] == b["status"], a["id"]
            if a["status"] == "assessed":
                assert a["statement"] == b["statement"], a["id"]
    assert [e.to_dict() for e in new.evidence] == [e.to_dict() for e in old.evidence]
    assert [c.to_dict() for c in new.coverage] == [c.to_dict() for c in old.coverage]
    old_studies = old.extensions["studies"]["studies"]
    new_studies = new.extensions["studies"]["studies"]
    assert [(s["study"], s["n_points"], s["fits"]) for s in new_studies] == [
        (s["study"], s["n_points"], s["fits"]) for s in old_studies
    ]
    assert new.extensions["inputs"] == old.extensions["inputs"]


def test_the_assessment_gives_the_verdicts_of_pkpdsim_0_19(
    credibility_model: PKPDModel,
    as_model: Callable[[PKPDModel], PKPDModel],
    tmp_path: Path,
) -> None:
    model = as_model(credibility_model)
    out = tmp_path / "assessed"
    args = ["-o", str(out), "--source-commit", "0123abc"]
    assert credibility_command(model, args) == out
    assert {p.name for p in out.iterdir()} == FILES_0_19 | RUN_RECORDS
    old = CredibilityAssessment.from_json(DATA / "toy-0.19.json")
    new = CredibilityAssessment.from_json(out / "credibility.json")
    _compare(old, new)
    assert (new.provenance.commit, new.provenance.dirty) == ("0123abc", None)
    assert new.verification is not None and new.verification.battery is not None
    battery = {r.check.id: r.status.value for r in new.verification.battery.results}
    assert list(battery) == list(PK_BATTERY.order)
    # the toy has no substance with a molecular weight and does not start at
    # its steady state; the generic non-negativity runs on every model
    assert battery == {
        "mass_balance": "not_applicable",
        "non_negativity": "passed",
        "zero_dose": "not_applicable",
        "steady_state": "not_applicable",
        "superposition": "not_applicable",
        "cross_simulator": "not_applicable",
        "plausibility": "passed",
        "rebuild": "passed",
    }
    assert [c.verified for c in new.contexts] == [c.verified for c in old.contexts]


def test_the_uncertainty_and_the_variability_of_pkpdsim_0_19(
    credibility_model: PKPDModel,
    as_model: Callable[[PKPDModel], PKPDModel],
    tmp_path: Path,
) -> None:
    model = as_model(
        credibility_model.replace(credibility=".credibility:CREDIBILITY_VARIABILITY")
    )
    out = tmp_path / "uncertain"
    args = [
        *("-o", str(out), "--no-verification", "--uncertainty", "--draws", "8"),
        *("--variability", "--subjects", "8", "--source-commit", "0123abc"),
    ]
    credibility_command(model, args)
    # without the battery the table of the fates has no rows, as before
    assert {p.name for p in out.iterdir()} == FILES_0_19 | RUN_RECORDS
    old = CredibilityAssessment.from_json(DATA / "toy-uncertainty-0.19.json")
    new = CredibilityAssessment.from_json(out / "credibility.json")
    _compare(old, new)
    for before, after in zip(old.contexts, new.contexts, strict=True):
        assert before.uncertainty is not None and after.uncertainty is not None
        assert [r["decision"] for r in after.uncertainty["rules"]] == [
            r["decision"] for r in before.uncertainty["rules"]
        ]
        assert before.variability is not None and after.variability is not None
        assert after.variability["agreement"] == before.variability["agreement"]


@pytest.mark.parametrize("record", ["toy-0.19.json", "toy-uncertainty-0.19.json"])
def test_a_record_of_pkpdsim_0_19_is_reported(record: str) -> None:
    assessment = CredibilityAssessment.from_json(DATA / record)
    html = render_html(assessment, EXTENSIONS, PK_BATTERY)
    assert "Temporal validation" in html
    assert 'id="studies"' in html
    assert "Study data" in html
