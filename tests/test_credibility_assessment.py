"""The assessment of a model by sbmlsim with the extensions and the battery of pkpdsim.

`tests/data/credibility/toy-0.19.json` and `toy-uncertainty-0.19.json` are the
records pkpdsim 0.19 wrote for the toy credibility model (`run_credibility`
before the credibility moved into sbmlsim, the second with the propagation of
the uncertainty and the variability, without the verification),
`toy-extensions-0.19.json` the record of the toy with the ratio pairs of its
renal context, a study collected after the fit and the volume `Vc` checked
against the groups of a study (`_data`), and `battery-drug-0.19.json` the
battery of the drug model of `test_credibility_battery`; the assessment by the
kind `credibility` gives the same verdicts and evidence, the battery of sbmlsim
with the PK battery the same checks.
"""

import dataclasses
import json
import math
import re
from collections.abc import Callable
from datetime import date, timedelta
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from sbmlsim.credibility import Battery, CredibilityAssessment
from sbmlsim.credibility.battery import run_battery
from sbmlsim.credibility.report import render_html
from sbmlsim.fit.result import fit_date
from test_credibility_battery import _model
from test_credibility_ratios import _renal

from pkpdsim.cli import credibility_command
from pkpdsim.credibility import DOSE_STATES, EXTENSIONS, PK_BATTERY
from pkpdsim.credibility.ratios import section_pairs
from pkpdsim.data import StudyRecord, study_record_path, table_digests
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

#: the factor whose status and statement change: the generic non-negativity of
#: the battery of sbmlsim assesses the numerical code verification of every
#: model, which pkpdsim 0.19 did only for a model with a substance
CODE_VERIFICATION = "numerical_code_verification"

#: the texts of the statements of the factors which sbmlsim words otherwise,
#: the text of pkpdsim 0.19 by the one of sbmlsim
SBMLSIM_TEXTS = {
    "pkpdsim and sbmlsim are released packages": (
        "sbmlsim and the libraries of the model are released packages"
    ),
    "The builds are valid SBML": "The variants are valid SBML",
    "assesses it by the pairs of its arms;": (
        "assesses it by the evidence extension of the pairing;"
    ),
    # the correlation of a propagation comes from the posterior of a fit as well
    # since sbmlsim 0.20, which added the Laplace posterior as a source
    "with the correlation of the Fisher information": (
        "with the correlation of their posterior or Fisher information"
    ),
}


#: the checks pkpdsim added after 0.19: a record of 0.19 has neither their
#: result nor their sentence in the statement of their factor, so the stored
#: records are compared without them and the rest exactly. `dose_states` is
#: the check of 0.22, the dose states of a model the battery does not know
#: (matthiaskoenig/pkpdsim#115)
CHECKS_AFTER_0_19 = (DOSE_STATES,)

#: the sentences of `CHECKS_AFTER_0_19` in the statement of a factor, as
#: sbmlsim words them: the label of the check, its cases, its deviation and
#: its problems, up to its result
AFTER_0_19 = re.compile(
    "|".join(
        rf"{c.label[0].upper()}{re.escape(c.label[1:])} of \d+ cases?.*?"
        r": (?:passed|failed)\. "
        for c in CHECKS_AFTER_0_19
    )
)


def _before(statement: str) -> str:
    """Get a statement without the sentences of the checks added after 0.19."""
    return AFTER_0_19.sub("", statement)


def _worded(statement: str) -> str:
    """Get a statement of pkpdsim 0.19 in the words of sbmlsim."""
    for old, new in SBMLSIM_TEXTS.items():
        statement = statement.replace(old, new)
    return statement


#: a number of a statement: an integer, a decimal or a number with an exponent,
#: not a digit inside a word such as `C1` or `within_2fold`
NUMBER = re.compile(r"(?<![\w.])[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?(?![\w])")

#: the point of the largest relative difference of a check of the verification,
#: e.g. `, at CredA.fm_renal at 6,`, which another platform may find at another
#: point of equal difference within the noise of the integrator
LARGEST_AT = re.compile(r"(Largest relative difference .*?, at )[^,]+,")

#: the relative tolerance of a decimal number of a statement: the numbers with
#: a fraction or an exponent are the largest relative differences of the
#: verification, printed with two significant digits, the maximum over the
#: decisive points of the difference between two integrations at tolerances
#: tightened 100-fold; that is the noise of the integrator, whose size the
#: tolerances set and whose bits a platform (its libm, its compiler) changes,
#: so a value within a factor of 2 (`math.isclose` relative to the larger of
#: the two) is the same; an integer, a count of points, cases or rules, is
#: compared exactly
NOISE_RTOL = 0.5


def _statement_parts(statement: str) -> tuple[str, list[str]]:
    """Get a statement with its numbers as `<n>` and the location of a difference as `<at>`, and the numbers."""
    masked = LARGEST_AT.sub(r"\1<at>,", statement)
    return NUMBER.sub("<n>", masked), NUMBER.findall(masked)


def _same_statement(after: str, before: str) -> bool:
    """Compare two statements: the text exactly, the decimal numbers within the noise."""
    text_after, numbers_after = _statement_parts(after)
    text_before, numbers_before = _statement_parts(before)
    if text_after != text_before or len(numbers_after) != len(numbers_before):
        return False
    for a, b in zip(numbers_after, numbers_before, strict=True):
        decimal = any(c in a + b for c in ".eE")
        if decimal and not math.isclose(float(a), float(b), rel_tol=NOISE_RTOL):
            return False
        if not decimal and a != b:
            return False
    return True


def test_the_comparison_of_statements() -> None:
    before = (
        "Largest relative difference 2.5e-06 of 5 decisive predictions with the "
        "tolerances of the fits divided by 100, at CredA.fm_renal at 6, "
        "acceptance unchanged: passed."
    )
    moved = before.replace("2.5e-06", "3.1e-06").replace("at 6,", "at 8,")
    assert _same_statement(moved, before)
    assert not _same_statement(before.replace("2.5e-06", "9e-06"), before)
    assert not _same_statement(before.replace("of 5 decisive", "of 6 decisive"), before)
    assert not _same_statement(before.replace("passed", "failed"), before)
    assert _same_statement("C1 within_2fold 10 of 10.", "C1 within_2fold 10 of 10.")


#: the tolerance of a decimal number of a record: the metrics, the values of
#: the rules and the inputs of the adopted fit of the toy, which the fit of
#: the session gives again; another platform fits a slightly different rate
#: (macOS differed in the last digit of the GMFE, about 1e-16 relative), and
#: the least squares stops at its tolerance of about 1e-8 of the parameter, so
#: a metric within 1e-6 relative (1e-12 absolute) is the same, far below any
#: limit of a rule (1.3 of a GMFE, a percentage of points); an integer, a text
#: and a bool are compared exactly
FIT_RTOL = 1e-6


def _same(a: object, b: object) -> bool:
    """Compare two values of a record, decimal numbers within the noise of the fit.

    Mappings and sequences are compared element by element; `None` and `NaN`
    are each only themselves.
    """
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(_same(a[k], b[k]) for k in a)
    if isinstance(a, list | tuple) and isinstance(b, list | tuple):
        return len(a) == len(b) and all(_same(x, y) for x, y in zip(a, b, strict=True))
    decimal = isinstance(a, float) or isinstance(b, float)
    if decimal and isinstance(a, int | float) and isinstance(b, int | float):
        if isinstance(a, bool) or isinstance(b, bool):
            return a == b
        x, y = float(a), float(b)
        return (math.isnan(x) and math.isnan(y)) or math.isclose(
            x, y, rel_tol=FIT_RTOL, abs_tol=1e-12
        )
    return a == b


def _rows(a: list[Any], b: list[Any]) -> bool:
    return _same(list(a), list(b))


def test_the_comparison_of_records() -> None:
    assert _same({"GMFE": 1.0022358937189946}, {"GMFE": 1.002235893718994})
    assert not _same({"GMFE": 1.0022}, {"GMFE": 1.0023})
    assert not _same({"n": 10}, {"n": 11})
    assert _same([{"v": None, "w": math.nan}], [{"v": None, "w": math.nan}])
    assert not _same({"v": None}, {"v": 0.0})
    assert not _same({"met": True}, {"met": False})


def _compare(
    old: CredibilityAssessment,
    new: CredibilityAssessment,
    dated: Callable[[CredibilityAssessment, Any], Any] = lambda record, data: data,
) -> None:
    """Compare the contexts, the evidence and the extensions of two records.

    `dated` gives the sections of a record with its days relative to the day
    of its fit, for records of fits of other days.
    """
    assert [c.id for c in new.contexts] == [c.id for c in old.contexts]
    for before, after in zip(old.contexts, new.contexts, strict=True):
        assert after.verdict == before.verdict, after.id
        assert (after.n_decisive, after.n_evaluable, after.minimum_points) == (
            before.n_decisive,
            before.n_evaluable,
            before.minimum_points,
        )
        assert _rows(list(after.acceptance), list(before.acceptance))
        assert _same(list(after.quantities), list(before.quantities))
        # the sections of the extensions, the fields of 0.19 read as sections;
        # a context without a ratio pairing has no section of the ratios, 0.19
        # wrote an empty list
        assert _same(
            section_pairs(after.section("ratios")),
            section_pairs(before.section("ratios")),
        )
        for name in ("temporal", "inputs"):
            assert _same(
                dated(new, after.section(name)), dated(old, before.section(name))
            ), (after.id, name)
        for a, b in zip(after.factors, before.factors, strict=True):
            assert a["id"] == b["id"]
            if a["id"] == CODE_VERIFICATION:
                continue
            assert a.keys() == b.keys()
            assert (a["id"], a["status"], a["name"], a["step"]) == (
                b["id"],
                b["status"],
                b["name"],
                b["step"],
            )
            assert _same_statement(_before(a["statement"]), _worded(b["statement"])), (
                a["id"],
                a["statement"],
                b["statement"],
            )
    assert _same(
        [e.to_dict() for e in new.evidence], [e.to_dict() for e in old.evidence]
    )
    assert _same(
        [c.to_dict() for c in new.coverage], [c.to_dict() for c in old.coverage]
    )
    # the folder of a study is the one of the run
    old_studies = dated(old, old.extensions["studies"])["studies"]
    new_studies = dated(new, new.extensions["studies"])["studies"]
    assert _same(
        [{**s, "folder": None} for s in new_studies],
        [{**s, "folder": None} for s in old_studies],
    )
    assert _same(new.extensions["inputs"], old.extensions["inputs"])


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
    # its steady state; the generic non-negativity runs on every model, and the
    # dose states need neither a substance nor a simulation
    assert battery == {
        "dose_states": "passed",
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


def _data(data: Path, day: date) -> Path:
    """Get a data folder of the studies of the toy: the groups of CredA, the records of both.

    CredA is collected the day after the fit, CredB ten days before it; the
    volume of the group of CredA is 1050 ml.
    """
    a = data / "CredA"
    a.mkdir(parents=True)
    pd.DataFrame(
        [("all", "NA", "volume", 1050.0, "ml")],
        columns=["name", "parent", "measurement_type", "mean", "unit"],
    ).to_csv(a / ".CredA_TabGroups.tsv", sep="\t", index=False)
    for study, folder, collected in (
        ("CredA", a, day + timedelta(days=1)),
        ("CredB", data / "CredB", day - timedelta(days=10)),
    ):
        folder.mkdir(parents=True, exist_ok=True)
        record = StudyRecord(
            study=study,
            collected=collected,
            dating="collection",
            tables=table_digests(folder),
        )
        study_record_path(folder).write_text(json.dumps(record.to_dict()))
    return data


def _days(record: CredibilityAssessment, data: Any) -> Any:  # noqa: ANN401
    """Get a section with the days in it as days after the day of the fit of the record."""
    (fit_id,) = record.provenance.fits
    day = fit_date(fit_id)

    def offset(value: Any) -> Any:  # noqa: ANN401
        if isinstance(value, dict):
            return {k: offset(v) for k, v in value.items()}
        if isinstance(value, list):
            return [offset(v) for v in value]
        if isinstance(value, str):
            try:
                return (date.fromisoformat(value) - day).days
            except ValueError:
                return value
        return value

    return offset(data)


def test_the_extensions_give_the_evidence_of_pkpdsim_0_19(
    credibility_model: PKPDModel,
    as_model: Callable[[PKPDModel], PKPDModel],
    tmp_path: Path,
) -> None:
    """The ratio pairs, the temporal validation, the inputs and the study data."""
    (fit_id,) = credibility_model.adopted_ids
    data = _data(tmp_path / "data", fit_date(fit_id))
    declaration = dataclasses.replace(
        _renal(), covariates={"Vc": "volume", "BW": "weight", "HEIGHT": "height"}
    )
    model = as_model(
        credibility_model.replace(credibility=declaration, data_path=(data,))
    )
    out = tmp_path / "extensions"
    credibility_command(model, ["-o", str(out), "--source-commit", "0123abc"])
    old = CredibilityAssessment.from_json(DATA / "toy-extensions-0.19.json")
    new = CredibilityAssessment.from_json(out / "credibility.json")
    _compare(old, new, _days)
    renal = next(c for c in new.contexts if c.id == "renal")
    assert section_pairs(renal.section("ratios"))
    assert renal.acceptance[-1]["extension"] == "ratios"
    for context in new.contexts:
        inputs = context.section("inputs")
        assert inputs is not None and inputs["inputs"]["n_checked"] > 0
    dose = next(c for c in new.contexts if c.id == "dose")
    temporal = dose.section("temporal")
    assert temporal is not None and temporal["temporal"]["studies"]


#: the largest deviation of a check of the battery which passed is rounding or
#: noise of the integrator (the mass balance about 5e-15, the superposition
#: about 7e-10 at rtol 1e-10), orders below the limits of the checks (1e-6,
#: 0.01), which the status compares; another platform gives other bits of it,
#: so a deviation within a factor of 10 is the same, and zero (nothing
#: negative, a state which does not move) stays zero
NOISE_FACTOR = 10.0

#: the absolute tolerance of a part of the fate of a dose, relative to the dose:
#: the parts after two days integrated at rtol 1e-10, far below the parts the
#: fate records (above 1e-6 of the dose)
PART_ATOL = 1e-8


def _same_noise(a: float, b: float) -> bool:
    """Compare two largest deviations of the battery within the noise, `NaN` with `NaN`."""
    if math.isnan(a) or math.isnan(b):
        return math.isnan(a) and math.isnan(b)
    if a == 0.0 or b == 0.0:
        return a == b
    return 1 / NOISE_FACTOR <= a / b <= NOISE_FACTOR


def test_the_battery_of_the_drug_model_of_pkpdsim_0_19(tmp_path: Path) -> None:
    old = Battery.from_dict(
        json.loads((DATA / "battery-drug-0.19.json").read_text(encoding="utf-8"))
    )
    new = run_battery(_model(tmp_path), rebuild=False, cross_simulator=False)
    added = {c.id for c in CHECKS_AFTER_0_19}
    assert added <= {r.check.id for r in new.results}
    results = [r for r in new.results if r.check.id not in added]
    assert [r.check.id for r in results] == [r.check.id for r in old.results]
    for before, after in zip(old.results, results, strict=True):
        assert (after.status, after.cases) == (before.status, before.cases), after
        assert _same_noise(after.largest, before.largest), (
            after.check.id,
            after.largest,
            before.largest,
        )
    old_fates = old.sections["fates"]
    new_fates = new.sections["fates"]
    assert [(f["build"], f["dose"], f["substance"]) for f in new_fates] == [
        (f["build"], f["dose"], f["substance"]) for f in old_fates
    ]
    for before, after in zip(old_fates, new_fates, strict=True):
        assert after["amount"] == pytest.approx(before["amount"])
        assert after["parts"].keys() == before["parts"].keys()
        for part, value in before["parts"].items():
            assert after["parts"][part] == pytest.approx(
                value, rel=FIT_RTOL, abs=PART_ATOL * before["amount"]
            )
