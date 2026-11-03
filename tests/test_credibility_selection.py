"""Tests of the kinds of the adopted fits and of the selection of the evidence."""

import dataclasses
import shutil
from pathlib import Path

import pandas as pd
import pytest
from pkpdsim_credibility import MODEL
from pkpdsim_credibility.credibility import CREDIBILITY
from pkpdsim_credibility.experiments import CredA, CredB
from pkpdsim_credibility.fitting import DOSE
from sbmlsim.fit import FitMapping
from sbmlsim.fit.helpers import FitMappings
from sbmlsim.fit.metrics import metadata_dict

from pkpdsim.credibility import (
    UNCLAIMED,
    ContextOfUse,
    CouLevel,
    Credibility,
    Domain,
    Evidence,
    EvidenceRole,
    FitTable,
    KindMismatchError,
    Level,
    QuantityOfInterest,
    SelectionError,
    Selections,
    Unclaimed,
    check_declaration,
    evaluate_fits,
    kind_differences,
    mapping_table,
    normalize_observable,
    recorded_kinds,
    select,
)
from pkpdsim.fit import metadata_filter, named_filter
from pkpdsim.metadata import Health, Route
from pkpdsim.model import PKPDModel

#: the kinds of the adopted fit of the test package, see `test_credibility_package`
KINDS = {
    ("CredA", "fm_train"): "training",
    ("CredA", "fm_dose"): "validation",
    ("CredA", "fm_renal"): "validation",
    ("CredA", "fm_limit"): "validation",
    ("CredA", "fm_subject"): "validation",
    ("CredA", "fm_outlier"): "outlier",
    ("CredB", "fm_dose"): "training",
    ("CredB", "fm_zero"): "validation",
}
#: the number of data points of a mapping of the hand-built table
POINTS = {("CredA", "fm_limit"): 2, ("CredA", "fm_outlier"): 3}
_A = (QuantityOfInterest("A", "the concentration of A"),)


@pytest.fixture(scope="module")
def fits() -> FitTable:
    """Get a table of the data points of the kinds of the test package."""
    mappings = FitMappings(
        experiment_classes=[CredA, CredB],
        base_path=MODEL.path,
        data_path=list(MODEL.data_paths),
    )

    def lookup(experiment: str, mapping: str) -> FitMapping:
        return mappings.experiments[experiment]._fit_mappings[mapping]

    rows = [
        {
            "fit": "DECAY",
            "experiment": experiment,
            "mapping": mapping,
            "kind": kind,
            "observable": "[A]",
            **metadata_dict(lookup(experiment, mapping).metadata, mapping),
            "DV": 1.0,
        }
        for (experiment, mapping), kind in KINDS.items()
        for _ in range(POINTS.get((experiment, mapping), 5))
    ]
    table = pd.DataFrame(rows)
    return FitTable(
        table=table,
        lookup=lookup,
        recorded_kinds={("DECAY", e, m): k for (e, m), k in KINDS.items()},
    )


def _credibility(
    *evidence: Evidence,
    limitations: tuple[Evidence, ...] = (),
    quantities: tuple[QuantityOfInterest, ...] = _A,
) -> Credibility:
    """Get a declaration of one context of use with the evidence."""
    return Credibility(
        contexts=(
            ContextOfUse(
                id="dose",
                title="Dose",
                question="What follows from a dose?",
                level=CouLevel.C2,
                influence=Level.LOW,
                consequence=Level.LOW,
                quantities=quantities,
                domain=Domain(),
                evidence=evidence,
                extrapolation="the decay is linear",
            ),
        ),
        limitations=limitations,
    )


def _pairs(fits: FitTable, rows: pd.Index) -> set[tuple[str, str]]:
    """Get the experiments and mappings of rows of the table."""
    selected = fits.table.loc[rows]
    return set(zip(selected.experiment, selected.mapping, strict=True))


def _n(fits: FitTable, experiment: str, mapping: str) -> int:
    """Get the number of rows of a mapping."""
    table = fits.table
    return int(((table.experiment == experiment) & (table.mapping == mapping)).sum())


def _problems(credibility: Credibility, fits: FitTable) -> str:
    """Get the message of the selection error of a declaration."""
    with pytest.raises(SelectionError) as err:
        select(credibility, fits)
    assert err.value.problems
    return str(err.value)


def test_normalize_observable() -> None:
    assert normalize_observable("[A]") == "A"
    assert normalize_observable("A") == "A"
    assert normalize_observable("pk.cmax") == "pk.cmax"


def test_a_validation_on_training_data_is_an_error(fits: FitTable) -> None:
    message = _problems(
        _credibility(
            Evidence(
                "fitted",
                EvidenceRole.VALIDATION,
                named_filter({"CredB": {"fm_dose"}}),
            )
        ),
        fits,
    )
    assert "CredB.fm_dose" in message


def test_the_key_alone_does_not_select(fits: FitTable) -> None:
    selections = select(
        _credibility(Evidence("doses", EvidenceRole.VALIDATION, named_filter(DOSE))),
        fits,
    )
    (dose,) = selections.evidence
    assert _pairs(fits, dose.rows) == {("CredA", "fm_dose")}
    assert len(dose.rows) == _n(fits, "CredA", "fm_dose")


def test_a_calibration_keeps_its_training_points(fits: FitTable) -> None:
    selections = select(
        _credibility(
            Evidence(
                "intravenous",
                EvidenceRole.CALIBRATION,
                metadata_filter(route={Route.IV}),
            )
        ),
        fits,
    )
    (calibration,) = selections.evidence
    assert _pairs(fits, calibration.rows) == {
        ("CredA", "fm_train"),
        ("CredB", "fm_dose"),
    }
    assert set(fits.table.loc[calibration.rows, "kind"]) == {"training"}


def test_outliers_are_dropped_and_counted(fits: FitTable) -> None:
    selections = select(
        _credibility(
            Evidence(
                "doses",
                EvidenceRole.VALIDATION,
                named_filter({"CredA": {"fm_dose", "fm_outlier"}}),
            )
        ),
        fits,
    )
    (dose,) = selections.evidence
    assert _pairs(fits, dose.rows) == {("CredA", "fm_dose")}
    assert dose.outliers_dropped == _n(fits, "CredA", "fm_outlier") == 3
    assert dose.limitation_subtracted == 0


def test_limitations_are_subtracted_and_counted(fits: FitTable) -> None:
    limitation = Evidence(
        "scaled", EvidenceRole.LIMITATION, named_filter({"CredA": {"fm_limit"}})
    )
    selections = select(
        _credibility(
            Evidence(
                "doses",
                EvidenceRole.VALIDATION,
                named_filter({"CredA": {"fm_dose", "fm_limit"}}),
            ),
            limitations=(limitation,),
        ),
        fits,
    )
    (dose,) = selections.evidence
    assert _pairs(fits, dose.rows) == {("CredA", "fm_dose")}
    assert dose.limitation_subtracted == _n(fits, "CredA", "fm_limit") == 2
    assert _pairs(fits, selections.limitations) == {("CredA", "fm_limit")}


def test_the_limitations_of_a_context_are_subtracted(fits: FitTable) -> None:
    credibility = _credibility(
        Evidence(
            "doses",
            EvidenceRole.VALIDATION,
            named_filter({"CredA": {"fm_dose", "fm_limit"}}),
        )
    )
    cou = dataclasses.replace(
        credibility.contexts[0],
        limitations=(
            Evidence(
                "scaled",
                EvidenceRole.LIMITATION,
                named_filter({"CredA": {"fm_limit"}}),
            ),
        ),
    )
    selections = select(Credibility(contexts=(cou,)), fits)
    (dose,) = selections.evidence
    assert _pairs(fits, dose.rows) == {("CredA", "fm_dose")}
    assert dose.limitation_subtracted == 2


def test_a_fitted_limitation_is_an_error(fits: FitTable) -> None:
    message = _problems(
        _credibility(
            Evidence(
                "doses", EvidenceRole.VALIDATION, named_filter({"CredA": {"fm_dose"}})
            ),
            limitations=(
                Evidence(
                    "fitted",
                    EvidenceRole.LIMITATION,
                    named_filter({"CredA": {"fm_train"}}),
                ),
            ),
        ),
        fits,
    )
    assert "CredA.fm_train" in message
    assert "fitted" in message


def test_unclaimed_selects_the_rest(fits: FitTable) -> None:
    selections = select(CREDIBILITY, fits)
    (unclaimed,) = (
        s
        for s in selections.evidence
        if s.context == "dose" and isinstance(s.evidence.selection, Unclaimed)
    )
    assert _pairs(fits, unclaimed.rows) == {("CredA", "fm_subject")}
    assert unclaimed.limitation_subtracted == _n(fits, "CredA", "fm_limit")

    dose = CREDIBILITY.context("dose")
    everything = dataclasses.replace(
        dose,
        evidence=(
            dose.evidence[0],
            Evidence("the rest", EvidenceRole.VALIDATION, UNCLAIMED),
        ),
    )
    credibility = dataclasses.replace(
        CREDIBILITY,
        contexts=tuple(
            everything if c.id == "dose" else c for c in CREDIBILITY.contexts
        ),
    )
    (unclaimed,) = (
        s for s in select(credibility, fits).evidence if s.evidence.name == "the rest"
    )
    assert _pairs(fits, unclaimed.rows) == {
        ("CredA", "fm_subject"),
        ("CredB", "fm_zero"),
    }


def test_the_unclaimed_validation_data(fits: FitTable) -> None:
    selections = select(CREDIBILITY, fits)
    assert isinstance(selections, Selections)
    assert selections.unclaimed.to_dict("records") == [
        {
            "fit": "DECAY",
            "experiment": "CredB",
            "mapping": "fm_zero",
            "n": _n(fits, "CredB", "fm_zero"),
        }
    ]
    assert _pairs(fits, selections.limitations) == {("CredA", "fm_limit")}


def test_quantities_match_without_brackets(fits: FitTable) -> None:
    doses = Evidence("doses", EvidenceRole.VALIDATION, named_filter(DOSE))
    for quantity in ("A", "[A]"):
        (dose,) = select(
            _credibility(
                doses, quantities=(QuantityOfInterest(quantity, "the quantity"),)
            ),
            fits,
        ).evidence
        assert _pairs(fits, dose.rows) == {("CredA", "fm_dose")}

    message = _problems(
        _credibility(doses, quantities=(*_A, QuantityOfInterest("B", "no observable"))),
        fits,
    )
    assert "'B'" in message


def test_an_unknown_fit_subset_is_an_error(fits: FitTable) -> None:
    message = _problems(
        _credibility(
            Evidence("doses", EvidenceRole.VALIDATION, named_filter(DOSE), fit="PD")
        ),
        fits,
    )
    assert "'PD'" in message
    assert "['DECAY']" in message
    assert "selects no point" not in message


def test_the_fit_of_an_evidence_restricts_it(fits: FitTable) -> None:
    other = fits.table.assign(fit="PK")
    table = pd.concat([fits.table, other], ignore_index=True)
    both = dataclasses.replace(fits, table=table)
    doses = Evidence("doses", EvidenceRole.VALIDATION, named_filter(DOSE), fit="PK")
    (dose,) = select(_credibility(doses), both).evidence
    assert set(table.loc[dose.rows, "fit"]) == {"PK"}
    assert len(dose.rows) == _n(fits, "CredA", "fm_dose")


def test_an_empty_evidence_is_an_error(fits: FitTable) -> None:
    message = _problems(
        _credibility(
            Evidence(
                "nothing",
                EvidenceRole.VALIDATION,
                named_filter({"CredA": {"fm_nothing"}}),
            )
        ),
        fits,
    )
    assert "'dose'" in message
    assert "'nothing'" in message
    assert "DECAY" in message


def test_every_problem_is_raised_at_once(fits: FitTable) -> None:
    with pytest.raises(SelectionError) as err:
        select(
            _credibility(
                Evidence(
                    "fitted",
                    EvidenceRole.VALIDATION,
                    named_filter({"CredB": {"fm_dose"}}),
                ),
                Evidence("pd", EvidenceRole.VALIDATION, named_filter(DOSE), fit="PD"),
            ),
            fits,
        )
    assert len(err.value.problems) == 2


def test_kind_differences(fits: FitTable) -> None:
    assert kind_differences(fits.table, fits.recorded_kinds) == []
    recorded = {**fits.recorded_kinds, ("DECAY", "CredA", "fm_dose"): "training"}
    (message,) = kind_differences(fits.table, recorded)
    assert "DECAY" in message
    assert "CredA" in message
    assert "fm_dose" in message

    missing = {k: v for k, v in fits.recorded_kinds.items() if k[2] != "fm_zero"} | {
        ("DECAY", "CredC", "fm_new"): "validation"
    }
    messages = kind_differences(fits.table, missing)
    assert len(messages) == 2
    assert any("CredB.fm_zero" in m for m in messages)
    assert any("CredC.fm_new" in m for m in messages)


def test_recorded_kinds(credibility_model: PKPDModel, tmp_path: Path) -> None:
    (fit_id,) = credibility_model.fits
    record = credibility_model.fits_dir / fit_id
    assert recorded_kinds(record, "DECAY") == {
        ("DECAY", e, m): k for (e, m), k in KINDS.items()
    }

    with pytest.raises(FileNotFoundError, match="adopt the fit again"):
        recorded_kinds(tmp_path, "DECAY")

    table = pd.read_csv(record / "metrics_mappings.tsv", sep="\t")
    twice = pd.concat([table, table.assign(kind="training")], ignore_index=True)
    twice.to_csv(tmp_path / "metrics_mappings.tsv", sep="\t", index=False)
    with pytest.raises(ValueError, match="fm_dose"):
        recorded_kinds(tmp_path, "DECAY")


def _edited_model(model: PKPDModel, fits_dir: Path) -> PKPDModel:
    """Get the model with a copy of its fits whose record has another kind."""
    shutil.copytree(model.fits_dir, fits_dir)
    (fit_id,) = model.fits
    path = fits_dir / fit_id / "metrics_mappings.tsv"
    table = pd.read_csv(path, sep="\t")
    edited = (table.experiment == "CredA") & (table.mapping == "fm_dose")
    table.loc[edited, "kind"] = "training"
    table.to_csv(path, sep="\t", index=False)
    return dataclasses.replace(model, fits_path=fits_dir)


def test_evaluate_fits(credibility_model: PKPDModel, tmp_path: Path) -> None:
    fits = evaluate_fits(credibility_model)
    table = fits.table
    assert {"fit", "DV", "IPRED", "DV_SD", "observable"} <= set(table.columns)
    assert set(table.fit) == {"DECAY"}
    pairs = zip(table.experiment, table.mapping, strict=True)
    kinds = dict(zip(pairs, table.kind, strict=True))
    assert kinds == KINDS
    assert kind_differences(table, fits.recorded_kinds) == []
    assert fits.lookup("CredB", "fm_zero").metadata.route is Route.PO  # ty: ignore[unresolved-attribute]
    assert fits.lookup("CredA", "fm_renal").metadata.health is Health.RENAL_IMPAIRMENT  # ty: ignore[unresolved-attribute]

    # the static table has the same mappings with the same observables
    static = mapping_table(credibility_model)
    columns = ["fit", "experiment", "mapping", "kind", "observable", "route", "health"]
    assert sorted(map(tuple, static.table[columns].to_numpy().tolist())) == sorted(
        map(tuple, table[columns].drop_duplicates().to_numpy().tolist())
    )
    assert static.recorded_kinds == fits.recorded_kinds

    edited = _edited_model(credibility_model, tmp_path / "fits")
    with pytest.raises(KindMismatchError, match="adopt") as err:
        evaluate_fits(edited)
    (problem,) = err.value.problems
    assert "CredA.fm_dose" in problem


def test_check_declaration_of_the_test_package(
    credibility_model: PKPDModel, tmp_path: Path
) -> None:
    assert check_declaration(credibility_model) == []

    broken = dataclasses.replace(credibility_model, credibility=".credibility:BROKEN")
    messages = check_declaration(broken)
    assert len(messages) == 1
    assert "CredB" in messages[0]

    edited = _edited_model(broken, tmp_path / "fits")
    messages = check_declaration(edited)
    assert any("CredA.fm_dose" in m and "training" in m for m in messages)
    assert any("CredB.fm_dose" in m for m in messages)

    missing = tmp_path / "missing"
    shutil.copytree(credibility_model.fits_dir, missing)
    (fit_id,) = credibility_model.fits
    (missing / fit_id / "metrics_mappings.tsv").unlink()
    (message,) = check_declaration(
        dataclasses.replace(credibility_model, fits_path=missing)
    )
    assert "adopt the fit again" in message
