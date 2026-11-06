"""The equivalency of the inputs (`pkpdsim.credibility.inputs`).

The functions on the tables of a study and the record of a context run
without a simulation; `input_checks` runs on the adopted fit of the test
package with the volume `Vc` as the covariate, reported by the group of the
study tables which the tests write.
"""

import json
import math
from pathlib import Path

import pandas as pd
import pytest
from credibility_tables import fit_table
from pkpdsim_credibility.credibility import CREDIBILITY
from sbmlsim import Q
from sbmlsim.uncertainty import Uncertain
from sbmlsim.units import ureg

from pkpdsim.credibility import (
    DEFAULT_COVARIATES,
    Credibility,
    CredibilityAssessment,
    DeclarationError,
    InputStatus,
    Provenance,
    assess_table,
    evaluate_fits,
    input_checks,
)
from pkpdsim.credibility.inputs import (
    InputCheck,
    Reported,
    StudyTables,
    _magnitude,
    check_covariate,
    check_dose,
    inputs_record,
    inputs_statement,
    read_study_tables,
    reported_doses,
    reported_subjects,
)
from pkpdsim.credibility.report import (
    input_ranges_table,
    inputs_table,
    write_outputs,
)
from pkpdsim.credibility.scorecard import scorecard
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

GROUPS = pd.DataFrame(
    [
        # all: weight with its range, height in m
        ("all", None, "weight", 70.0, None, None, None, None, 55.0, 90.0, None, "kg"),
        ("all", None, "height", 1.75, None, None, None, None, None, None, None, "m"),
        # A: weight with an SD, a child of all
        ("A", "all", "weight", 60.0, None, 5.0, None, 10, None, None, None, "kg"),
        # B: weight with an SE, its height from all
        ("B", "all", "weight", 80.0, None, None, 2.0, 25, None, None, None, "kg"),
        # C: a geometric mean
        ("C", "all", "weight", None, None, None, None, 9, None, None, 1.1, "kg"),
    ],
    columns=[
        "name",
        "parent",
        "measurement_type",
        "mean",
        "median",
        "sd",
        "se",
        "count",
        "min",
        "max",
        "gsd",
        "unit",
    ],
).assign(gmean=[None, None, None, None, 75.0])

INDIVIDUALS = pd.DataFrame(
    [
        ("S1", "A", "weight", 64.0, "kg"),
        ("S2", "B", "height", 180.0, "cm"),
    ],
    columns=["name", "group", "measurement_type", "value", "unit"],
)

INTERVENTIONS = pd.DataFrame(
    [
        ("D25", "dosing", 25.0, "mg"),
        ("D50", "dosing", 50.0, "mg"),
        ("K1", "dosing", 1.0, "mg/kg"),
        ("Q", "qualitative dosing", None, None),
    ],
    columns=["name", "measurement_type", "value", "unit"],
)

TABLES = StudyTables(
    groups=GROUPS, individuals=INDIVIDUALS, interventions=INTERVENTIONS
)


def test_the_reported_values_of_the_groups() -> None:
    (whole,) = reported_subjects(TABLES, ["all"], [], "weight", "kg", ureg)
    assert whole == Reported(by="group all", central=70.0, low=55.0, high=90.0)
    # the mean with two SD
    (a,) = reported_subjects(TABLES, ["A"], [], "weight", "kg", ureg)
    assert (a.central, a.low, a.high) == (60.0, 50.0, 70.0)
    # the SD of the SE and the count
    (b,) = reported_subjects(TABLES, ["B"], [], "weight", "kg", ureg)
    assert b.low == pytest.approx(60.0) and b.high == pytest.approx(100.0)
    # the geometric mean with the GSD
    (c,) = reported_subjects(TABLES, ["C"], [], "weight", "kg", ureg)
    assert c.central == 75.0
    assert c.low == pytest.approx(75.0 / 1.21) and c.high == pytest.approx(75.0 * 1.21)
    # a group without the measurement reads its parent, in the unit of the model
    (height,) = reported_subjects(TABLES, ["B"], [], "height", "cm", ureg)
    assert height.by == "group all"
    assert height.central == pytest.approx(175.0)
    # an unknown group and a measurement without a value report nothing
    assert reported_subjects(TABLES, ["X"], [], "weight", "kg", ureg) == []
    assert reported_subjects(TABLES, ["all"], [], "age", "yr", ureg) == []


def test_the_reported_values_of_the_individuals() -> None:
    (s1,) = reported_subjects(TABLES, [], ["S1"], "weight", "kg", ureg)
    assert s1 == Reported(by="individual S1", central=64.0)
    # an individual without the measurement reads its group
    (s2,) = reported_subjects(TABLES, [], ["S2"], "weight", "kg", ureg)
    assert s2.by == "group B"
    (h2,) = reported_subjects(TABLES, [], ["S2"], "height", "m", ureg)
    assert h2.by == "individual S2"
    assert h2.central == pytest.approx(1.8)
    # each subject once
    assert len(reported_subjects(TABLES, ["B"], ["S2"], "weight", "kg", ureg)) == 1


def test_a_unit_which_does_not_convert() -> None:
    with pytest.raises(ValueError, match="cannot be converted"):
        reported_subjects(TABLES, ["all"], [], "weight", "l", ureg)


def test_the_check_of_a_covariate() -> None:
    reported = [Reported(by="group all", central=70.0, low=55.0, high=90.0)]
    assert check_covariate(75.0, reported) == (InputStatus.EQUIVALENT, reported[0])
    assert check_covariate(88.0, reported) == (InputStatus.IN_RANGE, reported[0])
    assert check_covariate(95.0, reported) == (InputStatus.DIFFERS, reported[0])
    assert check_covariate(75.0, []) == (InputStatus.NOT_REPORTED, None)
    # the first equivalent one decides, else the nearest
    two = [Reported(by="a", central=50.0), Reported(by="b", central=100.0)]
    assert check_covariate(98.0, two) == (InputStatus.EQUIVALENT, two[1])
    assert check_covariate(75.0, two)[0] is InputStatus.DIFFERS
    assert check_covariate(60.0, two)[1] == two[0]
    # a value without a central value and a range is not comparable
    assert check_covariate(1.0, [Reported(by="a", central=None)])[0] is (
        InputStatus.NOT_REPORTED
    )


def test_the_doses() -> None:
    doses = reported_doses(TABLES, ["D25", "D50", "Q"], "mg", ureg)
    assert [(d.by, d.central) for d in doses] == [
        ("intervention D25", 25.0),
        ("intervention D50", 50.0),
    ]
    # a dose per body weight needs the body weight
    assert reported_doses(TABLES, ["K1"], "mg", ureg) == []
    (per_kg,) = reported_doses(TABLES, ["K1"], "mg", ureg, Q(70, "kg"))
    assert per_kg.central == pytest.approx(70.0)
    assert check_dose(25.4, doses) == (InputStatus.EQUIVALENT, doses[0])
    status, nearest = check_dose(45.8, doses)
    assert status is InputStatus.DIFFERS and nearest == doses[1]
    assert check_dose(25.0, []) == (InputStatus.NOT_REPORTED, None)
    assert (
        reported_doses(StudyTables(*(pd.DataFrame(),) * 3), ["D25"], "mg", ureg) == []
    )


def test_the_values_of_the_changes() -> None:
    assert _magnitude(Q(1, "g"), "mg") == pytest.approx(1000.0)
    assert _magnitude(25, "mg") == 25.0
    assert _magnitude(Uncertain(Q(5, "mg"), cv=0.1), "mg") == pytest.approx(5.0)
    # a dose added to its state, see add_doses
    assert _magnitude("PODOSE_a + 12.5", "mg") == 12.5
    assert _magnitude("-na_urine", "mg") is None
    assert _magnitude(Q(1, "l"), "mg") is None
    assert _magnitude(math.nan, "mg") is None


def test_the_tables_of_a_study(tmp_path: Path) -> None:
    folder = tmp_path / "S"
    folder.mkdir()
    GROUPS.to_csv(folder / ".S_TabGroups.tsv", sep="\t", index=False, na_rep="NA")
    tables = read_study_tables(folder, "S")
    assert len(tables.groups) == len(GROUPS)
    assert tables.individuals.empty and tables.interventions.empty
    assert read_study_tables(None, "S").groups.empty


def _check(mapping: str, target: str, value: float, status: str) -> InputCheck:
    """Get a check of a mapping of CredA of the fit DECAY."""
    return InputCheck(
        fit="DECAY",
        experiment="CredA",
        mapping=mapping,
        target=target,
        kind="dose" if target.startswith("PODOSE") else "covariate",
        measurement="dosing" if target.startswith("PODOSE") else "weight",
        unit="mg" if target.startswith("PODOSE") else "kg",
        simulated=value,
        source="simulation",
        reported=70.0 if status != "not_reported" else None,
        low=None,
        high=None,
        reported_by=["group all"] if status != "not_reported" else [],
        status=status,
        ratio=value / 70.0 if status != "not_reported" else None,
    )


def test_the_record_of_a_context() -> None:
    checks = [
        _check("fm_a", "BW", 70.0, "equivalent"),
        _check("fm_b", "BW", 90.0, "in_range"),
        _check("fm_c", "BW", 120.0, "differs"),
        _check("fm_c", "PODOSE_a", 25.0, "not_reported"),
        # not decisive
        _check("fm_x", "BW", 200.0, "differs"),
    ]
    decisive = pd.DataFrame(
        {
            "fit": "DECAY",
            "experiment": "CredA",
            "mapping": ["fm_a", "fm_a", "fm_b", "fm_c", "fm_d"],
        }
    )
    record = inputs_record(checks, decisive)
    assert record["n_mappings"] == 4
    assert record["n_checked"] == 3
    bw, dose = record["inputs"]
    assert (bw["target"], bw["low"], bw["high"], bw["n_mappings"]) == ("BW", 70, 120, 3)
    assert (bw["n_equivalent"], bw["n_in_range"], bw["n_differs"]) == (1, 1, 1)
    assert (dose["kind"], dose["n_not_reported"]) == ("dose", 1)
    assert [c["mapping"] for c in record["differs"]] == ["fm_c"]
    statement = inputs_statement(record)
    assert statement.startswith("The simulated inputs of 3 of the 4 decisive mappings")
    assert (
        "BW (1 equivalent, 1 in the reported range, 1 differs; simulated 70 to 120 kg)"
        in statement
    )
    assert "Differ: CredA/fm_c BW 120 kg against 70 kg (group all)." in statement
    assert inputs_record(checks, decisive.iloc[0:0])["n_mappings"] == 0


def test_the_covariates_of_a_declaration() -> None:
    assert Credibility(contexts=()).covariates == DEFAULT_COVARIATES
    declared = Credibility(contexts=(), covariates={"BW_kg": "weight"})
    assert dict(declared.covariates) == {"BW_kg": "weight"}
    with pytest.raises(DeclarationError, match="covariate 'BW'"):
        Credibility(contexts=(), covariates={"BW": " "})


def _data(data: Path) -> Path:
    """Get a data folder with the study tables of CredA."""
    folder = data / "CredA"
    folder.mkdir(parents=True)
    pd.DataFrame(
        [("all", "NA", "volume", 1050.0, "ml"), ("all", "NA", "rate", 0.2, "1/min")],
        columns=["name", "parent", "measurement_type", "mean", "unit"],
    ).to_csv(folder / ".CredA_TabGroups.tsv", sep="\t", index=False)
    return data


def test_the_checks_of_the_adopted_fit(
    credibility_model: PKPDModel, tmp_path: Path
) -> None:
    model = credibility_model.replace(data_path=(_data(tmp_path / "data"),))
    fits = evaluate_fits(model)
    # the toy has no body weight and no doses
    assert input_checks(model, fits, DEFAULT_COVARIATES) == []
    # the fitted rate is no input, a target the model has not is left out
    checks = input_checks(model, fits, {"Vc": "volume", "k": "rate", "X": "weight"})
    by_experiment = {(c["experiment"], c["mapping"]): c for c in checks}
    keys = set(fits.table[["experiment", "mapping"]].itertuples(index=False, name=None))
    assert set(by_experiment) == keys
    check = by_experiment["CredA", "fm_dose"]
    assert check["target"] == "Vc" and check["source"] == "model"
    assert (check["simulated"], check["unit"]) == (1.0, "l")
    assert check["status"] == "equivalent"
    assert check["reported"] == pytest.approx(1.05)
    assert check["reported_by"] == ["group all"]
    # CredB has no study tables
    assert by_experiment["CredB", "fm_dose"]["status"] == "not_reported"


def test_the_record_and_the_outputs(tmp_path: Path) -> None:
    checks = [
        _check("fm_dose", "BW", 72.0, "equivalent"),
        _check("fm_dose", "PODOSE_a", 25.0, "differs"),
        _check("fm_renal", "BW", 60.0, "in_range"),
    ]
    assessment = assess_table(CREDIBILITY, fit_table(), PROVENANCE, inputs=checks)
    contexts = {c.id: c for c in assessment.contexts}
    dose = contexts["dose"]
    assert dose.inputs is not None
    assert dose.inputs["n_checked"] == 1
    assert [s["target"] for s in dose.inputs["inputs"]] == ["BW", "PODOSE_a"]
    factor = next(f for f in dose.factors if f["id"] == "equivalency_of_inputs")
    assert factor["status"] == "assessed"
    assert "Differ: CredA/fm_dose PODOSE_a 25 mg against 70 mg" in factor["statement"]
    # a context without a checked decisive mapping keeps the factor not assessed
    describe = contexts["describe"]
    assert describe.inputs is not None and describe.inputs["n_checked"] == 0
    factor = next(f for f in describe.factors if f["id"] == "equivalency_of_inputs")
    assert factor["status"] == "not_assessed"
    path = tmp_path / "credibility.json"
    assessment.to_json(path)
    read = CredibilityAssessment.from_json(path)
    assert read.to_dict() == assessment.to_dict()
    assert read.inputs == checks
    # a record of pkpdsim 0.18 has no inputs
    data = json.loads(path.read_text(encoding="utf-8"))
    del data["inputs"]
    for context in data["contexts"]:
        del context["inputs"]
    old = CredibilityAssessment.from_dict(data)
    assert old.inputs == [] and old.contexts[0].inputs is None
    # without the checks no context has a record
    plain = assess_table(CREDIBILITY, fit_table(), PROVENANCE)
    assert all(c.inputs is None for c in plain.contexts)
    assert inputs_table(assessment).loc[1, "reported_by"] == "group all"
    ranges = input_ranges_table(assessment)
    assert list(ranges.loc[ranges["context"] == "dose", "target"]) == ["BW", "PODOSE_a"]
    assert (
        scorecard({"credibility": assessment})
        .set_index("context")
        .loc["dose", "inputs_differ"]
        == 1
    )
    output = write_outputs(assessment, tmp_path / "out")
    assert (
        (output / "inputs.tsv")
        .read_text(encoding="utf-8")
        .startswith("fit\texperiment\tmapping")
    )
    assert (
        (output / "input_ranges.tsv")
        .read_text(encoding="utf-8")
        .startswith("context\ttarget")
    )
    html = (output / "index.html").read_text(encoding="utf-8")
    assert "<h4>Inputs</h4>" in html
    assert "CredA/fm_dose" in html
