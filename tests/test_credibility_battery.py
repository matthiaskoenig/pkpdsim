"""The battery of the checks of the builds of a model (`pkpdsim.credibility.battery`).

The builds are small drug models in antimony which follow the conventions of
the PBPK models: a drug `d` with its doses and a metabolite `m` excreted into
the urine; every test breaks one property and expects its check to fail.
"""

import dataclasses
import importlib.util
import math
from pathlib import Path

import antimony
import pytest
from credibility_tables import fit_table, verification_of
from pkpdsim_credibility import MODEL as CREDIBILITY_MODEL
from pkpdsim_credibility.credibility import CREDIBILITY

from pkpdsim.credibility import (
    CredibilityAssessment,
    FactorStatus,
    Provenance,
    Verdict,
    assess_table,
)
from pkpdsim.credibility.battery import (
    DOSE,
    Battery,
    BatteryStatus,
    Fate,
    FlatBuild,
    run_battery,
    same_sbml,
)
from pkpdsim.credibility.report import battery_table, fates_table, write_outputs
from pkpdsim.credibility.verification import BATTERY_CHECKS, Check
from pkpdsim.model import CheckSettings, PKPDModel
from pkpdsim.testing import check_battery
from pkpdsim_toy import MODEL as TOY

#: a drug model of the conventions: the doses enter the venous blood, the drug
#: is metabolized in the liver and eliminated, the metabolite is excreted
DRUG = """
model drug()
  compartment Vve = 5, Vli = 1.5;
  species Cve_d in Vve = 0;
  species Cli_m in Vli = 0;
  substanceOnly species Aurine_m in Vve = 0;
  Mr_d = 300; Mr_m = 316;
  IVDOSE_d = 0; Ri_d = 0; PODOSE_d = 0;
  ka = 0.1; kabs = 0.02; k_met = 0.01; k_ex = 0.005; CL = 0.002;
  IVDOSE_d' = -ka*IVDOSE_d + Ri_d;
  PODOSE_d' = -kabs*PODOSE_d;
  iv_d: => Cve_d; ka*IVDOSE_d/Mr_d;
  po_d: => Cve_d; kabs*PODOSE_d/Mr_d;
  MET: Cve_d => Cli_m; k_met*Vve*Cve_d;
  EX: Cli_m => Aurine_m; k_ex*Vli*Cli_m;
  ELIM: Cve_d => ; CL*Cve_d;
  {extra}
end
"""

PROVENANCE = Provenance(
    model="credibility",
    title="Toy model with every kind of evidence",
    fits=["PK"],
    versions={"pkpdsim": "0.0.0"},
    commit=None,
    dirty=None,
    date="2026-11-03T12:00:00+01:00",
)


def _sbml(path: Path, text: str) -> Path:
    """Write the SBML of an antimony model."""
    antimony.clearPreviousLoads()
    if antimony.loadAntimonyString(text) < 0:
        raise ValueError(antimony.getLastError())
    path.write_text(antimony.getSBMLString("drug"), encoding="utf-8")
    return path


def _model(
    tmp_path: Path, extra: str = "", *, replace: tuple[str, str] = ("", "")
) -> PKPDModel:
    """Get a model of the drug model with an extra line or a replaced text."""
    text = DRUG.replace("{extra}", extra)
    if replace[0]:
        assert replace[0] in text
        text = text.replace(*replace)
    path = _sbml(tmp_path / "drug.xml", text)
    return PKPDModel(
        name="drug",
        package="drug",
        title="drug",
        path=tmp_path,
        model_path=path,
        factory="drug:create_drug_models",
    )


def _battery(model: PKPDModel) -> Battery:
    """Run the battery without the rebuild, the drug models have no factory."""
    return run_battery(model, rebuild=False, cross_simulator=False)


def test_a_sound_model_passes(tmp_path: Path) -> None:
    battery = _battery(_model(tmp_path))
    assert [r.check for r in battery.results] == list(BATTERY_CHECKS)
    status = {r.check: r.status for r in battery.results}
    for check in (
        Check.MASS_BALANCE,
        Check.NON_NEGATIVITY,
        Check.ZERO_DOSE,
        Check.STEADY_STATE,
        Check.SUPERPOSITION,
        Check.PLAUSIBILITY,
    ):
        assert status[check] is BatteryStatus.PASSED, battery.result(check)
    assert status[Check.CROSS_SIMULATOR] is BatteryStatus.NOT_APPLICABLE
    assert status[Check.REBUILD] is BatteryStatus.NOT_APPLICABLE
    mass = battery.result(Check.MASS_BALANCE)
    # the run without a dose and one per dose
    assert mass.cases == 4
    assert mass.largest < 1e-10


def test_the_fate_of_a_dose(tmp_path: Path) -> None:
    battery = _battery(_model(tmp_path))
    fates = {f.dose: f for f in battery.fates}
    assert set(fates) == {"IVDOSE_d", "Ri_d", "PODOSE_d"}
    fate = fates["IVDOSE_d"]
    assert fate.substance == "d"
    assert fate.amount == pytest.approx(DOSE / 300)
    parts = fate.parts
    # the metabolism moves the drug into the metabolite: no part of its own
    assert "MET" not in parts
    assert parts["ELIM"] < 0
    assert parts["urine"] > 0
    left = parts["body"] + parts["urine"] + parts.get("dose", 0.0)
    assert left - parts["ELIM"] == pytest.approx(fate.amount, rel=1e-8)


def test_a_dose_with_the_wrong_molecular_weight_breaks_the_mass_balance(
    tmp_path: Path,
) -> None:
    model = _model(
        tmp_path,
        replace=(
            "iv_d: => Cve_d; ka*IVDOSE_d/Mr_d",
            "iv_d: => Cve_d; ka*IVDOSE_d/(2*Mr_d)",
        ),
    )
    result = _battery(model).result(Check.MASS_BALANCE)
    assert result.status is BatteryStatus.FAILED
    assert any("male IVDOSE_d" in p for p in result.problems)
    assert any("male Ri_d" in p for p in result.problems)
    assert not any("PODOSE_d" in p for p in result.problems)


def test_a_rule_on_a_species_breaks_the_mass_balance(tmp_path: Path) -> None:
    model = _model(tmp_path, "species Cbo_d in Vve = 0; Cbo_d' = 0.001*Cve_d;")
    result = _battery(model).result(Check.MASS_BALANCE)
    assert result.status is BatteryStatus.FAILED
    assert (
        "male: the species 'Cbo_d' is changed by a rule or an event"
        in result.problems[0]
    )


def test_a_production_of_the_drug_fails_the_zero_dose_and_the_steady_state(
    tmp_path: Path,
) -> None:
    battery = _battery(_model(tmp_path, "PROD: => Cve_d; 1e-6;"))
    zero = battery.result(Check.ZERO_DOSE)
    assert zero.status is BatteryStatus.FAILED
    assert any("'Cve_d' of the substance 'd'" in p for p in zero.problems)
    assert battery.result(Check.STEADY_STATE).status is BatteryStatus.FAILED
    # the production is a flow of the balance, which still closes
    assert battery.result(Check.MASS_BALANCE).status is BatteryStatus.PASSED


def test_a_negative_species_fails_the_non_negativity(tmp_path: Path) -> None:
    battery = _battery(_model(tmp_path, "LOSS: Cli_m => ; 1e-6;"))
    result = battery.result(Check.NON_NEGATIVITY)
    assert result.status is BatteryStatus.FAILED
    assert any("'Cli_m' reaches" in p for p in result.problems)


def test_a_saturated_elimination_fails_the_superposition(tmp_path: Path) -> None:
    model = _model(
        tmp_path,
        "Vmax = 1e-6; Km = 1e-9;",
        replace=(
            "ELIM: Cve_d => ; CL*Cve_d",
            "ELIM: Cve_d => ; Vmax*Cve_d/(Km + Cve_d)",
        ),
    )
    result = _battery(model).result(Check.SUPERPOSITION)
    assert result.status is BatteryStatus.FAILED
    assert "to twice the dose differs from the double" in result.problems[0]


@pytest.mark.parametrize(
    ("k", "status"), [(0.002, BatteryStatus.PASSED), (0.003, BatteryStatus.FAILED)]
)
def test_an_endogenous_substance_at_its_steady_state(
    tmp_path: Path, k: float, status: BatteryStatus
) -> None:
    """The production 0.01 mmole/min equals the elimination at 1 mM only for k 0.002."""
    extra = (
        f"species Cve_e in Vve = 1; Mr_e = 100; ke = {k};"
        " PRODE: => Cve_e; 0.01; ELIME: Cve_e => ; ke*Vve*Cve_e;"
    )
    battery = _battery(_model(tmp_path, extra))
    steady = battery.result(Check.STEADY_STATE)
    assert steady.status is status
    if status is BatteryStatus.FAILED:
        assert "'Cve_e' changes by" in steady.problems[0]
    # the drug still starts at zero, the endogenous substance is no case
    assert battery.result(Check.ZERO_DOSE).cases == 3


def test_a_model_without_a_substance(tmp_path: Path) -> None:
    battery = run_battery(TOY, rebuild=False, cross_simulator=False)
    for check in (
        Check.MASS_BALANCE,
        Check.NON_NEGATIVITY,
        Check.ZERO_DOSE,
        Check.STEADY_STATE,
        Check.SUPERPOSITION,
        Check.CROSS_SIMULATOR,
    ):
        result = battery.result(check)
        assert result.status is BatteryStatus.NOT_APPLICABLE
        assert "Mr_<s>" in result.statement
    assert battery.result(Check.PLAUSIBILITY).status is BatteryStatus.PASSED
    assert battery.fates == ()


def test_a_compartment_without_a_volume_is_implausible(tmp_path: Path) -> None:
    result = _battery(_model(tmp_path, "compartment Vx = 0;")).result(
        Check.PLAUSIBILITY
    )
    assert result.status is BatteryStatus.FAILED
    assert result.problems == ("male: the compartment 'Vx' has the size 0",)


def test_a_fitted_value_at_its_bound_is_implausible(
    credibility_model: PKPDModel,
) -> None:
    records = credibility_model.fit_records()
    pid, parameter = next(iter(records.parameters.items()))
    at_bound = dataclasses.replace(parameter, value=parameter.upper_bound)
    records.parameters[pid] = at_bound

    class Model:
        fits = credibility_model.fits

        @staticmethod
        def fit_records() -> object:
            return records

    from pkpdsim.credibility.battery import _Collector, _plausibility

    collector = _Collector()
    _plausibility(Model(), [], collector)  # ty: ignore[invalid-argument-type]
    assert collector.problems == [
        f"the value {parameter.upper_bound:g} of '{pid}' of the fit "
        f"'{parameter.fit_id}' lies at its upper bound {parameter.upper_bound:g}"
    ]


def test_the_rebuild(tmp_path: Path) -> None:
    passed = run_battery(CREDIBILITY_MODEL, cross_simulator=False).result(Check.REBUILD)
    assert passed.status is BatteryStatus.PASSED
    assert passed.cases == 1
    copy = tmp_path / "models" / CREDIBILITY_MODEL.model_path.name
    copy.parent.mkdir()
    copy.write_text(
        CREDIBILITY_MODEL.model_path.read_text(encoding="utf-8").replace(
            'name="rate of the decay"', 'name="the rate of the decay"'
        ),
        encoding="utf-8",
    )
    model = dataclasses.replace(CREDIBILITY_MODEL, model_path=copy)
    failed = run_battery(model, cross_simulator=False).result(Check.REBUILD)
    assert failed.status is BatteryStatus.FAILED
    assert failed.problems == ("'toy.xml' differs from its rebuild",)


def test_a_factory_which_fails(tmp_path: Path) -> None:
    result = run_battery(_model(tmp_path), cross_simulator=False).result(Check.REBUILD)
    assert result.status is BatteryStatus.FAILED
    assert result.problems[0].startswith("the factory failed: ")


def test_same_sbml_ignores_the_dates() -> None:
    a = "<a>\n<dcterms:W3CDTF>2026-11-03T12:00:00Z</dcterms:W3CDTF>\n<b/>\n</a>"
    b = "<a>\n<dcterms:W3CDTF>2026-11-04T12:00:00Z</dcterms:W3CDTF>\n<b/>\n</a>"
    assert same_sbml(a, b)
    assert not same_sbml(a, b.replace("<b/>", "<c/>"))


def test_a_dose_through_the_stomach(tmp_path: Path) -> None:
    """A state of the dose in mg which a reaction reads, and a sum of the doses."""
    model = _model(
        tmp_path,
        "unit mg = 1e-3 gram; STOMACH_d = 0; cum_dose_d = 0; kst = 0.05;"
        " STOMACH_d has mg; cum_dose_d has mg;"
        " STOMACH_d' = kabs*PODOSE_d - kst*STOMACH_d; cum_dose_d' = Ri_d;",
        replace=(
            "po_d: => Cve_d; kabs*PODOSE_d/Mr_d",
            "po_d: => Cve_d; kst*STOMACH_d/Mr_d",
        ),
    )
    build = FlatBuild.read("male", model.model_path)
    assert dict(build.dose_states) == {
        "IVDOSE_d": "d",
        "PODOSE_d": "d",
        "STOMACH_d": "d",
    }
    assert build.inputs == {"iv_d", "po_d"}
    assert _battery(model).result(Check.MASS_BALANCE).status is BatteryStatus.PASSED


def test_a_species_by_its_prefix(tmp_path: Path) -> None:
    model = _model(
        tmp_path,
        "species d_lumen in Vve = 0; ABS: d_lumen => Cve_d; 0.01*d_lumen*Vve;",
    )
    build = FlatBuild.read("male", model.model_path)
    assert build.tracked["d_lumen"] == "d"
    assert build.tracked["Cve_d"] == "d"


def test_the_battery_round_trip(tmp_path: Path) -> None:
    battery = _battery(_model(tmp_path))
    read = Battery.from_dict(battery.to_dict())
    # equal but for NaN, which is not equal to itself
    assert read.to_dict() == battery.to_dict()
    assert read.fates == battery.fates
    assert math.isnan(read.result(Check.REBUILD).largest)


def _with_battery(battery: Battery) -> CredibilityAssessment:
    fits = fit_table()
    verification = verification_of(fits)
    summary = dataclasses.replace(verification.summary, battery=battery)
    return assess_table(
        CREDIBILITY,
        fits,
        PROVENANCE,
        dataclasses.replace(verification, summary=summary),
    )


def _failed(battery: Battery, check: Check) -> Battery:
    return dataclasses.replace(
        battery,
        results=tuple(
            dataclasses.replace(r, status=BatteryStatus.FAILED, problems=("broken",))
            if r.check is check
            else r
            for r in battery.results
        ),
    )


def test_the_battery_in_an_assessment(tmp_path: Path) -> None:
    battery = _battery(_model(tmp_path))
    assessment = _with_battery(battery)
    applicable = [r.check.value for r in battery.results if r.passed is not None]
    for context in assessment.contexts:
        checks = [r["check"] for r in context.verification]
        assert checks[3:] == applicable
        factors = {f["id"]: f for f in context.factors}
        code = factors["numerical_code_verification"]
        assert code["status"] == FactorStatus.ASSESSED.value
        assert code["statement"].startswith("sbmlsim runs the SBML Test Suite")
        assert "Mass balance of 4 cases" in code["statement"]
        # the rebuild did not run, its factor keeps the platform statement
        assert factors["software_quality_assurance"]["status"] == "platform"
    dose = next(c for c in assessment.contexts if c.id == "dose")
    assert dose.verdict is Verdict.SUPPORTED
    broken = _with_battery(_failed(battery, Check.STEADY_STATE))
    for context in broken.contexts:
        if context.n_decisive:
            assert context.verdict is Verdict.NOT_VERIFIED, context.id
    record = tmp_path / "credibility.json"
    broken.to_json(record)
    read = CredibilityAssessment.from_json(record).verification
    assert read is not None
    assert broken.verification is not None
    assert read.to_dict() == broken.verification.to_dict()


def test_the_outputs_of_the_battery(tmp_path: Path) -> None:
    battery = _battery(_model(tmp_path))
    assessment = _with_battery(battery)
    table = battery_table(assessment)
    assert list(table["check"]) == [c.value for c in BATTERY_CHECKS]
    assert set(table["status"]) == {"passed", "not_applicable"}
    fates = fates_table(assessment)
    assert set(fates["dose"]) == {"IVDOSE_d", "Ri_d", "PODOSE_d"}
    iv = fates[fates["dose"] == "IVDOSE_d"].set_index("part")
    assert iv.loc["ELIM", "share"] < 0
    output = write_outputs(assessment, tmp_path / "credibility")
    assert (
        (output / "battery.tsv")
        .read_text(encoding="utf-8")
        .startswith("check\tfactor\tstatus")
    )
    html = (output / "index.html").read_text(encoding="utf-8")
    assert "Battery of the builds" in html
    assert "Fate of the doses" in html
    assert "battery of the builds</a>, 6 checks" in html


def test_check_battery(tmp_path: Path) -> None:
    model = _model(tmp_path, "PROD: => Cve_d; 1e-6;")
    with pytest.raises(AssertionError, match="the check 'zero_dose' failed"):
        check_battery(model, rebuild=False, cross_simulator=False)
    known = dataclasses.replace(
        model,
        checks=CheckSettings(
            battery={"zero_dose": "a production", "steady_state": "a production"}
        ),
    )
    check_battery(known, rebuild=False, cross_simulator=False)
    stale = dataclasses.replace(
        _model(tmp_path), checks=CheckSettings(battery={"zero_dose": "fixed since"})
    )
    with pytest.raises(
        AssertionError, match="'zero_dose' no longer applies, it passed"
    ):
        check_battery(stale, rebuild=False, cross_simulator=False)


def test_the_known_failures_of_the_battery_are_checked() -> None:
    with pytest.raises(ValueError, match="'unknown' is no check of the battery"):
        CheckSettings(battery={"unknown": "reason"})
    with pytest.raises(ValueError, match="'zero_dose' needs a reason"):
        CheckSettings(battery={"zero_dose": " "})
    assert hash(CheckSettings(battery={"zero_dose": "x"})) == hash(
        CheckSettings(battery={"zero_dose": "x"})
    )


def test_without_amici_the_second_simulator_is_not_applicable(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import sys

    monkeypatch.setitem(sys.modules, "amici.sim.sundials", None)
    battery = run_battery(_model(tmp_path), rebuild=False)
    result = battery.result(Check.CROSS_SIMULATOR)
    assert result.status is BatteryStatus.NOT_APPLICABLE
    assert "AMICI is not installed" in result.statement


@pytest.mark.skipif(
    importlib.util.find_spec("amici") is None, reason="AMICI is not installed"
)
def test_the_second_simulator(tmp_path: Path) -> None:
    result = run_battery(_model(tmp_path), rebuild=False).result(Check.CROSS_SIMULATOR)
    assert result.status is BatteryStatus.PASSED, result
    # without a dose and one per dose
    assert result.cases == 4
    assert result.largest < 1e-4


def test_a_fate_round_trip() -> None:
    fate = Fate(
        build="male",
        dose="IVDOSE_d",
        substance="d",
        amount=0.1,
        parts={"body": 0.05, "ELIM": -0.05},
        imbalance=1e-15,
    )
    assert Fate.from_dict(fate.to_dict()) == fate


@pytest.mark.skipif(
    importlib.util.find_spec("amici") is None, reason="AMICI is not installed"
)
def test_a_compilation_which_takes_too_long(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from pkpdsim.credibility import battery

    monkeypatch.setattr(battery.tempfile, "gettempdir", lambda: str(tmp_path / "tmp"))
    monkeypatch.setattr(battery, "AMICI_TIMEOUT", 0.1)
    build = FlatBuild.read("male", _model(tmp_path).model_path)
    with pytest.raises(battery._CompileTimeout):
        battery._amici_model(build)
    model = _model(tmp_path)
    result = run_battery(model, rebuild=False).result(Check.CROSS_SIMULATOR)
    assert result.status is BatteryStatus.NOT_APPLICABLE
    assert (
        result.statement
        == "AMICI did not compile the male build within 0.00166667 min."
    )
