"""The PK battery of the checks of the builds of a model (`pkpdsim.credibility.battery`).

The PK battery is the battery hook of every `PKPDModel` in the battery of
sbmlsim (`sbmlsim.credibility.battery.run_battery`). The builds are small drug
models in antimony which follow the conventions of the PBPK models: a drug `d`
with its doses and a metabolite `m` excreted into the urine; every test breaks
one property and expects its check to fail.
"""

import dataclasses
import importlib.util
import math
from pathlib import Path

import antimony
import pytest
from credibility_tables import assess, fit_table, verification_of
from pkpdsim_credibility import MODEL as CREDIBILITY_MODEL
from pkpdsim_credibility.credibility import CREDIBILITY
from sbmlsim.credibility import (
    CHECKS,
    Battery,
    BatteryStatus,
    CredibilityAssessment,
    FactorStatus,
    Verdict,
)
from sbmlsim.credibility import battery as core
from sbmlsim.credibility.battery import run_battery
from sbmlsim.credibility.report import battery_table, write_outputs
from sbmlsim.project import CheckSettings, Project
from sbmlsim.project_checks import check_battery

from pkpdsim.credibility import EXTENSIONS, PK_BATTERY, PK_CHECKS
from pkpdsim.credibility.battery import (
    DOSE,
    HORIZON,
    STEPS,
    Fate,
    FlatBuild,
    dose_cases,
    fates_table,
)
from pkpdsim.model import PKPDModel
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
    assert [r.check.id for r in battery.results] == list(PK_BATTERY.order)
    status = {r.check.id: r.status for r in battery.results}
    for check in (
        "mass_balance",
        "non_negativity",
        "zero_dose",
        "steady_state",
        "superposition",
        "plausibility",
    ):
        assert status[check] is BatteryStatus.PASSED, battery.result(check)
    assert status["cross_simulator"] is BatteryStatus.NOT_APPLICABLE
    assert status["rebuild"] is BatteryStatus.NOT_APPLICABLE
    mass = battery.result("mass_balance")
    # the run without a dose and one per dose
    assert mass.cases == 4
    assert mass.largest < 1e-10


def _fates(battery: Battery) -> list[Fate]:
    return [Fate.from_dict(f) for f in battery.sections["fates"]]


def test_the_fate_of_a_dose(tmp_path: Path) -> None:
    battery = _battery(_model(tmp_path))
    fates = {f.dose: f for f in _fates(battery)}
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
    result = _battery(model).result("mass_balance")
    assert result.status is BatteryStatus.FAILED
    assert any("male IVDOSE_d" in p for p in result.problems)
    assert any("male Ri_d" in p for p in result.problems)
    assert not any("PODOSE_d" in p for p in result.problems)


def test_a_rule_on_a_species_breaks_the_mass_balance(tmp_path: Path) -> None:
    model = _model(tmp_path, "species Cbo_d in Vve = 0; Cbo_d' = 0.001*Cve_d;")
    result = _battery(model).result("mass_balance")
    assert result.status is BatteryStatus.FAILED
    assert (
        "male: the species 'Cbo_d' is changed by a rule or an event"
        in result.problems[0]
    )


def test_a_production_of_the_drug_fails_the_zero_dose_and_the_steady_state(
    tmp_path: Path,
) -> None:
    battery = _battery(_model(tmp_path, "PROD: => Cve_d; 1e-6;"))
    zero = battery.result("zero_dose")
    assert zero.status is BatteryStatus.FAILED
    assert any("'Cve_d' of the substance 'd'" in p for p in zero.problems)
    assert battery.result("steady_state").status is BatteryStatus.FAILED
    # the production is a flow of the balance, which still closes
    assert battery.result("mass_balance").status is BatteryStatus.PASSED


def test_a_negative_species_fails_the_non_negativity(tmp_path: Path) -> None:
    battery = _battery(_model(tmp_path, "LOSS: Cli_m => ; 1e-6;"))
    result = battery.result("non_negativity")
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
    result = _battery(model).result("superposition")
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
    steady = battery.result("steady_state")
    assert steady.status is status
    if status is BatteryStatus.FAILED:
        assert "'Cve_e' changes by" in steady.problems[0]
    # the drug still starts at zero, the endogenous substance is no case
    assert battery.result("zero_dose").cases == 3


def test_a_model_without_a_substance(tmp_path: Path) -> None:
    battery = run_battery(TOY, rebuild=False, cross_simulator=False)
    for check in ("mass_balance", "zero_dose", "superposition"):
        result = battery.result(check)
        assert result.status is BatteryStatus.NOT_APPLICABLE
        assert "Mr_<s>" in result.statement
    # the steady state leaves a build without a substance out
    steady = battery.result("steady_state")
    assert steady.status is BatteryStatus.NOT_APPLICABLE
    assert "Mr_<s>" in steady.statement
    # the other generic checks of sbmlsim run on every model
    assert battery.result("non_negativity").status is BatteryStatus.PASSED
    assert battery.result("plausibility").status is BatteryStatus.PASSED
    assert battery.result("cross_simulator").statement == "Not run."
    assert _fates(battery) == []


def test_a_compartment_without_a_volume_is_implausible(tmp_path: Path) -> None:
    result = _battery(_model(tmp_path, "compartment Vx = 0;")).result("plausibility")
    assert result.status is BatteryStatus.FAILED
    assert result.problems == ("male: the compartment 'Vx' has the size 0",)


def test_the_rebuild(tmp_path: Path) -> None:
    passed = run_battery(CREDIBILITY_MODEL, cross_simulator=False).result("rebuild")
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
    model = CREDIBILITY_MODEL.replace(
        root=tmp_path,
        base_path=CREDIBILITY_MODEL.root,
        variants={"male": copy},
    )
    failed = run_battery(model, cross_simulator=False).result("rebuild")
    assert failed.status is BatteryStatus.FAILED
    assert failed.problems == ("'models/toy.xml' differs from its rebuild",)


def test_a_factory_which_fails(tmp_path: Path) -> None:
    result = run_battery(_model(tmp_path), cross_simulator=False).result("rebuild")
    assert result.status is BatteryStatus.FAILED
    assert result.problems[0].startswith("the factory failed: ")


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
    assert _battery(model).result("mass_balance").status is BatteryStatus.PASSED


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
    assert read.sections == battery.sections
    assert math.isnan(read.result("rebuild").largest)


def _with_battery(battery: Battery) -> CredibilityAssessment:
    fits = fit_table()
    verification = verification_of(fits)
    summary = dataclasses.replace(verification.summary, battery=battery)
    return assess(CREDIBILITY, fits, dataclasses.replace(verification, summary=summary))


def _failed(battery: Battery, check: str) -> Battery:
    return dataclasses.replace(
        battery,
        results=tuple(
            dataclasses.replace(r, status=BatteryStatus.FAILED, problems=("broken",))
            if r.check.id == check
            else r
            for r in battery.results
        ),
    )


def test_the_battery_in_an_assessment(tmp_path: Path) -> None:
    battery = _battery(_model(tmp_path))
    assessment = _with_battery(battery)
    applicable = [r.check.id for r in battery.results if r.passed is not None]
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
    broken = _with_battery(_failed(battery, "steady_state"))
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
    assert list(table["check"]) == list(PK_BATTERY.order)
    assert set(table["status"]) == {"passed", "not_applicable"}
    assert assessment.verification is not None
    assert assessment.verification.battery is not None
    fates = fates_table(assessment.verification.battery.sections["fates"])
    assert set(fates["dose"]) == {"IVDOSE_d", "Ri_d", "PODOSE_d"}
    iv = fates[fates["dose"] == "IVDOSE_d"].set_index("part")
    assert iv.loc["ELIM", "share"] < 0
    output = write_outputs(assessment, tmp_path / "credibility", EXTENSIONS, PK_BATTERY)
    assert (
        (output / "battery.tsv")
        .read_text(encoding="utf-8")
        .startswith("check\tfactor\tstatus")
    )
    html = (output / "index.html").read_text(encoding="utf-8")
    assert "Battery of the variants" in html
    assert "Fate of the doses" in html
    assert "battery of the variants</a>, 6 checks" in html
    assert (output / "fates.tsv").read_text(encoding="utf-8").startswith("build\tdose")


def _pk_battery(project: Project) -> Battery:
    """Run the PK battery of a model, the hook of `check_battery` of sbmlsim."""
    assert isinstance(project, PKPDModel)
    return run_battery(project, rebuild=False, cross_simulator=False)


def test_check_battery(tmp_path: Path) -> None:
    model = _model(tmp_path, "PROD: => Cve_d; 1e-6;")
    with pytest.raises(AssertionError, match="the check 'zero_dose' failed"):
        check_battery(model, _pk_battery)
    known = model.replace(
        checks=CheckSettings(
            battery={"zero_dose": "a production", "steady_state": "a production"}
        ),
    )
    check_battery(known, _pk_battery)
    stale = _model(tmp_path).replace(
        checks=CheckSettings(battery={"zero_dose": "fixed since"})
    )
    with pytest.raises(
        AssertionError, match="'zero_dose' no longer applies, it passed"
    ):
        check_battery(stale, _pk_battery)


def test_the_known_failures_of_the_battery_are_checked(tmp_path: Path) -> None:
    model = _model(tmp_path)
    with pytest.raises(ValueError, match=r"\['unknown'\] are no checks of the battery"):
        model.replace(checks=CheckSettings(battery={"unknown": "reason"}))
    with pytest.raises(ValueError, match="'zero_dose' needs a name and a reason"):
        CheckSettings(battery={"zero_dose": " "})


def test_without_amici_the_second_simulator_is_not_applicable(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(core, "_installed", lambda engine: False)
    battery = run_battery(_model(tmp_path), rebuild=False)
    result = battery.result("cross_simulator")
    assert result.status is BatteryStatus.NOT_APPLICABLE
    assert "AMICI is not installed" in result.statement


@pytest.mark.skipif(
    importlib.util.find_spec("amici") is None, reason="AMICI is not installed"
)
def test_the_second_simulator(tmp_path: Path) -> None:
    result = run_battery(_model(tmp_path), rebuild=False).result("cross_simulator")
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


def test_the_battery_hook_of_a_model() -> None:
    """Every model has the PK battery: two days, a model at its steady state, AMICI."""
    assert CREDIBILITY_MODEL.battery_hook() is PK_BATTERY
    settings = core.battery_settings_of(TOY, PK_BATTERY)
    assert (settings.horizon, settings.steps) == (HORIZON, STEPS) == (2880.0, 288)
    assert settings.steady_state
    assert settings.engines == ("amici",)
    assert core.battery_settings_of(CREDIBILITY_MODEL, PK_BATTERY).steady_state
    # without a battery the table of the fates has no rows
    (fates,) = PK_BATTERY.report_tables({}).values()
    assert fates.empty
    assert list(fates.columns) == [
        "build",
        "dose",
        "substance",
        "amount",
        "part",
        "value",
        "share",
    ]
    for check in PK_CHECKS:
        assert CHECKS[check.id] == check
        assert check.factor == "numerical_code_verification"
    assert sorted(PK_BATTERY.order) == sorted(
        [*(c.id for c in core.BATTERY_CHECKS), *(c.id for c in PK_CHECKS)]
    )


def test_the_cases_of_the_doses(tmp_path: Path) -> None:
    build = FlatBuild.read("male", _model(tmp_path).model_path)
    cases = {case.name: case for case in dose_cases(build)}
    assert list(cases)[:5] == [
        "IVDOSE_d",
        "IVDOSE_d small at 0",
        "IVDOSE_d small twice at 0",
        "IVDOSE_d small at 720 min",
        "IVDOSE_d small at 0 and 720 min",
    ]
    # the second simulator repeats the doses, not the superposition
    assert [n for n, c in cases.items() if c.cross_simulator] == [
        "IVDOSE_d",
        "Ri_d",
        "PODOSE_d",
    ]
    assert cases["IVDOSE_d"].added == {"IVDOSE_d"}
    # an infusion over an hour, stopped
    assert cases["Ri_d"].changes == {0.0: {"Ri_d": DOSE / 60.0}, 60.0: {"Ri_d": 0.0}}
    assert cases["Ri_d"].added == frozenset()
