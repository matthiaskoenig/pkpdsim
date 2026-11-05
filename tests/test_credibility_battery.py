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

from pkpdsim.credibility import DOSE_STATES, EXTENSIONS, PK_BATTERY, PK_CHECKS
from pkpdsim.credibility.battery import (
    _INFUSION,
    _SINKS,
    DOSE,
    HORIZON,
    LEGACY_DOSE_STATES,
    STEPS,
    Fate,
    FlatBuild,
    _bolus,
    dose_cases,
    fates_table,
)
from pkpdsim.model import PKPDModel
from pkpdsim.pbpk import ids
from pkpdsim.pbpk.equilibration import DOSE_PREFIXES as LIBRARY_DOSES
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
        "dose_states",
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
    # the dose states of a build are read without a substance
    assert battery.result("dose_states").status is BatteryStatus.PASSED
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


def test_the_dose_states_of_the_library_are_doses_of_the_battery() -> None:
    """Every dose state the PBPK library creates is a dose of the battery.

    A subset and not an equality: the battery also knows the dose states of
    the models which are not built with the library
    (`LEGACY_DOSE_STATES`). A route the battery does not know would leave its
    dose out of the mass balance, the superposition and the fates.

    What each line pins: the first cannot fail while `dose_prefixes` is
    derived from the list of the library, and says that it must stay derived;
    the second pins the one id the battery spells itself, the infusion rate;
    the third and the fourth pin the legacy list, which is written by hand
    and is what a future entry is added to.
    """
    assert set(LIBRARY_DOSES) <= {f"{prefix}_" for prefix in (*_bolus(), _INFUSION)}
    # the battery spells the infusion rate itself, the library names it
    assert f"{_INFUSION}_" == ids.infusion_rate("")
    # a legacy state is a dose and names its model
    assert set(LEGACY_DOSE_STATES) <= {f"{prefix}_" for prefix in _bolus()}
    assert all(LEGACY_DOSE_STATES.values())


def test_the_sinks_of_the_library_are_sinks_of_a_fate() -> None:
    """Every sink the PBPK library writes is a part of a fate of its own.

    Without this, a rename of an id of the library would move the sink into
    the `body` of every fate with the suite green. The bile of the library
    (`Cbile_<s>`) is no sink: it transits into the intestinal lumen.
    """
    sinks = {
        ids.urine(""),
        ids.feces(""),
        ids.depot_loss("im", ""),
        ids.depot_loss("sc", ""),
    }
    assert all(
        any(sink.startswith(prefix) for prefix in _SINKS.values()) for sink in sinks
    ), (sinks, _SINKS)


def test_a_dose_state_the_battery_does_not_know_is_reported(tmp_path: Path) -> None:
    """An unknown state in mg which a reaction reads fails its own check.

    The failure mode of an unlisted dose state is a silent one: its dose is
    simulated by nothing and no fate is recorded. The check `dose_states`
    names it while the mass balance, whose arithmetic holds, passes. The case
    which asked for the split is `ethinylestradiol` of pkdb_models: its patch
    and ring depots are unknown and its amounts balance to 2.3e-15.
    """
    model = _model(
        tmp_path,
        "unit mg = 1e-3 gram; XDOSE_d = 0; XDOSE_d has mg; kx = 0.3;"
        " XDOSE_d' = -kx*XDOSE_d;"
        " x_d: => Cve_d; kx*XDOSE_d/Mr_d;",
    )
    build = FlatBuild.read("male", model.model_path)
    assert build.unknown_doses == ("XDOSE_d",)
    assert "XDOSE_d" not in build.doses
    battery = _battery(model)
    # the arithmetic of the mass balance holds and is asserted first, so the
    # test fails here with the behaviour of the round before the split
    mass = battery.result("mass_balance")
    assert mass.status is BatteryStatus.PASSED, mass.problems
    assert mass.largest < 1e-10
    result = battery.result("dose_states")
    assert result.status is BatteryStatus.FAILED
    assert result.cases == 1
    assert any("XDOSE_d" in problem for problem in result.problems), result.problems


def test_a_dose_which_disappears_fails_only_the_mass_balance(tmp_path: Path) -> None:
    """A dose of a known prefix which is not conserved fails the arithmetic.

    The other side of the split: the mass balance fails and the dose states
    of the model are all known. The case is the inhaled route of `albuterol`
    of pkdb_models, which drops the fraction of the dose that stays in the
    device or is exhaled without a species, 0.279 of the dose.
    """
    model = _model(
        tmp_path,
        "unit mg = 1e-3 gram; IMDOSE_d = 0; IMDOSE_d has mg; kim = 0.3;"
        " f_lost = 0.25;"
        " IMDOSE_d' = -kim*IMDOSE_d;"
        " im_d: => Cve_d; (1 - f_lost)*kim*IMDOSE_d/Mr_d;",
    )
    build = FlatBuild.read("male", model.model_path)
    assert build.unknown_doses == ()
    assert build.dose_states["IMDOSE_d"] == "d"
    battery = _battery(model)
    mass = battery.result("mass_balance")
    assert mass.status is BatteryStatus.FAILED
    assert mass.largest == pytest.approx(0.25, rel=1e-2)
    assert any("IMDOSE_d" in problem for problem in mass.problems), mass.problems
    assert battery.result("dose_states").status is BatteryStatus.PASSED


def test_a_dose_state_of_a_submodel(tmp_path: Path) -> None:
    """A dose state behind the qualifier of a flattened submodel is a dose.

    A flat model qualifies the ids of a submodel with `<SUBMODEL>__`, and a
    dose state of a submodel is a dose state of the model, as a species of a
    submodel belongs to its substance. Before the strip, `GU__PODOSE_d` with
    the prefix `PODOSE_` listed and `Mr_d` present was reported as unknown
    and the remedy the message named, adding the prefix, was already
    satisfied (`GU__PODOSE_ethanol` of chlorzoxazone is the live case).
    """
    model = _model(
        tmp_path,
        "unit mg = 1e-3 gram; GU__PODOSE_d = 0; GU__PODOSE_d has mg; kg = 0.4;"
        " GU__PODOSE_d' = -kg*GU__PODOSE_d;"
        " gu_d: => Cve_d; kg*GU__PODOSE_d/Mr_d;",
    )
    build = FlatBuild.read("male", model.model_path)
    assert build.unknown_doses == ()
    assert build.doses["GU__PODOSE_d"] == "d"
    assert build.dose_states["GU__PODOSE_d"] == "d"
    battery = _battery(model)
    assert battery.result("dose_states").status is BatteryStatus.PASSED
    mass = battery.result("mass_balance")
    assert mass.status is BatteryStatus.PASSED, mass.problems
    # the qualified dose has its own fate, and the run without a dose and one
    # per dose are the cases
    fates = {f.dose: f for f in _fates(battery)}
    assert "GU__PODOSE_d" in fates
    assert mass.cases == 5
    assert fates["GU__PODOSE_d"].parts["urine"] > 0


def test_a_qualified_dose_is_added_and_not_set(tmp_path: Path) -> None:
    """The superposition adds a dose to a qualified dose state.

    A dose of a dose state is added to its value, a dose of an infusion rate
    replaces it; the decision reads the prefix of the id, which is `GU` for a
    qualified id unless the qualifier is stripped, so the second dose of the
    superposition would have replaced the first.
    """
    model = _model(
        tmp_path,
        "unit mg = 1e-3 gram; GU__PODOSE_d = 0; GU__PODOSE_d has mg; kg = 0.4;"
        " GU__PODOSE_d' = -kg*GU__PODOSE_d;"
        " gu_d: => Cve_d; kg*GU__PODOSE_d/Mr_d;",
    )
    # imported here and not at the top of the module, so that the module
    # imports with the battery of the round before the strip and
    # `test_a_dose_state_of_a_submodel` fails on its behaviour there
    from pkpdsim.credibility.battery import _is_bolus

    assert _is_bolus("GU__PODOSE_d")
    assert not _is_bolus("GU__Ri_d")
    result = _battery(model).result("superposition")
    assert result.status is BatteryStatus.PASSED, result.problems


def test_an_injected_depot(tmp_path: Path) -> None:
    """The dose of a depot route is a dose, and its loss a part of the fate."""
    model = _model(
        tmp_path,
        "unit mg = 1e-3 gram; IMDOSE_d = 0; IMDOSE_d has mg;"
        " substanceOnly species Aim_d in Vve = 0;"
        " substanceOnly species Aloss_im_d in Vve = 0;"
        " ti_im = 10; ka_im = 0.5; f_im = 0.8;"
        " IMDOSE_d' = -ti_im*IMDOSE_d;"
        " im_app_d: => Aim_d; ti_im*IMDOSE_d/Mr_d;"
        " im_abs_d: Aim_d => Cve_d; f_im*ka_im*Aim_d;"
        " im_loss_d: Aim_d => Aloss_im_d; (1 - f_im)*ka_im*Aim_d;",
    )
    build = FlatBuild.read("male", model.model_path)
    assert build.dose_states["IMDOSE_d"] == "d"
    assert build.doses["IMDOSE_d"] == "d"
    battery = _battery(model)
    assert battery.result("mass_balance").status is BatteryStatus.PASSED
    # the run without a dose and one per dose, the intramuscular one included
    assert battery.result("mass_balance").cases == 5
    fate = {f.dose: f for f in _fates(battery)}["IMDOSE_d"]
    assert fate.parts["loss"] == pytest.approx(0.2 * fate.amount, rel=1e-6)
    assert fate.parts["urine"] > 0


def test_the_depot_routes_of_a_build_of_the_library(tmp_path: Path) -> None:
    """The fates of the depot doses of a build of the PBPK library.

    `test_an_injected_depot` imitates the ids of the library in antimony;
    this one builds them, so a rename of `ids.depot_loss` or of a dose state
    fails here instead of moving the loss of a depot back into the body.
    """
    from pkpdsim.pbpk import Binding, Excretion, Kp, PBPKDefinition, Route, Substance
    from pkpdsim.pbpk.builder import build

    definition = PBPKDefinition(
        sid="depots",
        substances=[
            Substance(
                sid="drg",
                name="drug",
                mr=300.0,
                binding=Binding(fu=0.1, protein="albumin"),
                distribution=Kp(2.0, ftissue=1.0),
                excretion=Excretion(urine=True),
            )
        ],
        chemistry=[],
        routes=[
            Route.iv("drg"),
            Route.im("drg", ka=0.5, fraction=0.8),
            Route.sc("drg", ka=0.3, fraction=0.7),
            Route.rectal("drg", ka=0.4, fraction=0.6, portal=0.3),
        ],
    )
    flat = build(definition, tmp_path).flat_path
    model = PKPDModel(
        name="depots",
        package="depots",
        title="depots",
        path=tmp_path,
        model_path=flat,
        factory="depots:create_depots_models",
    )
    build_read = FlatBuild.read("male", flat)
    assert build_read.unknown_doses == ()
    battery = _battery(model)
    mass = battery.result("mass_balance")
    assert mass.status is BatteryStatus.PASSED, mass.problems
    fates = {f.dose: f for f in _fates(battery)}
    assert {"IVDOSE_drg", "IMDOSE_drg", "SCDOSE_drg", "RECDOSE_drg"} <= set(fates)
    # the fraction which is not absorbed is the loss of an injected depot
    for dose, absorbed in (("IMDOSE_drg", 0.8), ("SCDOSE_drg", 0.7)):
        fate = fates[dose]
        assert fate.parts["loss"] == pytest.approx(
            (1 - absorbed) * fate.amount, rel=1e-4
        ), dose
    # the unabsorbed rectal dose leaves into the feces, it has no loss
    rectal = fates["RECDOSE_drg"]
    assert "loss" not in rectal.parts
    assert rectal.parts["feces"] == pytest.approx(0.4 * rectal.amount, rel=1e-4)
    assert "loss" not in fates["IVDOSE_drg"].parts
    # the amounts account for the whole dose: the model has no elimination
    # reaction, its urine is a species
    parts = {"body", "urine", "feces", "loss", "dose"}
    for dose, fate in fates.items():
        total = sum(v for k, v in fate.parts.items() if k in parts)
        assert total == pytest.approx(fate.amount, rel=1e-6), dose


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
    # the applicable checks of the drug model, the dose states included
    assert "battery of the variants</a>, 7 checks" in html
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
    assert PK_CHECKS[0] is DOSE_STATES
    for check in PK_CHECKS:
        assert CHECKS[check.id] == check
        # the dose states are a use error, the others verify the numerics
        expected = (
            "use_error" if check is DOSE_STATES else "numerical_code_verification"
        )
        assert check.factor == expected
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
