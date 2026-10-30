"""The drug-free equilibration of a condition (`PBPKSimulationExperiment.equilibrated`).

A PD component with a homeostasis, `X' = kin*f*(1 - E(C)) - kout*X`, starts at
the steady state of the healthy condition `f = 1` (`X = 100`); a condition
`f = 0.5` (an impairment) moves it to `X = 50` with the time constant
`1/kout`, 7 days. The second component also has an accumulating sink,
`V' = kout*X` (the urine of a PD model), which has no steady state.
"""

from collections.abc import Callable
from typing import cast

import numpy as np
import pytest
from sbmlsim import Q
from sbmlsim.experiment import ExperimentRunner
from sbmlsim.model import AbstractModel
from sbmlsim.result import TimecourseResult
from sbmlsim.simulation import (
    Change,
    Dimension,
    Scan,
    Simulation,
    SteadyState,
)
from sbmlsim.simulator import Simulator
from sbmlutils.factory import Model, Parameter, Port, PortType, RateRule, Sbase
from sbmlutils.metadata import SBO

from pkpdsim.model import PKPDModel
from pkpdsim.pbpk import MODEL_UNITS, PDComponent, Sex, U, build, ids
from pkpdsim.pbpk.builder import BuildResult
from pkpdsim.pbpk.experiment import PBPKSimulationExperiment

from .conftest import definition

#: the elimination of the homeostasis [1/min], a time constant of 7 days
KOUT = 1e-4
#: the condition, the steady state of the homeostasis in it
CONDITION = {"HOMEO__f": Q(0.5, "dimensionless")}
STEADY_X = 50.0


def homeostasis(*, sink: bool) -> PDComponent:
    """Get the PD component of the homeostasis, with or without the sink."""
    sbo = SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER
    objects: list[Sbase] = [
        Parameter(
            "C_drug",
            0.0,
            U.mM,
            constant=False,
            sboTerm=sbo,
            port=Port("C_drug_port", portType=PortType.INPUT_PORT),
        ),
        Parameter("f", 1.0, U.dimensionless, sboTerm=sbo),
        Parameter("kin", KOUT * 100, U.per_min, sboTerm=sbo),
        Parameter("kout", KOUT, U.per_min, sboTerm=sbo),
        Parameter("EC50", 1e-3, U.mM, sboTerm=sbo),
        Parameter("X", 100.0, U.dimensionless, constant=False, sboTerm=sbo),
        RateRule(
            "X",
            "kin * f * (1 dimensionless - C_drug / (EC50 + C_drug)) - kout * X",
            U.per_min,
        ),
    ]
    if sink:
        objects += [
            Parameter("V", 0.0, U.dimensionless, constant=False, sboTerm=sbo),
            RateRule("V", "kout * X", U.per_min),
        ]
    model = Model("HOMEO", units=U, model_units=MODEL_UNITS, objects=objects)
    return PDComponent(model=model, inputs={"C_drug": ids.arm_free("drg")}, outputs={})


@pytest.fixture(scope="module")
def builds(tmp_path_factory: pytest.TempPathFactory) -> dict[str, BuildResult]:
    """The homeostasis without and with the sink, male; with the sink female."""
    out = tmp_path_factory.mktemp("equilibration")
    return {
        "homeostasis": build(
            definition(sid="homeostasis", pd=[homeostasis(sink=False)]), out
        ),
        "sink": build(definition(sid="sink", pd=[homeostasis(sink=True)]), out),
        "sink_female": build(
            definition(sid="sink_female", sex=Sex.FEMALE, pd=[homeostasis(sink=True)]),
            out,
        ),
    }


def _model(result: BuildResult, female: BuildResult | None = None) -> PKPDModel:
    return PKPDModel(
        name="toy_equilibration",
        package="toy_equilibration",
        title="Toy equilibration",
        path=result.flat_path.parent,
        model_path=result.flat_path,
        female_model_path=None if female is None else female.flat_path,
        factory=".factory:create_models",
    )


def _run(
    model: PKPDModel,
    simulations: Callable[[PBPKSimulationExperiment], dict[str, Simulation | Scan]],
    sex: Sex = Sex.MALE,
) -> PBPKSimulationExperiment:
    """Initialize and run an experiment of the simulations, as a run does."""

    class Equilibration(PBPKSimulationExperiment):
        pkpd_model = model
        selections = ("time", "HOMEO__X", "Aurine_drg", "Carm_drg", "IVDOSE_drg")

        def simulations(self) -> dict[str, Simulation | Scan]:
            return simulations(self)

    Equilibration.sex = sex
    runner = ExperimentRunner(
        experiment_classes=[Equilibration],
        base_path=model.path,
        simulator=Simulator(
            n_workers=1, absolute_tolerance=1e-10, relative_tolerance=1e-10
        ),
        on_error="raise",
    )
    (experiment,) = runner.experiments.values()
    experiment = cast("PBPKSimulationExperiment", experiment)
    experiment.run(runner.simulator, output_path=None)
    return experiment


def _at_start(experiment: PBPKSimulationExperiment, key: str, sid: str) -> float:
    result = cast("TimecourseResult", experiment.results[f"task_{key}"])
    return float(np.asarray(result[sid]).ravel()[0])


def _dosed(**kwargs: object) -> Simulation:
    """A dose of the drug in the condition, given before the initialization."""
    return Simulation(
        end=600,
        preinit_changes={**CONDITION, "IVDOSE_drg": Q(10, "mg")},
        **kwargs,  # ty: ignore[invalid-argument-type]
    )


def test_the_steady_state_equilibrates_a_homeostasis_and_a_sink(
    builds: dict[str, BuildResult],
) -> None:
    """The choice: a presimulation to the steady state, not a pre-run of a duration.

    The presimulation of sbmlsim integrates until the rates of change are
    within its tolerances; an accumulating sink reaches them as well, its rate
    relative to its growing value, and is set back at the start. A hidden
    pre-run of 28 days from a negative start leaves the homeostasis of 7 days
    1.8 % off its steady state, and the default tolerances of sbmlsim, per
    minute, stop the slow homeostasis early.
    """
    for name in ("homeostasis", "sink"):
        model = _model(builds[name])
        resets = {"HOMEO__V": 0.0} if name == "sink" else None

        def simulations(
            experiment: PBPKSimulationExperiment,
            resets: dict[str, float] | None = resets,
        ) -> dict[str, Simulation | Scan]:
            return {
                "equilibrated": experiment.equilibrated(_dosed(), resets=resets),
                "prerun": Simulation(
                    start=-28 * 24 * 60,
                    end=600,
                    preinit_changes=CONDITION,
                    changes=[Change(0, {"IVDOSE_drg": Q(10, "mg")})],
                ),
                "default_tolerances": Simulation(
                    end=600,
                    preinit_changes=CONDITION,
                    changes=[Change(0, {"IVDOSE_drg": Q(10, "mg")})],
                    presimulation=SteadyState(),
                ),
            }

        experiment = _run(model, simulations)
        x = _at_start(experiment, "equilibrated", "HOMEO__X")
        assert x == pytest.approx(STEADY_X, rel=1e-6), name
        result = cast("TimecourseResult", experiment.results["task_prerun"])
        times = np.asarray(result["time"]).ravel()
        prerun = float(np.asarray(result["HOMEO__X"]).ravel()[np.argmin(np.abs(times))])
        assert abs(prerun / STEADY_X - 1) > 1e-2
        if name == "homeostasis":
            early = _at_start(experiment, "default_tolerances", "HOMEO__X")
            assert abs(early / STEADY_X - 1) > 1e-5


def test_the_doses_start_after_the_equilibration(
    builds: dict[str, BuildResult],
) -> None:
    """The doses move to the start, the condition stays, the sinks start empty."""
    simulation = _dosed(changes=[Change(60, {"IVDOSE_drg": Q(5, "mg")})])

    def simulations(
        experiment: PBPKSimulationExperiment,
    ) -> dict[str, Simulation | Scan]:
        equilibrated = experiment.equilibrated(simulation, resets={"HOMEO__V": 0.0})
        assert dict(equilibrated.preinit_changes) == CONDITION
        assert equilibrated.presimulation is not None
        first, later = equilibrated.changes
        assert first.times == (0.0,)
        assert dict(first.values) == {
            "IVDOSE_drg": Q(10, "mg"),
            "Aurine_drg": 0.0,
            "Afeces_drg": 0.0,
            "Aurine_drg_m": 0.0,
            "Afeces_drg_m": 0.0,
            "HOMEO__V": 0.0,
        }
        assert later is simulation.changes[0]
        return {"equilibrated": equilibrated}

    experiment = _run(_model(builds["sink"]), simulations)
    assert _at_start(experiment, "equilibrated", "IVDOSE_drg") > 0
    assert _at_start(experiment, "equilibrated", "Aurine_drg") == 0.0
    result = cast("TimecourseResult", experiment.results["task_equilibrated"])
    # the drug acts on the homeostasis after the dose only
    x = np.asarray(result["HOMEO__X"]).ravel()
    assert x[0] == pytest.approx(STEADY_X, rel=1e-6)
    assert x[-1] < x[0]


def test_a_stationary_condition_is_simulated_as_without_the_equilibration(
    builds: dict[str, BuildResult],
) -> None:
    """In the healthy condition the start is the steady state already."""
    times = list(np.linspace(0, 600, 61))

    def simulations(
        experiment: PBPKSimulationExperiment,
    ) -> dict[str, Simulation | Scan]:
        plain = Simulation(
            end=600, preinit_changes={"IVDOSE_drg": Q(10, "mg")}, times=times
        )
        return {"plain": plain, "equilibrated": experiment.equilibrated(plain)}

    experiment = _run(_model(builds["homeostasis"]), simulations)
    values = {
        key: np.asarray(experiment.results[f"task_{key}"]["Carm_drg"]).ravel()
        for key in ("plain", "equilibrated")
    }
    assert np.allclose(
        values["equilibrated"],
        values["plain"],
        rtol=1e-6,
        atol=1e-6 * values["plain"].max(),
    )


def test_a_scan_is_equilibrated_at_every_point(
    builds: dict[str, BuildResult],
) -> None:
    """A dimension of the condition is a condition of the equilibration."""

    def simulations(
        experiment: PBPKSimulationExperiment,
    ) -> dict[str, Simulation | Scan]:
        scan = Scan(
            simulation=Simulation(end=60, preinit_changes={"IVDOSE_drg": Q(10, "mg")}),
            dimensions=[
                Dimension("dim_f", values={"HOMEO__f": Q([0.5, 2.0], "dimensionless")}),
                Dimension("dim_dose", values={"IVDOSE_drg": Q([5.0, 10.0], "mg")}),
            ],
        )
        return {"scan": experiment.equilibrated(scan)}

    experiment = _run(_model(builds["homeostasis"]), simulations)
    from pkpdsim.scan import start_values

    x = np.asarray(start_values(experiment.results["task_scan"], "HOMEO__X"))
    assert np.allclose(x, [[50.0, 50.0], [200.0, 200.0]], rtol=1e-6)


def test_the_female_build_is_equilibrated(builds: dict[str, BuildResult]) -> None:
    """A study of women equilibrates on the female build."""
    model = _model(builds["sink"], builds["sink_female"])

    def simulations(
        experiment: PBPKSimulationExperiment,
    ) -> dict[str, Simulation | Scan]:
        return {"women": experiment.equilibrated(_dosed(), resets={"HOMEO__V": 0.0})}

    experiment = _run(model, simulations, sex=Sex.FEMALE)
    female = experiment.models()["model"]
    assert isinstance(female, AbstractModel)
    assert female.source.path == builds["sink_female"].flat_path
    assert _at_start(experiment, "women", "HOMEO__X") == pytest.approx(
        STEADY_X, rel=1e-6
    )


@pytest.mark.parametrize(
    ("simulation", "message"),
    [
        (
            _dosed(presimulation=SteadyState()),
            "has a presimulation",
        ),
        (
            _dosed(changes=[Change(0, {"IVDOSE_drg": "IVDOSE_drg + 1"})]),
            "before the initialization and by a change at its start",
        ),
    ],
)
def test_a_simulation_which_cannot_be_equilibrated_is_refused(
    builds: dict[str, BuildResult], simulation: Simulation, message: str
) -> None:
    def simulations(
        experiment: PBPKSimulationExperiment,
    ) -> dict[str, Simulation | Scan]:
        return {"refused": experiment.equilibrated(simulation)}

    with pytest.raises(ValueError, match=message):
        _run(_model(builds["homeostasis"]), simulations)


def test_a_dose_of_a_scan_must_be_a_dose_of_its_simulation(
    builds: dict[str, BuildResult],
) -> None:
    """A dimension over a dose which the simulation does not give is refused.

    The dimension would set the dose before the initialization, which the
    presimulation would absorb.
    """

    def simulations(
        experiment: PBPKSimulationExperiment,
    ) -> dict[str, Simulation | Scan]:
        scan = Scan(
            simulation=Simulation(end=60, preinit_changes=CONDITION),
            dimensions=[Dimension("dim_dose", values={"IVDOSE_drg": Q([5.0], "mg")})],
        )
        return {"scan": experiment.equilibrated(scan)}

    with pytest.raises(ValueError, match="IVDOSE_drg"):
        _run(_model(builds["homeostasis"]), simulations)


def test_a_model_must_be_loaded(builds: dict[str, BuildResult]) -> None:
    """The ids of the doses and the sinks are those of a loaded model."""

    def simulations(
        experiment: PBPKSimulationExperiment,
    ) -> dict[str, Simulation | Scan]:
        return {"other": experiment.equilibrated(_dosed(), model="other")}

    with pytest.raises(ValueError, match="not loaded"):
        _run(_model(builds["homeostasis"]), simulations)
