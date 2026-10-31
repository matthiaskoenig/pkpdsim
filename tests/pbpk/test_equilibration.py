"""The drug-free equilibration of a condition (`PBPKSimulationExperiment.equilibrated`).

A PD component with a homeostasis, `X' = kin*f*(1 - E(C)) - kout*X`, starts at
the steady state of the healthy condition `f = 1` (`X = 100`); a condition
`f = 0.5` (an impairment) moves it to `X = 50` with the time constant
`1/kout`, 7 days. The second component also has an accumulating sink,
`V' = kout*X` (the urine of a PD model), which has no steady state.
"""

import dataclasses
from collections.abc import Callable
from pathlib import Path
from typing import assert_type, cast

import libsbml
import numpy as np
import pandas as pd
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
from sbmlsim.simulator.executor import SteadyStateError
from sbmlsim.simulator.simulator import ScanError
from sbmlutils.factory import Model, Parameter, Port, PortType, RateRule, Sbase
from sbmlutils.metadata import SBO

from pkpdsim.model import PKPDModel, SimulatorSettings
from pkpdsim.pbpk import (
    MODEL_UNITS,
    CustomReaction,
    Excretion,
    FirstOrder,
    Kp,
    Metabolism,
    ModelParameter,
    Organ,
    Part,
    PBPKDefinition,
    PDComponent,
    Pool,
    Route,
    Sex,
    Substance,
    U,
    build,
    ids,
)
from pkpdsim.pbpk.builder import BuildResult
from pkpdsim.pbpk.equilibration import (
    accumulating_states,
    drug_free_zeros,
    equilibrated_simulation,
)
from pkpdsim.pbpk.experiment import PBPKSimulationExperiment

from .conftest import definition, endogenous_reference_definition

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


def oscillator() -> PDComponent:
    """Get a PD component which oscillates: it has no steady state."""
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
        Parameter("w", 0.01, U.per_min, sboTerm=sbo),
        Parameter("X", 1.0, U.dimensionless, constant=False, sboTerm=sbo),
        Parameter("Y", 0.0, U.dimensionless, constant=False, sboTerm=sbo),
        RateRule("X", "w * Y", U.per_min),
        RateRule("Y", "-w * X", U.per_min),
    ]
    model = Model("OSC", units=U, model_units=MODEL_UNITS, objects=objects)
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
        "oscillator": build(definition(sid="oscillator", pd=[oscillator()]), out),
        "endogenous": build(
            dataclasses.replace(
                endogenous_reference_definition(), deltas=(ids.arm("end"),)
            ),
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
    selections: tuple[str, ...] = (
        "time",
        "HOMEO__X",
        "Aurine_drg",
        "Carm_drg",
        "IVDOSE_drg",
    ),
) -> PBPKSimulationExperiment:
    """Initialize and run an experiment of the simulations, as a run does."""

    chosen = selections

    class Equilibration(PBPKSimulationExperiment):
        pkpd_model = model
        selections = chosen

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


def _concentrations(path: Path) -> set[str]:
    """The species of a model whose ids are their concentrations."""
    model = libsbml.readSBMLFromFile(str(path)).getModel()
    return {
        sp.getId()
        for sp in model.getListOfSpecies()
        if not sp.getHasOnlySubstanceUnits()
    }


def _zero_targets(path: Path) -> list[str]:
    """The targets of the drug-free zeros, `[sid]` for a concentration."""
    concentrations = _concentrations(path)
    return [
        f"[{sid}]" if sid in concentrations else sid for sid in drug_free_zeros(path)
    ]


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
        zeros = _zero_targets(builds["sink"].flat_path)
        assert dict(first.values) == {
            "IVDOSE_drg": Q(10, "mg"),
            "Aurine_drg": 0.0,
            "Afeces_drg": 0.0,
            "Aurine_drg_m": 0.0,
            "Afeces_drg_m": 0.0,
            **dict.fromkeys(zeros, 0.0),
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


def test_the_accumulating_states_of_a_model(builds: dict[str, BuildResult]) -> None:
    """The states which nothing drains: the sinks, not the homeostasis."""
    sink = accumulating_states(builds["sink"].flat_path.read_text(encoding="utf-8"))
    assert "HOMEO__V" in sink
    assert {"Aurine_drg", "Afeces_drg"} <= set(sink)
    assert "HOMEO__X" not in sink
    plain = builds["homeostasis"].flat_path.read_text(encoding="utf-8")
    assert "HOMEO__V" not in accumulating_states(plain)
    # two states which read each other are no sinks
    oscillating = builds["oscillator"].flat_path.read_text(encoding="utf-8")
    assert not {"OSC__X", "OSC__Y"} & set(accumulating_states(oscillating))


def test_a_growing_state_without_a_reset_is_refused(
    builds: dict[str, BuildResult],
) -> None:
    """The urine of the PD component grows; without a reset its start is meaningless."""

    def simulations(
        experiment: PBPKSimulationExperiment,
    ) -> dict[str, Simulation | Scan]:
        return {"forgotten": experiment.equilibrated(_dosed())}

    with pytest.raises(ValueError, match=r"\['HOMEO__V'\] grow without a steady"):
        _run(_model(builds["sink"]), simulations)


def test_a_model_without_a_steady_state_fails(builds: dict[str, BuildResult]) -> None:
    """An oscillation reaches no steady state: sbmlsim's SteadyStateError."""

    def simulations(
        experiment: PBPKSimulationExperiment,
    ) -> dict[str, Simulation | Scan]:
        dosed = Simulation(end=60, preinit_changes={"IVDOSE_drg": Q(10, "mg")})
        return {"oscillating": experiment.equilibrated(dosed, max_time=1e4)}

    with pytest.raises(ScanError, match="did not reach a steady state") as raised:
        _run(_model(builds["oscillator"]), simulations, selections=("time",))
    causes = []
    error: BaseException | None = raised.value
    while error is not None:
        causes.append(type(error))
        error = error.__cause__
    assert SteadyStateError in causes


def test_a_delta_starts_after_the_equilibration(builds: dict[str, BuildResult]) -> None:
    """The start value of a delta is the value after the equilibration, the delta 0."""
    arm = ids.arm("end")

    def simulations(
        experiment: PBPKSimulationExperiment,
    ) -> dict[str, Simulation | Scan]:
        condition = Simulation(
            end=60, preinit_changes={"f_renal_function": Q(0.5, "dimensionless")}
        )
        return {
            "equilibrated": experiment.equilibrated(condition),
            "plain": condition,
        }

    experiment = _run(
        _model(builds["endogenous"]),
        simulations,
        selections=("time", arm, ids.delta(arm)),
    )

    def at_start(key: str, sid: str) -> float:
        result = cast("TimecourseResult", experiment.results[f"task_{key}"])
        return float(np.asarray(result[sid]).ravel()[0])

    # the equilibration moves the arm vein concentration in the impairment
    assert at_start("equilibrated", arm) > 1.05 * at_start("plain", arm)
    assert at_start("equilibrated", ids.delta(arm)) == pytest.approx(
        0.0, abs=1e-12 * at_start("equilibrated", arm)
    )


# --- the species of the drugs start at exactly 0 (ruling R25) ---------------------


def _species(path: Path, *substances: str) -> set[str]:
    model = libsbml.readSBMLFromFile(str(path)).getModel()
    every = [s.getId() for s in model.getListOfSpecies()]
    return {sid for s in substances for sid in ids.species_of(s, every)}


def test_the_species_of_the_drugs_are_drug_free_zeros(
    builds: dict[str, BuildResult],
) -> None:
    """The drug and its metabolite: every species of them, not the sinks or the PD."""
    path = builds["sink"].flat_path
    zeros = set(drug_free_zeros(path))
    assert zeros == _species(path, "drg", "drg_m")
    assert len(zeros) > 20
    # a species which the condition sets holds a value, and so does what it feeds
    venous = ids.pool_species(Pool.VENOUS, "drg")
    assert venous in zeros
    fed = set(drug_free_zeros(path, {venous}))
    assert fed < zeros
    # the venous blood reaches every organ, not the lumen of the intestine
    assert all("lumen" in sid or "Asolid" in sid for sid in fed)
    assert drug_free_zeros(path, {"HOMEO__f"}) == tuple(sorted(zeros))


def test_an_endogenous_substance_and_its_products_are_no_zeros(
    tmp_path: Path,
) -> None:
    """A baseline starts the substance, its metabolism feeds the metabolite."""
    d = PBPKDefinition(
        sid="endogenous_products",
        substances=[
            Substance(
                "end",
                name="endogenous substance",
                mr=113.0,
                distribution=Kp(0.6, ftissue=6.0),
                excretion=Excretion(urine=True),
                baseline=0.02,
            ),
            Substance("end_m", name="its metabolite", mr=100.0),
            Substance("drg", name="drug", mr=300.0, distribution=Kp(2.0, ftissue=1.0)),
        ],
        chemistry=[
            Metabolism(
                "END2M",
                Organ.LIVER,
                "end",
                {"end_m": 1.0},
                FirstOrder(k=0.01),
                protein="CYP3A4",
            )
        ],
        routes=[Route.iv("drg"), Route.oral("end")],
    )
    path = build(d, tmp_path).flat_path
    zeros = set(drug_free_zeros(path))
    assert zeros == _species(path, "drg")
    assert not zeros & _species(path, "end", "end_m")


def test_the_drug_starts_at_exactly_zero_after_the_equilibration(
    builds: dict[str, BuildResult],
) -> None:
    """A placebo in the condition: every species of the drug is 0, not a residue."""
    path = builds["sink"].flat_path
    zeros = drug_free_zeros(path)

    def simulations(
        experiment: PBPKSimulationExperiment,
    ) -> dict[str, Simulation | Scan]:
        placebo = Simulation(end=600, preinit_changes=CONDITION)
        return {
            "placebo": experiment.equilibrated(placebo, resets={"HOMEO__V": 0.0}),
            "dosed": experiment.equilibrated(_dosed(), resets={"HOMEO__V": 0.0}),
        }

    experiment = _run(
        _model(builds["sink"]),
        simulations,
        selections=("time", "HOMEO__X", "Carm_drg", *(f"[{z}]" for z in zeros)),
    )
    for key in ("placebo", "dosed"):
        result = cast("TimecourseResult", experiment.results[f"task_{key}"])
        for sid in ("Carm_drg", *(f"[{z}]" for z in zeros)):
            assert float(np.asarray(result[sid]).ravel()[0]) == 0.0, (key, sid)
    placebo = cast("TimecourseResult", experiment.results["task_placebo"])
    assert np.all(np.asarray(placebo["Carm_drg"]).ravel() == 0.0)
    assert _at_start(experiment, "placebo", "HOMEO__X") == pytest.approx(
        STEADY_X, rel=1e-6
    )


def test_a_species_set_at_the_start_keeps_its_value(
    builds: dict[str, BuildResult],
) -> None:
    """A change of the simulation at its start sets a species of the drug, no zero."""
    venous = ids.pool_species(Pool.VENOUS, "drg")

    def simulations(
        experiment: PBPKSimulationExperiment,
    ) -> dict[str, Simulation | Scan]:
        start = Simulation(
            end=60, preinit_changes=CONDITION, changes=[Change(0, {venous: 0.01})]
        )
        equilibrated = experiment.equilibrated(start, resets={"HOMEO__V": 0.0})
        first, later = equilibrated.changes
        assert venous not in first.values
        assert later.values[venous] == 0.01
        return {"start": equilibrated, "plain": start}

    experiment = _run(
        _model(builds["sink"]), simulations, selections=("time", f"[{venous}]")
    )
    value = _at_start(experiment, "start", f"[{venous}]")
    assert value > 0.0
    assert value == pytest.approx(_at_start(experiment, "plain", f"[{venous}]"))


def test_an_equilibrated_simulation_exports_to_petab(
    builds: dict[str, BuildResult], tmp_path: Path
) -> None:
    """The start changes of the equilibration have a representation in PEtab.

    A condition of PEtab sets the concentration of a species whose id is its
    concentration (`hasOnlySubstanceUnits` false) and the amount of an amount
    species, so the zeros of the drug are written as `[sid]` for a
    concentration and the sinks (amounts) by their ids; the exporter of
    sbmlsim refused the zeros as amounts of concentration species once.
    """
    from sbmlsim.data import DataSet
    from sbmlsim.fit import FitData, FitMapping, FitMappingCollection, FitParameter
    from sbmlsim.fit.optimization import OptimizationProblem
    from sbmlsim.fit.petab_v2 import to_petab

    from pkpdsim.fit import default_fit_settings
    from pkpdsim.metadata import Dosing, Health, PKPDMappingMetaData, Route, Tissue

    path = builds["sink"].flat_path
    model = _model(builds["sink"])
    zeros = drug_free_zeros(path)
    assert set(zeros) & _concentrations(path)
    captured: dict[str, Simulation | Scan] = {}

    class Exported(PBPKSimulationExperiment):
        pkpd_model = model
        selections = ("time", "Carm_drg")

        def datasets(self) -> dict[str, DataSet]:
            df = pd.DataFrame(
                {
                    "time": [60.0, 120.0, 240.0],
                    "time_unit": "min",
                    "drg": [0.02, 0.01, 0.005],
                    "drg_unit": "mM",
                    "count": 6,
                }
            )
            return {"data_drg": DataSet.from_df(df, ureg=self.ureg)}

        def simulations(self) -> dict[str, Simulation | Scan]:
            captured["dosed"] = self.equilibrated(_dosed(), resets={"HOMEO__V": 0.0})
            return dict(captured)

        def fit_mappings(self) -> dict[str, FitMapping]:
            return {
                "fm_drg": FitMapping(
                    self,
                    reference=FitData(
                        self, dataset="data_drg", xid="time", yid="drg", count="count"
                    ),
                    observable=FitData(
                        self, task="task_dosed", xid="time", yid="Carm_drg"
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.IV,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                    ),
                )
            }

    problem = OptimizationProblem(
        opid="exported",
        mapping_collections=[FitMappingCollection(Exported, mappings=["fm_drg"])],
        fit_parameters=[
            FitParameter(
                pid="Kp_drg",
                start_value=2.0,
                lower_bound=0.5,
                upper_bound=5.0,
                unit="dimensionless",
            )
        ],
        base_path=path.parent,
        data_path=path.parent,
    )
    problem.initialize(default_fit_settings(absolute_tolerance=1e-10))
    yaml_file = to_petab(problem, tmp_path / "petab")
    assert yaml_file.exists()
    simulation = captured["dosed"]
    assert isinstance(simulation, Simulation)
    start = dict(simulation.changes[0].values)
    for target in [*_zero_targets(path), "Aurine_drg", "Afeces_drg"]:
        assert start[target] == 0.0, target


# --- the equilibration outside an experiment (ruling R23) -------------------------


def _choose(scan: bool) -> Simulation | Scan:
    """A simulation or a scan, typed as the union a caller with both holds."""
    simulation = _dosed()
    if scan:
        return Scan(simulation, [Dimension("dim_f", values={"HOMEO__f": [0.5, 2.0]})])
    return simulation


def test_the_equilibration_of_a_union(builds: dict[str, BuildResult]) -> None:
    """A caller with `Simulation | Scan` needs no match: the third overload."""

    def simulations(
        experiment: PBPKSimulationExperiment,
    ) -> dict[str, Simulation | Scan]:
        both = {
            key: experiment.equilibrated(_choose(scan), resets={"HOMEO__V": 0.0})
            for key, scan in (("simulation", False), ("scan", True))
        }
        assert_type(both["scan"], Simulation | Scan)
        assert isinstance(both["simulation"], Simulation)
        assert isinstance(both["scan"], Scan)
        return dict(both)

    experiment = _run(_model(builds["sink"]), simulations)
    assert _at_start(experiment, "simulation", "HOMEO__X") == pytest.approx(
        STEADY_X, rel=1e-6
    )


def test_a_simulation_of_a_model_file_is_equilibrated(
    builds: dict[str, BuildResult],
) -> None:
    """`equilibrated_simulation` is the equilibration of an experiment, from the file."""
    path = builds["sink"].flat_path
    captured: dict[str, Simulation | Scan] = {}

    def simulations(
        experiment: PBPKSimulationExperiment,
    ) -> dict[str, Simulation | Scan]:
        captured["experiment"] = experiment.equilibrated(
            _dosed(), resets={"HOMEO__V": 0.0}
        )
        return dict(captured)

    _run(_model(builds["sink"]), simulations)
    outside = equilibrated_simulation(_dosed(), path, resets={"HOMEO__V": 0.0})
    inside = captured["experiment"]
    assert isinstance(inside, Simulation)
    assert dict(outside.preinit_changes) == dict(inside.preinit_changes)
    assert dict(outside.changes[0].values) == dict(inside.changes[0].values)
    assert outside.presimulation == inside.presimulation
    result = Simulator(
        n_workers=1, absolute_tolerance=1e-10, relative_tolerance=1e-10
    ).run(path, outside)
    from pkpdsim.scan import start_values

    assert float(start_values(result, "HOMEO__X").item()) == pytest.approx(
        STEADY_X, rel=1e-6
    )
    with pytest.raises(ValueError, match=r"of the model 'sink_flat.xml'"):
        equilibrated_simulation(_dosed(), path)


def test_a_sensitivity_group_runs_from_the_equilibrated_state(
    builds: dict[str, BuildResult], tmp_path: Path
) -> None:
    """Every point of the design of an equilibrated group starts at its steady state.

    The group without the equilibration starts at the healthy steady state
    of the model (X 100), the equilibrated one at that of its condition (X
    50), whatever the EC50 of the drug, which the design varies.
    """
    from sbmlsim.simulation import Formula

    from pkpdsim.sensitivity import (
        ParameterType,
        SamplingAnalysis,
        SensitivityGroup,
        SensitivityOutput,
        SensitivityParameter,
        SensitivityRunner,
    )

    placebo = Simulation(end=60, preinit_changes=CONDITION)
    groups = [
        SensitivityGroup("plain", "Plain", placebo),
        SensitivityGroup(
            "equilibrated",
            "Equilibrated",
            placebo,
            equilibrate=True,
            resets={"HOMEO__V": 0.0},
        ),
    ]
    runner = SensitivityRunner(
        model=dataclasses.replace(_model(builds["sink"]), results_path=tmp_path),
        parameters=lambda: [
            SensitivityParameter(
                uid="HOMEO__EC50",
                name="EC50",
                value=1e-3,
                lower_bound=0.9e-3,
                upper_bound=1.1e-3,
                unit="mM",
                type=ParameterType.NA,
            )
        ],
        outputs=[SensitivityOutput("X_start", "X at the start")],
        observables=[Formula("X_start", "at(HOMEO__X, 0)", unit="dimensionless")],
        analyses={"sampling": SamplingAnalysis(groups=groups, n=8)},
        default=("sampling",),
    )
    path = runner.run(
        ["sampling"],
        simulator=SimulatorSettings(
            absolute_tolerance=1e-10, relative_tolerance=1e-10, n_workers=1
        ),
    )
    starts = {
        group.uid: pd.read_csv(path / "sampling" / group.uid / "samples.tsv", sep="\t")[
            "X_start"
        ].to_numpy()
        for group in groups
    }
    assert np.allclose(starts["plain"], 100.0, rtol=1e-9)
    assert np.allclose(starts["equilibrated"], STEADY_X, rtol=1e-6)


def test_a_sensitivity_group_resets_only_when_it_equilibrates() -> None:
    from pkpdsim.sensitivity import SensitivityGroup

    with pytest.raises(ValueError, match="has resets but no equilibration"):
        SensitivityGroup("g", "G", _dosed(), resets={"HOMEO__V": 0.0})


def test_a_reset_of_a_concentration_species_is_set_as_its_concentration(
    builds: dict[str, BuildResult],
) -> None:
    """A reset by the id of a concentration species is set as `[sid]`, once.

    A condition of PEtab sets the concentration of such a species, as the
    zeros of the drug are set; the reset of a species which is also a zero
    replaces the zero and is not set twice.
    """
    path = builds["sink"].flat_path
    venous = ids.pool_species(Pool.VENOUS, "drg")
    assert venous in _concentrations(path)
    assert venous in drug_free_zeros(path)
    for target in (venous, f"[{venous}]"):
        simulation = equilibrated_simulation(
            _dosed(), path, resets={"HOMEO__V": 0.0, target: 0.002}
        )
        values = simulation.changes[0].values
        assert values[f"[{venous}]"] == 0.002
        assert venous not in values
        assert values["HOMEO__V"] == 0.0


def test_a_build_is_parsed_and_compiled_once_for_both_checks(
    builds: dict[str, BuildResult], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The equilibration and the validity of the kidney outflow share one model.

    A copy of the build, new to the caches: its growing states and the
    flows of its kidney are read from one parse and one compilation.
    """
    import shutil

    from pkpdsim.pbpk import flat_model as flat_module
    from pkpdsim.pbpk.equilibration import growing_states
    from pkpdsim.pbpk.kidney_validity import kidney_outflow

    path = tmp_path / "copy_flat.xml"
    shutil.copyfile(builds["sink"].flat_path, path)
    loaded: list[Path] = []
    parsed: list[str] = []
    compile_model = flat_module.RoadrunnerSBMLModel
    read = flat_module.libsbml.readSBMLFromFile

    def counted_compile(**kwargs: object) -> object:
        loaded.append(Path(str(kwargs["source"])))
        return compile_model(**kwargs)  # ty: ignore[invalid-argument-type]

    def counted_read(name: str) -> object:
        parsed.append(name)
        return read(name)

    monkeypatch.setattr(flat_module, "RoadrunnerSBMLModel", counted_compile)
    monkeypatch.setattr(flat_module.libsbml, "readSBMLFromFile", counted_read)
    growing = growing_states(path, CONDITION, SteadyState())
    outflows = kidney_outflow(path, {ids.F_RENAL_FUNCTION: 0.5})
    drug_free_zeros(path)
    assert growing == ("HOMEO__V",)
    assert len(outflows) == 2
    assert loaded == [path]
    assert parsed == [str(path)]
    selections = flat_module.flat_model(path).selections
    assert {"HOMEO__V", ids.GFR, ids.HCT} <= set(selections)


def _custom(equation: str, rate: str, *, baseline: bool = True) -> PBPKDefinition:
    """`end` (with a baseline, or not) and `pre` (without) with a custom reaction in muscle."""
    return PBPKDefinition(
        sid="custom_direction",
        substances=[
            Substance(
                "end",
                name="endogenous substance",
                mr=113.0,
                distribution=Kp(0.6, ftissue=6.0),
                excretion=Excretion(urine=True),
                baseline=0.02 if baseline else None,
            ),
            Substance(
                "pre",
                name="its precursor",
                mr=100.0,
                distribution=Kp(0.6, ftissue=6.0),
                excretion=Excretion(urine=True),
            ),
        ],
        chemistry=[CustomReaction("CONV", Organ.MUSCLE, Part.TISSUE, equation, rate)],
        parameters=[ModelParameter("k_conv", 1e-3, U.l_per_min)],
        routes=[Route.iv("pre")],
    )


def test_a_custom_reaction_runs_backwards_only_when_reversible(
    tmp_path: Path,
) -> None:
    """The direction is the reversible flag: `A <-> B` feeds A from a started B.

    The exchange `k*(Cti_pre - Cti_end)` runs backwards from the endogenous
    `end` into `pre`. Written reversible, `pre` holds a value and is no zero;
    written `->` the analysis takes the reaction one way and zeroes `pre`,
    the convention the docs state.
    """
    rate = "k_conv * (Cti_pre - Cti_end)"
    both = build(_custom("Cti_pre <-> Cti_end", rate), tmp_path / "both").flat_path
    assert not set(drug_free_zeros(both)) & _species(both, "pre", "end")
    one = build(_custom("Cti_pre -> Cti_end", rate), tmp_path / "one").flat_path
    assert set(drug_free_zeros(one)) == _species(one, "pre")


def test_a_product_of_a_source_keeps_its_value(tmp_path: Path) -> None:
    """A substance without a baseline which a source (no dose) produces is no zero.

    `end` has no baseline and is produced in muscle by the source `-> Cti_end`
    (a production, not a dose): it starts at 0, reaches its steady state
    during the equilibration and keeps it at the start, while `pre`, which
    only its dose reaches, starts at exactly 0.
    """
    d = dataclasses.replace(
        _custom("-> Cti_end", "k_conv * 1 mM", baseline=False),
        parameters=[ModelParameter("k_conv", 1e-3, U.l_per_min)],
    )
    result = build(d, tmp_path)
    path = result.flat_path
    zeros = set(drug_free_zeros(path))
    assert zeros == _species(path, "pre")

    def simulations(
        experiment: PBPKSimulationExperiment,
    ) -> dict[str, Simulation | Scan]:
        return {"production": experiment.equilibrated(Simulation(end=60))}

    run = _run(
        _model(result),
        simulations,
        selections=("time", ids.arm("end"), ids.arm("pre")),
    )
    end = cast("TimecourseResult", run.results["task_production"])
    arm = np.asarray(end[ids.arm("end")]).ravel()
    assert arm[0] > 0.0
    assert arm[-1] == pytest.approx(arm[0], rel=1e-6)
    assert _at_start(run, "production", ids.arm("pre")) == 0.0
