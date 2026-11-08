"""The experiments of a machine perfusion: the condition and the medium changes.

A medium change replaces a fraction of the reservoir with fresh perfusate, a
`Change` of sbmlsim at its time; the tests run the changes on a built rat
liver (and kidney) circuit with sbmlsim's simulator, as an experiment does.
"""

import math
from pathlib import Path
from typing import ClassVar, cast

import libsbml
import numpy as np
import pytest
from sbmlsim import Q
from sbmlsim.experiment import ExperimentRunner
from sbmlsim.fit.options import DataAtChangeType
from sbmlsim.result import TimecourseResult
from sbmlsim.simulation import Change, Dimension, Scan, Simulation
from sbmlsim.simulator import Simulator
from sbmlsim.units import Quantity, ureg

from pkpdsim.fit import default_fit_settings
from pkpdsim.model import PKPDModel
from pkpdsim.perfusion import circuit, ids
from pkpdsim.perfusion.builder import build, compose, model_ids
from pkpdsim.perfusion.definition import PerfusionDefinition
from pkpdsim.perfusion.experiment import (
    PerfusionSimulationExperiment,
    medium_change,
    observable_labels,
    observable_units,
    perfusion_changes,
)

from .conftest import kidney_definition, liver_definition

TRACER = "trc"
#: the time of the first and the second medium change (min)
FIRST, SECOND = 60.0, 120.0


@pytest.fixture(scope="module")
def liver() -> PerfusionDefinition:
    """The rat liver circuit with two zones and a tracer."""
    return liver_definition(zones=2, tracers=(TRACER,), sid="mp_liver")


@pytest.fixture(scope="module")
def liver_path(
    liver: PerfusionDefinition, tmp_path_factory: pytest.TempPathFactory
) -> Path:
    return build(liver, tmp_path_factory.mktemp("liver"))


@pytest.fixture(scope="module")
def kidney() -> PerfusionDefinition:
    return kidney_definition(tracers=(TRACER,), sid="mp_kidney")


@pytest.fixture(scope="module")
def kidney_path(
    kidney: PerfusionDefinition, tmp_path_factory: pytest.TempPathFactory
) -> Path:
    return build(kidney, tmp_path_factory.mktemp("kidney"))


def _reservoir(definition: PerfusionDefinition) -> list[str]:
    """The quantities of the reservoir."""
    return list(circuit.quantities(definition.perfusate, definition.released))


def _simulate(
    path: Path,
    definition: PerfusionDefinition,
    changes: list[Change],
    times: list[float],
    extra: tuple[str, ...] = (),
) -> TimecourseResult:
    """Simulate the circuit for 3 h with the tracer at 1 mM after the start."""
    simulator = Simulator(
        n_workers=1, absolute_tolerance=1e-14, relative_tolerance=1e-10
    )
    model = simulator.load(path)
    quantities = _reservoir(definition)
    model.selections = [
        "time",
        *(f"[{ids.species(q, ids.RESERVOIR)}]" for q in quantities),
        *(ids.fresh(q) for q in quantities),
        *extra,
    ]
    start = Change(0.0, {f"[{ids.species(TRACER, ids.RESERVOIR)}]": 1.0})
    simulation = Simulation(
        end=180.0, time_unit="min", changes=[start, *changes], times=times
    )
    return simulator.simulate(model, simulation)


def _column(result: TimecourseResult, sid: str) -> np.ndarray:
    return np.asarray(result[sid], dtype=float).ravel()


def _rows(result: TimecourseResult, time: float) -> np.ndarray:
    """The rows of the output at a time."""
    return np.flatnonzero(np.isclose(_column(result, "time"), time, rtol=0, atol=1e-12))


def _distances(
    result: TimecourseResult, definition: PerfusionDefinition, row: int
) -> dict[str, float]:
    """The distance of every reservoir species to the fresh perfusate at a row."""
    return {
        q: _column(result, f"[{ids.species(q, ids.RESERVOIR)}]")[row]
        - _column(result, ids.fresh(q))[row]
        for q in _reservoir(definition)
    }


class _Liver(PerfusionSimulationExperiment):
    """An experiment of the circuit of the tests, bound in `experiment`."""


@pytest.fixture(scope="module")
def experiment(
    liver: PerfusionDefinition, liver_path: Path
) -> type[PerfusionSimulationExperiment]:
    """An experiment class of the liver circuit, bound to a model object."""
    model = PKPDModel(
        name="mp_liver",
        package="mp_liver",
        title="Perfused rat liver",
        path=liver_path.parent,
        model_path=liver_path,
        factory=".factory:create_models",
    )

    class Liver(PerfusionSimulationExperiment):
        pkpd_model = model
        definition = liver

    return Liver


def test_a_medium_change_halves_the_distance_to_the_fresh_perfusate(
    liver: PerfusionDefinition,
    liver_path: Path,
    experiment: type[PerfusionSimulationExperiment],
) -> None:
    changes = experiment.medium_changes([FIRST], 0.5)
    result = _simulate(liver_path, liver, changes, [0.0, 30.0, FIRST, 90.0])
    before, after = _rows(result, FIRST)
    d_before = _distances(result, liver, before)
    d_after = _distances(result, liver, after)
    assert 0.5 < d_before[TRACER] < 1.0, "the tracer mixes into the circuit"
    assert abs(d_before["lac"]) > 1e-4, "the liver changes the lactate"
    for q in d_before:
        assert d_after[q] == pytest.approx(0.5 * d_before[q], rel=1e-9, abs=1e-15), q


def test_the_value_at_the_time_of_a_change_is_the_value_before_it(
    liver: PerfusionDefinition,
    liver_path: Path,
    experiment: type[PerfusionSimulationExperiment],
) -> None:
    """sbmlsim's order: an output time at a change has two rows, before and after.

    The first row at the time is the state before the change (the state of
    the simulation without the change), the second the state after it, and
    a fit compares a data point at the time with the state before.
    """
    times = [0.0, 30.0, FIRST, 90.0]
    changed = _simulate(
        liver_path, liver, experiment.medium_changes([FIRST], 0.5), times
    )
    plain = _simulate(liver_path, liver, [], times)
    assert list(_column(changed, "time")) == [0.0, 30.0, FIRST, FIRST, 90.0]
    assert list(_column(plain, "time")) == times
    before, after = _rows(changed, FIRST)
    (unchanged,) = _rows(plain, FIRST)
    assert (before, after) == (2, 3)
    for q in _reservoir(liver):
        sid = f"[{ids.species(q, ids.RESERVOIR)}]"
        assert _column(changed, sid)[before] == pytest.approx(
            _column(plain, sid)[unchanged], rel=1e-12, abs=1e-18
        ), q
    lactate = "[lac_res]"
    assert _column(changed, lactate)[after] != pytest.approx(
        _column(changed, lactate)[before], rel=1e-6
    )
    assert default_fit_settings(1e-10).data_at_change is DataAtChangeType.BEFORE


def test_two_medium_changes_compose(
    liver: PerfusionDefinition,
    liver_path: Path,
    experiment: type[PerfusionSimulationExperiment],
) -> None:
    changes = experiment.medium_changes([FIRST, SECOND], 0.5)
    assert [c.times for c in changes] == [
        (ureg.Quantity(FIRST, "min"),),
        (ureg.Quantity(SECOND, "min"),),
    ]
    result = _simulate(liver_path, liver, changes, [0.0, FIRST, 90.0, SECOND, 150.0])
    for time in (FIRST, SECOND):
        before, after = _rows(result, time)
        d_before = _distances(result, liver, before)
        d_after = _distances(result, liver, after)
        for q in d_before:
            assert d_after[q] == pytest.approx(
                0.5 * d_before[q], rel=1e-9, abs=1e-15
            ), (time, q)
    # the tracer is inert: the second change halves what the first left
    (_, after_first) = _rows(result, FIRST)
    (before_second, after_second) = _rows(result, SECOND)
    trc = f"[{ids.species(TRACER, ids.RESERVOIR)}]"
    assert _column(result, trc)[after_second] == pytest.approx(
        0.5 * _column(result, trc)[before_second], rel=1e-9
    )
    assert _column(result, trc)[after_second] < _column(result, trc)[after_first]


def test_a_change_in_hours_is_at_its_time(
    liver: PerfusionDefinition,
    liver_path: Path,
    experiment: type[PerfusionSimulationExperiment],
) -> None:
    (change,) = experiment.medium_changes([1.0], 0.5, time_unit="hr")
    assert change.times == (ureg.Quantity(1.0, "hr"),)
    result = _simulate(liver_path, liver, [change], [0.0, FIRST, 90.0])
    assert len(_rows(result, FIRST)) == 2


def test_every_species_of_the_reservoir_is_changed(
    liver: PerfusionDefinition,
    liver_path: Path,
    kidney: PerfusionDefinition,
    kidney_path: Path,
) -> None:
    for definition, path in ((liver, liver_path), (kidney, kidney_path)):
        model = libsbml.readSBMLFromFile(str(path)).getModel()
        reservoir = {
            s.getId()
            for s in model.getListOfSpecies()
            if s.getCompartment() == ids.volume(ids.RESERVOIR)
        }
        values = medium_change(definition, 0.5)
        assert {t.strip("[]") for t in values} == reservoir
        q = "glc"
        assert values[f"[{ids.species(q, ids.RESERVOIR)}]"] == (
            "(1 - 0.5)*[glc_res] + 0.5*glc_fresh"
        )
        for formula in values.values():
            fresh = formula.rsplit("*", 1)[1]
            assert model.getParameter(fresh) is not None, fresh


@pytest.mark.parametrize("fraction", [0.0, -0.1, 1.5, math.nan, math.inf])
def test_the_fraction_of_a_change_is_checked(
    liver: PerfusionDefinition, fraction: float
) -> None:
    with pytest.raises(ValueError, match="fraction"):
        medium_change(liver, fraction)


def test_a_change_of_the_whole_reservoir_is_fresh_perfusate(
    liver: PerfusionDefinition,
    liver_path: Path,
    experiment: type[PerfusionSimulationExperiment],
) -> None:
    result = _simulate(
        liver_path, liver, experiment.medium_changes([FIRST], 1.0), [0.0, FIRST]
    )
    _, after = _rows(result, FIRST)
    for q, d in _distances(result, liver, after).items():
        assert d == pytest.approx(0.0, abs=1e-15), q


def test_a_change_keeps_the_volume_of_the_reservoir_of_a_kidney(
    kidney: PerfusionDefinition, kidney_path: Path
) -> None:
    """The reservoir of a kidney loses the urine; a change keeps its volume."""
    changes = [Change(FIRST, medium_change(kidney, 0.5))]
    v_res = ids.volume(ids.RESERVOIR)
    result = _simulate(kidney_path, kidney, changes, [0.0, FIRST, 90.0], extra=(v_res,))
    before, after = _rows(result, FIRST)
    volume = _column(result, v_res)
    assert volume[before] < volume[0], "the urine leaves the reservoir"
    assert volume[after] == volume[before]
    d_before = _distances(result, kidney, before)
    d_after = _distances(result, kidney, after)
    for q in d_before:
        assert d_after[q] == pytest.approx(0.5 * d_before[q], rel=1e-9, abs=1e-15), q


def test_the_perfusion_changes_set_the_condition(
    liver: PerfusionDefinition, liver_path: Path, kidney: PerfusionDefinition
) -> None:
    changes = perfusion_changes(liver, Q(2.0, "ml/min"), Q(15.0, "g"), Q(20.0, "degC"))
    assert set(changes) == {ids.Q_PUMP, ids.mass("li"), ids.TEMPERATURE}
    assert changes[ids.TEMPERATURE] == pytest.approx(20.0)
    kelvin = perfusion_changes(liver, Q(2.0, "ml/min"), Q(15.0, "g"), Q(293.15, "K"))
    assert kelvin[ids.TEMPERATURE] == pytest.approx(20.0)
    assert set(perfusion_changes(liver, Q(2.0, "ml/min"), Q(15.0, "g"))) == {
        ids.Q_PUMP,
        ids.mass("li"),
    }
    assert ids.mass("ki") in perfusion_changes(kidney, Q(1.0, "ml/min"), Q(1.0, "g"))

    simulator = Simulator(n_workers=1)
    model = simulator.load(liver_path)
    model.selections = ["time", ids.Q_PUMP, ids.mass("li"), ids.TEMPERATURE, ids.F_T]
    result = simulator.simulate(
        model, Simulation(end=1.0, steps=1, preinit_changes=changes)
    )
    assert _column(result, ids.Q_PUMP)[0] == pytest.approx(0.002)
    assert _column(result, ids.mass("li"))[0] == pytest.approx(0.015)
    assert _column(result, ids.TEMPERATURE)[0] == pytest.approx(20.0)
    assert 0.0 < _column(result, ids.F_T)[0] < 1.0


@pytest.mark.parametrize(
    ("flow", "mass", "temperature"),
    [
        (Q(-1.0, "ml/min"), Q(15.0, "g"), None),
        (Q(1.0, "ml/min"), Q(0.0, "g"), None),
        (Q(1.0, "ml"), Q(15.0, "g"), None),
        (Q(1.0, "ml/min"), Q(15.0, "ml"), None),
        (Q(1.0, "ml/min"), Q(15.0, "g"), Q(-300.0, "degC")),
    ],
)
def test_the_perfusion_changes_are_checked(
    liver: PerfusionDefinition,
    flow: Quantity,
    mass: Quantity,
    temperature: Quantity | None,
) -> None:
    with pytest.raises(ValueError):
        perfusion_changes(liver, flow, mass, temperature)


def test_the_observables_have_labels_and_units(
    liver: PerfusionDefinition, kidney: PerfusionDefinition
) -> None:
    for definition in (liver, kidney):
        labels = observable_labels(definition)
        units = observable_units(definition)
        assert set(labels) == set(units)
        model = compose(definition)
        assert {sid.strip("[]") for sid in labels} <= model_ids(model)
        species = {sp.sid for sp in model.species}
        for sid in labels:
            # a species by its concentration, an assignment rule by its id
            assert (sid.strip("[]") in species) == sid.startswith("["), sid
        for sid, unit in units.items():
            ureg(unit)  # a unit of sbmlsim
            assert labels[sid], sid
    labels, units = observable_labels(liver), observable_units(liver)
    assert labels["[glc_res]"] == "glucose (reservoir)"
    assert labels["glc_pp"] == "glucose (periportal)"
    assert labels["pO2_out"] == "pO2 (outflow)"
    assert labels["alt_act_res"] == "ALT activity (reservoir)"
    assert labels["[atp_z2_cell]"] == "ATP (cells, zone 2)"
    assert labels["atp_pc_cell"] == "ATP (cells, pericentral)"
    assert labels["f_alive_pc"] == "living fraction (pericentral)"
    assert units["alt_act_res"] == "µkat/l"
    assert units["bile"] == "µl"
    assert units[f"[{TRACER}_res]"] == "mM"
    assert "urine" in observable_units(kidney)
    assert "[glc_urine]" in observable_units(kidney)
    assert "bile" not in observable_units(kidney)


def test_an_experiment_has_the_labels_of_its_definition(
    liver: PerfusionDefinition, experiment: type[PerfusionSimulationExperiment]
) -> None:
    assert dict(experiment.labels) == observable_labels(liver)
    assert dict(experiment.units) == observable_units(liver)

    class Own(experiment):  # ty: ignore[unsupported-base]
        labels: ClassVar = {"[glc_res]": "glucose"}

    assert dict(Own.labels) == {"[glc_res]": "glucose"}
    assert dict(Own.units) == observable_units(liver)


def test_an_experiment_without_a_definition_is_refused() -> None:
    with pytest.raises(TypeError, match="definition"):
        _Liver.medium_changes([FIRST], 0.5)
    with pytest.raises(TypeError, match="definition"):
        _Liver.perfusion_changes(Q(1.0, "ml/min"), Q(10.0, "g"))


def test_an_experiment_runs_the_conditions_with_their_medium_changes(
    experiment: type[PerfusionSimulationExperiment],
) -> None:
    """A run of an experiment as the experiments of a study do."""

    class Flows(experiment):  # ty: ignore[unsupported-base]
        selections: ClassVar = ("time", "[glc_res]", "alt_act_res", "pO2_out")

        def simulations(self) -> dict[str, Simulation | Scan]:
            return {
                f"q{q:g}_mc{n}": Simulation(
                    end=360.0,
                    time_unit="min",
                    preinit_changes=self.perfusion_changes(
                        Q(q, "ml/min"), Q(15.0, "g"), Q(37.0, "degC")
                    ),
                    changes=self.medium_changes([120.0, 240.0][:n], 0.5),
                    times=[0.0, 60.0, 120.0, 180.0, 240.0, 300.0, 360.0],
                )
                for q, n in ((1.0, 0), (3.0, 0), (3.0, 2))
            }

    model = experiment.pkpd_model
    runner = ExperimentRunner(
        experiment_classes=[Flows],
        base_path=model.path,
        simulator=Simulator(n_workers=1),
        on_error="raise",
    )
    (run,) = runner.experiments.values()
    run = cast("PerfusionSimulationExperiment", run)
    run.run(runner.simulator, output_path=None)
    alt = {
        key: _column(
            cast("TimecourseResult", run.results[f"task_{key}"]), "alt_act_res"
        )
        for key in ("q1_mc0", "q3_mc0", "q3_mc2")
    }
    assert len(alt["q3_mc0"]) == 7
    assert len(alt["q3_mc2"]) == 9, "two rows at each change"
    assert alt["q3_mc2"][-1] < alt["q3_mc0"][-1], "the changes dilute the ALT"
    assert np.all(np.isfinite(alt["q1_mc0"]))


def test_a_species_of_the_map_is_its_concentration(
    liver: PerfusionDefinition, liver_path: Path
) -> None:
    """The selection of a species in the labels is its concentration (mM).

    roadrunner reads the plain id of a species as its amount: the inflow
    tubing of 10 ml holds 0.065 mmole of glucose at 6.5 mM.
    """
    units = observable_units(liver)
    selections = [s for s in units if s.startswith("[") and s.endswith("_in]")]
    selections += ["[glc_z1]", "[atp_z1_cell]"]
    assert set(selections) <= set(units)
    simulator = Simulator(n_workers=1)
    model = simulator.load(liver_path)
    v_in, v_z1, v_cell = ids.volume(ids.INFLOW), ids.vascular("z1"), ids.cells("z1")
    model.selections = [
        "time",
        *selections,
        *(s.strip("[]") for s in selections),
        v_in,
        v_z1,
        v_cell,
    ]
    result = simulator.simulate(model, Simulation(end=30.0, times=[0.0, 30.0]))
    volume = {"_in": v_in, "_z1": v_z1, "_cell": v_cell}
    for sid in selections:
        amount = _column(result, sid.strip("[]"))
        v = next(
            _column(result, x) for end, x in volume.items() if sid.endswith(f"{end}]")
        )
        assert _column(result, sid) == pytest.approx(amount / v, rel=1e-12), sid
    assert _column(result, "[glc_in]")[0] == pytest.approx(6.5)
    assert _column(result, "glc_in")[0] == pytest.approx(0.065)


def test_sbmlsim_converts_a_temperature_quantity_to_kelvin(liver_path: Path) -> None:
    """The reason to give `T` as a number: `Q(37, "degC")` becomes 310.15."""
    simulator = Simulator(n_workers=1)
    model = simulator.load(liver_path)
    model.selections = ["time", ids.TEMPERATURE, ids.F_T]
    result = simulator.simulate(
        model,
        Simulation(
            end=1.0, steps=1, preinit_changes={ids.TEMPERATURE: Q(37.0, "degC")}
        ),
    )
    assert _column(result, ids.TEMPERATURE)[0] == pytest.approx(310.15)
    assert _column(result, ids.F_T)[0] > 100.0


def test_an_experiment_refuses_a_temperature_quantity(
    experiment: type[PerfusionSimulationExperiment],
) -> None:
    flow, mass = Q(3.0, "ml/min"), Q(10.0, "g")

    def runner(simulation: Simulation | Scan) -> ExperimentRunner:
        class Warm(experiment):  # ty: ignore[unsupported-base]
            def simulations(self) -> dict[str, Simulation | Scan]:
                return {"warm": simulation}

        return ExperimentRunner(
            experiment_classes=[Warm],
            base_path=experiment.pkpd_model.path,
            simulator=Simulator(n_workers=1),
            on_error="raise",
        )

    changes = experiment.perfusion_changes(flow, mass)
    refused = [
        Simulation(end=10.0, preinit_changes={**changes, "T": Q(37.0, "degC")}),
        Simulation(end=10.0, changes=[Change(5.0, {"T": Q(20.0, "degC")})]),
        Scan(
            Simulation(end=10.0, preinit_changes=changes),
            [Dimension("temperature", values={"T": Q([4.0, 37.0], "degC")})],
        ),
    ]
    for simulation in refused:
        with pytest.raises(ValueError, match="perfusion_changes"):
            runner(simulation)
    # a number in degC, as perfusion_changes gives it
    good = experiment.perfusion_changes(flow, mass, Q(37.0, "degC"))
    runner(Simulation(end=10.0, preinit_changes=good))
