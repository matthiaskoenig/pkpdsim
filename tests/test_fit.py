"""The fit settings and the helpers of the mapping collections."""

import copy
import dataclasses
import multiprocessing
import os
from concurrent.futures import ProcessPoolExecutor

import pytest
from sbmlsim.experiment import ExperimentRunner
from sbmlsim.fit import FitMapping, FitSettings
from sbmlsim.fit.objects import MappingMetaData
from sbmlsim.fit.optimization import OptimizationProblem
from sbmlsim.fit.options import (
    DataAtChangeType,
    LossFunctionType,
    ResidualType,
    WeightingCurvesType,
    WeightingPointsType,
)
from sbmlsim.model.tolerances import AbsoluteTolerance
from sbmlsim.project import SimulatorSettings

from pkpdsim.fit import (
    check_sampling_sites,
    default_fit_settings,
    mapping_metadata,
    metadata_filter,
)
from pkpdsim.metadata import (
    Coadministration,
    Dosing,
    PKPDMappingMetaData,
    Route,
    SamplingSite,
    Tissue,
)
from pkpdsim_toy import DATA_PATHS, MODEL, TOY_PATH
from pkpdsim_toy.experiments.decay import ToyDecay


def test_the_default_settings_are_those_of_the_models() -> None:
    """The settings 40 models of pkdb_models share apart from the tolerances."""
    tolerance = AbsoluteTolerance(amount=1e-12, concentration=1e-12, other=1e-6)
    assert default_fit_settings(absolute_tolerance=tolerance) == FitSettings(
        residual=ResidualType.NORMALIZED,
        loss_function=LossFunctionType.LINEAR,
        weighting_curves=(WeightingCurvesType.MAPPING, WeightingCurvesType.POINTS),
        weighting_points=WeightingPointsType.ERROR_WEIGHTING,
        variable_step_size=True,
        relative_tolerance=1e-6,
        absolute_tolerance=tolerance,
    )
    assert (
        default_fit_settings(1e-6, relative_tolerance=1e-8).relative_tolerance == 1e-8
    )


def test_the_data_at_a_change_is_compared_with_the_state_before_it() -> None:
    """sbmlsim's default, written out so that a change of it does not change the fits."""
    assert default_fit_settings(1e-8).data_at_change is DataAtChangeType.BEFORE


@pytest.fixture(scope="module")
def fit_mapping() -> FitMapping:
    runner = ExperimentRunner(
        experiment_classes=[ToyDecay],
        base_path=TOY_PATH,
        data_path=DATA_PATHS,
        simulator=None,
    )
    return runner.experiments["ToyDecay"]._fit_mappings["fm_decay_A"]


def test_the_metadata_of_a_mapping(fit_mapping: FitMapping) -> None:
    metadata = mapping_metadata(fit_mapping, PKPDMappingMetaData)
    assert metadata.route == "iv"


def test_other_metadata_is_an_error(fit_mapping: FitMapping) -> None:
    class OtherMetaData(MappingMetaData):
        pass

    with pytest.raises(TypeError, match="OtherMetaData required on"):
        mapping_metadata(fit_mapping, OtherMetaData)


def with_metadata(fit_mapping: FitMapping, **changes: object) -> FitMapping:
    """A copy of the mapping whose metadata differs in some fields."""
    metadata = mapping_metadata(fit_mapping, PKPDMappingMetaData)
    other = copy.copy(fit_mapping)
    other.metadata = dataclasses.replace(metadata, **changes)
    return other


def test_metadata_filter_defaults_to_the_shared_metadata() -> None:
    f = metadata_filter(route={Route.IV})
    assert f.metadata_type is PKPDMappingMetaData
    assert repr(
        metadata_filter(
            route={Route.PO, Route.IV},
            coadministration={Coadministration.NONE},
            excluded={"dosing": {Dosing.MULTIPLE}},
        )
    ) == (
        "MetadataFilter(PKPDMappingMetaData: route in [IV, PO] "
        "and coadministration in [NONE] and dosing not in [MULTIPLE])"
    )


def test_metadata_filter_takes_a_local_subclass() -> None:
    @dataclasses.dataclass(kw_only=True)
    class LocalMetaData(PKPDMappingMetaData):
        diplotype: str | None = None

    f = metadata_filter(LocalMetaData, diplotype={"*1/*1"})
    assert f.metadata_type is LocalMetaData


def test_mapping_metadata_defaults_to_the_shared_metadata(
    fit_mapping: FitMapping,
) -> None:
    assert mapping_metadata(fit_mapping).route == "iv"


def test_check_sampling_sites_names_every_mismatch(fit_mapping: FitMapping) -> None:
    assert not check_sampling_sites([fit_mapping], {"A": SamplingSite.ARM_VEIN})
    wrong = check_sampling_sites([fit_mapping], {"A": SamplingSite.VENOUS})
    assert len(wrong) == 1
    assert "'A'" in wrong[0]
    assert "ARM_VEIN" in wrong[0]
    assert "VENOUS" in wrong[0]
    other = with_metadata(fit_mapping, sampling_site=SamplingSite.ARTERIAL)
    problems = check_sampling_sites([fit_mapping, other], {"A": SamplingSite.ARM_VEIN})
    assert len(problems) == 1
    unknown = check_sampling_sites([fit_mapping, other], {"B": SamplingSite.ARM_VEIN})
    assert len(unknown) == 2
    assert all("no sampling site" in problem for problem in unknown)


def test_a_filter_on_the_sampling_site_compares_the_resolved_site(
    fit_mapping: FitMapping,
) -> None:
    arm = metadata_filter(sampling_site={SamplingSite.ARM_VEIN})
    not_applicable = metadata_filter(sampling_site={SamplingSite.NOT_APPLICABLE})
    assert arm("k", fit_mapping)
    assert not not_applicable("k", fit_mapping)
    urine = with_metadata(fit_mapping, tissue=Tissue.URINE)
    assert not_applicable("k", urine)
    assert not arm("k", urine)


def test_the_filter_value_of_the_sampling_site_is_the_resolved_site(
    fit_mapping: FitMapping,
) -> None:
    metadata = mapping_metadata(fit_mapping)
    assert metadata.sampling_site is None
    assert metadata.filter_value("sampling_site") is SamplingSite.ARM_VEIN
    assert metadata.filter_value("route") is Route.IV


#: the environment variable of sbmlsim of the level of a model loaded without one
JIT_VARIABLE = "SBMLSIM_JIT_OPTIMIZATION"


def environment_level() -> str | None:
    """Get the level of the environment of a process."""
    return os.environ.get(JIT_VARIABLE)


def levels_in_worker(problem: OptimizationProblem, settings: FitSettings) -> tuple:
    """Get the levels a worker of a pool sees: of the problem, its environment, its models."""
    problem.initialize(settings)
    assert problem.runner is not None
    return (
        str(problem.jit_optimization),
        os.environ.get(JIT_VARIABLE),
        str(problem.runner.jit_optimization),
    )


# windows has no forkserver
@pytest.mark.parametrize(
    "start_method",
    [
        m
        for m in ("forkserver", "spawn")
        if m in multiprocessing.get_all_start_methods()
    ],
)
def test_a_worker_loads_the_models_at_the_level_of_the_model(
    monkeypatch: pytest.MonkeyPatch, start_method: str
) -> None:
    """The problem carries the level into a worker, whatever the environment of the pool.

    A forkserver which started before, without the level, starts the workers
    with its own environment, so the level cannot come from the environment
    of the parent; the problem is pickled into every task of a pool with it.
    """
    monkeypatch.delenv(JIT_VARIABLE, raising=False)
    context = multiprocessing.get_context(start_method)
    with ProcessPoolExecutor(1, mp_context=context) as pool:
        assert pool.submit(environment_level).result() is None
    leveled = MODEL.replace(simulator=SimulatorSettings(jit_optimization="O0"))
    definition = leveled.load_fit_definitions()["DECAY"]
    problem = definition.problem(opid="levels")
    assert problem.jit_optimization == "O0"
    with ProcessPoolExecutor(1, mp_context=context) as pool:
        levels = pool.submit(levels_in_worker, problem, definition.settings).result()
    assert levels == ("O0", None, "O0")
    # nothing of the level is left in the process
    assert JIT_VARIABLE not in os.environ
    # a model without a level gives its definitions without one
    plain = MODEL.load_fit_definitions()["DECAY"]
    assert plain.jit_optimization is None
    assert plain.problem(opid="default").jit_optimization == "O2"
