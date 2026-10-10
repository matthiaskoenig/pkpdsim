"""The fit settings and the helpers of the mapping collections."""

import pytest
from sbmlsim.experiment import ExperimentRunner
from sbmlsim.fit import FitMapping, FitSettings
from sbmlsim.fit.objects import MappingMetaData
from sbmlsim.fit.options import (
    LossFunctionType,
    ResidualType,
    WeightingCurvesType,
    WeightingPointsType,
)
from sbmlsim.model.tolerances import AbsoluteTolerance

from pkpdsim.fit import default_fit_settings, mapping_metadata, observable_id
from pkpdsim_toy import DATA_PATHS, TOY_PATH
from pkpdsim_toy.experiments.decay import ToyDecay, ToyMappingMetaData


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
    metadata = mapping_metadata(fit_mapping, ToyMappingMetaData)
    assert metadata.route == "iv"


def test_other_metadata_is_an_error(fit_mapping: FitMapping) -> None:
    class OtherMetaData(MappingMetaData):
        pass

    with pytest.raises(TypeError, match="OtherMetaData required on"):
        mapping_metadata(fit_mapping, OtherMetaData)


def test_the_observable_of_a_mapping(fit_mapping: FitMapping) -> None:
    assert observable_id(fit_mapping) == "A"
