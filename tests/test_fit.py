"""The fit settings and the helpers of the mapping collections."""

import copy
import dataclasses

import pytest
from sbmlsim.experiment import ExperimentRunner
from sbmlsim.fit import FitMapping, FitMappingCollection, FitSettings, MappingKind
from sbmlsim.fit.objects import MappingMetaData
from sbmlsim.fit.options import (
    DataAtChangeType,
    LossFunctionType,
    ResidualType,
    WeightingCurvesType,
    WeightingPointsType,
)
from sbmlsim.model.tolerances import AbsoluteTolerance

from pkpdsim.fit import (
    default_fit_settings,
    fit_subset,
    mapping_metadata,
    metadata_filter,
    observable_filter,
    observable_id,
    select_mappings,
)
from pkpdsim.metadata import (
    Coadministration,
    Dosing,
    Fasting,
    Health,
    PKPDMappingMetaData,
    Route,
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


def test_the_observable_of_a_mapping(fit_mapping: FitMapping) -> None:
    assert observable_id(fit_mapping) == "A"


def with_metadata(fit_mapping: FitMapping, **changes: object) -> FitMapping:
    """A copy of the mapping whose metadata differs in some fields."""
    metadata = mapping_metadata(fit_mapping, PKPDMappingMetaData)
    other = copy.copy(fit_mapping)
    other.metadata = dataclasses.replace(metadata, **changes)
    return other


def test_a_mapping_with_an_allowed_value_passes(fit_mapping: FitMapping) -> None:
    control = metadata_filter(PKPDMappingMetaData, route={Route.IV, Route.PO})
    assert control("fm_decay_A", fit_mapping)
    assert not control("fm_decay_A", with_metadata(fit_mapping, route=Route.SC))


def test_every_clause_must_hold(fit_mapping: FitMapping) -> None:
    control = metadata_filter(
        PKPDMappingMetaData, route={Route.IV}, health={Health.HEALTHY}
    )
    assert control("fm_decay_A", fit_mapping)
    assert not control("fm_decay_A", with_metadata(fit_mapping, route=Route.PO))


def test_an_excluded_value_fails(fit_mapping: FitMapping) -> None:
    single = metadata_filter(
        PKPDMappingMetaData, excluded={"dosing": {Dosing.MULTIPLE}}
    )
    assert single("fm_decay_A", fit_mapping)
    assert not single("fm_decay_A", with_metadata(fit_mapping, dosing=Dosing.MULTIPLE))


def test_a_constant_set_is_a_clause(fit_mapping: FitMapping) -> None:
    control_fasting = {Fasting.FASTED, Fasting.NR}
    fasted = metadata_filter(PKPDMappingMetaData, fasting=control_fasting)
    assert fasted("fm_decay_A", fit_mapping)
    assert not fasted("fm_decay_A", with_metadata(fit_mapping, fasting=Fasting.FED))


def test_an_unknown_field_fails_when_the_filter_is_built() -> None:
    with pytest.raises(TypeError, match=r"PKPDMappingMetaData has no field.*diplotype"):
        metadata_filter(PKPDMappingMetaData, diplotype={"*1/*1"})
    with pytest.raises(TypeError, match=r"has no field.*tisue"):
        metadata_filter(PKPDMappingMetaData, excluded={"tisue": {Tissue.URINE}})


def test_a_filter_needs_a_clause_and_values() -> None:
    with pytest.raises(ValueError, match="at least one clause"):
        metadata_filter(PKPDMappingMetaData)
    with pytest.raises(ValueError, match=r"selects nothing.*route"):
        metadata_filter(PKPDMappingMetaData, route=set())


def test_a_single_value_is_not_a_collection() -> None:
    with pytest.raises(TypeError, match=r"route.*collection"):
        metadata_filter(PKPDMappingMetaData, route=Route.IV)
    with pytest.raises(TypeError, match=r"dosing.*collection"):
        metadata_filter(PKPDMappingMetaData, excluded={"dosing": Dosing.SINGLE})
    with pytest.raises(TypeError, match=r"route.*collection"):
        metadata_filter(PKPDMappingMetaData, route=None)  # ty: ignore[invalid-argument-type]
    with pytest.raises(TypeError, match=r"observables.*collection"):
        observable_filter("Cve_hctz")


def test_an_empty_exclusion_or_observable_set_is_an_error() -> None:
    with pytest.raises(ValueError, match=r"excludes nothing.*dosing"):
        metadata_filter(PKPDMappingMetaData, excluded={"dosing": set()})
    with pytest.raises(ValueError, match=r"observables.*selects nothing"):
        observable_filter([])


def test_the_metadata_class_is_not_an_instance() -> None:
    metadata = PKPDMappingMetaData(
        tissue=Tissue.PLASMA,
        route=Route.IV,
        dosing=Dosing.SINGLE,
        health=Health.HEALTHY,
    )
    with pytest.raises(TypeError, match="class, not an instance"):
        metadata_filter(metadata, route={Route.IV})  # ty: ignore[invalid-argument-type]


def test_the_metadata_class_must_be_a_dataclass() -> None:
    class OtherMetaData:
        route: Route

    with pytest.raises(TypeError, match="OtherMetaData is not a dataclass"):
        metadata_filter(OtherMetaData, route={Route.IV})  # ty: ignore[invalid-argument-type]


def test_the_filter_checks_the_type_of_the_metadata(fit_mapping: FitMapping) -> None:
    other = dataclasses.make_dataclass(
        "Other", [("route", Route)], bases=(MappingMetaData,), kw_only=True
    )
    with pytest.raises(TypeError, match="Other required on"):
        metadata_filter(other, route={Route.IV})("fm_decay_A", fit_mapping)  # ty: ignore[invalid-argument-type]


def test_the_filter_says_what_it_selects() -> None:
    control = metadata_filter(
        PKPDMappingMetaData,
        route={Route.PO, Route.IV},
        excluded={"dosing": {Dosing.MULTIPLE}},
    )
    assert repr(control) == (
        "MetadataFilter(PKPDMappingMetaData: route in [IV, PO] "
        "and dosing not in [MULTIPLE])"
    )


def test_the_observable_filter(fit_mapping: FitMapping) -> None:
    pk_filter = observable_filter({"A", "B"})
    assert pk_filter("fm_decay_A", fit_mapping)
    assert not observable_filter({"B"})("fm_decay_A", fit_mapping)
    assert repr(pk_filter) == "ObservableFilter(A, B)"


def training(collections: dict[str, list[FitMappingCollection]]) -> list[list[str]]:
    return [
        c.mappings for c in collections["ToyDecay"] if c.kind is MappingKind.TRAINING
    ]


def test_the_mappings_are_selected_by_the_filters() -> None:
    selected = select_mappings(
        MODEL, [ToyDecay], [metadata_filter(PKPDMappingMetaData, route={Route.IV})]
    )
    assert training(selected) == [["fm_decay_A"]]
    excluded = select_mappings(
        MODEL, [ToyDecay], [metadata_filter(PKPDMappingMetaData, route={Route.PO})]
    )
    assert training(excluded) == []


def test_outliers_and_validation_are_passed_on() -> None:
    outlier = select_mappings(MODEL, [ToyDecay], outliers={"fm_decay_A"})
    assert training(outlier) == []
    with pytest.raises(ValueError, match="no_mapping"):
        select_mappings(MODEL, [ToyDecay], validation={"no_mapping"})


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


@pytest.mark.parametrize(
    ("fit_id", "subset"),
    [
        ("PK_20261016_185018__3efe3", "PK"),
        ("PK_X_20261016_185018__3efe3", "PK_X"),
        ("20261016_185018__3efe3", ""),
    ],
)
def test_the_subset_of_a_fit_id(fit_id: str, subset: str) -> None:
    if subset:
        assert fit_subset(fit_id) == subset
    else:
        with pytest.raises(ValueError, match="not the id of a fit"):
            fit_subset(fit_id)


@pytest.mark.parametrize(
    "fit_id",
    [
        "PK",
        "PK_2026_1",
        "PK_20261016_185018_3efe3",
        "../PK_20261016_185018__3efe3",
        "a/PK_20261016_185018__3efe3",
        "a\\PK_20261016_185018__3efe3",
        "..PK_20261016_185018__3efe3",
    ],
)
def test_an_invalid_fit_id_has_no_subset(fit_id: str) -> None:
    with pytest.raises(ValueError, match="not the id of a fit"):
        fit_subset(fit_id)
