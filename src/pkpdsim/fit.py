"""The fit settings, the metadata filters and the selection of the mappings.

The models of pkdb_models share their fit settings apart from the integrator
tolerances (`default_fit_settings`). Every `mapping_collections.py` checks the
type of the metadata of a mapping (`mapping_metadata`), reads the observable of
a mapping (`observable_id`) and selects the training data by the values of the
metadata (`metadata_filter`) or the observable (`observable_filter`).
`select_mappings` selects the data of a fit with the paths of the model object,
see `sbmlsim.fit.helpers.FitMappings.select` for the filters, outliers and
validation data.
"""

from __future__ import annotations

import dataclasses
from collections.abc import Collection, Hashable, Iterable, Mapping
from dataclasses import dataclass
from enum import Enum
from typing import TYPE_CHECKING, overload

from sbmlsim.fit import FitMapping, FitMappingCollection, FitSettings
from sbmlsim.fit.helpers import FitMappings, MappingFilter
from sbmlsim.fit.objects import MappingMetaData
from sbmlsim.fit.options import (
    DataAtChangeType,
    LossFunctionType,
    ResidualType,
    WeightingCurvesType,
    WeightingPointsType,
)

from pkpdsim.metadata import PKPDMappingMetaData

if TYPE_CHECKING:
    from sbmlsim.experiment import SimulationExperiment
    from sbmlsim.model.tolerances import AbsoluteTolerance

    from pkpdsim.model import PKPDModel


def default_fit_settings(
    absolute_tolerance: float | AbsoluteTolerance,
    relative_tolerance: float = 1e-6,
) -> FitSettings:
    """Get the settings of the fits of a model.

    Normalized residuals with a linear loss; the curves are weighted by the
    weights of their mappings and by their number of points, the points by their
    errors (a mapping without errors is weighted with a CV of 0.5); the
    integrator takes variable steps; a data point at the time of a change after
    the start is compared with the state before the change (`data_at_change`),
    sbmlsim's default, written out so that a change of the default does not
    change the fits.

    Args:
        absolute_tolerance: the absolute tolerance of the integrator, one value or
            an `AbsoluteTolerance` per kind and state.
        relative_tolerance: the relative tolerance of the integrator.

    Returns:
        The settings, stored with the result of a fit and read back by its report.
    """
    return FitSettings(
        residual=ResidualType.NORMALIZED,
        loss_function=LossFunctionType.LINEAR,
        weighting_curves=(
            WeightingCurvesType.MAPPING,  # user defined weights
            WeightingCurvesType.POINTS,  # number of points
        ),
        weighting_points=WeightingPointsType.ERROR_WEIGHTING,
        variable_step_size=True,
        data_at_change=DataAtChangeType.BEFORE,
        relative_tolerance=relative_tolerance,
        absolute_tolerance=absolute_tolerance,
    )


@overload
def mapping_metadata(fit_mapping: FitMapping) -> PKPDMappingMetaData: ...


@overload
def mapping_metadata[M: MappingMetaData](
    fit_mapping: FitMapping, metadata_type: type[M]
) -> M: ...


def mapping_metadata(
    fit_mapping: FitMapping, metadata_type: type[MappingMetaData] = PKPDMappingMetaData
) -> MappingMetaData:
    """Get the metadata of a fit mapping, checked against the type of the model.

    Args:
        fit_mapping: the fit mapping.
        metadata_type: the metadata class, `PKPDMappingMetaData` unless a model
            extends it.

    Returns:
        The metadata.

    Raises:
        TypeError: if the mapping has no metadata of the type.
    """
    metadata = fit_mapping.metadata
    if not isinstance(metadata, metadata_type):
        raise TypeError(f"{metadata_type.__name__} required on: '{fit_mapping}'")
    return metadata


def observable_id(fit_mapping: FitMapping) -> str:
    """Get the observable of a fit mapping, the index of its data.

    The index of the y data of the observable (`sbmlsim.data.Data.index`): the
    selection without the brackets of a concentration, `Cve_hctz` for
    `[Cve_hctz]`, `Aurine_hctz` for `Aurine_hctz`. sbmlsim 0.9 writes the
    concentration into the id of the data (`<task>__conc__Cve_hctz`), so the
    id is no longer read.

    Args:
        fit_mapping: the fit mapping.

    Returns:
        The id of the observable, e.g. `Cve_hctz`.
    """
    return fit_mapping.observable.y.index


def _show(values: frozenset[Hashable]) -> str:
    """Show values sorted, enum members by their name."""
    names = sorted(v.name if isinstance(v, Enum) else repr(v) for v in values)
    return f"[{', '.join(names)}]"


def _values(name: str, values: object) -> frozenset[Hashable]:
    """Check that the values of a clause are a collection and make a set of them.

    A single value would silently become a set of its characters when it is a
    string (every `StrEnum` member is one).

    Raises:
        TypeError: if the values are a string or not a collection.
    """
    if isinstance(values, str | bytes) or not isinstance(values, Collection):
        raise TypeError(
            f"{name} needs a collection of values, e.g. a set, not {values!r}"
        )
    return frozenset(values)


@dataclass(frozen=True)
class MetadataFilter:
    """A filter of fit mappings by the values of their metadata.

    A mapping passes when every allowed field has one of its values and no
    excluded field has one of its values. Build it with `metadata_filter`.

    Attributes:
        metadata_type: the metadata class, `PKPDMappingMetaData` unless a model
            extends it.
        allowed: the fields with the values a mapping must have.
        excluded: the fields with the values a mapping must not have.
    """

    metadata_type: type[MappingMetaData]
    allowed: tuple[tuple[str, frozenset[Hashable]], ...]
    excluded: tuple[tuple[str, frozenset[Hashable]], ...]

    def __call__(self, fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
        """Check a fit mapping, see `sbmlsim.fit.helpers.MappingFilter`."""
        metadata = mapping_metadata(fit_mapping, self.metadata_type)
        return all(
            getattr(metadata, name) in values for name, values in self.allowed
        ) and not any(
            getattr(metadata, name) in values for name, values in self.excluded
        )

    def __repr__(self) -> str:
        """Say what the filter selects."""
        clauses = [f"{name} in {_show(values)}" for name, values in self.allowed]
        clauses += [f"{name} not in {_show(values)}" for name, values in self.excluded]
        return f"MetadataFilter({self.metadata_type.__name__}: {' and '.join(clauses)})"


def metadata_filter(
    metadata_type: type[MappingMetaData] = PKPDMappingMetaData,
    /,
    *,
    excluded: Mapping[str, Collection[Hashable]] | None = None,
    **allowed: Collection[Hashable],
) -> MetadataFilter:
    """Build a filter of fit mappings by the values of their metadata.

    ```python
    # not executed
    filter_control = metadata_filter(
        route={Route.PO, Route.IV},
        coadministration={Coadministration.NONE},
        excluded={"dosing": {Dosing.MULTIPLE}},
    )
    ```

    Args:
        metadata_type: the metadata class, `PKPDMappingMetaData` unless a model
            extends it; a dataclass.
        excluded: the fields with the values a mapping must not have.
        **allowed: the fields with the values a mapping must have.

    Returns:
        The filter, called by the selection of the mappings of a fit.

    Raises:
        TypeError: if the metadata class is no dataclass (or an instance of
            one), has no field of a given name, or the values of a field are not
            a collection (a single value, also a string or `None`), i.e. when
            the module of the filters is imported.
        ValueError: for a filter without clauses, with an empty set of allowed
            values, which would select nothing, or with an empty set of excluded
            values, which would exclude nothing.
    """
    if not isinstance(metadata_type, type):
        raise TypeError(
            f"metadata_type must be the class, not an instance: {metadata_type!r}"
        )
    if not dataclasses.is_dataclass(metadata_type):
        raise TypeError(f"{metadata_type.__name__} is not a dataclass")
    excluded = excluded or {}
    names = {field.name for field in dataclasses.fields(metadata_type)}
    unknown = sorted((set(allowed) | set(excluded)) - names)
    if unknown:
        raise TypeError(f"{metadata_type.__name__} has no fields {unknown}")
    if not allowed and not excluded:
        raise ValueError("a metadata filter needs at least one clause")
    allowed_sets = {name: _values(name, values) for name, values in allowed.items()}
    excluded_sets = {name: _values(name, values) for name, values in excluded.items()}
    empty = sorted(name for name, values in allowed_sets.items() if not values)
    if empty:
        raise ValueError(f"an empty set of allowed values selects nothing: {empty}")
    empty = sorted(name for name, values in excluded_sets.items() if not values)
    if empty:
        raise ValueError(f"an empty set of excluded values excludes nothing: {empty}")
    return MetadataFilter(
        metadata_type=metadata_type,
        allowed=tuple(allowed_sets.items()),
        excluded=tuple(excluded_sets.items()),
    )


@dataclass(frozen=True)
class ObservableFilter:
    """A filter of fit mappings by their observable, see `observable_id`.

    Attributes:
        observables: the ids of the observables a mapping may have, the index
            of the data (`Data.index`).
    """

    observables: frozenset[str]

    def __call__(self, fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
        """Check a fit mapping, see `sbmlsim.fit.helpers.MappingFilter`."""
        return observable_id(fit_mapping) in self.observables

    def __repr__(self) -> str:
        """Say what the filter selects."""
        return f"ObservableFilter({', '.join(sorted(self.observables))})"


def observable_filter(observables: Collection[str]) -> ObservableFilter:
    """Build a filter of fit mappings by their observable, e.g. the PK data.

    Args:
        observables: the ids of the observables, the index of the data
            (`Data.index`), e.g. `Cve_hctz`.

    Returns:
        The filter.

    Raises:
        TypeError: if the observables are a single id or not a collection.
        ValueError: for no observables, which would select nothing.
    """
    ids = _values("observables", observables)
    if not ids:
        raise ValueError("an empty set of observables selects nothing")
    return ObservableFilter(frozenset(str(i) for i in ids))


def select_mappings(
    model: PKPDModel,
    experiment_classes: Iterable[type[SimulationExperiment]],
    filters: MappingFilter | Iterable[MappingFilter] = (),
    outliers: Iterable[str] = (),
    validation: Iterable[str] | MappingFilter = (),
) -> dict[str, list[FitMappingCollection]]:
    """Select the data of a fit of a model, see `sbmlsim.fit.helpers.FitMappings`.

    The experiments are created with the paths of the model object; the counts of
    the selected mappings per experiment and kind are printed.

    Args:
        model: the model object.
        experiment_classes: the experiments with the fit mappings.
        filters: the filters of the training data.
        outliers: the keys of the training data which is not usable.
        validation: the keys or a filter of the validation data.

    Returns:
        The fit mapping collections of all kinds by experiment id.
    """
    fit_mappings = FitMappings(
        experiment_classes=experiment_classes,
        base_path=model.path,
        data_path=list(model.data_paths),
    )
    selection = fit_mappings.select(
        filters=filters, outliers=outliers, validation=validation, print_info=False
    )
    selection.print(detail=False)
    return selection.collections
