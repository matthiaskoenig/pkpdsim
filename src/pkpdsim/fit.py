"""The fit settings, the metadata check and the sampling sites of the models.

The models of pkdb_models share their fit settings apart from the integrator
tolerances (`default_fit_settings`). Every `mapping_collections.py` checks the
type of the metadata of a mapping (`mapping_metadata`) and selects the training
data by the values of the metadata (`metadata_filter`, the filter of
`sbmlsim.fit.filters` for `PKPDMappingMetaData`). The other filters, the
selection of the mappings, the start of a fit, the fixed adopted fits and the
records of the adopted fits are those of `sbmlsim.fit` (`sbmlsim.fit.filters`,
`sbmlsim.fit.helpers.select_mappings`, `sbmlsim.fit.project_fits`,
`sbmlsim.fit.records`).
"""

from __future__ import annotations

from collections.abc import Collection, Hashable, Iterable, Mapping
from typing import TYPE_CHECKING, overload

from sbmlsim.fit import FitMapping, FitSettings, filters
from sbmlsim.fit.filters import MetadataFilter, observable_id
from sbmlsim.fit.helpers import FitMappings
from sbmlsim.fit.objects import MappingMetaData
from sbmlsim.fit.options import (
    DataAtChangeType,
    LossFunctionType,
    ResidualType,
    WeightingCurvesType,
    WeightingPointsType,
)

from pkpdsim.metadata import PKPDMappingMetaData, SamplingSite

if TYPE_CHECKING:
    from sbmlsim.model.tolerances import AbsoluteTolerance


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


def metadata_filter(
    metadata_type: type[PKPDMappingMetaData] = PKPDMappingMetaData,
    /,
    *,
    excluded: Mapping[str, Collection[Hashable]] | None = None,
    **allowed: Collection[Hashable],
) -> MetadataFilter:
    """Build a filter of fit mappings by the values of their PK metadata.

    The filter of `sbmlsim.fit.filters.metadata_filter` with the metadata of
    the models, `PKPDMappingMetaData`; the sampling site it compares is the
    resolved one (`PKPDMappingMetaData.filter_value`).

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
            extends it; positional only.
        excluded: the fields with the values a mapping must not have.
        **allowed: the fields with the values a mapping must have.

    Returns:
        The filter, called by the selection of the mappings of a fit.

    Raises:
        TypeError: if the metadata class is no dataclass (or an instance of
            one), has no field of a given name, or the values of a field are not
            a collection (a single value, also a string or `None`).
        ValueError: for a filter without clauses, with an empty set of allowed
            values, which would select nothing, or with an empty set of excluded
            values, which would exclude nothing.
    """
    return filters.metadata_filter(metadata_type, excluded=excluded, **allowed)


def check_sampling_sites(
    fit_mappings: FitMappings | Iterable[FitMapping],
    sites: Mapping[str, SamplingSite],
    metadata_type: type[PKPDMappingMetaData] = PKPDMappingMetaData,
) -> list[str]:
    """Check the sampling site of the mappings against their observable.

    Args:
        fit_mappings: the fit mappings, or the instantiated experiments of a
            model (`FitMappings`, every fit mapping of every experiment).
        sites: the sampling site of every observable by observable id, e.g.
            `Carm_hctz` is `ARM_VEIN` and `Aurine_hctz` is `NOT_APPLICABLE`.
        metadata_type: the metadata class, `PKPDMappingMetaData` unless a model
            extends it.

    Returns:
        One problem per mapping whose site differs from the site of its
        observable or whose observable is not in `sites`; empty if all fit.
    """
    if isinstance(fit_mappings, FitMappings):
        mappings: list[tuple[str, FitMapping]] = [
            (f"{experiment_id}.{key}", fit_mapping)
            for experiment_id, experiment in fit_mappings.experiments.items()
            for key, fit_mapping in experiment._fit_mappings.items()
        ]
    else:
        mappings = [(str(fm), fm) for fm in fit_mappings]
    problems: list[str] = []
    for name, fit_mapping in mappings:
        observable = observable_id(fit_mapping)
        site = mapping_metadata(fit_mapping, metadata_type).resolved_sampling_site()
        if observable not in sites:
            problems.append(f"{name}: observable '{observable}' has no sampling site")
        elif site is not sites[observable]:
            problems.append(
                f"{name}: sampling site {site.name} of '{observable}' "
                f"differs from {sites[observable].name}"
            )
    return problems
