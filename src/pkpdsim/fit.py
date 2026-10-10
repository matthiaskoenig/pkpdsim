"""The fit settings and the helpers of the mapping collections of a model.

The models of pkdb_models share their fit settings apart from the integrator
tolerances, and every `mapping_collections.py` checks the type of the metadata
of a mapping and reads the observable of a mapping the same way. The selection
of the mappings by filters, outliers and validation data is
`sbmlsim.fit.helpers.select_mapping_collections`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from sbmlsim.fit import FitMapping, FitSettings
from sbmlsim.fit.objects import MappingMetaData
from sbmlsim.fit.options import (
    LossFunctionType,
    ResidualType,
    WeightingCurvesType,
    WeightingPointsType,
)

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
    integrator takes variable steps.

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
        relative_tolerance=relative_tolerance,
        absolute_tolerance=absolute_tolerance,
    )


def mapping_metadata[M: MappingMetaData](
    fit_mapping: FitMapping, metadata_type: type[M]
) -> M:
    """Get the metadata of a fit mapping, checked against the type of the model.

    Args:
        fit_mapping: the fit mapping.
        metadata_type: the metadata class of the model.

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
    """Get the observable of a fit mapping, the id without the task.

    Args:
        fit_mapping: the fit mapping.

    Returns:
        The id of the observable, e.g. `Cve_hctz`.
    """
    return "__".join(fit_mapping.observable.y.sid.split("__")[1:])
