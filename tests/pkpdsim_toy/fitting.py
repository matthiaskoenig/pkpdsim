"""The fit problem of the toy model: the rate of the decay from the data of A."""

from sbmlsim.fit import FitMappingCollection, FitParameter
from sbmlsim.fit.cli import FitDefinition

from pkpdsim.fit import default_fit_settings, metadata_filter, select_mappings
from pkpdsim.metadata import Route
from pkpdsim_toy import DATA_PATHS, MODEL, TOY_PATH
from pkpdsim_toy.experiments.decay import ToyDecay, ToyMappingMetaData

#: the intravenous data
filter_iv = metadata_filter(ToyMappingMetaData, route={Route.IV})


def collections_decay() -> dict[str, list[FitMappingCollection]]:
    """Get the fit mappings of the decay."""
    return select_mappings(MODEL, [ToyDecay], [filter_iv])


FIT_DEFINITIONS: dict[str, FitDefinition] = {
    "DECAY": FitDefinition(
        mapping_collections=collections_decay,
        parameters=[
            FitParameter(
                pid="k",
                start_value=0.5,
                lower_bound=0.01,
                upper_bound=1.0,
                unit="1/min",
            )
        ],
        base_path=TOY_PATH,
        data_path=DATA_PATHS,
        settings=default_fit_settings(absolute_tolerance=1e-8),
    ),
}
