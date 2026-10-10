"""The fit problem of the toy model: the rate of the decay from the data of A."""

from sbmlsim.fit import FitMappingCollection, FitParameter
from sbmlsim.fit.cli import FitDefinition
from sbmlsim.fit.helpers import select_mapping_collections

from pkpdsim.fit import default_fit_settings
from pkpdsim_toy import DATA_PATHS, TOY_PATH
from pkpdsim_toy.experiments.decay import ToyDecay


def collections_decay() -> dict[str, list[FitMappingCollection]]:
    """Get the fit mappings of the decay."""
    return select_mapping_collections(
        experiment_classes=[ToyDecay],
        base_path=TOY_PATH,
        data_path=DATA_PATHS,
        print_info=False,
    )


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
