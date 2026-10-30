"""The fit problem of the toy model per sex: one rate from the data of both sexes."""

from sbmlsim.fit import FitMappingCollection, FitParameter
from sbmlsim.fit.cli import FitDefinition

from pkpdsim.fit import default_fit_settings, select_mappings
from pkpdsim_sexes import MODEL, SEXES_PATH
from pkpdsim_sexes.experiments.decay import SexesMen, SexesWomen


def collections_decay() -> dict[str, list[FitMappingCollection]]:
    """Get the fit mappings of the decay in men and in women."""
    return select_mappings(MODEL, [SexesMen, SexesWomen], [])


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
        base_path=SEXES_PATH,
        data_path=[SEXES_PATH],
        settings=default_fit_settings(absolute_tolerance=1e-8),
    ),
}
