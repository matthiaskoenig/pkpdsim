"""The fit of the credibility test model: the rate of the decay, with its named sets."""

from sbmlsim.fit import FitMappingCollection, FitParameter
from sbmlsim.fit.cli import FitDefinition

from pkpdsim.fit import (
    NamedMappings,
    default_fit_settings,
    named_filter,
    select_mappings,
)
from pkpdsim_credibility import MODEL
from pkpdsim_credibility.experiments import CredA, CredB
from pkpdsim_toy import DATA_PATHS, TOY_PATH

#: the validation of the dose, by study
DOSE: NamedMappings = {"CredA": {"fm_dose"}}
#: the validation of the population
POPULATION: NamedMappings = {"CredA": {"fm_renal"}}
#: the data the model does not describe
LIMITATION: NamedMappings = {"CredA": {"fm_limit"}}
#: data of zeros, whose fold errors are not evaluable
ZERO: NamedMappings = {"CredB": {"fm_zero"}}
#: validation data no evidence claims
SUBJECTS: NamedMappings = {"CredA": {"fm_subject"}}

#: every named set of validation data
VALIDATION: NamedMappings = {}
for _set in (DOSE, POPULATION, LIMITATION, ZERO, SUBJECTS):
    for _study, _keys in _set.items():
        VALIDATION.setdefault(_study, set()).update(_keys)


def collections_decay() -> dict[str, list[FitMappingCollection]]:
    """Get the fit mappings of the decay, with the validation and the outlier."""
    return select_mappings(
        MODEL,
        [CredA, CredB],
        filters=(),
        outliers=["fm_outlier"],
        validation=named_filter(VALIDATION),
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
