"""The evidence and the battery of the credibility of a PBPK model.

The credibility of a model for its contexts of use is `sbmlsim.credibility`:
the declaration, the selection of the evidence, the assessment, the
verification, the report and the scorecard. pkpdsim adds what is specific to
the models of pkdb_models, registered by every `pkpdsim.model.PKPDModel`
(`pkpdsim.model.CREDIBILITY_SETTINGS`: the metadata `PKPDMappingMetaData`, the
fields of the domain `DOMAIN_FIELDS`, the extensions and the battery hook):

- the evidence extensions (`EXTENSIONS`, `sbmlsim.credibility.extensions`):
  the ratio pairs of the ratio quantities with the limits of Guest et al.
  (`ratios`, `Ratio`, `GuestAcceptance`), the temporal validation on the
  studies collected after their fit (`temporal`), the equivalency of the
  simulated and the reported inputs of the subjects (`inputs`) and the
  provenance of the study data (`studies`);
- the battery hook of the PK checks (`PK_BATTERY`, `pkpdsim.credibility.battery`):
  the dose states of the model, mass balance, zero dose and superposition on
  the doses of every route of every substance, with the fates of the doses.
"""

from typing import Final

from pkpdsim.credibility.battery import (
    DOSE_STATES,
    MASS_BALANCE,
    PK_BATTERY,
    PK_CHECKS,
    SUPERPOSITION,
    ZERO_DOSE,
    Fate,
    FlatBuild,
    PKBattery,
)
from pkpdsim.credibility.inputs import (
    DEFAULT_COVARIATES,
    INPUTS,
    InputCheck,
    InputsExtension,
    InputsRecord,
    InputStatus,
    InputSummary,
    input_checks,
)
from pkpdsim.credibility.ratios import (
    DEFAULT_RATIO_ACCEPTANCE,
    GUEST_DELTA,
    RATIOS,
    GuestAcceptance,
    Ratio,
    RatioRecord,
    RatiosExtension,
    Summary,
    guest_limit,
)
from pkpdsim.credibility.studies import STUDIES, StudiesExtension, StudyData, study_data
from pkpdsim.credibility.temporal import (
    TEMPORAL,
    Chronology,
    TemporalExtension,
    TemporalRecord,
    TemporalStudy,
    chronology,
)

#: the evidence extensions of pkpdsim, in the order a `PKPDModel` registers them
EXTENSIONS: Final = (RATIOS, TEMPORAL, INPUTS, STUDIES)

__all__ = [
    "DEFAULT_COVARIATES",
    "DEFAULT_RATIO_ACCEPTANCE",
    "DOSE_STATES",
    "EXTENSIONS",
    "GUEST_DELTA",
    "INPUTS",
    "MASS_BALANCE",
    "PK_BATTERY",
    "PK_CHECKS",
    "RATIOS",
    "STUDIES",
    "SUPERPOSITION",
    "TEMPORAL",
    "ZERO_DOSE",
    "Chronology",
    "Fate",
    "FlatBuild",
    "GuestAcceptance",
    "InputCheck",
    "InputStatus",
    "InputSummary",
    "InputsExtension",
    "InputsRecord",
    "PKBattery",
    "Ratio",
    "RatioRecord",
    "RatiosExtension",
    "StudiesExtension",
    "StudyData",
    "Summary",
    "TemporalExtension",
    "TemporalRecord",
    "TemporalStudy",
    "chronology",
    "guest_limit",
    "input_checks",
    "study_data",
]
