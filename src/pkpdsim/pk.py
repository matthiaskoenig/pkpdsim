"""The labels and units of the pharmacokinetic parameters in the figures.

The keys are the parameters of the non-compartmental analysis of pkpdutils
(`nca_single`); the base of the experiments of a model has them as `pk_labels`
and `pk_units`.
"""

from collections.abc import Mapping

#: the label of a pharmacokinetic parameter
PK_LABELS: Mapping[str, str] = {
    "auc_last": "AUCend",
    "auc_inf_obs": "AUC",
    "cl_f": "Apparent clearance",
    "cl_renal": "Renal clearance",
    "cl_fecal": "Fecal clearance",
    "cmax": "Cmax",
    "thalf": "Half-life",
    "lambda_z": "Elimination rate",
    "vz_f": "Apparent volume",
}

#: the unit of a pharmacokinetic parameter in the figures
PK_UNITS: Mapping[str, str] = {
    "auc_last": "µmole/l*hr",
    "auc_inf_obs": "µmole/l*hr",
    "cl_f": "ml/min",
    "cl_renal": "ml/min",
    "cl_fecal": "ml/min",
    "cmax": "µmole/l",
    "thalf": "hr",
    "lambda_z": "1/hr",
    "vz_f": "l",
}
