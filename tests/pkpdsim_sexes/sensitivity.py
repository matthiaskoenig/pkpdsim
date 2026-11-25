"""The sensitivity analysis of the toy model per sex, `sensitivities` of its model object.

The PK observables of `A` in men and in women, a local analysis and a small
Latin hypercube, which `run_sensitivity pkpdsim_sexes` and the pipeline of the
package run.
"""

from sbmlsim.sensitivity.parameters import (
    SensitivityParameter,
    parameters_from_sbml,
    set_relative_bounds,
)
from sbmlsim.sensitivity.runner import (
    SamplingAnalysis,
    SensitivityGroup,
    SensitivityOutput,
    SensitivityRunner,
)
from sbmlsim.simulation import Simulation

from pkpdsim.sensitivity import NonzeroLocalAnalysis, pk_observables
from pkpdsim.sex import Sex
from pkpdsim_sexes import MODEL_PATH


def sensitivity_parameters() -> list[SensitivityParameter]:
    """Get the parameters of the toy model, read from its male build, +-15 %."""
    parameters = parameters_from_sbml(MODEL_PATH)
    set_relative_bounds(parameters, 0.15)
    return parameters


SIMULATION = Simulation(end=60, steps=600, time_unit="min")

GROUPS = [
    SensitivityGroup("men", "Men", SIMULATION),
    SensitivityGroup("women", "Women", SIMULATION, variant=Sex.FEMALE),
]

#: the analyses of the model, the local one by default
SENSITIVITY = SensitivityRunner(
    parameters=sensitivity_parameters,
    outputs=[
        SensitivityOutput("A_auc_inf_obs", "A AUC∞"),
        SensitivityOutput("A_thalf", "A half-life", observable="pk_A.thalf"),
    ],
    observables=lambda: pk_observables(
        "A", "[A]", {"auc_inf_obs": "mM*min", "thalf": "min"}
    ),
    analyses={
        "local": NonzeroLocalAnalysis(groups=GROUPS),
        "sampling": SamplingAnalysis(groups=GROUPS, n=8),
    },
    default=("local",),
)
