"""The credibility factors of ASME V&V 40 and their status in the assessment.

Every factor is reported for every context of use, also the ones the
assessment does not assess, so a reader sees what is missing: `assessed` from
the evidence of the context, `platform` for a statement which holds for the
whole stack and `not_assessed` with the step of the credibility work which will
assess it.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum, unique
from typing import Final, TypedDict

__all__ = ["FACTORS", "Factor", "FactorRecord", "FactorStatus"]


@unique
class FactorStatus(StrEnum):
    """How a credibility factor is assessed."""

    ASSESSED = "assessed"
    PLATFORM = "platform"
    NOT_ASSESSED = "not_assessed"


class FactorRecord(TypedDict):
    """A factor as plain data, see `Factor`."""

    id: str
    name: str
    status: str
    statement: str
    step: int | None


@dataclass(frozen=True, kw_only=True)
class Factor:
    """A credibility factor of V&V 40.

    Attributes:
        id: the id of the factor.
        name: the name of the factor in V&V 40.
        status: how the factor is assessed.
        statement: what the factor assesses, or the platform statement.
        step: the step of the credibility work which will assess a factor which
            is not assessed, `None` otherwise.
    """

    id: str
    name: str
    status: FactorStatus
    statement: str
    step: int | None = None

    def to_dict(self) -> FactorRecord:
        """Get the factor as plain data."""
        return FactorRecord(
            id=self.id,
            name=self.name,
            status=self.status.value,
            statement=self.statement,
            step=self.step,
        )


#: the credibility factors of V&V 40 in its order
FACTORS: Final[tuple[Factor, ...]] = (
    Factor(
        id="software_quality_assurance",
        name="Software quality assurance",
        status=FactorStatus.PLATFORM,
        statement=(
            "pkpdsim and sbmlsim are released packages with tests and continuous "
            "integration; the model is released reproducibly from its repository "
            "with the versions and the commit of the provenance."
        ),
    ),
    Factor(
        id="numerical_code_verification",
        name="Numerical code verification",
        status=FactorStatus.PLATFORM,
        statement="sbmlsim runs the SBML Test Suite and the PEtab test suite.",
    ),
    Factor(
        id="discretization_error",
        name="Discretization error",
        status=FactorStatus.NOT_ASSESSED,
        statement="The convergence of the results in the tolerances of the integrator.",
        step=2,
    ),
    Factor(
        id="numerical_solver_error",
        name="Numerical solver error",
        status=FactorStatus.NOT_ASSESSED,
        statement="The results against a second simulator and the mass balance.",
        step=2,
    ),
    Factor(
        id="use_error",
        name="Use error",
        status=FactorStatus.NOT_ASSESSED,
        statement="The checks of the setup of the simulations, e.g. doses and units.",
        step=2,
    ),
    Factor(
        id="model_form",
        name="Model form",
        status=FactorStatus.NOT_ASSESSED,
        statement=(
            "The structure of the model across studies and drugs: the ratios, the "
            "cross-validation and the qualification across models."
        ),
        step=3,
    ),
    Factor(
        id="model_inputs",
        name="Model inputs",
        status=FactorStatus.NOT_ASSESSED,
        statement="The uncertainty and the sensitivity of the parameters.",
        step=4,
    ),
    Factor(
        id="test_samples",
        name="Test samples",
        status=FactorStatus.ASSESSED,
        statement="The data points of the evidence: points, mappings and studies.",
    ),
    Factor(
        id="test_conditions",
        name="Test conditions",
        status=FactorStatus.ASSESSED,
        statement="The conditions of the evidence by the metadata of its mappings.",
    ),
    Factor(
        id="equivalency_of_inputs",
        name="Equivalency of input parameters",
        status=FactorStatus.NOT_ASSESSED,
        statement=(
            "Whether the inputs of the simulations of the evidence are the ones of "
            "the context of use."
        ),
        step=3,
    ),
    Factor(
        id="output_comparison",
        name="Output comparison",
        status=FactorStatus.ASSESSED,
        statement="The acceptance of the predictions of the decisive evidence.",
    ),
    Factor(
        id="relevance_of_quantities",
        name="Relevance of the quantities of interest",
        status=FactorStatus.NOT_ASSESSED,
        statement=(
            "Whether the quantities of the evidence answer the question, e.g. the "
            "ratios of a special population."
        ),
        step=3,
    ),
    Factor(
        id="relevance_of_validation",
        name="Relevance of the validation activities to the COU",
        status=FactorStatus.ASSESSED,
        statement="The coverage of the claimed domain by the evidence.",
    ),
)
