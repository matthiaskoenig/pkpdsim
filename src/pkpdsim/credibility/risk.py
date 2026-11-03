"""The model risk and the default acceptance of a context of use.

The model risk follows ICH M15 and V&V 40 as a lookup of the influence of the
model on the decision and the consequence of a wrong decision. It is a table
and not a formula, so it can be changed in one place.
"""

from collections.abc import Mapping
from types import MappingProxyType
from typing import Final

from sbmlsim.fit.evaluation import Fold

from pkpdsim.credibility.declaration import Acceptance, AcceptanceRule, Level, MaxGmfe

__all__ = ["DEFAULT_ACCEPTANCE", "RISK_TABLE", "default_acceptance", "model_risk"]

#: the model risk by influence and consequence
RISK_TABLE: Final[Mapping[tuple[Level, Level], Level]] = MappingProxyType(
    {
        (Level.LOW, Level.LOW): Level.LOW,
        (Level.LOW, Level.MEDIUM): Level.LOW,
        (Level.LOW, Level.HIGH): Level.MEDIUM,
        (Level.MEDIUM, Level.LOW): Level.LOW,
        (Level.MEDIUM, Level.MEDIUM): Level.MEDIUM,
        (Level.MEDIUM, Level.HIGH): Level.HIGH,
        (Level.HIGH, Level.LOW): Level.MEDIUM,
        (Level.HIGH, Level.MEDIUM): Level.HIGH,
        (Level.HIGH, Level.HIGH): Level.HIGH,
    }
)

#: the acceptance of a context of use which does not declare its own, by risk
DEFAULT_ACCEPTANCE: Final[Mapping[Level, tuple[AcceptanceRule, ...]]] = (
    MappingProxyType(
        {
            Level.LOW: (Acceptance(Fold(2.0), 70.0),),
            Level.MEDIUM: (Acceptance(Fold(2.0), 80.0), MaxGmfe(1.5)),
            Level.HIGH: (
                Acceptance(Fold(2.0), 90.0),
                Acceptance(Fold(1.5), 70.0),
                MaxGmfe(1.3),
            ),
        }
    )
)


def model_risk(influence: Level, consequence: Level) -> Level:
    """Get the model risk of an influence and a consequence.

    Args:
        influence: the influence of the model on the decision.
        consequence: the consequence of a wrong decision.

    Returns:
        The model risk.
    """
    return RISK_TABLE[(influence, consequence)]


def default_acceptance(risk: Level) -> tuple[AcceptanceRule, ...]:
    """Get the acceptance a context of use of a risk takes without its own.

    Args:
        risk: the model risk.

    Returns:
        The rules, all of which must be met.
    """
    return DEFAULT_ACCEPTANCE[risk]
