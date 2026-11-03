"""The credibility of a model for its contexts of use.

A model declares what it is used for (`ContextOfUse`), how risky a wrong answer
is, which evidence supports it and when it is accepted (`Credibility`). The
evidence is selected from the data points of the adopted fits (`evaluate_fits`,
`mapping_table`) with the kinds their records store (`select`,
`check_declaration`).
"""

from pkpdsim.credibility.check import check_declaration
from pkpdsim.credibility.declaration import (
    UNCLAIMED,
    Acceptance,
    AcceptanceResult,
    AcceptanceRule,
    AfeRange,
    Consequence,
    ContextOfUse,
    CouLevel,
    Credibility,
    DeclarationError,
    Domain,
    Evidence,
    EvidenceRole,
    Influence,
    Level,
    MaxGmfe,
    QuantityOfInterest,
    Risk,
    Selection,
    Unclaimed,
)
from pkpdsim.credibility.records import (
    FitTable,
    KindMismatchError,
    MappingLookup,
    evaluate_fits,
    kind_differences,
    mapping_table,
    recorded_kinds,
)
from pkpdsim.credibility.risk import (
    DEFAULT_ACCEPTANCE,
    RISK_TABLE,
    default_acceptance,
    model_risk,
)
from pkpdsim.credibility.selection import (
    EvidenceSelection,
    SelectionError,
    Selections,
    normalize_observable,
    select,
)

__all__ = [
    "DEFAULT_ACCEPTANCE",
    "RISK_TABLE",
    "UNCLAIMED",
    "Acceptance",
    "AcceptanceResult",
    "AcceptanceRule",
    "AfeRange",
    "Consequence",
    "ContextOfUse",
    "CouLevel",
    "Credibility",
    "DeclarationError",
    "Domain",
    "Evidence",
    "EvidenceRole",
    "EvidenceSelection",
    "FitTable",
    "Influence",
    "KindMismatchError",
    "Level",
    "MappingLookup",
    "MaxGmfe",
    "QuantityOfInterest",
    "Risk",
    "Selection",
    "SelectionError",
    "Selections",
    "Unclaimed",
    "check_declaration",
    "default_acceptance",
    "evaluate_fits",
    "kind_differences",
    "mapping_table",
    "model_risk",
    "normalize_observable",
    "recorded_kinds",
    "select",
]
