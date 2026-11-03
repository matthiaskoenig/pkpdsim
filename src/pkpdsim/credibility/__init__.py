"""The credibility of a model for its contexts of use.

A model declares what it is used for (`ContextOfUse`), how risky a wrong answer
is, which evidence supports it and when it is accepted (`Credibility`). The
evidence is selected from the data points of the adopted fits (`evaluate_fits`,
`mapping_table`) with the kinds their records store (`select`,
`check_declaration`). The assessment (`assess`, `CredibilityAssessment`) gives
every context of use its metrics, acceptance, domain coverage, verdict and the
factors of V&V 40 (`FACTORS`) as plain data, written as JSON. The verification
(`verify`) evaluates the fits again with tightened tolerances and a second
integration method and checks the builds, for the factors of the verification,
and runs the battery of the checks of the builds (`run_battery`): mass
balance, non-negativity, zero dose, steady state, superposition, a second
simulator, plausibility and rebuild.
"""

from pkpdsim.credibility.assessment import (
    AcceptanceRecord,
    ContextResult,
    Coverage,
    CoverageEntry,
    CredibilityAssessment,
    EvidenceResult,
    Provenance,
    QuantityRecord,
    UnclaimedRecord,
    Verdict,
    assess,
    assess_table,
    describe_selection,
)
from pkpdsim.credibility.battery import (
    Battery,
    BatteryResult,
    BatteryStatus,
    Fate,
    run_battery,
)
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
from pkpdsim.credibility.factors import FACTORS, Factor, FactorRecord, FactorStatus
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
    DEFAULT_MIN_POINTS,
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
from pkpdsim.credibility.verification import (
    BATTERY_CHECKS,
    Check,
    CheckRecord,
    Verification,
    VerificationSummary,
    verify,
)

__all__ = [
    "BATTERY_CHECKS",
    "DEFAULT_ACCEPTANCE",
    "DEFAULT_MIN_POINTS",
    "FACTORS",
    "RISK_TABLE",
    "UNCLAIMED",
    "Acceptance",
    "AcceptanceRecord",
    "AcceptanceResult",
    "AcceptanceRule",
    "AfeRange",
    "Battery",
    "BatteryResult",
    "BatteryStatus",
    "Check",
    "CheckRecord",
    "Consequence",
    "ContextOfUse",
    "ContextResult",
    "CouLevel",
    "Coverage",
    "CoverageEntry",
    "Credibility",
    "CredibilityAssessment",
    "DeclarationError",
    "Domain",
    "Evidence",
    "EvidenceResult",
    "EvidenceRole",
    "EvidenceSelection",
    "Factor",
    "FactorRecord",
    "FactorStatus",
    "Fate",
    "FitTable",
    "Influence",
    "KindMismatchError",
    "Level",
    "MappingLookup",
    "MaxGmfe",
    "Provenance",
    "QuantityOfInterest",
    "QuantityRecord",
    "Risk",
    "Selection",
    "SelectionError",
    "Selections",
    "Unclaimed",
    "UnclaimedRecord",
    "Verdict",
    "Verification",
    "VerificationSummary",
    "assess",
    "assess_table",
    "check_declaration",
    "default_acceptance",
    "describe_selection",
    "evaluate_fits",
    "kind_differences",
    "mapping_table",
    "model_risk",
    "normalize_observable",
    "recorded_kinds",
    "run_battery",
    "select",
    "verify",
]
