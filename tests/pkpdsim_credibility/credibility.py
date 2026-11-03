"""The contexts of use of the credibility test model."""

from pkpdsim.credibility import (
    ContextOfUse,
    CouLevel,
    Credibility,
    Domain,
    Evidence,
    EvidenceRole,
    Level,
    QuantityOfInterest,
    Unclaimed,
)
from pkpdsim.fit import metadata_filter, named_filter
from pkpdsim.metadata import Health, Route
from pkpdsim_credibility.fitting import DOSE, LIMITATION, POPULATION

_A = (QuantityOfInterest("[A]", "the concentration of A"),)

CREDIBILITY = Credibility(
    contexts=(
        ContextOfUse(
            id="describe",
            title="Describe the decay",
            question="Does the model describe the intravenous decay of A?",
            level=CouLevel.C1,
            influence=Level.LOW,
            consequence=Level.LOW,
            quantities=_A,
            domain=Domain(),
            evidence=(
                Evidence(
                    "intravenous data",
                    EvidenceRole.CALIBRATION,
                    metadata_filter(route={Route.IV}),
                ),
            ),
        ),
        ContextOfUse(
            id="dose",
            title="Decay after a dose",
            question="Does the model predict the decay after a dose of the data?",
            level=CouLevel.C2,
            influence=Level.MEDIUM,
            consequence=Level.MEDIUM,
            quantities=_A,
            domain=Domain(metadata={"route": (Route.IV,), "health": (Health.HEALTHY,)}),
            evidence=(
                Evidence("doses", EvidenceRole.VALIDATION, named_filter(DOSE)),
                Evidence(
                    "other healthy data",
                    EvidenceRole.VALIDATION,
                    Unclaimed(
                        within=metadata_filter(
                            route={Route.IV}, health={Health.HEALTHY}
                        )
                    ),
                ),
            ),
        ),
        ContextOfUse(
            id="renal",
            title="Decay in renal impairment",
            question="Does the model predict the decay in renal impairment?",
            level=CouLevel.C3,
            influence=Level.MEDIUM,
            consequence=Level.HIGH,
            quantities=_A,
            domain=Domain(
                metadata={"health": (Health.HEALTHY, Health.RENAL_IMPAIRMENT)}
            ),
            evidence=(
                Evidence(
                    "renal impairment",
                    EvidenceRole.VALIDATION,
                    named_filter(POPULATION),
                ),
            ),
        ),
    ),
    limitations=(
        Evidence(
            "scaled data",
            EvidenceRole.LIMITATION,
            named_filter(LIMITATION),
            "data at three times the decay which the model does not describe",
        ),
    ),
)

#: a declaration which validates on fitted data, CredB `fm_dose` is training data
BROKEN = Credibility(
    contexts=(
        ContextOfUse(
            id="fitted",
            title="Validation on fitted data",
            question="Does the model predict the data it was fitted on?",
            level=CouLevel.C2,
            influence=Level.LOW,
            consequence=Level.LOW,
            quantities=_A,
            domain=Domain(),
            evidence=(
                Evidence(
                    "fitted dose",
                    EvidenceRole.VALIDATION,
                    named_filter({"CredB": {"fm_dose"}}),
                ),
            ),
        ),
    ),
)
