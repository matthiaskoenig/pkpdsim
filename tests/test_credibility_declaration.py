"""Tests of the declaration of the contexts of use."""

import math

import numpy as np
import pandas as pd
import pytest
from sbmlsim.fit.evaluation import Fold

from pkpdsim.credibility import (
    UNCLAIMED,
    Acceptance,
    ContextOfUse,
    CouLevel,
    Credibility,
    DeclarationError,
    Domain,
    Evidence,
    EvidenceRole,
    Level,
    MaxGmfe,
    QuantityOfInterest,
    default_acceptance,
    model_risk,
)
from pkpdsim.metadata import Health, Route

L, M, H = Level.LOW, Level.MEDIUM, Level.HIGH


def _select(key: str, fit_mapping: object) -> bool:
    return True


def _cou(**kw: object) -> ContextOfUse:
    args: dict[str, object] = {
        "id": "dose",
        "title": "Dose",
        "question": "What follows from a dose?",
        "level": CouLevel.C2,
        "influence": M,
        "consequence": M,
        "quantities": (QuantityOfInterest("Cve", "plasma"),),
        "domain": Domain(),
        "evidence": (Evidence("v", EvidenceRole.VALIDATION, _select),),
    }
    args.update(kw)
    return ContextOfUse(**args)  # ty: ignore[invalid-argument-type]


@pytest.mark.parametrize(
    ("i", "c", "expected"),
    [
        (L, L, L), (L, M, L), (L, H, M),
        (M, L, L), (M, M, M), (M, H, H),
        (H, L, M), (H, M, H), (H, H, H),
    ],
)  # fmt: skip
def test_the_risk_table(i: Level, c: Level, expected: Level) -> None:
    assert model_risk(i, c) is expected


def test_the_default_acceptance() -> None:
    expected = {
        L: (["within_2fold"], [">= 70"]),
        M: (["within_2fold", "GMFE"], [">= 80", "<= 1.5"]),
        H: (["within_2fold", "within_1.5fold", "GMFE"], [">= 90", ">= 70", "<= 1.3"]),
    }
    for risk, (names, limits) in expected.items():
        rules = default_acceptance(risk)
        assert [r.name for r in rules] == names
        assert [r.limit for r in rules] == limits


def test_a_rule_evaluates_points() -> None:
    points = pd.DataFrame(
        {"DV": [1, 1, 1, 1], "IPRED": [1.1, 1.9, 2.5, 0.4], "DV_SD": np.nan}
    )
    result = Acceptance(Fold(2.0), 50.0).evaluate(points)
    assert result.value == 50.0
    assert result.met is True
    assert MaxGmfe(1.5).evaluate(points).met is False


def test_a_rule_without_applicable_points_is_not_evaluable() -> None:
    points = pd.DataFrame({"DV": [0, -1], "IPRED": [1, 1], "DV_SD": np.nan})
    result = Acceptance(Fold(2.0), 70.0).evaluate(points)
    assert math.isnan(result.value)
    assert result.met is None


def test_a_cou_takes_the_default_of_its_risk() -> None:
    cou = _cou(influence=M, consequence=L)
    assert cou.risk is L
    assert cou.rules == default_acceptance(L)


def test_the_declaration_reports_every_problem() -> None:
    good = _cou()
    bad = _cou(
        id="Bad-Id",
        question="Not a question.",
        quantities=(),
        evidence=(
            Evidence("u1", EvidenceRole.VALIDATION, UNCLAIMED),
            Evidence("u2", EvidenceRole.VALIDATION, UNCLAIMED),
            Evidence("u3", EvidenceRole.CALIBRATION, UNCLAIMED),
        ),
        domain=Domain(metadata={"dose": (Route.PO,), "health": (Route.PO,)}),
    )
    no_validation = _cou(
        id="nov", evidence=(Evidence("c", EvidenceRole.CALIBRATION, _select),)
    )
    no_calibration = _cou(
        id="noc",
        level=CouLevel.C1,
        evidence=(Evidence("v", EvidenceRole.VALIDATION, _select),),
    )
    with pytest.raises(DeclarationError) as err:
        Credibility(
            contexts=(good, good, bad, no_validation, no_calibration),
            limitations=(Evidence("l", EvidenceRole.VALIDATION, _select),),
        )
    text = "\n".join(err.value.problems)
    for needle in [
        "'dose': the id is not unique",
        "'Bad-Id': the id must match",
        "'Bad-Id': the question must be a question",
        "'Bad-Id': at least one quantity",
        "'Bad-Id': at most one unclaimed",
        "'Bad-Id': the unclaimed evidence 'u3' must have the role validation",
        "'Bad-Id': the domain field 'dose'",
        "'Bad-Id': the domain value",
        "'nov': a C2 context needs validation evidence",
        "'noc': a C1 context needs calibration evidence",
        "model limitation 'l' must have the role limitation",
    ]:
        assert needle in text, needle
    assert len(err.value.problems) == 11


def test_a_c4_cou_needs_only_its_extrapolation() -> None:
    cou = _cou(
        level=CouLevel.C4,
        evidence=(Evidence("c", EvidenceRole.CALIBRATION, _select),),
        extrapolation="platform qualification",
        domain=Domain(metadata={"route": [Route.PO], "health": [Health.HEALTHY]}),
    )
    assert Credibility(contexts=(cou,)).context("dose") is cou
