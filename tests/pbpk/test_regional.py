"""Tests of the derivation of a regional absorption profile with the model (R36b)."""

from pathlib import Path
from typing import TYPE_CHECKING

import pytest

from pkpdsim.pbpk import (
    Organ,
    RegionalAbsorption,
    Route,
    SegmentedIntestine,
    Variants,
    build,
    derive_regional_absorption,
    regional_key,
    release_auc_ratios,
)

from .conftest import definition, load

if TYPE_CHECKING:
    import roadrunner

S = "drg"
DOSE = 10.0
END, POINTS = 2 * 24 * 60.0, 2881
"""Two days at one point per minute: the AUC of the simple drug, fast enough for a test."""


def _model(
    tmp_path: Path, routes: list[Route] | None = None
) -> "roadrunner.RoadRunner":
    """The simple drug with an oral solution and a profile of the lower ileum and the colon."""
    d = definition(
        routes=routes or [Route.oral_solution(S)],
        variants=Variants(
            {
                Organ.SMALL_INTESTINE: SegmentedIntestine(
                    regional=RegionalAbsorption({"ile2": 0.5, "colon": 0.3})
                )
            }
        ),
    )
    return load(build(d, tmp_path))


def test_the_key_of_a_release_segment() -> None:
    """A segment of the small intestine has its own factor, the large intestine the group."""
    assert regional_key("ile2") == "ile2"
    assert regional_key("du") == "du"
    assert regional_key("asc") == "colon"
    assert regional_key("rec") == "colon"
    with pytest.raises(ValueError, match="'st' is no absorbing segment"):
        regional_key("st")


def test_the_derived_factors_give_the_ratios(tmp_path: Path) -> None:
    """The ratios of known factors, derived back from the factors of the model.

    The release in the lower ileum and in the colon ascendens give AUC ratios
    below the oral solution, the colon the lower one; from other factors of
    the model the derivation finds the known ones.
    """
    r = _model(tmp_path)
    known = {"ile2": 0.3, "colon": 0.1}
    ratios = release_auc_ratios(
        r, S, ["ile2", "asc"], DOSE, factors=known, end=END, points=POINTS
    )
    assert 0 < ratios["asc"] < ratios["ile2"] < 1
    model = release_auc_ratios(r, S, ["ile2", "asc"], DOSE, end=END, points=POINTS)
    assert model["ile2"] > ratios["ile2"]
    derived = derive_regional_absorption(
        r, S, ratios, DOSE, xtol=1e-8, end=END, points=POINTS
    )
    assert derived == pytest.approx(known, rel=1e-5)
    again = release_auc_ratios(
        r, S, ["ile2", "asc"], DOSE, factors=derived, end=END, points=POINTS
    )
    assert again == pytest.approx(ratios, rel=1e-6)
    # the model keeps its own factors, not the derived ones
    assert r["f_ka_ile2_drg"] == 0.5
    assert r["f_ka_colon_drg"] == 0.3


def test_what_the_derivation_cannot_do_is_refused(tmp_path: Path) -> None:
    """Two releases of one factor, an unreachable ratio, a model without an oral solution."""
    r = _model(tmp_path)
    with pytest.raises(ValueError, match=r"share the factors \['colon'\]"):
        derive_regional_absorption(
            r, S, {"asc": 0.2, "des": 0.1}, DOSE, end=END, points=POINTS
        )
    with pytest.raises(ValueError, match=r"No factor of 'colon' within"):
        derive_regional_absorption(
            r, S, {"asc": 0.999}, DOSE, bounds=(1e-4, 1.0), end=END, points=POINTS
        )
    tablet = _model(tmp_path / "tablet", routes=[Route.oral(S)])
    with pytest.raises(ValueError, match=r"has no \['SOLDOSE_drg'\]"):
        release_auc_ratios(tablet, S, ["ile2"], DOSE, end=END, points=POINTS)


def test_the_condition_of_the_study_is_taken_back(tmp_path: Path) -> None:
    """The ratios in a condition (`changes`, initial values); afterwards `r` is reset.

    Fed and at 100 kg (the initial assignments follow: `T_st_eff`, the
    volumes) the ratio differs from the fasted reference; afterwards the
    initial values are those of before and every value is reset, so the
    fasted ratio comes again.
    """
    r = _model(tmp_path)
    fasted = release_auc_ratios(r, S, ["ile2"], DOSE, end=END, points=POINTS)
    fed = release_auc_ratios(
        r,
        S,
        ["ile2"],
        DOSE,
        changes={"FED": 1.0, "BW": 100.0},
        end=END,
        points=POINTS,
    )
    assert fed["ile2"] != pytest.approx(fasted["ile2"], rel=1e-6)
    assert r["FED"] == 0.0
    assert r["BW"] == pytest.approx(73.0)
    assert r["T_st_eff"] == pytest.approx(22.79)
    assert r["SOLDOSE_drg"] == 0.0
    assert r["f_ka_ile2_drg"] == 0.5
    again = release_auc_ratios(r, S, ["ile2"], DOSE, end=END, points=POINTS)
    assert again == pytest.approx(fasted, rel=1e-12)
    with pytest.raises(ValueError, match=r"has no \['X_nowhere'\]"):
        release_auc_ratios(r, S, ["ile2"], DOSE, changes={"X_nowhere": 1.0})
