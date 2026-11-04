"""The leave-one-study-out cross-validation (`pkpdsim.credibility.cross_validation`).

The adopted fit of the test package fits the decay rate to the training
mappings of the studies `CredA` and `CredB`; the tests write into a copy of
its folder of fits, never into the one of the session.
"""

import dataclasses
import math
import shutil
from pathlib import Path

import pytest
from sbmlsim.fit.objects import MappingKind

from pkpdsim.cli import fit_command
from pkpdsim.credibility import CredibilityAssessment, FactorStatus, assess
from pkpdsim.credibility.cross_validation import (
    CROSS_VALIDATION_FILE,
    CrossValidation,
    HeldOutStudy,
    cross_validate,
    held_out,
    read_cross_validations,
    training_studies,
    write_cross_validation,
)
from pkpdsim.credibility.report import write_outputs
from pkpdsim.fit import fit_subset, leveled_definitions
from pkpdsim.model import PKPDModel


@pytest.fixture
def model(credibility_model: PKPDModel, tmp_path: Path) -> PKPDModel:
    """Get the credibility model with a copy of its folder of adopted fits."""
    fits = tmp_path / "fits"
    shutil.copytree(credibility_model.fits_dir, fits)
    return dataclasses.replace(credibility_model, fits_path=fits)


def _definition(model: PKPDModel):  # noqa: ANN202
    return leveled_definitions(model)[fit_subset(model.fits[0])]


def test_the_training_studies(model: PKPDModel) -> None:
    assert training_studies(_definition(model)) == ["CredA", "CredB"]


def test_held_out(model: PKPDModel) -> None:
    definition = _definition(model)
    held = held_out(definition, "CredA")
    kinds = {c.kind for c in held.mapping_collections()["CredA"]}
    assert MappingKind.TRAINING not in kinds
    assert training_studies(held) == ["CredB"]
    # the definition of the model keeps its kinds
    assert training_studies(definition) == ["CredA", "CredB"]
    with pytest.raises(KeyError, match="'Nothing' is no study"):
        held_out(definition, "Nothing").mapping_collections()


def test_cross_validate(model: PKPDModel, tmp_path: Path) -> None:
    fit_id = model.fits[0]
    result = cross_validate(model, fit_id, output_dir=tmp_path / "cv")
    assert result.fit == fit_id
    assert [s.study for s in result.studies] == ["CredA", "CredB"]
    for study in result.studies:
        assert study.n_points > 0
        assert math.isfinite(study.gmfe_held_out)
        assert study.gmfe_held_out >= 1.0
        assert math.isfinite(study.cost)
    assert result.n_points == sum(s.n_points for s in result.studies)
    assert result.statement().startswith(
        f"Leave-one-study-out cross-validation of the fit {result.subset} (2 studies, "
    )
    path = write_cross_validation(model, result)
    assert path == model.fits_dir / fit_id / CROSS_VALIDATION_FILE
    read = read_cross_validations(model)[result.subset]
    assert read.to_dict() == result.to_dict()

    assessment = assess(model, verification=False)
    assert [c.fit for c in assessment.cross_validations] == [fit_id]
    for context in assessment.contexts:
        factors = {f["id"]: f for f in context.factors}
        assert factors["model_form"]["status"] == FactorStatus.ASSESSED.value
        assert factors["model_form"]["statement"] == result.statement()
    record = tmp_path / "credibility.json"
    assessment.to_json(record)
    assert CredibilityAssessment.from_json(record).cross_validations[0].fit == fit_id
    output = write_outputs(assessment, tmp_path / "out")
    table = (output / "cross_validation.tsv").read_text().splitlines()
    assert table[0].startswith("fit\tsubset\tstudy\tn_points")
    assert len(table) == 3
    assert 'id="cross-validation"' in (output / "index.html").read_text()


def test_the_studies_must_be_training_studies(model: PKPDModel, tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="no training studies"):
        cross_validate(model, model.fits[0], ["CredC"], output_dir=tmp_path)
    with pytest.raises(ValueError, match="is not adopted"):
        cross_validate(model, "DECAY_20260101_000000__00000", output_dir=tmp_path)


def test_the_command(model: PKPDModel, tmp_path: Path) -> None:
    fit_id = model.fits[0]
    assert (
        fit_command(
            model,
            [
                "--cross-validation",
                fit_id,
                "--experiments",
                "CredB",
                "-o",
                str(tmp_path),
            ],
        )
        == {}
    )
    (result,) = read_cross_validations(model).values()
    assert [s.study for s in result.studies] == ["CredB"]
    with pytest.raises(SystemExit) as exc:
        fit_command(model, ["--cross-validation", fit_id, "--experiments", "CredC"])
    assert exc.value.code == 1


def test_a_round_trip_without_points() -> None:
    result = CrossValidation(
        fit="PK_20260101_000000__00000",
        subset="PK",
        studies=(
            HeldOutStudy(
                study="A",
                n_points=0,
                gmfe_held_out=math.nan,
                gmfe_fitted=math.nan,
                within_2fold_held_out=math.nan,
                within_2fold_fitted=math.nan,
                cost=1.0,
            ),
        ),
        gmfe_held_out=math.nan,
        gmfe_fitted=math.nan,
        seed=1,
        date="2026-11-04T12:00:00+01:00",
    )
    read = CrossValidation.from_dict(result.to_dict())
    assert read.to_dict() == result.to_dict()
    assert read.to_dict()["gmfe_held_out"] is None
