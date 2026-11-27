"""The command `run_petab`: the fit problems of a model as PEtab v2 problems."""

import json
from collections.abc import Callable
from pathlib import Path

import pytest

from pkpdsim.cli import petab_command, run_petab_main
from pkpdsim.model import PKPDModel
from pkpdsim_toy import MODEL

type AsModel = Callable[[PKPDModel], PKPDModel]

#: the records of a run of a kind of the run layer
RUN_FILES = {"run.json", "spec.json", "summary.json", "progress.jsonl"}


def _adopted(fits_dir: Path, as_model: AsModel) -> tuple[PKPDModel, str]:
    """Get the toy model with its adopted fit and the id of the fit."""
    (record,) = fits_dir.iterdir()
    model = as_model(MODEL.replace(adopted_ids=[record.name], fits_dir=fits_dir))
    return model, record.name


def test_run_petab_writes_the_problems_of_the_adopted_fit(
    fits_dir: Path, tmp_path: Path, as_model: AsModel
) -> None:
    model, fit_id = _adopted(fits_dir, as_model)
    out = petab_command(model, ["-o", str(tmp_path / "petab")])
    assert out == tmp_path / "petab"
    directory = out / "DECAY"
    # the training problem next to the records of the run, the toy fit has no
    # validation data
    assert {path.name for path in directory.iterdir()} >= RUN_FILES
    assert (directory / "training" / "problem.yaml").is_file()
    assert not (directory / "validation").exists()
    readme = (directory / "training" / "README.md").read_text(encoding="utf-8")
    assert f"| parameter values | the best parameters of the fit `{fit_id}` |" in readme
    spec = json.loads((directory / "spec.json").read_text(encoding="utf-8"))
    assert spec["fit"] == str((fits_dir / fit_id).resolve())
    assert spec["split"] and spec["unsupported"] == "skip"
    summary = json.loads((directory / "summary.json").read_text(encoding="utf-8"))
    assert summary["problems"] == {"training": "petab/training/problem.yaml"}

    # a second run replaces the first
    petab_command(model, ["-o", str(tmp_path / "petab")])
    assert (directory / "training" / "problem.yaml").is_file()


def test_run_petab_without_an_adopted_fit_writes_the_start_values(
    tmp_path: Path, as_model: AsModel
) -> None:
    model = as_model(MODEL.replace(adopted_ids=[], results_path=tmp_path / "results"))
    out = petab_command(model, ["--portable", "--no-split"])
    assert out == model.results_dir / "petab"
    directory = out / "DECAY"
    readme = (directory / "README.md").read_text(encoding="utf-8")
    assert "the start values of the fit problem `DECAY`, no fit" in readme
    assert "The extension `sbmlsim` is not required" in readme
    spec = json.loads((directory / "spec.json").read_text(encoding="utf-8"))
    assert "fit" not in spec or spec["fit"] is None
    assert spec["required_extension"] is False and spec["split"] is False


def test_run_petab_refuses_an_unknown_subset(
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit) as exit_info:
        run_petab_main(["pkpdsim_toy", "--subset", "PK"])
    assert exit_info.value.code == 2
    assert "invalid choice: 'PK'" in capsys.readouterr().err


def test_run_petab_refuses_an_output_in_the_adopted_fits(
    fits_dir: Path, as_model: AsModel, capsys: pytest.CaptureFixture[str]
) -> None:
    model, _ = _adopted(fits_dir, as_model)
    with pytest.raises(SystemExit) as exit_info:
        petab_command(model, ["-o", str(fits_dir)])
    assert exit_info.value.code == 2
    assert "is or contains the input" in capsys.readouterr().err
