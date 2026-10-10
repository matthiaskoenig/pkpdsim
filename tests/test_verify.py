"""The verification of a change of a model, run on the toy model."""

import json
import pickle
import re
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest
from PIL import Image

from pkpdsim.verify import _attribute, describe, verify
from pkpdsim_toy import DATA_PATHS, MODEL, MODEL_PATH, TOY_PATH

TESTS_DIR = Path(__file__).parent


def test_describe_a_model_object() -> None:
    info = describe("pkpdsim_toy")
    assert info.base_path == TOY_PATH
    assert info.model_path == MODEL_PATH
    assert info.data_path == DATA_PATHS
    assert info.experiments is not None
    assert set(info.experiments) == {"studies", "misc", "all"}
    assert info.absolute_tolerance == MODEL.simulator.absolute_tolerance
    assert info.fit_definitions is not None
    assert set(info.fit_definitions) == {"DECAY"}


def test_describe_a_package_by_the_conventions() -> None:
    info = describe("pkpdsim_legacy")
    assert info.base_path == TOY_PATH
    assert info.model_path == MODEL_PATH
    assert info.experiments is not None
    assert info.absolute_tolerance == 1e-10
    assert info.fit_definitions is not None


def test_dump_and_compare(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    common = ["--root", str(TESTS_DIR), "--atol", "model", "--rtol", "1e-10"]
    old, new = tmp_path / "old.pkl", tmp_path / "new.pkl"
    assert (
        verify(["dump", "--model", "pkpdsim_legacy", "--out", str(old), *common]) == 0
    )
    assert (
        verify(
            [
                "dump",
                "--model",
                "pkpdsim_toy",
                "--out",
                str(new),
                "--order",
                str(old),
                *common,
            ]
        )
        == 0
    )
    dumped = pickle.loads(new.read_bytes())
    assert dumped["results"]
    assert all(isinstance(v[1], np.ndarray) for v in dumped["results"].values())
    capsys.readouterr()
    assert verify(["compare", str(old), str(new)]) == 0
    match = re.search(r"simulations: (\d+) arrays compared", capsys.readouterr().out)
    assert match
    assert int(match.group(1)) == len(dumped["results"]) > 0


def test_definitions(tmp_path: Path) -> None:
    out = tmp_path / "definitions.json"
    assert (
        verify(
            [
                "definitions",
                "--root",
                str(TESTS_DIR),
                "--model",
                "pkpdsim_toy",
                "--out",
                str(out),
            ]
        )
        == 0
    )
    content = json.loads(out.read_text(encoding="utf-8"))
    assert set(content) == {"DECAY"}
    assert [p["pid"] for p in content["DECAY"]["parameters"]] == ["k"]
    legacy = tmp_path / "legacy.json"
    verify(
        [
            "definitions",
            "--root",
            str(TESTS_DIR),
            "--model",
            "pkpdsim_legacy",
            "--out",
            str(legacy),
        ]
    )
    assert legacy.read_text(encoding="utf-8") == out.read_text(encoding="utf-8")


def test_sbml(tmp_path: Path) -> None:
    copy = tmp_path / "copy.xml"
    shutil.copy(MODEL_PATH, copy)
    assert verify(["sbml", str(MODEL_PATH), str(copy)]) == 0
    changed = tmp_path / "changed.xml"
    changed.write_text(
        MODEL_PATH.read_text(encoding="utf-8").replace('value="0.2"', 'value="0.3"'),
        encoding="utf-8",
    )
    assert verify(["sbml", str(MODEL_PATH), str(changed)]) == 1


def write_png(path: Path, value: int = 0) -> None:
    """Write a small image, a different one for a different value."""
    path.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(np.full((4, 4, 3), value, dtype=np.uint8)).save(path)


def test_png(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    output = MODEL.run_experiments(
        MODEL.load_experiments()["studies"], tmp_path / "a", figure_formats=["png"]
    )
    figures = list(output.rglob("*.png"))
    assert figures
    # the figures are below the experiments, not in the directory itself
    assert not list(output.glob("*.png"))
    capsys.readouterr()
    assert verify(["png", str(output), str(output)]) == 0
    assert f"{len(figures)} identical, 0 differ" in capsys.readouterr().out


def test_png_compares_recursively_by_relative_path(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    old, new = tmp_path / "old", tmp_path / "new"
    for directory in (old, new):
        write_png(directory / "A" / "fig.png")
        write_png(directory / "_figures" / "fig.png", 7)
    capsys.readouterr()
    assert verify(["png", str(old), str(new)]) == 0
    assert "2 identical, 0 differ" in capsys.readouterr().out
    write_png(new / "A" / "fig.png", 200)
    assert verify(["png", str(old), str(new)]) == 1
    assert "differs: A" in capsys.readouterr().out


def test_png_reports_files_on_one_side(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    old, new = tmp_path / "old", tmp_path / "new"
    write_png(old / "A" / "fig.png")
    write_png(new / "A" / "fig.png")
    write_png(old / "B" / "only_old.png")
    write_png(new / "C" / "only_new.png")
    capsys.readouterr()
    assert verify(["png", str(old), str(new)]) == 1
    out = capsys.readouterr().out
    assert "only in old: B" in out
    assert "only in new: C" in out
    assert "1 identical, 2 differ" in out


def test_png_without_figures_or_directory_fails(tmp_path: Path) -> None:
    empty = tmp_path / "empty"
    empty.mkdir()
    assert verify(["png", str(empty), str(empty)]) == 2
    assert verify(["png", str(tmp_path / "missing"), str(empty)]) == 2
    assert verify(["png", str(empty), str(tmp_path / "missing")]) == 2


def test_sbml_compares_directories_recursively(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    old, new = tmp_path / "old", tmp_path / "new"
    for directory in (old, new):
        (directory / "sub").mkdir(parents=True)
        shutil.copy(MODEL_PATH, directory / "sub" / "toy.xml")
    capsys.readouterr()
    assert verify(["sbml", str(old), str(new)]) == 0
    assert "sub" in capsys.readouterr().out
    (new / "sub" / "toy.xml").write_text(
        MODEL_PATH.read_text(encoding="utf-8").replace('value="0.2"', 'value="0.3"'),
        encoding="utf-8",
    )
    assert verify(["sbml", str(old), str(new)]) == 1
    shutil.copy(MODEL_PATH, new / "extra.xml")
    capsys.readouterr()
    assert verify(["sbml", str(old), str(new)]) == 1
    assert "extra.xml: only in new" in capsys.readouterr().out


def test_sbml_without_files_or_directory_fails(tmp_path: Path) -> None:
    empty = tmp_path / "empty"
    empty.mkdir()
    assert verify(["sbml", str(empty), str(empty)]) == 2
    assert verify(["sbml", str(tmp_path / "missing"), str(empty)]) == 2


def test_a_package_from_another_root_fails(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # pkpdsim_toy is imported from TESTS_DIR already, not from this root
    assert "pkpdsim_toy" in sys.modules
    other = tmp_path / "other"
    other.mkdir()
    out = tmp_path / "definitions.json"
    capsys.readouterr()
    code = verify(
        [
            "definitions",
            "--root",
            str(other),
            "--model",
            "pkpdsim_toy",
            "--out",
            str(out),
        ]
    )
    assert code == 2
    assert "not from below" in capsys.readouterr().err
    assert not out.exists()


def test_the_console_script_is_registered() -> None:
    from importlib.metadata import entry_points

    (script,) = entry_points(group="console_scripts", name="pkpdsim-verify")
    assert script.value == "pkpdsim.verify:main"


def test_an_unknown_model_package_fails() -> None:
    with pytest.raises(ModuleNotFoundError):
        describe("pkpdsim_missing")


def test_a_package_without_fit_definitions_is_described() -> None:
    info = describe("pkpdsim_legacy_nofit")
    assert info.fit_definitions is None
    assert info.experiments is None
    assert info.absolute_tolerance is None


def test_a_missing_nested_module_has_no_attribute() -> None:
    assert _attribute("pkpdsim_legacy_nofit.fitting.fitting", "FIT_DEFINITIONS") is None
    assert _attribute("pkpdsim_legacy_nofit.helpers", "ABSOLUTE_TOLERANCE") is None


def test_the_reraise_path_of_attribute_for_a_broken_module() -> None:
    # the module exists but its own import fails: not "no such attribute"
    with pytest.raises(ModuleNotFoundError, match="pkpdsim_does_not_exist"):
        _attribute("pkpdsim_legacy_broken.helpers", "ABSOLUTE_TOLERANCE")
