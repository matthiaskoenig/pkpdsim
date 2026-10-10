"""The verification of a change of a model, run on the toy model."""

import json
import shutil
from pathlib import Path

import pytest

from pkpdsim.verify import describe, verify
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


def test_dump_and_compare(tmp_path: Path) -> None:
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
    assert verify(["compare", str(old), str(new)]) == 0


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


def test_png(tmp_path: Path) -> None:
    output = MODEL.run_experiments(
        MODEL.load_experiments()["studies"], tmp_path / "a", figure_formats=["png"]
    )
    assert verify(["png", str(output), str(output)]) == 0


def test_the_console_script_is_registered() -> None:
    from importlib.metadata import entry_points

    (script,) = entry_points(group="console_scripts", name="pkpdsim-verify")
    assert script.value == "pkpdsim.verify:main"


def test_an_unknown_model_package_fails() -> None:
    with pytest.raises(ModuleNotFoundError):
        describe("pkpdsim_missing")
