"""`pkpdsim-verify`, the wrapper of `sbmlsim.verify` for the models of pkpdsim."""

import json
import pickle
from pathlib import Path

import pytest
from PIL import Image

from pkpdsim.verify import MODEL_COMMANDS, UsageError, translate, verify

TESTS_DIR = Path(__file__).parent


def _dump(tmp_path: Path, package: str, name: str = "dump.pkl") -> Path:
    out = tmp_path / name
    arguments = ["dump", "--root", str(TESTS_DIR), "--model", package]
    assert verify([*arguments, "--out", str(out)]) == 0
    return out


def test_the_console_script_is_registered() -> None:
    from importlib.metadata import entry_points

    (script,) = entry_points(group="console_scripts", name="pkpdsim-verify")
    assert script.value == "pkpdsim.verify:main"


def test_the_model_option_is_the_model_object_of_the_package() -> None:
    assert translate(["dump", "--root", "r", "--model", "pkg", "--out", "o"]) == [
        "dump",
        "--root",
        "r",
        "--project",
        "pkg:MODEL",
        "--out",
        "o",
    ]
    assert "--project=pkg:MODEL" in translate(
        ["settings", "--root", "r", "--model=pkg", "--out", "o"]
    )


@pytest.mark.parametrize("command", MODEL_COMMANDS)
def test_a_command_of_a_model_needs_the_root_and_the_model(command: str) -> None:
    with pytest.raises(UsageError, match="--model"):
        translate([command, "--root", "r", "--out", "o"])
    with pytest.raises(UsageError, match="--root"):
        translate([command, "--model", "pkg", "--out", "o"])
    assert verify([command, "--model", "pkg"]) == 2


def test_an_abbreviated_option_is_refused() -> None:
    """An abbreviation would pass the parser but not the translation."""
    assert verify(["dump", "--root", "r", "--mod", "pkg", "--out", "o"]) == 2


def test_the_scorecard_pools_the_ratios_by_health_by_default() -> None:
    assert translate(["scorecard", "records"]) == [
        "scorecard",
        "records",
        "--extensions",
        "pkpdsim.credibility:EXTENSIONS",
        "--pool",
        "ratios",
        "--field",
        "health",
    ]
    assert "health" not in translate(["scorecard", "r", "--field", "age"])


def test_dump_compare_and_the_builds_per_sex(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """A dump records the sex of the build of every task and the versions of pkpdsim."""
    old = _dump(tmp_path, "pkpdsim_sexes")
    dumped = pickle.loads(old.read_bytes())
    assert dumped["builds"] == {
        "SexesMen": {"task_decay": "male"},
        "SexesWomen": {"task_decay": "female"},
        "SexesGroups": {"task_men": "male", "task_women": "female"},
    }
    capsys.readouterr()
    assert verify(["compare", str(old), str(old)]) == 0
    assert "builds: 4 old, 4 new, 0 differ" in capsys.readouterr().out

    dumped["builds"]["SexesWomen"]["task_decay"] = "male"
    new = tmp_path / "new.pkl"
    new.write_bytes(pickle.dumps(dumped))
    assert verify(["compare", str(old), str(new)]) == 1
    assert "'female' -> 'male'" in capsys.readouterr().out


def test_the_variant_of_a_simulation_is_the_value_of_its_sex() -> None:
    from pkpdsim_sexes.experiments.decay import SexesGroups, SexesMen, SexesWomen

    assert SexesMen().variant_of("decay") == "male"
    assert SexesWomen().variant_of("decay") == "female"
    assert SexesGroups().variant_of("women") == "female"
    assert SexesGroups().variant_of("men") == "male"


def test_ranges_and_the_fit_definitions_of_the_toy_model(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    dumped = _dump(tmp_path, "pkpdsim_toy")
    versions = pickle.loads(dumped.read_bytes())["versions"]
    assert {"sbmlsim", "pkpdsim", "sbmlutils"} <= set(versions)
    assert verify(["ranges", str(dumped)]) == 0
    definitions = tmp_path / "definitions.json"
    arguments = ["--root", str(TESTS_DIR), "--model", "pkpdsim_toy"]
    assert verify(["definitions", *arguments, "--out", str(definitions)]) == 0
    assert "DECAY" in json.loads(definitions.read_text(encoding="utf-8"))
    settings = tmp_path / "settings.json"
    assert verify(["settings", *arguments, "--out", str(settings)]) == 0
    assert settings.exists()


def test_a_minimal_fit_of_the_toy_model(tmp_path: Path) -> None:
    arguments = ["--root", str(TESTS_DIR), "--model", "pkpdsim_toy"]
    assert verify(["fits", *arguments, "--out", str(tmp_path / "fits")]) == 0


def test_png_and_sbml_compare_files(tmp_path: Path) -> None:
    from pkpdsim_toy import MODEL_PATH

    assert verify(["sbml", str(MODEL_PATH), str(MODEL_PATH)]) == 0
    (old, new) = (tmp_path / "old", tmp_path / "new")
    old.mkdir()
    new.mkdir()
    for directory in (old, new):
        Image.new("RGB", (4, 4), "white").save(directory / "figure.png")
    assert verify(["png", str(old), str(new)]) == 0


def test_a_package_without_a_model_object_is_a_configuration_error(
    tmp_path: Path,
) -> None:
    code = verify(
        [
            "dump",
            "--root",
            str(TESTS_DIR),
            "--model",
            "pkpdsim_toy.experiments",
            "--out",
            str(tmp_path / "x.pkl"),
        ]
    )
    assert code == 4


@pytest.mark.parametrize(
    "command", [None, "dump", "compare", "fits", "definitions", "settings", "scorecard"]
)
def test_the_help_shows_the_options_of_pkpdsim(
    command: str | None, capsys: pytest.CaptureFixture[str]
) -> None:
    assert verify([command, "--help"] if command else ["--help"]) == 0
    text = capsys.readouterr().out
    assert "pkpdsim-verify" in text
    assert "--project" not in text
    assert "sbmlsim.run" not in text
    if command in MODEL_COMMANDS:
        assert "--model PACKAGE" in text


def test_a_scorecard_without_pooled_rows_prints_no_qualification(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    from sbmlsim.credibility.scorecard import load_records
    from test_cli_credibility import _assessment, _write

    records = _write(tmp_path / "records", _assessment("a"))
    assert len(load_records([records])) == 1
    assert verify(["scorecard", str(records), "--field", "nothing"]) == 0
    assert "qualification by" not in capsys.readouterr().out
