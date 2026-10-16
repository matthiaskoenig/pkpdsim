"""The console scripts of the commands, which take the model package first."""

from collections.abc import Callable, Sequence
from importlib.metadata import distribution
from pathlib import Path

import pytest

from pkpdsim.cli import (
    COMMANDS,
    ModelNotFoundError,
    resolve_model,
    run_fit_main,
    run_identifiability_main,
    run_report_main,
    run_simulation_main,
)
from pkpdsim_toy import MODEL

Main = Callable[[Sequence[str] | None], None]

MAINS: dict[str, Main] = {
    "run_simulation": run_simulation_main,
    "run_fit": run_fit_main,
    "run_report": run_report_main,
    "run_identifiability": run_identifiability_main,
}


@pytest.fixture
def packages(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Packages which are no model: no `MODEL`, a `MODEL` of another type, ...."""
    root = tmp_path / "packages"
    files = {
        "nomodel/__init__.py": "",
        "wrongmodel/__init__.py": "MODEL = 42\n",
        "brokenmodel/__init__.py": "import missing_dependency_of_the_model\n",
        "aliasmodel/__init__.py": "from pkpdsim_toy import MODEL\n",
        "nofits/__init__.py": (
            "import dataclasses\n"
            "from pkpdsim_toy import MODEL as TOY\n"
            "MODEL = dataclasses.replace(TOY, package=__name__, fit_definitions=None)\n"
        ),
    }
    for name, text in files.items():
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
    monkeypatch.syspath_prepend(str(root))
    return root


def test_the_commands_are_console_scripts() -> None:
    scripts = {
        entry.name: entry.value
        for entry in distribution("pkpdsim").entry_points
        if entry.group == "console_scripts"
    }
    for name in MAINS:
        assert scripts[name] == f"pkpdsim.cli:{name}_main"
    assert set(COMMANDS) == set(MAINS)


def test_a_model_package_resolves_to_its_model_object() -> None:
    assert resolve_model("pkpdsim_toy") is MODEL


@pytest.mark.parametrize(
    ("package", "message"),
    [
        ("pkpdsim_no_such_package", "is not found"),
        ("pkpdsim_no_such_package.models.toy", "is not found"),
        ("pkpdsim_toy.no_such_module", "is not found"),
        ("nomodel", "has no model object MODEL"),
        ("wrongmodel", "is no PKPDModel"),
        ("aliasmodel", "belongs to the package 'pkpdsim_toy'"),
        ("", "is no dotted package"),
        ("pkpdsim_toy/", "is no dotted package"),
        (".pkpdsim_toy", "is no dotted package"),
    ],
)
def test_a_package_which_is_no_model_is_not_found(
    packages: Path, package: str, message: str
) -> None:
    with pytest.raises(ModelNotFoundError, match=message):
        resolve_model(package)


def test_an_import_error_inside_a_model_keeps_its_traceback(packages: Path) -> None:
    with pytest.raises(ModuleNotFoundError, match="missing_dependency_of_the_model"):
        resolve_model("brokenmodel")


@pytest.mark.parametrize("name", MAINS)
def test_the_help_without_a_model_names_the_package(
    name: str, capsys: pytest.CaptureFixture[str]
) -> None:
    with pytest.raises(SystemExit) as exit_info:
        MAINS[name](["--help"])
    assert exit_info.value.code == 0
    out = capsys.readouterr().out
    assert out.startswith(f"usage: {name} PACKAGE [options]")
    assert f"{name} PACKAGE --help" in out


@pytest.mark.parametrize("name", MAINS)
def test_the_help_of_a_model_is_the_help_of_its_command(
    name: str, capsys: pytest.CaptureFixture[str]
) -> None:
    with pytest.raises(SystemExit) as exit_info:
        MAINS[name](["pkpdsim_toy", "--help"])
    assert exit_info.value.code == 0
    assert capsys.readouterr().out.startswith(f"usage: {name} pkpdsim_toy ")


@pytest.mark.parametrize(
    ("args", "message"),
    [
        ([], "the model package is required"),
        (["-a", "factory", "pkpdsim_toy"], "the model package is the first argument"),
        (["pkpdsim_no_such_package"], "is not found"),
        (["nomodel"], "has no model object MODEL"),
    ],
)
@pytest.mark.parametrize("name", MAINS)
def test_a_command_without_a_model_is_a_usage_error(
    packages: Path,
    name: str,
    args: list[str],
    message: str,
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit) as exit_info:
        MAINS[name](args)
    assert exit_info.value.code == 2
    err = capsys.readouterr().err
    assert err.startswith(f"usage: {name} PACKAGE [options]")
    assert message in err


@pytest.mark.parametrize("name", ["run_fit", "run_report", "run_identifiability"])
def test_the_fit_commands_need_fit_definitions(
    packages: Path, name: str, capsys: pytest.CaptureFixture[str]
) -> None:
    with pytest.raises(SystemExit) as exit_info:
        MAINS[name](["nofits", "--help"])
    assert exit_info.value.code == 2
    assert "the model 'toy' (nofits) has no fit definitions" in capsys.readouterr().err


def test_the_simulation_of_a_model_without_fit_definitions_runs(
    packages: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    with pytest.raises(SystemExit) as exit_info:
        run_simulation_main(["nofits", "--help"])
    assert exit_info.value.code == 0


def test_run_simulation_runs_the_experiments(tmp_path: Path) -> None:
    run_simulation_main(
        ["pkpdsim_toy", "-a", "simulate", "-e", "ToyDecay", "-r", str(tmp_path)]
    )
    assert (tmp_path / "simulation" / "custom_selection").is_dir()


def test_the_fit_commands_run_the_fit_and_its_reports(tmp_path: Path) -> None:
    fit_dir = tmp_path / "fit"
    run_fit_main(
        [
            "pkpdsim_toy",
            "--runs",
            "1",
            "--no-mapping-figures",
            "--output_dir",
            str(fit_dir),
        ]
    )
    (parameters,) = fit_dir.glob("**/parameters.json")
    run_report_main(
        [
            "pkpdsim_toy",
            "--subset",
            "DECAY",
            "--no-mapping-figures",
            "--output_dir",
            str(tmp_path / "report"),
            str(parameters),
        ]
    )
    assert list((tmp_path / "report").glob("**/index.html"))
    run_identifiability_main(
        [
            "pkpdsim_toy",
            "--subset",
            "DECAY",
            "--max-points",
            "3",
            "--no-fisher",
            "--output_dir",
            str(tmp_path / "identifiability"),
            str(parameters),
        ]
    )
    assert list((tmp_path / "identifiability").glob("**/identifiability.json"))
