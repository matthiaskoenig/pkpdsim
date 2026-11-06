"""The run command of a model."""

from pathlib import Path

import pytest
from sbmlsim.project import SimulatorSettings

from pkpdsim.cli import (
    Action,
    parse_formats,
    parse_workers,
    resolve_experiments,
    run_command,
    run_parser,
)
from pkpdsim.model import RUN_ACTIONS, RUN_WORKERS, ExtraAction, PKPDModel
from pkpdsim_toy import MODEL
from pkpdsim_toy.experiments.decay import ToyDecay, ToyDecaySlow

GREET = ExtraAction(
    name="greet", help="write a greeting", function=".extras:greet", in_all=True
)
EXTRA = MODEL.replace(extra_actions=(GREET,), extra_arguments=".extras:add_arguments")
BROKEN = EXTRA.replace(experiments=".broken:BROKEN")


def building_in(model: PKPDModel, root: Path, factory_dir: Path) -> PKPDModel:
    """Get a model whose factory builds into a temporary directory.

    The directory of the factory lies inside the root of a project, so the
    model moves its root there and keeps the base path of its experiments.
    """
    return model.replace(root=root, base_path=model.root, factory_dir=factory_dir)


def test_the_actions_of_the_model_and_the_command_agree() -> None:
    assert {action.value for action in Action} == set(RUN_ACTIONS)


def test_the_options_of_the_command() -> None:
    text = run_parser(MODEL).format_help()
    for option in [
        "--action",
        "--results-dir",
        "--experiments",
        "--open",
        "--formats",
    ]:
        assert option in text
    assert "usage: run_simulation pkpdsim_toy " in text
    assert "svg,png" in text


def test_the_extra_options_and_actions() -> None:
    text = run_parser(EXTRA).format_help()
    assert "--greeting" in text
    assert "greet: write a greeting" in text


def test_resolve_groups_and_experiments() -> None:
    registry = MODEL.load_experiments()
    classes, not_found = resolve_experiments(registry, ["studies", "ToyDecaySlow"])
    assert classes == [ToyDecay, ToyDecaySlow]
    assert not_found == []


def test_resolve_drops_duplicates() -> None:
    registry = MODEL.load_experiments()
    classes, not_found = resolve_experiments(registry, ["all", "ToyDecay", "x"])
    assert classes == [ToyDecay, ToyDecaySlow]
    assert not_found == ["x"]


def test_formats_are_normalized() -> None:
    assert parse_formats(" PNG, svg ,") == ["png", "svg"]


def test_empty_formats_are_a_usage_error(tmp_path: Path) -> None:
    with pytest.raises(SystemExit) as error:
        run_command(MODEL, ["-a", "list_experiments", "--formats", ","])
    assert error.value.code == 2


def test_simulate_needs_experiments() -> None:
    for experiments in [[], ["--experiments", " , "]]:
        with pytest.raises(SystemExit) as error:
            run_command(MODEL, ["-a", "simulate", *experiments])
        assert error.value.code == 2


def test_unknown_experiments_only_exit_with_1(tmp_path: Path) -> None:
    with pytest.raises(SystemExit) as error:
        run_command(
            MODEL,
            ["-a", "simulate", "-e", "Unknown", "--results-dir", str(tmp_path)],
        )
    assert error.value.code == 1


def test_simulate_writes_into_the_results_dir(tmp_path: Path) -> None:
    run_command(
        MODEL,
        [
            "-a",
            "simulate",
            "-e",
            "studies",
            "--results-dir",
            str(tmp_path),
            "--formats",
            "png",
        ],
    )
    figures = tmp_path / "simulation" / "studies" / "_figures"
    assert list(figures.glob("ToyDecay_*.png"))
    assert not list(figures.glob("*.svg"))


def test_a_selection_is_a_custom_selection(tmp_path: Path) -> None:
    run_command(
        MODEL,
        ["-a", "simulate", "-e", "ToyDecay,misc", "--results-dir", str(tmp_path)],
    )
    assert (tmp_path / "simulation" / "custom_selection" / "_figures").is_dir()


def test_a_relative_results_dir_is_relative_to_the_working_directory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    run_command(MODEL, ["-a", "simulate", "-e", "ToyDecay", "-r", "out"])
    assert (tmp_path / "out" / "simulation" / "custom_selection").is_dir()


def test_the_factory_action(tmp_path: Path) -> None:
    run_command(building_in(MODEL, tmp_path, tmp_path), ["-a", "factory"])
    assert (tmp_path / "toy.xml").exists()


def test_list_experiments(capsys: pytest.CaptureFixture[str]) -> None:
    run_command(MODEL, ["-a", "list_experiments"])
    output = capsys.readouterr().out
    assert "ToyDecaySlow" in output


def test_an_extra_action(tmp_path: Path) -> None:
    run_command(
        EXTRA, ["-a", "greet", "--greeting", "hi", "--results-dir", str(tmp_path)]
    )
    assert (tmp_path / "greeting.txt").read_text(encoding="utf-8") == "hi toy"


WORKERS = ExtraAction(
    name="workers", help="write the processes", function=".extras:workers"
)


@pytest.mark.parametrize(
    ("options", "n_workers"), [([], "1"), (["-w", "3"], "3"), (["-w", "auto"], "None")]
)
def test_an_extra_action_gets_the_simulator_of_the_run(
    tmp_path: Path, options: list[str], n_workers: str
) -> None:
    model = MODEL.replace(extra_actions=(WORKERS,))
    run_command(model, ["-a", "workers", *options, "--results-dir", str(tmp_path)])
    text = (tmp_path / "workers.txt").read_text(encoding="utf-8")
    assert text == f"toy {n_workers}"


def test_an_extra_action_of_a_model_without_experiments_gets_its_simulator(
    tmp_path: Path,
) -> None:
    model = MODEL.replace(experiments=None, fits={}, extra_actions=(WORKERS,))
    run_command(model, ["-a", "workers", "--results-dir", str(tmp_path)])
    assert (tmp_path / "workers.txt").read_text(encoding="utf-8") == "toy 1"


def test_all_runs_the_factory_the_simulations_and_the_extras(tmp_path: Path) -> None:
    model = building_in(EXTRA, tmp_path, tmp_path / "models")
    run_command(model, ["-a", "all", "--results-dir", str(tmp_path)])
    assert (tmp_path / "models" / "toy.xml").exists()
    assert (tmp_path / "simulation" / "all" / "_figures").is_dir()
    assert (tmp_path / "greeting.txt").exists()


NO_EXPERIMENTS = MODEL.replace(experiments=None, fits={})


@pytest.mark.parametrize("action", ["simulate", "list_experiments"])
def test_a_model_without_experiments_cannot_simulate(action: str) -> None:
    with pytest.raises(SystemExit) as error:
        run_command(NO_EXPERIMENTS, ["-a", action])
    assert error.value.code == 2


def test_all_of_a_model_without_experiments_runs_the_factory(tmp_path: Path) -> None:
    model = building_in(NO_EXPERIMENTS, tmp_path, tmp_path / "models")
    run_command(model, ["-a", "all", "-r", str(tmp_path)])
    assert (tmp_path / "models" / "toy.xml").exists()
    assert not (tmp_path / "simulation").exists()


def test_the_factory_of_a_model_without_experiments(tmp_path: Path) -> None:
    model = building_in(NO_EXPERIMENTS, tmp_path, tmp_path)
    run_command(model, ["-a", "factory"])
    assert (tmp_path / "toy.xml").exists()


def test_the_help_of_a_model_without_experiments_offers_what_applies() -> None:
    text = " ".join(run_parser(NO_EXPERIMENTS).format_help().split())
    assert "Run the model factory of the" in text
    assert "simulation experiments" not in text
    for absent in [
        "--experiments",
        "--formats",
        "--open",
        "simulate",
        "list_experiments",
    ]:
        assert absent not in text
    assert "{factory,all}" in text


def test_the_help_of_a_model_without_experiments_names_its_extra_actions() -> None:
    model = NO_EXPERIMENTS.replace(extra_actions=(GREET,))
    text = " ".join(run_parser(model).format_help().split())
    assert "{factory,all,greet}" in text
    assert "greet: write a greeting" in text
    assert "Run the model factory of the" in text
    assert "and its extra actions." in text
    assert "--experiments" not in text


def test_all_of_a_model_without_experiments_runs_the_extras(tmp_path: Path) -> None:
    model = building_in(NO_EXPERIMENTS, tmp_path, tmp_path / "models").replace(
        extra_actions=(GREET,), extra_arguments=".extras:add_arguments"
    )
    run_command(model, ["-a", "all", "-r", str(tmp_path)])
    assert (tmp_path / "models" / "toy.xml").exists()
    assert (tmp_path / "greeting.txt").exists()
    assert not (tmp_path / "simulation").exists()


def test_the_workers_of_a_run() -> None:
    assert parse_workers("model") == RUN_WORKERS == 1
    assert parse_workers("auto") is None
    assert parse_workers("3") == 3


@pytest.mark.parametrize("text", ["0", "-2", "two", "1.5", ""])
def test_the_workers_are_a_positive_number(text: str) -> None:
    with pytest.raises(ValueError, match="--workers"):
        parse_workers(text)


def test_a_bad_number_of_workers_is_a_usage_error(tmp_path: Path) -> None:
    with pytest.raises(SystemExit) as info:
        run_command(
            MODEL, ["-a", "simulate", "-e", "ToyDecay", "-w", "0", "-r", str(tmp_path)]
        )
    assert info.value.code == 2


def test_the_workers_reach_the_run(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    seen: list[tuple[SimulatorSettings, int | None]] = []
    original = SimulatorSettings.simulator

    def create(settings: SimulatorSettings, n_workers: int | None) -> object:
        seen.append((settings, n_workers))
        return original(settings, n_workers)

    monkeypatch.setattr(SimulatorSettings, "simulator", create)
    run_command(
        MODEL, ["-a", "simulate", "-e", "ToyDecay", "-w", "2", "-r", str(tmp_path)]
    )
    assert [n for _, n in seen] == [2]
    assert seen[0][0] == MODEL.simulator


def test_a_failure_exits_with_1_after_every_action(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    model = building_in(BROKEN, tmp_path, tmp_path / "models")
    with pytest.raises(SystemExit) as info:
        run_command(model, ["-a", "all", "--results-dir", str(tmp_path)])
    assert info.value.code == 1
    # the report, the collected figures and the extra action are there
    assert list((tmp_path / "simulation" / "all").glob("**/*.html"))
    assert list((tmp_path / "simulation" / "all" / "_figures").glob("ToyDecay_*"))
    assert (tmp_path / "greeting.txt").exists()
    out = capsys.readouterr().out
    assert "ToyCannotInitialize" in out
    assert "ToyRefusedUnit" in out
    assert "figures_mpl()" in out


def test_a_run_without_failures_exits_normally(tmp_path: Path) -> None:
    run_command(MODEL, ["-a", "simulate", "-e", "all", "-r", str(tmp_path)])
    assert list((tmp_path / "simulation" / "all" / "_figures").glob("ToyDecay_*"))


def test_an_extra_option_named_simulator_is_an_error(tmp_path: Path) -> None:
    model = EXTRA.replace(extra_arguments=".extras:add_simulator_argument")
    with pytest.raises(ValueError, match="--simulator"):
        run_command(model, ["-a", "greet", "--results-dir", str(tmp_path)])
