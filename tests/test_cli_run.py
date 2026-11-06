"""The run command of a model, a wrapper of the kinds of `python -m sbmlsim.run`."""

import json
from collections.abc import Callable
from pathlib import Path

import pytest
from sbmlsim.run import execute

import pkpdsim.cli
from pkpdsim.cli import (
    Action,
    parse_formats,
    parse_workers,
    project_reference,
    resolve_experiments,
    run_command,
    run_parser,
    run_workers,
)
from pkpdsim.model import RUN_ACTIONS, RUN_WORKERS, ExtraAction, PKPDModel
from pkpdsim_toy import MODEL
from pkpdsim_toy.experiments.decay import ToyDecay, ToyDecaySlow

GREET = ExtraAction(
    name="greet", help="write a greeting", function=".extras:greet", in_all=True
)
EXTRA = MODEL.replace(extra_actions=(GREET,), extra_arguments=".extras:add_arguments")
BROKEN = EXTRA.replace(experiments=".broken:BROKEN")

#: a serial run, the toy experiments start no pool
SERIAL = ("-w", "1")

type AsModel = Callable[[PKPDModel], PKPDModel]


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
        "--rebuild",
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


def test_the_commands_run_the_model_object_of_its_package(as_model: AsModel) -> None:
    assert project_reference(MODEL) == "pkpdsim_toy:MODEL"
    with pytest.raises(ValueError, match="another one"):
        project_reference(EXTRA)
    assert project_reference(as_model(EXTRA)) == "pkpdsim_toy:MODEL"


def test_empty_formats_are_a_usage_error() -> None:
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


def test_simulate_runs_the_kind_into_the_folder_of_the_group(tmp_path: Path) -> None:
    """The end-to-end run of `run_simulation`: the files of the old command."""
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
    out = tmp_path / "simulation" / "studies"
    files = {p.relative_to(out).as_posix() for p in out.rglob("*") if p.is_file()}
    assert files == {
        "ToyDecay/ToyDecay.html",
        "ToyDecay/ToyDecay.json",
        "ToyDecay/ToyDecay_fig_decay.png",
        "ToyDecay/datasets/ToyDecay_data_decay.tsv",
        "_figures/ToyDecay_fig_decay.png",
        "index.html",
        # the records of the run of the kind `simulate`
        "run.json",
        "spec.json",
        "summary.json",
        "progress.jsonl",
    }
    spec = json.loads((out / "spec.json").read_text(encoding="utf-8"))
    assert spec["project"] == "pkpdsim_toy:MODEL"
    assert spec["experiments"] == ["ToyDecay"]
    assert spec["figure_formats"] == ["png"]
    run = json.loads((out / "run.json").read_text(encoding="utf-8"))
    assert (run["kind"], run["status"]) == ("simulate", "done")


def test_a_run_replaces_the_earlier_results(tmp_path: Path) -> None:
    arguments = ["-a", "simulate", "-e", "studies", *SERIAL, "-r", str(tmp_path)]
    run_command(MODEL, arguments)
    stale = tmp_path / "simulation" / "studies" / "stale.txt"
    stale.write_text("an earlier result", encoding="utf-8")
    run_command(MODEL, arguments)
    assert not stale.exists()
    assert (tmp_path / "simulation" / "studies" / "index.html").is_file()


def test_unknown_experiments_are_left_out(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    run_command(
        MODEL,
        ["-a", "simulate", "-e", "ToyDecay,Unknown", *SERIAL, "-r", str(tmp_path)],
    )
    assert "Not found: Unknown" in capsys.readouterr().out
    spec = tmp_path / "simulation" / "custom_selection" / "spec.json"
    assert json.loads(spec.read_text(encoding="utf-8"))["experiments"] == ["ToyDecay"]


def test_a_selection_is_a_custom_selection(tmp_path: Path) -> None:
    run_command(
        MODEL,
        ["-a", "simulate", "-e", "ToyDecay,misc", *SERIAL, "-r", str(tmp_path)],
    )
    assert (tmp_path / "simulation" / "custom_selection" / "_figures").is_dir()


def test_a_relative_results_dir_is_relative_to_the_working_directory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    run_command(MODEL, ["-a", "simulate", "-e", "ToyDecay", *SERIAL, "-r", "out"])
    assert (tmp_path / "out" / "simulation" / "custom_selection").is_dir()


def test_the_factory_action(tmp_path: Path, as_model: AsModel) -> None:
    run_command(as_model(building_in(MODEL, tmp_path, tmp_path)), ["-a", "factory"])
    assert (tmp_path / "toy.xml").exists()


def test_a_rebuild_compares_the_models_with_the_package(tmp_path: Path) -> None:
    run_command(MODEL, ["-a", "factory", "--rebuild", "-r", str(tmp_path)])
    out = tmp_path / "factory"
    summary = json.loads((out / "summary.json").read_text(encoding="utf-8"))
    assert summary["changed"] == []
    assert (out / summary["directory"] / "toy.xml").is_file()


def test_a_rebuild_which_differs_exits_with_1(
    tmp_path: Path, as_model: AsModel, capsys: pytest.CaptureFixture[str]
) -> None:
    package = tmp_path / "package"
    model = as_model(
        building_in(MODEL, package, package).replace(
            variants={"male": package / "toy.xml"}, factory_dir=package
        )
    )
    run_command(model, ["-a", "factory"])
    (package / "toy.xml").write_text("<sbml/>", encoding="utf-8")
    with pytest.raises(SystemExit) as info:
        run_command(model, ["-a", "factory", "--rebuild", "-r", str(tmp_path / "r")])
    assert info.value.code == 1
    assert "toy.xml" in capsys.readouterr().out


def test_list_experiments(capsys: pytest.CaptureFixture[str]) -> None:
    run_command(MODEL, ["-a", "list_experiments"])
    output = capsys.readouterr().out
    assert "ToyDecaySlow" in output
    assert "misc" in output


def test_an_extra_action(tmp_path: Path, as_model: AsModel) -> None:
    run_command(
        as_model(EXTRA),
        ["-a", "greet", "--greeting", "hi", "--results-dir", str(tmp_path)],
    )
    assert (tmp_path / "greeting.txt").read_text(encoding="utf-8") == "hi toy"


def test_the_extra_actions_are_actions_of_the_project(tmp_path: Path) -> None:
    assert list(EXTRA.actions) == ["greet"]
    assert EXTRA.actions["greet"].description == "write a greeting"
    assert EXTRA.replace(title="other").actions == EXTRA.actions
    function = EXTRA.load_action("greet")
    function(EXTRA, tmp_path, workers=1)
    assert (tmp_path / "greeting.txt").read_text(encoding="utf-8") == "hello toy"
    function(EXTRA, tmp_path, workers=1, greeting="hi")
    assert (tmp_path / "greeting.txt").read_text(encoding="utf-8") == "hi toy"
    with pytest.raises(ValueError, match="extra actions"):
        MODEL.replace(actions=EXTRA.actions)


def test_the_kind_action_runs_an_extra_action(
    tmp_path: Path, as_model: AsModel
) -> None:
    as_model(EXTRA)
    spec = {
        "project": "pkpdsim_toy:MODEL",
        "name": "greet",
        "arguments": {"greeting": "hey"},
    }
    summary = execute("action", spec, tmp_path / "out")
    assert summary == {"action": "greet"}
    assert (tmp_path / "out" / "greeting.txt").read_text(encoding="utf-8") == "hey toy"


WORKERS = ExtraAction(
    name="workers", help="write the processes", function=".extras:workers"
)


@pytest.mark.parametrize(
    ("options", "n_workers"), [([], "1"), (["-w", "3"], "3"), (["-w", "auto"], "None")]
)
def test_an_extra_action_gets_the_simulator_of_the_run(
    tmp_path: Path, as_model: AsModel, options: list[str], n_workers: str
) -> None:
    model = as_model(MODEL.replace(extra_actions=(WORKERS,)))
    run_command(model, ["-a", "workers", *options, "--results-dir", str(tmp_path)])
    text = (tmp_path / "workers.txt").read_text(encoding="utf-8")
    assert text == f"toy {n_workers}"


def test_an_extra_action_of_a_model_without_experiments_gets_its_simulator(
    tmp_path: Path, as_model: AsModel
) -> None:
    model = as_model(MODEL.replace(experiments=None, fits={}, extra_actions=(WORKERS,)))
    run_command(model, ["-a", "workers", "--results-dir", str(tmp_path)])
    assert (tmp_path / "workers.txt").read_text(encoding="utf-8") == "toy 1"


def test_all_runs_the_factory_the_simulations_and_the_extras(
    tmp_path: Path, as_model: AsModel
) -> None:
    model = as_model(building_in(EXTRA, tmp_path, tmp_path / "models"))
    run_command(model, ["-a", "all", *SERIAL, "--results-dir", str(tmp_path)])
    assert (tmp_path / "models" / "toy.xml").exists()
    assert (tmp_path / "simulation" / "all" / "_figures").is_dir()
    assert (tmp_path / "greeting.txt").exists()


NO_EXPERIMENTS = MODEL.replace(experiments=None, fits={})


@pytest.mark.parametrize("action", ["simulate", "list_experiments"])
def test_a_model_without_experiments_cannot_simulate(action: str) -> None:
    with pytest.raises(SystemExit) as error:
        run_command(NO_EXPERIMENTS, ["-a", action])
    assert error.value.code == 2


def test_all_of_a_model_without_experiments_runs_the_factory(
    tmp_path: Path, as_model: AsModel
) -> None:
    model = as_model(building_in(NO_EXPERIMENTS, tmp_path, tmp_path / "models"))
    run_command(model, ["-a", "all", "-r", str(tmp_path)])
    assert (tmp_path / "models" / "toy.xml").exists()
    assert not (tmp_path / "simulation").exists()


def test_the_factory_of_a_model_without_experiments(
    tmp_path: Path, as_model: AsModel
) -> None:
    model = as_model(building_in(NO_EXPERIMENTS, tmp_path, tmp_path))
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


def test_all_of_a_model_without_experiments_runs_the_extras(
    tmp_path: Path, as_model: AsModel
) -> None:
    model = building_in(NO_EXPERIMENTS, tmp_path, tmp_path / "models").replace(
        extra_actions=(GREET,), extra_arguments=".extras:add_arguments"
    )
    run_command(as_model(model), ["-a", "all", "-r", str(tmp_path)])
    assert (tmp_path / "models" / "toy.xml").exists()
    assert (tmp_path / "greeting.txt").exists()
    assert not (tmp_path / "simulation").exists()


def test_the_workers_of_a_run() -> None:
    assert parse_workers("model") == RUN_WORKERS == 1
    assert parse_workers("auto") is None
    assert parse_workers("3") == 3
    # the run of the kind: every CPU of the quota for `model` and `auto`
    assert run_workers("model") is None
    assert run_workers("auto") is None
    assert run_workers("3") == 3


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


@pytest.mark.parametrize(("options", "workers"), [([], None), (["-w", "2"], 2)])
def test_the_workers_reach_the_run(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    options: list[str],
    workers: int | None,
) -> None:
    seen: list[int | None] = []
    original = pkpdsim.cli.execute

    def recording(kind: str, spec: dict, out: Path, **kwargs: object) -> dict:
        workers = kwargs.pop("workers", None)
        assert workers is None or isinstance(workers, int)
        seen.append(workers)
        return original(kind, spec, out, workers=workers, **kwargs)  # ty: ignore[invalid-argument-type]

    monkeypatch.setattr(pkpdsim.cli, "execute", recording)
    run_command(
        MODEL, ["-a", "simulate", "-e", "ToyDecay", *options, "-r", str(tmp_path)]
    )
    assert seen == [workers]


def test_a_failure_exits_with_1_after_every_action(
    tmp_path: Path, as_model: AsModel, capsys: pytest.CaptureFixture[str]
) -> None:
    model = as_model(building_in(BROKEN, tmp_path, tmp_path / "models"))
    with pytest.raises(SystemExit) as info:
        run_command(model, ["-a", "all", *SERIAL, "--results-dir", str(tmp_path)])
    assert info.value.code == 1
    # the report, the collected figures and the extra action are there
    assert list((tmp_path / "simulation" / "all").glob("**/*.html"))
    assert list((tmp_path / "simulation" / "all" / "_figures").glob("ToyDecay_*"))
    assert (tmp_path / "greeting.txt").exists()
    run = json.loads(
        (tmp_path / "simulation" / "all" / "run.json").read_text(encoding="utf-8")
    )
    assert run["status"] == "failed"
    out = capsys.readouterr().out
    assert "ToyCannotInitialize" in out
    assert "ToyRefusedUnit" in out
    assert "figures_mpl()" in out


def test_a_run_without_failures_exits_normally(tmp_path: Path) -> None:
    run_command(MODEL, ["-a", "simulate", "-e", "all", *SERIAL, "-r", str(tmp_path)])
    assert list((tmp_path / "simulation" / "all" / "_figures").glob("ToyDecay_*"))


def test_an_extra_option_named_simulator_is_an_error(
    tmp_path: Path, as_model: AsModel
) -> None:
    model = as_model(EXTRA.replace(extra_arguments=".extras:add_simulator_argument"))
    with pytest.raises(ValueError, match="--simulator"):
        run_command(model, ["-a", "greet", "--results-dir", str(tmp_path)])
