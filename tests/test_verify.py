"""The verification of a change of a model, run on the toy model."""

import copy
import dataclasses
import importlib.metadata
import json
import math
import pickle
import re
import shutil
import sys
from collections.abc import Callable, Sequence
from pathlib import Path
from types import MappingProxyType

import numpy as np
import pytest
from figure_builders import FigureStudy, build_figures
from PIL import Image
from sbmlsim.experiment import SimulationExperiment
from sbmlsim.plot import Figure
from sbmlsim.simulation import Scan, Simulation
from sbmlsim.simulation.sampling import local

from pkpdsim.model import FigureSettings
from pkpdsim.verify import (
    FIGURE_SETTINGS_08,
    NEW_FIGURE_SETTINGS,
    Dump,
    align_figures,
    apply_renames,
    describe,
    figure_spec,
    load,
    load_renames,
    matched_rows,
    rename_mapping,
    verify,
)
from pkpdsim_toy import DATA_PATHS, MODEL, MODEL_PATH, TOY_PATH
from pkpdsim_toy.experiments.decay import ToyDecay

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


def test_dump_and_compare(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    common = ["--root", str(TESTS_DIR), "--atol", "model", "--rtol", "1e-10"]
    old, new = tmp_path / "old.pkl", tmp_path / "new.pkl"
    assert verify(["dump", "--model", "pkpdsim_toy", "--out", str(old), *common]) == 0
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


def test_capture_settings_of_a_model() -> None:
    from pkpdsim.verify import capture_settings

    settings = capture_settings(MODEL.run_experiments, ToyDecay, TESTS_DIR)
    assert settings["experiment"] == "ToyDecay"
    assert settings["simulator"]["relative_tolerance"] == 1e-10
    assert settings["run"]["figure_formats"] == ["svg", "png"]
    assert settings["run"]["reduced_selections"] is True
    assert settings["run"]["save_results"] is False
    assert settings["run"]["figure_settings"]["legend_fontsize"] == 7
    assert settings["rc"]["savefig.dpi"] == 30
    assert settings["calls"] == {"simulator": 1, "runner": 1, "run": 1}
    assert "<root>" in settings["runner"]["base_path"]
    assert settings["simulator"]["n_workers"] == 1
    assert "figure" not in settings


def test_capture_settings_sees_a_difference() -> None:
    from pkpdsim.verify import capture_settings

    other = dataclasses.replace(
        MODEL, figures=FigureSettings(formats=("png",), legend_fontsize=9)
    )
    a = capture_settings(MODEL.run_experiments, ToyDecay, TESTS_DIR)
    b = capture_settings(other.run_experiments, ToyDecay, TESTS_DIR)
    assert a["run"]["figure_settings"] != b["run"]["figure_settings"]
    assert a["run"]["figure_formats"] != b["run"]["figure_formats"]


def test_dump_with_the_relative_tolerance_of_the_model(tmp_path: Path) -> None:
    out = tmp_path / "dump.pkl"
    code = verify(
        [
            "dump",
            "--root",
            str(TESTS_DIR),
            "--model",
            "pkpdsim_toy",
            "--out",
            str(out),
            "--atol",
            "model",
            "--rtol",
            "model",
        ]
    )
    assert code == 0


def test_a_dump_records_an_experiment_which_cannot_be_initialized(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    import pkpdsim_toy

    monkeypatch.setattr(
        pkpdsim_toy, "MODEL", dataclasses.replace(MODEL, experiments=".broken:BROKEN")
    )
    # the experiment which cannot be initialized runs first
    order = tmp_path / "order.pkl"
    with order.open("wb") as f:
        pickle.dump({"order": ["ToyCannotInitialize", "ToyDecay", "ToyFigureFails"]}, f)
    out = tmp_path / "dump.pkl"
    code = verify(
        [
            "dump",
            "--root",
            str(TESTS_DIR),
            "--model",
            "pkpdsim_toy",
            "--out",
            str(out),
            "--order",
            str(order),
        ]
    )
    assert code == 1
    assert "ToyCannotInitialize: ok" not in capsys.readouterr().err
    dump = load(str(out))
    # an experiment which the order does not name runs last
    assert dump["order"] == [
        "ToyCannotInitialize",
        "ToyDecay",
        "ToyFigureFails",
        "ToyRefusedUnit",
    ]
    assert list(dump["errors"]) == ["ToyCannotInitialize", "ToyRefusedUnit"]
    assert (
        "RuntimeError: the simulations of ToyCannotInitialize fail"
        in dump["errors"]["ToyCannotInitialize"]
    )
    assert "ValueError: The unit 'IU_toy_refused'" in dump["errors"]["ToyRefusedUnit"]
    # the experiments after it ran and are stored
    assert {sid for sid, _, _ in dump["results"]} == {"ToyDecay", "ToyFigureFails"}
    assert ("ToyDecay", "task_decay", "[A]") in dump["results"]
    assert set(dump["mappings"]) == {"ToyDecay", "ToyFigureFails"}


def test_a_dump_stores_the_time_of_the_grid(tmp_path: Path) -> None:
    out = tmp_path / "dump.pkl"
    assert (
        verify(
            [
                "dump",
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
    dims, values = load(str(out))["results"]["ToyDecay", "task_decay", "time"]
    assert dims == ("time",)
    assert values[0] == 0.0
    assert values.size == 51
    dims, _ = load(str(out))["results"]["ToyDecay", "task_decay", "[A]"]
    assert dims == ("time",)


class ToyDesign(ToyDecay):
    """The decay with a local design of k, whose creation constructs simulators of its own."""

    # the toy experiments have no scan, this one adds its design
    def simulations(self) -> dict[str, Simulation | Scan]:  # ty: ignore[invalid-method-override]
        """The decay and its local design."""
        sim = Simulation(end=10, steps=50)
        return {
            "decay": sim,
            "decay_sensitivity": Scan(
                sim,
                [local(["k"], delta=0.5, model=self._models["model"], simulation=sim)],
            ),
        }


def test_capture_settings_of_an_experiment_with_a_design() -> None:
    from pkpdsim.verify import capture_settings

    settings = capture_settings(MODEL.run_experiments, ToyDesign, TESTS_DIR)
    assert settings["simulator"]["n_workers"] == 1
    assert settings["simulator"]["absolute_tolerance"] == 1e-10
    assert settings["calls"]["simulator"] == 1


def _fake_run(positional: bool) -> "Callable[..., object]":
    from sbmlsim.experiment import ExperimentRunner
    from sbmlsim.simulator import Simulator

    def run(experiment: type[SimulationExperiment], output_dir: Path) -> None:
        if positional:
            simulator = Simulator(1, absolute_tolerance=1e-10)
            runner = ExperimentRunner([experiment], TOY_PATH, DATA_PATHS, simulator)
            runner.run_experiments(output_dir, False, False, ["png"])
        else:
            simulator = Simulator(n_workers=1, absolute_tolerance=1e-10)
            runner = ExperimentRunner(
                experiment_classes=[experiment],
                simulator=simulator,
                data_path=DATA_PATHS,
                base_path=TOY_PATH,
            )
            runner.run_experiments(
                output_path=output_dir,
                show_figures=False,
                save_results=False,
                figure_formats=["png"],
            )

    return run


def test_capture_settings_binds_positional_arguments() -> None:
    from pkpdsim.verify import capture_settings

    a = capture_settings(_fake_run(True), ToyDecay, TESTS_DIR)
    b = capture_settings(_fake_run(False), ToyDecay, TESTS_DIR)
    assert a == b
    assert a["simulator"]["absolute_tolerance"] == 1e-10
    assert a["run"]["figure_formats"] == ["png"]
    assert a["runner"]["experiment_classes"] == ["ToyDecay"]


def test_capture_settings_requires_a_complete_run(tmp_path: Path) -> None:
    from pkpdsim.verify import capture_settings

    with pytest.raises(RuntimeError, match="did not reach sbmlsim"):
        capture_settings(lambda experiment, output_dir: None, ToyDecay, TESTS_DIR)


def test_capture_settings_counts_the_calls() -> None:
    from sbmlsim.experiment import ExperimentRunner

    from pkpdsim.verify import capture_settings

    inner = _fake_run(False)

    def twice(experiment: type[SimulationExperiment], output_dir: Path) -> None:
        inner(experiment, output_dir)
        ExperimentRunner([experiment], TOY_PATH, DATA_PATHS)

    settings = capture_settings(twice, ToyDecay, TESTS_DIR)
    assert settings["calls"] == {"simulator": 1, "runner": 2, "run": 1}


def test_settings_of_a_run_which_misses_sbmlsim_exit_with_2(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(MODEL.__class__, "run_experiments", lambda self, *a, **k: None)
    code = verify(
        [
            "settings",
            "--root",
            str(TESTS_DIR),
            "--model",
            "pkpdsim_toy",
            "--out",
            str(tmp_path / "s.json"),
        ]
    )
    assert code == 2


def test_capture_settings_records_the_matplotlib_settings() -> None:
    import matplotlib

    from pkpdsim.verify import capture_settings

    other = dataclasses.replace(MODEL, figures=FigureSettings(savefig_dpi=55))
    assert (
        capture_settings(other.run_experiments, ToyDecay, TESTS_DIR)["rc"][
            "savefig.dpi"
        ]
        == 55
    )

    def with_font(experiment: type[SimulationExperiment], output_dir: Path) -> None:
        with matplotlib.rc_context({"font.size": 3}):
            MODEL.run_experiments(experiment, output_dir)

    rc = capture_settings(with_font, ToyDecay, TESTS_DIR)["rc"]
    assert rc["font.size"] == 3
    assert rc["savefig.dpi"] == 30


def test_plain_converts_nested_paths() -> None:
    from pkpdsim.verify import _plain

    value = _plain({"a": [{"b": TESTS_DIR / "x"}], "c": (TESTS_DIR,)})
    assert value == {
        "a": [{"b": (TESTS_DIR / "x").as_posix()}],
        "c": [TESTS_DIR.as_posix()],
    }


def test_dump_with_a_bad_tolerance_exits_with_2(tmp_path: Path) -> None:
    for option in ("--atol", "--rtol"):
        code = verify(
            [
                "dump",
                "--root",
                str(TESTS_DIR),
                "--model",
                "pkpdsim_toy",
                "--out",
                str(tmp_path / "d.pkl"),
                option,
                "abc",
            ]
        )
        assert code == 2


def test_a_package_without_a_model_object_is_a_usage_error(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    with pytest.raises(ValueError, match="defines no model object"):
        describe("pkpdsim_toy.experiments")
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
    assert code == 2
    assert "pkpdsim_toy.experiments defines no model object" in capsys.readouterr().err


def test_the_metadata_of_a_mapping_is_renamed() -> None:
    description = (
        "fm|d|t|time|A|4|None|None|t|t|time|[A]|None|None|None|weight=None|"
        + str([("data", "MEAN"), ("dosing", "MULTI"), ("route", "IV")])
    )
    renamed = rename_mapping(
        description,
        {
            "fields": {"data": "data_type"},
            "members": {"dosing": {"MULTI": "MULTIPLE"}},
            "added": {"application_form": "NR"},
        },
    )
    assert renamed.endswith(
        str(
            [
                ("application_form", "NR"),
                ("data_type", "MEAN"),
                ("dosing", "MULTIPLE"),
                ("route", "IV"),
            ]
        )
    )


def test_compare_applies_the_expected_renames(tmp_path: Path) -> None:
    old, new = tmp_path / "old.pkl", tmp_path / "new.pkl"
    args = ["--root", str(TESTS_DIR), "--model", "pkpdsim_toy"]
    assert verify(["dump", *args, "--out", str(old)]) == 0
    dump = load(str(old))
    dump["mappings"] = {
        sid: [d.replace("('route', 'IV')", "('route', 'IVX')") for d in descriptions]
        for sid, descriptions in dump["mappings"].items()
    }
    dump["figures"] = {
        sid: {key: spec.replace('"A"', '"A renamed"') for key, spec in figures.items()}
        for sid, figures in dump["figures"].items()
    }
    new.write_bytes(pickle.dumps(dump))
    renames = tmp_path / "renames.json"
    assert verify(["compare", str(old), str(new)]) == 1
    renames.write_text(
        json.dumps(
            {
                "metadata": {"members": {"route": {"IV": "IVX"}}},
                "strings": {'"A"': '"A renamed"'},
            }
        )
    )
    assert verify(["compare", str(old), str(new), "--renames", str(renames)]) == 0


def test_a_renames_file_with_unknown_parts_is_an_error(tmp_path: Path) -> None:
    renames = tmp_path / "renames.json"
    renames.write_text(json.dumps({"metadata": {"member": {}}}))
    with pytest.raises(ValueError, match="member"):
        load_renames(str(renames))


@pytest.mark.parametrize(
    ("content", "part"),
    [
        ([], "object"),
        ({"strings": ["a"]}, "strings"),
        ({"strings": {"a": 1}}, "strings"),
        ({"metadata": []}, "metadata"),
        ({"metadata": {"fields": []}}, "fields"),
        ({"metadata": {"fields": {"a": 1}}}, "fields"),
        ({"metadata": {"members": {"route": "PO"}}}, "members"),
        ({"metadata": {"members": {"route": {"a": 1}}}}, "members"),
        ({"metadata": {"added": {"a": 1}}}, "added"),
    ],
)
def test_a_renames_file_of_the_wrong_type_is_an_error(
    tmp_path: Path, content: object, part: str
) -> None:
    renames = tmp_path / "renames.json"
    renames.write_text(json.dumps(content))
    with pytest.raises(ValueError, match=part):
        load_renames(str(renames))


def test_renamed_mappings_are_sorted_like_a_dump() -> None:
    old = {
        "mappings": {
            "E": [
                "fm|a|t|time|A|4|None|None|t|t|time|[A]|None|None|None|weight=None|"
                + str([("route", "B")]),
                "fm|a|t|time|A|4|None|None|t|t|time|[A]|None|None|None|weight=None|"
                + str([("route", "C")]),
            ]
        }
    }
    renames = {"metadata": {"members": {"route": {"B": "Z"}}}}
    renamed = apply_renames(old, renames)  # ty: ignore[invalid-argument-type]
    assert renamed["mappings"]["E"] == sorted(renamed["mappings"]["E"])
    assert renamed["mappings"]["E"][0].endswith(str([("route", "C")]))


T_OLD = np.array([0.0, 1.0, 2.0, 3.0])
#: a change at 2: the state before it, then the state after it
T_NEW = np.array([0.0, 1.0, 2.0, 2.0, 3.0])


def synthetic(
    results: dict | None = None,
    *,
    datasets: dict | None = None,
    mappings: dict | None = None,
    figures: dict | None = None,
    ranges: dict | None = None,
    versions: dict[str, str] | None = None,
    builds: dict[str, dict[str, str]] | None = None,
) -> Dump:
    dump = Dump(
        model="m",
        order=["E"],
        errors={},
        results=results or {},
        datasets=datasets or {},
        mappings=mappings or {"E": []},
        figures=figures or {"E": {}},
        ranges=ranges or {"E": {}},
    )
    if versions is not None:
        dump["versions"] = versions
    if builds is not None:
        dump["builds"] = builds
    return dump


def written(tmp_path: Path, name: str, dump: Dump) -> str:
    path = tmp_path / name
    path.write_bytes(pickle.dumps(dump))
    return str(path)


def old_scan() -> Dump:
    """A scan of sbmlsim 0.8: (_point, dim_k), the state after the change at 2."""
    values = np.array([[1.0, 10.0], [2.0, 20.0], [3.5, 35.0], [4.0, 40.0]])
    return synthetic(
        {
            ("E", "task_x", "time"): (
                ("_point", "dim_k"),
                np.stack([T_OLD, T_OLD], axis=1),
            ),
            ("E", "task_x", "[A]"): (("_point", "dim_k"), values),
        }
    )


def new_scan(before: float = 3.0) -> Dump:
    """The scan of sbmlsim 0.9: (dim_k, time), the states before and after the change."""
    values = np.array([[1.0, 2.0, before, 3.5, 4.0], [10.0, 20.0, 30.0, 35.0, 40.0]])
    return synthetic(
        {
            ("E", "task_x", "time"): (("time",), T_NEW),
            ("E", "task_x", "[A]"): (("dim_k", "time"), values),
        },
        versions={"sbmlsim": "0.9.3", "pkpdsim": "0.5.0"},
    )


def test_the_matched_rows() -> None:
    match = matched_rows(T_OLD, T_NEW)
    assert match.old.tolist() == [0, 1, 2, 3]
    assert match.new.tolist() == [0, 1, 3, 4]
    assert (match.ignored_old, match.ignored_new, match.unmatched) == (0, 1, 0)
    same = matched_rows(T_NEW, T_NEW)
    assert same.new.tolist() == [0, 1, 2, 3, 4]
    assert (same.ignored_old, same.ignored_new) == (0, 0)


def test_a_dump_of_sbmlsim_08_compares_with_one_of_09(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    old, new = (
        written(tmp_path, "old.pkl", old_scan()),
        written(tmp_path, "new.pkl", new_scan()),
    )
    code = verify(["compare", old, new, "--threshold", "1e-12"])
    out = capsys.readouterr().out
    assert code == 0, out
    assert (
        "rows at a change: 1 rows of the new dump and 0 of the old dump left out" in out
    )
    assert "changes: E 1 rows" in out
    assert "new: pkpdsim 0.5.0, sbmlsim 0.9.3" in out


def test_two_dumps_of_sbmlsim_09_compare_every_row(tmp_path: Path) -> None:
    a = written(tmp_path, "a.pkl", new_scan())
    b = written(tmp_path, "b.pkl", new_scan(before=3.1))
    assert verify(["compare", a, a]) == 0
    assert verify(["compare", a, b, "--threshold", "1e-12"]) == 1


def ragged(times: Sequence[Sequence[float]], values: Sequence[Sequence[float]]) -> Dump:
    """A scan of sbmlsim 0.9 in the ragged layout: every point has its own times.

    The rows of a point are padded with NaN to the longest point, `(dim_k, _point)`.
    """
    rows = max(len(t) for t in times)

    def padded(columns: Sequence[Sequence[float]]) -> np.ndarray:
        return np.array([[*c, *[np.nan] * (rows - len(c))] for c in columns])

    return synthetic(
        {
            ("E", "task_x", "time"): (("dim_k", "_point"), padded(times)),
            ("E", "task_x", "[A]"): (("dim_k", "_point"), padded(values)),
        },
        versions={"sbmlsim": "0.9.3", "pkpdsim": "0.5.0"},
    )


def test_the_rows_of_every_point_are_matched_by_its_own_times(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """A ragged result: the second point has a change at 2, the first has not."""
    old = ragged([[0, 1, 2, 3], [0, 1, 2, 3]], [[1, 2, 3, 4], [10, 20, 35, 40]])
    new = ragged([[0, 1, 2, 3], [0, 1, 2, 2, 3]], [[1, 2, 3, 4], [10, 20, 30, 35, 40]])
    code = verify(
        [
            "compare",
            written(tmp_path, "old.pkl", old),
            written(tmp_path, "new.pkl", new),
            "--threshold",
            "1e-12",
        ]
    )
    out = capsys.readouterr().out
    assert code == 0, out
    assert "changes: E 1 rows" in out
    # the value after the change of the second point is compared
    new["results"]["E", "task_x", "[A]"][1][1, 3] = 36.0
    code = verify(
        [
            "compare",
            written(tmp_path, "old.pkl", old),
            written(tmp_path, "new.pkl", new),
            "--threshold",
            "1e-12",
        ]
    )
    assert code == 1


def test_a_time_of_another_point_than_the_first_is_compared(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    old = ragged([[0, 1, 2, 3], [0, 1, 2, 3]], [[1, 2, 3, 4], [10, 20, 30, 40]])
    new = ragged([[0, 1, 2, 3], [0, 1, 2, 5]], [[1, 2, 3, 4], [10, 20, 30, 40]])
    code = verify(
        [
            "compare",
            written(tmp_path, "old.pkl", old),
            written(tmp_path, "new.pkl", new),
        ]
    )
    out = capsys.readouterr().out
    assert code == 1, out
    # the old row at 3 is compared with the new rows at 2 and 5 interpolated,
    # the new row at 5 has no time of the old dump around it
    assert "rows at other times: 1 rows of the old dump and 0 of the new" in out
    assert "1 rows of the new dump without a time in the old dump" in out
    assert "times of E/task_x: 1 rows only in new, at 5" in out


def test_an_old_row_without_a_new_time_is_a_difference(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    old = old_scan()
    times = T_OLD.copy()
    # after the last time of the new dump: no interpolation
    times[3] = 3.5
    old["results"]["E", "task_x", "time"] = (
        ("_point", "dim_k"),
        np.stack([times, times], axis=1),
    )
    code = verify(
        [
            "compare",
            written(tmp_path, "old.pkl", old),
            written(tmp_path, "new.pkl", new_scan()),
        ]
    )
    assert code == 1
    assert (
        "1 rows of the old dump without a time in the new dump"
        in capsys.readouterr().out
    )


def timecourse(
    times: Sequence[float] | np.ndarray, values: Sequence[float] | np.ndarray
) -> Dump:
    """A dump of one task with its output times and the values of A."""
    return synthetic(
        {
            ("E", "t", "time"): (("time",), np.array(times)),
            ("E", "t", "[A]"): (("time",), np.array(values)),
        }
    )


def test_repeated_times_are_matched_from_the_last_row() -> None:
    # two dumps of sbmlsim 0.9 with other output times: every row of the
    # common times is compared, the states before a change included
    match = matched_rows(np.array([0.0, 1, 2, 2, 3]), np.array([0.0, 1, 2, 2, 3, 4]))
    assert match.old.tolist() == [0, 1, 2, 3, 4]
    assert match.new.tolist() == [0, 1, 2, 3, 4]
    assert match[2:6] == (0, 0, 0, 1)
    assert match.times_old.tolist() == []
    assert match.times_new.tolist() == [4.0]
    # the old dump has the time twice: its state after the change is compared
    match = matched_rows(np.array([0.0, 1, 2, 2, 3]), T_OLD)
    assert match.old.tolist() == [0, 1, 3, 4]
    assert match.new.tolist() == [0, 1, 2, 3]
    assert match[2:6] == (1, 0, 0, 0)
    # a change at the last time
    match = matched_rows(T_OLD, np.array([0.0, 1, 2, 3, 3]))
    assert match.old.tolist() == [0, 1, 2, 3]
    assert match.new.tolist() == [0, 1, 2, 4]
    assert match[2:6] == (0, 1, 0, 0)


def test_a_state_before_a_change_differs_although_the_times_differ(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    old = written(tmp_path, "o.pkl", timecourse([0, 1, 2, 2, 3], [1, 2, 3, 3.5, 4]))
    times = [0, 1, 2, 2, 3, 4]
    same = written(tmp_path, "s.pkl", timecourse(times, [1, 2, 3, 3.5, 4, 5]))
    other = written(tmp_path, "n.pkl", timecourse(times, [1, 2, 99, 3.5, 4, 5]))
    # the new output time is expected, the states at 2 are compared
    renames = tmp_path / "renames.json"
    renames.write_text(json.dumps({"expected": {"results": ["E/t/time"]}}))
    args = ["--threshold", "1e-12", "--renames", str(renames)]
    assert verify(["compare", old, same, *args]) == 0
    assert verify(["compare", old, other, *args]) == 1


def test_a_new_output_time_is_a_difference(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    old = written(tmp_path, "o.pkl", timecourse([0, 1, 2, 3], [1, 2, 3, 4]))
    new = written(tmp_path, "n.pkl", timecourse([0, 1, 2, 3, 4], [1, 2, 3, 4, 5]))
    capsys.readouterr()
    assert verify(["compare", old, new, "--threshold", "1e-12"]) == 1
    out = capsys.readouterr().out
    assert "1 rows of the new dump without a time in the old dump (0 expected)" in out
    assert "  times of E/t: 1 rows only in new, at 4\n" in out
    # the pattern of the variable time of the task names the new output times
    renames = tmp_path / "renames.json"
    renames.write_text(json.dumps({"expected": {"results": ["E/t/time"]}}))
    assert verify(["compare", old, new, "--renames", str(renames)]) == 0
    out = capsys.readouterr().out
    assert "1 rows of the new dump without a time in the old dump (1 expected)" in out
    assert "  times of E/t: 1 rows only in new, at 4 (expected)" in out
    # a time only the old dump has is a difference, also with the pattern
    assert verify(["compare", new, old, "--renames", str(renames)]) == 1
    assert "  times of E/t: 1 rows only in old, at 4\n" in capsys.readouterr().out


def dosed(times: np.ndarray, rate: float = 0.5) -> np.ndarray:
    """A decay with a second dose at 2: a trajectory of an adaptive grid.

    `times` has the time of the dose twice, the state before and after it.
    """
    before = np.exp(-rate * times)
    after = (math.exp(-2 * rate) + 1.0) * np.exp(-rate * (times - 2.0))
    values = np.where(times < 2.0, before, after)
    values[int(np.flatnonzero(times == 2.0)[0])] = math.exp(-2 * rate)
    return values


#: two adaptive grids of one simulation: other steps, the dose at 2 in both;
#: their interpolants are within 4e-7 of the trajectory (the steps of an
#: integrator at a relative tolerance of about 1e-8)
GRID_A = np.concatenate([np.linspace(0, 2, 61), np.linspace(2, 10, 121)])
GRID_B = np.concatenate([[0.0], np.geomspace(0.01, 2, 50), np.geomspace(2, 10, 100)])


def test_two_grids_of_one_trajectory_do_not_differ(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """A row at a time only one dump has is compared with the other interpolated."""
    old = written(tmp_path, "o.pkl", timecourse(GRID_A, dosed(GRID_A)))
    new = written(tmp_path, "n.pkl", timecourse(GRID_B, dosed(GRID_B)))
    capsys.readouterr()
    assert verify(["compare", old, new]) == 0, capsys.readouterr().out
    assert verify(["compare", new, old]) == 0
    out = capsys.readouterr().out
    assert "0 rows of the old dump without a time in the new dump" in out
    # the common times are 0, the dose at 2 (twice) and 10
    assert (
        f"rows at other times: {GRID_B.size - 4} rows of the old dump and "
        f"{GRID_A.size - 4} of the new dump compared with the other dump "
        "interpolated at their times" in out
    )


def test_a_trajectory_on_another_grid_differs(tmp_path: Path) -> None:
    """The interpolation leaves a difference of the trajectory: a slower decay."""
    old = written(tmp_path, "o.pkl", timecourse(GRID_A, dosed(GRID_A)))
    new = written(tmp_path, "n.pkl", timecourse(GRID_B, dosed(GRID_B, rate=0.501)))
    assert verify(["compare", old, new]) == 1
    assert verify(["compare", new, old]) == 1


def test_a_difference_between_the_common_times_is_found(tmp_path: Path) -> None:
    """A change of the trajectory of 1e-4 which vanishes at the common times."""
    old = written(tmp_path, "o.pkl", timecourse(GRID_A, dosed(GRID_A)))
    for eps, code in ((1e-4, 1), (1e-7, 0)):
        bump = np.where(
            GRID_B <= 2,
            1 + eps * np.sin(np.pi * GRID_B / 2),
            1 + eps * np.sin(np.pi * (GRID_B - 2) / 8),
        )
        new = written(tmp_path, "n.pkl", timecourse(GRID_B, dosed(GRID_B) * bump))
        assert verify(["compare", old, new]) == code


def tasks(values: dict[str, Sequence[float]]) -> Dump:
    """A dump of tasks of one experiment with the output times 0 to 3."""
    results: dict = {}
    for task, task_values in values.items():
        results["E", task, "time"] = (("time",), np.array([0.0, 1.0, 2.0, 3.0]))
        results["E", task, "[A]"] = (("time",), np.array(task_values))
    return synthetic(results)


def test_the_noise_of_a_control_is_no_difference(tmp_path: Path) -> None:
    """A variable which is zero in a task is compared at the scale of the experiment."""
    dose = [1.0, 0.5, 0.25, 0.125]
    old = written(tmp_path, "o.pkl", tasks({"dose": dose, "control": [0, 1e-24, 0, 0]}))
    noise = tasks({"dose": dose, "control": [0, 3e-24, -1e-25, 0]})
    assert verify(["compare", old, written(tmp_path, "n.pkl", noise)]) == 0
    # a small value which is not noise is compared at its own scale
    small = written(
        tmp_path, "s.pkl", tasks({"dose": dose, "control": [0, 1e-6, 0, 0]})
    )
    other = tasks({"dose": dose, "control": [0, 1.1e-6, 0, 0]})
    assert verify(["compare", small, written(tmp_path, "m.pkl", other)]) == 1


def test_the_interpolation_does_not_cross_a_change() -> None:
    """A row is interpolated between the rows of the segment of its time."""
    old_times = np.array([0.0, 1.0, 2.0, 2.0, 3.0])
    new_times = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.0, 2.5, 3.0, 4.0])
    match = matched_rows(old_times, new_times)
    assert match.old.tolist() == [0, 1, 2, 3, 4]
    assert match.interpolated_old.tolist() == []
    # the new rows between the common times are interpolated in the old dump,
    # the row after its last time is not
    assert match.interpolated_new.tolist() == [1, 3, 6]
    assert match.unmatched_new == 1
    assert match.times_new.tolist() == [4.0]
    # a repeated time which the other dump does not have is not interpolated
    match = matched_rows(np.array([0.0, 1.0, 1.5, 1.5, 3.0]), np.array([0.0, 1, 3]))
    assert match.interpolated_old.tolist() == []
    assert match.unmatched == 2
    # nor a row after the last time of the other dump
    match = matched_rows(np.array([0.0, 1.0, 1.0, 2.0]), np.array([0.0, 0.5, 1.0]))
    assert match.interpolated_old.tolist() == []
    assert match.interpolated_new.tolist() == [1]
    assert match.unmatched == 1


def test_close_times_are_not_matched() -> None:
    """Two rows of times 1e-7 apart are interpolated, not compared as one time."""
    match = matched_rows(np.array([0.0, 1.0, 2.0]), np.array([0.0, 1.0 + 1e-7, 2.0]))
    assert match.old.tolist() == [0, 2]
    assert match.interpolated_old.tolist() == [1]
    assert match.interpolated_new.tolist() == [1]


def test_repeated_times_of_the_old_dump_are_left_out(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    old = written(tmp_path, "o.pkl", timecourse([0, 1, 2, 2, 3], [1, 2, 9, 3, 4]))
    new = written(tmp_path, "n.pkl", timecourse([0, 1, 2, 3], [1, 2, 3, 4]))
    capsys.readouterr()
    assert verify(["compare", old, new, "--threshold", "1e-12"]) == 0
    assert (
        "0 rows of the new dump and 1 of the old dump left out"
        in capsys.readouterr().out
    )


def test_the_state_after_a_change_at_the_last_time_is_compared(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # sbmlsim 0.8 gave the state after a change at the last output time
    old = written(tmp_path, "o.pkl", timecourse([0, 1, 2, 3], [1, 2, 3, 10]))
    times = [0, 1, 2, 3, 3]
    new = written(tmp_path, "n.pkl", timecourse(times, [1, 2, 3, 4, 10]))
    swapped = written(tmp_path, "w.pkl", timecourse(times, [1, 2, 3, 10, 4]))
    capsys.readouterr()
    assert verify(["compare", old, new, "--threshold", "1e-12"]) == 0
    assert "changes: E 1 rows" in capsys.readouterr().out
    assert verify(["compare", old, swapped, "--threshold", "1e-12"]) == 1


def test_a_scan_of_one_point_compares_with_a_simulation(tmp_path: Path) -> None:
    values = np.array([1.0, 2.0, 3.0, 4.0])
    old = synthetic(
        {
            ("E", "task_x", "time"): (("_point", "dim_ref"), T_OLD[:, np.newaxis]),
            ("E", "task_x", "[A]"): (("_point", "dim_ref"), values[:, np.newaxis]),
        }
    )
    new = synthetic(
        {
            ("E", "task_x", "time"): (("time",), T_OLD),
            ("E", "task_x", "[A]"): (("time",), values),
        }
    )
    assert (
        verify(
            [
                "compare",
                written(tmp_path, "o.pkl", old),
                written(tmp_path, "n.pkl", new),
            ]
        )
        == 0
    )


def test_an_array_without_a_time_is_compared_relative_to_its_scale(
    tmp_path: Path,
) -> None:
    # the values at the end of the points of a scan, e.g. a clearance per renal
    # function, have a dimension of the scan and no time
    def dump(values: list[float]) -> Dump:
        return synthetic(
            {
                ("E", "task_x", "time"): (("time",), T_OLD),
                ("E", "task_x", "cl"): (("dim_renal",), np.array(values)),
            }
        )

    old = written(tmp_path, "o.pkl", dump([2.0, 4.0, 10.0]))
    close = written(tmp_path, "c.pkl", dump([2.0 + 1e-10, 4.0, 10.0]))
    far = written(tmp_path, "f.pkl", dump([2.1, 4.0, 10.0]))
    nan = written(tmp_path, "n.pkl", dump([np.nan, 4.0, 10.0]))
    assert verify(["compare", old, close]) == 0
    assert verify(["compare", old, far]) == 1
    assert verify(["compare", old, nan]) == 1


def test_arrays_without_a_time_of_other_dimensions_differ(tmp_path: Path) -> None:
    def dump(dim: str) -> Dump:
        return synthetic(
            {
                ("E", "task_x", "time"): (("time",), T_OLD),
                ("E", "task_x", "cl"): ((dim,), np.array([2.0, 4.0])),
            }
        )

    old = written(tmp_path, "o.pkl", dump("dim_renal"))
    assert verify(["compare", old, written(tmp_path, "n.pkl", dump("dim_q"))]) == 1


def test_arrays_of_other_dimensions_differ(tmp_path: Path) -> None:
    new = new_scan()
    _, values = new["results"]["E", "task_x", "[A]"]
    new["results"]["E", "task_x", "[A]"] = (("dim_q", "time"), values)
    old = written(tmp_path, "o.pkl", old_scan())
    assert verify(["compare", old, written(tmp_path, "n.pkl", new)]) == 1


def test_arrays_of_other_sizes_differ(tmp_path: Path) -> None:
    new = new_scan()
    _, values = new["results"]["E", "task_x", "[A]"]
    new["results"]["E", "task_x", "[A]"] = (
        ("dim_k", "time"),
        np.vstack([values, values]),
    )
    old = written(tmp_path, "o.pkl", old_scan())
    assert verify(["compare", old, written(tmp_path, "n.pkl", new)]) == 1


def test_the_ids_of_the_old_dump_become_sids(tmp_path: Path) -> None:
    def figure(sid: str, task: str) -> str:
        spec = {"sid": sid, "curves": [{"task_id": task, "name": "po 2.5"}]}
        return json.dumps(spec, sort_keys=True)

    def dump(task: str, dataset: str, fig: str, mapping: str) -> Dump:
        return synthetic(
            {
                ("E", task, "time"): (("time",), T_OLD),
                ("E", task, "[A]"): (("time",), T_OLD),
            },
            datasets={("E", dataset, "mean"): np.array([1.0])},
            figures={"E": {fig: figure(fig, task)}},
            ranges={"E": {mapping: (0.0, 1.0, 0.0, 1.0, "hr")}},
        )

    old = written(tmp_path, "o.pkl", dump("task_po 2.5", "Fig1_A b", "fig 1", "fm_a b"))
    new = written(tmp_path, "n.pkl", dump("task_po_2p5", "Fig1_A_b", "fig_1", "fm_a_b"))
    renames = tmp_path / "renames.json"
    renames.write_text(json.dumps({"sids": True}))
    assert verify(["compare", old, new]) == 1
    assert verify(["compare", old, new, "--renames", str(renames)]) == 0


def test_the_ids_of_a_mapping_description_are_renamed() -> None:
    description = (
        "fm_a b|Fig1_A b|None|time|mean|4|mean_sd|None"
        "|None|task_po 2.5|time|A|None|None|None|weight=None|[('route', 'IV')]"
    )
    renamed = apply_renames(synthetic(mappings={"E": [description]}), {"sids": True})
    assert renamed["mappings"]["E"] == [
        "fm_a_b|Fig1_A_b|None|time|mean|4|mean_sd|None"
        "|None|task_po_2p5|time|A|None|None|None|weight=None|[('route', 'IV')]"
    ]


def test_a_key_with_a_bar_cannot_be_renamed() -> None:
    description = (
        "fm|a|x|Fig1|None|time|mean|4|mean_sd|None"
        "|None|task|time|A|None|None|None|weight=None|[]"
    )
    with pytest.raises(ValueError, match="cannot be renamed"):
        apply_renames(synthetic(mappings={"E": [description]}), {"sids": True})


def test_two_ids_with_one_new_id_are_an_error() -> None:
    dump = synthetic(
        datasets={
            ("E", "a b", "x"): np.array([1.0]),
            ("E", "a-b", "x"): np.array([2.0]),
        }
    )
    with pytest.raises(ValueError, match="same new id"):
        apply_renames(dump, {"sids": True})


@pytest.mark.parametrize(
    "parts",
    [
        {
            "results": {
                ("E", "task a b", "x"): (("time",), T_OLD),
                ("E", "task a-b", "x"): (("time",), T_OLD),
            }
        },
        {"figures": {"E": {"fig a": '{"sid": "fig a"}', "fig-a": '{"sid": "fig-a"}'}}},
        {"ranges": {"E": {"fm a": "no x values", "fm-a": "no x values"}}},
    ],
    ids=["results", "figures", "ranges"],
)
def test_two_ids_of_a_part_with_one_new_id_are_an_error(parts: dict) -> None:
    with pytest.raises(ValueError, match="same new id"):
        apply_renames(synthetic(**parts), {"sids": True})


def test_two_keys_with_one_key_after_the_strings_are_an_error() -> None:
    dump = synthetic(
        datasets={("E", "d", "x1"): np.array([1.0]), ("E", "d", "x2"): np.array([2.0])}
    )
    with pytest.raises(ValueError, match="same new key"):
        apply_renames(dump, {"strings": {"x1": "x2"}})


def test_the_strings_do_not_rename_the_parts_of_a_dump() -> None:
    dump = synthetic(datasets={("E", "a", "x"): np.array([1.0])})
    renamed = apply_renames(dump, {"strings": {"a": "b"}})
    assert renamed.keys() == dump.keys()
    assert list(renamed["datasets"]) == [("E", "b", "x")]


def test_the_renames_do_not_change_the_old_dump() -> None:
    descriptions = ["b", "a"]
    dump = synthetic(mappings={"E": descriptions})
    before = copy.deepcopy(dump)
    result = apply_renames(dump, {})
    assert result["mappings"]["E"] == ["a", "b"]
    assert dump["mappings"]["E"] is descriptions
    assert descriptions == ["b", "a"]
    assert dump.keys() == before.keys()


#: the figures of `build_figures` as a dump of pkpdsim 0.4.1 with sbmlsim 0.8.5
#: stored them (written by that function with sbmlsim 0.8.5, which sets a
#: setting as an attribute of the figure, `setattr`)
FIGURES_08: dict[str, str] = json.loads(
    (TESTS_DIR / "data" / "figures-sbmlsim-0.8.5.json").read_text(encoding="utf-8")
)["figures"]


def setting_of_09(figure: Figure, name: str, value: object) -> None:
    """Set a setting of a figure as sbmlsim 0.9 does (R1.5 of the port)."""
    figure.settings = {**figure.settings, name: value}


def figures_09(
    change: Callable[[dict[str, Figure]], object] | None = None,
) -> dict[str, str]:
    """Get the figures of `build_figures` as a dump of sbmlsim 0.9 stores them."""
    experiment = FigureStudy(base_path=TOY_PATH, data_path=DATA_PATHS)
    figures = build_figures(experiment, setting_of_09)
    if change is not None:
        change(figures)
    return {sid: figure_spec(figure, experiment) for sid, figure in figures.items()}


def compare_figures(
    tmp_path: Path, old: dict[str, str], new: dict[str, str], *args: str
) -> int:
    return verify(
        [
            "compare",
            written(tmp_path, "o.pkl", synthetic(figures={"E": old})),
            written(tmp_path, "n.pkl", synthetic(figures={"E": new})),
            *args,
        ]
    )


def test_the_settings_of_a_figure_of_sbmlsim_08_are_those_of_sbmlsim_09() -> None:
    from sbmlsim.plot import FigureSettings as SbmlsimFigureSettings

    names = tuple(f.name for f in dataclasses.fields(SbmlsimFigureSettings))
    # the settings only newer sbmlsim has (absent from 0.9.5) are not those of 0.8
    old = tuple(name for name in names if name not in NEW_FIGURE_SETTINGS)
    assert old == FIGURE_SETTINGS_08
    assert set(names) - set(old) <= set(NEW_FIGURE_SETTINGS)


def test_the_figures_of_sbmlsim_08_compare_with_those_of_09(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    assert compare_figures(tmp_path, FIGURES_08, figures_09()) == 0
    out = capsys.readouterr().out
    assert "figures: 3 old, 3 new, 0 differ" in out
    assert "figures across versions of sbmlsim: 3" in out


@pytest.mark.parametrize(
    ("change", "difference"),
    [
        (
            lambda f: setting_of_09(f["Fig1"], "legend_position", "inside"),
            "_settings.legend_position 'outside' -> 'inside'",
        ),
        (
            lambda f: setattr(f["Fig3"], "settings", None),
            "_settings.fig_dpi only in old: 300",
        ),
        (
            lambda f: setting_of_09(f["Fig2"], "legend_fontsize", 9),
            "_settings.legend_fontsize only in new: 9",
        ),
        (lambda f: setattr(f["Fig1"], "width", 12), "width 10 -> 12"),
        (lambda f: setattr(f["Fig2"], "height", None), "height 15 -> 10.0"),
        (
            lambda f: setting_of_09(f["Fig3"], "panel_width", 6.0),
            "width 14.0 -> 12.0",
        ),
        (
            lambda f: setattr(f["Fig3"].subplots[0].plot, "width", 9.0),
            "subplots[0].plot.width None -> 9.0",
        ),
    ],
)
def test_a_setting_or_a_size_across_versions_is_compared(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    change: Callable[[dict[str, Figure]], object],
    difference: str,
) -> None:
    assert compare_figures(tmp_path, FIGURES_08, figures_09(change)) == 1
    out = capsys.readouterr().out
    assert "figures: 3 old, 3 new, 1 differ (0 expected)" in out
    assert difference in out


def test_a_figure_of_one_dump_only_is_a_difference(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    new = figures_09()
    del new["Fig3"]
    assert compare_figures(tmp_path, FIGURES_08, new) == 1
    assert "figures: 3 old, 2 new, 1 differ (0 expected)" in capsys.readouterr().out
    old = dict(FIGURES_08)
    del old["Fig1"]
    assert compare_figures(tmp_path, old, figures_09()) == 1


def with_new_fields(spec: str) -> str:
    """Add a selection, a dimension of the lines and a band to a figure."""
    figure = json.loads(spec)
    plot = figure["subplots"][0]["plot"]
    curve = plot["_curves"][0]
    curve["over"] = ["dim_dose"]
    curve["y"]["sel"] = {"dim_bw": 0}
    plot["_bands"] = [{"__type__": "Band", "across": ["dim_sens"]}]
    return json.dumps(figure, sort_keys=True)


def test_the_fields_only_sbmlsim_09_has_are_listed(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    new = figures_09()
    new["Fig2"] = with_new_fields(new["Fig2"])
    assert compare_figures(tmp_path, FIGURES_08, new) == 0
    out = capsys.readouterr().out
    assert 'new fields: E/Fig2: subplots[0].plot._curves[0].over = ["dim_dose"]' in out
    assert (
        'new fields: E/Fig2: subplots[0].plot._curves[0].y.sel = {"dim_bw": 0}' in out
    )
    assert "new fields: E/Fig2: subplots[0].plot._bands = [{" in out
    # an empty selection is the default and not listed
    assert "E/Fig1" not in out


def test_two_dumps_of_sbmlsim_09_compare_the_new_fields(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    old = figures_09()
    new = dict(old)
    new["Fig2"] = with_new_fields(new["Fig2"])
    assert compare_figures(tmp_path, old, new) == 1
    assert "figures across versions" not in capsys.readouterr().out


def test_the_settings_and_the_size_a_figure_is_rendered_with() -> None:
    class Wide(FigureStudy):
        figure_settings = MappingProxyType(
            {"panel_width": 3.0, "legend_position": "outside"}
        )

    experiment = Wide(base_path=TOY_PATH, data_path=DATA_PATHS)
    figure = Figure(
        experiment, sid="f", num_cols=2, settings={"legend_position": "inside"}
    )
    resolved = json.loads(figure_spec(figure, experiment))["__resolved__"]
    assert resolved == {
        "settings": {"panel_width": 3.0, "legend_position": "inside"},
        "width": 6.0,
        "height": 5.0,
    }


def test_the_figure_fields_are_left_out_at_the_level_of_the_figure(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    def change(figures: dict[str, Figure]) -> None:
        figures["Fig2"].subplots[0].plot.width = 14.0

    renames = tmp_path / "renames.json"
    renames.write_text(json.dumps({"figure_fields": ["width", "height"]}))
    old, new = figures_09(), figures_09(change)
    assert compare_figures(tmp_path, old, new, "--renames", str(renames)) == 1
    assert "subplots[0].plot.width 8.333333333333334 -> 14.0" in capsys.readouterr().out
    # a field of the figure is left out
    extra = {**json.loads(old["Fig1"]), "extra": 1}
    new = {**old, "Fig1": json.dumps(extra, sort_keys=True)}
    renames.write_text(json.dumps({"figure_fields": ["extra"]}))
    assert compare_figures(tmp_path, old, new) == 1
    assert compare_figures(tmp_path, old, new, "--renames", str(renames)) == 0


def test_a_dump_stores_the_settings_and_the_size_of_a_figure(tmp_path: Path) -> None:
    out = tmp_path / "dump.pkl"
    args = ["dump", "--root", str(TESTS_DIR), "--model", "pkpdsim_toy"]
    assert verify([*args, "--out", str(out), "ToyDecay"]) == 0
    spec = json.loads(load(str(out))["figures"]["ToyDecay"]["fig_decay"])
    assert spec["__resolved__"] == {"settings": {}, "width": 7.0, "height": 5.0}


def test_a_figure_field_of_one_dump_only_is_named(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    old_spec = {"sid": "f", "plots": [{"x": 1}]}
    new_spec = {"sid": "f", "width": None, "plots": [{"x": 1, "_bands": []}]}
    old = written(
        tmp_path, "o.pkl", synthetic(figures={"E": {"f": json.dumps(old_spec)}})
    )
    new = written(
        tmp_path, "n.pkl", synthetic(figures={"E": {"f": json.dumps(new_spec)}})
    )
    capsys.readouterr()
    assert verify(["compare", old, new]) == 1
    out = capsys.readouterr().out
    # a field which is null in one dump and missing in the other is a difference
    assert "width only in new: None" in out
    assert "plots[0]._bands only in new: []" in out


def test_an_expected_difference_is_not_counted(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    old = written(
        tmp_path, "o.pkl", synthetic(datasets={("E", "d", "mean"): np.array([1.0])})
    )
    new = written(
        tmp_path, "n.pkl", synthetic(datasets={("E", "d", "mean"): np.array([2.0])})
    )
    renames = tmp_path / "renames.json"
    renames.write_text(json.dumps({"expected": {"datasets": ["E/d/*"]}}))
    capsys.readouterr()
    assert verify(["compare", old, new, "--renames", str(renames)]) == 0
    assert "datasets: 1 old, 1 new, 1 differ (1 expected)" in capsys.readouterr().out


def mapping(key: str, y: str, task: str) -> str:
    """The description of a fit mapping of the column y to [A] of a task."""
    return (
        f"{key}|d|None|time|{y}|4|None|None|None|{task}|time|[A]|None|None|None"
        "|weight=None|[]"
    )


def test_an_expected_mapping_difference_is_one_mapping(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    old = synthetic(
        mappings={"E": [mapping("fm1", "A", "t"), mapping("fm2", "A", "t")]}
    )
    new = synthetic(
        mappings={"E": [mapping("fm1", "B", "t"), mapping("fm2", "A", "u")]}
    )
    a, b = written(tmp_path, "o.pkl", old), written(tmp_path, "n.pkl", new)
    renames = tmp_path / "renames.json"
    capsys.readouterr()
    for patterns, code in ((["E"], 1), (["E/fm1"], 1), (["E/fm*"], 0)):
        renames.write_text(json.dumps({"expected": {"mappings": patterns}}))
        assert verify(["compare", a, b, "--renames", str(renames)]) == code, patterns
    out = capsys.readouterr().out
    assert "mappings: 2 old, 2 new, 2 differ (1 expected)" in out
    assert "expected difference which does not occur: mappings: E" in out
    assert (
        f"  E/fm2\n    - {mapping('fm2', 'A', 't')}\n    + {mapping('fm2', 'A', 'u')}"
        in out
    )


def test_two_mappings_of_one_key_are_an_error(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    dump = synthetic(mappings={"E": [mapping("fm", "A", "t"), mapping("fm", "B", "t")]})
    path = written(tmp_path, "d.pkl", dump)
    assert verify(["compare", path, path]) == 2
    assert "Two fit mappings of E have the key fm" in capsys.readouterr().err


def test_a_bracket_of_a_pattern_is_literal(tmp_path: Path) -> None:
    def dump(a: float, bracket: float) -> Dump:
        return synthetic(
            {
                ("E", "t", "time"): (("time",), T_OLD),
                ("E", "t", "A"): (("time",), np.full(4, a)),
                ("E", "t", "[A]"): (("time",), np.full(4, bracket)),
            }
        )

    old = written(tmp_path, "o.pkl", dump(1.0, 1.0))
    renames = tmp_path / "renames.json"
    renames.write_text(json.dumps({"expected": {"results": ["E/t/[A]"]}}))
    other_bracket = written(tmp_path, "b.pkl", dump(1.0, 2.0))
    other_plain = written(tmp_path, "p.pkl", dump(2.0, 1.0))
    assert verify(["compare", old, other_bracket, "--renames", str(renames)]) == 0
    assert verify(["compare", old, other_plain, "--renames", str(renames)]) == 1


def test_an_expected_difference_which_does_not_occur_is_counted(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    same = written(
        tmp_path, "s.pkl", synthetic(datasets={("E", "d", "mean"): np.array([1.0])})
    )
    renames = tmp_path / "renames.json"
    renames.write_text(json.dumps({"expected": {"figures": ["E/fig_x"]}}))
    capsys.readouterr()
    assert verify(["compare", same, same, "--renames", str(renames)]) == 1
    assert (
        "expected difference which does not occur: figures: E/fig_x"
        in capsys.readouterr().out
    )


@pytest.mark.parametrize(
    "renames",
    [
        {"ids": {"a": 1}},
        {"sids": "yes"},
        {"figure_fields": "width"},
        {"expected": {"results": "E/*"}},
        {"expected": {"plots": ["E"]}},
    ],
)
def test_a_renames_file_with_wrong_new_parts_is_an_error(
    tmp_path: Path, renames: dict
) -> None:
    path = tmp_path / "renames.json"
    path.write_text(json.dumps(renames))
    with pytest.raises(ValueError):
        load_renames(str(path))


def test_png_pairs_the_figures_of_renamed_keys(tmp_path: Path) -> None:
    old, new = tmp_path / "old", tmp_path / "new"
    write_png(old / "A" / "A_fig 2.5.png")
    write_png(new / "A" / "A_fig_2p5.png")
    assert verify(["png", str(old), str(new)]) == 1
    assert verify(["png", "--sids", str(old), str(new)]) == 0


def test_png_with_two_old_figures_of_one_sid_fails(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    old, new = tmp_path / "old", tmp_path / "new"
    write_png(old / "A" / "A_fig 2.5.png")
    write_png(old / "A" / "A_fig-2.5.png")
    write_png(new / "A" / "A_fig_2p5.png")
    assert verify(["png", "--sids", str(old), str(new)]) == 2
    message = f"two old figures have the same SId: {Path('A', 'A_fig_2p5.png')}"
    assert message in capsys.readouterr().err


def test_a_dump_names_the_versions(tmp_path: Path) -> None:
    import pkpdsim

    out = tmp_path / "dump.pkl"
    assert (
        verify(
            [
                "dump",
                "--root",
                str(TESTS_DIR),
                "--model",
                "pkpdsim_toy",
                "--out",
                str(out),
                "ToyDecay",
            ]
        )
        == 0
    )
    versions = load(str(out))["versions"]
    assert versions["sbmlsim"] == importlib.metadata.version("sbmlsim")
    assert versions["pkpdsim"] == pkpdsim.__version__


def _curve(name: str, **kwargs: object) -> str:
    """Get a figure of one curve as a dump stores it."""
    curve = {"__type__": "Curve", "name": name, "kwargs": kwargs}
    return json.dumps(
        {"__type__": "Figure", "subplots": [{"plot": {"_curves": [curve]}}]}
    )


def test_the_label_of_a_curve_across_versions_of_sbmlsim_is_its_name() -> None:
    """sbmlsim 0.14 keeps the label of a curve as its name and writes ` ± SD`."""
    old = synthetic(figures={"E": {"Fig1": _curve("d±SD", label="d±SD", color="k")}})
    new = synthetic(figures={"E": {"Fig1": _curve("d ± SD", color="k")}})
    aligned = align_figures(old, new)
    assert aligned.old["figures"]["E"]["Fig1"] == aligned.new["figures"]["E"]["Fig1"]
    # a label which is not the name of its curve stays a difference
    other = synthetic(figures={"E": {"Fig1": _curve("d", label="x", color="k")}})
    aligned = align_figures(other, new)
    assert aligned.old["figures"]["E"]["Fig1"] != aligned.new["figures"]["E"]["Fig1"]


# --- the builds per sex -----------------------------------------------------------


def test_a_dump_records_the_build_of_every_task(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """A study which moves to the build of another sex is a difference."""
    old = tmp_path / "old.pkl"
    common = ["--root", str(TESTS_DIR), "--model", "pkpdsim_sexes"]
    assert verify(["dump", *common, "--out", str(old)]) == 0
    dumped = load(str(old))
    assert dumped["builds"] == {
        "SexesMen": {"task_decay": "male"},
        "SexesWomen": {"task_decay": "female"},
        "SexesGroups": {"task_men": "male", "task_women": "female"},
    }
    capsys.readouterr()
    assert verify(["compare", str(old), str(old)]) == 0
    assert "builds: 4 old, 4 new, 0 differ" in capsys.readouterr().out

    moved = copy.deepcopy(dumped)
    moved["builds"]["SexesWomen"]["task_decay"] = "male"
    new = written(tmp_path, "new.pkl", moved)
    assert verify(["compare", str(old), new]) == 1
    out = capsys.readouterr().out
    assert "builds: 4 old, 4 new, 1 differ (0 expected)" in out
    assert "SexesWomen/task_decay" in out
    assert "'female' -> 'male'" in out
    renames = tmp_path / "renames.json"
    renames.write_text(json.dumps({"expected": {"builds": ["SexesWomen/*"]}}))
    assert verify(["compare", str(old), new, "--renames", str(renames)]) == 0


def test_the_tasks_of_the_builds_are_renamed(tmp_path: Path) -> None:
    dump = synthetic(builds={"E": {"task_a": "female"}})
    renamed = apply_renames(dump, {"ids": {"task_a": "task_b"}})
    assert renamed["builds"] == {"E": {"task_b": "female"}}


def test_a_dump_without_builds_is_compared_without_them(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    old = written(tmp_path, "old.pkl", synthetic())
    new = written(tmp_path, "new.pkl", synthetic(builds={"E": {"task": "male"}}))
    assert verify(["compare", old, new]) == 0
    assert "builds: in the new dump only, not compared" in capsys.readouterr().out
