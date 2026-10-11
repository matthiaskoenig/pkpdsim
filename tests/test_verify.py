"""The verification of a change of a model, run on the toy model."""

import dataclasses
import json
import pickle
import re
import shutil
import sys
from collections.abc import Callable
from pathlib import Path

import numpy as np
import pytest
from PIL import Image
from sbmlsim.experiment import SimulationExperiment
from sbmlsim.plot import Figure

from pkpdsim.model import FigureSettings
from pkpdsim.verify import (
    apply_renames,
    describe,
    load,
    load_renames,
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
    assert settings["figure"]["legend_fontsize"] == 7
    assert settings["rc"]["savefig.dpi"] == 30
    assert settings["calls"] == {"simulator": 1, "runner": 1, "run": 1}
    assert "<root>" in settings["simulator"]["model"]


def test_capture_settings_sees_a_difference() -> None:
    from pkpdsim.verify import capture_settings

    other = dataclasses.replace(
        MODEL, figures=FigureSettings(formats=("png",), legend_fontsize=9)
    )
    a = capture_settings(MODEL.run_experiments, ToyDecay, TESTS_DIR)
    b = capture_settings(other.run_experiments, ToyDecay, TESTS_DIR)
    assert a["figure"] != b["figure"]
    assert a["run"]["figure_formats"] != b["run"]["figure_formats"]


def test_settings_leave_the_figure_settings_as_they_were(tmp_path: Path) -> None:
    before = Figure.legend_fontsize
    verify(
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
    assert Figure.legend_fontsize == before


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


def _fake_run(positional: bool) -> "Callable[..., object]":
    from sbmlsim.experiment import ExperimentRunner
    from sbmlsim.simulator.simulation_serial import SimulatorSerial

    def run(experiment: type[SimulationExperiment], output_dir: Path) -> None:
        if positional:
            simulator = SimulatorSerial(MODEL_PATH, absolute_tolerance=1e-10)
            runner = ExperimentRunner([experiment], TOY_PATH, DATA_PATHS, simulator)
            runner.run_experiments(output_dir, False, False, ["png"])
        else:
            simulator = SimulatorSerial(model=MODEL_PATH, absolute_tolerance=1e-10)
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
