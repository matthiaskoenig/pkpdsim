"""The runner of the sensitivity analyses of a script, on the toy model.

The toy model is the decay of A with the rate `k` in the compartment `Vc`; its
concentration does not depend on `Vc`. The group `large` sets `Vc`, so an
analysis over both groups varies only `k`.
"""

import dataclasses
import logging
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from sbmlsim import sensitivity
from sbmlsim.simulation import PK, Formula, Scan, Simulation, sampling

from pkpdsim.model import SimulatorSettings
from pkpdsim.run import create_simulator
from pkpdsim.sensitivity import (
    FastAnalysis,
    LocalAnalysis,
    SamplingAnalysis,
    SensitivityGroup,
    SensitivityOutput,
    SensitivityParameter,
    SensitivityRunner,
    SobolAnalysis,
    parameters_from_sbml,
    set_relative_bounds,
)
from pkpdsim_toy import MODEL, MODEL_PATH

#: the simulations of the tests run in the process of the test
SERIAL = SimulatorSettings(
    absolute_tolerance=1e-10, relative_tolerance=1e-10, n_workers=1
)

OBSERVABLES = [
    PK("pk_A", "[A]", parameters=["auc_inf_obs", "thalf"]),
    Formula("A_auc", "pk_A.auc_inf_obs", unit="mM*min"),
    Formula("A_thalf", "pk_A.thalf", unit="min"),
    Formula("A_min", "min([A])", unit="mM"),
]

#: the third output is named in the tables by another id than its observable
OUTPUTS = [
    SensitivityOutput("A_auc", "A AUC∞"),
    SensitivityOutput("A_thalf", "A half-life"),
    SensitivityOutput("[A]_min", "A minimum", observable="A_min"),
]


def simulation(changes: dict[str, float]) -> Simulation:
    return Simulation(end=60, steps=600, time_unit="min", preinit_changes=changes)


CONTROL = SensitivityGroup("control", "Control", simulation({}))
LARGE = SensitivityGroup("large", "Large volume", simulation({"Vc": 2.0}))


def toy_parameters() -> list[SensitivityParameter]:
    parameters = parameters_from_sbml(MODEL_PATH)
    set_relative_bounds(parameters, 0.15)
    return parameters


def runner(tmp_path: Path, **kwargs: object) -> SensitivityRunner:
    options: dict[str, Any] = {
        "model": dataclasses.replace(MODEL, results_path=tmp_path),
        "parameters": toy_parameters,
        "outputs": OUTPUTS,
        "observables": OBSERVABLES,
        "analyses": {
            "local": LocalAnalysis(groups=[CONTROL, LARGE]),
            "sampling": SamplingAnalysis(groups=[CONTROL, LARGE], n=8),
            "fast": FastAnalysis(groups=[CONTROL], n=65),
            "sobol": SobolAnalysis(groups=[CONTROL], n=8),
        },
    }
    options.update(kwargs)
    return SensitivityRunner(**options)


def test_an_output_is_its_observable_unless_it_names_another() -> None:
    assert SensitivityOutput("cmax", "Cmax").observable == "cmax"
    output = SensitivityOutput("[Cve_api]_cmax", "API Cmax", observable="Cve_api_cmax")
    assert output.observable == "Cve_api_cmax"


def test_the_defaults_of_the_analyses() -> None:
    assert LocalAnalysis(groups=[CONTROL]).delta == 0.01
    assert SamplingAnalysis(groups=[CONTROL]).n == 1000
    assert SamplingAnalysis(groups=[CONTROL]).quantiles == (
        0.025,
        0.05,
        0.5,
        0.95,
        0.975,
    )
    assert FastAnalysis(groups=[CONTROL]).n == 1000
    assert SobolAnalysis(groups=[CONTROL]).n == 4096
    # the groups are kept as a tuple
    assert LocalAnalysis(groups=[CONTROL, LARGE]).groups == (CONTROL, LARGE)


def test_local_writes_the_normalized_sensitivities_of_every_group(
    tmp_path: Path,
) -> None:
    path = runner(tmp_path).run(["local"], simulator=SERIAL)
    assert path == tmp_path / "sensitivity"
    assert (path / "parameters.tsv").is_file()
    assert (path / "parameters.tex").is_file()
    for group in ["control", "large"]:
        group_path = path / "local" / group
        assert sorted(p.name for p in group_path.iterdir()) == sorted(
            [
                "A_auc.normalized.tsv",
                "A_thalf.normalized.tsv",
                "[A]_min.normalized.tsv",
                "normalized.png",
            ]
        )
    # the table of an output is named by its id, its column too
    df = pd.read_csv(path / "local" / "control" / "[A]_min.normalized.tsv", sep="\t")
    assert "[A]_min.normalized" in df.columns


def test_local_is_the_analysis_of_sbmlsim(tmp_path: Path) -> None:
    """The table is the local analysis of sbmlsim on the scan of the group."""
    path = runner(tmp_path).run(["local"], simulator=SERIAL)
    # `Vc` is set by the group `large`, the local analysis varies `k` only
    design = sampling.local(
        ["k"], delta=0.01, model=MODEL_PATH, simulation=CONTROL.simulation
    )
    result = sensitivity.local(
        create_simulator(SERIAL).run(
            MODEL_PATH,
            Scan(CONTROL.simulation, [design]),
            OBSERVABLES,
            keep=["A_auc", "A_thalf", "A_min"],
        )
    )
    expected = result.to_dataframe("A_thalf.normalized")
    df = pd.read_csv(
        path / "local" / "control" / "A_thalf.normalized.tsv", sep="\t", index_col=0
    )
    assert list(df.index) == list(expected.index)
    assert np.allclose(df["A_thalf.normalized"], expected["A_thalf.normalized"])
    # the half-life is inversely proportional to the rate
    assert df["A_thalf.normalized"].iloc[0] == pytest.approx(-1.0, rel=1e-2)


def test_an_analysis_varies_no_parameter_which_a_group_sets(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    with caplog.at_level(logging.INFO, logger="pkpdsim.sensitivity"):
        path = runner(tmp_path).run(["local"], simulator=SERIAL)
    assert "Removing fixed parameter: Vc" in caplog.text
    df = pd.read_csv(path / "local" / "control" / "A_auc.normalized.tsv", sep="\t")
    assert "Vc" not in df.to_string()


def test_sampling_writes_the_samples_the_summary_and_the_distributions(
    tmp_path: Path,
) -> None:
    path = runner(tmp_path).run(["sampling"], simulator=SERIAL)
    for group in ["control", "large"]:
        group_path = path / "sampling" / group
        assert sorted(p.name for p in group_path.iterdir()) == [
            "distributions.png",
            "samples.tsv",
            "summary.tsv",
        ]
        samples = pd.read_csv(group_path / "samples.tsv", sep="\t")
        assert len(samples) == 8
        assert {"A_auc", "A_thalf", "[A]_min"} <= set(samples.columns)
        # the values of the design, `Vc` is set by the group `large`
        assert "k" in samples.columns
        assert "Vc" not in samples.columns
        summary = pd.read_csv(group_path / "summary.tsv", sep="\t")
        assert summary["output"].tolist() == ["A_auc", "A_thalf", "[A]_min"]
        assert summary["name"].tolist() == ["A AUC∞", "A half-life", "A minimum"]
        assert list(summary.columns) == [
            "output",
            "name",
            "mean",
            "sd",
            "cv",
            "min",
            "max",
            "q0.025",
            "q0.05",
            "q0.5",
            "q0.95",
            "q0.975",
            "unit",
        ]
        assert summary["unit"].tolist()[1] == "minute"


@pytest.mark.parametrize(
    ("name", "n"), [("fast", 65), ("sobol", 8)], ids=["fast", "sobol"]
)
def test_a_global_analysis_writes_s1_and_st_of_every_output(
    tmp_path: Path, name: str, n: int
) -> None:
    path = runner(tmp_path).run([name], simulator=SERIAL)
    group_path = path / name / "control"
    expected = []
    for uid in ["A_auc", "A_thalf", "[A]_min"]:
        expected += [f"{uid}.S1.tsv", f"{uid}.ST.tsv", f"{uid}.png"]
    assert sorted(p.name for p in group_path.iterdir()) == sorted(expected)
    st = pd.read_csv(group_path / "A_thalf.ST.tsv", sep="\t", index_col=0)
    assert "A_thalf.ST" in st.columns
    assert not (path / name / "large").exists()


def test_the_sample_size_of_a_smoke_run(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("PKPDSIM_SENSITIVITY_N", "4")
    analyses = {"sampling": SamplingAnalysis(groups=[CONTROL])}
    path = runner(tmp_path, analyses=analyses, default=["sampling"]).run(
        simulator=SERIAL
    )
    assert (
        len(pd.read_csv(path / "sampling" / "control" / "samples.tsv", sep="\t")) == 4
    )


def test_the_analyses_run_by_default_and_into_another_folder(
    tmp_path: Path,
) -> None:
    analyses = {
        "local": LocalAnalysis(groups=[CONTROL]),
        "sampling": SamplingAnalysis(groups=[CONTROL], n=4),
    }
    path = runner(tmp_path, analyses=analyses).run(
        path=tmp_path / "single", simulator=SERIAL
    )
    assert path == tmp_path / "single"
    assert (path / "local" / "control").is_dir()
    assert not (path / "sampling").exists()


def test_an_unknown_analysis_is_refused(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match=r"\['morris'\].*'local', 'sampling'"):
        runner(tmp_path).run(["local", "morris"], simulator=SERIAL)
    assert not (tmp_path / "sensitivity").exists()


def test_the_simulator_of_the_model_with_a_pool(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A run uses the tolerances of the model and lets sbmlsim choose the processes."""
    settings: list[SimulatorSettings] = []

    def serial(s: SimulatorSettings) -> object:
        settings.append(s)
        return create_simulator(dataclasses.replace(s, n_workers=1))

    monkeypatch.setattr("pkpdsim.sensitivity.create_simulator", serial)
    analyses = {"local": LocalAnalysis(groups=[CONTROL])}
    runner(tmp_path, analyses=analyses).run()
    assert settings == [dataclasses.replace(MODEL.simulator, n_workers=None)]


def test_main_runs_the_analyses_of_the_command_line(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: list[tuple[str, ...]] = []

    def run(self: SensitivityRunner, names: tuple[str, ...]) -> Path:
        calls.append(tuple(names))
        return tmp_path

    monkeypatch.setattr(SensitivityRunner, "run", run)
    toy = runner(tmp_path, default=("local", "sampling"))
    toy.main([])
    toy.main(["sobol", "fast"])
    assert calls == [("local", "sampling"), ("sobol", "fast")]


def test_main_refuses_an_unknown_analysis(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    with pytest.raises(SystemExit) as info:
        runner(tmp_path).main(["morris"])
    assert info.value.code == 2
    assert "invalid choice: 'morris'" in capsys.readouterr().err
    assert not (tmp_path / "sensitivity").exists()


def test_main_lists_the_analyses(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("COLUMNS", "200")
    with pytest.raises(SystemExit) as info:
        runner(tmp_path).main(["--help"])
    assert info.value.code == 0
    out = capsys.readouterr().out
    assert "local, sampling, fast, sobol" in out
    assert "Toy PK model" in out
    assert "PKPDSIM_SENSITIVITY_N" in out


@pytest.mark.parametrize(
    ("make", "match"),
    [
        (lambda: LocalAnalysis(groups=[]), "at least one group"),
        (lambda: LocalAnalysis(groups=[CONTROL, CONTROL]), r"groups \['control'\]"),
        (lambda: LocalAnalysis(groups=[CONTROL], delta=0), "delta"),
        (lambda: LocalAnalysis(groups=[CONTROL], delta=1.0), r"delta.*\(0, 1\)"),
        (lambda: SamplingAnalysis(groups=[CONTROL], n=0), "number of samples"),
        (lambda: SamplingAnalysis(groups=[CONTROL], n=True), "number of samples"),
        (
            lambda: SamplingAnalysis(groups=[CONTROL], quantiles=(0.5, 1.5)),
            r"quantiles.*1\.5",
        ),
        (lambda: SobolAnalysis(groups=[CONTROL], n=-1), "number of samples"),
        (lambda: SobolAnalysis(groups=[CONTROL], n=1000), "power of two.*1000"),
        (lambda: SobolAnalysis(groups=[CONTROL], n=1), "power of two"),
        (lambda: FastAnalysis(groups=[CONTROL], n=64), "more than 64.*64"),
    ],
    ids=[
        "no-group",
        "same-group",
        "delta",
        "delta-one",
        "sampling-n",
        "sampling-bool",
        "quantiles",
        "sobol-n",
        "sobol-power",
        "sobol-one",
        "fast-n",
    ],
)
def test_an_analysis_is_checked(make: Callable[[], object], match: str) -> None:
    with pytest.raises(ValueError, match=match):
        make()


@pytest.mark.parametrize(
    ("options", "match"),
    [
        ({"outputs": []}, "at least one output"),
        (
            {"outputs": [*OUTPUTS, SensitivityOutput("A_auc", "again")]},
            r"outputs \['A_auc'\]",
        ),
        (
            {"outputs": [*OUTPUTS, SensitivityOutput("A_min2", "x", "A_min")]},
            r"observables \['A_min'\]",
        ),
        (
            {"outputs": [*OUTPUTS, SensitivityOutput("A_max", "A maximum")]},
            r"no observables \['A_max'\]",
        ),
        ({"analyses": {}}, "at least one analysis"),
        ({"default": ("local", "morris")}, r"default analyses \['morris'\]"),
        ({"default": ()}, "at least one default analysis"),
        ({"default": ("local", "local")}, r"default analyses \['local'\].*twice"),
    ],
    ids=[
        "no-output",
        "same-output",
        "same-observable",
        "no-observable",
        "no-analysis",
        "default",
        "no-default",
        "same-default",
    ],
)
def test_a_runner_is_checked(
    tmp_path: Path, options: dict[str, object], match: str
) -> None:
    with pytest.raises(ValueError, match=match):
        runner(tmp_path, **options)


def test_the_default_is_no_string(tmp_path: Path) -> None:
    with pytest.raises(TypeError, match="default"):
        runner(tmp_path, default="local")


def test_the_names_of_a_run_are_no_string(tmp_path: Path) -> None:
    with pytest.raises(TypeError, match="names"):
        runner(tmp_path).run("local", simulator=SERIAL)
    assert not (tmp_path / "sensitivity").exists()


def test_a_name_is_run_once(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match=r"\['local'\].*twice"):
        runner(tmp_path).run(["local", "sampling", "local"], simulator=SERIAL)
    assert not (tmp_path / "sensitivity").exists()


@pytest.mark.parametrize(
    ("value", "names", "match"),
    [
        ("16", ["local", "sampling", "fast"], "more than 64"),
        ("65", ["local", "sobol"], "power of two"),
        ("many", ["local", "sampling"], "PKPDSIM_SENSITIVITY_N"),
    ],
    ids=["fast", "sobol", "malformed"],
)
def test_a_smoke_size_is_checked_before_anything_runs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    value: str,
    names: list[str],
    match: str,
) -> None:
    monkeypatch.setenv("PKPDSIM_SENSITIVITY_N", value)
    with pytest.raises(ValueError, match=match) as info:
        runner(tmp_path).run(names, simulator=SERIAL)
    assert not (tmp_path / "sensitivity").exists()
    if value != "many":
        assert any("PKPDSIM_SENSITIVITY_N" in note for note in info.value.__notes__)


def test_the_distributions_are_checked_before_anything_runs(tmp_path: Path) -> None:
    def reversed_bounds() -> list[SensitivityParameter]:
        return [
            SensitivityParameter(
                uid="k", name="k", value=0.2, lower_bound=0.3, upper_bound=0.1
            ),
            SensitivityParameter(uid="Vc", name="Vc", value=1.0),
        ]

    with pytest.raises(ValueError, match=r"'k' \(0.3 > 0.1\)"):
        runner(tmp_path, parameters=reversed_bounds).run(
            ["local", "sampling"], simulator=SERIAL
        )
    assert not (tmp_path / "sensitivity").exists()


def test_an_analysis_which_varies_nothing_is_refused_before_anything_runs(
    tmp_path: Path,
) -> None:
    fixed = SensitivityGroup("fixed", "Fixed", simulation({"k": 0.1, "Vc": 2.0}))
    analyses = {
        "local": LocalAnalysis(groups=[CONTROL]),
        "fixed": LocalAnalysis(groups=[fixed]),
    }
    with pytest.raises(ValueError, match=r"varies no parameter.*\['fixed'\]") as info:
        runner(tmp_path, analyses=analyses).run(["local", "fixed"], simulator=SERIAL)
    assert info.value.__notes__ == ["The analysis 'fixed' cannot run."]
    assert not (tmp_path / "sensitivity").exists()


def test_the_seed_of_the_designs(tmp_path: Path) -> None:
    analyses = {"sampling": SamplingAnalysis(groups=[CONTROL], n=4)}

    def samples(seed: int, folder: str) -> pd.DataFrame:
        path = runner(tmp_path, analyses=analyses, default=["sampling"], seed=seed).run(
            path=tmp_path / folder, simulator=SERIAL
        )
        return pd.read_csv(path / "sampling" / "control" / "samples.tsv", sep="\t")

    first = samples(1234, "first")
    pd.testing.assert_frame_equal(first, samples(1234, "again"))
    assert not np.allclose(first["k"], samples(4321, "other")["k"])


def test_the_analyses_go_into_a_subfolder_the_table_stays(tmp_path: Path) -> None:
    """The analyses of one of several simulations, the one table of the parameters."""
    analyses = {"local": LocalAnalysis(groups=[CONTROL])}
    toy = runner(tmp_path, analyses=analyses)
    path = toy.run(["local"], subfolder="single", simulator=SERIAL)
    toy.run(["local"], subfolder="multi", simulator=SERIAL)
    assert path == tmp_path / "sensitivity"
    assert sorted(p.name for p in path.iterdir()) == [
        "multi",
        "parameters.tex",
        "parameters.tsv",
        "single",
    ]
    for key in ["single", "multi"]:
        assert (path / key / "local" / "control" / "A_auc.normalized.tsv").is_file()


@pytest.mark.parametrize("subfolder", ["/tmp/elsewhere", "../up", "a/../../up", ""])
def test_a_subfolder_lies_below_the_folder(tmp_path: Path, subfolder: str) -> None:
    with pytest.raises(ValueError, match="subfolder"):
        runner(tmp_path).run(["local"], subfolder=subfolder, simulator=SERIAL)
    assert not (tmp_path / "sensitivity").exists()


def test_the_command_line_with_a_choice_of_the_script(tmp_path: Path) -> None:
    toy = runner(tmp_path, default=("local", "sampling"))
    choices = ["single", "multi"]
    assert toy.parse_command_line([], choices=choices) == (
        ("local", "sampling"),
        (),
    )
    words = toy.parse_command_line(["sobol", "multi", "local"], choices=choices)
    assert words.analyses == ("sobol", "local")
    assert words.choices == ("multi",)
    assert toy.parse_command_line(["single"], choices=choices) == (
        ("local", "sampling"),
        ("single",),
    )


@pytest.mark.parametrize(
    ("argv", "message"),
    [
        (["morris"], "invalid choice: 'morris' (choose from local, sampling, fast, "),
        (["local", "local"], "'local' given twice"),
        (["multi", "multi"], "'multi' given twice"),
    ],
    ids=["unknown", "same-analysis", "same-choice"],
)
def test_the_command_line_refuses_a_word(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], argv: list[str], message: str
) -> None:
    with pytest.raises(SystemExit) as info:
        runner(tmp_path).parse_command_line(argv, choices=["single", "multi"])
    assert info.value.code == 2
    err = capsys.readouterr().err
    assert message in err
    if argv == ["morris"]:
        assert "single, multi" in err


def test_a_choice_of_the_script_is_no_analysis(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match=r"\['local'\]"):
        runner(tmp_path).parse_command_line([], choices=["single", "local"])
    with pytest.raises(TypeError, match="argv"):
        runner(tmp_path).parse_command_line("local")


def test_the_help_names_the_choices_of_the_script(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("COLUMNS", "200")
    with pytest.raises(SystemExit):
        runner(tmp_path).parse_command_line(["--help"], choices=["single", "multi"])
    out = capsys.readouterr().out
    assert "local, sampling, fast, sobol" in out
    assert "single, multi" in out
