"""The runner of the sensitivity analyses of a script, on the toy model.

The toy model is the decay of A with the rate `k` in the compartment `Vc`; its
concentration does not depend on `Vc`. The group `large` sets `Vc`, so an
analysis over both groups varies only `k`.
"""

import dataclasses
import logging
import math
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from matplotlib.figure import Figure
from sbmlsim import sensitivity
from sbmlsim.sensitivity import SensitivityResult
from sbmlsim.sensitivity.classification import (
    uncertainty_classification,
    uncertainty_classification_symbol,
)
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
    run_subfolders,
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
    assert FastAnalysis(groups=[CONTROL]).n == 1000
    assert SobolAnalysis(groups=[CONTROL]).n == 4096
    # the groups are kept as a tuple
    assert LocalAnalysis(groups=[CONTROL, LARGE]).groups == (CONTROL, LARGE)
    # the statistics of the sampling are those of sbmlsim 0.8.5, no quantiles
    with pytest.raises(TypeError, match="quantiles"):
        SamplingAnalysis(groups=[CONTROL], quantiles=(0.5,))  # ty: ignore[unknown-argument]


def test_local_writes_one_table_and_heatmap_per_group(tmp_path: Path) -> None:
    """The table and heatmap of sbmlsim 0.8.5 are the only result of the local analysis."""
    path = runner(tmp_path).run(["local"], simulator=SERIAL)
    assert path == tmp_path / "sensitivity"
    assert (path / "parameters.tsv").is_file()
    assert (path / "parameters.tex").is_file()
    assert sorted(p.name for p in (path / "local").iterdir()) == [
        "local_d0.01_nvar3_00_control.png",
        "local_d0.01_nvar3_00_control.tsv",
        "local_d0.01_nvar3_01_large.png",
        "local_d0.01_nvar3_01_large.tsv",
    ]


def test_an_analysis_varies_no_parameter_which_a_group_sets(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    with caplog.at_level(logging.INFO, logger="pkpdsim.sensitivity"):
        path = runner(tmp_path).run(["local"], simulator=SERIAL)
    assert "Removing fixed parameter: Vc" in caplog.text
    df = pd.read_csv(path / "local" / "local_d0.01_nvar3_00_control.tsv", sep="\t")
    assert set(df["parameter"]) == {"k"}


def test_sampling_writes_the_samples_of_every_group(tmp_path: Path) -> None:
    """The samples have no table of sbmlsim 0.8.5; its statistics and figures do."""
    path = runner(tmp_path).run(["sampling"], simulator=SERIAL)
    for group in ["control", "large"]:
        group_path = path / "sampling" / group
        assert sorted(p.name for p in group_path.iterdir()) == ["samples.tsv"]
        samples = pd.read_csv(group_path / "samples.tsv", sep="\t")
        assert len(samples) == 8
        assert {"A_auc", "A_thalf", "[A]_min"} <= set(samples.columns)
        # the values of the design, `Vc` is set by the group `large`
        assert "k" in samples.columns
        assert "Vc" not in samples.columns


@pytest.mark.parametrize(
    ("name", "n"), [("fast", 65), ("sobol", 8)], ids=["fast", "sobol"]
)
def test_a_global_analysis_writes_the_tables_of_s1_and_st(
    tmp_path: Path, name: str, n: int
) -> None:
    """The indices have no table of sbmlsim 0.8.5, its figures show them."""
    path = runner(tmp_path).run([name], simulator=SERIAL)
    group_path = path / name / "control"
    expected = []
    for uid in ["A_auc", "A_thalf", "[A]_min"]:
        expected += [f"{uid}.S1.tsv", f"{uid}.ST.tsv"]
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
    assert (path / "local" / "local_d0.01_nvar3_00_control.tsv").is_file()
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

    def run(
        self: SensitivityRunner,
        names: tuple[str, ...],
        *,
        simulator: SimulatorSettings | None = None,
    ) -> Path:
        assert simulator is None
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
        assert (path / key / "local" / "local_d0.01_nvar3_00_control.tsv").is_file()


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
        None,
    )
    words = toy.parse_command_line(["sobol", "multi", "local"], choices=choices)
    assert words.analyses == ("sobol", "local")
    assert words.choices == ("multi",)
    assert toy.parse_command_line(["single"], choices=choices) == (
        ("local", "sampling"),
        ("single",),
        None,
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


# the artifacts of sbmlsim 0.8.5, beside the tables and figures per group


def test_local_writes_the_table_and_heatmap_of_sbmlsim_085(tmp_path: Path) -> None:
    path = runner(tmp_path).run(["local"], simulator=SERIAL)
    local = path / "local"
    assert sorted(p.name for p in local.iterdir() if p.is_file()) == [
        "local_d0.01_nvar3_00_control.png",
        "local_d0.01_nvar3_00_control.tsv",
        "local_d0.01_nvar3_01_large.png",
        "local_d0.01_nvar3_01_large.tsv",
    ]
    df = pd.read_csv(local / "local_d0.01_nvar3_00_control.tsv", sep="\t")
    assert list(df.columns) == [
        "parameter",
        "output",
        "effect",
        "Sij",
        "|Sij|",
        "classification",
        "method",
        "difference",
        "n_var",
    ]
    # `Vc` is set by the group `large`: one parameter, three outputs
    assert sorted(df["output"]) == ["A_auc", "A_thalf", "[A]_min"]
    assert set(df["parameter"]) == {"k"}
    assert df["|Sij|"].tolist() == sorted(df["|Sij|"], reverse=True)
    assert set(df["method"]) == {"LocalSensitivity"}
    assert set(df["difference"]) == {0.01}
    assert set(df["n_var"]) == {3}
    # the mean of the normalized sensitivities at delta/3, 2 delta/3 and delta
    expected = []
    for delta in [0.01 / 3, 0.02 / 3, 0.01]:
        design = sampling.local(
            ["k"], delta=delta, model=MODEL_PATH, simulation=CONTROL.simulation
        )
        result = sensitivity.local(
            create_simulator(SERIAL).run(
                MODEL_PATH,
                Scan(CONTROL.simulation, [design]),
                OBSERVABLES,
                keep=["A_thalf"],
            )
        )
        expected.append(float(result["A_thalf.normalized"].sel(parameter="k")))
    row = df[df["output"] == "A_thalf"].iloc[0]
    assert row["Sij"] == pytest.approx(np.mean(expected), rel=1e-12)
    assert row["|Sij|"] == pytest.approx(abs(row["Sij"]))
    assert row["Sij"] == pytest.approx(-1.0, rel=1e-2)
    assert row["effect"] == "---"
    assert row["classification"] == "high"


def test_n_var_of_a_local_analysis() -> None:
    assert LocalAnalysis(groups=[CONTROL]).n_var == 3
    with pytest.raises(ValueError, match="n_var"):
        LocalAnalysis(groups=[CONTROL], n_var=0)


def test_sampling_writes_the_tables_and_figures_of_sbmlsim_085(tmp_path: Path) -> None:
    path = runner(tmp_path).run(["sampling"], simulator=SERIAL)
    folder = path / "sampling"
    assert sorted(p.name for p in folder.iterdir() if p.is_file()) == [
        "sampling_N8.tex",
        "sampling_N8.tsv",
        "sampling_N8_compact.tsv",
        "sampling_N8_sensitivity_outputs.png",
        "sampling_N8_sensitivity_samples.png",
    ]
    df = pd.read_csv(folder / "sampling_N8.tsv", sep="\t", keep_default_na=False)
    assert list(df.columns) == [
        "gid",
        "gname",
        "uid",
        "name",
        "N",
        "mean",
        "median",
        "std",
        "cv",
        "min",
        "q005",
        "q095",
        "max",
        "U",
        "symbol",
        "classification",
        "unit",
    ]
    assert df["gid"].tolist() == ["control"] * 3 + ["large"] * 3
    assert df["gname"].tolist() == ["Control"] * 3 + ["Large volume"] * 3
    assert df["uid"].tolist() == ["A_auc", "A_thalf", "[A]_min"] * 2
    assert df["name"].tolist()[:3] == ["A AUC∞", "A half-life", "A minimum"]
    assert set(df["N"]) == {8}
    # the units as the observables give them
    assert df["unit"].tolist()[:3] == ["mM*min", "min", "mM"]
    # the statistics of sbmlsim 0.8.5 over the samples of a group
    values = pd.read_csv(folder / "control" / "samples.tsv", sep="\t")["A_thalf"]
    v = values.to_numpy()
    row = df.iloc[1]
    assert row["mean"] == pytest.approx(np.mean(v))
    assert row["median"] == pytest.approx(np.median(v))
    assert row["std"] == pytest.approx(np.std(v))
    assert row["cv"] == pytest.approx(np.std(v) / np.mean(v))
    assert row["min"] == pytest.approx(np.min(v))
    assert row["q005"] == pytest.approx(np.quantile(v, 0.05))
    assert row["q095"] == pytest.approx(np.quantile(v, 0.95))
    assert row["max"] == pytest.approx(np.max(v))
    u = (np.percentile(v, 97.5) - np.percentile(v, 2.5)) / np.percentile(v, 50)
    assert row["U"] == pytest.approx(u)
    assert row["symbol"] == uncertainty_classification_symbol(u)
    assert row["classification"] == uncertainty_classification(u).value

    compact = pd.read_csv(folder / "sampling_N8_compact.tsv", sep="\t")
    assert list(compact.columns) == ["output", "control", "large", "unit"]
    assert compact["output"].tolist() == ["A AUC∞", "A half-life", "A minimum"]
    mean, cv = row["mean"], row["cv"]
    assert compact["control"].tolist()[1] == f"{mean:.3g} ({cv * 100:.1f})"
    tex = (folder / "sampling_N8.tex").read_text(encoding="utf-8")
    assert r"A AUC$\infty$" in tex
    assert r"mM$\cdot$min" in tex


@pytest.mark.parametrize(
    ("name", "prefix"),
    [("fast", "fast_M4_N65"), ("sobol", "sobol_N8")],
    ids=["fast", "sobol"],
)
def test_a_global_analysis_writes_the_figures_of_sbmlsim_085(
    tmp_path: Path, name: str, prefix: str
) -> None:
    path = runner(tmp_path).run([name], simulator=SERIAL)
    folder = path / name
    uids = ["A_auc", "A_thalf", "[A]_min"]
    expected = [f"{prefix}_sensitivity_00_control_{key}.png" for key in ["S1", "ST"]]
    expected += [
        f"{prefix}_sensitivity_00_control_{ko:03}_{uid}.png"
        for ko, uid in enumerate(uids)
    ]
    assert sorted(p.name for p in folder.iterdir() if p.is_file()) == sorted(expected)


def test_the_heatmaps_of_sbmlsim_085_are_labelled_with_the_names(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: dict[str, tuple[list[str], list[str], dict[str, object]]] = {}
    plot_heatmap = sensitivity.plot_heatmap

    def recorded(result: SensitivityResult, index: str, **kwargs: object) -> Figure:
        calls[Path(str(kwargs["path"])).name] = (
            result.parameters,
            result.observables,
            kwargs,
        )
        return plot_heatmap(result, index, **kwargs)  # ty: ignore[invalid-argument-type]

    monkeypatch.setattr(sensitivity, "plot_heatmap", recorded)
    analyses = {
        "local": LocalAnalysis(groups=[CONTROL]),
        "fast": FastAnalysis(groups=[CONTROL], n=65),
    }
    runner(tmp_path, analyses=analyses).run(["local", "fast"], simulator=SERIAL)
    names = ["A AUC∞", "A half-life", "A minimum"]
    parameters, observables, kwargs = calls["local_d0.01_nvar3_00_control.png"]
    assert parameters == ["k: rate of the decay", "Vc: central"]
    assert observables == names
    assert {key: kwargs[key] for key in ["cutoff", "cluster_rows", "cmap"]} == {
        "cutoff": 0.05,
        "cluster_rows": False,
        "cmap": "seismic",
    }
    assert (kwargs["vmin"], kwargs["vmax"]) == (-2.0, 2.0)
    parameters, observables, kwargs = calls["fast_M4_N65_sensitivity_00_control_ST.png"]
    assert observables == names
    assert (kwargs["cmap"], kwargs["vmin"], kwargs["vmax"]) == ("viridis", 0.0, 1.0)
    # the heatmaps of sbmlsim 0.8.5 are the only heatmaps
    assert sorted(calls) == [
        "fast_M4_N65_sensitivity_00_control_S1.png",
        "fast_M4_N65_sensitivity_00_control_ST.png",
        "local_d0.01_nvar3_00_control.png",
    ]


def test_a_heatmap_without_a_parameter_above_the_cutoff_shows_every_parameter(
    tmp_path: Path,
) -> None:
    """A group which sets `k` varies only `Vc`, which the concentrations ignore."""
    slow = SensitivityGroup("slow", "Slow decay", simulation({"k": 0.1}))
    analyses = {"local": LocalAnalysis(groups=[slow])}
    path = runner(tmp_path, analyses=analyses).run(["local"], simulator=SERIAL)
    assert (path / "local" / "local_d0.01_nvar3_00_slow.png").is_file()


def test_the_color_of_a_group() -> None:
    assert CONTROL.color is None
    group = SensitivityGroup("mildRI", "Mild", simulation({}), color="#66c2a4")
    assert group.color == "#66c2a4"


def test_the_names_of_the_outputs_are_unique(tmp_path: Path) -> None:
    outputs = [*OUTPUTS[:2], SensitivityOutput("[A]_min", "A half-life", "A_min")]
    with pytest.raises(ValueError, match=r"names \['A half-life'\]"):
        runner(tmp_path, outputs=outputs)


# several subfolders, one table of the parameters


def subfolder_runners(tmp_path: Path) -> dict[str, SensitivityRunner]:
    return {
        "single": runner(tmp_path, analyses={"local": LocalAnalysis(groups=[CONTROL])}),
        "multi": runner(
            tmp_path,
            outputs=OUTPUTS[:1],
            analyses={"local": LocalAnalysis(groups=[LARGE])},
        ),
    }


def test_the_table_of_several_subfolders_is_written_once(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    import pkpdsim.sensitivity as module

    tables: list[Path] = []
    write = module.write_parameters

    def counted(
        parameters: Sequence[SensitivityParameter], directory: Path
    ) -> pd.DataFrame:
        tables.append(directory)
        return write(parameters, directory)

    monkeypatch.setattr(module, "write_parameters", counted)
    path = run_subfolders(subfolder_runners(tmp_path), ["local"], simulator=SERIAL)
    assert path == tmp_path / "sensitivity"
    assert tables == [path]
    assert capsys.readouterr().out.count("rate of the decay") == 1
    assert sorted(p.name for p in path.iterdir()) == [
        "multi",
        "parameters.tex",
        "parameters.tsv",
        "single",
    ]
    single = pd.read_csv(
        path / "single" / "local" / "local_d0.01_nvar3_00_control.tsv", sep="\t"
    )
    assert sorted(set(single["output"])) == ["A_auc", "A_thalf", "[A]_min"]
    multi = pd.read_csv(
        path / "multi" / "local" / "local_d0.01_nvar3_00_large.tsv", sep="\t"
    )
    assert sorted(set(multi["output"])) == ["A_auc"]


def test_some_of_several_subfolders(tmp_path: Path) -> None:
    path = run_subfolders(
        subfolder_runners(tmp_path), subfolders=["multi"], simulator=SERIAL
    )
    assert sorted(p.name for p in path.iterdir()) == [
        "multi",
        "parameters.tex",
        "parameters.tsv",
    ]


def test_several_subfolders_are_checked_before_anything_runs(tmp_path: Path) -> None:
    runners = subfolder_runners(tmp_path)
    with pytest.raises(ValueError, match=r"\['other'\]"):
        run_subfolders(runners, subfolders=["single", "other"], simulator=SERIAL)
    # the analysis of the second runner varies nothing
    fixed = SensitivityGroup("fixed", "Fixed", simulation({"k": 0.1, "Vc": 2.0}))
    runners["multi"] = runner(
        tmp_path, analyses={"local": LocalAnalysis(groups=[fixed])}
    )
    with pytest.raises(ValueError, match="varies no parameter") as info:
        run_subfolders(runners, ["local"], simulator=SERIAL)
    assert "The analysis 'local' of the subfolder 'multi' cannot run." in (
        info.value.__notes__
    )
    with pytest.raises(ValueError, match="no runner"):
        run_subfolders({}, simulator=SERIAL)
    with pytest.raises(ValueError, match="at least one subfolder"):
        run_subfolders(runners, subfolders=(), simulator=SERIAL)
    assert not (tmp_path / "sensitivity").exists()


def test_several_subfolders_share_the_model_and_the_parameters(tmp_path: Path) -> None:
    runners = subfolder_runners(tmp_path)

    def other() -> list[SensitivityParameter]:
        parameters = toy_parameters()
        parameters[0].upper_bound *= 1.01
        return parameters

    runners["multi"] = runner(tmp_path, parameters=other)
    with pytest.raises(ValueError, match="parameters"):
        run_subfolders(runners, ["local"], simulator=SERIAL)
    runners["multi"] = dataclasses.replace(
        runners["single"], model=dataclasses.replace(MODEL, title="other")
    )
    with pytest.raises(ValueError, match="model"):
        run_subfolders(runners, ["local"], simulator=SERIAL)
    assert not (tmp_path / "sensitivity").exists()


def test_an_undefined_or_negative_index_is_classified(tmp_path: Path) -> None:
    """sbmlsim's classifications refuse `NaN` and a negative uncertainty.

    An output which is zero at the reference has an undefined normalized
    sensitivity and uncertainty, a negative output a negative uncertainty,
    which is classified by its magnitude.
    """
    observables = [
        *OBSERVABLES,
        Formula("A_zero", "max([A]) - max([A])", unit="mM"),
        Formula("A_negative", "-pk_A.thalf", unit="min"),
    ]
    outputs = [
        *OUTPUTS,
        SensitivityOutput("A_zero", "A zero"),
        SensitivityOutput("A_negative", "A negative half-life"),
    ]
    analyses = {
        "local": LocalAnalysis(groups=[CONTROL]),
        "sampling": SamplingAnalysis(groups=[CONTROL], n=8),
    }
    path = runner(
        tmp_path, outputs=outputs, observables=observables, analyses=analyses
    ).run(["local", "sampling"], simulator=SERIAL)
    local = pd.read_csv(
        path / "local" / "local_d0.01_nvar3_00_control.tsv",
        sep="\t",
        keep_default_na=False,
        na_values=[""],
    )
    zero = local[(local["output"] == "A_zero") & (local["parameter"] == "k")].iloc[0]
    assert math.isnan(zero["Sij"])
    assert math.isnan(zero["effect"])
    assert math.isnan(zero["classification"])
    df = pd.read_csv(
        path / "sampling" / "sampling_N8.tsv",
        sep="\t",
        keep_default_na=False,
        na_values=[""],
    ).set_index("uid")
    assert math.isnan(df.loc["A_zero", "U"])
    assert math.isnan(df.loc["A_zero", "classification"])
    u = df.loc["A_negative", "U"]
    assert u < 0
    assert df.loc["A_negative", "classification"] == (
        uncertainty_classification(abs(u)).value
    )
    assert df.loc["A_negative", "U"] == pytest.approx(-df.loc["A_thalf", "U"])


@pytest.mark.parametrize(
    ("argv", "workers"),
    [
        ([], None),
        (["--workers", "0"], None),
        (["--workers", "1"], 1),
        (["local", "--workers", "4", "multi"], 4),
        (["--workers=3", "sobol"], 3),
    ],
)
def test_the_command_line_reads_the_workers(
    tmp_path: Path, argv: list[str], workers: int | None
) -> None:
    words = runner(tmp_path).parse_command_line(argv, choices=["multi"])
    assert words.workers == workers


@pytest.mark.parametrize("value", ["-1", "two", "1.5", ""])
def test_the_command_line_refuses_workers(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], value: str
) -> None:
    with pytest.raises(SystemExit) as info:
        runner(tmp_path).parse_command_line(["--workers", value])
    assert info.value.code == 2
    assert "--workers" in capsys.readouterr().err


def test_the_help_names_the_workers(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("COLUMNS", "200")
    with pytest.raises(SystemExit):
        runner(tmp_path).main(["--help"])
    assert "--workers" in capsys.readouterr().out


def test_the_workers_set_the_processes_of_the_simulator(tmp_path: Path) -> None:
    toy = runner(tmp_path)
    assert toy.simulator_settings(None) is None
    assert toy.simulator_settings(4) == dataclasses.replace(
        MODEL.simulator, n_workers=4
    )


def test_main_passes_the_workers_to_the_run(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: list[tuple[tuple[str, ...], SimulatorSettings | None]] = []

    def run(
        self: SensitivityRunner,
        names: tuple[str, ...],
        *,
        simulator: SimulatorSettings | None = None,
    ) -> Path:
        calls.append((tuple(names), simulator))
        return tmp_path

    monkeypatch.setattr(SensitivityRunner, "run", run)
    toy = runner(tmp_path)
    toy.main(["local", "--workers", "2"])
    toy.main(["local", "--workers", "0"])
    toy.main(["local"])
    assert calls == [
        (("local",), dataclasses.replace(MODEL.simulator, n_workers=2)),
        (("local",), None),
        (("local",), None),
    ]
