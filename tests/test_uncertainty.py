"""Tests of the uncertainty of the fitted parameters: records, distrib, draws, bands."""

import dataclasses
import json
import math
import re
import shutil
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from sbmlsim.data import Data
from sbmlsim.experiment import ExperimentRunner
from sbmlsim.fit.options import ParameterScaleType
from sbmlsim.plot import Curve, Figure
from sbmlsim.plot.plotting import Band
from sbmlsim.report.render import compile_pdf
from sbmlsim.simulation import Dimension, Scan, Simulation
from sbmlsim.simulation.sampling import LogNormal, Normal, Truncated, Uniform
from sbmlsim.simulator import Simulator
from sbmlutils.factory import Parameter, create_model

from pkpdsim.cli import fit_command, run_command
from pkpdsim.experiment import PKPDSimulationExperiment
from pkpdsim.scan import point_band, scan_dimensions
from pkpdsim.uncertainty import (
    RECORD_FILES,
    UNCERTAINTY_DIMENSION,
    FitRecord,
    FitRecords,
    ParameterUncertainty,
    UncertaintyMethod,
    UncertaintyRun,
    UncertaintySettings,
    adopt_fit,
    check_fit_parameters,
    distribution_from_formula,
    parameter_distributions,
    uncertainty_dimension,
)
from pkpdsim_toy import MODEL, TOY_PATH
from pkpdsim_toy.experiments.decay import ToyDecay, ToyDecaySlow
from pkpdsim_toy.factory import U, toy_model


@pytest.fixture(scope="module")
def toy_fit(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Fit the toy model with its profile likelihood and Fisher information."""
    output = tmp_path_factory.mktemp("fit")
    runs = fit_command(
        MODEL,
        [
            "--runs",
            "1",
            "--identifiability",
            "--max-points",
            "5",
            "--no-mapping-figures",
            "--output_dir",
            str(output),
        ],
    )
    (opid,) = runs
    return output / opid


@pytest.fixture(scope="module")
def fits_dir(toy_fit: Path, tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Get a folder of fits with the adopted toy fit."""
    directory = tmp_path_factory.mktemp("fits")
    adopt_fit(toy_fit, directory)
    return directory


@pytest.fixture(scope="module")
def toy_with_uncertainty(
    toy_fit: Path, fits_dir: Path, tmp_path_factory: pytest.TempPathFactory
) -> Path:
    """Create the toy model whose rate is the value of the fit with its uncertainty."""
    records = FitRecords(fits_dir, [toy_fit.name])
    model = toy_model()
    model.parameters = [
        Parameter("k", **records.parameter("k", unit=U.per_min), name="rate")
    ]
    path = tmp_path_factory.mktemp("model") / "toy.xml"
    create_model(model=model, filepath=path, sbml_level=3, sbml_version=2)
    return path


def uncertainty(**changes: object) -> ParameterUncertainty:
    """Get the uncertainty of a parameter with a closed profile in log10."""
    fields: dict[str, object] = {
        "pid": "k",
        "value": 1.0,
        "unit": "1/min",
        "lower_bound": 0.01,
        "upper_bound": 10.0,
        "scale": ParameterScaleType.LOG10,
        "alpha": 0.95,
        "z": 1.959963984540054,
        "method": UncertaintyMethod.PROFILE,
        "ci_lower": 0.5,
        "ci_upper": 2.0,
        "fit_id": "FIT_1",
        "cost": 1.0,
        "record": "FIT_1/",
    }
    fields.update(changes)
    return ParameterUncertainty(**fields)  # ty: ignore[invalid-argument-type]


# ----------------------------------------------------------------------------
# adopting a fit
# ----------------------------------------------------------------------------
def test_adopt_copies_the_record(toy_fit: Path, tmp_path: Path) -> None:
    record = adopt_fit(toy_fit, tmp_path)
    assert record == tmp_path / toy_fit.name
    names = {p.name for p in record.iterdir()}
    assert {"parameters.json", "identifiability.json", "fisher.json"} <= names
    assert names <= set(RECORD_FILES) | {"plots"}


REPORT_FILES = [
    "report.pdf",
    "report.typ",
    "report.tex",
    "report_standalone.tex",
    "references.bib",
    "references.csl",
]


def test_adopt_copies_the_report(toy_fit: Path, tmp_path: Path) -> None:
    record = adopt_fit(toy_fit, tmp_path)
    for name in REPORT_FILES:
        assert (record / name).exists(), name
    referenced = set(
        re.findall(
            r"plots/[\w.\-]+\.(?:png|svg|pdf)",
            (record / "report.typ").read_text(encoding="utf-8")
            + (record / "report.tex").read_text(encoding="utf-8"),
        )
    )
    assert referenced
    kept = {p.relative_to(record).as_posix() for p in (record / "plots").iterdir()}
    # the figures the report shows, nothing else of the plots
    assert kept == referenced


def test_the_report_compiles_in_the_record(toy_fit: Path, tmp_path: Path) -> None:
    record = adopt_fit(toy_fit, tmp_path)
    pdf = compile_pdf(record / "report.typ", tmp_path / "recompiled.pdf")
    assert pdf.stat().st_size > 0


def test_adopt_skips_a_missing_report(toy_fit: Path, tmp_path: Path) -> None:
    fit = tmp_path / "old" / toy_fit.name
    shutil.copytree(toy_fit, fit)
    for name in REPORT_FILES:
        (fit / name).unlink(missing_ok=True)
    record = adopt_fit(fit, tmp_path / "fits")
    assert (record / "parameters.json").exists()
    assert not (record / "report.typ").exists()
    assert not (record / "plots").exists()


def test_adopt_refuses_an_adopted_fit_unless_forced(
    toy_fit: Path, tmp_path: Path
) -> None:
    adopt_fit(toy_fit, tmp_path)
    with pytest.raises(FileExistsError, match="adopted already"):
        adopt_fit(toy_fit, tmp_path)
    assert adopt_fit(toy_fit, tmp_path, force=True).exists()


def test_adopt_refuses_a_fit_without_uncertainty(toy_fit: Path, tmp_path: Path) -> None:
    fit = tmp_path / "FIT"
    fit.mkdir()
    for name in ("parameters.json", "optimization_result.json"):
        shutil.copy2(toy_fit / name, fit / name)
    with pytest.raises(FileNotFoundError, match="neither a profile"):
        adopt_fit(fit, tmp_path / "fits")
    with pytest.raises(FileNotFoundError, match="does not exist"):
        adopt_fit(tmp_path / "missing", tmp_path / "fits")
    (fit / "parameters.json").unlink()
    with pytest.raises(FileNotFoundError, match=r"parameters\.json"):
        adopt_fit(fit, tmp_path / "fits")


def test_fit_adopt_command(toy_fit: Path, tmp_path: Path) -> None:
    model = dataclasses.replace(MODEL, fits_path=tmp_path / "fits")
    runs = fit_command(
        model, ["--adopt", toy_fit.name, "--output_dir", str(toy_fit.parent)]
    )
    assert runs == {}
    assert (tmp_path / "fits" / toy_fit.name / "parameters.json").exists()
    with pytest.raises(SystemExit):
        fit_command(
            model, ["--adopt", toy_fit.name, "--output_dir", str(toy_fit.parent)]
        )


# ----------------------------------------------------------------------------
# the record and the distribution of a parameter
# ----------------------------------------------------------------------------
def test_the_record_of_the_toy_fit(fits_dir: Path, toy_fit: Path) -> None:
    record = FitRecord.load(fits_dir / toy_fit.name)
    k = record.parameters["k"]
    best = json.loads((toy_fit / "parameters.json").read_text())["sets"][0]
    assert k.value == best["values"]["k"]
    assert k.method == UncertaintyMethod.PROFILE
    assert (k.lower_bound, k.upper_bound) == (0.01, 1.0)
    assert k.standard_error is not None
    assert k.fisher_ci is not None
    assert k.n_data == 5
    assert k.identifiability is not None
    # the data of the toy model do not close the profile within its bounds
    expected = "uniform(" if k.sigma() is None else "lognormal("
    assert k.distribution().startswith(expected)
    distribution_from_formula(k.distribution())


def test_a_record_without_profiles_falls_back_on_the_fisher_information(
    fits_dir: Path, toy_fit: Path, tmp_path: Path
) -> None:
    record_dir = tmp_path / toy_fit.name
    shutil.copytree(fits_dir / toy_fit.name, record_dir)
    (record_dir / "identifiability.json").unlink()
    k = FitRecord.load(record_dir).parameters["k"]
    assert k.method == UncertaintyMethod.FISHER
    assert k.identifiability is None
    assert (k.ci_lower, k.ci_upper) == k.fisher_ci
    assert k.sigma() == k.standard_error
    assert k.distribution().startswith("lognormal(")
    (record_dir / "fisher.json").unlink()
    with pytest.raises(FileNotFoundError, match="neither a profile"):
        FitRecord.load(record_dir)


def test_the_distribution_of_a_closed_profile() -> None:
    p = uncertainty()
    sigma = math.log10(2.0) / p.z
    assert p.sigma() == pytest.approx(sigma)
    formula = p.distribution()
    assert formula == (f"lognormal(0.0, {sigma * math.log(10)!r}, 0.01, 10.0)")
    # the 95 % interval of the draws is the profile interval
    d = distribution_from_formula(formula)
    lower, upper = d.ppf(np.array([0.025, 0.975]))
    assert (lower, upper) == pytest.approx((0.5, 2.0), rel=1e-6)


def test_an_asymmetric_profile_takes_the_mean_of_its_half_widths() -> None:
    p = uncertainty(ci_lower=0.1, ci_upper=2.0)
    assert p.sigma() == pytest.approx((1.0 + math.log10(2.0)) / 2 / p.z)


def test_an_open_side_is_truncated_at_the_bound() -> None:
    p = uncertainty(value=3.0, upper_bound=3.0, ci_lower=2.6, ci_upper=None)
    assert p.sigma() == pytest.approx(math.log10(3.0 / 2.6) / p.z)
    assert p.interval() == (2.6, 3.0)
    d = distribution_from_formula(p.distribution())
    assert isinstance(d, Truncated)
    sigma = p.sigma()
    assert sigma is not None
    draws = np.asarray(d.ppf(np.linspace(0.01, 0.99, 50)), dtype=float)
    assert np.all(draws <= 3.0)
    # the upper half of the lognormal is cut off: the median of the draws is
    # its lower quartile
    assert float(np.median(draws)) == pytest.approx(
        3.0 * 10 ** (-0.6745 * sigma), rel=1e-2
    )


def test_a_profile_open_on_both_sides_is_uniform_in_the_bounds() -> None:
    p = uncertainty(ci_lower=None, ci_upper=None, identifiability="non_identifiable")
    assert p.sigma() is None
    assert p.distribution() == "uniform(0.01, 10.0)"
    assert isinstance(distribution_from_formula(p.distribution()), Uniform)


def test_a_linear_scale_gives_a_truncated_normal() -> None:
    p = uncertainty(scale=ParameterScaleType.LINEAR, ci_lower=0.8, ci_upper=1.2)
    assert p.distribution().startswith("normal(1.0, ")
    assert p.sigma() == pytest.approx(0.2 / p.z)
    d = distribution_from_formula(p.distribution())
    assert isinstance(d, Truncated)
    assert isinstance(d.distribution, Normal)


def test_the_fisher_fallback_takes_the_standard_error() -> None:
    p = uncertainty(method=UncertaintyMethod.FISHER, standard_error=0.1)
    assert p.sigma() == 0.1
    assert p.distribution() == (f"lognormal(0.0, {0.1 * math.log(10)!r}, 0.01, 10.0)")
    assert p.standard_error_model() == pytest.approx(math.log(10) * 0.1)


def test_the_notes_have_the_provenance() -> None:
    notes = uncertainty(ci_upper=None).notes()
    for text in ["FIT_1", "profile likelihood", "[0.5, 10]", "open sides", "lognormal"]:
        assert text in notes


@pytest.mark.parametrize(
    ("formula", "kind"),
    [
        ("lognormal(0.0, 0.5)", LogNormal),
        ("normal(1, 0.2)", Normal),
        ("uniform(1, 2)", Uniform),
        ("lognormal(0.0, 0.5, 0.1, 3)", Truncated),
    ],
)
def test_the_distributions_of_distrib(formula: str, kind: type) -> None:
    assert isinstance(distribution_from_formula(formula), kind)


@pytest.mark.parametrize("formula", ["gamma(1, 2)", "normal(1)", "normal(a, 1)", "x"])
def test_an_unsupported_distribution_raises(formula: str) -> None:
    with pytest.raises(ValueError, match="distribution"):
        distribution_from_formula(formula)


def test_the_lognormal_of_distrib_has_its_median_and_sigma() -> None:
    d = distribution_from_formula("lognormal(1.0, 0.5)")
    median, upper = d.ppf(np.array([0.5, 0.8413447460685429]))
    assert median == pytest.approx(math.e)
    assert upper == pytest.approx(math.exp(1.5))


# ----------------------------------------------------------------------------
# the records of a model
# ----------------------------------------------------------------------------
def test_the_records_give_the_arguments_of_a_parameter(
    fits_dir: Path, toy_fit: Path
) -> None:
    records = FitRecords(fits_dir, [toy_fit.name])
    kwargs = records.parameter("k", unit=U.per_min)
    assert set(kwargs) == {"value", "uncertainties", "unit"}
    assert kwargs["value"] == records["k"].value
    assert set(records.parameter("k")) == {"value", "uncertainties"}
    with pytest.raises(KeyError, match="No adopted fit has the parameter 'x'"):
        records.parameter("x")


def test_a_parameter_of_two_fits_raises(fits_dir: Path, toy_fit: Path) -> None:
    shutil.copytree(fits_dir / toy_fit.name, fits_dir / "COPY", dirs_exist_ok=True)
    try:
        with pytest.raises(ValueError, match="one adopted fit"):
            FitRecords(fits_dir, [toy_fit.name, "COPY"])
    finally:
        shutil.rmtree(fits_dir / "COPY")


def test_a_fit_which_is_not_adopted_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="--adopt FIT_X"):
        FitRecords(tmp_path, ["FIT_X"])


def test_the_model_has_the_uncertainty(toy_with_uncertainty: Path) -> None:
    text = toy_with_uncertainty.read_text()
    assert "distrib:uncertainty" in text
    assert 'distrib:type="confidenceInterval"' in text
    assert 'distrib:type="standardError"' in text
    distributions = parameter_distributions(toy_with_uncertainty)
    assert set(distributions) == {"k"}
    # the profile of the toy data is open: uniform in the bounds of the fit
    assert isinstance(distributions["k"], Truncated | Uniform)


def test_the_draws_of_the_model(toy_with_uncertainty: Path) -> None:
    dimension = uncertainty_dimension(toy_with_uncertainty, 8, seed=1)
    assert dimension.id == UNCERTAINTY_DIMENSION
    values = np.asarray(dimension.values["k"], dtype=float)
    assert values.shape == (8,)
    assert np.all((values >= 0.01) & (values <= 1.0))
    again = uncertainty_dimension(toy_with_uncertainty, 8, seed=1)
    assert np.array_equal(values, np.asarray(again.values["k"], dtype=float))


def test_a_model_without_distributions_cannot_be_drawn() -> None:
    with pytest.raises(ValueError, match="No parameter"):
        uncertainty_dimension(MODEL.model_path, 3, seed=1)


def test_the_check_of_the_fit_parameters(
    fits_dir: Path, toy_fit: Path, toy_with_uncertainty: Path
) -> None:
    assert check_fit_parameters(MODEL) == []
    model = dataclasses.replace(
        MODEL,
        fits=(toy_fit.name,),
        fits_path=fits_dir,
        model_path=toy_with_uncertainty,
    )
    assert check_fit_parameters(model) == []
    # the model of the factory has neither the value nor the distribution
    problems = check_fit_parameters(
        dataclasses.replace(model, model_path=MODEL.model_path)
    )
    assert any("has no distribution" in p for p in problems)
    assert any("take it from the record" in p for p in problems)
    missing = dataclasses.replace(model, fits=("FIT_X",))
    assert any("not adopted" in p for p in check_fit_parameters(missing))


def test_the_settings_of_a_run() -> None:
    assert UncertaintySettings().quantiles == (0.05, 0.95)
    with pytest.raises(ValueError, match="two draws"):
        UncertaintySettings(draws=1)
    with pytest.raises(ValueError, match="quantiles"):
        UncertaintySettings(quantiles=(0.9, 0.1))


# ----------------------------------------------------------------------------
# the experiments over the draws
# ----------------------------------------------------------------------------
def draws_run(path: Path, draws: int = 3) -> UncertaintyRun:
    """Get a run over the draws of a model."""
    settings = UncertaintySettings(draws=draws)
    return UncertaintyRun(settings, uncertainty_dimension(path, draws, seed=1))


def initialized(
    experiment: type[PKPDSimulationExperiment], run: UncertaintyRun | None
) -> PKPDSimulationExperiment:
    """Initialize an experiment, over the draws if given, as a run does."""
    settings: dict[str, Any] = {} if run is None else {"uncertainty": run}
    runner = ExperimentRunner(
        experiment_classes=[experiment],
        base_path=TOY_PATH,
        data_path=list(MODEL.data_paths),
        simulator=None,
        on_error="raise",
        **settings,
    )
    (exp,) = runner.experiments.values()
    assert isinstance(exp, PKPDSimulationExperiment)
    return exp


def test_an_experiment_runs_over_the_draws(toy_with_uncertainty: Path) -> None:
    run = draws_run(toy_with_uncertainty)
    exp = initialized(ToyDecay, run)
    assert exp.uncertainty is run
    simulation = exp._simulations["decay"]
    assert isinstance(simulation, Scan)
    assert [d.id for d in simulation.dimensions] == [UNCERTAINTY_DIMENSION]
    (plot,) = exp._figures["fig_decay"].get_plots()
    assert not [c for c in plot.curves if c.y.task_id is not None]
    (band,) = plot.bands
    assert isinstance(band, Band)
    assert band.across == UNCERTAINTY_DIMENSION
    assert band.name == "A"


def test_a_parameter_which_a_simulation_sets_is_not_drawn(
    toy_with_uncertainty: Path,
) -> None:
    exp = initialized(ToyDecaySlow, draws_run(toy_with_uncertainty))
    # the slow decay sets the rate itself, nothing is left to draw: it runs
    # once and keeps its curve
    assert isinstance(exp._simulations["decay"], Simulation)
    (plot,) = exp._figures["fig_decay"].get_plots()
    assert len(plot.curves) == 1
    assert not plot.bands


def test_without_draws_the_experiment_is_unchanged() -> None:
    exp = initialized(ToyDecay, None)
    assert exp.uncertainty is None
    assert isinstance(exp._simulations["decay"], Simulation)


def test_a_wrong_setting_raises() -> None:
    with pytest.raises(TypeError, match="no UncertaintyRun"):
        ExperimentRunner(
            experiment_classes=[ToyDecay],
            base_path=TOY_PATH,
            data_path=list(MODEL.data_paths),
            simulator=None,
            on_error="raise",
            uncertainty=3,
        )


def test_a_scan_gets_the_draws_appended(toy_with_uncertainty: Path) -> None:
    from pkpdsim.uncertainty import _over_draws

    run = draws_run(toy_with_uncertainty)
    scan = Scan(
        Simulation(end=10, steps=10),
        [Dimension("dim_volume", values={"Vc": np.array([1.0, 2.0])})],
    )
    over = _over_draws(scan, run.dimension, "scan", "exp")
    assert isinstance(over, Scan)
    assert [d.id for d in over.dimensions] == ["dim_volume", UNCERTAINTY_DIMENSION]
    clash = Scan(Simulation(end=10, steps=10), [run.dimension])
    with pytest.raises(ValueError, match="id of the draws"):
        _over_draws(clash, run.dimension, "scan", "exp")


def test_the_band_of_a_point(toy_with_uncertainty: Path) -> None:
    run = draws_run(toy_with_uncertainty, draws=5)
    scan = Scan(
        Simulation(end=10, steps=10),
        [
            Dimension("dim_volume", values={"Vc": np.array([1.0, 2.0])}),
            run.dimension,
        ],
    )
    simulator = Simulator()
    res = simulator.run(simulator.load(MODEL.model_path), scan)
    assert scan_dimensions(res, UNCERTAINTY_DIMENSION) == ["dim_volume"]
    time, lower, median, upper = point_band(
        res, "[A]", across=UNCERTAINTY_DIMENSION, dim_volume=0
    )
    assert time.shape == lower.shape == median.shape == upper.shape == (11,)
    assert np.all(lower <= median) and np.all(median <= upper)
    assert lower[-1] < upper[-1]
    with pytest.raises(ValueError, match="no dimension"):
        point_band(res, "[A]", across="draws", dim_volume=0)


def test_the_run_writes_into_the_results_of_the_draws(
    toy_with_uncertainty: Path, tmp_path: Path
) -> None:
    model = dataclasses.replace(MODEL, model_path=toy_with_uncertainty)
    run_command(
        model,
        [
            "-a",
            "simulate",
            "-e",
            "studies",
            "--uncertainty",
            "--draws",
            "3",
            "--results-dir",
            str(tmp_path),
            "--formats",
            "png",
        ],
    )
    figures = tmp_path / "simulation_uncertainty" / "studies" / "_figures"
    assert list(figures.glob("ToyDecay_*.png"))
    assert not (tmp_path / "simulation").exists()


def test_the_options_of_the_draws(tmp_path: Path) -> None:
    with pytest.raises(SystemExit) as error:
        run_command(
            MODEL,
            [
                "-a",
                "simulate",
                "-e",
                "studies",
                "-u",
                "--draws",
                "1",
                "-r",
                str(tmp_path),
            ],
        )
    assert error.value.code == 2


class ToyPhase(ToyDecay):
    """The decay with a curve of A over A, whose x is not the time."""

    def figures(self) -> dict[str, Figure]:
        """Get the figure of the decay and a curve over the concentration."""
        figures = super().figures()
        (plot,) = figures["fig_decay"].get_plots()
        plot.curves.append(
            Curve(
                x=Data(index="[A]", task="task_decay"),
                y=Data(index="[A]", task="task_decay"),
                name="A over A",
            )
        )
        return figures


def test_a_curve_over_another_variable_draws_the_fitted_values(
    toy_with_uncertainty: Path,
) -> None:
    exp = initialized(ToyPhase, draws_run(toy_with_uncertainty))
    (plot,) = exp._figures["fig_decay"].get_plots()
    (band,) = plot.bands
    assert band.y.task_id == "task_decay"
    (curve,) = [c for c in plot.curves if c.y.task_id is not None]
    assert curve.x.task_id == curve.y.task_id == "task_decay_fitted"
    assert curve.name == "A over A (fit)"
    assert isinstance(exp._simulations["decay_fitted"], Simulation)
    assert exp._tasks["task_decay_fitted"].simulation_id == "decay_fitted"


def test_a_parameter_which_the_model_has_not_is_not_drawn(
    toy_with_uncertainty: Path,
) -> None:
    from pkpdsim.uncertainty import _over_draws

    run = draws_run(toy_with_uncertainty)
    values = dict(run.dimension.values)
    dimension = Dimension(
        UNCERTAINTY_DIMENSION,
        values={**values, "missing": np.ones(3)},
    )
    over = _over_draws(Simulation(end=1, steps=1), dimension, "s", "e", {"k"})
    assert isinstance(over, Scan)
    assert set(over.dimensions[0].values) == {"k"}
    # a model with none of the drawn parameters runs the simulation once
    simulation = Simulation(end=1, steps=1)
    assert _over_draws(simulation, dimension, "s", "e", {"x"}) is simulation
