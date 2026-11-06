"""Tests of the uncertainty of the fitted parameters of a model of pkpdsim.

The records of the adopted fits, their distributions and the run of the
experiments over the draws are sbmlsim's (`sbmlsim.fit.records`,
`sbmlsim.simulation.sampling`, `sbmlsim.experiment.draws`), tested in sbmlsim;
these tests check them on the builds of the sexes, the experiments of pkpdsim and
the command `run_simulation --uncertainty`.
"""

import dataclasses
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from sbmlsim.distrib import UNCERTAINTY_DIMENSION
from sbmlsim.experiment import ExperimentRunner
from sbmlsim.experiment.draws import (
    UncertaintyRun,
    UncertaintySettings,
    check_fit_parameters,
)
from sbmlsim.fit.records import AdoptedRecords
from sbmlsim.plot.plotting import Band
from sbmlsim.result.scan import point_timecourse
from sbmlsim.simulation import Scan
from sbmlsim.simulation.sampling import uncertainty_dimension
from sbmlsim.simulator import Simulator
from sbmlutils.factory import Parameter, create_model

from pkpdsim.cli import fit_command, run_command
from pkpdsim.experiment import PKPDSimulationExperiment
from pkpdsim.model import PKPDModel
from pkpdsim_toy import MODEL, TOY_PATH
from pkpdsim_toy.experiments.decay import ToyDecay
from pkpdsim_toy.factory import U, toy_model


@pytest.fixture(scope="module")
def toy_with_uncertainty(
    toy_fit: Path, fits_dir: Path, tmp_path_factory: pytest.TempPathFactory
) -> Path:
    """Create the toy model whose rate is the value of the fit with its uncertainty."""
    records = AdoptedRecords(fits_dir, [toy_fit.name])
    model = toy_model()
    model.parameters = [
        Parameter("k", **records.parameter("k", unit=U.per_min), name="rate")
    ]
    path = tmp_path_factory.mktemp("model") / "toy.xml"
    create_model(model=model, filepath=path, sbml_level=3, sbml_version=2)
    return path


def test_fit_adopt_command(
    toy_fit: Path, tmp_path: Path, as_model: Callable[[PKPDModel], PKPDModel]
) -> None:
    model = as_model(MODEL.replace(fits_dir=tmp_path / "fits"))
    runs = fit_command(
        model, ["--adopt", toy_fit.name, "--output_dir", str(toy_fit.parent)]
    )
    assert runs == {}
    assert (tmp_path / "fits" / toy_fit.name / "parameters.json").exists()
    # the run of the kind `adopt` in the folder of the fit
    assert (toy_fit / "adopt" / "fits" / toy_fit.name / "parameters.json").exists()
    with pytest.raises(SystemExit):
        fit_command(
            model, ["--adopt", toy_fit.name, "--output_dir", str(toy_fit.parent)]
        )


def test_the_check_of_the_fit_parameters(
    fits_dir: Path, toy_fit: Path, toy_with_uncertainty: Path
) -> None:
    assert check_fit_parameters(MODEL) == []
    model = MODEL.replace(
        adopted_ids=(toy_fit.name,),
        fits_dir=fits_dir,
        variants={"male": toy_with_uncertainty},
        factory_dir=MODEL.factory_dir,
    )
    assert check_fit_parameters(model) == []
    # the model of the factory has neither the value nor the distribution
    problems = check_fit_parameters(model.replace(variants=MODEL.variants))
    assert any("has no distribution" in p for p in problems)
    assert any("take it from the record" in p for p in problems)
    missing = model.replace(adopted_ids=("DECAY_20260101_000000__abcde",))
    assert any("not adopted" in p for p in check_fit_parameters(missing))


def _female_build(
    path: Path,
    records: AdoptedRecords,
    *,
    upper: float | None = None,
    fitted: bool = True,
) -> Path:
    """Create a female build of the toy model with the fitted rate of the records.

    Args:
        path: the SBML file.
        records: the adopted fits.
        upper: another upper bound of the distribution of the draws.
        fitted: give the rate the value and the uncertainty of the fit.

    Returns:
        The SBML file.
    """
    model = toy_model()
    model.sid = "toy_female"
    if fitted:
        kwargs = records.parameter("k", unit=U.per_min)
        if upper is not None:
            uncertainty = dataclasses.replace(records["k"], upper_bound=upper)
            kwargs["uncertainties"] = [uncertainty.uncertainty(U.per_min)]
        model.parameters = [Parameter("k", **kwargs, name="rate")]
    create_model(model=model, filepath=path, sbml_level=3, sbml_version=2)
    return path


def test_the_check_of_the_fit_parameters_sees_every_build(
    fits_dir: Path, toy_fit: Path, toy_with_uncertainty: Path, tmp_path: Path
) -> None:
    """The female build, a variant, has the fitted values and distributions of the male one."""
    records = AdoptedRecords(fits_dir, [toy_fit.name])
    model = MODEL.replace(
        adopted_ids=(toy_fit.name,),
        fits_dir=fits_dir,
        variants={
            "male": toy_with_uncertainty,
            "female": _female_build(tmp_path / "female.xml", records),
        },
        factory_dir=MODEL.factory_dir,
    )
    assert check_fit_parameters(model) == []
    stale = model.replace(
        variants={
            **model.variants,
            "female": _female_build(tmp_path / "stale.xml", records, fitted=False),
        },
    )
    problems = check_fit_parameters(stale)
    assert problems
    assert all("(the variant 'female', stale.xml)" in p for p in problems)
    assert any("has no distribution" in p for p in problems)
    assert any("take it from the record" in p for p in problems)
    other = model.replace(
        variants={
            **model.variants,
            "female": _female_build(tmp_path / "other.xml", records, upper=5.0),
        },
    )
    (problem,) = check_fit_parameters(other)
    assert "distribution" in problem
    assert "differs from the default variant" in problem
    missing = model.replace(
        variants={**model.variants, "female": tmp_path / "missing.xml"}
    )
    assert check_fit_parameters(missing) == [
        f"the variant 'female' {tmp_path / 'missing.xml'} does not exist"
    ]


def test_the_draws_reach_a_female_experiment(toy_with_uncertainty: Path) -> None:
    """The draws of the male build run on the female build, which has the same ids.

    Every draw of the rate decays from the female start: A(t) = 0.8 exp(-k t).
    """
    from pkpdsim_sexes import FEMALE_A0, FEMALE_MODEL_PATH
    from pkpdsim_sexes.experiments.decay import SexesWomen

    run = draws_run(toy_with_uncertainty)
    runner = ExperimentRunner(
        experiment_classes=[SexesWomen],
        base_path=TOY_PATH,
        simulator=Simulator(n_workers=1),
        on_error="raise",
        uncertainty=run,
    )
    (exp,) = runner.experiments.values()
    assert exp._models["model"].source.path == FEMALE_MODEL_PATH
    simulation = exp._simulations["decay"]
    assert isinstance(simulation, Scan)
    assert [d.id for d in simulation.dimensions] == [UNCERTAINTY_DIMENSION]
    exp.run(runner.simulator, output_path=None)
    result = exp.results["task_decay"]
    rates = np.asarray(run.dimension.values["k"], dtype=float)
    assert len(set(rates)) == len(rates) == 3
    for index, k in enumerate(rates):
        times, values = point_timecourse(
            result, "[A]", **{UNCERTAINTY_DIMENSION: index}
        )
        np.testing.assert_allclose(
            np.asarray(values),
            FEMALE_A0 * np.exp(-k * np.asarray(times)),
            rtol=1e-6,
            atol=1e-8,
        )


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


def test_the_run_writes_into_the_results_of_the_draws(
    toy_with_uncertainty: Path,
    tmp_path: Path,
    as_model: Callable[[PKPDModel], PKPDModel],
) -> None:
    model = as_model(
        MODEL.replace(
            variants={"male": toy_with_uncertainty}, factory_dir=MODEL.factory_dir
        )
    )
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
            "-w",
            "1",
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
