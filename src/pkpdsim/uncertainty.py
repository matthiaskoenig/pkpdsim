"""The bands of the predictions of a model over the uncertainty of its fitted parameters.

A fit of a model becomes part of the model in three steps:

1. **Adopt** the fit: `run_fit <package> --adopt <fit id>` copies the record of a
   finished fit (its parameters, the profile likelihood and the Fisher
   information) from the results of the fits into the tracked folder of the
   fits next to the model (`PKPDModel.fits_dir`), see
   `sbmlsim.fit.records.adopt_fit`.
2. **Use** it in the factory: the model lists its adopted fits
   (`PKPDModel.adopted_ids`), and a model file takes the value of a fitted
   parameter from the record, together with its uncertainty, see
   `sbmlsim.fit.records.AdoptedRecords.parameter`:

   ```python
   # not executed
   FITS = MODEL.fit_records()
   Parameter("k", **FITS.parameter("k", unit=U.per_min), name="rate")
   ```

   `check_fit_parameters` checks that the model has the values and the
   distributions of the records.
3. **Draw** it: `run_simulation <package> --uncertainty` reads the distributions
   of the parameters of the model, draws them into the dimension `uncertainty`
   (`sbmlsim.simulation.sampling.uncertainty_dimension`) and runs every
   simulation over the draws; every curve of a simulation becomes the band of a
   quantile range with the median, see `apply_uncertainty`.

The draws are independent per parameter: distrib states one distribution per
parameter, the correlations of the Fisher information stay in the record.
"""

from __future__ import annotations

import logging
import math
from collections.abc import Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import libsbml
from sbmlsim.data import Data
from sbmlsim.plot.plotting import Curve
from sbmlsim.simulation import Dimension, Scan, Simulation
from sbmlsim.simulation.sampling import (
    Distribution,
    parameter_distributions,
    uncertainty_dimension,
)
from sbmlsim.task import Task

from pkpdsim.sex import Sex

if TYPE_CHECKING:
    from sbmlsim.experiment import SimulationExperiment

    from pkpdsim.model import PKPDModel

logger = logging.getLogger(__name__)

#: the relative tolerance of the value of a fitted parameter in the model
VALUE_RTOL = 1e-12


# ----------------------------------------------------------------------------
# the experiments over the draws
# ----------------------------------------------------------------------------
@dataclass(frozen=True, kw_only=True)
class UncertaintySettings:
    """The settings of a run over the draws of the fitted parameters.

    Attributes:
        draws: the number of draws.
        seed: the seed of the draws.
        quantiles: the lower and the upper quantile of the bands.
    """

    draws: int = 100
    seed: int = 1234
    quantiles: tuple[float, float] = (0.05, 0.95)

    def __post_init__(self) -> None:
        """Check the settings.

        Raises:
            ValueError: for less than two draws or quantiles which are not
                increasing within `(0, 1)`.
        """
        object.__setattr__(self, "quantiles", tuple(float(q) for q in self.quantiles))
        if self.draws < 2:
            raise ValueError(f"A band needs at least two draws, not {self.draws}.")
        low, high = self.quantiles
        if not 0.0 < low < high < 1.0:
            raise ValueError(
                f"The quantiles of a band are increasing within (0, 1), not "
                f"{self.quantiles}."
            )


@dataclass(frozen=True)
class UncertaintyRun:
    """The draws and the settings of a run, a setting of the experiments.

    Attributes:
        settings: the settings of the run.
        dimension: the draws of the fitted parameters.
    """

    settings: UncertaintySettings
    dimension: Dimension

    @property
    def quantiles(self) -> tuple[float, float]:
        """The quantiles of the bands."""
        return self.settings.quantiles

    @property
    def id(self) -> str:
        """The id of the dimension of the draws."""
        return self.dimension.id

    @classmethod
    def create(cls, model: PKPDModel, settings: UncertaintySettings) -> UncertaintyRun:
        """Draw the fitted parameters of a model.

        The draws are those of the male build (`model_path`); every build has
        the same fitted parameters with the same distributions
        (`check_fit_parameters`), so the draws run on the experiments of every
        build.

        Args:
            model: the model, whose male build (`model_path`) has the
                distributions.
            settings: the settings of the run.

        Returns:
            The run.
        """
        return cls(
            settings=settings,
            dimension=uncertainty_dimension(
                model.model_path, settings.draws, settings.seed
            ),
        )


def _set_targets(simulation: Simulation) -> set[str]:
    """Get the targets a simulation sets itself."""
    targets = set(simulation.preinit_changes)
    if simulation.presimulation is not None:
        targets |= set(simulation.presimulation.preinit_changes)
    for change in simulation.changes:
        targets |= set(change.values)
    return targets


def _without(dimension: Dimension, targets: set[str]) -> Dimension | None:
    """Get the dimension without the values of some targets, `None` if none is left."""
    shared = targets & set(dimension.values)
    if not shared:
        return dimension
    values = {k: v for k, v in dimension.values.items() if k not in shared}
    if not values:
        return None
    return Dimension(
        dimension.id, values=values, labels=dimension.labels, design=dimension.design
    )


def _over_draws(
    simulation: Simulation | Scan,
    dimension: Dimension,
    key: str,
    sid: str,
    available: set[str] | None = None,
) -> Simulation | Scan:
    """Get a simulation or a scan over the draws of the fitted parameters.

    A parameter which the simulation, or a dimension of the scan, sets itself
    keeps its value: it is left out of the draws of this simulation, as is a
    parameter which the model of the simulation does not have.

    Args:
        simulation: the simulation or scan.
        dimension: the draws.
        key: the key of the simulation, for the messages.
        sid: the id of the experiment, for the messages.
        available: the parameters of the model of the simulation, `None` for
            all of the draws.

    Returns:
        The scan over the draws, the simulation itself if no drawn parameter
        is left.

    Raises:
        ValueError: if the scan has a dimension of the id of the draws.
    """
    if isinstance(simulation, Scan):
        if any(d.id == dimension.id for d in simulation.dimensions):
            raise ValueError(
                f"The scan '{key}' of the experiment '{sid}' has a dimension "
                f"'{dimension.id}', the id of the draws of the fitted parameters."
            )
        base = simulation.simulation
        targets = _set_targets(base)
        for d in simulation.dimensions:
            targets |= set(d.values)
            for s in d.simulations.values():
                targets |= _set_targets(s)
        dimensions = list(simulation.dimensions)
    else:
        base = simulation
        targets = _set_targets(simulation)
        dimensions = []
    left_out = sorted(targets & set(dimension.values))
    if left_out:
        logger.info(
            "The simulation '%s' of '%s' sets the fitted parameters %s itself, "
            "they are not drawn.",
            key,
            sid,
            left_out,
        )
    excluded = set(targets)
    if available is not None:
        excluded |= set(dimension.values) - available
    draws = _without(dimension, excluded)
    if draws is None:
        return simulation
    return Scan(base, [*dimensions, draws])


def _has_dimension(simulation: Simulation | Scan, dimension_id: str) -> bool:
    """Check whether a simulation is a scan with a dimension of an id."""
    return isinstance(simulation, Scan) and any(
        d.id == dimension_id for d in simulation.dimensions
    )


def _curve_color(curve: Curve) -> str | None:
    """Get the color of a curve, `None` for the color cycle."""
    color = (curve.kwargs or {}).get("color")
    if color is None:
        line = getattr(curve.style, "line", None)
        line_color = getattr(line, "color", None)
        color = getattr(line_color, "color", line_color)
    return None if color is None else str(color)


def _model_parameters(
    experiment: SimulationExperiment, model_id: str
) -> set[str] | None:
    """Get the parameters of a model of an experiment, `None` if it is not loaded."""
    model = experiment._models.get(model_id)
    runner = getattr(model, "r", None)
    if runner is None:
        return None
    return set(runner.model.getGlobalParameterIds())


#: the suffix of the simulation and the task of the fitted values, which the
#: curves draw whose x is neither the time nor the values of a dimension of the scan
FITTED_SUFFIX = "_fitted"


def _fitted_task(
    experiment: SimulationExperiment,
    task_id: str,
    originals: Mapping[str, Simulation | Scan],
) -> str:
    """Get the task of the simulation of a task with the fitted values, created once.

    Returns:
        The id of the task.
    """
    task = experiment._tasks[task_id]
    fitted_task = f"{task_id}{FITTED_SUFFIX}"
    if fitted_task not in experiment._tasks:
        fitted_simulation = f"{task.simulation_id}{FITTED_SUFFIX}"
        experiment._simulations[fitted_simulation] = originals[task.simulation_id]
        experiment._tasks[fitted_task] = Task(
            model=task.model_id, simulation=fitted_simulation
        )
    return fitted_task


def _on_task(data: Data | None, task_id: str) -> Data | None:
    """Get the data of a task from another task."""
    if data is None or data.task_id is None:
        return data
    return Data(data.selection, task=task_id, sel=data.sel)


def _over_scan_dimension(data: Data, simulation: Simulation | Scan) -> bool:
    """Check whether data are the values of a dimension of a scan.

    The values of a dimension (`<dimension>.<target>`) are the same in every
    draw, a curve over them has a band across the draws.
    """
    if not isinstance(simulation, Scan):
        return False
    dimension = str(data.index).split(".", 1)[0]
    return any(d.id == dimension for d in simulation.dimensions)


def _draw_data(data: Data, run_id: str, k: int) -> Data:
    """Get the data of a draw, an id of its own (the id does not depend on `sel`)."""
    sid = Data(data.selection, task=data.task_id, sel=data.sel).sid
    return Data(
        data.selection,
        task=data.task_id,
        sel={**data.sel, run_id: k},
        sid=f"{sid}__{run_id}_{k}",
    )


def _draw_curves(curve: Curve, run: UncertaintyRun) -> list[Curve]:
    """Get a line of every draw of a curve over another variable than the time.

    x differs between the draws, so the quantiles of y have no common x: the
    lines of the draws, thin, translucent and without a legend entry, in the
    colour of the curve, give the area of the draws; their opacity falls with
    the number of draws, so that the area is denser where more draws lie.
    """
    draws = run.settings.draws
    kwargs: dict[str, Any] = {
        **(curve.kwargs or {}),
        "label": None,
        "alpha": min(0.3, 10 / draws),
        "linewidth": 0.5,
        "marker": "",
    }
    kwargs.pop("markeredgecolor", None)
    color = _curve_color(curve)
    if color is not None:
        kwargs["color"] = color
    return [
        Curve(
            x=_draw_data(curve.x, run.id, k),
            y=_draw_data(curve.y, run.id, k),
            yaxis_position=curve.yaxis_position,
            **kwargs,
        )
        for k in range(draws)
    ]


def apply_uncertainty(experiment: SimulationExperiment, run: UncertaintyRun) -> None:
    """Run an initialized experiment over the draws of the fitted parameters.

    Every simulation becomes a scan over the draws (a scan gets the dimension
    of the draws appended); a parameter which the simulation sets itself, or
    which the model of the simulation does not have, is not drawn, and a
    simulation without a drawn parameter runs once. Every curve of a figure
    over the time of a task over the draws, or over the values of a dimension
    of its scan (e.g. an amount at a time over the scanned renal function),
    becomes the band of the quantiles of the draws with the median, with the
    name and the colour of the curve. A curve whose x is another variable of
    the task (e.g. an effect over a concentration) has no common x across the
    draws: it draws the simulation with the fitted values, which runs once next
    to the draws (`<simulation>_fitted`), over a thin translucent line of every
    draw in its colour (the area of the draws, without a legend entry; not for
    a curve with a line per point of a dimension). The curves of the data stay.

    Args:
        experiment: the initialized experiment.
        run: the draws and the settings.

    Raises:
        ValueError: if a scan of the experiment has a dimension of the id of
            the draws, or a curve of a task draws a line per point of more
            than one dimension (a band is drawn per point of one).
    """
    sid = experiment.sid
    simulations = experiment._simulations
    originals = dict(simulations)
    available: dict[str, set[str] | None] = {}
    for task in experiment._tasks.values():
        parameters = _model_parameters(experiment, task.model_id)
        known = available.get(task.simulation_id)
        if task.simulation_id not in available or known is None:
            available[task.simulation_id] = parameters
        elif parameters is not None:
            available[task.simulation_id] = known & parameters
    for key, simulation in originals.items():
        simulations[key] = _over_draws(
            simulation, run.dimension, key, sid, available.get(key)
        )
    # the tasks whose simulation runs over the draws; a simulation without a
    # drawn parameter runs once and keeps its curves
    drawn = {
        task_id
        for task_id, task in experiment._tasks.items()
        if _has_dimension(simulations[task.simulation_id], run.id)
    }
    for figure in experiment._figures.values():
        for plot in figure.get_plots():
            curves: list[Curve] = []
            # the lines of the draws come first, below every curve of the plot
            draw_curves: list[Curve] = []
            for curve in plot.curves:
                task_id = None if curve.y.task_id is None else str(curve.y.task_id)
                if task_id is None or task_id not in drawn:
                    curves.append(curve)
                    continue
                original = originals[experiment._tasks[task_id].simulation_id]
                if curve.x.index != "time" and not _over_scan_dimension(
                    curve.x, original
                ):
                    if not curve.over:
                        draw_curves.extend(_draw_curves(curve, run))
                    fitted = _fitted_task(experiment, task_id, originals)
                    curve.x = _on_task(curve.x, fitted)  # ty: ignore[invalid-assignment]
                    curve.y = _on_task(curve.y, fitted)  # ty: ignore[invalid-assignment]
                    curve.xerr = _on_task(curve.xerr, fitted)
                    curve.yerr = _on_task(curve.yerr, fitted)
                    curve.name = f"{curve.name} (fit)" if curve.name else curve.name
                    curves.append(curve)
                    continue
                if len(curve.over) > 1:
                    raise ValueError(
                        f"The curve '{curve.sid}' of the figure '{figure.sid}' of "
                        f"'{sid}' draws a line per point of {list(curve.over)}; a "
                        f"band is drawn per point of one dimension."
                    )
                plot.band(
                    x=curve.x,
                    y=curve.y,
                    across=run.id,
                    quantiles=run.quantiles,
                    over=curve.over or None,
                    name=curve.name,
                    color=_curve_color(curve),
                    yaxis_position=curve.yaxis_position,
                )
            plot.curves = draw_curves + curves


# ----------------------------------------------------------------------------
# the check of a model
# ----------------------------------------------------------------------------
def check_fit_parameters(model: PKPDModel) -> list[str]:
    """Check that the model has the fitted values and the uncertainty of its adopted fits.

    Every parameter of an adopted fit must be a parameter of every build of the
    model (`model_paths`, the flat models) with the value of the fit and a
    distribution of its draws. The draws are those of the male build
    (`UncertaintyRun.create`) and run on every build, so the distribution of a
    parameter in another build is the one of the male build.

    Args:
        model: the model.

    Returns:
        The problems, one line each; empty for a model without adopted fits.
    """
    if not model.adopted_ids:
        return []
    try:
        records = model.fit_records()
    except (FileNotFoundError, ValueError) as err:
        return [str(err)]
    problems: list[str] = []
    male: dict[str, Distribution] | None = None
    for sex, path in model.model_paths.items():
        build = "model" if len(model.model_paths) == 1 else f"{sex} build"
        if not path.exists():
            problems.append(f"the {build} {path} does not exist")
            continue
        sbml_model = libsbml.readSBMLFromFile(str(path)).getModel()
        distributions = parameter_distributions(path)
        if sex is Sex.MALE:
            male = distributions
        where = "" if len(model.model_paths) == 1 else f" (the {build} {path.name})"
        for pid, parameter in records.parameters.items():
            p = sbml_model.getParameter(pid)
            if p is None:
                problems.append(
                    f"the parameter '{pid}' of the fit '{parameter.fit_id}' is no "
                    f"parameter of {path.name}{where}"
                )
                continue
            if not math.isclose(p.getValue(), parameter.value, rel_tol=VALUE_RTOL):
                problems.append(
                    f"the parameter '{pid}' is {p.getValue()!r} in the model, the "
                    f"fit '{parameter.fit_id}' gives {parameter.value!r}: take it "
                    f"from the record (`AdoptedRecords.parameter`){where}"
                )
            if pid not in distributions:
                problems.append(
                    f"the parameter '{pid}' of the fit '{parameter.fit_id}' has no "
                    f"distribution in the model{where}"
                )
            elif male is not None and pid in male and distributions[pid] != male[pid]:
                problems.append(
                    f"the distribution of the parameter '{pid}' of the fit "
                    f"'{parameter.fit_id}' differs from the male build{where}: "
                    f"{distributions[pid]} instead of {male[pid]}"
                )
    return problems
