"""The sensitivity analyses of a model: the table of their parameters and their runner.

sbmlsim 0.9 analyses the sensitivity of a model on the result of a scan: a
design of `sbmlsim.simulation.sampling` (`local`, `sobol`, `fast`, `lhs`) is a
dimension of a `Scan`, and `sbmlsim.sensitivity` computes the indices on its
result. It dropped the table of the parameters of an analysis
(`sbmlsim.sensitivity.parameters` of sbmlsim 0.8), which the sensitivity
scripts of the models write next to their results. This module is that table,
with the behavior of sbmlsim 0.8, so the tables of the scripts stay the same,
and the distributions of the global designs:

```python
# not executed
parameters = parameters_from_sbml(MODEL_PATH, exclude_ids={"Mr_hctz"})
set_bounds(parameters, [("f_renal_function", 0.85, 1.15, ParameterType.SCALING)])
set_relative_bounds(parameters, 0.15)
write_parameters(parameters, MODEL.results_dir / "sensitivity")
design = sobol(distributions(parameters), n=sample_size(4096), seed=1234, model=MODEL_PATH)
```

`SensitivityRunner` runs the analyses of a script (`LocalAnalysis`,
`SamplingAnalysis`, `FastAnalysis`, `SobolAnalysis`) over the simulations of
its groups (`SensitivityGroup`) and writes the tables and figures of its
outputs (`SensitivityOutput`); the script keeps only what is the model's own,
its parameters, groups, outputs and observables:

```python
# not executed
SENSITIVITY = SensitivityRunner(
    model=MODEL,
    parameters=sensitivity_parameters,
    outputs=OUTPUTS,
    observables=OBSERVABLES,
    analyses={
        "local": LocalAnalysis(groups=GROUPS[:1]),
        "sampling": SamplingAnalysis(groups=GROUPS),
    },
)

if __name__ == "__main__":
    SENSITIVITY.main()
```
"""

from __future__ import annotations

import argparse
import dataclasses
import logging
import math
import os
from abc import ABC, abstractmethod
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass
from enum import StrEnum
from pathlib import Path
from types import MappingProxyType

import libsbml
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import roadrunner
from sbmlsim import display, sensitivity
from sbmlsim.result import ScanResult
from sbmlsim.sensitivity import SensitivityResult
from sbmlsim.sensitivity.uncertainty import plot_distribution
from sbmlsim.simulation import Dimension, Observable, Scan, Simulation, sampling
from sbmlsim.simulation.sampling import Uniform
from sbmlsim.simulator import Simulator
from sbmlutils.console import console
from sbmlutils.report.units import udef_to_string

from pkpdsim.model import PKPDModel, SimulatorSettings
from pkpdsim.run import create_simulator

logger = logging.getLogger(__name__)

#: the columns of the table of the parameters, in the order of sbmlsim 0.8
COLUMNS: tuple[str, ...] = (
    "uid",
    "name",
    "value",
    "lower_bound",
    "upper_bound",
    "unit",
    "type",
    "reference",
)


class ParameterType(StrEnum):
    """Where the value and the bounds of a parameter of an analysis come from."""

    DATA = "data"
    SCALING = "scaling"
    NA = "na"
    FIT = "fitted"


@dataclass(kw_only=True)
class SensitivityParameter:
    """A parameter of a sensitivity analysis.

    Attributes:
        uid: the id in the model.
        name: the name in the model, the id if it has none.
        value: the value in the model after its initialization.
        lower_bound: the lower bound of the analysis.
        upper_bound: the upper bound of the analysis.
        unit: the unit in the model.
        type: where the value and the bounds come from.
        reference: the references of the value, e.g. a LaTeX citation.
    """

    uid: str
    name: str
    value: float = math.nan
    lower_bound: float = math.nan
    upper_bound: float = math.nan
    unit: str | None = None
    type: ParameterType = ParameterType.NA
    reference: str = ""


def parameters_from_sbml(
    sbml_path: Path,
    exclude_ids: Iterable[str] | None = None,
    exclude_na: bool = True,
    exclude_zero: bool = True,
) -> list[SensitivityParameter]:
    """Get the parameters of a model which an analysis varies, as sbmlsim 0.8 did.

    The constant parameters, the constant compartments and the constant or
    boundary species, in the order of the model, with their values after the
    initialization (the concentration of a species with only substance units).
    A parameter whose value is not a number (`exclude_na`) or is close to zero
    (`exclude_zero`), and every id of `exclude_ids`, is left out.

    Args:
        sbml_path: the SBML model.
        exclude_ids: the ids to leave out.
        exclude_na: leave out a parameter without a value.
        exclude_zero: leave out a parameter whose value is close to zero.

    Returns:
        The parameters, without bounds.
    """
    r = roadrunner.RoadRunner(str(sbml_path))
    model: libsbml.Model = libsbml.readSBMLFromFile(str(sbml_path)).getModel()
    excluded = set(exclude_ids or ())
    parameters: list[SensitivityParameter] = []

    def parameter(sbase: libsbml.SBase) -> SensitivityParameter:
        uid = sbase.getId()
        selection = uid
        if isinstance(sbase, libsbml.Species) and sbase.getHasOnlySubstanceUnits():
            selection = f"[{uid}]"
        return SensitivityParameter(
            uid=uid,
            name=sbase.getName() if sbase.isSetName() else uid,
            value=float(r.getValue(selection)),
            unit=udef_to_string(
                sbase.getDerivedUnitDefinition(), model=None, format="str"
            ),
        )

    for p in model.getListOfParameters():
        if p.getConstant():
            if exclude_na and np.isnan(p.getValue()):
                excluded.add(p.getId())
            if exclude_zero and np.isclose(r.getValue(p.getId()), 0.0):
                excluded.add(p.getId())
            parameters.append(parameter(p))
    for c in model.getListOfCompartments():
        if c.getConstant():
            if exclude_na and np.isnan(c.getSize()):
                excluded.add(c.getId())
            if exclude_zero and np.isclose(r.getValue(c.getId()), 0.0):
                excluded.add(c.getId())
            parameters.append(parameter(c))
    for s in model.getListOfSpecies():
        amount, concentration = s.isSetInitialAmount(), s.isSetInitialConcentration()
        if exclude_na and (
            (not amount and not concentration)
            or (amount and np.isnan(s.getInitialAmount()))
            or (concentration and np.isnan(s.getInitialConcentration()))
        ):
            excluded.add(s.getId())
        if exclude_zero and (
            (amount and np.isclose(s.getInitialAmount(), 0.0))
            or (concentration and np.isclose(s.getInitialConcentration(), 0.0))
        ):
            excluded.add(s.getId())
        if s.getConstant() or s.getBoundaryCondition():
            parameters.append(parameter(s))
    logger.info(
        "Excluded parameters: %s",
        sorted(p.uid for p in parameters if p.uid in excluded),
    )
    return [p for p in parameters if p.uid not in excluded]


def set_bounds(
    parameters: Iterable[SensitivityParameter],
    bounds: Iterable[tuple[str, float, float, ParameterType]],
) -> None:
    """Set the bounds and the type of parameters.

    Args:
        parameters: the parameters, changed in place.
        bounds: `(uid, lower bound, upper bound, type)`; a uid which is no
            parameter is logged and skipped, as sbmlsim 0.8 printed it.
    """
    by_uid = {p.uid: p for p in parameters}
    for uid, lower, upper, kind in bounds:
        if uid not in by_uid:
            logger.warning("unused bounds definition: %s = [%s, %s]", uid, lower, upper)
            continue
        p = by_uid[uid]
        p.lower_bound, p.upper_bound, p.type = lower, upper, kind


def set_relative_bounds(
    parameters: Iterable[SensitivityParameter], fraction: float
) -> None:
    """Let every parameter without bounds vary by `value * (1 +- fraction)`.

    The bounds of a negative value are reversed, the lower bound
    `value * (1 - fraction)` is above the upper bound `value * (1 + fraction)`,
    as in the table of sbmlsim 0.8, which stays the same; `distributions`
    refuses such a parameter, naming it, so a model sets its bounds with
    `set_bounds`. A value of 0 has the bounds 0 and 0.

    Args:
        parameters: the parameters, changed in place.
        fraction: the relative half width.
    """
    for p in parameters:
        if math.isnan(p.lower_bound) and math.isnan(p.upper_bound):
            p.lower_bound = p.value * (1 - fraction)
            p.upper_bound = p.value * (1 + fraction)


def parameters_to_df(
    parameters: Sequence[SensitivityParameter], sort: bool = True
) -> pd.DataFrame:
    """Get the table of the parameters.

    Args:
        parameters: the parameters.
        sort: sort the rows by the type and the id, as sbmlsim 0.8 did.

    Returns:
        One row per parameter, the columns `COLUMNS`, the type as its value.
    """
    rows = [{**asdict(p), "type": str(p.type)} for p in parameters]
    df = pd.DataFrame(rows, columns=list(COLUMNS))
    if sort:
        df = df.sort_values(by=["type", "uid"], ascending=True, ignore_index=True)
    return df


def write_parameters(
    parameters: Sequence[SensitivityParameter], directory: Path
) -> pd.DataFrame:
    """Write the table of the parameters as `parameters.tsv` and `parameters.tex`.

    Args:
        parameters: the parameters.
        directory: the directory, created if it does not exist.

    Returns:
        The table.
    """
    directory.mkdir(parents=True, exist_ok=True)
    df = parameters_to_df(parameters)
    df.to_csv(directory / "parameters.tsv", sep="\t", index=False)
    tex = df.to_latex(None, index=False, float_format="{:.3g}".format)
    (directory / "parameters.tex").write_text(tex.replace("_", r"\_"), encoding="utf-8")
    return df


def distributions(
    parameters: Sequence[SensitivityParameter],
) -> dict[str, Uniform]:
    """Get the distributions of the parameters of a global analysis, uniform in their bounds.

    Args:
        parameters: the parameters, with bounds.

    Returns:
        The distribution of every parameter, by its id.

    Raises:
        ValueError: naming the parameters without bounds, or with a lower bound
            above the upper one (the relative bounds of a negative value, see
            `set_relative_bounds`).
    """
    missing = sorted(
        p.uid
        for p in parameters
        if math.isnan(p.lower_bound) or math.isnan(p.upper_bound)
    )
    if missing:
        raise ValueError(f"The parameters {missing} have no bounds.")
    reversed_bounds = [
        f"'{p.uid}' ({p.lower_bound:g} > {p.upper_bound:g})"
        for p in sorted(parameters, key=lambda p: p.uid)
        if p.lower_bound > p.upper_bound
    ]
    if reversed_bounds:
        raise ValueError(
            f"The lower bound of the parameters {', '.join(reversed_bounds)} is "
            "above the upper one, e.g. the relative bounds of a negative value: "
            "set their bounds with `set_bounds`."
        )
    return {
        p.uid: Uniform(lower=p.lower_bound, upper=p.upper_bound) for p in parameters
    }


def free_parameters(
    parameters: Sequence[SensitivityParameter], simulations: Iterable[Simulation]
) -> list[SensitivityParameter]:
    """Get the parameters an analysis varies: those which none of its simulations sets.

    A design replaces its targets wherever a simulation sets them, so a
    parameter which the simulation of a group sets (in its `preinit_changes`, the
    `preinit_changes` of its presimulation or in a `Change`, e.g. `f_renal_function` of the degree of renal impairment)
    would take the values of the design in every group, and the groups would
    be one analysis. sbmlsim 0.8.5 left such a parameter out of an analysis,
    sbmlsim 0.9's designs do not, so a script takes the parameters of a design
    from here. The ids are compared as written, plain ids as sbmlsim 0.8.5
    did: a target given as a selection (`[S]`) does not match the parameter `S`.
    Every parameter left out is logged ("Removing fixed parameter"),
    as sbmlsim 0.8.5 printed it.

    Args:
        parameters: the parameters of the analysis.
        simulations: the simulations of the groups of the analysis.

    Returns:
        The parameters which no simulation sets, in their order.
    """
    fixed: set[str] = set()
    for simulation in simulations:
        fixed.update(simulation.preinit_changes)
        if simulation.presimulation is not None:
            fixed.update(simulation.presimulation.preinit_changes)
        for change in simulation.changes:
            fixed.update(change.values)
    for p in parameters:
        if p.uid in fixed:
            logger.info("Removing fixed parameter: %s", p.uid)
    return [p for p in parameters if p.uid not in fixed]


#: the environment variable which replaces the number of samples of a global
#: analysis, for a smoke run of the sensitivity scripts
SAMPLE_SIZE = "PKPDSIM_SENSITIVITY_N"


def sample_size(n: int) -> int:
    """Get the number of samples of a global analysis, `n` unless a smoke run sets another.

    The environment variable `PKPDSIM_SENSITIVITY_N` (a positive int) replaces
    `n`, e.g. `PKPDSIM_SENSITIVITY_N=16` for a quick run of the scripts in a
    verification.

    Args:
        n: the number of samples of the analysis.

    Returns:
        The number of samples.

    Raises:
        ValueError: if the variable is set and is no positive int.
    """
    value = os.environ.get(SAMPLE_SIZE, "")
    if not value:
        return n
    try:
        size = int(value)
    except ValueError:
        size = 0
    if size < 1:
        raise ValueError(
            f"{SAMPLE_SIZE} is a positive number of samples, not '{value}'."
        )
    return size


#: the seed of the designs of the analyses
SEED = 1234


@dataclass(frozen=True)
class SensitivityOutput:
    """An output of the analyses, an observable of the simulation of a group.

    Attributes:
        uid: the id of the output in the tables and figures, `<uid>.<index>.tsv`.
        name: the name of the output in the titles of the figures and in the
            summary of `sampling`.
        observable: the id of the observable in the result of the scan, the uid
            when empty; an id of an output which is no identifier, e.g.
            `[Cve_api]_cmax`, names an observable `Cve_api_cmax`.
    """

    uid: str
    name: str
    observable: str = ""

    def __post_init__(self) -> None:
        """Let the output be its own observable when it names none."""
        if not self.observable:
            object.__setattr__(self, "observable", self.uid)


@dataclass(frozen=True)
class SensitivityGroup:
    """A group of the analyses, e.g. a degree of renal impairment.

    Attributes:
        uid: the id of the group, the folder of its tables and figures.
        name: the name of the group in the titles of the figures.
        simulation: the simulation of the group, its changes in
            `preinit_changes`; an analysis varies no parameter which it sets,
            see `free_parameters`.
    """

    uid: str
    name: str
    simulation: Simulation


@dataclass(frozen=True)
class _Context:
    """What a run hands to its analyses."""

    simulator: Simulator
    model_path: Path
    parameters: Sequence[SensitivityParameter]
    outputs: Sequence[SensitivityOutput]
    observables: Sequence[Observable]
    seed: int

    def scan(self, group: SensitivityGroup, design: Dimension) -> ScanResult:
        """Run a design over the simulation of a group, the outputs per point."""
        return self.simulator.run(
            self.model_path,
            Scan(group.simulation, [design]),
            self.observables,
            keep=[o.observable for o in self.outputs],
        )

    def write_indices(
        self, result: SensitivityResult, indices: Sequence[str], path: Path
    ) -> None:
        """Write every index of every output as a table, `<output>.<index>.tsv`."""
        for output in self.outputs:
            for index in indices:
                column = f"{output.uid}.{index}"
                result.to_dataframe(f"{output.observable}.{index}").rename(
                    columns={f"{output.observable}.{index}": column}
                ).to_csv(path / f"{column}.tsv", sep="\t")


def _repeated(values: Iterable[str]) -> list[str]:
    """Get the values which are given more than once, sorted."""
    seen: set[str] = set()
    repeated: set[str] = set()
    for value in values:
        (repeated if value in seen else seen).add(value)
    return sorted(repeated)


def _check_samples(n: int) -> None:
    """Check a number of samples of an analysis.

    Raises:
        ValueError: if it is no positive int.
    """
    if n < 1:
        raise ValueError(f"The number of samples is a positive int, not {n}.")


@dataclass(frozen=True, kw_only=True)
class Analysis(ABC):
    """An analysis of a script over the simulations of its groups.

    The parameters it varies are the parameters of the script which no
    simulation of its groups sets (`free_parameters`); its tables and figures
    are written per group, into `<analysis>/<group>`.

    Attributes:
        groups: the groups of the analysis, at least one, the ids unique;
            stored as a tuple.
    """

    groups: Sequence[SensitivityGroup]

    def __post_init__(self) -> None:
        """Check the groups.

        Raises:
            ValueError: for no group or two groups of one id.
        """
        object.__setattr__(self, "groups", tuple(self.groups))
        if not self.groups:
            raise ValueError(f"{type(self).__name__} needs at least one group.")
        repeated = _repeated(g.uid for g in self.groups)
        if repeated:
            raise ValueError(
                f"The groups {repeated} of {type(self).__name__} are given twice."
            )

    def varied(
        self, parameters: Sequence[SensitivityParameter]
    ) -> list[SensitivityParameter]:
        """Get the parameters the analysis varies, which no group sets.

        Args:
            parameters: the parameters of the script.

        Returns:
            The parameters, see `free_parameters`.
        """
        return free_parameters(parameters, [g.simulation for g in self.groups])

    @staticmethod
    def _group_path(path: Path, group: SensitivityGroup) -> Path:
        """Get the folder of a group, created if it does not exist."""
        group_path = path / group.uid
        group_path.mkdir(parents=True, exist_ok=True)
        return group_path

    @abstractmethod
    def _run(self, context: _Context, path: Path) -> None:
        """Run the analysis and write its tables and figures.

        Args:
            context: the simulator, the model, the parameters, the outputs and
                observables and the seed of the run.
            path: the folder of the analysis.
        """


@dataclass(frozen=True, kw_only=True)
class LocalAnalysis(Analysis):
    """The local sensitivities of every group, each parameter by `+-delta`.

    The design `sampling.local` per group (its values are those of the
    simulation of the group), the indices `sensitivity.local`, the table of
    `normalized` per output and its heatmap, `normalized.png`.

    Attributes:
        delta: the relative change of a parameter.
    """

    delta: float = 0.01

    def __post_init__(self) -> None:
        """Check the groups and the change.

        Raises:
            ValueError: see `Analysis`, and for a change which is not positive.
        """
        super().__post_init__()
        if not self.delta > 0:
            raise ValueError(
                f"The delta of a local analysis is positive, not {self.delta}."
            )

    def _run(self, context: _Context, path: Path) -> None:
        """Run the analysis, see the class."""
        targets = [p.uid for p in self.varied(context.parameters)]
        for group in self.groups:
            group_path = self._group_path(path, group)
            design = sampling.local(
                targets,
                delta=self.delta,
                model=context.model_path,
                simulation=group.simulation,
            )
            result = sensitivity.local(context.scan(group, design))
            context.write_indices(result, ["normalized"], group_path)
            sensitivity.plot_heatmap(
                result,
                "normalized",
                title=f"Local sensitivity ({group.name})",
                path=group_path / "normalized.png",
            )
            plt.close("all")


@dataclass(frozen=True, kw_only=True)
class SamplingAnalysis(Analysis):
    """The uncertainty of the outputs over a Latin hypercube of the parameters.

    One design `sampling.lhs` of `sample_size(n)` points for every group;
    per group the outputs of every point, `samples.tsv`, their statistics,
    `summary.tsv` (a row per output: its id, name, the statistics of
    `ScanResult.summary` with the quantiles, the unit), and their
    distributions, `distributions.png`.

    Attributes:
        n: the number of samples, see `sample_size`.
        quantiles: the quantiles of the summary.
    """

    n: int = 1000
    quantiles: Sequence[float] = (0.025, 0.05, 0.5, 0.95, 0.975)

    def __post_init__(self) -> None:
        """Check the groups and the number of samples.

        Raises:
            ValueError: see `Analysis`, and for a number of samples which is
                not positive.
        """
        super().__post_init__()
        _check_samples(self.n)
        object.__setattr__(self, "quantiles", tuple(self.quantiles))

    def _run(self, context: _Context, path: Path) -> None:
        """Run the analysis, see the class."""
        design = sampling.lhs(
            distributions(self.varied(context.parameters)),
            n=sample_size(self.n),
            seed=context.seed,
        )
        observables = {o.observable: o.uid for o in context.outputs}
        for group in self.groups:
            group_path = self._group_path(path, group)
            res = context.scan(group, design)
            res.ds[list(observables)].to_dataframe().rename(columns=observables).to_csv(
                group_path / "samples.tsv", sep="\t"
            )
            summary = res.summary(design.id, quantiles=list(self.quantiles))
            rows = []
            for output in context.outputs:
                values = summary[output.observable]
                row: dict[str, object] = {"output": output.uid, "name": output.name}
                for statistic, value in zip(
                    values["statistic"].values, values.values, strict=True
                ):
                    row[str(statistic)] = float(value)
                row["unit"] = res.units[output.observable]
                rows.append(row)
            pd.DataFrame(rows).to_csv(group_path / "summary.tsv", sep="\t", index=False)
            self._plot(context.outputs, res, design, group, group_path)

    @staticmethod
    def _plot(
        outputs: Sequence[SensitivityOutput],
        res: ScanResult,
        design: Dimension,
        group: SensitivityGroup,
        path: Path,
    ) -> None:
        """Plot the distribution of every output, four in a row."""
        n_cols = 4
        n_rows = -(-len(outputs) // n_cols)
        fig, axes = plt.subplots(
            nrows=n_rows,
            ncols=n_cols,
            figsize=(5 * n_cols, 4 * n_rows),
            layout="constrained",
            squeeze=False,
        )
        for ax, output in zip(axes.flatten(), outputs, strict=False):
            plot_distribution(res, output.observable, dim=design.id, ax=ax)
            ax.set_title(output.name)
            ax.locator_params(axis="x", nbins=4)
        for ax in axes.flatten()[len(outputs) :]:
            ax.set_visible(False)
        fig.suptitle(f"Uncertainty ({group.name})")
        fig.savefig(path / "distributions.png", dpi=300, bbox_inches="tight")
        plt.close("all")


@dataclass(frozen=True, kw_only=True)
class _GlobalAnalysis(Analysis):
    """A variance based analysis: the indices `S1` and `ST` of every group.

    One design for every group; per group the tables of `S1` and `ST` per
    output and the figure of its indices, `<output>.png`.

    Attributes:
        n: the base number of samples, see `sample_size`.
    """

    n: int

    def __post_init__(self) -> None:
        """Check the groups and the number of samples.

        Raises:
            ValueError: see `Analysis`, and for a number of samples which is
                not positive.
        """
        super().__post_init__()
        _check_samples(self.n)

    @abstractmethod
    def _design(
        self, parameters: Sequence[SensitivityParameter], seed: int
    ) -> Dimension:
        """Get the design over the parameters the analysis varies."""

    @abstractmethod
    def _indices(self, result: ScanResult) -> SensitivityResult:
        """Get the indices of the result of the design."""

    def _run(self, context: _Context, path: Path) -> None:
        """Run the analysis, see the class."""
        design = self._design(self.varied(context.parameters), context.seed)
        for group in self.groups:
            group_path = self._group_path(path, group)
            result = self._indices(context.scan(group, design))
            context.write_indices(result, ["S1", "ST"], group_path)
            for output in context.outputs:
                sensitivity.plot_indices(
                    result,
                    output.observable,
                    title=f"{output.name} ({group.name})",
                    path=group_path / f"{output.uid}.png",
                )
                plt.close("all")


@dataclass(frozen=True, kw_only=True)
class FastAnalysis(_GlobalAnalysis):
    """The indices of the extended FAST, `sampling.fast` and `sensitivity.fast`.

    Attributes:
        n: the number of samples per parameter, more than 64, see `sample_size`.
    """

    n: int = 1000

    def _design(
        self, parameters: Sequence[SensitivityParameter], seed: int
    ) -> Dimension:
        """Get the design over the parameters, see `sampling.fast`."""
        return sampling.fast(
            distributions(parameters), n=sample_size(self.n), seed=seed
        )

    def _indices(self, result: ScanResult) -> SensitivityResult:
        """Get the indices, see `sensitivity.fast`."""
        return sensitivity.fast(result)


@dataclass(frozen=True, kw_only=True)
class SobolAnalysis(_GlobalAnalysis):
    """The Sobol indices, `sampling.sobol` and `sensitivity.sobol`.

    Attributes:
        n: the base number of samples, a power of 2, see `sample_size`.
    """

    n: int = 4096

    def _design(
        self, parameters: Sequence[SensitivityParameter], seed: int
    ) -> Dimension:
        """Get the design over the parameters, see `sampling.sobol`."""
        return sampling.sobol(
            distributions(parameters), n=sample_size(self.n), seed=seed
        )

    def _indices(self, result: ScanResult) -> SensitivityResult:
        """Get the indices, see `sensitivity.sobol`."""
        return sensitivity.sobol(result)


@dataclass(frozen=True, kw_only=True)
class SensitivityRunner:
    """The sensitivity analyses of a model, run by its script.

    The script gives what is the model's own (the parameters, the outputs with
    their observables, the analyses with their groups) and runs them with `main`
    from its command line, or with `run`. A run writes the table of the
    parameters (`write_parameters`) into its folder,
    `model.results_dir / "sensitivity"`, and the tables and figures of an
    analysis into `<analysis>/<group>` below it.

    Attributes:
        model: the model object: the model which is simulated, the folder of
            the results and the settings of the simulator.
        parameters: the function which gives the parameters of the analyses
            with their bounds; it reads the model, so it is called by a run.
        outputs: the outputs of the analyses, at least one, the ids and the
            observables unique; stored as a tuple.
        observables: the observables of the simulation of a group, which give
            the observables of the outputs.
        analyses: the analyses by name, the name on the command line and the
            folder; at least one, stored read-only.
        default: the analyses which run without a name on the command line.
        seed: the seed of the designs.
    """

    model: PKPDModel
    parameters: Callable[[], Sequence[SensitivityParameter]]
    outputs: Sequence[SensitivityOutput]
    observables: Sequence[Observable]
    analyses: Mapping[str, Analysis]
    default: Sequence[str] = ("local",)
    seed: int = SEED

    def __post_init__(self) -> None:
        """Check the outputs and the analyses.

        Raises:
            ValueError: for no output, two outputs of one id or one
                observable, an output whose observable is none of the
                observables, no analysis, or a default analysis which is none
                of the analyses.
        """
        object.__setattr__(self, "outputs", tuple(self.outputs))
        object.__setattr__(self, "observables", tuple(self.observables))
        object.__setattr__(self, "analyses", MappingProxyType(dict(self.analyses)))
        object.__setattr__(self, "default", tuple(self.default))
        if not self.outputs:
            raise ValueError("A sensitivity analysis needs at least one output.")
        repeated = _repeated(o.uid for o in self.outputs)
        if repeated:
            raise ValueError(f"The outputs {repeated} are given twice.")
        repeated = _repeated(o.observable for o in self.outputs)
        if repeated:
            raise ValueError(f"The observables {repeated} are outputs twice.")
        ids = {observable.id for observable in self.observables}
        missing = [o.observable for o in self.outputs if o.observable not in ids]
        if missing:
            raise ValueError(f"The outputs name no observables {missing}.")
        if not self.analyses:
            raise ValueError("A sensitivity analysis needs at least one analysis.")
        unknown = [name for name in self.default if name not in self.analyses]
        if unknown:
            raise ValueError(
                f"The default analyses {unknown} are none of {list(self.analyses)}."
            )

    def run(
        self,
        names: Sequence[str] | None = None,
        *,
        path: Path | None = None,
        simulator: SimulatorSettings | None = None,
    ) -> Path:
        """Run analyses and write the table of the parameters.

        Args:
            names: the names of the analyses, `default` if `None`.
            path: the folder of the results, `model.results_dir / "sensitivity"`
                if `None`; e.g. the folder of one of several simulations.
            simulator: the settings of the simulator, the settings of the model
                with `n_workers=None` (sbmlsim chooses the processes of the
                pool, the script has the guard `if __name__ == "__main__":`) if
                `None`.

        Returns:
            The folder of the results.

        Raises:
            ValueError: naming the names which are no analysis, before anything
                runs.
        """
        names = self.default if names is None else tuple(names)
        unknown = [name for name in names if name not in self.analyses]
        if unknown:
            raise ValueError(f"Unknown analyses {unknown}, use {list(self.analyses)}.")
        path = self.model.results_dir / "sensitivity" if path is None else path
        path.mkdir(parents=True, exist_ok=True)
        parameters = list(self.parameters())
        df_parameters = write_parameters(parameters, path)
        with pd.option_context("display.float_format", "{:.5g}".format):
            console.print(df_parameters)

        settings = (
            dataclasses.replace(self.model.simulator, n_workers=None)
            if simulator is None
            else simulator
        )
        context = _Context(
            simulator=create_simulator(settings),
            model_path=self.model.model_path,
            parameters=parameters,
            outputs=self.outputs,
            observables=self.observables,
            seed=self.seed,
        )
        for name in names:
            self.analyses[name]._run(context, path / name)

        display.section("Report", icon=display.ICON_REPORT)
        display.key_values({"analyses": ", ".join(names)})
        display.link("results", path)
        return path

    def main(self, argv: Sequence[str] | None = None) -> None:
        """Run the analyses named on the command line, see `run`.

        Args:
            argv: the arguments, those of the process if `None`: the names of
                the analyses, `default` without one; an unknown name exits with
                2.
        """
        parser = argparse.ArgumentParser(
            description=(
                f"The sensitivity analyses of the {self.model.title}, into "
                f"{self.model.results_dir / 'sensitivity'}."
            ),
            epilog=(
                f"The environment variable {SAMPLE_SIZE} replaces the numbers of "
                "samples of the sampling and the global analyses for a smoke run."
            ),
        )
        # the names are checked here: argparse of python 3.13 checks the
        # default of a positional `nargs="*"` as one choice
        parser.add_argument(
            "analyses",
            nargs="*",
            metavar="ANALYSIS",
            help=(
                f"the analyses, of {', '.join(self.analyses)} "
                f"(default: {' '.join(self.default)})"
            ),
        )
        names: list[str] = parser.parse_args(argv).analyses
        unknown = [name for name in names if name not in self.analyses]
        if unknown:
            parser.error(
                f"invalid choice: {', '.join(repr(name) for name in unknown)} "
                f"(choose from {', '.join(self.analyses)})"
            )
        self.run(tuple(names) or self.default)
