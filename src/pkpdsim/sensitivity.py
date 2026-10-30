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
from collections.abc import Callable, Collection, Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass
from enum import StrEnum
from pathlib import Path
from types import MappingProxyType
from typing import NamedTuple, override
from xml.etree import ElementTree

import libsbml
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from sbmlsim import display, sensitivity
from sbmlsim.model.model_roadrunner import JitOptimization, load_roadrunner
from sbmlsim.result import ScanResult
from sbmlsim.sensitivity import SensitivityResult
from sbmlsim.sensitivity.classification import (
    sensitivity_classification,
    sensitivity_classification_symbol,
    uncertainty_classification,
    uncertainty_classification_symbol,
)
from sbmlsim.sensitivity.result import PARAMETER
from sbmlsim.simulation import Dimension, Observable, Scan, Simulation, sampling
from sbmlsim.simulation.sampling import Uniform
from sbmlsim.simulator import Simulator
from sbmlutils.console import console
from sbmlutils.report.units import udef_to_string

from pkpdsim.model import PKPDModel, SimulatorSettings
from pkpdsim.run import create_simulator
from pkpdsim.sex import Sex

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
    *,
    exclude_derived: bool = False,
) -> list[SensitivityParameter]:
    """Get the parameters of a model which an analysis varies, as sbmlsim 0.8 did.

    The constant parameters, the constant compartments and the constant or
    boundary species, in the order of the model, with their values after the
    initialization (the concentration of a species with only substance units).
    A parameter whose value is not a number (`exclude_na`) or is close to zero
    (`exclude_zero`), and every id of `exclude_ids`, is left out; so is every
    derived parameter, compartment and species (`exclude_derived`).

    Args:
        sbml_path: the SBML model.
        exclude_ids: the ids to leave out.
        exclude_na: leave out a parameter without a value.
        exclude_zero: leave out a parameter whose value is close to zero.
        exclude_derived: leave out every parameter, compartment and species
            with an initial assignment, e.g. the organ volumes and the
            clearances which a library model derives from its covariates; the
            covariates and the fitted parameters have none.

    Returns:
        The parameters, without bounds.
    """
    # the values only: compiled with LLJIT at O0, which compiles fastest and
    # gives the values of roadrunner's default MCJIT bitwise (the creatinine
    # model of pkdb_models in 4.5 s instead of 48 s)
    r = load_roadrunner(sbml_path, JitOptimization.O0)
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
    if exclude_derived:
        excluded.update(
            p.uid for p in parameters if model.getInitialAssignmentBySymbol(p.uid)
        )
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
    parameters: Sequence[SensitivityParameter],
    sort: bool = True,
    builds: Mapping[Sex, Mapping[str, float]] | None = None,
) -> pd.DataFrame:
    """Get the table of the parameters.

    Args:
        parameters: the parameters.
        sort: sort the rows by the type and the id, as sbmlsim 0.8 did.
        builds: the values of the parameters in other builds of the model than
            the one of the table, by sex (`SensitivityRunner` gives those of the
            sexes of its groups other than male); each is the column
            `value_<sex>` after `value`.

    Returns:
        One row per parameter, the columns `COLUMNS` (and `value_<sex>`), the
        type as its value.
    """
    builds = builds or {}
    columns = list(COLUMNS)
    position = columns.index("value") + 1
    columns[position:position] = [f"value_{sex}" for sex in builds]
    rows = [
        {
            **asdict(p),
            "type": str(p.type),
            **{
                f"value_{sex}": values.get(p.uid, math.nan)
                for sex, values in builds.items()
            },
        }
        for p in parameters
    ]
    df = pd.DataFrame(rows, columns=columns)
    if sort:
        df = df.sort_values(by=["type", "uid"], ascending=True, ignore_index=True)
    return df


def write_parameters(
    parameters: Sequence[SensitivityParameter],
    directory: Path,
    builds: Mapping[Sex, Mapping[str, float]] | None = None,
) -> pd.DataFrame:
    """Write the table of the parameters as `parameters.tsv` and `parameters.tex`.

    Args:
        parameters: the parameters.
        directory: the directory, created if it does not exist.
        builds: the values in the other builds of the model, see
            `parameters_to_df`.

    Returns:
        The table.
    """
    directory.mkdir(parents=True, exist_ok=True)
    df = parameters_to_df(parameters, builds=builds)
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


def build_values(sbml_path: Path, uids: Iterable[str]) -> dict[str, float]:
    """Get the values of parameters in a build of a model, as `parameters_from_sbml` reads them.

    Args:
        sbml_path: the SBML model, e.g. the female build.
        uids: the ids of the parameters, compartments and species.

    Returns:
        The value after the initialization of every id the model has.
    """
    wanted = set(uids)
    return {
        p.uid: p.value
        for p in parameters_from_sbml(sbml_path, exclude_na=False, exclude_zero=False)
        if p.uid in wanted
    }


def design_factors(
    parameters: Sequence[SensitivityParameter], values: Mapping[str, float]
) -> dict[str, float]:
    """Get the factors which carry a design of the table to the values of a build.

    A global design draws the values of the parameters in their bounds around
    the values of the table, which are those of the male build of the model. A
    group of another build samples the same design relative to the values of
    its build: every value of a parameter is multiplied with the ratio of its
    value in the build to its value in the table, so a group of women varies
    the body weight around the female body weight by the same factors as the
    men around the male one. The bounds of the table are taken relative to its
    values.

    Args:
        parameters: the parameters which a design varies.
        values: the values of the parameters in the build of a group.

    Returns:
        The ratio of every parameter whose value differs, by its id.

    Raises:
        ValueError: for a parameter which has no value in the build, and for
            a parameter whose value differs and whose bounds cannot be taken
            relative to its value: a value of the table which is zero or not
            finite, or a ratio which is not positive (another sign); naming the
            parameter and both values.
    """
    missing = [p.uid for p in parameters if p.uid not in values]
    if missing:
        raise ValueError(
            f"The varied parameter {', '.join(repr(m) for m in missing)} has no "
            "value in the build, so the design has no counterpart in it."
        )
    factors: dict[str, float] = {}
    problems: list[str] = []
    for p in parameters:
        value = values[p.uid]
        if value == p.value:
            continue
        ratio = value / p.value if p.value != 0 else math.nan
        if not (math.isfinite(p.value) and math.isfinite(ratio) and ratio > 0):
            problems.append(
                f"'{p.uid}' ({p.value:g} in the table, {value:g} in the build)"
            )
            continue
        factors[p.uid] = ratio
    if problems:
        raise ValueError(
            f"The values of the parameters {', '.join(problems)} differ between "
            "the table and the build, and their bounds cannot be taken relative "
            "to the value of the table, so the design has no counterpart in the "
            "build: give them a value in the table which is neither zero nor of "
            "another sign, or leave them out of the analysis."
        )
    return factors


def relative_design(design: Dimension, factors: Mapping[str, float]) -> Dimension:
    """Get a design with the values of its parameters multiplied by factors.

    The design of a group of another build, see `design_factors`; the bounds
    of the uniform distributions of its record are multiplied as well.

    Args:
        design: the design over the values of the table.
        factors: the factor of a parameter, by its id; 1 for the others.

    Returns:
        The design of the group, the given one without a factor.
    """
    factors = {uid: f for uid, f in factors.items() if uid in design.values}
    if not factors:
        return design
    record = design.design
    if record is not None:
        distributions = {}
        for uid, d in record.distributions.items():
            if uid in factors:
                d = {
                    **d,
                    **{
                        key: d[key] * factors[uid]
                        for key in ("lower", "upper")
                        if isinstance(d.get(key), int | float)
                    },
                }
            distributions[uid] = d
        record = dataclasses.replace(record, distributions=distributions)
    return Dimension(
        design.id,
        values={
            uid: np.asarray(v) * factors.get(uid, 1.0)
            for uid, v in design.values.items()
        },
        labels=design.labels,
        design=record,
    )


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

#: the cutoff of the heatmaps of sbmlsim 0.8.5
CUTOFF_085 = 0.05


@dataclass(frozen=True)
class SensitivityOutput:
    """An output of the analyses, an observable of the simulation of a group.

    Attributes:
        uid: the id of the output in the tables and the names of the files:
            the columns `output` and `uid` of the tables of sbmlsim 0.8.5, the
            columns of `samples.tsv`, `<uid>.S1.tsv`, `<uid>.ST.tsv` and the
            figure of its indices, `..._<o>_<uid>.png`.
        name: the name of the output in the titles of the figures, in the
            heatmaps and in the column `name` of the tables of sbmlsim 0.8.5,
            unique among the outputs of a runner.
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
        color: the color of the group in the figures of `sampling` which
            compare the groups, a color of matplotlib; the color of the
            property cycle of matplotlib (`C0`, `C1`, ...) of its position if
            `None`.
        sex: the sex of the group, whose build of the model it runs
            (`PKPDModel.model_paths`), the male build by default; the design
            of an analysis is the same for every group.
    """

    uid: str
    name: str
    simulation: Simulation
    color: str | None = None
    sex: Sex = Sex.MALE


@dataclass(frozen=True)
class _Context:
    """What a run hands to its analyses."""

    simulator: Simulator
    model: PKPDModel
    outputs: Sequence[SensitivityOutput]
    observables: Sequence[Observable]
    seed: int
    #: the values of the parameters in the builds of the groups other than
    #: the male one, which the table has, by sex
    builds: Mapping[Sex, Mapping[str, float]] = MappingProxyType({})

    def model_path(self, group: SensitivityGroup) -> Path:
        """Get the build of the model which a group runs, the build of its sex."""
        return self.model.model_path_of(group.sex)

    def design(
        self,
        group: SensitivityGroup,
        design: Dimension,
        varied: Sequence[SensitivityParameter],
    ) -> Dimension:
        """Get the global design of a group, around the values of its build.

        The design of the table for a male group, the same relative design
        around the values of the build of the group otherwise
        (`design_factors`, `relative_design`).
        """
        values = self.builds.get(group.sex)
        if values is None:
            return design
        return relative_design(design, design_factors(varied, values))

    def scan(self, group: SensitivityGroup, design: Dimension) -> ScanResult:
        """Run a design over the simulation of a group, the outputs per point."""
        return self.simulator.run(
            self.model_path(group),
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

    def unit(self, result: ScanResult, output: SensitivityOutput) -> str:
        """Get the unit of an output: as its observable gives it, else of the result.

        The unit of a `Formula` or `Custom` observable is the one the script
        wrote (`mM*min`), as in the tables of sbmlsim 0.8.5; a `PK` observable
        has the unit of the result (`millimolar * minute`).
        """
        for observable in self.observables:
            if observable.id == output.observable:
                unit = getattr(observable, "unit", None)
                if unit:
                    return str(unit)
        return str(result.units.get(output.observable, ""))


def _repeated(values: Iterable[str]) -> list[str]:
    """Get the values which are given more than once, sorted."""
    seen: set[str] = set()
    repeated: set[str] = set()
    for value in values:
        (repeated if value in seen else seen).add(value)
    return sorted(repeated)


def _sequence_of_names(names: Sequence[str], what: str) -> tuple[str, ...]:
    """Get names as a tuple, refusing a single string, which is no sequence of names.

    Raises:
        TypeError: for a string, e.g. `run("local")` instead of `run(["local"])`.
    """
    if isinstance(names, str):
        raise TypeError(
            f"The {what} are a sequence of names, not the string {names!r}."
        )
    return tuple(names)


def _check_subfolder(subfolder: str) -> None:
    """Check that a subfolder lies below the folder of the results.

    Raises:
        ValueError: for an empty folder, one with a drive or a root (`/tmp`,
            which windows takes as relative to the drive) or with `..`.
    """
    if not subfolder or Path(subfolder).anchor or ".." in Path(subfolder).parts:
        raise ValueError(
            f"The subfolder {subfolder!r} is a relative folder below the folder "
            "of the results."
        )


def _reaches(result: SensitivityResult, index: str, cutoff: float) -> bool:
    """Check whether an index reaches a cutoff, as `plot_heatmap` keeps its rows.

    A row is kept if a value reaches the cutoff or one is undefined.
    """
    values = np.asarray(result.index(index).values, dtype=float)
    finite = np.isfinite(values)
    magnitude = np.abs(np.where(finite, values, 0.0))
    return bool((~finite).any() or (magnitude >= cutoff).any())


def _heatmap(
    result: SensitivityResult,
    index: str,
    *,
    cutoff: float,
    path: Path,
    **options: object,
) -> None:
    """Draw a heatmap of `plot_heatmap`, every parameter if none reaches the cutoff.

    `plot_heatmap` refuses an index which no parameter reaches, e.g. the local
    sensitivities of a group whose parameters do not move its outputs; the
    heatmap then shows every parameter, which a warning names.
    """
    if cutoff > 0 and not _reaches(result, index, cutoff):
        logger.warning(
            "No parameter reaches the cutoff %g of '%s' in '%s', the heatmap "
            "shows every parameter.",
            cutoff,
            index,
            path.name,
        )
        cutoff = 0.0
    sensitivity.plot_heatmap(result, index, cutoff=cutoff, path=path, **options)  # ty: ignore[invalid-argument-type]


def _labelled(
    result: SensitivityResult,
    parameters: Sequence[SensitivityParameter],
    outputs: Sequence[SensitivityOutput],
) -> SensitivityResult:
    """Get a result labelled as the heatmaps of sbmlsim 0.8.5.

    The parameters are `uid: name`, the observables the names of the outputs.
    """
    names = {o.observable: o.name for o in outputs}
    renames: dict[str, str] = {}
    for variable in result.ds.data_vars:
        observable, index = str(variable).rsplit(".", 1)
        if observable in names:
            renames[str(variable)] = f"{names[observable]}.{index}"
    labels = {p.uid: f"{p.uid}: {p.name}" for p in parameters}
    ds = result.ds.rename(renames).assign_coords(
        {PARAMETER: [labels.get(p, p) for p in result.parameters]}
    )
    ds.attrs = {
        **ds.attrs,
        "units": {renames.get(k, k): v for k, v in result.units.items()},
    }
    return SensitivityResult(ds)


def _classes(
    value: float,
    classify: Callable[[float], StrEnum | None],
    symbolize: Callable[[float], str],
) -> tuple[str, str]:
    """Get the symbol and the class of a value, empty if it is undefined.

    The classifications of sbmlsim give `None` for `NaN`, the sensitivity of
    an output which is zero at the reference or the uncertainty of one with the
    median zero, and classify a negative uncertainty by its magnitude.
    """
    classification = classify(value)
    if classification is None:
        return "", ""
    return symbolize(value), classification.value


def _figshape(n: int) -> tuple[int, int]:
    """Get the rows and columns of a figure of `n` panels, as sbmlsim 0.8.5 did."""
    if n <= 4:
        return 1, n
    ncols = np.ceil(np.sqrt(n))
    n_empty = ncols * ncols - n
    nrows = ncols - np.floor(n_empty / ncols)
    return int(nrows), int(ncols)


def _plot_groups(
    panels: Sequence[tuple[str, Sequence[np.ndarray]]],
    groups: Sequence[SensitivityGroup],
    path: Path,
    seed: int,
) -> None:
    """Draw the values of every group as box, violin and points, a panel per quantity.

    The figures of the sampling of sbmlsim 0.8.5: a box (without outliers), a
    violin (the mean and the median) and the points of every group side by
    side in the color of the group. A group with fewer than two distinct
    values has no violin, which needs a spread.

    Args:
        panels: the label of the axis and the values of every group.
        groups: the groups.
        path: the file of the figure.
        seed: the seed of the jitter of the points.
    """
    nrows, ncols = _figshape(len(panels))
    # a `Figure`, not a pyplot figure: pyplot binds it to the backend, which
    # aborts the process with TkAgg next to a pool of processes
    fig = Figure(figsize=(4 * ncols, 4 * nrows), layout="constrained")
    axes = fig.subplots(nrows=nrows, ncols=ncols, squeeze=False)
    colors = [g.color or f"C{k}" for k, g in enumerate(groups)]
    rng = np.random.default_rng(seed)
    for k, ax in enumerate(axes.flat):
        if k >= len(panels):
            ax.axis("off")
            continue
        label, values = panels[k]
        data = [v[np.isfinite(v)] for v in values]
        boxes = ax.boxplot(
            data,
            positions=range(len(groups)),
            tick_labels=[g.uid for g in groups],
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "black"},
            whiskerprops={"color": "black"},
            capprops={"color": "black"},
        )
        for box, color in zip(boxes["boxes"], colors, strict=True):
            box.set_facecolor(color)
        for kg, (d, color) in enumerate(zip(data, colors, strict=True)):
            if np.unique(d).size > 1:
                violin = ax.violinplot(
                    d,
                    positions=[kg + 0.3],
                    showmeans=True,
                    showmedians=True,
                    showextrema=False,
                )
                for body in violin["bodies"]:
                    body.set_facecolor(color)
                # the lines of the mean and the median in one color, as the one
                # call of `violinplot` of sbmlsim 0.8.5 drew them
                for key in ["cmeans", "cmedians"]:
                    violin[key].set_color("C0")
            ax.scatter(
                rng.normal(kg + 0.3, 0.02, d.size),
                d,
                alpha=0.7,
                s=30,
                color="white",
                edgecolors="black",
            )
        ax.set_ylabel(label, fontsize=13, fontweight="bold")
        for tick in ax.get_xticklabels() + ax.get_yticklabels():
            tick.set_fontweight("bold")
        ax.tick_params(axis="x", labelrotation=90)
    fig.savefig(path, dpi=300, bbox_inches="tight")


@dataclass(frozen=True, kw_only=True)
class Analysis(ABC):
    """An analysis of a script over the simulations of its groups.

    The parameters it varies are the parameters of the script which no
    simulation of its groups sets (`free_parameters`). It writes one result per
    quantity: the tables and figures of sbmlsim 0.8.5, with their names, into
    `<analysis>`, and what sbmlsim 0.8.5 wrote no file of (the samples, the
    tables of the global indices) per group into `<analysis>/<group>`. The analyses are the four of
    this module, `LocalAnalysis`, `SamplingAnalysis`, `FastAnalysis` and
    `SobolAnalysis`; a script with an analysis of its own runs it itself.

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

    def _check(
        self, parameters: Sequence[SensitivityParameter]
    ) -> list[SensitivityParameter]:
        """Check what the analysis needs, at the start of a run.

        Args:
            parameters: the parameters of the script.

        Returns:
            The parameters the analysis varies.

        Raises:
            ValueError: if it varies no parameter.
        """
        varied = self.varied(parameters)
        if not varied:
            raise ValueError(
                f"{type(self).__name__} varies no parameter: the simulations of "
                f"its groups {[g.uid for g in self.groups]} set every parameter."
            )
        return varied

    @staticmethod
    def _group_path(path: Path, group: SensitivityGroup) -> Path:
        """Get the folder of a group, created if it does not exist."""
        group_path = path / group.uid
        group_path.mkdir(parents=True, exist_ok=True)
        return group_path

    @abstractmethod
    def _run(
        self, context: _Context, varied: Sequence[SensitivityParameter], path: Path
    ) -> None:
        """Run the analysis and write its tables and figures.

        Args:
            context: the simulator, the model, the outputs and observables and
                the seed of the run.
            varied: the parameters the analysis varies, see `_check`.
            path: the folder of the analysis.
        """


@dataclass(frozen=True, kw_only=True)
class LocalAnalysis(Analysis):
    """The local sensitivities of every group, each parameter by `+-delta`, as sbmlsim 0.8.5.

    Per group the designs `sampling.local` of `delta * k / n_var`,
    `k = 1, ..., n_var` (their values are those of the simulation of the
    group), the indices `sensitivity.local`, and the table of sbmlsim 0.8.5,
    `local_d<delta>_nvar<n_var>_<k>_<group>.tsv` (a row per parameter and
    output: `parameter`, `output`, `effect`, `Sij`, `|Sij|`,
    `classification`, `method`, `difference`, `n_var`, sorted by `|Sij|`;
    `Sij` the mean of `normalized` over the designs, as sbmlsim 0.8.5 averaged
    its steps), and its heatmap (`.png`, the parameters with their names, the
    outputs by name, cutoff 0.05, no clustering), in the folder of the
    analysis.

    Attributes:
        delta: the relative change of a parameter, in (0, 1).
        n_var: the number of steps of sbmlsim 0.8.5 up to `delta`, at least 1.
    """

    delta: float = 0.01
    n_var: int = 3

    def __post_init__(self) -> None:
        """Check the groups, the change and its steps.

        Raises:
            ValueError: see `Analysis`, for a change outside of (0, 1), which
                `sampling.local` refuses, and for a number of steps which is no
                positive int.
        """
        super().__post_init__()
        if not 0 < self.delta < 1:
            raise ValueError(
                f"The delta of a local analysis is in (0, 1), not {self.delta}."
            )
        if (
            isinstance(self.n_var, bool)
            or not isinstance(self.n_var, int)
            or self.n_var < 1
        ):
            raise ValueError(
                f"The n_var of a local analysis is a positive int, not {self.n_var!r}."
            )

    def _run(
        self, context: _Context, varied: Sequence[SensitivityParameter], path: Path
    ) -> None:
        """Run the analysis, see the class."""
        targets = [p.uid for p in varied]
        prefix = f"local_d{self.delta}_nvar{self.n_var}"
        path.mkdir(parents=True, exist_ok=True)
        for kg, group in enumerate(self.groups):
            results = [
                sensitivity.local(
                    context.scan(
                        group,
                        sampling.local(
                            targets,
                            delta=kv / self.n_var * self.delta,
                            model=context.model_path(group),
                            simulation=group.simulation,
                        ),
                    )
                )
                for kv in range(1, self.n_var + 1)
            ]
            mean = self._mean(results, context.outputs)
            self._table(mean, results, varied, context.outputs).to_csv(
                path / f"{prefix}_{kg:>02}_{group.uid}.tsv", sep="\t", index=False
            )
            _heatmap(
                _labelled(mean, varied, context.outputs),
                "normalized",
                cutoff=CUTOFF_085,
                cluster_rows=False,
                cmap="seismic",
                vmin=-2.0,
                vmax=2.0,
                path=path / f"{prefix}_{kg:>02}_{group.uid}.png",
            )

    def _mean(
        self, results: Sequence[SensitivityResult], outputs: Sequence[SensitivityOutput]
    ) -> SensitivityResult:
        """Get the mean of `normalized` over the steps, the `Sij` of sbmlsim 0.8.5."""
        ds = results[0].ds.copy()
        for output in outputs:
            key = f"{output.observable}.normalized"
            ds[key] = sum(r.ds[key] for r in results) / len(results)
        return SensitivityResult(ds)

    def _table(
        self,
        mean: SensitivityResult,
        results: Sequence[SensitivityResult],
        varied: Sequence[SensitivityParameter],
        outputs: Sequence[SensitivityOutput],
    ) -> pd.DataFrame:
        """Get the table of sbmlsim 0.8.5 and log the nonlinear sensitivities.

        A sensitivity whose `raw` values of the steps differ by more than 10 %
        of their mean while it is larger than 0.01 is logged as an error, as
        sbmlsim 0.8.5 did: the output is nonlinear in the range or is not
        computed accurately.
        """
        rows = []
        for p in varied:
            for output in outputs:
                s = float(
                    mean[f"{output.observable}.normalized"].sel({PARAMETER: p.uid})
                )
                raw = np.array(
                    [
                        float(r[f"{output.observable}.raw"].sel({PARAMETER: p.uid}))
                        for r in results
                    ]
                )
                if not np.isclose(raw.mean(), 0.0) and s > 0.01:
                    max_diff = (raw.max() - raw.min()) / raw.mean()
                    if max_diff > 0.10:
                        logger.error(
                            "Large delta difference: %s%% for %s. Parameter '%s' on "
                            "output '%s'.",
                            format(max_diff * 100, ".1f"),
                            raw,
                            p.uid,
                            output.uid,
                        )
                effect, classification = _classes(
                    s, sensitivity_classification, sensitivity_classification_symbol
                )
                rows.append(
                    {
                        "parameter": p.uid,
                        "output": output.uid,
                        "effect": effect,
                        "Sij": s,
                        "|Sij|": np.abs(s),
                        "classification": classification,
                        "method": "LocalSensitivity",
                        "difference": self.delta,
                        "n_var": self.n_var,
                    }
                )
        df = pd.DataFrame(rows)
        return df.sort_values(
            by="|Sij|", ascending=False, na_position="last", ignore_index=True
        )


#: the severities of the impairments of a library model, 0 in the healthy
#: individual, by their largest value (Child-Pugh C, NYHA IV)
SEVERITIES: Mapping[str, float] = MappingProxyType(
    {"CP_severity": 3.0, "HF_severity": 4.0}
)


@dataclass(frozen=True, kw_only=True)
class NonzeroLocalAnalysis(LocalAnalysis):
    """The local analysis of the parameters whose value is not 0.

    A local design changes a value relative to itself, so a severity at 0
    (the healthy individual) is not changed and its sensitivity is 0/0; the
    sampled analyses vary the severities within their bounds.
    """

    @override
    def varied(
        self, parameters: Sequence[SensitivityParameter]
    ) -> list[SensitivityParameter]:
        return [p for p in super().varied(parameters) if p.value != 0]


def notes_sources(sbase: libsbml.SBase) -> list[str]:
    """Get the sources in the notes of an element of a library model.

    The library writes the sources of a value as a list after a paragraph
    `Sources:` of its notes.

    Args:
        sbase: the element.

    Returns:
        The sources, one text per list item; empty without notes or sources.
    """
    if not sbase.isSetNotes():
        return []
    root = ElementTree.fromstring(sbase.getNotesString())
    elements = list(root.iter())
    for k, element in enumerate(elements):
        if element.tag.endswith("}p") and "".join(element.itertext()).strip() == (
            "Sources:"
        ):
            lists = [e for e in elements[k + 1 :] if e.tag.endswith("}ul")]
            if lists:
                return [
                    " ".join("".join(item.itertext()).split())
                    for item in lists[0]
                    if item.tag.endswith("}li")
                ]
    return []


def library_parameters(
    model_path: Path,
    ids: Sequence[str],
    *,
    fit_ids: Collection[str],
    scaling_ids: Collection[str],
    references: Mapping[str, str] = MappingProxyType({}),
    relative: float = 0.15,
) -> list[SensitivityParameter]:
    """Get the parameters of the analyses of a model of the PBPK library.

    A library model has hundreds of constants, most of them initial
    assignments of the covariates, which an analysis of all constants would
    vary as well; the analysis names the parameters it varies. They are the
    `ids`, in this order, with these bounds:

    - fitted parameters (`fit_ids`, type fitted): `value * (1 +- relative)`;
    - scaling factors (`scaling_ids`, type scaling, 1 in the healthy
      individual): `1 +- relative`;
    - severities (`SEVERITIES`, type scaling): 0 to `relative * top`, since
      relative bounds around 0 would be 0;
    - any other parameter: `value * (1 +- relative)`.

    The reference of a parameter is the sources of its notes (`;` separated),
    replaced by the text in `references` (e.g. the guideline of the classes of
    the renal function, whose notes give no source). A parameter of the type
    not set with a reference is of the type data.

    Args:
        model_path: the SBML model.
        ids: the parameters of the analysis.
        fit_ids: the fitted parameters among the `ids`.
        scaling_ids: the scaling factors among the `ids`.
        references: references by id, which replace those of the notes.
        relative: the relative half width of the bounds.

    Returns:
        The parameters of the `ids` with bounds, types and references.

    Raises:
        ValueError: naming every id which is no parameter of the model, and
            every id of `fit_ids`, `scaling_ids` or `references` which is not
            in `ids`.
    """
    by_uid = {
        p.uid: p
        for p in parameters_from_sbml(model_path, exclude_na=False, exclude_zero=False)
    }
    missing = [uid for uid in ids if uid not in by_uid]
    if missing:
        raise ValueError(
            f"The parameters {missing} of the sensitivity analysis are no "
            f"parameters of the model {model_path.name}."
        )
    outside = sorted(
        {*fit_ids, *scaling_ids, *references}.difference(ids),
    )
    if outside:
        raise ValueError(
            f"The ids {outside} of the sensitivity analysis are not in its parameters."
        )
    parameters = [by_uid[uid] for uid in ids]

    # the fitted parameters, their bounds are set relative to the value below
    set_bounds(
        parameters, [(uid, np.nan, np.nan, ParameterType.FIT) for uid in fit_ids]
    )
    # the scaling factors of the organ functions and of the protein activities
    set_bounds(
        parameters,
        [
            (uid, 1 - relative, 1 + relative, ParameterType.SCALING)
            for uid in scaling_ids
        ],
    )
    # the severities, whose relative bounds around 0 would be 0
    set_bounds(
        parameters,
        [
            (uid, 0.0, relative * top, ParameterType.SCALING)
            for uid, top in SEVERITIES.items()
        ],
    )

    # the references of the values, as the sources of their notes
    model: libsbml.Model = libsbml.readSBMLFromFile(str(model_path)).getModel()
    texts = {
        uid: "; ".join(sources)
        for uid in ids
        if (sources := notes_sources(model.getParameter(uid)))
    }
    texts.update(references)
    for uid, reference in texts.items():
        p = by_uid[uid]
        p.reference = reference
        if p.type == ParameterType.NA:
            p.type = ParameterType.DATA

    # every parameter without bounds varies around its value
    set_relative_bounds(parameters, relative)
    return parameters


@dataclass(frozen=True, kw_only=True)
class _SampledAnalysis(Analysis):
    """An analysis of a design of `n` samples, see `sample_size`.

    Attributes:
        n: the number of samples of the design.
    """

    n: int

    def __post_init__(self) -> None:
        """Check the groups and the number of samples.

        Raises:
            ValueError: see `Analysis`, and for a number of samples which the
                design refuses.
        """
        super().__post_init__()
        self._check_size(self.n)

    def _check_size(self, n: int) -> None:
        """Check a number of samples, as the design of the analysis does.

        Raises:
            ValueError: if it is no positive int.
        """
        if isinstance(n, bool) or not isinstance(n, int) or n < 1:
            raise ValueError(
                f"The number of samples of {type(self).__name__} is a positive "
                f"int, not {n!r}."
            )

    def _samples(self) -> int:
        """Get the number of samples of a run, `n` or the one of a smoke run.

        Raises:
            ValueError: for a value of `PKPDSIM_SENSITIVITY_N` which is no
                positive int or which the design refuses.
        """
        n = sample_size(self.n)
        try:
            self._check_size(n)
        except ValueError as err:
            err.add_note(
                f"{SAMPLE_SIZE} replaces the number of samples {self.n} of "
                f"{type(self).__name__} by {n}."
            )
            raise
        return n

    def _check(
        self, parameters: Sequence[SensitivityParameter]
    ) -> list[SensitivityParameter]:
        """Check the parameters, their distributions and the number of samples.

        Raises:
            ValueError: see `Analysis`, `distributions` and `_samples`.
        """
        varied = super()._check(parameters)
        distributions(varied)
        self._samples()
        return varied


@dataclass(frozen=True, kw_only=True)
class SamplingAnalysis(_SampledAnalysis):
    """The uncertainty of the outputs over a Latin hypercube of the parameters.

    One design `sampling.lhs` of `sample_size(n)` points for every group; per
    group the parameters and outputs of every point, `samples.tsv` (sbmlsim
    0.8.5 kept them in its cache only). The tables and figures of sbmlsim 0.8.5
    compare the groups: `sampling_N<n>.tsv` (a row per group and output:
    `gid`, `gname`, `uid`, `name`, `N`, `mean`, `median`, `std`, `cv`, `min`,
    `q005`, `q095`, `max`, the uncertainty `U = (P97.5 - P2.5) / P50`, its
    `symbol` and `classification`, `unit`; `std` of the population),
    `sampling_N<n>_compact.tsv` and its LaTeX table `sampling_N<n>.tex` (an
    output per row, `mean (cv in %)` per group), and the figures of the samples
    of the parameters and of the outputs per group,
    `sampling_N<n>_sensitivity_samples.png` and `..._outputs.png`.

    Attributes:
        n: the number of samples, see `sample_size`.
    """

    n: int = 1000

    def _run(
        self, context: _Context, varied: Sequence[SensitivityParameter], path: Path
    ) -> None:
        """Run the analysis, see the class."""
        n = self._samples()
        design = sampling.lhs(distributions(varied), n=n, seed=context.seed)
        observables = {o.observable: o.uid for o in context.outputs}
        results: dict[str, ScanResult] = {}
        for group in self.groups:
            group_path = self._group_path(path, group)
            res = context.scan(group, context.design(group, design, varied))
            results[group.uid] = res
            res.ds[list(observables)].to_dataframe().rename(columns=observables).to_csv(
                group_path / "samples.tsv", sep="\t"
            )
        self._write_085(context, varied, results, n, path)

    def _write_085(
        self,
        context: _Context,
        varied: Sequence[SensitivityParameter],
        results: Mapping[str, ScanResult],
        n: int,
        path: Path,
    ) -> None:
        """Write the tables and figures of sbmlsim 0.8.5, see the class."""
        prefix = f"sampling_N{n}"
        units = {
            o.uid: context.unit(results[self.groups[0].uid], o) for o in context.outputs
        }
        rows = []
        stats: dict[tuple[str, str], dict[str, float]] = {}
        for group in self.groups:
            for output in context.outputs:
                values = np.asarray(
                    results[group.uid].ds[output.observable].values, dtype=float
                )
                stats[group.uid, output.uid] = _uncertainty(values)
                symbol, classification = _classes(
                    stats[group.uid, output.uid]["U"],
                    uncertainty_classification,
                    uncertainty_classification_symbol,
                )
                rows.append(
                    {
                        "gid": group.uid,
                        "gname": group.name,
                        "uid": output.uid,
                        "name": output.name,
                        "N": n,
                        **stats[group.uid, output.uid],
                        "symbol": symbol,
                        "classification": classification,
                        "unit": units[output.uid],
                    }
                )
        pd.DataFrame(rows).to_csv(path / f"{prefix}.tsv", sep="\t", index=False)

        compact = pd.DataFrame(
            [
                {
                    "output": output.name,
                    **{
                        g.uid: (
                            f"{stats[g.uid, output.uid]['mean']:.3g} "
                            f"({stats[g.uid, output.uid]['cv'] * 100:.1f})"
                        )
                        for g in self.groups
                    },
                    "unit": units[output.uid],
                }
                for output in context.outputs
            ]
        )
        compact.to_csv(path / f"{prefix}_compact.tsv", sep="\t", index=False)
        latex = compact.to_latex(None, index=False)
        latex = latex.replace("∞", r"$\infty$").replace("*", r"$\cdot$")
        (path / f"{prefix}.tex").write_text(latex, encoding="utf-8")

        samples = [
            (
                f"{p.uid}: {p.name} [{p.unit if p.unit else 'AU'}]",
                [
                    np.asarray(results[g.uid].ds[p.uid].values, dtype=float)
                    for g in self.groups
                ],
            )
            for p in varied
        ]
        _plot_groups(
            samples,
            self.groups,
            path / f"{prefix}_sensitivity_samples.png",
            context.seed,
        )
        outputs = [
            (
                f"{o.name} [{units[o.uid] if units[o.uid] else 'AU'}]",
                [
                    np.asarray(results[g.uid].ds[o.observable].values, dtype=float)
                    for g in self.groups
                ],
            )
            for o in context.outputs
        ]
        _plot_groups(
            outputs,
            self.groups,
            path / f"{prefix}_sensitivity_outputs.png",
            context.seed,
        )


def _uncertainty(values: np.ndarray) -> dict[str, float]:
    """Get the statistics of the samples of an output, as sbmlsim 0.8.5 did.

    `std` is the standard deviation of the population, `cv = std / mean`, the
    uncertainty `U = (P97.5 - P2.5) / P50`.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        return {
            "mean": float(np.mean(values)),
            "median": float(np.median(values)),
            "std": float(np.std(values)),
            "cv": float(np.std(values) / np.mean(values)),
            "min": float(np.min(values)),
            "q005": float(np.quantile(values, q=0.05)),
            "q095": float(np.quantile(values, q=0.95)),
            "max": float(np.max(values)),
            "U": float(
                (np.percentile(values, 97.5) - np.percentile(values, 2.5))
                / np.percentile(values, 50)
            ),
        }


@dataclass(frozen=True, kw_only=True)
class _GlobalAnalysis(_SampledAnalysis):
    """A variance based analysis: the indices `S1` and `ST` of every group.

    One design for every group; per group the tables of `S1` and `ST` per
    output, `<output>.S1.tsv` and `<output>.ST.tsv` (sbmlsim 0.8.5 wrote no
    tables of a global analysis), and the figures of sbmlsim 0.8.5: the
    heatmaps of `S1` and `ST` (`<prefix>_sensitivity_<k>_<group>_S1.png`,
    `..._ST.png`, the parameters with their names, the outputs by name, cutoff
    0.05, no clustering, viridis from 0 to 1) and the indices of every output
    with their intervals, `<prefix>_sensitivity_<k>_<group>_<o>_<output>.png`
    (`plot_indices`).

    Attributes:
        n: the base number of samples, see `sample_size`.
    """

    @abstractmethod
    def _prefix(self, n: int) -> str:
        """Get the prefix of the files of sbmlsim 0.8.5."""

    @abstractmethod
    def _design(
        self, parameters: Sequence[SensitivityParameter], n: int, seed: int
    ) -> Dimension:
        """Get the design over the parameters the analysis varies."""

    @abstractmethod
    def _indices(self, result: ScanResult) -> SensitivityResult:
        """Get the indices of the result of the design."""

    def _run(
        self, context: _Context, varied: Sequence[SensitivityParameter], path: Path
    ) -> None:
        """Run the analysis, see the class."""
        n = self._samples()
        prefix = self._prefix(n)
        design = self._design(varied, n, context.seed)
        for kg, group in enumerate(self.groups):
            group_path = self._group_path(path, group)
            result = self._indices(
                context.scan(group, context.design(group, design, varied))
            )
            context.write_indices(result, ["S1", "ST"], group_path)
            figures = f"{prefix}_sensitivity_{kg:>02}_{group.uid}"
            for ko, output in enumerate(context.outputs):
                sensitivity.plot_indices(
                    result,
                    output.observable,
                    title=f"{output.name} ({group.name})",
                    path=path / f"{figures}_{ko:>03}_{output.uid}.png",
                )
            labelled = _labelled(result, varied, context.outputs)
            for key in ["ST", "S1"]:
                _heatmap(
                    labelled,
                    key,
                    cutoff=CUTOFF_085,
                    cluster_rows=False,
                    cmap="viridis",
                    vmin=0.0,
                    vmax=1.0,
                    path=path / f"{figures}_{key}.png",
                )


@dataclass(frozen=True, kw_only=True)
class FastAnalysis(_GlobalAnalysis):
    """The indices of the extended FAST, `sampling.fast` and `sensitivity.fast`.

    The files of sbmlsim 0.8.5 have the prefix `fast_M4_N<n>`.

    Attributes:
        n: the number of samples per parameter, more than 64 (`4 m^2` of the
            `m = 4` harmonics of `sampling.fast`), see `sample_size`.
    """

    n: int = 1000

    def _check_size(self, n: int) -> None:
        """Check a number of samples, more than 64, as `sampling.fast` does.

        Raises:
            ValueError: if it is no int above 64.
        """
        super()._check_size(n)
        if n <= 64:
            raise ValueError(
                f"The number of samples of FastAnalysis is more than 64, not {n}."
            )

    def _prefix(self, n: int) -> str:
        """Get the prefix of the files of sbmlsim 0.8.5, `fast_M4_N<n>`."""
        return f"fast_M4_N{n}"

    def _design(
        self, parameters: Sequence[SensitivityParameter], n: int, seed: int
    ) -> Dimension:
        """Get the design over the parameters, see `sampling.fast`."""
        return sampling.fast(distributions(parameters), n=n, seed=seed)

    def _indices(self, result: ScanResult) -> SensitivityResult:
        """Get the indices, see `sensitivity.fast`."""
        return sensitivity.fast(result)


@dataclass(frozen=True, kw_only=True)
class SobolAnalysis(_GlobalAnalysis):
    """The Sobol indices, `sampling.sobol` and `sensitivity.sobol`.

    The files of sbmlsim 0.8.5 have the prefix `sobol_N<n>`.

    Attributes:
        n: the base number of samples, a power of two of at least 2, see
            `sample_size`.
    """

    n: int = 4096

    def _check_size(self, n: int) -> None:
        """Check a number of samples, a power of two, as `sampling.sobol` does.

        Raises:
            ValueError: if it is no power of two of at least 2.
        """
        super()._check_size(n)
        if n < 2 or n & (n - 1):
            raise ValueError(
                "The number of samples of SobolAnalysis is a power of two of at "
                f"least 2, not {n}."
            )

    def _prefix(self, n: int) -> str:
        """Get the prefix of the files of sbmlsim 0.8.5, `sobol_N<n>`."""
        return f"sobol_N{n}"

    def _design(
        self, parameters: Sequence[SensitivityParameter], n: int, seed: int
    ) -> Dimension:
        """Get the design over the parameters, see `sampling.sobol`."""
        return sampling.sobol(distributions(parameters), n=n, seed=seed)

    def _indices(self, result: ScanResult) -> SensitivityResult:
        """Get the indices, see `sensitivity.sobol`."""
        return sensitivity.sobol(result)


class CommandLine(NamedTuple):
    """The words of the command line of a sensitivity script.

    Attributes:
        analyses: the names of the analyses, the default analyses without one.
        choices: the words of the choices of the script, in their order, none
            if none is given.
        workers: the processes of `--workers`, `None` if not given or `0`
            (`n_workers=None`, sbmlsim chooses).
    """

    analyses: tuple[str, ...]
    choices: tuple[str, ...]
    workers: int | None = None


def _workers(text: str) -> int:
    """Read the value of `--workers`: `0` or a number of processes of at least 1."""
    try:
        value = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{text!r} is no whole number") from None
    if value < 0:
        raise argparse.ArgumentTypeError(f"{value} is negative (0 or at least 1)")
    return value


def _write_table(
    parameters: Sequence[SensitivityParameter],
    path: Path,
    builds: Mapping[Sex, Mapping[str, float]] | None = None,
) -> None:
    """Write the table of the parameters into a folder and print it."""
    path.mkdir(parents=True, exist_ok=True)
    df_parameters = write_parameters(parameters, path, builds)
    with pd.option_context("display.float_format", "{:.5g}".format):
        console.print(df_parameters)


def _builds(
    model: PKPDModel, sexes: Iterable[Sex], parameters: Sequence[SensitivityParameter]
) -> dict[Sex, dict[str, float]]:
    """Get the values of the parameters in the builds of the sexes, by sex."""
    uids = [p.uid for p in parameters]
    return {
        sex: build_values(model.model_path_of(sex), uids)
        for sex in sorted(sexes, key=list(Sex).index)
    }


def _report(items: Mapping[str, str], path: Path) -> None:
    """Print the report of a run."""
    display.section("Report", icon=display.ICON_REPORT)
    display.key_values(dict(items))
    display.link("results", path)


@dataclass(frozen=True, kw_only=True)
class SensitivityRunner:
    """The sensitivity analyses of a model, run by its script.

    The script gives what is the model's own (the parameters, the outputs with
    their observables, the analyses with their groups) and runs them with `main`
    from its command line, or with `run`. A run writes the table of the
    parameters (`write_parameters`) into its folder,
    `model.results_dir / "sensitivity"`, and the results of an analysis into
    `<analysis>` below it, or below a subfolder of it, see `Analysis`. A model with several
    simulations has a runner per simulation, see `run_subfolders`. What can be
    checked is checked when the runner is created, and the rest at the start
    of a run, before anything is written.

    Attributes:
        model: the model object: the models which are simulated (the build
            of the sex of every group, `SensitivityGroup.sex`), the folder of
            the results and the settings of the simulator.
        parameters: the function which gives the parameters of the analyses
            with their bounds; it reads the model, so it is called by a run.
        outputs: the outputs of the analyses, at least one, the ids, the names
            and the observables unique; stored as a tuple.
        observables: the observables of the simulation of a group, which give
            the observables of the outputs.
        analyses: the analyses by name, the name on the command line and the
            folder; at least one, stored read-only.
        default: the analyses which run without a name on the command line, at
            least one.
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
            ValueError: for no output, two outputs of one id, name or
                observable, an output whose observable is none of the
                observables, no analysis, no default analysis, or a default
                analysis which is none of the analyses or is given twice, or a
                group of a sex which the model has no build of.
            TypeError: for a default which is a string.
        """
        object.__setattr__(self, "outputs", tuple(self.outputs))
        object.__setattr__(self, "observables", tuple(self.observables))
        object.__setattr__(self, "analyses", MappingProxyType(dict(self.analyses)))
        object.__setattr__(
            self, "default", _sequence_of_names(self.default, "default analyses")
        )
        if not self.outputs:
            raise ValueError("A sensitivity analysis needs at least one output.")
        repeated = _repeated(o.uid for o in self.outputs)
        if repeated:
            raise ValueError(f"The outputs {repeated} are given twice.")
        repeated = _repeated(o.name for o in self.outputs)
        if repeated:
            raise ValueError(f"The names {repeated} of the outputs are given twice.")
        repeated = _repeated(o.observable for o in self.outputs)
        if repeated:
            raise ValueError(f"The observables {repeated} are outputs twice.")
        ids = {observable.id for observable in self.observables}
        missing = [o.observable for o in self.outputs if o.observable not in ids]
        if missing:
            raise ValueError(f"The outputs name no observables {missing}.")
        if not self.analyses:
            raise ValueError("A sensitivity analysis needs at least one analysis.")
        if not self.default:
            raise ValueError(
                "A sensitivity analysis needs at least one default analysis."
            )
        unknown = [name for name in self.default if name not in self.analyses]
        if unknown:
            raise ValueError(
                f"The default analyses {unknown} are none of {list(self.analyses)}."
            )
        repeated = _repeated(self.default)
        if repeated:
            raise ValueError(f"The default analyses {repeated} are given twice.")
        for name, analysis in self.analyses.items():
            for group in analysis.groups:
                try:
                    self.model.model_path_of(group.sex)
                except ValueError as err:
                    raise ValueError(
                        f"The group '{group.uid}' of the analysis '{name}' runs "
                        f"the {group.sex} build: {err}"
                    ) from err

    def _names(self, names: Sequence[str] | None) -> tuple[str, ...]:
        """Get the names of the analyses of a run, `default` if `None`.

        Raises:
            TypeError: for names which are a string.
            ValueError: for a name which is no analysis or is given twice.
        """
        names = self.default if names is None else _sequence_of_names(names, "names")
        unknown = [name for name in names if name not in self.analyses]
        if unknown:
            raise ValueError(f"Unknown analyses {unknown}, use {list(self.analyses)}.")
        repeated = _repeated(names)
        if repeated:
            raise ValueError(f"The analyses {repeated} are given twice.")
        return tuple(names)

    def _plan(
        self,
        names: Sequence[str],
        parameters: Sequence[SensitivityParameter],
        where: str = "",
    ) -> dict[str, list[SensitivityParameter]]:
        """Check the named analyses and get the parameters each varies.

        Raises:
            ValueError: see `Analysis._check`, with a note naming the analysis.
        """
        varied: dict[str, list[SensitivityParameter]] = {}
        for name in names:
            try:
                varied[name] = self.analyses[name]._check(parameters)
            except ValueError as err:
                err.add_note(f"The analysis '{name}'{where} cannot run.")
                raise
        return varied

    def _context(
        self, simulator: Simulator, builds: Mapping[Sex, Mapping[str, float]]
    ) -> _Context:
        """Get the context of the analyses of a run."""
        return _Context(
            simulator=simulator,
            model=self.model,
            outputs=self.outputs,
            observables=self.observables,
            seed=self.seed,
            builds=builds,
        )

    def _sexes(self, names: Sequence[str]) -> set[Sex]:
        """Get the sexes of the groups of the named analyses other than male."""
        return {
            group.sex
            for name in names
            for group in self.analyses[name].groups
            if group.sex is not Sex.MALE
        }

    def _check_designs(
        self,
        names: Sequence[str],
        varied: Mapping[str, Sequence[SensitivityParameter]],
        builds: Mapping[Sex, Mapping[str, float]],
        where: str = "",
    ) -> None:
        """Check that the global designs of the named analyses carry to every build.

        Raises:
            ValueError: see `design_factors`, with a note naming the analysis
                and the group.
        """
        for name in names:
            analysis = self.analyses[name]
            if isinstance(analysis, LocalAnalysis):
                # a local design is drawn around the values of each build
                continue
            for group in analysis.groups:
                if group.sex not in builds:
                    continue
                try:
                    design_factors(varied[name], builds[group.sex])
                except ValueError as err:
                    err.add_note(
                        f"The group '{group.uid}' of the analysis '{name}'{where} "
                        f"runs the {group.sex} build."
                    )
                    raise

    def simulator_settings(self, workers: int | None) -> SimulatorSettings | None:
        """Get the settings of the simulator of a run which caps the processes.

        Args:
            workers: the processes of `--workers` (`CommandLine.workers`).

        Returns:
            The settings of the model with `n_workers=workers`, `None` (the
            default of `run`, `n_workers=None`) if `workers` is `None`.
        """
        if workers is None:
            return None
        return dataclasses.replace(self.model.simulator, n_workers=workers)

    def _simulator(self, simulator: SimulatorSettings | None) -> Simulator:
        """Get the simulator of a run, of the model with `n_workers=None` by default."""
        settings = (
            dataclasses.replace(self.model.simulator, n_workers=None)
            if simulator is None
            else simulator
        )
        return create_simulator(settings)

    def _execute(
        self,
        names: Sequence[str],
        varied: Mapping[str, Sequence[SensitivityParameter]],
        path: Path,
        simulator: Simulator,
        builds: Mapping[Sex, Mapping[str, float]] | None = None,
    ) -> None:
        """Run the named analyses, every one into `path / <name>`."""
        context = self._context(simulator, builds or {})
        for name in names:
            self.analyses[name]._run(context, varied[name], path / name)

    def run(
        self,
        names: Sequence[str] | None = None,
        *,
        path: Path | None = None,
        subfolder: str | None = None,
        simulator: SimulatorSettings | None = None,
    ) -> Path:
        """Run analyses and write the table of the parameters.

        The names, the parameters (their distributions, whether an analysis
        varies one, whether a global design carries to the build of every
        group, `design_factors`) and the numbers of samples of a smoke run are
        checked before anything is written. The values of the table are those
        of the male build; with groups of another sex the table has their
        values as well (`value_<sex>`).

        Args:
            names: the names of the analyses, `default` if `None`.
            path: the folder of the results and of the table of the parameters,
                `model.results_dir / "sensitivity"` if `None`.
            subfolder: the folder of the analyses below `path`; the table of
                the parameters stays in `path`. Several simulations of a model
                run with `run_subfolders`, which writes the table once.
            simulator: the settings of the simulator, the settings of the model
                with `n_workers=None` (sbmlsim chooses the processes of the
                pool, the script has the guard `if __name__ == "__main__":`) if
                `None`.

        Returns:
            The folder of the results, `path`.

        Raises:
            TypeError: for names which are a string.
            ValueError: for a name which is no analysis or is given twice, a
                subfolder which does not lie below `path`, and see
                `distributions`, `sample_size` and `design_factors` (a global
                design of a group of another sex); before anything is written.
        """
        names = self._names(names)
        if subfolder is not None:
            _check_subfolder(subfolder)
        path = self.model.results_dir / "sensitivity" if path is None else path
        parameters = list(self.parameters())
        varied = self._plan(names, parameters)
        builds = _builds(self.model, self._sexes(names), parameters)
        self._check_designs(names, varied, builds)

        _write_table(parameters, path, builds)
        analyses_path = path if subfolder is None else path / subfolder
        self._execute(names, varied, analyses_path, self._simulator(simulator), builds)
        _report({"analyses": ", ".join(names)}, analyses_path)
        return path

    def parse_command_line(
        self, argv: Sequence[str] | None = None, *, choices: Sequence[str] = ()
    ) -> CommandLine:
        """Read the names of the analyses, and of choices of the script, from a command line.

        A script with a choice of its own, e.g. the simulation of a model with
        several, reads it here with the names of the analyses, in any order
        (`local sampling multi`), and runs the analyses itself.

        Args:
            argv: the arguments, those of the process if `None`.
            choices: the words of the choices of the script, none an analysis.

        Returns:
            The names of the analyses, the default analyses without one, the
            words of the choices and the `--workers` (`None` if absent or `0`),
            which a script passes as `simulator=self.simulator_settings(workers)`.

        Raises:
            TypeError: for arguments which are a string.
            ValueError: for a choice which is the name of an analysis.
            SystemExit: with 2 for a word which is neither an analysis nor a
                choice or which is given twice (argparse), with 0 for `--help`.
        """
        if argv is not None:
            argv = _sequence_of_names(argv, "argv")
        choices = _sequence_of_names(choices, "choices")
        clashes = sorted(set(choices) & set(self.analyses))
        if clashes:
            raise ValueError(f"The choices {clashes} are names of analyses.")
        words = ", ".join([*self.analyses, *choices])
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
        # the words are checked here: argparse of python 3.13 checks the
        # default of a positional `nargs="*"` with `choices` as one choice
        parser.add_argument(
            "words",
            nargs="*",
            metavar="ANALYSIS" if not choices else "WORD",
            help=(
                f"the analyses, of {', '.join(self.analyses)} "
                f"(default: {' '.join(self.default)})"
                + (
                    f"; the choices of the script, of {', '.join(choices)}"
                    if choices
                    else ""
                )
            ),
        )
        parser.add_argument(
            "--workers",
            type=_workers,
            default=0,
            metavar="N",
            help=(
                "the processes of the pool of a run, at least 1 (a pool needs "
                "about 300 MB per process for a PBPK model); 0 or absent: "
                "sbmlsim chooses (every CPU for large designs)"
            ),
        )
        args = parser.parse_intermixed_args(argv)
        given: list[str] = args.words
        unknown = [w for w in given if w not in self.analyses and w not in choices]
        if unknown:
            parser.error(
                f"invalid choice: {', '.join(repr(w) for w in unknown)} "
                f"(choose from {words})"
            )
        repeated = _repeated(given)
        if repeated:
            parser.error(f"{', '.join(repr(w) for w in repeated)} given twice")
        analyses = tuple(w for w in given if w in self.analyses)
        return CommandLine(
            analyses=analyses or tuple(self.default),
            choices=tuple(w for w in given if w in choices),
            workers=args.workers or None,
        )

    def main(self, argv: Sequence[str] | None = None) -> None:
        """Run the analyses named on the command line, see `parse_command_line`, `run`.

        Args:
            argv: the arguments, those of the process if `None`: the names of
                the analyses, `default` without one; an unknown name exits with
                2; `--workers N` caps the processes of the pool.
        """
        words = self.parse_command_line(argv)
        self.run(words.analyses, simulator=self.simulator_settings(words.workers))


def run_subfolders(
    runners: Mapping[str, SensitivityRunner],
    names: Sequence[str] | None = None,
    *,
    subfolders: Sequence[str] | None = None,
    path: Path | None = None,
    simulator: SimulatorSettings | None = None,
) -> Path:
    """Run the analyses of several simulations of a model, each into a subfolder.

    A model with several simulations (dulaglutide: a single dose and weekly
    doses) has a runner per simulation, with its outputs and groups, by the
    name of its subfolder. The table of the parameters is written and printed
    once into `path`, the analyses of a runner run into
    `path / <subfolder> / <analysis> / <group>`. Everything is checked before
    anything is written, for every runner: the names, the parameters and the
    sizes, as `SensitivityRunner.run` checks them.

    Args:
        runners: the runners by the name of their subfolder.
        names: the names of the analyses, the default analyses of every runner
            if `None`.
        subfolders: the subfolders which run, every one if `None`.
        path: the folder of the results, `model.results_dir / "sensitivity"` if
            `None`.
        simulator: the settings of the simulator, see `SensitivityRunner.run`.

    Returns:
        The folder of the results, `path`.

    Raises:
        TypeError: for names or subfolders which are a string.
        ValueError: for no runner, no subfolder, a subfolder which has no
            runner, is given twice or does not lie below `path`, runners of
            different models or with different parameters (they share one
            table), and see `SensitivityRunner.run`; before anything is
            written.
    """
    if not runners:
        raise ValueError(
            "A run of subfolders needs at least one runner, there is no runner."
        )
    keys = (
        tuple(runners)
        if subfolders is None
        else _sequence_of_names(subfolders, "subfolders")
    )
    if not keys:
        raise ValueError(
            f"A run of subfolders needs at least one subfolder, of {list(runners)}."
        )
    unknown = [key for key in keys if key not in runners]
    if unknown:
        raise ValueError(
            f"The subfolders {unknown} have no runner, use {list(runners)}."
        )
    repeated = _repeated(keys)
    if repeated:
        raise ValueError(f"The subfolders {repeated} are given twice.")
    for key in keys:
        _check_subfolder(key)
    first = runners[keys[0]]
    others = [key for key in keys if runners[key].model != first.model]
    if others:
        raise ValueError(
            f"The runners of the subfolders {others} have another model than the "
            f"one of '{keys[0]}'; the subfolders of a run share one model."
        )

    # the parameters of every runner, read once per function
    read: dict[int, list[SensitivityParameter]] = {}
    parameters: dict[str, list[SensitivityParameter]] = {}
    for key in keys:
        function = runners[key].parameters
        if id(function) not in read:
            read[id(function)] = list(function())
        parameters[key] = read[id(function)]
    table = parameters_to_df(parameters[keys[0]])
    different = [
        key for key in keys if not parameters_to_df(parameters[key]).equals(table)
    ]
    if different:
        raise ValueError(
            f"The runners of the subfolders {different} give other parameters than "
            f"the one of '{keys[0]}'; the subfolders of a run share one table of "
            "the parameters."
        )

    plans: dict[str, tuple[tuple[str, ...], dict[str, list[SensitivityParameter]]]] = {}
    for key in keys:
        runner = runners[key]
        try:
            run_names = runner._names(names)
        except ValueError as err:
            err.add_note(f"The subfolder '{key}' cannot run.")
            raise
        plans[key] = (
            run_names,
            runner._plan(run_names, parameters[key], f" of the subfolder '{key}'"),
        )

    sexes = set().union(*(runners[key]._sexes(plans[key][0]) for key in keys))
    # the values of every parameter of every subfolder (the tables are equal,
    # the union keeps a parameter of any subfolder in the values)
    union = list({p.uid: p for key in keys for p in parameters[key]}.values())
    builds = _builds(first.model, sexes, union)
    for key in keys:
        run_names, varied = plans[key]
        runners[key]._check_designs(
            run_names, varied, builds, f" of the subfolder '{key}'"
        )

    path = first.model.results_dir / "sensitivity" if path is None else path
    _write_table(parameters[keys[0]], path, builds)
    simulator_ = first._simulator(simulator)
    for key in keys:
        run_names, varied = plans[key]
        runners[key]._execute(run_names, varied, path / key, simulator_, builds)
    _report(
        {
            "subfolders": ", ".join(keys),
            "analyses": "; ".join(f"{key}: {', '.join(plans[key][0])}" for key in keys),
        },
        path,
    )
    return path
