"""The parameters of the sensitivity analyses of a model.

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
"""

from __future__ import annotations

import logging
import math
import os
from collections.abc import Iterable, Sequence
from dataclasses import asdict, dataclass
from enum import StrEnum
from pathlib import Path

import libsbml
import numpy as np
import pandas as pd
import roadrunner
from sbmlsim.simulation.sampling import Uniform
from sbmlutils.report.units import udef_to_string

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
        ValueError: naming the parameters without bounds.
    """
    missing = sorted(
        p.uid
        for p in parameters
        if math.isnan(p.lower_bound) or math.isnan(p.upper_bound)
    )
    if missing:
        raise ValueError(f"The parameters {missing} have no bounds.")
    return {
        p.uid: Uniform(lower=p.lower_bound, upper=p.upper_bound) for p in parameters
    }


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
