"""The uncertainty of the fitted parameters of a model and the bands of its predictions.

A fit of a model becomes part of the model in three steps:

1. **Adopt** the fit: `run_fit <package> --adopt <fit id>` copies the record of a
   finished fit (its parameters, the profile likelihood and the Fisher
   information) from the results of the fits into the tracked folder of the
   fits next to the model (`PKPDModel.fits_dir`), see `adopt_fit`.
2. **Use** it in the factory: the model lists its adopted fits
   (`PKPDModel.fits`), and a model file takes the value of a fitted parameter
   from the record, together with its uncertainty, see `FitRecords.parameter`:

   ```python
   # not executed
   FITS = MODEL.fit_records()
   Parameter("k", **FITS.parameter("k", unit=U.per_min), name="rate")
   ```

   The uncertainty is a distrib `Uncertainty` on the SBML parameter: the
   confidence interval of the profile likelihood (the Fisher information when
   the fit has no profile), the standard error of the Fisher information, the
   distribution the draws of the parameter come from and the provenance in its
   notes, see `ParameterUncertainty`. `check_fit_parameters` checks that the
   model has them.
3. **Draw** it: `run_simulation <package> --uncertainty` reads the distributions of the
   parameters of the model (`parameter_distributions`), draws them into the
   dimension `uncertainty` (`uncertainty_dimension`) and runs every simulation
   over the draws; every curve of a simulation becomes the band of a quantile
   range with the median, see `apply_uncertainty`.

The draws are independent per parameter: distrib states one distribution per
parameter, the correlations of the Fisher information stay in the record.
"""

from __future__ import annotations

import json
import logging
import math
import re
import shutil
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path
from typing import TYPE_CHECKING, Any

import libsbml
import numpy as np
from sbmlsim.data import Data
from sbmlsim.fit.fisher import FisherInformation
from sbmlsim.fit.options import ParameterScaleType
from sbmlsim.plot.plotting import Curve
from sbmlsim.simulation import Dimension, Scan, Simulation, sampling
from sbmlsim.simulation.sampling import (
    Distribution,
    LogNormal,
    Normal,
    Truncated,
    Uniform,
)
from sbmlsim.task import Task
from sbmlsim.units import ureg
from scipy import stats

if TYPE_CHECKING:
    from sbmlsim.experiment import SimulationExperiment
    from sbmlutils.factory import Uncertainty
    from sbmlutils.factory.units import UnitType

    from pkpdsim.model import PKPDModel

logger = logging.getLogger(__name__)

#: the id of the dimension of the draws of the fitted parameters
UNCERTAINTY_DIMENSION = "uncertainty"

#: the files of a fit which an adopted fit keeps; the first ones are required
RECORD_FILES: tuple[str, ...] = (
    "parameters.json",
    "optimization_result.json",
    "identifiability.json",
    "identifiability.tsv",
    "fisher.json",
    "fisher.tsv",
    "metrics.tsv",
    "metrics_mappings.tsv",
    "report.txt",
)
#: the files without which a fit cannot be adopted
REQUIRED_FILES: tuple[str, ...] = ("parameters.json", "optimization_result.json")
#: one of these files gives the uncertainty of the parameters
UNCERTAINTY_FILES: tuple[str, ...] = ("identifiability.json", "fisher.json")

#: the relative tolerance of the value of a fitted parameter in the model
VALUE_RTOL = 1e-12


class UncertaintyMethod(StrEnum):
    """Where the uncertainty of a fitted parameter comes from."""

    PROFILE = "profile likelihood"
    FISHER = "Fisher information"


# ----------------------------------------------------------------------------
# adopting a fit
# ----------------------------------------------------------------------------
def adopt_fit(fit_dir: Path, fits_dir: Path, *, force: bool = False) -> Path:
    """Copy the record of a finished fit into the folder of the fits of a model.

    The record is the part of the results of a fit which the model needs: the
    parameters, the settings and the problem (`optimization_result.json`), the
    profile likelihood and the Fisher information with their tables, the
    metrics and the text report (`RECORD_FILES`). The runs, the figures and
    the HTML report stay in the results of the fit.

    Args:
        fit_dir: the directory of the results of the fit, named by its id.
        fits_dir: the folder of the adopted fits of the model.
        force: replace a record of the same id.

    Returns:
        The directory of the record.

    Raises:
        FileNotFoundError: if the fit does not exist, lacks a required file or
            has neither a profile likelihood nor a Fisher information.
        FileExistsError: if the record exists and `force` is not set.
    """
    if not fit_dir.is_dir():
        raise FileNotFoundError(f"The fit '{fit_dir}' does not exist.")
    missing = [name for name in REQUIRED_FILES if not (fit_dir / name).exists()]
    if missing:
        raise FileNotFoundError(f"The fit '{fit_dir.name}' has no {missing}.")
    if not any((fit_dir / name).exists() for name in UNCERTAINTY_FILES):
        raise FileNotFoundError(
            f"The fit '{fit_dir.name}' has neither a profile likelihood nor a "
            f"Fisher information ({list(UNCERTAINTY_FILES)}): run the fit with "
            f"--identifiability."
        )
    target = fits_dir / fit_dir.name
    if target.exists():
        if not force:
            raise FileExistsError(
                f"The fit '{fit_dir.name}' is adopted already: {target}."
            )
        shutil.rmtree(target)
    target.mkdir(parents=True)
    for name in RECORD_FILES:
        source = fit_dir / name
        if source.exists():
            shutil.copy2(source, target / name)
    return target


# ----------------------------------------------------------------------------
# the uncertainty of a parameter
# ----------------------------------------------------------------------------
def _number(value: float) -> str:
    """Get a number as the shortest text which reads back as the same float."""
    return repr(float(value))


@dataclass(frozen=True, kw_only=True)
class ParameterUncertainty:
    """The fitted value of a parameter and its uncertainty, from the record of a fit.

    All numbers are in the unit of the model (`unit`); the scale is the space
    the fit searched the parameter in.

    Attributes:
        pid: the id of the parameter in the flat model.
        value: the fitted value.
        unit: the unit of the model, as the fit wrote it.
        lower_bound: the lower bound of the fit.
        upper_bound: the upper bound of the fit.
        scale: the space the fit searched the parameter in.
        alpha: the confidence level of the intervals.
        z: the quantile of the standard normal of the confidence level of the
            profile, `sqrt(chi2(alpha, df))`.
        method: where the confidence interval and the distribution come from.
        ci_lower: the lower end of the confidence interval, `None` if it is
            open (the profile reaches the bound).
        ci_upper: the upper end of the confidence interval, `None` if open.
        identifiability: the class of the profile likelihood, `None` without.
        standard_error: the standard error of the Fisher information in the
            space of the fit (`scale`), `None` without.
        fisher_ci: the confidence interval of the Fisher information.
        fit_id: the id of the fit.
        cost: the cost of the fit.
        n_data: the number of data points of the fit, `None` if unknown.
        record: the path of the record, relative to the folder of the fits.
    """

    pid: str
    value: float
    unit: str | None
    lower_bound: float
    upper_bound: float
    scale: ParameterScaleType
    alpha: float
    z: float
    method: UncertaintyMethod
    ci_lower: float | None
    ci_upper: float | None
    identifiability: str | None = None
    standard_error: float | None = None
    fisher_ci: tuple[float, float] | None = None
    fit_id: str
    cost: float
    n_data: int | None = None
    record: str

    def _to_scale(self, x: float) -> float:
        """Get a value in the space of the fit."""
        return float(self.scale.to_scale(x))

    def sigma(self) -> float | None:
        """Get the standard deviation of the distribution in the space of the fit.

        From the profile it is the mean of the closed half widths of the
        confidence interval divided by `z`; from the Fisher information its
        standard error. `None` if neither side of the profile is closed.
        """
        if self.method == UncertaintyMethod.FISHER:
            return self.standard_error
        center = self._to_scale(self.value)
        widths = [
            abs(center - self._to_scale(bound))
            for bound in (self.ci_lower, self.ci_upper)
            if bound is not None
        ]
        widths = [w for w in widths if w > 0.0]
        if not widths:
            return None
        return float(np.mean(widths)) / self.z

    def interval(self) -> tuple[float, float]:
        """Get the confidence interval, an open side at the bound of the fit."""
        lower = self.lower_bound if self.ci_lower is None else self.ci_lower
        upper = self.upper_bound if self.ci_upper is None else self.ci_upper
        return lower, upper

    def standard_error_model(self) -> float | None:
        """Get the standard error in the unit of the model (delta method)."""
        if self.standard_error is None:
            return None
        if self.scale is ParameterScaleType.LOG10:
            return self.value * math.log(10) * self.standard_error
        if self.scale is ParameterScaleType.LOG:
            return self.value * self.standard_error
        return self.standard_error

    def distribution(self) -> str:
        """Get the distribution of the draws as a formula of distrib.

        Returns:
            On a logarithmic scale `lognormal(ln(value), sigma, lb, ub)` with
            sigma in natural logarithms, on the linear scale `normal(value,
            sigma, lb, ub)`, both truncated to the bounds of the fit; without a
            closed side of the profile `uniform(lb, ub)`.
        """
        lb, ub = _number(self.lower_bound), _number(self.upper_bound)
        sigma = self.sigma()
        if sigma is None or not sigma > 0.0:
            return f"uniform({lb}, {ub})"
        if self.scale.is_log:
            sigma_ln = (
                sigma * math.log(10)
                if self.scale is ParameterScaleType.LOG10
                else sigma
            )
            mu = math.log(self.value)
            return f"lognormal({_number(mu)}, {_number(sigma_ln)}, {lb}, {ub})"
        return f"normal({_number(self.value)}, {_number(sigma)}, {lb}, {ub})"

    def notes(self) -> str:
        """Get the provenance of the value and of its uncertainty as markdown."""
        lower, upper = self.interval()
        rows = [
            ("fit", f"`{self.fit_id}`"),
            ("record", f"`{self.record}`"),
            ("cost", f"{self.cost:.6g}"),
            ("data points", "" if self.n_data is None else str(self.n_data)),
            ("value", f"{self.value:.6g} {self.unit or ''}".strip()),
            ("bounds", f"[{self.lower_bound:.6g}, {self.upper_bound:.6g}]"),
            ("scale", self.scale.name),
            ("method", self.method.value),
            ("identifiability", self.identifiability or ""),
            (
                f"{100 * self.alpha:g} % confidence interval",
                f"[{lower:.6g}, {upper:.6g}]"
                + (
                    " (open sides at the bounds)"
                    if self.ci_lower is None or self.ci_upper is None
                    else ""
                ),
            ),
        ]
        if self.fisher_ci is not None:
            rows.append(
                (
                    "Fisher confidence interval",
                    f"[{self.fisher_ci[0]:.6g}, {self.fisher_ci[1]:.6g}]",
                )
            )
        rows.append(("distribution of the draws", f"`{self.distribution()}`"))
        table = "\n".join(f"| {key} | {value} |" for key, value in rows)
        return (
            f"The value of `{self.pid}` is the result of the fit `{self.fit_id}`, "
            f"its uncertainty is the {self.method.value} of the fit.\n\n"
            f"| | |\n|---|---|\n{table}\n"
        )

    def uncertainty(self, unit: UnitType = None) -> Uncertainty:
        """Get the distrib uncertainty of the parameter.

        Args:
            unit: the unit of the values, the unit of the parameter in the
                model definition (`sbmlutils.factory` unit).

        Returns:
            The uncertainty: the confidence interval, the standard error of
            the Fisher information if there is one and the distribution of the
            draws, the provenance in its notes.
        """
        from sbmlutils.factory import Uncertainty, UncertParameter, UncertSpan

        lower, upper = self.interval()
        children: list[UncertParameter | UncertSpan] = [
            UncertSpan(
                type=libsbml.DISTRIB_UNCERTTYPE_CONFIDENCEINTERVAL,
                valueLower=lower,
                valueUpper=upper,
                unit=unit,
                name=f"{100 * self.alpha:g} % confidence interval ({self.method.value})",
            )
        ]
        se = self.standard_error_model()
        if se is not None and math.isfinite(se):
            children.append(
                UncertParameter(
                    type=libsbml.DISTRIB_UNCERTTYPE_STANDARDERROR,
                    value=se,
                    unit=unit,
                    name="standard error (Fisher information)",
                )
            )
        children.append(
            UncertParameter(
                type=libsbml.DISTRIB_UNCERTTYPE_DISTRIBUTION,
                math=self.distribution(),
                name="distribution of the draws",
            )
        )
        return Uncertainty(
            uncertParameters=children,
            name=f"fit {self.fit_id}",
            notes=self.notes(),
        )


# ----------------------------------------------------------------------------
# the records of the fits
# ----------------------------------------------------------------------------
def _read_json(path: Path) -> dict[str, Any]:
    """Read a JSON object."""
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def _convert(value: float, unit: str | None, target: str | None) -> float:
    """Convert a number from one unit into another, both units of pint."""
    if unit is None or target is None or unit == target:
        return float(value)
    return float(ureg.Quantity(value, unit).to(target).magnitude)


def _convert_optional(
    value: float | None, unit: str | None, target: str | None
) -> float | None:
    """Convert a number which may be missing, see `_convert`."""
    return None if value is None else _convert(value, unit, target)


def _scale(name: str | None, default: ParameterScaleType) -> ParameterScaleType:
    if name is None:
        return default
    return ParameterScaleType[str(name)]


def _fisher(path: Path) -> FisherInformation:
    """Read the Fisher information a fit wrote (`FisherInformation.to_dict`)."""
    d = _read_json(path)
    return FisherInformation(
        opid=d["opid"],
        sid=d["sid"],
        pids=list(d["pids"]),
        values=np.asarray(d["values"], dtype=float),
        scale=ParameterScaleType[d["scale"]],
        matrix=np.asarray(d["matrix"], dtype=float),
        cost=float(d["cost"]),
        n=int(d["n"]),
        alpha=float(d.get("alpha", 0.95)),
        units=list(d.get("units") or []),
        scales=[ParameterScaleType[s] for s in d.get("scales") or []],
        targets=list(d.get("targets") or []),
    )


@dataclass(frozen=True)
class FitRecord:
    """The record of an adopted fit: the fitted parameters with their uncertainty.

    Attributes:
        directory: the directory of the record, named by the id of the fit.
        parameters: the uncertainty of every fitted parameter by pid.
    """

    directory: Path
    parameters: Mapping[str, ParameterUncertainty] = field(repr=False)

    @property
    def fit_id(self) -> str:
        """The id of the fit."""
        return self.directory.name

    @classmethod
    def load(cls, directory: Path) -> FitRecord:
        """Read the record of an adopted fit.

        The value of a parameter is the best parameter set of the fit
        (`parameters.json`), the bounds and the scale are those of the problem
        (`optimization_result.json`); the confidence interval is the one of the
        profile likelihood (`identifiability.json`), the Fisher information
        (`fisher.json`) gives the standard error and is the fallback for a
        parameter without a profile.

        Args:
            directory: the directory of the record.

        Returns:
            The record.

        Raises:
            FileNotFoundError: if the record lacks a required file or has
                neither a profile nor a Fisher information.
            ValueError: if a fitted parameter has neither a profile nor a
                Fisher information.
        """
        fit_id = directory.name
        missing = [n for n in REQUIRED_FILES if not (directory / n).exists()]
        if missing:
            raise FileNotFoundError(f"The record '{directory}' has no {missing}.")
        best = _read_json(directory / "parameters.json")["sets"][0]
        values: dict[str, float] = {k: float(v) for k, v in best["values"].items()}
        units: dict[str, str | None] = dict(best.get("units") or {})
        result = _read_json(directory / "optimization_result.json")
        fit_parameters = {p["pid"]: p for p in result["parameters"]}
        default_scale = _scale(
            result.get("settings", {}).get("parameter_scale"),
            ParameterScaleType.LINEAR,
        )

        profiles: dict[str, Mapping[str, Any]] = {}
        alpha, df = 0.95, 1
        if (directory / "identifiability.json").exists():
            identifiability = _read_json(directory / "identifiability.json")
            profiles = identifiability.get("profiles") or {}
            settings = identifiability.get("settings") or {}
            alpha = float(settings.get("alpha", alpha))
            df = int(settings.get("degrees_of_freedom", df))
        fisher: FisherInformation | None = None
        if (directory / "fisher.json").exists():
            fisher = _fisher(directory / "fisher.json")
        if not profiles and fisher is None:
            raise FileNotFoundError(
                f"The record '{directory}' has neither a profile likelihood nor a "
                f"Fisher information."
            )
        z = math.sqrt(float(stats.chi2.ppf(alpha, df)))

        fisher_se: dict[str, float] = {}
        fisher_ci: dict[str, tuple[float, float]] = {}
        fisher_alpha = alpha
        if fisher is not None:
            fisher_alpha = fisher.alpha
            errors = fisher.standard_errors
            lower, upper = fisher.confidence_intervals()
            for k, pid in enumerate(fisher.pids):
                if math.isfinite(float(errors[k])):
                    fisher_se[pid] = float(errors[k])
                    fisher_ci[pid] = (float(lower[k]), float(upper[k]))

        parameters: dict[str, ParameterUncertainty] = {}
        for pid, value in values.items():
            fp = fit_parameters[pid]
            unit = units.get(pid)
            fp_unit = fp.get("unit")
            scale = _scale(fp.get("scale"), default_scale)
            lower_bound = _convert(float(fp["lower_bound"]), fp_unit, unit)
            upper_bound = _convert(float(fp["upper_bound"]), fp_unit, unit)
            profile = profiles.get(pid)
            if profile is not None:
                method = UncertaintyMethod.PROFILE
                ci_lower = _convert_optional(profile.get("ci_lower"), fp_unit, unit)
                ci_upper = _convert_optional(profile.get("ci_upper"), fp_unit, unit)
                identifiability_class = profile.get("identifiability")
                level = alpha
            elif pid in fisher_se:
                method = UncertaintyMethod.FISHER
                ci_lower, ci_upper = fisher_ci[pid]
                identifiability_class = None
                level = fisher_alpha
            else:
                raise ValueError(
                    f"The parameter '{pid}' of the fit '{fit_id}' has neither a "
                    f"profile likelihood nor a Fisher information."
                )
            parameters[pid] = ParameterUncertainty(
                pid=pid,
                value=value,
                unit=unit,
                lower_bound=lower_bound,
                upper_bound=upper_bound,
                scale=scale,
                alpha=level,
                z=z,
                method=method,
                ci_lower=ci_lower,
                ci_upper=ci_upper,
                identifiability=identifiability_class,
                standard_error=fisher_se.get(pid),
                fisher_ci=fisher_ci.get(pid),
                fit_id=fit_id,
                cost=float(best["cost"]),
                n_data=None if fisher is None else fisher.n,
                record=f"{fit_id}/",
            )
        return cls(directory=directory, parameters=parameters)


class FitRecords:
    """The adopted fits of a model: the fitted values and their uncertainty for the factory.

    A model file takes the value of a fitted parameter from the record of its
    fit, so the value, its uncertainty and its provenance come from one place:

    ```python
    # not executed
    FITS = MODEL.fit_records()
    Parameter("HCTZEX_f", **FITS.parameter("KI__HCTZEX_f", unit=U.dimensionless))
    ```
    """

    def __init__(self, directory: Path, fit_ids: Sequence[str]) -> None:
        """Read the records of the adopted fits.

        Args:
            directory: the folder of the adopted fits.
            fit_ids: the ids of the adopted fits.

        Raises:
            FileNotFoundError: if a record does not exist.
            ValueError: if a parameter is in two fits.
        """
        self.directory = directory
        self.records: dict[str, FitRecord] = {}
        for fit_id in fit_ids:
            path = directory / fit_id
            if not path.is_dir():
                raise FileNotFoundError(
                    f"The fit '{fit_id}' is not adopted: {path} does not exist, "
                    f"adopt it with `run_fit <package> --adopt {fit_id}`."
                )
            self.records[fit_id] = FitRecord.load(path)
        self.parameters: dict[str, ParameterUncertainty] = {}
        for record in self.records.values():
            for pid, parameter in record.parameters.items():
                if pid in self.parameters:
                    raise ValueError(
                        f"The parameter '{pid}' is a result of the fits "
                        f"'{self.parameters[pid].fit_id}' and '{record.fit_id}'; "
                        f"a parameter takes its value from one adopted fit."
                    )
                self.parameters[pid] = parameter

    def __getitem__(self, pid: str) -> ParameterUncertainty:
        """Get the uncertainty of a fitted parameter.

        Raises:
            KeyError: if no adopted fit has the parameter.
        """
        try:
            return self.parameters[pid]
        except KeyError:
            raise KeyError(
                f"No adopted fit has the parameter '{pid}', the fitted parameters "
                f"are {sorted(self.parameters)}."
            ) from None

    def parameter(self, pid: str, unit: UnitType = None) -> dict[str, Any]:
        """Get the keyword arguments of a `sbmlutils.factory.Parameter` of a fitted parameter.

        Args:
            pid: the id of the parameter in the flat model, as the fit names it.
            unit: the unit of the parameter in the model definition, set on the
                parameter and its uncertainty.

        Returns:
            `value` and `uncertainties`, and `unit` if given.

        Raises:
            KeyError: if no adopted fit has the parameter.
        """
        parameter = self[pid]
        kwargs: dict[str, Any] = {
            "value": parameter.value,
            "uncertainties": [parameter.uncertainty(unit)],
        }
        if unit is not None:
            kwargs["unit"] = unit
        return kwargs


# ----------------------------------------------------------------------------
# the distributions of a model
# ----------------------------------------------------------------------------
_CALL = re.compile(r"^\s*(\w+)\s*\((.*)\)\s*$")


def distribution_from_formula(formula: str) -> Distribution:
    """Get the sbmlsim distribution of a formula of distrib.

    Args:
        formula: `lognormal(mu, sigma[, lower, upper])` (mu and sigma of the
            natural logarithm), `normal(mean, sd[, lower, upper])` or
            `uniform(lower, upper)`, with numbers as arguments.

    Returns:
        The distribution, truncated if the formula has bounds.

    Raises:
        ValueError: for another formula.
    """
    match = _CALL.match(formula)
    if match is None:
        raise ValueError(f"The distribution '{formula}' is no function of numbers.")
    name = match.group(1)
    try:
        args = [float(a) for a in match.group(2).split(",")]
    except ValueError:
        raise ValueError(
            f"The distribution '{formula}' has arguments which are no numbers."
        ) from None
    inner: Distribution
    if name == "uniform" and len(args) == 2:
        return Uniform(args[0], args[1])
    if name == "lognormal" and len(args) in (2, 4):
        mu, sigma = args[0], args[1]
        inner = LogNormal(math.exp(mu), math.sqrt(math.expm1(sigma**2)))
    elif name == "normal" and len(args) in (2, 4):
        inner = Normal(args[0], args[1])
    else:
        raise ValueError(
            f"The distribution '{formula}' is not supported: lognormal and normal "
            f"with two or four arguments, uniform with two."
        )
    if len(args) == 4:
        return Truncated(inner, args[2], args[3])
    return inner


def parameter_distributions(sbml_path: Path) -> dict[str, Distribution]:
    """Get the distributions of the parameters of a model.

    Every parameter with a distrib uncertainty which states a distribution
    (an uncert parameter of the type `distribution`) has one.

    Args:
        sbml_path: the SBML model, the flat model a simulation runs.

    Returns:
        The distribution of every such parameter by its id.

    Raises:
        ValueError: if a distribution is not supported, see
            `distribution_from_formula`, or a parameter states two.
    """
    doc = libsbml.readSBMLFromFile(str(sbml_path))
    model = doc.getModel()
    if model is None:
        raise ValueError(f"'{sbml_path}' has no SBML model.")
    distributions: dict[str, Distribution] = {}
    for parameter in model.getListOfParameters():
        plugin = parameter.getPlugin("distrib")
        if plugin is None:
            continue
        for k in range(plugin.getNumUncertainties()):
            uncertainty = plugin.getUncertainty(k)
            for j in range(uncertainty.getNumUncertParameters()):
                child = uncertainty.getUncertParameter(j)
                if child.getType() != libsbml.DISTRIB_UNCERTTYPE_DISTRIBUTION:
                    continue
                if not child.isSetMath():
                    continue
                pid = parameter.getId()
                if pid in distributions:
                    raise ValueError(f"The parameter '{pid}' states two distributions.")
                formula = libsbml.formulaToL3String(child.getMath())
                distributions[pid] = distribution_from_formula(formula)
    return distributions


def uncertainty_dimension(sbml_path: Path, draws: int, seed: int) -> Dimension:
    """Draw the parameters of a model which have a distribution.

    Args:
        sbml_path: the SBML model.
        draws: the number of draws.
        seed: the seed of the Latin hypercube.

    Returns:
        The dimension `uncertainty` of the draws.

    Raises:
        ValueError: if no parameter of the model has a distribution.
    """
    distributions = parameter_distributions(sbml_path)
    if not distributions:
        raise ValueError(
            f"No parameter of '{sbml_path}' has a distribution: the uncertainty "
            f"of the fitted parameters comes from the adopted fits of the model."
        )
    return sampling.lhs(distributions, draws, seed=seed, id=UNCERTAINTY_DIMENSION)


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

        Args:
            model: the model, whose `model_path` has the distributions.
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
#: curves draw whose x is not the time
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


def apply_uncertainty(experiment: SimulationExperiment, run: UncertaintyRun) -> None:
    """Run an initialized experiment over the draws of the fitted parameters.

    Every simulation becomes a scan over the draws (a scan gets the dimension
    of the draws appended); a parameter which the simulation sets itself, or
    which the model of the simulation does not have, is not drawn, and a
    simulation without a drawn parameter runs once. Every curve of a figure
    over the time of a task over the draws becomes the band of the quantiles
    of the draws with the median, with the name and the colour of the curve.
    A curve whose x is another variable of the task (e.g. an effect over a
    concentration) has no band across the draws: it draws the simulation with
    the fitted values, which runs once next to the draws (`<simulation>_fitted`).
    The curves of the data stay.

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
            for curve in plot.curves:
                task_id = None if curve.y.task_id is None else str(curve.y.task_id)
                if task_id is None or task_id not in drawn:
                    curves.append(curve)
                    continue
                if curve.x.index != "time":
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
            plot.curves = curves


# ----------------------------------------------------------------------------
# the check of a model
# ----------------------------------------------------------------------------
def check_fit_parameters(model: PKPDModel) -> list[str]:
    """Check that the model has the fitted values and the uncertainty of its adopted fits.

    Every parameter of an adopted fit must be a parameter of the flat model
    (`model_path`) with the value of the fit and a distribution of its draws.

    Args:
        model: the model.

    Returns:
        The problems, one line each; empty for a model without adopted fits.
    """
    if not model.fits:
        return []
    try:
        records = model.fit_records()
    except (FileNotFoundError, ValueError) as err:
        return [str(err)]
    if not model.model_path.exists():
        return [f"the model {model.model_path} does not exist"]
    doc = libsbml.readSBMLFromFile(str(model.model_path))
    sbml_model = doc.getModel()
    distributions = parameter_distributions(model.model_path)
    problems: list[str] = []
    for pid, parameter in records.parameters.items():
        p = sbml_model.getParameter(pid)
        if p is None:
            problems.append(
                f"the parameter '{pid}' of the fit '{parameter.fit_id}' is no "
                f"parameter of {model.model_path.name}"
            )
            continue
        if not math.isclose(p.getValue(), parameter.value, rel_tol=VALUE_RTOL):
            problems.append(
                f"the parameter '{pid}' is {p.getValue()!r} in the model, the fit "
                f"'{parameter.fit_id}' gives {parameter.value!r}: take it from the "
                f"record (`FitRecords.parameter`)"
            )
        if pid not in distributions:
            problems.append(
                f"the parameter '{pid}' of the fit '{parameter.fit_id}' has no "
                f"distribution in the model"
            )
    return problems
