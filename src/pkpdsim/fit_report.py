"""The report of an adopted fit, recreated from its record.

The record of an adopted fit (`pkpdsim.uncertainty.adopt_fit`) holds only what
is needed to recreate the report of the fit: the parameters
(`parameters.json`), the profile likelihood (`identifiability.json`) and the
Fisher information (`fisher.json`). `report_adopted_fit` builds the fit problem
of the record again from the fit definitions of the model and writes the full
report with its figures, its tables and its identifiability section
(`run_report <package> --adopted <fit id>`).
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import TYPE_CHECKING

from sbmlsim.fit.identifiability import IdentifiabilityResult
from sbmlsim.fit.parameters import ParameterSets
from sbmlsim.fit.report import FitReport

from pkpdsim.uncertainty import read_fisher

if TYPE_CHECKING:
    from pkpdsim.model import PKPDModel

FIT_ID = re.compile(r"^(?P<subset>.+)_\d{8}_\d{6}__[0-9a-f]+$")
"""The id of a fit: the fit problem, the time and a short hash
(`sbmlsim.fit.result.fit_id`), e.g. `PK_20261016_185018__3efe3`."""


def fit_subset(fit_id: str) -> str:
    """Get the fit problem of a fit from its id.

    Args:
        fit_id: the id of the fit, e.g. `PK_20261016_185018__3efe3`.

    Returns:
        The name of the fit problem, the `--subset` of the fit, e.g. `PK`.

    Raises:
        ValueError: if the text is no id of a fit.
    """
    match = FIT_ID.match(fit_id)
    if match is None:
        raise ValueError(
            f"'{fit_id}' is no id of a fit, i.e. <problem>_<YYYYMMDD>_<HHMMSS>__<hash>."
        )
    return match["subset"]


def fit_record_dir(model: PKPDModel, fit_id: str) -> Path:
    """Get the record of an adopted fit of a model.

    Args:
        model: the model.
        fit_id: the id of the adopted fit.

    Returns:
        The directory of the record in the adopted fits of the model.

    Raises:
        ValueError: if the text is no id of a fit or its fit problem is no fit
            definition of the model.
        FileNotFoundError: if the model has no record of the fit or the record
            has no `parameters.json`.
    """
    subset = fit_subset(fit_id)
    definitions = model.load_fit_definitions()
    if subset not in definitions:
        raise ValueError(
            f"The fit '{fit_id}' is of the problem '{subset}', which is no fit "
            f"definition of the model '{model.name}': {sorted(definitions)}."
        )
    directory = model.fits_dir / fit_id
    if not directory.is_dir():
        raise FileNotFoundError(
            f"The model '{model.name}' has no record of the fit '{fit_id}' in "
            f"{model.fits_dir}."
        )
    if not (directory / "parameters.json").exists():
        raise FileNotFoundError(f"The record '{directory}' has no parameters.json.")
    return directory


def report_adopted_fit(model: PKPDModel, fit_id: str, output_dir: Path) -> Path:
    """Create the report of an adopted fit from its record.

    The fit problem is the fit definition of the model named by the id of the
    fit; the parameter sets come from `parameters.json`, the identifiability
    section from `identifiability.json` and `fisher.json` of the record, each
    when present. The report is written with all of its figures.

    Args:
        model: the model.
        fit_id: the id of the adopted fit, its record is in `model.fits_dir`.
        output_dir: the directory of the reports, the report is written into
            `output_dir / fit_id`.

    Returns:
        The directory of the report.

    Raises:
        ValueError: if the text is no id of a fit or its fit problem is no fit
            definition of the model.
        FileNotFoundError: if the model has no record of the fit or the record
            has no `parameters.json`.
    """
    directory = fit_record_dir(model, fit_id)
    definition = model.load_fit_definitions()[fit_subset(fit_id)]
    identifiability_path = directory / "identifiability.json"
    fisher_path = directory / "fisher.json"
    report = FitReport(
        problem=definition.problem(opid=fit_id),
        settings=definition.settings,
        parameter_sets=ParameterSets.from_json(directory / "parameters.json"),
        identifiability=(
            IdentifiabilityResult.from_json(identifiability_path)
            if identifiability_path.exists()
            else None
        ),
        fisher=read_fisher(fisher_path) if fisher_path.exists() else None,
        show_titles=False,
    )
    report.create(output_dir=output_dir, name=fit_id)
    return output_dir / fit_id
