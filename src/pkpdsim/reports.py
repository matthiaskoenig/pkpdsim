"""The reports of the adopted fits of a model.

A fit is adopted with `run_fit <package> --adopt <fit id>`, its record
(`PKPDModel.fits_dir`) holds everything its report is built from. The report of
every adopted fit is written separately, there is no report which merges the
fits of a model: the fits have different problems, data and parameters.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from sbmlsim.fit.cli import load_fit_records, load_parameter_sets
from sbmlsim.fit.report import FitReport

from pkpdsim.fit import fit_subset, report_context

if TYPE_CHECKING:
    from pathlib import Path

    from pkpdsim.model import PKPDModel


def adopted_reports(
    model: PKPDModel, output_dir: Path | None = None, *, mapping_figures: bool = True
) -> list[Path]:
    """Write the report of every adopted fit of a model, each from its record.

    No fit is run. The report of a record is the report which
    `run_report <package> --subset <subset> <record>/parameters.json` writes
    for it, with the context of the model (`report_context`): the title, the
    description of the fit and of its parameters, and the adopted fits whose
    parameters were fixed. A record without provenance gives a report which
    states that the provenance was not recorded.

    Args:
        model: the model, its `fits` are the fits to report.
        output_dir: the directory of the reports, `model.report_dir` by default;
            every report is written into the folder named by its fit id.
        mapping_figures: create the figures of the fit mappings.

    Returns:
        The directories of the reports, in the order of `model.fits`.

    Raises:
        ValueError: if the model has no adopted fits, or the subset of an
            adopted fit has no fit definition.
        FileNotFoundError: if the record of an adopted fit has no
            `parameters.json`.
    """
    if not model.fits:
        raise ValueError(f"The model '{model.name}' has no adopted fits.")
    definitions = model.load_fit_definitions()
    directories: list[Path] = []
    for fit_id in model.fits:
        subset = fit_subset(fit_id)
        if subset not in definitions:
            raise ValueError(
                f"The subset '{subset}' of the adopted fit '{fit_id}' has no fit "
                f"definition, the subsets are {', '.join(definitions)}."
            )
        definition = definitions[subset]
        record = model.fits_dir / fit_id
        parameters = record / "parameters.json"
        if not parameters.exists():
            raise FileNotFoundError(
                f"The record of the adopted fit '{fit_id}' has no {parameters}."
            )
        records = load_fit_records(record)
        # the settings the fit ran with, the definition may have changed since
        opt_result = records.opt_result
        settings = (
            opt_result.settings
            if opt_result is not None and opt_result.settings is not None
            else definition.settings
        )
        problem = definition.problem(opid=fit_id)
        problem.initialize(settings)
        report = FitReport(
            problem=problem,
            settings=settings,
            parameter_sets=load_parameter_sets([parameters]),
            opt_result=records.opt_result,
            identifiability=records.identifiability,
            fisher=records.fisher,
            show_titles=False,
            mapping_figures=mapping_figures,
            context=report_context(model, subset, definition),
        )
        directories.append(
            report.create(output_dir=output_dir or model.report_dir, name=fit_id)
        )
    return directories
