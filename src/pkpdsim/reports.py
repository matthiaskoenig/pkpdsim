"""The reports of the adopted fits of a model.

A fit is adopted with `run_fit <package> --adopt <fit id>`, its record
(`PKPDModel.fits_dir`) holds everything its report is built from. The report of
every adopted fit is written separately, there is no report which merges the
fits of a model: the fits have different problems, data and parameters.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from sbmlsim.fit.cli import load_fit_records, report_from_records

from pkpdsim.fit import fit_subset, records_context

if TYPE_CHECKING:
    from pathlib import Path

    from pkpdsim.model import PKPDModel


def adopted_record(model: PKPDModel, fit_id: str) -> Path:
    """Get the record of an adopted fit of a model, which its report is built from.

    Args:
        model: the model.
        fit_id: the id of the adopted fit, e.g. `PK_20261016_185018__3efe3`.

    Returns:
        The directory of the record in `model.fits_dir`.

    Raises:
        ValueError: if the id is not the id of a fit (`pkpdsim.fit.fit_subset`)
            or its subset has no fit definition.
        FileNotFoundError: if the record of the fit has no `parameters.json`.
    """
    subset = fit_subset(fit_id)
    definitions = model.load_fit_definitions()
    if subset not in definitions:
        raise ValueError(
            f"The subset '{subset}' of the adopted fit '{fit_id}' has no fit "
            f"definition, the subsets are {', '.join(definitions)}."
        )
    record = model.fits_dir / fit_id
    parameters = record / "parameters.json"
    if not parameters.exists():
        raise FileNotFoundError(
            f"The record of the adopted fit '{fit_id}' has no {parameters}."
        )
    return record


def adopted_reports(
    model: PKPDModel, output_dir: Path | None = None, *, mapping_figures: bool = True
) -> list[Path]:
    """Write the report of every adopted fit of a model, each from its record.

    No fit is run. The report of a record is the report which
    `run_report <package> --subset <subset> --name <fit id>
    <record>/parameters.json` writes for it: the same problem (its id is the
    subset), the settings the fit ran with, the records of the fit
    (`sbmlsim.fit.cli.report_from_records`) and the context of the model
    (`pkpdsim.fit.records_context`), i.e. the title, the description of the
    fit and of its parameters, and the adopted fits whose parameters were
    fixed as the provenance of the record names them. A record without
    provenance gives a report which states that the provenance was not
    recorded and names no fixed fits.

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
    context = records_context(model)
    directories: list[Path] = []
    for fit_id in model.fits:
        record = adopted_record(model, fit_id)
        subset = fit_subset(fit_id)
        definition = definitions[subset]
        records = load_fit_records(record)
        report = report_from_records(
            definition,
            record,
            opid=subset,
            records=records,
            context=context(subset, definition, records),
            mapping_figures=mapping_figures,
        )
        directories.append(
            report.create(output_dir=output_dir or model.report_dir, name=fit_id)
        )
    return directories
