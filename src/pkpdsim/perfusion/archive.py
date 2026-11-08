"""The factory of a machine perfusion model in one call: the files and the COMBINE archive.

A model package builds its `PerfusionDefinition` with
`pkpdsim.perfusion.builder.build` and writes next to the SBML the
differential equations as markdown, the antimony of the model, the figure of
the circuit, the two license texts and the COMBINE archive. `create_archive`
does all of it, as `pkpdsim.pbpk.archive.create_archive` does for a PBPK
model: it builds the definition, writes the figure
(`pkpdsim.perfusion.figure.circuit_figure`) and the rest with the archive
writer of sbmlsim, `sbmlsim.release.archive.write_archive`, which imports
`sbmlode` and `pymetadata` when it writes. A perfusion model has no
submodels, so its SBML is its own flat model, the master entry of the
archive.
"""

from pathlib import Path
from typing import Any

from sbmlsim.release.archive import write_archive

from pkpdsim.perfusion.builder import build
from pkpdsim.perfusion.definition import PerfusionDefinition
from pkpdsim.perfusion.figure import circuit_figure

__all__ = ["create_archive"]


def create_archive(
    definition: PerfusionDefinition,
    model_output_dir: Path,
    *,
    model: str,
    copyright: str,
    readme: Path | None = None,
) -> dict[str, dict[str, Any]]:
    """Build a definition and write its files and the COMBINE archive.

    Builds the definition into `model_output_dir` (`<sid>.xml`), writes the
    figure of the circuit there (`<sid>.svg`), the ODEs of the model as
    markdown (`<sid>.md`, `sbmlode`) and its antimony (`<sid>.ant`), the
    texts of CC BY-SA 4.0 and MIT (with the copyright) into the parent of the
    directory and the archive `<model>_model.omex` there, with the README,
    the licenses, the model (the master entry), the markdown, the antimony and
    the figure. Every file of the archive lies below the parent of
    `model_output_dir` and is located at its path relative to it
    (`./models/<sid>.xml` for the directory `models` of a model package, the
    licenses `./mit.txt`); the README alone is placed as `./README.md`.

    Args:
        definition: the definition of the circuit.
        model_output_dir: directory of the model; the licenses and the
            archive are written into its parent.
        model: name of the model, the archive is `<model>_model.omex`.
        copyright: copyright of the MIT license, e.g. `2026 Matthias König`.
        readme: the README of the model repository, left out of the archive if
            `None`.

    Returns:
        The files with their entries in the COMBINE archive, by the keys of a
        factory, each `{"path": Path, "entry": ManifestEntry}`: `README` (with
        a README), `cc-by`, `mit`, `<sid>` of the model, `<sid>_md` of its
        ODEs, `<sid>_ant` of its antimony and `<sid>_svg` of the figure.

    Raises:
        ValueError: for a definition which does not build (`build`).
    """
    model_output_dir.mkdir(parents=True, exist_ok=True)
    path = build(definition, model_output_dir)
    figure = circuit_figure(definition, model_output_dir / f"{definition.sid}.svg")
    return write_archive(
        [path],
        model_output_dir.parent,
        name=f"{model}_model",
        copyright=copyright,
        readme=readme,
        extra_files=[figure],
    )
