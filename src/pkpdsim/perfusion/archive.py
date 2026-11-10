"""The factory of a machine perfusion model in one call: the files and the COMBINE archive.

A model package builds its `PerfusionDefinition` with
`pkpdsim.perfusion.builder.build` and writes next to the SBML the
differential equations as markdown, the antimony of the model, the figures of
the circuit, the two license texts and the COMBINE archive. `create_archive`
does all of it, as `pkpdsim.pbpk.archive.create_archive` does for a PBPK
model: it builds the definition, writes the schematic
(`pkpdsim.perfusion.figure.circuit_figure`) and the organ in its circuit
(`pkpdsim.perfusion.anatomy.circuit_anatomy_figure`), and the rest with the
archive writer of sbmlsim, `sbmlsim.release.archive.write_archive`, which imports
`sbmlode` and `pymetadata` when it writes. The model of a liver is composed of
the cell model of its zones, so the archive ships the three SBML files the
build writes: the comp model, the cell model, which stays loadable alone (it is
the model a mesh or another architecture couples to), and the flat model, the
master entry, which is the model that is simulated. The model of an organ
without zones is its own flat model and the only entry.
"""

from pathlib import Path
from typing import Any

from sbmlsim.release.archive import ArchiveModel, write_archive

from pkpdsim.perfusion.anatomy import circuit_anatomy_figure
from pkpdsim.perfusion.builder import build, zone_cell_model
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

    Builds the definition into `model_output_dir` (`<sid>.xml`, and for a
    liver, whose zones are submodels, the cell model `<sid>_core.xml` and the
    flat model `<sid>_flat.xml`), writes the
    schematic of the circuit there (`<sid>.svg`, panel B), the anatomical
    figure of the organ in its circuit (`<sid>_anatomy.svg`, panel A), the
    ODEs of the flat model as markdown (`<sid>.md`, `<sid>_flat.md` for a liver,
    `sbmlode`) and its antimony (`<sid>.ant`, `<sid>_flat.ant`), the texts of CC BY-SA 4.0 and MIT (with the copyright)
    into the parent of the directory and the archive `<model>_model.omex`
    there, with the README, the licenses, every model the build writes (the
    flat model the master entry, and for a liver the comp model and the cell
    model, which stays loadable alone), the
    markdown, the antimony and the two figures (the schematic first). Every
    file of the archive lies below the parent of `model_output_dir` and is
    located at its path relative to it (`./models/<sid>.xml` for the
    directory `models` of a model package, the licenses `./mit.txt`); the
    README alone is placed as `./README.md`.

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
        a README), `cc-by`, `mit`, `<sid>` of the model (the comp model for a
        liver), `<sid>_md` of its
        ODEs, `<sid>_ant` of its antimony (for a liver `<sid>_flat` of the flat
        model, `<sid>_flat_md`, `<sid>_flat_ant` and `<sid>_core` of the cell
        model), `<sid>_svg` of the schematic and
        `<sid>_anatomy_svg` of the anatomical figure.

    Raises:
        ValueError: for a definition which does not build (`build`).
    """
    model_output_dir.mkdir(parents=True, exist_ok=True)
    path = build(definition, model_output_dir)
    sid = definition.sid
    cell = zone_cell_model(definition)
    # the flat model is the master entry: the comp model of a liver carries it
    models: list[Path | ArchiveModel] = (
        [path]
        if cell is None
        else [ArchiveModel(model_output_dir / f"{sid}.xml", flat=path)]
    )
    # the cell model of the zones, which is loadable and simulates alone
    cells = [] if cell is None else [model_output_dir / f"{cell.sid}.xml"]
    figures = [
        circuit_figure(definition, model_output_dir / f"{sid}.svg"),
        circuit_anatomy_figure(definition, model_output_dir / f"{sid}_anatomy.svg"),
    ]
    return write_archive(
        models,
        model_output_dir.parent,
        name=f"{model}_model",
        copyright=copyright,
        readme=readme,
        extra_models=cells,
        extra_files=figures,
    )
