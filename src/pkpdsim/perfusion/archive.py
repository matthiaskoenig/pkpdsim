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
without zones is its own flat model and the only entry. The interactive panels
of the model (`sbmlsim.panel.Panel`, `create_archive(..., panels=...)`) are an
annotation of that model, written before its ODEs, its antimony and the
archive are.
"""

from collections.abc import Sequence
from pathlib import Path
from typing import Any

from sbmlsim.panel import Panel
from sbmlsim.panel_sbml import write_panels
from sbmlsim.release.archive import ArchiveModel, write_archive

from pkpdsim.pbpk.archive import check_panels
from pkpdsim.perfusion.anatomy import circuit_anatomy_figure
from pkpdsim.perfusion.builder import PerfusionExtension, build, zone_cell_model
from pkpdsim.perfusion.definition import PerfusionDefinition
from pkpdsim.perfusion.figure import circuit_figure, interactive_circuit

__all__ = ["create_archive"]


def create_archive(
    definition: PerfusionDefinition,
    model_output_dir: Path,
    *,
    model: str,
    copyright: str,
    readme: Path | None = None,
    extension: PerfusionExtension | None = None,
    panels: Sequence[Panel] | None = None,
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
    markdown, the antimony and the two figures (the schematic first). The
    interactive snippet of each figure is written beside it, `<sid>.html` and
    `<sid>_anatomy.html` (`pkpdsim.perfusion.figure.interactive_circuit`),
    which a page of the documentation of the model embeds; the snippets are no
    entry of the archive, which carries the figures as SVG. Every
    file of the archive lies below the parent of `model_output_dir` and is
    located at its path relative to it (`./models/<sid>.xml` for the
    directory `models` of a model package, the licenses `./mit.txt`); the
    README alone is placed as `./README.md`. With `panels` the interactive
    panels of the model are written into the flat model, the master entry
    (`sbmlsim.panel_sbml.write_panels`), before its markdown, its antimony and
    the archive are written.

    Args:
        definition: the definition of the circuit.
        model_output_dir: directory of the model; the licenses and the
            archive are written into its parent.
        model: name of the model, the archive is `<model>_model.omex`.
        copyright: copyright of the MIT license, e.g. `2026 Matthias König`.
        readme: the README of the model repository, left out of the archive if
            `None`.
        extension: a mechanism of a package, in the cell model of the zones
            and in the composed model (`pkpdsim.perfusion.PerfusionExtension`).
        panels: the interactive panels of the model, the annotation of its
            flat model; `None` or none annotates nothing.

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
        ValueError: for a definition which does not build (`build`), for panels
            which repeat an id (before anything is written) and for panels
            which do not read back from the flat model as they were declared
            (`write_panels`).
        TypeError: for panels which are no `sbmlsim.panel.Panel` (before
            anything is written).
    """
    annotation = check_panels(panels or (), f"the model '{definition.sid}'")
    model_output_dir.mkdir(parents=True, exist_ok=True)
    path = build(definition, model_output_dir, extension)
    if annotation:
        write_panels(path, path, annotation)
    sid = definition.sid
    cell = zone_cell_model(definition, extension)
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
    # the interactive snippets beside the figures, for the documentation
    for svg, kind in zip(figures, ("schematic", "anatomy"), strict=True):
        interactive_circuit(
            definition, svg.with_suffix(".html"), figure=kind, snippet=True
        )
    return write_archive(
        models,
        model_output_dir.parent,
        name=f"{model}_model",
        copyright=copyright,
        readme=readme,
        extra_models=cells,
        extra_files=figures,
    )
