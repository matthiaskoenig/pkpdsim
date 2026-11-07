"""The factory of a library model in one call: the files and the COMBINE archive.

A model package builds its `PBPKDefinition` with `pkpdsim.pbpk.build` and
writes next to the SBML the differential equations as markdown, the antimony
of the flat model, the two license texts and the COMBINE archive.
`create_archive` does all of it, so the factory of a model is one call: it
builds the definition and writes the rest with the archive writer of sbmlsim,
`sbmlsim.release.archive.write_archive`, which imports `sbmlode` and
`pymetadata` when it writes.

A model with a build per sex builds the one definition for every sex
(`sexes`): the build of the sex of the definition keeps its id (`<sid>.xml`,
`<sid>_flat.xml`, the master entry of the archive) and the build of another
sex is the definition with that sex and the id `<sid>_<sex>`, e.g.
`<sid>_female_flat.xml` next to the male `<sid>_flat.xml`. The model object
holds both (`PKPDModel(model_path=..., female_model_path=...)`); its `model_path`
is the male build, so a model with builds per sex keeps its definition male
(the default), whose build is the master entry.
"""

import dataclasses
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from sbmlsim.release.archive import ArchiveModel, write_archive

from pkpdsim.pbpk.builder import build
from pkpdsim.pbpk.definition import PBPKDefinition
from pkpdsim.pbpk.figure import model_figure
from pkpdsim.sex import Sex


def sex_definition(definition: PBPKDefinition, sex: Sex) -> PBPKDefinition:
    """Get the definition of the build of a sex.

    Args:
        definition: the definition of the model.
        sex: the sex of the build.

    Returns:
        The definition itself for its own sex, else the definition with the
        sex, the id `<sid>_<sex>`, e.g. `hctz_female`, and the name with the
        sex, `<name> (female)`.
    """
    sex = Sex(sex)
    if sex is definition.sex:
        return definition
    name = None if definition.name is None else f"{definition.name} ({sex})"
    return dataclasses.replace(
        definition, sex=sex, sid=f"{definition.sid}_{sex}", name=name
    )


def create_archive(
    definition: PBPKDefinition,
    model_output_dir: Path,
    *,
    model: str,
    copyright: str,
    readme: Path | None = None,
    extra_models: Sequence[Path] = (),
    extra_files: Sequence[Path] = (),
    sexes: Sequence[Sex] | None = None,
    figure: bool = True,
) -> dict[str, dict[str, Any]]:
    """Build a definition and write its files and the COMBINE archive.

    Builds the definition into `model_output_dir` (`build`), writes the ODEs of
    the flat model as markdown (`sbmlode`) and its antimony, the texts of CC
    BY-SA 4.0 and MIT (with the copyright) into the parent of the directory and
    the archive `<model>_model.omex` there, with the README, the licenses, the
    comp model, the flat model (the master entry), the markdown, the antimony,
    the extra models and the extra files. With `sexes` every sex is built
    (`sex_definition`), each with its comp model, flat model, markdown and
    antimony in the archive; the flat model of the sex of the definition stays
    the master entry, the first model of the archive. With `figure` the figure of the definition
    (`model_figure`, an SVG schematic of the organs, routes, substances and
    processes) is written into `model_output_dir` as `<sid>.svg` and is an
    entry of the archive.

    Args:
        definition: the definition of the model.
        model_output_dir: directory of the models; the licenses and the archive
            are written into its parent.
        model: name of the model, the archive is `<model>_model.omex`.
        copyright: copyright of the MIT license, e.g. `2024-2026 Matthias König`.
        readme: the README of the model repository, left out of the archive if
            `None`.
        extra_models: further SBML files (in L3V2) which are entries of the
            archive, e.g. a submodel which is simulated on its own.
        extra_files: further files of other formats which are entries of the
            archive, e.g. the ODEs of an extra model (`.md`) or a figure
            (`.png`), with the format of their suffix
            (`sbmlsim.release.archive.EXTRA_FORMATS`, the antimony of `.ant`
            as the antimony of the flat model).
        sexes: the sexes to build, the sex of the definition among them and
            each once, e.g. `(Sex.MALE, Sex.FEMALE)`; `None` builds the sex of
            the definition only.
        figure: write the figure of the definition, `<sid>.svg` (the key
            `<sid>_svg`, the first extra file of the archive); one figure for
            all sexes, which have the same organs, substances and processes,
            its title names the builds.

    Every file of the archive has one location rule: it is below the parent
    of `model_output_dir`, and its location in the archive is its path
    relative to it. The built models, their markdown and antimony are in
    `model_output_dir` (`./models/<name>` for the directory `models` of a
    model package), the licenses in the parent (`./mit.txt`), an extra model
    or file where it lies (`./models/<sid>_fluid.xml`, `./figures/<x>.png`).
    The README alone is placed, as `./README.md` of the archive, wherever
    the file lies.

    Raises:
        ValueError: if `sexes` leaves out the sex of the definition or names
            a sex twice (before anything is written), or, after the builds and
            before the archive is written, if the stem of an extra model or the
            key of an extra file is the key of another entry, if an extra model
            or file is not below the parent of `model_output_dir` or the suffix
            of an extra file has no format.

    Returns:
        The files with their entries in the COMBINE archive, by the keys of a
        factory: `README`, `cc-by`, `mit`, the sid of the definition, `<sid>_flat`,
        `<sid>_flat_md`, `<sid>_flat_ant` (and the same keys of the id of the
        build of every other sex, `<sid>_female_flat`), the stem of every extra
        model, `<sid>_svg` of the figure and `<stem>_<suffix>` of every extra
        file (`fluid_md`), each `{"path": Path, "entry": ManifestEntry}`.
    """
    sexes = (definition.sex,) if sexes is None else tuple(Sex(s) for s in sexes)
    if definition.sex not in sexes:
        raise ValueError(
            f"The sexes {[str(s) for s in sexes]} of the archive of "
            f"'{definition.sid}' leave out the sex of the definition "
            f"'{definition.sex}', the build of the master entry."
        )
    repeated = sorted({str(s) for s in sexes if sexes.count(s) > 1})
    if repeated:
        raise ValueError(
            f"The archive of '{definition.sid}' builds the sexes {repeated} twice."
        )
    model_output_dir.mkdir(parents=True, exist_ok=True)
    # the build of the sex of the definition first, the master entry
    models = []
    for sex in (definition.sex, *(s for s in sexes if s is not definition.sex)):
        built = build(sex_definition(definition, sex), model_output_dir)
        models.append(ArchiveModel(built.comp_path, flat=built.flat_path))
    figures = []
    if figure:
        path = model_output_dir / f"{definition.sid}.svg"
        figures.append(model_figure(definition, path, sexes))
    return write_archive(
        models,
        model_output_dir.parent,
        name=f"{model}_model",
        copyright=copyright,
        readme=readme,
        extra_models=extra_models,
        extra_files=[*figures, *extra_files],
    )
