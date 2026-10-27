"""The factory of a library model in one call: the files and the COMBINE archive.

A model package builds its `PBPKDefinition` with `pkpdsim.pbpk.build` and
writes next to the SBML the differential equations as markdown, the antimony
of the flat model, the two license texts and the COMBINE archive.
`create_archive` does all of it, so the factory of a model is one call.
"""

from collections.abc import Sequence
from pathlib import Path
from typing import Any

from pkpdsim.pbpk.builder import build
from pkpdsim.pbpk.definition import PBPKDefinition


def create_archive(
    definition: PBPKDefinition,
    model_output_dir: Path,
    *,
    model: str,
    copyright: str,
    readme: Path | None = None,
    extra_models: Sequence[Path] = (),
) -> dict[str, dict[str, Any]]:
    """Build a definition and write its files and the COMBINE archive.

    Builds the definition into `model_output_dir` (`build`), writes the ODEs of
    the flat model as markdown (`sbmlode`) and its antimony, the texts of CC
    BY-SA 4.0 and MIT (with the copyright) into the parent of the directory and
    the archive `<model>_model.omex` there, with the README, the licenses, the
    comp model, the flat model (the master entry), the markdown, the antimony
    and the extra models.

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

    Raises:
        ValueError: if the stem of an extra model is the key of another entry.

    Returns:
        The files with their entries in the COMBINE archive, by the keys of a
        factory: `README`, `cc-by`, `mit`, the sid of the definition, `<sid>_flat`,
        `<sid>_flat_md`, `<sid>_flat_ant` and the stem of every extra model, each
        `{"path": Path, "entry": ManifestEntry}`.
    """
    from pymetadata.omex import EntryFormat, ManifestEntry, Omex
    from sbmlode import OdeSystem
    from sbmlutils.io import sbml_to_antimony

    from pkpdsim.licenses import License, write_license

    def entry(
        location: str, fmt: EntryFormat, *, master: bool = False
    ) -> ManifestEntry:
        return ManifestEntry(location=location, format=fmt, master=master)

    model_output_dir.mkdir(parents=True, exist_ok=True)
    built = build(definition, model_output_dir)
    flat = built.flat_path

    # differential equations and antimony of the flat model
    md_path = flat.with_suffix(".md")
    OdeSystem.from_sbml(flat).write(md_path)
    ant_path = flat.with_suffix(".ant")
    ant_path.write_text(sbml_to_antimony(flat), encoding="utf-8")

    parent = model_output_dir.parent
    results: dict[str, dict[str, Any]] = {}
    if readme is not None:
        results["README"] = {
            "path": readme,
            "entry": entry("./README.md", EntryFormat.MARKDOWN),
        }
    cc_by = write_license(parent, License.CC_BY_SA_4_0)
    mit = write_license(parent, License.MIT, copyright=copyright)
    results["cc-by"] = {
        "path": cc_by,
        "entry": entry(f"./{cc_by.name}", EntryFormat.TXT),
    }
    results["mit"] = {
        "path": mit,
        "entry": entry(f"./{mit.name}", EntryFormat.TXT),
    }
    results[definition.sid] = {
        "path": built.comp_path,
        "entry": entry(f"./models/{built.comp_path.name}", EntryFormat.SBML_L3V2),
    }
    results[f"{definition.sid}_flat"] = {
        "path": flat,
        "entry": entry(f"./models/{flat.name}", EntryFormat.SBML_L3V2, master=True),
    }
    results[f"{definition.sid}_flat_md"] = {
        "path": md_path,
        "entry": entry(f"./models/{md_path.name}", EntryFormat.MARKDOWN),
    }
    results[f"{definition.sid}_flat_ant"] = {
        "path": ant_path,
        "entry": entry(f"./models/{ant_path.name}", EntryFormat.ANTIMONY),
    }
    for extra in extra_models:
        if extra.stem in results:
            raise ValueError(
                f"The extra model {extra.name!r} has the key {extra.stem!r} of an "
                "entry of the archive."
            )
        results[extra.stem] = {
            "path": extra,
            "entry": entry(f"./models/{extra.name}", EntryFormat.SBML_L3V2),
        }

    omex = Omex()
    for info in results.values():
        omex.add_entry(entry_path=info["path"], entry=info["entry"])
    omex.to_omex(omex_path=parent / f"{model}_model.omex")
    return results
