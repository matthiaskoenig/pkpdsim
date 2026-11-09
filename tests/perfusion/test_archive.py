"""Tests of the archive of a machine perfusion: the model, its files and the COMBINE archive."""

from pathlib import Path

from pymetadata.omex import EntryFormat, Omex

from pkpdsim.perfusion.archive import create_archive

from .conftest import kidney_definition, liver_definition


def test_the_archive_has_the_model_its_files_and_the_figure(tmp_path: Path) -> None:
    """The model, its ODEs and antimony, the figure, the licenses and the README."""
    readme = tmp_path / "README.md"
    readme.write_text("# perfusion\n", encoding="utf-8")
    d = kidney_definition(sid="kidney_circuit")
    results = create_archive(
        d,
        tmp_path / "models",
        model="kidney_circuit",
        copyright="2026 Matthias König",
        readme=readme,
    )
    assert set(results) == {
        "README",
        "cc-by",
        "mit",
        "kidney_circuit",
        "kidney_circuit_md",
        "kidney_circuit_ant",
        "kidney_circuit_svg",
        "kidney_circuit_anatomy_svg",
    }
    models = tmp_path / "models"
    locations = {
        "kidney_circuit": (
            "./models/kidney_circuit.xml",
            models / "kidney_circuit.xml",
        ),
        "kidney_circuit_md": (
            "./models/kidney_circuit.md",
            models / "kidney_circuit.md",
        ),
        "kidney_circuit_ant": (
            "./models/kidney_circuit.ant",
            models / "kidney_circuit.ant",
        ),
        "kidney_circuit_svg": (
            "./models/kidney_circuit.svg",
            models / "kidney_circuit.svg",
        ),
        "kidney_circuit_anatomy_svg": (
            "./models/kidney_circuit_anatomy.svg",
            models / "kidney_circuit_anatomy.svg",
        ),
    }
    for key, (location, path) in locations.items():
        assert results[key]["entry"].location == location, key
        assert results[key]["path"] == path, key
        assert path.is_file(), key
    assert results["kidney_circuit"]["entry"].master
    for key in ("kidney_circuit_svg", "kidney_circuit_anatomy_svg"):
        assert results[key]["entry"].format == EntryFormat.SVG
    # the schematic first, then the anatomical figure
    keys = list(results)
    assert keys.index("kidney_circuit_svg") < keys.index("kidney_circuit_anatomy_svg")
    assert (tmp_path / "mit.txt").is_file() and (
        tmp_path / "cc-by-sa-4.0.txt"
    ).is_file()
    omex = Omex.from_omex(tmp_path / "kidney_circuit_model.omex")
    assert "./models/kidney_circuit.svg" in {e.location for e in omex.manifest.entries}


def test_the_archive_of_a_liver_has_the_flat_model_as_its_master(
    tmp_path: Path,
) -> None:
    """The liver is simulated as its flat model, which is the master entry; the
    composed model and the cell model are written beside it."""
    results = create_archive(
        liver_definition(2, sid="liver_circuit"),
        tmp_path / "models",
        model="liver_circuit",
        copyright="2026 Matthias König",
    )
    assert results["liver_circuit_flat"]["entry"].master
    assert results["liver_circuit_flat"]["path"] == (
        tmp_path / "models" / "liver_circuit_flat.xml"
    )
    for name in ("liver_circuit.xml", "liver_circuit_core.xml"):
        assert (tmp_path / "models" / name).is_file()
