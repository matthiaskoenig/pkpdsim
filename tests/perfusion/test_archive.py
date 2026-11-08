"""Tests of the archive of a machine perfusion: the model, its files and the COMBINE archive."""

from pathlib import Path

from pymetadata.omex import EntryFormat, Omex

from pkpdsim.perfusion.archive import create_archive

from .conftest import kidney_definition


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
    }
    for key, (location, path) in locations.items():
        assert results[key]["entry"].location == location, key
        assert results[key]["path"] == path, key
        assert path.is_file(), key
    assert results["kidney_circuit"]["entry"].master
    assert results["kidney_circuit_svg"]["entry"].format == EntryFormat.SVG
    assert (tmp_path / "mit.txt").is_file() and (
        tmp_path / "cc-by-sa-4.0.txt"
    ).is_file()
    omex = Omex.from_omex(tmp_path / "kidney_circuit_model.omex")
    assert "./models/kidney_circuit.svg" in {e.location for e in omex.manifest.entries}
