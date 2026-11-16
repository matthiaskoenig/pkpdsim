"""Tests of the archive of a machine perfusion: the model, its files and the COMBINE archive."""

import zipfile
from pathlib import Path

import pytest
from pymetadata.omex import EntryFormat, Omex
from sbmlsim.panel import Panel, PanelCurve, PanelSimulation
from sbmlsim.panel_sbml import PANEL_NS, panels_from_sbml

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
    # the interactive snippet of every figure lies beside it, no entry of the archive
    for stem, kind in (
        ("kidney_circuit", "schematic"),
        ("kidney_circuit_anatomy", "anatomy"),
    ):
        snippet = (models / f"{stem}.html").read_text(encoding="utf-8")
        assert snippet.startswith(f'<div id="pkf-kidney_circuit-{kind}"')
    assert not [key for key in results if key.endswith("_html")]
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
    # the three models the build writes are entries of the archive, the cell
    # model among them: it stays loadable alone, for a mesh or another coupling
    for key in ("liver_circuit", "liver_circuit_flat", "liver_circuit_core"):
        assert results[key]["path"].is_file(), key
    omex = Omex.from_omex(tmp_path / "liver_circuit_model.omex")
    locations = {e.location for e in omex.manifest.entries}
    assert {
        "./models/liver_circuit.xml",
        "./models/liver_circuit_flat.xml",
        "./models/liver_circuit_core.xml",
    } <= locations
    masters = {e.location for e in omex.manifest.entries if e.master}
    assert masters == {"./models/liver_circuit_flat.xml"}


def _panel(pid: str) -> Panel:
    """Get a panel of a perfusion model."""
    return Panel(
        id=pid,
        title=f"panel {pid}",
        curves=(PanelCurve(id="time", label="time"),),
        simulation=PanelSimulation(end=10.0),
    )


def test_the_flat_model_of_a_liver_carries_the_panels(tmp_path: Path) -> None:
    """The panels are the annotation of the flat model, the entry of the archive.

    The comp model and the cell model of the zones carry none, and every
    entry of the archive is the file on disk, so the archive is written once
    from the annotated model.
    """
    panels = [_panel("a"), _panel("b")]
    results = create_archive(
        liver_definition(2, sid="liver_circuit"),
        tmp_path / "models",
        model="liver_circuit",
        copyright="2026 Matthias König",
        panels=panels,
    )
    flat = results["liver_circuit_flat"]["path"]
    read, problems = panels_from_sbml(flat)
    assert problems == ()
    assert list(read.values()) == panels
    for key in ("liver_circuit", "liver_circuit_core"):
        assert PANEL_NS.encode() not in results[key]["path"].read_bytes()
    with zipfile.ZipFile(tmp_path / "liver_circuit_model.omex") as z:
        for info in results.values():
            entry = info["entry"].location.removeprefix("./")
            assert z.read(entry) == info["path"].read_bytes(), entry


def test_without_panels_the_model_has_no_annotation_of_panels(tmp_path: Path) -> None:
    """`panels=None`, the default, writes no annotation, and neither do no panels."""
    for panels in (None, []):
        out = tmp_path / str(panels)
        results = create_archive(
            kidney_definition(sid="kidney_circuit"),
            out / "models",
            model="kidney_circuit",
            copyright="2026 Matthias König",
            panels=panels,
        )
        model = results["kidney_circuit"]["path"].read_bytes()
        assert PANEL_NS.encode() not in model


def test_repeated_panels_are_refused_before_anything_is_written(
    tmp_path: Path,
) -> None:
    """Two panels of one id raise before the model is built."""
    with pytest.raises(ValueError, match=r"the model 'kidney_circuit' repeat"):
        create_archive(
            kidney_definition(sid="kidney_circuit"),
            tmp_path / "models",
            model="kidney_circuit",
            copyright="2026 Matthias König",
            panels=[_panel("a"), _panel("a")],
        )
    assert list(tmp_path.iterdir()) == []
