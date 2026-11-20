"""Tests of `create_archive`, the factory of a library model in one call."""

import json
import zipfile
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import pytest
from pymetadata.omex import EntryFormat
from sbmlsim.licenses import License, license_text
from sbmlsim.panel import Panel, PanelCurve, PanelSimulation
from sbmlsim.panel_sbml import PANEL_NS, panels_from_sbml

from pkpdsim.pbpk import Sex, build, create_archive, sex_definition

from .conftest import definition

COPYRIGHT = "2024-2026 Matthias König"
MINIMAL_SBML = (
    '<?xml version="1.0" encoding="UTF-8"?>\n'
    '<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">'
    '<model id="fluid"/></sbml>\n'
)


def _create(
    tmp_path: Path,
    *,
    extra: tuple[Path, ...] = (),
    files: tuple[Path, ...] = (),
    readme: bool = True,
) -> dict[str, dict[str, Any]]:
    """Create the archive of the simple definition below `tmp_path`."""
    readme_path = tmp_path / "README.md"
    readme_path.write_text("# simple\n", encoding="utf-8")
    return create_archive(
        definition(),
        tmp_path / "models",
        model="simple",
        copyright=COPYRIGHT,
        readme=readme_path if readme else None,
        extra_models=extra,
        extra_files=files,
    )


def test_the_archive_has_the_flat_model_as_master(tmp_path: Path) -> None:
    """The flat model is the master entry, the files and formats are those of a factory."""
    results = _create(tmp_path)
    assert set(results) == {
        "README",
        "cc-by",
        "mit",
        "simple",
        "simple_flat",
        "simple_flat_md",
        "simple_flat_ant",
        "simple_svg",
        "simple_body_svg",
        "simple_body.readouts_json",
    }
    masters = [k for k, v in results.items() if v["entry"].master]
    assert masters == ["simple_flat"]
    assert results["simple_flat"]["entry"].location == "./models/simple_flat.xml"
    assert results["simple"]["entry"].format == EntryFormat.SBML_L3V2
    assert results["simple_flat_md"]["entry"].format == EntryFormat.MARKDOWN
    assert results["simple_flat_ant"]["entry"].format == EntryFormat.ANTIMONY
    assert results["README"]["entry"].location == "./README.md"
    for info in results.values():
        assert info["path"].exists()
    omex = tmp_path / "simple_model.omex"
    with zipfile.ZipFile(omex) as z:
        names = set(z.namelist())
    assert {
        "README.md",
        "cc-by-sa-4.0.txt",
        "mit.txt",
        "models/simple.xml",
        "models/simple_flat.xml",
        "models/simple_flat.md",
        "models/simple_flat.ant",
        "models/simple.svg",
        "models/simple_body.svg",
        "models/simple_body.readouts.json",
        "manifest.xml",
    } <= names


def test_the_snippets_lie_beside_the_figures(tmp_path: Path) -> None:
    """The interactive snippet of every figure, for the docs, not in the archive."""
    results = _create(tmp_path)
    for stem, kind in (("simple", "schematic"), ("simple_body", "anatomy")):
        snippet = (tmp_path / "models" / f"{stem}.html").read_text(encoding="utf-8")
        assert snippet.startswith(f'<div id="pkf-simple-{kind}"')
    assert not [key for key in results if key.endswith("_html")]
    with zipfile.ZipFile(tmp_path / "simple_model.omex") as z:
        assert not [name for name in z.namelist() if name.endswith(".html")]


def test_the_readouts_of_the_body_figure_are_an_entry(tmp_path: Path) -> None:
    """The sidecar of the anatomical figure is an entry of the archive, as JSON.

    A reader of the archive has the places of the drawing and what each of
    them reads without an import of pkpdsim.
    """
    results = _create(tmp_path)
    info = results["simple_body.readouts_json"]
    assert info["path"] == tmp_path / "models" / "simple_body.readouts.json"
    assert info["entry"].location == "./models/simple_body.readouts.json"
    assert info["entry"].format == "https://purl.org/NET/mediatypes/application/json"
    assert not info["entry"].master
    payload = json.loads(info["path"].read_text(encoding="utf-8"))
    assert payload["contract"] == "figurereadouts/1"
    assert payload["keys"]
    with zipfile.ZipFile(tmp_path / "simple_model.omex") as z:
        assert json.loads(z.read("models/simple_body.readouts.json")) == payload


def test_the_archive_carries_the_extra_models(tmp_path: Path) -> None:
    """An extra SBML file is an entry by its stem, after the models.

    The figure and the extra files come after the extra models, the order of
    `sbmlsim.release.archive.write_archive`.
    """
    (tmp_path / "models").mkdir()
    extra = tmp_path / "models" / "simple_fluid.xml"
    extra.write_text(MINIMAL_SBML, encoding="utf-8")
    results = _create(tmp_path, extra=(extra,))
    assert list(results)[-4:] == [
        "simple_fluid",
        "simple_svg",
        "simple_body_svg",
        "simple_body.readouts_json",
    ]
    entry = results["simple_fluid"]["entry"]
    assert entry.location == "./models/simple_fluid.xml"
    assert entry.format == EntryFormat.SBML_L3V2
    assert not entry.master
    with zipfile.ZipFile(tmp_path / "simple_model.omex") as z:
        assert "models/simple_fluid.xml" in z.namelist()


def test_the_licenses_are_the_texts(tmp_path: Path) -> None:
    """The license files are the texts, the MIT one with the copyright."""
    results = _create(tmp_path)
    assert results["cc-by"]["path"] == tmp_path / "cc-by-sa-4.0.txt"
    assert results["cc-by"]["path"].read_bytes() == license_text(
        License.CC_BY_SA_4_0
    ).encode("utf-8")
    assert results["mit"]["path"].read_bytes() == license_text(
        License.MIT, copyright=COPYRIGHT
    ).encode("utf-8")


def test_the_archive_returns_the_entries_of_a_factory(tmp_path: Path) -> None:
    """Every value is a path and an entry, without a README the key is left out."""
    results = _create(tmp_path, readme=False)
    assert "README" not in results
    for info in results.values():
        assert set(info) == {"path", "entry"}
        assert isinstance(info["path"], Path)


def test_an_extra_model_may_not_take_a_key_of_the_archive(tmp_path: Path) -> None:
    """An extra model named like a built file is an error."""
    (tmp_path / "models").mkdir()
    extra = tmp_path / "models" / "simple_flat.xml"
    extra.write_text(MINIMAL_SBML, encoding="utf-8")
    with pytest.raises(ValueError, match="simple_flat"):
        _create(tmp_path, extra=(extra,))


def test_the_archive_carries_extra_files(tmp_path: Path) -> None:
    """A further file (the ODEs of an extra model, a figure) is an entry by its name.

    The key is the stem with the suffix (`simple_fluid_md`), the location the
    path below the parent of the directory of the models, the format the one
    of the suffix.
    """
    (tmp_path / "models").mkdir()
    (tmp_path / "figures").mkdir()
    extra = tmp_path / "models" / "simple_fluid.xml"
    extra.write_text(MINIMAL_SBML, encoding="utf-8")
    md = tmp_path / "models" / "simple_fluid.md"
    md.write_text("# fluid\n", encoding="utf-8")
    png = tmp_path / "figures" / "simple.png"
    png.write_bytes(b"\x89PNG")
    results = _create(tmp_path, extra=(extra,), files=(md, png))
    assert list(results)[-6:] == [
        "simple_fluid",
        "simple_svg",
        "simple_body_svg",
        "simple_body.readouts_json",
        "simple_fluid_md",
        "simple_png",
    ]
    assert results["simple_fluid_md"]["entry"].location == "./models/simple_fluid.md"
    assert results["simple_fluid_md"]["entry"].format == EntryFormat.MARKDOWN
    assert results["simple_png"]["entry"].location == "./figures/simple.png"
    assert results["simple_png"]["entry"].format == EntryFormat.PNG
    assert not any(
        results[k]["entry"].master for k in ("simple_fluid_md", "simple_png")
    )
    with zipfile.ZipFile(tmp_path / "simple_model.omex") as z:
        names = set(z.namelist())
    assert {"models/simple_fluid.md", "figures/simple.png"} <= names


@pytest.mark.parametrize(
    ("name", "match"),
    [
        ("simple_flat.md", "simple_flat_md"),
        ("simple.unknown", r"suffix '\.unknown'"),
    ],
    ids=["key of the archive", "unknown suffix"],
)
def test_an_extra_file_is_checked(tmp_path: Path, name: str, match: str) -> None:
    """An extra file with the key of an entry or a suffix of no format is an error."""
    (tmp_path / "models").mkdir()
    extra = tmp_path / "models" / name
    extra.write_text("x\n", encoding="utf-8")
    with pytest.raises(ValueError, match=match):
        _create(tmp_path, files=(extra,))


def test_an_extra_file_outside_the_archive_is_an_error(tmp_path: Path) -> None:
    """A file which is not below the parent of the directory of the models has no location."""
    outside = tmp_path.parent / f"{tmp_path.name}_outside.md"
    outside.write_text("x\n", encoding="utf-8")
    try:
        with pytest.raises(ValueError, match="not below"):
            _create(tmp_path, files=(outside,))
    finally:
        outside.unlink()


def test_an_extra_model_has_the_location_rule_of_an_extra_file(tmp_path: Path) -> None:
    """An extra model below the parent is located by its path, outside it is an error."""
    (tmp_path / "standalone").mkdir()
    extra = tmp_path / "standalone" / "simple_liver.xml"
    extra.write_text(MINIMAL_SBML, encoding="utf-8")
    results = _create(tmp_path, extra=(extra,))
    assert results["simple_liver"]["entry"].location == "./standalone/simple_liver.xml"
    with zipfile.ZipFile(tmp_path / "simple_model.omex") as z:
        assert "standalone/simple_liver.xml" in z.namelist()
    outside = tmp_path.parent / f"{tmp_path.name}_outside.xml"
    outside.write_text(MINIMAL_SBML, encoding="utf-8")
    try:
        with pytest.raises(ValueError, match=r"extra model .* not below"):
            _create(tmp_path, extra=(outside,))
    finally:
        outside.unlink()


def test_an_extra_antimony_file_is_antimony(tmp_path: Path) -> None:
    """The antimony of an extra model has the format of the antimony of the flat model."""
    (tmp_path / "models").mkdir()
    ant = tmp_path / "models" / "simple_fluid.ant"
    ant.write_text("model fluid()\nend\n", encoding="utf-8")
    results = _create(tmp_path, files=(ant,))
    assert results["simple_fluid_ant"]["entry"].format == EntryFormat.ANTIMONY


def test_the_built_files_have_the_location_rule_of_the_extra_files(
    tmp_path: Path,
) -> None:
    """Every file of the archive is located by its path below the parent of the models.

    The built models, their markdown, antimony and figure and an extra file
    in the directory of the models share one rule, also when the directory is not
    named `models`; the README is the `README.md` of the archive.
    """
    readme = tmp_path / "README.md"
    readme.write_text("# simple\n", encoding="utf-8")
    out = tmp_path / "sbml"
    out.mkdir()
    extra = out / "notes.md"
    extra.write_text("notes\n", encoding="utf-8")
    results = create_archive(
        definition(),
        out,
        model="simple",
        copyright=COPYRIGHT,
        readme=readme,
        extra_files=(extra,),
    )
    locations = {k: v["entry"].location for k, v in results.items()}
    assert locations == {
        "README": "./README.md",
        "cc-by": "./cc-by-sa-4.0.txt",
        "mit": "./mit.txt",
        "simple": "./sbml/simple.xml",
        "simple_flat": "./sbml/simple_flat.xml",
        "simple_flat_md": "./sbml/simple_flat.md",
        "simple_flat_ant": "./sbml/simple_flat.ant",
        "simple_svg": "./sbml/simple.svg",
        "simple_body_svg": "./sbml/simple_body.svg",
        "simple_body.readouts_json": "./sbml/simple_body.readouts.json",
        "notes_md": "./sbml/notes.md",
    }
    with zipfile.ZipFile(tmp_path / "simple_model.omex") as z:
        assert "sbml/simple_flat.xml" in z.namelist()


def _panel(pid: str) -> Panel:
    """Get a panel of the simple model."""
    return Panel(
        id=pid,
        title=f"panel {pid}",
        curves=(PanelCurve(id="[Cve_drg]", label="drug in venous blood"),),
        simulation=PanelSimulation(end=10.0),
    )


A, B, C = _panel("a"), _panel("b"), _panel("c")

#: the panels of the archives with builds per sex, by the name of the case
PANEL_CASES: dict[str, Sequence[Panel] | Mapping[Sex | str, Sequence[Panel]] | None] = {
    "none": None,
    "one sex": {Sex.MALE: [A, B]},
    "both sexes": {"male": [A], Sex.FEMALE: [C, A]},
    "every build": [A],
}

#: the panels of the flat model of every build, by case and key of the build
EXPECTED_PANELS: dict[str, dict[str, list[Panel]]] = {
    "none": {"simple_flat": [], "simple_female_flat": []},
    "one sex": {"simple_flat": [A, B], "simple_female_flat": []},
    "both sexes": {"simple_flat": [A], "simple_female_flat": [C, A]},
    "every build": {"simple_flat": [A], "simple_female_flat": [A]},
}


@pytest.fixture(scope="module")
def panel_archives(
    tmp_path_factory: pytest.TempPathFactory,
) -> dict[str, tuple[Path, dict[str, dict[str, Any]]]]:
    """The archives of the male and female build of every case of `PANEL_CASES`.

    Each with its directory, the parent of the models, and the files of
    `create_archive`.
    """
    archives = {}
    for case, panels in PANEL_CASES.items():
        root = tmp_path_factory.mktemp("panels")
        files = create_archive(
            definition(),
            root / "models",
            model="simple",
            copyright=COPYRIGHT,
            sexes=(Sex.MALE, Sex.FEMALE),
            figure=False,
            panels=panels,
        )
        archives[case] = (root, files)
    return archives


def _members(root: Path) -> dict[str, bytes]:
    """Get the members of the archive below `root` by their name."""
    with zipfile.ZipFile(root / "simple_model.omex") as z:
        return {name: z.read(name) for name in z.namelist()}


@pytest.mark.parametrize("case", list(PANEL_CASES))
def test_the_flat_models_of_the_archive_carry_the_panels_of_their_build(
    panel_archives: dict[str, tuple[Path, dict[str, dict[str, Any]]]], case: str
) -> None:
    """The flat model of a build in the archive has exactly the panels of its build.

    A sex without panels has no annotation of panels, and neither has a comp
    model.
    """
    root, files = panel_archives[case]
    members = _members(root)
    for key, expected in EXPECTED_PANELS[case].items():
        entry = files[key]["entry"].location.removeprefix("./")
        panels, problems = panels_from_sbml(members[entry].decode("utf-8"))
        assert problems == ()
        assert list(panels.values()) == expected
    for key in ("simple", "simple_female"):
        assert PANEL_NS.encode() not in files[key]["path"].read_bytes()


@pytest.mark.parametrize("case", list(PANEL_CASES))
def test_the_archive_is_written_once_from_the_annotated_files(
    panel_archives: dict[str, tuple[Path, dict[str, dict[str, Any]]]], case: str
) -> None:
    """Every entry of the archive is the file on disk, the annotated flat models among them."""
    root, files = panel_archives[case]
    members = _members(root)
    for info in files.values():
        entry = info["entry"].location.removeprefix("./")
        assert members[entry] == info["path"].read_bytes(), entry


@pytest.mark.parametrize("case", [c for c in PANEL_CASES if c != "none"])
def test_the_panels_change_only_the_annotated_flat_models(
    panel_archives: dict[str, tuple[Path, dict[str, dict[str, Any]]]], case: str
) -> None:
    """An archive with panels is the archive without them but for its annotated flat models.

    The comp models, the ODEs and the antimony of every build, the licenses
    and the flat model of a build without panels are the bytes of the archive
    without panels, so `panels=None` (the default) is the archive as before
    the keyword; the manifest and the metadata, which carry the dates of the
    archive, are left out of the comparison.
    """
    none = _members(panel_archives["none"][0])
    root, files = panel_archives[case]
    members = _members(root)
    assert set(members) == set(none)
    annotated = {
        files[key]["entry"].location.removeprefix("./")
        for key, expected in EXPECTED_PANELS[case].items()
        if expected
    }
    dated = {"manifest.xml", "metadata.rdf"}
    for name in set(members) - annotated - dated:
        assert members[name] == none[name], name
    for name in annotated:
        assert members[name] != none[name]


def test_without_panels_the_built_models_are_the_build(
    panel_archives: dict[str, tuple[Path, dict[str, dict[str, Any]]]],
    tmp_path: Path,
) -> None:
    """`panels=None` writes no annotation: the models are the bytes of `build`."""
    members = _members(panel_archives["none"][0])
    for sex in (Sex.MALE, Sex.FEMALE):
        built = build(sex_definition(definition(), sex), tmp_path)
        for path in (built.comp_path, built.flat_path):
            assert members[f"models/{path.name}"] == path.read_bytes()
            assert PANEL_NS.encode() not in path.read_bytes()


@pytest.mark.parametrize(
    ("panels", "sexes", "error", "match"),
    [
        ({"cp": [_panel("a")]}, None, ValueError, "'cp', which is no sex"),
        ({Sex.FEMALE: [_panel("a")]}, None, ValueError, "'female', which is not built"),
        ([_panel("a"), _panel("a")], None, ValueError, r"repeat the ids \['a'\]"),
        ([object()], None, TypeError, r"no panels .*\['object'\]"),
        (
            {Sex.FEMALE: [_panel("b"), _panel("b")]},
            (Sex.MALE, Sex.FEMALE),
            ValueError,
            "female build of 'simple' repeat",
        ),
    ],
    ids=["no sex", "not built", "repeated id", "no panel", "repeated in a build"],
)
def test_the_panels_are_checked_before_anything_is_written(
    tmp_path: Path,
    panels: Any,  # noqa: ANN401
    sexes: tuple[Sex, ...] | None,
    error: type[Exception],
    match: str,
) -> None:
    """A key of no built sex, a repeated id or no panel raises before a file is written."""
    with pytest.raises(error, match=match):
        create_archive(
            definition(),
            tmp_path / "models",
            model="simple",
            copyright=COPYRIGHT,
            sexes=sexes,
            panels=panels,
        )
    assert list(tmp_path.iterdir()) == []


def test_an_empty_sequence_of_panels_annotates_nothing(tmp_path: Path) -> None:
    """A build whose panels are empty is written as without panels."""
    files = create_archive(
        definition(),
        tmp_path / "models",
        model="simple",
        copyright=COPYRIGHT,
        figure=False,
        panels={Sex.MALE: []},
    )
    assert PANEL_NS.encode() not in files["simple_flat"]["path"].read_bytes()
