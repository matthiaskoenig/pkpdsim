"""Tests of `create_archive`, the factory of a library model in one call."""

import zipfile
from pathlib import Path

from pymetadata.omex import EntryFormat

from pkpdsim.licenses import License, license_text
from pkpdsim.pbpk import create_archive

from .conftest import definition

COPYRIGHT = "2024-2026 Matthias König"


def _create(
    tmp_path: Path, *, extra: tuple[Path, ...] = (), readme: bool = True
) -> dict:
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
        "manifest.xml",
    } <= names


def test_the_archive_carries_the_extra_models(tmp_path: Path) -> None:
    """An extra SBML file is an entry by its stem, after the others."""
    (tmp_path / "models").mkdir()
    extra = tmp_path / "models" / "simple_fluid.xml"
    extra.write_text("<sbml/>", encoding="utf-8")
    results = _create(tmp_path, extra=(extra,))
    assert list(results)[-1] == "simple_fluid"
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
