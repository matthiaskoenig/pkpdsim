"""Reuse only intact, validated PBPK builds."""

import os
from dataclasses import replace
from pathlib import Path
from zipfile import ZipFile

import libsbml
import pytest
from sbmlutils.factory import Creator, ExternalModelDefinition, Model, Submodel
from sbmlutils.factory import model as factory_model

from pkpdsim.pbpk import builder
from pkpdsim.pbpk.builder import BuildResult
from pkpdsim.pbpk.definition import PBPKDefinition
from pkpdsim.sbml import cache

from .conftest import definition


@pytest.fixture(autouse=True)
def local_cache(monkeypatch: pytest.MonkeyPatch) -> None:
    """Use the production default, regardless of the test session cache."""
    monkeypatch.delenv("PKPDSIM_BUILD_CACHE", raising=False)


def _outputs(result: BuildResult) -> tuple[Path, Path]:
    return result.comp_path, result.flat_path


def test_unchanged_generated_history_preserves_both_outputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    d = definition(
        creators=[
            Creator(
                familyName="Tester",
                givenName="Model",
                email="model@example.org",
                organization="Model Lab",
            )
        ]
    )
    monkeypatch.setattr(
        factory_model, "date_now", lambda: libsbml.Date("2025-01-01T00:00:00Z")
    )
    result = builder.build(d, tmp_path)
    contents = [path.read_bytes() for path in _outputs(result)]
    assert b"2025-01-01" in contents[0]
    for path in _outputs(result):
        os.utime(path, ns=(1_000_000_000, 1_000_000_000))

    monkeypatch.setattr(
        factory_model, "date_now", lambda: libsbml.Date("2025-01-02T00:00:00Z")
    )
    with monkeypatch.context() as patch:
        patch.setattr(builder, "_validate", lambda *args: pytest.fail("revalidated"))
        patch.setattr(builder, "_flatten", lambda *args: pytest.fail("reflattened"))
        assert builder.build(d, tmp_path) == result
        assert [path.read_bytes() for path in _outputs(result)] == contents
        assert [path.stat().st_mtime_ns for path in _outputs(result)] == [
            1_000_000_000
        ] * 2

    monkeypatch.setenv("PKPDSIM_BUILD_CACHE", "off")
    assert builder.build(d, tmp_path) == result
    assert [path.read_bytes() for path in _outputs(result)] == contents
    assert [path.stat().st_mtime_ns for path in _outputs(result)] == [1_000_000_000] * 2


def test_shared_cache_restores_the_validated_pair(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    shared = tmp_path / "shared"
    monkeypatch.setenv("PKPDSIM_BUILD_CACHE", str(shared))
    first = builder.build(definition(), tmp_path / "first")
    monkeypatch.setattr(builder, "flat_model", lambda *args: pytest.fail("cache miss"))
    second = builder.build(definition(), tmp_path / "second")
    assert len(list(shared.glob("*.zip"))) == 1
    assert [path.read_bytes() for path in _outputs(first)] == [
        path.read_bytes() for path in _outputs(second)
    ]
    second.flat_path.write_bytes(b"broken output")
    builder.build(definition(), tmp_path / "second")
    assert second.flat_path.read_bytes() == first.flat_path.read_bytes()


@pytest.fixture
def validations(monkeypatch: pytest.MonkeyPatch) -> list[bool]:
    recorded: list[bool] = []
    validate = builder._validate

    def record(
        doc: libsbml.SBMLDocument, path: Path, options: builder.ValidationOptions
    ) -> None:
        recorded.append(options.units_consistency)
        validate(doc, path, options)

    monkeypatch.setattr(builder, "_validate", record)
    return recorded


@pytest.mark.parametrize("change", ["definition", "version", "settings", "off"])
def test_build_revalidates_after_a_cache_key_change(
    change: str,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    validations: list[bool],
) -> None:
    d = definition()
    result = builder.build(d, tmp_path)
    validations.clear()

    match change:
        case "definition":
            d = replace(d, name="Renamed model")
        case "version":
            monkeypatch.setattr(cache, "__version__", "future")
        case "settings":
            fast = builder.ValidationOptions.fast
            monkeypatch.setattr(
                builder.ValidationOptions,
                "fast",
                lambda: replace(fast(), log_errors=False),
            )
        case "off":
            monkeypatch.setenv("PKPDSIM_BUILD_CACHE", "off")

    assert builder.build(d, tmp_path) == result
    assert validations == [False, True]


def test_corrupt_cache_member_is_rebuilt(
    tmp_path: Path, validations: list[bool]
) -> None:
    d = definition()
    result = builder.build(d, tmp_path)
    entry = tmp_path / ".pkpdsim-cache" / f"{d.sid}.zip"
    with ZipFile(entry) as archive:
        manifest = archive.read("manifest.json")
        comp = archive.read("comp.xml")
    with ZipFile(entry, "w") as archive:
        archive.writestr("manifest.json", manifest)
        archive.writestr("comp.xml", comp)
        archive.writestr("flat.xml", b"corrupt")
    original_flat = result.flat_path.read_bytes()
    validations.clear()

    assert builder.build(d, tmp_path) == result
    assert validations == [False, True]
    assert result.flat_path.read_bytes() == original_flat


def test_failed_validation_is_never_published(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    validate = builder._validate

    def fail_strict(
        doc: libsbml.SBMLDocument, path: Path, options: builder.ValidationOptions
    ) -> None:
        if options.units_consistency:
            raise ValueError("flat validation failed")
        validate(doc, path, options)

    d = definition()
    monkeypatch.setattr(builder, "_validate", fail_strict)
    with pytest.raises(ValueError, match="flat validation failed"):
        builder.build(d, tmp_path)
    assert not list(tmp_path.rglob("*.zip"))
    assert not list(tmp_path.glob("*.xml"))


def test_external_sources_are_resolved_relative_to_the_output_directory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, validations: list[bool]
) -> None:
    source = tmp_path / "component.xml"
    doc = libsbml.SBMLDocument(3, 2)
    doc.createModel().setId("component")
    libsbml.writeSBMLToFile(doc, str(source))
    compose = builder.compose

    def external_model(d: PBPKDefinition) -> Model:
        model = compose(d)
        model.external_model_definitions.append(
            ExternalModelDefinition(
                "component", source=source.name, modelRef="component"
            )
        )
        model.submodels.append(Submodel("external", modelRef="component"))
        return model

    monkeypatch.setattr(builder, "compose", external_model)
    d = definition()
    builder.build(d, tmp_path)
    builder.build(d, tmp_path)
    assert validations == [False, True, False, True]
    assert not list(tmp_path.rglob("*.zip"))
