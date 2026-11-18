"""Cache validated pairs of composed and flat SBML documents."""

from __future__ import annotations

import json
import os
import re
from dataclasses import asdict
from hashlib import sha256
from importlib.metadata import version
from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import BadZipFile, ZipFile
from zlib import error as ZlibError

import libsbml
from sbmlutils.validation import ValidationOptions

from pkpdsim import __version__

# Increment when the flattening or cache protocol changes.
_SCHEMA = 1
_HISTORY_DATE = re.compile(
    rb"(<dcterms:(?:created|modified)\b[^>]*>\s*<dcterms:W3CDTF>)"
    rb"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:Z|[+-]\d{2}:\d{2})?"
    rb"(</dcterms:W3CDTF>)"
)


def has_external_models(comp: bytes) -> bool:
    """Conservatively detect a model whose key would need dependency files."""
    return b"externalModelDefinition" in comp


def _digest(content: bytes) -> str:
    return sha256(content).hexdigest()


def _normalized(comp: bytes) -> bytes:
    return _HISTORY_DATE.sub(rb"\1\2", comp)


def cache_key(
    comp: bytes, comp_options: ValidationOptions, flat_options: ValidationOptions
) -> str:
    """Fingerprint the source and the software and options which validate it."""
    settings = json.dumps(
        {
            "schema": _SCHEMA,
            "pkpdsim": __version__,
            "sbmlutils": version("sbmlutils"),
            "libsbml": libsbml.getLibSBMLDottedVersion(),
            "comp": asdict(comp_options),
            "flat": asdict(flat_options),
        },
        sort_keys=True,
    )
    return _digest(settings.encode() + b"\0" + _normalized(comp))


def cache_path(output_dir: Path, sid: str, key: str) -> Path | None:
    """Use a per-output cache unless a shared directory or 'off' is set."""
    configured = os.environ.get("PKPDSIM_BUILD_CACHE")
    if configured == "off":
        return None
    if configured:
        return Path(configured).expanduser() / f"{key}.zip"
    return output_dir / ".pkpdsim-cache" / f"{sid}.zip"


def read_cache(path: Path, key: str, comp: bytes) -> tuple[bytes, bytes] | None:
    """Read a matching, intact cached pair, or return no hit."""
    try:
        with ZipFile(path) as archive:
            manifest = json.loads(archive.read("manifest.json"))
            cached_comp = archive.read("comp.xml")
            flat = archive.read("flat.xml")
    except OSError, BadZipFile, KeyError, ValueError, RuntimeError, ZlibError:
        return None
    if manifest != _manifest(key, cached_comp, flat):
        return None
    if _normalized(cached_comp) != _normalized(comp):
        return None
    return cached_comp, flat


def write_cache(path: Path, key: str, comp: bytes, flat: bytes) -> None:
    """Publish a validated pair atomically, with an integrity manifest."""
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        with TemporaryDirectory(dir=path.parent) as temporary:
            staged = Path(temporary) / path.name
            with ZipFile(staged, mode="w") as archive:
                archive.writestr(
                    "manifest.json", json.dumps(_manifest(key, comp, flat))
                )
                archive.writestr("comp.xml", comp)
                archive.writestr("flat.xml", flat)
            staged.replace(path)
    except OSError:
        # A build remains usable when its optional cache cannot be written.
        return


def _manifest(key: str, comp: bytes, flat: bytes) -> dict[str, str]:
    return {"key": key, "comp": _digest(comp), "flat": _digest(flat)}


def _read_output(path: Path) -> bytes | None:
    """Read an existing output if present."""
    try:
        return path.read_bytes()
    except FileNotFoundError:
        return None


def publish_outputs(output_dir: Path, staged_files: dict[Path, bytes]) -> None:
    """Keep the existing pair if unchanged, otherwise atomically publish each file."""
    outputs = [
        (staged, output_dir / staged.name, content)
        for staged, content in staged_files.items()
    ]
    previous = [_read_output(path) for _, path, _ in outputs]
    if all(
        existing is not None
        and (existing == content or _normalized(existing) == _normalized(content))
        for (_, _, content), existing in zip(outputs, previous, strict=True)
    ):
        return
    for (staged, path, content), existing in zip(outputs, previous, strict=True):
        if existing != content:
            staged.write_bytes(content)
            staged.replace(path)
