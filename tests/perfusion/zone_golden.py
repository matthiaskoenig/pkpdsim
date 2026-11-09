"""The objects of a zone of the cellular core as they were before the cell model.

`tests/perfusion/data/zone_objects.json` was written from `core.zone_objects`
at commit 37679df, before the emitter was shared with `cell_model`, by
`python -m perfusion.zone_golden` run from `tests/`. It is the one witness that
the move changed no id, value, unit, formula or SBO term, so it is never
regenerated from newer code: a regeneration would make the tests which read it
agree with whatever the code says.

Every public attribute of an object is described (the formula of a rule or a
kinetic law as text, units by their id, SBO terms by their value), the key of
an object is its id (the variable or symbol of a rule or an assignment).
"""

import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from pkpdsim.perfusion import core

GOLDEN = Path(__file__).parent / "data" / "zone_objects.json"
CELLS = "Vcell_z1"
VASCULAR = "Vvas_z1"
ZONE = "z1"

CONFIGS: dict[str, dict[str, Any]] = {
    "plain": {},
    "profile": {"profile": {"OXPHOS": 1.5, "BILE": 0.5}},
    "reduced": {"released": ("alt", "k"), "bile": False},
}
"""The configurations of the golden file, the keywords of `zone_objects`."""

KINDS = ("species", "reactions", "parameters", "rules", "rate_rules", "assignments")


def _plain(value: Any) -> Any:  # noqa: ANN401
    """Get a JSON value of an attribute of an object."""
    if value is None or isinstance(value, bool | int | float | str):
        return value
    if isinstance(value, list | tuple):
        return [_plain(v) for v in value]
    if hasattr(value, "to_string"):
        return str(value.to_string())
    if hasattr(value, "math"):
        return {
            "math": str(value.math),
            "unit": _plain(value.unit),
            "local_parameters": _plain(value.local_parameters),
        }
    sid = getattr(value, "sid", None)
    if sid is not None:
        return sid
    inner = getattr(value, "value", None)
    if isinstance(inner, str):
        return inner
    return str(value)


def describe(obj: object) -> dict[str, object]:
    """Describe an object by all its public attributes."""
    return {
        name: _plain(attribute)
        for name, attribute in sorted(vars(obj).items())
        if not name.startswith("_")
    }


def key(obj: object) -> str:
    """Get the id of an object, the variable or the symbol of a rule."""
    for name in ("sid", "variable", "symbol"):
        found = getattr(obj, name, None)
        if found:
            return str(found)
    raise ValueError(f"No id of {obj!r}.")


def describe_objects(objects: core.CoreObjects) -> dict[str, dict[str, Any]]:
    """Describe the objects by kind and id."""
    return {
        kind: {key(o): describe(o) for o in getattr(objects, kind)} for kind in KINDS
    }


def load() -> Mapping[str, dict[str, dict[str, Any]]]:
    """Get the golden file, by configuration."""
    return json.loads(GOLDEN.read_text())


def main() -> None:
    """Write the golden file from `core.zone_objects`."""
    params = core.CoreParameters.default()
    out = {
        name: describe_objects(
            core.zone_objects(ZONE, CELLS, VASCULAR, params, **options)
        )
        for name, options in CONFIGS.items()
    }
    GOLDEN.write_text(json.dumps(out, indent=1, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
