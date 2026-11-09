"""The map of the ids of the perfusion liver from before the zonation engine to after.

`tests/perfusion/data/organ_id_map.json` is the record of what moved when the
zones of the liver went onto the engine (the vascular space, the chain, the
cells as a submodel of every zone): every id of `data/organ_ids.json` (the
model before, commit 98aacca) which has another id in the model now, the ids
the model lost and the ids it gained. It is read by `test_organ_reference`,
which checks it against the model that is built and uses it to compare the
trajectories of the two models id by id.

The rules below are written from the id schemes (the ids of the engine,
`pkpdsim.zonation.ids`, and the flat ids of a comp submodel, `M<kk>__<id>`),
not read off the new model: an old id which a rule sends to an id the model
does not have fails the test. Run from `tests/` with
`python -m perfusion.organ_id_map` after a deliberate change of the ids.
"""

import json
import re
from pathlib import Path
from tempfile import TemporaryDirectory

import libsbml

from pkpdsim.perfusion.builder import build

from .organ_reference import DATA, IDS, definition

MAP = DATA / "organ_id_map.json"
PERFUSATE = {"o2", "glc", "lac", "h", "co2", "hco3", "k", "alt", "ast", "ldh", "tr"}
"""The quantities of the perfusate of the reference (the tracer `tr` included)."""
CELL_PROCESSES = re.compile(r"[A-Z0-9]+")


def rename(old: str) -> str | None:
    """Get the id in the model now of an id of the model before, `None` if it is gone."""
    if old == "V_li":
        return None
    if m := re.fullmatch(r"(.+)_z(\d+)_(cell|dead)", old):
        return f"M{int(m[2]):02}__{m[1]}_{m[3]}"
    if m := re.fullmatch(r"Vcell_z(\d+)", old):
        return f"V_cell_z{m[1]}"
    if m := re.fullmatch(r"Vvas_z(\d+)", old):
        return f"V_vas_z{m[1]}"
    if m := re.fullmatch(r"FLOW_(\w+?)_(res|in|out|z\d+)_(res|in|out|z\d+)", old):
        source, target = (
            f"vas_{site}" if site.startswith("z") else site for site in m.group(2, 3)
        )
        return f"FLOW_{m[1]}_{source}_{target}"
    if m := re.fullmatch(r"(.+)_z(\d+)", old):
        stem, k = m[1], int(m[2])
        if stem in PERFUSATE:
            return f"{stem}_vas_z{k}"
        if stem in {"BUF", "pH"}:
            return f"M{k:02}__{stem}_ext"
        if stem in {"f_alive", "v_nec", "Vbile"} or CELL_PROCESSES.fullmatch(stem):
            return f"M{k:02}__{stem}"
        return old
    return {"f_vas_li": "f_vas", "f_cell_li": "f_cell"}.get(old, old)


def model_ids(path: Path) -> set[str]:
    """The ids of the compartments, species, parameters and reactions of a file."""
    model = libsbml.readSBMLFromFile(str(path)).getModel()
    groups = (
        model.getListOfCompartments(),
        model.getListOfSpecies(),
        model.getListOfParameters(),
        model.getListOfReactions(),
    )
    return {g.get(i).getId() for g in groups for i in range(g.size())}


def main() -> None:
    """Write the map from the rules and the model which is built."""
    before = json.loads(IDS.read_text())
    with TemporaryDirectory() as tmp:
        now = model_ids(build(definition(), Path(tmp)))
    renamed = {old: new for old in before if (new := rename(old)) not in (old, None)}
    kept = {rename(old) or old for old in before}
    out = {
        "before": "98aacca",
        "renamed": renamed,
        "dropped": [old for old in before if rename(old) is None],
        "added": sorted(now - kept - set(renamed.values())),
    }
    MAP.write_text(json.dumps(out, indent=1, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
