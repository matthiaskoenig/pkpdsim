"""The map of the ids of the perfusion organs from before the zonation engine to after.

`tests/perfusion/data/organ_id_map.json` (the liver), `kidney_id_map.json` and
`generic_id_map.json` are the record of what moved when the zones of the liver
went onto the engine (the vascular space, the chain, the cells as a submodel of
every zone, whose states the liver exposes under their old ids) and the ids of the zones of every organ took the names of the engine
(`V_vas_<z>`, `V_cell_<z>`, `<q>_vas_<z>`, the flows between them): every id of
`data/<name>_ids.json` (the model before, commit 98aacca) which has another id in
the model now, the ids the model lost and the ids it gained. The kidney and the
generic organ emit their cells themselves, so the ids of their cells keep their
zone (`GLC2PYR_z1`) and carry no submodel prefix: only the liver has submodels, and the ids of its
cells which the liver observes are the old ones (`atp_z1_cell`, `f_alive_z1`). It is read by `test_organ_reference`,
which checks it against the model that is built and uses it to compare the
trajectories of the two models id by id.

The rules below are written from the id schemes (the ids of the engine,
`pkpdsim.zonation.ids`, and the flat ids of a comp submodel, `M<kk>__<id>` for the elements which the
liver does not expose),
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

from .organ_reference import DATA, DEFINITIONS, ids_file

PERFUSATE = {"o2", "glc", "lac", "h", "co2", "hco3", "k", "alt", "ast", "ldh", "tr"}
"""The quantities of the perfusate of the references (the tracer `tr` included)."""
SUBMODEL = {"organ": True, "kidney": False, "generic": False}
"""The cells of the zones are a submodel of the model (the liver) or are emitted by the builder."""


def map_file(name: str) -> Path:
    """Get the file of the map of a reference."""
    return DATA / f"{name}_id_map.json"


def rename(old: str, submodel: bool = True) -> str | None:
    """Get the id in the model now of an id of the model before, `None` if it is gone."""
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
        return old
    if submodel:
        return {"f_vas_li": "f_vas", "f_cell_li": "f_cell"}.get(old, old)
    return old


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
    """Write the maps from the rules and the models which are built."""
    for name, make in DEFINITIONS.items():
        submodel = SUBMODEL[name]
        before = json.loads(ids_file(name).read_text())
        with TemporaryDirectory() as tmp:
            now = model_ids(build(make(), Path(tmp)))
        renamed = {
            old: new
            for old in before
            if (new := rename(old, submodel)) not in (old, None)
        }
        kept = {rename(old, submodel) or old for old in before}
        out = {
            "before": "98aacca",
            "renamed": renamed,
            "dropped": [old for old in before if rename(old, submodel) is None],
            "added": sorted(now - kept - set(renamed.values())),
        }
        map_file(name).write_text(json.dumps(out, indent=1, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
