"""The reference of the perfusion organs before they moved onto the zonation engine.

`tests/perfusion/data/organ_ids.json` and `organ_trajectories.json` were written
from `builder.compose` and `builder.build` at commit 98aacca (the head of the
branch before `_organ` was rewritten), by `python -m perfusion.organ_reference`
run from `tests/`. They are the witness that the move of the organs onto the
engine changed no id and no trajectory, so they are never regenerated from newer
code: a regeneration would make the tests which read them agree with whatever the
code says.

The definition is `definition()`: the rat liver of ten zones with a profile of the
oxidative phosphorylation and of the bile secretion, a tracer, in the circuit of
the parameter database. The ids are those of the compartments, species,
parameters and reactions of `compose`, sorted. The trajectories are the
simulation of the written model from 0 to 60 min at 13 time points (O0, absolute
tolerance 1e-14, relative 1e-10) of every floating species, the compartment
volumes and the rule parameters.
"""

import json
from pathlib import Path
from tempfile import TemporaryDirectory

from perfusion.conftest import liver_definition, load_circuit
from pkpdsim.perfusion.builder import build, compose, model_ids
from pkpdsim.perfusion.definition import PerfusionDefinition

DATA = Path(__file__).parent / "data"
IDS = DATA / "organ_ids.json"
TRAJECTORIES = DATA / "organ_trajectories.json"
POINTS = 13
END = 60.0


def definition() -> PerfusionDefinition:
    """Get the definition of the reference."""
    return liver_definition(
        10,
        profiles={
            "OXPHOS": [1, 2, 3, 4, 5, 5, 4, 3, 2, 1],
            "BILE": [3, 3, 2, 2, 1, 1, 1, 1, 1, 1],
        },
        tracers=("tr",),
    )


def simulate(path: Path) -> dict[str, list[float]]:
    """Simulate a written model and get its trajectories by id."""
    r = load_circuit(path)
    r.timeCourseSelections = [
        "time",
        *[f"[{s}]" for s in r.model.getFloatingSpeciesIds()],
        *r.model.getCompartmentIds(),
        *r.model.getGlobalParameterIds(),
    ]
    result = r.simulate(0, END, POINTS)
    return {
        name.strip("[]"): [float(v) for v in result[:, k]]
        for k, name in enumerate(r.timeCourseSelections)
    }


def main() -> None:
    """Write the reference files."""
    d = definition()
    IDS.write_text(json.dumps(sorted(model_ids(compose(d))), indent=1) + "\n")
    with TemporaryDirectory() as tmp:
        data = simulate(build(d, Path(tmp)))
    TRAJECTORIES.write_text(json.dumps(data, indent=0) + "\n")


if __name__ == "__main__":
    main()
