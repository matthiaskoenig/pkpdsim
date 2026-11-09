"""The perfusion organs on the zonation engine are the organs from before it.

The ids and the trajectories of `data/<name>_ids.json` and
`data/<name>_trajectories.json` (the liver, the kidney and a generic organ at
commit 98aacca, before the zones of the liver were the engine's; see
`organ_reference`) are compared with the model that is built now through
`data/<name>_id_map.json`, the record of what moved (see `organ_id_map`). Both
claims are equalities: the ids of the model are exactly the images of the ids
before without the ones it lost, and the ones it gained are exactly the listed
ones; every trajectory of the model before, the concentrations, the volumes, the
parameters and the flux of every reaction, is the trajectory of its image.
"""

import json
from pathlib import Path

import libsbml
import numpy as np
import pytest

from pkpdsim.perfusion.builder import build, compose, model_ids

from .organ_id_map import SUBMODEL, map_file, rename
from .organ_id_map import model_ids as file_ids
from .organ_reference import DEFINITIONS, ids_file, simulate, trajectories_file

RTOL = 1e-6
"""The relative difference of a trajectory."""
ATOL = 1e-12
"""The absolute difference which is no difference, 1/100 of the absolute tolerance
of the integrator (1e-14, `load_circuit`) times 1e4: the pool of NADP is depleted to
1e-12 mM and differs by 1.6e-17 mM between the models, which is noise. A series
agrees when `|a - b| <= RTOL * max(|a|, |b|) + ATOL`."""

NAMES = list(DEFINITIONS)


@pytest.fixture(scope="module", params=NAMES)
def organ(
    request: pytest.FixtureRequest, tmp_path_factory: pytest.TempPathFactory
) -> tuple[str, Path]:
    """The name of a reference and its model as it is built now."""
    name: str = request.param
    path = build(DEFINITIONS[name](), tmp_path_factory.mktemp(name))
    return name, path


def _map(name: str) -> dict[str, object]:
    return json.loads(map_file(name).read_text())


def _renamed(name: str) -> dict[str, str]:
    renamed = _map(name)["renamed"]
    assert isinstance(renamed, dict)
    return renamed


def _dropped(name: str) -> set[str]:
    dropped = _map(name)["dropped"]
    assert isinstance(dropped, list)
    return set(dropped)


def _before(name: str) -> set[str]:
    return set(json.loads(ids_file(name).read_text()))


def _images(name: str) -> set[str]:
    return {_renamed(name).get(old, old) for old in _before(name) - _dropped(name)}


def test_the_map_is_what_the_rules_say() -> None:
    """The record of the map is the rules of `organ_id_map`, id by id."""
    for name in NAMES:
        before = _before(name)
        assert set(_renamed(name)) <= before
        for old in before:
            expected = rename(old, SUBMODEL[name])
            if expected is None:
                assert old in _dropped(name)
            else:
                assert _renamed(name).get(old, old) == expected, (name, old)
        assert not _dropped(name) - before


def test_the_ids_are_the_ids_before_through_the_map(
    organ: tuple[str, Path],
) -> None:
    """Equality, not inclusion: nothing is gained but the listed, nothing lost."""
    name, path = organ
    now = file_ids(path)
    added = set(_map(name)["added"])  # ty: ignore[invalid-argument-type]
    renamed = _renamed(name)
    assert len(set(renamed.values())) == len(renamed), "the map is injective"
    assert len(_images(name)) == len(_before(name) - _dropped(name)), (
        "no two ids became one"
    )
    assert added.isdisjoint(_images(name))
    assert _images(name) == now - added
    assert added <= now


def test_only_the_liver_has_the_ids_of_a_submodel(organ: tuple[str, Path]) -> None:
    """An id `M<kk>__<id>` is that of an element of a submodel: the kidney and a
    generic organ have none, and no comp plugin."""
    name, path = organ
    prefixed = {sid for sid in file_ids(path) if sid.startswith("M") and "__" in sid}
    assert bool(prefixed) is SUBMODEL[name]
    model = libsbml.readSBMLFromFile(str(path)).getModel()
    assert model.getPlugin("comp") is None, "the model which is simulated is flat"


def test_the_composed_liver_has_the_ids_of_the_engine_and_the_circuit() -> None:
    """`compose` is the model of the liver before the cells are flat: the ids
    of the circuit and of the zones' vascular space and cells, no cell model."""
    ids_ = model_ids(compose(DEFINITIONS["organ"]()))
    assert "V_vas_z3" in ids_
    assert "glc_vas_z3" in ids_
    assert "M03__atp_cell" not in ids_, "the cells are in the submodel"
    assert "Vvas_z3" not in ids_


def test_the_trajectories_are_the_trajectories_before_through_the_map(
    organ: tuple[str, Path],
) -> None:
    """Every series of the model before, its fluxes included, and its image agree."""
    name, path = organ
    reference = json.loads(trajectories_file(name).read_text())
    now = simulate(path)
    document = libsbml.readSBMLFromFile(str(path))  # keeps the model alive
    model = document.getModel()
    reactions = {model.getReaction(k).getId() for k in range(model.getNumReactions())}
    dropped = _dropped(name)
    worst = {"state": (0.0, ""), "flux": (0.0, "")}
    compared = {"state": 0, "flux": 0}
    absolute = 0.0
    for old, before in reference.items():
        if old in dropped:
            continue
        new = _renamed(name).get(old, old)
        assert new in now, f"{old} -> {new} is not simulated"
        a, b = np.asarray(before), np.asarray(now[new])
        scale = np.maximum(np.abs(a), np.abs(b))
        # |a - b| <= RTOL * scale + ATOL, as a relative difference under the criterion
        difference = float(np.max(np.abs(a - b) / (scale + ATOL / RTOL)))
        kind = "flux" if new in reactions else "state"
        compared[kind] += 1
        absolute = max(absolute, float(np.max(np.abs(a - b))))
        if difference > worst[kind][0]:
            worst[kind] = (difference, f"{old} -> {new}")
    print(
        f"{name}: maximum relative difference of the states {worst['state'][0]:.3e} "
        f"at {worst['state'][1]} ({compared['state']} series), of the fluxes "
        f"{worst['flux'][0]:.3e} at {worst['flux'][1]} ({compared['flux']} series), "
        f"maximum absolute difference {absolute:.3e}"
    )
    assert sum(compared.values()) == len(set(reference) - dropped)
    assert compared["state"] > 100
    assert compared["flux"] > 10, "the fluxes of the reactions are compared"
    for kind, (difference, where) in worst.items():
        assert difference <= RTOL, f"{kind}: {difference:.3e} at {where}"
