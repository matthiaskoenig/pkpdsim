"""The perfusion liver on the zonation engine is the liver from before it.

The ids and the trajectories of `data/organ_ids.json` and
`data/organ_trajectories.json` (the model at commit 98aacca, before the zones
of the liver were the engine's; see `organ_reference`) are compared with the
model that is built now through `data/organ_id_map.json`, the record of what
moved (see `organ_id_map`). Both claims are equalities: the ids of the model
are exactly the images of the ids before without the ones it lost, and the
ones it gained are exactly the listed ones; every trajectory of the model
before is the trajectory of its image.
"""

import json
from pathlib import Path

import numpy as np
import pytest

from pkpdsim.perfusion.builder import build, compose, model_ids

from .organ_id_map import MAP, rename
from .organ_id_map import model_ids as file_ids
from .organ_reference import IDS, TRAJECTORIES, definition, simulate

TOLERANCE = 1e-6
"""The largest relative difference of a trajectory."""
FLOOR = 1e-10
"""Below this a value is zero: a relative difference needs a denominator. It is
1e4 times the absolute tolerance of the integrator (1e-14, `load_circuit`): the
pool of NADP is depleted to 1e-12 mM and differs by 1.6e-17 mM between the two
models, a relative difference of 1.6e-5 of that value and 5e-16 of its maximum."""


@pytest.fixture(scope="module")
def built(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The flat model of the reference definition."""
    return build(definition(), tmp_path_factory.mktemp("organ"))


def _map() -> dict[str, object]:
    return json.loads(MAP.read_text())


def _before() -> set[str]:
    return set(json.loads(IDS.read_text()))


def _renamed() -> dict[str, str]:
    renamed = _map()["renamed"]
    assert isinstance(renamed, dict)
    return renamed


def _dropped() -> set[str]:
    dropped = _map()["dropped"]
    assert isinstance(dropped, list)
    return set(dropped)


def _images() -> set[str]:
    return {_renamed().get(old, old) for old in _before() - _dropped()}


def test_the_map_is_what_the_rules_say() -> None:
    """The record of the map is the rules of `organ_id_map`, id by id."""
    before = _before()
    assert set(_renamed()) <= before
    for old in before:
        expected = rename(old)
        if expected is None:
            assert old in _dropped()
        else:
            assert _renamed().get(old, old) == expected, old
    assert not _dropped() - before


def test_the_ids_are_the_ids_before_through_the_map(built: Path) -> None:
    """Equality, not inclusion: nothing is gained but the listed, nothing lost."""
    now = file_ids(built)
    added = set(_map()["added"])  # ty: ignore[invalid-argument-type]
    assert len(set(_renamed().values())) == len(_renamed()), "the map is injective"
    assert len(_images()) == len(_before() - _dropped()), "no two ids became one"
    assert added.isdisjoint(_images())
    assert _images() == now - added
    assert added <= now


def test_the_composed_model_has_the_ids_of_the_engine_and_the_circuit() -> None:
    """`compose` is the model of the liver before the cells are flat: the ids
    of the circuit and of the zones' vascular space and cells, no cell model."""
    ids_ = model_ids(compose(definition()))
    assert "V_vas_z3" in ids_
    assert "glc_vas_z3" in ids_
    assert "M03__atp_cell" not in ids_, "the cells are in the submodel"
    assert "Vvas_z3" not in ids_


def test_the_trajectories_are_the_trajectories_before_through_the_map(
    built: Path,
) -> None:
    """The relative difference of every series of the model before and of its image."""
    reference = json.loads(TRAJECTORIES.read_text())
    now = simulate(built)
    dropped = _dropped()
    worst, where = 0.0, ""
    absolute = 0.0
    compared = 0
    for old, before in reference.items():
        if old in dropped:
            continue
        new = _renamed().get(old, old)
        assert new in now, f"{old} -> {new} is not simulated"
        a, b = np.asarray(before), np.asarray(now[new])
        scale = np.maximum(np.maximum(np.abs(a), np.abs(b)), FLOOR)
        difference = float(np.max(np.abs(a - b) / scale))
        compared += 1
        absolute = max(absolute, float(np.max(np.abs(a - b))))
        if difference > worst:
            worst, where = difference, f"{old} -> {new}"
    print(
        f"maximum relative difference {worst:.3e} at {where}, maximum absolute "
        f"difference {absolute:.3e}, {compared} series"
    )
    assert compared == len(set(reference) - dropped) > 600
    assert worst < TOLERANCE, f"{worst:.3e} at {where}"
