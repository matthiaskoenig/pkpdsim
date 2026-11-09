"""The cellular core as a standalone cell model with ports.

The model is what `core.zone_objects` emits for a zone, once, without the zone
in its ids, with the compartments of a zone as ports. What "unchanged" means is
checked against `data/zone_objects.json`, the objects of `zone_objects` before
the emitter was shared (see `zone_golden`): every object is compared by all its
attributes, and the mapping of the objects is equal, not a superset.
"""

import dataclasses
import re
from pathlib import Path
from typing import Any

import libsbml
import numpy as np
import pytest
from pbpk.conftest import compiled, individual
from sbmlutils.factory import PORT_SUFFIX, Model, create_model
from zonation.conftest import mass_g

from pkpdsim.pbpk.substance import Substance
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Value
from pkpdsim.perfusion import core, ids
from pkpdsim.perfusion.cell_model import cell_model
from pkpdsim.zonation import ZonationBuilder
from pkpdsim.zonation.composition import CellPorts, compose, write
from pkpdsim.zonation.settings import ZonationSettings

from . import zone_golden
from .conftest import CELLS, VASCULAR, ZONE, load_zone

PARAMS = core.CoreParameters.default("liver_rat")
PORTS = {"Vext", "Vli", "Vbi", "position", "temperature"}
SITES = {"cell", "dead", "ext"}
SUFFIX = re.compile(r"\bz1\b|\bext\b")


def _written(model: Model, path: Path) -> Path:
    """Write the model, validated."""
    result = create_model(model, filepath=path, sbml_level=3, sbml_version=2)
    assert result.validation is not None
    assert result.validation.error_count == 0
    assert result.validation.warning_count == 0
    return path


def _document_ports(path: Path) -> set[str]:
    """The ids of the ports of the written model, without their suffix.

    Read from the document on purpose: `port=True` of sbmlutils creates the
    ports when the model is written and leaves `Model.ports` empty, so a test
    against `m.ports` finds no port.
    """
    document = libsbml.readSBMLFromFile(str(path))  # keeps the model alive
    plugin = document.getModel().getPlugin("comp")
    return {
        plugin.getPort(k).getId().removesuffix(PORT_SUFFIX)
        for k in range(plugin.getNumPorts())
    }


def test_the_core_is_a_standalone_model_with_the_cell_ports(tmp_path: Path) -> None:
    """The ports are those of `CellPorts` and one concentration port per
    quantity of the perfusate."""
    path = _written(cell_model(PARAMS), tmp_path / "core.xml")
    quantities = {
        ids.species(q, "ext") for q in ids.perfusate_quantities(ids.RELEASABLE)
    }
    assert _document_ports(path) == PORTS | quantities


def test_the_ports_are_the_cell_ports_of_the_engine() -> None:
    """The ids are what `compose` replaces, not a second spelling of them."""
    ports = CellPorts()
    assert {
        ports.external,
        ports.cell,
        ports.bile,
        ports.position,
        ports.temperature,
    } == PORTS
    assert ports.concentration.format(sid="glc") == ids.species("glc", "ext")


def test_the_standalone_core_simulates_on_its_own(tmp_path: Path) -> None:
    """The core runs without a circuit: finite, the living cells use their ATP
    and oxygen, bile is secreted."""
    r = load_zone(_written(cell_model(PARAMS), tmp_path / "core.xml"))
    result = r.simulate(0, 60, 61)
    assert np.all(np.isfinite(result))
    assert r["[glc_ext]"] > 0
    assert r["Vbile"] > 0
    assert 0 < r["f_alive"] <= 1


def test_the_core_is_per_temperature_of_its_port(tmp_path: Path) -> None:
    """The temperature port sets the Arrhenius factor, warmer is faster."""
    path = _written(cell_model(PARAMS), tmp_path / "core.xml")
    warm = individual(path, temperature=37.0 + 5.0)
    cold = individual(path, temperature=37.0 - 5.0)
    reference = individual(path, temperature=37.0)
    assert reference["f_T"] == pytest.approx(1.0, rel=1e-9)
    assert warm["f_T"] > 1.0 > cold["f_T"]


@pytest.mark.parametrize("config", ["plain", "reduced"])
def test_the_standalone_core_is_valid_sbml(tmp_path: Path, config: str) -> None:
    """libsbml consistency checks report no error and no warning."""
    model = cell_model(PARAMS, **zone_golden.CONFIGS[config])
    document = libsbml.readSBMLFromFile(str(_written(model, tmp_path / "core.xml")))
    assert document.checkConsistency() == 0
    assert document.getNumErrors() == 0


def _normalise_text(text: str) -> str:
    """Rename the ids of a zone in a text to those of the standalone core."""
    perfusate = set(ids.PERFUSATE) | {"pH", "BUF"}

    def rename(token: re.Match[str]) -> str:
        t = token.group(0)
        if t == CELLS:
            return "Vli"
        if t == VASCULAR:
            return "Vext"
        if t.endswith(("_z1_cell", "_z1_dead")):
            return t.replace("_z1", "")
        if t.endswith("_z1"):
            stem = t.removesuffix("_z1")
            return f"{stem}_ext" if stem in perfusate else stem
        return t

    return re.sub(r"[A-Za-z_][A-Za-z0-9_]*", rename, text)


def _normalise(value: Any, key: str = "") -> Any:  # noqa: ANN401
    """Rename the ids in a description of an object of a zone; a name loses
    its sites (`z1`, `ext`) on both sides of the comparison."""
    if isinstance(value, dict):
        return {k: _normalise(v, k) for k, v in value.items()}
    if isinstance(value, list):
        return [_normalise(v, key) for v in value]
    if isinstance(value, str):
        if key == "name":
            return _site_free(value)
        return _normalise_text(value)
    return value


def _site_free(name: str) -> str:
    """Get a name without its sites."""
    bare = re.sub(r"\s*\(\s*\)", "", SUFFIX.sub("", name))
    return " ".join(re.sub(r"\s+\)", ")", bare).split())


def _without_port(described: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """Drop the port of the descriptions (the species of the perfusate have one)."""
    return {
        kind: {
            sid: {a: v for a, v in d.items() if a != "port"} for sid, d in by.items()
        }
        for kind, by in described.items()
    }


def _named(described: dict[str, Any]) -> dict[str, Any]:
    """The description with a name without sites."""
    return {
        a: _site_free(v) if a == "name" and isinstance(v, str) else v
        for a, v in described.items()
    }


def _cell_objects(model: Model) -> dict[str, dict[str, Any]]:
    """The objects of the model, described as the golden file is."""
    return _without_port(
        {
            kind: {
                zone_golden.key(o): _named(zone_golden.describe(o))
                for o in getattr(model, kind)
            }
            for kind in zone_golden.KINDS
        }
    )


def _golden(config: str) -> dict[str, Any]:
    """The objects of a zone before the move, as a standalone core names them."""
    golden = _without_port(zone_golden.load()[config])
    return {
        kind: {_normalise_text(sid): _normalise(d) for sid, d in objects.items()}
        for kind, objects in golden.items()
    }


def _count(by_kind: dict[str, dict[str, Any]]) -> dict[str, int]:
    return {kind: len(objects) for kind, objects in by_kind.items()}


@pytest.mark.parametrize("config", list(zone_golden.CONFIGS))
def test_zone_objects_are_the_golden_objects(config: str) -> None:
    """The zone emits what it emitted before the emitter was shared."""
    objects = core.zone_objects(
        ZONE, CELLS, VASCULAR, PARAMS, **zone_golden.CONFIGS[config]
    )
    assert _without_port(zone_golden.describe_objects(objects)) == _without_port(
        zone_golden.load()[config]
    )


EXTRA_PARAMETERS = {"position", "temperature", "T_ref", "R_gas", "Ea", "f_T"}
"""The parameters of the standalone core which no zone emits: the ports of the
position and temperature and the Arrhenius factor, which the circuit defines
for the zones."""


@pytest.mark.parametrize("config", ["plain", "reduced"])
def test_the_core_is_the_objects_of_a_zone(config: str) -> None:
    """Every object of the core equals the object of the zone, by all its
    attributes, but its id; the exceptions are the ports and the constants."""
    options = zone_golden.CONFIGS[config]
    model = cell_model(PARAMS, **options)
    found = _cell_objects(model)
    expected = _golden(config)
    constants = {
        p.sid
        for p in (
            *core.core_parameters(PARAMS, **options),
            *core.buffer_parameters(),
        )
    }
    extra = EXTRA_PARAMETERS | constants
    assert extra <= set(found["parameters"])
    for kind in zone_golden.KINDS:
        got = {
            sid: d
            for sid, d in found[kind].items()
            if not (kind == "parameters" and sid in extra)
            and not (kind == "assignments" and sid == "f_T")
        }
        assert got == expected[kind], kind
    assert _count(expected) == _COUNTS[config]


_COUNTS = {
    "plain": {
        "species": 32,
        "reactions": 22,
        "parameters": 4,
        "rules": 2,
        "rate_rules": 2,
        "assignments": 20,
    },
    "reduced": {
        "species": 26,
        "reactions": 18,
        "parameters": 3,
        "rules": 2,
        "rate_rules": 1,
        "assignments": 18,
    },
}
"""The objects compared by kind (82 and 68), so an empty comparison fails."""


def test_without_bile_the_bile_compartment_is_still_a_port(tmp_path: Path) -> None:
    """`compose` wires `Vbi` only with a bile, the port is there either way."""
    model = cell_model(PARAMS, bile=False)
    path = _written(model, tmp_path / "core.xml")
    assert _document_ports(path) >= PORTS
    assert "Vbile" not in {p.sid for p in model.parameters}


def test_the_engine_places_the_core_once_per_zone(tmp_path: Path) -> None:
    """Composed into four zones the cell model is wired to the space of Disse
    and flattens to a model which simulates: one core per zone."""
    create_model(
        cell_model(PARAMS),
        filepath=tmp_path / "perfusion_core.xml",
        sbml_level=3,
        sbml_version=2,
    )
    glucose = Substance("glc", "glucose", Value(180.0, U.g_per_mole))
    builder = ZonationBuilder(
        ZonationSettings(n_zones=4, liver_mass=mass_g(1000.0)), [glucose], None
    )
    flat = write(
        compose(builder, "perfusion_core", CellPorts(), deletions=()), tmp_path
    )
    r = compiled(flat)
    atp = [k for k in r.model.getFloatingSpeciesIds() if k.endswith("__atp_cell")]
    assert len(atp) == 4
    local = set(r.model.getFloatingSpeciesIds())
    assert not {k for k in local if k.endswith("__glc_ext")}, "replaced by Disse"
    assert {f"M{k:02}__o2_ext" for k in range(1, 5)} <= local, "not a substance"
    assert np.all(np.isfinite(r.simulate(0, 10, 11)))


def test_potassium_which_is_not_released_stays_in_the_cells(tmp_path: Path) -> None:
    """The cells hold potassium whatever the organ releases (`k` is no enzyme),
    so its content is a constant of the model: valid SBML for `released=("alt",)`
    in a zone and in the standalone core (before, `K_content` was undefined)."""
    released = ("alt",)
    constants = {p.sid for p in core.core_parameters(PARAMS, released=released)}
    assert "K_content" in constants
    assert {"AST_content", "LDH_content"}.isdisjoint(constants)
    model = cell_model(PARAMS, released=released)
    document = libsbml.readSBMLFromFile(str(_written(model, tmp_path / "core.xml")))
    assert document.checkConsistency() == 0
    assert document.getNumErrors() == 0
    assert "k_dead" not in {s.sid for s in model.species}


def test_the_content_of_potassium_is_needed_without_its_release() -> None:
    """A core without the content of potassium is refused for any organ."""
    params = dataclasses.replace(PARAMS, K_content=None)
    assert any("K_content" in p for p in params.problems((), bile=False))
