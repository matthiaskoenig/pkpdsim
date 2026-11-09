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
from pkpdsim.perfusion.cell_model import constants as cell_constants
from pkpdsim.zonation import ZonationBuilder
from pkpdsim.zonation.composition import CellPorts, compose, write
from pkpdsim.zonation.settings import ZonationSettings

from . import zone_golden
from .conftest import ZONE, load_zone, zone_objects

PARAMS = core.CoreParameters.default("liver_rat")
GOLDEN_CELLS = "Vcell_z1"
GOLDEN_VASCULAR = "Vvas_z1"
"""The compartments of the zone in the golden file, which has the ids of the
zone before the engine named them (`V_cell_z1`, `V_vas_z1`)."""
PORTS = {"Vext", "Vli", "Vbi", "position", "temperature"}
SITES = {"cell", "dead", "ext"}


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
    constants = {p.sid for p in cell_constants(PARAMS)}
    assert _document_ports(path) == PORTS | quantities | constants


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
    """The closed core (1 ml of perfusate around 10 g of cells) does its work
    in an hour: the cells use up the oxygen of the perfusate (0.47 mM to 0), run
    short of ATP (3 mM to 0.9 mM), draw on their glycogen (200 to 192 mM) and
    release glucose (6.5 mM to 11.7 mM, the dilute perfusate takes it), a part of
    the cells dies (living fraction 0.38) and bile is secreted (0.21 ml)."""
    r = load_zone(_written(cell_model(PARAMS), tmp_path / "core.xml"))
    start = {k: r[k] for k in ("[o2_ext]", "[atp_cell]", "[glyc_cell]", "[glc_ext]")}
    assert np.all(np.isfinite(r.simulate(0, 60, 61)))
    assert r["[o2_ext]"] < 0.01 * start["[o2_ext]"]
    assert r["[atp_cell]"] < 0.5 * start["[atp_cell]"]
    assert r["[glyc_cell]"] < start["[glyc_cell]"] - 5.0
    assert r["[glc_ext]"] > start["[glc_ext]"] + 3.0
    assert 0.2 < r["f_alive"] < 0.6
    assert 1e-4 < r["Vbile"] < 4e-4


def test_the_temperature_port_sets_the_arrhenius_factor(tmp_path: Path) -> None:
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
        if t == GOLDEN_CELLS:
            return "Vli"
        if t == GOLDEN_VASCULAR:
            return "Vext"
        if t.endswith(("_z1_cell", "_z1_dead")):
            return t.replace("_z1", "")
        if t.endswith("_z1"):
            stem = t.removesuffix("_z1")
            return f"{stem}_ext" if stem in perfusate else stem
        return t

    return re.sub(r"[A-Za-z_][A-Za-z0-9_]*", rename, text)


def _cell_name(name: str, external: bool) -> str:
    """Get the name of an object of a zone as the standalone core names it.

    The site stays in the name: an object of the perfusate (id `<..>_ext`)
    says `ext` where the zone said `z1`, every other object loses its zone
    (`cells z1` to `cells`, `(z1)` to nothing). A name which says `z1` or `ext`
    in the wrong place then differs.
    """
    if external:
        return name.replace("z1", "ext")
    return name.replace(" (z1)", "").replace(" z1", "")


def _normalise(value: Any, key: str = "", external: bool = False) -> Any:  # noqa: ANN401
    """Rename the ids and names in a description of an object of a zone."""
    if isinstance(value, dict):
        return {k: _normalise(v, k, external) for k, v in value.items()}
    if isinstance(value, list):
        return [_normalise(v, key, external) for v in value]
    if isinstance(value, str):
        if key == "name":
            return _cell_name(value, external)
        return _normalise_text(value)
    return value


def _without_port(described: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """Drop the port of the descriptions (the species of the perfusate have one)."""
    return {
        kind: {
            sid: {a: v for a, v in d.items() if a != "port"} for sid, d in by.items()
        }
        for kind, by in described.items()
    }


def _cell_objects(model: Model) -> dict[str, dict[str, Any]]:
    """The objects of the model, described as the golden file is."""
    return _without_port(
        {
            kind: {
                zone_golden.key(o): zone_golden.describe(o)
                for o in getattr(model, kind)
            }
            for kind in zone_golden.KINDS
        }
    )


def _golden(config: str) -> dict[str, Any]:
    """The objects of a zone before the move, as a standalone core names them."""
    golden = _without_port(zone_golden.load()[config])
    return {
        kind: {
            _normalise_text(sid): _normalise(
                d, external=_normalise_text(sid).endswith("_ext")
            )
            for sid, d in objects.items()
        }
        for kind, objects in golden.items()
    }


def _count(by_kind: dict[str, dict[str, Any]]) -> dict[str, int]:
    return {kind: len(objects) for kind, objects in by_kind.items()}


def _flat_text(text: str) -> str:
    """Rename the ids of a zone in the golden file to those of the engine.

    `Vcell_z1` and `Vvas_z1` are `V_cell_z1` and `V_vas_z1`, a perfusate
    species `glc_z1` is `glc_vas_z1`, an element of the cell model is that of
    the submodel of the zone, `atp_z1_cell` is `M01__atp_cell`, `GLC2PYR_z1`
    is `M01__GLC2PYR` and the buffer and the pH of the zone are those of its
    cell model (`M01__BUF_ext`, `M01__pH_ext`); the factors of a profile
    `p_<process>_z1` are the ones of the engine, the same.
    """
    perfusate = set(ids.PERFUSATE)

    def rename(token: re.Match[str]) -> str:
        t = token.group(0)
        if t == GOLDEN_CELLS:
            return ids.cells(ZONE)
        if t == GOLDEN_VASCULAR:
            return ids.vascular(ZONE)
        for kind in ("cell", "dead"):
            if t.endswith(f"_z1_{kind}"):
                return ids.species(
                    t.removesuffix(f"_{kind}").removesuffix("_z1"), f"z1_{kind}"
                )
        if not t.endswith("_z1") or t.startswith("p_"):
            return t
        stem = t.removesuffix("_z1")
        if stem in perfusate:
            return ids.species(stem, ZONE)
        return {
            "BUF": ids.buffer(ZONE),
            "pH": ids.ph(ZONE),
            "f_alive": ids.alive(ZONE),
            "v_nec": ids.necrosis_rate(ZONE),
            "Vbile": ids.bile(ZONE),
        }.get(stem) or ids.reaction(stem, ZONE)

    return re.sub(r"[A-Za-z_][A-Za-z0-9_]*", rename, text)


def _flat_golden(value: Any, key: str = "") -> Any:  # noqa: ANN401
    """Rename the ids in a description of an object of the golden file."""
    if isinstance(value, dict):
        return {k: _flat_golden(v, k) for k, v in value.items()}
    if isinstance(value, list):
        return [_flat_golden(v, key) for v in value]
    if isinstance(value, str) and key != "name":
        return _flat_text(value)
    return value


@pytest.mark.parametrize("config", list(zone_golden.CONFIGS))
def test_zone_objects_are_the_golden_objects(config: str) -> None:
    """The zone emits what it emitted before the emitter was shared, with the
    ids of the engine (`_flat_text`): every object by all its attributes."""
    objects = zone_objects(
        ZONE,
        ids.cells(ZONE),
        ids.vascular(ZONE),
        PARAMS,
        **zone_golden.CONFIGS[config],
    )
    golden = _flat_golden(
        {
            kind: {_flat_text(sid): d for sid, d in by.items()}
            for kind, by in _without_port(zone_golden.load()[config]).items()
        }
    )
    assert _without_port(zone_golden.describe_objects(objects)) == golden
    if config in _COUNTS:
        assert _count(golden) == _COUNTS[config]


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


def test_the_constants_are_ports_which_the_zones_share(tmp_path: Path) -> None:
    """Every constant of the core and of the buffer is a port, and no other parameter is."""
    model = cell_model(PARAMS)
    path = _written(model, tmp_path / "core.xml")
    constants = {p.sid for p in cell_constants(PARAMS)}
    assert "GLC2PYR_Vmax" in constants
    assert ids.BUFFER_KA in constants
    assert constants <= _document_ports(path)
    assert _document_ports(path) - constants == PORTS | {
        f"{q}_ext" for q in ids.perfusate_quantities(ids.RELEASABLE)
    }


def test_a_process_with_a_zonal_factor_has_the_port_of_it(tmp_path: Path) -> None:
    """The factor `p_<process>` (1) is a port the engine replaces by the one of
    the zone; a process without a profile has neither factor nor port."""
    model = cell_model(PARAMS, processes=("OXPHOS", "BILE"))
    path = _written(model, tmp_path / "core.xml")
    constants = {p.sid for p in cell_constants(PARAMS)}
    assert _document_ports(path) - constants == PORTS | {
        "p_OXPHOS",
        "p_BILE",
    } | {f"{q}_ext" for q in ids.perfusate_quantities(ids.RELEASABLE)}
    factors = {p.sid: p.value for p in model.parameters if p.sid.startswith("p_")}
    assert factors == {"p_OXPHOS": 1.0, "p_BILE": 1.0}
    rates = {r.sid: str(r.formula.math) for r in model.reactions}  # ty: ignore[unresolved-attribute]
    assert "p_OXPHOS" in rates["OXPHOS"]
    assert "p_OXPHOS" not in rates["GLC2PYR"]
    (bile,) = [r for r in model.rate_rules if r.variable == "Vbile"]
    assert "p_BILE" in str(bile.value)


def test_a_process_of_no_core_is_refused() -> None:
    with pytest.raises(ValueError, match="CYP2E1"):
        cell_model(PARAMS, processes=("CYP2E1",))
