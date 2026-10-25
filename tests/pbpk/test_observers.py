"""Tests of the observers: assignment rules of the top model, no states.

The observers read the pool species and the outflow concentrations of the
organs: the unbound and the blood concentration of every pool, the hepatic
vein and the arm vein, the flow-weighted mix of the venous outflows of the
organs of a preset of `data/arm_vein.toml`.
"""

import dataclasses
from pathlib import Path
from typing import TYPE_CHECKING

import libsbml
import numpy as np
import pytest

from pkpdsim.pbpk import (
    Erythrocytes,
    Organ,
    Part,
    PBPKDefinition,
    Pool,
    Route,
    ids,
    load,
)
from pkpdsim.pbpk.builder import build
from pkpdsim.pbpk.observers import arm_vein_weights

from .conftest import definition, individual, simple_drug

if TYPE_CHECKING:
    import roadrunner

S = "drg"
"""The substance of the tests."""

TIMES = np.linspace(0.0, 240.0, 49)
"""The time points of the curves (min)."""


def _run(path: Path, cols: list[str], **changes: float) -> np.ndarray:
    """Simulate an IV bolus of 10 mg, the columns without time."""
    r = individual(path, **changes)
    r.integrator.absolute_tolerance = 1e-12
    r.integrator.relative_tolerance = 1e-10
    r[ids.iv_dose(S)] = 10.0
    return r.simulate(times=TIMES, selections=["time", *cols])[:, 1:]


def _build(tmp_path: Path, **overrides: object) -> Path:
    """Build the simple definition with IV route only."""
    d = definition(routes=[Route.iv(S)], **overrides)
    return build(d, tmp_path).flat_path


def _rr(path: Path) -> "roadrunner.RoadRunner":
    """Load a model without changes."""
    return individual(path)


def test_arm_vein_weighted_outflow(tmp_path: Path) -> None:
    """`Carm` is the PK-Sim default mix `0.713*Cout_sk + 0.287*Cout_mu`."""
    path = _build(tmp_path)
    sk, mu = ids.organ_outflow(Organ.SKIN, S), ids.organ_outflow(Organ.MUSCLE, S)
    carm, csk, cmu = _run(path, [ids.arm(S), sk, mu]).T
    np.testing.assert_allclose(carm, 0.713 * csk + 0.287 * cmu, rtol=1e-12, atol=0)
    assert np.max(carm) > 0


@pytest.mark.parametrize("preset", ["pksim_default", "levitt2004", "musther2015"])
def test_arm_vein_presets(tmp_path: Path, preset: str) -> None:
    """Every preset is a mix of outflows (weights sum to 1) which `Carm` follows."""
    weights = arm_vein_weights(preset)
    assert sum(w.value for w in weights.values()) == pytest.approx(1.0, abs=1e-12)
    path = _build(tmp_path, arm_vein=preset)
    cols = [
        f"[{ids.pool_species(site, S)}]"
        if isinstance(site, Pool)
        else ids.organ_outflow(site, S)
        for site in weights
    ]
    result = _run(path, [ids.arm(S), *cols])
    expected = sum(w.value * result[:, i + 1] for i, w in enumerate(weights.values()))
    np.testing.assert_allclose(result[:, 0], expected, rtol=1e-12, atol=0)


def test_arterialised_is_arterial(tmp_path: Path) -> None:
    """The arterialised hand vein is the arterial plasma."""
    path = _build(tmp_path, arm_vein="arterialised")
    carm, car = _run(path, [ids.arm(S), f"[{ids.pool_species(Pool.ARTERIAL, S)}]"]).T
    np.testing.assert_allclose(carm, car, rtol=1e-12, atol=0)


def test_arm_vein_weights_are_parameters(tmp_path: Path) -> None:
    """The weights are parameters with their source, which a model may change."""
    path = _build(tmp_path)
    r = _rr(path)
    assert r[ids.arm_weight(Organ.SKIN)] == 0.713
    assert r[ids.arm_weight(Organ.MUSCLE)] == 0.287
    doc = libsbml.readSBMLFromFile(str(path))
    notes = doc.getModel().getParameter(ids.arm_weight(Organ.SKIN)).getNotesString()
    assert "PK-Sim" in notes


def test_unknown_arm_vein_preset_raises() -> None:
    """The preset is validated with the definition."""
    with pytest.raises(ValueError, match="arm_vein 'nope'"):
        definition(arm_vein="nope").validate()


def test_arm_vein_is_observer(tmp_path: Path) -> None:
    """The observer changes no state: the curves are the same with every preset.

    No species, reaction or rate rule reads an observer; it is an assignment
    rule which only other observers read.
    """
    species = [
        f"[{ids.pool_species(Pool.VENOUS, S)}]",
        "[LI__Cti_drg]",
        "[MU__Cti_drg]",
    ]
    curves = [
        _run(_build(tmp_path / preset, arm_vein=preset), species)
        for preset in ("pksim_default", "arterialised")
    ]
    np.testing.assert_array_equal(curves[0], curves[1])

    path = _build(tmp_path / "check")
    model = libsbml.readSBMLFromFile(str(path)).getModel()
    observers = {ids.arm(S), ids.arm_free(S), ids.arm_blood(S), ids.hepatic_vein(S)}
    observers |= {f(p, S) for p in Pool for f in (ids.free, ids.blood)}
    assert all(model.getSpecies(o) is None for o in observers)
    readers = [
        libsbml.formulaToL3String(r.getKineticLaw().getMath())
        for r in model.getListOfReactions()
    ] + [
        libsbml.formulaToL3String(rule.getMath())
        for rule in model.getListOfRules()
        if rule.getVariable() not in observers
    ]
    for formula in readers:
        names = set(formula.replace("(", " ").replace(")", " ").split())
        assert not names & observers, formula


@pytest.mark.parametrize(
    "ery",
    [
        None,
        Erythrocytes("equilibrium", bp=1.3),
        Erythrocytes("kinetic", kbc=2.0, k_exchange=1.0),
    ],
    ids=["plasma", "equilibrium", "kinetic"],
)
def test_free_and_blood_observers(tmp_path: Path, ery: Erythrocytes | None) -> None:
    """Unbound `fu*C` and blood `(1-HCT)*C + HCT*Cbc` of the pools and the arm vein.

    Without kinetic blood cells the blood cells hold `Kbc*C`, the blood
    `BP*C` (`BP = 1 - HCT` without erythrocyte distribution); kinetic blood
    cells are their species, at the arm vein the flow-weighted mix of the
    blood cell outflows of skin and muscle.
    """
    drug = dataclasses.replace(simple_drug(), erythrocytes=ery)
    path = _build(tmp_path, substances=[drug, *definition().substances[1:]])
    r = _rr(path)
    fu, bp, hct = r[ids.fu(S)], r[ids.bp(S)], r[ids.HCT]
    for pool in Pool:
        c = f"[{ids.pool_species(pool, S)}]"
        cols = [c, ids.free(pool, S), ids.blood(pool, S)]
        if ery is not None and ery.mode == "kinetic":
            cols.append(f"[{ids.pool_species(pool, S, part=Part.BLOOD_CELLS)}]")
        result = _run(path, cols)
        np.testing.assert_allclose(result[:, 1], fu * result[:, 0], rtol=1e-12)
        if ery is not None and ery.mode == "kinetic":
            blood = (1 - hct) * result[:, 0] + hct * result[:, 3]
        else:
            blood = bp * result[:, 0]
        np.testing.assert_allclose(result[:, 2], blood, rtol=1e-12)

    cols = [ids.arm(S), ids.arm_free(S), ids.arm_blood(S)]
    if ery is not None and ery.mode == "kinetic":
        cols += [
            ids.organ_outflow_bc(Organ.SKIN, S),
            ids.organ_outflow_bc(Organ.MUSCLE, S),
        ]
    result = _run(path, cols)
    np.testing.assert_allclose(result[:, 1], fu * result[:, 0], rtol=1e-12)
    if ery is not None and ery.mode == "kinetic":
        cbc = 0.713 * result[:, 3] + 0.287 * result[:, 4]
        blood = (1 - hct) * result[:, 0] + hct * cbc
        assert np.max(cbc) > 0
    else:
        blood = bp * result[:, 0]
    np.testing.assert_allclose(result[:, 2], blood, rtol=1e-12)


def test_arm_vein_data_sources() -> None:
    """Every preset of the data file names its sources; the weights are fractions."""
    data = load("arm_vein")
    assert {"pksim_default", "levitt2004", "musther2015", "arterialised"} <= set(data)
    for preset in data:
        for w in arm_vein_weights(preset).values():
            assert 0 < w.value <= 1
            assert w.source


def test_definition_arm_vein_default() -> None:
    """The PK-Sim default is the preset of a definition which names none."""
    assert PBPKDefinition(sid="x", substances=[]).arm_vein == "pksim_default"
