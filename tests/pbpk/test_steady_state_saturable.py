"""The kidney at the steady state of a saturable secretion (`steady_state.Saturable`).

An endogenous substance is produced in the muscle tissue and excreted only by
the kidney, filtered with a reabsorbed fraction and secreted from the kidney
plasma with Michaelis-Menten kinetics (a custom reaction into the urine), as
creatinine is. The kidney plasma of the elimination of the production is the
root of a quadratic (`kidney_plasma(..., eliminated=..., secretion=...)`), the
baseline the arterial plasma of that kidney plasma (`kidney_arterial` with the
elimination, which reads the kidney plasma once): with
them as starts every species stays at its start in the simulation, for the
reference individual and changed ones.
"""

import dataclasses
import math
from pathlib import Path
from typing import TYPE_CHECKING

import libsbml
import numpy as np
import pytest

from pkpdsim.pbpk import (
    Binding,
    CustomReaction,
    Erythrocytes,
    Excretion,
    Formula,
    Kp,
    ModelParameter,
    Organ,
    Part,
    PBPKDefinition,
    Substance,
    U,
    Variants,
    build,
    ids,
)
from pkpdsim.pbpk import steady_state as ss

from .conftest import individual

if TYPE_CHECKING:
    import roadrunner

S = "end"
DAY = 24 * 60.0
STANDARD = Variants()
SECRETION = ss.Saturable(vmax="V_sec", km="km_sec")

PARAMETERS = [
    ModelParameter("k_prod", 0.009, U.mmole_per_min),
    ModelParameter("V_sec", 0.002, U.mmole_per_min),
    ModelParameter("km_sec", 0.08, U.mM),
]


def _substance(erythrocytes: Erythrocytes | None = None) -> Substance:
    """The substance, its baseline the arterial plasma of the kidney plasma `C_ki_0`."""
    s = Substance(
        S,
        name="endogenous substance",
        mr=113.1,
        binding=Binding(fu=0.8, protein="albumin"),
        erythrocytes=erythrocytes,
        distribution=Kp(0.6, ftissue=6.0),
        excretion=Excretion(urine=True, fraction_reabsorbed=0.1),
    )
    baseline = ss.kidney_arterial(s, "C_ki_0", eliminated="k_prod", variants=STANDARD)
    return dataclasses.replace(s, baseline=Formula(baseline, U.mM))


def _definition(s: Substance) -> PBPKDefinition:
    """The production of the muscle, the secretion of the kidney, their steady state."""
    kidney = ss.kidney_plasma(
        s, eliminated="k_prod", secretion=SECRETION, variants=STANDARD
    )
    muscle = ss.organ_plasma(s, Organ.MUSCLE, rate="k_prod", variants=STANDARD)
    return PBPKDefinition(
        sid="saturable",
        substances=[s],
        chemistry=[
            CustomReaction(
                "PROD", Organ.MUSCLE, Part.TISSUE, f"-> Cti_{S}", "k_prod", []
            ),
            CustomReaction(
                "SATSEC",
                Organ.KIDNEY,
                Part.PLASMA,
                f"Cpl_{S} -> Aurine_{S}",
                f"V_sec * Cpl_{S} / (km_sec + Cpl_{S})",
            ),
        ],
        parameters=[
            *PARAMETERS,
            ModelParameter("C_ki_0", Formula(kidney, U.mM), U.mM),
        ],
        starts={
            f"MU__Cpl_{S}": Formula(muscle, U.mM),
            f"MU__Cti_{S}": Formula(
                ss.organ_tissue(
                    s, Organ.MUSCLE, muscle, rate="k_prod", variants=STANDARD
                ),
                U.mM,
            ),
            f"KI__Cpl_{S}": Formula("C_ki_0", U.mM),
            f"KI__Cti_{S}": Formula(
                ss.organ_tissue(s, Organ.KIDNEY, "C_ki_0", variants=STANDARD), U.mM
            ),
        },
    )


def _drift(r: "roadrunner.RoadRunner") -> dict[str, float]:
    """The largest relative change of every species but the urine over a day."""
    sids = [s for s in r.model.getFloatingSpeciesIds() if s != ids.urine(S)]
    start = np.array([r[f"[{sid}]"] for sid in sids])
    values = np.asarray(r.simulate(0, DAY, 97, selections=[f"[{sid}]" for sid in sids]))
    scale = np.where(start > 0, start, 1.0)
    return dict(zip(sids, np.max(np.abs(values - start), axis=0) / scale, strict=True))


def _evaluate(formula: str, values: dict[str, float]) -> float:
    """Evaluate a formula with numbers for its ids."""
    ast = libsbml.parseL3Formula(formula.replace(" dimensionless", ""))
    names: set[str] = set()
    stack = [ast]
    while stack:
        node = stack.pop()
        if node.getType() == libsbml.AST_NAME:
            names.add(node.getName())
        stack.extend(node.getChild(k) for k in range(node.getNumChildren()))
    expression = libsbml.formulaToL3String(ast).replace("^", "**")
    namespace = {"sqrt": math.sqrt}
    return float(eval(expression, namespace, {n: values[n] for n in names}))


@pytest.mark.parametrize(
    "changes",
    [{}, {"BW": 110.0, "HEIGHT": 190.0}, {"f_renal_function": 0.3}],
    ids=["reference", "body", "renal function"],
)
@pytest.mark.parametrize(
    "erythrocytes",
    [None, Erythrocytes("equilibrium", bp=0.8)],
    ids=["plasma", "erythrocytes"],
)
def test_the_saturable_secretion_starts_stationary(
    tmp_path: Path, changes: dict[str, float], erythrocytes: Erythrocytes | None
) -> None:
    """Every species stays within 1e-9 of its start; the urine takes the production."""
    s = _substance(erythrocytes)
    r = individual(build(_definition(s), tmp_path).flat_path, **changes)
    r.integrator.absolute_tolerance = 1e-14
    r.integrator.relative_tolerance = 1e-12
    moved = {sid: d for sid, d in _drift(r).items() if d > 1e-9}
    assert not moved, moved
    r.reset()
    car, cki = r[f"[Car_{S}]"], r[f"[KI__Cpl_{S}]"]
    assert car == pytest.approx(r[ids.baseline(S)], rel=1e-12)
    # the secretion saturates: the kidney plasma is of the order of km
    assert 0.1 < cki / r["km_sec"] < 10
    assert car > cki
    # the kidney plasma of the arterial form is that of the elimination, the
    # arterial plasma of the secretion that of the elimination
    arterial = ss.kidney_plasma(
        s, arterial=ids.baseline(S), secretion=SECRETION, variants=STANDARD
    )
    names = [*r.model.getGlobalParameterIds(), *r.model.getCompartmentIds()]
    values = {n: r[n] for n in names}
    assert _evaluate(arterial, values) == pytest.approx(cki, rel=1e-12)
    inverse = ss.kidney_arterial(s, "C_ki_0", secretion=SECRETION, variants=STANDARD)
    assert _evaluate(inverse, values) == pytest.approx(car, rel=1e-12)
    # the excretion rate of the urine is the production
    assert r[ids.urine_rate(S)] == pytest.approx(r["k_prod"], rel=1e-10)


@pytest.mark.parametrize("load", [1e-5, 0.009, 0.5, 50.0])
def test_the_root_solves_the_kidney_balance(load: float) -> None:
    """The root is exact at a low and a high load (no cancellation, either sign of b)."""
    s = _substance()
    formula = ss.kidney_plasma(
        s, arterial="Car", secretion=SECRETION, clearance="CL", variants=STANDARD
    )
    values = {
        "Qki": 1.2,
        "HCT": 0.45,
        "fu_end": 0.8,
        "GFR": 0.12,
        "f_reabsorbed_end": 0.1,
        "Car": load,
        "CL": 0.01,
        "V_sec": 0.002,
        "km_sec": 0.08,
    }
    c = _evaluate(formula, values)
    rest = 1.2 * (1 - 0.45) - 0.8 * 0.12 * 0.9
    residual = rest * (load - c) - 0.01 * c - 0.002 * c / (0.08 + c)
    assert c > 0
    # the terms of the balance are of the order rest*load: a few ulp of them
    assert abs(residual) <= 1e-13 * rest * load


def test_without_a_saturable_secretion_the_form_is_linear() -> None:
    """The elimination of a filtered substance alone: `C = E/G`; both inputs refused."""
    s = _substance()
    plain = ss.kidney_plasma(s, eliminated="k_prod", variants=STANDARD)
    assert "sqrt" not in plain
    values = {"k_prod": 0.009, "fu_end": 0.8, "GFR": 0.12, "f_reabsorbed_end": 0.1}
    assert _evaluate(plain, values) == pytest.approx(0.009 / (0.8 * 0.12 * 0.9))
    with pytest.raises(ValueError, match="not both"):
        ss.kidney_plasma(s, arterial="Car", eliminated="k_prod", variants=STANDARD)
    unfiltered = dataclasses.replace(s, excretion=Excretion())
    with pytest.raises(ValueError, match="eliminates nothing"):
        ss.kidney_plasma(unfiltered, eliminated="k_prod", variants=STANDARD)
    alone = ss.kidney_plasma(
        unfiltered, eliminated="k_prod", secretion=SECRETION, variants=STANDARD
    )
    assert "sqrt" not in alone
    assert ss.kidney_arterial(s, "C", variants=STANDARD) == "C"
    with pytest.raises(ValueError, match="needs no rate"):
        ss.kidney_arterial(
            s, "C", eliminated="E", secretion=SECRETION, variants=STANDARD
        )
