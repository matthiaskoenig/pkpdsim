"""The steady state of the equations of the library (`pkpdsim.pbpk.steady_state`).

An endogenous substance produced in the muscle (tissue) and supplied in the
liver (plasma), filtered and secreted by the kidney: with the formulas of the
module as the starts of the organs every species stays at its start, for the
reference individual and a changed one, with and without erythrocytes in
equilibrium; the arm vein of the formulas is the observer at the start. The
variants which have no closed steady state are refused, naming the reason.
"""

import dataclasses
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pytest

from pkpdsim.pbpk import (
    Binding,
    CustomReaction,
    Erythrocytes,
    Excretion,
    Formula,
    Kp,
    Mode,
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

#: the production of the muscle, the supply of the liver (mmole/min), the
#: secretion clearance of the kidney plasma (l/min)
PARAMETERS = [
    ModelParameter("k_prod", 0.0074, U.mmole_per_min),
    ModelParameter("k_syn", 0.002, U.mmole_per_min),
    ModelParameter("CL_sec", 0.05, U.l_per_min),
]

CHEMISTRY = [
    CustomReaction("PROD", Organ.MUSCLE, Part.TISSUE, f"-> Cti_{S}", "k_prod", []),
    CustomReaction("SYN", Organ.LIVER, Part.PLASMA, f"-> Cpl_{S}", "k_syn", []),
    CustomReaction(
        "SEC",
        Organ.KIDNEY,
        Part.PLASMA,
        f"Cpl_{S} -> Aurine_{S}",
        f"CL_sec * Cpl_{S}",
        [],
    ),
]


def _substance(erythrocytes: Erythrocytes | None = None) -> Substance:
    """A bound substance, filtered with a reabsorbed fraction, at its steady state."""
    return Substance(
        S,
        name="endogenous substance",
        mr=113.1,
        binding=Binding(fu=0.5, protein="albumin"),
        erythrocytes=erythrocytes,
        distribution=Kp(0.6, ftissue=6.0, overrides={Organ.MUSCLE: 0.8}),
        excretion=Excretion(urine=True, fraction_reabsorbed=0.1),
        # production and supply over the renal clearance of the arterial plasma
        baseline=Formula("(k_prod + k_syn) / CL_ar", U.mM),
    )


def _definition(
    s: Substance, variants: Variants | None = None, **fields: object
) -> PBPKDefinition:
    """The definition of the substance with the starts of the module."""
    starts = _starts(s, variants)
    return PBPKDefinition(
        sid="steady",
        substances=[s],
        chemistry=CHEMISTRY,
        variants=variants or Variants(),
        parameters=[
            *PARAMETERS,
            ModelParameter(
                "CL_ar",
                Formula(
                    ss.renal_clearance(s, clearance="CL_sec", variants=variants),
                    U.l_per_min,
                ),
                U.l_per_min,
            ),
        ],
        starts=starts,
        **fields,  # ty: ignore[invalid-argument-type]
    )


def _muscle_plasma(s: Substance) -> str:
    """The plasma leaving the muscle with its production."""
    return ss.organ_plasma(s, Organ.MUSCLE, rate="k_prod")


def _starts(s: Substance, variants: Variants | None) -> dict[str, Formula]:
    """The steady state of the muscle, the kidney and the liver."""
    muscle = _muscle_plasma(s)
    kidney = ss.kidney_plasma(s, clearance="CL_sec", variants=variants)
    liver = ss.liver_plasma(s, rate="k_syn", variants=variants)
    return {
        f"MU__Cpl_{S}": Formula(muscle, U.mM),
        f"MU__Cti_{S}": Formula(
            ss.organ_tissue(s, Organ.MUSCLE, muscle, rate="k_prod"), U.mM
        ),
        f"KI__Cpl_{S}": Formula(kidney, U.mM),
        f"KI__Cti_{S}": Formula(ss.organ_tissue(s, Organ.KIDNEY, kidney), U.mM),
        f"LI__Cpl_{S}": Formula(liver, U.mM),
        f"LI__Cti_{S}": Formula(ss.organ_tissue(s, Organ.LIVER, liver), U.mM),
    }


def _drift(r: "roadrunner.RoadRunner") -> dict[str, float]:
    """The largest relative change of every species but the urine over a day."""
    sids = [s for s in r.model.getFloatingSpeciesIds() if s != ids.urine(S)]
    start = np.array([r[f"[{sid}]"] for sid in sids])
    values = np.asarray(r.simulate(0, DAY, 97, selections=[f"[{sid}]" for sid in sids]))
    scale = np.where(start > 0, start, 1.0)
    return dict(zip(sids, np.max(np.abs(values - start), axis=0) / scale, strict=True))


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
def test_the_steady_state_of_the_organs_is_stationary(
    tmp_path: Path, changes: dict[str, float], erythrocytes: Erythrocytes | None
) -> None:
    """The formulas as starts: every species stays within 1e-6 of its start.

    The muscle with its production (plasma flow, the share of ftissue of its
    tissue, its own Kp), the kidney with the filtration less the reabsorbed
    fraction and a secretion, the liver with its supply past the shunt; the
    baseline is the arterial plasma of the steady state (`renal_clearance`).
    """
    s = _substance(erythrocytes)
    r = individual(build(_definition(s), tmp_path).flat_path, **changes)
    r.integrator.absolute_tolerance = 1e-12
    r.integrator.relative_tolerance = 1e-10
    moved = {sid: d for sid, d in _drift(r).items() if d > 1e-6}
    assert not moved, moved
    assert r[f"[MU__Cpl_{S}]"] > r[f"[Car_{S}]"] > r[f"[KI__Cpl_{S}]"]
    assert r[f"[Car_{S}]"] == pytest.approx(r[ids.baseline(S)], rel=1e-12)
    for organ in (Organ.MUSCLE, Organ.LIVER):
        volume = r[ids.flat(organ, "Vti")]
        assert _evaluate(r, ss.tissue_volume(organ)) == pytest.approx(volume, rel=1e-12)


def test_the_arm_vein_of_the_steady_state_is_the_observer(tmp_path: Path) -> None:
    """`arm_vein` is the arm vein observer at the steady state, for every preset.

    The preset `musther2015` mixes the arterial plasma and adipose besides
    skin and muscle.
    """
    s = _substance()
    for preset in ("pksim_default", "musther2015"):
        d = _definition(s, arm_vein=preset)
        r = individual(build(d, tmp_path / preset).flat_path, BW=90.0)
        formula = ss.arm_vein(preset, s, {Organ.MUSCLE: _muscle_plasma(s)})
        for symbol in sorted(
            {"w_arm_sk", "w_arm_mu", "w_arm_ad", "w_arm_ar"} & set(formula.split())
        ):
            assert r[symbol] > 0, symbol
        assert r[ids.arm(S)] == pytest.approx(_evaluate(r, formula), rel=1e-12)


def _evaluate(r: "roadrunner.RoadRunner", formula: str) -> float:
    """Evaluate a formula over the ids of a loaded model."""
    import libsbml

    ast = libsbml.parseL3Formula(formula.replace(" dimensionless", ""))
    names: set[str] = set()
    stack = [ast]
    while stack:
        node = stack.pop()
        if node.getType() == libsbml.AST_NAME:
            names.add(node.getName())
        stack.extend(node.getChild(k) for k in range(node.getNumChildren()))
    expression = libsbml.formulaToL3String(ast).replace("^", "**")
    return float(eval(expression, {}, {n: r[n] for n in names}))


def test_the_zonated_liver_has_the_inflow_in_every_zone(tmp_path: Path) -> None:
    """In the zonated liver (no source) every zone starts at `liver_plasma`.

    The zonated liver refuses custom reactions, so the supply moves to the
    muscle; the plasma and the tissue of every zone stay at the inflow mix.
    A source in the zonated liver is refused.
    """
    variants = Variants(organs={Organ.LIVER: "zonated"})
    s = _substance()
    liver = ss.liver_plasma(s, variants=variants)
    with pytest.raises(NotImplementedError, match="variant of the liver"):
        ss.liver_plasma(s, rate="k_syn", variants=variants)
    muscle = ss.organ_plasma(s, Organ.MUSCLE, rate="k_prod + k_syn")
    kidney = ss.kidney_plasma(s, clearance="CL_sec")
    d = PBPKDefinition(
        sid="zonated",
        substances=[s],
        chemistry=[
            CustomReaction(
                "PROD", Organ.MUSCLE, Part.TISSUE, f"-> Cti_{S}", "k_prod + k_syn", []
            ),
            CHEMISTRY[2],
        ],
        variants=variants,
        parameters=[
            *PARAMETERS,
            ModelParameter(
                "CL_ar",
                Formula(ss.renal_clearance(s, clearance="CL_sec"), U.l_per_min),
                U.l_per_min,
            ),
        ],
        starts={
            f"MU__Cpl_{S}": Formula(muscle, U.mM),
            f"MU__Cti_{S}": Formula(
                ss.organ_tissue(s, Organ.MUSCLE, muscle, rate="k_prod + k_syn"), U.mM
            ),
            f"KI__Cpl_{S}": Formula(kidney, U.mM),
            f"KI__Cti_{S}": Formula(ss.organ_tissue(s, Organ.KIDNEY, kidney), U.mM),
        },
    )
    r = individual(build(d, tmp_path).flat_path, f_renal_function=0.5)
    r.integrator.absolute_tolerance = 1e-12
    r.integrator.relative_tolerance = 1e-10
    moved = {sid: x for sid, x in _drift(r).items() if x > 1e-6}
    assert not moved, moved
    zones = [
        sid
        for sid in r.model.getFloatingSpeciesIds()
        if sid.startswith("LI__") and sid.endswith(f"_{S}")
    ]
    assert len(zones) > 2
    plasma = _evaluate(r, liver)
    for sid in zones:
        expected = plasma if "Cpl" in sid else r["Kp_end"] * plasma
        assert r[f"[{sid}]"] == pytest.approx(expected, rel=1e-12), sid
    tissue = ss.organ_tissue(s, Organ.LIVER, liver, variants=variants)
    assert _evaluate(r, tissue) == pytest.approx(r["Kp_end"] * plasma, rel=1e-12)


def test_the_nephron_kidney_is_refused() -> None:
    """The tubule of the nephron kidney starts empty: no steady state of the kidney."""
    variants = Variants(organs={Organ.KIDNEY: "nephron"})
    s = _substance()
    for f in (ss.kidney_plasma, ss.renal_clearance):
        with pytest.raises(NotImplementedError, match="nephron kidney"):
            f(s, variants=variants)
    with pytest.raises(NotImplementedError, match="variant of the organ 'kidney'"):
        ss.organ_tissue(s, Organ.KIDNEY, "1 mM", variants=variants)


@pytest.mark.parametrize("organ", [Organ.SMALL_INTESTINE, Organ.LARGE_INTESTINE])
def test_the_intestine_is_refused(organ: Organ) -> None:
    """The segmented intestine refuses custom reactions: its parts start at the equilibrium."""
    s = _substance()
    with pytest.raises(NotImplementedError, match="intestine refuses custom"):
        ss.organ_plasma(s, organ, rate="k_prod")
    with pytest.raises(NotImplementedError, match="intestine refuses custom"):
        ss.organ_tissue(s, organ, "1 mM")


def test_the_forms_without_a_closed_steady_state_are_refused() -> None:
    """Kinetic erythrocytes, the full form, an exchange of the organ and the wrong helper."""
    kinetic = _substance(Erythrocytes("kinetic", kbc=1.5, k_exchange=0.1))
    with pytest.raises(NotImplementedError, match="kinetic erythrocytes"):
        ss.organ_plasma(kinetic, Organ.MUSCLE)
    full = dataclasses.replace(
        _substance(), distribution=Kp(0.6, ftissue=6.0, mode=Mode.PERMEABILITY)
    )
    with pytest.raises(NotImplementedError, match="permeability"):
        ss.organ_tissue(full, Organ.MUSCLE, "1 mM")
    own = dataclasses.replace(
        _substance(), distribution=Kp(0.6, ftissue=6.0, ps={Organ.MUSCLE: 1.0})
    )
    with pytest.raises(NotImplementedError, match="value of the organ"):
        ss.organ_tissue(own, Organ.MUSCLE, "1 mM", rate="k_prod")
    with pytest.raises(ValueError, match="kidney_plasma"):
        ss.organ_plasma(_substance(), Organ.KIDNEY)
    with pytest.raises(ValueError, match="liver_plasma"):
        ss.organ_plasma(_substance(), Organ.LIVER)
    unfiltered = dataclasses.replace(_substance(), excretion=Excretion())
    with pytest.raises(ValueError, match="not filtered"):
        ss.renal_filtration(unfiltered)
    # without filtration the kidney plasma is the arterial plasma with the
    # secretion alone, the renal clearance the secretion of the plasma flow
    assert "GFR" not in ss.kidney_plasma(unfiltered, clearance="CL_sec")
    assert "GFR" not in ss.renal_clearance(unfiltered, clearance="CL_sec")
