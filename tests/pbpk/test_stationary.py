"""The generated steady state (`pkpdsim.pbpk.stationary`, `CustomReaction.stationary`).

An endogenous substance whose custom reactions declare their rate at the
steady state starts stationary without a start of the model: every species
stays within 1e-12 (relative) of its start over a day, for the standard
organs (a production of the muscle tissue, a supply of the liver plasma, a
clearance of the liver tissue, a source of a portal organ, a saturable
secretion of the kidney), the nephron kidney, the zonated liver, a slow
pool, a lung which moves the substance, a threshold reabsorption, another
individual and an activity of an organ set before the simulation. The
baseline of a substance without one is solved by the library; a substance
with a baseline is pinned at it. What has no closed form is refused with
the reason.
"""

import dataclasses
import warnings
from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING

import libsbml
import numpy as np
import pytest
from sbmlutils.factory import Parameter

from pkpdsim.pbpk import (
    MM,
    Binding,
    CustomReaction,
    Erythrocytes,
    Excretion,
    Formula,
    Kp,
    Metabolism,
    Mode,
    ModelParameter,
    Organ,
    Part,
    PBPKDefinition,
    Route,
    Saturable,
    SlowPool,
    Substance,
    ThresholdReabsorption,
    U,
    Variants,
    build,
    ids,
    steady_state,
)
from pkpdsim.pbpk.stationary import (
    BASELINE_SIZE_LIMIT,
    SteadyStateWarning,
    expanded_size,
    linear_clearance,
    steady_plan,
)

from .conftest import (
    assert_reference_sbml,
    endogenous_reference_definition,
    individual,
)

if TYPE_CHECKING:
    import roadrunner

S = "end"
DAY = 24 * 60.0
STANDARD = Variants()

PARAMETERS = [
    ModelParameter("k_prod", 0.0074, U.mmole_per_min),
    ModelParameter("k_syn", 0.002, U.mmole_per_min),
    ModelParameter("k_pa", 0.001, U.mmole_per_min),
    ModelParameter("k_cl", 0.01, U.per_min),
    ModelParameter("k_pacl", 0.002, U.l_per_min),
    ModelParameter("Vm", 0.01, U.mmole_per_min),
    ModelParameter("Km", 0.05, U.mM),
]


def _rate(expression: str) -> Formula:
    """A rate at the steady state."""
    return Formula(expression, U.mmole_per_min)


PROD = CustomReaction(
    "PROD",
    Organ.MUSCLE,
    Part.TISSUE,
    f"-> Cti_{S}",
    "k_prod",
    stationary=_rate("k_prod"),
)
SYN = CustomReaction(
    "SYN", Organ.LIVER, Part.PLASMA, f"-> Cpl_{S}", "k_syn", stationary=_rate("k_syn")
)
LICL = CustomReaction(
    "LICL",
    Organ.LIVER,
    Part.TISSUE,
    f"Cti_{S} ->",
    f"f_act * k_cl * Vti * Cti_{S}",
    parameters=[Parameter("f_act", 1.0, U.dimensionless, constant=True)],
    stationary="linear",
)
PASRC = CustomReaction(
    "PASRC",
    Organ.PANCREAS,
    Part.PLASMA,
    f"-> Cpl_{S}",
    "k_pa",
    stationary=_rate("k_pa"),
)
SEC = CustomReaction(
    "SEC",
    Organ.KIDNEY,
    Part.PLASMA,
    f"Cpl_{S} -> Aurine_{S}",
    f"Vm * Cpl_{S} / (Km + Cpl_{S})",
    stationary=Saturable("Vm", "Km"),
)
PACL = CustomReaction(
    "PACL",
    Organ.PANCREAS,
    Part.PLASMA,
    f"Cpl_{S} ->",
    f"k_pacl * Cpl_{S}",
    stationary="linear",
)
CHEMISTRY = [PROD, SYN, LICL, PASRC, PACL]
"""The linear steady state: the baseline is the net source over the clearance
(a portal organ with a source and a clearance lowers the inflow of the liver)."""


def _substance(**fields: object) -> Substance:
    """A bound substance, filtered with a reabsorbed fraction, without baseline."""
    s = Substance(
        S,
        name="endogenous substance",
        mr=113.1,
        binding=Binding(fu=0.5, protein="albumin"),
        distribution=Kp(0.6, ftissue=6.0, overrides={Organ.MUSCLE: 0.8}),
        excretion=Excretion(urine=True, fraction_reabsorbed=0.1),
    )
    return dataclasses.replace(s, **fields)


def _definition(
    chemistry: Sequence[CustomReaction] = CHEMISTRY,
    substances: Sequence[Substance] | None = None,
    **fields: object,
) -> PBPKDefinition:
    """The definition of the generated steady state."""
    return PBPKDefinition(
        sid="stationary",
        substances=list(substances or [_substance()]),
        chemistry=list(chemistry),
        parameters=PARAMETERS + list(fields.pop("extra", [])),  # ty: ignore[invalid-argument-type]
        **fields,  # ty: ignore[invalid-argument-type]
    )


def _load(d: PBPKDefinition, path: Path, **changes: float) -> "roadrunner.RoadRunner":
    """Build and load a definition with changes before the simulation, tight tolerances."""
    r = individual(build(d, path).flat_path, **changes)
    r.integrator.absolute_tolerance = 1e-16
    r.integrator.relative_tolerance = 1e-13
    return r


def _drift(r: "roadrunner.RoadRunner") -> dict[str, float]:
    """The largest relative change of every species but the sinks and of the arm vein over a day."""
    sids = [
        sid
        for sid in r.model.getFloatingSpeciesIds()
        if not sid.startswith(("Aurine_", "Afeces_", "Abile_"))
    ]
    selections = [f"[{sid}]" for sid in sids]
    arm = [x for x in r.model.getGlobalParameterIds() if x.startswith("Carm_")]
    selections += arm
    start = np.array([r[x] for x in selections])
    # densely in the first hour, where a wrong start of a fast part relaxes,
    # then over the day
    early = np.asarray(r.simulate(0, 60, 61, selections=selections))
    late = np.asarray(r.simulate(60, DAY, 24, selections=selections))
    values = np.vstack([early, late])
    scale = np.where(np.abs(start) > 0, np.abs(start), 1.0)
    return dict(
        zip(selections, np.max(np.abs(values - start), axis=0) / scale, strict=True)
    )


def _assert_stationary(r: "roadrunner.RoadRunner") -> None:
    """Every species stays within 1e-12 of its start over a day."""
    moved = {sid: d for sid, d in _drift(r).items() if d > 1e-12}
    assert not moved, moved


# --- the standard organs ---------------------------------------------------------


@pytest.mark.parametrize(
    "changes",
    [
        {},
        {"BW": 110.0, "HEIGHT": 190.0},
        {"f_renal_function": 0.3},
        {"f_shunt": 0.3, "f_functional": 0.6},
    ],
    ids=["reference", "body", "renal function", "liver"],
)
@pytest.mark.parametrize(
    "erythrocytes",
    [None, Erythrocytes("equilibrium", bp=0.8)],
    ids=["plasma", "erythrocytes"],
)
def test_the_solved_steady_state_is_stationary(
    tmp_path: Path, changes: dict[str, float], erythrocytes: Erythrocytes | None
) -> None:
    """A substance without baseline starts at the steady state which the library solves.

    The production of the muscle tissue, the supply of the liver plasma, the
    clearance of the liver tissue, the source and the clearance of the
    pancreas (the portal vein above the arterial plasma, the portal
    eliminations lower the inflow of the liver), the filtration: every
    species and the arm vein within 1e-12 of their start over a day, for
    another individual, an impaired renal function and a portosystemic shunt
    with a reduced functional liver; the generated baseline is the arterial
    plasma. The saturable root is `test_a_saturable_secretion_starts_at_the_root`.
    """
    d = _definition(substances=[_substance(erythrocytes=erythrocytes)])
    r = _load(d, tmp_path, **changes)
    _assert_stationary(r)
    assert r[f"[Car_{S}]"] == pytest.approx(r[ids.baseline(S)], rel=1e-14)
    # the production of the muscle, the source of the pancreas: their plasma
    # above the arterial plasma
    assert r[f"[MU__Cpl_{S}]"] > r[f"[Car_{S}]"]
    assert r[f"[Cpo_{S}]"] != pytest.approx(r[f"[Car_{S}]"], rel=1e-6)


@pytest.mark.parametrize(
    "changes",
    [{}, {"BW": 60.0, "HEIGHT": 160.0, "f_renal_function": 0.3}],
    ids=["reference", "body"],
)
def test_a_saturable_secretion_starts_at_the_root(
    tmp_path: Path, changes: dict[str, float]
) -> None:
    """The saturable secretion of the kidney: the kidney at the root of the quadratic.

    The production of the muscle eliminated by the filtration, the liver
    clearance and the saturable secretion of the kidney plasma; the kidney
    starts at the root which the top model solves with the baseline
    (`Css_ki_<s>`).
    """
    d = _definition([PROD, LICL, SEC])
    r = _load(d, tmp_path, **changes)
    _assert_stationary(r)
    assert r[ids.steady_plasma(Organ.KIDNEY, S)] == pytest.approx(
        r[f"[KI__Cpl_{S}]"], rel=1e-14
    )
    assert r[f"[KI__Cpl_{S}]"] < r[f"[Car_{S}]"]


def test_the_baseline_follows_the_mass_balance(tmp_path: Path) -> None:
    """Without the saturable secretion the baseline is the sources over the clearance.

    The production and the supply are eliminated by the filtration
    (`fu*GFR*(1 - f_reabsorbed)`) and the liver; without the liver clearance
    the filtration alone: `Car = (k_prod + k_syn)/(fu*GFR*(1 - f_reabsorbed))`.
    """
    d = _definition([PROD, SYN])
    r = _load(d, tmp_path)
    _assert_stationary(r)
    filtered = r[ids.fu(S)] * r[ids.GFR] * (1 - r[ids.fraction_reabsorbed(S)])
    assert r[ids.baseline(S)] == pytest.approx((0.0074 + 0.002) / filtered, rel=1e-12)


def test_a_lung_which_clears_a_solved_substance(tmp_path: Path) -> None:
    """The lung plasma is the arterial plasma: it eliminates `CLss*Car - Rss` in the solve.

    The production of the muscle, a source and a linear clearance of the lung
    plasma; the venous pool starts at the inflow of the lung which gives the
    arterial plasma.
    """
    lung = [
        CustomReaction(
            "LUCL",
            Organ.LUNG,
            Part.PLASMA,
            f"Cpl_{S} ->",
            f"k_lu * Cpl_{S}",
            stationary="linear",
        ),
        CustomReaction(
            "LUSRC",
            Organ.LUNG,
            Part.TISSUE,
            f"-> Cti_{S}",
            "k_pa",
            stationary=_rate("k_pa"),
        ),
    ]
    d = _definition([PROD, *lung], extra=[ModelParameter("k_lu", 0.5, U.l_per_min)])
    r = _load(d, tmp_path, BW=80.0)
    _assert_stationary(r)
    assert r[f"[LU__Cpl_{S}]"] == pytest.approx(r[f"[Car_{S}]"], rel=1e-14)
    assert r[f"[Cve_{S}]"] > r[f"[Car_{S}]"]


def test_an_activity_set_before_the_simulation_moves_the_steady_state(
    tmp_path: Path,
) -> None:
    """A parameter of a custom reaction set before the simulation is read by the starts (D16).

    The activity `f_act` of the clearance of the liver tissue (a parameter
    of the organ, `LI__f_act`) at 2.5: the solved baseline falls, every
    species starts at the new steady state.
    """
    d = _definition([PROD, SYN, LICL])
    path = build(d, tmp_path).flat_path
    reference = individual(path)
    r = individual(path, LI__f_act=2.5)
    r.integrator.absolute_tolerance = 1e-16
    r.integrator.relative_tolerance = 1e-13
    _assert_stationary(r)
    assert r[ids.baseline(S)] < 0.99 * reference[ids.baseline(S)]
    assert r[f"[LI__Cti_{S}]"] < reference[f"[LI__Cti_{S}]"]


# --- a pinned substance: the baseline given ------------------------------------------


def _pinned(
    variants: Variants = STANDARD, lung: str = "linear", supply: str = ""
) -> PBPKDefinition:
    """A substance with a baseline and a threshold reabsorption, supplied to balance.

    The muscle takes up `k_up` from its plasma (saturable, its capacity from
    the steady state of the muscle plasma, `steady_state.organ_plasma`), the
    lung clears its plasma (`linear`: `k_lu`; `saturable`: `V_lu`, `Km_lu`;
    the lung plasma is the arterial plasma), and the supply of the heart
    balances the uptake, the lung clearance and the urinary excretion at the
    baseline: the balance of the model, with the formulas of
    `pkpdsim.pbpk.steady_state` (the supply sits in the heart, the zonated
    liver refuses custom reactions). `supply` is added to the supply (an
    imbalance).
    """
    s = _substance(
        baseline=0.03,
        excretion=Excretion(
            urine=True, reabsorption=ThresholdReabsorption(threshold=0.1, n=4.0)
        ),
    )
    excreted = f"{steady_state.renal_filtration(s, variants=STANDARD)} * baseline_{S}"
    muscle = steady_state.organ_plasma(s, Organ.MUSCLE, rate="-k_up", variants=STANDARD)
    cleared = {
        "linear": f" + k_lu * baseline_{S}",
        "saturable": f" + V_lu * baseline_{S} / (Km_lu + baseline_{S})",
    }.get(lung, "")
    rate = f"k_up + {excreted}{cleared}" + (f" + {supply}" if supply else "")
    chemistry = [
        CustomReaction(
            "UP",
            Organ.MUSCLE,
            Part.PLASMA,
            f"Cpl_{S} ->",
            f"V_up * Cpl_{S} / (Km_up + Cpl_{S})",
            stationary=_rate("k_up"),
        ),
        CustomReaction(
            "SUP",
            Organ.HEART,
            Part.PLASMA,
            f"-> Cpl_{S}",
            "k_sup",
            stationary=_rate("k_sup"),
        ),
    ]
    if lung == "saturable":
        chemistry.append(
            CustomReaction(
                "LUSAT",
                Organ.LUNG,
                Part.PLASMA,
                f"Cpl_{S} ->",
                f"V_lu * Cpl_{S} / (Km_lu + Cpl_{S})",
                stationary=Saturable("V_lu", "Km_lu"),
            )
        )
    elif lung == "linear":
        chemistry.append(
            CustomReaction(
                "LUCL",
                Organ.LUNG,
                Part.PLASMA,
                f"Cpl_{S} ->",
                f"k_lu * Cpl_{S}",
                stationary="linear",
            )
        )
    parameters = [
        ModelParameter("k_up", 0.004, U.mmole_per_min),
        ModelParameter("Km_up", 0.01, U.mM),
        ModelParameter(
            "V_up",
            Formula(f"k_up * (Km_up + {muscle}) / ({muscle})", U.mmole_per_min),
            U.mmole_per_min,
        ),
        ModelParameter("k_lu", 0.05, U.l_per_min),
        ModelParameter("V_lu", 0.002, U.mmole_per_min),
        ModelParameter("Km_lu", 0.02, U.mM),
        ModelParameter("k_sup", Formula(rate, U.mmole_per_min), U.mmole_per_min),
    ]
    return PBPKDefinition(
        sid="pinned",
        substances=[s],
        chemistry=chemistry,
        parameters=parameters,
        variants=variants,
    )


@pytest.mark.parametrize(
    "changes",
    [{}, {"BW": 60.0, "HEIGHT": 165.0}, {"f_renal_function": 0.4}],
    ids=["reference", "body", "renal function"],
)
@pytest.mark.parametrize("lung", ["linear", "saturable"])
def test_a_pinned_steady_state_is_stationary(
    tmp_path: Path, changes: dict[str, float], lung: str
) -> None:
    """A substance with a baseline starts at the steady state of the baseline.

    A saturable uptake of the muscle declared at its rate, a supply of the
    heart, a linear or saturable clearance of the lung (the venous pool at the inflow of
    the lung which gives the arterial plasma) and a threshold reabsorption
    (the kidney plasma at the arterial plasma).
    """
    r = _load(_pinned(lung=lung), tmp_path, **changes)
    _assert_stationary(r)
    assert r[f"[Car_{S}]"] == pytest.approx(0.03, rel=1e-14)
    assert r[f"[Cve_{S}]"] > r[f"[Car_{S}]"] > r[f"[MU__Cpl_{S}]"]
    assert r[f"[KI__Cpl_{S}]"] == pytest.approx(0.03, rel=1e-14)


# --- the variants ------------------------------------------------------------------


@pytest.mark.parametrize("changes", [{}, {"BW": 95.0, "f_renal_function": 0.5}])
def test_the_nephron_kidney_starts_with_the_filtrate_in_the_tubule(
    tmp_path: Path, changes: dict[str, float]
) -> None:
    """The tubule starts at the flow of the filtrate, `fu*Car/f_tub_<seg>`.

    The nephron kidney refuses custom reactions: the substance is produced
    in the muscle and eliminated by the filtration and the liver; a pinned
    substance with a threshold reabsorption in the nephron kidney as well.
    """
    d = _definition(
        [PROD, SYN, LICL], variants=Variants(organs={Organ.KIDNEY: "nephron"})
    )
    r = _load(d, tmp_path / "solved", **changes)
    _assert_stationary(r)
    car = r[f"[Car_{S}]"]
    for seg in ids.TUBULE_SEGMENTS:
        expected = r[ids.fu(S)] * car / r[f"KI__{ids.tubule_flow_fraction(seg)}"]
        assert r[f"[KI__{ids.tubule(S, seg)}]"] == pytest.approx(expected, rel=1e-12)
    pinned = _pinned(Variants(organs={Organ.KIDNEY: "nephron"}))
    _assert_stationary(_load(pinned, tmp_path / "pinned", **changes))


def test_the_zonated_liver_starts_at_the_mix_of_its_inflows(tmp_path: Path) -> None:
    """Every zone of the zonated liver at the mix of the hepatic artery and the portal vein.

    The pancreas raises the portal vein above the arterial plasma; the
    zonated liver refuses custom reactions, so the liver neither produces nor
    clears. A pinned substance in the zonated liver as well.
    """
    variants = Variants(organs={Organ.LIVER: "zonated"})
    d = _definition([PROD, PASRC, SEC], variants=variants)
    r = _load(d, tmp_path / "solved", f_shunt=0.2)
    _assert_stationary(r)
    mix = (r["Qha"] * r[f"[Car_{S}]"] + r["Qpo"] * r[f"[Cpo_{S}]"]) / r["Qh"]
    zones = [
        sid
        for sid in r.model.getFloatingSpeciesIds()
        if sid.startswith("LI__C") and sid.endswith(f"_{S}")
    ]
    assert len(zones) > 2
    for sid in zones:
        expected = mix if "Cpl" in sid else r[f"Kp_{S}"] * mix
        assert r[f"[{sid}]"] == pytest.approx(expected, rel=1e-12), sid
    _assert_stationary(_load(_pinned(variants), tmp_path / "pinned"))


def test_a_slow_pool_starts_at_its_tissue(tmp_path: Path) -> None:
    """A slow pool starts at `Kpslow` times the start of its tissue, not of the baseline."""
    d = _definition(
        slow_pools=[SlowPool(S, frozenset({Organ.MUSCLE}), k=0.001, kp=3.0)]
    )
    r = _load(d, tmp_path)
    _assert_stationary(r)
    assert r[f"[MU__Cslow_{S}]"] == pytest.approx(3.0 * r[f"[MU__Cti_{S}]"], rel=1e-14)
    assert r[f"[MU__Cti_{S}]"] > 0.8 * r[f"[MU__Cpl_{S}]"]


def test_parts_without_a_process_start_at_their_inflow(tmp_path: Path) -> None:
    """The full form and the flow-limited organ start at their inflow, a plasma-only substance.

    A substance which enters only some organs (the full form of the
    permeability mode, the flow-limited organs of the perfusion mode) with
    its processes in organs without tissue.
    """
    organs = frozenset({Organ.ADIPOSE, Organ.SKIN})
    for mode in (Mode.PERMEABILITY, Mode.PERFUSION):
        s = _substance(distribution=Kp(0.6, ftissue=6.0, organs=organs, mode=mode))
        d = _definition(
            [dataclasses.replace(PROD, part=Part.PLASMA, equation=f"-> Cpl_{S}"), SYN],
            substances=[s],
        )
        _assert_stationary(_load(d, tmp_path / mode.value))
    plain = _substance(distribution=None)
    d = _definition(
        [dataclasses.replace(PROD, part=Part.PLASMA, equation=f"-> Cpl_{S}"), SYN],
        substances=[plain],
    )
    _assert_stationary(_load(d, tmp_path / "plasma"))


def test_a_product_starts_at_the_steady_state_of_its_substrate(tmp_path: Path) -> None:
    """A linear conversion in the liver: the product's source reads the substrate's start.

    Two molecules of the product per molecule of the substrate: the source of
    the product is twice the rate at the start of the substrate.
    """
    product = _substance(sid="end_m", name="product")
    convert = CustomReaction(
        "CONV",
        Organ.LIVER,
        Part.PLASMA,
        f"Cpl_{S} -> 2 Cpl_end_m",
        f"k_conv * Cpl_{S}",
        stationary="linear",
    )
    d = _definition(
        [PROD, convert],
        substances=[_substance(), product],
        extra=[ModelParameter("k_conv", 0.2, U.l_per_min)],
    )
    r = _load(d, tmp_path)
    _assert_stationary(r)
    assert r[ids.baseline("end_m")] > 0


def test_an_explicit_start_overrides_the_generated_one(tmp_path: Path) -> None:
    """A start of the definition replaces the generated start of its species."""
    d = _definition(starts={f"MU__Cpl_{S}": 0.123})
    r = individual(build(d, tmp_path).flat_path)
    assert r[f"[MU__Cpl_{S}]"] == pytest.approx(0.123)


def test_an_oral_substance_keeps_its_lumen_empty(tmp_path: Path) -> None:
    """The segmented intestine starts at the arterial plasma, the lumen empty."""
    d = _definition(routes=[Route.oral(S)])
    r = _load(d, tmp_path)
    _assert_stationary(r)
    assert r[f"[GU__Cpl_{S}]"] == pytest.approx(r[f"[Car_{S}]"], rel=1e-14)


def test_the_flow_limited_zonated_liver_passes_its_inflow(tmp_path: Path) -> None:
    """The zones of a flow-limited zonated liver start at the inflow of the zone upstream.

    A substance of the perfusion mode which enters the liver and the adipose
    tissue, produced in the plasma of the muscle and supplied by the pancreas
    (the portal vein above the arterial plasma): the organ species of every
    zone at `K_or` times the plasma leaving the zone upstream.
    """
    s = _substance(
        distribution=Kp(
            0.6,
            ftissue=6.0,
            organs=frozenset({Organ.LIVER, Organ.ADIPOSE}),
            mode=Mode.PERFUSION,
        )
    )
    production = dataclasses.replace(PROD, part=Part.PLASMA, equation=f"-> Cpl_{S}")
    d = _definition(
        [production, PASRC],
        substances=[s],
        variants=Variants(organs={Organ.LIVER: "zonated"}),
    )
    r = _load(d, tmp_path)
    _assert_stationary(r)
    mix = (r["Qha"] * r[f"[Car_{S}]"] + r["Qpo"] * r[f"[Cpo_{S}]"]) / r["Qh"]
    zones = sorted(
        x
        for x in r.model.getFloatingSpeciesIds()
        if x.startswith("LI__Cor_z") and x.endswith(f"_{S}")
    )
    assert len(zones) > 2
    for sid in zones:
        assert r[f"[{sid}]"] == pytest.approx(r[f"LI__K_or_{S}"] * mix, rel=1e-12)


def test_the_single_lumen_starts_at_the_arterial_plasma(tmp_path: Path) -> None:
    """The intestine of a single lumen: the wall at the arterial plasma, the lumen empty."""
    d = _definition(
        routes=[Route.oral(S)],
        variants=Variants(organs={Organ.SMALL_INTESTINE: "single_lumen"}),
    )
    r = _load(d, tmp_path)
    _assert_stationary(r)
    assert r[f"[GU__Cpl_{S}]"] == pytest.approx(r[f"[Car_{S}]"], rel=1e-14)


# --- the checks of the build --------------------------------------------------------


def test_a_large_solved_baseline_warns(tmp_path: Path) -> None:
    """A solved saturable process of a body which also eliminates elsewhere: the build warns.

    The root repeats the clearance and the source of the rest of the body;
    above `BASELINE_SIZE_LIMIT` nodes every start compiles a slow copy.
    Without the saturable secretion the baseline is small and the build is
    quiet.
    """
    with pytest.warns(SteadyStateWarning, match="expands to"):
        build(_definition([*CHEMISTRY, SEC]), tmp_path / "large")
    with warnings.catch_warnings():
        warnings.simplefilter("error", SteadyStateWarning)
        result = build(_definition(), tmp_path / "small")
    model = libsbml.readSBMLFromFile(str(result.flat_path)).getModel()
    assert expanded_size(model, ids.baseline(S)) < BASELINE_SIZE_LIMIT


def test_a_pinned_imbalance_warns(tmp_path: Path) -> None:
    """A supply which does not balance the substance at its baseline: the build warns with the imbalance."""
    with pytest.warns(SteadyStateWarning, match="do not balance"):
        build(_pinned(supply="0.0001 mmole_per_min"), tmp_path / "off")
    with warnings.catch_warnings():
        warnings.simplefilter("error", SteadyStateWarning)
        build(_pinned(), tmp_path / "balanced")


# --- the refusals --------------------------------------------------------------------


def _refused(d: PBPKDefinition, match: str) -> None:
    """The plan refuses the definition with a message."""
    with pytest.raises(ValueError, match=match):
        steady_plan(d)


def test_an_undeclared_reaction_is_refused() -> None:
    """A custom reaction which moves a generated substance needs a declaration."""
    undeclared = dataclasses.replace(SYN, stationary=None)
    _refused(_definition([PROD, undeclared]), "custom reaction 'SYN' moves")


def test_a_nonlinear_rate_law_is_not_linear() -> None:
    """`linear` needs the substrate as one factor of the rate law."""
    for formula in (
        f"k_cl * Cpl_{S}^2",
        f"k_cl * Cpl_{S} / (Km + Cpl_{S})",
        f"k_cl * (Cpl_{S} + Km)",
        f"exp(Cpl_{S})",
    ):
        r = CustomReaction(
            "CL", Organ.LIVER, Part.PLASMA, f"Cpl_{S} ->", formula, stationary="linear"
        )
        _refused(_definition([PROD, r]), "not linear")


def test_a_library_process_is_refused() -> None:
    """A metabolism of the library moves the substance: no declared steady state."""
    m = Metabolism("MET", Organ.LIVER, S, {}, MM(vmax=0.01, km=0.1))
    d = PBPKDefinition(
        sid="lib", substances=[_substance()], chemistry=[PROD, m], parameters=PARAMETERS
    )
    _refused(d, "metabolism 'MET' of the library")


@pytest.mark.parametrize(
    ("fields", "chemistry", "match"),
    [
        (
            {"erythrocytes": Erythrocytes("kinetic", kbc=1.5, k_exchange=0.1)},
            [PROD, SEC],
            "kinetic erythrocytes",
        ),
        (
            {
                "excretion": Excretion(
                    urine=True, reabsorption=ThresholdReabsorption(0.1, 4.0)
                )
            },
            [PROD],
            "threshold reabsorption makes",
        ),
        (
            {
                "baseline": 0.01,
                "excretion": Excretion(
                    urine=True, reabsorption=ThresholdReabsorption(0.1, 4.0)
                ),
            },
            [PROD, SEC],
            "threshold reabsorption with a declared process",
        ),
        ({"excretion": Excretion()}, [PROD], "nothing eliminates"),
        (
            {"distribution": Kp(0.6, ftissue=6.0, mode=Mode.PERMEABILITY)},
            [PROD],
            "full form",
        ),
    ],
    ids=["kinetic", "threshold solved", "threshold kidney", "no elimination", "full"],
)
def test_a_substance_without_a_closed_form_is_refused(
    fields: dict[str, object], chemistry: list[CustomReaction], match: str
) -> None:
    """The properties of a substance whose steady state has no closed form."""
    _refused(_definition(chemistry, substances=[_substance(**fields)]), match)


def test_the_saturable_processes_of_a_solved_substance_are_refused() -> None:
    """Two saturable processes, or one of the liver: give the substance a baseline."""
    liver = dataclasses.replace(
        SEC, sid="SEC_LI", organ=Organ.LIVER, equation=f"Cpl_{S} ->"
    )
    _refused(_definition([PROD, SEC, liver]), "saturable processes in")
    _refused(_definition([PROD, liver]), "systemic organ")
    pinned = _definition([PROD, SEC, liver], substances=[_substance(baseline=0.01)])
    assert steady_plan(pinned).pinned == {S}


def test_two_saturable_rates_of_an_organ_are_refused() -> None:
    """The plasma of an organ has a closed steady state with one saturable rate."""
    second = dataclasses.replace(SEC, sid="SEC2")
    pinned = [_substance(baseline=0.01)]
    _refused(_definition([PROD, SEC, second], substances=pinned), "second saturable")


def test_a_declaration_reading_its_own_species_is_refused() -> None:
    """A stationary formula is the rate at the steady state: it reads no moved species."""
    own = dataclasses.replace(SYN, stationary=_rate(f"k_cl * Cpl_{S}"))
    _refused(_definition([PROD, own]), "reads species of")


@pytest.mark.parametrize(
    ("reaction", "fields", "match"),
    [
        (
            CustomReaction(
                "R",
                Organ.ADIPOSE,
                Part.ORGAN,
                f"-> Cor_{S}",
                "k_pa",
                stationary=_rate("k_pa"),
            ),
            {"distribution": Kp(0.6, ftissue=6.0, mode=Mode.PERFUSION)},
            "of the organ: a declared",
        ),
        (
            CustomReaction(
                "R",
                Organ.MUSCLE,
                Part.BLOOD_CELLS,
                f"-> Cbc_{S}",
                "k_pa",
                stationary=_rate("k_pa"),
            ),
            {"erythrocytes": Erythrocytes("kinetic", kbc=1.5, k_exchange=0.1)},
            "blood cells",
        ),
        (
            CustomReaction(
                "R",
                Organ.MUSCLE,
                Part.TISSUE,
                f"-> Cslow_{S}",
                "k_pa",
                stationary=_rate("k_pa"),
            ),
            {},
            "slow pool",
        ),
        (
            CustomReaction(
                "R",
                Organ.MUSCLE,
                Part.PLASMA,
                f"Car_{S} ->",
                f"k_cl * Car_{S}",
                stationary="linear",
            ),
            {},
            "blood pool",
        ),
        (
            CustomReaction(
                "R",
                Organ.LIVER,
                Part.PLASMA,
                f"Cpl_{S} + Cpl_b ->",
                f"k_cl * Vpl * Cpl_{S} * Cpl_b",
                stationary="linear",
            ),
            {},
            "needs one substrate",
        ),
        (
            CustomReaction(
                "R",
                Organ.LIVER,
                Part.TISSUE,
                f"Cti_{S} ->",
                f"Vm * Cti_{S} / (Km + Cti_{S})",
                stationary=Saturable("Vm", "Km"),
            ),
            {},
            "a rate of the plasma",
        ),
    ],
    ids=["flow-limited", "blood cells", "slow pool", "blood pool", "two", "tissue"],
)
def test_a_declared_process_without_a_closed_form_is_refused(
    reaction: CustomReaction, fields: dict[str, object], match: str
) -> None:
    """A declared process of a part or of substrates without a closed steady state."""
    substances = [_substance(**fields), _substance(sid="b")]
    _refused(_definition([PROD, reaction], substances=substances), match)


def test_a_transfer_between_parts_is_refused() -> None:
    """A linear transfer from the plasma to the tissue of the substance."""
    transfer = CustomReaction(
        "TR",
        Organ.MUSCLE,
        Part.PLASMA,
        f"Cpl_{S} -> Cti_{S}",
        f"k_cl * Vti * Cpl_{S}",
        stationary="linear",
    )
    _refused(_definition([PROD, transfer]), "transfer between two parts")


def test_a_cycle_of_products_is_refused() -> None:
    """Two substances converted into each other read each other's steady state."""
    there = CustomReaction(
        "AB",
        Organ.LIVER,
        Part.PLASMA,
        f"Cpl_{S} -> Cpl_b",
        f"k_cl * Vpl * Cpl_{S}",
        stationary="linear",
    )
    back = CustomReaction(
        "BA",
        Organ.MUSCLE,
        Part.PLASMA,
        f"Cpl_b -> Cpl_{S}",
        "k_cl * Vpl * Cpl_b",
        stationary="linear",
    )
    d = _definition([PROD, there, back], substances=[_substance(), _substance(sid="b")])
    _refused(d, "cycle of products")


def test_a_wrong_declaration_is_refused() -> None:
    """A stationary formula in mmole/min, a known kind of declaration."""
    wrong = dataclasses.replace(SYN, stationary=Formula("k_syn", U.mM))
    with pytest.raises(ValueError, match="stationary must be in"):
        _definition([PROD, wrong]).validate()
    other = dataclasses.replace(SYN, stationary="quadratic")
    with pytest.raises(ValueError, match="stationary is a str"):
        _definition([PROD, other]).validate()


# --- the linear rate law ----------------------------------------------------------


@pytest.mark.parametrize(
    ("formula", "clearance"),
    [
        ("k * C", "k"),
        ("C * k", "k"),
        ("k * V * C", "k * V"),
        ("k * C / V", "k / V"),
        ("C / V", "1 dimensionless / V"),
        ("C", None),
        ("-(k * C)", None),
        ("k * (C + 1 mM)", None),
        ("k * C * C", None),
        ("k / C", None),
        ("k * x", None),
    ],
)
def test_the_clearance_of_a_linear_rate_law(
    formula: str, clearance: str | None
) -> None:
    """The rate law without its substrate, `None` for a rate law not linear in it.

    A negated rate law is no clearance (a production of the substrate) and is
    refused, as the species alone (no clearance as a factor).
    """
    result = linear_clearance(formula, "C")
    if clearance is None:
        assert result is None
        return
    assert result is not None
    expected = libsbml.formulaToL3String(libsbml.parseL3Formula(clearance))
    assert result.replace(" ", "") == expected.replace(" ", "")


def test_the_endogenous_starts_without_a_declaration_are_those_of_87a5e79(
    tmp_path: Path,
) -> None:
    """Without a declaration the SBML of endogenous substances is that of 87a5e79.

    The generated steady state refactored the starts at the equilibrium with
    the baseline (the `Start` records, the zones, the inputs of the custom
    reactions): the merged tissue with a slow pool, an undeclared custom
    reaction, the mucosa and the lumen, the full form with kinetic
    erythrocytes, the flow-limited organ, the zonated liver, the nephron
    kidney and a start of the definition keep their SBML
    (`endogenous_reference_definition`).
    """
    assert_reference_sbml(build(endogenous_reference_definition(), tmp_path))


def test_without_a_declaration_there_is_no_plan() -> None:
    """A definition without a declared steady state has an empty plan: nothing changes."""
    d = _definition([dataclasses.replace(p, stationary=None) for p in CHEMISTRY])
    assert not steady_plan(d).substances


# --- the nephron kidney through the kink of the threshold ----------------------


def test_the_nephron_threshold_passes_a_transient(tmp_path: Path) -> None:
    """A dose through the steady state of the threshold reabsorption of the nephron kidney.

    At the steady state the flow out of the collecting duct is the filtered
    load, the kink of `min(leaving, fu*GFR*Car)*f(Car)`; a bolus of 2 g
    carries the arterial plasma above the renal threshold and back: the
    integrator passes the kink at the start and on the way back (no
    smoothing needed), the reabsorption never exceeds the filtration, and the
    body returns to its baseline.
    """
    d = dataclasses.replace(
        _pinned(Variants(organs={Organ.KIDNEY: "nephron"})),
        routes=[Route.iv(S)],
    )
    r = individual(build(d, tmp_path).flat_path, **{ids.iv_dose(S): 2000.0})
    r.integrator.absolute_tolerance = 1e-14
    r.integrator.relative_tolerance = 1e-12
    reab, gfr = (
        f"KI__{ids.tubule_reabsorption_reaction(S)}",
        f"KI__{ids.filtration_reaction(S)}",
    )
    values = np.asarray(
        r.simulate(0, 5 * DAY, 2001, selections=["time", f"[Car_{S}]", reab, gfr])
    )
    car = values[:, 1]
    assert np.all(np.isfinite(values))
    # above the renal threshold: fu*C above 0.1 mM
    assert car.max() * 0.5 > 0.1
    assert np.all(values[:, 2] <= values[:, 3] * (1 + 1e-12))
    assert car[-1] == pytest.approx(0.03, rel=1e-4)
