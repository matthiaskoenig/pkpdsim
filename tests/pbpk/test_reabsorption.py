"""The saturable reabsorption of the filtrate (`Excretion(reabsorption=ThresholdReabsorption(...))`).

The reabsorbed fraction of the filtrate at the plasma concentration `C` is
`1/(1 + (fu*C/threshold)^n)^(1/n)`: almost all of a filtrate below the renal
threshold, and the reabsorption approaches the transport maximum
`GFR*threshold` above it. The glomerulus filters the afferent plasma `Car`
and the plasma leaving the kidney lacks the net filtrate at `Cout`, as with
the constant fraction.
"""

import dataclasses
from pathlib import Path
from typing import TYPE_CHECKING, Any

import libsbml
import numpy as np
import pytest

from pkpdsim.pbpk import (
    Binding,
    CustomReaction,
    Direction,
    Excretion,
    FirstOrder,
    Formula,
    Kp,
    Membrane,
    ModelParameter,
    Organ,
    Part,
    PBPKDefinition,
    Pool,
    Route,
    Substance,
    ThresholdReabsorption,
    Transport,
    U,
    Value,
    Variants,
    build,
    ids,
)
from pkpdsim.pbpk import steady_state as ss

from .conftest import (
    assert_reference_sbml,
    chemistry_reference_definition,
    dose_amount,
    individual,
    nephron_reference_definition,
    total_amount,
)

if TYPE_CHECKING:
    import roadrunner

S = "sat"
FU = 0.5
THRESHOLD = 0.2  # mM, of the filtrate
N = 1.5  # not an integer: the power of a negative concentration needs abs
DAY = 24 * 60.0
STANDARD = Variants()
NEPHRON = Variants(organs={Organ.KIDNEY: "nephron"})


def _substance(**fields: Any) -> Substance:  # noqa: ANN401
    """A drug filtered with a saturable reabsorption, in the merged tissue."""
    values: dict[str, Any] = {
        "sid": S,
        "name": "reabsorbed drug",
        "mr": 180.0,
        "binding": Binding(fu=FU, protein="albumin"),
        "distribution": Kp(0.8, ftissue=2.0),
        "excretion": Excretion(
            urine=True, reabsorption=ThresholdReabsorption(threshold=THRESHOLD, n=N)
        ),
    }
    values.update(fields)
    return Substance(**values)


def _definition(variants: Variants = STANDARD, **fields: Any) -> PBPKDefinition:  # noqa: ANN401
    """The drug alone, given intravenously."""
    return PBPKDefinition(
        sid="reabsorption",
        substances=[_substance(**fields)],
        variants=variants,
        routes=[Route.iv(S)],
    )


@pytest.fixture(scope="module")
def standard(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The standard kidney with the saturable reabsorption, built once."""
    return build(_definition(), tmp_path_factory.mktemp("reab_standard")).flat_path


@pytest.fixture(scope="module")
def nephron(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The nephron kidney with the saturable reabsorption, built once."""
    return build(
        _definition(NEPHRON), tmp_path_factory.mktemp("reab_nephron")
    ).flat_path


def _fraction(c: float, fu: float = FU) -> float:
    """The reabsorbed fraction of the filtrate at the plasma concentration `c`."""
    return 1.0 / (1.0 + (fu * c / THRESHOLD) ** N) ** (1.0 / N)


def _load(path: Path, **changes: float) -> "roadrunner.RoadRunner":
    """Load a model with changes before its initialization, tight tolerances."""
    r = individual(path, **changes)
    r.integrator.absolute_tolerance = 1e-14
    r.integrator.relative_tolerance = 1e-12
    return r


CONCENTRATIONS = np.logspace(-6, 3, 46)
"""The plasma concentrations of the sweep (mM), far below to far above the threshold."""


def test_the_reabsorption_never_exceeds_the_filtration(standard: Path) -> None:
    """Over the sweep of the arterial plasma the reabsorption stays in `[0, filtration]`.

    It never exceeds the transport maximum `GFR*threshold` either, and it is
    the law of the threshold.
    """
    r = _load(standard)
    gfr = r[ids.GFR]
    arterial = ids.pool_species(Pool.ARTERIAL, S)
    filtration = ids.flat(Organ.KIDNEY, ids.filtration_reaction(S))
    for c in CONCENTRATIONS:
        r[f"[{arterial}]"] = c
        filtered = FU * gfr * c
        net = r[filtration]
        reabsorbed = filtered - net
        assert 0.0 <= net <= filtered
        assert 0.0 <= reabsorbed <= filtered
        assert reabsorbed <= gfr * THRESHOLD * (1 + 1e-12)
        assert reabsorbed == pytest.approx(filtered * _fraction(c), rel=1e-9)
        assert net == pytest.approx(filtered * (1 - _fraction(c)), rel=1e-9, abs=0)


def test_the_outflow_lacks_the_net_filtrate_of_the_outflow(standard: Path) -> None:
    """The plasma leaving the kidney lacks the net filtrate at its own concentration."""
    r = _load(standard)
    gfr, carrier = r[ids.GFR], r[ids.organ_flow(Organ.KIDNEY)] * (1 - r[ids.HCT])
    kidney = ids.flat(Organ.KIDNEY, ids.species(Part.PLASMA, S))
    outflow = ids.flat(Organ.KIDNEY, ids.outflow_reaction(Part.PLASMA, S))
    for c in CONCENTRATIONS:
        r[f"[{kidney}]"] = c
        cout = r[ids.organ_outflow(Organ.KIDNEY, S)]
        expected = carrier * cout - FU * gfr * cout * (1 - _fraction(cout))
        assert r[outflow] == pytest.approx(expected, rel=1e-12)
        assert r[outflow] > 0.0


def test_the_saturable_reabsorption_conserves_the_amount(
    standard: Path, nephron: Path
) -> None:
    """Body and urine keep the dose of an IV bolus, in both kidneys."""
    for path in (standard, nephron):
        r = _load(path)
        r[ids.iv_dose(S)] = 500.0
        r.simulate(0, 600, 2)
        amount = total_amount(r, S) + dose_amount(r, S)
        assert amount == pytest.approx(500.0 / r[ids.mr(S)], rel=1e-9)
        assert r[ids.urine(S)] > 0.0


def test_the_nephron_kidney_reabsorbs_the_fraction_of_the_filtered_load(
    nephron: Path,
) -> None:
    """The collecting duct returns the fraction of the afferent plasma, up to the filtered load.

    With little drug in the duct (the tubule filling) the fraction of the
    flow leaving it returns; with much (a secretion into the tubule) the
    fraction of the filtered load, never more than the filtration nor than
    the flow leaving the duct.
    """
    r = _load(nephron)
    gfr = r[ids.GFR]
    arterial = ids.pool_species(Pool.ARTERIAL, S)
    duct = ids.flat(Organ.KIDNEY, ids.tubule(S, "cd"))
    reab = ids.flat(Organ.KIDNEY, ids.tubule_reabsorption_reaction(S))
    urine = ids.flat(Organ.KIDNEY, ids.tubule_reaction(S, "cd"))
    f_cd = r[ids.flat(Organ.KIDNEY, ids.tubule_flow_fraction("cd"))]
    for duct_c in (1e-4, 1.0, 1e3):
        r[f"[{duct}]"] = duct_c
        leaving = gfr * f_cd * duct_c
        for c in CONCENTRATIONS:
            r[f"[{arterial}]"] = c
            filtered = FU * gfr * c
            expected = min(leaving, filtered) * _fraction(c)
            assert r[reab] == pytest.approx(expected, rel=1e-12)
            assert r[reab] + r[urine] == pytest.approx(leaving, rel=1e-12)
            assert 0.0 <= r[reab] <= min(leaving, filtered)


def test_a_secretion_is_not_reabsorbed_beyond_the_filtration(tmp_path: Path) -> None:
    """In the nephron kidney with a secretion the dose is conserved and the urine gets the secreted drug.

    At steady state of an infusion the reabsorption is the fraction of the
    filtered load, so the excretion is at least the secretion.
    """
    definition = dataclasses.replace(
        _definition(NEPHRON),
        chemistry=[
            Transport(
                "SATSEC",
                Organ.KIDNEY,
                S,
                Direction.EFFLUX,
                Membrane.APICAL,
                FirstOrder(k=0.2, per_volume=False),
            )
        ],
    )
    r = _load(build(definition, tmp_path).flat_path)
    r[ids.iv_dose(S)] = 500.0
    r.simulate(0, 600, 2)
    amount = total_amount(r, S) + dose_amount(r, S)
    assert amount == pytest.approx(500.0 / r[ids.mr(S)], rel=1e-9)
    gfr = r[ids.GFR]
    car = r[f"[{ids.pool_species(Pool.ARTERIAL, S)}]"]
    reabsorbed = r[ids.flat(Organ.KIDNEY, ids.tubule_reabsorption_reaction(S))]
    assert reabsorbed <= FU * gfr * car * _fraction(car) * (1 + 1e-12)


def test_the_law_is_defined_below_zero(standard: Path, nephron: Path) -> None:
    """A concentration a little below 0 (an undershoot) gives finite rates."""
    for path in (standard, nephron):
        r = _load(path)
        r[f"[{ids.pool_species(Pool.ARTERIAL, S)}]"] = -1e-12
        r[f"[{ids.flat(Organ.KIDNEY, ids.species(Part.PLASMA, S))}]"] = -1e-12
        rates = np.asarray(r.model.getReactionRates())
        assert np.all(np.isfinite(rates))


def test_without_a_threshold_reabsorption_the_nephron_sbml_is_unchanged(
    tmp_path: Path,
) -> None:
    """The nephron kidney with a constant fraction and a secretion is the SBML of 311895c."""
    assert_reference_sbml(build(nephron_reference_definition(), tmp_path))


def test_the_parameters_of_the_threshold(standard: Path, tmp_path: Path) -> None:
    """`reab_threshold_<s>` and `reab_n_<s>` replace `f_reabsorbed_<s>`; a fit sets them."""
    r = _load(standard)
    parameters = set(r.model.getGlobalParameterIds())
    assert {ids.reabsorption_threshold(S), ids.reabsorption_n(S)} <= parameters
    assert ids.fraction_reabsorbed(S) not in parameters
    assert r[ids.reabsorption_threshold(S)] == THRESHOLD
    assert r[ids.reabsorption_n(S)] == N
    fitted = ThresholdReabsorption(
        threshold=Value(0.3, U.mM, notes="The fitted threshold."),
        n=Value(3.0, U.dimensionless),
    )
    path = build(
        _definition(excretion=Excretion(urine=True, reabsorption=fitted)), tmp_path
    ).flat_path
    model = libsbml.readSBMLFromFile(str(path)).getModel()
    threshold = model.getParameter(ids.reabsorption_threshold(S))
    assert threshold.getValue() == 0.3
    assert "The fitted threshold." in threshold.getNotesString()
    assert model.getParameter(ids.reabsorption_n(S)).getValue() == 3.0


def test_without_a_saturable_reabsorption_the_sbml_is_unchanged(
    tmp_path: Path,
) -> None:
    """The default is no saturable reabsorption, and the SBML is that of 0c538ca."""
    assert Excretion().reabsorption is None
    assert_reference_sbml(build(chemistry_reference_definition(), tmp_path))


@pytest.mark.parametrize(
    ("excretion", "message"),
    [
        (
            Excretion(urine=False, reabsorption=ThresholdReabsorption(THRESHOLD, N)),
            r"sat \(saturable reabsorption\): no filtrate to reabsorb",
        ),
        (
            Excretion(
                urine=True,
                fraction_reabsorbed=0.2,
                reabsorption=ThresholdReabsorption(THRESHOLD, N),
            ),
            r"fraction_reabsorbed must stay 0",
        ),
        (
            Excretion(urine=True, reabsorption=ThresholdReabsorption(0.0, N)),
            r"threshold must be positive",
        ),
        (
            Excretion(urine=True, reabsorption=ThresholdReabsorption(THRESHOLD, -1.0)),
            r"n must be positive",
        ),
        (
            Excretion(
                urine=True, reabsorption=ThresholdReabsorption(Value(0.2, U.mmole), N)
            ),
            r"threshold must be in mmole/liter, got mmole",
        ),
        (
            Excretion(
                urine=True,
                reabsorption=ThresholdReabsorption(THRESHOLD, Value(2.0, U.mM)),
            ),
            r"n must be in dimensionless",
        ),
    ],
)
def test_an_invalid_reabsorption_is_refused(excretion: Excretion, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        _definition(excretion=excretion).validate()


# --- the steady state ---------------------------------------------------------

BASELINE = 0.15  # mM, near the threshold of the filtrate (fu*C = 0.075 mM)


def _endogenous(variants: Variants = STANDARD) -> PBPKDefinition:
    """An endogenous substance produced in the muscle and excreted only by the kidney.

    The production `k_prod` is the renal elimination at the baseline,
    `renal_clearance*baseline`; the muscle starts at its steady state with
    the production, the kidney at the steady state of the saturable
    reabsorption (`kidney_plasma`, the arterial plasma), the other organs at
    the equilibrium of the baseline.
    """
    s = _substance(baseline=BASELINE)
    clearance = ss.renal_clearance(s, variants=variants)
    kidney = ss.kidney_plasma(s, variants=variants)
    muscle = ss.organ_plasma(s, Organ.MUSCLE, rate="k_prod", variants=variants)
    return PBPKDefinition(
        sid="reabsorption_steady",
        substances=[s],
        chemistry=[
            CustomReaction(
                "PROD", Organ.MUSCLE, Part.TISSUE, f"-> Cti_{S}", "k_prod", []
            ),
        ],
        parameters=[
            ModelParameter(
                "k_prod",
                Formula(f"({clearance}) * {ids.baseline(S)}", U.mmole_per_min),
                U.mmole_per_min,
            ),
        ],
        starts={
            f"MU__Cpl_{S}": Formula(muscle, U.mM),
            f"MU__Cti_{S}": Formula(
                ss.organ_tissue(
                    s, Organ.MUSCLE, muscle, rate="k_prod", variants=variants
                ),
                U.mM,
            ),
            f"KI__Cpl_{S}": Formula(kidney, U.mM),
            f"KI__Cti_{S}": Formula(
                ss.organ_tissue(s, Organ.KIDNEY, kidney, variants=variants), U.mM
            ),
        },
        routes=[Route.iv(S)],
    )


def _drift(r: "roadrunner.RoadRunner") -> dict[str, float]:
    """The largest relative change of every species but the urine over a day."""
    sids = [s for s in r.model.getFloatingSpeciesIds() if s != ids.urine(S)]
    start = np.array([r[f"[{sid}]"] for sid in sids])
    values = np.asarray(r.simulate(0, DAY, 97, selections=[f"[{sid}]" for sid in sids]))
    scale = np.where(start > 0, start, 1.0)
    return dict(zip(sids, np.max(np.abs(values - start), axis=0) / scale, strict=True))


@pytest.mark.parametrize(
    "changes",
    [
        {},
        {"BW": 100.0, "HEIGHT": 190.0},
        {"f_renal_function": 0.3},
        {"f_renal_function": 1.8, "HF_severity": 2.0},
    ],
    ids=["reference", "large", "ckd", "hyperfiltration_hf"],
)
def test_the_steady_state_of_a_saturable_reabsorption(
    tmp_path_factory: pytest.TempPathFactory, changes: dict[str, float]
) -> None:
    """With the starts of the steady state every species stays at its start."""
    path = _steady_path(tmp_path_factory)
    r = _load(path, **changes)
    gfr, fu = r[ids.GFR], r[ids.fu(S)]
    assert r[f"[KI__Cpl_{S}]"] == r[ids.baseline(S)] == BASELINE
    # the renal clearance of the arterial plasma is the net filtered plasma flow
    assert r["k_prod"] == pytest.approx(
        fu * gfr * (1 - _fraction(BASELINE, fu)) * BASELINE, rel=1e-12
    )
    drift = _drift(r)
    assert max(drift.values()) < 1e-9, sorted(drift.items(), key=lambda x: -x[1])[:3]


_STEADY: dict[str, Path] = {}


def _steady_path(factory: pytest.TempPathFactory) -> Path:
    """The model of the steady state, built once per process."""
    if "path" not in _STEADY:
        _STEADY["path"] = build(_endogenous(), factory.mktemp("reab_steady")).flat_path
    return _STEADY["path"]


def test_the_steady_state_helpers_of_a_saturable_reabsorption() -> None:
    """The net filtrate reads the arterial plasma; what has no closed form is refused."""
    s = _substance(baseline=BASELINE)
    fraction = (
        "1 dimensionless / (1 dimensionless + (fu_sat * abs(baseline_sat) / "
        "reab_threshold_sat)^reab_n_sat)^(1 dimensionless / reab_n_sat)"
    )
    assert ss.renal_filtration(s, variants=STANDARD) == (
        f"fu_sat * GFR * (1 dimensionless - {fraction})"
    )
    assert "abs(Car_x)" in ss.renal_filtration(s, arterial="Car_x", variants=STANDARD)
    assert ss.renal_clearance(s, variants=STANDARD) == ss.renal_filtration(
        s, variants=STANDARD
    )
    assert ss.kidney_plasma(s, variants=STANDARD) == ids.baseline(S)
    assert ss.kidney_plasma(s, arterial="C_a", variants=STANDARD) == "C_a"
    assert ss.kidney_arterial(s, "C_k", variants=STANDARD) == "C_k"
    refused = [
        lambda: ss.kidney_plasma(s, rate="r", variants=STANDARD),
        lambda: ss.kidney_plasma(s, clearance="cl", variants=STANDARD),
        lambda: ss.kidney_plasma(
            s, secretion=ss.Saturable("v", "km"), variants=STANDARD
        ),
        lambda: ss.kidney_plasma(s, eliminated="e", variants=STANDARD),
        lambda: ss.kidney_arterial(s, "C_k", rate="r", variants=STANDARD),
        lambda: ss.kidney_arterial(s, "C_k", eliminated="e", variants=STANDARD),
        lambda: ss.renal_clearance(s, clearance="cl", variants=STANDARD),
    ]
    for call in refused:
        with pytest.raises(NotImplementedError, match="saturable reabsorption"):
            call()
    with pytest.raises(NotImplementedError, match="nephron"):
        ss.kidney_plasma(s, variants=NEPHRON)
    constant = dataclasses.replace(s, excretion=Excretion(urine=True))
    assert ss.renal_filtration(constant, arterial="Car_x", variants=STANDARD) == (
        "fu_sat * GFR * (1 dimensionless - f_reabsorbed_sat)"
    )
