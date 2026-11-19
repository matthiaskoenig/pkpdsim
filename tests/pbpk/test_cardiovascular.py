"""Pressure-driven circulation against PBPK, conservation and nodal benchmarks."""

import dataclasses
import warnings
from pathlib import Path
from typing import Literal

import libsbml
import numpy as np
import pytest
from sbmlsim.engine import Engine
from sbmlsim.simulation import Change

from pkpdsim.pbpk import (
    DynamicHeart,
    Erythrocytes,
    Formula,
    Kp,
    Mode,
    ModelParameter,
    NephronKidney,
    Organ,
    PBPKDefinition,
    Pool,
    PressureFlow,
    Sex,
    SingleLumen,
    Substance,
    U,
    Variants,
    ZonatedLiver,
    build,
    ids,
)
from pkpdsim.pbpk.builder import BuildResult
from pkpdsim.pbpk.experiment import healthy_gfr
from pkpdsim.pbpk.flat_model import loaded_model
from pkpdsim.pbpk.kidney_validity import individual_values
from pkpdsim.pbpk.organs import PORTAL_ORGANS
from pkpdsim.pbpk.stationary import SteadyStateWarning

from .cardiovascular_solver import simulate
from .conftest import cardiac_pd

TRACER = "trc"
FLOWS = ["CO", ids.QPO, ids.QH, *(ids.organ_flow(o) for o in Organ)]


def _definition(circulation: PressureFlow | DynamicHeart | None) -> PBPKDefinition:
    return PBPKDefinition(
        sid="cv_tracer",
        substances=[
            Substance(
                TRACER,
                "closed circulation tracer",
                100.0,
                distribution=Kp(1.0, mode=Mode.PERFUSION),
            )
        ],
        variants=Variants(circulation=circulation),
        starts={ids.pool_species(Pool.VENOUS, TRACER): 1.0},
    )


@pytest.fixture(scope="module")
def baseline(tmp_path_factory: pytest.TempPathFactory) -> BuildResult:
    return build(_definition(None), tmp_path_factory.mktemp("cv_baseline"))


@pytest.fixture(scope="module")
def pressure(tmp_path_factory: pytest.TempPathFactory) -> BuildResult:
    return build(_definition(PressureFlow()), tmp_path_factory.mktemp("cv_pressure"))


@pytest.fixture(scope="module")
def heart(tmp_path_factory: pytest.TempPathFactory) -> BuildResult:
    return build(_definition(DynamicHeart()), tmp_path_factory.mktemp("cv_heart"))


def _species(path: Path) -> list[str]:
    doc = libsbml.readSBMLFromFile(str(path))
    all_ids = [sp.getId() for sp in doc.getModel().getListOfSpecies()]
    return ids.species_of(TRACER, all_ids)


def test_cardiovascular_species_discovery_matches_the_exact_substance() -> None:
    expected = [
        f"CV_C_{node}{part}_{TRACER}"
        for node in ("pa", "pv", "lv", "rv")
        for part in ("", "_bc")
    ]
    other = [
        f"CV_C_pa_{TRACER}_m2",
        f"CV_C_pa_x{TRACER}",
        f"CV_C_invalid_{TRACER}",
        f"CV_P_pa_{TRACER}",
    ]
    assert ids.species_of(TRACER, [*expected, *other]) == expected


def test_fractional_flows_remain_the_default(baseline: BuildResult) -> None:
    doc = libsbml.readSBMLFromFile(str(baseline.flat_path))
    model = doc.getModel()
    assert not any(p.getId().startswith("CV_") for p in model.getListOfParameters())
    assert not any(rule.isAlgebraic() for rule in model.getListOfRules())


@pytest.mark.parametrize("changes", [{}, {"BW": 100.0, "HEIGHT": 180.0}])
def test_pressure_flows_reproduce_reference_perfusion(
    baseline: BuildResult,
    pressure: BuildResult,
    changes: dict[str, float],
) -> None:
    expected = simulate(baseline.flat_path, FLOWS, preinit_changes=changes)
    actual = simulate(pressure.flat_path, FLOWS, preinit_changes=changes)
    for sid in FLOWS:
        np.testing.assert_allclose(actual[sid], expected[sid], rtol=1e-7, atol=1e-9)


def test_reference_tracer_matches_fractional_flow_pharmacokinetics(
    baseline: BuildResult, pressure: BuildResult
) -> None:
    species = _species(baseline.flat_path)
    expected = simulate(baseline.flat_path, species, end=5.0)
    actual = simulate(pressure.flat_path, species, end=5.0)
    for sid in species:
        np.testing.assert_allclose(actual[sid], expected[sid], rtol=1e-6, atol=1e-9)
    total = np.sum([actual[sid] for sid in species], axis=0)
    assert total[0] > 0
    np.testing.assert_allclose(total, total[0], rtol=1e-8, atol=1e-9)


def test_local_resistance_changes_its_own_perfusion(pressure: BuildResult) -> None:
    reference = simulate(pressure.flat_path, FLOWS)
    changed = simulate(pressure.flat_path, FLOWS, preinit_changes={"CV_fR_mu": 2.0})
    np.testing.assert_allclose(changed["Qmu"], reference["Qmu"] / 2, rtol=1e-7)
    for sid in ("Qki", ids.QPO, ids.QH):
        np.testing.assert_allclose(changed[sid], reference[sid], rtol=1e-7)
    assert changed["CO"][-1] < reference["CO"][-1]


def test_resistance_change_during_simulation_preserves_tracer(
    pressure: BuildResult,
) -> None:
    species = _species(pressure.flat_path)
    result = simulate(
        pressure.flat_path,
        [*species, "CV_P_po", "CV_P_li", ids.QPO, ids.QH, "Qki"],
        end=1.0,
        changes=[Change(0.5, {"CV_fR_hepatic": 4.0})],
    )
    expected = simulate(
        pressure.flat_path,
        ["CV_P_po", "CV_P_li", ids.QPO, ids.QH],
        preinit_changes={"CV_fR_hepatic": 4.0},
    )
    assert result["CV_P_po"][-1] > result["CV_P_po"][0]
    assert result[ids.QPO][-1] < result[ids.QPO][0]
    for sid in ("CV_P_po", "CV_P_li", ids.QPO, ids.QH):
        np.testing.assert_allclose(result[sid][-1], expected[sid][0], rtol=1e-7)
    np.testing.assert_allclose(result["Qki"][-1], result["Qki"][0], rtol=1e-7)
    total = np.sum([result[sid] for sid in species], axis=0)
    np.testing.assert_allclose(total, total[0], rtol=1e-8, atol=1e-9)


@pytest.mark.parametrize("hepatic", [1.0, 2.0, 4.0])
def test_hepatic_network_matches_independent_nodal_solve(
    pressure: BuildResult, hepatic: float
) -> None:
    resistance_ids = [
        "CV_R_portal",
        "CV_R_hepatic",
        "CV_R_li",
        *(f"CV_R_{organ.value}" for organ in PORTAL_ORGANS),
    ]
    result = simulate(
        pressure.flat_path,
        [*resistance_ids, "CV_P_ar", "CV_P_ve", "CV_P_po", "CV_P_li", ids.QPO, ids.QH],
        preinit_changes={
            "CV_fR_hepatic": hepatic,
            "CV_fR_portal": 1.3,
            "CV_fR_si": 2.0,
            "CV_fR_li": 1.4,
        },
    )
    values = {sid: result[sid][0] for sid in result.columns}
    arterial, venous = values["CV_P_ar"], values["CV_P_ve"]
    splanchnic = sum(
        1.0
        / (
            values[f"CV_R_{organ.value}"]
            * (2.0 if organ is Organ.SMALL_INTESTINE else 1.0)
        )
        for organ in PORTAL_ORGANS
    )
    portal = 1.0 / (values["CV_R_portal"] * 1.3)
    artery = 1.0 / (values["CV_R_li"] * 1.4)
    outflow = 1.0 / (values["CV_R_hepatic"] * hepatic)
    matrix = np.array(
        [[splanchnic + portal, -portal], [-portal, portal + artery + outflow]]
    )
    rhs = np.array([splanchnic * arterial, artery * arterial + outflow * venous])
    portal_pressure, sinusoidal_pressure = np.linalg.solve(matrix, rhs)
    np.testing.assert_allclose(result["CV_P_po"], portal_pressure, rtol=1e-7)
    np.testing.assert_allclose(result["CV_P_li"], sinusoidal_pressure, rtol=1e-7)
    np.testing.assert_allclose(
        result[ids.QPO], portal * (portal_pressure - sinusoidal_pressure), rtol=1e-7
    )
    np.testing.assert_allclose(
        result[ids.QH], outflow * (sinusoidal_pressure - venous), rtol=1e-7
    )


def test_hepatic_resistance_raises_portal_pressure_and_reduces_flow(
    pressure: BuildResult,
) -> None:
    selections = ["CV_P_po", "CV_P_li", ids.QPO, ids.QH, "Qki"]
    reference = simulate(pressure.flat_path, selections)
    changed = simulate(
        pressure.flat_path, selections, preinit_changes={"CV_fR_hepatic": 4.0}
    )
    assert changed["CV_P_po"][-1] > reference["CV_P_po"][-1]
    assert changed["CV_P_li"][-1] > reference["CV_P_li"][-1]
    assert changed[ids.QPO][-1] < reference[ids.QPO][-1]
    assert changed[ids.QH][-1] < reference[ids.QH][-1]
    np.testing.assert_allclose(changed["Qki"], reference["Qki"], rtol=1e-7)


@pytest.mark.parametrize("value", [0.0, -1.0, float("inf"), float("nan")])
def test_invalid_resistances_are_refused(value: float) -> None:
    with pytest.raises(ValueError, match="resistance"):
        PressureFlow(hepatic_resistance=value)
    with pytest.raises(ValueError, match="resistance"):
        PressureFlow(resistance={Organ.KIDNEY: value})


def test_individual_values_and_healthy_gfr_use_the_dae_engine(
    pressure: BuildResult, baseline: BuildResult
) -> None:
    assert loaded_model(pressure.flat_path).engine is Engine.CORE
    assert loaded_model(baseline.flat_path).engine is Engine.ROADRUNNER
    changes = {"BW": 100.0, "HEIGHT": 180.0, "CV_fR_ki": 2.0}
    expected = simulate(
        pressure.flat_path, ["GFR_base", "Qki"], preinit_changes=changes
    )
    values = individual_values(pressure.flat_path, changes)
    assert values.number("Qki") == pytest.approx(expected["Qki"][0], rel=1e-8)
    assert values.number("GFR_base") == pytest.approx(expected["GFR_base"][0], rel=1e-8)
    gfr = healthy_gfr(pressure.flat_path, bw=100.0, height=180.0)
    assert float(gfr.m_as("ml/min")) == pytest.approx(
        expected["GFR_base"][0] * 1000, rel=1e-8
    )


def test_dynamic_heart_is_pulsatile_and_conserves_blood_and_tracer(
    heart: BuildResult,
) -> None:
    species = _species(heart.flat_path)
    stores = [*(ids.pool_volume(pool) for pool in Pool)] + [
        f"CV_V_{node}" for node in ("pa", "pv", "lv", "rv")
    ]
    result = simulate(
        heart.flat_path,
        [*species, *stores, "CV_V_total", "CV_V_blood_total", "Vblood", "CV_P_ar"],
        end=0.5,
        steps=1000,
        maximum_time_step=1 / 6000,
    )
    assert np.all(np.isfinite(result.values))
    assert np.ptp(result["CV_P_ar"][result.time > 0.4]) > 1.0
    for sid in stores:
        assert np.all(result[sid] > 0), sid
    volume = np.sum([result[sid] for sid in stores], axis=0)
    np.testing.assert_allclose(volume, volume[0], rtol=1e-8, atol=1e-9)
    np.testing.assert_allclose(volume, result["CV_V_total"], rtol=1e-8, atol=1e-9)
    np.testing.assert_allclose(
        result["CV_V_blood_total"], result["Vblood"], rtol=1e-8, atol=1e-9
    )
    amount = np.sum([result[sid] for sid in species], axis=0)
    assert amount[0] > 0
    np.testing.assert_allclose(amount, amount[0], rtol=1e-7, atol=1e-9)


def test_dynamic_heart_initial_renal_values_use_consistent_dae_state(
    heart: BuildResult,
) -> None:
    changes = {"BW": 85.0, "HEIGHT": 170.0, "CV_fR_ki": 2.0}
    expected = simulate(
        heart.flat_path,
        ["Qki", "GFR_base"],
        end=0.001,
        steps=1,
        preinit_changes=changes,
    )
    values = individual_values(heart.flat_path, changes)
    assert values.number("Qki") == pytest.approx(expected["Qki"][0], rel=1e-8)
    assert values.number("GFR_base") == pytest.approx(expected["GFR_base"][0], rel=1e-8)
    assert float(
        healthy_gfr(heart.flat_path, bw=85.0, height=170.0).m_as("ml/min")
    ) == (pytest.approx(expected["GFR_base"][0] * 1000, rel=1e-8))


def test_drug_feedback_changes_resistance_perfusion_and_pharmacokinetics(
    tmp_path: Path, pressure: BuildResult
) -> None:
    component = dataclasses.replace(
        cardiac_pd(emax=0.5, ec50=0.03, reads=TRACER),
        outputs={"f_cardiac": "CV_fR_systemic"},
    )
    definition = dataclasses.replace(_definition(PressureFlow()), pd=[component])
    feedback = build(definition, tmp_path)
    arm = ids.arm_free(TRACER)
    venous = f"[{ids.pool_species(Pool.VENOUS, TRACER)}]"
    selections = [arm, venous, "CV_fR_systemic", "Qmu"]
    reference = simulate(pressure.flat_path, selections, end=1.0)
    result = simulate(feedback.flat_path, selections, end=1.0)
    np.testing.assert_allclose(
        result["CV_fR_systemic"], 1.0 - 0.5 * result[arm] / (0.03 + result[arm])
    )
    assert np.ptp(result["CV_fR_systemic"]) > 0.1
    np.testing.assert_allclose(
        result["Qmu"] * result["CV_fR_systemic"], reference["Qmu"], rtol=1e-7
    )
    assert np.max(np.abs(result[venous] - reference[venous])) > 1e-4


def test_combine_merges_compatible_circulation_requirements() -> None:
    parameters = {"P_ar": 100.0}
    multipliers = {Organ.KIDNEY: 2.0}
    choice = PressureFlow(parameters=parameters, resistance=multipliers)
    parameters["P_ar"] = 110.0
    multipliers[Organ.KIDNEY] = 3.0
    equal = PressureFlow(parameters={"P_ar": 100.0}, resistance={Organ.KIDNEY: 2.0})
    assert choice == equal
    assert choice.parameters["P_ar"] == 100.0
    assert choice.resistance[Organ.KIDNEY] == 2.0
    for definitions in (
        (_definition(None), _definition(choice)),
        (_definition(choice), _definition(None)),
        (_definition(choice), _definition(equal)),
    ):
        combined = PBPKDefinition.combine(*definitions, sid="cv_combined")
        assert combined.variants.circulation == equal
    with pytest.raises(ValueError, match="circulation"):
        PBPKDefinition.combine(
            _definition(PressureFlow()), _definition(DynamicHeart()), sid="conflict"
        )


@pytest.mark.parametrize("mode", ["equilibrium", "kinetic"])
def test_dynamic_erythrocyte_transport_conserves_tracer(
    tmp_path: Path, mode: Literal["equilibrium", "kinetic"]
) -> None:
    definition = _definition(DynamicHeart())
    substance = dataclasses.replace(
        definition.substances[0],
        distribution=Kp(
            1.0,
            mode=Mode.TISSUE if mode == "kinetic" else Mode.PERFUSION,
            ftissue=1.0 if mode == "kinetic" else None,
        ),
        erythrocytes=Erythrocytes(
            mode=mode, kbc=0.7, k_exchange=5.0 if mode == "kinetic" else None
        ),
    )
    built = build(dataclasses.replace(definition, substances=[substance]), tmp_path)
    species = _species(built.flat_path)
    result = simulate(
        built.flat_path,
        species,
        end=0.1,
        steps=200,
        maximum_time_step=1 / 6000,
    )
    amount = np.sum([result[sid] for sid in species], axis=0)
    assert amount[0] > 0
    np.testing.assert_allclose(amount, amount[0], rtol=1e-7, atol=1e-9)


@pytest.mark.parametrize("imbalance", [False, True])
def test_pressure_model_checks_pinned_endogenous_balance(
    tmp_path: Path, imbalance: bool
) -> None:
    from pkpdsim.pbpk import CustomReaction, Part

    definition = PBPKDefinition(
        sid="cv_endogenous",
        substances=[
            Substance(
                "end",
                "endogenous",
                100.0,
                distribution=Kp(1.0, ftissue=1.0),
                baseline=0.03,
            )
        ],
        variants=Variants(circulation=PressureFlow()),
        parameters=[
            ModelParameter(
                "source_rate", 0.004 if imbalance else 0.003, U.mmole_per_min
            ),
            ModelParameter("clearance", 0.1, U.l_per_min),
        ],
        chemistry=[
            CustomReaction(
                "SOURCE",
                Organ.MUSCLE,
                Part.PLASMA,
                "-> Cpl_end",
                "source_rate",
                stationary=Formula("source_rate", U.mmole_per_min),
            ),
            CustomReaction(
                "LOSS",
                Organ.MUSCLE,
                Part.PLASMA,
                "Cpl_end ->",
                "clearance*Cpl_end",
                stationary="linear",
            ),
        ],
    )
    if imbalance:
        with pytest.warns(SteadyStateWarning, match="do not balance"):
            build(definition, tmp_path)
    else:
        with warnings.catch_warnings():
            warnings.simplefilter("error", SteadyStateWarning)
            built = build(definition, tmp_path)
        result = simulate(built.flat_path, ["[Cve_end]"], end=5.0)
        np.testing.assert_allclose(result["[Cve_end]"], 0.03, rtol=1e-8)


@pytest.mark.parametrize("circulation", [PressureFlow(), DynamicHeart()])
def test_pressure_circulation_supports_combined_organ_variants(
    tmp_path: Path, circulation: PressureFlow
) -> None:
    definition = _definition(circulation)
    definition = dataclasses.replace(
        definition,
        variants=Variants(
            circulation=circulation,
            organs={
                Organ.LIVER: ZonatedLiver(n_zones=2),
                Organ.KIDNEY: NephronKidney(),
                Organ.SMALL_INTESTINE: SingleLumen(),
            },
        ),
        substances=[
            dataclasses.replace(
                definition.substances[0], distribution=Kp(1.0, ftissue=1.0)
            )
        ],
    )
    built = build(definition, tmp_path)
    species = _species(built.flat_path)
    result = simulate(
        built.flat_path,
        [*species, ids.QPO, "CV_P_li"],
        end=0.03,
        steps=60,
        changes=[Change(0.015, {"CV_fR_hepatic": 4.0})],
        maximum_time_step=1 / 6000,
    )
    amount = np.sum([result[sid] for sid in species], axis=0)
    assert amount[0] > 0
    np.testing.assert_allclose(amount, amount[0], rtol=1e-7, atol=1e-9)
    event_rows = np.flatnonzero(np.isclose(result.time, 0.015, rtol=0, atol=1e-12))
    assert len(event_rows) == 2
    before, after = event_rows
    assert result[ids.QPO][after] < result[ids.QPO][before]
    assert result["CV_P_li"][after] > result["CV_P_li"][before]


def test_female_dynamic_heart_preserves_reference_blood_volume(tmp_path: Path) -> None:
    definition = dataclasses.replace(_definition(DynamicHeart()), sex=Sex.FEMALE)
    built = build(definition, tmp_path)
    result = simulate(
        built.flat_path,
        ["CV_V_blood_total", "Vblood", "CV_P_ar"],
        end=0.03,
        steps=60,
        maximum_time_step=1 / 6000,
    )
    np.testing.assert_allclose(
        result["CV_V_blood_total"], result["Vblood"], rtol=1e-8, atol=1e-9
    )
    assert np.ptp(result["CV_P_ar"]) > 1.0


@pytest.mark.parametrize("hepatic", [1.0, 4.0])
@pytest.mark.parametrize("saturable", [False, True])
def test_generated_pressure_baselines_solve_the_consistent_initial_flow_network(
    tmp_path: Path, hepatic: float, saturable: bool
) -> None:
    from .test_stationary import CHEMISTRY, LICL, PROD, SEC
    from .test_stationary import _definition as stationary_definition

    definition = stationary_definition(
        [PROD, SEC, LICL] if saturable else CHEMISTRY,
        variants=Variants(circulation=PressureFlow(hepatic_resistance=hepatic)),
    )
    built = build(definition, tmp_path)
    doc = libsbml.readSBMLFromFile(str(built.flat_path))
    model = doc.getModel()
    species = ids.species_of("end", [sp.getId() for sp in model.getListOfSpecies()])
    snapshots = ["baseline_end", "CLss_end", "Rss_end"]
    if saturable:
        snapshots.append("Css_ki_end")
    for sid in snapshots:
        assert model.getInitialAssignmentBySymbol(sid) is not None
        assert model.getRuleByVariable(sid) is None
    result = simulate(built.flat_path, [*species, *snapshots], end=5.0)
    assert result["baseline_end"][0] > 0
    for sid in species:
        np.testing.assert_allclose(result[sid], result[sid][0], rtol=1e-7, atol=1e-10)
    intervention = simulate(
        built.flat_path,
        [*snapshots, ids.QH],
        end=0.1,
        steps=20,
        changes=[Change(0.05, {"CV_fR_hepatic": 2.0 * hepatic})],
    )
    assert intervention[ids.QH][-1] < intervention[ids.QH][0]
    for sid in snapshots:
        np.testing.assert_allclose(
            intervention[sid], intervention[sid][0], rtol=1e-10, atol=1e-12
        )


@pytest.mark.parametrize("target", ["CV_fR_systemic", "CV_fR_hepatic"])
def test_generated_pressure_baselines_refuse_circular_pharmacodynamic_flow_feedback(
    target: str,
) -> None:
    from .test_stationary import LICL, PROD, SEC
    from .test_stationary import _definition as stationary_definition

    feedback = dataclasses.replace(
        cardiac_pd(reads="end"), outputs={"f_cardiac": target}
    )
    definition = stationary_definition(
        [PROD, SEC, LICL],
        variants=Variants(circulation=PressureFlow()),
        pd=[feedback],
    )
    from pkpdsim.pbpk import compose

    with pytest.raises(ValueError, match=r"change during|PD output"):
        compose(definition)


def test_user_formulas_still_cannot_snapshot_pressure_driven_flows() -> None:
    from pkpdsim.pbpk import compose

    definition = dataclasses.replace(
        _definition(PressureFlow()),
        parameters=[
            ModelParameter("user_flow", Formula("Qki", U.l_per_min), U.l_per_min)
        ],
    )
    with pytest.raises(ValueError, match="change during"):
        compose(definition)
