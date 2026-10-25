"""Tests of the impaired organ function: Child-Pugh, CKD and heart failure."""

import dataclasses
import math
from pathlib import Path
from typing import Any

import pytest
import roadrunner
from sbmlutils.factory import Model, create_model

from pkpdsim import impairment as presets
from pkpdsim.pbpk import (
    MODEL_UNITS,
    Binding,
    Metabolism,
    Organ,
    PBPKDefinition,
    Route,
    Substance,
    U,
    build,
    ids,
    load,
)
from pkpdsim.pbpk.builder import BuildResult
from pkpdsim.pbpk.chemistry import MM
from pkpdsim.pbpk.impairment import (
    impairment_objects,
    physiology_factors,
    unbound_fraction,
)
from pkpdsim.pbpk.organs import PORTAL_ORGANS
from pkpdsim.pbpk.physiology import Sex, physiology_objects

from .conftest import definition, individual, simple_drug
from .conftest import load as load_model

PARALLEL = [o for o in Organ if o is not Organ.LUNG]
DATA = load("impairment")
CP = DATA["child_pugh"]
CKD = DATA["ckd"]
HF = DATA["heart_failure"]
HEALTHY_HCT = CP["healthy_hematocrit"]["value"]


def _impaired_definition(**overrides: Any) -> PBPKDefinition:  # noqa: ANN401
    """The simple definition with a CYP2D6 metabolism of the liver besides CYP3A4."""
    d = definition(**overrides)
    extra = Metabolism(
        sid="DRG2D6",
        organ=Organ.LIVER,
        substrate="drg",
        products={"drg_m": 1.0},
        kinetics=MM(vmax=1e-3, km=0.01),
        protein="CYP2D6",
    )
    return dataclasses.replace(d, chemistry=[*d.chemistry, extra])


@pytest.fixture(scope="module")
def model(tmp_path_factory: pytest.TempPathFactory) -> BuildResult:
    """The simple definition with the liver enzymes CYP3A4 and CYP2D6, built once."""
    return build(_impaired_definition(), tmp_path_factory.mktemp("impairment"))


def _flows(r: roadrunner.RoadRunner) -> dict[Organ, float]:
    return {o: r[ids.organ_flow(o)] for o in PARALLEL}


def test_child_pugh_class_values(model: BuildResult) -> None:
    """At `CP_severity` 1, 2, 3 every factor of cirrhosis equals the data."""
    for k in (1, 2, 3):
        r = individual(model.flat_path, CP_severity=float(k))
        for key in (
            "functional_liver",
            "hepatic_artery_flow",
            "kidney_flow",
            "cardiac_index",
            "gfr",
            "albumin",
            "agp",
        ):
            assert r[ids.condition_factor("cp", key)] == pytest.approx(
                CP[key]["value"][k - 1], rel=1e-12
            ), key
        assert r[ids.condition_factor("cp", "hematocrit")] == pytest.approx(
            CP["hematocrit"]["value"][k - 1] / HEALTHY_HCT, rel=1e-12
        )
        # the MathML of the interpolation has 15 significant digits
        assert r[ids.F_SHUNT] == pytest.approx(CP["shunt"]["value"][k - 1], rel=1e-14)
        assert r[ids.F_FUNCTIONAL] == pytest.approx(
            CP["functional_liver"]["value"][k - 1], rel=1e-12
        )
        for protein in ("CYP3A4", "CYP2D6"):
            expected = CP["proteins"][protein]["value"][k - 1]
            assert r[ids.condition_factor("cp", protein)] == pytest.approx(expected)
            assert r[ids.impairment_activity(protein, Organ.LIVER)] == pytest.approx(
                expected
            )


def test_shunt_is_the_cirrhosis_of_the_old_models() -> None:
    """The shunt of the Child-Pugh classes is `pkpdsim.impairment.CIRRHOSIS`, bit-identical."""
    assert CP["shunt"]["value"] == list(presets.CIRRHOSIS.values())[1:]


def test_child_pugh_interpolates(model: BuildResult) -> None:
    """Between the classes every factor interpolates linearly, outside it is constant."""
    albumin = [1.0, *CP["albumin"]["value"]]
    shunt = [0.0, *CP["shunt"]["value"]]
    for severity in (0.5, 1.5, 2.25, 3.5):
        r = individual(model.flat_path, CP_severity=severity)
        i = min(int(severity), 2)
        t = min(severity, 3.0) - i
        for key, points in (("albumin", albumin), ("shunt", shunt)):
            expected = points[i] + t * (points[i + 1] - points[i])
            assert r[ids.condition_factor("cp", key)] == pytest.approx(expected), (
                key,
                severity,
            )


def _assert_flows_sum(r: roadrunner.RoadRunner) -> None:
    """The flows in parallel sum to the cardiac output, which perfuses the lung."""
    assert sum(_flows(r).values()) == pytest.approx(r[ids.CO], rel=1e-12)
    assert r[ids.organ_flow(Organ.LUNG)] == pytest.approx(r[ids.CO], rel=1e-12)
    assert r[ids.QPO] == pytest.approx(
        sum(r[ids.organ_flow(o)] for o in PORTAL_ORGANS), rel=1e-12
    )


def test_flows_renormalised_under_impairment(model: BuildResult) -> None:
    """The flows sum to the cardiac output for every class and condition.

    Cirrhosis: the cardiac output is the healthy one times the cardiac index,
    the hepatic artery and the kidney have their factors, portal organs and
    brain are unchanged.
    """
    healthy = individual(model.flat_path)
    _assert_flows_sum(healthy)
    q0 = _flows(healthy)
    for k in (1, 2, 3):
        r = individual(model.flat_path, CP_severity=float(k))
        _assert_flows_sum(r)
        assert r[ids.CO] == pytest.approx(
            healthy[ids.CO] * CP["cardiac_index"]["value"][k - 1], rel=1e-12
        )
        q = _flows(r)
        assert q[Organ.LIVER] / q0[Organ.LIVER] == pytest.approx(
            CP["hepatic_artery_flow"]["value"][k - 1], rel=1e-12
        )
        assert q[Organ.KIDNEY] / q0[Organ.KIDNEY] == pytest.approx(
            CP["kidney_flow"]["value"][k - 1], rel=1e-12
        )
        for o in (*PORTAL_ORGANS, Organ.BRAIN):
            assert q[o] == pytest.approx(q0[o], rel=1e-12)
        assert q[Organ.MUSCLE] / q0[Organ.MUSCLE] == pytest.approx(
            q[Organ.SKIN] / q0[Organ.SKIN], rel=1e-12
        )
        assert q[Organ.MUSCLE] > q0[Organ.MUSCLE]
    for changes in (
        {"HF_severity": 2.0},
        {"HF_severity": 4.0},
        {"f_renal_function": 0.1},
        {"CP_severity": 3.0, "HF_severity": 3.0, "f_renal_function": 0.3},
    ):
        _assert_flows_sum(individual(model.flat_path, **changes))


def _fu(fu_ref: float, factor: float) -> float:
    """The fraction unbound of McNamara 2002 with a binding factor."""
    return 1.0 / (1.0 + factor * (1.0 - fu_ref) / fu_ref)


def test_fu_follows_albumin(tmp_path: Path) -> None:
    """fu follows the factor of its binding protein (McNamara 2002).

    An unknown protein takes no factor, an unbound substance stays unbound.
    The CKD factor divides the binding factor: it multiplies fu of a strongly
    bound substance (Malik 2020) and fu stays below 1.
    """
    substances = [
        Substance(
            sid="alb", name="alb", mr=300.0, binding=Binding(fu=0.1, protein="albumin")
        ),
        Substance(
            sid="agp", name="agp", mr=300.0, binding=Binding(fu=0.2, protein="agp")
        ),
        Substance(sid="unk", name="unk", mr=300.0, binding=Binding(fu=0.3)),
        Substance(
            sid="high",
            name="high",
            mr=300.0,
            binding=Binding(fu=0.9, protein="albumin"),
        ),
        Substance(sid="free", name="free", mr=300.0),
    ]
    d = PBPKDefinition(sid="binding", substances=substances)
    path = build(d, tmp_path).flat_path
    r = individual(path, CP_severity=2.0)
    albumin, agp = CP["albumin"]["value"][1], CP["agp"]["value"][1]
    assert r[ids.fu("alb")] == pytest.approx(_fu(0.1, albumin), rel=1e-12)
    assert r[ids.fu("agp")] == pytest.approx(_fu(0.2, agp), rel=1e-12)
    assert r[ids.fu("unk")] == pytest.approx(0.3, rel=1e-12)
    assert r[ids.fu("free")] == 1.0
    assert r[ids.fu("alb")] > 0.1
    r = individual(path, f_renal_function=0.05)
    factor = CKD["plasma_protein_fu"]["value"][3]
    assert r[ids.fu("alb")] == pytest.approx(_fu(0.1, 1.0 / factor), rel=1e-12)
    assert r[ids.fu("unk")] == pytest.approx(_fu(0.3, 1.0 / factor), rel=1e-12)
    assert r[ids.fu("high")] == pytest.approx(_fu(0.9, 1.0 / factor), rel=1e-12)
    assert r[ids.fu("high")] < 1.0
    assert _fu(0.001, 1.0 / factor) == pytest.approx(0.001 * factor, rel=1e-3)
    assert unbound_fraction(substances[4]) == ids.fu_ref("free")


def _stage_function(r: roadrunner.RoadRunner, egfr: float) -> float:
    """The renal function which gives a GFR per 1.73 m2 (ml/min) of the model."""
    return egfr / (r[ids.EGFR] * 1000.0)


def test_ckd_stage_factors(model: BuildResult) -> None:
    """The renal function gives the CKD stage and its factors.

    Normal >= 90, mild 60-89, moderate 30-59, severe 15-29, end stage < 15
    ml/min/1.73 m2; the kidney volume and flow follow the regressions, the
    body weight is kept (adipose takes the kidney mass), the hepatic factors
    act with `CKD_hepatic` only.
    """
    healthy = individual(model.flat_path)
    assert healthy[ids.EGFR] * 1000.0 > 90.0
    for egfr, stage in ((100.0, None), (75.0, 0), (45.0, 1), (22.0, 2), (10.0, 3)):
        f = _stage_function(healthy, egfr)
        r = individual(model.flat_path, f_renal_function=f)
        assert r[ids.EGFR] * 1000.0 == pytest.approx(egfr, rel=1e-12)
        for key in ("plasma_protein_fu", "gastric_emptying_time", "si_transit_time"):
            expected = 1.0 if stage is None else CKD[key]["value"][stage]
            assert r[ids.condition_factor("ckd", key)] == expected, (key, egfr)
        hct = 1.0 if stage is None else CKD["hematocrit"]["value"][stage] / HEALTHY_HCT
        assert r[ids.HCT] == pytest.approx(healthy[ids.HCT] * hct, rel=1e-12)
        x0 = healthy[ids.EGFR] * 1000.0
        a, b, _ = CKD["kidney_volume"]["value"]
        volume = math.exp(a * (egfr**2 - x0**2) + b * (egfr - x0))
        aq, bq, _ = CKD["kidney_flow"]["value"]
        flow = math.exp(aq * (egfr**2 - x0**2) + bq * (egfr - x0)) * volume
        mki = ids.organ_mass(Organ.KIDNEY)
        mad = ids.organ_mass(Organ.ADIPOSE)
        assert r[mki] == pytest.approx(healthy[mki] * volume, rel=1e-12)
        assert r[mad] - healthy[mad] == pytest.approx(healthy[mki] - r[mki], rel=1e-9)
        assert r[ids.condition_factor("ckd", "kidney_flow")] == pytest.approx(
            flow, rel=1e-12
        )
        assert r[ids.GFR] == pytest.approx(healthy[ids.GFR] * f, rel=1e-12)
        _assert_flows_sum(r)
        cyp2d6 = 1.0 if stage is None else CKD["proteins"]["CYP2D6"]["value"][stage]
        liver_2d6 = ids.impairment_activity("CYP2D6", Organ.LIVER)
        assert r[liver_2d6] == 1.0
        on = individual(model.flat_path, f_renal_function=f, CKD_hepatic=1.0)
        assert on[liver_2d6] == pytest.approx(cyp2d6, rel=1e-12)


def test_heart_failure_redistribution(model: BuildResult) -> None:
    """NYHA classes reduce the flows of liver, portal organs, kidney and periphery.

    The cardiac output is the sum of the flows, heart and brain keep theirs,
    the GFR follows the renal blood flow.
    """
    healthy = individual(model.flat_path)
    q0 = _flows(healthy)
    for k in (1, 2, 3, 4):
        r = individual(model.flat_path, HF_severity=float(k))
        q = _flows(r)
        liver = HF["liver_flow"]["value"][k - 1]
        kidney = HF["kidney_flow"]["value"][k - 1]
        peripheral = HF["peripheral_flow"]["value"][k - 1]
        for o in (Organ.LIVER, *PORTAL_ORGANS):
            assert q[o] / q0[o] == pytest.approx(liver, rel=1e-12)
        assert q[Organ.KIDNEY] / q0[Organ.KIDNEY] == pytest.approx(kidney, rel=1e-12)
        for o in (Organ.SKIN, Organ.ADIPOSE, Organ.MUSCLE):
            assert q[o] / q0[o] == pytest.approx(peripheral, rel=1e-12)
        for o in (Organ.HEART, Organ.BRAIN, Organ.BONE, Organ.GONADS):
            assert q[o] == pytest.approx(q0[o], rel=1e-12)
        assert r[ids.CO] == pytest.approx(sum(q.values()), rel=1e-12)
        assert r[ids.GFR] / healthy[ids.GFR] == pytest.approx(kidney, rel=1e-12)
    assert individual(model.flat_path, HF_severity=4.0)[ids.CO] < healthy[ids.CO]


def _values(path: Path) -> dict[str, float]:
    """The values of the parameters and compartments of a model."""
    r = roadrunner.RoadRunner(str(path))
    return {
        pid: r[pid]
        for pid in [
            *r.model.getGlobalParameterIds(),
            *r.model.getCompartmentIds(),
        ]
    }


def test_healthy_is_identity(tmp_path: Path, model: BuildResult) -> None:
    """All severities 0 give the healthy model.

    The physiology with the impairment has the values of the physiology
    without it (to the last digits: the products with the factors 1 are
    evaluated in another order), and every factor of the model is 1: fu is
    the reference fraction, the liver has no shunt and is fully functional.
    """
    d = _impaired_definition()
    for sex in Sex:
        paths = {}
        for name, objects in (
            ("healthy", physiology_objects(sex)),
            (
                "impaired",
                [
                    *physiology_objects(sex, physiology_factors()),
                    *impairment_objects(d),
                ],
            ),
        ):
            paths[name] = tmp_path / f"{name}_{sex.value}.xml"
            create_model(
                Model(
                    sid=f"{name}_{sex.value}",
                    units=U,
                    model_units=MODEL_UNITS,
                    objects=objects,
                ),
                filepath=paths[name],
                validate=False,
            )
        healthy, impaired = _values(paths["healthy"]), _values(paths["impaired"])
        assert healthy == pytest.approx({k: impaired[k] for k in healthy}, rel=1e-15)
        assert impaired[ids.HCT] == impaired[ids.ref(ids.HCT)]
    r = load_model(model)
    factors = [
        p
        for p in r.model.getGlobalParameterIds()
        if p.startswith(("fimp_", "fcp_", "fckd_", "fhf_"))
        and p != ids.condition_factor("cp", "shunt")
    ]
    assert factors
    for p in factors:
        assert r[p] == 1.0, p
    assert r[ids.fu("drg")] == r[ids.fu_ref("drg")]
    assert (r[ids.F_SHUNT], r[ids.F_FUNCTIONAL]) == (0.0, 1.0)


def test_enzyme_factor_multiplies_activity(model: BuildResult) -> None:
    """The impairment multiplies the organ-local activity `LI__f_cyp3a4`.

    A genotype or an interaction which sets the local activity acts on top of
    the impairment: the rate law has both factors, and halving the local
    activity under cirrhosis halves the rate.
    """
    import libsbml

    doc = libsbml.readSBMLFromFile(str(model.flat_path))
    law = libsbml.formulaToL3String(
        doc.getModel().getReaction("LI__DRG2M").getKineticLaw().getMath()
    )
    assert "LI__f_cyp3a4" in law
    assert ids.impairment_activity("CYP3A4", Organ.LIVER) in law
    rates = {}
    for f in (1.0, 0.5):
        r = individual(model.flat_path, CP_severity=3.0, LI__f_cyp3a4=f)
        r["[LI__Cti_drg]"] = 0.001
        rates[f] = r["LI__DRG2M"]
        assert r[ids.impairment_activity("CYP3A4", Organ.LIVER)] == pytest.approx(
            CP["proteins"]["CYP3A4"]["value"][2]
        )
    assert rates[0.5] == pytest.approx(0.5 * rates[1.0], rel=1e-12)


def test_gut_transit_follows_ckd(tmp_path: Path) -> None:
    """The stomach and small intestinal transit of the intestine read the CKD factors."""
    r = individual(build(definition(), tmp_path).flat_path, f_renal_function=0.1)
    assert r["GU__T_st"] > 0
    assert (
        r[ids.impairment_factor("t_stomach")]
        == CKD["gastric_emptying_time"]["value"][3]
    )
    assert r[ids.impairment_factor("t_si")] == CKD["si_transit_time"]["value"][3]


@pytest.mark.parametrize("dynamic", [True, False])
def test_dynamic_renal_function(tmp_path: Path, dynamic: bool) -> None:
    """A dynamic renal function changes GFR and urinary excretion during a simulation.

    Without `dynamic` the change of `f_renal_function` in the middle of the
    simulation changes nothing (the GFR is an initial assignment), as
    documented.
    """
    d = PBPKDefinition(
        sid="renal",
        substances=[simple_drug()],
        routes=[Route.iv("drg")],
        dynamic=("f_renal_function",) if dynamic else (),
    )
    r = load_model(build(d, tmp_path))
    r[ids.iv_dose("drg")] = 10.0
    r.simulate(0, 120, 2)
    filtration = ids.flat(Organ.KIDNEY, ids.filtration_reaction("drg"))
    cpl = "[KI__Cpl_drg]"
    gfr, clearance, fu = r[ids.GFR], r[filtration] / r[cpl], r[ids.fu("drg")]
    r[ids.F_RENAL_FUNCTION] = 0.5
    factor = 0.5 if dynamic else 1.0
    assert r[ids.GFR] == pytest.approx(factor * gfr, rel=1e-12)
    # moderate CKD: fu follows the plasma protein factor as well, the
    # volumes the hematocrit and the kidney mass (so the concentration)
    assert (r[ids.fu("drg")] > fu) is dynamic
    assert r[filtration] / r[cpl] == pytest.approx(
        factor * clearance * r[ids.fu("drg")] / fu, rel=1e-12
    )
    r.simulate(120, 180, 2)
    assert r[ids.GFR] == pytest.approx(factor * gfr, rel=1e-12)


def test_dynamic_names_parameters(tmp_path: Path) -> None:
    """A dynamic name must be a parameter of the top model."""
    with pytest.raises(ValueError, match="nope"):
        build(definition(dynamic=("nope",)), tmp_path)
    with pytest.raises(ValueError, match="SBML SId"):
        definition(dynamic=("1x",)).validate()
