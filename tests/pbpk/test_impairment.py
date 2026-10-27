"""Tests of the impaired organ function: Child-Pugh, CKD and heart failure."""

import dataclasses
import itertools
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
            "hematocrit",
        ):
            assert r[ids.condition_factor("cp", key)] == pytest.approx(
                CP[key]["value"][k - 1], rel=1e-12
            ), key
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


@pytest.mark.parametrize("fu_ref", [0.001, 0.05, 0.3, 0.9])
def test_fu_follows_albumin(tmp_path: Path, fu_ref: float) -> None:
    """fu follows the factor of its binding protein (McNamara 2002).

    An unknown protein takes no factor, an unbound substance stays unbound.
    The CKD factor is a factor of the binding capacity: it divides the
    binding factor, so fu of a strongly bound substance rises by it, and fu
    stays below 1.
    """
    substances = [
        Substance(
            sid="alb",
            name="alb",
            mr=300.0,
            binding=Binding(fu=fu_ref, protein="albumin"),
        ),
        Substance(
            sid="agp", name="agp", mr=300.0, binding=Binding(fu=fu_ref, protein="agp")
        ),
        Substance(sid="unk", name="unk", mr=300.0, binding=Binding(fu=fu_ref)),
        Substance(sid="free", name="free", mr=300.0),
    ]
    d = PBPKDefinition(sid="binding", substances=substances)
    path = build(d, tmp_path).flat_path
    healthy = individual(path)
    for sid in ("alb", "agp", "unk"):
        assert healthy[ids.fu(sid)] == pytest.approx(fu_ref, rel=1e-14)
    r = individual(path, CP_severity=2.0)
    albumin, agp = CP["albumin"]["value"][1], CP["agp"]["value"][1]
    assert r[ids.fu("alb")] == pytest.approx(_fu(fu_ref, albumin), rel=1e-12)
    assert r[ids.fu("agp")] == pytest.approx(_fu(fu_ref, agp), rel=1e-12)
    assert r[ids.fu("unk")] == pytest.approx(fu_ref, rel=1e-12)
    assert r[ids.fu("free")] == 1.0
    assert r[ids.fu("alb")] > fu_ref
    r = individual(path, f_renal_function=0.05)
    factor = CKD["plasma_protein"]["value"][3]
    for sid in ("alb", "unk"):
        assert r[ids.fu(sid)] == pytest.approx(_fu(fu_ref, 1.0 / factor), rel=1e-12)
        assert fu_ref < r[ids.fu(sid)] < 1.0
    if fu_ref <= 0.001:
        assert r[ids.fu("alb")] == pytest.approx(fu_ref * factor, rel=2e-3)
    assert unbound_fraction(substances[3]) == ids.fu_ref("free")


def test_an_unbound_fraction_without_factor_is_the_reference() -> None:
    """No protein and no CKD factor leaves no factor: fu is `fu_ref` itself."""
    s = Substance(sid="x", name="x", mr=300.0, binding=Binding(0.4, ckd=False))
    formula = unbound_fraction(s)
    assert formula == ids.fu_ref("x")
    assert "dimensionless" not in formula


def test_binding_free_of_ckd(tmp_path: Path) -> None:
    """A binding with `ckd=False` ignores the CKD factor of the binding capacity.

    fu stays fu_ref in renal impairment, while the factor of the binding
    protein still acts in cirrhosis; the default builds the formula as before.
    """
    fu_ref = 0.07
    substances = [
        Substance(
            sid="keep",
            name="keep",
            mr=300.0,
            binding=Binding(fu=fu_ref, protein="albumin", ckd=False),
        ),
        Substance(
            sid="rise",
            name="rise",
            mr=300.0,
            binding=Binding(fu=fu_ref, protein="albumin"),
        ),
    ]
    assert Binding(fu=fu_ref).ckd is True
    d = PBPKDefinition(sid="binding", substances=substances)
    path = build(d, tmp_path).flat_path
    r = individual(path, f_renal_function=0.2)
    assert r[ids.fu("keep")] == pytest.approx(fu_ref, rel=1e-12)
    assert r[ids.fu("rise")] > fu_ref * 1.2
    r = individual(path, CP_severity=2.0)
    albumin = CP["albumin"]["value"][1]
    assert r[ids.fu("keep")] == pytest.approx(_fu(fu_ref, albumin), rel=1e-12)
    assert r[ids.fu("keep")] > fu_ref
    r = individual(path, CP_severity=2.0, f_renal_function=0.2)
    assert r[ids.fu("keep")] == pytest.approx(_fu(fu_ref, albumin), rel=1e-12)
    keep, rise = (unbound_fraction(s) for s in substances)
    assert "fimp_plasma_protein" not in keep
    assert rise.replace(" ", "") == (
        "1dimensionless/(1dimensionless+fimp_albumin/fimp_plasma_protein*"
        f"(1dimensionless-{ids.fu_ref('rise')})/{ids.fu_ref('rise')})"
    )


def _stage_function(r: roadrunner.RoadRunner, egfr: float) -> float:
    """The renal function which gives a GFR per 1.73 m2 (ml/min) of the model."""
    return egfr / (r[ids.EGFR] * 1000.0)


ANCHORS = CKD["stage_anchor"]["value"]
"""The eGFR (ml/min/1.73 m2) of the values of mild, moderate, severe and end stage."""
LIMIT = CKD["stage_gfr"]["value"][0]
"""The eGFR above which the factors are 1."""
PROTEIN_KEYS = (
    "plasma_protein",
    "hematocrit",
    "kidney_volume",
    "kidney_flow",
    "gastric_emptying_time",
    "si_transit_time",
)


def _expected(values: list[float], egfr: float) -> float:
    """The factor at an eGFR: piecewise linear through 1 at 90 and the anchors.

    Constant (1) above 90 and (end stage value) below the last anchor.
    """
    xs = [LIMIT, *ANCHORS]
    ys = [1.0, *values]
    if egfr >= xs[0]:
        return ys[0]
    if egfr <= xs[-1]:
        return ys[-1]
    for i in range(len(xs) - 1):
        if xs[i + 1] <= egfr <= xs[i]:
            t = (xs[i] - egfr) / (xs[i] - xs[i + 1])
            return ys[i] + t * (ys[i + 1] - ys[i])
    raise AssertionError(egfr)  # pragma: no cover


def test_ckd_anchors_are_the_stage_midpoints() -> None:
    """The anchors are the midpoints of the stages, the end stage one at 7.5."""
    assert ANCHORS == [75.0, 45.0, 22.5, 7.5]
    assert CKD["stage_anchor"]["source"] == "ModelChoice"
    assert "continuous covariate" in CKD["stage_anchor"]["notes"]


def test_ckd_factors_at_the_anchors(model: BuildResult) -> None:
    """At each anchor eGFR the factor is the value of the stage, 1 from 90."""
    healthy = individual(model.flat_path)
    assert healthy[ids.EGFR] * 1000.0 > LIMIT
    for egfr, stage in (
        (100.0, None),
        (LIMIT, None),
        *zip(ANCHORS, range(4), strict=True),
    ):
        f = _stage_function(healthy, egfr)
        r = individual(model.flat_path, f_renal_function=f)
        assert r[ids.EGFR] * 1000.0 == pytest.approx(egfr, rel=1e-12)
        for key in PROTEIN_KEYS:
            expected = 1.0 if stage is None else CKD[key]["value"][stage]
            assert r[ids.condition_factor("ckd", key)] == pytest.approx(
                expected, rel=1e-12
            ), (key, egfr)
        assert r[ids.HCT] == healthy[ids.HCT]
        mki = ids.organ_mass(Organ.KIDNEY)
        assert r[mki] == healthy[mki]
        assert r[ids.GFR] == pytest.approx(healthy[ids.GFR] * f, rel=1e-12)
        _assert_flows_sum(r)
        cyp2d6 = 1.0 if stage is None else CKD["proteins"]["CYP2D6"]["value"][stage]
        liver_2d6 = ids.impairment_activity("CYP2D6", Organ.LIVER)
        assert r[liver_2d6] == 1.0
        on = individual(model.flat_path, f_renal_function=f, CKD_hepatic=1.0)
        assert on[liver_2d6] == pytest.approx(cyp2d6, rel=1e-12)
    # constant below the end stage anchor and above 90
    for egfr in (0.0, 1.0, 5.0):
        r = individual(model.flat_path, f_renal_function=_stage_function(healthy, egfr))
        assert r[ids.condition_factor("ckd", "plasma_protein")] == pytest.approx(
            CKD["plasma_protein"]["value"][3], rel=1e-12
        )


def test_ckd_factors_interpolate_between_anchors(model: BuildResult) -> None:
    """Between two anchors (and from 90 to the mild one) the factor is linear."""
    healthy = individual(model.flat_path)
    liver_2d6 = ids.impairment_activity("CYP2D6", Organ.LIVER)
    for egfr in (89.0, 82.5, 70.0, 60.0, 59.9, 45.1, 33.0, 15.0, 10.0, 7.6):
        f = _stage_function(healthy, egfr)
        r = individual(model.flat_path, f_renal_function=f, CKD_hepatic=1.0)
        for key in ("plasma_protein", "hematocrit"):
            assert r[ids.condition_factor("ckd", key)] == pytest.approx(
                _expected(CKD[key]["value"], egfr), rel=1e-12
            ), (key, egfr)
        assert r[liver_2d6] == pytest.approx(
            _expected(CKD["proteins"]["CYP2D6"]["value"], egfr), rel=1e-12
        )
    # no step at the former limit of moderate CKD (1 below, 1.45 above)
    below = individual(model.flat_path, f_renal_function=_stage_function(healthy, 59.9))
    above = individual(model.flat_path, f_renal_function=_stage_function(healthy, 60.1))
    pp = ids.condition_factor("ckd", "plasma_protein")
    assert above[pp] == pytest.approx(_expected(CKD["plasma_protein"]["value"], 60.1))
    assert below[pp] == pytest.approx(_expected(CKD["plasma_protein"]["value"], 59.9))
    assert abs(above[pp] - below[pp]) < 5e-3


def test_ckd_fu_scan_is_continuous(tmp_path: Path) -> None:
    """A fine scan of `f_renal_function` from 0 to 1.2 gives a continuous fu.

    fu = c/(c + k) of the capacity factor c (k = (1 - fu_ref)/fu_ref), so the
    jump between neighbouring points is bounded by the steepest slope of the
    factor between the anchors times the step of the eGFR times the largest
    slope of fu in c, `k/(c + k)^2` at c = 1. The former staircase jumped by
    0.021 at the limit of moderate CKD, 20 times more.
    """
    fu_ref = 0.05
    k = (1.0 - fu_ref) / fu_ref
    substance = Substance(
        sid="alb", name="alb", mr=300.0, binding=Binding(fu=fu_ref, protein="albumin")
    )
    d = PBPKDefinition(sid="fuscan", substances=[substance])
    path = build(d, tmp_path).flat_path
    n = 240
    steps = [i * 1.2 / n for i in range(n + 1)]
    egfr_step = 0.0
    fus = []
    for f in steps:
        r = individual(path, f_renal_function=f)
        fus.append(r[ids.fu("alb")])
        egfr_step = r[ids.EGFR] * 1000.0 / f if f else egfr_step
    egfr_step *= 1.2 / n
    xs = [LIMIT, *ANCHORS]
    ys = [1.0, *CKD["plasma_protein"]["value"]]
    slope = max(abs(ys[i + 1] - ys[i]) / (xs[i] - xs[i + 1]) for i in range(4))
    bound = slope * egfr_step * k / (1.0 + k) ** 2
    jumps = [abs(b - a) for a, b in itertools.pairwise(fus)]
    assert max(jumps) <= bound * (1.0 + 1e-9)
    assert max(jumps) < 0.001
    # no step: the scan spans the factors 1 to the end stage value
    assert fus[0] == pytest.approx(_fu(fu_ref, 1.0 / ys[-1]), rel=1e-12)
    assert fus[-1] == pytest.approx(fu_ref, rel=1e-12)


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
    the reference fraction (to the last digit, McNamara with the factor 1),
    the liver has no shunt and is fully functional.
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
    assert r[ids.fu("drg")] == pytest.approx(r[ids.fu_ref("drg")], rel=1e-14)
    assert (r[ids.F_SHUNT], r[ids.F_FUNCTIONAL]) == (0.0, 1.0)


def _without_impairment(path: Path, out: Path) -> Path:
    """The flat model with the impairment removed.

    Every factor of impairment is a constant 1 (the shunt 0), `f_shunt` and
    `f_functional` keep their healthy values and fu is the reference fraction.
    """
    import libsbml

    doc = libsbml.readSBMLFromFile(str(path))
    model = doc.getModel()
    for ia in list(model.getListOfInitialAssignments()):
        symbol = ia.getSymbol()
        if symbol.startswith(("fimp_", "fcp_", "fckd_", "fhf_")):
            value = 0.0 if symbol == ids.condition_factor("cp", "shunt") else 1.0
            model.getParameter(symbol).setValue(value)
            model.removeInitialAssignment(symbol)
        elif symbol in (ids.F_SHUNT, ids.F_FUNCTIONAL):
            model.removeInitialAssignment(symbol)
        elif symbol == ids.fu("drg"):
            ia.setMath(libsbml.parseL3Formula(ids.fu_ref("drg")))
    libsbml.writeSBMLToFile(doc, str(out))
    return out


def test_healthy_model_output_is_identity(tmp_path: Path, model: BuildResult) -> None:
    """The arm-vein concentration of an IV dose with all severities 0 is that without impairment."""
    curves = []
    for path in (
        model.flat_path,
        _without_impairment(model.flat_path, tmp_path / "x.xml"),
    ):
        r = load_model(BuildResult(comp_path=path, flat_path=path))
        r[ids.iv_dose("drg")] = 10.0
        result = r.simulate(0, 600, 61, selections=["time", ids.arm("drg")])
        curves.append(result[:, 1])
    assert max(curves[0]) > 0
    assert list(curves[0]) == pytest.approx(list(curves[1]), rel=1e-10, abs=1e-15)


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


@pytest.mark.parametrize("cp_severity", [0.0, 2.0])
def test_gfr_base_is_the_gfr_at_full_renal_function(
    model: BuildResult, cp_severity: float
) -> None:
    """`GFR_base` is the GFR at renal function 1, the renal function its factor.

    It does not change with `f_renal_function` (defined without renal
    function, where the GFR is 0) and carries the other impairments
    (`fimp_GFR`) and the body surface area.
    """
    bases = []
    for f in (1.0, 0.3, 0.0):
        r = individual(
            model.flat_path, f_renal_function=f, CP_severity=cp_severity, BW=100.0
        )
        base = r[ids.GFR_BASE]
        assert r[ids.GFR] == pytest.approx(base * f, rel=1e-12, abs=1e-15)
        assert base == pytest.approx(
            r[ids.ref(ids.GFR)]
            * r[ids.BSA]
            / r[ids.ref(ids.BSA)]
            * r[ids.impairment_factor(ids.GFR)],
            rel=1e-12,
        )
        bases.append(base)
    assert bases == pytest.approx([bases[0]] * 3, rel=1e-12)
    assert (r[ids.impairment_factor(ids.GFR)] < 1.0) is (cp_severity > 0)


def test_dynamic_names_parameters(tmp_path: Path) -> None:
    """A dynamic name must be a parameter of the top model."""
    with pytest.raises(ValueError, match="nope"):
        build(definition(dynamic=("nope",)), tmp_path)
    with pytest.raises(ValueError, match="SBML SId"):
        definition(dynamic=("1x",)).validate()


def test_unused_factor_raises() -> None:
    """A factor for a quantity which the physiology does not compute raises."""
    with pytest.raises(ValueError, match="nope"):
        physiology_objects(Sex.MALE, {"nope": "fimp_nope"})


def test_dynamic_must_be_a_parameter(tmp_path: Path) -> None:
    """A compartment or species of the top model is no dynamic parameter."""
    for sid in ("Var_pl", "Cve_drg"):
        with pytest.raises(ValueError, match=sid):
            build(definition(dynamic=(sid,)), tmp_path)


@pytest.mark.parametrize(
    ("severity", "value"), [("CP_severity", 3.0), ("HF_severity", 4.0)]
)
def test_dynamic_severity(tmp_path: Path, severity: str, value: float) -> None:
    """A dynamic severity changes the flows and the liver during a simulation."""
    d = definition(routes=[Route.iv("drg")], dynamic=(severity,))
    r = load_model(build(d, tmp_path))
    r[ids.iv_dose("drg")] = 10.0
    r.simulate(0, 60, 2)
    co, qli, shunt = r[ids.CO], r[ids.QHA], r[ids.F_SHUNT]
    r[severity] = value
    _assert_flows_sum(r)
    if severity == "CP_severity":
        assert r[ids.CO] == pytest.approx(
            co * CP["cardiac_index"]["value"][2], rel=1e-12
        )
        assert r[ids.F_SHUNT] == pytest.approx(CP["shunt"]["value"][2], rel=1e-14)
        assert shunt == 0.0
    else:
        assert r[ids.CO] < co
        assert r[ids.QHA] == pytest.approx(
            qli * HF["liver_flow"]["value"][3], rel=1e-12
        )
    r.simulate(60, 120, 2)
    _assert_flows_sum(r)


def test_model_choices_are_cited_as_such(model: BuildResult) -> None:
    """A vector with model-choice classes cites ModelChoice in the SBML notes."""
    import libsbml

    doc = libsbml.readSBMLFromFile(str(model.flat_path))
    notes = (
        doc.getModel()
        .getParameter(ids.condition_factor("cp", "CYP3A4"))
        .getNotesString()
    )
    assert "El-Khateeb" in notes
    assert "A: a model choice" in notes
    assert "Choice of the pkpdsim PBPK library" in notes
    for key in ("albumin", "shunt"):
        notes = (
            doc.getModel()
            .getParameter(ids.condition_factor("cp", key))
            .getNotesString()
        )
        assert "not verified" in notes
        assert "Edginton AN" not in notes
