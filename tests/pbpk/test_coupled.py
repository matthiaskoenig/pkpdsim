"""Tests of the coupled physiology: the definitions of several packages on one body.

`PBPKDefinition.combine` unites the substances, chemistry, routes, PD
components and dynamic parameters of the definitions on one physiology; an
inhibitor of one package acts on the processes of the other (a drug-drug
interaction), since the organs place the inhibitors of every substance in
the rates of the processes of their protein.
"""

import dataclasses
from pathlib import Path

import numpy as np
import pytest
from sbmlutils.factory import Creator

from pkpdsim.pbpk import (
    MM,
    Binding,
    Inhibition,
    Kp,
    Metabolism,
    Mode,
    Observer,
    Organ,
    Part,
    PBPKDefinition,
    Pool,
    Route,
    Substance,
    U,
    Variants,
    ids,
)
from pkpdsim.pbpk.builder import build, compose
from pkpdsim.pbpk.physiology import Sex

from .conftest import cardiac_pd, load

DAY = 24 * 60.0
"""One day (min)."""

KI = 1e-3
"""The competitive inhibition constant of the inhibitor on CYP3A4 (mM)."""


def _drg_package() -> PBPKDefinition:
    """The package of the drug: a low-extraction CYP3A4 substrate of the liver.

    Flow-limited (the hepatic clearance is the well-stirred one), no other
    elimination; `CLint_u = Vmax*Vli/Km` about 0.02 l/min with fu 0.5, an
    extraction of about 1 % of the hepatic plasma flow, and a dose far below
    Km (linear kinetics).
    """
    drug = Substance(
        sid="drg",
        name="drug",
        mr=300.0,
        binding=Binding(fu=0.5, protein="albumin"),
        distribution=Kp(1.0, mode=Mode.PERFUSION),
    )
    return PBPKDefinition(
        sid="drg_body",
        substances=[drug],
        chemistry=[
            Metabolism(
                sid="DRG3A4",
                organ=Organ.LIVER,
                substrate="drg",
                products={},
                kinetics=MM(vmax=0.011, km=1.0),
                protein="CYP3A4",
            )
        ],
        routes=[Route.iv("drg")],
    )


def _inh_package(ki: float = KI) -> PBPKDefinition:
    """The package of the inhibitor: competitive on CYP3A4, not eliminated.

    Without elimination its amount stays the dose, so after the
    distribution its unbound concentration is constant.
    """
    inhibitor = Substance(
        sid="inh",
        name="inhibitor",
        mr=400.0,
        distribution=Kp(1.0, mode=Mode.PERFUSION),
        inhibits=(Inhibition(protein="CYP3A4", kind="competitive", ki=ki),),
    )
    return PBPKDefinition(
        sid="inh_body", substances=[inhibitor], routes=[Route.iv("inh")]
    )


def test_combine_two_definitions() -> None:
    """The union of two packages on one physiology, in the order of the packages."""
    drg, inh = _drg_package(), _inh_package()
    pd = cardiac_pd()
    drg = dataclasses.replace(drg, pd=[pd], dynamic=(ids.F_RENAL_FUNCTION,))
    inh = dataclasses.replace(
        inh,
        dynamic=(ids.F_RENAL_FUNCTION,),
        variants=Variants({Organ.KIDNEY: "nephron"}),
    )
    d = PBPKDefinition.combine(drg, inh, sid="coupled")
    assert d.sid == "coupled"
    assert [s.sid for s in d.substances] == ["drg", "inh"]
    assert [p.sid for p in d.chemistry] == ["DRG3A4"]
    assert [(r.kind, r.substance) for r in d.routes] == [("iv", "drg"), ("iv", "inh")]
    assert list(d.pd) == [pd]
    assert d.dynamic == (ids.F_RENAL_FUNCTION,)
    assert dict(d.variants.organs) == {Organ.KIDNEY: "nephron"}
    assert d.sex is Sex.MALE
    assert d.population == drg.population
    d.validate()
    model = compose(d)
    rate = next(
        r
        for md in model.model_definitions or []
        if md.sid.endswith("_li")
        for r in md.reactions or []
        if str(r.sid) == "DRG3A4"
    )
    # the inhibitor of the other package is in the rate of the CYP3A4 process
    assert rate.formula is not None
    assert ids.ki("CYP3A4", "inh") in str(rate.formula.math)


def test_combine_shared_object_accepted() -> None:
    """A substance which two packages define as the same object is taken once."""
    shared = Substance(sid="shared_m", name="shared metabolite", mr=200.0)
    a = dataclasses.replace(
        _drg_package(), substances=[*_drg_package().substances, shared]
    )
    b = dataclasses.replace(
        _inh_package(), substances=[*_inh_package().substances, shared]
    )
    d = PBPKDefinition.combine(a, b, sid="coupled")
    assert [s.sid for s in d.substances] == ["drg", "shared_m", "inh"]


@pytest.mark.parametrize(
    ("change", "match"),
    [
        (
            {"substances": [Substance(sid="drg", name="other drug", mr=1.0)]},
            "substance 'drg'",
        ),
        (
            {
                "chemistry": [
                    dataclasses.replace(
                        _drg_package().chemistry[0], kinetics=MM(vmax=1.0, km=1.0)
                    )
                ]
            },
            "process 'DRG3A4'",
        ),
        (
            {"pd": [cardiac_pd(reads="inh")]},
            "PD component 'cardio' has different inputs or outputs",
        ),
        ({"variants": Variants({Organ.LIVER: "zonated"})}, "organ 'li'"),
        ({"sex": Sex.FEMALE}, "sex"),
        ({"population": "east_asian_tanaka"}, "population"),
    ],
    ids=["substance", "process", "pd", "variant", "sex", "population"],
)
def test_combine_duplicate_ids_raises(change: dict[str, object], match: str) -> None:
    """The same id defined differently by two packages raises, naming it.

    So do different variants of an organ and a different sex or population:
    a coupled model has one physiology.
    """
    base = dataclasses.replace(
        _drg_package(), pd=[cardiac_pd()], variants=Variants({Organ.LIVER: "standard"})
    )
    other = dataclasses.replace(
        _inh_package(), pd=[cardiac_pd()], variants=Variants({Organ.LIVER: "standard"})
    )
    other = dataclasses.replace(
        other,
        **{
            k: (
                [*other.substances, *v]
                if k == "substances" and isinstance(v, list)
                else v
            )
            for k, v in change.items()
        },
    )
    with pytest.raises(ValueError, match=match):
        PBPKDefinition.combine(base, other, sid="coupled")


def test_combine_pd_equal_by_model_id_and_ports() -> None:
    """The same PD component, constructed by each package for itself, is taken once.

    PD components are equal by the id of their model, their inputs and their
    outputs (R50): two separately built `cardiac_pd()` merge, a component of
    the same id with other inputs raises (the "pd" case of
    `test_combine_duplicate_ids_raises`).
    """
    a = dataclasses.replace(_drg_package(), pd=[cardiac_pd()])
    b = dataclasses.replace(_inh_package(), pd=[cardiac_pd()])
    assert a.pd[0].model is not b.pd[0].model
    d = PBPKDefinition.combine(a, b, sid="coupled")
    assert len(d.pd) == 1
    assert d.pd[0] == cardiac_pd()
    different = dataclasses.replace(b, pd=[cardiac_pd(reads="inh")])
    with pytest.raises(ValueError, match="different inputs or outputs"):
        PBPKDefinition.combine(a, different, sid="coupled")


def test_combine_omitted_variant_is_no_requirement() -> None:
    """An organ which a package leaves out takes the variant of the other package.

    The inhibitor package does not choose a liver: the coupled model has the
    zonated liver of the drug package, in either order.
    """
    drg = dataclasses.replace(
        _drg_package(), variants=Variants({Organ.LIVER: "zonated"})
    )
    inh = _inh_package()
    assert Organ.LIVER not in inh.variants.organs
    for order in ((drg, inh), (inh, drg)):
        d = PBPKDefinition.combine(*order, sid="coupled")
        assert dict(d.variants.organs) == {Organ.LIVER: "zonated"}
    standard = dataclasses.replace(inh, variants=Variants({Organ.LIVER: "standard"}))
    with pytest.raises(ValueError, match="organ 'li'"):
        PBPKDefinition.combine(drg, standard, sid="coupled")


def test_combine_nothing_raises() -> None:
    """Combining no definition raises."""
    with pytest.raises(ValueError, match="no definitions"):
        PBPKDefinition.combine(sid="coupled")


def _auc(
    r: object, sid: str, start: float, end: float, dose: float
) -> tuple[float, float, np.ndarray]:
    """Dose the drug at `start` and integrate its venous plasma concentration.

    A fine grid over the distribution (10 h), a coarser one after it; the
    trapezoidal rule, and the terminal tail `C_end/lambda_z` from the last
    two points. Returns the AUC to infinity (mM*min), the share of the tail
    in it and the unbound concentration of the inhibitor in the liver over
    the time.
    """
    import roadrunner

    assert isinstance(r, roadrunner.RoadRunner)
    c = f"[{ids.pool_species(Pool.VENOUS, sid)}]"
    iu = ids.flat(Organ.LIVER, ids.unbound(Part.ORGAN, "inh"))
    r[ids.iv_dose(sid)] = dose
    early = r.simulate(start, start + 600, 6001, selections=["time", c, iu])
    late = r.simulate(start + 600, end, 4001, selections=["time", c, iu])
    auc = float(np.trapezoid(early[c], early["time"]))
    auc += float(np.trapezoid(late[c], late["time"]))
    t, y = late["time"], late[c]
    lambda_z = float(np.log(y[-2] / y[-1]) / (t[-1] - t[-2]))
    tail = float(y[-1]) / lambda_z
    return auc + tail, tail / (auc + tail), np.concatenate([early[iu], late[iu]])


def test_cross_package_inhibition(tmp_path: Path) -> None:
    """The inhibitor of one package raises the AUC of the drug of the other by `1 + I_u/Ki`.

    The drug is a low-extraction CYP3A4 substrate of the flow-limited liver
    without other elimination, dosed intravenously far below Km (linear);
    the inhibitor is competitive on CYP3A4, dosed three days before the drug
    and not eliminated, so its unbound concentration in the liver `I_u` is
    constant over the AUC of the drug (checked to 1e-6). The intrinsic
    clearance falls by `1 + I_u/Ki` and the AUC rises by it.

    The expectation is the well-stirred liver of the model, not the
    intrinsic-clearance limit: the flow-limited liver is perfused by the
    plasma flow `Qp = Qh*(1 - HCT)` (the drug does not enter the
    erythrocytes, no shunt) and metabolizes `CLint_u*fu*Cout/(1 + r)` with
    `CLint_u = Vmax*Vli/Km`, so the plasma clearance is
    `Qp*fu*CLint_u/(1 + r)/(Qp + fu*CLint_u/(1 + r))` and the AUC ratio
    `(Qp*(1 + r) + fu*CLint_u)/(Qp + fu*CLint_u)`, about 0.6 % below
    `1 + r` at an extraction of about 1 %. No other organ eliminates, so the
    AUC is that of any pool. Tolerance 1e-3: the saturation of the drug
    (`C/Km` < 1e-4), the trapezoidal rule and the extrapolated tail
    (< 1e-4 of the AUC) are below it.
    """
    d = PBPKDefinition.combine(_drg_package(), _inh_package(), sid="ddi")
    result = build(d, tmp_path)
    t_dose, end = 3 * DAY, 3 * DAY + 60 * DAY

    control = load(result)
    control.simulate(0, t_dose, 11)
    auc_control, tail_control, _ = _auc(control, "drg", t_dose, end, dose=1.0)

    inhibited = load(result)
    inhibited[ids.iv_dose("inh")] = 24.0
    inhibited.simulate(0, t_dose, 101)
    auc_inhibited, tail_inhibited, iu = _auc(inhibited, "drg", t_dose, end, dose=1.0)

    assert iu.max() == pytest.approx(iu.min(), rel=1e-6)
    ratio = 1 + float(iu.mean()) / KI
    assert 1.5 < ratio < 3.0
    # the extrapolated terminal tail is a small part of the AUC
    assert tail_control < 1e-3
    assert tail_inhibited < 1e-2
    r = ratio - 1
    qp = float(control[ids.QH]) * (1 - float(control[ids.HCT]))
    clint = (
        float(control[ids.flat(Organ.LIVER, "DRG3A4_Vmax")])
        * float(control[ids.organ_volume(Organ.LIVER)])
        / float(control[ids.flat(Organ.LIVER, "DRG3A4_Km")])
        * float(control[ids.fu("drg")])
    )
    assert float(control[ids.F_SHUNT]) == 0.0
    assert float(control[ids.F_FUNCTIONAL]) == 1.0
    assert float(control[ids.flat(Organ.LIVER, ids.activity("CYP3A4"))]) == 1.0
    expected = (qp * (1 + r) + clint) / (qp + clint)
    assert 0.003 < 1 - expected / ratio < 0.02
    assert auc_inhibited / auc_control == pytest.approx(expected, rel=1e-3)


def _creator(given: str = "Matthias") -> Creator:
    """A creator of a package (a placeholder email, no personal data)."""
    return Creator(
        familyName="König",
        givenName=given,
        email="creator@example.org",
        organization="Humboldt-University Berlin",
    )


def test_combine_observers_and_creators() -> None:
    """The observers and the creators of the packages are united.

    An observer which two packages define equally is taken once, one which
    they define differently raises; a creator is the same person when all its
    fields are equal (a sbmlutils `Creator` compares by identity).
    """
    total = Observer("Carm_total", "Carm_drg + Carm_inh", U.mM, "drug and inhibitor")
    a = dataclasses.replace(
        _drg_package(),
        observers=[total, Observer("Carm_drg_uM", "Carm_drg", U.mM, "drug")],
        creators=[_creator()],
    )
    b = dataclasses.replace(
        _inh_package(),
        observers=[
            Observer("Carm_total", "Carm_drg + Carm_inh", U.mM, "drug and inhibitor")
        ],
        creators=[_creator(), _creator("Michelle")],
    )
    d = PBPKDefinition.combine(a, b, sid="coupled")
    assert [o.sid for o in d.observers] == ["Carm_total", "Carm_drg_uM"]
    assert [c.givenName for c in d.creators] == ["Matthias", "Michelle"]
    d.validate()
    other = dataclasses.replace(
        b, observers=[Observer("Carm_total", "Carm_inh", U.mM, "inhibitor")]
    )
    with pytest.raises(ValueError, match="observer 'Carm_total'"):
        PBPKDefinition.combine(a, other, sid="coupled")
