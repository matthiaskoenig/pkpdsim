"""Tests of the `Formula` values: a value tied to parameters of the top model."""

from pathlib import Path

import libsbml
import pytest
import roadrunner
from sbmlutils.factory import Uncertainty, UncertParameter

from pkpdsim.pbpk import (
    Binding,
    Direction,
    Excretion,
    FirstOrder,
    Formula,
    Kp,
    Membrane,
    Organ,
    PBPKDefinition,
    Route,
    Substance,
    Transport,
    U,
    Value,
    build,
    ids,
)
from pkpdsim.pbpk.placement import formula_ids
from pkpdsim.pbpk.values import as_value

#: a fitted Kp of the parent with its uncertainty
FITTED_KP = Value(
    0.8,
    U.dimensionless,
    "fit",
    uncertainties=[
        Uncertainty(
            uncertParameters=[
                UncertParameter(
                    type=libsbml.DISTRIB_UNCERTTYPE_DISTRIBUTION,
                    math="normal(0.8, 0.1)",
                )
            ]
        )
    ],
)


def _definition(
    k: object = 0.01, kp_m: object = 1.0, dynamic: tuple[str, ...] = ()
) -> PBPKDefinition:
    """A parent and a metabolite whose biliary excretion and Kp may be formulas."""
    parent = Substance(
        "p",
        "parent",
        300.0,
        binding=Binding(fu=0.1),
        distribution=Kp(FITTED_KP, ftissue=1.0),
    )
    metabolite = Substance(
        "pm",
        "metabolite",
        310.0,
        binding=Binding(fu=0.5),
        distribution=Kp(kp_m, ftissue=1.0),  # ty: ignore[invalid-argument-type]
        excretion=Excretion(urine=True),
    )
    bile = Transport(
        "PMBIL",
        Organ.LIVER,
        "pm",
        Direction.EFFLUX,
        Membrane.APICAL,
        FirstOrder(k),  # ty: ignore[invalid-argument-type]
    )
    return PBPKDefinition(
        sid="formula",
        substances=[parent, metabolite],
        chemistry=[bile],
        routes=[Route.iv("p"), Route.iv("pm")],
        dynamic=dynamic,
    )


def _model(path: Path) -> libsbml.Model:
    """The flat model of a build."""
    return libsbml.readSBMLFromFile(str(path)).getModel()


def _formula(model: libsbml.Model, sid: str) -> str | None:
    """The formula of the initial assignment of an id, `None` without one."""
    ia = model.getInitialAssignmentBySymbol(sid)
    return None if ia is None else libsbml.formulaToL3String(ia.getMath())


def test_formula_is_a_value() -> None:
    """A formula keeps expression, unit and notes; `as_value` returns it as it is."""
    f = Formula("Kp_p", U.dimensionless, notes="The Kp of the parent.")
    assert f.expression == "Kp_p"
    assert f.unit == U.dimensionless
    assert f.notes == "The Kp of the parent."
    assert as_value(f, U.mM) is f
    assert Formula("Kp_p", U.dimensionless).notes == ""


def test_an_organ_value_follows_a_parameter_of_the_top_model(tmp_path: Path) -> None:
    """A formula of an organ reads ids of the top model, before the simulation.

    The clearance per volume of the biliary excretion in the liver is tied to
    the body weight; a change of `BW` before the simulation (an initial value,
    as `preinit_changes`) changes it.
    """
    k = Formula(
        "0.0001 per_min * BW / 1 kg", U.per_min, notes="Tied to the body weight."
    )
    flat = build(_definition(k=k), tmp_path).flat_path
    model = _model(flat)
    # the organ reads BW of the top model (the units are those of the organ)
    ia = model.getInitialAssignmentBySymbol("LI__PMBIL_k")
    assert formula_ids(libsbml.formulaToL3String(ia.getMath())) == {"BW"}
    notes = model.getParameter("LI__PMBIL_k").getNotesString()
    assert "Tied to the body weight." in notes
    assert "0.0001 per_min * BW / 1 kg" in notes
    r = roadrunner.RoadRunner(str(flat))
    assert r["LI__PMBIL_k"] == pytest.approx(0.0001 * r["BW"])
    r.setValue("init(BW)", 100.0)
    r.reset()
    assert r["LI__PMBIL_k"] == pytest.approx(0.01)


def test_a_dynamic_id_keeps_the_formula_dynamic(tmp_path: Path) -> None:
    """A formula which reads a dynamic id becomes an assignment rule (R21)."""
    k = Formula("0.0001 per_min * BW / 1 kg", U.per_min)
    model = _model(build(_definition(k=k, dynamic=("BW",)), tmp_path).flat_path)
    assert model.getAssignmentRuleByVariable("LI__PMBIL_k") is not None
    assert model.getInitialAssignmentBySymbol("LI__PMBIL_k") is None


def test_a_fitted_value_tied_by_a_formula_appears_once(tmp_path: Path) -> None:
    """The metabolite takes the Kp of the parent: one fitted parameter, one uncertainty."""
    kp_m = Formula("Kp_p", U.dimensionless, notes="The Kp of the parent.")
    flat = build(_definition(kp_m=kp_m), tmp_path).flat_path
    model = _model(flat)
    assert _formula(model, ids.kp("pm")) == "Kp_p"
    assert model.getInitialAssignmentBySymbol(ids.kp("p")) is None
    uncertain = [
        p.getId()
        for p in model.getListOfParameters()
        if p.getPlugin("distrib") is not None
        and p.getPlugin("distrib").getNumUncertainties()
    ]
    assert uncertain == [ids.kp("p")]
    r = roadrunner.RoadRunner(str(flat))
    assert r[ids.kp("pm")] == pytest.approx(0.8)
    r.setValue(f"init({ids.kp('p')})", 0.5)
    r.reset()
    assert r[ids.kp("pm")] == pytest.approx(0.5)


@pytest.mark.parametrize(
    ("k", "kp_m", "missing"),
    [
        (0.01, Formula("Kp_unknown", U.dimensionless), "Kp_unknown"),
        (Formula("0.1 per_min * f_unknown", U.per_min), 1.0, "f_unknown"),
    ],
    ids=["top", "organ"],
)
def test_an_unknown_id_raises(
    k: object, kp_m: object, missing: str, tmp_path: Path
) -> None:
    """Every id of a formula is an id of the top model, the build names the others."""
    with pytest.raises(ValueError, match=rf"unknown in the top model: \['{missing}'\]"):
        build(_definition(k=k, kp_m=kp_m), tmp_path)


def test_a_formula_is_no_number_for_a_kp_method() -> None:
    """The Kp method computes with fu: a formula there is a problem of the definition."""
    from pkpdsim.pbpk import KpMethod

    s = Substance(
        "q",
        "drug",
        300.0,
        binding=Binding(fu=Formula("fu_ref_p", U.dimensionless)),
        distribution=KpMethod(
            method="rodgers_rowland",
            logp=2.0,
            pka=(),
            compound="neutral",
            ftissue=1.0,
        ),
    )
    d = PBPKDefinition(sid="m", substances=[s], routes=[Route.iv("q")])
    with pytest.raises(ValueError, match=r"q: .*Kp method.*number"):
        d.validate()
