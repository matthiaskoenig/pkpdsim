"""Tests of the `Formula` values: a value tied to parameters of the top model."""

import dataclasses
from pathlib import Path

import libsbml
import pytest
import roadrunner
from sbmlutils.factory import Uncertainty, UncertParameter

from pkpdsim.pbpk import (
    MM,
    Binding,
    CustomReaction,
    Direction,
    Excretion,
    FirstOrder,
    Formula,
    Kp,
    Membrane,
    Metabolism,
    ModelParameter,
    Organ,
    Part,
    PBPKDefinition,
    Pool,
    Route,
    SegmentedIntestine,
    Substance,
    Transport,
    U,
    Value,
    Variants,
    ZonatedLiver,
    build,
    ids,
)
from pkpdsim.pbpk.placement import formula_ids
from pkpdsim.pbpk.values import as_value

from .conftest import definition, individual, simple_drug, total_amount

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
    assert Formula("Kp_p", U.dimensionless).notes is None


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


def _liver_k(k: Formula) -> PBPKDefinition:
    """The simple definition with a first order CYP3A4 metabolism of the liver."""
    chemistry = [
        Metabolism(
            sid="DRG2M",
            organ=Organ.LIVER,
            substrate="drg",
            products={"drg_m": 1.0},
            kinetics=FirstOrder(k=k, per_volume=False),
            protein="CYP3A4",
        )
    ]
    return definition(chemistry=chemistry)


def test_an_organ_formula_reads_the_id_of_the_top_model(tmp_path: Path) -> None:
    """An id of a formula is an id of the top model, not a local id of the organ.

    The liver reads the impairment factor of CYP3A4 of the top model,
    `fimp_li_cyp3a4`, under its local id `fimp_cyp3a4`: a formula which reads
    `fimp_cyp3a4` reads no id of the top model and raises, one which reads
    `fimp_li_cyp3a4` reads it.
    """
    with pytest.raises(
        ValueError, match=r"unknown in the top model: \['fimp_cyp3a4'\]"
    ):
        build(_liver_k(Formula("0.01 l_per_min * fimp_cyp3a4", U.l_per_min)), tmp_path)
    flat = build(
        _liver_k(Formula("0.01 l_per_min * fimp_li_cyp3a4", U.l_per_min)),
        tmp_path / "top",
    ).flat_path
    r = roadrunner.RoadRunner(str(flat))
    r.setValue("init(CP_severity)", 2.0)
    r.reset()
    assert r["fimp_li_cyp3a4"] < 1.0
    assert r["LI__DRG2M_k"] == pytest.approx(0.01 * r["fimp_li_cyp3a4"])


@pytest.mark.parametrize("sid", ["Cve_drg", "Carm_drg"], ids=["species", "observer"])
def test_a_formula_reads_no_state_of_the_top_model(sid: str, tmp_path: Path) -> None:
    """A species or an observer changes during the simulation, a formula would freeze it.

    The organs are built before the observers of the top model, so an organ
    refuses an observer as unknown.
    """
    kp = Formula(f"2 dimensionless + {sid} / 1 mM", U.dimensionless)
    drug = dataclasses.replace(simple_drug(), distribution=Kp(kp, ftissue=1.0))
    d = definition(
        substances=[drug, Substance(sid="drg_m", name="metabolite", mr=316.0)]
    )
    with pytest.raises(ValueError, match=rf"change during the simulation: \['{sid}'\]"):
        build(d, tmp_path)
    with pytest.raises(ValueError, match=rf"\['{sid}'\]"):
        build(
            _liver_k(Formula(f"0.01 l_per_min * {sid} / 1 mM", U.l_per_min)),
            tmp_path / "o",
        )


def test_an_organ_formula_reads_a_compartment_of_the_top_model(tmp_path: Path) -> None:
    """The liver reads the urine compartment of the top model through an input port."""
    k = Formula("0.01 l_per_min * Vurine / 1 litre", U.l_per_min)
    result = build(_liver_k(k), tmp_path)
    doc = libsbml.readSBMLFromFile(str(result.comp_path))
    liver = doc.getPlugin("comp").getModelDefinition("standard_li")
    assert liver.getCompartment("Vurine") is not None
    r = roadrunner.RoadRunner(str(result.flat_path))
    r.setValue("init(Vurine)", 2.0)
    r.reset()
    assert r["LI__DRG2M_k"] == pytest.approx(0.02)


def test_a_formula_in_the_zonated_liver(tmp_path: Path) -> None:
    """The Vmax of a zonated liver follows the body weight it is tied to."""
    vmax = Formula("0.001 mmole_per_min_l * BW / 70 kg", U.mmole_per_min_l)
    chemistry = [
        Metabolism(
            sid="DRG2M",
            organ=Organ.LIVER,
            substrate="drg",
            products={"drg_m": 1.0},
            kinetics=MM(vmax=vmax, km=0.01),
            protein="CYP3A4",
        )
    ]
    d = definition(
        chemistry=chemistry,
        variants=Variants(organs={Organ.LIVER: ZonatedLiver(n_zones=3)}),
    )
    r = roadrunner.RoadRunner(str(build(d, tmp_path).flat_path))
    r.setValue("init(BW)", 140.0)
    r.reset()
    assert r["LI__DRG2M_Vmax"] == pytest.approx(0.002)


def test_a_formula_in_the_segmented_intestine(tmp_path: Path) -> None:
    """The absorption and a kinetic constant of the intestinal wall follow the body weight."""
    k = Formula("0.5 l_per_min * BW / 70 kg", U.l_per_min)
    chemistry = [
        Metabolism(
            sid="DRG2M_GU",
            organ=Organ.SMALL_INTESTINE,
            substrate="drg",
            products={"drg_m": 1.0},
            kinetics=FirstOrder(k=k, per_volume=False),
            protein="CYP3A4",
        )
    ]
    ka = Formula("1 per_hr * BW / 70 kg", U.per_hr)
    d = definition(
        chemistry=chemistry,
        variants=Variants(organs={Organ.SMALL_INTESTINE: SegmentedIntestine(ka=ka)}),
    )
    r = roadrunner.RoadRunner(str(build(d, tmp_path).flat_path))
    r.setValue("init(BW)", 140.0)
    r.reset()
    assert r["GU__DRG2M_GU_k"] == pytest.approx(1.0)
    assert r["ka_drg"] == pytest.approx(2.0)


def test_the_validation_skips_the_numbers_of_a_formula() -> None:
    """A formula has no number before the simulation: no positivity or interval check.

    Its declared unit is checked as that of a `Value`.
    """
    drug = dataclasses.replace(
        simple_drug(),
        binding=Binding(fu=Formula("fu_ref_drg_m", U.dimensionless)),
        distribution=Kp(
            Formula("Kp_drg_m", U.dimensionless),
            ftissue=1.0,
            ps={Organ.LIVER: Formula("ftissue_drg / 1 litre", U.per_min)},
        ),
        excretion=Excretion(
            urine=True, fraction_reabsorbed=Formula("fu_ref_drg_m", U.dimensionless)
        ),
    )
    metabolite = Substance(
        sid="drg_m",
        name="metabolite",
        mr=316.0,
        binding=Binding(fu=0.5),
        distribution=Kp(0.7, ftissue=1.0),
    )
    definition(substances=[drug, metabolite]).validate()
    wrong = dataclasses.replace(
        drug,
        distribution=Kp(
            1.0, ftissue=1.0, ps={Organ.LIVER: Formula("ftissue_drg", U.l_per_min)}
        ),
    )
    with pytest.raises(ValueError, match=r"ps 'li' must be in 1/min"):
        definition(substances=[wrong, metabolite]).validate()


def test_a_formula_per_substance_in_the_segmented_intestine(tmp_path: Path) -> None:
    """The `ka` of each oral substance is its own formula."""
    d = definition(
        routes=[Route.oral("drg"), Route.oral("drg_m")],
        variants=Variants(
            organs={
                Organ.SMALL_INTESTINE: SegmentedIntestine(
                    ka={
                        "drg": Formula("1 per_hr * BW / 70 kg", U.per_hr),
                        "drg_m": 0.5,
                    }
                )
            }
        ),
    )
    r = roadrunner.RoadRunner(str(build(d, tmp_path).flat_path))
    r.setValue("init(BW)", 140.0)
    r.reset()
    assert r["ka_drg"] == pytest.approx(2.0)
    assert r["ka_drg_m"] == pytest.approx(0.5)


def _k_prod(value: float = 0.025) -> ModelParameter:
    """The constant `k_prod` of a definition (mM)."""
    return ModelParameter("k_prod", value, U.mM, name="production constant")


def test_a_formula_baseline_reads_a_definition_parameter(tmp_path: Path) -> None:
    """A baseline `2 * k_prod` reads the parameter of the definition.

    A change of `k_prod` before the simulation (`preinit_changes`) moves the
    baseline and every start with it.
    """
    s = Substance(
        "end",
        "endogenous substance",
        150.0,
        distribution=Kp(2.0, ftissue=1.0),
        baseline=Formula("2 dimensionless * k_prod", U.mM),
    )
    d = PBPKDefinition(sid="bl_parameter", substances=[s], parameters=[_k_prod()])
    path = build(d, tmp_path).flat_path
    plain = individual(path)
    doubled = individual(path, k_prod=0.05)
    venous = f"[{ids.pool_species(Pool.VENOUS, 'end')}]"
    assert plain[ids.baseline("end")] == pytest.approx(0.05)
    assert plain[venous] == pytest.approx(0.05)
    assert plain["[MU__Cti_end]"] == pytest.approx(0.1)
    assert doubled[ids.baseline("end")] == pytest.approx(0.1)
    assert doubled[venous] == pytest.approx(0.1)
    assert doubled["[MU__Cti_end]"] == pytest.approx(0.2)


def test_a_custom_reaction_reads_a_definition_parameter(tmp_path: Path) -> None:
    """A production in the muscle reads `k_prod` of the definition through an input port.

    The substance is not eliminated, so the body holds `k_prod * t`; a
    change of `k_prod` before the simulation doubles it. The muscle reads
    the parameter under its id, the flat model has no `MU__k_prod`.
    """
    s = Substance("prd", "produced substance", 150.0, distribution=Kp(2.0, ftissue=1.0))
    production = CustomReaction(
        "PRDPROD",
        Organ.MUSCLE,
        Part.TISSUE,
        "-> Cti_prd",
        "k_prod",
        parameters=[],
    )
    d = PBPKDefinition(
        sid="custom_parameter",
        substances=[s],
        chemistry=[production],
        parameters=[
            ModelParameter("k_prod", 0.002, U.mmole_per_min, name="production rate")
        ],
    )
    result = build(d, tmp_path)
    comp = libsbml.readSBMLFromFile(str(result.comp_path))
    muscle = comp.getPlugin("comp").getModelDefinition("standard_mu")
    port = muscle.getPlugin("comp").getPort("k_prod_port")
    assert port is not None and port.getIdRef() == "k_prod"
    flat = _model(result.flat_path)
    assert flat.getParameter("MU__k_prod") is None
    reaction = flat.getReaction("MU__PRDPROD")
    assert libsbml.formulaToL3String(reaction.getKineticLaw().getMath()) == "k_prod"
    for k_prod in (0.002, 0.004):
        r = individual(result.flat_path, k_prod=k_prod)
        r.integrator.absolute_tolerance = 1e-12
        r.integrator.relative_tolerance = 1e-10
        r.simulate(0, 100)
        assert total_amount(r, "prd") == pytest.approx(k_prod * 100, rel=1e-6)


def test_a_custom_reaction_reads_a_covariate_and_a_compartment(tmp_path: Path) -> None:
    """A custom reaction reads `BW` and the compartment `Vurine` of the top model.

    The muscle reads the covariate through a parameter input and the
    compartment through a compartment input, each under its id; a change of
    `BW` before the simulation doubles the production.
    """
    s = Substance("prd", "produced substance", 150.0, distribution=Kp(2.0, ftissue=1.0))
    production = CustomReaction(
        "PRDPROD",
        Organ.MUSCLE,
        Part.TISSUE,
        "-> Cti_prd",
        "0.002 mmole_per_min * BW / 70 kg * Vurine / Vurine",
        parameters=[],
    )
    d = PBPKDefinition(sid="custom_top", substances=[s], chemistry=[production])
    result = build(d, tmp_path)
    comp = libsbml.readSBMLFromFile(str(result.comp_path))
    muscle = comp.getPlugin("comp").getModelDefinition("standard_mu")
    assert muscle.getCompartment("Vurine") is not None
    assert muscle.getParameter("BW") is not None
    ports = muscle.getPlugin("comp")
    assert ports.getPort("Vurine_port").getIdRef() == "Vurine"
    assert ports.getPort("BW_port").getIdRef() == "BW"
    flat = _model(result.flat_path)
    assert flat.getCompartment("MU__Vurine") is None
    assert flat.getParameter("MU__BW") is None
    for bw in (70.0, 140.0):
        r = individual(result.flat_path, BW=bw)
        r.integrator.absolute_tolerance = 1e-12
        r.integrator.relative_tolerance = 1e-10
        r.simulate(0, 100)
        expected = 0.002 * bw / 70.0 * 100
        assert total_amount(r, "prd") == pytest.approx(expected, rel=1e-6)


def test_an_organ_formula_reads_a_definition_parameter(tmp_path: Path) -> None:
    """A kinetic constant of the liver is a formula over a parameter of the definition.

    The liver reads `k_bl` through an input port under its id; a change
    before the simulation is followed.
    """
    k = Formula("2 dimensionless * k_bl", U.per_min)
    d = dataclasses.replace(
        _definition(k=k),
        parameters=[ModelParameter("k_bl", 0.005, U.per_min, name="bile constant")],
    )
    result = build(d, tmp_path)
    comp = libsbml.readSBMLFromFile(str(result.comp_path))
    liver = comp.getPlugin("comp").getModelDefinition("standard_li")
    assert liver.getPlugin("comp").getPort("k_bl_port").getIdRef() == "k_bl"
    model = _model(result.flat_path)
    assert formula_ids(str(_formula(model, "LI__PMBIL_k"))) == {"k_bl"}
    r = individual(result.flat_path)
    assert r["LI__PMBIL_k"] == pytest.approx(0.01)
    r = individual(result.flat_path, k_bl=0.01)
    assert r["LI__PMBIL_k"] == pytest.approx(0.02)
