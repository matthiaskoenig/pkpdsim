"""Tests of the standard organ: collapsing of its parts, boundary and placement."""

import dataclasses
from collections.abc import Sequence
from pathlib import Path

import pytest
from sbmlutils.factory import Model, ModelDefinition, Parameter, Reaction, create_model

from pkpdsim.pbpk import (
    MM,
    Binding,
    CustomReaction,
    Direction,
    Erythrocytes,
    FirstOrder,
    Hill,
    Kp,
    Membrane,
    Metabolism,
    Mode,
    Organ,
    OrganContext,
    Substance,
    Transport,
    Value,
    parts_for,
)
from pkpdsim.pbpk.organ import Boundary, StandardOrgan
from pkpdsim.pbpk.placement import place
from pkpdsim.pbpk.units import MODEL_UNITS, U

from .conftest import simple_drug

Process = Metabolism | Transport | CustomReaction

LIVER = Organ.LIVER


def _metabolite() -> Substance:
    """A metabolite without tissue distribution."""
    return Substance(sid="drg_m", name="metabolite", mr=316.0)


def _drug(**changes: object) -> Substance:
    """The simple drug with changed fields."""
    return dataclasses.replace(simple_drug(), **changes)


def _metabolism(organ: Organ = LIVER, **changes: object) -> Metabolism:
    """The metabolism `DRG2M` of the drug into the metabolite (CYP3A4)."""
    m = Metabolism(
        sid="DRG2M",
        organ=organ,
        substrate="drg",
        products={"drg_m": 1.0},
        kinetics=MM(vmax=1e-3, km=0.01),
        protein="CYP3A4",
    )
    return dataclasses.replace(m, **changes)


def _transport(
    direction: Direction = Direction.INFLUX,
    membrane: Membrane = Membrane.BASOLATERAL,
    organ: Organ = LIVER,
) -> Transport:
    """A transport of the drug (OATP1B1)."""
    return Transport(
        sid="DRGT",
        organ=organ,
        substance="drg",
        direction=direction,
        membrane=membrane,
        kinetics=MM(vmax=1e-3, km=0.01),
        protein="OATP1B1",
    )


def _ctx(
    organ: Organ,
    substances: Sequence[Substance] = (),
    chemistry: Sequence[Process] = (),
) -> OrganContext:
    """The context of an organ, the parts by `parts_for`."""
    subs = tuple(substances) or (simple_drug(), _metabolite())
    return OrganContext(
        organ=organ,
        substances=subs,
        parts={s.sid: parts_for(s, organ, chemistry) for s in subs},
        chemistry=tuple(chemistry),
        inhibitions=(),
        inductions=(),
    )


def _build(
    organ: Organ = LIVER,
    substances: Sequence[Substance] = (),
    chemistry: Sequence[Process] = (),
) -> tuple[ModelDefinition, Boundary]:
    """Build the standard organ."""
    return StandardOrgan().build(_ctx(organ, substances, chemistry))


def _species(md: ModelDefinition) -> dict[str, str]:
    """The species of a model definition with their compartment."""
    return {str(s.sid): s.compartment for s in md.species or []}


def _reaction(md: ModelDefinition, sid: str) -> Reaction:
    """A reaction of a model definition by id."""
    return next(r for r in md.reactions or [] if r.sid == sid)


def _formula(reaction: Reaction) -> str:
    """The math of the kinetic law of a reaction."""
    assert reaction.formula is not None
    return str(reaction.formula.math)


def _equation(
    reaction: Reaction,
) -> tuple[list[tuple[str, float | None]], list[tuple[str, float | None]]]:
    """The reactants and products of a reaction with their stoichiometry."""
    eq = reaction.equation
    return (
        [(p.species, p.stoichiometry) for p in eq.reactants],
        [(p.species, p.stoichiometry) for p in eq.products],
    )


def _rules(md: ModelDefinition) -> dict[str, str]:
    """The assignment rules and initial assignments of a model definition."""
    rules = {r.variable: str(r.value) for r in md.rules or []}
    return rules | {a.symbol: str(a.value) for a in md.assignments or []}


def _parameters(md: ModelDefinition) -> dict[str, Parameter]:
    """The parameters of a model definition, also those of the reactions."""
    pars = {str(p.sid): p for p in md.parameters or []}
    for r in md.reactions or []:
        pars |= {str(p.sid): p for p in r.pars or []}
    return pars


def _validate(md: ModelDefinition, tmp_path: Path) -> None:
    """The model definition is valid SBML (with comp) and its units are consistent.

    libsbml checks the units of the main model only, so the content of the
    definition is validated as a main model as well.
    """
    container = Model(
        sid="organ_test",
        name="organ test",
        units=U,
        model_units=MODEL_UNITS,
        model_definitions=[md],
    )
    main = Model(
        sid=md.sid,
        name=md.name,
        units=U,
        model_units=MODEL_UNITS,
        compartments=md.compartments,
        species=md.species,
        parameters=md.parameters,
        rules=md.rules,
        assignments=md.assignments,
        reactions=md.reactions,
    )
    for model in (container, main):
        result = create_model(model, filepath=tmp_path / f"{model.sid}.xml")
        assert result.validation is not None
        assert result.validation.error_count == 0, model.sid
        assert result.validation.warning_count == 0, model.sid


# --- collapsing --------------------------------------------------------------


def test_merged_parts_species() -> None:
    """Merged tissue: plasma and tissue species, no blood cells, no interstitium."""
    md, _ = _build()
    species = _species(md)
    assert species["Cpl_drg"] == "Vpl"
    assert species["Cti_drg"] == "Vti"
    for sid in ("Cbc_drg", "Cis_drg", "Cic_drg", "Cor_drg"):
        assert sid not in species
    rules = _rules(md)
    assert rules["Vpl"] == "Vli * f_vas_li * (1 dimensionless - HCT)"
    assert rules["Vti"] == "Vli * (1 dimensionless - f_vas_li) * f_functional"
    assert _formula(_reaction(md, "J_ti_drg")) == (
        "PS_drg * fu_drg * (Cpl_drg - Cti_drg / Kp_drg)"
    )
    assert rules["PS_drg"] == (
        "ftissue_drg * Vli * (1 dimensionless - f_vas_li) * f_functional / "
        "Vti_total_ref_drg"
    )


def test_full_parts_species() -> None:
    """Permeability-limited: plasma, interstitial and intracellular species."""
    drug = _drug(distribution=Kp(2.0, mode=Mode.PERMEABILITY, ftissue=1.0))
    md, _ = _build(substances=[drug, _metabolite()])
    species = _species(md)
    assert species["Cpl_drg"] == "Vpl"
    assert species["Cis_drg"] == "Vis"
    assert species["Cic_drg"] == "Vic"
    assert "Cti_drg" not in species
    rules = _rules(md)
    assert rules["Vis"] == "Vli * f_int_li * f_functional"
    assert rules["Vic"] == (
        "Vli * (1 dimensionless - f_vas_li - f_int_li) * f_functional"
    )
    assert rules["K_int_pls_drg"] == (
        "(f_water_int + r_protein_int_pls * (1 dimensionless / fu_drg - "
        "f_water_pls)) * fu_drg"
    )
    f_cell = "(1 dimensionless - f_vas_li - f_int_li)"
    k_cell = f"(Kp_drg * (f_int_li + {f_cell}) - f_int_li * K_int_pls_drg) / {f_cell}"
    floor = "1e-06 dimensionless"
    assert rules["K_cell_pls_drg"] == (
        f"piecewise({k_cell}, {k_cell} > {floor}, {floor})"
    )
    notes = str(_parameters(md)["K_cell_pls_drg"].notes)
    assert "must exceed the interstitial share" in notes
    assert _formula(_reaction(md, "J_is_drg")) == (
        "PS_vas * fu_drg * (Cpl_drg - Cis_drg / K_int_pls_drg)"
    )
    assert _formula(_reaction(md, "J_ic_drg")) == (
        "PS_drg * (fu_drg / K_int_pls_drg * Cis_drg - fu_drg / K_cell_pls_drg * "
        "Cic_drg)"
    )
    pars = _parameters(md)
    assert pars["f_water_int"].value == 0.935
    assert "Schmitt" in str(pars["f_water_int"].notes)


def test_perfusion_single_species() -> None:
    """Flow-limited: one species for the whole organ, outflow `Cor/K_or`."""
    drug = _drug(distribution=Kp(2.0, mode=Mode.PERFUSION))
    md, boundary = _build(substances=[drug, _metabolite()])
    species = _species(md)
    assert species["Cor_drg"] == "Vor"
    assert "Cpl_drg" not in species
    assert "Cti_drg" not in species
    rules = _rules(md)
    assert rules["Vor"] == "Vli"
    assert rules["K_or_drg"] == (
        "f_vas_li * (1 dimensionless - HCT) + (1 dimensionless - f_vas_li) * "
        "Kp_drg * f_functional"
    )
    assert rules["Cout_drg"] == "Cor_drg / K_or_drg"
    assert _formula(_reaction(md, "J_out_pl_drg")) == (
        "(1 dimensionless - f_shunt) * Qh * (1 dimensionless - HCT) * Cout_drg"
    )
    assert boundary.outputs["Cout_drg"] == "Cout_li_drg"


def test_no_distribution_plasma_only() -> None:
    """A substance without distribution has its plasma species only."""
    md, _ = _build()
    species = _species(md)
    assert species["Cpl_drg_m"] == "Vpl"
    assert not {"Cti_drg_m", "Cis_drg_m", "Cic_drg_m", "Cor_drg_m"} & set(species)


def test_equilibrium_erythrocytes_flows_carry_blood() -> None:
    """Erythrocytes in equilibrium: `Q*BP*C`, capacity `V*f_vas*BP`."""
    drug = _drug(erythrocytes=Erythrocytes(mode="equilibrium", bp=1.5))
    md, boundary = _build(substances=[drug, _metabolite()])
    species = _species(md)
    assert species["Cpl_drg"] == "Vpl_drg"
    assert "Cbc_drg" not in species
    assert _rules(md)["Vpl_drg"] == "Vli * f_vas_li * BP_drg"
    assert _formula(_reaction(md, "J_in_ar_pl_drg")) == (
        "(1 dimensionless - f_shunt) * Qha * BP_drg * Car_drg"
    )
    assert _formula(_reaction(md, "J_out_pl_drg")) == (
        "(1 dimensionless - f_shunt) * Qh * BP_drg * Cout_drg"
    )
    assert species["Car_drg"] == "Var_pl_drg"
    assert boundary.replacements["Var_pl_drg"] == "Var_pl_drg"
    assert boundary.replacements["BP_drg"] == "BP_drg"


def test_equilibrium_erythrocytes_perfusion() -> None:
    """Flow-limited with erythrocytes in equilibrium: the blood cells hold drug."""
    drug = _drug(
        erythrocytes=Erythrocytes(mode="equilibrium", kbc=2.0),
        distribution=Kp(2.0, mode=Mode.PERFUSION),
    )
    md, _ = _build(substances=[drug, _metabolite()])
    assert _rules(md)["K_or_drg"] == (
        "f_vas_li * BP_drg + (1 dimensionless - f_vas_li) * Kp_drg * f_functional"
    )


def test_kinetic_erythrocytes_blood_cells() -> None:
    """Kinetic erythrocytes: blood cell species, flows with HCT and exchange."""
    drug = _drug(erythrocytes=Erythrocytes(mode="kinetic", kbc=2.0, k_exchange=0.5))
    md, boundary = _build(substances=[drug, _metabolite()])
    species = _species(md)
    assert species["Cbc_drg"] == "Vbc"
    assert species["Car_bc_drg"] == "Var_bc"
    assert species["Cpo_bc_drg"] == "Vpo_bc"
    assert species["Cve_bc_drg"] == "Vve_bc"
    assert _rules(md)["Vbc"] == "Vli * f_vas_li * HCT"
    assert _formula(_reaction(md, "J_in_ar_bc_drg")) == (
        "(1 dimensionless - f_shunt) * Qha * HCT * Car_bc_drg"
    )
    assert _formula(_reaction(md, "J_out_bc_drg")) == (
        "(1 dimensionless - f_shunt) * Qh * HCT * Cout_bc_drg"
    )
    assert _formula(_reaction(md, "J_shunt_po_bc_drg")) == (
        "f_shunt * Qpo * HCT * Cpo_bc_drg"
    )
    assert _rules(md)["Cout_bc_drg"] == "Cbc_drg"
    assert boundary.outputs["Cout_bc_drg"] == "Cout_li_bc_drg"
    assert _formula(_reaction(md, "J_bc_drg")) == (
        "kex_bc_drg * Vbc * (Cpl_drg - Cbc_drg / Kbc_drg)"
    )
    assert {"Kbc_drg", "kex_bc_drg", "Var_bc"} <= set(boundary.replacements)


def test_kp_override() -> None:
    """An organ with a Kp override reads `Kp_<o>_<s>`."""
    drug = _drug(distribution=Kp(2.0, ftissue=1.0, overrides={LIVER: 5.0}))
    md, boundary = _build(substances=[drug, _metabolite()])
    assert "Kp_li_drg" in _formula(_reaction(md, "J_ti_drg"))
    assert boundary.replacements["Kp_li_drg"] == "Kp_li_drg"
    assert "Kp_drg" not in boundary.replacements
    md_mu, boundary_mu = _build(Organ.MUSCLE, substances=[drug, _metabolite()])
    assert "Kp_drg" in _formula(_reaction(md_mu, "J_ti_drg"))
    assert "Kp_drg" in boundary_mu.replacements


def test_ps_of_organ() -> None:
    """A value of the organ replaces the share of `ftissue`."""
    ps = Value(0.2, U.l_per_min, source="Schmitt2008")
    drug = _drug(distribution=Kp(2.0, ftissue=1.0, ps={LIVER: ps}))
    md, boundary = _build(substances=[drug, _metabolite()])
    par = _parameters(md)["PS_drg"]
    assert par.value == 0.2
    assert "doi:10.1016/j.tiv.2007.09.010" in str(par.notes)
    assert "PS_drg" not in _rules(md)
    assert "ftissue_drg" not in boundary.replacements
    assert "Vti_total_ref_drg" not in boundary.replacements


def test_organ_not_entered() -> None:
    """An organ outside the organs of the distribution has plasma only."""
    drug = _drug(distribution=Kp(2.0, ftissue=1.0, organs=frozenset({Organ.MUSCLE})))
    md, _ = _build(substances=[drug, _metabolite()])
    assert "Cti_drg" not in _species(md)


def test_kp_below_interstitial_share_raises() -> None:
    """A Kp below the interstitial share has no positive `K_cell_pls`."""
    drug = _drug(distribution=Kp(0.01, mode=Mode.PERMEABILITY, ftissue=1.0))
    with pytest.raises(ValueError, match=r"drg.*'li'.*K_cell_pls"):
        _build(substances=[drug, _metabolite()])


def test_full_without_ftissue_raises() -> None:
    """A flow-limited substance forced into the full form needs a permeability."""
    drug = _drug(distribution=Kp(2.0, mode=Mode.PERFUSION))
    with pytest.raises(ValueError, match="ftissue"):
        _build(substances=[drug, _metabolite()], chemistry=[_transport()])


def test_transport_without_distribution() -> None:
    """A transporter gives a vascular substance cells without passive exchange."""
    md, _ = _build(
        substances=[_metabolite(), simple_drug()],
        chemistry=[
            Transport(
                sid="MT",
                organ=LIVER,
                substance="drg_m",
                direction=Direction.INFLUX,
                membrane=Membrane.BASOLATERAL,
                kinetics=FirstOrder(k=0.1),
            )
        ],
    )
    assert _parameters(md)["PS_drg_m"].value == 0.0
    assert _rules(md)["K_cell_pls_drg_m"] == "K_int_pls_drg_m"


# --- boundary ----------------------------------------------------------------


def test_boundary_names_inflow_and_outflow() -> None:
    """The lung reads venous and writes arterial blood, the portal organs the portal vein."""
    md, boundary = _build(Organ.LUNG)
    assert _equation(_reaction(md, "J_in_ve_pl_drg")) == (
        [("Cve_drg", 1.0)],
        [("Cpl_drg", 1.0)],
    )
    assert _equation(_reaction(md, "J_out_pl_drg"))[1] == [("Car_drg", 1.0)]
    assert _formula(_reaction(md, "J_out_pl_drg")) == (
        "Qlu * (1 dimensionless - HCT) * Cout_drg"
    )
    assert {"Cve_drg", "Car_drg", "Vve_pl", "Var_pl", "Qlu"} <= set(
        boundary.replacements
    )

    md, boundary = _build(LIVER)
    assert _formula(_reaction(md, "J_in_ar_pl_drg")) == (
        "(1 dimensionless - f_shunt) * Qha * (1 dimensionless - HCT) * Car_drg"
    )
    assert _formula(_reaction(md, "J_in_po_pl_drg")) == (
        "(1 dimensionless - f_shunt) * Qpo * (1 dimensionless - HCT) * Cpo_drg"
    )
    assert _equation(_reaction(md, "J_out_pl_drg"))[1] == [("Cve_drg", 1.0)]
    assert {"Car_drg", "Cpo_drg", "Cve_drg", "Qha", "Qpo", "Qh"} <= set(
        boundary.replacements
    )

    md, boundary = _build(Organ.SPLEEN)
    assert _equation(_reaction(md, "J_in_ar_pl_drg"))[0] == [("Car_drg", 1.0)]
    assert _equation(_reaction(md, "J_out_pl_drg"))[1] == [("Cpo_drg", 1.0)]
    assert "Cve_drg" not in boundary.replacements


def test_boundary_physiology_and_substance_parameters() -> None:
    """Organ volume, flow, HCT, fractions and substance parameters are replaced."""
    md, boundary = _build(Organ.MUSCLE)
    for pid in (
        "Vmu",
        "Qmu",
        "HCT",
        "f_vas_mu",
        "fu_drg",
        "Kp_drg",
        "ftissue_drg",
        "Vti_total_ref_drg",
        "fu_drg_m",
    ):
        assert boundary.replacements[pid] == pid
    assert boundary.outputs == {
        "Cout_drg": "Cout_mu_drg",
        "Cout_drg_m": "Cout_mu_drg_m",
    }
    pars = _parameters(md)
    for pid in boundary.replacements:
        element = pars.get(pid) or next(
            e for e in (md.compartments or []) + (md.species or []) if e.sid == pid
        )
        assert element.port is not None, pid
    assert pars["Cout_drg"].port is not None


# --- placement ---------------------------------------------------------------


def test_metabolism_on_unbound_with_activity() -> None:
    """The rate acts on the unbound tissue concentration and has the activity.

    The activity is the organ-local base activity times the impairment of the
    liver, an input of the body (`fimp_li_cyp3a4`).
    """
    md, boundary = _build(chemistry=[_metabolism()])
    reaction = _reaction(md, "DRG2M")
    formula = _formula(reaction)
    assert "f_cyp3a4" in formula
    assert "fu_drg" in formula
    assert formula == (
        "f_cyp3a4 * fimp_cyp3a4 * f_functional * DRG2M_Vmax * Vli * "
        "(fu_drg / Kp_drg * Cti_drg) / (DRG2M_Km + fu_drg / Kp_drg * Cti_drg)"
    )
    assert reaction.compartment == "Vti"
    assert _parameters(md)["f_cyp3a4"].port is True
    assert "f_cyp3a4" not in boundary.replacements
    assert boundary.replacements["fimp_cyp3a4"] == "fimp_li_cyp3a4"


@pytest.mark.parametrize(
    ("distribution", "part", "unbound"),
    [
        (Kp(2.0, mode=Mode.PERMEABILITY, ftissue=1.0), "ic", "fu_drg / K_cell_pls_drg"),
        (Kp(2.0, ftissue=1.0), "ti", "fu_drg / Kp_drg"),
        (Kp(2.0, mode=Mode.PERFUSION), "or", "fu_drg / K_or_drg"),
        (None, "pl", "fu_drg"),
    ],
)
def test_metabolism_compartment(
    distribution: Kp | None, part: str, unbound: str
) -> None:
    """The metabolism acts in the part of the tissue of the substrate."""
    drug = _drug(distribution=distribution)
    md, _ = _build(substances=[drug, _metabolite()], chemistry=[_metabolism()])
    reaction = _reaction(md, "DRG2M")
    assert reaction.compartment == f"V{part}"
    assert f"({unbound} * C{part}_drg)" in _formula(reaction)


def test_product_without_tissue_goes_to_plasma() -> None:
    """A product without species in the part of the substrate goes into plasma."""
    md, _ = _build(chemistry=[_metabolism()])
    assert _equation(_reaction(md, "DRG2M")) == (
        [("Cti_drg", 1.0)],
        [("Cpl_drg_m", 1.0)],
    )


def test_product_with_tissue_stays_in_part() -> None:
    """A product with a species in the part of the substrate stays there."""
    metabolite = dataclasses.replace(_metabolite(), distribution=Kp(1.0, ftissue=1.0))
    other = Substance(sid="drg_n", name="other", mr=300.0)
    md, _ = _build(
        substances=[simple_drug(), metabolite, other],
        chemistry=[_metabolism(products={"drg_m": 0.25, "drg_n": 0.75})],
    )
    assert _equation(_reaction(md, "DRG2M"))[1] == [
        ("Cti_drg_m", 0.25),
        ("Cpl_drg_n", 0.75),
    ]


def test_product_into_equivalent_tissue_part() -> None:
    """A product goes into its equivalent tissue part before plasma (R19)."""
    full = Kp(2.0, mode=Mode.PERMEABILITY, ftissue=1.0)
    merged = Kp(1.0, ftissue=1.0)
    perfusion = Kp(1.0, mode=Mode.PERFUSION)
    cases = [
        (full, merged, "Cic_drg", "Cti_drg_m"),
        (merged, full, "Cti_drg", "Cic_drg_m"),
        (full, perfusion, "Cic_drg", "Cor_drg_m"),
        (perfusion, full, "Cor_drg", "Cic_drg_m"),
        (None, perfusion, "Cpl_drg", "Cor_drg_m"),
    ]
    for parent, child, substrate, product in cases:
        drug = _drug(distribution=parent)
        metabolite = dataclasses.replace(_metabolite(), distribution=child)
        md, _ = _build(substances=[drug, metabolite], chemistry=[_metabolism()])
        assert _equation(_reaction(md, "DRG2M")) == (
            [(substrate, 1.0)],
            [(product, 1.0)],
        ), (parent, child)


def test_metabolism_kinetics() -> None:
    """First order and Hill kinetics, an elimination without products."""
    first = _metabolism(sid="E1", products={}, kinetics=FirstOrder(k=0.5), protein=None)
    hill = _metabolism(sid="E2", products={}, kinetics=Hill(vmax=1.0, km=0.1, n=2.0))
    md, _ = _build(chemistry=[first, hill])
    assert _formula(_reaction(md, "E1")) == (
        "f_functional * E1_k * fu_drg / Kp_drg * Cti_drg"
    )
    assert _equation(_reaction(md, "E1")) == ([("Cti_drg", 1.0)], [])
    assert _formula(_reaction(md, "E2")) == (
        "f_cyp3a4 * fimp_cyp3a4 * f_functional * E2_Vmax * Vli * ((fu_drg / Kp_drg * "
        "Cti_drg) / E2_Km)^E2_n / (1 dimensionless + ((fu_drg / Kp_drg * Cti_drg) / "
        "E2_Km)^E2_n)"
    )


def test_first_order_per_volume() -> None:
    """A first order clearance per litre of organ: the rate `k*V<o>*Cu`, `k` in 1/min."""
    first = _metabolism(
        sid="E1",
        products={},
        kinetics=FirstOrder(k=0.5, per_volume=True),
        protein=None,
    )
    whole = _metabolism(sid="E2", products={}, kinetics=FirstOrder(k=0.5), protein=None)
    md, _ = _build(chemistry=[first, whole])
    assert _formula(_reaction(md, "E1")) == (
        "f_functional * E1_k * Vli * fu_drg / Kp_drg * Cti_drg"
    )
    pars = _parameters(md)
    assert pars["E1_k"].unit is U.per_min
    assert pars["E2_k"].unit is U.l_per_min
    assert "per litre of organ" in str(pars["E1_k"].notes)


def test_kidney_processes_follow_the_renal_function() -> None:
    """Every metabolism and transport of the kidney runs with the tubular capacity.

    `1 - f_tubular + f_tubular*f_renal_function`, the share `f_tubular` 1 of
    the intact nephron hypothesis with its provenance (a model choice); the
    renal function is an input of the body. A kidney without processes reads
    no renal function.
    """
    elimination = _metabolism(
        Organ.KIDNEY, products={}, kinetics=FirstOrder(k=0.5), protein=None
    )
    secretion = _transport(Direction.EFFLUX, Membrane.APICAL, Organ.KIDNEY)
    md, boundary = _build(Organ.KIDNEY, chemistry=[elimination, secretion])
    capacity = "(1 dimensionless - f_tubular + f_tubular * f_renal_function)"
    # the transport gives the drug the full form: the metabolism acts in the cell
    assert _formula(_reaction(md, "DRG2M")) == (
        f"{capacity} * DRG2M_k * fu_drg / K_cell_pls_drg * Cic_drg"
    )
    assert _formula(_reaction(md, "DRGT")).startswith(f"f_oatp1b1 * {capacity} * ")
    assert boundary.replacements["f_renal_function"] == "f_renal_function"
    share = _parameters(md)["f_tubular"]
    assert share.value == 1.0
    assert "intact nephron hypothesis (Bricker 1960), not verified" in str(share.notes)
    assert "Choice of the pkpdsim PBPK library" in str(share.notes)
    md, boundary = _build(Organ.KIDNEY)
    assert "f_renal_function" not in boundary.replacements
    assert "f_tubular" not in _parameters(md)


def test_basolateral_transport() -> None:
    """Basolateral influx from the interstitium, efflux into it."""
    influx = _transport()
    efflux = dataclasses.replace(_transport(Direction.EFFLUX), sid="DRGE")
    md, _ = _build(chemistry=[influx, efflux])
    assert _equation(_reaction(md, "DRGT")) == ([("Cis_drg", 1.0)], [("Cic_drg", 1.0)])
    assert "fu_drg / K_int_pls_drg * Cis_drg" in _formula(_reaction(md, "DRGT"))
    assert "f_oatp1b1" in _formula(_reaction(md, "DRGT"))
    assert _equation(_reaction(md, "DRGE")) == ([("Cic_drg", 1.0)], [("Cis_drg", 1.0)])
    assert "fu_drg / K_cell_pls_drg * Cic_drg" in _formula(_reaction(md, "DRGE"))


@pytest.mark.parametrize(
    ("organ", "sink", "compartment"),
    [
        (LIVER, "Cbile_drg", "Vbile"),
        (Organ.KIDNEY, "Aurine_drg", "Vurine"),
        (Organ.SMALL_INTESTINE, "Clumen_drg", "Vlumen"),
    ],
)
def test_apical_efflux_writes_to_sink_port(
    organ: Organ, sink: str, compartment: str
) -> None:
    """Apical efflux moves the drug from the cell into the sink of the organ (a port)."""
    md, boundary = _build(
        organ, chemistry=[_transport(Direction.EFFLUX, Membrane.APICAL, organ)]
    )
    assert _equation(_reaction(md, "DRGT")) == ([("Cic_drg", 1.0)], [(sink, 1.0)])
    species = {s.sid: s for s in md.species or []}
    assert species[sink].compartment == compartment
    assert species[sink].port is True
    assert species[sink].hasOnlySubstanceUnits is (organ is Organ.KIDNEY)
    assert boundary.replacements[sink] == sink
    assert boundary.replacements[compartment] == compartment


def test_apical_influx_from_sink() -> None:
    """Apical influx takes up the drug from the bile (a transit compartment)."""
    md, _ = _build(chemistry=[_transport(Direction.INFLUX, Membrane.APICAL)])
    reaction = _reaction(md, "DRGT")
    assert _equation(reaction) == ([("Cbile_drg", 1.0)], [("Cic_drg", 1.0)])
    assert "(Cbile_drg) / (DRGT_Km + Cbile_drg)" in _formula(reaction)


def test_apical_influx_from_urine_raises() -> None:
    """The standard kidney takes up nothing from the cumulative urine (R18)."""
    with pytest.raises(ValueError, match=r"'DRGT'.*fraction_reabsorbed.*nephron"):
        _build(
            Organ.KIDNEY,
            chemistry=[_transport(Direction.INFLUX, Membrane.APICAL, Organ.KIDNEY)],
        )


def test_apical_transport_without_sink_raises() -> None:
    """An organ without apical sink rejects apical transport."""
    with pytest.raises(ValueError, match="'mu' has no apical sink"):
        _build(
            Organ.MUSCLE,
            chemistry=[_transport(Direction.EFFLUX, Membrane.APICAL, Organ.MUSCLE)],
        )


def test_custom_reaction_in_part() -> None:
    """A custom reaction is placed as given in its part.

    Its rate is its own: in the liver `f_functional` does not scale it.
    """
    from pkpdsim.pbpk import Part

    custom = CustomReaction(
        sid="CUSTOM",
        organ=LIVER,
        part=Part.TISSUE,
        equation="Cti_drg -> Cpl_drg_m",
        formula="k_custom * Vti * Cti_drg",
        parameters=[Parameter("k_custom", 0.1, U.per_min)],
    )
    md, _ = _build(chemistry=[custom])
    reaction = _reaction(md, "CUSTOM")
    assert reaction.compartment == "Vti"
    assert _formula(reaction) == "k_custom * Vti * Cti_drg"
    bad = dataclasses.replace(custom, part=Part.INTRACELLULAR)
    with pytest.raises(ValueError, match="CUSTOM"):
        _build(chemistry=[bad])


def test_place_rejects_other_organ() -> None:
    """A process of another organ is not placed."""
    ctx = _ctx(LIVER)
    with pytest.raises(ValueError, match="'ki'"):
        place(_metabolism(Organ.KIDNEY), ctx, ctx.parts)


def test_chemistry_of_other_organs_ignored() -> None:
    """The organ builds only its own processes."""
    md, _ = _build(Organ.MUSCLE, chemistry=[_metabolism()])
    assert all(r.sid != "DRG2M" for r in md.reactions or [])


# --- SBML ----------------------------------------------------------------------


@pytest.mark.parametrize(
    ("organ", "drug", "chemistry"),
    [
        (LIVER, simple_drug(), [_metabolism()]),
        (
            LIVER,
            _drug(
                distribution=Kp(2.0, mode=Mode.PERMEABILITY, ftissue=1.0),
                erythrocytes=Erythrocytes(mode="kinetic", kbc=2.0, k_exchange=0.5),
            ),
            [
                _metabolism(kinetics=Hill(vmax=1.0, km=0.1, n=2.0)),
                _transport(),
                dataclasses.replace(
                    _transport(Direction.EFFLUX, Membrane.APICAL), sid="DRGB"
                ),
            ],
        ),
        (
            Organ.KIDNEY,
            _drug(binding=Binding(fu=0.5)),
            [
                _transport(Direction.EFFLUX, Membrane.APICAL, Organ.KIDNEY),
                _metabolism(Organ.KIDNEY, kinetics=FirstOrder(k=0.1)),
            ],
        ),
        (
            Organ.LUNG,
            _drug(
                distribution=Kp(2.0, mode=Mode.PERFUSION),
                erythrocytes=Erythrocytes(mode="equilibrium", bp=1.2),
            ),
            [],
        ),
        (
            Organ.SMALL_INTESTINE,
            _drug(erythrocytes=Erythrocytes(mode="equilibrium", bp=1.2)),
            [_transport(Direction.EFFLUX, Membrane.APICAL, Organ.SMALL_INTESTINE)],
        ),
    ],
    ids=["merged", "full", "kidney", "perfusion", "intestine"],
)
def test_valid_sbml(
    organ: Organ, drug: Substance, chemistry: list[Process], tmp_path: Path
) -> None:
    """The model definitions are valid SBML with consistent units."""
    md, _ = _build(organ, substances=[drug, _metabolite()], chemistry=chemistry)
    _validate(md, tmp_path)


def test_custom_species_collision_raises() -> None:
    """A species of a custom reaction must not replace a species of the organ."""
    from sbmlutils.factory import Species

    from pkpdsim.pbpk import Part

    custom = CustomReaction(
        sid="CUSTOM",
        organ=LIVER,
        part=Part.TISSUE,
        equation="Cti_drg -> Cpl_drg_m",
        formula="k_custom * Vti * Cti_drg",
        parameters=[Parameter("k_custom", 0.1, U.per_min)],
        species=[Species("Cti_drg", "Vti", initialConcentration=0.0)],
    )
    with pytest.raises(ValueError, match=r"CUSTOM.*Cti_drg"):
        _build(chemistry=[custom])


def test_blood_capacities_not_physical() -> None:
    """The capacities of the blood with BP are implicit, not physical compartments."""
    from sbmlutils.metadata import SBO

    drug = _drug(erythrocytes=Erythrocytes(mode="equilibrium", bp=1.5))
    md, _ = _build(substances=[drug, _metabolite()])
    compartments = {c.sid: c for c in md.compartments or []}
    for cid in ("Vpl_drg", "Var_pl_drg", "Vpo_pl_drg", "Vve_pl_drg"):
        assert compartments[cid].sboTerm == SBO.IMPLICIT_COMPARTMENT, cid
        assert "plasma-equivalent capacity" in str(compartments[cid].name), cid
        assert "Not a physical volume" in str(compartments[cid].notes), cid
    assert compartments["Vpl"].sboTerm == SBO.PHYSICAL_COMPARTMENT


def test_ftissue_meaning_in_notes() -> None:
    """The notes say that ftissue multiplies the unbound driving force (R20)."""
    md, _ = _build()
    pars = _parameters(md)
    for pid in ("ftissue_drg", "PS_drg"):
        assert "fu*(Cpl - Cti/Kp)" in str(pars[pid].notes), pid


def test_liver_shunt_past_the_liver() -> None:
    """The shunt flows from the arterial and the portal pool into the venous pool."""
    md, boundary = _build(LIVER)
    for pool, flow in (("ar", "Qha"), ("po", "Qpo")):
        reaction = _reaction(md, f"J_shunt_{pool}_pl_drg")
        assert _equation(reaction) == (
            [(f"C{pool}_drg", 1.0)],
            [("Cve_drg", 1.0)],
        )
        assert _formula(reaction) == (
            f"f_shunt * {flow} * (1 dimensionless - HCT) * C{pool}_drg"
        )
    assert {"f_shunt", "f_functional"} <= set(boundary.replacements)
    md, boundary = _build(Organ.MUSCLE)
    assert not [r for r in md.reactions or [] if str(r.sid).startswith("J_shunt")]
    assert not {"f_shunt", "f_functional"} & set(boundary.replacements)


@pytest.mark.parametrize(
    ("distribution", "source", "plasma"),
    [
        (Kp(2.0, ftissue=1.0), "Cpl_drg", "Cpl_drg"),
        (Kp(2.0, mode=Mode.PERFUSION), "Cor_drg", "Cor_drg / K_or_drg"),
    ],
    ids=["merged", "perfusion"],
)
def test_kidney_filtration(distribution: Kp, source: str, plasma: str) -> None:
    """The kidney filters `fu*GFR*Cpl*(1 - f_reabsorbed)` into the urine."""
    md, boundary = _build(Organ.KIDNEY, substances=[_drug(distribution=distribution)])
    reaction = _reaction(md, "J_gfr_drg")
    assert _equation(reaction) == ([(source, 1.0)], [("Aurine_drg", 1.0)])
    assert _formula(reaction) == (
        f"fu_drg * GFR * {plasma} * (1 dimensionless - f_reabsorbed_drg)"
    )
    assert {"GFR", "f_reabsorbed_drg", "Aurine_drg", "Vurine"} <= set(
        boundary.replacements
    )
    md, _ = _build(Organ.LIVER, substances=[_drug(distribution=distribution)])
    assert not [r for r in md.reactions or [] if r.sid == "J_gfr_drg"]
