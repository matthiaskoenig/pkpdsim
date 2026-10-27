"""Tests of the id scheme."""

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.organs import Organ, Part, Pool
from pkpdsim.pbpk.physiology import physiology_objects


def test_examples() -> None:
    """The documented examples."""
    assert ids.submodel(Organ.LIVER) == "LI"
    assert ids.submodel(Organ.LARGE_INTESTINE) == "LA"
    assert ids.flat(Organ.LIVER, "Cpl_api") == "LI__Cpl_api"
    assert ids.compartment(Part.INTERSTITIAL) == "Vis"
    assert ids.species(Part.PLASMA, "api") == "Cpl_api"
    assert ids.pool_compartment(Pool.VENOUS, Part.PLASMA) == "Vve_pl"
    assert ids.pool_species(Pool.VENOUS, "api") == "Cve_api"
    assert ids.pool_species(Pool.VENOUS, "api", Part.BLOOD_CELLS) == "Cve_bc_api"
    assert ids.organ_volume(Organ.LIVER) == "Vli"
    assert ids.organ_flow(Organ.LIVER) == "Qli"
    assert ids.outflow_conc("api") == "Cout_api"
    assert ids.activity("CYP3A4") == "f_cyp3a4"


def test_physiology_ids() -> None:
    """The ids of the physiology: reference values, fractions, masses and pools."""
    assert ids.ref(ids.FFM) == "FFM_ref"
    assert ids.ref(ids.organ_volume(Organ.LIVER)) == "Vli_ref"
    assert ids.organ_mass(Organ.ADIPOSE) == "Mad"
    assert ids.flow_fraction(Organ.LIVER) == "fQli"
    assert ids.pool_volume(Pool.VENOUS) == "Vve"
    assert ids.pool_fraction(Pool.PORTAL) == "fVpo"
    assert ids.VBLOOD == "Vblood"
    assert ids.SV_TIS == "SV_tis"
    assert ids.F_CARDIAC_FUNCTION == "f_cardiac_function"


def test_substance_parameters() -> None:
    """Parameter, sink and observer ids."""
    assert ids.mr("api") == "Mr_api"
    assert ids.fu("api") == "fu_api"
    assert ids.fu_ref("api") == "fu_ref_api"
    assert ids.kp("api") == "Kp_api"
    assert ids.kp_organ(Organ.LIVER, "api") == "Kp_li_api"
    assert ids.ftissue("api") == "ftissue_api"
    assert ids.bp("api") == "BP_api"
    assert ids.kbc("api") == "Kbc_api"
    assert ids.urine("api") == "Aurine_api"
    assert ids.feces("api") == "Afeces_api"
    assert ids.bile("api") == "Cbile_api"
    assert ids.arm("api") == "Carm_api"
    assert ids.arm_free("api") == "Carm_free_api"
    assert ids.arm_blood("api") == "Carm_blood_api"
    assert ids.free(Pool.VENOUS, "api") == "Cve_free_api"
    assert ids.blood(Pool.ARTERIAL, "api") == "Car_blood_api"
    assert ids.hepatic_vein("api") == "Chv_api"


def test_species_of_exact_substance() -> None:
    """A substance does not match the species of another one with its id as suffix."""
    all_ids = [
        "LI__Cpl_api",
        "LI__Cti_api",
        "LI__Cpl_z1_api",
        "LI__Cbc_z12_api",
        "LI__Cu_ti_z1_api",
        "Cve_api",
        "Cve_bc_api",
        "Car_bc_api",
        "LI__Cpl_api_m2",
        "LI__Cti_z2_api_m2",
        "Cve_api_m2",
        "Cpl_xdrg",
        "Cpl_drg",
        "Cve_free_api",
        "Carm_api",
        "Mr_api",
    ]
    assert ids.species_of("api", all_ids) == [
        "LI__Cpl_api",
        "LI__Cti_api",
        "LI__Cpl_z1_api",
        "LI__Cbc_z12_api",
        "Cve_api",
        "Cve_bc_api",
        "Car_bc_api",
    ]
    assert ids.species_of("api_m2", all_ids) == [
        "LI__Cpl_api_m2",
        "LI__Cti_z2_api_m2",
        "Cve_api_m2",
    ]
    assert ids.species_of("drg", all_ids) == ["Cpl_drg"]
    assert ids.species_of("m2", all_ids) == []


def test_organ_ids() -> None:
    """The ids of the standard organ, its boundary and the sinks."""
    assert ids.vascular_fraction(Organ.LIVER) == "f_vas_li"
    assert ids.interstitial_fraction(Organ.LIVER) == "f_int_li"
    assert ids.tissue_volume_total_ref("api") == "Vti_total_ref_api"
    assert ids.organ_outflow(Organ.LIVER, "api") == "Cout_li_api"
    assert ids.plasma_capacity("api") == "Vpl_api"
    assert ids.pool_capacity(Pool.VENOUS, "api") == "Vve_pl_api"
    assert ids.ps("api") == "PS_api"
    assert ids.PS_VAS == "PS_vas"
    assert ids.k_int_pls("api") == "K_int_pls_api"
    assert ids.k_cell_pls("api") == "K_cell_pls_api"
    assert ids.k_organ("api") == "K_or_api"
    assert ids.kbc_exchange("api") == "kex_bc_api"
    assert ids.lumen("api") == "Clumen_api"
    assert (ids.VBILE, ids.VURINE, ids.VLUMEN) == ("Vbile", "Vurine", "Vlumen")
    assert ids.inflow_reaction(Pool.PORTAL, Part.PLASMA, "api") == "J_in_po_pl_api"
    assert ids.outflow_reaction(Part.BLOOD_CELLS, "api") == "J_out_bc_api"
    assert ids.exchange_reaction(Part.TISSUE, "api") == "J_ti_api"


def test_zone_ids() -> None:
    """The ids of the parts of a zone carry the zone between part and substance."""
    z = ids.zone(1)
    assert z == "z1"
    assert ids.compartment(Part.TISSUE, z) == "Vti_z1"
    assert ids.species(Part.PLASMA, "api", z) == "Cpl_z1_api"
    assert ids.unbound(Part.TISSUE, "api", z) == "Cu_ti_z1_api"
    assert ids.plasma_capacity("api", z) == "Vpl_z1_api"
    assert ids.exchange_reaction(Part.TISSUE, "api", z) == "J_ti_z1_api"
    assert ids.sinusoid_reaction(Part.PLASMA, "api", z) == "J_sin_z1_pl_api"
    assert ids.segment_process("DRG2M", z) == "DRG2M_z1"
    assert ids.placement_weight("cyp3a4", z) == "w_cyp3a4_z1"
    assert ids.placement_profile("cyp3a4", z) == "p_cyp3a4_z1"
    assert ids.F_ZONE == "f_zone"


def test_species_of_zones_exact_substance() -> None:
    """Zone species belong to exactly their substance.

    A substance id does not start with a zone token (`z1_`, refused by the
    validation of the definition), so `Cpl_z1_drg` is a zone species of `drg`;
    a substance `z1` is allowed and has the species `Cpl_z1` and `Cpl_z2_z1`.
    """
    all_ids = ["LI__Cpl_z1", "LI__Cpl_z2_z1", "LI__Cpl_z1_drg", "LI__Cu_ti_z1_drg"]
    assert ids.species_of("z1", all_ids) == ["LI__Cpl_z1", "LI__Cpl_z2_z1"]
    assert ids.species_of("drg", all_ids) == ["LI__Cpl_z1_drg"]


def test_submodels_unique_and_not_physiology() -> None:
    """The submodel ids are unique and no id of the physiology (`CO`)."""
    submodels = [ids.submodel(o) for o in Organ]
    assert len(set(submodels)) == len(submodels)
    taken = {str(getattr(o, "sid", "")) for o in physiology_objects()}
    assert not taken & set(submodels)


def test_tubule_ids() -> None:
    """The ids of the tubule of the nephron kidney: the segment is always in the id."""
    assert ids.TUBULE_SEGMENTS == ("pt", "lh", "dt", "cd")
    assert ids.tubule("api", "pt") == "Ctub_pt_api"
    assert ids.tubule_volume("cd") == "Vtub_cd"
    assert ids.tubule_flow_fraction("lh") == "f_tub_lh"
    assert ids.tubule_reaction("api", "cd") == "J_tub_cd_api"
    assert ids.tubule_reabsorption_reaction("api") == "J_reab_api"


def test_species_of_tubule_exact_substance() -> None:
    """The tubule species match their substance only, also for ids which start with a segment token.

    The segment token is mandatory after `Ctub_`, so `Ctub_pt_pt_x` is the
    substance `pt_x` and never `x`; no substance token is reserved for it.
    """
    all_ids = [
        "KI__Ctub_pt_x",
        "KI__Ctub_cd_x",
        "KI__Ctub_pt_pt_x",
        "KI__Ctub_lh_x_m",
        "KI__Ctub_x",
        "KI__Ctub_ic_x",
        "KI__Vtub_pt",
    ]
    assert ids.species_of("x", all_ids) == ["KI__Ctub_pt_x", "KI__Ctub_cd_x"]
    assert ids.species_of("pt_x", all_ids) == ["KI__Ctub_pt_pt_x"]
    assert ids.species_of("x_m", all_ids) == ["KI__Ctub_lh_x_m"]
