"""Tests of the id scheme."""

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.organs import Organ, Part, Pool


def test_examples() -> None:
    """The documented examples."""
    assert ids.submodel(Organ.LIVER) == "LI"
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
        "Cve_api",
        "Cve_bc_api",
        "Car_bc_api",
        "LI__Cpl_api_m2",
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
        "Cve_api",
        "Cve_bc_api",
        "Car_bc_api",
    ]
    assert ids.species_of("api_m2", all_ids) == ["LI__Cpl_api_m2", "Cve_api_m2"]
    assert ids.species_of("drg", all_ids) == ["Cpl_drg"]
    assert ids.species_of("m2", all_ids) == []
