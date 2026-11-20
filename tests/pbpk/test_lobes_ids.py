"""Tests of the lobe token in the id scheme of the liver."""

import pytest

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.organs import Organ, Part, Pool


def test_site_composes_lobe_and_zone() -> None:
    assert ids.site(None, None) is None
    assert ids.site("lll", None) == "lll"
    assert ids.site(None, "z3") == "z3"
    assert ids.site("lll", "z3") == "lll_z3"


def test_compartment_and_species_take_the_composite_site() -> None:
    assert ids.compartment(Part.PLASMA, ids.site("lll", "z3")) == "Vpl_lll_z3"
    assert ids.species(Part.PLASMA, "apap", ids.site("lll", "z3")) == "Cpl_lll_z3_apap"
    assert ids.compartment(Part.TISSUE, ids.site("sg8", None)) == "Vti_sg8"


def test_lobe_anatomy_ids() -> None:
    assert ids.lobe_volume("lll") == "Vli_lll"
    assert ids.lobe_volume_fraction("lll") == "f_vol_lll"
    assert ids.lobe_present("lll") == "f_present_lll"


def test_lobe_flow_ids() -> None:
    assert ids.lobe_flow(ids.QHA, "lll") == "Qha_lll"
    assert ids.lobe_flow(ids.QPO, "lll") == "Qpo_lll"
    assert ids.lobe_flow(ids.QH, "lll") == "Qh_lll"
    assert ids.lobe_perfusion("lll") == "q_lll"


def test_lobe_conductance_ids() -> None:
    assert ids.arterial_conductance("lll") == "g_ha_lll"
    assert ids.portal_conductance("lll") == "g_po_lll"


def test_the_whole_liver_ids_of_the_lobes() -> None:
    assert ids.VLI_REF == "Vli_ref"
    assert ids.F_REMNANT == "f_remnant"


def test_reaction_ids_without_a_site_are_unchanged() -> None:
    assert ids.inflow_reaction(Pool.ARTERIAL, Part.PLASMA, "apap") == "J_in_ar_pl_apap"
    assert ids.outflow_reaction(Part.PLASMA, "apap") == "J_out_pl_apap"
    assert ids.shunt_reaction(Pool.PORTAL, Part.PLASMA, "apap") == "J_shunt_po_pl_apap"
    assert ids.outflow_conc("apap") == "Cout_apap"
    assert ids.outflow_conc_bc("apap") == "Cout_bc_apap"


def test_reaction_ids_with_a_site_carry_it() -> None:
    site = ids.site("lll", "z1")
    assert (
        ids.inflow_reaction(Pool.ARTERIAL, Part.PLASMA, "apap", site)
        == "J_in_ar_pl_apap_lll_z1"
    )
    assert ids.outflow_reaction(Part.PLASMA, "apap", site) == "J_out_pl_apap_lll_z1"
    assert (
        ids.shunt_reaction(Pool.PORTAL, Part.PLASMA, "apap", site)
        == "J_shunt_po_pl_apap_lll_z1"
    )


def test_the_outflow_concentration_of_a_lobe() -> None:
    assert ids.outflow_conc("apap", "lll") == "Cout_lll_apap"
    assert ids.outflow_conc_bc("apap", "lll") == "Cout_bc_lll_apap"


def test_the_zonal_mean_of_a_lobe() -> None:
    assert ids.zonal_mean(Part.PLASMA, "apap", True) == "Cpl_pp_apap"
    assert ids.zonal_mean(Part.PLASMA, "apap", True, "lll") == "Cpl_lll_pp_apap"
    assert ids.zonal_mean(Part.PLASMA, "apap", False, "lll") == "Cpl_lll_pc_apap"


def test_every_lobe_id_is_a_valid_sid() -> None:
    made = [
        ids.lobe_volume("sg8"),
        ids.lobe_volume_fraction("sg8"),
        ids.lobe_present("sg8"),
        ids.lobe_flow(ids.QH, "sg8"),
        ids.lobe_perfusion("sg8"),
        ids.arterial_conductance("sg8"),
        ids.portal_conductance("sg8"),
        ids.outflow_conc("apap", "sg8"),
        ids.inflow_reaction(Pool.ARTERIAL, Part.PLASMA, "apap", ids.site("sg8", "z2")),
    ]
    for sid in made:
        assert sid.isidentifier(), sid


def test_the_flat_id_of_a_lobe_element() -> None:
    assert ids.flat(Organ.LIVER, ids.lobe_volume("lll")) == "LI__Vli_lll"


@pytest.mark.parametrize("flow", ["Qli", "nonsense"])
def test_lobe_flow_refuses_a_flow_which_is_no_hepatic_flow(flow: str) -> None:
    with pytest.raises(ValueError, match=flow):
        ids.lobe_flow(flow, "lll")


def test_lobe_supply_fraction_ids() -> None:
    assert ids.arterial_fraction("lll") == "f_ha_lll"
    assert ids.portal_fraction("lll") == "f_po_lll"


def test_lobe_remodelling_ids() -> None:
    assert ids.lobe_reference_fraction("lll") == "f_vol_ref_lll"
    assert ids.lobe_growth("lll") == "dVli_lll"
    assert ids.VLI_LOBES == "Vli_lobes"


def test_remodelling_reaction_id() -> None:
    assert (
        ids.remodelling_reaction(Part.TISSUE, "apap", "lll") == "J_remodel_lll_ti_apap"
    )
