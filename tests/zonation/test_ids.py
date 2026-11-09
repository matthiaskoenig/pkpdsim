"""The id scheme of the zonated liver."""

import pytest

from pkpdsim.pbpk import ids as pbpk_ids
from pkpdsim.perfusion import ids as perf_ids
from pkpdsim.zonation import ids
from pkpdsim.zonation.geometry import Layer


def test_zone_ids_agree_with_both_existing_conventions() -> None:
    """The zone token is the one of the PBPK library and of the perfusion."""
    assert ids.zone(3) == pbpk_ids.zone(3) == perf_ids.zone_id(3) == "z3"


def test_profile_id_matches_the_perfusion_pattern() -> None:
    """A profile is the factor `p_<process>_<zone>` of the perfusion."""
    assert ids.profile("cyp2e1", ids.zone(1)) == "p_cyp2e1_z1"
    assert ids.profile("cyp2e1", ids.zone(1)) == perf_ids.profile("cyp2e1", "z1")


def test_the_ids_of_a_layer() -> None:
    """The site, the volume and the species of a layer of a zone."""
    assert ids.site(Layer.DISSE, "z3") == "disse_z3"
    assert ids.volume(Layer.DISSE, "z3") == "V_disse_z3"
    assert ids.species("glc", Layer.HEPATOCYTE, "z2") == "glc_hepatocyte_z2"


def test_the_submodel_id_is_zero_padded() -> None:
    """The submodel of zone k is `M<kk>`."""
    assert ids.submodel(3) == "M03"
    assert ids.submodel(12) == "M12"


def test_a_zone_below_one_raises() -> None:
    """Zones are numbered from 1."""
    with pytest.raises(ValueError, match="from 1"):
        ids.zone(0)
    with pytest.raises(ValueError, match="from 1"):
        ids.submodel(0)


def test_the_flow_and_diffusion_ids() -> None:
    """The ids of the chain."""
    assert ids.layer_flow(Layer.PLASMA) == "Q_plasma_li"
    assert ids.axial("glc", "disse_z1", "disse_z2") == "D_glc_disse_z1_disse_z2"


def test_a_flat_id_is_prefixed_with_its_submodel() -> None:
    """What the flattening of a comp model makes of the id of a submodel element."""
    assert ids.flat(3, "glc_cell") == "M03__glc_cell"
    assert ids.flat(12, "x") == "M12__x"
