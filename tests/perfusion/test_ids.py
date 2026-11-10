"""The ids of the perfusion library which the zonation engine names.

The vascular space, the chain and the cells of a zone of a liver are the
engine's (`pkpdsim.zonation.ids`); the ids of the perfusion library are the names
of the flat model, for the cells of a zone as well (the liver exposes the states
of its cell submodels under them).
"""

import pytest

from pkpdsim.perfusion import ids as perf_ids
from pkpdsim.zonation import ids
from pkpdsim.zonation.geometry import Layer


def test_the_perfusion_ids_of_a_zone_are_those_of_the_engine() -> None:
    """The vascular space and the cells of a zone are named by the engine."""
    assert perf_ids.vascular("z3") == ids.volume(ids.VASCULAR, "z3") == "V_vas_z3"
    assert perf_ids.cells("z3") == ids.volume(Layer.CELL, "z3")
    assert perf_ids.species("glc", "z3") == ids.species("glc", ids.VASCULAR, "z3")
    assert perf_ids.flow("glc", "z2", "z3") == "FLOW_glc_vas_z2_vas_z3"
    assert perf_ids.flow("glc", "in", "z1") == "FLOW_glc_in_vas_z1"
    assert perf_ids.vascular_fraction("li") == ids.fraction(ids.VASCULAR)
    assert perf_ids.cell_fraction("li") == ids.fraction(Layer.CELL)
    assert perf_ids.vascular_fraction("ki") == "f_vas_ki"


@pytest.mark.parametrize(
    ("made", "expected"),
    [
        (perf_ids.species("glc", "z3_cell"), "glc_z3_cell"),
        (perf_ids.species("alt", "z12_dead"), "alt_z12_dead"),
        (perf_ids.reaction("GLC2PYR", "z3"), "GLC2PYR_z3"),
        (perf_ids.buffer("z3"), "BUF_z3"),
        (perf_ids.ph("z3"), "pH_z3"),
        (perf_ids.alive("z3"), "f_alive_z3"),
        (perf_ids.necrosis_rate("z3"), "v_nec_z3"),
        (perf_ids.bile("z3"), "Vbile_z3"),
        (perf_ids.profile("OXPHOS", "z3"), "p_OXPHOS_z3"),
        # the sites which are no zone are not changed
        (perf_ids.species("glc", "res"), "glc_res"),
        (perf_ids.species("glc", "pp"), "glc_pp"),
        (perf_ids.alive("pp"), "f_alive_pp"),
        (perf_ids.ph("pc"), "pH_pc"),
        (perf_ids.buffer("ext"), "BUF_ext"),
        (perf_ids.reaction("GLC2PYR", "ext"), "GLC2PYR_ext"),
    ],
)
def test_the_perfusion_ids_of_the_cells_are_those_of_the_flat_model(
    made: str, expected: str
) -> None:
    """The ids of a zone are the ids of the flat model of every organ: the cells of
    the liver are a submodel, whose states the liver exposes under these names, the
    cells of the kidney and a generic organ are emitted by the builder."""
    assert made == expected
