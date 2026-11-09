"""The seams of the engine for an organ with its own boundary, fractions and cell model.

The perfusion organs call the engine with their own fractions of the organ
(`ZonationSettings.fractions`), the circuit as the boundary of the chain
(`FlowBoundary.EXTERNAL`), the cell species in the cell model
(`cell_species=False`) and the processes of the cell model which a profile
scales (`cell_processes`).
"""

from dataclasses import dataclass
from pathlib import Path

import libsbml
import pytest
import roadrunner
from sbmlutils.factory import PORT_SUFFIX, Model, Parameter, create_model

from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Value
from pkpdsim.perfusion.units import PU
from pkpdsim.zonation import ZonationBuilder, ids
from pkpdsim.zonation.composition import CellPorts, compose, write
from pkpdsim.zonation.flow import ImposedFlow
from pkpdsim.zonation.geometry import Layer, geometry
from pkpdsim.zonation.profiles import Weights
from pkpdsim.zonation.settings import FlowBoundary, Geometry, ZonationSettings
from pkpdsim.zonation.transport import Solute

from .conftest import mass_g
from .data.toy_cell import TOY_ID, UPTAKE, create_toy_cell, toy_cell_model


@dataclass(frozen=True)
class _Gas:
    sid: str = "o2"
    name: str = "oxygen"
    mr: None = None
    diffusion: None = None


DRUG = _Gas("drg", "drug")


def _settings(**kwargs: object) -> ZonationSettings:
    base: dict[str, object] = {
        "n_zones": 4,
        "geometry": Geometry.WELL_STIRRED,
        "bile": None,
        "flow": ImposedFlow(Value(0.02, U.l_per_min)),
        "liver_mass": mass_g(10.0),
    }
    return ZonationSettings(**{**base, **kwargs})  # ty: ignore[invalid-argument-type]


def _builder(
    solutes: tuple[Solute, ...] = (_Gas(), DRUG),
    processes: tuple[str, ...] = (),
    **kwargs: object,
) -> ZonationBuilder:
    return ZonationBuilder(_settings(**kwargs), solutes, None, cell_processes=processes)


def _constants(model: Model) -> dict[str, float]:
    """The values of the parameters of a model by id."""
    return {str(p.sid): float(p.value) for p in model.parameters}  # ty: ignore[invalid-argument-type]


# fractions of the organ


def test_the_fractions_of_the_organ_set_the_volumes() -> None:
    """A vascular fraction 0.15 and a cell fraction 0.6 of 10 g at 1.1 kg/l, 4 zones."""
    geo = geometry(_settings().species)
    import dataclasses

    own = dataclasses.replace(geo, f_par=1.0, density=1.1)
    builder = _builder(
        fractions={Layer.PLASMA: 0.15, Layer.CELL: 0.6}, geometry_constants=own
    )
    model = builder.build()
    volume = {str(c.sid): c.value for c in model.compartments}
    litre = 10e-3 / 1.1
    assert volume[ids.volume(ids.VASCULAR, "z2")] == pytest.approx(0.15 * litre / 4)
    assert volume[ids.volume(Layer.CELL, "z2")] == pytest.approx(0.6 * litre / 4)
    constants = _constants(model)
    assert constants[ids.fraction(ids.VASCULAR)] == pytest.approx(0.15)
    assert constants[ids.fraction(Layer.CELL)] == pytest.approx(0.6)


def test_the_volumes_follow_a_change_of_the_mass(tmp_path: Path) -> None:
    """The initial assignments carry the organ fractions, so the mass sets the zones."""
    import dataclasses

    own = dataclasses.replace(geometry(_settings().species), f_par=1.0, density=1.0)
    model = _builder(
        fractions={Layer.PLASMA: 0.2, Layer.CELL: 0.5}, geometry_constants=own
    ).build()
    path = tmp_path / "m.xml"
    create_model(model, filepath=path, sbml_level=3, sbml_version=2, validate=False)
    r = roadrunner.RoadRunner(str(path))
    r[f"init({ids.MASS})"] = 0.02
    r.reset()
    assert r[ids.volume(ids.VASCULAR, "z1")] == pytest.approx(0.2 * 0.02 / 4)


def test_without_fractions_the_table_gives_them() -> None:
    constants = _constants(_builder().build())
    assert 0 < constants[ids.fraction(Layer.CELL)] < 1
    assert constants[ids.fraction(ids.VASCULAR)] != 0.15


# the boundary of the circuit


def test_an_external_boundary_has_neither_inflow_nor_outflow() -> None:
    model = _builder(boundary=FlowBoundary.EXTERNAL).build()
    reactions = {str(r.sid) for r in model.reactions}
    assert reactions == {f"FLOW_o2_vas_z{k}_vas_z{k + 1}" for k in (1, 2, 3)} | {
        f"FLOW_drg_vas_z{k}_vas_z{k + 1}" for k in (1, 2, 3)
    }
    assert not [p for p in model.parameters if str(p.sid).endswith("_in")]


def test_an_open_boundary_has_both() -> None:
    reactions = {str(r.sid) for r in _builder().build().reactions}
    assert "FLOW_o2_vas_in_vas_z1" in reactions
    assert "FLOW_o2_vas_z4_vas_out" in reactions


# the cell species


def test_the_cell_model_holds_the_cell_species() -> None:
    species = {str(s.sid) for s in _builder(cell_species=False).build().species}
    assert species == {f"{q}_vas_z{k}" for q in ("o2", "drg") for k in (1, 2, 3, 4)}


def test_the_cell_species_are_the_default() -> None:
    species = {str(s.sid) for s in _builder().build().species}
    assert "o2_cell_z1" in species


def test_a_sinusoid_needs_the_cell_species() -> None:
    with pytest.raises(ValueError, match="cell_species"):
        _settings(geometry=Geometry.SINUSOID, cell_species=False)


# the processes of the cell model


def test_a_cell_process_has_a_zonal_factor_which_normalizes_to_the_mean() -> None:
    model = _builder(
        processes=("OXPHOS",), profiles={"OXPHOS": Weights((1, 2, 3, 6))}
    ).build()
    p = _constants(model)
    factors = [p[ids.profile("OXPHOS", ids.zone(k))] for k in (1, 2, 3, 4)]
    assert factors == pytest.approx([4 * w / 12 for w in (1, 2, 3, 6)])
    assert sum(factors) == pytest.approx(4.0)


def test_a_profile_of_a_process_nothing_has_is_refused() -> None:
    with pytest.raises(ValueError, match="OXPHOS"):
        _builder(profiles={"OXPHOS": Weights((1, 2, 3, 6))})


# composition


def _toy(processes: tuple[str, ...] = ()) -> Model:
    return toy_cell_model(processes=processes)


def test_the_zonal_factors_replace_the_profile_ports() -> None:
    builder = _builder(
        processes=("OXPHOS",), profiles={"OXPHOS": Weights((1, 2, 3, 6))}
    )
    m = compose(builder, TOY_ID, CellPorts(), exchange=("drg",))
    found = {
        (r.submodelRef, r.elementRef)
        for r in m.replaced_elements
        if r.portRef == f"p_OXPHOS{PORT_SUFFIX}"
    }
    assert found == {(f"M0{k}", f"p_OXPHOS_z{k}") for k in (1, 2, 3, 4)}


def test_a_process_without_a_profile_has_no_port() -> None:
    m = compose(_builder(processes=("OXPHOS",)), TOY_ID, CellPorts(), exchange=("drg",))
    assert not [r for r in m.replaced_elements if "OXPHOS" in str(r.portRef)]


def test_the_zonal_factor_scales_the_rate_of_its_zone(tmp_path: Path) -> None:
    builder = _builder(processes=("UP",), profiles={"UP": Weights((1, 2, 3, 6))})
    flat = write(
        compose(builder, TOY_ID, CellPorts(), exchange=("drg",)),
        tmp_path,
        cell_model=toy_cell_model(processes=("UP",)),
    )
    r = roadrunner.RoadRunner(str(flat))
    for k in (1, 2, 3, 4):
        r[f"[{ids.species('drg', ids.VASCULAR, ids.zone(k))}]"] = 1.0
    for k, w in zip((1, 2, 3, 4), (1, 2, 3, 6), strict=True):
        assert r[f"M0{k}__UP"] == pytest.approx(
            UPTAKE * (4 * w / 12) * r[ids.volume(ids.VASCULAR, ids.zone(k))]
        )


def test_only_the_exchanged_solutes_are_wired() -> None:
    builder = _builder(solutes=(DRUG, _Gas("tr", "tracer")))
    m = compose(builder, TOY_ID, CellPorts(), exchange=("drg",))
    ports = {r.portRef for r in m.replaced_elements if r.submodelRef == "M01"}
    assert f"drg_ext{PORT_SUFFIX}" in ports
    assert f"tr_ext{PORT_SUFFIX}" not in ports


def test_a_given_temperature_parameter_replaces_the_default() -> None:
    mine = Parameter(
        "T",
        25.0,
        unit=PU.degC,
        name="temperature of the circuit",
        constant=False,
    )
    m = compose(_builder(), TOY_ID, CellPorts(), exchange=("drg",), temperature=mine)
    (found,) = [p for p in m.parameters if p.sid in ("T", ids.TEMPERATURE)]
    assert found is mine
    wired = {
        r.elementRef
        for r in m.replaced_elements
        if r.portRef == f"temperature{PORT_SUFFIX}"
    }
    assert wired == {"T"}


# the file of the cell model


def test_write_writes_the_cell_model_beside_the_composed_model(tmp_path: Path) -> None:
    flat = write(
        compose(_builder(), TOY_ID, CellPorts(), exchange=("drg",)),
        tmp_path,
        cell_model=_toy(),
    )
    assert (tmp_path / f"{TOY_ID}.xml").is_file()
    assert flat == tmp_path / f"{_builder().build().sid}_flat.xml"
    model = libsbml.readSBMLFromFile(str(flat)).getModel()
    assert model.getReaction("M01__UP") is not None


def test_write_writes_the_same_cell_model_again(tmp_path: Path) -> None:
    for _ in range(2):
        write(
            compose(_builder(), TOY_ID, CellPorts(), exchange=("drg",)),
            tmp_path,
            cell_model=_toy(),
        )


def test_write_replaces_an_older_version_of_the_cell_model(tmp_path: Path) -> None:
    """A rebuilt model overwrites its own file: the version with a process is the new one."""
    create_toy_cell(tmp_path)
    write(
        compose(_builder(), TOY_ID, CellPorts(), exchange=("drg",)),
        tmp_path,
        cell_model=_toy(),
    )
    older = (tmp_path / f"{TOY_ID}.xml").read_text()
    newer = toy_cell_model(processes=("UP",))
    write(
        compose(_builder(), TOY_ID, CellPorts(), exchange=("drg",)),
        tmp_path,
        cell_model=newer,
    )
    assert (tmp_path / f"{TOY_ID}.xml").read_text() != older
    assert "p_UP" in (tmp_path / f"{TOY_ID}.xml").read_text()


def test_write_refuses_to_overwrite_another_model(tmp_path: Path) -> None:
    other = toy_cell_model()
    other.sid = "something_else"
    create_model(
        other,
        filepath=tmp_path / f"{TOY_ID}.xml",
        sbml_level=3,
        sbml_version=2,
        validate=False,
    )
    before = (tmp_path / f"{TOY_ID}.xml").read_text()
    with pytest.raises(FileExistsError, match="something_else"):
        write(
            compose(_builder(), TOY_ID, CellPorts(), exchange=("drg",)),
            tmp_path,
            cell_model=_toy(),
        )
    assert (tmp_path / f"{TOY_ID}.xml").read_text() == before


def test_write_refuses_a_cell_model_of_another_id(tmp_path: Path) -> None:
    other = _toy()
    other.sid = "other"
    with pytest.raises(ValueError, match="other"):
        write(
            compose(_builder(), TOY_ID, CellPorts(), exchange=("drg",)),
            tmp_path,
            cell_model=other,
        )


# the constants the zones share


def test_a_shared_constant_is_one_parameter_of_the_liver(tmp_path: Path) -> None:
    """The constant of the cell model is a port, which the liver replaces in
    every zone by its own: changing it changes every zone."""
    model = compose(
        _builder(), TOY_ID, CellPorts(), exchange=("drg",), shared=("k_up",)
    )
    twice = Parameter(
        "k_up",
        2 * UPTAKE,
        unit=U.per_min,
        name="uptake rate of all zones",
        constant=True,
    )
    model.parameters.append(twice)
    flat = write(model, tmp_path, cell_model=toy_cell_model(shared=True))
    doc = libsbml.readSBMLFromFile(str(flat)).getModel()
    assert not [
        p.getId() for p in doc.getListOfParameters() if p.getId().endswith("__k_up")
    ]
    r = roadrunner.RoadRunner(str(flat))
    for k in (1, 2, 3, 4):
        r[f"[{ids.species('drg', ids.VASCULAR, ids.zone(k))}]"] = 1.0
    volume = {k: r[ids.volume(ids.VASCULAR, ids.zone(k))] for k in (1, 4)}
    for k in (1, 4):
        assert r[f"M0{k}__UP"] == pytest.approx(2 * UPTAKE * volume[k])
    r["k_up"] = 3 * UPTAKE
    assert r["M04__UP"] == pytest.approx(3 * UPTAKE * volume[4])


def test_a_shared_constant_is_wired_in_every_zone() -> None:
    m = compose(_builder(), TOY_ID, CellPorts(), exchange=("drg",), shared=("k_up",))
    found = {
        (r.submodelRef, r.elementRef)
        for r in m.replaced_elements
        if r.portRef == f"k_up{PORT_SUFFIX}"
    }
    assert found == {(f"M0{k}", "k_up") for k in (1, 2, 3, 4)}
