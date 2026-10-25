"""Tests of the collapsing of the parts of a substance in an organ."""

import dataclasses

from pkpdsim.pbpk import (
    MM,
    Direction,
    Erythrocytes,
    Kp,
    Membrane,
    Metabolism,
    Mode,
    Organ,
    Parts,
    Transport,
    parts_for,
)

from .conftest import definition, simple_drug

LIVER = Organ.LIVER


def _transport(
    organ: Organ = LIVER,
    substance: str = "drg",
    membrane: Membrane = Membrane.BASOLATERAL,
) -> Transport:
    """Get a transport."""
    return Transport(
        sid="T1",
        organ=organ,
        substance=substance,
        direction=Direction.INFLUX,
        membrane=membrane,
        kinetics=MM(vmax=1e-3, km=0.01),
        protein="OATP1B1",
    )


def test_definition_validates() -> None:
    """The definition of the helpers is valid."""
    definition().validate()


def test_no_erythrocytes() -> None:
    """Without erythrocytes the flows carry plasma only."""
    parts = parts_for(simple_drug(), LIVER, [])
    assert not parts.blood_cells
    assert parts.flow_carrier == "plasma"


def test_erythrocytes_equilibrium() -> None:
    """Erythrocytes in equilibrium need no blood cell species."""
    s = dataclasses.replace(
        simple_drug(), erythrocytes=Erythrocytes(mode="equilibrium", bp=1.2)
    )
    parts = parts_for(s, LIVER, [])
    assert not parts.blood_cells
    assert parts.flow_carrier == "blood_bp"


def test_erythrocytes_kinetic() -> None:
    """Kinetic erythrocytes build a blood cell species, flows carry plasma."""
    s = dataclasses.replace(
        simple_drug(),
        erythrocytes=Erythrocytes(mode="kinetic", kbc=0.5, k_exchange=1.0),
    )
    parts = parts_for(s, LIVER, [])
    assert parts.blood_cells
    assert parts.flow_carrier == "plasma"


def test_perfusion_mode() -> None:
    """Flow-limited distribution is one species for the organ."""
    s = dataclasses.replace(simple_drug(), distribution=Kp(2.0, mode=Mode.PERFUSION))
    assert parts_for(s, LIVER, []).tissue == "organ"


def test_merged_default() -> None:
    """A Kp gives the merged tissue by default."""
    assert parts_for(simple_drug(), LIVER, []).tissue == "merged"


def test_permeability_mode() -> None:
    """Permeability-limited distribution builds the full form."""
    s = dataclasses.replace(
        simple_drug(), distribution=Kp(2.0, mode=Mode.PERMEABILITY, ftissue=1.0)
    )
    assert parts_for(s, LIVER, []).tissue == "full"


def test_transport_forces_full() -> None:
    """A transport of the substance in the organ forces the full form."""
    parts = parts_for(simple_drug(), LIVER, [_transport()])
    assert parts == Parts(blood_cells=False, flow_carrier="plasma", tissue="full")
    perfusion = dataclasses.replace(
        simple_drug(), distribution=Kp(2.0, mode=Mode.PERFUSION)
    )
    assert parts_for(perfusion, LIVER, [_transport()]).tissue == "full"


def test_transport_of_other_substance_or_organ_does_not_force_full() -> None:
    """Only a transport of this substance in this organ counts."""
    chemistry = [_transport(substance="other"), _transport(organ=Organ.KIDNEY)]
    assert parts_for(simple_drug(), LIVER, chemistry).tissue == "merged"


def test_apical_process_forces_full() -> None:
    """An apical process of the substance in the organ forces the full form."""
    chemistry = [_transport(membrane=Membrane.APICAL)]
    assert parts_for(simple_drug(), LIVER, chemistry).tissue == "full"


def test_metabolism_does_not_force_full() -> None:
    """Metabolism acts in the plasma part and does not change the collapsing."""
    chemistry = definition().chemistry
    assert isinstance(chemistry[0], Metabolism)
    assert parts_for(simple_drug(), LIVER, chemistry).tissue == "merged"


def test_not_entered_organ_has_no_tissue() -> None:
    """An organ which the substance does not enter has the vascular parts only."""
    s = dataclasses.replace(
        simple_drug(),
        distribution=Kp(2.0, ftissue=1.0, organs=frozenset({Organ.LIVER})),
    )
    assert parts_for(s, LIVER, []).tissue == "merged"
    assert parts_for(s, Organ.KIDNEY, []).tissue == "none"


def test_no_distribution_has_no_tissue() -> None:
    """A substance without distribution is vascular only."""
    s = dataclasses.replace(simple_drug(), distribution=None)
    assert parts_for(s, LIVER, []).tissue == "none"


def test_transport_without_distribution_is_full() -> None:
    """A transport forces the full form without a distribution, and where not entered."""
    s = dataclasses.replace(simple_drug(), distribution=None)
    assert parts_for(s, LIVER, [_transport()]).tissue == "full"
    entered = dataclasses.replace(
        simple_drug(),
        distribution=Kp(2.0, ftissue=1.0, organs=frozenset({Organ.HEART})),
    )
    assert parts_for(entered, LIVER, [_transport()]).tissue == "full"
