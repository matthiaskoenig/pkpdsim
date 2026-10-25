"""Tests of the definition dataclasses and their validation."""

import dataclasses

import pytest

from pkpdsim.pbpk import (
    MM,
    PROTEINS,
    Binding,
    Direction,
    Erythrocytes,
    Excretion,
    FirstOrder,
    Kp,
    Membrane,
    Metabolism,
    Mode,
    Organ,
    PBPKDefinition,
    Route,
    Substance,
    Transport,
    load,
    sources,
)


def _api(**kwargs: object) -> Substance:
    """Get a valid substance."""
    values: dict[str, object] = {
        "sid": "api",
        "name": "apixaban",
        "mr": 459.5,
        "binding": Binding(fu=0.13, protein="albumin"),
        "distribution": Kp(2.0, ftissue=0.5),
        "excretion": Excretion(urine=True),
    }
    values.update(kwargs)
    return Substance(**values)  # ty: ignore[invalid-argument-type]


def _definition(**kwargs: object) -> PBPKDefinition:
    """Get a valid definition."""
    values: dict[str, object] = {
        "sid": "m",
        "substances": [_api(), Substance("api_m2", "metabolite", 445.5)],
        "chemistry": [
            Metabolism(
                "met",
                Organ.LIVER,
                "api",
                {"api_m2": 1.0},
                MM(vmax=1.0, km=0.1),
                protein="CYP3A4",
            )
        ],
        "routes": [Route.iv("api"), Route.oral("api"), Route.oral_solution("api")],
    }
    values.update(kwargs)
    return PBPKDefinition(**values)  # ty: ignore[invalid-argument-type]


def test_valid_definition_passes() -> None:
    """A valid definition validates."""
    _definition().validate()


def test_duplicate_substance_ids() -> None:
    """Duplicate substance ids are named."""
    with pytest.raises(ValueError, match="api_m2"):
        _definition(substances=[_api(sid="api_m2"), _api(sid="api_m2")]).validate()


@pytest.mark.parametrize("sid", ["BW", "CO", "HCT", "GFR", "Vli", "Qki", "Vblood"])
def test_substance_id_collides_with_physiology_id(sid: str) -> None:
    """A substance id which is a physiology id is rejected."""
    with pytest.raises(ValueError, match=f"'{sid}'"):
        _definition(substances=[_api(sid=sid)], chemistry=[], routes=[]).validate()


@pytest.mark.parametrize(
    "sid", ["1abc", "a-b", "a b", "a__b", "pl_x", "bc_x", "free_x", "blood_x", "ve_x"]
)
def test_invalid_sid(sid: str) -> None:
    """Ids which are no SId, contain '__' or take a token of the id scheme."""
    with pytest.raises(ValueError, match=sid):
        _definition(substances=[_api(sid=sid)], chemistry=[], routes=[]).validate()


def test_unknown_substance_in_chemistry() -> None:
    """A process on an unknown substance names it."""
    t = Transport(
        "t",
        Organ.LIVER,
        "nope",
        Direction.INFLUX,
        Membrane.BASOLATERAL,
        FirstOrder(1.0),
    )
    with pytest.raises(ValueError, match="nope"):
        _definition(chemistry=[t]).validate()
    m = Metabolism("m", Organ.LIVER, "api", {"ghost": 1.0}, FirstOrder(1.0))
    with pytest.raises(ValueError, match="ghost"):
        _definition(chemistry=[m]).validate()


def test_unknown_protein() -> None:
    """An unknown protein is named."""
    m = Metabolism(
        "m", Organ.LIVER, "api", {"api_m2": 1.0}, FirstOrder(1.0), protein="CYP9Z9"
    )
    with pytest.raises(ValueError, match="CYP9Z9"):
        _definition(chemistry=[m]).validate()


def test_product_fractions_sum_to_one() -> None:
    """Product fractions which do not sum to one are rejected."""
    m = Metabolism("m", Organ.LIVER, "api", {"api_m2": 0.7}, FirstOrder(1.0))
    with pytest.raises(ValueError, match=r"sum to 0\.7"):
        _definition(chemistry=[m]).validate()


def test_perfusion_with_kinetic_erythrocytes_rejected() -> None:
    """Kinetic erythrocytes need an exchange mode."""
    s = _api(
        erythrocytes=Erythrocytes("kinetic", kbc=1.0, k_exchange=0.1),
        distribution=Kp(2.0, mode=Mode.PERFUSION),
    )
    with pytest.raises(ValueError, match="api: kinetic erythrocytes"):
        _definition(substances=[s], chemistry=[], routes=[]).validate()


def test_erythrocyte_parameters() -> None:
    """Equilibrium takes exactly one of bp and kbc, kinetic kbc and k_exchange."""
    for ery in (
        Erythrocytes("equilibrium"),
        Erythrocytes("equilibrium", bp=1.0, kbc=1.0),
        Erythrocytes("kinetic", kbc=1.0),
    ):
        with pytest.raises(ValueError, match="erythrocytes"):
            _definition(
                substances=[_api(erythrocytes=ery)], chemistry=[], routes=[]
            ).validate()
    ok = [
        Erythrocytes("equilibrium", bp=1.0),
        Erythrocytes("equilibrium", kbc=1.0),
        Erythrocytes("kinetic", kbc=1.0, k_exchange=0.1),
    ]
    for ery in ok:
        _definition(
            substances=[_api(erythrocytes=ery)], chemistry=[], routes=[]
        ).validate()


@pytest.mark.parametrize("mode", [Mode.TISSUE, Mode.PERMEABILITY])
def test_tissue_mode_requires_ftissue(mode: Mode) -> None:
    """The exchange modes need the whole-body exchange."""
    s = _api(distribution=Kp(2.0, mode=mode))
    with pytest.raises(ValueError, match="requires ftissue"):
        _definition(substances=[s], chemistry=[], routes=[]).validate()


def test_override_for_organ_not_entered() -> None:
    """An override or exchange of an organ outside `organs` is rejected."""
    s = _api(
        distribution=Kp(
            2.0,
            overrides={Organ.BRAIN: 0.1},
            organs=frozenset({Organ.LIVER}),
            ftissue=0.5,
        )
    )
    with pytest.raises(ValueError, match=r"override.*'br'"):
        _definition(substances=[s], chemistry=[], routes=[]).validate()


def test_transport_in_organ_without_tissue() -> None:
    """A transport in an organ the substance does not enter is valid."""
    s = _api(distribution=Kp(2.0, organs=frozenset({Organ.LIVER}), ftissue=0.5))
    t = Transport(
        "t",
        Organ.KIDNEY,
        "api",
        Direction.EFFLUX,
        Membrane.APICAL,
        MM(1.0, 0.1),
        protein="ABCB1",
    )
    _definition(substances=[s], chemistry=[t], routes=[]).validate()


def test_unknown_route_substance() -> None:
    """A route of an unknown substance is rejected."""
    with pytest.raises(ValueError, match="zzz"):
        _definition(routes=[Route.iv("zzz")]).validate()


def test_all_problems_reported() -> None:
    """The message lists every problem."""
    with pytest.raises(ValueError, match=r"(?s)BW.*zzz"):
        _definition(
            substances=[_api(sid="BW")], chemistry=[], routes=[Route.iv("zzz")]
        ).validate()


def test_dataclasses_are_frozen() -> None:
    """Definitions cannot be changed."""
    with pytest.raises(dataclasses.FrozenInstanceError):
        _api().sid = "x"  # ty: ignore[invalid-assignment]


def test_proteins_table() -> None:
    """The table has the required proteins with verified-form accessions."""
    required = [
        "CYP1A2",
        "CYP2C8",
        "CYP2C9",
        "CYP2C19",
        "CYP2D6",
        "CYP3A4",
        "UGT1A1",
        "UGT1A9",
        "UGT2B7",
        "SULT1A1",
        "SLCO1B1",
        "SLCO1B3",
        "SLC22A1",
        "SLC22A2",
        "SLC22A6",
        "SLC22A8",
        "SLC47A1",
        "ABCB1",
        "ABCG2",
        "ABCC2",
    ]
    assert set(required) <= set(PROTEINS)
    assert PROTEINS["CYP3A4"].uniprot == "P08684"
    assert PROTEINS["SLC22A6"].ncit is None
    assert load("proteins")["source"] in sources()
