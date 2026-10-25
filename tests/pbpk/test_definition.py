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
    Hill,
    Kp,
    Membrane,
    Metabolism,
    Mode,
    Observer,
    Organ,
    PBPKDefinition,
    Route,
    Substance,
    Transport,
    U,
    Value,
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
    "sid",
    [
        "1abc",
        "a-b",
        "a b",
        "a__b",
        "pl_x",
        "bc_x",
        "free_x",
        "blood_x",
        "ve_x",
        "du_x",
        "z1_x",
        "z12_x",
    ],
)
def test_invalid_sid(sid: str) -> None:
    """Ids which are no SId, contain '__' or take a token of the id scheme.

    The tokens of the parts, pools, observers, lumen segments (`du_`) and the
    zones of a zonated organ (`z1_`), which precede the substance in an id.
    """
    with pytest.raises(ValueError, match=sid):
        _definition(substances=[_api(sid=sid)], chemistry=[], routes=[]).validate()


@pytest.mark.parametrize("sid", ["z1", "zx_drg", "zinc", "drg_z1"])
def test_sid_like_zone_token_valid(sid: str) -> None:
    """Only a zone token followed by `_` at the start of a substance id is reserved."""
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


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"binding": Binding(fu=0.0)}, "api: fu"),
        ({"binding": Binding(fu=1.5)}, "api: fu"),
        ({"binding": Binding(fu=Value(-1.0, U.dimensionless))}, "api: fu"),
        ({"mr": 0.0}, "api: mr"),
        ({"distribution": Kp(-1.0, ftissue=0.5)}, "api: Kp"),
        (
            {"distribution": Kp(1.0, overrides={Organ.LIVER: 0.0}, ftissue=0.5)},
            "Kp override 'li'",
        ),
        ({"distribution": Kp(1.0, ftissue=0.0)}, "api: ftissue"),
        ({"distribution": Kp(1.0, ftissue=0.5, ps={Organ.LIVER: -1.0})}, "ps 'li'"),
        (
            {"excretion": Excretion(urine=True, fraction_reabsorbed=1.2)},
            "fraction_reabsorbed",
        ),
        ({"erythrocytes": Erythrocytes("equilibrium", bp=0.0)}, "api: bp"),
        ({"erythrocytes": Erythrocytes("equilibrium", kbc=-1.0)}, "api: kbc"),
    ],
)
def test_substance_numbers_checked(kwargs: dict[str, object], message: str) -> None:
    """Plain numbers and the numbers of values are range checked."""
    with pytest.raises(ValueError, match=message):
        _definition(substances=[_api(**kwargs)], chemistry=[], routes=[]).validate()


@pytest.mark.parametrize(
    "kinetics",
    [FirstOrder(0.0), MM(vmax=-1.0, km=1.0), MM(vmax=1.0, km=0.0), Hill(1.0, 1.0, 0.0)],
)
def test_kinetics_numbers_checked(kinetics: FirstOrder | MM | Hill) -> None:
    """Kinetic parameters must be positive."""
    m = Metabolism("met1", Organ.LIVER, "api", {"api_m2": 1.0}, kinetics)
    with pytest.raises(ValueError, match=r"met1: .* must be positive"):
        _definition(chemistry=[m]).validate()


def test_product_fraction_must_be_positive() -> None:
    """A zero or negative product fraction is rejected."""
    m = Metabolism(
        "met1", Organ.LIVER, "api", {"api_m2": 1.0, "api_x": 0.0}, FirstOrder(1.0)
    )
    with pytest.raises(ValueError, match="product fraction of 'api_x'"):
        _definition(
            substances=[
                _api(),
                Substance("api_m2", "m", 1.0),
                Substance("api_x", "x", 1.0),
            ],
            chemistry=[m],
        ).validate()


def test_observer_with_unknown_symbol_is_rejected() -> None:
    """A symbol of an observer which is no id of the top model is named."""
    with pytest.raises(ValueError, match=r"'x' reads \['Cnope'\]"):
        _definition(observers=[Observer("x", "Cnope + 1 mM", U.mM, "x")]).validate()


def test_observer_with_an_organ_id_is_rejected() -> None:
    """An id of an organ (`LI__Cpl_api`) is not visible in the top model."""
    observer = Observer("x", "LI__Cpl_api", U.mM, "x")
    with pytest.raises(ValueError, match="LI__Cpl_api"):
        _definition(observers=[observer]).validate()


def test_observer_with_a_body_id_is_rejected() -> None:
    """An observer id which is an id of the body is named."""
    observer = Observer("Cve_api", "Car_api", U.mM, "x")
    with pytest.raises(ValueError, match="observer id 'Cve_api' is an id of the body"):
        _definition(observers=[observer]).validate()


def test_observer_ids_and_formulas_checked() -> None:
    """Observer ids are SIds without `__`, unique; a formula must parse."""
    observers = [
        Observer("a", "Carm_api", U.mM, "a"),
        Observer("a", "Carm_api", U.mM, "a again"),
        Observer("1b", "Carm_api", U.mM, "b"),
        Observer("c__d", "Carm_api", U.mM, "c"),
        Observer("e", "Carm_api +* 2", U.mM, "e"),
    ]
    with pytest.raises(ValueError) as info:
        _definition(observers=observers).validate()
    message = str(info.value)
    assert "duplicate observer ids: ['a']" in message
    assert "observer id '1b' is not an SBML SId" in message
    assert "observer id 'c__d' contains '__'" in message
    assert "observer 'e': formula 'Carm_api +* 2' cannot be parsed" in message


def test_observer_symbols_checked_without_compose() -> None:
    """`validate(symbols=False)` leaves the symbols to `compose`."""
    observer = Observer("x", "Cnope", U.mM, "x")
    _definition(observers=[observer]).validate(symbols=False)


def test_observer_formula_symbols() -> None:
    """The symbols of a formula by libsbml's parser: no units, no functions, `time`."""
    from pkpdsim.pbpk.observers import formula_symbols

    assert formula_symbols("exp(-k * time) * 2 mM + piecewise(a, b > 1 mM, c)") == {
        "k",
        "time",
        "a",
        "b",
        "c",
    }
