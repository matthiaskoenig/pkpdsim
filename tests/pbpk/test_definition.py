"""Tests of the definition dataclasses and their validation."""

import dataclasses
import re
from collections.abc import Callable
from pathlib import Path

import libsbml
import pytest
import roadrunner
from sbmlutils.factory import Parameter, Uncertainty, UncertParameter

from pkpdsim.pbpk import (
    MM,
    PROTEINS,
    Binding,
    CustomReaction,
    Direction,
    Erythrocytes,
    Excretion,
    FirstOrder,
    FirstOrderDissolution,
    Formula,
    Hill,
    Inhibition,
    Kp,
    KpMethod,
    Membrane,
    Metabolism,
    Mode,
    ModelParameter,
    Observer,
    Organ,
    Part,
    PBPKDefinition,
    Route,
    Sex,
    Substance,
    Transport,
    U,
    Value,
    Variants,
    ZonatedLiver,
    build,
    compose,
    ids,
    load,
    sources,
)
from pkpdsim.pbpk.builder import organ_component

from .conftest import (
    assert_reference_sbml,
    dateless,
    definition,
    individual,
    reference_definition,
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


def _liver_only(mode: Mode, organs: frozenset[Organ]) -> PBPKDefinition:
    """A metabolite which enters the given organs, the liver with an exchange of its own."""
    m = Substance(
        "m",
        "metabolite",
        300.0,
        distribution=Kp(1.0, mode=mode, organs=organs, ps={Organ.LIVER: 1.0}),
        excretion=Excretion(urine=True),
    )
    bile = Transport(
        "MBIL",
        Organ.LIVER,
        "m",
        Direction.EFFLUX,
        Membrane.APICAL,
        FirstOrder(0.01),
    )
    return PBPKDefinition(
        sid="liver_only", substances=[m], chemistry=[bile], routes=[Route.iv("m")]
    )


@pytest.mark.parametrize("mode", [Mode.TISSUE, Mode.PERMEABILITY])
def test_no_ftissue_when_every_entered_organ_has_its_own_ps(
    mode: Mode, tmp_path: Path
) -> None:
    """`ftissue` is needed only for an entered organ without an exchange of its own.

    A metabolite which enters only the liver, with the liver exchange `ps` and
    a biliary excretion, is built without `ftissue`: the liver has
    `PS = kps*Vli*f_functional`, no organ has a share of a whole-body exchange.
    """
    d = _liver_only(mode, frozenset({Organ.LIVER}))
    d.validate()
    flat = build(d, tmp_path).flat_path
    model = libsbml.readSBMLFromFile(str(flat)).getModel()
    assert model.getParameter(ids.ftissue("m")) is None
    assert model.getParameter(ids.tissue_volume_total_ref("m")) is None
    assert model.getParameter("LI__kps_m").getValue() == 1.0
    ia = model.getInitialAssignmentBySymbol("LI__PS_m")
    assert libsbml.formulaToL3String(ia.getMath()) == "LI__kps_m * Vli * f_functional"
    r = roadrunner.RoadRunner(str(flat))
    r["IVDOSE_m"] = 1.0
    r.simulate(0, 600, 11)
    assert r["LI__PS_m"] == pytest.approx(r["Vli"] * r["f_functional"])
    assert r["Cbile_m"] > 0


@pytest.mark.parametrize("mode", [Mode.TISSUE, Mode.PERMEABILITY])
def test_ftissue_for_an_entered_organ_without_ps(mode: Mode) -> None:
    """An entered organ without an exchange of its own needs `ftissue`."""
    d = _liver_only(mode, frozenset({Organ.LIVER, Organ.MUSCLE}))
    with pytest.raises(ValueError, match=r"requires ftissue.*\['mu'\]"):
        d.validate()


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
        "ABCC3",
        "ABCC4",
        "SLC6A8",
    ]
    assert set(required) <= set(PROTEINS)
    assert PROTEINS["CYP3A4"].uniprot == "P08684"
    assert PROTEINS["SLC22A6"].ncit is None
    assert PROTEINS["ABCC3"].uniprot == "O15438"
    assert PROTEINS["ABCC4"].uniprot == "O15439"
    assert PROTEINS["SLC6A8"].uniprot == "P48029"
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
            {
                "distribution": Kp(
                    1.0, ftissue=0.5, ps={Organ.LIVER: Value(0.2, U.l_per_min)}
                )
            },
            "ps 'li' must be in 1/min, got l/min",
        ),
        (
            {
                "distribution": Kp(
                    1.0,
                    ftissue=0.5,
                    ps={Organ.LIVER: Value(0.2, U.per_min)},
                    ps_per_volume=False,
                )
            },
            "ps 'li' must be in l/min, got 1/min",
        ),
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


def test_observer_with_an_organ_id_which_is_not_exposed_is_rejected() -> None:
    """An organ exposes its species (`LI__Cpl_api`) to an observer, no other element (`LI__Vpl`)."""
    _definition(observers=[Observer("x", "LI__Cpl_api", U.mM, "x")]).validate()
    observer = Observer("x", "LI__Vpl", U.mM, "x")
    with pytest.raises(ValueError, match="'LI__Vpl', which is no species"):
        _definition(observers=[observer]).validate()


def test_observer_with_a_body_id_is_rejected() -> None:
    """An observer id which is an id of the body is named."""
    observer = Observer("Cve_api", "Car_api", U.mM, "x")
    with pytest.raises(ValueError, match="observer id 'Cve_api' is an id of the body"):
        _definition(observers=[observer]).validate()


@pytest.mark.parametrize("sid", ["LI", "standard_li"])
def test_observer_reading_a_submodel_id_is_rejected(sid: str) -> None:
    """A submodel or model definition id is no value: a formula cannot read it.

    It stays an id of the body, which an observer must not take.
    """
    observer = Observer("x", f"{sid} + Car_api", U.mM, "x")
    with pytest.raises(ValueError, match=rf"'x' reads \['{sid}'\]"):
        _definition(observers=[observer]).validate()
    observer = Observer(sid, "Car_api", U.mM, "x")
    with pytest.raises(ValueError, match=f"observer id '{sid}' is an id of the body"):
        _definition(observers=[observer]).validate()


def test_observer_calling_an_unknown_function_is_rejected() -> None:
    """A call of a function which SBML L3 infix does not define is named.

    The models of the library have no function definitions; the built-in
    functions pass.
    """
    from pkpdsim.pbpk.observers import formula_functions

    assert formula_functions("exp(a) * piecewise(b, c > 1, 0) + foo(d, bar(e))") == {
        "foo",
        "bar",
    }
    observer = Observer("x", "foo(Car_api)", U.mM, "x")
    with pytest.raises(ValueError, match=r"observer 'x' calls \['foo'\]"):
        _definition(observers=[observer]).validate(symbols=False)
    _definition(
        observers=[Observer("y", "exp(-Car_api / 1 mM) * 1 mM", U.mM, "y")]
    ).validate()


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


def _untyped[T](cls: Callable[..., T], *args: object, **kwargs: object) -> T:
    """Construct a dataclass from untyped arguments, as a script or a notebook does."""
    return cls(*args, **kwargs)


@pytest.mark.parametrize(
    "organ", [Organ.ADIPOSE, "ad", "adipose", "Adipose", "ADIPOSE"]
)
def test_kp_overrides_accept_organ_texts(organ: object) -> None:
    """The keys of `Kp.overrides` are organs, their values or their names (R60)."""
    kp = _untyped(Kp, 2.0, ftissue=0.5, overrides={organ: 0.4})
    assert kp.overrides == {Organ.ADIPOSE: 0.4}
    assert all(type(k) is Organ for k in kp.overrides)


def test_organ_fields_are_normalised() -> None:
    """`Kp.organs`, `Kp.ps`, `KpMethod` and `Variants.organs` hold `Organ` members."""
    for organs in ({"liver", "ki"}, frozenset({"liver", "ki"}), ["li", "kidney"]):
        kp = _untyped(Kp, 1.0, organs=organs, ftissue=0.5)
        assert kp.organs == frozenset({Organ.LIVER, Organ.KIDNEY})
        assert isinstance(kp.organs, frozenset)
    assert Kp(1.0).organs is None
    kp = _untyped(Kp, 1.0, ps={"kidney": 0.2})
    assert kp.ps == {Organ.KIDNEY: 0.2}
    method = _untyped(
        KpMethod,
        method="rodgers_rowland",
        logp=1.0,
        pka=(),
        compound="neutral",
        organs=("liver",),
        ps={"li": 0.1},
    )
    assert method.organs == frozenset({Organ.LIVER})
    assert method.ps == {Organ.LIVER: 0.1}
    variants = _untyped(Variants, organs={"liver": "zonated"})
    assert variants.organs == {Organ.LIVER: "zonated"}
    replaced = _untyped(dataclasses.replace, kp, overrides={"adipose": 0.4})
    assert replaced.overrides == {Organ.ADIPOSE: 0.4}


_METHOD = {"method": "rodgers_rowland", "logp": 1.0, "pka": (), "compound": "neutral"}
"""The arguments of a Kp method but the organs."""


@pytest.mark.parametrize(
    ("cls", "kwargs", "field"),
    [
        (Kp, {"value": 2.0, "overrides": {"adipos": 0.4}}, "Kp.overrides"),
        (Kp, {"value": 2.0, "overrides": [("ad", 0.4)]}, "Kp.overrides"),
        (Kp, {"value": 2.0, "ps": {3: 0.2}}, "Kp.ps"),
        (Kp, {"value": 2.0, "organs": {"lever"}}, "Kp.organs"),
        (Kp, {"value": 2.0, "organs": "liver"}, "Kp.organs"),
        (Kp, {"value": 2.0, "organs": {"li": 1.0}}, "Kp.organs"),
        (KpMethod, {**_METHOD, "organs": {"brian"}}, "KpMethod.organs"),
        (KpMethod, {**_METHOD, "ps": {"brian": 0.1}}, "KpMethod.ps"),
        (Variants, {"organs": {"lever": "zonated"}}, "Variants.organs"),
        (Variants, {"organs": {"li": "zonated", "liver": "x"}}, "Variants.organs"),
    ],
)
def test_organ_fields_reject_other_values(
    cls: Callable[..., object], kwargs: dict[str, object], field: str
) -> None:
    """Anything but an organ, its value or its name raises, naming the field (R60).

    The message lists the allowed values, or says what is wrong with the
    collection (no mapping, no set, an organ given twice).
    """
    with pytest.raises(ValueError, match=re.escape(field)) as err:
        _untyped(cls, **kwargs)
    message = str(err.value)
    assert any(s in message for s in ("'li' (liver)", "given twice", "e.g."))


@pytest.mark.parametrize("sex", [Sex.FEMALE, "female", "FEMALE", "Female"])
def test_sex_accepts_its_value(sex: object) -> None:
    """The sex is a `Sex` member or its value (R60), also through `replace`."""
    assert _definition(sex=sex).sex is Sex.FEMALE
    assert _untyped(dataclasses.replace, _definition(), sex=sex).sex is Sex.FEMALE


def test_sex_rejects_other_values() -> None:
    """Another sex raises with the field and the allowed values."""
    with pytest.raises(
        ValueError, match=r"PBPKDefinition\.sex.*'male' \(male\), 'female'"
    ):
        _definition(sex="f")


def test_definition_with_texts_composes() -> None:
    """The examples of the spec compose: the texts select the organs and the sex.

    `overrides={"adipose": 0.4}`, `sex="female"` and a variant by the name of
    the organ give the override of adipose, the female reference and the
    zonated liver.
    """
    kp = _untyped(Kp, 2.0, ftissue=0.5, overrides={"adipose": 0.4})
    d = _definition(
        substances=[_api(distribution=kp), Substance("api_m2", "metabolite", 445.5)],
        sex="female",
        variants=_untyped(Variants, organs={"liver": "zonated"}),
    )
    d.validate()
    assert isinstance(organ_component(d, Organ.LIVER), ZonatedLiver)
    parameters = {str(p.sid): p.value for p in compose(d).parameters or []}
    assert parameters[ids.kp_organ(Organ.ADIPOSE, "api")] == 0.4
    assert parameters[ids.BW] == 60.0


def _text_processes() -> tuple[Metabolism, Transport, CustomReaction]:
    """The processes of the chemistry with every enum field given as text."""
    return (
        _untyped(Metabolism, "EL", "liver", "api", {"api_m2": 1.0}, MM(1e-3, 0.01)),
        _untyped(Transport, "UP", "li", "api", "influx", "BASOLATERAL", MM(1e-3, 0.01)),
        _untyped(
            CustomReaction,
            sid="KX",
            organ="kidney",
            part="tissue",
            equation="Cti_api ->",
            formula="kx * Cti_api",
            parameters=[Parameter("kx", 0.01, U.l_per_min)],
        ),
    )


def test_chemistry_enum_fields_accept_texts(tmp_path: Path) -> None:
    """Organ, direction, membrane, part and mode given as text build the reactions (R61).

    Before, a process of a text organ was dropped from the model.
    """
    metabolism, transport, custom = _text_processes()
    assert metabolism.organ is Organ.LIVER
    assert (transport.organ, transport.direction, transport.membrane) == (
        Organ.LIVER,
        Direction.INFLUX,
        Membrane.BASOLATERAL,
    )
    assert type(transport.direction) is Direction
    assert (custom.organ, custom.part) == (Organ.KIDNEY, Part.TISSUE)
    kp = _untyped(Kp, 2.0, mode="TISSUE", ftissue=0.5)
    assert kp.mode is Mode.TISSUE
    method = _untyped(KpMethod, **_METHOD, mode="perfusion")
    assert method.mode is Mode.PERFUSION
    d = _definition(
        substances=[_api(distribution=kp), Substance("api_m2", "metabolite", 445.5)],
        chemistry=list(_text_processes()),
    )
    doc = libsbml.readSBMLFromFile(str(build(d, tmp_path).flat_path))
    reactions = {r.getId() for r in doc.getModel().getListOfReactions()}
    assert {"LI__EL", "LI__UP", "KI__KX"} <= reactions


@pytest.mark.parametrize(
    ("cls", "args", "field"),
    [
        (Metabolism, ("EL", "lever", "api", {}, MM(1e-3, 0.01)), "Metabolism.organ"),
        (
            Transport,
            ("T", "lever", "api", "influx", "apical", MM(1e-3, 0.01)),
            "Transport.organ",
        ),
        (
            Transport,
            ("T", "li", "api", "inflow", "apical", MM(1e-3, 0.01)),
            "Transport.direction",
        ),
        (
            Transport,
            ("T", "li", "api", "influx", "apikal", MM(1e-3, 0.01)),
            "Transport.membrane",
        ),
        (CustomReaction, ("C", "brian", "ti", "", "", []), "CustomReaction.organ"),
        (CustomReaction, ("C", "br", "cell", "", "", []), "CustomReaction.part"),
        (Kp, (2.0,), "Kp.mode"),
    ],
)
def test_chemistry_enum_fields_reject_other_texts(
    cls: Callable[..., object], args: tuple[object, ...], field: str
) -> None:
    """An unknown text raises with the field and the allowed values (R61)."""
    kwargs = {"mode": "flow"} if cls is Kp else {}
    with pytest.raises(ValueError, match=re.escape(field)) as err:
        _untyped(cls, *args, **kwargs)
    assert "allowed values are" in str(err.value)


def test_basolateral_text_hits_the_gut_refusal() -> None:
    """A basolateral transport of the intestine given as text is refused as with the member."""
    transport = _untyped(
        Transport,
        "GUTUP",
        "small_intestine",
        "api",
        "influx",
        "basolateral",
        FirstOrder(0.1),
    )
    d = _definition(routes=[Route.iv("api"), Route.oral("api")], chemistry=[transport])
    with pytest.raises(ValueError, match="Transport 'GUTUP' is basolateral"):
        compose(d)


@pytest.mark.parametrize(
    ("make", "field"),
    [
        (lambda: _untyped(Binding, 0.1, "albumine"), "Binding.protein"),
        (lambda: _untyped(Erythrocytes, "kinetik", kbc=2.0), "Erythrocytes.mode"),
        (lambda: _untyped(Route, "ivv", "api"), "Route.kind"),
        (lambda: _untyped(Inhibition, "CYP3A4", "competitve", 1.0), "Inhibition.kind"),
        (
            lambda: _untyped(KpMethod, **{**_METHOD, "compound": "bas"}),
            "KpMethod.compound",
        ),
        (lambda: _untyped(KpMethod, **{**_METHOD, "method": "rr"}), "KpMethod.method"),
    ],
)
def test_choices_reject_other_texts(make: Callable[[], object], field: str) -> None:
    """A text which is none of the choices of a `Literal` field raises (R61).

    Before, a route of an unknown kind composed no doses and an unknown
    erythrocyte mode was taken for another.
    """
    with pytest.raises(ValueError, match=re.escape(field) + ".*allowed values"):
        make()


# --- parameters of the definition ------------------------------------------

#: a fitted production rate with its uncertainty, as `fitted(FITS.parameter(...))`
K_PROD = Value(
    0.002,
    U.mmole_per_min,
    "fit",
    notes="The fitted production.",
    uncertainties=[
        Uncertainty(
            uncertParameters=[
                UncertParameter(
                    type=libsbml.DISTRIB_UNCERTTYPE_DISTRIBUTION,
                    math="normal(0.002, 0.0002)",
                )
            ]
        )
    ],
)


def _k_prod(value: float | Value = K_PROD) -> ModelParameter:
    """The production rate `k_prod` of a definition."""
    return ModelParameter(
        "k_prod",
        value,
        U.mmole_per_min,
        name="production rate",
        notes="The production shared by the baseline and the muscle.",
    )


def test_a_definition_parameter_is_a_top_parameter(tmp_path: Path) -> None:
    """A parameter of the definition is a constant parameter of the top model.

    It keeps value, unit, name, notes and the uncertainty of a fitted value
    (the distribution of `check_model`); a `Formula` parameter is the initial
    assignment over top ids (another parameter of the definition, `BW`) and
    follows a change before the simulation; an observer of the definition
    reads the parameters.
    """
    d = _definition(
        parameters=[
            _k_prod(),
            ModelParameter(
                "k_scaled",
                Formula("2 dimensionless * k_prod * BW / 70 kg", U.mmole_per_min),
                U.mmole_per_min,
            ),
        ],
        observers=[Observer("A_prod", "k_prod * time", U.mmole, "produced amount")],
    )
    d.validate()
    result = build(d, tmp_path)
    model = libsbml.readSBMLFromFile(str(result.flat_path)).getModel()
    p = model.getParameter("k_prod")
    assert p.getConstant()
    assert p.getValue() == 0.002
    assert p.getUnits() == "mmole_per_min"
    assert p.getName() == "production rate"
    notes = p.getNotesString()
    assert "The production shared by the baseline and the muscle." in notes
    assert "The fitted production." in notes
    assert "Fitted value." in notes
    distrib = p.getPlugin("distrib")
    assert distrib is not None and distrib.getNumUncertainties() == 1
    assert p.getSBOTermID() == "SBO:0000002"  # quantitative systems description
    assert p.isSetMetaId() and p.getCVTerms().getSize() > 0
    scaled = model.getParameter("k_scaled")
    assert scaled.getName() == "k_scaled"
    assert scaled.getConstant()
    ia = model.getInitialAssignmentBySymbol("k_scaled")
    assert "k_prod" in libsbml.formulaToL3String(ia.getMath())
    r = individual(result.flat_path, k_prod=0.004, BW=70.0)
    assert r["k_scaled"] == pytest.approx(0.008)
    r.simulate(0, 10)
    assert r["A_prod"] == pytest.approx(0.04)
    # the parameter is no input of an organ which does not read it
    comp = libsbml.readSBMLFromFile(str(result.comp_path))
    plugin = comp.getPlugin("comp")
    for k in range(plugin.getNumModelDefinitions()):
        md = plugin.getModelDefinition(k)
        assert md.getParameter("k_prod") is None, md.getId()


@pytest.mark.parametrize(
    ("sid", "where"),
    [
        ("GFR", "physiology"),
        ("BW", "physiology"),
        ("fu_api", "top model"),
        ("Kp_api", "top model"),
        ("Carm_api", "top model"),
        ("Aurine_api", "top model"),
        ("LI", "organ submodel"),
        ("GU", "library"),
        ("Vti", "organ"),
        ("met_Vmax", "organ"),
    ],
)
def test_a_definition_parameter_may_not_shadow_a_library_id(
    sid: str, where: str
) -> None:
    """An id of the library is refused, naming it.

    The ids of the physiology are refused by `validate`, every other id of
    the top model (the parameters of the substances, the observers, the
    sinks, the submodels) and the local ids of an organ (which a custom
    reaction of the organ would read instead) when the model is composed.
    """
    d = _definition(parameters=[ModelParameter(sid, 1.0, U.dimensionless)])
    with pytest.raises(ValueError, match=rf"parameter '{sid}'.*{where}"):
        compose(d)


def test_the_parameters_of_a_definition_are_validated() -> None:
    """SIds without `__`, unique, a `Value` or `Formula` in the declared unit."""
    d = _definition(
        parameters=[
            ModelParameter("a", 1.0, U.dimensionless),
            ModelParameter("a", 2.0, U.dimensionless),
            ModelParameter("1b", 1.0, U.dimensionless),
            ModelParameter("c__d", 1.0, U.dimensionless),
            ModelParameter("e", Value(1.0, U.mM), U.per_min),
            ModelParameter("f", Formula("2 * a", U.mM), U.per_min),
        ]
    )
    with pytest.raises(ValueError) as info:
        d.validate()
    message = str(info.value)
    assert "duplicate parameter ids: ['a']" in message
    assert "parameter id '1b' is not an SBML SId" in message
    assert "parameter id 'c__d' contains '__'" in message
    assert "parameter 'e': value must be in 1/min, got mmole/liter" in message
    assert "parameter 'f': value must be in 1/min, got mmole/liter" in message


def test_the_formulas_of_parameters_may_not_read_each_other_in_a_cycle() -> None:
    """Two formulas which read each other are refused, naming them."""
    d = _definition(
        parameters=[
            ModelParameter("a", Formula("2 * b", U.dimensionless), U.dimensionless),
            ModelParameter("b", Formula("a / 2", U.dimensionless), U.dimensionless),
        ]
    )
    with pytest.raises(ValueError, match=r"\['a', 'b'\] read each other in a cycle"):
        compose(d)


def test_a_custom_reaction_reads_no_unknown_id() -> None:
    """An id of a custom reaction which neither the organ nor the top model has is named."""
    custom = CustomReaction(
        "PROD",
        Organ.MUSCLE,
        Part.TISSUE,
        "-> Cti_api",
        "k_unknown * k_prod",
        parameters=[],
    )
    d = _definition(chemistry=[custom], parameters=[_k_prod()])
    with pytest.raises(ValueError, match=r"'PROD'.*\['k_unknown'\]"):
        compose(d)
    observer = CustomReaction(
        "PROD", Organ.MUSCLE, Part.TISSUE, "-> Cti_api", "Carm_api", parameters=[]
    )
    with pytest.raises(ValueError, match=r"'PROD'.*\['Carm_api'\]"):
        compose(_definition(chemistry=[observer]))


def test_combine_unites_the_parameters() -> None:
    """Equal parameters of two definitions are taken once, different ones raise."""
    a = _definition(parameters=[_k_prod()])
    b = _definition(parameters=[_k_prod(), ModelParameter("x", 1.0, U.dimensionless)])
    d = PBPKDefinition.combine(a, b, sid="coupled")
    assert [p.sid for p in d.parameters] == ["k_prod", "x"]
    with pytest.raises(ValueError, match="parameter 'k_prod'"):
        PBPKDefinition.combine(a, _definition(parameters=[_k_prod(1.0)]), sid="c")


def test_a_definition_without_parameters_is_unchanged(tmp_path: Path) -> None:
    """Without parameters the SBML is that of f008d56, line by line.

    With a parameter which nothing reads the models differ only by its
    element in the top model.
    """
    plain = reference_definition()
    assert plain.parameters == ()
    plain_result = build(plain, tmp_path / "plain")
    assert_reference_sbml(plain_result)
    with_parameter = dataclasses.replace(
        plain, parameters=[ModelParameter("k_prod", 0.002, U.mmole_per_min)]
    )
    result = build(with_parameter, tmp_path / "parameter")
    for path, reference in (
        (result.comp_path, plain_result.comp_path),
        (result.flat_path, plain_result.flat_path),
    ):
        doc = libsbml.readSBMLFromFile(str(path))
        assert doc.getModel().removeParameter("k_prod") is not None
        assert dateless(libsbml.writeSBMLToString(doc)) == dateless(
            reference.read_text(encoding="utf-8")
        )


def test_a_definition_parameter_may_not_take_an_id_of_the_definition() -> None:
    """An observer of the definition or a PD component under the id is named as such."""
    from .conftest import cardiac_pd

    observer = Observer("A_prod", "1 mmole_per_min * time", U.mmole, "produced")
    d = _definition(
        parameters=[ModelParameter("A_prod", 1.0, U.mmole)], observers=[observer]
    )
    with pytest.raises(
        ValueError,
        match="parameter 'A_prod' is the id of an observer of the definition",
    ):
        compose(d)
    d = _definition(
        parameters=[ModelParameter("cardio", 1.0, U.dimensionless)],
        pd=[cardiac_pd(reads="api")],
    )
    with pytest.raises(ValueError, match="parameter 'cardio' is an id of the PD"):
        compose(d)


def test_a_constant_of_a_formulation_is_read_from_the_definition(
    tmp_path: Path,
) -> None:
    """A constant of a formulation lives in the definition and the model (D14).

    The experiments of a prolonged release set the factor of the dissolution
    from the constant of the definition, which the archive carries.
    """
    released = definition(
        routes=[Route.oral("drg", FirstOrderDissolution(k=2.0, release=1.0))],
        parameters=[
            ModelParameter(
                "f_dis_pr_drg",
                Value(0.27, U.dimensionless, notes="calibrated"),
                U.dimensionless,
                notes="The factor of the dissolution of the prolonged release.",
            ),
            ModelParameter(
                "f_twice",
                Formula("2 dimensionless * f_dis_pr_drg", U.dimensionless),
                U.dimensionless,
            ),
        ],
    )
    assert released.parameter_value("f_dis_pr_drg") == 0.27
    result = build(released, tmp_path)
    r = roadrunner.RoadRunner(str(result.flat_path))
    assert r["f_dis_pr_drg"] == 0.27
    with pytest.raises(ValueError, match=r"no parameter 'f_dis_ir_drg'.*f_dis_pr_drg"):
        released.parameter_value("f_dis_ir_drg")
    with pytest.raises(ValueError, match="formula"):
        released.parameter_value("f_twice")
    plain = definition(parameters=[ModelParameter("f", 0.5, U.dimensionless)])
    assert plain.parameter_value("f") == 0.5
