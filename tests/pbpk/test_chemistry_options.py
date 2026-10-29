"""The options of the chemistry: custom parameters, shared rates, product fractions.

- the parameters of a custom reaction as `ModelParameter`s, local ids of the
  organ with value, provenance and the uncertainty of a fit (D3);
- the capacity of a process as a ratio of another process of the organ,
  `RateOf` (D4);
- the product fractions of a metabolism as values and fitted values with the
  remainder `REST`, the parameters `f_<metabolism>_<product>` (D5);
- a definition which uses none of them builds the SBML of 0c538ca.
"""

import dataclasses
from collections.abc import Mapping
from pathlib import Path

import libsbml
import pytest
from sbmlutils.factory import Parameter, Uncertainty, UncertParameter

from pkpdsim.pbpk import (
    MM,
    REST,
    CustomReaction,
    Direction,
    FirstOrder,
    Formula,
    Hill,
    Membrane,
    Metabolism,
    ModelParameter,
    Organ,
    Part,
    PBPKDefinition,
    RateOf,
    Route,
    Substance,
    Transport,
    U,
    Value,
    Variants,
    build,
    ids,
)
from pkpdsim.pbpk.builder import BuildResult

from .conftest import (
    assert_reference_sbml,
    chemistry_reference_definition,
    dose_amount,
    individual,
    load,
    simple_drug,
    total_amount,
)

DAY = 24 * 60.0


def _uncertainty(distribution: str) -> list[Uncertainty]:
    """The uncertainty of a fitted value with the distribution of its draws."""
    return [
        Uncertainty(
            uncertParameters=[
                UncertParameter(
                    type=libsbml.DISTRIB_UNCERTTYPE_DISTRIBUTION,
                    math=distribution,
                    name="distribution of the draws",
                )
            ],
            name="fit F1",
            notes="The value is the result of the fit `F1`.",
        )
    ]


def _fitted(value: float, distribution: str, unit: object = U.dimensionless) -> Value:
    """A fitted value with the uncertainty of the fit."""
    return Value(
        value,
        unit,  # ty: ignore[invalid-argument-type]
        source="fit",
        notes="Fitted in a test.",
        uncertainties=_uncertainty(distribution),
    )


def _substances() -> list[Substance]:
    """The drug and three metabolites, which stay in the body."""
    return [
        simple_drug(),
        Substance("drg_m1", name="metabolite 1", mr=316.0),
        Substance("drg_m2", name="metabolite 2", mr=332.0),
        Substance("drg_m3", name="metabolite 3", mr=348.0),
    ]


def _definition(
    chemistry: list[Metabolism | Transport | CustomReaction],
    sid: str = "options",
    **fields: object,
) -> PBPKDefinition:
    """A definition of the drug and its metabolites with an IV route."""
    return PBPKDefinition(
        sid=sid,
        substances=_substances(),
        chemistry=chemistry,
        routes=[Route.iv("drg")],
        **fields,  # ty: ignore[invalid-argument-type]
    )


_DOCUMENTS: list[libsbml.SBMLDocument] = []
"""The documents read by the tests: a libsbml model belongs to its document,
which must outlive it."""


def _read(path: Path) -> libsbml.SBMLDocument:
    """Read an SBML document and keep it alive."""
    doc = libsbml.readSBMLFromFile(str(path))
    _DOCUMENTS.append(doc)
    return doc


def _flat(result: BuildResult) -> libsbml.Model:
    """The flat model of a build."""
    model = _read(result.flat_path).getModel()
    assert model is not None
    return model


def _comp_definition(result: BuildResult, organ: Organ) -> libsbml.Model:
    """The model definition of an organ in the comp model."""
    comp = _read(result.comp_path).getPlugin("comp")
    for k in range(comp.getNumModelDefinitions()):
        md = comp.getModelDefinition(k)
        if md.getId().endswith(f"_{organ.value}"):
            return md
    raise AssertionError(f"no model definition of {organ.value}")


def _assignment(model: libsbml.Model, sid: str) -> str:
    """The formula of the initial assignment of an id."""
    ia = model.getInitialAssignmentBySymbol(sid)
    assert ia is not None, sid
    return libsbml.formulaToL3String(ia.getMath())


def _notes(model: libsbml.Model, sid: str) -> str:
    """The notes of a parameter as text."""
    p = model.getParameter(sid)
    assert p is not None, sid
    return p.getNotesString()


def _distribution(model: libsbml.Model, sid: str) -> str | None:
    """The distribution of the draws of a parameter, `None` without."""
    p = model.getParameter(sid)
    assert p is not None, sid
    plugin = p.getPlugin("distrib")
    for k in range(plugin.getNumUncertainties() if plugin else 0):
        u = plugin.getUncertainty(k)
        for j in range(u.getNumUncertParameters()):
            child = u.getUncertParameter(j)
            if child.getType() == libsbml.DISTRIB_UNCERTTYPE_DISTRIBUTION:
                return libsbml.formulaToL3String(child.getMath())
    return None


def _problems(d: PBPKDefinition) -> str:
    """The message of the validation of a definition."""
    with pytest.raises(ValueError) as err:
        d.validate()
    return str(err.value)


# --- custom parameters (D3) ---------------------------------------------------


def _deep(*parameters: ModelParameter | Parameter) -> CustomReaction:
    """A deep pool of the muscle tissue with the exchange constant `k_deep`."""
    from sbmlutils.factory import Species

    return CustomReaction(
        "DRGDEEP",
        Organ.MUSCLE,
        Part.TISSUE,
        "Cti_drg <-> Cdeep_drg",
        "k_deep * Vti * (Kp_deep * Cti_drg - Cdeep_drg)",
        parameters=list(parameters),
        species=[
            Species(
                "Cdeep_drg",
                "Vti",
                initialConcentration=0.0,
                substanceUnit=U.mmole,
                hasOnlySubstanceUnits=False,
            )
        ],
    )


_KP_DEEP = ModelParameter(
    "Kp_deep",
    Value(1.5, U.dimensionless, source="ModelChoice", notes="The deep pool to tissue."),
    U.dimensionless,
    name="deep pool to tissue",
)


def test_a_custom_parameter_is_a_local_id_of_the_organ(tmp_path: Path) -> None:
    """A `ModelParameter` of a custom reaction is a constant of its organ.

    Local id in the model definition of the organ, `<ORGAN>__<id>` in the flat
    model, with value, unit, name, notes and the citation of its source; no
    parameter of the top model under the id.
    """
    k_deep = ModelParameter(
        "k_deep",
        2e-3,
        U.per_min,
        name="exchange of the deep pool",
        notes="Per litre of tissue.",
    )
    result = build(_definition([_deep(k_deep, _KP_DEEP)]), tmp_path)
    md = _comp_definition(result, Organ.MUSCLE)
    assert md.getParameter("k_deep") is not None
    flat = _flat(result)
    p = flat.getParameter(ids.flat(Organ.MUSCLE, "k_deep"))
    assert p is not None
    assert p.getValue() == 2e-3
    assert p.getConstant()
    assert md.getParameter("k_deep").getUnits() == "per_min"
    assert p.getName() == "exchange of the deep pool"
    assert "Per litre of tissue." in p.getNotesString()
    assert flat.getParameter("k_deep") is None
    kp = _notes(flat, ids.flat(Organ.MUSCLE, "Kp_deep"))
    assert "The deep pool to tissue." in kp
    assert "Choice of the pkpdsim PBPK library" in kp


def test_a_fitted_custom_parameter_carries_its_distribution_and_provenance(
    tmp_path: Path,
) -> None:
    """A fitted custom parameter has the uncertainty of the fit and its notes.

    The distribution of the draws (`run_simulation --uncertainty`, the check of
    `check_model`) is that of the fit record, the notes say it is fitted.
    """
    k_deep = ModelParameter(
        "k_deep",
        _fitted(2e-3, "lognormal(-6.2, 0.3, 1e-5, 0.1)", U.per_min),
        U.per_min,
        name="exchange of the deep pool",
    )
    flat = _flat(build(_definition([_deep(k_deep, _KP_DEEP)]), tmp_path))
    sid = ids.flat(Organ.MUSCLE, "k_deep")
    assert flat.getParameter(sid).getValue() == 2e-3
    assert _distribution(flat, sid) == "lognormal(-6.2, 0.3, 1e-5, 0.1)"
    notes = _notes(flat, sid)
    assert "Fitted in a test." in notes
    assert "Fitted value." in notes


def test_the_distribution_of_a_custom_parameter_is_drawn(tmp_path: Path) -> None:
    """The draws of the uncertainty runs find the distribution of a custom parameter."""
    from pkpdsim.uncertainty import parameter_distributions

    k_deep = ModelParameter(
        "k_deep",
        _fitted(2e-3, "normal(0.002, 0.0002, 0.001, 0.003)", U.per_min),
        U.per_min,
    )
    result = build(_definition([_deep(k_deep, _KP_DEEP)]), tmp_path)
    assert ids.flat(Organ.MUSCLE, "k_deep") in parameter_distributions(result.flat_path)


def test_a_custom_parameter_of_a_formula_reads_the_top_model(tmp_path: Path) -> None:
    """A `Formula` of a custom parameter is its initial assignment over top ids.

    The organ reads the ids of the formula through input ports, the value
    follows a change of a covariate before the simulation.
    """
    k_deep = ModelParameter(
        "k_deep", Formula("2e-3 per_min * BW / 70 kg", U.per_min), U.per_min
    )
    result = build(_definition([_deep(k_deep, _KP_DEEP)]), tmp_path)
    md = _comp_definition(result, Organ.MUSCLE)
    assert "BW" in _assignment(md, "k_deep")
    sid = ids.flat(Organ.MUSCLE, "k_deep")
    r = individual(result.flat_path, BW=140.0)
    assert r[sid] == pytest.approx(4e-3)


def test_a_custom_parameter_simulates_as_a_sbmlutils_parameter(
    tmp_path: Path,
) -> None:
    """A `ModelParameter` and a sbmlutils `Parameter` of one value give one model."""
    as_model = build(
        _definition([_deep(ModelParameter("k_deep", 2e-3, U.per_min), _KP_DEEP)]),
        tmp_path / "model",
    )
    as_sbml = build(
        _definition(
            [
                _deep(
                    Parameter("k_deep", 2e-3, U.per_min, constant=True),
                    Parameter("Kp_deep", 1.5, U.dimensionless, constant=True),
                )
            ]
        ),
        tmp_path / "sbml",
    )
    values = []
    for result in (as_model, as_sbml):
        r = load(result)
        r[ids.iv_dose("drg")] = 10.0
        r.simulate(0, DAY, 20)
        values.append(r["[MU__Cdeep_drg]"])
    assert values[0] > 0
    assert values[0] == pytest.approx(values[1], rel=1e-9)


@pytest.mark.parametrize(
    ("parameters", "match"),
    [
        (
            [ModelParameter("k__deep", 2e-3, U.per_min)],
            "DRGDEEP: parameter id 'k__deep' contains '__'",
        ),
        (
            [ModelParameter("1k", 2e-3, U.per_min)],
            "DRGDEEP: parameter id '1k' is not an SBML SId",
        ),
        (
            [
                ModelParameter("k_deep", 2e-3, U.per_min),
                ModelParameter("k_deep", 3e-3, U.per_min),
            ],
            "parameter 'k_deep' of organ 'mu' is defined twice",
        ),
        (
            [ModelParameter("k_deep", Value(2e-3, U.per_hr), U.per_min)],
            "DRGDEEP: parameter 'k_deep': value must be in",
        ),
    ],
)
def test_wrong_custom_parameters_are_refused(
    parameters: list[ModelParameter], match: str
) -> None:
    """An id which is no SId or has `__`, an id twice in an organ, a value in another unit."""
    assert match in _problems(_definition([_deep(*parameters, _KP_DEEP)]))


def test_a_custom_parameter_twice_in_two_reactions_of_an_organ_is_refused() -> None:
    """Two custom reactions of an organ define a parameter once; the other reads it."""
    first = _deep(ModelParameter("k_deep", 2e-3, U.per_min), _KP_DEEP)
    second = dataclasses.replace(
        first,
        sid="DRGDEEP2",
        species=[],
        equation="Cti_drg -> Cdeep_drg",
        formula="k_deep * Vti * Cti_drg",
        parameters=[ModelParameter("k_deep", 2e-3, U.per_min)],
    )
    problems = _problems(_definition([first, second]))
    assert "parameter 'k_deep' of organ 'mu' is defined twice" in problems
    reads = dataclasses.replace(second, parameters=[])
    _definition([first, reads]).validate()


def test_custom_parameters_of_two_organs_are_two_parameters(tmp_path: Path) -> None:
    """The same id in two organs is a constant of each organ, `<ORGAN>__<id>`."""
    muscle = _deep(ModelParameter("k_deep", 2e-3, U.per_min), _KP_DEEP)
    skin = dataclasses.replace(
        muscle,
        sid="DRGDEEPSK",
        organ=Organ.SKIN,
        parameters=[ModelParameter("k_deep", 5e-3, U.per_min), _KP_DEEP],
    )
    flat = _flat(build(_definition([muscle, skin]), tmp_path))
    assert flat.getParameter("MU__k_deep").getValue() == 2e-3
    assert flat.getParameter("SK__k_deep").getValue() == 5e-3


# --- shared rates (D4) --------------------------------------------------------


def _uptake(sid: str, vmax: float | Value | RateOf, km: float = 0.05) -> Transport:
    """A basolateral uptake of the drug into the liver."""
    return Transport(
        sid,
        Organ.LIVER,
        "drg",
        Direction.INFLUX,
        Membrane.BASOLATERAL,
        MM(vmax=vmax, km=km),
        protein="SLCO1B1" if sid.endswith("1B1") else "SLCO1B3",
    )


_RATIO = Value(
    0.6701,
    U.dimensionless,
    source="ModelChoice",
    notes="The abundance times the in vitro Vmax ratio.",
)


def test_rate_of_writes_the_capacity_as_a_ratio_of_another_process(
    tmp_path: Path,
) -> None:
    """`RateOf` gives the capacity `Vmax_b = factor*Vmax_a`, local ids of the organ.

    The ratio is the constant `<b>_Vmax_ratio` with its provenance, the
    capacity its initial assignment, so a change of the capacity of the other
    process before the simulation (a fit) moves both.
    """
    d = _definition(
        [
            _uptake(
                "UP_1B1",
                _fitted(0.02, "lognormal(-3.9, 0.2, 0.001, 1.0)", U.mmole_per_min_l),
            ),
            _uptake("UP_1B3", RateOf("UP_1B1", _RATIO)),
        ]
    )
    result = build(d, tmp_path)
    md = _comp_definition(result, Organ.LIVER)
    assert _assignment(md, "UP_1B3_Vmax") == "UP_1B3_Vmax_ratio * UP_1B1_Vmax"
    flat = _flat(result)
    ratio = ids.flat(Organ.LIVER, ids.rate_ratio("UP_1B3", "Vmax"))
    assert ratio == "LI__UP_1B3_Vmax_ratio"
    assert flat.getParameter(ratio).getValue() == pytest.approx(0.6701)
    assert "The abundance times the in vitro Vmax ratio." in _notes(flat, ratio)
    assert _distribution(flat, "LI__UP_1B3_Vmax") is None
    r = individual(result.flat_path)
    assert r["LI__UP_1B3_Vmax"] == pytest.approx(0.6701 * 0.02)
    r = individual(result.flat_path, LI__UP_1B1_Vmax=0.05)
    assert r["LI__UP_1B3_Vmax"] == pytest.approx(0.6701 * 0.05)


def test_rate_of_shares_a_first_order_clearance(tmp_path: Path) -> None:
    """`RateOf` of the `k` of a first order process: `k_b = factor*k_a`."""
    d = _definition(
        [
            Metabolism("M1", Organ.LIVER, "drg", {"drg_m1": 1.0}, FirstOrder(k=0.2)),
            Metabolism(
                "M2", Organ.LIVER, "drg", {"drg_m2": 1.0}, FirstOrder(RateOf("M1", 0.5))
            ),
        ]
    )
    r = individual(build(d, tmp_path).flat_path)
    assert r["LI__M2_k"] == pytest.approx(0.1)
    assert r["LI__M2_k_ratio"] == pytest.approx(0.5)


@pytest.mark.parametrize("liver", ["zonated", "standard"])
def test_rate_of_in_the_variants(tmp_path: Path, liver: str) -> None:
    """The zonated liver and the intestine share a capacity as the standard organ does."""
    gut = Metabolism(
        "GU1",
        Organ.SMALL_INTESTINE,
        "drg",
        {"drg_m1": 1.0},
        Hill(vmax=1e-3, km=0.01, n=2.0),
        protein="CYP3A4",
    )
    d = _definition(
        [
            _uptake("UP_1B1", 0.02),
            _uptake("UP_1B3", RateOf("UP_1B1", 0.5)),
            gut,
            dataclasses.replace(
                gut, sid="GU2", kinetics=Hill(vmax=RateOf("GU1", 2.0), km=0.02, n=1.0)
            ),
        ],
        variants=Variants(organs={Organ.LIVER: liver}),
    )
    r = individual(build(d, tmp_path).flat_path)
    assert r["LI__UP_1B3_Vmax"] == pytest.approx(0.01)
    assert r["GU__GU2_Vmax"] == pytest.approx(2e-3)


def _two_processes(second: Metabolism | Transport, *others: object) -> PBPKDefinition:
    """A definition of an uptake of the liver and a second process."""
    return _definition([_uptake("UP_1B1", 0.02), second, *others])  # ty: ignore[invalid-argument-type]


@pytest.mark.parametrize(
    ("second", "match"),
    [
        (_uptake("UP_1B3", RateOf("UP_X", 0.5)), "UP_1B3: RateOf 'UP_X' names no"),
        (
            Transport(
                "UP_KI",
                Organ.KIDNEY,
                "drg",
                Direction.INFLUX,
                Membrane.BASOLATERAL,
                MM(vmax=RateOf("UP_1B1", 0.5), km=0.05),
            ),
            "UP_KI: RateOf 'UP_1B1' is a process of organ 'li', not of 'ki'",
        ),
        (
            Metabolism("M1", Organ.LIVER, "drg", {}, FirstOrder(RateOf("UP_1B1", 0.5))),
            "M1: RateOf 'UP_1B1' has no first order constant k",
        ),
        (
            _uptake("UP_1B3", RateOf("UP_1B1", -0.5)),
            "UP_1B3: RateOf factor must be positive",
        ),
        (_uptake("UP_1B3", RateOf("UP_1B3", 0.5)), "UP_1B3: RateOf 'UP_1B3' is"),
        (
            _uptake("UP_1B3", RateOf("UP_1B1", Value(0.5, U.per_min))),
            "UP_1B3: RateOf factor must be in dimensionless",
        ),
    ],
)
def test_rate_of_is_refused_for_another_organ_or_a_missing_process(
    second: Metabolism | Transport, match: str
) -> None:
    """The other process must be a metabolism or transport of the organ with the same constant.

    A missing process, a process of another organ, another kind of constant
    (a Vmax for a k, a k of the whole organ for one per litre), a process
    whose constant is itself a ratio (or the process itself) and a factor
    which is not positive or not dimensionless are refused.
    """
    assert match in _problems(_two_processes(second))


def test_rate_of_a_constant_of_another_scaling_is_refused() -> None:
    """A clearance of the whole organ is no ratio of one per litre of organ."""
    d = _definition(
        [
            Metabolism("M1", Organ.LIVER, "drg", {}, FirstOrder(k=0.2)),
            Metabolism(
                "M2",
                Organ.LIVER,
                "drg",
                {},
                FirstOrder(k=RateOf("M1", 0.5), per_volume=False),
            ),
        ]
    )
    problems = _problems(d)
    assert (
        "M2: RateOf 'M1' has no first order constant k of the whole organ" in problems
    )


def test_rate_of_a_ratio_is_refused() -> None:
    """A chain of ratios is refused: the other process states its own constant."""
    d = _two_processes(
        _uptake("UP_1B3", RateOf("UP_1B1", 0.5)),
        _uptake("UP_2B1", RateOf("UP_1B3", 0.5)),
    )
    assert "UP_2B1: RateOf 'UP_1B3' is itself a ratio" in _problems(d)


def test_rate_of_a_custom_reaction_is_refused() -> None:
    """A custom reaction has no constant of the library to share."""
    d = _definition(
        [
            _deep(ModelParameter("k_deep", 2e-3, U.per_min), _KP_DEEP),
            Transport(
                "UP_MU",
                Organ.MUSCLE,
                "drg",
                Direction.INFLUX,
                Membrane.BASOLATERAL,
                MM(vmax=RateOf("DRGDEEP", 0.5), km=0.05),
            ),
        ]
    )
    assert "UP_MU: RateOf 'DRGDEEP' is no metabolism or transport" in _problems(d)


# --- product fractions (D5) ---------------------------------------------------


def _split(
    products: Mapping[str, object], kinetics: FirstOrder | MM | Hill | None = None
) -> Metabolism:
    """The metabolism of the drug in the liver into the products."""
    return Metabolism(
        "DRG2M",
        Organ.LIVER,
        "drg",
        products,  # ty: ignore[invalid-argument-type]
        kinetics or FirstOrder(k=0.5),
        protein="CYP3A4",
    )


_F_M1 = _fitted(0.3, "normal(0.3, 0.05, 0.01, 0.6)")


@pytest.mark.parametrize(
    ("products", "match"),
    [
        (
            {
                "drg_m1": _fitted(0.5, "normal(0.5, 0.1, 0.1, 0.7)"),
                "drg_m2": _fitted(0.3, "normal(0.3, 0.1, 0.1, 0.4)"),
                "drg_m3": REST,
            },
            "DRG2M: the product fractions can sum to 1.1",
        ),
        (
            {"drg_m1": _fitted(0.5, "normal(0.5, 0.1)"), "drg_m3": REST},
            "DRG2M: the product fractions can sum to inf > 1",
        ),
        (
            {"drg_m1": Value(0.7, U.dimensionless), "drg_m2": 0.5, "drg_m3": REST},
            "DRG2M: the product fractions can sum to 1.2",
        ),
        (
            {"drg_m1": _fitted(0.5, "uniform(0.0, 1.5)"), "drg_m3": REST},
            "DRG2M: the product fractions can sum to 1.5",
        ),
        (
            {"drg_m1": _fitted(0.5, "normal(0.5, 0.1, -0.1, 0.6)"), "drg_m3": REST},
            "DRG2M: product fraction of 'drg_m1' can be -0.1",
        ),
        (
            {"drg_m1": 1.0, "drg_m3": REST},
            "DRG2M: the other fractions sum to 1, the remainder REST of 'drg_m3' is 0",
        ),
        (
            {"drg_m1": REST, "drg_m3": REST},
            "DRG2M: more than one product is the remainder REST",
        ),
        ({"drg_m1": REST}, "DRG2M: the remainder REST needs another product"),
        (
            {"drg_m1": _F_M1, "drg_m2": 0.7},
            "DRG2M: the fitted product fraction of 'drg_m1' needs a product which "
            "is the remainder REST",
        ),
        (
            {"drg_m1": Formula("0.3 dimensionless", U.dimensionless), "drg_m3": REST},
            "DRG2M: product fraction of 'drg_m1' is a Formula",
        ),
        (
            {"drg_m1": Value(0.3, U.per_min), "drg_m3": REST},
            "DRG2M: product fraction of 'drg_m1' must be in dimensionless",
        ),
        (
            {"drg_m1": Value(0.3, U.dimensionless), "drg_m2": 0.6},
            "DRG2M: product fractions sum to 0.9, not 1",
        ),
    ],
)
def test_fractions_which_can_exceed_one_are_refused(
    products: Mapping[str, object], match: str
) -> None:
    """The fractions of a metabolism never sum above 1, the remainder is never negative.

    The range of a fitted fraction is its value and the support of the
    distribution of its draws (the bounds of the fit); without a remainder
    the fractions are numbers or values which sum to 1, a fitted fraction
    needs the remainder `REST`; a formula cannot be bounded.
    """
    assert match in _problems(_definition([_split(products)]))


def test_the_fractions_of_a_metabolism_are_parameters(tmp_path: Path) -> None:
    """Values and the remainder are the constants `f_<metabolism>_<product>` of the organ.

    A fitted fraction keeps the uncertainty of its fit, the remainder is the
    initial assignment `1 - f_a - f_b`; the metabolism consumes the substrate
    and a reaction per product forms the product with its fraction of the
    rate.
    """
    d = _definition(
        [
            _split(
                {
                    "drg_m1": _F_M1,
                    "drg_m2": Value(0.2, U.dimensionless, source="ModelChoice"),
                    "drg_m3": REST,
                }
            )
        ]
    )
    result = build(d, tmp_path)
    md = _comp_definition(result, Organ.LIVER)
    f1, f2, f3 = (ids.product_fraction("DRG2M", f"drg_m{k}") for k in (1, 2, 3))
    assert f1 == "f_DRG2M_drg_m1"
    assert _assignment(md, f3) == f"1 dimensionless - {f1} - {f2}"
    reaction = md.getReaction("DRG2M")
    assert [s.getSpecies() for s in reaction.getListOfReactants()] == ["Cti_drg"]
    assert reaction.getNumProducts() == 0
    formed = md.getReaction(ids.product_reaction("DRG2M", "drg_m1"))
    assert formed.getId() == "DRG2M_to_drg_m1"
    assert [s.getSpecies() for s in formed.getListOfProducts()] == ["Cpl_drg_m1"]
    assert formed.getNumReactants() == 0
    assert (
        libsbml.formulaToL3String(formed.getKineticLaw().getMath()) == f"{f1} * DRG2M"
    )
    flat = _flat(result)
    assert flat.getParameter(f"LI__{f1}").getValue() == 0.3
    assert _distribution(flat, f"LI__{f1}") == "normal(0.3, 0.05, 0.01, 0.6)"
    assert flat.getParameter(f"LI__{f2}").getValue() == 0.2
    r = individual(result.flat_path)
    assert r[f"LI__{f3}"] == pytest.approx(0.5)


def _formed(r: object, sid: str) -> float:
    """The amount of a metabolite in the body (species of the id scheme)."""
    return total_amount(r, sid)  # ty: ignore[invalid-argument-type]


@pytest.mark.parametrize(
    "kinetics", [FirstOrder(k=50.0), MM(vmax=5.0, km=0.01)], ids=["first", "mm"]
)
def test_a_metabolism_with_fitted_fractions_keeps_the_mass_balance(
    tmp_path: Path, kinetics: FirstOrder | MM
) -> None:
    """The drug and its products keep the amount of the dose, in the ratio of the fractions.

    Also after a change of a fitted fraction before the simulation (the
    remainder follows).
    """
    d = _definition(
        [
            _split(
                {
                    "drg_m1": _F_M1,
                    "drg_m2": Value(0.2, U.dimensionless),
                    "drg_m3": REST,
                },
                kinetics,
            )
        ]
    )
    result = build(d, tmp_path)
    for f1, expected in ((None, 0.3), (0.5, 0.5)):
        changes = {} if f1 is None else {"LI__f_DRG2M_drg_m1": f1}
        r = individual(result.flat_path, **changes)
        r.integrator.absolute_tolerance = 1e-14
        r.integrator.relative_tolerance = 1e-10
        r[ids.iv_dose("drg")] = 10.0
        mmole = dose_amount(r, "drg")
        r.simulate(0, DAY, 50)
        substances = ("drg", "drg_m1", "drg_m2", "drg_m3")
        total = sum(total_amount(r, s) + dose_amount(r, s) for s in substances)
        assert total == pytest.approx(mmole, rel=1e-6)
        metabolized = sum(_formed(r, s) for s in substances[1:])
        assert metabolized > 0.01 * mmole
        assert _formed(r, "drg_m1") / metabolized == pytest.approx(expected, rel=1e-6)
        assert _formed(r, "drg_m2") / metabolized == pytest.approx(0.2, rel=1e-6)
        assert _formed(r, "drg_m3") / metabolized == pytest.approx(
            0.8 - expected, rel=1e-6
        )


@pytest.mark.parametrize("liver", ["zonated", "standard"])
def test_fitted_fractions_in_the_variants(tmp_path: Path, liver: str) -> None:
    """The zonated liver and the mucosa of the intestine form the products by the fractions."""
    products = {"drg_m1": _F_M1, "drg_m3": REST}
    d = dataclasses.replace(
        _definition(
            [
                _split(products),
                Metabolism(
                    "DRG2MGU",
                    Organ.SMALL_INTESTINE,
                    "drg",
                    products,
                    FirstOrder(k=0.5),
                    protein="CYP3A4",
                ),
            ],
            variants=Variants(organs={Organ.LIVER: liver}),
        ),
        routes=[Route.iv("drg"), Route.oral("drg")],
    )
    result = build(d, tmp_path)
    r = load(result)
    r.integrator.absolute_tolerance = 1e-14
    r[ids.iv_dose("drg")] = 10.0
    r[ids.po_dose("drg")] = 10.0
    mmole = dose_amount(r, "drg")
    r.simulate(0, DAY, 50)
    substances = ("drg", "drg_m1", "drg_m3")
    total = sum(total_amount(r, s) + dose_amount(r, s) for s in substances)
    assert total == pytest.approx(mmole, rel=1e-6)
    metabolized = _formed(r, "drg_m1") + _formed(r, "drg_m3")
    assert _formed(r, "drg_m1") / metabolized == pytest.approx(0.3, rel=1e-6)
    flat = _flat(result)
    assert flat.getParameter("GU__f_DRG2MGU_drg_m1") is not None
    zone = "_z1" if liver == "zonated" else ""
    assert flat.getReaction(f"LI__DRG2M_to_drg_m1{zone}") is not None


# --- byte identity ------------------------------------------------------------


@pytest.mark.parametrize("liver", ["standard", "zonated"])
def test_the_chemistry_without_the_options_is_that_of_0c538ca(
    tmp_path: Path, liver: str
) -> None:
    """Numbers as products, own capacities and sbmlutils parameters build the SBML of 0c538ca."""
    assert_reference_sbml(build(chemistry_reference_definition(liver), tmp_path))
