"""Tests of the conversions in the lumen and the regional absorption profile.

`LumenConversion` converts a substance into another in chosen lumen segments
of the segmented intestine, first order on the dissolved concentration;
`RegionalAbsorption` scales the absorption `ka*f_area` of chosen segments
of a substance. Neither creates an element of the model unless it is set.
"""

import dataclasses
import math
from pathlib import Path
from typing import TYPE_CHECKING

import libsbml
import numpy as np
import pytest

from pkpdsim.pbpk import (
    Excretion,
    FirstOrder,
    FirstOrderDissolution,
    LumenConversion,
    Organ,
    PBPKDefinition,
    Pool,
    RegionalAbsorption,
    Route,
    SegmentedIntestine,
    Substance,
    U,
    Value,
    Variants,
    build,
    ids,
)
from pkpdsim.pbpk.chemistry import Direction, Membrane, Transport
from pkpdsim.pbpk.data import value
from pkpdsim.pbpk.variants.segments import DATA, SEGMENTS

from .conftest import (
    assert_reference_sbml,
    chemistry_reference_definition,
    dose_amount,
    load,
    reference_definition,
    simple_drug,
    total_amount,
)

if TYPE_CHECKING:
    from collections.abc import Callable

    import roadrunner

S, M = "drg", "drg_m"
DAY = 24 * 60.0
SMALL = ("du", "jej1", "jej2", "ile1", "ile2")
COLON = ("cae", "asc", "tra", "des", "sig", "rec")


def _gu(local: str) -> str:
    """The flat id of a local id of the intestine."""
    return ids.flat(Organ.SMALL_INTESTINE, local)


def _transit(seg: str) -> float:
    """The transit time (min, fasted) of a lumen segment of the data."""
    return value(DATA, SEGMENTS[seg], "transit_time").value


def _definition(
    sid: str,
    intestine: SegmentedIntestine,
    routes: list[Route] | None = None,
    chemistry: tuple[Transport, ...] = (),
) -> PBPKDefinition:
    """The drug (no elimination) and its product (renal excretion), the oral solution of the drug."""
    return PBPKDefinition(
        sid=sid,
        substances=[
            dataclasses.replace(simple_drug(), excretion=Excretion()),
            Substance(M, name="product", mr=316.0, excretion=Excretion(urine=True)),
        ],
        chemistry=list(chemistry),
        variants=Variants(organs={Organ.SMALL_INTESTINE: intestine}),
        routes=routes if routes is not None else [Route.oral_solution(S)],
    )


def _sbml_ids(path: Path) -> set[str]:
    """The ids of the parameters, species and reactions of a flat model."""
    model = libsbml.readSBMLFromFile(str(path)).getModel()
    return {
        x.getId()
        for xs in (
            model.getListOfParameters(),
            model.getListOfSpecies(),
            model.getListOfReactions(),
        )
        for x in xs
    }


def _law(path: Path, rid: str) -> str:
    """The rate law of a reaction of the flat model as L3 text."""
    model = libsbml.readSBMLFromFile(str(path)).getModel()
    return libsbml.formulaToL3String(model.getReaction(rid).getKineticLaw().getMath())


# --- the conversion in the lumen --------------------------------------------


def _converted(r: "roadrunner.RoadRunner", route: str, dose_mg: float) -> float:
    """The fraction of an oral dose of the drug converted, with no absorption.

    Neither substance is absorbed (`ka = 0`), so after a week both are in
    the feces: the product there per dose (in mmole).
    """
    r[route] = dose_mg
    r.simulate(0, 7 * DAY, 50)
    return r[ids.feces(M)] / (dose_mg / r[ids.mr(S)])


_NO_ABSORPTION = {S: 0.0, M: 0.0}
"""The absorption of neither substance, the conversion is all that acts in the lumen."""


def test_the_converted_fraction_of_a_solution_in_the_stomach(tmp_path: Path) -> None:
    """A solution in the stomach: `f = k*T/(1 + k*T)` with the gastric transit `T`.

    The stomach does not absorb, the dissolved drug leaves it with the
    transit `V*C/T` and converts with `k/60*V*C` (`k` in 1/hr), so the
    converted fraction is `k/60*T/(1 + k/60*T)`; in the fed state with the
    transit time of the fed stomach.
    """
    k = 0.5
    intestine = SegmentedIntestine(
        ka=_NO_ABSORPTION, conversions=[LumenConversion(S, M, k, ["st"])]
    )
    result = build(_definition("conv_sol", intestine), tmp_path)
    fed_transit = value(DATA, "stomach", "transit_time_fed").value
    for fed, t in ((0.0, _transit("st")), (1.0, fed_transit)):
        r = load(result)
        r[ids.FED] = fed
        assert r[ids.lumen_conversion_rate(S, M)] == k
        kt = k / 60 * t
        assert _converted(r, ids.sol_dose(S), 10.0) == pytest.approx(
            kt / (1 + kt), rel=1e-5
        )


def test_the_converted_fraction_over_several_segments(tmp_path: Path) -> None:
    """Over several segments the drug survives each with `1/(1 + k*T_seg)`.

    The converted fraction of a solution is `1 - prod(1/(1 + k/60*T_seg))`
    over the segments of the conversion (no absorption).
    """
    k, segments = 0.2, ("st", "du", "jej1", "cae")
    intestine = SegmentedIntestine(
        ka=_NO_ABSORPTION, conversions=[LumenConversion(S, M, k, segments)]
    )
    r = load(build(_definition("conv_segments", intestine), tmp_path))
    survives = math.prod(1 / (1 + k / 60 * _transit(seg)) for seg in segments)
    assert _converted(r, ids.sol_dose(S), 10.0) == pytest.approx(1 - survives, rel=1e-5)


def test_the_converted_fraction_of_a_tablet_in_the_stomach(tmp_path: Path) -> None:
    """A tablet: only the drug dissolved in the stomach converts there.

    With the first order dissolution `kdis` (1/hr) the fraction which
    dissolves in the stomach is `g = kdis/60*T/(1 + kdis/60*T)`, the
    undissolved rest moves on; the converted fraction is
    `g * k/60*T/(1 + k/60*T)`.
    """
    k, kdis = 0.5, 6.0
    intestine = SegmentedIntestine(
        ka=_NO_ABSORPTION, conversions=[LumenConversion(S, M, k, ("st",))]
    )
    d = _definition(
        "conv_tablet", intestine, routes=[Route.oral(S, FirstOrderDissolution(k=kdis))]
    )
    r = load(build(d, tmp_path))
    t = _transit("st")
    g = kdis / 60 * t / (1 + kdis / 60 * t)
    kt = k / 60 * t
    assert _converted(r, ids.po_dose(S), 10.0) == pytest.approx(
        g * kt / (1 + kt), rel=1e-5
    )


def test_the_conversion_conserves_the_amount(tmp_path: Path) -> None:
    """Dose = both substances in body, lumen, mucosa, feces and urine.

    Both absorbed, the conversion in every segment, tablet and solution; the
    product is excreted into the urine.
    """
    intestine = SegmentedIntestine(
        ka={S: 1.0, M: 0.5},
        conversions=[LumenConversion(S, M, 0.3, ids.GUT_SEGMENTS)],
    )
    d = _definition(
        "conv_balance", intestine, routes=[Route.oral(S), Route.oral_solution(S)]
    )
    r = load(build(d, tmp_path))
    mr = r[ids.mr(S)]
    r[ids.po_dose(S)] = 10.0
    r[ids.sol_dose(S)] = 5.0
    r.simulate(0, 120, 10)
    assert r[f"[{_gu(ids.lumen(M, 'du'))}]"] > 0
    r.simulate(120, 3 * DAY, 50)
    total = sum(total_amount(r, s) + dose_amount(r, s) for s in (S, M))
    assert total == pytest.approx(15.0 / mr, rel=1e-6)
    assert r[ids.urine(M)] > 0
    assert r[ids.feces(M)] > 0


def test_the_conversion_reactions(tmp_path: Path) -> None:
    """One reaction per segment, `Clumen_<seg>_<s> -> Clumen_<seg>_<p>`, `k/60*V*C`.

    The rate constant `kconv_<s>_<p>` (1/hr) is a parameter of the top model
    with its provenance; the product without an oral route gets lumen
    species (and no absorption, its `ka` is 0).
    """
    k = Value(0.02, U.per_hr, source="ModelChoice", notes="a test value")
    intestine = SegmentedIntestine(conversions=[LumenConversion(S, M, k, ("st", "du"))])
    result = build(_definition("conv_ids", intestine), tmp_path)
    found = _sbml_ids(result.flat_path)
    rate = ids.lumen_conversion_rate(S, M)
    assert rate == "kconv_drg_drg_m"
    assert rate in found
    reactions = {
        seg: _gu(ids.lumen_conversion_reaction(S, M, seg)) for seg in ("st", "du")
    }
    assert reactions["st"] == "GU__J_conv_st_drg_drg_m"
    assert set(reactions.values()) <= found
    assert not {i for i in found if "J_conv_" in i} - set(reactions.values())
    assert _gu(ids.lumen(M, "rec")) in found
    model = libsbml.readSBMLFromFile(str(result.flat_path)).getModel()
    for seg, rid in reactions.items():
        reaction = model.getReaction(rid)
        assert [x.getSpecies() for x in reaction.getListOfReactants()] == [
            _gu(ids.lumen(S, seg))
        ]
        assert [x.getSpecies() for x in reaction.getListOfProducts()] == [
            _gu(ids.lumen(M, seg))
        ]
        assert reaction.getSBOTermID() == "SBO:0000176"
    parameter = model.getParameter(rate)
    assert parameter.getUnits() == "per_hr"
    assert "a test value" in parameter.getNotesString()
    r = load(result)
    r[ids.sol_dose(S)] = 10.0
    r.simulate(0, 30, 10)
    for seg, rid in reactions.items():
        expected = (
            0.02 / 60 * r[_gu(ids.lumen_volume(seg))] * r[f"[{_gu(ids.lumen(S, seg))}]"]
        )
        assert expected > 0
        assert r[rid] == pytest.approx(expected, rel=1e-9)
    assert r[ids.absorption(M)] == 0.0


def test_a_substrate_which_reaches_the_lumen_by_the_bile_converts(
    tmp_path: Path,
) -> None:
    """The substrate of a conversion may reach the lumen by its biliary excretion.

    A chain: the product of the first conversion is the substrate of the
    second, which reaches the lumen only through the first.
    """
    bile = Transport(
        "DRGMBILE",
        Organ.LIVER,
        M,
        Direction.EFFLUX,
        Membrane.APICAL,
        FirstOrder(k=0.1),
    )
    intestine = SegmentedIntestine(conversions=[LumenConversion(M, S, 0.1, ("du",))])
    d = _definition("conv_bile", intestine, routes=[Route.iv(M)], chemistry=(bile,))
    found = _sbml_ids(build(d, tmp_path).flat_path)
    assert _gu(ids.lumen_conversion_reaction(M, S, "du")) in found
    assert _gu(ids.lumen(S, "du")) in found


@pytest.mark.parametrize(
    ("conversions", "message"),
    [
        (lambda: [LumenConversion(S, S, 0.1, ("st",))], "converts 'drg' into itself"),
        (lambda: [LumenConversion(S, "xyz", 0.1, ("st",))], "'xyz' is no substance"),
        (lambda: [LumenConversion(S, M, 0.1, ())], "no segments"),
        (lambda: [LumenConversion(S, M, 0.1, ("st", "st"))], "more than once"),
        (lambda: [LumenConversion(S, M, 0.1, ("colon",))], "unknown segments"),
        (lambda: [LumenConversion(S, M, -0.1, ("st",))], "must not be negative"),
        (
            lambda: [LumenConversion(S, M, Value(0.1, U.per_min), ("st",))],
            "must be in 1/hr",
        ),
        (
            lambda: [
                LumenConversion(S, M, 0.1, ("st",)),
                LumenConversion(S, M, 0.2, ("du",)),
            ],
            "twice",
        ),
        (
            lambda: [LumenConversion(M, S, 0.1, ("st",))],
            "'drg_m' does not reach the lumen",
        ),
    ],
)
def test_an_invalid_conversion_is_refused(
    conversions: "Callable[[], list[LumenConversion]]",
    message: str,
    tmp_path: Path,
) -> None:
    """A conversion of the same substance, of an unknown substance or segment, a negative k, a twin, or a substrate which never reaches the lumen."""
    with pytest.raises(ValueError, match=message):
        intestine = SegmentedIntestine(conversions=conversions())
        build(_definition("conv_bad", intestine), tmp_path)


# --- the regional absorption ------------------------------------------------


def _absorption_rate(r: "roadrunner.RoadRunner", sid: str, seg: str) -> float:
    """`ka/60*f_area*V*C` of a segment without factors."""
    return (
        r[ids.absorption(sid)]
        / 60
        * r[_gu(ids.area_factor(seg))]
        * r[_gu(ids.lumen_volume(seg))]
        * r[f"[{_gu(ids.lumen(sid, seg))}]"]
    )


def test_the_regional_profile_scales_the_absorption_of_its_segments(
    tmp_path: Path,
) -> None:
    """A factor `f_ka_<seg>_<s>` multiplies the absorption of its segment only."""
    profile = RegionalAbsorption({"ile2": 0.5, "asc": 0.2, "rec": 0.0})
    intestine = SegmentedIntestine(regional=profile)
    result = build(_definition("regional", intestine), tmp_path)
    found = _sbml_ids(result.flat_path)
    factors = {i for i in found if i.startswith("f_ka_")}
    assert factors == {
        ids.regional_absorption_factor(S, seg) for seg in ("ile2", "asc", "rec")
    }
    assert ids.regional_absorption_factor(S, "asc") == "f_ka_asc_drg"
    r = load(result)
    r[ids.sol_dose(S)] = 10.0
    r.simulate(0, 8 * 60.0, 20)
    for seg, factor in (("jej1", 1.0), ("ile2", 0.5), ("asc", 0.2), ("rec", 0.0)):
        rate = _absorption_rate(r, S, seg)
        assert r[_gu(ids.absorption_reaction(S, seg))] == pytest.approx(
            factor * rate, rel=1e-9, abs=1e-18
        )
    assert r[f"[{_gu(ids.lumen(S, 'asc'))}]"] > 0


def test_a_profile_of_ones_equals_the_default(tmp_path: Path) -> None:
    """With every factor 1 the model simulates as without the profile."""
    selections = [
        "time",
        *(f"[{_gu(ids.lumen(S, seg))}]" for seg in SMALL + COLON),
        f"[{ids.pool_species(Pool.VENOUS, S)}]",
    ]
    runs = []
    for name, regional in (
        ("default", None),
        ("ones", RegionalAbsorption(dict.fromkeys(SMALL + COLON, 1.0))),
    ):
        d = _definition(f"regional_{name}", SegmentedIntestine(regional=regional))
        r = load(build(d, tmp_path / name))
        r[ids.sol_dose(S)] = 10.0
        runs.append(np.asarray(r.simulate(0, DAY, 100, selections=selections)))
    np.testing.assert_allclose(runs[1], runs[0], rtol=1e-8, atol=1e-16)


def test_a_profile_per_substance(tmp_path: Path) -> None:
    """A mapping gives the profile to the substances in it, the others keep `ka*f_area`."""
    intestine = SegmentedIntestine(
        regional={S: RegionalAbsorption({"cae": Value(0.3, U.dimensionless)})}
    )
    d = _definition(
        "regional_map", intestine, routes=[Route.oral_solution(S), Route.oral(M)]
    )
    result = build(d, tmp_path)
    found = _sbml_ids(result.flat_path)
    assert ids.regional_absorption_factor(S, "cae") in found
    assert not {i for i in found if i.startswith("f_ka_") and i.endswith(f"_{M}")}
    assert "f_ka_cae_drg" in _law(
        result.flat_path, _gu(ids.absorption_reaction(S, "cae"))
    )
    assert "f_ka_" not in _law(result.flat_path, _gu(ids.absorption_reaction(M, "cae")))


def test_a_profile_and_the_colon_factor_multiply_outside_the_colon(
    tmp_path: Path,
) -> None:
    """`f_ka_colon` with a profile of the small intestine: both act, each on its segments."""
    intestine = SegmentedIntestine(
        f_ka_colon=0.4, regional=RegionalAbsorption({"ile1": 0.5})
    )
    path = build(_definition("regional_colon", intestine), tmp_path).flat_path
    assert "f_ka_ile1_drg" in _law(path, _gu(ids.absorption_reaction(S, "ile1")))
    law = _law(path, _gu(ids.absorption_reaction(S, "asc")))
    assert "f_ka_colon_drg" in law
    assert "f_ka_asc" not in law


@pytest.mark.parametrize(
    ("make", "message"),
    [
        (lambda: RegionalAbsorption({"st": 0.5}), "stomach does not absorb"),
        (lambda: RegionalAbsorption({"colon": 0.5}), "unknown segments"),
        (lambda: RegionalAbsorption({}), "no segments"),
        (lambda: RegionalAbsorption({"asc": -0.1}), "must not be negative"),
        (
            lambda: RegionalAbsorption({"asc": Value(0.1, U.per_hr)}),
            "must be in dimensionless",
        ),
    ],
)
def test_an_invalid_profile_is_refused(
    make: "Callable[[], RegionalAbsorption]", message: str
) -> None:
    """A segment which does not absorb or is unknown, no segment, a negative factor or a unit."""
    with pytest.raises(ValueError, match=message):
        make()


@pytest.mark.parametrize(
    ("intestine", "message"),
    [
        (
            SegmentedIntestine(regional={M: RegionalAbsorption({"asc": 0.5})}),
            "'drg_m' has no oral route",
        ),
        (
            SegmentedIntestine(regional={"xyz": RegionalAbsorption({"asc": 0.5})}),
            "'xyz' is no substance",
        ),
        (
            SegmentedIntestine(
                f_ka_colon=0.4, regional=RegionalAbsorption({"asc": 0.5})
            ),
            "both f_ka_colon and a regional factor of the large intestine",
        ),
    ],
)
def test_an_invalid_regional_mapping_is_refused(
    intestine: SegmentedIntestine, message: str, tmp_path: Path
) -> None:
    """A key without oral route or no substance, the colon factor twice."""
    with pytest.raises(ValueError, match=message):
        build(_definition("regional_bad", intestine), tmp_path)


# --- the defaults ----------------------------------------------------------


def test_the_options_are_hashable_and_compare_by_value() -> None:
    """The intestine with a profile and conversions keeps its hash and equality."""

    def make() -> SegmentedIntestine:
        return SegmentedIntestine(
            regional={S: RegionalAbsorption({"asc": 0.2, "cae": 0.3})},
            conversions=[LumenConversion(S, M, 0.1, ["st"])],
        )

    assert make() == make()
    assert hash(make()) == hash(make())
    assert make().conversions == (LumenConversion(S, M, 0.1, ("st",)),)
    assert RegionalAbsorption({"asc": 0.2, "cae": 0.3}) == RegionalAbsorption(
        {"cae": 0.3, "asc": 0.2}
    )


def test_unset_options_create_nothing(tmp_path: Path) -> None:
    """Without a profile and a conversion the intestine has no factor, rate or reaction of them."""
    found = _sbml_ids(
        build(_definition("unset", SegmentedIntestine()), tmp_path).flat_path
    )
    assert not {i for i in found if i.startswith(("f_ka_", "kconv_")) or "J_conv_" in i}


def test_unset_options_leave_the_reference_sbml(tmp_path: Path) -> None:
    """The stored references, with the intestine options given explicitly unset."""
    explicit = Variants(
        organs={
            Organ.SMALL_INTESTINE: SegmentedIntestine(regional=None, conversions=())
        }
    )
    assert_reference_sbml(
        build(
            dataclasses.replace(reference_definition(), variants=explicit),
            tmp_path / "base",
        )
    )
    chemistry = chemistry_reference_definition()
    organs = {**chemistry.variants.organs, **explicit.organs}
    assert_reference_sbml(
        build(
            dataclasses.replace(chemistry, variants=Variants(organs=organs)),
            tmp_path / "chemistry",
        )
    )
