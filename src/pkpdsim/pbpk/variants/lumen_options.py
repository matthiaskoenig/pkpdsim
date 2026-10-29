"""The options of the lumen of the segmented intestine: regional absorption and conversions.

`RegionalAbsorption` is the regional absorption profile of a substance, a
factor of the absorption `ka_<s>*f_area_<seg>` per segment or per group of
segments (`colon`, caecum to rectum, the parameter `f_ka_colon_<s>`, of which
`SegmentedIntestine.f_ka_colon` is the shorthand); `LumenConversion` converts
a substance into another in chosen lumen segments, first order. Both are
options of `pkpdsim.pbpk.variants.intestine.SegmentedIntestine`; the builders
of the intestine (`gut`, `gut_lumen`) read them here.
"""

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

from sbmlutils.factory import Sbase
from sbmlutils.metadata import SBO

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.annotation import annotate
from pkpdsim.pbpk.organs import Organ
from pkpdsim.pbpk.placement import value_parameter
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import FrozenMap, Q, check_unit, number
from pkpdsim.pbpk.variants.segments import SEGMENTS, SMALL

if TYPE_CHECKING:
    from pkpdsim.pbpk.definition import PBPKDefinition

__all__ = [
    "ABSORBING",
    "COLON",
    "GROUPS",
    "LARGE",
    "LumenConversion",
    "RegionalAbsorption",
    "check_conversions",
    "conversion_objects",
    "converted",
    "segment_factor_objects",
]

ABSORBING = tuple(t for t in ids.GUT_SEGMENTS if t != "st")
"""The lumen segments which absorb: every segment but the stomach."""

LARGE = tuple(t for t in ABSORBING if t not in SMALL)
"""The lumen segments of the large intestine, caecum to rectum."""

COLON = "colon"
"""The key of the group of the segments of the large intestine in a profile."""

GROUPS: dict[str, tuple[str, ...]] = {COLON: LARGE}
"""The groups of segments with one shared factor of a profile, by key."""

_LABELS = {t: key.replace("_", " ") for t, key in SEGMENTS.items()}
"""The segments in names by token, `upper jejunum`."""

_SBO_PARAMETER = SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER


def _segment_problems(
    where: str, tokens: Sequence[str], allowed: Sequence[str]
) -> list[str]:
    """The problems of the segment tokens of an option: none, twice or unknown."""
    if not tokens:
        return [f"{where}: no segments, give the tokens {list(allowed)}"]
    problems = []
    twice = sorted({t for t in tokens if list(tokens).count(t) > 1})
    if twice:
        problems.append(f"{where}: the segments {twice} more than once")
    unknown = [t for t in tokens if t not in allowed and t not in ids.GUT_SEGMENTS]
    if unknown:
        problems.append(
            f"{where}: unknown segments {unknown}, the segments are {list(allowed)}"
        )
    return problems


def _raise(problems: list[str]) -> None:
    """Raise the problems of an option, if any.

    Raises:
        ValueError: with every problem.
    """
    if problems:
        raise ValueError("; ".join(problems))


def _region(tokens: Iterable[str]) -> Organ | None:
    """The organ of the segments, `None` for the stomach or segments of both regions."""
    tokens = set(tokens)
    if tokens <= set(SMALL):
        return Organ.SMALL_INTESTINE
    if tokens <= set(LARGE):
        return Organ.LARGE_INTESTINE
    return None


@dataclass(frozen=True)
class RegionalAbsorption:
    """The regional absorption profile of a substance in the segmented intestine.

    A factor of a segment multiplies its absorption `ka_<s>*f_area_<seg>`:
    the permeability of the segment relative to the small intestine, which
    `ka_<s>` describes (the area factor holds the mucosal surface per lumen
    volume, Helander 2014). A factor is a parameter of the top model: of a
    segment `f_ka_<seg>_<s>`; of the group `colon` (`GROUPS`, the segments of
    the large intestine, caecum to rectum) one shared `f_ka_colon_<s>`, which
    `SegmentedIntestine.f_ka_colon` is the shorthand of. A segment missing
    from the profile has none (1).

    There is no generic default profile: the colonic absorption relative to
    the small intestine depends on the substance, high for a highly
    permeable drug (BCS class I, relative colonic bioavailability above 70
    %) and low for a drug of low permeability (BCS III and IV, below 50 %;
    Tannergren 2009, 42 drugs), so the default is no factor (1). A substance
    whose regional absorption is measured (an intubation or a capsule
    released in a segment, `Byon2018` for apixaban) takes its own profile.

    Attributes:
        factors: the factor of a segment by its token (`du` to `rec`, not
            the stomach, which does not absorb) or of a group by its key
            (`colon`): a number or a `Value` (dimensionless), not negative.

    Raises:
        ValueError: for no segment, the stomach, an unknown token, a segment
            given by itself and in its group, a negative factor or a `Value`
            which is not dimensionless.
    """

    factors: Mapping[str, Q]

    def __post_init__(self) -> None:
        """Store the factors as a hashable read-only mapping and check them."""
        object.__setattr__(self, "factors", FrozenMap(self.factors))
        where = "RegionalAbsorption"
        allowed = [*ABSORBING, *GROUPS]
        problems = _segment_problems(where, list(self.factors), allowed)
        if "st" in self.factors:
            problems.append(
                f"{where}: the stomach does not absorb, the segments are {allowed}"
            )
        for key, members in GROUPS.items():
            both = [t for t in members if t in self.factors]
            if key in self.factors and both:
                problems.append(
                    f"{where}: the segments {both} have a factor of their own and "
                    f"of the group '{key}', which scale the same absorption"
                )
        for seg, q in self.factors.items():
            x = number(q)
            if x is not None and x < 0:
                problems.append(
                    f"{where}: the factor of '{seg}' must not be negative, got {x}"
                )
            problems += check_unit(where, f"the factor of '{seg}'", q, U.dimensionless)
        _raise(problems)

    @property
    def segments(self) -> frozenset[str]:
        """The segments with a factor of their own (no group)."""
        return frozenset(t for t in self.factors if t not in GROUPS)

    def covers(self, tokens: Iterable[str]) -> bool:
        """Whether the profile has a factor of one of the segments, its own or of its group."""
        covered = set(self.segments)
        for key, members in GROUPS.items():
            if key in self.factors:
                covered |= set(members)
        return bool(covered & set(tokens))


@dataclass(frozen=True)
class LumenConversion:
    """A first order conversion of a substance into another in the lumen.

    In every segment of the conversion the dissolved substrate converts into
    the product, `J_conv_<seg>_<s>_<m> = kconv_<s>_<m>/60*Vlumen_<seg>*Clumen_<seg>_<s>`,
    one mole of product per mole of substrate: a chemical conversion (the
    cyclization of creatine to creatinine in the acid of the stomach), a
    hydrolysis of a prodrug or a bacterial conversion in the colon. The
    undissolved drug of a tablet does not convert. The product gets the
    lumen species, transit and absorption of the intestine (its `ka_<m>` is
    0 without an oral route, unless `SegmentedIntestine.ka` gives it one).

    A solution in the stomach, a conversion in the stomach only: the
    stomach does not absorb, so the converted fraction of the dose is
    `f = k*T/(1 + k*T)` with `k = kconv/60` (1/min) and the transit time of
    the stomach `T` (`T_st`, fasted 22.79 min), and `kconv = 60*f/((1 - f)*T)`.
    An immediate-release tablet with the first order dissolution `Ka_dis_<s>`
    (1/hr, without the factors `f_dis_fed_<s>` and `f_dis_release_<s>`)
    converts only the part which dissolves in the stomach,
    `g = kd*T/(1 + kd*T)` with `kd = Ka_dis/60` (1/min), so
    `f = g*k*T/(1 + k*T)`; the fasted transit time, the fed state has
    `T_st_fed`. Over several segments without absorption the drug survives
    each with `1/(1 + k*T_seg)`.

    Attributes:
        substrate: the id of the converted substance, which reaches the lumen
            (an oral route, biliary excretion, a transport of the intestine
            or the product of another conversion).
        product: the id of the formed substance.
        k: the rate constant (1/hr, the unit of `ka_<s>`), the parameter
            `kconv_<s>_<m>` of the top model, not negative.
        segments: the lumen segments it acts in, tokens of `ids.GUT_SEGMENTS`
            (`st` the stomach to `rec` the rectum).

    Raises:
        ValueError: for a substrate which is the product, no segment, an
            unknown token or one twice, a negative `k` or a `Value` which is
            not in 1/hr.
    """

    substrate: str
    product: str
    k: Q
    segments: Sequence[str]

    def __post_init__(self) -> None:
        """Store the segments as a tuple and check the conversion."""
        object.__setattr__(self, "segments", tuple(self.segments))
        where = f"LumenConversion {self.substrate} -> {self.product}"
        problems = []
        if self.substrate == self.product:
            problems.append(f"{where}: converts '{self.substrate}' into itself")
        problems += _segment_problems(where, self.segments, ids.GUT_SEGMENTS)
        x = number(self.k)
        if x is not None and x < 0:
            problems.append(f"{where}: k must not be negative, got {x}")
        problems += check_unit(where, "k", self.k, U.per_hr)
        _raise(problems)


def check_conversions(
    conversions: Sequence[LumenConversion], sids: Iterable[str]
) -> None:
    """Check the substances of the conversions and that none is given twice.

    Raises:
        ValueError: for a substance which is no substance of the definition,
            or two conversions of the same ids.
    """
    known = frozenset(sids)
    problems = []
    seen: set[str] = set()
    for c in conversions:
        where = f"LumenConversion {c.substrate} -> {c.product}"
        for sid in (c.substrate, c.product):
            if sid not in known:
                problems.append(
                    f"{where}: '{sid}' is no substance of the definition "
                    f"{sorted(known)}"
                )
        rate = ids.lumen_conversion_rate(c.substrate, c.product)
        if rate in seen:
            problems.append(f"{where}: given twice (the id '{rate}')")
        seen.add(rate)
    _raise(problems)


def converted(conversions: Sequence[LumenConversion], luminal: set[str]) -> set[str]:
    """The products of the conversions, which reach the lumen by them.

    A product is the substrate of a further conversion in a chain.

    Raises:
        ValueError: for a conversion whose substrate never reaches the lumen.
    """
    reached = set(luminal)
    changed = True
    while changed:
        changed = False
        for c in conversions:
            if c.substrate in reached and c.product not in reached:
                reached.add(c.product)
                changed = True
    missing = [c for c in conversions if c.substrate not in reached]
    if missing:
        c = missing[0]
        raise ValueError(
            f"LumenConversion {c.substrate} -> {c.product}: '{c.substrate}' does "
            "not reach the lumen (no oral route, biliary excretion, transport of "
            "the intestine or conversion into it), it cannot convert there."
        )
    return reached - luminal


def segment_factor_objects(
    definition: "PBPKDefinition", profiles: Mapping[str, RegionalAbsorption]
) -> list[Sbase]:
    """The factors `f_ka_<seg>_<s>` of the segments of the profiles (not the groups).

    In the order of the substances and the segments.
    """
    return [
        annotate(
            value_parameter(
                ids.regional_absorption_factor(s.sid, seg),
                profiles[s.sid].factors[seg],
                U.dimensionless,
                f"absorption factor {_LABELS[seg]} {s.sid}",
                _SBO_PARAMETER,
                "Factor of the absorption of the segment relative to the "
                f"small intestine, `ka*f_area_{seg}*f_ka_{seg}`: the regional "
                "absorption profile of the substance.",
            ),
            "regional_absorption",
            organ=_region([seg]),
        )
        for s in definition.substances
        if s.sid in profiles
        for seg in ABSORBING
        if seg in profiles[s.sid].segments
    ]


def conversion_objects(conversions: Sequence[LumenConversion]) -> list[Sbase]:
    """The rate constants `kconv_<s>_<m>` of the conversions.

    Annotated with the organ of their segments (none for the stomach or both
    regions).
    """
    return [
        annotate(
            value_parameter(
                ids.lumen_conversion_rate(c.substrate, c.product),
                c.k,
                U.per_hr,
                f"conversion {c.substrate} to {c.product} in the lumen [1/hr]",
                SBO.FORWARD_RATE_CONSTANT,
                "First order rate constant of the conversion in the lumen "
                f"segments {list(c.segments)}, `kconv/60*V*C`.",
            ),
            "lumen_conversion_rate",
            organ=_region(c.segments),
        )
        for c in conversions
    ]
