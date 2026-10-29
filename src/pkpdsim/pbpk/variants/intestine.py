"""The intestine variants: the segmented intestine and the single lumen.

Both replace small and large intestine with one submodel (`GU`, the id of
the small intestine for every intestine variant, `ids.submodel`) with the
boundary of the standard organs on the blood side: the arterial inflow of
both organs, `Qgu = Qsi + Qla`, drains into the portal vein. The output
`Cout_si_<s>` is the plasma leaving the wall; the venous blood of the mucosa
enters the portal vein directly. Their parts:

- the intestinal wall, the standard organ (`StandardBuilder`) of the volume
  `Vgu = Vsi + Vla` with the volume-weighted vascular and interstitial
  fractions, the distribution of the small intestine (its Kp) and the share
  of `ftissue` of its tissue volume;
- the lumen segments: stomach, duodenum, upper and lower jejunum, upper and
  lower ileum, caecum, colon ascendens, transversum, descendens, sigmoid and
  rectum (`segmented`, Thelen 2011), or one lumen, the fluid of the small
  intestine (`single_lumen`); a substance which reaches the lumen (oral
  route, biliary excretion, transport into the lumen) has the dissolved
  concentration `Clumen_<seg>_<s>` per segment, which moves on with the
  transit `V*C/T_<seg>` and from the last segment into the feces
  `Afeces_<s>`; an oral tablet has the undissolved amount `Asolid_<seg>_<s>`
  (in the first segment the dose state `PODOSE_<s>`), which dissolves with
  the dissolution of its formulation (`pkpdsim.pbpk.routes`) and moves with
  the transit too; the bile drains into the duodenum (the single lumen);
  the fed state `FED = 1` slows the gastric emptying (the stomach transit
  time of the fed state);
- the absorption: `ka_<s>` (1/hr) is the absorption rate constant of the
  small intestine, a segment absorbs with `ka_<s> * f_area_<seg>` (the
  mucosal surface per lumen volume relative to the small intestine; the
  stomach does not absorb), a segment of the large intestine optionally
  times `f_ka_colon_<s>` (`SegmentedIntestine.f_ka_colon`, a substance
  with an oral route), a segment of the regional absorption profile of a
  substance times `f_ka_<seg>_<s>` (`RegionalAbsorption`); the single
  lumen absorbs with `ka_<s>`; of the
  absorbed flux the fraction `fabs_<s>` enters the body and the rest the
  feces (the single lumen for every substance, the segments with
  `SegmentedIntestine.fabs` for a substance with an oral route);
- the mucosa (enterocytes) per intestinal segment, `Cmu_<seg>_<s>`, only for
  a substance with chemistry of the intestine (substrate, product,
  transported substance, inhibitor or inducer of a protein of a process of
  the small intestine): flow-limited with the mucosa to plasma coefficient
  `K_mu_<s>` (the Kp of the small intestine, 1 without), perfused by
  arterial blood with its flow `Qmu_<seg> = f_mucosa_flow * Q<region> *
  share` (the mucosa fraction of the blood flow of its region, small or
  large intestine, times the transit share of the segment in the region;
  the wall gets the rest) and draining into the portal vein; the absorbed
  drug enters it; its volume `Vmu_<seg>` is carved out of the tissue of the
  wall of that substance (`V<part>_wall`), so wall and mucosa together have
  the capacity of the wall of a substance without mucosa; the unbound
  concentration `Cu_mu_<seg>_<s>` is an observer. Without such chemistry
  the substance has no mucosa species and is absorbed into the plasma of
  the wall.

The chemistry of the small intestine is placed into the mucosa segments by
the relative expression of its protein (CYP3A4 of the data) times the
mucosa volume, the sum kept: the Vmax of a process stays per litre of small
intestine (`Vmax * Vsi` for the whole intestine, as in the standard organ),
a segment has the weight `w_<protein>_<seg>`; a protein without expression
data (or a process without protein) is placed by mucosa volume over the
segments of the small intestine, the organ of the process. Apical efflux moves the drug from the
mucosa into the lumen of the segment, apical influx back. A protein which a
substance inactivates or induces has one turnover state `E_<protein>` for
the intestine, driven by the unbound mucosa concentration weighted by the
placement weights of the protein.

The bile of every substance drains: of a substance with lumen species into
the duodenum (the single lumen), of any other into the feces.

A conversion in the lumen (`LumenConversion`, the segmented intestine)
converts a dissolved substance into another in chosen segments, first order
with `kconv_<s>_<m>` (1/hr); the product gets lumen species.
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

from sbmlutils.factory import ModelDefinition, Parameter, Sbase
from sbmlutils.metadata import SBO

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.annotation import annotate
from pkpdsim.pbpk.data import value
from pkpdsim.pbpk.organ import Boundary, OrganContext
from pkpdsim.pbpk.organs import Organ
from pkpdsim.pbpk.placement import value_parameter
from pkpdsim.pbpk.routes import ORAL
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import FrozenMap, Q, check_unit, number
from pkpdsim.pbpk.variants.gut import GUT_ORGANS, GutBuilder
from pkpdsim.pbpk.variants.segments import (
    DATA,
    SEGMENTS,
    SMALL,
    segments,
    single_lumen,
)

if TYPE_CHECKING:
    from collections.abc import Iterable

    from pkpdsim.pbpk.definition import PBPKDefinition
    from pkpdsim.pbpk.routes import Route

__all__ = [
    "LumenConversion",
    "PerSubstance",
    "RegionalAbsorption",
    "SegmentedIntestine",
    "SingleLumen",
]

_SBO_PARAMETER = SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER

type PerSubstance = Q | Mapping[str, Q]
"""One value for every oral substance, or a mapping from substance id to value."""


_ABSORBING = tuple(t for t in ids.GUT_SEGMENTS if t != "st")
"""The lumen segments which absorb: every segment but the stomach."""

_LARGE = tuple(t for t in ids.GUT_SEGMENTS if t != "st" and t not in SMALL)
"""The lumen segments of the large intestine, caecum to rectum."""

_LABELS = {t: key.replace("_", " ") for t, key in SEGMENTS.items()}
"""The segments in names by token, `upper jejunum`."""

_REGION = {
    t: Organ.LARGE_INTESTINE if t in _LARGE else Organ.SMALL_INTESTINE
    for t in _ABSORBING
}
"""The organ of an absorbing segment, for the annotation of its factor."""


def _segment_problems(where: str, tokens: "Sequence[str]") -> list[str]:
    """The problems of the segment tokens of an option: none, twice or unknown."""
    if not tokens:
        return [f"{where}: no segments, give the tokens {list(ids.GUT_SEGMENTS)}"]
    problems = []
    twice = sorted({t for t in tokens if list(tokens).count(t) > 1})
    if twice:
        problems.append(f"{where}: the segments {twice} more than once")
    unknown = [t for t in tokens if t not in ids.GUT_SEGMENTS]
    if unknown:
        problems.append(
            f"{where}: unknown segments {unknown}, the segments are "
            f"{list(ids.GUT_SEGMENTS)}"
        )
    return problems


def _raise(problems: list[str]) -> None:
    """Raise the problems of an option, if any.

    Raises:
        ValueError: with every problem.
    """
    if problems:
        raise ValueError("; ".join(problems))


@dataclass(frozen=True)
class RegionalAbsorption:
    """The regional absorption profile of a substance in the segmented intestine.

    A factor of a segment multiplies its absorption `ka_<s>*f_area_<seg>`:
    the permeability of the segment relative to the small intestine, which
    `ka_<s>` describes (the area factor holds the mucosal surface per lumen
    volume, Helander 2014). Every factor is a parameter `f_ka_<seg>_<s>` of
    the top model; a segment missing from the profile has none (1).

    There is no generic default profile: the colonic absorption relative to
    the small intestine depends on the substance, high for a highly
    permeable drug (BCS class I, relative colonic bioavailability above 70
    %) and low for a drug of low permeability (BCS III and IV, below 50 %;
    Tannergren 2009, 42 drugs), so the default is no factor (1). A substance
    whose regional absorption is measured (an intubation or a capsule
    released in a segment, `Byon2018` for apixaban) takes its own profile.

    Attributes:
        factors: the factor of a segment by its token (`du` to `rec`, not
            the stomach, which does not absorb): a number or a `Value`
            (dimensionless), not negative.

    Raises:
        ValueError: for no segment, the stomach, an unknown token, a negative
            factor or a `Value` which is not dimensionless.
    """

    factors: Mapping[str, Q]

    def __post_init__(self) -> None:
        """Store the factors as a hashable read-only mapping and check them."""
        object.__setattr__(self, "factors", FrozenMap(self.factors))
        where = "RegionalAbsorption"
        problems = _segment_problems(where, list(self.factors))
        if "st" in self.factors:
            problems.append(
                f"{where}: the stomach does not absorb, the segments are "
                f"{list(_ABSORBING)}"
            )
        for seg, q in self.factors.items():
            x = number(q)
            if x is not None and x < 0:
                problems.append(
                    f"{where}: the factor of '{seg}' must not be negative, got {x}"
                )
            problems += check_unit(where, f"the factor of '{seg}'", q, U.dimensionless)
        _raise(problems)


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
    the stomach `T` (`T_st`, fasted 22.79 min), and `kconv = 60*f/((1 - f)*T)`;
    a tablet with the first order dissolution `kdis` (1/hr) converts only
    the part which dissolves in the stomach, `g = kd*T/(1 + kd*T)` with
    `kd = kdis/60`, so `f = g*k*T/(1 + k*T)`. Over several segments without
    absorption the drug survives each with `1/(1 + k*T_seg)`.

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
        problems += _segment_problems(where, self.segments)
        x = number(self.k)
        if x is not None and x < 0:
            problems.append(f"{where}: k must not be negative, got {x}")
        problems += check_unit(where, "k", self.k, U.per_hr)
        _raise(problems)


def _oral(routes: "Sequence[Route]") -> frozenset[str]:
    """The substances with an oral route (tablet or solution)."""
    return frozenset(r.substance for r in routes if r.kind in ORAL)


def _freeze(spec: "PerSubstance | None") -> "PerSubstance | None":
    """Store a mapping as a hashable read-only mapping."""
    if isinstance(spec, Mapping):
        return FrozenMap(spec)
    return spec


def _per_substance(spec: PerSubstance, sid: str, default: Q) -> Q:
    """Get the value of a substance: the mapping entry (else `default`) or the scalar."""
    if isinstance(spec, Mapping):
        return spec.get(sid, default)
    return spec


def _check_keys(
    attribute: str,
    spec: object,
    known: frozenset[str],
    allowed: frozenset[str] | None,
) -> None:
    """Check that the keys of a per-substance mapping are substances (of `allowed` if set).

    Raises:
        ValueError: for a key which is no substance of the definition, or, if
            `allowed` is set, a substance which is not in it (no oral route).
    """
    if not isinstance(spec, Mapping):
        return
    for sid in spec:
        if sid not in known:
            raise ValueError(
                f"{attribute}: '{sid}' is no substance of the definition "
                f"{sorted(known)}"
            )
        if allowed is not None and sid not in allowed:
            raise ValueError(
                f"{attribute}: '{sid}' has no oral route, only the substances "
                f"{sorted(allowed)} can have it"
            )


def _only(spec: "PerSubstance | None", oral: frozenset[str]) -> frozenset[str]:
    """The substances which get the parameter: the oral ones, the keys of a mapping among them."""
    if spec is None:
        return frozenset()
    if isinstance(spec, Mapping):
        return frozenset(spec) & oral
    return oral


def _fabs_objects(
    definition: "PBPKDefinition",
    fabs: PerSubstance,
    substances: frozenset[str] | None = None,
    default: Q = 1.0,
) -> list[Sbase]:
    """The fraction absorbed `fabs_<s>` of the substances (all if `None`).

    A substance missing from a mapping has `default`.
    """
    return [
        annotate(
            value_parameter(
                ids.fraction_absorbed(s.sid),
                _per_substance(fabs, s.sid, default),
                U.dimensionless,
                f"fraction absorbed {s.sid}",
                _SBO_PARAMETER,
                "Fraction of the absorbed flux which enters the body, the "
                "rest goes into the feces.",
            ),
            "fraction_absorbed",
        )
        for s in definition.substances
        if substances is None or s.sid in substances
    ]


def _colon_objects(
    definition: "PBPKDefinition", f_ka_colon: PerSubstance, substances: frozenset[str]
) -> list[Sbase]:
    """The factor of the absorption in the large intestine `f_ka_colon_<s>` of the substances."""
    return [
        annotate(
            value_parameter(
                ids.colon_absorption_factor(s.sid),
                _per_substance(f_ka_colon, s.sid, 1.0),
                U.dimensionless,
                f"absorption factor large intestine {s.sid}",
                _SBO_PARAMETER,
                "Factor of the absorption of the segments of the large "
                "intestine, `ka*f_area*f_ka_colon`; a model fits it.",
            ),
            "colon_absorption",
            organ=Organ.LARGE_INTESTINE,
        )
        for s in definition.substances
        if s.sid in substances
    ]


def _ka_value(ka: "PerSubstance | None", sid: str, oral: frozenset[str]) -> Q:
    """The `ka` of a substance: its entry or the scalar, else the default (0 without oral route)."""
    default = value(DATA, "absorption", "ka")
    if isinstance(ka, Mapping) and sid in ka:
        return ka[sid]
    if sid not in oral:
        return 0.0
    return default if ka is None or isinstance(ka, Mapping) else ka


def _ka_objects(definition: "PBPKDefinition", ka: "PerSubstance | None") -> list[Sbase]:
    """The absorption `ka_<s>` of every substance.

    A scalar `ka` is the value of every oral substance; a mapping gives its
    entries (also to a substance without oral route: the reabsorption after
    biliary excretion), the other oral substances the default of the data.
    Every other substance has 0.
    """
    oral = _oral(definition.routes)
    return [
        annotate(
            value_parameter(
                ids.absorption(s.sid),
                _ka_value(ka, s.sid, oral),
                U.per_hr,
                f"absorption {s.sid} [1/hr]",
                _SBO_PARAMETER,
                "Absorption rate constant of the small intestine; a segment absorbs "
                "with the factor of its area, `ka*f_area`; a model fits it."
                if s.sid in oral or (isinstance(ka, Mapping) and s.sid in ka)
                else "No absorption: the substance has no oral route, it is not "
                "reabsorbed after biliary excretion unless a model sets it.",
            ),
            "absorption_rate",
            organ=Organ.SMALL_INTESTINE,
        )
        for s in definition.substances
    ]


@dataclass(frozen=True)
class SegmentedIntestine:
    """The segmented intestine (Thelen 2011), the default intestine.

    Every segment but the stomach absorbs with `ka_<s> * f_area_<seg>`.
    Two optional factors change it, for every substance with an oral route
    (a model fits the parameter of its substance): `fabs` sends the fraction
    `1 - fabs_<s>` of the absorbed flux of every segment into the feces
    (incomplete absorption with a fast onset), `f_ka_colon` scales the
    absorption of the segments of the large intestine (a substance absorbed
    less in the colon). Without them the model has neither parameter, nor
    has a substance without an oral route (its `ka_<s>` is 0): it keeps the
    absorption `ka_<s> * f_area_<seg>` and has no reaction
    `J_feces_<seg>_<s>`. Unlike `SingleLumen`, whose `fabs=None` is
    `fabs_<s> = 1` for every substance, `None` creates no parameter.

    `regional` gives a substance with an oral route the regional absorption
    profile (`RegionalAbsorption`), a factor `f_ka_<seg>_<s>` of the
    absorption of a segment; `conversions` convert a substance into another
    in the lumen (`LumenConversion`). Neither creates an element unless set.

    Attributes:
        name: the name of the component.
        ka: the absorption rate constant of a substance with an oral route,
            the default of the data if `None`; a mapping from substance id
            gives each substance its own (a substance missing from it has the
            default, a substance without oral route may be in it).
        fabs: the fraction of the absorbed flux which enters the body,
            `fabs_<s>` of a substance with an oral route; none (all of it)
            if `None`; a mapping gives it only to the substances in it.
        f_ka_colon: the factor of the absorption of the segments of the large
            intestine (caecum to rectum), `f_ka_colon_<s>` of a substance
            with an oral route; none (1) if `None`; a mapping gives it only
            to the substances in it.
        regional: the regional absorption profile of every substance with an
            oral route, or a mapping from substance id to the profile of the
            substances in it; none (every factor 1) if `None`, the default
            (no generic profile, see `RegionalAbsorption`). A substance with
            `f_ka_colon` may have a profile of the small intestine only.
        conversions: the conversions in the lumen.

    Raises:
        ValueError: (on build) for a key of a mapping which is no substance of
            the definition, or of `fabs`, `f_ka_colon` or `regional` a
            substance without an oral route; for a substance with both
            `f_ka_colon` and a factor of the large intestine; for a
            conversion of a substance which is no substance of the
            definition or never reaches the lumen, or two conversions of the
            same substrate and product.
    """

    name: str = "segmented"
    ka: "PerSubstance | None" = None
    fabs: "PerSubstance | None" = None
    f_ka_colon: "PerSubstance | None" = None
    regional: "RegionalAbsorption | Mapping[str, RegionalAbsorption] | None" = None
    conversions: Sequence[LumenConversion] = ()

    def __post_init__(self) -> None:
        """Store a mapping as a hashable read-only mapping, the conversions as a tuple."""
        for name in ("ka", "fabs", "f_ka_colon"):
            object.__setattr__(self, name, _freeze(getattr(self, name)))
        if isinstance(self.regional, Mapping):
            object.__setattr__(
                self, "regional", FrozenMap[RegionalAbsorption](self.regional)
            )
        object.__setattr__(self, "conversions", tuple(self.conversions))

    @property
    def replaces(self) -> tuple[Organ, ...]:
        """The organs the intestine replaces: small and large intestine."""
        return GUT_ORGANS

    def _oral_sets(
        self, sids: "Iterable[str]", routes: "Sequence[Route]"
    ) -> tuple[frozenset[str], frozenset[str]]:
        """The substances with `fabs_<s>` and with `f_ka_colon_<s>`, checked.

        The one place which body_objects and build use, so that neither can
        drop a key silently.

        Raises:
            ValueError: for a key which is no substance or (`fabs`,
                `f_ka_colon`) has no oral route.
        """
        known, oral = frozenset(sids), _oral(routes)
        _check_keys("ka", self.ka, known, None)
        _check_keys("fabs", self.fabs, known, oral)
        _check_keys("f_ka_colon", self.f_ka_colon, known, oral)
        return _only(self.fabs, oral), _only(self.f_ka_colon, oral)

    def _profiles(
        self, sids: "Iterable[str]", routes: "Sequence[Route]"
    ) -> dict[str, RegionalAbsorption]:
        """The regional absorption profiles by substance, checked.

        Raises:
            ValueError: for a key which is no substance or has no oral route,
                or a substance with `f_ka_colon` and a factor of the large
                intestine.
        """
        known, oral = frozenset(sids), _oral(routes)
        _check_keys("regional", self.regional, known, oral)
        if self.regional is None:
            return {}
        if isinstance(self.regional, Mapping):
            profiles = {sid: p for sid, p in self.regional.items() if sid in oral}
        else:
            profiles = dict.fromkeys(sorted(oral), self.regional)
        colon = _only(self.f_ka_colon, oral)
        twice = sorted(
            sid
            for sid, p in profiles.items()
            if sid in colon and any(seg in _LARGE for seg in p.factors)
        )
        if twice:
            raise ValueError(
                f"regional: {twice} have both f_ka_colon and a regional factor of "
                "the large intestine, which scale the same absorption; give the "
                "large intestine by one of them"
            )
        return profiles

    def _check_conversions(self, sids: "Iterable[str]") -> None:
        """Check the substances of the conversions and that none is given twice.

        Raises:
            ValueError: for a substance which is no substance of the
                definition, or two conversions of the same ids.
        """
        known = frozenset(sids)
        problems = []
        seen: set[str] = set()
        for c in self.conversions:
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

    def _regional_objects(
        self, definition: "PBPKDefinition", profiles: Mapping[str, RegionalAbsorption]
    ) -> list[Sbase]:
        """The factors `f_ka_<seg>_<s>` of the profiles, in the order of the substances and segments."""
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
                organ=_REGION[seg],
            )
            for s in definition.substances
            if s.sid in profiles
            for seg in _ABSORBING
            if seg in profiles[s.sid].factors
        ]

    def _conversion_objects(self) -> list[Sbase]:
        """The rate constants `kconv_<s>_<m>` of the conversions."""
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
                organ=Organ.SMALL_INTESTINE,
            )
            for c in self.conversions
        ]

    def body_objects(self, definition: "PBPKDefinition") -> list[Sbase]:
        """Get `ka_<s>`, `fabs_<s>`, `f_ka_colon_<s>`, the regional factors, the conversion rates, `FED`."""
        sids = [s.sid for s in definition.substances]
        fabs, colon = self._oral_sets(sids, definition.routes)
        profiles = self._profiles(sids, definition.routes)
        self._check_conversions(sids)
        return [
            *_ka_objects(definition, self.ka),
            *(
                _fabs_objects(definition, self.fabs, fabs)
                if self.fabs is not None
                else []
            ),
            *(
                _colon_objects(definition, self.f_ka_colon, colon)
                if self.f_ka_colon is not None
                else []
            ),
            *self._regional_objects(definition, profiles),
            *self._conversion_objects(),
            annotate(
                Parameter(
                    ids.FED,
                    0.0,
                    U.dimensionless,
                    constant=True,
                    name="fed state (0 fasted, 1 fed)",
                    sboTerm=_SBO_PARAMETER,
                    notes="The prandial state: 1 slows the gastric emptying to the "
                    "stomach transit time of the fed state (an initial value, an "
                    "experiment sets it before the simulation).",
                ),
                "fed",
            ),
        ]

    def build(self, ctx: OrganContext) -> tuple[ModelDefinition, Boundary]:
        """Build the model of the intestine and its boundary.

        Raises:
            ValueError: for chemistry the intestine cannot place: a process
                of the large intestine, a custom reaction, a basolateral
                transport of the mucosa.
        """
        sids = [s.sid for s in ctx.substances]
        fabs, colon = self._oral_sets(sids, ctx.routes)
        profiles = self._profiles(sids, ctx.routes)
        self._check_conversions(sids)
        return GutBuilder(
            ctx,
            self.name,
            segments(),
            single=False,
            fabs=fabs,
            colon=colon,
            regional={sid: frozenset(p.factors) for sid, p in profiles.items()},
            conversions=tuple(self.conversions),
        ).build()


@dataclass(frozen=True)
class SingleLumen:
    """The single intestinal lumen with first order absorption.

    The first order absorption of most current models: one lumen (the fluid
    of the small intestine) with the transit into the feces, absorption with
    `ka_<s>`, of which the fraction `fabs_<s>` enters the body. The lumen
    keeps a physiological transit `t_lumen` into the feces, so the fraction
    of a solution which enters the body is `fabs*ka/(ka + 1/t_lumen)`: the
    fitted absorption and fraction absorbed of the hand-written models,
    which had no transit, do not port unchanged (a refit is expected).

    `fabs_<s>` is a parameter of every substance, 1 if `fabs=None` (unlike
    `SegmentedIntestine`, which creates it only when set and only for the
    substances with an oral route).

    Attributes:
        name: the name of the component.
        ka: the absorption rate constant of a substance with an oral route,
            the default of the data if `None`; a mapping from substance id
            gives each substance its own (the others have the default).
        fabs: the fraction of the absorbed flux which enters the body, 1 if
            `None`; in a mapping a substance missing from it has 1.

    Raises:
        ValueError: (on build) for a key of a mapping which is no substance of
            the definition.
    """

    name: str = "single_lumen"
    ka: "PerSubstance | None" = None
    fabs: "PerSubstance | None" = None

    def __post_init__(self) -> None:
        """Store a mapping as a hashable read-only mapping."""
        for name in ("ka", "fabs"):
            object.__setattr__(self, name, _freeze(getattr(self, name)))

    @property
    def replaces(self) -> tuple[Organ, ...]:
        """The organs the intestine replaces: small and large intestine."""
        return GUT_ORGANS

    def _check(self, sids: "Iterable[str]") -> None:
        """Check the keys of the mappings (any substance has `fabs_<s>`).

        Raises:
            ValueError: for a key which is no substance.
        """
        known = frozenset(sids)
        _check_keys("ka", self.ka, known, None)
        _check_keys("fabs", self.fabs, known, None)

    def body_objects(self, definition: "PBPKDefinition") -> list[Sbase]:
        """Get `ka_<s>` and `fabs_<s>` of every substance."""
        self._check(s.sid for s in definition.substances)
        return [
            *_ka_objects(definition, self.ka),
            *_fabs_objects(
                definition,
                self.fabs if self.fabs is not None else 1.0,
                default=1.0,
            ),
        ]

    def build(self, ctx: OrganContext) -> tuple[ModelDefinition, Boundary]:
        """Build the model of the intestine and its boundary.

        Raises:
            ValueError: for chemistry the intestine cannot place, or a key of
                a mapping which is no substance.
        """
        self._check(s.sid for s in ctx.substances)
        return GutBuilder(
            ctx,
            self.name,
            single_lumen(),
            single=True,
            fabs=frozenset(s.sid for s in ctx.substances),
        ).build()
