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
  the transit too; an enteric-coated tablet does not dissolve in the
  stomach, it empties undissolved and dissolves from the duodenum on; the
  bile drains into the duodenum (the single lumen);
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
from typing import TYPE_CHECKING, NamedTuple

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
from pkpdsim.pbpk.values import FrozenMap, Q, Value, check_positive, check_unit
from pkpdsim.pbpk.variants.gut import GUT_ORGANS, GutBuilder, mucosal_substances
from pkpdsim.pbpk.variants.lumen_options import (
    COLON,
    LARGE,
    LumenConversion,
    RegionalAbsorption,
    check_conversions,
    conversion_objects,
    segment_factor_objects,
)
from pkpdsim.pbpk.variants.segments import (
    intestine_table,
    segments,
    single_lumen,
)

if TYPE_CHECKING:
    from pkpdsim.pbpk.chemistry import CustomReaction, Metabolism, Transport
    from pkpdsim.pbpk.definition import PBPKDefinition
    from pkpdsim.pbpk.routes import Route
    from pkpdsim.pbpk.substance import Substance

__all__ = ["PerSubstance", "SegmentedIntestine", "SingleLumen", "absorbed_substances"]

_SBO_PARAMETER = SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER

type PerSubstance = Q | Mapping[str, Q]
"""One value for every oral substance, or a mapping from substance id to value."""


def _oral(routes: Sequence[Route]) -> frozenset[str]:
    """The substances with an oral route (tablet or solution)."""
    return frozenset(r.substance for r in routes if r.kind in ORAL)


def _freeze(spec: PerSubstance | None) -> PerSubstance | None:
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
    spec: PerSubstance | RegionalAbsorption | Mapping[str, RegionalAbsorption] | None,
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


def _mucosa_map(spec: Mapping[str, Q] | None, per_volume: bool) -> FrozenMap[Q] | None:
    """Store `mucosa_ps` as a hashable read-only mapping, its values checked.

    The values are checked as the `ps` of an organ is
    (`PBPKDefinition.validate`): positive, and a `Value` in the unit of the
    form (1/min per litre of mucosa, l/min of the mucosa of the small
    intestine). A negative exchange gives a negative hybrid flow, so the flows
    of the mucosa run backwards; zero cuts the enterocytes off from the blood.
    A `Formula` passes both checks: its value is only known when the model is
    built, as for the `ps` of an organ, so a formula which can turn negative
    is the model's responsibility.

    Args:
        spec: the exchange of the mucosa by substance id, or `None`.
        per_volume: the values are per litre of mucosa
            (`mucosa_ps_per_volume`).

    Returns:
        The read-only copy, `None` for `None`.

    Raises:
        ValueError: for a form which is no bool or is given without
            `mucosa_ps`, for anything but a mapping (the permeability-limited
            mucosa is opt-in per substance, there is no value for all of
            them), and for a value which is not positive or is a `Value` in
            another unit.
    """
    if not isinstance(per_volume, bool):
        raise ValueError(
            f"mucosa_ps_per_volume: must be a bool, got {per_volume!r}: the exchange "
            "is given per litre of mucosa (True) or of the mucosa of the small "
            "intestine (False)"
        )
    if spec is None:
        if not per_volume:
            raise ValueError(
                "mucosa_ps_per_volume: without mucosa_ps it changes nothing, the "
                "mucosa of every substance is flow-limited; name the substances "
                "with the exchange of the mucosa of their small intestine in "
                "mucosa_ps"
            )
        return None
    if not isinstance(spec, Mapping):
        raise ValueError(
            f"mucosa_ps: a mapping from substance id to its exchange, not {spec!r}: "
            "the permeability-limited mucosa is given per substance"
        )
    unit = U.per_min if per_volume else U.l_per_min
    problems: list[str] = []
    for sid, q in spec.items():
        problems += check_positive("mucosa_ps", f"the exchange of '{sid}'", q)
        problems += check_unit("mucosa_ps", f"the exchange of '{sid}'", q, unit)
    if problems:
        raise ValueError("; ".join(problems))
    return FrozenMap(spec)


def _check_mucosa_ps(
    spec: Mapping[str, Q] | None,
    substances: Sequence[Substance],
    chemistry: Sequence[Metabolism | Transport | CustomReaction],
) -> dict[str, Q]:
    """Check `mucosa_ps` against the substances with a mucosa.

    Args:
        spec: the exchange of the mucosa by substance id, `None` for none.
        substances: the substances of the definition.
        chemistry: the processes of the definition.

    Returns:
        The exchange by substance id, empty without the option.

    Raises:
        ValueError: for a key which is no substance of the definition or has
            no mucosa species (no chemistry of the small intestine).
    """
    if spec is None:
        return {}
    known = frozenset(s.sid for s in substances)
    _check_keys("mucosa_ps", spec, known, None)
    mucosal = mucosal_substances(substances, chemistry)
    for sid in spec:
        if sid not in mucosal:
            raise ValueError(
                f"mucosa_ps: '{sid}' has no mucosa: only a substance with "
                "chemistry of the small intestine (a substrate, a product, a "
                "transported substance, or an inhibitor, inducer or stabiliser "
                "of a protein of such a process) has mucosa species; the "
                f"substances with a mucosa are {sorted(mucosal)}"
            )
    return dict(spec)


def _only(spec: PerSubstance | None, oral: frozenset[str]) -> frozenset[str]:
    """The substances which get the parameter: the oral ones, the keys of a mapping among them."""
    if spec is None:
        return frozenset()
    if isinstance(spec, Mapping):
        return frozenset(spec) & oral
    return oral


def _fabs_objects(
    definition: PBPKDefinition,
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


def _colon_text(q: Q) -> str:
    """The description of `f_ka_colon_<s>`: a sourced value names its source, else a fit.

    A `Value` with a source (not `fit`) is a measured regional absorption
    (`RegionalAbsorption({"colon": Value(..., source="Byon2018")})`), which
    no fit sets; the notes of the value follow with the citation.
    """
    text = (
        "Factor of the absorption of the segments of the large intestine, "
        "`ka*f_area*f_ka_colon`"
    )
    if isinstance(q, Value) and q.source not in (None, "fit"):
        return (
            f"{text}: the regional absorption profile of the substance, the value "
            f"of its source {q.source}, not fitted."
        )
    return f"{text}; a model fits it."


def _colon_objects(definition: PBPKDefinition, values: Mapping[str, Q]) -> list[Sbase]:
    """The factor of the absorption in the large intestine `f_ka_colon_<s>` of the substances.

    The values by substance: of `f_ka_colon` or of the group `colon` of a
    regional absorption profile, one parameter either way; the notes name the
    source of a sourced value (`_colon_text`).
    """
    return [
        annotate(
            value_parameter(
                ids.colon_absorption_factor(s.sid),
                values[s.sid],
                U.dimensionless,
                f"absorption factor large intestine {s.sid}",
                _SBO_PARAMETER,
                _colon_text(values[s.sid]),
            ),
            "colon_absorption",
            organ=Organ.LARGE_INTESTINE,
        )
        for s in definition.substances
        if s.sid in values
    ]


def _ka_value(ka: PerSubstance | None, sid: str, oral: frozenset[str], table: str) -> Q:
    """The `ka` of a substance: its entry or the scalar, else the default (0 without oral route).

    The default is the one of the data file of the intestine of the species.
    """
    default = value(table, "absorption", "ka")
    if isinstance(ka, Mapping) and sid in ka:
        return ka[sid]
    if sid not in oral:
        return 0.0
    return default if ka is None or isinstance(ka, Mapping) else ka


def absorbed_substances(
    definition: PBPKDefinition, ka: PerSubstance | None
) -> list[str]:
    """Get the substances which the intestine absorbs, `ka_<s>` of the definition.

    Every substance with an oral route, and a substance which a mapping
    `ka` names (also without oral route: the reabsorption after biliary
    excretion); every other substance has `ka_<s> = 0`.

    Args:
        definition: the definition.
        ka: the `ka` of the intestine variant.

    Returns:
        The ids, in the order of the substances.
    """
    oral = _oral(definition.routes)
    named = frozenset(ka) if isinstance(ka, Mapping) else frozenset()
    return [s.sid for s in definition.substances if s.sid in oral or s.sid in named]


def _ka_objects(definition: PBPKDefinition, ka: PerSubstance | None) -> list[Sbase]:
    """The absorption `ka_<s>` of every substance.

    A scalar `ka` is the value of every oral substance; a mapping gives its
    entries (also to a substance without oral route: the reabsorption after
    biliary excretion), the other oral substances the default of the data.
    Every other substance has 0.
    """
    oral = _oral(definition.routes)
    absorbed = absorbed_substances(definition, ka)
    return [
        annotate(
            value_parameter(
                ids.absorption(s.sid),
                _ka_value(ka, s.sid, oral, intestine_table(definition.species)),
                U.per_hr,
                f"absorption {s.sid} [1/hr]",
                _SBO_PARAMETER,
                "Absorption rate constant of the small intestine; a segment absorbs "
                "with the factor of its area, `ka*f_area`; a model fits it."
                if s.sid in absorbed
                else "No absorption: the substance has no oral route, it is not "
                "reabsorbed after biliary excretion unless a model sets it.",
            ),
            "absorption_rate",
            organ=Organ.SMALL_INTESTINE,
        )
        for s in definition.substances
    ]


class _Absorption(NamedTuple):
    """The factors of the absorption of the segmented intestine by substance.

    Attributes:
        fabs: the substances with `fabs_<s>`.
        colon: the value of `f_ka_colon_<s>` by substance (of `f_ka_colon`
            or of the group `colon` of a profile).
        profiles: the regional absorption profiles by substance.
        mucosa_ps: the exchange of the permeability-limited mucosa by
            substance (not a factor of the absorption, checked in the same
            place).
    """

    fabs: frozenset[str]
    colon: dict[str, Q]
    profiles: dict[str, RegionalAbsorption]
    mucosa_ps: dict[str, Q]


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
            to the substances in it. The shorthand of the group `colon` of a
            regional profile, `RegionalAbsorption({"colon": x})`, which makes
            the same parameter.
        regional: the regional absorption profile of every substance with an
            oral route, or a mapping from substance id to the profile of the
            substances in it; none (every factor 1) if `None`, the default
            (no generic profile, see `RegionalAbsorption`). A substance with
            `f_ka_colon` may have a profile of the small intestine only.
        conversions: the conversions in the lumen.
        mucosa_ps: the exchange between the mucosal blood and the enterocytes
            of the substances in it (`PS_mu_<seg>_<s>`), which makes their
            mucosa permeability-limited instead of flow-limited (the Qgut
            model, see the intestine page); the mucosa of every other
            substance is flow-limited, and `None` (the default) creates
            nothing. A substance without mucosa species is refused.
        mucosa_ps_per_volume: the values of `mucosa_ps` are per litre of
            mucosa (`kps_mu_<s>`, 1/min), so the exchange of a segment follows
            the mucosa of the individual; `False` takes them as the exchange
            of the mucosa of the small intestine (`PS_mu_ref_<s>`, l/min, the
            `CLperm` of the Qgut model, which the literature measures for that
            organ), of which a segment has its share by mucosa volume relative
            to the mucosa of the small intestine: the segments of the small
            intestine sum to the given value and the colon has its own share
            at the same exchange per litre of mucosa. Both forms are the one
            quantity `PS_mu_ref_<s>/Vmu_si`.

    Raises:
        ValueError: (on build) for a key of a mapping which is no substance of
            the definition, or of `fabs`, `f_ka_colon` or `regional` a
            substance without an oral route; for a substance with both
            `f_ka_colon` and a factor of the large intestine; for a
            conversion of a substance which is no substance of the
            definition or never reaches the lumen, or two conversions of the
            same substrate and product; for a `mucosa_ps` of a substance
            without mucosa species.
    """

    name: str = "segmented"
    ka: PerSubstance | None = None
    fabs: PerSubstance | None = None
    f_ka_colon: PerSubstance | None = None
    regional: RegionalAbsorption | Mapping[str, RegionalAbsorption] | None = None
    conversions: Sequence[LumenConversion] = ()
    mucosa_ps: Mapping[str, Q] | None = None
    mucosa_ps_per_volume: bool = True

    def __post_init__(self) -> None:
        """Store a mapping as a hashable read-only mapping, the conversions as a tuple."""
        for name in ("ka", "fabs", "f_ka_colon"):
            object.__setattr__(self, name, _freeze(getattr(self, name)))
        if isinstance(self.regional, Mapping):
            object.__setattr__(
                self, "regional", FrozenMap[RegionalAbsorption](self.regional)
            )
        object.__setattr__(self, "conversions", tuple(self.conversions))
        object.__setattr__(
            self, "mucosa_ps", _mucosa_map(self.mucosa_ps, self.mucosa_ps_per_volume)
        )

    @property
    def replaces(self) -> tuple[Organ, ...]:
        """The organs the intestine replaces: small and large intestine."""
        return GUT_ORGANS

    def _resolve(
        self,
        substances: Sequence[Substance],
        routes: Sequence[Route],
        chemistry: Sequence[Metabolism | Transport | CustomReaction],
    ) -> _Absorption:
        """The factors of the absorption by substance, and the mucosa, checked.

        The one place which body_objects and build use, so that neither can
        drop a key silently. `f_ka_colon` is the shorthand of the group
        `colon` of a profile: both give `f_ka_colon_<s>`.

        Raises:
            ValueError: for a key which is no substance or (`fabs`,
                `f_ka_colon`, `regional`) has no oral route, a substance with
                `f_ka_colon` and a factor of the large intestine in its
                profile, a conversion of an unknown substance or twice, or a
                `mucosa_ps` of a substance without mucosa species.
        """
        sids = [s.sid for s in substances]
        known, oral = frozenset(sids), _oral(routes)
        _check_keys("ka", self.ka, known, None)
        _check_keys("fabs", self.fabs, known, oral)
        _check_keys("f_ka_colon", self.f_ka_colon, known, oral)
        _check_keys("regional", self.regional, known, oral)
        check_conversions(self.conversions, sids)
        if isinstance(self.regional, Mapping):
            profiles = {sid: p for sid, p in self.regional.items() if sid in oral}
        elif self.regional is not None:
            profiles = dict.fromkeys(sorted(oral), self.regional)
        else:
            profiles = {}
        colon: dict[str, Q] = {}
        if self.f_ka_colon is not None:
            colon = {
                sid: _per_substance(self.f_ka_colon, sid, 1.0)
                for sid in _only(self.f_ka_colon, oral)
            }
        twice = sorted(
            sid for sid, p in profiles.items() if sid in colon and p.covers(LARGE)
        )
        if twice:
            raise ValueError(
                f"regional: {twice} have both f_ka_colon and a regional factor of "
                "the large intestine, which scale the same absorption; give the "
                "large intestine by one of them (f_ka_colon is the shorthand of "
                "the group 'colon' of the profile)"
            )
        for sid, p in profiles.items():
            if COLON in p.factors:
                colon[sid] = p.factors[COLON]
        return _Absorption(
            _only(self.fabs, oral),
            colon,
            profiles,
            _check_mucosa_ps(self.mucosa_ps, substances, chemistry),
        )

    def absorption(self, definition: PBPKDefinition) -> dict[str, tuple[str, ...]]:
        """Get the absorbed substances with the options of their absorption.

        The substances of `absorbed_substances` with what the model gives
        each of them (`_resolve`, as the build): `fraction absorbed`
        (`fabs_<s>`), `colon factor` (`f_ka_colon_<s>`) and the regional
        factors of segments (`f_ka_<seg>_<s>`), `factor of ile2` for one
        segment, `factors of ile2, asc` for several.

        Args:
            definition: the definition.

        Returns:
            The options by substance id, in the order of the substances.
        """
        a = self._resolve(
            definition.substances, definition.routes, definition.chemistry
        )
        result = {}
        for sid in absorbed_substances(definition, self.ka):
            options = []
            if sid in a.fabs:
                options.append("fraction absorbed")
            if sid in a.colon:
                options.append("colon factor")
            if sid in a.profiles and a.profiles[sid].segments:
                own = [t for t in ids.GUT_SEGMENTS if t in a.profiles[sid].segments]
                noun = "factors" if len(own) > 1 else "factor"
                options.append(f"{noun} of " + ", ".join(own))
            result[sid] = tuple(options)
        return result

    def body_objects(self, definition: PBPKDefinition) -> list[Sbase]:
        """Get `ka_<s>`, `fabs_<s>`, `f_ka_colon_<s>`, the regional factors, the conversion rates, `FED`."""
        a = self._resolve(
            definition.substances, definition.routes, definition.chemistry
        )
        return [
            *_ka_objects(definition, self.ka),
            *(
                _fabs_objects(definition, self.fabs, a.fabs)
                if self.fabs is not None
                else []
            ),
            *_colon_objects(definition, a.colon),
            *segment_factor_objects(definition, a.profiles),
            *conversion_objects(self.conversions),
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
                transport of the mucosa; a conversion of a substance which
                never reaches the lumen.
        """
        a = self._resolve(ctx.substances, ctx.routes, ctx.chemistry)
        return GutBuilder(
            ctx,
            self.name,
            segments(intestine_table(ctx.species)),
            single=False,
            fabs=a.fabs,
            colon=frozenset(a.colon),
            regional={sid: p.segments for sid, p in a.profiles.items()},
            conversions=tuple(self.conversions),
            mucosa_ps=a.mucosa_ps,
            mucosa_ps_per_volume=self.mucosa_ps_per_volume,
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

    The single lumen has no segments: a regional absorption profile and the
    conversions in the lumen need the segmented intestine, `regional` and
    `conversions` are refused with that message; it has no stomach and no
    fed state, so an enteric-coated tablet and a factor of the dissolution
    in the fed state are refused on build, naming the substance.

    Attributes:
        name: the name of the component.
        ka: the absorption rate constant of a substance with an oral route,
            the default of the data if `None`; a mapping from substance id
            gives each substance its own (the others have the default).
        fabs: the fraction of the absorbed flux which enters the body, 1 if
            `None`; in a mapping a substance missing from it has 1.
        regional: refused if set (no segments).
        conversions: refused if set (no segments).
        mucosa_ps: the exchange between the mucosal blood and the enterocytes
            of the substances in it, which makes their mucosa
            permeability-limited (as of `SegmentedIntestine`; the one mucosa
            of the single lumen is that of the small intestine).
        mucosa_ps_per_volume: the values of `mucosa_ps` are per litre of
            mucosa (1/min), else the exchange of the mucosa of the small
            intestine (l/min), which the one mucosa of the single lumen is.

    Raises:
        ValueError: for `regional` or `conversions`; (on build) for a key of
            a mapping which is no substance of the definition, an
            enteric-coated tablet or a fed factor of a dissolution, or a
            `mucosa_ps` of a substance without mucosa species.
    """

    name: str = "single_lumen"
    ka: PerSubstance | None = None
    fabs: PerSubstance | None = None
    regional: RegionalAbsorption | Mapping[str, RegionalAbsorption] | None = None
    conversions: Sequence[LumenConversion] = ()
    mucosa_ps: Mapping[str, Q] | None = None
    mucosa_ps_per_volume: bool = True

    def __post_init__(self) -> None:
        """Store a mapping as a hashable read-only mapping, refuse the options of the segments.

        Raises:
            ValueError: for `regional` or `conversions`.
        """
        for name in ("ka", "fabs"):
            object.__setattr__(self, name, _freeze(getattr(self, name)))
        given = [n for n in ("regional", "conversions") if getattr(self, n)]
        if given:
            raise ValueError(
                f"SingleLumen: {given} need the lumen segments of the segmented "
                "intestine (SegmentedIntestine), the single lumen has one lumen"
            )
        object.__setattr__(self, "conversions", ())
        object.__setattr__(
            self, "mucosa_ps", _mucosa_map(self.mucosa_ps, self.mucosa_ps_per_volume)
        )

    @property
    def replaces(self) -> tuple[Organ, ...]:
        """The organs the intestine replaces: small and large intestine."""
        return GUT_ORGANS

    def _check(
        self,
        substances: Sequence[Substance],
        chemistry: Sequence[Metabolism | Transport | CustomReaction],
    ) -> dict[str, Q]:
        """Check the keys of the mappings (any substance has `fabs_<s>`).

        Returns:
            The exchange of the permeability-limited mucosa by substance.

        Raises:
            ValueError: for a key which is no substance, or a `mucosa_ps` of
                a substance without mucosa species.
        """
        known = frozenset(s.sid for s in substances)
        _check_keys("ka", self.ka, known, None)
        _check_keys("fabs", self.fabs, known, None)
        return _check_mucosa_ps(self.mucosa_ps, substances, chemistry)

    def absorption(self, definition: PBPKDefinition) -> dict[str, tuple[str, ...]]:
        """Get the absorbed substances with the options of their absorption.

        The substances of `absorbed_substances`, `fraction absorbed` for a
        substance whose `fabs_<s>` the model gives (every substance has one,
        1 unless `fabs` sets it).

        Args:
            definition: the definition.

        Returns:
            The options by substance id, in the order of the substances.
        """
        self._check(definition.substances, definition.chemistry)
        fabs = self.fabs
        return {
            sid: (
                ("fraction absorbed",)
                if fabs is not None and (not isinstance(fabs, Mapping) or sid in fabs)
                else ()
            )
            for sid in absorbed_substances(definition, self.ka)
        }

    def body_objects(self, definition: PBPKDefinition) -> list[Sbase]:
        """Get `ka_<s>` and `fabs_<s>` of every substance."""
        self._check(definition.substances, definition.chemistry)
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
            ValueError: for chemistry the intestine cannot place, a key of
                a mapping which is no substance, an enteric-coated tablet or
                a fed factor of a dissolution.
        """
        mucosa_ps = self._check(ctx.substances, ctx.chemistry)
        return GutBuilder(
            ctx,
            self.name,
            single_lumen(intestine_table(ctx.species)),
            single=True,
            fabs=frozenset(s.sid for s in ctx.substances),
            mucosa_ps=mucosa_ps,
            mucosa_ps_per_volume=self.mucosa_ps_per_volume,
        ).build()
