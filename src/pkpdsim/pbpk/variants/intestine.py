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
  times `f_ka_colon_<s>` (`SegmentedIntestine.f_ka_colon`); the single
  lumen absorbs with `ka_<s>`; of the absorbed flux the fraction `fabs_<s>`
  enters the body and the rest the feces (the single lumen, the segments
  with `SegmentedIntestine.fabs`);
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
"""

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
from pkpdsim.pbpk.values import Q
from pkpdsim.pbpk.variants.gut import GUT_ORGANS, GutBuilder
from pkpdsim.pbpk.variants.segments import DATA, segments, single_lumen

if TYPE_CHECKING:
    from pkpdsim.pbpk.definition import PBPKDefinition

__all__ = ["SegmentedIntestine", "SingleLumen"]

_SBO_PARAMETER = SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER


def _fabs_objects(definition: "PBPKDefinition", fabs: Q) -> list[Sbase]:
    """The fraction absorbed `fabs_<s>` of every substance."""
    return [
        annotate(
            value_parameter(
                ids.fraction_absorbed(s.sid),
                fabs,
                U.dimensionless,
                f"fraction absorbed {s.sid}",
                _SBO_PARAMETER,
                "Fraction of the absorbed flux which enters the body, the "
                "rest goes into the feces.",
            ),
            "fraction_absorbed",
        )
        for s in definition.substances
    ]


def _colon_objects(definition: "PBPKDefinition", f_ka_colon: Q) -> list[Sbase]:
    """The factor of the absorption in the large intestine `f_ka_colon_<s>` of every substance."""
    return [
        annotate(
            value_parameter(
                ids.colon_absorption_factor(s.sid),
                f_ka_colon,
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
    ]


def _ka_objects(definition: "PBPKDefinition", ka: Q | None) -> list[Sbase]:
    """The absorption `ka_<s>` of every substance: of the oral ones the default, else 0."""
    oral = {r.substance for r in definition.routes if r.kind in ORAL}
    default = ka if ka is not None else value(DATA, "absorption", "ka")
    return [
        annotate(
            value_parameter(
                ids.absorption(s.sid),
                default if s.sid in oral else 0.0,
                U.per_hr,
                f"absorption {s.sid} [1/hr]",
                _SBO_PARAMETER,
                "Absorption rate constant of the small intestine; a segment absorbs "
                "with the factor of its area, `ka*f_area`; a model fits it."
                if s.sid in oral
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
    Two optional factors change it, for every substance of the definition
    (a model fits the parameter of its substance): `fabs` sends the fraction
    `1 - fabs_<s>` of the absorbed flux of every segment into the feces
    (incomplete absorption with a fast onset), `f_ka_colon` scales the
    absorption of the segments of the large intestine (a substance absorbed
    less in the colon). Without them the model has neither parameter.

    Attributes:
        name: the name of the component.
        ka: the absorption rate constant of a substance with an oral route,
            the default of the data if `None`.
        fabs: the fraction of the absorbed flux which enters the body,
            `fabs_<s>`; none (all of it) if `None`.
        f_ka_colon: the factor of the absorption of the segments of the large
            intestine (caecum to rectum), `f_ka_colon_<s>`; none (1) if
            `None`.
    """

    name: str = "segmented"
    ka: Q | None = None
    fabs: Q | None = None
    f_ka_colon: Q | None = None

    @property
    def replaces(self) -> tuple[Organ, ...]:
        """The organs the intestine replaces: small and large intestine."""
        return GUT_ORGANS

    def body_objects(self, definition: "PBPKDefinition") -> list[Sbase]:
        """Get `ka_<s>` (`fabs_<s>`, `f_ka_colon_<s>`) of every substance and `FED`."""
        return [
            *_ka_objects(definition, self.ka),
            *(_fabs_objects(definition, self.fabs) if self.fabs is not None else []),
            *(
                _colon_objects(definition, self.f_ka_colon)
                if self.f_ka_colon is not None
                else []
            ),
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
        return GutBuilder(
            ctx,
            self.name,
            segments(),
            single=False,
            fabs=self.fabs is not None,
            colon=self.f_ka_colon is not None,
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

    Attributes:
        name: the name of the component.
        ka: the absorption rate constant of a substance with an oral route,
            the default of the data if `None`.
        fabs: the fraction of the absorbed flux which enters the body, 1 if
            `None`.
    """

    name: str = "single_lumen"
    ka: Q | None = None
    fabs: Q | None = None

    @property
    def replaces(self) -> tuple[Organ, ...]:
        """The organs the intestine replaces: small and large intestine."""
        return GUT_ORGANS

    def body_objects(self, definition: "PBPKDefinition") -> list[Sbase]:
        """Get `ka_<s>` and `fabs_<s>` of every substance."""
        return [
            *_ka_objects(definition, self.ka),
            *_fabs_objects(definition, self.fabs if self.fabs is not None else 1.0),
        ]

    def build(self, ctx: OrganContext) -> tuple[ModelDefinition, Boundary]:
        """Build the model of the intestine and its boundary.

        Raises:
            ValueError: for chemistry the intestine cannot place.
        """
        return GutBuilder(
            ctx, self.name, single_lumen(), single=True, fabs=True
        ).build()
