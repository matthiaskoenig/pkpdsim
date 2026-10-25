"""The intestine variants: the segmented intestine and the single lumen.

Both replace small and large intestine with one submodel (`GU`) with the
boundary of the standard organs on the blood side: the arterial inflow of
both organs, `Qgu = Qsi + Qla`, drains into the portal vein. Their parts:

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
  stomach does not absorb); the single lumen absorbs with `ka_<s>`, of which
  the fraction `fabs_<s>` enters the body and the rest the feces;
- the mucosa (enterocytes) per intestinal segment, `Cmu_<seg>_<s>`, only for
  a substance with chemistry of the intestine (substrate, product,
  transported substance, inhibitor or inducer of a protein of a process of
  the small intestine): flow-limited with the mucosa to plasma coefficient
  `K_mu_<s>` (the Kp of the small intestine, 1 without), perfused by
  arterial blood with its share `Qmu_<seg>` of the intestinal blood flow
  (the wall gets the rest) and draining into the portal vein; the absorbed
  drug enters it. Without such chemistry the substance has no mucosa
  species and is absorbed into the plasma of the wall.

The chemistry of the small intestine is placed into the mucosa segments by
the relative expression of its protein (CYP3A4, ABCB1, ABCG2 of the data)
times the mucosa volume, the sum kept: the Vmax of a process stays per
litre of small intestine (`Vmax * Vsi` for the whole intestine, as in the
standard organ), a segment has the weight `w_<protein>_<seg>`; a protein
without expression data is placed by mucosa volume. Apical efflux moves the
drug from the mucosa into the lumen of the segment, apical influx back.
"""

from dataclasses import dataclass
from typing import TYPE_CHECKING

from sbmlutils.factory import ModelDefinition, Parameter, Sbase
from sbmlutils.metadata import SBO

from pkpdsim.pbpk import ids
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


def _ka_objects(definition: "PBPKDefinition", ka: Q | None) -> list[Sbase]:
    """The absorption `ka_<s>` of every substance: of the oral ones the default, else 0."""
    oral = {r.substance for r in definition.routes if r.kind in ORAL}
    default = ka if ka is not None else value(DATA, "absorption", "ka")
    return [
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
        )
        for s in definition.substances
    ]


@dataclass(frozen=True)
class SegmentedIntestine:
    """The segmented intestine (Thelen 2011), the default intestine.

    Attributes:
        name: the name of the component.
        ka: the absorption rate constant of a substance with an oral route,
            the default of the data if `None`.
    """

    name: str = "segmented"
    ka: Q | None = None

    @property
    def replaces(self) -> tuple[Organ, ...]:
        """The organs the intestine replaces: small and large intestine."""
        return GUT_ORGANS

    def body_objects(self, definition: "PBPKDefinition") -> list[Sbase]:
        """Get `ka_<s>` of every substance and the prandial state `FED`."""
        return [
            *_ka_objects(definition, self.ka),
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
        ]

    def build(self, ctx: OrganContext) -> tuple[ModelDefinition, Boundary]:
        """Build the model of the intestine and its boundary.

        Raises:
            ValueError: for chemistry the intestine cannot place: a process
                of the large intestine, a custom reaction, a basolateral
                transport of the mucosa.
        """
        return GutBuilder(ctx, self.name, segments(), single=False).build()


@dataclass(frozen=True)
class SingleLumen:
    """The single intestinal lumen with first order absorption.

    The first order absorption of most current models: one lumen (the fluid
    of the small intestine) with the transit into the feces, absorption with
    `ka_<s>`, of which the fraction `fabs_<s>` enters the body.

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
            *(
                value_parameter(
                    ids.fraction_absorbed(s.sid),
                    self.fabs if self.fabs is not None else 1.0,
                    U.dimensionless,
                    f"fraction absorbed {s.sid}",
                    _SBO_PARAMETER,
                    "Fraction of the absorbed flux which enters the body, the "
                    "rest goes into the feces.",
                )
                for s in definition.substances
            ),
        ]

    def build(self, ctx: OrganContext) -> tuple[ModelDefinition, Boundary]:
        """Build the model of the intestine and its boundary.

        Raises:
            ValueError: for chemistry the intestine cannot place.
        """
        return GutBuilder(ctx, self.name, single_lumen(), single=True).build()
