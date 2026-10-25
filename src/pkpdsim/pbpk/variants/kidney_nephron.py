"""The nephron kidney: glomerular filtration into the tubule and its segments.

`NephronKidney` builds the kidney with the tubule of the nephron. The unbound
plasma concentration is filtered into the lumen of the proximal tubule,
`fu*GFR*Cpl` (of the plasma leaving the kidney, `Cor/K_or` for the
flow-limited kidney), and the tubular fluid flows through the lumens of the
segments in series, proximal tubule `pt`, loop of Henle `lh`, distal tubule
`dt` and collecting duct `cd`; the collecting duct drains into the urine
`Aurine_<s>`. Every segment reabsorbs water: the flow leaving a segment is
the fraction `f_tub_<seg>` of the GFR (Scotcher 2016), so the tubular flows
follow the GFR of the individual and the renal function. The drug moves with
the fluid, `GFR*f_tub_<seg>*Ctub_<seg>_<s>` out of a segment, and
concentrates as the water is reabsorbed (the filtrate has no proteins, the
lumen concentration is unbound).

The lumens are well-stirred compartments of constant volume (`Vtub_<seg>`,
Huang 2018, both kidneys). They are lumen, not tissue: they are not carved
out of the kidney volume, the blood and the tissue of the kidney keep the
volumes of the standard kidney and the tubule adds its fluid (about 0.12 l)
to the kidney.

Without secretion and reabsorption every filtered molecule reaches the
urine: at steady state the flow out of every segment equals the flow into
it, so the urinary excretion is the filtration `fu*GFR*Cpl` of the standard
kidney, whatever the volumes and the flows of the segments; the tubule
delays the excretion by minutes.

The fraction reabsorbed of a filtered substance (`f_reabsorbed_<s>` of
`Excretion.fraction_reabsorbed`) is an empirical passive reabsorption at the
end of the collecting duct: that fraction of the drug leaving it returns
into the plasma of the kidney (`J_reab_<s>`), the rest is urine, so the
renal clearance without transport is `fu*GFR*(1 - f_reabsorbed)` as in the
standard kidney. With secretion the fraction acts on the secreted drug as
well, the lumped `(fu*GFR + CL_sec)*(1 - F_reabs)` of the literature. A
reabsorption by a transporter is an apical influx instead.

The chemistry of the kidney is placed by membrane: a metabolism in the cells
of the kidney and a basolateral transport between interstitium and cells as
in the standard kidney; an apical transport between the cells and the lumen
of the proximal tubule, an efflux secretes into it, an influx reabsorbs from
it on the lumen concentration (the inhibitors of its protein act from the
lumen of the proximal tubule). A custom reaction states its rate for a part
of the standard kidney: the nephron kidney refuses it (R34). Every
metabolism and transport runs with the tubular capacity of the standard
kidney, `1 - f_tubular + f_tubular*f_renal_function`. A substance has
the species of the tubule when it is filtered, transported apically in the
kidney or inhibits the protein of an apical influx of the kidney.

The blood side (perfusion, tissue distribution, blood cells, the unbound
observers, the outflow concentrations `Cout_<s>`) is the code of the
standard kidney (`StandardBuilder`), and the boundary is that of the
standard kidney: the same inputs replaced by the same elements of the body
and the same outputs. A substance which is secreted but not filtered reads
the GFR in addition, for the flow of the tubule.
"""

import itertools
from dataclasses import dataclass

from sbmlutils.factory import Compartment, ModelDefinition, Species
from sbmlutils.metadata import SBO

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.chemistry import (
    CustomReaction,
    Direction,
    Membrane,
    Metabolism,
    Transport,
)
from pkpdsim.pbpk.data import cite, value
from pkpdsim.pbpk.organ import Boundary, OrganContext, StandardBuilder
from pkpdsim.pbpk.organ_builder import ONE
from pkpdsim.pbpk.organs import Organ, Part
from pkpdsim.pbpk.parts import Parts
from pkpdsim.pbpk.placement import Apical, place
from pkpdsim.pbpk.substance import Substance
from pkpdsim.pbpk.units import U

__all__ = ["NephronKidney"]

_SEGMENTS: tuple[tuple[str, str, str], ...] = tuple(
    zip(
        ids.TUBULE_SEGMENTS,
        ("proximal_tubule", "loop_of_henle", "distal_tubule", "collecting_duct"),
        ("proximal tubule", "loop of Henle", "distal tubule", "collecting duct"),
        strict=True,
    )
)
"""The tubule segments in the order of the flow: token, table of the data, name."""

_DATA = "kidney_human"
"""The data file of the nephron."""

_LUMEN = (
    "Lumen of the tubule, not tissue: not part of the kidney volume, the blood "
    "and the tissue of the kidney keep the volumes of the standard kidney. The "
    "volume sets the delay of the excretion, not the renal clearance at "
    "steady state."
)
"""The notes of a lumen volume."""

_APICAL = Apical(
    name="lumen of the proximal tubule", species=lambda s: ids.tubule(s, "pt")
)
"""The apical side of the cells of the kidney."""


@dataclass(frozen=True)
class NephronKidney:
    """The kidney with glomerular filtration into the segments of the tubule.

    Attributes:
        name: the name of the component, the prefix of the id of the model
            definition (`nephron_ki`).
    """

    name: str = "nephron"

    def build(self, ctx: OrganContext) -> tuple[ModelDefinition, Boundary]:
        """Build the model of the kidney and its boundary.

        Args:
            ctx: the context of the kidney.

        Returns:
            The model definition and the boundary of the standard kidney.

        Raises:
            ValueError: if the organ is not the kidney, a process cannot be
                placed (a custom reaction), or the standard kidney refuses a
                substance.
        """
        return NephronBuilder(ctx, self.name).build()


class NephronBuilder(StandardBuilder):
    """The SBML objects of the nephron kidney: the standard kidney, the tubule and its chemistry."""

    def __init__(self, ctx: OrganContext, name: str) -> None:
        """Start the collection of the elements of the nephron kidney.

        Args:
            ctx: the context of the kidney.
            name: the name of the component.

        Raises:
            ValueError: for another organ than the kidney, or a custom
                reaction of the kidney.
        """
        if ctx.organ is not Organ.KIDNEY:
            raise ValueError(
                f"The nephron kidney is a variant of the kidney, not of organ "
                f"'{ctx.organ.value}'."
            )
        super().__init__(ctx, name)
        self.processes = self._placed()
        self.tubular = self._tubular()

    def _notes(self) -> str:
        """The notes of the model definition."""
        return (
            f"The {self.name} kidney: the standard kidney with glomerular filtration "
            "into the lumen of the proximal tubule, the tubular flow through "
            "proximal tubule, loop of Henle, distal tubule and collecting duct with "
            "the reabsorption of water, apical transport between the cells and the "
            "lumen of the proximal tubule; the collecting duct drains into the "
            "urine; the boundary of the standard kidney.\n\nSources:\n\n"
            f"- {cite('Scotcher2016')}\n- {cite('Huang2018')}"
        )

    def _placed(self) -> list[Metabolism | Transport]:
        """The processes of the kidney.

        Raises:
            ValueError: for a custom reaction.
        """
        processes: list[Metabolism | Transport] = []
        for p in self._processes():
            if isinstance(p, CustomReaction):
                raise ValueError(
                    f"Custom reaction '{p.sid}' of the kidney: the nephron kidney "
                    "places metabolism and transport only, a custom reaction "
                    "states its rate for a part of the standard kidney."
                )
            processes.append(p)
        return processes

    def _tubular(self) -> set[str]:
        """The substances with species in the tubule.

        The filtered substances, the substances of an apical transport and
        the inhibitors of the protein of an apical influx (which act from the
        lumen of the proximal tubule).
        """
        apical = [
            p
            for p in self.processes
            if isinstance(p, Transport) and p.membrane is Membrane.APICAL
        ]
        influx = {
            p.protein
            for p in apical
            if p.direction is Direction.INFLUX and p.protein is not None
        }
        return (
            {s.sid for s in self.ctx.substances if s.excretion.urine}
            | {p.substance for p in apical}
            | {sid for sid, inh in self.ctx.inhibitions if inh.protein in influx}
        )

    @property
    def _gfr(self) -> str:
        """The glomerular filtration rate, which the tubular flows are fractions of."""
        return self._input(ids.GFR, U.l_per_min, "glomerular filtration rate [l/min]")

    def _compartments(self) -> None:
        """The compartments of the standard kidney, the lumens of the tubule and its flows."""
        super()._compartments()
        if not self.tubular:
            return
        for seg, key, name in _SEGMENTS:
            cid = ids.tubule_volume(seg)
            v = value(_DATA, key, "volume")
            self.compartments[cid] = self._annotate(
                Compartment(
                    cid,
                    v.value,
                    U.litre,
                    constant=True,
                    name=f"lumen {name} kidney",
                    sboTerm=SBO.PHYSICAL_COMPARTMENT,
                    notes=f"{_LUMEN} {v.notes}\n\nSource:\n\n- {cite(str(v.source))}",
                    port=True,
                ),
                "tubule",
                segment=seg,
            )
            self._constant(
                ids.tubule_flow_fraction(seg),
                value(_DATA, key, "flow_fraction"),
                U.dimensionless,
                f"tubular flow leaving the {name} relative to the GFR",
                f"The tubular fluid flow leaving the {name} as fraction of the "
                "glomerular filtration rate; the water reabsorbed upstream is "
                "the rest. An experiment may change it (e.g. the urine flow out "
                "of the collecting duct).",
                kind="tubule_fraction",
                segment=seg,
            )

    def _source(self, parts: Parts) -> Part:
        """The part of the kidney which is filtered: plasma, the organ when flow-limited."""
        return Part.ORGAN if parts.tissue == "organ" else Part.PLASMA

    def _substance(self, s: Substance) -> None:
        """The species, flows and exchanges of a substance of the standard kidney and its tubule."""
        super()._substance(s)
        if s.sid in self.tubular:
            self._tubule(s)

    def _filtration(self, s: Substance, parts: Parts, fu: str, cplasma: str) -> None:
        """The glomerular filtration into the lumen of the proximal tubule, `fu*GFR*Cpl`.

        The unbound concentration of the plasma leaving the kidney is
        filtered (of the flow-limited kidney `Cor/K_or`); the fraction
        reabsorbed acts at the end of the collecting duct (`_tubule`).
        """
        source = ids.species(self._source(parts), s.sid)
        self._reaction(
            ids.filtration_reaction(s.sid),
            f"{source} -> {ids.tubule(s.sid, 'pt')}",
            f"{fu} * {self._gfr} * {cplasma}",
            f"glomerular filtration {s.sid} into the proximal tubule",
            kind="filtration",
            segment="pt",
        )

    def _tubule(self, s: Substance) -> None:
        """The species of a substance in the lumens of the tubule and the tubular flows.

        `GFR*f_tub_<seg>*Ctub_<seg>_<s>` out of every segment into the next,
        out of the collecting duct into the urine; of a filtered substance
        the fraction `f_reabsorbed_<s>` of the flow out of the collecting
        duct returns into the plasma of the kidney instead.
        """
        gfr = self._gfr
        self._input_compartment(ids.VURINE, "urine")
        urine = self._input_species(
            ids.urine(s.sid), ids.VURINE, f"{s.name} (urine)", amount=True
        )
        for seg, _, name in _SEGMENTS:
            sid = ids.tubule(s.sid, seg)
            self.species[sid] = self._annotate(
                Species(
                    sid,
                    ids.tubule_volume(seg),
                    initialConcentration=0.0,
                    substanceUnit=U.mmole,
                    hasOnlySubstanceUnits=False,
                    name=f"{s.name} (lumen {name}, kidney)",
                    sboTerm=SBO.SIMPLE_CHEMICAL,
                    port=True,
                ),
                "concentration_tubule",
                segment=seg,
                substance=s,
            )
        for (seg, _, name), (downstream, _, downstream_name) in itertools.pairwise(
            _SEGMENTS
        ):
            c = ids.tubule(s.sid, seg)
            self._reaction(
                ids.tubule_reaction(s.sid, seg),
                f"{c} -> {ids.tubule(s.sid, downstream)}",
                f"{gfr} * {ids.tubule_flow_fraction(seg)} * {c}",
                f"tubular flow {s.sid} from the {name} into the {downstream_name}",
                kind="tubular_flow",
                segment=seg,
            )
        cd = ids.tubule(s.sid, "cd")
        leaving = f"{gfr} * {ids.tubule_flow_fraction('cd')} * {cd}"
        if not s.excretion.urine:
            self._reaction(
                ids.tubule_reaction(s.sid, "cd"),
                f"{cd} -> {urine}",
                leaving,
                f"urinary excretion {s.sid} from the collecting duct",
                kind="tubular_flow",
                segment="cd",
            )
            return
        reabsorbed = self._input(
            ids.fraction_reabsorbed(s.sid),
            U.dimensionless,
            f"fraction reabsorbed {s.sid}",
            "The fraction of the drug leaving the collecting duct which is "
            "reabsorbed into the plasma of the kidney (passive, empirical).",
        )
        self._reaction(
            ids.tubule_reaction(s.sid, "cd"),
            f"{cd} -> {urine}",
            f"{leaving} * ({ONE} - {reabsorbed})",
            f"urinary excretion {s.sid} from the collecting duct",
            kind="tubular_flow",
            segment="cd",
        )
        source = ids.species(self._source(self.ctx.parts[s.sid]), s.sid)
        self._reaction(
            ids.tubule_reabsorption_reaction(s.sid),
            f"{cd} -> {source}",
            f"{leaving} * {reabsorbed}",
            f"reabsorbed fraction {s.sid} from the collecting duct into the plasma",
            kind="reabsorption",
            segment="cd",
        )

    def _chemistry(self) -> None:
        """The processes of the kidney by membrane, apical between the cells and the proximal tubule."""
        for p in self.processes:
            if p.protein is not None:
                self._activity(p.protein)
        capacity = self._tubular_capacity() if self.processes else None
        for p in self.processes:
            self.reactions += place(
                p, self.ctx, self.ctx.parts, capacity, apical=_APICAL
            )
