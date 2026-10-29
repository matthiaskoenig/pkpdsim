"""The organ components: their interface and the standard organ.

An organ component builds the comp model definition of one organ and its
boundary to the body. The standard organ (PK-Sim 4-compartment structure,
Willmann 2005, Schmitt 2008) has the parts

```
plasma          V*f_vas*(1-HCT)    Cpl
blood cells     V*f_vas*HCT        Cbc   (kinetic erythrocytes)
interstitial    V*f_int            Cis   (full form)
intracellular   V*(1-f_vas-f_int)  Cic   (full form)
tissue          V*(1-f_vas)        Cti   (merged tissue)
organ           V                  Cor   (flow-limited)
```

of which only the parts a substance needs are built (`parts_for`). Every
species is a concentration, every reaction an amount per time. The organ owns
its inflow from the arterial pool (and the portal vein for the liver) and its
outflow into the venous pool (the portal vein for the portal organs, the
arterial pool for the lung); it reads every element of the body (volumes,
flows, hematocrit, fractions, substance parameters, pool and sink species)
through a port, and `Boundary.replacements` names the element of the body
which replaces it. The Kp of a substance is the whole-tissue (interstitial +
intracellular) to plasma coefficient in every form.

The liver has a portosystemic shunt, the fraction `f_shunt` of its arterial
and portal inflow which flows past it into the venous pool, and a functional
fraction `f_functional` of its tissue, which scales the tissue volumes and
the capacity of its processes (both inputs of the top model, healthy 0 and 1).
The kidney filters the unbound afferent (arterial) plasma concentration of a
substance with urinary excretion into the urine,
`J_gfr = fu*GFR*Car*(1 - f_reabsorbed)`; the plasma leaving the kidney lacks
the net filtrate at the outflow concentration,
`Qki*(1-HCT)*Cout - fu*GFR*(1 - f_reabsorbed)*Cout`, so the filtration
clearance of the arterial plasma is `fu*GFR*(1 - f_reabsorbed)` and a
secretion clears the plasma flow less the filtrate. A saturable reabsorption
(`ThresholdReabsorption`) replaces the constant fraction by the fraction
`1/(1 + (fu*C/reab_threshold)^reab_n)^(1/reab_n)` at the concentration of
each term (`Car` in the filtration, `Cout` in the outflow). Its
metabolism and transport (the tubular processes) run with the tubular
capacity `1 - f_tubular + f_tubular*f_renal_function`, which follows the
renal function (the intact nephron hypothesis, `f_tubular` 1). The kidney
exposes the urinary excretion rate `Rurine_<s>`, the sum of its fluxes into
the urine, as an output.

The volumes of the parts, the coefficients and the exchanges depend on
parameters only and are constants with initial assignments (they follow a
change of the body before a simulation, not during it); only what reads a
species is an assignment rule: the outflow concentrations and the unbound
concentrations of the tissue parts `Cu_<part>_<s>` (observers).

The elements (inputs with their ports, compartments, species, rules,
constants) are collected by `pkpdsim.pbpk.organ_builder.OrganBuilder`, which
the variants reuse. The standard organ is one well-stirred zone; the zonated
liver (`pkpdsim.pbpk.variants.liver_zonated`) builds the same parts, flows and
exchanges in each of its zones, which the blood flows through in series.
"""

import copy
import itertools
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Protocol, runtime_checkable

from sbmlutils.factory import ModelDefinition, Port, PortType, Sbase, Species
from sbmlutils.factory.units import UnitType
from sbmlutils.metadata import SBO

from pkpdsim.pbpk import distribution, ids
from pkpdsim.pbpk.annotation import annotate
from pkpdsim.pbpk.chemistry import (
    CustomReaction,
    Direction,
    Membrane,
    Metabolism,
    Transport,
)
from pkpdsim.pbpk.data import value
from pkpdsim.pbpk.organ_builder import (
    CAPACITY,
    COMPARTMENT_KINDS,
    ONE,
    OrganBuilder,
    label,
)
from pkpdsim.pbpk.organs import Organ, Part, Pool, inflow, outflow
from pkpdsim.pbpk.parts import Parts
from pkpdsim.pbpk.physiology import Sex
from pkpdsim.pbpk.placement import (
    apical_sink,
    built_parts,
    in_zone,
    kp_id,
    place,
    unbound,
    with_modifiers,
)
from pkpdsim.pbpk.proteins import Induction, Inhibition
from pkpdsim.pbpk.routes import Route
from pkpdsim.pbpk.slow_pool import SlowPool
from pkpdsim.pbpk.substance import Substance
from pkpdsim.pbpk.units import U

if TYPE_CHECKING:
    from pkpdsim.pbpk.definition import PBPKDefinition

__all__ = [
    "Boundary",
    "GutComponent",
    "OrganComponent",
    "OrganContext",
    "StandardBuilder",
    "StandardOrgan",
]


@dataclass(frozen=True)
class OrganContext:
    """What an organ component is built from.

    Attributes:
        organ: the organ.
        substances: the substances of the definition.
        parts: the parts of the organ by substance id.
        chemistry: the processes of the definition.
        inhibitions: the inhibitions as (inhibitor id, inhibition).
        inductions: the inductions as (inducer id, induction).
        sex: the sex of the build, whose reference individual the checks of
            the organ use (`distribution.check_kp`).
        routes: the routes of administration (the intestine applies the
            oral doses).
        top: the parameters and compartments of the top model which exist
            when the organs are built, id -> (`parameter` or `compartment`,
            unit): the ids a `Formula` of the organ may read, through an
            input port each.
        slow_pools: the slow pools of the definition (`SlowPool`); the
            standard organ builds those of its organ.
    """

    organ: Organ
    substances: tuple[Substance, ...]
    parts: Mapping[str, Parts]
    chemistry: tuple[Metabolism | Transport | CustomReaction, ...]
    inhibitions: tuple[tuple[str, Inhibition], ...]
    inductions: tuple[tuple[str, Induction], ...]
    sex: Sex
    routes: tuple[Route, ...] = ()
    top: Mapping[str, tuple[str, UnitType]] = field(default_factory=dict)
    slow_pools: tuple[SlowPool, ...] = ()


@dataclass(frozen=True)
class Boundary:
    """What an organ component exposes to the body.

    Attributes:
        replacements: the elements of the body which replace an element of the
            organ, local id -> id in the body: organ volume, flows, `HCT`,
            fractions, substance parameters, pool and sink species and their
            compartments.
        outputs: the outputs of the organ, local id -> id in the body
            (`Cout_<s>` -> `Cout_<o>_<s>`).
    """

    replacements: Mapping[str, str]
    outputs: Mapping[str, str]


class OrganComponent(Protocol):
    """An organ of the body, the standard one or a variant."""

    @property
    def name(self) -> str:
        """The name of the component (variant)."""
        ...

    def build(self, ctx: OrganContext) -> tuple[ModelDefinition, Boundary]:
        """Build the model of the organ and its boundary."""
        ...


@runtime_checkable
class GutComponent(OrganComponent, Protocol):
    """An intestine which replaces several organs with one submodel (`GU`).

    The segmented and the single lumen intestine replace small and large
    intestine: they receive the arterial inflow of both and drain into the
    portal vein, own the intestinal lumen, apply the oral doses and take the
    bile; the body gets the parameters they read from it (`body_objects`).
    """

    @property
    def replaces(self) -> tuple[Organ, ...]:
        """The organs the component replaces, the small intestine first."""
        ...

    def body_objects(self, definition: "PBPKDefinition") -> list[Sbase]:
        """Get the parameters of the top model which the component reads."""
        ...


@dataclass(frozen=True)
class StandardOrgan:
    """The standard organ with the collapsing of its parts.

    Attributes:
        name: the name of the component, the prefix of the id of the model
            definition (`standard_li`).
    """

    name: str = "standard"

    def build(self, ctx: OrganContext) -> tuple[ModelDefinition, Boundary]:
        """Build the model of the organ and its boundary.

        Args:
            ctx: the context of the organ.

        Returns:
            The model definition and the boundary.

        Raises:
            ValueError: if a substance cannot be built in the organ: a Kp not
                above the interstitial share in the full form, an exchange
                without `ftissue`, a process which cannot be placed.
        """
        return StandardBuilder(ctx, self.name).build()


_IMPAIRMENT = "impairment"
"""The data file of the impairment, which holds the tubular capacity."""


_FTISSUE = (
    "ftissue multiplies `fu*(Cpl - Cti/Kp)` (the unbound plasma concentration "
    "and the tissue concentration as plasma equivalent), so its fitted values "
    "differ from those of the hand-written models, where it multiplied the "
    "total concentrations."
)
"""The meaning of ftissue in the notes."""


_TISSUE_PARTS = {
    Part.TISSUE: "tissue",
    Part.INTERSTITIAL: "interstitium",
    Part.INTRACELLULAR: "intracellular space",
    Part.ORGAN: "organ",
}
"""The parts of the tissue with an unbound concentration observer, in names."""


class StandardBuilder(OrganBuilder):
    """The SBML objects of one standard organ: its parts, flows, exchanges and chemistry.

    The organ is one well-stirred zone with the ids without zone. A zonated
    organ gives its zones (`zones`): the blood flows from the pools into the
    first zone, from each zone into the next and from the last one into the
    pool (the outflow concentration is that of the last zone); every zone has
    the parts of the organ with the share `f_zone` (`1/n`) of their volumes
    and of the passive exchanges, with the coefficients of the organ. The
    intestine variants derive from it for the intestinal wall.

    Attributes:
        zones: the zones of a zonated organ in the order of the blood flow,
            none for the organ without zones.
    """

    zones: tuple[str, ...] = ()

    def build(self) -> tuple[ModelDefinition, Boundary]:
        """Build the model definition and the boundary."""
        self._compartments()
        for s in self.ctx.substances:
            self._substance(s)
        self._chemistry()
        self._slow_pools()
        if self.organ is Organ.KIDNEY:
            self._urine_rates()
        md = self.model_definition(self._notes())
        boundary = Boundary(
            replacements=dict(self.replacements), outputs=dict(self.outputs)
        )
        return md, boundary

    def _notes(self) -> str:
        """The notes of the model definition."""
        return (
            f"The {self.name} organ {label(self.organ)}: the parts which the "
            "substances need, their blood flows, exchanges and chemistry."
        )

    @property
    def _each_zone(self) -> tuple[str | None, ...]:
        """The zones in the order of the blood flow, `(None,)` for the organ without zones."""
        return self.zones or (None,)

    def _zone_fraction(self) -> str:
        """The volume fraction of a zone, `f_zone = 1/n` of `n` zones of equal volume."""
        return self._constant(
            ids.F_ZONE,
            1.0 / len(self.zones),
            U.dimensionless,
            "volume fraction of a zone",
            f"The organ has {len(self.zones)} zones of equal volume in series; a "
            "zone has this share of the volumes and of the passive exchanges of "
            "the organ. A constant of the structure of the model, not a "
            "parameter to change.",
            kind="zone_fraction",
        )

    def _share(self, quantity: str, zone: str | None) -> str:
        """The share of a zone of a volume or an exchange of the organ, without zone all of it."""
        if zone is None:
            return quantity
        return f"{quantity} * {self._zone_fraction()}"

    def _flows(self, s: Substance) -> tuple[list[tuple[Pool, str]], str]:
        """The inflows (pool, flow) of a substance into the organ and the flow of its outflow."""
        del s
        if self.organ is Organ.LIVER:
            return [(Pool.ARTERIAL, ids.QHA), (Pool.PORTAL, ids.QPO)], ids.QH
        flow = ids.organ_flow(self.organ)
        return [(pool, flow) for pool in inflow(self.organ)], flow

    def _capacity(self) -> str | None:
        """The functional fraction of the tissue: `f_functional` in the liver.

        The non-functional tissue of the liver (fibrosis, necrosis) holds no
        substance and has no capacity: `f_functional` multiplies the tissue
        volumes, the tissue coefficient of the flow-limited liver, the
        exchange share of the liver and the rate of every metabolism and
        transport of the liver (`place`). The other organs have no such
        fraction.

        Returns:
            The id of the input `f_functional` in the liver, `None` elsewhere.
        """
        if self.organ is not Organ.LIVER:
            return None
        return self._input(
            ids.F_FUNCTIONAL,
            U.dimensionless,
            "functional liver fraction",
            "The functional fraction of the liver tissue.",
        )

    def _process_capacity(self) -> str | None:
        """The factor of the capacity of every metabolism and transport of the organ.

        The functional fraction `f_functional` in the liver (`_capacity`),
        the tubular capacity in the kidney (`_tubular_capacity`), none
        elsewhere.

        Returns:
            The factor (formula), `None` for none.
        """
        if self.organ is Organ.KIDNEY:
            return self._tubular_capacity()
        return self._capacity()

    def _tubular_capacity(self) -> str:
        """The capacity of the tubular processes: `1 - f_tubular + f_tubular*f_renal_function`.

        Every metabolism and transport of the kidney runs with the share
        `f_tubular` of its capacity which follows the renal function, the
        intact nephron hypothesis (a model choice of `data/impairment.toml`,
        `[ckd] tubular_capacity`): with the share 1 a tubular secretion falls
        with the GFR, as the secretion of a fraction of the GFR of the
        hand-written models did. The filtration `fu*GFR` follows the GFR
        itself.

        Returns:
            The factor (formula).
        """
        f = self._input(
            ids.F_RENAL_FUNCTION,
            U.dimensionless,
            "renal function",
            "The GFR relative to the normal GFR of the individual.",
        )
        share = self._constant(
            ids.F_TUBULAR,
            value(_IMPAIRMENT, "ckd", "tubular_capacity"),
            U.dimensionless,
            "share of the tubular capacity which follows the renal function",
            "The capacity of every metabolism and transport of the kidney is "
            f"`1 - {ids.F_TUBULAR} + {ids.F_TUBULAR}*{ids.F_RENAL_FUNCTION}` of "
            "the healthy kidney.",
            kind="tubular_share",
        )
        return f"({ONE} - {share} + {share} * {f})"

    def _functional(self, tissue: str) -> str:
        """A quantity of the tissue of the organ, of the liver its functional share."""
        capacity = self._capacity()
        return tissue if capacity is None else f"{tissue} * {capacity}"

    def _functional_text(self) -> str:
        """The notes of a quantity of the tissue scaled by the functional fraction."""
        if self.organ is not Organ.LIVER:
            return ""
        return f"The functional share `{ids.F_FUNCTIONAL}` of the tissue of the liver."

    def _compartments(self) -> None:
        """The compartments of the parts which hold a species of any substance (in every zone)."""
        built = {p for parts in self.ctx.parts.values() for p in built_parts(parts)}
        v, f_vas = self.v, self.f_vas
        volumes: list[tuple[Part, str, str, str]] = []
        if Part.PLASMA in built:
            volumes.append(
                (Part.PLASMA, f"{v} * {f_vas} * ({ONE} - {self.hct})", "plasma", "")
            )
        if Part.BLOOD_CELLS in built:
            volumes.append(
                (Part.BLOOD_CELLS, f"{v} * {f_vas} * {self.hct}", "blood cells", "")
            )
        if Part.INTERSTITIAL in built:
            volumes.append(
                (
                    Part.INTERSTITIAL,
                    self._functional(f"{v} * {self.f_int}"),
                    "interstitium",
                    self._functional_text(),
                )
            )
            volumes.append(
                (
                    Part.INTRACELLULAR,
                    self._functional(f"{v} * ({ONE} - {f_vas} - {self.f_int})"),
                    "intracellular space",
                    self._functional_text(),
                )
            )
        if Part.TISSUE in built:
            volumes.append(
                (
                    Part.TISSUE,
                    self._functional(f"{v} * ({ONE} - {f_vas})"),
                    "tissue",
                    self._functional_text(),
                )
            )
        if Part.ORGAN in built:
            volumes.append((Part.ORGAN, v, "organ", ""))
        for zone in self._each_zone:
            for part, formula, name, text in volumes:
                self._compartment(
                    ids.compartment(part, zone),
                    self._share(formula, zone),
                    name,
                    text,
                    zone=zone,
                    kind=COMPARTMENT_KINDS[part],
                )

    def _substance(self, s: Substance) -> None:
        """The species, outflow concentration, flows and exchanges of a substance."""
        parts = self.ctx.parts[s.sid]
        if parts.tissue == "organ" and parts.blood_cells:
            raise ValueError(
                f"{s.sid}: kinetic erythrocytes need a distribution with exchange, "
                f"not the flow-limited organ '{self.organ.value}'."
            )
        fu = self._fu(s)
        for zone in self._each_zone:
            if parts.flow_carrier == "blood_bp" and parts.tissue != "organ":
                self._compartment(
                    ids.plasma_capacity(s.sid, zone),
                    self._share(f"{self.v} * {self.f_vas} * {self._bp(s)}", zone),
                    f"plasma-equivalent capacity V*f_vas*BP of {s.sid},",
                    CAPACITY,
                    zone=zone,
                    kind="capacity",
                )
            for part in built_parts(parts):
                self._species(part, s, parts, zone)
        k_organ = self._k_organ(s, parts) if parts.tissue == "organ" else None
        last = self._each_zone[-1]
        cout = self._leaving(s, Part.PLASMA, last, k_organ)
        out = ids.outflow_conc(s.sid)
        self._rule(
            out,
            cout,
            U.mM,
            f"outflow plasma concentration {s.sid} {label(self.organ)} [mM]",
            "Plasma concentration of the blood leaving the organ.",
            port=Port(f"{out}_port", portType=PortType.OUTPUT_PORT),
            kind="outflow",
        )
        self.outputs[out] = ids.organ_outflow(self.organ, s.sid)
        if parts.blood_cells:
            out_bc = ids.outflow_conc_bc(s.sid)
            self._rule(
                out_bc,
                self._leaving(s, Part.BLOOD_CELLS, last, k_organ),
                U.mM,
                f"outflow blood cell concentration {s.sid} {label(self.organ)} [mM]",
                "Blood cell concentration of the blood leaving the organ.",
                port=Port(f"{out_bc}_port", portType=PortType.OUTPUT_PORT),
                kind="outflow_blood_cells",
            )
            self.outputs[out_bc] = ids.organ_outflow_bc(self.organ, s.sid)
        filters = self.organ is Organ.KIDNEY and s.excretion.urine
        self._blood_flows(
            s, parts, k_organ, self._filtrate(s, fu, out) if filters else None
        )
        if parts.blood_cells:
            for zone in self._each_zone:
                self._blood_cells(s, zone)
        if filters:
            self._filtration(s, parts, fu)
        for zone in self._each_zone:
            if parts.tissue == "merged":
                self._merged(s, fu, zone)
            elif parts.tissue == "full":
                self._full(s, fu, zone)
            self._unbound_observers(s, parts, zone)

    def _leaving(
        self, s: Substance, part: Part, zone: str | None, k_organ: str | None
    ) -> str:
        """The concentration of the blood leaving a zone (the organ), plasma or blood cells.

        The plasma concentration, of the flow-limited organ `Cor/K_or` (in
        equilibrium with the plasma leaving it, `k_organ` the id of its organ
        to plasma coefficient, `None` for an organ with a plasma species); the
        blood cell concentration.
        """
        if part is Part.BLOOD_CELLS:
            return ids.species(Part.BLOOD_CELLS, s.sid, zone)
        if k_organ is not None:
            return f"{ids.species(Part.ORGAN, s.sid, zone)} / {k_organ}"
        return ids.species(Part.PLASMA, s.sid, zone)

    def _unbound_observers(
        self, s: Substance, parts: Parts, zone: str | None = None
    ) -> None:
        """The unbound concentrations of the tissue parts of a substance, `Cu_<part>_<s>`.

        `fu/K*C` with the coefficient of the part to plasma (`Kp`,
        `K_int_pls`, `K_cell_pls`, `K_or`): observers, assignment rules of
        the organ (they read species), the concentration which the processes
        of the part act on; in a zone `Cu_<part>_<zone>_<s>`.
        """
        for part in built_parts(parts):
            if part not in _TISSUE_PARTS:
                continue
            self._rule(
                ids.unbound(part, s.sid, zone),
                unbound(part, s, self.organ, zone),
                U.mM,
                in_zone(
                    f"unbound concentration {s.sid} {_TISSUE_PARTS[part]} "
                    f"{label(self.organ)}",
                    zone,
                )
                + " [mM]",
                f"Unbound concentration of the {_TISSUE_PARTS[part]}; an observer.",
                kind="unbound_concentration",
            )

    def _k_organ(self, s: Substance, parts: Parts) -> str:
        """The organ to plasma coefficient of the flow-limited organ.

        `K_or = (Vpl + Vbc*Kbc + Vti*Kp)/V`; the blood cells hold the substance
        only with erythrocytes in equilibrium, `Vpl + Vbc*Kbc = V*f_vas*BP`.
        """
        kp = self._kp(s)
        assert kp is not None
        if parts.flow_carrier == "blood_bp":
            vascular = self._bp(s)
        else:
            vascular = f"({ONE} - {self.hct})"
        return self._initial(
            ids.k_organ(s.sid),
            distribution.k_organ(self.f_vas, vascular, self._functional(kp)),
            U.dimensionless,
            f"organ to plasma coefficient {s.sid}",
            "Organ to plasma coefficient of the flow-limited organ, "
            "`(Vpl + Vbc*Kbc + Vti*Kp)/V`: the organ is in equilibrium with the "
            f"plasma leaving it. {self._functional_text()}".rstrip(),
            kind="partition",
        )

    def _blood_flows(
        self,
        s: Substance,
        parts: Parts,
        k_organ: str | None,
        filtrate: str | None = None,
    ) -> None:
        """The inflow from the pools and the outflow into the pool, plasma and blood cells.

        The plasma carries `Q*(1-HCT)*C`, or `Q*BP*C` with erythrocytes in
        equilibrium; kinetic blood cells carry `Q*HCT*Cbc`. The plasma leaving
        the kidney carries the flow of its carrier less the net filtrate at
        the outflow concentration, `filtrate` (`_filtrate`):
        `Q*(1-HCT)*Cout - fu*GFR*(1 - f_reabsorbed)*Cout` (a saturable
        reabsorption at `Cout`). The outflow is valid while the net filtered
        plasma flow `fu*GFR*(1 - f_reabsorbed)` (`fu*GFR` for a saturable
        reabsorption, whose fraction falls to 0) is below the flow of the
        carrier `Q*(1-HCT)` (`Q*BP`), which holds for every physiological GFR
        (about 0.125 against 0.76 l/min in the reference individual),
        hyperfiltration and the classes of impairment included; the GFR, the
        flow and the hematocrit are parameters which an experiment changes, so
        `PBPKSimulationExperiment` checks it for every simulated individual
        (`pkpdsim.pbpk.experiment.kidney_outflow`). The fraction
        `f_shunt` of the inflow of the liver bypasses it (portosystemic shunt)
        from the arterial and the portal pool directly into the venous pool;
        the liver is perfused by and drains `(1 - f_shunt)*Q`. The zones of a
        zonated organ are perfused in series: the inflow enters the first
        zone, the flow leaving a zone enters the next, the last one drains;
        `k_organ` is the organ to plasma coefficient of the flow-limited
        organ (`_leaving`).
        """
        inflows, out_flow = self._flows(s)
        out_pool = outflow(self.organ)
        if parts.flow_carrier == "blood_bp":
            carrier = self._bp(s)
        else:
            carrier = f"({ONE} - {self.hct})"
        target = Part.ORGAN if parts.tissue == "organ" else Part.PLASMA
        carried = [(Part.PLASMA, carrier, target)]
        if parts.blood_cells:
            carried.append((Part.BLOOD_CELLS, self.hct, Part.BLOOD_CELLS))
        perfused, f_shunt = "", None
        if self.organ is Organ.LIVER:
            f_shunt = self._input(
                ids.F_SHUNT,
                U.dimensionless,
                "portosystemic shunt",
                "The fraction of the liver inflow which bypasses the liver.",
            )
            perfused = f"({ONE} - {f_shunt}) * "
        first, last = self._each_zone[0], self._each_zone[-1]
        for part, share, local_part in carried:
            what = "" if part is Part.PLASMA else " blood cells"
            local = ids.species(local_part, s.sid, first)
            out = self._pool_species(out_pool, part, s, parts)
            for pool, flow in inflows:
                c = self._pool_species(pool, part, s, parts)
                rate = f"{self._flow(flow)} * {share} * {c}"
                self._reaction(
                    ids.inflow_reaction(pool, part, s.sid),
                    f"{c} -> {local}",
                    f"{perfused}{rate}",
                    f"inflow {s.sid}{what} from {pool.name.lower()} pool",
                    kind="blood_flow",
                    pool=pool,
                )
                if f_shunt is not None:
                    self._reaction(
                        ids.shunt_reaction(pool, part, s.sid),
                        f"{c} -> {out}",
                        f"{f_shunt} * {rate}",
                        f"shunt {s.sid}{what} from {pool.name.lower()} pool past "
                        f"the liver into {out_pool.name.lower()} pool",
                        kind="blood_flow",
                        pool=pool,
                    )
            for upstream, downstream in itertools.pairwise(self.zones):
                self._reaction(
                    ids.sinusoid_reaction(part, s.sid, upstream),
                    f"{ids.species(local_part, s.sid, upstream)} -> "
                    f"{ids.species(local_part, s.sid, downstream)}",
                    f"{perfused}{self._flow(out_flow)} * {share} * "
                    f"{self._leaving(s, part, upstream, k_organ)}",
                    f"flow {s.sid}{what} from zone {upstream} into zone {downstream}",
                    kind="blood_flow",
                )
            leaving = (
                ids.outflow_conc(s.sid)
                if part is Part.PLASMA
                else ids.outflow_conc_bc(s.sid)
            )
            rate = f"{perfused}{self._flow(out_flow)} * {share} * {leaving}"
            if part is Part.PLASMA and filtrate is not None:
                rate = f"{rate} - {filtrate}"
            self._reaction(
                ids.outflow_reaction(part, s.sid),
                f"{ids.species(local_part, s.sid, last)} -> {out}",
                rate,
                f"outflow {s.sid}{what} into {out_pool.name.lower()} pool",
                kind="blood_flow",
                pool=out_pool,
            )

    def _fraction_reabsorbed(self, s: Substance) -> str:
        """The constant fraction `f_reabsorbed_<s>` of the filtrate which is reabsorbed, an input."""
        return self._input(
            ids.fraction_reabsorbed(s.sid),
            U.dimensionless,
            f"fraction reabsorbed {s.sid}",
            "The fraction of the filtrate which is reabsorbed.",
        )

    def _reabsorbed(self, s: Substance, fu: str, c: str) -> str:
        """The reabsorbed fraction of the filtrate of the plasma concentration `c`.

        The constant fraction `f_reabsorbed_<s>` (`Excretion.fraction_reabsorbed`),
        or of a saturable reabsorption (`Excretion.reabsorption`, `ThresholdReabsorption`)
        the fraction `1/(1 + (fu*|c|/reab_threshold_<s>)^reab_n_<s>)^(1/reab_n_<s>)`
        of the unbound concentration `fu*c` of the filtrate, at most 1, so the
        reabsorption never exceeds the filtration (`|c|` keeps the power
        defined for a concentration a little below 0, an undershoot of the
        integrator).
        """
        saturable = s.excretion.reabsorption
        if saturable is None:
            return self._fraction_reabsorbed(s)
        threshold = self._input(
            ids.reabsorption_threshold(s.sid),
            U.mM,
            f"renal threshold {s.sid} [mM]",
            "The renal threshold of the saturable reabsorption, a concentration "
            "of the filtrate: above it the reabsorption approaches the transport "
            "maximum `GFR` times the threshold.",
        )
        n = self._input(
            ids.reabsorption_n(s.sid),
            U.dimensionless,
            f"sharpness of the renal threshold {s.sid}",
            "How sharp the renal threshold of the saturable reabsorption is.",
        )
        # abs: a concentration a little below 0 (an undershoot of the
        # integrator) keeps the power of a non-integer n defined
        return f"{ONE} / ({ONE} + ({fu} * abs({c}) / {threshold})^{n})^({ONE} / {n})"

    def _filtrate(self, s: Substance, fu: str, c: str) -> str:
        """The net filtrate of the plasma concentration `c`, `fu*GFR*c*(1 - f_reabsorbed)`.

        The unbound plasma is filtered with the GFR, the reabsorbed fraction
        (`_reabsorbed`, constant or saturable at `c`) returns into the plasma:
        the glomerular filtration of the afferent plasma (`_filtration`) and
        the filtrate which the plasma leaving the kidney lacks at the outflow
        concentration (`_blood_flows`).
        """
        gfr = self._input(ids.GFR, U.l_per_min, "glomerular filtration rate [l/min]")
        return f"{fu} * {gfr} * {c} * ({ONE} - {self._reabsorbed(s, fu, c)})"

    def _filtration(self, s: Substance, parts: Parts, fu: str) -> None:
        """The glomerular filtration into the urine, `fu*GFR*Car*(1 - f_reabsorbed)`.

        The glomerulus is upstream of the peritubular capillaries: it filters
        the afferent plasma, the arterial plasma concentration `Car_<s>` of
        the inflow (in every form and erythrocyte mode, a modifier), out of
        the plasma of the kidney (its organ species when flow-limited); the
        reabsorbed fraction of the filtrate (`f_reabsorbed_<s>`, or of a
        saturable reabsorption at `Car_<s>`) returns into the plasma (it is
        not taken out). The plasma leaving the kidney lacks the net
        filtrate (`_blood_flows`): at pure filtration the kidney plasma is the
        arterial plasma and the clearance is `fu*GFR*(1 - f_reabsorbed)`, a
        secretion clears the plasma flow less the filtrate.
        """
        self._input_compartment(ids.VURINE, "urine")
        urine = self._input_species(
            ids.urine(s.sid), ids.VURINE, f"{s.name} (urine)", amount=True
        )
        source = ids.species(
            Part.ORGAN if parts.tissue == "organ" else Part.PLASMA, s.sid
        )
        arterial = self._pool_species(Pool.ARTERIAL, Part.PLASMA, s, parts)
        self._reaction(
            ids.filtration_reaction(s.sid),
            with_modifiers(f"{source} -> {urine}", (arterial,)),
            self._filtrate(s, fu, arterial),
            f"glomerular filtration {s.sid}",
            kind="filtration",
        )

    def _blood_cells(self, s: Substance, zone: str | None = None) -> None:
        """The exchange between plasma and kinetic blood cells (of a zone).

        `kex_bc * Vbc * (Cpl - Cbc/Kbc)`: the blood cells approach `Kbc*Cpl`
        with the rate `kex_bc/Kbc`.
        """
        kbc = self._input(
            ids.kbc(s.sid), U.dimensionless, f"blood cell to plasma coefficient {s.sid}"
        )
        kex = self._input(
            ids.kbc_exchange(s.sid),
            U.per_min,
            f"exchange plasma and blood cells {s.sid} [1/min]",
        )
        cpl = ids.species(Part.PLASMA, s.sid, zone)
        cbc = ids.species(Part.BLOOD_CELLS, s.sid, zone)
        vbc = ids.compartment(Part.BLOOD_CELLS, zone)
        self._reaction(
            ids.exchange_reaction(Part.BLOOD_CELLS, s.sid, zone),
            f"{cpl} <-> {cbc}",
            f"{kex} * {vbc} * ({cpl} - {cbc} / {kbc})",
            in_zone(f"exchange {s.sid} plasma and blood cells", zone),
            kind="exchange",
        )

    def _ps(self, s: Substance) -> str:
        """The passive exchange between plasma and tissue of a substance in the organ.

        `PS = ftissue * Vti/Vti_total_ref`, the share of the whole-body exchange
        by the tissue volume relative to the tissue volume of the reference
        individual (so the exchange scales with the individual), or the value
        of the organ (`ps` of the Kp): per litre of organ `PS = kps * V<o>`
        (the default, it scales with the individual as well), else the
        exchange of the whole organ in l/min; no
        passive exchange for a substance which does not enter the tissue of
        the organ (it has the full form by its transporters).

        Raises:
            ValueError: if the substance enters the tissue but has neither
                `ftissue` nor a value of the organ.
        """
        pid = ids.ps(s.sid)
        if pid in self.parameters:
            return pid
        name = f"exchange plasma and tissue {s.sid} [l/min]"
        dist = s.distribution
        if dist is None or kp_id(s, self.organ) is None:
            return self._constant(
                pid,
                0.0,
                U.l_per_min,
                name,
                "No passive exchange: the substance does not distribute into the "
                "tissue of the organ, it enters the cells by its transporters only.",
                kind="permeability",
            )
        if self.organ in dist.ps:
            if not dist.ps_per_volume:
                whole = (
                    "Passive exchange of the whole organ, given instead of the "
                    "share of the whole-body exchange; it does not scale with the "
                    "organ volume."
                )
                if self.organ is not Organ.LIVER:
                    return self._constant(
                        pid,
                        dist.ps[self.organ],
                        U.l_per_min,
                        name,
                        whole,
                        kind="permeability",
                    )
                ref = self._constant(
                    ids.ps_ref(s.sid),
                    dist.ps[self.organ],
                    U.l_per_min,
                    f"exchange plasma and tissue {s.sid} of the healthy liver [l/min]",
                    whole,
                    kind="permeability",
                )
                return self._initial(
                    pid,
                    self._functional(ref),
                    U.l_per_min,
                    name,
                    f"{whole} {self._functional_text()}",
                    kind="permeability",
                )
            kps = self._constant(
                ids.kps(s.sid),
                dist.ps[self.organ],
                U.per_min,
                f"exchange per volume plasma and tissue {s.sid} [1/min]",
                "Passive exchange per litre of organ, given instead of the share "
                f"of the whole-body exchange; the organ has `{ids.kps(s.sid)} * "
                f"{self.v}`.",
                kind="permeability",
            )
            return self._initial(
                pid,
                f"{kps} * " + self._functional(self.v),
                U.l_per_min,
                name,
                "Passive exchange between plasma and tissue: the exchange per "
                "litre of organ times the organ volume, so it follows the organ "
                f"of the individual. {self._functional_text()}".strip(),
                kind="permeability",
            )
        if dist.ftissue is None:
            raise ValueError(
                f"{s.sid}: the exchange between plasma and tissue in organ "
                f"'{self.organ.value}' needs ftissue or a value of the organ "
                f"(ps), mode '{dist.mode.value}' has neither."
            )
        ftissue = self._input(
            ids.ftissue(s.sid),
            U.l_per_min,
            f"whole-body exchange {s.sid} [l/min]",
            _FTISSUE,
        )
        total = self._input(
            ids.tissue_volume_total_ref(s.sid),
            U.litre,
            f"reference tissue volume of the organs {s.sid} enters [l]",
            "Tissue volume of the reference individual summed over the organs "
            "whose exchange is a share of ftissue.",
        )
        return self._initial(
            pid,
            f"{ftissue} * "
            + self._functional(f"{self.v} * ({ONE} - {self.f_vas})")
            + f" / {total}",
            U.l_per_min,
            name,
            "Passive exchange between plasma and tissue: the share of the "
            "whole-body exchange `ftissue` by the tissue volume of the organ "
            f"relative to the reference tissue volume. {_FTISSUE}",
            kind="permeability",
        )

    def _merged(self, s: Substance, fu: str, zone: str | None = None) -> None:
        """The exchange between plasma and the merged tissue, `PS*fu*(Cpl - Cti/Kp)`.

        A zone has the share `f_zone` of the exchange.
        """
        kp = self._kp(s)
        assert kp is not None
        ps = self._ps(s)
        cpl = ids.species(Part.PLASMA, s.sid, zone)
        cti = ids.species(Part.TISSUE, s.sid, zone)
        self._reaction(
            ids.exchange_reaction(Part.TISSUE, s.sid, zone),
            f"{cpl} <-> {cti}",
            f"{self._share(ps, zone)} * {fu} * ({cpl} - {cti} / {kp})",
            in_zone(f"exchange {s.sid} plasma and tissue", zone),
            kind="exchange",
        )

    def _full(self, s: Substance, fu: str, zone: str | None = None) -> None:
        """The exchanges plasma - interstitium - intracellular space.

        `PS_vas*fu*(Cpl - Cis/K_int_pls)` and
        `PS*(fu/K_int_pls*Cis - fu/K_cell_pls*Cic)`, with `K_int_pls` of
        Schmitt 2008 and `K_cell_pls` such that the whole tissue has the Kp;
        a zone has the share `f_zone` of both exchanges.
        """
        k_int, k_cell = self._cell_coefficients(s, fu)
        if ids.PS_VAS not in self.parameters:
            self.parameters[distribution.K_VAS] = distribution.constant(
                distribution.K_VAS
            )
            self._initial(
                ids.PS_VAS,
                f"{distribution.K_VAS} * {self.v}",
                U.l_per_min,
                "exchange plasma and interstitium [l/min]",
                "Exchange between plasma and interstitium across the endothelium.",
                kind="permeability",
            )
        ps = self._ps(s)
        cpl = ids.species(Part.PLASMA, s.sid, zone)
        cis = ids.species(Part.INTERSTITIAL, s.sid, zone)
        cic = ids.species(Part.INTRACELLULAR, s.sid, zone)
        self._reaction(
            ids.exchange_reaction(Part.INTERSTITIAL, s.sid, zone),
            f"{cpl} <-> {cis}",
            f"{self._share(ids.PS_VAS, zone)} * {fu} * ({cpl} - {cis} / {k_int})",
            in_zone(f"exchange {s.sid} plasma and interstitium", zone),
            kind="exchange",
        )
        self._reaction(
            ids.exchange_reaction(Part.INTRACELLULAR, s.sid, zone),
            f"{cis} <-> {cic}",
            f"{self._share(ps, zone)} * ({fu} / {k_int} * {cis} - {fu} / {k_cell} "
            f"* {cic})",
            in_zone(f"passive exchange {s.sid} interstitium and cells", zone),
            kind="exchange",
        )

    def _cell_coefficients(self, s: Substance, fu: str) -> tuple[str, str]:
        """The interstitial and the intracellular to plasma coefficients of a substance.

        `K_int_pls` of Schmitt 2008 and `K_cell_pls` such that the whole
        tissue has the Kp, once for the organ (every zone has them).
        """
        if ids.k_cell_pls(s.sid) in self.parameters:
            return ids.k_int_pls(s.sid), ids.k_cell_pls(s.sid)
        for key in distribution.INTERSTITIAL:
            if key not in self.parameters:
                self.parameters[key] = distribution.constant(key)
        k_int = self._initial(
            ids.k_int_pls(s.sid),
            distribution.k_int_pls(fu),
            U.dimensionless,
            f"interstitial to plasma coefficient {s.sid}",
            "Interstitial to plasma coefficient (Schmitt 2008): the water of the "
            "interstitium and the plasma proteins of the interstitial fluid.",
            kind="partition",
        )
        kp = self._kp(s)
        if kp is None:
            k_cell = self._initial(
                ids.k_cell_pls(s.sid),
                k_int,
                U.dimensionless,
                f"intracellular to plasma coefficient {s.sid}",
                "The substance does not distribute into the tissue of the organ: "
                "without a Kp the cells bind it as the interstitial fluid.",
                kind="partition",
            )
        else:
            distribution.check_kp(s, self.organ, self.ctx.sex)
            k_cell = self._initial(
                ids.k_cell_pls(s.sid),
                distribution.k_cell_pls(kp, k_int, self.f_vas, self.f_int),
                U.dimensionless,
                f"intracellular to plasma coefficient {s.sid}",
                "Intracellular to plasma coefficient such that interstitium and "
                "cells together have the whole-tissue Kp, "
                "`(f_int*K_int_pls + f_cell*K_cell_pls)/(f_int + f_cell) = Kp`. "
                "The Kp must exceed the interstitial share "
                "`f_int*K_int_pls/(f_int + f_cell)` (the lower bound of a fitted "
                f"Kp); below it the coefficient is floored at "
                f"{distribution.K_CELL_FLOOR}.",
                kind="partition",
            )
        return k_int, k_cell

    def _chemistry(self) -> None:
        """The processes of the organ with the activities of their proteins and their sinks.

        Raises:
            ValueError: if a species of a custom reaction has the id of a
                species of the organ.
        """
        processes = self._processes()
        self._proteins_and_sinks(processes)
        for p in processes:
            if not isinstance(p, CustomReaction):
                continue
            for sp in p.species:
                sid = str(sp.sid)
                if sid in self.species:
                    raise ValueError(
                        f"Custom reaction '{p.sid}': species '{sid}' is a species "
                        f"of the organ '{self.organ.value}' already."
                    )
                self.species[sid] = annotate(
                    copy.copy(sp),
                    "custom_amount" if sp.hasOnlySubstanceUnits else "custom_species",
                )
        capacity = (
            self._process_capacity()
            if any(not isinstance(p, CustomReaction) for p in processes)
            else None
        )
        for p in processes:
            self.reactions += place(p, self.ctx, self.ctx.parts, capacity)

    def _slow_pools(self) -> None:
        """The slow pools of the organ: `Cslow_<s>` in the tissue and its exchange.

        `kslow_<s>*Vti*(Kpslow_<s>*Cti_<s> - Cslow_<s>)` with the constants of
        the top model (`SlowPool`); the pool of an endogenous substance starts
        at its equilibrium, `Kpslow_<s>*Kp*baseline_<s>`.

        Raises:
            ValueError: if the substance has no merged tissue in the organ or
                the organ has zones (`PBPKDefinition.validate` refuses both).
        """
        by_sid = {s.sid: s for s in self.ctx.substances}
        for pool in self.ctx.slow_pools:
            if self.organ not in pool.organs:
                continue
            s = by_sid[pool.substance]
            if self.ctx.parts[s.sid].tissue != "merged" or self.zones:
                raise ValueError(
                    f"{s.sid}: the slow pool of organ '{self.organ.value}' needs "
                    "the merged tissue of the standard organ."
                )
            k = self._input(
                ids.kslow(s.sid),
                U.per_min,
                f"exchange of the slow pool {s.sid} per litre of tissue [1/min]",
            )
            kp = self._input(
                ids.kp_slow(s.sid),
                U.dimensionless,
                f"slow pool to tissue coefficient {s.sid}",
            )
            vti = ids.compartment(Part.TISSUE)
            cti = ids.species(Part.TISSUE, s.sid)
            cslow = ids.slow_species(s.sid)
            self.species[cslow] = self._annotate(
                Species(
                    cslow,
                    vti,
                    initialConcentration=0.0,
                    substanceUnit=U.mmole,
                    hasOnlySubstanceUnits=False,
                    name=f"{s.name} (slow pool, {label(self.organ)})",
                    sboTerm=SBO.SIMPLE_CHEMICAL,
                    notes="The slow (deep) pool of the tissue, a concentration in "
                    f"the tissue volume: at equilibrium `{cslow} = {kp}*{cti}`.",
                    port=True,
                ),
                "concentration_tissue",
                substance=s,
            )
            tissue_kp = self._kp(s)
            self._defer_start(cslow, s, coefficient=f"{kp} * {tissue_kp}")
            self._reaction(
                ids.slow_reaction(s.sid),
                f"{cti} <-> {cslow}",
                f"{k} * {vti} * ({kp} * {cti} - {cslow})",
                f"exchange {s.sid} tissue and slow pool",
                kind="exchange",
            )

    def _urine_rates(self) -> None:
        """The urinary excretion rate `Rurine_<s>` of every substance which reaches the urine.

        A substance with urinary excretion or an apical transport of the
        kidney: the net flux of the reactions of the kidney into the urine
        `Aurine_<s>`, the sum of the fluxes into it (the filtration and the
        apical efflux of the standard kidney, the outflow of the collecting
        duct of the nephron kidney) less those out of it (a custom reaction
        which returns a part of the filtrate), an output which replaces
        `Rurine_<s>` of the body.
        """
        apical = {
            p.substance
            for p in self.ctx.chemistry
            if isinstance(p, Transport)
            and p.organ is self.organ
            and p.membrane is Membrane.APICAL
        }
        for s in self.ctx.substances:
            if not (s.excretion.urine or s.sid in apical):
                continue
            urine = ids.urine(s.sid)
            fluxes = [
                str(r.sid)
                for r in self.reactions
                if any(p.species == urine for p in r.equation.products)
            ]
            if not fluxes:
                continue
            returned = [
                str(r.sid)
                for r in self.reactions
                if any(p.species == urine for p in r.equation.reactants)
            ]
            pid = ids.urine_rate(s.sid)
            self._rule(
                pid,
                " + ".join(fluxes) + "".join(f" - {r}" for r in returned),
                U.mmole_per_min,
                f"urinary excretion rate {s.sid} [mmole/min]",
                "The rate of the urinary excretion: the sum of the fluxes of the "
                "kidney into the urine"
                + (" less those out of it." if returned else "."),
                port=Port(f"{pid}_port", portType=PortType.OUTPUT_PORT),
                kind="urine_rate",
            )
            self.outputs[pid] = pid

    def _processes(self) -> list[Metabolism | Transport | CustomReaction]:
        """The processes of the chemistry of the organ."""
        return [p for p in self.ctx.chemistry if p.organ is self.organ]

    def _proteins_and_sinks(
        self, processes: Sequence[Metabolism | Transport | CustomReaction]
    ) -> None:
        """The activities of the proteins of the processes and the sinks of their apical transport.

        An apical influx reads the inhibitors of its protein from the sink.
        """
        for p in processes:
            if isinstance(p, Metabolism | Transport) and p.protein is not None:
                self._activity(p.protein)
            if isinstance(p, Transport) and p.membrane is Membrane.APICAL:
                self._sink(p.substance)
                if p.direction is Direction.INFLUX:
                    for sid, inh in self.ctx.inhibitions:
                        if inh.protein == p.protein:
                            self._sink(sid)

    def _sink(self, s: str) -> None:
        """The apical sink species of a substance in the organ, an element of the body.

        Read by an apical transport of the substance and, for an apical
        influx, by the inhibitors of its protein, which act from the sink.
        """
        sink = apical_sink(self.organ, s)
        self._input_compartment(sink.compartment, f"sink {sink.compartment}")
        self._input_species(
            sink.species,
            sink.compartment,
            f"{s} ({sink.compartment})",
            amount=sink.amount,
        )
