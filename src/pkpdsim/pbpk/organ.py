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
The kidney filters the unbound plasma concentration of a substance with
urinary excretion into the urine, `fu*GFR*Cpl*(1 - f_reabsorbed)`.

The volumes of the parts, the coefficients and the exchanges depend on
parameters only and are constants with initial assignments (they follow a
change of the body before a simulation, not during it); only what reads a
species is an assignment rule: the outflow concentrations and the unbound
concentrations of the tissue parts `Cu_<part>_<s>` (observers).

The elements (inputs with their ports, compartments, species, rules,
constants) are collected by `pkpdsim.pbpk.organ_builder.OrganBuilder`, which
the variants reuse.
"""

from collections.abc import Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING, Protocol, runtime_checkable

from sbmlutils.factory import ModelDefinition, Port, PortType, Sbase
from sbmlutils.metadata import SBO

from pkpdsim.pbpk import distribution, ids
from pkpdsim.pbpk.chemistry import (
    CustomReaction,
    Direction,
    Membrane,
    Metabolism,
    Transport,
)
from pkpdsim.pbpk.organ_builder import CAPACITY, ONE, OrganBuilder, label
from pkpdsim.pbpk.organs import Organ, Part, Pool, inflow, outflow
from pkpdsim.pbpk.parts import Parts
from pkpdsim.pbpk.placement import apical_sink, built_parts, kp_id, place, unbound
from pkpdsim.pbpk.proteins import Induction, Inhibition
from pkpdsim.pbpk.routes import Route
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
        routes: the routes of administration (the intestine applies the
            oral doses).
    """

    organ: Organ
    substances: tuple[Substance, ...]
    parts: Mapping[str, Parts]
    chemistry: tuple[Metabolism | Transport | CustomReaction, ...]
    inhibitions: tuple[tuple[str, Inhibition], ...]
    inductions: tuple[tuple[str, Induction], ...]
    routes: tuple[Route, ...] = ()


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

    The intestine variants derive from it for the intestinal wall.
    """

    def build(self) -> tuple[ModelDefinition, Boundary]:
        """Build the model definition and the boundary."""
        self._compartments()
        for s in self.ctx.substances:
            self._substance(s)
        self._chemistry()
        md = self.model_definition(
            f"The {self.name} organ {label(self.organ)}: the parts which the "
            "substances need, their blood flows, exchanges and chemistry."
        )
        boundary = Boundary(
            replacements=dict(self.replacements), outputs=dict(self.outputs)
        )
        return md, boundary

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
        """The compartments of the parts which hold a species of any substance."""
        built = {p for parts in self.ctx.parts.values() for p in built_parts(parts)}
        v, f_vas = self.v, self.f_vas
        if Part.PLASMA in built:
            self._compartment(
                ids.compartment(Part.PLASMA),
                f"{v} * {f_vas} * ({ONE} - {self.hct})",
                "plasma",
            )
        if Part.BLOOD_CELLS in built:
            self._compartment(
                ids.compartment(Part.BLOOD_CELLS),
                f"{v} * {f_vas} * {self.hct}",
                "blood cells",
            )
        if Part.INTERSTITIAL in built:
            self._compartment(
                ids.compartment(Part.INTERSTITIAL),
                self._functional(f"{v} * {self.f_int}"),
                "interstitium",
                self._functional_text(),
            )
            self._compartment(
                ids.compartment(Part.INTRACELLULAR),
                self._functional(f"{v} * ({ONE} - {f_vas} - {self.f_int})"),
                "intracellular space",
                self._functional_text(),
            )
        if Part.TISSUE in built:
            self._compartment(
                ids.compartment(Part.TISSUE),
                self._functional(f"{v} * ({ONE} - {f_vas})"),
                "tissue",
                self._functional_text(),
            )
        if Part.ORGAN in built:
            self._compartment(ids.compartment(Part.ORGAN), v, "organ")

    def _substance(self, s: Substance) -> None:
        """The species, outflow concentration, flows and exchanges of a substance."""
        parts = self.ctx.parts[s.sid]
        if parts.tissue == "organ" and parts.blood_cells:
            raise ValueError(
                f"{s.sid}: kinetic erythrocytes need a distribution with exchange, "
                f"not the flow-limited organ '{self.organ.value}'."
            )
        fu = self._fu(s)
        if parts.flow_carrier == "blood_bp" and parts.tissue != "organ":
            self._compartment(
                ids.plasma_capacity(s.sid),
                f"{self.v} * {self.f_vas} * {self._bp(s)}",
                f"plasma-equivalent capacity V*f_vas*BP of {s.sid},",
                CAPACITY,
                SBO.IMPLICIT_COMPARTMENT,
            )
        for part in built_parts(parts):
            self._species(part, s, parts)
        if parts.tissue == "organ":
            cout = f"{ids.species(Part.ORGAN, s.sid)} / {self._k_organ(s, parts)}"
        else:
            cout = ids.species(Part.PLASMA, s.sid)
        out = ids.outflow_conc(s.sid)
        self._rule(
            out,
            cout,
            U.mM,
            f"outflow plasma concentration {s.sid} {label(self.organ)} [mM]",
            "Plasma concentration of the blood leaving the organ.",
            port=Port(f"{out}_port", portType=PortType.OUTPUT_PORT),
        )
        self.outputs[out] = ids.organ_outflow(self.organ, s.sid)
        if parts.blood_cells:
            out_bc = ids.outflow_conc_bc(s.sid)
            self._rule(
                out_bc,
                ids.species(Part.BLOOD_CELLS, s.sid),
                U.mM,
                f"outflow blood cell concentration {s.sid} {label(self.organ)} [mM]",
                "Blood cell concentration of the blood leaving the organ.",
                port=Port(f"{out_bc}_port", portType=PortType.OUTPUT_PORT),
            )
            self.outputs[out_bc] = ids.organ_outflow_bc(self.organ, s.sid)
        self._blood_flows(s, parts)
        if parts.blood_cells:
            self._blood_cells(s)
        if self.organ is Organ.KIDNEY and s.excretion.urine:
            self._filtration(s, parts, fu, cout)
        if parts.tissue == "merged":
            self._merged(s, fu)
        elif parts.tissue == "full":
            self._full(s, fu)
        self._unbound_observers(s, parts)

    def _unbound_observers(self, s: Substance, parts: Parts) -> None:
        """The unbound concentrations of the tissue parts of a substance, `Cu_<part>_<s>`.

        `fu/K*C` with the coefficient of the part to plasma (`Kp`,
        `K_int_pls`, `K_cell_pls`, `K_or`): observers, assignment rules of
        the organ (they read species), the concentration which the processes
        of the part act on.
        """
        for part in built_parts(parts):
            if part not in _TISSUE_PARTS:
                continue
            self._rule(
                ids.unbound(part, s.sid),
                unbound(part, s, self.organ),
                U.mM,
                f"unbound concentration {s.sid} {_TISSUE_PARTS[part]} "
                f"{label(self.organ)} [mM]",
                f"Unbound concentration of the {_TISSUE_PARTS[part]}; an observer.",
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
        )

    def _blood_flows(self, s: Substance, parts: Parts) -> None:
        """The inflow from the pools and the outflow into the pool, plasma and blood cells.

        The plasma carries `Q*(1-HCT)*C`, or `Q*BP*C` with erythrocytes in
        equilibrium; kinetic blood cells carry `Q*HCT*Cbc`. The fraction
        `f_shunt` of the inflow of the liver bypasses it (portosystemic shunt)
        from the arterial and the portal pool directly into the venous pool;
        the liver is perfused by and drains `(1 - f_shunt)*Q`.
        """
        inflows, out_flow = self._flows(s)
        out_pool = outflow(self.organ)
        if parts.flow_carrier == "blood_bp":
            carrier = self._bp(s)
        else:
            carrier = f"({ONE} - {self.hct})"
        target = Part.ORGAN if parts.tissue == "organ" else Part.PLASMA
        carried = [(Part.PLASMA, carrier, ids.species(target, s.sid))]
        if parts.blood_cells:
            carried.append(
                (Part.BLOOD_CELLS, self.hct, ids.species(Part.BLOOD_CELLS, s.sid))
            )
        perfused, f_shunt = "", None
        if self.organ is Organ.LIVER:
            f_shunt = self._input(
                ids.F_SHUNT,
                U.dimensionless,
                "portosystemic shunt",
                "The fraction of the liver inflow which bypasses the liver.",
            )
            perfused = f"({ONE} - {f_shunt}) * "
        for part, share, local in carried:
            what = "" if part is Part.PLASMA else " blood cells"
            out = self._pool_species(out_pool, part, s, parts)
            for pool, flow in inflows:
                c = self._pool_species(pool, part, s, parts)
                rate = f"{self._flow(flow)} * {share} * {c}"
                self._reaction(
                    ids.inflow_reaction(pool, part, s.sid),
                    f"{c} -> {local}",
                    f"{perfused}{rate}",
                    f"inflow {s.sid}{what} from {pool.name.lower()} pool",
                )
                if f_shunt is not None:
                    self._reaction(
                        ids.shunt_reaction(pool, part, s.sid),
                        f"{c} -> {out}",
                        f"{f_shunt} * {rate}",
                        f"shunt {s.sid}{what} from {pool.name.lower()} pool past "
                        f"the liver into {out_pool.name.lower()} pool",
                    )
            leaving = (
                ids.outflow_conc(s.sid)
                if part is Part.PLASMA
                else ids.outflow_conc_bc(s.sid)
            )
            self._reaction(
                ids.outflow_reaction(part, s.sid),
                f"{local} -> {out}",
                f"{perfused}{self._flow(out_flow)} * {share} * {leaving}",
                f"outflow {s.sid}{what} into {out_pool.name.lower()} pool",
            )

    def _filtration(self, s: Substance, parts: Parts, fu: str, cplasma: str) -> None:
        """The glomerular filtration into the urine, `fu*GFR*Cpl*(1 - f_reabsorbed)`.

        The unbound plasma concentration is filtered, the fraction
        `f_reabsorbed_<s>` of the filtrate returns into the plasma (it is not
        taken out); the flow-limited kidney filters the plasma concentration of
        its outflow.
        """
        gfr = self._input(ids.GFR, U.l_per_min, "glomerular filtration rate [l/min]")
        reabsorbed = self._input(
            ids.fraction_reabsorbed(s.sid),
            U.dimensionless,
            f"fraction reabsorbed {s.sid}",
            "The fraction of the filtrate which is reabsorbed.",
        )
        self._input_compartment(ids.VURINE, "urine")
        urine = self._input_species(
            ids.urine(s.sid), ids.VURINE, f"{s.name} (urine)", amount=True
        )
        source = ids.species(
            Part.ORGAN if parts.tissue == "organ" else Part.PLASMA, s.sid
        )
        self._reaction(
            ids.filtration_reaction(s.sid),
            f"{source} -> {urine}",
            f"{fu} * {gfr} * {cplasma} * ({ONE} - {reabsorbed})",
            f"glomerular filtration {s.sid}",
        )

    def _blood_cells(self, s: Substance) -> None:
        """The exchange between plasma and kinetic blood cells.

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
        cpl, cbc = ids.species(Part.PLASMA, s.sid), ids.species(Part.BLOOD_CELLS, s.sid)
        self._reaction(
            ids.exchange_reaction(Part.BLOOD_CELLS, s.sid),
            f"{cpl} <-> {cbc}",
            f"{kex} * {ids.compartment(Part.BLOOD_CELLS)} * ({cpl} - {cbc} / {kbc})",
            f"exchange {s.sid} plasma and blood cells",
        )

    def _ps(self, s: Substance) -> str:
        """The passive exchange between plasma and tissue of a substance in the organ.

        `PS = ftissue * Vti/Vti_total_ref`, the share of the whole-body exchange
        by the tissue volume relative to the tissue volume of the reference
        individual (so the exchange scales with the individual), or the value
        of the organ (`ps` of the Kp); no
        passive exchange for a substance which does not enter the tissue of
        the organ (it has the full form by its transporters).

        Raises:
            ValueError: if the substance enters the tissue but has neither
                `ftissue` nor a value of the organ.
        """
        pid = ids.ps(s.sid)
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
            )
        if self.organ in dist.ps:
            return self._constant(
                pid,
                dist.ps[self.organ],
                U.l_per_min,
                name,
                "Passive exchange of the organ, given instead of the share of the "
                "whole-body exchange.",
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
        )

    def _merged(self, s: Substance, fu: str) -> None:
        """The exchange between plasma and the merged tissue, `PS*fu*(Cpl - Cti/Kp)`."""
        kp = self._kp(s)
        assert kp is not None
        ps = self._ps(s)
        cpl, cti = ids.species(Part.PLASMA, s.sid), ids.species(Part.TISSUE, s.sid)
        self._reaction(
            ids.exchange_reaction(Part.TISSUE, s.sid),
            f"{cpl} <-> {cti}",
            f"{ps} * {fu} * ({cpl} - {cti} / {kp})",
            f"exchange {s.sid} plasma and tissue",
        )

    def _full(self, s: Substance, fu: str) -> None:
        """The exchanges plasma - interstitium - intracellular space.

        `PS_vas*fu*(Cpl - Cis/K_int_pls)` and
        `PS*(fu/K_int_pls*Cis - fu/K_cell_pls*Cic)`, with `K_int_pls` of
        Schmitt 2008 and `K_cell_pls` such that the whole tissue has the Kp.
        """
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
            )
        else:
            distribution.check_kp(s, self.organ)
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
            )
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
            )
        ps = self._ps(s)
        cpl = ids.species(Part.PLASMA, s.sid)
        cis = ids.species(Part.INTERSTITIAL, s.sid)
        cic = ids.species(Part.INTRACELLULAR, s.sid)
        self._reaction(
            ids.exchange_reaction(Part.INTERSTITIAL, s.sid),
            f"{cpl} <-> {cis}",
            f"{ids.PS_VAS} * {fu} * ({cpl} - {cis} / {k_int})",
            f"exchange {s.sid} plasma and interstitium",
        )
        self._reaction(
            ids.exchange_reaction(Part.INTRACELLULAR, s.sid),
            f"{cis} <-> {cic}",
            f"{ps} * ({fu} / {k_int} * {cis} - {fu} / {k_cell} * {cic})",
            f"passive exchange {s.sid} interstitium and cells",
        )

    def _chemistry(self) -> None:
        """The processes of the organ with the activities of their proteins and their sinks.

        Raises:
            ValueError: if a species of a custom reaction has the id of a
                species of the organ.
        """
        processes = [p for p in self.ctx.chemistry if p.organ is self.organ]
        for p in processes:
            if isinstance(p, Metabolism | Transport) and p.protein is not None:
                self._activity(p.protein)
            if isinstance(p, Transport) and p.membrane is Membrane.APICAL:
                self._sink(p.substance)
                if p.direction is Direction.INFLUX:
                    for sid, inh in self.ctx.inhibitions:
                        if inh.protein == p.protein:
                            self._sink(sid)
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
                self.species[sid] = sp
        capacity = self._capacity()
        for p in processes:
            self.reactions += place(p, self.ctx, self.ctx.parts, capacity)

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
