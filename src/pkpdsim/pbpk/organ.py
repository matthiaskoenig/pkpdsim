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

The volumes of the parts, the coefficients and the exchanges depend on
parameters only and are constants with initial assignments (they follow a
change of the body before a simulation, not during it); only the outflow
concentration, which reads a species, is an assignment rule.
"""

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Protocol

from sbmlutils.factory import (
    AssignmentRule,
    Compartment,
    InitialAssignment,
    ModelDefinition,
    Parameter,
    Port,
    PortType,
    Reaction,
    Species,
    UnitDefinition,
)
from sbmlutils.metadata import SBO

from pkpdsim.pbpk import distribution, ids
from pkpdsim.pbpk.chemistry import CustomReaction, Membrane, Metabolism, Transport
from pkpdsim.pbpk.organs import Organ, Part, Pool, inflow, outflow
from pkpdsim.pbpk.parts import Parts
from pkpdsim.pbpk.placement import (
    apical_sink,
    built_parts,
    kp_id,
    place,
    species_compartment,
    value_parameter,
)
from pkpdsim.pbpk.proteins import Induction, Inhibition
from pkpdsim.pbpk.substance import Substance
from pkpdsim.pbpk.units import MODEL_UNITS, U
from pkpdsim.pbpk.values import Q

__all__ = ["Boundary", "OrganComponent", "OrganContext", "StandardOrgan"]


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
    """

    organ: Organ
    substances: tuple[Substance, ...]
    parts: Mapping[str, Parts]
    chemistry: tuple[Metabolism | Transport | CustomReaction, ...]
    inhibitions: tuple[tuple[str, Inhibition], ...]
    inductions: tuple[tuple[str, Induction], ...]


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
        return _Builder(ctx, self.name).build()


_INPUT = (
    "Input of the organ, replaced by the element of the body; the value is a "
    "placeholder."
)
"""The notes of an element of the body in the organ."""

_ONE = "1 dimensionless"

_POOL_PART = {Part.PLASMA: "plasma", Part.BLOOD_CELLS: "blood cells"}
"""The parts of a pool in names."""

_CAPACITY = (
    "Not a physical volume: the plasma-equivalent capacity of the blood of a "
    "substance with erythrocytes in equilibrium, the vascular volume times the "
    "blood to plasma ratio, so that the plasma concentration in it gives the "
    "amount in the blood."
)
"""The notes of a capacity of the blood."""

_FTISSUE = (
    "ftissue multiplies `fu*(Cpl - Cti/Kp)` (the unbound plasma concentration "
    "and the tissue concentration as plasma equivalent), so its fitted values "
    "differ from those of the hand-written models, where it multiplied the "
    "total concentrations."
)
"""The meaning of ftissue in the notes."""


def _label(organ: Organ) -> str:
    """The organ in a name, `small intestine`."""
    return organ.name.lower().replace("_", " ")


def _flows(organ: Organ) -> tuple[list[tuple[Pool, str]], str]:
    """The inflows (pool, flow) of an organ and the flow of its outflow."""
    if organ is Organ.LIVER:
        return [(Pool.ARTERIAL, ids.QHA), (Pool.PORTAL, ids.QPO)], ids.QH
    flow = ids.organ_flow(organ)
    return [(pool, flow) for pool in inflow(organ)], flow


class _Builder:
    """The SBML objects of one standard organ, collected by id."""

    def __init__(self, ctx: OrganContext, name: str) -> None:
        self.ctx = ctx
        self.organ = ctx.organ
        self.name = name
        self.compartments: dict[str, Compartment] = {}
        self.species: dict[str, Species] = {}
        self.parameters: dict[str, Parameter] = {}
        self.rules: list[AssignmentRule] = []
        self.assignments: list[InitialAssignment] = []
        self.reactions: list[Reaction] = []
        self.replacements: dict[str, str] = {}
        self.outputs: dict[str, str] = {}

    def build(self) -> tuple[ModelDefinition, Boundary]:
        """Build the model definition and the boundary."""
        self._compartments()
        for s in self.ctx.substances:
            self._substance(s)
        self._chemistry()
        label = _label(self.organ)
        md = ModelDefinition(
            sid=f"{self.name}_{self.organ.value}",
            name=f"{label} ({self.name} organ)",
            notes=f"The {self.name} organ {label}: the parts which the substances "
            "need, their blood flows, exchanges and chemistry.",
            units=U,
            model_units=MODEL_UNITS,
            compartments=list(self.compartments.values()),
            species=list(self.species.values()),
            parameters=list(self.parameters.values()),
            rules=self.rules,
            assignments=self.assignments,
            reactions=self.reactions,
        )
        boundary = Boundary(
            replacements=dict(self.replacements), outputs=dict(self.outputs)
        )
        return md, boundary

    # --- elements --------------------------------------------------------

    def _input(self, pid: str, unit: UnitDefinition, name: str, text: str = "") -> str:
        """A parameter of the body in the organ, with an input port."""
        if pid not in self.parameters:
            self.parameters[pid] = Parameter(
                pid,
                1.0,
                unit,
                constant=True,
                name=name,
                sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                notes=f"{text}\n\n{_INPUT}" if text else _INPUT,
                port=Port(f"{pid}_port", portType=PortType.INPUT_PORT),
            )
            self.replacements[pid] = pid
        return pid

    def _input_compartment(
        self, cid: str, name: str, text: str = "", sbo: str = SBO.PHYSICAL_COMPARTMENT
    ) -> str:
        """A compartment of the body in the organ, with an input port."""
        if cid not in self.compartments:
            self.compartments[cid] = Compartment(
                cid,
                1.0,
                U.litre,
                constant=True,
                name=name,
                sboTerm=sbo,
                notes=f"{text}\n\n{_INPUT}" if text else _INPUT,
                port=Port(f"{cid}_port", portType=PortType.INPUT_PORT),
            )
            self.replacements[cid] = cid
        return cid

    def _input_species(
        self, sid: str, compartment: str, name: str, amount: bool = False
    ) -> str:
        """A species of the body (pool or sink) which the organ reads and writes."""
        if sid not in self.species:
            self.species[sid] = Species(
                sid,
                compartment,
                initialAmount=0.0 if amount else None,
                initialConcentration=None if amount else 0.0,
                substanceUnit=U.mmole,
                hasOnlySubstanceUnits=amount,
                name=name,
                sboTerm=SBO.SIMPLE_CHEMICAL,
                notes=_INPUT,
                port=True,
            )
            self.replacements[sid] = sid
        return sid

    def _rule(
        self,
        pid: str,
        formula: str,
        unit: UnitDefinition,
        name: str,
        text: str,
        port: bool | Port = False,
    ) -> str:
        """A parameter of the organ with its assignment rule (it reads species)."""
        self.parameters[pid] = Parameter(
            pid,
            None,
            unit,
            constant=False,
            name=name,
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes=f"{text}\n\n`{pid} = {formula}`",
            port=port,
        )
        self.rules.append(AssignmentRule(pid, formula, unit))
        return pid

    def _initial(
        self,
        pid: str,
        formula: str,
        unit: UnitDefinition,
        name: str,
        text: str,
    ) -> str:
        """A constant parameter of the organ with its initial assignment.

        For a quantity of the parameters of the body (volumes, fractions,
        substance parameters), which changes with them before a simulation.
        """
        self.parameters[pid] = Parameter(
            pid,
            None,
            unit,
            constant=True,
            name=name,
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes=f"{text}\n\n`{pid} = {formula}`",
        )
        self.assignments.append(InitialAssignment(pid, formula, unit))
        return pid

    def _constant(
        self, pid: str, q: Q, unit: UnitDefinition, name: str, text: str = ""
    ) -> str:
        """A constant parameter of the organ with its provenance."""
        if pid not in self.parameters:
            self.parameters[pid] = value_parameter(
                pid, q, unit, name, SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER, text
            )
        return pid

    def _compartment(
        self,
        cid: str,
        formula: str,
        name: str,
        text: str = "",
        sbo: str = SBO.PHYSICAL_COMPARTMENT,
    ) -> None:
        """A constant compartment of the organ with its volume as initial assignment."""
        if cid in self.compartments:
            return
        self.compartments[cid] = Compartment(
            cid,
            None,
            U.litre,
            constant=True,
            name=f"{name} {_label(self.organ)}",
            sboTerm=sbo,
            notes=f"{text}\n\n`{cid} = {formula}`" if text else f"`{cid} = {formula}`",
            port=True,
        )
        self.assignments.append(InitialAssignment(cid, formula, U.litre))

    def _species(self, part: Part, s: Substance, parts: Parts) -> None:
        """The species of a substance in a part of the organ."""
        sid = ids.species(part, s.sid)
        self.species[sid] = Species(
            sid,
            species_compartment(part, s.sid, parts),
            initialConcentration=0.0,
            substanceUnit=U.mmole,
            hasOnlySubstanceUnits=False,
            name=f"{s.name} ({part.name.lower().replace('_', ' ')}, "
            f"{_label(self.organ)})",
            sboTerm=SBO.SIMPLE_CHEMICAL,
            port=True,
        )

    def _reaction(self, sid: str, equation: str, rate: str, name: str) -> None:
        """A transport of the organ, the rate in mmole/min."""
        self.reactions.append(
            Reaction(
                sid=sid,
                equation=equation,
                formula=(rate, U.mmole_per_min),
                name=name,
                sboTerm=SBO.TRANSPORT_REACTION,
            )
        )

    # --- physiology ------------------------------------------------------

    @property
    def v(self) -> str:
        """The organ volume."""
        return self._input(
            ids.organ_volume(self.organ), U.litre, f"volume {_label(self.organ)} [l]"
        )

    @property
    def hct(self) -> str:
        """The hematocrit."""
        return self._input(ids.HCT, U.dimensionless, "hematocrit")

    @property
    def f_vas(self) -> str:
        """The vascular fraction of the organ volume."""
        return self._input(
            ids.vascular_fraction(self.organ),
            U.dimensionless,
            f"vascular fraction {_label(self.organ)}",
        )

    @property
    def f_int(self) -> str:
        """The interstitial fraction of the organ volume."""
        return self._input(
            ids.interstitial_fraction(self.organ),
            U.dimensionless,
            f"interstitial fraction {_label(self.organ)}",
        )

    def _flow(self, fid: str) -> str:
        """A blood flow."""
        return self._input(fid, U.l_per_min, f"blood flow {fid} [l/min]")

    def _compartments(self) -> None:
        """The compartments of the parts which hold a species of any substance."""
        built = {p for parts in self.ctx.parts.values() for p in built_parts(parts)}
        v, f_vas = self.v, self.f_vas
        if Part.PLASMA in built:
            self._compartment(
                ids.compartment(Part.PLASMA),
                f"{v} * {f_vas} * ({_ONE} - {self.hct})",
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
                f"{v} * {self.f_int}",
                "interstitium",
            )
            self._compartment(
                ids.compartment(Part.INTRACELLULAR),
                f"{v} * ({_ONE} - {f_vas} - {self.f_int})",
                "intracellular space",
            )
        if Part.TISSUE in built:
            self._compartment(
                ids.compartment(Part.TISSUE), f"{v} * ({_ONE} - {f_vas})", "tissue"
            )
        if Part.ORGAN in built:
            self._compartment(ids.compartment(Part.ORGAN), v, "organ")

    # --- substances ------------------------------------------------------

    def _substance(self, s: Substance) -> None:
        """The species, outflow concentration, flows and exchanges of a substance."""
        parts = self.ctx.parts[s.sid]
        if parts.tissue == "organ" and parts.blood_cells:
            raise ValueError(
                f"{s.sid}: kinetic erythrocytes need a distribution with exchange, "
                f"not the flow-limited organ '{self.organ.value}'."
            )
        fu = self._input(ids.fu(s.sid), U.dimensionless, f"fraction unbound {s.sid}")
        if parts.flow_carrier == "blood_bp" and parts.tissue != "organ":
            self._compartment(
                ids.plasma_capacity(s.sid),
                f"{self.v} * {self.f_vas} * {self._bp(s)}",
                f"plasma-equivalent capacity V*f_vas*BP of {s.sid},",
                _CAPACITY,
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
            f"outflow plasma concentration {s.sid} {_label(self.organ)} [mM]",
            "Plasma concentration of the blood leaving the organ.",
            port=Port(f"{out}_port", portType=PortType.OUTPUT_PORT),
        )
        self.outputs[out] = ids.organ_outflow(self.organ, s.sid)
        self._blood_flows(s, parts)
        if parts.blood_cells:
            self._blood_cells(s)
        if parts.tissue == "merged":
            self._merged(s, fu)
        elif parts.tissue == "full":
            self._full(s, fu)

    def _bp(self, s: Substance) -> str:
        """The blood to plasma ratio of a substance."""
        return self._input(
            ids.bp(s.sid), U.dimensionless, f"blood to plasma ratio {s.sid}"
        )

    def _kp(self, s: Substance) -> str | None:
        """The whole-tissue to plasma coefficient of a substance in the organ."""
        kp = kp_id(s, self.organ)
        if kp is None:
            return None
        return self._input(kp, U.dimensionless, f"tissue to plasma coefficient {s.sid}")

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
            vascular = f"({_ONE} - {self.hct})"
        return self._initial(
            ids.k_organ(s.sid),
            distribution.k_organ(self.f_vas, vascular, kp),
            U.dimensionless,
            f"organ to plasma coefficient {s.sid}",
            "Organ to plasma coefficient of the flow-limited organ, "
            "`(Vpl + Vbc*Kbc + Vti*Kp)/V`: the organ is in equilibrium with the "
            "plasma leaving it.",
        )

    def _pool_species(self, pool: Pool, part: Part, s: Substance, parts: Parts) -> str:
        """The species of a substance in a pool, an element of the body."""
        if part is Part.PLASMA and parts.flow_carrier == "blood_bp":
            compartment = ids.pool_capacity(pool, s.sid)
            self._input_compartment(
                compartment,
                f"plasma-equivalent capacity V*BP of {s.sid}, {pool.name.lower()} pool",
                _CAPACITY,
                SBO.IMPLICIT_COMPARTMENT,
            )
        else:
            compartment = ids.pool_compartment(pool, part)
            self._input_compartment(
                compartment, f"{pool.name.lower()} {_POOL_PART[part]}"
            )
        return self._input_species(
            ids.pool_species(pool, s.sid, part),
            compartment,
            f"{s.name} ({pool.name.lower()} {_POOL_PART[part]})",
        )

    def _blood_flows(self, s: Substance, parts: Parts) -> None:
        """The inflow from the pools and the outflow into the pool, plasma and blood cells.

        The plasma carries `Q*(1-HCT)*C`, or `Q*BP*C` with erythrocytes in
        equilibrium; kinetic blood cells carry `Q*HCT*Cbc`.
        """
        inflows, out_flow = _flows(self.organ)
        out_pool = outflow(self.organ)
        if parts.flow_carrier == "blood_bp":
            carrier = self._bp(s)
        else:
            carrier = f"({_ONE} - {self.hct})"
        target = Part.ORGAN if parts.tissue == "organ" else Part.PLASMA
        organ_species = ids.species(target, s.sid)
        for pool, flow in inflows:
            c = self._pool_species(pool, Part.PLASMA, s, parts)
            self._reaction(
                ids.inflow_reaction(pool, Part.PLASMA, s.sid),
                f"{c} -> {organ_species}",
                f"{self._flow(flow)} * {carrier} * {c}",
                f"inflow {s.sid} from {pool.name.lower()} pool",
            )
        c = self._pool_species(out_pool, Part.PLASMA, s, parts)
        self._reaction(
            ids.outflow_reaction(Part.PLASMA, s.sid),
            f"{organ_species} -> {c}",
            f"{self._flow(out_flow)} * {carrier} * {ids.outflow_conc(s.sid)}",
            f"outflow {s.sid} into {out_pool.name.lower()} pool",
        )
        if not parts.blood_cells:
            return
        cbc = ids.species(Part.BLOOD_CELLS, s.sid)
        for pool, flow in inflows:
            c = self._pool_species(pool, Part.BLOOD_CELLS, s, parts)
            self._reaction(
                ids.inflow_reaction(pool, Part.BLOOD_CELLS, s.sid),
                f"{c} -> {cbc}",
                f"{self._flow(flow)} * {self.hct} * {c}",
                f"inflow {s.sid} blood cells from {pool.name.lower()} pool",
            )
        c = self._pool_species(out_pool, Part.BLOOD_CELLS, s, parts)
        self._reaction(
            ids.outflow_reaction(Part.BLOOD_CELLS, s.sid),
            f"{cbc} -> {c}",
            f"{self._flow(out_flow)} * {self.hct} * {cbc}",
            f"outflow {s.sid} blood cells into {out_pool.name.lower()} pool",
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
            f"{ftissue} * {self.v} * ({_ONE} - {self.f_vas}) / {total}",
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

    # --- chemistry -------------------------------------------------------

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
                sink = apical_sink(self.organ, p.substance)
                self._input_compartment(sink.compartment, f"sink {sink.compartment}")
                self._input_species(
                    sink.species,
                    sink.compartment,
                    f"{p.substance} ({sink.compartment})",
                    amount=sink.amount,
                )
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
        for p in processes:
            self.reactions += place(p, self.ctx, self.ctx.parts)

    def _activity(self, protein: str) -> None:
        """The activity `f_<protein>` of a protein in the organ."""
        pid = ids.activity(protein)
        if pid in self.parameters:
            return
        self.parameters[pid] = Parameter(
            pid,
            1.0,
            U.dimensionless,
            constant=True,
            name=f"activity {protein}",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes=f"Activity of {protein} in the organ relative to normal (1.0); "
            "it multiplies the rates of its processes.",
            port=True,
        )
