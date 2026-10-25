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
"""

from collections.abc import Mapping
from dataclasses import dataclass
from functools import cache
from typing import Protocol

from sbmlutils.factory import (
    AssignmentRule,
    Compartment,
    ModelDefinition,
    Parameter,
    Port,
    PortType,
    Reaction,
    Species,
    UnitDefinition,
)
from sbmlutils.metadata import SBO

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.chemistry import CustomReaction, Membrane, Metabolism, Transport
from pkpdsim.pbpk.data import load
from pkpdsim.pbpk.organs import Organ, Part, Pool, inflow, outflow
from pkpdsim.pbpk.parts import Parts
from pkpdsim.pbpk.physiology import reference_fractions
from pkpdsim.pbpk.placement import (
    apical_sink,
    built_parts,
    kp_id,
    place,
    species_compartment,
    value_parameter,
)
from pkpdsim.pbpk.proteins import Induction, Inhibition
from pkpdsim.pbpk.substance import Kp, Substance
from pkpdsim.pbpk.units import MODEL_UNITS, U
from pkpdsim.pbpk.values import Q, Value

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

    name: str

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


_DATA_UNITS: dict[str, UnitDefinition] = {
    "dimensionless": U.dimensionless,
    "1/min": U.per_min,
}
"""The units of the model for the units of `data/distribution.toml`."""


@cache
def _distribution() -> dict[str, dict[str, Value]]:
    """The constants of `data/distribution.toml` as values with their source."""
    return {
        group: {
            key: Value(e["value"], _DATA_UNITS[e["unit"]], e["source"], e.get("notes"))
            for key, e in table.items()
        }
        for group, table in load("distribution").items()
    }


def _number(q: Q) -> float:
    """The number of a plain float or a `Value`."""
    return q.value if isinstance(q, Value) else q


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
            reactions=self.reactions,
        )
        boundary = Boundary(
            replacements=dict(self.replacements), outputs=dict(self.outputs)
        )
        return md, boundary

    # --- elements --------------------------------------------------------

    def _input(self, pid: str, unit: UnitDefinition, name: str) -> str:
        """A parameter of the body in the organ, with an input port."""
        if pid not in self.parameters:
            self.parameters[pid] = Parameter(
                pid,
                1.0,
                unit,
                constant=False,
                name=name,
                sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                notes=_INPUT,
                port=Port(f"{pid}_port", portType=PortType.INPUT_PORT),
            )
            self.replacements[pid] = pid
        return pid

    def _input_compartment(self, cid: str, name: str) -> str:
        """A compartment of the body in the organ, with an input port."""
        if cid not in self.compartments:
            self.compartments[cid] = Compartment(
                cid,
                1.0,
                U.litre,
                constant=False,
                name=name,
                sboTerm=SBO.PHYSICAL_COMPARTMENT,
                notes=_INPUT,
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
        """A parameter of the organ with its assignment rule."""
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

    def _constant(
        self, pid: str, q: Q, unit: UnitDefinition, name: str, text: str = ""
    ) -> str:
        """A constant parameter of the organ with its provenance."""
        if pid not in self.parameters:
            self.parameters[pid] = value_parameter(
                pid, q, unit, name, SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER, text
            )
        return pid

    def _compartment(self, cid: str, formula: str, name: str) -> None:
        """A compartment of the organ with its volume as assignment rule."""
        if cid in self.compartments:
            return
        self.compartments[cid] = Compartment(
            cid,
            None,
            U.litre,
            constant=False,
            name=f"{name} {_label(self.organ)}",
            sboTerm=SBO.PHYSICAL_COMPARTMENT,
            notes=f"`{cid} = {formula}`",
            port=True,
        )
        self.rules.append(AssignmentRule(cid, formula, U.litre))

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
                f"blood (plasma volume of {s.sid} with the blood to plasma ratio)",
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
            vascular = f"{self.f_vas} * {self._bp(s)}"
        else:
            vascular = f"{self.f_vas} * ({_ONE} - {self.hct})"
        return self._rule(
            ids.k_organ(s.sid),
            f"{vascular} + ({_ONE} - {self.f_vas}) * {kp}",
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
            text = f"{pool.name.lower()} blood of {s.sid} (blood to plasma ratio)"
        else:
            compartment = ids.pool_compartment(pool, part)
            text = f"{pool.name.lower()} {_POOL_PART[part]}"
        self._input_compartment(compartment, text)
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

        `PS = ftissue * Vti/Vti_total`, the share of the whole-body exchange by
        the tissue volume, or the value of the organ (`ps` of the Kp); no
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
            ids.ftissue(s.sid), U.l_per_min, f"whole-body exchange {s.sid} [l/min]"
        )
        total = self._input(
            ids.tissue_volume_total(s.sid),
            U.litre,
            f"tissue volume of the organs {s.sid} enters [l]",
        )
        return self._rule(
            pid,
            f"{ftissue} * {self.v} * ({_ONE} - {self.f_vas}) / {total}",
            U.l_per_min,
            name,
            "Passive exchange between plasma and tissue: the share of the "
            "whole-body exchange `ftissue` by the tissue volume of the organ.",
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
        constants = _distribution()["interstitial"]
        f_water_int = self._constant(
            "f_water_int",
            constants["f_water_int"],
            U.dimensionless,
            "water fraction interstitial fluid",
        )
        f_water_pls = self._constant(
            "f_water_pls",
            constants["f_water_pls"],
            U.dimensionless,
            "water fraction plasma",
        )
        r_protein = self._constant(
            "r_protein_int_pls",
            constants["r_protein_int_pls"],
            U.dimensionless,
            "protein of the interstitial fluid relative to plasma",
        )
        k_int = self._rule(
            ids.k_int_pls(s.sid),
            f"({f_water_int} + {r_protein} * ({_ONE} / {fu} - {f_water_pls})) * {fu}",
            U.dimensionless,
            f"interstitial to plasma coefficient {s.sid}",
            "Interstitial to plasma coefficient (Schmitt 2008): the water of the "
            "interstitium and the plasma proteins of the interstitial fluid.",
        )
        kp = self._kp(s)
        if kp is None:
            k_cell = self._rule(
                ids.k_cell_pls(s.sid),
                k_int,
                U.dimensionless,
                f"intracellular to plasma coefficient {s.sid}",
                "The substance does not distribute into the tissue of the organ: "
                "without a Kp the cells bind it as the interstitial fluid.",
            )
        else:
            self._check_k_cell(s)
            f_cell = f"({_ONE} - {self.f_vas} - {self.f_int})"
            k_cell = self._rule(
                ids.k_cell_pls(s.sid),
                f"({kp} * ({self.f_int} + {f_cell}) - {self.f_int} * {k_int}) / "
                f"{f_cell}",
                U.dimensionless,
                f"intracellular to plasma coefficient {s.sid}",
                "Intracellular to plasma coefficient such that interstitium and "
                "cells together have the whole-tissue Kp, "
                "`(f_int*K_int_pls + f_cell*K_cell_pls)/(f_int + f_cell) = Kp`.",
            )
        if ids.PS_VAS not in self.parameters:
            k_vas = self._constant(
                "k_vas",
                _distribution()["vascular"]["k_vas"],
                U.per_min,
                "exchange plasma and interstitium per organ volume [1/min]",
            )
            self._rule(
                ids.PS_VAS,
                f"{k_vas} * {self.v}",
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

    def _check_k_cell(self, s: Substance) -> None:
        """Check that the Kp leaves a positive `K_cell_pls` (reference individuals).

        Raises:
            ValueError: if the Kp is not above the interstitial share of the
                whole tissue, `f_int*K_int_pls/(f_int + f_cell)`.
        """
        dist = s.distribution
        if not isinstance(dist, Kp):
            return  # the Kp of a method is computed per organ
        kp = _number(dist.overrides.get(self.organ, dist.value))
        fu = _number(s.binding.fu) if s.binding is not None else 1.0
        c = {key: v.value for key, v in _distribution()["interstitial"].items()}
        k_int = (
            c["f_water_int"] + c["r_protein_int_pls"] * (1.0 / fu - c["f_water_pls"])
        ) * fu
        for sex, (f_vas, f_int) in reference_fractions(self.organ).items():
            f_cell = 1.0 - f_vas - f_int
            share = f_int * k_int / (f_int + f_cell)
            if kp <= share:
                raise ValueError(
                    f"{s.sid}: Kp {kp} in organ '{self.organ.value}' is not above "
                    f"the interstitial share {share:.4g} of the tissue ({sex.value} "
                    f"reference individual), the full form would need "
                    f"K_cell_pls <= 0."
                )

    # --- chemistry -------------------------------------------------------

    def _chemistry(self) -> None:
        """The processes of the organ with the activities of their proteins and their sinks."""
        for p in self.ctx.chemistry:
            if p.organ is not self.organ:
                continue
            if isinstance(p, Metabolism | Transport) and p.protein is not None:
                pid = ids.activity(p.protein)
                if pid not in self.parameters:
                    self.parameters[pid] = Parameter(
                        pid,
                        1.0,
                        U.dimensionless,
                        constant=True,
                        name=f"activity {p.protein}",
                        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                        notes=f"Activity of {p.protein} in the organ relative to "
                        "normal (1.0); it multiplies the rates of its processes.",
                        port=True,
                    )
            if isinstance(p, Transport) and p.membrane is Membrane.APICAL:
                sink = apical_sink(self.organ, p.substance)
                self._input_compartment(sink.compartment, f"sink {sink.compartment}")
                self._input_species(
                    sink.species,
                    sink.compartment,
                    f"{p.substance} ({sink.compartment})",
                    amount=sink.amount,
                )
            if isinstance(p, CustomReaction):
                for sp in p.species:
                    self.species[str(sp.sid)] = sp
            self.reactions += place(p, self.ctx, self.ctx.parts)
