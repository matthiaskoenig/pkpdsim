"""The elements of an organ component: the base of the standard organ and its variants.

An organ component builds a comp model definition whose elements of the body
(volumes, flows, hematocrit, fractions, substance parameters, pool and sink
species) are inputs with a port, which the body replaces
(`replacements`), and whose outputs replace parameters of the top model
(`outputs`). `OrganBuilder` collects the elements of one organ by id and
creates each element once: the inputs, the compartments and species of the
organ, the rules (for what reads a species), the initial assignments (for
what depends on parameters only), the constants with their provenance, the
reactions and the activity of the proteins of its processes. The standard
organ and the variants derive from it and add their kinetics.

A protein of a process of the organ has the base activity `f_<protein>`, a
constant of the organ (1, set by genotype and the experiments); in the liver
the impairment multiplies it by the input `fimp_<protein>`
(`pkpdsim.pbpk.impairment`).
When a substance of the definition inactivates it (mechanism-based) or
induces it, the organ has the turnover state `E_<protein>`, the amount of the
protein relative to its baseline (initially 1),

```
dE/dt = kdeg*(1 + sum Emax*I/(EC50 + I)) - kdeg*E - sum kinact*I/(KI + I)*E
```

(Mayhew 2000, Almond 2009, Fahmi 2008), with the unbound concentration `I`
of the inhibitor or inducer in the cell (intracellular, else the equivalent
tissue part, else plasma) and the degradation rate `kdeg_<protein>` of the
organ; the rates of the protein are multiplied by `E_<protein>`. The
constants of the interactions (`Ki`, `Kiu`, `kinact`, `Emax`, `EC50`) are
parameters of the body, read through a port.
"""

from typing import TYPE_CHECKING

from sbmlutils.factory import (
    AssignmentRule,
    Compartment,
    InitialAssignment,
    ModelDefinition,
    Parameter,
    Port,
    PortType,
    RateRule,
    Reaction,
    Species,
    UnitDefinition,
)
from sbmlutils.metadata import SBO

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.annotation import annotate
from pkpdsim.pbpk.data import cite
from pkpdsim.pbpk.impairment import organ_activity
from pkpdsim.pbpk.organs import Organ, Part, Pool
from pkpdsim.pbpk.parts import Parts
from pkpdsim.pbpk.placement import (
    ONE,
    has_turnover,
    in_zone,
    inhibition_parameters,
    inhibitor_unbound,
    kp_id,
    species_compartment,
    value_parameter,
)
from pkpdsim.pbpk.proteins import turnover
from pkpdsim.pbpk.substance import Substance
from pkpdsim.pbpk.units import MODEL_UNITS, U
from pkpdsim.pbpk.values import Q

if TYPE_CHECKING:
    from pkpdsim.pbpk.organ import OrganContext

__all__ = ["CAPACITY", "COMPARTMENT_KINDS", "ONE", "OrganBuilder", "label"]

CAPACITY = (
    "Not a physical volume: the plasma-equivalent capacity of the blood of a "
    "substance with erythrocytes in equilibrium, the vascular volume times the "
    "blood to plasma ratio, so that the plasma concentration in it gives the "
    "amount in the blood."
)
"""The notes of a capacity of the blood."""

_INPUT = (
    "Input of the organ, replaced by the element of the body; the value is a "
    "placeholder."
)
"""The notes of an element of the body in the organ."""

_POOL_PART = {Part.PLASMA: "plasma", Part.BLOOD_CELLS: "blood cells"}

_SPECIES_KINDS = {
    Part.PLASMA: "concentration_plasma",
    Part.BLOOD_CELLS: "concentration",
    Part.INTERSTITIAL: "concentration_tissue",
    Part.INTRACELLULAR: "concentration_tissue",
    Part.TISSUE: "concentration_tissue",
    Part.ORGAN: "concentration_tissue",
}
"""The kind (`pkpdsim.pbpk.annotation`) of the species of a part of an organ."""

COMPARTMENT_KINDS = {
    Part.PLASMA: "plasma",
    Part.BLOOD_CELLS: "blood_cells",
    Part.INTERSTITIAL: "interstitial",
    Part.INTRACELLULAR: "intracellular",
    Part.TISSUE: "tissue",
    Part.ORGAN: "organ",
}
"""The kind (`pkpdsim.pbpk.annotation`) of the compartment of a part of an organ."""

"""The parts of a pool in names."""


def label(organ: Organ) -> str:
    """Get the organ in a name, `small intestine`."""
    return organ.name.lower().replace("_", " ")


class OrganBuilder:
    """The SBML elements of one organ, collected by id, and its boundary.

    Attributes:
        ctx: the context of the organ.
        organ: the organ.
        name: the name of the component, the prefix of the id of the model
            definition (`standard_li`).
        replacements: the elements of the body which replace an input, local
            id -> id in the body.
        outputs: the outputs of the organ, local id -> id in the body.
    """

    def __init__(self, ctx: "OrganContext", name: str) -> None:
        """Start the collection of the elements of an organ.

        Args:
            ctx: the context of the organ.
            name: the name of the component.
        """
        self.ctx = ctx
        self.organ = ctx.organ
        self.name = name
        self.compartments: dict[str, Compartment] = {}
        self.species: dict[str, Species] = {}
        self.parameters: dict[str, Parameter] = {}
        self.rules: list[AssignmentRule] = []
        self.assignments: list[InitialAssignment] = []
        self.rate_rules: list[RateRule] = []
        self.reactions: list[Reaction] = []
        self.replacements: dict[str, str] = {}
        self.outputs: dict[str, str] = {}

    def _annotate[E: (Compartment, Parameter, Reaction, Species)](
        self,
        element: E,
        kind: str,
        pool: Pool | None = None,
        segment: str | None = None,
        protein: str | None = None,
        substance: Substance | None = None,
    ) -> E:
        """Annotate an element of the organ as `kind` (`pkpdsim.pbpk.annotation`).

        An input of the organ is annotated when the organ is composed, with
        the annotation of the element of the body which replaces it.
        """
        return annotate(
            element,
            kind,
            organ=self.organ,
            pool=pool,
            segment=segment,
            protein=protein,
            substance=substance,
        )

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
        *,
        kind: str,
        pool: Pool | None = None,
        segment: str | None = None,
    ) -> str:
        """A parameter of the organ with its assignment rule (it reads species)."""
        self.parameters[pid] = self._annotate(
            Parameter(
                pid,
                None,
                unit,
                constant=False,
                name=name,
                sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                notes=f"{text}\n\n`{pid} = {formula}`",
                port=port,
            ),
            kind,
            pool=pool,
            segment=segment,
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
        *,
        kind: str,
        pool: Pool | None = None,
        segment: str | None = None,
        protein: str | None = None,
    ) -> str:
        """A constant parameter of the organ with its initial assignment.

        For a quantity of the parameters of the body (volumes, fractions,
        substance parameters), which changes with them before a simulation.
        """
        self.parameters[pid] = self._annotate(
            Parameter(
                pid,
                None,
                unit,
                constant=True,
                name=name,
                sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                notes=f"{text}\n\n`{pid} = {formula}`",
            ),
            kind,
            pool=pool,
            segment=segment,
            protein=protein,
        )
        self.assignments.append(InitialAssignment(pid, formula, unit))
        return pid

    def _constant(
        self,
        pid: str,
        q: Q,
        unit: UnitDefinition,
        name: str,
        text: str = "",
        *,
        kind: str,
        segment: str | None = None,
    ) -> str:
        """A constant parameter of the organ with its provenance."""
        if pid not in self.parameters:
            self.parameters[pid] = self._annotate(
                value_parameter(
                    pid,
                    q,
                    unit,
                    name,
                    SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                    text,
                ),
                kind,
                segment=segment,
            )
        return pid

    def _compartment(
        self,
        cid: str,
        formula: str,
        name: str,
        text: str = "",
        zone: str | None = None,
        *,
        kind: str,
        segment: str | None = None,
    ) -> None:
        """A constant compartment of the organ (of a zone) with its volume as initial assignment.

        It is annotated as `kind` (`pkpdsim.pbpk.annotation`), which gives its
        SBO term.
        """
        if cid in self.compartments:
            return
        self.compartments[cid] = self._annotate(
            Compartment(
                cid,
                None,
                U.litre,
                constant=True,
                name=in_zone(f"{name} {label(self.organ)}", zone),
                notes=f"{text}\n\n`{cid} = {formula}`"
                if text
                else f"`{cid} = {formula}`",
                port=True,
            ),
            kind,
            segment=segment,
        )
        self.assignments.append(InitialAssignment(cid, formula, U.litre))

    def _species(
        self, part: Part, s: Substance, parts: Parts, zone: str | None = None
    ) -> None:
        """The species of a substance in a part of the organ (of a zone of it)."""
        sid = ids.species(part, s.sid, zone)
        self.species[sid] = self._annotate(
            Species(
                sid,
                species_compartment(part, s.sid, parts, zone),
                initialConcentration=0.0,
                substanceUnit=U.mmole,
                hasOnlySubstanceUnits=False,
                name=f"{s.name} ({part.name.lower().replace('_', ' ')}, "
                + in_zone(label(self.organ), zone)
                + ")",
                sboTerm=SBO.SIMPLE_CHEMICAL,
                port=True,
            ),
            _SPECIES_KINDS[part],
            substance=s,
        )

    def _reaction(
        self,
        sid: str,
        equation: str,
        rate: str,
        name: str,
        *,
        kind: str,
        pool: Pool | None = None,
        segment: str | None = None,
    ) -> None:
        """A transport of the organ, the rate in mmole/min, annotated as `kind`."""
        self.reactions.append(
            self._annotate(
                Reaction(
                    sid=sid,
                    equation=equation,
                    formula=(rate, U.mmole_per_min),
                    name=name,
                    sboTerm=SBO.TRANSPORT_REACTION,
                ),
                kind,
                pool=pool,
                segment=segment,
            )
        )

    @property
    def v(self) -> str:
        """The organ volume."""
        return self._input(
            ids.organ_volume(self.organ), U.litre, f"volume {label(self.organ)} [l]"
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
            f"vascular fraction {label(self.organ)}",
        )

    @property
    def f_int(self) -> str:
        """The interstitial fraction of the organ volume."""
        return self._input(
            ids.interstitial_fraction(self.organ),
            U.dimensionless,
            f"interstitial fraction {label(self.organ)}",
        )

    def _flow(self, fid: str) -> str:
        """A blood flow."""
        return self._input(fid, U.l_per_min, f"blood flow {fid} [l/min]")

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

    def _pool_species(self, pool: Pool, part: Part, s: Substance, parts: Parts) -> str:
        """The species of a substance in a pool, an element of the body."""
        if part is Part.PLASMA and parts.flow_carrier == "blood_bp":
            compartment = ids.pool_capacity(pool, s.sid)
            self._input_compartment(
                compartment,
                f"plasma-equivalent capacity V*BP of {s.sid}, {pool.name.lower()} pool",
                CAPACITY,
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

    def _fu(self, s: Substance) -> str:
        """The fraction unbound in plasma of a substance."""
        return self._input(ids.fu(s.sid), U.dimensionless, f"fraction unbound {s.sid}")

    def _activity(self, protein: str) -> None:
        """The base activity `f_<protein>` of a protein in the organ and its interactions.

        The constants of the inhibitions and inductions of the protein are
        inputs of the organ; a protein which a substance inactivates or
        induces has the turnover state `E_<protein>`.
        """
        pid = ids.activity(protein)
        if pid in self.parameters:
            return
        self.parameters[pid] = self._annotate(
            Parameter(
                pid,
                1.0,
                U.dimensionless,
                constant=True,
                name=f"activity {protein}",
                sboTerm=SBO.RELATIVE_ACTIVITY,
                notes=f"Base activity of {protein} in the organ relative to normal "
                "(1.0), set by genotype, impairment and the experiments; it "
                "multiplies the rates of its processes. Inhibition and induction "
                "act on top of it (the turnover state and the inhibitor terms of "
                "the rates).",
                port=True,
            ),
            "activity",
            protein=protein,
        )
        impaired = organ_activity(self.organ, protein)
        if impaired is not None:
            self._input(
                impaired,
                U.dimensionless,
                f"impairment of {protein}",
                f"The factor of the impaired organ function on the activity of "
                f"{protein} in the organ (`pkpdsim.pbpk.impairment`); the rates "
                f"run with `{pid}` times it.",
            )
            self.replacements[impaired] = ids.impairment_activity(protein, self.organ)
        for sid, inh in self.ctx.inhibitions:
            if inh.protein != protein:
                continue
            for par, field in inhibition_parameters(sid, inh):
                unit = U.per_min if field == "kinact" else U.mM
                self._input(par, unit, f"{field} {protein} {sid}")
        for sid, ind in self.ctx.inductions:
            if ind.protein == protein:
                self._input(
                    ids.emax(protein, sid), U.dimensionless, f"Emax {protein} {sid}"
                )
                self._input(ids.ec50(protein, sid), U.mM, f"EC50 {protein} {sid}")
        if has_turnover(protein, self.ctx):
            self._turnover(protein)

    def _turnover(self, protein: str) -> None:
        """The turnover state `E_<protein>` under mechanism-based inhibition and induction."""
        kdeg = ids.kdeg(protein)
        self.parameters[kdeg] = self._annotate(
            value_parameter(
                kdeg,
                turnover(self.organ),
                U.per_min,
                f"degradation rate {protein} [1/min]",
                SBO.DECAY_CONSTANT,
                f"First order degradation of {protein} in the organ; the synthesis "
                "of the baseline amount balances it.",
            ),
            "kdeg",
            protein=protein,
        )
        synthesis = [ONE]
        inactivation: list[str] = []
        for sid, ind in self.ctx.inductions:
            if ind.protein != protein:
                continue
            i = self._unbound_in_cell(sid, protein)
            synthesis.append(
                f"{ids.emax(protein, sid)} * ({i}) / ({ids.ec50(protein, sid)} + {i})"
            )
        for sid, inh in self.ctx.inhibitions:
            if inh.protein != protein or inh.kind != "mechanism":
                continue
            i = self._unbound_in_cell(sid, protein)
            inactivation.append(
                f"{ids.kinact(protein, sid)} * ({i}) / ({ids.ki(protein, sid)} + {i})"
            )
        e = ids.enzyme(protein)
        formula = f"{kdeg} * ({' + '.join(synthesis)}) - {kdeg} * {e}"
        if inactivation:
            formula += f" - ({' + '.join(inactivation)}) * {e}"
        self.parameters[e] = self._annotate(
            Parameter(
                e,
                1.0,
                U.dimensionless,
                constant=False,
                name=f"amount {protein} relative to baseline",
                sboTerm=SBO.RELATIVE_ACTIVITY,
                notes=f"Turnover of {protein} in the organ (initially 1): synthesis "
                "`kdeg*(1 + sum Emax*I/(EC50 + I))` with the induction of Almond "
                "2009, degradation `kdeg*E` and the mechanism-based inactivation "
                "`kinact*I/(KI + I)*E` of Mayhew 2000, combined as in Fahmi 2008; "
                "`I` is the unbound concentration in the cell. It multiplies the "
                f"rates of {protein}.\n\n`d{e}/dt = {formula}`\n\nSources:\n\n"
                f"- {cite('Mayhew2000')}\n- {cite('Almond2009')}\n- {cite('Fahmi2008')}",
            ),
            "enzyme_amount",
            protein=protein,
        )
        self.rate_rules.append(RateRule(e, formula, U.per_min))

    def _unbound_in_cell(self, sid: str, protein: str) -> str:
        """The unbound concentration of a substance in the cell of the organ.

        Intracellular, else the equivalent tissue part (merged tissue, whole
        organ), else plasma; the same for every protein (a variant with
        several cells weights them by the protein).
        """
        del protein
        s = next(s for s in self.ctx.substances if s.sid == sid)
        return inhibitor_unbound(Part.INTRACELLULAR, s, self.organ, self.ctx.parts)

    def model_definition(self, notes: str) -> ModelDefinition:
        """Get the model definition of the collected elements.

        Args:
            notes: the notes of the model definition.

        Returns:
            The model definition `<name>_<organ>`.
        """
        return ModelDefinition(
            sid=f"{self.name}_{self.organ.value}",
            name=f"{label(self.organ)} ({self.name} organ)",
            notes=notes,
            units=U,
            model_units=MODEL_UNITS,
            compartments=list(self.compartments.values()),
            species=list(self.species.values()),
            parameters=list(self.parameters.values()),
            rules=self.rules,
            assignments=self.assignments,
            rate_rules=self.rate_rules,
            reactions=self.reactions,
        )
