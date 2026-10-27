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
from sbmlutils.factory.units import UnitType
from sbmlutils.metadata import SBO

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.annotation import annotate, initial_assignment
from pkpdsim.pbpk.chemistry import CustomReaction
from pkpdsim.pbpk.data import cite
from pkpdsim.pbpk.impairment import organ_activity
from pkpdsim.pbpk.organs import Organ, Part, Pool
from pkpdsim.pbpk.parts import Parts
from pkpdsim.pbpk.placement import (
    ONE,
    FormulaParameter,
    formula_assignment,
    formula_ids,
    formula_problems,
    has_turnover,
    in_zone,
    inhibition_parameters,
    inhibitor_unbound,
    kp_id,
    species_compartment,
    species_start,
    start_formula,
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
        self.starts: list[tuple[str, Substance, Part | None, str | None]] = []

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

    def _input(self, pid: str, unit: UnitType, name: str, text: str = "") -> str:
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
        self.assignments.append(initial_assignment(pid, formula, unit, name))
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
        self.assignments.append(
            initial_assignment(
                cid, formula, U.litre, in_zone(f"{name} {label(self.organ)}", zone)
            )
        )

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
        self._defer_start(sid, s, part)

    def _defer_start(
        self,
        sid: str,
        s: Substance,
        part: Part | None = None,
        coefficient: str | None = None,
    ) -> None:
        """Record the start of a species of an endogenous substance, built at the end.

        The starts are built when the model definition is created
        (`model_definition`), after every exchange: the elements they add
        (the input `baseline_<s>`, a coefficient, the initial assignments)
        come after those of the organ, which keep their order, so the model
        of a substance without baseline is unchanged and that of one with a
        baseline only gains elements.

        Args:
            sid: the id of the species.
            s: the substance; nothing is recorded without baseline.
            part: the part of the species, whose coefficient the start reads
                (`_equilibrium`), `None` for `coefficient`.
            coefficient: the coefficient of a species of no part (the
                mucosa).
        """
        if s.baseline is not None:
            self.starts.append((sid, s, part, coefficient))

    def _equilibrium(self, part: Part, s: Substance) -> str | None:
        """The coefficient of a part to plasma at the equilibrium of its exchange.

        `Kbc` of the kinetic blood cells (`kex_bc*Vbc*(Cpl - Cbc/Kbc)`), the
        Kp of the merged tissue (`PS*fu*(Cpl - Cti/Kp)`), `K_int_pls` of the
        interstitium (`PS_vas*fu*(Cpl - Cis/K_int_pls)`), `K_cell_pls` of the
        cells (`PS*(fu/K_int_pls*Cis - fu/K_cell_pls*Cic)`), `K_or` of the
        flow-limited organ (it releases `Cor/K_or`); `None` for plasma.
        """
        if part is Part.BLOOD_CELLS:
            return self._input(
                ids.kbc(s.sid),
                U.dimensionless,
                f"blood cell to plasma coefficient {s.sid}",
            )
        if part is Part.TISSUE:
            return self._kp(s)
        if part is Part.INTERSTITIAL:
            return ids.k_int_pls(s.sid)
        if part is Part.INTRACELLULAR:
            return ids.k_cell_pls(s.sid)
        if part is Part.ORGAN:
            return ids.k_organ(s.sid)
        return None

    def _start(self, sid: str, s: Substance, coefficient: str | None) -> None:
        """The start of a species of an endogenous substance at its equilibrium with plasma.

        The initial assignment `coefficient * baseline_<s>` (`baseline_<s>`
        for plasma), the baseline an input of the organ; nothing for a
        substance without baseline, whose species start at 0.
        """
        start = start_formula(s, coefficient)
        if start is None:
            return
        self._input(ids.baseline(s.sid), U.mM, f"baseline {s.sid} [mM]")
        species = self.species[sid]
        species.initialConcentration = None
        self.assignments.append(species_start(sid, start, str(species.name)))

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
        self.rate_rules.append(
            RateRule(e, formula, U.per_min, name=f"turnover of {protein}")
        )

    def _unbound_in_cell(self, sid: str, protein: str) -> str:
        """The unbound concentration of a substance in the cell of the organ.

        Intracellular, else the equivalent tissue part (merged tissue, whole
        organ), else plasma; the same for every protein (a variant with
        several cells weights them by the protein).
        """
        del protein
        s = next(s for s in self.ctx.substances if s.sid == sid)
        return inhibitor_unbound(Part.INTRACELLULAR, s, self.organ, self.ctx.parts)

    def _formula(self, p: FormulaParameter) -> None:
        """The initial assignment of a `Formula` of the organ and its inputs.

        Every id of the formula is a constant parameter or compartment of the
        top model (`OrganContext.top`), read through an input port of the
        organ whose local id is the id of the top model: an input the organ
        already has under that id is used, else the input is added; in the
        flat model the formula reads the elements of the top model. An input
        whose local id differs from the id of the top model (the liver reads
        `fimp_li_cyp3a4` as `fimp_cyp3a4`) is no id of a formula, and a
        local element of the organ under the id of the top model is refused.
        """
        pid = str(p.sid)
        where = (
            f"The formula of '{pid}' in organ '{self.organ.value}' "
            f"(`{p.formula.expression}`)"
        )
        problems = formula_problems(where, p.formula.expression, self.ctx.top)
        clash: list[str] = []
        for symbol in sorted(formula_ids(p.formula.expression)):
            kind = self.ctx.top.get(symbol, ("", None))[0]
            if kind not in ("parameter", "compartment"):
                continue
            if self.replacements.get(symbol) == symbol:
                continue
            if (
                symbol in self.parameters
                or symbol in self.compartments
                or symbol in self.species
            ):
                clash.append(symbol)
            elif kind == "compartment":
                self._input_compartment(symbol, f"{symbol} (formula of {pid})")
            else:
                self._input(
                    symbol, self.ctx.top[symbol][1], f"{symbol} (formula of {pid})"
                )
        if clash:
            problems.append(
                f"{where} reads ids of the top model which are local elements of "
                f"the organ: {clash}"
            )
        if problems:
            raise ValueError("; ".join(problems) + ".")
        self.assignments.append(formula_assignment(p))

    def _custom_inputs(self) -> None:
        """The inputs of the ids of the top model which the custom reactions read.

        An id of the formula of a custom reaction is an element of the organ
        (its species, compartments, parameters, the parameters of the
        reaction), `time`, or a constant parameter or compartment of the top
        model (`OrganContext.top`, a parameter of the definition, a
        covariate), which the organ reads through an input port under the
        same id, as for the ids of a `Formula`.

        Raises:
            ValueError: if an id is neither, naming the reaction and the ids.
        """
        placed = {str(r.sid) for r in self.reactions}
        local = {
            *self.parameters,
            *self.compartments,
            *self.species,
            *placed,
            *(str(p.sid) for r in self.reactions for p in r.pars or []),
            "time",
        }
        for p in self.ctx.chemistry:
            if not isinstance(p, CustomReaction) or p.sid not in placed:
                continue
            refused: list[str] = []
            for symbol in sorted(formula_ids(p.formula) - local):
                kind, unit = self.ctx.top.get(symbol, ("", None))
                if kind == "compartment":
                    self._input_compartment(symbol, f"{symbol} (read by {p.sid})")
                elif kind == "parameter":
                    self._input(symbol, unit, f"{symbol} (read by {p.sid})")
                else:
                    refused.append(symbol)
            if refused:
                raise ValueError(
                    f"Custom reaction '{p.sid}' of organ '{self.organ.value}' reads "
                    f"{refused}, which are no elements of the organ nor constant "
                    "parameters or compartments of the top model (a parameter of "
                    "the definition, a covariate); a species or an observer of the "
                    "top model is not read by an organ"
                )
            local |= formula_ids(p.formula)

    def model_definition(self, notes: str) -> ModelDefinition:
        """Get the model definition of the collected elements.

        Args:
            notes: the notes of the model definition.

        Returns:
            The model definition `<name>_<organ>`.

        Raises:
            ValueError: if a `Formula` of the organ reads an id which is no
                parameter or compartment of the top model, or an element of
                the organ which is no input, or a custom reaction reads an
                id which is neither an element of the organ nor a parameter
                or compartment of the top model.
        """
        for sid, s, part, coefficient in self.starts:
            if part is not None:
                coefficient = self._equilibrium(part, s)
            self._start(sid, s, coefficient)
        self.starts = []
        self._custom_inputs()
        parameters = [
            *self.parameters.values(),
            *(p for r in self.reactions for p in (r.pars or [])),
        ]
        for p in parameters:
            if isinstance(p, FormulaParameter):
                self._formula(p)
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
