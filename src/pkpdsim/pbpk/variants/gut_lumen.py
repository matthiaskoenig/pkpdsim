"""The lumen of the intestine variants: transit, dissolution, absorption, bile, doses.

`LumenBuilder` adds the lumen segments to the standard organ of the
intestinal wall; `pkpdsim.pbpk.variants.gut.GutBuilder` derives from it and
sets the attributes (the segments, the substances which reach the lumen,
the mucosa) before it builds.
"""

from sbmlutils.factory import (
    Compartment,
    Parameter,
    PortType,
    RateRule,
    Species,
)
from sbmlutils.factory import Port as _Port
from sbmlutils.metadata import SBO

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.data import cite, value
from pkpdsim.pbpk.impairment import GASTRIC_EMPTYING, SI_TRANSIT
from pkpdsim.pbpk.organ import StandardBuilder
from pkpdsim.pbpk.organ_builder import ONE
from pkpdsim.pbpk.organs import Part
from pkpdsim.pbpk.parts import Parts
from pkpdsim.pbpk.routes import (
    Route,
    dissolution_inputs,
    dissolution_rate,
    is_weibull,
)
from pkpdsim.pbpk.substance import Substance
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Value
from pkpdsim.pbpk.variants.segments import DATA, SMALL, Segment

__all__ = ["LumenBuilder"]

_SBO_PARAMETER = SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER


class LumenBuilder(StandardBuilder):
    """The lumen segments of the intestine on top of the standard organ of the wall.

    Attributes:
        segments: the lumen segments in the order of the transit.
        single: the single lumen (absorption with `fabs`, no fed state).
        oral: the oral routes by (kind, substance).
        bile: the substances excreted into the bile.
        luminal: the substances which reach the lumen.
        mucosal: the substances with mucosa species.
    """

    segments: tuple[Segment, ...]
    single: bool
    oral: dict[tuple[str, str], Route]
    bile: set[str]
    luminal: set[str]
    mucosal: set[str]

    def _lumen_segments(self) -> None:
        """The compartments, transit times and area factors of the lumen segments."""
        for seg in self.segments:
            v = value(DATA, seg.key, "volume")
            cid = ids.lumen_volume(seg.token)
            self.compartments[cid] = self._annotate(
                Compartment(
                    cid,
                    v.value,
                    U.litre,
                    constant=True,
                    name=f"lumen {seg.label}",
                    sboTerm=SBO.PHYSICAL_COMPARTMENT,
                    notes=_notes("Fluid volume of the lumen segment.", v),
                    port=True,
                ),
                "lumen",
                segment=seg.token,
            )
            self._constant(
                ids.transit_time(seg.token),
                value(DATA, seg.key, "transit_time"),
                U.min,
                f"transit time {seg.label} [min]",
                "Mean transit time of the lumen segment (fasted), the transit "
                "into the next segment is `V*C/T`.",
                kind="gastric_transit_time" if seg.token == "st" else "transit_time",
                segment=seg.token,
            )
            if seg.absorbs:
                self._constant(
                    ids.area_factor(str(seg.token)),
                    value(DATA, seg.key, "area_factor"),
                    U.dimensionless,
                    f"absorption area factor {seg.label}",
                    "Mucosal surface per lumen volume relative to the small "
                    "intestine: the absorption of the segment is `ka*f_area`.",
                    kind="intestine_geometry",
                    segment=seg.token,
                )
        if not self.single:
            self._constant(
                ids.T_STOMACH_FED,
                value(DATA, "stomach", "transit_time_fed"),
                U.min,
                "transit time stomach fed [min]",
                "Mean transit time of the stomach in the fed state (`FED = 1`).",
                kind="gastric_transit_time",
                segment="st",
            )

    def _transit(self, seg: Segment) -> str:
        """The transit time of a segment, of the stomach by the prandial state.

        The impairment (chronic kidney disease) multiplies the transit time of
        the stomach and of the small intestine (its segments, the single
        lumen) by its factors, inputs of the body. The time is returned in
        parentheses, so that a caller divides by the whole product.
        """
        t = ids.transit_time(seg.token)
        if seg.token == "st":
            fed = self._input(ids.FED, U.dimensionless, "fed state (0 fasted, 1 fed)")
            factor = self._input(
                GASTRIC_EMPTYING,
                U.dimensionless,
                "impairment factor gastric emptying time",
            )
            return f"(({t} + {fed} * ({ids.T_STOMACH_FED} - {t})) * {factor})"
        if seg.token is None or seg.token in SMALL:
            factor = self._input(
                SI_TRANSIT,
                U.dimensionless,
                "impairment factor small intestinal transit time",
            )
            return f"({t} * {factor})"
        return t

    def _gut_species(
        self,
        sid: str,
        compartment: str,
        name: str,
        amount: bool = False,
        *,
        kind: str,
        substance: Substance,
        segment: str | None,
    ) -> str:
        """A species of a substance in the intestine, initially zero."""
        self.species[sid] = self._annotate(
            Species(
                sid,
                compartment,
                initialAmount=0.0 if amount else None,
                initialConcentration=None if amount else 0.0,
                substanceUnit=U.mmole,
                hasOnlySubstanceUnits=amount,
                name=name,
                sboTerm=SBO.SIMPLE_CHEMICAL,
                port=True,
            ),
            kind,
            segment=segment,
            substance=substance,
        )
        return sid

    def _dose_input(self, pid: str, name: str) -> str:
        """A dose state of the body, which the rate rule of the intestine changes."""
        self.parameters[pid] = Parameter(
            pid,
            0.0,
            U.mg,
            constant=False,
            name=name,
            sboTerm=_SBO_PARAMETER,
            notes="Dose state of the body, replaced by it; the intestine applies "
            "it and its rate rule changes it.",
            port=_Port(f"{pid}_port", portType=PortType.INPUT_PORT),
        )
        self.replacements[pid] = pid
        return pid

    def _absorbed_into(self, s: Substance, parts: Parts, seg: Segment) -> str:
        """The species the absorbed drug enters: the mucosa, else the plasma of the wall."""
        if s.sid in self.mucosal:
            return ids.mucosa(s.sid, seg.token)
        part = Part.ORGAN if parts.tissue == "organ" else Part.PLASMA
        return ids.species(part, s.sid)

    def _lumen(self, s: Substance, parts: Parts) -> None:
        """The lumen of a substance: transit, absorption, bile and oral doses."""
        sid = s.sid
        lumen = [
            self._gut_species(
                ids.lumen(sid, seg.token),
                ids.lumen_volume(seg.token),
                f"{s.name} (lumen {seg.label})",
                kind="concentration_lumen",
                substance=s,
                segment=seg.token,
            )
            for seg in self.segments
        ]
        self._input_compartment(ids.VFECES, "feces")
        feces = self._input_species(
            ids.feces(sid), ids.VFECES, f"{s.name} (feces)", amount=True
        )
        ka = self._input(ids.absorption(sid), U.per_hr, f"absorption {sid} [1/hr]")
        for i, seg in enumerate(self.segments):
            c, v = lumen[i], ids.lumen_volume(seg.token)
            target = lumen[i + 1] if i + 1 < len(lumen) else feces
            self._reaction(
                ids.transit_reaction(sid, seg.token),
                f"{c} -> {target}",
                f"{v} * {c} / {self._transit(seg)}",
                f"transit {sid} {seg.label}",
                kind="transit",
                segment=seg.token,
            )
            if seg.absorbs:
                into = self._absorbed_into(s, parts, seg)
                self._reaction(
                    ids.absorption_reaction(sid, seg.token),
                    f"{c} -> {into}",
                    f"{ka} / 60 min_per_hr * {ids.area_factor(str(seg.token))} * "
                    f"{v} * {c}",
                    f"absorption {sid} {seg.label}",
                    kind="absorption",
                    segment=seg.token,
                )
            elif self.single:
                into = self._absorbed_into(s, parts, seg)
                fabs = self._input(
                    ids.fraction_absorbed(sid),
                    U.dimensionless,
                    f"fraction absorbed {sid}",
                )
                rate = f"{ka} / 60 min_per_hr * {v} * {c}"
                self._reaction(
                    ids.absorption_reaction(sid),
                    f"{c} -> {into}",
                    f"{fabs} * {rate}",
                    f"absorption {sid}",
                    kind="absorption",
                )
                self._reaction(
                    ids.feces_reaction(sid),
                    f"{c} -> {feces}",
                    f"({ONE} - {fabs}) * {rate}",
                    f"unabsorbed fraction {sid} into the feces",
                    kind="fecal_excretion",
                )
        self._bile(s, lumen[0 if self.single else 1])
        tablet = self.oral.get(("oral", sid))
        if tablet is not None:
            self._tablet(s, tablet, lumen, feces)
        if ("oral_solution", sid) in self.oral:
            self._solution(s, lumen[0])

    def _bile(self, s: Substance, into: str) -> None:
        """The transit of the bile of a substance into the duodenum (or the feces).

        Every substance drains its bile: one with lumen species into the
        duodenum (the single lumen), any other into the feces, so nothing
        accumulates in the bile.
        """
        sid = s.sid
        self._input_compartment(ids.VBILE, "bile")
        bile = self._input_species(ids.bile(sid), ids.VBILE, f"{s.name} (bile)")
        k_bile = self._input(
            ids.K_BILE, U.per_min, "transit bile to intestinal lumen [1/min]"
        )
        into_feces = into == ids.feces(sid)
        self._reaction(
            ids.bile_reaction(sid),
            f"{bile} -> {into}",
            f"{k_bile} * {ids.VBILE} * {bile}",
            f"transit {sid} bile to "
            + ("feces (no lumen species)" if into_feces else "intestinal lumen"),
            kind="fecal_excretion" if into_feces else "bile_flow",
        )

    def _bile_to_feces(self, s: Substance) -> None:
        """The bile of a substance without lumen species drains into the feces."""
        self._input_compartment(ids.VFECES, "feces")
        feces = self._input_species(
            ids.feces(s.sid), ids.VFECES, f"{s.name} (feces)", amount=True
        )
        self._bile(s, feces)

    def _solution(self, s: Substance, stomach: str) -> None:
        """The oral solution: swallowed into the lumen of the stomach, dissolved."""
        sid = s.sid
        dose = self._dose_input(ids.sol_dose(sid), f"oral dose {sid} (solution) [mg]")
        mr = self._input(ids.mr(sid), U.g_per_mole, f"molecular weight {sid} [g/mole]")
        k = self._input(
            ids.solution(sid), U.per_hr, f"application rate {sid} (solution) [1/hr]"
        )
        reaction = ids.dose_reaction("sol", sid)
        self._reaction(
            reaction,
            f"-> {stomach}",
            f"{k} / 60 min_per_hr * {dose} / {mr}",
            f"oral application {sid} (solution)",
            kind="oral_application",
            segment=self.segments[0].token,
        )
        self.rate_rules.append(
            RateRule(dose, f"-{reaction} * {mr}", U.mg_per_min, name=f"oral dose {sid}")
        )

    def _tablet(self, s: Substance, route: Route, lumen: list[str], feces: str) -> None:
        """The oral tablet: dissolution and transit of the undissolved drug.

        The dose state is the undissolved drug of the first segment; the
        other segments have the undissolved amount `Asolid_<seg>_<s>`.
        """
        sid = s.sid
        dose = self._dose_input(ids.po_dose(sid), f"oral dose {sid} (tablet) [mg]")
        mr = self._input(ids.mr(sid), U.g_per_mole, f"molecular weight {sid} [g/mole]")
        for pid, unit, name in dissolution_inputs(route):
            self._input(pid, unit, name)
        solids = [
            self._gut_species(
                ids.solid(sid, str(seg.token)),
                ids.lumen_volume(seg.token),
                f"{s.name} (undissolved, lumen {seg.label})",
                amount=True,
                kind="amount_undissolved",
                substance=s,
                segment=seg.token,
            )
            for seg in self.segments[1:]
        ]
        undissolved = f"({dose} + {mr} * ({' + '.join(solids)}))" if solids else dose
        age = ids.dissolution_age(sid)
        kdis, reads = dissolution_rate(route, undissolved, age)
        if reads:
            kdis = self._rule(
                ids.dissolution_rate(sid),
                kdis,
                U.per_min,
                f"dissolution rate {sid} [1/min]",
                "Dissolution rate of the undissolved drug in every segment.",
                kind="dissolution",
            )
        else:
            kdis = f"({kdis})"
        first = ids.dose_reaction("po", sid)
        self._reaction(
            first,
            f"-> {lumen[0]}",
            f"{kdis} * {dose} / {mr}",
            f"dissolution {sid} (tablet) {self.segments[0].label}",
            kind="dissolution_reaction",
            segment=self.segments[0].token,
        )
        empty = ids.transit_reaction(sid, self.segments[0].token, solid=True)
        self._reaction(
            empty,
            f"-> {solids[0] if solids else feces}",
            f"{dose} / {mr} / {self._transit(self.segments[0])}",
            f"transit {sid} undissolved {self.segments[0].label}",
            kind="transit",
            segment=self.segments[0].token,
        )
        self.rate_rules.append(
            RateRule(
                dose,
                f"-({first} + {empty}) * {mr}",
                U.mg_per_min,
                name=f"oral dose {sid} (tablet)",
            )
        )
        dissolved, last = [first], empty
        for i, seg in enumerate(self.segments[1:], start=1):
            a = solids[i - 1]
            dis = ids.dissolution_reaction(sid, str(seg.token))
            self._reaction(
                dis,
                f"{a} -> {lumen[i]}",
                f"{kdis} * {a}",
                f"dissolution {sid} {seg.label}",
                kind="dissolution_reaction",
                segment=seg.token,
            )
            last = ids.transit_reaction(sid, seg.token, solid=True)
            self._reaction(
                last,
                f"{a} -> {solids[i] if i < len(solids) else feces}",
                f"{a} / {self._transit(seg)}",
                f"transit {sid} undissolved {seg.label}",
                kind="transit",
                segment=seg.token,
            )
            dissolved.append(dis)
        if is_weibull(route):
            self._age(
                sid,
                f"{dose} / {mr}" + "".join(f" + {a}" for a in solids),
                [
                    *dissolved,
                    last,
                ],
            )

    def _age(self, sid: str, pool: str, outflows: list[str]) -> None:
        """The mean time since the dose of the undissolved drug.

        `dW/dt = A - (sum J_out)*age`, `age = W/A` with the undissolved
        amount `A`: every undissolved particle ages, the outflows (dissolution
        and the transit into the feces) take the mean age with them.
        """
        w, age = ids.dissolution_clock(sid), ids.dissolution_age(sid)
        self.parameters[w] = self._annotate(
            Parameter(
                w,
                0.0,
                U.mmole_min,
                constant=False,
                name=f"age-weighted undissolved amount {sid} [mmole*min]",
                sboTerm=_SBO_PARAMETER,
                notes="The undissolved amount weighted by its time since the dose; "
                f"`{age} = {w}/A` is the mean age of the undissolved amount `A`.",
            ),
            "dissolution",
        )
        self.rate_rules.append(
            RateRule(
                w,
                f"{pool} - ({' + '.join(outflows)}) * {age}",
                U.mmole,
                name=f"age-weighted undissolved amount {sid}",
            )
        )
        self._rule(
            age,
            f"piecewise({w} / ({pool}), ({pool}) > 1e-18 mmole, 0 min)",
            U.min,
            f"mean time since the dose of the undissolved drug {sid} [min]",
            "The time since the dose of a Weibull dissolution: exact for one "
            "dose, the amount-weighted mean of doses which overlap.",
            kind="dissolution",
        )


def _notes(text: str, v: Value) -> str:
    """The notes of a value of the data with its source."""
    notes = f"{text} {v.notes}" if v.notes else text
    return f"{notes}\n\nSource:\n\n- {cite(str(v.source))}"
