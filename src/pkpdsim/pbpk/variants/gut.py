"""The builder of the intestine variants: wall, mucosa and chemistry.

`GutBuilder` builds the submodel `GU` of the segmented intestine and of the
single lumen (`pkpdsim.pbpk.variants.intestine`): the wall is the standard
organ of small and large intestine, the lumen that of `LumenBuilder`, the
mucosa and the chemistry of the small intestine are placed here.
"""

from collections.abc import Callable
from dataclasses import replace

from sbmlutils.factory import ModelDefinition, Reaction
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
from pkpdsim.pbpk.organ import Boundary, OrganContext
from pkpdsim.pbpk.organ_builder import ONE
from pkpdsim.pbpk.organs import Organ, Part, Pool
from pkpdsim.pbpk.parts import Parts, parts_for
from pkpdsim.pbpk.placement import (
    inhibitor_terms,
    rate_law,
    stoichiometry,
    with_modifiers,
)
from pkpdsim.pbpk.routes import ORAL
from pkpdsim.pbpk.substance import Kp, Substance
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Value
from pkpdsim.pbpk.variants.gut_lumen import LumenBuilder
from pkpdsim.pbpk.variants.segments import DATA, Segment

__all__ = ["GUT_ORGANS", "GutBuilder"]

GUT_ORGANS = (Organ.SMALL_INTESTINE, Organ.LARGE_INTESTINE)
"""The organs the intestine variants replace."""

_PROFILED = {"CYP3A4": "cyp3a4", "ABCB1": "abcb1", "ABCG2": "abcg2"}
"""The proteins with an expression of the data, by their key."""

_UNIFORM = "uniform"
"""The key of the placement by mucosa volume."""


def _gut_chemistry(
    ctx: OrganContext,
) -> tuple[
    list[Metabolism | Transport], tuple[Metabolism | Transport | CustomReaction, ...]
]:
    """The processes of the intestine and the others.

    Raises:
        ValueError: for a process of the large intestine, a custom reaction or
            a basolateral transport of the intestine.
    """
    gut: list[Metabolism | Transport] = []
    rest = []
    for p in ctx.chemistry:
        if p.organ not in GUT_ORGANS:
            rest.append(p)
            continue
        if p.organ is Organ.LARGE_INTESTINE:
            raise ValueError(
                f"Process '{p.sid}' is of the large intestine: the intestine "
                "places the chemistry of the small intestine (organ "
                "SMALL_INTESTINE) into the mucosa of all segments."
            )
        if isinstance(p, CustomReaction):
            raise ValueError(
                f"Custom reaction '{p.sid}' of the small intestine: the "
                "intestine variants place metabolism and transport only."
            )
        if isinstance(p, Transport) and p.membrane is Membrane.BASOLATERAL:
            raise ValueError(
                f"Transport '{p.sid}' is basolateral: the mucosa of the intestine "
                "is flow-limited, it has apical transport (lumen) only."
            )
        gut.append(p)
    return gut, tuple(rest)


def _check_wall(substances: tuple[Substance, ...]) -> None:
    """Check that no substance gives the large intestine a distribution of its own.

    Raises:
        ValueError: for a Kp or an exchange of the large intestine, which the
            intestine (with the distribution of the small intestine) ignores.
    """
    for s in substances:
        dist = s.distribution
        if dist is None:
            continue
        own = isinstance(dist, Kp) and Organ.LARGE_INTESTINE in dist.overrides
        if own or Organ.LARGE_INTESTINE in dist.ps:
            raise ValueError(
                f"{s.sid}: a Kp or exchange of the large intestine is not used, the "
                "intestine wall has the distribution of the small intestine."
            )


class GutBuilder(LumenBuilder):
    """The SBML objects of the intestine: wall, lumen segments, mucosa and chemistry."""

    def __init__(
        self,
        ctx: OrganContext,
        name: str,
        segments: tuple[Segment, ...],
        *,
        single: bool,
    ) -> None:
        """Start the collection of the elements of the intestine.

        The wall is the standard organ without the chemistry of the
        intestine, which goes into the mucosa.

        Args:
            ctx: the context of the small intestine.
            name: the name of the component.
            segments: the lumen segments.
            single: the single lumen.

        Raises:
            ValueError: for chemistry the intestine cannot place, or a
                distribution of the large intestine.
        """
        processes, rest = _gut_chemistry(ctx)
        _check_wall(ctx.substances)
        wall = replace(
            ctx,
            parts={
                s.sid: parts_for(s, Organ.SMALL_INTESTINE, rest) for s in ctx.substances
            },
            chemistry=rest,
        )
        super().__init__(wall, name)
        self.processes = processes
        self.segments = segments
        self.single = single
        self.mucosa_segments = tuple(seg for seg in segments if seg.region)
        proteins = {p.protein for p in processes if p.protein is not None}
        mucosal = {
            sid
            for p in processes
            for sid in (
                [p.substrate, *p.products]
                if isinstance(p, Metabolism)
                else [p.substance]
            )
        }
        mucosal |= {sid for sid, i in ctx.inhibitions if i.protein in proteins}
        mucosal |= {sid for sid, i in ctx.inductions if i.protein in proteins}
        self.mucosal = mucosal
        self.bile = {
            p.substance
            for p in ctx.chemistry
            if isinstance(p, Transport)
            and p.organ is Organ.LIVER
            and p.membrane is Membrane.APICAL
            and p.direction is Direction.EFFLUX
        }
        influx = {
            p.protein
            for p in processes
            if isinstance(p, Transport) and p.direction is Direction.INFLUX
        }
        self.oral = {(r.kind, r.substance): r for r in ctx.routes if r.kind in ORAL}
        self.luminal = (
            {s for _, s in self.oral}
            | self.bile
            | {p.substance for p in processes if isinstance(p, Transport)}
            | {sid for sid, i in ctx.inhibitions if i.protein in influx}
        )

    def build(self) -> tuple[ModelDefinition, Boundary]:
        """Build the model definition and the boundary."""
        self._compartments()
        self._lumen_segments()
        if self.mucosal:
            self._mucosa()
        for s in self.ctx.substances:
            parts = self.ctx.parts[s.sid]
            self._substance(s)
            if s.sid in self.mucosal:
                self._mucosa_substance(s, parts)
            if s.sid in self.luminal:
                self._lumen(s, parts)
        self._chemistry()
        kind = "single lumen" if self.single else "segmented"
        md = self.model_definition(
            f"The {kind} intestine, which replaces small and large intestine: "
            "the intestinal wall, the lumen with the transit, dissolution and "
            "absorption, the mucosa and the chemistry of the intestine."
        )
        md.name = f"intestine ({self.name})"
        return md, Boundary(
            replacements=dict(self.replacements), outputs=dict(self.outputs)
        )

    # --- the wall: the standard organ of small and large intestine ---------

    @property
    def v(self) -> str:
        """The volume of the intestine, `Vsi + Vla`."""
        if ids.GUT_VOLUME not in self.parameters:
            vsi, vla = (
                self._input(
                    ids.organ_volume(o), U.litre, f"volume {o.name.lower()} [l]"
                )
                for o in GUT_ORGANS
            )
            self._initial(
                ids.GUT_VOLUME,
                f"{vsi} + {vla}",
                U.litre,
                "volume intestine [l]",
                "Volume of small and large intestine.",
            )
        return ids.GUT_VOLUME

    def _fraction(self, pid: str, fraction: Callable[[Organ], str], name: str) -> str:
        """A fraction of the intestine, the volume-weighted fraction of both organs."""
        if pid not in self.parameters:
            v = self.v
            terms = [
                f"{ids.organ_volume(o)} * "
                + self._input(fraction(o), U.dimensionless, f"{name} {o.name.lower()}")
                for o in GUT_ORGANS
            ]
            self._initial(
                pid,
                f"({' + '.join(terms)}) / {v}",
                U.dimensionless,
                f"{name} intestine",
                f"The {name} of small and large intestine, weighted by volume.",
            )
        return pid

    @property
    def f_vas(self) -> str:
        """The vascular fraction of the intestine."""
        return self._fraction(
            ids.GUT_VASCULAR_FRACTION, ids.vascular_fraction, "vascular fraction"
        )

    @property
    def f_int(self) -> str:
        """The interstitial fraction of the intestine."""
        return self._fraction(
            ids.GUT_INTERSTITIAL_FRACTION,
            ids.interstitial_fraction,
            "interstitial fraction",
        )

    def _flow(self, fid: str) -> str:
        """A blood flow: the local flows of the intestine or an input."""
        if fid in (ids.GUT_FLOW, ids.GUT_WALL_FLOW):
            return fid
        return super()._flow(fid)

    def _gut_flow(self) -> str:
        """The blood flow of the intestine, `Qsi + Qla`."""
        if ids.GUT_FLOW not in self.parameters:
            flows = [
                super(GutBuilder, self)._flow(ids.organ_flow(o)) for o in GUT_ORGANS
            ]
            self._initial(
                ids.GUT_FLOW,
                " + ".join(flows),
                U.l_per_min,
                "blood flow intestine [l/min]",
                "Blood flow of small and large intestine.",
            )
        return ids.GUT_FLOW

    def _flows(self, s: Substance) -> tuple[list[tuple[Pool, str]], str]:
        """The arterial inflow of the wall and its outflow into the portal vein.

        A substance with a mucosa perfuses the wall with the flow past the
        mucosa, any other the wall with the whole intestinal flow.
        """
        flow = self._gut_flow()
        if s.sid in self.mucosal:
            flow = ids.GUT_WALL_FLOW
        return [(Pool.ARTERIAL, flow)], flow

    # --- the mucosa -----------------------------------------------------

    def _mucosa(self) -> None:
        """The volumes and blood flows of the mucosa segments and the flow of the wall."""
        f_vol = self._constant(
            ids.F_MUCOSA_VOLUME,
            value(DATA, "mucosa", "volume_fraction"),
            U.dimensionless,
            "mucosa fraction of the intestinal tissue",
        )
        f_flow = self._constant(
            ids.F_MUCOSA_FLOW,
            value(DATA, "mucosa", "flow_fraction"),
            U.dimensionless,
            "mucosa fraction of the intestinal blood flow",
        )
        volumes, flows = [], []
        for seg in self.mucosa_segments:
            region = seg.region
            assert region is not None
            share = f"{seg.share:.6g} dimensionless"
            vo = self._input(
                ids.organ_volume(region), U.litre, f"volume {region.name.lower()} [l]"
            )
            fv = self._input(
                ids.vascular_fraction(region),
                U.dimensionless,
                f"vascular fraction {region.name.lower()}",
            )
            cid = ids.mucosa_volume(seg.token)
            self._compartment(
                cid,
                f"{f_vol} * {vo} * ({ONE} - {fv}) * {share}",
                f"mucosa {seg.key.replace('_', ' ')}",
                "Mucosa (enterocytes) of the segment: its share of the tissue of "
                "the region by the transit time.",
            )
            q = self._initial(
                ids.mucosa_flow(seg.token),
                f"{f_flow} * {super()._flow(ids.organ_flow(region))} * {share}",
                U.l_per_min,
                f"blood flow mucosa {seg.key.replace('_', ' ')} [l/min]",
                "Blood flow of the mucosa of the segment: its share of the blood "
                "flow of the region by the transit time.",
            )
            volumes.append(cid)
            flows.append(q)
        self._initial(
            ids.MUCOSA_VOLUME_TOTAL,
            " + ".join(volumes),
            U.litre,
            "volume mucosa [l]",
            "Volume of the mucosa of all segments.",
        )
        self._initial(
            ids.MUCOSA_FLOW_TOTAL,
            " + ".join(flows),
            U.l_per_min,
            "blood flow mucosa [l/min]",
            "Blood flow of the mucosa of all segments.",
        )
        self._initial(
            ids.GUT_WALL_FLOW,
            f"{self._gut_flow()} - {ids.MUCOSA_FLOW_TOTAL}",
            U.l_per_min,
            "blood flow intestinal wall [l/min]",
            "Blood flow of the wall past the mucosa, for a substance with a mucosa.",
        )

    def _mucosa_substance(self, s: Substance, parts: Parts) -> None:
        """The mucosa species of a substance and their arterial perfusion.

        Flow-limited: `Qmu*carrier*Car` in, `Qmu*carrier*Cmu/K_mu` out into
        the portal vein; the kinetic blood cells of the mucosal blood pass it.
        """
        sid = s.sid
        kp = self._kp(s)
        k = self._initial(
            ids.mucosa_coefficient(sid),
            kp if kp is not None else ONE,
            U.dimensionless,
            f"mucosa to plasma coefficient {sid}",
            "Mucosa to plasma coefficient of the flow-limited mucosa: the Kp of "
            "the small intestine, 1 for a substance without tissue distribution.",
        )
        if parts.flow_carrier == "blood_bp":
            carrier = self._bp(s)
        else:
            carrier = f"({ONE} - {self.hct})"
        car = self._pool_species(Pool.ARTERIAL, Part.PLASMA, s, parts)
        cpo = self._pool_species(Pool.PORTAL, Part.PLASMA, s, parts)
        for seg in self.mucosa_segments:
            q = ids.mucosa_flow(seg.token)
            c = self._gut_species(
                ids.mucosa(sid, seg.token),
                ids.mucosa_volume(seg.token),
                f"{s.name} (mucosa {seg.key.replace('_', ' ')})",
            )
            self._reaction(
                ids.mucosa_reaction("in", sid, seg.token),
                f"{car} -> {c}",
                f"{q} * {carrier} * {car}",
                f"inflow {sid} into mucosa {seg.key.replace('_', ' ')}",
            )
            self._reaction(
                ids.mucosa_reaction("out", sid, seg.token),
                f"{c} -> {cpo}",
                f"{q} * {carrier} * {c} / {k}",
                f"outflow {sid} from mucosa {seg.key.replace('_', ' ')}",
            )
            if parts.blood_cells:
                car_bc = self._pool_species(Pool.ARTERIAL, Part.BLOOD_CELLS, s, parts)
                cpo_bc = self._pool_species(Pool.PORTAL, Part.BLOOD_CELLS, s, parts)
                self._reaction(
                    ids.mucosa_reaction("bc", sid, seg.token),
                    f"{car_bc} -> {cpo_bc}",
                    f"{q} * {self.hct} * {car_bc}",
                    f"blood cells {sid} past mucosa {seg.key.replace('_', ' ')}",
                )

    def _mucosa_unbound(self, s: Substance, seg: Segment) -> str:
        """The unbound concentration of a substance in the mucosa of a segment."""
        fu = self._fu(s)
        return (
            f"{fu} / {ids.mucosa_coefficient(s.sid)} * {ids.mucosa(s.sid, seg.token)}"
        )

    def _unbound_in_cell(self, sid: str) -> str:
        """The unbound concentration in the cell, the mean over the mucosa segments."""
        if sid not in self.mucosal:
            return super()._unbound_in_cell(sid)
        s = next(s for s in self.ctx.substances if s.sid == sid)
        terms = " + ".join(
            f"{ids.mucosa_volume(seg.token)} * {self._mucosa_unbound(s, seg)}"
            for seg in self.mucosa_segments
        )
        return f"({terms}) / {ids.MUCOSA_VOLUME_TOTAL}"

    # --- the chemistry --------------------------------------------------

    def _weight(self, protein: str | None, seg: Segment) -> str | None:
        """The share of a process of a protein in the mucosa of a segment.

        `w = e*Vmu/sum(e*Vmu)` with the expression `e` of the protein of the
        data, by mucosa volume (`e = 1`) for a protein without; `None` where
        the protein is not expressed.
        """
        key = _UNIFORM if self.single else _PROFILED.get(protein or "", _UNIFORM)

        def expression(other: Segment) -> Value | None:
            if key == _UNIFORM:
                return None
            return value(DATA, other.key, key)

        e = expression(seg)
        if e is not None and e.value == 0.0:
            return None
        pid = ids.placement_weight(key, seg.token)
        if pid in self.parameters:
            return pid

        def term(other: Segment) -> str:
            ex = expression(other)
            v = ids.mucosa_volume(other.token)
            return v if ex is None else f"{ex.value!r} dimensionless * {v}"

        total = " + ".join(
            term(o)
            for o in self.mucosa_segments
            if (ex := expression(o)) is None or ex.value != 0.0
        )
        if key == _UNIFORM:
            text = "by mucosa volume (the protein has no expression data)"
        else:
            sources = sorted(
                {str(x.source) for o in self.mucosa_segments if (x := expression(o))}
            )
            text = (
                "by the relative expression times the mucosa volume\n\nSources:\n\n"
                + "\n".join(f"- {cite(src)}" for src in sources)
            )
        self._initial(
            pid,
            f"{term(seg)} / ({total})",
            U.dimensionless,
            f"share {key} mucosa {seg.key.replace('_', ' ')}",
            f"Share of a process of the protein in the mucosa of the segment, {text}",
        )
        return pid

    def _chemistry(self) -> None:
        """The processes of the small intestine in the mucosa segments."""
        for p in self.processes:
            if p.protein is not None:
                self._activity(p.protein)
            for seg in self.mucosa_segments:
                weight = self._weight(p.protein, seg)
                if weight is not None:
                    self._place(p, seg, weight)

    def _site(
        self, seg: Segment, lumen: bool
    ) -> Callable[[Substance], tuple[str, str]]:
        """The species and unbound concentration of an inhibitor in a segment."""

        def site(inhibitor: Substance) -> tuple[str, str]:
            if lumen:
                sp = ids.lumen(inhibitor.sid, seg.token)
                return sp, sp
            return ids.mucosa(inhibitor.sid, seg.token), self._mucosa_unbound(
                inhibitor, seg
            )

        return site

    def _place(self, p: Metabolism | Transport, seg: Segment, weight: str) -> None:
        """A process of the chemistry in the mucosa of a segment."""
        by_id = {s.sid: s for s in self.ctx.substances}
        name = seg.key.replace("_", " ")
        if isinstance(p, Metabolism):
            s = by_id[p.substrate]
            cu = self._mucosa_unbound(s, seg)
            products = " + ".join(
                stoichiometry(f, ids.mucosa(product, seg.token))
                for product, f in p.products.items()
            )
            equation = f"{ids.mucosa(s.sid, seg.token)} -> {products}"
            lumen, sbo = False, SBO.BIOCHEMICAL_REACTION
            what = f"{p.substrate} -> {', '.join(p.products) or 'elimination'}"
        else:
            s = by_id[p.substance]
            cell, outside = ids.mucosa(s.sid, seg.token), ids.lumen(s.sid, seg.token)
            lumen = p.direction is Direction.INFLUX
            if lumen:
                equation, cu = f"{outside} -> {cell}", outside
            else:
                equation, cu = f"{cell} -> {outside}", self._mucosa_unbound(s, seg)
            sbo = SBO.TRANSPORT_REACTION
            what = f"{p.direction.value} {s.sid} apical"
        inhibitors = inhibitor_terms(p, self.ctx, self._site(seg, lumen))
        rate, pars = rate_law(p, self.ctx, cu, weight, inhibitors)
        for par in pars:
            self.parameters.setdefault(str(par.sid), par)
        self.reactions.append(
            Reaction(
                sid=ids.segment_process(p.sid, seg.token),
                name=f"{what} mucosa {name}" + (f" ({p.protein})" if p.protein else ""),
                equation=with_modifiers(equation, inhibitors.modifiers),
                formula=(rate, U.mmole_per_min),
                compartment=ids.mucosa_volume(seg.token),
                sboTerm=sbo,
                notes=f"Process `{p.sid}` of the small intestine in the mucosa of "
                f"the {name}, its share `{weight}` of the capacity, on `{cu}`.",
            )
        )
