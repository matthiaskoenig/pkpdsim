"""The placement of the chemistry in the parts of the standard organ.

A process of the chemistry acts in one part of the organ: a metabolism in the
part which holds the tissue of its substrate (intracellular in the full form,
the merged tissue, the whole organ when flow-limited, the plasma without
tissue), a transport across the membrane of the cell (basolateral between
interstitium and cell, apical between cell and the sink of the organ: bile,
urine or the intestinal lumen; the standard kidney has no uptake from the
cumulative urine), a custom reaction in the part it names.

Every rate acts on the unbound concentration of its part and is an amount per
time: the Vmax of a process is per litre of organ, so the whole organ has the
capacity `Vmax * V<o>`; the clearance `k` of a first order process is per
litre of organ as well (`FirstOrder.per_volume`, the default, 1/min), so the
whole organ has `k * V<o>` as with a Vmax, or that of the whole organ (l/min). The activity
`f_<protein>` of its protein multiplies the rate (times the turnover state
`E_<protein>` when a substance of the definition inactivates or induces the
protein). The organ may give a factor of its capacity, which multiplies the
rate of every metabolism and transport as well: the functional fraction
`f_functional` of the liver, the share of the capacity in the functional
tissue, and the tubular capacity of the kidney, which follows the renal
function `f_renal_function` (the intact nephron hypothesis); a custom
reaction states its rate itself and is not scaled.

The reversible inhibitors of the protein enter the rate of every process of
the protein with their unbound concentration `I` in the compartment from
which the process takes its substrate (the part of the substrate of a
metabolism, the interstitium of a basolateral influx, the cell of an efflux,
the sink of an apical influx): the inhibitor species of the same part, else
of the equivalent tissue part, else plasma. With the competitive sum
`c = sum I/Ki` (competitive, mixed with its competitive constant,
noncompetitive, and mechanism-based inhibition with `Ki = KI`) and the
uncompetitive sum `u = sum I/Kiu` (uncompetitive, mixed, noncompetitive), the
rate laws of mutually exclusive inhibitors (Segel 1975) are

```
MM           Vmax*V*C / (Km*(1 + c) + C*(1 + u))
Hill         Vmax*V/(1 + u) * x^n / (1 + x^n),  x = C/(Km*(1 + c)/(1 + u))
first order  k*C / (1 + c), per volume k*V*C / (1 + c)
```

i.e. `Vmax_app = Vmax/(1 + u)` and `Km_app = Km*(1 + c)/(1 + u)`; the first
order rate is the limit `C << Km`, which an uncompetitive inhibitor does not
change.

A zonated organ (`pkpdsim.pbpk.variants.liver_zonated`) places every
metabolism and transport into the parts of each of its zones (`zone`), with
the ids of the zone (`Cti_z1_<s>`, `DRG2M_z1`) and the share of the zone of
the capacity of the organ (`capacity`); the apical sink is that of the organ.

A variant with a lumen of its own on the apical side of the cell gives it
(`Apical`): the nephron kidney (`pkpdsim.pbpk.variants.kidney_nephron`)
places the apical transport between the cells of the kidney and the lumen of
the proximal tubule, an efflux secretes into it, an influx reabsorbs from it
(with the inhibitors of its protein in the lumen).
"""

import copy
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING

from sbmlutils.factory import Parameter, Reaction, UnitDefinition
from sbmlutils.metadata import SBO

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.annotation import annotate
from pkpdsim.pbpk.chemistry import (
    MM,
    CustomReaction,
    Direction,
    FirstOrder,
    Hill,
    Membrane,
    Metabolism,
    Transport,
)
from pkpdsim.pbpk.data import cite
from pkpdsim.pbpk.impairment import organ_activity
from pkpdsim.pbpk.organs import Organ, Part
from pkpdsim.pbpk.parts import Parts
from pkpdsim.pbpk.proteins import Inhibition
from pkpdsim.pbpk.substance import KpMethod, Substance
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Q, as_value

if TYPE_CHECKING:
    from pkpdsim.pbpk.organ import OrganContext

__all__ = [
    "ONE",
    "Apical",
    "Inhibitors",
    "Sink",
    "apical_sink",
    "built_parts",
    "has_turnover",
    "in_zone",
    "inhibition_parameters",
    "inhibitor_terms",
    "inhibitor_unbound",
    "kp_id",
    "metabolism_part",
    "place",
    "process_source",
    "rate_law",
    "species_compartment",
    "stoichiometry",
    "unbound",
    "value_parameter",
    "with_modifiers",
]

_PART_NAMES = {
    Part.PLASMA: "plasma",
    Part.BLOOD_CELLS: "blood cells",
    Part.INTERSTITIAL: "interstitial",
    Part.INTRACELLULAR: "intracellular",
    Part.TISSUE: "tissue",
    Part.ORGAN: "organ",
}
"""The parts in names."""


ONE = "1 dimensionless"
"""The number one in a formula."""


_EQUIVALENT_PARTS: dict[Part, tuple[Part, ...]] = {
    Part.INTRACELLULAR: (Part.INTRACELLULAR, Part.TISSUE, Part.ORGAN, Part.PLASMA),
    Part.INTERSTITIAL: (Part.INTERSTITIAL, Part.TISSUE, Part.ORGAN, Part.PLASMA),
    Part.TISSUE: (Part.TISSUE, Part.INTRACELLULAR, Part.ORGAN, Part.PLASMA),
    Part.ORGAN: (Part.ORGAN, Part.TISSUE, Part.INTRACELLULAR, Part.PLASMA),
    Part.PLASMA: (Part.PLASMA, Part.ORGAN),
}
"""The part of another substance which corresponds to a part, in order: the
same part, the equivalent tissue part, plasma (the whole organ for a
flow-limited substance, which has no plasma species); where a product of a
metabolism goes and where an inhibitor acts."""


@dataclass(frozen=True)
class Sink:
    """The apical sink of an organ for a substance.

    Attributes:
        species: the id of the sink species (top model).
        compartment: the id of its compartment (top model).
        amount: the species is an amount (`hasOnlySubstanceUnits`), not a
            concentration.
    """

    species: str
    compartment: str
    amount: bool


@dataclass(frozen=True)
class Apical:
    """The apical side of the cell of a variant, instead of the sink of the organ.

    Attributes:
        name: the apical compartment in words, for the notes
            (`lumen of the proximal tubule`).
        species: the local id of the species of a substance (by its id) on
            the apical side; its concentration is unbound (a lumen without
            proteins).
    """

    name: str
    species: Callable[[str], str]


def apical_sink(organ: Organ, s: str) -> Sink:
    """Get the apical sink of an organ: bile, urine or the intestinal lumen.

    The lumen `Clumen_<s>` of a standard small intestine is an element of the
    body which the body has not: the intestine variants, which own the
    lumen, place the apical transport of the small intestine themselves, and
    the builder refuses it with the standard intestine.

    Args:
        organ: the organ.
        s: the substance id.

    Returns:
        The sink species of the substance and its compartment.

    Raises:
        ValueError: if the organ has no apical sink.
    """
    if organ is Organ.LIVER:
        return Sink(ids.bile(s), ids.VBILE, amount=False)
    if organ is Organ.KIDNEY:
        return Sink(ids.urine(s), ids.VURINE, amount=True)
    if organ is Organ.SMALL_INTESTINE:
        return Sink(ids.lumen(s), ids.VLUMEN, amount=False)
    raise ValueError(
        f"Organ '{organ.value}' has no apical sink: apical transport is defined "
        f"for the liver (bile), the kidney (urine) and the small intestine (lumen)."
    )


def built_parts(parts: Parts) -> tuple[Part, ...]:
    """Get the parts of an organ which hold a species of a substance.

    Args:
        parts: the parts of the substance in the organ.

    Returns:
        The parts with a species of the substance.
    """
    if parts.tissue == "organ":
        return (Part.ORGAN,)
    built = [Part.PLASMA]
    if parts.blood_cells:
        built.append(Part.BLOOD_CELLS)
    if parts.tissue == "merged":
        built.append(Part.TISSUE)
    elif parts.tissue == "full":
        built += [Part.INTERSTITIAL, Part.INTRACELLULAR]
    return tuple(built)


def metabolism_part(parts: Parts) -> Part:
    """Get the part in which a substance is metabolized: the part of its tissue.

    Args:
        parts: the parts of the substance in the organ.

    Returns:
        Intracellular (full form), the merged tissue, the whole organ
        (flow-limited) or plasma (no tissue).
    """
    return {
        "full": Part.INTRACELLULAR,
        "merged": Part.TISSUE,
        "organ": Part.ORGAN,
        "none": Part.PLASMA,
    }[parts.tissue]


def species_compartment(
    part: Part, s: str, parts: Parts, zone: str | None = None
) -> str:
    """Get the local compartment of the species of a substance in a part.

    The plasma species of a substance with erythrocytes in equilibrium lives
    in the vascular volume times the blood to plasma ratio (`Vpl_<s>`), so
    that it holds the amount of the blood; every other species lives in the
    compartment of its part.

    Args:
        part: the part.
        s: the substance id.
        parts: the parts of the substance in the organ.
        zone: the zone of a zonated organ, `None` for the organ.

    Returns:
        The local compartment id.
    """
    if part is Part.PLASMA and parts.flow_carrier == "blood_bp":
        return ids.plasma_capacity(s, zone)
    return ids.compartment(part, zone)


def kp_id(substance: Substance, organ: Organ) -> str | None:
    """Get the id of the whole tissue to plasma coefficient of a substance in an organ.

    Args:
        substance: the substance.
        organ: the organ.

    Returns:
        `Kp_<o>_<s>` for an organ with an override or a Kp method, `Kp_<s>`
        otherwise, `None` if the substance does not enter the tissue of the
        organ.
    """
    dist = substance.distribution
    if dist is None or (dist.organs is not None and organ not in dist.organs):
        return None
    if isinstance(dist, KpMethod) or organ in dist.overrides:
        return ids.kp_organ(organ, substance.sid)
    return ids.kp(substance.sid)


def unbound(
    part: Part, substance: Substance, organ: Organ, zone: str | None = None
) -> str:
    """Get the unbound concentration of a substance in a part, as formula.

    Args:
        part: the part.
        substance: the substance.
        organ: the organ.
        zone: the zone of a zonated organ, `None` for the organ.

    Returns:
        `fu*C` in plasma, `fu/K*C` with the coefficient of the part to plasma
        elsewhere (`K_int_pls`, `K_cell_pls`, `Kp`, `K_or`).

    Raises:
        ValueError: for the blood cells, or the merged tissue of a substance
            which does not enter the tissue of the organ.
    """
    s = substance.sid
    fu, c = ids.fu(s), ids.species(part, s, zone)
    if part is Part.PLASMA:
        return f"{fu} * {c}"
    if part is Part.INTERSTITIAL:
        return f"{fu} / {ids.k_int_pls(s)} * {c}"
    if part is Part.INTRACELLULAR:
        return f"{fu} / {ids.k_cell_pls(s)} * {c}"
    if part is Part.ORGAN:
        return f"{fu} / {ids.k_organ(s)} * {c}"
    kp = kp_id(substance, organ)
    if part is Part.TISSUE and kp is not None:
        return f"{fu} / {kp} * {c}"
    raise ValueError(
        f"{s}: no unbound concentration in part '{part.value}' of organ "
        f"'{organ.value}'."
    )


def process_source(
    process: Metabolism | Transport, parts: Mapping[str, Parts]
) -> Part | None:
    """Get the part from which a process takes its substrate.

    Args:
        process: the metabolism or transport.
        parts: the parts of the substances in the organ by substance id.

    Returns:
        The part of the substrate of a metabolism, the interstitium of a
        basolateral influx, the cell of an efflux, `None` for the sink of an
        apical influx.
    """
    if isinstance(process, Metabolism):
        return metabolism_part(parts[process.substrate])
    if process.direction is Direction.EFFLUX:
        return Part.INTRACELLULAR
    if process.membrane is Membrane.BASOLATERAL:
        return Part.INTERSTITIAL
    return None


def inhibitor_unbound(
    source: Part | None,
    inhibitor: Substance,
    organ: Organ,
    parts: Mapping[str, Parts],
    zone: str | None = None,
) -> str:
    """Get the unbound concentration of an inhibitor where a process takes its substrate.

    Args:
        source: the part of the substrate (`process_source`), `None` for the
            apical sink.
        inhibitor: the inhibiting (or inducing) substance.
        organ: the organ.
        parts: the parts of the substances in the organ by substance id.
        zone: the zone of a zonated organ, `None` for the organ.

    Returns:
        The unbound concentration (formula) of the inhibitor species of the
        same part, else of the equivalent tissue part, else of plasma; the
        concentration of the inhibitor in the sink of an apical influx.
    """
    return _inhibitor_site(source, inhibitor, organ, parts, zone)[1]


def _inhibitor_site(
    source: Part | None,
    inhibitor: Substance,
    organ: Organ,
    parts: Mapping[str, Parts],
    zone: str | None = None,
    apical: Apical | None = None,
) -> tuple[str, str]:
    """The species of an inhibitor where a process takes its substrate and its unbound concentration.

    `source` `None` is the apical side: the sink of the organ, or the
    apical compartment of a variant.
    """
    if source is None:
        if apical is not None:
            species = apical.species(inhibitor.sid)
        else:
            species = apical_sink(organ, inhibitor.sid).species
        return species, species
    built = built_parts(parts[inhibitor.sid])
    part = next(p for p in _EQUIVALENT_PARTS[source] if p in built)
    return (
        ids.species(part, inhibitor.sid, zone),
        unbound(part, inhibitor, organ, zone),
    )


def has_turnover(protein: str, ctx: "OrganContext") -> bool:
    """Check whether a protein has a turnover state.

    A protein has one when a substance of the definition inactivates it
    (mechanism-based) or induces it.

    Args:
        protein: the protein id.
        ctx: the context of the organ.

    Returns:
        `True` if the protein needs the state `E_<protein>`.
    """
    return any(
        i.protein == protein and i.kind == "mechanism" for _, i in ctx.inhibitions
    ) or any(i.protein == protein for _, i in ctx.inductions)


def _activity(protein: str | None, ctx: "OrganContext") -> str:
    """The factor of the protein on a rate: `f_<protein>`, times `E_<protein>`.

    In the liver times the factor of the impairment `fimp_<protein>` of a
    protein with data of impairment (`pkpdsim.pbpk.impairment`).
    """
    if protein is None:
        return ""
    factor = f"{ids.activity(protein)} * "
    impaired = organ_activity(ctx.organ, protein)
    if impaired is not None:
        factor += f"{impaired} * "
    if has_turnover(protein, ctx):
        factor += f"{ids.enzyme(protein)} * "
    return factor


@dataclass(frozen=True)
class Inhibitors:
    """The inhibitors of the protein of a process.

    Attributes:
        competitive: `(1 + sum I/Ki)` of the competitive part (competitive,
            mixed, noncompetitive, mechanism-based with `Ki = KI`), `None`
            without one.
        uncompetitive: `(1 + sum I/Kiu)` of the uncompetitive part
            (uncompetitive, mixed with `Kiu`, noncompetitive with
            `Kiu = Ki`), `None` without one.
        modifiers: the inhibitor species which the rate reads.
    """

    competitive: str | None = None
    uncompetitive: str | None = None
    modifiers: tuple[str, ...] = ()


def _inhibitors(
    process: Metabolism | Transport,
    ctx: "OrganContext",
    parts: Mapping[str, Parts],
    zone: str | None,
    apical: Apical | None = None,
) -> Inhibitors:
    """The competitive and the uncompetitive terms of the inhibitors of the protein."""
    source = process_source(process, parts) if process.protein is not None else None
    return inhibitor_terms(
        process,
        ctx,
        lambda s: _inhibitor_site(source, s, ctx.organ, parts, zone, apical),
    )


def inhibitor_terms(
    process: Metabolism | Transport,
    ctx: "OrganContext",
    site: Callable[[Substance], tuple[str, str]],
) -> Inhibitors:
    """Get the competitive and the uncompetitive terms of the inhibitors of a process.

    Args:
        process: the metabolism or transport.
        ctx: the context of the organ.
        site: the species of an inhibitor where the process takes its
            substrate and its unbound concentration (formula).

    Returns:
        The terms of the inhibitors of the protein of the process.
    """
    if process.protein is None:
        return Inhibitors()
    competitive: list[str] = []
    uncompetitive: list[str] = []
    modifiers: list[str] = []
    for sid, inh in ctx.inhibitions:
        if inh.protein != process.protein:
            continue
        species, i = site(_substance(ctx, sid))
        if species not in modifiers:
            modifiers.append(species)
        ki = ids.ki(inh.protein, sid)
        if inh.kind in ("competitive", "mixed", "noncompetitive", "mechanism"):
            competitive.append(f"({i}) / {ki}")
        if inh.kind in ("uncompetitive", "noncompetitive"):
            uncompetitive.append(f"({i}) / {ki}")
        elif inh.kind == "mixed":
            uncompetitive.append(f"({i}) / {ids.kiu(inh.protein, sid)}")
    return Inhibitors(
        competitive=f"({ONE} + {' + '.join(competitive)})" if competitive else None,
        uncompetitive=(
            f"({ONE} + {' + '.join(uncompetitive)})" if uncompetitive else None
        ),
        modifiers=tuple(modifiers),
    )


def inhibition_parameters(sid: str, inh: Inhibition) -> list[tuple[str, str]]:
    """Get the ids of the constants of an inhibition and their quantity.

    Args:
        sid: the inhibitor id.
        inh: the inhibition.

    Returns:
        (id, field of `Inhibition`) of `Ki`, and `Kiu` of a mixed, `kinact`
        of a mechanism-based inhibition.
    """
    pars = [(ids.ki(inh.protein, sid), "ki")]
    if inh.kind == "mixed":
        pars.append((ids.kiu(inh.protein, sid), "kiu"))
    if inh.kind == "mechanism":
        pars.append((ids.kinact(inh.protein, sid), "kinact"))
    return pars


def value_parameter(
    pid: str,
    q: Q,
    unit: UnitDefinition,
    name: str,
    sbo: str,
    text: str = "",
) -> Parameter:
    """Get a constant parameter of a number with its provenance.

    Args:
        pid: the id of the parameter.
        q: the number, plain (in `unit`) or a `Value` with unit and source.
        unit: the unit of a plain number.
        name: the name of the parameter.
        sbo: the SBO term.
        text: the description of the parameter.

    Returns:
        The parameter, the citation of the source in the notes and the
        uncertainty of a fitted value.
    """
    v = as_value(q, unit)
    lines = [text] if text else []
    if v.notes:
        lines.append(v.notes)
    if v.source == "fit":
        lines.append("Fitted value.")
    elif v.source is not None:
        lines.append(f"Source:\n\n- {cite(v.source)}")
    return Parameter(
        pid,
        v.value,
        v.unit,
        constant=True,
        name=name,
        sboTerm=sbo,
        notes="\n\n".join(lines) if lines else None,
        uncertainties=v.uncertainties,
    )


def place(
    process: Metabolism | Transport | CustomReaction,
    ctx: "OrganContext",
    parts: Mapping[str, Parts],
    capacity: str | None = None,
    zone: str | None = None,
    apical: Apical | None = None,
) -> list[Reaction]:
    """Place a process of the chemistry in the parts of the standard organ.

    The reactions reference the ids of the standard organ (species, `fu_<s>`,
    the coefficients of the parts, `V<o>`, `f_<protein>`, the sink species),
    which the organ creates; in a zone the ids of the zone.

    Args:
        process: the process, of the organ of `ctx`.
        ctx: the context of the organ.
        parts: the parts of the substances in the organ by substance id.
        capacity: the factor (formula) of the capacity of the organ which
            multiplies the rate of a metabolism or a transport: the
            functional fraction (`f_functional` of the liver), in a zone times
            the share of the zone, the tubular capacity of the kidney (with
            `f_renal_function`); `None` for none. A custom reaction states
            its rate itself and is not scaled.
        zone: the zone of a zonated organ in which the process is placed,
            `None` for the organ.
        apical: the apical side of the cell of a variant (the lumen of the
            proximal tubule of the nephron kidney), `None` for the sink of
            the organ.

    Returns:
        The reactions of the process.

    Raises:
        ValueError: if the process is of another organ, transports a substance
            without interstitial and intracellular species, is apical in an
            organ without sink, is an apical influx from the urine of the
            standard kidney, acts in a part which is not built, or is a
            custom reaction in a zone (its rate is that of the whole organ).
    """
    if process.organ is not ctx.organ:
        raise ValueError(
            f"Process '{process.sid}' is of organ '{process.organ.value}', not "
            f"of '{ctx.organ.value}'."
        )
    if isinstance(process, Metabolism):
        return [_metabolism(process, ctx, parts, capacity, zone)]
    if isinstance(process, Transport):
        return [_transport(process, ctx, parts, capacity, zone, apical)]
    if zone is not None:
        raise ValueError(
            f"Custom reaction '{process.sid}' cannot be placed in zone '{zone}': "
            "its rate is stated for the whole organ."
        )
    return [_custom(process, parts)]


def _substance(ctx: "OrganContext", sid: str) -> Substance:
    """Get a substance of the context by id."""
    for s in ctx.substances:
        if s.sid == sid:
            return s
    raise ValueError(f"Unknown substance '{sid}' in organ '{ctx.organ.value}'.")


def rate_law(
    process: Metabolism | Transport,
    ctx: "OrganContext",
    cu: str,
    capacity: str | None,
    inhibitors: Inhibitors,
) -> tuple[str, list[Parameter]]:
    """Get the rate law of a process on the unbound concentration `cu`.

    The activity of the protein (with its turnover state) and the factor of
    the capacity of the organ `capacity` multiply it; the inhibitors of the
    protein enter as the competitive and uncompetitive sums (module
    docstring). A first order clearance per volume is multiplied by the
    organ volume, as a Vmax.

    Args:
        process: the metabolism or transport.
        ctx: the context of the organ.
        cu: the unbound concentration (formula) the process acts on.
        capacity: the factor (formula) of the capacity, `None` for none.
        inhibitors: the terms of the inhibitors of the protein.

    Returns:
        The rate (mmole/min) and the parameters of the process.
    """
    pid, k = process.sid, process.kinetics
    kinds = _parameter_kinds(process)
    factor = _activity(process.protein, ctx)
    if capacity is not None:
        factor += f"{capacity} * "
    c, u = inhibitors.competitive, inhibitors.uncompetitive
    v = ids.organ_volume(ctx.organ)
    if isinstance(k, FirstOrder):
        p_k = f"{pid}_k"
        if k.per_volume:
            rate = f"{factor}{p_k} * {v} * {cu}"
            unit, name, text, kind = (
                U.per_min,
                f"clearance per volume {pid} [1/min]",
                "Clearance per litre of organ on the unbound concentration; the "
                f"organ has `{p_k} * {v}`.",
                kinds["k_volume"],
            )
        else:
            rate = f"{factor}{p_k} * {cu}"
            unit, name, text, kind = (
                U.l_per_min,
                f"clearance {pid} [l/min]",
                "Clearance of the whole organ on the unbound concentration.",
                kinds["k"],
            )
        if c:
            rate += f" / {c}"
        return rate, [
            annotate(
                value_parameter(p_k, k.k, unit, name, SBO.KINETIC_CONSTANT, text),
                kind,
                organ=ctx.organ,
                protein=process.protein,
            )
        ]
    vmax, km = f"{pid}_Vmax", f"{pid}_Km"
    pars = [
        annotate(
            value_parameter(
                vmax,
                k.vmax,
                U.mmole_per_min_l,
                f"Vmax {pid} [mmole/min/l]",
                SBO.MAXIMAL_VELOCITY,
                f"Maximal rate per litre of organ; the organ has `{vmax} * {v}`.",
            ),
            kinds["vmax"],
            organ=ctx.organ,
            protein=process.protein,
        ),
        annotate(
            value_parameter(
                km,
                k.km,
                U.mM,
                f"Km {pid} [mM]",
                SBO.MICHAELIS_CONSTANT,
                "On the unbound concentration.",
            ),
            kinds["km"],
            organ=ctx.organ,
            protein=process.protein,
        ),
    ]
    km_c = f"{km} * {c}" if c else km
    if isinstance(k, MM):
        s_u = f"({cu}) * {u}" if u else cu
        return f"{factor}{vmax} * {v} * ({cu}) / ({km_c} + {s_u})", pars
    assert isinstance(k, Hill)
    n = f"{pid}_n"
    pars.append(
        annotate(
            value_parameter(
                n, k.n, U.dimensionless, f"Hill coefficient {pid}", SBO.HILL_COEFFICIENT
            ),
            "hill",
            organ=ctx.organ,
            protein=process.protein,
        )
    )
    if u:
        vmax_app, km_app = f"{vmax} / {u}", f"({km_c} / {u})"
    else:
        vmax_app, km_app = vmax, f"({km_c})" if c else km
    # the power of the dimensionless ratio, so that the units of the rate are
    # checkable: `C^n/(Km^n + C^n) = x^n/(1 + x^n)` with `x = C/Km`
    x = f"(({cu}) / {km_app})"
    return (
        f"{factor}{vmax_app} * {v} * {x}^{n} / ({ONE} + {x}^{n})",
        pars,
    )


def _parameter_kinds(process: Metabolism | Transport) -> dict[str, str]:
    """The kinds (`pkpdsim.pbpk.annotation`) of the parameters of a process."""
    if isinstance(process, Metabolism):
        return {
            "k": "k_metabolism",
            "k_volume": "k_metabolism_volume",
            "vmax": "vmax",
            "km": "km",
        }
    side = "basolateral" if process.membrane is Membrane.BASOLATERAL else "apical"
    return {
        "k": "k_influx" if process.direction is Direction.INFLUX else "k_efflux",
        "k_volume": "k_transport_volume",
        "vmax": f"vmax_{side}",
        "km": f"km_{side}",
    }


def process_kind(process: Metabolism | Transport, organ: Organ) -> str:
    """Get the kind (`pkpdsim.pbpk.annotation`) of the reaction of a process.

    Args:
        process: the metabolism or transport.
        organ: the organ of the process.

    Returns:
        The kind: a metabolism (in the intestine the gut metabolism, else
        linear or saturable), or a transport (the apical efflux of the liver
        into the bile is the biliary excretion, of the kidney into the urine
        the renal secretion, the apical influx of the kidney the
        reabsorption and of the intestine the absorption).
    """
    if isinstance(process, Metabolism):
        if organ in {Organ.SMALL_INTESTINE, Organ.LARGE_INTESTINE}:
            return "gut_metabolism"
        return (
            "metabolism"
            if isinstance(process.kinetics, FirstOrder)
            else ("metabolism_saturable")
        )
    if process.membrane is Membrane.APICAL:
        efflux = process.direction is Direction.EFFLUX
        if organ is Organ.LIVER and efflux:
            return "biliary_excretion"
        if organ is Organ.KIDNEY:
            return "renal_secretion" if efflux else "reabsorption"
        if organ is Organ.SMALL_INTESTINE and not efflux:
            return "absorption"
    return "transport"


def with_modifiers(equation: str, modifiers: tuple[str, ...]) -> str:
    """An equation with the inhibitor species which its rate reads as modifiers.

    A species of the equation itself (an inhibitor of its own process) is not
    repeated.
    """
    sides = equation.replace("<->", "->").split("->")
    present = {
        part.strip().split(" ")[-1] for side in sides for part in side.split("+")
    }
    extra = [m for m in modifiers if m not in present]
    return f"{equation} [{', '.join(extra)}]" if extra else equation


def stoichiometry(fraction: float, species: str) -> str:
    """A product of an equation, the fraction left out if it is 1."""
    return species if fraction == 1.0 else f"{fraction!r} {species}"


def in_zone(text: str, zone: str | None) -> str:
    """Get a name or a note of an element, in a zone with the zone.

    Args:
        text: the name or note.
        zone: the zone of a zonated organ, `None` for the organ.

    Returns:
        The text, in a zone followed by `, zone <zone>`.
    """
    return text if zone is None else f"{text}, zone {zone}"


def _metabolism(
    process: Metabolism,
    ctx: "OrganContext",
    parts: Mapping[str, Parts],
    capacity: str | None,
    zone: str | None,
) -> Reaction:
    """The reaction of a metabolism in the part of the tissue of its substrate.

    A product goes into the first part of `_EQUIVALENT_PARTS` of the part of
    the substrate in which it has a species: the same part, else the
    equivalent tissue part (intracellular, merged tissue, whole organ), else
    plasma.
    """
    substrate = _substance(ctx, process.substrate)
    part = metabolism_part(parts[substrate.sid])
    cu = unbound(part, substrate, ctx.organ, zone)
    products = []
    for product, fraction in process.products.items():
        built = built_parts(parts[product])
        target = next(t for t in _EQUIVALENT_PARTS[part] if t in built)
        products.append(stoichiometry(fraction, ids.species(target, product, zone)))
    inhibitors = _inhibitors(process, ctx, parts, zone)
    rate, pars = rate_law(process, ctx, cu, capacity, inhibitors)
    names = ", ".join(process.products) or "elimination"
    reaction = Reaction(
        sid=ids.segment_process(process.sid, zone),
        name=in_zone(
            f"{process.substrate} -> {names}"
            + (f" ({process.protein})" if process.protein else ""),
            zone,
        ),
        equation=with_modifiers(
            f"{ids.species(part, substrate.sid, zone)} -> {' + '.join(products)}",
            inhibitors.modifiers,
        ),
        formula=(rate, U.mmole_per_min),
        pars=pars,
        compartment=species_compartment(
            part, substrate.sid, parts[substrate.sid], zone
        ),
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        notes=in_zone(f"Metabolism in the {_PART_NAMES[part]} part", zone)
        + f" on the unbound concentration `{cu}`.",
    )
    return annotate(
        reaction,
        process_kind(process, ctx.organ),
        organ=ctx.organ,
        protein=process.protein,
    )


def _transport(
    process: Transport,
    ctx: "OrganContext",
    parts: Mapping[str, Parts],
    capacity: str | None,
    zone: str | None,
    apical: Apical | None = None,
) -> Reaction:
    """The reaction of an active transport across the membrane of the cell."""
    substance = _substance(ctx, process.substance)
    s = substance.sid
    if parts[s].tissue != "full":
        raise ValueError(
            f"Transport '{process.sid}' needs the interstitial and intracellular "
            f"species of '{s}' in organ '{ctx.organ.value}', the parts are "
            f"'{parts[s].tissue}'."
        )
    cell = ids.species(Part.INTRACELLULAR, s, zone)
    if process.membrane is Membrane.BASOLATERAL:
        other = ids.species(Part.INTERSTITIAL, s, zone)
        outside = unbound(Part.INTERSTITIAL, substance, ctx.organ, zone)
        where = "interstitium"
    elif apical is not None:
        other = outside = apical.species(s)
        where = apical.name
    else:
        sink = apical_sink(ctx.organ, s)
        if process.direction is Direction.INFLUX and ctx.organ is Organ.KIDNEY:
            raise ValueError(
                f"Transport '{process.sid}': apical influx from the urine is not "
                "part of the standard kidney (the urine is cumulative); "
                "reabsorption is `Excretion.fraction_reabsorbed` or the nephron "
                "variant."
            )
        other = outside = sink.species
        where = {
            Organ.LIVER: "bile",
            Organ.KIDNEY: "urine",
            Organ.SMALL_INTESTINE: "lumen",
        }[ctx.organ]
    if process.direction is Direction.INFLUX:
        equation, cu = f"{other} -> {cell}", outside
    else:
        equation, cu = (
            f"{cell} -> {other}",
            unbound(Part.INTRACELLULAR, substance, ctx.organ, zone),
        )
    inhibitors = _inhibitors(process, ctx, parts, zone, apical)
    rate, pars = rate_law(process, ctx, cu, capacity, inhibitors)
    reaction = Reaction(
        sid=ids.segment_process(process.sid, zone),
        name=in_zone(
            f"{process.direction.value} {s} {process.membrane.value}"
            + (f" ({process.protein})" if process.protein else ""),
            zone,
        ),
        equation=with_modifiers(equation, inhibitors.modifiers),
        formula=(rate, U.mmole_per_min),
        pars=pars,
        compartment=ids.compartment(Part.INTRACELLULAR, zone),
        sboTerm=SBO.TRANSPORT_REACTION,
        notes=in_zone(
            f"{process.membrane.value.capitalize()} {process.direction.value} "
            f"between the cell and the {where}",
            zone,
        )
        + f" on the unbound concentration `{cu}`.",
    )
    return annotate(
        reaction,
        process_kind(process, ctx.organ),
        organ=ctx.organ,
        protein=process.protein,
    )


def _custom(process: CustomReaction, parts: Mapping[str, Parts]) -> Reaction:
    """The reaction of a custom process in the part it names."""
    if process.part not in _built(parts.values()):
        raise ValueError(
            f"Custom reaction '{process.sid}' acts in part '{process.part.value}' "
            f"of organ '{process.organ.value}', which no substance has there."
        )
    return annotate(
        Reaction(
            sid=process.sid,
            equation=process.equation,
            formula=(process.formula, U.mmole_per_min),
            pars=[
                annotate(copy.copy(par), "custom_parameter", organ=process.organ)
                for par in process.parameters
            ],
            compartment=ids.compartment(process.part),
            notes=f"Custom reaction in the {_PART_NAMES[process.part]} part.",
        ),
        "custom_reaction",
        organ=process.organ,
    )


def _built(parts: Iterable[Parts]) -> set[Part]:
    """The parts which hold a species of any substance."""
    return {part for p in parts for part in built_parts(p)}
