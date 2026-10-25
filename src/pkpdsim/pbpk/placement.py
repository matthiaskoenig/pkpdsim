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
capacity `Vmax * V<o>`, and the activity `f_<protein>` of its protein
multiplies it. The organ may give a functional fraction (`f_functional` of
the liver), which multiplies the rate of every metabolism and transport as
well, the share of the capacity in the functional tissue; a custom reaction
states its rate itself and is not scaled.
"""

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING

from sbmlutils.factory import Parameter, Reaction, UnitDefinition
from sbmlutils.metadata import SBO

from pkpdsim.pbpk import ids
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
from pkpdsim.pbpk.organs import Organ, Part
from pkpdsim.pbpk.parts import Parts
from pkpdsim.pbpk.substance import KpMethod, Substance
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Q, as_value

if TYPE_CHECKING:
    from pkpdsim.pbpk.organ import OrganContext

__all__ = [
    "Sink",
    "apical_sink",
    "built_parts",
    "kp_id",
    "metabolism_part",
    "place",
    "species_compartment",
    "unbound",
    "value_parameter",
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


_PRODUCT_PARTS: dict[Part, tuple[Part, ...]] = {
    Part.INTRACELLULAR: (Part.INTRACELLULAR, Part.TISSUE, Part.ORGAN, Part.PLASMA),
    Part.TISSUE: (Part.TISSUE, Part.INTRACELLULAR, Part.ORGAN, Part.PLASMA),
    Part.ORGAN: (Part.ORGAN, Part.TISSUE, Part.INTRACELLULAR, Part.PLASMA),
    Part.PLASMA: (Part.PLASMA, Part.ORGAN),
}
"""The parts a product of a metabolism goes into by the part of the substrate,
in order: the same part, the equivalent tissue part, plasma (the whole organ
for a flow-limited product, which has no plasma species)."""


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


def apical_sink(organ: Organ, s: str) -> Sink:
    """Get the apical sink of an organ: bile, urine or the intestinal lumen.

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


def species_compartment(part: Part, s: str, parts: Parts) -> str:
    """Get the local compartment of the species of a substance in a part.

    The plasma species of a substance with erythrocytes in equilibrium lives
    in the vascular volume times the blood to plasma ratio (`Vpl_<s>`), so
    that it holds the amount of the blood; every other species lives in the
    compartment of its part.

    Args:
        part: the part.
        s: the substance id.
        parts: the parts of the substance in the organ.

    Returns:
        The local compartment id.
    """
    if part is Part.PLASMA and parts.flow_carrier == "blood_bp":
        return ids.plasma_capacity(s)
    return ids.compartment(part)


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


def unbound(part: Part, substance: Substance, organ: Organ) -> str:
    """Get the unbound concentration of a substance in a part, as formula.

    Args:
        part: the part.
        substance: the substance.
        organ: the organ.

    Returns:
        `fu*C` in plasma, `fu/K*C` with the coefficient of the part to plasma
        elsewhere (`K_int_pls`, `K_cell_pls`, `Kp`, `K_or`).

    Raises:
        ValueError: for the blood cells, or the merged tissue of a substance
            which does not enter the tissue of the organ.
    """
    s = substance.sid
    fu, c = ids.fu(s), ids.species(part, s)
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
) -> list[Reaction]:
    """Place a process of the chemistry in the parts of the standard organ.

    The reactions reference the ids of the standard organ (species, `fu_<s>`,
    the coefficients of the parts, `V<o>`, `f_<protein>`, the sink species),
    which the organ creates.

    Args:
        process: the process, of the organ of `ctx`.
        ctx: the context of the organ.
        parts: the parts of the substances in the organ by substance id.
        capacity: the id of the functional fraction of the organ which
            multiplies the rate of a metabolism or a transport (`f_functional`
            of the liver), `None` for none; a custom reaction states its
            rate itself and is not scaled.

    Returns:
        The reactions of the process.

    Raises:
        ValueError: if the process is of another organ, transports a substance
            without interstitial and intracellular species, is apical in an
            organ without sink, or acts in a part which is not built.
    """
    if process.organ is not ctx.organ:
        raise ValueError(
            f"Process '{process.sid}' is of organ '{process.organ.value}', not "
            f"of '{ctx.organ.value}'."
        )
    if isinstance(process, Metabolism):
        return [_metabolism(process, ctx, parts, capacity)]
    if isinstance(process, Transport):
        return [_transport(process, ctx, parts, capacity)]
    return [_custom(process, parts)]


def _substance(ctx: "OrganContext", sid: str) -> Substance:
    """Get a substance of the context by id."""
    for s in ctx.substances:
        if s.sid == sid:
            return s
    raise ValueError(f"Unknown substance '{sid}' in organ '{ctx.organ.value}'.")


def _rate(
    process: Metabolism | Transport, organ: Organ, cu: str, capacity: str | None
) -> tuple[str, list[Parameter]]:
    """Get the rate law of a process on the unbound concentration `cu`.

    The activity of the protein and the functional fraction `capacity`
    multiply it.

    Returns:
        The rate (mmole/min) and the parameters of the process.
    """
    pid, k = process.sid, process.kinetics
    factor = f"{ids.activity(process.protein)} * " if process.protein else ""
    if capacity is not None:
        factor += f"{capacity} * "
    v = ids.organ_volume(organ)
    if isinstance(k, FirstOrder):
        p_k = f"{pid}_k"
        return f"{factor}{p_k} * {cu}", [
            value_parameter(
                p_k,
                k.k,
                U.l_per_min,
                f"clearance {pid} [l/min]",
                SBO.KINETIC_CONSTANT,
                "Clearance of the whole organ on the unbound concentration.",
            )
        ]
    vmax, km = f"{pid}_Vmax", f"{pid}_Km"
    pars = [
        value_parameter(
            vmax,
            k.vmax,
            U.mmole_per_min_l,
            f"Vmax {pid} [mmole/min/l]",
            SBO.MAXIMAL_VELOCITY,
            f"Maximal rate per litre of organ; the organ has `{vmax} * {v}`.",
        ),
        value_parameter(
            km,
            k.km,
            U.mM,
            f"Km {pid} [mM]",
            SBO.MICHAELIS_CONSTANT,
            "On the unbound concentration.",
        ),
    ]
    if isinstance(k, MM):
        return f"{factor}{vmax} * {v} * ({cu}) / ({km} + {cu})", pars
    assert isinstance(k, Hill)
    n = f"{pid}_n"
    pars.append(
        value_parameter(
            n, k.n, U.dimensionless, f"Hill coefficient {pid}", SBO.HILL_COEFFICIENT
        )
    )
    return (
        f"{factor}{vmax} * {v} * ({cu})^{n} / ({km}^{n} + ({cu})^{n})",
        pars,
    )


def _stoichiometry(fraction: float, species: str) -> str:
    """A product of an equation, the fraction left out if it is 1."""
    return species if fraction == 1.0 else f"{fraction!r} {species}"


def _metabolism(
    process: Metabolism,
    ctx: "OrganContext",
    parts: Mapping[str, Parts],
    capacity: str | None,
) -> Reaction:
    """The reaction of a metabolism in the part of the tissue of its substrate.

    A product goes into the first part of `_PRODUCT_PARTS` of the part of
    the substrate in which it has a species: the same part, else the
    equivalent tissue part (intracellular, merged tissue, whole organ), else
    plasma.
    """
    substrate = _substance(ctx, process.substrate)
    part = metabolism_part(parts[substrate.sid])
    cu = unbound(part, substrate, ctx.organ)
    products = []
    for product, fraction in process.products.items():
        built = built_parts(parts[product])
        target = next(t for t in _PRODUCT_PARTS[part] if t in built)
        products.append(_stoichiometry(fraction, ids.species(target, product)))
    rate, pars = _rate(process, ctx.organ, cu, capacity)
    names = ", ".join(process.products) or "elimination"
    return Reaction(
        sid=process.sid,
        name=f"{process.substrate} -> {names}"
        + (f" ({process.protein})" if process.protein else ""),
        equation=f"{ids.species(part, substrate.sid)} -> {' + '.join(products)}",
        formula=(rate, U.mmole_per_min),
        pars=pars,
        compartment=species_compartment(part, substrate.sid, parts[substrate.sid]),
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        notes=f"Metabolism in the {_PART_NAMES[part]} part on the unbound "
        f"concentration `{cu}`.",
    )


def _transport(
    process: Transport,
    ctx: "OrganContext",
    parts: Mapping[str, Parts],
    capacity: str | None,
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
    cell = ids.species(Part.INTRACELLULAR, s)
    if process.membrane is Membrane.BASOLATERAL:
        other = ids.species(Part.INTERSTITIAL, s)
        outside = unbound(Part.INTERSTITIAL, substance, ctx.organ)
        where = "interstitium"
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
            unbound(Part.INTRACELLULAR, substance, ctx.organ),
        )
    rate, pars = _rate(process, ctx.organ, cu, capacity)
    return Reaction(
        sid=process.sid,
        name=f"{process.direction.value} {s} {process.membrane.value}"
        + (f" ({process.protein})" if process.protein else ""),
        equation=equation,
        formula=(rate, U.mmole_per_min),
        pars=pars,
        compartment=ids.compartment(Part.INTRACELLULAR),
        sboTerm=SBO.TRANSPORT_REACTION,
        notes=f"{process.membrane.value.capitalize()} {process.direction.value} "
        f"between the cell and the {where} on the unbound concentration `{cu}`.",
    )


def _custom(process: CustomReaction, parts: Mapping[str, Parts]) -> Reaction:
    """The reaction of a custom process in the part it names."""
    if process.part not in _built(parts.values()):
        raise ValueError(
            f"Custom reaction '{process.sid}' acts in part '{process.part.value}' "
            f"of organ '{process.organ.value}', which no substance has there."
        )
    return Reaction(
        sid=process.sid,
        equation=process.equation,
        formula=(process.formula, U.mmole_per_min),
        pars=list(process.parameters),
        compartment=ids.compartment(process.part),
        notes=f"Custom reaction in the {_PART_NAMES[process.part]} part.",
    )


def _built(parts: Iterable[Parts]) -> set[Part]:
    """The parts which hold a species of any substance."""
    return {part for p in parts for part in built_parts(p)}
