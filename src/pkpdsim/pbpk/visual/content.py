"""The content of a definition which the figures show.

Read once from a `PBPKDefinition` through the readers of the schematic
(`pkpdsim.pbpk.figure`), so the anatomical figure, its interactive page and the
schematic show the same model: the substances with their colours and
mechanisms, per organ the chips (the substance enters the tissue or stays in
the plasma) and the process lines, the parts of every substance in every organ
(`parts_for`), the processes as records, the routes, the lumen, the bile, the
urine and the PD components.
"""

from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType

from pkpdsim.pbpk.chemistry import Metabolism, Transport
from pkpdsim.pbpk.definition import PBPKDefinition
from pkpdsim.pbpk.figure import _Figure
from pkpdsim.pbpk.organs import Organ
from pkpdsim.pbpk.parts import parts_for

__all__ = ["Chip", "ModelContent", "OrganContent", "Process", "SubstanceInfo"]


@dataclass(frozen=True)
class SubstanceInfo:
    """A substance of the definition: its id, name, colour and mechanisms in words."""

    sid: str
    name: str
    color: str
    details: str


@dataclass(frozen=True)
class Chip:
    """A substance in an organ: its id, colour and whether it enters the tissue."""

    sid: str
    color: str
    tissue: bool


@dataclass(frozen=True)
class OrganContent:
    """An organ of the definition.

    Attributes:
        key: the `Organ` value.
        chips: the substances of the organ.
        lines: the processes as (badge, text) of the schematic.
        variant: the variant of the organ in words, `None` for the standard organ.
        parts: the tissue form of every substance (`none`, `organ`, `merged`,
            `full`, see `pkpdsim.pbpk.parts.Parts`).
    """

    key: str
    chips: tuple[Chip, ...]
    lines: tuple[tuple[str, str], ...]
    variant: str | None
    parts: Mapping[str, str]


@dataclass(frozen=True)
class Process:
    """A metabolism or transport of an organ (a custom reaction has `kind` `CUS`).

    Attributes:
        kind: `MET`, `TRA` or `CUS`.
        sid: the id of the process.
        organ: the `Organ` value.
        substrate: the substrate of a metabolism, the substance of a transport.
        products: the products of a metabolism.
        protein: the enzyme or transporter.
        membrane: `basolateral` or `apical` of a transport.
        direction: `influx` or `efflux` of a transport.
    """

    kind: str
    sid: str
    organ: str
    substrate: str = ""
    products: tuple[str, ...] = ()
    protein: str | None = None
    membrane: str = ""
    direction: str = ""


@dataclass(frozen=True)
class ModelContent:
    """What the figures of a definition show.

    Attributes:
        sid: the id of the model.
        name: the name of the model.
        substances: the substances in the order of the definition.
        organs: the organs by their key.
        processes: the processes of the chemistry.
        lumen_lines: the process lines of the lumen, empty without a lumen.
        routes: the routes as (kind, substance), kind `iv`, `oral` or
            `oral_solution`.
        route_text: the routes in words.
        luminal: the substances which reach the lumen.
        bile: the substances excreted into the bile.
        urine: the substances excreted into the urine.
        pd: the PD components as (id, reads, sets).
    """

    sid: str
    name: str
    substances: tuple[SubstanceInfo, ...]
    organs: Mapping[str, OrganContent]
    processes: tuple[Process, ...]
    lumen_lines: tuple[tuple[str, str], ...]
    routes: tuple[tuple[str, str], ...]
    route_text: str
    luminal: tuple[str, ...]
    bile: tuple[str, ...]
    urine: tuple[str, ...]
    pd: tuple[tuple[str, tuple[str, ...], tuple[str, ...]], ...]

    @classmethod
    def of(cls, definition: PBPKDefinition) -> "ModelContent":
        """Read the content of a definition through the readers of the schematic."""
        f = _Figure(definition, (definition.sex,))
        organs = {}
        for o in Organ:
            organs[o.value] = OrganContent(
                key=o.value,
                chips=tuple(Chip(s, c, t) for s, c, t in f._chips(o)),
                lines=tuple(f._organ_lines([o])),
                variant=f._variant(o),
                parts=MappingProxyType(
                    {
                        s.sid: parts_for(s, o, definition.chemistry).tissue
                        for s in definition.substances
                    }
                ),
            )
        processes = []
        for p in definition.chemistry:
            if isinstance(p, Metabolism):
                processes.append(
                    Process(
                        "MET",
                        p.sid,
                        Organ(p.organ).value,
                        p.substrate,
                        tuple(p.products),
                        p.protein,
                    )
                )
            elif isinstance(p, Transport):
                processes.append(
                    Process(
                        "TRA",
                        p.sid,
                        Organ(p.organ).value,
                        p.substance,
                        (),
                        p.protein,
                        str(p.membrane),
                        str(p.direction),
                    )
                )
            else:
                processes.append(Process("CUS", p.sid, Organ(p.organ).value))
        lumen = f._lumen_row(f.gut).lines if (f.gut is not None and f.luminal) else []
        return cls(
            sid=definition.sid,
            name=definition.name or definition.sid,
            substances=tuple(
                SubstanceInfo(s.sid, s.name, f.colors[s.sid], f._details(s))
                for s in definition.substances
            ),
            organs=MappingProxyType(organs),
            processes=tuple(processes),
            lumen_lines=tuple(lumen),
            routes=tuple((r.kind, r.substance) for r in definition.routes),
            route_text=f._route_text(),
            luminal=tuple(f.luminal),
            bile=tuple(f.bile),
            urine=tuple(s.sid for s in definition.substances if s.excretion.urine),
            pd=tuple(
                (pd.model.sid, tuple(pd.inputs.values()), tuple(pd.outputs.values()))
                for pd in definition.pd
            ),
        )

    def color(self, sid: str) -> str:
        """The colour of a substance."""
        return next(s.color for s in self.substances if s.sid == sid)

    def has_route(self, kind: str) -> bool:
        """Whether a route of a kind (`iv`, `oral`) exists; `oral` includes the solutions."""
        if kind == "oral":
            return any(k.startswith("oral") for k, _ in self.routes)
        return any(k == kind for k, _ in self.routes)
