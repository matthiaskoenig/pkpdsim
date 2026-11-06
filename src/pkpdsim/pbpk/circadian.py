"""The circadian rhythms of the physiology: the changes over the day and the night.

A model with circadian rhythms (`PBPKDefinition.circadian`, `Circadian`)
knows the clock time: `CLOCK` is the clock time (hours) at the start of a
simulation, time 0, an initial value which an experiment sets (a dose in the
morning or in the evening, `PBPKSimulationExperiment.clock_changes`), and
`TOD` the clock time during the simulation (0 to 24 hours). Every rhythm of
`database/tables/circadian_human.toml` is a cosinor factor relative to the
mean over 24 hours,

    fcirc_<key> = 1 + A_circ_<key> * cos(2*pi*(CLOCK + time - phi_circ_<key>)/24 h),

with the relative amplitude `A_circ_<key>` and the acrophase
`phi_circ_<key>` (the clock time of the maximum) as constants of the top
model, so a fit or an experiment may change them. A factor is an assignment
rule of the time and multiplies its quantity, as the factors of the
impairment do (`pkpdsim.pbpk.impairment`): the GFR (`fimp_GFR`), the transit
time of the stomach in the fed state (`T_st_fed`, an input of the
intestine), the plasma albumin of the binding (`fimp_albumin`, not a default
rhythm) and the activity of a protein in an organ (`fimp_<o>_<protein>`).
Every quantity which reads a factor follows it during the simulation (an
assignment rule, as a parameter of `PBPKDefinition.dynamic`); without
`circadian` the model has none of these elements.
"""

import math
from collections.abc import Mapping
from dataclasses import dataclass
from functools import cache
from typing import TYPE_CHECKING, Any

from sbmlutils.factory import AssignmentRule, Parameter, Sbase
from sbmlutils.factory.units import UnitType
from sbmlutils.metadata import SBO

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.age import provenance
from pkpdsim.pbpk.annotation import annotate
from pkpdsim.pbpk.data import cite, load
from pkpdsim.pbpk.organs import Organ
from pkpdsim.pbpk.units import HOUR, U

if TYPE_CHECKING:
    from pkpdsim.pbpk.definition import PBPKDefinition

__all__ = [
    "CIRCADIAN",
    "CLOCK_DEFAULT",
    "PHYSIOLOGY_RHYTHMS",
    "Circadian",
    "Rhythm",
    "acrophase",
    "amplitude",
    "circadian_factor",
    "circadian_ids",
    "circadian_keys",
    "circadian_objects",
    "protein_rhythms",
    "rhythms",
]

CIRCADIAN = "circadian_human"
"""The table of the circadian rhythms in the database."""

CLOCK_DEFAULT = 8.0
"""The default clock time (hours) at the start of a simulation: 8:00 in the
morning, the dose of most studies."""

PHYSIOLOGY_RHYTHMS: Mapping[str, str] = {
    "gfr": ids.GFR,
    "gastric_emptying_fed": "t_stomach_fed",
    "albumin": "albumin",
}
"""The keys of the factors of the rhythms of the physiology by their entry of
the table: `fcirc_GFR` (on `fimp_GFR`), `fcirc_t_stomach_fed` (on the transit
time of the stomach in the fed state `T_st_fed`, an input of the intestine)
and `fcirc_albumin` (on `fimp_albumin`)."""

_SBO_PARAMETER = SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER


@dataclass(frozen=True)
class Rhythm:
    """A circadian rhythm, a cosinor of the clock time.

    Attributes:
        key: the key of its factor (`CO`, `li_cyp3a4`).
        amplitude: the amplitude relative to the mean over 24 hours.
        acrophase: the clock time of the maximum (hours).
        sources: the keys of the sources.
        notes: where the values come from.
        default: a rhythm of every model with `Circadian()`; a rhythm which
            is not is chosen by its key (`Circadian(rhythms=(...))`).
    """

    key: str
    amplitude: float
    acrophase: float
    sources: tuple[str, ...]
    notes: str
    default: bool = True

    def value(self, clock: float) -> float:
        """The factor at a clock time (hours)."""
        return 1.0 + self.amplitude * math.cos(
            2.0 * math.pi * (clock - self.acrophase) / 24.0
        )


def circadian_factor(key: str) -> str:
    """The id of the factor of a rhythm, `fcirc_CO`, `fcirc_li_cyp3a4`."""
    return f"fcirc_{key}"


def amplitude(key: str) -> str:
    """The id of the relative amplitude of a rhythm, `A_circ_CO`."""
    return f"A_circ_{key}"


def acrophase(key: str) -> str:
    """The id of the acrophase of a rhythm, `phi_circ_CO`."""
    return f"phi_circ_{key}"


def _rhythm(key: str, entry: Mapping[str, Any]) -> Rhythm:
    """The rhythm of an entry of the table.

    Raises:
        ValueError: for an amplitude not in [0, 1) or an acrophase not in
            [0, 24) or not in hours.
    """
    a = float(entry["amplitude"]["value"])
    phi = float(entry["acrophase"]["value"])
    if not 0.0 <= a < 1.0:
        raise ValueError(
            f"{CIRCADIAN}.toml: the amplitude {a} of '{key}' is not in [0, 1)."
        )
    if entry["acrophase"]["unit"] != "hr" or not 0.0 <= phi < 24.0:
        raise ValueError(
            f"{CIRCADIAN}.toml: the acrophase of '{key}' is no clock time in hours."
        )
    sources, notes = provenance(entry["amplitude"], entry["acrophase"])
    return Rhythm(key, a, phi, sources, notes, bool(entry.get("default", True)))


@cache
def rhythms() -> dict[str, Rhythm]:
    """The rhythms of the physiology of the table, by the key of their factor."""
    data = load(CIRCADIAN)
    return {
        key: _rhythm(key, data[name])
        for name, key in PHYSIOLOGY_RHYTHMS.items()
        if name in data
    }


@cache
def protein_rhythms() -> dict[tuple[Organ, str], Rhythm]:
    """The rhythms of the activity of the proteins, by organ and protein."""
    data = load(CIRCADIAN).get("protein", {})
    return {
        (Organ(o), p): _rhythm(f"{o}_{p.lower()}", entry)
        for o, proteins in data.items()
        for p, entry in proteins.items()
    }


@dataclass(frozen=True)
class Circadian:
    """The circadian rhythms of a model.

    Attributes:
        rhythms: the keys of the factors of the rhythms of the model (`GFR`,
            `t_stomach_fed`, `albumin`, `li_cyp3a4`); `None` for the default
            rhythms of the table (of the proteins those which the processes of
            the model use).
        clock: the default clock time (hours) at the start of a
            simulation, the value of `CLOCK`.

    Raises:
        ValueError: for an unknown rhythm or a clock time not in [0, 24).
    """

    rhythms: tuple[str, ...] | None = None
    clock: float = CLOCK_DEFAULT

    def __post_init__(self) -> None:
        """Check the rhythms and the clock time."""
        if not 0.0 <= self.clock < 24.0:
            raise ValueError(
                f"Circadian: the clock time {self.clock} is not in [0, 24)."
            )
        if self.rhythms is not None:
            known = set(rhythms()) | {r.key for r in protein_rhythms().values()}
            unknown = sorted(set(self.rhythms) - known)
            if unknown:
                raise ValueError(
                    f"Circadian: unknown rhythms {unknown}, the table has {sorted(known)}."
                )


def circadian_keys(definition: "PBPKDefinition") -> dict[str, Rhythm]:
    """The rhythms of a definition by the key of their factor.

    The rhythms of the physiology and of the proteins which a metabolism or
    a transport of the definition uses in their organ; none without
    `PBPKDefinition.circadian`.
    """
    # chemistry imports the definition of the organs, not this module
    from pkpdsim.pbpk.chemistry import Metabolism, Transport

    c = definition.circadian
    if c is None:
        return {}
    used = {
        (p.organ, p.protein)
        for p in definition.chemistry
        if isinstance(p, Metabolism | Transport) and p.protein is not None
    }
    found = dict(rhythms())
    found |= {r.key: r for pair, r in protein_rhythms().items() if pair in used}
    if c.rhythms is None:
        return {k: r for k, r in found.items() if r.default}
    return {k: r for k, r in found.items() if k in c.rhythms}


def circadian_ids(definition: "PBPKDefinition") -> tuple[str, ...]:
    """The ids of the factors of a definition, which change during a simulation."""
    return tuple(circadian_factor(k) for k in circadian_keys(definition))


def _parameter(
    pid: str,
    value: float | None,
    unit: UnitType,
    name: str,
    notes: str,
    kind: str,
    *,
    constant: bool = True,
) -> Parameter:
    """A parameter of the top model with its annotation."""
    return annotate(
        Parameter(
            pid,
            value,
            unit,
            constant=constant,
            name=name,
            sboTerm=_SBO_PARAMETER,
            notes=notes,
        ),
        kind,
    )


def circadian_objects(definition: "PBPKDefinition") -> list[Sbase]:
    """The clock, the time of day and the factors of the rhythms of a definition.

    Args:
        definition: the definition with `circadian`.

    Returns:
        `CLOCK`, the rule of `TOD`, the amplitude and acrophase of every
        rhythm and its factor, an assignment rule of the time; empty without
        `circadian`.
    """
    c = definition.circadian
    if c is None:
        return []
    # the clock time of the simulation in hours, a formula of the time
    clock = f"({ids.CLOCK} + time / 60 {U.min_per_hr.sid})"
    objects: list[Sbase] = [
        _parameter(
            ids.CLOCK,
            c.clock,
            HOUR,
            "clock time at the start [hr]",
            "The clock time at the start of the simulation, time 0 (hours, 0 to "
            "24): an initial value, which an experiment sets for a dose in the "
            "morning or in the evening (`pkpdsim.pbpk.circadian`).",
            "clock_time",
        ),
        _parameter(
            ids.TIME_OF_DAY,
            None,
            HOUR,
            "time of day [hr]",
            "The clock time during the simulation (hours, 0 to 24), `CLOCK` plus "
            "the time of the simulation, modulo 24 hours.",
            "clock_time",
            constant=False,
        ),
        AssignmentRule(
            ids.TIME_OF_DAY,
            f"{clock} - 24 {HOUR.sid} * floor({clock} / 24 {HOUR.sid})",
            HOUR,
            name="time of day",
        ),
    ]
    for key, r in circadian_keys(definition).items():
        a, phi, f = amplitude(key), acrophase(key), circadian_factor(key)
        citations = "\n".join(f"- {cite(s)}" for s in r.sources)
        objects += [
            _parameter(
                a,
                r.amplitude,
                U.dimensionless,
                f"circadian amplitude {key}",
                "The amplitude of the circadian rhythm relative to its mean over "
                f"24 hours: {r.notes}\n\n{citations}",
                "circadian",
            ),
            _parameter(
                phi,
                r.acrophase,
                HOUR,
                f"circadian acrophase {key} [hr]",
                f"The clock time of the maximum of the circadian rhythm: {r.notes}"
                f"\n\n{citations}",
                "circadian",
            ),
            _parameter(
                f,
                None,
                U.dimensionless,
                f"circadian factor {key}",
                f"The factor of the circadian rhythm on `{key}`, a cosinor of the "
                f"clock time relative to the mean over 24 hours: `1 + {a}*cos(2*pi*"
                f"(CLOCK + time - {phi})/24 hr)`.",
                "circadian",
                constant=False,
            ),
            AssignmentRule(
                f,
                f"1 dimensionless + {a} * cos(2 dimensionless * pi * "
                f"({clock} - {phi}) / 24 {HOUR.sid})",
                U.dimensionless,
                name=f"circadian factor {key}",
            ),
        ]
    return objects
