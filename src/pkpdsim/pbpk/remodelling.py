"""The remodelling of the lobes of the liver: the anatomy follows the simulation.

A lobe which loses its portal inflow atrophies and the lobes which keep theirs
hypertrophy, over days: the volume of a lobe (`Vli_<lobe>`) is therefore a
state of the model and not a constant of the individual, and its share of the
liver `f_vol_<lobe>` follows it. The driver of a `LobedLiver` says what moves
it, and every driver moves it through the one quantity `dVli_<lobe>` (l/min),
the growth of the lobe, so the liver always knows how fast it is losing
tissue and what that tissue carries away.

Two kinds of driver, both of which a study needs:

- `FlowDriven` predicts the remodelling. The stimulus is the perfusion per
  volume of the lobe against the reference perfusion of the liver, which is
  what the hyperperfusion of the remnant after a resection and the
  hypoperfusion of a deportalised lobe after a ligation are: a lobe which is
  perfused above the reference grows, one below it shrinks, and the liver
  comes to rest when every lobe is at the reference again.
- `PrescribedVolume` and `PrescribedFlow` map a measured course onto the
  model: the volume of a lobe over time, as an imaging study reports it, or
  its growth rate.

The tissue a lobe loses carries its content away. The state of a species of
the liver is its amount, so a lobe which grows dilutes what it holds by
itself, which is what growing tissue does; a lobe which shrinks would
concentrate it instead, which no tissue does. The remodelling therefore emits
one flux per species and part of a shrinking lobe into the sink
`Aremodelled_<s>`, `A*max(0, -dVli_<lobe>)/Vli_<lobe>`, so that the
concentration of a shrinking lobe is left alone and the amount which left the
body is auditable.
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from sbmlutils.factory import (
    AssignmentRule,
    Compartment,
    InitialAssignment,
    Parameter,
    RateRule,
    Sbase,
    Species,
)
from sbmlutils.metadata import SBO

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.annotation import annotate
from pkpdsim.pbpk.data import value
from pkpdsim.pbpk.lobes import Lobe
from pkpdsim.pbpk.organs import Organ
from pkpdsim.pbpk.units import U
from pkpdsim.sbml.parameters import value_parameter
from pkpdsim.sbml.values import Formula, Value

if TYPE_CHECKING:
    from pkpdsim.pbpk.definition import PBPKDefinition

__all__ = [
    "REMODELLING_SINK",
    "VREMODELLED",
    "FlowDriven",
    "PrescribedFlow",
    "PrescribedVolume",
    "lobes_dynamic",
    "remodelling_sink",
]

REMODELLING_SINK = "Aremodelled"
"""The prefix of the sink of the tissue a lobe has lost, `Aremodelled_<s>`."""

VREMODELLED = "Vremodelled"
"""The nominal compartment of the sink of the remodelling."""

_SBO_PARAMETER = SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER
_TABLE = "liver_remodelling"


def remodelling_sink(s: str) -> str:
    """The id of the sink of the tissue a lobe has lost, `Aremodelled_<s>` (mmole)."""
    return f"{REMODELLING_SINK}_{s}"


def lobe_target(lobe: str) -> str:
    """The id of the prescribed volume a lobe follows, `Vli_target_<lobe>` (l)."""
    return f"Vli_target_{lobe}"


K_GROWTH = "k_growth"
"""The id of the rate constant of the flow driven growth of a lobe (1/min)."""

Q_TARGET = "q_li_target"
"""The id of the perfusion per volume a flow driven lobe grows towards (1/min)."""

TAU_TRACK = "tau_track"
"""The id of the tracking time of a prescribed lobe volume (min)."""


def _default(key: str) -> Value:
    """A default of the remodelling table."""
    return value(_TABLE, "defaults", key)


def _k_growth() -> Value:
    return _default(K_GROWTH)


def _tau_track() -> Value:
    return _default(TAU_TRACK)


def lobes_dynamic(definition: PBPKDefinition) -> tuple[str, ...]:
    """The lobe volumes which change during a simulation.

    The quantities of the top model and of the organs which read them (the
    share of a lobe, the compartments of its sites, its perfusion) must
    follow them, so they are the seeds of `pkpdsim.pbpk.builder._dynamic`.

    Args:
        definition: the definition.

    Returns:
        The ids of the lobe volumes, empty for a liver without lobes or
        without a driver of the remodelling.

    Raises:
        ValueError: for a species without a lobe set.
    """
    from pkpdsim.pbpk.builder import organ_component
    from pkpdsim.pbpk.variants.liver_lobes import LobedLiver

    liver = organ_component(definition, Organ.LIVER)
    if not isinstance(liver, LobedLiver) or liver.remodelling is None:
        return ()
    lobes = liver.lobe_set(definition.species)
    if len(lobes) < 2:
        return ()
    return tuple(ids.lobe_volume(lobe.token) for lobe in lobes)


def _parameter(pid: str, v: float, unit: object, name: str, text: str) -> Parameter:
    """A variable parameter of the remodelling, which a rule follows."""
    return Parameter(
        pid,
        v,
        unit,  # ty: ignore[invalid-argument-type]
        constant=False,
        name=name,
        sboTerm=_SBO_PARAMETER,
        notes=text,
    )


def _growth(lobe: Lobe, formula: str) -> list[Sbase]:
    """The growth of a lobe and the rate rule of its volume.

    Every driver emits these two: the named growth `dVli_<lobe>`, which the
    liver reads for the tissue a shrinking lobe carries away, and the rate
    rule which moves the volume by it. A lobe which is not in place
    (`f_present_<lobe>` 0, resected) does not move.
    """
    growth = ids.lobe_growth(lobe.token)
    return [
        annotate(
            _parameter(
                growth,
                0.0,
                U.l_per_min,
                f"growth of {lobe.name} [l/min]",
                "Rate of change of the volume of the lobe: positive is "
                "hypertrophy, negative atrophy. The liver reads it for the "
                "amount which atrophied tissue carries away.",
            ),
            "lobe_flow",
            organ=Organ.LIVER,
        ),
        AssignmentRule(
            growth, f"{ids.lobe_present(lobe.token)} * ({formula})", U.l_per_min
        ),
        RateRule(ids.lobe_volume(lobe.token), growth, U.l_per_min),
    ]


@dataclass(frozen=True)
class FlowDriven:
    """The predictive driver: a lobe grows towards the reference perfusion per volume.

    `dVli_<lobe> = f_present_<lobe>*k_growth*Vli_<lobe>*(q_<lobe>/q_target - 1)`,
    the stimulus being the perfusion per volume of the lobe against the
    reference perfusion of the liver. The remnant after a resection carries
    the whole hepatic inflow in a fraction of the volume, so it is perfused
    above the reference and grows; a lobe whose portal branch is ligated
    keeps its arterial inflow only, so it is perfused below the reference and
    atrophies. The liver comes to rest when every lobe is at the reference
    perfusion again, which is the intact liver, so the law needs no ceiling
    of its own; `total_cap` adds one for a study which wants a hard limit.

    Attributes:
        k: the rate constant of the growth (1/min), the relative growth of a
            lobe at twice the reference perfusion. A model choice to be
            fitted, not a measured constant
            (`pkpdsim.database.tables.liver_remodelling`).
        target: the perfusion per volume a lobe grows towards (1/min),
            `None` for the reference perfusion of the liver of the
            individual, `Qh/Vli`.
        total_cap: the growth of a lobe is scaled by
            `max(0, 1 - Vli_lobes/Vli_ref)`, so the liver cannot grow past
            the volume of the intact organ whatever the perfusion is. It is
            off by default: the perfusion stimulus already brings the liver
            to rest at the intact volume, and the cap slows the whole
            approach to it.
        name: the name of the driver.
    """

    k: Value = field(default_factory=_k_growth)
    target: Value | None = None
    total_cap: bool = False
    name: str = "flow_driven"

    def __post_init__(self) -> None:
        """Check the constants.

        Raises:
            ValueError: for a negative rate constant or a target which is
                not positive.
        """
        if self.k.value < 0.0:
            raise ValueError(
                f"The growth constant of the remodelling must be >= 0, not {self.k.value}."
            )
        if self.target is not None and self.target.value <= 0.0:
            raise ValueError(
                f"The perfusion target of the remodelling must be > 0, not "
                f"{self.target.value}."
            )

    def objects(self, lobes: Sequence[Lobe], species: object) -> list[Sbase]:
        """Get the constants, the growth and the rate rules of the flow driven growth."""
        objects: list[Sbase] = [
            annotate(
                value_parameter(
                    K_GROWTH,
                    self.k,
                    U.per_min,
                    "growth constant of the remodelling [1/min]",
                    _SBO_PARAMETER,
                    "Rate constant of the flow driven growth of a lobe: the "
                    "relative growth at twice the reference perfusion per "
                    "volume.",
                ),
                "lobe_perfusion",
                organ=Organ.LIVER,
            ),
        ]
        if self.target is None:
            objects += [
                annotate(
                    Parameter(
                        Q_TARGET,
                        1.0,
                        U.per_min,
                        constant=True,
                        name="perfusion target of the remodelling [1/min]",
                        sboTerm=_SBO_PARAMETER,
                        notes="The perfusion per volume a lobe grows towards: "
                        "the reference perfusion of the liver of the "
                        "individual, the total hepatic blood flow per liver "
                        "volume.",
                    ),
                    "lobe_perfusion",
                    organ=Organ.LIVER,
                ),
                InitialAssignment(
                    Q_TARGET,
                    f"{ids.QH} / {ids.organ_volume(Organ.LIVER)}",
                    U.per_min,
                ),
            ]
        else:
            objects.append(
                annotate(
                    value_parameter(
                        Q_TARGET,
                        self.target,
                        U.per_min,
                        "perfusion target of the remodelling [1/min]",
                        _SBO_PARAMETER,
                        "The perfusion per volume a lobe grows towards.",
                    ),
                    "lobe_perfusion",
                    organ=Organ.LIVER,
                )
            )
        cap = (
            f" * max(0 dimensionless, 1 dimensionless - "
            f"{ids.VLI_LOBES} / {ids.VLI_REF})"
            if self.total_cap
            else ""
        )
        for lobe in lobes:
            excess = (
                f"({ids.lobe_perfusion(lobe.token)} / {Q_TARGET} - 1 dimensionless)"
            )
            objects += _growth(
                lobe,
                f"{K_GROWTH} * {ids.lobe_volume(lobe.token)} * {excess}{cap}",
            )
        return objects


def _curves(
    curves: Mapping[str, Formula], lobes: Sequence[Lobe], what: str
) -> dict[str, Formula]:
    """Check that every prescribed curve names a lobe of the liver.

    Raises:
        ValueError: for a key which is no lobe of the liver, naming the
            lobes it has.
    """
    tokens = {lobe.token for lobe in lobes}
    unknown = sorted(set(curves) - tokens)
    if unknown:
        raise ValueError(
            f"The prescribed {what} of the remodelling names {unknown}, which "
            f"are no lobes of the liver; its lobes are {sorted(tokens)}."
        )
    return dict(curves)


@dataclass(frozen=True)
class PrescribedVolume:
    """The data driven driver: a lobe follows a measured volume over time.

    The volume of a lobe is a state which tracks the given curve with the
    time constant `tau_track`, `dVli_<lobe> = (Vli_target_<lobe> -
    Vli_<lobe>)/tau_track`, and not an assignment rule of the curve itself:
    the growth of a lobe has to be a quantity of the model, since the liver
    reads it for the tissue a shrinking lobe carries away, and the rate of
    change of an assignment rule is not available to a rate law. One hour is
    fast against the days of liver remodelling, and the state lags the curve
    by the tracking time times its slope.

    Attributes:
        curves: the volume of a lobe over time (l), by lobe token; a lobe
            without a curve keeps its volume. The formula reads the constants
            of the top model and `time`, as every `Formula` of the library
            does (`Vli_ref` is the volume of the intact liver).
        tau: the tracking time (min).
        name: the name of the driver.
    """

    curves: Mapping[str, Formula] = field(default_factory=dict)
    tau: Value = field(default_factory=_tau_track)
    name: str = "prescribed_volume"

    def __post_init__(self) -> None:
        """Check the tracking time.

        Raises:
            ValueError: for a tracking time which is not positive.
        """
        if self.tau.value <= 0.0:
            raise ValueError(
                f"The tracking time of the remodelling must be > 0, not {self.tau.value}."
            )

    def objects(self, lobes: Sequence[Lobe], species: object) -> list[Sbase]:
        """Get the targets, the tracking time, the growth and the rate rules.

        Raises:
            ValueError: for a curve which names no lobe of the liver.
        """
        curves = _curves(self.curves, lobes, "volumes")
        objects: list[Sbase] = [
            annotate(
                value_parameter(
                    TAU_TRACK,
                    self.tau,
                    U.min,
                    "tracking time of the prescribed volume [min]",
                    _SBO_PARAMETER,
                    "The time constant with which the volume of a lobe "
                    "follows its prescribed course.",
                ),
                "lobe_perfusion",
                organ=Organ.LIVER,
            )
        ]
        for lobe in lobes:
            curve = curves.get(lobe.token)
            if curve is None:
                objects += _growth(lobe, "0 l_per_min")
                continue
            target = lobe_target(lobe.token)
            objects += [
                annotate(
                    _parameter(
                        target,
                        1.0,
                        U.litre,
                        f"prescribed volume of {lobe.name} [l]",
                        "The volume of the lobe over time, the course which "
                        "the model is given.",
                    ),
                    "lobe_volume",
                    organ=Organ.LIVER,
                ),
                AssignmentRule(target, curve.expression, U.litre),
            ]
            objects += _growth(
                lobe, f"({target} - {ids.lobe_volume(lobe.token)}) / {TAU_TRACK}"
            )
        return objects


@dataclass(frozen=True)
class PrescribedFlow:
    """The data driven driver: a lobe grows at a measured rate over time.

    `dVli_<lobe>` is the given curve, so the volume of the lobe is its
    integral. It is the driver of a study which reports the growth rate of a
    lobe instead of its volume, and the exact form of `PrescribedVolume`,
    which has to track its course.

    Attributes:
        curves: the growth of a lobe over time (l/min), by lobe token; a lobe
            without a curve keeps its volume.
        name: the name of the driver.
    """

    curves: Mapping[str, Formula] = field(default_factory=dict)
    name: str = "prescribed_flow"

    def objects(self, lobes: Sequence[Lobe], species: object) -> list[Sbase]:
        """Get the growth and the rate rules of the prescribed growth.

        Raises:
            ValueError: for a curve which names no lobe of the liver.
        """
        curves = _curves(self.curves, lobes, "flows")
        objects: list[Sbase] = []
        for lobe in lobes:
            curve = curves.get(lobe.token)
            objects += _growth(
                lobe, curve.expression if curve is not None else "0 l_per_min"
            )
        return objects


def sink_objects(definition: PBPKDefinition) -> list[Sbase]:
    """Get the sink of the remodelling: the tissue the lobes have lost.

    Args:
        definition: the definition.

    Returns:
        The nominal compartment and one amount species per substance.
    """
    objects: list[Sbase] = [
        annotate(
            Compartment(
                VREMODELLED,
                1.0,
                U.litre,
                constant=True,
                name="remodelled tissue",
                notes="Nominal volume of a compartment of amounts, which "
                "enters no rate: the tissue the lobes of the liver have lost.",
            ),
            "remodelled_tissue",
        )
    ]
    for s in definition.substances:
        objects.append(
            annotate(
                Species(
                    remodelling_sink(s.sid),
                    VREMODELLED,
                    initialAmount=0.0,
                    substanceUnit=U.mmole,
                    hasOnlySubstanceUnits=True,
                    name=f"{s.name} (remodelled tissue)",
                    sboTerm=SBO.SIMPLE_CHEMICAL,
                    notes=f"Cumulative amount of {s.sid} which the tissue of "
                    "an atrophying lobe of the liver carried away.",
                ),
                "amount_remodelled",
                substance=s,
            )
        )
    return objects
