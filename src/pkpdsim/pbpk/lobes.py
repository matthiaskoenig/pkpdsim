"""The lobes of the liver: the sets of a species, their shares and the configuration.

The liver of a model is one organ with one hepatic artery inflow, one portal
vein inflow and one hepatic vein outflow, so everything it says is a statement
about the whole organ. A lobed liver splits it along a second spatial axis,
parallel to the zones along the sinusoid: one chain per lobe, each with its own
mass, its own arterial and portal inflow and its own drainage, so that the
perfusion of a deportalised lobe, the exposure of a lobe which lives on
arterial blood alone and the function left after a resection are quantities of
the model.

A `Lobe` carries the share of the liver mass of a lobe and, when a source
gives them, its shares of the arterial and of the portal inflow. No supply
fraction was found for any species, so a lobe without one is supplied in
proportion to its mass: a Couinaud segment is defined as a portal territory
and the hepatic artery branch runs with the portal branch in the portal
triad, which makes an equal perfusion per mass the defensible default. The
measured perfusion per lobe of a portal vein ligation study is the test of it.

The sets of the species are the tables
`database/tables/liver_lobes_{human,rat,mouse}.toml`: the anatomical lobes for
the rodents (the regions a lobe resolved study reports) and the eight Couinaud
segments for the human (the units of a planned resection or a portal vein
embolisation). The pig has a physiology table but no lobe data was sought, so
it has no lobes.

This module is configuration: it holds no SBML. `HABR` configures the hepatic
arterial buffer response, whose objects `pkpdsim.pbpk.cardiovascular` writes,
and `Remodelling` is the protocol of a driver of the lobe masses, whose
implementations are in `pkpdsim.pbpk.remodelling`.
"""

from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Protocol, runtime_checkable

from pkpdsim.database.tables import load
from pkpdsim.pbpk.data import value
from pkpdsim.sbml.values import Value
from pkpdsim.species import Species

if TYPE_CHECKING:
    from sbmlutils.factory import Sbase

__all__ = [
    "HABR",
    "LOBE_SPECIES",
    "Lobe",
    "Remodelling",
    "check_lobes",
    "lobe_table",
    "lobe_table_name",
]

LOBE_SPECIES = frozenset({Species.HUMAN, Species.RAT, Species.MOUSE})
"""The species with a lobe set; the pig has a physiology but no lobe data."""

TOLERANCE = 1e-6
"""How far the shares of a lobe set may miss 1."""


def lobe_table_name(species: Species) -> str:
    """The name of the table of the lobe set of a species.

    Args:
        species: the species.

    Returns:
        The name of the table without extension, `liver_lobes_<species>`.

    Raises:
        ValueError: for a species without a lobe set, naming it and the
            species which have one.
    """
    if species not in LOBE_SPECIES:
        have = ", ".join(sorted(s.value for s in LOBE_SPECIES))
        raise ValueError(
            f"Species '{species.value}' has no liver lobe set; the species with "
            f"a lobe set are {have}."
        )
    return f"liver_lobes_{species.value}"


@dataclass(frozen=True)
class Lobe:
    """A lobe of the liver with its share of the mass and of the blood supply.

    Attributes:
        token: the token of the lobe in the ids of the model, e.g. `lll` (the
            left lateral lobe of the rat) or `sg8` (Couinaud segment VIII). A
            lower case identifier; it joins the zone token of a zonated liver
            into the site `<lobe>_<zone>`.
        name: the name of the lobe, for the name of its elements.
        mass_fraction: the share of the liver mass of the lobe; the shares of
            a set sum to 1.
        arterial_fraction: the share of the hepatic artery inflow, `None` for
            the share of the mass.
        portal_fraction: the share of the portal vein inflow, `None` for the
            share of the mass.
        annotation: the term of an ontology for the lobe, `None` for a lobe
            whose term was not researched.
    """

    token: str
    name: str
    mass_fraction: Value
    arterial_fraction: Value | None = None
    portal_fraction: Value | None = None
    annotation: str | None = None

    def arterial_share(self) -> float:
        """The share of the hepatic artery inflow, the share of the mass without one."""
        if self.arterial_fraction is None:
            return self.mass_fraction.value
        return self.arterial_fraction.value

    def portal_share(self) -> float:
        """The share of the portal vein inflow, the share of the mass without one."""
        if self.portal_fraction is None:
            return self.mass_fraction.value
        return self.portal_fraction.value


def check_lobes(lobes: Sequence[Lobe]) -> tuple[Lobe, ...]:
    """Check a lobe set and return it.

    Args:
        lobes: the lobes.

    Returns:
        The lobes as a tuple, in the order they were given.

    Raises:
        ValueError: for an empty set, a token which is no lower case
            identifier, a repeated token, a negative share, mass fractions
            which do not sum to 1, or supply fractions which do not sum to 1
            when any lobe gives them.
    """
    if not lobes:
        raise ValueError("A liver needs at least one lobe.")
    seen: set[str] = set()
    for lobe in lobes:
        if not lobe.token.isidentifier():
            raise ValueError(f"Lobe '{lobe.token}': the token must be an identifier.")
        if lobe.token != lobe.token.lower():
            raise ValueError(f"Lobe '{lobe.token}': the token must be lower case.")
        if lobe.token in seen:
            raise ValueError(f"Lobe '{lobe.token}' is given twice.")
        seen.add(lobe.token)
        if lobe.mass_fraction.value < 0.0:
            raise ValueError(
                f"Lobe '{lobe.token}': the mass fraction "
                f"{lobe.mass_fraction.value} is negative."
            )
    _check_sum("mass", [lobe.mass_fraction.value for lobe in lobes])
    if any(lobe.arterial_fraction is not None for lobe in lobes):
        _check_sum("arterial", [lobe.arterial_share() for lobe in lobes])
    if any(lobe.portal_fraction is not None for lobe in lobes):
        _check_sum("portal", [lobe.portal_share() for lobe in lobes])
    return tuple(lobes)


def _check_sum(what: str, shares: Sequence[float]) -> None:
    """Check that the shares of a lobe set sum to 1.

    Args:
        what: what the shares are of, for the message.
        shares: the shares.

    Raises:
        ValueError: if the shares miss 1 by more than `TOLERANCE`.
    """
    total = sum(shares)
    if abs(total - 1.0) > TOLERANCE:
        raise ValueError(
            f"The {what} fractions of the liver lobes sum to {total}, not 1."
        )


def lobe_table(species: Species) -> tuple[Lobe, ...]:
    """Get the lobe set of a species from its table.

    Args:
        species: the species.

    Returns:
        The lobes in the order of the table.

    Raises:
        ValueError: for a species without a lobe set, or for a table whose
            shares are not a lobe set (`check_lobes`).
    """
    name = lobe_table_name(species)
    table = load(name)
    lobes = [
        Lobe(
            token=token,
            name=str(section["name"]),
            mass_fraction=value(name, token, "mass_fraction"),
            arterial_fraction=(
                value(name, token, "arterial_fraction")
                if "arterial_fraction" in section
                else None
            ),
            portal_fraction=(
                value(name, token, "portal_fraction")
                if "portal_fraction" in section
                else None
            ),
            annotation=section.get("annotation"),
        )
        for token, section in table.items()
    ]
    return check_lobes(lobes)


def _habr_gain() -> Value:
    """The measured gain of the hepatic arterial buffer response."""
    return value("cardiovascular", "defaults", "habr_gain")


@dataclass(frozen=True)
class HABR:
    """The hepatic arterial buffer response of a lobe.

    The hepatic artery of a territory whose portal inflow falls dilates, so a
    deportalised lobe is not simply ischaemic. The conductance of the arterial
    branch of a lobe rises with the share of its portal flow which is missing,
    `1 + gain*max(0, 1 - Qpo_<lobe>/Qpo_<lobe>_ref)`, capped at `maximum`.
    Without a buffer response (`None` in the lobed liver, the default) the
    arterial branches do not react, which is what a model without lobes does
    today.

    Attributes:
        gain: the fraction of a lost portal flow the artery compensates,
            `None` for the measured default of the table
            (`cardiovascular.toml`, 0.255 of Lautt 1985; Lautt 1990 reports
            0.44 with the passive pressure effect included).
        maximum: the largest factor the arterial conductance reaches, so that
            a fully deportalised lobe does not take an unbounded arterial
            flow.
    """

    gain: Value = field(default_factory=_habr_gain)
    maximum: float = 2.0

    def __post_init__(self) -> None:
        """Check the gain and the maximum.

        Raises:
            ValueError: for a negative gain or a maximum below 1.
        """
        if self.gain.value < 0.0:
            raise ValueError(
                f"The gain of the buffer response is {self.gain.value}, it must "
                "not be negative."
            )
        if self.maximum < 1.0:
            raise ValueError(
                f"The maximum of the buffer response is {self.maximum}, it must "
                "be at least 1 (no response)."
            )


@runtime_checkable
class Remodelling(Protocol):
    """A driver of the lobe masses after an intervention.

    A lobe which loses its portal inflow atrophies and its neighbours
    hypertrophy. A driver says how: `pkpdsim.pbpk.remodelling.FlowDriven`
    grows a lobe towards the mass at which its perfusion per mass is back at
    the target, which predicts a ligation pattern it was not fitted to;
    `PrescribedVolume` and `PrescribedFlow` follow a measured time course,
    which replays the data. Both write the same ids.
    """

    @property
    def name(self) -> str:
        """The name of the driver, for the notes of its elements."""
        ...

    def objects(self, lobes: Sequence[Lobe], species: Species) -> list[Sbase]:
        """Get the SBML objects which drive the masses of the lobes."""
        ...
