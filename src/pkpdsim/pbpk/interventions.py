"""The interventions on the lobes of the liver: ligation and resection.

An intervention is a set of changes of parameters of the model, as the
impairment maps of the library are (`pkpdsim.pbpk.impairment`): the models of
the library have no events, so a study applies the changes before a simulation
or between the phases of one, which is what a surgery at a time point is.
Three procedures:

- `PortalLigation` closes the portal branch of its lobes (`g_po_<lobe>` 0).
  The lobes keep their hepatic artery inflow, so they stay perfused below the
  reference and atrophy, and the lobes which keep their portal inflow take
  over the portal flow and hypertrophy. It is the portal vein ligation or
  embolisation which grows the future liver remnant before a hepatectomy.
- `ArterialLigation` closes the hepatic artery branch instead; both together
  devascularise a lobe.
- `Resection` removes its lobes: both branches close, `f_present_<lobe>` goes
  to 0, so the lobe has no capacity and no growth, and the remnant fraction of
  the liver falls to the share which is left.

An occlusion acts on the conductance of a branch and never on its resistance,
so no rate law divides by zero and a liver whose every branch is closed is
simply a liver without inflow.

The named procedures of the literature are the scenarios of
`pkpdsim.database.tables.liver_interventions`, one set per species:
`scenarios(species)` gives them as interventions and `scenario_changes(species)`
as the parameter maps which an experiment reads.
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.data import load
from pkpdsim.pbpk.lobes import lobe_table
from pkpdsim.species import Species

__all__ = [
    "ArterialLigation",
    "Intervention",
    "PortalLigation",
    "Resection",
    "scenario_changes",
    "scenarios",
]

_TABLE = "liver_interventions"


def _checked(lobes: Sequence[str], of_liver: Sequence[str], what: str) -> None:
    """Check that an intervention names lobes of the liver.

    Raises:
        ValueError: for a lobe which the liver does not have, naming the
            lobes it has.
    """
    unknown = sorted(set(lobes) - set(of_liver))
    if unknown:
        raise ValueError(
            f"The {what} names {unknown}, which are no lobes of the liver; its "
            f"lobes are {sorted(of_liver)}."
        )
    if not lobes:
        raise ValueError(f"The {what} names no lobe.")


@dataclass(frozen=True)
class PortalLigation:
    """The portal vein ligation of a set of lobes.

    Attributes:
        lobes: the tokens of the lobes whose portal branch is closed.
        name: the name of the intervention.
    """

    lobes: tuple[str, ...]
    name: str = "portal_ligation"

    def __init__(self, lobes: Sequence[str], name: str = "portal_ligation") -> None:
        """Take the lobes in the order they are given."""
        object.__setattr__(self, "lobes", tuple(lobes))
        object.__setattr__(self, "name", name)

    def changes(self, of_liver: Sequence[str]) -> dict[str, float]:
        """The parameter changes of the ligation.

        Args:
            of_liver: the lobe tokens of the liver of the model.

        Returns:
            The conductance of the portal branch of every ligated lobe, 0.

        Raises:
            ValueError: for a lobe which the liver does not have, or none.
        """
        _checked(self.lobes, of_liver, "portal vein ligation")
        return {ids.portal_conductance(lobe): 0.0 for lobe in self.lobes}


@dataclass(frozen=True)
class ArterialLigation:
    """The hepatic artery ligation of a set of lobes.

    Attributes:
        lobes: the tokens of the lobes whose hepatic artery branch is closed.
        name: the name of the intervention.
    """

    lobes: tuple[str, ...]
    name: str = "arterial_ligation"

    def __init__(self, lobes: Sequence[str], name: str = "arterial_ligation") -> None:
        """Take the lobes in the order they are given."""
        object.__setattr__(self, "lobes", tuple(lobes))
        object.__setattr__(self, "name", name)

    def changes(self, of_liver: Sequence[str]) -> dict[str, float]:
        """The parameter changes of the ligation.

        Args:
            of_liver: the lobe tokens of the liver of the model.

        Returns:
            The conductance of the arterial branch of every ligated lobe, 0.

        Raises:
            ValueError: for a lobe which the liver does not have, or none.
        """
        _checked(self.lobes, of_liver, "hepatic artery ligation")
        return {ids.arterial_conductance(lobe): 0.0 for lobe in self.lobes}


@dataclass(frozen=True)
class Resection:
    """The resection of a set of lobes: a hepatectomy.

    Both branches of the lobe close and the lobe is no longer in place, so it
    holds no capacity, does not grow and leaves the remnant fraction of the
    liver. What the lobe still held stays in it, cut off from the circulation;
    a study which wants the content of the resected tissue removed as well
    sets the species of the lobe to 0 in the same step.

    Attributes:
        lobes: the tokens of the lobes which are removed.
        name: the name of the intervention.
    """

    lobes: tuple[str, ...]
    name: str = "resection"

    def __init__(self, lobes: Sequence[str], name: str = "resection") -> None:
        """Take the lobes in the order they are given."""
        object.__setattr__(self, "lobes", tuple(lobes))
        object.__setattr__(self, "name", name)

    def changes(self, of_liver: Sequence[str]) -> dict[str, float]:
        """The parameter changes of the resection.

        Args:
            of_liver: the lobe tokens of the liver of the model.

        Returns:
            Both conductances and the presence of every removed lobe, 0.

        Raises:
            ValueError: for a lobe which the liver does not have, or none.
        """
        _checked(self.lobes, of_liver, "resection")
        changes: dict[str, float] = {}
        for lobe in self.lobes:
            changes[ids.arterial_conductance(lobe)] = 0.0
            changes[ids.portal_conductance(lobe)] = 0.0
            changes[ids.lobe_present(lobe)] = 0.0
        return changes


Intervention = PortalLigation | ArterialLigation | Resection
"""An intervention on the lobes of the liver."""

_KINDS: Mapping[
    str, type[PortalLigation] | type[ArterialLigation] | type[Resection]
] = {
    "portal_ligation": PortalLigation,
    "arterial_ligation": ArterialLigation,
    "resection": Resection,
}


def scenarios(species: Species) -> dict[str, Intervention]:
    """Get the named procedures of a species.

    Args:
        species: the species.

    Returns:
        The interventions by the name of the procedure
        (`pkpdsim.database.tables.liver_interventions`); a species without a
        procedure gets an empty mapping. The mapping refuses a name of
        another species, naming the names it has.

    Raises:
        ValueError: for a scenario whose kind is no intervention, or which
            names a lobe the species does not have.
    """
    table = load(_TABLE)
    entries = table.get(species.value, {})
    found: dict[str, Intervention] = {}
    tokens = [lobe.token for lobe in lobe_table(species)] if entries else []
    for name, entry in entries.items():
        kind = str(entry["kind"])
        if kind not in _KINDS:
            raise ValueError(
                f"Scenario '{name}' of {species.value} has the kind '{kind}', "
                f"which is no intervention; the kinds are {sorted(_KINDS)}."
            )
        lobes = [str(lobe) for lobe in entry["lobes"]]
        _checked(lobes, tokens, f"scenario '{name}' of {species.value}")
        found[name] = _KINDS[kind](lobes, name=name)
    return _Scenarios(species, found)


def scenario_changes(species: Species) -> dict[str, dict[str, float]]:
    """Get the parameter changes of the named procedures of a species.

    The map of an experiment (`PBPKSimulationExperiment.lobe_map`): the name
    of a procedure to the parameters it changes and their values.

    Args:
        species: the species.

    Returns:
        The changes by the name of the procedure.

    Raises:
        ValueError: whatever `scenarios` refuses.
    """
    tokens = [lobe.token for lobe in lobe_table(species)]
    return {
        name: intervention.changes(tokens)
        for name, intervention in scenarios(species).items()
    }


class _Scenarios(dict[str, Intervention]):
    """The procedures of a species, which refuse the name of another species."""

    def __init__(self, species: Species, found: Mapping[str, Intervention]) -> None:
        """Keep the species for the message of an unknown name."""
        super().__init__(found)
        self.species = species

    def __missing__(self, name: str) -> Intervention:
        """Refuse a name which the species does not have.

        Raises:
            ValueError: always, naming the procedures the species has.
        """
        raise ValueError(
            f"'{name}' is no procedure of the liver of the {self.species.value}; "
            f"its procedures are {sorted(self)}."
        )
