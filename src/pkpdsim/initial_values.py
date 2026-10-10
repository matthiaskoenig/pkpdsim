"""The initial values of a substance in every tissue of a model.

A substance which is transported in the body is initialized in every tissue
to the same value, which removes the distribution kinetics: the simulation
starts in (or close to) the steady state of the distribution. The functions
return the changes of a simulation, `[sid]` for the concentration and `sid`
for the amount of a species, and leave the urine alone.

The species of a substance are found by the naming pattern of the ids of the
models, `C<tissue>_<substance>` for a concentration and `A<tissue>_<substance>`
for an amount (`Cve_len`, `Cli_plasma_len`, `Aurine_len`), see
`get_species_keys`.
"""

import logging
import re

import roadrunner

logger = logging.getLogger(__name__)

# -------------------------------------------------------------------------------------------------
# Initial values
# -------------------------------------------------------------------------------------------------
# Helper functions for setting initial values in the model.
# Substances which are transported in the body can hereby be initialized in all
# tissues to identical values (which removes the distribution kinetics).
# -------------------------------------------------------------------------------------------------


def init_concentrations_changes(
    r: roadrunner.RoadRunner, skey: str, value: float
) -> dict[str, float]:
    """Get changes to set initial concentrations for skey.

    Args:
        r: roadrunner model.
        skey: substance key, e.g., `len` for the species `Cve_len`.
        value: initial concentration in model units.

    Returns:
        Changes for the model, species selection to value.
    """
    return _set_initial_values(r, skey, value, method="concentration")


def init_amounts_changes(
    r: roadrunner.RoadRunner, skey: str, value: float
) -> dict[str, float]:
    """Get changes to set initial amounts for skey.

    Args:
        r: roadrunner model.
        skey: substance key, e.g., `len` for the species `Cve_len`.
        value: initial amount in model units.

    Returns:
        Changes for the model, species selection to value.
    """
    return _set_initial_values(r, skey, value, method="amount")


def _set_initial_values(
    r: roadrunner.RoadRunner, skey: str, value: float, method: str = "concentration"
) -> dict[str, float]:
    """Set the initial value of a distributing substance.

    Takes care of all the compartment values so starting close/in steady state.
    Units are in model units.

    Args:
        r: roadrunner model.
        skey: substance key.
        value: initial value in model units.
        method: `concentration` or `amount`.

    Returns:
        Changes for the model, species selection to value.

    Raises:
        ValueError: for an unknown method.
    """
    if method not in ["amount", "concentration"]:
        raise ValueError(f"method must be 'amount' or 'concentration': '{method}'")

    species_ids = r.model.getFloatingSpeciesIds() + r.model.getBoundarySpeciesIds()
    changeset: dict[str, float] = {}

    for key in get_species_keys(skey, species_ids):
        if "urine" in key:
            logger.debug("urinary values are not set")
            continue

        # `[sid]` selects the concentration, `sid` the amount
        rkey = f"[{key}]" if method == "concentration" else key
        changeset[rkey] = value

    return changeset


def get_species_keys(skey: str, species_ids: list[str]) -> list[str]:
    """Get keys of substance in given list of ids.

    Relies on naming patterns of ids. This does not get the species ids of the
    submodels, but only of the top model.

    Args:
        skey: substance key.
        species_ids: list of species ids to filter from.

    Returns:
        Species ids of the substance.
    """
    # The pattern is not very robust, e.g., blood vs. plasma.
    pattern = re.compile(rf"^[AC][a-z]+(_plasma)*_{skey}$")
    return [species_id for species_id in species_ids if pattern.search(species_id)]
