"""The populations: the reference individual and the distributions of the covariates.

A population is a set of data of the library (`data/populations/<key>.toml`):
the reference covariates of each sex (body weight, height and, where the data
give it, the age), which become the default individual of a model, and the
distributions of the covariates for a later sampler of virtual individuals.
The reference anatomy stays ICRP 89 (`data/physiology_human.toml`): the organ
volumes and flows of a population's reference individual follow from the
covariate scaling of `pkpdsim.pbpk.physiology` (Huisinga 2012) at its body
weight and height. The population is chosen when the model is built
(`PBPKDefinition.population`), the European reference individuals of ICRP 89
are the default.

The distributions are data of the library and never SBML `distrib`
uncertainties: `run_simulation --uncertainty` draws every parameter with a
`distrib` uncertainty and must only draw the uncertainty of a fit.

The populations:

- `european_icrp`: the ICRP 89 reference individuals (the default); ICRP 89
  gives reference values, not distributions, so the population has none;
- `white_american_nhanes`, `black_american_nhanes`, `mexican_american_nhanes`:
  adults of 20 to 49 years of the NHANES examinations 2013-2018 as
  preprocessed by httk (HTTK-Pop, Ring 2017), survey-weighted means and
  standard deviations;
- `east_asian_tanaka`: the Asian Reference Man of Tanaka 1996 (Tables 9 and
  10, 20 to 49 years) with the distributions of the Japanese of the same
  tables.

PK-Sim also offers a Japanese population (based on Japanese national
statistics of 2015). Its primary data were not read, so the library has no
`japanese` population; the Japanese data of Tanaka 1996 (1976 to 1988) are
the distributions of `east_asian_tanaka`.
"""

import tomllib
from collections.abc import Mapping
from dataclasses import dataclass
from functools import cache
from importlib import resources
from typing import Any, Literal

from sbmlutils.factory import UnitDefinition

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.data import load
from pkpdsim.pbpk.physiology import AGE_DEFAULT, MODEL_DEFAULT, Sex
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Value

__all__ = [
    "COVARIATES",
    "DEFAULT_POPULATION",
    "POPULATIONS",
    "Distribution",
    "Population",
    "population",
]

DEFAULT_POPULATION = "european_icrp"
"""The population of a model which does not choose one."""

COVARIATES: tuple[str, ...] = (ids.BW, ids.HEIGHT, ids.AGE)
"""The covariates a population gives reference values for."""

_UNITS: dict[str, UnitDefinition] = {"kg": U.kg, "cm": U.cm, "yr": U.yr}
"""The units of the model for the units of the population files."""

_COVARIATE_UNITS: dict[str, UnitDefinition] = {
    ids.BW: U.kg,
    ids.HEIGHT: U.cm,
    ids.AGE: U.yr,
}


@dataclass(frozen=True)
class Distribution:
    """The distribution of a covariate in a population.

    Attributes:
        kind: `normal` or `lognormal`; the mean and the standard deviation are
            those of the covariate itself (not of its logarithm) for both.
        mean: the mean.
        sd: the standard deviation.
        unit: the unit of mean and standard deviation.
        source: the key of the source in `data/sources.toml`.
        notes: where the numbers come from.
    """

    kind: Literal["normal", "lognormal"]
    mean: float
    sd: float
    unit: UnitDefinition
    source: str
    notes: str | None = None


@dataclass(frozen=True)
class Population:
    """A population: the reference individual and the covariate distributions per sex.

    Attributes:
        key: the key of the population, `PBPKDefinition.population`.
        name: the name.
        source: the key of the source of the population in `data/sources.toml`.
        reference: the reference covariates (`BW`, `HEIGHT`, `AGE`) of each
            sex, the default individual of a model.
        covariates: the distributions of the covariates of each sex.
    """

    key: str
    name: str
    source: str
    reference: Mapping[Sex, Mapping[str, Value]]
    covariates: Mapping[Sex, Mapping[str, Distribution]]


def _unit(key: str, entry: Mapping[str, Any], covariate: str) -> UnitDefinition:
    """The unit of an entry, which must be the unit of the covariate in the model.

    Raises:
        ValueError: if the covariate or its unit is not one of the model.
    """
    if covariate not in _COVARIATE_UNITS:
        raise ValueError(
            f"population '{key}': unknown covariate '{covariate}', the "
            f"covariates are {list(COVARIATES)}."
        )
    unit = _UNITS.get(entry["unit"])
    if unit is not _COVARIATE_UNITS[covariate]:
        raise ValueError(
            f"population '{key}': {covariate} is given in '{entry['unit']}', the "
            f"model needs '{_COVARIATE_UNITS[covariate].sid}'."
        )
    return unit


def _icrp_reference(sex: Sex) -> dict[str, Value]:
    """The reference covariates of the ICRP 89 reference individual of a sex."""
    p = load("physiology_human")[sex.value]
    return {
        ids.BW: Value(
            float(p["bodyweight"]["value"]),
            U.kg,
            p["bodyweight"]["source"],
            p["bodyweight"].get("notes"),
        ),
        ids.HEIGHT: Value(
            float(p["height"]["value"]),
            U.cm,
            p["height"]["source"],
            p["height"].get("notes"),
        ),
        ids.AGE: Value(
            AGE_DEFAULT, U.yr, None, f"{MODEL_DEFAULT} No rule uses AGE yet."
        ),
    }


def _population(key: str, data: Mapping[str, Any]) -> Population:
    """A population of its data file."""
    reference: dict[Sex, Mapping[str, Value]] = {}
    covariates: dict[Sex, Mapping[str, Distribution]] = {}
    for sex in Sex:
        table = data.get(sex.value, {})
        values = _icrp_reference(sex)
        for covariate, entry in table.get("reference", {}).items():
            values[covariate] = Value(
                float(entry["value"]),
                _unit(key, entry, covariate),
                entry["source"],
                entry.get("notes"),
            )
        reference[sex] = values
        distributions: dict[str, Distribution] = {}
        for covariate, entry in table.get("covariates", {}).items():
            if entry["kind"] not in ("normal", "lognormal"):
                raise ValueError(
                    f"population '{key}': {covariate} has the unknown kind "
                    f"'{entry['kind']}', not normal or lognormal."
                )
            distributions[covariate] = Distribution(
                kind=entry["kind"],
                mean=float(entry["mean"]),
                sd=float(entry["sd"]),
                unit=_unit(key, entry, covariate),
                source=entry["source"],
                notes=entry.get("notes"),
            )
        covariates[sex] = distributions
    return Population(
        key=key,
        name=data["name"],
        source=data["source"],
        reference=reference,
        covariates=covariates,
    )


@cache
def _populations() -> dict[str, Population]:
    """All populations of `data/populations`, by key."""
    folder = resources.files("pkpdsim.pbpk").joinpath("data", "populations")
    found: dict[str, Population] = {}
    for path in sorted(folder.iterdir(), key=lambda p: p.name):
        if path.name.endswith(".toml"):
            key = path.name.removesuffix(".toml")
            found[key] = _population(key, tomllib.loads(path.read_text("utf-8")))
    return found


POPULATIONS: dict[str, Population] = _populations()
"""The populations by key."""


def population(key: str) -> Population:
    """Get a population.

    Args:
        key: the key of the population.

    Returns:
        The population.

    Raises:
        ValueError: if there is no population of the key.
    """
    try:
        return POPULATIONS[key]
    except KeyError:
        raise ValueError(
            f"Unknown population '{key}', the populations are {sorted(POPULATIONS)}."
        ) from None
