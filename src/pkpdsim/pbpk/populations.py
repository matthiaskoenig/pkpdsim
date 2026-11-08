"""The populations: the reference individual and the covariates over age.

A population is a set of data of the library (`database/tables/populations/<key>.toml`):
the reference covariates of each sex (body weight, height and, where the data
give it, the age), which become the default individual of a model, and the
distribution of the body weight and the height of each sex over age
(`Anthropometry`), from which `pkpdsim.pbpk.virtual` draws virtual
individuals. The reference anatomy stays ICRP 89 (`database/tables/physiology_human.toml`,
with its children and the changes of the adults over age of
`pkpdsim.pbpk.age`): the organ volumes and flows of an individual follow
from the covariate scaling of `pkpdsim.pbpk.physiology` (Huisinga 2012) at
its body weight, height and age. The population is chosen when the model is
built (`PBPKDefinition.population`), the European reference individuals of
ICRP 89 are the default.

The distributions are data of the library and never SBML `distrib`
uncertainties: `run_simulation --uncertainty` draws every parameter with a
`distrib` uncertainty and must only draw the uncertainty of a fit.

The populations:

- `european_icrp`: the ICRP 89 reference individuals (the default); ICRP 89
  gives reference values, not distributions, so the population has none;
- `white_american_nhanes`, `black_american_nhanes`, `mexican_american_nhanes`:
  the NHANES examinations 2013-2018 (`scripts/database/nhanes_populations.py`
  from the public files of the survey): the reference individual the
  survey-weighted mean of the adults of 20 to 49 years, the distributions of
  every age from birth to 80 years and over, survey-weighted kernel moments;
- `east_asian_tanaka`: the Asian Reference Man of Tanaka 1996 with the
  distributions of the Japanese of the same tables.

PK-Sim also offers a Japanese population (based on Japanese national
statistics of 2015). Its primary data were not read, so the library has no
`japanese` population; the Japanese data of Tanaka 1996 (1976 to 1988) are
the distributions of `east_asian_tanaka`.
"""

from collections.abc import Mapping
from dataclasses import dataclass
from functools import cache
from typing import Any

import numpy as np
from sbmlutils.factory import UnitDefinition

from pkpdsim.database.tables import load_folder
from pkpdsim.pbpk import ids
from pkpdsim.pbpk.data import load
from pkpdsim.pbpk.physiology import AGE_REFERENCE, DEFAULT_POPULATION, Sex
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Value

__all__ = [
    "COVARIATES",
    "DEFAULT_POPULATION",
    "POPULATIONS",
    "Anthropometry",
    "CovariateDistribution",
    "Population",
    "population",
]

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
class CovariateDistribution:
    """The distribution of height and body weight of the individuals of one age and sex.

    The height (cm) and the logarithm of the body weight (kg) are bivariate
    normal: the body weight is lognormal and correlated with the height.

    Attributes:
        height_mean: the mean height (cm).
        height_sd: the standard deviation of the height (cm).
        ln_bw_mean: the mean of the logarithm of the body weight (ln kg).
        ln_bw_sd: the standard deviation of the logarithm of the body weight.
        rho: the correlation of height and logarithm of the body weight;
            `None` where the data give none, the two are drawn independently.
    """

    height_mean: float
    height_sd: float
    ln_bw_mean: float
    ln_bw_sd: float
    rho: float | None

    @property
    def bw_mean(self) -> float:
        """The mean body weight (kg) of the lognormal distribution."""
        return float(np.exp(self.ln_bw_mean + self.ln_bw_sd**2 / 2))


@dataclass(frozen=True)
class Anthropometry:
    """The distribution of height and body weight of one sex of a population over age.

    The rows are the distributions at ages (`CovariateDistribution`); between
    two ages every moment is interpolated linearly, below the first and above
    the last age of the rows the nearest row holds, within `age_range`, the
    ages the data describe.

    Attributes:
        ages: the ages of the rows (years), increasing.
        height_mean: the mean height of each row (cm).
        height_sd: the standard deviation of the height of each row (cm).
        ln_bw_mean: the mean logarithm of the body weight of each row (ln kg).
        ln_bw_sd: the standard deviation of the logarithm of the body weight.
        rho: the correlation of height and log body weight of each row, or
            `None`: the data give no correlation, the two are independent.
        age_range: the ages (years) the distributions describe.
        source: the key of the source in `database/sources.toml`.
        notes: where the numbers come from.
    """

    ages: tuple[float, ...]
    height_mean: tuple[float, ...]
    height_sd: tuple[float, ...]
    ln_bw_mean: tuple[float, ...]
    ln_bw_sd: tuple[float, ...]
    rho: tuple[float, ...] | None
    age_range: tuple[float, float]
    source: str
    notes: str | None = None

    def __post_init__(self) -> None:
        """Check the rows.

        Raises:
            ValueError: for rows of different lengths, ages which do not
                increase or lie outside the age range, a standard deviation
                which is not positive or a correlation outside of (-1, 1).
        """
        columns = [self.height_mean, self.height_sd, self.ln_bw_mean, self.ln_bw_sd]
        if self.rho is not None:
            columns.append(self.rho)
        if not self.ages or any(len(c) != len(self.ages) for c in columns):
            raise ValueError("Anthropometry: the columns differ in length.")
        if any(b <= a for a, b in zip(self.ages, self.ages[1:], strict=False)):
            raise ValueError(f"Anthropometry: the ages {self.ages} do not increase.")
        lo, hi = self.age_range
        if not lo <= self.ages[0] or not self.ages[-1] <= hi:
            raise ValueError(
                f"Anthropometry: the ages {self.ages} lie outside of the age "
                f"range {self.age_range}."
            )
        if min(self.height_sd) <= 0 or min(self.ln_bw_sd) <= 0:
            raise ValueError("Anthropometry: a standard deviation is not positive.")
        if self.rho is not None and any(not -1 < r < 1 for r in self.rho):
            raise ValueError("Anthropometry: a correlation is outside of (-1, 1).")

    def at(self, age: float) -> CovariateDistribution:
        """The distribution of an age.

        Args:
            age: the age (years).

        Returns:
            The moments interpolated linearly in age.

        Raises:
            ValueError: for an age outside of `age_range`.
        """
        lo, hi = self.age_range
        if not lo <= age <= hi:
            raise ValueError(
                f"The age {age} is outside of the ages {lo} to {hi} of the data."
            )

        def interp(column: tuple[float, ...]) -> float:
            return float(np.interp(age, self.ages, column))

        return CovariateDistribution(
            height_mean=interp(self.height_mean),
            height_sd=interp(self.height_sd),
            ln_bw_mean=interp(self.ln_bw_mean),
            ln_bw_sd=interp(self.ln_bw_sd),
            rho=None if self.rho is None else interp(self.rho),
        )

    def sample(
        self, ages: np.ndarray, rng: np.random.Generator, truncation: float = 3.0
    ) -> tuple[np.ndarray, np.ndarray]:
        """Draw the height and body weight of individuals of given ages.

        Every individual is drawn from the bivariate normal distribution of
        height and log body weight of its age, truncated at `truncation`
        standard deviations of each (the draws beyond are drawn again).

        Args:
            ages: the ages of the individuals (years).
            rng: the random generator.
            truncation: the truncation of both marginals in standard
                deviations.

        Returns:
            The heights (cm) and the body weights (kg).

        Raises:
            ValueError: for an age outside of `age_range` or a truncation
                which is not positive.
        """
        if not truncation > 0:
            raise ValueError(f"The truncation {truncation} is not positive.")
        rows = [self.at(float(a)) for a in ages]
        n = len(rows)
        z = np.empty((n, 2))
        todo = np.arange(n)
        while todo.size:
            draw = rng.standard_normal((todo.size, 2))
            ok = np.all(np.abs(draw) <= truncation, axis=1)
            z[todo[ok]] = draw[ok]
            todo = todo[~ok]
        heights = np.empty(n)
        bws = np.empty(n)
        for i, row in enumerate(rows):
            rho = 0.0 if row.rho is None else row.rho
            z1 = z[i, 0]
            z2 = rho * z1 + np.sqrt(1.0 - rho**2) * z[i, 1]
            heights[i] = row.height_mean + row.height_sd * z1
            bws[i] = np.exp(row.ln_bw_mean + row.ln_bw_sd * z2)
        return heights, bws


@dataclass(frozen=True)
class Population:
    """A population: the reference individual and the covariates over age per sex.

    Attributes:
        key: the key of the population, `PBPKDefinition.population`.
        name: the name.
        source: the key of the source of the population in `database/sources.toml`.
        reference: the reference covariates (`BW`, `HEIGHT`, `AGE`) of each
            sex, the default individual of a model.
        anthropometry: the distribution of height and body weight over age
            of each sex; a sex without data has none.
    """

    key: str
    name: str
    source: str
    reference: Mapping[Sex, Mapping[str, Value]]
    anthropometry: Mapping[Sex, Anthropometry]


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
            AGE_REFERENCE,
            U.yr,
            "ICRP89",
            "Section 1.4.2 para 14: adults are 20-50 years, the adult reference "
            "values are tabulated at the average age of 35 years",
        ),
    }


def _floats(table: Mapping[str, Any], key: str) -> tuple[float, ...]:
    """A column of the rows of an anthropometry."""
    return tuple(float(v) for v in table[key])


def _anthropometry(key: str, table: Mapping[str, Any]) -> Anthropometry:
    """The anthropometry of a sex of a population file.

    The log body weight is given by its moments (`ln_bw_mean`, `ln_bw_sd`) or
    by the mean and standard deviation of the lognormal body weight
    (`bw_mean`, `bw_sd`).

    Raises:
        ValueError: for a table with neither form of the body weight.
    """
    if "ln_bw_mean" in table:
        ln_mean, ln_sd = _floats(table, "ln_bw_mean"), _floats(table, "ln_bw_sd")
    elif "bw_mean" in table:
        mean, sd = (
            np.array(_floats(table, "bw_mean")),
            np.array(_floats(table, "bw_sd")),
        )
        var = np.log1p((sd / mean) ** 2)
        ln_mean = tuple(float(v) for v in np.log(mean) - var / 2)
        ln_sd = tuple(float(v) for v in np.sqrt(var))
    else:
        raise ValueError(
            f"population '{key}': `by_age` gives neither `ln_bw_mean` nor `bw_mean`."
        )
    lo, hi = table["age_range"]
    return Anthropometry(
        ages=_floats(table, "age"),
        height_mean=_floats(table, "height_mean"),
        height_sd=_floats(table, "height_sd"),
        ln_bw_mean=ln_mean,
        ln_bw_sd=ln_sd,
        rho=_floats(table, "rho") if "rho" in table else None,
        age_range=(float(lo), float(hi)),
        source=table["source"],
        notes=table.get("notes"),
    )


def _population(key: str, data: Mapping[str, Any]) -> Population:
    """A population of its data file."""
    reference: dict[Sex, Mapping[str, Value]] = {}
    anthropometry: dict[Sex, Anthropometry] = {}
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
        if "by_age" in table:
            anthropometry[sex] = _anthropometry(key, table["by_age"])
    return Population(
        key=key,
        name=data["name"],
        source=data["source"],
        reference=reference,
        anthropometry=anthropometry,
    )


@cache
def _populations() -> dict[str, Population]:
    """All populations of the tables `populations` of the database, by key."""
    return {key: _population(key, data) for key, data in load_folder("populations")}


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
