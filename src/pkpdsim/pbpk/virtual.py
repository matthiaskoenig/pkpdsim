"""Virtual populations: individuals drawn from the factors of a population.

A virtual individual is a combination of factors, each of which the library
knows: the population (the distribution of body weight and height over age of
`pkpdsim.pbpk.populations`), the sex (the build of the model), the age (the
physiology over age of `pkpdsim.pbpk.age`) and the impaired organ function
(the severities of `pkpdsim.pbpk.impairment`). `virtual_population` draws
them: the sex by its share, the age uniformly within a range or an age class
(`pkpdsim.impairment.AGE_CLASSES`), the height and the body weight from the
bivariate distribution of the age and sex of the population, and a severity
as a fixed value or uniformly within a range. The changes of an individual
(`AGE`, `BW`, `HEIGHT` and the severities) are initial values of the model:
set before the initialization, the physiology of the build of the sex follows
(`VirtualPopulation.values`, `PBPKSimulationExperiment.population_dimension`).

The draws are reproducible by the seed; the covariates are library data and
never SBML `distrib` uncertainties.
"""

from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType

import numpy as np

from pkpdsim.impairment import AGE_CLASSES
from pkpdsim.pbpk import ids
from pkpdsim.pbpk.populations import population as get_population
from pkpdsim.sex import Sex

__all__ = [
    "MAX_AGE",
    "SEVERITIES",
    "VirtualIndividual",
    "VirtualPopulation",
    "virtual_population",
]

MAX_AGE = 100.0
"""The oldest age (years) of an age class without an upper limit."""

SEVERITIES: Mapping[str, tuple[float, float]] = MappingProxyType(
    {
        ids.CP_SEVERITY: (0.0, 3.0),
        ids.HF_SEVERITY: (0.0, 4.0),
        ids.F_RENAL_FUNCTION: (0.0, 1.5),
    }
)
"""The conditions a virtual individual may have, with the range of their value."""


@dataclass(frozen=True)
class VirtualIndividual:
    """An individual of a virtual population.

    Attributes:
        index: the number of the individual in its population.
        sex: the sex, whose build simulates it.
        age: the age (years).
        bw: the body weight (kg).
        height: the height (cm).
        conditions: the severities of the impairments by parameter id
            (`CP_severity`, `f_renal_function`, `HF_severity`).
    """

    index: int
    sex: Sex
    age: float
    bw: float
    height: float
    conditions: Mapping[str, float] = field(default_factory=dict)

    @property
    def bmi(self) -> float:
        """The body mass index (kg/m2)."""
        return self.bw / (self.height / 100.0) ** 2

    def changes(self) -> dict[str, float]:
        """The initial values of the individual by parameter id, in the units of the model."""
        return {
            ids.AGE: self.age,
            ids.BW: self.bw,
            ids.HEIGHT: self.height,
            **self.conditions,
        }


@dataclass(frozen=True)
class VirtualPopulation:
    """A virtual population: the individuals and how they were drawn.

    Attributes:
        population: the key of the population of the covariates.
        individuals: the individuals.
        seed: the seed of the draws.
    """

    population: str
    individuals: tuple[VirtualIndividual, ...]
    seed: int

    def __len__(self) -> int:
        """The number of individuals."""
        return len(self.individuals)

    def of_sex(self, sex: Sex) -> tuple[VirtualIndividual, ...]:
        """The individuals of a sex, simulated with the build of the sex."""
        return tuple(i for i in self.individuals if i.sex is Sex(sex))

    def values(self, sex: Sex) -> dict[str, np.ndarray]:
        """The initial values of the individuals of a sex, an array per parameter id.

        The arrays are the values of a dimension of a scan (`AGE` in years,
        `BW` in kg, `HEIGHT` in cm, the severities dimensionless); a
        condition which only some individuals have is the value of the
        healthy organ function for the others (`f_renal_function` 1, the
        severities 0).

        Raises:
            ValueError: if the population has no individual of the sex.
        """
        members = self.of_sex(sex)
        if not members:
            raise ValueError(
                f"The virtual population has no individual of the sex {Sex(sex).value}."
            )
        keys = sorted({k for i in members for k in i.conditions})
        healthy = {ids.F_RENAL_FUNCTION: 1.0}
        values = {
            ids.AGE: np.array([i.age for i in members]),
            ids.BW: np.array([i.bw for i in members]),
            ids.HEIGHT: np.array([i.height for i in members]),
        }
        for k in keys:
            values[k] = np.array(
                [i.conditions.get(k, healthy.get(k, 0.0)) for i in members]
            )
        return values


def _ages(
    ages: tuple[float, float] | str, n: int, rng: np.random.Generator
) -> np.ndarray:
    """Draw ages uniformly within a range or an age class.

    Raises:
        ValueError: for an empty range or an unknown class.
    """
    if isinstance(ages, str):
        limits = AGE_CLASSES[ages].limits
        if limits is None or limits[0] is None:
            raise ValueError(f"The age class '{ages}' has no lower limit.")
        lo = limits[0]
        hi = MAX_AGE if limits[1] is None else limits[1]
    else:
        lo, hi = ages
    if not 0.0 <= lo < hi:
        raise ValueError(f"The ages {lo} to {hi} are no range of ages.")
    return rng.uniform(lo, hi, n)


def _conditions(
    conditions: Mapping[str, float | tuple[float, float]],
    n: int,
    rng: np.random.Generator,
) -> dict[str, np.ndarray]:
    """Draw the severities: a fixed value or uniformly within a range.

    Raises:
        ValueError: for an unknown condition or a value outside of its range.
    """
    out: dict[str, np.ndarray] = {}
    for key, spec in conditions.items():
        if key not in SEVERITIES:
            raise ValueError(
                f"Unknown condition '{key}', the conditions are {sorted(SEVERITIES)}."
            )
        lo, hi = SEVERITIES[key]
        a, b = (spec, spec) if isinstance(spec, int | float) else spec
        if not lo <= a <= b <= hi:
            raise ValueError(f"'{key}' {spec} is outside of {lo} to {hi}.")
        out[key] = np.full(n, float(a)) if a == b else rng.uniform(a, b, n)
    return out


def virtual_population(
    population: str,
    n: int,
    *,
    ages: tuple[float, float] | str = (20.0, 50.0),
    female: float = 0.5,
    conditions: Mapping[str, float | tuple[float, float]] | None = None,
    seed: int = 1234,
    truncation: float = 3.0,
) -> VirtualPopulation:
    """Draw a virtual population.

    Args:
        population: the key of the population of the covariates
            (`pkpdsim.pbpk.populations`), which must have distributions of
            both sexes drawn (`european_icrp` has none).
        n: the number of individuals.
        ages: the range of the ages (years), drawn uniformly, or the key of an
            age class of `pkpdsim.impairment.AGE_CLASSES` (`child`,
            `old_old`, ...); the population must describe the ages.
        female: the share of women, 0 (men only) to 1 (women only).
        conditions: the severities of the impairments by parameter id, a
            fixed value or a range drawn uniformly (`{"CP_severity": (1, 3)}`).
        seed: the seed of the random generator.
        truncation: the truncation of height and log body weight in standard
            deviations of the age.

    Returns:
        The population.

    Raises:
        ValueError: for no individual, a share outside of 0 to 1, a
            population without the distributions of a drawn sex, ages it does
            not describe, or an unknown condition.
    """
    if n < 1:
        raise ValueError(f"A virtual population needs individuals, not {n}.")
    if not 0.0 <= female <= 1.0:
        raise ValueError(f"The share of women {female} is outside of 0 to 1.")
    rng = np.random.default_rng(seed)
    p = get_population(population)
    sexes = np.where(rng.uniform(size=n) < female, 1, 0)
    age = _ages(ages, n, rng)
    drawn = _conditions(conditions or {}, n, rng)
    height = np.empty(n)
    bw = np.empty(n)
    for code, sex in ((0, Sex.MALE), (1, Sex.FEMALE)):
        members = np.flatnonzero(sexes == code)
        if not members.size:
            continue
        anthropometry = p.anthropometry.get(sex)
        if anthropometry is None:
            raise ValueError(
                f"The population '{population}' has no distribution of the "
                f"covariates of the sex {sex.value}."
            )
        height[members], bw[members] = anthropometry.sample(
            age[members], rng, truncation
        )
    individuals = tuple(
        VirtualIndividual(
            index=i,
            sex=Sex.FEMALE if sexes[i] else Sex.MALE,
            age=float(age[i]),
            bw=float(bw[i]),
            height=float(height[i]),
            conditions=MappingProxyType({k: float(v[i]) for k, v in drawn.items()}),
        )
        for i in range(n)
    )
    return VirtualPopulation(population=population, individuals=individuals, seed=seed)
