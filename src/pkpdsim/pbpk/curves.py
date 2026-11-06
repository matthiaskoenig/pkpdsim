"""Curves of the age: functions of the age as numbers and as SBML formulas.

The physiology over age (`pkpdsim.pbpk.age`) is a set of curves of the
covariate `AGE`, each of which is written once and evaluated twice: as a
number by `Curve.value` (the reference individuals of an age, the tables of
the documentation, the sampler of virtual individuals) and as the SBML
formula of an initial assignment by `Curve.formula` (the model). The tests
compare the two for every curve of the library, so the numbers of the
documentation are those of the model.

A curve is dimensionless: the value of a quantity of the physiology is its
reference value times a curve, or a curve of values in the unit of the
quantity whose formula carries the unit (`Curve.formula(..., unit=)`). The
age enters the formulas as the dimensionless number of years
(`AGE / 1 yr`), and the postmenstrual age of the maturation functions as
weeks after the conception, `PMA = 40 + AGE * 365.25 / 7` for a birth at
term (40 weeks).
"""

import bisect
import math
from abc import ABC, abstractmethod
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Literal

__all__ = [
    "DAYS_PER_YEAR",
    "PMA_AT_BIRTH",
    "WEEKS_PER_YEAR",
    "Affine",
    "Constant",
    "Curve",
    "Hill",
    "Interpolation",
    "Linear",
    "Logarithmic",
    "Product",
    "Quotient",
    "Relative",
    "Spliced",
    "check_curve",
    "pma",
]

PMA_AT_BIRTH = 40.0
"""The postmenstrual age at a birth at term (weeks)."""

DAYS_PER_YEAR = 365.25
"""The days of a year."""

WEEKS_PER_YEAR = DAYS_PER_YEAR / 7.0
"""The weeks of a year."""


def pma(age: float) -> float:
    """The postmenstrual age (weeks) of a postnatal age (years), born at term."""
    return PMA_AT_BIRTH + age * WEEKS_PER_YEAR


def _n(x: float) -> str:
    """A number of a formula."""
    return f"{float(x)!r} dimensionless"


def _at_least(x: str, low: str) -> str:
    """The larger of a formula and a bound, as a piecewise (SBML L3V1 has no `max`)."""
    return f"piecewise({low}, {x} < {low}, {x})"


def _variable(age: str, variable: Literal["age", "pma"]) -> str:
    """The variable of a curve in a formula, of the age in years `age`."""
    if variable == "age":
        return age
    return f"({_n(PMA_AT_BIRTH)} + {age} * {_n(WEEKS_PER_YEAR)})"


class Curve(ABC):
    """A dimensionless function of the age."""

    @abstractmethod
    def value(self, age: float) -> float:
        """The value at an age.

        Args:
            age: the postnatal age (years).
        """

    @abstractmethod
    def math(self, age: str) -> str:
        """The dimensionless formula of the curve.

        Args:
            age: the dimensionless formula of the age in years, `(AGE / 1 yr)`.
        """

    def formula(self, age: str, unit: str = "dimensionless") -> str:
        """The formula of the curve in a unit.

        Args:
            age: the id of the age (yr), `AGE`.
            unit: the id of the unit of the values of the curve.

        Returns:
            The formula, `(<math>) * 1 <unit>`; a dimensionless curve without
            the factor.
        """
        math_ = self.math(f"({age} / 1 yr)")
        if unit == "dimensionless":
            return math_
        return f"({math_}) * 1 {unit}"


@dataclass(frozen=True)
class Constant(Curve):
    """A constant.

    Attributes:
        c: the value.
    """

    c: float

    def value(self, age: float) -> float:
        """The constant."""
        return self.c

    def math(self, age: str) -> str:
        """The constant as a number."""
        return _n(self.c)


@dataclass(frozen=True)
class Interpolation(Curve):
    """The piecewise linear interpolation of values at ages, constant beyond them.

    Attributes:
        ages: the ages (years), increasing.
        values: the values at the ages.
    """

    ages: tuple[float, ...]
    values: tuple[float, ...]

    def __post_init__(self) -> None:
        """Check the points.

        Raises:
            ValueError: for no point, unequal lengths or ages which do not
                increase.
        """
        if not self.ages or len(self.ages) != len(self.values):
            raise ValueError(
                f"Interpolation: {len(self.ages)} ages and {len(self.values)} values."
            )
        if any(b <= a for a, b in zip(self.ages, self.ages[1:], strict=False)):
            raise ValueError(f"Interpolation: the ages {self.ages} do not increase.")

    def value(self, age: float) -> float:
        """The interpolated value."""
        if age <= self.ages[0]:
            return self.values[0]
        if age >= self.ages[-1]:
            return self.values[-1]
        i = bisect.bisect_right(self.ages, age) - 1
        a0, a1 = self.ages[i], self.ages[i + 1]
        v0, v1 = self.values[i], self.values[i + 1]
        return v0 + (v1 - v0) * (age - a0) / (a1 - a0)

    def math(self, age: str) -> str:
        """The piecewise formula."""
        if len(self.ages) == 1:
            return _n(self.values[0])
        terms = [_n(self.values[0]), f"{age} <= {_n(self.ages[0])}"]
        for i in range(len(self.ages) - 1):
            a0, a1 = self.ages[i], self.ages[i + 1]
            v0, v1 = self.values[i], self.values[i + 1]
            slope = (v1 - v0) / (a1 - a0)
            terms += [
                f"{_n(v0)} + {_n(slope)} * ({age} - {_n(a0)})",
                f"{age} < {_n(a1)}",
            ]
        terms.append(_n(self.values[-1]))
        return f"piecewise({', '.join(terms)})"


@dataclass(frozen=True)
class Linear(Curve):
    """A polynomial of the age, `c0 + c1*age + c2*age^2 + ...`.

    Attributes:
        coefficients: the coefficients from the constant on.
    """

    coefficients: tuple[float, ...]

    def value(self, age: float) -> float:
        """The polynomial at the age."""
        return math.fsum(c * age**k for k, c in enumerate(self.coefficients))

    def math(self, age: str) -> str:
        """The polynomial."""
        terms = []
        for k, c in enumerate(self.coefficients):
            if k == 0:
                terms.append(_n(c))
            elif k == 1:
                terms.append(f"{_n(c)} * {age}")
            else:
                terms.append(f"{_n(c)} * {age}^{k}")
        return f"({' + '.join(terms)})"


@dataclass(frozen=True)
class Hill(Curve):
    """A sigmoid maturation, `base + (top - base) * x^n / (x50^n + x^n)`.

    Attributes:
        x50: the variable at half of the maturation (years of age, or weeks
            of postmenstrual age).
        n: the Hill coefficient.
        base: the value at the variable 0.
        top: the value of the mature organ.
        variable: `age` (postnatal years) or `pma` (postmenstrual weeks, a
            birth at term).
        x_min: the smallest variable, below which the curve is constant (a
            function of the days which starts at the first day).
    """

    x50: float
    n: float
    base: float = 0.0
    top: float = 1.0
    variable: Literal["age", "pma"] = "age"
    x_min: float = 0.0

    def value(self, age: float) -> float:
        """The maturation at the age."""
        x = age if self.variable == "age" else pma(age)
        x = max(x, self.x_min)
        return self.base + (self.top - self.base) * x**self.n / (
            self.x50**self.n + x**self.n
        )

    def math(self, age: str) -> str:
        """The sigmoid."""
        x = _variable(age, self.variable)
        if self.variable == "age" or self.x_min > 0:
            x = _at_least(x, _n(self.x_min))
        p = f"{x}^{_n(self.n)}"
        return (
            f"({_n(self.base)} + {_n(self.top - self.base)} * {p} / "
            f"({_n(self.x50**self.n)} + {p}))"
        )


@dataclass(frozen=True)
class Logarithmic(Curve):
    """A logarithm of the postnatal age in days, `a + b*ln(max(days, 1))`.

    Attributes:
        a: the value at the first day.
        b: the slope per e-fold of the days.
    """

    a: float
    b: float

    def value(self, age: float) -> float:
        """The value at the age."""
        return self.a + self.b * math.log(max(age * DAYS_PER_YEAR, 1.0))

    def math(self, age: str) -> str:
        """The logarithm of the days, at least one."""
        days = _at_least(f"{age} * {_n(DAYS_PER_YEAR)}", "1 dimensionless")
        return f"({_n(self.a)} + {_n(self.b)} * ln({days}))"


@dataclass(frozen=True)
class Affine(Curve):
    """A constant plus multiples of curves, `c + sum(k_i * curve_i)`.

    Attributes:
        c: the constant.
        terms: the factors and their curves.
    """

    c: float
    terms: tuple[tuple[float, Curve], ...]

    def value(self, age: float) -> float:
        """The sum at the age."""
        return self.c + math.fsum(k * curve.value(age) for k, curve in self.terms)

    def math(self, age: str) -> str:
        """The sum."""
        terms = [_n(self.c), *(f"{_n(k)} * ({c.math(age)})" for k, c in self.terms)]
        return f"({' + '.join(terms)})"


@dataclass(frozen=True)
class Quotient(Curve):
    """The quotient of two curves.

    Attributes:
        numerator: the numerator.
        denominator: the denominator.
    """

    numerator: Curve
    denominator: Curve

    def value(self, age: float) -> float:
        """The quotient at the age."""
        return self.numerator.value(age) / self.denominator.value(age)

    def math(self, age: str) -> str:
        """The quotient."""
        return f"({self.numerator.math(age)}) / ({self.denominator.math(age)})"


@dataclass(frozen=True)
class Product(Curve):
    """The product of curves.

    Attributes:
        curves: the factors.
    """

    curves: tuple[Curve, ...]

    def value(self, age: float) -> float:
        """The product of the values."""
        return math.prod(c.value(age) for c in self.curves)

    def math(self, age: str) -> str:
        """The product of the formulas."""
        return " * ".join(f"({c.math(age)})" for c in self.curves)


@dataclass(frozen=True)
class Relative(Curve):
    """A curve relative to its value at an age, 1 at that age.

    Attributes:
        curve: the curve.
        at: the age (years) of the value 1.
    """

    curve: Curve
    at: float
    norm: float = field(init=False)

    def __post_init__(self) -> None:
        """The value of the curve at the age.

        Raises:
            ValueError: if it is not positive.
        """
        norm = self.curve.value(self.at)
        if not norm > 0:
            raise ValueError(f"Relative: the curve is {norm} at the age {self.at}.")
        object.__setattr__(self, "norm", norm)

    def value(self, age: float) -> float:
        """The value relative to the value at `at`."""
        return self.curve.value(age) / self.norm

    def math(self, age: str) -> str:
        """The formula over the formula at `at`, so the curve is exactly 1 at `at`."""
        return f"({self.curve.math(age)}) / ({self.curve.math(_n(self.at))})"


@dataclass(frozen=True)
class Spliced(Curve):
    """Curves of consecutive age ranges, continuous at the splices.

    The first curve holds below the first age of `ages`, the second from it
    to the second, and so on; the values are not adjusted, so the curves must
    meet at the splices (checked).

    Attributes:
        ages: the ages of the splices (years), increasing.
        curves: one more curve than splices.
    """

    ages: tuple[float, ...]
    curves: tuple[Curve, ...]

    def __post_init__(self) -> None:
        """Check the splices.

        Raises:
            ValueError: for a wrong number of curves or curves which do not
                meet at a splice (relative difference above 1e-9).
        """
        if len(self.curves) != len(self.ages) + 1:
            raise ValueError(
                f"Spliced: {len(self.curves)} curves for {len(self.ages)} splices."
            )
        for a, left, right in zip(
            self.ages, self.curves, self.curves[1:], strict=False
        ):
            lv, rv = left.value(a), right.value(a)
            if not math.isclose(lv, rv, rel_tol=1e-9, abs_tol=1e-12):
                raise ValueError(
                    f"Spliced: the curves do not meet at the age {a}: {lv} and {rv}."
                )

    def _index(self, age: float) -> int:
        """The curve of an age."""
        return bisect.bisect_right(self.ages, age)

    def value(self, age: float) -> float:
        """The value of the curve of the age."""
        return self.curves[self._index(age)].value(age)

    def math(self, age: str) -> str:
        """The piecewise formula of the curves."""
        terms: list[str] = []
        for a, c in zip(self.ages, self.curves, strict=False):
            terms += [c.math(age), f"{age} < {_n(a)}"]
        terms.append(self.curves[-1].math(age))
        return f"piecewise({', '.join(terms)})"


def check_curve(curve: Curve, ages: Sequence[float]) -> list[str]:
    """Check that a curve is finite and positive at ages.

    Args:
        curve: the curve.
        ages: the ages (years).

    Returns:
        The problems.
    """
    problems = []
    for a in ages:
        v = curve.value(a)
        if not (math.isfinite(v) and v > 0):
            problems.append(f"the curve is {v} at the age {a}")
    return problems
