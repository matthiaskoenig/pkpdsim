"""The adipose tissue, the body composition and the physiology of obesity.

The adipose tissue of the library is the remainder of the body weight after
the lean organs and the blood (`pkpdsim.pbpk.physiology`: the fat-free mass of
Janmahasatian 2005 and the scaling of Huisinga 2012), so the masses of an
individual always sum to its body weight and the body composition needs no
equation of its own. `PBPKDefinition.adipose` (`Adipose`) adds what the
remainder does not give; it is a choice of the build, and a build without it
is unchanged:

- **The body composition** of every build with the option: the lean body
  weight `LBW = BW - Mad`, the body fat mass `FM = BW - FFM` with the body fat
  fraction `BFF` and the fat mass index `FMI`, and the split of the adipose
  tissue into the visceral `VAT = f_vat*Mad` and the subcutaneous
  `SAT = Mad - VAT` with the visceral share `f_vat` of Shen 2009 (a function
  of the age of the sex). The classes of the body mass index are
  `pkpdsim.impairment.BMI_CLASSES`, which an experiment scans and draws as the
  classes of an impairment.
- **The physiology of overweight and obesity**: the factors `fobes_<key>` of
  `database/tables/adipose_human.toml` on the cardiac output, on the flow
  fractions of the adipose tissue, the muscle, the brain, the skin, the heart,
  the gut and the pancreas, on the plasma albumin and alpha-1-acid
  glycoprotein and on the hepatic CYP3A4 and CYP2E1, each the residual of the
  repository of Berton 2022 over what the scaling of the library already does
  and each 1 up to a body mass index of 25 kg/m2, so a normal weight
  individual is unchanged (`Adipose(obesity=False)` leaves them out).
- **The change of the body weight** over time for a model which gives an
  effect (`Adipose(effect=Formula(...))`, the GLP-1 agonists): the effect is
  the fractional change of the body weight at the set point (`f_bw`), the
  weight approaches `BW_set = BW0*(1 + f_bw)` with the rate `k_bw`, and `BW`
  is a state of the model (a rate rule). The physiology follows the weight:
  `BW` is a dynamic parameter, so every organ mass, volume and flow, the blood
  volume, the cardiac output, the GFR and the factors of obesity become
  assignment rules of it, as they are of a covariate which an experiment
  changes before a simulation. The fat and the lean mass of the change are
  those of the fat-free mass of Janmahasatian 2005: 0.69 to 0.73 of a change
  of the weight of an obese individual is adipose tissue, which is the fat
  fraction the DXA substudies of the GLP-1 trials measure (0.74 in
  SURMOUNT-1, 0.709 with an SD of 0.19 over 40 reports) and which the Forbes
  relation predicts, see `docs/pbpk/adipose.md`.

The body weight of the covariates is the weight at the start, `BW0`, and the
model reports `BW_change` and `BW_ratio` over it.
"""

import math
from collections.abc import Mapping
from dataclasses import dataclass
from functools import cache
from types import MappingProxyType
from typing import TYPE_CHECKING, Any

from sbmlutils.factory import (
    AssignmentRule,
    InitialAssignment,
    Parameter,
    RateRule,
    Sbase,
)
from sbmlutils.factory.units import UnitType
from sbmlutils.metadata import SBO

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.age import AGE_REFERENCE, AgeCurve
from pkpdsim.pbpk.annotation import annotate, initial_assignment
from pkpdsim.pbpk.curves import Curve, Interpolation
from pkpdsim.pbpk.data import cite, load, value
from pkpdsim.pbpk.organs import Organ
from pkpdsim.pbpk.physiology import body_composition
from pkpdsim.pbpk.units import KG_PER_MIN, U
from pkpdsim.pbpk.values import Formula, Value, number
from pkpdsim.sex import Sex
from pkpdsim.species import Species as BuildSpecies

if TYPE_CHECKING:
    from pkpdsim.pbpk.definition import PBPKDefinition

__all__ = [
    "ADIPOSE",
    "BINDING_FACTORS",
    "BMI_NORMAL",
    "CARDIAC_OUTPUT",
    "FLOW_FACTORS",
    "Adipose",
    "ObesityFactor",
    "StartValue",
    "VisceralComposition",
    "WaistModel",
    "adipose_dynamic",
    "adipose_objects",
    "obesity_factor",
    "obesity_factors",
    "obesity_protein_factors",
    "visceral_composition",
    "visceral_share",
    "waist_factor",
    "waist_model",
]

ADIPOSE = "adipose_human"
"""The table of the adipose tissue and the body composition in the database."""

BMI_NORMAL = 25.0
"""The body mass index (kg/m2) up to which the factors of obesity are 1: the
threshold of overweight of the WHO classification, at which the factors of the
table are normalised."""

CARDIAC_OUTPUT = "cardiac_output"
"""The key of the factor of obesity on the cardiac output, which carries every
organ flow, so the cardiac output of a model stays the sum of its flows."""

FLOW_FACTORS: Mapping[str, tuple[Organ, ...]] = MappingProxyType(
    {
        "adipose_flow": (Organ.ADIPOSE,),
        "muscle_flow": (Organ.MUSCLE,),
        "brain_flow": (Organ.BRAIN,),
        "skin_flow": (Organ.SKIN,),
        "heart_flow": (Organ.HEART,),
        "gut_flow": (Organ.STOMACH, Organ.SMALL_INTESTINE, Organ.LARGE_INTESTINE),
        "pancreas_flow": (Organ.PANCREAS,),
    }
)
"""The organs whose flow fraction a factor of obesity multiplies, by the key of
the factor in the table; the kidney, the spleen, the bone, the gonads and the
hepatic artery keep their fraction (`docs/pbpk/adipose.md`)."""

BINDING_FACTORS = ("albumin", "agp")
"""The plasma proteins with a factor of obesity on their concentration."""

_SBO_PARAMETER = SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER


class StartValue(InitialAssignment):
    """The initial assignment of a constant which records a value at the start.

    It stays an initial assignment when it reads a dynamic quantity
    (`pkpdsim.pbpk.builder`), as the start of a species of an endogenous
    substance does: `BW0` is the body weight at the start of the simulation,
    of which the state `BW` then departs.
    """


def obesity_factor(pid: str) -> str:
    """The id of the factor of obesity on a quantity, `fobes_CO`, `fobes_Qad`."""
    return f"fobes_{pid}"


@dataclass(frozen=True)
class ObesityFactor:
    """A factor of overweight and obesity over the body mass index.

    Attributes:
        curve: the dimensionless factor over the body mass index (kg/m2), 1 up
            to `BMI_NORMAL`.
        sources: the keys of the sources.
        notes: where the values come from.
    """

    curve: Curve
    sources: tuple[str, ...]
    notes: str

    def value(self, bmi: float) -> float:
        """The factor at a body mass index (kg/m2)."""
        return self.curve.value(bmi)

    def formula(self) -> str:
        """The formula of the factor over `BMI`."""
        return self.curve.math(f"({ids.BMI} / 1 kg_per_m2)")


@cache
def _table() -> dict[str, Any]:
    """The table of the adipose tissue and the body composition."""
    return load(ADIPOSE)


def _provenance(entry: Mapping[str, Any]) -> tuple[tuple[str, ...], str]:
    """The keys of the sources and the notes of an entry of the table."""
    sources = tuple(dict.fromkeys((entry["source"], *entry.get("sources", ()))))
    return sources, str(entry["notes"])


def _interpolation(entry: Mapping[str, Any], x: str, sex: Sex) -> Interpolation:
    """The interpolation of an entry over `x`, of the sex where it differs.

    Raises:
        ValueError: for an entry with neither `value` nor the key of the sex.
    """
    values = entry.get("value", entry.get(sex.value))
    if values is None:
        raise ValueError(
            f"tables/{ADIPOSE}.toml: an entry over '{x}' has neither `value` nor "
            f"'{sex.value}'."
        )
    return Interpolation(
        tuple(float(a) for a in entry[x]), tuple(float(v) for v in values)
    )


@cache
def visceral_share(sex: Sex) -> AgeCurve:
    """The visceral share of the adipose tissue of a sex over age.

    Args:
        sex: the sex of the build.

    Returns:
        The share, a piecewise linear function of the postnatal age (years).
    """
    entry = _table()["visceral_share"][sex.value]
    sources, notes = _provenance(entry)
    return AgeCurve(_interpolation(entry, "age", sex), sources, notes)


@dataclass(frozen=True)
class WaistModel:
    """The waist circumference of an individual of the covariates of a build.

    The prediction a model without a measured waist uses, and the reference of
    the covariate `f_waist` of a model which has one: an individual whose
    waist is the prediction keeps the visceral share of its age and sex.

    Attributes:
        intercept: the intercept of the prediction [cm].
        bw: the coefficient of the body weight [cm/kg].
        height: the coefficient of the height [cm/cm].
        age: the coefficient of the age [cm/yr].
        sensitivity: the relative change of the visceral adipose tissue per cm
            of a deviation of the waist from the prediction [1/cm], at a fixed
            age, height and body weight.
        sources: the keys of the sources.
        notes: where the coefficients come from.
    """

    intercept: float
    bw: float
    height: float
    age: float
    sensitivity: float
    sources: tuple[str, ...]
    notes: str

    def value(self, bw: float, height: float, age: float = AGE_REFERENCE) -> float:
        """The predicted waist circumference [cm] of an individual.

        Args:
            bw: the body weight [kg].
            height: the height [cm].
            age: the age [years].

        Returns:
            The waist circumference [cm].
        """
        return self.intercept + self.bw * bw + self.height * height + self.age * age

    def formula(self, age: str) -> str:
        """The formula of the predicted waist [cm] over the covariates.

        Args:
            age: the age of the individual, the id `AGE` over age and the
                number of the reference age otherwise, dimensionless.

        Returns:
            The formula of `WAIST_base`.
        """
        out = repr(self.intercept)
        for coefficient, factor in (
            (self.bw, f"({ids.BW} / 1 kg)"),
            (self.height, f"({ids.HEIGHT} / 1 cm)"),
            (self.age, age),
        ):
            sign = "-" if coefficient < 0 else "+"
            out = f"{out} {sign} {abs(coefficient)!r} * {factor}"
        return f"({out}) * 1 cm"

    def factor(self) -> str:
        """The formula of the factor of a deviation of the waist on `f_vat`."""
        deviation = f"(({ids.WAIST} - {ids.WAIST_BASE}) / 1 cm)"
        return f"exp({self.sensitivity!r} * {deviation})"


@cache
def waist_model(sex: Sex) -> WaistModel:
    """The waist circumference of a sex and its effect on the visceral depot.

    Args:
        sex: the sex of the build.

    Returns:
        The prediction of the waist from the covariates of the build and the
        sensitivity of the visceral depot to a deviation from it.
    """
    entry = _table()["waist"][sex.value]
    coefficients = dict(
        zip(entry["coefficient"], (float(v) for v in entry["value"]), strict=True)
    )
    sensitivity = _table()["visceral_waist"][sex.value]
    sources, notes = _provenance(entry)
    more, sensitivity_notes = _provenance(sensitivity)
    return WaistModel(
        intercept=coefficients["intercept"],
        bw=coefficients["bw"],
        height=coefficients["height"],
        age=coefficients["age"],
        sensitivity=float(sensitivity["value"]),
        sources=tuple(dict.fromkeys((*sources, *more))),
        notes=f"{notes}\n\nThe sensitivity of the visceral depot: {sensitivity_notes}",
    )


def waist_factor(
    sex: Sex, waist: float, bw: float, height: float, age: float = AGE_REFERENCE
) -> float:
    """The covariate `f_waist` of a measured waist circumference.

    The waist of an individual relative to the waist of the population at its
    sex, body weight, height and age, which is what a model sets: it keeps the
    ratio over a change of the body weight. The waist is the circumference at
    the site of the WHO protocol, the midpoint between the lowest rib and the
    iliac crest.

    Args:
        sex: the sex of the individual.
        waist: the measured waist circumference [cm].
        bw: the body weight [kg].
        height: the height [cm].
        age: the age [years].

    Returns:
        The waist over the predicted waist of the covariates.

    Raises:
        ValueError: for a waist which is not positive.
    """
    if waist <= 0:
        raise ValueError(f"The waist circumference {waist} cm is not positive.")
    return waist / waist_model(sex).value(bw, height, age)


@dataclass(frozen=True)
class VisceralComposition:
    """The waist circumference and the depots of an individual, as numbers.

    What the build computes from the covariates and `f_waist`, so that a table
    of the documentation, a test and an experiment restate no equation of the
    body composition (`pkpdsim.pbpk.physiology.body_composition` gives the
    masses it builds on).

    Attributes:
        waist_base: the predicted waist of the covariates [cm], `WAIST_base`.
        waist: the waist circumference [cm], `WAIST`.
        whtr: the waist-to-height ratio, `WHtR`.
        f_visceral: the visceral share of the adipose tissue, `f_vat`.
        visceral: the visceral adipose tissue [kg], `VAT`.
        subcutaneous: the subcutaneous adipose tissue [kg], `SAT`.
    """

    waist_base: float
    waist: float
    whtr: float
    f_visceral: float
    visceral: float
    subcutaneous: float


def visceral_composition(
    sex: Sex,
    bw: float,
    height: float,
    age: float = AGE_REFERENCE,
    f_waist: float = 1.0,
) -> VisceralComposition:
    """The waist and the visceral and subcutaneous depot of an individual.

    Args:
        sex: the sex of the build.
        bw: the body weight [kg].
        height: the height [cm].
        age: the age [years]; the reference age of a build which is not over
            age, whose share is that of the reference age.
        f_waist: the waist over the predicted waist of the covariates
            (`waist_factor` for a measured waist).

    Returns:
        The waist circumference and the depots.
    """
    model = waist_model(sex)
    waist_base = model.value(bw, height, age)
    waist = f_waist * waist_base
    share = visceral_share(sex).curve.value(age) * math.exp(
        model.sensitivity * (waist - waist_base)
    )
    adipose = body_composition(sex, bw, height, age).adipose
    return VisceralComposition(
        waist_base=waist_base,
        waist=waist,
        whtr=waist / height,
        f_visceral=share,
        visceral=share * adipose,
        subcutaneous=(1.0 - share) * adipose,
    )


@cache
def obesity_factors(sex: Sex) -> Mapping[str, ObesityFactor]:
    """The factors of overweight and obesity on the physiology of a sex.

    Args:
        sex: the sex of the build.

    Returns:
        The factor by its key in the table (`cardiac_output`, the keys of
        `FLOW_FACTORS` and `BINDING_FACTORS`), a function of the body mass
        index which is 1 up to `BMI_NORMAL`.
    """
    factors = {}
    for key, entry in _table()["factor"].items():
        sources, notes = _provenance(entry)
        factors[key] = ObesityFactor(_interpolation(entry, "bmi", sex), sources, notes)
    return MappingProxyType(factors)


@cache
def obesity_protein_factors(sex: Sex) -> Mapping[tuple[Organ, str], ObesityFactor]:
    """The factors of obesity on the activity of the proteins of the organs.

    Args:
        sex: the sex of the build.

    Returns:
        The factor by the organ and the protein, per amount of tissue, as the
        factors of cirrhosis are.
    """
    factors = {}
    for organ, proteins in _table()["protein"].items():
        for protein, entry in proteins.items():
            sources, notes = _provenance(entry)
            factors[Organ(organ), protein] = ObesityFactor(
                _interpolation(entry, "bmi", sex), sources, notes
            )
    return MappingProxyType(factors)


def _default_rate() -> Value:
    """The default rate of the change of the body weight."""
    return value(ADIPOSE, "weight_change", ids.K_BW)


@dataclass(frozen=True)
class Adipose:
    """The body composition, the physiology of obesity and a change of the weight.

    Attributes:
        effect: the fractional change of the body weight at the set point
            (`f_bw`, negative for a loss), a dimensionless `Formula` over the
            ids of the top model, e.g. the Emax of the concentration of a
            GLP-1 agonist; it is an assignment rule and may therefore read the
            species and the observers of the top model, unlike the other
            formulas of a definition. `None` for a model whose body weight
            does not change, which has the body composition alone.
        k_bw: the rate (1/min) at which the body weight approaches its set
            point, a number or a fitted `Value`; the default of the table (a
            half time of 35 weeks, the measured time scale of a change of the
            body weight).
        obesity: whether the factors of overweight and obesity enter the
            physiology (the default); `False` leaves them out and keeps the
            body composition and the change of the weight alone.

    Raises:
        ValueError: for an effect which is not dimensionless, a rate which is
            not positive, or a rate without an effect.
    """

    effect: Formula | None = None
    k_bw: float | Value | None = None
    obesity: bool = True

    def __post_init__(self) -> None:
        """Check the rate and the effect.

        Raises:
            ValueError: for an effect which is not dimensionless, a rate
                without an effect or a rate which is not positive.
        """
        if self.effect is not None:
            unit = self.effect.unit
            sid = unit if isinstance(unit, str) else getattr(unit, "sid", unit)
            if str(sid) != "dimensionless":
                raise ValueError(
                    f"Adipose: the effect is in '{sid}', the fractional change "
                    "of the body weight is dimensionless."
                )
        if self.k_bw is None:
            return
        if self.effect is None:
            raise ValueError(
                "Adipose: the rate `k_bw` needs an effect which changes the weight."
            )
        rate = number(self.k_bw)
        if rate is None or rate <= 0:
            raise ValueError(f"Adipose: the rate {self.k_bw} is not positive.")

    @property
    def rate(self) -> float | Value:
        """The rate of the change of the body weight, the default of the table."""
        return _default_rate() if self.k_bw is None else self.k_bw


def _parameter(
    pid: str,
    unit: UnitType,
    name: str,
    notes: str,
    kind: str,
    *,
    constant: bool = True,
    value: float | None = None,
) -> Parameter:
    """A parameter of the top model with its annotation.

    A parameter without a `value` takes it from its initial assignment; a
    covariate states its default.
    """
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


def _constant(
    pid: str, formula: str, unit: UnitType, name: str, notes: str, kind: str
) -> list[Sbase]:
    """A constant with its initial assignment, as every derived quantity of the body.

    A quantity which depends on parameters only is an initial assignment; it
    becomes an assignment rule of the build whose body weight changes
    (`adipose_dynamic`, `pkpdsim.pbpk.builder`).
    """
    return [
        _parameter(pid, unit, name, f"{notes}\n\n`{pid} = {formula}`", kind),
        initial_assignment(pid, formula, unit, name),
    ]


def _composition(definition: "PBPKDefinition") -> list[Sbase]:
    """The lean body weight, the fat mass, the waist and the two depots."""
    sex = definition.sex
    share = visceral_share(sex)
    waist = waist_model(sex)
    citations = "\n".join(f"- {cite(s)}" for s in share.sources)
    waist_citations = "\n".join(f"- {cite(s)}" for s in waist.sources)
    fat = ids.organ_mass(Organ.ADIPOSE)
    height = f"({ids.HEIGHT} / 100 cm)"
    if definition.age:
        f_vat = share.curve.formula(ids.AGE)
        age = f"({ids.AGE} / 1 yr)"
        over = f"of the age of the {sex.value} individual"
    else:
        f_vat = f"{share.curve.value(AGE_REFERENCE)!r} dimensionless"
        age = repr(float(AGE_REFERENCE))
        over = (
            f"of the {sex.value} individual of {AGE_REFERENCE:g} years "
            "(the build is not over age)"
        )
    return [
        *_constant(
            ids.LBW,
            f"{ids.BW} - {fat}",
            U.kg,
            "lean body weight [kg]",
            "The body weight without the adipose tissue of the library, which is "
            "the remainder of the body weight after the lean organs and the "
            "blood and which therefore holds the separable adipose tissue of "
            "ICRP 89 and its unmapped tissues.",
            "body_composition",
        ),
        *_constant(
            ids.FM,
            f"{ids.BW} - {ids.FFM}",
            U.kg,
            "body fat mass [kg]",
            "The body fat mass, the body weight without the fat-free mass of "
            "Janmahasatian 2005. It is the fat of a body composition measurement "
            f"and not the mass of the adipose tissue compartment `{fat}`, which "
            "exceeds it by 2 to 4 kg over the whole range of the body mass index "
            "(the neutral lipid is about 85 % of the adipose tissue, and the "
            "compartment also holds the tissues of ICRP 89 which the library does "
            "not map).",
            "fat_mass",
        ),
        *_constant(
            ids.BFF,
            f"{ids.FM} / {ids.BW}",
            U.dimensionless,
            "body fat fraction",
            "The body fat mass over the body weight, the percent body fat of a "
            "body composition measurement divided by 100.",
            "body_composition",
        ),
        *_constant(
            ids.FMI,
            f"({ids.FM} / 1 kg) / {height}^2 * 1 kg_per_m2",
            U.kg_per_m2,
            "fat mass index [kg/m2]",
            "The body fat mass over the squared height, the fat analogue of the "
            "body mass index; the reference values and the classes of a fat "
            f"deficit and of an excess of fat are those of Kelly 2009 (DXA, "
            f"NHANES).\n\n- {cite('Kelly2009')}",
            "body_composition",
        ),
        *_constant(
            ids.WAIST_BASE,
            waist.formula(age),
            U.cm,
            "predicted waist circumference [cm]",
            "The waist circumference of the population at the sex, the body "
            "weight, the height and the age of the individual, the reference "
            f"of the covariate `{ids.F_WAIST}`: {waist.notes}\n\n"
            f"{waist_citations}",
            "body_composition",
        ),
        _parameter(
            ids.F_WAIST,
            U.dimensionless,
            "waist circumference over the predicted waist",
            "Covariate, an initial value; the waist circumference of the "
            "individual over the waist which its body weight, height and age "
            f"predict (`{ids.WAIST_BASE}`), 1.0 for the individual of the "
            "population. It is the one covariate of the distribution of the "
            "adipose tissue between the visceral and the subcutaneous depot, "
            "and the individual keeps it over a change of the body weight, so "
            f"that `{ids.WAIST}` follows the weight. "
            f"`pkpdsim.pbpk.waist_factor` gives it for a measured waist "
            "circumference (the site of the WHO protocol, the midpoint "
            "between the lowest rib and the iliac crest).",
            "body_composition",
            value=1.0,
        ),
        *_constant(
            ids.WAIST,
            f"{ids.F_WAIST} * {ids.WAIST_BASE}",
            U.cm,
            "waist circumference [cm]",
            "The waist circumference of the individual at the site of the WHO "
            "protocol, the midpoint between the lowest rib and the iliac "
            "crest. The action points of the WHO are 94 and 102 cm in men and "
            "80 and 88 cm in women.",
            "body_composition",
        ),
        *_constant(
            ids.WHTR,
            f"{ids.WAIST} / {ids.HEIGHT}",
            U.dimensionless,
            "waist-to-height ratio",
            "The waist circumference over the height, the index of central "
            "adiposity which needs no classification by sex: the boundary of "
            "0.5 is the one every guideline shares.",
            "body_composition",
        ),
        *_constant(
            ids.F_VISCERAL,
            f"({f_vat}) * {waist.factor()}",
            U.dimensionless,
            "visceral share of the adipose tissue",
            f"The share of the adipose tissue which is visceral: the share "
            f"{over}, times the factor of a deviation of the waist "
            f"circumference from the waist of the population "
            f"(`{ids.WAIST_BASE}`), which is 1 for the individual of the "
            f"population. The share of the age: {share.notes}\n\n"
            f"{citations}\n{waist_citations}",
            "body_composition",
        ),
        *_constant(
            ids.VAT,
            f"{ids.F_VISCERAL} * {fat}",
            U.kg,
            "visceral adipose tissue [kg]",
            "The visceral adipose tissue, the share `f_vat` of the adipose "
            "tissue. It is a mass of the body composition and not an organ of "
            "its own: the adipose tissue of the library is one well-stirred "
            "organ which drains into the venous blood, while the visceral depot "
            "drains through the portal vein.",
            "adipose_depot",
        ),
        *_constant(
            ids.SAT,
            f"{fat} - {ids.VAT}",
            U.kg,
            "subcutaneous adipose tissue [kg]",
            "The subcutaneous adipose tissue, the adipose tissue which is not "
            "visceral.",
            "adipose_depot",
        ),
    ]


def _weight_change(adipose: Adipose) -> list[Sbase]:
    """The set point of the body weight, its approach and the observers."""
    effect = adipose.effect
    if effect is None:  # pragma: no cover - guarded by the caller
        raise ValueError("Adipose: no effect on the body weight.")
    rate = adipose.rate
    return [
        _parameter(
            ids.BW_BASE,
            U.kg,
            "body weight at the start [kg]",
            "The body weight of the covariates, the weight at the start of the "
            "simulation, which the effect of the model moves.\n\n"
            f"`{ids.BW_BASE} = {ids.BW}`",
            "body_weight",
        ),
        StartValue(
            ids.BW_BASE,
            ids.BW,
            U.kg,
            name="body weight at the start",
            sboTerm=SBO.MATHEMATICAL_EXPRESSION,
        ),
        _parameter(
            ids.BW_EFFECT,
            U.dimensionless,
            "effect on the body weight",
            "The fractional change of the body weight at the set point, the "
            "effect of the model (negative for a loss).\n\n"
            f"`{ids.BW_EFFECT} = {effect.expression}`"
            + (f"\n\n{effect.notes}" if effect.notes else ""),
            "body_composition",
            constant=False,
        ),
        AssignmentRule(
            ids.BW_EFFECT,
            effect.expression,
            U.dimensionless,
            name="effect on the body weight",
        ),
        *_constant(
            ids.BW_SET,
            f"{ids.BW_BASE} * (1 dimensionless + {ids.BW_EFFECT})",
            U.kg,
            "body weight set point [kg]",
            "The body weight which the individual approaches: the weight at the "
            "start and the effect of the model.",
            "body_weight",
        ),
        annotate(
            Parameter(
                ids.K_BW,
                number(rate),
                U.per_min,
                constant=True,
                name="rate of the change of the body weight [1/min]",
                sboTerm=_SBO_PARAMETER,
                notes="The rate at which the body weight approaches its set "
                "point; a model fits it."
                + (
                    f"\n\n- {cite(rate.source)}"
                    if isinstance(rate, Value) and rate.source
                    else ""
                ),
            ),
            "weight_change",
        ),
        RateRule(
            ids.BW,
            f"{ids.K_BW} * ({ids.BW_SET} - {ids.BW})",
            KG_PER_MIN,
            name="change of the body weight",
        ),
        *_constant(
            "BW_change",
            f"{ids.BW} - {ids.BW_BASE}",
            U.kg,
            "body weight change [kg]",
            "The change of the body weight over the weight at the start.",
            "body_weight",
        ),
        *_constant(
            "BW_ratio",
            f"{ids.BW} / {ids.BW_BASE}",
            U.dimensionless,
            "body weight ratio",
            "The body weight over the weight at the start.",
            "body_weight",
        ),
        _parameter(
            "WAIST0",
            U.cm,
            "waist circumference at the start [cm]",
            "The waist circumference at the start of the simulation, from "
            f"which the waist of the individual departs with its weight.\n\n"
            f"`WAIST0 = {ids.WAIST}`",
            "body_composition",
        ),
        StartValue(
            "WAIST0",
            ids.WAIST,
            U.cm,
            name="waist circumference at the start",
            sboTerm=SBO.MATHEMATICAL_EXPRESSION,
        ),
        *_constant(
            "WAIST_change",
            f"{ids.WAIST} - WAIST0",
            U.cm,
            "waist circumference change [cm]",
            "The change of the waist circumference over the waist at the "
            "start, the endpoint which the trials of a weight loss report "
            "beside the weight. It follows the body weight through the "
            f"prediction `{ids.WAIST_BASE}`, at the ratio `{ids.F_WAIST}` of "
            "the individual, which is held: the waist of the prediction falls "
            "by 0.766 cm per kg in a man and by 0.790 cm per kg in a woman, "
            "while the trials of the GLP-1 receptor agonists report about "
            "0.85 cm per kg (SURMOUNT-1: about 18 cm at a loss of 20.9 % of "
            "104.8 kg), so "
            "the model under-predicts the waist of a large weight loss by "
            "about a tenth: the weight these drugs remove is more central "
            "than the weight by which the individuals of a cross section "
            "differ.",
            "body_composition",
        ),
    ]


def _check(definition: "PBPKDefinition") -> None:
    """Check that a definition may have the adipose model.

    Raises:
        ValueError: for a laboratory animal (the body composition and the
            obesity of the library are those of the human).
    """
    if definition.species is not BuildSpecies.HUMAN:
        raise ValueError(
            f"Invalid adipose model of '{definition.sid}': the species is "
            f"{definition.species.value}, the body composition and the physiology "
            "of obesity of the library are those of the human."
        )


def adipose_objects(definition: "PBPKDefinition") -> list[Sbase]:
    """The body composition and the change of the body weight of a definition.

    The factors of obesity are not here: they enter the products `fimp_<id>`
    of all conditions and are created by `pkpdsim.pbpk.impairment`.

    Args:
        definition: the definition with `adipose`.

    Returns:
        The lean body weight, the fat mass with its fraction and index and the
        visceral and subcutaneous adipose tissue, and for a model with an
        effect the set point of the body weight, its rate and the observers of
        the change; empty without `adipose`.

    Raises:
        ValueError: for a laboratory animal.
    """
    adipose = definition.adipose
    if adipose is None:
        return []
    _check(definition)
    objects = _composition(definition)
    if adipose.effect is not None:
        objects += _weight_change(adipose)
    return objects


def adipose_dynamic(definition: "PBPKDefinition") -> tuple[str, ...]:
    """The ids which change during a simulation of a model whose weight changes.

    The body weight is a state (a rate rule), so every quantity of the
    physiology, of the body composition and of the factors of obesity which
    reads it follows it (`pkpdsim.pbpk.builder`).
    """
    adipose = definition.adipose
    if adipose is None or adipose.effect is None:
        return ()
    return (ids.BW, ids.BW_EFFECT)
