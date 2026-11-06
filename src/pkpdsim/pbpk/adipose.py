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
    "adipose_dynamic",
    "adipose_objects",
    "obesity_factor",
    "obesity_factors",
    "obesity_protein_factors",
    "visceral_share",
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
) -> Parameter:
    """A parameter of the top model with its annotation."""
    return annotate(
        Parameter(
            pid,
            None,
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
    """The lean body weight, the fat mass and the visceral and subcutaneous depot."""
    share = visceral_share(definition.sex)
    citations = "\n".join(f"- {cite(s)}" for s in share.sources)
    fat = ids.organ_mass(Organ.ADIPOSE)
    height = f"({ids.HEIGHT} / 100 cm)"
    if definition.age:
        f_vat = share.curve.formula(ids.AGE)
        over = f"of the age of the {definition.sex.value} individual"
    else:
        f_vat = f"{share.curve.value(AGE_REFERENCE)!r} dimensionless"
        over = (
            f"of the {definition.sex.value} individual of {AGE_REFERENCE:g} years "
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
            ids.F_VISCERAL,
            f_vat,
            U.dimensionless,
            "visceral share of the adipose tissue",
            f"The share of the adipose tissue which is visceral, a function "
            f"{over}: {share.notes}\n\n{citations}",
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
