"""The physiology over age: growth, maturation and aging of a sex.

The reference values and the factors of the physiology are curves of the
covariate `AGE` (`pkpdsim.pbpk.curves`), written once and evaluated as a
number (`Curve.value`, the reference individual of an age, the tables of the
documentation, the virtual individuals) and as the formula of an initial
assignment of the model (`Curve.formula`). An individual of any age from
birth to 100 years follows from the scaling of the physiology to its body
weight and height (Huisinga 2012) relative to the reference individual of its
age and sex:

- growth (`reference_curves`): the reference children of ICRP 89 at birth and
  1, 5, 10 and 15 years (`database/tables/children_human.toml`) and the
  adult of ICRP 89 (`physiology_human.toml`) at `ADULT_AGE`, interpolated
  linearly in age: the reference body weight and height, the tissue masses
  of the lean organs, the blood volume and the cardiac output; Huisinga 2012
  take the reference individual of the age class, the interpolation between
  the classes is a choice of the library;
- aging of the adults: the adult reference value times an aging factor of
  the age (`database/tables/age_human.toml`) at a fixed body size, 1 at
  `AGE_REFERENCE` (35 years, the adult of ICRP 89), so the default individual
  of ICRP 89 keeps its values: the liver and the kidney mass, the cardiac
  output;
- the factors of the functions (`factor_curves`), 1 at 35 years: the GFR
  (maturation of Rhodin 2009 and aging), the hematocrit (neonates and NHANES),
  the albumin and the alpha-1-acid glycoprotein of the binding, the fat-free
  mass of the children (Al-Sallami 2015) and the blood flows (the brain of
  the children, the liver and the kidney of the elderly, `flow_factors`);
- the ontogeny of the enzymes and transporters (`protein_factor`,
  `database/tables/ontogeny.toml`).

The differences of the sexes beyond the reference anatomy of each sex are
`pkpdsim.pbpk.sex_differences`.

The age enters the formulas as `AGE` in years; every curve is constant
beyond the ages of its data.
"""

from collections.abc import Mapping
from dataclasses import dataclass
from functools import cache
from typing import Any

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.curves import (
    Affine,
    Constant,
    Curve,
    Hill,
    Interpolation,
    Logarithmic,
    Product,
    Quotient,
    Relative,
    Spliced,
)
from pkpdsim.pbpk.data import load
from pkpdsim.pbpk.organs import PORTAL_ORGANS, Organ
from pkpdsim.sex import Sex

__all__ = [
    "ADULT_AGE",
    "AGE_REFERENCE",
    "LIVER_GROUP",
    "AgeCurve",
    "FlowFactors",
    "age_factor",
    "entry_curve",
    "factor_curves",
    "flow_factors",
    "flow_fraction_curves",
    "protein_factor",
    "protein_factors",
    "provenance",
    "reference_curves",
    "rest_denominator",
]

ADULT_AGE = 20.0
"""The age (years) at which the growth ends: the reference values of the
children reach the adult of ICRP 89 at this age (ICRP 89 Section 1.4.2 para
14: adults are 20 to 50 years)."""

AGE_REFERENCE = 35.0
"""The age (years) of the adult reference values of ICRP 89 (Section 1.4.2
para 14: "an average age of 35 years is used"), the default age of the
individual and the age at which every factor of age is 1."""

CHILDREN = "children_human"
"""The table of the reference children in the database."""

AGE = "age_human"
"""The table of the maturation and the aging in the database."""

ONTOGENY = "ontogeny"
"""The table of the ontogeny of the proteins in the database."""


@dataclass(frozen=True)
class AgeCurve:
    """A quantity over age with its provenance.

    Attributes:
        curve: the value over age, in the unit of the quantity or
            dimensionless for a factor.
        sources: the keys of the sources.
        notes: where the values come from.
    """

    curve: Curve
    sources: tuple[str, ...]
    notes: str


def age_factor(pid: str) -> str:
    """The id of the factor of the age on a quantity, `fage_GFR`, `fage_li_cyp3a4`."""
    return f"fage_{pid}"


def _table(name: str) -> dict[str, Any]:
    """A table of the database."""
    return load(name)


def _points(entry: Mapping[str, Any]) -> Interpolation:
    """The interpolation of an entry of the form `points`."""
    return Interpolation(
        tuple(float(a) for a in entry["age"]),
        tuple(float(v) for v in entry["value"]),
    )


def entry_curve(entry: Mapping[str, Any]) -> Curve:
    """The curve of an entry of age_human.toml, ontogeny.toml or sex_human.toml.

    Args:
        entry: the entry with its `form` (`points`, `hill`, `log` or
            `constant`) and the numbers of the form.

    Returns:
        The curve over age.

    Raises:
        ValueError: for an unknown form or variable.
    """
    form = entry["form"]
    if form == "constant":
        return Constant(float(entry["value"]))
    if form == "points":
        return _points(entry)
    if form == "hill":
        variable = entry["variable"]
        if variable not in ("age", "pma"):
            raise ValueError(f"Unknown variable '{variable}' of a hill curve.")
        return Hill(
            x50=float(entry["x50"]),
            n=float(entry["n"]),
            base=float(entry.get("base", 0.0)),
            top=float(entry.get("top", 1.0)),
            variable=variable,
            x_min=float(entry.get("x_min", 0.0)),
        )
    if form == "log":
        return Logarithmic(float(entry["a"]), float(entry["b"]))
    raise ValueError(f"Unknown form '{form}' of a curve.")


def provenance(*entries: Mapping[str, Any]) -> tuple[tuple[str, ...], str]:
    """The sources and the notes of entries of the database.

    Args:
        entries: entries with a `source`, optional further `sources` and
            `notes`.

    Returns:
        The keys of the sources, each once in their order, and the notes
        joined.
    """
    sources = tuple(
        dict.fromkeys(k for e in entries for k in (e["source"], *e.get("sources", ())))
    )
    return sources, "; ".join(e.get("notes", "") for e in entries)


def _nhanes(entry: Mapping[str, Any], sex: Sex) -> Interpolation:
    """The values of a quantity of the NHANES table of a sex, the top-coded row moved.

    The row of 80 years and over is placed at `top_coded_age` of the entry.
    """
    table = _table(entry["table"])[entry["quantity"]][sex.value]
    ages = [float(a) for a in table["age"]]
    if ages[-1] == 80.0:
        ages[-1] = float(entry["top_coded_age"])
    return Interpolation(tuple(ages), tuple(float(v) for v in table["value"]))


def _aging(sex: Sex, key: str | None) -> Curve:
    """The aging factor of a reference value of the adults, 1 at `AGE_REFERENCE`."""
    if key is None:
        return Constant(1.0)
    entry = _table(AGE)[key]["aging"]
    if "form" not in entry:
        entry = entry[sex.value]
    return Relative(entry_curve(entry), AGE_REFERENCE)


_AGING_KEYS: dict[str, str] = {
    ids.ref(ids.CO): "cardiac_output",
    ids.ref(ids.organ_mass(Organ.LIVER)): "liver_mass",
    ids.ref(ids.organ_mass(Organ.KIDNEY)): "kidney_mass",
}
"""The reference values with an aging factor, by their id."""


def _growth(child: Mapping[str, Any], adult: float, aging: Curve, unit: str) -> Curve:
    """The curve of a reference value: the children, the adult at `ADULT_AGE`, the aging.

    Args:
        child: the entry of the children (values at the ages of the table).
        adult: the adult reference value (at `AGE_REFERENCE`).
        aging: the aging factor of the adult, 1 at `AGE_REFERENCE`.
        unit: the unit of the quantity, which the children's entry must have.

    Raises:
        ValueError: if the unit of the children is another.
    """
    if child["unit"] != unit:
        raise ValueError(
            f"{CHILDREN}.toml: the unit '{child['unit']}' is not the unit '{unit}'."
        )
    ages = tuple(float(a) for a in _table(CHILDREN)["ages"])
    values = tuple(float(v) for v in child["value"])
    at_adult = adult * aging.value(ADULT_AGE)
    growth = Interpolation((*ages, ADULT_AGE), (*values, at_adult))
    if isinstance(aging, Constant):
        return growth
    return Spliced((ADULT_AGE,), (growth, Product((Constant(adult), aging))))


@cache
def reference_curves(sex: Sex) -> dict[str, AgeCurve]:
    """The reference values over age of a sex, by the id of their parameter.

    The reference body weight and height (`BW_ref`, `HEIGHT_ref`), the
    tissue masses of the lean organs (`M<o>_ref`), the blood volume
    (`Vblood_ref`) and the cardiac output (`CO_ref`).

    Args:
        sex: the sex.

    Returns:
        The curves in the units of the model (kg, cm, l, l/min).
    """
    children = _table(CHILDREN)[sex.value]
    adult = load("physiology_human")[sex.value]
    curves: dict[str, AgeCurve] = {}

    def add(
        pid: str, child: Mapping[str, Any], adult_entry: Mapping[str, Any], unit: str
    ) -> None:
        key = _AGING_KEYS.get(pid)
        aging = _aging(sex, key)
        sources = [child["source"], adult_entry["source"]]
        notes = (
            f"children (newborn, 1, 5, 10, 15 years): {child['notes']}; the adult "
            f"from {ADULT_AGE:g} years: {adult_entry.get('notes', '')}; linear "
            "between the ages"
        )
        if key is not None:
            entry = _table(AGE)[key]["aging"]
            entry = entry if "form" in entry else entry[sex.value]
            sources.append(entry["source"])
            notes += f"; times the aging of the adults: {entry['notes']}"
        curves[pid] = AgeCurve(
            _growth(child, float(adult_entry["value"]), aging, unit),
            tuple(dict.fromkeys(sources)),
            notes,
        )

    add(ids.ref(ids.BW), children["bodyweight"], adult["bodyweight"], "kg")
    add(ids.ref(ids.HEIGHT), children["height"], adult["height"], "cm")
    add(ids.ref(ids.VBLOOD), children["blood_volume"], adult["blood"]["volume"], "l")
    add(ids.ref(ids.CO), children["cardiac_output"], adult["cardiac_output"], "l/min")
    for o in Organ:
        if o is Organ.ADIPOSE:
            continue
        add(
            ids.ref(ids.organ_mass(o)),
            children["organs"][o],
            adult["organs"][o]["mass"],
            "kg",
        )
    return curves


def _hematocrit(sex: Sex) -> AgeCurve:
    """The hematocrit relative to 35 years: neonates, then NHANES."""
    data = _table(AGE)["hematocrit"]
    newborn, adult = data["newborn"], data["adult"]
    nhanes = _nhanes(adult, sex)
    points = _points(newborn)
    curve = Interpolation(
        (*points.ages, *nhanes.ages), (*points.values, *nhanes.values)
    )
    sources, notes = provenance(newborn, adult)
    return AgeCurve(Relative(curve, AGE_REFERENCE), sources, notes)


def _albumin(sex: Sex) -> AgeCurve:
    """The albumin relative to 35 years: Johnson 2006 below 12 years, then NHANES."""
    data = _table(AGE)["albumin"]
    children, adult = data["children"], data["adult"]
    nhanes = _nhanes(adult, sex)
    young = nhanes.ages[0]
    johnson = entry_curve(children)
    scaled = Product((Constant(nhanes.value(young) / johnson.value(young)), johnson))
    curve = Spliced((young,), (scaled, nhanes))
    sources, notes = provenance(children, adult)
    return AgeCurve(Relative(curve, AGE_REFERENCE), sources, notes)


def _relative(entry: Mapping[str, Any], at: float = AGE_REFERENCE) -> AgeCurve:
    """The curve of an entry relative to its value at an age."""
    sources, notes = provenance(entry)
    return AgeCurve(Relative(entry_curve(entry), at), sources, notes)


def _gfr(sex: Sex) -> AgeCurve:
    """The GFR at a fixed body surface area: maturation times aging."""
    data = _table(AGE)["gfr"]
    maturation, aging = data["maturation"], data["aging"][sex.value]
    sources, _notes = provenance(maturation, aging)
    curve = Product(
        (
            Relative(entry_curve(maturation), AGE_REFERENCE),
            Relative(entry_curve(aging), AGE_REFERENCE),
        )
    )
    return AgeCurve(
        curve, sources, f"maturation: {maturation['notes']}; aging: {aging['notes']}"
    )


def _fat_free_mass(sex: Sex) -> AgeCurve:
    """The multiplier of the fat-free mass of the children, 1 from `ADULT_AGE`."""
    entry = _table(AGE)["fat_free_mass"][sex.value]
    sources, notes = provenance(entry)
    curve = Spliced(
        (ADULT_AGE,), (Relative(entry_curve(entry), ADULT_AGE), Constant(1.0))
    )
    return AgeCurve(curve, sources, notes)


@cache
def factor_curves(sex: Sex) -> dict[str, AgeCurve]:
    """The factors of the age on the functions of a sex, by the id they multiply.

    `GFR` (the GFR at renal function 1 at a fixed body surface area), `HCT`
    (the hematocrit), `FFM` (the fat-free mass at a fixed body weight and
    height), `albumin` and `agp` (the binding proteins); every factor is 1 at
    `AGE_REFERENCE`.

    Args:
        sex: the sex.

    Returns:
        The dimensionless curves.
    """
    data = _table(AGE)
    return {
        ids.GFR: _gfr(sex),
        ids.HCT: _hematocrit(sex),
        ids.FFM: _fat_free_mass(sex),
        "albumin": _albumin(sex),
        "agp": _relative(data["agp"]["maturation"]),
    }


LIVER_GROUP: tuple[Organ, ...] = (Organ.LIVER, *sorted(PORTAL_ORGANS))
"""The organs whose flows follow the hepatic blood flow: the hepatic artery and
the portal organs."""


@dataclass(frozen=True)
class FlowFactors:
    """The factors of the age on the fractions of the cardiac output of a sex.

    The flow fractions of the adults of ICRP 89 (`fQ<o>`) are multiplied by a
    factor of the age, so that they sum to 1 at every age:

    - the brain by `brain`, its flow per mass of the age (Chiron 1992) times
      its mass over the cardiac output, both relative to 35 years (the
      children's brain takes up to half of the cardiac output);
    - every other organ shares the rest in the proportions of the adults,
      `rest = (1 - fQbr*brain)/(1 - fQbr)`;
    - of these, the liver and the portal organs by `liver` (the hepatic flow
      of the elderly follows the liver volume) and the kidney by `kidney`
      (the renal plasma flow over the cardiac output of the elderly);
    - the remaining organs by `other`, which closes the sum:
      `(1 - fQbr - sum fQ_s*h_s)/(1 - fQbr - sum fQ_s)` over the liver group
      and the kidney.

    Attributes:
        brain: the factor of the brain.
        liver: the factor of the liver group (before `rest`).
        kidney: the factor of the kidney (before `rest`).
        sources: the keys of the sources.
        notes: where the factors come from.
    """

    brain: Curve
    liver: Curve
    kidney: Curve
    sources: tuple[str, ...]
    notes: str


@cache
def flow_factors(sex: Sex) -> FlowFactors:
    """The factors of the age on the flow fractions of a sex.

    Args:
        sex: the sex.

    Returns:
        The curves of the brain, the liver group and the kidney.
    """
    data = _table(AGE)
    refs = reference_curves(sex)
    perfusion = data["brain_perfusion"]["children"]
    mbr = refs[ids.ref(ids.organ_mass(Organ.BRAIN))].curve
    co = refs[ids.ref(ids.CO)].curve
    brain = Product(
        (
            Relative(entry_curve(perfusion), AGE_REFERENCE),
            Relative(Quotient(mbr, co), AGE_REFERENCE),
        )
    )
    liver = Relative(entry_curve(data["liver_flow"]["aging"]), AGE_REFERENCE)
    kidney = Quotient(
        Relative(entry_curve(data["kidney_flow"]["aging"]), AGE_REFERENCE),
        _aging(sex, "cardiac_output"),
    )
    sources, notes = provenance(
        perfusion,
        data["liver_flow"]["aging"],
        data["kidney_flow"]["aging"],
        data["cardiac_output"]["aging"],
    )
    return FlowFactors(brain, liver, kidney, sources, notes)


def rest_denominator(brain: float, fractions: list[float]) -> float:
    """`1 - fQbr - sum fQ_s`, subtracted in order as the formula of the model does."""
    d = 1.0 - brain
    for f in fractions:
        d -= f
    return d


def flow_fraction_curves(
    sex: Sex, fractions: Mapping[Organ, float]
) -> dict[Organ, Curve]:
    """The flow fractions over age of the organs in parallel.

    Args:
        sex: the sex.
        fractions: the fractions of the adults, `fQ<o>` without the lung,
            summing to 1.

    Returns:
        The fraction of every organ over age, the brain of the female build
        with its factor of the sex (`pkpdsim.pbpk.sex_differences`); they sum
        to 1 at every age.
    """
    # sex_differences imports this module (the curves of the table)
    from pkpdsim.pbpk.sex_differences import BRAIN_FLOW, sex_curves

    f = flow_factors(sex)
    fbr = fractions[Organ.BRAIN]
    brain = f.brain
    if (sex_brain := sex_curves(sex).get(BRAIN_FLOW)) is not None:
        brain = Product((f.brain, sex_brain.curve))
    rest = Affine(1.0 / (1.0 - fbr), ((-fbr / (1.0 - fbr), brain),))
    special = dict.fromkeys(LIVER_GROUP, f.liver) | {Organ.KIDNEY: f.kidney}
    denominator = rest_denominator(fbr, [fractions[o] for o in special])
    other = Quotient(
        Affine(
            1.0 - fbr,
            tuple((-fractions[o], h) for o, h in special.items()),
        ),
        Constant(denominator),
    )
    out: dict[Organ, Curve] = {}
    for o, fq in fractions.items():
        if o is Organ.BRAIN:
            out[o] = Product((Constant(fq), brain))
        elif o in special:
            out[o] = Product((Constant(fq), rest, special[o]))
        else:
            out[o] = Product((Constant(fq), rest, other))
    return out


@cache
def _ontogeny() -> dict[str, Any]:
    """The table of the ontogeny of the proteins."""
    return _table(ONTOGENY)


def protein_factor(organ: Organ, protein: str) -> AgeCurve | None:
    """The ontogeny of a protein in an organ relative to 35 years.

    Args:
        organ: the organ.
        protein: the protein.

    Returns:
        The curve, `None` for a protein without ontogeny in the organ.
    """
    entry = _ontogeny()["ontogeny"].get(organ.value, {}).get(protein)
    if entry is None:
        return None
    return _relative(entry)


def protein_factors() -> frozenset[tuple[Organ, str]]:
    """The organs and proteins with a factor of age (their ontogeny)."""
    return frozenset(
        (Organ(o), p)
        for o, proteins in _ontogeny()["ontogeny"].items()
        for p in proteins
    )
