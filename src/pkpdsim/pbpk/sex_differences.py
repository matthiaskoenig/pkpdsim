"""The physiology of the women relative to the men of the same age and size.

The sex is a choice of the build. The female build has the reference anatomy
of the women of ICRP 89 (organ masses, blood volume, flow fractions,
hematocrit; `pkpdsim.pbpk.physiology`), the fat-free mass of the women
(Janmahasatian 2005) and the curves of the women over age where they differ
(`pkpdsim.pbpk.age`). What differs beyond these between a woman and a man of
the same age, body weight and height are the factors of this module, each a
function of the age (`database/tables/sex_human.toml`, the forms of
`pkpdsim.pbpk.curves`), an initial assignment `fsex_<key>` of the female
build; the male build has none, its physiology is the reference of the
parameters of a substance:

- the cardiac output (`fsex_CO`): the same cardiac index of both sexes at the
  same size, which the scaling with the lean tissue misses;
- the brain fraction of the cardiac output (`fsex_Qbr`): the higher cerebral
  perfusion of the women from puberty to about 65 years;
- the plasma albumin of the binding (`fsex_albumin`): lower in the women
  from puberty to about 55 years;
- the transit time of the stomach in the fed state (`fsex_t_stomach_fed`):
  the slower emptying of a meal of the women from puberty to the menopause;
- the transit times of the large intestine (`fsex_t_colon`);
- the activity of a protein per amount of tissue (`fsex_<o>_<protein>`).

The differences arise at puberty: every factor is 1 in the children.
"""

from collections.abc import Mapping
from functools import cache
from typing import Any

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.age import AGE_REFERENCE, AgeCurve, entry_curve, provenance
from pkpdsim.pbpk.curves import Constant, Interpolation
from pkpdsim.pbpk.data import load
from pkpdsim.pbpk.organs import Organ
from pkpdsim.sex import Sex

__all__ = [
    "ALBUMIN",
    "BRAIN_FLOW",
    "CARDIAC_OUTPUT",
    "COLON_TRANSIT",
    "GASTRIC_EMPTYING_FED",
    "SEX",
    "protein_sex_factor",
    "protein_sex_factors",
    "same_cardiac_index",
    "sex_curves",
    "sex_factor",
]

SEX = "sex_human"
"""The table of the differences of the sexes in the database."""

CARDIAC_OUTPUT = ids.CO
"""The key of the factor on the cardiac output, `fsex_CO`."""

BRAIN_FLOW = ids.organ_flow(Organ.BRAIN)
"""The key of the factor on the brain fraction of the cardiac output, `fsex_Qbr`."""

ALBUMIN = "albumin"
"""The key of the factor on the plasma albumin, `fsex_albumin`."""

GASTRIC_EMPTYING_FED = "t_stomach_fed"
"""The key of the factor on the transit time of the stomach in the fed state,
`fsex_t_stomach_fed`."""

COLON_TRANSIT = "t_colon"
"""The key of the factor on the transit times of the large intestine, `fsex_t_colon`."""

_ENTRIES: dict[str, str] = {
    CARDIAC_OUTPUT: "cardiac_output",
    BRAIN_FLOW: "brain_flow",
    ALBUMIN: "albumin",
    GASTRIC_EMPTYING_FED: "gastric_emptying_fed",
    COLON_TRANSIT: "colon_transit",
}
"""The entries of the table by the key of their factor."""


def sex_factor(pid: str) -> str:
    """The id of the factor of the sex on a quantity, `fsex_CO`, `fsex_li_cyp3a4`."""
    return f"fsex_{pid}"


def _ratio(entry: Mapping[str, Any]) -> float:
    """The ratio of the women to the men of a table of the database at 35 years."""
    table = load(entry["table"])[entry["quantity"]]

    def at(sex: Sex) -> float:
        rows = table[sex.value]
        curve = Interpolation(
            tuple(float(a) for a in rows["age"]),
            tuple(float(v) for v in rows["value"]),
        )
        return curve.value(AGE_REFERENCE)

    return at(Sex.FEMALE) / at(Sex.MALE)


def same_cardiac_index(age: float) -> float:
    """The factor on the female cardiac output which makes it that of a man of the same size.

    The ratio of the cardiac output of the scaling of the male build to that
    of the female build (`pkpdsim.pbpk.physiology.scaled_cardiac_output`) at
    the body weight and height between the reference individuals of the
    sexes at the age.

    Args:
        age: the age (years).

    Returns:
        The factor (dimensionless).
    """
    # physiology imports this module (the factors of the sex)
    from pkpdsim.pbpk.physiology import reference_individual, scaled_cardiac_output

    refs = [reference_individual(sex, age) for sex in Sex]
    bw = sum(r.bw for r in refs) / len(refs)
    height = sum(r.height for r in refs) / len(refs)
    return scaled_cardiac_output(Sex.MALE, bw, height, age) / scaled_cardiac_output(
        Sex.FEMALE, bw, height, age
    )


def _factor(entry: Mapping[str, Any]) -> AgeCurve:
    """The curve of an entry of the table with its provenance."""
    sources, notes = provenance(entry)
    form = entry["form"]
    if form == "ratio":
        curve = Constant(_ratio(entry))
    elif form == "same_cardiac_index":
        ages = tuple(float(a) for a in entry["age"])
        curve = Interpolation(ages, tuple(same_cardiac_index(a) for a in ages))
    else:
        curve = entry_curve(entry)
    return AgeCurve(curve, sources, notes)


@cache
def sex_curves(sex: Sex) -> dict[str, AgeCurve]:
    """The factors of the sex on the physiology of a build, by the key of their id.

    Args:
        sex: the sex of the build.

    Returns:
        The dimensionless curves over age of the female build relative to
        the male of the same age and size (`CARDIAC_OUTPUT`, `BRAIN_FLOW`,
        `ALBUMIN`, `GASTRIC_EMPTYING_FED`, `COLON_TRANSIT`); none for the male
        build.
    """
    if sex is not Sex.FEMALE:
        return {}
    data = load(SEX)
    return {key: _factor(data[name]) for key, name in _ENTRIES.items()}


@cache
def _proteins() -> dict[str, dict[str, Any]]:
    """The entries of the proteins of the women, by organ token and protein."""
    return load(SEX).get("protein", {})


def protein_sex_factor(sex: Sex, organ: Organ, protein: str) -> AgeCurve | None:
    """The factor of the sex on the activity of a protein per amount of tissue.

    Args:
        sex: the sex of the build.
        organ: the organ.
        protein: the protein.

    Returns:
        The curve over age of the women relative to the men; `None` for the
        male build and for a protein which the sex does not change.
    """
    if sex is not Sex.FEMALE:
        return None
    entry = _proteins().get(organ.value, {}).get(protein)
    return None if entry is None else _factor(entry)


def protein_sex_factors() -> frozenset[tuple[Organ, str]]:
    """The organs and proteins with a factor of the sex."""
    return frozenset(
        (Organ(o), p) for o, proteins in _proteins().items() for p in proteins
    )
