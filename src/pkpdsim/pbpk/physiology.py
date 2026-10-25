"""The standard human physiology and its scaling to an individual.

The physiology is a set of SBML parameters and assignment rules of the top
model. The covariates of the individual (`BW`, `HEIGHT`, `SEX`, `AGE`) are
constant parameters which an experiment changes before the simulation, as it
changes `BW` today; everything else follows from them by assignment rules:

- the reference values of the sex (ICRP 89, `data/physiology_human.toml`),
  `piecewise(<female>, SEX > 0.5, <male>)`, so that one model covers both
  sexes;
- body mass index, body surface area (Du Bois 1916) and fat-free mass
  (Janmahasatian 2005);
- the volumes of the lean organs and the blood volume scale with the fat-free
  mass, adipose is the remainder of the body weight (Huisinga 2012);
- the cardiac output scales with the fat-free mass to the power 0.75
  (Huisinga 2012) and the organ flows are fractions of it, renormalised over
  the organs so that they sum to the cardiac output;
- the GFR scales with the body surface area.

An organ volume is the volume of its tissue and its blood,
`mass/density + blood*Vblood`; the body weight is the mass of the tissues
without blood and of the blood, `sum(M<o>) + Vblood*rho_blood`.
"""

from collections.abc import Iterable
from dataclasses import dataclass
from enum import StrEnum
from typing import Any

from sbmlutils.factory import AssignmentRule, Parameter, Sbase, UnitDefinition
from sbmlutils.metadata import SBO

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.data import cite, load
from pkpdsim.pbpk.organs import PORTAL_ORGANS, Organ, Pool
from pkpdsim.pbpk.units import U

__all__ = ["Sex", "physiology_objects"]


class Sex(StrEnum):
    """The sex of the reference individual."""

    MALE = "male"
    FEMALE = "female"


LEAN_ORGANS: tuple[Organ, ...] = tuple(o for o in Organ if o is not Organ.ADIPOSE)
"""The organs whose volume scales with the fat-free mass."""

AGE_DEFAULT = 40.0
"""The age of an individual if an experiment does not set it (yr)."""

_POOL_KEYS = {Pool.ARTERIAL: "arterial", Pool.VENOUS: "venous", Pool.PORTAL: "portal"}


@dataclass(frozen=True)
class _Ref:
    """A reference value of one sex.

    Attributes:
        value: the value in the unit of the model.
        sources: the keys of the sources of the value.
        notes: where the value comes from.
    """

    value: float
    sources: tuple[str, ...]
    notes: str


def _entry(table: dict[str, Any], key: str) -> _Ref:
    """Get a value of the data file as it is."""
    entry = table[key]
    return _Ref(float(entry["value"]), (entry["source"],), entry.get("notes", ""))


def _number(table: dict[str, Any], key: str) -> float:
    return float(table[key]["value"])


def _density(table: dict[str, Any]) -> float:
    """The density of an organ, 1.0 kg/l if the data state none."""
    return _number(table, "density") if "density" in table else 1.0


def _sources(table: dict[str, Any], keys: Iterable[str]) -> tuple[str, ...]:
    return tuple(table[key]["source"] for key in keys if key in table)


def _dubois(bw: float, height: float) -> float:
    """Body surface area of Du Bois 1916 (m2) of body weight (kg) and height (cm)."""
    return 0.007184 * bw**0.425 * height**0.725


def _janmahasatian(bw: float, height: float, sex: Sex) -> float:
    """Fat-free mass of Janmahasatian 2005 (kg) of body weight (kg) and height (cm)."""
    bmi = bw / (height / 100) ** 2
    if sex is Sex.MALE:
        return 9270 * bw / (6680 + 216 * bmi)
    return 9270 * bw / (8780 + 244 * bmi)


def _references(sex: Sex, data: dict[str, Any]) -> dict[str, _Ref]:
    """The reference values of one sex, by parameter id.

    Besides the ids of parameters, `rho_blood`, `rho_ad` and `blood_ad` are the
    densities of blood and adipose and the blood fraction of adipose, which
    enter the formulas of adipose as numbers.

    Raises:
        ValueError: if the lean organs and the blood of the reference
            individual leave no mass for adipose.
    """
    p = data[sex.value]
    organs, blood = p["organs"], p["blood"]
    bw, height = _number(p, "bodyweight"), _number(p, "height")
    vblood = _number(blood, "volume")
    individual = f"{sex.value} reference individual, {bw} kg, {height} cm"

    refs: dict[str, _Ref] = {
        ids.BW: _entry(p, "bodyweight"),
        ids.HEIGHT: _entry(p, "height"),
        ids.ref(ids.BSA): _Ref(
            _dubois(bw, height),
            ("DuBois1916", p["bodyweight"]["source"]),
            f"Du Bois formula of the {individual}",
        ),
        ids.ref(ids.FFM): _Ref(
            _janmahasatian(bw, height, sex),
            ("Janmahasatian2005", p["bodyweight"]["source"]),
            f"Janmahasatian formula of the {individual}",
        ),
        ids.ref(ids.CO): _entry(p, "cardiac_output"),
        ids.HCT: _entry(p, "hematocrit"),
        ids.ref(ids.VBLOOD): _entry(blood, "volume"),
        "rho_blood": _entry(blood, "density"),
    }
    gfr = _entry(p, "gfr")
    refs[ids.ref(ids.GFR)] = _Ref(gfr.value / 1000, gfr.sources, gfr.notes)

    lean_mass = 0.0
    for o in LEAN_ORGANS:
        t = organs[o]
        mass = _number(t, "mass")
        lean_mass += mass
        refs[ids.ref(ids.organ_volume(o))] = _Ref(
            mass / _density(t) + _number(t, "blood") * vblood,
            (*_sources(t, ("mass", "density", "blood")), blood["volume"]["source"]),
            f"mass/density + blood*Vblood: {mass} kg / {_density(t)} kg/l + "
            f"{_number(t, 'blood')} * {vblood} l; {t['mass']['notes']}",
        )
    refs[ids.ref(ids.MLEAN)] = _Ref(
        lean_mass,
        tuple(dict.fromkeys(organs[o]["mass"]["source"] for o in LEAN_ORGANS)),
        "sum of the tissue masses (without blood) of the lean organs",
    )

    adipose = bw - lean_mass - vblood * _number(blood, "density")
    if adipose <= 0:
        raise ValueError(
            f"The {sex.value} reference individual leaves no mass for adipose: "
            f"the lean organs and the blood exceed the body weight by "
            f"{-adipose:.3f} kg."
        )
    ad = organs[Organ.ADIPOSE]
    refs["rho_ad"] = (
        _entry(ad, "density")
        if "density" in ad
        else _Ref(1.0, ("ICRP89",), "the density of an organ without stated density")
    )
    refs["blood_ad"] = _entry(ad, "blood")

    parallel = [o for o in Organ if o is not Organ.LUNG]
    total = sum(_number(organs[o], "flow") for o in parallel)
    for o in parallel:
        flow = _entry(organs[o], "flow")
        refs[ids.flow_fraction(o)] = _Ref(
            flow.value / total,
            flow.sources,
            f"{flow.value} / {total:.4f}, the fraction renormalised over the "
            f"organs (the flows of the tissues outside the organs are "
            f"distributed proportionally); {flow.notes}",
        )

    organ_blood = sum(_number(organs[o], "blood") for o in Organ)
    for pool in (Pool.ARTERIAL, Pool.PORTAL):
        refs[ids.pool_fraction(pool)] = _entry(blood, _POOL_KEYS[pool])
    arterial, portal = _number(blood, "arterial"), _number(blood, "portal")
    unmapped = 1.0 - organ_blood - arterial - portal - _number(blood, "venous")
    refs[ids.pool_fraction(Pool.VENOUS)] = _Ref(
        1.0 - organ_blood - arterial - portal,
        (blood["venous"]["source"],),
        f"the venous pool {_number(blood, 'venous')} and the blood of the tissues "
        f"outside the organs {unmapped:.4f}, so that pools and organ blood are "
        f"the blood volume; {blood['venous']['notes']}",
    )
    return refs


def _notes(text: str, sources: Iterable[str]) -> str:
    """Notes with the citations of the sources."""
    citations = "\n".join(f"- {cite(key)}" for key in dict.fromkeys(sources))
    return f"{text}\n\nSources:\n\n{citations}"


def _parameter(
    pid: str,
    unit: UnitDefinition,
    name: str,
    text: str,
    sources: Iterable[str],
    value: float | None = None,
) -> Parameter:
    """A parameter; constant if it has a value, else the variable of a rule."""
    return Parameter(
        pid,
        value,
        unit,
        constant=value is not None,
        name=name,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        notes=_notes(text, sources),
    )


def _number_of(value: float, unit: UnitDefinition) -> str:
    """A number with its unit in a formula."""
    return f"{value!r} {unit.sid}"


def _sexed_formula(
    refs: dict[Sex, dict[str, _Ref]], pid: str, unit: UnitDefinition
) -> str:
    """The value of the sex, `piecewise(<female>, SEX > 0.5, <male>)`."""
    male = refs[Sex.MALE][pid].value
    female = refs[Sex.FEMALE][pid].value
    if male == female:
        return _number_of(male, unit)
    return (
        f"piecewise({_number_of(female, unit)}, {ids.SEX} > 0.5 dimensionless, "
        f"{_number_of(male, unit)})"
    )


def _rule(
    pid: str,
    formula: str,
    unit: UnitDefinition,
    name: str,
    text: str,
    sources: Iterable[str],
) -> list[Sbase]:
    """A parameter with its assignment rule."""
    return [
        _parameter(pid, unit, name, f"{text}\n\n`{pid} = {formula}`", sources),
        AssignmentRule(pid, formula, unit),
    ]


def _sexed(
    refs: dict[Sex, dict[str, _Ref]],
    pid: str,
    unit: UnitDefinition,
    name: str,
    text: str,
) -> list[Sbase]:
    """A reference parameter of the sex with the provenance of both values."""
    notes, sources = _provenance(refs, pid)
    return _rule(
        pid, _sexed_formula(refs, pid, unit), unit, name, f"{text}\n\n{notes}", sources
    )


def _provenance(
    refs: dict[Sex, dict[str, _Ref]], key: str
) -> tuple[str, tuple[str, ...]]:
    """The notes and sources of the male and female value of a reference."""
    male, female = refs[Sex.MALE][key], refs[Sex.FEMALE][key]
    return (
        f"- male: {male.notes}\n- female: {female.notes}",
        (*male.sources, *female.sources),
    )


def physiology_objects(sex: Sex = Sex.MALE) -> list[Sbase]:
    """The parameters and rules of the standard human physiology.

    Args:
        sex: the sex of the default individual: the default of `SEX`, `BW` and
            `HEIGHT`. The model covers both sexes, an experiment switches by
            setting `SEX`.

    Returns:
        The parameters and assignment rules of the top model.

    Raises:
        ValueError: if a reference individual leaves no mass for adipose.
    """
    data = load("physiology_human")
    refs = {s: _references(s, data) for s in Sex}
    default = refs[sex]
    ffm_ratio = f"({ids.FFM} / {ids.ref(ids.FFM)})"
    objects: list[Sbase] = []

    # covariates and factors
    objects += [
        _parameter(
            ids.BW,
            U.kg,
            "body weight [kg]",
            "Covariate; the default is the reference individual of the sex. "
            f"{default[ids.BW].notes}",
            default[ids.BW].sources,
            value=default[ids.BW].value,
        ),
        _parameter(
            ids.HEIGHT,
            U.cm,
            "height [cm]",
            "Covariate; the default is the reference individual of the sex. "
            f"{default[ids.HEIGHT].notes}",
            default[ids.HEIGHT].sources,
            value=default[ids.HEIGHT].value,
        ),
        _parameter(
            ids.SEX,
            U.dimensionless,
            "sex (0 male, 1 female)",
            "Covariate; selects the male (0) or female (1) reference values of "
            "ICRP 89 in every reference parameter.",
            ("ICRP89",),
            value=0.0 if sex is Sex.MALE else 1.0,
        ),
        _parameter(
            ids.AGE,
            U.yr,
            "age [yr]",
            "Covariate of the adult individual. The reference values are those "
            "of the adult of ICRP 89; no rule of the standard physiology "
            "depends on the age yet.",
            ("ICRP89",),
            value=AGE_DEFAULT,
        ),
        _parameter(
            ids.F_RENAL_FUNCTION,
            U.dimensionless,
            "renal function",
            "The GFR relative to the normal GFR of the individual "
            "(ICRP 89 reference scaled with the body surface area); 1.0 normal, "
            "< 1.0 reduced.",
            ("ICRP89",),
            value=1.0,
        ),
        _parameter(
            ids.F_CARDIAC_FUNCTION,
            U.dimensionless,
            "cardiac function",
            "The cardiac output relative to the normal cardiac output of the "
            "individual; it multiplies `CO` and so every blood flow; 1.0 "
            "normal, < 1.0 reduced.",
            ("ICRP89",),
            value=1.0,
        ),
    ]

    # reference values of the sex
    objects += _sexed(
        refs,
        ids.ref(ids.BSA),
        U.m2,
        "reference body surface area [m2]",
        "Body surface area of the reference individual, the reference of `GFR`.",
    )
    objects += _sexed(
        refs,
        ids.ref(ids.FFM),
        U.kg,
        "reference fat-free mass [kg]",
        "Fat-free mass of the reference individual, the reference of the "
        "lean organ volumes, the blood volume and the cardiac output.",
    )
    objects += _sexed(
        refs,
        ids.ref(ids.CO),
        U.l_per_min,
        "reference cardiac output [l/min]",
        "Cardiac output of the reference individual.",
    )
    objects += _sexed(
        refs,
        ids.ref(ids.GFR),
        U.l_per_min,
        "reference glomerular filtration rate [l/min]",
        "GFR of the reference individual.",
    )
    objects += _sexed(
        refs,
        ids.ref(ids.VBLOOD),
        U.litre,
        "reference blood volume [l]",
        "Total blood volume of the reference individual.",
    )
    objects += _sexed(
        refs,
        ids.ref(ids.MLEAN),
        U.kg,
        "reference mass of the lean organ tissues [kg]",
        "Mass of the tissues of the lean organs (without their blood) of the "
        "reference individual.",
    )
    objects += _sexed(
        refs,
        ids.HCT,
        U.dimensionless,
        "hematocrit",
        "The volume fraction of the blood cells in blood.",
    )
    for o in LEAN_ORGANS:
        objects += _sexed(
            refs,
            ids.ref(ids.organ_volume(o)),
            U.litre,
            f"reference volume {o.name.lower().replace('_', ' ')} [l]",
            "Volume of tissue and blood of the reference individual.",
        )
    for o in Organ:
        if o is Organ.LUNG:
            continue
        objects += _sexed(
            refs,
            ids.flow_fraction(o),
            U.dimensionless,
            f"blood flow fraction {o.name.lower().replace('_', ' ')}",
            "Blood flow as fraction of cardiac output"
            + (" (hepatic artery)." if o is Organ.LIVER else "."),
        )
    for pool in Pool:
        objects += _sexed(
            refs,
            ids.pool_fraction(pool),
            U.dimensionless,
            f"blood fraction {pool.name.lower()} pool",
            "Volume of the blood pool as fraction of the blood volume.",
        )

    # the individual
    bmi = f"({ids.BW} / 1 kg) / ({ids.HEIGHT} / 100 cm)^2 * 1 kg_per_m2"
    objects += _rule(
        ids.BMI,
        bmi,
        U.kg_per_m2,
        "body mass index [kg/m2]",
        "Body mass index of the individual.",
        ("Janmahasatian2005",),
    )
    objects += _rule(
        ids.BSA,
        f"0.007184 m2 * ({ids.BW} / 1 kg)^0.425 * ({ids.HEIGHT} / 1 cm)^0.725",
        U.m2,
        "body surface area [m2]",
        "Body surface area of the individual, Du Bois formula.",
        ("DuBois1916",),
    )

    def ffm(k: int, a: int) -> str:
        """The FFM per body weight of Janmahasatian, `9270 / (k + a*BMI)`."""
        return (
            f"9270 dimensionless / ({k} dimensionless + {a} dimensionless * "
            f"{ids.BMI} / 1 kg_per_m2)"
        )

    objects += _rule(
        ids.FFM,
        f"{ids.BW} * piecewise({ffm(8780, 244)}, {ids.SEX} > 0.5 dimensionless, "
        f"{ffm(6680, 216)})",
        U.kg,
        "fat-free mass [kg]",
        "Fat-free mass of the individual, the sex specific formula of "
        "Janmahasatian 2005; Huisinga 2012 use it for the lean body weight.",
        ("Janmahasatian2005", "Huisinga2012"),
    )
    objects += _rule(
        ids.VBLOOD,
        f"{ids.ref(ids.VBLOOD)} * {ffm_ratio}",
        U.litre,
        "blood volume [l]",
        "Blood volume of the individual; it scales with the fat-free mass as "
        "the lean tissues do (Huisinga 2012 scale all tissues but adipose, "
        "brain and skin with the lean body weight).",
        ("Huisinga2012",),
    )
    for o in LEAN_ORGANS:
        vid = ids.organ_volume(o)
        objects += _rule(
            vid,
            f"{ids.ref(vid)} * {ffm_ratio}",
            U.litre,
            f"volume {o.name.lower().replace('_', ' ')} [l]",
            "Volume of tissue and blood of the organ; it scales with the "
            "fat-free mass.",
            ("Huisinga2012",),
        )
    objects += _rule(
        ids.MLEAN,
        f"{ids.ref(ids.MLEAN)} * {ffm_ratio}",
        U.kg,
        "mass of the lean organ tissues [kg]",
        "Mass of the tissues of the lean organs without their blood; it scales "
        "with the fat-free mass.",
        ("Huisinga2012",),
    )
    mad = ids.organ_mass(Organ.ADIPOSE)
    rho_blood, rho_blood_sources = _provenance(refs, "rho_blood")
    objects += _rule(
        mad,
        f"{ids.BW} - {ids.MLEAN} - {ids.VBLOOD} * "
        f"{_sexed_formula(refs, 'rho_blood', U.kg_per_l)}",
        U.kg,
        "mass of adipose tissue [kg]",
        "Mass of adipose tissue without its blood, the remainder of the body "
        "weight (Huisinga 2012); it includes the tissues of ICRP 89 which are "
        f"not organs of the model. Density of blood:\n\n{rho_blood}",
        ("Huisinga2012", *rho_blood_sources),
    )
    rho_ad, rho_ad_sources = _provenance(refs, "rho_ad")
    blood_ad, blood_ad_sources = _provenance(refs, "blood_ad")
    objects += _rule(
        ids.organ_volume(Organ.ADIPOSE),
        f"{mad} / {_sexed_formula(refs, 'rho_ad', U.kg_per_l)} + "
        f"{_sexed_formula(refs, 'blood_ad', U.dimensionless)} * {ids.VBLOOD}",
        U.litre,
        "volume adipose [l]",
        f"Volume of adipose tissue and its blood. Density:\n\n{rho_ad}\n\n"
        f"Blood as fraction of the blood volume:\n\n{blood_ad}",
        (*rho_ad_sources, *blood_ad_sources),
    )
    for pool in Pool:
        objects += _rule(
            ids.pool_volume(pool),
            f"{ids.pool_fraction(pool)} * {ids.VBLOOD}",
            U.litre,
            f"volume {pool.name.lower()} blood [l]",
            "Volume of the blood pool.",
            ("ICRP89",),
        )
    objects += _rule(
        ids.CO,
        f"{ids.ref(ids.CO)} * {ffm_ratio}^0.75 * {ids.F_CARDIAC_FUNCTION}",
        U.l_per_min,
        "cardiac output [l/min]",
        "Cardiac output of the individual; it scales with the fat-free mass to "
        "the power 0.75 and with the cardiac function.",
        ("Huisinga2012",),
    )
    objects += _rule(
        ids.organ_flow(Organ.LUNG),
        ids.CO,
        U.l_per_min,
        "blood flow lung [l/min]",
        "The lung is in series and receives the cardiac output.",
        ("ICRP89",),
    )
    for o in Organ:
        if o is Organ.LUNG:
            continue
        objects += _rule(
            ids.organ_flow(o),
            f"{ids.flow_fraction(o)} * {ids.CO}",
            U.l_per_min,
            f"blood flow {o.name.lower().replace('_', ' ')} [l/min]",
            "Blood flow of the organ"
            + (", the hepatic artery." if o is Organ.LIVER else "."),
            ("ICRP89",),
        )
    objects += _rule(
        ids.QHA,
        ids.organ_flow(Organ.LIVER),
        U.l_per_min,
        "hepatic arterial blood flow [l/min]",
        "Blood flow of the hepatic artery, the flow of the liver of the flow "
        "fractions.",
        ("ICRP89",),
    )
    objects += _rule(
        ids.QPO,
        " + ".join(ids.organ_flow(o) for o in Organ if o in PORTAL_ORGANS),
        U.l_per_min,
        "portal venous blood flow [l/min]",
        "Blood flow of the portal vein, the flows of the organs which drain into it.",
        ("ICRP89",),
    )
    objects += _rule(
        ids.QH,
        f"{ids.QHA} + {ids.QPO}",
        U.l_per_min,
        "liver blood flow [l/min]",
        "Total blood flow of the liver, hepatic artery and portal vein.",
        ("ICRP89",),
    )
    objects += _rule(
        ids.GFR,
        f"{ids.ref(ids.GFR)} * {ids.BSA} / {ids.ref(ids.BSA)} * {ids.F_RENAL_FUNCTION}",
        U.l_per_min,
        "glomerular filtration rate [l/min]",
        "GFR of the individual: the reference GFR, which is roughly "
        "proportional to body size (ICRP 89 Section 8.3.1), scaled with the "
        "body surface area, i.e. the same GFR per 1.73 m2, times the renal "
        "function.",
        ("ICRP89",),
    )
    return objects
