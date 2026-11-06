"""The standard human physiology and its scaling to an individual.

The physiology is a set of SBML parameters of the top model. The covariates
of the individual (`BW`, `HEIGHT`, `AGE`) and the factors of organ function
are constant parameters which an experiment changes before the simulation,
as it changes `BW` today; everything else follows from them by initial
assignments on constant parameters, which sbmlsim evaluates again after the
changes of an experiment (a value which no rule recomputes during the
integration keeps the compiled model small and fast to load):

- the reference values of the sex of the build (ICRP 89,
  `database/tables/physiology_human.toml`) and of the age: the reference
  values are curves of `AGE` (`pkpdsim.pbpk.age`, the reference children of
  ICRP 89 and the aging of the adults, 1 at 35 years); the sex is chosen when
  the model is built, a study of both sexes uses one model per sex;
- the default individual (`BW`, `HEIGHT`, `AGE`) is the reference individual
  of the sex of the population of the build (`pkpdsim.pbpk.populations`,
  ICRP 89 by default); the reference anatomy stays ICRP 89 for every
  population, the organs of another population's individual follow from the
  scaling;
- body mass index, body surface area (Du Bois 1916) and fat-free mass
  (Janmahasatian 2005);
- the tissues scale as in Huisinga 2012 (eqs. 9, 10, 12): brain is constant,
  skin scales with the body surface area, the other lean organs and the blood
  with the fat-free mass less brain and skin (`SV_tis`); adipose is the
  remainder of the body weight;
- the cardiac output scales as the tissues, `SV_tis` (Huisinga 2012 eq. 13),
  and the organ flows are fractions of it, renormalised over the organs so
  that they sum to the cardiac output;
- the GFR scales with the body surface area and matures and declines with
  the age (`fage_GFR`); the hematocrit, the fat-free mass and the flow
  fractions carry their factors of the age (`fage_<id>`).

An organ volume is the volume of its tissue and its blood,
`M<o>/density + blood*Vblood`; the body weight is the mass of the tissues
without blood and of the blood, `sum(M<o>) + Vblood*rho_blood`. The vascular
fraction of an organ volume is its blood, `blood*Vblood/V<o>`, so that the
vascular spaces of the organs and the pools are the blood volume for every
individual; the interstitial fraction is the interstitial fraction of the
tissue without blood (Schmitt 2008) times the tissue share `1 - f_vas`.

A covariate or a factor of organ function is an initial value: changed before
the simulation (a change of an experiment, `init(...)` in roadrunner), every
quantity which depends on it follows; changed during a simulation, nothing
follows (unless the definition names it `dynamic`, see
`pkpdsim.pbpk.builder`).

The impaired organ function (`pkpdsim.pbpk.impairment`) multiplies cardiac
output, the blood flows, the GFR, the hematocrit and the kidney mass by its
factors (`factors` of `physiology_objects`); without them the physiology is
the healthy one.
"""

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from functools import cache
from types import MappingProxyType
from typing import TYPE_CHECKING, Any

from sbmlutils.factory import Parameter, Sbase, UnitDefinition
from sbmlutils.metadata import SBO

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.age import (
    ADULT_AGE,
    AGE_REFERENCE,
    LIVER_GROUP,
    age_factor,
    factor_curves,
    flow_factors,
    reference_curves,
    rest_denominator,
)
from pkpdsim.pbpk.annotation import annotate, initial_assignment
from pkpdsim.pbpk.data import cite, load
from pkpdsim.pbpk.organs import PORTAL_ORGANS, Organ, Pool
from pkpdsim.pbpk.sex_differences import (
    BRAIN_FLOW,
    CARDIAC_OUTPUT,
    COLON_TRANSIT,
    GASTRIC_EMPTYING_FED,
    sex_curves,
    sex_factor,
)
from pkpdsim.pbpk.units import U
from pkpdsim.sex import Sex

if TYPE_CHECKING:
    from pkpdsim.pbpk.populations import Population

__all__ = [
    "ADULT_AGE",
    "AGE_REFERENCE",
    "DEFAULT_POPULATION",
    "ReferenceIndividual",
    "Sex",
    "normal_gfr",
    "physiology_objects",
    "reference_fractions",
    "reference_hematocrit",
    "reference_individual",
    "reference_tissue_volumes",
]


LEAN_ORGANS: tuple[Organ, ...] = tuple(o for o in Organ if o is not Organ.ADIPOSE)
"""The organs whose tissue mass scales with the individual; adipose is the rest."""

DEFAULT_POPULATION = "european_icrp"
"""The population of a model which does not choose one (`pkpdsim.pbpk.populations`)."""

MODEL_DEFAULT = "Model default without a source."
"""The notes of a parameter whose value is a choice of the model, not data."""

LIVER_DEFAULTS: dict[str, float] = {ids.F_SHUNT: 0.0, ids.F_FUNCTIONAL: 1.0}
"""The healthy liver: no portosystemic shunt, all tissue functional. The
values at runtime are those of the hepatic impairment (`pkpdsim.impairment`),
which sets them as an experiment sets a covariate."""

_POOL_KEYS = {Pool.ARTERIAL: "arterial", Pool.VENOUS: "venous", Pool.PORTAL: "portal"}

_CONVERSIONS: dict[tuple[str, str], float] = {("ml/min", "l/min"): 1e-3}
"""Factors from a unit of the data to the unit of the model."""


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


def _number(table: dict[str, Any], key: str, unit: str) -> float:
    """Get a number of the data file in the unit of the model.

    Raises:
        ValueError: if the unit of the data cannot be converted to `unit`.
    """
    entry = table[key]
    if entry["unit"] == unit:
        return float(entry["value"])
    try:
        return float(entry["value"]) * _CONVERSIONS[(entry["unit"], unit)]
    except KeyError:
        raise ValueError(
            f"'{key}' is given in '{entry['unit']}', the model needs '{unit}'."
        ) from None


def _entry(table: dict[str, Any], key: str, unit: str) -> _Ref:
    """Get a value of the data file in the unit of the model."""
    entry = table[key]
    return _Ref(_number(table, key, unit), (entry["source"],), entry.get("notes", ""))


def _density(table: dict[str, Any]) -> _Ref:
    """The density of an organ, 1.0 kg/l if the data state none."""
    if "density" in table:
        return _entry(table, "density", "kg/l")
    return _Ref(
        1.0,
        (),
        "1.0 kg/l, the convention of physiology_human.toml for an organ "
        "without stated density",
    )


def _dubois(bw: float, height: float) -> float:
    """Body surface area of Du Bois 1916 (m2) of body weight (kg) and height (cm)."""
    return 0.007184 * bw**0.425 * height**0.725


def _janmahasatian(bw: float, height: float, sex: Sex) -> float:
    """Fat-free mass of Janmahasatian 2005 (kg) of body weight (kg) and height (cm)."""
    bmi = bw / (height / 100) ** 2
    if sex is Sex.MALE:
        return 9270 * bw / (6680 + 216 * bmi)
    return 9270 * bw / (8780 + 244 * bmi)


def _dubois_formula(bw: str, height: str) -> str:
    """The Du Bois body surface area of the ids of body weight (kg) and height (cm)."""
    return f"0.007184 m2 * ({bw} / 1 kg)^0.425 * ({height} / 1 cm)^0.725"


def _ffm_formula(bw: str, height: str, sex: Sex, bmi: str | None = None) -> str:
    """The fat-free mass of Janmahasatian 2005, `BW * 9270 / (k + a*BMI)`.

    Args:
        bw: the id of the body weight (kg).
        height: the id of the height (cm).
        sex: the sex of the coefficients.
        bmi: the id of the body mass index (kg/m2), else computed of both.
    """
    k, a = (6680, 216) if sex is Sex.MALE else (8780, 244)
    index = (
        f"{bmi} / 1 kg_per_m2"
        if bmi is not None
        else f"({bw} / 1 kg) / ({height} / 100 cm)^2"
    )
    return (
        f"{bw} * 9270 dimensionless / ({k} dimensionless + {a} dimensionless * {index})"
    )


def _references(sex: Sex, data: dict[str, Any]) -> dict[str, _Ref]:
    """The reference values of one sex, by parameter id.

    Besides the ids of parameters, `rho_blood`, `rho_<o>` and `blood_<o>` are
    the densities of blood and of the organs and the blood of the organs as
    fraction of the blood volume, which enter the formulas as numbers.

    Raises:
        ValueError: if the lean organs and the blood of the reference
            individual leave no mass for adipose, or a unit of the data is not
            the unit of the model.
    """
    p = data[sex.value]
    organs, blood = p["organs"], p["blood"]
    bw, height = _number(p, "bodyweight", "kg"), _number(p, "height", "cm")
    individual = f"{sex.value} reference individual, {bw} kg, {height} cm"
    refs: dict[str, _Ref] = {
        ids.BW: _entry(p, "bodyweight", "kg"),
        ids.HEIGHT: _entry(p, "height", "cm"),
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
        ids.ref(ids.CO): _entry(p, "cardiac_output", "l/min"),
        ids.ref(ids.GFR): _entry(p, "gfr", "l/min"),
        ids.HCT: _entry(p, "hematocrit", "dimensionless"),
        ids.ref(ids.VBLOOD): _entry(blood, "volume", "l"),
        "rho_blood": _entry(blood, "density", "kg/l"),
    }
    for o in Organ:
        refs[f"rho_{o}"] = _density(organs[o])
        refs[f"blood_{o}"] = _entry(organs[o], "blood", "dimensionless")
        f_vas = _number(organs[o], "f_vas", "dimensionless")
        f_int = _entry(organs[o], "f_int", "dimensionless")
        refs[f"fint_{o}"] = _Ref(
            f_int.value / (1.0 - f_vas),
            f_int.sources,
            f"f_int / (1 - f_vas) = {f_int.value} / (1 - {f_vas}), the "
            f"interstitial fraction of the tissue without blood; {f_int.notes}",
        )
    for o in LEAN_ORGANS:
        refs[ids.ref(ids.organ_mass(o))] = _entry(organs[o], "mass", "kg")

    lean = sum(refs[ids.ref(ids.organ_mass(o))].value for o in LEAN_ORGANS)
    adipose = bw - lean - refs[ids.ref(ids.VBLOOD)].value * refs["rho_blood"].value
    if adipose <= 0:
        raise ValueError(
            f"The {sex.value} reference individual leaves no mass for adipose: "
            f"the lean organs and the blood exceed the body weight by "
            f"{-adipose:.3f} kg."
        )

    parallel = [o for o in Organ if o is not Organ.LUNG]
    total = sum(_number(organs[o], "flow", "dimensionless") for o in parallel)
    for o in parallel:
        flow = _entry(organs[o], "flow", "dimensionless")
        refs[ids.flow_fraction(o)] = _Ref(
            flow.value / total,
            flow.sources,
            f"{flow.value} / {total:.4f}, the fraction renormalised over the "
            f"organs (the flows of the tissues outside the organs are "
            f"distributed proportionally); {flow.notes}",
        )

    for pool in (Pool.ARTERIAL, Pool.PORTAL):
        refs[ids.pool_fraction(pool)] = _entry(blood, _POOL_KEYS[pool], "dimensionless")
    organ_blood = sum(refs[f"blood_{o}"].value for o in Organ)
    arterial = refs[ids.pool_fraction(Pool.ARTERIAL)].value
    portal = refs[ids.pool_fraction(Pool.PORTAL)].value
    venous = _entry(blood, "venous", "dimensionless")
    refs[ids.pool_fraction(Pool.VENOUS)] = _Ref(
        1.0 - organ_blood - arterial - portal,
        venous.sources,
        f"the venous pool {venous.value} and the blood of the tissues outside "
        f"the organs {1.0 - organ_blood - arterial - portal - venous.value:.4f}, "
        f"so that pools and organ blood are the blood volume; {venous.notes}",
    )
    return refs


def _fraction(pid: str, value: float) -> float:
    """Check that a default of the physiology is a fraction in [0, 1].

    Raises:
        ValueError: if it is not.
    """
    if not 0.0 <= value <= 1.0:
        raise ValueError(f"The default of '{pid}' must be in [0, 1], got {value}.")
    return value


def _notes(text: str, sources: Iterable[str]) -> str:
    """Notes with the citations of the sources."""
    keys = list(dict.fromkeys(sources))
    if not keys:
        return text
    citations = "\n".join(f"- {cite(key)}" for key in keys)
    return f"{text}\n\nSources:\n\n{citations}"


def _parameter(
    pid: str,
    unit: UnitDefinition,
    name: str,
    text: str,
    sources: Iterable[str],
    value: float | None = None,
    *,
    kind: str,
    organ: Organ | None = None,
    pool: Pool | None = None,
) -> Parameter:
    """A constant parameter, of a value or of its initial assignment.

    It is annotated as `kind` (`pkpdsim.pbpk.annotation`) of the organ or pool.
    """
    return annotate(
        Parameter(
            pid,
            value,
            unit,
            constant=True,
            name=name,
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes=_notes(text, sources),
        ),
        kind,
        organ=organ,
        pool=pool,
    )


def _initial(
    pid: str,
    formula: str,
    unit: UnitDefinition,
    name: str,
    text: str,
    sources: Iterable[str],
    *,
    kind: str,
    organ: Organ | None = None,
    pool: Pool | None = None,
) -> list[Sbase]:
    """A constant parameter with its initial assignment, annotated as `kind`."""
    return [
        _parameter(
            pid,
            unit,
            name,
            f"{text}\n\n`{pid} = {formula}`",
            sources,
            kind=kind,
            organ=organ,
            pool=pool,
        ),
        initial_assignment(pid, formula, unit, name),
    ]


class _Physiology:
    """The SBML objects of the physiology, from the reference values of one sex."""

    def __init__(
        self,
        sex: Sex,
        factors: Mapping[str, str],
        population: "Population",
        default: bool,
    ) -> None:
        self.sex = sex
        self.refs = _references(sex, load("physiology_human"))
        self.curves = reference_curves(sex)
        self.factors_of_age = factor_curves(sex)
        self.factors_of_sex = sex_curves(sex)
        self.population = population
        self.default = default
        self.factors = factors
        self.applied: set[str] = set()

    def initial(
        self,
        pid: str,
        formula: str,
        unit: UnitDefinition,
        name: str,
        text: str,
        sources: Iterable[str],
        *,
        kind: str,
        organ: Organ | None = None,
        pool: Pool | None = None,
    ) -> list[Sbase]:
        """A constant parameter with its initial assignment, times its factor of impairment."""
        factor = self.factors.get(pid)
        if factor is not None:
            formula = f"({formula}) * {factor}"
            text += (
                f" Times `{factor}`, the factor of the impaired organ function "
                "(`pkpdsim.pbpk.impairment`)."
            )
            self.applied.add(pid)
        return _initial(
            pid, formula, unit, name, text, sources, kind=kind, organ=organ, pool=pool
        )

    def number(self, key: str, unit: UnitDefinition) -> str:
        """The reference value of the sex in a formula, `1.04 kg_per_l`."""
        return f"{self.refs[key].value!r} {unit.sid}"

    def provenance(self, *keys: str) -> tuple[str, tuple[str, ...]]:
        """The notes and sources of reference values of the sex."""
        lines, sources = [], []
        for key in keys:
            ref = self.refs[key]
            label = f"{key}, " if len(keys) > 1 else ""
            lines.append(f"- {label}{self.sex.value}: {ref.notes}")
            sources += ref.sources
        return "\n".join(lines), tuple(sources)

    def reference(
        self,
        pid: str,
        unit: UnitDefinition,
        name: str,
        text: str,
        key: str | None = None,
        *,
        kind: str,
        organ: Organ | None = None,
        pool: Pool | None = None,
    ) -> list[Sbase]:
        """A constant reference parameter of the sex with its provenance.

        `key` is the reference value of the parameter, `pid` by default.
        """
        key = key or pid
        notes, sources = self.provenance(key)
        return [
            _parameter(
                pid,
                unit,
                name,
                f"{text}\n\n{notes}",
                sources,
                value=self.refs[key].value,
                kind=kind,
                organ=organ,
                pool=pool,
            )
        ]

    def covariates(self) -> list[Sbase]:
        """The covariates of the individual and the factors of organ function."""
        individual = self.population.reference[self.sex]
        default = (
            "the reference individual of the sex"
            if self.default
            else f"the reference individual of the population {self.population.name}"
        )

        def covariate(
            pid: str, unit: UnitDefinition, name: str, kind: str
        ) -> Parameter:
            """A covariate with the value of the reference individual."""
            v = individual[pid]
            origin = (
                f"the default is {default}"
                if v.source is not None
                else "the population gives no value, the default is a choice of "
                "the model"
            )
            return _parameter(
                pid,
                unit,
                name,
                f"Covariate, an initial value; {origin}. {v.notes or ''}".rstrip(),
                (v.source,) if v.source is not None else (),
                value=v.value,
                kind=kind,
            )

        return [
            covariate(ids.BW, U.kg, "body weight [kg]", kind="body_weight"),
            covariate(ids.HEIGHT, U.cm, "height [cm]", kind="height"),
            covariate(ids.AGE, U.yr, "age [yr]", kind="age"),
            _parameter(
                ids.F_RENAL_FUNCTION,
                U.dimensionless,
                "renal function",
                "The GFR relative to the normal GFR of the individual; 1.0 "
                "normal, < 1.0 reduced. An initial value, as the covariates. "
                f"{MODEL_DEFAULT}",
                (),
                value=1.0,
                kind="organ_function",
                organ=Organ.KIDNEY,
            ),
            _parameter(
                ids.F_CARDIAC_FUNCTION,
                U.dimensionless,
                "cardiac function",
                "The cardiac output relative to the normal cardiac output of the "
                "individual; it multiplies `CO` and so every blood flow; 1.0 "
                "normal, < 1.0 reduced. An initial value, as the covariates. "
                f"{MODEL_DEFAULT}",
                (),
                value=1.0,
                kind="organ_function",
                organ=Organ.HEART,
            ),
            _parameter(
                ids.F_SHUNT,
                U.dimensionless,
                "portosystemic shunt",
                "The fraction in [0, 1] of the hepatic arterial and portal venous "
                "inflow which bypasses the liver into the venous blood; 0.0 "
                "healthy, the hepatic impairment sets it. An initial value, as "
                f"the covariates. {MODEL_DEFAULT}",
                (),
                value=_fraction(ids.F_SHUNT, LIVER_DEFAULTS[ids.F_SHUNT]),
                kind="shunt",
                organ=Organ.LIVER,
            ),
            _parameter(
                ids.F_FUNCTIONAL,
                U.dimensionless,
                "functional liver fraction",
                "The functional fraction in [0, 1] of the liver tissue, which "
                "scales the tissue volume of the liver and the capacity of its "
                "processes; 1.0 healthy, the hepatic impairment sets it. An "
                f"initial value, as the covariates. {MODEL_DEFAULT}",
                (),
                value=_fraction(ids.F_FUNCTIONAL, LIVER_DEFAULTS[ids.F_FUNCTIONAL]),
                kind="organ_function",
                organ=Organ.LIVER,
            ),
        ]

    def curve(
        self,
        pid: str,
        unit: UnitDefinition,
        name: str,
        text: str,
        *,
        kind: str,
        organ: Organ | None = None,
    ) -> list[Sbase]:
        """A reference value of the sex over age, the initial assignment of its curve."""
        c = self.curves[pid]
        return _initial(
            pid,
            c.curve.formula(ids.AGE, str(unit.sid)),
            unit,
            name,
            f"{text} A function of the age (`pkpdsim.pbpk.age`): {c.notes}",
            c.sources,
            kind=kind,
            organ=organ,
        )

    def factor(
        self,
        key: str,
        name: str,
        *,
        kind: str,
        organ: Organ | None = None,
    ) -> list[Sbase]:
        """The factor of the age on a quantity, `fage_<key>`, 1 at 35 years."""
        c = self.factors_of_age[key]
        return _initial(
            age_factor(key),
            c.curve.formula(ids.AGE),
            U.dimensionless,
            f"age factor {name}",
            f"The factor of the age on the {name} of the {self.sex.value} "
            f"individual at a fixed body size, 1 at {AGE_REFERENCE:g} years "
            f"(`pkpdsim.pbpk.age`): {c.notes}",
            c.sources,
            kind="age_factor",
            organ=organ,
        )

    def with_sex(self, key: str, formula: str) -> str:
        """A formula times the factor of the sex on a quantity, if the build has one."""
        if key in self.factors_of_sex:
            return f"{formula} * {sex_factor(key)}"
        return formula

    def sex_factors(self) -> list[Sbase]:
        """The factors of the sex of the female build, `fsex_<key>`, curves of the age.

        The cardiac output, the brain fraction of the cardiac output and the
        transit times of the stomach in the fed state and of the large
        intestine (inputs of the intestine); the factors of the albumin and of
        the proteins are those of the impairment
        (`pkpdsim.pbpk.sex_differences`). None in the male build.
        """
        names = {
            CARDIAC_OUTPUT: "cardiac output",
            BRAIN_FLOW: "brain fraction of the cardiac output",
            GASTRIC_EMPTYING_FED: "transit time of the stomach in the fed state",
            COLON_TRANSIT: "transit times of the large intestine",
        }
        objects: list[Sbase] = []
        for key, name in names.items():
            c = self.factors_of_sex.get(key)
            if c is None:
                continue
            objects += _initial(
                sex_factor(key),
                c.curve.formula(ids.AGE),
                U.dimensionless,
                f"sex factor {name}",
                f"The factor of the sex on the {name} of the {self.sex.value} "
                "individual relative to a man of the same age, body weight and "
                f"height (`pkpdsim.pbpk.sex_differences`): {c.notes}",
                c.sources,
                kind="sex_factor",
            )
        return objects

    def flow_fractions(self) -> list[Sbase]:
        """The flow fractions over age: the adults of ICRP 89 times the factors of the age.

        The brain by `fage_Qbr` (and in the female build by its factor of the
        sex `fsex_Qbr`), the other organs share the rest in the proportions of
        the adults (`fage_Qrest`), the liver group and the kidney with their
        factors `fage_Qli` and `fage_Qki`, the remaining organs close the sum
        (`fage_Qother`); `pkpdsim.pbpk.age.FlowFactors`.
        """
        f = flow_factors(self.sex)
        text = (
            f"of the {self.sex.value} individual at a fixed body size, 1 at {AGE_REFERENCE:g} years "
            f"(`pkpdsim.pbpk.age.FlowFactors`): {f.notes}"
        )
        brain, _kidney = ids.flow_fraction(Organ.BRAIN), ids.flow_fraction(Organ.KIDNEY)
        qbr, qli, qki = (
            age_factor(ids.organ_flow(o))
            for o in (Organ.BRAIN, Organ.LIVER, Organ.KIDNEY)
        )
        rest, other = age_factor("Qrest"), age_factor("Qother")
        objects: list[Sbase] = []
        for pid, curve, name in (
            (qbr, f.brain, "brain flow"),
            (qli, f.liver, "liver flow"),
            (qki, f.kidney, "kidney flow"),
        ):
            objects += _initial(
                pid,
                curve.formula(ids.AGE),
                U.dimensionless,
                f"age factor {name}",
                f"The factor of the age on the {name} fraction {text}",
                f.sources,
                kind="age_factor",
            )
        fbr = self.refs[brain].value
        qbr = self.with_sex(BRAIN_FLOW, qbr)
        special = dict.fromkeys(LIVER_GROUP, qli) | {Organ.KIDNEY: qki}
        denominator = rest_denominator(
            fbr, [self.refs[ids.flow_fraction(o)].value for o in special]
        )
        objects += _initial(
            rest,
            f"(1 dimensionless - {fbr!r} dimensionless * {qbr}) / "
            f"{1.0 - fbr!r} dimensionless",
            U.dimensionless,
            "age factor of the organs besides the brain",
            "The share of the cardiac output of the organs besides the brain "
            f"relative to the adults, `(1 - fQbr*{qbr})/(1 - fQbr)`.",
            f.sources,
            kind="age_factor",
        )
        terms = " - ".join(
            f"{self.refs[ids.flow_fraction(o)].value!r} dimensionless * {h}"
            for o, h in special.items()
        )
        objects += _initial(
            other,
            f"(1 dimensionless - {fbr!r} dimensionless - {terms}) / "
            f"{denominator!r} dimensionless",
            U.dimensionless,
            "age factor of the organs without a factor of their own",
            "The factor which closes the sum of the fractions: the share of "
            "the organs besides the brain, the liver group and the kidney, "
            "relative to the adults.",
            f.sources,
            kind="age_factor",
        )
        for o in Organ:
            if o is Organ.LUNG:
                continue
            pid = ids.flow_fraction(o)
            if o is Organ.BRAIN:
                factor = qbr
            elif o in special:
                factor = f"{rest} * {special[o]}"
            else:
                factor = f"{rest} * {other}"
            notes, sources = self.provenance(pid)
            objects += _initial(
                pid,
                f"{self.number(pid, U.dimensionless)} * {factor}",
                U.dimensionless,
                f"blood flow fraction {_label(o)}",
                "Blood flow as fraction of cardiac output"
                + (" (hepatic artery)" if o is Organ.LIVER else "")
                + f": the fraction of the adults times `{factor}`.\n\n{notes}",
                sources,
                kind="flow_fraction",
                organ=o,
            )
        return objects

    def references(self) -> list[Sbase]:
        """The reference values of the sex, those of the size over age."""
        objects = self.curve(
            ids.ref(ids.BW),
            U.kg,
            "reference body weight [kg]",
            "Body weight of the reference individual of the age.",
            kind="body_weight",
        )
        objects += self.curve(
            ids.ref(ids.HEIGHT),
            U.cm,
            "reference height [cm]",
            "Height of the reference individual of the age.",
            kind="height",
        )
        bw, height = ids.ref(ids.BW), ids.ref(ids.HEIGHT)
        objects += _initial(
            ids.ref(ids.BSA),
            _dubois_formula(bw, height),
            U.m2,
            "reference body surface area [m2]",
            "Body surface area of the reference individual of the age, Du Bois "
            "formula; the reference of the skin.",
            ("DuBois1916",),
            kind="bsa",
        )
        objects += self.factor(ids.FFM, "fat-free mass", kind="lean_mass")
        objects += _initial(
            ids.ref(ids.FFM),
            f"{_ffm_formula(bw, height, self.sex)} * {age_factor(ids.FFM)}",
            U.kg,
            "reference fat-free mass [kg]",
            f"Fat-free mass of the reference individual of the age, the "
            f"{self.sex.value} formula of Janmahasatian 2005 times the factor of "
            "the children (Al-Sallami 2015).",
            ("Janmahasatian2005", "AlSallami2015"),
            kind="lean_mass",
        )
        objects += self.curve(
            ids.ref(ids.CO),
            U.l_per_min,
            "reference cardiac output [l/min]",
            "Cardiac output of the reference individual of the age.",
            kind="cardiac_output",
        )
        objects += self.reference(
            ids.ref(ids.GFR),
            U.l_per_min,
            "reference glomerular filtration rate [l/min]",
            "GFR of the adult reference individual.",
            kind="gfr",
            organ=Organ.KIDNEY,
        )
        objects += self.curve(
            ids.ref(ids.VBLOOD),
            U.litre,
            "reference blood volume [l]",
            "Total blood volume of the reference individual of the age.",
            kind="blood_volume",
        )
        objects += self.factor(ids.HCT, "hematocrit", kind="hematocrit")
        notes, sources = self.provenance(ids.HCT)
        objects += _initial(
            ids.ref(ids.HCT),
            f"{self.number(ids.HCT, U.dimensionless)} * {age_factor(ids.HCT)}",
            U.dimensionless,
            "reference hematocrit",
            "The volume fraction of the blood cells in blood of the healthy "
            f"individual of the age: the adult value times `{age_factor(ids.HCT)}`."
            f"\n\n{notes}",
            sources,
            kind="hematocrit",
        )
        objects += self.initial(
            ids.HCT,
            ids.ref(ids.HCT),
            U.dimensionless,
            "hematocrit",
            "The volume fraction of the blood cells in blood, the reference of "
            "the age.",
            sources,
            kind="hematocrit",
        )
        for o in LEAN_ORGANS:
            objects += self.curve(
                ids.ref(ids.organ_mass(o)),
                U.kg,
                f"reference tissue mass {_label(o)} [kg]",
                "Mass of the organ tissue without its blood of the reference "
                "individual of the age.",
                kind="organ_mass",
                organ=o,
            )
        objects += self.flow_fractions()
        for pool in Pool:
            objects += self.reference(
                ids.pool_fraction(pool),
                U.dimensionless,
                f"blood fraction {pool.name.lower()} pool",
                "Volume of the blood pool as fraction of the blood volume.",
                kind="pool_fraction",
                pool=pool,
            )
        return objects

    def body(self) -> list[Sbase]:
        """Body mass index, body surface area, fat-free mass and `SV_tis`."""
        objects = self.initial(
            ids.BMI,
            f"({ids.BW} / 1 kg) / ({ids.HEIGHT} / 100 cm)^2 * 1 kg_per_m2",
            U.kg_per_m2,
            "body mass index [kg/m2]",
            "Body mass index of the individual.",
            ("Janmahasatian2005",),
            kind="bmi",
        )
        objects += self.initial(
            ids.BSA,
            _dubois_formula(ids.BW, ids.HEIGHT),
            U.m2,
            "body surface area [m2]",
            "Body surface area of the individual, Du Bois formula.",
            ("DuBois1916",),
            kind="bsa",
        )

        objects += self.initial(
            ids.FFM,
            f"{_ffm_formula(ids.BW, ids.HEIGHT, self.sex, bmi=ids.BMI)} * "
            f"{age_factor(ids.FFM)}",
            U.kg,
            "fat-free mass [kg]",
            f"Fat-free mass of the individual, the {self.sex.value} formula of "
            "Janmahasatian 2005 times the factor of the children (Al-Sallami "
            "2015, 1 from 20 years); Huisinga 2012 (eq. 6) use it for the lean "
            "body weight.",
            ("Janmahasatian2005", "AlSallami2015", "Huisinga2012"),
            kind="lean_mass",
        )
        br, sk = ids.organ_mass(Organ.BRAIN), ids.organ_mass(Organ.SKIN)
        objects += self.initial(
            ids.SV_TIS,
            f"({ids.FFM} - {br} - {sk}) / "
            f"({ids.ref(ids.FFM)} - {ids.ref(br)} - {ids.ref(sk)})",
            U.dimensionless,
            "scaling factor of the lean tissues",
            "Scaling factor of the lean organs but brain and skin, of the blood "
            "volume and of the cardiac output: the fat-free mass less brain and "
            "skin relative to the reference individual (Huisinga 2012 eq. 12, "
            "with the tissue masses without blood).",
            ("Huisinga2012",),
            kind="scaling_factor",
        )
        return objects

    def volumes(self) -> list[Sbase]:
        """Tissue masses, organ volumes, blood volume and blood pools."""
        objects = self.initial(
            ids.VBLOOD,
            f"{ids.ref(ids.VBLOOD)} * {ids.SV_TIS}",
            U.litre,
            "blood volume [l]",
            "Blood volume of the individual; it scales with the lean tissues "
            "(Huisinga 2012 eq. 12, blood is one of the remaining tissues).",
            ("Huisinga2012",),
            kind="blood_volume",
        )
        for o in LEAN_ORGANS:
            mid = ids.organ_mass(o)
            if o is Organ.BRAIN:
                formula = ids.ref(mid)
                text = "it is constant (Huisinga 2012 eq. 9)."
            elif o is Organ.SKIN:
                formula = f"{ids.ref(mid)} * {ids.BSA} / {ids.ref(ids.BSA)}"
                text = (
                    "it scales with the body surface area (Huisinga 2012 eq. 10; "
                    "Huisinga use the Mosteller formula of the body surface area, "
                    "the model the Du Bois formula)."
                )
            else:
                formula = f"{ids.ref(mid)} * {ids.SV_TIS}"
                text = "it scales with the lean tissues (Huisinga 2012 eq. 12)."
            objects += self.initial(
                mid,
                formula,
                U.kg,
                f"tissue mass {_label(o)} [kg]",
                f"Mass of the organ tissue without its blood; {text}",
                ("Huisinga2012",),
                kind="organ_mass",
                organ=o,
            )
        objects += self.initial(
            ids.MLEAN,
            " + ".join(ids.organ_mass(o) for o in LEAN_ORGANS),
            U.kg,
            "mass of the lean organ tissues [kg]",
            "Mass of the tissues of the lean organs without their blood.",
            ("ICRP89",),
            kind="lean_mass",
        )
        mad = ids.organ_mass(Organ.ADIPOSE)
        notes, sources = self.provenance("rho_blood")
        objects += self.initial(
            mad,
            f"{ids.BW} - {ids.MLEAN} - {ids.VBLOOD} * "
            f"{self.number('rho_blood', U.kg_per_l)}",
            U.kg,
            "tissue mass adipose [kg]",
            "Mass of adipose tissue without its blood, the remainder of the body "
            "weight; it includes the tissues of ICRP 89 which are not organs of "
            f"the model. Density of blood:\n\n{notes}",
            sources,
            kind="organ_mass",
            organ=Organ.ADIPOSE,
        )
        for o in Organ:
            notes, sources = self.provenance(f"rho_{o}", f"blood_{o}")
            objects += self.initial(
                ids.organ_volume(o),
                f"{ids.organ_mass(o)} / {self.number(f'rho_{o}', U.kg_per_l)} + "
                f"{self.number(f'blood_{o}', U.dimensionless)} * {ids.VBLOOD}",
                U.litre,
                f"volume {_label(o)} [l]",
                "Volume of the organ tissue and its blood, with the density of "
                "the tissue and the blood of the organ as fraction of the blood "
                f"volume:\n\n{notes}",
                sources,
                kind="organ_volume",
                organ=o,
            )
        for pool in Pool:
            objects += self.initial(
                ids.pool_volume(pool),
                f"{ids.pool_fraction(pool)} * {ids.VBLOOD}",
                U.litre,
                f"volume {pool.name.lower()} blood [l]",
                "Volume of the blood pool.",
                ("ICRP89",),
                kind="pool_volume",
                pool=pool,
            )
        return objects

    def fractions(self) -> list[Sbase]:
        """The vascular and interstitial fractions of the organ volumes."""
        objects: list[Sbase] = []
        for o in Organ:
            f_vas, f_int = ids.vascular_fraction(o), ids.interstitial_fraction(o)
            notes, sources = self.provenance(f"blood_{o}")
            objects += self.initial(
                f_vas,
                f"{self.number(f'blood_{o}', U.dimensionless)} * {ids.VBLOOD} / "
                f"{ids.organ_volume(o)}",
                U.dimensionless,
                f"vascular fraction {_label(o)}",
                "Vascular fraction of the organ volume: the blood of the organ "
                "(fraction of the blood volume) relative to the organ volume; it "
                "is the `f_vas` of physiology_human.toml for the reference "
                f"individual. Blood of the organ:\n\n{notes}",
                sources,
                kind="vascular_fraction",
                organ=o,
            )
            notes, sources = self.provenance(f"fint_{o}")
            objects += self.initial(
                f_int,
                f"{self.number(f'fint_{o}', U.dimensionless)} * (1 dimensionless - "
                f"{f_vas})",
                U.dimensionless,
                f"interstitial fraction {_label(o)}",
                "Interstitial fraction of the organ volume: the interstitial "
                "fraction of the tissue without blood times the tissue share of "
                "the organ volume; it is the `f_int` of physiology_human.toml for "
                f"the reference individual.\n\n{notes}",
                sources,
                kind="interstitial_fraction",
                organ=o,
            )
        return objects

    def flows(self) -> list[Sbase]:
        """Cardiac output, blood flows and GFR.

        An organ flow is its fraction of the cardiac output without
        impairment `CO_base`, so that a factor of impairment on the cardiac
        output (the sum of the impaired flows) is not a part of every flow;
        without impairment `CO = CO_base`.
        """
        objects = _initial(
            ids.CO_BASE,
            self.with_sex(
                CARDIAC_OUTPUT,
                f"{ids.ref(ids.CO)} * {ids.SV_TIS} * {ids.F_CARDIAC_FUNCTION}",
            ),
            U.l_per_min,
            "cardiac output without impairment [l/min]",
            "Cardiac output of the individual without the impaired organ "
            "function; it scales as the lean tissues (Huisinga 2012 eq. 13, "
            "SQ_co = SV_tis) and with the cardiac function"
            + (
                ", and with the factor of the sex (the same cardiac index of "
                "women and men of the same size)"
                if CARDIAC_OUTPUT in self.factors_of_sex
                else ""
            )
            + ". The organ flows are fractions of it.",
            ("Huisinga2012",),
            kind="cardiac_output",
        )
        objects += self.initial(
            ids.CO,
            ids.CO_BASE,
            U.l_per_min,
            "cardiac output [l/min]",
            "Cardiac output of the individual, the sum of the organ flows.",
            ("Huisinga2012",),
            kind="cardiac_output",
        )
        objects += self.initial(
            ids.organ_flow(Organ.LUNG),
            ids.CO,
            U.l_per_min,
            "blood flow lung [l/min]",
            "The lung is in series and receives the cardiac output.",
            ("ICRP89",),
            kind="organ_flow",
            organ=Organ.LUNG,
        )
        for o in Organ:
            if o is Organ.LUNG:
                continue
            objects += self.initial(
                ids.organ_flow(o),
                f"{ids.flow_fraction(o)} * {ids.CO_BASE}",
                U.l_per_min,
                f"blood flow {_label(o)} [l/min]",
                "Blood flow of the organ"
                + (", the hepatic artery" if o is Organ.LIVER else "")
                + f": its fraction of the cardiac output without impairment "
                f"`{ids.CO_BASE}`.",
                ("ICRP89",),
                kind="organ_flow",
                organ=o,
            )
        objects += self.initial(
            ids.QHA,
            ids.organ_flow(Organ.LIVER),
            U.l_per_min,
            "hepatic arterial blood flow [l/min]",
            "Blood flow of the hepatic artery, the flow of the liver of the flow "
            "fractions.",
            ("ICRP89",),
            kind="hepatic_artery_flow",
            organ=Organ.LIVER,
        )
        objects += self.initial(
            ids.QPO,
            " + ".join(ids.organ_flow(o) for o in Organ if o in PORTAL_ORGANS),
            U.l_per_min,
            "portal venous blood flow [l/min]",
            "Blood flow of the portal vein, the flows of the organs which drain "
            "into it.",
            ("ICRP89",),
            kind="portal_flow",
            pool=Pool.PORTAL,
        )
        objects += self.initial(
            ids.QH,
            f"{ids.QHA} + {ids.QPO}",
            U.l_per_min,
            "liver blood flow [l/min]",
            "Total blood flow of the liver, hepatic artery and portal vein.",
            ("ICRP89",),
            kind="hepatic_flow",
            organ=Organ.LIVER,
        )
        objects += self.factor(ids.GFR, "GFR", kind="age_factor", organ=Organ.KIDNEY)
        objects += self.initial(
            ids.GFR_BASE,
            f"{ids.ref(ids.GFR)} * {ids.BSA} / {self.number(ids.ref(ids.BSA), U.m2)} * "
            f"{age_factor(ids.GFR)}",
            U.l_per_min,
            "glomerular filtration rate at renal function 1 [l/min]",
            "GFR of the individual at renal function 1: the reference GFR, which "
            "is roughly proportional to body size (ICRP 89 Section 8.3.1), scaled "
            "with the body surface area, i.e. the same GFR per 1.73 m2, times the "
            f"factor of the age `{age_factor(ids.GFR)}` (the maturation of the "
            "children and the decline of the adults). A model "
            "which needs the GFR of the kidney without chronic kidney disease "
            "reads it (defined at renal function 0, where the GFR is 0).",
            ("ICRP89",),
            kind="gfr_base",
            organ=Organ.KIDNEY,
        )
        objects += _initial(
            ids.GFR,
            f"{ids.GFR_BASE} * {ids.F_RENAL_FUNCTION}",
            U.l_per_min,
            "glomerular filtration rate [l/min]",
            f"GFR of the individual: the GFR at renal function 1 `{ids.GFR_BASE}` "
            "times the renal function.",
            ("ICRP89",),
            kind="gfr",
            organ=Organ.KIDNEY,
        )
        return objects


def _label(organ: Organ) -> str:
    """The organ in a name, `small intestine`."""
    return organ.name.lower().replace("_", " ")


def physiology_objects(
    sex: Sex = Sex.MALE,
    factors: Mapping[str, str] | None = None,
    population: str = DEFAULT_POPULATION,
) -> list[Sbase]:
    """The parameters and initial assignments of the standard human physiology.

    Args:
        sex: the sex of the model: its reference values and the default
            individual (`BW`, `HEIGHT`). A study of both sexes uses one model
            per sex.
        factors: the factors of the impaired organ function by the id of the
            quantity they multiply (`pkpdsim.pbpk.impairment.physiology_factors`):
            the formula of the quantity is multiplied by the factor (a formula);
            the hematocrit becomes `HCT = HCT_ref * factor`. None: the healthy
            physiology.
        population: the population (`pkpdsim.pbpk.populations`) whose reference
            individual of the sex is the default individual (`BW`, `HEIGHT`,
            `AGE`); the reference anatomy stays ICRP 89.

    Returns:
        The constant parameters and their initial assignments of the top model.

    Raises:
        ValueError: if the reference individual leaves no mass for adipose, a
            unit of the data is not the unit of the model, a factor names
            no quantity of the physiology with an initial assignment, or the
            population is unknown.
    """
    # populations imports this module (Sex, the defaults), so it is imported here
    from pkpdsim.pbpk.populations import population as get_population

    physiology = _Physiology(
        sex,
        factors or {},
        get_population(population),
        population == DEFAULT_POPULATION,
    )
    objects = [
        *physiology.covariates(),
        *physiology.references(),
        *physiology.sex_factors(),
        *physiology.body(),
        *physiology.volumes(),
        *physiology.fractions(),
        *physiology.flows(),
    ]
    unused = sorted(set(physiology.factors) - physiology.applied)
    if unused:
        raise ValueError(
            f"Factors of impairment for quantities which the physiology does not "
            f"compute: {unused}."
        )
    return objects


def reference_fractions(organ: Organ) -> dict[Sex, tuple[float, float]]:
    """The vascular and interstitial fraction of an organ volume of the reference individuals.

    Args:
        organ: the organ.

    Returns:
        `(f_vas, f_int)` of physiology_human.toml by sex.
    """
    data = load("physiology_human")
    return {
        sex: (
            _number(data[sex.value]["organs"][organ], "f_vas", "dimensionless"),
            _number(data[sex.value]["organs"][organ], "f_int", "dimensionless"),
        )
        for sex in Sex
    }


def normal_gfr(sex: Sex, age: float = AGE_REFERENCE) -> float:
    """The GFR per 1.73 m2 body surface area of a healthy individual.

    Args:
        sex: the sex of the individual.
        age: the age (years); the GFR per 1.73 m2 matures in the children and
            declines in the adults (`fage_GFR`), 35 years the adult of ICRP 89.

    Returns:
        The reference GFR scaled to 1.73 m2 with the reference body surface
        area (l/min) times the factor of the age, the GFR of the normal renal
        function of the age.
    """
    refs = _references(sex, load("physiology_human"))
    factor = factor_curves(sex)[ids.GFR].curve.value(age)
    return refs[ids.ref(ids.GFR)].value * 1.73 / refs[ids.ref(ids.BSA)].value * factor


def _at_age(sex: Sex, age: float) -> dict[str, float]:
    """The reference values of a sex at an age, by the id of their parameter."""
    return {pid: c.curve.value(age) for pid, c in reference_curves(sex).items()}


def reference_tissue_volumes(
    sex: Sex, age: float = AGE_REFERENCE
) -> dict[Organ, float]:
    """The tissue volumes of the organs of a reference individual.

    The tissue volume of an organ is its volume without its blood,
    `V<o>*(1 - f_vas_<o>) = M<o>/density`; adipose is the remainder of the
    body weight, as in the model.

    Args:
        sex: the sex of the reference individual.
        age: the age of the reference individual (years).

    Returns:
        The tissue volume (l) by organ.
    """
    refs = _references(sex, load("physiology_human"))
    at = _at_age(sex, age)
    masses = {o: at[ids.ref(ids.organ_mass(o))] for o in LEAN_ORGANS}
    masses[Organ.ADIPOSE] = (
        at[ids.ref(ids.BW)]
        - sum(masses.values())
        - at[ids.ref(ids.VBLOOD)] * refs["rho_blood"].value
    )
    return {o: masses[o] / refs[f"rho_{o}"].value for o in Organ}


def reference_hematocrit(sex: Sex, age: float = AGE_REFERENCE) -> float:
    """The hematocrit of a reference individual.

    Args:
        sex: the sex of the reference individual.
        age: the age of the reference individual (years).

    Returns:
        The hematocrit of physiology_human.toml times the factor of the age.
    """
    factor = factor_curves(sex)[ids.HCT].curve.value(age)
    return _references(sex, load("physiology_human"))[ids.HCT].value * factor


@dataclass(frozen=True)
class ReferenceIndividual:
    """The reference individual of a sex (ICRP 89, `database/tables/physiology_human.toml`).

    The individual whose anatomy the physiology scales, the default
    individual of the populations of ICRP 89: a model which states a value of
    the hand-written models per litre of an organ or per kg of body weight
    converts it with these numbers instead of copies of them.

    Attributes:
        sex: the sex.
        bw: the body weight [kg], `BW` of the default individual.
        height: the height [cm], `HEIGHT` of the default individual.
        gfr: the healthy GFR [l/min], the GFR at renal function 1 of the
            default individual (`GFR_base`, the reference GFR `GFR_ref` of
            the sex at its own body surface area), which `healthy_gfr` of
            the build of the sex gives in ml/min; the GFR of a component
            which a model keeps outside the library (the fluid model of
            hctz) reads it instead of a copy.
        volumes: the organ volumes with their blood `V<o>` [l],
            `M<o>/density + blood*Vblood` as in the model.
        tissue_volumes: the tissue volumes without blood [l]
            (`reference_tissue_volumes`).
        intracellular_volumes: the volumes of the cells of the organs [l],
            `V<o>*(1 - f_vas_<o> - f_int_<o>)` as in the model (the tissue
            volume without its interstitial fraction), the intracellular space
            `Vic` of the full form (the hepatocytes of the liver, 1.3859 l of a
            man).
        lumen_volumes: the fluid volumes of the lumen segments of the
            segmented intestine by their token (`ids.GUT_SEGMENTS`) [l],
            `Vlumen_<seg>`, the same for both sexes (the small intestine
            `du` to `ile2` 0.105 l).
    """

    sex: Sex
    bw: float
    height: float
    gfr: float
    volumes: Mapping[Organ, float]
    tissue_volumes: Mapping[Organ, float]
    intracellular_volumes: Mapping[Organ, float]
    lumen_volumes: Mapping[str, float]


@cache
def reference_individual(sex: Sex, age: float = AGE_REFERENCE) -> ReferenceIndividual:
    """The reference individual of a sex and an age.

    Args:
        sex: the sex of the reference individual.
        age: the age (years), the adult of ICRP 89 by default; a child is the
            interpolation of the reference children of ICRP 89, an older
            adult the adult with the factors of the age (`pkpdsim.pbpk.age`).

    Returns:
        Its body weight, height, healthy GFR, organ, cell and lumen volumes.
    """
    from pkpdsim.pbpk.data import value
    from pkpdsim.pbpk.variants.segments import DATA, SEGMENTS

    refs = _references(sex, load("physiology_human"))
    at = _at_age(sex, age)
    tissue = reference_tissue_volumes(sex, age)
    vblood = at[ids.ref(ids.VBLOOD)]
    volumes = {o: tissue[o] + refs[f"blood_{o}"].value * vblood for o in Organ}
    bw, height = at[ids.ref(ids.BW)], at[ids.ref(ids.HEIGHT)]
    # GFR_base = GFR_ref*BSA/BSA_ref(adult)*fage_GFR, at the body surface area of the age
    gfr = (
        refs[ids.ref(ids.GFR)].value
        * _dubois(bw, height)
        / refs[ids.ref(ids.BSA)].value
        * factor_curves(sex)[ids.GFR].curve.value(age)
    )
    return ReferenceIndividual(
        sex=sex,
        bw=bw,
        height=height,
        gfr=gfr,
        volumes=MappingProxyType(volumes),
        tissue_volumes=MappingProxyType(tissue),
        intracellular_volumes=MappingProxyType(
            {o: tissue[o] * (1.0 - refs[f"fint_{o}"].value) for o in Organ}
        ),
        lumen_volumes=MappingProxyType(
            {t: value(DATA, key, "volume").value for t, key in SEGMENTS.items()}
        ),
    )
