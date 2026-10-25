"""The standard human physiology and its scaling to an individual.

The physiology is a set of SBML parameters of the top model. The covariates
of the individual (`BW`, `HEIGHT`, `AGE`) and the factors of organ function
are constant parameters which an experiment changes before the simulation,
as it changes `BW` today; everything else follows from them by initial
assignments on constant parameters, which sbmlsim evaluates again after the
changes of an experiment (a value which no rule recomputes during the
integration keeps the compiled model small and fast to load):

- the reference values of the sex of the build (ICRP 89,
  `data/physiology_human.toml`); the sex is chosen when the model is built,
  a study of both sexes uses one model per sex;
- body mass index, body surface area (Du Bois 1916) and fat-free mass
  (Janmahasatian 2005);
- the tissues scale as in Huisinga 2012 (eqs. 9, 10, 12): brain is constant,
  skin scales with the body surface area, the other lean organs and the blood
  with the fat-free mass less brain and skin (`SV_tis`); adipose is the
  remainder of the body weight;
- the cardiac output scales as the tissues, `SV_tis` (Huisinga 2012 eq. 13),
  and the organ flows are fractions of it, renormalised over the organs so
  that they sum to the cardiac output;
- the GFR scales with the body surface area.

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
follows.
"""

from collections.abc import Iterable
from dataclasses import dataclass
from enum import StrEnum
from typing import Any

from sbmlutils.factory import InitialAssignment, Parameter, Sbase, UnitDefinition
from sbmlutils.metadata import SBO

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.data import cite, load
from pkpdsim.pbpk.organs import PORTAL_ORGANS, Organ, Pool
from pkpdsim.pbpk.units import U

__all__ = [
    "Sex",
    "physiology_objects",
    "reference_fractions",
    "reference_tissue_volumes",
]


class Sex(StrEnum):
    """The sex of the reference individual."""

    MALE = "male"
    FEMALE = "female"


LEAN_ORGANS: tuple[Organ, ...] = tuple(o for o in Organ if o is not Organ.ADIPOSE)
"""The organs whose tissue mass scales with the individual; adipose is the rest."""

AGE_DEFAULT = 40.0
"""The age of an individual if an experiment does not set it (yr)."""

MODEL_DEFAULT = "Model default without a source."
"""The notes of a parameter whose value is a choice of the model, not data."""

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
) -> Parameter:
    """A constant parameter, of a value or of its initial assignment."""
    return Parameter(
        pid,
        value,
        unit,
        constant=True,
        name=name,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        notes=_notes(text, sources),
    )


def _initial(
    pid: str,
    formula: str,
    unit: UnitDefinition,
    name: str,
    text: str,
    sources: Iterable[str],
) -> list[Sbase]:
    """A constant parameter with its initial assignment."""
    return [
        _parameter(pid, unit, name, f"{text}\n\n`{pid} = {formula}`", sources),
        InitialAssignment(pid, formula, unit),
    ]


class _Physiology:
    """The SBML objects of the physiology, from the reference values of one sex."""

    def __init__(self, sex: Sex) -> None:
        self.sex = sex
        self.refs = _references(sex, load("physiology_human"))

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
        self, pid: str, unit: UnitDefinition, name: str, text: str
    ) -> list[Sbase]:
        """A constant reference parameter of the sex with its provenance."""
        notes, sources = self.provenance(pid)
        return [
            _parameter(
                pid,
                unit,
                name,
                f"{text}\n\n{notes}",
                sources,
                value=self.refs[pid].value,
            )
        ]

    def covariates(self) -> list[Sbase]:
        """The covariates of the individual and the factors of organ function."""
        default = self.refs
        return [
            _parameter(
                ids.BW,
                U.kg,
                "body weight [kg]",
                "Covariate, an initial value; the default is the reference "
                f"individual of the sex. {default[ids.BW].notes}",
                default[ids.BW].sources,
                value=default[ids.BW].value,
            ),
            _parameter(
                ids.HEIGHT,
                U.cm,
                "height [cm]",
                "Covariate, an initial value; the default is the reference "
                f"individual of the sex. {default[ids.HEIGHT].notes}",
                default[ids.HEIGHT].sources,
                value=default[ids.HEIGHT].value,
            ),
            _parameter(
                ids.AGE,
                U.yr,
                "age [yr]",
                f"Covariate. {MODEL_DEFAULT} No rule uses AGE yet.",
                (),
                value=AGE_DEFAULT,
            ),
            _parameter(
                ids.F_RENAL_FUNCTION,
                U.dimensionless,
                "renal function",
                "The GFR relative to the normal GFR of the individual; 1.0 "
                "normal, < 1.0 reduced. An initial value, as the covariates. "
                f"{MODEL_DEFAULT}",
                (),
                value=1.0,
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
            ),
        ]

    def references(self) -> list[Sbase]:
        """The reference values of the sex."""
        objects = self.reference(
            ids.ref(ids.BSA),
            U.m2,
            "reference body surface area [m2]",
            "Body surface area of the reference individual, the reference of "
            "skin and GFR.",
        )
        objects += self.reference(
            ids.ref(ids.FFM),
            U.kg,
            "reference fat-free mass [kg]",
            "Fat-free mass of the reference individual.",
        )
        objects += self.reference(
            ids.ref(ids.CO),
            U.l_per_min,
            "reference cardiac output [l/min]",
            "Cardiac output of the reference individual.",
        )
        objects += self.reference(
            ids.ref(ids.GFR),
            U.l_per_min,
            "reference glomerular filtration rate [l/min]",
            "GFR of the reference individual.",
        )
        objects += self.reference(
            ids.ref(ids.VBLOOD),
            U.litre,
            "reference blood volume [l]",
            "Total blood volume of the reference individual.",
        )
        objects += self.reference(
            ids.HCT,
            U.dimensionless,
            "hematocrit",
            "The volume fraction of the blood cells in blood.",
        )
        for o in LEAN_ORGANS:
            objects += self.reference(
                ids.ref(ids.organ_mass(o)),
                U.kg,
                f"reference tissue mass {_label(o)} [kg]",
                "Mass of the organ tissue without its blood of the reference "
                "individual.",
            )
        for o in Organ:
            if o is Organ.LUNG:
                continue
            objects += self.reference(
                ids.flow_fraction(o),
                U.dimensionless,
                f"blood flow fraction {_label(o)}",
                "Blood flow as fraction of cardiac output"
                + (" (hepatic artery)." if o is Organ.LIVER else "."),
            )
        for pool in Pool:
            objects += self.reference(
                ids.pool_fraction(pool),
                U.dimensionless,
                f"blood fraction {pool.name.lower()} pool",
                "Volume of the blood pool as fraction of the blood volume.",
            )
        return objects

    def body(self) -> list[Sbase]:
        """Body mass index, body surface area, fat-free mass and `SV_tis`."""
        objects = _initial(
            ids.BMI,
            f"({ids.BW} / 1 kg) / ({ids.HEIGHT} / 100 cm)^2 * 1 kg_per_m2",
            U.kg_per_m2,
            "body mass index [kg/m2]",
            "Body mass index of the individual.",
            ("Janmahasatian2005",),
        )
        objects += _initial(
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

        k, a = (6680, 216) if self.sex is Sex.MALE else (8780, 244)
        objects += _initial(
            ids.FFM,
            f"{ids.BW} * {ffm(k, a)}",
            U.kg,
            "fat-free mass [kg]",
            f"Fat-free mass of the individual, the {self.sex.value} formula of "
            "Janmahasatian 2005; Huisinga 2012 (eq. 6) use it for the lean body "
            "weight.",
            ("Janmahasatian2005", "Huisinga2012"),
        )
        br, sk = ids.organ_mass(Organ.BRAIN), ids.organ_mass(Organ.SKIN)
        objects += _initial(
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
        )
        return objects

    def volumes(self) -> list[Sbase]:
        """Tissue masses, organ volumes, blood volume and blood pools."""
        objects = _initial(
            ids.VBLOOD,
            f"{ids.ref(ids.VBLOOD)} * {ids.SV_TIS}",
            U.litre,
            "blood volume [l]",
            "Blood volume of the individual; it scales with the lean tissues "
            "(Huisinga 2012 eq. 12, blood is one of the remaining tissues).",
            ("Huisinga2012",),
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
            objects += _initial(
                mid,
                formula,
                U.kg,
                f"tissue mass {_label(o)} [kg]",
                f"Mass of the organ tissue without its blood; {text}",
                ("Huisinga2012",),
            )
        objects += _initial(
            ids.MLEAN,
            " + ".join(ids.organ_mass(o) for o in LEAN_ORGANS),
            U.kg,
            "mass of the lean organ tissues [kg]",
            "Mass of the tissues of the lean organs without their blood.",
            ("ICRP89",),
        )
        mad = ids.organ_mass(Organ.ADIPOSE)
        notes, sources = self.provenance("rho_blood")
        objects += _initial(
            mad,
            f"{ids.BW} - {ids.MLEAN} - {ids.VBLOOD} * "
            f"{self.number('rho_blood', U.kg_per_l)}",
            U.kg,
            "tissue mass adipose [kg]",
            "Mass of adipose tissue without its blood, the remainder of the body "
            "weight; it includes the tissues of ICRP 89 which are not organs of "
            f"the model. Density of blood:\n\n{notes}",
            sources,
        )
        for o in Organ:
            notes, sources = self.provenance(f"rho_{o}", f"blood_{o}")
            objects += _initial(
                ids.organ_volume(o),
                f"{ids.organ_mass(o)} / {self.number(f'rho_{o}', U.kg_per_l)} + "
                f"{self.number(f'blood_{o}', U.dimensionless)} * {ids.VBLOOD}",
                U.litre,
                f"volume {_label(o)} [l]",
                "Volume of the organ tissue and its blood, with the density of "
                "the tissue and the blood of the organ as fraction of the blood "
                f"volume:\n\n{notes}",
                sources,
            )
        for pool in Pool:
            objects += _initial(
                ids.pool_volume(pool),
                f"{ids.pool_fraction(pool)} * {ids.VBLOOD}",
                U.litre,
                f"volume {pool.name.lower()} blood [l]",
                "Volume of the blood pool.",
                ("ICRP89",),
            )
        return objects

    def fractions(self) -> list[Sbase]:
        """The vascular and interstitial fractions of the organ volumes."""
        objects: list[Sbase] = []
        for o in Organ:
            f_vas, f_int = ids.vascular_fraction(o), ids.interstitial_fraction(o)
            notes, sources = self.provenance(f"blood_{o}")
            objects += _initial(
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
            )
            notes, sources = self.provenance(f"fint_{o}")
            objects += _initial(
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
            )
        return objects

    def flows(self) -> list[Sbase]:
        """Cardiac output, blood flows and GFR."""
        objects = _initial(
            ids.CO,
            f"{ids.ref(ids.CO)} * {ids.SV_TIS} * {ids.F_CARDIAC_FUNCTION}",
            U.l_per_min,
            "cardiac output [l/min]",
            "Cardiac output of the individual; it scales as the lean tissues "
            "(Huisinga 2012 eq. 13, SQ_co = SV_tis) and with the cardiac function.",
            ("Huisinga2012",),
        )
        objects += _initial(
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
            objects += _initial(
                ids.organ_flow(o),
                f"{ids.flow_fraction(o)} * {ids.CO}",
                U.l_per_min,
                f"blood flow {_label(o)} [l/min]",
                "Blood flow of the organ"
                + (", the hepatic artery." if o is Organ.LIVER else "."),
                ("ICRP89",),
            )
        objects += _initial(
            ids.QHA,
            ids.organ_flow(Organ.LIVER),
            U.l_per_min,
            "hepatic arterial blood flow [l/min]",
            "Blood flow of the hepatic artery, the flow of the liver of the flow "
            "fractions.",
            ("ICRP89",),
        )
        objects += _initial(
            ids.QPO,
            " + ".join(ids.organ_flow(o) for o in Organ if o in PORTAL_ORGANS),
            U.l_per_min,
            "portal venous blood flow [l/min]",
            "Blood flow of the portal vein, the flows of the organs which drain "
            "into it.",
            ("ICRP89",),
        )
        objects += _initial(
            ids.QH,
            f"{ids.QHA} + {ids.QPO}",
            U.l_per_min,
            "liver blood flow [l/min]",
            "Total blood flow of the liver, hepatic artery and portal vein.",
            ("ICRP89",),
        )
        objects += _initial(
            ids.GFR,
            f"{ids.ref(ids.GFR)} * {ids.BSA} / {ids.ref(ids.BSA)} * "
            f"{ids.F_RENAL_FUNCTION}",
            U.l_per_min,
            "glomerular filtration rate [l/min]",
            "GFR of the individual: the reference GFR, which is roughly "
            "proportional to body size (ICRP 89 Section 8.3.1), scaled with the "
            "body surface area, i.e. the same GFR per 1.73 m2, times the renal "
            "function.",
            ("ICRP89",),
        )
        return objects


def _label(organ: Organ) -> str:
    """The organ in a name, `small intestine`."""
    return organ.name.lower().replace("_", " ")


def physiology_objects(sex: Sex = Sex.MALE) -> list[Sbase]:
    """The parameters and initial assignments of the standard human physiology.

    Args:
        sex: the sex of the model: its reference values and the default
            individual (`BW`, `HEIGHT`). A study of both sexes uses one model
            per sex.

    Returns:
        The constant parameters and their initial assignments of the top model.

    Raises:
        ValueError: if the reference individual leaves no mass for adipose, or
            a unit of the data is not the unit of the model.
    """
    physiology = _Physiology(sex)
    return [
        *physiology.covariates(),
        *physiology.references(),
        *physiology.body(),
        *physiology.volumes(),
        *physiology.fractions(),
        *physiology.flows(),
    ]


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


def reference_tissue_volumes(sex: Sex) -> dict[Organ, float]:
    """The tissue volumes of the organs of a reference individual.

    The tissue volume of an organ is its volume without its blood,
    `V<o>*(1 - f_vas_<o>) = M<o>/density`; adipose is the remainder of the
    body weight, as in the model.

    Args:
        sex: the sex of the reference individual.

    Returns:
        The tissue volume (l) by organ.
    """
    refs = _references(sex, load("physiology_human"))
    masses = {o: refs[ids.ref(ids.organ_mass(o))].value for o in LEAN_ORGANS}
    masses[Organ.ADIPOSE] = (
        refs[ids.BW].value
        - sum(masses.values())
        - refs[ids.ref(ids.VBLOOD)].value * refs["rho_blood"].value
    )
    return {o: masses[o] / refs[f"rho_{o}"].value for o in Organ}
