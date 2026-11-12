"""The pregnancy: the maternal physiology over the gestational age and the fetus.

A model with a pregnancy (`PBPKDefinition.pregnancy`, `Pregnancy`) knows the
gestational age: `GA` (weeks since the first day of the last menstrual period,
0 not pregnant, 40 term) is an initial value which an experiment sets, as the
covariates. Every change of the maternal physiology is a factor
`fpreg_<key>` of `database/tables/pregnancy_human.toml`, the published fits
relative to the non-pregnant woman (1 at GA 0), a piecewise linear function of
`GA` which multiplies its quantity as the factors of the impairment do: the
cardiac output and with it every organ flow, the blood volume, the hematocrit,
the GFR, the plasma albumin and alpha-1-acid glycoprotein, the small
intestinal transit and the activity of the enzymes and transporters of the
liver and the kidney.

The body weight `BW` stays the body weight before the pregnancy, from which
the organs of the woman are scaled; the masses of the pregnancy (fetus,
placenta, amniotic fluid, the gain of the uterus and of the breasts) are added
to `BW_total`, the weight of the pregnant woman.

The fetoplacental unit (`Pregnancy(fetus=True)`, the default) is the placenta
and the fetus as two compartments of the top model:

- the placenta `Vpla` is perfused with `Qpla`, the gain of the uterine
  arterial blood flow (`fpreg_upla` of the cardiac output), from the arterial
  blood into the venous blood, and holds the substance at `Kp_pla_<s>`;
- the fetus `Vfet` exchanges with the placenta with the transfer clearance
  `PS_pla_<s>` on the unbound concentrations and holds the substance at
  `Kp_fet_<s>`; `Cfet_<s>` is the fetal concentration, the model of the cord
  blood, and the fetal to maternal ratio follows from it.

The fetus neither metabolises nor excretes the substance in this version.
"""

import bisect
from collections.abc import Mapping
from dataclasses import dataclass
from functools import cache
from typing import TYPE_CHECKING, Any

from sbmlutils.factory import (
    Compartment,
    Parameter,
    Reaction,
    Sbase,
    Species,
)
from sbmlutils.metadata import SBO

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.annotation import annotate, initial_assignment
from pkpdsim.pbpk.data import cite, load
from pkpdsim.pbpk.organs import Organ
from pkpdsim.pbpk.units import WEEK, U
from pkpdsim.sex import Sex
from pkpdsim.species import Species as BuildSpecies

if TYPE_CHECKING:
    from pkpdsim.pbpk.definition import PBPKDefinition

__all__ = [
    "GESTATION",
    "ORGAN_FLOW",
    "PLACENTA_SHARE",
    "TERM",
    "Factor",
    "Pregnancy",
    "factors",
    "pregnancy_factor",
    "pregnancy_objects",
    "pregnancy_physiology_factors",
    "protein_factors",
]

GESTATION = "pregnancy_human"
"""The table of the pregnancy in the database."""

TERM = 40.0
"""The gestational age (weeks) of a birth at term."""

ORGAN_FLOW = "Qorgans"
"""The key of the factor of the pregnancy on every organ flow: the cardiac
output less the share of the placenta."""

PLACENTA_SHARE = "Qpla"
"""The key of the blood flow of the placenta as fraction of the cardiac output
of the non-pregnant woman."""

_V_MIN = 1e-9
"""The volume (l) of the placenta and the fetus before the pregnancy, so that
their concentrations are defined at GA 0 (they are not perfused there)."""

_SBO_PARAMETER = SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER

_KEYS: Mapping[str, str] = {
    "cardiac_output": ids.CO,
    "placental_fraction": "upla",
    "blood_volume": ids.VBLOOD,
    "hematocrit": ids.HCT,
    "gfr": ids.GFR,
    "kidney_flow": ids.organ_flow(Organ.KIDNEY),
    "albumin": "albumin",
    "agp": "agp",
    "t_si": "t_si",
}
"""The key of the factor of every entry of the table."""


@dataclass(frozen=True)
class Factor:
    """A quantity of the pregnancy over the gestational age.

    Attributes:
        key: the key of its id (`CO`, `li_cyp3a4`, `mass_fetus`).
        ga: the gestational ages (weeks) of the values, increasing.
        values: the factor (or the mass in kg) at the gestational ages,
            linear between them and constant beyond them.
        unit: the unit of the values.
        sources: the keys of the sources.
        notes: where the values come from.
    """

    key: str
    ga: tuple[float, ...]
    values: tuple[float, ...]
    unit: str
    sources: tuple[str, ...]
    notes: str

    def value(self, ga: float) -> float:
        """The factor at a gestational age (weeks)."""
        if ga <= self.ga[0]:
            return self.values[0]
        if ga >= self.ga[-1]:
            return self.values[-1]
        i = bisect.bisect_right(self.ga, ga) - 1
        a0, a1 = self.ga[i], self.ga[i + 1]
        v0, v1 = self.values[i], self.values[i + 1]
        return v0 + (v1 - v0) * (ga - a0) / (a1 - a0)

    def formula(self) -> str:
        """The dimensionless piecewise formula of the factor over `GA`."""

        def n(x: float) -> str:
            return f"{float(x)!r} dimensionless"

        week = f"({ids.GA} / 1 week)"
        terms = [n(self.values[0]), f"{week} <= {n(self.ga[0])}"]
        for i in range(len(self.ga) - 1):
            a0, a1 = self.ga[i], self.ga[i + 1]
            v0, v1 = self.values[i], self.values[i + 1]
            slope = (v1 - v0) / (a1 - a0)
            terms += [
                f"{n(v0)} + {n(slope)} * ({week} - {n(a0)})",
                f"{week} < {n(a1)}",
            ]
        terms.append(n(self.values[-1]))
        return f"piecewise({', '.join(terms)})"


def pregnancy_factor(key: str) -> str:
    """The id of a factor of the pregnancy, `fpreg_CO`, `fpreg_li_cyp3a4`."""
    return f"fpreg_{key}"


def _factor(key: str, entry: Mapping[str, Any], ga: tuple[float, ...]) -> Factor:
    """The factor of an entry of the table."""
    return Factor(
        key,
        ga,
        tuple(float(v) for v in entry["value"]),
        str(entry["unit"]),
        (entry["source"], *entry.get("sources", ())),
        str(entry["notes"]),
    )


@cache
def _table() -> dict[str, Any]:
    """The table of the pregnancy."""
    return load(GESTATION)


@cache
def factors() -> dict[str, Factor]:
    """The factors of the maternal physiology and the masses, by the key of their id."""
    data = _table()
    ga = tuple(float(g) for g in data["ga"])
    found = {key: _factor(key, data[name], ga) for name, key in _KEYS.items()}
    found |= {
        f"mass_{name}": _factor(f"mass_{name}", entry, ga)
        for name, entry in data["mass"].items()
    }
    return found


@cache
def protein_factors() -> dict[tuple[Organ, str], Factor]:
    """The factors of the enzymes and transporters, by organ and protein."""
    data = _table()
    ga = tuple(float(g) for g in data["ga"])
    return {
        (Organ(organ), protein): _factor(f"{organ}_{protein.lower()}", entry, ga)
        for organ, proteins in data["protein"].items()
        for protein, entry in proteins.items()
    }


@dataclass(frozen=True)
class Pregnancy:
    """The pregnancy of a model.

    Attributes:
        ga: the default gestational age (weeks, 0 to `TERM`), the value of
            `GA`; 0 (not pregnant) by default, so the model of a pregnancy
            which an experiment does not set is the non-pregnant woman.
        fetus: build the fetoplacental unit (the placenta and the fetus as
            compartments of the top model); `False` gives the maternal
            physiology alone.

    Raises:
        ValueError: for a gestational age outside of `[0, TERM]`.
    """

    ga: float = 0.0
    fetus: bool = True

    def __post_init__(self) -> None:
        """Check the gestational age."""
        if not 0.0 <= self.ga <= TERM:
            raise ValueError(
                f"Pregnancy: the gestational age {self.ga} is not in [0, {TERM:g}] weeks."
            )


def _check(definition: PBPKDefinition) -> None:
    """Check that a definition may be pregnant.

    Raises:
        ValueError: for a male build, a laboratory animal, or a substance
            whose erythrocytes are kinetic (the placenta and the fetus carry
            the plasma alone).
    """
    problems = []
    if definition.sex is not Sex.FEMALE:
        problems.append(
            f"the build is {definition.sex.value}: a pregnancy needs the female build"
        )
    if definition.species is not BuildSpecies.HUMAN:
        problems.append(
            f"the species is {definition.species.value}: the pregnancy of the "
            "library is the human one"
        )
    kinetic = [
        s.sid
        for s in definition.substances
        if s.erythrocytes is not None and s.erythrocytes.mode == "kinetic"
    ]
    if kinetic:
        problems.append(
            f"the substances {kinetic} have kinetic erythrocytes, which the "
            "placenta and the fetus do not carry"
        )
    if problems:
        raise ValueError(
            f"Invalid pregnancy of '{definition.sid}': " + "; ".join(problems) + "."
        )


def pregnancy_physiology_factors(definition: PBPKDefinition) -> dict[str, str]:
    """The factors of the pregnancy on the quantities of the physiology.

    Args:
        definition: the definition.

    Returns:
        The factor by the id of the quantity it multiplies (the blood
        volume); empty without a pregnancy. The other quantities (the
        hematocrit, the GFR, the flows) carry their factor through the
        product of the impairment (`pkpdsim.pbpk.impairment`).
    """
    if definition.pregnancy is None:
        return {}
    return {ids.VBLOOD: pregnancy_factor(ids.VBLOOD)}


def _parameter(
    pid: str,
    value: float | None,
    unit: Any,  # noqa: ANN401
    name: str,
    notes: str,
    kind: str,
) -> Parameter:
    """A constant parameter of the top model with its annotation."""
    return annotate(
        Parameter(
            pid,
            value,
            unit,
            constant=True,
            name=name,
            sboTerm=_SBO_PARAMETER,
            notes=notes,
        ),
        kind,
    )


def _factor_objects(f: Factor) -> list[Sbase]:
    """The parameter and the initial assignment of a factor of the pregnancy."""
    pid = pregnancy_factor(f.key)
    unit = U.kg if f.unit == "kg" else U.dimensionless
    citations = "\n".join(f"- {cite(s)}" for s in f.sources)
    return [
        _parameter(
            pid,
            None,
            unit,
            f"pregnancy {f.key.replace('_', ' ')}",
            f"The {'mass' if f.unit == 'kg' else 'factor'} of the pregnancy over "
            f"the gestational age `{ids.GA}`, linear between the gestational ages "
            f"of the table: {f.notes}\n\n{citations}",
            "pregnancy",
        ),
        initial_assignment(
            pid,
            f.formula() if f.unit != "kg" else f"({f.formula()}) * 1 kg",
            unit,
            f"pregnancy {f.key}",
        ),
    ]


def _unit_objects() -> list[Sbase]:
    """The compartments and parameters of the fetoplacental unit."""
    placenta = pregnancy_factor("mass_placenta")
    fetus = pregnancy_factor("mass_fetus")
    return [
        annotate(
            Compartment(
                ids.V_PLACENTA,
                None,
                U.litre,
                constant=True,
                name="volume placenta [l]",
                notes="The volume of the placenta of the gestational age (the mass "
                f"of 1 kg/l), at least {_V_MIN:g} l so that its concentrations are "
                "defined before the pregnancy, where it has no blood flow.",
            ),
            "placenta",
        ),
        initial_assignment(
            ids.V_PLACENTA,
            f"piecewise({_V_MIN!r} litre, {placenta} * 1 l_per_kg < {_V_MIN!r} litre, "
            f"{placenta} * 1 l_per_kg)",
            U.litre,
            "volume placenta",
        ),
        annotate(
            Compartment(
                ids.V_FETUS,
                None,
                U.litre,
                constant=True,
                name="volume fetus [l]",
                notes="The volume of the fetus of the gestational age (the mass of "
                f"1 kg/l), at least {_V_MIN:g} l so that its concentrations are "
                "defined before the pregnancy.",
            ),
            "fetus",
        ),
        initial_assignment(
            ids.V_FETUS,
            f"piecewise({_V_MIN!r} litre, {fetus} * 1 l_per_kg < {_V_MIN!r} litre, "
            f"{fetus} * 1 l_per_kg)",
            U.litre,
            "volume fetus",
        ),
        _parameter(
            ids.Q_PLACENTA,
            None,
            U.l_per_min,
            "blood flow placenta [l/min]",
            "The blood flow of the placenta: the gain of the bilateral uterine "
            "arterial flow of the pregnancy, the fraction `fpreg_upla` of the "
            "cardiac output; 0 before the pregnancy.",
            "placenta",
        ),
        initial_assignment(
            ids.Q_PLACENTA,
            f"{ids.CO_BASE} * {pregnancy_factor(PLACENTA_SHARE)}",
            U.l_per_min,
            "blood flow placenta",
        ),
    ]


def _substance_objects(definition: PBPKDefinition) -> list[Sbase]:
    """The species of the placenta and the fetus and their exchange."""
    objects: list[Sbase] = []
    plasma_flow = f"({ids.Q_PLACENTA} * (1 dimensionless - {ids.HCT}))"
    for s in definition.substances:
        sid = s.sid
        pla, fet = ids.placenta_species(sid), ids.fetus_species(sid)
        kp_pla, kp_fet = ids.kp_placenta(sid), ids.kp_fetus(sid)
        ps = ids.ps_placenta(sid)
        objects += [
            _parameter(
                kp_pla,
                1.0,
                U.dimensionless,
                f"placenta to plasma coefficient {sid}",
                "The concentration in the placenta over the plasma concentration "
                "leaving it; 1 by default (the placenta distributes like plasma), "
                "a model fits it or sets it from a perfusion study.",
                "kp",
            ),
            _parameter(
                kp_fet,
                1.0,
                U.dimensionless,
                f"fetus to fetal plasma coefficient {sid}",
                "The concentration in the fetus over its plasma concentration; 1 "
                "by default, a model sets it from the Kp of the fetal tissues.",
                "kp",
            ),
            _parameter(
                ps,
                None,
                U.l_per_min,
                f"transfer clearance placenta {sid} [l/min]",
                "The transfer clearance of the placenta on the unbound "
                "concentrations; the plasma flow of the placenta by default, so "
                "the transfer is limited by the flow (a choice of the library), a "
                "model fits it or sets it from a placenta perfusion study.",
                "kp",
            ),
            initial_assignment(
                ps, plasma_flow, U.l_per_min, f"transfer clearance placenta {sid}"
            ),
            annotate(
                Species(
                    pla,
                    ids.V_PLACENTA,
                    initialConcentration=0.0,
                    substanceUnit=U.mmole,
                    hasOnlySubstanceUnits=False,
                    name=f"{s.name} (placenta)",
                ),
                "concentration",
                substance=s,
            ),
            annotate(
                Species(
                    fet,
                    ids.V_FETUS,
                    initialConcentration=0.0,
                    substanceUnit=U.mmole,
                    hasOnlySubstanceUnits=False,
                    name=f"{s.name} (fetus)",
                ),
                "concentration",
                substance=s,
            ),
            annotate(
                Reaction(
                    sid=f"J_pla_in_{sid}",
                    equation=f"{ids.pool_species(ids.Pool.ARTERIAL, sid)} -> {pla}",
                    formula=(
                        f"{plasma_flow} * {ids.pool_species(ids.Pool.ARTERIAL, sid)}",
                        U.mmole_per_min,
                    ),
                    name=f"inflow placenta {sid}",
                ),
                "exchange",
            ),
            annotate(
                Reaction(
                    sid=f"J_pla_out_{sid}",
                    equation=f"{pla} -> {ids.pool_species(ids.Pool.VENOUS, sid)}",
                    formula=(
                        f"{plasma_flow} * {pla} / {kp_pla}",
                        U.mmole_per_min,
                    ),
                    name=f"outflow placenta {sid}",
                ),
                "exchange",
            ),
            annotate(
                Reaction(
                    sid=f"J_fet_{sid}",
                    equation=f"{pla} <-> {fet}",
                    formula=(
                        f"{ps} * ({ids.fu(sid)} * {pla} / {kp_pla} - "
                        f"{ids.fu(sid)} * {fet} / {kp_fet})",
                        U.mmole_per_min,
                    ),
                    name=f"placental transfer {sid}",
                ),
                "exchange",
            ),
        ]
    return objects


def pregnancy_objects(definition: PBPKDefinition) -> list[Sbase]:
    """The gestational age, the factors, the masses and the fetoplacental unit.

    Args:
        definition: the definition with `pregnancy`.

    Returns:
        `GA`, every factor `fpreg_<key>` of the physiology and of the proteins
        the definition uses, the masses of the pregnancy with `BW_total`, and
        the placenta and the fetus with their species and exchanges; empty
        without a pregnancy.

    Raises:
        ValueError: for a male build, a laboratory animal or kinetic
            erythrocytes (`_check`).
    """
    from pkpdsim.pbpk.chemistry import Metabolism, Transport

    p = definition.pregnancy
    if p is None:
        return []
    _check(definition)
    objects: list[Sbase] = [
        _parameter(
            ids.GA,
            p.ga,
            WEEK,
            "gestational age [week]",
            "The gestational age, the weeks since the first day of the last "
            "menstrual period: 0 not pregnant, 40 at term. An initial value, as "
            "the covariates, which an experiment sets (`pkpdsim.pbpk.pregnancy`).",
            "gestational_age",
        )
    ]
    found = factors()
    for f in found.values():
        objects += _factor_objects(f)
    used = {
        (x.organ, x.protein)
        for x in definition.chemistry
        if isinstance(x, Metabolism | Transport) and x.protein is not None
    }
    for pair, f in protein_factors().items():
        if pair in used:
            objects += _factor_objects(f)
    # the flow of the organs and the share of the placenta
    co, upla = pregnancy_factor(ids.CO), pregnancy_factor("upla")
    objects += [
        _parameter(
            pregnancy_factor(ORGAN_FLOW),
            None,
            U.dimensionless,
            "pregnancy factor organ flows",
            "The factor of the pregnancy on the blood flow of every organ: the "
            "cardiac output less the share of the placenta, which the organs do "
            f"not receive, `{co} * (1 - {upla})`.",
            "pregnancy",
        ),
        initial_assignment(
            pregnancy_factor(ORGAN_FLOW),
            f"{co} * (1 dimensionless - {upla})",
            U.dimensionless,
            "pregnancy factor organ flows",
        ),
        _parameter(
            pregnancy_factor(PLACENTA_SHARE),
            None,
            U.dimensionless,
            "pregnancy blood flow placenta",
            "The blood flow of the placenta as fraction of the cardiac output of "
            f"the non-pregnant woman, `{co} * {upla}`; it enters the cardiac "
            "output next to the flows of the organs.",
            "pregnancy",
        ),
        initial_assignment(
            pregnancy_factor(PLACENTA_SHARE),
            f"{co} * {upla}",
            U.dimensionless,
            "pregnancy blood flow placenta",
        ),
        _parameter(
            ids.BW_TOTAL,
            None,
            U.kg,
            "total body weight [kg]",
            "The body weight of the pregnant woman: the body weight `BW` of the "
            "covariates (before the pregnancy, from which the organs are scaled) "
            "and the masses of the pregnancy (fetus, placenta, amniotic fluid and "
            "the gain of the uterus and of the breasts).",
            "body_weight",
        ),
        initial_assignment(
            ids.BW_TOTAL,
            " + ".join(
                [ids.BW, *(pregnancy_factor(k) for k in found if k.startswith("mass_"))]
            ),
            U.kg,
            "total body weight",
        ),
    ]
    if p.fetus:
        objects += _unit_objects()
        objects += _substance_objects(definition)
    return objects
