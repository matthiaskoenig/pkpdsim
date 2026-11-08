"""The observers of a machine perfusion: every quantity at every site.

The perfusate quantities are species of every site, `<q>_<site>` (the
reservoir `res`, the inflow `in`, the outflow `out`, the vascular space of
every zone `z<k>`); the cells of a zone are the site `z<k>_cell`. The
observers add what the species do not give, as assignment rules:

- `pO2_<site>` (kPa), the free oxygen over its solubility at the
  temperature, `o2_<site>/alpha_o2`;
- `<e>_act_<site>` (µkat/l), the catalytic activity of an enzyme of
  `ids.ENZYMES`, its concentration times its specific activity `a_<e>`
  (µkat/mmole, `[enzymes]` of the perfusion table), the unit of the data;
- the periportal and pericentral means, the sites `pp` and `pc` (`pp_cell`
  for the cells): every concentration of a zone averaged over the first and
  the last third of the length of the sinusoid, each zone weighted by its
  overlap with that third (the zones have the same volume), and the pO2,
  the pH and the activities of these means as at any site. With one zone
  both means are the zone; with two zones the periportal mean is the first
  zone and the pericentral the second; with three zones the first and the
  last zone; with four zones the periportal mean is `3/4` of the first and
  `1/4` of the second zone;
- `pCO2_<site>` (kPa), the dissolved CO2 over its solubility;
- the bile of the organ, the rate `bile_rate` (µl/min, the sum of the
  secretion of the zones) and the volume `bile` (µl, the sum of the bile
  volumes of the zones);
- the oxygen uptake of the organ `VO2` (mmole/min), the oxygen which the
  cells take up from the perfusate, the sum of the zones.
"""

from collections.abc import Sequence
from fractions import Fraction

from sbmlutils.factory import AssignmentRule, Parameter
from sbmlutils.factory.units import UnitDefinition
from sbmlutils.metadata import SBO

from pkpdsim.pbpk.placement import value_parameter
from pkpdsim.pbpk.units import U
from pkpdsim.perfusion import core, ids
from pkpdsim.perfusion.definition import table_value
from pkpdsim.perfusion.units import PU

__all__ = [
    "bile_objects",
    "enzyme_parameters",
    "mean_objects",
    "site_objects",
    "third_weights",
    "uptake_objects",
]

_RULE = SBO.MATHEMATICAL_EXPRESSION


def _observer(
    sid: str, formula: str, unit: UnitDefinition, name: str, sbo: str
) -> tuple[Parameter, AssignmentRule]:
    """Get the parameter of an observer and its assignment rule."""
    return (
        Parameter(sid, None, unit, constant=False, name=name, sboTerm=sbo),
        AssignmentRule(sid, formula, unit, name=name, sboTerm=_RULE),
    )


def _add(objects: core.CoreObjects, observer: tuple[Parameter, AssignmentRule]) -> None:
    """Add an observer to the objects."""
    objects.parameters.append(observer[0])
    objects.rules.append(observer[1])


def enzyme_parameters() -> list[Parameter]:
    """Get the specific activity of every enzyme, `a_<e>` (µkat/mmole)."""
    return [
        value_parameter(
            ids.specific_activity(e),
            table_value("enzymes", e.upper()),
            PU.ukat_per_mmole,
            f"specific activity of {core.quantity_name(e)}",
            SBO.SPECIFIC_ACTIVITY,
        )
        for e in ids.ENZYMES
    ]


def site_objects(site: str) -> core.CoreObjects:
    """Get the pO2, the pCO2 and the enzyme activities of the perfusate at a site."""
    out = core.CoreObjects()
    _add(
        out,
        _observer(
            ids.po2(site),
            f"{ids.species('o2', site)}/{ids.SOLUBILITY}",
            PU.kPa,
            f"oxygen partial pressure ({site})",
            SBO.PRESSURE,
        ),
    )
    _add(
        out,
        _observer(
            ids.pco2(site),
            f"{ids.species('co2', site)}/{ids.CO2_SOLUBILITY}",
            PU.kPa,
            f"CO2 partial pressure ({site})",
            SBO.PRESSURE,
        ),
    )
    for e in ids.ENZYMES:
        _add(
            out,
            _observer(
                ids.activity(e, site),
                f"{ids.species(e, site)} * {ids.specific_activity(e)}",
                PU.ukat_per_l,
                f"{core.quantity_name(e)} activity ({site})",
                SBO.ENZYME_ACTIVITY,
            ),
        )
    return out


def third_weights(n: int, periportal: bool) -> list[tuple[int, Fraction]]:
    """Get the weights of the zones in the mean of a third of the sinusoid.

    Zone `k` of `n` covers `[(k - 1)/n, k/n]` of the length of the sinusoid;
    its weight is its overlap with the first (periportal) or the last
    (pericentral) third over the length of the third, so the weights sum to 1.

    Args:
        n: the number of zones.
        periportal: the first third, else the last.

    Returns:
        The zones (from 1) with a weight above 0 and their weights.
    """
    third = Fraction(1, 3)
    out = []
    for k in range(1, n + 1):
        lower, upper = Fraction(k - 1, n), Fraction(k, n)
        if not periportal:
            lower, upper = 1 - upper, 1 - lower
        overlap = min(upper, third) - lower
        if overlap > 0:
            out.append((k, overlap / third))
    return out


def _mean(of_zones: Sequence[str], weights: list[tuple[int, Fraction]]) -> str:
    """Get the weighted mean of an id of the zones (periportal first)."""
    terms = [
        of_zones[k - 1] if w == 1 else f"{float(w)!r} dimensionless * {of_zones[k - 1]}"
        for k, w in weights
    ]
    return " + ".join(terms)


def mean_objects(zones: Sequence[str], perfusate: Sequence[str]) -> core.CoreObjects:
    """Get the periportal and pericentral means of the zones.

    Args:
        zones: the zones, periportal first.
        perfusate: the quantities of the perfusate (with the tracers).

    Returns:
        The means of the perfusate concentrations, the cell concentrations
        and the living fraction at the sites `pp` and `pc`, and their pO2,
        pH and enzyme activities.
    """
    out = core.CoreObjects()
    for site, periportal in ((ids.PERIPORTAL, True), (ids.PERICENTRAL, False)):
        weights = third_weights(len(zones), periportal)
        where = "periportal" if periportal else "pericentral"

        for q in perfusate:
            _add(
                out,
                _observer(
                    ids.species(q, site),
                    _mean([ids.species(q, z) for z in zones], weights),
                    U.mM,
                    f"{q} (perfusate, {where} mean)",
                    SBO.MOLAR_CONCENTRATION_OF_AN_ENTITY,
                ),
            )
        cell = ids.cell_site(site)
        for q in ids.CELL:
            _add(
                out,
                _observer(
                    ids.species(q, cell),
                    _mean([ids.species(q, ids.cell_site(z)) for z in zones], weights),
                    U.mM,
                    f"{core.quantity_name(q)} (cells, {where} mean)",
                    SBO.MOLAR_CONCENTRATION_OF_AN_ENTITY,
                ),
            )
        _add(
            out,
            _observer(
                ids.alive(site),
                _mean([ids.alive(z) for z in zones], weights),
                U.dimensionless,
                f"living fraction of the cells ({where} mean)",
                SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            ),
        )
        _add(
            out,
            _observer(
                ids.ph(site),
                f"3 dimensionless - log10({ids.species('h', site)}/(1 mM))",
                U.dimensionless,
                f"pH ({where} mean)",
                SBO.PH,
            ),
        )
        out += site_objects(site)
    return out


def bile_objects(rates: Sequence[str], zones: Sequence[str]) -> core.CoreObjects:
    """Get the bile rate and the bile volume of the organ.

    Args:
        rates: the formulas of the bile secretion of the zones (l/min), the
            rate rules of `Vbile_<zone>` of the core.
        zones: the zones.

    Returns:
        `bile_rate` (µl/min) and `bile` (µl), the sum of the bile volumes of
        the zones, both assignment rules.
    """
    out = core.CoreObjects()
    ul = f"1000000 {PU.ul_per_l.sid}"
    total = " + ".join(f"({r})" for r in rates)
    _add(
        out,
        _observer(
            ids.BILE_RATE,
            f"{ul} * ({total})",
            PU.ul_per_min,
            "bile secretion of the organ",
            SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
    )
    _add(
        out,
        _observer(
            ids.BILE,
            f"{ul} * ({' + '.join(ids.bile(z) for z in zones)})",
            PU.ul,
            "bile volume of the organ",
            SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
    )
    return out


def uptake_objects(zones: Sequence[str]) -> core.CoreObjects:
    """Get the oxygen uptake of the organ, `VO2` (mmole/min)."""
    out = core.CoreObjects()
    _add(
        out,
        _observer(
            ids.VO2,
            " + ".join(ids.reaction("O2IM", z) for z in zones),
            U.mmole_per_min,
            "oxygen uptake of the organ",
            SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
    )
    return out
