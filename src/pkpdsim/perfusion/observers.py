"""The observers of a machine perfusion: every quantity at every site.

The perfusate quantities are species of every site, `<q>_<site>` (the
reservoir `res`, the inflow `in`, the outflow `out`, the vascular space of
every zone `z<k>`); the cells of a zone are the site `z<k>_cell`. The
observers add what the species do not give, as assignment rules:

- `pO2_<site>` (kPa), the free oxygen over its solubility at the
  temperature, `o2_<site>/alpha_o2`;
- `<e>_act_<site>` (µkat/l), the catalytic activity of an enzyme of
  `ids.ENZYMES` which the organ releases, its concentration times its
  specific activity `a_<e>` (µkat/mmole, `[enzymes]` of the perfusion
  table), the unit of the data;
- the periportal and pericentral means, the sites `pp` and `pc` (`pp_cell`
  for the cells): every concentration of a zone averaged over the first and
  the last third of the length of the sinusoid, each zone weighted by its
  overlap with that third (`third_weights` and `weighted_mean` of
  `pkpdsim.zones`, which every zonated model of pkpdsim shares; the zones
  have the same volume), and the pO2,
  the pH and the activities of these means as at any site. With one zone
  both means are the zone; with two zones the periportal mean is the first
  zone and the pericentral the second; with three zones the first and the
  last zone; with four zones the periportal mean is `3/4` of the first and
  `1/4` of the second zone;
- `pCO2_<site>` (kPa), the dissolved CO2 over its solubility;
- the bile of an organ with bile, the rate `bile_rate` (µl/min, the sum of
  the secretion of the zones), the volume `bile` (µl, the sum of the bile
  volumes of the zones) and the volume per gram of the organ `bile_per_mass`
  (µl/g, `bile/(1000 g/kg * M_<o>)`, the unit of the data);
- the urine of an organ with urine, the rate `urine_rate` (µl/min, the net
  urine flow `Q_urine`) and the volume `urine` (µl, the urine volume
  `V_urine`);
- the oxygen uptake of the organ `VO2` (mmole/min), the oxygen which the
  cells take up from the perfusate, the sum of the zones.

The parameter of an observer with its assignment rule (`observer`, `RULE`) is
the SBML writing kit's (`pkpdsim.sbml.parameters`), which the zonation engine
writes its observers with as well; this module re-exports them unchanged.
"""

from collections.abc import Sequence

from sbmlutils.factory import AssignmentRule, Parameter
from sbmlutils.metadata import SBO

from pkpdsim.pbpk.placement import value_parameter
from pkpdsim.pbpk.units import U
from pkpdsim.perfusion import circuit, core, ids
from pkpdsim.perfusion.definition import table_value
from pkpdsim.perfusion.units import PU
from pkpdsim.sbml.parameters import RULE, observer
from pkpdsim.zones import third_weights, weighted_mean

__all__ = [
    "RULE",
    "bile_objects",
    "enzyme_parameters",
    "mean_objects",
    "observer",
    "site_objects",
    "third_weights",
    "uptake_objects",
    "urine_objects",
    "weighted_mean",
]


def _add(objects: core.CoreObjects, observer: tuple[Parameter, AssignmentRule]) -> None:
    """Add an observer to the objects."""
    objects.parameters.append(observer[0])
    objects.rules.append(observer[1])


def enzyme_parameters(enzymes: Sequence[str] = ids.ENZYMES) -> list[Parameter]:
    """Get the specific activity of every enzyme, `a_<e>` (µkat/mmole).

    Args:
        enzymes: the enzymes of the perfusate, those the organ releases.
    """
    return [
        value_parameter(
            ids.specific_activity(e),
            table_value("enzymes", e.upper()),
            PU.ukat_per_mmole,
            f"specific activity of {core.quantity_name(e)}",
            SBO.SPECIFIC_ACTIVITY,
        )
        for e in enzymes
    ]


def site_objects(site: str, enzymes: Sequence[str] = ids.ENZYMES) -> core.CoreObjects:
    """Get the pO2, the pCO2 and the enzyme activities of the perfusate at a site.

    Args:
        site: the site.
        enzymes: the enzymes of the perfusate, those the organ releases.
    """
    out = core.CoreObjects()
    _add(
        out,
        observer(
            ids.po2(site),
            f"{ids.species('o2', site)}/{ids.SOLUBILITY}",
            PU.kPa,
            f"oxygen partial pressure ({site})",
            SBO.PRESSURE,
        ),
    )
    _add(
        out,
        observer(
            ids.pco2(site),
            f"{ids.species('co2', site)}/{ids.CO2_SOLUBILITY}",
            PU.kPa,
            f"CO2 partial pressure ({site})",
            SBO.PRESSURE,
        ),
    )
    for e in enzymes:
        _add(
            out,
            observer(
                ids.activity(e, site),
                f"{ids.species(e, site)} * {ids.specific_activity(e)}",
                PU.ukat_per_l,
                f"{core.quantity_name(e)} activity ({site})",
                SBO.ENZYME_ACTIVITY,
            ),
        )
    return out


def mean_objects(
    zones: Sequence[str],
    perfusate: Sequence[str],
    cells: Sequence[str] = ids.CELL,
    enzymes: Sequence[str] = ids.ENZYMES,
    *,
    submodel: bool = False,
) -> core.CoreObjects:
    """Get the periportal and pericentral means of the zones.

    The means of the first and the last third of the zones of any organ;
    the periportal and the pericentral mean in the liver.

    Args:
        zones: the zones, periportal first.
        perfusate: the quantities of the perfusate (with the tracers).
        cells: the quantities of the cells.
        enzymes: the enzymes of the perfusate.
        submodel: the cells of the zones are submodels (the liver).

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
                observer(
                    ids.species(q, site),
                    weighted_mean([ids.species(q, z) for z in zones], weights),
                    U.mM,
                    f"{circuit.quantity_label(q)} (perfusate, {where} mean)",
                    SBO.MOLAR_CONCENTRATION_OF_AN_ENTITY,
                ),
            )
        cell = ids.cell_site(site)
        for q in cells:
            _add(
                out,
                observer(
                    ids.species(q, cell),
                    weighted_mean(
                        [
                            ids.species(q, ids.cell_site(z), submodel=submodel)
                            for z in zones
                        ],
                        weights,
                    ),
                    U.mM,
                    f"{core.quantity_name(q)} (cells, {where} mean)",
                    SBO.MOLAR_CONCENTRATION_OF_AN_ENTITY,
                ),
            )
        _add(
            out,
            observer(
                ids.alive(site),
                weighted_mean(
                    [ids.alive(z, submodel=submodel) for z in zones], weights
                ),
                U.dimensionless,
                f"living fraction of the cells ({where} mean)",
                SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            ),
        )
        _add(
            out,
            observer(
                ids.ph(site),
                f"3 dimensionless - log10({ids.species('h', site)}/(1 mM))",
                U.dimensionless,
                f"pH ({where} mean)",
                SBO.PH,
            ),
        )
        out += site_objects(site, enzymes)
    return out


def bile_objects(
    rates: Sequence[str], zones: Sequence[str], mass: str
) -> core.CoreObjects:
    """Get the bile rate, the bile volume and the bile per mass of the organ.

    Args:
        rates: the formulas of the bile secretion of the zones (l/min), the
            rate rules of `Vbile_<zone>` of the core.
        zones: the zones.
        mass: the id of the mass of the organ (kg), `M_<o>`.

    Returns:
        `bile_rate` (µl/min), `bile` (µl), the sum of the bile volumes of
        the zones, and `bile_per_mass` (µl/g), the bile volume per gram of
        the organ, all assignment rules.
    """
    out = core.CoreObjects()
    ul = f"1000000 {PU.ul_per_l.sid}"
    total = " + ".join(f"({r})" for r in rates)
    _add(
        out,
        observer(
            ids.BILE_RATE,
            f"{ul} * ({total})",
            PU.ul_per_min,
            "bile secretion of the organ",
            SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
    )
    _add(
        out,
        observer(
            ids.BILE,
            f"{ul} * ({' + '.join(ids.bile(z, submodel=True) for z in zones)})",
            PU.ul,
            "bile volume of the organ",
            SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
    )
    _add(
        out,
        observer(
            ids.BILE_PER_MASS,
            f"{ids.BILE}/(1000 {PU.g_per_kg.sid} * {mass})",
            PU.ul_per_g,
            "bile volume per mass of the organ",
            SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
    )
    return out


def uptake_objects(zones: Sequence[str], *, submodel: bool = False) -> core.CoreObjects:
    """Get the oxygen uptake of the organ, `VO2` (mmole/min).

    Args:
        zones: the zones.
        submodel: the cells of the zones are submodels (the liver).
    """
    out = core.CoreObjects()
    _add(
        out,
        observer(
            ids.VO2,
            " + ".join(ids.reaction("O2IM", z, submodel=submodel) for z in zones),
            U.mmole_per_min,
            "oxygen uptake of the organ",
            SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
    )
    return out


def urine_objects() -> core.CoreObjects:
    """Get the urine rate and the urine volume of the organ.

    Returns:
        `urine_rate` (µl/min), the net urine flow `Q_urine`, and `urine`
        (µl), the urine volume `V_urine`, both assignment rules.
    """
    out = core.CoreObjects()
    ul = f"1000000 {PU.ul_per_l.sid}"
    _add(
        out,
        observer(
            ids.URINE_RATE,
            f"{ul} * {ids.Q_URINE}",
            PU.ul_per_min,
            "urine flow of the organ",
            SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
    )
    _add(
        out,
        observer(
            ids.URINE,
            f"{ul} * {ids.volume(ids.URINE_SITE)}",
            PU.ul,
            "urine volume of the organ",
            SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
    )
    return out
