"""The observers of a PBPK model: assignment rules of the top model, no states.

An observer reads the pool species and the outflow concentrations of the
organs and changes nothing, so it conserves mass. Per substance `<s>`:

- `C<pool>_free_<s> = fu*C<pool>_<s>`, the unbound plasma concentration of a
  pool (`ar`, `ve`, `po`);
- `C<pool>_blood_<s>`, the blood concentration of a pool,
  `(1-HCT)*C + HCT*C_bc`: `BP*C` unless the blood cells are kinetic;
- `Chv_<s>`, the hepatic vein: the outflow of the liver mixed with the
  portosystemic shunt by their flows,
  `((1-f_shunt)*Qh*Cout_li + f_shunt*(Qha*Car + Qpo*Cpo))/Qh`, with
  `Chv_free_<s>` and `Chv_blood_<s>`;
- `Carm_<s> = sum w*C`, the arm vein (peripheral venous) sampling site: the
  venous plasma outflows of organs and the arterial plasma mixed with the
  weights `w_arm_<site>` of a preset of `data/arm_vein.toml`
  (`PBPKDefinition.arm_vein`, by default the mix of skin and muscle of the
  PK-Sim documentation); `Carm_free_<s>` and `Carm_blood_<s>` as for the
  pools.

The tissue concentrations are the species of the organs themselves, their
unbound concentrations `Cu_<part>_<s>` observers of the organs, the amounts
in urine and feces the sink species.
"""

from sbmlutils.factory import AssignmentRule, Parameter, Sbase
from sbmlutils.metadata import SBO

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.annotation import annotate
from pkpdsim.pbpk.data import load, value
from pkpdsim.pbpk.definition import PBPKDefinition
from pkpdsim.pbpk.organs import Organ, Part, Pool
from pkpdsim.pbpk.placement import value_parameter
from pkpdsim.pbpk.substance import Substance
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Value

__all__ = ["arm_vein_weights", "observer_objects"]

_ONE = "1 dimensionless"
_SBO_PARAMETER = SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER
_OBSERVER = "An observer: an assignment rule which changes no state."


def arm_vein_weights(preset: str) -> dict[Organ | Pool, Value]:
    """Get the weights of a preset of the arm vein.

    Args:
        preset: the name of the preset, e.g. `pksim_default`.

    Returns:
        The weight of every organ outflow or pool of the mix with its source.

    Raises:
        ValueError: if the preset is unknown or names neither an organ nor a
            pool.
    """
    data = load("arm_vein")
    if preset not in data:
        raise ValueError(
            f"Unknown arm_vein '{preset}', the presets are {sorted(data)}."
        )
    return {_site(preset, key): value("arm_vein", preset, key) for key in data[preset]}


def _site(preset: str, key: str) -> Organ | Pool:
    """Get the organ or the pool of a key of a preset, `skin`, `arterial`."""
    for enum in (Organ, Pool):
        if key.upper() in enum.__members__:
            return enum[key.upper()]
    raise ValueError(
        f"arm_vein.toml [{preset}] {key}: neither an organ nor a pool of the body."
    )


def observer_objects(definition: PBPKDefinition) -> list[Sbase]:
    """Get the observers of a definition.

    Args:
        definition: the definition.

    Returns:
        The weights of the arm vein and the parameters with their assignment
        rules of every substance.
    """
    weights = arm_vein_weights(definition.arm_vein)
    objects: list[Sbase] = [
        annotate(
            value_parameter(
                ids.arm_weight(site),
                w,
                U.dimensionless,
                f"weight of {_site_label(site)} in the arm vein",
                _SBO_PARAMETER,
                f"Share of the {_site_label(site)} in the arm vein sample "
                f"(preset '{definition.arm_vein}').",
            ),
            "arm_weight",
            organ=site if isinstance(site, Organ) else None,
            pool=site if isinstance(site, Pool) else None,
        )
        for site, w in weights.items()
    ]
    for s in definition.substances:
        objects += _substance(s, weights)
    return objects


def _site_label(site: Organ | Pool) -> str:
    """The venous outflow of an organ or the plasma of a pool in a name."""
    if isinstance(site, Pool):
        return f"{site.name.lower()} plasma"
    return f"venous outflow of the {site.name.lower().replace('_', ' ')}"


def _observer(
    pid: str,
    formula: str,
    name: str,
    text: str,
    kind: str,
    pool: Pool | None = None,
    site: str | None = None,
) -> list[Sbase]:
    """A concentration (mM) with its assignment rule, annotated as `kind`."""
    return [
        annotate(
            Parameter(
                pid,
                None,
                U.mM,
                constant=False,
                name=name,
                sboTerm=_SBO_PARAMETER,
                notes=f"{text} {_OBSERVER}\n\n`{pid} = {formula}`",
            ),
            kind,
            pool=pool,
            site=site,
        ),
        AssignmentRule(pid, formula, U.mM),
    ]


def _blood(s: Substance, plasma: str, cells: str) -> str:
    """The blood concentration of the plasma and blood cell concentrations, as formula."""
    ery = s.erythrocytes
    if ery is not None and ery.mode == "kinetic":
        return f"({_ONE} - {ids.HCT}) * {plasma} + {ids.HCT} * {cells}"
    return f"{ids.bp(s.sid)} * {plasma}"


def _hepatic(s: Substance, cells: bool) -> str:
    """The hepatic vein mix of the liver outflow and the shunt, as formula.

    `((1-f_shunt)*Qh*Cout_li + f_shunt*(Qha*Car + Qpo*Cpo))/Qh` of the plasma,
    or of the kinetic blood cells.
    """
    sid = s.sid
    part = Part.BLOOD_CELLS if cells else Part.PLASMA
    liver = (
        ids.organ_outflow_bc(Organ.LIVER, sid)
        if cells
        else ids.organ_outflow(Organ.LIVER, sid)
    )
    car = ids.pool_species(Pool.ARTERIAL, sid, part)
    cpo = ids.pool_species(Pool.PORTAL, sid, part)
    return (
        f"(({_ONE} - {ids.F_SHUNT}) * {ids.QH} * {liver} + {ids.F_SHUNT} * "
        f"({ids.QHA} * {car} + {ids.QPO} * {cpo})) / {ids.QH}"
    )


def _substance(s: Substance, weights: dict[Organ | Pool, Value]) -> list[Sbase]:
    """The observers of a substance."""
    sid, fu = s.sid, ids.fu(s.sid)
    objects: list[Sbase] = []
    for pool in Pool:
        c = ids.pool_species(pool, sid)
        label = pool.name.lower()
        objects += _observer(
            ids.free(pool, sid),
            f"{fu} * {c}",
            f"unbound plasma concentration {sid} {label} [mM]",
            f"Unbound plasma concentration of the {label} pool.",
            "free_concentration",
            pool=pool,
        )
        objects += _observer(
            ids.blood(pool, sid),
            _blood(s, c, ids.pool_species(pool, sid, Part.BLOOD_CELLS)),
            f"blood concentration {sid} {label} [mM]",
            f"Blood concentration of the {label} pool.",
            "blood_concentration",
            pool=pool,
        )

    hv = ids.hepatic_vein(sid)
    objects += _observer(
        hv,
        _hepatic(s, cells=False),
        f"hepatic vein plasma concentration {sid} [mM]",
        "Plasma concentration of the hepatic vein: the outflow of the liver "
        "mixed with the portosystemic shunt by their flows.",
        "plasma_concentration",
        site="hepatic_vein",
    )
    objects += _observer(
        ids.hepatic_vein_free(sid),
        f"{fu} * {hv}",
        f"hepatic vein unbound plasma concentration {sid} [mM]",
        "Unbound plasma concentration of the hepatic vein.",
        "free_concentration",
        site="hepatic_vein",
    )
    ery = s.erythrocytes
    kinetic = ery is not None and ery.mode == "kinetic"
    objects += _observer(
        ids.hepatic_vein_blood(sid),
        _blood(s, hv, f"({_hepatic(s, cells=True)})" if kinetic else ""),
        f"hepatic vein blood concentration {sid} [mM]",
        "Blood concentration of the hepatic vein.",
        "blood_concentration",
        site="hepatic_vein",
    )

    def mix(cells: bool) -> str:
        terms = []
        for site in weights:
            if isinstance(site, Pool):
                c = ids.pool_species(
                    site, sid, Part.BLOOD_CELLS if cells else Part.PLASMA
                )
            elif cells:
                c = ids.organ_outflow_bc(site, sid)
            else:
                c = ids.organ_outflow(site, sid)
            terms.append(f"{ids.arm_weight(site)} * {c}")
        return " + ".join(terms)

    arm = ids.arm(sid)
    objects += _observer(
        arm,
        mix(cells=False),
        f"arm vein plasma concentration {sid} [mM]",
        "Plasma concentration of an arm vein sample: the venous outflows of the "
        "organs and the arterial plasma mixed with the weights of the preset. "
        "Plasma data without a stated sampling site map to it.",
        "plasma_concentration",
        site="arm_vein",
    )
    objects += _observer(
        ids.arm_free(sid),
        f"{fu} * {arm}",
        f"arm vein unbound plasma concentration {sid} [mM]",
        "Unbound plasma concentration of an arm vein sample.",
        "free_concentration",
        site="arm_vein",
    )
    cells = f"({mix(cells=True)})" if kinetic else ""
    objects += _observer(
        ids.arm_blood(sid),
        _blood(s, arm, cells),
        f"arm vein blood concentration {sid} [mM]",
        "Blood concentration of an arm vein sample.",
        "blood_concentration",
        site="arm_vein",
    )
    return objects
