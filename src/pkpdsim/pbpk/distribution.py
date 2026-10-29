"""The coefficients of the distribution between plasma and tissue in the standard organ.

The Kp of a substance is the whole-tissue (interstitial + intracellular) to
plasma coefficient in every form of the organ:

- full form: the interstitial to plasma coefficient follows from fu (Schmitt
  2008), `K_int_pls = (f_water_int + r_protein_int_pls*(1/fu - f_water_pls))*fu`,
  and the intracellular one is derived so that both together have the Kp,
  `K_cell_pls = (Kp*(f_int + f_cell) - f_int*K_int_pls)/f_cell`. A Kp at or
  below the interstitial share `f_int*K_int_pls/(f_int + f_cell)` has no
  positive `K_cell_pls`: the build refuses it for the reference individual
  of its sex, and the rule floors `K_cell_pls` at `K_CELL_FLOOR` should a fit
  or another individual reach it (the interstitial share is the lower bound of
  a fitted Kp);
- flow-limited organ: `K_or = (Vpl + Vbc*Kbc + Vti*Kp)/V`, the organ to plasma
  coefficient of the well-stirred organ.

The interstitial share `f_int*K_int_pls/(1 - f_vas)` of the reference
individual of the sex of the build (`interstitial_share`) is the Kp of an
extracellular substance, `Kp(interstitial=...)` scales it in every organ.

Every formula is written once as a template, which gives the SBML rule (with
ids) and the number (with values) alike.
"""

from functools import cache

import libsbml
from sbmlutils.factory import Parameter, UnitDefinition
from sbmlutils.metadata import SBO

from pkpdsim.pbpk.annotation import annotate
from pkpdsim.pbpk.data import load
from pkpdsim.pbpk.kp import kp_values
from pkpdsim.pbpk.organs import Organ
from pkpdsim.pbpk.physiology import Sex, reference_fractions
from pkpdsim.pbpk.placement import value_parameter
from pkpdsim.pbpk.substance import Kp, Substance
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Formula, Q, Value

__all__ = [
    "INTERSTITIAL",
    "K_CELL_FLOOR",
    "K_VAS",
    "check_kp",
    "constant",
    "constants",
    "interstitial_share",
    "k_cell_pls",
    "k_cell_pls_value",
    "k_int_pls",
    "k_int_pls_value",
    "k_organ",
    "kp_of_organ",
    "number",
]

INTERSTITIAL = ("f_water_int", "f_water_pls", "r_protein_int_pls")
"""The local ids of the constants of `K_int_pls`, keys of `data/distribution.toml`."""

K_VAS = "k_vas"
"""The local id of the exchange between plasma and interstitium per organ volume."""

K_CELL_FLOOR = 1e-6
"""The lower bound of `K_cell_pls` in the rule, against a Kp below the interstitial share."""

_ONE = "1 dimensionless"

_K_INT = "({f_water_int} + {r_protein_int_pls} * ({one} / {fu} - {f_water_pls})) * {fu}"
_K_CELL = "({kp} * ({f_int} + {f_cell}) - {f_int} * {k_int}) / {f_cell}"
_F_CELL = "({one} - {f_vas} - {f_int})"

_DATA_UNITS: dict[str, UnitDefinition] = {
    "dimensionless": U.dimensionless,
    "1/min": U.per_min,
}
"""The units of the model for the units of `data/distribution.toml`."""

_NAMES = {
    "f_water_int": "water fraction interstitial fluid",
    "f_water_pls": "water fraction plasma",
    "r_protein_int_pls": "protein of the interstitial fluid relative to plasma",
    K_VAS: "exchange plasma and interstitium per organ volume [1/min]",
}


@cache
def constants() -> dict[str, Value]:
    """The constants of `data/distribution.toml` by key, values with their source."""
    return {
        key: Value(e["value"], _DATA_UNITS[e["unit"]], e["source"], e.get("notes"))
        for table in load("distribution").values()
        for key, e in table.items()
    }


def constant(key: str) -> Parameter:
    """The parameter of a constant of `data/distribution.toml`, with its source."""
    value = constants()[key]
    assert isinstance(value.unit, UnitDefinition)
    return annotate(
        value_parameter(
            key,
            value,
            value.unit,
            _NAMES[key],
            SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        "endothelial_exchange" if key == K_VAS else "tissue_composition",
    )


def number(q: Q) -> float:
    """The number of a plain float or a `Value`.

    Raises:
        TypeError: for a `Formula`, which has no number before the simulation.
    """
    if isinstance(q, Formula):
        raise TypeError(f"The formula `{q.expression}` has no number.")
    return q.value if isinstance(q, Value) else q


def _evaluate(formula: str) -> float:
    """The value of a formula of numbers."""
    ast = libsbml.parseL3Formula(formula)
    if ast is None:
        raise ValueError(f"Cannot parse '{formula}': {libsbml.getLastParseL3Error()}")
    return float(libsbml.SBMLTransforms.evaluateASTNode(ast))


def k_int_pls(fu: str) -> str:
    """The rule of the interstitial to plasma coefficient (Schmitt 2008).

    Args:
        fu: the id of the unbound fraction in plasma.

    Returns:
        The formula with the ids of `INTERSTITIAL`.
    """
    return _K_INT.format(fu=fu, one=_ONE, **{k: k for k in INTERSTITIAL})


def k_int_pls_value(fu: float) -> float:
    """The interstitial to plasma coefficient of an unbound fraction."""
    values = {k: repr(constants()[k].value) for k in INTERSTITIAL}
    return _evaluate(_K_INT.format(fu=repr(fu), one="1", **values))


def interstitial_share(organ: Organ, fu: float, sex: Sex = Sex.MALE) -> float:
    """The interstitial share of the whole tissue of an organ.

    `f_int*K_int_pls(fu)/(1 - f_vas)`, the Kp of a substance which reaches the
    whole interstitial fluid of the organ and not its cells, with the vascular
    and interstitial fractions of the reference individual of the sex
    (`physiology.reference_fractions`): the lower bound of the Kp of the full
    form (`check_kp`) and the coefficients of `Kp(interstitial=...)` of a
    build of the sex.

    Args:
        organ: the organ.
        fu: the fraction unbound in plasma.
        sex: the sex of the build, male by default (the default of
            `PBPKDefinition`).

    Returns:
        The interstitial share (dimensionless).
    """
    f_vas, f_int = reference_fractions(organ)[sex]
    return f_int * k_int_pls_value(fu) / (1.0 - f_vas)


def kp_of_organ(substance: Substance, organ: Organ, sex: Sex) -> Q | None:
    """The whole-tissue Kp which a `Kp` gives an organ, before the simulation.

    Args:
        substance: the substance, with a `Kp` as distribution.
        organ: the organ.
        sex: the sex of the build, whose interstitial share a scale multiplies.

    Returns:
        The override of the organ, else the scale times the interstitial share
        (a number, `None` for a scale or an fu which is a `Formula`), else the
        value of the substance.

    Raises:
        ValueError: if the distribution of the substance is no `Kp`.
    """
    dist = substance.distribution
    if not isinstance(dist, Kp):
        raise ValueError(f"{substance.sid}: the distribution is no Kp.")
    if organ in dist.overrides:
        return dist.overrides[organ]
    if dist.interstitial is None:
        return dist.value
    fu = substance.binding.fu if substance.binding is not None else 1.0
    if isinstance(fu, Formula) or isinstance(dist.interstitial, Formula):
        return None
    return number(dist.interstitial) * interstitial_share(organ, number(fu), sex)


def k_cell_pls(kp: str, k_int: str, f_vas: str, f_int: str) -> str:
    """The rule of the intracellular to plasma coefficient, floored at `K_CELL_FLOOR`.

    Args:
        kp: the id of the whole-tissue to plasma coefficient.
        k_int: the id of the interstitial to plasma coefficient.
        f_vas: the id of the vascular fraction of the organ volume.
        f_int: the id of the interstitial fraction of the organ volume.

    Returns:
        `piecewise(K, K > floor, floor)` with `K` such that interstitium and
        cells together have the Kp.
    """
    f_cell = _F_CELL.format(one=_ONE, f_vas=f_vas, f_int=f_int)
    k = _K_CELL.format(kp=kp, k_int=k_int, f_int=f_int, f_cell=f_cell)
    floor = f"{K_CELL_FLOOR!r} dimensionless"
    return f"piecewise({k}, {k} > {floor}, {floor})"


def k_cell_pls_value(kp: float, fu: float, f_vas: float, f_int: float) -> float:
    """The intracellular to plasma coefficient of the numbers, without the floor."""
    f_cell = _F_CELL.format(one="1", f_vas=repr(f_vas), f_int=repr(f_int))
    return _evaluate(
        _K_CELL.format(
            kp=repr(kp),
            k_int=repr(k_int_pls_value(fu)),
            f_int=repr(f_int),
            f_cell=f_cell,
        )
    )


def k_organ(f_vas: str, vascular: str, kp: str) -> str:
    """The rule of the organ to plasma coefficient of the flow-limited organ.

    Args:
        f_vas: the id of the vascular fraction of the organ volume.
        vascular: the capacity of the blood relative to its plasma, `1 - HCT`
            (no erythrocyte distribution) or `BP` (erythrocytes in
            equilibrium), as formula.
        kp: the id of the whole-tissue to plasma coefficient.

    Returns:
        `f_vas*vascular + (1 - f_vas)*Kp`, i.e. `(Vpl + Vbc*Kbc + Vti*Kp)/V`.
    """
    return f"{f_vas} * {vascular} + ({_ONE} - {f_vas}) * {kp}"


def check_kp(substance: Substance, organ: Organ, sex: Sex) -> None:
    """Check that the Kp of a substance leaves a positive `K_cell_pls` in an organ.

    The check uses the fractions of the reference individual of the sex of
    the build and the fu of the definition; the Kp of a method is computed
    for the sex (its reference hematocrit) and checked with the factor of the
    method. A build of the other sex is checked when it is built.

    Args:
        substance: the substance.
        organ: the organ.
        sex: the sex of the build.

    Raises:
        ValueError: if the Kp is not above the interstitial share of the
            whole tissue, `f_int*K_int_pls/(f_int + f_cell)`.
    """
    dist = substance.distribution
    if dist is None:
        return
    fu_q = substance.binding.fu if substance.binding is not None else 1.0
    kp_q = kp_of_organ(substance, organ, sex) if isinstance(dist, Kp) else dist.factor
    if isinstance(fu_q, Formula) or kp_q is None or isinstance(kp_q, Formula):
        # a formula has no number before the simulation; `K_cell_pls` has the
        # lower bound 1e-6 in the model
        return
    fu = number(fu_q)
    f_vas, f_int = reference_fractions(organ)[sex]
    if isinstance(dist, Kp):
        kp = number(kp_q)
    else:
        kp = kp_values(substance, sex)[organ].value * number(dist.factor)
    if k_cell_pls_value(kp, fu, f_vas, f_int) <= 0:
        share = f_int * k_int_pls_value(fu) / (1.0 - f_vas)
        raise ValueError(
            f"{substance.sid}: Kp {kp:.4g} in organ '{organ.value}' is not "
            f"above the interstitial share {share:.4g} of the tissue "
            f"({sex.value} reference individual), the full form would need "
            f"K_cell_pls <= 0."
        )
