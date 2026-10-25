"""The tissue to plasma coefficients of the Kp methods.

A Kp method computes the whole-tissue (interstitial + intracellular) to
plasma coefficient of every organ when the model is built (`KpMethod` of a
substance). The builder writes one `Kp_<o>_<s>` per organ, an initial
assignment `fKp_<s> * value` with the factor `fKp_<s>` of the method, so that a
fit scales every coefficient at once.

Rodgers & Rowland (Rodgers, Leahy & Rowland 2005 for moderate-to-strong bases,
Rodgers & Rowland 2006 for acids, very weak bases, neutrals and zwitterions)
predict the tissue to unbound plasma coefficient `Kpu` (`Kp = Kpu * fu`) in the
form stated by Yau 2020 (eqs. 1 to 5), from the composition of
`data/tissue_composition.toml`:

- moderate-to-strong bases (a basic pKa >= 7), electrostatic binding to the
  acidic phospholipids `AP` of the tissue (eq. 1),
  `Kpu = X*f_iw/Y + f_ew + (P*f_nl + (0.3*P + 0.7)*f_np)/Y + Ka_AP*AP*(X - 1)/Y`,
  with the affinity from the erythrocytes (eqs. 3, 5),
  `Kpu_bc = (HCT - 1 + BP)/(fu*HCT)`,
  `Ka_AP = (Kpu_bc - X_bc*f_iw_bc/Y - (P*f_nl_bc + (0.3*P + 0.7)*f_np_bc)/Y)
  * Y/(AP_bc*(X_bc - 1))`;
- acids, very weak bases and zwitterions (albumin) and neutrals
  (lipoproteins), binding to the extracellular protein `PR` (eqs. 2, 4),
  `Kpu = X*f_iw/Y + f_ew + (P*f_nl + (0.3*P + 0.7)*f_np)/Y + Ka_PR*[PR]_T`,
  `Ka_PR*[PR]_T = (1/fu - 1 - (P*f_nl_p + (0.3*P + 0.7)*f_np_p)/Y) * [PR]_T/[PR]_P`;
- `P` is the octanol:water partition coefficient, the vegetable oil:water one
  `log P_vo:w = 1.115*log P - 1.35` in the lipid terms of adipose; `X` and `Y`
  are the ionisation terms of intracellular water and plasma, `1 + 10^(pKa -
  pH)` for a base, `1 + 10^(pH - pKa)` for an acid, 1 for a neutral; a
  negative affinity is set to zero (Yau 2020, citing Rodgers & Rowland 2006).

Choices of the library (not verified against the 2005/2006 papers, which were
not accessible): several pKa of an acid or a base give the ionisation of
successive dissociation, `1 + 10^(pH - pKa1) + 10^(2*pH - pKa1 - pKa2)` for an
acid (pKa ascending) and `1 + 10^(pKa1 - pH) + 10^(pKa1 + pKa2 - 2*pH)` for a
base (pKa descending); a zwitterion (`pka = (acidic, basic)`) has `1 + 10^(pH -
pKa_acid) + 10^(pKa_base - pH)` and binds to the acidic phospholipids if its
basic pKa is >= 7, to albumin otherwise.

The method applies the rat tissue composition to human tissues (Rodgers &
Rowland 2007, Coutinho 2024); the inputs are the substance's reference `fu`
(`Binding.fu`), its blood to plasma ratio (`Erythrocytes.bp`, or `1 - HCT +
HCT*kbc`) and the reference hematocrit of the sex. The coefficients are
numbers of the build: a change of `fu` by an impairment (the binding factor of
`pkpdsim.pbpk.impairment`) or of the hematocrit does not recompute them.
"""

from collections.abc import Mapping
from dataclasses import dataclass
from functools import cache
from typing import Any, Literal

from pkpdsim.pbpk.data import load
from pkpdsim.pbpk.organs import Organ
from pkpdsim.pbpk.physiology import Sex, reference_hematocrit
from pkpdsim.pbpk.substance import KpMethod, Substance
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Q, Value

__all__ = [
    "AP_PKA",
    "Composition",
    "binds_phospholipids",
    "blood_to_plasma",
    "check_kp_method",
    "composition",
    "kp_sources",
    "kp_values",
    "rodgers_rowland",
]

AP_PKA = 7.0
"""A base with a pKa at or above it binds to the acidic phospholipids (Yau 2020)."""

_KEYS = ("f_ew", "f_iw", "f_nl", "f_np", "albumin", "lipoprotein", "ap")
_UNITS = {"dimensionless", "mg/g"}


@dataclass(frozen=True)
class Composition:
    """The composition of a tissue of the Rodgers & Rowland method.

    Attributes:
        tissue: the tissue of `data/tissue_composition.toml`.
        f_ew: the fractional volume of extracellular water.
        f_iw: the fractional volume of intracellular water.
        f_nl: the fractional volume of neutral lipids.
        f_np: the fractional volume of neutral phospholipids.
        albumin: the tissue to plasma ratio of albumin.
        lipoprotein: the tissue to plasma ratio of lipoproteins.
        ap: the acidic phospholipids (mg/g tissue).
        sources: the keys of the sources of the numbers.
    """

    tissue: str
    f_ew: float
    f_iw: float
    f_nl: float
    f_np: float
    albumin: float
    lipoprotein: float
    ap: float
    sources: tuple[str, ...] = ()


def _number(table: Mapping[str, Any], key: str) -> float:
    """A number of the composition file.

    Raises:
        ValueError: if its unit is not one of the file.
    """
    entry = table[key]
    if entry["unit"] not in _UNITS:
        raise ValueError(
            f"tissue_composition.toml: '{key}' has the unit '{entry['unit']}', "
            f"not one of {sorted(_UNITS)}."
        )
    return float(entry["value"])


@cache
def _data() -> dict[str, Any]:
    """The content of `data/tissue_composition.toml`."""
    return load("tissue_composition")


@cache
def _compositions() -> dict[str, Composition]:
    """The compositions of the tissues by name."""
    found: dict[str, Composition] = {}
    for tissue, table in _data()["tissues"].items():
        found[tissue] = Composition(
            tissue,
            **{key: _number(table, key) for key in _KEYS},
            sources=tuple(dict.fromkeys(e["source"] for e in table.values())),
        )
    return found


def composition(organ: Organ) -> Composition:
    """Get the composition of the tissue of an organ.

    Args:
        organ: the organ.

    Returns:
        The composition of its tissue, a surrogate tissue for an organ
        without data of its own (`[organs]` of the data file).
    """
    return _compositions()[_data()["organs"][organ.value]["tissue"]]


def kp_sources(organ: Organ) -> tuple[str, ...]:
    """Get the keys of the sources of the Rodgers & Rowland Kp of an organ.

    Args:
        organ: the organ.

    Returns:
        The equations (Yau 2020), the composition of the tissue, the choice of
        the tissue, plasma and erythrocytes.
    """
    data = _data()
    keys = [
        "Yau2020",
        *composition(organ).sources,
        data["organs"][organ.value]["source"],
    ]
    for section in ("plasma", "intracellular", "erythrocytes", "adipose"):
        keys += [e["source"] for e in data[section].values()]
    return tuple(dict.fromkeys(keys))


def _constant(section: str, key: str) -> float:
    """A constant of plasma, intracellular water, erythrocytes or adipose."""
    return _number(_data()[section], key)


def _ionisation(
    compound: Literal["neutral", "acid", "base", "zwitterion"],
    pka: tuple[float, ...],
    ph: float,
) -> float:
    """The total to neutral concentration of a compound at a pH, `X` and `Y`."""
    if compound == "neutral":
        return 1.0
    if compound == "zwitterion":
        acidic, basic = pka
        return 1.0 + 10 ** (ph - acidic) + 10 ** (basic - ph)
    total, term = 1.0, 0.0
    ordered = sorted(pka) if compound == "acid" else sorted(pka, reverse=True)
    for k, p in enumerate(ordered, start=1):
        term += p
        total += 10 ** (k * ph - term) if compound == "acid" else 10 ** (term - k * ph)
    return total


def binds_phospholipids(method: KpMethod) -> bool:
    """Whether a substance binds to the acidic phospholipids (eq. 1 of Yau 2020).

    Args:
        method: the Kp method of the substance.

    Returns:
        True for a base with a pKa >= 7 and a zwitterion with a basic pKa >= 7.
    """
    if method.compound == "base":
        return bool(method.pka) and max(method.pka) >= AP_PKA
    if method.compound == "zwitterion":
        return len(method.pka) == 2 and method.pka[1] >= AP_PKA
    return False


def rodgers_rowland(
    substance: KpMethod,
    organ: Organ,
    fu: float,
    bp: float | None,
    hct: float,
    *,
    tissue: Composition | None = None,
) -> float:
    """Get the whole-tissue to plasma coefficient of Rodgers & Rowland.

    Args:
        substance: the Kp method of the substance (logP, pKa, compound).
        organ: the organ.
        fu: the fraction unbound in plasma.
        bp: the blood to plasma ratio, required for a substance which binds to
            the acidic phospholipids.
        hct: the hematocrit.
        tissue: the composition, the tissue of the organ by default.

    Returns:
        `Kp = Kpu * fu`, without the factor of the method.

    Raises:
        ValueError: if a substance which binds to the acidic phospholipids has
            no blood to plasma ratio.
    """
    c = tissue or composition(organ)
    compound, pka = substance.compound, tuple(substance.pka)
    p = 10**substance.logp
    p_tissue = p
    if c.tissue == "adipose":
        p_tissue = 10 ** (
            _constant("adipose", "vo_slope") * substance.logp
            + _constant("adipose", "vo_intercept")
        )
    x = _ionisation(compound, pka, _constant("intracellular", "ph"))
    y = _ionisation(compound, pka, _constant("plasma", "ph"))
    kpu = (
        x * c.f_iw / y
        + c.f_ew
        + (p_tissue * c.f_nl + (0.3 * p_tissue + 0.7) * c.f_np) / y
    )
    if binds_phospholipids(substance):
        if bp is None:
            raise ValueError(
                "Rodgers & Rowland: a substance which binds to the acidic "
                "phospholipids (a base with pKa >= 7) needs its blood to plasma "
                "ratio."
            )
        kpu_bc = (hct - 1.0 + bp) / (hct * fu)
        x_bc = _ionisation(compound, pka, _constant("erythrocytes", "ph"))
        lipids_bc = (
            p * _constant("erythrocytes", "f_nl")
            + (0.3 * p + 0.7) * _constant("erythrocytes", "f_np")
        ) / y
        ka = (
            (kpu_bc - x_bc * _constant("erythrocytes", "f_iw") / y - lipids_bc)
            * y
            / (_constant("erythrocytes", "ap") * (x_bc - 1.0))
        )
        kpu += max(ka, 0.0) * c.ap * (x - 1.0) / y
    else:
        ratio = c.lipoprotein if compound == "neutral" else c.albumin
        lipids_p = (
            p * _constant("plasma", "f_nl")
            + (0.3 * p + 0.7) * _constant("plasma", "f_np")
        ) / y
        kpu += max(1.0 / fu - 1.0 - lipids_p, 0.0) * ratio
    return kpu * fu


def _number_of(q: Q) -> float:
    """The number of a plain float or a `Value`."""
    return q.value if isinstance(q, Value) else q


def blood_to_plasma(substance: Substance, hct: float) -> float | None:
    """Get the blood to plasma ratio of a substance at a hematocrit.

    Args:
        substance: the substance.
        hct: the hematocrit.

    Returns:
        `bp` of the erythrocytes, `1 - HCT + HCT*kbc` of a blood cell to plasma
        coefficient, `None` without erythrocyte distribution.
    """
    ery = substance.erythrocytes
    if ery is None:
        return None
    if ery.bp is not None:
        return _number_of(ery.bp)
    if ery.kbc is not None:
        return 1.0 - hct + hct * _number_of(ery.kbc)
    return None


def kp_values(substance: Substance, sex: Sex) -> dict[Organ, Value]:
    """Get the coefficients of the Kp method of a substance for the organs it enters.

    Args:
        substance: the substance, with a `KpMethod` as distribution.
        sex: the sex of the model, whose reference hematocrit enters the
            affinity of a base to the acidic phospholipids.

    Returns:
        The whole-tissue to plasma coefficient of every organ the substance
        enters, without the factor of the method, with its inputs as notes.

    Raises:
        ValueError: if the distribution of the substance is no Kp method.
    """
    method = substance.distribution
    if not isinstance(method, KpMethod):
        raise ValueError(f"{substance.sid}: the distribution is no Kp method.")
    fu = _number_of(substance.binding.fu) if substance.binding is not None else 1.0
    hct = reference_hematocrit(sex)
    bp = blood_to_plasma(substance, hct)
    organs = [o for o in Organ if method.organs is None or o in method.organs]
    binding = (
        "acidic phospholipids"
        if binds_phospholipids(method)
        else ("lipoproteins" if method.compound == "neutral" else "albumin")
    )
    values: dict[Organ, Value] = {}
    for organ in organs:
        tissue = composition(organ).tissue
        kp = rodgers_rowland(method, organ, fu, bp, hct)
        values[organ] = Value(
            kp,
            U.dimensionless,
            "Yau2020",
            f"Rodgers & Rowland Kp of the tissue '{tissue}' ({method.compound}, "
            f"binding to {binding}): logP {method.logp}, pKa {list(method.pka)}, "
            f"fu {fu}, BP {bp}, HCT {hct} ({sex.value} reference individual). "
            "The reference fu: an impairment which changes fu or HCT does not "
            "recompute the coefficient.",
        )
    return values


def check_kp_method(substance: Substance) -> list[str]:
    """Check the Kp method of a substance.

    Args:
        substance: the substance.

    Returns:
        The problems: the number of pKa values of the compound type, a missing
        blood to plasma ratio of a substance which binds to the acidic
        phospholipids.
    """
    method = substance.distribution
    if not isinstance(method, KpMethod):
        return []
    sid, n = substance.sid, len(method.pka)
    problems: list[str] = []
    if method.compound == "neutral" and n:
        problems.append(f"{sid}: a neutral compound has no pKa, got {list(method.pka)}")
    if method.compound in ("acid", "base") and not n:
        problems.append(f"{sid}: an acid or a base needs at least one pKa")
    if method.compound == "zwitterion" and n != 2:
        problems.append(
            f"{sid}: a zwitterion needs two pKa values (acidic, basic), got "
            f"{list(method.pka)}"
        )
    ery = substance.erythrocytes
    if (
        not problems
        and binds_phospholipids(method)
        and (ery is None or (ery.bp is None and ery.kbc is None))
    ):
        problems.append(
            f"{sid}: Rodgers & Rowland needs the blood to plasma ratio of a base "
            "with pKa >= 7 (the affinity to the acidic phospholipids follows from "
            "the erythrocytes, Rodgers, Leahy & Rowland 2005): give "
            "Erythrocytes(bp=...) or a kbc"
        )
    return problems
