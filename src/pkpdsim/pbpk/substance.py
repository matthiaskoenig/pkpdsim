"""Substances and how the body treats them: binding, erythrocytes, distribution."""

from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import StrEnum
from typing import Literal

from sbmlutils.metadata import BQB

from pkpdsim.pbpk.organs import Organ
from pkpdsim.pbpk.proteins import Induction, Inhibition
from pkpdsim.pbpk.values import Q

__all__ = [
    "Binding",
    "Erythrocytes",
    "Excretion",
    "Kp",
    "KpMethod",
    "Mode",
    "Substance",
]


@dataclass(frozen=True)
class Binding:
    """The binding to plasma proteins.

    Attributes:
        fu: the fraction unbound in plasma.
        protein: the binding protein.
    """

    fu: Q
    protein: Literal["albumin", "agp", "unknown"] = "unknown"


@dataclass(frozen=True)
class Erythrocytes:
    """The distribution into the erythrocytes.

    Attributes:
        mode: `equilibrium` (instantaneous, given by `bp` or `kbc`) or
            `kinetic` (explicit blood cell species with exchange rates).
        bp: the blood to plasma ratio (equilibrium).
        kbc: the blood cell to plasma coefficient (equilibrium or kinetic).
        k_exchange: the exchange rate between plasma and blood cells (kinetic).
    """

    mode: Literal["equilibrium", "kinetic"]
    bp: Q | None = None
    kbc: Q | None = None
    k_exchange: Q | None = None


class Mode(StrEnum):
    """How a substance is exchanged between plasma and tissue."""

    PERFUSION = "perfusion"
    """Flow-limited: the tissue is in equilibrium with the blood leaving it."""
    TISSUE = "tissue"
    """Merged tissue with finite exchange `PS * fu * (Cpl - Cti/Kp)`."""
    PERMEABILITY = "permeability"
    """Interstitial and intracellular space with a permeability each."""


@dataclass(frozen=True, kw_only=True)
class _Exchange:
    """What the Kp and the Kp method share.

    Attributes:
        organs: the organs the substance enters, all if `None`.
        mode: the exchange between plasma and tissue.
        ftissue: the whole-body exchange (l/min), required for the modes
            `TISSUE` and `PERMEABILITY`.
        ps: the exchange of single organs (l/min) instead of the share of
            `ftissue`.
    """

    organs: frozenset[Organ] | None = None
    mode: Mode = Mode.TISSUE
    ftissue: Q | None = None
    ps: Mapping[Organ, Q] = field(default_factory=dict)


@dataclass(frozen=True)
class Kp(_Exchange):
    """A tissue to plasma coefficient of the whole tissue.

    Attributes:
        value: the coefficient of the substance.
        overrides: the coefficients of single organs.
    """

    value: Q
    overrides: Mapping[Organ, Q] = field(default_factory=dict)


@dataclass(frozen=True)
class KpMethod(_Exchange):
    """The tissue to plasma coefficients from a calculation method.

    Attributes:
        method: the method.
        logp: the log10 of the octanol:water partition coefficient.
        pka: the pKa values.
        compound: the type of the compound.
        factor: a factor on all coefficients.
    """

    method: Literal["rodgers_rowland"]
    logp: float
    pka: tuple[float, ...]
    compound: Literal["neutral", "acid", "base", "zwitterion"]
    factor: Q = 1.0


@dataclass(frozen=True)
class Excretion:
    """The excretion which is not a process of the chemistry.

    Attributes:
        urine: glomerular filtration `fu*GFR`.
        fraction_reabsorbed: the fraction of the filtrate which is reabsorbed.
    """

    urine: bool = False
    fraction_reabsorbed: Q = 0.0


@dataclass(frozen=True)
class Substance:
    """A substance of the model.

    Attributes:
        sid: the short id, unique over all models (metabolites take the
            prefix of the parent, `api_m2`).
        name: the name.
        mr: the molecular weight (g/mole).
        annotations: the annotations (ChEBI, PubChem, NCIT).
        binding: the plasma protein binding, `None` is fu = 1.
        erythrocytes: the distribution into erythrocytes, `None` is none.
        distribution: the distribution into tissues, `None` is vascular only.
        excretion: the excretion.
        inhibits: the inhibitions of proteins by the substance.
        induces: the inductions of proteins by the substance.
    """

    sid: str
    name: str
    mr: Q
    annotations: list[tuple[BQB, str]] = field(default_factory=list)
    binding: Binding | None = None
    erythrocytes: Erythrocytes | None = None
    distribution: Kp | KpMethod | None = None
    excretion: Excretion = field(default_factory=Excretion)
    inhibits: tuple[Inhibition, ...] = ()
    induces: tuple[Induction, ...] = ()
