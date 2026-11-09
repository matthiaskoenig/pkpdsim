"""The transport coefficients of a substance in the zonated liver.

A new substance needs only its molar mass: the aqueous diffusion coefficient is
derived from it, the cytosolic one is a fraction of it and the sieving of the
fenestrae follows from the radius of the molecule. A `Substance.diffusion`
overrides the derived value. The units are SI (m, m**2/s).

Correlation
-----------
The aqueous diffusion coefficient at 37 C is the Stokes-Einstein scaling with
a hydrodynamic radius which grows with a power of the molar mass,

    D = D_REF * (Mr / MR_REF) ** -EXPONENT,    r = kT / (6 pi eta D),

with T = 310.15 K and the viscosity of water at 37 C, `ETA` = 0.692 mPa s
(Kestin, Sokolov and Wakeham 1978). `EXPONENT` = 0.4 lies between the
Stokes-Einstein scaling of a compact sphere (1/3, Young, Carroad and Bell
1980) and the 1/2 of an extended chain (Polson 1950); a single exponent of 1/3
misses the anchors by 1.15 (the Gd chelate) to 1.41 (insulin), albumin by
1.35; 0.4 keeps every anchor within 1.15. `D_REF` and `MR_REF` are the glucose
anchor. `test_transport.py` checks that the four sourced anchors are reproduced
within a factor of 1.5 (insulin within 2).

Anchors, D in 1e-9 m**2/s at 37 C. Where the measurement is at another
temperature it is scaled with `T / eta` of water (Stokes-Einstein: x 1.338 from
25 C with eta 0.890 and 0.692 mPa s, x 1.532 from 20 C with eta 1.002 mPa s):

- glucose, 180.16 g/mole: 6.7e-10 m2/s at 25 C (table of diffusion
  coefficients in water at 25 C, biofilmbook Table 4-1, hypertextbookshop.com),
  0.90 at 37 C.
- sucrose, 342.30 g/mole: 5.2e-10 m2/s at 25 C (same table), 0.70 at 37 C.
- Gd-DOTA, 558.6 g/mole, the comparable Gd chelate: Ga-DOTA in D2O at 310 K,
  4.38 +- 0.04e-10 m2/s by NMR diffusometry, valid for Gd-DOTA as the two are
  hydrodynamically similar (ISMRM 2014, abstract 2618). Corrected from D2O to
  H2O by the assumed viscosity ratio of about 1.23: 0.54 at 37 C. Consistency
  check, not the anchor: 2.92 +- 0.25e-10 in blood plasma, same abstract. No
  measured diffusion coefficient of gadoxetate itself was found.
- human serum albumin, 66.5 kDa: D(20,w) = 6.06 +- 0.07e-11 m2/s (bovine serum
  albumin, dynamic light scattering, extrapolated to water at 20 C and
  infinite dilution; Bio-protocol e4195, doi 10.21769/BioProtoc.4195), 0.093 at
  37 C; the radius it gives, 3.5 nm, is the accepted one of albumin.
- insulin, 5808 g/mole: 0.20 at 37 C, from 1.5e-10 at 25 C. Not confirmed by a
  source: insulin self-associates, so its measured D depends on concentration,
  pH and zinc, monomer and hexamer differ by a factor of about 1.6. It is kept
  as an indicative anchor with a factor of 2 in the test.

The exponent is an interpolation between the two cited scalings, not a fit.

The cytosolic diffusion of a small solute is about a quarter of the aqueous
one: Seksek, Biwersi and Verkman 1997 (J Cell Biol 138:131) measured
D/D0 = 0.26 +- 0.01 in MDCK cells and 0.27 +- 0.01 in fibroblasts by FRAP,
independent of the size of FITC-dextran and Ficoll; `F_CYT` = 0.27. Dix and
Verkman 2008 (Annu Rev Biophys 37:247) found the slowing less marked than
generally assumed.

The sieving of the fenestrae is the Renkin 1954 (J Gen Physiol 38:225)
partition and hindrance of a spherical solute in a cylindrical pore.
"""

import math
from typing import Protocol, runtime_checkable

from sbmlutils.factory import Parameter
from sbmlutils.factory.units import UnitDefinition
from sbmlutils.metadata import SBO

from pkpdsim.sbml.parameters import value_parameter
from pkpdsim.sbml.units import U
from pkpdsim.sbml.values import DIFFUSION_UNITS, Q, Value, number, unit_text
from pkpdsim.zonation import ids
from pkpdsim.zonation.geometry import SinusoidGeometry

__all__ = [
    "CM2_PER_S",
    "ETA",
    "F_CYT",
    "L_PER_M3",
    "M2_PER_S",
    "M_PER_S",
    "Solute",
    "aqueous_diffusion",
    "axial_velocity",
    "cytosolic_diffusion",
    "diffuses",
    "diffusion_of",
    "diffusion_parameter",
    "dispersion",
    "fenestral_permeability",
    "molecular_radius",
    "peclet",
    "renkin",
]

M2_PER_S = UnitDefinition("m2_per_s", "meter^2/second")
"""The unit of a diffusion coefficient."""

M_PER_S = UnitDefinition("m_per_s", "meter/second")
"""The unit of a permeability."""

CM2_PER_S = UnitDefinition("cm2_per_s", "cm^2/second")
"""A diffusion coefficient as it is often published."""

L_PER_M3 = UnitDefinition("l_per_m3", "litre/metre^3", name="litre per cubic metre")
"""The conversion of a volume in cubic metres into litres, which every
conductance of a permeability times an area needs (`P*A*1000 l_per_m3`)."""

T_BODY = 310.15
"""The body temperature (K)."""

ETA = 0.692e-3
"""The viscosity of water at 37 C (Pa s), Kestin, Sokolov and Wakeham 1978."""

K_B = 1.380649e-23
"""The Boltzmann constant (J/K)."""

D_REF = 0.90e-9
"""The aqueous diffusion of the reference, glucose, at 37 C (m**2/s)."""

MR_REF = 180.16
"""The molar mass of the reference, glucose (g/mole)."""

EXPONENT = 0.4
"""The exponent of the molar mass of the diffusion coefficient."""

F_CYT = 0.27
"""The cytosolic diffusion as a fraction of the aqueous one, Seksek, Biwersi and Verkman 1997."""

CORRELATION = (
    f"D = {D_REF:g} m2/s * (Mr / {MR_REF:g} g/mole)^-{EXPONENT:g} at 37 C "
    "(Stokes-Einstein scaling, anchors in pkpdsim.zonation.transport)"
)


MR_UNITS: dict[str, float] = {
    "g/mole": 1.0,
    "g/mol": 1.0,
    "kg/mole": 1e3,
    "kg/mol": 1e3,
}
"""The units of a molar mass, the factor to g/mole, by `unit_text`."""

LENGTH_UNITS: dict[str, float] = {
    "metre": 1.0,
    "meter": 1.0,
    "m": 1.0,
    "nanometer": 1e-9,
    "nm": 1e-9,
    "micrometer": 1e-6,
}
"""The units of a radius, the factor to m, by `unit_text`."""


def _mr_number(mr: Q) -> float:
    """The molar mass (g/mole) of a plain number (g/mole) or a `Value`, converted."""
    x = number(mr)
    if x is None or not x > 0.0:
        raise ValueError(f"The molar mass must be a positive number, got {mr!r}.")
    if isinstance(mr, Value):
        try:
            return x * MR_UNITS[unit_text(mr.unit)]
        except KeyError:
            raise ValueError(
                f"A molar mass needs a mass per amount, got {unit_text(mr.unit)}."
            ) from None
    return x


def aqueous_diffusion(mr: Q) -> Value:
    """Get the diffusion coefficient in water at 37 C derived from the molar mass.

    Args:
        mr: the molar mass (g/mole), a plain number or a `Value`.

    Returns:
        The diffusion coefficient (m**2/s).
    """
    d = D_REF * (_mr_number(mr) / MR_REF) ** -EXPONENT
    return Value(d, M2_PER_S, notes=CORRELATION)


def molecular_radius(mr: Q) -> Value:
    """Get the Stokes-Einstein radius of a molecule (m) from its aqueous diffusion."""
    d = aqueous_diffusion(mr).value
    return Value(K_B * T_BODY / (6.0 * math.pi * ETA * d), U.meter)


def cytosolic_diffusion(d_aq: Value, f_cyt: float = F_CYT) -> Value:
    """Get the diffusion in the cytosol as a fraction of the aqueous one (m**2/s)."""
    return Value(f_cyt * _as_m2_per_s(d_aq), M2_PER_S)


def axial_velocity(length: float, lumen: float, flow: float) -> float:
    """Get the velocity of the blood along a lumen (m/s).

    The length over the residence time, `L / (lumen / flow)`; the one
    expression the engine and the zonated liver of the PBPK library share,
    each with its own lumen.

    Args:
        length: the length of the lumen (m).
        lumen: the volume of the lumen (l).
        flow: the flow through the lumen (l/min).

    Returns:
        The velocity (m/s); 0 without a flow.

    Raises:
        ValueError: for a lumen of 0, which has no velocity.
    """
    if flow == 0.0:
        return 0.0
    if lumen <= 0.0:
        raise ValueError(f"The lumen {lumen} l has no volume, so it has no velocity.")
    return length / (lumen / flow * 60.0)


def dispersion(velocity: float, length: float, n_zones: int) -> float:
    """Get the numerical dispersion of an upwind chain of zones (m**2/s).

    `v * dx / 2` with the length `dx = L / n` of a zone: the dispersion the
    upwind discretization of the convection adds to the physical diffusion.

    Args:
        velocity: the velocity of the blood (m/s).
        length: the length of the chain (m).
        n_zones: the number of zones.

    Returns:
        The dispersion (m**2/s).
    """
    return velocity * length / n_zones / 2.0


def peclet(velocity: float, length: float, n_zones: int, d_eff: float) -> float:
    """Get the Peclet number an upwind chain of zones realizes.

    `v * L / (d_eff + v * dx / 2)`, which is `2 n` when the physical
    diffusion is negligible against the dispersion of the chain.

    Args:
        velocity: the velocity of the blood (m/s).
        length: the length of the chain (m).
        n_zones: the number of zones.
        d_eff: the effective axial diffusion (m**2/s).

    Returns:
        The Peclet number; 0 without a flow.
    """
    if velocity == 0.0:
        return 0.0
    return velocity * length / (d_eff + dispersion(velocity, length, n_zones))


def renkin(lambda_: float) -> float:
    """Get the Renkin hindrance of a spherical solute in a cylindrical pore.

    `(1 - l)**2 * (1 - 2.104 l + 2.09 l**3 - 0.95 l**5)` for the ratio `l` of
    the radius of the solute to that of the pore; 1 for a point and 0 once the
    solute no longer fits (`l >= 1`). The polynomial stays positive on [0, 1]
    (0.036 at 1), so the hindrance is never negative; the clamp is a guard.

    Args:
        lambda_: the ratio of the radius of the solute to that of the pore.

    Returns:
        The hindrance in [0, 1].

    Raises:
        ValueError: for a negative ratio, which is no radius.
    """
    if lambda_ < 0.0:
        raise ValueError(f"The radius ratio {lambda_} is negative.")
    if lambda_ >= 1.0:
        return 0.0
    poly = 1.0 - 2.104 * lambda_ + 2.09 * lambda_**3 - 0.95 * lambda_**5
    return max(0.0, (1.0 - lambda_) ** 2 * poly)


def fenestral_permeability(
    d_aq: Value, radius: Value, geo: SinusoidGeometry, porosity: float
) -> Value:
    """Get the permeability of the endothelium through its fenestrae (m/s).

    `porosity * d_aq / t_endo * renkin(radius / r_fen)`.

    Args:
        d_aq: the aqueous diffusion coefficient (an area per time).
        radius: the radius of the molecule (a length, converted to m).
        geo: the geometry with the thickness of the endothelium and the radius
            of a fenestra.
        porosity: the fenestral porosity (fraction).

    Returns:
        The permeability (m/s).

    Raises:
        ValueError: for a radius which is no length.
    """
    try:
        r = radius.value * LENGTH_UNITS[unit_text(radius.unit)]
    except KeyError:
        raise ValueError(
            f"A radius needs a length, got {unit_text(radius.unit)}."
        ) from None
    sieve = renkin(r / geo.r_fen)
    return Value(porosity * _as_m2_per_s(d_aq) / geo.t_endo * sieve, M_PER_S)


def _as_m2_per_s(d: Value) -> float:
    """The number of a diffusion coefficient in m**2/s."""
    try:
        return d.value * DIFFUSION_UNITS[unit_text(d.unit)]
    except KeyError:
        raise ValueError(
            f"A diffusion coefficient needs an area per time, got {unit_text(d.unit)}."
        ) from None


@runtime_checkable
class Solute(Protocol):
    """What rides the chain of zones: an id and a name, nothing more is needed.

    A perfusate quantity (oxygen, glucose, an enzyme) is no PBPK substance and
    has no molar mass; the engine carries it through the zones by convection.
    The optional attributes add the transport which needs them: a `diffusion`
    (or a molar mass `mr` to derive it from) adds the axial diffusion, a molar
    mass the path through the fenestrae (the radius of the molecule). A
    `pkpdsim.pbpk.substance.Substance` satisfies it.

    Attributes:
        sid: the id of the solute.
        name: the name of the solute.
        mr: the molar mass, `None` for none.
        diffusion: the aqueous diffusion coefficient, `None` derives it from
            the molar mass.
    """

    @property
    def sid(self) -> str:
        """The id."""
        ...

    @property
    def name(self) -> str:
        """The name."""
        ...

    @property
    def mr(self) -> Q | None:
        """The molar mass."""
        ...

    @property
    def diffusion(self) -> Value | None:
        """The aqueous diffusion coefficient."""
        ...


def diffuses(substance: Solute) -> bool:
    """Get whether a solute has a diffusion: a given one or a molar mass to derive it."""
    return substance.diffusion is not None or substance.mr is not None


def diffusion_of(substance: Solute) -> Value:
    """Get the aqueous diffusion of a substance (m**2/s), the override or derived.

    Raises:
        ValueError: for a solute with neither a diffusion nor a molar mass.
    """
    if substance.diffusion is not None:
        return Value(
            _as_m2_per_s(substance.diffusion),
            M2_PER_S,
            substance.diffusion.source,
            substance.diffusion.notes,
        )
    if substance.mr is None:
        raise ValueError(
            f"The solute '{substance.sid}' has neither a diffusion nor a molar mass."
        )
    return aqueous_diffusion(substance.mr)


def diffusion_parameter(substance: Solute, note: str = "") -> Parameter:
    """Get the diffusion of a substance as a parameter `D_<sid>` with its provenance.

    Args:
        substance: the substance.
        note: a paragraph which is appended to the notes, what the model of
            the caller does with the coefficient (the zonated liver of the
            PBPK library states the Peclet number of its chain).

    Returns:
        The parameter.
    """
    d = diffusion_of(substance)
    if substance.diffusion is not None:
        text = (
            f"The diffusion coefficient of {substance.name} in water, given "
            f"(not derived) as {substance.diffusion.value:g} "
            f"{unit_text(substance.diffusion.unit)}."
        )
    else:
        assert substance.mr is not None  # diffusion_of raised otherwise
        text = (
            f"The diffusion coefficient of {substance.name} in water at 37 C, "
            f"derived from its molar mass {_mr_number(substance.mr):g} g/mole: "
            f"{CORRELATION}."
        )
    return value_parameter(
        ids.diffusion(substance.sid),
        Value(d.value, M2_PER_S, d.source, d.notes),
        M2_PER_S,
        f"diffusion {substance.name}",
        SBO.DIFFUSION_COEFFICIENT,
        f"{text}\n\n{note}" if note else text,
    )
