"""The flow through the sinusoids: imposed, or driven by a pressure drop.

A zonated liver has the blood flow of its sinusoids either **imposed** (a
number, `ImposedFlow`: the flow of the liver, or the flow per g of liver as the
machine perfusion states it) or **driven by the pressure** between the portal
and the central end of the sinusoid (`PressureFlow`).

The pressure-driven sinusoid is the leaky tube of the spt model: a tube of the
radius `r_sin` with the specific hydraulic resistance `W = 8 nu / (pi r_sin**4)`
(Pa s/m**4) of the lumen, whose wall is crossed by the fenestrae with the
specific hydraulic resistance `w` (Pa s/m**2) of all pores in parallel. The
pressure `P(x)` along the tube follows `P'' = (P - P0) / lambda**2` with
`lambda = sqrt(w / W)` and the pressures `P(0) = p_pp`, `P(L) = p_pc` at its
ends, so the flow `Q(x) = -P'(x) / W` is, with `l = L / lambda` and `s = x /
lambda`,

    Q(x) = (p_pp' cosh(l - s) - p_pc' cosh(s)) / (W lambda sinh(l)),

`p' = p - P0`. This is the exponential solution of spt (`PP_Q` of
`spt_liver_zonated.py`) written with hyperbolic functions, which stays finite
when the leak vanishes (`lambda` diverges, and the flow tends to the
Poiseuille flow `(p_pp - p_pc) / (W L)` of a tube without a leak).

**The flow depends on the position, and the chain evaluates it at the faces of
its zones.** What leaks through the fenestrae leaves the lumen, so a
pressure-driven sinusoid has a position-dependent flow and the convection of
the chain no longer carries one `Q` through all zones. A flow crosses a face,
not a volume, so the engine is a cell-centred finite volume chain: the `n + 1`
faces of the zones (`profiles.face_positions`) each have the flow
`Q_sin_<site>` (m**3/s of one sinusoid) and `Q_li_<site>` (l/min, the liver has
`N_sin` sinusoids, which follows the liver mass), where the site is `in` for
the periportal face of the chain (`x = 0`) and the zone for its pericentral
face (`x = k L / n`, `ids.zone_flow`). The inflow of the chain reads the flow
of the periportal face, the convection out of a zone the flow of the face it
crosses, and the outflow the flow at `x = L`, so the organ receives the flow at
the portal face which the calibration targets, at every zone count. The
volumes and the zonal profiles stay at the centres of the zones
(`profiles.zone_positions`), and a zone is read at its centre, but the upwind
scheme gives it the value of its pericentral face wherever the flow varies
along the chain: the steady concentration of a solute the fenestrae stop is
`Q(0) / Q(x_k^face)`, +4.9 % in the first zone to -5.0 % in the last against
the average over the zone at ten zones. That bias is the first order error of
the scheme, which the zone count lowers and a higher order reconstruction
would remove (`docs/zonation.md`, "Numerical dispersion"). spt evaluated the flow at `x = 0` only (its
`FLOW: Q dependent on position`); this finishes it. The chain does not move the
leaking fluid out of the lumen: the convection carries the solute with the flow
of the face from zone to zone, and the leak of the solute is the fenestral
exchange of the radial barrier.

**The count of the fenestrae is derived from the porosity.** `w` needs the
number of fenestrae per area, `N_fen`. The table gives `n_fen = 1e13 /m**2`
which with the radius `r_fen = 53.5 nm` is a porosity of 8.99 %, while the
measured porosity of the same table (Braet and Wisse 2002) is 6 % periportal
and 8 % pericentral and the sieving of the barrier uses it. Using both would
give one endothelium two porosities, so the count is `N_fen = eps / (pi
r_fen**2)` (6.7e12 periportal, 8.9e12 pericentral), and the table value is a
cross-check (`tests/zonation/test_flow.py`, within 35 %). The leaky tube has
one `w`, so it is evaluated at the mean porosity of the sinusoid, the mean of
the zones, which is exact for a tube of uniform porosity and the first order
for the linear gradient of the table. The analytic solution assumes a uniform
leak, so the hydraulic resistance uses the mean porosity and the per-zone
`fenestra_counts` enter only the cross-check; the zonal `porosity` profile
scales the solute path of the barrier and not the hydraulic resistance.

**The oncotic pressure** `P0` (`p_onc`, `Poc - Pot` of spt) defaults to the mean
of the two pressures, as spt sets `P0 = 0.5 (Pa + Pb)`. The net filtration is
zero then: what leaks out of the first half of the tube returns in the second
half, and the flow is symmetric about the middle. A default which filtered
would change the mass balance of every model which switches to the pressure
flow, so a lower oncotic pressure is a choice of the model.

**The flow is not monotonic.** With `P0` at the mean of the two pressures the
flow falls from the portal end to a minimum in the middle and rises again to
the central end, though the pressure falls monotonically: fluid filters out of
the first half of the tube and is reabsorbed in the second. This is the
signature of `P0` and no defect.

**The resistance factor is calibrated.** `nu = resistance_factor * viscosity`
and the default factor is derived per species so that the periportal inflow of
the **reference** liver at the default pressures is the reference flow of the
settings (`ZonationSettings.reference_flow`, which the PBPK library fills with
the total liver flow of its reference individual,
`pkpdsim.pbpk.physiology.reference_flow`), the flow an imposed flow
takes by default: switching between the two modes does not change the flow the
organ receives. The chain evaluates its flows at the faces of its zones, so
the flow it receives is the one at the portal face the factor is calibrated to
and the guarantee holds at every zone count, which
`test_the_chain_receives_the_reference_flow_in_both_modes` pins to 1e-12.
Both modes follow the liver mass from there - the pressure
flow because the number of sinusoids in parallel does
(`sinusoid_count`), the imposed one because it scales the reference flow with
`mass / reference_liver_mass` (`imposed_flow`) - so they agree at every liver
mass and not only at the reference one. The 10 of the spt model gives about
five times that flow for the rat. The guarantee holds at the default pressures,
`p_onc`, viscosity and the table geometry; the factor is calibrated at those
and does not follow a change of them.

**The inflow is calibrated, the mean flow is lower.** The target is the
periportal inflow, what the portal vein and the hepatic artery deliver. `Q_li`,
the mean over the zones, sets the velocity and the Peclet number and is not the
reference flow: with `p_onc` at the mean pressure the fluid filters out of the
first half and is reabsorbed in the second, so the mean throughput is 17 %
below the inflow at 10 zones (`Q_li` is 0.8315 of the reference flow, 0.8298 at
40 zones and 0.8297 in the continuum), and the residence time is
correspondingly longer than `V / QH`. It is the mean over the zones of the mean
of the two faces of a zone, the trapezoidal rule over the faces (`mean_flow`),
which is the inflow itself for a single zone, whose two faces are the two ends
of the chain and carry the same flow.

A pressure which gives a negative flow in a zone is refused: the chain does
not run backwards.
"""

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

from sbmlutils.factory import AssignmentRule, InitialAssignment, Parameter
from sbmlutils.factory.units import UnitDefinition
from sbmlutils.metadata import SBO

from pkpdsim.sbml.annotation import initial_assignment
from pkpdsim.sbml.parameters import value_parameter
from pkpdsim.sbml.units import U
from pkpdsim.sbml.values import (
    Q,
    Value,
    derived,
    model_choice,
    number,
    unit_text,
)
from pkpdsim.species import Species
from pkpdsim.zonation import ids
from pkpdsim.zonation.geometry import (
    Layer,
    SinusoidGeometry,
    lumen_fraction,
    mass_kg,
    porosities,
    reference_liver_mass,
    sinusoid_count,
)
from pkpdsim.zonation.profiles import face_positions

if TYPE_CHECKING:
    from pkpdsim.zonation.settings import ZonationSettings

__all__ = [
    "ImposedFlow",
    "LeakyTube",
    "PressureFlow",
    "calibrated_resistance_factor",
    "face_flows",
    "fenestra_counts",
    "flow_objects",
    "imposed_flow",
    "leaky_tube",
    "mean_flow",
]

PA = UnitDefinition("Pa", "pascal", name="pascal")
"""The unit of a pressure."""

PA_S = UnitDefinition("Pa_s", "pascal*second", name="pascal second")
"""The unit of a viscosity."""

PA_S_PER_M4 = UnitDefinition(
    "Pa_s_per_m4", "pascal*second/meter^4", name="pascal second per metre^4"
)
"""The unit of the specific hydraulic resistance of a tube."""

PA_S_PER_M2 = UnitDefinition(
    "Pa_s_per_m2", "pascal*second/meter^2", name="pascal second per metre^2"
)
"""The unit of the specific hydraulic resistance of the pores."""

M3_PER_S = UnitDefinition("m3_per_s", "meter^3/second", name="cubic metre per second")
"""The unit of the flow of one sinusoid."""

PER_M2 = UnitDefinition("per_m2", "1/meter^2", name="per square metre")
"""The unit of the count of fenestrae per area."""

_P = SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER
_FLOW_TO_L_PER_MIN = {"l/min": 1.0, "ml/min": 1e-3}
_PER_MASS_TO_L_PER_MIN_KG = {"ml/min/g": 1.0, "l/min/kg": 1.0, "ml/min/kg": 1e-3}
_PRESSURE_TO_PA = {"Pa": 1.0, "kPa": 1e3, "mmHg": 133.322}
_VISCOSITY_TO_PA_S = {"Pa*s": 1.0, "Pa s": 1.0, "mPa*s": 1e-3, "mPa s": 1e-3}


_PORTAL_PRESSURE = Value(1333.22, "Pa", notes="10 mmHg, the portal pressure of spt")
_CENTRAL_PRESSURE = Value(266.64, "Pa", notes="2 mmHg, the perivenous pressure of spt")
_PLASMA_VISCOSITY = Value(0.0018, "Pa*s", notes="the plasma viscosity of spt")

_NO_REFERENCE = (
    "The flow of the organ is its reference flow, and the settings state none. "
    "Give the number: `ZonationSettings(..., reference_flow=Value(...))`, the "
    "blood flow of the whole reference organ in l/min or ml/min; for a species of "
    "the PBPK library it is `pkpdsim.pbpk.physiology.reference_flow(species)`, "
    "which returns that `Value`. The engine reads no physiology of its own. Or "
    "state the flow itself instead of a reference: `ImposedFlow(value=...)` or "
    "`ImposedFlow(per_mass=...)` for a measured flow, "
    "`PressureFlow(resistance_factor=...)` for a flow driven by the pressures."
)
"""The message of a flow which takes a reference flow the settings do not state."""


@dataclass(frozen=True)
class ImposedFlow:
    """A volumetric flow imposed on the organ.

    Attributes:
        value: the flow through the whole liver, in a volume per time
            (l/min or ml/min); a plain number is l/min.
        per_mass: the flow per mass of liver (ml/min/g, l/min/kg or
            ml/min/kg), which the flow of the liver follows, as the machine
            perfusion states it. `value` and `per_mass` exclude each other.
            Giving neither takes the reference flow of the settings
            (`ZonationSettings.reference_flow`, for a species of the PBPK
            library `pkpdsim.pbpk.physiology.reference_flow(species)`), and is
            an error when the settings state none.
    """

    value: Q | None = None
    per_mass: Q | None = None

    def __post_init__(self) -> None:
        """Check that the flow is given once.

        Raises:
            ValueError: for a `value` and a `per_mass` together.
        """
        if self.value is not None and self.per_mass is not None:
            raise ValueError(
                "An imposed flow is the flow of the liver (value) or a flow per "
                "mass of liver (per_mass), not both."
            )


@dataclass(frozen=True)
class PressureFlow:
    """A flow driven by the pressure drop over the sinusoid.

    Attributes:
        p_pp: the pressure at the periportal end; the default 1333.22 Pa is
            10 mmHg (the portal pressure of the spt model).
        p_pc: the pressure at the pericentral end; the default 266.64 Pa is
            2 mmHg (the perivenous pressure of the spt model).
        viscosity: the viscosity of the plasma.
        resistance_factor: the factor of the viscosity of the blood in the
            sinusoid over the plasma viscosity, `nu = resistance_factor *
            viscosity`. `None` (the default) is calibrated per species
            (`calibrated_resistance_factor`) so that the periportal inflow at
            the default pressures is the reference flow of the settings
            (`ZonationSettings.reference_flow`, for a species of the PBPK
            library `pkpdsim.pbpk.physiology.reference_flow(species)`), which
            is what an `ImposedFlow` without a value gives, and is an error
            when the settings state none. The 10 of the spt
            model is an uncalibrated stand-in for the microvascular geometry
            and gives about five times the reference flow for the rat. The
            agreement with an imposed flow holds only with the default
            `p_onc`, `viscosity` and pressures: an explicit one breaks it.
        p_onc: the oncotic pressure `P0 = Poc - Pot` the pressure relaxes to
            along the leaky tube; `None` (the default) is the mean of the two
            pressures, as spt sets it, which makes the net filtration zero.
    """

    p_pp: Q = _PORTAL_PRESSURE
    p_pc: Q = _CENTRAL_PRESSURE
    viscosity: Q = _PLASMA_VISCOSITY
    resistance_factor: float | None = None
    p_onc: Q | None = None


def _converted(q: Q, table: dict[str, float], what: str) -> float:
    """The number of `q` in the base unit of `table`, a plain number is in the base unit."""
    x = number(q)
    if x is None:
        raise ValueError(f"The {what} needs a number, not a formula.")
    if not isinstance(q, Value):
        return x
    unit = unit_text(q.unit)
    try:
        return x * table[unit]
    except KeyError:
        raise ValueError(
            f"The {what} needs a unit of {sorted(table)}, got {unit}."
        ) from None


def imposed_flow(
    flow: ImposedFlow, reference: Value | None, mass: float, reference_mass: float
) -> float:
    """Get an imposed flow of the liver in l/min.

    A flow which takes the reference flow of the settings is scaled with the
    liver mass, `reference * mass / reference_mass`: the reference flow is
    the flow of the reference liver of the species, and a liver of another
    mass receives its share of it, as a flow per mass does and as the flow of
    a `PressureFlow` does (its sinusoid count follows the mass). The two
    modes therefore give the same flow at every liver mass and not only at
    the reference one.

    Args:
        flow: the imposed flow.
        reference: the reference flow of the organ
            (`ZonationSettings.reference_flow`), which a flow without a value
            and without a flow per mass takes; `None` for none.
        mass: the mass of the liver (kg).
        reference_mass: the mass of the reference liver of the species (kg,
            `geometry.reference_liver_mass`), which the reference flow belongs
            to.

    Returns:
        The flow (l/min), as given, not checked for its sign.

    Raises:
        ValueError: for a flow which is a formula, has no flow unit, a flow
            per mass without a unit, or a flow which takes the reference flow
            of settings which state none.
    """
    if flow.per_mass is not None:
        if not isinstance(flow.per_mass, Value):
            raise ValueError(
                "A flow per mass needs a unit (ml/min/g, l/min/kg or ml/min/kg)."
            )
        return mass * _converted(
            flow.per_mass, _PER_MASS_TO_L_PER_MIN_KG, "flow per mass"
        )
    if flow.value is None:
        if reference is None:
            raise ValueError(_NO_REFERENCE)
        q = _converted(reference, _FLOW_TO_L_PER_MIN, "reference flow")
        return q * mass / reference_mass
    if number(flow.value) is None:
        raise ValueError(
            "The imposed flow needs a number, not a formula: a flow computed from "
            "the pressures is a PressureFlow."
        )
    return _converted(flow.value, _FLOW_TO_L_PER_MIN, "flow")


@dataclass(frozen=True)
class LeakyTube:
    """The leaky tube of a sinusoid, in SI units.

    Attributes:
        r_sin: the radius of the lumen (m).
        length: the length (m).
        t_endo: the thickness of the endothelium (m).
        r_fen: the radius of a fenestra (m).
        n_fen: the number of fenestrae per endothelial area (1/m**2); 0 is a
            tube without a leak.
        viscosity: the viscosity of the blood in the sinusoid (Pa s).
        p_pp: the periportal pressure (Pa).
        p_pc: the pericentral pressure (Pa).
        p_onc: the oncotic pressure (Pa).
    """

    r_sin: float
    length: float
    t_endo: float
    r_fen: float
    n_fen: float
    viscosity: float
    p_pp: float
    p_pc: float
    p_onc: float

    def resistance_lumen(self) -> float:
        """Get the specific hydraulic resistance `W` of the lumen (Pa s/m**4)."""
        return 8.0 * self.viscosity / (math.pi * self.r_sin**4)

    def resistance_pores(self) -> float:
        """Get the specific hydraulic resistance `w` of all pores (Pa s/m**2).

        `4 nu t_endo / (pi**2 r_fen**4 r_sin N_fen)`, infinite without
        fenestrae.
        """
        if self.n_fen == 0.0:
            return math.inf
        return (
            4.0
            * self.viscosity
            * self.t_endo
            / (math.pi**2 * self.r_fen**4 * self.r_sin * self.n_fen)
        )

    def decay_length(self) -> float:
        """Get `lambda = sqrt(w / W)` (m), infinite without fenestrae."""
        return math.sqrt(self.resistance_pores() / self.resistance_lumen())

    def flow(self, position: float) -> float:
        """Get the flow (m**3/s) of one sinusoid at a position.

        Args:
            position: the position along the tube as a fraction of its
                length, 0 periportal.

        Returns:
            The flow, the Poiseuille flow `(p_pp - p_pc) / (W L)` of a tube
            without fenestrae.
        """
        w_cap = self.resistance_lumen()
        if self.n_fen == 0.0:
            return (self.p_pp - self.p_pc) / (w_cap * self.length)
        lam = self.decay_length()
        big, s = self.length / lam, position * self.length / lam
        a, b = self.p_pp - self.p_onc, self.p_pc - self.p_onc
        return (a * math.cosh(big - s) - b * math.cosh(s)) / (
            w_cap * lam * math.sinh(big)
        )


def fenestra_counts(geo: SinusoidGeometry, n_zones: int) -> tuple[float, ...]:
    """Get the count of fenestrae per endothelial area of every zone (1/m**2).

    The porosity of the zone over the area of a fenestra, `eps / (pi r_fen**2)`:
    the porosity of the table, interpolated over the positions of the zones,
    is the one number of the endothelium (see the module docstring). The
    gradient itself is `geometry.porosities`, the one the barrier reads, so
    the hydraulics and the sieving cannot drift apart.

    Args:
        geo: the geometry.
        n_zones: the number of zones.

    Returns:
        The counts from periportal to pericentral.
    """
    area = math.pi * geo.r_fen.value**2
    return tuple(eps.value / area for eps in porosities(n_zones, geo))


def _mean_porosity(geo: SinusoidGeometry) -> float:
    """The porosity of the leaky tube, the mean of the periportal and pericentral porosity."""
    return 0.5 * (geo.porosity_pp.value + geo.porosity_pc.value)


def calibrated_resistance_factor(
    species: Species,
    geo: SinusoidGeometry,
    reference: Value,
    fr: Mapping[Layer, Value],
) -> float:
    """Get the resistance factor for which the inflow is the reference flow.

    The flow of the leaky tube is inversely proportional to the viscosity (`W`
    and `w` both scale with it, `lambda` does not), so the factor is the
    periportal inflow of the whole reference liver at the default pressures
    and the plasma viscosity (factor 1) over the reference flow. The
    entrance flow is what the organ receives; the flows of the zones do not
    add up in a leaky tube. The factor is a property of the species, its
    geometry and its volume fractions, which give the number of sinusoids the
    inflow is spread over: it does not follow the pressures, the oncotic
    pressure or the viscosity of a model, so the agreement with an imposed flow
    holds at their defaults only. The fractions are those of the build
    (`ZonationSettings.full_fractions`) and not the table of the species, or
    the calibration would contradict the geometry it calibrates for a caller
    which states its own fractions.

    The inflow is calibrated, not the mean over the zones: `Q_li`, which sets
    the velocity and the Peclet number, is 0.8315 of the reference flow at 10
    zones and 0.8298 at 40, because with `P0` at the mean pressure fluid
    filters out of the first half of the tube and returns in the second, so
    the residence time is longer than `V / QH`. What the chain receives is
    this inflow, the flow of the periportal face of its first zone, at every
    zone count.

    Args:
        species: the species.
        geo: the geometry.
        reference: the reference flow of the organ
            (`ZonationSettings.reference_flow`), the target of the calibration.
        fr: the volume fractions of the layers of the build
            (`ZonationSettings.full_fractions`), which give the lumen the
            sinusoids carry.

    Returns:
        The factor (about 50 for the rat).
    """
    tube = leaky_tube(PressureFlow(resistance_factor=1.0), geo, species, reference, fr)
    n_sin = sinusoid_count(reference_liver_mass(species), geo, fr)
    inflow = tube.flow(0.0) * n_sin * 1e3 * 60.0
    return inflow / _converted(reference, _FLOW_TO_L_PER_MIN, "reference flow")


def resolved_resistance_factor(
    flow: PressureFlow,
    geo: SinusoidGeometry,
    species: Species,
    reference: Value | None,
    fr: Mapping[Layer, Value],
) -> float:
    """Get the resistance factor of a pressure flow: the given or the calibrated one.

    Args:
        flow: the pressure flow.
        geo: the geometry.
        species: the species.
        reference: the reference flow of the organ, the target of a
            calibration.
        fr: the volume fractions of the layers of the build, which a
            calibration needs.

    Returns:
        The factor.

    Raises:
        ValueError: for a factor which is calibrated and settings which state
            no reference flow.
    """
    if flow.resistance_factor is not None:
        return flow.resistance_factor
    if reference is None:
        raise ValueError(_NO_REFERENCE)
    return calibrated_resistance_factor(species, geo, reference, fr)


def leaky_tube(
    flow: PressureFlow,
    geo: SinusoidGeometry,
    species: Species,
    reference: Value | None,
    fr: Mapping[Layer, Value],
) -> LeakyTube:
    """Get the leaky tube of a pressure flow and a geometry.

    Args:
        flow: the pressure flow.
        geo: the geometry.
        species: the species, which calibrates a resistance factor of `None`.
        reference: the reference flow of the organ
            (`ZonationSettings.reference_flow`), the target of the calibration.
        fr: the volume fractions of the layers of the build, which a
            calibration needs.

    Returns:
        The tube, at the mean porosity of the sinusoid.
    """
    p_pp = _converted(flow.p_pp, _PRESSURE_TO_PA, "periportal pressure")
    p_pc = _converted(flow.p_pc, _PRESSURE_TO_PA, "pericentral pressure")
    p_onc = (
        0.5 * (p_pp + p_pc)
        if flow.p_onc is None
        else _converted(flow.p_onc, _PRESSURE_TO_PA, "oncotic pressure")
    )
    nu = resolved_resistance_factor(flow, geo, species, reference, fr) * _converted(
        flow.viscosity, _VISCOSITY_TO_PA_S, "viscosity"
    )
    return LeakyTube(
        r_sin=geo.r_sin.value,
        length=geo.L_sin.value,
        t_endo=geo.t_endo.value,
        r_fen=geo.r_fen.value,
        n_fen=_mean_porosity(geo) / (math.pi * geo.r_fen.value**2),
        viscosity=nu,
        p_pp=p_pp,
        p_pc=p_pc,
        p_onc=p_onc,
    )


def face_flows(
    settings: ZonationSettings, geo: SinusoidGeometry, mass: Q
) -> tuple[float, ...]:
    """Get the flow of the liver at every face of the zones of a pressure flow (l/min).

    The `n + 1` faces of `n` zones (`profiles.face_positions`), periportal
    first: the chain carries the flow of a face across it, so the first entry
    is the flow the organ receives (`x = 0`, the target of the calibration)
    and the last the flow it releases (`x = L`).

    Args:
        settings: the settings, with a `PressureFlow`.
        geo: the geometry.
        mass: the liver mass.

    Returns:
        The flow of the sinusoids of the liver at every face, `n + 1` values.

    Raises:
        ValueError: for a pressure drop which gives a negative flow.
    """
    flow = settings.flow
    if not isinstance(flow, PressureFlow):
        raise TypeError("Only a pressure flow has a flow per face.")
    fr = settings.full_fractions()
    tube = leaky_tube(flow, geo, settings.species, settings.reference_flow, fr)
    n_sin = sinusoid_count(mass, geo, fr)
    out = tuple(
        tube.flow(x) * n_sin * 1e3 * 60.0 for x in face_positions(settings.n_zones)
    )
    if min(out) < 0.0:
        raise ValueError(
            f"The pressures {tube.p_pp} Pa (periportal) and {tube.p_pc} Pa "
            f"(pericentral) give a negative flow ({min(out):.3g} l/min): the "
            "blood would run from the pericentral to the periportal end, which "
            "is not a flow of the chain."
        )
    return out


def mean_flow(flows: Sequence[float]) -> float:
    """Get the mean flow of the chain over its zones from the flows of its faces.

    The mean over the zones of the mean of the two faces of a zone, which is
    the trapezoidal rule over the faces: `(Q(0)/2 + Q(1/n) + ... + Q(1 -
    1/n) + Q(1)/2) / n`. This is `Q_li` of a pressure flow, the flow which
    sets the velocity and the Peclet number of the chain, and it is below the
    inflow `Q(0)` at the default oncotic pressure (0.8315 of it at 10 zones,
    0.8298 at 40), and equal to it for a single zone, whose two faces are the
    two ends of the chain.

    Args:
        flows: the flows of the `n + 1` faces (`face_flows`), periportal
            first.

    Returns:
        The mean flow, in the unit of `flows`.
    """
    n = len(flows) - 1
    return (0.5 * (flows[0] + flows[-1]) + sum(flows[1:-1])) / n


def _constant(pid: str, q: Q, unit: UnitDefinition, name: str) -> Parameter:
    """A constant parameter of a SI number, with the provenance of a `Value`.

    A number of the geometry is the `Value` of its row of the table, whose
    source and notes reach the notes of the parameter; a plain number carries
    its unit alone, as the pressures and the viscosity of the leaky tube do.
    """
    return value_parameter(pid, q, unit, name, _P)


def _derived(pid: str, q: Q, unit: UnitDefinition | object, name: str) -> Parameter:
    """A constant parameter whose value an initial assignment sets again.

    It carries the provenance of a `Value` as a constant does: the assignment
    recomputes the number, not where it comes from.
    """
    if isinstance(q, Value):
        return value_parameter(pid, q, unit, name, _P)  # ty: ignore[invalid-argument-type]
    return Parameter(pid, q, unit=unit, name=name, constant=True, sboTerm=_P)  # ty: ignore[invalid-argument-type]


def _ruled(
    pid: str, value: float, unit: UnitDefinition | object, name: str
) -> Parameter:
    """A parameter of the flow chain, which an assignment rule computes at any time."""
    return Parameter(pid, value, unit=unit, name=name, constant=False, sboTerm=_P)  # ty: ignore[invalid-argument-type]


def assignment_rule(
    pid: str, formula: str, unit: UnitDefinition | object, name: str
) -> AssignmentRule:
    """Get the assignment rule of a quantity of the flow chain.

    The flow chain is computed by rules and not by initial assignments: a
    simulator applies a change after the initialization, so a change of the
    organ flow (or of a pressure, of the reference of the bile flow) reaches
    every flow which derives from it only through a rule.

    Args:
        pid: the id of the target.
        formula: the formula.
        unit: the unit of the formula.
        name: the name of the target.

    Returns:
        The assignment rule.
    """
    return AssignmentRule(
        pid,
        formula,
        unit,  # ty: ignore[invalid-argument-type]
        name=name,
        sboTerm=SBO.MATHEMATICAL_EXPRESSION,
    )


def flow_objects(
    settings: ZonationSettings, geo: SinusoidGeometry, mass: Q
) -> list[Parameter | InitialAssignment | AssignmentRule]:
    """Get the parameters and initial assignments of the flow of the liver.

    An imposed flow is the parameter `Q_li` (l/min). A pressure flow is the
    leaky tube: the pressures, the viscosity, the geometry, the resistances
    `W_sin` and `w_sin`, `lambda_sin`, the flow of one sinusoid `Q_sin_<zone>`
    (m**3/s) and of the liver `Q_li_<zone>` (l/min) at every zone, and `Q_li`,
    the mean of the zones (the flow which sets the velocity of the chain).
    The flow follows the pressures and the liver mass (the initial
    assignments are evaluated again when they change before a simulation).

    Args:
        settings: the settings.
        geo: the geometry.
        mass: the liver mass.

    Returns:
        The parameters, initial assignments and assignment rules.

    Raises:
        ValueError: for a negative flow, or a flow which is a formula.
    """
    flow = settings.flow
    if isinstance(flow, ImposedFlow):
        reference_mass = reference_liver_mass(settings.species)
        q = imposed_flow(
            flow, settings.reference_flow, mass_kg(mass), reference_mass.value
        )
        if q < 0.0:
            raise ValueError(
                f"The imposed flow {q} l/min is negative: it would run the "
                "blood from the pericentral to the periportal end, which is "
                "not a flow of the chain."
            )
        return [
            value_parameter(
                ids.FLOW,
                _imposed_flow_value(q, flow, settings.reference_flow, reference_mass),
                U.l_per_min,
                "blood flow through the liver",
                _P,
            )
        ]
    return _pressure_objects(settings, flow, geo, mass)


def _imposed_flow_value(
    q: float, flow: ImposedFlow, reference: Value | None, reference_mass: Value
) -> Value:
    """The imposed flow of the liver with its provenance.

    A flow which takes the reference flow of the settings rests on that flow
    and on the mass of the reference liver, both sourced, so it is a derived
    value naming them; a flow or a flow per mass the caller states is the
    model choice of the caller.
    """
    if flow.per_mass is not None:
        return model_choice(
            q,
            U.l_per_min,
            "the flow per g of the liver which the caller of the engine "
            "imposed, times the mass of the liver",
        )
    if flow.value is not None:
        return model_choice(
            q,
            U.l_per_min,
            "the blood flow of the liver which the caller of the engine imposed",
        )
    assert reference is not None  # imposed_flow raises without it
    return derived(
        q,
        U.l_per_min,
        (reference, reference_mass),
        "the reference flow of the whole organ scaled with the mass of the "
        "liver, `reference * mass / reference_mass`, so that the imposed and "
        "the pressure-driven mode agree at every liver mass. Reference flow: "
        f"{reference.notes}. Reference liver mass ({reference_mass.value:g} "
        f"kg): {reference_mass.notes}",
    )


def _pressure_objects(
    settings: ZonationSettings,
    flow: PressureFlow,
    geo: SinusoidGeometry,
    mass: Q,
) -> list[Parameter | InitialAssignment | AssignmentRule]:
    """The leaky tube: parameters and initial assignments, the flows as rules."""
    fr = settings.full_fractions()
    tube = leaky_tube(flow, geo, settings.species, settings.reference_flow, fr)
    factor = resolved_resistance_factor(
        flow, geo, settings.species, settings.reference_flow, fr
    )
    flows = face_flows(settings, geo, mass)
    n = settings.n_zones
    f_lumen = derived(
        lumen_fraction(fr),
        U.dimensionless,
        [fr[Layer.PLASMA], *([fr[Layer.RBC]] if Layer.RBC in fr else [])],
        "the volume fraction of the lumen of the sinusoids in the liver, the "
        "plasma and the erythrocytes of it",
    )
    n_sin = sinusoid_count(mass, geo, fr)
    w_cap, w_pore = tube.resistance_lumen(), tube.resistance_pores()
    flat = tube.n_fen == 0.0
    params: list[Parameter | InitialAssignment] = [
        _constant(ids.RADIUS, geo.r_sin, U.meter, "radius of the sinusoid"),
        _constant(ids.T_ENDO, geo.t_endo, U.meter, "thickness of the endothelium"),
        _constant(ids.LENGTH, geo.L_sin, U.meter, "length of the sinusoid"),
        _constant(ids.FENESTRA_RADIUS, geo.r_fen, U.meter, "radius of a fenestra"),
        _constant(
            ids.PRESSURE_PORTAL,
            tube.p_pp,
            PA,
            "pressure at the periportal end of the sinusoid",
        ),
        _constant(
            ids.PRESSURE_CENTRAL,
            tube.p_pc,
            PA,
            "pressure at the pericentral end of the sinusoid",
        ),
        _derived(
            ids.PRESSURE_ONCOTIC,
            tube.p_onc,
            PA,
            "oncotic pressure the pressure relaxes to along the sinusoid",
        ),
        _constant(
            ids.VISCOSITY,
            tube.viscosity / factor,
            PA_S,
            "viscosity of the plasma",
        ),
        _constant(
            ids.VISCOSITY_FACTOR,
            factor,
            U.dimensionless,
            "factor of the hepatic viscosity over the plasma viscosity",
        ),
        _derived(
            ids.VISCOSITY_LIVER,
            tube.viscosity,
            PA_S,
            "viscosity of the blood in the sinusoid",
        ),
        _constant(
            ids.POROSITY_MEAN,
            derived(
                _mean_porosity(geo),
                U.dimensionless,
                (geo.porosity_pp, geo.porosity_pc),
                "the fenestral porosity the leaky tube is evaluated at, the "
                "mean of the periportal and the pericentral porosity of the "
                f"sinusoid. Periportal: {geo.porosity_pp.notes}. Pericentral: "
                f"{geo.porosity_pc.notes}",
            ),
            U.dimensionless,
            "fenestral porosity the leaky tube is evaluated at (mean of the sinusoid)",
        ),
        _derived(
            ids.FENESTRA_COUNT,
            derived(
                tube.n_fen,
                PER_M2,
                (geo.porosity_pp, geo.porosity_pc, geo.r_fen),
                "the number of fenestrae per endothelial area, the mean "
                "porosity of the sinusoid over the area of a fenestra "
                f"(`eps / (pi r_fen**2)`); the count {geo.n_fen.value:g} /m**2 "
                "of the table is a cross-check and no input of the hydraulics "
                f"(see the module docstring). Fenestra radius: {geo.r_fen.notes}. "
                f"Count of the table: {geo.n_fen.notes}",
            ),
            PER_M2,
            "number of fenestrae per endothelial area (porosity over the area of a fenestra)",
        ),
        _constant(
            ids.LUMEN_FRACTION,
            f_lumen,
            U.dimensionless,
            "volume fraction of the lumen of the sinusoids in the liver",
        ),
        _derived(
            ids.N_SINUSOIDS, n_sin, U.dimensionless, "number of sinusoids of the liver"
        ),
        _derived(
            ids.RESISTANCE_LUMEN,
            w_cap,
            PA_S_PER_M4,
            "specific hydraulic resistance of the lumen of the sinusoid",
        ),
    ]
    assignments: list[Parameter | InitialAssignment] = []
    if flow.p_onc is None:
        assignments.append(
            initial_assignment(
                ids.PRESSURE_ONCOTIC,
                f"0.5 dimensionless * ({ids.PRESSURE_PORTAL} + {ids.PRESSURE_CENTRAL})",
                PA,
                "oncotic pressure the pressure relaxes to along the sinusoid",
            )
        )
    assignments += [
        initial_assignment(
            ids.VISCOSITY_LIVER,
            f"{ids.VISCOSITY_FACTOR} * {ids.VISCOSITY}",
            PA_S,
            "viscosity of the blood in the sinusoid",
        ),
        initial_assignment(
            ids.FENESTRA_COUNT,
            f"{ids.POROSITY_MEAN} / (pi * {ids.FENESTRA_RADIUS}^2)",
            PER_M2,
            "number of fenestrae per endothelial area",
        ),
        initial_assignment(
            ids.N_SINUSOIDS,
            f"{ids.LUMEN_FRACTION} * {ids.PARENCHYMA} * {ids.MASS} / {ids.DENSITY} "
            f"/ 1000 l_per_m3 / (pi * {ids.RADIUS}^2 * {ids.LENGTH})",
            U.dimensionless,
            "number of sinusoids of the liver",
        ),
        initial_assignment(
            ids.RESISTANCE_LUMEN,
            f"8 dimensionless * {ids.VISCOSITY_LIVER} / (pi * {ids.RADIUS}^4)",
            PA_S_PER_M4,
            "specific hydraulic resistance of the lumen of the sinusoid",
        ),
    ]
    out: list[Parameter | InitialAssignment | AssignmentRule] = [
        *params,
        *assignments,
    ]
    if not flat:
        out += [
            _derived(
                ids.RESISTANCE_PORES,
                w_pore,
                PA_S_PER_M2,
                "specific hydraulic resistance of all pores of the endothelium",
            ),
            _derived(
                ids.DECAY_LENGTH,
                tube.decay_length(),
                U.meter,
                "length over which the pressure of the leaky tube relaxes",
            ),
            initial_assignment(
                ids.RESISTANCE_PORES,
                f"4 dimensionless * {ids.VISCOSITY_LIVER} * {ids.T_ENDO} / "
                f"(pi^2 * {ids.FENESTRA_RADIUS}^4 * {ids.RADIUS} * {ids.FENESTRA_COUNT})",
                PA_S_PER_M2,
                "specific hydraulic resistance of all pores of the endothelium",
            ),
            initial_assignment(
                ids.DECAY_LENGTH,
                f"sqrt({ids.RESISTANCE_PORES} / {ids.RESISTANCE_LUMEN})",
                U.meter,
                "length over which the pressure of the leaky tube relaxes",
            ),
        ]
    sites = (ids.INFLOW, *(ids.zone(k) for k in range(1, n + 1)))
    for site, x, q_li in zip(sites, face_positions(n), flows, strict=True):
        out += _face_objects(site, x, tube.flow(x), q_li, flat)
    out += _mean_flow(n, flows)
    return out


def _face_objects(
    site: str, x: float, q_sin: float, q_li: float, flat: bool
) -> list[Parameter | AssignmentRule]:
    """The flow of one sinusoid and of the liver at the face of a site.

    The face of a site is the one the chain carries its flow across: the
    periportal face of the chain for `ids.INFLOW` (`x = 0`), the pericentral
    face of the zone for a zone.
    """
    if flat:
        formula = (
            f"({ids.PRESSURE_PORTAL} - {ids.PRESSURE_CENTRAL}) / "
            f"({ids.RESISTANCE_LUMEN} * {ids.LENGTH})"
        )
    else:
        big = f"({ids.LENGTH} / {ids.DECAY_LENGTH})"
        s = f"({x!r} dimensionless * {ids.LENGTH} / {ids.DECAY_LENGTH})"
        formula = (
            f"(({ids.PRESSURE_PORTAL} - {ids.PRESSURE_ONCOTIC}) * cosh({big} - {s}) - "
            f"({ids.PRESSURE_CENTRAL} - {ids.PRESSURE_ONCOTIC}) * cosh({s})) / "
            f"({ids.RESISTANCE_LUMEN} * {ids.DECAY_LENGTH} * sinh({big}))"
        )
    sin_id, li_id = ids.sinusoid_flow(site), ids.zone_flow(None, site)
    where = (
        f"at the periportal face of the chain (position {x:.3g})"
        if site == ids.INFLOW
        else f"at the pericentral face of zone {site[1:]} (position {x:.3g})"
    )
    return [
        _ruled(sin_id, q_sin, M3_PER_S, f"flow of one sinusoid {where}"),
        assignment_rule(sin_id, formula, M3_PER_S, f"flow of one sinusoid {where}"),
        _ruled(li_id, q_li, U.l_per_min, f"blood flow through the liver {where}"),
        assignment_rule(
            li_id,
            f"{sin_id} * {ids.N_SINUSOIDS} * 1000 l_per_m3 * 60 s_per_min",
            U.l_per_min,
            f"blood flow through the liver {where}",
        ),
    ]


def _mean_flow(n: int, flows: Sequence[float]) -> list[Parameter | AssignmentRule]:
    """The flow of the liver, the mean over the zones of the flows of its faces.

    The mean of the two faces of a zone, averaged over the zones, which is the
    trapezoidal rule over the faces (`mean_flow`): the ends count half,
    because each belongs to one zone only.
    """
    ends = (
        f"0.5 dimensionless * ({ids.zone_flow(None, ids.INFLOW)} + "
        f"{ids.zone_flow(None, ids.zone(n))})"
    )
    inner = [ids.zone_flow(None, ids.zone(k)) for k in range(1, n)]
    terms = " + ".join([ends, *inner])
    name = "blood flow through the liver (mean over the zones)"
    return [
        _ruled(ids.FLOW, mean_flow(flows), U.l_per_min, name),
        assignment_rule(ids.FLOW, f"({terms}) / {n} dimensionless", U.l_per_min, name),
    ]
