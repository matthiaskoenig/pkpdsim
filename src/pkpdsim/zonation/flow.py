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

**The flow depends on the position.** What leaks through the fenestrae leaves
the lumen, so a pressure-driven sinusoid has a position-dependent flow and the
convection of the chain no longer carries one `Q` through all zones: every
zone `k` has the flow `Q_sin_<zone>` (m**3/s of one sinusoid) at its own
position `zone_positions(n)[k - 1]`, the same position the zonal profiles
use, and `Q_li_<zone>` (l/min, the liver has `N_sin` sinusoids, which follows
the liver mass). The convection out of a zone, and the inflow into the first
zone, read the flow of that zone (`ids.zone_flow`). spt evaluated the flow at
`x = 0` only (its `FLOW: Q dependent on position`); this finishes it. The
chain does not move the leaking fluid out of the lumen: the convection carries
the solute with the local flow from zone to zone, and the leak of the solute
is the fenestral exchange of the radial barrier.

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
for the linear gradient of the table; the zonal `porosity` profile scales the
solute path of the barrier and not the hydraulic resistance.

**The oncotic pressure** `P0` (`p_onc`, `Poc - Pot` of spt) defaults to the mean
of the two pressures, as spt sets `P0 = 0.5 (Pa + Pb)`. The net filtration is
zero then: what leaks out of the first half of the tube returns in the second
half, and the flow is symmetric about the middle. A default which filtered
would change the mass balance of every model which switches to the pressure
flow, so a lower oncotic pressure is a choice of the model.

A pressure which gives a negative flow in a zone is refused: the chain does
not run backwards.
"""

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

from sbmlutils.factory import InitialAssignment, Parameter
from sbmlutils.factory.units import UnitDefinition
from sbmlutils.metadata import SBO

from pkpdsim.pbpk.annotation import initial_assignment
from pkpdsim.pbpk.placement import value_parameter
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Q, Value, number, unit_text
from pkpdsim.zonation import ids
from pkpdsim.zonation.geometry import (
    Layer,
    SinusoidGeometry,
    fractions,
    reference_flow,
    sinusoid_count,
)
from pkpdsim.zonation.profiles import zone_positions

if TYPE_CHECKING:
    from pkpdsim.zonation.settings import ZonationSettings

__all__ = [
    "ImposedFlow",
    "LeakyTube",
    "PressureFlow",
    "fenestra_counts",
    "flow_objects",
    "imposed_flow",
    "leaky_tube",
    "zone_flows",
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


@dataclass(frozen=True)
class ImposedFlow:
    """A volumetric flow imposed on the organ.

    Attributes:
        value: the flow through the whole liver, in a volume per time
            (l/min or ml/min); a plain number is l/min.
        per_mass: the flow per mass of liver (ml/min/g, l/min/kg or
            ml/min/kg), which the flow of the liver follows, as the machine
            perfusion states it. `value` and `per_mass` exclude each other;
            neither takes the reference flow of the species.
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
            viscosity` (10 in the spt model).
        p_onc: the oncotic pressure `P0 = Poc - Pot` the pressure relaxes to
            along the leaky tube; `None` (the default) is the mean of the two
            pressures, as spt sets it, which makes the net filtration zero.
    """

    p_pp: Q = _PORTAL_PRESSURE
    p_pc: Q = _CENTRAL_PRESSURE
    viscosity: Q = _PLASMA_VISCOSITY
    resistance_factor: float = 10.0
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


def _mass_kg(mass: Q) -> float:
    """The mass of the liver in kg."""
    return _converted(mass, {"g": 1e-3, "kg": 1.0, "mg": 1e-6}, "liver mass")


def imposed_flow(flow: ImposedFlow, species: str, mass_kg: float) -> float:
    """Get an imposed flow of the liver in l/min.

    Args:
        flow: the imposed flow.
        species: the species of the reference flow.
        mass_kg: the mass of the liver (kg).

    Returns:
        The flow (l/min), as given, not checked for its sign.

    Raises:
        ValueError: for a flow which is a formula, has no flow unit, or a flow
            per mass without a unit.
    """
    if flow.per_mass is not None:
        if not isinstance(flow.per_mass, Value):
            raise ValueError(
                "A flow per mass needs a unit (ml/min/g, l/min/kg or ml/min/kg)."
            )
        return mass_kg * _converted(
            flow.per_mass, _PER_MASS_TO_L_PER_MIN_KG, "flow per mass"
        )
    if flow.value is None:
        return reference_flow(species).value  # ty: ignore[invalid-argument-type]
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
    is the one number of the endothelium (see the module docstring).

    Args:
        geo: the geometry.
        n_zones: the number of zones.

    Returns:
        The counts from periportal to pericentral.
    """
    area = math.pi * geo.r_fen**2
    return tuple(
        (geo.porosity_pp + (geo.porosity_pc - geo.porosity_pp) * x) / area
        for x in zone_positions(n_zones)
    )


def _mean_porosity(geo: SinusoidGeometry) -> float:
    """The porosity of the leaky tube, the mean of the periportal and pericentral porosity."""
    return 0.5 * (geo.porosity_pp + geo.porosity_pc)


def leaky_tube(flow: PressureFlow, geo: SinusoidGeometry) -> LeakyTube:
    """Get the leaky tube of a pressure flow and a geometry.

    Args:
        flow: the pressure flow.
        geo: the geometry.

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
    nu = flow.resistance_factor * _converted(
        flow.viscosity, _VISCOSITY_TO_PA_S, "viscosity"
    )
    return LeakyTube(
        r_sin=geo.r_sin,
        length=geo.L_sin,
        t_endo=geo.t_endo,
        r_fen=geo.r_fen,
        n_fen=_mean_porosity(geo) / (math.pi * geo.r_fen**2),
        viscosity=nu,
        p_pp=p_pp,
        p_pc=p_pc,
        p_onc=p_onc,
    )


def zone_flows(
    settings: "ZonationSettings", geo: SinusoidGeometry, mass: Q
) -> tuple[float, ...]:
    """Get the flow of the liver at every zone of a pressure flow (l/min).

    Args:
        settings: the settings, with a `PressureFlow`.
        geo: the geometry.
        mass: the liver mass.

    Returns:
        The flow of the sinusoids of the liver at the position of each zone.

    Raises:
        ValueError: for a pressure drop which gives a negative flow.
    """
    flow = settings.flow
    if not isinstance(flow, PressureFlow):
        raise TypeError("Only a pressure flow has a flow per zone.")
    tube = leaky_tube(flow, geo)
    n_sin = sinusoid_count(mass, geo, fractions(settings.species, settings.haematocrit))
    out = tuple(
        tube.flow(x) * n_sin * 1e3 * 60.0 for x in zone_positions(settings.n_zones)
    )
    if min(out) < 0.0:
        raise ValueError(
            f"The pressures {tube.p_pp} Pa (periportal) and {tube.p_pc} Pa "
            f"(pericentral) give a negative flow ({min(out):.3g} l/min): the "
            "blood would run from the pericentral to the periportal end, which "
            "is not a flow of the chain."
        )
    return out


def _constant(pid: str, value: float, unit: UnitDefinition, name: str) -> Parameter:
    """A constant parameter of a SI number."""
    return value_parameter(pid, Value(value, unit), unit, name, _P)


def _derived(
    pid: str, value: float, unit: UnitDefinition | object, name: str
) -> Parameter:
    """A constant parameter whose value an initial assignment sets again."""
    return Parameter(pid, value, unit=unit, name=name, constant=True, sboTerm=_P)  # ty: ignore[invalid-argument-type]


def flow_objects(
    settings: "ZonationSettings", geo: SinusoidGeometry, mass: Q
) -> list[Parameter | InitialAssignment]:
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
        The parameters and initial assignments.

    Raises:
        ValueError: for a negative flow, or a flow which is a formula.
    """
    flow = settings.flow
    if isinstance(flow, ImposedFlow):
        q = imposed_flow(flow, settings.species, _mass_kg(mass))
        if q < 0.0:
            raise ValueError(
                f"The imposed flow {q} l/min is negative: it would run the "
                "blood from the pericentral to the periportal end, which is "
                "not a flow of the chain."
            )
        return [
            value_parameter(
                ids.FLOW, q, U.l_per_min, "blood flow through the liver", _P
            )
        ]
    return _pressure_objects(settings, flow, geo, mass)


def _pressure_objects(
    settings: "ZonationSettings",
    flow: PressureFlow,
    geo: SinusoidGeometry,
    mass: Q,
) -> list[Parameter | InitialAssignment]:
    """The leaky tube as parameters and initial assignments."""
    tube = leaky_tube(flow, geo)
    flows = zone_flows(settings, geo, mass)
    n = settings.n_zones
    fr = fractions(settings.species, settings.haematocrit)
    f_lumen = fr[Layer.PLASMA] + fr.get(Layer.RBC, 0.0)
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
            tube.viscosity / flow.resistance_factor,
            PA_S,
            "viscosity of the plasma",
        ),
        _constant(
            ids.VISCOSITY_FACTOR,
            flow.resistance_factor,
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
            _mean_porosity(geo),
            U.dimensionless,
            "fenestral porosity the leaky tube is evaluated at (mean of the sinusoid)",
        ),
        _derived(
            ids.FENESTRA_COUNT,
            tube.n_fen,
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
    out: list[Parameter | InitialAssignment] = [*params, *assignments]
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
    for k, (x, q_li) in enumerate(zip(zone_positions(n), flows, strict=True), start=1):
        zone = ids.zone(k)
        out += _zone_objects(zone, x, tube.flow(x), q_li, flat)
    out += _mean_flow(n, flows)
    return out


def _zone_objects(
    zone: str, x: float, q_sin: float, q_li: float, flat: bool
) -> list[Parameter | InitialAssignment]:
    """The flow of one sinusoid and of the liver at the position of a zone."""
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
    sin_id, li_id = ids.sinusoid_flow(zone), ids.zone_flow(None, zone)
    where = f"at the position {x:.3g} of zone {zone[1:]}"
    return [
        _derived(sin_id, q_sin, M3_PER_S, f"flow of one sinusoid {where}"),
        initial_assignment(sin_id, formula, M3_PER_S, f"flow of one sinusoid {where}"),
        _derived(li_id, q_li, U.l_per_min, f"blood flow through the liver {where}"),
        initial_assignment(
            li_id,
            f"{sin_id} * {ids.N_SINUSOIDS} * 1000 l_per_m3 * 60 s_per_min",
            U.l_per_min,
            f"blood flow through the liver {where}",
        ),
    ]


def _mean_flow(n: int, flows: Sequence[float]) -> list[Parameter | InitialAssignment]:
    """The flow of the liver, the mean of the flows of its zones."""
    terms = " + ".join(ids.zone_flow(None, ids.zone(k)) for k in range(1, n + 1))
    return [
        _derived(
            ids.FLOW,
            sum(flows) / n,
            U.l_per_min,
            "blood flow through the liver (mean of the zones)",
        ),
        initial_assignment(
            ids.FLOW,
            f"({terms}) / {n} dimensionless",
            U.l_per_min,
            "blood flow through the liver (mean of the zones)",
        ),
    ]
