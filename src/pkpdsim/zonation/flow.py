"""The flow through the sinusoids: the data of the two ways to drive it.

This module holds the settings only. The mathematics of the flow along a
leaky tube and the SBML objects it emits are added with the flow of the
engine.
"""

from dataclasses import dataclass

from pkpdsim.pbpk.values import Q

__all__ = ["ImposedFlow", "PressureFlow"]


@dataclass(frozen=True)
class ImposedFlow:
    """A volumetric flow imposed on the organ.

    Attributes:
        flow: the flow through the whole liver, in a volume per time
            (l/min or ml/min); `None` takes the reference flow of the
            species.
    """

    flow: Q | None = None


@dataclass(frozen=True)
class PressureFlow:
    """A flow driven by the pressure drop over the sinusoid.

    Attributes:
        inlet_pressure: the pressure at the portal inlet, in kPa.
        outlet_pressure: the pressure at the central outlet, in kPa.
        conductance: the hydraulic conductance of the fenestrated wall, in
            ml/min/kPa per g of liver.
    """

    inlet_pressure: Q
    outlet_pressure: Q
    conductance: Q
