"""The settings of the zonated liver."""

from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import StrEnum

from pkpdsim.pbpk.values import Value
from pkpdsim.species import Species
from pkpdsim.zonation.flow import ImposedFlow, PressureFlow
from pkpdsim.zonation.geometry import SinusoidGeometry
from pkpdsim.zonation.profiles import Profile

__all__ = ["BileTransport", "Geometry", "ZonationSettings"]


class Geometry(StrEnum):
    """The geometry of a zone."""

    SINUSOID = "sinusoid"
    """Six layers from the lumen to the bile canaliculus, with transport between them."""

    WELL_STIRRED = "well_stirred"
    """A vascular, a cell and a bile compartment, the vascular one the lumen, endothelium and Disse."""


class BileTransport(StrEnum):
    """The transport along the bile canaliculus (not implemented before the bile of the engine)."""

    ADVECTIVE = "advective"
    """The bile flows towards the portal tract with the bile secretion."""

    DIFFUSIVE = "diffusive"
    """The bile mixes by diffusion along the canaliculus."""


@dataclass(frozen=True)
class ZonationSettings:
    """The settings of a zonated liver.

    Attributes:
        n_zones: the number of zones from periportal to pericentral, at least 1
            (one zone is the unzonated organ).
        geometry: the geometry of a zone.
        flow: how the flow through the sinusoids is given.
        bile: the transport in the bile canaliculus; `None` omits the bile
            compartment in both geometries.
        haematocrit: the haematocrit of the lumen; `None` takes the one of the
            species, 0 leaves out the erythrocyte layer.
        profiles: the zonal profile of a process by the name of the process.
        species: the species of the geometry and of the reference liver.
        liver_mass: the liver mass (kg, g or mg); `None` takes the reference
            liver of the species.
        geometry_constants: the geometry replacing the one the table gives
            for the species; `None` takes `geometry(species)`.
    """

    n_zones: int = 10
    geometry: Geometry = Geometry.SINUSOID
    flow: ImposedFlow | PressureFlow = field(default_factory=ImposedFlow)
    bile: BileTransport | None = BileTransport.ADVECTIVE
    haematocrit: float | None = None
    profiles: Mapping[str, Profile] = field(default_factory=dict)
    species: Species = Species.RAT
    liver_mass: Value | None = None
    geometry_constants: SinusoidGeometry | None = None

    def __post_init__(self) -> None:
        """Check the number of zones.

        Raises:
            ValueError: for less than one zone.
        """
        if self.n_zones < 1:
            raise ValueError(f"A liver has at least one zone, not {self.n_zones}.")
