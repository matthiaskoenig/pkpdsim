"""The settings of the zonated liver."""

from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import StrEnum

from pkpdsim.sbml.units import U
from pkpdsim.sbml.values import Value, model_choice
from pkpdsim.species import Species
from pkpdsim.zonation.flow import ImposedFlow, PressureFlow
from pkpdsim.zonation.geometry import Layer, SinusoidGeometry, fractions
from pkpdsim.zonation.profiles import Profile

__all__ = ["BileTransport", "FlowBoundary", "Geometry", "ZonationSettings"]


class Geometry(StrEnum):
    """The geometry of a zone."""

    SINUSOID = "sinusoid"
    """Six layers from the lumen to the bile canaliculus, with transport between them."""

    WELL_STIRRED = "well_stirred"
    """A vascular, a cell and a bile compartment, the vascular one the lumen, endothelium and Disse."""


class BileTransport(StrEnum):
    """The transport along the bile canaliculus (not implemented before the bile of the engine)."""

    ADVECTIVE = "advective"
    """The bile flows towards the portal tract, collecting the secretion of the zones."""

    DIFFUSIVE = "diffusive"
    """The bile mixes by diffusion along the canaliculus."""


class FlowBoundary(StrEnum):
    """Who supplies the ends of the chain of zones."""

    OPEN = "open"
    """The engine: the concentration `<s>_<layer>_in` (a parameter, 0) flows into
    the first zone and the last zone flows out into nothing."""

    EXTERNAL = "external"
    """The caller (a perfusion circuit): the chain has the flows between the zones
    and neither an inflow nor an outflow, which the caller adds."""


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
        haematocrit: the haematocrit of the lumen entering the chain; `None`
            takes the one of the species, 0 leaves out the erythrocyte layer.
            The erythrocytes of a substance exchange with its plasma when the
            substance states how it distributes into them
            (`pkpdsim.pbpk.substance.Substance.erythrocytes`, read as
            `pkpdsim.zonation.transport.ErythrocyteDistribution`); a substance
            which states none keeps an erythrocyte species which exchanges
            with nothing and carries what a caller puts into it, the behaviour
            of every substance before the exchange existed.
        profiles: the zonal profile of a process by the name of the process.
        species: the species of the geometry and of the reference liver.
        liver_mass: the liver mass (kg, g or mg); `None` takes the reference
            liver of the species.
        reference_flow: the blood flow of the whole reference organ (l/min or
            ml/min), which an `ImposedFlow` without a value takes and against
            which a `PressureFlow` without a `resistance_factor` is
            calibrated; `None` (the default) states none, which makes those
            two an error naming what to give. The engine reads no physiology
            of its own: the caller gives the number, and for a species of the
            PBPK library it is the total liver flow of the reference
            individual, `pkpdsim.pbpk.physiology.reference_flow(species)`,
            which returns that `Value`.
        geometry_constants: the geometry replacing the one the table gives
            for the species; `None` takes `geometry(species)`.
        fractions: the volume fractions of the layers of the organ, replacing
            the table of the species (a perfused organ states its own); the
            fraction of a layer is of the parenchyma of the geometry
            (`f_par` of `geometry_constants`). `None` takes the table. The
            well-stirred geometry takes the vascular fraction as the sum of
            the lumen layers (`Layer.PLASMA`) and the cells as
            `Layer.CELL`. A fraction is a `Value` with its provenance, which
            reaches the notes of the parameter the builder writes
            (`f_vas`, `f_cell`); a plain number is wrapped as the model
            choice it is (`layer_fractions`) and never looks sourced.
        boundary: who supplies the ends of the chain.
        cell_species: every solute has a species in the cell layer; `False`
            when the cell model holds the species of the cell (well-stirred
            geometry only).
    """

    n_zones: int = 10
    geometry: Geometry = Geometry.SINUSOID
    flow: ImposedFlow | PressureFlow = field(default_factory=ImposedFlow)
    bile: BileTransport | None = BileTransport.ADVECTIVE
    haematocrit: float | None = None
    profiles: Mapping[str, Profile] = field(default_factory=dict)
    species: Species = Species.RAT
    liver_mass: Value | None = None
    reference_flow: Value | None = None
    geometry_constants: SinusoidGeometry | None = None
    fractions: Mapping[Layer, Value | float] | None = None
    boundary: FlowBoundary = FlowBoundary.OPEN
    cell_species: bool = True

    def layer_fractions(self) -> dict[Layer, Value] | None:
        """Get the fractions the caller states, every one with its provenance.

        Returns:
            The fraction of every layer the caller states as a `Value`, a
            plain number as the model choice of the caller it is; `None`
            where the caller states none and the table of the species
            applies.
        """
        if self.fractions is None:
            return None
        return {
            layer: f
            if isinstance(f, Value)
            else model_choice(
                f,
                U.dimensionless,
                f"the volume fraction of the {layer} in the organ, supplied "
                "by the caller of the engine and no value of a table",
            )
            for layer, f in self.fractions.items()
        }

    def full_fractions(self) -> dict[Layer, Value]:
        """Get the volume fractions of the layers of the build, with their provenance.

        The values of the table (`geometry.fractions`), or the fractions the
        caller states (`fractions`, a perfused organ its own), a plain number
        of which is the model choice it is (`layer_fractions`).

        `bile=None` drops the bile volume and does not renormalise the other
        fractions: the organ loses that fraction of its volume (0.4 % in the
        rat) instead of redistributing it. This is deliberate, and the same in
        both geometries.

        These are the fractions of the whole build: the compartments, the
        sinusoid count of the flow chain and the calibration of a resistance
        factor read them here, so a caller who states its own fractions gets
        one geometry and not two.

        Returns:
            The fraction of every layer of the build.
        """
        given = self.layer_fractions()
        fr = given if given is not None else fractions(self.species, self.haematocrit)
        return {
            layer: f
            for layer, f in fr.items()
            if layer is not Layer.BILE or self.bile is not None
        }

    def __post_init__(self) -> None:
        """Check the number of zones and the haematocrit.

        Raises:
            ValueError: for less than one zone, or a haematocrit outside [0, 1).
        """
        if not self.cell_species and self.geometry is not Geometry.WELL_STIRRED:
            raise ValueError(
                "cell_species=False needs the well-stirred geometry: the radial "
                "exchange of the sinusoid reads the species of the cell layer."
            )
        if self.n_zones < 1:
            raise ValueError(f"A liver has at least one zone, not {self.n_zones}.")
        h = self.haematocrit
        if h is not None and not 0.0 <= h < 1.0:
            raise ValueError(
                f"The haematocrit {h} is outside [0, 1): a haematocrit of 1 leaves "
                "no plasma, and without plasma there is no plasma chain."
            )
