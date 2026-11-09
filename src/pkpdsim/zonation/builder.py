"""The build of a zonated liver: one SBML model of the settings.

`ZonationBuilder` turns `ZonationSettings`, the substances and a cell model
into a model. The model has the liver mass `M_li` (kg), the density `rho_li`,
the parenchymal fraction `f_par` and the volume fraction `f_<layer>` of every
layer as constants, and the layers of the zones as compartments
`V_<layer>_<zone>` (l) with the initial assignment `f_<layer> * f_par * M_li /
rho_li / n`, so a change of the mass before a simulation sets all volumes, as
the zones of `pkpdsim.perfusion` follow their organ. The sinusoid has a
compartment per layer of `pkpdsim.zonation.geometry.fractions` (the
erythrocytes only with a haematocrit above 0); the well-stirred geometry a
vascular compartment (lumen, endothelium, Disse, the sum of their volumes), a cell compartment
(hepatocytes) and the bile compartment per zone; `bile=None` omits the bile
in both geometries.

The build grows by parts, the transport between the layers, the flow and the
bile are added to `build` by the next steps of the engine.
"""

from collections.abc import Sequence

from sbmlutils.factory import Compartment, InitialAssignment, Model, Parameter
from sbmlutils.metadata import BQB, SBO

from pkpdsim.database.tables import entry
from pkpdsim.pbpk.annotation import TAXA, initial_assignment
from pkpdsim.pbpk.placement import value_parameter
from pkpdsim.pbpk.substance import Substance
from pkpdsim.pbpk.units import MODEL_UNITS, U
from pkpdsim.pbpk.values import Value, number
from pkpdsim.zonation import ids
from pkpdsim.zonation.geometry import (
    Layer,
    SinusoidGeometry,
    fractions,
    geometry,
    layer_volumes,
)
from pkpdsim.zonation.settings import Geometry, ZonationSettings

__all__ = ["ZonationBuilder"]

_MASS_TO_KG = {"g": 1e-3, "kg": 1.0, "mg": 1e-6}
_VASCULAR = (Layer.RBC, Layer.PLASMA, Layer.ENDOTHELIUM, Layer.DISSE)
_P = SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER
_VASCULAR_NAME = "vascular space"
_NAMES: dict[Layer | str, str] = {
    ids.VASCULAR: _VASCULAR_NAME,
    Layer.RBC: "erythrocytes",
    Layer.PLASMA: "plasma",
    Layer.ENDOTHELIUM: "endothelium",
    Layer.DISSE: "space of Disse",
    Layer.HEPATOCYTE: "hepatocytes",
    Layer.BILE: "bile canaliculus",
}


class ZonationBuilder:
    """Build the SBML model of a zonated liver.

    Attributes:
        settings: the settings.
        substances: the substances of the model.
        cell_model: the cell model of the hepatocytes, supplied by the port of
            the perfusion core (`None` until then; it has no interface yet).
    """

    def __init__(
        self,
        settings: ZonationSettings,
        substances: Sequence[Substance],
        cell_model: object | None,
    ) -> None:
        """Check the settings.

        Args:
            settings: the settings.
            substances: the substances.
            cell_model: the cell model of the hepatocytes or `None`.

        Raises:
            ValueError: for a profile of a process nothing uses.
        """
        self.settings = settings
        self.substances = tuple(substances)
        self.cell_model = cell_model
        unknown = sorted(set(settings.profiles) - self.processes())
        if unknown:
            raise ValueError(
                f"The profiles {unknown} name no process of the model "
                f"(processes: {sorted(self.processes())})."
            )

    def processes(self) -> frozenset[str]:
        """Get the names of the processes a profile can scale.

        Returns:
            The processes; none before the transport and the cell model add
            theirs. This is the registry a later task extends with the
            processes it adds, and the validation of the profiles is only as
            strong as its contents.
        """
        return frozenset()

    def _geo(self) -> SinusoidGeometry:
        """The geometry: the given one or that of the species."""
        given = self.settings.geometry_constants
        return given if given is not None else geometry(self.settings.species)

    def _full_fractions(self) -> dict[Layer, float]:
        """The volume fractions of the layers of the sinusoid, without the bile if omitted."""
        fr = fractions(self.settings.species, self.settings.haematocrit)
        return {
            layer: f
            for layer, f in fr.items()
            if layer is not Layer.BILE or self.settings.bile is not None
        }

    def _fractions(self) -> dict[Layer | str, float]:
        """The volume fractions of the layers of a zone, by the geometry."""
        fr = self._full_fractions()
        out: dict[Layer | str, float] = {}
        if self.settings.geometry is Geometry.SINUSOID:
            for layer, f in fr.items():
                out[layer] = f
            return out
        out[ids.VASCULAR] = sum(fr.get(layer, 0.0) for layer in _VASCULAR)
        for layer, f in fr.items():
            if layer not in _VASCULAR:
                out[layer] = f
        return out

    def _layers(self) -> list[Layer | str]:
        """The tokens of the layers of a zone, from the lumen outwards."""
        return list(self._fractions())

    def _mass(self) -> Value:
        """The liver mass, in kg."""
        mass = self.settings.liver_mass
        if mass is not None:
            unit = getattr(mass.unit, "sid", mass.unit)
            value = number(mass)
            assert value is not None
            return Value(value * _MASS_TO_KG[str(unit)], U.kg, mass.source, mass.notes)
        section = "male.organs.li"
        table = f"physiology_{self.settings.species}"
        e = entry(table, section, "mass")
        return Value(e.value, U.kg, e.source, e.notes)  # ty: ignore[invalid-argument-type]

    def compartments(self) -> list[Compartment]:
        """Get the layer compartments of all zones.

        Returns:
            A compartment `V_<layer>_<zone>` (l) per layer and zone, with the
            volume of the reference liver.
        """
        n = self.settings.n_zones
        volumes = layer_volumes(self._mass(), n, self._geo(), self._full_fractions())
        out: list[Compartment] = []
        for k in range(1, n + 1):
            zone = ids.zone(k)
            for layer in self._layers():
                label = _NAMES[layer]
                if self.settings.geometry is Geometry.WELL_STIRRED:
                    # one source of truth: the vascular volume is the sum of its layers
                    volume = (
                        sum(
                            volumes[(v, k)].value
                            for v in _VASCULAR
                            if (v, k) in volumes
                        )
                        if layer == ids.VASCULAR
                        else volumes[(Layer(layer), k)].value
                    )
                else:
                    volume = volumes[(Layer(layer), k)].value
                out.append(
                    Compartment(
                        ids.volume(layer, zone),
                        volume,
                        unit=U.litre,
                        name=f"{label} of zone {k} of {n}",
                        sboTerm=SBO.PHYSICAL_COMPARTMENT,
                    )
                )
        return out

    def _constants(self) -> list[Parameter]:
        """The liver constants the volumes are computed from."""
        geo = self._geo()
        params = [
            value_parameter(ids.MASS, self._mass(), U.kg, "mass of the liver", _P),
            value_parameter(
                ids.DENSITY, geo.density, U.kg_per_l, "density of the liver", _P
            ),
            value_parameter(
                ids.PARENCHYMA,
                geo.f_par,
                U.dimensionless,
                "parenchymal fraction of the liver",
                _P,
            ),
        ]
        for layer, f in self._fractions().items():
            params.append(
                value_parameter(
                    ids.fraction(layer),
                    f,
                    U.dimensionless,
                    f"volume fraction of the {_NAMES[layer]} in the liver",
                    _P,
                )
            )
        return params

    def _assignments(self) -> list[InitialAssignment]:
        """The initial assignments of the compartments, which follow the mass."""
        n = self.settings.n_zones
        out: list[InitialAssignment] = []
        for k in range(1, n + 1):
            zone = ids.zone(k)
            for layer in self._layers():
                out.append(
                    initial_assignment(
                        ids.volume(layer, zone),
                        f"{ids.fraction(layer)} * {ids.PARENCHYMA} * {ids.MASS} / "
                        f"{ids.DENSITY} / {n} dimensionless",
                        U.litre,
                        f"{_NAMES[layer]} of zone {k} of {n}",
                    )
                )
        return out

    def build(self) -> Model:
        """Build the model of the parts that exist: the layer compartments.

        Returns:
            The model.
        """
        s = self.settings
        return Model(
            sid=f"zonated_liver_{s.species}",
            name=f"Zonated liver of the {s.species}",
            notes=(
                f"A liver of {s.n_zones} zone(s) from periportal to pericentral, "
                f"{s.geometry} geometry (pkpdsim.zonation)."
            ),
            annotations=[(BQB.HAS_TAXON, TAXA[s.species])],
            units=U,
            model_units=MODEL_UNITS,
            compartments=self.compartments(),
            parameters=self._constants(),
            assignments=self._assignments(),
        )
