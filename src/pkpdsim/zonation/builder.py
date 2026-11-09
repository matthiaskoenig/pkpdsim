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

Every substance has a species `<s>_<layer>_<zone>` (mM) in every layer
compartment. The zones are a chain from periportal to pericentral:

- convection `FLOW_<s>_<layer>_z<k>_<layer>_z<k+1>` (`Q * c` of the upstream
  zone) carries the plasma and the erythrocytes (the vascular site of the
  well-stirred geometry) from one zone to the next; the plasma and the
  erythrocytes have the shares `1 - Hct` and `Hct` of the blood flow `Q_li`
  (`Q_plasma_li`, `Q_rbc_li`), so they have one velocity. The inflow
  `FLOW_<s>_<layer>_in_<layer>_z1` brings the concentration `<s>_<layer>_in` (a
  parameter, 0) into the first zone and the outflow
  `FLOW_<s>_<layer>_z<n>_<layer>_out` takes the last zone's concentration out.
- axial diffusion `D_<s>_<layer>_z<k>_<layer>_z<k+1>` (reversible, `D * A / dx *
  (c_k - c_k+1)`) in the plasma and the space of Disse of the sinusoid
  geometry. `A / dx` of zone `k` is `n**2 * V_k / L_sin**2` (the layer
  volume of the zone over the squared sinusoid length, the cross section
  `axial_area` of the whole liver over the zone length), so the conductance
  follows the liver mass as the volumes do (the two forms are equal because
  `V_k = A * L / n`, the layer volume of a zone is its cross section times
  its length). A zone of the well-stirred
  geometry is mixed and has no axial diffusion.

A flow of 0 and a diffusion of 0 emit no reaction, so the chain is closed or
has no diffusion by construction and no rate law is zero or divides by zero.
Convection is a first order upwind scheme: a chain of `n` zones disperses like
`n` tanks in series, a Peclet number of `2 n` (`effective_peclet`).
"""

import itertools
import warnings
from collections.abc import Sequence

from sbmlutils.factory import (
    Compartment,
    InitialAssignment,
    Model,
    Parameter,
    Reaction,
    Species,
)
from sbmlutils.metadata import BQB, SBO

from pkpdsim.database.tables import entry
from pkpdsim.pbpk.annotation import TAXA, initial_assignment
from pkpdsim.pbpk.placement import value_parameter
from pkpdsim.pbpk.substance import Substance
from pkpdsim.pbpk.units import MODEL_UNITS, U
from pkpdsim.pbpk.values import Value, number, unit_text
from pkpdsim.perfusion.ids import flow as flow_id
from pkpdsim.zonation import ids
from pkpdsim.zonation.flow import PressureFlow
from pkpdsim.zonation.geometry import (
    Layer,
    SinusoidGeometry,
    fractions,
    geometry,
    layer_volumes,
    reference_flow,
)
from pkpdsim.zonation.profiles import Profile
from pkpdsim.zonation.settings import Geometry, ZonationSettings
from pkpdsim.zonation.transport import diffusion_parameter

__all__ = ["ZonationBuilder", "effective_peclet"]

AXIAL_DIFFUSION = "axial_diffusion"
"""The process of the axial diffusion, which a profile scales by zone."""

_FLOW_TO_L_PER_MIN = {"l/min": 1.0, "ml/min": 1e-3}
_CARRIERS = (Layer.RBC, Layer.PLASMA)
_DIFFUSING = (Layer.PLASMA, Layer.DISSE)

_MASS_TO_KG = {"g": 1e-3, "kg": 1.0, "mg": 1e-6}
_VASCULAR = (Layer.RBC, Layer.PLASMA, Layer.ENDOTHELIUM, Layer.DISSE)
_P = SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER
_NAMES: dict[Layer | str, str] = {
    ids.VASCULAR: "vascular space",
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
            The processes: the axial diffusion (the zonal factor
            `p_axial_diffusion_<zone>` scales the conductance of an interface
            by the arithmetic mean of its two zones, a choice and not a result); the cell model and the barrier add
            theirs. This is the registry a later task extends with the
            processes it adds, and the validation of the profiles is only as
            strong as its contents.
        """
        return frozenset({AXIAL_DIFFUSION})

    def _geo(self) -> SinusoidGeometry:
        """The geometry: the given one or that of the species."""
        given = self.settings.geometry_constants
        return given if given is not None else geometry(self.settings.species)

    def _full_fractions(self) -> dict[Layer, float]:
        """The volume fractions of the layers of the sinusoid.

        `bile=None` drops the bile volume and does not renormalise the other
        fractions: the organ loses that fraction of its volume (0.4 % in the
        rat) instead of redistributing it. This is deliberate, and the same
        in both geometries.
        """
        fr = fractions(self.settings.species, self.settings.haematocrit)
        return {
            layer: f
            for layer, f in fr.items()
            if layer is not Layer.BILE or self.settings.bile is not None
        }

    def _fractions(self) -> dict[Layer | str, float]:
        """The volume fractions of the layers of a zone, by the geometry."""
        fr = self._full_fractions()
        if self.settings.geometry is Geometry.SINUSOID:
            return {**fr}
        out: dict[Layer | str, float] = {
            ids.VASCULAR: sum(fr.get(layer, 0.0) for layer in _VASCULAR)
        }
        for layer, f in fr.items():
            if layer not in _VASCULAR:
                out[layer] = f
        return out

    def _mass(self) -> Value:
        """The liver mass, in kg."""
        mass = self.settings.liver_mass
        if mass is not None:
            unit = getattr(mass.unit, "sid", mass.unit)
            value = number(mass)
            assert value is not None
            try:
                factor = _MASS_TO_KG[str(unit)]
            except KeyError:
                raise ValueError(
                    f"A liver mass needs a mass in {sorted(_MASS_TO_KG)}, got {unit}."
                ) from None
            return Value(value * factor, U.kg, mass.source, mass.notes)
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
            for layer in self._fractions():
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
            for layer in self._fractions():
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

    def _flow(self) -> float:
        """The blood flow through the whole liver (l/min).

        Raises:
            NotImplementedError: for a flow driven by a pressure, which the
                leaky tube of Task 9 solves.
            ValueError: for a flow which is no number or has no flow unit.
        """
        flow = self.settings.flow
        if isinstance(flow, PressureFlow):
            raise NotImplementedError(
                "A flow driven by the pressure needs the leaky tube of Task 9; "
                "only an imposed flow is built."
            )
        if flow.flow is None:
            return reference_flow(self.settings.species).value
        x = number(flow.flow)
        if x is None:
            raise ValueError(
                "The imposed flow needs a number, not a formula: a computed "
                "flow arrives with the pressure-driven flow of Task 9."
            )
        if not isinstance(flow.flow, Value):
            return x
        unit = unit_text(flow.flow.unit)
        try:
            return x * _FLOW_TO_L_PER_MIN[unit]
        except KeyError:
            raise ValueError(
                f"A flow needs a volume per time in {sorted(_FLOW_TO_L_PER_MIN)}, "
                f"got {unit}."
            ) from None

    def _haematocrit(self) -> float:
        """The haematocrit of the lumen, 0 without an erythrocyte layer."""
        fr = self._full_fractions()
        rbc = fr.get(Layer.RBC, 0.0)
        return rbc / (rbc + fr[Layer.PLASMA])

    def _velocity(self) -> float:
        """The velocity (m/s) of the blood in the sinusoid, 0 without a flow."""
        q = self._flow()
        if q == 0.0:
            return 0.0
        fr = self._full_fractions()
        geo = self._geo()
        lumen = fr[Layer.PLASMA] + fr.get(Layer.RBC, 0.0)
        litre = lumen * geo.f_par * self._mass().value / geo.density
        return geo.L_sin / (litre / q * 60.0)

    def _dispersion(self) -> float:
        """The numerical dispersion (m**2/s) of the upwind chain, `v * dx / 2`."""
        return self._velocity() * self._geo().L_sin / self.settings.n_zones / 2.0

    def _diffusion(self, substance: Substance) -> Parameter | None:
        """Get the axial diffusion of a substance, `None` for none.

        The aqueous diffusion of the substance; a subclass may override it
        (a tracer without diffusion).
        """
        return diffusion_parameter(substance)

    def _diffusions(self) -> dict[str, Parameter]:
        """The diffusion parameters of the substances which diffuse."""
        out: dict[str, Parameter] = {}
        for substance in self.substances:
            d = self._diffusion(substance)
            if d is not None:
                out[substance.sid] = d
        return out

    def _carriers(self) -> list[Layer | str]:
        """The layers which the blood flows through."""
        if self.settings.geometry is Geometry.WELL_STIRRED:
            return [ids.VASCULAR]
        return [layer for layer in _CARRIERS if layer in self._fractions()]

    def _carrier_flow(self, layer: Layer | str) -> str:
        """The id of the flow of a carrier layer."""
        return ids.FLOW if layer == ids.VASCULAR else ids.layer_flow(layer)

    def _species(self) -> list[Species]:
        """The species of every substance in every layer of every zone."""
        out: list[Species] = []
        for substance in self.substances:
            for k in range(1, self.settings.n_zones + 1):
                zone = ids.zone(k)
                for layer in self._fractions():
                    out.append(
                        Species(
                            ids.species(substance.sid, layer, zone),
                            ids.volume(layer, zone),
                            initialConcentration=0.0,
                            substanceUnit=U.mmole,
                            hasOnlySubstanceUnits=False,
                            name=f"{substance.name} ({_NAMES[layer]} of zone {k})",
                            sboTerm=SBO.SIMPLE_CHEMICAL,
                        )
                    )
        return out

    def _flow_parameters(self) -> tuple[list[Parameter], list[InitialAssignment]]:
        """The flow of the liver, its split between the carriers, the boundary."""
        q = self._flow()
        params = [
            value_parameter(
                ids.FLOW, q, U.l_per_min, "blood flow through the liver", _P
            )
        ]
        assignments: list[InitialAssignment] = []
        if self.settings.geometry is Geometry.SINUSOID:
            params.append(
                value_parameter(
                    ids.HAEMATOCRIT,
                    self._haematocrit(),
                    U.dimensionless,
                    "haematocrit of the sinusoidal lumen",
                    _P,
                )
            )
            share = {
                Layer.PLASMA: f"(1 - {ids.HAEMATOCRIT})",
                Layer.RBC: ids.HAEMATOCRIT,
            }
            for layer in self._carriers():
                pid = ids.layer_flow(layer)
                params.append(
                    Parameter(
                        pid,
                        q
                        * (
                            self._haematocrit()
                            if layer is Layer.RBC
                            else 1.0 - self._haematocrit()
                        ),
                        unit=U.l_per_min,
                        name=f"flow of the {_NAMES[layer]} through the liver",
                        constant=True,
                        sboTerm=_P,
                    )
                )
                assignments.append(
                    initial_assignment(
                        pid,
                        f"{ids.FLOW} * {share[Layer(layer)]}",
                        U.l_per_min,
                        f"flow of the {_NAMES[layer]} through the liver",
                    )
                )
        if q > 0.0:
            for substance in self.substances:
                for layer in self._carriers():
                    params.append(
                        Parameter(
                            ids.species(substance.sid, layer, ids.INFLOW),
                            0.0,
                            unit=U.mM,
                            name=f"{substance.name} in the {_NAMES[layer]} entering the liver",
                            constant=True,
                            sboTerm=_P,
                        )
                    )
        return params, assignments

    def _diffusion_parameters(self) -> list[Parameter]:
        """The sinusoid length, the diffusions and the profiles of the diffusion."""
        s = self.settings
        if s.geometry is Geometry.WELL_STIRRED or s.n_zones < 2:
            return []
        diffusions = self._diffusions()
        if not diffusions:
            return []
        params = [
            value_parameter(
                ids.LENGTH,
                Value(self._geo().L_sin, U.meter),
                U.meter,
                "length of the sinusoid",
                _P,
            ),
            *diffusions.values(),
        ]
        profile = s.profiles.get(AXIAL_DIFFUSION)
        if profile is not None:
            params += self._profile_parameters(AXIAL_DIFFUSION, profile)
        return params

    def _profile_parameters(self, process: str, profile: Profile) -> list[Parameter]:
        """The zonal factors `p_<process>_<zone>` of a profile."""
        values = profile.values(self.settings.n_zones)
        return [
            Parameter(
                ids.profile(process, ids.zone(k)),
                v,
                unit=U.dimensionless,
                name=f"zonal factor of the {process.replace('_', ' ')} in zone {k}",
                constant=True,
                sboTerm=_P,
            )
            for k, v in enumerate(values, start=1)
        ]

    def _flow_reaction(
        self, substance: Substance, layer: Layer | str, source: str, target: str
    ) -> Reaction:
        """The flow of a substance in a layer from one site to the next."""
        a = ids.species(substance.sid, layer, source)
        b = ids.species(substance.sid, layer, target)
        upstream = a
        edge = target if source == ids.INFLOW else source
        left = "" if source == ids.INFLOW else a
        right = "" if target == ids.OUTFLOW else b
        return Reaction(
            flow_id(substance.sid, ids.site(layer, source), ids.site(layer, target)),
            f"{left} -> {right}",
            formula=(f"{self._carrier_flow(layer)} * {upstream}", U.mmole_per_min),
            compartment=ids.volume(layer, edge),
            name=f"flow of {substance.name} in the {_NAMES[layer]} from {source} to {target}",
            sboTerm=SBO.TRANSPORT_REACTION,
        )

    def convection(self) -> list[Reaction]:
        """Get the convection along the chain, the inflow and the outflow.

        Returns:
            The reactions `FLOW_<s>_<layer>_<zone>_<layer>_<next>` of every
            substance and carrier layer (plasma and erythrocytes; the
            vascular site of the well-stirred geometry), the inflow
            `FLOW_<s>_<layer>_in_<layer>_z1` of the periportal boundary and
            the outflow `FLOW_<s>_<layer>_z<n>_<layer>_out` of the pericentral
            end; none for a flow of 0.

        Raises:
            NotImplementedError: for a flow driven by a pressure.
        """
        n = self.settings.n_zones
        if self._flow() == 0.0:
            return []
        out: list[Reaction] = []
        zones = [ids.zone(k) for k in range(1, n + 1)]
        chain = [ids.INFLOW, *zones, ids.OUTFLOW]
        for substance in self.substances:
            for layer in self._carriers():
                out += [
                    self._flow_reaction(substance, layer, a, b)
                    for a, b in itertools.pairwise(chain)
                ]
        return out

    def axial_diffusion(self) -> list[Reaction]:
        """Get the axial diffusion between consecutive zones.

        Returns:
            The reversible reactions `D_<s>_<layer>_<zone>_<layer>_<next>` of
            every diffusing substance in the plasma and the space of Disse of
            the sinusoid; none for a diffusion of 0, a single zone and the
            well-stirred geometry.
        """
        n = self.settings.n_zones
        if self.settings.geometry is Geometry.WELL_STIRRED or n < 2:
            return []
        has_profile = AXIAL_DIFFUSION in self.settings.profiles
        out: list[Reaction] = []
        for substance in self.substances:
            d = self._diffusion(substance)
            if d is None:
                continue
            for layer in _DIFFUSING:
                if layer not in self._fractions():
                    continue
                for k in range(1, n):
                    za, zb = ids.zone(k), ids.zone(k + 1)
                    a = ids.species(substance.sid, layer, za)
                    b = ids.species(substance.sid, layer, zb)
                    factor = (
                        f" * ({ids.profile(AXIAL_DIFFUSION, za)} + "
                        f"{ids.profile(AXIAL_DIFFUSION, zb)}) / 2 dimensionless"
                        if has_profile
                        else ""
                    )
                    out.append(
                        Reaction(
                            ids.axial(
                                substance.sid,
                                ids.site(layer, za),
                                ids.site(layer, zb),
                            ),
                            f"{a} <-> {b}",
                            formula=(
                                f"{d.sid} * {ids.volume(layer, za)} * {n**2} "
                                f"dimensionless / ({ids.LENGTH} * {ids.LENGTH}) * 60 "
                                f"s_per_min * ({a} - {b}){factor}",
                                U.mmole_per_min,
                            ),
                            compartment=ids.volume(layer, za),
                            name=(
                                f"axial diffusion of {substance.name} in the "
                                f"{_NAMES[layer]} between zones {k} and {k + 1}"
                            ),
                            sboTerm=SBO.TRANSPORT_REACTION,
                        )
                    )
        return out

    def _warn_dispersion(self) -> None:
        """Warn once when a diffusion lies below the numerical dispersion."""
        limit = self._dispersion()
        low = sorted(
            substance.sid
            for substance in self.substances
            if (d := self._diffusion(substance)) is not None
            and d.value is not None
            and 0.0 < float(d.value) < limit
        )
        if low:
            warnings.warn(
                f"The diffusion of {low} lies below the numerical dispersion "
                f"{limit:.3g} m**2/s (v * dx / 2) of {self.settings.n_zones} "
                "zones: the zone count and not the physics sets the gradient "
                "along the sinusoid.",
                UserWarning,
                stacklevel=3,
            )

    def build(self) -> Model:
        """Build the model of the parts that exist: the chain of layers and zones.

        Returns:
            The model.
        """
        s = self.settings
        flow_parameters, flow_assignments = self._flow_parameters()
        self._warn_dispersion()
        dispersion = self._dispersion()
        velocity = self._velocity()
        peclet = velocity * self._geo().L_sin / dispersion if dispersion > 0.0 else 0.0
        return Model(
            sid=f"zonated_liver_{s.species}",
            name=f"Zonated liver of the {s.species}",
            notes=(
                f"A liver of {s.n_zones} zone(s) from periportal to pericentral, "
                f"{s.geometry} geometry (pkpdsim.zonation). The upwind convection "
                f"disperses numerically with v * dx / 2 = {dispersion:.3g} m**2/s "
                f"(v = {velocity:.3g} m/s), a Peclet number of {peclet:.3g} "
                "without physical diffusion."
            ),
            annotations=[(BQB.HAS_TAXON, TAXA[s.species])],
            units=U,
            model_units=MODEL_UNITS,
            compartments=self.compartments(),
            species=self._species(),
            parameters=[
                *self._constants(),
                *flow_parameters,
                *self._diffusion_parameters(),
            ],
            assignments=[*self._assignments(), *flow_assignments],
            reactions=[*self.convection(), *self.axial_diffusion()],
        )


def effective_peclet(settings: ZonationSettings, d_eff: float) -> float:
    """Get the Peclet number a chain of zones realizes.

    The upwind convection of `n` zones adds the numerical dispersion
    `v * dx / 2` to the physical diffusion, `Pe = v * L / (d_eff + v * L / (2
    n))`, which is `2 n` when the physical diffusion is negligible.

    Args:
        settings: the settings (the zones, the flow, the species).
        d_eff: the effective axial diffusion (m**2/s).

    Returns:
        The Peclet number; 0 without a flow.

    Note:
        A rat liver at 10 zones with an aqueous diffusion coefficient is
        expected to trigger the warning of the build: `v * dx / 2` is of the
        order of 1e-9 m**2/s, the same as the diffusion. That is a property of
        the discretization and not a defect of the settings.
    """
    builder = ZonationBuilder(settings, (), None)
    v = builder._velocity()
    if v == 0.0:
        return 0.0
    length = builder._geo().L_sin
    return v * length / (d_eff + builder._dispersion())
