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

- radial exchange (`radial_exchange`, sinusoid geometry only, the well-stirred
  vascular compartment contains the Disse space already): `FEN_<s>_<zone>`
  (plasma to Disse through the fenestrae), `END_<s>_<zone>_in` and `_out`
  (plasma to endothelial cell to Disse) and `PS_<s>_<zone>` (Disse to
  hepatocyte). The fenestral path is physics the engine knows, the paracellular
  permeability `P_fen_<s>` of an open pore (`fenestral_permeability`, from the
  aqueous diffusion, the radius of the molecule and the geometry) times the
  porosity of the zone `eps_fen_<zone>`, which rises from `porosity_pp` to
  `porosity_pc` of the table over the zone positions (a profile
  `porosity` multiplies it). The transcellular paths are model input of a
  transporter set and not of the molar mass: `P_endo_in_<s>`, `P_endo_out_<s>`
  and `P_hep_<s>` (m/s of the lateral area) are 0, the path is off, until a
  model sets them. The lateral areas `A_pd_<zone>` (lumen), `A_pe_<zone>` (abluminal face of the
  endothelial cell, at `r_sin + t_endo`) and `A_dh_<zone>` follow the
  mass: `n_sin * pi * 2 * r_sin * dx` is `2 V_lumen / r_sin` of the zone
  (`n_sin * pi * r_sin**2 * L_sin` is the volume of the whole lumen) and the
  Disse to hepatocyte interface is that area times `(r_sin + t_endo + t_disse)
  / r_sin`.

A flow of 0 and a diffusion of 0 emit no reaction, so the chain is closed or
has no diffusion by construction (a solute larger than a fenestra has no fenestral
reaction either) and no rate law is zero or divides by zero.
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
    bile_flow,
    fractions,
    geometry,
    layer_volumes,
    reference_flow,
    reference_liver_mass,
)
from pkpdsim.zonation.profiles import Profile, zone_positions
from pkpdsim.zonation.settings import BileTransport, Geometry, ZonationSettings
from pkpdsim.zonation.transport import (
    M_PER_S,
    diffusion_of,
    diffusion_parameter,
    fenestral_permeability,
    molecular_radius,
)

__all__ = ["ZonationBuilder", "effective_peclet"]

AXIAL_DIFFUSION = "axial_diffusion"
"""The process of the axial diffusion, which a profile scales by zone."""

POROSITY = "porosity"
"""The process of the fenestral porosity, which a profile scales by zone."""

ENDOTHELIAL_PERMEABILITY = "endothelial_permeability"
"""The process of the transcellular path of the endothelium."""

BASOLATERAL_PERMEABILITY = "basolateral_permeability"
"""The process of the exchange of the space of Disse with the hepatocyte."""

CANALICULAR_FLOW = "canalicular_flow"
"""The process of the bile flow through the canaliculus (advective bile)."""

CANALICULAR_DIFFUSION = "canalicular_diffusion"
"""The process of the diffusion along the canaliculus (diffusive bile)."""

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
        unused = sorted(set(settings.profiles) & self._unused_bile_processes())
        if unused:
            raise ValueError(
                f"The profiles {unused} scale a process the bile of the settings "
                f"({settings.bile}) does not have."
            )
        unknown = sorted(set(settings.profiles) - self.processes())
        if unknown:
            raise ValueError(
                f"The profiles {unknown} name no process of the model "
                f"(processes: {sorted(self.processes())})."
            )

    def _unused_bile_processes(self) -> frozenset[str]:
        """The processes of the bile the transport of the settings does not have."""
        match self.settings.bile:
            case BileTransport.ADVECTIVE:
                return frozenset({CANALICULAR_DIFFUSION})
            case BileTransport.DIFFUSIVE:
                return frozenset({CANALICULAR_FLOW})
            case _:
                return frozenset({CANALICULAR_FLOW, CANALICULAR_DIFFUSION})

    def processes(self) -> frozenset[str]:
        """Get the names of the processes a profile can scale.

        Returns:
            The processes: the axial diffusion (the zonal factor
            `p_axial_diffusion_<zone>` scales the conductance of an interface
            by the arithmetic mean of its two zones, a choice and not a result); the cell model and the barrier add
            theirs. The barrier adds the `porosity` (a factor on the gradient
            of the table, `p_porosity_<zone>`), the `endothelial_permeability`
            and the `basolateral_permeability`. The bile adds the
            `canalicular_flow` (the factor `p_canalicular_flow_<zone>` of the
            zone the bile leaves scales the flow of an interface) and the
            `canalicular_diffusion` (the mean of the two zones, as for the
            axial diffusion; the duct interface takes the factor of `z1`); a
            profile of the one the bile of the settings does not have is
            refused. This is the registry a later task extends with the
            processes it adds, and the validation of the profiles is only as
            strong as its contents.
        """
        return frozenset(
            {
                AXIAL_DIFFUSION,
                POROSITY,
                ENDOTHELIAL_PERMEABILITY,
                BASOLATERAL_PERMEABILITY,
                CANALICULAR_FLOW,
                CANALICULAR_DIFFUSION,
            }
        )

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
        """The blood flow through the whole liver (l/min), not negative.

        Raises:
            NotImplementedError: for a flow driven by a pressure, which the
                leaky tube of Task 9 solves.
            ValueError: for a flow which is no number, has no flow unit or is
                negative (it would reverse the chain).
        """
        q = self._imposed_flow()
        if q < 0.0:
            raise ValueError(
                f"The imposed flow {q} l/min is negative: it would run the "
                "blood from the pericentral to the periportal end, which is "
                "not a flow of the chain."
            )
        return q

    def _imposed_flow(self) -> float:
        """The imposed flow (l/min) as given, in l/min."""
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

    def velocity(self) -> float:
        """The velocity (m/s) of the blood in the sinusoid, 0 without a flow."""
        q = self._flow()
        if q == 0.0:
            return 0.0
        fr = self._full_fractions()
        geo = self._geo()
        lumen = fr[Layer.PLASMA] + fr.get(Layer.RBC, 0.0)
        litre = lumen * geo.f_par * self._mass().value / geo.density
        return geo.L_sin / (litre / q * 60.0)

    def dispersion(self) -> float:
        """The numerical dispersion (m**2/s) of the upwind chain, `v * dx / 2`."""
        return self.velocity() * self._geo().L_sin / self.settings.n_zones / 2.0

    def peclet(self, d_eff: float) -> float:
        """The Peclet number the chain realizes with a physical diffusion `d_eff`.

        `v * L / (d_eff + v * dx / 2)`: the one expression which the notes of
        the model, the warning and `effective_peclet` report; 0 without a flow.
        """
        v = self.velocity()
        if v == 0.0:
            return 0.0
        return v * self._geo().L_sin / (d_eff + self.dispersion())

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
        """The sinusoid length, the diffusions and the profiles of the diffusion.

        The axial diffusion of the blood and the diffusion along the
        canaliculus share the length and the diffusions of the substances.
        """
        s = self.settings
        axial = s.geometry is Geometry.SINUSOID and s.n_zones >= 2
        canalicular = s.bile is BileTransport.DIFFUSIVE
        if not (axial or canalicular):
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
        if axial and profile is not None:
            params += self._profile_parameters(AXIAL_DIFFUSION, profile)
        return params

    def _bile_parameters(self) -> list[Parameter]:
        """The bile flow, its reference and the profile of the canalicular transport."""
        s = self.settings
        if s.bile is None:
            return []
        if s.bile is BileTransport.DIFFUSIVE:
            profile = s.profiles.get(CANALICULAR_DIFFUSION)
            return (
                self._profile_parameters(CANALICULAR_DIFFUSION, profile)
                if profile is not None
                else []
            )
        flow = bile_flow(s.species)
        reference_mass = reference_liver_mass(s.species)
        params = [
            value_parameter(
                ids.BILE_FLOW_REFERENCE,
                flow,
                U.l_per_min,
                "bile flow of the reference liver",
                _P,
            ),
            value_parameter(
                ids.MASS_REFERENCE,
                reference_mass,
                U.kg,
                "mass of the reference liver of the bile flow",
                _P,
            ),
            Parameter(
                ids.BILE_FLOW,
                flow.value * self._mass().value / reference_mass.value,
                unit=U.l_per_min,
                name="bile flow through the liver",
                constant=True,
                sboTerm=_P,
            ),
        ]
        profile = s.profiles.get(CANALICULAR_FLOW)
        if profile is not None:
            params += self._profile_parameters(CANALICULAR_FLOW, profile)
        return params

    def _bile_assignments(self) -> list[InitialAssignment]:
        """The bile flow follows the liver mass: the reference flow per reference mass."""
        if self.settings.bile is not BileTransport.ADVECTIVE:
            return []
        return [
            initial_assignment(
                ids.BILE_FLOW,
                f"{ids.BILE_FLOW_REFERENCE} * {ids.MASS} / {ids.MASS_REFERENCE}",
                U.l_per_min,
                "bile flow through the liver",
            )
        ]

    def _duct_compartments(self) -> list[Compartment]:
        """The compartment of the duct, which holds amounts and has no volume of its own."""
        if self.settings.bile is None:
            return []
        return [
            Compartment(
                ids.volume(Layer.BILE, ids.DUCT),
                1.0,
                unit=U.litre,
                name="bile duct (a collector of amounts, no physical volume)",
                constant=True,
                sboTerm=SBO.PHYSICAL_COMPARTMENT,
            )
        ]

    def _duct_species(self) -> list[Species]:
        """The amount of every substance which has left the liver with the bile."""
        if self.settings.bile is None:
            return []
        return [
            Species(
                ids.species(substance.sid, Layer.BILE, ids.DUCT),
                ids.volume(Layer.BILE, ids.DUCT),
                initialAmount=0.0,
                substanceUnit=U.mmole,
                hasOnlySubstanceUnits=True,
                name=f"{substance.name} in the bile duct (amount)",
                sboTerm=SBO.SIMPLE_CHEMICAL,
            )
            for substance in self.substances
        ]

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

    def bile(self) -> list[Reaction]:
        """Get the transport along the canaliculus, from `z<n>` to `z1` and the duct.

        The bile runs the other way than the blood: pericentral to periportal,
        and what leaves the canaliculus of `z1` is collected in the duct
        (`<s>_bile_duct`, an amount), which gives nothing back.

        - `ADVECTIVE`: `FLOW_<s>_bile_z<k+1>_bile_z<k>` is `Q_bile * c` of the
          upstream (pericentral) zone, and `FLOW_<s>_bile_z1_bile_duct` takes
          the first zone out; a profile `canalicular_flow` scales every one
          by the factor of the zone it leaves.
        - `DIFFUSIVE`: `D_<s>_bile_z<k+1>_bile_z<k>` is reversible, `D * A / dx *
          (c_k+1 - c_k)`, with `A / dx = n**2 * V_k / L_sin**2` (the layer
          volume of a zone is its cross-section times its length, `V_k = A *
          L_sin / n`, as for the axial diffusion of the blood). The duct is
          an absorbing boundary (concentration 0) half a zone from the centre
          of `z1`, so `D_<s>_bile_z1_bile_duct` is `2 D A / dx * c_1`, without
          return. It needs no bile flow.

        Returns:
            The reactions of every substance; none for `bile=None`.
        """
        s = self.settings
        if s.bile is None:
            return []
        n = s.n_zones
        chain = [ids.zone(k) for k in range(n, 0, -1)] + [ids.DUCT]
        out: list[Reaction] = []
        for substance in self.substances:
            d = self._diffusion(substance)
            for source, target in itertools.pairwise(chain):
                if s.bile is BileTransport.ADVECTIVE:
                    out.append(self._bile_flow_reaction(substance, source, target))
                elif d is not None:
                    out.append(
                        self._bile_diffusion_reaction(substance, d, source, target)
                    )
        return out

    def _bile_flow_reaction(
        self, substance: Substance, source: str, target: str
    ) -> Reaction:
        """The advection of a substance in the canaliculus from one site to the next."""
        a = ids.species(substance.sid, Layer.BILE, source)
        b = ids.species(substance.sid, Layer.BILE, target)
        factor = (
            f" * {ids.profile(CANALICULAR_FLOW, source)}"
            if CANALICULAR_FLOW in self.settings.profiles
            else ""
        )
        return Reaction(
            flow_id(
                substance.sid,
                ids.site(Layer.BILE, source),
                ids.site(Layer.BILE, target),
            ),
            f"{a} -> {b}",
            formula=(f"{ids.BILE_FLOW} * {a}{factor}", U.mmole_per_min),
            compartment=ids.volume(Layer.BILE, source),
            name=(
                f"bile flow of {substance.name} in the canaliculus from {source} "
                f"to {target}"
            ),
            sboTerm=SBO.TRANSPORT_REACTION,
        )

    def _bile_diffusion_reaction(
        self, substance: Substance, d: Parameter, source: str, target: str
    ) -> Reaction:
        """The diffusion of a substance in the canaliculus between two sites."""
        n = self.settings.n_zones
        a = ids.species(substance.sid, Layer.BILE, source)
        b = ids.species(substance.sid, Layer.BILE, target)
        to_duct = target == ids.DUCT
        # the duct is an absorbing boundary half a zone away: twice the conductance
        gain = "2 dimensionless * " if to_duct else ""
        if CANALICULAR_DIFFUSION in self.settings.profiles:
            pa = ids.profile(CANALICULAR_DIFFUSION, source)
            pb = ids.profile(CANALICULAR_DIFFUSION, target)
            factor = f" * {pa}" if to_duct else f" * ({pa} + {pb}) / 2 dimensionless"
        else:
            factor = ""
        gradient = a if to_duct else f"({a} - {b})"
        return Reaction(
            ids.axial(
                substance.sid,
                ids.site(Layer.BILE, source),
                ids.site(Layer.BILE, target),
            ),
            f"{a} -> {b}" if to_duct else f"{a} <-> {b}",
            formula=(
                f"{gain}{d.sid} * {ids.volume(Layer.BILE, source)} * {n**2} "
                f"dimensionless / ({ids.LENGTH} * {ids.LENGTH}) * 60 s_per_min * "
                f"{gradient}{factor}",
                U.mmole_per_min,
            ),
            compartment=ids.volume(Layer.BILE, source),
            name=(
                f"diffusion of {substance.name} in the canaliculus between "
                f"{source} and {target}"
            ),
            sboTerm=SBO.TRANSPORT_REACTION,
        )

    def _radius(self, substance: Substance) -> Value:
        """Get the molecular radius of a substance (a subclass may override it)."""
        return molecular_radius(substance.mr)

    def _has_radial(self) -> bool:
        """Whether the zones have a barrier: a sinusoid with its radial layers."""
        fr = self._fractions()
        return self.settings.geometry is Geometry.SINUSOID and all(
            layer in fr
            for layer in (
                Layer.PLASMA,
                Layer.ENDOTHELIUM,
                Layer.DISSE,
                Layer.HEPATOCYTE,
            )
        )

    def _porosities(self) -> tuple[float, ...]:
        """The fenestral porosity of every zone: the gradient of the table."""
        geo = self._geo()
        return tuple(
            geo.porosity_pp + (geo.porosity_pc - geo.porosity_pp) * x
            for x in zone_positions(self.settings.n_zones)
        )

    def _lumen(self, zone: str) -> str:
        """The expression of the lumen volume (l) of a zone."""
        rbc = Layer.RBC in self._fractions()
        plasma = ids.volume(Layer.PLASMA, zone)
        return f"({plasma} + {ids.volume(Layer.RBC, zone)})" if rbc else plasma

    def _lumen_litre(self, k: int) -> float:
        """The lumen volume (l) of zone `k`."""
        n = self.settings.n_zones
        volumes = layer_volumes(self._mass(), n, self._geo(), self._full_fractions())
        return sum(
            volumes[(layer, k)].value
            for layer in (Layer.PLASMA, Layer.RBC)
            if (layer, k) in volumes
        )

    def _radial_parameters(self) -> list[Parameter]:
        """The geometry, the porosities, the lateral areas and the permeabilities."""
        if not self._has_radial():
            return []
        geo = self._geo()
        params = [
            value_parameter(
                ids.RADIUS,
                Value(geo.r_sin, U.meter),
                U.meter,
                "radius of the sinusoid",
                _P,
            ),
            value_parameter(
                ids.T_ENDO,
                Value(geo.t_endo, U.meter),
                U.meter,
                "thickness of the endothelium",
                _P,
            ),
            value_parameter(
                ids.T_DISSE,
                Value(geo.t_disse, U.meter),
                U.meter,
                "thickness of the space of Disse",
                _P,
            ),
        ]
        outer = (geo.r_sin + geo.t_endo + geo.t_disse) / geo.r_sin
        for k, eps in enumerate(self._porosities(), start=1):
            zone = ids.zone(k)
            area = 2.0 * self._lumen_litre(k) * 1e-3 / geo.r_sin
            params += [
                Parameter(
                    ids.porosity(zone),
                    eps,
                    unit=U.dimensionless,
                    name=f"fenestral porosity of the endothelium in zone {k}",
                    constant=True,
                    sboTerm=_P,
                ),
                Parameter(
                    ids.area("pd", zone),
                    area,
                    unit=U.m2,
                    name=f"lateral area between the plasma and the space of Disse in zone {k}",
                    constant=False,
                    sboTerm=_P,
                ),
                Parameter(
                    ids.area("pe", zone),
                    area * (geo.r_sin + geo.t_endo) / geo.r_sin,
                    unit=U.m2,
                    name=f"lateral area of the abluminal membrane of the endothelial cell in zone {k}",
                    constant=False,
                    sboTerm=_P,
                ),
                Parameter(
                    ids.area("dh", zone),
                    area * outer,
                    unit=U.m2,
                    name=f"lateral area between the space of Disse and the hepatocytes in zone {k}",
                    constant=False,
                    sboTerm=_P,
                ),
            ]
        for substance in self.substances:
            sid = substance.sid
            fen = fenestral_permeability(
                diffusion_of(substance), self._radius(substance), geo, 1.0
            )
            params.append(
                Parameter(
                    ids.permeability("fen", sid),
                    fen.value,
                    unit=M_PER_S,
                    name=(
                        f"permeability of {substance.name} through an open "
                        "fenestra (the porosity of the zone scales it)"
                    ),
                    constant=True,
                    sboTerm=_P,
                )
            )
            for path, text in (
                ("endo_in", "the luminal membrane of the endothelial cell"),
                ("endo_out", "the abluminal membrane of the endothelial cell"),
                ("hep", "the basolateral membrane of the hepatocyte"),
            ):
                params.append(
                    Parameter(
                        ids.permeability(path, sid),
                        0.0,
                        unit=M_PER_S,
                        name=f"permeability of {substance.name} through {text} (model input, 0: off)",
                        constant=True,
                        sboTerm=_P,
                    )
                )
        for process in (POROSITY, ENDOTHELIAL_PERMEABILITY, BASOLATERAL_PERMEABILITY):
            profile = self.settings.profiles.get(process)
            if profile is not None:
                params += self._profile_parameters(process, profile)
        return params

    def _radial_assignments(self) -> list[InitialAssignment]:
        """The lateral areas, which follow the volumes of the lumen."""
        if not self._has_radial():
            return []
        out: list[InitialAssignment] = []
        for k in range(1, self.settings.n_zones + 1):
            zone = ids.zone(k)
            area = f"2 dimensionless * {self._lumen(zone)} / {ids.RADIUS} / 1000 dimensionless"
            out.append(
                initial_assignment(
                    ids.area("pd", zone),
                    area,
                    U.m2,
                    f"lateral area between the plasma and the space of Disse in zone {k}",
                )
            )
            out.append(
                initial_assignment(
                    ids.area("pe", zone),
                    f"{ids.area('pd', zone)} * ({ids.RADIUS} + {ids.T_ENDO}) / "
                    f"{ids.RADIUS}",
                    U.m2,
                    f"lateral area of the abluminal membrane of the endothelial cell in zone {k}",
                )
            )
            out.append(
                initial_assignment(
                    ids.area("dh", zone),
                    f"{ids.area('pd', zone)} * ({ids.RADIUS} + {ids.T_ENDO} + "
                    f"{ids.T_DISSE}) / {ids.RADIUS}",
                    U.m2,
                    f"lateral area between the space of Disse and the hepatocytes in zone {k}",
                )
            )
        return out

    def _exchange(
        self,
        substance: Substance,
        kind: str,
        suffix: str,
        zone: str,
        layers: tuple[Layer, Layer],
        permeability: str,
        interface: str,
        text: str,
        factor: str = "",
    ) -> Reaction:
        """A reversible exchange of a substance between two layers of a zone."""
        a = ids.species(substance.sid, layers[0], zone)
        b = ids.species(substance.sid, layers[1], zone)
        porosity = f"{ids.porosity(zone)} * " if kind == "FEN" else ""
        return Reaction(
            ids.exchange(kind, substance.sid, zone, suffix),
            f"{a} <-> {b}",
            formula=(
                f"{permeability} * {porosity}{ids.area(interface, zone)} * 1000 "
                f"dimensionless * 60 s_per_min * ({a} - {b}){factor}",
                U.mmole_per_min,
            ),
            compartment=ids.volume(layers[0], zone),
            name=f"{text} of {substance.name} in zone {zone[1:]}",
            sboTerm=SBO.TRANSPORT_REACTION,
        )

    def radial_exchange(self) -> list[Reaction]:
        """Get the radial exchange between the layers of the zones.

        Returns:
            The reversible reactions per zone and substance: `FEN_<s>_<zone>`
            (plasma to the space of Disse through the fenestrae, the derived
            paracellular path `P_fen_<s> * eps_fen_<zone> * A_pd_<zone>`; none
            for a solute whose radius reaches the radius of a fenestra, which
            has no fenestral flux by construction), `END_<s>_<zone>_in`
            (plasma to endothelial cell) and `END_<s>_<zone>_out` (endothelial
            cell to Disse, over the abluminal area `A_pe_<zone>`), with `P_endo_in_<s>` and `P_endo_out_<s>`, and
            `PS_<s>_<zone>` (Disse to hepatocyte, `P_hep_<s>` over the area
            `A_dh_<zone>`). The transcellular permeabilities are model input
            and 0 (off) by default. None for the well-stirred geometry.

        Note:
            The fenestral path is a conductance. The transcellular path is two
            conductances in series through a compartment with volume, the
            endothelial cell, so the two paths add as parallel paths only in
            the quasi-steady state and not instantaneously: early after a
            change the endothelium is a sink (and drains the Disse space
            through its reversible abluminal membrane), so with a weak cell
            path both paths can deliver less than the fenestrae alone.

            A flow of 0 and a diffusion of 0 emit no reaction, but the
            transcellular reactions exist with a permeability of 0. That is no
            inconsistency: `Q` and the diffusion coefficient are structural
            settings fixed at the build, while `P_endo_in_<s>`,
            `P_endo_out_<s>` and `P_hep_<s>` are model input a caller sets on
            the built model, and a parameter which can change after the build
            needs its reaction present to act on. The fenestral reaction is
            structural (the radius of the solute is fixed at the build) and
            absent for a solute which cannot pass.
        """
        if not self._has_radial():
            return []
        profiles = self.settings.profiles
        geo = self._geo()
        out: list[Reaction] = []
        for substance in self.substances:
            sieve = fenestral_permeability(
                diffusion_of(substance), self._radius(substance), geo, 1.0
            ).value
            for k in range(1, self.settings.n_zones + 1):
                zone = ids.zone(k)

                def scaled(process: str, zone: str = zone) -> str:
                    return (
                        f" * {ids.profile(process, zone)}"
                        if process in profiles
                        else ""
                    )

                if sieve > 0.0:
                    out.append(
                        self._exchange(
                            substance,
                            "FEN",
                            "",
                            zone,
                            (Layer.PLASMA, Layer.DISSE),
                            ids.permeability("fen", substance.sid),
                            "pd",
                            "fenestral exchange",
                            scaled(POROSITY),
                        )
                    )
                endo = scaled(ENDOTHELIAL_PERMEABILITY)
                out.append(
                    self._exchange(
                        substance,
                        "END",
                        "in",
                        zone,
                        (Layer.PLASMA, Layer.ENDOTHELIUM),
                        ids.permeability("endo_in", substance.sid),
                        "pd",
                        "uptake into the endothelial cell",
                        endo,
                    )
                )
                out.append(
                    self._exchange(
                        substance,
                        "END",
                        "out",
                        zone,
                        (Layer.ENDOTHELIUM, Layer.DISSE),
                        ids.permeability("endo_out", substance.sid),
                        "pe",
                        "release of the endothelial cell into the space of Disse",
                        endo,
                    )
                )
                out.append(
                    self._exchange(
                        substance,
                        "PS",
                        "",
                        zone,
                        (Layer.DISSE, Layer.HEPATOCYTE),
                        ids.permeability("hep", substance.sid),
                        "dh",
                        "basolateral exchange",
                        scaled(BASOLATERAL_PERMEABILITY),
                    )
                )
        return out

    def _warn_dispersion(self) -> None:
        """Warn once when a diffusion lies below the numerical dispersion."""
        limit = self.dispersion()
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
        dispersion = self.dispersion()
        velocity = self.velocity()
        peclet = self.peclet(0.0)
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
            compartments=[*self.compartments(), *self._duct_compartments()],
            species=[*self._species(), *self._duct_species()],
            parameters=[
                *self._constants(),
                *flow_parameters,
                *self._diffusion_parameters(),
                *self._bile_parameters(),
                *self._radial_parameters(),
            ],
            assignments=[
                *self._assignments(),
                *flow_assignments,
                *self._radial_assignments(),
                *self._bile_assignments(),
            ],
            reactions=[
                *self.convection(),
                *self.axial_diffusion(),
                *self.radial_exchange(),
                *self.bile(),
            ],
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
    return ZonationBuilder(settings, (), None).peclet(d_eff)
