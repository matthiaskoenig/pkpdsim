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
(the cells, the hepatocytes of the liver) and the bile compartment per zone; `bile=None` omits the bile
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
  model sets them. The molar mass gates the fenestral path alone, because it is
  the radius of the molecule that it gives: a solute without one (a perfusate
  quantity, oxygen) has the transcellular and the basolateral path and no
  fenestral one. The lateral areas `A_pd_<zone>` (lumen), `A_pe_<zone>` (abluminal face of the
  endothelial cell, at `r_sin + t_endo`) and `A_dh_<zone>` follow the
  mass: `n_sin * pi * 2 * r_sin * dx` is `2 V_lumen / r_sin` of the zone
  (`n_sin * pi * r_sin**2 * L_sin` is the volume of the whole lumen) and the
  Disse to hepatocyte interface is that area times `(r_sin + t_endo + t_disse)
  / r_sin`.

- the filtrate (`filtrate`, a pressure-driven sinusoid only): the plasma which
  leaves the lumen through the fenestrae where the flow of the leaky tube falls
  takes solute with it (`SD_<s>_<zone>`, the solvent drag, scaled by
  `1 - sigma_<s>` of the reflection coefficient of the fenestrae), flows along
  the space of Disse (`FLOW_<s>_disse_z<k>_disse_z<k+1>` of `Q_fil_<zone>`) and
  is reabsorbed where the lumen regains flow. The water is conserved term by
  term, so a freely sieved solute keeps the concentration of the inflow through
  the whole chain and a solute no fenestra lets through is enriched by
  `Q(z1)/Q(z<k>)`, which is the physics of the two. The transmural flux is the
  leak of the lumen, `Q_li_<z<k-1>> - Q_li_<z<k>>`: an erythrocyte does not
  pass a fenestra, so the whole leak of the tube is the plasma's. The
  erythrocytes therefore keep their flow, the share `Hct` of the blood which
  enters the chain, in every zone (`Q_rbc_li_<zone>`), and the plasma carries
  the rest of the local flow of the tube (`Q_plasma_li_<zone> = Q_li_<zone> -
  Hct * Q_li_z1`); the three add up to the inflow in every zone. The static
  haematocrit, the ratio of the two volumes, stays the one of the inflow while
  the flow haematocrit rises towards the pericentral end: the engine does not
  zone the volumes of the two conduits.

- the exchange of the plasma with the erythrocytes (`erythrocyte_exchange`,
  sinusoid geometry at a haematocrit above 0): `RBC_<s>_<zone>`,
  `kex_bc_<s> * V_rbc_<zone> * (c_plasma - c_rbc / Kbc_<s>)` of every
  substance which states its distribution into the erythrocytes, the rate law
  and the parameters of the PBPK library. A substance in equilibrium gets a
  rate which is fast against the transit of a zone (`EQUILIBRATION`) and the
  inflow of its erythrocytes follows that of its plasma; a kinetic one gets
  its own rate and keeps its own boundary. The chain then carries the blood
  concentration of the library, `Q * BP * c_plasma` with
  `BP = 1 - Hct + Hct*Kbc`, as the sum of its two conduits. A substance which
  states no distribution keeps an erythrocyte species which exchanges with
  nothing and carries what a caller puts into it.

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
    AssignmentRule,
    Compartment,
    InitialAssignment,
    Model,
    Parameter,
    Reaction,
    Species,
)
from sbmlutils.metadata import BQB, SBO

from pkpdsim.database.tables import entry
from pkpdsim.sbml.annotation import TAXA, initial_assignment
from pkpdsim.sbml.parameters import value_parameter
from pkpdsim.sbml.units import MODEL_UNITS, OrganUnits, U
from pkpdsim.sbml.values import Value, check_unit, number
from pkpdsim.zonation import ids
from pkpdsim.zonation.flow import (
    M3_PER_S,
    PA,
    PA_S,
    PA_S_PER_M2,
    PA_S_PER_M4,
    PER_M2,
    PressureFlow,
    assignment_rule,
    flow_objects,
    imposed_flow,
    zone_flows,
)
from pkpdsim.zonation.geometry import (
    Layer,
    SinusoidGeometry,
    bile_flow,
    fractions,
    geometry,
    layer_volumes,
    mass_kg,
    porosities,
    reference_liver_mass,
)
from pkpdsim.zonation.ids import flow as flow_id
from pkpdsim.zonation.profiles import Profile
from pkpdsim.zonation.settings import (
    BileTransport,
    FlowBoundary,
    Geometry,
    ZonationSettings,
)
from pkpdsim.zonation.transport import (
    L_PER_M3,
    M2_PER_S,
    M_PER_S,
    ErythrocyteDistribution,
    Solute,
    axial_velocity,
    blood_cell_coefficient,
    diffuses,
    diffusion_of,
    diffusion_parameter,
    dispersion,
    erythrocytes_of,
    fenestral_permeability,
    molecular_radius,
    peclet,
    radius_ratio,
    reflection,
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

EQUILIBRATION = 100.0
"""How much faster than the transit of a zone the erythrocytes of a substance
in equilibrium exchange with the plasma.

The mode `equilibrium` of an erythrocyte distribution states that the
erythrocytes follow the plasma instantaneously. The engine has an erythrocyte
species per zone, because the haematocrit splits the lumen into two conduits,
so it cannot collapse the two into one and realizes the equilibrium as an
exchange which is fast against the residence time of the blood in a zone: the
relaxation of the difference `c_plasma - c_rbc/Kbc` has the rate
`kex_bc * (V_rbc/V_plasma + 1/Kbc)`, and `kex_bc` is chosen so that this rate
is `EQUILIBRATION` times the reciprocal of the shortest residence time of a
zone. A factor of 100 leaves the erythrocytes of a zone within about 1 % of
their equilibrium; a model which wants the departure itself states the mode
`kinetic` and its rate."""

_CARRIERS = (Layer.RBC, Layer.PLASMA)
_DIFFUSING = (Layer.PLASMA, Layer.DISSE)

_VASCULAR = (Layer.RBC, Layer.PLASMA, Layer.ENDOTHELIUM, Layer.DISSE)


class ZonationUnits(OrganUnits):
    """The unit definitions of the model of a zonated organ.

    Those of a model of an organ (`pkpdsim.sbml.units.OrganUnits`: the mass
    of an organ, the temperature, the oxygen partial pressure and the other
    units which `U` lacks) and the units of the transport and the pressure
    flow of the engine.
    """

    l_per_m3 = L_PER_M3
    m2_per_s = M2_PER_S
    m_per_s = M_PER_S
    m3_per_s = M3_PER_S
    Pa = PA
    Pa_s = PA_S
    Pa_s_per_m4 = PA_S_PER_M4
    Pa_s_per_m2 = PA_S_PER_M2
    per_m2 = PER_M2


_P = SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER
_NAMES: dict[Layer | str, str] = {
    ids.VASCULAR: "vascular space",
    Layer.RBC: "erythrocytes",
    Layer.PLASMA: "plasma",
    Layer.ENDOTHELIUM: "endothelium",
    Layer.DISSE: "space of Disse",
    Layer.CELL: "cells",
    Layer.BILE: "bile canaliculus",
}


class ZonationBuilder:
    """Build the SBML model of a zonated liver.

    Attributes:
        settings: the settings.
        substances: the substances of the model.
        cell_model: the cell model of the cells, supplied by the port of
            the perfusion core (`None` until then; it has no interface yet).
    """

    def __init__(
        self,
        settings: ZonationSettings,
        substances: Sequence[Solute],
        cell_model: object | None,
        cell_processes: Sequence[str] = (),
    ) -> None:
        """Check the settings.

        Args:
            settings: the settings.
            substances: the substances.
            cell_model: the cell model of the cells or `None`.
            cell_processes: the processes of the cell model which a profile
                can scale; the zonal factor `p_<process>_<zone>` of a profile
                replaces the port `p_<process>` of the cell model of the zone
                (`pkpdsim.zonation.composition.compose`).

        Raises:
            TypeError: for a substance which is no solute.
        """
        for substance in substances:
            if not isinstance(substance, Solute):
                raise TypeError(
                    f"{substance!r} is no solute: it needs an id (sid), a name, "
                    "a molar mass (mr) and a diffusion (diffusion, both may be None)."
                )
        self.settings = settings
        self.substances = tuple(substances)
        self.cell_model = cell_model
        self.cell_processes = tuple(cell_processes)

    def _check_profiles(self) -> None:
        """Check the profiles of the settings against the processes of the model.

        The check belongs to `build()`, which is where the profiles are
        consumed, and not to the constructor: a builder which is made to read
        one number of the settings (`observers.observer_objects`,
        `effective_peclet`) does not know the processes of the cell model, and
        a profile of such a process is the feature the engine exists for.

        Raises:
            ValueError: for a profile of a process nothing uses.
        """
        settings = self.settings
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
            theirs (`cell_processes`). The barrier adds the `porosity` (a factor on the gradient
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
                *self.cell_processes,
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
        given = self.settings.fractions
        fr = (
            dict(given)
            if given is not None
            else fractions(self.settings.species, self.settings.haematocrit)
        )
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
        """The liver mass, in kg (`geometry.mass_kg` converts it).

        Raises:
            ValueError: for a mass whose unit is no mass.
        """
        mass = self.settings.liver_mass
        if mass is not None:
            return Value(mass_kg(mass), U.kg, mass.source, mass.notes)
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
                        metaId=ids.meta_id(ids.volume(layer, zone)),
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

        A pressure flow has a flow which depends on the position: this is
        the mean of the flows of the zones, the flow which sets the velocity.

        Raises:
            ValueError: for a flow which is no number, has no flow unit or is
                negative (it would reverse the chain).
        """
        flow = self.settings.flow
        if isinstance(flow, PressureFlow):
            flows = self._zone_flows()
            return sum(flows) / len(flows)
        q = imposed_flow(
            flow,
            self.settings.reference_flow,
            self._mass().value,
            reference_liver_mass(self.settings.species).value,
        )
        if q < 0.0:
            raise ValueError(
                f"The imposed flow {q} l/min is negative: it would run the "
                "blood from the pericentral to the periportal end, which is "
                "not a flow of the chain."
            )
        return q

    def _zone_flows(self) -> tuple[float, ...]:
        """The flow of the liver (l/min) at every zone of a pressure flow."""
        return zone_flows(self.settings, self._geo(), self._mass())

    def _haematocrit(self) -> float:
        """The haematocrit of the lumen entering the chain, 0 without an erythrocyte layer.

        The volume fractions of the two conduits give it. It is the
        haematocrit of the blood which enters the first zone: the plasma of a
        pressure-driven sinusoid filters out along the lobule and the
        erythrocytes do not, so the flow haematocrit rises towards the
        pericentral end (`_plasma_flows`) while the static one, the ratio of
        the two volumes, is this number in every zone.
        """
        fr = self._full_fractions()
        rbc = fr.get(Layer.RBC, 0.0)
        return rbc / (rbc + fr[Layer.PLASMA])

    def _blood_flows(self) -> tuple[float, ...]:
        """The blood flow of the lumen (l/min) at every zone.

        The flow of the position of the zone for a pressure flow, the one
        imposed flow in every zone otherwise.
        """
        if isinstance(self.settings.flow, PressureFlow):
            return self._zone_flows()
        return (self._flow(),) * self.settings.n_zones

    def _erythrocyte_flow(self) -> float:
        """The flow of the erythrocytes (l/min), the same in every zone.

        An erythrocyte does not pass a fenestra, so the erythrocyte conduit
        neither loses nor gains water along the sinusoid: its flow is the
        share `Hct` of the blood which enters the chain, in every zone. The
        plasma carries the whole leak of the leaky tube, which is what makes
        the chain conservative; the share `Hct` of the local flow, which the
        engine carried before, let the erythrocytes lose water to nowhere and
        enriched a solute in them by `Q(z1)/Q(z<k>)`.
        """
        return self._haematocrit() * self._blood_flows()[0]

    def _plasma_flows(self) -> tuple[float, ...]:
        """The flow of the plasma (l/min) at every zone.

        The lumen of a zone carries the flow of the leaky tube at its
        position, of which the erythrocytes carry a constant share, so the
        plasma carries the rest: `Q_li_<zone> - Hct * Q_li_z1`. The three
        add up to the inflow in every zone, with the filtrate of the space of
        Disse (`_filtrate_flows`), so no water is lost.

        Raises:
            ValueError: for a flow which leaves the plasma conduit nothing:
                the tube would filter more plasma than enters it.
        """
        rbc = self._erythrocyte_flow()
        flows = tuple(q - rbc for q in self._blood_flows())
        if min(flows) <= 0.0 < self._blood_flows()[0]:
            raise ValueError(
                f"The flow of the lumen falls to {min(self._blood_flows()):.3g} "
                f"l/min, at or below the flow of the erythrocytes {rbc:.3g} "
                f"l/min (a haematocrit of {self._haematocrit():g}): the "
                "sinusoid would filter more plasma than enters it. Lower the "
                "haematocrit or the pressure drop across the endothelium."
            )
        return flows

    def velocity(self) -> float:
        """The velocity (m/s) of the blood in the sinusoid, 0 without a flow."""
        # TODO(geometry): both numbers are outputs of this geometry paired with
        # the PBPK reference flow (default settings), none is a measurement:
        # 110 um/s for the rat (the lower end of the intravital literature,
        # roughly 100 to 500 um/s) and 47 um/s for the human (below the roughly
        # 100 to 300 um/s of the literature). The velocity comes from pairing
        # the reference flow with the table geometry, not from the calibration
        # of the pressure flow.
        # Candidates: the derived sinusoid count (lumen fraction over
        # pi r_sin**2 L_sin), the sinusoid radius, or the sinusoid length; the
        # human entries of the table are `fallback: rat`. Do not adjust the
        # geometry to hit a velocity; open question for the maintainer.
        q = self._flow()
        if q == 0.0:
            return 0.0
        fr = self._full_fractions()
        geo = self._geo()
        lumen = fr[Layer.PLASMA] + fr.get(Layer.RBC, 0.0)
        litre = lumen * geo.f_par * self._mass().value / geo.density
        return axial_velocity(geo.L_sin, litre, q)

    def dispersion(self) -> float:
        """The numerical dispersion (m**2/s) of the upwind chain, `v * dx / 2`."""
        return dispersion(self.velocity(), self._geo().L_sin, self.settings.n_zones)

    def peclet(self, d_eff: float) -> float:
        """The Peclet number the chain realizes with a physical diffusion `d_eff`.

        `v * L / (d_eff + v * dx / 2)` of `pkpdsim.zonation.transport.peclet`,
        the one expression which the notes of the model, the warning,
        `effective_peclet` and the zonated liver of the PBPK library report;
        0 without a flow.
        """
        return peclet(self.velocity(), self._geo().L_sin, self.settings.n_zones, d_eff)

    def _diffusion(self, substance: Solute) -> Parameter | None:
        """Get the axial diffusion of a substance, `None` for none.

        The aqueous diffusion of the substance, none for a solute with
        neither a diffusion nor a molar mass; a subclass may override it (a
        tracer without diffusion).
        """
        return diffusion_parameter(substance) if diffuses(substance) else None

    def _diffusions(self) -> dict[str, Parameter]:
        """The diffusion parameters of the substances which diffuse."""
        out: dict[str, Parameter] = {}
        for substance in self.substances:
            d = self._diffusion(substance)
            if d is not None:
                out[substance.sid] = d
        return out

    def layers(self) -> list[Layer | str]:
        """Get the layers of a zone the geometry produced, from the lumen to the bile.

        Returns:
            The layers of the sinusoid geometry without the erythrocytes at
            a haematocrit of 0 and without the bile for `bile=None`; the
            vascular site, the hepatocyte and the bile of the well-stirred
            geometry.
        """
        return list(self._fractions())

    def has_flow(self) -> bool:
        """Get whether the blood flows through the chain (a flow above 0)."""
        return self._flow() > 0.0

    def carriers(self) -> list[Layer | str]:
        """The layers which the blood flows through."""
        if self.settings.geometry is Geometry.WELL_STIRRED:
            return [ids.VASCULAR]
        return [layer for layer in _CARRIERS if layer in self._fractions()]

    def carrier_flow(self, layer: Layer | str, zone: str) -> str:
        """The id of the flow of a carrier layer out of a zone.

        An imposed flow is one number for all zones; a pressure flow has the
        flow of the position of the zone.

        A **well-stirred** zone of a pressure flow is the exception: its one
        vascular compartment holds the lumen and the space of Disse, so the
        plasma which filters out of the lumen stays in it and flows on with
        it. The axial water flux of the compartment is the sum of the two,
        which is the flow entering the chain, the same in every zone
        (`Q_li_z1`); the position-dependent `Q_li_<zone>` of the lumen alone
        would make the zones lose water to nowhere and enrich every interior
        zone by `Q(z1)/Q(zk)`. The sinusoid geometry has the two conduits
        apart and carries the filtrate explicitly (`filtrate`).
        """
        if isinstance(self.settings.flow, PressureFlow):
            if self.settings.geometry is Geometry.WELL_STIRRED:
                return ids.zone_flow(layer, ids.zone(1))
            return ids.zone_flow(layer, zone)
        return ids.FLOW if layer == ids.VASCULAR else ids.layer_flow(layer)

    def _species(self) -> list[Species]:
        """The species of every substance in every layer of every zone."""
        out: list[Species] = []
        for substance in self.substances:
            for k in range(1, self.settings.n_zones + 1):
                zone = ids.zone(k)
                for layer in self._fractions():
                    if layer == Layer.CELL and not self.settings.cell_species:
                        continue
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

    def _flow_parameters(
        self,
    ) -> tuple[list[Parameter], list[InitialAssignment], list[AssignmentRule]]:
        """The flow of the liver, its split between the carriers, the boundary.

        A pressure flow adds the leaky tube and the flow of every zone
        (`flow_objects`), and the carriers split the flow of each zone
        (`_carrier_flows`). The
        flow chain (the flow of a zone, of a carrier) is assignment rules, so a
        change of the flow of the liver reaches all of it at any time.
        """
        q = self._flow()
        objects = flow_objects(self.settings, self._geo(), self._mass())
        params = [o for o in objects if isinstance(o, Parameter)]
        assignments = [o for o in objects if isinstance(o, InitialAssignment)]
        rules = [o for o in objects if isinstance(o, AssignmentRule)]
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
            for layer in self.carriers():
                for zone, flow, formula in self._carrier_flows(Layer(layer)):
                    pid = (
                        ids.layer_flow(layer)
                        if zone is None
                        else ids.zone_flow(layer, zone)
                    )
                    params.append(
                        Parameter(
                            pid,
                            flow,
                            unit=U.l_per_min,
                            name=f"flow of the {_NAMES[layer]} through the liver",
                            constant=False,
                            sboTerm=_P,
                        )
                    )
                    rules.append(
                        assignment_rule(
                            pid,
                            formula,
                            U.l_per_min,
                            f"flow of the {_NAMES[layer]} through the liver",
                        )
                    )
        if q > 0.0 and self.settings.boundary is FlowBoundary.OPEN:
            for substance in self.substances:
                for layer in self.carriers():
                    ruled = self._ruled_inflow(substance, layer)
                    params.append(
                        Parameter(
                            ids.species(substance.sid, layer, ids.INFLOW),
                            0.0,
                            unit=U.mM,
                            name=f"{substance.name} in the {_NAMES[layer]} entering the liver",
                            constant=not ruled,
                            sboTerm=_P,
                        )
                    )
        return params, assignments, rules

    def _carrier_flows(self, layer: Layer) -> list[tuple[str | None, float, str]]:
        """The flow of a carrier layer: the zone, the value and the formula.

        One entry without a zone for an imposed flow, which is the same in
        every zone, and one per zone for a pressure flow. The erythrocytes
        carry the share `Hct` of the blood entering the chain, which no zone
        changes (`_erythrocyte_flow`), and the plasma the rest of the local
        flow of the tube (`_plasma_flows`); the two shares of one flow are
        the same thing where the flow is the same in every zone.

        Args:
            layer: the carrier, the plasma or the erythrocytes.

        Returns:
            The entries, each the zone (`None` for all of them), the value
            (l/min) and the formula of its assignment rule.
        """
        if not isinstance(self.settings.flow, PressureFlow):
            q = self._flow()
            h = self._haematocrit()
            if layer is Layer.RBC:
                return [(None, h * q, f"{ids.FLOW} * {ids.HAEMATOCRIT}")]
            return [
                (
                    None,
                    (1.0 - h) * q,
                    f"{ids.FLOW} * (1 dimensionless - {ids.HAEMATOCRIT})",
                )
            ]
        rbc = f"{ids.zone_flow(None, ids.zone(1))} * {ids.HAEMATOCRIT}"
        flows = self._plasma_flows() if layer is Layer.PLASMA else None
        out: list[tuple[str | None, float, str]] = []
        for k in range(1, self.settings.n_zones + 1):
            zone = ids.zone(k)
            if flows is None:
                out.append((zone, self._erythrocyte_flow(), rbc))
            else:
                out.append((zone, flows[k - 1], f"{ids.zone_flow(None, zone)} - {rbc}"))
        return out

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
                constant=False,
                sboTerm=_P,
            ),
        ]
        profile = s.profiles.get(CANALICULAR_FLOW)
        if profile is not None:
            params += self._profile_parameters(CANALICULAR_FLOW, profile)
        return params

    def _bile_rules(self) -> list[AssignmentRule]:
        """The bile flow follows its reference and the liver mass, at any time."""
        if self.settings.bile is not BileTransport.ADVECTIVE:
            return []
        return [
            assignment_rule(
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
                metaId=ids.meta_id(ids.volume(Layer.BILE, ids.DUCT)),
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

    def _cell_profile_parameters(self) -> list[Parameter]:
        """The zonal factors of the processes of the cell model which have a profile."""
        out: list[Parameter] = []
        for process in self.cell_processes:
            profile = self.settings.profiles.get(process)
            if profile is not None:
                out += self._profile_parameters(process, profile)
        return out

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
        self, substance: Solute, layer: Layer | str, source: str, target: str
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
            formula=(
                f"{self.carrier_flow(layer, edge)} * {upstream}",
                U.mmole_per_min,
            ),
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
        chain = (
            [ids.INFLOW, *zones, ids.OUTFLOW]
            if self.settings.boundary is FlowBoundary.OPEN
            else zones
        )
        for substance in self.substances:
            for layer in self.carriers():
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
          L_sin / n`, as for the axial diffusion of the blood; the volume of
          an interface is the harmonic mean of its two zones). The duct is
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
        self, substance: Solute, source: str, target: str
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
        self, substance: Solute, d: Parameter, source: str, target: str
    ) -> Reaction:
        """The diffusion of a substance in the canaliculus between two sites."""
        n = self.settings.n_zones
        a = ids.species(substance.sid, Layer.BILE, source)
        b = ids.species(substance.sid, Layer.BILE, target)
        to_duct = target == ids.DUCT
        # The duct is an absorbing boundary (Dirichlet, c = 0) at the face of
        # z1, so the distance from the centre of z1 is dx / 2 and the
        # conductance D * A / (dx / 2) is twice that of a zone-to-zone
        # interface by construction (the standard treatment of a boundary at a
        # cell face); do not "fix" it to match the interior interfaces.
        gain = "2 dimensionless * " if to_duct else ""
        if CANALICULAR_DIFFUSION in self.settings.profiles:
            pa = ids.profile(CANALICULAR_DIFFUSION, source)
            pb = ids.profile(CANALICULAR_DIFFUSION, target)
            factor = f" * {pa}" if to_duct else f" * ({pa} + {pb}) / 2 dimensionless"
        else:
            factor = ""
        gradient = a if to_duct else f"({a} - {b})"
        va = ids.volume(Layer.BILE, source)
        # The face area of an interior interface is the harmonic mean of the two
        # zone areas, `2 Va Vb / (Va + Vb)` in volumes (`A_k = V_k n / L_sin`):
        # the two half-cells `dx / 2` are conductances in series, `1 / G = dx / 2
        # / (D A_a) + dx / 2 / (D A_b)`, so the conductance is one number of the
        # face and the diffusion is reciprocal. With the uniform bile volume of
        # today it is `V`, identical to the volume of the source zone alone, so
        # this changes no present number. The duct has one zone only.
        volume = (
            va
            if to_duct
            else f"(2 dimensionless * {va} * {ids.volume(Layer.BILE, target)} / "
            f"({va} + {ids.volume(Layer.BILE, target)}))"
        )
        return Reaction(
            ids.axial(
                substance.sid,
                ids.site(Layer.BILE, source),
                ids.site(Layer.BILE, target),
            ),
            f"{a} -> {b}" if to_duct else f"{a} <-> {b}",
            formula=(
                f"{gain}{d.sid} * {volume} * {n**2} "
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

    def _radius(self, substance: Solute) -> Value:
        """Get the molecular radius of a substance (a subclass may override it)."""
        assert substance.mr is not None  # a solute without a molar mass has no radius
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
                Layer.CELL,
            )
        )

    def _porosities(self) -> tuple[float, ...]:
        """The fenestral porosity of every zone: the gradient of the table."""
        return porosities(self.settings.n_zones, self._geo())

    def _reflection(self, substance: Solute) -> float:
        """The reflection coefficient of the endothelium for a substance.

        `1 - renkin(r / r_fen)` of `transport.reflection`; 0 (the filtered
        water takes the solute with it, unhindered) for a solute without a
        molar mass, which has no radius the engine could derive a hindrance
        from, as 0 is the permeability of a path which is model input.
        """
        if substance.mr is None:
            return 0.0
        return reflection(radius_ratio(self._radius(substance), self._geo()))

    def _has_filtrate(self) -> bool:
        """Whether the chain carries a filtrate in the space of Disse.

        A pressure-driven sinusoid loses plasma through the fenestrae where its
        flow falls and regains it where the flow rises; an imposed flow is the
        same in every zone and filters nothing, and the well-stirred geometry
        holds the lumen and the space of Disse in one compartment, so the
        filtrate never leaves it (and the chain carries the whole axial water
        flux, `carrier_flow`).
        """
        return (
            isinstance(self.settings.flow, PressureFlow)
            and self._has_radial()
            and self.has_flow()
        )

    def _filtrate_flows(self) -> tuple[float, ...]:
        """The filtrate flow (l/min) the space of Disse carries out of every zone.

        The plasma which has left the lumen up to the pericentral face of the
        zone, `Q_plasma(z1) - Q_plasma(zk)`: the lumen of the first zone
        receives the inflow at its own flow and so neither filters nor
        reabsorbs, which makes the first entry 0, and every further zone
        filters the flow its lumen loses. What filters is plasma: an
        erythrocyte does not pass a fenestra, so the whole leak of the tube
        is the plasma conduit's (`_plasma_flows`) and the filtrate is
        `Q_li_z1 - Q_li_<zone>`, the leak of the lumen.
        """
        flows = self._plasma_flows()
        return tuple(flows[0] - q for q in flows)

    def _filtrate_objects(self) -> tuple[list[Parameter], list[AssignmentRule]]:
        """The filtrate flow of every zone and the reflection of every substance.

        `Q_fil_<zone>` is an assignment rule of the plasma flows, so it follows
        a change of the pressures as the flows do; the first zone has none (it
        is 0 by construction).
        """
        if not self._has_filtrate():
            return [], []
        n = self.settings.n_zones
        flows = self._filtrate_flows()
        first = ids.zone_flow(Layer.PLASMA, ids.zone(1))
        params: list[Parameter] = []
        rules: list[AssignmentRule] = []
        for k in range(2, n + 1):
            zone = ids.zone(k)
            pid = ids.filtrate_flow(zone)
            name = (
                f"filtrate the space of Disse carries out of zone {k} of {n} "
                "(the plasma the lumen has lost through the fenestrae)"
            )
            params.append(
                Parameter(
                    pid,
                    flows[k - 1],
                    unit=U.l_per_min,
                    name=name,
                    constant=False,
                    sboTerm=_P,
                )
            )
            rules.append(
                assignment_rule(
                    pid,
                    f"{first} - {ids.zone_flow(Layer.PLASMA, zone)}",
                    U.l_per_min,
                    name,
                )
            )
        for substance in self.substances:
            derived = substance.mr is not None
            params.append(
                Parameter(
                    ids.reflection(substance.sid),
                    self._reflection(substance),
                    unit=U.dimensionless,
                    name=(
                        f"reflection coefficient of the endothelium for "
                        f"{substance.name}"
                        + ("" if derived else " (no molar mass: freely sieved)")
                    ),
                    constant=True,
                    sboTerm=_P,
                )
            )
        return params, rules

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
            # the molar mass gives the radius of the molecule, which the
            # fenestral path needs and no other path does: a solute without one
            # (a perfusate quantity) keeps the transcellular and the
            # basolateral path, whose permeabilities are model input in m/s
            if substance.mr is not None:
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
            area = (
                f"2 dimensionless * {self._lumen(zone)} / {ids.RADIUS} / 1000 l_per_m3"
            )
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
        substance: Solute,
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
                f"l_per_m3 * 60 s_per_min * ({a} - {b}){factor}",
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
            has no fenestral flux by construction, and none for a solute
            without a molar mass, which has no radius to derive the path
            from), `END_<s>_<zone>_in`
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
            sieve = (
                fenestral_permeability(
                    diffusion_of(substance), self._radius(substance), geo, 1.0
                ).value
                if substance.mr is not None
                else 0.0
            )
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
                        (Layer.DISSE, Layer.CELL),
                        ids.permeability("hep", substance.sid),
                        "dh",
                        "basolateral exchange",
                        scaled(BASOLATERAL_PERMEABILITY),
                    )
                )
        return out

    def _has_erythrocytes(self) -> bool:
        """Whether the zones have an erythrocyte layer which exchanges with the plasma.

        The sinusoid geometry at a haematocrit above 0. The well-stirred
        geometry holds the erythrocytes and the plasma in one compartment, so
        its erythrocytes are in equilibrium with its plasma by construction
        and there is nothing to exchange.
        """
        return (
            self.settings.geometry is Geometry.SINUSOID
            and Layer.RBC in self._fractions()
        )

    def _erythrocytes(self, substance: Solute) -> ErythrocyteDistribution | None:
        """The erythrocyte distribution of a substance, `None` for none.

        `None` for a substance which states none and for every substance
        where the zones have no erythrocyte layer to exchange with.
        """
        if not self._has_erythrocytes():
            return None
        return erythrocytes_of(substance)

    def _kbc(self, ery: ErythrocyteDistribution) -> float:
        """The erythrocyte to plasma coefficient of a distribution at the haematocrit."""
        return blood_cell_coefficient(ery, self._haematocrit())

    def _transit(self) -> float:
        """The shortest residence time of the blood in a zone (min), 0 without a flow.

        The volume of the lumen of a zone over the largest flow any zone
        carries: the time in which the erythrocytes of the fastest zone have
        to reach their equilibrium.
        """
        q = max(self._blood_flows())
        if q <= 0.0:
            return 0.0
        return self._lumen_litre(1) / q

    def _exchange_rate(
        self, sid: str, ery: ErythrocyteDistribution, kbc: float
    ) -> float:
        """The exchange rate of the erythrocytes of a distribution (1/min).

        `k_exchange` of the kinetic mode; `EQUILIBRATION` relaxations per
        residence time of a zone for the equilibrium mode, which the engine
        realizes as a fast exchange because its erythrocytes are a species of
        their own.

        Args:
            sid: the id of the substance, which an error names.
            ery: the distribution.
            kbc: its coefficient, which the derived rate reads.

        Returns:
            The rate (1/min).

        Raises:
            ValueError: for a kinetic mode whose rate is no number or is given
                in another unit than 1/min, and for the equilibrium mode of a
                chain without a flow, which gives the engine no time scale to
                be fast against.
        """
        if ery.mode != "equilibrium":
            k = number(ery.k_exchange)
            if k is None or k <= 0.0:
                raise ValueError(
                    f"The kinetic erythrocytes of '{sid}' need a positive "
                    f"exchange rate k_exchange (1/min), got {ery.k_exchange!r}."
                )
            problems = check_unit(sid, "k_exchange", ery.k_exchange, U.per_min)
            if problems:
                raise ValueError(problems[0])
            return k
        transit = self._transit()
        if transit <= 0.0:
            raise ValueError(
                "Erythrocytes in equilibrium need a flow: the engine makes "
                "the exchange fast against the residence time of a zone, and "
                "a chain without a flow has none. Give the mode 'kinetic' "
                "with its rate instead."
            )
        h = self._haematocrit()
        return EQUILIBRATION / transit / (h / (1.0 - h) + 1.0 / kbc)

    def _erythrocyte_parameters(self) -> list[Parameter]:
        """The coefficient and the exchange rate of every substance with erythrocytes."""
        out: list[Parameter] = []
        for substance in self.substances:
            ery = self._erythrocytes(substance)
            if ery is None:
                continue
            kbc = self._kbc(ery)
            equilibrium = ery.mode == "equilibrium"
            if number(ery.kbc) is not None:
                source = " (kbc of the substance)."
            else:
                source = (
                    f", inverted from the blood to plasma ratio {number(ery.bp):g} "
                    f"at the haematocrit {self._haematocrit():g} of the lumen "
                    "(BP = 1 - Hct + Hct*Kbc)."
                )
            boundary = (
                " The blood entering a chain whose erythrocytes are in "
                "equilibrium carries them in equilibrium with its plasma, so "
                f"the inflow {ids.species(substance.sid, Layer.RBC, ids.INFLOW)} "
                "is this coefficient times the one of the plasma and no "
                "independent input."
                if equilibrium
                else ""
            )
            out.append(
                Parameter(
                    ids.kbc(substance.sid),
                    kbc,
                    unit=U.dimensionless,
                    name=f"erythrocyte to plasma coefficient of {substance.name}",
                    constant=True,
                    sboTerm=_P,
                    notes=(
                        "The concentration of the erythrocytes over the one of "
                        f"the plasma at equilibrium{source}{boundary}"
                    ),
                )
            )
            rate = self._exchange_rate(substance.sid, ery, kbc)
            if equilibrium:
                how = (
                    f"Derived: {EQUILIBRATION:g} relaxations of "
                    "c_plasma - c_rbc/Kbc per residence time "
                    f"{self._transit():.3g} min of the fastest zone, which the "
                    "mode 'equilibrium' asks for (the engine has an erythrocyte "
                    "species of its own and cannot collapse it into the plasma)."
                )
            else:
                how = "Given by the substance (k_exchange)."
            out.append(
                Parameter(
                    ids.kbc_exchange(substance.sid),
                    rate,
                    unit=U.per_min,
                    name=(
                        f"exchange of {substance.name} between the plasma and "
                        "the erythrocytes"
                    ),
                    constant=True,
                    sboTerm=_P,
                    notes=f"The rate of {ery.mode} erythrocytes. {how}",
                )
            )
        return out

    def _ruled_inflow(self, substance: Solute, layer: Layer | str) -> bool:
        """Whether an assignment rule sets what a carrier brings into the first zone.

        The erythrocytes of a substance in equilibrium, whose inflow follows
        that of its plasma (`_erythrocyte_rules`) and is therefore no constant
        of the model.
        """
        if layer is not Layer.RBC:
            return False
        ery = self._erythrocytes(substance)
        return ery is not None and ery.mode == "equilibrium"

    def _erythrocyte_rules(self) -> list[AssignmentRule]:
        """The inflow of the erythrocytes of a substance in equilibrium.

        `<s>_rbc_in = Kbc_<s> * <s>_plasma_in`: the mode `equilibrium` says
        that the erythrocytes hold no state of their own, so the blood which
        enters the chain carries them in equilibrium with its plasma and the
        boundary is no independent input. Kinetic erythrocytes keep the free
        parameter, whose value is the history of the blood before the liver.
        Nothing without the inflow parameters (an external boundary, no flow).
        """
        if self.settings.boundary is not FlowBoundary.OPEN or not self.has_flow():
            return []
        out: list[AssignmentRule] = []
        for substance in self.substances:
            if not self._ruled_inflow(substance, Layer.RBC):
                continue
            pid = ids.species(substance.sid, Layer.RBC, ids.INFLOW)
            out.append(
                assignment_rule(
                    pid,
                    f"{ids.kbc(substance.sid)} * "
                    f"{ids.species(substance.sid, Layer.PLASMA, ids.INFLOW)}",
                    U.mM,
                    f"{substance.name} in the erythrocytes entering the liver",
                )
            )
        return out

    def erythrocyte_exchange(self) -> list[Reaction]:
        """Get the exchange of the plasma with the erythrocytes of every zone.

        `RBC_<s>_<zone>`: `kex_bc_<s> * V_rbc_<zone> * (c_plasma - c_rbc /
        Kbc_<s>)`, the rate law of the kinetic blood cells of the PBPK library
        (`pkpdsim.pbpk.organ`), so the two libraries describe the
        erythrocytes with one vocabulary. The erythrocytes approach
        `Kbc_<s> * c_plasma` with the rate `kex_bc_<s> / Kbc_<s>`.

        A substance in equilibrium gets the same reaction with a rate which
        is fast against the transit of a zone (`EQUILIBRATION`), because the
        erythrocytes of the engine are a species of their own: the chain then
        carries the blood concentration `Q * BP * c_plasma` of the library,
        `BP = 1 - Hct + Hct*Kbc`, as the sum of its two conduits.

        Returns:
            The reactions, one per substance with an erythrocyte distribution
            and zone; none for a substance without one (its erythrocyte
            species exchanges with nothing and keeps what a caller puts into
            it), none at a haematocrit of 0 (there is no erythrocyte layer)
            and none for the well-stirred geometry (one compartment holds
            both).
        """
        out: list[Reaction] = []
        for substance in self.substances:
            ery = self._erythrocytes(substance)
            if ery is None:
                continue
            for k in range(1, self.settings.n_zones + 1):
                zone = ids.zone(k)
                plasma = ids.species(substance.sid, Layer.PLASMA, zone)
                rbc = ids.species(substance.sid, Layer.RBC, zone)
                out.append(
                    Reaction(
                        ids.exchange("RBC", substance.sid, zone),
                        f"{plasma} <-> {rbc}",
                        formula=(
                            f"{ids.kbc_exchange(substance.sid)} * "
                            f"{ids.volume(Layer.RBC, zone)} * "
                            f"({plasma} - {rbc} / {ids.kbc(substance.sid)})",
                            U.mmole_per_min,
                        ),
                        compartment=ids.volume(Layer.PLASMA, zone),
                        name=(
                            f"exchange of {substance.name} between the plasma "
                            f"and the erythrocytes in zone {k}"
                        ),
                        sboTerm=SBO.TRANSPORT_REACTION,
                    )
                )
        return out

    def _upwind(self, flux: str, source: str, target: str) -> str:
        """The concentration a flux carries: the source's where it is positive.

        The first order upwind scheme of the chain, written so that the
        direction follows the sign of the flux at any time: a flux which the
        pressures of the model reverse carries the concentration of the site
        it then comes from.
        """
        return f"piecewise({source}, {flux} > 0 l_per_min, {target})"

    def filtrate(self) -> list[Reaction]:
        """Get the transport of the filtrate: the solvent drag and its chain.

        The plasma which leaves the lumen of a zone through the fenestrae
        takes solute with it (the solvent drag), flows along the space of
        Disse and is reabsorbed where the lumen regains flow. Both are the
        water balance of the pressure-driven sinusoid, written for the solute:

        - `SD_<s>_<zone>` (plasma to the space of Disse of the zone where its
          lumen loses flow, the reverse where it gains): `(1 - sigma_<s>) *
          (Q_fil_<zone> - Q_fil_<previous>) * c`, with the concentration of
          the side the water comes from (upwind). `sigma_<s>` is the
          reflection coefficient of the fenestrae
          (`transport.reflection`, the Renkin ratio of the fenestral
          permeability): a solute which no fenestra lets through
          (`sigma = 1`) stays in the lumen and is enriched by the filtration,
          which is the physics; a freely sieved one (`sigma = 0`) is carried
          with the water and keeps the concentration of the inflow.
        - `FLOW_<s>_disse_z<k>_disse_z<k+1>`: `Q_fil_<z<k>> * c` of the space
          of Disse, upwind, the filtrate flowing along the barrier. Without it
          the space of Disse of a zone would be a dead end, the net transmural
          flux would be zero at steady state whatever the drag, and the lumen
          would keep the whole enrichment (the defect this fixes).
        - `FLOW_<s>_disse_z<n>_disse_out` at the pericentral end, the filtrate
          which the lumen has not taken back (the hepatic lymph): the net
          filtration of the whole sinusoid, `Q_fil_<z<n>>`, which is 0 at the
          default oncotic pressure (the mean of the two pressures, where what
          filters out of the first half returns in the second) and is then not
          emitted. The reaction is the zero-gradient boundary of the space of
          Disse: a net filtration leaves at the concentration of the last
          zone, and a net reabsorption draws water in at that same
          concentration, so no gradient is invented at a boundary the engine
          has no interstitium behind. The observers count it in `cum_out_<s>`
          (`observers.observer_objects`), which keeps the mass balance exact.

        Returns:
            The reactions; none without a filtrate (`_has_filtrate`: an
            imposed flow, the well-stirred geometry, no flow) and none of the
            drag for a solute which is fully reflected.
        """
        if not self._has_filtrate():
            return []
        n = self.settings.n_zones
        flows = self._filtrate_flows()
        out: list[Reaction] = []
        for substance in self.substances:
            sid = substance.sid
            if self._reflection(substance) < 1.0:
                for k in range(2, n + 1):
                    out.append(self._drag_reaction(substance, k))
            for k in range(2, n):
                out.append(self._filtrate_reaction(sid, ids.zone(k), ids.zone(k + 1)))
            if flows[-1] != 0.0:
                out.append(
                    self._filtrate_reaction(sid, ids.zone(n), ids.OUTFLOW, lymph=True)
                )
        return out

    def _drag_reaction(self, substance: Solute, k: int) -> Reaction:
        """The solvent drag of a substance across the endothelium of zone `k`."""
        zone = ids.zone(k)
        a = ids.species(substance.sid, Layer.PLASMA, zone)
        b = ids.species(substance.sid, Layer.DISSE, zone)
        # the filtrate of the first zone is 0 by construction, so the flux of
        # the second zone is its own filtrate and not a difference
        flux = (
            ids.filtrate_flow(zone)
            if k == 2
            else f"({ids.filtrate_flow(zone)} - {ids.filtrate_flow(ids.zone(k - 1))})"
        )
        return Reaction(
            ids.exchange("SD", substance.sid, zone),
            f"{a} <-> {b}",
            formula=(
                f"(1 dimensionless - {ids.reflection(substance.sid)}) * {flux} * "
                f"{self._upwind(flux, a, b)}",
                U.mmole_per_min,
            ),
            compartment=ids.volume(Layer.PLASMA, zone),
            name=(
                f"solvent drag of {substance.name} with the filtered plasma of zone {k}"
            ),
            sboTerm=SBO.TRANSPORT_REACTION,
        )

    def _filtrate_reaction(
        self, sid: str, source: str, target: str, lymph: bool = False
    ) -> Reaction:
        """The filtrate carrying a substance along the space of Disse."""
        a = ids.species(sid, Layer.DISSE, source)
        b = "" if lymph else ids.species(sid, Layer.DISSE, target)
        flux = ids.filtrate_flow(source)
        where = (
            "out of the liver (the lymph of the net filtration)"
            if lymph
            else f"from {source} to {target}"
        )
        return Reaction(
            flow_id(sid, ids.site(Layer.DISSE, source), ids.site(Layer.DISSE, target)),
            f"{a} -> " if lymph else f"{a} <-> {b}",
            formula=(
                f"{flux} * {a}" if lymph else f"{flux} * {self._upwind(flux, a, b)}",
                U.mmole_per_min,
            ),
            compartment=ids.volume(Layer.DISSE, source),
            name=f"filtrate flow of {sid} in the space of Disse {where}",
            sboTerm=SBO.TRANSPORT_REACTION,
        )

    def filtrate_efflux(self, substance: Solute) -> str | None:
        """Get the flux (mmole/min) the filtrate takes out of the organ, `None` for none.

        The lymph of the net filtration at the pericentral end
        (`filtrate`), negative for a net reabsorption, which draws water in
        there; `None` when the sinusoid filters and reabsorbs the same amount
        (the default oncotic pressure) or carries no filtrate at all. The
        observers add it to the outflow ledger.
        """
        if not self._has_filtrate() or self._filtrate_flows()[-1] == 0.0:
            return None
        zone = ids.zone(self.settings.n_zones)
        return (
            f"{ids.filtrate_flow(zone)} * "
            f"{ids.species(substance.sid, Layer.DISSE, zone)}"
        )

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

        Raises:
            ValueError: for a profile of a process nothing uses
                (`_check_profiles`).
        """
        self._check_profiles()
        s = self.settings
        flow_parameters, flow_assignments, flow_rules = self._flow_parameters()
        filtrate_parameters, filtrate_rules = self._filtrate_objects()
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
            units=ZonationUnits,
            model_units=MODEL_UNITS,
            compartments=[*self.compartments(), *self._duct_compartments()],
            species=[*self._species(), *self._duct_species()],
            parameters=_unique(
                [
                    *self._constants(),
                    *flow_parameters,
                    *self._diffusion_parameters(),
                    *self._bile_parameters(),
                    *self._radial_parameters(),
                    *self._erythrocyte_parameters(),
                    *filtrate_parameters,
                    *self._cell_profile_parameters(),
                ]
            ),
            assignments=[
                *self._assignments(),
                *flow_assignments,
                *self._radial_assignments(),
            ],
            rules=[
                *flow_rules,
                *self._erythrocyte_rules(),
                *filtrate_rules,
                *self._bile_rules(),
            ],
            reactions=[
                *self.convection(),
                *self.axial_diffusion(),
                *self.radial_exchange(),
                *self.erythrocyte_exchange(),
                *self.filtrate(),
                *self.bile(),
            ],
        )


def _unique(parameters: Sequence[Parameter]) -> list[Parameter]:
    """The parameters without a repeated id: the geometry the flow and the barrier share.

    The leaky tube and the barrier (and the axial diffusion) read the same
    radius, thickness and length; every part emits what it reads, and the first
    of an id is kept (they carry the same value, which the test of the flow
    checks).
    """
    seen: dict[str, Parameter] = {}
    for p in parameters:
        seen.setdefault(str(p.sid), p)
    return list(seen.values())


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
