"""The zonated liver: the sinusoid from periportal to pericentral in zones.

`ZonatedLiver` builds the liver as `n` zones in series along the sinusoid,
from the periportal zone `z1` to the pericentral zone `z<n>`, the pattern of
the sinusoidal unit of the zonation models of pkdb_models (`misc/zonation`):
a chain of the same compartments, the chemistry placed into the cells of
every zone. Its boundary is that of the standard liver (`StandardOrgan`):
the same inputs (volumes, flows, `HCT`, fractions, substance parameters,
pool and sink species) and the same outputs, the concentrations of the
plasma (and the kinetic blood cells) leaving the liver `Cout_<s>`, so the
hepatic and the arm vein observers are those of the standard liver.

Every zone has the parts of the standard liver which a substance needs
(`pkpdsim.pbpk.parts`), with the share `f_zone = 1/n` of their volumes and of
the passive exchanges of the liver and the coefficients of the liver (`Kp`,
`K_int_pls`, `K_cell_pls`, `K_or`), so the zones together have the capacity
and the exchange of the standard liver. The blood (plasma and kinetic blood
cells) flows from the hepatic artery and the portal vein into the first
zone, from each zone into the next (`J_sin_<zone>_<part>_<s>`) and from the
last one into the venous pool; the portosystemic shunt `f_shunt` bypasses
all zones and the functional fraction `f_functional` scales the tissue of
every zone and the capacity of every process, as in the standard liver. The
species, compartments, unbound concentrations (observers), exchanges and
processes of a zone carry the zone in their ids: `Cpl_z1_<s>`, `Vti_z1`,
`Cu_ti_z1_<s>`, `J_ti_z1_<s>`, `DRG2M_z1`.

Every metabolism and transport of the liver is placed into every zone
(`pkpdsim.pbpk.placement.place`), with the share `w_<zone>` of the capacity
of the liver: by volume `f_zone`, or by the profile of its protein
(`profiles`, the relative abundance `p_<key>_<zone>` per zone, e.g. a
pericentral CYP3A4), `w = p*V_zone/sum(p*V_zone) = p/sum(p)`, so that the
total capacity `Vmax*V` (`k` of a first order process) is that of the
standard liver. A metabolism with a `profile` takes the profile of that
name, any other process the profile of its protein, if there is one. The
inhibitors act inline in every zone with their unbound concentration in the
zone; the bile is the sink of all zones. A protein which a substance
inactivates or induces has one turnover state `E_<protein>` for the liver
(as the intestine), driven by the unbound concentration in the cells of the
zones weighted by the placement of the protein: exact for the same
concentration in every zone, the mean of the zones along a gradient. A
custom reaction is placed in every zone with its rate law over the ids of
the zone (`CUSTOM_z1`), so it reads the species of its zone; the constants
of its rate law and the species it brings itself belong to the liver and are
created once, the species with the zone in their ids.

The zones are well-stirred tanks in series (`Geometry.WELL_STIRRED`). The
sinusoid geometry (`Geometry.SINUSOID`) adds the two transport steps of the
sinusoid which the well-stirred zone leaves out, with the geometry of the
species of the zonation engine (`pkpdsim.zonation`): the diffusion along the
sinusoid between the zones (`D_z1_pl_<s>`, in the plasma of the lumen and in
the interstitium, the space of Disse) and the fenestrated endothelium
between them, whose conductance `PS_vas_<zone>_<s> =
P_fen_<s>*eps_fen_<zone>*A_pd_<zone>` replaces the generic endothelium
`PS_vas` of the standard organ (the fenestral porosity rises from the
periportal to the pericentral end, 0.06 to 0.08 of the table of the engine,
Braet and Wisse 2002; the lateral area follows the lumen of the zone). The
layers of the engine are the parts of the liver: the plasma of the lumen is
the plasma, the erythrocytes the blood cells, the space of Disse the
interstitium and the hepatocyte the intracellular space; a substance whose
form has no interstitium has no barrier either (the space of Disse is inside
its tissue), and the endothelial cell and the bile canaliculus are no part
of the PBPK liver.

Neither step is a small correction of the chain, and the axial diffusion is
the larger of the two. For a small molecule its conductance between two
zones, `D*V_<part>_<zone>*n^2/L_sin^2`, is of the order of the plasma flow
through the chain (the Peclet number of the notes is about 6 for a drug of
0.3 kDa in ten zones, against the `2n = 20` of a substance which does not
diffuse), so it back-mixes the zones and partly undoes the zonation: the
venous concentration of such a drug deviates from the well-stirred liver by
about 7e-03 with the diffusion on and by 2e-05 with it off. A modeller who
chooses a zone count, or who wants the gradient along the lobule and not its
mean, needs that number; `axial_diffusion=False` leaves the term out. The
fenestral path, in contrast, limits nothing for a small molecule: its
conductance is hundreds of times the plasma flow, and the basolateral `PS`
of the substance is the step which limits the transfer into the tissue.

With one zone the zonated liver is the standard liver. With `n` zones and
linear kinetics the availability of a flow-limited liver at steady state is
`(1 + x/n)^-n` with `x = fu*CLint/(Q*BP)` instead of the well-stirred
`1/(1 + x)`: with the same intrinsic clearance the zonated liver extracts
more, the difference of the extraction is about `x/2*(1 - 1/n)` of it (second
order in `x`), and for large `n` it approaches the parallel tube, `exp(-x)`.
"""

import copy
import itertools
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

from sbmlutils.factory import ModelDefinition
from sbmlutils.factory.units import Units

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.annotation import annotate
from pkpdsim.pbpk.chemistry import CustomReaction, Metabolism, Transport
from pkpdsim.pbpk.organ import Boundary, OrganContext, StandardBuilder
from pkpdsim.pbpk.organ_builder import stationary_ids
from pkpdsim.pbpk.organs import Organ, Part
from pkpdsim.pbpk.physiology import reference_fractions, reference_liver_flow
from pkpdsim.pbpk.placement import (
    built_parts,
    formula_ids,
    in_zone,
    inhibitor_unbound,
    place,
    zone_ids,
)
from pkpdsim.pbpk.substance import Substance
from pkpdsim.pbpk.units import U
from pkpdsim.sex import Sex
from pkpdsim.zonation import ids as zids
from pkpdsim.zonation.geometry import (
    SinusoidGeometry,
    geometry,
    porosities,
    reference_liver_mass,
)
from pkpdsim.zonation.profiles import Profile
from pkpdsim.zonation.settings import Geometry
from pkpdsim.zonation.transport import (
    L_PER_M3,
    M2_PER_S,
    M_PER_S,
    axial_velocity,
    diffusion_of,
    diffusion_parameter,
    dispersion,
    fenestral_permeability,
    molecular_radius,
    peclet,
)
from pkpdsim.zones import third_weights, weighted_mean

__all__ = ["ZonatedLiver"]

_DIFFUSING: Mapping[Part, str] = {
    Part.PLASMA: "plasma",
    Part.INTERSTITIAL: "interstitium",
}
"""The parts a substance diffuses in along the sinusoid, in names.

The plasma of the sinusoidal lumen and the space of Disse, the two layers the
engine diffuses in; the erythrocytes and the cells do not diffuse axially."""


class SinusoidUnits(U):
    """The units of the library and the three of the sinusoid geometry.

    The unit definitions are those of the transport of the engine
    (`pkpdsim.zonation.transport`), so that the model of the liver and the
    model of the engine write one definition of each; the units of the
    pressure flow of the engine are none of the liver, and the module of its
    unit class imports the PBPK library, which imports this variant.
    """

    l_per_m3 = L_PER_M3
    m_per_s = M_PER_S
    m2_per_s = M2_PER_S


_REFERENCE_SEX = Sex.MALE
"""The sex of the one individual the numbers of the notes are of.

The mass of the reference liver of the geometry table of the engine is that
of the reference male of `physiology_<species>.toml`
(`pkpdsim.zonation.geometry.reference_liver_mass`), so the vascular fraction
and the hepatic blood flow of the notes are the male ones too: the velocity,
the dispersion and the Peclet number describe one individual, not a mix of
a male liver and a female flow. Every species of the geometry table has a
reference male."""


_PART_NAMES: Mapping[Part, str] = {
    Part.PLASMA: "plasma",
    Part.BLOOD_CELLS: "blood cells",
    Part.INTERSTITIAL: "interstitium",
    Part.INTRACELLULAR: "cells",
    Part.TISSUE: "tissue",
}
"""The parts of a zone in words, for the names of the means."""


_TABLE = (
    "the geometry of the liver of the species in the table of the zonation "
    "engine (`pkpdsim.zonation.geometry`, "
    "`pkpdsim/database/tables/zonation.toml`)"
)
"""Where a number of the sinusoid geometry comes from, for the notes."""


@dataclass(frozen=True)
class ZonatedLiver:
    """The liver as zones in series along the sinusoid.

    Attributes:
        n_zones: the number of zones, periportal to pericentral.
        profiles: the relative abundance of a protein in each zone, by the id
            of the protein or the `profile` of a metabolism; the processes of
            a protein without a profile are placed by volume. A sequence
            gives the abundance of every zone from the periportal to the
            pericentral one (`n_zones` values, not all 0); a
            `pkpdsim.zonation.profiles.Profile` is a shape of the lobular
            position which the engine evaluates at the zone count
            (`Linear(1.0, 4.0)` is the four-fold pericentral gradient at any
            number of zones), normalized to the mean 1.
        geometry: the geometry of a zone
            (`pkpdsim.zonation.settings.Geometry`): `WELL_STIRRED` (the
            default) is one mixed zone of the parts of the standard liver,
            `SINUSOID` adds the axial diffusion along the sinusoid and the
            fenestrated barrier between the plasma and the space of Disse of
            a substance in the full form.
        axial_diffusion: whether the zones of the sinusoid geometry diffuse
            into each other along the sinusoid (the default). `False` is the
            limit of a chain which only convects: no reaction between two
            zones, no diffusion coefficient and no length of the sinusoid
            are emitted, as the engine emits none for a substance without a
            diffusion (`ZonationBuilder._diffusion` returning `None`). A
            coefficient of 0 is no way to say it: a `Substance` needs a
            positive diffusion. Without the sinusoid geometry there is no
            axial diffusion either way.
        means: whether every zonated part gets the periportal and the
            pericentral mean of its concentration, `C<part>_pp_<s>` and
            `C<part>_pc_<s>` (mM): the mean of the first and the last third
            of the sinusoid, every zone weighted by its overlap with the
            third, the observers of the engine
            (`pkpdsim.zonation.observers`, the weights of
            `pkpdsim.zones.third_weights`). They are read-only
            and on request, as the observers of the engine are added on
            request and the readouts of the intestine
            (`PBPKDefinition.readouts`) are.
        name: the name of the component, the prefix of the id of the model
            definition (`zonated_li`).
    """

    n_zones: int = 10
    profiles: Mapping[str, Profile | Sequence[float]] = field(default_factory=dict)
    geometry: Geometry = Geometry.WELL_STIRRED
    axial_diffusion: bool = True
    means: bool = False
    name: str = "zonated"

    def __post_init__(self) -> None:
        """Check the zones and the profiles.

        Raises:
            ValueError: for less than one zone, a profile key which is no
                id, two keys which differ in case only, or a profile which is
                a sequence without a value per zone, with a negative value or
                only zeros.
        """
        if isinstance(self.n_zones, bool) or self.n_zones < 1:
            raise ValueError(
                f"The zonated liver needs at least one zone, not {self.n_zones}."
            )
        lower: dict[str, str] = {}
        resolved: dict[str, Profile | Sequence[float]] = {}
        for key, values in self.profiles.items():
            if not key.isidentifier():
                raise ValueError(f"Profile '{key}': the key must be an id.")
            if key.lower() in lower:
                raise ValueError(
                    f"Profiles '{lower[key.lower()]}' and '{key}' have the same id."
                )
            lower[key.lower()] = key
            if isinstance(values, Profile):
                resolved[key] = values
                continue
            if len(values) != self.n_zones:
                raise ValueError(
                    f"Profile '{key}' has {len(values)} values, the liver has "
                    f"{self.n_zones} zones."
                )
            if any(not math.isfinite(v) or v < 0.0 for v in values) or not any(values):
                raise ValueError(
                    f"Profile '{key}' needs finite values >= 0, not all 0: "
                    f"{list(values)}."
                )
            resolved[key] = tuple(float(v) for v in values)
        object.__setattr__(self, "profiles", resolved)

    def build(self, ctx: OrganContext) -> tuple[ModelDefinition, Boundary]:
        """Build the model of the liver and its boundary.

        Args:
            ctx: the context of the liver.

        Returns:
            The model definition and the boundary of the standard liver.

        Raises:
            ValueError: if the organ is not the liver, a process cannot be
                placed (a metabolism with a profile the liver does not give),
                the processes of a protein have different profiles, a profile
                is used by no process, the standard liver refuses a
                substance, or the sinusoid geometry has no geometry of the
                species or no fenestra which lets a substance through.
        """
        return ZonatedBuilder(ctx, self).build()


class ZonatedBuilder(StandardBuilder):
    """The SBML objects of the zonated liver: the standard liver in zones and its chemistry."""

    def __init__(self, ctx: OrganContext, liver: ZonatedLiver) -> None:
        """Start the collection of the elements of the zonated liver.

        Args:
            ctx: the context of the liver.
            liver: the component.

        Raises:
            ValueError: for another organ than the liver, a process the
                zonated liver cannot place, or a species whose liver the
                geometry table of the engine does not have (the sinusoid
                geometry only).
        """
        if ctx.organ is not Organ.LIVER:
            raise ValueError(
                f"The zonated liver is a variant of the liver, not of organ "
                f"'{ctx.organ.value}'."
            )
        super().__init__(ctx, liver.name)
        # the chain is the engine's: its zone tokens, its count and its geometry
        self.zones = tuple(zids.zone(i) for i in range(1, liver.n_zones + 1))
        self.profiles = liver.profiles
        self.geometry = liver.geometry
        self.axial = liver.axial_diffusion
        self.means = liver.means
        self.geo = geometry(ctx.species) if self._is_sinusoid else None
        self.customs = [p for p in self._processes() if isinstance(p, CustomReaction)]
        self.processes, self.profile_of = self._placed()
        self.protein_profile = {
            p.protein: self.profile_of[p.sid]
            for p in self.processes
            if p.protein is not None
        }

    def _notes(self) -> str:
        """The notes of the model definition.

        The well-stirred zones keep the notes they had before the sinusoid
        geometry, which the stored references of `tests/pbpk` pin line by
        line; the sinusoid names its two transport steps.
        """
        n = len(self.zones)
        sinusoid = f" {self._sinusoid_notes()}" if self._is_sinusoid else ""
        return (
            f"The {self.name} liver: {n} zones in series along the sinusoid, from "
            f"the periportal zone {self.zones[0]} to the pericentral zone "
            f"{self.zones[-1]}, each with the parts which the substances need and "
            "its share of the volumes, the exchanges and the chemistry of the "
            "liver; the boundary of the standard liver." + sinusoid
        )

    @property
    def _is_sinusoid(self) -> bool:
        """Whether the zones have the geometry of the sinusoid."""
        return self.geometry is Geometry.SINUSOID

    def _sinusoid_notes(self) -> str:
        """The two transport steps of the sinusoid and the numbers of its chain.

        The dispersion of the upwind chain and its Peclet number at the
        reference male of the species: a reader who chooses a zone count
        needs the order of magnitude of both, and the flow and the volume of
        the simulated individual are parameters of the model and no numbers
        at build time (`_reference_velocity`).
        """
        n = len(self.zones)
        geo = self._geometry()
        v = self._reference_velocity()
        axial = (
            "the plasma and the interstitium of neighbouring zones diffuse into "
            "each other and "
            if self.axial
            else "the zones do not diffuse into each other (`axial_diffusion` is "
            "off) and "
        )
        return (
            "The zones are tanks in series along the sinusoid: "
            + axial
            + "the plasma exchanges with the interstitium through the fenestrae "
            "of the endothelium, with the geometry of the liver of the species. "
            "The chain is upwind, so the convection from zone to zone disperses "
            "numerically with `v*L_sin/(2n)` of the velocity v of the plasma: at "
            f"the reference male of the species (the liver of "
            f"{reference_liver_mass(self.ctx.species).value:.3g} kg and its blood "
            f"flow {self._reference_flow():.3g} l/min) the blood moves at "
            f"{v * 1e6:.3g} um/s, the chain of {n} zones disperses with "
            f"{dispersion(v, geo.L_sin, n):.3g} m**2/s and its Peclet number is "
            f"{peclet(v, geo.L_sin, n, 0.0):.3g} = 2n for a substance whose "
            "axial diffusion is negligible against the dispersion "
            "(`pkpdsim.zonation.transport.peclet`, the `effective_peclet` of the "
            "engine)"
            + (
                "; the notes of the diffusion coefficient of a substance give "
                "its own Peclet number."
                if self.axial
                else ", which is the Peclet number of this chain: it has no "
                "axial diffusion."
            )
        )

    def _reference_flow(self) -> float:
        """The hepatic blood flow of the reference male of the species (l/min)."""
        return reference_liver_flow(self.ctx.species, _REFERENCE_SEX)

    def _reference_velocity(self) -> float:
        """The velocity of the blood in the sinusoid of the reference individual (m/s).

        The length of the sinusoid over the residence time of its lumen, the
        velocity of the engine (`ZonationBuilder.velocity`,
        `pkpdsim.zonation.transport.axial_velocity`) with the lumen and the
        flow of one individual, the reference male of the species: the
        vascular fraction of the liver of `physiology_<species>.toml` times
        the mass of its reference liver over the density of the liver, and
        the hepatic blood flow of the same individual. The mass is that of
        the reference male (`reference_liver_mass` of the engine reads
        `male.organs.li`, the liver the geometry table belongs to), so the
        fraction and the flow are the male ones as well and not those of the
        sex of the build: a number of a mixed individual would belong to
        nobody. The volume and the flow of the simulated individual are
        parameters of the model, so the velocity of the notes is that of the
        reference male and not of the simulation.

        It is about half of the velocity of the engine (22 against 47 um/s
        for the human): the lumen of a zone of the PBPK liver is the whole
        vascular space of the liver (`f_vas`, the sinusoids and the larger
        vessels), the lumen of the engine only the sinusoidal one. The
        number describes the chain of this model, whose zones hold that
        volume.
        """
        geo = self._geometry()
        f_vas = reference_fractions(self.organ, self.ctx.species)[_REFERENCE_SEX][0]
        lumen = f_vas * reference_liver_mass(self.ctx.species).value / geo.density
        return axial_velocity(geo.L_sin, lumen, self._reference_flow())

    def _units(self) -> type[Units]:
        """The unit definitions of the model of the liver.

        The sinusoid needs the units of its transport (`SinusoidUnits`:
        m/s, m**2/s, l/m**3); the well-stirred liver has the units of the
        library, so that its model is that of the standard liver.
        """
        return SinusoidUnits if self._is_sinusoid else super()._units()

    def _geometry(self) -> SinusoidGeometry:
        """The geometry of the sinusoid of the species of the build."""
        assert self.geo is not None
        return self.geo

    def _placed(
        self,
    ) -> tuple[list[Metabolism | Transport], dict[str, str | None]]:
        """The metabolism and transport of the liver and their profiles (`None`: by volume).

        A custom reaction has no profile: it states its rate itself and is
        placed in every zone with the ids of the zone (`self.customs`).

        Raises:
            ValueError: for a metabolism with a profile the liver does not
                give, the processes of a protein with different profiles, or
                a profile which no process uses.
        """
        processes: list[Metabolism | Transport] = []
        profile_of: dict[str, str | None] = {}
        for p in self._processes():
            if isinstance(p, CustomReaction):
                continue
            processes.append(p)
            profile_of[p.sid] = self._profile(p)
        by_protein: dict[str, set[str | None]] = {}
        for p in processes:
            if p.protein is not None:
                by_protein.setdefault(p.protein, set()).add(profile_of[p.sid])
        for protein, keys in by_protein.items():
            if len(keys) > 1:
                sids = sorted(p.sid for p in processes if p.protein == protein)
                raise ValueError(
                    f"The processes {sids} of protein '{protein}' have different "
                    f"profiles {sorted(str(k) for k in keys)}: a protein has one "
                    "abundance profile in the zonated liver."
                )
        unused = sorted(set(self.profiles) - set(profile_of.values()))
        if unused:
            raise ValueError(
                f"Profiles {unused} of the zonated liver are used by no process "
                "of the liver (by its protein or the profile of a metabolism)."
            )
        return processes, profile_of

    def _profile(self, p: Metabolism | Transport) -> str | None:
        """The profile of a process: of its name (metabolism) or of its protein.

        Raises:
            ValueError: for a metabolism with a profile the liver does not give.
        """
        if isinstance(p, Metabolism) and p.profile is not None:
            if p.profile not in self.profiles:
                raise ValueError(
                    f"Process '{p.sid}' has the profile '{p.profile}', which the "
                    f"zonated liver does not give; its profiles are "
                    f"{sorted(self.profiles)}."
                )
            return p.profile
        if p.protein is not None and p.protein in self.profiles:
            return p.protein
        return None

    def _levels(self, key: str) -> tuple[float, ...]:
        """The relative abundance of a profile in every zone, periportal first.

        A sequence is taken as it is given, a `Profile` is evaluated by the
        engine at the number of zones (`pkpdsim.zonation.profiles`), which
        normalizes the levels to the mean 1 over the zones of equal volume.
        """
        profile = self.profiles[key]
        if isinstance(profile, Profile):
            return profile.values(len(self.zones))
        return tuple(float(v) for v in profile)

    def _weight(self, key: str | None, zone: str) -> str:
        """The share of a zone of the capacity of a process with a profile, by volume without.

        `w = p*V_zone/sum(p*V_zone) = p/sum(p)` with the relative abundance
        `p` of the profile in the zones (the zones have the same volume), so
        that the shares sum to 1; by volume the share is `f_zone`.
        """
        if key is None:
            return self._zone_fraction()
        token = key.lower()
        weight = ids.placement_weight(token, zone)
        if weight in self.parameters:
            return weight
        abundances = [ids.placement_profile(token, z) for z in self.zones]
        for z, pid, value in zip(
            self.zones, abundances, self._levels(key), strict=True
        ):
            self._constant(
                pid,
                value,
                U.dimensionless,
                f"relative abundance {key} zone {z}",
                f"Relative abundance of the profile `{key}` in the zone, given "
                "with the zonated liver (periportal zone first). A value may be "
                "changed (the shares `p/sum(p)` follow and keep the capacity of "
                "the liver), but not every value of the profile set to 0: the "
                "shares are then undefined.",
                kind="placement_weight",
            )
        total = " + ".join(abundances)
        for z, pid in zip(self.zones, abundances, strict=True):
            self._initial(
                ids.placement_weight(token, z),
                f"{pid} / ({total})",
                U.dimensionless,
                f"share {key} zone {z}",
                f"Share of the zone of the capacity of a process of the profile "
                f"`{key}`: its relative abundance times the zone volume over the "
                "sum of the zones (the zones have the same volume); the shares "
                "sum to 1, the capacity of the liver is kept.",
                kind="placement_weight",
            )
        return weight

    def _chemistry(self) -> None:
        """The processes of the liver in every zone, each with its share of the capacity.

        A metabolism and a transport take the share of the capacity of their
        zone; a custom reaction is placed in every zone with its own rate law
        over the ids of the zone, and the constants of the rate law and the
        species it brings itself are created once for the liver.
        """
        self._proteins_and_sinks([*self.processes, *self.customs])
        capacity = self._capacity()
        assert capacity is not None
        for p in self.processes:
            for zone in self.zones:
                share = f"{capacity} * {self._weight(self.profile_of[p.sid], zone)}"
                for reaction in place(p, self.ctx, self.ctx.parts, share, zone):
                    for par in reaction.pars:
                        self.parameters.setdefault(str(par.sid), par)
                    reaction.pars = []
                    self.reactions.append(reaction)
        for p in self.customs:
            self._custom_species(p)
            for k, zone in enumerate(self.zones):
                for reaction in place(p, self.ctx, self.ctx.parts, None, zone):
                    if k:
                        # the constants of the rate law belong to the liver,
                        # not to a zone: they are created with the reaction
                        # of the first zone and read by the others
                        reaction.pars = []
                    self.reactions.append(reaction)

    def _custom_species(self, p: CustomReaction) -> None:
        """The species a custom reaction brings itself, one per zone.

        Raises:
            ValueError: if a species of the zone has the id of a species of
                the liver already.
        """
        for zone in self.zones:
            mapping = zone_ids(self.ctx, self.ctx.parts, zone)
            for sp in p.species:
                sid = ids.segment_process(str(sp.sid), zone)
                if sid in self.species:
                    raise ValueError(
                        f"Custom reaction '{p.sid}': species '{sid}' is a species "
                        f"of the zone {zone} of the liver already."
                    )
                species = copy.copy(sp)
                species.sid = sid
                compartment = str(sp.compartment)
                species.compartment = mapping.get(compartment, compartment)
                self.species[sid] = annotate(
                    species,
                    "custom_amount" if sp.hasOnlySubstanceUnits else "custom_species",
                )

    def _custom_placed(
        self, p: CustomReaction, placed: set[str]
    ) -> frozenset[str] | None:
        """The ids of the reactions of a custom process, one per zone."""
        reactions = frozenset(ids.segment_process(p.sid, zone) for zone in self.zones)
        return reactions if reactions <= placed else None

    def _custom_read(self, p: CustomReaction) -> set[str]:
        """The ids which the reaction of a custom process reads in every zone.

        The ids of its rate law and of its declared steady state, each the id
        of the zone where the liver has one (`zone_ids`), so that the inputs
        of the ids of the top model are those of the rewritten rate laws.
        """
        read = formula_ids(p.formula) | stationary_ids(p)
        out: set[str] = set()
        for zone in self.zones:
            mapping = zone_ids(
                self.ctx, self.ctx.parts, zone, [str(sp.sid) for sp in p.species]
            )
            out |= {mapping.get(i, i) for i in read}
        return out

    def _radius(self) -> str:
        """The radius of the lumen of a sinusoid (m)."""
        return self._constant(
            zids.RADIUS,
            self._geometry().r_sin,
            U.meter,
            "radius of the sinusoid [m]",
            f"Radius of the lumen of a sinusoid, {_TABLE}. The lateral area of "
            "the zone follows it.",
            kind="sinusoid_geometry",
        )

    def _length(self) -> str:
        """The length of a sinusoid (m), which the zones divide."""
        return self._constant(
            zids.LENGTH,
            self._geometry().L_sin,
            U.meter,
            "length of the sinusoid [m]",
            f"Length of a sinusoid, {_TABLE}; the zones divide it, so two zones "
            "are `L_sin/n` apart.",
            kind="sinusoid_geometry",
        )

    def _porosity(self, zone: str) -> str:
        """The fenestral porosity of the endothelium of a zone (fraction).

        The gradient of the table over the positions of the zones (the
        engine, `pkpdsim.zonation.geometry.porosities`), which rises towards
        the pericentral end (0.06 to 0.08 for the human, Braet and Wisse
        2002): the pericentral endothelium is the more open one.
        """
        values = porosities(len(self.zones), self._geometry())
        return self._constant(
            zids.porosity(zone),
            values[self.zones.index(zone)],
            U.dimensionless,
            in_zone("fenestral porosity of the endothelium", zone),
            "Fraction of the endothelium of the zone which the fenestrae open, "
            f"{_TABLE}: the gradient at the position of the zone, which rises "
            "from the periportal to the pericentral end. It scales the "
            "fenestral exchange of the zone.",
            kind="sinusoid_geometry",
        )

    def _lateral_area(self, zone: str) -> str:
        """The lateral area between the plasma and the space of Disse of a zone (m**2).

        `2*V_lumen/r_sin`, the surface of the cylinder of the lumen of the
        zone (the engine, `ZonationBuilder._radial_assignments`), from the
        vascular volume of the zone, so that it follows the liver volume of
        the individual.

        The lumen of a zone of this model is the whole vascular space of the
        liver (`f_vas`, the sinusoids and the larger vessels), where the
        engine takes the sinusoidal lumen only, so the area and with it
        `PS_vas_<zone>_<s>` are about 2.1 times the engine's for the same
        liver (the factor of `_reference_velocity`, which reports the lower
        velocity of the same volume). It changes no behaviour, since neither
        conductance limits the transfer of a small molecule.
        """
        lumen = self._share(f"{self.v} * {self.f_vas}", zone)
        return self._initial(
            zids.area("pd", zone),
            f"2 dimensionless * {lumen} / {self._radius()} / 1000 l_per_m3",
            U.m2,
            in_zone("lateral area plasma and interstitium [m2]", zone),
            "Lateral area between the plasma and the space of Disse of the zone: "
            "the surface of the cylinder of its lumen, `2*V_lumen/r_sin` with the "
            "vascular volume of the zone (plasma and blood cells).",
            kind="sinusoid_geometry",
        )

    def _fenestral(self, s: Substance) -> str:
        """The permeability of a substance through an open fenestra (m/s).

        `d_aq/t_endo * renkin(r/r_fen)` of the engine
        (`pkpdsim.zonation.transport.fenestral_permeability`) at the porosity
        1, which the porosity of the zone scales; a model input, to be raised
        where the endothelium is known not to limit the substance.

        Raises:
            ValueError: for a molecule which no fenestra lets through (its
                radius reaches the radius of a fenestra), which would cut it
                off from the tissue of the liver.
        """
        geo = self._geometry()
        radius = molecular_radius(s.mr)
        p = fenestral_permeability(diffusion_of(s), radius, geo, 1.0)
        if p.value <= 0.0:
            raise ValueError(
                f"{s.sid}: the molecular radius {radius.value:.3g} m reaches the "
                f"radius of a fenestra ({geo.r_fen:.3g} m), so no fenestra lets "
                "the substance through and the sinusoid geometry of the zonated "
                "liver would cut it off from the tissue. It needs a transcellular "
                "path through the endothelium, which the liver does not have "
                "(it has no endothelial compartment); the geometry "
                "`Geometry.WELL_STIRRED` holds such a molecule, across the "
                "generic endothelium of the standard organ."
            )
        return self._constant(
            zids.permeability("fen", s.sid),
            p.value,
            SinusoidUnits.m_per_s,
            f"fenestral permeability {s.sid} [m/s]",
            f"Permeability of {s.sid} through an open fenestra: its diffusion "
            "coefficient in water over the thickness of the endothelium times "
            "the hindrance of its molecular radius "
            f"({radius.value:.3g} m) in a fenestra ({geo.r_fen:.3g} m, "
            f"Renkin 1954); {_TABLE}.",
            kind="pore_permeability",
        )

    def _ps_vas(self, s: Substance, zone: str | None) -> str:
        """The exchange between the plasma and the interstitium of a zone (l/min).

        The well-stirred zone exchanges across the generic endothelium of the
        standard organ, the one `PS_vas` of the organ whose share it takes;
        the sinusoid exchanges through the fenestrae of its zone, the
        conductance of the engine (`ZonationBuilder.radial_exchange`) over
        the lateral area of the lumen of the zone, as the named parameter
        `PS_vas_<zone>_<s> = P_fen_<s>*eps_fen_<zone>*A_pd_<zone>` of the
        zone and the substance: a conductance a sensitivity analysis scans
        and a fit changes is a parameter of the model and not a formula in a
        kinetic law. It is a constant with an initial assignment, as
        `PS_vas` of the standard organ is, so it follows the liver volume of
        the individual. The transcellular path of the endothelium is no path
        of the PBPK liver, which has no endothelial compartment.

        It is the conductance of the engine over the lateral area of this
        model, which is that of the whole vascular space of the zone and
        about 2.1 times the engine's sinusoidal wall (`_lateral_area`). For
        a small molecule neither number limits anything: both are hundreds
        of times the plasma flow through the chain, and the basolateral
        `PS_<s>` of the substance is the limiting step.
        """
        if not self._is_sinusoid:
            return super()._ps_vas(s, zone)
        assert zone is not None
        return self._initial(
            ids.ps_vas(s.sid, zone),
            f"{self._fenestral(s)} * {self._porosity(zone)} * "
            f"{self._lateral_area(zone)} * 1000 l_per_m3 * 60 s_per_min",
            U.l_per_min,
            in_zone(f"exchange {s.sid} plasma and interstitium [l/min]", zone),
            f"Exchange of {s.sid} between the plasma of the sinusoid and the "
            "space of Disse of the zone, through the fenestrae of its "
            "endothelium: its permeability through an open fenestra times the "
            "fenestral porosity of the zone and the lateral area of its lumen.",
            kind="permeability",
        )

    def _diffusion(self, s: Substance) -> str:
        """The diffusion coefficient of a substance in water (m**2/s).

        The coefficient of the engine
        (`pkpdsim.zonation.transport.diffusion_parameter`): the value given
        with the substance or the one derived from its molar mass, with its
        provenance in the notes.
        """
        pid = zids.diffusion(s.sid)
        if pid not in self.parameters:
            parameter = diffusion_parameter(s, self._peclet_notes(s))
            self.parameters[pid] = self._annotate(parameter, "diffusion")
        return pid

    def _peclet_notes(self, s: Substance) -> str:
        """What the chain of zones does to the axial diffusion of a substance.

        The Peclet number of the chain with this coefficient at the
        reference male of the species, against the numerical dispersion
        of the upwind convection: where the dispersion exceeds the
        coefficient, the zone count and not the physics sets the gradient
        along the sinusoid (the warning of the engine,
        `ZonationBuilder._warn_dispersion`).
        """
        n = len(self.zones)
        geo = self._geometry()
        v = self._reference_velocity()
        numerical = dispersion(v, geo.L_sin, n)
        d = float(diffusion_of(s).value)
        above = "above" if d > numerical else "below"
        return (
            f"At the reference male of the species the chain of {n} zones "
            f"realizes the Peclet number {peclet(v, geo.L_sin, n, d):.3g} with "
            f"this coefficient, against {peclet(v, geo.L_sin, n, 0.0):.3g} = 2n "
            f"for a substance which does not diffuse: the coefficient lies "
            f"{above} the numerical dispersion {numerical:.3g} m**2/s "
            "(`v*L_sin/(2n)`) of the upwind chain, which the zone count and not "
            "the physics sets."
        )

    def _substance(self, s: Substance) -> None:
        """The species, flows and exchanges of a substance, its diffusion and its means."""
        super()._substance(s)
        if self._is_sinusoid and self.axial:
            self._axial_diffusion(s)
        if self.means:
            self._zonal_means(s)

    def _zonal_means(self, s: Substance) -> None:
        """The periportal and the pericentral mean of a substance in every zonated part.

        The observers of the engine (`pkpdsim.zonation.observers`): the mean
        of the concentration over the first (periportal) and the last
        (pericentral) third of the sinusoid, every zone weighted by its
        overlap with the third
        (`pkpdsim.zones.third_weights`, the weighting every zonated model of
        pkpdsim shares, so a zonated model of the engine and one of this
        library report the same two means). They read the species
        of the zones and nothing reads them, so they change no state of the
        model.
        """
        n = len(self.zones)
        for part in built_parts(self.ctx.parts[s.sid]):
            of_zones = [ids.species(part, s.sid, zone) for zone in self.zones]
            for periportal in (True, False):
                third = "periportal" if periportal else "pericentral"
                self._rule(
                    ids.zonal_mean(part, s.sid, periportal),
                    weighted_mean(of_zones, third_weights(n, periportal)),
                    U.mM,
                    f"{third} mean {s.sid} {_PART_NAMES[part]}",
                    f"Mean concentration of {s.sid} in the {_PART_NAMES[part]} of "
                    f"the {third} third of the sinusoid, the zones of the third "
                    "weighted by their overlap with it (the observers of the "
                    "zonation engine, `pkpdsim.zonation.observers`). It reads the "
                    "zones and nothing reads it.",
                    kind="organ_concentration",
                )

    def _axial_diffusion(self, s: Substance) -> None:
        """The diffusion of a substance along the sinusoid, between consecutive zones.

        `D_<s>*V<part>_<zone>*n^2/L_sin^2*(C_<zone> - C_<next>)` in the
        plasma and in the interstitium (the space of Disse), the reactions of
        the engine (`ZonationBuilder.axial_diffusion`) over the parts of the
        liver: the cross-section of a part is its volume in all zones over
        the length of the sinusoid and two zones are `L_sin/n` apart. None
        with one zone, and none in a part the substance does not have.
        """
        n = len(self.zones)
        if n < 2:
            return
        parts = [p for p in built_parts(self.ctx.parts[s.sid]) if p in _DIFFUSING]
        if not parts:
            return
        d = self._diffusion(s)
        length = self._length()
        for part in parts:
            for k, (upstream, downstream) in enumerate(
                itertools.pairwise(self.zones), start=1
            ):
                a = ids.species(part, s.sid, upstream)
                b = ids.species(part, s.sid, downstream)
                self._reaction(
                    ids.axial_reaction(part, s.sid, upstream),
                    f"{a} <-> {b}",
                    f"{d} * {ids.compartment(part, upstream)} * {n * n} "
                    f"dimensionless / ({length} * {length}) * 60 s_per_min * "
                    f"({a} - {b})",
                    f"axial diffusion {s.sid} {_DIFFUSING[part]}, zones {k} "
                    f"and {k + 1}",
                    kind="exchange",
                )

    def _unbound_in_cell(self, sid: str, protein: str) -> str:
        """The unbound concentration in the cells which drives the turnover of a protein.

        One turnover state `E_<protein>` for the liver, driven by the
        unbound concentration in the cells of the zones (intracellular, else
        the equivalent tissue part, else plasma) weighted by the placement of
        the protein, `sum w_zone*Cu_zone`.
        """
        s = next(s for s in self.ctx.substances if s.sid == sid)
        key = self.protein_profile.get(protein)
        terms = [
            f"{self._weight(key, zone)} * "
            + inhibitor_unbound(Part.INTRACELLULAR, s, self.organ, self.ctx.parts, zone)
            for zone in self.zones
        ]
        return f"({' + '.join(terms)})"
