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
custom reaction states its rate for the whole liver and is not split over
the zones: the zonated liver refuses it.

The zones are well-stirred tanks in series. With one zone the zonated liver
is the standard liver. With `n` zones and linear kinetics the availability
of a flow-limited liver at steady state is `(1 + x/n)^-n` with
`x = fu*CLint/(Q*BP)` instead of the well-stirred `1/(1 + x)`: with the same
intrinsic clearance the zonated liver extracts more, the difference of the
extraction is about `x/2*(1 - 1/n)` of it (second order in `x`), and for
large `n` it approaches the parallel tube, `exp(-x)`.
"""

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

from sbmlutils.factory import ModelDefinition

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.chemistry import CustomReaction, Metabolism, Transport
from pkpdsim.pbpk.organ import Boundary, OrganContext, StandardBuilder
from pkpdsim.pbpk.organs import Organ, Part
from pkpdsim.pbpk.placement import inhibitor_unbound, place
from pkpdsim.pbpk.units import U

__all__ = ["ZonatedLiver"]


@dataclass(frozen=True)
class ZonatedLiver:
    """The liver as zones in series along the sinusoid.

    Attributes:
        n_zones: the number of zones, periportal to pericentral.
        profiles: the relative abundance of a protein in each zone, from the
            periportal to the pericentral zone (`n_zones` values, not all 0),
            by the id of the protein or the `profile` of a metabolism; the
            processes of a protein without a profile are placed by volume.
        name: the name of the component, the prefix of the id of the model
            definition (`zonated_li`).
    """

    n_zones: int = 8
    profiles: Mapping[str, Sequence[float]] = field(default_factory=dict)
    name: str = "zonated"

    def __post_init__(self) -> None:
        """Check the zones and the profiles.

        Raises:
            ValueError: for less than one zone, a profile key which is no
                id, two keys which differ in case only, or a profile without
                a value per zone, with a negative value or only zeros.
        """
        if isinstance(self.n_zones, bool) or self.n_zones < 1:
            raise ValueError(
                f"The zonated liver needs at least one zone, not {self.n_zones}."
            )
        lower: dict[str, str] = {}
        for key, values in self.profiles.items():
            if not key.isidentifier():
                raise ValueError(f"Profile '{key}': the key must be an id.")
            if key.lower() in lower:
                raise ValueError(
                    f"Profiles '{lower[key.lower()]}' and '{key}' have the same id."
                )
            lower[key.lower()] = key
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
        object.__setattr__(
            self,
            "profiles",
            {
                key: tuple(float(v) for v in values)
                for key, values in self.profiles.items()
            },
        )

    def build(self, ctx: OrganContext) -> tuple[ModelDefinition, Boundary]:
        """Build the model of the liver and its boundary.

        Args:
            ctx: the context of the liver.

        Returns:
            The model definition and the boundary of the standard liver.

        Raises:
            ValueError: if the organ is not the liver, a process cannot be
                placed (a custom reaction, a metabolism with a profile the
                liver does not give), the processes of a protein have
                different profiles, a profile is used by no process, or the
                standard liver refuses a substance.
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
            ValueError: for another organ than the liver, or a process the
                zonated liver cannot place.
        """
        if ctx.organ is not Organ.LIVER:
            raise ValueError(
                f"The zonated liver is a variant of the liver, not of organ "
                f"'{ctx.organ.value}'."
            )
        super().__init__(ctx, liver.name)
        self.zones = tuple(ids.zone(i) for i in range(1, liver.n_zones + 1))
        self.profiles = liver.profiles
        self.processes, self.profile_of = self._placed()
        self.protein_profile = {
            p.protein: self.profile_of[p.sid]
            for p in self.processes
            if p.protein is not None
        }

    def _notes(self) -> str:
        """The notes of the model definition."""
        n = len(self.zones)
        return (
            f"The {self.name} liver: {n} zones in series along the sinusoid, from "
            f"the periportal zone {self.zones[0]} to the pericentral zone "
            f"{self.zones[-1]}, each with the parts which the substances need and "
            "its share of the volumes, the exchanges and the chemistry of the "
            "liver; the boundary of the standard liver."
        )

    def _placed(
        self,
    ) -> tuple[list[Metabolism | Transport], dict[str, str | None]]:
        """The processes of the liver and their profiles (`None`: by volume).

        Raises:
            ValueError: for a custom reaction, a metabolism with a profile
                the liver does not give, the processes of a protein with
                different profiles, or a profile which no process uses.
        """
        processes: list[Metabolism | Transport] = []
        profile_of: dict[str, str | None] = {}
        for p in self._processes():
            if isinstance(p, CustomReaction):
                raise ValueError(
                    f"Custom reaction '{p.sid}' of the liver: the zonated liver "
                    "places metabolism and transport only, a custom reaction "
                    "states its rate for the whole liver."
                )
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
            self.zones, abundances, self.profiles[key], strict=True
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
            )
        return weight

    def _chemistry(self) -> None:
        """The processes of the liver in every zone, each with its share of the capacity."""
        self._proteins_and_sinks(self.processes)
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
