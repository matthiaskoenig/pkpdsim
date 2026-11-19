"""The lobed liver: the liver as lobes perfused in parallel.

The liver of the standard organ is one chain fed by the hepatic artery and the
portal vein. The lobed liver is one chain per lobe: each lobe holds its share
of the liver volume (`f_vol_<lobe>`), is fed at its own flows (`Qha_<lobe>`,
`Qpo_<lobe>`) and drains at its own concentration (`Cout_<lobe>_<s>`), and the
body receives the flow weighted mix of the lobes under the unchanged id
`Cout_<s>`. The whole organ flows `Qha`, `Qpo` and `Qh` are the sums of the
lobe flows, so every consumer of them (the observers, the steady state, the PK
and the pressure driven circulation) reads what it read before.

The lobe and the zone are the two spatial axes of the liver and they compose:
the lobes are perfused in parallel, the zones of a lobe in series along its
sinusoid, and a site (`ids.site`) carries both. The single lobe liver is the
liver without lobes, exactly as the single zone liver is the liver without
zones: it emits no lobe token and no lobe parameter, so a model which asks for
one lobe is the standard liver.

Each lobe is supplied in proportion to its reference supply share, which the
lobe table gives or which is the share of its mass, renormalized over the
patent lobes: a lobe whose branch is closed (`g_po_<lobe>` 0, portal vein
ligation) passes its flow to the lobes which are open, which is what happens
when the splanchnic bed keeps delivering its blood to a liver with one branch
less. With the pressure driven circulation the flows come from per lobe
vascular branches instead, through the same ids.
"""

from collections.abc import Sequence
from dataclasses import dataclass, replace
from typing import TYPE_CHECKING

from sbmlutils.factory import (
    AssignmentRule,
    InitialAssignment,
    ModelDefinition,
    Parameter,
    Sbase,
)
from sbmlutils.metadata import SBO

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.annotation import annotate
from pkpdsim.pbpk.lobes import HABR, Lobe, Remodelling, check_lobes, lobe_table
from pkpdsim.pbpk.organ import Boundary, OrganContext, StandardOrgan
from pkpdsim.pbpk.organs import Organ
from pkpdsim.pbpk.placement import built_parts
from pkpdsim.pbpk.remodelling import VREMODELLED, remodelling_sink, sink_objects
from pkpdsim.pbpk.substance import Substance
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.variants.liver_zonated import ZonatedBuilder, ZonatedLiver
from pkpdsim.sbml.parameters import value_parameter
from pkpdsim.species import Species

if TYPE_CHECKING:
    from pkpdsim.pbpk.definition import PBPKDefinition

__all__ = ["MAX_SITES", "LobedBuilder", "LobedLiver"]

MAX_SITES = 40
"""How many sites (lobes times zones) a liver may have.

A model choice: a site carries the parts of the liver for every substance, and
the tests of the PBPK library compile their models at `O0`, so a liver of 40
sites with the six layers of the sinusoid geometry is the largest the suite
builds in its budget. The eight Couinaud segments therefore take at most five
zones each; raise it knowingly.
"""

_SBO_PARAMETER = SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER
"""The SBO term of a parameter of the body."""


@dataclass(frozen=True)
class LobedLiver:
    """The liver as lobes perfused in parallel.

    Attributes:
        lobes: the lobes, `None` for the set of the species
            (`pkpdsim.pbpk.lobes.lobe_table`: the anatomical lobes of rat and
            mouse, the eight Couinaud segments of the human). The shares of
            the mass sum to 1, the tokens are unique lower case identifiers.
        zonation: the zonation of every lobe, `None` for lumped lobes. Each
            lobe is then the zonated liver of this setting, with its share of
            the volume and of the flow.
        remodelling: the driver of the lobe volumes during a simulation,
            `None` for a constant anatomy. A deportalised lobe atrophies and
            its neighbours hypertrophy; the driver says how
            (`pkpdsim.pbpk.remodelling`).
        habr: the hepatic arterial buffer response, `None` for arterial
            branches which do not react to a lost portal inflow.
        name: the name of the component, the prefix of the id of the model
            definition (`lobes_li`).
    """

    lobes: Sequence[Lobe] | None = None
    zonation: ZonatedLiver | None = None
    remodelling: Remodelling | None = None
    habr: HABR | None = None
    name: str = "lobes"

    def __post_init__(self) -> None:
        """Check the lobes which are given and the number of sites.

        Raises:
            ValueError: for a lobe set which is no set
                (`pkpdsim.pbpk.lobes.check_lobes`), or for more sites than
                `MAX_SITES`.
        """
        if self.lobes is not None:
            object.__setattr__(self, "lobes", check_lobes(self.lobes))
            _check_sites(len(self.lobes), self._n_zones)

    @property
    def _n_zones(self) -> int:
        """The number of zones of every lobe, 1 for lumped lobes."""
        return 1 if self.zonation is None else self.zonation.n_zones

    def lobe_set(self, species: Species) -> tuple[Lobe, ...]:
        """The lobes of the component, the set of a species where it gives none.

        Args:
            species: the species of the build.

        Returns:
            The lobes.

        Raises:
            ValueError: for a species without a lobe set, or for more sites
                than `MAX_SITES`.
        """
        if self.lobes is not None:
            return tuple(self.lobes)
        lobes = lobe_table(species)
        _check_sites(len(lobes), self._n_zones)
        return lobes

    def resolved(self, ctx: OrganContext) -> tuple[Lobe, ...]:
        """The lobes of the component in the context of a build."""
        return self.lobe_set(ctx.species)

    def body_objects(self, definition: PBPKDefinition) -> list[Sbase]:
        """Get the anatomy, the supply and the flows of the lobes in the top model.

        Args:
            definition: the definition of the model.

        Returns:
            The objects of the top model, none for a liver of one lobe (it is
            the liver without lobes).

        Raises:
            ValueError: for a species without a lobe set.
        """
        lobes = self.lobe_set(definition.species)
        if len(lobes) < 2:
            return []
        objects = lobe_objects(lobes)
        if self.remodelling is not None:
            objects += self.remodelling.objects(lobes, definition.species)
            objects += sink_objects(definition)
        return objects

    def build(self, ctx: OrganContext) -> tuple[ModelDefinition, Boundary]:
        """Build the model of the lobed liver and its boundary.

        A liver of one lobe is the liver without lobes: the standard liver,
        or the zonated liver of the zonation of the component.

        Args:
            ctx: the context of the liver.

        Returns:
            The model definition and the boundary.

        Raises:
            ValueError: if the organ is not the liver, the species has no
                lobe set, there are more sites than `MAX_SITES`, or the
                liver refuses a substance or a process (the standard and the
                zonated liver).
        """
        if ctx.organ is not Organ.LIVER:
            raise ValueError(
                f"The lobed liver is a variant of the liver, not of organ "
                f"'{ctx.organ.value}'."
            )
        lobes = self.resolved(ctx)
        if len(lobes) < 2:
            if self.zonation is not None:
                return self.zonation.build(ctx)
            return StandardOrgan().build(ctx)
        return LobedBuilder(ctx, self, lobes).build()


def _check_sites(n_lobes: int, n_zones: int) -> None:
    """Check the number of sites of a lobed liver.

    Args:
        n_lobes: the number of lobes.
        n_zones: the number of zones of every lobe.

    Raises:
        ValueError: for more than `MAX_SITES` sites, naming the product and
            the limit.
    """
    sites = n_lobes * n_zones
    if sites > MAX_SITES:
        raise ValueError(
            f"The lobed liver would have {n_lobes} lobes of {n_zones} zones, "
            f"{sites} sites; a liver has at most {MAX_SITES} sites."
        )


def _parameter(
    pid: str,
    v: float,
    unit: object,
    name: str,
    text: str,
    kind: str,
    *,
    constant: bool = True,
) -> Parameter:
    """A parameter of the body of the lobes, annotated by its kind.

    `constant=False` for a parameter an assignment rule follows during the
    simulation (the perfusion of a lobe, the remnant fraction of the liver).
    """
    return annotate(
        Parameter(
            pid,
            v,
            unit,  # ty: ignore[invalid-argument-type]
            constant=constant,
            name=name,
            sboTerm=_SBO_PARAMETER,
            notes=text,
        ),
        kind,
        organ=Organ.LIVER,
    )


def _patency(lobe: Lobe) -> list[Sbase]:
    """The patency of the two inflow branches of a lobe and its presence."""
    return [
        _parameter(
            ids.arterial_conductance(lobe.token),
            1.0,
            U.dimensionless,
            f"patency of the hepatic artery branch of {lobe.name}",
            "1 is the open branch, 0 the ligated one. An occlusion acts on "
            "the conductance of a branch and never on its resistance, so no "
            "rate law divides by zero. An initial value; an intervention of "
            "the model sets it during a simulation.",
            "lobe_patency",
        ),
        _parameter(
            ids.portal_conductance(lobe.token),
            1.0,
            U.dimensionless,
            f"patency of the portal vein branch of {lobe.name}",
            "1 is the open branch, 0 the ligated one (portal vein ligation). "
            "The portal flow of a closed branch goes to the lobes which are "
            "open. An initial value; an intervention of the model sets it "
            "during a simulation.",
            "lobe_patency",
        ),
        _parameter(
            ids.lobe_present(lobe.token),
            1.0,
            U.dimensionless,
            f"presence of {lobe.name}",
            "1 is the lobe in place, 0 the resected one, which holds no "
            "substance and has no function.",
            "lobe_present",
        ),
    ]


def _shares(lobe: Lobe) -> list[Sbase]:
    """The reference shares of the volume and of the two inflows of a lobe."""
    objects: list[Sbase] = [
        annotate(
            value_parameter(
                ids.lobe_reference_fraction(lobe.token),
                lobe.mass_fraction,
                U.dimensionless,
                f"reference volume fraction of {lobe.name}",
                _SBO_PARAMETER,
                "The share of the volume of the intact liver which this lobe "
                "holds, the constant of the lobe table; the share which it "
                "holds during a simulation follows its volume.",
            ),
            "lobe_fraction",
            organ=Organ.LIVER,
        ),
        _parameter(
            ids.lobe_volume_fraction(lobe.token),
            lobe.mass_fraction.value,
            U.dimensionless,
            f"volume fraction of {lobe.name}",
            "The share of the volume of the liver which this lobe holds, "
            "`Vli_<lobe>/Vli` against the reference volume: it follows the "
            "volume of the lobe when the anatomy of the liver changes during "
            "a simulation, and the shares then sum to the remnant fraction "
            "of the liver instead of 1.",
            "lobe_fraction",
        ),
        InitialAssignment(
            ids.lobe_volume_fraction(lobe.token),
            f"{ids.lobe_volume(lobe.token)} / {ids.organ_volume(Organ.LIVER)}",
            U.dimensionless,
        ),
    ]
    supply = (
        (
            ids.arterial_fraction,
            lobe.arterial_fraction,
            lobe.arterial_share,
            "hepatic artery",
        ),
        (ids.portal_fraction, lobe.portal_fraction, lobe.portal_share, "portal vein"),
    )
    for name, given, share, what in supply:
        text = (
            f"The share of the {what} inflow of the liver which this lobe "
            "receives, relative to the other lobes."
        )
        if given is None:
            text += (
                " No source gives it, so it is the share of the mass of the "
                "lobe: a Couinaud segment is defined as a portal territory "
                "and the hepatic artery branch runs with the portal branch "
                "in the portal triad, so the perfusion per volume is the "
                "same in every lobe."
            )
        objects.append(
            annotate(
                value_parameter(
                    name(lobe.token),
                    given if given is not None else share(),
                    U.dimensionless,
                    f"{what} fraction of {lobe.name}",
                    _SBO_PARAMETER,
                    text,
                ),
                "lobe_supply_fraction",
                organ=Organ.LIVER,
            )
        )
    return objects


def _flows(lobes: Sequence[Lobe]) -> list[Sbase]:
    """The blood flows of the lobes and their perfusion per volume.

    The inflow of a lobe is its share of the inflow of the liver among the
    lobes whose branch is open: `Qpo*g_po*f_po/sum(g_po*f_po)`. An occlusion
    therefore passes the flow of the closed branch to the lobes which are
    open, with no pressure model present, and the sum of the lobe flows is
    the flow of the liver whatever is open.
    """
    objects: list[Sbase] = []
    branches = (
        (ids.QHA, ids.arterial_conductance, ids.arterial_fraction),
        (ids.QPO, ids.portal_conductance, ids.portal_fraction),
    )
    for flow, conductance, fraction in branches:
        total = " + ".join(
            f"{conductance(lobe.token)} * {fraction(lobe.token)}" for lobe in lobes
        )
        for lobe in lobes:
            objects.append(
                InitialAssignment(
                    ids.lobe_flow(flow, lobe.token),
                    f"{flow} * {conductance(lobe.token)} * "
                    f"{fraction(lobe.token)} / ({total})",
                    U.l_per_min,
                )
            )
    for lobe in lobes:
        objects.append(
            InitialAssignment(
                ids.lobe_flow(ids.QH, lobe.token),
                f"{ids.lobe_flow(ids.QHA, lobe.token)} + "
                f"{ids.lobe_flow(ids.QPO, lobe.token)}",
                U.l_per_min,
            )
        )
    return objects


def lobe_objects(lobes: Sequence[Lobe]) -> list[Sbase]:
    """Get the objects of the top model of a lobe set.

    The shares of the volume and of the inflows, the patency of the branches,
    the volume, the three blood flows and the perfusion per volume of every
    lobe, and the remnant fraction of the liver.

    Args:
        lobes: the lobes, a checked set.

    Returns:
        The objects.
    """
    objects: list[Sbase] = []
    for lobe in lobes:
        objects += _shares(lobe)
        objects += _patency(lobe)
    objects += [
        _parameter(
            ids.VLI_REF,
            1.0,
            U.litre,
            "reference volume of the liver [l]",
            "The volume of the intact liver of the individual, which the "
            "shares of the lobes are measured against and which a remnant "
            "grows back towards. It is taken at the start of a simulation "
            "and does not change during it.",
            "lobe_volume",
        ),
        InitialAssignment(ids.VLI_REF, ids.organ_volume(Organ.LIVER), U.litre),
        _parameter(
            ids.VLI_LOBES,
            1.0,
            U.litre,
            "volume of the lobes of the liver [l]",
            "The sum of the volumes of the lobes: the liver which is there, "
            "the reference volume of the intact organ and less after a "
            "resection or an atrophy.",
            "lobe_volume",
            constant=False,
        ),
        AssignmentRule(
            ids.VLI_LOBES,
            " + ".join(ids.lobe_volume(lobe.token) for lobe in lobes),
            U.litre,
        ),
    ]
    for lobe in lobes:
        objects.append(
            _parameter(
                ids.lobe_volume(lobe.token),
                1.0,
                U.litre,
                f"volume of {lobe.name} [l]",
                "The volume of the lobe; the volumes of the lobes sum to the "
                "volume of the liver.",
                "lobe_volume",
            )
        )
        objects.append(
            InitialAssignment(
                ids.lobe_volume(lobe.token),
                f"{ids.organ_volume(Organ.LIVER)} * "
                f"{ids.lobe_reference_fraction(lobe.token)}",
                U.litre,
            )
        )
    for flow in ids.HEPATIC_FLOWS:
        for lobe in lobes:
            objects.append(
                _parameter(
                    ids.lobe_flow(flow, lobe.token),
                    1.0,
                    U.l_per_min,
                    f"blood flow {flow} of {lobe.name} [l/min]",
                    "The share of the blood flow of the liver which this "
                    "lobe receives; the lobe flows sum to the flow of the "
                    "liver.",
                    "lobe_flow",
                )
            )
    objects += _flows(lobes)
    for lobe in lobes:
        objects.append(
            _parameter(
                ids.lobe_perfusion(lobe.token),
                1.0,
                U.per_min,
                f"perfusion per volume of {lobe.name} [1/min]",
                "The total blood flow of the lobe per volume of the lobe: "
                "the quantity a lobe is hyperperfused or hypoperfused in, "
                "and the stimulus of the remodelling.",
                "lobe_perfusion",
                constant=False,
            )
        )
        objects.append(
            AssignmentRule(
                ids.lobe_perfusion(lobe.token),
                f"{ids.lobe_flow(ids.QH, lobe.token)} / {ids.lobe_volume(lobe.token)}",
                U.per_min,
            )
        )
    objects.append(
        _parameter(
            ids.F_REMNANT,
            1.0,
            U.dimensionless,
            "remnant fraction of the liver",
            "The share of the liver which is left: 1 is the intact liver, "
            "and a resection of the lobes of half the mass leaves 0.5. It is "
            "the quantity the hepatectomy studies report as the remnant.",
            "remnant_fraction",
            constant=False,
        )
    )
    objects.append(
        AssignmentRule(
            ids.F_REMNANT,
            f"{ids.VLI_LOBES} / {ids.VLI_REF}",
            U.dimensionless,
        )
    )
    return objects


class LobedBuilder(ZonatedBuilder):
    """The SBML objects of the lobed liver: the chains of the lobes in parallel.

    It is the builder of the zonated liver with a second axis: the lobes are
    the parallel chains, the zones the steps along each chain, and a site
    (`ids.site`) carries both. A lumped lobed liver has no zone token, as the
    zonated liver of one lobe has no lobe token.
    """

    def __init__(
        self, ctx: OrganContext, liver: LobedLiver, lobes: Sequence[Lobe]
    ) -> None:
        """Start the collection of the elements of the lobed liver.

        Args:
            ctx: the context of the liver.
            liver: the component.
            lobes: the resolved lobes.

        Raises:
            ValueError: for another organ than the liver, or whatever the
                zonated liver refuses (a process it cannot place, a species
                without the geometry of the sinusoid).
        """
        zonation = liver.zonation if liver.zonation is not None else ZonatedLiver()
        super().__init__(ctx, replace(zonation, name=liver.name))
        if liver.zonation is None:
            # the lumped lobes are the chains of one site, which carries the
            # lobe and no zone, as the liver of one zone carries none
            self.zones = ()
        self.liver = liver
        self.lobe_set = tuple(lobes)
        self.lobes = tuple(lobe.token for lobe in lobes)

    def _notes(self) -> str:
        """The notes of the model definition."""
        where = (
            "lumped"
            if not self.zones
            else f"in the {len(self.zones)} zones of its sinusoid"
        )
        return (
            f"The liver in the lobes {list(self.lobes)}, perfused in parallel, "
            f"each {where}: a lobe holds its share of the volume and of the "
            "capacity of the liver, is fed at its own hepatic artery and "
            "portal vein flows and drains at its own concentration; the organ "
            "drains the flow weighted mix of the lobes."
        )

    def _substance(self, s: Substance) -> None:
        """The species, flows and exchanges of a substance, and what atrophy carries away."""
        super()._substance(s)
        if self.liver.remodelling is not None:
            self._remodelled(s)

    def _remodelled(self, s: Substance) -> None:
        """The amount of a substance which the tissue of a shrinking lobe carries away.

        The state of a species of the liver is its amount, so a lobe which
        grows dilutes what it holds by itself; a lobe which shrinks would
        concentrate it instead, which no tissue does. Every part of every
        site of a lobe therefore loses the share of its amount which the
        lobe loses, `A*max(0, -dVli_<lobe>)/Vli_<lobe>`, into the sink of
        the body: the concentration of a shrinking lobe is left alone and
        the amount which left the body is auditable.
        """
        self._input_compartment(VREMODELLED, "remodelled tissue")
        sink = self._input_species(
            remodelling_sink(s.sid),
            VREMODELLED,
            f"{s.name} (remodelled tissue)",
            amount=True,
        )
        for lobe, sites in self._chains:
            assert lobe is not None
            loss = (
                f"max(0 l_per_min, -{self._growth(lobe)}) / {self._lobe_volume(lobe)}"
            )
            for site in sites:
                for part in built_parts(self.ctx.parts[s.sid]):
                    conc = ids.species(part, s.sid, site)
                    self._reaction(
                        ids.remodelling_reaction(part, s.sid, str(site)),
                        f"{conc} -> {sink}",
                        f"{ids.compartment(part, site)} * {conc} * {loss}",
                        f"loss {s.sid} {part.value} of lobe {lobe} to atrophy",
                        kind="remodelling",
                    )

    def _growth(self, lobe: str) -> str:
        """The growth of a lobe, the input `dVli_<lobe>` (l/min)."""
        return self._input(
            ids.lobe_growth(lobe),
            U.l_per_min,
            f"growth of lobe {lobe} [l/min]",
            "The rate of change of the volume of the lobe, which the "
            "remodelling of the body drives.",
        )

    def _lobe_volume(self, lobe: str) -> str:
        """The volume of a lobe, the input `Vli_<lobe>` (l)."""
        return self._input(
            ids.lobe_volume(lobe),
            U.litre,
            f"volume of lobe {lobe} [l]",
            "The volume of the lobe of the body.",
        )

    def _zone_of(self, site: str | None) -> str:
        """The zone of a site, the lobe token stripped.

        The porosity of the endothelium and the relative abundance of a
        profile are quantities of the zone and of no lobe, so the lobes share
        one set of them.
        """
        assert site is not None
        lobe = self._lobe_of(site)
        assert lobe is not None
        return site[len(lobe) + 1 :]
