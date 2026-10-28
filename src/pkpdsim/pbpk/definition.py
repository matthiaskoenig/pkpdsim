"""The definition of a PBPK model: the substances, the chemistry and the choices."""

import re
from collections.abc import Callable, Hashable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field

from sbmlutils.factory import Creator
from sbmlutils.factory.units import UnitDefinition, UnitType

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.chemistry import (
    MM,
    CustomReaction,
    FirstOrder,
    Hill,
    Metabolism,
    Transport,
)
from pkpdsim.pbpk.data import load
from pkpdsim.pbpk.kp import check_kp_method
from pkpdsim.pbpk.observers import Observer, formula_functions
from pkpdsim.pbpk.organ import OrganComponent
from pkpdsim.pbpk.organs import Organ, Part, Pool, as_member, as_organ_keys
from pkpdsim.pbpk.pd import PDComponent
from pkpdsim.pbpk.physiology import DEFAULT_POPULATION, Sex
from pkpdsim.pbpk.populations import POPULATIONS
from pkpdsim.pbpk.proteins import PROTEINS
from pkpdsim.pbpk.routes import Route
from pkpdsim.pbpk.substance import Erythrocytes, Kp, KpMethod, Mode, Substance
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Formula, Q, Value

__all__ = ["ModelParameter", "PBPKDefinition", "Variants"]


_SID = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
_RESERVED = (
    *(p.value for p in Part),
    *(p.value for p in Pool),
    "free",
    "blood",
    *ids.GUT_SEGMENTS,
)
"""The tokens of the id scheme which a substance id must not start with (the
segment tokens follow `Clumen_`, `Cmu_` and `Asolid_`)."""

_ZONE = re.compile(r"z\d+_")
"""The zone tokens (`z1_`) of the id scheme, which follow the part in the ids
of a zonated organ (`Cpl_z1_<s>`); a substance id must not start with one."""


def _physiology_ids() -> frozenset[str]:
    """Get all ids of the physiology which a substance id must not take."""
    names = {
        value
        for name, value in vars(ids).items()
        if name.isupper() and name != "__all__" and isinstance(value, str)
    }
    for organ in Organ:
        names |= {
            ids.organ_volume(organ),
            ids.organ_flow(organ),
            ids.organ_mass(organ),
            ids.flow_fraction(organ),
        }
    for pool in Pool:
        names |= {ids.pool_volume(pool), ids.pool_fraction(pool)}
    return frozenset(names)


@dataclass(frozen=True)
class ModelParameter:
    """A constant parameter of the top model which the model package states.

    A constant of the model which no field of the library holds (a
    production rate which the baseline of a substance and a custom reaction
    of the muscle share): the formulas of the definition, its observers and
    the custom reactions of the organs read it under its id, an organ through
    an input port under the same id.

    Attributes:
        sid: the id of the parameter, an id of the top model which is no id
            of the library.
        value: the number in `unit`, a `Value` (a fitted value with its
            uncertainty, `fitted(FITS.parameter(pid, unit=...))`) or a
            `Formula` over ids of the top model (its initial assignment); a
            `Value` or a `Formula` is in `unit`.
        unit: the unit of the parameter.
        name: the name of the parameter, the id by default.
        notes: the notes of the parameter, e.g. its meaning and source.
    """

    sid: str
    value: Q
    unit: UnitDefinition
    name: str | None = None
    notes: str | None = None


@dataclass(frozen=True)
class Variants:
    """The organ variants of a definition, the standard organ everywhere else.

    Attributes:
        organs: the variant of an organ, by name or as component; the names
            are resolved from the variant registry when the model is built.
            The organs may be given as `Organ` members, their values or their
            names (`"li"`, `"liver"`), normalised to `Organ` members.
    """

    organs: Mapping[Organ, str | OrganComponent] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalise the organs to `Organ` members."""
        object.__setattr__(
            self, "organs", as_organ_keys(self.organs, "Variants.organs")
        )


@dataclass(frozen=True)
class PBPKDefinition:
    """The definition of a PBPK model.

    Attributes:
        sid: the id of the model.
        substances: the substances.
        chemistry: the processes, each defined once per organ.
        variants: the organ variants.
        routes: the routes of administration.
        pd: the pharmacodynamic components.
        sex: the sex of the model, chosen when it is built: its reference
            values and default individual; a study of both sexes uses one
            model per sex. A `Sex` member or its value (`"female"`),
            normalised to the member.
        population: the population (`pkpdsim.pbpk.populations`), chosen when
            the model is built: its reference individual of the sex is the
            default individual (`BW`, `HEIGHT`, `AGE`); the reference anatomy
            stays ICRP 89.
        arm_vein: the preset of the weights of the arm vein sampling site
            (`data/arm_vein.toml`), the PK-Sim default skin and muscle mix.
        dynamic: the parameters of the top model which change during a
            simulation (an impairment changed in the middle of a study, a PD
            output): every quantity which depends on them is an assignment
            rule instead of an initial assignment, so it follows the change;
            every other quantity follows a change before the simulation only.
        name: the name of the model.
        notes: the notes of the model.
        observers: the observers of the definition, assignment rules of the
            top model after the observers of the library (a sum of
            metabolites, a ratio); each reads a state.
        creators: the creators of the model, its model history (sbmlutils
            `Creator`).
        parameters: the constant parameters of the top model which the model
            package states (`ModelParameter`), read by the formulas, the
            observers and the custom reactions of the definition.
        starts: the initial concentrations (mM) of species by their flat id
            (`MU__Cpl_crtn`, a pool `Cve_crtn`, a species of a custom
            reaction): a number, a `Value` or a `Formula` over the constants
            of the top model, the initial assignment of the species. For the
            organs which a production or an elimination of the chemistry
            moves away from the equilibrium with the baseline of an
            endogenous substance (`Substance.baseline`): their steady state,
            so that the body starts stationary.
    """

    sid: str
    substances: Sequence[Substance]
    chemistry: Sequence[Metabolism | Transport | CustomReaction] = ()
    variants: Variants = field(default_factory=Variants)
    routes: Sequence[Route] = ()
    pd: Sequence[PDComponent] = ()
    sex: Sex = Sex.MALE
    population: str = DEFAULT_POPULATION
    arm_vein: str = "pksim_default"
    dynamic: tuple[str, ...] = ()
    name: str | None = None
    notes: str | None = None
    observers: Sequence[Observer] = ()
    creators: Sequence[Creator] = ()
    parameters: Sequence[ModelParameter] = ()
    starts: Mapping[str, Q] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalise the sex to a `Sex` member."""
        object.__setattr__(self, "sex", as_member(Sex, self.sex, "PBPKDefinition.sex"))

    @classmethod
    def combine(cls, *definitions: "PBPKDefinition", sid: str) -> "PBPKDefinition":
        """Couple the definitions of several packages into one physiology.

        The substances, the chemistry, the routes, the PD components, the
        dynamic parameters, the observers and parameters of the definitions
        and the creators are the unions of the definitions (the observers of the library follow
        from the substances), the organ variants are merged; there is
        one physiology, so the definitions must agree on sex, population and
        the arm vein preset. Substance ids are unique over all packages, so
        nothing is renamed: an id which two definitions define is accepted
        only if both define equal objects (a shared metabolite; PD components
        are equal by the id of their model, their inputs and outputs; a
        creator is the same person when all its fields are equal). An
        organ variant is a requirement: an organ which a definition leaves
        out (the default organ) means no requirement, so the coupled model
        takes the variant another definition chooses (e.g. its zonated
        liver); two definitions which choose different variants of an organ
        raise.

        Args:
            *definitions: the definitions, e.g. of a drug and an inhibitor.
            sid: the id of the coupled model.

        Returns:
            The coupled definition, in the order of the definitions.

        Raises:
            ValueError: if no definition is given, the definitions differ in
                sex, population or arm vein preset, two definitions define an
                id (substance, process, route, PD component, observer,
                parameter) differently, or choose different variants of an organ; naming
                the id.
        """
        if not definitions:
            raise ValueError(f"Combine '{sid}': no definitions given.")
        for attribute in ("sex", "population", "arm_vein"):
            values = {getattr(d, attribute) for d in definitions}
            if len(values) > 1:
                raise ValueError(
                    f"Combine '{sid}': the definitions differ in {attribute} "
                    f"{sorted(str(v) for v in values)}, a coupled model has one "
                    "physiology."
                )
        first = definitions[0]
        sids = ", ".join(d.sid for d in definitions)
        return cls(
            sid=sid,
            substances=_union(
                sid,
                "substance",
                (x for d in definitions for x in d.substances),
                lambda x: x.sid,
            ),
            chemistry=_union(
                sid,
                "process",
                (x for d in definitions for x in d.chemistry),
                lambda x: x.sid,
            ),
            variants=Variants(
                organs=_merge_variants(sid, [d.variants for d in definitions])
            ),
            routes=_union(
                sid,
                "route",
                (x for d in definitions for x in d.routes),
                lambda x: (x.kind, x.substance),
            ),
            pd=_union(
                sid,
                "PD component",
                (x for d in definitions for x in d.pd),
                lambda x: str(x.model.sid),
                "has different inputs or outputs in two definitions",
            ),
            sex=first.sex,
            population=first.population,
            arm_vein=first.arm_vein,
            dynamic=tuple(dict.fromkeys(p for d in definitions for p in d.dynamic)),
            notes=f"The coupled physiology of the definitions {sids}.",
            observers=_union(
                sid,
                "observer",
                (x for d in definitions for x in d.observers),
                lambda x: x.sid,
            ),
            creators=_creators(definitions),
            parameters=_union(
                sid,
                "parameter",
                (x for d in definitions for x in d.parameters),
                lambda x: x.sid,
            ),
            starts=dict(
                _union(
                    sid,
                    "start of the species",
                    (x for d in definitions for x in d.starts.items()),
                    lambda x: x[0],
                )
            ),
        )

    def validate(self, *, symbols: bool = True) -> None:
        """Check the definition.

        The observers of the definition are checked against the top model:
        `validate` composes the model without them (`pkpdsim.pbpk.builder`)
        and checks that an observer id is no id of the body nor of another
        observer, that every symbol of its formula is an id of the top model
        (or another observer) and that it reads a state; `compose` makes the
        same check on the model it composes and calls `validate` without it.

        Args:
            symbols: check the observers against the top model.

        Raises:
            ValueError: with all problems found, naming the offending ids.
        """
        problems: list[str] = []
        problems += self._check_substances()
        known = {s.sid for s in self.substances}
        problems += self._check_chemistry(known)
        problems += self._check_routes(known)
        problems += self._check_dynamic()
        problems += self._check_observers()
        problems += self._check_parameters()
        problems += self._check_starts()
        if self.population not in POPULATIONS:
            problems.append(
                f"population '{self.population}' is unknown, the populations are "
                f"{sorted(POPULATIONS)}"
            )
        presets = sorted(load("arm_vein"))
        if self.arm_vein not in presets:
            problems.append(
                f"arm_vein '{self.arm_vein}' is no preset of data/arm_vein.toml, "
                f"the presets are {presets}"
            )
        if not problems and symbols and self.observers:
            # the builder imports the definition: imported on use
            from pkpdsim.pbpk.builder import observer_problems

            problems += observer_problems(self)
        if problems:
            raise ValueError(
                f"Invalid PBPK definition '{self.sid}':\n- " + "\n- ".join(problems)
            )

    def _check_observers(self) -> list[str]:
        """Check the ids and the formulas of the observers of the definition.

        Their symbols are checked against the top model by
        `pkpdsim.pbpk.builder.observer_problems`.
        """
        problems: list[str] = []
        sids = [o.sid for o in self.observers]
        twice = sorted({s for s in sids if sids.count(s) > 1})
        if twice:
            problems.append(f"duplicate observer ids: {twice}")
        for o in self.observers:
            if not _SID.fullmatch(o.sid):
                problems.append(f"observer id '{o.sid}' is not an SBML SId")
            elif "__" in o.sid:
                problems.append(
                    f"observer id '{o.sid}' contains '__', the separator of the "
                    "ids of the organs in the flat model"
                )
            try:
                functions = formula_functions(o.formula)
            except ValueError as err:
                problems.append(f"observer '{o.sid}': {err}")
                continue
            if functions:
                problems.append(
                    f"observer '{o.sid}' calls {sorted(functions)}, which are no "
                    "functions of SBML L3 infix (exp, ln, piecewise, ...); the "
                    "model defines no function definitions"
                )
        return problems

    def _check_parameters(self) -> list[str]:
        """Check the ids and the units of the parameters of the definition.

        That an id is no other id of the top model nor a local id of an organ
        is checked by `compose` (`pkpdsim.pbpk.builder`).
        """
        problems: list[str] = []
        sids = [p.sid for p in self.parameters]
        twice = sorted({s for s in sids if sids.count(s) > 1})
        if twice:
            problems.append(f"duplicate parameter ids: {twice}")
        taken = _physiology_ids()
        submodels = {ids.submodel(organ) for organ in Organ}
        for p in self.parameters:
            if not _SID.fullmatch(p.sid):
                problems.append(f"parameter id '{p.sid}' is not an SBML SId")
            elif "__" in p.sid:
                problems.append(
                    f"parameter id '{p.sid}' contains '__', the separator of the "
                    "ids of the organs in the flat model"
                )
            elif p.sid in taken:
                problems.append(
                    f"parameter '{p.sid}' is an id of the physiology of the library"
                )
            elif p.sid in submodels:
                problems.append(
                    f"parameter '{p.sid}' is the id of an organ submodel of the library"
                )
            problems += _unit(f"parameter '{p.sid}'", "value", p.value, p.unit)
        return problems

    def _check_starts(self) -> list[str]:
        """Check the ids and the values of the starts of the species.

        That an id is a species of the model is checked by `compose`
        (`pkpdsim.pbpk.builder`).
        """
        problems: list[str] = []
        for sid, q in self.starts.items():
            if not _SID.fullmatch(sid):
                problems.append(f"start '{sid}': the id is not an SBML SId")
                continue
            x = _number(q)
            if x is not None and not x >= 0:
                problems.append(f"start '{sid}': must not be negative, got {x}")
            problems += _unit(f"start '{sid}'", "value", q, U.mM)
        return problems

    def _check_substances(self) -> list[str]:
        """Check the ids and the mechanisms of the substances."""
        problems: list[str] = []
        sids = [s.sid for s in self.substances]
        duplicates = sorted({s for s in sids if sids.count(s) > 1})
        if duplicates:
            problems.append(f"duplicate substance ids: {duplicates}")
        taken = _physiology_ids()
        for s in self.substances:
            if not _SID.fullmatch(s.sid):
                problems.append(f"substance id '{s.sid}' is not an SBML SId")
            if "__" in s.sid:
                problems.append(f"substance id '{s.sid}' contains '__'")
            for token in _RESERVED:
                if s.sid.startswith(f"{token}_"):
                    problems.append(
                        f"substance id '{s.sid}' starts with the reserved token "
                        f"'{token}_' of the id scheme"
                    )
            if zone := _ZONE.match(s.sid):
                problems.append(
                    f"substance id '{s.sid}' starts with the reserved zone token "
                    f"'{zone.group()}' of the id scheme"
                )
            if s.sid in taken:
                problems.append(f"substance id '{s.sid}' is an id of the physiology")
            problems += self._check_substance(s)
        return problems

    def _check_substance(self, s: Substance) -> list[str]:
        """Check the mechanisms of a substance."""
        problems: list[str] = []
        for inh in s.inhibits:
            if inh.protein not in PROTEINS:
                problems.append(
                    f"{s.sid}: unknown protein '{inh.protein}' (inhibition)"
                )
        for ind in s.induces:
            if ind.protein not in PROTEINS:
                problems.append(f"{s.sid}: unknown protein '{ind.protein}' (induction)")
        problems += _check_interactions(s)
        problems += _positive(s.sid, "mr", s.mr)
        if s.binding is not None:
            problems += _in_unit_interval(s.sid, "fu", s.binding.fu, open_low=True)
        problems += _in_unit_interval(
            s.sid, "fraction_reabsorbed", s.excretion.fraction_reabsorbed
        )
        if s.baseline is not None:
            problems += _positive(s.sid, "baseline", s.baseline)
            problems += _unit(s.sid, "baseline", s.baseline, U.mM)
        ery = s.erythrocytes
        if ery is not None:
            problems += _check_erythrocytes(s.sid, ery)
        dist = s.distribution
        if dist is not None:
            if isinstance(dist, Kp | KpMethod):
                problems += _check_distribution(s.sid, dist)
            problems += check_kp_method(s)
            if (
                ery is not None
                and ery.mode == "kinetic"
                and dist.mode is Mode.PERFUSION
            ):
                problems.append(
                    f"{s.sid}: kinetic erythrocytes need a distribution with "
                    "exchange, not mode 'perfusion'"
                )
        return problems

    def _check_chemistry(self, known: set[str]) -> list[str]:
        """Check the processes of the chemistry."""
        problems: list[str] = []
        sids = [p.sid for p in self.chemistry]
        duplicates = sorted({s for s in sids if sids.count(s) > 1})
        if duplicates:
            problems.append(f"duplicate process ids: {duplicates}")
        for p in self.chemistry:
            if not _SID.fullmatch(p.sid):
                problems.append(f"process id '{p.sid}' is not an SBML SId")
            used: list[str] = []
            if isinstance(p, Metabolism):
                used = [p.substrate, *p.products]
                for product, fraction in p.products.items():
                    if fraction <= 0:
                        problems.append(
                            f"{p.sid}: product fraction of '{product}' must be "
                            f"positive, got {fraction}"
                        )
                total = sum(p.products.values())
                if p.products and abs(total - 1.0) > 1e-9:
                    problems.append(f"{p.sid}: product fractions sum to {total}, not 1")
            elif isinstance(p, Transport):
                used = [p.substance]
            if isinstance(p, Metabolism | Transport):
                problems += _check_kinetics(p.sid, p.kinetics)
            unknown = sorted(u for u in used if u not in known)
            if unknown:
                problems.append(f"{p.sid}: unknown substances {unknown}")
            protein = getattr(p, "protein", None)
            if protein is not None and protein not in PROTEINS:
                problems.append(f"{p.sid}: unknown protein '{protein}'")
        return problems

    def _check_dynamic(self) -> list[str]:
        """Check the names of the dynamic parameters.

        That they are parameters of the top model is checked by `compose`.
        """
        problems = [
            f"dynamic '{pid}' is not an SBML SId"
            for pid in self.dynamic
            if not _SID.fullmatch(pid)
        ]
        twice = sorted({p for p in self.dynamic if self.dynamic.count(p) > 1})
        if twice:
            problems.append(f"dynamic names {twice} more than once")
        return problems

    def _check_routes(self, known: set[str]) -> list[str]:
        """Check the routes."""
        return [
            f"route {r.kind}: unknown substance '{r.substance}'"
            for r in self.routes
            if r.substance not in known
        ]


def _union[T](
    sid: str,
    what: str,
    items: Iterable[T],
    key: Callable[[T], Hashable],
    differs: str = "is defined differently by two definitions",
) -> list[T]:
    """The union of the items of several definitions by their id, in order.

    Raises:
        ValueError: if two items of one id differ, naming the id and how
            they differ (`differs`).
    """
    result: dict[Hashable, T] = {}
    for item in items:
        k = key(item)
        if k in result and result[k] != item:
            name = " ".join(k) if isinstance(k, tuple) else k
            raise ValueError(f"Combine '{sid}': the {what} '{name}' {differs}.")
        result.setdefault(k, item)
    return list(result.values())


def _creators(definitions: Iterable["PBPKDefinition"]) -> list[Creator]:
    """The creators of several definitions, each person once, in order.

    A sbmlutils `Creator` compares and hashes by value: two creators are the
    same person when all their fields are equal.
    """
    return list(dict.fromkeys(c for d in definitions for c in d.creators))


def _merge_variants(
    sid: str, variants: Sequence[Variants]
) -> dict[Organ, "str | OrganComponent"]:
    """The organ variants of several definitions.

    Raises:
        ValueError: if two definitions choose different variants of an organ.
    """
    merged: dict[Organ, str | OrganComponent] = {}
    for v in variants:
        for organ, choice in v.organs.items():
            if organ in merged and merged[organ] != choice:
                raise ValueError(
                    f"Combine '{sid}': the definitions choose different variants "
                    f"of organ '{organ.value}'."
                )
            merged.setdefault(organ, choice)
    return merged


def _check_erythrocytes(sid: str, ery: Erythrocytes) -> list[str]:
    """Check the parameters of the erythrocytes of a substance."""
    problems = [
        *_positive(sid, "bp", ery.bp),
        *_positive(sid, "kbc", ery.kbc),
        *_positive(sid, "k_exchange", ery.k_exchange),
    ]
    if problems:
        return problems
    if ery.mode == "equilibrium":
        if (ery.bp is None) == (ery.kbc is None):
            return [f"{sid}: equilibrium erythrocytes need exactly one of bp and kbc"]
        if ery.k_exchange is not None:
            return [f"{sid}: equilibrium erythrocytes have no k_exchange"]
    elif ery.kbc is None or ery.k_exchange is None:
        return [f"{sid}: kinetic erythrocytes need kbc and k_exchange"]
    elif ery.bp is not None:
        return [f"{sid}: kinetic erythrocytes have no bp"]
    return []


def _check_distribution(sid: str, dist: Kp | KpMethod) -> list[str]:
    """Check the distribution of a substance."""
    problems: list[str] = []
    entered = dist.organs
    if dist.mode in {Mode.TISSUE, Mode.PERMEABILITY} and dist.ftissue is None:
        # an organ with an exchange of its own (`ps`) has no share of ftissue
        shared = sorted(o.value for o in (entered or Organ) if o not in dist.ps)
        if shared:
            problems.append(
                f"{sid}: mode '{dist.mode.value}' requires ftissue for the "
                f"organs without ps: {shared}"
            )
    named: dict[str, Mapping[Organ, Q]] = {"ps": dist.ps}
    if isinstance(dist, Kp):
        named["override"] = dist.overrides
        problems += _positive(sid, "Kp", dist.value)
        for organ, value in dist.overrides.items():
            problems += _positive(sid, f"Kp override '{organ.value}'", value)
    else:
        problems += _positive(sid, "factor", dist.factor)
    if dist.ftissue is not None:
        problems += _positive(sid, "ftissue", dist.ftissue)
    ps_unit = U.per_min if dist.ps_per_volume else U.l_per_min
    for organ, value in dist.ps.items():
        problems += _positive(sid, f"ps '{organ.value}'", value)
        problems += _unit(sid, f"ps '{organ.value}'", value, ps_unit)
    if entered is not None:
        for what, mapping in named.items():
            outside = sorted(o.value for o in mapping if o not in entered)
            if outside:
                problems.append(f"{sid}: {what} for organs not entered: {outside}")
    return problems


def _number(q: Q | None) -> float | None:
    """Get the number of a plain float or a `Value`, `None` for none or a `Formula`."""
    if isinstance(q, Formula):
        return None
    return q.value if isinstance(q, Value) else q


def _positive(where: str, name: str, q: Q | None) -> list[str]:
    """Check that a number is positive."""
    x = _number(q)
    if x is not None and not x > 0:
        return [f"{where}: {name} must be positive, got {x}"]
    return []


def _unit(where: str, name: str, q: Q | None, unit: UnitType) -> list[str]:
    """Check that a `Value` is given in the unit of its parameter.

    The numbers of a definition are in the units of the model; a `Value` in
    another unit (an exchange in l/min where one per litre of organ is meant)
    would be taken as it is.
    """
    if not isinstance(q, Value | Formula) or q.unit is None:
        return []
    given, expected = _unit_text(q.unit), _unit_text(unit)
    if given != expected:
        return [f"{where}: {name} must be in {expected}, got {given}"]
    return []


def _unit_text(unit: UnitType) -> str:
    """The definition of a unit, e.g. `1/min`."""
    return str(getattr(unit, "definition", unit))


def _in_unit_interval(
    where: str, name: str, q: Q | None, *, open_low: bool = False
) -> list[str]:
    """Check that a number is a fraction, in (0, 1] or [0, 1]."""
    x = _number(q)
    if x is None:
        return []
    low_ok = x > 0 if open_low else x >= 0
    if not (low_ok and x <= 1):
        interval = "(0, 1]" if open_low else "[0, 1]"
        return [f"{where}: {name} must be in {interval}, got {x}"]
    return []


def _check_interactions(s: Substance) -> list[str]:
    """Check the inhibitions and inductions of a substance.

    One inhibition and one induction per protein (a mechanism-based
    inhibition includes the competitive part), the constants of the kind of
    the inhibition and none other, positive constants.
    """
    problems: list[str] = []
    for what, proteins in (
        ("inhibitions", [i.protein for i in s.inhibits]),
        ("inductions", [i.protein for i in s.induces]),
    ):
        twice = sorted({p for p in proteins if proteins.count(p) > 1})
        if twice:
            problems.append(f"{s.sid}: several {what} of {twice}")
    for inh in s.inhibits:
        where = f"{s.sid} ({inh.kind} inhibition of {inh.protein})"
        if (inh.kiu is None) is (inh.kind == "mixed"):
            problems.append(
                f"{where}: kiu is required for a mixed inhibition and only there"
            )
        if (inh.kinact is None) is (inh.kind == "mechanism"):
            problems.append(
                f"{where}: kinact is required for a mechanism-based inhibition "
                "and only there"
            )
        for name in ("ki", "kiu", "kinact"):
            problems += _positive(where, name, getattr(inh, name))
    for ind in s.induces:
        where = f"{s.sid} (induction of {ind.protein})"
        problems += _positive(where, "emax", ind.emax)
        problems += _positive(where, "ec50", ind.ec50)
    return problems


def _check_kinetics(where: str, k: FirstOrder | MM | Hill) -> list[str]:
    """Check that the parameters of a rate law are positive."""
    problems: list[str] = []
    for name in ("k", "vmax", "km", "n"):
        problems += _positive(where, name, getattr(k, name, None))
    return problems
