"""Tests that every value of a built model is inside its domain.

The domain of a quantity is written into the SBML by its kind (`KINDS` of
`pkpdsim.pbpk.annotation`), classified by hand. This is the check of that
classification: it builds models of the shapes in which a kind exists, reads
the value of every quantity which has a domain at the start and along a time
course, and requires it inside the domain. A domain which is too wide revives
the bug the contract exists to fix (a fraction absorbed bounded by 1.127), one
which is too narrow clips a legitimate value; no other test looks at either.

A failure names a kind: either the domain of the row is wrong, or the value is
wrong, which is a bug of the library.

The limit of the check: the sweep shows that no domain is too narrow for the
values the models actually reach, and it cannot show that an upper bound is too
wide when no value ever approaches it. A green file is no proof of more.
The check reads sample points (every 30 min of the course), not an envelope, so
a brief excursion between two points is not seen; a monotone quantity like the
fraction absorbed of pkpdsim#207 is caught either way.
"""

import math
from collections.abc import Callable
from pathlib import Path

import pytest
from sbmlutils.factory import Parameter
from sbmlutils.factory import Species as SpeciesSBML

from pkpdsim.pbpk import (
    REST,
    Adipose,
    Binding,
    Circadian,
    CustomReaction,
    Direction,
    Excretion,
    FirstOrder,
    Formula,
    Hill,
    Kp,
    LumenConversion,
    Membrane,
    Metabolism,
    Mode,
    ModelParameter,
    Observer,
    Organ,
    Part,
    PBPKDefinition,
    Pregnancy,
    PrescribedVolume,
    RateOf,
    RegionalAbsorption,
    Route,
    Saturable,
    SegmentedIntestine,
    Sex,
    Species,
    Substance,
    ThresholdReabsorption,
    Transport,
    U,
    Value,
    Variants,
    ZonatedLiver,
    ids,
)
from pkpdsim.pbpk.builder import build
from pkpdsim.pbpk.interventions import scenarios
from pkpdsim.pbpk.pregnancy import TERM
from pkpdsim.pbpk.remodelling import FlowDriven
from pkpdsim.pbpk.variants.liver_lobes import LobedLiver
from pkpdsim.sbml.domains import (
    BOOLEAN,
    FRACTION,
    NONNEGATIVE,
    REAL,
    Domain,
    Values,
    domains,
    integers,
)
from pkpdsim.zonation.settings import Geometry

from .conftest import (
    chemistry_reference_definition,
    definition,
    endogenous_reference_definition,
    individual,
    nephron_reference_definition,
    reference_definition,
    turnover_reference_definition,
)

ABS_TOL = 1e-12
"""The tolerance of a bound (the units of the quantity), the absolute tolerance
of the integrator (`_check`).

It absorbs the noise of the solver on a quantity which is zero in exact
arithmetic and nothing else: a state which decays to 0 (an amount in a lumen
after the dose has passed, `IVDOSE_drg`) dips to -3e-17 at this tolerance and to
-3e-12 at the default one of roadrunner. The integrator runs at abs 1e-12 and
rel 1e-10; with it the worst excursion beyond any bound across every shape,
measured with tolerance 0, is 2.6e-17, all of it states decaying to 0, so the
tolerance is noise and no slack for a real violation. It is absolute only: no
quantity of the sweep left its upper bound at all, so a relative slack would
only widen the domains (the bug of pkpdsim#207 is a bound off by 13 percent)."""

RAT_TOKENS = ("cl", "lll", "lml", "rl", "rml")
"""The lobes of the liver of the rat."""

DOSE = 10.0
"""The dose of every route the shapes give (the units of the model)."""

DAY = 24 * 60.0
"""A day (min)."""

SAMPLE = 30.0
"""The interval between the points of a time course (min)."""

MIN_IDS = 350
"""The floor of the ids a shape checks: a sweep which checks nothing looks like a
green run. The shapes check 385 (the steady state, no gut) to 995 (endogenous)."""

Individual = tuple[dict[str, float], float]
"""The changes of an individual before the initialization and the span (min)
of its time course."""

Shape = tuple[PBPKDefinition, list[Individual]]
"""A definition and the individuals which are built from it and simulated."""


def _outside(
    values: dict[str, float], bounds: dict[str, Domain], when: str
) -> list[tuple[str, str, float, str]]:
    """The quantities whose value is outside of its domain beyond the tolerance.

    Args:
        values: the values by id.
        bounds: the domains by id.
        when: the label of the point of the time course.

    Returns:
        The id, the label, the value and the domain of each, sorted.
    """
    outside = []
    for uid, d in bounds.items():
        x = values[uid]
        if d.values is not Values.CONTINUOUS and abs(x - round(x)) <= ABS_TOL:
            x = float(round(x))
        if not Domain(d.lower - ABS_TOL, d.upper + ABS_TOL, d.values).contains(x):
            outside.append((uid, when, values[uid], str(d)))
    return sorted(outside)


def _check(path: Path, changes: dict[str, float], span: float) -> int:
    """Check every quantity with a domain of an individual, at its start and over time.

    Args:
        path: the flat model.
        changes: the values of parameters by id before the initialization.
        span: the span of the time course (min).

    Returns:
        The number of quantities checked.

    Raises:
        AssertionError: if fewer than `MIN_IDS` quantities have a domain.
    """
    bounds = domains(path)
    r = individual(path, **changes)
    r.integrator.absolute_tolerance = 1e-12
    r.integrator.relative_tolerance = 1e-10
    names = sorted(bounds)
    assert len(names) >= MIN_IDS, f"only {len(names)} quantities have a domain"
    start = {uid: float(r.getValue(uid)) for uid in names}
    outside = _outside(start, bounds, "t=0")
    r.timeCourseSelections = ["time", *names]
    result = r.simulate(0.0, span, round(span / SAMPLE) + 1)
    for row in result:
        values = dict(zip(names, map(float, row[1:]), strict=True))
        outside += _outside(values, bounds, f"t={row[0]:g}")
    # one entry per quantity: the first time it left its domain, the value there
    first = {}
    for entry in outside:
        first.setdefault(entry[0], entry)
    assert sorted(first.values()) == [], f"individual {changes}"
    return len(names)


def _iv(sid: str) -> dict[str, float]:
    """An intravenous dose of a substance."""
    return {ids.iv_dose(sid): DOSE}


def _both(sid: str) -> dict[str, float]:
    """An intravenous and an oral dose of a substance."""
    return _iv(sid) | {ids.po_dose(sid): DOSE}


def _dose(route: Route) -> str:
    """The id of the dose of a route."""
    match route.kind:
        case "iv":
            return ids.iv_dose(route.substance)
        case "oral":
            return ids.po_dose(route.substance)
        case "oral_solution":
            return ids.sol_dose(route.substance)
        case kind:
            return ids.depot_dose(kind, route.substance)


def _reference(make: Callable[[], PBPKDefinition]) -> Shape:
    """A reference definition with a dose on every route it has."""
    d = make()
    changes: dict[str, float] = {}
    for route in d.routes:
        changes[_dose(route)] = DOSE
    return d, [(changes, DAY)]


def _shapes() -> dict[str, Callable[[], Shape]]:
    """The shapes: a definition, the individuals and the span of each."""

    def zonated() -> Shape:
        liver = ZonatedLiver(n_zones=3, profiles={"CYP3A4": [1, 2, 3]})
        d = definition(sid="zonated", variants=Variants(organs={Organ.LIVER: liver}))
        return d, [(_both("drg"), DAY)]

    def segmented() -> Shape:
        # the options which create parameters: the default intestine has none
        intestine = SegmentedIntestine(
            ka={"drg": 2.0},
            fabs=0.6,
            f_ka_colon=0.5,
            regional=RegionalAbsorption({"ile2": 0.5}),
            conversions=[LumenConversion("drg", "drg_m", 0.5, ["st", "cae"])],
        )
        d = definition(
            sid="segmented",
            variants=Variants(organs={Organ.SMALL_INTESTINE: intestine}),
        )
        return d, [({ids.po_dose("drg"): DOSE}, 2 * DAY)]

    def adipose() -> Shape:
        effect = Formula(
            "-0.2 dimensionless * Cve_drg / (Cve_drg + 0.0002 mM)", U.dimensionless
        )
        d = definition(sid="adipose", adipose=Adipose(effect=effect))
        # a day is too short for the weight to follow the effect (k_bw)
        return d, [({"BW": 105.0, "HEIGHT": 176.0} | _both("drg"), 30 * DAY)]

    def pregnancy() -> Shape:
        d = definition(sid="pregnant", sex=Sex.FEMALE, pregnancy=Pregnancy())
        return d, [({"GA": ga} | _both("drg"), DAY) for ga in (0.0, 13.0, 26.0, TERM)]

    def hepatectomy() -> Shape:
        d = definition(
            sid="lobes",
            species=Species.RAT,
            chemistry=[
                Metabolism(
                    sid="DRG2M",
                    organ=Organ.LIVER,
                    substrate="drg",
                    products={"drg_m": 1.0},
                    kinetics=FirstOrder(k=0.5),
                    protein="CYP3A4",
                )
            ],
            routes=[Route.iv("drg")],
            variants=Variants(
                organs={Organ.LIVER: LobedLiver(remodelling=FlowDriven())}
            ),
        )
        changes = scenarios(Species.RAT)["hepatectomy90"].changes(RAT_TOKENS)
        # the remodelling of the remnant takes days
        return d, [(changes | _iv("drg"), 7 * DAY), (_iv("drg"), DAY)]

    def impaired() -> Shape:
        d = definition(sid="impaired")
        worst = {ids.CP_SEVERITY: 3.0, ids.HF_SEVERITY: 4.0, ids.F_RENAL_FUNCTION: 0.05}
        return d, [(worst | _both("drg"), DAY)]

    def readouts() -> Shape:
        # the readouts of the intestine and an observer of a species of an organ,
        # on an intestine whose fabs is a fraction of the dose (pkpdsim#207)
        flat = ids.flat(Organ.MUSCLE, ids.species(Part.TISSUE, "drg"))
        d = definition(
            sid="readouts",
            variants=Variants(
                organs={Organ.SMALL_INTESTINE: SegmentedIntestine(fabs=0.94)}
            ),
            chemistry=[
                Transport(
                    "DRGEFF",
                    Organ.SMALL_INTESTINE,
                    "drg",
                    Direction.EFFLUX,
                    Membrane.APICAL,
                    FirstOrder(k=0.1, per_volume=False),
                )
            ],
            readouts=tuple(
                f"{r}_drg" for r in ("Alumen", "Rabs", "Reff", "Rfeces", "Aabs", "Fa")
            ),
            observers=[Observer("Cmu_twice", f"2 dimensionless * {flat}", U.mM, "x2")],
        )
        return d, [({ids.po_dose("drg"): DOSE}, 2 * DAY)]

    def depots() -> Shape:
        d = definition(
            sid="depots",
            routes=[
                Route.iv("drg"),
                Route.oral("drg"),
                Route.oral_solution("drg"),
                Route.im("drg"),
                Route.sc("drg"),
                Route.rectal("drg"),
                Route.inhaled("drg"),
            ],
        )
        return d, [({_dose(r): DOSE for r in d.routes}, DAY)]

    def circadian() -> Shape:
        d = definition(sid="circadian", circadian=Circadian())
        return d, [({ids.CLOCK: c} | _both("drg"), DAY) for c in (0.0, 8.0, 20.0)]

    def steady() -> Shape:
        # the generated steady state of an endogenous substance: a source, a
        # linear clearance and a saturable secretion
        chemistry = [
            CustomReaction(
                "PROD",
                Organ.MUSCLE,
                Part.TISSUE,
                "-> Cti_end",
                "k_prod",
                stationary=Formula("k_prod", U.mmole_per_min),
            ),
            CustomReaction(
                "LICL",
                Organ.LIVER,
                Part.TISSUE,
                "Cti_end ->",
                "k_cl * Vti * Cti_end",
                stationary="linear",
            ),
            CustomReaction(
                "SEC",
                Organ.KIDNEY,
                Part.PLASMA,
                "Cpl_end -> Aurine_end",
                "Vm * Cpl_end / (Km + Cpl_end)",
                stationary=Saturable("Vm", "Km"),
            ),
        ]
        d = PBPKDefinition(
            sid="steady",
            substances=[
                Substance(
                    "end",
                    name="endogenous substance",
                    mr=113.1,
                    binding=Binding(fu=0.5, protein="albumin"),
                    distribution=Kp(0.6, ftissue=6.0),
                    excretion=Excretion(urine=True),
                )
            ],
            chemistry=chemistry,
            parameters=[
                ModelParameter("k_prod", 0.0074, U.mmole_per_min, domain=NONNEGATIVE),
                ModelParameter("k_cl", 0.01, U.per_min, domain=NONNEGATIVE),
                ModelParameter("Vm", 0.01, U.mmole_per_min, domain=NONNEGATIVE),
                ModelParameter("Km", 0.05, U.mM, domain=NONNEGATIVE),
            ],
        )
        return d, [({}, DAY)]

    def processes() -> Shape:
        # the options of a process which create parameters: a clearance of the
        # whole organ and one per litre, a Hill law, a ratio of two capacities, a
        # product fraction by a value and by the remainder, and an interstitial
        # scale of the distribution
        d = definition(
            sid="processes",
            substances=[
                Substance(
                    "drg",
                    name="drug",
                    mr=300.0,
                    binding=Binding(fu=0.1, protein="albumin"),
                    distribution=Kp(
                        interstitial=1.0, overrides={Organ.LIVER: 2.0}, ftissue=1.0
                    ),
                    excretion=Excretion(urine=True),
                ),
                Substance(
                    "drg_m1",
                    name="metabolite 1",
                    mr=316.0,
                    excretion=Excretion(
                        urine=True, reabsorption=ThresholdReabsorption(0.1, 2.0)
                    ),
                ),
                Substance("drg_m2", name="metabolite 2", mr=332.0),
            ],
            chemistry=[
                Metabolism(
                    "DRG2M",
                    Organ.LIVER,
                    "drg",
                    products={
                        "drg_m1": Value(0.25, U.dimensionless),
                        "drg_m2": REST,
                    },
                    kinetics=FirstOrder(k=0.05, per_volume=False),
                    protein="CYP3A4",
                ),
                Metabolism(
                    "DRG2MB",
                    Organ.LIVER,
                    "drg",
                    products={"drg_m1": 1.0},
                    kinetics=FirstOrder(RateOf("DRG2M", 0.5), per_volume=False),
                    protein="CYP3A4",
                ),
                Transport(
                    "DRGIM",
                    Organ.LIVER,
                    "drg",
                    Direction.INFLUX,
                    Membrane.BASOLATERAL,
                    FirstOrder(k=0.1, per_volume=False),
                ),
                Transport(
                    "DRGIM2",
                    Organ.LIVER,
                    "drg",
                    Direction.INFLUX,
                    Membrane.BASOLATERAL,
                    FirstOrder(k=0.1),
                ),
                Transport(
                    "DRGSEC",
                    Organ.LIVER,
                    "drg",
                    Direction.EFFLUX,
                    Membrane.APICAL,
                    Hill(vmax=1e-3, km=0.01, n=2.0),
                ),
            ],
            routes=[Route.iv("drg"), Route.oral("drg")],
        )
        return d, [(_both("drg"), DAY)]

    def sinusoid() -> Shape:
        # the lobes of the rat with the remodelling of their volumes and the
        # zones of a sinusoid (the resistance network has algebraic rules which
        # roadrunner cannot simulate, so its kinds stay unreached)
        liver = LobedLiver(
            zonation=ZonatedLiver(n_zones=3, geometry=Geometry.SINUSOID),
            remodelling=FlowDriven(),
        )
        d = definition(
            sid="sinusoid",
            species=Species.RAT,
            substances=[
                Substance(
                    "drg",
                    name="drug",
                    mr=300.0,
                    binding=Binding(fu=0.1, protein="albumin"),
                    distribution=Kp(2.0, mode=Mode.PERMEABILITY, ftissue=1.0),
                    excretion=Excretion(urine=True),
                ),
                *definition().substances[1:],
            ],
            routes=[Route.iv("drg")],
            variants=Variants(organs={Organ.LIVER: liver}),
        )
        changes = scenarios(Species.RAT)["hepatectomy90"].changes(RAT_TOKENS)
        return d, [(changes | _iv("drg"), DAY), (_iv("drg"), DAY)]

    def aged() -> Shape:
        # the covariates over the age
        d = definition(sid="aged", adipose=Adipose(), age=True)
        return d, [({ids.AGE: age} | _both("drg"), DAY) for age in (20.0, 60.0)]

    def amounts() -> Shape:
        # an amount of a custom reaction and an observer which reads it
        reaction = CustomReaction(
            "DRGDEEP",
            Organ.MUSCLE,
            Part.TISSUE,
            "Cti_drg -> Adeep_drg",
            "k_deep * Vti * Cti_drg",
            parameters=[Parameter("k_deep", 1e-3, U.per_min)],
            species=[
                SpeciesSBML(
                    "Adeep_drg",
                    "Vti",
                    initialAmount=0.0,
                    substanceUnit=U.mmole,
                    hasOnlySubstanceUnits=True,
                )
            ],
        )
        deep = ids.flat(Organ.MUSCLE, "Adeep_drg")
        d = definition(
            sid="amounts",
            chemistry=[reaction],
            observers=[Observer("Adeep_total", deep, U.mmole, "deep amount")],
        )
        return d, [(_iv("drg"), DAY)]

    def prescribed() -> Shape:
        # the volume of a lobe which follows a prescribed course
        curve = Formula(
            f"{ids.VLI_REF} * 0.3 dimensionless * "
            "(1 dimensionless + 1 dimensionless * time / 7200 min)",
            "l",
        )
        d = definition(
            sid="prescribed",
            species=Species.RAT,
            routes=[Route.iv("drg")],
            variants=Variants(
                organs={
                    Organ.LIVER: LobedLiver(
                        remodelling=PrescribedVolume({"lll": curve})
                    )
                }
            ),
        )
        return d, [(_iv("drg"), DAY)]

    return {
        "reference": lambda: _reference(reference_definition),
        "chemistry": lambda: _reference(chemistry_reference_definition),
        "nephron": lambda: _reference(nephron_reference_definition),
        "endogenous": lambda: _reference(endogenous_reference_definition),
        "turnover": lambda: _reference(turnover_reference_definition),
        "zonated": zonated,
        "segmented": segmented,
        "adipose": adipose,
        "pregnancy": pregnancy,
        "hepatectomy": hepatectomy,
        "impaired": impaired,
        "readouts": readouts,
        "depots": depots,
        "circadian": circadian,
        "steady": steady,
        "processes": processes,
        "sinusoid": sinusoid,
        "prescribed": prescribed,
        "aged": aged,
        "amounts": amounts,
    }


SHAPES = _shapes()


@pytest.mark.parametrize("name", list(SHAPES))
def test_every_value_of_a_build_is_inside_its_domain(name: str, tmp_path: Path) -> None:
    """Every quantity with a domain stays in it, at the start and along the time course."""
    d, individuals = SHAPES[name]()
    path = build(d, tmp_path).flat_path
    for changes, span in individuals:
        _check(path, changes, span)


def test_the_check_names_a_value_outside_and_absorbs_noise() -> None:
    """A fraction of 1.127 (pkpdsim#207) is named, the noise of a zero is not."""
    bounds = {
        "f": FRACTION,
        "v": NONNEGATIVE,
        "n": integers(0, 3),
        "b": BOOLEAN,
        "x": REAL,
    }
    noise = {"f": 1.0 + 1e-15, "v": -1e-17, "n": 2.0, "b": 1.0, "x": -5.0}
    assert _outside(noise, bounds, "t=0") == []
    wrong = {"f": 1.127, "v": -1e-6, "n": 1.5, "b": 2.0, "x": math.nan}
    assert [uid for uid, *_ in _outside(wrong, bounds, "t=0")] == sorted(bounds)
