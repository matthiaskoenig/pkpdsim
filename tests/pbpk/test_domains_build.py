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

from pkpdsim.pbpk import (
    Adipose,
    FirstOrder,
    LumenConversion,
    Metabolism,
    Organ,
    PBPKDefinition,
    Pregnancy,
    RegionalAbsorption,
    Route,
    SegmentedIntestine,
    Sex,
    Species,
    Variants,
    ZonatedLiver,
    ids,
)
from pkpdsim.pbpk.builder import build
from pkpdsim.pbpk.interventions import scenarios
from pkpdsim.pbpk.pregnancy import TERM
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Formula
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

MIN_IDS = 400
"""The floor of the ids a shape checks: a sweep which checks nothing looks like a
green run. The shapes check 459 (the impaired default) to 995 (endogenous)."""

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
            variants=Variants(organs={Organ.LIVER: LobedLiver()}),
        )
        changes = scenarios(Species.RAT)["hepatectomy90"].changes(RAT_TOKENS)
        return d, [(changes | _iv("drg"), DAY), (_iv("drg"), DAY)]

    def impaired() -> Shape:
        d = definition(sid="impaired")
        worst = {ids.CP_SEVERITY: 3.0, ids.HF_SEVERITY: 4.0, ids.F_RENAL_FUNCTION: 0.05}
        return d, [(worst | _both("drg"), DAY)]

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
