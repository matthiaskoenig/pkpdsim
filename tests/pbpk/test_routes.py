"""Tests of the routes: the depot routes and the options of a tablet.

A depot route applies its dose state (`IMDOSE_<s>`, `SCDOSE_<s>`,
`RECDOSE_<s>`, mg) into a depot amount of the top model (`Aim_<s>`,
`Asc_<s>`, `Arec_<s>`), which is absorbed with first order into the venous
blood; the fraction `f_<k>_<s>` is absorbed, the rest leaves into a sink, and
the rectal route sends the fraction `f_rec_portal_<s>` of the absorbed dose
through the portal vein and the liver. A tablet without the options of its
dissolution (`fed`, `enteric`) builds the SBML it built before.
"""

import dataclasses
import math
from pathlib import Path
from typing import get_args, get_type_hints

import libsbml
import numpy as np
import pytest

from pkpdsim.pbpk import (
    Depot,
    DoseDependentDissolution,
    Excretion,
    FirstOrder,
    FirstOrderDissolution,
    Kp,
    Metabolism,
    Mode,
    Organ,
    Part,
    PBPKDefinition,
    Pool,
    Route,
    Substance,
    ids,
)
from pkpdsim.pbpk.builder import BuildResult, build
from pkpdsim.pbpk.equilibration import DOSE_PREFIXES
from pkpdsim.pbpk.routes import DEPOT

from .conftest import (
    assert_reference_sbml,
    chemistry_reference_definition,
    dateless,
    definition,
    dose_amount,
    endogenous_reference_definition,
    load,
    reference_definition,
    simple_drug,
    total_amount,
)

DAY = 24 * 60.0
"""One day (min)."""

S = "drg"
"""The substance of the tests."""


def _closed_drug(**changes: object) -> Substance:
    """The simple drug without excretion."""
    return dataclasses.replace(simple_drug(), excretion=Excretion(), **changes)


def _depot(routes: list[Route], sid: str, **changes: object) -> PBPKDefinition:
    """A definition of the closed drug and its metabolite with the routes."""
    return definition(
        sid=sid,
        substances=[_closed_drug(), Substance("drg_m", name="metabolite", mr=316.0)],
        chemistry=[],
        routes=routes,
        **changes,
    )


def _ids(path: Path) -> set[str]:
    """The ids of the parameters, species, compartments and reactions of a model."""
    model = libsbml.readSBMLFromFile(str(path)).getModel()
    found: set[str] = set()
    for elements in (
        model.getListOfParameters(),
        model.getListOfSpecies(),
        model.getListOfCompartments(),
        model.getListOfReactions(),
    ):
        found |= {e.getId() for e in elements}
    return found


def test_every_route_kind_is_a_dose_of_the_equilibration() -> None:
    """Every kind of a route has its dose state in `DOSE_PREFIXES`.

    `DOSE_PREFIXES` names the kinds itself, so a kind added to `Route` would
    otherwise be invisible to the equilibration, which moves a dose given
    before the initialization to the start, and to the battery of the
    credibility assessment, which reads the dose states of a model from that
    list. The kinds come from the `Literal` of the field, so a route of any
    kind, not only a depot one, fails this until its dose state is listed.
    """
    prefixes: dict[str, str] = {
        "iv": ids.iv_dose(""),
        "oral": ids.po_dose(""),
        "oral_solution": ids.sol_dose(""),
        "im": ids.depot_dose("im", ""),
        "sc": ids.depot_dose("sc", ""),
        "rectal": ids.depot_dose("rectal", ""),
    }
    kinds = set(get_args(get_type_hints(Route)["kind"]))
    assert kinds == set(prefixes)
    assert kinds >= DEPOT
    assert set(prefixes.values()) <= set(DOSE_PREFIXES)
    # the infusion rate belongs to the intravenous route
    assert ids.infusion_rate("") in DOSE_PREFIXES


def test_the_ids_of_the_depot_routes() -> None:
    """The ids of the dose states, depots, rates and fractions of the routes."""
    assert [ids.depot_dose(k, S) for k in ("im", "sc", "rectal")] == [
        "IMDOSE_drg",
        "SCDOSE_drg",
        "RECDOSE_drg",
    ]
    assert [ids.depot(k, S) for k in ("im", "sc", "rectal")] == [
        "Aim_drg",
        "Asc_drg",
        "Arec_drg",
    ]
    assert [ids.depot_absorption_rate(k, S) for k in ("im", "sc", "rectal")] == [
        "ka_im_drg",
        "ka_sc_drg",
        "ka_rec_drg",
    ]
    assert [ids.depot_fraction(k, S) for k in ("im", "sc", "rectal")] == [
        "f_im_drg",
        "f_sc_drg",
        "f_rec_drg",
    ]
    assert ids.rectal_portal_fraction(S) == "f_rec_portal_drg"


def test_a_depot_route_creates_its_dose_depot_and_absorption(tmp_path: Path) -> None:
    """Every depot route: dose state, depot, rate, fraction; the rectal portal fraction."""
    routes = [Route.im(S), Route.sc(S), Route.rectal(S)]
    found = _ids(build(_depot(routes, "depots"), tmp_path).flat_path)
    for kind in ("im", "sc", "rectal"):
        assert {
            ids.depot_dose(kind, S),
            ids.depot(kind, S),
            ids.depot_absorption_rate(kind, S),
            ids.depot_fraction(kind, S),
            ids.depot_application_time(kind, S),
            ids.depot_compartment(kind),
            ids.depot_application(kind, S),
            ids.depot_absorption(kind, S),
            ids.depot_loss_reaction(kind, S),
        } <= found, kind
    assert ids.rectal_portal_fraction(S) in found
    assert ids.depot_absorption("rectal", S, portal=True) in found
    assert {ids.depot_loss("im", S), ids.depot_loss("sc", S)} <= found
    # the substance without a route has no depot
    assert not {ids.depot(k, "drg_m") for k in ("im", "sc", "rectal")} & found


def test_the_defaults_of_the_depot_routes_state_their_source(tmp_path: Path) -> None:
    """The rates are choices of the library, the portal fraction is de Boer 1979."""
    routes = [Route.im(S), Route.rectal(S)]
    path = build(_depot(routes, "defaults"), tmp_path).flat_path
    model = libsbml.readSBMLFromFile(str(path)).getModel()

    def notes(pid: str) -> str:
        return model.getParameter(pid).getNotesString()

    assert "not verified" in notes(ids.depot_absorption_rate("im", S))
    assert "de Boer" in notes(ids.rectal_portal_fraction(S))
    r = load(build(_depot(routes, "defaults_values"), tmp_path))
    assert r[ids.depot_fraction("im", S)] == 1.0
    assert r[ids.rectal_portal_fraction(S)] == pytest.approx(0.45)


def test_a_model_sets_the_parameters_of_its_depot(tmp_path: Path) -> None:
    """The rate, the fraction, the portal fraction and the application of a route."""
    route = Route.rectal(S, ka=0.5, fraction=0.6, portal=0.3, application=600.0)
    r = load(build(_depot([route], "set"), tmp_path))
    assert r[ids.depot_absorption_rate("rectal", S)] == 0.5
    assert r[ids.depot_fraction("rectal", S)] == 0.6
    assert r[ids.rectal_portal_fraction(S)] == 0.3
    assert r[ids.depot_application_time("rectal", S)] == 600.0


def test_a_depot_of_the_wrong_route_is_refused() -> None:
    """A depot belongs to a depot route, a portal fraction to the rectal route."""
    with pytest.raises(ValueError, match="depot"):
        Route("iv", S, depot=Depot(ka=1.0))
    with pytest.raises(ValueError, match="portal"):
        Route("sc", S, depot=Depot(portal=0.5))


def test_an_intramuscular_dose_is_absorbed_into_the_venous_blood(
    tmp_path: Path,
) -> None:
    """The fraction `f_im` of the dose is absorbed, the rest leaves into the sink."""
    r = load(build(_depot([Route.im(S, fraction=0.8)], "im"), tmp_path))
    mmole = 10.0 / r[ids.mr(S)]
    r[ids.depot_dose("im", S)] = 10.0
    assert dose_amount(r, S) == pytest.approx(mmole, rel=1e-12)
    r.simulate(0, 60, 10)
    assert r[f"[{ids.pool_species(Pool.VENOUS, S)}]"] > 0
    r.simulate(60, 48 * 60, 50)
    depot, loss = r[ids.depot("im", S)], r[ids.depot_loss("im", S)]
    assert depot < 1e-12 * mmole
    assert dose_amount(r, S) < 1e-12 * mmole
    total = total_amount(r, S)
    assert total + dose_amount(r, S) == pytest.approx(mmole, rel=1e-9)
    absorbed = total - depot - loss
    assert absorbed == pytest.approx(0.8 * mmole, rel=1e-9)
    assert loss == pytest.approx(0.2 * mmole, rel=1e-9)


def test_a_subcutaneous_dose_is_absorbed_with_its_rate(tmp_path: Path) -> None:
    """After the application the depot falls with `exp(-ka_sc*t)`."""
    r = load(build(_depot([Route.sc(S, ka=2.0)], "sc"), tmp_path))
    r[ids.depot_dose("sc", S)] = 10.0
    r.timeCourseSelections = ["time", ids.depot("sc", S)]
    result = r.simulate(0, 240, 241)
    t, a = result[:, 0], result[:, 1]
    ka = r[ids.depot_absorption_rate("sc", S)] / 60.0
    # from 10 min on, the application (half-life 10 s) is over
    later = t >= 10
    expected = a[later][0] * np.exp(-ka * (t[later] - t[later][0]))
    assert np.allclose(a[later], expected, rtol=1e-6, atol=0)
    assert ka == pytest.approx(2.0 / 60.0)


def _rectal_auc(result: BuildResult, portal: float) -> tuple[float, float]:
    """The systemic AUC of a rectal dose and the hepatic availability during it."""
    r = load(result)
    r[ids.rectal_portal_fraction(S)] = portal
    r[ids.depot_dose("rectal", S)] = 10.0
    venous = f"[{ids.pool_species(Pool.VENOUS, S)}]"
    out = ids.flat(Organ.LIVER, ids.outflow_reaction(Part.PLASMA, S))
    metabolism = ids.flat(Organ.LIVER, "DRGLI")
    r.timeCourseSelections = ["time", venous, out, metabolism]
    result = r.simulate(0, 3 * DAY, 30001)
    t, c = result[:, 0], result[:, 1]
    # the liver is one well-stirred compartment (perfusion-limited, first
    # order): the fraction of its inflow which leaves it unmetabolized is the
    # same at every time, the hepatic availability
    j_out, j_met = result[-20000, 2], result[-20000, 3]
    return float(np.trapezoid(c, t)), float(j_out / (j_out + j_met))


def test_a_rectal_dose_passes_the_liver_by_its_portal_fraction(
    tmp_path: Path,
) -> None:
    """Through the portal vein the systemic AUC is the AUC times the hepatic availability."""
    drug = _closed_drug(distribution=Kp(2.0, mode=Mode.PERFUSION))
    d = definition(
        sid="rectal",
        substances=[drug, Substance("drg_m", name="metabolite", mr=316.0)],
        chemistry=[
            Metabolism(
                sid="DRGLI",
                organ=Organ.LIVER,
                substrate=S,
                products={"drg_m": 1.0},
                kinetics=FirstOrder(k=10.0),
            )
        ],
        routes=[Route.rectal(S)],
    )
    result = build(d, tmp_path)
    systemic, fh = _rectal_auc(result, 0.0)
    portal, fh_portal = _rectal_auc(result, 1.0)
    assert fh == pytest.approx(fh_portal, rel=1e-6)
    assert 0.2 < fh < 0.8
    assert portal == pytest.approx(systemic * fh, rel=1e-3)


def test_a_second_depot_dose_adds_to_the_depot(tmp_path: Path) -> None:
    """A change of the dose state adds its amount, as `add_doses` expects of a dose."""
    r = load(build(_depot([Route.im(S, ka=0.5)], "im_twice"), tmp_path))
    mr = r[ids.mr(S)]
    dose = ids.depot_dose("im", S)
    r[dose] = 10.0
    r.simulate(0, 60, 10)
    before = r[ids.depot("im", S)]
    assert 0 < before < 10.0 / mr
    # `PKPDSimulationExperiment.add_doses`: `IMDOSE_drg + 5.0`
    r[dose] = r[dose] + 5.0
    r.simulate(60, 61, 10)
    # within the minute the depot takes the new dose (half-life 10 s) and
    # absorbs less than 1 % of it
    added = r[ids.depot("im", S)]
    decay = math.exp(-0.5 / 60.0)
    assert added == pytest.approx(before * decay + 5.0 / mr, rel=1e-2)
    r.simulate(61, 3 * DAY, 50)
    assert total_amount(r, S) + dose_amount(r, S) == pytest.approx(15.0 / mr, rel=1e-9)


def test_a_definition_without_depot_routes_is_unchanged(tmp_path: Path) -> None:
    """The stored references: the depot routes add nothing to a definition without them."""
    assert_reference_sbml(build(reference_definition(), tmp_path / "base"))
    assert_reference_sbml(
        build(chemistry_reference_definition(), tmp_path / "chemistry")
    )


def test_a_depot_route_adds_only_its_own_elements(tmp_path: Path) -> None:
    """The intramuscular route of a definition adds its elements and nothing else."""
    base = _ids(build(_depot([Route.iv(S)], "iv_only"), tmp_path).flat_path)
    im = _ids(build(_depot([Route.iv(S), Route.im(S)], "iv_im"), tmp_path).flat_path)
    assert base <= im
    assert im - base == {
        ids.depot_dose("im", S),
        ids.depot("im", S),
        ids.depot_loss("im", S),
        ids.depot_absorption_rate("im", S),
        ids.depot_fraction("im", S),
        ids.depot_application_time("im", S),
        ids.depot_compartment("im"),
        ids.depot_application("im", S),
        ids.depot_absorption("im", S),
        ids.depot_loss_reaction("im", S),
    }


def _flat_text(d: PBPKDefinition, path: Path) -> list[str]:
    """The lines of the flat SBML of a definition, without the dates."""
    return dateless(build(d, path).flat_path.read_text(encoding="utf-8"))


def test_a_dissolution_without_the_options_is_unchanged(tmp_path: Path) -> None:
    """Without `fed` and `enteric` a tablet builds the SBML it built before.

    The stored references with the default tablet are unchanged; a
    dose-dependent dissolution has no factor of the fed state, its rate is
    `Ka_dis/60*D50/(D50 + D)` and the tablet dissolves in the stomach; the
    defaults stated explicitly build the same SBML.
    """
    assert_reference_sbml(build(reference_definition(), tmp_path / "base"))
    assert_reference_sbml(
        build(endogenous_reference_definition(), tmp_path / "endogenous")
    )
    for name, plain, explicit in (
        (
            "first_order",
            FirstOrderDissolution(k=2.0),
            FirstOrderDissolution(k=2.0, fed=None, release=None, enteric=False),
        ),
        (
            "dose_dependent",
            DoseDependentDissolution(k=2.0, dose50=50.0),
            DoseDependentDissolution(k=2.0, dose50=50.0, fed=None, enteric=False),
        ),
    ):
        d = _depot([Route.oral(S, plain)], f"tablet_{name}")
        same = dataclasses.replace(d, routes=[Route.oral(S, explicit)])
        assert _flat_text(d, tmp_path / name) == _flat_text(
            same, tmp_path / f"{name}_explicit"
        )
    result = build(
        _depot([Route.oral(S, DoseDependentDissolution(k=2.0, dose50=50.0))], "dd"),
        tmp_path / "dd",
    )
    found = _ids(result.flat_path)
    assert not {i for i in found if i.startswith("f_dis_")}
    gu = Organ.SMALL_INTESTINE
    assert ids.flat(gu, ids.dose_reaction("po", S)) in found
    model = libsbml.readSBMLFromFile(str(result.flat_path)).getModel()
    rule = model.getAssignmentRuleByVariable(ids.flat(gu, ids.dissolution_rate(S)))
    solids = " + ".join(ids.flat(gu, ids.solid(S, seg)) for seg in ids.GUT_SEGMENTS[1:])
    k, d50 = ids.dissolution(S), ids.dissolution_dose50(S)
    assert libsbml.formulaToL3String(rule.getMath()) == (
        f"({k} / 60 GU__min_per_hr) * {d50} / ({d50} + {ids.po_dose(S)} + "
        f"{ids.mr(S)} * ({solids}))"
    )
