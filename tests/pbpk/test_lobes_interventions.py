"""Tests of the interventions on the lobes of the liver: ligation and resection.

An intervention is a set of changes of parameters of the model, as the
impairment maps of the library are: a portal vein ligation closes the portal
branch of its lobes (`g_po_<lobe>` 0), an arterial ligation their hepatic
artery branch, and a resection removes them (`f_present_<lobe>` 0 as well).
The model has no events, so a study applies the changes before a simulation or
between the phases of one, which is what a surgery at a time point is.

The flows of the lobes which stay open take over the flow of the branches
which closed, so the hepatic inflow of the liver is kept whatever is ligated;
that invariant and the two degenerate cases, every branch closed and every
lobe resected, are the tests of the structure.
"""

from pathlib import Path

import numpy as np
import pytest
import roadrunner

from pkpdsim.pbpk import (
    FirstOrder,
    Metabolism,
    Organ,
    PBPKDefinition,
    Pool,
    Route,
    Variants,
    ids,
)
from pkpdsim.pbpk.builder import build
from pkpdsim.pbpk.interventions import (
    ArterialLigation,
    Intervention,
    PortalLigation,
    Resection,
    scenario_changes,
    scenarios,
)
from pkpdsim.pbpk.organ import OrganComponent
from pkpdsim.pbpk.remodelling import FlowDriven
from pkpdsim.pbpk.variants.liver_lobes import LobedLiver
from pkpdsim.species import Species

from .conftest import definition, individual
from .test_lobes_build import RAT_TOKENS

S = "drg"
DAY = 24 * 60.0


def _definition(liver: OrganComponent | None) -> PBPKDefinition:
    return definition(
        sid="m",
        species=Species.RAT,
        chemistry=[
            Metabolism(
                sid="DRG2M",
                organ=Organ.LIVER,
                substrate=S,
                products={"drg_m": 1.0},
                kinetics=FirstOrder(k=0.5),
                protein="CYP3A4",
            )
        ],
        routes=[Route.iv(S)],
        variants=Variants(organs={Organ.LIVER: liver})
        if liver is not None
        else Variants(),
    )


def _built(tmp_path: Path, liver: OrganComponent | None) -> Path:
    return build(_definition(liver), tmp_path).flat_path


def _after(
    tmp_path: Path, *interventions: Intervention, **kwargs: object
) -> roadrunner.RoadRunner:
    """A model of the rat lobes with the changes of the interventions applied."""
    liver = LobedLiver(**kwargs)  # ty: ignore[invalid-argument-type]
    changes: dict[str, float] = {}
    for intervention in interventions:
        changes |= intervention.changes(RAT_TOKENS)
    return individual(_built(tmp_path, liver), **changes)


def test_a_portal_ligation_closes_only_its_lobes(tmp_path: Path) -> None:
    r = _after(tmp_path, PortalLigation(["lll", "lml"]))
    assert r[ids.portal_conductance("lll")] == 0.0
    assert r[ids.portal_conductance("lml")] == 0.0
    assert r[ids.portal_conductance("rl")] == 1.0
    # the artery of a deportalised lobe stays open
    assert r[ids.arterial_conductance("lll")] == 1.0


def test_a_portal_ligation_redistributes_the_portal_flow(tmp_path: Path) -> None:
    """The lobes which stay open take over the flow of the branch which closed."""
    r = _after(tmp_path, PortalLigation(["lll", "lml"]))
    assert r[ids.lobe_flow(ids.QPO, "lll")] == pytest.approx(0.0)
    assert r[ids.lobe_flow(ids.QPO, "lml")] == pytest.approx(0.0)
    total = sum(r[ids.lobe_flow(ids.QPO, t)] for t in RAT_TOKENS)
    assert total == pytest.approx(r[ids.QPO], rel=1e-9)
    # the artery of the deportalised lobe keeps its share
    assert r[ids.lobe_flow(ids.QHA, "lll")] > 0.0


def test_a_deportalised_lobe_is_hypoperfused(tmp_path: Path) -> None:
    r = _after(tmp_path, PortalLigation(["lll", "lml"]))
    assert r[ids.lobe_perfusion("lll")] < r[ids.lobe_perfusion("rl")]


def test_an_arterial_ligation_closes_the_artery_only(tmp_path: Path) -> None:
    r = _after(tmp_path, ArterialLigation(["lll"]))
    assert r[ids.arterial_conductance("lll")] == 0.0
    assert r[ids.portal_conductance("lll")] == 1.0
    assert r[ids.lobe_flow(ids.QHA, "lll")] == pytest.approx(0.0)
    total = sum(r[ids.lobe_flow(ids.QHA, t)] for t in RAT_TOKENS)
    assert total == pytest.approx(r[ids.QHA], rel=1e-9)


def test_a_resection_removes_the_lobe_and_leaves_the_remnant(tmp_path: Path) -> None:
    r = _after(tmp_path, Resection(["lll", "lml"]))
    assert r[ids.lobe_present("lll")] == 0.0
    assert r[ids.lobe_flow(ids.QH, "lll")] == pytest.approx(0.0)
    assert r[ids.F_REMNANT] == pytest.approx(
        1.0 - r[ids.lobe_volume_fraction("lll")] - r[ids.lobe_volume_fraction("lml")],
        rel=1e-6,
    )


def test_a_resected_lobe_metabolizes_nothing(tmp_path: Path) -> None:
    """The capacity of a lobe which is gone is 0, so the liver loses its share."""
    before = _after(tmp_path / "a")
    after = _after(tmp_path / "b", Resection(["lll", "lml"]))
    rate = f"{ids.flat(Organ.LIVER, 'DRG2M_lll')}"
    before[ids.iv_dose(S)] = 10.0
    after[ids.iv_dose(S)] = 10.0
    for r in (before, after):
        r.simulate(0.0, 60.0, 11)
    assert before[rate] > 0.0
    assert after[rate] == pytest.approx(0.0)


def test_a_resection_with_remodelling_stays_finite(tmp_path: Path) -> None:
    r = _after(tmp_path, Resection(["lll", "lml"]), remodelling=FlowDriven())
    r[ids.iv_dose(S)] = 10.0
    cve = f"[{ids.pool_species(Pool.VENOUS, S)}]"
    result = r.simulate(0.0, 5 * DAY, 201, selections=["time", cve, ids.VLI_LOBES])
    assert np.isfinite(result[:, 1]).all()
    assert np.isfinite(result[:, 2]).all()


def test_the_ligated_lobes_shrink_and_the_others_grow(tmp_path: Path) -> None:
    """The deportalised lobes atrophy, the lobes which keep their portal inflow hypertrophy."""
    r = _after(tmp_path, PortalLigation(["lll", "lml"]), remodelling=FlowDriven())
    result = r.simulate(
        0.0,
        5 * DAY,
        201,
        selections=["time", ids.lobe_volume("lll"), ids.lobe_volume("rl")],
    )
    assert result[:, 1][-1] < result[:, 1][0]
    assert result[:, 2][-1] > result[:, 2][0]


def test_every_branch_closed_leaves_the_model_finite(tmp_path: Path) -> None:
    """A liver with no inflow at all: the flows are 0 and no rate is undefined."""
    r = _after(tmp_path, PortalLigation(RAT_TOKENS), ArterialLigation(RAT_TOKENS))
    for token in RAT_TOKENS:
        assert r[ids.lobe_flow(ids.QH, token)] == pytest.approx(0.0)
    r[ids.iv_dose(S)] = 10.0
    cve = f"[{ids.pool_species(Pool.VENOUS, S)}]"
    result = r.simulate(0.0, 240.0, 61, selections=["time", cve])
    assert np.isfinite(result[:, 1]).all()


def test_a_total_portal_occlusion_strands_the_drug_in_the_portal_vein(
    tmp_path: Path,
) -> None:
    """Every portal branch closed is a liver which draws nothing from the portal vein.

    The portal pool has one outlet in the body, the inflow of the liver, and
    the splanchnic organs keep delivering into it: with the supply fractions
    the model is finite and stationary in the arithmetic sense and physically
    a dead end, the whole dose parks in the portal vein and the systemic
    exposure goes to 0. This is what the model says, so the page says it as
    well (`docs/pbpk/lobes.md`): a total portal occlusion asks for the
    portosystemic collaterals the gaps of the page name, and the model of one
    lobe without a portal inflow does not.
    """
    r = _after(tmp_path, PortalLigation(RAT_TOKENS))
    assert r[ids.QPO] > 0.0, "the portal vein of the body still carries its flow"
    for token in RAT_TOKENS:
        assert r[ids.lobe_flow(ids.QPO, token)] == pytest.approx(0.0)
    r[ids.iv_dose(S)] = 10.0
    cpo = f"[{ids.pool_species(Pool.PORTAL, S)}]"
    cve = f"[{ids.pool_species(Pool.VENOUS, S)}]"
    result = r.simulate(0.0, 2880.0, 97, selections=["time", cpo, cve])
    portal, venous = result[:, 1], result[:, 2]
    assert np.isfinite(result).all()
    # the portal pool only fills (the plateau drifts by the tolerance of the
    # integrator alone) and the systemic concentration dies away
    assert (np.diff(portal) >= -1e-6 * portal[-1]).all()
    # a 10 mg dose in the portal plasma of a rat, which is where it stays
    assert portal[-1] > 1000.0
    assert venous[-1] < 1e-6 * venous[1]


def test_the_scenarios_of_the_rat_are_those_of_the_table() -> None:
    names = scenarios(Species.RAT)
    assert set(names) == {"pvl70", "pvl90", "hepatectomy70", "hepatectomy90"}
    # the name of an intervention of the table is the name of its scenario
    assert names["pvl70"] == PortalLigation(("lll", "lml", "rml"), name="pvl70")
    assert names["hepatectomy90"] == Resection(
        ("lll", "lml", "rml", "rl"), name="hepatectomy90"
    )


def test_the_changes_of_a_scenario_name_its_lobes() -> None:
    changes = scenario_changes(Species.RAT)
    assert set(changes["pvl70"]) == {
        ids.portal_conductance(t) for t in ("lll", "lml", "rml")
    }
    assert set(changes["hepatectomy70"]) >= {
        ids.lobe_present(t) for t in ("lll", "lml", "rml")
    }


def test_the_ninety_percent_hepatectomy_of_the_rat_is_ninety_percent(
    tmp_path: Path,
) -> None:
    """The scenario of the table and the lobe table agree on the remnant."""
    r = _after(tmp_path, scenarios(Species.RAT)["hepatectomy90"])
    assert r[ids.F_REMNANT] == pytest.approx(0.10, abs=1e-9)


def test_a_scenario_of_another_species_is_refused() -> None:
    with pytest.raises(ValueError, match="right_hepatectomy"):
        scenarios(Species.RAT)["right_hepatectomy"]


def test_an_intervention_of_an_unknown_lobe_is_refused() -> None:
    with pytest.raises(ValueError, match="nope"):
        PortalLigation(["nope"]).changes(RAT_TOKENS)


def test_the_human_right_hepatectomy_leaves_the_left_lobe() -> None:
    right = scenarios(Species.HUMAN)["right_hepatectomy"]
    assert right.lobes == ("sg5", "sg6", "sg7", "sg8")


def test_a_surgery_between_two_phases_takes_effect(tmp_path: Path) -> None:
    """A ligation at a time point: the changes are applied to the running model.

    The patency of a branch is a variable of the model and no initial
    assignment reads it as a constant, so a study which ligates between two
    phases of a timecourse gets the new flows from the point it made the
    change, which is what a surgery at a time point is.
    """
    r = individual(_built(tmp_path, LobedLiver()))
    before = r[ids.lobe_flow(ids.QPO, "rl")]
    r.simulate(0.0, 60.0, 11)
    for pid, value in PortalLigation(["lll", "lml"]).changes(RAT_TOKENS).items():
        r[pid] = value
    r.simulate(60.0, 120.0, 11)
    assert r[ids.lobe_flow(ids.QPO, "lll")] == pytest.approx(0.0)
    assert r[ids.lobe_flow(ids.QPO, "rl")] > before


def test_a_species_without_lobes_has_no_changes_either() -> None:
    """`scenario_changes` answers a species without a lobe table as `scenarios` does.

    A species the table of the lobes does not cover has no procedures, which
    is an empty mapping and not an error: the two functions say the same.
    """
    assert scenarios(Species.PIG) == {}
    assert scenario_changes(Species.PIG) == {}
