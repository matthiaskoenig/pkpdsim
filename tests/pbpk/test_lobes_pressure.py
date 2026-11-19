"""Tests of the lobes of the liver on the pressure driven circulation.

With `PressureFlow` the hepatic flows come from a resistance network and not
from the flow fractions of the individual, so a lobed liver gets a sinusoidal
node of its own per lobe and three branches to it: the hepatic artery, the
portal vein and the hepatic vein. The balance of the three at every node
solves the pressure of the lobe, so the flows of the lobes follow the anatomy
of the tree instead of being shares of a prescribed total, which is what makes
a ligation raise the portal pressure and the hepatic arterial buffer response
possible.

Two invariants hold whatever the network does: the flows of the lobes sum to
the flows of the liver, which every consumer of the liver reads, and at
baseline the network agrees with the supply fraction mode of a liver without
pressures. Both are tests here.

The models carry algebraic rules, so they are solved by SUNDIALS IDA of the
core engine of sbmlsim and not by roadrunner (`cardiovascular_solver`); the
first row of a result is the solved network of the initialized individual,
which is what an intervention applied before the initialization is read from.
"""

from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pytest

from pkpdsim.pbpk import (
    DynamicHeart,
    FirstOrder,
    Metabolism,
    Organ,
    PBPKDefinition,
    PressureFlow,
    Route,
    Variants,
    ids,
)
from pkpdsim.pbpk.builder import build
from pkpdsim.pbpk.interventions import ArterialLigation, Intervention, PortalLigation
from pkpdsim.pbpk.lobes import HABR
from pkpdsim.pbpk.remodelling import FlowDriven
from pkpdsim.pbpk.variants.liver_lobes import LobedLiver
from pkpdsim.species import Species

from .cardiovascular_solver import simulate
from .conftest import definition
from .test_lobes_build import RAT_TOKENS

S = "drg"
PRESSURES = ("CV_P_ar", "CV_P_po", "CV_P_ve", "CV_P_li")
LOBE_FLOWS = tuple(
    ids.lobe_flow(flow, token) for flow in ids.HEPATIC_FLOWS for token in RAT_TOKENS
)
NODES = tuple(ids.lobe_pressure(token) for token in RAT_TOKENS)
BUFFERS = tuple(ids.buffer_factor(token) for token in RAT_TOKENS)
EVERYTHING = (*PRESSURES, *ids.HEPATIC_FLOWS, *LOBE_FLOWS, *NODES, *BUFFERS)


def _definition(
    liver: LobedLiver | None, circulation: PressureFlow | None
) -> PBPKDefinition:
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
        variants=Variants(
            organs={Organ.LIVER: liver} if liver is not None else {},
            circulation=circulation,
        ),
    )


def _model(
    tmp_path: Path, liver: LobedLiver | None, circulation: PressureFlow | None
) -> Path:
    return build(_definition(liver, circulation), tmp_path).flat_path


def _solved(
    path: Path,
    selections: Sequence[str],
    interventions: Sequence[Intervention] = (),
) -> dict[str, float]:
    """The solved network of the individual, with the changes of the interventions."""
    changes: dict[str, float] = {}
    for intervention in interventions:
        changes |= intervention.changes(RAT_TOKENS)
    result = simulate(path, list(selections), preinit_changes=changes or None)
    return {sid: float(result[sid][0]) for sid in selections}


def _pressure(
    tmp_path: Path,
    *interventions: Intervention,
    selections: Sequence[str] = EVERYTHING,
    habr: HABR | None = None,
) -> dict[str, float]:
    """A lobed liver on the pressure network, solved after the interventions."""
    path = _model(tmp_path, LobedLiver(habr=habr), PressureFlow())
    return _solved(path, selections, interventions)


def test_the_lobe_flows_sum_to_the_liver_flows(tmp_path: Path) -> None:
    v = _pressure(tmp_path)
    for flow in ids.HEPATIC_FLOWS:
        total = sum(v[ids.lobe_flow(flow, t)] for t in RAT_TOKENS)
        assert total == pytest.approx(v[flow], rel=1e-7), flow


def test_every_lobe_has_its_own_sinusoidal_node(tmp_path: Path) -> None:
    v = _pressure(tmp_path)
    for token in RAT_TOKENS:
        assert v["CV_P_ve"] < v[ids.lobe_pressure(token)] < v["CV_P_po"]


def test_the_balance_of_every_node_holds(tmp_path: Path) -> None:
    """The three branches of a lobe balance, which is what solves its pressure."""
    v = _pressure(tmp_path)
    for token in RAT_TOKENS:
        inflow = v[ids.lobe_flow(ids.QHA, token)] + v[ids.lobe_flow(ids.QPO, token)]
        assert inflow == pytest.approx(v[ids.lobe_flow(ids.QH, token)], rel=1e-7)


def test_the_pressure_and_the_supply_fraction_modes_agree_at_baseline(
    tmp_path: Path,
) -> None:
    """A network whose branches are scaled by the lobe shares gives those shares."""
    selections = (*ids.HEPATIC_FLOWS, *LOBE_FLOWS)
    network = _pressure(tmp_path / "a", selections=selections)
    shares = _solved(_model(tmp_path / "b", LobedLiver(), None), selections)
    for flow in ids.HEPATIC_FLOWS:
        for token in RAT_TOKENS:
            of_network = network[ids.lobe_flow(flow, token)] / network[flow]
            of_shares = shares[ids.lobe_flow(flow, token)] / shares[flow]
            assert of_network == pytest.approx(of_shares, rel=1e-6), (flow, token)


def test_a_portal_ligation_raises_the_portal_pressure(tmp_path: Path) -> None:
    """The splanchnic bed delivers its blood through fewer branches."""
    path = _model(tmp_path, LobedLiver(), PressureFlow())
    before = _solved(path, PRESSURES)
    after = _solved(path, PRESSURES, [PortalLigation(["lll", "lml"])])
    assert after["CV_P_po"] > before["CV_P_po"]


def test_the_total_portal_flow_falls_only_a_little_after_a_ligation(
    tmp_path: Path,
) -> None:
    """The portal bed is a pressure source with a resistance, not a prescribed flow."""
    path = _model(tmp_path, LobedLiver(), PressureFlow())
    before = _solved(path, (ids.QPO,))
    after = _solved(path, (ids.QPO,), [PortalLigation(["lll", "lml"])])
    assert after[ids.QPO] < before[ids.QPO]
    assert after[ids.QPO] == pytest.approx(before[ids.QPO], rel=0.4)


def test_a_ligated_lobe_carries_no_portal_flow(tmp_path: Path) -> None:
    v = _pressure(tmp_path, PortalLigation(["lll"]))
    assert v[ids.lobe_flow(ids.QPO, "lll")] == pytest.approx(0.0, abs=1e-12)
    assert v[ids.lobe_flow(ids.QHA, "lll")] > 0.0


def test_a_devascularised_lobe_sits_at_the_venous_pressure(tmp_path: Path) -> None:
    """Both branches closed: the node has no inflow, so it is at the venous pressure."""
    v = _pressure(tmp_path, PortalLigation(["lll"]), ArterialLigation(["lll"]))
    assert v[ids.lobe_pressure("lll")] == pytest.approx(v["CV_P_ve"], rel=1e-6)
    assert v[ids.lobe_flow(ids.QH, "lll")] == pytest.approx(0.0, abs=1e-12)


def test_a_devascularised_lobe_leaves_the_solution_finite(tmp_path: Path) -> None:
    path = _model(tmp_path, LobedLiver(), PressureFlow())
    result = simulate(
        path,
        [ids.QH, *NODES],
        end=60.0,
        steps=60,
        preinit_changes=PortalLigation(["lll"]).changes(RAT_TOKENS)
        | ArterialLigation(["lll"]).changes(RAT_TOKENS),
    )
    for sid in (ids.QH, *NODES):
        assert np.isfinite(result[sid]).all(), sid


def test_the_buffer_response_raises_the_arterial_flow_of_a_ligated_lobe(
    tmp_path: Path,
) -> None:
    path = _model(tmp_path, LobedLiver(habr=HABR()), PressureFlow())
    arterial = (ids.lobe_flow(ids.QHA, "lll"), ids.buffer_factor("lll"))
    before = _solved(path, arterial)
    after = _solved(path, arterial, [PortalLigation(["lll"])])
    assert after[ids.buffer_factor("lll")] > before[ids.buffer_factor("lll")]
    assert after[arterial[0]] > before[arterial[0]]


def test_the_buffer_response_compensates_the_measured_fraction(tmp_path: Path) -> None:
    """Lautt 1985: the hepatic artery compensates 25.5 +- 2.7 % of a lost portal flow.

    The gain is the compensation of the conductance of the branch, so the flow
    which arrives also follows the pressure of the node, which the lost portal
    inflow lowers; the compensated fraction is therefore of the order of the
    gain and not the gain itself.
    """
    path = _model(tmp_path, LobedLiver(habr=HABR()), PressureFlow())
    selections = (ids.lobe_flow(ids.QHA, "lll"), ids.lobe_flow(ids.QPO, "lll"))
    before = _solved(path, selections)
    after = _solved(path, selections, [PortalLigation(["lll"])])
    lost = before[selections[1]]
    gained = after[selections[0]] - before[selections[0]]
    assert 0.5 * HABR().gain.value < gained / lost < 2.0 * HABR().gain.value


def test_without_the_buffer_response_the_artery_only_follows_its_node(
    tmp_path: Path,
) -> None:
    """Without it a deportalised lobe gains far less than the measured fraction."""
    path = _model(tmp_path, LobedLiver(), PressureFlow())
    selections = (ids.lobe_flow(ids.QHA, "lll"), ids.lobe_flow(ids.QPO, "lll"))
    before = _solved(path, selections)
    after = _solved(path, selections, [PortalLigation(["lll"])])
    gained = after[selections[0]] - before[selections[0]]
    assert 0.0 < gained / before[selections[1]] < 0.5 * HABR().gain.value


def test_the_buffer_response_is_capped(tmp_path: Path) -> None:
    """The artery of a lobe cannot grow beyond the ceiling, whatever it would compensate.

    Deportalising every lobe but the caudate is the strongest ligation the
    pressure network has a solution for: the portal inflow of a lobe is a
    multiple of its arterial inflow, so an uncapped response would raise the
    conductance of a deportalised artery by far more than the ceiling.
    """
    ligated = [t for t in RAT_TOKENS if t != "cl"]
    path = _model(tmp_path, LobedLiver(habr=HABR(maximum=1.1)), PressureFlow())
    v = _solved(path, BUFFERS, [PortalLigation(ligated)])
    for token in ligated:
        assert v[ids.buffer_factor(token)] == pytest.approx(1.1, rel=1e-9)
    assert v[ids.buffer_factor("cl")] == pytest.approx(1.0, rel=1e-6)


def test_the_uncapped_response_is_the_gain_on_the_lost_conductance(
    tmp_path: Path,
) -> None:
    """The factor is 1 + gain * (lost portal flow)/(reference arterial flow)."""
    path = _model(tmp_path, LobedLiver(habr=HABR(maximum=10.0)), PressureFlow())
    selections = (
        ids.buffer_factor("lll"),
        ids.lobe_reference_flow(ids.QPO, "lll"),
        ids.lobe_reference_flow(ids.QHA, "lll"),
    )
    v = _solved(path, selections, [PortalLigation(["lll"])])
    expected = 1.0 + HABR().gain.value * v[selections[1]] / v[selections[2]]
    assert v[selections[0]] == pytest.approx(expected, rel=1e-7)


def test_every_portal_branch_closed_has_no_solution_on_the_network(
    tmp_path: Path,
) -> None:
    """The splanchnic bed has nowhere to drain, which the solver cannot solve.

    With a prescribed flow the same occlusion is simply a liver without
    portal inflow (`test_lobes_interventions`), but the pressure network has
    to send the flow of the portal organs somewhere: the portal junction
    balance has no solution once every branch of the liver is closed, which
    is what a total portal occlusion without a portosystemic shunt is.
    """
    path = _model(tmp_path, LobedLiver(), PressureFlow())
    with pytest.raises(RuntimeError, match="integration failed"):
        _solved(path, (ids.QPO,), [PortalLigation(RAT_TOKENS)])


def test_a_liver_without_lobes_keeps_the_network_it_had(tmp_path: Path) -> None:
    """The one node, the two hepatic resistances and the one balance of the liver."""
    text = _model(tmp_path, None, PressureFlow()).read_text()
    for sid in ("CV_hepatic_balance", "CV_R_portal", "CV_R_hepatic"):
        assert f'id="{sid}"' in text or f'"{sid}"' in text, sid
    for token in RAT_TOKENS:
        assert ids.lobe_pressure(token) not in text
        assert ids.lobe_balance(token) not in text


def test_a_lobed_liver_replaces_the_one_node_by_the_nodes_of_the_lobes(
    tmp_path: Path,
) -> None:
    """The whole liver resistances and balance are the lobes' now."""
    text = _model(tmp_path, LobedLiver(), PressureFlow()).read_text()
    assert ids.lobe_balance("lll") in text
    assert '"CV_hepatic_balance"' not in text
    assert '"CV_R_portal"' not in text
    assert '"CV_R_hepatic"' not in text


def test_a_dynamic_heart_refuses_a_remodelling(tmp_path: Path) -> None:
    """A beating heart needs a constant anatomy, which a remodelling is not."""
    broken = definition(
        sid="m",
        species=Species.HUMAN,
        routes=[Route.iv(S)],
        variants=Variants(
            organs={Organ.LIVER: LobedLiver(remodelling=FlowDriven())},
            circulation=DynamicHeart(),
        ),
    )
    with pytest.raises(ValueError, match="constant anatomy"):
        build(broken, tmp_path)


def test_the_sinusoidal_pressure_of_the_liver_is_the_mean_of_the_lobes(
    tmp_path: Path,
) -> None:
    """`CV_P_li` keeps its id and is the flow weighted mean of the lobe nodes."""
    v = _pressure(tmp_path)
    weighted = sum(
        v[ids.lobe_flow(ids.QH, t)] * v[ids.lobe_pressure(t)] for t in RAT_TOKENS
    )
    assert v["CV_P_li"] == pytest.approx(weighted / v[ids.QH], rel=1e-7)
