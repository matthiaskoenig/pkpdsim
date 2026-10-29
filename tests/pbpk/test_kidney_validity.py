"""The validity of the kidney outflow (`pkpdsim.pbpk.experiment.kidney_outflow`).

The plasma leaving the kidney lacks the net filtrate, so the outflow is valid
while the net filtered plasma flow of a substance, `fu*GFR*(1 - f_reabsorbed)`
(`fu*GFR` with a saturable reabsorption), is below the plasma flow of the
kidney, `Qki*(1 - HCT)` (`Qki*BP` with erythrocytes in equilibrium).
`PBPKSimulationExperiment` refuses an experiment which simulates an individual
outside of it and a renal function beyond it, and derives the range of the
renal function from it.
"""

from functools import cache
from pathlib import Path
from typing import cast

import numpy as np
import pytest
from sbmlsim import Q
from sbmlsim.model.model_roadrunner import RoadrunnerSBMLModel
from sbmlsim.simulation import Dimension, Scan, Simulation
from sbmlsim.simulation.definition import Value as Change
from sbmlsim.simulator import Simulator

from pkpdsim.model import PKPDModel
from pkpdsim.pbpk import (
    Binding,
    Erythrocytes,
    Excretion,
    Organ,
    PBPKDefinition,
    Route,
    Saturable,
    Substance,
    build,
    ids,
)
from pkpdsim.pbpk.builder import BuildResult
from pkpdsim.pbpk.experiment import (
    KidneyOutflow,
    PBPKSimulationExperiment,
    kidney_outflow,
)
from pkpdsim.uncertainty import UNCERTAINTY_DIMENSION

from .conftest import individual

FU = 0.9


def _definition() -> PBPKDefinition:
    """Four filtered substances, one for each carrier and reabsorption.

    `pls` without erythrocytes (the plasma carries it, `BP = 1 - HCT`), `kin`
    with kinetic erythrocytes (the plasma carries it), `equ` in equilibrium
    with the erythrocytes (`BP` carries it), `rea` with a constant fraction
    reabsorbed and `sat` with a saturable reabsorption (`fu*GFR` counts).
    """
    return PBPKDefinition(
        sid="validity",
        substances=[
            Substance(
                "pls",
                name="plasma",
                mr=100.0,
                binding=Binding(fu=FU),
                excretion=Excretion(urine=True),
            ),
            Substance(
                "kin",
                name="kinetic",
                mr=100.0,
                binding=Binding(fu=FU),
                erythrocytes=Erythrocytes("kinetic", kbc=2.0, k_exchange=1.0),
                excretion=Excretion(urine=True),
            ),
            Substance(
                "equ",
                name="equilibrium",
                mr=100.0,
                binding=Binding(fu=FU),
                erythrocytes=Erythrocytes("equilibrium", bp=0.8),
                excretion=Excretion(urine=True),
            ),
            Substance(
                "rea",
                name="reabsorbed",
                mr=100.0,
                binding=Binding(fu=FU),
                excretion=Excretion(urine=True, fraction_reabsorbed=0.5),
            ),
            Substance(
                "sat",
                name="saturable",
                mr=100.0,
                binding=Binding(fu=FU),
                excretion=Excretion(urine=True, reabsorption=Saturable(0.1, 2.0)),
            ),
            Substance("non", name="not filtered", mr=100.0),
        ],
        routes=[Route.iv("pls")],
    )


@pytest.fixture(scope="module")
def built(tmp_path_factory: pytest.TempPathFactory) -> BuildResult:
    """The definition of the four carriers, built once."""
    return build(_definition(), tmp_path_factory.mktemp("validity"))


def _experiment(result: BuildResult) -> type[PBPKSimulationExperiment]:
    """The base of the experiments of the model."""

    class ValidityExperiment(PBPKSimulationExperiment):
        pkpd_model = PKPDModel(
            name="validity",
            package="validity",
            title="Validity",
            path=result.flat_path.parent,
            model_path=result.flat_path,
            factory=".factory:create_models",
        )

    return ValidityExperiment


def _study(
    base: type[PBPKSimulationExperiment],
    name: str,
    simulations: dict[str, Simulation | Scan],
) -> PBPKSimulationExperiment:
    """An experiment of the model with the given simulations, not initialized."""

    def simulations_of(self: PBPKSimulationExperiment) -> dict[str, Simulation | Scan]:
        return simulations

    cls = cast(
        "type[PBPKSimulationExperiment]",
        type(name, (base,), {"simulations": simulations_of}),
    )
    experiment = cls(base_path=base.pkpd_model.path)
    # the runner resolves the models before the initialization
    experiment._models = {"model": _model(base.pkpd_model.model_path)}
    return experiment


@cache
def _model(path: Path) -> RoadrunnerSBMLModel:
    """The loaded model of the experiments, once."""
    return RoadrunnerSBMLModel(source=path)


def _renal(f: float) -> dict[str, Change]:
    return {ids.F_RENAL_FUNCTION: Q(f, "dimensionless")}


def test_the_flows_of_every_carrier(built: BuildResult) -> None:
    """The net filtered plasma flow and the plasma flow of each substance, in ml/min."""
    changes = {"BW": 90.0, "HF_severity": 2.0, ids.F_RENAL_FUNCTION: 1.4}
    r = individual(built.flat_path, **changes)
    gfr = r[ids.GFR] * 1000
    flow = r[ids.organ_flow(Organ.KIDNEY)] * 1000
    plasma = flow * (1 - r[ids.HCT])
    quantities = {
        "BW": Q(90.0, "kg"),
        "HF_severity": Q(2.0, "dimensionless"),
        ids.F_RENAL_FUNCTION: Q(1.4, "dimensionless"),
    }
    outflows = {o.substance: o for o in kidney_outflow(built.flat_path, quantities)}
    assert set(outflows) == {"pls", "kin", "equ", "rea", "sat"}
    expected = {
        "pls": (FU * gfr, plasma),
        "kin": (FU * gfr, plasma),
        "equ": (FU * gfr, flow * 0.8),
        "rea": (FU * gfr * 0.5, plasma),
        "sat": (FU * gfr, plasma),
    }
    for s, (filtered, plasma_flow) in expected.items():
        assert outflows[s].gfr == pytest.approx(gfr, rel=1e-12)
        assert outflows[s].filtered == pytest.approx(filtered, rel=1e-12)
        assert outflows[s].plasma_flow == pytest.approx(plasma_flow, rel=1e-12)
        assert outflows[s].valid


def test_an_outflow_is_valid_below_the_plasma_flow() -> None:
    assert KidneyOutflow("x", gfr=100.0, filtered=99.0, plasma_flow=100.0).valid
    assert not KidneyOutflow("x", gfr=100.0, filtered=100.0, plasma_flow=100.0).valid


def test_the_range_of_the_renal_function_is_derived(built: BuildResult) -> None:
    """The highest renal function is the one at which `fu*GFR` reaches the plasma flow.

    The substances which the plasma carries (`pls`, `kin`, `sat`, `fu*GFR`
    against `Qki*(1 - HCT)`) limit it; the blood to plasma ratio 0.8 of `equ`
    is above the plasma share `1 - HCT`.
    """
    cls = _experiment(built)
    r = individual(built.flat_path, BW=60.0, HEIGHT=165.0)
    plasma = r[ids.organ_flow(Organ.KIDNEY)] * min(1 - r[ids.HCT], 0.8)
    expected = plasma / (FU * r[ids.GFR_BASE])
    low, high = cls.renal_function_limits(bw=60.0, height=165.0)
    assert low == 0.0
    assert high == pytest.approx(expected, rel=1e-12)
    assert 2.0 < high < 10.0
    healthy = float(cls.healthy_gfr(bw=60.0, height=165.0).m_as("ml/min"))
    below = cls.renal_function(0.999 * high * healthy, bw=60.0, height=165.0)
    assert below.magnitude == pytest.approx(0.999 * high)
    with pytest.raises(
        ValueError,
        match=(
            r"outside the validity of the kidney outflow: the net filtered plasma "
            r"flow of '\w+' fu\*GFR\*\(1 - f_reabsorbed\) [\d.]+ ml/min at the GFR "
            r"[\d.]+ ml/min against the plasma flow of the kidney [\d.]+ ml/min"
        ),
    ):
        cls.renal_function(1.001 * high * healthy, bw=60.0, height=165.0)


def test_an_individual_outside_the_validity_fails_the_initialization(
    built: BuildResult,
) -> None:
    """A simulation, a point of a scan: the message names experiment, GFR and plasma flow."""
    base = _experiment(built)
    high = base.renal_function_limits()[1]
    valid = _study(
        base,
        "ValidStudy",
        {
            "high": Simulation(end=10, preinit_changes=_renal(0.99 * high)),
            "scan": Scan(
                Simulation(end=10),
                [
                    Dimension(
                        "renal", values={ids.F_RENAL_FUNCTION: [0.1, 1.0, 0.99 * high]}
                    )
                ],
            ),
        },
    )
    valid.initialize()
    simulation = _study(
        base,
        "HyperfiltrationStudy",
        {"high": Simulation(end=10, preinit_changes=_renal(1.01 * high))},
    )
    with pytest.raises(
        ValueError,
        match=(
            r"The experiment 'HyperfiltrationStudy' simulates in 'high' an individual "
            r"outside the validity of the kidney outflow: the net filtered plasma flow "
            r"of '\w+' .* at the GFR [\d.]+ ml/min against the plasma flow of the "
            r"kidney [\d.]+ ml/min"
        ),
    ):
        simulation.initialize()
    scan = _study(
        base,
        "ScanStudy",
        {
            "scan": Scan(
                Simulation(end=10),
                [
                    Dimension(
                        "renal",
                        values={ids.F_RENAL_FUNCTION: np.linspace(0.5, 1.01 * high, 5)},
                    )
                ],
            )
        },
    )
    with pytest.raises(ValueError, match=r"'ScanStudy' simulates in 'scan'"):
        scan.initialize()


def test_doses_and_draws_make_no_individual(
    built: BuildResult, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Simulations of one individual are simulated once; the draws are not checked."""
    base = _experiment(built)
    high = base.renal_function_limits()[1]
    runs: list[object] = []
    run = Simulator.run

    def counted(self: Simulator, *args: object, **kwargs: object) -> object:
        runs.append(args)
        return run(self, *args, **kwargs)  # ty: ignore[invalid-argument-type]

    monkeypatch.setattr(Simulator, "run", counted)
    bw = {"BW": Q(81.2345, "kg")}
    study = _study(
        base,
        "DoseStudy",
        {
            f"dose{k}": Simulation(
                end=10, preinit_changes={**bw, ids.iv_dose("pls"): Q(k, "mg")}
            )
            for k in (1, 2, 3)
        }
        | {
            "doses": Scan(
                Simulation(end=10, preinit_changes=bw),
                [Dimension("dose", values={ids.iv_dose("pls"): [5.0, 6.0, 7.0]})],
            ),
            "draws": Scan(
                Simulation(end=10, preinit_changes=bw),
                [
                    Dimension(
                        UNCERTAINTY_DIMENSION,
                        values={ids.F_RENAL_FUNCTION: [1.0, 2 * high]},
                    )
                ],
            ),
        },
    )
    study.initialize()
    assert len(runs) == 1


def test_a_model_without_filtration_is_always_valid(tmp_path: Path) -> None:
    """Without a filtered substance there is no limit."""
    definition = PBPKDefinition(
        sid="no_urine",
        substances=[Substance("non", name="not filtered", mr=100.0)],
        routes=[Route.iv("non")],
    )
    result = build(definition, tmp_path)
    assert kidney_outflow(result.flat_path) == ()
    assert _experiment(result).renal_function_limits() == (0.0, float("inf"))
