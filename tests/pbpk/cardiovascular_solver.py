"""Simulate cardiovascular SBML with the public sbmlsim DAE engine."""

from collections.abc import Mapping, Sequence
from pathlib import Path

from sbmlsim.model import AbstractModel
from sbmlsim.result.timecourse import TimecourseResult
from sbmlsim.simulation import Change, Simulation
from sbmlsim.simulator import Simulator


def simulate(
    path: Path,
    selections: Sequence[str],
    *,
    end: float = 1.0,
    steps: int = 100,
    preinit_changes: Mapping[str, float] | None = None,
    changes: Sequence[Change] | None = None,
    maximum_time_step: float = 0.0,
) -> TimecourseResult:
    """Solve a cardiovascular model in its native time unit, the minute.

    Algebraic rules select SUNDIALS IDA automatically in the core engine.
    Selecting columns explicitly avoids unrelated pharmacokinetic outputs.
    """
    return Simulator(
        engine="core",
        n_workers=1,
        absolute_tolerance=1e-10,
        relative_tolerance=1e-8,
        maximum_time_step=maximum_time_step,
    ).simulate(
        AbstractModel(source=path, selections=["time", *selections]),
        Simulation(
            end=end,
            steps=steps,
            preinit_changes=preinit_changes,
            changes=changes,
        ),
    )
