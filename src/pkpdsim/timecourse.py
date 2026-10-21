"""The timecourse of a roadrunner model which is simulated in parts.

A model of a tissue (the tissue and perfusion models of pkdb_models, a
zonated liver) is often simulated in parts: a part sets values of the model at
its start (a perfusion, a substrate, a switch) and runs for its duration, the
next part continues from the state in which the part before ends. A
`Simulation` is the list of its parts, `run_simulation` runs them one after the
other and returns one timecourse whose time continues over the parts.

The parts run on a roadrunner instance of the model directly, without the
experiments of sbmlsim: the model is loaded, its tolerances and its selections
are set by the script which runs it.
"""

import logging
import time
from dataclasses import dataclass

import pandas as pd
import roadrunner

logger = logging.getLogger(__name__)


@dataclass
class SimPart:
    """A part of a simulation.

    Attributes:
        tend: the duration of the part, in the time unit of the model.
        changes: the values the part sets at its start, by the selection of
            roadrunner (`"[S_ext]"` for a concentration, `"S_ext"` for an
            amount or a parameter).
        steps: the number of steps of the timecourse of the part; the
            timecourse has `steps + 1` points, its start and its end included.
    """

    tend: float
    changes: dict[str, float]
    steps: int


@dataclass
class Simulation:
    """A simulation of parts which run one after the other.

    Attributes:
        sid: the id of the simulation, e.g. the key of its results.
        title: the title of the simulation, e.g. of its figures.
        parts: the parts, in the order in which they run.
    """

    sid: str
    title: str
    parts: list[SimPart]


def run_simulation(r: roadrunner.RoadRunner, sim: Simulation) -> pd.DataFrame:
    """Run the parts of a simulation and concatenate their timecourses.

    The model is reset (`resetAll`) before the first part, so that a
    simulation does not depend on the simulations before it. A part sets its
    changes and simulates from the state in which the part before ends; the
    time of its timecourse continues the time of the part before, so the last
    point of a part and the first point of the next part have the same time
    (the state before and after the changes). The columns are the
    `timeCourseSelections` of `r`.

    Args:
        r: the roadrunner instance of the model.
        sim: the simulation.

    Returns:
        The timecourse of the simulation, the timecourses of its parts one
        after the other (the index of each part starts at 0).
    """
    start = time.perf_counter()
    r.resetAll()
    dfs: list[pd.DataFrame] = []
    t_offset: float = 0.0
    for part in sim.parts:
        for key, value in part.changes.items():
            r.setValue(key, value)
        s = r.simulate(start=0, end=part.tend, steps=part.steps)
        df = pd.DataFrame(s, columns=s.colnames)
        df["time"] = df["time"] + t_offset
        t_offset = t_offset + part.tend
        dfs.append(df)

    df = pd.concat(dfs)
    logger.info(
        "simulation '%s': %d parts in %.3f s",
        sim.sid,
        len(sim.parts),
        time.perf_counter() - start,
    )
    return df
