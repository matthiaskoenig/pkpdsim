"""Write `tests/data/figures-sbmlsim-0.8.5.json`, the figures of sbmlsim 0.8.5.

The fixture is the serialization (`pkpdsim.verify.serialize`) of the figures of
`build_figures` of `tests/test_verify.py`, as pkpdsim 0.4.1 with sbmlsim 0.8.5
stored them; the tests compare it with the figures sbmlsim 0.9 builds. sbmlsim
0.8.5 cannot be installed next to the environment of the repository, so this
script is not part of the test suite: it needs an environment with sbmlsim
0.8.5 and pkpdsim 0.4.1 and is run by hand, only when `build_figures` changes:

    uv venv venv085 && uv pip install --python venv085 "sbmlsim==0.8.5" "pkpdsim==0.4.1"
    venv085/bin/python -I scripts/figures_sbmlsim_085.py tests/data/figures-sbmlsim-0.8.5.json

`build_figures` is a copy of the function of `tests/test_verify.py` (that
module imports sbmlsim 0.9); change both together. Run again with the same
versions, the output is the committed file byte for byte.
"""

import json
import logging
import sys
from collections.abc import Callable

import sbmlsim
from sbmlsim.experiment import SimulationExperiment
from sbmlsim.plot import Axis, Figure

import pkpdsim
from pkpdsim.verify import serialize

logging.disable(logging.CRITICAL)


class FigureStudy(SimulationExperiment):
    """An experiment the figures belong to."""


def build_figures(
    experiment: SimulationExperiment, set_setting: Callable[[Figure, str, object], None]
) -> dict[str, Figure]:
    figures: dict[str, Figure] = {}
    # a setting and an explicit width set after the creation (Harris1992 of creatinine)
    fig = Figure(experiment=experiment, sid="Fig1", name="FigureStudy Fig1")
    set_setting(fig, "legend_position", "outside")
    fig.width = 10
    plots = fig.create_plots(xaxis=Axis("time", unit="hr"), legend=True)
    plots[0].set_yaxis("creatine", unit="mM")
    plots[0].add_data(task="task_cr", xid="time", yid="[Cve_cr]", label="Sim")
    plots[0].add_data(dataset="d", xid="time", yid="mean", yid_sd="mean_sd", label="d")
    figures["Fig1"] = fig
    # an explicit size of 2 x 3 panels (Kharasch1997 of midazolam)
    fig = Figure(experiment, sid="Fig2", num_rows=2, num_cols=3, height=15, width=25)
    plots = fig.create_plots(xaxis=Axis("time", unit="min"), legend=True)
    plots[0].add_data(task="task_a", xid="time", yid="[A]", label="A")
    figures["Fig2"] = fig
    # the size of the panels, 1 x 2, and a resolution (Bedada2016 of chlorzoxazone)
    fig = Figure(experiment, sid="Fig3", num_cols=2)
    set_setting(fig, "fig_dpi", 300)
    plots = fig.create_plots(xaxis=Axis("time", unit="min"), legend=True)
    plots[1].add_data(task="task_a", xid="time", yid="[A]", label="A")
    figures["Fig3"] = fig
    return figures


figures = build_figures(FigureStudy(), setattr)
content = {
    "sbmlsim": sbmlsim.__version__,
    "pkpdsim": pkpdsim.__version__,
    "figures": {
        sid: json.dumps(serialize(figure), sort_keys=True, default=repr)
        for sid, figure in figures.items()
    },
}
with open(sys.argv[1], "w", encoding="utf-8") as f:
    json.dump(content, f, indent=1, sort_keys=True)
    f.write("\n")
