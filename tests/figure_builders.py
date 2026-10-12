"""The figures of the fixture of sbmlsim 0.8.5 and of the tests built against it.

`tests/test_verify.py` builds the figures with sbmlsim 0.9, and
`scripts/figures_sbmlsim_085.py` builds them with sbmlsim 0.8.5 to write
`tests/data/figures-sbmlsim-0.8.5.json`, so this module is their one source and
imports nothing but sbmlsim (it must load in both environments).
"""

from collections.abc import Callable

from sbmlsim.experiment import SimulationExperiment
from sbmlsim.plot import Axis, Figure


class FigureStudy(SimulationExperiment):
    """An experiment the figures of `build_figures` belong to."""


def build_figures(
    experiment: SimulationExperiment,
    set_setting: Callable[[Figure, str, object], None],
) -> dict[str, Figure]:
    """Build figures as the models do, with sbmlsim 0.8.5 or 0.9.

    The fixture of sbmlsim 0.8.5 was written by this function, so it runs with
    both versions; only the setting of a figure is set differently.
    """
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
