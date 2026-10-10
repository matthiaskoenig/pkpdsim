"""The commands of a model.

A model package registers its commands in `[project.scripts]` with its model
object, e.g. for hctz:

```toml
run_hctz = "pkdb_models.models.hctz:MODEL.run_main"
fit_hctz = "pkdb_models.models.hctz:MODEL.fit_main"
report_hctz = "pkdb_models.models.hctz:MODEL.report_main"
identifiability_hctz = "pkdb_models.models.hctz:MODEL.identifiability_main"
```

`run_<model>` creates the models with the factory and runs the simulation
experiments:

```bash
run_hctz --action factory
run_hctz --action list_experiments
run_hctz --action simulate --experiments "misc,Anderson1961" --open
run_hctz --action all --results-dir /path/to/results
```

A single group is written into the directory of the group, any other selection
into `custom_selection`, both in the `simulation` folder of the results. The
output ends with the links to the report and to the figures.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Mapping, Sequence
from enum import StrEnum
from pathlib import Path
from typing import TYPE_CHECKING

from sbmlsim import display
from sbmlutils.console import console

from pkpdsim.model import resolve

if TYPE_CHECKING:
    from sbmlsim.experiment import SimulationExperiment

    from pkpdsim.model import PKPDModel


class Action(StrEnum):
    """The actions of the `run` command, `pkpdsim.model.RUN_ACTIONS`."""

    SIMULATE = "simulate"
    LIST_EXPERIMENTS = "list_experiments"
    FACTORY = "factory"
    ALL = "all"


def parse_formats(text: str) -> list[str]:
    """Get the figure formats of a comma separated list.

    Args:
        text: the formats, e.g. `"svg, PNG"`.

    Returns:
        The formats in lower case, without blanks and empty entries.
    """
    return [part.strip().lower() for part in text.split(",") if part.strip()]


def resolve_experiments(
    registry: Mapping[str, Sequence[type[SimulationExperiment]]],
    names: Sequence[str],
) -> tuple[list[type[SimulationExperiment]], list[str]]:
    """Resolve the names of groups and experiments to the experiment classes.

    Args:
        registry: the groups of experiments of a model.
        names: names of groups or of experiment classes.

    Returns:
        The experiment classes, every class once in the order of the names, and
        the names which are neither a group nor an experiment.
    """
    experiments = {
        experiment.__name__: experiment
        for group in registry.values()
        for experiment in group
    }
    classes: dict[str, type[SimulationExperiment]] = {}
    not_found: list[str] = []
    for name in names:
        if name in registry:
            classes.update({e.__name__: e for e in registry[name]})
        elif name in experiments:
            classes[name] = experiments[name]
        else:
            not_found.append(name)
    return list(classes.values()), not_found


def run_parser(model: PKPDModel) -> argparse.ArgumentParser:
    """Get the parser of the `run` command of a model."""
    actions = [action.value for action in Action]
    actions += [extra.name for extra in model.extra_actions]
    action_help = "the action to perform" + "".join(
        f"; {extra.name}: {extra.help}" for extra in model.extra_actions
    )
    p = argparse.ArgumentParser(
        prog=f"run_{model.name}",
        description="Run the model factory and the simulation experiments of the "
        f"{model.title}.",
    )
    p.add_argument("-a", "--action", required=True, choices=actions, help=action_help)
    p.add_argument(
        "-r",
        "--results-dir",
        type=Path,
        default=model.results_dir,
        help="directory of the results, the simulations are written to its "
        "`simulation` folder (default: the results of the model)",
    )
    p.add_argument(
        "-e",
        "--experiments",
        help="comma separated groups and experiments to simulate, see "
        "`--action list_experiments`",
    )
    p.add_argument(
        "-o",
        "--open",
        action="store_true",
        help="open the report of the simulations in the web browser",
    )
    p.add_argument(
        "-f",
        "--formats",
        default=",".join(model.figures.formats),
        help="comma separated formats of the figures: static images drawn by "
        "matplotlib (e.g. 'svg', 'png') and interactive pages drawn by plotly "
        "('html'). The report shows the SVG and links the interactive page; "
        "without 'svg' the interactive page is embedded in the report "
        "(default: %(default)s)",
    )
    if model.extra_arguments is not None:
        resolve(model.extra_arguments, model.package)(p)
    return p


def _run_factory(model: PKPDModel) -> None:
    """Create the SBML models, the ODEs and the COMBINE archive."""
    display.section("Model factory", icon=":factory:")
    model.create_models()
    display.link("models", model.factory_dir)
    if model.archive_path is not None:
        display.link("archive", model.archive_path)


def _list_experiments(model: PKPDModel) -> None:
    """Print the groups and the experiments which can be simulated."""
    registry = model.load_experiments()
    display.section("Simulation experiments", icon=":test_tube:")
    display.key_values(
        {
            "groups": ", ".join(registry),
            **{
                group: ", ".join(e.__name__ for e in experiments)
                for group, experiments in registry.items()
                if group != Action.ALL
            },
        }
    )


def _simulate(
    model: PKPDModel,
    names: Sequence[str],
    simulation_dir: Path,
    show_report: bool,
    figure_formats: Sequence[str],
) -> None:
    """Simulate the given groups and experiments.

    Args:
        model: the model.
        names: names of groups or experiments.
        simulation_dir: directory the results of the simulations are written to.
        show_report: open the report in a web browser.
        figure_formats: the formats of the figures.
    """
    registry = model.load_experiments()
    classes, not_found = resolve_experiments(registry, names)
    if not_found:
        console.print(f"Not found: {', '.join(not_found)}", style="warning")
    if not classes:
        console.print("No experiments to simulate.", style="error")
        sys.exit(1)

    single_group = len(names) == 1 and names[0] in registry
    output_dir = simulation_dir / (names[0] if single_group else "custom_selection")
    display.section("Simulations", icon=":chart_increasing:")
    display.key_values({"experiments": len(classes)})
    display.link("output", output_dir)
    model.run_simulation_experiments(
        classes,
        output_dir,
        show_report=show_report,
        figure_formats=figure_formats,
    )


def run_command(model: PKPDModel, args: Sequence[str] | None = None) -> None:
    """Run the action of the command line of `run_<model>`.

    Args:
        model: the model.
        args: the command line arguments, `sys.argv[1:]` by default.
    """
    p = run_parser(model)
    options = p.parse_args(args)
    formats = parse_formats(options.formats)
    if not formats:
        p.error("'--formats' requires at least one format")
    names = [
        name.strip() for name in (options.experiments or "").split(",") if name.strip()
    ]
    if options.action == Action.SIMULATE and not names:
        p.error("'--action simulate' requires '--experiments'")

    display.section(model.title, icon=":pill:")
    display.key_values({"action": options.action})
    display.link("results", options.results_dir)
    results_dir = options.results_dir.resolve()
    simulation_dir = results_dir / "simulation"
    extras = {extra.name: extra for extra in model.extra_actions}

    if options.action == Action.FACTORY:
        _run_factory(model)
    elif options.action == Action.LIST_EXPERIMENTS:
        _list_experiments(model)
    elif options.action == Action.SIMULATE:
        _simulate(model, names, simulation_dir, options.open, formats)
    elif options.action == Action.ALL:
        _run_factory(model)
        _simulate(model, ["all"], simulation_dir, options.open, formats)
        for extra in model.extra_actions:
            if extra.in_all:
                resolve(extra.function, model.package)(model, options, results_dir)
    else:
        extra = extras[options.action]
        resolve(extra.function, model.package)(model, options, results_dir)
