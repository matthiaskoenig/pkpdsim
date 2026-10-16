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

`--uncertainty` runs the experiments over draws of the fitted parameters of the
model and draws the curves of the simulations as bands of the quantiles of the
draws, into the `simulation_uncertainty` folder of the results (see
`pkpdsim.uncertainty`):

```bash
run_hctz --action simulate --experiments all --uncertainty --draws 100
```

`fit_<model>` fits the parameters of a fit problem and, with
`--identifiability` (the default of a model with `fit_identifiability`),
analyses the identifiability of the fitted parameters: the profile likelihood
of every fitted parameter and the Fisher information around the best parameter
set, which end up in the identifiability section of the report of the fit.
`fit_<model> --adopt <fit id>` copies the record of a finished fit into the
adopted fits of the model, whose values and uncertainty the factory takes
(`pkpdsim.uncertainty.adopt_fit`).
`report_<model>` and `identifiability_<model>` work on the `parameters.json`
of an earlier fit (`sbmlsim.fit.cli.report_cli`, `identifiability_cli`). All
write into the results of the model unless `--output_dir` is given.
"""

from __future__ import annotations

import argparse
import dataclasses
import sys
from collections.abc import Mapping, Sequence
from enum import StrEnum
from pathlib import Path
from typing import TYPE_CHECKING

from sbmlsim import display, log
from sbmlsim.experiment import ExperimentRunError
from sbmlsim.fit import display as fit_display
from sbmlsim.fit.cli import (
    ALGORITHMS,
    STRATEGY_INFO,
    FitDefinition,
    FitRun,
    identifiability_cli,
    report_cli,
    run_fit,
)
from sbmlsim.fit.fisher import FisherInformation
from sbmlsim.fit.identifiability import IdentifiabilityResult, ProfileSettings
from sbmlsim.fit.options import OptimizationStrategy
from sbmlsim.fit.result import fit_id
from sbmlsim.utils import paths_text
from sbmlutils.console import console

from pkpdsim.model import resolve
from pkpdsim.uncertainty import UncertaintySettings, adopt_fit

if TYPE_CHECKING:
    from sbmlsim.experiment import SimulationExperiment

    from pkpdsim.model import PKPDModel, SimulatorSettings


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


def parse_workers(text: str, model: PKPDModel) -> int | None:
    """Get the processes of a run from the option `--workers`.

    Args:
        text: `model` for the setting of the model, `auto` for the choice of
            sbmlsim (`None`: a pool of every CPU for a scan of 256 points or
            more), or a positive number.
        model: the model.

    Returns:
        The number of processes, `None` for the choice of sbmlsim.

    Raises:
        ValueError: for any other text.
    """
    if text == "model":
        return model.simulator.n_workers
    if text == "auto":
        return None
    try:
        number = int(text)
    except ValueError:
        number = 0
    if number < 1:
        raise ValueError(
            f"--workers is a positive number, 'auto' or 'model', not '{text}'"
        )
    return number


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
    with_experiments = model.experiments is not None
    if with_experiments:
        actions = [action.value for action in Action]
    else:
        actions = [Action.FACTORY.value, Action.ALL.value]
    actions += [extra.name for extra in model.extra_actions]
    action_help = "the action to perform" + "".join(
        f"; {extra.name}: {extra.help}" for extra in model.extra_actions
    )
    if with_experiments:
        description = f"Run the model factory and the simulation experiments of the {model.title}."
    else:
        extras_text = " and its extra actions" if model.extra_actions else ""
        description = f"Run the model factory of the {model.title}{extras_text}."
    p = argparse.ArgumentParser(prog=f"run_{model.name}", description=description)
    p.add_argument("-a", "--action", required=True, choices=actions, help=action_help)
    p.add_argument(
        "-r",
        "--results-dir",
        type=Path,
        default=model.results_dir,
        help="directory of the results, the simulations are written to its "
        "`simulation` folder (default: the results of the model)",
    )
    if with_experiments:
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
        p.add_argument(
            "-u",
            "--uncertainty",
            action="store_true",
            help="run the experiments over draws of the fitted parameters of the "
            "model and draw the simulations as bands, into the "
            "`simulation_uncertainty` folder of the results; an extra action gets "
            "the settings as `options.uncertainty`",
        )
        p.add_argument(
            "--draws",
            type=int,
            default=UncertaintySettings.draws,
            help="the number of draws of --uncertainty (default: %(default)s)",
        )
        p.add_argument(
            "--seed",
            type=int,
            default=UncertaintySettings.seed,
            help="the seed of the draws of --uncertainty (default: %(default)s)",
        )
        p.add_argument(
            "--quantiles",
            type=float,
            nargs=2,
            default=UncertaintySettings.quantiles,
            metavar=("LOWER", "UPPER"),
            help="the quantiles of the bands of --uncertainty (default: %(default)s)",
        )
        p.add_argument(
            "-w",
            "--workers",
            default="model",
            help="the processes of the scans: a number, 'auto' for a pool of "
            "every CPU for a scan of 256 points or more, or 'model' for the "
            "setting of the model (default: %(default)s); an extra action gets "
            "the settings of the simulator as `options.simulator`",
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
    simulator: SimulatorSettings,
    uncertainty: UncertaintySettings | None = None,
) -> list[str]:
    """Simulate the given groups and experiments.

    Args:
        model: the model.
        names: names of groups or experiments.
        simulation_dir: directory the results of the simulations are written to.
        show_report: open the report in a web browser.
        figure_formats: the formats of the figures.
        simulator: the settings of the simulator.
        uncertainty: run over draws of the fitted parameters, `None` for once.

    Returns:
        The failures of the run, one line each.
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
    try:
        model.run_simulation_experiments(
            classes,
            output_dir,
            show_report=show_report,
            figure_formats=figure_formats,
            simulator=simulator,
            uncertainty=uncertainty,
        )
    except ExperimentRunError as err:
        return list(err.failures)
    return []


def run_command(model: PKPDModel, args: Sequence[str] | None = None) -> None:
    """Run the action of the command line of `run_<model>`.

    Args:
        model: the model.
        args: the command line arguments, `sys.argv[1:]` by default.

    A failure of the simulations does not stop the command: the remaining
    actions run, then the failures are printed and the command exits with 1.
    An extra action gets the settings of the simulator of the run as
    `options.simulator`: those of the model with the processes of `--workers`,
    and the settings of `--uncertainty` as `options.uncertainty` (`None`
    without).
    """
    p = run_parser(model)
    options = p.parse_args(args)
    formats: list[str] = []
    names: list[str] = []
    simulator = model.simulator
    uncertainty: UncertaintySettings | None = None
    if model.experiments is not None:
        if options.uncertainty:
            try:
                uncertainty = UncertaintySettings(
                    draws=options.draws,
                    seed=options.seed,
                    quantiles=tuple(options.quantiles),
                )
            except ValueError as err:
                p.error(str(err))
        formats = parse_formats(options.formats)
        if not formats:
            p.error("'--formats' requires at least one format")
        try:
            n_workers = parse_workers(options.workers, model)
        except ValueError as err:
            p.error(str(err))
        simulator = dataclasses.replace(model.simulator, n_workers=n_workers)
        names = [
            name.strip()
            for name in (options.experiments or "").split(",")
            if name.strip()
        ]
        if options.action == Action.SIMULATE and not names:
            p.error("'--action simulate' requires '--experiments'")

    # an extra action runs its simulations with the settings of the run, which
    # `--workers` changed; an extra option of that name would be replaced
    if hasattr(options, "simulator"):
        msg = (
            f"the extra option '--simulator' of {model.name} clashes with "
            "`options.simulator`, the settings of the simulator of the run"
        )
        raise ValueError(msg)
    if model.experiments is None and hasattr(options, "uncertainty"):
        msg = (
            f"the extra option '--uncertainty' of {model.name} clashes with "
            "`options.uncertainty`, the settings of the draws of the run"
        )
        raise ValueError(msg)
    options.simulator = simulator
    options.uncertainty = uncertainty
    display.section(model.title, icon=":pill:")
    display.key_values({"action": options.action})
    display.link("results", options.results_dir)
    results_dir = options.results_dir.resolve()
    simulation_dir = results_dir / (
        "simulation" if uncertainty is None else "simulation_uncertainty"
    )
    extras = {extra.name: extra for extra in model.extra_actions}
    failures: list[str] = []

    if options.action == Action.FACTORY:
        _run_factory(model)
    elif options.action == Action.LIST_EXPERIMENTS:
        _list_experiments(model)
    elif options.action == Action.SIMULATE:
        failures = _simulate(
            model, names, simulation_dir, options.open, formats, simulator, uncertainty
        )
    elif options.action == Action.ALL:
        _run_factory(model)
        if model.experiments is not None:
            failures = _simulate(
                model,
                ["all"],
                simulation_dir,
                options.open,
                formats,
                simulator,
                uncertainty,
            )
        for extra in model.extra_actions:
            if extra.in_all:
                resolve(extra.function, model.package)(model, options, results_dir)
    else:
        extra = extras[options.action]
        resolve(extra.function, model.package)(model, options, results_dir)

    if failures:
        console.print(
            f"{len(failures)} failures of the simulation experiments, their "
            "tracebacks are in the log above:",
            style="error",
        )
        for line in failures:
            console.print(f"  {line}", style="error", markup=False)
        sys.exit(1)


def with_output_dir(args: Sequence[str], output_dir: Path) -> list[str]:
    """Add the output directory of the model to command line arguments.

    The commands of `sbmlsim.fit.cli` write into the working directory by
    default, a model writes into its results instead.

    Args:
        args: command line arguments without the program name.
        output_dir: output directory used if the arguments do not name one.

    Returns:
        The arguments with the output directory in front of them. argparse keeps
        the last value of an option, so an output directory in the arguments,
        also an abbreviation like `--out`, wins.
    """
    return ["--output_dir", str(output_dir), *args]


def fit_parser(
    model: PKPDModel, definitions: Mapping[str, FitDefinition]
) -> argparse.ArgumentParser:
    """Get the parser of the `fit` command of a model.

    The fit arguments are those of `sbmlsim.fit.cli.fit_cli`, the
    identifiability arguments are the settings of the profile likelihood which
    follows the fit.
    """
    keys = list(definitions)
    if not keys:
        raise ValueError(
            f"The model '{model.name}' has no fit definitions: "
            "at least one FitDefinition is required."
        )
    defaults = ProfileSettings()
    identifiability_text = (
        "The identifiability of the fitted parameters is analysed after the fit, "
        "--no-identifiability skips it."
        if model.fit_identifiability
        else "--identifiability analyses the identifiability of the fitted "
        "parameters after the fit."
    )
    p = argparse.ArgumentParser(
        prog=f"fit_{model.name}",
        description=f"Parameter fitting of the {model.title}. {identifiability_text}",
    )
    fit = p.add_argument_group("fit")
    fit.add_argument(
        "-x", "--subset", choices=keys, default=keys[0], help="fit problem to run"
    )
    fit.add_argument(
        "-n", "--name", default=None, help="name of the report, the fit id by default"
    )
    fit.add_argument(
        "--no-mapping-figures",
        action="store_true",
        help="leave the two figures of every fit mapping out of the report. "
        "They are almost the whole cost of it, so a report which is read for "
        "its tables and its overview figures is created much faster",
    )
    fit.add_argument(
        "-c",
        "--cores",
        type=int,
        default=1,
        help="number of cores for the fitting and the profile scans",
    )
    fit.add_argument(
        "-r", "--runs", type=int, default=4, help="number of optimization runs"
    )
    fit.add_argument(
        "-s", "--seed", type=int, default=1234, help="seed of the optimization"
    )
    fit.add_argument(
        "-m",
        "--method",
        choices=list(ALGORITHMS),
        default="LSQ",
        help="optimization algorithm",
    )
    fit.add_argument(
        "-t",
        "--strategy",
        type=OptimizationStrategy,
        choices=list(OptimizationStrategy),
        default=OptimizationStrategy.ALL,
        help="fit the experiments together or every experiment on its own",
    )
    fit.add_argument(
        "-e",
        "--experiments",
        nargs="+",
        default=None,
        help="experiments to fit, all experiments of the problem by default",
    )
    fit.add_argument(
        "--timeout",
        type=float,
        default=None,
        help="seconds a single optimization may run, no limit by default",
    )
    fit.add_argument(
        "-o",
        "--output_dir",
        type=Path,
        default=model.fit_dir,
        help="directory for the results of the fit",
    )

    adopt = p.add_argument_group("adopt")
    adopt.add_argument(
        "--adopt",
        metavar="FIT_ID",
        default=None,
        help="copy the record of the finished fit FIT_ID from the output "
        "directory into the adopted fits of the model "
        f"({model.fits_dir}) instead of fitting; add its id to `fits` of the "
        "model object",
    )
    adopt.add_argument(
        "--force",
        action="store_true",
        help="replace the record of an adopted fit of the same id",
    )

    identifiability = p.add_argument_group("identifiability")
    identifiability.add_argument(
        "--identifiability",
        action=argparse.BooleanOptionalAction,
        default=model.fit_identifiability,
        help="analyse the identifiability of the fitted parameters after the fit "
        "(default: %(default)s)",
    )
    identifiability.add_argument(
        "-p",
        "--parameter",
        action="append",
        default=None,
        help="parameter to profile, all fitted parameters by default; repeatable",
    )
    identifiability.add_argument(
        "-a",
        "--alpha",
        type=float,
        default=defaults.alpha,
        help="confidence level of the intervals",
    )
    identifiability.add_argument(
        "--df",
        type=int,
        default=defaults.degrees_of_freedom,
        help="degrees of freedom of the threshold, 1 for pointwise intervals",
    )
    identifiability.add_argument(
        "--max-points",
        type=int,
        default=defaults.max_points,
        help="largest number of points of a profile in one direction",
    )
    identifiability.add_argument(
        "--no-reoptimize",
        action="store_true",
        help="keep the other parameters at the optimum, i.e., scan the cost "
        "instead of computing the profile likelihood",
    )
    identifiability.add_argument(
        "--no-fisher",
        action="store_true",
        help="skip the Fisher information, i.e., report the profiles alone",
    )
    return p


def analyse_identifiability(
    run: FitRun, options: argparse.Namespace
) -> tuple[IdentifiabilityResult | None, FisherInformation | None]:
    """Analyse the identifiability of the parameters of a finished fit.

    The profile likelihood is the exact answer and the Fisher information the
    cheap local one next to it; both are computed around the best parameter set
    of the fit and end up in the identifiability section of its report.

    Args:
        run: finished fit whose best parameter set is analysed.
        options: command line arguments of `fit_parser`.

    Returns:
        The profiles and the Fisher information, `None` for what was skipped.
    """
    identifiability = run.identifiability(
        profile_settings=ProfileSettings(
            alpha=options.alpha,
            degrees_of_freedom=options.df,
            max_points=options.max_points,
            reoptimize=not options.no_reoptimize,
        ),
        pids=options.parameter,
        n_cores=options.cores,
    )
    if options.no_fisher:
        return identifiability, None

    fisher = run.fisher(alpha=options.alpha)
    fit_display.section("Fisher information", icon=fit_display.ICON_IDENTIFIABILITY)
    fit_display.key_values(
        {
            "rank": f"{fisher.rank} of {fisher.k}",
            "condition number": f"{fisher.condition_number:.4g}",
        }
    )
    return identifiability, fisher


def _adopt(model: PKPDModel, fit_dir: Path, *, force: bool) -> Path:
    """Adopt a finished fit and print what the model object needs.

    Raises:
        SystemExit: if the fit cannot be adopted.
    """
    fit_display.section("Adopt fit", icon=fit_display.ICON_FIT)
    try:
        record = adopt_fit(fit_dir, model.fits_dir, force=force)
    except (FileNotFoundError, FileExistsError) as err:
        console.print(str(err), style="error")
        sys.exit(1)
    fit_display.key_values({"fit": fit_dir.name, "record": record})
    if fit_dir.name not in model.fits:
        console.print(
            f"Add '{fit_dir.name}' to `fits` of the model object, the factory takes "
            f"the fitted values and their uncertainty from the record "
            f"(`FitRecords.parameter`).",
            style="warning",
        )
    return record


def fit_command(
    model: PKPDModel, args: Sequence[str] | None = None
) -> dict[str, FitRun]:
    """Run a fit of a model and analyse the identifiability of its parameters.

    Args:
        model: the model.
        args: command line arguments, `sys.argv[1:]` by default.

    Returns:
        The finished fits by optimization id, none for `--adopt`.
    """
    definitions = model.load_fit_definitions()
    options = fit_parser(model, definitions).parse_args(args)
    log.enable_rich_logging()
    if options.adopt is not None:
        _adopt(model, options.output_dir / options.adopt, force=options.force)
        return {}

    definition = definitions[options.subset]
    algorithm = ALGORITHMS[options.method]
    # the id of the fit, created before it starts, so that the problem, its
    # result and the directory of its report carry the same key
    opid = fit_id(options.subset)

    fit_display.section("Fit", icon=fit_display.ICON_FIT)
    fit_display.key_values(
        {
            "fit": opid,
            "problem": options.subset,
            "strategy": f"{options.strategy.value} ({STRATEGY_INFO[options.strategy]})",
            "algorithm": f"{options.method} ({algorithm.name})",
            "experiments": options.experiments or "all of the problem",
            "runs": f"{options.runs} on {options.cores} core(s)",
            "seed": options.seed,
            "timeout": (f"{options.timeout} s per run" if options.timeout else "none"),
            "identifiability": (
                "profile likelihood"
                + ("" if options.no_fisher else " and Fisher information")
                if options.identifiability
                else "none"
            ),
            "base path": paths_text(definition.base_path),
            "data path": paths_text(definition.data_path),
            "output": options.output_dir,
        }
    )
    # the parameters are printed per problem by `run_fit`, with their coverage
    fit_display.print_settings(definition.settings)

    runs = run_fit(
        definition=definition,
        opid=opid,
        strategy=options.strategy,
        algorithm=algorithm,
        size=options.runs,
        n_cores=options.cores,
        seed=options.seed,
        study_ids=options.experiments,
        timeout=options.timeout,
        output_dir=options.output_dir,
    )

    # the fit only optimizes, the report is created from its parameters and
    # carries the identifiability of the parameters it reports
    for run in runs.values():
        identifiability, fisher = (
            analyse_identifiability(run, options)
            if options.identifiability
            else (None, None)
        )
        run.report(
            output_dir=options.output_dir,
            name=options.name,
            identifiability=identifiability,
            fisher=fisher,
            show_titles=False,
            mapping_figures=not options.no_mapping_figures,
        )
    return runs


def report_command(model: PKPDModel, args: Sequence[str] | None = None) -> Path:
    """Report stored parameter sets of a model (`sbmlsim.fit.cli.report_cli`).

    Args:
        model: the model.
        args: command line arguments, `sys.argv[1:]` by default.

    Returns:
        The directory of the report.
    """
    return report_cli(
        dict(model.load_fit_definitions()),
        args=with_output_dir(sys.argv[1:] if args is None else args, model.report_dir),
        prog=f"report_{model.name}",
    )


def identifiability_command(
    model: PKPDModel, args: Sequence[str] | None = None
) -> Path:
    """Analyse the identifiability of stored parameters of a model.

    See `sbmlsim.fit.cli.identifiability_cli`.

    Args:
        model: the model.
        args: command line arguments, `sys.argv[1:]` by default.

    Returns:
        The directory of the report.
    """
    return identifiability_cli(
        dict(model.load_fit_definitions()),
        args=with_output_dir(
            sys.argv[1:] if args is None else args, model.identifiability_dir
        ),
        prog=f"identifiability_{model.name}",
    )
