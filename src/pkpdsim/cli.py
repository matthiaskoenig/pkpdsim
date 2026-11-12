"""The commands of the models, wrappers of the run layer of sbmlsim.

pkpdsim installs five console scripts which take the dotted name of a model
package first and the options of the command for that model after it
(`main`, `resolve_model`):

```bash
run_simulation pkdb_models.models.hctz --action factory
run_simulation pkdb_models.models.hctz --action list_experiments
run_simulation pkdb_models.models.hctz --action simulate --experiments "misc,Anderson1961" --open
run_simulation pkdb_models.models.hctz --action all --results-dir /path/to/results
run_simulation pkdb_models.models.hctz --action simulate --experiments all --uncertainty --draws 100
run_fit pkdb_models.models.hctz --subset PK --identifiability
run_fit pkdb_models.models.hctz --adopt PK_<fit id>
run_report pkdb_models.models.hctz --subset PK results/fit/PK_<fit id>/parameters.json
run_report pkdb_models.models.hctz --adopted
run_identifiability pkdb_models.models.hctz --subset PK results/fit/PK_<fit id>/parameters.json
run_credibility pkdb_models.models.hctz
run_credibility pkdb_models.models.hctz --from-json results/credibility/credibility.json -o /tmp/credibility
```

The package is imported and its model object `MODEL` gives everything which
differs between two models; a model package registers no commands of its own.

The commands `run_simulation`, `run_fit`, `run_report` and
`run_identifiability` keep their options and defaults and translate them into
the spec of a kind of the run layer of sbmlsim (`python -m sbmlsim.run`, see
`sbmlsim.run`), which they run with `sbmlsim.run.execute` on the model object
of the package (the project `<package>:MODEL`, `project_reference`). A kind
runs into an empty output directory and writes the records of the run into it
(`run.json`, `spec.json`, `summary.json`, `progress.jsonl`), so a command runs
it into the folder of its result below the directory the command wrote into
before. The run writes into a staging directory next to the folder, and an
earlier result is replaced only after the run succeeded, a failed run or plan
leaves it untouched (`_staged`); an output which is an input of the command
(the records of a fit, the adopted fits, the factory, a build), contains one
or lies inside one is refused before anything happens (`_check_out`):

- `run_simulation --action simulate` (`run_command`) runs the kind `simulate`
  into the folder of the group (`simulation/<group>` of the results), any
  other selection into `simulation/custom_selection`; the figures are
  collected into `_figures` of it (`sbmlsim.project.collect_figure_files`).
  `--uncertainty` runs the experiments over draws of the fitted parameters of
  the model (the field `uncertainty` of the spec) into
  `simulation_uncertainty/<group>`.
- `run_fit` (`fit_command`) runs the kind `fit` into `<output_dir>/<fit id>`,
  with the report of the fit, the identifiability of its parameters with
  `--identifiability` (the default of a model with `fit_identifiability`) and
  the adopted fits of the other subsets recorded (`record_adopted`); the
  records of the fit, which the kind writes into its folder `<fit id>`, are
  moved into the output of the run, so they are where they were.
  A fit which fails keeps the records it wrote in `<output_dir>/<fit id>`.
  `run_fit <package> --adopt <fit id>` checks the id, runs the kind `adopt`
  into a temporary output and copies the record into the adopted fits of the
  model (`sbmlsim.fit.records.adopt_fit`); nothing is written into the fit.
- `run_report` (`report_command`) runs the kind `report` on the records next
  to a `parameters.json` into `<output_dir>/<name>` (`report` by default);
  `run_report <package> --adopted` runs it with `adopted` for every adopted fit
  into `<output_dir>/<fit id>`.
- `run_identifiability` (`identifiability_command`) runs the kind
  `identifiability` into `analysis/` of `<output_dir>/<name>`
  (`identifiability` by default) and the kind `report` of the analysis into
  it, the report in the folder itself and the records of its run in `report/`.

Three parts call the library of sbmlsim directly instead: `--action factory`
builds the models of the package in place with `create_models`, as before;
`--action factory --rebuild` runs the kind `factory`, which builds them into
`factory/` of the results and reports the outputs which differ. The extra
actions of a model write into the results of the model, which are no empty
output; they are the `actions` of the project as well
(`pkpdsim.model.ExtraActionFunction`), which the kind `action` runs. `run_fit
--strategy single` fits every experiment on its own, which no kind does, and
`run_report` and `run_identifiability` of several parameter files or of a file
without the records of its fit next to it use `sbmlsim.fit.cli.report_cli` and
`identifiability_cli`.

`run_credibility` (`credibility_command`) runs the kind `credibility`, which
assesses the contexts of use of the model against its adopted fits
(`sbmlsim.credibility.assess`, with the evidence extensions and the PK battery
of `pkpdsim.credibility`) and writes the record, the tables and the report
(`sbmlsim.credibility.report`) into `credibility/` of its output; the command
moves them up into its `--output_dir`, `--from-json` is the input `record`
(the outputs written again from a written record) and the lock of `--lock`,
which the kind writes into its output, is moved into the folder of the model,
where the assessment compares the evidence with it. All write into the
results of the model unless `--output_dir` is given.

`--profile` anywhere after the model package runs a command under the
sampling profiler of python 3.15 and writes its flame graph and statistics
into `profiling/<command>_<time>/` of the results of the model, see
`sbmlsim.run.profiling`; the commands also run as
`python -m pkpdsim.cli <command> <package> ...`.

`--progress MODE` anywhere after the model package says how the progress of
the runs of a command is shown: `auto` (the default: the bars of a terminal,
a line of the log every 30 s otherwise), `bar`, `log` or `none`, see
`sbmlsim.progress`. The option sets `SBMLSIM_PROGRESS`, so every run of the
command and every process it starts reads it, and the events of the progress
are written into the `progress.jsonl` of a run whatever is chosen: the
experiments of a simulation, the points of a scan, the runs of a fit, the
scans of an identifiability and the figures of a report then say where they
are, how long they still take and what they use of the machine.
"""

from __future__ import annotations

import argparse
import dataclasses
import datetime
import importlib
import json
import logging
import math
import os
import re
import shutil
import sys
import tempfile
import uuid
import webbrowser
from collections.abc import Callable, Mapping, Sequence
from enum import StrEnum
from pathlib import Path
from typing import TYPE_CHECKING, Any, NoReturn

from sbmlsim import display, log
from sbmlsim.atomic import write_text_atomic
from sbmlsim.cgroup import available_cpus
from sbmlsim.experiment import ExperimentRunError
from sbmlsim.experiment.draws import UncertaintySettings
from sbmlsim.experiment.runner import EXPERIMENTS_PER_WORKER
from sbmlsim.fit import display as fit_display
from sbmlsim.fit.cli import (
    ALGORITHMS,
    STRATEGY_INFO,
    FitDefinition,
    FitRun,
    identifiability_cli,
    records_context,
    report_cli,
    report_context,
    run_fit,
)
from sbmlsim.fit.fisher import FisherInformation
from sbmlsim.fit.identifiability import IdentifiabilityResult, ProfileSettings
from sbmlsim.fit.options import OptimizationAlgorithmType, OptimizationStrategy
from sbmlsim.fit.project_fits import (
    FIXED_FITS,
    fixed_values,
    remainder_problems,
    start_point,
    with_start,
)
from sbmlsim.fit.records import adopt_fit
from sbmlsim.fit.result import FIT_ID, fit_id, fit_subset
from sbmlsim.fit.sampling import SamplingType
from sbmlsim.progress import ENV_NAME as PROGRESS_ENV
from sbmlsim.progress import NAMES as PROGRESS_NAMES
from sbmlsim.project import Project, collect_figure_files, load_project, resolve
from sbmlsim.run import ConfigurationError, execute, list_analyses, plan
from sbmlsim.run.context import SPEC_FILE
from sbmlsim.run.profiling import (
    PROFILE_DIR,
    PROFILE_OPTION,
    profiler_available,
    record_exit_code,
    run_profiled,
)
from sbmlsim.utils import paths_text
from sbmlutils.console import console

from pkpdsim.model import RUN_WORKERS, PKPDModel

if TYPE_CHECKING:
    from sbmlsim.experiment import SimulationExperiment


logger = logging.getLogger(__name__)


class Action(StrEnum):
    """The actions of the `run` command, `pkpdsim.model.RUN_ACTIONS`."""

    SIMULATE = "simulate"
    LIST_EXPERIMENTS = "list_experiments"
    FACTORY = "factory"
    ALL = "all"


#: the option of the progress of the runs of a command, which `main` takes out
#: of the arguments and sets as `SBMLSIM_PROGRESS`, as `--profile` of
#: `sbmlsim.run.profiling` is taken out of them
PROGRESS_OPTION = "--progress"

#: the attribute of the model object in its package
MODEL_ATTRIBUTE = "MODEL"

#: the folder of the results of `--action factory --rebuild`
FACTORY_DIR = "factory"

#: the folder of the run of the kind `identifiability` in the folder of an
#: identifiability analysis, next to its report
ANALYSIS_DIR = "analysis"

#: the folder of the run of the kind `report` in the folder of an
#: identifiability analysis, whose report is moved into the folder
REPORT_RUN_DIR = "report"

#: the folder of the report in the output of the kind `report`
#: (`sbmlsim.run.fit_kinds.REPORT_NAME`)
REPORT_NAME = "report"

#: the folder of the adopted fit in the output of the kind `adopt`
#: (`sbmlsim.run.fit_kinds.FITS_DIR`)
FITS_DIR = "fits"

#: the file of the best parameters of a fit, which the records of a fit have
PARAMETERS_FILE = "parameters.json"


#: the folder of the outputs of the kind `credibility` in its output
CREDIBILITY_DIR = "credibility"

#: the record of an assessment, which marks a directory of its outputs
CREDIBILITY_RECORD = "credibility.json"


def parse_formats(text: str) -> list[str]:
    """Get the figure formats of a comma separated list.

    Args:
        text: the formats, e.g. `"svg, PNG"`.

    Returns:
        The formats in lower case, without blanks and empty entries.
    """
    return [part.strip().lower() for part in text.split(",") if part.strip()]


def parse_workers(text: str) -> int | None:
    """Get the processes of the option `--workers` which an extra action gets.

    Args:
        text: `model` for the default of a model (`pkpdsim.model.RUN_WORKERS`),
            `auto` for the choice of sbmlsim (`None`, see
            `sbmlsim.project.Project.run_experiments`), or a positive number.

    Returns:
        The number of processes, `None` for the choice of sbmlsim.

    Raises:
        ValueError: for any other text.
    """
    if text == "model":
        return RUN_WORKERS
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


def run_workers(text: str, n_experiments: int) -> int | None:
    """Get the workers of a run of experiments from `--workers`.

    Args:
        text: the option, see `parse_workers`.
        n_experiments: the number of experiments of the run.

    Returns:
        The processes of the kind `simulate`: for `model`, the default of a
        model, a worker per four experiments, at most one per CPU, serially
        below (the rule of `sbmlsim.experiment.ExperimentRunner`), every
        experiment with its scans serially in its worker; the number as given;
        `None` for `auto`, the choice of sbmlsim, which the kind cannot take
        (it resolves its workers), so the experiments run with
        `Project.run_experiments` itself.

    Raises:
        ValueError: for an invalid option, see `parse_workers`.
    """
    number = parse_workers(text)
    if text == "model":
        return max(
            1, min(available_cpus(), math.ceil(n_experiments / EXPERIMENTS_PER_WORKER))
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


# ----------------------------------------------------------------------------
# the run layer
# ----------------------------------------------------------------------------
def project_reference(model: PKPDModel) -> str:
    """Get the reference of a model as the project of the run layer.

    The kinds of `sbmlsim.run` load their project by its reference, so a
    command runs the model object of the package, `<package>:MODEL`.

    Args:
        model: the model.

    Returns:
        The reference `<package>:MODEL`.

    Raises:
        ValueError: if the reference does not give the model, e.g. for a copy
            of the model object with other fields.
    """
    reference = f"{model.package}:{MODEL_ATTRIBUTE}"
    try:
        project = load_project(reference)
    except (ImportError, ValueError) as err:
        raise ValueError(
            f"The model '{model.name}' is not loaded by its reference "
            f"'{reference}': {err}"
        ) from err
    if project is not model and project != model:
        raise ValueError(
            f"The commands run the model object of the package, '{reference}', "
            f"which the run layer of sbmlsim loads by its reference; the model "
            f"'{model.name}' is another one, e.g. a copy with other fields."
        )
    return reference


def _same(a: Path, b: Path) -> bool:
    """Check whether two resolved paths are one, also on a case-insensitive file system."""
    if a == b:
        return True
    try:
        return os.path.samefile(a, b)
    except OSError:
        return False


def _overlaps(out: Path, path: Path) -> bool:
    """Check whether an output is a path, contains it or lies inside it."""
    target, other = out.resolve(), path.resolve()
    return (
        _same(target, other)
        or any(_same(target, parent) for parent in other.parents)
        or any(_same(other, parent) for parent in target.parents)
    )


def _check_out(
    parser: argparse.ArgumentParser,
    model: PKPDModel,
    out: Path,
    inputs: Sequence[Path] = (),
) -> None:
    """Refuse an output which would replace an input of a command.

    A command replaces its output, so the output must not be an input of the
    command (the records of a fit), the adopted fits, the factory or a build
    of the model, contain one of them or lie inside one, nor contain the
    package or a folder of its data.

    Raises:
        SystemExit: a usage error naming the input.
    """
    for path in [*inputs, model.fits_dir, model.factory_dir, *model.variants.values()]:
        if _overlaps(out, path):
            parser.error(
                f"The output '{out}' is or contains the input '{path}' of the "
                "command, or lies inside it, which the output would replace: "
                "give another output."
            )
    for path in [model.root, *model.data_paths]:
        target, other = out.resolve(), path.resolve()
        if _same(target, other) or any(_same(target, p) for p in other.parents):
            parser.error(
                f"The output '{out}' contains the input '{path}' of the command, "
                "which the output would replace: give another output."
            )


def _staged[T](
    target: Path, run: Callable[[Path], T], *, inner: str | None = None
) -> T:
    """Run into a staging directory and replace the result only when the run succeeded.

    The run writes into `.<name>.<uuid>` next to the target, a kind of the
    run layer into an empty output; the result folder `inner` of the kind
    (e.g. the records of a fit in `<fit id>`) is moved up into the staging
    directory, then an earlier result is renamed aside, the staging directory
    renamed into its place and the earlier result removed (the pattern of
    `sbmlsim.fit.records.adopt_fit`). A failed run, also a failed plan,
    removes the staging directory and leaves the earlier result untouched.

    Args:
        target: the folder of the result.
        run: runs into the staging directory, which it creates.
        inner: the folder of the result in the staging directory, relative.

    Returns:
        What the run returns.

    Raises:
        SystemExit: exit 1 if the result is not moved into its place.
    """
    token = uuid.uuid4().hex
    staging = target.parent / f".{target.name}.{token}"
    target.parent.mkdir(parents=True, exist_ok=True)
    _sweep(target)
    try:
        result = run(staging)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    try:
        if inner is not None:
            _moved_up(staging / inner, staging)
    except OSError as err:
        shutil.rmtree(staging, ignore_errors=True)
        _failed(f"The result '{target}' is not written: {err}")
    earlier = target.parent / f".{target.name}.old.{token}"
    if target.exists():
        console.print(f"The earlier results in '{target}' are replaced.", markup=False)
        try:
            target.rename(earlier)
        except OSError as err:
            shutil.rmtree(staging, ignore_errors=True)
            _failed(f"The result '{target}' is not written: {err}")
    try:
        staging.rename(target)
    except OSError as err:
        shutil.rmtree(staging, ignore_errors=True)
        if earlier.exists():
            try:
                earlier.rename(target)
            except OSError as restore:
                _failed(
                    f"The result '{target}' is not written ({err}) and the earlier "
                    f"result is not restored ({restore}), it is in '{earlier}'."
                )
        _failed(f"The result '{target}' is not written: {err}")
    shutil.rmtree(earlier, ignore_errors=True)
    return result


def _sweep(target: Path) -> None:
    """Remove the stagings and earlier results of a target which a killed run left.

    A run which was killed leaves `.<name>.<uuid>` or `.<name>.old.<uuid>`
    next to the target, see `_staged`; each removal is logged.
    """
    pattern = re.compile(rf"^\.{re.escape(target.name)}\.(old\.)?[0-9a-f]{{32}}$")
    for entry in sorted(target.parent.iterdir()):
        if entry.is_dir() and pattern.match(entry.name):
            logger.warning("A run which was killed left '%s', it is removed.", entry)
            shutil.rmtree(entry, ignore_errors=True)


def _relocated_input(spec_file: Path, staged: Path, final: Path) -> None:
    """Name an input of a run where it is after the staged swap, see `_staged`.

    The spec of a run records the absolute paths of its inputs; an input in
    the staging directory is renamed with it.

    Args:
        spec_file: the planned spec of the run, `spec.json`.
        staged: the input in the staging directory.
        final: the input after the swap.
    """
    spec = json.loads(spec_file.read_text(encoding="utf-8"))
    old = str(staged.resolve())
    changed = {
        key: str(final.resolve())
        for key, value in spec.items()
        if isinstance(value, str) and value == old
    }
    if changed:
        spec.update(changed)
        write_text_atomic(spec_file, json.dumps(spec, indent=2) + "\n")


def _plain_name(parser: argparse.ArgumentParser, name: str, option: str) -> str:
    """Check that a name of a folder of the results is one plain name.

    Raises:
        SystemExit: a usage error for a name with a separator, `.` or `..`.
    """
    if not name or Path(name).name != name or name in (".", ".."):
        parser.error(f"{option} is the name of a folder, not '{name}'")
    return name


def _moved_up(source: Path, target: Path) -> None:
    """Move the entries of a directory into another one and remove it.

    The output of a kind has its result in a folder of its own, e.g. the
    records of a fit in `<fit id>` of the output; a command moves them up into
    the folder of its result, where they were before.

    Args:
        source: the folder of the result.
        target: the folder the entries are moved into.

    Raises:
        FileExistsError: if an entry exists in the target.
    """
    entries = sorted(source.iterdir())
    clashes = [entry.name for entry in entries if (target / entry.name).exists()]
    if clashes:
        raise FileExistsError(
            f"The results of '{source}' are not moved into '{target}', which "
            f"has {clashes}."
        )
    for entry in entries:
        entry.rename(target / entry.name)
    source.rmdir()


def _execute(
    parser: argparse.ArgumentParser,
    kind: str,
    spec: dict[str, Any],
    out: Path,
    *,
    workers: int | None = None,
    inputs: Mapping[str, Path] | None = None,
) -> dict[str, Any]:
    """Run a kind of the run layer, see `sbmlsim.run.execute`.

    The spec is planned first: a spec which does not plan, e.g. an unknown
    experiment, an invalid start of a fit or a remainder of the model which
    can become negative, is a usage error (exit 2); a configuration error of
    the run itself, e.g. a fit adopted already, exits with 1.

    Returns:
        The summary of the run.

    Raises:
        SystemExit: for a spec which does not plan or a run which is
            configured wrongly (`ConfigurationError`).
    """
    given = {
        **spec,
        **{name: str(Path(path).resolve()) for name, path in (inputs or {}).items()},
    }
    try:
        plan(kind, given)
    except ConfigurationError as err:
        parser.error(str(err))
    try:
        return execute(kind, spec, out, workers=workers, inputs=inputs)
    except ConfigurationError as err:
        _failed(str(err))


def _failed(message: str) -> NoReturn:
    """Print an error of a command and exit with 1."""
    console.print(message, style="error")
    sys.exit(1)


# ----------------------------------------------------------------------------
# run_simulation
# ----------------------------------------------------------------------------
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
    p = argparse.ArgumentParser(
        prog=f"run_simulation {model.package}", description=description
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
        "--rebuild",
        action="store_true",
        help="with the factory (--action factory or all): build the models "
        "into the `factory` folder of the results instead of the package and "
        "report the outputs of the factory which differ from the ones of the "
        "package (the kind `factory` of `python -m sbmlsim.run`); the command "
        "exits with 1 if one differs",
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
            help="the processes of the run: a number (the experiments in a pool, "
            "every experiment serially in its worker, a single experiment with "
            "its scans in the pool), 'model' for the default of a model (a "
            "worker per four experiments, at most one per CPU, serially below, "
            "the scans of every experiment serially in its worker) or 'auto' "
            "for the choice of sbmlsim (a worker per four experiments, the "
            "scans in a pool of every CPU when one has 256 points or more) "
            "(default: %(default)s); an extra action gets the settings of the "
            "simulator as `options.simulator` and the processes as "
            "`options.n_workers`: the number, 1 for 'model' and None for 'auto'",
        )
    if model.extra_arguments is not None:
        resolve(model.extra_arguments, model.package)(p)
    return p


def _run_factory(
    model: PKPDModel,
    parser: argparse.ArgumentParser,
    results_dir: Path,
    *,
    rebuild: bool,
) -> list[str]:
    """Create the SBML models, the ODEs and the COMBINE archive.

    Without `rebuild` the factory builds the models of the package in place,
    with it the kind `factory` builds them into the results and compares them
    with the ones of the package.

    Returns:
        The outputs of a rebuild which differ from the ones of the package.
    """
    display.section("Model factory", icon=":factory:")
    if not rebuild:
        model.create_models()
        display.link("models", model.factory_dir)
        if model.archive_path is not None:
            display.link("archive", model.archive_path)
        return []
    out = results_dir / FACTORY_DIR
    _check_out(parser, model, out)
    spec = {"project": project_reference(model)}
    summary = _staged(out, lambda staging: _execute(parser, "factory", spec, staging))
    changed: list[str] = list(summary["changed"])
    display.key_values({"compared": summary["outputs"], "changed": len(changed)})
    display.link("models", out / summary["directory"])
    return changed


def _list_experiments(model: PKPDModel) -> None:
    """Print the groups and the experiments which can be simulated."""
    entries = list_analyses(project_reference(model), experiments=True)
    display.section("Simulation experiments", icon=":test_tube:")
    display.key_values(
        {
            entry["name"]: entry["description"]
            for entry in entries
            if entry["kind"] == "simulate"
        }
    )


def _simulate(
    model: PKPDModel,
    parser: argparse.ArgumentParser,
    names: Sequence[str],
    simulation_dir: Path,
    *,
    show_report: bool,
    figure_formats: Sequence[str],
    workers: str,
    uncertainty: UncertaintySettings | None = None,
) -> list[str]:
    """Simulate the given groups and experiments with the kind `simulate`.

    A name which is neither a group nor an experiment is reported and left
    out, the run of no experiment exits with 1. The run replaces an earlier
    result only when it ran (`_staged`), also with failed experiments, whose
    report is written. `--workers auto` runs the experiments with
    `Project.run_experiments` and the choice of sbmlsim, without the records
    of a run, see `run_workers`.

    Args:
        model: the model.
        parser: the parser of the command, for the errors of the spec.
        names: names of groups or experiments.
        simulation_dir: directory the results of the simulations are written to.
        show_report: open the report in a web browser.
        figure_formats: the formats of the figures.
        workers: the option `--workers`, see `run_workers`.
        uncertainty: run over draws of the fitted parameters, `None` for once.

    Returns:
        The failures of the run, one line each.
    """
    reference = project_reference(model)
    registry = model.load_experiments()
    classes, not_found = resolve_experiments(registry, names)
    if not_found:
        console.print(f"Not found: {', '.join(not_found)}", style="warning")
    if not classes:
        _failed("No experiments to simulate.")

    single_group = len(names) == 1 and names[0] in registry
    output_dir = simulation_dir / (names[0] if single_group else "custom_selection")
    display.section("Simulations", icon=":chart_increasing:")
    display.key_values({"experiments": len(classes)})
    display.link("output", output_dir)
    _check_out(parser, model, output_dir)
    n_workers = run_workers(workers, len(classes))
    if n_workers is not None and (model.save_results or not model.reduced_selections):
        console.print(
            "The kind `simulate` runs the experiments with the selections they "
            "use and without their results; `save_results` and "
            "`reduced_selections` of the model apply to `run_experiments`.",
            style="warning",
        )
    spec: dict[str, Any] = {
        "project": reference,
        "experiments": [name for name in names if name not in not_found],
        "figure_formats": list(figure_formats),
    }
    if uncertainty is not None:
        spec["uncertainty"] = {
            "draws": uncertainty.draws,
            "seed": uncertainty.seed,
            "quantiles": list(uncertainty.quantiles),
        }

    def run(staging: Path) -> list[str]:
        try:
            if n_workers is None:
                Project.run_experiments(
                    model,
                    classes,
                    staging,
                    n_workers=None,
                    figure_formats=figure_formats,
                    save_results=model.save_results,
                    reduced_selections=model.reduced_selections,
                    uncertainty=uncertainty,
                )
            else:
                _execute(parser, "simulate", spec, staging, workers=n_workers)
        except ExperimentRunError as err:
            collect_figure_files(staging)
            return list(err.failures)
        collect_figure_files(staging)
        return []

    failures = _staged(output_dir, run)
    display.link("results", output_dir / "index.html")
    if show_report:
        webbrowser.open((output_dir / "index.html").resolve().as_uri())
    return failures


def run_command(model: PKPDModel, args: Sequence[str] | None = None) -> None:
    """Run the action of the command line of `run_simulation`.

    Args:
        model: the model, the model object of its package.
        args: the command line arguments, `sys.argv[1:]` by default.

    A failure of the simulations or a rebuild which differs from the models of
    the package does not stop the command: the remaining actions run, then the
    failures are printed and the command exits with 1. An extra action gets
    the settings of the simulator of the model as `options.simulator`, the
    processes of `--workers` as `options.n_workers` (`parse_workers`), and the
    settings of `--uncertainty` as `options.uncertainty` (`None` without).
    """
    p = run_parser(model)
    options = p.parse_args(args)
    formats: list[str] = []
    names: list[str] = []
    simulator = model.simulator
    n_workers: int | None = RUN_WORKERS
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
            n_workers = parse_workers(options.workers)
        except ValueError as err:
            p.error(str(err))
        names = [
            name.strip()
            for name in (options.experiments or "").split(",")
            if name.strip()
        ]
        if options.action == Action.SIMULATE and not names:
            p.error("'--action simulate' requires '--experiments'")

    # an extra action runs its simulations with the settings of the run and the
    # processes of `--workers`; an extra option of these names would be replaced
    if hasattr(options, "simulator"):
        msg = (
            f"the extra option '--simulator' of {model.name} clashes with "
            "`options.simulator`, the settings of the simulator of the run"
        )
        raise ValueError(msg)
    if hasattr(options, "n_workers"):
        msg = (
            f"the extra option '--n-workers' of {model.name} clashes with "
            "`options.n_workers`, the processes of the run"
        )
        raise ValueError(msg)
    if model.experiments is None and hasattr(options, "uncertainty"):
        msg = (
            f"the extra option '--uncertainty' of {model.name} clashes with "
            "`options.uncertainty`, the settings of the draws of the run"
        )
        raise ValueError(msg)
    options.simulator = simulator
    options.n_workers = n_workers
    options.uncertainty = uncertainty
    log.enable_rich_logging()
    display.section(model.title, icon=":pill:")
    display.key_values({"action": options.action})
    display.link("results", options.results_dir)
    results_dir = options.results_dir.resolve()
    simulation_dir = results_dir / (
        "simulation" if uncertainty is None else "simulation_uncertainty"
    )
    extras = {extra.name: extra for extra in model.extra_actions}
    failures: list[str] = []
    changed: list[str] = []

    def simulate(selection: Sequence[str]) -> list[str]:
        return _simulate(
            model,
            p,
            selection,
            simulation_dir,
            show_report=options.open,
            figure_formats=formats,
            workers=options.workers,
            uncertainty=uncertainty,
        )

    if options.action == Action.FACTORY:
        changed = _run_factory(model, p, results_dir, rebuild=options.rebuild)
    elif options.action == Action.LIST_EXPERIMENTS:
        _list_experiments(model)
    elif options.action == Action.SIMULATE:
        failures = simulate(names)
    elif options.action == Action.ALL:
        changed = _run_factory(model, p, results_dir, rebuild=options.rebuild)
        if model.experiments is not None:
            failures = simulate(["all"])
        for extra in model.extra_actions:
            if extra.in_all:
                resolve(extra.function, model.package)(model, options, results_dir)
    else:
        extra = extras[options.action]
        resolve(extra.function, model.package)(model, options, results_dir)

    if changed:
        console.print(
            f"{len(changed)} outputs of the factory differ from the ones of the "
            "package:",
            style="error",
        )
        for path in changed:
            console.print(f"  {path}", style="error", markup=False)
    if failures:
        console.print(
            f"{len(failures)} failures of the simulation experiments, their "
            "tracebacks are in the log above:",
            style="error",
        )
        for line in failures:
            console.print(f"  {line}", style="error", markup=False)
    if changed or failures:
        sys.exit(1)


# ----------------------------------------------------------------------------
# run_fit
# ----------------------------------------------------------------------------
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
        prog=f"run_fit {model.package}",
        description=f"Parameter fitting of the {model.title}. {identifiability_text}",
    )
    fit = p.add_argument_group("fit")
    fit.add_argument(
        "-x", "--subset", choices=keys, default=keys[0], help="fit problem to run"
    )
    fit.add_argument(
        "-n",
        "--name",
        default=None,
        help="name of the reports of --strategy single, the fit id by default; "
        "ignored otherwise, the report of a fit is in the folder of its id",
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
        "--start",
        metavar="JSON_OR_FIT_ID",
        default=None,
        help="start least squares from these parameters instead of samples of "
        "the bounds: a json file of the values by parameter id in the units of "
        "the fit parameters, or the id of a fit (adopted, or finished in the "
        "output directory) whose best parameters are taken; every run starts "
        "from the point, so the fit has one run (--runs 1); a parameter the "
        "start does not give starts from its start value",
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
    adopt.add_argument(
        "--cross-validation",
        metavar="FIT_ID",
        default=None,
        help="cross-validate the adopted fit FIT_ID instead of fitting: refit "
        "its problem once per training study without the study, from the "
        "adopted parameters, and store the predictions of the held-out studies "
        "in its record (cross_validation.json); --experiments limits the "
        "studies, --cores and --seed are those of the refits",
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


def _adopt(
    model: PKPDModel,
    parser: argparse.ArgumentParser,
    output_dir: Path,
    fit_id: str,
    *,
    force: bool,
) -> Path:
    """Adopt a finished fit and print what the model object needs.

    The argument is checked first: the id of a fit, whose records are in the
    output directory. The kind `adopt` checks the fit and writes its record
    into a temporary output, which is copied into the adopted fits of the
    model (`adopt_fit`); nothing is written into the fit.

    Raises:
        SystemExit: a usage error for an argument which is no id of a fit,
            exit 1 if the fit is not found or cannot be adopted.
    """
    fit_display.section("Adopt fit", icon=fit_display.ICON_FIT)
    if Path(fit_id).name != fit_id or FIT_ID.match(fit_id) is None:
        parser.error(
            f"--adopt is the id of a fit, <subset>_<date>_<time>__<hash>, "
            f"not '{fit_id}'"
        )
    fit_dir = output_dir / fit_id
    if not (fit_dir / PARAMETERS_FILE).is_file():
        _failed(f"The fit '{fit_dir}' does not exist, it has no {PARAMETERS_FILE}.")
    reference = project_reference(model)
    with tempfile.TemporaryDirectory(prefix="pkpdsim-adopt-") as tmp:
        out = Path(tmp) / "adopt"
        summary = _execute(
            parser,
            "adopt",
            {"project": reference, "force": force},
            out,
            inputs={"fit": fit_dir},
        )
        try:
            record = adopt_fit(
                out / FITS_DIR / summary["opid"], model.fits_dir, force=force
            )
        except (FileNotFoundError, FileExistsError) as err:
            _failed(str(err))
    fit_display.key_values({"fit": record.name, "record": record})
    if record.name not in model.adopted_ids:
        console.print(
            f"Add '{record.name}' to `fits` of the model object, the factory takes "
            f"the fitted values and their uncertainty from the record "
            f"(`AdoptedRecords.parameter`).",
            style="warning",
        )
    return record


def _start(model: PKPDModel, start: str | None, output_dir: Path) -> str | None:
    """Get the start of the spec of a fit from `--start`.

    The id of a fit is an adopted fit of the model, which the plan of the fit
    finds by its id, else a finished fit in the output directory, whose
    records are given by their path; any other start is a path.

    Returns:
        The start of the spec, `None` without.
    """
    if start is None:
        return None
    path = Path(start)
    if len(path.parts) == 1 and FIT_ID.match(start) is not None:
        if (model.fits_dir / start / PARAMETERS_FILE).is_file():
            return start
        if (output_dir / start / PARAMETERS_FILE).is_file():
            return str((output_dir / start).resolve())
    return str(path.resolve()) if path.exists() else start


def _show_fit(
    options: argparse.Namespace, opid: str, definition: FitDefinition, output: Path
) -> None:
    """Show the settings of a fit before it runs."""
    algorithm = ALGORITHMS[options.method]
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
            "start": options.start or "samples of the bounds",
            "identifiability": (
                "profile likelihood"
                + ("" if options.no_fisher else " and Fisher information")
                if options.identifiability
                else "none"
            ),
            "base path": paths_text(definition.base_path),
            "data path": paths_text(definition.data_path),
            "output": output,
        }
    )
    # the parameters are printed per problem by `run_fit`, with their coverage
    fit_display.print_settings(definition.settings)


def _fit(
    model: PKPDModel,
    parser: argparse.ArgumentParser,
    options: argparse.Namespace,
    definition: FitDefinition,
) -> dict[str, Path]:
    """Fit the experiments of a problem together with the kind `fit`.

    Returns:
        The records of the fit by its id.
    """
    reference = project_reference(model)
    if options.name is not None:
        console.print(
            "--name is ignored: the report of a fit is in the folder of its id.",
            style="warning",
        )
    opid = fit_id(options.subset)
    out = options.output_dir / opid
    _check_out(parser, model, out)
    _show_fit(options, opid, definition, out)
    spec: dict[str, Any] = {
        "project": reference,
        "subset": options.subset,
        "opid": opid,
        "repeats": options.runs,
        "algorithm": options.method,
        "seed": options.seed,
        "studies": options.experiments,
        "timeout": options.timeout,
        "start": _start(model, options.start, options.output_dir),
        "record_adopted": True,
        "identifiability": options.identifiability,
        "profile": {
            "alpha": options.alpha,
            "degrees_of_freedom": options.df,
            "max_points": options.max_points,
            "reoptimize": not options.no_reoptimize,
            "pids": options.parameter,
        },
        "fisher": not options.no_fisher,
        "report": True,
        "mapping_figures": not options.no_mapping_figures,
    }
    try:
        _execute(parser, "fit", spec, out, workers=options.cores)
    finally:
        # the records of the fit where they were, in the folder of the fit id,
        # also of a fit which failed, so `--start <fit id>` finds them
        if (out / opid).is_dir():
            try:
                _moved_up(out / opid, out)
            except OSError as err:
                console.print(
                    f"The records of the fit stay in '{out / opid}': {err}",
                    style="error",
                    markup=False,
                )
    # the report of the kind named its folder before it was moved
    display.link("results", out / "index.html")
    return {opid: out}


def _fit_each(
    model: PKPDModel,
    parser: argparse.ArgumentParser,
    options: argparse.Namespace,
    definition: FitDefinition,
) -> dict[str, Path]:
    """Fit every experiment of a problem on its own, `--strategy single`.

    No kind fits the experiments one at a time, so the fit runs with
    `sbmlsim.fit.cli.run_fit`, as the command did before.

    Returns:
        The records of the fits by their ids.
    """
    algorithm = ALGORITHMS[options.method]
    # the start of least squares: every run starts from the point (the sampling
    # START of sbmlsim), the global optimizer draws its population itself
    # keyword arguments of `run_fit`, whose other keywords have other types
    sampling: dict[str, Any] = {}
    if options.start is not None:
        if algorithm != OptimizationAlgorithmType.LEAST_SQUARE:
            parser.error(
                "--start is the start of least squares; differential evolution "
                "draws its population from the bounds."
            )
        if options.runs != 1:
            parser.error(
                "--start starts every run from the point, the runs would be "
                "identical: give --runs 1."
            )
        try:
            point = start_point(model, definition, options.start, options.output_dir)
        except ValueError as err:
            parser.error(str(err))
        definition = with_start(definition, point)
        sampling["sampling"] = SamplingType.START
    # the id of the fit, created before it starts, so that the problem, its
    # result and the directory of its report carry the same key
    opid = fit_id(options.subset)

    _show_fit(options, opid, definition, options.output_dir)

    # the adopted fits of the other subsets and the values of their parameters,
    # in the units of the model, which this fit does not optimize itself
    fixed_fits, fixed = fixed_values(model, options.subset, definition)
    # every build which the experiments of the fit run has the remainders
    remainders = list(
        dict.fromkeys(
            problem
            for path in model.model_paths.values()
            for problem in remainder_problems(path, definition.parameters, fixed)
        )
    )
    if remainders:
        parser.error("; ".join(remainders))
    context = report_context(model, options.subset, definition, fixed_fits=fixed_fits)

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
        arguments={
            "subset": options.subset,
            "runs": options.runs,
            "cores": options.cores,
            "method": options.method,
            "command": f"run_fit {model.package}",
            "start": options.start,
            FIXED_FITS: fixed_fits,
        },
        fixed_values=fixed,
        **sampling,
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
            context=context,
        )
    return {key: options.output_dir / key for key in runs}


def fit_command(model: PKPDModel, args: Sequence[str] | None = None) -> dict[str, Path]:
    """Run a fit of a model and analyse the identifiability of its parameters.

    The experiments of the problem are fitted together with the kind `fit`
    of the run layer into `<output_dir>/<fit id>`, every experiment on its own
    (`--strategy single`) with `sbmlsim.fit.cli.run_fit`; `--adopt` adopts a
    finished fit with the kind `adopt`.

    Args:
        model: the model, the model object of its package.
        args: command line arguments, `sys.argv[1:]` by default.

    Returns:
        The records of the finished fits by their ids, none for `--adopt` and
        `--cross-validation`.
    """
    definitions = model.load_fit_definitions()
    parser = fit_parser(model, definitions)
    options = parser.parse_args(args)
    log.enable_rich_logging()
    if options.adopt is not None:
        _adopt(model, parser, options.output_dir, options.adopt, force=options.force)
        return {}
    if options.cross_validation is not None:
        _cross_validate(model, options)
        return {}
    # the problems of the fit load the models at the level of the model, in the
    # process and in every worker of a pool (`PKPDModel.load_fit_definitions`)
    definition = definitions[options.subset]
    if options.strategy == OptimizationStrategy.ALL:
        return _fit(model, parser, options, definition)
    return _fit_each(model, parser, options, definition)


def _cross_validate(model: PKPDModel, options: argparse.Namespace) -> None:
    """Cross-validate an adopted fit and store the result in its record."""
    from sbmlsim.credibility.cross_validation import (
        CROSS_VALIDATION_FILE,
        cross_validate,
    )

    try:
        result = cross_validate(
            model,
            options.cross_validation,
            options.experiments,
            output_dir=options.output_dir / "cross_validation",
            seed=options.seed if options.seed is not None else 1234,
            n_cores=options.cores,
        )
    except (ValueError, KeyError, FileNotFoundError) as err:
        console.print(str(err), style="error")
        sys.exit(1)
    # stored in the record of its adopted fit, which the assessment reads
    # (`sbmlsim.credibility.read_cross_validations`)
    path = model.adopted_record(result.fit) / CROSS_VALIDATION_FILE
    write_text_atomic(path, json.dumps(result.to_dict(), indent=2) + "\n")
    fit_display.section("Cross-validation", icon=fit_display.ICON_FIT)
    for study in result.studies:
        console.print(
            f"{study.study}: GMFE {study.gmfe_held_out:.3g} held out, "
            f"{study.gmfe_fitted:.3g} fitted ({study.n_points} points)"
        )
    console.print(result.statement())
    display.key_values({"cross-validation": path})


# ----------------------------------------------------------------------------
# run_report and run_identifiability
# ----------------------------------------------------------------------------
def _add_report_arguments(
    parser: argparse.ArgumentParser, definitions: Mapping[str, FitDefinition]
) -> None:
    """Add the options of a report, those of `sbmlsim.fit.cli`."""
    keys = list(definitions)
    parser.add_argument(
        "-x", "--subset", choices=keys, default=keys[0], help="fit problem to run"
    )
    parser.add_argument(
        "-n", "--name", default=None, help="name of the report, the fit id by default"
    )
    parser.add_argument(
        "--no-mapping-figures",
        action="store_true",
        help="leave the figure of every fit mapping out of the report. "
        "They are almost the whole cost of it, so a report which is read for "
        "its tables and its overview figures is created much faster",
    )
    parser.add_argument(
        "--group-by",
        action="append",
        default=None,
        metavar="COLUMNS",
        help="evaluate the report by the values of columns of its data points, "
        "a comma separated list per grouping, e.g. tissue or tissue,route; "
        "repeatable, replaces the groupings of the definition",
    )


def _group_by(options: argparse.Namespace) -> list[list[str]] | None:
    """Get the groupings of a report spec from `--group-by`."""
    if not options.group_by:
        return None
    return [
        [column.strip() for column in value.split(",")] for value in options.group_by
    ]


def report_parser(
    model: PKPDModel, definitions: Mapping[str, FitDefinition]
) -> argparse.ArgumentParser:
    """Get the parser of the `report` command of a model, of parameter files."""
    p = argparse.ArgumentParser(
        prog=f"run_report {model.package}",
        description=(
            "Report of a fit for stored parameters. `--adopted` instead "
            "writes the report of every adopted fit of the model, without "
            "a parameters file."
        ),
    )
    p.add_argument(
        "parameters",
        type=Path,
        nargs="+",
        help="JSON files with the parameter sets to report",
    )
    _add_report_arguments(p, definitions)
    p.add_argument(
        "-o",
        "--output_dir",
        type=Path,
        default=model.report_dir,
        help="directory for the report",
    )
    return p


def _has_records(parameters: Path) -> bool:
    """Check whether a parameter file is the `parameters.json` of the records of a fit."""
    return parameters.name == PARAMETERS_FILE and parameters.is_file()


def _adopted_reports(model: PKPDModel, args: Sequence[str]) -> Path:
    """Write the report of every adopted fit of a model, `run_report --adopted`.

    No fit is run: the kind `report` with `adopted` writes the report of an
    adopted fit from its record into `<output_dir>/<fit id>`, with the
    problem of its subset, the settings the fit ran with and the context of
    the model.

    Returns:
        The directory of the reports.
    """
    p = argparse.ArgumentParser(
        prog=f"run_report {model.package}",
        description="Report of every adopted fit of the model, from its record.",
    )
    p.add_argument("--adopted", action="store_true", required=True)
    p.add_argument(
        "-o",
        "--output_dir",
        type=Path,
        default=model.report_dir,
        help="directory of the reports, one folder per fit",
    )
    p.add_argument("--no-mapping-figures", action="store_true")
    options = p.parse_args(args)
    log.enable_rich_logging()
    if not model.adopted_ids:
        _failed(f"The model '{model.name}' has no adopted fits.")
    reference = project_reference(model)
    directories = [options.output_dir / adopted_id for adopted_id in model.adopted_ids]
    for out in directories:
        _check_out(p, model, out)
    for adopted_id, out in zip(model.adopted_ids, directories, strict=True):
        spec = {
            "project": reference,
            "subset": fit_subset(adopted_id),
            "adopted": True,
            "mapping_figures": not options.no_mapping_figures,
        }
        _staged(
            out,
            lambda staging, spec=spec: _execute(p, "report", spec, staging),
            inner=adopted_id,
        )
    for directory in directories:
        display.link(directory.name, directory / "index.html")
    return options.output_dir


def report_command(model: PKPDModel, args: Sequence[str] | None = None) -> Path:
    """Report stored parameter sets of a model.

    The kind `report` writes the report of the records next to a
    `parameters.json` into `<output_dir>/<name>`, `report` by default, with
    the context of the model, which names the fixed adopted fits of the record
    (`sbmlsim.fit.cli.records_context`). Several parameter files, or a file
    which is no `parameters.json` of the records of a fit, are reported with
    `sbmlsim.fit.cli.report_cli`. With `--adopted` the report of every
    adopted fit of the model is written from its record, see
    `_adopted_reports`.

    Args:
        model: the model, the model object of its package.
        args: command line arguments, `sys.argv[1:]` by default.

    Returns:
        The directory of the report, with `--adopted` the directory of the
        reports.
    """
    arguments = sys.argv[1:] if args is None else list(args)
    if "--adopted" in arguments:
        return _adopted_reports(model, arguments)
    definitions = model.load_fit_definitions()
    p = report_parser(model, definitions)
    options = p.parse_args(arguments)
    if len(options.parameters) > 1 or not _has_records(options.parameters[0]):
        return report_cli(
            definitions,
            args=with_output_dir(arguments, model.report_dir),
            prog=f"run_report {model.package}",
            description=p.description or "",
            context=records_context(model),
        )
    log.enable_rich_logging()
    reference = project_reference(model)
    name = _plain_name(p, options.name or REPORT_NAME, "--name")
    out = options.output_dir / name
    (parameters,) = options.parameters
    spec = {
        "project": reference,
        "subset": options.subset,
        "mapping_figures": not options.no_mapping_figures,
        "group_by": _group_by(options),
    }
    _check_out(p, model, out, [parameters.parent])
    inputs = {"fit": parameters.parent}
    _staged(
        out,
        lambda staging: _execute(p, "report", spec, staging, inputs=inputs),
        inner=REPORT_NAME,
    )
    # the report of the kind named its folder before it was moved
    display.link("results", out / "index.html")
    return out


def identifiability_parser(
    model: PKPDModel, definitions: Mapping[str, FitDefinition]
) -> argparse.ArgumentParser:
    """Get the parser of the `identifiability` command of a model.

    The options of `sbmlsim.fit.cli.identifiability_cli`.
    """
    defaults = ProfileSettings()
    p = argparse.ArgumentParser(
        prog=f"run_identifiability {model.package}",
        description="Identifiability of stored parameters by profile likelihood.",
    )
    p.add_argument(
        "parameters",
        type=Path,
        help="JSON file with the parameter set to analyse, the first set is used",
    )
    _add_report_arguments(p, definitions)
    p.add_argument(
        "-c", "--cores", type=int, default=1, help="number of cores for the scans"
    )
    p.add_argument(
        "-p",
        "--parameter",
        action="append",
        default=None,
        help=(
            "parameter to profile, all parameters which are no elements of a "
            "network by default; repeatable"
        ),
    )
    p.add_argument(
        "-a",
        "--alpha",
        type=float,
        default=defaults.alpha,
        help="confidence level of the intervals",
    )
    p.add_argument(
        "--df",
        type=int,
        default=defaults.degrees_of_freedom,
        help="degrees of freedom of the threshold, 1 for pointwise intervals",
    )
    p.add_argument(
        "--max-points",
        type=int,
        default=defaults.max_points,
        help="largest number of points of a profile in one direction",
    )
    p.add_argument(
        "--initial-step",
        type=float,
        default=defaults.initial_step,
        help="first step of a scan in decades of the parameter",
    )
    p.add_argument(
        "--min-step",
        type=float,
        default=defaults.min_step,
        help="smallest step of a scan in decades",
    )
    p.add_argument(
        "--max-step",
        type=float,
        default=defaults.max_step,
        help="largest step of a scan in decades",
    )
    p.add_argument(
        "--no-reoptimize",
        action="store_true",
        help="keep the other parameters at the optimum, i.e., scan the cost "
        "instead of computing the profile likelihood",
    )
    p.add_argument(
        "--no-fisher",
        action="store_true",
        help="skip the Fisher information, i.e., report the profiles alone",
    )
    p.add_argument(
        "-o",
        "--output_dir",
        type=Path,
        default=model.identifiability_dir,
        help="directory for the report",
    )
    return p


def identifiability_command(
    model: PKPDModel, args: Sequence[str] | None = None
) -> Path:
    """Analyse the identifiability of stored parameters of a model.

    The kind `identifiability` profiles the parameters of the records next to
    a `parameters.json` into `analysis/` of `<output_dir>/<name>`
    (`identifiability` by default), with the Fisher information unless
    `--no-fisher` is given, and the kind `report` writes the report of the
    analysis into the folder, the records of its run in `report/`. A file
    which is no `parameters.json` of the records of a fit is analysed with
    `sbmlsim.fit.cli.identifiability_cli`.

    Args:
        model: the model, the model object of its package.
        args: command line arguments, `sys.argv[1:]` by default.

    Returns:
        The directory of the report.
    """
    arguments = sys.argv[1:] if args is None else list(args)
    definitions = model.load_fit_definitions()
    p = identifiability_parser(model, definitions)
    options = p.parse_args(arguments)
    if not _has_records(options.parameters):
        return identifiability_cli(
            definitions,
            args=with_output_dir(arguments, model.identifiability_dir),
            prog=f"run_identifiability {model.package}",
        )
    log.enable_rich_logging()
    reference = project_reference(model)
    name = _plain_name(p, options.name or "identifiability", "--name")
    out = options.output_dir / name
    _check_out(p, model, out, [options.parameters.parent])
    spec = {
        "project": reference,
        "subset": options.subset,
        "alpha": options.alpha,
        "degrees_of_freedom": options.df,
        "max_points": options.max_points,
        "initial_step": options.initial_step,
        "min_step": options.min_step,
        "max_step": options.max_step,
        "reoptimize": not options.no_reoptimize,
        "pids": options.parameter,
        "fisher": not options.no_fisher,
    }
    report_spec = {
        "project": reference,
        "subset": options.subset,
        "mapping_figures": not options.no_mapping_figures,
        "group_by": _group_by(options),
    }

    def run(staging: Path) -> Path:
        analysis = staging / ANALYSIS_DIR
        _execute(
            p,
            "identifiability",
            spec,
            analysis,
            workers=options.cores,
            inputs={"fit": options.parameters.parent},
        )
        report_run = staging / REPORT_RUN_DIR
        _execute(p, "report", report_spec, report_run, inputs={"fit": analysis})
        return analysis

    staged = _staged(out, run, inner=f"{REPORT_RUN_DIR}/{REPORT_NAME}")
    # the report read the analysis in the staging directory: its spec names
    # the analysis where it is now
    _relocated_input(out / REPORT_RUN_DIR / SPEC_FILE, staged, out / ANALYSIS_DIR)
    # the report of the kind named its folder before it was moved
    display.link("results", out / "index.html")
    return out


def credibility_command(model: PKPDModel, args: Sequence[str] | None = None) -> Path:
    """Assess the contexts of use of a model against its adopted fits.

    Runs the kind `credibility` of `sbmlsim.run` (`sbmlsim.credibility.assess`
    with the evidence extensions and the PK battery of the model, see
    `pkpdsim.model.CREDIBILITY_SETTINGS`) into the output directory through a
    staged swap: the outputs of the kind (`credibility/` of its output) are
    moved up into the directory, next to the records of the run. The adopted
    fits are evaluated on the current model and data, `--from-json` instead
    reads a written record (the input `record`) without evaluating a fit. The
    decisive predictions are verified (`sbmlsim.credibility.verification`)
    unless `--no-verification` is given. `--lock` writes the lock of the
    evidence, which the kind writes into its output, into the folder of the
    model. `--source-commit` records the commit of the source repository in
    the provenance, which a release passes, since it assesses the model in its
    release repository. The outputs replace an earlier assessment in the
    directory, a directory with other files is refused; an error of the
    declaration, the evidence, the records or the directory is printed and
    exits with 1.

    Args:
        model: the model.
        args: command line arguments, `sys.argv[1:]` by default.

    Returns:
        The directory of the outputs.
    """
    p = argparse.ArgumentParser(
        prog=f"run_credibility {model.package}",
        description=(
            f"Assess the contexts of use of the {model.title} against its "
            "adopted fits, with the verification of the decisive predictions: "
            "the record credibility.json, the tables contexts.tsv, evidence.tsv, "
            "domain.tsv and verification.tsv and the report index.html."
        ),
    )
    p.add_argument(
        "-o",
        "--output_dir",
        type=Path,
        default=model.credibility_dir,
        help="directory of the outputs, replaced if it holds an earlier "
        "assessment (default: %(default)s)",
    )
    p.add_argument(
        "--from-json",
        type=Path,
        metavar="PATH",
        help="write the outputs again from a written credibility.json, without "
        "evaluating the fits",
    )
    p.add_argument(
        "--no-verification",
        action="store_true",
        help="do not verify the decisive predictions with tightened tolerances, "
        "a second integration method and the checks of the builds; the factors "
        "of the verification are then not assessed",
    )
    p.add_argument(
        "--lock",
        action="store_true",
        help="write the lock of the evidence, credibility.lock.json in the folder "
        "of the model, instead of comparing the evidence with it: the splits "
        "fixed for review",
    )
    p.add_argument(
        "--uncertainty",
        action="store_true",
        help="propagate the parameter uncertainty of the adopted fits to the "
        "decisive predictions: every fit is evaluated at the draws of its "
        "parameters, the factor model inputs is then assessed",
    )
    p.add_argument(
        "--draws",
        type=int,
        default=None,
        metavar="N",
        help="the draws of every fit of --uncertainty (default 128)",
    )
    p.add_argument(
        "--seed",
        type=int,
        default=None,
        metavar="SEED",
        help="the seed of the draws of --uncertainty (default 1234)",
    )
    p.add_argument(
        "--variability",
        action="store_true",
        help="check the decisive points with the variability of the subjects the "
        "declaration has (Credibility.variability): the visual predictive check, "
        "a mean of a group against the means of virtual groups of its size, "
        "nested in the draws of --uncertainty",
    )
    p.add_argument(
        "--subjects",
        type=int,
        default=None,
        metavar="N",
        help="the subjects of every study of --variability (default 64)",
    )
    p.add_argument(
        "--studies",
        type=int,
        default=None,
        metavar="N",
        help="the studies of every parameter vector of --variability with a "
        "variability between the studies (default 16)",
    )
    p.add_argument(
        "--source-commit",
        metavar="SHA",
        help="the commit of the source repository of the model, recorded in the "
        "provenance instead of the commit of the working directory (whether its "
        "tree was clean is not known); a release passes it",
    )
    options = p.parse_args(sys.argv[1:] if args is None else list(args))
    if options.from_json is not None and options.source_commit is not None:
        p.error("--source-commit is the provenance of an assessment, not of a record")
    if options.from_json is not None and options.no_verification:
        p.error("--no-verification is an option of an assessment, not of a record")
    if options.from_json is not None and options.lock:
        p.error("--lock is an option of an assessment, not of a record")
    if options.from_json is not None and options.uncertainty:
        p.error("--uncertainty is an option of an assessment, not of a record")
    if options.from_json is not None and options.variability:
        p.error("--variability is an option of an assessment, not of a record")
    if not options.variability and (
        options.subjects is not None or options.studies is not None
    ):
        p.error("--subjects and --studies are options of --variability")
    if not options.uncertainty and (
        options.draws is not None or options.seed is not None
    ):
        p.error("--draws and --seed are options of --uncertainty")
    out = options.output_dir
    _check_out(p, model, out)
    _check_credibility_out(out)
    spec: dict[str, Any] = {
        "project": project_reference(model),
        "verification": not options.no_verification,
        "lock": options.lock,
        "uncertainty": options.uncertainty,
        "draws": options.draws,
        "seed": options.seed,
        "variability": options.variability,
        "subjects": options.subjects,
        "studies": options.studies,
        "source_commit": options.source_commit,
    }
    inputs = None if options.from_json is None else {"record": options.from_json}

    try:
        _staged(
            out,
            lambda staging: _execute(p, "credibility", spec, staging, inputs=inputs),
            inner=CREDIBILITY_DIR,
        )
    except (LookupError, ValueError, TypeError, OSError) as err:
        # the errors of the declaration, the evidence, the records and the
        # lock, e.g. `DeclarationError`, `SelectionError`, `KindMismatchError`
        _failed(str(err))
    if options.lock:
        from sbmlsim.credibility import lock as locks

        # the lock of the package, which the next assessment compares with
        lock = locks.lock_path(model)
        try:
            shutil.move(out / locks.LOCK_FILE, lock)
        except OSError as err:
            _failed(
                f"The lock '{lock}' is not written ({err}); the assessment is in "
                f"'{out}'."
            )
        display.key_values({"lock": lock})
    display.key_values({"report": out / "index.html"})
    return out


def _check_credibility_out(out: Path) -> None:
    """Refuse an output of an assessment which holds other files than an assessment.

    The outputs replace an earlier assessment, a directory with
    `credibility.json`, as `sbmlsim.credibility.report.write_outputs` does.

    Raises:
        SystemExit: exit 1 for a file or a directory with other files.
    """
    if out.exists() and not out.is_dir():
        _failed(f"The output directory '{out}' is a file.")
    if out.is_dir() and any(out.iterdir()) and not (out / CREDIBILITY_RECORD).is_file():
        _failed(
            f"The output directory '{out}' is not empty and holds no "
            f"{CREDIBILITY_RECORD}; it is not replaced."
        )


DOTTED_PACKAGE = re.compile(r"^[A-Za-z_]\w*(\.[A-Za-z_]\w*)*$")
"""The form of a dotted package, e.g. `pkdb_models.models.hctz`."""


class ModelNotFoundError(LookupError):
    """A package which is not found or is no model package."""


def resolve_model(package: str) -> PKPDModel:
    """Get the model object `MODEL` of a model package.

    Args:
        package: the dotted name of the model package, e.g.
            `pkdb_models.models.hctz`.

    Returns:
        The model object of the package.

    Raises:
        ModelNotFoundError: if the name is no dotted package, the package is not
            found, it has no `MODEL`, its `MODEL` is no `PKPDModel` or belongs
            to another package. Any other error of the import of the package,
            e.g. a missing module which the package imports, is raised as it is,
            with its traceback.
    """
    if not DOTTED_PACKAGE.match(package):
        raise ModelNotFoundError(f"'{package}' is no dotted package.")
    try:
        module = importlib.import_module(package)
    except ModuleNotFoundError as err:
        # only the package itself or one of its parents is not found
        if err.name is not None and (
            package == err.name or package.startswith(f"{err.name}.")
        ):
            raise ModelNotFoundError(
                f"The model package '{package}' is not found."
            ) from err
        raise
    model = getattr(module, "MODEL", None)
    if model is None:
        raise ModelNotFoundError(f"The package '{package}' has no model object MODEL.")
    if not isinstance(model, PKPDModel):
        raise ModelNotFoundError(
            f"MODEL of the package '{package}' is no PKPDModel but a "
            f"{type(model).__name__}."
        )
    if model.package != package:
        raise ModelNotFoundError(
            f"MODEL of the package '{package}' belongs to the package "
            f"'{model.package}'."
        )
    return model


@dataclasses.dataclass(frozen=True)
class Command:
    """A command of the models, a console script of pkpdsim.

    Attributes:
        name: the name of the console script.
        description: what the command does, the description of its help.
        function: the command as a function of the model and the arguments
            after the model package.
        fits: the command needs the fit definitions of the model.
    """

    name: str
    description: str
    function: Callable[[PKPDModel, Sequence[str]], object]
    fits: bool


COMMANDS: dict[str, Command] = {
    command.name: command
    for command in (
        Command(
            name="run_simulation",
            description="Run the model factory, the simulation experiments and "
            "the extra actions of a model.",
            function=run_command,
            fits=False,
        ),
        Command(
            name="run_fit",
            description="Fit the parameters of a model and analyse their "
            "identifiability, or adopt a finished fit.",
            function=fit_command,
            fits=True,
        ),
        Command(
            name="run_report",
            description="Report stored parameters of a model.",
            function=report_command,
            fits=True,
        ),
        Command(
            name="run_identifiability",
            description="Analyse the identifiability of stored parameters of a model.",
            function=identifiability_command,
            fits=True,
        ),
        Command(
            name="run_credibility",
            description="Assess the contexts of use of a model against its "
            "adopted fits.",
            function=credibility_command,
            fits=True,
        ),
    )
}
"""The commands by the names of their console scripts."""


def profile_dir(model: PKPDModel, command: str) -> Path:
    """Get a new folder of the profile of a command in the results of a model.

    Args:
        model: the model.
        command: the name of the command, a key of `COMMANDS`.

    Returns:
        `profiling/<command>_<time>` of the results of the model, not created.
    """
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    return model.results_dir / PROFILE_DIR / f"{command}_{stamp}"


def command_parser(command: Command) -> argparse.ArgumentParser:
    """Get the parser of a command before its model is known.

    The options of a command depend on the model, so this parser only prints
    the general usage and the errors of the model package.
    """
    p = argparse.ArgumentParser(
        prog=command.name,
        usage=(
            f"{command.name} PACKAGE [{PROFILE_OPTION}] "
            f"[{PROGRESS_OPTION} MODE] [options]"
        ),
        description=command.description,
        epilog=f"`{command.name} PACKAGE --help` shows the options of the model.",
    )
    p.add_argument(
        "package",
        metavar="PACKAGE",
        help="the dotted name of the model package, which defines the model "
        "object MODEL, e.g. pkdb_models.models.hctz",
    )
    p.add_argument(
        PROFILE_OPTION,
        action="store_true",
        help="run the command under the sampling profiler of python 3.15 and "
        "write its flame graph and statistics into profiling/ of the results "
        "of the model (`sbmlsim.run.profiling`)",
    )
    p.add_argument(
        PROGRESS_OPTION,
        metavar="MODE",
        choices=PROGRESS_NAMES,
        help="how the progress of the runs of the command is shown: auto (the "
        "bars of a terminal, a line of the log otherwise), bar, log or none; "
        "the events of the progress are written whatever is chosen "
        "(`sbmlsim.progress`)",
    )
    return p


def _progress_setting(args: Sequence[str]) -> list[str]:
    """Take the progress out of the arguments of a command and set it.

    `--progress MODE` and `--progress=MODE` anywhere after the model package
    set the environment variable `SBMLSIM_PROGRESS`, which every run of the
    command and every process it starts reads, see `sbmlsim.progress.resolve`.
    The option is taken out of the arguments, so the parser of the command
    does not see it, as `--profile` is.

    Args:
        args: the arguments after the model package.

    Returns:
        The arguments without the option.

    Raises:
        ValueError: for a mode which is no rendering of the progress and for
            the option without one.
    """
    rest: list[str] = []
    args = list(args)
    index = 0
    while index < len(args):
        argument = args[index]
        if argument == PROGRESS_OPTION:
            if index + 1 == len(args):
                raise ValueError(f"{PROGRESS_OPTION} needs one of {PROGRESS_NAMES}")
            value, index = args[index + 1], index + 2
        elif argument.startswith(f"{PROGRESS_OPTION}="):
            value, index = argument.split("=", 1)[1], index + 1
        else:
            rest.append(argument)
            index += 1
            continue
        if value.strip().lower() not in PROGRESS_NAMES:
            raise ValueError(
                f"the progress '{value}' is none of {PROGRESS_NAMES} "
                f"({PROGRESS_OPTION})"
            )
        os.environ[PROGRESS_ENV] = value.strip().lower()
    return rest


def main(name: str, argv: Sequence[str] | None = None) -> None:
    """Run a command with the model package as the first argument.

    The model is resolved first, since the options of the command depend on it;
    the arguments after the package are the arguments of the command.

    Args:
        name: the name of the command, a key of `COMMANDS`.
        argv: the command line arguments, `sys.argv[1:]` by default.
    """
    command = COMMANDS[name]
    args = list(sys.argv[1:] if argv is None else argv)
    p = command_parser(command)
    if not args:
        p.error("the model package is required")
    if args[0] in ("-h", "--help"):
        p.print_help()
        p.exit()
    if args[0].startswith("-"):
        p.error(f"the model package is the first argument, not '{args[0]}'")
    package, *rest = args
    try:
        rest = _progress_setting(rest)
    except ValueError as err:
        p.error(str(err))
    try:
        model = resolve_model(package)
    except ModelNotFoundError as err:
        p.error(str(err))
    if command.fits and not model.fits:
        p.error(f"the model '{model.name}' ({package}) has no fit definitions")
    if PROFILE_OPTION in rest:
        rest = [arg for arg in rest if arg != PROFILE_OPTION]
        if not profiler_available():
            p.error(
                f"{PROFILE_OPTION} needs the sampling profiler of python 3.15, "
                f"this is python {sys.version.split()[0]}"
            )
        code = run_profiled(
            [name, package, *rest], profile_dir(model, name), module="pkpdsim.cli"
        )
        if code:
            sys.exit(code)
        return
    command.function(model, rest)


def run_simulation_main(argv: Sequence[str] | None = None) -> None:
    """The console script `run_simulation`, see `run_command`."""
    main("run_simulation", argv)


def run_fit_main(argv: Sequence[str] | None = None) -> None:
    """The console script `run_fit`, see `fit_command`.

    Returns nothing: a console script exits with the return value of the
    function it calls, and 1 for anything which is not `None` or an int.
    """
    main("run_fit", argv)


def run_report_main(argv: Sequence[str] | None = None) -> None:
    """The console script `run_report`, see `report_command`."""
    main("run_report", argv)


def run_identifiability_main(argv: Sequence[str] | None = None) -> None:
    """The console script `run_identifiability`, see `identifiability_command`."""
    main("run_identifiability", argv)


def run_credibility_main(argv: Sequence[str] | None = None) -> None:
    """The console script `run_credibility`, see `credibility_command`."""
    main("run_credibility", argv)


if __name__ == "__main__":
    # `python -m pkpdsim.cli <command> <package> ...`, which the profiler runs;
    # its exit code goes to the profiler through a file (`record_exit_code`)
    try:
        main(sys.argv[1], sys.argv[2:])
    except SystemExit as stop:
        record_exit_code(stop.code)
        raise
    except BaseException:
        record_exit_code(1)
        raise
    record_exit_code(0)
