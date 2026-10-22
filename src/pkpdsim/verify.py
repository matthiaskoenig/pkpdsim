"""Verify that a change of a model changes only what it is meant to change.

A refactoring of a model, an update of sbmlsim or a refresh of the study data
should leave the simulations of the model as they were, and a fix should change
them only where it is meant to. The subcommands make this checkable:

- `dump` runs the simulation experiments of a model in a checkout and stores the
  results, the datasets, the fit mappings and the figures in a pickle,
- `compare` compares two dumps: the simulations on the matched rows of their
  output times, relative to the scale of every array (a time which only one
  dump has is a difference), everything else exactly; with `--renames FILE` the
  expected renames of a migration (metadata fields and members, added fields,
  texts such as legend keys, ids) and the expected differences are applied to
  the old dump first, see `load_renames`,
- `ranges` lists the fit mappings of a dump whose data lie outside of the
  simulated range, e.g. data times beyond the end of the simulation,
- `fits` builds the optimization problem of every fit definition of a model and
  runs a minimal fit,
- `definitions` writes the fit definitions of a model as json: the settings,
  the parameters and the mappings of the collections, so that two checkouts
  are compared with `diff`,
- `settings` writes what a run of the experiments of a model hands to sbmlsim, so
  that the settings of two checkouts are compared with `diff`,
- `png` compares the PNG figures of two directories pixel by pixel,
- `sbml` compares SBML models by their elements, ignoring the order and the
  annotations.

A model before and after a change, with tight tolerances so that the solver does
not hide a difference:

```bash
git worktree add --detach ../base develop
# the experiments simulate the model at `MODEL_PATH`, which is generated:
# create it with the factory of each checkout right before its dump
(cd ../base && uv run python pkdb_models/models/hctz/models/factory.py)
uv run python pkdb_models/models/hctz/models/factory.py
uv run pkpdsim-verify dump --root ../base --model pkdb_models.models.hctz \
    --out base.pkl --atol 1e-16 --rtol 1e-12
uv run pkpdsim-verify dump --root . --model pkdb_models.models.hctz \
    --out new.pkl --order base.pkl --atol 1e-16 --rtol 1e-12
uv run pkpdsim-verify compare base.pkl new.pkl
```

`--root` is the directory which is put first on the python path, e.g. a checkout
of the repository of the model, so that its packages are imported. `--model` is
the dotted name of a model package, which defines its model object `MODEL`
(`pkpdsim.model.PKPDModel`) which describes it; a package without a model object
`MODEL` is a usage error (exit code 2). The command only needs the
python packages of the environment it runs in, so a checkout which needs another
version of sbmlsim is dumped with the python of an environment with that
version, as long as it is one which pkpdsim supports.

Across versions of sbmlsim (a dump of pkpdsim 0.4 with sbmlsim 0.8 and one of
0.5 with 0.9) `compare` moves the time of every array first and matches the
other dimensions by name, matches the rows of a repeated time from the last
one, so that the state after a change, which is all sbmlsim 0.8 gave, is
compared, and counts the rows of a time which only one dump has as differences
(a new output time of a task is expected with the pattern of its variable
`time`), and compares a figure of sbmlsim 0.8 with the one of sbmlsim 0.9 in the
form of sbmlsim 0.9: the settings and the size it is rendered with, without the
fields which only sbmlsim 0.9 has, which it lists (`align_figures`). The renames
of a port are the ids which became SIds (`sids`) and the differences which are
expected (`expected`).
"""

import argparse
import ast
import fnmatch
import importlib
import importlib.metadata
import inspect
import json
import math
import os
import pickle
import pkgutil
import sys
import tempfile
import traceback
from collections.abc import Callable, Collection, Iterable, Mapping, Sequence
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import (
    TYPE_CHECKING,
    Any,
    Literal,
    NamedTuple,
    NotRequired,
    TypedDict,
    cast,
)

import libsbml
import numpy as np

if TYPE_CHECKING:
    from sbmlsim.experiment import SimulationExperiment
    from sbmlsim.fit import FitData, FitMapping
    from sbmlsim.fit.cli import FitDefinition
    from sbmlsim.plot import Figure

ROOT_HELP = (
    "directory put first on the python path, e.g. a checkout of the repository "
    "of the model"
)
MODEL_HELP = "dotted name of the model package, e.g. pkdb_models.models.hctz"

#: dimensions along which the simulations of a task are compared: the time of a
#: grid (sbmlsim 0.9), the interpolated time and the points of a timecourse
#: (sbmlsim 0.8)
TIME_DIMENSIONS = ("time", "_time", "_point")

#: experiment groups of the registry which are run first, in this order
GROUPS = ("studies", "misc", "scans")

#: the settings of a figure of sbmlsim 0.8: attributes of the class `Figure`,
#: which a figure overrode as attributes of its own; sbmlsim 0.9 has them as
#: the fields of `FigureSettings`, which a figure overrides in its `_settings`
FIGURE_SETTINGS_08: tuple[str, ...] = (
    "fig_dpi",
    "fig_facecolor",
    "fig_subplots_wspace",
    "fig_subplots_hspace",
    "panel_width",
    "panel_height",
    "fig_titlesize",
    "fig_titleweight",
    "axes_titlesize",
    "axes_titleweight",
    "axes_labelsize",
    "axes_labelweight",
    "xtick_labelsize",
    "ytick_labelsize",
    "legend_fontsize",
    "legend_position",
    "legend_loc",
)

#: the settings of a figure which only sbmlsim from 0.9.6 has (the scope of a
#: legend, `legend_scope`): sbmlsim 0.8 has no counterpart, and a figure which
#: does not override them has the same settings in both
NEW_FIGURE_SETTINGS: tuple[str, ...] = ("legend_scope",)

#: the key of the specification of a figure (sbmlsim 0.9) with the settings of
#: its experiment and its own and the size it is rendered with, without the
#: settings of a run, see `figure_spec`
RESOLVED = "__resolved__"

#: the fields of the specification of a figure which only sbmlsim 0.9 has: the
#: selection of a `Data`, the dimensions of the lines of a `Curve`, the bands of
#: a `Plot`
NEW_FIGURE_FIELDS: tuple[str, ...] = ("sel", "over", "_bands")


class Dump(TypedDict):
    """Everything `dump` stores of the experiments of a model."""

    model: str
    #: names of the experiment classes in the order they were run
    order: list[str]
    #: traceback per experiment (or `import:<module>`) which failed
    errors: dict[str, str]
    #: (experiment, task, variable) -> (dimensions, values)
    results: dict[tuple[str, str, str], tuple[tuple[str, ...], np.ndarray]]
    #: (experiment, dataset, column) -> values
    datasets: dict[tuple[str, str, str], np.ndarray]
    #: experiment -> description of every fit mapping
    mappings: dict[str, list[str]]
    #: experiment -> figure -> specification of the figure (json)
    figures: dict[str, dict[str, str]]
    #: experiment -> fit mapping -> (data min, data max, simulation min, simulation
    #: max, unit) of the x values, or the error message
    ranges: dict[str, dict[str, tuple[float, float, float, float, str] | str]]
    #: the versions of the packages which made the dump (since pkpdsim 0.5)
    versions: NotRequired[dict[str, str]]


# --- the experiments of a model -------------------------------------------------


@dataclass
class ModelInfo:
    """What the verification needs of a model package."""

    package: str
    base_path: Path
    model_path: Path
    data_path: Path | list[Path] | None
    experiments: Mapping[str, Sequence[type["SimulationExperiment"]]] | None
    absolute_tolerance: object | None
    relative_tolerance: object | None
    run_experiments: Callable[..., object] | None
    fit_definitions: Mapping[str, "FitDefinition"] | None


def _outside_root(package: str, root: Path) -> str | None:
    """Check that a package was imported from below `root`.

    `--root` is put first on the python path, but a package which is already
    imported (a second call in the same process) or installed in an environment
    in which the root does not take effect would silently be dumped from
    somewhere else.

    Returns:
        A message if the package comes from elsewhere, else `None`.
    """
    module = importlib.import_module(package)
    file = getattr(module, "__file__", None)
    if file is not None and Path(file).resolve().is_relative_to(root):
        return None
    return f"{package} is imported from {file}, not from below {root}"


def describe(package: str) -> ModelInfo:
    """Describe a model package by its model object `MODEL`.

    Args:
        package: the dotted name of the model package.

    Returns:
        What the verification needs of the model.

    Raises:
        ValueError: if the package defines no model object.
    """
    from pkpdsim.model import PKPDModel

    module = importlib.import_module(package)
    model = getattr(module, "MODEL", None)
    if not isinstance(model, PKPDModel):
        raise ValueError(
            f"{package} defines no model object MODEL (pkpdsim.model.PKPDModel)"
        )
    return ModelInfo(
        package=package,
        base_path=model.path,
        model_path=model.model_path,
        data_path=list(model.data_paths) or None,
        experiments=model.load_experiments() if model.experiments else None,
        absolute_tolerance=model.simulator.absolute_tolerance,
        relative_tolerance=model.simulator.relative_tolerance,
        run_experiments=model.run_experiments,
        fit_definitions=(
            dict(model.load_fit_definitions()) if model.fit_definitions else None
        ),
    )


def _describe(args: argparse.Namespace) -> ModelInfo | int:
    """Describe the model of the arguments, or print the problem and return 2."""
    try:
        return describe(args.model)
    except ValueError as err:
        print(err, file=sys.stderr)
        return 2


def registered_experiments(
    info: ModelInfo,
) -> list[type["SimulationExperiment"]] | None:
    """Get the experiments of the registry of a model.

    Returns:
        The experiment classes, those of `GROUPS` first, without duplicates,
        or `None` if the model has no registry.
    """
    registry = info.experiments
    if registry is None:
        return None
    keys = [key for key in GROUPS if key in registry]
    keys += [key for key in registry if key not in GROUPS]
    return list(dict.fromkeys(cls for key in keys for cls in registry[key]))


def defined_experiments(
    package: str, errors: dict[str, str]
) -> list[type["SimulationExperiment"]]:
    """Get the experiments defined in the modules of `experiments/` of a model.

    Every subclass of `SimulationExperiment` defined in a module below
    `<package>.experiments` counts, except the modules `base*`
    and `metadata` and the classes from which another found class derives (the
    base experiments).

    Args:
        package: the dotted name of the model package.
        errors: the traceback of every module which does not import is added
            here, as `import:<module>`.

    Returns:
        The experiment classes in the order of module name and definition.
    """
    from sbmlsim.experiment import SimulationExperiment

    def on_error(name: str) -> None:
        errors[f"import:{name}"] = traceback.format_exc()

    experiments = importlib.import_module(f"{package}.experiments")
    found: list[tuple[str, int, type[SimulationExperiment]]] = []
    for info in pkgutil.walk_packages(
        experiments.__path__, prefix=f"{experiments.__name__}.", onerror=on_error
    ):
        last = info.name.rsplit(".", 1)[-1]
        if info.ispkg or last.startswith("base") or last == "metadata":
            continue
        try:
            module = importlib.import_module(info.name)
        except Exception:  # every error of a module is reported
            on_error(info.name)
            continue
        for _, cls in inspect.getmembers(module, inspect.isclass):
            if cls.__module__ == module.__name__ and issubclass(
                cls, SimulationExperiment
            ):
                try:
                    line = inspect.getsourcelines(cls)[1]
                except (OSError, TypeError):
                    line = 0
                found.append((info.name, line, cls))
    found.sort(key=lambda item: (item[0], item[1]))
    classes = [cls for _, _, cls in found]
    return [
        cls
        for cls in classes
        if not any(other is not cls and issubclass(other, cls) for other in classes)
    ]


# --- dump -----------------------------------------------------------------------


def serialize(obj: object, depth: int = 0) -> object:
    """Serialize the attributes of an object (a figure) into json values.

    References back to the experiment and to parents are replaced by their type,
    so that the specification of a figure only contains the figure.
    """
    from sbmlsim.experiment import SimulationExperiment

    if depth > 12:
        return "<depth>"
    if obj is None or isinstance(obj, bool | int | float | str):
        return obj
    if isinstance(obj, Enum):
        return f"{type(obj).__name__}.{obj.name}"
    if isinstance(obj, SimulationExperiment):
        return f"<experiment {type(obj).__name__}>"
    if isinstance(obj, dict):
        return {str(k): serialize(v, depth + 1) for k, v in obj.items()}
    if isinstance(obj, list | tuple):
        return [serialize(v, depth + 1) for v in obj]
    if isinstance(obj, type):
        return f"<class {obj.__qualname__}>"
    if hasattr(obj, "__dict__"):
        skip = {"experiment", "_experiment", "figure", "_figure", "parent", "_parent"}
        return {"__type__": type(obj).__name__} | {
            k: f"<{type(v).__name__}>" if k in skip else serialize(v, depth + 1)
            for k, v in vars(obj).items()
        }
    return repr(obj)


def metadata_names(metadata: object) -> dict[str, str]:
    """Get the metadata of a fit mapping as {field: name}.

    Enum members are stored by name, so a change of their values or of the enum
    type (`StrEnum`) does not count as a change.
    """
    if metadata is None:
        return {}
    fields = getattr(metadata, "__dataclass_fields__", None)
    items = (
        {k: getattr(metadata, k) for k in fields} if fields else dict(vars(metadata))
    )
    return {
        k: v.name if isinstance(v, Enum) else str(v)
        for k, v in items.items()
        if not k.startswith("_")
    }


def describe_mapping(key: str, mapping: "FitMapping") -> str:
    """Describe a fit mapping by its data, its observable, weight and metadata."""

    def fit_data(d: "FitData") -> str:
        # the data of a mapping without x values (sbmlsim 0.9) have no x
        x = d.x.index if d.x else None
        sd = d.y_sd.index if d.y_sd else None
        se = d.y_se.index if d.y_se else None
        return f"{d.dset_id}|{d.task_id}|{x}|{d.y.index}|{d.count}|{sd}|{se}"

    metadata = sorted(metadata_names(mapping.metadata).items())
    weight = getattr(mapping, "weight", None)
    return (
        f"{key}|{fit_data(mapping.reference)}|{fit_data(mapping.observable)}"
        f"|weight={weight}|{metadata}"
    )


def mapping_range(
    mapping: "FitMapping",
) -> tuple[float, float, float, float, str] | str:
    """Get the range of the x values of the data and of the simulation of a mapping.

    The data are converted into the unit of the simulation.
    """
    try:
        data = mapping.reference.get_data().x
        simulation = mapping.observable.get_data().x
        if data is None or simulation is None:
            return "no x values"
        x = np.asarray(data.to(simulation.units).magnitude, dtype=float)
        s = np.asarray(simulation.magnitude, dtype=float)
        return (
            float(np.nanmin(x)),
            float(np.nanmax(x)),
            float(np.nanmin(s)),
            float(np.nanmax(s)),
            str(simulation.units),
        )
    except Exception as err:  # reported per mapping
        return f"{type(err).__name__}: {err}"


def dump(args: argparse.Namespace) -> int:
    """Run the experiments of a model and store everything in a pickle."""
    root = Path(args.root).resolve()
    sys.path.insert(0, str(root))
    from sbmlsim.experiment import ExperimentRunner
    from sbmlsim.simulator import Simulator

    if message := _outside_root(args.model, root):
        print(message, file=sys.stderr)
        return 2
    info = _describe(args)
    if isinstance(info, int):
        return info
    errors: dict[str, str] = {}
    classes = registered_experiments(info)
    if classes is None:
        classes = defined_experiments(info.package, errors)
    elif args.all:
        classes += [
            c for c in defined_experiments(info.package, errors) if c not in classes
        ]
    if args.order:
        with Path(args.order).open("rb") as f:
            rank = {name: k for k, name in enumerate(pickle.load(f)["order"])}
        classes.sort(key=lambda c: rank.get(c.__name__, len(rank)))
    if args.experiments:
        unknown = set(args.experiments) - {c.__name__ for c in classes}
        if unknown:
            print(f"unknown experiments: {sorted(unknown)}", file=sys.stderr)
            return 2
        classes = [c for c in classes if c.__name__ in args.experiments]
    classes = [c for c in classes if c.__name__ not in args.skip]
    print(f"{len(classes)} experiments", file=sys.stderr)

    tolerances: dict[str, Any] = {}
    if args.atol == "model":
        if info.absolute_tolerance is None:
            print(f"{args.model} has no tolerance of its own", file=sys.stderr)
            return 2
        tolerances["absolute_tolerance"] = info.absolute_tolerance
    elif args.atol is not None:
        try:
            tolerances["absolute_tolerance"] = float(args.atol)
        except ValueError:
            print(f"--atol is a number or `model`, not {args.atol}", file=sys.stderr)
            return 2
    if args.rtol == "model":
        if info.relative_tolerance is None:
            print(f"{args.model} has no relative tolerance of its own", file=sys.stderr)
            return 2
        tolerances["relative_tolerance"] = info.relative_tolerance
    elif args.rtol is not None:
        try:
            tolerances["relative_tolerance"] = float(args.rtol)
        except ValueError:
            print(f"--rtol is a number or `model`, not {args.rtol}", file=sys.stderr)
            return 2

    import pkpdsim

    result = Dump(
        model=args.model,
        order=[c.__name__ for c in classes],
        errors=errors,
        results={},
        datasets={},
        mappings={},
        figures={},
        ranges={},
        versions={
            "sbmlsim": importlib.metadata.version("sbmlsim"),
            "pkpdsim": pkpdsim.__version__,
        },
    )
    for cls in classes:
        try:
            runner = ExperimentRunner(
                experiment_classes=[cls],
                data_path=info.data_path,
                base_path=info.base_path,
                # a dump runs serially, the order of its experiments is the
                # order of the run
                simulator=Simulator(n_workers=1, **tolerances),
                # an experiment which cannot be initialized raises and is
                # recorded below; sbmlsim 0.9 only logs it by default
                on_error="raise",
            )
            for sid, experiment in runner.experiments.items():
                experiment.run(
                    runner.simulator,
                    output_path=None,
                    reduced_selections=args.reduced_selections,
                )
                store_experiment(result, sid, experiment)
        except Exception:  # every error of an experiment is reported
            result["errors"][cls.__name__] = traceback.format_exc()
            print(
                f"ERROR {cls.__name__}\n{result['errors'][cls.__name__]}",
                file=sys.stderr,
            )
            continue
        print(f"{cls.__name__}: ok", file=sys.stderr)

    with Path(args.out).open("wb") as f:
        pickle.dump(result, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(
        f"{len(classes)} experiments, {len(result['results'])} arrays, "
        f"{len(result['errors'])} errors -> {args.out}"
    )
    return 1 if result["errors"] else 0


def store_experiment(
    result: Dump, sid: str, experiment: "SimulationExperiment"
) -> None:
    """Store the results, datasets, fit mappings and figures of a run experiment.

    The time of a grid is a coordinate of a result of sbmlsim 0.9; it is stored
    as the variable `time` over `("time",)`, which the comparison reads.
    """
    for task, res in experiment.results.items():
        ds = res.ds
        for var in ds.data_vars:
            result["results"][sid, task, str(var)] = (
                tuple(str(d) for d in ds[var].dims),
                np.asarray(ds[var].values),
            )
        if "time" in ds.coords and "time" not in ds.data_vars:
            result["results"][sid, task, "time"] = (
                ("time",),
                np.asarray(ds["time"].values),
            )
    for key, dset in experiment.datasets().items():
        for column in dset.columns:
            values = np.asarray(dset[column].values)
            if values.dtype.kind not in "fiub":
                values = values.astype(str)
            result["datasets"][sid, key, str(column)] = values
    mappings = experiment.fit_mappings()
    result["mappings"][sid] = sorted(
        describe_mapping(k, m) for k, m in mappings.items()
    )
    result["ranges"][sid] = {k: mapping_range(m) for k, m in mappings.items()}
    result["figures"][sid] = {
        fid: figure_spec(fig, experiment) for fid, fig in experiment.figures().items()
    }


def figure_spec(figure: "Figure", experiment: "SimulationExperiment") -> str:
    """Get the specification of a figure, its attributes as json.

    `RESOLVED` adds the settings of the experiment and of the figure (the
    `figure_settings` of the experiment overridden by the `settings` of the
    figure) and the width and height the figure is rendered with, as sbmlsim
    resolves them without the settings of a run. A dump of sbmlsim 0.8 has
    them as attributes of the figure (`_width`, `_height` and a setting such as
    `legend_position`), so that `compare` compares them across versions.
    """
    spec = cast("dict[str, object]", serialize(figure))
    width, height = figure.size(figure.resolved_settings(experiment))
    spec[RESOLVED] = {
        "settings": serialize({**experiment.figure_settings, **figure.settings}),
        "width": width,
        "height": height,
    }
    return json.dumps(spec, sort_keys=True, default=repr)


# --- compare --------------------------------------------------------------------


def load(path: str) -> Dump:
    """Load a dump."""
    with Path(path).open("rb") as f:
        return pickle.load(f)


def time_axis(dims: Sequence[str]) -> int | None:
    """Get the position of the time (or scan point) dimension."""
    for k, dim in enumerate(dims):
        if dim in TIME_DIMENSIONS:
            return k
    return None


def point_times(
    dump: Dump, sid: str, task: str, dims: Sequence[str], shape: tuple[int, ...]
) -> np.ndarray | None:
    """Get the output time of every value of a variable of a task, in its layout.

    A grid of sbmlsim 0.9 has one vector of times (`("time",)`) for every
    point of a scan; sbmlsim 0.8 and the ragged layout of sbmlsim 0.9 (`_point`)
    have the times of every point, over the dimensions of the variable, padded
    with NaN to the longest point.

    Args:
        dump: the dump.
        sid: the experiment.
        task: the task.
        dims: the dimensions of the variable.
        shape: the shape of the variable.

    Returns:
        The times in the shape of the variable, or `None` if the task has no
        times or they do not fit the variable.
    """
    entry = dump["results"].get((sid, task, "time"))
    axis = time_axis(dims)
    if entry is None or axis is None:
        return None
    tdims, values = entry
    times = values.astype(float)
    if len(tdims) == 1 and time_axis(tdims) == 0 and times.size == shape[axis]:
        expand = [1] * len(shape)
        expand[axis] = -1
        return np.broadcast_to(times.reshape(expand), shape)
    if sorted(tdims) != sorted(dims):
        return None
    times = np.transpose(times, [list(tdims).index(d) for d in dims])
    return times if times.shape == shape else None


#: the fields of the specification of a figure with the id of a task or a dataset
ID_FIELDS: tuple[str, ...] = ("task_id", "dset_id")

#: the parts of a dump whose differences can be expected, and their keys
EXPECTED_PARTS: tuple[str, ...] = (
    "results",
    "datasets",
    "mappings",
    "figures",
    "ranges",
)


class RowMatch(NamedTuple):
    """The rows of the output of a task in two dumps which are compared."""

    #: the rows of the old and of the new output which are compared
    old: np.ndarray
    new: np.ndarray
    #: rows left out at a time which the other output has fewer times (the
    #: states before a change)
    ignored_old: int
    ignored_new: int
    #: rows of the old output without a time in the new output
    unmatched: int
    #: rows of the new output without a time in the old output
    unmatched_new: int
    #: the times of the old output which the new output does not have, and
    #: the other way round
    times_old: np.ndarray
    times_new: np.ndarray


def _runs(times: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Get the first row and the number of rows of every time of an output."""
    starts = np.flatnonzero(np.append(True, np.diff(times) != 0))
    return starts, np.diff(np.append(starts, times.size))


def _last_rows(ends: np.ndarray, counts: np.ndarray) -> np.ndarray:
    """Get the last `counts[k]` rows before `ends[k]` of every time, in order."""
    offsets = np.arange(counts.sum()) - np.repeat(np.cumsum(counts) - counts, counts)
    return np.repeat(ends - counts, counts) + offsets


def matched_rows(old_times: np.ndarray, new_times: np.ndarray) -> RowMatch:
    """Match the rows of the output times of a task in two dumps.

    The times are matched relative to the largest time 1e-6. The rows of a time
    which an output has more than once (a change: the state before it, then the
    state after it) are matched from the last one, the k-th last old row with
    the k-th last new row, and the surplus rows of the output with more of them
    are left out. So the state after a change is compared across versions
    (sbmlsim 0.8 had only this one, sbmlsim 0.9 has the states before and after
    it), and two dumps of sbmlsim 0.9 compare every row of their common times.

    Args:
        old_times: the output times of the old dump.
        new_times: the output times of the new dump.

    Returns:
        The rows which are compared and the rows which are not.
    """
    if old_times.shape == new_times.shape and np.array_equal(old_times, new_times):
        rows = np.arange(old_times.size)
        none = np.zeros(0)
        return RowMatch(rows, rows, 0, 0, 0, 0, none, none)
    if old_times.size == 0 or new_times.size == 0:
        empty = np.arange(0)
        return RowMatch(
            empty,
            empty,
            0,
            0,
            int(old_times.size),
            int(new_times.size),
            np.unique(old_times),
            np.unique(new_times),
        )
    starts_a, counts_a = _runs(old_times)
    starts_b, counts_b = _runs(new_times)
    a, b = old_times[starts_a], new_times[starts_b]
    j = np.clip(np.searchsorted(b, a), 0, b.size - 1)
    previous = np.clip(j - 1, 0, b.size - 1)
    j = np.where(np.abs(b[previous] - a) < np.abs(b[j] - a), previous, j)
    common = np.abs(b[j] - a) < 1e-6 * max(1.0, float(np.max(np.abs(a))))
    # a new time is matched with one old time only, the first close to it
    first = np.zeros(a.size, dtype=bool)
    first[np.unique(np.where(common, j, -1), return_index=True)[1]] = True
    ka = np.flatnonzero(common & first)
    kb = j[ka]
    n = np.minimum(counts_a[ka], counts_b[kb])
    only_a, only_b = np.ones(a.size, dtype=bool), np.ones(b.size, dtype=bool)
    only_a[ka], only_b[kb] = False, False
    return RowMatch(
        _last_rows(starts_a[ka] + counts_a[ka], n),
        _last_rows(starts_b[kb] + counts_b[kb], n),
        int(np.sum(counts_a[ka] - n)),
        int(np.sum(counts_b[kb] - n)),
        int(old_times.size - np.sum(counts_a[ka])),
        int(new_times.size - np.sum(counts_b[kb])),
        a[only_a],
        b[only_b],
    )


def _squeezed(
    values: np.ndarray, dims: list[str], other: list[str]
) -> tuple[np.ndarray, list[str]]:
    """Leave out the dimensions of size 1 which the other array does not have."""
    drop = [
        k for k, d in enumerate(dims) if d not in other and values.shape[1 + k] == 1
    ]
    if not drop:
        return values, dims
    kept = [d for k, d in enumerate(dims) if k not in drop]
    return np.squeeze(values, axis=tuple(1 + k for k in drop)), kept


def _aligned(
    dims_a: Sequence[str], a: np.ndarray, dims_b: Sequence[str], b: np.ndarray
) -> tuple[np.ndarray, np.ndarray] | None:
    """Get two arrays with the time first and the other dimensions in the order of `a`.

    Both arrays have a time dimension. A dimension of size 1 which only one array
    has is left out (a scan of one point is a simulation).

    Returns:
        The arrays as floats, or `None` if their other dimensions differ.
    """
    axis_a, axis_b = time_axis(dims_a), time_axis(dims_b)
    assert axis_a is not None and axis_b is not None
    x = np.moveaxis(a.astype(float), axis_a, 0)
    y = np.moveaxis(b.astype(float), axis_b, 0)
    rest_a = [str(d) for k, d in enumerate(dims_a) if k != axis_a]
    rest_b = [str(d) for k, d in enumerate(dims_b) if k != axis_b]
    x, rest_a = _squeezed(x, rest_a, rest_b)
    y, rest_b = _squeezed(y, rest_b, rest_a)
    if sorted(rest_a) != sorted(rest_b):
        return None
    return x, np.transpose(y, [0, *(1 + rest_b.index(d) for d in rest_a)])


class Columns(NamedTuple):
    """A variable of a task in two dumps, one column per point of a scan.

    The time is the first dimension, the points are the other dimensions of the
    variable, matched by name and flattened in the order of the old dump.
    """

    #: the values of the old and of the new dump, (rows, points)
    old: np.ndarray
    new: np.ndarray
    #: the output time of every value
    old_times: np.ndarray
    new_times: np.ndarray


def columns(old: Dump, new: Dump, key: tuple[str, str, str]) -> Columns | None:
    """Get a variable of a task in two dumps with the times of every value.

    Returns:
        The columns, or `None` if the dimensions of the variables differ or a
        dump has no times of the task which fit its variable.
    """
    (dims_a, a), (dims_b, b) = old["results"][key], new["results"][key]
    sid, task, _ = key
    ta = point_times(old, sid, task, dims_a, a.shape)
    tb = point_times(new, sid, task, dims_b, b.shape)
    if ta is None or tb is None:
        return None
    values = _aligned(dims_a, a, dims_b, b)
    times = _aligned(dims_a, ta, dims_b, tb)
    if values is None or times is None or values[0].shape[1:] != values[1].shape[1:]:
        return None
    x, y, tx, ty = (v.reshape(v.shape[0], -1) for v in (*values, *times))
    return Columns(x, y, tx, ty)


def task_columns(old: Dump, new: Dump, sid: str, task: str) -> Columns | None:
    """Get the output times of every point of a task in two dumps.

    The times are those of a variable of the task (the first by name), so
    that the points are matched by the names of their dimensions as the
    variables are; a task with the times only compares their arrays.

    Returns:
        The times as the columns (values and times), or `None`.
    """
    variables = sorted(
        v
        for (s, t, v) in old["results"]
        if s == sid and t == task and v != "time" and (s, t, v) in new["results"]
    )
    for variable in variables:
        found = columns(old, new, (sid, task, variable))
        if found is not None:
            return Columns(found.old_times, found.new_times, *found[2:])
    key = (sid, task, "time")
    if key in old["results"] and key in new["results"]:
        return columns(old, new, key)
    return None


class PointMatch(NamedTuple):
    """The rows of points of a task which are compared, see `point_matches`."""

    #: the points (columns) with these times
    points: np.ndarray
    #: the rows of the old and of the new values which are compared
    old: np.ndarray
    new: np.ndarray
    #: the match of the times of the points
    match: RowMatch


def _same_columns(times: np.ndarray) -> bool:
    """Check that every point has the times of the first, NaN equal to NaN."""
    first = times[:, :1]
    return bool(np.all((times == first) | (np.isnan(times) & np.isnan(first))))


def point_matches(old_times: np.ndarray, new_times: np.ndarray) -> list[PointMatch]:
    """Match the rows of every point of a task by its own times.

    The points of a grid have the same times and one match; a point of a
    ragged result is matched by its own times, without the NaN which pad it to
    the longest point. Points with the same times in both dumps share a match.

    Args:
        old_times: the times of the old dump, (rows, points).
        new_times: the times of the new dump, (rows, points).

    Returns:
        The matches.
    """
    points = old_times.shape[1]
    if _same_columns(old_times) and _same_columns(new_times):
        groups = [list(range(points))]
    else:
        by_times: dict[bytes, list[int]] = {}
        for p in range(points):
            key = old_times[:, p].tobytes() + new_times[:, p].tobytes()
            by_times.setdefault(key, []).append(p)
        groups = list(by_times.values())
    matches: list[PointMatch] = []
    for group in groups:
        ta, tb = old_times[:, group[0]], new_times[:, group[0]]
        rows_a, rows_b = np.flatnonzero(~np.isnan(ta)), np.flatnonzero(~np.isnan(tb))
        match = matched_rows(ta[rows_a], tb[rows_b])
        matches.append(
            PointMatch(np.array(group), rows_a[match.old], rows_b[match.new], match)
        )
    return matches


def relative_difference(old: Dump, new: Dump, key: tuple[str, str, str]) -> float:
    """Get the difference of a result array on the matched rows.

    The time is moved first, the other dimensions are matched by name (a
    dimension of size 1 which only one array has is left out), and the rows of
    every point are matched by its own times (`point_matches`). The difference
    is max|new - old| / max|old|; other dimensions and NaN at different
    positions are an infinite difference. A numeric array without a time, e.g.
    the values at the end of the points of a scan, is compared on the same
    dimensions relative to its scale, any other array exactly; the variable
    `time` is compared as the output times of the points of the task.
    """
    (dims_a, a), (dims_b, b) = old["results"][key], new["results"][key]
    sid, task, variable = key
    if variable == "time":
        # sbmlsim 0.8 repeated the time of a scan for every point, (_point, *dims);
        # a grid of sbmlsim 0.9 has it once, ("time",): compare the output times
        found = task_columns(old, new, sid, task)
    else:
        numeric = a.dtype.kind in "fiub" and b.dtype.kind in "fiub"
        timeless = time_axis(dims_a) is None and time_axis(dims_b) is None
        if numeric and timeless:
            if tuple(dims_a) != tuple(dims_b) or a.shape != b.shape:
                return math.inf
            return scaled_difference([(a, b)])
        if not numeric or time_axis(dims_a) is None or time_axis(dims_b) is None:
            same = a.shape == b.shape and np.array_equal(a, b)
            return 0.0 if same else math.inf
        found = columns(old, new, key)
    if found is None:
        return math.inf
    return scaled_difference(
        (
            found.old[np.ix_(match.old, match.points)],
            found.new[np.ix_(match.new, match.points)],
        )
        for match in point_matches(found.old_times, found.new_times)
    )


def scaled_difference(pairs: Iterable[tuple[np.ndarray, np.ndarray]]) -> float:
    """Get max|new - old| / max|old| over pairs of old and new arrays.

    NaN at different positions is an infinite difference.
    """
    scale = difference = 0.0
    for x, y in pairs:
        if not np.array_equal(np.isnan(x), np.isnan(y)):
            return math.inf
        x, y = np.nan_to_num(x), np.nan_to_num(y)
        if x.size:
            scale = max(scale, float(np.max(np.abs(x))))
            difference = max(difference, float(np.max(np.abs(x - y))))
    return difference / scale if scale > 1e-300 else difference


def equal(a: object, b: object) -> bool:
    """Check that two values (arrays, lists, dicts) are equal, NaN equal to NaN."""
    if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        x, y = np.asarray(a), np.asarray(b)
        if x.shape != y.shape or x.dtype.kind != y.dtype.kind:
            return False
        return bool(np.array_equal(x, y, equal_nan=x.dtype.kind == "f"))
    return a == b


class MetadataRenames(TypedDict, total=False):
    """The expected changes of the metadata of the fit mappings."""

    #: old field name -> new field name
    fields: dict[str, str]
    #: field (new name) -> old member name -> new member name
    members: dict[str, dict[str, str]]
    #: field -> member name of a field the old metadata did not have
    added: dict[str, str]


class Renames(TypedDict, total=False):
    """The expected renames and differences between two dumps."""

    metadata: MetadataRenames
    #: old text -> new text, in every key and string of the old dump
    strings: dict[str, str]
    #: old id -> new id, of a task, a dataset, a figure or a fit mapping
    ids: dict[str, str]
    #: every id of the old dump becomes its SId (`pkpdsim.experiment.to_sid`)
    sids: bool
    #: the keys of the specification of a figure which are left out of both
    #: dumps, at the level of the figure, after `align_figures`: a field of a
    #: figure which `compare` cannot map
    figure_fields: list[str]
    #: part of the dump -> patterns of the keys whose differences are expected
    #: (`*` and `?` are wildcards, a bracket is literal, e.g. `E/task_po/[A]`):
    #: `experiment/task/variable` (results; `experiment/task/time` also names
    #: the output times which only the new dump has),
    #: `experiment/dataset/column` (datasets), `experiment/figure` (figures),
    #: `experiment/mapping` (mappings, ranges)
    expected: dict[str, list[str]]


#: the parts of a renames file
RENAME_PARTS: tuple[str, ...] = (
    "metadata",
    "strings",
    "ids",
    "sids",
    "figure_fields",
    "expected",
)


def load_renames(path: str) -> Renames:
    """Load the expected renames from a json file.

    Raises:
        ValueError: for a part which is not known.
    """
    renames = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(renames, dict):
        raise ValueError(f"the renames {path} must be a json object")
    unknown = sorted(set(renames) - set(RENAME_PARTS)) + sorted(
        set(renames.get("metadata", {})) - {"fields", "members", "added"}
    )
    if unknown:
        raise ValueError(f"unknown parts of the renames {path}: {unknown}")
    _check_renames(path, renames)
    return cast("Renames", renames)


def _is_strings(value: object) -> bool:
    """Check for a json object of strings."""
    return isinstance(value, dict) and all(
        isinstance(k, str) and isinstance(v, str) for k, v in value.items()
    )


def _is_texts(value: object) -> bool:
    """Check for a json list of strings."""
    return isinstance(value, list) and all(isinstance(v, str) for v in value)


def _check_renames(path: str, renames: dict[str, object]) -> None:
    """Check the types of the parts of the renames.

    Raises:
        ValueError: naming the part which has the wrong type.
    """
    if "strings" in renames and not _is_strings(renames["strings"]):
        raise ValueError(
            f"the part strings of the renames {path} must map text to text"
        )
    if "ids" in renames and not _is_strings(renames["ids"]):
        raise ValueError(f"the part ids of the renames {path} must map text to text")
    if "sids" in renames and not isinstance(renames["sids"], bool):
        raise ValueError(f"the part sids of the renames {path} must be true or false")
    if "figure_fields" in renames and not _is_texts(renames["figure_fields"]):
        raise ValueError(
            f"the part figure_fields of the renames {path} must be a list of texts"
        )
    if "expected" in renames:
        expected = renames["expected"]
        if (
            not isinstance(expected, dict)
            or set(expected) - set(EXPECTED_PARTS)
            or not all(_is_texts(v) for v in expected.values())
        ):
            raise ValueError(
                f"the part expected of the renames {path} must map parts of "
                f"{list(EXPECTED_PARTS)} to lists of patterns"
            )
    if "metadata" not in renames:
        return
    metadata = renames["metadata"]
    if not isinstance(metadata, dict):
        raise ValueError(f"the part metadata of the renames {path} must be an object")
    for part in ("fields", "added"):
        if part in metadata and not _is_strings(metadata[part]):
            raise ValueError(
                f"the part {part} of the metadata of the renames {path} "
                "must map text to text"
            )
    members = metadata.get("members", {})
    if not isinstance(members, dict) or not all(
        _is_strings(v) for v in members.values()
    ):
        raise ValueError(
            f"the part members of the metadata of the renames {path} must map "
            "a field to an object of text to text"
        )


def rename_mapping(description: str, renames: MetadataRenames) -> str:
    """Apply the expected renames to the metadata of a mapping description."""
    head, metadata = description.rsplit("|", 1)
    fields = renames.get("fields", {})
    items = {fields.get(k, k): v for k, v in ast.literal_eval(metadata)}
    for name, members in renames.get("members", {}).items():
        if name in items:
            items[name] = members.get(items[name], items[name])
    for name, member in renames.get("added", {}).items():
        items.setdefault(name, member)
    return f"{head}|{sorted(items.items())}"


def _rename_strings(value: object, strings: Mapping[str, str]) -> object:
    """Replace texts in every key and string of a value, arrays untouched.

    Raises:
        ValueError: if two keys of a dict have the same new key.
    """
    if isinstance(value, str):
        for old, new in strings.items():
            value = value.replace(old, new)
        return value
    if isinstance(value, dict):
        renamed: dict[object, object] = {}
        for k, v in value.items():
            key = _rename_strings(k, strings)
            if key in renamed:
                raise ValueError(
                    f"Two keys of the old dump have the same new key {key!r} "
                    "after the renames of strings."
                )
            renamed[key] = _rename_strings(v, strings)
        return renamed
    if isinstance(value, list):
        return [_rename_strings(v, strings) for v in value]
    if isinstance(value, tuple):
        return tuple(_rename_strings(v, strings) for v in value)
    return value


def _rename_description(description: str, rename: Callable[[str], str]) -> str:
    """Rename the key, the datasets and the tasks of a mapping description.

    Raises:
        ValueError: if the head of the description is not its 16 fields, i.e. a
            key has a '|'.
    """
    head, metadata = description.rsplit("|", 1)
    fields = head.split("|")
    if len(fields) != 16:
        raise ValueError(
            f"The mapping description cannot be renamed, a key has a '|': {description}"
        )
    for k in (0, 1, 2, 8, 9):
        fields[k] = rename(fields[k])
    return "|".join([*fields, metadata])


def _rename_figure(spec: str, rename: Callable[[str], str]) -> str:
    """Rename the id of a figure and the ids of its tasks and datasets in its specification."""

    def walk(value: object, top: bool) -> object:
        if isinstance(value, dict):
            return {
                k: rename(v)
                if isinstance(v, str) and (k in ID_FIELDS or (top and k == "sid"))
                else walk(v, False)
                for k, v in value.items()
            }
        if isinstance(value, list):
            return [walk(v, False) for v in value]
        return value

    return json.dumps(walk(json.loads(spec), True), sort_keys=True)


def _unique[K, V](renamed: dict[K, V], before: Mapping[K, V], part: str) -> dict[K, V]:
    """Check that no two keys of a part got the same new key.

    Raises:
        ValueError: if two keys have the same new id.
    """
    if len(renamed) != len(before):
        raise ValueError(f"Two ids of the {part} of the old dump have the same new id.")
    return renamed


def rename_ids(dump: Dump, rename: Callable[[str], str]) -> Dump:
    """Rename the ids of the tasks, datasets, figures and fit mappings of a dump.

    Returns:
        A new dump; the given one is not changed.

    Raises:
        ValueError: if two ids of a part have the same new id, or a mapping
            description cannot be renamed.
    """
    result = cast("Dump", dict(dump))
    result["results"] = _unique(
        {(s, rename(t), v): x for (s, t, v), x in dump["results"].items()},
        dump["results"],
        "results",
    )
    result["datasets"] = _unique(
        {(s, rename(d), c): x for (s, d, c), x in dump["datasets"].items()},
        dump["datasets"],
        "datasets",
    )
    result["mappings"] = {
        sid: [_rename_description(d, rename) for d in descriptions]
        for sid, descriptions in dump["mappings"].items()
    }
    result["ranges"] = {
        sid: _unique({rename(k): v for k, v in items.items()}, items, "ranges")
        for sid, items in dump["ranges"].items()
    }
    result["figures"] = {
        sid: _unique(
            {rename(f): _rename_figure(spec, rename) for f, spec in figures.items()},
            figures,
            "figures",
        )
        for sid, figures in dump["figures"].items()
    }
    return result


def drop_figure_fields(dump: Dump, fields: Collection[str]) -> Dump:
    """Leave out keys of the specifications of the figures, of the figure only.

    A key is left out at the level of the figure, not inside its plots, curves
    or axes: a plot has a `width` and a `height` of its own.

    Returns:
        A new dump; the given one is not changed.
    """
    names = set(fields)

    def without(spec: str) -> str:
        figure = json.loads(spec)
        if isinstance(figure, dict):
            figure = {k: v for k, v in figure.items() if k not in names}
        return json.dumps(figure, sort_keys=True)

    result = cast("Dump", dict(dump))
    result["figures"] = {
        sid: {f: without(spec) for f, spec in figures.items()}
        for sid, figures in dump["figures"].items()
    }
    return result


def _is_figure_of_08(spec: object) -> bool:
    """Check whether a specification is one of a figure of sbmlsim 0.8.

    A figure of sbmlsim 0.8 always has a size (`_width`, `_height`), one of
    sbmlsim 0.9 has its settings (`_settings`).
    """
    return isinstance(spec, dict) and "_width" in spec and "_settings" not in spec


def _is_number(value: object) -> bool:
    """Check for a json number."""
    return isinstance(value, int | float) and not isinstance(value, bool)


def _without_derived_plot_sizes(figure: dict[str, Any]) -> None:
    """Leave out the size of a plot which its figure gave it, in place.

    A plot which is added to a figure of a size gets the width of the figure
    divided by its columns times the columns of the plot, likewise the height:
    in sbmlsim 0.8 always, a figure had a size, in sbmlsim 0.9 only for a
    figure of an explicit size. The size of the figure is compared, so a plot
    size which equals the one of the figure is `None`; a size of its own stays.
    """
    for subplot in figure.get("subplots") or []:
        plot = subplot.get("plot") if isinstance(subplot, dict) else None
        if not isinstance(plot, dict):
            continue
        for size, count, span in (
            ("width", "num_cols", "col_span"),
            ("height", "num_rows", "row_span"),
        ):
            total, n, k = figure.get(size), figure.get(count), subplot.get(span)
            value = plot.get(size)
            if (
                all(_is_number(x) for x in (total, n, k, value))
                and n
                and math.isclose(value, total / n * k, rel_tol=1e-12)
            ):
                plot[size] = None


def _figure_of_08(spec: dict[str, Any]) -> dict[str, Any]:
    """Bring the specification of a figure of sbmlsim 0.8 into the compared form.

    The settings the figure set as attributes of its own go into `_settings`,
    the size (`_width`, `_height`, which sbmlsim 0.8 computed from the panels
    when the figure was created) into `width` and `height`.
    """
    figure = dict(spec)
    figure["_settings"] = {
        name: figure.pop(name) for name in FIGURE_SETTINGS_08 if name in figure
    }
    figure["width"] = figure.pop("_width")
    figure["height"] = figure.pop("_height", None)
    _without_derived_plot_sizes(figure)
    return figure


def _without_new_fields(
    value: object, path: str, listed: list[tuple[str, object]]
) -> object:
    """Leave out the fields of the objects which only sbmlsim 0.9 has.

    A field which is not empty is added to `listed` with its path.
    """
    if isinstance(value, dict):
        result: dict[str, object] = {}
        for k, v in value.items():
            if "__type__" in value and k in NEW_FIGURE_FIELDS:
                if v:
                    listed.append((f"{path}.{k}".lstrip("."), v))
                continue
            result[k] = _without_new_fields(v, f"{path}.{k}", listed)
        return result
    if isinstance(value, list):
        return [
            _without_new_fields(v, f"{path}[{k}]", listed) for k, v in enumerate(value)
        ]
    return value


def _figure_of_09(
    spec: dict[str, Any], listed: list[tuple[str, object]]
) -> dict[str, Any]:
    """Bring the specification of a figure of sbmlsim 0.9 into the compared form.

    `_settings`, `width` and `height` are those the figure is rendered with
    without the settings of the run (`RESOLVED`): the settings of its
    experiment and its own, the explicit size or the one of its panels. The
    fields which only sbmlsim 0.9 has (`NEW_FIGURE_FIELDS`) are left out and
    the ones which are not empty added to `listed`.
    """
    figure = cast("dict[str, Any]", _without_new_fields(spec, "", listed))
    resolved = figure.pop(RESOLVED, None)
    if isinstance(resolved, dict):
        figure["_settings"] = resolved.get("settings")
        figure["width"] = resolved.get("width")
        figure["height"] = resolved.get("height")
    _without_derived_plot_sizes(figure)
    return figure


class AlignedFigures(NamedTuple):
    """The figures of two dumps in the form in which they are compared."""

    old: Dump
    new: Dump
    #: the number of figures of sbmlsim 0.8 compared with one of sbmlsim 0.9
    across: int
    #: experiment, figure, path and value of a field which only sbmlsim 0.9 has
    new_fields: list[tuple[str, str, str, object]]


def align_figures(old: Dump, new: Dump) -> AlignedFigures:
    """Bring a figure of sbmlsim 0.8 and the one of sbmlsim 0.9 into one form.

    sbmlsim 0.9 moved the settings of a figure from its attributes
    (`legend_position`) into `_settings` and its size from `_width` and
    `_height`, which sbmlsim 0.8 computed when the figure was created, into
    `width` and `height`, which are `None` for the size of the panels at the
    time it is rendered; a plot gets a size from its figure only for an
    explicit one. A figure which is in both dumps, of sbmlsim 0.8 in the old and
    of sbmlsim 0.9 in the new one, is compared in the form of sbmlsim 0.9, with
    the settings and the size it is rendered with apart from the settings of
    the run (`RESOLVED` of the dump, see `figure_spec`), and without the size
    of a plot which its figure gave it. The selection of a `Data`, the
    dimensions of the lines of a `Curve` and the bands of a `Plot`
    (`NEW_FIGURE_FIELDS`) have no counterpart in sbmlsim 0.8: they are left
    out and listed. Two figures of one version are compared as they are.

    Returns:
        The two dumps (new ones, the given ones are not changed), the number of
        figures compared across versions and the listed fields.
    """
    figures_old = {sid: dict(figures) for sid, figures in old["figures"].items()}
    figures_new = {sid: dict(figures) for sid, figures in new["figures"].items()}
    across = 0
    new_fields: list[tuple[str, str, str, object]] = []
    for sid in sorted(figures_old.keys() & figures_new.keys()):
        a_figures, b_figures = figures_old[sid], figures_new[sid]
        for fid in sorted(a_figures.keys() & b_figures.keys()):
            a, b = json.loads(a_figures[fid]), json.loads(b_figures[fid])
            if (
                not _is_figure_of_08(a)
                or not isinstance(b, dict)
                or "_settings" not in b
            ):
                continue
            listed: list[tuple[str, object]] = []
            a_figures[fid] = json.dumps(_figure_of_08(a), sort_keys=True)
            b_figures[fid] = json.dumps(_figure_of_09(b, listed), sort_keys=True)
            across += 1
            new_fields += [(sid, fid, path, value) for path, value in listed]
    return AlignedFigures(
        old={**old, "figures": figures_old},
        new={**new, "figures": figures_new},
        across=across,
        new_fields=new_fields,
    )


def _renamer(renames: Renames) -> Callable[[str], str]:
    """Get the rename of an id: `ids`, then its SId with `sids`."""
    from pkpdsim.experiment import to_sid

    ids = renames.get("ids", {})
    sids = renames.get("sids", False)

    def rename(text: str) -> str:
        text = ids.get(text, text)
        return to_sid(text) if sids else text

    return rename


def apply_renames(dump: Dump, renames: Renames) -> Dump:
    """Apply the expected renames to a dump, the old one of a comparison.

    The strings first, then the ids (`ids`, `sids`), then the metadata.

    Returns:
        A new dump; the given one is not changed.

    Raises:
        ValueError: if two keys or ids get the same new key or id, or a mapping
            description cannot be renamed.
    """
    result = cast("Dump", dict(dump))
    if strings := renames.get("strings"):
        # the texts of the parts, not the names of the parts
        result = cast(
            "Dump", {part: _rename_strings(v, strings) for part, v in result.items()}
        )
    if renames.get("ids") or renames.get("sids"):
        result = rename_ids(result, _renamer(renames))
    if metadata := renames.get("metadata"):
        result["mappings"] = {
            sid: [rename_mapping(d, metadata) for d in descriptions]
            for sid, descriptions in result["mappings"].items()
        }
    # a dump stores the mappings sorted, and a rename can change the order
    result["mappings"] = {
        sid: sorted(descriptions) for sid, descriptions in result["mappings"].items()
    }
    return result


class _Expected:
    """The expected differences of a renames file and which of them occurred."""

    def __init__(self, patterns: Mapping[str, Sequence[str]]) -> None:
        self.patterns = {part: list(p) for part, p in patterns.items()}
        self.used: set[tuple[str, str]] = set()

    def match(self, part: str, key: str | tuple[str, ...]) -> bool:
        """Check whether the difference of a key is expected, recording the patterns.

        `*` and `?` of a pattern are the wildcards of `fnmatch`; a bracket is
        literal, so that `E/task/[A]` is the variable `[A]`.
        """
        text = "/".join(key) if isinstance(key, tuple) else key
        found = False
        for pattern in self.patterns.get(part, []):
            if fnmatch.fnmatchcase(text, pattern.replace("[", "[[]")):
                self.used.add((part, pattern))
                found = True
        return found

    def unused(self) -> list[tuple[str, str]]:
        """Get the patterns which matched no difference."""
        return [
            (part, pattern)
            for part, patterns in self.patterns.items()
            for pattern in patterns
            if (part, pattern) not in self.used
        ]


def _flat(
    dump: Dump, part: Literal["datasets", "mappings", "figures", "ranges"]
) -> dict[str | tuple[str, ...], object]:
    """Get a part of a dump by the keys a comparison reports.

    A fit mapping is keyed by its experiment and its key, the first field of its
    description.

    Raises:
        ValueError: if two fit mappings of an experiment have one key (a key with
            a '|').
    """
    if part == "mappings":
        mappings: dict[str | tuple[str, ...], object] = {}
        for sid, descriptions in dump["mappings"].items():
            for description in descriptions:
                key = (sid, description.split("|", 1)[0])
                if key in mappings:
                    raise ValueError(
                        f"Two fit mappings of {sid} have the key {key[1]}, a key "
                        "has a '|'."
                    )
                mappings[key] = description
        return mappings
    if part in ("figures", "ranges"):
        return {
            (sid, key): value
            for sid, items in dump[part].items()
            for key, value in items.items()
        }
    return dict(dump[part].items())


def _listed(times: np.ndarray, count: int = 5) -> str:
    """List the first times of a vector."""
    listed = ", ".join(f"{t:g}" for t in times[:count])
    return listed + (", ..." if times.size > count else "")


def _detail(part: str, a: object, b: object) -> list[str]:
    """Describe the difference of a changed key."""
    if part == "mappings":
        return diff_lines([a] if a else [], [b] if b else [])
    if part == "figures":
        return json_paths(_json(a), _json(b))[:5]
    return [f"{a!r} -> {b!r}"]


def compare(args: argparse.Namespace) -> int:
    """Compare two dumps; the exit code is 1 if they differ beyond the threshold.

    An expected difference (`expected` of the renames) is printed and not
    counted; an expected difference which does not occur is counted. A row at a
    time which only one dump has is a difference: one of the new dump is
    expected with the pattern of the variable `time` of its task
    (`experiment/task/time`, a changed output grid), one of the old dump always
    counts.
    """
    old, new = load(args.old), load(args.new)
    renames: Renames = {}
    if args.renames:
        try:
            renames = load_renames(args.renames)
            old = apply_renames(old, renames)
        except ValueError as err:
            print(err, file=sys.stderr)
            return 2
    aligned = align_figures(old, new)
    old, new = aligned.old, aligned.new
    if fields := renames.get("figure_fields"):
        old, new = drop_figure_fields(old, fields), drop_figure_fields(new, fields)
    try:
        parts = {
            part: (_flat(old, part), _flat(new, part))
            for part in ("datasets", "mappings", "figures", "ranges")
        }
    except ValueError as err:
        print(err, file=sys.stderr)
        return 2
    expected = _Expected(renames.get("expected", {}))
    for label, dump_ in (("old", old), ("new", new)):
        if versions := dump_.get("versions"):
            print(
                f"{label}: "
                + ", ".join(f"{k} {v}" for k, v in sorted(versions.items()))
            )
    names: set[str] = set(args.experiments)

    def selected(sid: str) -> bool:
        return not names or sid in names

    differing = 0
    for label, dump_ in (("old", old), ("new", new)):
        for sid in sorted(filter(selected, dump_["errors"])):
            print(f"error in {label}: {sid}")
            differing += 1
    only_old = [s for s in old["order"] if s not in new["order"] and selected(s)]
    only_new = [s for s in new["order"] if s not in old["order"] and selected(s)]
    for sid in only_old:
        print(f"only in old: {sid}")
    for sid in only_new:
        print(f"only in new: {sid}")
    differing += len(only_old) + len(only_new)

    # simulations
    keys_old = {k for k in old["results"] if selected(k[0])}
    keys_new = {k for k in new["results"] if selected(k[0])}
    one_side = sorted(keys_old ^ keys_new)
    unexpected_side = [k for k in one_side if not expected.match("results", k)]
    for key in one_side[: args.top]:
        mark = "" if key in unexpected_side else "expected: "
        print(
            f"{mark}array only in {'old' if key in keys_old else 'new'}: {'/'.join(key)}"
        )
    rows = sorted(
        ((relative_difference(old, new, key), key) for key in keys_old & keys_new),
        reverse=True,
    )
    per_experiment: dict[str, float] = {}
    for difference, key in rows:
        per_experiment[key[0]] = max(per_experiment.get(key[0], 0.0), difference)
    above = [row for row in rows if row[0] > args.threshold]
    unexpected = [row for row in above if not expected.match("results", row[1])]
    print(
        f"simulations: {len(rows)} arrays compared, {len(above)} differ by more "
        f"than {args.threshold:g} ({len(above) - len(unexpected)} expected), "
        f"{len(one_side)} in one dump only"
    )
    for sid, difference in sorted(per_experiment.items(), key=lambda t: -t[1]):
        if difference > 0:
            print(f"  {difference:.3e} {sid}")
    if rows and rows[0][0] > 0:
        print("largest differences:")
        for difference, key in rows[: args.top]:
            if difference > 0:
                print(f"  {difference:.3e} {'/'.join(key)}")
    differing += len(unexpected) + len(unexpected_side)

    # the rows at a change
    ignored_old = ignored_new = unmatched = unmatched_new = expected_new = 0
    changes: dict[str, int] = {}
    other_times: list[str] = []
    for sid, task in sorted({(k[0], k[1]) for k in keys_old & keys_new}):
        times = task_columns(old, new, sid, task)
        matches = [] if times is None else point_matches(times.old, times.new)
        for match in (m.match for m in matches):
            ignored_old += match.ignored_old
            ignored_new += match.ignored_new
            unmatched += match.unmatched
            unmatched_new += match.unmatched_new
            if match.ignored_new:
                changes[sid] = changes.get(sid, 0) + match.ignored_new
            if match.unmatched:
                other_times.append(
                    f"  times of {sid}/{task}: {match.unmatched} rows only in old, "
                    f"at {_listed(match.times_old)}"
                )
            if match.unmatched_new:
                # a changed output grid of a task is named by its variable `time`
                named = expected.match("results", (sid, task, "time"))
                expected_new += match.unmatched_new if named else 0
                other_times.append(
                    f"  times of {sid}/{task}: {match.unmatched_new} rows only in new, "
                    f"at {_listed(match.times_new)}{' (expected)' if named else ''}"
                )
    print(
        f"rows at a change: {ignored_new} rows of the new dump and {ignored_old} of "
        f"the old dump left out (the states before a change), {unmatched} rows of "
        f"the old dump without a time in the new dump, {unmatched_new} rows of the "
        f"new dump without a time in the old dump ({expected_new} expected)"
    )
    for sid, count in sorted(changes.items()):
        print(f"changes: {sid} {count} rows")
    for line in other_times[: args.top]:
        print(line)
    # a time which only one dump has is a difference; one of the new dump can be
    # expected, the old dump is the reference
    differing += unmatched + unmatched_new - expected_new

    # the figures of sbmlsim 0.8 compared with those of sbmlsim 0.9
    if aligned.across:
        print(
            f"figures across versions of sbmlsim: {aligned.across} (the settings "
            "and the size as rendered without the settings of the run; "
            f"{', '.join(NEW_FIGURE_FIELDS)} of the new dump listed, not compared)"
        )
    for sid, fid, path, value in aligned.new_fields:
        if selected(sid):
            print(f"new fields: {sid}/{fid}: {path} = {json.dumps(value)}")

    # everything else exactly
    for part, (flat_old, flat_new) in parts.items():
        a = {k: v for k, v in flat_old.items() if selected(_sid(k))}
        b = {k: v for k, v in flat_new.items() if selected(_sid(k))}
        changed = sorted(
            (k for k in a.keys() | b.keys() if not equal(a.get(k), b.get(k))), key=str
        )
        unexpected_keys = [k for k in changed if not expected.match(part, k)]
        print(
            f"{part}: {len(a)} old, {len(b)} new, {len(changed)} differ "
            f"({len(changed) - len(unexpected_keys)} expected)"
        )
        for k in changed[: args.top]:
            mark = "" if k in unexpected_keys else " (expected)"
            print(f"  {'/'.join(k) if isinstance(k, tuple) else k}{mark}")
            for line in _detail(part, a.get(k), b.get(k))[: args.top]:
                print(f"    {line}")
        differing += len(unexpected_keys)

    for part, pattern in expected.unused():
        print(f"expected difference which does not occur: {part}: {pattern}")
        differing += 1
    return 1 if differing else 0


def _sid(key: str | tuple[str, ...]) -> str:
    """Get the experiment of a key of a dump."""
    return key[0] if isinstance(key, tuple) else key


def diff_lines(a: object, b: object) -> list[str]:
    """Describe the difference of two lists or dicts line by line.

    The values of a dict are compared exactly, json values (the specification of
    a figure) by the paths of their elements which differ.
    """
    if isinstance(a, dict) and isinstance(b, dict):
        lines: list[str] = []
        for k in sorted(a.keys() | b.keys(), key=str):
            if k not in a or k not in b:
                lines.append(f"{k}: only {'old' if k in a else 'new'}")
            elif not equal(a[k], b[k]):
                paths = json_paths(_json(a[k]), _json(b[k]))
                lines.append(f"{k}: {', '.join(paths[:5]) or 'changed'}")
        return lines
    x = list(a) if isinstance(a, list) else []
    y = list(b) if isinstance(b, list) else []
    return [f"- {v}" for v in x if v not in y] + [f"+ {v}" for v in y if v not in x]


def _json(value: object) -> object:
    """Parse a json string, any other value is returned as it is."""
    if isinstance(value, str):
        try:
            return json.loads(value)
        except json.JSONDecodeError:
            return value
    return value


def json_paths(a: object, b: object, path: str = "") -> list[str]:
    """Get the paths of the elements in which two json values differ.

    A key which only one object has is named with its value, also a `null` one
    (a field which one version of sbmlsim added to a figure).
    """
    if isinstance(a, dict) and isinstance(b, dict):
        paths: list[str] = []
        for k in sorted(a.keys() | b.keys()):
            if k in a and k in b:
                paths += json_paths(a[k], b[k], f"{path}.{k}")
            else:
                side, value = ("old", a[k]) if k in a else ("new", b[k])
                where = f"{path}.{k}".lstrip(".")
                paths.append(f"{where} only in {side}: {value!r}")
        return paths
    if isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        return [
            p
            for k, (x, y) in enumerate(zip(a, b, strict=True))
            for p in json_paths(x, y, f"{path}[{k}]")
        ]
    if a != b:
        return [f"{path.lstrip('.')} {a!r} -> {b!r}"]
    return []


# --- ranges ---------------------------------------------------------------------


def ranges(args: argparse.Namespace) -> int:
    """List the fit mappings whose data lie outside of the simulated range."""
    data = load(args.dump)
    outside = 0
    for sid, mappings in sorted(data["ranges"].items()):
        for key, value in sorted(mappings.items()):
            if isinstance(value, str):
                print(f"{sid} {key}: {value}")
                outside += 1
                continue
            dmin, dmax, smin, smax, unit = value
            tolerance = 1e-9 * max(1.0, abs(smin), abs(smax))
            ok = dmin >= smin - tolerance and dmax <= smax + tolerance
            if not ok or args.verbose:
                status = "ok" if ok else "OUTSIDE"
                print(
                    f"{status} {sid} {key}: data [{dmin:.4g}, {dmax:.4g}], "
                    f"simulation [{smin:.4g}, {smax:.4g}] {unit}"
                )
            outside += not ok
    print(f"{outside} fit mappings with data outside of the simulation")
    return 1 if outside else 0


# --- fits -----------------------------------------------------------------------


def fits(args: argparse.Namespace) -> int:
    """Run a minimal fit (one start, a few evaluations) of every fit definition."""
    root = Path(args.root).resolve()
    sys.path.insert(0, str(root))
    from sbmlsim.fit.cli import run_fit

    if message := _outside_root(args.model, root):
        print(message, file=sys.stderr)
        return 2
    info = _describe(args)
    if isinstance(info, int):
        return info
    if info.fit_definitions is None:
        print(f"{args.model} has no fit definitions", file=sys.stderr)
        return 2
    definitions = info.fit_definitions
    status: dict[str, str] = {}
    for key in args.subsets or list(definitions):
        try:
            run_fit(
                definitions[key],
                opid=f"smoke_{key}",
                size=1,
                seed=1234,
                output_dir=Path(args.out),
                max_nfev=args.max_nfev,
            )
            status[key] = "ok"
        except Exception:  # every failing subset is reported
            traceback.print_exc()
            status[key] = "FAILED"
    print(status)
    return 1 if "FAILED" in status.values() else 0


def definitions(args: argparse.Namespace) -> int:
    """Write the fit definitions of a model as json.

    The settings, the parameters and the mappings of the collections of
    every fit definition, so that the definitions of two checkouts are
    compared with `diff`. The checkout directory is replaced by `<root>`.
    """
    root_path = Path(args.root).resolve()
    root = str(root_path)
    sys.path.insert(0, root)
    if message := _outside_root(args.model, root_path):
        print(message, file=sys.stderr)
        return 2
    info = _describe(args)
    if isinstance(info, int):
        return info
    if info.fit_definitions is None:
        print(f"{args.model} has no fit definitions", file=sys.stderr)
        return 2
    content: dict[str, object] = {}
    for key, definition in info.fit_definitions.items():
        collections = definition.mapping_collections()
        content[key] = {
            "settings": serialize(definition.settings),
            "parameters": [serialize(p) for p in definition.parameters],
            "collections": {
                experiment: [serialize(c) for c in items]
                for experiment, items in sorted(collections.items())
            },
        }
    text = json.dumps(content, indent=2, sort_keys=True, default=str)
    Path(args.out).write_text(
        text.replace(json.dumps(root)[1:-1], "<root>").replace(root, "<root>"),
        encoding="utf-8",
    )
    return 0


# --- settings -------------------------------------------------------------------


def _plain(value: object) -> object:
    """A json value of a setting: paths as posix strings, classes by name."""
    if isinstance(value, Path):
        return value.as_posix()
    if isinstance(value, type):
        return value.__name__
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, list | tuple):
        return [_plain(item) for item in value]
    return serialize(value)


def _bound(
    original: Callable[..., object],
    args: tuple[object, ...],
    kwargs: dict[str, object],
    skip: Collection[str] = (),
) -> dict[str, object]:
    """The arguments of a call bound to the signature of `original`, defaults applied.

    A positional and a keyword call of the same values give the same result.
    `self` and the names in `skip` are left out.
    """
    signature = inspect.signature(original)
    bound = signature.bind(*args, **kwargs)
    bound.apply_defaults()
    arguments: dict[str, object] = {}
    for name, value in bound.arguments.items():
        if signature.parameters[name].kind is inspect.Parameter.VAR_KEYWORD:
            arguments.update(value)
        else:
            arguments[name] = value
    return {
        name: _plain(value)
        for name, value in sorted(arguments.items())
        if name != "self" and name not in skip
    }


def capture_settings(
    run_experiments: Callable[..., object],
    experiment: type["SimulationExperiment"],
    root: Path,
) -> dict[str, Any]:
    """Capture the settings a run of an experiment hands to sbmlsim.

    The run is intercepted where it reaches sbmlsim: the arguments of the
    simulator (`Simulator`), of the runner (without the simulator) and of the
    run (without the output path; with the figure settings of the run), bound
    to the signatures, so positional and keyword calls agree, and every
    matplotlib `rcParams` entry which differs from the value before the run, at
    the time of the run. `calls` counts the intercepted calls. The simulator is
    the one created before the runner: the designs of an experiment
    (`sbmlsim.simulation.sampling.parameters_of`, `local`) create simulators of
    their own while the runner initializes it, which are neither recorded nor
    counted. Nothing is simulated and no report is written.

    Args:
        run_experiments: the run of a model, `MODEL.run_experiments`.
        experiment: the experiment class the run gets.
        root: the checkout; it is replaced by `<root>` in the paths.

    Returns:
        The settings as json values.

    Raises:
        RuntimeError: if the run does not reach the simulator, the runner and
            the run of the experiments of sbmlsim.
    """
    from unittest import mock

    import matplotlib
    from sbmlsim.experiment import ExperimentRunner
    from sbmlsim.report.experiment_report import ExperimentReport
    from sbmlsim.simulator import Simulator

    captured: dict[str, Any] = {"experiment": experiment.__name__}
    calls = dict.fromkeys(("simulator", "runner", "run"), 0)
    simulator_init = Simulator.__init__
    runner_init = ExperimentRunner.__init__
    run_original = ExperimentRunner.run_experiments
    baseline = dict(matplotlib.rcParams)

    def record_simulator(*args: object, **kwargs: object) -> None:
        # a design of an experiment creates simulators of its own while the
        # runner initializes it
        if calls["runner"] == 0:
            calls["simulator"] += 1
            captured["simulator"] = _bound(simulator_init, args, kwargs)
        simulator_init(*args, **kwargs)  # ty: ignore[invalid-argument-type]

    def record_runner(*args: object, **kwargs: object) -> None:
        calls["runner"] += 1
        # the simulator is captured on its own; as an object it carries memory
        # addresses, which differ between any two runs
        captured["runner"] = _bound(runner_init, args, kwargs, skip={"simulator"})
        runner_init(*args, **kwargs)  # ty: ignore[invalid-argument-type]

    def record_run(*args: object, **kwargs: object) -> list[object]:
        calls["run"] += 1
        captured["run"] = _bound(run_original, args, kwargs, skip={"output_path"})
        captured["rc"] = {
            key: _plain(value)
            for key, value in sorted(dict(matplotlib.rcParams).items())
            if key not in baseline or bool(value != baseline[key])
        }
        return []

    def no_report(
        self: object, output_path: Path, *args: object, **kwargs: object
    ) -> Path:
        return output_path

    with (
        tempfile.TemporaryDirectory() as tmp,
        mock.patch.object(Simulator, "__init__", record_simulator),
        mock.patch.object(ExperimentRunner, "__init__", record_runner),
        mock.patch.object(ExperimentRunner, "run_experiments", record_run),
        mock.patch.object(ExperimentReport, "create_report", no_report),
    ):
        run_experiments(experiment, output_dir=Path(tmp) / "settings")
    missing = [
        key for key in ("simulator", "runner", "run", "rc") if key not in captured
    ]
    if missing:
        raise RuntimeError(
            f"the run did not reach sbmlsim, nothing captured for {missing}"
        )
    captured["calls"] = calls
    text = json.dumps(captured, sort_keys=True, default=str)
    for path in (root.resolve(), root.absolute()):
        for form in (path.as_posix(), str(path), json.dumps(str(path))[1:-1]):
            text = text.replace(form, "<root>")
    return json.loads(text)


def settings(args: argparse.Namespace) -> int:
    """Write the settings of a run of the experiments of a model as json."""
    root = Path(args.root).resolve()
    sys.path.insert(0, str(root))
    if message := _outside_root(args.model, root):
        print(message, file=sys.stderr)
        return 2
    info = _describe(args)
    if isinstance(info, int):
        return info
    if info.run_experiments is None:
        print(f"{args.model} has no run of its experiments", file=sys.stderr)
        return 2
    classes = registered_experiments(info) or []
    if args.experiment:
        experiment = {c.__name__: c for c in classes}.get(args.experiment)
    else:
        experiment = next(iter(classes), None)
    if experiment is None:
        print(
            f"no experiment {args.experiment or ''} in the registry of {args.model}",
            file=sys.stderr,
        )
        return 2
    try:
        captured = capture_settings(
            info.run_experiments, experiment, Path(args.root).absolute()
        )
    except RuntimeError as err:
        print(f"{args.model}: {err}", file=sys.stderr)
        return 2
    Path(args.out).write_text(
        json.dumps(captured, indent=2, sort_keys=True), encoding="utf-8"
    )
    return 0


# --- png ------------------------------------------------------------------------


def _files(directory: Path, pattern: str) -> dict[Path, Path]:
    """Get the files below a directory by their path relative to it."""
    return {p.relative_to(directory): p for p in sorted(directory.rglob(pattern))}


def _directories(old: Path, new: Path) -> str | None:
    """Get the message about a directory which does not exist, or `None`."""
    for directory in (old, new):
        if not directory.is_dir():
            return f"{directory} is not a directory"
    return None


def png(args: argparse.Namespace) -> int:
    """Compare the PNG figures of two directories pixel by pixel.

    The directories are searched recursively and the figures are paired by the
    path relative to the directory; with `--sids` an old figure is paired with
    the new figure of the SId of its file name (`pkpdsim.experiment.to_sid`).

    Returns:
        0 if all figures are identical, 1 if some differ or exist on one side
        only, 2 if a directory does not exist, no figure was found or two old
        figures have the same SId.
    """
    from PIL import Image

    old, new = Path(args.old), Path(args.new)
    if message := _directories(old, new):
        print(message, file=sys.stderr)
        return 2

    def selected(path: Path) -> bool:
        return not args.names or any(s in path.name for s in args.names)

    old_files = {k: v for k, v in _files(old, "*.png").items() if selected(v)}
    new_files = {k: v for k, v in _files(new, "*.png").items() if selected(v)}
    if args.sids:
        from pkpdsim.experiment import to_sid

        by_sid: dict[Path, Path] = {}
        for key, path in old_files.items():
            renamed = key.with_name(to_sid(key.stem) + key.suffix)
            if renamed in by_sid:
                print(f"two old figures have the same SId: {renamed}", file=sys.stderr)
                return 2
            by_sid[renamed] = path
        old_files = by_sid
    if not old_files and not new_files:
        print(f"no PNG figures in {old} and {new}", file=sys.stderr)
        return 2
    same = differing = 0
    for key in sorted(old_files.keys() | new_files.keys()):
        if key not in new_files:
            print(f"only in old: {key}")
            differing += 1
            continue
        if key not in old_files:
            print(f"only in new: {key}")
            differing += 1
            continue
        a = np.asarray(Image.open(old_files[key]))
        b = np.asarray(Image.open(new_files[key]))
        if a.shape == b.shape and np.array_equal(a, b):
            same += 1
            continue
        fraction = (a != b).any(axis=-1).mean() if a.shape == b.shape else 1.0
        print(f"differs: {key} {a.shape} {b.shape}, pixels {fraction:.2e}")
        differing += 1
    print(f"{same} identical, {differing} differ")
    return 1 if differing else 0


# --- sbml -----------------------------------------------------------------------


def sbml_elements(path: Path) -> dict[str, object] | None:
    """Get the elements of an SBML model which determine its simulation.

    Returns:
        The elements by kind and id, `None` if the file has no SBML model (the
        layouts the factories write next to the models are no SBML).
    """

    def formula(math_: libsbml.ASTNode | None) -> str | None:
        return libsbml.formulaToL3String(math_) if math_ is not None else None

    def references(items: Iterable[libsbml.SpeciesReference]) -> tuple:
        return tuple(sorted((r.getSpecies(), r.getStoichiometry()) for r in items))

    # the document owns the model, it has to stay referenced while the model is read
    document = libsbml.readSBMLFromFile(str(path))
    model = document.getModel()
    if model is None:
        return None
    elements: dict[str, object] = {}
    for p in model.getListOfParameters():
        # an unset value is NaN, as a value of NaN
        value = p.getValue()
        elements[f"parameter {p.getId()}"] = (value, p.getConstant(), p.getUnits())
    for s in model.getListOfSpecies():
        elements[f"species {s.getId()}"] = (
            s.getInitialAmount() if s.isSetInitialAmount() else None,
            s.getInitialConcentration() if s.isSetInitialConcentration() else None,
            s.getCompartment(),
            s.getHasOnlySubstanceUnits(),
            s.getBoundaryCondition(),
            s.getConstant(),
            s.getSubstanceUnits(),
        )
    for c in model.getListOfCompartments():
        size = c.getSize()
        elements[f"compartment {c.getId()}"] = (size, c.getConstant(), c.getUnits())
    for r in model.getListOfRules():
        elements[f"rule {r.getVariable()} {r.getTypeCode()}"] = formula(r.getMath())
    for a in model.getListOfInitialAssignments():
        elements[f"initial assignment {a.getSymbol()}"] = formula(a.getMath())
    for r in model.getListOfReactions():
        law = r.getKineticLaw()
        elements[f"reaction {r.getId()}"] = (
            formula(law.getMath()) if law else None,
            references(r.getListOfReactants()),
            references(r.getListOfProducts()),
            tuple(sorted(m.getSpecies() for m in r.getListOfModifiers())),
        )
    for e in model.getListOfEvents():
        elements[f"event {e.getId()}"] = (
            formula(e.getTrigger().getMath()),
            tuple(
                (a.getVariable(), formula(a.getMath()))
                for a in e.getListOfEventAssignments()
            ),
        )
    for u in model.getListOfUnitDefinitions():
        elements[f"unit {u.getId()}"] = tuple(
            sorted(
                (x.getKind(), x.getExponent(), x.getScale(), x.getMultiplier())
                for x in u.getListOfUnits()
            )
        )
    return {k: _nan_equal(v) for k, v in elements.items()}


def _nan_equal(value: object) -> object:
    """Replace NaN by a string, so that NaN equals NaN."""
    if isinstance(value, tuple):
        return tuple(_nan_equal(v) for v in value)
    if isinstance(value, float) and math.isnan(value):
        return "nan"
    return value


def sbml(args: argparse.Namespace) -> int:
    """Compare SBML models: two files, or the files of the same path in two dirs.

    The directories are searched recursively and the files are paired by the
    path relative to the directory.

    Returns:
        0 if the models are equal, 1 if they differ or exist on one side only,
        2 if a directory does not exist or no file was found.
    """
    old, new = Path(args.old), Path(args.new)
    if old.is_dir() or new.is_dir():
        if message := _directories(old, new):
            print(message, file=sys.stderr)
            return 2
        names = sorted(_files(old, "*.xml").keys() | _files(new, "*.xml").keys())
        if not names:
            print(f"no XML files in {old} and {new}", file=sys.stderr)
            return 2
        pairs = [(old / name, new / name) for name in names]
    else:
        pairs = [(old, new)]
    differing = 0
    for a_path, b_path in pairs:
        label = str(a_path.relative_to(old)) if old.is_dir() else a_path.name
        if not a_path.exists() or not b_path.exists():
            print(f"{label}: only in {'new' if b_path.exists() else 'old'}")
            differing += 1
            continue
        a, b = sbml_elements(a_path), sbml_elements(b_path)
        if a is None or b is None:
            print(f"{label}: no SBML model")
            continue
        changed = [k for k in sorted(a.keys() | b.keys()) if a.get(k) != b.get(k)]
        print(f"{label}: {len(changed)} of {len(a.keys() | b.keys())} differ")
        for k in changed:
            print(f"  {k}\n    old: {a.get(k)}\n    new: {b.get(k)}")
        differing += len(changed)
    return 1 if differing else 0


# --- command line ---------------------------------------------------------------


def parser() -> argparse.ArgumentParser:
    """Create the parser of the command line."""
    p = argparse.ArgumentParser(
        prog="pkpdsim-verify", description=__doc__.split("\n\n")[0] if __doc__ else None
    )
    sub = p.add_subparsers(required=True)

    def command(
        name: str, function: Callable[[argparse.Namespace], int], help_: str
    ) -> argparse.ArgumentParser:
        s = sub.add_parser(name, help=help_, description=help_)
        s.set_defaults(function=function)
        return s

    s = command("dump", dump, "run the experiments of a model and store everything")
    s.add_argument("--root", required=True, help=ROOT_HELP)
    s.add_argument("--model", required=True, help=MODEL_HELP)
    s.add_argument("--out", required=True, help="pickle the dump is written to")
    s.add_argument(
        "--all",
        action="store_true",
        help="also the experiments of experiments/ which are not in the registry",
    )
    s.add_argument("--order", help="run in the order of a previous dump")
    s.add_argument(
        "--skip", default=[], type=lambda s: s.split(","), help="A,B: leave these out"
    )
    s.add_argument(
        "--atol", help="absolute tolerance, or `model` for that of the model"
    )
    s.add_argument(
        "--rtol", help="relative tolerance, or `model` for that of the model"
    )
    s.add_argument(
        "--reduced-selections",
        action="store_true",
        help="only the selections the experiments need, as a run of the model does",
    )
    s.add_argument("experiments", nargs="*", help="only these experiment classes")

    s = command("compare", compare, "compare two dumps")
    s.add_argument("old")
    s.add_argument("new")
    s.add_argument(
        "--renames",
        help=(
            "json file of the expected renames and differences (metadata, strings, "
            "ids, sids, figure_fields, expected), applied to the old dump; in a "
            "pattern of expected, * and ? are wildcards and a bracket is literal, "
            "experiment/task/time names the output times only the new dump has"
        ),
    )
    s.add_argument("--threshold", type=float, default=1e-6, help="default 1e-6")
    s.add_argument("--top", type=int, default=20, help="lines per listing")
    s.add_argument("experiments", nargs="*", help="only these experiments")

    s = command("ranges", ranges, "fit mappings with data outside of the simulation")
    s.add_argument("dump")
    s.add_argument("--verbose", action="store_true", help="every mapping")

    s = command("fits", fits, "a minimal fit of every fit definition of a model")
    s.add_argument("--root", required=True, help=ROOT_HELP)
    s.add_argument("--model", required=True, help=MODEL_HELP)
    s.add_argument("--out", required=True, help="output directory of the fits")
    s.add_argument("--max-nfev", type=int, default=3, help="evaluations per start")
    s.add_argument("subsets", nargs="*", help="only these fit definitions")

    s = command("definitions", definitions, "the fit definitions of a model as json")
    s.add_argument("--root", required=True, help=ROOT_HELP)
    s.add_argument("--model", required=True, help=MODEL_HELP)
    s.add_argument(
        "--out", required=True, help="json file the definitions are written to"
    )

    s = command(
        "settings", settings, "the settings a run of a model hands to sbmlsim as json"
    )
    s.add_argument("--root", required=True, help=ROOT_HELP)
    s.add_argument("--model", required=True, help=MODEL_HELP)
    s.add_argument("--out", required=True, help="json file the settings are written to")
    s.add_argument(
        "--experiment",
        help="the experiment class to run, the first of the registry by default",
    )

    s = command("png", png, "compare the PNG figures of two directories")
    s.add_argument("old")
    s.add_argument("new")
    s.add_argument(
        "--sids",
        action="store_true",
        help=(
            "pair an old figure with the new figure of the SId of its file name "
            "(pkpdsim.experiment.to_sid), for a port which renamed the keys"
        ),
    )
    s.add_argument("names", nargs="*", help="only files containing one of these")

    s = command("sbml", sbml, "compare two SBML files or the SBML files of two dirs")
    s.add_argument("old")
    s.add_argument("new")
    return p


def verify(args: Sequence[str] | None = None) -> int:
    """Run a subcommand.

    Args:
        args: the command line arguments, `sys.argv[1:]` by default.

    Returns:
        The exit code, 1 if the subcommand found a difference.
    """
    os.environ.setdefault("MPLBACKEND", "Agg")
    options = parser().parse_args(args)
    return options.function(options)


def main() -> None:
    """The command `pkpdsim-verify`."""
    sys.exit(verify())


if __name__ == "__main__":
    # the workers of a fit start with forkserver/spawn and import this module
    main()
