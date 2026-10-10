"""Verify that a change of a model changes only what it is meant to change.

A refactoring of a model, an update of sbmlsim or a refresh of the study data
should leave the simulations of the model as they were, and a fix should change
them only where it is meant to. The subcommands make this checkable:

- `dump` runs the simulation experiments of a model in a checkout and stores the
  results, the datasets, the fit mappings and the figures in a pickle,
- `compare` compares two dumps: the simulations on their common output times,
  relative to the scale of every array, everything else exactly,
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
the dotted name of a model package. A package with a `MODEL`
(`pkpdsim.model.PKPDModel`) is described by it; a package without is described
by the conventions of pkdb_models before pkpdsim 0.3.0: `<NAME>_PATH` (else the
directory of the package), `MODEL_PATH`, `DATA_PATHS` (or `DATA_PATH`), the
registry `simulations.EXPERIMENTS`, the tolerance `helpers.ABSOLUTE_TOLERANCE`
and the fit definitions `fitting.fitting.FIT_DEFINITIONS`, so that a checkout
from before the switch is the base of a comparison. The command only needs the
python packages of the environment it runs in, so a checkout which needs another
version of sbmlsim is dumped with the python of an environment with that
version, as long as it is one which pkpdsim supports.
"""

import argparse
import importlib
import inspect
import json
import math
import os
import pickle
import pkgutil
import sys
import tempfile
import traceback
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import TYPE_CHECKING, Any, TypedDict

import libsbml
import numpy as np

if TYPE_CHECKING:
    from sbmlsim.experiment import SimulationExperiment
    from sbmlsim.fit import FitData, FitMapping
    from sbmlsim.fit.cli import FitDefinition

ROOT_HELP = (
    "directory put first on the python path, e.g. a checkout of the repository "
    "of the model"
)
MODEL_HELP = "dotted name of the model package, e.g. pkdb_models.models.hctz"

#: dimensions along which the simulations of a task are compared, the time of a
#: timecourse and the points of a scan
TIME_DIMENSIONS = ("_time", "_point")

#: experiment groups of the registry which are run first, in this order
GROUPS = ("studies", "misc", "scans")


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


def _attribute(module_name: str, attribute: str) -> object | None:
    """Get an attribute of a module, `None` if the module or a parent does not exist."""
    try:
        module = importlib.import_module(module_name)
    except ModuleNotFoundError as err:
        if err.name != module_name and not module_name.startswith(f"{err.name}."):
            # the module exists but does not import
            raise
        return None
    return getattr(module, attribute, None)


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
    """Describe a model package by its `MODEL` or by the conventions.

    Args:
        package: the dotted name of the model package.

    Returns:
        What the verification needs of the model.
    """
    from pkpdsim.model import PKPDModel

    module = importlib.import_module(package)
    model = getattr(module, "MODEL", None)
    if isinstance(model, PKPDModel):
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
    name = package.rsplit(".", 1)[-1]
    base_path = (
        getattr(module, f"{name.upper()}_PATH", None)
        or Path(str(module.__file__)).parent
    )
    registry = _attribute(f"{package}.simulations", "EXPERIMENTS")
    definitions = _attribute(f"{package}.fitting.fitting", "FIT_DEFINITIONS")
    run = _attribute(f"{package}.helpers", "run_experiments")
    return ModelInfo(
        package=package,
        base_path=Path(base_path),
        model_path=Path(module.MODEL_PATH),
        data_path=getattr(module, "DATA_PATHS", None)
        or getattr(module, "DATA_PATH", None),
        experiments=registry if isinstance(registry, dict) else None,
        absolute_tolerance=_attribute(f"{package}.helpers", "ABSOLUTE_TOLERANCE"),
        relative_tolerance=_attribute(f"{package}.helpers", "RELATIVE_TOLERANCE"),
        run_experiments=run if callable(run) else None,
        fit_definitions=definitions if isinstance(definitions, dict) else None,
    )


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
        sd = d.y_sd.index if d.y_sd else None
        se = d.y_se.index if d.y_se else None
        return f"{d.dset_id}|{d.task_id}|{d.x.index}|{d.y.index}|{d.count}|{sd}|{se}"

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
    from sbmlsim.simulator.simulation_serial import SimulatorSerial

    if message := _outside_root(args.model, root):
        print(message, file=sys.stderr)
        return 2
    info = describe(args.model)
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

    tolerances: dict[str, object] = {}
    if args.atol == "model":
        if info.absolute_tolerance is None:
            print(f"{args.model} has no tolerance of its own", file=sys.stderr)
            return 2
        tolerances["absolute_tolerance"] = info.absolute_tolerance
    elif args.atol is not None:
        tolerances["absolute_tolerance"] = float(args.atol)
    if args.rtol == "model":
        if info.relative_tolerance is None:
            print(f"{args.model} has no relative tolerance of its own", file=sys.stderr)
            return 2
        tolerances["relative_tolerance"] = info.relative_tolerance
    elif args.rtol is not None:
        tolerances["relative_tolerance"] = float(args.rtol)

    result = Dump(
        model=args.model,
        order=[c.__name__ for c in classes],
        errors=errors,
        results={},
        datasets={},
        mappings={},
        figures={},
        ranges={},
    )
    for cls in classes:
        try:
            runner = ExperimentRunner(
                experiment_classes=[cls],
                data_path=info.data_path,
                base_path=info.base_path,
                simulator=SimulatorSerial(model=info.model_path, **tolerances),
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
    """Store the results, datasets, fit mappings and figures of a run experiment."""
    for task, xres in experiment.results.items():
        # `xdataset` up to sbmlsim 0.8.3
        ds = getattr(xres, "xds", None)
        if ds is None:
            ds = getattr(xres, "xdataset")  # noqa: B009
        for var in ds.data_vars:
            result["results"][sid, task, str(var)] = (
                tuple(str(d) for d in ds[var].dims),
                np.asarray(ds[var].values),
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
        fid: json.dumps(serialize(fig), sort_keys=True, default=repr)
        for fid, fig in experiment.figures().items()
    }


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


def time_values(dump: Dump, sid: str, task: str) -> np.ndarray | None:
    """Get the output times of a task as a vector."""
    entry = dump["results"].get((sid, task, "time"))
    if entry is None:
        return None
    dims, values = entry
    axis = time_axis(dims)
    if axis is None:
        return None
    values = np.moveaxis(values.astype(float), axis, 0)
    return values.reshape(values.shape[0], -1)[:, 0]


def relative_difference(old: Dump, new: Dump, key: tuple[str, str, str]) -> float:
    """Get the difference of a result array on the common output times.

    At a repeated output time (the end of a timecourse and the start of the next)
    the last value is used. The difference is max|new - old| / max|old| on the
    common times; NaN at different positions is an infinite difference.
    """
    (dims_a, a), (dims_b, b) = old["results"][key], new["results"][key]
    axis_a, axis_b = time_axis(dims_a), time_axis(dims_b)
    if a.dtype.kind not in "fiub" or axis_a is None or axis_b is None:
        same = a.shape == b.shape and np.array_equal(a, b)
        return 0.0 if same else math.inf
    sid, task, _ = key
    ta, tb = time_values(old, sid, task), time_values(new, sid, task)
    if ta is None or tb is None:
        return math.inf
    a = np.moveaxis(a.astype(float), axis_a, 0)
    b = np.moveaxis(b.astype(float), axis_b, 0)
    keep = np.append(np.diff(ta) != 0, True)
    a, ta = a[keep], ta[keep]
    j = np.clip(np.searchsorted(tb, ta), 0, len(tb) - 1)
    previous = np.clip(j - 1, 0, len(tb) - 1)
    j = np.where(np.abs(tb[previous] - ta) < np.abs(tb[j] - ta), previous, j)
    common = np.abs(tb[j] - ta) < 1e-6 * max(1.0, float(np.max(np.abs(ta))))
    x, y = a[common], b[j[common]]
    if x.shape[1:] != y.shape[1:]:
        return math.inf
    if not np.array_equal(np.isnan(x), np.isnan(y)):
        return math.inf
    x, y = np.nan_to_num(x), np.nan_to_num(y)
    if x.size == 0:
        return 0.0
    scale = float(np.max(np.abs(x)))
    difference = float(np.max(np.abs(x - y)))
    return difference / scale if scale > 1e-300 else difference


def equal(a: object, b: object) -> bool:
    """Check that two values (arrays, lists, dicts) are equal, NaN equal to NaN."""
    if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        x, y = np.asarray(a), np.asarray(b)
        if x.shape != y.shape or x.dtype.kind != y.dtype.kind:
            return False
        return bool(np.array_equal(x, y, equal_nan=x.dtype.kind == "f"))
    return a == b


def compare(args: argparse.Namespace) -> int:
    """Compare two dumps; the exit code is 1 if they differ beyond the threshold."""
    old, new = load(args.old), load(args.new)
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
    for key in sorted(keys_old ^ keys_new)[: args.top]:
        print(f"array only in {'old' if key in keys_old else 'new'}: {'/'.join(key)}")
    rows = sorted(
        ((relative_difference(old, new, key), key) for key in keys_old & keys_new),
        reverse=True,
    )
    per_experiment: dict[str, float] = {}
    for difference, key in rows:
        per_experiment[key[0]] = max(per_experiment.get(key[0], 0.0), difference)
    above = [row for row in rows if row[0] > args.threshold]
    print(
        f"simulations: {len(rows)} arrays compared, {len(above)} differ by more "
        f"than {args.threshold:g}, {len(keys_old ^ keys_new)} in one dump only"
    )
    # the experiments and arrays which differ at all, the largest first
    for sid, difference in sorted(per_experiment.items(), key=lambda t: -t[1]):
        if difference > 0:
            print(f"  {difference:.3e} {sid}")
    if rows and rows[0][0] > 0:
        print("largest differences:")
        for difference, key in rows[: args.top]:
            if difference > 0:
                print(f"  {difference:.3e} {'/'.join(key)}")
    differing += len(above) + len(keys_old ^ keys_new)

    # everything else exactly
    for part in ("datasets", "mappings", "figures", "ranges"):
        a = {k: v for k, v in old[part].items() if selected(_sid(k))}
        b = {k: v for k, v in new[part].items() if selected(_sid(k))}
        changed = sorted(
            (k for k in a.keys() | b.keys() if not equal(a.get(k), b.get(k))), key=str
        )
        print(f"{part}: {len(a)} old, {len(b)} new, {len(changed)} differ")
        for k in changed[: args.top]:
            print(f"  {'/'.join(k) if isinstance(k, tuple) else k}")
            if isinstance(a.get(k), list | dict) or isinstance(b.get(k), list | dict):
                for line in diff_lines(a.get(k), b.get(k))[: args.top]:
                    print(f"    {line}")
        differing += len(changed)
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
    """Get the paths of the elements in which two json values differ."""
    if isinstance(a, dict) and isinstance(b, dict):
        return [
            p
            for k in sorted(a.keys() | b.keys())
            for p in json_paths(a.get(k), b.get(k), f"{path}.{k}")
        ]
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
    info = describe(args.model)
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
    info = describe(args.model)
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
    if isinstance(value, list | tuple):
        return [_plain(item) for item in value]
    return serialize(value)


def capture_settings(
    run_experiments: Callable[..., object],
    experiment: type["SimulationExperiment"],
    root: Path,
) -> dict[str, Any]:
    """Capture the settings a run of an experiment hands to sbmlsim.

    The run is intercepted where it reaches sbmlsim: the arguments of the
    simulator, of the runner and of the run of the experiments, and the
    settings of the figures and `savefig.dpi` of matplotlib at the time of the
    run. Nothing is simulated and no report is written. The figure settings of
    sbmlsim are restored afterwards, also those a legacy run leaves set.

    Args:
        run_experiments: the run of a model, `MODEL.run_experiments` or the
            `run_experiments` of a legacy `helpers` module.
        experiment: the experiment class the run gets.
        root: the checkout; it is replaced by `<root>` in the paths.

    Returns:
        The settings as json values.
    """
    from unittest import mock

    import matplotlib
    from sbmlsim.experiment import ExperimentRunner
    from sbmlsim.plot import Figure
    from sbmlsim.report.experiment_report import ExperimentReport
    from sbmlsim.simulator.simulation_serial import SimulatorSerial

    from pkpdsim.model import FIGURE_ATTRIBUTES
    from pkpdsim.testing import restored_figure_settings

    captured: dict[str, object] = {"experiment": experiment.__name__}
    simulator_init = SimulatorSerial.__init__
    runner_init = ExperimentRunner.__init__

    def record_simulator(self: object, *args: object, **kwargs: object) -> None:
        captured["simulator"] = {k: _plain(v) for k, v in sorted(kwargs.items())}
        simulator_init(self, *args, **kwargs)  # ty: ignore[invalid-argument-type]

    def record_runner(self: object, *args: object, **kwargs: object) -> None:
        # the simulator is captured on its own; as an object it carries memory
        # addresses, which differ between any two runs
        captured["runner"] = {
            k: _plain(v) for k, v in sorted(kwargs.items()) if k != "simulator"
        }
        runner_init(self, *args, **kwargs)  # ty: ignore[invalid-argument-type]

    def record_run(self: object, *args: object, **kwargs: object) -> list[object]:
        captured["run"] = {
            k: _plain(v) for k, v in sorted(kwargs.items()) if k != "output_path"
        }
        captured["figure"] = {
            name: _plain(getattr(Figure, name)) for name in FIGURE_ATTRIBUTES
        }
        captured["savefig_dpi"] = _plain(matplotlib.rcParams["savefig.dpi"])
        return []

    def no_report(
        self: object, output_path: Path, *args: object, **kwargs: object
    ) -> Path:
        return output_path

    with (
        tempfile.TemporaryDirectory() as tmp,
        restored_figure_settings(),
        mock.patch.object(SimulatorSerial, "__init__", record_simulator),
        mock.patch.object(ExperimentRunner, "__init__", record_runner),
        mock.patch.object(ExperimentRunner, "run_experiments", record_run),
        mock.patch.object(ExperimentReport, "create_report", no_report),
    ):
        run_experiments(experiment, output_dir=Path(tmp) / "settings")
    text = json.dumps(captured, sort_keys=True, default=str)
    resolved = root.resolve().as_posix()
    return json.loads(text.replace(resolved, "<root>"))


def settings(args: argparse.Namespace) -> int:
    """Write the settings of a run of the experiments of a model as json."""
    root = Path(args.root).resolve()
    sys.path.insert(0, str(root))
    if message := _outside_root(args.model, root):
        print(message, file=sys.stderr)
        return 2
    info = describe(args.model)
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
    captured = capture_settings(info.run_experiments, experiment, root)
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
    path relative to the directory.

    Returns:
        0 if all figures are identical, 1 if some differ or exist on one side
        only, 2 if a directory does not exist or no figure was found.
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
