"""The pipeline of a model: its commands in order, timed and logged.

The analysis of a model is a sequence of the commands of pkpdsim: the factory,
the simulation experiments, their bands over the draws of the fitted
parameters, the sensitivity analysis, the assessment of the credibility, the
refits of the adopted fits. Which of them a model runs, in which order and with
which options differs between models, so a model package declares its pipeline
(`Pipeline`, usually `PIPELINE` in its `pipeline.py`) and its model object
references it (`pipeline=".pipeline:PIPELINE"`); `run_pipeline <package>` runs
it (`pkpdsim.cli.pipeline_command`):

```python
PIPELINE = Pipeline(
    description="The complete analysis of the model from its adopted fits.",
    steps=(
        Step("factory", "run_simulation", ("--action", "factory")),
        Step("simulate", "run_simulation", ("--action", "simulate", "--experiments", "all", "--workers", "{cores}")),
        Step("sensitivity", "run_sensitivity", ("--workers", "{cores}")),
        Step("credibility", "run_credibility", ("--uncertainty",)),
    ),
)
```

Every step runs its command in a process of its own
(`python -m pkpdsim.cli <command> <package> ...`), one after the other, so a
step which fails or runs out of memory leaves the runner and the results of the
steps before it intact, and every step has the machine to itself: a command
parallelizes internally (the runs of a fit, the scans of a profile, the
experiments of a simulation, the points of a sensitivity design), and two of
them side by side would compete for the same cores. `{cores}` in the arguments
of a step is the number of processes the pipeline gives every step
(`--cores`, the CPUs available to the process by default) and `{results}` the
results of the model, so a pipeline states where it parallelizes and the
machine decides how wide.

A run writes into its own folder below `pipeline/` of the results of the model:
the output of every step in `logs/<step>.log`, the times in `timings.tsv`
(wall clock, CPU time of the step and of every process it started, the cores
this is on average, the peak memory of its processes) and the versions and the
machine in `pipeline.json`, and shows the steps with their times and the last
line of progress of the running one on the console.
"""

from __future__ import annotations

import datetime
import json
import os
import platform
import re
import subprocess
import sys
import time
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from rich.live import Live
from rich.table import Table
from sbmlsim.cgroup import available_cpus
from sbmlutils.console import console

__all__ = [
    "PLACEHOLDERS",
    "Pipeline",
    "Step",
    "StepResult",
    "run_pipeline",
]

#: the placeholders of the arguments of a step: the processes of a step and the
#: results of the model
PLACEHOLDERS = ("cores", "results")

#: the packages whose versions a run records
VERSIONS = (
    "pkpdsim",
    "sbmlsim",
    "sbmlsim-core",
    "sbmlutils",
    "sbmlode",
    "pymetadata",
    "pkpdutils",
    "libroadrunner",
    "numpy",
    "scipy",
)

#: the name of a step: a word of letters, digits and underscores
STEP_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")

#: the seconds between two refreshes of the console and samples of the memory
REFRESH = 1.0

#: the characters of the line of progress of a running step on the console
PROGRESS_WIDTH = 72

#: the source location rich appends to a line of the log of a step
_LOCATION = re.compile(r"\s+[\w.]+\.py:\d+\s*$")

#: the level rich writes before a line of the log of a step
_LEVEL = re.compile(r"^(DEBUG|INFO)\s+")


@dataclass(frozen=True)
class Step:
    """A step of a pipeline: a command of pkpdsim on the model of the package.

    Attributes:
        name: the name of the step, unique in its pipeline, a word.
        command: the command, a console script of pkpdsim, e.g. `run_fit`.
        args: the arguments of the command after the package; `{cores}` is the
            processes of a step, `{results}` the results of the model.
        description: what the step does, listed by `run_pipeline --list`.
    """

    name: str
    command: str
    args: tuple[str, ...] = ()
    description: str = ""

    def __post_init__(self) -> None:
        """Check the name and the placeholders of the arguments.

        Raises:
            ValueError: for a name which is no word, arguments which are a
                string, or a placeholder which is none of `PLACEHOLDERS`.
        """
        if not STEP_NAME.match(self.name):
            raise ValueError(f"The name of the step '{self.name}' is no word.")
        if isinstance(self.args, str):
            raise ValueError(
                f"The arguments of the step '{self.name}' are a string, not a tuple."
            )
        for arg in self.args:
            for placeholder in re.findall(r"\{(\w*)\}", arg):
                if placeholder not in PLACEHOLDERS:
                    raise ValueError(
                        f"The argument '{arg}' of the step '{self.name}' has the "
                        f"placeholder '{{{placeholder}}}', none of {PLACEHOLDERS}."
                    )

    def arguments(self, values: Mapping[str, object]) -> list[str]:
        """Get the arguments with the placeholders replaced by their values."""
        return [
            re.sub(r"\{(\w+)\}", lambda m: str(values[m.group(1)]), arg)
            for arg in self.args
        ]


@dataclass(frozen=True)
class Pipeline:
    """The steps of the analysis of a model, run one after the other.

    Attributes:
        steps: the steps in the order they run, at least one.
        description: what the pipeline does.
    """

    steps: tuple[Step, ...]
    description: str = ""

    def __post_init__(self) -> None:
        """Check the steps.

        Raises:
            ValueError: for no step or two steps of one name.
        """
        if not self.steps:
            raise ValueError("A pipeline has at least one step.")
        names = [step.name for step in self.steps]
        twice = sorted({name for name in names if names.count(name) > 1})
        if twice:
            raise ValueError(f"The pipeline has two steps of the names {twice}.")

    @property
    def names(self) -> list[str]:
        """The names of the steps in their order."""
        return [step.name for step in self.steps]

    def select(
        self,
        only: Sequence[str] = (),
        skip: Sequence[str] = (),
        start: str | None = None,
    ) -> list[Step]:
        """Get the steps of a run, in the order of the pipeline.

        Args:
            only: the names of the steps which run, every step without.
            skip: the names of the steps which do not run.
            start: the name of the first step which runs, the steps before it
                do not.

        Raises:
            ValueError: for a name which is no step of the pipeline.
        """
        for name in [*only, *skip, *([start] if start else [])]:
            if name not in self.names:
                raise ValueError(
                    f"'{name}' is no step of the pipeline, the steps are "
                    f"{', '.join(self.names)}."
                )
        steps = list(self.steps)
        if start is not None:
            steps = steps[self.names.index(start) :]
        if only:
            steps = [step for step in steps if step.name in only]
        return [step for step in steps if step.name not in skip]


@dataclass
class StepResult:
    """What a step of a run did.

    Attributes:
        step: the step.
        argv: the command line of its process.
        log: the file of its output.
        status: `pending`, `running`, `ok`, `failed` or `skipped`.
        exit_code: the exit code of its process, `None` before it ended.
        wall: the seconds it ran.
        cpu: the CPU seconds of its process and of every process it started.
        peak_rss: the peak of the resident bytes of its processes, sampled.
        progress: the last line of progress of its output.
    """

    step: Step
    argv: list[str]
    log: Path
    status: str = "pending"
    exit_code: int | None = None
    wall: float = 0.0
    cpu: float = 0.0
    peak_rss: int = 0
    progress: str = ""
    _started: float = field(default=0.0, repr=False)

    @property
    def cores(self) -> float:
        """The cores the step used on average, its CPU time over its wall time."""
        return self.cpu / self.wall if self.wall > 0 else 0.0

    def row(self) -> dict[str, object]:
        """Get the row of the step in `timings.tsv`."""
        return {
            "step": self.step.name,
            "command": self.step.command,
            "status": self.status,
            "exit_code": self.exit_code,
            "wall_s": round(self.wall, 1),
            "cpu_s": round(self.cpu, 1),
            "cores": round(self.cores, 2),
            "peak_rss_gb": round(self.peak_rss / 1024**3, 2),
            "log": self.log.name,
        }


def _duration(seconds: float) -> str:
    """Get a duration as `1 h 02 m`, `3 m 05 s` or `12.3 s`."""
    if seconds >= 3600:
        hours, rest = divmod(round(seconds), 3600)
        return f"{hours} h {rest // 60:02d} m"
    if seconds >= 60:
        minutes, rest = divmod(round(seconds), 60)
        return f"{minutes} m {rest:02d} s"
    return f"{seconds:.1f} s"


#: the symbols of the status of a step on the console
_STATUS = {
    "pending": "[dim]·[/dim]",
    "running": "[yellow]▶[/yellow]",
    "ok": "[green]✔[/green]",
    "failed": "[red]✘[/red]",
    "skipped": "[dim]–[/dim]",
}


def _table(title: str, results: Sequence[StepResult], total: float) -> Table:
    """Get the table of the steps of a run for the console."""
    table = Table(title=title, title_justify="left", expand=False)
    table.add_column("", width=1)
    table.add_column("step")
    table.add_column("command")
    table.add_column("time", justify="right")
    table.add_column("cores", justify="right")
    table.add_column("memory", justify="right")
    active = any(result.status == "running" for result in results)
    if active:
        table.add_column("progress", overflow="ellipsis", no_wrap=True)
    for result in results:
        if result.status == "running":
            wall = time.monotonic() - result._started
        else:
            wall = result.wall
        done = result.status in ("ok", "failed")
        running = result.status == "running"
        cells = [
            _STATUS[result.status],
            result.step.name,
            result.step.command,
            _duration(wall) if result.status != "pending" else "",
            f"{result.cpu / wall:.1f}" if (done or running) and wall > 0 else "",
            f"{result.peak_rss / 1024**3:.1f} GB" if result.peak_rss else "",
        ]
        if active:
            cells.append(result.progress[:PROGRESS_WIDTH] if running else "")
        table.add_row(*cells)
    table.caption = f"total {_duration(total)}"
    table.caption_justify = "left"
    return table


def _last_line(log: Path, offset: int) -> tuple[str, int]:
    """Get the last line of output written to a log since an offset.

    Returns:
        The line without the source location rich appends, empty if there is
        none, and the new offset.
    """
    try:
        with log.open("rb") as stream:
            stream.seek(offset)
            data = stream.read()
    except OSError:
        return "", offset
    lines = [
        _LEVEL.sub("", _LOCATION.sub("", line).strip())
        for line in data.decode("utf-8", "replace").splitlines()
    ]
    lines = [line for line in lines if line]
    return (lines[-1] if lines else ""), offset + len(data)


def _children_cpu() -> float:
    """Get the CPU seconds of the processes of this process which ended.

    The usage of the processes a process waited for includes the usage of the
    processes they waited for, so it is the CPU of a step and of its pool. On
    a platform without `resource` (Windows) it is 0.
    """
    try:
        import resource
    except ImportError:  # pragma: no cover - Windows
        return 0.0
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    return usage.ru_utime + usage.ru_stime


#: where the processes are read on linux
PROC = Path("/proc")


class _TreeUsage:
    """The CPU time and the memory of the processes of a step, sampled.

    The usage of the processes a process waited for reaches its parent
    (`_children_cpu`), but not the usage of a process whose parent did not
    wait for it, e.g. the workers of a pool started by a fork server, which
    ends with the step. So every process of the tree of the step is read from
    `/proc` once per sample, and the CPU time of a process is the last one
    seen, which misses at most the last interval of a process. Without `/proc`
    nothing is sampled.

    Attributes:
        cpu: the CPU seconds of the processes seen, by process and start.
        peak_rss: the peak of the resident bytes of the tree.
    """

    def __init__(self, pid: int) -> None:
        """Create the usage of the tree of a process."""
        self.pid = pid
        self.cpu: dict[tuple[int, int], float] = {}
        self.peak_rss = 0
        self._ticks = os.sysconf("SC_CLK_TCK") if hasattr(os, "sysconf") else 100
        self._page = os.sysconf("SC_PAGE_SIZE") if hasattr(os, "sysconf") else 4096

    def sample(self) -> None:
        """Read the processes of the tree of the step."""
        if not PROC.is_dir():
            return
        children: dict[int, list[int]] = {}
        stats: dict[int, tuple[int, float, int]] = {}
        for entry in PROC.iterdir():
            if not entry.name.isdigit():
                continue
            try:
                text = (entry / "stat").read_text(encoding="utf-8", errors="replace")
            except OSError:
                continue
            fields = text[text.rfind(")") + 1 :].split()
            try:
                ppid = int(fields[1])
                seconds = (int(fields[11]) + int(fields[12])) / self._ticks
                start = int(fields[19])
                rss = int(fields[21]) * self._page
            except IndexError, ValueError:
                continue
            pid = int(entry.name)
            stats[pid] = (start, seconds, rss)
            children.setdefault(ppid, []).append(pid)
        if self.pid not in stats:
            return
        resident = 0
        stack = [self.pid]
        while stack:
            pid = stack.pop()
            start, seconds, rss = stats[pid]
            self.cpu[(pid, start)] = seconds
            resident += rss
            stack.extend(children.get(pid, ()))
        self.peak_rss = max(self.peak_rss, resident)

    @property
    def cpu_seconds(self) -> float:
        """The CPU seconds of the processes of the tree which were seen."""
        return sum(self.cpu.values())


def _versions() -> dict[str, str | None]:
    """Get the versions of the packages of the run."""
    versions: dict[str, str | None] = {}
    for name in VERSIONS:
        try:
            versions[name] = version(name)
        except PackageNotFoundError:
            versions[name] = None
    return versions


def _run_step(
    result: StepResult, live: Live | None, update: object, env: Mapping[str, str]
) -> None:
    """Run the process of a step to its end, sampling its usage and progress.

    The CPU time of the step is the larger of the usage of the processes the
    step waited for and the sampled usage of its tree (`_TreeUsage`): the
    first is exact for the processes it covers, the second covers the
    processes whose parent did not wait for them.
    """
    result.status = "running"
    result._started = time.monotonic()
    cpu_before = _children_cpu()
    offset = 0
    progress = ""
    with result.log.open("wb") as stream:
        process = subprocess.Popen(
            result.argv, stdout=stream, stderr=subprocess.STDOUT, env=dict(env)
        )
        tree = _TreeUsage(process.pid)
        while True:
            try:
                process.wait(timeout=REFRESH)
                break
            except subprocess.TimeoutExpired:
                pass
            tree.sample()
            result.peak_rss = tree.peak_rss
            result.cpu = tree.cpu_seconds
            line, offset = _last_line(result.log, offset)
            progress = line or progress
            result.progress = progress
            if live is not None:
                live.update(update())  # ty: ignore[call-non-callable]
    result.wall = time.monotonic() - result._started
    result.cpu = max(_children_cpu() - cpu_before, tree.cpu_seconds)
    result.exit_code = process.returncode
    result.status = "ok" if process.returncode == 0 else "failed"


def run_pipeline(
    package: str,
    title: str,
    steps: Sequence[Step],
    results_dir: Path,
    out: Path,
    *,
    cores: int,
    keep_going: bool = False,
) -> list[StepResult]:
    """Run the steps of a pipeline one after the other, timed and logged.

    Args:
        package: the model package, the first argument of every command.
        title: the title of the run on the console.
        steps: the steps which run, in their order.
        results_dir: the results of the model, `{results}` of a step.
        out: the folder of the run, created; it holds `logs/`, `timings.tsv`
            and `pipeline.json`.
        cores: the processes of a step, `{cores}` of a step.
        keep_going: run the steps after a step which failed; they are skipped
            without.

    Returns:
        The results of the steps.
    """
    logs = out / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    values = {"cores": cores, "results": results_dir}
    results = [
        StepResult(
            step=step,
            argv=[
                sys.executable,
                "-m",
                "pkpdsim.cli",
                step.command,
                package,
                *step.arguments(values),
            ],
            log=logs / f"{k + 1:02d}_{step.name}.log",
        )
        for k, step in enumerate(steps)
    ]
    # the steps write a line of their progress every 30 s into their log
    env = {**os.environ, "SBMLSIM_PROGRESS": "log", "COLUMNS": "120"}
    started = datetime.datetime.now(datetime.UTC)
    begin = time.monotonic()
    record: dict[str, object] = {
        "package": package,
        "started": started.isoformat(timespec="seconds"),
        "cores": cores,
        "machine": {
            "platform": platform.platform(),
            "processor": platform.processor() or platform.machine(),
            "cpus": os.cpu_count(),
            "available_cpus": available_cpus(),
            "python": platform.python_version(),
        },
        "versions": _versions(),
        "steps": [
            {"name": result.step.name, "argv": result.argv[3:]} for result in results
        ],
    }

    def table() -> Table:
        return _table(title, results, time.monotonic() - begin)

    def write() -> None:
        record["ended"] = datetime.datetime.now(datetime.UTC).isoformat(
            timespec="seconds"
        )
        record["wall_s"] = round(time.monotonic() - begin, 1)
        record["results"] = [result.row() for result in results]
        (out / "pipeline.json").write_text(
            json.dumps(record, indent=2) + "\n", encoding="utf-8"
        )
        rows = [result.row() for result in results]
        header = list(rows[0])
        text = "\t".join(header) + "\n"
        for row in rows:
            text += "\t".join("" if row[k] is None else str(row[k]) for k in header)
            text += "\n"
        (out / "timings.tsv").write_text(text, encoding="utf-8")

    failed = False
    live = Live(table(), console=console, refresh_per_second=2, transient=False)
    interactive = console.is_terminal
    with live if interactive else _NoLive():
        for result in results:
            if failed and not keep_going:
                result.status = "skipped"
                continue
            if not interactive:
                console.print(
                    f"▶ {result.step.name}: {' '.join(result.argv[3:])}",
                    soft_wrap=True,
                    highlight=False,
                )
            _run_step(result, live if interactive else None, table, env)
            if not interactive:
                console.print(
                    f"{'✔' if result.status == 'ok' else '✘'} {result.step.name} "
                    f"{_duration(result.wall)}, {result.cores:.1f} cores, "
                    f"exit {result.exit_code}"
                )
            failed = failed or result.status == "failed"
            write()
            if interactive:
                live.update(table())
    write()
    if not interactive:
        console.print(table())
    return results


class _NoLive:
    """The context of a run without a live display, on no terminal."""

    def __enter__(self) -> None:
        return None

    def __exit__(self, *args: object) -> None:
        return None
