"""The pipeline of a model: its commands, side by side where they may, timed and logged.

The analysis of a model is a set of the commands of pkpdsim: the factory, the
simulation experiments, their bands over the draws of the fitted parameters,
the sensitivity analysis, the assessment of the credibility, the refits of the
adopted fits. Which of them a model runs, with which options and on how much of
the machine differs between models, so a model package declares its pipeline
(`Pipeline`, usually `PIPELINE` in its `pipeline.py`) and its model object
references it (`pipeline=".pipeline:PIPELINE"`); `run_pipeline <package>` runs
it (`pkpdsim.cli.pipeline_command`):

```python
PIPELINE = Pipeline(
    description="The complete analysis of the model from its adopted fits.",
    steps=(
        Step("factory", "run_simulation", ("--action", "factory"), cores=1),
        Step("fit_PK", "run_fit", ("--subset", "PK", "--cores", "{cores}"), cores=10, after=("factory",)),
        Step("simulate", "run_simulation", ("--action", "simulate", "--experiments", "all"), cores=3, after=("factory",)),
        Step("sensitivity", "run_sensitivity", ("--workers", "{cores}"), cores=6, memory_gb=4, after=("factory",)),
        Step("credibility", "run_credibility", ("--uncertainty",), after=("factory",)),
    ),
)
```

Every step runs its command in a process of its own
(`python -m pkpdsim.cli <command> <package> ...`), so a step which fails or
runs out of memory leaves the runner and the results of the other steps
intact. A command parallelizes internally (the runs of a fit, the scans of a
profile, the experiments of a simulation, the points of a sensitivity design),
but rarely over the whole machine all the time: it has serial phases, and a
pool waits for its longest task. So the runner starts every step whose steps
`after` ended as soon as its cores and its memory fit into what the running
steps leave of the budget (`--cores`, the CPUs available to the process by
default, and `--memory`, the memory available at the start):

- `cores` of a step is the cores it is given, `{cores}` of its arguments and
  the CPUs its processes see (`PYTHON_CPU_COUNT`, so a pool which sizes itself
  by the CPUs takes these), and the running steps are given no more cores
  than the run has: pools which together ask for more processes than there
  are CPUs slow the step which takes longest, and with it the run. A step
  without `cores` is given the whole budget and runs alone, so a pipeline
  which declares no `cores` runs its steps one after the other in their
  order; running steps side by side is a declaration of the pipeline, since
  only the pipeline knows how much of the machine a step uses.
- `memory_gb` of a step is its peak memory, an estimate (`peak_rss_gb` of a
  run in `timings.tsv`), which the running steps do not exceed together.
- `after` of a step is the steps whose output it reads: it starts after they
  ended and is skipped if one of them failed. A step without `after` reads the
  output of no other step of the pipeline.

A step whose needs exceed the budget runs when no other step runs, so every
pipeline runs on every machine. The steps start in their order and a step which
fits starts before an earlier one which does not, so the longest steps first is
the order which packs the machine best. `--jobs` bounds the steps which run at
the same time: `--jobs 1` runs one step at a time, each on the cores it is
given, which times the steps alone.

A run writes into its own folder below `pipeline/` of the results of the model:
the output of every step in `logs/<step>.log`, the times in `timings.tsv` (the
start of the step in the run, its wall clock, the CPU time of the step and of
every process it started, the cores it was given and the cores it used on
average, the peak memory of its processes) and the versions, the machine and
the budget in `pipeline.json`, and shows the steps with their times and the
last line of progress of the running ones on the console, with the wall clock
of the run against the sum of the times of the steps.
"""

from __future__ import annotations

import datetime
import json
import math
import os
import platform
import re
import signal
import subprocess
import sys
import threading
import time
from collections.abc import Callable, Generator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import dataclass, field
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import IO

from rich.live import Live
from rich.table import Table
from sbmlsim.cgroup import available_cpus
from sbmlutils.console import console

__all__ = [
    "PLACEHOLDERS",
    "Pipeline",
    "Step",
    "StepResult",
    "available_memory_gb",
    "run_pipeline",
]

#: the placeholders of the arguments of a step: the cores of a step and the
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

#: the seconds between two refreshes of the console and samples of the usage
REFRESH = 1.0

#: the seconds between two looks at the processes of the running steps
POLL = 0.2

#: the seconds the processes of a step have to end after an interrupt, before
#: they are killed
GRACE = 10.0

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
            cores of the step, `{results}` the results of the model.
        description: what the step does, listed by `run_pipeline --list`.
        after: the names of the steps whose output the step reads, declared
            before it: it starts after they ended and is skipped if one of
            them failed.
        cores: the cores the step is given, `{cores}` of its arguments and the
            CPUs its processes see, which the runner counts against the
            budget of the run; `None` (the default) is the whole budget, so
            the step runs alone.
        memory_gb: the peak memory of the step in GB, an estimate the running
            steps do not exceed together; 0 (the default) counts no memory.
    """

    name: str
    command: str
    args: tuple[str, ...] = ()
    description: str = ""
    after: tuple[str, ...] = ()
    cores: int | None = None
    memory_gb: float = 0.0

    def __post_init__(self) -> None:
        """Check the name, the placeholders of the arguments and the needs.

        Raises:
            ValueError: for a name which is no word, arguments or steps `after`
                which are a string, a placeholder which is none of
                `PLACEHOLDERS`, cores below 1 or a negative memory.
        """
        if not STEP_NAME.match(self.name):
            raise ValueError(f"The name of the step '{self.name}' is no word.")
        for attribute in ("args", "after"):
            if isinstance(getattr(self, attribute), str):
                raise ValueError(
                    f"The {attribute} of the step '{self.name}' are a string, "
                    "not a tuple."
                )
        for arg in self.args:
            for placeholder in re.findall(r"\{(\w*)\}", arg):
                if placeholder not in PLACEHOLDERS:
                    raise ValueError(
                        f"The argument '{arg}' of the step '{self.name}' has the "
                        f"placeholder '{{{placeholder}}}', none of {PLACEHOLDERS}."
                    )
        if self.cores is not None and self.cores < 1:
            raise ValueError(
                f"The step '{self.name}' has {self.cores} cores, at least 1."
            )
        if not self.memory_gb >= 0:
            raise ValueError(
                f"The step '{self.name}' has the memory {self.memory_gb} GB, "
                "at least 0."
            )

    def given_cores(self, budget: int) -> int:
        """Get the cores the step is given on a budget, `{cores}` of its arguments.

        A step without `cores` is given the budget, a step with more cores than
        the budget the budget.
        """
        return budget if self.cores is None else min(self.cores, budget)

    def arguments(self, values: Mapping[str, object]) -> list[str]:
        """Get the arguments with the placeholders replaced by their values.

        A path is written with forward slashes, which every platform reads, so
        `{results}/fit` is one path on Windows as well.
        """

        def text(value: object) -> str:
            return value.as_posix() if isinstance(value, Path) else str(value)

        return [
            re.sub(r"\{(\w+)\}", lambda m: text(values[m.group(1)]), arg)
            for arg in self.args
        ]


@dataclass(frozen=True)
class Pipeline:
    """The steps of the analysis of a model.

    Attributes:
        steps: the steps in their order, at least one; a step starts before
            the steps after it which wait for the same cores.
        description: what the pipeline does.
    """

    steps: tuple[Step, ...]
    description: str = ""

    def __post_init__(self) -> None:
        """Check the steps.

        Raises:
            ValueError: for no step, two steps of one name, or a step `after`
                which is not declared before the step.
        """
        if not self.steps:
            raise ValueError("A pipeline has at least one step.")
        names = [step.name for step in self.steps]
        twice = sorted({name for name in names if names.count(name) > 1})
        if twice:
            raise ValueError(f"The pipeline has two steps of the names {twice}.")
        for k, step in enumerate(self.steps):
            for name in step.after:
                if name not in names[:k]:
                    raise ValueError(
                        f"The step '{step.name}' runs after '{name}', which is "
                        "no step declared before it."
                    )

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

        A step `after` which does not run counts as ended: its output is the
        one of an earlier run.

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
        given_cores: the cores it was given, `{cores}` of its arguments.
        status: `pending`, `running`, `ok`, `failed`, `skipped` or
            `interrupted`.
        exit_code: the exit code of its process, `None` before it ended.
        start: the seconds from the start of the run to the start of the step.
        wall: the seconds it ran.
        cpu: the CPU seconds of its process and of every process it started.
        peak_rss: the peak of the resident bytes of its processes, sampled.
        progress: the last line of progress of its output.
    """

    step: Step
    argv: list[str]
    log: Path
    given_cores: int = 1
    status: str = "pending"
    exit_code: int | None = None
    start: float = 0.0
    wall: float = 0.0
    cpu: float = 0.0
    peak_rss: int = 0
    progress: str = ""
    _started: float = field(default=0.0, repr=False)

    @property
    def cores(self) -> float:
        """The cores the step used on average, its CPU time over its wall time."""
        return self.cpu / self.wall if self.wall > 0 else 0.0

    @property
    def ran(self) -> bool:
        """Whether the step started."""
        return self.status not in ("pending", "skipped")

    def elapsed(self) -> float:
        """Get the seconds the step ran, so far for a running step."""
        if self.status == "running":
            return time.monotonic() - self._started
        return self.wall

    def row(self) -> dict[str, object]:
        """Get the row of the step in `timings.tsv`."""
        return {
            "step": self.step.name,
            "command": self.step.command,
            "status": self.status,
            "exit_code": self.exit_code,
            "start_s": round(self.start, 1) if self.ran else None,
            "wall_s": round(self.wall, 1),
            "cpu_s": round(self.cpu, 1),
            "cores_given": self.given_cores,
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
    "interrupted": "[red]■[/red]",
}


def _summary(results: Sequence[StepResult], total: float) -> str:
    """Get the wall clock of a run against the sum of the times of its steps."""
    steps = sum(result.elapsed() for result in results)
    text = f"total {_duration(total)}, the steps {_duration(steps)}"
    if total > 0 and steps > 0:
        text += f" ({steps / total:.1f}x)"
    return text


def _table(title: str, results: Sequence[StepResult], total: float) -> Table:
    """Get the table of the steps of a run for the console.

    The cores of a step which ran are the cores it used on average over the
    cores it was given, of a waiting step the cores it will be given.
    """
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
        wall = result.elapsed()
        if result.ran and wall > 0:
            cores = f"{result.cpu / wall:.1f}/{result.given_cores}"
        elif result.status == "pending":
            cores = f"[dim]{result.given_cores}[/dim]"
        else:
            cores = ""
        cells = [
            _STATUS[result.status],
            result.step.name,
            result.step.command,
            _duration(wall) if result.ran else "",
            cores,
            f"{result.peak_rss / 1024**3:.1f} GB" if result.peak_rss else "",
        ]
        if active:
            running = result.status == "running"
            cells.append(result.progress[:PROGRESS_WIDTH] if running else "")
        table.add_row(*cells)
    table.caption = _summary(results, total)
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


#: where the processes and the memory are read on linux
PROC = Path("/proc")


def available_memory_gb() -> float:
    """Get the memory available for new processes in GB, `MemAvailable` of linux.

    Returns:
        The memory, infinite on a platform without `/proc/meminfo`.
    """
    try:
        text = (PROC / "meminfo").read_text(encoding="utf-8")
    except OSError:
        return math.inf
    match = re.search(r"^MemAvailable:\s+(\d+)\s+kB", text, re.MULTILINE)
    return int(match.group(1)) / 1024**2 if match else math.inf


@dataclass
class _Processes:
    """The processes of the machine at a moment, read from `/proc`.

    Attributes:
        children: the processes by their parent.
        stats: the start, the CPU seconds and the resident bytes by process.
    """

    children: dict[int, list[int]] = field(default_factory=dict)
    stats: dict[int, tuple[int, float, int]] = field(default_factory=dict)

    @classmethod
    def read(cls) -> _Processes:
        """Read the processes, none without `/proc`."""
        processes = cls()
        if not PROC.is_dir():
            return processes
        ticks = os.sysconf("SC_CLK_TCK") if hasattr(os, "sysconf") else 100
        page = os.sysconf("SC_PAGE_SIZE") if hasattr(os, "sysconf") else 4096
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
                seconds = (int(fields[11]) + int(fields[12])) / ticks
                start = int(fields[19])
                rss = int(fields[21]) * page
            except IndexError, ValueError:
                continue
            pid = int(entry.name)
            processes.stats[pid] = (start, seconds, rss)
            processes.children.setdefault(ppid, []).append(pid)
        return processes


class _TreeUsage:
    """The CPU time and the memory of the processes of a step, sampled.

    The usage of the processes a process waited for reaches its parent (and
    the runner, which waits for the process of the step), but not the usage of
    a process whose parent did not wait for it, e.g. the workers of a pool
    started by a fork server, which ends with the step. So every process of the
    tree of the step is read from `/proc` once per sample, and the CPU time of
    a process is the last one seen, which misses at most the last interval of
    a process. The tree of a step holds its own processes only, so the steps
    which run side by side are counted apart. Without `/proc` nothing is
    sampled.

    Attributes:
        cpu: the CPU seconds of the processes seen, by process and start.
        peak_rss: the peak of the resident bytes of the tree.
    """

    def __init__(self, pid: int) -> None:
        """Create the usage of the tree of a process."""
        self.pid = pid
        self.cpu: dict[tuple[int, int], float] = {}
        self.peak_rss = 0

    def sample(self, processes: _Processes | None = None) -> None:
        """Read the processes of the tree of the step.

        Args:
            processes: the processes of the machine, read now without.
        """
        processes = _Processes.read() if processes is None else processes
        if self.pid not in processes.stats:
            return
        resident = 0
        stack = [self.pid]
        while stack:
            pid = stack.pop()
            start, seconds, rss = processes.stats[pid]
            self.cpu[(pid, start)] = seconds
            resident += rss
            stack.extend(processes.children.get(pid, ()))
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


class _Schedule:
    """Which steps of a run start, by their steps `after` and the budget.

    A step `after` which is no step of the run counts as ended `ok`.

    Attributes:
        results: the results of the steps of the run, in their order.
        cores: the cores the running steps share.
        memory_gb: the memory the running steps share, in GB.
        jobs: the most steps which run at the same time, no bound for `None`.
        keep_going: start the steps which do not depend on a step which
            failed; without, no step starts after a step failed.
    """

    def __init__(
        self,
        results: Sequence[StepResult],
        cores: int,
        memory_gb: float = math.inf,
        jobs: int | None = None,
        keep_going: bool = False,
    ) -> None:
        """Create the schedule of the results of a run."""
        self.results = list(results)
        self.cores = cores
        self.memory_gb = memory_gb
        self.jobs = jobs
        self.keep_going = keep_going
        self._by_name = {result.step.name: result for result in self.results}

    @property
    def running(self) -> list[StepResult]:
        """The steps which run."""
        return [result for result in self.results if result.status == "running"]

    @property
    def pending(self) -> list[StepResult]:
        """The steps which wait to run."""
        return [result for result in self.results if result.status == "pending"]

    def _before(self, result: StepResult) -> list[str]:
        """Get the statuses of the steps `after` of a step which run."""
        return [
            self._by_name[name].status
            for name in result.step.after
            if name in self._by_name
        ]

    def skip(self) -> list[StepResult]:
        """Skip the waiting steps which will not run.

        Without `keep_going` every waiting step after a step failed, with it
        the steps after a step which did not end `ok`, transitively, since the
        steps `after` of a step precede it.

        Returns:
            The steps which were skipped now.
        """
        failed = any(result.status == "failed" for result in self.results)
        skipped = []
        for result in self.pending:
            if (failed and not self.keep_going) or any(
                status in ("failed", "skipped", "interrupted")
                for status in self._before(result)
            ):
                result.status = "skipped"
                skipped.append(result)
        return skipped

    def ready(self) -> list[StepResult]:
        """Get the waiting steps which start now, in their order.

        A step starts when its steps `after` ended `ok` and its cores and its
        memory fit beside the running steps and the steps which start before
        it, or when no step runs; a step which does not fit lets the steps
        after it start.
        """
        running = self.running
        cores = sum(result.given_cores for result in running)
        memory = sum(result.step.memory_gb for result in running)
        count = len(running)
        starts: list[StepResult] = []
        for result in self.pending:
            if self.jobs is not None and count >= self.jobs:
                break
            if any(status != "ok" for status in self._before(result)):
                continue
            fits = (
                cores + result.given_cores <= self.cores
                and memory + result.step.memory_gb <= self.memory_gb
            )
            if count and not fits:
                continue
            starts.append(result)
            cores += result.given_cores
            memory += result.step.memory_gb
            count += 1
        return starts


@dataclass
class _Process:
    """The process of a running step.

    Attributes:
        result: the result of the step.
        popen: the process.
        stream: the log the process writes to.
        tree: the usage of the processes of the step.
        offset: the bytes of the log read for the progress.
    """

    result: StepResult
    popen: subprocess.Popen[bytes]
    stream: IO[bytes]
    tree: _TreeUsage
    offset: int = 0

    def ended(self, block: bool = False) -> bool:
        """Reap the process if it ended, with the usage of its processes.

        With `os.wait4` the CPU time of the step is the larger of the usage of
        the processes the step waited for, which `wait4` gives for the process
        of the step alone, and the sampled usage of its tree: the first is
        exact for the processes it covers, the second covers the processes
        whose parent did not wait for them.

        Args:
            block: wait for the process to end.
        """
        cpu = 0.0
        if hasattr(os, "wait4"):
            pid, status, usage = os.wait4(self.popen.pid, 0 if block else os.WNOHANG)
            if pid == 0:
                return False
            self.popen.returncode = os.waitstatus_to_exitcode(status)
            cpu = usage.ru_utime + usage.ru_stime
        elif block:  # pragma: no cover - Windows
            self.popen.wait()
        elif self.popen.poll() is None:  # pragma: no cover - Windows
            return False
        self.stream.close()
        result = self.result
        result.wall = time.monotonic() - result._started
        result.cpu = max(cpu, self.tree.cpu_seconds)
        result.peak_rss = self.tree.peak_rss
        result.exit_code = self.popen.returncode
        result.status = "ok" if self.popen.returncode == 0 else "failed"
        return True

    def signal(self, kill: bool) -> None:
        """Send an interrupt, or the kill, to the processes of the step.

        On posix the step is a process group of its own, which the processes
        it starts join, so the signal reaches its pools as well.
        """
        try:
            if hasattr(os, "killpg"):
                os.killpg(self.popen.pid, signal.SIGKILL if kill else signal.SIGINT)
            elif kill:  # pragma: no cover - Windows
                self.popen.kill()
            else:  # pragma: no cover - Windows
                self.popen.terminate()
        except ProcessLookupError, PermissionError:
            pass


def _start(result: StepResult, env: Mapping[str, str], begin: float) -> _Process:
    """Start the process of a step.

    The processes of the step see the cores it was given as the CPUs of the
    machine (`PYTHON_CPU_COUNT`), so a pool which sizes itself by the CPUs
    takes these and the running steps do not ask for more cores than the run
    has.
    """
    result.status = "running"
    result._started = time.monotonic()
    result.start = result._started - begin
    stream = result.log.open("wb")
    popen = subprocess.Popen(
        result.argv,
        stdout=stream,
        stderr=subprocess.STDOUT,
        env={**env, "PYTHON_CPU_COUNT": str(result.given_cores)},
        start_new_session=hasattr(os, "killpg"),
    )
    return _Process(result, popen, stream, _TreeUsage(popen.pid))


def _wait(seconds: float) -> None:
    """Wait for the processes of the running steps."""
    time.sleep(seconds)


def _stop(processes: Sequence[_Process]) -> None:
    """Stop the processes of the running steps, which are `interrupted` then.

    The processes get an interrupt and have `GRACE` seconds to end, a second
    interrupt of the runner ends the grace; then every process of the steps
    is killed.
    """
    for process in processes:
        process.signal(kill=False)
    alive = list(processes)
    deadline = time.monotonic() + GRACE
    try:
        while alive and time.monotonic() < deadline:
            time.sleep(POLL)
            alive = [process for process in alive if not process.ended()]
    except KeyboardInterrupt:
        pass
    for process in processes:
        process.signal(kill=True)
    for process in alive:
        if process.result.status == "running":
            process.ended(block=True)
    for process in processes:
        process.result.status = "interrupted"


def _raise_interrupt(signum: int, frame: object) -> None:
    """Raise `KeyboardInterrupt` for a signal which ends the runner."""
    raise KeyboardInterrupt


@contextmanager
def _interrupts() -> Generator[None]:
    """Turn the signals which end the runner into `KeyboardInterrupt`.

    The steps are process groups of their own, which the interrupt of a
    terminal (Ctrl-C) or the end of its session does not reach, so the runner
    stops them on `SIGTERM` and `SIGHUP` as on `SIGINT`.
    """
    if threading.current_thread() is not threading.main_thread():
        yield
        return
    names = ("SIGTERM", "SIGHUP")
    numbers = [getattr(signal, name) for name in names if hasattr(signal, name)]
    previous = {number: signal.signal(number, _raise_interrupt) for number in numbers}
    try:
        yield
    finally:
        for number, handler in previous.items():
            signal.signal(number, handler)


def run_pipeline(
    package: str,
    title: str,
    steps: Sequence[Step],
    results_dir: Path,
    out: Path,
    *,
    cores: int,
    memory_gb: float | None = None,
    jobs: int | None = None,
    keep_going: bool = False,
    program: Sequence[str] | None = None,
) -> list[StepResult]:
    """Run the steps of a pipeline, side by side where they may, timed and logged.

    Args:
        package: the model package, the first argument of every command.
        title: the title of the run on the console.
        steps: the steps which run, in their order.
        results_dir: the results of the model, `{results}` of a step.
        out: the folder of the run, created; it holds `logs/`, `timings.tsv`
            and `pipeline.json`.
        cores: the cores the running steps share, `{cores}` of a step without
            `cores`.
        memory_gb: the memory the running steps share in GB, the memory
            available at the start without (`available_memory_gb`).
        jobs: the most steps which run at the same time, no bound without.
        keep_going: run the steps which do not depend on a step which failed;
            without, the steps which did not start are skipped.
        program: the start of the command line of a step, before the command,
            the package and the arguments (`python -m pkpdsim.cli` without).

    Returns:
        The results of the steps.

    Raises:
        KeyboardInterrupt: after the processes of the running steps were
            stopped on an interrupt (`SIGINT`, `SIGTERM`, `SIGHUP`).
    """
    logs = out / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    memory_gb = available_memory_gb() if memory_gb is None else memory_gb
    program = [sys.executable, "-m", "pkpdsim.cli"] if program is None else program
    results = []
    for k, step in enumerate(steps):
        given = step.given_cores(cores)
        values = {"cores": given, "results": results_dir}
        results.append(
            StepResult(
                step=step,
                argv=[*program, step.command, package, *step.arguments(values)],
                log=logs / f"{k + 1:02d}_{step.name}.log",
                given_cores=given,
            )
        )
    schedule = _Schedule(results, cores, memory_gb, jobs, keep_going)
    # the steps write a line of their progress every 30 s into their log
    env = {**os.environ, "SBMLSIM_PROGRESS": "log", "COLUMNS": "120"}
    started = datetime.datetime.now(datetime.UTC)
    begin = time.monotonic()
    record: dict[str, object] = {
        "package": package,
        "started": started.isoformat(timespec="seconds"),
        "cores": cores,
        "memory_gb": round(memory_gb, 1) if math.isfinite(memory_gb) else None,
        "jobs": jobs,
        "machine": {
            "platform": platform.platform(),
            "processor": platform.processor() or platform.machine(),
            "cpus": os.cpu_count(),
            "available_cpus": available_cpus(),
            "python": platform.python_version(),
        },
        "versions": _versions(),
        "steps": [
            {
                "name": result.step.name,
                "argv": result.argv[len(program) :],
                "after": list(result.step.after),
                "cores": result.given_cores,
                "memory_gb": result.step.memory_gb,
            }
            for result in results
        ],
    }

    def table() -> Table:
        return _table(title, results, time.monotonic() - begin)

    def write() -> None:
        record["ended"] = datetime.datetime.now(datetime.UTC).isoformat(
            timespec="seconds"
        )
        record["wall_s"] = round(time.monotonic() - begin, 1)
        record["steps_wall_s"] = round(sum(result.wall for result in results), 1)
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

    interactive = console.is_terminal

    def say(text: str) -> None:
        if not interactive:
            console.print(text, soft_wrap=True, highlight=False)

    live = Live(table(), console=console, refresh_per_second=2, transient=False)
    update: Callable[[], None] = (
        (lambda: live.update(table())) if interactive else (lambda: None)
    )
    processes: list[_Process] = []
    with live if interactive else _NoLive(), _interrupts():
        try:
            sampled = 0.0
            while True:
                for result in schedule.skip():
                    say(f"– {result.step.name} skipped")
                for result in schedule.ready():
                    say(
                        f"▶ {result.step.name}: {' '.join(result.argv[len(program) :])}"
                    )
                    processes.append(_start(result, env, begin))
                if not processes:
                    break
                _wait(POLL)
                ended = [process for process in processes if process.ended()]
                for process in ended:
                    processes.remove(process)
                    result = process.result
                    say(
                        f"{'✔' if result.status == 'ok' else '✘'} "
                        f"{result.step.name} {_duration(result.wall)}, "
                        f"{result.cores:.1f} cores, exit {result.exit_code}"
                    )
                if ended:
                    write()
                if ended or time.monotonic() - sampled >= REFRESH:
                    sampled = time.monotonic()
                    snapshot = _Processes.read()
                    for process in processes:
                        process.tree.sample(snapshot)
                        result = process.result
                        result.peak_rss = process.tree.peak_rss
                        result.cpu = process.tree.cpu_seconds
                        line, process.offset = _last_line(result.log, process.offset)
                        result.progress = line or result.progress
                    update()
        except KeyboardInterrupt:
            _stop(processes)
            for result in schedule.pending:
                result.status = "skipped"
            write()
            update()
            raise
        finally:
            for process in processes:
                process.stream.close()
        write()
        update()
    if not interactive:
        console.print(table())
    return results


class _NoLive:
    """The context of a run without a live display, on no terminal."""

    def __enter__(self) -> None:
        return None

    def __exit__(self, *args: object) -> None:
        return None
