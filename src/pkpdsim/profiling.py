"""The profile of a command of a model, with the sampling profiler of python 3.15.

`--profile` after the model package runs a command of the models under the
statistical sampling profiler of python 3.15 (`profiling.sampling`, PEP 799)
and writes the profile into the results of the model:

```bash
run_simulation pkdb_models.models.hctz --profile --action simulate --experiments all
```

The profiler samples the stack of the command at a fixed rate from a process of
its own, so the command runs at nearly its normal speed. The command runs as a
child of the profiler (`python -m profiling.sampling run -m pkpdsim.cli ...`),
whose binary profile is converted afterwards (`replay`) into the files of the
folder `profiling/<command>_<time>/` of the results of the model:

- `flamegraph.html`: the interactive flame graph of the sampled stacks;
- `heatmap/`: the source code of the sampled lines with their sample counts;
- `profile.pstats`: the statistics of the functions, for `pstats` and its
  viewers;
- `summary.txt`: the functions with the most samples;
- `profile.bin`: the samples, from which other formats are made and against
  which a later profile is compared (`replay --diff-flamegraph profile.bin`).

The profiler is part of python 3.15; with an older python `--profile` is an
error and the command does not run.
"""

from __future__ import annotations

import datetime
import importlib.util
import os
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING

from sbmlutils.console import console

if TYPE_CHECKING:
    from pkpdsim.model import PKPDModel

#: the option of the commands which runs a command under the profiler
PROFILE_OPTION = "--profile"

#: the python which has the sampling profiler
PROFILER_PYTHON = (3, 15)

#: the environment variable with the file into which the profiled command
#: writes its exit code: the profiler exits with 0 whatever the command did
EXIT_CODE_FILE = "PKPDSIM_PROFILE_EXIT_CODE"

#: the formats which the binary profile is converted into, as the option of
#: `replay` and the file in the folder of the profile
FORMATS: tuple[tuple[str, str], ...] = (
    ("--flamegraph", "flamegraph.html"),
    ("--heatmap", "heatmap"),
    ("--pstats", "profile.pstats"),
)


def profiler_available() -> bool:
    """Check whether the python has the sampling profiler `profiling.sampling`."""
    return (
        sys.version_info >= PROFILER_PYTHON
        and importlib.util.find_spec("profiling.sampling") is not None
    )


def record_exit_code(code: int | str | None) -> None:
    """Write the exit code of a profiled command, see `EXIT_CODE_FILE`.

    Args:
        code: the code of `SystemExit`: `None` is 0, a message is 1.
    """
    path = os.environ.get(EXIT_CODE_FILE)
    if not path:
        return
    value = 0 if code is None else code if isinstance(code, int) else 1
    Path(path).write_text(str(value), encoding="utf-8")


def profile_dir(model: PKPDModel, command: str) -> Path:
    """Get a new folder of the profile of a command in the results of a model.

    Args:
        model: the model.
        command: the name of the command, e.g. `run_simulation`.

    Returns:
        `profiling/<command>_<time>` in the results of the model, not created.
    """
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    return model.results_dir / "profiling" / f"{command}_{stamp}"


def run_profiled(
    model: PKPDModel, command: str, args: Sequence[str], output: Path | None = None
) -> int:
    """Run a command of a model under the sampling profiler, see the module.

    Args:
        model: the model.
        command: the name of the command, a key of `pkpdsim.cli.COMMANDS`.
        args: the arguments of the command after the model package, without
            `--profile`.
        output: the folder of the profile, `profile_dir` by default.

    Returns:
        The exit code of the command.

    Raises:
        RuntimeError: if the python has no sampling profiler.
    """
    if not profiler_available():
        version = ".".join(map(str, PROFILER_PYTHON))
        raise RuntimeError(
            f"{PROFILE_OPTION} needs the sampling profiler of python {version} "
            f"(profiling.sampling), this is python {sys.version.split()[0]}."
        )
    folder = output if output is not None else profile_dir(model, command)
    folder.mkdir(parents=True, exist_ok=True)
    binary = folder / "profile.bin"
    exit_code = folder / "exit_code"
    exit_code.unlink(missing_ok=True)
    profiler = [sys.executable, "-m", "profiling.sampling"]
    returned = subprocess.run(
        [
            *profiler,
            *("run", "--binary", "-o", str(binary)),
            *("-m", "pkpdsim.cli", command, model.package, *args),
        ],
        check=False,
        env={**os.environ, EXIT_CODE_FILE: str(exit_code)},
    ).returncode
    # the code of the command, which the profiler does not pass on
    code = int(exit_code.read_text()) if exit_code.is_file() else returned or 1
    exit_code.unlink(missing_ok=True)
    if not binary.is_file():
        console.print(f"The profiler wrote no profile into '{folder}'.", style="error")
        return code or 1
    for option, name in FORMATS:
        replay = subprocess.run(
            [*profiler, "replay", option, "-o", str(folder / name), str(binary)],
            check=False,
            capture_output=True,
            text=True,
        )
        if replay.returncode:
            raise RuntimeError(
                f"replay {option} failed:\n{replay.stderr}{replay.stdout}"
            )
    summary = subprocess.run(
        [*profiler, "replay", "--pstats", "--sort", "cumtime", "-l", "40", str(binary)],
        check=True,
        capture_output=True,
        text=True,
    )
    # the statistics, without the progress of the replay before them
    text = summary.stdout
    start = text.find("Profile Stats:")
    (folder / "summary.txt").write_text(
        text[start:] if start >= 0 else text, encoding="utf-8"
    )
    console.print(f"profile: file://{(folder / 'flamegraph.html').resolve()}")
    return code
