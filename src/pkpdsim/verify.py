"""Verify that a change of a model changes only what it is meant to change.

`pkpdsim-verify` is a wrapper of `python -m sbmlsim.run verify`, see
`sbmlsim.verify` for the subcommands (`dump`, `compare`, `ranges`, `fits`,
`definitions`, `settings`, `png`, `sbml`, `scorecard`) and the exit codes. It
keeps the options of pkpdsim and adds what a model of pkpdsim is:

- `--model PACKAGE` is the dotted name of a model package which defines its model
  object `MODEL` (`pkpdsim.model.PKPDModel`), i.e. `--project PACKAGE:MODEL`,
- a dump records the sex of the build of every task (`builds`), from
  `PKPDSimulationExperiment.variant_of`, and the versions of pkpdsim, sbmlsim and
  sbmlutils (the core libraries of the release settings of the model),
- `scorecard` pools with the evidence extensions of `pkpdsim.credibility`
  (`EXTENSIONS`), the qualification of the ratio pairs by `--field` (default
  `health`) in the wording of the limits of Guest et al.

A model before and after a change, with tight tolerances so that the solver does
not hide a difference:

```bash
git worktree add --detach ../base develop
(cd ../base && uv run python pkdb_models/models/hctz/models/factory.py)
uv run python pkdb_models/models/hctz/models/factory.py
uv run pkpdsim-verify dump --root ../base --model pkdb_models.models.hctz \
    --out base.pkl --atol 1e-16 --rtol 1e-12 --grid 1000
uv run pkpdsim-verify dump --root . --model pkdb_models.models.hctz \
    --out new.pkl --order base.pkl --atol 1e-16 --rtol 1e-12 --grid 1000
uv run pkpdsim-verify compare base.pkl new.pkl
```
"""

import os
import sys
from collections.abc import Iterator, Sequence

#: the subcommands which load a model
MODEL_COMMANDS = ("dump", "fits", "definitions", "settings")

#: the attribute of a model package which holds its model
MODEL_ATTRIBUTE = "MODEL"

#: the evidence extensions the scorecard judges with
EXTENSIONS = "pkpdsim.credibility:EXTENSIONS"

#: the extension which pools the qualification
POOL = "ratios"

#: the field of the metadata the qualification is pooled by without `--field`
DEFAULT_FIELD = "health"


class UsageError(ValueError):
    """The command line is no command of `pkpdsim-verify`."""


def _model_arguments(arguments: Sequence[str]) -> Iterator[str]:
    """Replace `--model PACKAGE` by `--project PACKAGE:MODEL`."""
    pending = False
    for argument in arguments:
        if pending:
            pending = False
            yield f"{argument}:{MODEL_ATTRIBUTE}"
        elif argument == "--model":
            pending = True
            yield "--project"
        elif argument.startswith("--model="):
            yield f"--project={argument.removeprefix('--model=')}:{MODEL_ATTRIBUTE}"
        else:
            yield argument
    if pending:
        # sbmlsim reports the missing value
        yield from ()


def translate(arguments: Sequence[str]) -> list[str]:
    """Get the arguments of `python -m sbmlsim.run verify` for those of `pkpdsim-verify`.

    Args:
        arguments: the arguments after `pkpdsim-verify`.

    Returns:
        The arguments after `verify`.

    Raises:
        UsageError: if a subcommand which loads a model gets no `--root` or no
            `--model`.
    """
    translated = list(_model_arguments(arguments))
    command = next((a for a in translated if not a.startswith("-")), None)
    helping = any(a in ("-h", "--help") for a in translated)
    if command in MODEL_COMMANDS and not helping:
        for option, given in (("--root", "--root"), ("--model", "--project")):
            if not any(a == given or a.startswith(f"{given}=") for a in translated):
                raise UsageError(f"the following arguments are required: {option}")
    if command == "scorecard" and not helping:
        translated += ["--extensions", EXTENSIONS, "--pool", POOL]
        if not any(a == "--field" or a.startswith("--field=") for a in translated):
            translated += ["--field", DEFAULT_FIELD]
    return translated


def verify(args: Sequence[str] | None = None) -> int:
    """Run a subcommand.

    Args:
        args: the command line arguments, `sys.argv[1:]` by default.

    Returns:
        The exit code of `sbmlsim.verify`: 0 nothing differs, 1 a difference, 2 a
        command which cannot be run, 3 a failure.
    """
    os.environ.setdefault("MPLBACKEND", "Agg")
    try:
        translated = translate(sys.argv[1:] if args is None else args)
    except UsageError as err:
        print(f"pkpdsim-verify: error: {err}", file=sys.stderr)
        return 2

    from sbmlsim.run.cli import main as run_main

    return run_main(["verify", *translated])


def main() -> None:
    """The command `pkpdsim-verify`."""
    from rich.console import Console
    from sbmlsim.log import enable_rich_logging

    enable_rich_logging(console=Console(stderr=True))
    sys.exit(verify())


if __name__ == "__main__":
    # the workers of a fit start with forkserver/spawn and import this module
    main()
