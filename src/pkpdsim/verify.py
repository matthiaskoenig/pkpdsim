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

import argparse
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


MODEL_HELP = "dotted name of the model package, e.g. pkdb_models.models.hctz"

ROOT_HELP = (
    "directory put first on the python path, e.g. a checkout of the repository "
    "of the model"
)


def parser() -> argparse.ArgumentParser:
    """Create the parser of the command line, which checks the arguments and shows the help.

    The arguments themselves are translated (`translate`) for `sbmlsim.verify`.
    """
    p = argparse.ArgumentParser(
        prog="pkpdsim-verify",
        description=__doc__.split("\n\n")[0] if __doc__ else None,
    )
    sub = p.add_subparsers(required=True)

    def command(name: str, help_: str) -> argparse.ArgumentParser:
        return sub.add_parser(name, help=help_, description=help_)

    s = command("dump", "run the experiments of a model and store everything")
    s.add_argument("--root", required=True, help=ROOT_HELP)
    s.add_argument("--model", required=True, metavar="PACKAGE", help=MODEL_HELP)
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
    s.add_argument(
        "--grid",
        type=_positive_int,
        default=None,
        metavar="N",
        help="output every simulation at N + 1 equidistant times from its start to "
        "its end and the times of its changes instead of the steps of the "
        "integrator, so two dumps share their times; a migration compare dumps "
        "both models with the same grid",
    )
    s.add_argument("experiments", nargs="*", help="only these experiment classes")

    s = command("compare", "compare two dumps")
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

    s = command("ranges", "fit mappings with data outside of the simulation")
    s.add_argument("dump")
    s.add_argument("--verbose", action="store_true", help="every mapping")

    s = command("fits", "a minimal fit of every fit definition of a model")
    s.add_argument("--root", required=True, help=ROOT_HELP)
    s.add_argument("--model", required=True, metavar="PACKAGE", help=MODEL_HELP)
    s.add_argument("--out", required=True, help="output directory of the fits")
    s.add_argument("--max-nfev", type=int, default=3, help="evaluations per start")
    s.add_argument("subsets", nargs="*", help="only these fit definitions")

    s = command("definitions", "the fit definitions of a model as json")
    s.add_argument("--root", required=True, help=ROOT_HELP)
    s.add_argument("--model", required=True, metavar="PACKAGE", help=MODEL_HELP)
    s.add_argument(
        "--out", required=True, help="json file the definitions are written to"
    )

    s = command("settings", "the settings a run of a model hands to sbmlsim as json")
    s.add_argument("--root", required=True, help=ROOT_HELP)
    s.add_argument("--model", required=True, metavar="PACKAGE", help=MODEL_HELP)
    s.add_argument("--out", required=True, help="json file the settings are written to")
    s.add_argument(
        "--experiment",
        help="the experiment class to run, the first of the registry by default",
    )

    s = command("png", "compare the PNG figures of two directories")
    s.add_argument("old")
    s.add_argument("new")
    s.add_argument(
        "--sids",
        action="store_true",
        help=(
            "pair an old figure with the new figure of the SId of its file name "
            "(sbmlsim.experiment.to_sid), for a port which renamed the keys"
        ),
    )
    s.add_argument("names", nargs="*", help="only files containing one of these")

    s = command("sbml", "compare two SBML files or the SBML files of two dirs")
    s.add_argument("old")
    s.add_argument("new")

    s = command(
        "scorecard",
        "the credibility of many models: the scorecard, the qualification of the "
        "ratios across models and the regressions against a baseline",
    )
    s.add_argument(
        "records", nargs="+", help="credibility.json files or directories with them"
    )
    s.add_argument(
        "--baseline", help="a scorecard TSV to compare with, a regression exits 1"
    )
    s.add_argument("--out", help="TSV file the scorecard is written to")
    s.add_argument("--qualification", help="TSV file the qualification is written to")
    s.add_argument(
        "--field",
        default="health",
        help="the field of the metadata of the numerator the ratios are pooled by",
    )
    s.add_argument(
        "--tolerance",
        type=float,
        default=0.05,
        help="the change of a metric which is no regression, default 0.05",
    )
    return p


def _positive_int(text: str) -> int:
    """Parse a positive integer of the command line."""
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError(f"a positive number, not {text}")
    return value


def verify(args: Sequence[str] | None = None) -> int:
    """Run a subcommand.

    Args:
        args: the command line arguments, `sys.argv[1:]` by default.

    Returns:
        The exit code of `sbmlsim.verify`: 0 nothing differs, 1 a difference, 2 a
        command which cannot be run, 3 a failure.
    """
    os.environ.setdefault("MPLBACKEND", "Agg")
    arguments = sys.argv[1:] if args is None else args
    try:
        parser().parse_args(arguments)
        translated = translate(arguments)
    except SystemExit as stop:
        return stop.code if isinstance(stop.code, int) else 2
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
