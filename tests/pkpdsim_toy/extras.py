"""Extra actions of the toy model: they write into the results."""

import argparse
from pathlib import Path

from pkpdsim.model import PKPDModel


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Add the option of the greeting."""
    parser.add_argument("--greeting", default="hello", help="the text of the greeting")


def greet(model: PKPDModel, options: argparse.Namespace, results_dir: Path) -> None:
    """Write the greeting of the model into the results."""
    results_dir.mkdir(parents=True, exist_ok=True)
    (results_dir / "greeting.txt").write_text(
        f"{options.greeting} {model.name}", encoding="utf-8"
    )


def workers(model: PKPDModel, options: argparse.Namespace, results_dir: Path) -> None:
    """Write the processes of the simulator of the run into the results."""
    results_dir.mkdir(parents=True, exist_ok=True)
    (results_dir / "workers.txt").write_text(
        f"{model.name} {options.simulator.n_workers}", encoding="utf-8"
    )


def add_simulator_argument(parser: argparse.ArgumentParser) -> None:
    """Add an option which clashes with `options.simulator` of the run."""
    parser.add_argument("--simulator", default="roadrunner")
