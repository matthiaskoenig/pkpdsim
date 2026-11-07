"""The HTML of the evidence extensions and of the battery hook in the report.

The report of an assessment is the template of sbmlsim
(`sbmlsim.credibility.report`); the extensions of pkpdsim give their blocks in
the card of a context of use, their sections and the fates of the doses of the
battery as HTML, rendered from the macros of `credibility_extensions.html` of
the templates of pkpdsim with the environment of the reports of sbmlsim.
"""

from __future__ import annotations

import math
from functools import cache
from pathlib import Path
from typing import TYPE_CHECKING, Final

from sbmlsim.report.templates import template_environment

if TYPE_CHECKING:
    import jinja2

__all__ = ["TEMPLATE", "TEMPLATE_DIR", "render_macro"]

#: the templates of pkpdsim
TEMPLATE_DIR: Final = Path(__file__).resolve().parents[1] / "resources" / "templates"

#: the macros of the blocks of the extensions
TEMPLATE: Final = "credibility_extensions.html"


def _number(value: float | None) -> str:
    """Show a number with three significant digits, `n/a` if not evaluable."""
    if value is None or math.isnan(value):
        return "n/a"
    return f"{value:.3g}" if abs(value) < 1000 else f"{value:.0f}"


@cache
def _environment() -> jinja2.Environment:
    """Get the environment of the macros, the one of the reports of sbmlsim."""
    env = template_environment(TEMPLATE_DIR)
    env.filters["number"] = _number
    return env


def render_macro(name: str, *args: object) -> str:
    """Render a macro of the blocks of the extensions.

    Args:
        name: the macro of `TEMPLATE`, e.g. `ratio_pairs`.
        *args: its arguments.

    Returns:
        The HTML.
    """
    module = _environment().get_template(TEMPLATE).module
    return str(getattr(module, name)(*args))
