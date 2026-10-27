"""The numbers which the page of the substances states for its creatinine example.

`docs/pbpk/substances.md` builds a creatinine model whose plasma species
start at the baseline `CRTNPROD_k / GFR` and describes how the plasma moves
from that start. The test executes the snippets of the page up to the
definition, builds it and checks the stated numbers, so that a change of the
library which changes them fails here and the page is rewritten.
"""

from pathlib import Path

import numpy as np
import pytest

from pkpdsim.pbpk import PBPKDefinition
from pkpdsim.pbpk.builder import build

from .test_snippets import DOCS_DIR, executable_blocks

DAY = 24 * 60.0
"""One day (min)."""


def _creatinine() -> PBPKDefinition:
    """The definition `creatinine` of the snippets of the page."""
    namespace: dict[str, object] = {"__name__": "__main__"}
    for number, code in executable_blocks(DOCS_DIR / "pbpk" / "substances.md"):
        exec(compile(code, f"<block {number}>", "exec"), namespace)
        definition = namespace.get("creatinine")
        if isinstance(definition, PBPKDefinition):
            return definition
    raise AssertionError("no snippet of substances.md defines `creatinine`")


def test_the_creatinine_example_settles_at_its_baseline(tmp_path: Path) -> None:
    """The arterial plasma dips, returns to the baseline in two days; the arm vein settles above.

    The page: the arterial plasma drops about 17 % below the baseline within
    the first 20 minutes and is back within 0.2 % after two days; the arm vein
    plasma settles about 5.5 % above the baseline.
    """
    import roadrunner

    result = build(_creatinine(), tmp_path)
    r = roadrunner.RoadRunner(str(result.flat_path))
    r.integrator.absolute_tolerance = 1e-12
    r.integrator.relative_tolerance = 1e-10
    baseline = r["baseline_crtn"]
    r.timeCourseSelections = ["time", "[Car_crtn]", "Carm_crtn"]
    s = r.simulate(0, 4 * DAY, int(4 * DAY) + 1)
    time, arterial, arm = s[:, 0], s[:, 1] / baseline, s[:, 2] / baseline
    lowest = int(np.argmin(arterial))
    assert arterial[lowest] == pytest.approx(0.83, abs=0.01)
    assert time[lowest] < 20.0
    assert arterial[int(2 * DAY)] == pytest.approx(0.998, abs=0.001)
    assert arterial[-1] == pytest.approx(1.0, abs=1e-4)
    assert arm[-1] == pytest.approx(1.055, abs=0.001)
