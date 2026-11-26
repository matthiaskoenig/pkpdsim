"""What pkpdsim imports: no matplotlib and no seaborn, also when it draws a figure.

The figures of the zonation are charts which the renderer of sbmlsim-core draws,
so a process which imports the package and renders them never loads either.
"""

import subprocess
import sys

CODE = """
import sys
import numpy as np
import pandas as pd

import pkpdsim
from pkpdsim.zonation import observers
from pkpdsim.zonation.figures import (
    ZoneProfile, layer_profile, lobule_map, space_time, zone_profiles,
)
from pkpdsim.zonation.geometry import Layer

n = 6
t = np.linspace(0.0, 60.0, 13)
cols = {"time": t}
for layer in (Layer.PLASMA, Layer.CELL):
    for k, sid in enumerate(observers.profile_ids("drg", layer, n), start=1):
        cols[sid] = k * t / 60.0
cols[observers.outflow_concentration("drg")] = t / 60.0
result = pd.DataFrame(cols)

profile = ZoneProfile("p", "mM", tuple(observers.profile_ids("drg", Layer.PLASMA, n)))
figures = [
    zone_profiles(result, [profile], times=(10.0, 60.0)),
    space_time(result, "drg", layers=(Layer.PLASMA, Layer.CELL)),
    layer_profile(result, "drg", time=30.0),
    lobule_map(result, "drg", time=30.0, zones=n, layer=Layer.CELL),
]
for figure in figures:
    assert figure.svg().lstrip().startswith("<svg"), "not an SVG"
print(sorted(m for m in ("matplotlib", "seaborn") if m in sys.modules))
"""


def test_a_zonation_figure_imports_no_matplotlib_or_seaborn() -> None:
    """`import pkpdsim` and the four figures, rendered to SVG, leave both out."""
    result = subprocess.run(
        [sys.executable, "-c", CODE],
        capture_output=True,
        text=True,
        check=False,
        timeout=300,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip().splitlines()[-1] == "[]"


def test_importing_pkpdsim_alone_imports_no_matplotlib_or_seaborn() -> None:
    """The package itself, before any figure is built."""
    code = (
        "import sys, pkpdsim, pkpdsim.zonation.figures; "
        "print(sorted(m for m in ('matplotlib', 'seaborn') if m in sys.modules))"
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        check=False,
        timeout=300,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip().splitlines()[-1] == "[]"
