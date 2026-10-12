"""Write `tests/data/figures-sbmlsim-0.8.5.json`, the figures of sbmlsim 0.8.5.

The fixture is the serialization (`pkpdsim.verify.serialize`) of the figures of
`build_figures` of `tests/figure_builders.py`, as pkpdsim 0.4.1 with sbmlsim 0.8.5
stored them; the tests compare it with the figures sbmlsim 0.9 builds. sbmlsim
0.8.5 cannot be installed next to the environment of the repository, so this
script is not part of the test suite: it needs an environment with sbmlsim
0.8.5 and pkpdsim 0.4.1 and is run by hand, only when `build_figures` changes:

    uv venv venv085 && uv pip install --python venv085 "sbmlsim==0.8.5" "pkpdsim==0.4.1"
    venv085/bin/python -I scripts/figures_sbmlsim_085.py tests/data/figures-sbmlsim-0.8.5.json

`tests/figure_builders.py` is the one source of the figures (it imports nothing
but sbmlsim); the script loads it by path, since `-I` leaves the folder of the
script out of `sys.path`. Run again with the same versions, the output is the
committed file byte for byte.
"""

import json
import logging
import sys
from pathlib import Path

import sbmlsim

import pkpdsim
from pkpdsim.verify import serialize

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tests"))
from figure_builders import FigureStudy, build_figures

logging.disable(logging.CRITICAL)


figures = build_figures(FigureStudy(), setattr)
content = {
    "sbmlsim": sbmlsim.__version__,
    "pkpdsim": pkpdsim.__version__,
    "figures": {
        sid: json.dumps(serialize(figure), sort_keys=True, default=repr)
        for sid, figure in figures.items()
    },
}
with open(sys.argv[1], "w", encoding="utf-8") as f:
    json.dump(content, f, indent=1, sort_keys=True)
    f.write("\n")
