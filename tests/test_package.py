import subprocess
import sys

import pkpdsim


def test_version() -> None:
    assert pkpdsim.__version__ == "0.10.1"


def test_all() -> None:
    assert pkpdsim.__all__ == ["__version__"]


def test_import_is_light() -> None:
    """`import pkpdsim` imports none of the modules and none of their dependencies."""
    code = (
        "import sys, pkpdsim; "
        "print(sorted(m for m in sys.modules "
        "if m.split('.')[0] in {'roadrunner', 'sbmlutils', 'libsbml'} "
        "or m.startswith('pkpdsim.')))"
    )
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    assert result.stdout.strip() == "[]"
