"""The inline plugin of tox: the libpython of the interpreter on the search path.

The extension of libroadrunner links against `libpython`. An interpreter of uv
ships it in its `lib` directory, where the loader of Linux does not look, so
`import roadrunner` fails with `cannot open shared object file`, or, worse,
loads the `libpython` of another interpreter of the same version from the
system, a second copy of the runtime in the process. Every test environment
puts the `lib` directory of its own interpreter first on `LD_LIBRARY_PATH`, as
`.github/actions/setup` does for the workflows.
"""

import subprocess
import sys
from pathlib import Path

from tox.plugin import impl
from tox.tox_env.api import ToxEnv
from tox.tox_env.python.api import Python


@impl
def tox_before_run_commands(tox_env: ToxEnv) -> None:
    """Put the `lib` directory of the interpreter of the environment on the path.

    Args:
        tox_env: the environment whose commands run next.
    """
    if not sys.platform.startswith("linux") or not isinstance(tox_env, Python):
        return
    base_prefix = subprocess.run(
        [str(tox_env.env_python()), "-c", "import sys; print(sys.base_prefix)"],
        capture_output=True,
        check=True,
        text=True,
    ).stdout.strip()
    lib = Path(base_prefix) / "lib"
    if not any(lib.glob("libpython*.so*")):
        return
    paths = [str(lib)]
    if current := tox_env.environment_variables.get("LD_LIBRARY_PATH"):
        paths.append(current)
    tox_env.conf["set_env"].update({"LD_LIBRARY_PATH": ":".join(paths)})
