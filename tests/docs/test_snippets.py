"""Run the python snippets of the documentation.

Every ```python fence of a page of `docs/` is executed, in the order it
appears and in one namespace per page, so that a fragment which uses the
objects of the snippet above it runs as a reader would run it. The page runs
in a subprocess with warnings as errors, in a temporary working directory, so
a snippet which writes a file writes nothing into the repository.

A fence whose first line is the comment `# not executed` is skipped: it is a
fragment which needs what a page cannot provide (a checkout of pkdb_data, a
model of a reader). The convention is documented in `docs/development.md`.
"""

import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

#: the root of the repository
REPO_DIR = Path(__file__).parent.parent.parent

#: the markdown sources of the documentation
DOCS_DIR = REPO_DIR / "docs"

#: fences which start with this comment are fragments and are not executed
SKIP_MARKER = "# not executed"

#: a fenced python block of a markdown page
BLOCK = re.compile(r"^```python\n(.*?)^```", re.S | re.M)

#: runs the blocks of one page in one namespace, `blocks.json` beside it
RUNNER = """
import json
import sys

blocks = json.loads(open("blocks.json", encoding="utf-8").read())
namespace = {"__name__": "__main__"}
for number, code in blocks:
    try:
        exec(compile(code, f"<block {number}>", "exec"), namespace)
    except BaseException:
        print(f"snippet {number} failed", file=sys.stderr)
        raise
"""


def executable_blocks(page: Path) -> list[tuple[int, str]]:
    """The python fences of a page which are meant to run, numbered from 1.

    Args:
        page: the markdown page.

    Returns:
        The number and the source of every fence without the skip marker.
    """
    blocks = BLOCK.findall(page.read_text(encoding="utf-8"))
    return [
        (number, code)
        for number, code in enumerate(blocks, start=1)
        if not code.lstrip().startswith(SKIP_MARKER)
    ]


#: the pages with at least one snippet to run, `docs/api/` holds only the
#: directives of the API reference and is left out
PAGES: list[str] = sorted(
    page.name for page in DOCS_DIR.glob("*.md") if executable_blocks(page)
)


def test_pages_are_found() -> None:
    """The pages with snippets are discovered."""
    assert "index.md" in PAGES


@pytest.mark.parametrize("name", PAGES)
def test_snippets(name: str, tmp_path: Path) -> None:
    """Every snippet of a page runs, in the order and the namespace of the page."""
    blocks = executable_blocks(DOCS_DIR / name)
    (tmp_path / "blocks.json").write_text(json.dumps(blocks), encoding="utf-8")
    (tmp_path / "run_page.py").write_text(RUNNER, encoding="utf-8")
    env = dict(os.environ, PYTHONUTF8="1")
    result = subprocess.run(
        [sys.executable, "-W", "error", "run_page.py"],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        encoding="utf-8",
        check=False,
    )
    assert result.returncode == 0, f"{name}\n{result.stdout}\n{result.stderr}"
