"""Git repositories for the tests of the release of a model repository.

The tests run real git on temporary repositories: a source repository (the
monorepo) with a copy of the toy model package, and a release repository with
its bare remote. The fixtures are in `conftest.py`; the helpers here are imported
by the tests.
"""

from pathlib import Path

from pkpdsim.release.source import run

#: a git without the configuration of the user, with a fixed author
GIT_ENV = {
    "GIT_AUTHOR_NAME": "Test",
    "GIT_AUTHOR_EMAIL": "test@example.org",
    "GIT_COMMITTER_NAME": "Test",
    "GIT_COMMITTER_EMAIL": "test@example.org",
    "GIT_CONFIG_GLOBAL": "/dev/null",
    "GIT_CONFIG_SYSTEM": "/dev/null",
}


def git(cwd: Path, *args: str) -> str:
    """Run git in a folder and return its output."""
    return run(["git", *args], cwd)


def init_repo(path: Path) -> None:
    """Create a git repository on the branch `main`."""
    path.mkdir(parents=True, exist_ok=True)
    git(path, "init", "-b", "main")
