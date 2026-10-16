"""Tests of the release notes, from the previous release and the new build."""

import dataclasses
import json
import shutil
from collections.abc import Sequence
from pathlib import Path

import pytest
from release_git import git

from pkpdsim import __version__
from pkpdsim.model import PKPDModel
from pkpdsim.release.notes import (
    CoreSummary,
    DataChanges,
    FitChange,
    LibraryChange,
    core_summaries,
    data_changes,
    experiment_changes,
    fit_changes,
    github_releases,
    library_changes,
    model_changes,
    write_release_notes,
)
from pkpdsim.release.project import ProjectContext
from pkpdsim.release.source import ReleaseError, ReleaseRepo, SourceRepo
from pkpdsim_toy import MODEL

OLD_FIT = "PK_20260101_000000__abcde"
NEW_FIT = "PK_20261016_185018__3efe3"
PD_FIT = "PD_20261016_185018__12345"
DATA_DIR = "src/pkpdsim_toy/data"
FITS_DIR = "src/pkpdsim_toy/fits"
COMMIT = "0123456789abcdef0123456789abcdef01234567"


def lock(*packages: tuple[str, str]) -> str:
    text = 'version = 1\n\n[[package]]\nname = "toy-model"\nversion = "0.1.0"\n'
    text += 'source = { editable = "." }\n'
    for name, version in packages:
        text += f'\n[[package]]\nname = "{name}"\nversion = "{version}"\n'
        text += 'source = { registry = "https://pypi.org/simple" }\n'
    return text


OLD_LOCK = lock(
    ("numpy", "2.0.0"),
    ("pkpdsim", "0.6.0"),
    ("removed-package", "1.0.0"),
    ("sbmlsim", "0.9.6"),
    ("sbmlutils", "0.16.0"),
)
NEW_LOCK = lock(
    ("added-package", "0.1.0"),
    ("numpy", "2.1.0"),
    ("pkpdsim", "0.7.0"),
    ("sbmlsim", "0.9.7"),
    ("sbmlutils", "0.16.0"),
)


pytestmark = pytest.mark.usefixtures("git_env")


def write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def parameters(values: dict[str, float]) -> str:
    return json.dumps({"sets": [{"sid": "best", "values": values, "cost": 1.0}]})


def toy_sbml() -> str:
    return (MODEL.path / "models" / "toy.xml").read_text()


def commit_release(path: Path, tag: str) -> ReleaseRepo:
    git(path, "add", "-A")
    git(path, "commit", "-m", f"release {tag}")
    git(path, "tag", tag)
    return ReleaseRepo(path=path, last_tag=tag)


@pytest.fixture
def model() -> PKPDModel:
    return dataclasses.replace(MODEL, fits=[NEW_FIT])


@pytest.fixture
def release(tmp_path: Path) -> ReleaseRepo:
    """A previous release 0.1.0 and the new build in the working tree."""
    path = tmp_path / "toy-model"
    path.mkdir()
    git(path, "init", "-b", "main")
    write(path / "uv.lock", OLD_LOCK)
    write(path / DATA_DIR / "drug" / "Same2000" / "Same2000.tsv", "a\n1\n")
    write(path / DATA_DIR / "drug" / "Changed2001" / "Changed2001.tsv", "a\n1\n")
    write(path / DATA_DIR / "drug" / "Removed2002" / "Removed2002.tsv", "a\n1\n")
    write(path / DATA_DIR / "collect_data.py", "print()\n")
    write(path / "simulations" / "A" / "index.html", "<html/>")
    write(path / "simulations" / "B" / "index.html", "<html/>")
    write(path / "simulations" / "_figures" / "a.svg", "<svg/>")
    write(path / "simulations" / "index.html", "<html/>")
    unset = '<parameter id="r" value="NaN" constant="false"/>\n<parameter '
    write(path / "models" / "toy.xml", toy_sbml().replace("<parameter ", unset, 1))
    write(path / "models" / "removed.xml", toy_sbml())
    write(path / "models" / "toy.md", "# toy")
    write(path / FITS_DIR / OLD_FIT / "parameters.json", parameters({"k": 0.2}))
    write(path / FITS_DIR / PD_FIT / "parameters.json", parameters({"e": 1.0}))
    repo = commit_release(path, "0.1.0")

    write(path / "uv.lock", NEW_LOCK)
    write(path / DATA_DIR / "drug" / "Changed2001" / "Changed2001.tsv", "a\n2\n")
    shutil.rmtree(path / DATA_DIR / "drug" / "Removed2002")
    write(path / DATA_DIR / "other" / "Added2003" / "Added2003.tsv", "a\n1\n")
    shutil.rmtree(path / "simulations" / "B")
    write(path / "simulations" / "C" / "index.html", "<html/>")
    changed = toy_sbml().replace('value="0.2"', 'value="0.3"')
    changed = changed.replace(
        "<parameter ",
        '<parameter id="k2" value="1" constant="true"/>\n'
        '<parameter id="r" constant="false"/>\n<parameter ',
        1,
    )
    write(path / "models" / "toy.xml", changed)
    (path / "models" / "removed.xml").unlink()
    write(path / "models" / "added.xml", toy_sbml())
    shutil.rmtree(path / FITS_DIR)
    write(path / FITS_DIR / NEW_FIT / "parameters.json", parameters({"k": 0.25}))
    return repo


@pytest.fixture
def old_layout(tmp_path: Path) -> ReleaseRepo:
    """A previous release of the old layout: no lock file, no fit records."""
    path = tmp_path / "old-model"
    path.mkdir()
    git(path, "init", "-b", "main")
    write(path / DATA_DIR / "drug" / "Same2000" / "Same2000.tsv", "a\n1\n")
    write(path / "simulations" / "A" / "index.html", "<html/>")
    write(path / "simulations" / "_figures" / "a.svg", "<svg/>")
    write(path / "models" / "toy.xml", toy_sbml())
    write(path / "fit" / "20260421_191243__251fc" / "index.html", "<html/>")
    repo = commit_release(path, "0.7.1")

    write(path / "uv.lock", NEW_LOCK)
    write(path / FITS_DIR / NEW_FIT / "parameters.json", parameters({"k": 0.25}))
    return repo


def context(
    model: PKPDModel, version: str, dependencies: tuple[str, ...] = ()
) -> ProjectContext:
    assert model.release is not None
    return ProjectContext(
        model=model,
        settings=model.release,
        version=version,
        python="3.14",
        dependencies=dependencies,
        constraints=(),
        pkpdsim_path=None,
        year=2026,
        requires_python=">=3.13",
    )


RELEASES = {
    "pkpdsim": {
        "0.5.0": "- too old. Not shown.\n",
        "0.6.0": "- the old version. Not shown.\n",
        "0.6.1": (
            "# Release notes\n\nIntro text.\n\n## New features\n"
            "- First of 0.6.1. Second sentence.\n"
            "  - nested, not shown\n"
            "- sixth bullet, over the limit.\n"
        ),
        "0.7.0": (
            "## Breaking changes\r\n- third bullet! More.\r\n- fourth bullet\r\n"
            "- fifth bullet? More.\r\n"
            "- `run_report` recreates the report, e.g. of a fit. More.\r\n"
        ),
        "0.8.0": "- newer than the release. Not shown.\n",
    },
    "sbmlsim": {
        "0.9.7": (
            "- **A bold sentence.** More.\n"
            "- **Bold** ([#1](https://x.org/1)). lowercase more.\n"
            "- Versions 1.5 and 2.0, i.e. two. More.\n"
        )
    },
}


def fake_fetch(repository: str) -> dict[str, str]:
    return RELEASES[repository]


def failing_fetch(repository: str) -> dict[str, str]:
    raise ReleaseError("'gh api' could not run")


def test_library_changes() -> None:
    assert library_changes(OLD_LOCK, NEW_LOCK) == [
        LibraryChange("added-package", None, "0.1.0"),
        LibraryChange("numpy", "2.0.0", "2.1.0"),
        LibraryChange("pkpdsim", "0.6.0", "0.7.0"),
        LibraryChange("removed-package", "1.0.0", None),
        LibraryChange("sbmlsim", "0.9.6", "0.9.7"),
    ]
    assert library_changes(None, NEW_LOCK) == []


LOCAL_PKPDSIM_LOCK = """\
version = 1

[[package]]
name = "pkpdsim"
source = { editable = "/home/user/pkpdsim" }

[[package]]
name = "sbmlsim"
version = "0.9.7"
source = { registry = "https://pypi.org/simple" }
"""


def test_library_changes_local_package_without_version() -> None:
    """uv locks a local package with a dynamic version without a version."""
    changes = library_changes(OLD_LOCK, LOCAL_PKPDSIM_LOCK)
    assert LibraryChange("pkpdsim", "0.6.0", "editable /home/user/pkpdsim") in changes


def test_data_changes(release: ReleaseRepo) -> None:
    assert data_changes(release, DATA_DIR) == DataChanges(
        added=["other/Added2003"],
        removed=["drug/Removed2002"],
        changed=["drug/Changed2001"],
    )


def test_experiment_changes(release: ReleaseRepo) -> None:
    assert experiment_changes(release) == (["C"], ["B"])


def test_model_changes(release: ReleaseRepo) -> None:
    changes = {change.file: change for change in model_changes(release)}
    assert sorted(changes) == ["added.xml", "removed.xml", "toy.xml"]
    toy = changes["toy.xml"]
    assert toy.status == "changed"
    assert toy.added == {"species": 0, "reactions": 0, "parameters": 1}
    assert toy.removed == {"species": 0, "reactions": 0, "parameters": 0}
    assert toy.changed_parameters == [("k", 0.2, 0.3)]
    assert changes["added.xml"].status == "added"
    assert changes["added.xml"].added == {"species": 1, "reactions": 1, "parameters": 1}
    assert changes["removed.xml"].status == "removed"
    assert changes["removed.xml"].removed["reactions"] == 1


def test_fit_changes(release: ReleaseRepo, model: PKPDModel) -> None:
    assert fit_changes(release, model) == [
        FitChange(PD_FIT, None, []),
        FitChange(OLD_FIT, NEW_FIT, [("k", 0.2, 0.25)]),
    ]


def test_core_summaries() -> None:
    changes = library_changes(OLD_LOCK, NEW_LOCK)
    assert core_summaries(changes, fetch=fake_fetch) == [
        CoreSummary(
            "pkpdsim",
            "0.6.0",
            "0.7.0",
            [
                "third bullet!",
                "fourth bullet",
                "fifth bullet?",
                "`run_report` recreates the report, e.g. of a fit.",
                "First of 0.6.1.",
            ],
            total=6,
        ),
        CoreSummary(
            "sbmlsim",
            "0.9.6",
            "0.9.7",
            [
                "**A bold sentence.**",
                "**Bold** ([#1](https://x.org/1)).",
                "Versions 1.5 and 2.0, i.e. two.",
            ],
            total=3,
        ),
    ]


def test_core_summaries_added_library() -> None:
    changes = [LibraryChange("sbmlutils", None, "0.16.0")]
    assert core_summaries(changes, fetch=failing_fetch) == [
        CoreSummary("sbmlutils", None, "0.16.0", [])
    ]


def test_core_summaries_fetch_fails() -> None:
    changes = library_changes(OLD_LOCK, NEW_LOCK)
    summaries = core_summaries(changes, fetch=failing_fetch)
    assert [(s.name, s.fetched, s.bullets) for s in summaries] == [
        ("pkpdsim", False, []),
        ("sbmlsim", False, []),
    ]


def test_github_releases() -> None:
    calls: list[list[str]] = []

    def runner(args: Sequence[str], cwd: Path) -> str:
        calls.append([str(a) for a in args])
        return '{"tag":"0.7.0","body":"- a"}\n{"tag":"v0.6.0","body":null}\n'

    assert github_releases("pkpdsim", runner) == {"0.7.0": "- a", "0.6.0": ""}
    assert calls[0][:3] == ["gh", "api", "repos/matthiaskoenig/pkpdsim/releases"]


def test_write_release_notes(release: ReleaseRepo, model: PKPDModel) -> None:
    source = SourceRepo(root=release.path, commit=COMMIT)
    path = write_release_notes(
        context(model, "0.2.0"), release, model, source, fetch=fake_fetch
    )
    assert path == release.path / "release-notes" / "0.2.0.md"
    text = path.read_text()
    assert text.startswith("# Release notes for toy-model v0.2.0\n")
    assert f"Made from pkdb_models {COMMIT[:10]} with pkpdsim {__version__}." in text
    for heading in ("Core libraries", "Libraries", "Data", "Analysis"):
        assert f"## {heading}\n" in text
    assert "pkpdsim 0.6.0 → 0.7.0 (latest 5 of 6):\n" in text
    assert "sbmlsim 0.9.6 → 0.9.7:\n" in text
    assert "third bullet!" in text
    assert "sixth bullet" not in text
    assert "| `numpy` | 2.0.0 | 2.1.0 |" in text
    assert "`other/Added2003`" in text
    assert "`drug/Removed2002`" in text
    assert "`drug/Changed2001`" in text
    assert "`k` 0.2 → 0.3" in text
    assert "`k` 0.2 → 0.25" in text
    assert "\u2014" not in text
    assert "\n\n\n" not in text


def test_write_release_notes_fetch_fails(
    release: ReleaseRepo, model: PKPDModel
) -> None:
    source = SourceRepo(root=release.path, commit=COMMIT)
    text = write_release_notes(
        context(model, "0.2.0"), release, model, source, fetch=failing_fetch
    ).read_text()
    assert "pkpdsim 0.6.0 → 0.7.0: the release notes could not be fetched." in text


def minimal_release(tmp_path: Path, old_lock: str) -> ReleaseRepo:
    """A previous release 0.1.0 with a lock file, a model and a fit record."""
    path = tmp_path / "toy-model"
    path.mkdir()
    git(path, "init", "-b", "main")
    write(path / "uv.lock", old_lock)
    write(path / "models" / "toy.xml", toy_sbml())
    write(path / FITS_DIR / NEW_FIT / "parameters.json", parameters({"k": 0.25}))
    return commit_release(path, "0.1.0")


def test_notes_no_changes(tmp_path: Path, model: PKPDModel) -> None:
    release = minimal_release(tmp_path, OLD_LOCK)
    source = SourceRepo(root=release.path, commit=COMMIT)
    text = write_release_notes(
        context(model, "0.1.1"), release, model, source, fetch=fake_fetch
    ).read_text()
    assert text.count("No changes.") == 6


def test_notes_direct_dependencies(tmp_path: Path, model: PKPDModel) -> None:
    release = minimal_release(
        tmp_path, lock(("numpy", "2.0.0"), ("pint", "0.25"), ("idna", "3.11"))
    )
    write(
        release.path / "uv.lock",
        lock(("numpy", "2.1.0"), ("pint", "0.26"), ("idna", "3.20")),
    )
    source = SourceRepo(root=release.path, commit=COMMIT)
    dependencies = ("NumPy>=2.0", "pkpdsim[docs]>=0.7")
    text = write_release_notes(
        context(model, "0.1.1", dependencies), release, model, source, fetch=fake_fetch
    ).read_text()
    libraries = text.split("## Libraries\n\n")[1].split("\n## Data")[0]
    assert libraries == (
        "| Package | Previous | New |\n"
        "| --- | --- | --- |\n"
        "| `numpy` | 2.0.0 | 2.1.0 |\n"
        "\n"
        "2 other packages changed.\n"
        "\n"
        "<details>\n"
        "<summary>All changed packages</summary>\n"
        "\n"
        "| Package | Previous | New |\n"
        "| --- | --- | --- |\n"
        "| `idna` | 3.11 | 3.20 |\n"
        "| `numpy` | 2.0.0 | 2.1.0 |\n"
        "| `pint` | 0.25 | 0.26 |\n"
        "\n"
        "</details>\n"
    )


@pytest.mark.parametrize(
    "content", ["{", '{"sets": []}', '{"sets": [{"sid": "best"}]}']
)
@pytest.mark.parametrize("where", ["tag", "worktree"])
def test_fit_changes_malformed_record(
    release: ReleaseRepo, model: PKPDModel, content: str, where: str
) -> None:
    if where == "tag":
        shutil.rmtree(release.path / FITS_DIR)
        write(release.path / FITS_DIR / OLD_FIT / "parameters.json", content)
        git(release.path, "add", "-A")
        git(release.path, "commit", "-m", "broken record")
        git(release.path, "tag", "0.1.1")
        release = ReleaseRepo(path=release.path, last_tag="0.1.1")
        write(release.path / FITS_DIR / NEW_FIT / "parameters.json", parameters({}))
        name = f"0.1.1:{FITS_DIR}/{OLD_FIT}/parameters.json"
    else:
        write(release.path / FITS_DIR / NEW_FIT / "parameters.json", content)
        name = str(release.path / FITS_DIR / NEW_FIT / "parameters.json")
    with pytest.raises(ReleaseError, match="malformed") as err:
        fit_changes(release, model)
    assert name in str(err.value)


def test_notes_old_layout(old_layout: ReleaseRepo, model: PKPDModel) -> None:
    source = SourceRepo(root=old_layout.path, commit=COMMIT)
    text = write_release_notes(
        context(model, "0.8.0"), old_layout, model, source, fetch=failing_fetch
    ).read_text()
    assert "no lock file in the previous release" in text
    assert "no fit record in the previous release" in text
    assert f"`{NEW_FIT}`" in text


def test_notes_no_previous_release(release: ReleaseRepo, model: PKPDModel) -> None:
    first = ReleaseRepo(path=release.path, last_tag=None)
    source = SourceRepo(root=release.path, commit=COMMIT)
    text = write_release_notes(
        context(model, "0.1.0"), first, model, source, fetch=fake_fetch
    ).read_text()
    assert "First release." in text
    for item in ("drug/Same2000", "other/Added2003", "`A`", "`C`", "toy.xml", NEW_FIT):
        assert item in text
    assert "drug/Removed2002" not in text
    assert "no fit record in the previous release" not in text
