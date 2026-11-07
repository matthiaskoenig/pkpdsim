"""Tests of the release hook of pkpdsim and of the release settings of a model.

The release itself is sbmlsim's (`sbmlsim.release`, tested there); pkpdsim adds
the framework `pkpdsim` with its core libraries, the PK parts of the page of
the credibility and the studies of the data as pkpdsim 0.19 named them.
"""

import dataclasses
import subprocess
import sys
from datetime import timedelta
from pathlib import Path, PurePosixPath

import pytest
from credibility_tables import assess, fit_table
from pkpdsim_credibility.credibility import CREDIBILITY
from release_git import git, init_repo
from sbmlsim.credibility import CredibilityAssessment
from sbmlsim.release.build import BuildPaths
from sbmlsim.release.docs import write_docs
from sbmlsim.release.hooks import ReleaseHook
from sbmlsim.release.notes import data_changes, data_dirs
from sbmlsim.release.project import ProjectContext
from sbmlsim.release.settings import ReleaseSettings
from sbmlsim.release.source import ReleaseRepo

from pkpdsim.model import PKPDModel, PKPDReleaseSettings
from pkpdsim.release.hook import HOOK, PKReleaseHook
from pkpdsim_toy import MODEL

assert MODEL.release is not None
SETTINGS: ReleaseSettings = MODEL.release


def test_a_model_is_released_with_pkpdsim() -> None:
    assert isinstance(SETTINGS, PKPDReleaseSettings)
    assert SETTINGS.framework == "pkpdsim"
    assert SETTINGS.core_libraries == ("pkpdsim", "sbmlsim", "sbmlutils")
    assert SETTINGS.library_owner == "matthiaskoenig"
    assert SETTINGS.hook == "pkpdsim.release.hook:HOOK"
    assert SETTINGS.archive is None
    assert MODEL.release_hook() is HOOK
    assert isinstance(HOOK, PKReleaseHook)
    assert SETTINGS.missing() == []


def _archive(model: PKPDModel) -> str | None:
    assert model.release is not None
    return model.release.archive


def test_the_archive_of_a_model_is_the_archive_of_its_release(tmp_path: Path) -> None:
    model = MODEL.replace(archive_path=MODEL.path / "models" / "toy_model.omex")
    assert model.release is not None
    assert model.release.archive == "models/toy_model.omex"
    assert _archive(MODEL.replace(archive_path="toy_model.omex")) == "toy_model.omex"
    # the archive of the release is the one of the model, also of a copy
    assert _archive(model.replace(archive_path=None)) is None
    with pytest.raises(ValueError, match="lies outside of the root"):
        MODEL.replace(archive_path=tmp_path / "toy_model.omex")
    # a model without release settings has no archive of a release
    assert MODEL.replace(release=None, archive_path=tmp_path / "a.omex").release is None


def test_the_settings_a_package_gives_are_kept() -> None:
    """A field the package gives keeps its value, also the default of sbmlsim."""
    required = ["title", "repository", "doi", "authors", "keywords"]
    required += ["docker_image", "funding", "description", "first_year"]
    fields = {name: getattr(SETTINGS, name) for name in required}
    defaults = PKPDReleaseSettings(**fields)
    assert (defaults.framework, defaults.hook) == (
        "pkpdsim",
        "pkpdsim.release.hook:HOOK",
    )
    own = PKPDReleaseSettings(
        **fields,
        framework="sbmlsim",
        core_libraries=("sbmlsim", "sbmlutils"),
        hook=None,
    )
    model = MODEL.replace(release=own)
    assert model.release is not None
    assert model.release.framework == "sbmlsim"
    assert model.release.core_libraries == ("sbmlsim", "sbmlutils")
    assert model.release.hook is None
    assert type(model.release_hook()) is ReleaseHook
    other = dataclasses.replace(SETTINGS, framework="other", hook="other:HOOK")
    model = MODEL.replace(release=other)
    assert model.release is not None
    assert (model.release.framework, model.release.hook) == ("other", "other:HOOK")
    assert model.release.core_libraries == ("pkpdsim", "sbmlsim", "sbmlutils")


def test_the_release_settings_of_sbmlsim_are_refused() -> None:
    """A model is released with pkpdsim, the settings of sbmlsim would not say so."""
    fields = {f.name: getattr(SETTINGS, f.name) for f in dataclasses.fields(SETTINGS)}
    with pytest.raises(ValueError, match="are a PKPDReleaseSettings"):
        MODEL.replace(release=ReleaseSettings(**fields))


def test_the_model_imports_no_release() -> None:
    code = (
        "import sys\n"
        "import pkpdsim.model\n"
        "print(sorted(m for m in sys.modules if m.startswith(('pkpdsim.release', "
        "'sbmlsim.release.workflow', 'jinja2'))))\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, encoding="utf-8", check=True
    )
    assert result.stdout.strip() == "[]"


@pytest.mark.parametrize(
    ("data_dir", "path", "study"),
    [
        ("data/hctz", "Beermann1976/Beermann1976_Tab1.tsv", "hctz/Beermann1976"),
        ("data", "hctz/Beermann1976/Beermann1976_Tab1.tsv", "hctz/Beermann1976"),
        ("data/hctz", "Beermann1976/raw/Fig1.png", "hctz/Beermann1976"),
        ("data", "hctz/Beermann1976.tsv", None),
        ("data/hctz", "collect.py", None),
        ("data", "README.md", None),
        ("experiments/data", "hctz/Beermann1976/a.tsv", None),
    ],
)
def test_the_study_of_a_file_of_the_data(
    data_dir: str, path: str, study: str | None
) -> None:
    assert HOOK.data_study(PurePosixPath(data_dir), PurePosixPath(path)) == study


def _data(package: Path, files: dict[str, str]) -> None:
    for relative, text in files.items():
        file = package / relative
        file.parent.mkdir(parents=True, exist_ok=True)
        file.write_text(text, encoding="utf-8")


def test_the_release_notes_name_the_studies_as_pkpdsim_0_19(
    tmp_path: Path, git_env: None
) -> None:
    """The studies of the data are `<substance>/<study>`, the names of 0.19.

    The previous release, made with pkpdsim 0.19, and the new build have the
    same layout; an unchanged study is no change, so the first release notes
    after the release moved into sbmlsim rename no study.
    """
    repo = tmp_path / "release"
    init_repo(repo)
    package = repo / "src" / "pkdb_models" / "models" / "apixaban"
    _data(
        package,
        {
            "data/apixaban/Frost2013/Frost2013_Tab1.tsv": "a\n",
            "data/apixaban/Kreutz2017/Kreutz2017_Fig1.tsv": "b\n",
            "data/edoxaban/Mendell2011/Mendell2011_Tab1.tsv": "c\n",
            "data/apixaban/collect.py": "x\n",
        },
    )
    git(repo, "add", "-A")
    git(repo, "commit", "-m", "0.1.0")
    git(repo, "tag", "0.1.0")
    _data(
        package,
        {
            "data/apixaban/Kreutz2017/Kreutz2017_Fig1.tsv": "changed\n",
            "data/apixaban/Upreti2013/Upreti2013_Fig2.tsv": "d\n",
            "data/apixaban/collect.py": "y\n",
        },
    )
    (package / "data/edoxaban/Mendell2011/Mendell2011_Tab1.tsv").unlink()
    model = MODEL.replace(
        root=package,
        variants={"male": package / "models" / "apixaban.xml"},
        data_path=(package / "data" / "apixaban", package / "data" / "edoxaban"),
    )
    assert data_dirs(model) == [
        PurePosixPath("data/apixaban"),
        PurePosixPath("data/edoxaban"),
    ]
    changes = data_changes(
        ReleaseRepo(path=repo, last_tag="0.1.0", ahead=0),
        "src/pkdb_models/models/apixaban",
        data_dirs(model),
        study=HOOK.data_study,
    )
    assert changes.added == ["apixaban/Upreti2013"]
    assert changes.removed == ["edoxaban/Mendell2011"]
    assert changes.changed == ["apixaban/Kreutz2017"]


# ---------------------------------------------------------------------------
# the page of the credibility
# ---------------------------------------------------------------------------
@pytest.fixture
def release(tmp_path: Path) -> tuple[ProjectContext, Path]:
    """A release repository of the toy model, with its credibility declared."""
    model = MODEL.replace(credibility="pkpdsim_toy.credibility:CREDIBILITY")
    assert model.release is not None
    context = ProjectContext(
        project=model,
        reference="pkpdsim_toy:MODEL",
        settings=model.release,
        version="0.10.0",
        python="3.14",
        pythons=("3.14", "3.15"),
        dependencies=(),
        constraints=(),
        framework_path=None,
        year=2026,
        requires_python=">=3.14",
        hook=model.release_hook(),
        built_steps=(),
    )
    return context, tmp_path / "toy-model"


def _page(
    release: tuple[ProjectContext, Path], assessment: CredibilityAssessment
) -> str:
    """Write the docs of the release with the assessment, get the credibility page."""
    context, path = release
    (path / BuildPaths.CREDIBILITY).mkdir(parents=True)
    assessment.to_json(path / BuildPaths.CREDIBILITY / "credibility.json")
    write_docs(context, path)
    return (path / "docs" / "credibility.md").read_text(encoding="utf-8")


def _section(page: str, start: str, end: str | None = None) -> str:
    return page[page.index(start) : page.index(end) if end else None]


def test_the_credibility_page_shows_the_ratio_pairs(
    release: tuple[ProjectContext, Path],
) -> None:
    from test_credibility_ratios import _renal

    fits = fit_table(factors={("CredA", "fm_renal"): 1.1}, times=True)
    page = _page(release, assess(_renal(), fits))
    renal = _section(page, "## renal")
    assert "| `A_ratio` | the concentration in renal impairment" in renal
    assert "| ratio | covered | 1 |" in renal
    # the rule of the ratio pairs is the last row of the acceptance table
    acceptance = _section(renal, "| Acceptance criterion |", "| Ratio |")
    assert acceptance.rstrip().endswith("| `Guest` (ratio pairs) | 100 | >= 80 | yes |")
    assert acceptance.count("(ratio pairs)") == 1
    assert (
        "| `A_ratio` | CredA | `fm_renal` / `fm_dose` | 1 | 1.1 | 1.25 | yes |" in renal
    )


def test_the_credibility_page_shows_the_temporal_validation(
    release: tuple[ProjectContext, Path],
) -> None:
    from test_credibility_temporal import FIT_DAY, _chronology

    assessment = assess(
        CREDIBILITY,
        fit_table(),
        prepared={"temporal": _chronology(FIT_DAY + timedelta(days=1))},
    )
    page = _page(release, assessment)
    dose = _section(page, "## dose", "## renal")
    assert "Temporal validation, the decisive studies the model collected" in dose
    assert "| `GMFE` | 1 | <= 1.5 | yes |" in dose
    assert "| CredA | `DECAY` | 2026-10-01 | 2026-10-02 | 2026-09-01 |" in dose
    # after the acceptance table, before the parameter uncertainty and the
    # verification, as in pkpdsim 0.19
    assert dose.index("| Acceptance criterion |") < dose.index("Temporal validation")
    # the rules of the temporal validation are no rows of the acceptance table
    acceptance = _section(dose, "| Acceptance criterion |", "Temporal validation")
    assert "(temporal)" not in acceptance
    describe = _section(page, "## describe", "## dose")
    assert "Temporal validation" not in describe


def test_the_credibility_page_shows_the_study_data(
    release: tuple[ProjectContext, Path], tmp_path: Path
) -> None:
    from test_credibility_studies import _model

    from pkpdsim.credibility.studies import study_data

    assessment = assess(
        CREDIBILITY,
        fit_table(),
        prepared={"studies": study_data(_model(tmp_path / "data"), fit_table().table)},
    )
    page = _page(release, assessment)
    section = _section(page, "## Study data")
    assert (
        "| CredA | 25 | 2026-09-01 (git) | 2026-09-15 | 2026-09-20 `PKDB09999` | "
        "`abcdef012345` | 1 as collected |"
    ) in section
    assert "| CredB | 10 | no data folder of the study | | | | |" in section
    # after the sections of every context of use, before the battery
    assert page.index("## Study data") > page.index("## renal")


def test_the_credibility_page_shows_the_inputs(
    release: tuple[ProjectContext, Path],
) -> None:
    from test_credibility_inputs import _check

    checks = [
        _check("fm_dose", "BW", 72.0, "equivalent"),
        _check("fm_dose", "PODOSE_a", 25.0, "differs"),
    ]
    page = _page(release, assess(CREDIBILITY, fit_table(), prepared={"inputs": checks}))
    dose = _section(page, "## dose", "## renal")
    assert (
        "Inputs, the simulated covariates of the subjects and the doses of 1 of" in dose
    )
    assert "| `BW` | 1 | 72 kg | 1 | 0 | 0 | 0 |" in dose
    assert "| `PODOSE_a` (dose) | 1 | 25 mg | 0 | 0 | 1 | 0 |" in dose
    assert "| CredA/fm_dose | `PODOSE_a` | 25 mg | 70 mg | group all |" in dose
    # between the parameter uncertainty and the variability, as in pkpdsim 0.19
    assert dose.index("Inputs, ") > dose.index("| Acceptance criterion |")


def test_the_inputs_come_after_the_uncertainty_before_the_variability(
    release: tuple[ProjectContext, Path],
) -> None:
    from credibility_tables import predictive_checks_of, propagation_of
    from test_credibility_inputs import _check

    checks = [_check("fm_dose", "BW", 72.0, "equivalent")]
    assessment = assess(
        CREDIBILITY,
        fit_table(),
        prepared={"inputs": checks},
        propagation=propagation_of(0.01),
        variability=predictive_checks_of(0.5),
    )
    dose = _section(_page(release, assessment), "## dose", "## renal")
    uncertainty = dose.index("Parameter uncertainty, ")
    inputs = dose.index("Inputs, ")
    variability = dose.index("Variability, ")
    assert uncertainty < inputs < variability


def test_the_credibility_page_without_the_evidence_of_pkpdsim(
    release: tuple[ProjectContext, Path],
) -> None:
    page = _page(release, assess(CREDIBILITY, fit_table()))
    for text in ("| Ratio |", "Temporal validation", "Inputs, ", "## Study data"):
        assert text not in page
    # the parts of pkpdsim leave no comment and no blank lines of their own
    assert "pkpdsim:" not in page
    assert "\n\n\n" not in page
