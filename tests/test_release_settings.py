"""Tests of the release metadata of a model."""

import pytest

from pkpdsim.model import PKPDModel
from pkpdsim.release.settings import Author, ReleaseSettings
from pkpdsim_toy import MODEL


def make_settings(**changes: object) -> ReleaseSettings:
    values: dict = {
        "repository": "matthiaskoenig/hctz-model",
        "doi": "10.5281/zenodo.15427418",
        "authors": [Author("König, Matthias", orcid="0000-0003-1725-179X")],
        "keywords": ["pharmacokinetics", "PBPK"],
        "docker_image": "matthiaskoenig/hctz",
        "funding": "Funded by someone.",
        "description": "A model.",
        "first_year": 2024,
    }
    values.update(changes)
    return ReleaseSettings(**values)


def test_settings_tuples() -> None:
    settings = make_settings()
    assert isinstance(settings.authors, tuple)
    assert isinstance(settings.keywords, tuple)
    assert isinstance(settings.exclude, tuple)
    assert hash(settings) == hash(make_settings())


def test_repository_format() -> None:
    assert make_settings().repository == "matthiaskoenig/hctz-model"
    for bad in ("hctz-model", "https://github.com/x/y"):
        with pytest.raises(ValueError, match="owner/name"):
            make_settings(repository=bad)


def test_urls() -> None:
    settings = make_settings()
    assert settings.owner == "matthiaskoenig"
    assert settings.repo_name == "hctz-model"
    assert settings.site_url == "https://matthiaskoenig.github.io/hctz-model/"
    assert settings.raw_url == (
        "https://raw.githubusercontent.com/matthiaskoenig/hctz-model/main/"
    )


def test_missing() -> None:
    assert make_settings().missing() == []
    assert make_settings(funding="", keywords=()).missing() == ["funding", "keywords"]
    assert make_settings(exclude=()).missing() == []
    assert make_settings(first_year=0).missing() == ["first_year"]
    assert make_settings(first_year=-1).missing() == ["first_year"]
    assert make_settings(first_year=True).missing() == ["first_year"]


def test_default_exclude() -> None:
    assert make_settings().exclude == ("TODO.md", "literature", "thesis")


def test_citation_name() -> None:
    assert Author("König, Matthias").citation_name == "König, M."


def test_model_release_default_none() -> None:
    assert (
        PKPDModel(
            name="x",
            package="x",
            title="x",
            path=MODEL.path,
            model_path=MODEL.model_path,
            factory="x:y",
        ).release
        is None
    )


def test_toy_has_release() -> None:
    assert MODEL.release is not None
    assert MODEL.release.missing() == []
