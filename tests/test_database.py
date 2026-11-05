"""Tests of the parameter database: its tables, studies and their integration."""

import math
from collections import defaultdict

import libsbml
import pytest

from pkpdsim.database import integration, studies, tables
from pkpdsim.database.integration import Estimate, log_location, random_effects
from pkpdsim.database.studies import Group, Record, Summary
from pkpdsim.pbpk.data import sbml_uncertainty
from pkpdsim.pbpk.units import U

CIRRHOSIS_C = {"condition": "cirrhosis", "child_pugh": "C"}


def record(
    summary: Summary, *, study: str = "ElKhateeb2021", group: str = "x"
) -> Record:
    """A record of a group with statistics."""
    return Record(
        study=study,
        measurement="protein",
        kind="protein",
        tissue="liver",
        unit="pmol/g",
        entity="CYP3A4",
        name="CYP3A4",
        group=Group(id=group, condition="cirrhosis", child_pugh="C"),
        summary=summary,
    )


def estimate(y: float, se: float, study: str) -> Estimate:
    """An estimate of a study."""
    return Estimate(
        study=study,
        kind="protein",
        groups=("g",),
        control="c",
        log_ratio=y,
        se=se,
        statistic="",
        imputed=False,
        n=10,
        n_control=10,
    )


def test_every_study_is_a_source_and_reads() -> None:
    """Every study folder is a source of sources.toml and its records read."""
    found = studies.studies()
    assert {"Drozdzik2021", "ElKhateeb2021", "Prasad2018"} <= set(found)
    for key, study in found.items():
        assert key in tables.sources()
        assert study.measurements
        for m in study.measurements.values():
            assert m.records
            assert {r.group.id for r in m.records} <= set(study.groups)


def test_every_inconsistent_record_has_a_note() -> None:
    """A record whose statistics disagree says why, the numbers stay as printed."""
    assert studies.check() == []


def test_the_classes_of_drozdzik2021_sum_to_their_groups() -> None:
    """The n of the classes A, B and C of a group is the n of its A+B+C column.

    A check of the extraction of Tables S2 and S3: every cell is in the row and
    column of its group and class. A sum which the source prints otherwise has
    a note at its A+B+C record (the mRNA of CYP2D6 in AIH).
    """
    for m in studies.studies()["Drozdzik2021"].measurements.values():
        n: dict[tuple[str, str], dict[str, int]] = defaultdict(dict)
        notes: dict[tuple[str, str], str] = {}
        for r in m.records:
            if r.group.id == "controls":
                continue
            prefix, cls = r.group.id.rsplit("_", 1)
            assert r.summary.n is not None
            n[(r.entity, prefix)][cls] = r.summary.n
            if cls == "all":
                notes[(r.entity, prefix)] = r.notes
        for key, by_class in n.items():
            parts = sum(v for k, v in by_class.items() if k != "all")
            assert by_class["all"] == parts or "sum to" in notes[key], (m.key, key)


def test_a_median_with_its_interval_is_its_own_location() -> None:
    """The median with its 95 % CI gives the log median and the SE of its interval."""
    loc = log_location(record(Summary(n=10, median=100.0, ci=(50.0, 200.0))))
    assert not isinstance(loc, str)
    assert loc.y == pytest.approx(math.log(100))
    assert loc.se == pytest.approx(math.log(4) / (2 * integration.Z))
    assert loc.statistic == "median and its 95 % CI"


def test_a_median_takes_the_spread_of_mean_and_sd() -> None:
    """The median with the mean and SD: the log-normal SD of the log, the SE of a median."""
    loc = log_location(record(Summary(n=16, mean=120.0, sd=60.0, median=100.0)))
    assert not isinstance(loc, str)
    sigma = math.sqrt(math.log1p(0.25))
    assert loc.y == pytest.approx(math.log(100))
    assert loc.se == pytest.approx(math.sqrt(math.pi / 2) * sigma / 4)
    assert not loc.imputed


def test_a_mean_and_sd_alone_give_the_geometric_mean() -> None:
    """The mean and SD alone: the geometric mean of the log-normal distribution."""
    loc = log_location(record(Summary(n=25, mean=309.0, sd=240.2)))
    assert not isinstance(loc, str)
    sigma2 = math.log1p((240.2 / 309.0) ** 2)
    assert loc.y == pytest.approx(math.log(309.0) - sigma2 / 2)
    assert loc.se == pytest.approx(math.sqrt(sigma2) / 5)


def test_a_group_without_spread_is_imputed_and_without_location_excluded() -> None:
    """A group without a spread takes the one of its study; a median of 0 has no log."""
    loc = log_location(record(Summary(n=6, median=4.7)))
    assert not isinstance(loc, str)
    assert loc.imputed
    assert isinstance(log_location(record(Summary(n=6, median=0.0, mean=0.0))), str)
    assert isinstance(log_location(record(Summary(median=4.7))), str)


def test_random_effects_of_one_study_is_the_study() -> None:
    """A single study has no between-study variance."""
    mu, se, tau2, q, i2 = random_effects([estimate(-1.0, 0.2, "a")])
    assert (mu, se, tau2, q, i2) == (-1.0, pytest.approx(0.2), 0.0, None, None)


def test_random_effects_is_dersimonian_laird() -> None:
    """The pooled estimate, its SE and tau^2 of DerSimonian and Laird."""
    data = [(-1.6, 0.30), (-1.2, 0.25), (-2.1, 0.40)]
    mu, se, tau2, q, i2 = random_effects(
        [estimate(y, s, str(k)) for k, (y, s) in enumerate(data)]
    )
    w = [1 / s**2 for _, s in data]
    fixed = sum(wi * y for wi, (y, _) in zip(w, data, strict=True)) / sum(w)
    q_ref = sum(wi * (y - fixed) ** 2 for wi, (y, _) in zip(w, data, strict=True))
    tau2_ref = max(0.0, (q_ref - 2) / (sum(w) - sum(wi**2 for wi in w) / sum(w)))
    ws = [1 / (s**2 + tau2_ref) for _, s in data]
    assert q == pytest.approx(q_ref)
    assert tau2 == pytest.approx(tau2_ref)
    assert mu == pytest.approx(
        sum(wi * y for wi, (y, _) in zip(ws, data, strict=True)) / sum(ws)
    )
    assert se == pytest.approx(math.sqrt(1 / sum(ws)))
    assert i2 == pytest.approx(max(0.0, (q_ref - 2) / q_ref))


def test_homogeneous_studies_pool_as_fixed_effect() -> None:
    """Without heterogeneity (Q below its degrees of freedom) tau^2 is 0."""
    mu, se, tau2, _, i2 = random_effects(
        [estimate(-1.0, 0.3, "a"), estimate(-1.01, 0.3, "b")]
    )
    assert tau2 == 0.0
    assert i2 == 0.0
    assert mu == pytest.approx(-1.005)
    assert se == pytest.approx(0.3 / math.sqrt(2))


def test_a_study_uses_its_aggregate_group_not_its_strata() -> None:
    """Drozdzik 2021 gives the class C of all diseases: its etiologies are not counted again."""
    study = studies.studies()["Drozdzik2021"]
    found, _ = integration.study_estimate(study, "protein", "CYP3A4", CIRRHOSIS_C)
    assert found is not None
    assert found.groups == ("cirrhosis_C",)
    assert found.n == 28


def test_a_study_without_an_aggregate_pools_its_strata() -> None:
    """Prasad 2018 has alcoholic and HCV cirrhosis, both assumed class C: pooled, one control."""
    study = studies.studies()["Prasad2018"]
    found, _ = integration.study_estimate(study, "protein", "CYP3A4", CIRRHOSIS_C)
    assert found is not None
    assert found.groups == ("ald", "hcv")
    assert found.n == 57
    assert found.n_control == 25
    assert found.assumptions


def test_the_consensus_pools_every_study_of_the_class() -> None:
    """CYP3A4 of class C: the three studies, a factor below 1 within their range."""
    c = integration.consensus("CYP3A4", CIRRHOSIS_C)
    assert c is not None
    assert set(c.sources) == {"Drozdzik2021", "ElKhateeb2021", "Prasad2018"}
    ratios = [e.ratio for e in c.estimates]
    assert min(ratios) <= c.value <= max(ratios)
    lo, hi = c.ci
    assert lo < c.value < hi
    assert c.prediction_se >= c.se
    assert "3 studies" in c.text()


def test_the_factor_prefers_the_protein() -> None:
    """The factor is the protein abundance where it exists, the mRNA only without."""
    protein = integration.factor("CYP3A4", CIRRHOSIS_C)
    assert protein is not None
    assert protein.kind == "protein"
    mrna = integration.factor("CYP3A4", CIRRHOSIS_C, kinds=("mrna",))
    assert mrna is not None
    assert mrna.kind == "mrna"
    assert mrna.sources == ("Drozdzik2021",)
    assert "mRNA" in mrna.text()


def test_an_unknown_kind_is_refused() -> None:
    """The kinds of measurement are a closed set."""
    with pytest.raises(ValueError, match="Unknown kind"):
        integration.consensus("CYP3A4", CIRRHOSIS_C, kind="activity")


def test_a_table_number_reads_its_uncertainty() -> None:
    """A number of a table states its uncertainty, read into `Uncertainty`."""
    found = tables._uncertainty(
        {"value": 1.0, "sd": 0.2, "n": 12, "ci": [0.8, 1.3]}, "x"
    )
    assert found == tables.Uncertainty(sd=0.2, n=12, ci=(0.8, 1.3))
    assert tables._uncertainty({"value": 1.0}, "x") is None
    with pytest.raises(ValueError, match="negative"):
        tables._uncertainty({"value": 1.0, "sd": -1.0}, "x")
    with pytest.raises(ValueError, match="lower, upper"):
        tables._uncertainty({"value": 1.0, "ci": [1.0]}, "x")


def test_an_uncertainty_is_written_as_distrib() -> None:
    """The uncertainty of a number becomes the distrib uncertainty of its parameter."""
    u = sbml_uncertainty(
        tables.Uncertainty(sd=0.2, se=0.05, cv=0.1, n=12, ci=(0.8, 1.3)),
        U.dimensionless,
        "ICRP89",
    )
    types = [p.type for p in u.uncertParameters]
    assert types == [
        libsbml.DISTRIB_UNCERTTYPE_STANDARDDEVIATION,
        libsbml.DISTRIB_UNCERTTYPE_STANDARDERROR,
        libsbml.DISTRIB_UNCERTTYPE_COEFFIENTOFVARIATION,
        libsbml.DISTRIB_UNCERTTYPE_SAMPLESIZE,
        libsbml.DISTRIB_UNCERTTYPE_CONFIDENCEINTERVAL,
    ]


def test_an_entry_names_an_unknown_source() -> None:
    """Every number of a table names a source of sources.toml."""
    e = tables.entry(
        "routes",
        *next(
            (section, key)
            for section, content in tables.load("routes").items()
            for key in content
        ),
    )
    assert e.source in tables.sources()
