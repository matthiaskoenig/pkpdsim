"""Tests of the parameter database: its tables, studies and their integration."""

import math
from collections import defaultdict
from pathlib import Path

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


def test_a_median_with_its_quartiles_takes_their_spread() -> None:
    """The median of a box plot with its quartiles: the SD of the log of the IQR."""
    loc = log_location(record(Summary(n=16, median=100.0, q1=50.0, q3=200.0)))
    assert not isinstance(loc, str)
    sigma = math.log(4) / (2 * 0.6744897501960817)
    assert loc.y == pytest.approx(math.log(100))
    assert loc.se == pytest.approx(math.sqrt(math.pi / 2) * sigma / 4)


def test_groups_below_the_limit_of_quantification_compare_their_means() -> None:
    """P-gp of Drozdzik 2020: the controls have a first quartile of 0, so the means compare.

    The ratio of the means of class C to the controls is the 260 % of the
    text; the ratio of the medians (0.52 of the controls) would be 16-fold.
    """
    study = studies.studies()["Drozdzik2020"]
    found, _ = integration.study_estimate(study, "protein", "ABCB1", CIRRHOSIS_C)
    assert found is not None
    assert found.statistic.startswith("arithmetic mean")
    assert found.ratio == pytest.approx(2.6, abs=0.1)


def test_a_mostly_undetected_group_has_no_ratio() -> None:
    """BCRP of Drozdzik 2020 is below the limit in most livers: no ratio."""
    study = studies.studies()["Drozdzik2020"]
    found, excluded = integration.study_estimate(study, "protein", "ABCG2", CIRRHOSIS_C)
    assert found is None
    assert excluded


@pytest.mark.parametrize(
    ("entity", "group", "text"),
    [
        ("ABCB1", "cirrhosis_C", 2.60),
        ("ABCC2", "cirrhosis_C", 0.30),
        ("SLC10A1", "cirrhosis_C", 0.34),
        ("SLC22A1", "cirrhosis_C", 0.35),
        ("SLCO1B1", "cirrhosis_C", 0.46),
        ("SLCO2B1", "cirrhosis_C", 0.27),
        ("ABCC2", "ald", 0.30),
        ("ABCC2", "aih", 0.18),
        ("ABCB1", "pbc", 3.50),
        ("ABCB1", "psc", 2.87),
        ("ABCB1", "aih", 4.10),
        ("SLC22A1", "hcv", 0.30),
        ("SLC22A1", "pbc", 0.38),
        ("SLC22A7", "ald", 0.26),
    ],
)
def test_the_digitized_means_of_drozdzik2020_are_those_of_its_text(
    entity: str, group: str, text: float
) -> None:
    """The ratio of the digitized means to the controls is the change the text states."""
    rows = {
        r.group.id: r.summary
        for r in studies.records(kind="protein", entity=entity)
        if r.study == "Drozdzik2020"
    }
    mean, control = rows[group].mean, rows["controls"].mean
    assert mean is not None and control is not None
    assert mean / control == pytest.approx(text, abs=0.03 * max(1.0, text))


def test_physiology_records_are_in_the_unit_of_their_quantity() -> None:
    """Every record of a quantity of the subjects is in the unit of the quantity."""
    found = studies.records(kind="physiology")
    assert {r.entity for r in found} == set(studies.QUANTITIES)
    for r in found:
        assert r.unit == studies.QUANTITIES[r.entity].unit
        assert r.tissue == "body"
    via = {s.via for s in studies.studies().values() if s.via}
    assert via == {"Li2015"}


def _write_study(folder: Path, unit: str) -> None:
    """A study of the albumin of one group, its record in the unit given."""
    folder.mkdir()
    folder.joinpath("study.toml").write_text(
        'source = "Viani1992"\ntissue = "body"\nnotes = ""\n\n'
        '[measurements.physiology]\nfile = "physiology.tsv"\nkind = "physiology"\n\n'
        '[groups.cirrhosis_A]\ncondition = "cirrhosis"\nchild_pugh = "A"\n',
        encoding="utf-8",
    )
    folder.joinpath("physiology.tsv").write_text(
        f"entity\tgroup\tn\tmean\tsd\tunit\nalbumin\tcirrhosis_A\t10\t4.0\t0.5\t{unit}\n",
        encoding="utf-8",
    )


def test_a_physiology_record_in_another_unit_is_refused(tmp_path: Path) -> None:
    """A source which prints another unit is converted when its numbers are entered."""
    _write_study(tmp_path / "Viani1992", "g/dL")
    with pytest.raises(ValueError, match="the unit of the quantity is 'g/L'"):
        studies._study(tmp_path / "Viani1992")


def test_a_physiology_record_takes_the_unit_of_its_quantity(tmp_path: Path) -> None:
    """A record without a unit is in the unit of its quantity."""
    _write_study(tmp_path / "Viani1992", "")
    study = studies._study(tmp_path / "Viani1992")
    (r,) = study.measurements["physiology"].records
    assert r.unit == "g/L"
    assert study.control is None


def test_a_level_pools_the_means_and_the_medians_of_the_studies() -> None:
    """The MELD of a class pools the mean of one cohort and the medians of the others."""
    meld = integration.level("meld", {"condition": "cirrhosis", "child_pugh": "B"})
    assert meld is not None
    assert meld.unit == "1"
    by_study = {e.study: e for e in meld.estimates}
    assert set(by_study) == {
        "Hong2025",
        "Kim2020",
        "Mesropyan2022",
        "Xu2026",
        "Yang2020",
    }
    assert by_study["Mesropyan2022"].statistic == "mean and SD"
    assert by_study["Mesropyan2022"].se == pytest.approx(4.2 / math.sqrt(32))
    kim = by_study["Kim2020"]
    assert kim.statistic == "median and quartiles"
    assert kim.value == 13
    assert kim.se == pytest.approx(
        math.sqrt(math.pi / 2) * (15 - 11) / (2 * 0.6744897501960817) / math.sqrt(51)
    )
    assert min(e.value for e in meld.estimates) < meld.value
    assert meld.value < max(e.value for e in meld.estimates)
    assert meld.tau2 > 0
    lo, hi = meld.ci
    assert lo < meld.value < hi
    assert (
        integration.level("meld", {"condition": "cirrhosis", "etiology": "nope"})
        is None
    )


def test_the_physiology_of_cirrhosis_worsens_with_the_class() -> None:
    """Albumin falls and the MELD and the cardiac index rise from class A to C."""
    for quantity, sign in (("albumin", -1), ("meld", 1), ("cardiac_index", 1)):
        values = []
        for cls in "ABC":
            c = integration.level(
                quantity, {"condition": "cirrhosis", "child_pugh": cls}
            )
            assert c is not None
            values.append(sign * c.value)
        assert values == sorted(values), quantity
