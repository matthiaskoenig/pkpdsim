"""The builds per sex of a library model: the archive, the experiments, the GFR.

One definition is built male and female (`create_archive(..., sexes=...)`),
the model object holds both builds and the experiments choose one by the sex
of their subjects. A study of women must run on the female build and take the
renal function of a woman from the female GFR: a female experiment on the male
build, or the GFR of a man for a woman, would pass every other test.
"""

import zipfile
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

import numpy as np
import pytest
from sbmlsim.experiment import ExperimentRunner
from sbmlsim.simulation import Simulation
from sbmlsim.simulator import Simulator

from pkpdsim import impairment
from pkpdsim.model import PKPDModel
from pkpdsim.pbpk import Sex, build, create_archive, ids, sex_definition
from pkpdsim.pbpk import experiment as pbpk_experiment
from pkpdsim.pbpk.experiment import PBPKSimulationExperiment, healthy_gfr

from .conftest import definition, individual

if TYPE_CHECKING:
    from sbmlsim.units import Quantity

COPYRIGHT = "2024-2026 Matthias König"


@pytest.fixture(scope="module")
def archive(tmp_path_factory: pytest.TempPathFactory) -> dict[str, dict[str, Any]]:
    """The archive of the simple definition with a male and a female build."""
    return create_archive(
        definition(),
        tmp_path_factory.mktemp("sexes") / "models",
        model="simple",
        copyright=COPYRIGHT,
        sexes=(Sex.MALE, Sex.FEMALE),
    )


@pytest.fixture(scope="module")
def model(archive: dict[str, dict[str, Any]]) -> PKPDModel:
    """The model object with both builds."""
    male = archive["simple_flat"]["path"]
    return PKPDModel(
        name="simple",
        package="simple",
        title="Simple PBPK model",
        path=male.parent.parent,
        model_path=male,
        female_model_path=archive["simple_female_flat"]["path"],
        factory=".factory:create_models",
    )


def _base(model: PKPDModel) -> type[PBPKSimulationExperiment]:
    """The base of the experiments of the model: the individual of a simulation."""

    class Base(PBPKSimulationExperiment):
        pkpd_model = model
        selections = ("time", "BW", ids.GFR_BASE, ids.GFR)

        def simulations(self) -> dict[str, Simulation]:
            return {"iv": Simulation(end=10)}

    return Base


def _subclass(
    base: type[PBPKSimulationExperiment], name: str, **attributes: object
) -> type[PBPKSimulationExperiment]:
    """A study of a base of experiments."""
    return cast("type[PBPKSimulationExperiment]", type(name, (base,), attributes))


def _run(cls: type[PBPKSimulationExperiment]) -> PBPKSimulationExperiment:
    """Initialize and run an experiment, as a run of the commands does."""
    runner = ExperimentRunner(
        experiment_classes=[cls],
        base_path=cls.pkpd_model.path,
        simulator=Simulator(n_workers=1),
        on_error="raise",
    )
    (experiment,) = runner.experiments.values()
    assert isinstance(experiment, PBPKSimulationExperiment)
    experiment.run(runner.simulator, output_path=None)
    return experiment


def _value(experiment: PBPKSimulationExperiment, task: str, sid: str) -> float:
    """Get the value of a quantity at the start of a task."""
    return float(np.asarray(experiment.results[task].ds[sid]).ravel()[0])


def _ml(gfr: "Quantity") -> float:
    """Get a GFR in ml/min."""
    return float(gfr.m_as("ml/min"))


# ----------------------------------------------------------------------------
# the archive
# ----------------------------------------------------------------------------
def test_the_archive_has_both_builds_with_the_male_one_as_master(
    archive: dict[str, dict[str, Any]],
) -> None:
    assert set(archive) == {
        "cc-by",
        "mit",
        "simple",
        "simple_flat",
        "simple_flat_md",
        "simple_flat_ant",
        "simple_female",
        "simple_female_flat",
        "simple_female_flat_md",
        "simple_female_flat_ant",
    }
    assert [k for k, v in archive.items() if v["entry"].master] == ["simple_flat"]
    female = archive["simple_female_flat"]
    assert female["path"].name == "simple_female_flat.xml"
    assert female["entry"].location == "./models/simple_female_flat.xml"
    for info in archive.values():
        assert info["path"].exists()
    omex = archive["simple_flat"]["path"].parent.parent / "simple_model.omex"
    with zipfile.ZipFile(omex) as z:
        names = set(z.namelist())
    assert {
        "models/simple_flat.xml",
        "models/simple_female.xml",
        "models/simple_female_flat.xml",
        "models/simple_female_flat.md",
        "models/simple_female_flat.ant",
    } <= names


def test_the_builds_have_the_reference_individual_of_their_sex(
    archive: dict[str, dict[str, Any]],
) -> None:
    male = individual(archive["simple_flat"]["path"])
    female = individual(archive["simple_female_flat"]["path"])
    assert male["BW"] == pytest.approx(73.0)
    assert female["BW"] == pytest.approx(60.0)
    assert female[ids.GFR_BASE] < male[ids.GFR_BASE]


def test_the_male_build_of_an_archive_is_the_build_of_the_definition(
    archive: dict[str, dict[str, Any]], tmp_path: Path
) -> None:
    """A second sex leaves the files of the build of the definition as they were."""
    single = create_archive(
        definition(), tmp_path / "models", model="simple", copyright=COPYRIGHT
    )
    assert set(single) == {k for k in archive if "female" not in k}
    for key in ("simple", "simple_flat", "simple_flat_md", "simple_flat_ant"):
        assert single[key]["path"].read_text() == archive[key]["path"].read_text()
    female = build(definition(sex=Sex.FEMALE, sid="simple_female"), tmp_path / "f")
    assert (
        female.flat_path.read_text()
        == archive["simple_female_flat"]["path"].read_text()
    )


@pytest.mark.parametrize(
    ("sexes", "match"),
    [
        ((Sex.FEMALE,), "the sex of the definition"),
        ((Sex.MALE, Sex.FEMALE, Sex.FEMALE), "twice"),
        ((), "the sex of the definition"),
    ],
    ids=["without-the-definition", "twice", "none"],
)
def test_the_sexes_of_an_archive_are_checked(
    tmp_path: Path, sexes: tuple[Sex, ...], match: str
) -> None:
    with pytest.raises(ValueError, match=match):
        create_archive(
            definition(),
            tmp_path / "models",
            model="simple",
            copyright=COPYRIGHT,
            sexes=sexes,
        )
    assert not (tmp_path / "models").exists()


# ----------------------------------------------------------------------------
# the experiments
# ----------------------------------------------------------------------------
def test_a_female_experiment_runs_the_female_model(
    model: PKPDModel, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The study of women simulates the female build and checks its kidney there."""
    checked: list[Path] = []
    outflow = pbpk_experiment.kidney_outflow

    def recorded(path: Path, changes: object = None) -> object:
        checked.append(Path(path))
        return outflow(path, changes)  # ty: ignore[invalid-argument-type]

    monkeypatch.setattr(pbpk_experiment, "kidney_outflow", recorded)
    base = _base(model)
    women = _run(_subclass(base, "Women", sex=Sex.FEMALE))
    female_path = model.model_path_of(Sex.FEMALE)
    assert {k: m.source.path for k, m in women._models.items()} == {
        "model": female_path
    }
    assert checked == [female_path]
    assert _value(women, "task_iv", "BW") == pytest.approx(60.0)
    assert _value(women, "task_iv", ids.GFR_BASE) * 1000 == pytest.approx(
        _ml(healthy_gfr(female_path)), rel=1e-9
    )

    checked.clear()
    men = _run(_subclass(base, "Men"))
    assert checked == [model.model_path]
    assert _value(men, "task_iv", "BW") == pytest.approx(73.0)


def test_a_study_of_both_sexes_runs_and_checks_each_build(
    model: PKPDModel, monkeypatch: pytest.MonkeyPatch
) -> None:
    checked: list[Path] = []
    outflow = pbpk_experiment.kidney_outflow

    def recorded(path: Path, changes: object = None) -> object:
        checked.append(Path(path))
        return outflow(path, changes)  # ty: ignore[invalid-argument-type]

    monkeypatch.setattr(pbpk_experiment, "kidney_outflow", recorded)

    def simulations(self: PBPKSimulationExperiment) -> dict[str, Simulation]:
        return {"men": Simulation(end=10), "women": Simulation(end=10)}

    groups_cls = _subclass(
        _base(model),
        "Groups",
        simulation_sexes={"women": Sex.FEMALE},
        simulations=simulations,
    )
    groups = _run(groups_cls)
    assert sorted(checked) == sorted(model.model_paths.values())
    assert _value(groups, "task_men", "BW") == pytest.approx(73.0)
    assert _value(groups, "task_women", "BW") == pytest.approx(60.0)


def test_the_renal_function_of_a_woman_uses_the_female_gfr(model: PKPDModel) -> None:
    """The healthy GFR, the renal function and the renal classes of a woman."""
    base = _base(model)
    women = _subclass(base, "Women", sex=Sex.FEMALE)
    female_path = model.model_path_of(Sex.FEMALE)
    female = _ml(healthy_gfr(female_path))
    male = _ml(healthy_gfr(model.model_path))
    assert female < male
    assert _ml(women.healthy_gfr()) == female
    assert _ml(base.healthy_gfr()) == male
    assert _ml(base.healthy_gfr(sex=Sex.FEMALE)) == female
    assert _ml(women.healthy_gfr(sex=Sex.MALE)) == male
    assert _ml(women.healthy_gfr(bw=50.0)) == _ml(healthy_gfr(female_path, bw=50.0))

    # the renal function of a measured clearance, the simulated GFR is the clearance
    f = float(women.renal_function(60.0).m_as("dimensionless"))
    assert f == pytest.approx(60.0 / female, rel=1e-12)
    assert float(base.renal_function(60.0, sex=Sex.FEMALE).m_as("dimensionless")) == f
    assert float(base.renal_function(60.0).m_as("dimensionless")) == pytest.approx(
        60.0 / male
    )
    gfr = individual(female_path, f_renal_function=f)[ids.GFR] * 1000
    assert gfr == pytest.approx(60.0, rel=1e-9)

    # the classes of renal impairment keep their GFR in ml/min in women
    for key, preset in impairment.RENAL_FUNCTION_GFR.items():
        assert women.renal_map[key] * female == pytest.approx(101.0 * preset)
        assert base.renal_map[key] * male == pytest.approx(101.0 * preset)
    assert dict(base.renal_map_of(Sex.FEMALE)) == pytest.approx(dict(women.renal_map))
    assert base.renal_map_of(Sex.MALE) is base.renal_map
    assert base.renal_map_of(Sex.FEMALE) is base.renal_map_of(Sex.FEMALE)

    # the range of the renal functions of a woman
    low, high = women.renal_function_limits()
    assert (low, high) == base.renal_function_limits(sex=Sex.FEMALE)
    assert (low, high) != base.renal_function_limits()


def test_the_renal_classes_of_a_sex_need_a_map_of_the_class() -> None:
    """A map with a GFR of its own is of one sex, a plain mapping of none."""
    from sbmlsim import Q

    from pkpdsim.pbpk.experiment import ScaledRenalMap

    own = ScaledRenalMap(impairment.RENAL_FUNCTION, lambda: Q(120.0, "ml/min"))
    with pytest.raises(TypeError, match="its own GFR"):
        own.of_sex(Sex.FEMALE)
    unbound = ScaledRenalMap(impairment.RENAL_FUNCTION)
    with pytest.raises(TypeError, match="no class"):
        unbound.of_sex(Sex.FEMALE)


def test_a_female_study_of_a_model_without_a_female_build_fails(
    tmp_path: Path,
) -> None:
    built = build(definition(), tmp_path)
    male_only = PKPDModel(
        name="simple",
        package="simple",
        title="Simple PBPK model",
        path=tmp_path,
        model_path=built.flat_path,
        factory=".factory:create_models",
    )
    women = _subclass(_base(male_only), "Women", sex=Sex.FEMALE)
    with pytest.raises(ValueError, match="'simple' has no female build"):
        women.healthy_gfr()
    with pytest.raises(ValueError, match="'simple' has no female build"):
        _run(women)


def test_the_definition_of_a_build_of_another_sex() -> None:
    named = definition(name="Simple model")
    assert sex_definition(named, Sex.MALE) is named
    female = sex_definition(named, Sex.FEMALE)
    assert (female.sid, female.sex, female.name) == (
        "simple_female",
        Sex.FEMALE,
        "Simple model (female)",
    )
    assert sex_definition(definition(), Sex.FEMALE).name is None
    assert female.substances == named.substances
