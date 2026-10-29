"""The base experiment of the models built on the library (`pkpdsim.pbpk.experiment`)."""

import os
import pickle
import shutil
import subprocess
import sys
from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING, cast

import pytest
from sbmlsim import Q
from sbmlsim.simulation import Simulation
from sbmlsim.simulator import Simulator

from pkpdsim import impairment
from pkpdsim.model import PKPDModel
from pkpdsim.pbpk import Sex, build, load
from pkpdsim.pbpk import experiment as pbpk_experiment
from pkpdsim.pbpk.builder import BuildResult
from pkpdsim.pbpk.experiment import (
    PBPKSimulationExperiment,
    ScaledRenalMap,
    healthy_gfr,
)
from pkpdsim.scan import start_values

from .conftest import definition, individual

if TYPE_CHECKING:
    from sbmlsim.units import Quantity


def _model(path: Path, model_path: Path) -> PKPDModel:
    """A model object around a flat model of the library."""
    return PKPDModel(
        name="toy_pbpk",
        package="toy_pbpk",
        title="Toy PBPK model",
        path=path,
        model_path=model_path,
        factory=".factory:create_models",
    )


@pytest.fixture(scope="module")
def male(tmp_path_factory: pytest.TempPathFactory) -> BuildResult:
    """The simple definition built male, once."""
    return build(definition(), tmp_path_factory.mktemp("experiment_male"))


@pytest.fixture(scope="module")
def female(tmp_path_factory: pytest.TempPathFactory) -> BuildResult:
    """The simple definition built female, once."""
    return build(
        definition(sex=Sex.FEMALE), tmp_path_factory.mktemp("experiment_female")
    )


def _experiment(result: BuildResult) -> type[PBPKSimulationExperiment]:
    """The base of the experiments of a toy model of the library."""

    class ToyExperiment(PBPKSimulationExperiment):
        pkpd_model = _model(result.flat_path.parent, result.flat_path)

    return ToyExperiment


def _subclass(
    base: type[PBPKSimulationExperiment], name: str, **attributes: object
) -> type[PBPKSimulationExperiment]:
    """A subclass of a base of experiments, e.g. a study or another model."""
    return cast("type[PBPKSimulationExperiment]", type(name, (base,), attributes))


def _ml(gfr: "Quantity") -> float:
    """Get a GFR in ml/min, the unit of `healthy_gfr`."""
    assert str(gfr.units) == "milliliter / minute"
    return float(gfr.m_as("ml/min"))


def test_the_healthy_gfr_is_the_gfr_base(male: BuildResult) -> None:
    """The healthy GFR is `GFR_base` [l/min] of the flat model in ml/min."""
    path = male.flat_path
    gfr = {bw: _ml(healthy_gfr(path, bw=bw)) for bw in (50.0, 100.0)}
    for bw, value in gfr.items():
        expected = individual(path, BW=bw)["GFR_base"] * 1000
        assert value == pytest.approx(expected, rel=1e-12)
    assert gfr[50.0] < gfr[100.0]
    assert _ml(healthy_gfr(path)) == pytest.approx(
        individual(path)["GFR_base"] * 1000, rel=1e-12
    )
    assert _ml(healthy_gfr(path, bw=60.0, height=165.0)) == pytest.approx(
        individual(path, BW=60.0, HEIGHT=165.0)["GFR_base"] * 1000, rel=1e-12
    )
    assert _ml(healthy_gfr(path, hf_severity=3.0)) == pytest.approx(
        individual(path, HF_severity=3.0)["GFR_base"] * 1000, rel=1e-12
    )
    assert _ml(healthy_gfr(path, hf_severity=3.0)) < _ml(healthy_gfr(path))
    cls = _experiment(male)
    assert _ml(cls.healthy_gfr(bw=50.0)) == gfr[50.0]
    assert _ml(cls.healthy_gfr(hf_severity=3.0)) == _ml(
        healthy_gfr(path, hf_severity=3.0)
    )


def test_the_healthy_gfr_is_cached(
    male: BuildResult, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An individual is simulated once, however its arguments are given."""
    path = male.flat_path
    _ml(healthy_gfr(path, bw=70.0))
    runs: list[object] = []
    run = Simulator.run

    def counted(self: Simulator, *args: object, **kwargs: object) -> object:
        runs.append(args)
        return run(self, *args, **kwargs)  # ty: ignore[invalid-argument-type]

    monkeypatch.setattr(Simulator, "run", counted)
    _ml(healthy_gfr(path, bw=70.0))
    _ml(healthy_gfr(Path(str(path)), bw=70.0, height=None, hf_severity=0.0))
    assert runs == []


def test_the_healthy_gfr_falls_with_the_cirrhosis(male: BuildResult) -> None:
    """At `CP_severity` 3 the healthy GFR is lower by the hepatorenal factor."""
    path = male.flat_path
    factor = load("impairment")["child_pugh"]["gfr"]["value"][2]
    cirrhosis = _ml(healthy_gfr(path, cp_severity=3.0))
    assert cirrhosis == pytest.approx(
        individual(path, CP_severity=3.0)["GFR_base"] * 1000, rel=1e-12
    )
    assert cirrhosis == pytest.approx(factor * _ml(healthy_gfr(path)), rel=1e-12)
    cls = _experiment(male)
    assert _ml(cls.healthy_gfr(cp_severity=3.0)) == cirrhosis
    f = cls.renal_function(30.0, bw=60.0, height=165.0, cp_severity=3.0)
    result = Simulator(n_workers=1).run(
        path,
        Simulation(
            end=1,
            time_unit="min",
            preinit_changes={
                "BW": Q(60.0, "kg"),
                "HEIGHT": Q(165.0, "cm"),
                "CP_severity": Q(3.0, "dimensionless"),
                "f_renal_function": f,
            },
        ),
    )
    gfr = float(start_values(result, "GFR").item())
    assert Q(gfr, result.units["GFR"]).to("ml/min").magnitude == pytest.approx(
        30.0, rel=1e-9
    )


def test_a_rebuilt_model_is_simulated_again(
    male: BuildResult, female: BuildResult, tmp_path: Path
) -> None:
    """The cache knows the modification time of the model file."""
    path = tmp_path / "model_flat.xml"
    shutil.copyfile(male.flat_path, path)
    assert _ml(healthy_gfr(path)) == pytest.approx(_ml(healthy_gfr(male.flat_path)))
    mtime = path.stat().st_mtime_ns
    shutil.copyfile(female.flat_path, path)
    os.utime(path, ns=(mtime + 1_000_000_000, mtime + 1_000_000_000))
    assert _ml(healthy_gfr(path)) == pytest.approx(_ml(healthy_gfr(female.flat_path)))


def test_the_renal_function_gives_the_measured_clearance(male: BuildResult) -> None:
    """The renal function of a clearance gives that GFR in a simulation."""
    cls = _experiment(male)
    f = cls.renal_function(40.0, bw=60.0, height=165.0)
    assert str(f.units) == "dimensionless"
    result = Simulator(n_workers=1).run(
        male.flat_path,
        Simulation(
            end=1,
            time_unit="min",
            preinit_changes={
                "BW": Q(60.0, "kg"),
                "HEIGHT": Q(165.0, "cm"),
                "f_renal_function": f,
            },
        ),
    )
    gfr = float(start_values(result, "GFR").item())
    assert Q(gfr, result.units["GFR"]).to("ml/min").magnitude == pytest.approx(
        40.0, rel=1e-9
    )
    heart_failure = cls.renal_function(40.0, hf_severity=3.0)
    assert heart_failure.magnitude == pytest.approx(
        40.0 / _ml(healthy_gfr(male.flat_path, hf_severity=3.0)), rel=1e-12
    )


def test_a_renal_function_outside_the_range_is_refused(male: BuildResult) -> None:
    """A clearance whose renal function leaves the range raises a ValueError.

    The range is that of the validity of the kidney outflow by default
    (`renal_function_limits`), a model narrows it with `renal_function_range`.
    """
    cls = _experiment(male)
    assert cls.renal_function_range is None
    assert cls.renal_function(0.0).magnitude == 0.0
    healthy = _ml(cls.healthy_gfr())
    assert cls.renal_function(2.1 * healthy).magnitude == pytest.approx(2.1)
    with pytest.raises(ValueError, match=r"outside of the supported range 0\.0 to inf"):
        cls.renal_function(-1.0)

    no_anuria = _subclass(cls, "NoAnuria", renal_function_range=(0.05, 2.0))

    with pytest.raises(ValueError, match=r"outside of the supported range 0\.05"):
        no_anuria.renal_function(0.0)
    with pytest.raises(ValueError, match=r"supported range 0\.05 to 2\.0"):
        no_anuria.renal_function(2.1 * healthy)
    assert no_anuria.renal_function_limits() == (0.05, 2.0)


def test_the_renal_map_is_lazy(
    male: BuildResult, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Defining a class loads no model, the first value does, once."""
    calls: list[Path] = []
    original = pbpk_experiment.healthy_gfr

    def counted(model_path: Path, **kwargs: float | None) -> "Quantity":
        calls.append(model_path)
        return original(model_path, **kwargs)  # ty: ignore[invalid-argument-type]

    monkeypatch.setattr(pbpk_experiment, "healthy_gfr", counted)

    class Missing(PBPKSimulationExperiment):
        pkpd_model = _model(tmp_path, tmp_path / "missing_flat.xml")

    assert list(Missing.renal_map) == list(impairment.RENAL_FUNCTION_GFR)
    assert len(Missing.renal_map) == len(impairment.RENAL_FUNCTION_GFR)
    assert "Normal renal function" in Missing.renal_map
    assert calls == []

    cls = _experiment(male)
    assert list(cls.renal_map) == list(impairment.RENAL_FUNCTION_GFR)
    assert calls == []
    value = cls.renal_map["Normal renal function"]
    assert calls == [male.flat_path]
    assert cls.renal_map["Mild renal impairment"] < value
    assert dict(cls.renal_map)
    assert calls == [male.flat_path]


def test_importing_the_library_imports_no_simulator() -> None:
    """`import pkpdsim.pbpk` stays light, the base experiment is exported lazily."""
    code = (
        "import sys, pkpdsim.pbpk; "
        "print(sorted(m for m in sys.modules "
        "if m.split('.')[0] in {'sbmlsim', 'roadrunner'} "
        "or m in {'pkpdsim.experiment', 'pkpdsim.pbpk.experiment'})); "
        "from pkpdsim.pbpk import REFERENCE_GFR, PBPKSimulationExperiment, "
        "ScaledRenalMap, healthy_gfr; "
        "from pkpdsim.pbpk import experiment; "
        "print(PBPKSimulationExperiment is experiment.PBPKSimulationExperiment "
        "and ScaledRenalMap is experiment.ScaledRenalMap "
        "and healthy_gfr is experiment.healthy_gfr "
        "and REFERENCE_GFR == experiment.REFERENCE_GFR == 101.0)"
    )
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    assert result.stdout.split() == ["[]", "True"]
    import pkpdsim.pbpk

    lazy = {
        "REFERENCE_GFR",
        "PBPKSimulationExperiment",
        "ScaledRenalMap",
        "healthy_gfr",
    }
    assert lazy <= set(pkpdsim.pbpk.__all__)
    assert lazy <= set(dir(pkpdsim.pbpk))
    assert {"build", "Substance", "ids"} <= set(dir(pkpdsim.pbpk))
    with pytest.raises(AttributeError, match="no_such_name"):
        _ = pkpdsim.pbpk.no_such_name


def test_the_renal_classes_keep_their_gfr(
    male: BuildResult, female: BuildResult
) -> None:
    """A class has the GFR of its preset in ml/min in every model, per class."""
    men = _experiment(male)
    women = _experiment(female)
    assert _ml(men.healthy_gfr()) != pytest.approx(_ml(women.healthy_gfr()), rel=1e-3)
    for cls in (men, women):
        gfr = _ml(cls.healthy_gfr())
        assert list(cls.renal_map) == list(impairment.RENAL_FUNCTION_GFR)
        for key, preset in impairment.RENAL_FUNCTION_GFR.items():
            assert cls.renal_map[key] * gfr == pytest.approx(101.0 * preset, rel=1e-12)
    assert men.renal_map["Normal renal function"] != pytest.approx(
        women.renal_map["Normal renal function"], rel=1e-3
    )

    study = _subclass(men, "Study")
    assert dict(study.renal_map) == dict(men.renal_map)
    assert study.renal_map is not men.renal_map

    midpoints = _subclass(
        men, "Midpoints", renal_map=ScaledRenalMap(impairment.RENAL_FUNCTION)
    )

    gfr = _ml(men.healthy_gfr())
    for key, preset in impairment.RENAL_FUNCTION.items():
        assert midpoints.renal_map[key] * gfr == pytest.approx(101.0 * preset)
    assert men.renal_map is not midpoints.renal_map


def test_the_maps_of_the_library(male: BuildResult) -> None:
    """The base has the maps onto the severities of the library."""
    cls = _experiment(male)
    assert cls.renal_colors is impairment.RENAL_COLORS_END_STAGE
    assert cls.cirrhosis_map is impairment.CHILD_PUGH
    assert cls.cirrhosis_colors is impairment.CIRRHOSIS_COLORS
    assert cls.cardiac_map is impairment.HEART_FAILURE
    assert cls.cardiac_colors is impairment.HEART_FAILURE_COLORS
    assert set(cls.renal_colors) == set(cls.renal_map)


def test_a_scaled_renal_map_of_a_function() -> None:
    """A map with its own GFR function needs no class."""
    gfr: Callable[[], Quantity] = lambda: Q(120.0, "ml/min")  # noqa: E731
    classes = ScaledRenalMap(impairment.RENAL_FUNCTION, gfr)
    for key, preset in impairment.RENAL_FUNCTION.items():
        assert classes[key] == pytest.approx(preset * 101.0 / 120.0, rel=1e-15)
    relative_to_90 = ScaledRenalMap(impairment.RENAL_FUNCTION_90, gfr, reference=90.0)
    assert relative_to_90["Normal renal function"] == pytest.approx(0.75)

    class Holder:
        renal_map = classes

    assert Holder.renal_map is classes
    unbound = ScaledRenalMap(impairment.RENAL_FUNCTION)
    assert list(unbound) == list(impairment.RENAL_FUNCTION)
    with pytest.raises(TypeError, match="healthy GFR"):
        _ = unbound["Normal renal function"]


def test_a_renal_map_is_picklable(male: BuildResult) -> None:
    """A map pickles without the maps of the classes which accessed it."""
    cls = _experiment(male)
    _ = cls.renal_map["Normal renal function"]
    unbound = PBPKSimulationExperiment.__dict__["renal_map"]
    copy = pickle.loads(pickle.dumps(unbound))
    assert isinstance(copy, ScaledRenalMap)
    assert list(copy) == list(impairment.RENAL_FUNCTION_GFR)
    assert copy.__dict__["_bound"] == {}
    with pytest.raises(TypeError, match="healthy GFR"):
        _ = copy["Normal renal function"]


def test_the_healthy_gfr_needs_a_model_object() -> None:
    """The base itself has no model object."""
    with pytest.raises(TypeError, match="has no model object"):
        _ml(PBPKSimulationExperiment.healthy_gfr())
