"""The browser examples are genuine, current, conservative IDA trajectories."""

import hashlib
import itertools
import json
import re
import runpy
import shutil
import subprocess
from pathlib import Path
from typing import Any
from xml.etree import ElementTree

import numpy as np
import pytest
from sbmlsim.model import AbstractModel
from sbmlsim.simulation import Simulation
from sbmlsim.simulator import Simulator

from pkpdsim.pbpk import (
    DynamicHeart,
    PBPKDefinition,
    PressureFlow,
    Sex,
    Variants,
    build,
)

SCRIPT = runpy.run_path(
    str(Path(__file__).parents[2] / "scripts/cardiovascular_examples.py")
)
EXAMPLES = SCRIPT["EXAMPLES"]
INTEGRATOR = SCRIPT["INTEGRATOR"]
OUTPUT = SCRIPT["OUTPUT"]


def test_javascript_network_and_app_interactions() -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node.js is required to test the documentation JavaScript.")
    for filename in ("cardiovascular_model.test.mjs", "cardiovascular_app.test.mjs"):
        result = subprocess.run(
            [node, str(Path(__file__).with_name(filename))],
            capture_output=True,
            text=True,
            timeout=50,
            check=False,
        )
        assert result.returncode == 0, result.stdout + result.stderr


@pytest.fixture(scope="module", params=["male", "female"])
def sex(request: pytest.FixtureRequest) -> Sex:
    return Sex(request.param)


@pytest.fixture(scope="module")
def data(sex: Sex) -> dict[str, Any]:
    path = OUTPUT if sex is Sex.MALE else OUTPUT.parent / "female" / OUTPUT.name
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def sbml(tmp_path_factory: pytest.TempPathFactory, sex: Sex) -> dict[str, Path]:
    folder = tmp_path_factory.mktemp("cardiovascular_examples")
    return {
        name: build(
            PBPKDefinition(
                f"{name}_examples", [], sex=sex, variants=Variants(circulation=variant)
            ),
            folder / name,
        ).flat_path
        for name, variant in (("pressure", PressureFlow()), ("heart", DynamicHeart()))
    }


def test_data_names_its_solver_units_and_reproducible_source(
    data: dict[str, Any], sbml: dict[str, Path]
) -> None:
    assert data["schema_version"] == 1
    provenance = data["provenance"]
    assert provenance["solver"] == "SUNDIALS IDA"
    assert provenance["engine"] == "core"
    assert provenance["generator"] == "scripts/cardiovascular_examples.py"
    assert provenance["units"]["time"] == "s"
    assert provenance["units"]["CV_V_lv"] == "ml"
    fingerprint = hashlib.sha256(
        sbml["pressure"].read_text(encoding="utf-8").encode("utf-8")
        + sbml["heart"].read_text(encoding="utf-8").encode("utf-8")
    ).hexdigest()
    assert provenance["sbml_sha256"] == fingerprint, (
        "Reference models changed: run `uv run python scripts/cardiovascular_examples.py`"
    )


def test_pressure_calibration_matches_current_consistent_sbml(
    data: dict[str, Any], sbml: dict[str, Path]
) -> None:
    pressure = data["pressure"]
    selections = ["CO", "CV_P_po", "CV_P_li"]
    for organ in pressure["organs"]:
        selections += [f"CV_R_{organ['id']}", f"Q{organ['id']}"]
    result = Simulator(engine="core", n_workers=1).simulate(
        AbstractModel(source=sbml["pressure"], selections=selections),
        Simulation(end=1e-9, times=[0.0]),
    )
    assert pressure["cardiac_output"] == pytest.approx(result["CO"][0])
    assert pressure["portal_pressure"] == pytest.approx(result["CV_P_po"][0])
    assert pressure["hepatic_pressure"] == pytest.approx(result["CV_P_li"][0])
    for organ in pressure["organs"]:
        assert organ["resistance"] == pytest.approx(result[f"CV_R_{organ['id']}"][0])
        assert organ["flow"] == pytest.approx(result[f"Q{organ['id']}"][0])
    portal_flow = sum(organ["flow"] for organ in pressure["organs"] if organ["portal"])
    artery_flow = next(
        organ["flow"] for organ in pressure["organs"] if organ["id"] == "li"
    )
    assert portal_flow == pytest.approx(
        (pressure["portal_pressure"] - pressure["hepatic_pressure"])
        / pressure["portal_resistance"]
    )
    assert portal_flow + artery_flow == pytest.approx(
        (pressure["hepatic_pressure"] - pressure["venous_pressure"])
        / pressure["hepatic_resistance"]
    )


def test_live_javascript_perturbations_match_ida(
    data: dict[str, Any], sbml: dict[str, Path]
) -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node.js is required to test the documentation JavaScript.")
    defaults = {
        "hepatic": 1,
        "systemic": 1,
        "portal": 1,
        "renal": 1,
        "arterial": 95,
        "venous": 5,
    }
    controls = [
        defaults | changes
        for changes in (
            {"hepatic": 4},
            {"systemic": 0.7, "portal": 3, "renal": 2},
            {"hepatic": 8, "arterial": 85, "venous": 12},
        )
    ]
    module = Path(__file__).parents[2] / "docs/pbpk/examples/cardiovascular-model.mjs"
    result = subprocess.run(
        [
            node,
            "--input-type=module",
            "-e",
            "const {solvePressure} = await import(process.argv[1]); "
            "const fs = await import('node:fs'); "
            "const {reference, controls} = JSON.parse(fs.readFileSync(0, 'utf8')); "
            "console.log(JSON.stringify(controls.map(c => solvePressure(reference, c))));",
            module.as_uri(),
        ],
        input=json.dumps({"reference": data["pressure"], "controls": controls}),
        capture_output=True,
        text=True,
        timeout=30,
        check=True,
    )
    browser_values = json.loads(result.stdout)
    simulator = Simulator(engine="core", n_workers=1, **INTEGRATOR)
    quantities = {
        "ppo": "CV_P_po",
        "pli": "CV_P_li",
        "qpo": "Qpo",
        "qh": "Qh",
        "qha": "Qha",
        "co": "CO",
    }
    for settings, values in zip(controls, browser_values, strict=True):
        actual = simulator.simulate(
            AbstractModel(
                source=sbml["pressure"], selections=list(quantities.values())
            ),
            Simulation(
                end=1e-9,
                times=[0.0],
                preinit_changes={
                    "CV_fR_hepatic": settings["hepatic"],
                    "CV_fR_systemic": settings["systemic"],
                    "CV_fR_portal": settings["portal"],
                    "CV_fR_ki": settings["renal"],
                    "CV_P_ar": settings["arterial"],
                    "CV_P_ve": settings["venous"],
                },
            ),
        )
        for key, sid in quantities.items():
            assert values[key] == pytest.approx(float(actual[sid][0]), rel=1e-7)


def test_every_case_is_finite_and_preserves_blood_volume(data: dict[str, Any]) -> None:
    assert [case["id"] for case in data["cases"]] == [
        example.sid for example in EXAMPLES
    ]
    for case in data["cases"]:
        assert case["summary"]["maximum_blood_volume_error"] < 1e-8
        for table in (case["series"], case["beats"]):
            values = np.array(table["values"])
            assert values.shape[1] == len(table["columns"])
            assert np.isfinite(values).all()
            assert (np.diff(values[:, 0]) >= 0).all()
            volume = values[:, table["columns"].index("CV_V_blood_total")]
            expected_volume = 5.3 if data["provenance"]["sex"] == "male" else 3.9
            np.testing.assert_allclose(volume, expected_volume, atol=1e-5)
        duration = case["beats"]["values"][-1][0] - case["beats"]["values"][0][0]
        assert duration == pytest.approx(3 * 60 / case["heart_rate"])


def test_increased_hepatic_resistance_raises_portal_pressure(
    data: dict[str, Any],
) -> None:
    cases = {case["id"]: case for case in data["cases"]}
    ordered = [
        cases[key]["summary"]["after"]
        for key in ("baseline", "hepatic-2", "hepatic-4", "hepatic-8")
    ]
    assert all(
        first["mean_portal_pressure"] < second["mean_portal_pressure"]
        for first, second in itertools.pairwise(ordered)
    )
    assert all(
        first["mean_hepatic_flow"] > second["mean_hepatic_flow"]
        for first, second in itertools.pairwise(ordered)
    )


def test_case_summaries_use_complete_beats_and_mean_pump_flow(
    data: dict[str, Any],
) -> None:
    assert data["provenance"]["summary_time_step_seconds"] <= 0.0002
    for case in data["cases"]:
        after = case["summary"]["after"]
        assert after["stroke_volume"] == pytest.approx(
            1000 * after["mean_cardiac_output"] / case["heart_rate"]
        )
        assert after["arterial_systolic"] > after["arterial_diastolic"]
        assert after["mean_cardiac_output"] > 0
    cases = {case["id"]: case for case in data["cases"]}
    baseline = cases["baseline"]["summary"]["after"]
    assert baseline["mean_cardiac_output"] == pytest.approx(
        baseline["mean_lung_flow"], abs=1e-4
    )
    assert (
        cases["contractility-low"]["summary"]["after"]["mean_cardiac_output"]
        < cases["baseline"]["summary"]["after"]["mean_cardiac_output"]
    )


def test_preview_samples_match_a_fresh_ida_integration(
    data: dict[str, Any], sbml: dict[str, Path]
) -> None:
    table = data["cases"][0]["series"]
    reference = np.array(table["values"])
    reference = reference[reference[:, 0] <= 1.2]
    columns = table["columns"]
    result = Simulator(engine="core", n_workers=1, **INTEGRATOR).simulate(
        AbstractModel(source=sbml["heart"], selections=columns),
        Simulation(
            end=float(reference[-1, 0]) / 60, times=(reference[:, 0] / 60).tolist()
        ),
    )
    for key in ("CV_P_ar", "CV_P_po", "CV_P_li", "CV_P_lv", "CV_V_lv", "CO"):
        factor = 1000 if key == "CV_V_lv" else 1
        np.testing.assert_allclose(
            factor * result[key],
            reference[:, columns.index(key)],
            rtol=1e-5,
            atol=1e-4,
        )


def test_anatomy_reuses_the_physiology_figure_bodies(tmp_path: Path) -> None:
    output = tmp_path / "anatomy.json"
    SCRIPT["export_anatomy"](output)
    source = OUTPUT.parent / "cardiovascular-anatomy.json"
    regenerated = json.loads(output.read_text(encoding="utf-8"))
    exported = json.loads(source.read_text(encoding="utf-8"))
    assert regenerated.keys() == exported.keys()
    for key in regenerated.keys() - {"bodies"}:
        assert regenerated[key] == exported[key]
    bodies = regenerated["bodies"]
    assert bodies.keys() == exported["bodies"].keys()
    assert set(bodies) == {"male", "female"}
    assert bodies["male"]["height"] == pytest.approx(176)
    assert bodies["female"]["height"] == pytest.approx(163)
    for sex, body in bodies.items():
        saved = exported["bodies"][sex]
        assert body["height"] == pytest.approx(saved["height"], abs=1e-8)
        assert body["anchors"].keys() == saved["anchors"].keys()
        for key, anchor in body["anchors"].items():
            np.testing.assert_allclose(anchor, saved["anchors"][key], atol=1e-4, rtol=0)
        assert 'data-k="li"' in body["svg"]
        assert set(body["anchors"]) == {organ.value for organ in SCRIPT["Organ"]}
        actual_paths = list(ElementTree.fromstring(f"<svg>{body['svg']}</svg>"))
        saved_paths = list(ElementTree.fromstring(f"<svg>{saved['svg']}</svg>"))
        assert len(actual_paths) == len(saved_paths)
        for actual, expected in zip(actual_paths, saved_paths, strict=True):
            assert actual.tag == expected.tag
            actual_d = actual.attrib.pop("d")
            expected_d = expected.attrib.pop("d")
            assert actual.attrib == expected.attrib
            numbers = r"-?\d+(?:\.\d+)?"
            assert re.sub(numbers, "#", actual_d) == re.sub(numbers, "#", expected_d)
            # visual.body.path quantizes centimetres to one decimal place (1 mm).
            # Platform rounding can place a half-way coordinate in either bin.
            np.testing.assert_allclose(
                [float(v) for v in re.findall(numbers, actual_d)],
                [float(v) for v in re.findall(numbers, expected_d)],
                atol=0.100001,
                rtol=0,
            )
