"""The dynamic sliders run the same SBML/IDA model in JavaScript and Python."""

import hashlib
import json
import os
import runpy
import shutil
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import sbmlsim_core
from sbmlsim.model import AbstractModel
from sbmlsim.simulation import Change, Simulation
from sbmlsim.simulator import Simulator
from sbmlsim_core import web

ROOT = Path(__file__).resolve().parents[2]
ASSETS = ROOT / "docs/pbpk/examples/solver"
GENERATOR = runpy.run_path(str(ROOT / "scripts/cardiovascular_examples.py"))


@pytest.fixture(scope="module")
def browser_result(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node.js is required for the WebAssembly worker comparison.")
    output = tmp_path_factory.mktemp("heart_worker") / "result.json"
    result = subprocess.run(
        [node, str(Path(__file__).with_name("cardiovascular_solver.test.mjs"))],
        env=dict(os.environ, PKPDSIM_HEART_TEST_OUTPUT=str(output)),
        text=True,
        capture_output=True,
        timeout=60,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    return json.loads(output.read_text(encoding="utf-8"))


def test_browser_runtime_is_exported_from_the_installed_core_release() -> None:
    config = json.loads((ASSETS / "config.json").read_text(encoding="utf-8"))
    description = json.loads((ASSETS / "model.json").read_text(encoding="utf-8"))
    assert description["n_algebraic"] == 1
    assert config["engine"] == "SUNDIALS IDA"
    assert (
        config["sbml_sha256"]
        == hashlib.sha256((ASSETS / "model.xml").read_bytes()).hexdigest()
    )
    for name, wasm in sbmlsim_core.runtime_wasm().items():
        assert (ASSETS / f"{name}.wasm").read_bytes() == wasm
    assert (ASSETS / "sbmlsim-core.js").read_text(
        encoding="utf-8"
    ) == web.script().replace("\r\n", "\n")
    assert (ASSETS / "THIRD_PARTY.md").read_text(
        encoding="utf-8"
    ) == sbmlsim_core.THIRD_PARTY.replace("\r\n", "\n")
    assert (ASSETS / "LICENSE").is_file()


@pytest.mark.parametrize("newline", ["\n", "\r\n"], ids=["LF", "CRLF"])
def test_browser_export_has_portable_text_and_fingerprints(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, newline: str
) -> None:
    source = tmp_path / "heart.xml"
    source.write_bytes(
        (ASSETS / "model.xml")
        .read_text(encoding="utf-8")
        .replace("\n", newline)
        .encode()
    )
    script = web.script().replace("\r\n", "\n").replace("\n", newline)
    monkeypatch.setattr(web, "script", lambda: script)
    monkeypatch.setattr(
        sbmlsim_core,
        "THIRD_PARTY",
        sbmlsim_core.THIRD_PARTY.replace("\r\n", "\n").replace("\n", newline),
    )
    output = tmp_path / "export"
    GENERATOR["export_browser_model"](source, output)
    for name in ("model.xml", "config.json", "sbmlsim-core.js", "THIRD_PARTY.md"):
        assert (output / name).read_bytes() == (ASSETS / name).read_bytes(), name


def test_slider_worker_matches_a_fresh_python_ida_solution(
    browser_result: dict[str, Any],
) -> None:
    changed = browser_result["changed"]
    after = changed["summary"]["after"]
    heart_rate = changed["heart_rate"]
    end = changed["provenance"]["end_time"] / 60
    start = end - 3 / heart_rate
    count = int(np.ceil((end - start) * 60 / 0.0002)) + 1
    times = np.unique(np.concatenate((np.linspace(start, end, count), [start, end])))
    columns = GENERATOR["COLUMNS"]
    native = Simulator(engine="core", n_workers=1, **GENERATOR["INTEGRATOR"]).simulate(
        AbstractModel(source=ASSETS / "model.xml", selections=list(columns)),
        Simulation(
            end=end,
            times=times.tolist(),
            changes=[Change(0.4, changed["changes"])],
        ),
    )
    summary = GENERATOR["_window_summary"](native, start, end, heart_rate)
    for key in (
        "mean_portal_pressure",
        "mean_hepatic_pressure",
        "mean_cardiac_output",
        "mean_lung_flow",
        "stroke_volume",
        "arterial_systolic",
        "arterial_diastolic",
    ):
        assert after[key] == pytest.approx(summary[key], rel=2e-5, abs=1e-5), key
    for name in ("Qha", "Qki", "Qpo", "Qh", "CV_Q_tv", "CV_Q_mv"):
        assert after["mean_flows"][name] == pytest.approx(
            summary["mean_flows"][name], rel=2e-5, abs=1e-5
        )
