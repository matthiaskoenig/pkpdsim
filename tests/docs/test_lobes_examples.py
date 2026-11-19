"""The liver lobe examples of the documentation are genuine, current IDA runs.

The page shows precomputed courses before a reader moves a control and runs
the same compiled model in the browser afterwards, so two things have to hold:
the committed dataset is what the generator produces now, and the committed
browser model was compiled from the SBML the library builds now. Both are
compared against a fresh run here, which is what keeps a stale number or a
stale `model.wasm` out of the published documentation.
"""

import itertools
import json
import runpy
import shutil
import subprocess
from importlib.metadata import version
from pathlib import Path
from typing import Any

import numpy as np
import pytest

ROOT = Path(__file__).parents[2]
SCRIPT = runpy.run_path(str(ROOT / "scripts/lobes_examples.py"))
OUTPUT: Path = SCRIPT["OUTPUT"]
SOLVERS: Path = SCRIPT["SOLVERS"]
RUNTIME: str = SCRIPT["RUNTIME"]
SPECIES = SCRIPT["SPECIES"]
EXAMPLES = ROOT / "docs/pbpk/examples"


@pytest.fixture(scope="module")
def committed() -> dict[str, Any]:
    """The dataset of the page as it is committed."""
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def fresh() -> dict[str, Any]:
    """The dataset of a fresh run of the generator, without the browser models."""
    return SCRIPT["generate"](browser=None)


def test_the_committed_dataset_is_a_current_run(
    committed: dict[str, Any], fresh: dict[str, Any]
) -> None:
    """Every number of the page comes from the model as the library builds it now."""
    assert committed["schema_version"] == fresh["schema_version"]
    assert set(committed["species"]) == set(fresh["species"])
    for species, of_species in fresh["species"].items():
        was = committed["species"][species]
        assert was["sbml_sha256"] == of_species["sbml_sha256"], species
        assert was["columns"] == of_species["columns"], species
        assert was["lobes"] == of_species["lobes"], species
        assert set(was["cases"]) == set(of_species["cases"]), species
        for name, case in of_species["cases"].items():
            assert was["cases"][name]["changes"] == case["changes"], (species, name)
            for sid, series in case["series"].items():
                np.testing.assert_allclose(
                    was["cases"][name]["series"][sid],
                    series,
                    rtol=1e-6,
                    atol=1e-12,
                    err_msg=f"{species} {name} {sid}",
                )


def test_the_provenance_names_the_installed_releases(committed: dict[str, Any]) -> None:
    provenance = committed["provenance"]
    assert provenance["generator"] == "scripts/lobes_examples.py"
    assert provenance["sbmlsim_version"] == version("sbmlsim")
    assert provenance["sbmlsim_core_version"] == version("sbmlsim-core")
    assert provenance["time_unit"] == "min"
    assert "not a reproduction" in provenance["notice"]


def test_every_species_has_a_compiled_browser_model(committed: dict[str, Any]) -> None:
    """The committed `model.wasm` was compiled from the SBML of the dataset."""
    for species in SPECIES:
        config = json.loads(
            (SOLVERS / species.value / "config.json").read_text(encoding="utf-8")
        )
        assert (
            config["sbml_sha256"] == committed["species"][species.value]["sbml_sha256"]
        )
        assert config["sbmlsim_core_version"] == version("sbmlsim-core")
        assert config["columns"] == committed["species"][species.value]["columns"]
        assert config["runtime"] == f"../{RUNTIME}"
        assert (SOLVERS / species.value / "model.wasm").stat().st_size > 0
        assert (SOLVERS / species.value / "model.json").stat().st_size > 0


def test_the_lobe_models_share_the_runtime_of_the_cardiovascular_export() -> None:
    """The engine is committed once, so the two exports must be the same release."""
    shared = json.loads(
        (EXAMPLES / RUNTIME / "config.json").read_text(encoding="utf-8")
    )
    assert shared["sbmlsim_core_version"] == version("sbmlsim-core")
    for name in ("engine.wasm", "mathlib.wasm", "sbmlsim-core.js"):
        assert (EXAMPLES / RUNTIME / name).stat().st_size > 0, name
    for species in SPECIES:
        assert not (SOLVERS / species.value / "engine.wasm").exists()


def test_a_resection_regenerates_the_liver(committed: dict[str, Any]) -> None:
    """The remnant of a 70 % hepatectomy of the rat grows back over the five days."""
    case = committed["species"]["rat"]["cases"]["hepatectomy70"]
    remnant = case["series"]["f_remnant"]
    assert remnant[0] == pytest.approx(0.30, abs=0.01)
    assert remnant[-1] > 0.9
    assert all(b >= a - 1e-9 for a, b in itertools.pairwise(remnant))


def test_a_deportalised_lobe_atrophies_and_its_neighbours_grow(
    committed: dict[str, Any],
) -> None:
    """A portal vein ligation moves mass from the ligated lobes to the others."""
    case = committed["species"]["rat"]["cases"]["pvl70"]
    for token in ("lll", "lml", "rml"):
        volume = case["series"][f"Vli_{token}"]
        assert volume[-1] < 0.5 * volume[0], token
    for token in ("cl", "rl"):
        volume = case["series"][f"Vli_{token}"]
        assert volume[-1] > 1.5 * volume[0], token


def test_an_intact_liver_does_not_move(committed: dict[str, Any]) -> None:
    """The remodelling is at rest when every lobe is at the reference perfusion."""
    for of_species in committed["species"].values():
        series = of_species["cases"]["intact"]["series"]
        for sid, values in series.items():
            if sid.startswith("Vli_") or sid == "f_remnant":
                np.testing.assert_allclose(values, values[0], rtol=1e-3, err_msg=sid)


def test_every_case_is_finite_and_the_lobe_flows_sum(
    committed: dict[str, Any],
) -> None:
    """No case of the page diverges and the lobes are the liver in every one."""
    for species, of_species in committed["species"].items():
        tokens = [lobe["token"] for lobe in of_species["lobes"]]
        for name, case in of_species["cases"].items():
            for sid, values in case["series"].items():
                assert np.isfinite(values).all(), (species, name, sid)
            for flow in ("Qha", "Qpo", "Qh"):
                total = np.sum([case["series"][f"{flow}_{t}"] for t in tokens], axis=0)
                np.testing.assert_allclose(
                    total, case["series"][flow], rtol=1e-6, err_msg=f"{name} {flow}"
                )


def test_the_javascript_of_the_page(tmp_path: Path) -> None:
    """The model helpers and the app of the page, under `node:test`."""
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not installed")
    for name in ("lobes_model.test.mjs", "lobes_app.test.mjs"):
        test = Path(__file__).parent / name
        result = subprocess.run(
            [node, str(test)],
            cwd=test.parent,
            capture_output=True,
            text=True,
            timeout=60,
            check=False,
        )
        assert result.returncode == 0, result.stdout + result.stderr
