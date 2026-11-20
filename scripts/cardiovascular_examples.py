"""Generate the reference data of the browser cardiovascular examples.

Run ``uv run python scripts/cardiovascular_examples.py`` after changing the
cardiovascular model or its defaults. The browser displays these IDA results;
it does not integrate a second, JavaScript approximation of the beating heart.
The pressure-network calibration is read from the same SBML model so the live
resistance explorer starts from the PBPK individual's actual reference flows.
"""

import argparse
import hashlib
import json
import shutil
from dataclasses import dataclass
from importlib.metadata import distribution, version
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

import numpy as np
import sbmlsim_core
from numpy.typing import NDArray
from sbmlsim.model import AbstractModel
from sbmlsim.model.model_sbml import IntegratorSettings
from sbmlsim.result.timecourse import TimecourseResult
from sbmlsim.simulation import Change, Simulation
from sbmlsim.simulator import Simulator
from sbmlsim_core import web as core_web

from pkpdsim.pbpk import (
    DynamicHeart,
    PBPKDefinition,
    PressureFlow,
    Sex,
    Variants,
    build,
    ids,
)
from pkpdsim.pbpk.organs import PORTAL_ORGANS, Organ

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs/pbpk/examples/cardiovascular-data.json"
END = 0.8
INTERVENTION = 0.4
INITIAL_HEART_RATE = 75.0
BEATS = 3
SUMMARY_STEP_SECONDS = 0.0002
INTEGRATOR: IntegratorSettings = {
    "absolute_tolerance": 1e-10,
    "relative_tolerance": 1e-8,
    "maximum_time_step": 1 / (100 * 100),
}
COLUMNS = (
    "time",
    "CV_P_ar",
    "CV_P_ve",
    "CV_P_po",
    "CV_P_li",
    "CV_P_pa",
    "CV_P_pv",
    "CV_P_lv",
    "CV_P_rv",
    "CV_V_lv",
    "CV_V_rv",
    "CO",
    "Qlu",
    ids.QPO,
    ids.QH,
    "CV_V_blood_total",
    "CV_Q_tv",
    "CV_Q_pv",
    "CV_Q_mv",
    ids.QHA,
    *(ids.organ_flow(organ) for organ in Organ if organ is not Organ.LUNG),
)


@dataclass(frozen=True)
class Example:
    """A controlled intervention after a shared drug-free settling interval."""

    sid: str
    label: str
    description: str
    parameter: str | None = None
    value: float = 1.0

    @property
    def changes(self) -> dict[str, float]:
        """The model parameters changed at 24 seconds."""
        return {} if self.parameter is None else {self.parameter: self.value}

    @property
    def heart_rate(self) -> float:
        """The rate after the intervention, in cycles per minute."""
        return self.value if self.parameter == "CV_HR" else INITIAL_HEART_RATE


EXAMPLES = (
    Example("baseline", "Reference circulation", "Unchanged reference circulation."),
    Example(
        "hepatic-2",
        "Hepatic resistance 2x",
        "Double hepatic outflow resistance.",
        "CV_fR_hepatic",
        2.0,
    ),
    Example(
        "hepatic-4",
        "Hepatic resistance 4x",
        "Increase hepatic outflow resistance fourfold.",
        "CV_fR_hepatic",
        4.0,
    ),
    Example(
        "hepatic-8",
        "Hepatic resistance 8x",
        "Increase hepatic outflow resistance eightfold.",
        "CV_fR_hepatic",
        8.0,
    ),
    Example(
        "systemic-low",
        "Systemic resistance 0.7x",
        "Lower systemic organ resistance by 30 percent.",
        "CV_fR_systemic",
        0.7,
    ),
    Example(
        "systemic-high",
        "Systemic resistance 1.4x",
        "Raise systemic organ resistance by 40 percent.",
        "CV_fR_systemic",
        1.4,
    ),
    Example(
        "contractility-low",
        "Contractility 0.65x",
        "Reduce active ventricular elastance by 35 percent.",
        "CV_contractility",
        0.65,
    ),
    Example(
        "heart-rate-high",
        "Heart rate 100/min",
        "Increase heart rate from 75 to 100 cycles per minute.",
        "CV_HR",
        100.0,
    ),
)


def _simulator() -> Simulator:
    return Simulator(engine="core", n_workers=1, **INTEGRATOR)


def pressure_reference(path: Path) -> dict[str, Any]:
    """Read the resistance calibration and baseline flows from consistent SBML."""
    organs = [organ for organ in Organ if organ is not Organ.LUNG]
    selections = [
        "CV_P_ar",
        "CV_P_ve",
        "CV_P_po",
        "CV_P_li",
        "CO",
        "CV_R_portal",
        "CV_R_hepatic",
        *(f"CV_R_{organ.value}" for organ in organs),
        *(ids.organ_flow(organ) for organ in organs),
    ]
    result = _simulator().simulate(
        AbstractModel(source=path, selections=selections),
        Simulation(end=1e-9, times=[0.0]),
    )
    value = {key: float(result[key][0]) for key in selections}
    return {
        "arterial_pressure": value["CV_P_ar"],
        "venous_pressure": value["CV_P_ve"],
        "portal_pressure": value["CV_P_po"],
        "hepatic_pressure": value["CV_P_li"],
        "cardiac_output": value["CO"],
        "portal_resistance": value["CV_R_portal"],
        "hepatic_resistance": value["CV_R_hepatic"],
        "organs": [
            {
                "id": organ.value,
                "name": organ.name.replace("_", " ").capitalize(),
                "portal": organ in PORTAL_ORGANS,
                "resistance": value[f"CV_R_{organ.value}"],
                "flow": value[ids.organ_flow(organ)],
            }
            for organ in organs
        ],
    }


def _window_summary(
    result: TimecourseResult, start: float, end: float, heart_rate: float
) -> dict[str, Any]:
    time = result.time
    selected = (time >= start - 1e-12) & (time <= end + 1e-12)
    indices = np.flatnonzero(selected)
    # The event has one row on each side. The pre-intervention interval ends
    # at the first row of its boundary, before any parameter is changed.
    if end == INTERVENTION:
        boundary = np.flatnonzero(np.isclose(time, end, rtol=0, atol=1e-12))
        indices = indices[indices <= boundary[0]]
    t = time[indices]

    def mean(key: str) -> float:
        return float(np.trapezoid(result[key][indices], t) / (end - start))

    return {
        "arterial_systolic": float(result["CV_P_ar"][indices].max()),
        "arterial_diastolic": float(result["CV_P_ar"][indices].min()),
        "mean_arterial_pressure": mean("CV_P_ar"),
        "mean_portal_pressure": mean("CV_P_po"),
        "mean_hepatic_pressure": mean("CV_P_li"),
        "mean_cardiac_output": mean("CO"),
        "mean_lung_flow": mean("Qlu"),
        "stroke_volume": 1000 * mean("CO") / heart_rate,
        "heart_rate": heart_rate,
        "mean_hepatic_flow": mean(ids.QH),
        "mean_portal_flow": mean(ids.QPO),
        "mean_flows": {
            key: mean(key)
            for key in COLUMNS
            if key == "CO" or key.startswith(("Q", "CV_Q_"))
        },
    }


def _table(result: TimecourseResult, indices: NDArray[np.int64]) -> dict[str, Any]:
    data = np.column_stack([result[column][indices] for column in COLUMNS])
    data[:, 0] *= 60
    for column in ("CV_V_lv", "CV_V_rv"):
        data[:, COLUMNS.index(column)] *= 1000
    return {"columns": list(COLUMNS), "values": np.round(data, 5).tolist()}


def heart_example(path: Path, example: Example) -> dict[str, Any]:
    """Solve a shared baseline and a controlled change with SUNDIALS IDA."""
    beat_start = END - BEATS / example.heart_rate
    beat_times = np.linspace(beat_start, END, 601)
    before_start = INTERVENTION - BEATS / INITIAL_HEART_RATE
    before_samples = (
        int(np.ceil((INTERVENTION - before_start) * 60 / SUMMARY_STEP_SECONDS)) + 1
    )
    after_samples = int(np.ceil((END - beat_start) * 60 / SUMMARY_STEP_SECONDS)) + 1
    times = np.unique(
        np.concatenate(
            (
                np.linspace(0, END, 15001),
                np.linspace(before_start, INTERVENTION, before_samples),
                np.linspace(beat_start, END, after_samples),
                beat_times,
                np.array([before_start, INTERVENTION, beat_start, END]),
            )
        )
    )
    result = _simulator().simulate(
        AbstractModel(source=path, selections=list(COLUMNS)),
        Simulation(
            end=END,
            times=times.tolist(),
            changes=[Change(INTERVENTION, example.changes)] if example.changes else [],
        ),
    )
    time = result.time
    overview = np.searchsorted(time, np.linspace(0, END, 501))
    event = np.flatnonzero(np.isclose(time, INTERVENTION, rtol=0, atol=1e-12))
    overview = np.unique(np.concatenate((overview, event)))
    beats = np.searchsorted(time, beat_times)
    return {
        "id": example.sid,
        "label": example.label,
        "description": example.description,
        "changes": example.changes,
        "heart_rate": example.heart_rate,
        "series": _table(result, overview),
        "beats": _table(result, beats),
        "summary": {
            "before": _window_summary(
                result, before_start, INTERVENTION, INITIAL_HEART_RATE
            ),
            "after": _window_summary(result, beat_start, END, example.heart_rate),
            "maximum_blood_volume_error": float(
                np.max(
                    np.abs(result["CV_V_blood_total"] - result["CV_V_blood_total"][0])
                )
            ),
        },
    }


def export_browser_model(path: Path, output: Path, *, runtime: bool = True) -> None:
    """Export the installed core host, runtime and compiled heart for Web Workers.

    All runtime bytes and their notices come from the same installed release
    that generated the reference trajectories, not a separately pinned copy.
    """
    output.mkdir(parents=True, exist_ok=True)
    # Canonical LF text keeps exports and fingerprints independent of the host OS.
    sbml = path.read_text(encoding="utf-8")
    compiled = sbmlsim_core.compile_sbml(sbml)
    (output / "model.wasm").write_bytes(compiled.wasm)
    (output / "model.json").write_text(
        json.dumps(compiled.description, separators=(",", ":"), allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    (output / "model.xml").write_text(sbml, encoding="utf-8", newline="\n")
    if runtime:
        for name, wasm in sbmlsim_core.runtime_wasm().items():
            (output / f"{name}.wasm").write_bytes(wasm)
        (output / "sbmlsim-core.js").write_text(
            core_web.script().replace("\r\n", "\n"), encoding="utf-8", newline="\n"
        )
        (output / "THIRD_PARTY.md").write_text(
            sbmlsim_core.THIRD_PARTY.replace("\r\n", "\n"),
            encoding="utf-8",
            newline="\n",
        )
        package = distribution("sbmlsim-core")
        license_file = next(
            file
            for file in package.files or ()
            if str(file).endswith("licenses/LICENSE")
        )
        shutil.copyfile(str(package.locate_file(license_file)), output / "LICENSE")
    config = {
        "generated_by": "scripts/cardiovascular_examples.py",
        "engine": "SUNDIALS IDA",
        "sbmlsim_core_version": version("sbmlsim-core"),
        "sbml_sha256": hashlib.sha256(sbml.encode("utf-8")).hexdigest(),
        "columns": list(COLUMNS),
        "start": 0.0,
        "end": END,
        "intervention": INTERVENTION,
        "initial_heart_rate": INITIAL_HEART_RATE,
        "summary_beats": BEATS,
        "summary_step_seconds": SUMMARY_STEP_SECONDS,
        "settings": {
            **INTEGRATOR,
            "absolute_tolerance": {"scalar": INTEGRATOR["absolute_tolerance"]},
            "maximum_num_steps": 1000000,
            "initial_time_step": 0,
            "stiff": True,
            "linear_solver": "auto",
            "jacobian": "auto",
        },
    }
    (output / "config.json").write_text(
        json.dumps(config, separators=(",", ":"), allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def generate(
    *, sex: Sex = Sex.MALE, browser_output: Path | None = None
) -> dict[str, Any]:
    """Build both SBML variants and generate the documentation dataset."""
    with TemporaryDirectory(prefix="pkpdsim-cardiovascular-examples-") as directory:
        folder = Path(directory)
        pressure = build(
            PBPKDefinition(
                "pressure_examples",
                [],
                sex=sex,
                variants=Variants(circulation=PressureFlow()),
            ),
            folder / "pressure",
        )
        heart = build(
            PBPKDefinition(
                "heart_examples",
                [],
                sex=sex,
                variants=Variants(circulation=DynamicHeart()),
            ),
            folder / "heart",
        )
        if browser_output is not None:
            export_browser_model(
                heart.flat_path, browser_output, runtime=sex is Sex.MALE
            )
        fingerprint = hashlib.sha256(
            pressure.flat_path.read_text(encoding="utf-8").encode("utf-8")
            + heart.flat_path.read_text(encoding="utf-8").encode("utf-8")
        ).hexdigest()
        return {
            "schema_version": 1,
            "provenance": {
                "generator": "scripts/cardiovascular_examples.py",
                "engine": "core",
                "solver": "SUNDIALS IDA",
                "model": f"illustrative adult human {sex.value}",
                "sex": sex.value,
                "sbmlsim_version": version("sbmlsim"),
                "sbmlsim_core_version": version("sbmlsim-core"),
                "sbml_sha256": fingerprint,
                "intervention_time": INTERVENTION * 60,
                "settling_time": INTERVENTION * 60,
                "end_time": END * 60,
                "time_unit": "s",
                "summary_beats": BEATS,
                "summary_time_step_seconds": SUMMARY_STEP_SECONDS,
                "integrator": INTEGRATOR,
                "units": {
                    key: "s"
                    if key == "time"
                    else "ml"
                    if key in ("CV_V_lv", "CV_V_rv")
                    else "l"
                    if key == "CV_V_blood_total"
                    else "mmHg"
                    if key.startswith("CV_P_")
                    else "l/min"
                    for key in COLUMNS
                },
                "notice": "Illustrative model parameters, not a calibrated clinical model. The browser displays IDA trajectories without interpolating between scenarios.",
            },
            "pressure": pressure_reference(pressure.flat_path),
            "cases": [heart_example(heart.flat_path, example) for example in EXAMPLES],
        }


def export_anatomy(output: Path) -> None:
    """Reuse the physiology figures' scaled bodies and organ anchors for the app."""
    from pkpdsim.pbpk.visual.body import reference_body

    bodies = {}
    for sex in Sex:
        body = reference_body(sex.value)
        anchors = {}
        for organ in Organ:
            key = organ.value
            if key in ("mu", "bo", "sk", "ad"):
                lm, cx = body.landmarks, body.anatomy.midline
                fraction = {"mu": 0.55, "bo": 0.7, "sk": 0.25, "ad": 0.45}[key]
                point = body.at(
                    cx - 8, lm["armpit"] + fraction * (lm["bottom"] - lm["armpit"])
                )
            elif key == "go" and sex is Sex.FEMALE:
                point = body.centroid("uterus")
            else:
                point = body.centroid(key)
            anchors[key] = [round(float(v), 4) for v in point]
        bodies[sex.value] = {
            "svg": body.svg(arteries=False, tags=True),
            "height": body.bounds()[2],
            "anchors": anchors,
        }
    output.write_text(
        json.dumps(
            {
                "generated_by": "scripts/cardiovascular_examples.py",
                "attribution": "Anatomy adapted from EMBL-EBI Expression Atlas anatomograms (CC BY 4.0), using pkpdsim physiology figures.",
                "source": "https://github.com/ebi-gene-expression-group/anatomogram",
                "license": "https://creativecommons.org/licenses/by/4.0/",
                "bodies": bodies,
            },
            separators=(",", ":"),
        )
        + "\n",
        encoding="utf-8",
        newline="\n",
    )


def write(output: Path = OUTPUT) -> None:
    """Write male/female examples and physiology anatomy with portable text."""
    output.parent.mkdir(parents=True, exist_ok=True)
    for sex in Sex:
        target = output if sex is Sex.MALE else output.parent / "female" / output.name
        target.parent.mkdir(parents=True, exist_ok=True)
        data = generate(sex=sex, browser_output=target.parent / "solver")
        target.write_text(
            json.dumps(data, separators=(",", ":"), ensure_ascii=False, allow_nan=False)
            + "\n",
            encoding="utf-8",
            newline="\n",
        )
    export_anatomy(output.parent / "cardiovascular-anatomy.json")


def main() -> None:
    """Run the reference generation from a checkout of pkpdsim."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    options = parser.parse_args()
    write(options.output)


if __name__ == "__main__":
    main()
