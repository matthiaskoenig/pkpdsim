"""Verify the SBML heart against an independent explicit circuit integration.

The reference assembles a directed hydraulic network with an incidence matrix,
eliminates the storage-free hepatic pressure numerically, and integrates the
remaining volumes with SciPy DOP853. It does not evaluate SBML formulas or use
the production cardiovascular builder to calculate a derivative. This checks
the translation of standard resistance/compliance/elastance components into
the coupled SBML DAE; it is not external clinical validation of the defaults.
"""

import numpy as np
import pytest
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from pkpdsim.pbpk import DynamicHeart, PBPKDefinition, Substance, Variants, build, ids
from pkpdsim.pbpk.builder import BuildResult
from pkpdsim.pbpk.organs import PORTAL_ORGANS, SYSTEMIC_ORGANS, Organ

from .cardiovascular_solver import simulate

STORAGE = ("ar", "ve", "po", "pa", "pv", "lv", "rv")
VOLUME_IDS = ("Var", "Vve", "Vpo", "CV_V_pa", "CV_V_pv", "CV_V_lv", "CV_V_rv")
VASCULAR = STORAGE[:5]


@pytest.fixture(scope="module")
def reference_heart(tmp_path_factory: pytest.TempPathFactory) -> BuildResult:
    return build(
        PBPKDefinition(
            sid="independent_heart_reference",
            substances=[Substance("trc", "transport probe", 100.0)],
            variants=Variants(circulation=DynamicHeart()),
        ),
        tmp_path_factory.mktemp("independent_heart_reference"),
    )


def _circuit(
    values: dict[str, float],
) -> tuple[
    NDArray[np.float64], NDArray[np.int64], NDArray[np.int64], NDArray[np.float64]
]:
    """Assemble a network independently of the production SBML rate rules."""
    nodes = (*STORAGE, "li")
    edges: list[tuple[str, str, float]] = []
    for organs, destination in ((SYSTEMIC_ORGANS, "ve"), (PORTAL_ORGANS, "po")):
        for organ in sorted(organs):
            conductance = values[ids.impairment_factor(ids.organ_flow(organ))] / (
                values[f"CV_R_{organ.value}"]
                * values[f"CV_fR_{organ.value}"]
                * values["CV_fR_systemic"]
            )
            edges.append(("ar", destination, conductance))
    hepatic_artery = values[ids.impairment_factor(ids.organ_flow(Organ.LIVER))] / (
        values["CV_R_li"] * values["CV_fR_li"] * values["CV_fR_systemic"]
    )
    edges += [
        ("ar", "li", hepatic_artery),
        ("po", "li", 1 / (values["CV_R_portal"] * values["CV_fR_portal"])),
        ("li", "ve", 1 / (values["CV_R_hepatic"] * values["CV_fR_hepatic"])),
        ("pa", "pv", 1 / (values["CV_R_lu"] * values["CV_fR_lu"])),
    ]
    valve_conductance = values["CV_scale"] / values["CV_R_valve"]
    edges += [
        (source, destination, valve_conductance)
        for source, destination in (
            ("ve", "rv"),
            ("rv", "pa"),
            ("pv", "lv"),
            ("lv", "ar"),
        )
    ]
    source = np.array([nodes.index(edge[0]) for edge in edges], dtype=np.int64)
    destination = np.array([nodes.index(edge[1]) for edge in edges], dtype=np.int64)
    conductance = np.array([edge[2] for edge in edges])
    incidence = np.zeros((len(nodes), len(edges)))
    for column, (start, end, _) in enumerate(edges):
        incidence[nodes.index(start), column] = -1.0
        incidence[nodes.index(end), column] = 1.0
    return incidence, source, destination, conductance


def _pressures(
    time: float, volume: NDArray[np.float64], values: dict[str, float]
) -> NDArray[np.float64]:
    pressure = np.empty(8)
    scale = values["CV_scale"]
    for index, node in enumerate(VASCULAR):
        pressure[index] = values[f"CV_P_{node}_ref"] + (
            volume[index] - values[f"CV_V_{node}_ref"]
        ) / (values[f"CV_C_{node}"] * scale)
    phase = values["CV_HR"] * time % 1.0
    systole = values["CV_systole"]
    activation = np.sin(np.pi * phase / systole) ** 2 if phase < systole else 0.0
    for index, node in ((5, "lv"), (6, "rv")):
        minimum, maximum = values[f"CV_E_{node}_min"], values[f"CV_E_{node}_max"]
        elastance = (
            minimum + activation * values["CV_contractility"] * (maximum - minimum)
        ) / scale
        pressure[index] = elastance * (volume[index] - values[f"CV_V0_{node}"] * scale)
    return pressure


@pytest.mark.parametrize("hepatic_resistance", [1.0, 3.0])
def test_dynamic_heart_matches_independent_explicit_circuit(
    reference_heart: BuildResult, hepatic_resistance: float
) -> None:
    parameters = [
        "CV_scale",
        "CV_HR",
        "CV_systole",
        "CV_contractility",
        "CV_R_valve",
        "CV_R_portal",
        "CV_R_hepatic",
        "CV_R_lu",
        "CV_fR_systemic",
        "CV_fR_portal",
        "CV_fR_hepatic",
        "CV_fR_lu",
        *(f"CV_R_{organ.value}" for organ in Organ if organ is not Organ.LUNG),
        *(f"CV_fR_{organ.value}" for organ in Organ if organ is not Organ.LUNG),
        *(
            ids.impairment_factor(ids.organ_flow(organ))
            for organ in Organ
            if organ is not Organ.LUNG
        ),
        *(f"CV_{kind}_{node}_ref" for kind in ("P", "V") for node in VASCULAR),
        *(f"CV_C_{node}" for node in VASCULAR),
        *(f"CV_E_{node}_{bound}" for node in ("lv", "rv") for bound in ("min", "max")),
        *(f"CV_V0_{node}" for node in ("lv", "rv")),
    ]
    pressure_ids = [f"CV_P_{node}" for node in (*STORAGE, "li")]
    result = simulate(
        reference_heart.flat_path,
        [*parameters, *VOLUME_IDS, *pressure_ids],
        end=3.0 / 75.0,
        steps=600,
        maximum_time_step=1.0 / (75 * 100),
        preinit_changes={"CV_fR_hepatic": hepatic_resistance},
    )
    values = {sid: float(result[sid][0]) for sid in parameters}
    incidence, source, destination, conductance = _circuit(values)

    def pressures(time: float, volume: NDArray[np.float64]) -> NDArray[np.float64]:
        pressure = _pressures(time, volume, values)

        def hepatic_balance(candidate: float) -> float:
            pressure[7] = candidate
            flow = np.maximum(pressure[source] - pressure[destination], 0) * conductance
            return float(incidence[7] @ flow)

        boundaries = pressure[[0, 1, 2]]
        pressure[7] = brentq(
            hepatic_balance,
            float(boundaries.min()),
            float(boundaries.max()),
            xtol=1e-12,
        )
        return pressure

    def rate(time: float, volume: NDArray[np.float64]) -> NDArray[np.float64]:
        pressure = pressures(time, volume)
        flow = np.maximum(pressure[source] - pressure[destination], 0) * conductance
        return incidence[:7] @ flow

    times = result["time"]
    solution = solve_ivp(
        rate,
        (float(times[0]), float(times[-1])),
        [float(result[sid][0]) for sid in VOLUME_IDS],
        method="DOP853",
        t_eval=times,
        max_step=1.0 / (values["CV_HR"] * 100),
        rtol=1e-10,
        atol=1e-12,
    )
    assert solution.success, solution.message
    expected_pressure = np.stack(
        [
            pressures(float(time), volume)
            for time, volume in zip(times, solution.y.T, strict=True)
        ]
    )
    for index, sid in enumerate(VOLUME_IDS):
        np.testing.assert_allclose(
            result[sid], solution.y[index], rtol=2e-5, atol=2e-7, err_msg=sid
        )
    for index, sid in enumerate(pressure_ids):
        np.testing.assert_allclose(
            result[sid], expected_pressure[:, index], rtol=2e-4, atol=2e-3, err_msg=sid
        )
