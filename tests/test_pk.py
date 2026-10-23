"""The non-compartmental analysis of the scans of a model."""

import numpy as np
import pytest
from pkpdutils import Dose, NCAOptions, Route, Timecourse, nca_single
from sbmlsim import Q
from sbmlsim.model import AbstractModel
from sbmlsim.result import ScanResult
from sbmlsim.simulation import Change, Dimension, Scan, Simulation
from sbmlsim.simulator import Simulator

from pkpdsim.pk import nca_frame, scan_nca
from pkpdsim.scan import point_timecourse
from pkpdsim_toy import MODEL_PATH

TOY = AbstractModel(
    source=MODEL_PATH,
    language_type=AbstractModel.LanguageType.SBML,
    selections=["time", "[A]"],
)
RATES = Dimension("dim_k", values={"k": Q([0.1, 0.2, 0.4], "1/min")})
PARAMETERS = ["auc_last", "auc_inf_obs", "cmax", "tmax", "thalf", "lambda_z", "cl"]


def run(simulation: Simulation, dimensions: list[Dimension]) -> ScanResult:
    return Simulator(n_workers=1).run(TOY, Scan(simulation, dimensions))


def per_point(res: ScanResult, dose: float, lloq_fraction: float | None) -> dict:
    """The analysis of every point on its own, `nca_single`."""
    rows: dict[str, list[float]] = {name: [] for name in PARAMETERS}
    for k in range(res.ds.sizes["dim_k"]):
        time, values = point_timecourse(res, "[A]", dim_k=k)
        cmax = float(values.max())
        result = nca_single(
            Timecourse(
                time=time,
                value=values,
                time_unit=res.units["time"],
                unit=res.units["[A]"],
                lloq=None if lloq_fraction is None else lloq_fraction * cmax,
                dose=Dose(amount=dose, unit="mmol", route=Route.IV_BOLUS),
                substance="A",
            ),
            options=NCAOptions(),
        )
        for name in PARAMETERS:
            rows[name].append(float(result.ds[name].values))
    return rows


@pytest.mark.parametrize("lloq_fraction", [None, 1e-3])
@pytest.mark.parametrize(
    "simulation",
    [Simulation(end=60, steps=120), Simulation(end=60)],
    ids=["grid", "steps of the integrator"],
)
def test_the_analysis_of_a_scan_is_the_one_of_every_point(
    simulation: Simulation, lloq_fraction: float | None
) -> None:
    res = run(simulation, [RATES])
    result = scan_nca(
        res,
        "[A]",
        dose=Q(np.full(3, 2.0), "mmole"),
        route=Route.IV_BOLUS,
        substance="A",
        lloq_fraction=lloq_fraction,
    )
    df = result.to_dataframe()
    assert df["dim_k"].tolist() == [0, 1, 2]
    expected = per_point(res, 2.0, lloq_fraction)
    for name in PARAMETERS:
        np.testing.assert_allclose(df[name], expected[name], rtol=1e-12)
    # the faster decay has the shorter half-life
    assert df["thalf"].is_monotonic_decreasing


def test_a_change_gives_the_state_after_it() -> None:
    """The state before a change at its time is left out, as per point."""
    redose = Simulation(end=60, changes=[Change([30.0], {"[A]": Q(1.0, "mM")})])
    res = run(redose, [RATES])
    result = scan_nca(
        res, "[A]", dose=Q(1.0, "mmole"), route=Route.IV_BOLUS, substance="A"
    )
    expected = per_point(res, 1.0, None)
    np.testing.assert_allclose(result.to_dataframe()["cmax"], expected["cmax"])


def test_the_frame_has_the_rows_of_every_point() -> None:
    """The quantities and the flags of `nca_single` of every point, as rows."""
    res = run(Simulation(end=60), [RATES])
    df = nca_frame(
        scan_nca(res, "[A]", dose=Q(2.0, "mmole"), route=Route.IV_BOLUS, substance="A")
    )
    assert df["dim_k"].tolist() == [0, 1, 2]
    time, values = point_timecourse(res, "[A]", dim_k=1)
    single = nca_single(
        Timecourse(
            time=time,
            value=values,
            time_unit=res.units["time"],
            unit=res.units["[A]"],
            dose=Dose(amount=2.0, unit="mmol", route=Route.IV_BOLUS),
            substance="A",
        )
    )
    for name, quantity in single.to_quantities().items():
        np.testing.assert_allclose(df[name].iloc[1], quantity.magnitude, rtol=1e-12)
        assert df[f"{name}_unit"].iloc[1] == str(quantity.units)
    assert df["nca_flags"].iloc[1] == ", ".join(single.flags())


def test_a_substance_without_a_dose() -> None:
    """A metabolite has the parameters of its curve, none of a dose."""
    res = run(Simulation(end=60), [RATES])
    df = nca_frame(
        scan_nca(res, "[A]", dose=None, route=None, substance="A", lloq_fraction=1e-3)
    )
    assert df["dim_k"].tolist() == [0, 1, 2]
    assert df["auc_last"].notna().all()
    assert "cl" not in df.columns
