"""The parameters of the sensitivity analyses, on the toy model."""

import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from sbmlsim.simulation.sampling import Uniform

from pkpdsim.sensitivity import (
    ParameterType,
    SensitivityParameter,
    distributions,
    parameters_from_sbml,
    parameters_to_df,
    sample_size,
    set_bounds,
    set_relative_bounds,
    write_parameters,
)
from pkpdsim_toy import MODEL_PATH

COLUMNS = [
    "uid",
    "name",
    "value",
    "lower_bound",
    "upper_bound",
    "unit",
    "type",
    "reference",
]


def test_the_constant_parameters_and_compartments_of_the_model() -> None:
    parameters = {p.uid: p for p in parameters_from_sbml(MODEL_PATH)}
    assert set(parameters) == {"k", "Vc"}
    assert parameters["k"].value == pytest.approx(0.2)
    assert parameters["k"].name == "rate of the decay"
    assert parameters["Vc"].value == pytest.approx(1.0)
    assert math.isnan(parameters["k"].lower_bound)
    assert parameters["k"].type is ParameterType.NA


def test_parameters_are_excluded() -> None:
    assert [p.uid for p in parameters_from_sbml(MODEL_PATH, exclude_ids={"Vc"})] == [
        "k"
    ]


def test_the_bounds_of_the_analysis() -> None:
    parameters = parameters_from_sbml(MODEL_PATH)
    set_bounds(
        parameters,
        [("k", 0.1, 0.3, ParameterType.FIT), ("unknown", 0, 1, ParameterType.DATA)],
    )
    set_relative_bounds(parameters, 0.15)
    by_uid = {p.uid: p for p in parameters}
    assert (by_uid["k"].lower_bound, by_uid["k"].upper_bound) == (0.1, 0.3)
    assert by_uid["k"].type is ParameterType.FIT
    assert by_uid["Vc"].lower_bound == pytest.approx(0.85)
    assert by_uid["Vc"].upper_bound == pytest.approx(1.15)


def test_the_table_is_sorted_by_type_and_id() -> None:
    parameters = [
        SensitivityParameter(uid="b", name="b", type=ParameterType.SCALING),
        SensitivityParameter(uid="a", name="a", type=ParameterType.SCALING),
        SensitivityParameter(uid="c", name="c", type=ParameterType.FIT),
    ]
    df = parameters_to_df(parameters)
    assert list(df.columns) == COLUMNS
    assert df["uid"].tolist() == ["c", "a", "b"]
    assert df["type"].tolist() == ["fitted", "scaling", "scaling"]
    assert parameters_to_df(parameters, sort=False)["uid"].tolist() == ["b", "a", "c"]


def test_the_table_is_written_as_tsv_and_latex(tmp_path: Path) -> None:
    parameters = parameters_from_sbml(MODEL_PATH)
    set_relative_bounds(parameters, 0.15)
    df = write_parameters(parameters, tmp_path / "sensitivity")
    read = pd.read_csv(tmp_path / "sensitivity" / "parameters.tsv", sep="\t")
    assert list(read.columns) == COLUMNS
    assert np.allclose(read["value"], df["value"])
    tex = (tmp_path / "sensitivity" / "parameters.tex").read_text(encoding="utf-8")
    assert r"lower\_bound" in tex
    assert "_bound" not in tex.replace(r"\_bound", "")


def test_the_distributions_of_a_global_analysis() -> None:
    parameters = parameters_from_sbml(MODEL_PATH)
    set_relative_bounds(parameters, 0.1)
    result = distributions(parameters)
    assert set(result) == {"k", "Vc"}
    assert isinstance(result["k"], Uniform)
    assert result["k"].lower == pytest.approx(0.18)


def test_a_parameter_without_bounds_has_no_distribution() -> None:
    with pytest.raises(ValueError, match="k"):
        distributions(parameters_from_sbml(MODEL_PATH))


def test_a_negative_value_has_the_bounds_of_sbmlsim_08() -> None:
    """The relative bounds of a negative value are reversed, as in sbmlsim 0.8."""
    parameters = [SensitivityParameter(uid="dE", name="dE", value=-2.0)]
    set_relative_bounds(parameters, 0.15)
    assert parameters[0].lower_bound == pytest.approx(-1.7)
    assert parameters[0].upper_bound == pytest.approx(-2.3)


def test_reversed_bounds_have_no_distribution() -> None:
    parameters = [
        SensitivityParameter(
            uid="k", name="k", value=1.0, lower_bound=0.9, upper_bound=1.1
        ),
        SensitivityParameter(uid="dE", name="dE", value=-2.0),
    ]
    set_relative_bounds(parameters, 0.15)
    with pytest.raises(ValueError, match=r"'dE' \(-1.7 > -2.3\)"):
        distributions(parameters)


def test_the_sample_size_of_a_global_analysis(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("PKPDSIM_SENSITIVITY_N", raising=False)
    assert sample_size(4096) == 4096
    monkeypatch.setenv("PKPDSIM_SENSITIVITY_N", "16")
    assert sample_size(4096) == 16


@pytest.mark.parametrize("value", ["0", "-1", "many"])
def test_a_smoke_run_needs_a_positive_sample_size(
    monkeypatch: pytest.MonkeyPatch, value: str
) -> None:
    monkeypatch.setenv("PKPDSIM_SENSITIVITY_N", value)
    with pytest.raises(ValueError, match="PKPDSIM_SENSITIVITY_N"):
        sample_size(4096)
