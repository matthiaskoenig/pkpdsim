"""The parameters of the sensitivity analyses, on the toy model."""

import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from pbpk.conftest import definition
from sbmlsim import Q
from sbmlsim.simulation import Change, Simulation, SteadyState
from sbmlsim.simulation.sampling import Uniform

from pkpdsim.pbpk.builder import build
from pkpdsim.sensitivity import (
    SEVERITIES,
    NonzeroLocalAnalysis,
    ParameterType,
    SensitivityGroup,
    SensitivityParameter,
    distributions,
    free_parameters,
    library_parameters,
    notes_sources,
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


def _parameters(*uids: str) -> list[SensitivityParameter]:
    return [SensitivityParameter(uid=uid, name=uid, value=1.0) for uid in uids]


def test_free_parameters_without_set_targets_keeps_all() -> None:
    parameters = _parameters("Ka", "CL", "f_renal_function")
    simulation = Simulation(end=10, steps=5, preinit_changes={"PODOSE": Q(1, "mg")})
    assert free_parameters(parameters, [simulation]) == parameters
    assert free_parameters(parameters, []) == parameters


def test_free_parameters_leaves_out_a_preinit_change() -> None:
    parameters = _parameters("Ka", "CL", "f_renal_function")
    simulation = Simulation(
        end=10, steps=5, preinit_changes={"f_renal_function": 0.5, "PODOSE": Q(1, "mg")}
    )
    free = free_parameters(parameters, [simulation])
    assert [p.uid for p in free] == ["Ka", "CL"]


def test_free_parameters_leaves_out_a_change() -> None:
    parameters = _parameters("Ka", "CL", "f_renal_function")
    simulation = Simulation(
        end=10, steps=5, changes=[Change([5.0], {"CL": Q(2.0, "1/hr")})]
    )
    free = free_parameters(parameters, [simulation])
    assert [p.uid for p in free] == ["Ka", "f_renal_function"]


def test_free_parameters_unites_the_simulations_and_keeps_the_order() -> None:
    parameters = _parameters("d", "c", "b", "a")
    simulations = [
        Simulation(end=10, steps=5, preinit_changes={"c": 1.0}),
        Simulation(end=10, steps=5, changes=[Change([1.0], {"a": 1.0})]),
    ]
    assert [p.uid for p in free_parameters(parameters, simulations)] == ["d", "b"]
    # a generator of simulations is consumed once
    assert [p.uid for p in free_parameters(parameters, iter(simulations))] == ["d", "b"]


def test_free_parameters_logs_the_removed_parameters(
    caplog: pytest.LogCaptureFixture,
) -> None:
    simulation = Simulation(end=10, steps=5, preinit_changes={"CL": 1.0})
    with caplog.at_level("INFO", logger="pkpdsim.sensitivity"):
        free_parameters(_parameters("Ka", "CL"), [simulation])
    assert "Removing fixed parameter: CL" in caplog.text
    assert "Ka" not in caplog.text


def test_free_parameters_leaves_out_a_target_of_the_presimulation() -> None:
    simulation = Simulation(
        end=10, steps=5, presimulation=SteadyState(preinit_changes={"CL": 1.0})
    )
    free = free_parameters(_parameters("Ka", "CL"), [simulation])
    assert [p.uid for p in free] == ["Ka"]


def test_free_parameters_leaves_out_a_change_at_several_times() -> None:
    simulation = Simulation(
        end=48, steps=5, changes=[Change([0.0, 24.0], {"CL": Q(2.0, "1/hr")})]
    )
    free = free_parameters(_parameters("Ka", "CL"), [simulation])
    assert [p.uid for p in free] == ["Ka"]


def test_free_parameters_compares_the_ids_as_written() -> None:
    simulation = Simulation(end=10, steps=5, changes=[Change([1.0], {"[S]": 1.0})])
    free = free_parameters(_parameters("S"), [simulation])
    assert [p.uid for p in free] == ["S"]


def derived_model(path: Path) -> Path:
    """Write the toy model with the derived parameter `k2 = 2 * k` into `path`."""
    import libsbml

    document = libsbml.readSBMLFromFile(str(MODEL_PATH))
    model = document.getModel()
    k2 = model.createParameter()
    k2.setId("k2")
    k2.setConstant(True)
    k2.setValue(0.0)
    assignment = model.createInitialAssignment()
    assignment.setSymbol("k2")
    assignment.setMath(libsbml.parseL3Formula("2 * k"))
    target = path / "derived.xml"
    libsbml.writeSBMLToFile(document, str(target))
    return target


def test_derived_parameters_can_be_excluded(tmp_path: Path) -> None:
    path = derived_model(tmp_path)
    assert {p.uid for p in parameters_from_sbml(path)} == {"k", "Vc", "k2"}
    kept = parameters_from_sbml(path, exclude_derived=True)
    assert {p.uid for p in kept} == {"k", "Vc"}
    both = parameters_from_sbml(path, exclude_ids={"Vc"}, exclude_derived=True)
    assert [p.uid for p in both] == ["k"]


def library_model(tmp_path: Path) -> Path:
    """Build the toy library model and get its flat SBML."""
    return build(definition(), tmp_path).flat_path


IDS = ("BW", "f_renal_function", "CP_severity", "HF_severity", "HCT_ref", "fu_ref_drg")


def test_library_parameters_refuse_a_missing_id(tmp_path: Path) -> None:
    path = library_model(tmp_path)
    with pytest.raises(ValueError, match=r"nope1.*nope2|nope2.*nope1"):
        library_parameters(
            path, ["BW", "nope1", "nope2"], fit_ids=["BW"], scaling_ids=[]
        )


def test_library_parameters_bound_severities_and_scalings(tmp_path: Path) -> None:
    path = library_model(tmp_path)
    parameters = library_parameters(
        path,
        IDS,
        fit_ids={"fu_ref_drg"},
        scaling_ids={"f_renal_function"},
        references={"f_renal_function": "KDIGO"},
        relative=0.2,
    )
    assert [p.uid for p in parameters] == list(IDS)
    by_uid = {p.uid: p for p in parameters}
    assert SEVERITIES == {"CP_severity": 3.0, "HF_severity": 4.0}
    for uid, top in SEVERITIES.items():
        p = by_uid[uid]
        assert (p.lower_bound, p.upper_bound) == (0.0, pytest.approx(0.2 * top))
        assert p.type is ParameterType.SCALING
    p = by_uid["f_renal_function"]
    assert (p.lower_bound, p.upper_bound) == (pytest.approx(0.8), pytest.approx(1.2))
    assert p.type is ParameterType.SCALING
    assert p.reference == "KDIGO"
    p = by_uid["fu_ref_drg"]
    assert p.type is ParameterType.FIT
    assert p.lower_bound == pytest.approx(p.value * 0.8)
    assert p.upper_bound == pytest.approx(p.value * 1.2)
    # a parameter with sources in its notes is data with its references
    p = by_uid["HCT_ref"]
    assert p.type is ParameterType.DATA
    assert p.reference
    assert p.lower_bound == pytest.approx(p.value * 0.8)
    assert p.upper_bound == pytest.approx(p.value * 1.2)
    assert by_uid["BW"].reference


def test_library_parameters_refuse_ids_outside_the_list(tmp_path: Path) -> None:
    path = library_model(tmp_path)
    with pytest.raises(ValueError, match="GFR_ref"):
        library_parameters(path, ["BW"], fit_ids=["GFR_ref"], scaling_ids=[])


def test_the_sources_of_the_notes(tmp_path: Path) -> None:
    import libsbml

    model = libsbml.readSBMLFromFile(str(library_model(tmp_path))).getModel()
    assert notes_sources(model.getParameter("BW"))
    assert notes_sources(model.getParameter("CP_severity")) == []


def test_the_nonzero_local_analysis_varies_nonzero_values_only(tmp_path: Path) -> None:
    parameters = library_parameters(
        library_model(tmp_path), IDS, fit_ids=[], scaling_ids=[]
    )
    varied = {
        p.uid
        for p in NonzeroLocalAnalysis(
            groups=[SensitivityGroup("g", "G", Simulation(end=1.0))]
        ).varied(parameters)
    }
    assert "CP_severity" not in varied
    assert "BW" in varied


def test_the_pk_observables_of_a_substance() -> None:
    """The analysis of a concentration and its outputs, as the copies of the models had them."""
    from pkpdutils import NCAOptions
    from sbmlsim.simulation import PK, Formula

    from pkpdsim.sensitivity import pk_observables

    outputs = {"auc_inf_obs": "mM*min", "cmax": "mM", "cl_f": "ml/min"}
    options = NCAOptions(lloq=1e-6)
    observables = pk_observables(
        "gli",
        "Carm_gli",
        outputs,
        dose="PODOSE_gli",
        route="oral",
        options=options,
        per_dose_mass={"cl_f", "vz_f"},
        mr="Mr_gli",
    )
    assert observables == [
        PK(
            "pk_gli",
            "Carm_gli",
            dose="PODOSE_gli",
            route="oral",
            options=options,
            parameters=list(outputs),
        ),
        Formula("gli_auc_inf_obs", "pk_gli.auc_inf_obs", unit="mM*min"),
        Formula("gli_cmax", "pk_gli.cmax", unit="mM"),
        Formula("gli_cl_f", "pk_gli.cl_f / at(Mr_gli, 0)", unit="ml/min"),
    ]
    # a metabolite without a dose, its outputs named by the concentration
    metabolite = pk_observables("Carm_m1", "Carm_m1", {"cmax": "mM"}, pk_id="pk_m1")
    assert metabolite == [
        PK("pk_m1", "Carm_m1", parameters=["cmax"]),
        Formula("Carm_m1_cmax", "pk_m1.cmax", unit="mM"),
    ]


def test_an_output_per_dose_mass_needs_the_molecular_weight() -> None:
    from pkpdsim.sensitivity import pk_observables

    with pytest.raises(ValueError, match="molecular weight"):
        pk_observables("gli", "Carm_gli", {"cl_f": "ml/min"}, per_dose_mass={"cl_f"})


def test_the_maximum_after_the_dose() -> None:
    """The maximum of the timecourse from the time of the dose, 0, on."""
    from pkpdsim.sensitivity import max_after_dose

    time = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
    values = {"FLUID__diuresis": np.array([9.0, 5.0, 1.0, 3.0, 2.0])}
    assert max_after_dose(time, values) == 3.0
    with pytest.raises(ValueError, match="one selection"):
        max_after_dose(time, {"a": time, "b": time})
