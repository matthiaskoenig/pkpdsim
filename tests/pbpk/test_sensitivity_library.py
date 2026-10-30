"""The parameters of the sensitivity analyses of a library model.

`library_parameters` takes the parameters of a model built by the library
with their bounds from the notes and the list of the model; the toy library
model is the definition of `tests/pbpk/conftest.py`.
"""

from pathlib import Path

import pytest
from sbmlsim.simulation import Simulation

from pkpdsim.pbpk.builder import build
from pkpdsim.sensitivity import (
    SEVERITIES,
    NonzeroLocalAnalysis,
    ParameterType,
    SensitivityGroup,
    library_parameters,
    notes_sources,
)

from .conftest import definition


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
