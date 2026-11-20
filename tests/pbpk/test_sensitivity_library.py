"""The parameters of the sensitivity analyses of a library model.

`library_parameters` takes the parameters of a model built by the library
with their bounds from the notes and the list of the model; the toy library
model is the definition of `tests/pbpk/conftest.py`.
"""

import logging
import math
from collections.abc import Collection
from pathlib import Path

import libsbml
import pytest
from sbmlsim.sensitivity.parameters import ParameterType, SensitivityParameter
from sbmlsim.sensitivity.runner import SensitivityGroup
from sbmlsim.simulation import Simulation

from pkpdsim.pbpk import ModelParameter, Organ, U, Variants
from pkpdsim.pbpk.builder import build
from pkpdsim.pbpk.variants.intestine import SegmentedIntestine
from pkpdsim.sbml.domains import Domain, domains
from pkpdsim.sensitivity import (
    SEVERITY_IDS,
    NonzeroLocalAnalysis,
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
    tops = {"CP_severity": 3.0, "HF_severity": 4.0}
    for uid, top in tops.items():
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


FABS = 0.94
RELATIVE = 0.15
LOWER = FABS * (1 - RELATIVE)  # 0.799
UPPER = FABS * (1 + RELATIVE)  # 1.081


def fraction_model(tmp_path: Path, value: float = FABS) -> Path:
    """Build the toy model with an oral route and a fraction absorbed of `value`."""
    definition_ = definition(
        variants=Variants(
            organs={Organ.SMALL_INTESTINE: SegmentedIntestine(fabs=value)}
        )
    )
    return build(definition_, tmp_path).flat_path


def fraction(
    path: Path, *, clip: bool = True, unclipped: Collection[str] = ()
) -> SensitivityParameter:
    """Get the parameter `fabs_drg` of the analysis of the model."""
    parameters = library_parameters(
        path,
        ["fabs_drg"],
        fit_ids=[],
        scaling_ids=[],
        relative=RELATIVE,
        clip=clip,
        unclipped=unclipped,
    )
    return parameters[0]


def test_severities_are_ids_in_a_fixed_order() -> None:
    assert SEVERITY_IDS == ("CP_severity", "HF_severity")


def test_a_severity_is_varied_to_the_top_of_its_domain(tmp_path: Path) -> None:
    # relative=0.15, the domain of CP_severity is [0, 3]
    p = library_parameters(
        library_model(tmp_path),
        ["CP_severity", "HF_severity"],
        fit_ids=[],
        scaling_ids=[],
        relative=RELATIVE,
    )
    assert (p[0].lower_bound, p[0].upper_bound) == (0.0, pytest.approx(0.45))
    assert (p[1].lower_bound, p[1].upper_bound) == (0.0, pytest.approx(0.6))


def test_a_severity_without_a_finite_top_raises(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        "pkpdsim.sensitivity.domains",
        lambda path: {"CP_severity": Domain(0.0, math.inf)},
    )
    with pytest.raises(ValueError, match="CP_severity"):
        library_parameters(
            library_model(tmp_path), ["CP_severity"], fit_ids=[], scaling_ids=[]
        )


def test_a_fraction_is_clipped_to_its_domain(tmp_path: Path) -> None:
    p = fraction(fraction_model(tmp_path))
    assert p.value == pytest.approx(FABS)
    assert (p.lower_bound, p.upper_bound) == (pytest.approx(LOWER), 1.0)


def test_the_clip_is_logged(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    with caplog.at_level(logging.INFO, logger="pkpdsim.sensitivity"):
        fraction(fraction_model(tmp_path))
    # the upper bound before the clip is 0.94 * 1.15 = 1.081
    assert "fabs_drg" in caplog.text
    assert "1.081" in caplog.text
    assert "[0, 1]" in caplog.text


def test_unclipped_keeps_the_span(tmp_path: Path) -> None:
    p = fraction(fraction_model(tmp_path), unclipped={"fabs_drg"})
    assert p.upper_bound == pytest.approx(UPPER)


def test_clip_false_keeps_every_span(tmp_path: Path) -> None:
    p = fraction(fraction_model(tmp_path), clip=False)
    assert (p.lower_bound, p.upper_bound) == pytest.approx((LOWER, UPPER))


def test_a_parameter_inside_its_domain_is_untouched(tmp_path: Path) -> None:
    p = library_parameters(
        library_model(tmp_path), ["BW"], fit_ids=[], scaling_ids=[], relative=0.15
    )[0]
    assert (p.lower_bound, p.upper_bound) == pytest.approx(
        (0.85 * p.value, 1.15 * p.value)
    )


def _without_domains(path: Path, tmp_path: Path) -> Path:
    """Write a copy of the model whose elements state no domain."""
    document = libsbml.readSBMLFromFile(str(path))
    model = document.getModel()
    for elements in (
        model.getListOfParameters(),
        model.getListOfCompartments(),
        model.getListOfSpecies(),
    ):
        for element in elements:
            plugin = element.getPlugin("distrib")
            while plugin.getNumUncertainties():
                plugin.removeUncertainty(0)
    stripped = tmp_path / "stripped_all.xml"
    libsbml.writeSBMLToFile(document, str(stripped))
    return stripped


def test_a_model_without_domains_changes_no_bound_and_says_so(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    stripped = _without_domains(fraction_model(tmp_path), tmp_path)
    with caplog.at_level(logging.INFO, logger="pkpdsim.sensitivity"):
        p = fraction(stripped)
    assert "states no domain for any quantity" in caplog.text
    assert (p.lower_bound, p.upper_bound) == pytest.approx((LOWER, UPPER))


def test_a_model_without_domains_varies_a_severity_by_the_relative_rule(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    path = _without_domains(library_model(tmp_path), tmp_path)
    assert not domains(path)
    with caplog.at_level(logging.WARNING, logger="pkpdsim.sensitivity"):
        p = library_parameters(
            path, ["CP_severity", "BW"], fit_ids=[], scaling_ids=[], relative=RELATIVE
        )
    assert "built before the domain contract" in caplog.text
    assert "no severity is varied over its domain" in caplog.text
    assert len([r for r in caplog.records if r.levelno >= logging.WARNING]) == 1
    severity, weight = p
    assert (severity.lower_bound, severity.upper_bound) == (0.0, 0.0)
    assert weight.lower_bound == pytest.approx(0.85 * weight.value)


def test_a_severity_without_a_domain_in_a_model_with_domains_raises(
    tmp_path: Path,
) -> None:
    document = libsbml.readSBMLFromFile(str(library_model(tmp_path)))
    plugin = document.getModel().getParameter("CP_severity").getPlugin("distrib")
    while plugin.getNumUncertainties():
        plugin.removeUncertainty(0)
    partial = tmp_path / "partial.xml"
    libsbml.writeSBMLToFile(document, str(partial))
    assert domains(partial)
    with pytest.raises(ValueError, match="CP_severity") as error:
        library_parameters(partial, ["CP_severity"], fit_ids=[], scaling_ids=[])
    assert "has to be rebuilt" in str(error.value)


def test_an_unclipped_id_outside_the_parameters_raises(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="fabs_typo"):
        library_parameters(
            fraction_model(tmp_path),
            ["fabs_drg"],
            fit_ids=[],
            scaling_ids=[],
            unclipped=["fabs_typo"],
        )


def test_a_value_outside_its_domain_raises(tmp_path: Path) -> None:
    document = libsbml.readSBMLFromFile(str(fraction_model(tmp_path)))
    document.getModel().getParameter("fabs_drg").setValue(1.2)
    broken = tmp_path / "broken.xml"
    libsbml.writeSBMLToFile(document, str(broken))
    with pytest.raises(ValueError, match="fabs_drg"):
        fraction(broken)


def test_an_inverted_pair_of_a_negative_value_is_left_alone(tmp_path: Path) -> None:
    """The reversed bounds of a negative value survive the clipping unchanged.

    `set_relative_bounds` of sbmlsim reverses the bounds of a negative value on
    purpose and `distributions` refuses such a parameter by name, so that a
    model sets explicit bounds with `set_bounds`; the clipping must not turn the
    pair into one which looks valid. The value -0.9 in the domain [-1, 1] gives
    the pair (-0.765, -1.035). The outcome follows from `Domain.clip` limiting
    each bound on its own side, so this test fails if that ever changes.
    """
    definition_ = definition(
        parameters=[ModelParameter("neg", -0.9, U.dimensionless, domain=Domain(-1, 1))]
    )
    path = build(definition_, tmp_path).flat_path
    p = library_parameters(
        path, ["neg"], fit_ids=[], scaling_ids=[], relative=RELATIVE
    )[0]
    assert p.value == pytest.approx(-0.9)
    assert (p.lower_bound, p.upper_bound) == pytest.approx((-0.765, -1.035))


def test_values_which_are_not_a_number_are_counted(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    document = libsbml.readSBMLFromFile(str(fraction_model(tmp_path)))
    document.getModel().getParameter("fabs_drg").setValue(math.nan)
    broken = tmp_path / "nan.xml"
    libsbml.writeSBMLToFile(document, str(broken))
    with caplog.at_level(logging.INFO, logger="pkpdsim.sensitivity"):
        p = fraction(broken)
    assert "1 of the 1 parameters" in caplog.text
    assert "not a number" in caplog.text
    assert "fabs_drg" in caplog.text
    assert math.isnan(p.lower_bound)


def test_a_nan_value_with_finite_bounds_is_not_clipped(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """The warning says "not clipped", so a NaN value keeps its finite bounds."""
    path = fraction_model(tmp_path)
    p = fraction(path)
    p.value = math.nan
    p.lower_bound, p.upper_bound = -0.5, 2.0
    from pkpdsim.sensitivity import _clip_to_domains

    with caplog.at_level(logging.WARNING, logger="pkpdsim.sensitivity"):
        _clip_to_domains([p], domains(path), (), path.name)
    assert (p.lower_bound, p.upper_bound) == (-0.5, 2.0)
    assert "not clipped" in caplog.text
