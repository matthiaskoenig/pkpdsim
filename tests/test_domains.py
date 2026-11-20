"""Tests of the domain of a quantity, written and read as a distrib range span."""

import logging
import math
from pathlib import Path

import libsbml
import pytest
from sbmlutils.factory import (
    Compartment,
    Model,
    Parameter,
    Species,
    Uncertainty,
    UncertParameter,
    UncertSpan,
    create_model,
)
from sbmlutils.io.antimony import sbml_to_antimony

from pkpdsim.sbml.domains import (
    BOOLEAN,
    DOMAIN,
    FRACTION,
    NONNEGATIVE,
    REAL,
    Domain,
    Values,
    domain_of,
    domain_uncertainty,
    domains,
    integers,
)


def _write(tmp_path: Path, objects: list) -> Path:
    path = tmp_path / "m.xml"
    create_model(model=Model(sid="m", objects=objects), filepath=path)
    return path


def _read(path: Path) -> libsbml.Model:
    return libsbml.readSBMLFromFile(str(path)).getModel()


def test_named_domains() -> None:
    assert (FRACTION.lower, FRACTION.upper) == (0.0, 1.0)
    assert (NONNEGATIVE.lower, NONNEGATIVE.upper) == (0.0, math.inf)
    assert (REAL.lower, REAL.upper) == (-math.inf, math.inf)


def test_clip_moves_only_the_bound_outside() -> None:
    assert FRACTION.clip(0.799, 1.127) == (0.799, 1.0)
    assert FRACTION.clip(-0.2, 0.5) == (0.0, 0.5)
    assert FRACTION.clip(0.2, 0.5) == (0.2, 0.5)
    assert NONNEGATIVE.clip(0.85, 1.15) == (0.85, 1.15)


def test_clip_keeps_a_bound_which_is_the_bound() -> None:
    assert FRACTION.clip(1.0, 1.0) == (1.0, 1.0)


def test_contains_rejects_what_is_not_finite() -> None:
    """An infinite bound is open and a state which blows up is no value."""
    for d in (NONNEGATIVE, REAL):
        for x in (math.inf, -math.inf, math.nan):
            assert not d.contains(x)
    assert REAL.contains(-1e300)


def test_contains_accepts_the_bounds() -> None:
    assert FRACTION.contains(0.0) and FRACTION.contains(1.0)
    assert not FRACTION.contains(1.0000001)
    assert NONNEGATIVE.contains(0.0) and not NONNEGATIVE.contains(-1e-9)


def test_text() -> None:
    assert str(FRACTION) == "[0, 1]"
    assert str(NONNEGATIVE) == "[0, inf)"
    assert str(REAL) == "(-inf, inf)"


def test_round_trip_of_a_bounded_domain(tmp_path: Path) -> None:
    path = _write(
        tmp_path,
        [
            Parameter(
                sid="f",
                value=0.5,
                unit="dimensionless",
                uncertainties=[domain_uncertainty(FRACTION)],
            )
        ],
    )
    assert domains(path) == {"f": FRACTION}


def test_round_trip_of_an_infinite_domain(tmp_path: Path) -> None:
    path = _write(
        tmp_path,
        [
            Parameter(
                sid="v",
                value=1.0,
                unit="dimensionless",
                uncertainties=[domain_uncertainty(NONNEGATIVE)],
            ),
            Parameter(
                sid="x",
                value=1.0,
                unit="dimensionless",
                uncertainties=[domain_uncertainty(REAL)],
            ),
        ],
    )
    assert domains(path) == {"v": NONNEGATIVE, "x": REAL}


def test_a_span_without_an_upper_value_reads_as_infinite(tmp_path: Path) -> None:
    namespaces = libsbml.SBMLNamespaces(3, 1, "distrib", 1)
    doc = libsbml.SBMLDocument(namespaces)
    doc.setPackageRequired("distrib", True)
    model = doc.createModel()
    parameter = model.createParameter()
    parameter.setId("v")
    uncertainty = parameter.getPlugin("distrib").createUncertainty()
    uncertainty.setName(DOMAIN)
    span = uncertainty.createUncertSpan()
    span.setType(libsbml.DISTRIB_UNCERTTYPE_RANGE)
    span.setValueLower(0.0)
    assert not span.isSetValueUpper()
    path = tmp_path / "span.xml"
    libsbml.writeSBMLToFile(doc, str(path))
    assert domain_of(_read(path).getParameter("v")) == NONNEGATIVE


def test_the_domain_is_read_next_to_the_uncertainty_of_a_source(tmp_path: Path) -> None:
    source = Uncertainty(
        uncertParameters=[
            UncertParameter(
                type=libsbml.DISTRIB_UNCERTTYPE_STANDARDDEVIATION, value=0.1
            ),
            UncertSpan(
                type=libsbml.DISTRIB_UNCERTTYPE_CONFIDENCEINTERVAL,
                valueLower=0.3,
                valueUpper=0.7,
            ),
        ],
        name="source",
    )
    path = _write(
        tmp_path,
        [
            Parameter(
                sid="f",
                value=0.5,
                unit="dimensionless",
                uncertainties=[source, domain_uncertainty(FRACTION)],
            )
        ],
    )
    assert domain_of(_read(path).getParameter("f")) == FRACTION


def test_an_element_without_a_domain_reads_as_none(tmp_path: Path) -> None:
    source = Uncertainty(
        uncertParameters=[
            UncertSpan(
                type=libsbml.DISTRIB_UNCERTTYPE_CONFIDENCEINTERVAL,
                valueLower=0.3,
                valueUpper=0.7,
            )
        ],
        name="source",
    )
    path = _write(
        tmp_path,
        [
            Parameter(sid="p", value=0.5, unit="dimensionless", uncertainties=[source]),
            Parameter(sid="q", value=0.5, unit="dimensionless"),
        ],
    )
    model = _read(path)
    assert domain_of(model.getParameter("p")) is None
    assert domain_of(model.getParameter("q")) is None
    assert domains(path) == {}


def test_a_document_without_the_distrib_plugin_reads_as_empty(tmp_path: Path) -> None:
    doc = libsbml.SBMLDocument(3, 2)
    model = doc.createModel()
    parameter = model.createParameter()
    parameter.setId("p")
    path = tmp_path / "plain.xml"
    libsbml.writeSBMLToFile(doc, str(path))
    assert (
        domain_of(libsbml.readSBMLFromFile(str(path)).getModel().getParameter("p"))
        is None
    )
    assert domains(path) == {}


def test_compartments_and_species_carry_a_domain(tmp_path: Path) -> None:
    path = _write(
        tmp_path,
        [
            Compartment(
                sid="C",
                value=1.0,
                unit="litre",
                uncertainties=[domain_uncertainty(NONNEGATIVE)],
            ),
            Species(
                sid="S",
                compartment="C",
                substanceUnit="mole",
                initialConcentration=1.0,
                uncertainties=[domain_uncertainty(NONNEGATIVE)],
            ),
        ],
    )
    assert domains(path) == {"C": NONNEGATIVE, "S": NONNEGATIVE}


def test_an_uncertainty_named_domain_without_a_range_span_reads_as_none(
    tmp_path: Path,
) -> None:
    interval = Uncertainty(
        uncertParameters=[
            UncertSpan(
                type=libsbml.DISTRIB_UNCERTTYPE_CONFIDENCEINTERVAL,
                valueLower=0.3,
                valueUpper=0.7,
            )
        ],
        name=DOMAIN,
    )
    spread = Uncertainty(
        uncertParameters=[
            UncertParameter(
                type=libsbml.DISTRIB_UNCERTTYPE_STANDARDDEVIATION, value=0.1
            )
        ],
        name=DOMAIN,
    )
    path = _write(
        tmp_path,
        [
            Parameter(
                sid="a", value=0.5, unit="dimensionless", uncertainties=[interval]
            ),
            Parameter(sid="b", value=0.5, unit="dimensionless", uncertainties=[spread]),
        ],
    )
    model = _read(path)
    assert domain_of(model.getParameter("a")) is None
    assert domain_of(model.getParameter("b")) is None
    assert domains(path) == {}


def test_text_of_a_non_finite_bound_does_not_raise() -> None:
    assert str(Domain(math.inf, math.inf)) == "(inf, inf)"
    assert str(Domain(-math.inf, -math.inf)) == "(-inf, -inf)"
    assert str(Domain(math.nan, 1.0)) == "(nan, 1]"


def test_writing_a_domain_logs_no_sbmlutils_warning(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """A domain is metadata of a value: it asks for no SBO term (no hint)."""
    p = Parameter(
        sid="f",
        value=0.5,
        unit="dimensionless",
        name="fraction",
        sboTerm="SBO:0000002",
        uncertainties=[domain_uncertainty(FRACTION)],
    )
    with caplog.at_level(logging.WARNING, logger="sbmlutils"):
        create_model(
            model=Model(sid="m", name="m", objects=[p]),
            filepath=tmp_path / "m.xml",
        )
    assert [r.getMessage() for r in caplog.records] == []


def test_a_boolean_domain_is_the_integers_of_zero_and_one() -> None:
    assert BOOLEAN.contains(0) and BOOLEAN.contains(1)
    assert not BOOLEAN.contains(0.5)
    assert (BOOLEAN.lower, BOOLEAN.upper) == (0, 1)


def test_an_integer_domain_clips_onto_its_integers_across_zero() -> None:
    d = integers(-2, 2)
    for value, expected in [(-2.6, -2.0), (-0.4, 0.0), (0.5, 0.0), (7.0, 2.0)]:
        clipped = d.clip(value)
        assert clipped == expected
        assert d.contains(clipped)


def test_an_integer_domain_without_an_integer_is_rejected() -> None:
    with pytest.raises(ValueError, match=r"\[0\.2, 0\.8\]"):
        integers(0.2, 0.8)
    with pytest.raises(ValueError, match="no integer"):
        Domain(0.2, 0.8, Values.BOOLEAN)


def test_clipping_an_infinite_value_returns_the_infinite_bound() -> None:
    assert integers(0, math.inf).clip(math.inf) == math.inf
    assert integers(-math.inf, 0).clip(-math.inf) == -math.inf


def test_clipping_not_a_number_raises() -> None:
    with pytest.raises(ValueError, match="not a number"):
        integers(0, 3).clip(math.nan)


def test_an_integer_domain_does_not_contain_a_non_integral_value() -> None:
    assert not integers(0, 3).contains(2.5)
    assert integers(0, 3).contains(2)


def test_a_continuous_domain_contains_a_non_integral_value() -> None:
    assert FRACTION.contains(0.94)


def test_clip_of_an_integer_domain_returns_a_value_it_contains() -> None:
    d = integers(0, 3)
    for x in (-1.0, 0.4, 2.5, 7.0):
        assert d.contains(d.clip(x))


def test_the_text_of_an_integer_domain_states_its_value_set() -> None:
    assert str(integers(0, 3)) == "[0, 3] integer"
    assert str(BOOLEAN) == "[0, 1] boolean"
    assert str(FRACTION) == "[0, 1]"


def _one_parameter(tmp_path: Path, d: Domain) -> Path:
    return _write(
        tmp_path,
        [
            Parameter(
                sid="x",
                value=1,
                unit="dimensionless",
                uncertainties=[domain_uncertainty(d)],
            )
        ],
    )


def test_an_integer_domain_round_trips_through_sbml(tmp_path: Path) -> None:
    path = _one_parameter(tmp_path, integers(0, 3))
    assert domains(path) == {"x": integers(0, 3)}
    assert domains(path)["x"].values is Values.INTEGER


def test_a_boolean_domain_round_trips_through_sbml(tmp_path: Path) -> None:
    path = _one_parameter(tmp_path, BOOLEAN)
    assert domains(path) == {"x": BOOLEAN}


def test_a_continuous_domain_writes_no_value_set(tmp_path: Path) -> None:
    path = _one_parameter(tmp_path, FRACTION)
    model = _read(path)  # the model owns the elements, keep it alive
    uncertainty = model.getParameter("x").getPlugin("distrib").getUncertainty(0)
    assert uncertainty.getNumUncertParameters() == 1
    assert isinstance(uncertainty.getUncertParameter(0), libsbml.UncertSpan)


@pytest.mark.parametrize("d", [integers(0, 3), BOOLEAN])
def test_a_document_with_an_integer_domain_validates(tmp_path: Path, d: Domain) -> None:
    path = _one_parameter(tmp_path, d)
    doc = libsbml.readSBMLFromFile(str(path))
    doc.checkConsistency()
    assert doc.getNumErrors() == 0


@pytest.mark.parametrize("d", [integers(0, 3), BOOLEAN])
def test_the_value_one_of_the_value_set_child_is_needed_by_antimony(
    tmp_path: Path, d: Domain
) -> None:
    """Keep the `value` of the value set child: libAntimony fails without it.

    libsbml accepts the child without a value, libAntimony does not and blames
    the sibling `uncertSpan` (sys-bio/antimony#190). Remove the value only when
    that issue is fixed.
    """
    # the archive of a model is written through antimony, which refuses an
    # externalParameter without a value
    path = _one_parameter(tmp_path, d)
    assert sbml_to_antimony(path)
