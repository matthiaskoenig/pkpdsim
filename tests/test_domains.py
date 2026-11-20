"""Tests of the domain of a quantity, written and read as a distrib range span."""

import math
from pathlib import Path

import libsbml
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

from pkpdsim.sbml.domains import (
    DOMAIN,
    FRACTION,
    NONNEGATIVE,
    REAL,
    domain_of,
    domain_uncertainty,
    domains,
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
