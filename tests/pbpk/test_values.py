"""Tests of the values with provenance and the units."""

import libsbml
from sbmlutils.factory import Uncertainty, UncertParameter

from pkpdsim.pbpk import MODEL_UNITS, Formula, U, Value, as_value, fitted
from pkpdsim.pbpk.values import check_domain, equal_by_value
from pkpdsim.sbml.domains import BOOLEAN, FRACTION


def test_fitted_takes_fit_record_kwargs() -> None:
    """A fit record gives a value with the source fit and its uncertainties."""
    v = fitted({"value": 2.5, "uncertainties": [], "unit": U.dimensionless})
    assert v.value == 2.5
    assert v.source == "fit"
    assert v.uncertainties == []
    assert v.unit == U.dimensionless
    assert v.notes is None


def test_fitted_keeps_notes() -> None:
    """A fitted value takes a remark, e.g. the literature of the process."""
    v = fitted(
        {"value": 0.6, "uncertainties": [], "unit": U.per_min},
        notes="The uptake by OAT1.",
    )
    assert v.notes == "The uptake by OAT1."
    assert v.source == "fit"
    assert v.value == 0.6


def test_as_value_wraps_float() -> None:
    """A number gets the unit, a value is kept."""
    assert as_value(3.0, U.mM) == Value(3.0, U.mM)
    v = Value(1.0, U.kg, source="ICRP89")
    assert as_value(v, U.mM) is v


def test_model_units() -> None:
    """The model units are those of the models."""
    assert MODEL_UNITS.time == U.min
    assert MODEL_UNITS.substance == U.mmole
    assert MODEL_UNITS.volume == U.litre
    assert MODEL_UNITS.area == U.m2


def _uncertainty(math: str) -> Uncertainty:
    """An uncertainty with a distribution, a new object at every call."""
    return Uncertainty(
        uncertParameters=[
            UncertParameter(type=libsbml.DISTRIB_UNCERTTYPE_DISTRIBUTION, math=math)
        ]
    )


def test_values_with_equal_uncertainties_are_equal_by_value() -> None:
    """sbmlutils compares uncertainties by identity, `equal_by_value` by their fields."""
    a = Value(
        0.5, U.dimensionless, "fit", uncertainties=[_uncertainty("normal(0.5, 0.1)")]
    )
    b = Value(
        0.5, U.dimensionless, "fit", uncertainties=[_uncertainty("normal(0.5, 0.1)")]
    )
    assert a != b
    assert equal_by_value(a, b)
    c = Value(
        0.5, U.dimensionless, "fit", uncertainties=[_uncertainty("normal(0.5, 0.2)")]
    )
    assert not equal_by_value(a, c)
    assert not equal_by_value(a, Value(0.5, U.dimensionless, "fit"))
    # containers compare their elements, a list is no tuple
    assert equal_by_value({"x": [a]}, {"x": [b]})
    assert not equal_by_value([a], (b,))
    assert equal_by_value(1, 1.0)


def test_check_domain_requires_a_domain() -> None:
    """A parameter which states no domain is a problem."""
    assert check_domain("parameter 'k'", "value", 0.5, None) == [
        "parameter 'k': value states no domain"
    ]


def test_check_domain_reports_a_value_outside_it() -> None:
    """A number outside the domain is named with the domain."""
    (problem,) = check_domain("parameter 'f'", "value", 1.3, FRACTION)
    assert "1.3" in problem and "[0, 1]" in problem


def test_check_domain_passes_a_value_inside_it() -> None:
    """A number inside the domain is no problem."""
    assert check_domain("parameter 'f'", "value", 0.94, FRACTION) == []


def test_check_domain_passes_a_formula() -> None:
    """A formula has no number to check."""
    formula = Formula("1 - fabs_drg", U.dimensionless)
    assert check_domain("parameter 'f'", "value", formula, FRACTION) == []


def test_check_domain_reports_a_non_integral_value_of_a_boolean() -> None:
    """A flag of 0.5 is neither off nor on."""
    (problem,) = check_domain("parameter 'ckd'", "value", 0.5, BOOLEAN)
    assert (
        "0.5" in problem and "not admitted" in problem and "[0, 1] boolean" in problem
    )
