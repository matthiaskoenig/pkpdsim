"""Tests of the values with provenance and the units."""

import libsbml
from sbmlutils.factory import Uncertainty, UncertParameter

from pkpdsim.pbpk import MODEL_UNITS, U, Value, as_value, fitted
from pkpdsim.pbpk.values import equal_by_value


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
