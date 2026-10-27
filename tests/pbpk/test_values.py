"""Tests of the values with provenance and the units."""

from pkpdsim.pbpk import MODEL_UNITS, U, Value, as_value, fitted


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
