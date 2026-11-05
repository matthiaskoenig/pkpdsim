"""The classes of impaired organ function and the labels of the PK parameters.

The values are those of the models of pkdb_models, written with the same
expressions, so that a model which takes a preset simulates exactly as before.
"""

import itertools

import pytest

from pkpdsim import impairment, pk


def test_the_renal_function_by_the_midpoints_of_the_classes() -> None:
    assert list(impairment.RENAL_FUNCTION.items()) == [
        ("Normal renal function", 101.0 / 101.0),
        ("Mild renal impairment", 69.5 / 101.0),
        ("Moderate renal impairment", 32.5 / 101.0),
        ("Severe renal impairment", 19.5 / 101.0),
    ]


def test_the_renal_function_by_the_limits_of_the_classes() -> None:
    assert list(impairment.RENAL_FUNCTION_GFR.items()) == [
        ("Normal renal function", 101.0 / 101.0),
        ("Mild renal impairment", 50.0 / 101.0),
        ("Moderate renal impairment", 35.0 / 101.0),
        ("Severe renal impairment", 20.0 / 101.0),
        ("End stage renal disease", 10.5 / 101.0),
    ]


def test_the_renal_function_on_a_base_of_90_ml_per_min() -> None:
    assert list(impairment.RENAL_FUNCTION_90.items()) == [
        ("Normal renal function", 90.0 / 90.0),
        ("Mild renal impairment", (60 + 89) / 2 / 90.0),
        ("Moderate renal impairment", (30 + 59) / 2 / 90.0),
        ("Severe renal impairment", (15 + 29) / 2 / 90.0),
        ("End stage renal disease", 9.0 / 90.0),
    ]


def test_the_renal_colors() -> None:
    assert list(impairment.RENAL_COLORS.items()) == [
        ("Normal renal function", "black"),
        ("Mild renal impairment", "#66c2a4"),
        ("Moderate renal impairment", "#2ca25f"),
        ("Severe renal impairment", "#006d2c"),
    ]
    assert list(impairment.RENAL_COLORS_END_STAGE.items()) == [
        *impairment.RENAL_COLORS.items(),
        ("End stage renal disease", "#006d5e"),
    ]


def test_the_cirrhosis_by_the_child_pugh_classes() -> None:
    assert list(impairment.CIRRHOSIS.items()) == [
        ("Control", 0),
        ("Mild cirrhosis", 0.3994897959183674),
        ("Moderate cirrhosis", 0.6979591836734694),
        ("Severe cirrhosis", 0.8127551020408164),
    ]
    assert list(impairment.CIRRHOSIS_COLORS.items()) == [
        ("Control", "black"),
        ("Mild cirrhosis", "#74a9cf"),
        ("Moderate cirrhosis", "#2b8cbe"),
        ("Severe cirrhosis", "#045a8d"),
    ]


def test_the_cardiac_function() -> None:
    assert list(impairment.CARDIAC.items()) == [
        ("Normal cardiac function", 1.0),
        ("Mild cardiac impairment", 4.75 / 5.25),
        ("Moderate cardiac impairment", 3.5 / 5.25),
        ("Severe cardiac impairment", 3 / 5.25),
        ("Cardiac failure", 0.5),
    ]
    assert list(impairment.CARDIAC_COLORS.items()) == [
        ("Normal cardiac function", "black"),
        ("Mild cardiac impairment", "#FDDBC7"),
        ("Moderate cardiac impairment", "#D6604D"),
        ("Severe cardiac impairment", "#67001F"),
        ("Cardiac failure", "#67001F"),
        ("CHF", "#67001F"),
    ]


def test_every_class_has_a_color() -> None:
    assert set(impairment.RENAL_FUNCTION) <= set(impairment.RENAL_COLORS)
    for renal in (impairment.RENAL_FUNCTION_GFR, impairment.RENAL_FUNCTION_90):
        assert set(renal) <= set(impairment.RENAL_COLORS_END_STAGE)
    assert set(impairment.CIRRHOSIS) <= set(impairment.CIRRHOSIS_COLORS)
    assert set(impairment.CARDIAC) <= set(impairment.CARDIAC_COLORS)
    assert set(impairment.CHILD_PUGH) <= set(impairment.CIRRHOSIS_COLORS)
    assert set(impairment.HEART_FAILURE) <= set(impairment.HEART_FAILURE_COLORS)


def test_the_child_pugh_severity_of_the_library() -> None:
    """The Child-Pugh classes are the severities `CP_severity` 0 to 3 of the library."""
    assert list(impairment.CHILD_PUGH.items()) == [
        ("Control", 0),
        ("Mild cirrhosis", 1),
        ("Moderate cirrhosis", 2),
        ("Severe cirrhosis", 3),
    ]
    assert list(impairment.CHILD_PUGH) == list(impairment.CIRRHOSIS)


def test_the_heart_failure_severity_of_the_library() -> None:
    """The NYHA classes are the severities `HF_severity` 0 to 4 of the library."""
    assert list(impairment.HEART_FAILURE.items()) == [
        ("Normal cardiac function", 0),
        ("NYHA I", 1),
        ("NYHA II", 2),
        ("NYHA III", 3),
        ("NYHA IV", 4),
    ]
    assert impairment.HEART_FAILURE_COLORS["Normal cardiac function"] == "black"


def test_old_maps_unchanged() -> None:
    """The maps of the models not built on the library keep their values."""
    assert dict(impairment.CIRRHOSIS) == {
        "Control": 0,
        "Mild cirrhosis": 0.3994897959183674,
        "Moderate cirrhosis": 0.6979591836734694,
        "Severe cirrhosis": 0.8127551020408164,
    }
    assert dict(impairment.CARDIAC) == {
        "Normal cardiac function": 1.0,
        "Mild cardiac impairment": 4.75 / 5.25,
        "Moderate cardiac impairment": 3.5 / 5.25,
        "Severe cardiac impairment": 3 / 5.25,
        "Cardiac failure": 0.5,
    }
    assert impairment.RENAL_FUNCTION["Moderate renal impairment"] == 32.5 / 101.0
    assert impairment.RENAL_FUNCTION_GFR["End stage renal disease"] == 10.5 / 101.0
    assert impairment.RENAL_FUNCTION_90["End stage renal disease"] == 9.0 / 90.0


def test_the_labels_and_units_of_the_pk_parameters() -> None:
    assert list(pk.PK_LABELS.items()) == [
        ("auc_last", "AUCend"),
        ("auc_inf_obs", "AUC"),
        ("cl_f", "Apparent clearance"),
        ("cl_renal", "Renal clearance"),
        ("cl_fecal", "Fecal clearance"),
        ("cmax", "Cmax"),
        ("thalf", "Half-life"),
        ("lambda_z", "Elimination rate"),
        ("vz_f", "Apparent volume"),
    ]
    assert list(pk.PK_UNITS.items()) == [
        ("auc_last", "µmole/l*hr"),
        ("auc_inf_obs", "µmole/l*hr"),
        ("cl_f", "ml/min"),
        ("cl_renal", "ml/min"),
        ("cl_fecal", "ml/min"),
        ("cmax", "µmole/l"),
        ("thalf", "hr"),
        ("lambda_z", "1/hr"),
        ("vz_f", "l"),
    ]


@pytest.mark.parametrize("classification", list(impairment.CLASSIFICATIONS.values()))
def test_the_classes_of_a_classification(
    classification: impairment.Classification,
) -> None:
    """The classes are unique, ordered by their severity and look themselves up."""
    keys = classification.keys()
    assert len(set(keys)) == len(keys) == len(classification)
    assert len({c.label for c in classification}) == len(keys)
    values = [c.value for c in classification if c.value is not None]
    assert values == sorted(values, reverse=classification.parameter == "eGFR")
    for c in classification:
        assert classification[c.key] is c
        assert c.definition
    assert impairment.CLASSIFICATIONS[classification.key] is classification
    assert list(classification.colors().values()) == [c.color for c in classification]
    with pytest.raises(KeyError, match="nope"):
        classification["nope"]


def test_the_classes_keep_the_colors_of_the_maps() -> None:
    """The Child-Pugh and NYHA classes have the colors of the severity maps."""
    for classification, severity, colors in (
        (
            impairment.CHILD_PUGH_CLASSES,
            impairment.CHILD_PUGH,
            impairment.CIRRHOSIS_COLORS,
        ),
        (
            impairment.NYHA_CLASSES,
            impairment.HEART_FAILURE,
            impairment.HEART_FAILURE_COLORS,
        ),
    ):
        assert [c.value for c in classification] == list(severity.values())
        assert [c.color for c in classification] == list(colors.values())


def test_the_gfr_categories_tile_the_egfr() -> None:
    """The limits of the KDIGO categories meet, the value of a class within them."""
    classes = list(impairment.KDIGO_GFR_CLASSES)
    assert classes[0].limits == (90.0, None)
    assert classes[-1].limits == (None, 15.0)
    for upper, lower in itertools.pairwise(classes):
        assert upper.limits is not None and lower.limits is not None
        assert lower.limits[1] == upper.limits[0]
    for c in classes[1:]:
        assert c.limits is not None and c.value is not None
        lo, hi = c.limits
        assert (lo or 0.0) <= c.value < hi  # ty: ignore[unsupported-operator]
