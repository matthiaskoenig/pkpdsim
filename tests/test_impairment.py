"""The classes of impaired organ function and the labels of the PK parameters.

The values are those of the models of pkdb_models, written with the same
expressions, so that a model which takes a preset simulates exactly as before.
"""

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
