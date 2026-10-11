"""The classes of impaired organ function: the factors of their parameters and colors.

The models scale a parameter of an organ by these factors, the renal function
by `KI__f_renal_function`, the cirrhosis by `f_cirrhosis` (the degree of the
liver which is cirrhotic) and the cardiac function by the cardiac output; the
colors draw a class in every figure of the model alike. The base of the
experiments of a model (`pkpdsim.experiment.PKPDSimulationExperiment`) has
`RENAL_FUNCTION`, `CIRRHOSIS` and `CARDIAC` as its maps; a model which defines
the classes otherwise assigns another preset in its base.
"""

from collections.abc import Mapping

#: the renal function by the midpoints of the classes of the glomerular
#: filtration rate, relative to a normal rate of 101 ml/min
RENAL_FUNCTION: Mapping[str, float] = {
    "Normal renal function": 101.0 / 101.0,
    "Mild renal impairment": 69.5 / 101.0,
    "Moderate renal impairment": 32.5 / 101.0,
    "Severe renal impairment": 19.5 / 101.0,
}

#: the renal function by the limits of the classes of the glomerular filtration
#: rate (50, 35 and 20 ml/min) and the end stage renal disease, relative to a
#: normal rate of 101 ml/min
RENAL_FUNCTION_GFR: Mapping[str, float] = {
    "Normal renal function": 101.0 / 101.0,
    "Mild renal impairment": 50.0 / 101.0,
    "Moderate renal impairment": 35.0 / 101.0,
    "Severe renal impairment": 20.0 / 101.0,
    "End stage renal disease": 10.5 / 101.0,
}

#: the renal function by the midpoints of the classes of the glomerular
#: filtration rate and the end stage renal disease, relative to a normal rate
#: of 90 ml/min
RENAL_FUNCTION_90: Mapping[str, float] = {
    "Normal renal function": 90.0 / 90.0,
    "Mild renal impairment": (60 + 89) / 2 / 90.0,
    "Moderate renal impairment": (30 + 59) / 2 / 90.0,
    "Severe renal impairment": (15 + 29) / 2 / 90.0,
    "End stage renal disease": 9.0 / 90.0,
}

#: the colors of the classes of `RENAL_FUNCTION`
RENAL_COLORS: Mapping[str, str] = {
    "Normal renal function": "black",
    "Mild renal impairment": "#66c2a4",
    "Moderate renal impairment": "#2ca25f",
    "Severe renal impairment": "#006d2c",
}

#: the colors of the classes of `RENAL_FUNCTION_GFR` and `RENAL_FUNCTION_90`
RENAL_COLORS_END_STAGE: Mapping[str, str] = {
    **RENAL_COLORS,
    "End stage renal disease": "#006d5e",
}

#: the degree of cirrhosis of the Child-Pugh classes A, B and C
CIRRHOSIS: Mapping[str, float] = {
    "Control": 0,
    "Mild cirrhosis": 0.3994897959183674,  # CPT A
    "Moderate cirrhosis": 0.6979591836734694,  # CPT B
    "Severe cirrhosis": 0.8127551020408164,  # CPT C
}

#: the colors of the classes of `CIRRHOSIS`
CIRRHOSIS_COLORS: Mapping[str, str] = {
    "Control": "black",
    "Mild cirrhosis": "#74a9cf",
    "Moderate cirrhosis": "#2b8cbe",
    "Severe cirrhosis": "#045a8d",
}

#: the cardiac output of the classes of heart failure, relative to a normal
#: output of 5.25 l/min
CARDIAC: Mapping[str, float] = {
    "Normal cardiac function": 1.0,
    "Mild cardiac impairment": 4.75 / 5.25,
    "Moderate cardiac impairment": 3.5 / 5.25,
    "Severe cardiac impairment": 3 / 5.25,
    "Cardiac failure": 0.5,
}

#: the colors of the classes of `CARDIAC`, and of the congestive heart failure
#: of the data
CARDIAC_COLORS: Mapping[str, str] = {
    "Normal cardiac function": "black",
    "Mild cardiac impairment": "#FDDBC7",
    "Moderate cardiac impairment": "#D6604D",
    "Severe cardiac impairment": "#67001F",
    "Cardiac failure": "#67001F",
    "CHF": "#67001F",
}
