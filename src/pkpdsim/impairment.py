"""The classes of impaired organ function: the factors of their parameters and colors.

The models scale a parameter of an organ by these factors, the renal function
by `KI__f_renal_function`, the cirrhosis by `f_cirrhosis` (the degree of the
liver which is cirrhotic) and the cardiac function by the cardiac output; the
colors draw a class in every figure of the model alike. The base of the
experiments of a model (`pkpdsim.experiment.PKPDSimulationExperiment`) has
`RENAL_FUNCTION`, `CIRRHOSIS` and `CARDIAC` as its maps; a model which defines
the classes otherwise assigns another preset in its base.

The models built on the PBPK library (`pkpdsim.pbpk.impairment`) have other
parameters: `CP_severity` (liver cirrhosis, 0 healthy and 1, 2, 3 the
Child-Pugh classes A, B, C), `f_renal_function` (the GFR relative to healthy,
which gives the stage of chronic kidney disease) and `HF_severity` (heart
failure, 0 healthy and 1 to 4 the NYHA classes I to IV). Their maps are
`CHILD_PUGH` (colors `CIRRHOSIS_COLORS`), the `RENAL_FUNCTION` presets and
`HEART_FAILURE` (colors `HEART_FAILURE_COLORS`); the maps of the other models
keep their values and meaning.

The clinical classes of the conditions are `Classification`s with their
definitions, colors and the value of the library parameter of every class:
`CHILD_PUGH_CLASSES` (Child-Pugh, `CP_severity`), `KDIGO_GFR_CLASSES` (the GFR
categories of KDIGO 2012, the eGFR of the library), `NYHA_CLASSES` (the NYHA
functional classes, `HF_severity`), `AGE_CLASSES` (the paediatric and
geriatric age classes of ICH E11 and ICH E7, `AGE`) and `BMI_CLASSES` (the
body mass index classes of the WHO, `BMI`), together `CLASSIFICATIONS`; a
model scans and draws them with `class_simulations` and `add_classes` of
`pkpdsim.pbpk.PBPKSimulationExperiment` (docs page `pbpk/classes.md`).
"""

import itertools
from collections.abc import Iterator, Mapping
from dataclasses import dataclass

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

#: the severity `CP_severity` of the Child-Pugh classes of the models built on
#: the PBPK library (`pkpdsim.pbpk.impairment`); colors `CIRRHOSIS_COLORS`
CHILD_PUGH: Mapping[str, float] = {
    "Control": 0,
    "Mild cirrhosis": 1,  # CPT A
    "Moderate cirrhosis": 2,  # CPT B
    "Severe cirrhosis": 3,  # CPT C
}

#: the severity `HF_severity` of the NYHA classes of heart failure of the models
#: built on the PBPK library (`pkpdsim.pbpk.impairment`)
HEART_FAILURE: Mapping[str, float] = {
    "Normal cardiac function": 0,
    "NYHA I": 1,
    "NYHA II": 2,
    "NYHA III": 3,
    "NYHA IV": 4,
}

#: the colors of the classes of `HEART_FAILURE`
HEART_FAILURE_COLORS: Mapping[str, str] = {
    "Normal cardiac function": "black",
    "NYHA I": "#FDDBC7",
    "NYHA II": "#F4A582",
    "NYHA III": "#D6604D",
    "NYHA IV": "#67001F",
}


@dataclass(frozen=True)
class ImpairmentClass:
    """A clinical class of an impaired organ function.

    Attributes:
        key: the key of the class, e.g. `A`, `G3a`, `III`; `healthy` for the
            reference without the condition.
        label: the label of the class in a legend.
        color: the color of the class in every figure.
        definition: the clinical definition of the class.
        value: the value of the parameter of its classification which
            represents the class in a model (`CP_severity`, `HF_severity`, or
            the eGFR in ml/min/1.73 m2 for the GFR categories); `None` for the
            reference of the model (the renal function 1 of the build).
        limits: the range of the quantity which defines the class (the
            Child-Pugh score, the eGFR), `(lower, upper)` with an open end
            `None`; `None` for a class defined by symptoms.
    """

    key: str
    label: str
    color: str
    definition: str
    value: float | None
    limits: tuple[float | None, float | None] | None = None


@dataclass(frozen=True)
class Classification:
    """The clinical classes of a condition with the parameter which sets them.

    Attributes:
        key: the key of the classification, e.g. `child_pugh`.
        name: the name, e.g. `Child-Pugh classification`.
        condition: the condition the classes grade.
        parameter: the parameter of the models built on the PBPK library whose
            value represents a class (`CP_severity`, `HF_severity`, or `eGFR`
            for the GFR categories, which `f_renal_function` sets).
        quantity: the quantity of the `limits` of the classes and its unit.
        reference: the key of the publication of the definition in the sources
            of the parameter database (`pkpdsim/database/sources.toml`).
        classes: the classes in the order of the severity, the reference
            without the condition first.
    """

    key: str
    name: str
    condition: str
    parameter: str
    quantity: str
    reference: str
    classes: tuple[ImpairmentClass, ...]

    def __iter__(self) -> Iterator[ImpairmentClass]:
        """Iterate over the classes in the order of the severity."""
        return iter(self.classes)

    def __len__(self) -> int:
        """The number of classes."""
        return len(self.classes)

    def __getitem__(self, key: str) -> ImpairmentClass:
        """Get a class by its key.

        Raises:
            KeyError: for an unknown key, naming the keys.
        """
        for c in self.classes:
            if c.key == key:
                return c
        raise KeyError(
            f"No class '{key}' of {self.name}, the classes are {self.keys()}."
        )

    def keys(self) -> list[str]:
        """The keys of the classes in the order of the severity."""
        return [c.key for c in self.classes]

    def colors(self) -> dict[str, str]:
        """The colors of the classes by label, for a legend."""
        return {c.label: c.color for c in self.classes}

    def color_of(self, value: float) -> str:
        """Get the color of a continuous value of the parameter of the classes.

        The color is interpolated in RGB between the colors of the two classes
        whose values enclose it, so a continuous scan is drawn in the colors
        of its classes; a class without a value (G1, the reference of the
        model) is at the lower limit of its range, a value beyond the first
        or the last class has its color.

        Args:
            value: the value, e.g. `CP_severity` 1.5 or an eGFR of 40.

        Returns:
            The color as `#rrggbb`.
        """
        points = []
        for c in self.classes:
            x = c.value
            if x is None and c.limits is not None:
                x = c.limits[0] if c.limits[0] is not None else c.limits[1]
            if x is not None:
                points.append((x, _rgb(c.color)))
        points.sort(key=lambda p: p[0])
        if value <= points[0][0]:
            return _hex(points[0][1])
        for (x0, c0), (x1, c1) in itertools.pairwise(points):
            if value <= x1:
                t = (value - x0) / (x1 - x0)
                return _hex(tuple(a + t * (b - a) for a, b in zip(c0, c1, strict=True)))
        return _hex(points[-1][1])


_NAMED = {"black": "#000000"}


def _rgb(color: str) -> tuple[float, float, float]:
    """The RGB of a color `#rrggbb` (or `black`), 0 to 255."""
    code = _NAMED.get(color, color).lstrip("#")
    return (int(code[0:2], 16), int(code[2:4], 16), int(code[4:6], 16))


def _hex(rgb: tuple[float, ...]) -> str:
    """The `#rrggbb` of an RGB."""
    return "#" + "".join(f"{round(v):02x}" for v in rgb)


CHILD_PUGH_CLASSES = Classification(
    key="child_pugh",
    name="Child-Pugh classification",
    condition="liver cirrhosis",
    parameter="CP_severity",
    quantity="Child-Pugh score (points)",
    reference="Pugh1973",
    classes=(
        ImpairmentClass(
            "healthy",
            "Healthy",
            "black",
            "no liver disease",
            0.0,
        ),
        ImpairmentClass(
            "A",
            "Child-Pugh A",
            "#74a9cf",
            "well-compensated cirrhosis, score 5-6",
            1.0,
            (5.0, 6.0),
        ),
        ImpairmentClass(
            "B",
            "Child-Pugh B",
            "#2b8cbe",
            "significant functional compromise, score 7-9",
            2.0,
            (7.0, 9.0),
        ),
        ImpairmentClass(
            "C",
            "Child-Pugh C",
            "#045a8d",
            "decompensated cirrhosis, score 10-15",
            3.0,
            (10.0, 15.0),
        ),
    ),
)
"""The Child-Pugh classes of liver cirrhosis (`CP_severity` 0 to 3)."""

KDIGO_GFR_CLASSES = Classification(
    key="kdigo_gfr",
    name="KDIGO GFR categories",
    condition="chronic kidney disease",
    parameter="eGFR",
    quantity="eGFR (ml/min/1.73 m2)",
    reference="KDIGO2012",
    classes=(
        ImpairmentClass(
            "G1",
            "G1 normal",
            "black",
            "normal or high, eGFR >= 90",
            None,
            (90.0, None),
        ),
        ImpairmentClass(
            "G2",
            "G2 mild",
            "#66c2a4",
            "mildly decreased, eGFR 60-89",
            75.0,
            (60.0, 90.0),
        ),
        ImpairmentClass(
            "G3a",
            "G3a mild-moderate",
            "#41ae76",
            "mildly to moderately decreased, eGFR 45-59",
            52.0,
            (45.0, 60.0),
        ),
        ImpairmentClass(
            "G3b",
            "G3b moderate-severe",
            "#238b45",
            "moderately to severely decreased, eGFR 30-44",
            37.0,
            (30.0, 45.0),
        ),
        ImpairmentClass(
            "G4",
            "G4 severe",
            "#006d2c",
            "severely decreased, eGFR 15-29",
            22.5,
            (15.0, 30.0),
        ),
        ImpairmentClass(
            "G5",
            "G5 kidney failure",
            "#00441b",
            "kidney failure, eGFR < 15 (not on dialysis in the library)",
            7.5,
            (None, 15.0),
        ),
    ),
)
"""The GFR categories of chronic kidney disease (KDIGO 2012), by the eGFR in
ml/min/1.73 m2 at the midpoint of each category (G5 at 7.5, the end-stage
anchor of the library); G1 is the renal function of the reference individual."""

NYHA_CLASSES = Classification(
    key="nyha",
    name="NYHA functional classification",
    condition="chronic heart failure",
    parameter="HF_severity",
    quantity="symptoms at physical activity",
    reference="NYHA1994",
    classes=(
        ImpairmentClass("healthy", "Healthy", "black", "no heart failure", 0.0),
        ImpairmentClass(
            "I",
            "NYHA I",
            "#FDDBC7",
            "no limitation of physical activity; ordinary activity causes no "
            "undue fatigue, palpitation or dyspnoea",
            1.0,
        ),
        ImpairmentClass(
            "II",
            "NYHA II",
            "#F4A582",
            "slight limitation; comfortable at rest, ordinary activity causes "
            "fatigue, palpitation or dyspnoea",
            2.0,
        ),
        ImpairmentClass(
            "III",
            "NYHA III",
            "#D6604D",
            "marked limitation; comfortable at rest, less than ordinary activity "
            "causes symptoms",
            3.0,
        ),
        ImpairmentClass(
            "IV",
            "NYHA IV",
            "#67001F",
            "unable to carry on any physical activity without discomfort; "
            "symptoms at rest",
            4.0,
        ),
    ),
)
"""The NYHA functional classes of chronic heart failure (`HF_severity` 0 to 4)."""

_DAY = 1.0 / 365.25
"""A day in years."""

AGE_CLASSES = Classification(
    key="age",
    name="Age classes (ICH E11, ICH E7)",
    condition="age",
    parameter="AGE",
    quantity="age (years)",
    reference="ICHE11",
    classes=(
        ImpairmentClass(
            "newborn",
            "Term newborn",
            "#c6dbef",
            "term newborn infants, 0 to 27 days (ICH E11); preterm newborn "
            "infants need the gestational age, which the library does not model",
            13.5 * _DAY,
            (0.0, 28 * _DAY),
        ),
        ImpairmentClass(
            "infant",
            "Infant and toddler",
            "#6baed6",
            "infants and toddlers, 28 days to 23 months (ICH E11)",
            1.0,
            (28 * _DAY, 2.0),
        ),
        ImpairmentClass(
            "child",
            "Child",
            "#2171b5",
            "children, 2 to 11 years (ICH E11)",
            6.5,
            (2.0, 12.0),
        ),
        ImpairmentClass(
            "adolescent",
            "Adolescent",
            "#08306b",
            "adolescents, 12 to 17 years (ICH E11: 12 to 16-18 years by "
            "region; the EMA guideline on paediatric pharmacokinetics 12 to 17)",
            15.0,
            (12.0, 18.0),
        ),
        ImpairmentClass(
            "adult",
            "Adult",
            "black",
            "adults, 18 to 64 years; the value of the class is the adult of "
            "ICRP 89 (35 years)",
            35.0,
            (18.0, 65.0),
        ),
        ImpairmentClass(
            "young_old",
            "65-74 years",
            "#fdae6b",
            "geriatric patients of 65 to 74 years (ICH E7 Questions and "
            "Answers 2010; ICH E7: the geriatric population is 65 years or older)",
            70.0,
            (65.0, 75.0),
        ),
        ImpairmentClass(
            "old_old",
            "75-84 years",
            "#e6550d",
            "geriatric patients of 75 to 84 years (ICH E7 Questions and Answers 2010)",
            80.0,
            (75.0, 85.0),
        ),
        ImpairmentClass(
            "oldest_old",
            "85 years and older",
            "#7f2704",
            "geriatric patients of 85 years and older (ICH E7 Questions and "
            "Answers 2010); the physiology beyond 85 years is extrapolated",
            90.0,
            (85.0, None),
        ),
    ),
)
"""The age classes of the paediatric (ICH E11) and geriatric (ICH E7)
populations (`AGE` in years); the value of a class is a representative age,
the physiology is continuous in age."""

BMI_CLASSES = Classification(
    key="bmi",
    name="Body mass index classes (WHO)",
    condition="overweight and obesity",
    parameter="BMI",
    quantity="body mass index (kg/m2)",
    reference="WHO2000",
    classes=(
        ImpairmentClass(
            "normal",
            "Normal weight",
            "black",
            "a body mass index of 18.5 to 24.9 kg/m2 (WHO Technical Report "
            "Series 894); the class has no value of its own, it is the "
            "reference individual of the build (a body mass index of 23.6 "
            "kg/m2 for the man and 22.6 for the woman of ICRP 89)",
            None,
            (18.5, 25.0),
        ),
        ImpairmentClass(
            "underweight",
            "Underweight",
            "#6baed6",
            "a body mass index below 18.5 kg/m2 (WHO: severe below 16, "
            "moderate 16 to 16.9, mild 17 to 18.4); the value of the class is "
            "the middle of the mild and moderate range",
            17.0,
            (None, 18.5),
        ),
        ImpairmentClass(
            "overweight",
            "Overweight",
            "#fdd0a2",
            "a body mass index of 25.0 to 29.9 kg/m2, pre-obese in WHO "
            "Technical Report Series 894; the value of the class is its middle",
            27.5,
            (25.0, 30.0),
        ),
        ImpairmentClass(
            "obesity_1",
            "Obesity class I",
            "#fdae6b",
            "a body mass index of 30.0 to 34.9 kg/m2; the value of the class "
            "is its middle",
            32.5,
            (30.0, 35.0),
        ),
        ImpairmentClass(
            "obesity_2",
            "Obesity class II",
            "#e6550d",
            "a body mass index of 35.0 to 39.9 kg/m2; the value of the class "
            "is its middle",
            37.5,
            (35.0, 40.0),
        ),
        ImpairmentClass(
            "obesity_3",
            "Obesity class III",
            "#7f2704",
            "a body mass index of 40.0 kg/m2 and above, severe or morbid "
            "obesity; the value of the class is 45 kg/m2, and the adipose "
            "model of the library is valid to about 40 (the equations of the "
            "fat-free mass and of the obese repository both saturate above it)",
            45.0,
            (40.0, None),
        ),
    ),
)
"""The classes of the body mass index of the WHO classification (`BMI`, the
body weight of the class at the height of the reference individual). The WHO
Expert Consultation 2004 retained these cut-offs as the international
classification and added the public health action points 23.0, 27.5, 32.5 and
37.5 kg/m2 for Asian populations, which the percent body fat of Gallagher 2000
reproduces (a body mass index of 22.7 in Asian men has the body fat of 25.0 in
white men); the library carries the action points in `docs/pbpk/adipose.md`
and not as a classification of its own."""

CLASSIFICATIONS: Mapping[str, Classification] = {
    c.key: c
    for c in (
        CHILD_PUGH_CLASSES,
        KDIGO_GFR_CLASSES,
        NYHA_CLASSES,
        AGE_CLASSES,
        BMI_CLASSES,
    )
}
"""The classifications by key."""
