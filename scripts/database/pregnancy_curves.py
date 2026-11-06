"""Write `src/pkpdsim/database/tables/pregnancy_human.toml` from the published fits.

The maternal physiology over the gestational age GA (weeks since the first day
of the last menstrual period; 0 not pregnant, 40 term) as factors relative to
the non-pregnant woman, evaluated on the grid `GRID` and interpolated linearly
by the library (`pkpdsim.pbpk.pregnancy`). The fits of Dallmann 2017 and 2018
are of the fertilization age `GA - 2` weeks (Dallmann 2017, ESM Tab. 1
footnote), so every factor is 1 up to GA 2; the factors of the enzymes and
transporters are the in vivo studies against postpartum, placed at the middle
of their windows.

Run from the repository: `uv run python scripts/database/pregnancy_curves.py`.
"""

import math
from collections.abc import Callable
from itertools import pairwise
from pathlib import Path

GRID: tuple[float, ...] = tuple(float(g) for g in range(0, 42, 2))
"""The gestational ages (weeks) of the table."""

OUT = (
    Path(__file__).parents[2]
    / "src"
    / "pkpdsim"
    / "database"
    / "tables"
    / "pregnancy_human.toml"
)

Entry = tuple[Callable[[float], float], str, tuple[str, ...], str]
"""A row of the table: the curve of GA, its source, further sources, the notes."""


def fa(ga: float) -> float:
    """The fertilization age (weeks), 0 before fertilization."""
    return max(ga - 2.0, 0.0)


def points(pairs: list[tuple[float, float]]) -> Callable[[float], float]:
    """The piecewise linear interpolation of points of GA, constant outside them."""

    def value(ga: float) -> float:
        if ga <= pairs[0][0]:
            return pairs[0][1]
        for (x0, y0), (x1, y1) in pairwise(pairs):
            if ga <= x1:
                return y0 + (y1 - y0) * (ga - x0) / (x1 - x0)
        return pairs[-1][1]

    return value


def ratio(f: Callable[[float], float]) -> Callable[[float], float]:
    """A fit of the fertilization age as the factor relative to GA 0."""
    return lambda ga: f(fa(ga)) / f(0.0)


# --- the fits of Dallmann 2017 (ESM 1, Tab. 1), x the fertilization age ----


def cardiac_output(x: float) -> float:
    """Cardiac output (l/min), Eq. 19 (natural logarithm)."""
    if x <= 0:
        return 5.14
    return 6.09 * math.exp(-math.exp(-0.352 * math.log(x) + 1.36)) + 5.14


def hematocrit(x: float) -> float:
    """Hematocrit (%), Eq. 21."""
    return 0.000401 * x**3 - 0.0180 * x**2 + 0.0299 * x + 40.1


def plasma_volume(x: float) -> float:
    """Plasma volume (l), Eq. 29."""
    return -0.0000817 * x**3 + 0.00395 * x**2 + 0.0000841 * x + 2.36


def blood_volume(x: float) -> float:
    """Blood volume (l), the plasma volume over `1 - HCT`."""
    return plasma_volume(x) / (1 - hematocrit(x) / 100)


def gfr(x: float) -> float:
    """GFR (ml/min), Eq. 27 with its misprints read as 0.0245, 113 and 150."""
    if x <= 14:
        return -0.000294 * x**4 + 0.0245 * x**3 - 0.716 * x**2 + 8.73 * x + 113
    return 150.0


def erpf(x: float) -> float:
    """Effective renal plasma flow (l/min), Eq. 25."""
    return -0.00000274 * x**4 + 0.000223 * x**3 - 0.00644 * x**2 + 0.0746 * x + 0.534


def renal_blood_flow(x: float) -> float:
    """Renal blood flow (l/min), the plasma flow over `1 - HCT`."""
    return erpf(x) / (1 - hematocrit(x) / 100)


def uterine_flow(x: float) -> float:
    """Bilateral uterine arterial blood flow (l/min), Eq. 64."""
    return 0.954 / (1 + 21.8 * math.exp(-0.153 * x)) - 0.0197


def albumin(x: float) -> float:
    """Serum albumin (g/l), Eq. 9."""
    return 14.7 * math.exp(-0.0454 * x) + 31.7


def agp(x: float) -> float:
    """Alpha-1-acid glycoprotein (g/l), Eq. 11."""
    return 0.0000768 * x**2 - 0.00573 * x + 0.701


def fetus(x: float) -> float:
    """Fetal volume (l), ESM Tab. 1."""
    return 0.0 if x <= 0 else 5.90 / (1 + math.exp(-4.45 * math.log(x) + 15.8))


def placenta(x: float) -> float:
    """Placental volume (l), ESM Tab. 1."""
    return max(0.937 * math.exp(-6.10 * math.exp(-0.0813 * x)) - 0.00211, 0.0)


def amniotic_fluid(x: float) -> float:
    """Amniotic fluid (l), ESM Tab. 1 (the geometric mean)."""
    if x <= 0:
        return 0.0
    return math.exp(-0.0000169 * x**4 + 0.00169 * x**3 - 0.0653 * x**2 + 1.22 * x - 10)


def uterus(x: float) -> float:
    """Uterus (l), Eq. 58."""
    return 0.959 / (1 + 28.9 * math.exp(-0.177 * x)) + 0.048


def breasts(x: float) -> float:
    """Breast tissue (l), Eq. 35."""
    return math.exp(0.0110 * x - 0.348)


# --- the derived factors of the library ------------------------------------


def placental_fraction(ga: float) -> float:
    """The blood flow of the placenta as fraction of the cardiac output.

    The gain of the uterine arterial flow over the non-pregnant one; the
    organs of the body share the rest of the cardiac output.
    """
    x = fa(ga)
    return (uterine_flow(x) - uterine_flow(0.0)) / cardiac_output(x)


def kidney_flow(ga: float) -> float:
    """The factor on the flow fraction of the kidney.

    The kidney keeps the measured renal blood flow while the organs share the
    cardiac output less the flow of the placenta.
    """
    organs = ratio(cardiac_output)(ga) * (1 - placental_fraction(ga))
    return ratio(renal_blood_flow)(ga) / organs


def cyp1a2(ga: float) -> float:
    """CYP1A2, Dallmann 2018 Eq. 2 of the fertilization age."""
    x = fa(ga)
    return 1 + (0.0291 * x**2 - 2.77 * x) / 100


def ugt1a4(ga: float) -> float:
    """UGT1A4, the lamotrigine clearance of Karanam 2025 (Table 2)."""
    return 1 + (4.9 / 2.79) * ga / (12 + ga) if ga > 2 else 1.0


FACTORS: dict[str, Entry] = {
    "cardiac_output": (
        ratio(cardiac_output),
        "Dallmann2017",
        (),
        "Eq. 19, the meta-analysed cardiac output (+40 % at term); the total "
        "cardiac output, of which the placenta takes `placental_fraction`",
    ),
    "placental_fraction": (
        placental_fraction,
        "Dallmann2017",
        (),
        "not a factor but a fraction, 0 at GA 0: the gain of the bilateral "
        "uterine arterial flow over the non-pregnant one (Eq. 64, 0.88 l/min at "
        "term) as fraction of the cardiac output (Eq. 19), the blood flow of the "
        "placenta",
    ),
    "blood_volume": (
        ratio(blood_volume),
        "Dallmann2017",
        (),
        "derived: the plasma volume (Eq. 29) over 1 - HCT (Eq. 21), +45 % at term",
    ),
    "hematocrit": (
        ratio(hematocrit),
        "Dallmann2017",
        (),
        "Eq. 21, the haemodilution of the pregnancy (-7 % at term)",
    ),
    "gfr": (
        ratio(gfr),
        "Dallmann2017",
        (),
        "Eq. 27, its misprints read as 0.0245, 113 and 150 ml/min (the text of "
        "the article and the tabulated curve of the PK-Sim database): +33 % from "
        "the end of the first trimester",
    ),
    "kidney_flow": (
        kidney_flow,
        "Dallmann2017",
        (),
        "derived: the renal blood flow (Eq. 25 over 1 - HCT) relative to the flow "
        "of the organs (the cardiac output less the placenta), so the kidney "
        "keeps its measured flow and the other organs share the rest; the renal "
        "flow peaks at the end of the first trimester and falls towards term",
    ),
    "albumin": (
        ratio(albumin),
        "Dallmann2017",
        (),
        "Eq. 9, the serum albumin (-26 % at term), which raises fu of an "
        "albumin-bound substance",
    ),
    "agp": (
        ratio(agp),
        "Dallmann2017",
        (),
        "Eq. 11, the alpha-1-acid glycoprotein (-15 % at term)",
    ),
    "t_si": (
        points([(0.0, 1.0), (8.0, 1.0), (20.0, 1.667), (34.0, 1.60)]),
        "Lawson1985",
        ("Chiloiro2001", "Szilagyi1996"),
        "the small intestinal transit time of the orocaecal transit studies "
        "against postpartum: unchanged in the first trimester, 125/75 min in the "
        "second (Lawson 1985), the weighted mean 1.60 of Lawson 1985, Chiloiro "
        "2001 and Szilagyi 1996 in the third; the gastric emptying does not change",
    ),
}
"""The factors of the maternal physiology by the key of their factor."""

PROTEINS: dict[tuple[str, str], Entry] = {
    ("li", "CYP1A2"): (
        cyp1a2,
        "Dallmann2018",
        ("Tracy2005",),
        "Eq. 2 of the meta-analysed caffeine and theophylline data (percent "
        "change 0.0291*FA^2 - 2.77*FA): -37 % at term; the dextromethorphan of "
        "Tracy 2005 (n = 25, within-subject against postpartum) gives 0.35",
    ),
    ("li", "CYP3A4"): (
        points([(0.0, 1.0), (2.0, 1.0), (14.0, 1.6), (40.0, 1.6)]),
        "Dallmann2018",
        ("Tracy2005", "Hebert2008"),
        "the weighted mean of the in vivo studies, +60 % from the first trimester "
        "on (Dallmann 2018 Tab. 1; the midazolam of Hebert 2008 gives 1.72, the "
        "dextromethorphan of Tracy 2005 is constant from 14 weeks); the rise "
        "between fertilization and 14 weeks is a choice of the library",
    ),
    ("li", "CYP2D6"): (
        points([(0.0, 1.0), (2.0, 1.0), (16.0, 1.256), (26.0, 1.348), (38.0, 1.478)]),
        "Tracy2005",
        (),
        "the dextromethorphan of Tracy 2005 at the middle of its windows: +48 % "
        "at term; the metoprolol of Ryu 2016 gives 1.8 to 3.0, not taken here",
    ),
    ("li", "CYP2C9"): (
        points([(0.0, 1.0), (2.0, 1.0), (8.0, 1.4), (20.0, 1.5), (34.0, 1.6)]),
        "Ke2014",
        (),
        "the phenytoin of Tomson 1994 as Ke 2014 use it, 1.4, 1.5 and 1.6 in the "
        "three trimesters; low confidence",
    ),
    ("li", "CYP2C19"): (
        points([(0.0, 1.0), (2.0, 1.0), (20.0, 0.38), (34.0, 0.32)]),
        "Ke2014",
        (),
        "the proguanil to cycloguanil of McGready 2003 as Ke 2014 use it: -62 % "
        "in the second and -68 % in the third trimester, the first not measured; "
        "low confidence",
    ),
    ("li", "CYP2B6"): (
        points([(0.0, 1.0), (2.0, 1.0), (8.0, 1.1), (20.0, 1.4), (34.0, 1.9)]),
        "Ke2014",
        (),
        "the in vitro induction by estradiol scaled by Ke 2014: 1.1, 1.4 and 1.9 "
        "in the three trimesters; low confidence (no probe drug)",
    ),
    ("li", "CYP2A6"): (
        points([(0.0, 1.0), (2.0, 1.0), (34.0, 1.82)]),
        "Dallmann2018",
        ("Dempsey2002",),
        "the nicotine and cotinine of the third trimester (Dallmann 2018 Tab. 1, "
        "the mean of +109 % and +54 %); only the third trimester is measured",
    ),
    ("li", "UGT1A4"): (
        ugt1a4,
        "Karanam2025",
        (),
        "the lamotrigine clearance of Karanam 2025 (the Emax model of its Table 2 "
        "of the 91 % subpopulation): +135 % at term",
    ),
    ("li", "UGT1A1"): (
        points([(0.0, 1.0), (2.0, 1.0), (12.0, 1.4), (40.0, 1.6)]),
        "Fischer2014",
        (),
        "the labetalol clearance against postpartum (Fischer 2014): 1.4 at 12 "
        "weeks and 1.6 at term; its glucuronidation is UGT1A1 and UGT2B7 and the "
        "factor is given to UGT1A1 alone (low confidence)",
    ),
    ("ki", "ABCB1"): (
        points([(0.0, 1.0), (2.0, 1.0), (30.0, 1.879)]),
        "Hebert2008",
        (),
        "the unbound renal secretion clearance of digoxin at 28-32 weeks against "
        "6-10 weeks postpartum (Hebert 2008): +88 %",
    ),
    ("ki", "SLC22A2"): (
        points([(0.0, 1.0), (2.0, 1.0), (20.0, 1.534), (34.0, 1.339)]),
        "Eyal2010",
        (),
        "the net renal secretion clearance of metformin against postpartum (Eyal "
        "2010): +53 % in the second and +34 % in the third trimester (OCT2 with "
        "MATE)",
    ),
    ("ki", "SLC22A6"): (
        points([(0.0, 1.0), (2.0, 1.0), (27.0, 1.3)]),
        "Thakur2026",
        (),
        "the biomarkers of OAT1 and OAT3 at 25-29 weeks (Thakur 2026): +30 %",
    ),
    ("ki", "SLC22A8"): (
        points([(0.0, 1.0), (2.0, 1.0), (27.0, 1.3)]),
        "Thakur2026",
        (),
        "the biomarkers of OAT1 and OAT3 at 25-29 weeks (Thakur 2026): +30 %",
    ),
}
"""The factors of the enzymes and transporters by organ token and protein."""

MASSES: dict[str, Entry] = {
    "fetus": (
        lambda ga: fetus(fa(ga)),
        "Dallmann2017",
        (),
        "the fetal volume of the ESM Tab. 1, 3.5 kg at term: the compartment of "
        "the fetus",
    ),
    "placenta": (
        lambda ga: placenta(fa(ga)),
        "Dallmann2017",
        (),
        "the placenta of the ESM Tab. 1, 0.71 kg at term: the compartment of the "
        "placenta",
    ),
    "amniotic_fluid": (
        lambda ga: amniotic_fluid(fa(ga)),
        "Dallmann2017",
        (),
        "the amniotic fluid of the ESM Tab. 1 (the geometric mean), 0.65 l at "
        "term; a mass of the body which takes no drug",
    ),
    "uterus": (
        lambda ga: uterus(fa(ga)) - uterus(0.0),
        "Dallmann2017",
        (),
        "the gain of the uterus over the non-pregnant one (Eq. 58), 0.90 kg at "
        "term; a mass of the body which takes no drug",
    ),
    "breasts": (
        lambda ga: breasts(fa(ga)) - breasts(0.0),
        "Dallmann2017",
        (),
        "the gain of the breast tissue over the non-pregnant one (Eq. 35), 0.37 "
        "kg at term; a mass of the body which takes no drug",
    ),
}
"""The masses of the pregnancy (kg), 0 at GA 0."""

HEADER = """\
# The maternal physiology over the gestational age GA (weeks since the first day
# of the last menstrual period; 0 not pregnant, 40 term), generated by
# scripts/database/pregnancy_curves.py from the published fits: do not edit.
#
# Every entry is a factor of the female physiology relative to the non-pregnant
# woman, 1 at GA 0 (`[placental_fraction]` a fraction and `[mass.*]` a mass in
# kg, both 0 at GA 0), at the gestational ages of `ga` and linear between them
# (`pkpdsim.pbpk.pregnancy`). The fits of Dallmann 2017 and Dallmann 2018 are of
# the fertilization age GA - 2 weeks, so every factor is 1 up to GA 2; the
# factors of the enzymes and transporters are the in vivo studies against
# postpartum, at the middle of their windows.
#
# No change (documented in docs/pbpk/pregnancy.md): the gastric emptying, the
# gastric pH, UGT2B7 and the intestinal P-gp; no data: the colonic transit and
# the other transporters; not applied: the liver and the kidney volume and the
# body water of Abduljalil 2012 (it counts the fetus).
"""


def section(name: str, entry: Entry, unit: str) -> str:
    """The TOML section of an entry: its values on the grid with their provenance."""
    curve, source, more, notes = entry
    values = ", ".join(repr(round(curve(ga), 6)) for ga in GRID)
    sources = f"sources = {list(more)!r}\n" if more else ""
    return (
        f"[{name}]\n"
        f'unit = "{unit}"\n'
        f'source = "{source}"\n'
        f"{sources}"
        f'notes = "{notes}"\n'
        f"value = [{values}]\n"
    )


def table() -> str:
    """The content of the table of the pregnancy."""
    parts = [HEADER, f"ga = [{', '.join(repr(g) for g in GRID)}]\n"]
    parts += [section(key, entry, "dimensionless") for key, entry in FACTORS.items()]
    parts += [
        section(f"protein.{organ}.{protein}", entry, "dimensionless")
        for (organ, protein), entry in PROTEINS.items()
    ]
    parts += [section(f"mass.{key}", entry, "kg") for key, entry in MASSES.items()]
    return "\n".join(parts)


if __name__ == "__main__":
    OUT.write_text(table(), encoding="utf-8")
    print(f"written {OUT}")
