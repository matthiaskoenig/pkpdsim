"""Write the data tables of the pages of the PBPK library and the database.

The pages of `docs/pbpk/` show the reference physiology, the populations, the
impairment vectors, the arm vein presets, the intestine, the nephron and the
proteins as tables, `docs/database.md` the studies of the database and their
consensus. The numbers are those of the parameter database
`pkpdsim/database` (its tables and the consensus of its studies), never typed
by hand: a page marks the region of a table with

```
<!-- generated: <name> (scripts/pbpk_tables.py), do not edit -->
...
<!-- end generated -->
```

and this script writes the table `<name>` (of `TABLES`) into it, with a line
which says that it is generated and from which file, and the sources the
table cites. Run it after a change of a data file:

```bash
uv run python scripts/pbpk_tables.py
```

`tests/docs/test_pbpk_tables.py` fails while a page is out of date.
"""

import html
import itertools
import re
from collections.abc import Callable, Iterable, Mapping
from pathlib import Path
from typing import Any

from sbmlutils.factory import UnitDefinition

from pkpdsim.database import integration
from pkpdsim.database.studies import KINDS, QUANTITIES, studies
from pkpdsim.impairment import (
    AGE_CLASSES,
    BMI_CLASSES,
    CHILD_PUGH_CLASSES,
    KDIGO_GFR_CLASSES,
    NYHA_CLASSES,
    ImpairmentClass,
)
from pkpdsim.pbpk import ids
from pkpdsim.pbpk.adipose import (
    BINDING_FACTORS,
    CARDIAC_OUTPUT,
    FLOW_FACTORS,
    ObesityFactor,
    obesity_factor,
    obesity_factors,
    obesity_protein_factors,
    visceral_composition,
    visceral_share,
    waist_model,
)
from pkpdsim.pbpk.age import (
    AGE_REFERENCE,
    factor_curves,
    flow_factors,
    flow_fraction_curves,
    protein_factor,
    protein_factors,
    reference_curves,
)
from pkpdsim.pbpk.circadian import circadian_factor, protein_rhythms, rhythms
from pkpdsim.pbpk.data import load, sources
from pkpdsim.pbpk.experiment import PBPKSimulationExperiment
from pkpdsim.pbpk.impairment import (
    _HF_ORGANS,
    child_pugh_protein,
    child_pugh_proteins,
    class_factors,
    meld_anchors,
)
from pkpdsim.pbpk.organs import Organ
from pkpdsim.pbpk.physiology import (
    Sex,
    _references,
    body_composition,
    physiology_table,
    reference_hematocrit,
    reference_individual,
    reference_tissue_volumes,
)
from pkpdsim.pbpk.populations import POPULATIONS
from pkpdsim.pbpk.pregnancy import factors as pregnancy_factors
from pkpdsim.pbpk.pregnancy import pregnancy_factor
from pkpdsim.pbpk.pregnancy import protein_factors as pregnancy_proteins
from pkpdsim.pbpk.proteins import PROTEINS
from pkpdsim.pbpk.sex_differences import (
    protein_sex_factor,
    protein_sex_factors,
    sex_curves,
    sex_factor,
)
from pkpdsim.pbpk.variants.segments import intestine_table
from pkpdsim.pbpk.visual import icons
from pkpdsim.species import Species

REPO_DIR: Path = Path(__file__).parent.parent
PBPK_DIR: Path = REPO_DIR / "docs" / "pbpk"
PAGES: list[Path] = [*sorted(PBPK_DIR.glob("*.md")), REPO_DIR / "docs" / "database.md"]
"""The pages with generated tables."""

END = "<!-- end generated -->"
_REGION = re.compile(
    r"^<!-- generated: (?P<name>[\w-]+) \(scripts/pbpk_tables\.py\), do not edit -->\n"
    r".*?^<!-- end generated -->$",
    re.S | re.M,
)

_ORGAN_NAMES: dict[Organ, str] = {
    Organ.LUNG: "lung",
    Organ.HEART: "heart",
    Organ.BRAIN: "brain",
    Organ.LIVER: "liver",
    Organ.KIDNEY: "kidney",
    Organ.STOMACH: "stomach",
    Organ.SMALL_INTESTINE: "small intestine",
    Organ.LARGE_INTESTINE: "large intestine",
    Organ.SPLEEN: "spleen",
    Organ.PANCREAS: "pancreas",
    Organ.MUSCLE: "muscle",
    Organ.SKIN: "skin",
    Organ.ADIPOSE: "adipose",
    Organ.BONE: "bone",
    Organ.GONADS: "gonads",
}


def begin_marker(name: str) -> str:
    """Get the line which opens the region of a table on a page."""
    return f"<!-- generated: {name} (scripts/pbpk_tables.py), do not edit -->"


def _number(x: float) -> str:
    """Format a number with up to four significant digits."""
    return f"{x:.4g}"


def _row(cells: Iterable[str]) -> str:
    """Format a row of a markdown table."""
    return "| " + " | ".join(cells) + " |"


def _numeric(cell: str) -> bool:
    """Check whether a cell is a number (a model choice marked by `*`) or empty."""
    try:
        float(cell.removesuffix("*") or "0")
    except ValueError:
        return False
    return True


def _table(header: list[str], rows: Iterable[list[str]]) -> list[str]:
    """Format a markdown table, the columns of numbers right aligned."""
    body = list(rows)
    align = [
        "---:" if body and all(_numeric(r[i]) for r in body) else "---"
        for i in range(len(header))
    ]
    return [_row(header), _row(align), *(_row(r) for r in body)]


def _generated(files: str, lines: list[str], cited: Iterable[str]) -> str:
    """Format the content of a region: the notice, the table and its sources."""
    keys = sorted(set(cited))
    known = sources()
    references = []
    for key in keys:
        source = known[key]
        doi = f" [doi:{source.doi}](https://doi.org/{source.doi})" if source.doi else ""
        references.append(f"- `{key}`: {source.citation}{doi}")
    notice = (
        f"*Generated from {files} by `scripts/pbpk_tables.py`; "
        "the notes of every number are in the database.*"
    )
    parts = [notice, "", *lines]
    if references:
        parts += ["", '??? note "Sources"', ""]
        parts += [f"    {r}" for r in references]
    return "\n".join(parts)


def _value(entry: Mapping[str, Any] | None) -> str:
    """Format the value of an entry of a data file, empty if there is none."""
    if entry is None:
        return ""
    value = entry["value"]
    if isinstance(value, list):
        return ", ".join(_number(v) for v in value)
    return _number(value)


def _unit(entries: Iterable[Mapping[str, Any]], key: str) -> str:
    """Get the one unit of a column of a data file, `-` for dimensionless.

    Raises:
        ValueError: if the rows of the column state different units.
    """
    units = {row[key]["unit"] for row in entries if key in row}
    if len(units) != 1:
        raise ValueError(f"column '{key}' has the units {sorted(units)}, not one")
    unit = units.pop()
    return "-" if unit == "dimensionless" else unit


def _unit_id(unit: object) -> str:
    """Get the id of a unit of the model (a `UnitDefinition` or its id)."""
    return str(unit.sid) if isinstance(unit, UnitDefinition) else str(unit)


def _header(label: str, entries: Iterable[Mapping[str, Any]], key: str) -> str:
    """Get the header of a column with the unit of the data file."""
    return f"{label} [{_unit(entries, key)}]"


def _column_sources(
    entries: Iterable[Mapping[str, Any]], columns: Mapping[str, str]
) -> tuple[str, list[str]]:
    """Format the sources of the columns of a table, the columns of one set together.

    Args:
        entries: the tables of the rows of the data file.
        columns: the label of a column by the key of its entry.

    Returns:
        The line which names the sources of every column, and the sources.
    """
    rows = list(entries)
    by_column: dict[str, tuple[str, ...]] = {}
    for key, label in columns.items():
        found = sorted({row[key]["source"] for row in rows if key in row})
        if found:
            by_column[label] = tuple(found)
    groups: dict[tuple[str, ...], list[str]] = {}
    for label, found in by_column.items():
        groups.setdefault(found, []).append(label)
    parts = [
        f"{', '.join(labels)}: {', '.join(f'`{s}`' for s in found)}"
        for found, labels in groups.items()
    ]
    cited = [s for found in groups for s in found]
    return "Sources of the columns: " + "; ".join(parts) + ".", cited


def physiology_general() -> str:
    """The covariates and the whole-body values of the reference individuals."""
    data = load("physiology_human")
    rows: list[list[str]] = []
    cited: list[str] = []
    entries = [
        ("body weight", (), "bodyweight"),
        ("height", (), "height"),
        ("body surface area", (), "body_surface_area"),
        ("cardiac output", (), "cardiac_output"),
        ("hematocrit", (), "hematocrit"),
        ("GFR", (), "gfr"),
        ("hepatic artery share of the liver flow", (), "hepatic_artery_fraction"),
        ("blood volume", ("blood",), "volume"),
        ("blood density", ("blood",), "density"),
        ("arterial pool (fraction of blood)", ("blood",), "arterial"),
        ("venous pool (fraction of blood)", ("blood",), "venous"),
        ("portal vein (fraction of blood)", ("blood",), "portal"),
    ]
    for label, path, key in entries:
        tables = [data[sex.value] for sex in (Sex.MALE, Sex.FEMALE)]
        for part in path:
            tables = [table[part] for table in tables]
        male, female = (table[key] for table in tables)
        rows.append([label, _unit(tables, key), _value(male), _value(female)])
        cited += [male["source"], female["source"]]
    lines = _table(["quantity", "unit", "male", "female"], rows)
    return _generated("`database/tables/physiology_human.toml`", lines, cited)


def _organs(sex: Sex, table: str = "physiology_human") -> str:
    """The organs of the reference individual of a sex of a physiology table."""
    organs = load(table)[sex.value]["organs"]
    rows = [
        [
            f"`{organ.value}` {_ORGAN_NAMES[organ]}",
            _value(organs[organ.value].get("mass")),
            _value(organs[organ.value].get("density", {"value": 1.0})),
            _value(organs[organ.value].get("blood")),
            _value(organs[organ.value].get("f_vas")),
            _value(organs[organ.value].get("f_int")),
            _value(organs[organ.value].get("flow")),
        ]
        for organ in Organ
    ]
    entries = list(organs.values())
    keys = ("mass", "density", "blood", "f_vas", "f_int", "flow")
    header = ["organ", *(_header(k, entries, k) for k in keys)]
    columns = {k: k for k in ("mass", "density", "blood", "f_vas", "f_int", "flow")}
    line, cited = _column_sources(organs.values(), columns)
    lines = [*_table(header, rows), "", line]
    return _generated(f"`database/tables/{table}.toml`", lines, cited)


def physiology_organs_male() -> str:
    """The organs of the male reference individual."""
    return _organs(Sex.MALE)


def physiology_organs_female() -> str:
    """The organs of the female reference individual."""
    return _organs(Sex.FEMALE)


ANIMALS: tuple[Species, ...] = (Species.RAT, Species.MOUSE, Species.PIG)
"""The laboratory animals of the library."""


def species_overview() -> str:
    """The reference animals of the laboratory species against the human."""
    rows: list[list[str]] = []
    cited: list[str] = []
    for species in (Species.HUMAN, *ANIMALS):
        data = load(physiology_table(species))
        for sex in Sex:
            if sex.value not in data:
                continue
            p = data[sex.value]
            refs = _references(sex, data, species)
            adipose = reference_tissue_volumes(sex, species=species)[Organ.ADIPOSE]
            ad = p["organs"]["ad"]
            rho = float(ad.get("density", {"value": 1.0})["value"])
            exponents = (
                ("BSA", "BSA")
                if species is Species.HUMAN
                else (
                    _number(refs[ids.CO_EXPONENT].value),
                    _number(refs[ids.GFR_EXPONENT].value),
                )
            )
            for key in ("bodyweight", "cardiac_output", "gfr", "hematocrit"):
                cited.append(p[key]["source"])
            cited.append(p["blood"]["volume"]["source"])
            rows.append(
                [
                    species.value,
                    sex.value,
                    _number(refs[ids.BW].value),
                    _number(refs[ids.ref(ids.CO)].value * 1000),
                    exponents[0],
                    _number(refs[ids.ref(ids.GFR)].value * 1000),
                    exponents[1],
                    _number(
                        refs[ids.ref(ids.VBLOOD)].value * 1000 / refs[ids.BW].value
                    ),
                    _number(refs[ids.HCT].value),
                    f"{_number(100 * adipose * rho / refs[ids.BW].value)} "
                    f"({_number(100 * float(ad['mass']['value']) / refs[ids.BW].value)})",
                ]
            )
    header = [
        "species",
        "sex",
        "BW [kg]",
        "CO [ml/min]",
        "CO scales with",
        "GFR [ml/min]",
        "GFR scales with",
        "blood [ml/kg]",
        "HCT",
        "adipose [% BW] (reference)",
    ]
    return _generated(
        "`database/tables/physiology_<species>.toml`", _table(header, rows), cited
    )


def species_intestine() -> str:
    """The transit times and fluid volumes of the lumen of the species."""
    rows: list[list[str]] = []
    cited: list[str] = []
    for species in (Species.HUMAN, *ANIMALS):
        data = load(intestine_table(species))
        small = [
            data[k]
            for k in (
                "duodenum",
                "upper_jejunum",
                "lower_jejunum",
                "upper_ileum",
                "lower_ileum",
            )
        ]
        large = [
            data[k]
            for k in (
                "caecum",
                "colon_ascendens",
                "colon_transversum",
                "colon_descendens",
                "colon_sigmoid",
                "rectum",
            )
        ]
        for e in (data["stomach"], *small, *large):
            cited += [e["transit_time"]["source"], e["volume"]["source"]]
        fed = data["stomach"].get("transit_time_fed")
        rows.append(
            [
                species.value,
                _number(data["stomach"]["transit_time"]["value"]),
                _number(fed["value"]) if fed else "",
                _number(sum(e["transit_time"]["value"] for e in small)),
                _number(sum(e["transit_time"]["value"] for e in large) / 60),
                _number(1000 * data["stomach"]["volume"]["value"]),
                _number(1000 * sum(e["volume"]["value"] for e in small)),
                _number(1000 * sum(e["volume"]["value"] for e in large)),
            ]
        )
    header = [
        "species",
        "stomach fasted [min]",
        "stomach fed [min]",
        "small intestine [min]",
        "large intestine [h]",
        "stomach [ml]",
        "small intestine [ml]",
        "large intestine [ml]",
    ]
    return _generated(
        "`database/tables/intestine_<species>.toml`", _table(header, rows), cited
    )


def species_organs_rat() -> str:
    """The organs of the reference rat."""
    return _organs(Sex.MALE, physiology_table(Species.RAT))


def species_organs_mouse() -> str:
    """The organs of the reference mouse."""
    return _organs(Sex.MALE, physiology_table(Species.MOUSE))


def species_organs_pig() -> str:
    """The organs of the reference pig (male)."""
    return _organs(Sex.MALE, physiology_table(Species.PIG))


def circadian_rhythms() -> str:
    """The circadian rhythms of the physiology and the proteins."""
    rows: list[list[str]] = []
    cited: list[str] = []
    entries = [(r, "") for r in rhythms().values()]
    entries += [
        (r, f"{protein} {_ORGAN_NAMES[organ]}")
        for (organ, protein), r in sorted(
            protein_rhythms().items(), key=lambda x: (x[0][0].value, x[0][1])
        )
    ]
    for r, protein in entries:
        cited += list(r.sources)
        hours = int(r.acrophase)
        minutes = round(60 * (r.acrophase - hours))
        rows.append(
            [
                protein or r.key,
                f"`{circadian_factor(r.key)}`",
                _number(r.amplitude),
                f"{hours:02d}:{minutes:02d}",
                _number(r.value(8.0)),
                _number(r.value(20.0)),
                "yes" if r.default else "no",
                ", ".join(f"`{k}`" for k in r.sources),
            ]
        )
    header = [
        "rhythm",
        "factor",
        "amplitude",
        "acrophase",
        "08:00",
        "20:00",
        "default",
        "sources",
    ]
    return _generated(
        "`database/tables/circadian_human.toml`", _table(header, rows), cited
    )


PREGNANCY_AGES: tuple[float, ...] = (0.0, 8.0, 13.0, 20.0, 26.0, 33.0, 40.0)
"""The gestational ages (weeks) of the table of the pregnancy."""


def pregnancy() -> str:
    """The factors of the pregnancy over the gestational age."""
    rows: list[list[str]] = []
    cited: list[str] = []
    entries = [(f.key, f) for f in pregnancy_factors().values()]
    entries += [
        (f"{protein} {_ORGAN_NAMES[organ]}", f)
        for (organ, protein), f in sorted(
            pregnancy_proteins().items(), key=lambda x: (x[0][0].value, x[0][1])
        )
    ]
    for name, f in entries:
        cited += list(f.sources)
        rows.append(
            [
                name.replace("mass_", "mass "),
                f"`{pregnancy_factor(f.key)}`",
                "kg" if f.unit == "kg" else "-",
                *[_number(f.value(ga)) for ga in PREGNANCY_AGES],
                ", ".join(f"`{k}`" for k in f.sources),
            ]
        )
    header = [
        "quantity",
        "factor",
        "unit",
        *[f"{a:g} wk" for a in PREGNANCY_AGES],
        "sources",
    ]
    return _generated(
        "`database/tables/pregnancy_human.toml`", _table(header, rows), cited
    )


ADIPOSE_BMI: tuple[float, ...] = (25.0, 30.0, 35.0, 40.0, 50.0)
"""The body mass indices (kg/m2) of the table of the factors of obesity."""

ADIPOSE_AGES: tuple[float, ...] = (20.0, 30.0, 40.0, 50.0, 60.0, 70.0)
"""The ages (years) of the table of the visceral share of the adipose tissue."""


def adipose_factors() -> str:
    """The factors of overweight and obesity over the body mass index."""
    rows: list[list[str]] = []
    cited: list[str] = []
    male, female = obesity_factors(Sex.MALE), obesity_factors(Sex.FEMALE)
    entries: list[tuple[str, str, str, ObesityFactor, ObesityFactor]] = [
        (
            CARDIAC_OUTPUT.replace("_", " "),
            f"`{obesity_factor(ids.CO)}`",
            "every organ flow",
            male[CARDIAC_OUTPUT],
            female[CARDIAC_OUTPUT],
        )
    ]
    for key, organs in FLOW_FACTORS.items():
        entries.append(
            (
                key.replace("_", " "),
                ", ".join(f"`{obesity_factor(ids.organ_flow(o))}`" for o in organs),
                ", ".join(f"`{ids.organ_flow(o)}`" for o in organs),
                male[key],
                female[key],
            )
        )
    for key in BINDING_FACTORS:
        entries.append(
            (
                key,
                f"`{obesity_factor(key)}`",
                f"`{ids.impairment_factor(key)}`",
                male[key],
                female[key],
            )
        )
    for (organ, protein), f in sorted(
        obesity_protein_factors(Sex.MALE).items(),
        key=lambda x: (x[0][0].value, x[0][1]),
    ):
        entries.append(
            (
                f"{protein} {_ORGAN_NAMES[organ]}",
                f"`{obesity_factor(f'{organ.value}_{protein.lower()}')}`",
                f"`{ids.impairment_activity(protein, organ)}`",
                f,
                obesity_protein_factors(Sex.FEMALE)[organ, protein],
            )
        )
    for name, factor, target, m, f in entries:
        cited += list(m.sources)
        same = all(abs(m.value(b) - f.value(b)) < 5e-4 for b in (*ADIPOSE_BMI, 60.0))
        values = [
            _number(m.value(b))
            if same
            else f"{_number(m.value(b))} / {_number(f.value(b))}"
            for b in ADIPOSE_BMI
        ]
        rows.append(
            [name, factor, target, *values, ", ".join(f"`{k}`" for k in m.sources)]
        )
    header = [
        "quantity",
        "factor",
        "on",
        *[f"BMI {b:g}" for b in ADIPOSE_BMI],
        "sources",
    ]
    return _generated(
        "`database/tables/adipose_human.toml`", _table(header, rows), cited
    )


def adipose_visceral() -> str:
    """The visceral share of the adipose tissue over age."""
    rows: list[list[str]] = []
    cited: list[str] = []
    for sex in Sex:
        share = visceral_share(sex)
        cited += list(share.sources)
        rows.append(
            [
                sex.value,
                *[_number(share.curve.value(a)) for a in ADIPOSE_AGES],
                ", ".join(f"`{k}`" for k in share.sources),
            ]
        )
    header = ["sex", *[f"{a:g} y" for a in ADIPOSE_AGES], "sources"]
    return _generated(
        "`database/tables/adipose_human.toml`", _table(header, rows), cited
    )


def adipose_waist() -> str:
    """The prediction of the waist circumference and its effect on the depot."""
    rows: list[list[str]] = []
    cited: list[str] = []
    for sex in Sex:
        model = waist_model(sex)
        reference = reference_individual(sex)
        cited += list(model.sources)
        rows.append(
            [
                sex.value,
                _number(model.intercept),
                _number(model.bw),
                _number(model.height),
                _number(model.age),
                _number(model.value(reference.bw, reference.height)),
                _number(100.0 * model.sensitivity),
                ", ".join(f"`{k}`" for k in model.sources),
            ]
        )
    header = [
        "sex",
        "intercept [cm]",
        "per kg",
        "per cm",
        "per year",
        "WAIST of the reference individual [cm]",
        "VAT per cm [%]",
        "sources",
    ]
    return _generated(
        "`database/tables/adipose_human.toml`", _table(header, rows), cited
    )


def adipose_classes() -> str:
    """The body composition of the classes of the body mass index."""
    rows: list[list[str]] = []
    order = (
        "underweight",
        "normal",
        "overweight",
        "obesity_1",
        "obesity_2",
        "obesity_3",
    )
    for sex in Sex:
        reference = reference_individual(sex)
        for key in order:
            c = BMI_CLASSES[key]
            bw = (
                reference.bw
                if c.value is None
                else c.value * (reference.height / 100.0) ** 2
            )
            comp = body_composition(sex, bw, reference.height)
            depots = visceral_composition(sex, bw, reference.height)
            rows.append(
                [
                    c.label,
                    sex.value,
                    _number(comp.bmi),
                    _number(comp.bw),
                    _number(comp.ffm),
                    _number(comp.fm),
                    _number(100.0 * comp.fm / comp.bw),
                    _number(comp.fm / (comp.height / 100.0) ** 2),
                    _number(comp.adipose),
                    _number(depots.waist),
                    _number(depots.whtr),
                    _number(depots.visceral),
                    _number(depots.subcutaneous),
                ]
            )
    header = [
        "class",
        "sex",
        "BMI",
        "BW [kg]",
        "FFM [kg]",
        "FM [kg]",
        "%fat",
        "FMI",
        "Mad [kg]",
        "WAIST [cm]",
        "WHtR",
        "VAT [kg]",
        "SAT [kg]",
    ]
    return _generated(
        "`pkpdsim.pbpk.physiology.body_composition` at the height of the "
        "reference individual of the sex (ICRP 89, 35 years)",
        _table(header, rows),
        (),
    )


def populations() -> str:
    """The reference individuals and the covariate distributions of the populations."""
    rows: list[list[str]] = []
    cited: list[str] = []
    default = "european_icrp"
    keys = [default, *sorted(k for k in POPULATIONS if k != default)]
    for key in keys:
        population = POPULATIONS[key]
        for sex in Sex:
            reference = population.reference[sex]
            cited.append(population.source)
            cited += [v.source for v in reference.values() if v.source]
            a = population.anthropometry.get(sex)
            if a is None:
                spread = "none"
            else:
                cited.append(a.source)
                lo, hi = a.age_range
                adult = a.at(min(max(AGE_REFERENCE, lo), hi))
                corr = (
                    "independent"
                    if adult.rho is None
                    else f"correlation {adult.rho:.2f} at 35 y"
                )
                spread = (
                    f"{_number(lo)}-{_number(hi)} y, {len(a.ages)} ages; "
                    f"HEIGHT &plusmn; {_number(adult.height_sd)} cm, "
                    f"ln BW &plusmn; {_number(adult.ln_bw_sd)}, {corr}"
                )
            rows.append(
                [
                    f"{population.name} `{key}`",
                    sex.value,
                    _number(reference["BW"].value),
                    _number(reference["HEIGHT"].value),
                    _number(reference["AGE"].value),
                    spread,
                ]
            )
    reference = POPULATIONS[default].reference[Sex.MALE]
    header = [
        "population",
        "sex",
        *(f"{c} [{_unit_id(reference[c].unit)}]" for c in ("BW", "HEIGHT", "AGE")),
        "distribution over age",
    ]
    lines = _table(header, rows)
    return _generated("`database/tables/populations/*.toml`", lines, cited)


AGES: tuple[float, ...] = (
    0.0,
    1.0,
    5.0,
    10.0,
    15.0,
    20.0,
    35.0,
    50.0,
    65.0,
    80.0,
    100.0,
)
"""The ages (years) of the tables of the physiology over age."""


def _age(a: float) -> str:
    """An age in a table, the newborn at 0."""
    return "newborn" if a == 0 else _number(a)


def _physiology_age(sex: Sex) -> str:
    """The reference individual of a sex over age."""
    curves = reference_curves(sex)
    cited = ["ICRP89", *(s for c in curves.values() for s in c.sources)]
    cited += list(factor_curves(sex)[ids.HCT].sources)
    cited += list(factor_curves(sex)[ids.GFR].sources)
    rows = []
    organs = (Organ.BRAIN, Organ.LIVER, Organ.KIDNEY, Organ.MUSCLE)
    for a in AGES:
        ref = reference_individual(sex, a)
        bsa = 0.007184 * ref.bw**0.425 * ref.height**0.725
        masses = [curves[ids.ref(ids.organ_mass(o))].curve.value(a) for o in organs]
        rows.append(
            [
                _age(a),
                _number(ref.bw),
                _number(ref.height),
                _number(bsa),
                *(_number(m) for m in masses),
                _number(curves[ids.ref(ids.VBLOOD)].curve.value(a)),
                _number(curves[ids.ref(ids.CO)].curve.value(a)),
                _number(ref.gfr * 1e3),
                _number(ref.gfr * 1e3 * 1.73 / bsa),
                _number(reference_hematocrit(sex, a)),
            ]
        )
    header = [
        "age [yr]",
        "BW [kg]",
        "HEIGHT [cm]",
        "BSA [m2]",
        *(f"M{o.value} [kg]" for o in organs),
        "Vblood [l]",
        "CO [l/min]",
        "GFR [ml/min]",
        "GFR [ml/min/1.73 m2]",
        "HCT [-]",
    ]
    return _generated(
        "`database/tables/children_human.toml`, `physiology_human.toml` and "
        "`age_human.toml`",
        _table(header, rows),
        cited,
    )


def physiology_age_male() -> str:
    """The male reference individual over age."""
    return _physiology_age(Sex.MALE)


def physiology_age_female() -> str:
    """The female reference individual over age."""
    return _physiology_age(Sex.FEMALE)


def _age_factors(sex: Sex) -> str:
    """The factors of the age of a sex."""
    factors = factor_curves(sex)
    refs = _references(sex, load("physiology_human"))
    fractions = {
        o: refs[ids.flow_fraction(o)].value for o in Organ if o is not Organ.LUNG
    }
    flows = flow_fraction_curves(sex, fractions)
    keys = (ids.GFR, ids.HCT, ids.FFM, "albumin", "agp")
    cited = [s for k in keys for s in factors[k].sources]
    cited += list(flow_factors(sex).sources)
    rows = []
    for a in AGES:
        rows.append(
            [
                _age(a),
                *(_number(factors[k].curve.value(a)) for k in keys),
                _number(flows[Organ.BRAIN].value(a)),
                _number(flows[Organ.LIVER].value(a) / fractions[Organ.LIVER]),
                _number(flows[Organ.KIDNEY].value(a) / fractions[Organ.KIDNEY]),
                _number(flows[Organ.MUSCLE].value(a) / fractions[Organ.MUSCLE]),
            ]
        )
    header = [
        "age [yr]",
        "fage_GFR",
        "fage_HCT",
        "fage_FFM",
        "fage_albumin",
        "fage_agp",
        "brain share of CO",
        "liver group",
        "kidney",
        "other organs",
    ]
    return _generated("`database/tables/age_human.toml`", _table(header, rows), cited)


def age_factors_male() -> str:
    """The factors of the age of men."""
    return _age_factors(Sex.MALE)


def age_factors_female() -> str:
    """The factors of the age of women."""
    return _age_factors(Sex.FEMALE)


ONTOGENY_AGES: tuple[float, ...] = (0.0, 1 / 12, 0.5, 1.0, 2.0, 5.0, 18.0)
"""The ages (years) of the table of the ontogeny."""


def ontogeny() -> str:
    """The ontogeny of the proteins."""
    rows = []
    cited: list[str] = []
    for organ, protein in sorted(protein_factors(), key=lambda x: (x[0].value, x[1])):
        curve = protein_factor(organ, protein)
        assert curve is not None
        cited += list(curve.sources)
        rows.append(
            [
                f"`{organ.value}` {_ORGAN_NAMES[organ]}",
                protein,
                *[_number(curve.curve.value(a)) for a in ONTOGENY_AGES],
                ", ".join(f"`{k}`" for k in curve.sources),
            ]
        )
    header = [
        "organ",
        "protein",
        "newborn",
        "1 month",
        "6 months",
        "1 y",
        "2 y",
        "5 y",
        "18 y",
        "sources",
    ]
    return _generated("`database/tables/ontogeny.toml`", _table(header, rows), cited)


SEX_AGES: tuple[float, ...] = (5.0, 12.0, 15.0, 18.0, 35.0, 50.0, 60.0, 80.0)
"""The ages (years) of the table of the factors of the sex."""

_SEX_NAMES: dict[str, str] = {
    "CO": "cardiac output",
    "Qbr": "brain fraction of the cardiac output",
    "albumin": "plasma albumin",
    "t_stomach_fed": "transit time of the stomach, fed",
    "t_colon": "transit times of the large intestine",
}
"""The quantities of the factors of the sex in the table."""


def sex_factors() -> str:
    """The factors of the female build relative to the male, over age."""
    rows = []
    cited: list[str] = []
    entries = [
        (_SEX_NAMES[key], sex_factor(key), c)
        for key, c in sex_curves(Sex.FEMALE).items()
    ]
    for organ, protein in sorted(
        protein_sex_factors(), key=lambda x: (x[0].value, x[1])
    ):
        c = protein_sex_factor(Sex.FEMALE, organ, protein)
        assert c is not None
        entries.append(
            (
                f"{protein} {_ORGAN_NAMES[organ]}",
                sex_factor(f"{organ.value}_{protein.lower()}"),
                c,
            )
        )
    for name, pid, c in entries:
        cited += list(c.sources)
        rows.append(
            [
                name,
                f"`{pid}`",
                *[_number(c.curve.value(a)) for a in SEX_AGES],
                ", ".join(f"`{k}`" for k in c.sources),
            ]
        )
    header = ["quantity", "factor", *[f"{a:g} y" for a in SEX_AGES], "sources"]
    return _generated("`database/tables/sex_human.toml`", _table(header, rows), cited)


def _class_table(
    section: str,
    classes: list[str],
    exclude: Iterable[str] = (),
    first: str = "factor",
) -> tuple[list[str], list[str]]:
    """Format a table of the factors of the classes of a section of the impairment.

    A factor states its values or is the consensus of the studies on a
    quantity (`class_factors`); a model choice is marked with `*`.
    """
    table: Mapping[str, Any] = load("impairment")
    for part in section.split("."):
        table = table[part]
    rows: list[list[str]] = []
    cited: list[str] = []
    for key, entry in table.items():
        if not isinstance(entry, dict) or key in exclude:
            continue
        if "value" not in entry and "consensus" not in entry:
            continue
        factors = class_factors(section, key)
        cells = [
            _number(v) + ("*" if c in factors.choice else "")
            for v, c in zip(factors.values, classes, strict=True)
        ]
        cited += factors.sources
        source = (
            f"consensus of `{entry['consensus']}`: " if "consensus" in entry else ""
        ) + ", ".join(f"`{s}`" for s in factors.sources)
        rows.append([f"`{key}`", *cells, source])
    return _table([first, *classes, "source"], rows), cited


_CHOICE = "A value marked `*` is a choice of the library, not a value of the source."


_CONSENSUS = (
    "A factor with a consensus is that of the studies of the database on the "
    "quantity, the ratio of the class to the controls of each study pooled over "
    "the studies, see [the database](../database.md); a value marked `*` is a "
    "choice of the library (for a consensus the interpolation of a class without "
    "a study)."
)
_SOURCES = "`database/tables/impairment.toml` and the studies of `pkpdsim/database`"


def impairment_child_pugh() -> str:
    """The factors of the Child-Pugh classes on the physiology."""
    lines, cited = _class_table("child_pugh", ["A", "B", "C"])
    return _generated(_SOURCES, [*lines, "", _CONSENSUS], cited)


def impairment_child_pugh_comparison() -> str:
    """The factors of the library against the Simcyp cirrhosis populations."""
    table = load("comparisons")["johnson2010"]["child_pugh"]
    cp = load("impairment")["child_pugh"]
    volume = class_factors("child_pugh", "functional_liver").values
    johnson_volume = [
        v / table["liver_volume"]["value"][0]
        for v in table["liver_volume"]["value"][1:]
    ]
    rows = []
    for entry in table.values():
        control, *classes = entry["value"]
        ratios = [v / control for v in classes]
        library = entry["library"]
        if library in cp:
            factors = list(class_factors("child_pugh", library).values)
            what = "physiology"
        else:
            try:
                factors = list(child_pugh_protein(library).values)
            except ValueError:
                continue
            # the capacity of the liver: abundance times the liver tissue
            factors = [f * v for f, v in zip(factors, volume, strict=True)]
            ratios = [r * v for r, v in zip(ratios, johnson_volume, strict=True)]
            what = "capacity"
        rows.append(
            [
                f"`{library}`",
                what,
                " / ".join(_number(round(f, 3)) for f in factors),
                " / ".join(_number(round(r, 3)) for r in ratios),
            ]
        )
    lines = _table(
        ["factor", "kind", "library A / B / C", "Johnson 2010 A / B / C"], rows
    )
    note = (
        "The physiological factors relative to healthy; for an enzyme the "
        "capacity of the liver relative to healthy, the abundance per amount of "
        "tissue times the liver volume in the library (`functional_liver`), the "
        "abundance per mg of microsomal protein times the functional hepatocyte "
        "volume in Johnson 2010 (Table III, the parameters of the Simcyp "
        "cirrhosis populations)."
    )
    return _generated(
        "`database/tables/comparisons.toml` and the factors of the library",
        [*lines, "", note],
        ["Johnson2010"],
    )


def impairment_meld() -> str:
    """The MELD scores of the severities of cirrhosis, the anchors of the mapping."""
    rows = []
    cited: list[str] = []
    for anchor in meld_anchors():
        if anchor.level is None:
            rows.append(
                [_number(anchor.severity), "healthy", _number(anchor.meld), "", ""]
            )
            continue
        lo, hi = anchor.level.ci
        cited += anchor.level.sources
        rows.append(
            [
                _number(anchor.severity),
                "ABC"[int(anchor.severity) - 1],
                _number(anchor.meld),
                f"{_number(lo)}-{_number(hi)}",
                ", ".join(
                    f"`{e.study}` {_number(e.value)}" for e in anchor.level.estimates
                ),
            ]
        )
    lines = _table(
        ["`CP_severity`", "class", "MELD", "95 % CI", "studies (MELD)"], rows
    )
    return _generated("the studies of `pkpdsim/database`", lines, cited)


def _consensus_cell(c: integration.Consensus | None) -> str:
    """A cell of a consensus: the ratio, its 95 % CI and the number of studies."""
    if c is None:
        return ""
    lo, hi = c.ci
    mark = " (mRNA)" if c.kind == "mrna" else ""
    return (
        f"{_number(c.value)} [{_number(lo)}-{_number(hi)}] ({len(c.estimates)}){mark}"
    )


def impairment_child_pugh_proteins() -> str:
    """The factors of the Child-Pugh classes on the proteins of the library in the liver."""
    rows = []
    cited: list[str] = []
    for protein in sorted(child_pugh_proteins() & set(PROTEINS)):
        entry = child_pugh_protein(protein)
        cells = [
            _number(v) + ("*" if c in entry.choice else "")
            for v, c in zip(entry.values, ("A", "B", "C"), strict=True)
        ]
        cited += entry.sources
        rows.append(
            [f"`{protein}`", *cells, ", ".join(f"`{s}`" for s in entry.sources)]
        )
    lines = _table(["protein", "A", "B", "C", "studies"], rows)
    note = (
        "The consensus of the studies of the database on the abundance per amount "
        "of liver tissue over the controls, see [the database](../database.md); "
        "a value marked `*` is the interpolation of a class without a study, a "
        "choice of the library."
    )
    return _generated("the studies of `pkpdsim/database`", [*lines, "", note], cited)


def database_studies() -> str:
    """The studies of the database with their measurements and groups."""
    rows = []
    for key, study in studies().items():
        for m in study.measurements.values():
            groups = sorted({r.group.id for r in m.records} - {study.control})
            rows.append(
                [
                    f"`{key}`",
                    KINDS[m.kind],
                    m.unit,
                    str(len({r.entity for r in m.records})),
                    str(len(m.records)),
                    str(len(groups) + 1),
                    m.location,
                ]
            )
    lines = _table(
        ["study", "measurement", "unit", "entities", "records", "groups", "location"],
        rows,
    )
    return _generated("`pkpdsim/database/studies`", lines, list(studies()))


def database_cirrhosis() -> str:
    """The consensus of every entity of the liver in the Child-Pugh classes."""
    rows = []
    cited: list[str] = []
    for entity in sorted(
        set(integration.entities("protein")) | set(integration.entities("mrna"))
    ):
        cells = []
        for cls in ("A", "B", "C"):
            c = integration.factor(
                entity, {"condition": "cirrhosis", "child_pugh": cls}
            )
            cells.append(_consensus_cell(c))
            if c is not None:
                cited += c.sources
        if any(cells):
            rows.append([f"`{entity}`", *cells])
    lines = _table(["entity", "A", "B", "C"], rows)
    note = (
        "The ratio of the class to the controls with its 95 % confidence "
        "interval and the number of studies, the protein abundance per amount "
        "of liver tissue, `(mRNA)` the mRNA expression of an entity without a "
        "protein abundance."
    )
    return _generated(
        "the studies of `pkpdsim/database`",
        [*lines, "", note],
        cited,
    )


def _level_cell(level: integration.Level | None) -> str:
    """A cell of a level: the value, its 95 % CI and the number of studies."""
    if level is None:
        return ""
    lo, hi = level.ci
    return (
        f"{_number(level.value)} [{_number(lo)}-{_number(hi)}] ({len(level.estimates)})"
    )


def _physiology_table(condition: str, field: str | None, classes: list[str]) -> str:
    """The levels and the ratios of the quantities of the subjects in a condition.

    Args:
        condition: the condition of the groups, e.g. `cirrhosis`.
        field: the field of the class of the groups, `None` for the condition
            as a whole.
        classes: the classes, the columns of the table.
    """
    rows = []
    cited: list[str] = []
    for quantity, q in QUANTITIES.items():
        controls = integration.level(quantity, {"condition": "control"})
        cells = []
        for cls in classes:
            fields = {"condition": condition}
            if field is not None:
                fields[field] = cls
            lv = integration.level(quantity, fields)
            c = integration.consensus(
                quantity, fields, kind="physiology", tissue="body"
            )
            ratio = _consensus_cell(c)
            cells.append(f"{_level_cell(lv)}<br>{ratio}" if ratio else _level_cell(lv))
            cited += [*(lv.sources if lv else ()), *(c.sources if c else ())]
        if not any(cells):
            continue
        if controls is not None:
            cited += controls.sources
        rows.append([f"`{quantity}` [{q.unit}]", _level_cell(controls), *cells])
    lines = _table(["quantity", "controls", *classes], rows)
    note = (
        "Per cell the level, the mean or median of the groups pooled over the "
        "studies on the linear scale, and below it the ratio to the controls of "
        "the studies with a control group, pooled on the log scale; each with "
        "its 95 % confidence interval and the number of studies. The controls "
        "are the control groups of every study of the quantity."
    )
    return _generated("the studies of `pkpdsim/database`", [*lines, "", note], cited)


def database_physiology() -> str:
    """The levels and ratios of the quantities of the subjects in cirrhosis."""
    return _physiology_table("cirrhosis", "child_pugh", ["A", "B", "C"])


def database_physiology_ckd() -> str:
    """The levels and ratios of the quantities of the subjects by the CKD stage."""
    return _physiology_table(
        "ckd", "ckd_stage", ["mild", "moderate", "severe", "end stage", "dialysis"]
    )


def database_physiology_heart_failure() -> str:
    """The levels and ratios of the quantities of the subjects by the NYHA class."""
    by_class = _physiology_table("heart_failure", "nyha", ["I", "II", "III", "IV"])
    overall = _physiology_table("heart_failure", None, ["heart failure"])
    return by_class + "\n\n" + overall.split("\n", 2)[2]


def impairment_ckd() -> str:
    """The factors of the CKD stages on the physiology."""
    section = dict(load("impairment")["ckd"])
    stages = section.pop("stage_gfr")
    anchors = section.pop("stage_anchor")
    tubular = section.pop("tubular_capacity")
    lines, cited = _class_table(
        "ckd",
        ["mild", "moderate", "severe", "end stage"],
        exclude=("stage_gfr", "stage_anchor", "tubular_capacity"),
    )
    limits = ", ".join(_number(v) for v in stages["value"])
    points = ", ".join(_number(v) for v in anchors["value"])
    head = [
        f"Stage limits of the GFR per 1.73 m2 ({stages['unit']}): {limits} "
        f"(`{stages['source']}`).",
        "",
        f"Anchors of the interpolation, the eGFR ({anchors['unit']}) of the "
        f"values of mild, moderate, severe and end stage: {points} "
        f"(`{anchors['source']}`).",
        "",
    ]
    tail = [
        f"Share of the tubular capacity which follows the renal function "
        f"(`f_tubular`, not per stage): {_number(tubular['value'])} "
        f"(`{tubular['source']}`).",
    ]
    return _generated(
        _SOURCES,
        [*head, *lines, "", _CONSENSUS, "", *tail],
        [*cited, stages["source"], anchors["source"], tubular["source"]],
    )


def impairment_ckd_proteins() -> str:
    """The factors of the CKD stages on the proteins of the liver (`CKD_hepatic`)."""
    stages = ["mild", "moderate", "severe", "end stage"]
    lines, cited = _class_table("ckd.proteins", stages, first="protein")
    return _generated("`database/tables/impairment.toml`", [*lines, "", _CHOICE], cited)


def impairment_ckd_binding() -> str:
    """The CKD profile of the binding per protein."""
    section = load("impairment")["ckd"]["binding"]
    lines, cited = _class_table(
        "ckd.binding", ["mild", "moderate", "severe", "end stage"]
    )
    rows = "; ".join(
        f"`{key}` the {entry['quantity']} (protein `{entry['protein']}`)"
        for key, entry in section.items()
    )
    tail = [f"The rows: {rows}. The factor of a protein is the product of its rows."]
    return _generated(_SOURCES, [*lines, "", _CONSENSUS, "", *tail], cited)


def impairment_heart_failure() -> str:
    """The factors of the NYHA classes on the blood flows."""
    lines, cited = _class_table("heart_failure", ["I", "II", "III", "IV"])
    return _generated(_SOURCES, [*lines, "", _CONSENSUS], cited)


def arm_vein() -> str:
    """The presets of the arm vein sampling site."""
    data = load("arm_vein")
    sites = sorted({site for preset in data.values() for site in preset})
    rows: list[list[str]] = []
    cited: list[str] = []
    for name, preset in data.items():
        preset_sources = sorted({e["source"] for e in preset.values()})
        cited += preset_sources
        rows.append(
            [
                f"`{name}`",
                *(_value(preset.get(site)) for site in sites),
                ", ".join(f"`{s}`" for s in preset_sources),
            ]
        )
    lines = _table(["preset", *sites, "source"], rows)
    return _generated("`database/tables/arm_vein.toml`", lines, cited)


def intestine_segments() -> str:
    """The lumen segments of the segmented intestine."""
    data = load("intestine_human")
    segments = {
        name: entry
        for name, entry in data.items()
        if "transit_time" in entry and name != "single_lumen"
    }
    columns = {
        "volume": "volume",
        "transit_time": "transit",
        "transit_time_fed": "transit fed",
        "ph": "pH",
        "area_factor": "area factor",
        "cyp3a4": "CYP3A4",
    }
    rows = [
        [name.replace("_", " "), *(_value(entry.get(k)) for k in columns)]
        for name, entry in segments.items()
    ]
    entries = list(segments.values())
    header = ["segment", *(_header(label, entries, k) for k, label in columns.items())]
    line, cited = _column_sources(segments.values(), columns)
    lines = [*_table(header, rows), "", line]
    return _generated("`database/tables/intestine_human.toml`", lines, cited)


def nephron_segments() -> str:
    """The tubule segments of the nephron kidney."""
    data = load("kidney_human")
    rows = [
        [
            name.replace("_", " ").replace("henle", "Henle"),
            _value(entry["flow_fraction"]),
            _value(entry["volume"]),
        ]
        for name, entry in data.items()
    ]
    entries = list(data.values())
    header = [
        "segment",
        "flow leaving, fraction of GFR " + f"[{_unit(entries, 'flow_fraction')}]",
        _header("volume", entries, "volume"),
    ]
    columns = {"flow_fraction": "flow", "volume": "volume"}
    line, cited = _column_sources(data.values(), columns)
    lines = [*_table(header, rows), "", line]
    return _generated("`database/tables/kidney_human.toml`", lines, cited)


def proteins() -> str:
    """The protein table: the enzymes and transporters with their annotation."""
    data = load("proteins")
    rows = [
        [
            f"`{key}`",
            entry["name"],
            f"[{entry['uniprot']}](https://www.uniprot.org/uniprotkb/{entry['uniprot']})",
            f"`{entry['ncit']}`" if "ncit" in entry else "",
        ]
        for key, entry in data["proteins"].items()
    ]
    turnover = data["turnover"]
    lines = _table(["protein", "name", "UniProt", "NCIT"], rows)
    lines += [
        "",
        "Degradation rate `kdeg_<protein>` of the turnover state "
        + f"({_unit([turnover], next(iter(turnover)))}): "
        + ", ".join(f"{organ} {_value(e)}" for organ, e in turnover.items())
        + ".",
    ]
    cited = [data["source"], *(e["source"] for e in turnover.values())]
    return _generated("`database/tables/proteins.toml`", lines, cited)


# --- the clinical classes ---------------------------------------------------

_CP_SCORE = [
    [
        "total bilirubin",
        "< 2 mg/dL (< 34 umol/L)",
        "2-3 mg/dL (34-50 umol/L)",
        "> 3 mg/dL (> 50 umol/L)",
    ],
    [
        "serum albumin",
        "> 3.5 g/dL (> 35 g/L)",
        "2.8-3.5 g/dL (28-35 g/L)",
        "< 2.8 g/dL (< 28 g/L)",
    ],
    ["prothrombin time: INR", "< 1.7", "1.7-2.3", "> 2.3"],
    ["(or its prolongation)", "< 4 s", "4-6 s", "> 6 s"],
    [
        "ascites",
        "none",
        "mild (controlled by diuretics)",
        "moderate to severe (refractory)",
    ],
    [
        "hepatic encephalopathy",
        "none",
        "grade I-II (or controlled by medication)",
        "grade III-IV (refractory)",
    ],
]
"""The criteria of the Child-Pugh score and their points 1, 2, 3."""


def _limits(c: ImpairmentClass, *, closed: bool) -> str:
    """The range of a class, e.g. `45-59` or `>= 90`.

    A score has closed integer limits (`5-6`), an eGFR the range `[lo, hi)`,
    printed as the integers in it (`45-59`).
    """
    if c.limits is None:
        return ""
    lo, hi = c.limits
    if lo is None:
        return "" if hi is None else f"< {_number(hi)}"
    if hi is None:
        return f">= {_number(lo)}"
    return f"{_number(lo)}-{_number(hi if closed else hi - 1)}"


def _swatch(color: str) -> str:
    """A square of the color of a class in a table cell, its code on hover."""
    return (
        f'<span class="pbpk-swatch" style="background:{color}" title="{color}"></span>'
    )


def _class_cell(c: ImpairmentClass) -> str:
    """The color, the key and the label of a class in a table cell."""
    return f"{_swatch(c.color)} `{c.key}` {c.label}"


def classes_child_pugh() -> str:
    """The Child-Pugh score and its classes with their severity."""
    score = _table(["criterion", "1 point", "2 points", "3 points"], _CP_SCORE)
    rows = [
        [
            _class_cell(c),
            c.definition,
            _limits(c, closed=True),
            _number(c.value) if c.value is not None else "",
        ]
        for c in CHILD_PUGH_CLASSES
    ]
    classes = _table(["class", "definition", "score", "CP_severity"], rows)
    return "\n".join(
        [
            "*Generated from `pkpdsim.impairment.CHILD_PUGH_CLASSES` by "
            "`scripts/pbpk_tables.py`.*",
            "",
            *score,
            "",
            *classes,
        ]
    )


def classes_age() -> str:
    """The age classes of ICH E11 and ICH E7 with the age of the scans."""
    rows = [
        [
            _class_cell(c),
            c.definition,
            _limits(c, closed=False) if c.key != "newborn" else "0-27 days",
            _number(c.value) if c.value is not None else "",
        ]
        for c in AGE_CLASSES
    ]
    return "\n".join(
        [
            "*Generated from `pkpdsim.impairment.AGE_CLASSES` by "
            "`scripts/pbpk_tables.py`.*",
            "",
            *_table(["class", "definition", "age [yr]", "AGE"], rows),
        ]
    )


def classes_bmi() -> str:
    """The WHO classes of the body mass index with the body weight of the scans."""
    rows = []
    for c in BMI_CLASSES:
        weights = ", ".join(
            f"{sex.value} "
            + _number(
                reference_individual(sex).bw
                if c.value is None
                else c.value * (reference_individual(sex).height / 100.0) ** 2
            )
            for sex in Sex
        )
        rows.append(
            [
                _class_cell(c),
                c.definition,
                _limits(c, closed=False),
                _number(c.value) if c.value is not None else "reference",
                weights,
            ]
        )
    return "\n".join(
        [
            "*Generated from `pkpdsim.impairment.BMI_CLASSES` by "
            "`scripts/pbpk_tables.py`.*",
            "",
            *_table(
                ["class", "definition", "BMI [kg/m2]", "BMI", "BW [kg]"],
                rows,
            ),
        ]
    )


def classes_kdigo() -> str:
    """The GFR categories of KDIGO with the eGFR and renal function of the scans."""
    stages = {
        "G2": "mild",
        "G3a": "moderate",
        "G3b": "moderate",
        "G4": "severe",
        "G5": "end stage",
    }
    rows = []
    for c in KDIGO_GFR_CLASSES:
        fs = []
        for sex in Sex:
            changes = PBPKSimulationExperiment.class_changes(
                KDIGO_GFR_CLASSES, c.key, sex=sex
            )
            f = changes["f_renal_function"].m_as("dimensionless")
            fs.append(_number(float(f"{f:.4g}")))
        rows.append(
            [
                _class_cell(c),
                c.definition,
                stages.get(c.key, "normal"),
                _number(c.value) if c.value is not None else "reference",
                " / ".join(fs),
            ]
        )
    lines = _table(
        ["class", "definition", "FDA stage", "eGFR of the scan", "f_renal_function"],
        rows,
    )
    return "\n".join(
        [
            "*Generated from `pkpdsim.impairment.KDIGO_GFR_CLASSES` by "
            "`scripts/pbpk_tables.py`; the eGFR in ml/min/1.73 m2, "
            "`f_renal_function` of the male / female build.*",
            "",
            *lines,
        ]
    )


def classes_nyha() -> str:
    """The NYHA classes with their severity."""
    rows = [
        [
            _class_cell(c),
            c.definition,
            _number(c.value) if c.value is not None else "",
        ]
        for c in NYHA_CLASSES
    ]
    lines = _table(["class", "definition", "HF_severity"], rows)
    return "\n".join(
        [
            "*Generated from `pkpdsim.impairment.NYHA_CLASSES` by "
            "`scripts/pbpk_tables.py`.*",
            "",
            *lines,
        ]
    )


_FONT = 'font-family="system-ui, sans-serif"'


def _svg(width: float, height: float, title: str, body: list[str]) -> str:
    """An inline SVG of the theme: text and axes in `currentColor`."""
    return "\n".join(
        [
            '<div class="pbpk-figure pbpk-wide" markdown="0">',
            f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {_number(width)} '
            f'{_number(height)}" role="img" aria-label="{title}" {_FONT} '
            'font-size="14" fill="currentColor">',
            f"<title>{title}</title>",
            *body,
            "</svg>",
            "</div>",
        ]
    )


def _band_color(c: ImpairmentClass) -> tuple[str, str]:
    """The fill of a class and its opacity; the healthy black is a neutral gray."""
    return ("#9e9e9e", "0.35") if c.color == "black" else (c.color, "0.9")


_PICTOGRAM = 72.0
"""The height of the row of the pictograms of the classes figure."""


def _scale_row(
    y: float,
    name: str,
    bands: list[tuple[ImpairmentClass, float, float]],
    marks: list[tuple[float, str]],
    note: str,
    organ: str,
    losses: list[tuple[float, str]],
) -> list[str]:
    """A row of the classes figure: the pictograms, the bands of the classes and the marks.

    Args:
        y: the top of the pictograms.
        name: the title of the row.
        bands: the classes with their extent on the axis.
        marks: the values of the scan, `(x, label)`.
        note: the text under the row.
        organ: the organ of the pictograms (`li`, `ki`, `he`).
        losses: the lost fraction of the organ of each class and its hover text.
    """
    out = [f'<text x="10" y="{_number(y - 10)}" font-weight="600">{name}</text>']
    width = min(x1 - x0 for _, x0, x1 in bands) - 10
    for (c, x0, x1), (loss, text) in zip(bands, losses, strict=True):
        fill = "#9e9e9e" if c.color == "black" else c.color
        out.append(
            icons.impaired_organ(
                organ,
                (x0 + x1) / 2,
                y + (_PICTOGRAM - 8) / 2,
                width,
                _PICTOGRAM - 14,
                loss,
                fill,
                title=f"{c.label}: {text}",
            )
        )
    y += _PICTOGRAM
    for c, x0, x1 in bands:
        fill, opacity = _band_color(c)
        out.append(
            f'<rect x="{_number(x0 + 1)}" y="{_number(y)}" width="{_number(x1 - x0 - 2)}" '
            f'height="22" rx="3" fill="{fill}" fill-opacity="{opacity}">'
            f"<title>{html.escape(f'{c.label}: {c.definition}')}</title></rect>"
        )
        out.append(
            f'<text x="{_number((x0 + x1) / 2)}" y="{_number(y + 41)}" '
            f'text-anchor="middle">{c.key}</text>'
        )
    for x, label in marks:
        out.append(
            f'<line x1="{_number(x)}" x2="{_number(x)}" y1="{_number(y - 3)}" '
            f'y2="{_number(y + 25)}" stroke="currentColor" stroke-width="2"/>'
        )
        out.append(
            f'<text x="{_number(x)}" y="{_number(y + 58)}" text-anchor="middle" '
            f'font-size="12" opacity="0.8">{label}</text>'
        )
    out.append(
        f'<text x="10" y="{_number(y + 80)}" font-size="12" opacity="0.8">{note}</text>'
    )
    return out


_EGFR_NORMAL = 114.5
"""The eGFR of the male build (ml/min/1.73 m2), the renal function 1."""


def _hf_cardiac_output(severity: float) -> float:
    """The cardiac output of heart failure relative to healthy, from its flows.

    The factors of heart failure multiply the flows of their organs, the
    cardiac output is the sum of the flows: `1 + sum f_o (factor_o - 1)` with the
    flow fractions `f_o` of the reference man.
    """
    refs = _references(Sex.MALE, load("physiology_human"))
    co = 1.0
    for key, organs in _HF_ORGANS.items():
        share = sum(refs[ids.flow_fraction(o)].value for o in organs)
        co += share * (_hf_factor(key)(severity) - 1.0)
    return co


def classes_figure() -> str:
    """The classes of the three conditions on their scales, with the scan values.

    Above every class a pictogram of its organ, the lost part pale: the liver by
    the portosystemic shunt of the class, the kidney by the eGFR of the class
    over that of the build, the heart by the cardiac output of the class.
    """
    left, right = 10.0, 710.0
    healthy = 90.0
    row = 118.0 + _PICTOGRAM
    body: list[str] = []
    # Child-Pugh: the score 5 to 15, healthy before it
    span = (right - left - healthy) / 11
    cp = list(CHILD_PUGH_CLASSES)
    bands = [(cp[0], left, left + healthy)]
    for c in cp[1:]:
        assert c.limits is not None
        lo, hi = c.limits
        assert lo is not None and hi is not None
        bands.append(
            (c, left + healthy + (lo - 5) * span, left + healthy + (hi - 4) * span)
        )
    marks = [((x0 + x1) / 2, _number(c.value or 0)) for c, x0, x1 in bands]
    melds = ", ".join(
        f"{'ABC'[int(a.severity) - 1]} {a.meld:.3g}" for a in meld_anchors()[1:]
    )
    shunt = _cp_factor("shunt")
    losses = [
        (
            shunt(c.value or 0.0),
            f"shunt {shunt(c.value or 0.0):.2f}, the portal blood which bypasses the liver",
        )
        for c in cp
    ]
    body += _scale_row(
        36,
        "Child-Pugh classes by the score 5-15, the CP_severity of each below",
        bands,
        marks,
        f"MELD of the classes (consensus of the database): healthy 6.43, {melds}",
        "li",
        losses,
    )

    # KDIGO: the eGFR from 120 down to 0
    def gx(egfr: float) -> float:
        return left + (120 - min(egfr, 120)) / 120 * (right - left)

    bands = []
    losses = []
    for c in KDIGO_GFR_CLASSES:
        assert c.limits is not None
        lo, hi = c.limits
        bands.append((c, gx(hi if hi is not None else 120), gx(lo or 0)))
        egfr = min(c.value or _EGFR_NORMAL, _EGFR_NORMAL)
        losses.append(
            (
                1.0 - egfr / _EGFR_NORMAL,
                f"eGFR {_number(egfr)}, {100 * egfr / _EGFR_NORMAL:.0f} % of the "
                "renal function of the build",
            )
        )
    marks = [
        (gx(c.value), f"{_number(c.value)}")
        for c in KDIGO_GFR_CLASSES
        if c.value is not None
    ]
    body += _scale_row(
        36 + row,
        "KDIGO GFR categories by the eGFR (120 to 0 ml/min/1.73 m2), the eGFR of each below",
        bands,
        marks,
        "the eGFR of the scan of a category; G1 is the renal function 1 of the "
        f"build (eGFR {_number(_EGFR_NORMAL)} male, 114.7 female)",
        "ki",
        losses,
    )
    # NYHA: equal bands
    ny = list(NYHA_CLASSES)
    width = (right - left) / len(ny)
    bands = [(c, left + i * width, left + (i + 1) * width) for i, c in enumerate(ny)]
    marks = [((x0 + x1) / 2, _number(c.value or 0)) for c, x0, x1 in bands]
    losses = []
    for c in ny:
        co = _hf_cardiac_output(c.value or 0.0)
        losses.append(
            (1.0 - co, f"cardiac output {co:.2f} of healthy (the sum of its flows)")
        )
    body += _scale_row(
        36 + 2 * row,
        "NYHA functional classes, the HF_severity of each below",
        bands,
        marks,
        "the classes grade the symptoms at physical activity",
        "he",
        losses,
    )
    return _svg(
        720,
        36 + 2 * row + _PICTOGRAM + 90,
        "The clinical classes of the impairments and their scan values",
        body,
    )


def _factor_panel(
    x0: float,
    y0: float,
    title: str,
    points: list[tuple[float, float]],
    marks: list[tuple[float, float, str]],
    bands: list[tuple[ImpairmentClass, float, float]],
    xscale: Callable[[float], float],
) -> list[str]:
    """A small panel of a factor over the severity: the class bands and the curve.

    Args:
        x0: the left of the panel.
        y0: the top of the panel.
        title: the name of the factor.
        points: the curve, `(severity, factor)` in the order of the x axis.
        marks: the classes on the curve, `(severity, factor, label)`.
        bands: the classes with their range of the severity, for the strip
            under the curve.
        xscale: the severity to the fraction of the width of the panel.
    """
    w, h = 136.0, 70.0
    values = [v for _, v in points]
    top = max(1.0, *values) * 1.08

    def px(x: float) -> float:
        return x0 + xscale(x) * w

    def py(v: float) -> float:
        return y0 + 18 + h - v / top * h

    out = [
        f'<text x="{_number(x0)}" y="{_number(y0 + 11)}" font-size="13">{title}</text>'
    ]
    for c, a, b in bands:
        fill, opacity = _band_color(c)
        out.append(
            f'<rect x="{_number(px(a))}" y="{_number(y0 + 18 + h + 3)}" '
            f'width="{_number(px(b) - px(a))}" height="5" fill="{fill}" '
            f'fill-opacity="{opacity}"><title>{c.label}</title></rect>'
        )
    out.append(
        f'<line x1="{_number(x0)}" x2="{_number(x0 + w)}" y1="{_number(py(1.0))}" '
        f'y2="{_number(py(1.0))}" stroke="currentColor" stroke-opacity="0.35" '
        'stroke-dasharray="3 3"/>'
    )
    out.append(
        f'<text x="{_number(x0 + w + 3)}" y="{_number(py(1.0) + 3)}" font-size="11" '
        'opacity="0.7">1</text>'
    )
    line = " ".join(f"{_number(px(x))},{_number(py(v))}" for x, v in points)
    out.append(
        f'<polyline points="{line}" fill="none" stroke="currentColor" '
        'stroke-width="2" stroke-linejoin="round"/>'
    )
    for x, v, label in marks:
        out.append(
            f'<circle cx="{_number(px(x))}" cy="{_number(py(v))}" r="3.5" '
            f'fill="currentColor"><title>{label}: {_number(v)}</title></circle>'
        )
    hi, lo = max(values), min(values)
    for v, show in ((hi, hi > 1.05), (lo, lo < 0.95)):
        # no label which would overlap the label of 1
        if show and abs(py(v) - py(1.0)) > 11:
            out.append(
                f'<text x="{_number(x0 + w + 3)}" y="{_number(py(v) + 3)}" '
                f'font-size="11" opacity="0.7">{_number(v)}</text>'
            )
    return out


def _factors_figure(
    title: str,
    panels: list[tuple[str, list[tuple[float, float]], list[tuple[float, float, str]]]],
    bands: list[tuple[ImpairmentClass, float, float]],
    xscale: Callable[[float], float],
) -> str:
    """The factors of a condition over its severity, a panel per factor."""
    columns = 4
    rows = (len(panels) + columns - 1) // columns
    body: list[str] = [
        f'<text x="10" y="18" font-weight="600">{title}</text>',
    ]
    for i, (key, points, marks) in enumerate(panels):
        x0 = 10 + (i % columns) * 180
        y0 = 34 + (i // columns) * 110
        body += _factor_panel(x0, y0, key, points, marks, bands, xscale)
    return _svg(730, 34 + rows * 110, title, body)


def _class_panels(
    section: str, keys: Iterable[str], severities: list[float], labels: list[str]
) -> list[tuple[str, list[tuple[float, float]], list[tuple[float, float, str]]]]:
    """The panels of the factors of a section whose classes are at `severities`."""
    panels = []
    for key in keys:
        f = class_factors(section, key)
        healthy = 0.0 if key == "shunt" else 1.0
        marks = [
            (x, v, label)
            for x, v, label in zip(severities, f.values, labels, strict=True)
        ]
        points = [(0.0, healthy), *((x, v) for x, v, _ in marks)]
        name = key if "." not in section else f"{section.split('.')[1]}: {key}"
        panels.append((name, points, marks))
    return panels


def _entries(section: Mapping[str, Any], skip: Iterable[str] = ()) -> list[str]:
    """The keys of the factors of a section of the impairment."""
    return [
        k
        for k, v in section.items()
        if isinstance(v, dict) and k not in skip and ("value" in v or "consensus" in v)
    ]


def factors_child_pugh() -> str:
    """The physiological factors of the Child-Pugh classes over `CP_severity`."""
    keys = _entries(load("impairment")["child_pugh"])
    bands = [(c, i - 0.5, i + 0.5) for i, c in enumerate(CHILD_PUGH_CLASSES)]
    panels = _class_panels("child_pugh", keys, [1.0, 2.0, 3.0], ["A", "B", "C"])
    return _factors_figure(
        "Liver cirrhosis: the factors over CP_severity (healthy, A, B, C)",
        panels,
        bands,
        lambda x: (x + 0.5) / 4,
    )


def factors_ckd() -> str:
    """The physiological factors of the CKD stages over the eGFR."""
    table = load("impairment")["ckd"]
    stages = ["mild", "moderate", "severe", "end stage"]
    anchors = [float(v) for v in table["stage_anchor"]["value"]]
    normal = float(table["stage_gfr"]["value"][0])
    panels = []
    for section, keys in (
        ("ckd", _entries(table, ("stage_gfr", "stage_anchor", "tubular_capacity"))),
        ("ckd.binding", _entries(table["binding"])),
    ):
        for name, _, marks in _class_panels(section, keys, anchors, stages):
            # 1 from the normal function up, constant below the end stage
            points = [(120.0, 1.0), (normal, 1.0), *((x, v) for x, v, _ in marks)]
            points.append((0.0, marks[-1][1]))
            panels.append((name, points, marks))
    bands = []
    for c in KDIGO_GFR_CLASSES:
        assert c.limits is not None
        lo, hi = c.limits
        bands.append((c, hi if hi is not None else 120.0, lo or 0.0))
    return _factors_figure(
        "Chronic kidney disease: the factors over the eGFR (120 to 0 ml/min/1.73 m2)",
        panels,
        bands,
        lambda x: (120 - x) / 120,
    )


def factors_heart_failure() -> str:
    """The factors of the NYHA classes over `HF_severity`."""
    keys = _entries(load("impairment")["heart_failure"])
    bands = [(c, i - 0.5, i + 0.5) for i, c in enumerate(NYHA_CLASSES)]
    panels = _class_panels(
        "heart_failure", keys, [1.0, 2.0, 3.0, 4.0], ["I", "II", "III", "IV"]
    )
    return _factors_figure(
        "Heart failure: the factors over HF_severity (healthy, NYHA I to IV)",
        panels,
        bands,
        lambda x: (x + 0.5) / 5,
    )


def _piecewise(points: list[tuple[float, float]], x: float) -> float:
    """The linear interpolation of points `(x, y)` sorted by x, constant outside."""
    if x <= points[0][0]:
        return points[0][1]
    for (x0, y0), (x1, y1) in itertools.pairwise(points):
        if x <= x1:
            return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
    return points[-1][1]


def _cp_factor(key: str, section: str = "child_pugh") -> Callable[[float], float]:
    """A factor of cirrhosis as a function of `CP_severity`."""
    healthy = 0.0 if key == "shunt" else 1.0
    values = class_factors(section, key).values
    points = [(0.0, healthy), *((i + 1.0, v) for i, v in enumerate(values))]
    return lambda cp: _piecewise(points, cp)


def _ckd_factor(key: str, section: str = "ckd") -> Callable[[float], float]:
    """A factor of CKD as a function of the eGFR (ml/min/1.73 m2)."""
    table = load("impairment")["ckd"]
    anchors = [float(v) for v in table["stage_anchor"]["value"]]
    normal = float(table["stage_gfr"]["value"][0])
    values = class_factors(section, key).values
    points = sorted([(normal, 1.0), *zip(anchors, values, strict=True)])
    return lambda egfr: _piecewise(points, egfr)


_RAMP = ("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b")
"""The sequential blue of the heat maps, light (low) to dark (high)."""


def _ramp(t: float) -> str:
    """The color of a fraction 0 to 1 of the ramp."""
    t = min(max(t, 0.0), 1.0) * (len(_RAMP) - 1)
    i = min(int(t), len(_RAMP) - 2)
    a, b = _RAMP[i].lstrip("#"), _RAMP[i + 1].lstrip("#")
    f = t - i
    rgb = [
        round(
            int(a[k : k + 2], 16) + f * (int(b[k : k + 2], 16) - int(a[k : k + 2], 16))
        )
        for k in (0, 2, 4)
    ]
    return "#" + "".join(f"{v:02x}" for v in rgb)


def _heatmap(
    x0: float,
    title: str,
    xlabel: str,
    xs: list[float],
    xticks: list[tuple[float, str]],
    value: Callable[[float, float], float],
) -> list[str]:
    """A heat map of a quantity over a severity (x) and the eGFR (y, 120 to 0)."""
    w, h, top = 200.0, 200.0, 40.0
    egfrs = [120.0 - 120.0 * (j + 0.5) / 24 for j in range(24)]
    grid = [[value(x, e) for x in xs] for e in egfrs]
    hi = max(max(row) for row in grid)
    lo = min(min(row) for row in grid)
    cw, ch = w / len(xs), h / len(egfrs)
    out = [
        f'<text x="{_number(x0)}" y="18" font-weight="600">{title}</text>',
        f'<text transform="translate({_number(x0 - 28)} {_number(top + h / 2)}) '
        'rotate(-90)" text-anchor="middle" font-size="12">eGFR (ml/min/1.73 m2)'
        "</text>",
    ]
    for j, (e, row) in enumerate(zip(egfrs, grid, strict=True)):
        for i, (x, v) in enumerate(zip(xs, row, strict=True)):
            t = (v - lo) / (hi - lo) if hi > lo else 0.0
            out.append(
                f'<rect x="{_number(x0 + i * cw)}" y="{_number(top + j * ch)}" '
                f'width="{_number(cw + 0.3)}" height="{_number(ch + 0.3)}" '
                f'fill="{_ramp(t)}"><title>{xlabel} {_number(round(x, 3))}, eGFR '
                f"{_number(round(e, 1))}: {v:.3g}</title></rect>"
            )
    for egfr in (90, 60, 30, 15):
        y = top + (120 - egfr) / 120 * h
        out.append(
            f'<line x1="{_number(x0 - 4)}" x2="{_number(x0)}" y1="{_number(y)}" '
            f'y2="{_number(y)}" stroke="currentColor"/>'
        )
        out.append(
            f'<text x="{_number(x0 - 6)}" y="{_number(y + 4)}" text-anchor="end" '
            f'font-size="11">{egfr}</text>'
        )
    span = xs[-1] - xs[0]
    for x, label in xticks:
        px = x0 + (x - xs[0]) / span * (w - cw) + cw / 2
        out.append(
            f'<line x1="{_number(px)}" x2="{_number(px)}" y1="{_number(top + h)}" '
            f'y2="{_number(top + h + 4)}" stroke="currentColor"/>'
        )
        out.append(
            f'<text x="{_number(px)}" y="{_number(top + h + 17)}" text-anchor="middle" '
            f'font-size="11">{label}</text>'
        )
    out.append(
        f'<text x="{_number(x0 + w / 2)}" y="{_number(top + h + 33)}" '
        f'text-anchor="middle" font-size="12">{xlabel}</text>'
    )
    # the color bar
    for k in range(20):
        out.append(
            f'<rect x="{_number(x0 + k * w / 20)}" y="{_number(top + h + 42)}" '
            f'width="{_number(w / 20 + 0.3)}" height="8" fill="{_ramp(k / 19)}"/>'
        )
    out.append(
        f'<text x="{_number(x0)}" y="{_number(top + h + 64)}" font-size="11">'
        f"{lo:.2g}</text>"
    )
    out.append(
        f'<text x="{_number(x0 + w)}" y="{_number(top + h + 64)}" font-size="11" '
        f'text-anchor="end">{hi:.2g}</text>'
    )
    return out


def combinations_figure() -> str:
    """Quantities of combined impairments over the severities, heat maps."""
    normal = _EGFR_NORMAL
    cp_gfr = _cp_factor("gfr")
    cp_alb = _cp_factor("albumin")
    hf_kf = _hf_factor("kidney_flow")
    ckd_kf = _ckd_factor("kidney_flow")
    ckd_alb = _ckd_factor("albumin", "ckd.binding")
    ckd_bind = _ckd_factor("albumin_binding", "ckd.binding")
    fu_ref = 0.05

    def fu(cp: float, egfr: float) -> float:
        f = cp_alb(cp) * ckd_alb(egfr) * ckd_bind(egfr)
        return 1 / (1 + f * (1 - fu_ref) / fu_ref) / fu_ref

    cps = [3.0 * i / 23 for i in range(24)]
    hfs = [4.0 * i / 23 for i in range(24)]
    cp_ticks = [(0.0, "0"), (1.0, "A"), (2.0, "B"), (3.0, "C")]
    hf_ticks = [(0.0, "0"), (1.0, "I"), (2.0, "II"), (3.0, "III"), (4.0, "IV")]
    body: list[str] = []
    body += _heatmap(
        50,
        "GFR, hepatorenal",
        "CP_severity",
        cps,
        cp_ticks,
        lambda cp, e: cp_gfr(cp) * e / normal,
    )
    body += _heatmap(
        300,
        "renal blood flow, cardiorenal",
        "HF_severity",
        hfs,
        hf_ticks,
        lambda h, e: hf_kf(h) * ckd_kf(e),
    )
    body += _heatmap(
        550,
        "fu of an albumin-bound drug",
        "CP_severity",
        cps,
        cp_ticks,
        fu,
    )
    return _svg(
        760,
        310,
        "Combined impairments: GFR, renal blood flow and fu relative to healthy",
        body,
    )


def _hf_factor(key: str) -> Callable[[float], float]:
    """A factor of heart failure as a function of `HF_severity`."""
    values = class_factors("heart_failure", key).values
    points = [(0.0, 1.0), *((i + 1.0, v) for i, v in enumerate(values))]
    return lambda hf: _piecewise(points, hf)


TABLES: dict[str, Callable[[], str]] = {
    "physiology-general": physiology_general,
    "physiology-organs-male": physiology_organs_male,
    "species-overview": species_overview,
    "circadian-rhythms": circadian_rhythms,
    "adipose-factors": adipose_factors,
    "adipose-visceral": adipose_visceral,
    "adipose-waist": adipose_waist,
    "adipose-classes": adipose_classes,
    "pregnancy-factors": pregnancy,
    "species-intestine": species_intestine,
    "species-organs-rat": species_organs_rat,
    "species-organs-mouse": species_organs_mouse,
    "species-organs-pig": species_organs_pig,
    "physiology-organs-female": physiology_organs_female,
    "populations": populations,
    "physiology-age-male": physiology_age_male,
    "physiology-age-female": physiology_age_female,
    "age-factors-male": age_factors_male,
    "age-factors-female": age_factors_female,
    "ontogeny": ontogeny,
    "sex-factors": sex_factors,
    "classes-age": classes_age,
    "classes-bmi": classes_bmi,
    "impairment-child-pugh": impairment_child_pugh,
    "impairment-child-pugh-proteins": impairment_child_pugh_proteins,
    "impairment-meld": impairment_meld,
    "impairment-child-pugh-comparison": impairment_child_pugh_comparison,
    "database-physiology": database_physiology,
    "database-physiology-ckd": database_physiology_ckd,
    "database-physiology-heart-failure": database_physiology_heart_failure,
    "impairment-ckd": impairment_ckd,
    "impairment-ckd-binding": impairment_ckd_binding,
    "impairment-ckd-proteins": impairment_ckd_proteins,
    "impairment-heart-failure": impairment_heart_failure,
    "arm-vein": arm_vein,
    "intestine-segments": intestine_segments,
    "nephron-segments": nephron_segments,
    "proteins": proteins,
    "database-studies": database_studies,
    "database-cirrhosis": database_cirrhosis,
    "classes-figure": classes_figure,
    "classes-child-pugh": classes_child_pugh,
    "classes-kdigo": classes_kdigo,
    "classes-nyha": classes_nyha,
    "factors-child-pugh": factors_child_pugh,
    "factors-ckd": factors_ckd,
    "factors-heart-failure": factors_heart_failure,
    "combinations-figure": combinations_figure,
}
"""The tables by the name of their region."""


def render_page(text: str) -> str:
    """Write the tables into the regions of a page.

    Args:
        text: the markdown of the page.

    Returns:
        The page with every region rewritten.

    Raises:
        KeyError: if a region names no table of `TABLES`.
    """

    def replace(match: re.Match[str]) -> str:
        name = match["name"]
        if name == REFERENCES:
            return f"{begin_marker(name)}\n{references(text)}\n{END}"
        if name not in TABLES:
            raise KeyError(f"No table '{name}', the tables are {sorted(TABLES)}.")
        return f"{begin_marker(name)}\n{TABLES[name]()}\n{END}"

    return _REGION.sub(replace, text)


REFERENCES = "references"
"""The region of the references of a page, the sources its text cites."""

_CITATION = re.compile(r"\]\(#ref-(?P<key>[A-Za-z0-9_]+)\)")
"""A citation in the text of a page, a link `[Pugh 1973](#ref-Pugh1973)` to the
anchor of a source of the parameter database."""


def references(text: str) -> str:
    """The references of a page: every source its text cites, by key.

    A citation is a link to `#ref-<key>` of a source of
    `pkpdsim/database/sources.toml`; the list gives each source once with the
    anchor of its key, in the order of the keys.

    Args:
        text: the markdown of the page.

    Returns:
        The list of the references.

    Raises:
        KeyError: for a citation of a key which is no source.
    """
    known = sources()
    keys = sorted(set(_CITATION.findall(text)), key=str.lower)
    unknown = [k for k in keys if k not in known]
    if unknown:
        raise KeyError(f"The page cites the unknown sources {unknown}.")
    lines = []
    for key in keys:
        source = known[key]
        doi = f" [doi:{source.doi}](https://doi.org/{source.doi})" if source.doi else ""
        lines.append(f'- <span id="ref-{key}"></span>{source.citation}{doi}')
    return "\n".join(lines)


def main() -> None:
    """Write the tables into the pages of `docs/pbpk/` and `docs/database.md`."""
    for page in PAGES:
        text = page.read_text(encoding="utf-8")
        rendered = render_page(text)
        if rendered != text:
            page.write_text(rendered, encoding="utf-8")
            print(f"written {page.relative_to(REPO_DIR)}")


if __name__ == "__main__":
    main()
