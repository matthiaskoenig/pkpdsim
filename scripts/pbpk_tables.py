"""Write the data tables of the pages of the PBPK library from its data files.

The pages of `docs/pbpk/` show the reference physiology, the populations, the
impairment vectors, the arm vein presets, the intestine, the nephron and the
proteins as tables. The numbers are those of the TOML files of
`pkpdsim/pbpk/data`, never typed by hand: a page marks the region of a table
with

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

import re
from collections.abc import Callable, Iterable, Mapping
from pathlib import Path
from typing import Any

from sbmlutils.factory import UnitDefinition

from pkpdsim.pbpk.data import load, sources
from pkpdsim.pbpk.organs import Organ
from pkpdsim.pbpk.physiology import Sex
from pkpdsim.pbpk.populations import POPULATIONS

REPO_DIR: Path = Path(__file__).parent.parent
PBPK_DIR: Path = REPO_DIR / "docs" / "pbpk"

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
        "the notes of every number are in the data file.*"
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
    return _generated("`data/physiology_human.toml`", lines, cited)


def _organs(sex: Sex) -> str:
    """The organs of the reference individual of a sex."""
    organs = load("physiology_human")[sex.value]["organs"]
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
    return _generated("`data/physiology_human.toml`", lines, cited)


def physiology_organs_male() -> str:
    """The organs of the male reference individual."""
    return _organs(Sex.MALE)


def physiology_organs_female() -> str:
    """The organs of the female reference individual."""
    return _organs(Sex.FEMALE)


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
            covariates = population.covariates[sex]
            spread = [
                f"{name} {d.kind} {_number(d.mean)} &plusmn; {_number(d.sd)}"
                for name, d in covariates.items()
            ]
            cited.append(population.source)
            cited += [v.source for v in reference.values() if v.source]
            cited += [d.source for d in covariates.values()]
            rows.append(
                [
                    f"{population.name} `{key}`",
                    sex.value,
                    _number(reference["BW"].value),
                    _number(reference["HEIGHT"].value),
                    _number(reference["AGE"].value),
                    "<br>".join(spread) if spread else "none",
                ]
            )
    reference = POPULATIONS[default].reference[Sex.MALE]
    header = [
        "population",
        "sex",
        *(f"{c} [{_unit_id(reference[c].unit)}]" for c in ("BW", "HEIGHT", "AGE")),
        "covariates",
    ]
    lines = _table(header, rows)
    return _generated("`data/populations/*.toml`", lines, cited)


def _vector_table(
    section: Mapping[str, Any], classes: list[str], first: str = "factor"
) -> tuple[list[str], list[str]]:
    """Format a table of vectors by class; a model choice is marked with `*`."""
    rows: list[list[str]] = []
    cited: list[str] = []
    for key, entry in section.items():
        if not isinstance(entry, dict) or "value" not in entry:
            continue
        choice = set(entry.get("choice", ()))
        values = [
            _number(v) + ("*" if c in choice else "")
            for v, c in zip(entry["value"], classes, strict=True)
        ]
        cited.append(entry["source"])
        rows.append([f"`{key}`", *values, f"`{entry['source']}`"])
    return _table([first, *classes, "source"], rows), cited


_CHOICE = "A value marked `*` is a choice of the library, not a value of the source."


def impairment_child_pugh() -> str:
    """The factors of the Child-Pugh classes on the physiology."""
    lines, cited = _vector_table(load("impairment")["child_pugh"], ["A", "B", "C"])
    return _generated("`data/impairment.toml`", lines, cited)


def impairment_child_pugh_proteins() -> str:
    """The factors of the Child-Pugh classes on the proteins of the liver."""
    section = load("impairment")["child_pugh"]["proteins"]
    lines, cited = _vector_table(section, ["A", "B", "C"], first="protein")
    return _generated("`data/impairment.toml`", [*lines, "", _CHOICE], cited)


def impairment_ckd() -> str:
    """The factors of the CKD stages on the physiology."""
    section = dict(load("impairment")["ckd"])
    stages = section.pop("stage_gfr")
    lines, cited = _vector_table(section, ["mild", "moderate", "severe", "end stage"])
    limits = ", ".join(_number(v) for v in stages["value"])
    head = [
        f"Stage limits of the GFR per 1.73 m2 ({stages['unit']}): {limits} "
        f"(`{stages['source']}`).",
        "",
    ]
    return _generated(
        "`data/impairment.toml`",
        [*head, *lines, "", _CHOICE],
        [*cited, stages["source"]],
    )


def impairment_ckd_proteins() -> str:
    """The factors of the CKD stages on the proteins of the liver (`CKD_hepatic`)."""
    section = load("impairment")["ckd"]["proteins"]
    stages = ["mild", "moderate", "severe", "end stage"]
    lines, cited = _vector_table(section, stages, first="protein")
    return _generated("`data/impairment.toml`", [*lines, "", _CHOICE], cited)


def impairment_heart_failure() -> str:
    """The factors of the NYHA classes on the blood flows."""
    section = load("impairment")["heart_failure"]
    lines, cited = _vector_table(section, ["I", "II", "III", "IV"])
    return _generated("`data/impairment.toml`", [*lines, "", _CHOICE], cited)


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
    return _generated("`data/arm_vein.toml`", lines, cited)


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
    return _generated("`data/intestine_human.toml`", lines, cited)


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
    return _generated("`data/kidney_human.toml`", lines, cited)


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
    return _generated("`data/proteins.toml`", lines, cited)


TABLES: dict[str, Callable[[], str]] = {
    "physiology-general": physiology_general,
    "physiology-organs-male": physiology_organs_male,
    "physiology-organs-female": physiology_organs_female,
    "populations": populations,
    "impairment-child-pugh": impairment_child_pugh,
    "impairment-child-pugh-proteins": impairment_child_pugh_proteins,
    "impairment-ckd": impairment_ckd,
    "impairment-ckd-proteins": impairment_ckd_proteins,
    "impairment-heart-failure": impairment_heart_failure,
    "arm-vein": arm_vein,
    "intestine-segments": intestine_segments,
    "nephron-segments": nephron_segments,
    "proteins": proteins,
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
        if name not in TABLES:
            raise KeyError(f"No table '{name}', the tables are {sorted(TABLES)}.")
        return f"{begin_marker(name)}\n{TABLES[name]()}\n{END}"

    return _REGION.sub(replace, text)


def main() -> None:
    """Write the tables into the pages of `docs/pbpk/`."""
    for page in sorted(PBPK_DIR.glob("*.md")):
        text = page.read_text(encoding="utf-8")
        rendered = render_page(text)
        if rendered != text:
            page.write_text(rendered, encoding="utf-8")
            print(f"written {page.relative_to(REPO_DIR)}")


if __name__ == "__main__":
    main()
