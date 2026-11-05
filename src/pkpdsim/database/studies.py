"""The data of the studies of the database, as their sources report them.

A study is a folder of `pkpdsim/database/studies` named by the key of its
source in `sources.toml`: `study.toml` describes the study (the tissue, the
control group, the groups and the measurements) and every measurement is a
TSV file of records, one row per entity (a gene symbol of HGNC, e.g.
`CYP3A4`, `SLCO1B1`) and group with the statistics the source reports. The
records are the numbers of the source, never derived ones: a statistic the
source does not report is empty, a number which disagrees with the others of
its row (a mean outside its range, a CV which is not the SD over the mean)
is kept as printed and has a note. `check` finds such a row without a note.

`study.toml`:

```toml
source = "Drozdzik2021"           # the key of sources.toml, the folder name
tissue = "liver"
control = "controls"             # the control group of every ratio
notes = "..."

[groups.controls]
condition = "control"
description = "histologically normal liver ..."

[groups.cirrhosis_C]
condition = "cirrhosis"
child_pugh = "C"
etiology = "mixed"
description = "..."

[measurements.protein]
file = "protein.tsv"
kind = "protein"                 # protein abundance, or `mrna` expression
unit = "fmol/mg tissue"
method = "..."
location = "Table S3"
```

The columns of a TSV file are `entity` and `group` (required), `name` (the
name in the source), `n`, `n_detected` (the samples above the limit of
quantification), `mean`, `sd`, `se`, `cv_percent`, `median`, `q1`, `q3`
(the quartiles), `min`, `max`, `unit` (of a record whose unit is not that
of its measurement), `ci_low`, `ci_high` (the 95 % confidence interval of the
median) and `notes`. The entity of a protein or an mRNA is its HGNC symbol,
that of the kind `physiology` a key of `QUANTITIES` (e.g. `albumin`, `meld`)
in the unit of the quantity (a measurement of the kind has no `unit`).
A study without a control group (`control` left out) gives levels and no
ratios; a study whose numbers are read from a review states the review as
`via`. Numbers digitized from a figure say so
in their notes, with the resolution of a pixel.
"""

import csv
import io
import math
import re
import tomllib
from dataclasses import dataclass
from functools import cache
from importlib import resources
from importlib.resources.abc import Traversable
from typing import Any

from pkpdsim.database.tables import sources

__all__ = [
    "CHILD_PUGH",
    "CONDITIONS",
    "ETIOLOGIES",
    "KINDS",
    "QUANTITIES",
    "Group",
    "Measurement",
    "Quantity",
    "Record",
    "Study",
    "Summary",
    "check",
    "records",
    "studies",
]

KINDS: dict[str, str] = {
    "protein": "protein abundance per amount of tissue",
    "mrna": "mRNA expression",
    "physiology": "a physiological or clinical quantity of the subjects",
}
"""The kinds of measurement."""


@dataclass(frozen=True)
class Quantity:
    """A quantity of the kind `physiology`.

    Attributes:
        description: what the quantity is.
        unit: the unit of every record of the quantity; a source which reports
            another unit is converted when its numbers are entered, the note of
            the record gives the printed numbers.
    """

    description: str
    unit: str


QUANTITIES: dict[str, Quantity] = {
    "albumin": Quantity("plasma (serum) albumin concentration", "g/L"),
    "hematocrit": Quantity("hematocrit", "%"),
    "bilirubin": Quantity("serum total bilirubin", "mg/dL"),
    "creatinine": Quantity("serum creatinine", "mg/dL"),
    "inr": Quantity("international normalized ratio of the prothrombin time", "1"),
    "meld": Quantity("MELD score (model for end-stage liver disease)", "1"),
    "cardiac_index": Quantity("cardiac output per body surface area", "L/min/m2"),
    "portal_vein_flow": Quantity("blood flow of the portal vein", "mL/min"),
    "portal_perfusion": Quantity(
        "portal perfusion of the liver per volume of liver", "mL/min/dL"
    ),
    "splenic_vein_flow": Quantity("blood flow of the splenic vein", "mL/min"),
    "hepatic_artery_flow": Quantity("blood flow of the hepatic artery", "mL/min"),
    "hepatic_arterial_perfusion": Quantity(
        "arterial perfusion of the liver per volume of liver", "mL/min/dL"
    ),
    "liver_volume_per_bw": Quantity("liver volume per body weight", "cm3/kg"),
    "liver_volume_per_height": Quantity("liver volume per body height", "cm3/m"),
    "hepatocyte_volume_per_bw": Quantity("hepatocyte volume per body weight", "cm3/kg"),
    "functional_liver_volume_per_bw": Quantity(
        "functional liver volume per body weight (99mTc-GSA scintigraphy)", "cm3/kg"
    ),
    "hepatic_ecv": Quantity(
        "extracellular volume fraction of the liver (MRI or CT)", "%"
    ),
}
"""The quantities of the kind `physiology`, the entities of its records, in
the unit of each (`1` dimensionless)."""

CONDITIONS: dict[str, str] = {
    "control": "the control group of a study, the healthy reference of its ratios",
    "cirrhosis": "liver cirrhosis or chronic liver disease",
}
"""The conditions of the groups."""

CHILD_PUGH = ("A", "B", "C")
"""The Child-Pugh classes of cirrhosis."""

ETIOLOGIES: dict[str, str] = {
    "mixed": "the groups of several etiologies together",
    "alcoholic liver disease": "alcoholic liver disease (ALD), alcoholic cirrhosis",
    "hepatitis C": "chronic hepatitis C, HCV cirrhosis",
    "primary biliary cholangitis": "primary biliary cholangitis (PBC)",
    "primary sclerosing cholangitis": "primary sclerosing cholangitis (PSC)",
    "autoimmune hepatitis": "autoimmune hepatitis (AIH)",
}
"""The etiologies of the liver diseases."""

_ENTITY = re.compile(r"^[A-Z][A-Z0-9-]*$")
"""A gene symbol of HGNC, the entity of a protein or mRNA."""

_COLUMNS = (
    "entity",
    "name",
    "group",
    "n",
    "n_detected",
    "mean",
    "sd",
    "se",
    "cv_percent",
    "median",
    "q1",
    "q3",
    "min",
    "max",
    "unit",
    "ci_low",
    "ci_high",
    "notes",
)
"""The columns a TSV file of records may have."""

_GROUP_FIELDS = {"condition", "child_pugh", "etiology", "description", "assumption"}
_MEASUREMENT_FIELDS = {"file", "kind", "unit", "method", "location", "notes"}
_STUDY_FIELDS = {
    "source",
    "via",
    "tissue",
    "control",
    "notes",
    "groups",
    "measurements",
}


@dataclass(frozen=True)
class Group:
    """A group of subjects or samples of a study.

    Attributes:
        id: the key of the group in its study.
        condition: the condition, a key of `CONDITIONS`.
        child_pugh: the Child-Pugh class of a cirrhosis group, `None` without.
        etiology: the etiology of a liver disease, a key of `ETIOLOGIES`.
        description: the description of the group.
        assumption: the assumption on the group which its source does not
            state (e.g. the Child-Pugh class of end-stage livers), `None`
            without.
    """

    id: str
    condition: str
    child_pugh: str | None = None
    etiology: str | None = None
    description: str = ""
    assumption: str | None = None

    def fields(self) -> dict[str, str]:
        """The fields which describe the subjects: condition, class and etiology."""
        out = {"condition": self.condition}
        if self.child_pugh is not None:
            out["child_pugh"] = self.child_pugh
        if self.etiology is not None:
            out["etiology"] = self.etiology
        return out


@dataclass(frozen=True)
class Summary:
    """The statistics of a group as its source reports them.

    Attributes:
        n: the number of samples.
        n_detected: the samples above the limit of quantification.
        mean: the arithmetic mean.
        sd: the standard deviation.
        se: the standard error of the mean.
        cv: the coefficient of variation, a fraction.
        median: the median.
        q1: the first quartile.
        q3: the third quartile.
        min: the smallest value.
        max: the largest value.
        ci: the 95 % confidence interval of the median.
    """

    n: int | None = None
    n_detected: int | None = None
    mean: float | None = None
    sd: float | None = None
    se: float | None = None
    cv: float | None = None
    median: float | None = None
    q1: float | None = None
    q3: float | None = None
    min: float | None = None
    max: float | None = None
    ci: tuple[float, float] | None = None

    def inconsistencies(self) -> list[str]:
        """The statistics which disagree with each other.

        Returns:
            A text per disagreement: a negative statistic, a CV which is not
            the SD over the mean (within the rounding of the source, 2
            percentage points and 5 %), a median or mean outside the range, a
            standard deviation of 0 with several samples, more samples
            detected than measured.
        """
        problems = []
        for name in (
            "n",
            "n_detected",
            "mean",
            "sd",
            "se",
            "cv",
            "median",
            "q1",
            "q3",
            "min",
            "max",
        ):
            v = getattr(self, name)
            if v is not None and v < 0:
                problems.append(f"{name} is negative")
        if self.mean and self.sd is not None and self.cv is not None:
            cv = self.sd / self.mean
            if abs(cv - self.cv) > max(0.02, 0.05 * self.cv):
                problems.append(
                    f"the CV {100 * self.cv:g} % is not SD/mean {100 * cv:.1f} %"
                )
        if (
            self.q1 is not None
            and self.q3 is not None
            and self.median is not None
            and not self.q1 - 1e-9 <= self.median <= self.q3 + 1e-9
        ):
            problems.append(f"the median {self.median:g} is outside the quartiles")
        if self.min is not None and self.max is not None:
            for name in ("mean", "median", "q1", "q3"):
                v = getattr(self, name)
                if v is not None and not self.min - 1e-9 <= v <= self.max + 1e-9:
                    problems.append(f"the {name} {v:g} is outside the range")
        if self.sd == 0 and (self.n or 0) > 1 and self.mean:
            problems.append("the SD is 0 with several samples")
        if (
            self.n is not None
            and self.n_detected is not None
            and self.n_detected > self.n
        ):
            problems.append("more samples detected than measured")
        return problems


@dataclass(frozen=True)
class Record:
    """A record of a study: the statistics of an entity in a group.

    Attributes:
        study: the key of the study, the key of its source.
        measurement: the key of the measurement in its study.
        kind: the kind of the measurement, a key of `KINDS`.
        tissue: the tissue.
        unit: the unit of the values.
        entity: the gene symbol of the entity.
        name: the name of the entity in the source.
        group: the group.
        summary: the statistics.
        notes: a remark on the record, e.g. why a statistic disagrees.
    """

    study: str
    measurement: str
    kind: str
    tissue: str
    unit: str
    entity: str
    name: str
    group: Group
    summary: Summary
    notes: str = ""


@dataclass(frozen=True)
class Measurement:
    """A measurement of a study: a kind of quantity in a tissue, one TSV file.

    Attributes:
        key: the key of the measurement in its study.
        kind: the kind, a key of `KINDS`.
        unit: the unit of the values (a record may state its own); empty for
            the kind `physiology`, whose records are in the unit of their
            quantity.
        method: the method of the measurement.
        location: where the values are in the source, e.g. `Table S3`.
        notes: a remark on the measurement.
        records: the records.
    """

    key: str
    kind: str
    unit: str
    method: str
    location: str
    notes: str
    records: tuple[Record, ...]


@dataclass(frozen=True)
class Study:
    """A study: its groups and measurements.

    Attributes:
        key: the key of the study, the key of its source in `sources.toml`.
        tissue: the tissue of the measurements, `body` for the quantities of
            the subjects.
        control: the control group, the reference of every ratio, `None` for a
            study without (its records give levels, no ratios).
        groups: the groups by key.
        measurements: the measurements by key.
        notes: a remark on the study.
        via: the key of the source the numbers are read from when it is not
            the study itself (a table of a review), `None` for the study.
    """

    key: str
    tissue: str
    control: str | None
    groups: dict[str, Group]
    measurements: dict[str, Measurement]
    notes: str
    via: str | None = None


def _root() -> Traversable:
    """The folder of the studies."""
    return resources.files("pkpdsim.database").joinpath("studies")


def _number(text: str, where: str) -> float | None:
    """A number of a TSV file, `None` for an empty cell."""
    text = text.strip()
    if not text:
        return None
    try:
        value = float(text)
    except ValueError:
        raise ValueError(f"{where}: '{text}' is no number.") from None
    if not math.isfinite(value):
        raise ValueError(f"{where}: '{text}' is not finite.")
    return value


def _integer(text: str, where: str) -> int | None:
    """An integer of a TSV file, `None` for an empty cell."""
    value = _number(text, where)
    if value is None:
        return None
    if value != int(value):
        raise ValueError(f"{where}: '{text}' is no integer.")
    return int(value)


def _group(key: str, data: dict[str, Any], where: str) -> Group:
    """A group of `study.toml`."""
    unknown = set(data) - _GROUP_FIELDS
    if unknown:
        raise ValueError(f"{where}: unknown fields {sorted(unknown)} of group '{key}'.")
    condition = data.get("condition")
    if condition not in CONDITIONS:
        raise ValueError(
            f"{where}: group '{key}' has the unknown condition {condition!r}."
        )
    child_pugh = data.get("child_pugh")
    if child_pugh is not None and child_pugh not in CHILD_PUGH:
        raise ValueError(
            f"{where}: group '{key}' has the unknown class {child_pugh!r}."
        )
    etiology = data.get("etiology")
    if etiology is not None and etiology not in ETIOLOGIES:
        raise ValueError(
            f"{where}: group '{key}' has the unknown etiology {etiology!r}."
        )
    if condition == "control" and (child_pugh or etiology):
        raise ValueError(f"{where}: the control group '{key}' has a class or etiology.")
    return Group(
        id=key,
        condition=condition,
        child_pugh=child_pugh,
        etiology=etiology,
        description=data.get("description", ""),
        assumption=data.get("assumption"),
    )


def _records(
    text: str,
    study: str,
    key: str,
    data: dict[str, Any],
    tissue: str,
    groups: dict[str, Group],
) -> tuple[Record, ...]:
    """The records of a TSV file."""
    where = f"studies/{study}/{data['file']}"
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    columns = reader.fieldnames or []
    unknown = set(columns) - set(_COLUMNS)
    if unknown or not {"entity", "group"} <= set(columns):
        raise ValueError(
            f"{where}: the columns {columns} are not columns of the records "
            f"{list(_COLUMNS)} with entity and group."
        )
    out = []
    seen = set()
    for line, row in enumerate(reader, start=2):
        at = f"{where}:{line}"
        entity, group = row["entity"].strip(), row["group"].strip()
        unit = (row.get("unit") or "").strip() or data.get("unit", "")
        if data["kind"] == "physiology":
            if entity not in QUANTITIES:
                raise ValueError(f"{at}: '{entity}' is no quantity of QUANTITIES.")
            expected = QUANTITIES[entity].unit
            if unit and unit != expected:
                raise ValueError(
                    f"{at}: {entity} in '{unit}', the unit of the quantity is "
                    f"'{expected}' (convert the numbers, the note gives the printed)."
                )
            unit = expected
        elif not unit:
            raise ValueError(f"{at}: the record has no unit.")
        elif not _ENTITY.match(entity):
            raise ValueError(f"{at}: '{entity}' is no gene symbol.")
        if group not in groups:
            raise ValueError(f"{at}: the group '{group}' is not a group of the study.")
        if (entity, group) in seen:
            raise ValueError(f"{at}: a second record of {entity} in '{group}'.")
        seen.add((entity, group))
        ci_low = _number(row.get("ci_low", "") or "", at)
        ci_high = _number(row.get("ci_high", "") or "", at)
        if (ci_low is None) != (ci_high is None):
            raise ValueError(f"{at}: an interval needs both bounds.")
        cv = _number(row.get("cv_percent", "") or "", at)
        summary = Summary(
            n=_integer(row.get("n", "") or "", at),
            n_detected=_integer(row.get("n_detected", "") or "", at),
            mean=_number(row.get("mean", "") or "", at),
            sd=_number(row.get("sd", "") or "", at),
            se=_number(row.get("se", "") or "", at),
            cv=cv / 100 if cv is not None else None,
            median=_number(row.get("median", "") or "", at),
            q1=_number(row.get("q1", "") or "", at),
            q3=_number(row.get("q3", "") or "", at),
            min=_number(row.get("min", "") or "", at),
            max=_number(row.get("max", "") or "", at),
            ci=(ci_low, ci_high)
            if ci_low is not None and ci_high is not None
            else None,
        )
        out.append(
            Record(
                study=study,
                measurement=key,
                kind=data["kind"],
                tissue=tissue,
                unit=unit,
                entity=entity,
                name=(row.get("name") or "").strip() or entity,
                group=groups[group],
                summary=summary,
                notes=(row.get("notes") or "").strip(),
            )
        )
    return tuple(out)


def _study(folder: Traversable) -> Study:
    """A study of its folder."""
    key = folder.name
    where = f"studies/{key}/study.toml"
    data = tomllib.loads(folder.joinpath("study.toml").read_text("utf-8"))
    unknown = set(data) - _STUDY_FIELDS
    if unknown:
        raise ValueError(f"{where}: unknown fields {sorted(unknown)}.")
    if data.get("source") != key:
        raise ValueError(
            f"{where}: the source {data.get('source')!r} is not the folder '{key}'."
        )
    if key not in sources():
        raise ValueError(f"{where}: '{key}' is no source of sources.toml.")
    groups = {k: _group(k, g, where) for k, g in data.get("groups", {}).items()}
    control = data.get("control")
    if control is not None and (
        control not in groups or groups[control].condition != "control"
    ):
        raise ValueError(f"{where}: the control {control!r} is no control group.")
    via = data.get("via")
    if via is not None and via not in sources():
        raise ValueError(f"{where}: the source {via!r} of `via` is unknown.")
    measurements = {}
    for k, m in data.get("measurements", {}).items():
        unknown = set(m) - _MEASUREMENT_FIELDS
        if unknown:
            raise ValueError(
                f"{where}: unknown fields {sorted(unknown)} of measurement '{k}'."
            )
        if m.get("kind") not in KINDS:
            raise ValueError(
                f"{where}: measurement '{k}' has the unknown kind {m.get('kind')!r}."
            )
        text = folder.joinpath(m["file"]).read_text("utf-8")
        measurements[k] = Measurement(
            key=k,
            kind=m["kind"],
            unit=m.get("unit", ""),
            method=m.get("method", ""),
            location=m.get("location", ""),
            notes=m.get("notes", ""),
            records=_records(text, key, k, m, data["tissue"], groups),
        )
    return Study(
        key=key,
        tissue=data["tissue"],
        control=control,
        groups=groups,
        measurements=measurements,
        notes=data.get("notes", ""),
        via=via,
    )


@cache
def studies() -> dict[str, Study]:
    """Get the studies of the database by key, read once per process."""
    return {
        folder.name: _study(folder)
        for folder in sorted(_root().iterdir(), key=lambda p: p.name)
        if folder.is_dir() and folder.joinpath("study.toml").is_file()
    }


def records(
    kind: str | None = None, entity: str | None = None, tissue: str | None = None
) -> list[Record]:
    """Get the records of the studies.

    Args:
        kind: only the records of this kind of measurement.
        entity: only the records of this entity.
        tissue: only the records of this tissue.

    Returns:
        The records in the order of the studies and their files.
    """
    return [
        r
        for study in studies().values()
        for m in study.measurements.values()
        for r in m.records
        if (kind is None or r.kind == kind)
        and (entity is None or r.entity == entity)
        and (tissue is None or r.tissue == tissue)
    ]


def check() -> list[str]:
    """Check the records of the studies.

    Returns:
        A problem per record whose statistics disagree and which has no note
        that says why, empty if every record is consistent or explained.
    """
    problems = []
    for r in records():
        found = r.summary.inconsistencies()
        if found and not r.notes:
            problems.append(
                f"studies/{r.study} {r.measurement} {r.entity} {r.group.id}: "
                f"{'; '.join(found)}, without a note"
            )
    return problems
