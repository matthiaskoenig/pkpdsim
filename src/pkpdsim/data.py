"""The study data of a model, collected from the study folders of PK-DB.

The studies of [PK-DB](https://pk-db.com) are curated in a checkout of
[pkdb_data](https://github.com/matthiaskoenig/pkdb_data): the folder `studies`
holds a folder per substance and in it a folder per study. A model keeps the
subset of the studies it uses in its own `data/` folder, so its simulations
run without the checkout; `collect_tsv_files` writes them there from the
folder `studies_dir()` finds.

The experiments of the models read the tables of study format 1, one hidden
`.tsv` file per source of a study (`.Boulton1996a_Fig1.tsv`, see
`pkpdsim.experiment.PKPDSimulationExperiment.load_dataframe`). pkdb_data is
curated in study format 2 (`"format": 2` of the `study.json`, the tables
`subjects.tsv`, `interventions.tsv`, `timecourses_<source>.tsv`,
`outputs_<source>.tsv`, ...). A study of format 2 is read with the library of
PK-DB (`pkdb.prepare`, the same reading and validation as `pkdb check` and the
upload to the database), and its data are written as the tables of format 1
(`study_tables`), so that the experiments read them unchanged: a geometric
mean in `gmean` (calculation `geometric mean`), a value below the limit of
quantification as a row without a central value, an output with a label
generated from its columns which is unique across the sources of the study,
the interventions in the natural order of format 2 (`GEM,SIM`). Format 2 has
no collection intervals (no `tstart` and `tend` of format 1) and no columns
of its own of a table, so neither have the tables. A study of format 1 is
copied as it is. `database_differences` compares the studies with the ones
uploaded to a PK-DB server.
"""

import json
import os
import re
import shutil
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd
from sbmlutils.console import console

if TYPE_CHECKING:
    from pkdb import Client
    from pkdb.schemas.study import (
        CanonicalStudy,
        Measurement,
        Notes,
        Observation,
        Subject,
    )

PKDB_DATA = "PKDB_DATA"
"""The environment variable with the path of a checkout of pkdb_data."""


def studies_dir() -> Path:
    """Get the folder of the PK-DB studies, `studies` of a checkout of pkdb_data.

    The checkout is the path of the environment variable `PKDB_DATA`; without
    the variable it is `../pkdb_data` relative to the working directory, i.e.
    a checkout next to the repository a script runs in.

    Returns:
        The absolute path of the folder `studies` of the checkout.

    Raises:
        FileNotFoundError: if the folder does not exist; the message names the
            environment variable and the default location.
    """
    default = (Path.cwd().parent / "pkdb_data").resolve()
    value = os.environ.get(PKDB_DATA)
    if value:
        path = Path(value).expanduser().resolve() / "studies"
        if path.is_dir():
            return path
        raise FileNotFoundError(
            f"The PK-DB studies '{path}' do not exist: the environment variable "
            f"{PKDB_DATA}='{value}' has to be the path of a checkout of pkdb_data "
            f"(https://github.com/matthiaskoenig/pkdb_data). Without "
            f"{PKDB_DATA} the checkout '{default}' is used."
        )

    path = default / "studies"
    if path.is_dir():
        return path
    raise FileNotFoundError(
        f"The PK-DB studies '{path}' do not exist: set the environment variable "
        f"{PKDB_DATA} to the path of a checkout of pkdb_data "
        f"(https://github.com/matthiaskoenig/pkdb_data), or clone it next to the "
        f"working directory, i.e. to '{default}'."
    )


def study_format(study_folder: Path) -> int:
    """Get the study format of a study folder of pkdb_data.

    Args:
        study_folder: the folder of a study.

    Returns:
        The `format` of its `study.json`; 1 for a study without one, or
        without the field, which is study format 1.

    Raises:
        ValueError: if the `study.json` is no JSON object or its format is no
            integer.
    """
    path = study_folder / "study.json"
    if not path.is_file():
        return 1
    try:
        study = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as err:
        raise ValueError(f"The study '{path}' is no JSON: {err}") from err
    if not isinstance(study, dict):
        raise ValueError(f"The study '{path}' is no JSON object.")
    value = study.get("format", 1)
    if not isinstance(value, int) or isinstance(value, bool):
        raise ValueError(f"The format of the study '{path}' is no integer: {value!r}.")
    return value


#: the columns of a table of measurements of study format 1 (`.<study>_<source>.tsv`)
MEASUREMENT_COLUMNS: tuple[str, ...] = (
    "study",
    "group",
    "individual",
    "count",
    "intervention",
    "measurement_type",
    "calculation_type",
    "choice",
    "substance",
    "tissue",
    "method",
    "label",
    "time",
    "time_unit",
    "value",
    "mean",
    "median",
    "sd",
    "se",
    "cv",
    "min",
    "max",
    "gmean",
    "gsd",
    "gcv",
    "unit",
    "comments",
)

#: the columns of the groups of study format 1 (`.<study>_TabGroups.tsv`)
GROUP_COLUMNS: tuple[str, ...] = (
    "study",
    "name",
    "group_count",
    "parent",
    "measurement_type",
    "count",
    "substance",
    "choice",
    "mean",
    "median",
    "sd",
    "se",
    "min",
    "max",
    "unit",
    "comments",
)

#: the columns of the individuals of study format 1 (`.<study>_TabIndividuals.tsv`)
INDIVIDUAL_COLUMNS: tuple[str, ...] = (
    "study",
    "name",
    "group",
    "measurement_type",
    "count",
    "substance",
    "choice",
    "value",
    "unit",
    "comments",
)

#: the columns of the interventions of study format 1 (`.<study>_TabInterventions.tsv`),
#: with the `interval` and the number of `doses` of a repeated dose of format 2
INTERVENTION_COLUMNS: tuple[str, ...] = (
    "study",
    "name",
    "time",
    "time_end",
    "time_unit",
    "interval",
    "doses",
    "route",
    "application",
    "measurement_type",
    "form",
    "substance",
    "value",
    "unit",
    "comments",
)

#: the statistics of a measurement which are columns of the tables of format 1
STATISTICS: tuple[str, ...] = (
    "mean",
    "median",
    "sd",
    "se",
    "cv",
    "min",
    "max",
    "gmean",
    "gsd",
    "gcv",
)

#: the prefixes of the data tables of study format 2, `<prefix>_<source>.tsv`
DATA_PREFIXES: tuple[str, ...] = ("timecourses_", "outputs_", "scatters_")

#: the columns of a dimension of a scatter of study format 1, `x_<column>` and
#: `y_<column>`; `""` is the column `x` (`y`) of the value, which `x_unit` is the unit of
SCATTER_DIMENSION_COLUMNS: tuple[str, ...] = (
    "measurement_type",
    "substance",
    "tissue",
    "method",
    "label",
    "value",
    "",
    "unit",
)

#: the columns of a scatter of study format 1 (`.<study>_<source>.tsv` of a scatter)
SCATTER_COLUMNS: tuple[str, ...] = (
    "study",
    "individual",
    "group",
    "count",
    "intervention",
    "label",
    "time",
    "time_unit",
    *("_".join(("x", column)).rstrip("_") for column in SCATTER_DIMENSION_COLUMNS),
    *("_".join(("y", column)).rstrip("_") for column in SCATTER_DIMENSION_COLUMNS),
    "comments",
)


def _comments(record: Notes) -> str | None:
    """Get the comments of a record of PK-DB as one text."""
    texts = [comment.text for comment in record.comments if comment.text]
    return " | ".join(texts) if texts else None


def _source(file: str) -> str:
    """Get the source of a data table of format 2, `Fig1` of `timecourses_Fig1.tsv`."""
    stem = Path(file).stem
    for prefix in DATA_PREFIXES:
        if stem.startswith(prefix):
            return stem.removeprefix(prefix)
    return stem


def _statistics(
    observation: Observation, *, complete: bool = False, count: int | None = None
) -> dict[str, Any]:
    """Get the statistics of an observation as the columns of format 1.

    A digitized error bar of format 2 (`error_bar`, the position of the end of
    the bar, with its `error_type`) gives the statistic it names, which format
    1 stored as the size of the bar: `sd` or `se` is the distance of the end
    from the mean, `gsd` the ratio of the end and the geometric mean. With
    `complete`, the missing `sd`, `se` and `cv` are calculated from the others
    and the count as PK-DB does (`pkdb.domain.statistics.complete_statistics`,
    `cv` is a fraction of the mean). A reported value is always kept.

    Args:
        observation: a measurement or characteristic of PK-DB.
        complete: whether to calculate the missing statistics.
        count: the count of the subjects of the observation, for `complete`.

    Returns:
        The statistics by the columns of format 1.
    """
    from pkdb.domain.statistics import complete_statistics

    statistics = observation.statistics
    if complete:
        statistics = complete_statistics(statistics, count)
    values = {key: getattr(statistics, key) for key in STATISTICS}
    error_bar, error_type = statistics.error_bar, statistics.error_type
    if error_bar is not None and error_type is not None and values[error_type] is None:
        if error_type == "gsd":
            gmean = values["gmean"]
            if gmean is not None and gmean > 0 and error_bar > 0:
                values["gsd"] = max(error_bar / gmean, gmean / error_bar)
        elif values["mean"] is not None:
            values[error_type] = abs(error_bar - values["mean"])
    return values


def _label(*parts: str | None) -> str:
    """Get the label of an output without one from the columns which identify it.

    The parts which are given are joined by `_`, every other character than a
    letter or a digit becomes `_` (`auc_hydrochlorothiazide_plasma_HCTZ25_all`).
    """
    text = "_".join(part for part in parts if part)
    return re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")


#: the calculation PK-DB gives every measurement of a group without one
SAMPLE_MEAN = "sample mean"
#: the calculation of a measurement of a group with a geometric mean but no mean
GEOMETRIC_MEAN = "geometric mean"


def _calculation_type(measurement: Measurement) -> str | None:
    """Get the calculation of a measurement as the column of format 1.

    PK-DB gives every measurement of a group without a `calculation` the
    calculation `sample mean` (`pkdb.prepare`), also one with a geometric mean
    and no mean: such a measurement is a `geometric mean`.
    """
    statistics = measurement.statistics
    if (
        measurement.calculation_type == SAMPLE_MEAN
        and statistics.mean is None
        and statistics.gmean is not None
    ):
        return GEOMETRIC_MEAN
    return measurement.calculation_type


def _measurement_tables(study: CanonicalStudy, name: str) -> dict[str, pd.DataFrame]:
    """Get the measurements of a study as the tables of format 1, by source.

    A table holds the reported measurements of a source, the outputs and the
    points of the timecourses of `outputs_<source>.tsv` and
    `timecourses_<source>.tsv`; the records PK-DB normalizes into other units
    or calculates from them are left out, and the missing `sd`, `se` and `cv`
    of a row are calculated from its statistics (see `_statistics`). A
    measurement of an individual has its mean in the column `value` and its
    group in `group`, the count of a row without one is the count of its
    subjects. A row without a central value (a value below the limit of
    quantification, curated with only its `max`) has none in `value`, `mean`,
    `median` and `gmean`. A row with a geometric mean but no mean has the
    calculation `geometric mean` (`_calculation_type`). The interventions are
    the names of the cell of format 2, comma separated in the order of the
    cell, which is the natural order of `pkdb format` (`GEM,SIM`).

    An output without a label (the outputs of format 2 have none) gets one
    from its measurement, substance, tissue, interventions and subject
    (`_label`), so that the experiments can group the rows by `label` as for
    the timecourses. A label is unique in the study: a generated label which
    another source of the study has as well gets the source appended
    (`concentration_sorafenib_plasma_SOF400_1_Fig1` and `..._Tab3`), in every
    source which has it.
    """
    group_counts = {group.name: group.count for group in study.groups}
    individual_groups = {
        individual.name: individual.group for individual in study.individuals
    }
    rows: dict[str, list[dict[str, Any]]] = {}
    generated_rows: list[tuple[str, dict[str, Any]]] = []
    for measurement in study.measurements:
        if (
            measurement.origin != "reported"
            or measurement.source is None
            or measurement.source.file.startswith("scatters_")
        ):
            continue
        individual = measurement.individual
        group = measurement.group
        count = measurement.statistics.count
        if individual is not None:
            group = group or individual_groups.get(individual)
            count = 1 if count is None else count
        elif count is None and group is not None:
            count = group_counts.get(group)
        statistics = _statistics(measurement, complete=True, count=count)
        value = None
        if individual is not None:
            value = statistics["mean"]
            statistics["mean"] = None
        interventions = ",".join(measurement.interventions) or None
        generated = not measurement.label
        label = measurement.label or _label(
            measurement.measurement_type,
            measurement.substance,
            measurement.tissue,
            "_".join(measurement.interventions),
            individual or group,
        )
        row = {
            "study": name,
            "group": group,
            "individual": individual,
            "count": count,
            "intervention": interventions,
            "measurement_type": measurement.measurement_type,
            "calculation_type": _calculation_type(measurement),
            "choice": measurement.choice,
            "substance": measurement.substance,
            "tissue": measurement.tissue,
            "method": measurement.method,
            "label": label,
            "time": None if measurement.time_not_reported else measurement.time,
            "time_unit": measurement.time_unit,
            "value": value,
            **statistics,
            "unit": measurement.unit,
            "comments": _comments(measurement),
        }
        source = _source(measurement.source.file)
        rows.setdefault(source, []).append(row)
        if generated:
            generated_rows.append((source, row))

    # the sources of a label, before a generated label gets its source
    sources: dict[str, set[str]] = {}
    for source, source_rows in rows.items():
        for row in source_rows:
            sources.setdefault(row["label"], set()).add(source)
    for source, row in generated_rows:
        if len(sources[row["label"]]) > 1:
            row["label"] = _label(row["label"], source)
    return {
        source: pd.DataFrame(source_rows, columns=list(MEASUREMENT_COLUMNS))
        for source, source_rows in rows.items()
    }


def _scatter_tables(study: CanonicalStudy, name: str) -> dict[str, pd.DataFrame]:
    """Get the scatters of a study as the tables of format 1, by source.

    A point of a scatter (a subject, the row of `scatters_<source>.tsv`) is a
    row with its two measurements in the columns `x_<column>` and
    `y_<column>`, the value in `x_value` and `x` (the column of sbmlsim which
    `x_unit` is the unit of) and in `y_value` and `y`, and the count 1; the
    label is the name of the scatter, the interventions are those of both
    measurements and the time is the one of `y` (else of `x`).
    """
    reported = {
        measurement.key: measurement
        for measurement in study.measurements
        if measurement.origin == "reported"
    }
    every = {measurement.key: measurement for measurement in study.measurements}

    def point(key: str) -> Measurement:
        measurement = every[key]
        if measurement.origin != "reported" and measurement.derived_from:
            measurement = reported[measurement.derived_from]
        return measurement

    rows: dict[str, list[dict[str, Any]]] = {}
    for scatter in study.scatters:
        if scatter.data_type != "scatter" or scatter.source is None:
            continue
        for subset in scatter.subsets:
            for keys in subset.points:
                x, y = (point(key) for key in keys)
                interventions = list(
                    dict.fromkeys([*x.interventions, *y.interventions])
                )
                row: dict[str, Any] = {
                    "study": name,
                    "individual": y.individual,
                    "group": y.group,
                    "count": 1,
                    "intervention": ",".join(interventions) or None,
                    "label": scatter.name,
                    "time": y.time if y.time is not None else x.time,
                    "time_unit": y.time_unit or x.time_unit,
                    "comments": _comments(y) or _comments(x),
                }
                for axis, measurement in (("x", x), ("y", y)):
                    row |= {
                        f"{axis}_measurement_type": measurement.measurement_type,
                        f"{axis}_substance": measurement.substance,
                        f"{axis}_tissue": measurement.tissue,
                        f"{axis}_method": measurement.method,
                        f"{axis}_label": measurement.label,
                        f"{axis}_value": measurement.statistics.mean,
                        axis: measurement.statistics.mean,
                        f"{axis}_unit": measurement.unit,
                    }
                rows.setdefault(_source(scatter.source.file), []).append(row)
    return {
        source: pd.DataFrame(source_rows, columns=list(SCATTER_COLUMNS))
        for source, source_rows in rows.items()
    }


def _characteristica_rows(subject: Subject) -> list[Observation | None]:
    """Get the reported characteristics of a subject, one `None` for none."""
    observations = [
        observation
        for observation in subject.characteristica
        if observation.origin == "reported"
    ]
    return list(observations) if observations else [None]


def _group_table(study: CanonicalStudy, name: str) -> pd.DataFrame:
    """Get the groups of a study and their characteristics as the table of format 1."""
    rows = []
    for group in study.groups:
        for observation in _characteristica_rows(group):
            row: dict[str, Any] = {
                "study": name,
                "name": group.name,
                "group_count": group.count,
                "parent": group.parent,
            }
            if observation is not None:
                statistics = _statistics(observation)
                row |= {
                    "measurement_type": observation.measurement_type,
                    "count": observation.statistics.count,
                    "substance": observation.substance,
                    "choice": observation.choice,
                    **{
                        key: statistics[key]
                        for key in GROUP_COLUMNS
                        if key in statistics
                    },
                    "unit": observation.unit,
                    "comments": _comments(observation),
                }
            rows.append(row)
    return pd.DataFrame(rows, columns=list(GROUP_COLUMNS))


def _individual_table(study: CanonicalStudy, name: str) -> pd.DataFrame:
    """Get the individuals of a study and their characteristics as the table of format 1."""
    rows = []
    for individual in study.individuals:
        for observation in _characteristica_rows(individual):
            row: dict[str, Any] = {
                "study": name,
                "name": individual.name,
                "group": individual.group,
            }
            if observation is not None:
                row |= {
                    "measurement_type": observation.measurement_type,
                    "count": observation.statistics.count,
                    "substance": observation.substance,
                    "choice": observation.choice,
                    "value": observation.statistics.mean,
                    "unit": observation.unit,
                    "comments": _comments(observation),
                }
            rows.append(row)
    return pd.DataFrame(rows, columns=list(INDIVIDUAL_COLUMNS))


def _intervention_table(study: CanonicalStudy, name: str) -> pd.DataFrame:
    """Get the reported interventions of a study as the table of format 1.

    The doses PK-DB normalizes into other units are left out.
    """
    rows = []
    for intervention in study.interventions:
        if intervention.origin != "reported":
            continue
        time = intervention.time
        if isinstance(time, list):
            time = ",".join(str(value) for value in time)
        rows.append(
            {
                "study": name,
                "name": intervention.name,
                "time": None if intervention.time_not_reported else time,
                "time_end": intervention.time_end,
                "time_unit": intervention.time_unit,
                "interval": intervention.interval,
                "doses": intervention.doses,
                "route": intervention.route,
                "application": intervention.application,
                "measurement_type": intervention.measurement_type,
                "form": intervention.form,
                "substance": intervention.substance,
                "value": intervention.statistics.mean,
                "unit": intervention.unit,
                "comments": _comments(intervention),
            }
        )
    return pd.DataFrame(rows, columns=list(INTERVENTION_COLUMNS))


def prepare_study(study_folder: Path) -> CanonicalStudy:
    """Read a study folder of format 2 with the library of PK-DB.

    The folder is read and validated as `pkdb check` and the upload to the
    database do (`pkdb.prepare`), with the vocabulary of the checkout
    (`vocabulary.lock.json` of pkdb_data), else the one bundled with pkdb.

    Args:
        study_folder: the folder of a study of format 2.

    Returns:
        The study as PK-DB reads it.
    """
    from pkdb import prepare
    from pkdb.cache import load_vocabulary, lock_file

    lock = lock_file(study_folder)
    vocabulary = load_vocabulary(lock) if lock is not None else None
    return prepare(study_folder, vocabulary=vocabulary).study


def study_tables(study_folder: Path) -> dict[str, pd.DataFrame]:
    """Get the data of a study of format 2 as the tables of format 1.

    The tables are keyed by the source of format 1, the part of the file name
    after the study (`Fig1` of `.Boulton1996a_Fig1.tsv`): a table per source
    of the measurements (the reported outputs and timecourses of
    `outputs_<source>.tsv` and `timecourses_<source>.tsv`, `MEASUREMENT_COLUMNS`),
    a table per source of the scatters of `scatters_<source>.tsv`
    (`SCATTER_COLUMNS`; `<source>Scatter` when the source also has
    measurements), and `TabGroups` (`GROUP_COLUMNS`), `TabIndividuals`
    (`INDIVIDUAL_COLUMNS`) and `TabInterventions` (`INTERVENTION_COLUMNS`)
    with the subjects, their characteristics and the interventions. The
    columns are those of format 1, so the experiments select their rows as
    before.

    Args:
        study_folder: the folder of a study of format 2.

    Returns:
        The tables by source.
    """
    name = study_folder.name
    study = prepare_study(study_folder)
    tables = _measurement_tables(study, name)
    for source, df in _scatter_tables(study, name).items():
        tables[f"{source}Scatter" if source in tables else source] = df
    tables["TabGroups"] = _group_table(study, name)
    tables["TabIndividuals"] = _individual_table(study, name)
    tables["TabInterventions"] = _intervention_table(study, name)
    return {source: tables[source] for source in sorted(tables)}


def write_study_tables(study_folder: Path, target_dir: Path) -> list[Path]:
    """Write the data of a study of format 2 into the data folder of a model.

    Every table of `study_tables` becomes `target_dir/<study>/.<study>_<source>.tsv`,
    the file `PKPDSimulationExperiment.load_dataframe` reads; the files of an
    earlier collection of the study are replaced, so a source which the study
    no longer has leaves no file behind.

    Args:
        study_folder: the folder of a study of format 2.
        target_dir: the data folder of the model.

    Returns:
        The written files.
    """
    name = study_folder.name
    tables = study_tables(study_folder)
    folder = target_dir / name
    if folder.is_dir():
        for old in folder.glob(f".{name}_*.tsv"):
            old.unlink()
    folder.mkdir(parents=True, exist_ok=True)
    paths = []
    for source, df in tables.items():
        path = folder / f".{name}_{source}.tsv"
        df.to_csv(path, sep="\t", index=False, na_rep="NA")
        paths.append(path)
    return paths


def collect_tsv_files(
    source_dir: Path,
    target_dir: Path,
    filter_study: Callable[[str], bool] | None = None,
) -> None:
    """Collect the data of the studies of a substance into the data folder of a model.

    Every study folder of `source_dir` which passes `filter_study` becomes a
    folder of the same name in `target_dir` with the data of the study as the
    tables of study format 1, and every written file is reported on the
    console. A study of format 2 is read with the library of PK-DB and
    written with `write_study_tables`; the hidden `.tsv` files of a study of
    format 1 are copied. The figures, the tables and the publication are left
    behind.

    Args:
        source_dir: the folder of the studies of a substance, e.g.
            `studies_dir() / "albuterol"`.
        target_dir: the data folder of the model, created if it does not exist
            (its parent has to).
        filter_study: decides by the name of a study folder whether the study
            is collected; all studies by default.

    Raises:
        ValueError: if a study has a format which is neither 1 nor 2, see
            `study_format`; nothing is written then.
    """
    if not filter_study:

        def filter_study(study_name: str) -> bool:
            return True

    study_folders = sorted(
        folder
        for folder in source_dir.iterdir()
        if folder.is_dir()
        and not folder.name.startswith(".")
        and filter_study(folder.name)
    )
    formats = {folder.name: study_format(folder) for folder in study_folders}
    other = {name: value for name, value in formats.items() if value not in (1, 2)}
    if other:
        names = ", ".join(f"{name} (format {value})" for name, value in other.items())
        raise ValueError(
            f"The studies {names} of '{source_dir}' are neither in study format 1 "
            f"nor 2. Nothing is written, the data of '{target_dir}' are left as "
            f"they are."
        )

    target_dir.mkdir(exist_ok=True)
    for study_folder in study_folders:
        dest_study_folder = target_dir / study_folder.name
        if formats[study_folder.name] == 2:
            console.rule(f"{dest_study_folder}", style="white", align="left")
            for path in write_study_tables(study_folder, target_dir):
                console.print(f"'{study_folder}' -> '{path}'")
            continue

        tsv_files = sorted(study_folder.glob(".*.tsv"))
        if tsv_files:
            dest_study_folder.mkdir(exist_ok=True)
            console.rule(f"{dest_study_folder}", style="white", align="left")

            for tsv_file in tsv_files:
                dest_file = dest_study_folder / tsv_file.name
                shutil.copy2(tsv_file, dest_file)
                console.print(f"'{tsv_file}' -> '{dest_file}'")


@dataclass(frozen=True)
class DatabaseDifference:
    """A study whose version in the database differs from its study folder.

    Attributes:
        study: the folder of the study.
        sid: the study id of PK-DB, `<substance>/<name>`.
        reason: why they differ, e.g. that the study is not in the database or
            that the digests of the source files differ.
    """

    study: Path
    sid: str
    reason: str


def database_differences(
    study_folders: Iterable[Path], client: Client | None = None
) -> list[DatabaseDifference]:
    """Compare studies of format 2 with their versions in a PK-DB database.

    A study is read with `prepare_study` and fetched from the server
    (`pkdb.Client().studies.get`); the two agree when the digests of their
    source files (`source_digest`) are equal, i.e. the study folder is the one
    that was uploaded. The server and the credentials are those of the
    environment of `pkdb.Client` (`PKDB_ENDPOINT`, `PKDB_API_KEY`) unless a
    client is given.

    Args:
        study_folders: the folders of the studies.
        client: the client of the server, `pkdb.Client()` by default.

    Returns:
        The studies which differ, in the order of `study_folders`.

    Raises:
        pkdb.errors.ClientError: if the server cannot be asked, e.g. without
            the credentials of a private study; a study the server does not
            have (404) is a difference.
    """
    from pkdb import Client
    from pkdb.errors import ClientError

    client = client if client is not None else Client()
    differences = []
    for folder in study_folders:
        local = prepare_study(folder)
        sid = f"{folder.parent.name}/{folder.name}"
        try:
            remote = client.studies.get(sid)
        except ClientError as err:
            if err.status_code != 404:
                raise
            differences.append(DatabaseDifference(folder, sid, "not in the database"))
            continue
        if remote.source_digest != local.source_digest:
            differences.append(
                DatabaseDifference(
                    folder,
                    sid,
                    f"the database has another version (source digest "
                    f"{remote.source_digest[:12]}, the folder {local.source_digest[:12]})",
                )
            )
    return differences
