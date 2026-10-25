# Study data

The studies against which the models are fitted and validated are curated in [PK-DB](https://pk-db.com), in a checkout of [pkdb_data](https://github.com/matthiaskoenig/pkdb_data): its folder `studies` holds a folder per substance and in it a folder per study, with the publication, its figures and tables and the digitized data as hidden `.tsv` files (`.Boulton1996a_Fig1.tsv`). A model keeps the studies it uses in its own `data/` folder, so that its simulations run without the checkout; `pkpdsim.data` copies them there.

## Finding the studies

`studies_dir()` returns the folder `studies` of the checkout. The checkout is the path of the environment variable `PKDB_DATA`; without it the checkout is `../pkdb_data` relative to the working directory, i.e. a clone next to the repository a script runs in. A folder which does not exist is an error which names both:

```python
import os
from pathlib import Path

from pkpdsim.data import PKDB_DATA, studies_dir

# a checkout of pkdb_data with one study of albuterol
study = Path("checkouts/pkdb_data/studies/albuterol/Boulton1996a")
study.mkdir(parents=True)
(study / ".Boulton1996a_Fig1.tsv").write_text("time\tmean\n0.0\t0.0\n")
(study / "Boulton1996a_Fig1.png").write_bytes(b"")

os.environ[PKDB_DATA] = "checkouts/pkdb_data"
print(studies_dir().relative_to(Path.cwd()).as_posix())

os.environ[PKDB_DATA] = "missing/pkdb_data"
try:
    studies_dir()
except FileNotFoundError as error:
    print(type(error).__name__)
```

```text
checkouts/pkdb_data/studies
FileNotFoundError
```

## Collecting the data of a model

The experiments read the tables of study format 1, one hidden `.tsv` file per source of a study (`.Boulton1996a_Fig1.tsv`). pkdb_data is curated in study format 2 (`"format": 2` in the `study.json`, with the tables `subjects.tsv`, `interventions.tsv`, `timecourses_<source>.tsv`, `outputs_<source>.tsv`, ...). `collect_tsv_files(source_dir, target_dir, filter_study=None, *, studies=(), keep_files=False)` writes the data of every study folder of `source_dir` (and of the further study folders `studies`) into a folder of the same name in `target_dir`, as the tables of format 1, so that the experiments read them unchanged:

- a study of format 2 is read with the library of PK-DB, `pkdb.prepare`, with the vocabulary of the checkout (`vocabulary.lock.json`); this is the reading and the validation of `pkdb check` and of the upload to the database, so a study which `pkdb check` refuses is an error. Its data become one table per source (the reported outputs and timecourses of `outputs_<source>.tsv` and `timecourses_<source>.tsv`, the statistics PK-DB calculates or normalizes from them are left out) one table per source of the scatters of `scatters_<source>.tsv` (a point per row with the columns `x_...` and `y_...` of format 1, the name of the scatter as `label`, `<source>Scatter` when the source also has measurements), and the tables `TabGroups`, `TabIndividuals` and `TabInterventions` of the subjects, their characteristics and the interventions (`study_tables`, `write_study_tables`). A digitized error bar gives the statistic it names (`sd`, `se` or `gsd`), the missing `sd`, `se` and `cv` of a row are calculated as PK-DB does (`cv` is a fraction of the mean), an output (which has no label in format 2) gets the label of its measurement, substance, tissue, interventions and subject (`cumulative_amount_change_chloride_urine_HCTZ25_sus_all`), a point of a scatter has its values also in `x` and `y`, and the earlier tables of the study in `target_dir` are replaced; with `keep_files=True` the tables of the study folder in `target_dir` which the collection does not write stay, e.g. tables a model keeps beside the ones of PK-DB (see [the rows and columns](#the-rows-and-columns));
- the hidden `.tsv` files of a study of format 1 are copied; a study without them is left out;
- the figures, the tables and the publication stay in pkdb_data, and `filter_study` selects studies of `source_dir` by the name of their folder;
- `studies` takes the study folders of other substances which a model uses (`studies=[studies_dir() / "midazolam" / "Chung2006"]` beside the studies of its own substance); two studies of the same name are an error.

What the curation should look at is reported as a `CollectionWarning` and the collection goes on: a curated label which is no SId (`corr_creatine_cl-pravastatin_cl`), the same output in two sources of a study (the same measurement, substance, tissue, subject, interventions, time and value), a curator 0 in a curve of geometric means (see below).

Every written file is reported on the console. Continuing the snippet above:

```python
from pkpdsim.data import collect_tsv_files

os.environ[PKDB_DATA] = "checkouts/pkdb_data"
collect_tsv_files(source_dir=studies_dir() / "albuterol", target_dir=Path("data"))
print(sorted(path.as_posix() for path in Path("data").rglob("*")))
```

The console shows the study folder and every copy, from the absolute path of the source to the target, then the snippet prints the data folder:

```text
['data/Boulton1996a', 'data/Boulton1996a/.Boulton1996a_Fig1.tsv']
```

A model calls it from the script `data/collect_data.py` of its folder:

```python
# not executed: data/collect_data.py of a model, which needs a checkout of pkdb_data
from pkdb_models.models.albuterol import DATA_PATH_ALBUTEROL

from pkpdsim.data import collect_tsv_files, studies_dir

collect_tsv_files(
    source_dir=studies_dir() / "albuterol", target_dir=DATA_PATH_ALBUTEROL
)
```

## The rows and columns

The tables of a study of format 2 have the columns of format 1 (`MEASUREMENT_COLUMNS`, `SCATTER_COLUMNS`, `GROUP_COLUMNS`, `INDIVIDUAL_COLUMNS`, `INTERVENTION_COLUMNS`), so an experiment selects its rows as before. A missing value is `NA`. The ports of the models to format 2 showed what the experiments have to know about the rows:

- **Geometric means.** A row reported as a geometric mean has its value in `gmean` and none in `mean`, its spread in `gsd` (a factor, from a geometric error bar or calculated from the `gcv`) and `gcv` (a fraction), and the `calculation_type` `geometric mean` (PK-DB itself gives every measurement of a group without a calculation the calculation `sample mean`). An experiment reads such a curve from `gmean` (`yid="gmean"`). `PKPDSimulationExperiment.load_dataframe` gives the geometric means their unit (`add_gmean_unit`, the column `gmean_unit`), which sbmlsim does not read from `unit` for them. A curve of geometric means has no curator 0: the mean of 0 at the time of the dose (aliskiren/Tapaninen2010 Fig1, the urine at 0 hr of apixaban/Chang2016 Tab2) cannot be a geometric mean and is left out with a `CollectionWarning`. The groups have the `cv`, `gmean`, `gsd` and `gcv` of their characteristics as well (the geometric mean of the renin activity of captopril/MacFadyen1991).
- **Values below the limit of quantification.** A point below the LOQ is a row without a central value, with the LOQ in `max` and a comment: it has `NA` in `value`, `mean`, `median` and `gmean`, so an experiment drops it with `df.dropna(subset=["mean"])` (`"value"` for individuals, `"gmean"` for geometric means). It has the `calculation_type` of the other rows of its curve, a point below the LOQ of a curve of geometric means is a `geometric mean`.
- **Errors.** `cv` and `gcv` are fractions of the mean, not percent as in format 2 and in some tables of format 1 (an experiment which divides them by 100 is wrong by that factor). The missing `sd`, `se` and `cv` of a row are calculated from the others and the count. A table always has the columns of the errors, a curve without any error has `NA` in all of them, and `DataSet.from_df` of sbmlsim renames `sd` to `mean_sd`: whether a curve has errors is `df["sd"].notna().any()` (`dset["mean_sd"].notna().any()`), not whether the column exists. Format 2 has the median with `min` and `max` but no quartiles or confidence intervals, so neither do the tables.
- **Counts.** `count` (and `group_count` of the groups, `doses` of the interventions) are integers, a missing count is `NA` (`53`, not `53.0`). A row of a group without one has the count of the group, a row of an individual counts 1.
- **Labels.** A timecourse has the label of its curation, an output the label generated from its measurement, substance, tissue, interventions and subject. A generated label which another source of the study has as well (the trough of visit 8 of sorafenib/Fucile2015 in `outputs_Fig1.tsv` and the troughs of the same patient in `outputs_Tab3.tsv`) gets its source appended in every source which has it (`concentration_sorafenib_plasma_SOF400_1_Fig1`, `..._Tab3`), so a label is unique across the sources of a study and the datasets of an experiment can be keyed by it. The labels which only one source has stay as they are; a source added to a study can change the label of another source. A generated label is an SId (a leading digit gets a `_`); a curated label which is no SId (`montelukast-M1_het`) stays as curated with a `CollectionWarning`, an experiment keys its datasets by `to_sid(label)`.
- **Interventions.** `intervention` holds the names of the `interventions` cell of format 2, comma separated, in the natural order which `pkdb format` gives the cell (`GEM,SIM`, `SEM0.5,SEM0.25,SEM1`); a generated label joins them with `_` (`..._GEM_SIM_...`). A selection of format 1 in another order (`SIM,GEM`) selects nothing.
- **Scatters.** A point of a scatter has the measurement, substance, tissue, method, value, unit, time and interventions of each axis in `x_...` and `y_...` (`x_time`, `x_intervention`); `time` is the one of `y` (else of `x`) and `intervention` the interventions of both axes, as before. A measurement without a substance (the eGFR of empagliflozin/Sarashina2014 Fig2) has `NA` in `x_substance`, so a selection by the substance of `x` drops it. A point of an individual counts 1, a point of a group the count of the group (the means of the volunteers of inulin/Prescott1991 Fig4); several points of a scatter with the same group, interventions and times are anonymous subjects of the group and count 1 (Sarashina2014 Fig2).
- **Individuals.** `TabIndividuals` is long, a row per characteristic of an individual (`measurement_type`, `choice`, `value`, `unit`), and `mean_of_individuals` calculates the mean, SD and count of the individuals of a table at every time, e.g. for a group whose means the publication does not report. A table with a column per characteristic is a pivot:

```python
import pandas as pd

individuals = pd.DataFrame(
    {
        "name": ["S1", "S1", "S2", "S2"],
        "measurement_type": ["age", "sex", "age", "sex"],
        "choice": [None, "M", None, "F"],
        "value": [24.0, None, 31.0, None],
    }
)
wide = individuals.pivot(index="name", columns="measurement_type", values="value")
wide["sex"] = individuals.pivot(
    index="name", columns="measurement_type", values="choice"
)["sex"]
print(wide.reset_index().to_string(index=False))
```

```text
name  age sex
  S1 24.0   M
  S2 31.0   F
```

A model whose digitized times scatter around the sampling times of the protocol moves them there with `sampling_times(df, times, start=None, by="label")`: within a curve the points keep their order and go to distinct sampling times, the assignment with the least total shift.

What format 2 does not store, the tables do not have either. The tables have no collection intervals: format 2 has only the `time` of a measurement (for a urine collection the end of the interval), PK-DB has no start and end of a collection yet ([pkdb#1032](https://github.com/matthiaskoenig/pkdb/issues/1032)), so the experiments which read `tstart` and `tend` (or `time_start` and `time_end`) of format 1 take the intervals from the publication until then. Format 2 has no digitized actual times beside the nominal ones, no SD of a single subject, no quartiles of a box plot and no confidence intervals either ([pkdb#1044](https://github.com/matthiaskoenig/pkdb/issues/1044)). Neither do the tables have the columns which the curators of format 1 added to a single table (`mean_se`, `dose`, `wlow`, `wup`, ...): an experiment which drops or reads such a column changes with the collection from format 2.

## Comparing with the database

The models read the study folders of the checkout, which are the studies of a revision of pkdb_data. `database_differences(study_folders, client=None)` compares them with the studies uploaded to a PK-DB server: a study agrees when the digest of its source files (`source_digest` of `pkdb.prepare`) is the one of the server, a study the server does not have and a study with another digest are a `DatabaseDifference` with the reason. The server and the credentials are those of `pkdb.Client` (`PKDB_ENDPOINT`, `PKDB_API_KEY`); a private study needs an API key.

```python
# not executed: needs a checkout of pkdb_data and a PK-DB server
from pkpdsim.data import database_differences, studies_dir

folders = sorted((studies_dir() / "albuterol").iterdir())
for difference in database_differences(folders):
    print(difference.sid, difference.reason)
```
