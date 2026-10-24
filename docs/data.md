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

The experiments read the tables of study format 1, one hidden `.tsv` file per source of a study (`.Boulton1996a_Fig1.tsv`). pkdb_data is curated in study format 2 (`"format": 2` in the `study.json`, with the tables `subjects.tsv`, `interventions.tsv`, `timecourses_<source>.tsv`, `outputs_<source>.tsv`, ...). `collect_tsv_files(source_dir, target_dir, filter_study=None)` writes the data of every study folder of `source_dir` into a folder of the same name in `target_dir`, as the tables of format 1, so that the experiments read them unchanged:

- a study of format 2 is read with the library of PK-DB, `pkdb.prepare`, with the vocabulary of the checkout (`vocabulary.lock.json`); this is the reading and the validation of `pkdb check` and of the upload to the database, so a study which `pkdb check` refuses is an error. Its data become one table per source (the reported outputs and timecourses of `outputs_<source>.tsv` and `timecourses_<source>.tsv`, the statistics PK-DB calculates or normalizes from them are left out) one table per source of the scatters of `scatters_<source>.tsv` (a point per row with the columns `x_...` and `y_...` of format 1, the name of the scatter as `label`, `<source>Scatter` when the source also has measurements), and the tables `TabGroups`, `TabIndividuals` and `TabInterventions` of the subjects, their characteristics and the interventions (`study_tables`, `write_study_tables`). A digitized error bar gives the statistic it names (`sd`, `se` or `gsd`), the missing `sd`, `se` and `cv` of a row are calculated as PK-DB does (`cv` is a fraction of the mean), an output (which has no label in format 2) gets the label of its measurement, substance, tissue, interventions and subject (`cumulative_amount_change_chloride_urine_HCTZ25_sus_all`), a point of a scatter has its values also in `x` and `y` and counts 1, and the earlier tables of the study in `target_dir` are replaced;
- the hidden `.tsv` files of a study of format 1 are copied; a study without them is left out;
- the figures, the tables and the publication stay in pkdb_data, and `filter_study` selects studies by the name of their folder.

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

## Comparing with the database

The models read the study folders of the checkout, which are the studies of a revision of pkdb_data. `database_differences(study_folders, client=None)` compares them with the studies uploaded to a PK-DB server: a study agrees when the digest of its source files (`source_digest` of `pkdb.prepare`) is the one of the server, a study the server does not have and a study with another digest are a `DatabaseDifference` with the reason. The server and the credentials are those of `pkdb.Client` (`PKDB_ENDPOINT`, `PKDB_API_KEY`); a private study needs an API key.

```python
# not executed: needs a checkout of pkdb_data and a PK-DB server
from pkpdsim.data import database_differences, studies_dir

folders = sorted((studies_dir() / "albuterol").iterdir())
for difference in database_differences(folders):
    print(difference.sid, difference.reason)
```
