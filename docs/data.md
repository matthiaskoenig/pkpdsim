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

The data files are the hidden `.tsv` files of study format 1. pkdb_data is curated in study format 2 (`"format": 2` in the `study.json`, with the tables `timecourses_<source>.tsv`, `outputs_<source>.tsv`, `subjects.tsv`, ...), which the models do not read yet: `collect_tsv_files` raises a `ValueError` which names the studies of another format (`study_format`) and copies nothing, instead of copying nothing silently and leaving the data of the model as they were.

`collect_tsv_files(source_dir, target_dir, filter_study=None)` copies the data files of every study folder of `source_dir` into a folder of the same name in `target_dir`; a study without data files is left out, as are the figures, the tables and the publication, and `filter_study` selects studies by the name of their folder. Every copy is reported on the console. Continuing the snippet above:

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
