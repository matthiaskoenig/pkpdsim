"""The study data of a model, collected from the study folders of PK-DB.

The studies of [PK-DB](https://pk-db.com) are curated in a checkout of
[pkdb_data](https://github.com/matthiaskoenig/pkdb_data): the folder `studies`
holds a folder per substance and in it a folder per study, with the digitized
data as hidden `.tsv` files next to the figures and tables they were read
from. A model keeps the subset of the studies it uses in its own `data/`
folder, so its simulations run without the checkout; `collect_tsv_files`
copies them there from the folder `studies_dir()` finds.

The hidden `.tsv` files are the data of study format 1. pkdb_data is curated
in study format 2 (`"format": 2` of the `study.json`, the tables
`timecourses_<source>.tsv`, `outputs_<source>.tsv`, ...), which the models do
not read yet; `collect_tsv_files` refuses a study of format 2 instead of
copying nothing (pkdb_models#109).
"""

import json
import os
import shutil
from collections.abc import Callable
from pathlib import Path

from sbmlutils.console import console

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


def collect_tsv_files(
    source_dir: Path,
    target_dir: Path,
    filter_study: Callable[[str], bool] | None = None,
) -> None:
    """Copy the data of the studies of a substance into the data folder of a model.

    The data files of a study folder of pkdb_data are its hidden `.tsv` files
    (`.Boulton1996a_Fig1.tsv`); the figures, the tables and the publication
    are left behind. Every study folder of `source_dir` which passes
    `filter_study` and has data files becomes a folder of the same name in
    `target_dir` with copies of them, and every copy is reported on the
    console.

    The hidden `.tsv` files are the data of study format 1. A study of format 2
    has none, so its collection would copy nothing and leave the data of the
    model as they were; the studies are therefore checked first, and nothing is
    copied if one of them is not in format 1.

    Args:
        source_dir: the folder of the studies of a substance, e.g.
            `studies_dir() / "albuterol"`.
        target_dir: the data folder of the model, created if it does not exist
            (its parent has to).
        filter_study: decides by the name of a study folder whether the study
            is collected; all studies by default.

    Raises:
        ValueError: if a study of `source_dir` which passes `filter_study` is
            not in study format 1, see `study_format`.
    """
    if not filter_study:

        def filter_study(study_name: str) -> bool:
            return True

    study_folders = sorted(
        folder
        for folder in source_dir.iterdir()
        if folder.is_dir() and filter_study(folder.name)
    )
    formats = {folder.name: study_format(folder) for folder in study_folders}
    other = {name: value for name, value in formats.items() if value != 1}
    if other:
        names = ", ".join(f"{name} (format {value})" for name, value in other.items())
        raise ValueError(
            f"The studies {names} of '{source_dir}' are not in study format 1: "
            f"collect_tsv_files copies the hidden .tsv files of format 1, which "
            f"a study of format 2 does not have. Nothing is copied, the data of "
            f"'{target_dir}' are left as they are."
        )

    target_dir.mkdir(exist_ok=True)
    for study_folder in study_folders:
        tsv_files = sorted(study_folder.glob(".*.tsv"))
        if tsv_files:
            dest_study_folder = target_dir / study_folder.name
            dest_study_folder.mkdir(exist_ok=True)
            console.rule(f"{dest_study_folder}", style="white", align="left")

            for tsv_file in tsv_files:
                dest_file = dest_study_folder / tsv_file.name
                shutil.copy2(tsv_file, dest_file)
                console.print(f"'{tsv_file}' -> '{dest_file}'")
