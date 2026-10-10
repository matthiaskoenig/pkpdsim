"""The study data of a model, collected from the study folders of PK-DB.

The studies of [PK-DB](https://pk-db.com) are curated in a checkout of
[pkdb_data](https://github.com/matthiaskoenig/pkdb_data): the folder `studies`
holds a folder per substance and in it a folder per study, with the digitized
data as hidden `.tsv` files next to the figures and tables they were read
from. A model keeps the subset of the studies it uses in its own `data/`
folder, so its simulations run without the checkout; `collect_tsv_files`
copies them there from the folder `studies_dir()` finds.
"""

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

    Args:
        source_dir: the folder of the studies of a substance, e.g.
            `studies_dir() / "albuterol"`.
        target_dir: the data folder of the model, created if it does not exist
            (its parent has to).
        filter_study: decides by the name of a study folder whether the study
            is collected; all studies by default.
    """
    if not filter_study:

        def filter_study(study_name: str) -> bool:
            return True

    target_dir.mkdir(exist_ok=True)
    for study_folder in source_dir.iterdir():
        if study_folder.is_dir() and filter_study(study_folder.name):
            tsv_files = list(study_folder.glob(".*.tsv"))

            if tsv_files:
                dest_study_folder = target_dir / study_folder.name
                dest_study_folder.mkdir(exist_ok=True)
                console.rule(f"{dest_study_folder}", style="white", align="left")

                for tsv_file in tsv_files:
                    dest_file = dest_study_folder / tsv_file.name
                    shutil.copy2(tsv_file, dest_file)
                    console.print(f"'{tsv_file}' -> '{dest_file}'")

    return
