"""Cross-check the reference physiology with the one of PK-Sim.

Prints the relative difference of the organ volumes and blood flows of the
reference physiology (`pkpdsim/pbpk/data/physiology_human.toml`, ICRP 89) to
those of the PK-Sim database for the population `European_ICRP_2002` at age 30.

PK-Sim (GPLv2) is a check only: no number of the database enters pkpdsim. The
database is not downloaded or bundled; the developer passes the path of a local
`PKSimDB.sqlite` (https://github.com/Open-Systems-Pharmacology/PK-Sim, folder
`src/Db`). It is opened read-only.

```bash
uv run python scripts/pbpk_crosscheck_pksim.py /path/to/PKSimDB.sqlite
```

The organ volume of the reference physiology is the tissue and its blood,
`mass/density + blood*V_blood`; the flow is `flow*cardiac_output`. Adipose is
compared with the ICRP 89 separable adipose mass (in the model adipose is the
remainder of the body weight). The flow of the liver is its arterial flow, as
in PK-Sim.
"""

import argparse
import itertools
import sqlite3
from pathlib import Path
from typing import Any

from pkpdsim.pbpk import Organ, load

POPULATION = "European_ICRP_2002"
AGE = 30.0

PKSIM_ORGANS: dict[str, Organ] = {
    "Lung": Organ.LUNG,
    "Heart": Organ.HEART,
    "Brain": Organ.BRAIN,
    "Liver": Organ.LIVER,
    "Kidney": Organ.KIDNEY,
    "Stomach": Organ.STOMACH,
    "SmallIntestine": Organ.SMALL_INTESTINE,
    "LargeIntestine": Organ.LARGE_INTESTINE,
    "Spleen": Organ.SPLEEN,
    "Pancreas": Organ.PANCREAS,
    "Muscle": Organ.MUSCLE,
    "Skin": Organ.SKIN,
    "Fat": Organ.ADIPOSE,
    "Bone": Organ.BONE,
    "Gonads": Organ.GONADS,
}
PKSIM_POOLS: dict[str, str] = {
    "ArterialBlood": "arterial",
    "VenousBlood": "venous",
    "PortalVein": "portal",
}


def _v(table: dict[str, Any], key: str) -> float:
    return float(table[key]["value"])


def pksim_values(db: Path, gender: str) -> dict[tuple[str, str], float]:
    """Read the volumes and specific blood flows of PK-Sim (l, l/min/l) at `AGE`.

    PK-Sim stores a curve over age per parameter and interpolates linearly
    between its ages, so does this function.
    """
    connection = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    try:
        rows = connection.execute(
            "SELECT container_name, parameter_name, age, mean_value "
            "FROM tab_container_parameter_curves "
            "WHERE population = ? AND gender = ? "
            "AND parameter_name IN ('Volume', 'Specific blood flow rate') "
            "ORDER BY container_name, parameter_name, age",
            (POPULATION, gender),
        ).fetchall()
    finally:
        connection.close()
    curves: dict[tuple[str, str], list[tuple[float, float]]] = {}
    for name, parameter, age, value in rows:
        curves.setdefault((name, parameter), []).append((age, value))
    return {key: _interpolate(curve, AGE) for key, curve in curves.items()}


def _interpolate(curve: list[tuple[float, float]], age: float) -> float:
    for (a0, v0), (a1, v1) in itertools.pairwise(curve):
        if a0 <= age <= a1:
            return v0 + (v1 - v0) * (age - a0) / (a1 - a0)
    raise ValueError(f"age {age} outside of the curve {curve}")


def our_values(sex: str) -> dict[str, tuple[float, float | None]]:
    """The volume (l) and flow (l/min) of every organ and pool."""
    p = load("physiology_human")[sex]
    blood_volume = _v(p["blood"], "volume")
    values: dict[str, tuple[float, float | None]] = {}
    for name, organ in PKSIM_ORGANS.items():
        table = p["organs"][organ]
        density = _v(table, "density") if "density" in table else 1.0
        volume = _v(table, "mass") / density + _v(table, "blood") * blood_volume
        flow = None
        if organ is not Organ.LUNG:
            flow = _v(table, "flow") * _v(p, "cardiac_output")
        values[name] = (volume, flow)
    for name, pool in PKSIM_POOLS.items():
        values[name] = (_v(p["blood"], pool) * blood_volume, None)
    return values


def _rel(ours: float, theirs: float) -> str:
    return f"{(ours - theirs) / theirs:+8.1%}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("db", type=Path, help="path of a local PKSimDB.sqlite")
    args = parser.parse_args()
    if not args.db.is_file():
        parser.error(f"no such file: {args.db}")

    for sex, gender in (("male", "MALE"), ("female", "FEMALE")):
        pksim = pksim_values(args.db, gender)
        print(f"\n{sex} ({POPULATION}, age {AGE:g}); relative difference ours/PK-Sim")
        print(
            f"{'organ':<16}{'V ours':>9}{'V PK-Sim':>10}{'dV':>9}"
            f"{'Q ours':>9}{'Q PK-Sim':>10}{'dQ':>9}"
        )
        for name, (volume, flow) in our_values(sex).items():
            their_volume = pksim[(name, "Volume")]
            line = f"{name:<16}{volume:9.4f}{their_volume:10.4f}"
            line += _rel(volume, their_volume)
            specific = pksim.get((name, "Specific blood flow rate"))
            if flow is not None and specific is not None:
                their_flow = specific * their_volume
                line += f"{flow:9.4f}{their_flow:10.4f}{_rel(flow, their_flow)}"
            print(line)


if __name__ == "__main__":
    main()
