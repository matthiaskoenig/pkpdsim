"""Fit the waist circumference of NHANES 2013-2018 on the covariates of a build.

The waist circumference is the covariate of the visceral depot of the body
composition (`pkpdsim.pbpk.adipose`), and a model which does not measure it
predicts it from the covariates it has. This script fits

    WC [cm] = b0 + b_bw * BW [kg] + b_height * HEIGHT [cm] + b_age * AGE [yr]

per sex by weighted least squares on the public files of the examinations
2013-2014, 2015-2016 and 2017-2018 (CDC/NCHS, public domain): demographics
(`DEMO`) and body measures (`BMX`, the waist `BMXWAIST` measured at the
uppermost lateral border of the ilium), all races and Hispanic origins, the
women pregnant at the examination left out, weighted with the 6-year
examination weights. The respondents of 80 years and over are left out
(`RIDAGEYR` is top-coded at 80, so their age is not known) and so are the
subjects below `ADULT`, whose waist the regression of the adults does not
describe.

The intercept is then re-anchored on the cohort of the sensitivity of the
visceral depot (`NAKO`, Haueise 2023): the shift is the difference between the
prediction at the mean age, height and body weight of that cohort and its mean
waist, which absorbs the measurement site (NHANES measures at the iliac crest,
the German National Cohort at the midpoint between the lowest rib and the
iliac crest, the site of the WHO protocol) together with the difference of the
populations, which these data cannot separate.

    uv run python scripts/database/nhanes_waist.py [--cache DIR]

prints the coefficients, the fit and the checks as the entries of the
`[waist.<sex>]` sections of `tables/adipose_human.toml`, which hold them (the
table is written by hand, as the digitized figures are).
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from nhanes_populations import SEXES, TOP_CODED, fetch, fmt

ADULT = 20.0
"""The age (years) from which a respondent is in the fit."""

COLUMNS = ("BMXWT", "BMXHT", "age")
"""The predictors: body weight (kg), height (cm) and the age (years)."""

KEYS = ("intercept", "bw", "height", "age")
"""The keys of the coefficients in the table."""

NAKO: dict[str, dict[str, float]] = {
    "male": {"n": 5697, "age": 52.1, "height": 178.2, "bw": 86.9, "waist": 97.4},
    "female": {"n": 5420, "age": 51.7, "height": 164.9, "bw": 71.2, "waist": 86.1},
}
"""The cohort of the sensitivity of the visceral depot, on whose level the
prediction is anchored: the means of Table 2 of Haueise 2023 (the German
National Cohort, whole-trunk MRI of 11,141 participants of 20 to 72 years)."""

WHITE = 3
"""`RIDRETH1` of the non-Hispanic white respondents, the check of the slopes."""


def fit(data: pd.DataFrame) -> tuple[np.ndarray, float, float]:
    """The weighted least squares of the waist on the covariates.

    Args:
        data: the respondents with the predictors, the waist `BMXWAIST` and
            the examination weight `w`.

    Returns:
        The coefficients in the order of `KEYS`, the weighted R2 and the
        weighted root mean squared residual (cm).
    """
    x = np.column_stack(
        [np.ones(len(data))] + [data[c].to_numpy(dtype=float) for c in COLUMNS]
    )
    y = data["BMXWAIST"].to_numpy(dtype=float)
    w = data["w"].to_numpy(dtype=float)
    root = np.sqrt(w)
    beta, *_ = np.linalg.lstsq(x * root[:, None], y * root, rcond=None)
    residual = y - x @ beta
    sse = float(np.average(residual**2, weights=w))
    sst = float(np.average((y - np.average(y, weights=w)) ** 2, weights=w))
    return beta, 1.0 - sse / sst, sse**0.5


def main() -> None:
    """Fit the waist of every sex and print the entries of the table."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--cache", type=Path, default=Path.home() / ".cache" / "pkpdsim" / "nhanes"
    )
    args = parser.parse_args()
    data = fetch(args.cache, ("BMX",))
    data = data[
        data["BMXWAIST"].notna()
        & data["BMXWT"].notna()
        & data["BMXHT"].notna()
        & (data["RIDAGEYR"] >= ADULT)
        & (data["RIDAGEYR"] < TOP_CODED)
    ]
    for sex, code in SEXES.items():
        part = data[data["RIAGENDR"] == code]
        beta, r2, rmse = fit(part)
        cohort = NAKO[sex]
        predicted = float(
            beta[0]
            + beta[1] * cohort["bw"]
            + beta[2] * cohort["height"]
            + beta[3] * cohort["age"]
        )
        shift = cohort["waist"] - predicted
        white, r2_white, _ = fit(part[part["RIDRETH1"] == WHITE])
        means = {c: float(np.average(part[c], weights=part["w"])) for c in COLUMNS}
        print(f"\n# [waist.{sex}]")
        print(f"# n = {len(part)}, R2 = {r2:.3f}, residual SD = {rmse:.2f} cm")
        print(
            "# weighted means: "
            + ", ".join(f"{c} {v:.1f}" for c, v in means.items())
            + f", waist {np.average(part['BMXWAIST'], weights=part['w']):.1f} cm"
        )
        print(
            f"# the prediction at the {sex} mean of the anchor cohort is "
            f"{predicted:.2f} cm against its {cohort['waist']:.1f} cm, "
            f"the shift {shift:+.2f} cm"
        )
        print(
            "# non-Hispanic white respondents only (the check of the slopes): "
            + ", ".join(f"{k} {fmt(b, 5)}" for k, b in zip(KEYS, white, strict=True))
            + f", R2 = {r2_white:.3f}"
        )
        values = (beta[0] + shift, *beta[1:])
        for key, b in zip(KEYS, values, strict=True):
            print(f"{key} = {fmt(b, 5)}")


if __name__ == "__main__":
    main()
