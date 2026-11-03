"""Tests of the credibility test model and its fit."""

import pandas as pd

from pkpdsim.credibility import CouLevel
from pkpdsim.model import PKPDModel


def test_load_credibility(credibility_model: PKPDModel) -> None:
    assert credibility_model.load_credibility().context("dose").level is CouLevel.C2


def test_the_fit_has_every_kind(credibility_model: PKPDModel) -> None:
    (fit,) = credibility_model.fits
    table = pd.read_csv(
        credibility_model.fits_dir / fit / "metrics_mappings.tsv", sep="\t"
    )
    kinds = {
        (e, m): k
        for e, m, k in zip(table.experiment, table.mapping, table.kind, strict=True)
    }
    assert kinds == {
        ("CredA", "fm_train"): "training",
        ("CredA", "fm_dose"): "validation",
        ("CredA", "fm_renal"): "validation",
        ("CredA", "fm_limit"): "validation",
        ("CredA", "fm_subject"): "validation",
        ("CredA", "fm_outlier"): "outlier",
        ("CredB", "fm_dose"): "training",
        ("CredB", "fm_zero"): "validation",
    }
