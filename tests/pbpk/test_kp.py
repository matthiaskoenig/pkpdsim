"""Tests of the Kp method of Rodgers & Rowland and its build."""

import copy
import dataclasses
from pathlib import Path

import pytest

from pkpdsim.pbpk import (
    Binding,
    Erythrocytes,
    Excretion,
    KpMethod,
    Mode,
    Organ,
    Route,
    Substance,
    U,
    Value,
    build,
    ids,
    load,
)
from pkpdsim.pbpk import kp as kp_module
from pkpdsim.pbpk.kp import (
    Composition,
    composition,
    kp_sources,
    kp_values,
    rodgers_rowland,
)
from pkpdsim.pbpk.physiology import Sex, reference_hematocrit

from .conftest import definition, dose_amount, individual, simple_drug, total_amount
from .conftest import load as load_model

DAY = 24 * 60.0


def _method(**changes: object) -> KpMethod:
    """A neutral compound of logP 2, merged tissue."""
    fields: dict[str, object] = {
        "method": "rodgers_rowland",
        "logp": 2.0,
        "pka": (),
        "compound": "neutral",
        "ftissue": 1.0,
    }
    fields.update(changes)
    return KpMethod(**fields)  # ty: ignore[invalid-argument-type]


def _drug(method: KpMethod, **changes: object) -> Substance:
    """The simple drug (fu 0.1) with a Kp method."""
    return dataclasses.replace(simple_drug(), distribution=method, **changes)


# Coutinho 2024 (Pharm Res 41:1121, open access, CC BY 4.0), Supplementary
# Tables S5 to S7: Kpu of the Rodgers & Rowland equation (their Eq. 3) for 13
# tissues and every logP column, computed with the rat composition of their
# Table III, which is the composition of data/tissue_composition.toml. Drug
# data of their Table II (weak bases, binding to albumin). The rows were
# transcribed by a script from the supplementary file (11095_2024_3703_MOESM1).
# Left out: adipose, whose two lipid fractions their Table III exchanges (the
# Kpu of adipose are computed with the exchanged values), and thymus, which is
# no organ of the model. Rodgers & Rowland 2005/2006 were not accessible, so
# their own examples are not reproduced (R44).
_PUBLISHED: list[tuple[str, float, float, float, Organ, float]] = [
    # drug, pKa, fup, logP, organ, Kpu
    # itraconazole (Table S5), pKa 4.57, fup 0.002
    ("itraconazole", 4.57, 0.002, 4.0, Organ.BONE, 226.2),
    ("itraconazole", 4.57, 0.002, 4.0, Organ.BRAIN, 418.3),
    ("itraconazole", 4.57, 0.002, 4.0, Organ.SMALL_INTESTINE, 486.9),
    ("itraconazole", 4.57, 0.002, 4.0, Organ.HEART, 241.3),
    ("itraconazole", 4.57, 0.002, 4.0, Organ.KIDNEY, 254.8),
    ("itraconazole", 4.57, 0.002, 4.0, Organ.LIVER, 247.4),
    ("itraconazole", 4.57, 0.002, 4.0, Organ.LUNG, 352.2),
    ("itraconazole", 4.57, 0.002, 4.0, Organ.MUSCLE, 152.4),
    ("itraconazole", 4.57, 0.002, 4.0, Organ.PANCREAS, 458.3),
    ("itraconazole", 4.57, 0.002, 4.0, Organ.SKIN, 746.7),
    ("itraconazole", 4.57, 0.002, 4.0, Organ.SPLEEN, 149.3),
    ("itraconazole", 4.57, 0.002, 4.5, Organ.BONE, 606.4),
    ("itraconazole", 4.57, 0.002, 4.5, Organ.BRAIN, 1269.5),
    ("itraconazole", 4.57, 0.002, 4.5, Organ.SMALL_INTESTINE, 1367.7),
    ("itraconazole", 4.57, 0.002, 4.5, Organ.HEART, 592.3),
    ("itraconazole", 4.57, 0.002, 4.5, Organ.KIDNEY, 663.9),
    ("itraconazole", 4.57, 0.002, 4.5, Organ.LIVER, 688.0),
    ("itraconazole", 4.57, 0.002, 4.5, Organ.LUNG, 883.8),
    ("itraconazole", 4.57, 0.002, 4.5, Organ.MUSCLE, 411.2),
    ("itraconazole", 4.57, 0.002, 4.5, Organ.PANCREAS, 1383.3),
    ("itraconazole", 4.57, 0.002, 4.5, Organ.SKIN, 2061.1),
    ("itraconazole", 4.57, 0.002, 4.5, Organ.SPLEEN, 366.3),
    ("itraconazole", 4.57, 0.002, 5.0, Organ.BONE, 1808.8),
    ("itraconazole", 4.57, 0.002, 5.0, Organ.BRAIN, 3961.0),
    ("itraconazole", 4.57, 0.002, 5.0, Organ.SMALL_INTESTINE, 4153.1),
    ("itraconazole", 4.57, 0.002, 5.0, Organ.HEART, 1702.3),
    ("itraconazole", 4.57, 0.002, 5.0, Organ.KIDNEY, 1957.8),
    ("itraconazole", 4.57, 0.002, 5.0, Organ.LIVER, 2081.4),
    ("itraconazole", 4.57, 0.002, 5.0, Organ.LUNG, 2564.7),
    ("itraconazole", 4.57, 0.002, 5.0, Organ.MUSCLE, 1229.7),
    ("itraconazole", 4.57, 0.002, 5.0, Organ.PANCREAS, 4308.1),
    ("itraconazole", 4.57, 0.002, 5.0, Organ.SKIN, 6217.4),
    ("itraconazole", 4.57, 0.002, 5.0, Organ.SPLEEN, 1052.4),
    ("itraconazole", 4.57, 0.002, 5.5, Organ.BONE, 5611.2),
    ("itraconazole", 4.57, 0.002, 5.5, Organ.BRAIN, 12472.3),
    ("itraconazole", 4.57, 0.002, 5.5, Organ.SMALL_INTESTINE, 12961.0),
    ("itraconazole", 4.57, 0.002, 5.5, Organ.HEART, 5212.5),
    ("itraconazole", 4.57, 0.002, 5.5, Organ.KIDNEY, 6049.3),
    ("itraconazole", 4.57, 0.002, 5.5, Organ.LIVER, 6487.8),
    ("itraconazole", 4.57, 0.002, 5.5, Organ.LUNG, 7880.3),
    ("itraconazole", 4.57, 0.002, 5.5, Organ.MUSCLE, 3818.0),
    ("itraconazole", 4.57, 0.002, 5.5, Organ.PANCREAS, 13557.3),
    ("itraconazole", 4.57, 0.002, 5.5, Organ.SKIN, 19360.8),
    ("itraconazole", 4.57, 0.002, 5.5, Organ.SPLEEN, 3222.1),
    ("itraconazole", 4.57, 0.002, 6.0, Organ.BONE, 17635.3),
    ("itraconazole", 4.57, 0.002, 6.0, Organ.BRAIN, 39387.4),
    ("itraconazole", 4.57, 0.002, 6.0, Organ.SMALL_INTESTINE, 40814.3),
    ("itraconazole", 4.57, 0.002, 6.0, Organ.HEART, 16312.6),
    ("itraconazole", 4.57, 0.002, 6.0, Organ.KIDNEY, 18987.9),
    ("itraconazole", 4.57, 0.002, 6.0, Organ.LIVER, 20422.1),
    ("itraconazole", 4.57, 0.002, 6.0, Organ.LUNG, 24689.7),
    ("itraconazole", 4.57, 0.002, 6.0, Organ.MUSCLE, 12002.8),
    ("itraconazole", 4.57, 0.002, 6.0, Organ.PANCREAS, 42805.9),
    ("itraconazole", 4.57, 0.002, 6.0, Organ.SKIN, 60923.8),
    ("itraconazole", 4.57, 0.002, 6.0, Organ.SPLEEN, 10083.2),
    ("itraconazole", 4.57, 0.002, 6.5, Organ.BONE, 55658.9),
    ("itraconazole", 4.57, 0.002, 6.5, Organ.BRAIN, 124500.4),
    ("itraconazole", 4.57, 0.002, 6.5, Organ.SMALL_INTESTINE, 128894.1),
    ("itraconazole", 4.57, 0.002, 6.5, Organ.HEART, 51414.3),
    ("itraconazole", 4.57, 0.002, 6.5, Organ.KIDNEY, 59903.2),
    ("itraconazole", 4.57, 0.002, 6.5, Organ.LIVER, 64486.1),
    ("itraconazole", 4.57, 0.002, 6.5, Organ.LUNG, 77845.8),
    ("itraconazole", 4.57, 0.002, 6.5, Organ.MUSCLE, 37885.6),
    ("itraconazole", 4.57, 0.002, 6.5, Organ.PANCREAS, 135298.1),
    ("itraconazole", 4.57, 0.002, 6.5, Organ.SKIN, 192357.8),
    ("itraconazole", 4.57, 0.002, 6.5, Organ.SPLEEN, 31780.0),
    ("itraconazole", 4.57, 0.002, 7.0, Organ.BONE, 175900.2),
    ("itraconazole", 4.57, 0.002, 7.0, Organ.BRAIN, 393651.3),
    ("itraconazole", 4.57, 0.002, 7.0, Organ.SMALL_INTESTINE, 407426.9),
    ("itraconazole", 4.57, 0.002, 7.0, Organ.HEART, 162415.5),
    ("itraconazole", 4.57, 0.002, 7.0, Organ.KIDNEY, 189288.7),
    ("itraconazole", 4.57, 0.002, 7.0, Organ.LIVER, 203828.8),
    ("itraconazole", 4.57, 0.002, 7.0, Organ.LUNG, 245940.0),
    ("itraconazole", 4.57, 0.002, 7.0, Organ.MUSCLE, 119734.0),
    ("itraconazole", 4.57, 0.002, 7.0, Organ.PANCREAS, 427783.9),
    ("itraconazole", 4.57, 0.002, 7.0, Organ.SKIN, 607988.5),
    ("itraconazole", 4.57, 0.002, 7.0, Organ.SPLEEN, 100391.3),
    # posaconazole (Table S6), pKa 4.67, fup 0.02
    ("posaconazole", 4.67, 0.02, 4.0, Organ.BONE, 181.1),
    ("posaconazole", 4.67, 0.02, 4.0, Organ.BRAIN, 396.6),
    ("posaconazole", 4.67, 0.02, 4.0, Organ.SMALL_INTESTINE, 415.7),
    ("posaconazole", 4.67, 0.02, 4.0, Organ.HEART, 170.5),
    ("posaconazole", 4.67, 0.02, 4.0, Organ.KIDNEY, 196.2),
    ("posaconazole", 4.67, 0.02, 4.0, Organ.LIVER, 208.6),
    ("posaconazole", 4.67, 0.02, 4.0, Organ.LUNG, 256.7),
    ("posaconazole", 4.67, 0.02, 4.0, Organ.MUSCLE, 123.5),
    ("posaconazole", 4.67, 0.02, 4.0, Organ.PANCREAS, 431.2),
    ("posaconazole", 4.67, 0.02, 4.0, Organ.SKIN, 621.9),
    ("posaconazole", 4.67, 0.02, 4.0, Organ.SPLEEN, 105.6),
    ("posaconazole", 4.67, 0.02, 4.5, Organ.BONE, 561.2),
    ("posaconazole", 4.67, 0.02, 4.5, Organ.BRAIN, 1247.4),
    ("posaconazole", 4.67, 0.02, 4.5, Organ.SMALL_INTESTINE, 1296.1),
    ("posaconazole", 4.67, 0.02, 4.5, Organ.HEART, 521.4),
    ("posaconazole", 4.67, 0.02, 4.5, Organ.KIDNEY, 605.2),
    ("posaconazole", 4.67, 0.02, 4.5, Organ.LIVER, 649.1),
    ("posaconazole", 4.67, 0.02, 4.5, Organ.LUNG, 788.1),
    ("posaconazole", 4.67, 0.02, 4.5, Organ.MUSCLE, 382.3),
    ("posaconazole", 4.67, 0.02, 4.5, Organ.PANCREAS, 1355.7),
    ("posaconazole", 4.67, 0.02, 4.5, Organ.SKIN, 1935.7),
    ("posaconazole", 4.67, 0.02, 4.5, Organ.SPLEEN, 322.5),
    ("posaconazole", 4.67, 0.02, 5.0, Organ.BONE, 1763.1),
    ("posaconazole", 4.67, 0.02, 5.0, Organ.BRAIN, 3937.9),
    ("posaconazole", 4.67, 0.02, 5.0, Organ.SMALL_INTESTINE, 4080.4),
    ("posaconazole", 4.67, 0.02, 5.0, Organ.HEART, 1631.0),
    ("posaconazole", 4.67, 0.02, 5.0, Organ.KIDNEY, 1898.6),
    ("posaconazole", 4.67, 0.02, 5.0, Organ.LIVER, 2041.9),
    ("posaconazole", 4.67, 0.02, 5.0, Organ.LUNG, 2468.4),
    ("posaconazole", 4.67, 0.02, 5.0, Organ.MUSCLE, 1200.4),
    ("posaconazole", 4.67, 0.02, 5.0, Organ.PANCREAS, 4279.5),
    ("posaconazole", 4.67, 0.02, 5.0, Organ.SKIN, 6090.4),
    ("posaconazole", 4.67, 0.02, 5.0, Organ.SPLEEN, 1008.4),
    ("posaconazole", 4.67, 0.02, 5.5, Organ.BONE, 5564.1),
    ("posaconazole", 4.67, 0.02, 5.5, Organ.BRAIN, 12445.9),
    ("posaconazole", 4.67, 0.02, 5.5, Organ.SMALL_INTESTINE, 12885.0),
    ("posaconazole", 4.67, 0.02, 5.5, Organ.HEART, 5139.8),
    ("posaconazole", 4.67, 0.02, 5.5, Organ.KIDNEY, 5988.5),
    ("posaconazole", 4.67, 0.02, 5.5, Organ.LIVER, 6446.7),
    ("posaconazole", 4.67, 0.02, 5.5, Organ.LUNG, 7781.9),
    ("posaconazole", 4.67, 0.02, 5.5, Organ.MUSCLE, 3787.7),
    ("posaconazole", 4.67, 0.02, 5.5, Organ.PANCREAS, 13525.2),
    ("posaconazole", 4.67, 0.02, 5.5, Organ.SKIN, 19228.8),
    ("posaconazole", 4.67, 0.02, 5.5, Organ.SPLEEN, 3177.2),
    ("posaconazole", 4.67, 0.02, 6.0, Organ.BONE, 17583.6),
    ("posaconazole", 4.67, 0.02, 6.0, Organ.BRAIN, 39350.7),
    ("posaconazole", 4.67, 0.02, 6.0, Organ.SMALL_INTESTINE, 40727.6),
    ("posaconazole", 4.67, 0.02, 6.0, Organ.HEART, 16235.7),
    ("posaconazole", 4.67, 0.02, 6.0, Organ.KIDNEY, 18922.1),
    ("posaconazole", 4.67, 0.02, 6.0, Organ.LIVER, 20375.6),
    ("posaconazole", 4.67, 0.02, 6.0, Organ.LUNG, 24584.9),
    ("posaconazole", 4.67, 0.02, 6.0, Organ.MUSCLE, 11969.4),
    ("posaconazole", 4.67, 0.02, 6.0, Organ.PANCREAS, 42762.6),
    ("posaconazole", 4.67, 0.02, 6.0, Organ.SKIN, 60776.0),
    ("posaconazole", 4.67, 0.02, 6.0, Organ.SPLEEN, 10035.7),
    ("posaconazole", 4.67, 0.02, 6.5, Organ.BONE, 55592.7),
    ("posaconazole", 4.67, 0.02, 6.5, Organ.BRAIN, 124431.2),
    ("posaconazole", 4.67, 0.02, 6.5, Organ.SMALL_INTESTINE, 128773.8),
    ("posaconazole", 4.67, 0.02, 6.5, Organ.HEART, 51324.0),
    ("posaconazole", 4.67, 0.02, 6.5, Organ.KIDNEY, 59821.8),
    ("posaconazole", 4.67, 0.02, 6.5, Organ.LIVER, 64422.8),
    ("posaconazole", 4.67, 0.02, 6.5, Organ.LUNG, 77720.7),
    ("posaconazole", 4.67, 0.02, 6.5, Organ.MUSCLE, 37842.3),
    ("posaconazole", 4.67, 0.02, 6.5, Organ.PANCREAS, 135219.4),
    ("posaconazole", 4.67, 0.02, 6.5, Organ.SKIN, 192159.7),
    ("posaconazole", 4.67, 0.02, 6.5, Organ.SPLEEN, 31724.2),
    ("posaconazole", 4.67, 0.02, 7.0, Organ.BONE, 175788.0),
    ("posaconazole", 4.67, 0.02, 7.0, Organ.BRAIN, 393479.2),
    ("posaconazole", 4.67, 0.02, 7.0, Organ.SMALL_INTESTINE, 407200.1),
    ("posaconazole", 4.67, 0.02, 7.0, Organ.HEART, 162282.8),
    ("posaconazole", 4.67, 0.02, 7.0, Organ.KIDNEY, 189157.8),
    ("posaconazole", 4.67, 0.02, 7.0, Organ.LIVER, 203712.2),
    ("posaconazole", 4.67, 0.02, 7.0, Organ.LUNG, 245750.6),
    ("posaconazole", 4.67, 0.02, 7.0, Organ.MUSCLE, 119659.5),
    ("posaconazole", 4.67, 0.02, 7.0, Organ.PANCREAS, 427593.4),
    ("posaconazole", 4.67, 0.02, 7.0, Organ.SKIN, 607631.5),
    ("posaconazole", 4.67, 0.02, 7.0, Organ.SPLEEN, 100309.2),
    # isavuconazole (Table S7), pKa 3.28, fup 0.01
    ("isavuconazole", 3.28, 0.01, 3.0, Organ.BONE, 27.9),
    ("isavuconazole", 3.28, 0.01, 3.0, Organ.BRAIN, 44.9),
    ("isavuconazole", 3.28, 0.01, 3.0, Organ.SMALL_INTESTINE, 57.2),
    ("isavuconazole", 3.28, 0.01, 3.0, Organ.HEART, 32.4),
    ("isavuconazole", 3.28, 0.01, 3.0, Organ.KIDNEY, 32.5),
    ("isavuconazole", 3.28, 0.01, 3.0, Organ.LIVER, 29.6),
    ("isavuconazole", 3.28, 0.01, 3.0, Organ.LUNG, 46.2),
    ("isavuconazole", 3.28, 0.01, 3.0, Organ.MUSCLE, 19.1),
    ("isavuconazole", 3.28, 0.01, 3.0, Organ.PANCREAS, 49.4),
    ("isavuconazole", 3.28, 0.01, 3.0, Organ.SKIN, 89.0),
    ("isavuconazole", 3.28, 0.01, 3.0, Organ.SPLEEN, 20.2),
    ("isavuconazole", 3.28, 0.01, 3.5, Organ.BONE, 66.0),
    ("isavuconazole", 3.28, 0.01, 3.5, Organ.BRAIN, 130.2),
    ("isavuconazole", 3.28, 0.01, 3.5, Organ.SMALL_INTESTINE, 145.4),
    ("isavuconazole", 3.28, 0.01, 3.5, Organ.HEART, 67.5),
    ("isavuconazole", 3.28, 0.01, 3.5, Organ.KIDNEY, 73.5),
    ("isavuconazole", 3.28, 0.01, 3.5, Organ.LIVER, 73.7),
    ("isavuconazole", 3.28, 0.01, 3.5, Organ.LUNG, 99.4),
    ("isavuconazole", 3.28, 0.01, 3.5, Organ.MUSCLE, 45.0),
    ("isavuconazole", 3.28, 0.01, 3.5, Organ.PANCREAS, 142.0),
    ("isavuconazole", 3.28, 0.01, 3.5, Organ.SKIN, 220.6),
    ("isavuconazole", 3.28, 0.01, 3.5, Organ.SPLEEN, 41.9),
    ("isavuconazole", 3.28, 0.01, 4.0, Organ.BONE, 186.4),
    ("isavuconazole", 3.28, 0.01, 4.0, Organ.BRAIN, 399.7),
    ("isavuconazole", 3.28, 0.01, 4.0, Organ.SMALL_INTESTINE, 424.3),
    ("isavuconazole", 3.28, 0.01, 4.0, Organ.HEART, 178.7),
    ("isavuconazole", 3.28, 0.01, 4.0, Organ.KIDNEY, 203.0),
    ("isavuconazole", 3.28, 0.01, 4.0, Organ.LIVER, 213.2),
    ("isavuconazole", 3.28, 0.01, 4.0, Organ.LUNG, 267.7),
    ("isavuconazole", 3.28, 0.01, 4.0, Organ.MUSCLE, 126.9),
    ("isavuconazole", 3.28, 0.01, 4.0, Organ.PANCREAS, 434.9),
    ("isavuconazole", 3.28, 0.01, 4.0, Organ.SKIN, 636.8),
    ("isavuconazole", 3.28, 0.01, 4.0, Organ.SPLEEN, 110.7),
    ("isavuconazole", 3.28, 0.01, 4.5, Organ.BONE, 567.2),
    ("isavuconazole", 3.28, 0.01, 4.5, Organ.BRAIN, 1252.0),
    ("isavuconazole", 3.28, 0.01, 4.5, Organ.SMALL_INTESTINE, 1306.3),
    ("isavuconazole", 3.28, 0.01, 4.5, Organ.HEART, 530.2),
    ("isavuconazole", 3.28, 0.01, 4.5, Organ.KIDNEY, 612.8),
    ("isavuconazole", 3.28, 0.01, 4.5, Organ.LIVER, 654.5),
    ("isavuconazole", 3.28, 0.01, 4.5, Organ.LUNG, 800.1),
    ("isavuconazole", 3.28, 0.01, 4.5, Organ.MUSCLE, 386.1),
    ("isavuconazole", 3.28, 0.01, 4.5, Organ.PANCREAS, 1361.2),
    ("isavuconazole", 3.28, 0.01, 4.5, Organ.SKIN, 1953.0),
    ("isavuconazole", 3.28, 0.01, 4.5, Organ.SPLEEN, 327.9),
    ("isavuconazole", 3.28, 0.01, 5.0, Organ.BONE, 1771.3),
    ("isavuconazole", 3.28, 0.01, 5.0, Organ.BRAIN, 3947.3),
    ("isavuconazole", 3.28, 0.01, 5.0, Organ.SMALL_INTESTINE, 4095.6),
    ("isavuconazole", 3.28, 0.01, 5.0, Organ.HEART, 1641.8),
    ("isavuconazole", 3.28, 0.01, 5.0, Organ.KIDNEY, 1908.4),
    ("isavuconazole", 3.28, 0.01, 5.0, Organ.LIVER, 2049.9),
    ("isavuconazole", 3.28, 0.01, 5.0, Organ.LUNG, 2483.4),
    ("isavuconazole", 3.28, 0.01, 5.0, Organ.MUSCLE, 1205.8),
    ("isavuconazole", 3.28, 0.01, 5.0, Organ.PANCREAS, 4290.1),
    ("isavuconazole", 3.28, 0.01, 5.0, Organ.SKIN, 6115.1),
    ("isavuconazole", 3.28, 0.01, 5.0, Organ.SPLEEN, 1015.0),
    ("isavuconazole", 3.28, 0.01, 5.5, Organ.BONE, 5579.0),
    ("isavuconazole", 3.28, 0.01, 5.5, Organ.BRAIN, 12470.5),
    ("isavuconazole", 3.28, 0.01, 5.5, Organ.SMALL_INTESTINE, 12915.9),
    ("isavuconazole", 3.28, 0.01, 5.5, Organ.HEART, 5156.9),
    ("isavuconazole", 3.28, 0.01, 5.5, Organ.KIDNEY, 6005.7),
    ("isavuconazole", 3.28, 0.01, 5.5, Organ.LIVER, 6462.5),
    ("isavuconazole", 3.28, 0.01, 5.5, Organ.LUNG, 7806.4),
    ("isavuconazole", 3.28, 0.01, 5.5, Organ.MUSCLE, 3797.7),
    ("isavuconazole", 3.28, 0.01, 5.5, Organ.PANCREAS, 13552.3),
    ("isavuconazole", 3.28, 0.01, 5.5, Organ.SKIN, 19276.9),
    ("isavuconazole", 3.28, 0.01, 5.5, Organ.SPLEEN, 3187.7),
]

_PLASMA_LIPIDS = (0.0023, 0.0013)
"""f_NL,P and f_NP,P of Coutinho 2024 (Materials and Methods, Eq. 3)."""


@pytest.mark.parametrize(("drug", "pka", "fu", "logp", "organ", "kpu"), _PUBLISHED)
def test_rr_published_values(
    drug: str, pka: float, fu: float, logp: float, organ: Organ, kpu: float
) -> None:
    """The weak bases of Coutinho 2024 within 0.3 (absolute Kpu).

    The tolerance is the one decimal of their tables (0.05) and the offset of
    up to 0.27 of heart, lung, spleen and pancreas, which is constant over the
    logP columns and so lies in the water terms: their spreadsheet must use
    water fractions which differ from Table III in the third digit (it is not
    published). The other tissues agree within 0.12.

    Coutinho 2024 do not set a negative plasma protein term
    `1/fu - 1 - lipids_p/Y` to zero, as the method does (Yau 2020, line 119):
    where it is negative (the high logP columns) their Kpu is lower by
    `term * albumin ratio`, which the expected value adds back.
    """
    method = _method(logp=logp, pka=(pka,), compound="base")
    p, y = 10**logp, 1 + 10 ** (pka - 7.4)
    f_nl_p, f_np_p = _PLASMA_LIPIDS
    term = 1 / fu - 1 - (p * f_nl_p + (0.3 * p + 0.7) * f_np_p) / y
    expected = kpu - min(term, 0.0) * composition(organ).albumin
    kp = rodgers_rowland(method, organ, fu, None, 0.45)
    assert kp / fu == pytest.approx(expected, abs=0.3), drug


def test_rr_neutral_without_lipophilicity() -> None:
    """A hydrophilic neutral compound without binding distributes into the water."""
    method = _method(logp=-10.0)
    c = composition(Organ.MUSCLE)
    kp = rodgers_rowland(method, Organ.MUSCLE, 1.0, None, 0.45)
    assert kp == pytest.approx(c.f_ew + c.f_iw + 0.7 * c.f_np, rel=1e-6)


def test_rr_neutral_binds_lipoproteins() -> None:
    """A bound neutral compound in muscle, computed by hand (eqs. 2, 4 of Yau 2020).

    Muscle of tissue_composition.toml: f_ew 0.118, f_iw 0.63, f_nl 0.01, f_np
    0.0072, lipoprotein ratio 0.059 (albumin ratio 0.064); plasma f_nl 0.0023,
    f_np 0.0013. logP 1 (P = 10, 0.3*P + 0.7 = 3.7), fu 0.2, X = Y = 1:

        Kpu = 0.118 + 0.63 + (10*0.01 + 3.7*0.0072)
              + (1/0.2 - 1 - (10*0.0023 + 3.7*0.0013)) * 0.059
            = 0.748 + 0.12664 + 3.97219*0.059 = 1.10899921
        Kp  = 0.2 * Kpu = 0.221799842

    With the albumin ratio the protein term would be 3.97219*0.064, Kp 0.22577.
    """
    method = _method(logp=1.0)
    kp = rodgers_rowland(method, Organ.MUSCLE, 0.2, None, 0.45)
    assert kp == pytest.approx(0.221799842, rel=1e-8)


def test_rr_strong_base_by_hand() -> None:
    """A strong base in muscle, computed by hand (eqs. 1, 3, 5 of Yau 2020).

    Base pKa 9, logP 1 (P = 10), fu 0.5, BP 1.0, HCT 0.45. Ionisation:
    X = 1 + 10^(9 - 7.0) = 101, Y = 1 + 10^(9 - 7.4) = 40.810717,
    X_bc = 1 + 10^(9 - 7.22) = 61.255959.
    Erythrocytes (f_iw 0.63, f_nl 0.0033, f_np 0.0012, AP 0.5 mg/g):

        Kpu_bc = (0.45 - 1 + 1.0)/(0.45*0.5) = 2.0
        Ka     = (2.0 - 61.255959*0.63/40.810717
                  - (10*0.0033 + 3.7*0.0012)/40.810717)
                 * 40.810717/(0.5*60.255959) = 1.4270038

    Muscle (f_ew 0.118, f_iw 0.63, f_nl 0.01, f_np 0.0072, AP 1.53 mg/g):

        Kpu = 101*0.63/40.810717 + 0.118 + (10*0.01 + 3.7*0.0072)/40.810717
              + 1.4270038*1.53*100/40.810717 = 7.0301113
        Kp  = 0.5 * Kpu = 3.5150556
    """
    method = _method(logp=1.0, pka=(9.0,), compound="base")
    kp = rodgers_rowland(method, Organ.MUSCLE, 0.5, 1.0, 0.45)
    assert kp == pytest.approx(3.5150556, rel=1e-7)


def test_rr_base_reproduces_erythrocytes(monkeypatch: pytest.MonkeyPatch) -> None:
    """The affinity of a base follows from the erythrocytes (eqs. 3, 5 of Yau 2020).

    The equation of the tissues applied to the composition of the erythrocytes
    of the data file (at the pH of the erythrocytes) gives back the erythrocyte
    to unbound plasma coefficient `Kpu_bc = (HCT - 1 + BP)/(fu*HCT)`.
    """
    data = copy.deepcopy(kp_module._data())
    data["intracellular"]["ph"] = data["erythrocytes"]["ph"]
    monkeypatch.setattr(kp_module, "_data", lambda: data)
    ery = {k: e["value"] for k, e in data["erythrocytes"].items()}
    method = _method(logp=1.5, pka=(9.2,), compound="base")
    fu, bp, hct = 0.4, 1.3, 0.45
    erythrocytes = Composition(
        "erythrocytes",
        f_ew=0.0,
        f_iw=ery["f_iw"],
        f_nl=ery["f_nl"],
        f_np=ery["f_np"],
        albumin=0.0,
        lipoprotein=0.0,
        ap=ery["ap"],
    )
    kp = rodgers_rowland(method, Organ.MUSCLE, fu, bp, hct, tissue=erythrocytes)
    assert kp / fu == pytest.approx((hct - 1 + bp) / (fu * hct), rel=1e-9)


def test_rr_strong_base_binds_phospholipids() -> None:
    """A strong base is concentrated in the tissues rich in acidic phospholipids."""
    method = _method(logp=1.0, pka=(9.5,), compound="base")
    kp = {o: rodgers_rowland(method, o, 0.5, 2.0, 0.45) for o in Organ}
    assert kp[Organ.KIDNEY] > kp[Organ.MUSCLE] > kp[Organ.BRAIN]
    assert kp[Organ.LUNG] > 10.0
    # a low BP gives a negative affinity, which is zero (Yau 2020)
    low = rodgers_rowland(method, Organ.KIDNEY, 0.5, 0.55, 0.45)
    assert low < kp[Organ.KIDNEY]
    assert low > 0


def test_rr_acid_ionisation() -> None:
    """An acid is excluded from the cells, which are more acidic than plasma."""
    neutral = rodgers_rowland(_method(logp=0.0), Organ.MUSCLE, 1.0, None, 0.45)
    acid = rodgers_rowland(
        _method(logp=0.0, pka=(4.0,), compound="acid"), Organ.MUSCLE, 1.0, None, 0.45
    )
    c = composition(Organ.MUSCLE)
    # X/Y = (1 + 10^(7.0 - 4))/(1 + 10^(7.4 - 4)) for the intracellular water
    x_y = (1 + 10**3.0) / (1 + 10**3.4)
    lipids = (c.f_nl + c.f_np) / (1 + 10**3.4)
    assert acid == pytest.approx(c.f_ew + x_y * c.f_iw + lipids, rel=1e-9)
    assert acid < neutral


def test_rr_strong_base_needs_bp() -> None:
    """Without a blood to plasma ratio a strong base cannot be computed."""
    method = _method(logp=1.0, pka=(9.5,), compound="base")
    with pytest.raises(ValueError, match="blood to plasma"):
        rodgers_rowland(method, Organ.LIVER, 0.5, None, 0.45)


def test_surrogate_tissues() -> None:
    """Organs without a composition of their own take a documented surrogate."""
    assert composition(Organ.STOMACH).tissue == "gut"
    assert composition(Organ.LARGE_INTESTINE).tissue == "gut"
    assert composition(Organ.GONADS).tissue == "muscle"
    # the acidic phospholipids and erythrocytes are ModelChoice (R43): only
    # the tissue choice tells a surrogate from an organ of its own
    data = load("tissue_composition")["organs"]
    assert data["go"]["source"] == "ModelChoice"
    assert data["li"]["source"] == "Coutinho2024"
    assert {"Yau2020", "Coutinho2024", "ModelChoice"} <= set(kp_sources(Organ.LIVER))
    assert "Huisinga2012Supplement" not in kp_sources(Organ.LIVER)
    assert {composition(o).tissue for o in Organ} >= {"adipose", "liver", "skin"}


def test_kp_values_use_sex_hematocrit() -> None:
    """The coefficients of a base depend on the reference hematocrit of the sex."""
    drug = _drug(
        _method(logp=1.0, pka=(9.5,), compound="base"),
        erythrocytes=Erythrocytes(mode="equilibrium", bp=1.5),
    )
    male, female = kp_values(drug, Sex.MALE), kp_values(drug, Sex.FEMALE)
    assert set(male) == set(Organ)
    assert reference_hematocrit(Sex.MALE) != reference_hematocrit(Sex.FEMALE)
    assert male[Organ.LIVER].value != female[Organ.LIVER].value
    hct = reference_hematocrit(Sex.MALE)
    expected = rodgers_rowland(drug.distribution, Organ.LIVER, 0.1, 1.5, hct)  # ty: ignore[invalid-argument-type]
    assert male[Organ.LIVER].value == pytest.approx(expected, rel=1e-12)


def test_kp_values_of_entered_organs() -> None:
    """Only the organs the substance enters get a coefficient."""
    organs = frozenset({Organ.LIVER, Organ.MUSCLE})
    values = kp_values(_drug(_method(organs=organs)), Sex.MALE)
    assert set(values) == organs


def test_kp_method_factor(tmp_path: Path) -> None:
    """Every `Kp_<o>_<s>` is the factor `fKp_<s>` times the value of the method.

    The factor is a parameter of the top model, a fit value carries its
    uncertainty; a change of the factor before the simulation scales every
    coefficient.
    """
    method = _method(factor=Value(1.5, U.dimensionless, "fit"))
    drug = _drug(method)
    result = build(definition(substances=[drug], chemistry=[]), tmp_path)
    values = kp_values(drug, Sex.MALE)
    r = load_model(result)
    assert r[ids.fkp("drg")] == 1.5
    for organ, v in values.items():
        assert r[ids.kp_organ(organ, "drg")] == pytest.approx(1.5 * v.value, rel=1e-12)
    assert ids.kp("drg") not in r.model.getGlobalParameterIds()

    changed = individual(result.flat_path, **{ids.fkp("drg"): 3.0})
    for organ, v in values.items():
        assert changed[ids.kp_organ(organ, "drg")] == pytest.approx(
            3.0 * v.value, rel=1e-12
        )


@pytest.mark.parametrize("mode", list(Mode))
def test_kp_method_modes_conserve_mass(tmp_path: Path, mode: Mode) -> None:
    """Perfusion, merged and full tissue build with a Kp method and keep the dose."""
    drug = _drug(
        _method(logp=1.0, pka=(9.0,), compound="base", mode=mode),
        erythrocytes=Erythrocytes(mode="equilibrium", bp=1.2),
        excretion=Excretion(),
    )
    d = definition(substances=[drug], chemistry=[], routes=[Route.iv("drg")])
    r = load_model(build(d, tmp_path))
    r[ids.iv_dose("drg")] = 10.0
    start = total_amount(r, "drg") + dose_amount(r, "drg")
    r.simulate(0, DAY, 50)
    assert start > 0
    assert total_amount(r, "drg") + dose_amount(r, "drg") == pytest.approx(
        start, rel=1e-6
    )
    tissue = {
        Mode.TISSUE: "MU__Cti_drg",
        Mode.PERMEABILITY: "MU__Cic_drg",
        Mode.PERFUSION: "MU__Cor_drg",
    }[mode]
    assert r[f"[{tissue}]"] > 0


def test_kp_method_validation() -> None:
    """The definition refuses a strong base without BP and wrong pKa counts."""
    base = _drug(_method(logp=1.0, pka=(9.0,), compound="base"))
    with pytest.raises(ValueError, match=r"drg: Rodgers & Rowland needs the blood"):
        definition(substances=[base], chemistry=[]).validate()
    weak = _drug(_method(logp=1.0, pka=(5.0,), compound="base"))
    definition(substances=[weak], chemistry=[]).validate()
    for pka, compound, match in [
        ((4.0,), "neutral", "neutral compound has no pKa"),
        ((), "acid", "at least one pKa"),
        ((3.0,), "zwitterion", "two pKa"),
        ((9.0, 3.0), "zwitterion", r"\(acidic, basic\).*acidic < basic"),
    ]:
        drug = _drug(_method(pka=pka, compound=compound))
        with pytest.raises(ValueError, match=match):
            definition(substances=[drug], chemistry=[]).validate()


def test_kp_notes_cite_sources(tmp_path: Path) -> None:
    """The coefficients cite the equations and the composition of their tissue."""
    import libsbml

    drug = _drug(_method())
    result = build(definition(substances=[drug], chemistry=[]), tmp_path)
    doc = libsbml.readSBMLFromFile(str(result.flat_path))
    notes = (
        doc.getModel().getParameter(ids.kp_organ(Organ.LIVER, "drg")).getNotesString()
    )
    for text in ("10.1208/s12248-020-0418-7", "10.1007/s11095-024-03703-4", "fKp_drg"):
        assert text in notes
    assert load("tissue_composition")["organs"]["li"]["tissue"] == "liver"


def test_rr_zwitterion_order() -> None:
    """A zwitterion takes `pka = (acidic, basic)`; a swapped tuple is refused."""
    method = _method(logp=1.0, pka=(9.0, 3.0), compound="zwitterion")
    with pytest.raises(ValueError, match=r"\(acidic, basic\)"):
        rodgers_rowland(method, Organ.MUSCLE, 0.5, 1.0, 0.45)
    ordered = _method(logp=1.0, pka=(3.0, 9.0), compound="zwitterion")
    assert rodgers_rowland(ordered, Organ.MUSCLE, 0.5, 1.0, 0.45) > 0


def test_full_mode_refuses_bound_acid_in_adipose(tmp_path: Path) -> None:
    """A highly bound hydrophilic acid is refused in the full form (documented).

    Its Rodgers & Rowland Kp of adipose (rat albumin ratio 0.049) lies below
    the interstitial share of the Schmitt `K_int_pls` of the library, so the
    full form would need `K_cell_pls <= 0`; merged and perfusion build.
    """
    acid = _method(logp=1.0, pka=(4.0,), compound="acid", mode=Mode.PERMEABILITY)
    drug = _drug(acid, binding=Binding(fu=0.01, protein="albumin"))
    d = definition(substances=[drug], chemistry=[], routes=[Route.iv("drg")])
    with pytest.raises(ValueError, match=r"drg: Kp .* organ 'ad' .*interstitial share"):
        build(d, tmp_path / "full")
    merged = dataclasses.replace(
        drug, distribution=dataclasses.replace(acid, mode=Mode.TISSUE)
    )
    build(
        definition(substances=[merged], chemistry=[], routes=[Route.iv("drg")]),
        tmp_path / "merged",
    )
