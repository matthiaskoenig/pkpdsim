"""Tests of the temperature of a machine perfusion."""

import math
from pathlib import Path

import numpy as np
import pytest
from pbpk.conftest import individual

from pkpdsim.database.tables import entry
from pkpdsim.perfusion import ids
from pkpdsim.perfusion.builder import build

from .conftest import liver_definition

R_GAS = 8.314462618
"""J/(mol*K), the exact molar gas constant of the SI since 2019."""


def _number(section: str, key: str) -> float:
    value = entry("perfusion", section, key).value
    assert isinstance(value, float)
    return value


@pytest.fixture
def path(tmp_path: Path) -> Path:
    return build(liver_definition(2), tmp_path)


def test_f_t_is_the_arrhenius_factor(path: Path) -> None:
    """`f_T` is the Arrhenius factor of the activation energy of the
    database, 1 at the reference temperature."""
    t_ref = _number("temperature", "reference") + 273.15
    ea = _number("temperature", "activation_energy_factor") * R_GAS * t_ref
    for celsius in (4.0, 20.0, 37.0):
        r = individual(path, **{ids.TEMPERATURE: celsius})
        expected = math.exp(ea / R_GAS * (1 / t_ref - 1 / (celsius + 273.15)))
        assert r[ids.F_T] == pytest.approx(expected, rel=1e-12), celsius
    assert r[ids.F_T] == pytest.approx(1.0, rel=1e-12)


def test_the_oxygen_consumption_falls_with_the_temperature(path: Path) -> None:
    """The oxygen uptake of the organ after an hour falls from 37 to 20 to
    4 °C, the solubility of oxygen rises."""
    uptake, solubility = [], []
    for celsius in (37.0, 20.0, 4.0):
        r = individual(path, **{ids.TEMPERATURE: celsius})
        solubility.append(r[ids.SOLUBILITY])
        r.simulate(0, 60, 61)
        uptake.append(r[ids.VO2])
    assert uptake[0] > uptake[1] > uptake[2] > 0
    assert solubility[0] < solubility[1] < solubility[2]


def test_a_cold_perfusion_runs_six_hours(path: Path) -> None:
    r = individual(path, **{ids.TEMPERATURE: 4.0})
    s = r.simulate(0, 360, 361)
    assert np.isfinite(np.asarray(s)).all()
    assert r[ids.alive("z1", submodel=True)] > 0
