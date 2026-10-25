"""Tests of the distribution coefficients of the standard organ."""

import dataclasses

import pytest

from pkpdsim.pbpk import Kp, Mode, Organ, distribution
from pkpdsim.pbpk.physiology import Sex, reference_fractions

from .conftest import simple_drug


def test_k_int_pls_value_and_rule_agree() -> None:
    """The number is the rule evaluated with the constants (one template)."""
    fu = 0.1
    expected = (0.935 + 0.37 * (1 / fu - 0.926)) * fu
    assert distribution.k_int_pls_value(fu) == pytest.approx(expected, rel=1e-12)
    assert distribution.k_int_pls("fu_drg") == (
        "(f_water_int + r_protein_int_pls * (1 dimensionless / fu_drg - "
        "f_water_pls)) * fu_drg"
    )


def test_k_cell_pls_gives_whole_tissue_kp() -> None:
    """Interstitium and cells together have the Kp."""
    f_vas, f_int = reference_fractions(Organ.LIVER)[Sex.MALE]
    f_cell = 1 - f_vas - f_int
    k_int = distribution.k_int_pls_value(0.1)
    k_cell = distribution.k_cell_pls_value(2.0, 0.1, f_vas, f_int)
    assert (f_int * k_int + f_cell * k_cell) / (f_int + f_cell) == pytest.approx(2.0)


def test_check_kp() -> None:
    """A Kp below the interstitial share is refused, above it accepted."""
    drug = simple_drug()
    distribution.check_kp(drug, Organ.SKIN)
    low = Kp(0.01, mode=Mode.PERMEABILITY, ftissue=1.0)
    with pytest.raises(ValueError, match=r"drg.*'sk'.*interstitial share"):
        distribution.check_kp(dataclasses.replace(drug, distribution=low), Organ.SKIN)
