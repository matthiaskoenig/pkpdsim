"""The PK parts of the sensitivity analyses, on the toy models.

The runner and the table of the parameters are sbmlsim's
(`sbmlsim.sensitivity.runner`, `sbmlsim.sensitivity.parameters`); a
`PKPDModel` is the project of a runner, with its builds per sex as the
variants of the groups.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from sbmlsim.sensitivity.parameters import (
    SensitivityParameter,
    parameters_from_sbml,
    set_relative_bounds,
)
from sbmlsim.sensitivity.runner import (
    SamplingAnalysis,
    SensitivityGroup,
    SensitivityOutput,
    SensitivityRunner,
)
from sbmlsim.simulation import Simulation

from pkpdsim.sensitivity import NonzeroLocalAnalysis, pk_observables
from pkpdsim.sex import Sex


def test_the_pk_observables_of_a_substance() -> None:
    """The analysis of a concentration and its outputs, as the copies of the models had them."""
    from pkpdutils import NCAOptions
    from sbmlsim.simulation import PK, Formula

    from pkpdsim.sensitivity import pk_observables

    outputs = {"auc_inf_obs": "mM*min", "cmax": "mM", "cl_f": "ml/min"}
    options = NCAOptions(lloq=1e-6)
    observables = pk_observables(
        "gli",
        "Carm_gli",
        outputs,
        dose="PODOSE_gli",
        route="oral",
        options=options,
        per_dose_mass={"cl_f", "vz_f"},
        mr="Mr_gli",
    )
    assert observables == [
        PK(
            "pk_gli",
            "Carm_gli",
            dose="PODOSE_gli",
            route="oral",
            options=options,
            parameters=list(outputs),
        ),
        Formula("gli_auc_inf_obs", "pk_gli.auc_inf_obs", unit="mM*min"),
        Formula("gli_cmax", "pk_gli.cmax", unit="mM"),
        Formula("gli_cl_f", "pk_gli.cl_f / at(Mr_gli, 0)", unit="ml/min"),
    ]
    # a metabolite without a dose, its outputs named by the concentration
    metabolite = pk_observables("Carm_m1", "Carm_m1", {"cmax": "mM"}, pk_id="pk_m1")
    assert metabolite == [
        PK("pk_m1", "Carm_m1", parameters=["cmax"]),
        Formula("Carm_m1_cmax", "pk_m1.cmax", unit="mM"),
    ]


def test_an_output_per_dose_mass_needs_the_molecular_weight() -> None:
    from pkpdsim.sensitivity import pk_observables

    with pytest.raises(ValueError, match="molecular weight"):
        pk_observables("gli", "Carm_gli", {"cl_f": "ml/min"}, per_dose_mass={"cl_f"})


def test_the_maximum_after_the_dose() -> None:
    """The maximum of the timecourse from the time of the dose, 0, on."""
    from pkpdsim.sensitivity import max_after_dose

    time = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
    values = {"FLUID__diuresis": np.array([9.0, 5.0, 1.0, 3.0, 2.0])}
    assert max_after_dose(time, values) == 3.0
    with pytest.raises(ValueError, match="one selection"):
        max_after_dose(time, {"a": time, "b": time})


def _sexes_parameters() -> list[SensitivityParameter]:
    """The parameters of the toy model per sex, read from its male build."""
    from pkpdsim_sexes import MODEL_PATH

    parameters = parameters_from_sbml(MODEL_PATH)
    set_relative_bounds(parameters, 0.15)
    return parameters


def _sexes_runner(out: Path | None) -> SensitivityRunner:
    """The analyses of the toy model per sex with the PK observables of `A`."""
    simulation = Simulation(end=60, steps=600, time_unit="min")
    groups = [
        SensitivityGroup("men", "Men", simulation),
        SensitivityGroup("women", "Women", simulation, variant=Sex.FEMALE),
    ]
    return SensitivityRunner(
        parameters=_sexes_parameters,
        outputs=[
            SensitivityOutput("A_auc_inf_obs", "A AUC∞"),
            SensitivityOutput("A_thalf", "A half-life", observable="pk_A.thalf"),
        ],
        observables=lambda: pk_observables(
            "A", "[A]", {"auc_inf_obs": "mM*min", "thalf": "min"}
        ),
        analyses={
            "local": NonzeroLocalAnalysis(groups=groups),
            "sampling": SamplingAnalysis(groups=groups, n=8),
        },
        out=out,
    )


def test_a_pkpd_model_runs_the_pk_observables_per_sex(tmp_path: Path) -> None:
    """The runner of sbmlsim on a `PKPDModel`: its builds are the variants."""
    from pkpdsim_sexes import FEMALE_A0
    from pkpdsim_sexes import MODEL as SEXES

    runner = _sexes_runner(tmp_path / "sensitivity")
    model = SEXES.replace(sensitivities={"decay": runner})
    bound = model.sensitivity("decay")
    assert bound.project is model
    path = bound.run(["local", "sampling"], workers=1)
    assert path == tmp_path / "sensitivity"
    for k, group in enumerate(("men", "women")):
        assert (path / "local" / f"local_d0.01_nvar3_{k:02d}_{group}.tsv").is_file()
    samples = {
        group: pd.read_csv(path / "sampling" / group / "samples.tsv", sep="\t")
        for group in ("men", "women")
    }
    men, women = samples["men"], samples["women"]
    # A0 differs by sex, the same relative design around the values of a build
    np.testing.assert_allclose(women["A0"], FEMALE_A0 * men["A0"], rtol=1e-12)
    np.testing.assert_allclose(
        women["A_auc_inf_obs"], FEMALE_A0 * men["A_auc_inf_obs"], rtol=1e-6
    )
    # the half-life of the PK analysis, named by the output directly
    np.testing.assert_allclose(women["A_thalf"], men["A_thalf"], rtol=1e-6)
    table = pd.read_csv(path / "parameters.tsv", sep="\t").set_index("uid")
    assert table.loc["A0", "value_female"] == pytest.approx(FEMALE_A0)
