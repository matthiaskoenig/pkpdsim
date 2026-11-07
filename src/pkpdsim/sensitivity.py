"""The PK and PBPK parts of the sensitivity analyses of the models.

The runner of the sensitivity analyses of a model and the table of their
parameters are sbmlsim's: `sbmlsim.sensitivity.runner` (`SensitivityRunner`,
`SensitivityGroup`, `SensitivityOutput`, `LocalAnalysis`, `SamplingAnalysis`,
`FastAnalysis`, `SobolAnalysis`, `run_subfolders`, `sample_size` with
`SBMLSIM_SENSITIVITY_N`) and `sbmlsim.sensitivity.parameters`
(`parameters_from_sbml`, `set_bounds`, `set_relative_bounds`,
`write_parameters`, `free_parameters`, ...). A `PKPDModel` is the project of
a runner, its builds are the variants of the groups (`variant=Sex.FEMALE`).
This module holds what is pharmacokinetic or of the PBPK library: the
observables of the non-compartmental analysis of a substance
(`pk_observables`, the `observables` of a runner, `max_after_dose`), the
parameters of an analysis of a library model (`library_parameters`,
`notes_sources`, `SEVERITIES`) and the local analysis of its nonzero
parameters (`NonzeroLocalAnalysis`):

```python
# not executed
SENSITIVITY = SensitivityRunner(
    project=MODEL,
    parameters=sensitivity_parameters,
    outputs=OUTPUTS,
    observables=lambda: pk_observables(
        "hctz", "[Cve_hctz]", {"cmax": "mM"}, dose="PODOSE_hctz", pk_id="pk_hctz"
    ),
    analyses={
        "local": NonzeroLocalAnalysis(groups=GROUPS[:1]),
        "sampling": SamplingAnalysis(groups=GROUPS),
    },
    out=MODEL.results_dir / "sensitivity",
)

if __name__ == "__main__":
    SENSITIVITY.main()
```
"""

from __future__ import annotations

from collections.abc import Collection, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import override
from xml.etree import ElementTree

import libsbml
import numpy as np
from pkpdutils import NCAOptions
from sbmlsim.sensitivity.parameters import (
    ParameterType,
    SensitivityParameter,
    parameters_from_sbml,
    set_bounds,
    set_relative_bounds,
)
from sbmlsim.sensitivity.runner import LocalAnalysis
from sbmlsim.simulation import PK, Formula, Observable


def pk_observables(
    prefix: str,
    selection: str,
    outputs: Mapping[str, str],
    *,
    dose: str | None = None,
    route: str | None = None,
    options: NCAOptions | None = None,
    per_dose_mass: Collection[str] = (),
    mr: str | None = None,
    pk_id: str | None = None,
) -> list[Observable]:
    """Get the non-compartmental analysis of a concentration and its outputs.

    The observables of a sensitivity analysis of the pharmacokinetics of a
    substance: the analysis `pk_id` (sbmlsim `PK`) of the concentration
    `selection` with the parameters of the outputs, then an observable
    `<prefix>_<key>` of every output in its unit (sbmlsim `Formula`). The
    dose of an analysis is a mass (`PODOSE_<s>` in mg), so an output whose
    unit has the amount of the dose (`per_dose_mass`, e.g. `cl_f` and `vz_f`)
    is divided by the molecular weight `mr` at the start (`at(Mr_<s>, 0)`, a
    value per simulation like the parameter of the analysis); the molecular
    weight of the dosed substance, also for a metabolite of it. The
    observables are the `observables` of a
    `sbmlsim.sensitivity.runner.SensitivityRunner`, as a list or as a function
    which builds them when a run starts, whose outputs
    (`SensitivityOutput(observable=...)`) name the ids `<prefix>_<key>` or the
    parameters `<pk_id>.<key>` of the analysis directly.

    ```python
    # not executed
    OBSERVABLES = [
        *pk_observables(
            "tor",
            "Carm_tor",
            {"auc_inf_obs": "mM*min", "cl_f": "ml/min"},
            dose="PODOSE_tor",
            route="oral",
            options=NCAOptions(lloq=LLOQ),
            per_dose_mass={"cl_f", "vz_f"},
            mr="Mr_tor",
        ),
        *pk_observables("Carm_tor_m5", "Carm_tor_m5", {"cmax": "mM"}, pk_id="pk_m5"),
    ]
    ```

    Args:
        prefix: the prefix of the ids of the outputs, e.g. the substance
            `tor` (`tor_cmax`) or the concentration `Carm_api`
            (`Carm_api_cmax`).
        selection: the concentration, e.g. `Carm_tor`.
        outputs: the unit of every output by the key of its parameter of the
            analysis (`pkpdutils`, e.g. `auc_inf_obs`, `cmax`, `cl_f`).
        dose: the dose state of the analysis, `None` for an analysis without
            a dose (a metabolite), see `sbmlsim.simulation.PK`.
        route: the route of the dose, e.g. `oral` or `iv_bolus`.
        options: the options of the analysis, e.g. the limit of
            quantification relative to the maximum (`NCAOptions(lloq=...)`).
        per_dose_mass: the outputs whose unit has the amount of the dose.
        mr: the id of the molecular weight of the dosed substance, which an
            output of `per_dose_mass` is divided by, e.g. `Mr_tor`.
        pk_id: the id of the analysis, `pk_<prefix>` by default; an id of the
            model, e.g. a species `hctz`, is no id of an observable.

    Returns:
        The analysis and the outputs.

    Raises:
        ValueError: if an output per dose mass has no molecular weight.
    """
    pk = f"pk_{prefix}" if pk_id is None else pk_id
    divided = [key for key in outputs if key in per_dose_mass]
    if divided and mr is None:
        raise ValueError(
            f"The outputs {divided} of '{prefix}' have the amount of the dose in "
            "their unit: give the molecular weight of the dosed substance, `mr`."
        )
    return [
        PK(
            pk,
            selection,
            dose=dose,
            route=route,
            options=options,
            parameters=list(outputs),
        ),
        *(
            Formula(
                f"{prefix}_{key}",
                f"{pk}.{key} / at({mr}, 0)" if key in per_dose_mass else f"{pk}.{key}",
                unit=unit,
            )
            for key, unit in outputs.items()
        ),
    ]


def max_after_dose(time: np.ndarray, values: Mapping[str, np.ndarray]) -> float:
    """Get the maximum of a timecourse after the dose at the time 0.

    The function of a `Custom` observable of one selection, e.g. the largest
    diuresis after a dose of a simulation which starts before it:
    `Custom("max_diuresis", max_after_dose, "ml/min", symbols=["FLUID__diuresis"])`.

    Args:
        time: the times of the timecourse.
        values: the timecourse of the one selection of the observable.

    Returns:
        The maximum from the time 0 on.

    Raises:
        ValueError: for more or less than one selection.
    """
    if len(values) != 1:
        raise ValueError(
            f"The maximum after the dose is of one selection, not of {sorted(values)}."
        )
    (x,) = values.values()
    return float(np.max(x[time >= 0]))


#: the severities of the impairments of a library model, 0 in the healthy
#: individual, by their largest value (Child-Pugh C, NYHA IV)
SEVERITIES: Mapping[str, float] = MappingProxyType(
    {"CP_severity": 3.0, "HF_severity": 4.0}
)


@dataclass(frozen=True, kw_only=True)
class NonzeroLocalAnalysis(LocalAnalysis):
    """The local analysis of the parameters whose value is not 0.

    A local design changes a value relative to itself, so a severity at 0
    (the healthy individual) is not changed and its sensitivity is 0/0; the
    sampled analyses vary the severities within their bounds.
    """

    @override
    def varied(
        self, parameters: Sequence[SensitivityParameter]
    ) -> list[SensitivityParameter]:
        return [p for p in super().varied(parameters) if p.value != 0]


def notes_sources(sbase: libsbml.SBase) -> list[str]:
    """Get the sources in the notes of an element of a library model.

    The library writes the sources of a value as a list after a paragraph
    `Sources:` of its notes.

    Args:
        sbase: the element.

    Returns:
        The sources, one text per list item; empty without notes or sources.
    """
    if not sbase.isSetNotes():
        return []
    root = ElementTree.fromstring(sbase.getNotesString())
    elements = list(root.iter())
    for k, element in enumerate(elements):
        if element.tag.endswith("}p") and "".join(element.itertext()).strip() == (
            "Sources:"
        ):
            lists = [e for e in elements[k + 1 :] if e.tag.endswith("}ul")]
            if lists:
                return [
                    " ".join("".join(item.itertext()).split())
                    for item in lists[0]
                    if item.tag.endswith("}li")
                ]
    return []


def library_parameters(
    model_path: Path,
    ids: Sequence[str],
    *,
    fit_ids: Collection[str],
    scaling_ids: Collection[str],
    references: Mapping[str, str] = MappingProxyType({}),
    relative: float = 0.15,
) -> list[SensitivityParameter]:
    """Get the parameters of the analyses of a model of the PBPK library.

    A library model has hundreds of constants, most of them initial
    assignments of the covariates, which an analysis of all constants would
    vary as well; the analysis names the parameters it varies. They are the
    `ids`, in this order, with these bounds:

    - fitted parameters (`fit_ids`, type fitted): `value * (1 +- relative)`;
    - scaling factors (`scaling_ids`, type scaling, 1 in the healthy
      individual): `1 +- relative`;
    - severities (`SEVERITIES`, type scaling): 0 to `relative * top`, since
      relative bounds around 0 would be 0;
    - any other parameter: `value * (1 +- relative)`.

    The reference of a parameter is the sources of its notes (`;` separated),
    replaced by the text in `references` (e.g. the guideline of the classes of
    the renal function, whose notes give no source). A parameter of the type
    not set with a reference is of the type data.

    Args:
        model_path: the SBML model.
        ids: the parameters of the analysis.
        fit_ids: the fitted parameters among the `ids`.
        scaling_ids: the scaling factors among the `ids`.
        references: references by id, which replace those of the notes.
        relative: the relative half width of the bounds.

    Returns:
        The parameters of the `ids` with bounds, types and references.

    Raises:
        ValueError: naming every id which is no parameter of the model, and
            every id of `fit_ids`, `scaling_ids` or `references` which is not
            in `ids`.
    """
    by_uid = {
        p.uid: p
        for p in parameters_from_sbml(model_path, exclude_na=False, exclude_zero=False)
    }
    missing = [uid for uid in ids if uid not in by_uid]
    if missing:
        raise ValueError(
            f"The parameters {missing} of the sensitivity analysis are no "
            f"parameters of the model {model_path.name}."
        )
    outside = sorted(
        {*fit_ids, *scaling_ids, *references}.difference(ids),
    )
    if outside:
        raise ValueError(
            f"The ids {outside} of the sensitivity analysis are not in its parameters."
        )
    parameters = [by_uid[uid] for uid in ids]

    # the fitted parameters, their bounds are set relative to the value below
    set_bounds(
        parameters, [(uid, np.nan, np.nan, ParameterType.FIT) for uid in fit_ids]
    )
    # the scaling factors of the organ functions and of the protein activities
    set_bounds(
        parameters,
        [
            (uid, 1 - relative, 1 + relative, ParameterType.SCALING)
            for uid in scaling_ids
        ],
    )
    # the severities, whose relative bounds around 0 would be 0
    set_bounds(
        parameters,
        [
            (uid, 0.0, relative * top, ParameterType.SCALING)
            for uid, top in SEVERITIES.items()
        ],
    )

    # the references of the values, as the sources of their notes
    model: libsbml.Model = libsbml.readSBMLFromFile(str(model_path)).getModel()
    texts = {
        uid: "; ".join(sources)
        for uid in ids
        if (sources := notes_sources(model.getParameter(uid)))
    }
    texts.update(references)
    for uid, reference in texts.items():
        p = by_uid[uid]
        p.reference = reference
        if p.type == ParameterType.NA:
            p.type = ParameterType.DATA

    # every parameter without bounds varies around its value
    set_relative_bounds(parameters, relative)
    return parameters
