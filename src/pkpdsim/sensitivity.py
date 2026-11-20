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
`notes_sources`, `SEVERITY_IDS`) and the local analysis of its nonzero
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

import logging
import math
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

from pkpdsim.sbml.domains import Domain, domains

logger = logging.getLogger(__name__)


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


#: the ids of the severities of the impairments of a library model, 0 in the
#: healthy individual; their largest value (Child-Pugh C 3, NYHA IV 4) is the
#: upper bound of their domain in the model
SEVERITY_IDS: tuple[str, ...] = ("CP_severity", "HF_severity")


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
    clip: bool = True,
    unclipped: Collection[str] = (),
) -> list[SensitivityParameter]:
    """Get the parameters of the analyses of a model of the PBPK library.

    A library model has hundreds of constants, most of them initial
    assignments of the covariates, which an analysis of all constants would
    vary as well; the analysis names the parameters it varies. They are the
    `ids`, in this order, with these bounds:

    - fitted parameters (`fit_ids`, type fitted): `value * (1 +- relative)`;
    - scaling factors (`scaling_ids`, type scaling, 1 in the healthy
      individual): `1 +- relative`;
    - severities (`SEVERITY_IDS`, type scaling): 0 to `relative * top`, the top
      being the upper bound of the domain of the severity in the model, since
      relative bounds around 0 would be 0;
    - any other parameter: `value * (1 +- relative)`.

    A relative span does not know what a value means: 0.94 of a fraction
    absorbed `+- 0.15` reaches 1.081, more than a dose. With `clip` the bounds
    of every parameter are intersected with the domain the model states for it
    (`pkpdsim.sbml.domains.domains`), and every bound which moves is logged as
    a warning. A parameter which states no domain keeps its bounds, which one
    summary warning counts. A model which states no domain for any quantity
    was built before the domain contract: it gets the plain relative span for
    every parameter, severities included, and one warning, and raises nothing.
    A pair of bounds whose lower bound is above its
    upper bound (the relative bounds of a negative value) is left as it is. The
    two ways out of the clip are `clip=False` for the whole analysis and
    `unclipped` for single parameters, which keep the span they were given.

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
        clip: intersect the bounds with the domains of the model.
        unclipped: the ids whose bounds are not intersected with their domain.

    Returns:
        The parameters of the `ids` with bounds, types and references.

    Raises:
        ValueError: naming every id which is no parameter of the model, and
            every id of `fit_ids`, `scaling_ids`, `references` or `unclipped`
            which is not in `ids`, a severity without a finite upper bound of
            its domain in a model which states domains, and a clipped
            parameter whose value is outside its domain.
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
        {*fit_ids, *scaling_ids, *references, *unclipped}.difference(ids),
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
    # the severities, whose relative bounds around 0 would be 0, up to the top
    # of their domain
    domain_of = domains(model_path)
    # a model which states no domain at all predates the contract: the plain
    # relative span for every parameter, no clip
    predates = not domain_of
    if predates:
        logger.warning(
            "The model %s states no domain for any quantity, so it was built "
            "before the domain contract: the bounds are the plain relative span, "
            "no severity is varied over its domain and nothing is clipped. "
            "Rebuilding the model gives the clip.",
            model_path.name,
        )
    severities = [] if predates else [uid for uid in ids if uid in SEVERITY_IDS]
    for uid in severities:
        if uid not in domain_of or not math.isfinite(domain_of[uid].upper):
            raise ValueError(
                f"The severity {uid} of the sensitivity analysis has no finite "
                f"upper bound of its domain in the model {model_path.name}: the "
                "model states domains but none for this quantity, or one without "
                "an upper bound. A model built before the domain contract has to "
                "be rebuilt."
            )
    set_bounds(
        parameters,
        [
            (uid, 0.0, relative * domain_of[uid].upper, ParameterType.SCALING)
            for uid in severities
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
    if clip and not predates:
        _clip_to_domains(parameters, domain_of, set(unclipped), model_path.name)
    return parameters


def _clip_to_domains(
    parameters: Sequence[SensitivityParameter],
    domain_of: Mapping[str, Domain],
    unclipped: Collection[str],
    name: str,
) -> None:
    """Intersect the bounds of the parameters with their domains, in place.

    Args:
        parameters: the parameters with bounds.
        domain_of: the domains of the model by id.
        unclipped: the ids which keep their bounds.
        name: the name of the model for the messages.

    Raises:
        ValueError: for a parameter whose value is outside its domain.
    """
    without = []
    not_a_number: list[str] = []
    for p in parameters:
        if p.uid in unclipped:
            continue
        domain = domain_of.get(p.uid)
        if domain is None:
            without.append(p.uid)
            continue
        if math.isnan(p.value):
            not_a_number.append(p.uid)
            continue
        if not domain.contains(p.value):
            raise ValueError(
                f"The value {p.value:g} of the parameter {p.uid} is outside its "
                f"domain {domain} in the model {name}."
            )
        lower, upper = p.lower_bound, p.upper_bound
        if math.isnan(lower) or math.isnan(upper):
            not_a_number.append(p.uid)
            continue
        clipped = domain.clip(lower, upper)
        if clipped != (lower, upper):
            logger.warning(
                "The bounds [%g, %g] of the sensitivity parameter %s leave its "
                "domain %s and are clipped to [%g, %g].",
                lower,
                upper,
                p.uid,
                domain,
                *clipped,
            )
            p.lower_bound, p.upper_bound = clipped
    if without:
        logger.warning(
            "%d of the %d parameters of the sensitivity analysis state no domain "
            "in the model %s, their bounds are not clipped: %s.",
            len(without),
            len(parameters),
            name,
            ", ".join(without),
        )
    if not_a_number:
        logger.warning(
            "%d of the %d parameters of the sensitivity analysis have a value or a "
            "bound which is not a number in the model %s, they are not clipped: %s.",
            len(set(not_a_number)),
            len(parameters),
            name,
            ", ".join(dict.fromkeys(not_a_number)),
        )
