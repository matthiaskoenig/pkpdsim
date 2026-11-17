"""The default web app of a model: pages of its doses, its diseases and its scans.

The factory of a model writes its web app with its models (`sbmlsim.webapp`,
`sbmlsim.project.Project.create_models`). A model which declares its own
panels and app shows them; every other model of pkpdsim shows the app of
`default_app`, which is built from what the model and its simulation
experiments state, so it needs nothing written by hand:

- the readouts: one small panel per output which the experiments with data
  compare with their studies, the observables of their fit mappings
  (`READOUTS` at most, the most compared first), in their units of the model
  (a concentration in uM and an amount in umol);
- the page `dose`: a knob of the dose of every route and substance the
  experiments with data dose (`PODOSE_<s>`, `IVDOSE_<s>`, ... of
  `DOSE_TARGETS`, `DOSE_KNOBS` at most, the most studied first), a range from
  0 to the largest studied dose with the studied doses as its stops, the most
  studied dose of the first route as its default;
- the page `protocol`: the doses repeated by a dosing protocol
  (`sbmlsim.panel.PanelDosing`), a choice of the interval (a single dose,
  once, twice or three times daily, `INTERVALS`) and the days of the dosing,
  over `PROTOCOL_DAYS` days;
- a page per disease the model is studied in, the subjects of its studies
  (`pkpdsim.metadata.Health` of the metadata of the fit mappings) or a scan
  of its parameter: the classes of the renal function, of liver cirrhosis and
  of heart failure as a choice knob, the classes its experiments simulate (the
  maps `renal_map`, `cirrhosis_map` and `cardiac_map` of the base of its
  experiments, e.g. the renal functions of the classes of the GFR at the
  healthy GFR of the model) on the parameter of the model
  (`f_renal_function` or `KI__f_renal_function`, `CP_severity` or
  `f_cirrhosis`, `HF_severity` or `f_cardiac_function`);
- a page per other scan of the scan experiments of the model (their
  `scan_map`, e.g. the body weight) as a range knob over the scanned range.

Every page shows the readouts as small panels (`AppPage.size`) under the
knobs it shares, the dose knobs on every page, and the app names the variants
of the model (the builds of the sexes), which a reader switches between.
"""

from __future__ import annotations

import contextlib
import logging
import math
import re
import statistics
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Final

import libsbml
from pint.errors import PintError
from sbmlsim.panel import (
    MAX_END,
    App,
    AppPage,
    AppVariant,
    Choice,
    Knob,
    Panel,
    PanelCurve,
    PanelDosing,
    PanelSimulation,
)
from sbmlsim.units import Quantity, ureg

from pkpdsim.sex import Sex

if TYPE_CHECKING:
    from sbmlsim.experiment import SimulationExperiment

    from pkpdsim.model import PKPDModel

logger = logging.getLogger(__name__)

__all__ = [
    "DOSE_KNOBS",
    "DOSE_TARGETS",
    "INTERVALS",
    "PROTOCOL_DAYS",
    "READOUTS",
    "default_app",
]

#: the targets of a dose of the PBPK library: the route and the substance
DOSE_TARGETS: Final = re.compile(r"^(PO|IV|SOL|IM|SC|REC|INH)DOSE_(\w+)$")
#: the routes of the dose targets as they are shown
ROUTE_LABELS: Final = {
    "PO": "oral",
    "IV": "intravenous",
    "SOL": "oral solution",
    "IM": "intramuscular",
    "SC": "subcutaneous",
    "REC": "rectal",
    "INH": "inhaled",
}
#: the most readouts of a page
READOUTS: Final = 8
#: the most dose knobs of a page
DOSE_KNOBS: Final = 3
#: the intervals of the dosing protocols, in hours, 0 a single dose
INTERVALS: Final = (
    ("single dose", 0.0),
    ("once daily", 24.0),
    ("twice daily", 12.0),
    ("three times daily", 8.0),
)
#: the days of the page `protocol`
PROTOCOL_DAYS: Final = 6
#: the output points of a page of single doses and of the page `protocol`
POINTS: Final = 300
PROTOCOL_POINTS: Final = 600


@dataclass(frozen=True)
class _Disease:
    """A disease with the map of its classes and the parameter they set.

    Attributes:
        key: the id of its page.
        title: the title of its page.
        label: the label of its knob.
        map: the attribute of the base of the experiments of a model which
            maps a class to the value of its parameter (`renal_map`,
            `cirrhosis_map`, `cardiac_map`), the classes the experiments of
            the model use.
        targets: the parameters a map sets, the first a model has: the one of
            the PBPK library and the one of the older models.
        words: the words of the health of the subjects of a study which have it.
    """

    key: str
    title: str
    label: str
    map: str
    targets: tuple[str, ...]
    words: tuple[str, ...]


#: the diseases of the models, in the order of their pages
DISEASES: Final = (
    _Disease(
        "renal",
        "Renal impairment",
        "renal function",
        "renal_map",
        ("f_renal_function", "KI__f_renal_function"),
        ("renal",),
    ),
    _Disease(
        "hepatic",
        "Liver cirrhosis",
        "liver cirrhosis",
        "cirrhosis_map",
        ("CP_severity", "f_cirrhosis"),
        ("cirrhosis", "hepatic", "hepato"),
    ),
    _Disease(
        "cardiac",
        "Heart failure",
        "cardiac function",
        "cardiac_map",
        ("HF_severity", "f_cardiac_function"),
        ("heart failure", "cardiac"),
    ),
)


@dataclass(frozen=True)
class _Studied:
    """What the simulation experiments with data of a model state.

    Attributes:
        readouts: the selections the fit mappings compare with data, the most
            compared first.
        doses: the doses in mg of every dose target, by the target, the most
            studied target first.
        ends: the ends of the simulations which dose, in hours.
        health: the health of the subjects of the studies.
        scans: the scans of the scan experiments, the `scan_map` entries.
        base: the base class of the experiments of the model, whose maps give
            the classes of the diseases, `None` for none.
    """

    readouts: tuple[str, ...]
    doses: Mapping[str, tuple[float, ...]]
    ends: tuple[float, ...]
    health: frozenset[str]
    scans: tuple[Mapping[str, Any], ...]
    base: type | None


def default_app(model: PKPDModel, sbml: Path) -> tuple[dict[str, Panel], App | None]:
    """Get the default panels and app of a model, see the module.

    Args:
        model: the model.
        sbml: the SBML of its default variant, which the factory just wrote.

    Returns:
        The panels by their id and the app; the panels of `Project.default_panels`
        and no app for a model whose experiments state no readout.
    """
    ids = _entities(sbml)
    studied = _studied(model, ids)
    readouts = [r for r in studied.readouts if r in ids][:READOUTS]
    if not readouts:
        logger.info(
            "The model '%s' compares no output of it with data, its app shows its species.",
            model.name,
        )
        return model.default_panels(sbml), None
    units = _units(sbml, readouts)
    labels = {r: _label(sbml, r) for r in readouts}
    time_unit = _time_unit(sbml)
    dose_knobs = _dose_knobs(
        studied.doses,
        ids,
        {t: _label(sbml, t) for t in studied.doses},
        _substance_names(sbml),
    )
    zero = {
        target: 0.0
        for target in sorted(ids)
        if DOSE_TARGETS.fullmatch(target)
        and target not in {k.target for k in dose_knobs}
    }
    sex = _sex(next(iter(model.variants), ""))
    end = _end(studied.ends, time_unit)
    single = PanelSimulation(end=end, points=POINTS, x_unit="hr", preinit=zero)

    pages: list[tuple[str, str, tuple[Knob, ...], PanelSimulation]] = []
    if dose_knobs:
        pages.append(("dose", "Dose", tuple(dose_knobs), single))
        pages.append(_protocol_page(dose_knobs, zero, time_unit))
    first = tuple(dose_knobs[:1])
    covered: set[str] = {str(k.target) for k in dose_knobs}
    for disease in DISEASES:
        target = next((t for t in disease.targets if t in ids), None)
        if target is None or not _studies(studied, disease, target):
            continue
        knob = _disease_knob(disease, target, studied.base, sex)
        if knob is not None:
            covered.add(target)
            pages.append((disease.key, disease.title, (knob, *first), single))
    for spec in studied.scans:
        knob = _scan_knob(spec, ids, covered)
        if knob is not None:
            covered.add(str(knob.target))
            title = re.sub(r"\s*\[[^\]]*\]$", "", str(spec.get("label") or knob.label))
            pages.append((f"scan_{knob.id}", title or knob.id, (knob, *first), single))

    if not pages:
        # neither a dose nor a disease nor a scan: the readouts as the model is
        pages.append(("model", "Model", (), single))
    panels: dict[str, Panel] = {}
    app_pages = []
    for page, title, knobs, simulation in pages:
        ids_of_page = []
        for k, readout in enumerate(readouts):
            panel_id = f"{page}_{k}"
            panels[panel_id] = Panel(
                id=panel_id,
                title=labels[readout],
                knobs=knobs,
                curves=(
                    PanelCurve(id=readout, label=labels[readout], unit=units[readout]),
                ),
                simulation=simulation,
            )
            ids_of_page.append(panel_id)
        app_pages.append(
            AppPage(id=page, title=title, panels=tuple(ids_of_page), size="small")
        )
    variants = tuple(
        AppVariant(id=name, label=name.replace("_", " ").capitalize())
        for name in model.variants
    )
    app = App(
        title=model.title or model.name,
        subtitle="the readouts of the model under its doses, diseases and scans",
        pages=tuple(app_pages),
        variants=variants if len(variants) > 1 else (),
        citation=_citation(model),
        doi=(model.release.doi or "") if model.release is not None else "",
        licenses=("Models: CC BY-SA 4.0", "Source code: MIT"),
    )
    return panels, app


def default_variant_panels(
    model: PKPDModel, panels: Mapping[str, Panel], variants: Mapping[str, Path]
) -> dict[str, dict[str, Panel]]:
    """Get the default panels of the builds of the other sexes, see the module.

    The panels of `default_app` are the ones of the default build. The classes
    of a renal impairment are the renal functions at the healthy GFR of a
    build (`ScaledRenalMap`), so the build of another sex gets the panels of
    the page of the disease with the classes at its own healthy GFR; every
    other panel is the same in every build.

    Args:
        model: the model.
        panels: the panels of `default_app` by their id.
        variants: the SBML file of every variant, the default first.

    Returns:
        The panels of a variant by its name and the id of the panel, only the
        panels which differ from the ones of the default build.
    """
    names = list(variants)
    if len(names) < 2:
        return {}
    base = _base(_classes(model))
    default = _sex(names[0])
    keys = {disease.key: disease for disease in DISEASES}
    shown: dict[str, dict[str, Panel]] = {}
    for name in names[1:]:
        sex = _sex(name)
        if sex is None or sex == default:
            continue
        changed: dict[str, Panel] = {}
        for panel_id, panel in panels.items():
            knobs = list(panel.knobs)
            for k, knob in enumerate(knobs):
                disease = keys.get(knob.id)
                if disease is None or knob.kind != "choice" or knob.target is None:
                    continue
                other = _disease_knob(disease, knob.target, base, sex)
                if other is not None and other != knob:
                    knobs[k] = other
            if knobs != list(panel.knobs):
                changed[panel_id] = panel.model_copy(update={"knobs": tuple(knobs)})
        if changed:
            shown[name] = changed
    return shown


def _entities(sbml: Path) -> set[str]:
    """Get the ids a panel sets or selects in a model: its entities and `[S]`."""
    doc: libsbml.SBMLDocument = libsbml.readSBMLFromFile(str(sbml))
    model: libsbml.Model = doc.getModel()
    ids = {p.getId() for p in model.getListOfParameters()}
    ids |= {c.getId() for c in model.getListOfCompartments()}
    for species in model.getListOfSpecies():
        ids.add(species.getId())
        ids.add(f"[{species.getId()}]")
    return ids


def _classes(model: PKPDModel) -> list[type[SimulationExperiment]]:
    """Get the classes of the experiments of a model, none for a model without."""
    try:
        return list(model.experiment_classes())
    except ValueError:
        # a model without simulation experiments
        return []


def _base(classes: Iterable[type]) -> type | None:
    """Get the base of the experiments of a model, whose maps give the classes."""
    from pkpdsim.experiment import PKPDSimulationExperiment

    return next((c for c in classes if issubclass(c, PKPDSimulationExperiment)), None)


def _sex(variant: str) -> Sex | None:
    """Get the sex of the build of a variant, `None` for a variant of no sex."""
    try:
        return Sex(variant)
    except ValueError:
        return None


def _studied(model: PKPDModel, ids: set[str]) -> _Studied:
    """Collect what the experiments with data of a model state, see `_Studied`."""
    from sbmlsim.fit.helpers import FitMappings

    classes = _classes(model)
    base = _base(classes)
    scans = tuple(
        spec
        for c in classes
        # a map of a class, not a property of its instances
        if isinstance(scan_map := getattr(c, "scan_map", None), Mapping)
        for spec in scan_map.values()
        if isinstance(spec, Mapping) and "parameter" in spec
    )
    try:
        mappings = FitMappings.for_project(model, classes)
    except Exception:  # an experiment which cannot be read loses its readouts only
        logger.warning(
            "The experiments of the model '%s' cannot be read, its app has no "
            "readouts of them.",
            model.name,
            exc_info=True,
        )
        return _Studied((), {}, (), frozenset(), scans, base)
    readouts: Counter[str] = Counter()
    health: set[str] = set()
    doses: dict[str, Counter[float]] = {}
    ends: list[float] = []
    for experiment in mappings.experiments.values():
        fit_mappings = getattr(experiment, "_fit_mappings", {}) or {}
        if not fit_mappings:
            continue
        for mapping in fit_mappings.values():
            # the selection as given, `[S]` the concentration and `S` the amount
            selection = getattr(
                getattr(mapping.observable, "y", None), "selection", None
            )
            if mapping.observable.task_id is not None and isinstance(selection, str):
                readouts[selection] += 1
            value = getattr(getattr(mapping, "metadata", None), "health", None)
            if value is not None:
                health.add(str(getattr(value, "value", value)))
        for simulation in _simulations(experiment):
            dosed = False
            for target, value in (simulation.preinit_changes or {}).items():
                if not DOSE_TARGETS.fullmatch(target) or target not in ids:
                    continue
                mg = _milligram(value)
                if mg is not None and mg > 0.0:
                    doses.setdefault(target, Counter())[round(mg, 6)] += 1
                    dosed = True
            if dosed:
                ends.append(_hours(simulation.end, simulation.time_unit))
    ordered = sorted(doses, key=lambda t: -sum(doses[t].values()))
    return _Studied(
        readouts=tuple(r for r, _ in readouts.most_common()),
        doses={t: tuple(d for d, _ in doses[t].most_common()) for t in ordered},
        ends=tuple(ends),
        health=frozenset(health),
        scans=scans,
        base=base,
    )


def _simulations(experiment: SimulationExperiment) -> Iterable[Any]:
    """The simulations of an experiment, the one of a scan for a scan."""
    for simulation in (getattr(experiment, "_simulations", None) or {}).values():
        yield getattr(simulation, "simulation", simulation)


def _milligram(value: object) -> float | None:
    """A dose in mg, a number of the model taken as mg; `None` for no number."""
    if isinstance(value, Quantity):
        try:
            return float(value.to("mg").magnitude)
        except PintError:
            return None
    if isinstance(value, (int, float)):
        return float(value)
    return None


def _hours(end: float, unit: object) -> float:
    """An end of a simulation in hours."""
    try:
        return float(
            ureg.Quantity(end, str(unit) if unit else "min").to("hr").magnitude
        )
    except PintError, ValueError, TypeError:
        return float(end) / 60.0


def _time_unit(sbml: Path) -> str:
    """The time unit of a model, from the SBML the factory wrote."""
    from sbmlsim.units import UnitsInformation

    uinfo = UnitsInformation.from_sbml(sbml)
    return str(uinfo.get("time") or "min")


def _end(ends: Sequence[float], time_unit: str) -> float:
    """The end of a page of single doses: the median end of the studies, 24 h without."""
    hours = statistics.median(ends) if ends else 24.0
    hours = max(hours, 1.0)
    end = float(ureg.Quantity(hours, "hr").to(time_unit).magnitude)
    return min(end, MAX_END)


def _units(sbml: Path, readouts: Sequence[str]) -> dict[str, str | None]:
    """The units the readouts are drawn in: uM for a concentration, umol for an amount."""
    from sbmlsim.units import UnitsInformation

    uinfo = UnitsInformation.from_sbml(sbml)
    units: dict[str, str | None] = {}
    for readout in readouts:
        unit = uinfo.get(readout)
        shown: str | None = None
        if unit:
            try:
                quantity = ureg.Quantity(1.0, unit)
                if quantity.check("[substance]/[length]**3"):
                    shown = "uM"
                elif quantity.check("[substance]"):
                    shown = "umol"
            except PintError, ValueError, TypeError:
                shown = None
        units[readout] = shown
    return units


def _label(sbml: Path, selection: str) -> str:
    """The name of the entity of a selection, its id without one."""
    doc: libsbml.SBMLDocument = libsbml.readSBMLFromFile(str(sbml))
    sid = selection.strip("[]")
    element = doc.getModel().getElementBySId(sid)
    name = element.getName() if element is not None else ""
    # the unit a name of the library carries, `... [mM]`, is the one of the axis
    return re.sub(r"\s*\[[^\]]*\]$", "", name) or sid


def _substance_names(sbml: Path) -> dict[str, str]:
    """The names of the substances of a model by their id, from its species.

    A species `<site>_<substance>` of the PBPK library is named `<substance
    name> (<site>)`, e.g. `Aurine_edo` `edoxaban (urine)`, so the name of a
    substance is the name of such a species before its site.
    """
    doc: libsbml.SBMLDocument = libsbml.readSBMLFromFile(str(sbml))
    names: dict[str, str] = {}
    for species in doc.getModel().getListOfSpecies():
        sid, name = species.getId(), species.getName()
        if "_" not in sid or " (" not in name:
            continue
        substance = sid.split("_", 1)[1]
        names.setdefault(substance, name.split(" (", 1)[0])
    return names


def _dose_knobs(
    doses: Mapping[str, Sequence[float]],
    ids: set[str],
    names: Mapping[str, str],
    substances: Mapping[str, str],
) -> list[Knob]:
    """The dose knobs of the most studied routes, the first one dosed by default.

    A knob is labelled by the name of its target in the model, e.g. `oral dose
    caffeine`, without its unit, which the knob states, else by its route and
    the name of its substance (`_substance_names`).
    """
    knobs = []
    for k, (target, values) in enumerate(doses.items()):
        if target not in ids or len(knobs) == DOSE_KNOBS:
            continue
        match = DOSE_TARGETS.fullmatch(target)
        assert match is not None
        route, substance = match.groups()
        largest = max(values)
        stops = tuple(
            Choice(label=f"{_number(v)} mg", value=v) for v in sorted(set(values))
        )
        substance_name = substances.get(substance, substance)
        label = names.get(target) or ""
        if not label or label == target:
            label = f"{ROUTE_LABELS.get(route, route)} dose {substance_name}"
        # a name of the id of its substance, `oral dose edo`, names the substance
        label = re.sub(rf"\b{re.escape(substance)}\b", substance_name, label)
        knobs.append(
            Knob(
                id=target,
                label=re.sub(r"\s*\[mg\]$", "", label),
                target=target,
                unit="mg",
                kind="range",
                min=0.0,
                max=largest,
                steps=5,
                choices=stops,
                default=values[0] if k == 0 else 0.0,
            )
        )
    return knobs


def _number(value: float) -> str:
    """A dose as it is shown, without a trailing zero."""
    return f"{value:g}"


def _protocol_page(
    dose_knobs: Sequence[Knob], zero: Mapping[str, float], time_unit: str
) -> tuple[str, str, tuple[Knob, ...], PanelSimulation]:
    """The page of the dosing protocols: the dose knobs repeated at an interval."""
    interval = Knob(
        id="interval",
        label="dosing",
        unit="hr",
        kind="choice",
        choices=tuple(Choice(label=label, value=value) for label, value in INTERVALS),
        default=24.0,
    )
    end = min(
        float(ureg.Quantity(PROTOCOL_DAYS, "day").to(time_unit).magnitude), MAX_END
    )
    # whole days, robust against the rounding of the conversion (6 days are 8640
    # min, which convert back to 6.000000000000001 or 5.999999999999999 days)
    days_max = math.floor(
        float(ureg.Quantity(end, time_unit).to("day").magnitude) + 1e-9
    )
    days = Knob(
        id="days",
        label="days of dosing",
        unit="day",
        kind="range",
        min=1.0,
        max=float(days_max),
        steps=days_max,
        default=float(days_max),
    )
    simulation = PanelSimulation(
        end=end,
        points=PROTOCOL_POINTS,
        x_unit="day",
        preinit=dict(zero),
        dosing=tuple(
            PanelDosing(
                target=str(knob.target),
                dose=knob.id,
                interval="interval",
                duration="days",
            )
            for knob in dose_knobs
        ),
    )
    return ("protocol", "Dosing protocol", (*dose_knobs, interval, days), simulation)


def _studies(studied: _Studied, disease: _Disease, target: str) -> bool:
    """Whether the model is studied in a disease: subjects with it or a scan of it."""
    if any(spec.get("parameter") == target for spec in studied.scans):
        return True
    return any(
        word in health.lower() for health in studied.health for word in disease.words
    )


def _disease_knob(
    disease: _Disease, target: str, base: type | None, sex: Sex | None
) -> Knob | None:
    """The choice knob of the classes of a disease, the first class its default.

    The classes are the map of the base of the experiments of the model
    (`PKPDSimulationExperiment.renal_map` and its siblings), the classes its
    experiments simulate, e.g. the renal functions of the classes of the GFR at
    the healthy GFR of the model; two classes of one value are one choice. A
    map of the classes of a sex (`ScaledRenalMap.of_sex`) gives the classes of
    the build of `sex`, at its own healthy GFR.
    """
    classes = getattr(base, disease.map, None) if base is not None else None
    of_sex = getattr(classes, "of_sex", None)
    if sex is not None and callable(of_sex):
        # a map of one build raises, its classes are the ones of every sex
        with contextlib.suppress(TypeError):
            classes = of_sex(sex)
    if not classes:
        return None
    choices: list[Choice] = []
    for label, value in classes.items():
        number = float(getattr(value, "magnitude", value))
        if any(math.isclose(number, other.value) for other in choices):
            continue
        choices.append(Choice(label=str(label), value=number))
    if len(choices) < 2:
        return None
    return Knob(
        id=disease.key,
        label=disease.label,
        target=target,
        unit="dimensionless",
        kind="choice",
        choices=tuple(choices),
        default=choices[0].value,
    )


def _scan_knob(
    spec: Mapping[str, Any], ids: set[str], covered: set[str]
) -> Knob | None:
    """The range knob of a scan which no other page covers, `None` for one it does."""
    target = str(spec.get("parameter", ""))
    if target not in ids or target in covered or DOSE_TARGETS.fullmatch(target):
        return None
    values = [float(v) for v in spec.get("range", ())]
    if len(values) < 2:
        return None
    low, high = min(values), max(values)
    default = float(spec.get("default", low))
    if not low <= default <= high or low == high:
        return None
    unit = str(spec.get("units") or "dimensionless")
    try:
        ureg.Quantity(1.0, unit)
    except PintError, ValueError, TypeError:
        unit = "dimensionless"
    key = re.sub(r"\W", "_", target)
    return Knob(
        id=key,
        label=re.sub(r"\s*\[[^\]]*\]$", "", str(spec.get("label") or target)) or target,
        target=target,
        unit=unit,
        kind="range",
        min=low,
        max=high,
        steps=min(len(set(values)), 9),
        default=default,
    )


def _citation(model: PKPDModel) -> str:
    """The citation of the release of the model, empty without a release."""
    if model.release is None:
        return ""
    return model.release.citation(markdown=False)
