"""The base of the simulation experiments of a model built on the library.

A model of the library has the impairments as the severities and the renal
function of the library (`CP_severity`, `HF_severity`, `f_renal_function`, see
`pkpdsim.pbpk.impairment`), so the classes of a study map onto other values than
those of the hand-written models. The base of its experiments derives from
`PBPKSimulationExperiment` instead of `PKPDSimulationExperiment` and binds the
model object as usual:

```python
# not executed
class GlimepirideSimulationExperiment(PBPKSimulationExperiment):
    pkpd_model = MODEL
```

The GFR of the library is `GFR = GFR_base * f_renal_function`, and `GFR_base`,
the GFR of the individual at renal function 1, scales with its body surface
area (and falls with heart failure and cirrhosis). `healthy_gfr` reads it from
the flat model, `PBPKSimulationExperiment.renal_function` gives the renal
function of a measured clearance (the simulated GFR is the measured one), and
the classes of renal impairment keep their GFR in ml/min (`ScaledRenalMap`).

The outflow of the kidney is valid while the net filtered plasma flow of a
substance is below the plasma flow of the kidney, `fu*GFR*(1 - f_reabsorbed)
< Qki*(1 - HCT)` (`Qki*BP` with erythrocytes in equilibrium, `fu*GFR` with a
saturable reabsorption): the plasma leaving the kidney lacks the net
filtrate. GFR, flow and hematocrit are parameters of the individual, so
`kidney_outflow` evaluates it for an individual (the changes of a simulation
before its initialization), and `PBPKSimulationExperiment` refuses an
experiment which simulates an individual outside of it and a renal function
beyond it (`renal_function`, `renal_function_limits`).

The module imports sbmlsim; `pkpdsim.pbpk` exports its names lazily, so
`import pkpdsim.pbpk` (a factory) imports no simulator.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Callable, Iterable, Iterator, Mapping
from dataclasses import dataclass
from functools import cache, lru_cache
from pathlib import Path
from typing import ClassVar, override

import libsbml
from sbmlsim import Q
from sbmlsim.model.model_roadrunner import RoadrunnerSBMLModel
from sbmlsim.simulation import Scan, Simulation
from sbmlsim.simulation.definition import Value as Change
from sbmlsim.simulation.scan import DimensionKind
from sbmlsim.simulator import Simulator
from sbmlsim.units import Quantity

from pkpdsim import impairment
from pkpdsim.experiment import PKPDSimulationExperiment
from pkpdsim.pbpk import ids
from pkpdsim.pbpk.organs import Organ, Part
from pkpdsim.scan import start_values
from pkpdsim.uncertainty import UNCERTAINTY_DIMENSION

#: the healthy GFR the presets of `pkpdsim.impairment` are relative to [ml/min];
#: import it from `pkpdsim.pbpk` (`from pkpdsim.pbpk import REFERENCE_GFR`), the
#: one path of the names of this module
REFERENCE_GFR = 101.0


def healthy_gfr(
    model_path: Path,
    *,
    bw: float | None = None,
    height: float | None = None,
    hf_severity: float = 0.0,
    cp_severity: float = 0.0,
) -> Quantity:
    """Get the healthy glomerular filtration rate of an individual (ml/min).

    The GFR of an individual at renal function 1 (`GFR_base` of the library,
    `GFR = GFR_base * f_renal_function`): it scales with the body surface area
    of the body weight and the height, the severity of heart failure lowers
    the renal blood flow and with it the GFR, and the cirrhosis lowers it by
    the hepatorenal factor of its Child-Pugh class. The GFR and its covariates
    are initial assignments, so it is read from a simulation of one minute
    with the covariates as changes before the initialization, which evaluates
    the initial assignments again (sbmlsim). It is cached per model file, its
    modification time and the individual, once for all experiments of a
    process; a model file which is created again is simulated again. The sex
    is no parameter of the model (the library builds the model of one sex).

    Args:
        model_path: the flat SBML model of the library.
        bw: the body weight of the individual [kg], `None` for that of the
            model.
        height: the height of the individual [cm], `None` for that of the
            model.
        hf_severity: the severity of heart failure `HF_severity` (0 to 4).
        cp_severity: the severity of cirrhosis `CP_severity` (0 to 3).

    Returns:
        The GFR at renal function 1 in ml/min, a quantity as the renal
        function of `PBPKSimulationExperiment.renal_function`.
    """
    path = Path(model_path)
    gfr = _healthy_gfr(
        path,
        path.stat().st_mtime_ns,
        bw,
        height,
        float(hf_severity),
        float(cp_severity),
    )
    return Q(gfr, "ml/min")


def _healthy_gfr(
    model_path: Path,
    mtime_ns: int,
    bw: float | None,
    height: float | None,
    hf_severity: float,
    cp_severity: float,
) -> float:
    """Simulate the healthy GFR of an individual once, see `healthy_gfr`.

    `mtime_ns`, the modification time of the model file, is a key of the
    cache: a model which is created again is simulated again.
    """
    return _individual(
        model_path, mtime_ns, _individual_changes(bw, height, hf_severity, cp_severity)
    ).ml_per_min(ids.GFR_BASE)


def _individual_changes(
    bw: float | None,
    height: float | None,
    hf_severity: float,
    cp_severity: float,
    renal_function: float | None = None,
) -> dict[str, Quantity]:
    """The changes of an individual before the initialization, by its covariates."""
    changes = {
        ids.HF_SEVERITY: Q(float(hf_severity), "dimensionless"),
        ids.CP_SEVERITY: Q(float(cp_severity), "dimensionless"),
    }
    if bw is not None:
        changes[ids.BW] = Q(bw, "kg")
    if height is not None:
        changes[ids.HEIGHT] = Q(height, "cm")
    if renal_function is not None:
        changes[ids.F_RENAL_FUNCTION] = Q(float(renal_function), "dimensionless")
    return changes


@dataclass(frozen=True)
class KidneyOutflow:
    """The net filtered plasma flow of a substance and the plasma flow of the kidney.

    The outflow of the kidney is valid while the net filtered plasma flow is
    below the plasma flow which carries the substance through the kidney.

    Attributes:
        substance: the substance, filtered into the urine.
        gfr: the GFR of the individual [ml/min].
        filtered: the net filtered plasma flow `fu*GFR*(1 - f_reabsorbed)`,
            `fu*GFR` with a saturable reabsorption (its fraction falls to 0
            at a high concentration) [ml/min].
        plasma_flow: the flow which carries the substance through the kidney,
            `Qki*(1 - HCT)`, `Qki*BP` with erythrocytes in equilibrium
            [ml/min].
    """

    substance: str
    gfr: float
    filtered: float
    plasma_flow: float

    @property
    def valid(self) -> bool:
        """Whether the net filtered plasma flow is below the plasma flow."""
        return self.filtered < self.plasma_flow

    def text(self) -> str:
        """Describe the flows, for a message."""
        return (
            f"the net filtered plasma flow of '{self.substance}' "
            f"fu*GFR*(1 - f_reabsorbed) {self.filtered:.4g} ml/min at the GFR "
            f"{self.gfr:.4g} ml/min against the plasma flow of the kidney "
            f"{self.plasma_flow:.4g} ml/min (Qki*(1 - HCT), Qki*BP with "
            "erythrocytes in equilibrium)"
        )


def kidney_outflow(
    model_path: Path, changes: Mapping[str, Change] | None = None
) -> tuple[KidneyOutflow, ...]:
    """Get the net filtered plasma flows and the plasma flow of the kidney of an individual.

    One `KidneyOutflow` for every substance which the kidney filters (its
    filtration reaction is in the model). The flows are initial assignments,
    so they are read from a simulation of one minute with the changes before
    the initialization, which evaluates the initial assignments again, as
    `healthy_gfr` reads the GFR. Only the changes of the parameters which the
    flows depend on make the individual (a dose does not), and the values are
    cached per model file, its modification time and the individual; the
    model is loaded once per process.

    Args:
        model_path: the flat SBML model of the library.
        changes: the changes of the individual before the initialization
            (`Simulation.preinit_changes`), numbers in the units of the model
            or quantities.

    Returns:
        The flows of every filtered substance, in the order of the model.
    """
    path = Path(model_path)
    mtime = path.stat().st_mtime_ns
    model = _kidney_model(path, mtime)
    values = _individual(path, mtime, _relevant(changes or {}, model.inputs))
    flow = values.ml_per_min(ids.organ_flow(Organ.KIDNEY))
    hct = values.number(ids.HCT)
    gfr = values.ml_per_min(ids.GFR)
    outflows = []
    for s, kinetic in model.filtered:
        reabsorbed = values.numbers.get(ids.fraction_reabsorbed(s), 0.0)
        carrier = 1.0 - hct if kinetic else values.number(ids.bp(s))
        outflows.append(
            KidneyOutflow(
                substance=s,
                gfr=gfr,
                filtered=values.number(ids.fu(s)) * gfr * (1.0 - reabsorbed),
                plasma_flow=flow * carrier,
            )
        )
    return tuple(outflows)


@dataclass(frozen=True)
class _KidneyModel:
    """What the validity of the outflow reads from a flat model.

    Attributes:
        filtered: every filtered substance with whether its erythrocytes are
            kinetic (the plasma carries it, else `BP` does, which is
            `1 - HCT` without erythrocytes).
        quantities: the ids of the flows of the validity and `GFR_base`
            which the model has, the selections of its simulations.
        inputs: the ids which the quantities depend on, through initial
            assignments and assignment rules (themselves included): the
            changes of an individual.
    """

    filtered: tuple[tuple[str, bool], ...]
    quantities: tuple[str, ...]
    inputs: frozenset[str]


@lru_cache(maxsize=8)
def _kidney_model(model_path: Path, mtime_ns: int) -> _KidneyModel:
    """Read the filtered substances and the inputs of the validity from a flat model.

    `mtime_ns`, the modification time of the model file, is only a key of the
    cache.
    """
    model = libsbml.readSBMLFromFile(str(model_path)).getModel()
    prefix = ids.flat(Organ.KIDNEY, ids.filtration_reaction(""))
    species = {sp.getId() for sp in model.getListOfSpecies()}
    filtered = tuple(
        (s, ids.flat(Organ.KIDNEY, ids.species(Part.BLOOD_CELLS, s)) in species)
        for s in (
            reaction.getId().removeprefix(prefix)
            for reaction in model.getListOfReactions()
            if reaction.getId().startswith(prefix)
        )
    )
    quantities = [ids.GFR, ids.GFR_BASE, ids.organ_flow(Organ.KIDNEY), ids.HCT]
    for s, _ in filtered:
        quantities += [ids.fu(s), ids.bp(s), ids.fraction_reabsorbed(s)]
    present = tuple(q for q in quantities if model.getParameter(q) is not None)
    return _KidneyModel(filtered, present, _dependencies(model, present))


def _dependencies(model: libsbml.Model, roots: Iterable[str]) -> frozenset[str]:
    """The ids which the roots depend on through initial assignments and rules."""
    reads: dict[str, set[str]] = {}
    for assignment in model.getListOfInitialAssignments():
        reads.setdefault(assignment.getSymbol(), set()).update(
            _names(assignment.getMath())
        )
    for rule in model.getListOfRules():
        if rule.isAssignment():
            reads.setdefault(rule.getVariable(), set()).update(_names(rule.getMath()))
    seen: set[str] = set()
    stack = list(roots)
    while stack:
        sid = stack.pop()
        if sid not in seen:
            seen.add(sid)
            stack.extend(reads.get(sid, ()))
    return frozenset(seen)


def _names(math: libsbml.ASTNode | None) -> set[str]:
    """The names a formula reads."""
    names: set[str] = set()
    stack = [] if math is None else [math]
    while stack:
        node = stack.pop()
        if node.getType() == libsbml.AST_NAME:
            names.add(node.getName())
        stack.extend(node.getChild(k) for k in range(node.getNumChildren()))
    return names


#: a change of an individual, hashable: the target and its value, a number, a
#: formula or the magnitude and the unit of a quantity
_Frozen = tuple[str, float | str | tuple[float, str]]


def _relevant(
    changes: Mapping[str, Change], inputs: frozenset[str]
) -> dict[str, Change]:
    """The changes of the ids which the validity of the outflow depends on."""
    return {k: v for k, v in changes.items() if k in inputs}


@dataclass(frozen=True)
class _Values:
    """The start values of an individual, in the units of the model.

    Attributes:
        numbers: id -> the value at the start.
        units: id -> its unit.
    """

    numbers: Mapping[str, float]
    units: Mapping[str, str]

    def number(self, sid: str) -> float:
        """The value of an id in the unit of the model."""
        return self.numbers[sid]

    def ml_per_min(self, sid: str) -> float:
        """The value of a flow in ml/min."""
        return float(Q(self.numbers[sid], self.units[sid]).to("ml/min").magnitude)


def _individual(
    model_path: Path, mtime_ns: int, changes: Mapping[str, Change]
) -> _Values:
    """Simulate an individual once, see `kidney_outflow`, cached by its changes."""
    return _individual_values(model_path, mtime_ns, _frozen(changes))


def _frozen(changes: Mapping[str, Change]) -> tuple[_Frozen, ...]:
    """The changes as a key of the cache, sorted by target."""
    frozen: list[_Frozen] = []
    for target, v in sorted(changes.items()):
        if isinstance(v, str):
            frozen.append((target, v))
        elif isinstance(v, Quantity):
            frozen.append((target, (float(v.magnitude), str(v.units))))
        else:
            frozen.append((target, float(v)))
    return tuple(frozen)


def _thawed(changes: tuple[_Frozen, ...]) -> dict[str, Change]:
    """The changes of a key of the cache."""
    return {
        target: Q(v[0], v[1]) if isinstance(v, tuple) else v for target, v in changes
    }


@cache
def _individual_values(
    model_path: Path, mtime_ns: int, changes: tuple[_Frozen, ...]
) -> _Values:
    """Simulate an individual for one minute and read its start values, once."""
    quantities = _kidney_model(model_path, mtime_ns).quantities
    result = Simulator(n_workers=1).run(
        _loaded(model_path, mtime_ns),
        Simulation(end=1, time_unit="min", preinit_changes=_thawed(changes)),
    )
    numbers = {sid: float(start_values(result, sid).item()) for sid in quantities}
    return _Values(numbers, {sid: str(result.units[sid]) for sid in quantities})


@lru_cache(maxsize=2)
def _loaded(model_path: Path, mtime_ns: int) -> RoadrunnerSBMLModel:
    """Load a model once per process for the simulations of individuals.

    It selects only the quantities of the validity of the kidney outflow and
    `GFR_base`, and it is compiled without optimization (`O0`), the fastest
    compilation, since its simulations of one minute only evaluate the
    initial assignments; the two models used last stay loaded.
    """
    quantities = _kidney_model(model_path, mtime_ns).quantities
    return RoadrunnerSBMLModel(
        source=model_path, selections=["time", *quantities], jit_optimization="O0"
    )


class ScaledRenalMap(Mapping[str, float]):
    """The renal function of the classes of renal impairment at a healthy GFR.

    A preset of `pkpdsim.impairment` gives the renal function of a class
    relative to a healthy GFR of `reference` ml/min (101); the classes keep
    their GFR in ml/min, so a value is `preset[class] * reference / gfr()`.
    The GFR is evaluated on the first access to a value, so importing the
    experiments loads no model; the keys are those of the preset at once.

    Without `gfr` the map is bound to the class it is an attribute of: a
    class attribute is a (non-data) descriptor, and the access through a
    class (or an instance) gives the map of that class, on its classmethod
    `healthy_gfr`, one map per class, created on the first access. So the base
    of the experiments of every model, and every subclass with another model
    object, gets the GFR of its own model, and a model overrides the preset in
    its base:

    ```python
    # not executed
    class HCTZSimulationExperiment(PBPKSimulationExperiment):
        pkpd_model = MODEL
        renal_map: ClassVar[Mapping[str, float]] = ScaledRenalMap(
            impairment.RENAL_FUNCTION
        )
    ```
    """

    def __init__(
        self,
        preset: Mapping[str, float],
        gfr: Callable[[], Quantity] | None = None,
        reference: float = REFERENCE_GFR,
    ) -> None:
        """Get the classes of a preset at a healthy GFR.

        Args:
            preset: the renal function of the classes relative to `reference`.
            gfr: the healthy GFR (a quantity of a flow), called on the first
                access to a value; `None` for the `healthy_gfr` of the class the
                map is an attribute of.
            reference: the healthy GFR the preset is relative to [ml/min].
        """
        self._preset = dict(preset)
        self._gfr = gfr
        self._reference = reference
        self._values: dict[str, float] | None = None
        #: the maps bound to the classes which accessed this unbound map; a
        #: plain dict, since a bound map holds the classmethod `healthy_gfr` of
        #: its class, which a weak key would not free either (the classes of
        #: experiments live as long as their module), a cache which `pickle`
        #: leaves out (`__getstate__`)
        self._bound: dict[type, ScaledRenalMap] = {}

    def __getstate__(self) -> dict[str, object]:
        """Get the state to pickle: without the maps bound to classes."""
        return {**self.__dict__, "_bound": {}}

    def __get__(
        self, instance: object | None, owner: type | None = None
    ) -> Mapping[str, float]:
        """Get the map of the class of an access, see the class.

        Args:
            instance: the instance of the access, `None` for a class.
            owner: the class of the access.

        Returns:
            This map if it has its own GFR (or the class has no `healthy_gfr`),
            else the map of the class on its `healthy_gfr`.
        """
        if owner is None:
            owner = type(instance)
        gfr = getattr(owner, "healthy_gfr", None)
        if self._gfr is not None or gfr is None:
            return self
        bound = self._bound.get(owner)
        if bound is None:
            bound = ScaledRenalMap(self._preset, gfr, self._reference)
            self._bound[owner] = bound
        return bound

    def _scaled(self) -> dict[str, float]:
        """Get the renal functions, evaluating the healthy GFR once.

        Raises:
            TypeError: for a map without a GFR which no class bound.
        """
        if self._values is None:
            if self._gfr is None:
                raise TypeError(
                    "The map has no healthy GFR: give `gfr` or make it an "
                    "attribute of a class with the classmethod `healthy_gfr`."
                )
            factor = self._reference / float(self._gfr().m_as("ml/min"))
            self._values = {key: v * factor for key, v in self._preset.items()}
        return self._values

    @override
    def __getitem__(self, key: str) -> float:
        return self._scaled()[key]

    @override
    def __iter__(self) -> Iterator[str]:
        return iter(self._preset)

    @override
    def __len__(self) -> int:
        return len(self._preset)

    @override
    def __contains__(self, key: object) -> bool:
        return key in self._preset

    @override
    def __repr__(self) -> str:
        return f"{type(self).__name__}({self._preset!r}, reference={self._reference})"


class PBPKSimulationExperiment(PKPDSimulationExperiment):
    """The base of the simulation experiments of a model built on the library.

    The maps of the classes of a study onto the library: the renal function of
    the classes of `impairment.RENAL_FUNCTION_GFR` (the limits of the GFR and
    the end stage renal disease) at the healthy GFR of the model, the
    Child-Pugh classes as `CP_severity` (`impairment.CHILD_PUGH`) and the NYHA
    classes as `HF_severity` (`impairment.HEART_FAILURE`), with their colors.
    The classes of renal impairment keep their GFR in ml/min: the renal
    function of a class is its GFR over the healthy GFR of the model
    (`ScaledRenalMap`), unlike the eGFR of the library, which sets the factors
    of chronic kidney disease and is relative to the GFR per 1.73 m2 of the
    reference individual (`pkpdsim.pbpk.impairment`). A model overrides any
    attribute, e.g. the classes by the midpoints of the GFR
    (`ScaledRenalMap(impairment.RENAL_FUNCTION)`), its own map of the heart
    failure or a narrower `renal_function_range`.

    The initialization refuses an experiment which simulates an individual
    outside the validity of the kidney outflow (`kidney_outflow`): the net
    filtered plasma flow of a substance must stay below the plasma flow of
    the kidney, at the start of every simulation and every point of a scan
    (the draws of the fitted parameters of a run over them are not checked).

    Attributes:
        renal_function_range: a narrower range of the renal functions the
            model supports, `None` (the default) for the range of the validity
            of the kidney outflow, from 0 (anuria, the substances leave by the
            other routes alone) to the renal function at which the net
            filtered plasma flow of a substance reaches the plasma flow of
            the kidney (`renal_function_limits`); a model whose urine is not
            defined without filtration starts above 0. The validity holds
            within the range as well.
    """

    renal_map: ClassVar[Mapping[str, float]] = ScaledRenalMap(
        impairment.RENAL_FUNCTION_GFR
    )
    renal_colors: ClassVar[Mapping[str, str]] = impairment.RENAL_COLORS_END_STAGE
    cirrhosis_map: ClassVar[Mapping[str, float]] = impairment.CHILD_PUGH
    cirrhosis_colors: ClassVar[Mapping[str, str]] = impairment.CIRRHOSIS_COLORS
    cardiac_map: ClassVar[Mapping[str, float]] = impairment.HEART_FAILURE
    cardiac_colors: ClassVar[Mapping[str, str]] = impairment.HEART_FAILURE_COLORS
    renal_function_range: ClassVar[tuple[float, float] | None] = None

    @override
    def initialize(self) -> None:
        """Initialize the experiment and check the kidney outflow of its individuals.

        Raises:
            ValueError: for a simulation of an individual outside the validity
                of the kidney outflow, naming the experiment, the simulation,
                the GFR and the plasma flow of the kidney.
        """
        super().initialize()
        try:
            self._check_kidney_outflow()
        except ValueError as err:
            err.add_note(f"The experiment '{self.sid}' cannot be initialized.")
            raise

    def _check_kidney_outflow(self) -> None:
        """Check the kidney outflow of every individual the experiment simulates.

        Raises:
            ValueError: for an individual outside the validity.
        """
        path = Path(self.pkpd_model.model_path)
        inputs = _kidney_model(path, path.stat().st_mtime_ns).inputs
        for key, simulation in self._simulations.items():
            for changes in _individuals(simulation, inputs):
                for outflow in kidney_outflow(path, changes):
                    if not outflow.valid:
                        raise ValueError(
                            f"The experiment '{self.sid}' simulates in '{key}' an "
                            "individual outside the validity of the kidney "
                            f"outflow: {outflow.text()}; the net filtered plasma "
                            "flow must stay below the plasma flow. The changes of "
                            f"the individual: {_changes_text(changes)}."
                        )

    @classmethod
    def healthy_gfr(
        cls,
        *,
        bw: float | None = None,
        height: float | None = None,
        hf_severity: float = 0.0,
        cp_severity: float = 0.0,
    ) -> Quantity:
        """Get the healthy glomerular filtration rate of an individual (ml/min).

        `healthy_gfr` of the module for the model of the model object, cached
        for all experiments.

        Args:
            bw: the body weight of the individual [kg], `None` for that of the
                model.
            height: the height of the individual [cm], `None` for that of the
                model.
            hf_severity: the severity of heart failure `HF_severity` (0 to 4).
            cp_severity: the severity of cirrhosis `CP_severity` (0 to 3).

        Returns:
            The GFR at renal function 1 in ml/min, a quantity.

        Raises:
            TypeError: if the class has no model object (`pkpd_model`).
        """
        if not hasattr(cls, "pkpd_model"):
            raise TypeError(
                f"{cls.__name__} has no model object: set `pkpd_model` in the "
                "base of the experiments of the model."
            )
        return healthy_gfr(
            cls.pkpd_model.model_path,
            bw=bw,
            height=height,
            hf_severity=hf_severity,
            cp_severity=cp_severity,
        )

    @classmethod
    def renal_function(
        cls,
        clearance: float,
        *,
        bw: float | None = None,
        height: float | None = None,
        hf_severity: float = 0.0,
        cp_severity: float = 0.0,
    ) -> Quantity:
        """Get the renal function of a measured clearance.

        The glomerular filtration rate of the model is `f_renal_function *
        GFR_base`, and `GFR_base`, the healthy rate of the individual of the
        simulation (`healthy_gfr`), scales with its body surface area and
        falls with the severity of heart failure and of cirrhosis. A measured
        clearance includes all of them, so the renal function is the
        clearance relative to the healthy rate of the individual of the
        study, and the simulated rate is the measured one; the simulation
        sets the same covariates.

        Args:
            clearance: measured clearance (creatinine clearance or GFR)
                [ml/min], not indexed to the body surface area.
            bw: the body weight of the simulation [kg], `None` for that of the
                model.
            height: the height of the simulation [cm], `None` for that of the
                model.
            hf_severity: the severity of heart failure of the simulation.
            cp_severity: the severity of cirrhosis of the simulation.

        Returns:
            The renal function `f_renal_function`.

        Raises:
            ValueError: for a renal function outside of the validity of the
                kidney outflow (the message names the GFR and the plasma flow
                of the kidney) or of the `renal_function_range` of the model.
        """
        healthy = cls.healthy_gfr(
            bw=bw, height=height, hf_severity=hf_severity, cp_severity=cp_severity
        )
        value = clearance / float(healthy.m_as("ml/min"))
        changes = _individual_changes(bw, height, hf_severity, cp_severity, value)
        for outflow in kidney_outflow(cls.pkpd_model.model_path, changes):
            if not outflow.valid:
                raise ValueError(
                    f"The renal function {value:.3g} of a clearance of {clearance} "
                    f"ml/min is outside the validity of the kidney outflow: "
                    f"{outflow.text()}; the net filtered plasma flow must stay "
                    "below the plasma flow."
                )
        low, high = cls.renal_function_range or (0.0, math.inf)
        if not low <= value <= high:
            raise ValueError(
                f"The renal function {value:.3g} of a clearance of {clearance} "
                f"ml/min is outside of the supported range {low} to {high} "
                "(`renal_function_range`)."
            )
        return Q(value, "dimensionless")

    @classmethod
    def renal_function_limits(
        cls,
        *,
        bw: float | None = None,
        height: float | None = None,
        hf_severity: float = 0.0,
        cp_severity: float = 0.0,
    ) -> tuple[float, float]:
        """Get the range of the renal functions of an individual which the model supports.

        From 0 to the renal function at which the net filtered plasma flow of
        a filtered substance reaches the plasma flow of the kidney
        (`kidney_outflow`), the flows of the individual at renal function 1
        (the net filtered plasma flow is the GFR, `GFR_base*f_renal_function`,
        times the factors of the substance, which do not change above renal
        function 1, where the factors of chronic kidney disease are 1), the
        upper end excluded; within the `renal_function_range` of the model if
        it has one. Infinite without a filtered substance.

        Args:
            bw: the body weight of the individual [kg], `None` for that of the
                model.
            height: the height of the individual [cm], `None` for that of the
                model.
            hf_severity: the severity of heart failure of the individual.
            cp_severity: the severity of cirrhosis of the individual.

        Returns:
            The lowest and the highest renal function.

        Raises:
            TypeError: if the class has no model object (`pkpd_model`).
        """
        if not hasattr(cls, "pkpd_model"):
            raise TypeError(
                f"{cls.__name__} has no model object: set `pkpd_model` in the "
                "base of the experiments of the model."
            )
        changes = _individual_changes(bw, height, hf_severity, cp_severity, 1.0)
        outflows = kidney_outflow(cls.pkpd_model.model_path, changes)
        high = min(
            (o.plasma_flow / o.filtered for o in outflows if o.filtered > 0),
            default=math.inf,
        )
        if cls.renal_function_range is None:
            return 0.0, high
        low, own = cls.renal_function_range
        return low, min(own, high)


def _individuals(
    simulation: Simulation | Scan, inputs: frozenset[str]
) -> list[dict[str, Change]]:
    """The individuals a simulation or scan simulates, each once.

    An individual is the changes before the initialization of the ids which
    the validity of the kidney outflow depends on (`inputs`), at the start of
    the simulation: of its presimulation and itself, of every simulation of a
    dimension of simulations and every point of the dimensions of values
    without a time (`Simulation.with_values`), not the draws of the fitted
    parameters (`pkpdsim.uncertainty`). A scan over models is not checked.
    """
    scan = Scan.of(simulation)
    if scan.dimension(DimensionKind.MODELS) is not None:
        return []
    dimensions = [
        d
        for d in scan.dimensions
        if d.kind is DimensionKind.VALUES
        and d.at is None
        and d.id != UNCERTAINTY_DIMENSION
        and inputs & set(d.values)
    ]
    individuals: dict[tuple[_Frozen, ...], dict[str, Change]] = {}
    for sim in scan.simulations():
        for point in itertools.product(*(range(len(d)) for d in dimensions)):
            values = {
                target: array[k]
                for d, k in zip(dimensions, point, strict=True)
                for target, array in d.values.items()
            }
            start = sim.with_values(values) if values else sim
            changes: dict[str, Change] = {}
            if start.presimulation is not None:
                changes.update(
                    (k, v)
                    for k, v in start.presimulation.preinit_changes.items()
                    if isinstance(v, float | int | str | Quantity)
                )
            changes.update(start.preinit_changes)
            relevant = _relevant(changes, inputs)
            individuals.setdefault(_frozen(relevant), relevant)
    return list(individuals.values())


def _changes_text(changes: Mapping[str, Change]) -> str:
    """The changes of an individual, for a message."""
    if not changes:
        return "none"
    return ", ".join(f"{k} = {v}" for k, v in sorted(changes.items()))
