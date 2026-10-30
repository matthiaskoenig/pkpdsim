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
All of them are of the build of the sex of the experiment (`sex`, the female
build for a study of women): the GFR of a woman is the GFR of the female
build.

The outflow of the kidney is valid while the net filtered plasma flow of a
substance is below the plasma flow of the kidney
(`pkpdsim.pbpk.kidney_validity`): `PBPKSimulationExperiment` refuses an
experiment which simulates an individual outside of it, checked on the model
of each task, and a renal function beyond it (`renal_function`,
`renal_function_limits`).

The module imports sbmlsim; `pkpdsim.pbpk` exports its names lazily, so
`import pkpdsim.pbpk` (a factory) imports no simulator.
"""

from __future__ import annotations

import functools
import math
from collections.abc import Callable, Iterator, Mapping
from functools import lru_cache
from pathlib import Path
from types import MappingProxyType
from typing import TYPE_CHECKING, ClassVar, overload, override

import libsbml
from sbmlsim import Q
from sbmlsim.model import AbstractModel
from sbmlsim.model.model_roadrunner import RoadrunnerSBMLModel
from sbmlsim.simulation import Change, Scan, Simulation, SteadyState
from sbmlsim.simulation.definition import Value
from sbmlsim.simulation.scan import DimensionKind
from sbmlsim.simulator import Simulator
from sbmlsim.uncertainty import Uncertain
from sbmlsim.units import Quantity

from pkpdsim import impairment
from pkpdsim.experiment import PKPDSimulationExperiment
from pkpdsim.pbpk import ids
from pkpdsim.pbpk.kidney_validity import (
    changes_text,
    covariate_changes,
    individual_values,
    kidney_outflow,
    simulated_individuals,
)
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import number, unit_text
from pkpdsim.scan import start_values
from pkpdsim.sex import Sex

if TYPE_CHECKING:
    from pkpdsim.pbpk.definition import PBPKDefinition

#: the absolute and the relative tolerance of the rates of change (per minute)
#: and the longest time (minutes) of the equilibration of a condition, see
#: `PBPKSimulationExperiment.equilibrated`
EQUILIBRATION_TOLERANCES: tuple[float, float, float] = (1e-12, 1e-10, 1e14)

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
    is no parameter of the model: the library builds the model of one sex, and
    the GFR of a woman is the GFR of the female build.

    Args:
        model_path: the flat SBML model of the library, the build of the sex
            of the individual.
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
    changes = covariate_changes(bw, height, hf_severity, cp_severity)
    gfr = individual_values(Path(model_path), changes).ml_per_min(ids.GFR_BASE)
    return Q(gfr, "ml/min")


def accumulating_states(sbml: str) -> tuple[str, ...]:
    """Get the states of a model which nothing drains: they only accumulate.

    A state (the variable of a rate rule, a species which reactions change)
    which no rate reads, neither a rate rule nor a kinetic law, directly or
    through the assignment rules, and which no reaction consumes, e.g. the
    urine and the feces of the library or the urine of a PD component: with an
    inflow at the steady state of the rest of the model it grows without a
    steady state. A species which a rate reads is read with its compartment.

    Args:
        sbml: the SBML of the model, flat.

    Returns:
        The ids, sorted.
    """
    model = libsbml.readSBMLFromString(sbml).getModel()

    def names(ast: libsbml.ASTNode | None) -> set[str]:
        found: set[str] = set()
        stack = [] if ast is None else [ast]
        while stack:
            node = stack.pop()
            if node.getType() == libsbml.AST_NAME:
                found.add(node.getName())
            stack.extend(node.getChild(k) for k in range(node.getNumChildren()))
        return found

    assigned: dict[str, set[str]] = {}
    rated: dict[str, set[str]] = {}
    for rule in model.getListOfRules():
        if rule.isAssignment():
            assigned[rule.getVariable()] = names(rule.getMath())
        elif rule.isRate():
            rated[rule.getVariable()] = names(rule.getMath())
    read = set().union(*rated.values())
    consumed: set[str] = set()
    changed: set[str] = set()
    for reaction in model.getListOfReactions():
        law = reaction.getKineticLaw()
        if law is not None:
            local = {p.getId() for p in law.getListOfLocalParameters()}
            read |= names(law.getMath()) - local
        consumed |= {r.getSpecies() for r in reaction.getListOfReactants()}
        changed |= {r.getSpecies() for r in reaction.getListOfProducts()}
    # what the rates read through the assignment rules, and the compartments of
    # the species they read
    compartments = {s.getId(): s.getCompartment() for s in model.getListOfSpecies()}
    stack = list(read)
    while stack:
        sid = stack.pop()
        for further in (assigned.get(sid, set()) | {compartments.get(sid, "")}) - read:
            if further:
                read.add(further)
                stack.append(further)
    states = set(rated) | (changed - consumed)
    return tuple(sorted(states - read - consumed))


@lru_cache(maxsize=4)
def _accumulation_model(
    model_path: Path, mtime_ns: int
) -> tuple[RoadrunnerSBMLModel, tuple[str, ...]]:
    """Load a model once per process with its accumulating states as selections."""
    states = accumulating_states(model_path.read_text(encoding="utf-8"))
    model = RoadrunnerSBMLModel(
        source=model_path, selections=["time", *states], jit_optimization="O0"
    )
    return model, states


#: the growing states of a model file by its path, modification time and the
#: tolerances of the equilibration, see `growing_states`
_GROWING: dict[tuple[Path, int, float, float, float], tuple[str, ...]] = {}


def growing_states(
    model_path: Path,
    condition: Mapping[str, Value],
    steady_state: SteadyState,
) -> tuple[str, ...]:
    """Get the states which grow without a steady state during an equilibration.

    The accumulating states of the model (`accumulating_states`) whose value
    after the equilibration of the condition differs from their value at the
    initialization by more than the absolute tolerance (and 1e-9 of the
    value): an inflow at the steady state, which the presimulation reaches
    only relative to their growing amount (after about `1/relative_tolerance`
    minutes), so their value at the start of the simulation is meaningless
    unless they are set back. They are found once per model file (and its
    modification time) and tolerances in a process, in the condition of the
    first equilibration: whether a state has no outflow is a property of the
    model.

    Args:
        model_path: the flat SBML model.
        condition: the changes before the initialization.
        steady_state: the presimulation.

    Returns:
        The ids of the growing states.

    Raises:
        SteadyStateError: if the model reaches no steady state in the
            condition by `max_time` (sbmlsim).
    """
    path = Path(model_path)
    key = (
        path,
        path.stat().st_mtime_ns,
        steady_state.absolute_tolerance,
        steady_state.relative_tolerance,
        steady_state.max_time,
    )
    if key not in _GROWING:
        model, states = _accumulation_model(path, key[1])
        if not states:
            _GROWING[key] = ()
            return ()
        simulator = Simulator(n_workers=1)
        initial = simulator.run(
            model, Simulation(end=1, preinit_changes=dict(condition))
        )
        equilibrated = simulator.run(
            model,
            Simulation(
                end=1, preinit_changes=dict(condition), presimulation=steady_state
            ),
        )
        growing = []
        for sid in states:
            before = float(start_values(initial, sid).item())
            after = float(start_values(equilibrated, sid).item())
            if abs(after - before) > steady_state.absolute_tolerance + 1e-9 * abs(
                before
            ):
                growing.append(sid)
        _GROWING[key] = tuple(growing)
    return _GROWING[key]


def molecular_weights(definition: PBPKDefinition) -> Mapping[str, Quantity]:
    """Get the molecular weights of a definition, for the conversion of the data.

    The molecular weight `Mr_<s>` of every substance of the definition, and
    of a substance of the data which is no substance of the model (the
    substance of a measurement which the model reads in another form, e.g. a
    protein of a PD readout) as a parameter of the definition `Mr_<x>` in
    g/mole; the same numbers the model has. The base of the experiments of a
    model reads them instead of a copy:

    ```python
    # not executed
    class TorasemideSimulationExperiment(PBPKSimulationExperiment):
        pkpd_model = MODEL
        Mr: ClassVar[Mapping[str, Quantity]] = molecular_weights(DEFINITION)

        def to_molar(self, dset: DataSet) -> DataSet:
            dset.unit_conversion("mean", 1 / self.Mr["tor"])
    ```

    Args:
        definition: the definition of the model.

    Returns:
        The molecular weight by the id of the substance (`tor`, without
        `Mr_`), a read-only mapping of quantities in g/mole.

    Raises:
        ValueError: for a molecular weight which is a formula, or a parameter
            `Mr_<x>` of another unit.
    """
    weights: dict[str, Quantity] = {}
    entries = [
        (s.sid, ids.mr(s.sid), s.mr, U.g_per_mole) for s in definition.substances
    ]
    entries += [
        (p.sid.removeprefix(ids.mr("")), p.sid, p.value, p.unit)
        for p in definition.parameters
        if p.sid.startswith(ids.mr(""))
    ]
    for key, sid, q, unit in entries:
        if unit_text(unit) != unit_text(U.g_per_mole):
            raise ValueError(
                f"The molecular weight '{sid}' is in {unit_text(unit)}, not in g/mole."
            )
        value = number(q)
        if value is None:
            raise ValueError(
                f"The molecular weight '{sid}' is a formula, whose number only the "
                "model has."
            )
        weights[key] = Q(float(value), "g/mole")
    return MappingProxyType(weights)


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
    object or another sex (a study of women gets the GFR of the female build),
    gets the GFR of its own model, and a model overrides the preset in its
    base:

    ```python
    # not executed
    class HCTZSimulationExperiment(PBPKSimulationExperiment):
        pkpd_model = MODEL
        renal_map: ClassVar[Mapping[str, float]] = ScaledRenalMap(
            impairment.RENAL_FUNCTION
        )
    ```

    The map of a class for another sex than the sex of the class, e.g. the
    group of women of a study of men, is `of_sex` (see
    `PBPKSimulationExperiment.renal_map_of`).
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
        #: the classmethod `healthy_gfr` of the class a bound map is bound to
        #: and the sex of the class, `None` for an unbound map or a map with
        #: its own GFR
        self._class_gfr: Callable[..., Quantity] | None = None
        self._class_sex: Sex = Sex.MALE
        #: the maps of the other sexes of a bound map, see `of_sex`
        self._sexes: dict[Sex, ScaledRenalMap] = {}

    def __getstate__(self) -> dict[str, object]:
        """Get the state to pickle: without the maps bound to classes or of a sex."""
        return {**self.__dict__, "_bound": {}, "_sexes": {}}

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
            bound._class_gfr = gfr
            bound._class_sex = Sex(getattr(owner, "sex", Sex.MALE))
            self._bound[owner] = bound
        return bound

    def of_sex(self, sex: Sex) -> ScaledRenalMap:
        """Get the map of a class at the healthy GFR of the build of a sex.

        Args:
            sex: the sex.

        Returns:
            The map of the class this map is bound to, on its `healthy_gfr` of
            the sex: this map for the sex of the class, else a map created on
            the first call.

        Raises:
            TypeError: for a map with its own GFR (of one build) or a map which
                no class bound.
        """
        class_gfr = self._class_gfr
        if class_gfr is None:
            if self._gfr is not None:
                raise TypeError(
                    "A map with its own GFR is of one build, it has no map of "
                    "another sex."
                )
            raise TypeError(
                "The map has no class: get the map of a sex from the map of a "
                "class (`renal_map_of`)."
            )
        sex = Sex(sex)
        if sex == self._class_sex:
            return self
        other = self._sexes.get(sex)
        if other is None:
            gfr = functools.partial(class_gfr, sex=sex)
            other = ScaledRenalMap(self._preset, gfr, self._reference)
            self._sexes[sex] = other
        return other

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

    The healthy GFR, the renal function of a clearance, the range of the renal
    functions and the renal classes are those of the build of the sex of the
    experiment (`sex`, see `PKPDSimulationExperiment`); the methods take another
    sex with `sex`, e.g. for the group of women of a study of men, whose renal
    classes are `renal_map_of(Sex.FEMALE)`.

    The initialization refuses an experiment which simulates an individual
    outside the validity of the kidney outflow (`kidney_outflow`): the net
    filtered plasma flow of a substance must stay below the plasma flow of
    the kidney, at the start of every simulation and every point of a scan,
    each on the model its task runs (`simulated_individuals`; not checked: a
    change during a simulation, also of a `dynamic` id at its start, a
    dimension with a time, a scan over models and the draws of the fitted
    parameters of a run over them).

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

    @overload
    def equilibrated(
        self,
        simulation: Simulation,
        *,
        resets: Mapping[str, Quantity | float] | None = None,
        model: str = "model",
        absolute_tolerance: float = EQUILIBRATION_TOLERANCES[0],
        relative_tolerance: float = EQUILIBRATION_TOLERANCES[1],
        max_time: float = EQUILIBRATION_TOLERANCES[2],
    ) -> Simulation: ...

    @overload
    def equilibrated(
        self,
        simulation: Scan,
        *,
        resets: Mapping[str, Quantity | float] | None = None,
        model: str = "model",
        absolute_tolerance: float = EQUILIBRATION_TOLERANCES[0],
        relative_tolerance: float = EQUILIBRATION_TOLERANCES[1],
        max_time: float = EQUILIBRATION_TOLERANCES[2],
    ) -> Scan: ...

    def equilibrated(
        self,
        simulation: Simulation | Scan,
        *,
        resets: Mapping[str, Quantity | float] | None = None,
        model: str = "model",
        absolute_tolerance: float = EQUILIBRATION_TOLERANCES[0],
        relative_tolerance: float = EQUILIBRATION_TOLERANCES[1],
        max_time: float = EQUILIBRATION_TOLERANCES[2],
    ) -> Simulation | Scan:
        """Get a simulation which starts from the drug-free steady state of its condition.

        The state of a model at its start is the steady state of the healthy
        reference individual; a condition of the simulation (an impairment, a
        covariate, a PD component with a homeostasis) moves a slow state away
        from it, which a study measures as the basal state of its subjects
        before the first dose. The simulation equilibrates first: the
        presimulation of sbmlsim (`SteadyState`) integrates the model with the
        changes before the initialization of the simulation, the condition,
        without the doses until the rates of change are within the
        tolerances, and the simulation starts from that state.

        - The doses of the changes before the initialization, the dose
          states and infusion rates of the library (`PODOSE_<s>`,
          `IVDOSE_<s>`, `SOLDOSE_<s>`, `Ri_<s>`), are given by a change at
          the start instead, after the equilibration; the doses of the changes
          at later times stay. A dose of zero moves as well, so a dimension of
          a scan over it replaces it at the start.
        - The sinks, whose amounts grow without a steady state (the urine and
          the feces of the library, `Aurine_<s>`, `Afeces_<s>`, and what a PD
          component accumulates, `resets`), are set back at the start: they
          are in steady state when their rate is within the tolerances
          relative to their amount, which they reach as well. A state without
          an outflow which grows during the equilibration and is neither a
          sink of the library nor in `resets` is refused, named
          (`growing_states`, checked once per model file in a process).
        - The base `<id>_0` of a delta observer (`PBPKDefinition.deltas`)
          is set at the start to the value after the equilibration, so the
          delta `<id>_delta` is 0 at the start of the simulation.
        - A model which reaches no steady state by `max_time` fails with the
          `SteadyStateError` of sbmlsim (within the `ScanError` of the run).

        The presimulation is the choice over a hidden pre-run of a duration
        from a negative start (tested on a homeostasis of 7 days: the pre-run
        of 28 days leaves it 1.8 % off the steady state, the presimulation is
        exact), and its default tolerances are tighter than those of sbmlsim:
        the rates are per minute, the time unit of the model, so the relative
        tolerance of sbmlsim, 1e-6, stops a homeostasis of 7 days 5e-5 off;
        at 1e-10 a sink reaches the tolerance after 1e10 minutes, within
        `max_time`.

        Every point of a scan equilibrates in its own condition; a dimension
        over a dose which the simulation does not give is refused (it would
        be a dose before the equilibration). The changes of a fit (the fitted
        parameters) and of a run over the draws of the fitted parameters are
        changes before the initialization as well, so they apply during the
        equilibration. The simulation runs on the build of its sex as every
        simulation of the experiment; the builds have the same ids.

        Args:
            simulation: the simulation or the scan, in the condition.
            resets: the values at the start of the further states which grow
                without a steady state, e.g. the urine of a PD component.
            model: the key of the loaded model whose ids are read, see
                `add_doses`.
            absolute_tolerance: the absolute tolerance of the rates of the
                equilibration (per minute).
            relative_tolerance: the relative tolerance of the rates.
            max_time: the longest equilibration (minutes).

        Returns:
            The equilibrated simulation or scan, a new object.

        Raises:
            ValueError: if the model is not loaded, the simulation has a
                presimulation, a dose is set before the initialization and by a
                change at the start, a state grows during the equilibration
                which is not set back, or a dimension of a scan sets a dose which
                the simulation does not give or varies simulations.
            ScanError: if the model reaches no steady state in the condition
                by `max_time` (a `SteadyStateError` of sbmlsim).
        """
        loaded = self._models.get(model)
        symbols = getattr(loaded, "symbols", None)
        if symbols is None:
            raise ValueError(
                f"The model '{model}' of the experiment '{self.sid}' is not loaded: "
                "the equilibration is created in `simulations()` or later."
            )
        dose_prefixes = tuple(
            f("") for f in (ids.po_dose, ids.iv_dose, ids.sol_dose, ids.infusion_rate)
        )
        sink_prefixes = (ids.urine(""), ids.feces(""))
        sinks = sorted(s for s in symbols.species if s.startswith(sink_prefixes))
        # the bases of the deltas, evaluated after the equilibration
        deltas = sorted(
            sid
            for sid in symbols.parameters
            if ids.delta_base(sid) in symbols.parameters
            and ids.delta(sid) in symbols.parameters
        )
        at_start: dict[str, float | Quantity | str] = {
            **dict.fromkeys(sinks, 0.0),
            **(resets or {}),
            **{ids.delta_base(sid): sid for sid in deltas},
        }
        steady_state = SteadyState(
            absolute_tolerance=absolute_tolerance,
            relative_tolerance=relative_tolerance,
            max_time=max_time,
        )
        base = simulation.simulation if isinstance(simulation, Scan) else simulation
        equilibrated = _equilibrated(base, dose_prefixes, at_start, steady_state)
        path = _model_path(self.models().get(model))
        if path is not None:
            growing = sorted(
                set(growing_states(path, equilibrated.preinit_changes, steady_state))
                - set(at_start)
            )
            if growing:
                raise ValueError(
                    f"The states {growing} grow without a steady state during the "
                    f"equilibration in the experiment '{self.sid}' (they have no "
                    "outflow): give their values at the start in `resets`."
                )
        if not isinstance(simulation, Scan):
            return equilibrated
        given = {target for change in equilibrated.changes for target in change.values}
        for dimension in simulation.dimensions:
            if dimension.kind != DimensionKind.VALUES:
                raise ValueError(
                    f"The dimension '{dimension.id}' varies {dimension.kind}: "
                    "equilibrate every simulation of it instead."
                )
            missing = sorted(
                target
                for target in dimension.values
                if target.startswith(dose_prefixes)
                and dimension.at is None
                and target not in given
            )
            if missing:
                raise ValueError(
                    f"The dimension '{dimension.id}' sets the doses {missing}, "
                    "which the simulation does not give: they would be doses "
                    "before the equilibration; give them in the simulation."
                )
        return Scan(simulation=equilibrated, dimensions=simulation.dimensions)

    def _check_kidney_outflow(self) -> None:
        """Check the kidney outflow of every individual the experiment simulates.

        Every task is checked on the model it runs (`models`), a model of the
        library, whose filtered substances `kidney_outflow` reads; a model
        given as SBML text is not checked.

        Raises:
            ValueError: for an individual outside the validity.
        """
        models = self.models()
        for task in self._tasks.values():
            path = _model_path(models.get(task.model_id))
            if path is None:
                continue
            key = task.simulation_id
            for changes in simulated_individuals(self._simulations[key], path):
                for outflow in kidney_outflow(path, changes):
                    if not outflow.valid:
                        raise ValueError(
                            f"The experiment '{self.sid}' simulates in '{key}' an "
                            "individual outside the validity of the kidney "
                            f"outflow: {outflow.text()}; the net filtered plasma "
                            "flow must stay below the plasma flow. The changes of "
                            f"the individual: {changes_text(changes)}."
                        )

    @classmethod
    def _build(cls, sex: Sex | None) -> Path:
        """Get the build of a sex of the model object, of `sex` of the class by default.

        Raises:
            TypeError: if the class has no model object (`pkpd_model`).
            ValueError: if the model object has no build of the sex.
        """
        if not hasattr(cls, "pkpd_model"):
            raise TypeError(
                f"{cls.__name__} has no model object: set `pkpd_model` in the "
                "base of the experiments of the model."
            )
        return cls.pkpd_model.model_path_of(cls.sex if sex is None else sex)

    @classmethod
    def healthy_gfr(
        cls,
        *,
        bw: float | None = None,
        height: float | None = None,
        hf_severity: float = 0.0,
        cp_severity: float = 0.0,
        sex: Sex | None = None,
    ) -> Quantity:
        """Get the healthy glomerular filtration rate of an individual (ml/min).

        `healthy_gfr` of the module for the build of the sex of the model
        object, cached for all experiments.

        Args:
            bw: the body weight of the individual [kg], `None` for that of the
                model.
            height: the height of the individual [cm], `None` for that of the
                model.
            hf_severity: the severity of heart failure `HF_severity` (0 to 4).
            cp_severity: the severity of cirrhosis `CP_severity` (0 to 3).
            sex: the sex of the individual, `None` for the sex of the
                experiment (`sex`).

        Returns:
            The GFR at renal function 1 in ml/min, a quantity.

        Raises:
            TypeError: if the class has no model object (`pkpd_model`).
            ValueError: if the model object has no build of the sex.
        """
        return healthy_gfr(
            cls._build(sex),
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
        sex: Sex | None = None,
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
            sex: the sex of the subjects, `None` for the sex of the experiment
                (`sex`); the GFR of a woman is the GFR of the female build.

        Returns:
            The renal function `f_renal_function`.

        Raises:
            ValueError: for a renal function outside of the validity of the
                kidney outflow (the message names the GFR and the plasma flow
                of the kidney) or of the `renal_function_range` of the model,
                or if the model object has no build of the sex.
        """
        healthy = cls.healthy_gfr(
            bw=bw,
            height=height,
            hf_severity=hf_severity,
            cp_severity=cp_severity,
            sex=sex,
        )
        value = clearance / float(healthy.m_as("ml/min"))
        changes = covariate_changes(bw, height, hf_severity, cp_severity, value)
        for outflow in kidney_outflow(cls._build(sex), changes):
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
        sex: Sex | None = None,
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
            sex: the sex of the individual, `None` for the sex of the
                experiment (`sex`).

        Returns:
            The lowest and the highest renal function.

        Raises:
            TypeError: if the class has no model object (`pkpd_model`).
            ValueError: if the model object has no build of the sex.
        """
        changes = covariate_changes(bw, height, hf_severity, cp_severity, 1.0)
        outflows = kidney_outflow(cls._build(sex), changes)
        high = min(
            (o.plasma_flow / o.filtered for o in outflows if o.filtered > 0),
            default=math.inf,
        )
        if cls.renal_function_range is None:
            return 0.0, high
        low, own = cls.renal_function_range
        return low, min(own, high)

    @classmethod
    def renal_map_of(cls, sex: Sex) -> Mapping[str, float]:
        """Get the renal classes of the subjects of a sex.

        `renal_map` is of the sex of the experiment; the classes of the other
        sex, e.g. of the group of women of a study of men, are those of the
        same preset at the healthy GFR of the build of that sex
        (`ScaledRenalMap.of_sex`). A map which is no `ScaledRenalMap` is not
        relative to a healthy GFR and is the map of every sex.

        Args:
            sex: the sex.

        Returns:
            The renal function of every class.

        Raises:
            TypeError: for a `ScaledRenalMap` with its own GFR (of one build)
                and another sex than the sex of the experiment.
        """
        renal_map = cls.renal_map
        if Sex(sex) == Sex(cls.sex) or not isinstance(renal_map, ScaledRenalMap):
            return renal_map
        return renal_map.of_sex(sex)


def _equilibrated(
    simulation: Simulation,
    dose_prefixes: tuple[str, ...],
    at_start: Mapping[str, float | Quantity | str],
    steady_state: SteadyState,
) -> Simulation:
    """Equilibrate a simulation, see `PBPKSimulationExperiment.equilibrated`.

    Args:
        simulation: the simulation.
        dose_prefixes: the prefixes of the ids of the doses.
        at_start: the further values at the start, after the doses.
        steady_state: the presimulation.

    Raises:
        ValueError: for a simulation with a presimulation or a dose before the
            initialization and in a change at the start.
    """
    if simulation.presimulation is not None:
        raise ValueError(
            "The simulation has a presimulation already; it is equilibrated once."
        )
    preinit = {**simulation.preinit_changes, **simulation.uncertain_preinit}
    doses = {k: v for k, v in preinit.items() if k.startswith(dose_prefixes)}
    condition = {k: v for k, v in preinit.items() if k not in doses}
    start = simulation.start
    at_the_start = {
        target
        for change in simulation.changes
        if any(
            math.isclose(
                _time_value(simulation, t),
                _time_value(simulation, start),
                abs_tol=1e-12,
            )
            for t in change.times
        )
        for target in change.values
    }
    twice = sorted(set(doses) & at_the_start)
    if twice:
        raise ValueError(
            f"The doses {twice} are set before the initialization and by a change at "
            "its start; give each dose once, the equilibration moves the doses "
            "before the initialization to the start."
        )
    values: dict[str, float | Quantity | str | Uncertain] = {**doses, **at_start}
    return Simulation(
        time_unit=simulation.time_unit,
        start=start,
        end=simulation.end,
        preinit_changes=condition,
        changes=[Change(start, values), *simulation.changes],
        presimulation=steady_state,
        times=None if simulation.times is None else list(simulation.times),
        steps=simulation.steps,
        time_shift=simulation.time_shift,
    )


def _time_value(simulation: Simulation, time: object) -> float:
    """Get a time of a simulation as a number in its time unit."""
    if isinstance(time, Quantity):
        unit = simulation.time_unit
        return float(time.m_as(unit) if unit else time.magnitude)
    return float(time)  # ty: ignore[invalid-argument-type]


def _model_path(model: AbstractModel | Path | str | None) -> Path | None:
    """The file of a model of an experiment, `None` for SBML text or no model."""
    if isinstance(model, AbstractModel):
        return model.source.path
    if model is None or (isinstance(model, str) and model.lstrip().startswith("<")):
        return None
    return Path(model)
