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
All of them are of the build of the sex of the subjects (`sex` of the
experiment, the female build for a study of women, or the sex of a
simulation, `simulation=`, in a study of a group of each sex): the GFR of a
woman is the GFR of the female build.

The outflow of the kidney is valid while the filtered plasma water (the
GFR) and the net filtered plasma flow of every substance are below the plasma
flow of the kidney (`pkpdsim.pbpk.kidney_validity`):
`PBPKSimulationExperiment` refuses an
experiment which simulates an individual outside of it, checked on the model
of each task, and a renal function beyond it (`renal_function`,
`renal_function_limits`).

The module imports sbmlsim; `pkpdsim.pbpk` exports its names lazily, so
`import pkpdsim.pbpk` (a factory) imports no simulator.
"""

from __future__ import annotations

import functools
import math
from collections.abc import Callable, Iterator, Mapping, Sequence
from pathlib import Path
from types import MappingProxyType
from typing import TYPE_CHECKING, ClassVar, overload, override

import numpy as np
from sbmlsim import Q
from sbmlsim.model import AbstractModel
from sbmlsim.simulation import Dimension, Scan, Simulation
from sbmlsim.units import Quantity

from pkpdsim import impairment
from pkpdsim.experiment import PKPDSimulationExperiment
from pkpdsim.pbpk import ids
from pkpdsim.pbpk.equilibration import (
    EQUILIBRATION_TOLERANCES,
    equilibrated_simulation,
)
from pkpdsim.pbpk.kidney_validity import (
    changes_text,
    covariate_changes,
    individual_values,
    kidney_outflow,
    simulated_individuals,
)
from pkpdsim.pbpk.physiology import AGE_REFERENCE, normal_gfr, reference_individual
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import number, unit_text
from pkpdsim.pbpk.virtual import VirtualPopulation
from pkpdsim.sex import Sex

if TYPE_CHECKING:
    from sbmlsim.plot import Plot

    from pkpdsim.pbpk.definition import PBPKDefinition

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
    age: float | None = None,
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
        age: the age of the individual [yr], `None` for that of the model;
            the GFR matures in the children and declines in the adults.

    Returns:
        The GFR at renal function 1 in ml/min, a quantity as the renal
        function of `PBPKSimulationExperiment.renal_function`.
    """
    changes = covariate_changes(bw, height, hf_severity, cp_severity, age=age)
    gfr = individual_values(Path(model_path), changes).ml_per_min(ids.GFR_BASE)
    return Q(gfr, "ml/min")


def _age(age: float | None) -> float:
    """The age of an individual, the adult of ICRP 89 (35 years) if not given."""
    return AGE_REFERENCE if age is None else float(age)


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

    The map of a class is of the sex of the class (`sex`); the map for
    another sex, e.g. the group of women of a study of men, is `of_sex` (see
    `PBPKSimulationExperiment.renal_map_of`, which takes the key of the
    simulation).
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
            sex = Sex(getattr(owner, "sex", Sex.MALE))
            bound = ScaledRenalMap(
                self._preset, functools.partial(gfr, sex=sex), self._reference
            )
            bound._class_gfr = gfr
            bound._class_sex = sex
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
    subjects: the sex of the experiment (`sex`, see
    `PKPDSimulationExperiment`), or of a simulation by its key, `simulation=`
    (`sex_of`, the sex of `simulation_sexes`), e.g.
    `renal_function(60.0, simulation="po_women")` and
    `renal_map_of(simulation="po_women")` for the group of women of a study of
    men; `sex=` gives a sex directly. An experiment whose `simulation_sexes`
    names another sex than `sex` refuses a call with neither, since the sex of
    the subjects is not the sex of the experiment for every simulation (the
    attribute `renal_map` stays the map of `sex`).

    The initialization refuses an experiment which simulates an individual
    outside the validity of the kidney outflow (`kidney_outflow`): the
    filtered plasma water (the GFR) and the net filtered plasma flow of a
    substance must stay below the plasma flow of the kidney, at the start of every simulation and every point of a scan,
    each on the model its task runs (`simulated_individuals`; not checked: a
    change during a simulation, also of a `dynamic` id at its start, a
    dimension with a time, a scan over models and the draws of the fitted
    parameters of a run over them).

    Attributes:
        renal_function_range: a narrower range of the renal functions the
            model supports, `None` (the default) for the range of the validity
            of the kidney outflow, from 0 (anuria, the substances leave by the
            other routes alone) to the renal function at which the GFR or the
            net filtered plasma flow of a substance reaches the plasma flow of
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

    @overload
    def equilibrated(
        self,
        simulation: Simulation | Scan,
        *,
        resets: Mapping[str, Quantity | float] | None = None,
        model: str = "model",
        absolute_tolerance: float = EQUILIBRATION_TOLERANCES[0],
        relative_tolerance: float = EQUILIBRATION_TOLERANCES[1],
        max_time: float = EQUILIBRATION_TOLERANCES[2],
    ) -> Simulation | Scan: ...

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
          the feces of the library, `Aurine_<s>`, `Afeces_<s>`, the readouts
          of the intestine which nothing drains, the amount absorbed
          `Aabs_<s>` and the oral dose applied `Aapp_<s>` (`GU__Aabs_<s>`,
          `GU__Aapp_<s>` with the readout `Fa_<s>` alone), and what a PD
          component accumulates, `resets`), are set back at the start: they
          are in steady state when their rate is within the tolerances
          relative to their amount, which they reach as well. A state without
          an outflow which grows during the equilibration and is neither a
          sink of the library nor in `resets` is refused, named
          (`growing_states`, checked once per model file in a process). A
          reset names a species by its id; a species whose id is its
          concentration is set as `[sid]`, as the zeros are.
        - The species which are 0 without a dose, of the substances without a
          baseline (the drugs and their metabolites) which start at 0 and
          which no reaction feeds from a species with a start, are set to
          exactly 0 at the start (`drug_free_zeros`): the presimulation leaves
          a residue of the rounding in them (1.5e-29 mM of hctz), which a
          pre-dose point of zero would count in the fold metrics. A species
          which a change of the condition or at the start sets keeps its
          value; an endogenous substance (a baseline, e.g. creatinine) keeps
          its equilibrated values.
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
        simulation of the experiment; the builds have the same ids. It is
        `equilibrated_simulation(simulation, model_path)` of
        `pkpdsim.pbpk.equilibration` on the file of the loaded model, the
        equilibration outside an experiment as well (the groups of a
        sensitivity analysis).

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
            ValueError: if the model is not loaded or not a file, the
                simulation has a presimulation, a dose is set before the
                initialization and by a change at the start, a state grows
                during the equilibration which is not set back, or a dimension
                of a scan sets a dose which the simulation does not give or
                varies simulations; the message is noted with the experiment.
            ScanError: if the model reaches no steady state in the condition
                by `max_time` (a `SteadyStateError` of sbmlsim).
        """
        if getattr(self._models.get(model), "symbols", None) is None:
            raise ValueError(
                f"The model '{model}' of the experiment '{self.sid}' is not loaded: "
                "the equilibration is created in `simulations()` or later."
            )
        path = _model_path(self.models().get(model))
        if path is None:
            raise ValueError(
                f"The model '{model}' of the experiment '{self.sid}' is no file: "
                "the equilibration reads the ids of the model from its file."
            )
        try:
            return equilibrated_simulation(
                simulation,
                path,
                resets=resets,
                absolute_tolerance=absolute_tolerance,
                relative_tolerance=relative_tolerance,
                max_time=max_time,
            )
        except ValueError as err:
            err.add_note(f"The equilibration in the experiment '{self.sid}'.")
            raise

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
                            f"outflow: {outflow.text()}; the filtered flow must "
                            "stay below the plasma flow. The changes of the "
                            f"individual: {changes_text(changes)}."
                        )

    @classmethod
    def subjects_sex(cls, sex: Sex | None = None, simulation: str | None = None) -> Sex:
        """Get the sex of the subjects of a call of the renal functions.

        Args:
            sex: the sex, `None` for the sex of `simulation`.
            simulation: the key of the simulation whose subjects are meant
                (`sex_of`), `None` for the sex of the experiment (`sex`).

        Returns:
            `sex`, else the sex of the simulation, else the sex of the experiment.

        Raises:
            ValueError: for both `sex` and `simulation`, and for neither in an
                experiment whose `simulation_sexes` names another sex than
                `sex` (its subjects have no one sex).
        """
        if sex is not None and simulation is not None:
            raise ValueError(
                f"{cls.__name__}: give the sex of the subjects (`sex`) or the key "
                "of their simulation (`simulation`), not both."
            )
        if simulation is not None:
            return cls.sex_of(simulation)
        if sex is not None:
            return Sex(sex)
        own = Sex(cls.sex)
        others = sorted(
            {str(Sex(s)) for s in cls.simulation_sexes.values()} - {str(own)}
        )
        if others:
            raise ValueError(
                f"{cls.__name__} simulates subjects of another sex than its sex "
                f"'{own}' ({', '.join(others)}, `simulation_sexes`): give the key of "
                "the simulation (`simulation=`) or the sex (`sex=`)."
            )
        return own

    @classmethod
    def _build(cls, sex: Sex | None, simulation: str | None = None) -> Path:
        """Get the build of the sex of the subjects of a call (`subjects_sex`).

        Raises:
            TypeError: if the class has no model object (`pkpd_model`).
            ValueError: if the model object has no build of the sex, or the sex
                is not given in an experiment of both sexes.
        """
        if not hasattr(cls, "pkpd_model"):
            raise TypeError(
                f"{cls.__name__} has no model object: set `pkpd_model` in the "
                "base of the experiments of the model."
            )
        return cls.pkpd_model.model_path_of(cls.subjects_sex(sex, simulation))

    @classmethod
    def healthy_gfr(
        cls,
        *,
        bw: float | None = None,
        height: float | None = None,
        hf_severity: float = 0.0,
        cp_severity: float = 0.0,
        age: float | None = None,
        sex: Sex | None = None,
        simulation: str | None = None,
    ) -> Quantity:
        """Get the healthy glomerular filtration rate of an individual (ml/min).

        `healthy_gfr` of the module for the build of the sex of the subjects
        (`subjects_sex`), cached for all experiments.

        Args:
            bw: the body weight of the individual [kg], `None` for that of the
                model.
            height: the height of the individual [cm], `None` for that of the
                model.
            hf_severity: the severity of heart failure `HF_severity` (0 to 4).
            cp_severity: the severity of cirrhosis `CP_severity` (0 to 3).
            age: the age of the individual [yr], `None` for that of the model.
            sex: the sex of the individual, `None` for the sex of `simulation`.
            simulation: the key of the simulation of the individual, `None`
                for the sex of the experiment (`sex`).

        Returns:
            The GFR at renal function 1 in ml/min, a quantity.

        Raises:
            TypeError: if the class has no model object (`pkpd_model`).
            ValueError: if the model object has no build of the sex, or the sex
                is not given in an experiment of both sexes (`subjects_sex`).
        """
        return healthy_gfr(
            cls._build(sex, simulation),
            bw=bw,
            height=height,
            hf_severity=hf_severity,
            cp_severity=cp_severity,
            age=age,
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
        age: float | None = None,
        sex: Sex | None = None,
        simulation: str | None = None,
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
            age: the age of the simulation [yr], `None` for that of the model.
            sex: the sex of the subjects, `None` for the sex of `simulation`;
                the GFR of a woman is the GFR of the female build.
            simulation: the key of the simulation of the subjects, whose sex
                (`sex_of`) gives the build, `None` for the sex of the
                experiment (`sex`).

        Returns:
            The renal function `f_renal_function`.

        Raises:
            ValueError: for a renal function outside of the validity of the
                kidney outflow (the message names the GFR and the plasma flow
                of the kidney) or of the `renal_function_range` of the model,
                if the model object has no build of the sex, or if the sex is
                not given in an experiment of both sexes (`subjects_sex`).
        """
        sex = cls.subjects_sex(sex, simulation)
        healthy = cls.healthy_gfr(
            bw=bw,
            height=height,
            hf_severity=hf_severity,
            cp_severity=cp_severity,
            age=age,
            sex=sex,
        )
        value = clearance / float(healthy.m_as("ml/min"))
        changes = covariate_changes(
            bw, height, hf_severity, cp_severity, value, age=age
        )
        for outflow in kidney_outflow(cls._build(sex), changes):
            if not outflow.valid:
                raise ValueError(
                    f"The renal function {value:.3g} of a clearance of {clearance} "
                    f"ml/min is outside the validity of the kidney outflow: "
                    f"{outflow.text()}; the filtered flow must stay below the "
                    "plasma flow."
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
        age: float | None = None,
        sex: Sex | None = None,
        simulation: str | None = None,
    ) -> tuple[float, float]:
        """Get the range of the renal functions of an individual which the model supports.

        From 0 to the smaller of the renal functions at which the filtered
        plasma water, the GFR, reaches the plasma flow of the kidney
        `Qki*(1 - HCT)` and at which the net filtered plasma flow of a
        filtered substance reaches the flow which carries it
        (`kidney_outflow`), the flows of the individual at renal function 1
        (the GFR is `GFR_base*f_renal_function`, the net filtered plasma flow
        the GFR times the factors of the substance, which do not change above
        renal function 1, where the factors of chronic kidney disease are 1),
        the upper end excluded; within the `renal_function_range` of the
        model if it has one. The water is the tighter bound (`fu <= 1`,
        `BP >= 1 - HCT`); infinite for a model without a filtered substance.

        Args:
            bw: the body weight of the individual [kg], `None` for that of the
                model.
            height: the height of the individual [cm], `None` for that of the
                model.
            hf_severity: the severity of heart failure of the individual.
            cp_severity: the severity of cirrhosis of the individual.
            age: the age of the individual [yr], `None` for that of the model.
            sex: the sex of the individual, `None` for the sex of `simulation`.
            simulation: the key of the simulation of the individual, `None`
                for the sex of the experiment (`sex`).

        Returns:
            The lowest and the highest renal function.

        Raises:
            TypeError: if the class has no model object (`pkpd_model`).
            ValueError: if the model object has no build of the sex, or the sex
                is not given in an experiment of both sexes (`subjects_sex`).
        """
        changes = covariate_changes(bw, height, hf_severity, cp_severity, 1.0, age=age)
        outflows = kidney_outflow(cls._build(sex, simulation), changes)
        high = min(
            (o.plasma_flow / o.filtered for o in outflows if o.filtered > 0),
            default=math.inf,
        )
        if cls.renal_function_range is None:
            return 0.0, high
        low, own = cls.renal_function_range
        return low, min(own, high)

    @classmethod
    def renal_map_of(
        cls, sex: Sex | None = None, *, simulation: str | None = None
    ) -> Mapping[str, float]:
        """Get the renal classes of the subjects of a sex or of a simulation.

        `renal_map` is of the sex of the experiment; the classes of the other
        sex, e.g. of the group of women of a study of men
        (`renal_map_of(simulation="po_women")`), are those of the same preset
        at the healthy GFR of the build of that sex (`ScaledRenalMap.of_sex`).
        A map which is no `ScaledRenalMap` is not relative to a healthy GFR and
        is the map of every sex.

        Args:
            sex: the sex, `None` for the sex of `simulation`.
            simulation: the key of the simulation of the subjects (`sex_of`),
                `None` for the sex of the experiment (`sex`).

        Returns:
            The renal function of every class.

        Raises:
            TypeError: for a `ScaledRenalMap` with its own GFR (of one build)
                and another sex than the sex of the experiment.
            ValueError: if the sex is not given in an experiment of both sexes
                (`subjects_sex`).
        """
        sex = cls.subjects_sex(sex, simulation)
        renal_map = cls.renal_map
        if sex == Sex(cls.sex) or not isinstance(renal_map, ScaledRenalMap):
            return renal_map
        return renal_map.of_sex(sex)

    @classmethod
    def class_changes(
        cls,
        classification: impairment.Classification,
        key: str,
        *,
        sex: Sex | None = None,
        simulation: str | None = None,
        age: float | None = None,
    ) -> dict[str, Quantity]:
        """Get the changes of a clinical class of a condition.

        The changes set the parameter of the class (`Classification.parameter`)
        before the initialization: `CP_severity` and `HF_severity` to the value
        of the class, and for a GFR category (`eGFR`) `f_renal_function` to the
        eGFR of the class over the eGFR of the normal renal function of the
        build of the subjects (`normal_gfr`, the GFR per 1.73 m2 of its
        reference individual), so the CKD factors of the library take the
        values of the class; the category of the normal function (G1) is the
        renal function 1. An age class and a class of the body mass index set
        the covariates instead: the age class its age and the body weight and
        height of the reference individual of that age, the class of the body
        mass index the height of the reference individual and the body weight
        at the body mass index of the class (the class of normal weight, which
        has no value of its own, the reference individual).

        Args:
            classification: the classification, e.g.
                `pkpdsim.impairment.NYHA_CLASSES`.
            key: the key of the class, e.g. `III`.
            sex: the sex of the subjects of a GFR category or of a class of
                the body mass index (`subjects_sex`).
            simulation: the key of the simulation of the subjects.
            age: the age of the subjects of a GFR category or of a class of
                the body mass index (years), `None` for the adult of ICRP 89
                (35 years).

        Returns:
            The changes, `{parameter: value}`, for `preinit_changes`.

        Raises:
            KeyError: for a key which is no class of the classification.
        """
        c = classification[key]
        if classification.parameter == ids.AGE:
            ref = reference_individual(cls.subjects_sex(sex, simulation), c.value)
            return {
                ids.AGE: Q(c.value, "yr"),
                ids.BW: Q(ref.bw, "kg"),
                ids.HEIGHT: Q(ref.height, "cm"),
            }
        if classification.parameter == ids.BMI:
            ref = reference_individual(cls.subjects_sex(sex, simulation), _age(age))
            bw = ref.bw if c.value is None else c.value * (ref.height / 100.0) ** 2
            return {ids.BW: Q(bw, "kg"), ids.HEIGHT: Q(ref.height, "cm")}
        if classification.parameter != "eGFR":
            return {classification.parameter: Q(c.value, "dimensionless")}
        if c.value is None:
            return {ids.F_RENAL_FUNCTION: Q(1.0, "dimensionless")}
        normal = normal_gfr(cls.subjects_sex(sex, simulation), _age(age)) * 1e3
        return {ids.F_RENAL_FUNCTION: Q(c.value / normal, "dimensionless")}

    @classmethod
    def severity_values(
        cls,
        *,
        cp_severity: Sequence[float] | float | None = None,
        meld: Sequence[float] | float | None = None,
        egfr: Sequence[float] | float | None = None,
        renal_function: Sequence[float] | float | None = None,
        hf_severity: Sequence[float] | float | None = None,
        sex: Sex | None = None,
        simulation: str | None = None,
        age: float | None = None,
    ) -> dict[str, np.ndarray]:
        """Get the values of the severities of a continuous course of the impairments.

        The conditions are continuous: every factor of the library is
        piecewise linear in its severity between the classes (`CP_severity`,
        `HF_severity`) or the stage anchors of the eGFR (`f_renal_function`),
        and the conditions combine, the factor of a quantity is the product
        of the factors of the conditions (`fimp_<id>`), e.g. the GFR of
        cirrhosis and CKD together (hepatorenal) or of heart failure and CKD
        (cardiorenal). A severity is given directly or by its clinical
        measure: `meld` the MELD score (`cp_severity_of_meld`), `egfr` the
        eGFR in ml/min/1.73 m2 of the intrinsic renal function, which
        `f_renal_function` sets relative to the eGFR of the normal function
        of the build of the subjects (`normal_gfr`).

        Args:
            cp_severity: the severity of cirrhosis, 0 to 3.
            meld: the MELD score instead of `cp_severity`.
            egfr: the eGFR (ml/min/1.73 m2) of the renal function.
            renal_function: `f_renal_function` instead of `egfr`.
            hf_severity: the severity of heart failure, 0 to 4.
            sex: the sex of the subjects of an eGFR (`subjects_sex`).
            simulation: the key of the simulation of the subjects.
            age: the age of the subjects of an eGFR (years), whose normal eGFR
                the renal function is relative to; `None` for 35 years.

        Returns:
            The values by parameter (`CP_severity`, `f_renal_function`,
            `HF_severity`), arrays of one length, the points of a course.

        Raises:
            ValueError: for no severity, two ways of one severity, sequences
                of different lengths, a value outside of its range or not
                finite.
        """
        given = {
            k: v
            for k, v in {
                "cp_severity": cp_severity,
                "meld": meld,
                "egfr": egfr,
                "renal_function": renal_function,
                "hf_severity": hf_severity,
            }.items()
            if v is not None
        }
        if not given:
            raise ValueError(f"{cls.__name__}: give at least one severity.")
        for a, b in (("cp_severity", "meld"), ("egfr", "renal_function")):
            if a in given and b in given:
                raise ValueError(f"{cls.__name__}: give `{a}` or `{b}`, not both.")
        arrays = {
            k: np.atleast_1d(np.asarray(v, dtype=float)) for k, v in given.items()
        }
        lengths = {len(v) for v in arrays.values()}
        if len(lengths) != 1 or 0 in lengths:
            raise ValueError(
                f"{cls.__name__}: the severities of a course have one length > 0, "
                f"not {sorted(lengths)}."
            )
        for key, values in arrays.items():
            if not np.all(np.isfinite(values)):
                raise ValueError(f"{cls.__name__}: `{key}` is not finite: {values}.")
        ranges = {
            "cp_severity": (0.0, 3.0),
            "hf_severity": (0.0, 4.0),
            "egfr": (0.0, math.inf),
            "renal_function": (0.0, math.inf),
        }
        for key, (lo, hi) in ranges.items():
            if key in arrays and (arrays[key].min() < lo or arrays[key].max() > hi):
                raise ValueError(
                    f"{cls.__name__}: `{key}` {arrays[key]} is outside of {lo} to {hi}."
                )
        out: dict[str, np.ndarray] = {}
        if "meld" in arrays:
            from pkpdsim.pbpk.impairment import cp_severity_of_meld

            out[ids.CP_SEVERITY] = np.array(
                [cp_severity_of_meld(m) for m in arrays["meld"]]
            )
        elif "cp_severity" in arrays:
            out[ids.CP_SEVERITY] = arrays["cp_severity"]
        if "egfr" in arrays:
            normal = normal_gfr(cls.subjects_sex(sex, simulation), _age(age)) * 1e3
            out[ids.F_RENAL_FUNCTION] = arrays["egfr"] / normal
        elif "renal_function" in arrays:
            out[ids.F_RENAL_FUNCTION] = arrays["renal_function"]
        if "hf_severity" in arrays:
            out[ids.HF_SEVERITY] = arrays["hf_severity"]
        return out

    @classmethod
    def severity_changes(
        cls,
        *,
        cp_severity: float | None = None,
        meld: float | None = None,
        egfr: float | None = None,
        renal_function: float | None = None,
        hf_severity: float | None = None,
        sex: Sex | None = None,
        simulation: str | None = None,
        age: float | None = None,
    ) -> dict[str, Quantity]:
        """Get the changes of one point of the continuous severities of the impairments.

        The continuous counterpart of `class_changes`, for `preinit_changes`;
        the arguments are those of `severity_values`, one number each.

        Args:
            cp_severity: the severity of cirrhosis, 0 to 3.
            meld: the MELD score instead of `cp_severity`.
            egfr: the eGFR (ml/min/1.73 m2) of the renal function.
            renal_function: `f_renal_function` instead of `egfr`.
            hf_severity: the severity of heart failure, 0 to 4.
            sex: the sex of the subjects of an eGFR (`subjects_sex`).
            simulation: the key of the simulation of the subjects.
            age: the age of the subjects of an eGFR (years), `None` for 35 years.

        Returns:
            The changes, `{parameter: value}`.
        """
        values = cls.severity_values(
            cp_severity=cp_severity,
            meld=meld,
            egfr=egfr,
            renal_function=renal_function,
            hf_severity=hf_severity,
            sex=sex,
            simulation=simulation,
            age=age,
        )
        return {k: Q(float(v[0]), "dimensionless") for k, v in values.items()}

    @classmethod
    def impairment_dimension(
        cls,
        dimension: str,
        *,
        cp_severity: Sequence[float] | None = None,
        meld: Sequence[float] | None = None,
        egfr: Sequence[float] | None = None,
        renal_function: Sequence[float] | None = None,
        hf_severity: Sequence[float] | None = None,
        sex: Sex | None = None,
        simulation: str | None = None,
        age: float | None = None,
    ) -> Dimension:
        """Get the dimension of a scan over a continuous course of the impairments.

        The severities given change together along the dimension (a course,
        e.g. cirrhosis and its hepatorenal decline of the renal function); a
        scan of two dimensions (`Scan(simulation, [d1, d2])`) is the grid of
        both, e.g. every severity of cirrhosis with every eGFR. The values are
        set before the initialization of every point (a dimension without a
        time), the given measures (`meld`, `egfr`, ...) are coordinates of the
        result.

        Args:
            dimension: the id of the dimension.
            cp_severity: the severity of cirrhosis, 0 to 3.
            meld: the MELD score instead of `cp_severity`.
            egfr: the eGFR (ml/min/1.73 m2) of the renal function.
            renal_function: `f_renal_function` instead of `egfr`.
            hf_severity: the severity of heart failure, 0 to 4.
            sex: the sex of the subjects of an eGFR (`subjects_sex`).
            simulation: the key of the simulation of the subjects.
            age: the age of the subjects of an eGFR (years), `None` for 35 years.

        Returns:
            The dimension of values.
        """
        given = {
            "cp_severity": cp_severity,
            "meld": meld,
            "egfr": egfr,
            "renal_function": renal_function,
            "hf_severity": hf_severity,
        }
        values = cls.severity_values(**given, sex=sex, simulation=simulation, age=age)
        coordinates = {
            k: np.asarray(v, dtype=float) for k, v in given.items() if v is not None
        }
        return Dimension(dimension, values=values, coordinates=coordinates)

    @staticmethod
    def population_dimension(
        population: VirtualPopulation,
        sex: Sex,
        dimension: str = "dim_individual",
    ) -> Dimension:
        """Get the dimension of a scan over the individuals of a virtual population of a sex.

        Every point is an individual of the sex (`VirtualPopulation.of_sex`):
        its age, body weight, height and severities set before the
        initialization, so the physiology of the build of the sex follows.
        A study of both sexes scans each sex on its build, the simulation of
        the other sex declared in `simulation_sexes`.

        Args:
            population: the virtual population (`pkpdsim.pbpk.virtual`).
            sex: the sex of the build of the simulation.
            dimension: the id of the dimension.

        Returns:
            The dimension of values.

        Raises:
            ValueError: if the population has no individual of the sex.
        """
        units = {ids.AGE: "yr", ids.BW: "kg", ids.HEIGHT: "cm"}
        values = {
            k: Q(v, units.get(k, "dimensionless"))
            for k, v in population.values(sex).items()
        }
        return Dimension(dimension, values=values)

    @staticmethod
    def clock_changes(hour: float) -> dict[str, Quantity]:
        """Get the change which starts a simulation at a clock time.

        For a model with circadian rhythms (`PBPKDefinition.circadian`): the
        clock time `CLOCK` at time 0 of the simulation, set before the
        initialization (`preinit_changes`), e.g. a dose at 8:00 or at 20:00.

        Args:
            hour: the clock time (hours, 0 to 24).

        Returns:
            The change of `CLOCK` in hours.

        Raises:
            ValueError: for a clock time not in [0, 24).
        """
        if not 0.0 <= hour < 24.0:
            raise ValueError(f"The clock time {hour} is not in [0, 24) hours.")
        return {ids.CLOCK: Q(hour, "hr")}

    def class_simulations(
        self,
        simulation: Simulation,
        classification: impairment.Classification,
        *,
        prefix: str,
        keys: Sequence[str] | None = None,
    ) -> dict[str, Simulation]:
        """Get a simulation per clinical class of a condition.

        Every simulation is `simulation` with the changes of its class
        (`class_changes`) before the initialization; its key is
        `<prefix>_<class key>` (`to_sid`), its task `task_<key>`, which
        `add_classes` draws. A condition with a homeostasis of a PD component
        equilibrates each simulation (`equilibrated`).

        Args:
            simulation: the simulation of the healthy subjects, e.g. a dose.
            classification: the classification, e.g.
                `pkpdsim.impairment.CHILD_PUGH_CLASSES`.
            prefix: the prefix of the keys of the simulations.
            keys: the keys of the classes, `None` for every class.

        Returns:
            The simulations by key, in the order of the classes.

        Raises:
            ValueError: if `simulation` sets the parameter of the classes
                itself.
        """
        parameter = (
            ids.F_RENAL_FUNCTION
            if classification.parameter == "eGFR"
            else classification.parameter
        )
        if parameter in simulation.preinit_changes:
            raise ValueError(
                f"{type(self).__name__}: the simulation of the classes of "
                f"{classification.name} sets '{parameter}' itself, which the "
                "classes set."
            )
        out = {}
        for key in classification.keys() if keys is None else keys:
            changes = self.class_changes(classification, key)
            out[self.to_sid(f"{prefix}_{key}")] = Simulation(
                time_unit=simulation.time_unit,
                start=simulation.start,
                end=simulation.end,
                preinit_changes={**simulation.preinit_changes, **changes},
                changes=list(simulation.changes),
                presimulation=simulation.presimulation,
                times=None if simulation.times is None else list(simulation.times),
                steps=simulation.steps,
                time_shift=simulation.time_shift,
            )
        return out

    def add_classes(
        self,
        plot: Plot,
        classification: impairment.Classification,
        *,
        prefix: str,
        yid: str,
        xid: str = "time",
        keys: Sequence[str] | None = None,
    ) -> None:
        """Draw a selection of the simulations of the classes of a condition.

        One curve per class from the task of its simulation of
        `class_simulations`, labelled and colored by the class.

        Args:
            plot: the plot.
            classification: the classification of the simulations.
            prefix: the prefix of the keys of the simulations.
            yid: the selection on the y axis, e.g. `[Cve_hctz]`.
            xid: the selection on the x axis.
            keys: the keys of the classes, `None` for every class.
        """
        for key in classification.keys() if keys is None else keys:
            c = classification[key]
            plot.add_data(
                task=f"task_{self.to_sid(f'{prefix}_{key}')}",
                xid=xid,
                yid=yid,
                label=c.label,
                color=c.color,
            )


def _model_path(model: AbstractModel | Path | str | None) -> Path | None:
    """The file of a model of an experiment, `None` for SBML text or no model."""
    if isinstance(model, AbstractModel):
        return model.source.path
    if model is None or (isinstance(model, str) and model.lstrip().startswith("<")):
        return None
    return Path(model)
