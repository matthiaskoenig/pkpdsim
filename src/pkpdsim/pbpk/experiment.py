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
from pathlib import Path
from typing import TYPE_CHECKING, ClassVar, override

from sbmlsim import Q
from sbmlsim.model import AbstractModel

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
from pkpdsim.sex import Sex

if TYPE_CHECKING:
    from sbmlsim.units import Quantity

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


def _model_path(model: AbstractModel | Path | str | None) -> Path | None:
    """The file of a model of an experiment, `None` for SBML text or no model."""
    if isinstance(model, AbstractModel):
        return model.source.path
    if model is None or (isinstance(model, str) and model.lstrip().startswith("<")):
        return None
    return Path(model)
