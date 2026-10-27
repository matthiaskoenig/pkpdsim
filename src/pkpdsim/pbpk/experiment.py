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

The module imports sbmlsim; `pkpdsim.pbpk` exports its names lazily, so
`import pkpdsim.pbpk` (a factory) imports no simulator.
"""

from __future__ import annotations

from collections.abc import Callable, Iterator, Mapping
from functools import cache
from pathlib import Path
from typing import TYPE_CHECKING, ClassVar, override

from sbmlsim import Q
from sbmlsim.simulation import Simulation
from sbmlsim.simulator import Simulator

from pkpdsim import impairment
from pkpdsim.experiment import PKPDSimulationExperiment
from pkpdsim.pbpk import ids
from pkpdsim.scan import start_values

if TYPE_CHECKING:
    from sbmlsim.units import Quantity

#: the healthy GFR the presets of `pkpdsim.impairment` are relative to [ml/min]
REFERENCE_GFR = 101.0


def healthy_gfr(
    model_path: Path,
    *,
    bw: float | None = None,
    height: float | None = None,
    hf_severity: float = 0.0,
    cp_severity: float = 0.0,
) -> float:
    """Get the healthy glomerular filtration rate of an individual [ml/min].

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
        The GFR at renal function 1 [ml/min].
    """
    path = Path(model_path)
    return _healthy_gfr(
        path,
        path.stat().st_mtime_ns,
        bw,
        height,
        float(hf_severity),
        float(cp_severity),
    )


@cache
def _healthy_gfr(
    model_path: Path,
    mtime_ns: int,  # a key of the cache only
    bw: float | None,
    height: float | None,
    hf_severity: float,
    cp_severity: float,
) -> float:
    """Simulate the healthy GFR of an individual once, see `healthy_gfr`.

    `mtime_ns`, the modification time of the model file, is only a key of the
    cache: a model which is created again is simulated again.
    """
    changes = {
        ids.HF_SEVERITY: Q(hf_severity, "dimensionless"),
        ids.CP_SEVERITY: Q(cp_severity, "dimensionless"),
    }
    if bw is not None:
        changes[ids.BW] = Q(bw, "kg")
    if height is not None:
        changes[ids.HEIGHT] = Q(height, "cm")
    result = Simulator(n_workers=1).run(
        model_path,
        Simulation(end=1, time_unit="min", preinit_changes=changes),
    )
    gfr = float(start_values(result, ids.GFR_BASE).item())
    return float(Q(gfr, result.units[ids.GFR_BASE]).to("ml/min").magnitude)


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
        gfr: Callable[[], float] | None = None,
        reference: float = REFERENCE_GFR,
    ) -> None:
        """Get the classes of a preset at a healthy GFR.

        Args:
            preset: the renal function of the classes relative to `reference`.
            gfr: the healthy GFR [ml/min], called on the first access to a
                value; `None` for the `healthy_gfr` of the class the map is an
                attribute of.
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
            factor = self._reference / self._gfr()
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

    Attributes:
        renal_function_range: the renal functions the model supports, see
            `renal_function`: up to 2 the filtration stays a physiological
            fraction of the renal plasma flow of the library, at 0 (anuria)
            the substances leave by the other routes alone; a model whose
            urine is not defined without filtration starts above 0.
    """

    renal_map: ClassVar[Mapping[str, float]] = ScaledRenalMap(
        impairment.RENAL_FUNCTION_GFR
    )
    renal_colors: ClassVar[Mapping[str, str]] = impairment.RENAL_COLORS_END_STAGE
    cirrhosis_map: ClassVar[Mapping[str, float]] = impairment.CHILD_PUGH
    cirrhosis_colors: ClassVar[Mapping[str, str]] = impairment.CIRRHOSIS_COLORS
    cardiac_map: ClassVar[Mapping[str, float]] = impairment.HEART_FAILURE
    cardiac_colors: ClassVar[Mapping[str, str]] = impairment.HEART_FAILURE_COLORS
    renal_function_range: ClassVar[tuple[float, float]] = (0.0, 2.0)

    @classmethod
    def healthy_gfr(
        cls,
        *,
        bw: float | None = None,
        height: float | None = None,
        hf_severity: float = 0.0,
        cp_severity: float = 0.0,
    ) -> float:
        """Get the healthy glomerular filtration rate of an individual [ml/min].

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
            The GFR at renal function 1 [ml/min].

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
            ValueError: for a renal function outside of `renal_function_range`.
        """
        value = clearance / cls.healthy_gfr(
            bw=bw, height=height, hf_severity=hf_severity, cp_severity=cp_severity
        )
        low, high = cls.renal_function_range
        if not low <= value <= high:
            raise ValueError(
                f"The renal function {value:.3g} of a clearance of {clearance} "
                f"ml/min is outside of the supported range {low} to {high}."
            )
        return Q(value, "dimensionless")
