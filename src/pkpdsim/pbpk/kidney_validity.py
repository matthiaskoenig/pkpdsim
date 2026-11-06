"""The simulated individual and the validity of the kidney outflow.

The plasma leaving the kidney lacks the net filtrate, so the outflow of the
kidney is valid while the filtered plasma water is below the plasma flow of
the kidney, `GFR < Qki*(1 - HCT)` (the glomeruli filter the water of the
plasma which flows through them, independent of the substances), and the net
filtered plasma flow of every substance is below the flow which carries it,
`fu*GFR*(1 - f_reabsorbed) < Qki*(1 - HCT)` (`Qki*BP` with erythrocytes in
equilibrium, `fu*GFR` with a threshold reabsorption, whose fraction falls to
0 at a high concentration). The water bounds the renal function of every
model which filters a substance (a PD component which filters the water at
the GFR, the fluid of hctz, as well); with `fu <= 1` and `BP >= 1 - HCT` it is
the tighter bound of the two. A model without a filtered substance is not
simulated and has no bound. GFR, flow,
hematocrit and fu are initial assignments of the parameters of the
individual, so `kidney_outflow` reads them for an individual, the changes of
a simulation before its initialization, from a simulation of one minute
(`individual_values`, which `healthy_gfr` reads the GFR with as well), and
`simulated_individuals` gives the individuals a simulation or a scan of an
experiment simulates. `PBPKSimulationExperiment`
(`pkpdsim.pbpk.experiment`) checks them all when it is initialized.

What makes an individual is the changes of the ids which the flows depend on
through initial assignments and assignment rules of the flat model, so a dose
makes no individual. The values are cached per model file, its modification
time and the individual, and the model is parsed and loaded once per process
(`pkpdsim.pbpk.flat_model`, shared with the equilibration). The check
reads the start of a simulation: its changes before the initialization (of
its presimulation as well) and the values of the dimensions of a scan without
a time. It does not check a change during a simulation (`Simulation.changes`,
also at its start, of an id the definition names `dynamic`), a dimension with
a time (`at`), a scan over models and the draws of the fitted parameters of a
run over them (`pkpdsim.uncertainty`).

The module imports sbmlsim; `pkpdsim.pbpk` exports its names
lazily, so `import pkpdsim.pbpk` (a factory) imports no simulator.
"""

from __future__ import annotations

import itertools
from collections.abc import Mapping
from dataclasses import dataclass
from functools import cache
from pathlib import Path

from sbmlsim import Q
from sbmlsim.result.scan import start_values
from sbmlsim.simulation import Scan, Simulation
from sbmlsim.simulation.definition import Value as Change
from sbmlsim.simulation.scan import DimensionKind
from sbmlsim.simulator import Simulator
from sbmlsim.units import Quantity

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.flat_model import flat_model, loaded_model
from pkpdsim.pbpk.organs import Organ
from pkpdsim.uncertainty import UNCERTAINTY_DIMENSION

__all__ = [
    "IndividualValues",
    "KidneyOutflow",
    "changes_text",
    "covariate_changes",
    "individual_values",
    "kidney_outflow",
    "simulated_individuals",
]


def covariate_changes(
    bw: float | None,
    height: float | None,
    hf_severity: float,
    cp_severity: float,
    renal_function: float | None = None,
) -> dict[str, Change]:
    """Get the changes before the initialization of an individual by its covariates.

    Args:
        bw: the body weight [kg], `None` for that of the model.
        height: the height [cm], `None` for that of the model.
        hf_severity: the severity of heart failure `HF_severity` (0 to 4).
        cp_severity: the severity of cirrhosis `CP_severity` (0 to 3).
        renal_function: the renal function `f_renal_function`, `None` for
            that of the model.

    Returns:
        The changes, quantities.
    """
    changes: dict[str, Change] = {
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
    """The net filtered plasma flow of a substance or the plasma water, and the plasma flow of the kidney.

    The outflow of the kidney is valid while the net filtered plasma flow is
    below the plasma flow which carries the substance through the kidney;
    the filtered plasma water (`substance` `None`, the `formula` `GFR`) is
    valid while the GFR is below the plasma flow of the kidney.

    Attributes:
        substance: the substance, filtered into the urine; `None` for the
            water of the plasma, which the glomeruli filter at the GFR.
        gfr: the GFR of the individual [ml/min].
        filtered: the net filtered plasma flow [ml/min], see `formula`.
        plasma_flow: the flow which carries the substance through the kidney,
            `Qki*(1 - HCT)`, `Qki*BP` with erythrocytes in equilibrium
            [ml/min].
        formula: the formula of the net filtered plasma flow:
            `fu*GFR*(1 - f_reabsorbed)`, `fu*GFR` with a threshold
            reabsorption (its fraction falls to 0 at a high concentration),
            `GFR` of the water.
    """

    substance: str | None
    gfr: float
    filtered: float
    plasma_flow: float
    formula: str = "fu*GFR*(1 - f_reabsorbed)"

    @property
    def valid(self) -> bool:
        """Whether the net filtered plasma flow is below the plasma flow."""
        return self.filtered < self.plasma_flow

    def text(self) -> str:
        """Describe the flows, for a message."""
        if self.substance is None:
            return (
                f"the filtered plasma water {self.formula} {self.filtered:.4g} "
                f"ml/min against the plasma flow of the kidney "
                f"{self.plasma_flow:.4g} ml/min (Qki*(1 - HCT))"
            )
        return (
            f"the net filtered plasma flow of '{self.substance}' {self.formula} "
            f"{self.filtered:.4g} ml/min at the GFR {self.gfr:.4g} ml/min against "
            f"the plasma flow of the kidney {self.plasma_flow:.4g} ml/min "
            "(Qki*(1 - HCT), Qki*BP with erythrocytes in equilibrium)"
        )


def kidney_outflow(
    model_path: Path, changes: Mapping[str, Change] | None = None
) -> tuple[KidneyOutflow, ...]:
    """Get the net filtered plasma flows and the plasma flow of the kidney of an individual.

    The filtered plasma water first (`substance` `None`, `GFR` against
    `Qki*(1 - HCT)`), then one `KidneyOutflow` for every substance which the
    kidney filters (its filtration reaction is in the model), read from the
    start values of the individual (`individual_values`). A model without a
    filtered substance has no outflow of its kidney to check and is not
    simulated (its GFR is no flow of the library), nor is a model without the
    GFR, the flow of the kidney and the hematocrit (no model of the library).

    Args:
        model_path: the flat SBML model of the library.
        changes: the changes of the individual before the initialization
            (`Simulation.preinit_changes`), numbers in the units of the model
            or quantities.

    Returns:
        The flows of the water and of every filtered substance, in the order
        of the model.
    """
    path = Path(model_path)
    model = flat_model(path)
    if not model.water or not model.filtered:
        return ()
    values = individual_values(path, changes or {})
    flow = values.ml_per_min(ids.organ_flow(Organ.KIDNEY))
    hct = values.number(ids.HCT)
    gfr = values.ml_per_min(ids.GFR)
    outflows = [KidneyOutflow(None, gfr, gfr, flow * (1.0 - hct), formula="GFR")]
    for f in model.filtered:
        carrier = 1.0 - hct if f.kinetic else values.number(ids.bp(f.sid))
        fu_gfr = values.number(ids.fu(f.sid)) * gfr
        if f.threshold:
            filtered, formula = fu_gfr, "fu*GFR (threshold reabsorption)"
        else:
            reabsorbed = values.number(ids.fraction_reabsorbed(f.sid))
            filtered, formula = fu_gfr * (1.0 - reabsorbed), KidneyOutflow.formula
        outflows.append(
            KidneyOutflow(f.sid, gfr, filtered, flow * carrier, formula=formula)
        )
    return tuple(outflows)


@dataclass(frozen=True)
class IndividualValues:
    """The start values of an individual, in the units of the model.

    Attributes:
        numbers: id -> the value at the start.
        units: id -> its unit.
    """

    numbers: Mapping[str, float]
    units: Mapping[str, str]

    def number(self, sid: str) -> float:
        """Get the value of an id in the unit of the model."""
        return self.numbers[sid]

    def ml_per_min(self, sid: str) -> float:
        """Get the value of a flow in ml/min."""
        return float(Q(self.numbers[sid], self.units[sid]).to("ml/min").magnitude)


def individual_values(
    model_path: Path, changes: Mapping[str, Change]
) -> IndividualValues:
    """Get the start values of the flows of the kidney outflow and `GFR_base` of an individual.

    The values are initial assignments, so they are read from a simulation
    of one minute with the changes before the initialization, which
    evaluates the initial assignments again. Only the changes of the ids
    which the values depend on make the individual (`simulated_individuals`),
    and the values are cached per model file, its modification time and the
    individual: a model file which is created again is simulated again.

    Args:
        model_path: the flat SBML model of the library.
        changes: the changes of the individual before the initialization.

    Returns:
        The values of `GFR`, `GFR_base`, `Qki`, `HCT` and of the filtered
        substances `fu_<s>`, `BP_<s>` and `f_reabsorbed_<s>`.
    """
    path = Path(model_path)
    mtime = path.stat().st_mtime_ns
    inputs = flat_model(path).inputs
    relevant = {k: v for k, v in _plain(changes).items() if k in inputs}
    return _individual_values(path, mtime, _frozen(relevant))


def simulated_individuals(
    simulation: Simulation | Scan, model_path: Path
) -> list[dict[str, Change]]:
    """Get the individuals a simulation or a scan simulates on a model, each once.

    An individual is the changes before the initialization of the ids which
    the flows of the validity depend on, at the start of the simulation: of
    its presimulation and of itself, of every simulation of a dimension of
    simulations and of every point of the dimensions of values without a time
    (`Simulation.with_values`). A dimension with a time, the draws of the
    fitted parameters (`pkpdsim.uncertainty`) and the changes during a
    simulation are not read, and a scan over models has no individual of this
    model.

    Args:
        simulation: a simulation or a scan of an experiment.
        model_path: the flat SBML model the simulation runs.

    Returns:
        The changes of every individual.
    """
    scan = Scan.of(simulation)
    if scan.dimension(DimensionKind.MODELS) is not None:
        return []
    path = Path(model_path)
    inputs = flat_model(path).inputs
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
                changes.update(_plain(start.presimulation.preinit_changes))
            changes.update(_plain(start.preinit_changes))
            relevant = {k: v for k, v in changes.items() if k in inputs}
            individuals.setdefault(_frozen(relevant), relevant)
    return list(individuals.values())


def changes_text(changes: Mapping[str, Change]) -> str:
    """Describe the changes of an individual, for a message."""
    if not changes:
        return "none"
    return ", ".join(f"{k} = {v}" for k, v in sorted(changes.items()))


def _plain(changes: Mapping[str, object]) -> dict[str, Change]:
    """The plain values of changes: numbers, quantities and formulas.

    sbmlsim keeps an uncertain value apart from the plain ones of a
    simulation, but the type of the changes of a presimulation admits it; it
    is a draw, not a value of the individual.
    """
    return {
        k: v for k, v in changes.items() if isinstance(v, float | int | str | Quantity)
    }


#: a change of an individual, hashable: the target and its value, a number, a
#: formula or the magnitude and the unit of a quantity
_Frozen = tuple[str, float | str | tuple[float, str]]


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
) -> IndividualValues:
    """Simulate an individual for one minute and read its start values, once."""
    quantities = flat_model(model_path).quantities
    result = Simulator(n_workers=1).run(
        loaded_model(model_path),
        Simulation(end=1, time_unit="min", preinit_changes=_thawed(changes)),
    )
    numbers = {sid: float(start_values(result, sid).item()) for sid in quantities}
    return IndividualValues(
        numbers, {sid: str(result.units[sid]) for sid in quantities}
    )
