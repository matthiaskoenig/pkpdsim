"""The starts of the species of a definition (`PBPKDefinition.starts`).

An endogenous substance starts at the equilibrium of every part with its
baseline (`Substance.baseline`), which is stationary without production and
elimination. A production in an organ moves the organ away from it: the
plasma leaving the muscle carries the production, its tissue holds more than
`Kp` times the plasma. A model states the steady state of such an organ in
`starts`, the initial concentration of a species by its flat id, and the body
starts stationary.
"""

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pytest
from sbmlutils.factory import Port, PortType

from pkpdsim.pbpk import (
    CustomReaction,
    Excretion,
    Formula,
    Kp,
    ModelParameter,
    Organ,
    Part,
    PBPKDefinition,
    Substance,
    U,
    Value,
    build,
    compose,
    ids,
)

from .conftest import dateless, individual, load

if TYPE_CHECKING:
    import roadrunner

S = "end"
DAY = 24 * 60.0
PRODUCTION = 0.0074
"""The production of the muscle (mmole/min)."""

#: the plasma flow of the muscle and the exchange of its tissue (l/min)
_QPL_MU = "Qmu * (1 dimensionless - HCT)"
_PS_MU = f"ftissue_{S} * Vmu * (1 dimensionless - f_vas_mu) / Vti_total_ref_{S}"
#: the steady state of the muscle with the production in its tissue (mM): the
#: plasma leaving the muscle carries the production, the tissue is in
#: equilibrium with the plasma plus the production over the exchange
MUSCLE_PLASMA = f"baseline_{S} + k_prod / ({_QPL_MU})"
MUSCLE_TISSUE = f"Kp_{S} * ({MUSCLE_PLASMA} + k_prod / ({_PS_MU}))"


def _definition(starts: dict[str, object] | None = None) -> PBPKDefinition:
    """A substance produced in the muscle and filtered, at the baseline of its steady state.

    The baseline `k_prod/GFR` is the arterial plasma of the steady state (the
    filtration clears `GFR*Car`); without `starts` every species starts at
    the equilibrium with it.
    """
    fields: dict[str, object] = {} if starts is None else {"starts": starts}
    return PBPKDefinition(
        sid="produced",
        substances=[
            Substance(
                S,
                name="endogenous substance",
                mr=113.1,
                distribution=Kp(0.6, ftissue=6.0),
                excretion=Excretion(urine=True),
                baseline=Formula("k_prod / GFR", U.mM),
            )
        ],
        chemistry=[
            CustomReaction(
                "PROD", Organ.MUSCLE, Part.TISSUE, f"-> Cti_{S}", "k_prod", []
            )
        ],
        parameters=[ModelParameter("k_prod", PRODUCTION, U.mmole_per_min)],
        **fields,  # ty: ignore[invalid-argument-type]
    )


STARTS: dict[str, object] = {
    f"MU__Cpl_{S}": Formula(MUSCLE_PLASMA, U.mM),
    f"MU__Cti_{S}": Formula(MUSCLE_TISSUE, U.mM),
}
"""The steady state of the muscle."""


def _drift(r: "roadrunner.RoadRunner", end: float = DAY) -> dict[str, float]:
    """The largest relative change of every floating species over `end` (min)."""
    sids = list(r.model.getFloatingSpeciesIds())
    start = np.array([r[f"[{sid}]"] for sid in sids])
    values = np.asarray(r.simulate(0, end, 97, selections=[f"[{sid}]" for sid in sids]))
    scale = np.where(start > 0, start, 1.0)
    return dict(zip(sids, np.max(np.abs(values - start), axis=0) / scale, strict=True))


def _sinks(sid: str) -> bool:
    """A sink accumulates: the urine of the substance."""
    return sid == ids.urine(S)


@pytest.mark.parametrize(
    "changes", [{}, {"BW": 110.0}, {"f_renal_function": 0.3}], ids=str
)
def test_the_starts_make_a_production_stationary(
    tmp_path: Path, changes: dict[str, float]
) -> None:
    """With the steady state of the muscle as starts every species stays at its start.

    For the default individual, a heavier one and an impaired kidney: the
    starts are initial assignments which follow a change before the
    simulation.
    """
    r = individual(build(_definition(STARTS), tmp_path).flat_path, **changes)
    drift = _drift(r)
    moved = {sid: d for sid, d in drift.items() if d > 1e-6 and not _sinks(sid)}
    assert not moved, moved
    assert r[f"[MU__Cpl_{S}]"] > r[f"[Car_{S}]"]


def test_without_starts_a_production_is_not_stationary(tmp_path: Path) -> None:
    """At the equilibrium with the baseline the plasma moves (the reason of `starts`)."""
    r = load(build(_definition(), tmp_path))
    drift = _drift(r)
    assert drift[f"MU__Cti_{S}"] > 0.01
    assert drift[f"Car_{S}"] > 0.01


def test_a_start_is_the_initial_assignment_of_its_species(tmp_path: Path) -> None:
    """A start replaces the start at the equilibrium; the others stay.

    The muscle species start at the formulas, read through input ports of the
    muscle (`k_prod`, `Qmu`, `HCT`); the plasma of the liver at the baseline.
    """
    r = load(build(_definition(STARTS), tmp_path))
    b = r[ids.baseline(S)]
    qpl = r["Qmu"] * (1 - r["HCT"])
    assert r[f"[MU__Cpl_{S}]"] == pytest.approx(b + PRODUCTION / qpl, rel=1e-12)
    assert r[f"[LI__Cpl_{S}]"] == pytest.approx(b, rel=1e-12)
    assert r[f"[LI__Cti_{S}]"] == pytest.approx(0.6 * b, rel=1e-12)
    model = compose(_definition(STARTS))
    muscle = next(md for md in model.model_definitions if md.sid.endswith("_mu"))
    inputs = {
        str(p.sid)
        for p in muscle.parameters or []
        if isinstance(p.port, Port) and p.port.portType == PortType.INPUT_PORT
    }
    assert {"k_prod", "Qmu", "HCT", ids.baseline(S)} <= inputs


def test_a_start_may_be_a_number(tmp_path: Path) -> None:
    """A number is the initial concentration of the species (mM)."""
    r = load(build(_definition({f"MU__Cti_{S}": 0.2}), tmp_path))
    assert r[f"[MU__Cti_{S}]"] == pytest.approx(0.2, rel=1e-12)
    r = load(
        build(
            _definition({f"MU__Cti_{S}": Value(0.3, U.mM, notes="a test")}),
            tmp_path / "value",
        )
    )
    assert r[f"[MU__Cti_{S}]"] == pytest.approx(0.3, rel=1e-12)


def test_a_start_of_a_pool(tmp_path: Path) -> None:
    """A species of the top model (a pool) takes its start as well."""
    r = load(
        build(
            _definition({f"Cve_{S}": Formula(f"2 dimensionless * baseline_{S}", U.mM)}),
            tmp_path,
        )
    )
    assert r[f"[Cve_{S}]"] == pytest.approx(2 * r[ids.baseline(S)], rel=1e-12)


def test_a_start_of_no_species_is_refused() -> None:
    """A start which names no species of the model raises, naming all at once.

    The starts of the organs and of the top model are collected before the
    error: a start of the muscle, of the liver and of the top model.
    """
    with pytest.raises(ValueError, match="MU__Cxx_end") as error:
        compose(_definition({"MU__Cxx_end": 0.1, "LI__Czz_end": 0.1, "Cyy_end": 0.1}))
    assert "LI__Czz_end" in str(error.value)
    assert "'Cyy_end'" in str(error.value)


def test_a_start_is_validated() -> None:
    """A start in another unit than mM, a negative start and a key which is no SId."""
    with pytest.raises(ValueError, match="MU__Cti_end") as error:
        _definition(
            {
                f"MU__Cti_{S}": Value(0.1, U.mmole_per_min),
                f"MU__Cpl_{S}": -1.0,
                "1x": 0.1,
            }
        ).validate()
    text = str(error.value)
    assert "MU__Cpl_end" in text
    assert "1x" in text


def test_a_start_formula_reads_the_top_model() -> None:
    """A formula of a start reads constants of the top model; an unknown id raises."""
    with pytest.raises(ValueError, match="unknown_id"):
        compose(_definition({f"MU__Cti_{S}": Formula("unknown_id * 1 mM", U.mM)}))


@pytest.mark.parametrize(
    "key", [f"MU__Cti_{S}", f"Cve_{S}"], ids=["organ species", "pool"]
)
def test_a_start_formula_refuses_a_state(key: str) -> None:
    """A formula of a start which reads what changes during the simulation raises.

    A species (`Car_<s>`) or an observer (`Carm_<s>`) would be read at its
    value before the species start.
    """
    for state in (f"Car_{S}", f"Carm_{S}"):
        with pytest.raises(ValueError, match=state):
            compose(_definition({key: Formula(f"2 dimensionless * {state}", U.mM)}))


def test_without_starts_the_model_is_unchanged(tmp_path: Path) -> None:
    """Empty starts give the SBML of a definition without the field."""
    a = build(_definition(), tmp_path / "a").flat_path
    b = build(_definition({}), tmp_path / "b").flat_path
    assert dateless(a.read_text()) == dateless(b.read_text())
