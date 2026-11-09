"""The observers of a zonated organ: means, outflow, amounts and the mass balance.

The tests use a tracer which does not react in the engine itself; where the
cells consume it, the test adds a first-order sink to the model and declares it
to the observers as the consumption, as the cell model of a composed liver does.
"""

import re
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest
import roadrunner
from pbpk.conftest import individual
from sbmlutils.factory import Parameter, Reaction, create_model

from pkpdsim.pbpk.substance import Substance
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Value
from pkpdsim.zonation import ZonationBuilder, ids, observers
from pkpdsim.zonation.flow import ImposedFlow, PressureFlow
from pkpdsim.zonation.geometry import Layer
from pkpdsim.zonation.observers import (
    add_observers,
    observer_objects,
    profile_ids,
)
from pkpdsim.zonation.settings import BileTransport, Geometry, ZonationSettings
from pkpdsim.zones import third_weights

from .conftest import mass_g

pytestmark = pytest.mark.filterwarnings("ignore:.*numerical dispersion:UserWarning")

DRG = Substance("drg", "drug", Value(180.0, U.g_per_mole))
SINK = 0.5
"""The first-order rate (1/min) with which the cells consume the drug."""
PERMEABILITY = 2e-6
"""The permeability (m/s) of the endothelial cell and the hepatocyte."""


def _settings(n: int = 10, **kwargs: object) -> ZonationSettings:
    """The settings of a chain of `n` zones of a 1 kg liver."""
    kwargs.setdefault("flow", ImposedFlow(1.0))
    return ZonationSettings(
        n_zones=n,
        liver_mass=mass_g(1000.0),
        **kwargs,  # ty: ignore[invalid-argument-type]
    )


def _sink_layer(settings: ZonationSettings) -> Layer | str:
    """The layer the cells consume from: the hepatocyte, the vascular site of the well-stirred."""
    if settings.geometry is Geometry.WELL_STIRRED:
        return ids.VASCULAR
    return Layer.HEPATOCYTE


def _consumed(settings: ZonationSettings) -> str:
    """The formula (mmole/min) of the consumption of the drug, the sum of the zone sinks."""
    return " + ".join(f"R_sink_{ids.zone(k)}" for k in range(1, settings.n_zones + 1))


def _model_path(
    tmp_path: Path,
    settings: ZonationSettings,
    *,
    watched: bool = True,
    sink: float = SINK,
    declared: bool = True,
) -> Path:
    """Write the model of the settings with the sink of the cells and the observers."""
    builder = ZonationBuilder(settings, [DRG], None)
    model = builder.build()
    layer = _sink_layer(settings)
    if sink:
        model.parameters.append(
            Parameter(
                "k_sink",
                sink,
                unit=U.per_min,
                name="first-order consumption by the cells",
                constant=True,
            )
        )
        for k in range(1, settings.n_zones + 1):
            zone = ids.zone(k)
            c = ids.species("drg", layer, zone)
            v = ids.volume(layer, zone)
            model.reactions.append(
                Reaction(
                    f"R_sink_{zone}",
                    f"{c} ->",
                    formula=(f"k_sink * {c} * {v}", U.mmole_per_min),
                    compartment=v,
                    name=f"consumption of the drug in zone {k}",
                )
            )
    if watched:
        consumption = {"drg": _consumed(settings)} if sink and declared else None
        add_observers(model, observer_objects(settings, [DRG], consumption))
    path = tmp_path / "observed.xml"
    create_model(model, filepath=path, sbml_level=3, sbml_version=2, validate=False)
    return path


def _runner(
    tmp_path: Path, settings: ZonationSettings, **kwargs: object
) -> roadrunner.RoadRunner:
    """The roadrunner model with the transcellular paths open and tight tolerances."""
    path = _model_path(tmp_path, settings, **kwargs)  # ty: ignore[invalid-argument-type]
    changes = {}
    if settings.geometry is Geometry.SINUSOID:
        changes = {
            "P_endo_in_drg": PERMEABILITY,
            "P_endo_out_drg": PERMEABILITY,
            "P_hep_drg": PERMEABILITY,
        }
    r = individual(path, **changes)
    r.integrator.absolute_tolerance = 1e-14
    r.integrator.relative_tolerance = 1e-10
    return r


def _layer_tokens(r: roadrunner.RoadRunner) -> list[str]:
    """The layers of the drug in the first zone, read from the species of the model."""
    return [
        sid.removeprefix("drg_").removesuffix("_z1")
        for sid in r.model.getFloatingSpeciesIds()
        if sid.startswith("drg_") and sid.endswith("_z1")
    ]


def _load(
    r: roadrunner.RoadRunner, settings: ZonationSettings, c: float = 1.0
) -> float:
    """Put `c` mM into every layer of every zone and record it as the dose, return mmole.

    The dose is computed from the compartment volumes of the model and added
    to the ledger the observers hold for what is put into the organ by hand.
    """
    dose = 0.0
    for layer in _layer_tokens(r):
        for k in range(1, settings.n_zones + 1):
            zone = ids.zone(k)
            r[f"[{ids.species('drg', layer, zone)}]"] = c
            dose += c * r[ids.volume(layer, zone)]
    observers.sync_dose(r, [DRG])
    assert r[observers.ledger("dose", "drg")] == pytest.approx(dose, rel=1e-12, abs=0)
    return dose


def _plasma(settings: ZonationSettings) -> str:
    """The layer of the blood the means are of."""
    return ids.VASCULAR if settings.geometry is Geometry.WELL_STIRRED else Layer.PLASMA


def _zone_values(r: roadrunner.RoadRunner, settings: ZonationSettings) -> list[float]:
    """The concentration in the blood of every zone, periportal first."""
    return [
        r[f"[{ids.species('drg', _plasma(settings), ids.zone(k))}]"]
        for k in range(1, settings.n_zones + 1)
    ]


def _run(
    tmp_path: Path, n: int = 10, horizon: float = 5.0, **kwargs: object
) -> tuple[roadrunner.RoadRunner, ZonationSettings]:
    """A chain with an inflow of 1 mM and a drug which the cells consume."""
    settings = _settings(n, **kwargs)
    r = _runner(tmp_path, settings)
    r["drg_" + str(_plasma(settings)) + "_in"] = 1.0
    r.simulate(0, horizon, 200)
    return r, settings


def _weights_in(objects: object, site: str) -> list[tuple[int, float]]:
    """The weights of the zones in the rule of the mean of the drug at a site."""
    rule = next(rule for rule in objects.rules if rule.variable == f"drg_{site}")  # ty: ignore[unresolved-attribute]
    terms = re.findall(
        r"(?:([0-9.e-]+) dimensionless \* )?drg_plasma_z(\d+)", rule.value
    )
    return [(int(k), float(w) if w else 1.0) for w, k in terms]


# --- the thirds are the perfusion library's ---------------------------------


@pytest.mark.parametrize("n", [1, 2, 3, 7, 10])
@pytest.mark.parametrize(("site", "periportal"), [("pp", True), ("pc", False)])
def test_the_thirds_reuse_the_perfusion_weighting(
    n: int, site: str, periportal: bool
) -> None:
    """The weights in the emitted rule are those of `third_weights`, at awkward n as well."""
    found = _weights_in(observer_objects(_settings(n=n), (DRG,)), site)
    expected = [(k, float(w)) for k, w in third_weights(n, periportal=periportal)]
    assert found == expected
    assert sum(w for _, w in found) == pytest.approx(1.0, abs=1e-15, rel=0)


def test_the_weights_come_from_the_function_of_the_library(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Replacing `third_weights` changes the rule: the module does not carry its own."""
    monkeypatch.setattr(
        observers, "third_weights", lambda n, periportal: [(2, Fraction(1))]
    )
    found = _weights_in(observer_objects(_settings(n=10), (DRG,)), "pp")
    assert found == [(2, 1.0)]


def _exact_third_mean(values: list[float], periportal: bool) -> float:
    """The mean of a piecewise constant profile over a third, by integration.

    The profile is `values[k - 1]` on `[(k - 1)/n, k/n]`; the third is
    integrated piece by piece, without the weights of the library.
    """
    n = len(values)
    lo, hi = (0.0, 1 / 3) if periportal else (2 / 3, 1.0)
    total = 0.0
    for k, v in enumerate(values):
        a, b = max(lo, k / n), min(hi, (k + 1) / n)
        total += v * max(0.0, b - a)
    return total * 3.0


@pytest.mark.parametrize("n", [3, 7, 10])
def test_the_means_integrate_the_profile_over_the_thirds(
    tmp_path: Path, n: int
) -> None:
    """A profile of 1 to n is averaged as its integral over the first and last third."""
    settings = _settings(n)
    r = _runner(tmp_path, settings)
    for k in range(1, n + 1):
        r[f"[drg_plasma_z{k}]"] = float(k)
    values = [float(k) for k in range(1, n + 1)]
    assert r["drg_pp"] == pytest.approx(_exact_third_mean(values, True), rel=1e-12)
    assert r["drg_pc"] == pytest.approx(_exact_third_mean(values, False), rel=1e-12)


# --- the means under a gradient ---------------------------------------------


def test_the_cells_consuming_the_drug_lower_the_pericentral_mean(
    tmp_path: Path,
) -> None:
    """With consumption the pericentral mean is below the periportal one, both bracketed."""
    r, settings = _run(tmp_path)
    zones = _zone_values(r, settings)
    pp, pc = r["drg_pp"], r["drg_pc"]
    assert pc < pp
    assert pp - pc > 0.01 * pp  # a gradient, not a rounding difference
    assert min(zones) <= pc <= max(zones)
    assert min(zones) <= pp <= max(zones)
    assert zones[0] > zones[-1]  # periportal first


def test_without_consumption_the_two_means_are_equal(tmp_path: Path) -> None:
    """Nothing consumed: the steady state is flat, so any weighting which sums to 1 gives one number."""
    settings = _settings(10)
    r = _runner(tmp_path, settings, sink=0.0)
    r["drg_plasma_in"] = 1.0
    r["drg_rbc_in"] = 1.0
    r.simulate(0, 200, 400)
    zones = _zone_values(r, settings)
    assert max(zones) - min(zones) < 1e-6
    assert r["drg_pp"] == pytest.approx(r["drg_pc"], rel=1e-6, abs=0)
    assert r["drg_pp"] == pytest.approx(1.0, rel=1e-6, abs=0)


# --- the outflow is not the pericentral mean ---------------------------------


def test_the_outflow_concentration_is_not_the_pericentral_mean(
    tmp_path: Path,
) -> None:
    """Under a gradient the concentration leaving the last zone differs from the mean of its third."""
    r, _ = _run(tmp_path)
    assert r["c_outflow_drg"] != pytest.approx(r["drg_pc"], rel=1e-3, abs=0)
    assert (
        r["c_outflow_drg"] < r["drg_pc"]
    )  # the last zone is the lowest, the mean includes zone 8 and 9


@pytest.mark.parametrize(
    "flow",
    [ImposedFlow(1.5), PressureFlow(), PressureFlow(p_onc=2000.0)],
    ids=["imposed", "pressure", "pressure_asymmetric"],
)
def test_the_outflow_concentration_is_the_flux_over_the_flow(
    tmp_path: Path, flow: ImposedFlow | PressureFlow
) -> None:
    """The concentration leaving is the flux over the flow of the last zone, flows included."""
    n = 5
    r = _runner(tmp_path, _settings(n, flow=flow))
    r["drg_plasma_in"] = 1.0
    r["drg_rbc_in"] = 1.0
    r.simulate(0, 5, 100)
    z = ids.zone(n)

    def carrier_flow(layer: Layer) -> float:
        sid = (
            ids.layer_flow(layer)
            if isinstance(flow, ImposedFlow)
            else ids.zone_flow(layer, z)
        )
        return float(r[sid])

    flux = sum(
        carrier_flow(layer) * r[f"[drg_{layer}_{z}]"]
        for layer in (Layer.PLASMA, Layer.RBC)
    )
    q = sum(carrier_flow(layer) for layer in (Layer.PLASMA, Layer.RBC))
    print(f"MEASURED {flow!r:.40} flux {flux!r} flow {q!r}")
    assert r["c_outflow_drg"] == pytest.approx(flux / q, rel=1e-12, abs=0)
    assert r[f"[drg_plasma_{z}]"] != pytest.approx(r[f"[drg_rbc_{z}]"], rel=1e-3)


def test_the_amounts_follow_the_volumes_of_their_zones(tmp_path: Path) -> None:
    """With a different volume in every zone, a layer amount sums c_k * V_k of the same zone."""
    n = 4
    r = _runner(tmp_path, _settings(n))
    for layer in _layer_tokens(r):
        for k in range(1, n + 1):
            r[f"V_{layer}_z{k}"] = 0.01 * (k + 1) * (1 + len(layer))
            r[f"[drg_{layer}_z{k}]"] = 0.2 * k + 0.1 * len(layer)
        expected = sum(
            r[f"[drg_{layer}_z{k}]"] * r[f"V_{layer}_z{k}"] for k in range(1, n + 1)
        )
        assert r[observers.layer_amount("drg", layer)] == pytest.approx(
            expected, rel=1e-12, abs=0
        )


def test_the_dose_of_the_initial_concentrations_follows_the_zones(
    tmp_path: Path,
) -> None:
    """The ledger starts as the amount of the initial concentrations, zone by zone and layer by layer."""
    n = 4
    settings = _settings(n)
    builder = ZonationBuilder(settings, [DRG], None)
    model = builder.build()
    expected = 0.0
    cells = [sp for sp in model.species if "_z" in str(sp.sid)]
    for sp in cells:
        sid = str(sp.sid)
        zone = int(sid.rsplit("_z", 1)[1])
        sp.initialConcentration = 0.1 * zone + 0.01 * len(sid)
    add_observers(model, observer_objects(settings, [DRG]))
    path = tmp_path / "initial.xml"
    create_model(model, filepath=path, sbml_level=3, sbml_version=2, validate=False)
    sizes = {
        str(c.sid): 0.01 * (i + 1)
        for i, c in enumerate(model.compartments)
        if str(c.sid).startswith("V_") and c.sid != "V_bile_duct"
    }
    r = individual(path, **sizes)  # a volume of its own in every zone and layer
    for sp in cells:
        expected += r[f"[{sp.sid}]"] * r[str(sp.compartment)]
    assert len({r[v] for v in sizes}) == len(sizes)
    assert r[observers.ledger("dose", "drg")] == pytest.approx(
        expected, rel=1e-12, abs=0
    )
    assert r[observers.amount("drg")] == pytest.approx(expected, rel=1e-12, abs=0)


def test_a_chain_without_flow_has_no_outflow_concentration(tmp_path: Path) -> None:
    """With no flow nothing leaves, and the concentration of nothing is not defined."""
    settings = _settings(3, flow=ImposedFlow(0.0))
    r = _runner(tmp_path, settings, sink=0.0)
    assert "c_outflow_drg" not in r.model.getGlobalParameterIds()


# --- the amounts ------------------------------------------------------------


@pytest.mark.parametrize(
    "kwargs",
    [
        {},
        {"haematocrit": 0.0},
        {"bile": None},
        {"geometry": Geometry.WELL_STIRRED},
    ],
    ids=["sinusoid", "hct0", "nobile", "well_stirred"],
)
def test_the_amounts_of_the_layers_sum_to_the_organ(
    tmp_path: Path, kwargs: dict[str, object]
) -> None:
    """The zones of a layer sum to its amount, the layers to the organ, over the layers made."""
    settings = _settings(4, **kwargs)
    r = _runner(tmp_path, settings)
    for layer in _layer_tokens(r):
        for k in range(1, 5):
            r[f"[drg_{layer}_z{k}]"] = 0.3 + 0.1 * k + 0.05 * len(layer)
    layers = _layer_tokens(r)
    by_layer = {
        layer: sum(r[f"[drg_{layer}_z{k}]"] * r[f"V_{layer}_z{k}"] for k in range(1, 5))
        for layer in layers
    }
    for layer, amount in by_layer.items():
        assert r[observers.layer_amount("drg", layer)] == pytest.approx(
            amount, rel=1e-12, abs=0
        )
        assert amount > 0.0
    assert r[observers.amount("drg")] == pytest.approx(
        sum(by_layer.values()), rel=1e-12, abs=0
    )
    # the layer set is the one of the model, so no layer is silently missing
    reported = {
        sid.removeprefix("amount_drg_")
        for sid in r.model.getGlobalParameterIds()
        if sid.startswith("amount_drg_")
        and sid.removeprefix("amount_drg_") in {*(str(x) for x in Layer), ids.VASCULAR}
    }
    assert reported == set(layers)
    assert ("rbc" in layers) == (
        kwargs.get("haematocrit", 0.45) != 0.0 and "geometry" not in kwargs
    )
    assert ("vas" in layers) == (kwargs.get("geometry") is Geometry.WELL_STIRRED)
    assert ("bile" in layers) == (kwargs.get("bile", 1) is not None)


# --- the mass balance -------------------------------------------------------

_CONFIGS = {
    "advective": {"bile": BileTransport.ADVECTIVE},
    "diffusive": {"bile": BileTransport.DIFFUSIVE},
    "nobile": {"bile": None},
    "hct0": {"haematocrit": 0.0},
    "well_stirred": {"geometry": Geometry.WELL_STIRRED},
    "imposed": {"flow": ImposedFlow(2.0)},
    "noflow": {"flow": ImposedFlow(0.0)},
    "pressure": {"flow": PressureFlow()},
    "pressure_asymmetric": {"flow": PressureFlow(p_onc=2000.0)},
    "pressure_well_stirred": {
        "flow": PressureFlow(),
        "geometry": Geometry.WELL_STIRRED,
    },
}


def _balance_run(
    tmp_path: Path, name: str, **model: object
) -> tuple[roadrunner.RoadRunner, float]:
    """A loaded chain with an inflow, run for a while; the dose and the final state."""
    settings = _settings(5, **_CONFIGS[name])
    r = _runner(tmp_path, settings, **model)
    dose = _load(r, settings, 0.5)
    for layer in ("plasma", "rbc", "vas"):
        if f"drg_{layer}_in" in r.model.getGlobalParameterIds():
            r[f"drg_{layer}_in"] = 1.0
    r.simulate(0, 3, 300)
    return r, dose


@pytest.mark.parametrize("name", list(_CONFIGS))
def test_the_balance_residual_stays_at_zero(tmp_path: Path, name: str) -> None:
    """The residual of the conservation law is zero to integrator tolerance."""
    r, dose = _balance_run(tmp_path, name)
    consumed = r[observers.ledger("consumed", "drg")]
    scale = dose + r[observers.ledger("in", "drg")]
    assert dose > 0.0
    assert consumed > 0.01 * scale  # the sink matters, so the balance is not vacuous
    assert abs(r["balance_drg"]) <= 1e-11 * scale


@pytest.mark.parametrize("name", ["advective", "diffusive"])
def test_the_balance_includes_what_the_duct_collected(
    tmp_path: Path, name: str
) -> None:
    """The bile configurations move drug into the duct, and the residual still vanishes."""
    r, _ = _balance_run(tmp_path, name)
    assert r["drg_bile_duct"] > 1e-4


def test_the_balance_sees_the_outflow_and_the_inflow(tmp_path: Path) -> None:
    """Both boundary ledgers carry drug in the run which stays at zero."""
    r, _ = _balance_run(tmp_path, "imposed")
    assert r[observers.ledger("in", "drg")] > 0.1
    assert r[observers.ledger("out", "drg")] > 0.1


def test_the_balance_fails_when_the_consumption_is_not_declared(
    tmp_path: Path,
) -> None:
    """A sink the observers do not know about is a residual of its size."""
    declared, dose = _balance_run(tmp_path, "imposed")
    consumed = declared[observers.ledger("consumed", "drg")]
    undeclared, _ = _balance_run(tmp_path, "imposed", declared=False)
    assert abs(undeclared["balance_drg"]) == pytest.approx(consumed, rel=1e-6)
    assert abs(declared["balance_drg"]) < 1e-9 * dose


def test_the_balance_fails_for_a_dose_not_recorded(tmp_path: Path) -> None:
    """A drug added by hand without the ledger is a residual of the dose, and with it zero."""
    settings = _settings(5)
    r = _runner(tmp_path, settings)
    r["drg_plasma_in"] = 1.0
    r.simulate(0, 1, 50)
    added = 0.25
    r["[drg_plasma_z3]"] += added / r["V_plasma_z3"]
    r.simulate(1, 2, 50)
    assert r["balance_drg"] == pytest.approx(added, rel=1e-6)
    r[observers.ledger("dose", "drg")] += added
    assert abs(r["balance_drg"]) < 1e-9


def test_the_balance_fails_for_a_leak(tmp_path: Path) -> None:
    """A reaction which removes drug without being declared is seen."""
    settings = _settings(5)
    builder = ZonationBuilder(settings, [DRG], None)
    model = builder.build()
    model.reactions.append(
        Reaction(
            "R_leak",
            "drg_disse_z2 ->",
            formula=("5 * drg_disse_z2 * V_disse_z2 * 1 per_min", U.mmole_per_min),
            compartment="V_disse_z2",
        )
    )
    add_observers(model, observer_objects(settings, [DRG]))
    path = tmp_path / "leak.xml"
    create_model(model, filepath=path, sbml_level=3, sbml_version=2, validate=False)
    r = individual(path, P_endo_in_drg=PERMEABILITY, P_endo_out_drg=PERMEABILITY)
    r["[drg_disse_z2]"] = 1.0
    r[observers.ledger("dose", "drg")] = r["V_disse_z2"]
    r.simulate(0, 1, 50)
    assert abs(r["balance_drg"]) > 1e-4


# --- the observers read, they do not act ------------------------------------


def test_the_observers_do_not_change_the_trajectories(tmp_path: Path) -> None:
    """The states of a model with and without the observers are the same."""
    settings = _settings(5)
    runs = {}
    for watched in (True, False):
        directory = tmp_path / str(watched)
        directory.mkdir()
        r = _runner(directory, settings, watched=watched)
        _ = [
            r.__setitem__(f"[drg_{layer}_z{k}]", 0.4)
            for layer in _layer_tokens(r)
            for k in range(1, 6)
        ]
        r["drg_plasma_in"] = 1.0
        r["drg_rbc_in"] = 1.0
        names = [
            f"[{sid}]"
            for sid in r.model.getFloatingSpeciesIds()
            if sid.startswith("drg_") and "duct" not in sid
        ]
        runs[watched] = r.simulate(0, 4, 81, selections=["time", *names])
    np.testing.assert_allclose(runs[True], runs[False], rtol=1e-7, atol=1e-10)
    assert np.ptp(runs[False][:, 1:]) > 0.1  # the states actually move


def test_the_observers_are_not_needed_to_build() -> None:
    """`build()` has no observers; they are added to the model on request."""
    model = ZonationBuilder(_settings(3), [DRG], None).build()
    assert not any(p.sid.startswith("balance_") for p in model.parameters)
    added = ZonationBuilder(_settings(3), [DRG], None).build()
    add_observers(added, observer_objects(_settings(3), [DRG]))
    assert any(p.sid == "balance_drg" for p in added.parameters)


# --- the axial profile ------------------------------------------------------


def test_the_profile_ids_run_from_periportal_to_pericentral(tmp_path: Path) -> None:
    """Reading `profile_ids` in order gives the zones periportal first."""
    settings = _settings(6)
    r = _runner(tmp_path, settings)
    listed = profile_ids("drg", Layer.PLASMA, 6)
    assert listed == [f"drg_plasma_z{k}" for k in range(1, 7)]
    for k, sid in enumerate(listed, start=1):
        r[f"[{sid}]"] = float(k)
    assert [r[f"[{sid}]"] for sid in profile_ids("drg", Layer.PLASMA, 6)] == [
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
    ]  # fmt: skip


def test_sync_dose_refuses_a_started_simulation(tmp_path: Path) -> None:
    """After the start the amount present is no longer the dose."""
    r = _runner(tmp_path, _settings(3))
    r.simulate(0, 1, 10)
    with pytest.raises(ValueError, match="only the dose at the start"):
        observers.sync_dose(r, [DRG])


def test_without_sync_dose_the_balance_shows_the_amount_set(tmp_path: Path) -> None:
    """Writing species after load leaves the ledger at 0: the residual is the amount set."""
    settings = _settings(3, flow=ImposedFlow(0.0), bile=None)
    r = _runner(tmp_path, settings, sink=0.0)
    r["[drg_plasma_z2]"] = 1.0
    assert r["balance_drg"] == pytest.approx(r["V_plasma_z2"], rel=1e-12, abs=0)
    observers.sync_dose(r, ["drg"])
    assert abs(r["balance_drg"]) < 1e-15


def test_the_observers_of_a_composed_model_close_the_balance(tmp_path: Path) -> None:
    """After `compose`, flattening and libsbml validation, the cell's uptake is declared and the balance holds."""
    import libsbml

    from pkpdsim.zonation.composition import CellPorts, compose, write

    from .data.toy_cell import TOY_ID, create_toy_cell

    n = 4
    settings = _settings(n)
    create_toy_cell(tmp_path)
    builder = ZonationBuilder(settings, [DRG], None)
    model = compose(builder, TOY_ID, CellPorts())
    uptake = " + ".join(f"M{k:02}__UP" for k in range(1, n + 1))
    add_observers(model, observer_objects(settings, [DRG], {"drg": uptake}))
    flat = write(model, tmp_path)
    document = libsbml.readSBMLFromFile(str(flat))
    document.checkConsistency()
    problems = [
        document.getError(i).getMessage()[:400] for i in range(document.getNumErrors())
    ]
    assert problems == []  # errors and warnings both count
    r = roadrunner.RoadRunner(str(flat))
    r.integrator.absolute_tolerance = 1e-14
    r.integrator.relative_tolerance = 1e-10
    r["drg_plasma_in"] = 1.0
    r["drg_rbc_in"] = 1.0
    for layer in _layer_tokens(r):
        for k in range(1, n + 1):
            r[f"[drg_{layer}_z{k}]"] = 0.5
    observers.sync_dose(r, [DRG])
    dose = r[observers.ledger("dose", "drg")]
    r.simulate(0, 3, 300)
    consumed = r[observers.ledger("consumed", "drg")]
    in_cells = sum(r[f"M{k:02}__drg"] for k in range(1, n + 1))  # amounts, mmole
    assert consumed > 0.01 * dose
    assert consumed == pytest.approx(
        in_cells, rel=1e-6, abs=0
    )  # the toy keeps what it takes up
    scale = dose + r[observers.ledger("in", "drg")]
    assert abs(r["balance_drg"]) <= 1e-11 * scale


@pytest.mark.parametrize(
    "flow",
    [ImposedFlow(1.5), PressureFlow(), PressureFlow(p_onc=2000.0)],
    ids=["imposed", "pressure", "pressure_asymmetric"],
)
def test_the_boundary_ledgers_integrate_the_flows_of_their_zones(
    tmp_path: Path, flow: ImposedFlow | PressureFlow
) -> None:
    """`cum_in` grows with the carrier flows of z1 and `cum_out` with those of the last zone."""
    n = 5
    settings = _settings(n, flow=flow)
    r = _runner(tmp_path, settings)
    c_in = {"plasma": 1.0, "rbc": 2.0}
    for layer, c in c_in.items():
        r[f"drg_{layer}_in"] = c
        for k in range(1, n + 1):
            r[f"[drg_{layer}_z{k}]"] = c * (1.0 + 0.1 * k)

    def carrier_flow(layer: Layer, k: int) -> float:
        sid = (
            ids.layer_flow(layer)
            if isinstance(flow, ImposedFlow)
            else ids.zone_flow(layer, ids.zone(k))
        )
        return float(r[sid])

    expected_in = sum(
        carrier_flow(layer, 1) * c
        for layer, c in ((Layer.PLASMA, 1.0), (Layer.RBC, 2.0))
    )
    expected_out = sum(
        carrier_flow(layer, n) * r[f"[drg_{layer}_z{n}]"]
        for layer in (Layer.PLASMA, Layer.RBC)
    )
    assert r[observers.ledger("in", "drg") + "'"] == pytest.approx(
        expected_in, rel=1e-12, abs=0
    )
    assert r[observers.ledger("out", "drg") + "'"] == pytest.approx(
        expected_out, rel=1e-12, abs=0
    )
    if not isinstance(flow, ImposedFlow):
        # the first and the last zone differ, so a swapped zone is seen
        first, last = carrier_flow(Layer.PLASMA, 1), carrier_flow(Layer.PLASMA, n)
        print(f"Q_plasma z1 {first!r} z{n} {last!r}")
        if flow.p_onc is not None:
            assert first != pytest.approx(last, rel=1e-3, abs=0)
