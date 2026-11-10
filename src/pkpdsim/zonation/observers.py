"""The observers of a zonated organ: the quantities a modeller plots and fits.

The observers are read-only: assignment rules and the rate rules of cumulative
ledgers, which no reaction of the model reads, so a model with them has the same
states as a model without them. They are not part of `ZonationBuilder.build()`;
`observer_objects` creates them and `add_observers` puts them into a model (for a
composed model after `compose`), so a caller who does not want them does not pay
for them.

Per substance `<s>` (layers are the ones whose species the model has: `rbc` only
at a haematocrit above 0, `bile` only with a bile, the single `vas` of the
well-stirred geometry instead of `rbc`, `plasma`, `endo` and `disse` (see
`ZonationBuilder.layers`), and no `cell` with `cell_species=False`, whose cells
are the species of a cell model):

- `<s>_pp`, `<s>_pc` (mM): the periportal and the pericentral mean of the
  concentration in the blood side of the lumen (the plasma, `vas` in the
  well-stirred geometry), and `<s>_<layer>_pp`, `<s>_<layer>_pc` of every other
  layer. The mean is that of the first (periportal) and the last (pericentral)
  third of the sinusoid, every zone weighted by `pkpdsim.zones.third_weights`
  (the weighting every zonated model of pkpdsim shares, not a copy).
- `c_outflow_<s>` (mM): the concentration of the blood leaving the last zone, the
  flux out over the flow out. It is **not** the pericentral mean (`<s>_pc`), which
  averages the last third of the length: with a gradient the two differ. It is not
  defined without a flow and not emitted then.
- `amount_<s>_<layer>` (mmole): the amount in a layer, the sum over the zones, and
  `amount_<s>` the sum over the layers, the amount in the organ (the bile duct is
  outside of the organ).
- `cum_in_<s>`, `cum_out_<s>` (mmole): what entered through the inflow and left
  through the outflow since the start; the outflow ledger holds the filtrate of
  a net filtration as well (`ZonationBuilder.filtrate_efflux`, which leaves the
  space of Disse at the pericentral end and is negative for a net
  reabsorption). The two and `balance_<s>` are of the ends of the chain and are
  written only for the boundary of the engine (`FlowBoundary.OPEN`): with an
  external boundary they belong to the caller, which holds both ends.
  `cum_consumed_<s>` (mmole): what the
  reactions the caller declares consumed. `cum_dose_<s>` (mmole): what was put into
  the organ by hand, the initial amounts and the doses; it starts as the amount of
  the initial concentrations of the model and a dose or a changed initial amount
  adds to it.
- `balance_<s>` (mmole), the mass balance residual, see `balance_expression`.

The axial profile is the species of the zones themselves, `<s>_<layer>_z<k>`;
`profile_ids` lists them periportal first, so no copy of them is made.
"""

from collections.abc import Mapping, Sequence

import roadrunner
from sbmlutils.factory import Model, Parameter, RateRule
from sbmlutils.metadata import SBO

from pkpdsim.sbml.annotation import initial_assignment
from pkpdsim.sbml.objects import SbmlObjects
from pkpdsim.sbml.parameters import RULE, observer
from pkpdsim.sbml.units import U
from pkpdsim.zonation import ids
from pkpdsim.zonation.builder import ZonationBuilder
from pkpdsim.zonation.geometry import Layer
from pkpdsim.zonation.settings import FlowBoundary, ZonationSettings
from pkpdsim.zonation.transport import Solute
from pkpdsim.zones import third_weights, weighted_mean

__all__ = [
    "add_observers",
    "amount",
    "balance",
    "balance_expression",
    "layer_amount",
    "ledger",
    "mean",
    "observer_objects",
    "outflow_concentration",
    "profile_ids",
    "sync_dose",
]

PERIPORTAL = "pp"
"""The suffix of the periportal mean."""
PERICENTRAL = "pc"
"""The suffix of the pericentral mean."""
LEDGERS = ("in", "out", "consumed", "dose")
"""The kinds of the cumulative amounts."""


def mean(
    substance: str, layer: Layer | str, periportal: bool, blood: Layer | str
) -> str:
    """Get the id of the mean of a substance in a layer over the first or last third.

    Args:
        substance: the id of the substance.
        layer: the layer.
        periportal: the first third, else the last.
        blood: the layer of the blood side of the lumen (`plasma`, or `vas`),
            whose means carry no layer in their id.

    Returns:
        `<s>_pp` / `<s>_pc` for the blood side, `<s>_<layer>_pp` / `_pc` else.
    """
    site = PERIPORTAL if periportal else PERICENTRAL
    return f"{substance}_{site}" if layer == blood else f"{substance}_{layer}_{site}"


def outflow_concentration(substance: str) -> str:
    """Get the id of the concentration of the blood leaving the organ, `c_outflow_<s>`."""
    return f"c_outflow_{substance}"


def amount(substance: str) -> str:
    """Get the id of the amount of a substance in the organ, `amount_<s>`."""
    return f"amount_{substance}"


def layer_amount(substance: str, layer: Layer | str) -> str:
    """Get the id of the amount of a substance in a layer, `amount_<s>_<layer>`."""
    return f"amount_{substance}_{layer}"


def ledger(kind: str, substance: str) -> str:
    """Get the id of a cumulative amount, `cum_<kind>_<s>` (kinds: `LEDGERS`).

    Raises:
        ValueError: for a kind which is no ledger.
    """
    if kind not in LEDGERS:
        raise ValueError(f"A ledger is one of {LEDGERS}, not '{kind}'.")
    return f"cum_{kind}_{substance}"


def balance(substance: str) -> str:
    """Get the id of the mass balance residual of a substance, `balance_<s>`."""
    return f"balance_{substance}"


def profile_ids(substance: str, layer: Layer | str, n_zones: int) -> list[str]:
    """Get the species of the axial profile of a substance in a layer.

    The profile is an array observable: one species per zone, concentrations
    in mM, which a caller reads in this order. The order is **periportal first**:
    index 0 is zone `z1` and the last is the pericentral `z<n>`, the order of
    `third_weights` and of `zone_positions` (0 periportal, 1 pericentral). The
    species are read as `[<id>]` in roadrunner. The choice of the species over
    a parameter or an assignment per zone is that no copy of a state has to
    be kept equal to it.

    Args:
        substance: the id of the substance.
        layer: the layer.
        n_zones: the number of zones.

    Returns:
        The ids from periportal to pericentral.
    """
    return [ids.species(substance, layer, ids.zone(k)) for k in range(1, n_zones + 1)]


def balance_expression(substance: str, duct: bool) -> str:
    """Get the formula of the mass balance residual of a substance.

    The conservation law is: what is in the organ now, plus what left through
    the outflow and into the bile duct and what the reactions consumed, equals
    what entered through the inflow plus what was put in by hand:

        balance = amount + cum_out + duct + cum_consumed - cum_in - cum_dose

    It stays at 0 (to the tolerance of the integrator) whenever the model
    moves the substance only between its compartments, the boundaries and the
    reactions declared as consumption: the axial convection and diffusion, the
    radial exchange and the bile in both transports, with and without a bile,
    at a haematocrit of 0, in the well-stirred geometry and with an imposed or a
    pressure driven flow. The solute is moved only by the reactions, each of
    which takes from one species what it gives to another, and the boundary
    ledgers integrate the same flows
    the first and last reaction use, so the balance holds there as well. It
    fails for a consumption or a dose which is not declared.

    The residual is written for the boundary of the engine
    (`FlowBoundary.OPEN`); a model whose ends a caller supplies has no inflow
    and no outflow ledger of its own and therefore no residual of its own.

    A caller who sets species after loading the model calls `sync_dose` once
    before the simulation: the initial assignment of `cum_dose_<s>` is
    evaluated when the model is loaded and reset, not when a caller writes a
    species (the behaviour of roadrunner), so without the call the residual
    equals the amount set by hand.

    The observers are read-only, but the ledgers are states the integrator
    controls the error of, so a model with observers and one without agree to
    the tolerance of the integrator and not bit for bit.

    Args:
        substance: the id of the substance.
        duct: whether the model has a bile duct.
    """
    terms = [
        amount(substance),
        ledger("out", substance),
        ledger("consumed", substance),
    ]
    if duct:
        terms.append(ids.species(substance, Layer.BILE, ids.DUCT))
    negative = [ledger("in", substance), ledger("dose", substance)]
    return " + ".join(terms) + " - " + " - ".join(negative)


def _layer_name(layer: Layer | str) -> str:
    """Get the name of a layer for the names of the observers."""
    return "vascular space" if layer == ids.VASCULAR else str(layer)


def _ledger(
    out: SbmlObjects,
    kind: str,
    substance: Solute,
    rate: str | None,
    start: str | None = None,
) -> None:
    """Add a cumulative amount, integrated from a rate (mmole/min), or constant 0."""
    sid = ledger(kind, substance.sid)
    names = {
        "in": "entered through the inflow",
        "out": "left through the outflow",
        "consumed": "consumed by the cells",
        "dose": "put into the organ by hand (initial amounts and doses)",
    }
    name = f"amount of {substance.name} {names[kind]}"
    out.parameters.append(
        Parameter(
            sid,
            0.0,
            U.mmole,
            constant=rate is None,
            name=name,
            sboTerm=SBO.AMOUNT_OF_AN_ENTITY_POOL,
        )
    )
    if rate is not None:
        out.rate_rules.append(
            RateRule(sid, rate, U.mmole_per_min, name=name, sboTerm=RULE)
        )
    if start is not None:
        out.assignments.append(initial_assignment(sid, start, U.mmole, name))


def observer_objects(
    settings: ZonationSettings,
    substances: Sequence[Solute],
    consumption: Mapping[str, str] | None = None,
) -> SbmlObjects:
    """Get the observers of a zonated organ.

    The weighting of the thirds is `pkpdsim.zones.third_weights`, which every
    zonated model of pkpdsim shares. `mean_objects` of
    `pkpdsim.perfusion.observers` is not reused: it writes the species ids of the
    perfusion core (`<q>_<site>`, the cells, the living fraction, the pH and the
    activities of the sites), which this model does not have; its weights are
    the part which carries over.

    Args:
        settings: the settings of the model.
        substances: the substances of the model.
        consumption: per substance the formula (mmole/min) of the net rate at
            which the reactions of the cells consume it, summed over the zones,
            e.g. the ids of the reactions of the cell models; a substance
            without one is not consumed.

    Returns:
        The parameters, the assignment rules, the rate rules and the initial
        assignments of the observers of every substance.

    Note:
        The layers are the ones whose species the model has: a model with
        `cell_species=False` keeps the species of its cells in its cell model,
        under the ids of that model, so the engine has no amount and no mean
        of the cell layer to write and writes none.

        A model with an external boundary (`FlowBoundary.EXTERNAL`, a
        perfusion circuit) has no inflow and no outflow of its own: the
        ledgers `cum_in_<s>` and `cum_out_<s>` and the residual `balance_<s>`
        belong to the caller, which holds the two ends of the chain, and are
        not written. `cum_consumed_<s>`, `cum_dose_<s>`, the amounts, the
        means and `c_outflow_<s>` (the concentration of the blood leaving the
        last zone, which needs no boundary) are written as always.
    """
    consumption = consumption or {}
    builder = ZonationBuilder(settings, substances, None)
    n = settings.n_zones
    layers = [
        layer
        for layer in builder.layers()
        if layer is not Layer.CELL or settings.cell_species
    ]
    closed = settings.boundary is FlowBoundary.OPEN
    carriers = builder.carriers()
    blood = ids.VASCULAR if ids.VASCULAR in layers else Layer.PLASMA
    last = ids.zone(n)
    flowing = builder.has_flow()
    out = SbmlObjects()
    for substance in substances:
        s = substance.sid
        for layer in layers:
            of_zones = profile_ids(s, layer, n)
            for periportal in (True, False):
                where = "periportal" if periportal else "pericentral"
                p, r = observer(
                    mean(s, layer, periportal, blood),
                    weighted_mean(of_zones, third_weights(n, periportal)),
                    U.mM,
                    f"{substance.name} ({_layer_name(layer)}, {where} mean)",
                    SBO.MOLAR_CONCENTRATION_OF_AN_ENTITY,
                )
                out.parameters.append(p)
                out.rules.append(r)
        for layer in layers:
            p, r = observer(
                layer_amount(s, layer),
                " + ".join(
                    f"{ids.species(s, layer, ids.zone(k))} * {ids.volume(layer, ids.zone(k))}"
                    for k in range(1, n + 1)
                ),
                U.mmole,
                f"amount of {substance.name} in the {_layer_name(layer)}",
                SBO.AMOUNT_OF_AN_ENTITY_POOL,
            )
            out.parameters.append(p)
            out.rules.append(r)
        p, r = observer(
            amount(s),
            " + ".join(layer_amount(s, layer) for layer in layers),
            U.mmole,
            f"amount of {substance.name} in the organ",
            SBO.AMOUNT_OF_AN_ENTITY_POOL,
        )
        out.parameters.append(p)
        out.rules.append(r)

        influx = efflux = None
        if flowing:
            # the flux of the blood out of the last zone, and the flow which
            # carries it: `blood` is the layer of the blood and stays it
            flux_out = " + ".join(
                f"{builder.carrier_flow(layer, last)} * {ids.species(s, layer, last)}"
                for layer in carriers
            )
            flow_out = " + ".join(
                builder.carrier_flow(layer, last) for layer in carriers
            )
            p, r = observer(
                outflow_concentration(s),
                f"({flux_out}) / ({flow_out})",
                U.mM,
                f"{substance.name} in the blood leaving the organ (not the pericentral mean)",
                SBO.MOLAR_CONCENTRATION_OF_AN_ENTITY,
            )
            out.parameters.append(p)
            out.rules.append(r)
        if flowing and closed:
            # the inflow crosses the periportal face of the chain, where a
            # calibrated pressure flow carries the reference flow; the flow of
            # the first zone is that of its pericentral face
            influx = " + ".join(
                f"{builder.carrier_flow(layer, ids.INFLOW)} * "
                f"{ids.species(s, layer, ids.INFLOW)}"
                for layer in carriers
            )
            # the filtrate which a net filtration takes out of the space of
            # Disse at the pericentral end leaves the organ as well, so the
            # outflow ledger integrates it (negative for a net reabsorption,
            # which draws water in there); it is not part of `c_outflow`,
            # which is the concentration of the blood
            lymph = builder.filtrate_efflux(substance)
            efflux = flux_out if lymph is None else f"{flux_out} + {lymph}"
        if closed:
            _ledger(out, "in", substance, influx)
            _ledger(out, "out", substance, efflux)
        _ledger(out, "consumed", substance, consumption.get(s))
        _ledger(
            out,
            "dose",
            substance,
            None,
            start=" + ".join(
                f"{ids.species(s, layer, ids.zone(k))} * {ids.volume(layer, ids.zone(k))}"
                for layer in layers
                for k in range(1, n + 1)
            ),
        )
        if closed:
            p, r = observer(
                balance(s),
                balance_expression(s, settings.bile is not None),
                U.mmole,
                f"mass balance residual of {substance.name} (0 if conserved)",
                SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            )
            out.parameters.append(p)
            out.rules.append(r)
    return out


def add_observers(model: Model, objects: SbmlObjects) -> None:
    """Add the observers to a model, which they do not change.

    Args:
        model: the model of `ZonationBuilder.build()` or of `compose`.
        objects: the observers (`observer_objects`).
    """
    model.parameters.extend(objects.parameters)
    model.rules.extend(objects.rules)
    model.rate_rules.extend(objects.rate_rules)
    model.assignments.extend(objects.assignments)


def sync_dose(
    runner: roadrunner.RoadRunner, substances: Sequence[Solute | str]
) -> None:
    """Bring the dose ledger in line with the amounts present at the start.

    Sets `cum_dose_<s>` to `amount_<s>`, the amount of the substance in the organ,
    for a model whose species a caller wrote after loading it. Call it once,
    after setting the species and before the first simulation; a dose given
    during a run is added to `cum_dose_<s>` by the caller. Amounts which a
    cell model holds in its own species are not in `amount_<s>`; the cell
    model's uptake into them is declared as consumption.

    Args:
        runner: the loaded roadrunner model with the observers.
        substances: the substances (or their ids).

    Raises:
        ValueError: if the simulation has started (time above 0), when the
            amount present is no longer the dose.
    """
    if runner["time"] != 0.0:
        raise ValueError(
            "sync_dose sets the dose to the amount present, which is only the dose "
            f"at the start; the model is at time {runner['time']}."
        )
    for substance in substances:
        s = substance.sid if isinstance(substance, Solute) else substance
        runner[ledger("dose", s)] = runner[amount(s)]
