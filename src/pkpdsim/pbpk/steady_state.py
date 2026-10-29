"""The steady state of the equations of the library, as formulas.

An endogenous substance starts at the equilibrium of every part with its
baseline (`Substance.baseline`), which is stationary without production and
elimination. An organ with a production or an elimination of the chemistry
(a custom reaction) is at another steady state: the plasma leaving it carries
the production, its tissue holds more than `Kp` times its plasma, the kidney
plasma lacks the secretion. A model gives such an organ its steady state as
the starts of its species (`PBPKDefinition.starts`); these functions return
the formulas of the steady state of the equations the library builds, over
ids of the top model, so that the model does not restate them:

- `plasma_flow`: the flow which carries the substance through an organ,
  `Q<o>*(1 - HCT)`, `Q<o>*BP` with erythrocytes in equilibrium;
- `tissue_volume`: the volume of the merged tissue of an organ, `Vti`,
  where a custom reaction states its rate per litre of tissue;
- `tissue_exchange`: the exchange of the merged tissue with the plasma,
  `PS*fu`, the share of `ftissue` of the organ;
- `organ_plasma`, `organ_tissue`: the plasma leaving a standard organ and its
  merged tissue with a net source `rate` (mmole/min) and a clearance of the
  plasma;
- `renal_filtration`, `kidney_plasma`, `kidney_arterial`, `renal_clearance`:
  the net filtrate of the kidney, the kidney plasma with a secretion (a
  clearance of the kidney plasma and a saturable secretion, `Saturable`, the
  root of a quadratic) at a given arterial plasma or a given renal
  elimination, the arterial plasma of a kidney plasma, and the renal
  clearance of the arterial plasma (glomerulus and tubule in series); a
  saturable reabsorption of the filtrate (`Excretion.reabsorption`, the
  `pkpdsim.pbpk.Saturable` of the substance, not this `Saturable` of a
  secretion) has the net filtrate of the arterial plasma and, without
  another process of the kidney, the kidney plasma at the arterial plasma;
- `liver_plasma`: the liver plasma, fed by the hepatic artery and the portal
  vein past the shunt `f_shunt`, with a source of the liver;
- `arm_vein`: the mix of the venous outflows of the arm vein of a preset.

Every formula is the steady state of an organ at a given inflow (the
arterial plasma, the baseline by default): the standard organ is linear in
its species at a fixed inflow, so it is a closed form. The formulas read the
same ids as the equations of the organs (the flows, `HCT`, `BP_<s>`, `fu_<s>`,
`GFR`, `f_reabsorbed_<s>` (`reab_threshold_<s>`, `reab_n_<s>`), the coefficients `Kp`, `ftissue_<s>`,
`Vti_total_ref_<s>`, the organ volume and its vascular fraction, `f_shunt`,
`f_functional`), so they follow a change of the individual before the
simulation, and the tests of the library check them against the equations of
every organ variant they apply to. Where no closed form applies, a function
raises `NotImplementedError` naming the reason: kinetic erythrocytes, the
full form and the flow-limited organ, an exchange of the organ of its own
(`ps`), the intestine (it refuses custom reactions, its parts start at the
equilibrium), the nephron kidney (its tubule starts empty, not at its steady
state) and a source in the zonated liver (it refuses custom reactions). The
functions of an organ which has variants take the variants of the definition
(`variants=`, required), so a model with a variant never gets the formula of
the standard organ silently. Every function which takes an inflow names it
`arterial` (and `portal` for the liver).
"""

from collections.abc import Mapping
from dataclasses import dataclass

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.definition import Variants
from pkpdsim.pbpk.observers import arm_vein_weights
from pkpdsim.pbpk.organs import Organ, Pool
from pkpdsim.pbpk.placement import ONE, kp_id
from pkpdsim.pbpk.substance import Kp, KpMethod, Mode, Substance

__all__ = [
    "Saturable",
    "arm_vein",
    "kidney_arterial",
    "kidney_plasma",
    "liver_plasma",
    "organ_plasma",
    "organ_tissue",
    "plasma_flow",
    "renal_clearance",
    "renal_filtration",
    "tissue_exchange",
    "tissue_volume",
]

_GUT = frozenset({Organ.SMALL_INTESTINE, Organ.LARGE_INTESTINE})
"""The organs of the intestine."""


def _variant(variants: Variants, organ: Organ) -> object | None:
    """The variant of an organ, `None` for the standard organ."""
    chosen = variants.organs.get(organ)
    name = getattr(chosen, "name", chosen)
    return None if name in (None, "standard") else chosen


def _carrier(substance: Substance) -> str:
    """The share of the blood flow which carries the plasma concentration."""
    ery = substance.erythrocytes
    if ery is None:
        return f"({ONE} - {ids.HCT})"
    if ery.mode == "equilibrium":
        return ids.bp(substance.sid)
    raise NotImplementedError(
        f"{substance.sid}: no closed steady state with kinetic erythrocytes "
        "(the blood cells exchange with the plasma of every organ)."
    )


def _flow_of(substance: Substance, flow: str) -> str:
    """A blood flow times the share which carries the plasma concentration."""
    return f"{flow} * {_carrier(substance)}"


def _terms(*terms: str) -> str:
    """The sum of the terms which are not `0`."""
    kept = [t for t in terms if t != "0"]
    return " + ".join(f"({t})" for t in kept) if kept else "0"


def _check_gut(organ: Organ) -> None:
    """Refuse the organs of the intestine."""
    if organ in _GUT:
        raise NotImplementedError(
            f"No steady state of the organ '{organ.name.lower()}': the intestine refuses "
            "custom reactions, its parts start at the equilibrium with the baseline."
        )


def _check_variant(organ: Organ, variants: Variants) -> None:
    """Refuse an organ with a variant."""
    if _variant(variants, organ) is not None:
        raise NotImplementedError(
            f"No steady state of the variant of the organ '{organ.name.lower()}'."
        )


def plasma_flow(substance: Substance, organ: Organ, *, variants: Variants) -> str:
    """Get the flow which carries the substance through an organ [l/min].

    `Q<o>*(1 - HCT)`, the plasma flow, or `Q<o>*BP_<s>` with erythrocytes in
    equilibrium (the plasma species carries the blood); the liver is perfused
    by `(1 - f_shunt)*Qh`. The flows of the variants of the liver and the
    kidney are those of the standard organ.

    Args:
        substance: the substance.
        organ: the organ, not the intestine.
        variants: the variants of the definition (`PBPKDefinition.variants`),
            required.

    Returns:
        The formula.

    Raises:
        NotImplementedError: with kinetic erythrocytes, for the intestine.
    """
    _check_gut(organ)
    if organ is Organ.LIVER:
        return f"({ONE} - {ids.F_SHUNT}) * " + _flow_of(substance, ids.QH)
    return _flow_of(substance, ids.organ_flow(organ))


def tissue_volume(organ: Organ, *, variants: Variants) -> str:
    """Get the volume of the merged tissue of a standard organ, `Vti` [l].

    `V<o>*(1 - f_vas_<o>)`, of the liver times `f_functional`: the local
    compartment `Vti` of the organ as a formula of the top model, e.g. for the
    rate of a custom reaction per litre of tissue at the steady state.

    Args:
        organ: the organ, a standard organ, not the intestine.
        variants: the variants of the definition (`PBPKDefinition.variants`),
            required: a variant is refused.

    Returns:
        The formula.

    Raises:
        NotImplementedError: for the intestine or a variant of the organ.
    """
    _check_gut(organ)
    _check_variant(organ, variants)
    tissue = f"{ids.organ_volume(organ)} * ({ONE} - {ids.vascular_fraction(organ)})"
    if organ is Organ.LIVER:
        tissue = f"{tissue} * {ids.F_FUNCTIONAL}"
    return tissue


def _check_merged(substance: Substance, organ: Organ) -> Kp | KpMethod:
    """Get the distribution of a substance which enters the merged tissue of an organ.

    Raises:
        ValueError: if the substance does not enter the tissue of the organ.
        NotImplementedError: for another mode than `Mode.TISSUE`.
    """
    dist = substance.distribution
    if dist is None or kp_id(substance, organ) is None:
        raise ValueError(
            f"{substance.sid} does not enter the tissue of the organ '{organ.name.lower()}'."
        )
    if dist.mode is not Mode.TISSUE:
        raise NotImplementedError(
            f"{substance.sid}: no closed steady state of the tissue in mode "
            f"'{dist.mode.value}', only of the merged tissue (Mode.TISSUE)."
        )
    return dist


def tissue_exchange(substance: Substance, organ: Organ, *, variants: Variants) -> str:
    """Get the exchange of the merged tissue of an organ with its plasma, `PS*fu` [l/min].

    The share of the whole-body exchange `ftissue_<s>` by the tissue volume of
    the organ (of the liver times `f_functional`) relative to the reference
    tissue volume, times the fraction unbound: the exchange
    `PS*fu*(Cpl - Cti/Kp)` of the merged tissue (`Mode.TISSUE`).

    Args:
        substance: the substance.
        organ: the organ, a standard organ, not the intestine.
        variants: the variants of the definition (`PBPKDefinition.variants`),
            required: a variant is refused.

    Returns:
        The formula.

    Raises:
        ValueError: if the substance does not enter the tissue of the organ.
        NotImplementedError: for another mode than `Mode.TISSUE`, an
            exchange of the organ of its own (`ps`), the intestine or a variant
            of the organ.
    """
    dist = _check_merged(substance, organ)
    if organ in dist.ps:
        raise NotImplementedError(
            f"{substance.sid}: the exchange of the organ '{organ.name.lower()}' is a "
            "value of the organ (ps), no constant of the top model."
        )
    return (
        f"{ids.ftissue(substance.sid)} * {tissue_volume(organ, variants=variants)} / "
        f"{ids.tissue_volume_total_ref(substance.sid)} * {ids.fu(substance.sid)}"
    )


def _check_standard(organ: Organ, variants: Variants) -> None:
    """Refuse the organs whose steady state the plasma of a standard organ is not."""
    _check_gut(organ)
    if organ is Organ.KIDNEY:
        raise ValueError("The kidney filters: use `kidney_plasma`.")
    if organ is Organ.LIVER:
        raise ValueError("The liver has two inflows and a shunt: use `liver_plasma`.")
    _check_variant(organ, variants)


def organ_plasma(
    substance: Substance,
    organ: Organ,
    *,
    arterial: str | None = None,
    rate: str = "0",
    clearance: str = "0",
    variants: Variants,
) -> str:
    """Get the plasma leaving a standard organ at its steady state [mM].

    The organ takes `F*Cin` (`plasma_flow`) and the net source `rate` of its
    chemistry (mmole/min, into the plasma or the tissue: at the steady state
    the tissue passes it on), and loses `F*Cpl` and `clearance*Cpl` (a first
    order elimination of its plasma): `Cpl = (F*Cin + rate)/(F + clearance)`.

    Args:
        substance: the substance.
        organ: the organ, not the kidney or the liver (`kidney_plasma`,
            `liver_plasma`) or the intestine.
        arterial: the arterial plasma flowing in, the baseline
            `baseline_<s>` by default.
        rate: the net source of the organ at the steady state (mmole/min),
            negative for an uptake.
        clearance: the clearance of the plasma of the organ (l/min).
        variants: the variants of the definition (`PBPKDefinition.variants`),
            required: a variant without a closed steady state is refused.

    Returns:
        The formula.
    """
    _check_standard(organ, variants)
    flow = plasma_flow(substance, organ, variants=variants)
    c_in = ids.baseline(substance.sid) if arterial is None else arterial
    return f"({_terms(f'{flow} * ({c_in})', rate)}) / ({_terms(flow, clearance)})"


def organ_tissue(
    substance: Substance,
    organ: Organ,
    plasma: str,
    *,
    rate: str = "0",
    variants: Variants,
) -> str:
    """Get the merged tissue of an organ at its steady state [mM].

    With the source `rate` of the tissue (mmole/min) the exchange
    `PS*fu*(Cpl - Cti/Kp)` passes it to the plasma:
    `Cti = Kp*(Cpl + rate/(PS*fu))`, `Kp*Cpl` without.

    Args:
        substance: the substance.
        organ: the organ.
        plasma: the plasma of the organ at the steady state (`organ_plasma`,
            `kidney_plasma`, `liver_plasma`).
        rate: the net source of the tissue (mmole/min).
        variants: the variants of the definition (`PBPKDefinition.variants`),
            required; a tissue of the zonated liver (every zone) without a
            source.

    Returns:
        The formula.
    """
    _check_gut(organ)
    variant = _variant(variants, organ)
    if variant is not None and not (organ is Organ.LIVER and rate == "0"):
        raise NotImplementedError(
            f"No steady state of a source in the variant of the organ '{organ.name.lower()}'."
        )
    _check_merged(substance, organ)
    kp = kp_id(substance, organ)
    if rate == "0":
        return f"{kp} * ({plasma})"
    exchange = tissue_exchange(substance, organ, variants=variants)
    return f"{kp} * ({plasma} + ({rate}) / ({exchange}))"


def renal_filtration(
    substance: Substance, *, arterial: str | None = None, variants: Variants
) -> str:
    """Get the net filtered plasma flow of the kidney, `fu*GFR*(1 - f_reabsorbed)` [l/min].

    With a saturable reabsorption (`Excretion.reabsorption`) the reabsorbed
    fraction is that of the arterial plasma `Car`,
    `fu*GFR*(1 - 1/(1 + (fu*Car/reab_threshold)^reab_n)^(1/reab_n))`.

    Args:
        substance: the substance, with urinary excretion.
        arterial: the arterial plasma of a saturable reabsorption, the
            baseline by default; a constant fraction does not read it.
        variants: the variants of the definition (`PBPKDefinition.variants`),
            required: the nephron kidney is refused.

    Returns:
        The formula.

    Raises:
        ValueError: for a substance without urinary excretion.
        NotImplementedError: for the nephron kidney.
    """
    _check_kidney(variants)
    s = substance.sid
    if not substance.excretion.urine:
        raise ValueError(f"{s} is not filtered (Excretion(urine=False)).")
    if substance.excretion.reabsorption is None:
        return f"{ids.fu(s)} * {ids.GFR} * ({ONE} - {ids.fraction_reabsorbed(s)})"
    car = ids.baseline(s) if arterial is None else arterial
    return f"{ids.fu(s)} * {ids.GFR} * ({ONE} - {_saturable_fraction(substance, car)})"


def _saturable_fraction(substance: Substance, plasma: str) -> str:
    """The reabsorbed fraction of a saturable reabsorption at a plasma concentration."""
    s = substance.sid
    n = ids.reabsorption_n(s)
    threshold = ids.reabsorption_threshold(s)
    return (
        f"{ONE} / ({ONE} + ({ids.fu(s)} * ({plasma}) / {threshold})^{n})^({ONE} / {n})"
    )


def _check_saturable(substance: Substance, closed: bool) -> None:
    """Refuse a saturable reabsorption where the steady state has no closed form.

    The net filtrate of a saturable reabsorption depends on the concentration
    (`Car` in the filtration, the kidney plasma in the outflow): without
    another process the kidney plasma is the arterial plasma, with one the
    balance of the kidney plasma is not linear in it.

    Args:
        substance: the substance.
        closed: whether the steady state asked for is the one without
            another process of the kidney, at a given arterial plasma.

    Raises:
        NotImplementedError: for a saturable reabsorption which is not closed.
    """
    if substance.excretion.reabsorption is not None and not closed:
        raise NotImplementedError(
            f"{substance.sid}: no closed steady state of a saturable reabsorption "
            "with another process of the kidney (a rate, a clearance, a "
            "secretion) or at a given elimination: its net filtrate depends on "
            "the concentration, the balance of the kidney plasma is not linear."
        )


def _check_kidney(variants: Variants) -> None:
    """Refuse the nephron kidney."""
    if _variant(variants, Organ.KIDNEY) is not None:
        raise NotImplementedError(
            "No steady state of the nephron kidney: its tubule starts empty, not "
            "at its steady state."
        )


@dataclass(frozen=True)
class Saturable:
    """A saturable secretion of the kidney plasma at the steady state.

    The rate `vmax*C/(km + C)` of the kidney plasma `C` (mmole/min), a
    Michaelis-Menten secretion of a custom reaction of the kidney into the
    urine, as formulas over ids of the top model.

    Attributes:
        vmax: the capacity of the whole kidney (mmole/min).
        km: the Michaelis constant (mM).
    """

    vmax: str
    km: str


def _kidney_flows(substance: Substance, variants: Variants) -> tuple[str, str, str]:
    """The plasma flow `F`, the net filtrate `G` and the flow past it `R = F - G`."""
    _check_kidney(variants)
    flow = plasma_flow(substance, Organ.KIDNEY, variants=variants)
    if not substance.excretion.urine:
        return flow, "0", flow
    gfr = renal_filtration(substance, variants=variants)
    return flow, gfr, f"{flow} - {gfr}"


def _saturable_root(a: str, vmax: str, km: str, d: str) -> str:
    """The positive root `C` of `a*C + vmax*C/(km + C) = d`.

    `a*C^2 + b*C - d*km = 0` with `b = a*km + vmax - d`, in the form
    `2*d*km/(b + sqrt(b^2 + 4*a*d*km))`, which has no cancellation for
    `b >= 0` (the linear elimination and the capacity at least the load) and
    no division by `a`.
    """
    b = f"({a}) * ({km}) + ({vmax}) - ({d})"
    root = f"sqrt(({b})^2 + 4 dimensionless * ({a}) * ({d}) * ({km}))"
    return f"2 dimensionless * ({d}) * ({km}) / ({b} + {root})"


def kidney_plasma(
    substance: Substance,
    *,
    arterial: str | None = None,
    eliminated: str | None = None,
    rate: str = "0",
    clearance: str = "0",
    secretion: Saturable | None = None,
    variants: Variants,
) -> str:
    """Get the plasma of the kidney at its steady state [mM].

    The glomerulus filters the arterial plasma, the plasma leaving the
    kidney lacks the net filtrate `G` (`renal_filtration`, 0 without urinary
    excretion), a secretion clears the kidney plasma (`clearance`, l/min, and
    a saturable `secretion`, `S(C) = vmax*C/(km + C)`) and a source returns
    `rate` into it (a reabsorption): with `R = F - G`

    `R*(Car - C) + rate = clearance*C + S(C)`.

    With a saturable reabsorption (`Excretion.reabsorption`) the net
    filtrate depends on the concentration: without another process of the
    kidney its plasma is the arterial plasma (the filtration of `Car` and the
    outflow lacking the net filtrate of `C` balance at `C = Car`); a rate, a
    clearance, a secretion or a given elimination have no closed form and
    are refused.

    The kidney plasma `C` follows from the arterial plasma `Car`, or from
    the amount the kidney eliminates into the urine (`eliminated`,
    mmole/min): `E = G*Car + clearance*C + S(C) - rate`, the production of an
    endogenous substance which only the kidney excretes, whose arterial plasma
    is unknown (`kidney_arterial` gives it). Without `secretion` the formula
    is linear in `C`; with it `C` is the positive root of a quadratic.

    Args:
        substance: the substance.
        arterial: the arterial plasma, the baseline by default.
        eliminated: the net excretion of the kidney into the urine at the
            steady state (mmole/min), instead of `arterial`.
        rate: the net source of the kidney (mmole/min).
        clearance: the first order clearance of the kidney plasma (l/min).
        secretion: a saturable secretion of the kidney plasma.
        variants: the variants of the definition (`PBPKDefinition.variants`),
            required: a variant without a closed steady state is refused.

    Returns:
        The formula.

    Raises:
        ValueError: for both `arterial` and `eliminated`, or `eliminated`
            without an elimination of the kidney.
        NotImplementedError: for a saturable reabsorption with a rate, a
            clearance, a secretion or `eliminated`.
    """
    _check_saturable(
        substance,
        eliminated is None and rate == clearance == "0" and secretion is None,
    )
    flow, gfr, rest = _kidney_flows(substance, variants)
    if substance.excretion.reabsorption is not None:
        # the outflow lacks the net filtrate of the kidney plasma: it is the
        # arterial plasma
        return ids.baseline(substance.sid) if arterial is None else arterial
    if eliminated is None:
        car = ids.baseline(substance.sid) if arterial is None else arterial
        load = _terms(f"({rest}) * ({car})", rate)
        linear = _terms(rest, clearance)
        vmax = None if secretion is None else secretion.vmax
    else:
        if arterial is not None:
            raise ValueError("Give the arterial plasma or the elimination, not both.")
        if gfr == "0" and clearance == "0" and secretion is None:
            raise ValueError(f"{substance.sid}: the kidney eliminates nothing.")
        # the kidney plasma passes the arterial plasma on: Car = C + (CL*C +
        # S(C) - rate)/R, so E = G*C + F/R*(CL*C + S(C) - rate)
        ratio = f"({flow}) / ({rest})"
        load = _terms(eliminated, f"{ratio} * ({rate})" if rate != "0" else "0")
        linear = _terms(gfr, f"{ratio} * ({clearance})" if clearance != "0" else "0")
        vmax = None if secretion is None else f"{ratio} * ({secretion.vmax})"
    if secretion is None or vmax is None:
        return f"({load}) / ({linear})"
    if linear == "0":
        # the saturable secretion alone: S(C) = load
        return f"({load}) * ({secretion.km}) / ({vmax} - ({load}))"
    return _saturable_root(linear, vmax, secretion.km, load)


def kidney_arterial(
    substance: Substance,
    plasma: str,
    *,
    eliminated: str | None = None,
    rate: str = "0",
    clearance: str = "0",
    secretion: Saturable | None = None,
    variants: Variants,
) -> str:
    """Get the arterial plasma whose kidney plasma at the steady state is `plasma` [mM].

    The inverse of `kidney_plasma`: `Car = C + (clearance*C + S(C) - rate)/R`
    with the flow past the filtrate `R = F - G`. With the net excretion of
    the kidney into the urine (`eliminated`, the `E` of `kidney_plasma`) it is
    `Car = (E + R*C)/F` whatever the secretion, which reads the kidney plasma
    once: the baseline of an endogenous substance which only the kidney
    excretes, from the kidney plasma of its elimination, whose formula the
    many starts of the substance read (the load time of a model grows with
    every copy of a large formula).

    Args:
        substance: the substance.
        plasma: the kidney plasma at the steady state.
        eliminated: the net excretion of the kidney into the urine at the
            steady state (mmole/min); `rate`, `clearance` and `secretion` are
            then not needed.
        rate: the net source of the kidney (mmole/min).
        clearance: the first order clearance of the kidney plasma (l/min).
        secretion: a saturable secretion of the kidney plasma.
        variants: the variants of the definition (`PBPKDefinition.variants`),
            required: a variant without a closed steady state is refused.

    Returns:
        The formula.

    With a saturable reabsorption and no other process of the kidney the
    arterial plasma is the kidney plasma; a rate, a clearance, a secretion or
    `eliminated` are refused (see `kidney_plasma`).

    Raises:
        ValueError: for `eliminated` together with a rate, a clearance or a
            secretion.
        NotImplementedError: for a saturable reabsorption with a rate, a
            clearance, a secretion or `eliminated`.
    """
    _check_saturable(
        substance,
        eliminated is None and rate == clearance == "0" and secretion is None,
    )
    flow, _, rest = _kidney_flows(substance, variants)
    if substance.excretion.reabsorption is not None:
        return plasma
    if eliminated is not None:
        if rate != "0" or clearance != "0" or secretion is not None:
            raise ValueError(
                "With the elimination the arterial plasma needs no rate, "
                "clearance or secretion: `(E + R*C)/F`."
            )
        return f"({eliminated} + ({rest}) * ({plasma})) / ({flow})"
    saturable = (
        "0"
        if secretion is None
        else f"({secretion.vmax}) * ({plasma}) / ({secretion.km} + {plasma})"
    )
    linear = f"({clearance}) * ({plasma})" if clearance != "0" else "0"
    removed = _terms(linear, saturable)
    if rate != "0":
        removed = f"{removed} - ({rate})" if removed != "0" else f"-({rate})"
    if removed == "0":
        return plasma
    return f"{plasma} + ({removed}) / ({rest})"


def renal_clearance(
    substance: Substance,
    *,
    clearance: str = "0",
    arterial: str | None = None,
    variants: Variants,
) -> str:
    """Get the renal clearance of the arterial plasma at the steady state [l/min].

    The filtration and a secretion of the kidney plasma (`clearance`) in
    series: `G + (F - G)*CL/(F - G + CL)`. With a saturable reabsorption the
    clearance is the net filtered plasma flow at the arterial plasma
    (`renal_filtration`), without a secretion.

    Args:
        substance: the substance.
        clearance: the clearance of the kidney plasma (l/min).
        arterial: the arterial plasma of a saturable reabsorption, the
            baseline by default; a constant fraction does not read it.
        variants: the variants of the definition (`PBPKDefinition.variants`),
            required: a variant without a closed steady state is refused.

    Returns:
        The formula.

    Raises:
        NotImplementedError: for a saturable reabsorption with a clearance.
    """
    _check_saturable(substance, clearance == "0")
    if substance.excretion.reabsorption is not None:
        return renal_filtration(substance, arterial=arterial, variants=variants)
    _, gfr, rest = _kidney_flows(substance, variants)
    if clearance == "0":
        return gfr
    secreted = f"({rest}) * ({clearance}) / ({rest} + {clearance})"
    return _terms(gfr, secreted)


def liver_plasma(
    substance: Substance,
    *,
    arterial: str | None = None,
    portal: str | None = None,
    rate: str = "0",
    clearance: str = "0",
    variants: Variants,
) -> str:
    """Get the plasma of the liver at its steady state [mM].

    The liver is perfused by `(1 - f_shunt)` of the hepatic artery and the
    portal vein and drains `(1 - f_shunt)*Qh`; a source `rate` of the liver
    (mmole/min) leaves with it and a clearance of its plasma (l/min) removes
    a part: `Cpl = ((1 - f_shunt)*(Qha*Car + Qpo*Cpo)*c + rate)/((1 -
    f_shunt)*Qh*c + clearance)` with the carrier `c` (`plasma_flow`). The
    zonated liver (no source and no clearance, it refuses custom reactions)
    has the same plasma in every zone.

    Args:
        substance: the substance.
        arterial: the arterial plasma, the baseline by default.
        portal: the portal plasma, the baseline by default (the portal organs
            without a source).
        rate: the net source of the liver (mmole/min).
        clearance: the clearance of the liver plasma (l/min).
        variants: the variants of the definition (`PBPKDefinition.variants`),
            required: a variant without a closed steady state is refused.

    Returns:
        The formula.
    """
    if (rate != "0" or clearance != "0") and _variant(
        variants, Organ.LIVER
    ) is not None:
        raise NotImplementedError(
            "No steady state of a source or a clearance in the variant of the liver."
        )
    b = ids.baseline(substance.sid)
    car = b if arterial is None else arterial
    cpo = b if portal is None else portal
    perfused = f"({ONE} - {ids.F_SHUNT})"
    carrier = _carrier(substance)
    inflow = f"{perfused} * ({ids.QHA} * ({car}) + {ids.QPO} * ({cpo})) * {carrier}"
    flow = plasma_flow(substance, Organ.LIVER, variants=variants)
    return f"({_terms(inflow, rate)}) / ({_terms(flow, clearance)})"


def arm_vein(
    substance: Substance,
    preset: str,
    outflows: Mapping[Organ, str],
    *,
    arterial: str | None = None,
) -> str:
    """Get the plasma of the arm vein of a preset at the steady state [mM].

    The mix `sum w_arm_<site>*C_<site>` of the venous outflows of the organs
    (and the arterial plasma) of the preset (`PBPKDefinition.arm_vein`); an
    organ not in `outflows` and the arterial plasma are at the baseline.

    Args:
        substance: the substance.
        preset: the preset of the arm vein, e.g. `pksim_default`.
        outflows: the plasma leaving the organs at the steady state.
        arterial: the arterial plasma, the baseline by default.

    Returns:
        The formula.
    """
    b = ids.baseline(substance.sid)
    terms = []
    for site in arm_vein_weights(preset):
        if site is Pool.ARTERIAL:
            c = b if arterial is None else arterial
        else:
            assert isinstance(site, Organ)
            c = outflows.get(site, b)
        terms.append(f"{ids.arm_weight(site)} * ({c})")
    return " + ".join(terms)
