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
- `renal_filtration`, `kidney_plasma`, `renal_clearance`: the net filtrate of
  the kidney, the kidney plasma with a secretion (a clearance of the kidney
  plasma) and the renal clearance of the arterial plasma (glomerulus and
  tubule in series);
- `liver_plasma`: the liver plasma, fed by the hepatic artery and the portal
  vein past the shunt `f_shunt`, with a source of the liver;
- `arm_vein`: the mix of the venous outflows of the arm vein of a preset.

Every formula is the steady state of an organ at a given inflow (the
arterial plasma, the baseline by default): the standard organ is linear in
its species at a fixed inflow, so it is a closed form. The formulas read the
same ids as the equations of the organs (the flows, `HCT`, `BP_<s>`, `fu_<s>`,
`GFR`, `f_reabsorbed_<s>`, the coefficients `Kp`, `ftissue_<s>`,
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
the standard organ silently.
"""

from collections.abc import Mapping

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.definition import Variants
from pkpdsim.pbpk.observers import arm_vein_weights
from pkpdsim.pbpk.organs import Organ, Pool
from pkpdsim.pbpk.placement import kp_id
from pkpdsim.pbpk.substance import Mode, Substance

__all__ = [
    "arm_vein",
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

ONE = "1 dimensionless"
"""The number one in a formula, with its unit (SBML unit checks)."""

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


def plasma_flow(substance: Substance, organ: Organ) -> str:
    """Get the flow which carries the substance through a standard organ [l/min].

    `Q<o>*(1 - HCT)`, the plasma flow, or `Q<o>*BP_<s>` with erythrocytes in
    equilibrium (the plasma species carries the blood); the liver is perfused
    by `(1 - f_shunt)*Qh`.

    Args:
        substance: the substance.
        organ: the organ.

    Returns:
        The formula.

    Raises:
        NotImplementedError: with kinetic erythrocytes.
    """
    if organ is Organ.LIVER:
        return f"({ONE} - {ids.F_SHUNT}) * " + _flow_of(substance, ids.QH)
    return _flow_of(substance, ids.organ_flow(organ))


def tissue_volume(organ: Organ) -> str:
    """Get the volume of the merged tissue of a standard organ, `Vti` [l].

    `V<o>*(1 - f_vas_<o>)`, of the liver times `f_functional`: the local
    compartment `Vti` of the organ as a formula of the top model, e.g. for the
    rate of a custom reaction per litre of tissue at the steady state.

    Args:
        organ: the organ.

    Returns:
        The formula.
    """
    tissue = f"{ids.organ_volume(organ)} * ({ONE} - {ids.vascular_fraction(organ)})"
    if organ is Organ.LIVER:
        tissue = f"{tissue} * {ids.F_FUNCTIONAL}"
    return tissue


def tissue_exchange(substance: Substance, organ: Organ) -> str:
    """Get the exchange of the merged tissue of an organ with its plasma, `PS*fu` [l/min].

    The share of the whole-body exchange `ftissue_<s>` by the tissue volume of
    the organ (of the liver times `f_functional`) relative to the reference
    tissue volume, times the fraction unbound: the exchange
    `PS*fu*(Cpl - Cti/Kp)` of the merged tissue (`Mode.TISSUE`).

    Args:
        substance: the substance.
        organ: the organ.

    Returns:
        The formula.

    Raises:
        ValueError: if the substance does not enter the tissue of the organ.
        NotImplementedError: for another mode than `Mode.TISSUE` or an
            exchange of the organ of its own (`ps`).
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
    if organ in dist.ps:
        raise NotImplementedError(
            f"{substance.sid}: the exchange of the organ '{organ.name.lower()}' is a "
            "value of the organ (ps), no constant of the top model."
        )
    return (
        f"{ids.ftissue(substance.sid)} * {tissue_volume(organ)} / "
        f"{ids.tissue_volume_total_ref(substance.sid)} * {ids.fu(substance.sid)}"
    )


def _check_standard(organ: Organ, variants: Variants) -> None:
    """Refuse the organs whose steady state the plasma of a standard organ is not."""
    if organ in _GUT:
        raise NotImplementedError(
            f"No steady state of the organ '{organ.name.lower()}': the intestine refuses "
            "custom reactions, its parts start at the equilibrium with the baseline."
        )
    if organ is Organ.KIDNEY:
        raise ValueError("The kidney filters: use `kidney_plasma`.")
    if organ is Organ.LIVER:
        raise ValueError("The liver has two inflows and a shunt: use `liver_plasma`.")
    if _variant(variants, organ) is not None:
        raise NotImplementedError(
            f"No steady state of the variant of the organ '{organ.name.lower()}'."
        )


def organ_plasma(
    substance: Substance,
    organ: Organ,
    *,
    inflow: str | None = None,
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
        inflow: the plasma flowing in, the baseline `baseline_<s>` (the
            arterial plasma of the steady state) by default.
        rate: the net source of the organ at the steady state (mmole/min),
            negative for an uptake.
        clearance: the clearance of the plasma of the organ (l/min).
        variants: the variants of the definition (`PBPKDefinition.variants`),
            required: a variant without a closed steady state is refused.

    Returns:
        The formula.
    """
    _check_standard(organ, variants)
    flow = plasma_flow(substance, organ)
    c_in = ids.baseline(substance.sid) if inflow is None else inflow
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
    if organ in _GUT:
        _check_standard(organ, variants)
    variant = _variant(variants, organ)
    if variant is not None and not (organ is Organ.LIVER and rate == "0"):
        raise NotImplementedError(
            f"No steady state of a source in the variant of the organ '{organ.name.lower()}'."
        )
    exchange = tissue_exchange(substance, organ)
    kp = kp_id(substance, organ)
    if rate == "0":
        return f"{kp} * ({plasma})"
    return f"{kp} * ({plasma} + ({rate}) / ({exchange}))"


def renal_filtration(substance: Substance) -> str:
    """Get the net filtered plasma flow of the kidney, `fu*GFR*(1 - f_reabsorbed)` [l/min].

    Args:
        substance: the substance, with urinary excretion.

    Returns:
        The formula.

    Raises:
        ValueError: for a substance without urinary excretion.
    """
    s = substance.sid
    if not substance.excretion.urine:
        raise ValueError(f"{s} is not filtered (Excretion(urine=False)).")
    return f"{ids.fu(s)} * {ids.GFR} * ({ONE} - {ids.fraction_reabsorbed(s)})"


def _check_kidney(variants: Variants) -> None:
    """Refuse the nephron kidney."""
    if _variant(variants, Organ.KIDNEY) is not None:
        raise NotImplementedError(
            "No steady state of the nephron kidney: its tubule starts empty, not "
            "at its steady state."
        )


def kidney_plasma(
    substance: Substance,
    *,
    arterial: str | None = None,
    rate: str = "0",
    clearance: str = "0",
    variants: Variants,
) -> str:
    """Get the plasma of the kidney at its steady state [mM].

    The glomerulus filters the arterial plasma, the plasma leaving the
    kidney lacks the net filtrate `CL_gfr` (`renal_filtration`, 0 without
    urinary excretion), a secretion clears the kidney plasma (`clearance`) and
    a source returns `rate` into it (a reabsorption): `Cpl = ((F - CL_gfr)*Car
    + rate)/(F - CL_gfr + clearance)`.

    Args:
        substance: the substance.
        arterial: the arterial plasma, the baseline by default.
        rate: the net source of the kidney (mmole/min).
        clearance: the clearance of the kidney plasma (l/min).
        variants: the variants of the definition (`PBPKDefinition.variants`),
            required: a variant without a closed steady state is refused.

    Returns:
        The formula.
    """
    _check_kidney(variants)
    flow = plasma_flow(substance, Organ.KIDNEY)
    if substance.excretion.urine:
        flow = f"{flow} - {renal_filtration(substance)}"
    car = ids.baseline(substance.sid) if arterial is None else arterial
    return f"({_terms(f'({flow}) * ({car})', rate)}) / ({_terms(flow, clearance)})"


def renal_clearance(
    substance: Substance,
    *,
    clearance: str = "0",
    variants: Variants,
) -> str:
    """Get the renal clearance of the arterial plasma at the steady state [l/min].

    The filtration and a secretion of the kidney plasma (`clearance`) in
    series: `CL_gfr + (F - CL_gfr)*CL/(F - CL_gfr + CL)`.

    Args:
        substance: the substance.
        clearance: the clearance of the kidney plasma (l/min).
        variants: the variants of the definition (`PBPKDefinition.variants`),
            required: a variant without a closed steady state is refused.

    Returns:
        The formula.
    """
    _check_kidney(variants)
    flow = plasma_flow(substance, Organ.KIDNEY)
    gfr = renal_filtration(substance) if substance.excretion.urine else "0"
    if clearance == "0":
        return gfr
    rest = f"{flow} - ({gfr})" if gfr != "0" else flow
    secreted = f"({rest}) * ({clearance}) / ({rest} + {clearance})"
    return _terms(gfr, secreted)


def liver_plasma(
    substance: Substance,
    *,
    arterial: str | None = None,
    portal: str | None = None,
    rate: str = "0",
    variants: Variants,
) -> str:
    """Get the plasma of the liver at its steady state [mM].

    The liver is perfused by `(1 - f_shunt)` of the hepatic artery and the
    portal vein and drains `(1 - f_shunt)*Qh`; a source `rate` of the liver
    (mmole/min) leaves with it:
    `Cpl = ((1 - f_shunt)*(Qha*Car + Qpo*Cpo)*c + rate)/((1 - f_shunt)*Qh*c)`
    with the carrier `c` (`plasma_flow`). The zonated liver (no source, it
    refuses custom reactions) has the same plasma in every zone.

    Args:
        substance: the substance.
        arterial: the arterial plasma, the baseline by default.
        portal: the portal plasma, the baseline by default (the portal organs
            without a source).
        rate: the net source of the liver (mmole/min).
        variants: the variants of the definition (`PBPKDefinition.variants`),
            required: a variant without a closed steady state is refused.

    Returns:
        The formula.
    """
    if rate != "0" and _variant(variants, Organ.LIVER) is not None:
        raise NotImplementedError(
            "No steady state of a source in the variant of the liver."
        )
    b = ids.baseline(substance.sid)
    car = b if arterial is None else arterial
    cpo = b if portal is None else portal
    perfused = f"({ONE} - {ids.F_SHUNT})"
    carrier = _carrier(substance)
    inflow = f"{perfused} * ({ids.QHA} * ({car}) + {ids.QPO} * ({cpo})) * {carrier}"
    return f"({_terms(inflow, rate)}) / ({plasma_flow(substance, Organ.LIVER)})"


def arm_vein(
    preset: str,
    substance: Substance,
    outflows: Mapping[Organ, str],
    *,
    arterial: str | None = None,
) -> str:
    """Get the plasma of the arm vein of a preset at the steady state [mM].

    The mix `sum w_arm_<site>*C_<site>` of the venous outflows of the organs
    (and the arterial plasma) of the preset (`PBPKDefinition.arm_vein`); an
    organ not in `outflows` and the arterial plasma are at the baseline.

    Args:
        preset: the preset of the arm vein, e.g. `pksim_default`.
        substance: the substance.
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
