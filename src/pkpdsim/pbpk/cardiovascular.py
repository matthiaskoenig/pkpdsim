"""Optional pressure-driven whole-body circulation in SBML.

The resistance network follows the PBPK organ topology, including the portal
organs and the two hepatic supplies. ``PressureFlow`` solves its zero-storage
junctions at prescribed arterial and venous pressures. ``DynamicHeart`` closes
the loop with compliant blood reservoirs and time-varying elastance ventricles.
Algebraic junction balances are solved by sbmlsim's core engine (SUNDIALS IDA).

The equations are an independent implementation of hydraulic circuit and
elastance methods, not a translation of an external library's parameter set.
See ``docs/pbpk/cardiovascular.md`` for assumptions, sources and examples.
"""

import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from sbmlutils.factory import (
    AlgebraicRule,
    AssignmentRule,
    Compartment,
    InitialAssignment,
    Parameter,
    RateRule,
    Reaction,
    Sbase,
    UnitDefinition,
)
from sbmlutils.factory import (
    Species as SBMLSpecies,
)
from sbmlutils.metadata import BQB, SBO

from pkpdsim.database.tables import load
from pkpdsim.pbpk import ids
from pkpdsim.pbpk.annotation import annotate
from pkpdsim.pbpk.organs import (
    PORTAL_ORGANS,
    SYSTEMIC_ORGANS,
    Organ,
    Part,
    Pool,
    as_organ_keys,
)
from pkpdsim.pbpk.placement import species_start, start_formula
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import FrozenMap
from pkpdsim.species import Species

if TYPE_CHECKING:
    from pkpdsim.pbpk.definition import PBPKDefinition

__all__ = ["CVU", "DynamicHeart", "PressureFlow"]


class CVU(U):
    """Pressure, resistance, compliance and elastance in PBPK time/volume units."""

    mmHg = UnitDefinition("mmHg", "133.322387415 pascal", name="mmHg")
    resistance = UnitDefinition("cv_resistance", "133.322387415 pascal*min/l")
    compliance = UnitDefinition("cv_compliance", "l/(133.322387415 pascal)")
    elastance = UnitDefinition("cv_elastance", "133.322387415 pascal/l")


_DEFAULTS = load("cardiovascular")["defaults"]
_UNITS = {
    "mmHg": CVU.mmHg,
    "mmHg*min/l": CVU.resistance,
    "l/mmHg": CVU.compliance,
    "mmHg/l": CVU.elastance,
    "l": U.litre,
    "1/min": U.per_min,
    "dimensionless": U.dimensionless,
}
_EXTRA_POOLS = ("pa", "pv", "lv", "rv")
_FLOW_IDS = (
    ids.CO,
    *(ids.organ_flow(o) for o in Organ),
    ids.QHA,
    ids.QPO,
    ids.QH,
)


@dataclass(frozen=True)
class PressureFlow:
    """A resistance network at prescribed systemic arterial/venous pressures.

    Attributes:
        parameters: overrides of the short keys of ``cardiovascular.toml``
            (e.g. ``P_ar`` in mmHg). Each is an illustrative model choice.
        resistance: positive resistance multipliers by systemic organ; liver
            means the hepatic artery. Change hepatic outflow with
            ``hepatic_resistance``. The lung multiplier is used by DynamicHeart.
        systemic_resistance: multiplier of the arterial organ resistances.
        hepatic_resistance: multiplier of the sinusoidal-to-venous resistance.
        portal_resistance: multiplier of the portal-to-sinusoidal resistance.
    """

    parameters: Mapping[str, float] = field(default_factory=dict)
    resistance: Mapping[Organ, float] = field(default_factory=dict)
    systemic_resistance: float = 1.0
    hepatic_resistance: float = 1.0
    portal_resistance: float = 1.0

    def __post_init__(self) -> None:
        """Freeze and validate the configuration and reference pressure ordering."""
        unknown = set(self.parameters) - set(_DEFAULTS)
        if unknown:
            raise ValueError(f"Unknown cardiovascular parameters: {sorted(unknown)}.")
        parameters = {k: float(v) for k, v in self.parameters.items()}
        resistance = {
            k: float(v)
            for k, v in as_organ_keys(
                self.resistance, "PressureFlow.resistance"
            ).items()
        }
        for key, value in {
            **parameters,
            **{f"resistance {k}": v for k, v in resistance.items()},
            "systemic_resistance": self.systemic_resistance,
            "hepatic_resistance": self.hepatic_resistance,
            "portal_resistance": self.portal_resistance,
        }.items():
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"Cardiovascular {key} must be finite and positive.")
        values = {k: float(v["value"]) for k, v in _DEFAULTS.items()} | parameters
        if not values["P_ar"] > values["P_po"] > values["P_li"] > values["P_ve"]:
            raise ValueError(
                "Reference pressures must satisfy P_ar > P_po > P_li > P_ve."
            )
        if values["P_pa"] <= values["P_pv"]:
            raise ValueError("Reference pressures must satisfy P_pa > P_pv.")
        if not 0 < values["systole"] < 1:
            raise ValueError("The systolic fraction must be between zero and one.")
        if values["f_pa"] + values["f_pv"] >= 1:
            raise ValueError("Pulmonary volume fractions must sum to less than one.")
        for ventricle in ("lv", "rv"):
            if values[f"E_{ventricle}_max"] <= values[f"E_{ventricle}_min"]:
                raise ValueError(f"E_{ventricle}_max must exceed E_{ventricle}_min.")
            if values[f"V_{ventricle}"] <= values[f"V0_{ventricle}"]:
                raise ValueError(f"V_{ventricle} must exceed V0_{ventricle}.")
        object.__setattr__(self, "parameters", FrozenMap(parameters))
        object.__setattr__(
            self, "resistance", FrozenMap({k.value: v for k, v in resistance.items()})
        )


@dataclass(frozen=True)
class DynamicHeart(PressureFlow):
    """A closed circulation with compliant vessels and two elastance ventricles.

    Parameters of ``PressureFlow`` also apply. ``contractility`` scales the
    active ventricular elastance and is exposed as ``CV_contractility`` for PD.
    The adult-human defaults are illustrative; a periodic state is obtained by
    simulating several cycles before dosing, not by a stationary solver.
    """

    contractility: float = 1.0

    def __post_init__(self) -> None:
        """Validate the inherited network and the contractility multiplier."""
        super().__post_init__()
        if not math.isfinite(self.contractility) or self.contractility <= 0:
            raise ValueError(
                "Cardiovascular contractility must be finite and positive."
            )


def check_circulation(definition: PBPKDefinition) -> None:
    """Reject topologies and initialization contracts the circulation cannot satisfy."""
    circulation = definition.variants.circulation
    if circulation is None:
        return
    if not isinstance(circulation, PressureFlow):
        raise ValueError("Variants.circulation must be PressureFlow or DynamicHeart.")
    if definition.pregnancy is not None:
        raise ValueError(
            "Cardiovascular variants do not yet include placental circulation."
        )
    if isinstance(circulation, DynamicHeart):
        if definition.species is not Species.HUMAN or definition.age:
            raise ValueError(
                "DynamicHeart currently requires an adult human definition."
            )
        if definition.adipose is not None and definition.adipose.effect is not None:
            raise ValueError("DynamicHeart requires a constant body size.")
        changing = set(definition.dynamic) | {
            output for pd in definition.pd for output in pd.outputs.values()
        }
        if changing & {ids.BW, ids.HEIGHT, ids.VBLOOD, ids.HCT}:
            raise ValueError("DynamicHeart requires constant anatomy and hematocrit.")
        if any(
            getattr(p, "stationary", None) is not None for p in definition.chemistry
        ):
            raise ValueError(
                "DynamicHeart cannot use generated stationary chemistry; "
                "initialize by a finite periodic presimulation."
            )
        if _remodelling_lobes(definition):
            raise ValueError(
                "DynamicHeart requires a constant anatomy, which a lobed "
                "liver with a remodelling is not; use PressureFlow for the "
                "days of a remodelling or drop the remodelling."
            )


def _remodelling_lobes(definition: PBPKDefinition) -> bool:
    """Whether the liver of a definition is a lobed liver whose anatomy moves."""
    from pkpdsim.pbpk.builder import organ_component
    from pkpdsim.pbpk.variants.liver_lobes import LobedLiver

    liver = organ_component(definition, Organ.LIVER)
    return isinstance(liver, LobedLiver) and liver.remodelling is not None


def _parameter(sid: str, value: float, unit: UnitDefinition, notes: str) -> Parameter:
    """A named cardiovascular parameter with explicit model-choice provenance."""
    return Parameter(
        sid,
        value,
        unit,
        name=sid.replace("CV_", "cardiovascular "),
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        annotations=[(BQB.IS_DESCRIBED_BY, "doi/10.1186/1475-925X-10-33")],
        notes=notes,
    )


def _derived(
    sid: str, formula: str, unit: UnitDefinition, *, initial: bool = False
) -> list[Sbase]:
    """A cardiovascular quantity with an initial assignment or assignment rule."""
    return [
        _parameter(sid, 0.0, unit, f"`{sid} = {formula}`"),
        InitialAssignment(sid, formula, unit)
        if initial
        else AssignmentRule(sid, formula, unit),
    ]


def lobed_liver(definition: PBPKDefinition) -> bool:
    """Whether the liver of a definition is a lobed liver of more than one lobe.

    A lobed liver brings its own vascular network, one node and three
    branches per lobe (`pkpdsim.pbpk.variants.liver_lobes.pressure_objects`),
    so the circulation leaves the hepatic part of the whole liver network to
    it: the single sinusoidal node, the two hepatic resistances, the three
    hepatic flows and the junction balance.

    Args:
        definition: the definition.

    Returns:
        Whether the liver has lobes.

    Raises:
        ValueError: for a species without a lobe set.
    """
    from pkpdsim.pbpk.builder import organ_component
    from pkpdsim.pbpk.variants.liver_lobes import LobedLiver

    liver = organ_component(definition, Organ.LIVER)
    if not isinstance(liver, LobedLiver):
        return False
    return len(liver.lobe_set(definition.species)) > 1


def _sum_flows(organs: frozenset[Organ]) -> str:
    """A stable sum of the flows of a set of organs."""
    return " + ".join(ids.organ_flow(o) for o in Organ if o in organs)


def circulation_physiology(definition: PBPKDefinition, objects: list[Sbase]) -> None:
    """Replace prescribed flows with a pressure network, preserving public PBPK ids."""
    circulation = definition.variants.circulation
    if circulation is None:
        return
    check_circulation(definition)
    dynamic = isinstance(circulation, DynamicHeart)
    # the lobes of a lobed liver bring their own nodes and branches
    lobed = lobed_liver(definition)
    remove = set(_FLOW_IDS)
    if dynamic:
        remove |= {ids.pool_volume(p) for p in Pool}
    objects[:] = [
        o
        for o in objects
        if not isinstance(o, InitialAssignment) or o.symbol not in remove
    ]
    for obj in objects:
        if isinstance(obj, Parameter) and obj.sid in remove:
            obj.constant = False
            obj.value = 0.0
            if obj.sid in _FLOW_IDS:
                obj.notes = (
                    "Pressure-driven flow of the optional cardiovascular network."
                )
    for key, entry in _DEFAULTS.items():
        value = circulation.parameters.get(key, entry["value"])
        # Reference pressures stay fixed while actual pressure boundaries may be changed.
        suffix = "_ref" if key.startswith("P_") else ""
        objects.append(
            _parameter(
                f"CV_{key}{suffix}",
                value,
                _UNITS[entry["unit"]],
                f"Model choice (cardiovascular.toml): {entry['notes']}",
            )
        )
    for key, value in (
        ("systemic", circulation.systemic_resistance),
        ("hepatic", circulation.hepatic_resistance),
        ("portal", circulation.portal_resistance),
        *((o.value, circulation.resistance.get(o, 1.0)) for o in Organ),
    ):
        objects.append(
            _parameter(
                f"CV_fR_{key}",
                value,
                U.dimensionless,
                "Positive resistance multiplier; may be replaced by a PD output.",
            )
        )
    for node in ("ar", "ve", "po", "li"):
        if dynamic and node != "li":
            continue
        entry = _DEFAULTS[f"P_{node}"]
        par = _parameter(
            f"CV_P_{node}",
            circulation.parameters.get(f"P_{node}", entry["value"]),
            CVU.mmHg,
            "Systemic pressure boundary or algebraic junction pressure.",
        )
        par.constant = node in ("ar", "ve")
        objects.append(par)
    for organ in Organ:
        if organ is Organ.LUNG or (organ is Organ.LIVER and lobed):
            continue
        down = (
            "po" if organ in PORTAL_ORGANS else "li" if organ is Organ.LIVER else "ve"
        )
        reference = f"{ids.flow_fraction(organ)} * {ids.CO_BASE}"
        objects += _derived(
            f"CV_R_{organ.value}",
            f"(CV_P_ar_ref - CV_P_{down}_ref)/({reference})",
            CVU.resistance,
            initial=True,
        )
        objects.append(
            AssignmentRule(
                ids.organ_flow(organ),
                f"max(CV_P_ar - CV_P_{down}, 0 mmHg) * "
                f"{ids.impairment_factor(ids.organ_flow(organ))} / "
                f"(CV_R_{organ.value} * CV_fR_{organ.value} * CV_fR_systemic)",
                U.l_per_min,
            )
        )
    portal_ref = " + ".join(
        f"{ids.flow_fraction(o)} * {ids.CO_BASE}" for o in Organ if o in PORTAL_ORGANS
    )
    if not lobed:
        objects += _derived(
            "CV_R_portal",
            f"(CV_P_po_ref - CV_P_li_ref)/({portal_ref})",
            CVU.resistance,
            initial=True,
        )
        objects += _derived(
            "CV_R_hepatic",
            f"(CV_P_li_ref - CV_P_ve_ref)/"
            f"({portal_ref} + {ids.flow_fraction(Organ.LIVER)} * {ids.CO_BASE})",
            CVU.resistance,
            initial=True,
        )
        objects += [
            AssignmentRule(ids.QHA, ids.organ_flow(Organ.LIVER), U.l_per_min),
            AssignmentRule(
                ids.QPO,
                "max(CV_P_po - CV_P_li, 0 mmHg)/(CV_R_portal * CV_fR_portal)",
                U.l_per_min,
            ),
            AssignmentRule(
                ids.QH,
                "max(CV_P_li - CV_P_ve, 0 mmHg)/(CV_R_hepatic * CV_fR_hepatic)",
                U.l_per_min,
            ),
            AlgebraicRule(
                "CV_hepatic_balance",
                f"{ids.QHA} + {ids.QPO} - {ids.QH}",
                U.l_per_min,
                name="hepatic junction volume balance",
            ),
        ]
    if dynamic:
        _heart_physiology(circulation, objects)
    else:
        objects += [
            AlgebraicRule(
                "CV_portal_balance",
                f"{_sum_flows(PORTAL_ORGANS)} - {ids.QPO}",
                U.l_per_min,
                name="portal junction volume balance",
            ),
            AssignmentRule(
                ids.CO, f"{_sum_flows(SYSTEMIC_ORGANS)} + {ids.QH}", U.l_per_min
            ),
            AssignmentRule(ids.organ_flow(Organ.LUNG), ids.CO, U.l_per_min),
        ]
    # sbmlutils does not infer parameter.constant from an assignment rule.
    ruled = {o.variable for o in objects if isinstance(o, AssignmentRule | RateRule)}
    for obj in objects:
        if isinstance(obj, Parameter) and obj.sid in ruled:
            obj.constant = False
        if isinstance(
            obj, InitialAssignment | AssignmentRule | RateRule | AlgebraicRule
        ):
            obj.name = (
                obj.name
                or f"cardiovascular {getattr(obj, 'variable', None) or obj.sid}"
            )
            obj.sboTerm = obj.sboTerm or SBO.MATHEMATICAL_EXPRESSION


def _heart_physiology(circulation: DynamicHeart, objects: list[Sbase]) -> None:
    """Compliant central reservoirs, ventricular storage and ideal diode valves."""
    objects += _derived(
        "CV_scale", f"{ids.VBLOOD}/CV_Vblood_ref", U.dimensionless, initial=True
    )
    # Reserve chamber and pulmonary volumes from existing central pools, so
    # adding a heart adds neither blood nor drug capacity to the whole body.
    for pool in Pool:
        objects += _derived(
            f"CV_V_{pool.value}_original",
            f"{ids.pool_fraction(pool)} * {ids.VBLOOD}",
            U.litre,
            initial=True,
        )
    for node in ("lv", "rv"):
        # The table's V_lv/V_rv are configurable initial volumes, not states.
        par = next(
            o for o in objects if isinstance(o, Parameter) and o.sid == f"CV_V_{node}"
        )
        par.sid = f"CV_V_{node}_initial"
        objects += _derived(
            f"CV_V_{node}", f"CV_V_{node}_initial * CV_scale", U.litre, initial=True
        )
    for node in ("pa", "pv"):
        objects += _derived(
            f"CV_V_{node}", f"CV_f_{node} * CV_V_ve_original", U.litre, initial=True
        )
    initial_volumes = {
        "ar": "CV_V_ar_original - CV_V_lv_initial * CV_scale",
        "ve": "CV_V_ve_original * (1 dimensionless - CV_f_pa - CV_f_pv) - CV_V_rv_initial * CV_scale",
        "po": "CV_V_po_original",
    }
    for node, formula in initial_volumes.items():
        objects.append(InitialAssignment(f"V{node}", formula, U.litre))
    derivatives = {
        "ar": f"CV_Q_av - ({_sum_flows(SYSTEMIC_ORGANS)} + {ids.QHA} + {_sum_flows(PORTAL_ORGANS)})",
        "ve": f"{_sum_flows(SYSTEMIC_ORGANS)} + {ids.QH} - CV_Q_tv",
        "po": f"{_sum_flows(PORTAL_ORGANS)} - {ids.QPO}",
        "pa": "CV_Q_pv - Qlu",
        "pv": "Qlu - CV_Q_mv",
        "lv": "CV_Q_mv - CV_Q_av",
        "rv": "CV_Q_tv - CV_Q_pv",
    }
    for node, derivative in derivatives.items():
        volume = f"V{node}" if node in initial_volumes else f"CV_V_{node}"
        objects.append(RateRule(volume, derivative, U.l_per_min))
        if node in ("lv", "rv"):
            continue
        initial = initial_volumes.get(node, f"CV_f_{node} * CV_V_ve_original")
        objects += _derived(f"CV_V_{node}_ref", initial, U.litre, initial=True)
        objects += _derived(
            f"CV_P_{node}",
            f"CV_P_{node}_ref + ({volume} - CV_V_{node}_ref)/(CV_C_{node} * CV_scale)",
            CVU.mmHg,
        )
    objects += [
        _parameter(
            "CV_phase", 0.0, U.dimensionless, "Integrated cardiac phase in cycles."
        ),
        RateRule("CV_phase", "CV_HR", U.per_min),
        _parameter(
            "CV_contractility",
            circulation.contractility,
            U.dimensionless,
            "Multiplier of active ventricular elastance; a target for PD.",
        ),
    ]
    objects += _derived(
        "CV_activation",
        "piecewise((1 dimensionless - cos(2 dimensionless * pi * "
        "(CV_phase - floor(CV_phase))/CV_systole))/2 dimensionless, "
        "CV_phase - floor(CV_phase) < CV_systole, 0 dimensionless)",
        U.dimensionless,
    )
    for node in ("lv", "rv"):
        objects += _derived(
            f"CV_E_{node}",
            f"(CV_E_{node}_min + (CV_E_{node}_max - CV_E_{node}_min) * "
            f"CV_contractility * CV_activation)/CV_scale",
            CVU.elastance,
        )
        objects += _derived(
            f"CV_P_{node}",
            f"CV_E_{node} * (CV_V_{node} - CV_V0_{node} * CV_scale)",
            CVU.mmHg,
        )
    for valve, source, target in (
        ("av", "lv", "ar"),
        ("tv", "ve", "rv"),
        ("pv", "rv", "pa"),
        ("mv", "pv", "lv"),
    ):
        objects += _derived(
            f"CV_Q_{valve}",
            f"max(CV_P_{source} - CV_P_{target}, 0 mmHg) * CV_scale/CV_R_valve",
            U.l_per_min,
        )
    objects += _derived(
        "CV_R_lu",
        "(CV_P_pa_ref - CV_P_pv_ref)/CO_base",
        CVU.resistance,
        initial=True,
    )
    objects += [
        AssignmentRule(
            "Qlu", "max(CV_P_pa - CV_P_pv, 0 mmHg)/(CV_R_lu * CV_fR_lu)", U.l_per_min
        ),
        AssignmentRule(ids.CO, "CV_Q_av", U.l_per_min),
    ]
    stores = "Var + Vve + Vpo + " + " + ".join(f"CV_V_{n}" for n in _EXTRA_POOLS)
    original = "CV_V_ar_original + CV_V_ve_original + CV_V_po_original"
    objects += _derived("CV_V_total", stores, U.litre)
    objects += _derived(
        "CV_V_blood_total", f"CV_V_total + Vblood - ({original})", U.litre
    )


def circulation_dynamic(definition: PBPKDefinition) -> tuple[str, ...]:
    """The seeds whose dependent organ quantities must follow pressure-driven flows."""
    if definition.variants.circulation is None:
        return ()
    volumes = (
        tuple(ids.pool_volume(p) for p in Pool)
        if isinstance(definition.variants.circulation, DynamicHeart)
        else ()
    )
    return (*_FLOW_IDS, *volumes)


def _pool_species(node: str, substance: str, cells: bool = False) -> str:
    """The concentration id of an additional cardiovascular blood pool."""
    return f"CV_C_{node}_{'bc_' if cells else ''}{substance}"


def circulation_transport(definition: PBPKDefinition) -> list[Sbase]:
    """Drug amounts move through the same valves and variable volumes as blood."""
    if not isinstance(definition.variants.circulation, DynamicHeart):
        return []
    objects: list[Sbase] = []
    for node in _EXTRA_POOLS:
        for part, factor in (("pl", "(1 dimensionless - HCT)"), ("bc", "HCT")):
            cid = f"CV_V_{node}_{part}"
            objects += [
                Compartment(
                    cid, None, U.litre, constant=False, name=f"{node} {part} blood"
                ),
                AssignmentRule(cid, f"CV_V_{node} * {factor}", U.litre),
            ]
        for s in definition.substances:
            ery = s.erythrocytes
            equilibrium = ery is not None and ery.mode == "equilibrium"
            kinetic = ery is not None and ery.mode == "kinetic"
            compartment = f"CV_V_{node}_pl"
            if equilibrium:
                compartment = f"CV_V_{node}_pl_{s.sid}"
                objects += [
                    Compartment(compartment, None, U.litre, constant=False),
                    AssignmentRule(
                        compartment, f"CV_V_{node} * {ids.bp(s.sid)}", U.litre
                    ),
                ]
            for cells in (False, True) if kinetic else (False,):
                sid = _pool_species(node, s.sid, cells)
                comp = f"CV_V_{node}_bc" if cells else compartment
                start = start_formula(s, ids.kbc(s.sid)) if cells else start_formula(s)
                given = definition.starts.get(sid)
                if given is not None:
                    raise ValueError(
                        f"Set cardiovascular start '{sid}' through a simulator preinitialization change."
                    )
                objects.append(
                    SBMLSpecies(
                        sid,
                        comp,
                        initialConcentration=0.0 if start is None else None,
                        substanceUnit=U.mmole,
                        hasOnlySubstanceUnits=False,
                        name=f"{s.name} in cardiovascular {node}",
                        sboTerm=SBO.SIMPLE_CHEMICAL,
                    )
                )
                if start is not None:
                    objects.append(species_start(sid, start, f"initial {sid}"))
            if kinetic:
                pl, bc = _pool_species(node, s.sid), _pool_species(node, s.sid, True)
                objects.append(
                    Reaction(
                        f"CV_J_exchange_{node}_bc_{s.sid}",
                        f"{pl} <-> {bc}",
                        f"{ids.kbc_exchange(s.sid)} * CV_V_{node}_bc * ({pl} - {bc}/{ids.kbc(s.sid)})",
                    )
                )
    for s in definition.substances:
        ery = s.erythrocytes
        equilibrium = ery is not None and ery.mode == "equilibrium"
        kinetic = ery is not None and ery.mode == "kinetic"
        for valve, source, target in (
            ("av", "lv", "ar"),
            ("tv", "ve", "rv"),
            ("pv", "rv", "pa"),
            ("mv", "pv", "lv"),
        ):
            for cells in (False, True) if kinetic else (False,):
                part = Part.BLOOD_CELLS if cells else Part.PLASMA

                src, dst = (
                    ids.pool_species(Pool(node), s.sid, part)
                    if node in ("ar", "ve")
                    else _pool_species(node, s.sid, cells)
                    for node in (source, target)
                )
                fraction = (
                    "HCT"
                    if cells
                    else ids.bp(s.sid)
                    if equilibrium
                    else "(1 dimensionless - HCT)"
                )
                objects.append(
                    Reaction(
                        f"CV_J_{valve}_{'bc_' if cells else ''}{s.sid}",
                        f"{src} -> {dst}",
                        f"CV_Q_{valve} * {fraction} * {src}",
                    )
                )
    for obj in objects:
        if isinstance(obj, Compartment | SBMLSpecies):
            obj.metaId = f"meta_{obj.sid}"
        if isinstance(obj, Compartment):
            obj.name = obj.name or f"cardiovascular {obj.sid}"
            annotate(
                obj,
                "blood_cells"
                if str(obj.sid).endswith("_bc")
                else "plasma"
                if str(obj.sid).endswith("_pl")
                else "capacity",
            )
        if isinstance(obj, Reaction):
            obj.name = f"cardiovascular {obj.sid}"
            annotate(
                obj,
                "exchange"
                if str(obj.sid).startswith("CV_J_exchange_")
                else "blood_flow",
            )
        if isinstance(obj, InitialAssignment | AssignmentRule):
            obj.name = f"cardiovascular {getattr(obj, 'variable', None) or obj.sid}"
            obj.sboTerm = SBO.MATHEMATICAL_EXPRESSION
    return objects


def circulation_boundary(
    definition: PBPKDefinition, organ: Organ, replacements: Mapping[str, str]
) -> Mapping[str, str]:
    """Connect the lung to pulmonary reservoirs when the heart is resolved."""
    if organ is not Organ.LUNG or not isinstance(
        definition.variants.circulation, DynamicHeart
    ):
        return replacements
    mapping = {}
    for node, pool in (("pa", Pool.VENOUS), ("pv", Pool.ARTERIAL)):
        for part in (Part.PLASMA, Part.BLOOD_CELLS):
            mapping[ids.pool_compartment(pool, part)] = f"CV_V_{node}_{part.value}"
            for s in definition.substances:
                mapping[ids.pool_species(pool, s.sid, part)] = _pool_species(
                    node, s.sid, part is Part.BLOOD_CELLS
                )
                mapping[ids.pool_capacity(pool, s.sid)] = f"CV_V_{node}_pl_{s.sid}"
    return {local: mapping.get(top, top) for local, top in replacements.items()}
