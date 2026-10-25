"""The impaired organ function: parameter vectors on the physiology.

Impairment is never a structurally different model: three severities of the
top model set factors on the quantities of the healthy physiology, the
liver and the binding to plasma proteins (`data/impairment.toml`, every
number with its source; a value which no publication read gives is a model
choice):

- `CP_severity`: liver cirrhosis, 0 healthy, 1, 2, 3 the Child-Pugh classes
  A, B, C; a severity between the classes interpolates every factor linearly
  (scans and fits use intermediate values). It sets the functional fraction
  `f_functional` and the portosystemic shunt `f_shunt` of the liver, the
  hepatic arterial and the renal blood flow, the cardiac output (cardiac
  index; the organs without a factor of their own close the cardiac output to
  it), the GFR, albumin and alpha-1-acid glycoprotein (and so fu) and the
  activity of the enzymes and transporters of the liver.
- `f_renal_function`: chronic kidney disease by the GFR relative to healthy
  (the factor of the physiology); its GFR per 1.73 m2 (`eGFR`) gives the
  stage (normal, mild, moderate, severe, end stage), which sets the binding
  capacity of the plasma proteins (and so fu); the factors are steps at the
  stage limits, a scan or fit of `f_renal_function` across a limit sees a
  discontinuity. The switch `CKD_hepatic` (0 off, 1 on) adds the reduced
  activity of hepatic CYP2D6 and OATP1B.
- `HF_severity`: heart failure, 0 healthy, 1 to 4 the NYHA classes I to IV,
  interpolated: the blood flows of liver and portal organs, kidney and skin,
  adipose and muscle are reduced, the cardiac output is the sum of the flows,
  the GFR follows the renal blood flow (which is not monotonic in the class,
  NYHA IV above NYHA III, as measured by Leithe 1984).

The hematocrit (anaemia of cirrhosis and CKD, combined multiplicatively on
the hematocrit of the sex of the build), the kidney volume and blood flow and
the gastric emptying and small intestinal transit of CKD have their factors,
which are 1 until a publication with the healthy reference is read (the data
say why).

Each condition has its factors `fcp_<key>`, `fckd_<key>`, `fhf_<key>`; the
factor of all conditions on a quantity is `fimp_<id>` (`fimp_GFR`,
`fimp_Qki`, `fimp_HCT`), their product. Blood flows: every organ flow is its
fraction of the cardiac output without impairment (`CO_base`) times
`fimp_Q<o>`, and the cardiac output is `CO_base` times
`fimp_CO = sum fQ<o>*fimp_Q<o>`, so the flows sum to it.
Plasma binding (McNamara 2002): `fu = 1/(1 + f*(1 - fu_ref)/fu_ref)` with the
binding factor `f`, the factor of the binding protein (`fimp_albumin`,
`fimp_agp`, 1 for an unknown protein) over the CKD factor of the binding
capacity `fimp_plasma_protein`: fu of a strongly bound substance rises by the
CKD factor, and fu stays below 1. Enzymes and transporters: a rate of a
protein in the liver runs with `f_<protein>*fimp_<protein>`, the organ-local
base activity (genotype, interactions of an experiment, `LI__f_cyp3a4`)
times the input `fimp_<protein>` of the impairment (`fimp_li_<protein>` of
the body), so both act multiplicatively.

All are initial assignments (`pkpdsim.pbpk.physiology`): set the severities
before the simulation, as a covariate; a severity which changes during a
simulation is named in `PBPKDefinition.dynamic`. `f_shunt` and
`f_functional` stay parameters an experiment may set directly. All
severities 0 give the healthy model: every factor is 1.
"""

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from functools import cache
from typing import TYPE_CHECKING, Any

from sbmlutils.factory import Parameter, Sbase, UnitDefinition
from sbmlutils.metadata import SBO

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.annotation import annotate, initial_assignment
from pkpdsim.pbpk.chemistry import Metabolism, Transport
from pkpdsim.pbpk.data import cite, load
from pkpdsim.pbpk.organs import PORTAL_ORGANS, Organ
from pkpdsim.pbpk.physiology import normal_gfr
from pkpdsim.pbpk.proteins import PROTEINS
from pkpdsim.pbpk.substance import Substance
from pkpdsim.pbpk.units import U

if TYPE_CHECKING:
    from pkpdsim.pbpk.definition import PBPKDefinition

__all__ = [
    "GASTRIC_EMPTYING",
    "SI_TRANSIT",
    "impairment_objects",
    "organ_activity",
    "physiology_factors",
    "unbound_fraction",
]

DATA = "impairment"
"""The data file of the impairment."""

GASTRIC_EMPTYING = ids.impairment_factor("t_stomach")
"""The factor on the transit time of the stomach, an input of the intestine."""

SI_TRANSIT = ids.impairment_factor("t_si")
"""The factor on the transit times of the small intestine, an input of the intestine."""

_ONE = "1 dimensionless"
_SBO = SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER
_STAGES = ("mild", "moderate", "severe", "end stage")
"""The CKD stages of the vectors of `[ckd]`; normal renal function is 1."""

_CP_FIXED = frozenset({*PORTAL_ORGANS, Organ.BRAIN})
"""The organs whose flow cirrhosis does not change: the portal organs (the
shunt represents the portal blood which does not reach the liver, see the
data) and the brain (autoregulation of the cerebral blood flow, a model
choice)."""

_HF_ORGANS: dict[str, frozenset[Organ]] = {
    "liver_flow": frozenset({Organ.LIVER, *PORTAL_ORGANS}),
    "kidney_flow": frozenset({Organ.KIDNEY}),
    "peripheral_flow": frozenset({Organ.SKIN, Organ.ADIPOSE, Organ.MUSCLE}),
}
"""The organs whose flow a factor of heart failure multiplies (Gu 2025: the
flows of the gastrointestinal organs and the spleen follow the liver; the
pancreas, which Gu 2025 do not model, is a portal organ as well)."""


@cache
def _data() -> dict[str, Any]:
    """The data of the impairment."""
    return load(DATA)


def _proteins() -> frozenset[str]:
    """The proteins with a factor of impairment in the liver."""
    data = _data()
    return frozenset(data["child_pugh"]["proteins"]) | frozenset(
        data["ckd"]["proteins"]
    )


def organ_activity(organ: Organ, protein: str) -> str | None:
    """Get the local input of the impairment on a protein in an organ.

    Args:
        organ: the organ.
        protein: the protein of a process of the organ.

    Returns:
        `fimp_<protein>` in the liver for a protein with data of impairment,
        replaced by `fimp_li_<protein>` of the body; `None` elsewhere.
    """
    if organ is Organ.LIVER and protein in _proteins():
        return ids.impairment_activity(protein)
    return None


def _flow_organs() -> list[Organ]:
    """The organs in parallel, whose flows are fractions of cardiac output."""
    return [o for o in Organ if o is not Organ.LUNG]


def physiology_factors() -> dict[str, str]:
    """Get the factors of the impairment on the quantities of the physiology.

    Returns:
        The factor (formula) by the id of the quantity it multiplies, for
        `pkpdsim.pbpk.physiology.physiology_objects`: the cardiac output (the
        sum of the impaired flows), the flows, the GFR, the hematocrit and the
        kidney mass.
    """
    co = ids.impairment_factor(ids.CO)
    factors = {
        ids.CO: co,
        ids.GFR: ids.impairment_factor(ids.GFR),
        ids.HCT: ids.impairment_factor(ids.HCT),
        ids.organ_mass(Organ.KIDNEY): ids.impairment_factor(
            ids.organ_mass(Organ.KIDNEY)
        ),
    }
    for o in _flow_organs():
        q = ids.organ_flow(o)
        factors[q] = ids.impairment_factor(q)
    return factors


def unbound_fraction(s: Substance) -> str:
    """Get the formula of the fraction unbound in plasma of a substance under impairment.

    Args:
        s: the substance.

    Returns:
        `fu_ref_<s>` for a substance without binding, else the McNamara 2002
        scaling `1/(1 + f*(1 - fu_ref)/fu_ref)` with the binding factor `f`:
        the factor of the binding protein (1 for an unknown protein) over the
        CKD factor of the binding capacity `fimp_plasma_protein`. The CKD
        factor raises fu of a strongly bound substance by itself
        (`fu ~ fimp_plasma_protein*fu_ref`), fu stays below 1; fu of a weakly
        bound substance rises less (fu_ref 0.5 at a factor 2.5: 0.71, not
        1.25).
    """
    fu_ref = ids.fu_ref(s.sid)
    if s.binding is None:
        return fu_ref
    protein = {
        "albumin": ids.impairment_factor("albumin"),
        "agp": ids.impairment_factor("agp"),
    }.get(s.binding.protein, _ONE)
    binding = f"{protein} / {ids.impairment_factor('plasma_protein')}"
    return f"{_ONE} / ({_ONE} + {binding} * ({_ONE} - {fu_ref}) / {fu_ref})"


def _number(x: float) -> str:
    """A number in a formula."""
    return f"{x!r} dimensionless"


def _interpolation(severity: str, points: Sequence[float]) -> str:
    """The piecewise linear interpolation of points at the severities 0, 1, 2, ...

    Constant below 0 and above the last severity.
    """
    terms = [_number(points[0]), f"{severity} <= 0 dimensionless"]
    for i in range(len(points) - 1):
        a, b = points[i], points[i + 1]
        terms += [
            f"{_number(a)} + {_number(b - a)} * ({severity} - {i} dimensionless)",
            f"{severity} < {i + 1} dimensionless",
        ]
    terms.append(_number(points[-1]))
    return f"piecewise({', '.join(terms)})"


def _notes(text: str, sources: Iterable[str]) -> str:
    """Notes with the citations of the sources."""
    keys = list(dict.fromkeys(sources))
    if not keys:
        return text
    return text + "\n\nSources:\n\n" + "\n".join(f"- {cite(k)}" for k in keys)


def _initial(
    pid: str,
    formula: str,
    name: str,
    text: str,
    sources: Iterable[str] = (),
    unit: UnitDefinition = U.dimensionless,
    *,
    kind: str,
    organ: Organ | None = None,
    protein: str | None = None,
) -> list[Sbase]:
    """A constant parameter with its initial assignment, annotated as `kind`."""
    return [
        annotate(
            Parameter(
                pid,
                None,
                unit,
                constant=True,
                name=name,
                sboTerm=_SBO,
                notes=_notes(f"{text}\n\n`{pid} = {formula}`", sources),
            ),
            kind,
            organ=organ,
            protein=protein,
        ),
        initial_assignment(pid, formula, unit, name),
    ]


def _input(pid: str, value: float, name: str, text: str, kind: str) -> Parameter:
    """A constant parameter of a severity, an initial value, annotated as `kind`."""
    return annotate(
        Parameter(
            pid,
            value,
            U.dimensionless,
            constant=True,
            name=name,
            sboTerm=_SBO,
            notes=text,
        ),
        kind,
    )


def _protein(key: str) -> str | None:
    """The protein of a factor by its key, `None` for a factor of the physiology."""
    return key if key in PROTEINS else None


@dataclass(frozen=True)
class _Entry:
    """A vector of the data with its provenance.

    Attributes:
        values: the values of the classes in order.
        source: the source of the values.
        notes: the notes of the values.
        choice: the classes whose value is a choice of the library, not a
            value of the source.
    """

    values: tuple[float, ...]
    source: str
    notes: str
    choice: tuple[str, ...]

    def provenance(self, classes: Sequence[str]) -> tuple[str, tuple[str, ...]]:
        """The text on the values and the sources of a parameter.

        Args:
            classes: the labels of the classes of the values.

        Returns:
            The values by class with the classes which are model choices, and
            the sources (ModelChoice besides the source for those classes).

        Raises:
            ValueError: if the vector does not have a value per class or a
                choice names no class.
        """
        if len(self.values) != len(classes):
            raise ValueError(
                f"{DATA}.toml: {len(self.values)} values for the classes {classes}."
            )
        unknown = sorted(set(self.choice) - set(classes))
        if unknown:
            raise ValueError(f"{DATA}.toml: choice of unknown classes {unknown}.")
        values = ", ".join(
            f"{c} {v}" for c, v in zip(classes, self.values, strict=True)
        )
        text = f"the values of the classes: {values}"
        sources = (self.source,)
        if self.choice:
            text += f" ({', '.join(self.choice)}: a model choice, not a value of the source)"
            sources += ("ModelChoice",)
        return f"{text}; {self.notes}", sources


def _entry(section: dict[str, Any], key: str) -> _Entry:
    """An entry of the data."""
    entry = section[key]
    return _Entry(
        values=tuple(float(v) for v in entry["value"]),
        source=entry["source"],
        notes=entry.get("notes", ""),
        choice=tuple(entry.get("choice", ())),
    )


_CP_CLASSES = ("A", "B", "C")
_HF_CLASSES = ("I", "II", "III", "IV")


class _Impairment:
    """The SBML objects of the impairment of a definition."""

    def __init__(self, definition: "PBPKDefinition") -> None:
        self.definition = definition
        self.data = _data()
        self.objects: list[Sbase] = []

    def add(self, objects: Iterable[Sbase]) -> None:
        """Add objects."""
        self.objects += objects

    def severities(self) -> None:
        """The severities of cirrhosis and heart failure and the switch of CKD."""
        self.add(
            [
                _input(
                    ids.CP_SEVERITY,
                    0.0,
                    "severity of liver cirrhosis",
                    "0 healthy, 1, 2, 3 the Child-Pugh classes A, B, C; between "
                    "the classes every factor interpolates linearly. An initial "
                    "value, as the covariates.",
                    "child_pugh",
                ),
                _input(
                    ids.HF_SEVERITY,
                    0.0,
                    "severity of heart failure",
                    "0 healthy, 1 to 4 the NYHA classes I to IV; between the "
                    "classes every factor interpolates linearly. An initial "
                    "value, as the covariates.",
                    "heart_failure",
                ),
                _input(
                    ids.CKD_HEPATIC,
                    0.0,
                    "hepatic changes of chronic kidney disease",
                    "Switch (0 off, 1 on) of the reduced activity of the hepatic "
                    "CYP2D6 and OATP1B of chronic kidney disease (by the stage of "
                    f"`{ids.F_RENAL_FUNCTION}`).",
                    "ckd",
                ),
            ]
        )

    # --- Child-Pugh -------------------------------------------------------

    def _cp(self, key: str, entry: _Entry, text: str, healthy: float = 1.0) -> str:
        """A factor of cirrhosis interpolated over the severity (healthy at 0)."""
        pid = ids.condition_factor("cp", key)
        provenance, sources = entry.provenance(_CP_CLASSES)
        self.add(
            _initial(
                pid,
                _interpolation(ids.CP_SEVERITY, [healthy, *entry.values]),
                f"cirrhosis factor {key}",
                f"{text} (CP_severity 1, 2, 3 the classes A, B, C, linear in "
                f"between), {provenance}",
                sources,
                kind="cirrhosis_factor",
                protein=_protein(key),
            )
        )
        return pid

    def child_pugh(self) -> None:
        """The factors of the Child-Pugh classes, the shunt and the functional fraction."""
        cp = self.data["child_pugh"]
        for key in (
            "functional_liver",
            "hepatic_artery_flow",
            "kidney_flow",
            "cardiac_index",
            "gfr",
            "albumin",
            "agp",
            "hematocrit",
        ):
            self._cp(
                key,
                _entry(cp, key),
                f"Factor of liver cirrhosis on the {key.replace('_', ' ')}",
            )
        shunt = self._cp(
            "shunt",
            _entry(cp, "shunt"),
            "Portosystemic shunt of liver cirrhosis",
            healthy=0.0,
        )
        self.add(
            [
                initial_assignment(
                    ids.F_SHUNT, shunt, U.dimensionless, "portosystemic shunt"
                ),
                initial_assignment(
                    ids.F_FUNCTIONAL,
                    ids.condition_factor("cp", "functional_liver"),
                    U.dimensionless,
                    "functional fraction of the liver",
                ),
            ]
        )
        for protein in self.liver_proteins(cp["proteins"]):
            self._cp(
                protein,
                _entry(cp["proteins"], protein),
                f"Factor of liver cirrhosis on the activity of {protein} in the "
                "liver (abundance per g liver)",
            )
        self._closing_flow()

    def _closing_flow(self) -> None:
        """The factor of cirrhosis on the flows of the organs without a factor.

        The organs without a factor of their own (`_closing_organs`) close the
        cardiac output to the cardiac index: their common factor `r` gives
        `sum fQ<o>*m<o> = ci` with the factors `m<o>` of the other organs,
        written as `r = 1 + (ci - 1 - sum fQ<o>*(m<o> - 1))/sum_rest fQ`, which
        is 1 exactly in health (the fractions sum to 1).
        """
        fixed = []
        for o in _flow_organs():
            factor = self._cp_flow(o)
            if factor is not None and factor != _ONE:
                fixed.append(f"{ids.flow_fraction(o)} * ({factor} - {_ONE})")
        rest = " + ".join(ids.flow_fraction(o) for o in self._closing_organs())
        index = ids.condition_factor("cp", "cardiac_index")
        self.add(
            _initial(
                ids.condition_factor("cp", "rest_flow"),
                f"{_ONE} + ({index} - {_ONE} - ({' + '.join(fixed)})) / ({rest})",
                "cirrhosis factor flow of the other organs",
                "Factor of liver cirrhosis on the blood flow of the organs without "
                "a factor of their own (heart, muscle, skin, adipose, bone, "
                "gonads): the cardiac index less the flows of the other organs "
                "(the portal organs and the brain unchanged, the hepatic artery "
                "and the kidney by their factors) over their share of cardiac "
                "output, so that the flows sum to the cardiac output of the "
                "cardiac index.",
                kind="cirrhosis_factor",
            )
        )

    def _cp_flow(self, organ: Organ) -> str | None:
        """The factor of cirrhosis on the flow of an organ, `None` for a closing organ."""
        if organ is Organ.LIVER:
            return ids.condition_factor("cp", "hepatic_artery_flow")
        if organ is Organ.KIDNEY:
            return ids.condition_factor("cp", "kidney_flow")
        if organ in _CP_FIXED:
            return _ONE
        return None

    def _closing_organs(self) -> list[Organ]:
        """The organs whose flows close the cardiac output under cirrhosis."""
        return [o for o in _flow_organs() if self._cp_flow(o) is None]

    def liver_proteins(self, section: dict[str, Any]) -> list[str]:
        """The proteins of a section which a process of the liver uses."""
        used = {
            p.protein
            for p in self.definition.chemistry
            if isinstance(p, Metabolism | Transport)
            and p.organ is Organ.LIVER
            and p.protein is not None
        }
        return sorted(p for p in section if p in used)

    # --- chronic kidney disease -----------------------------------------

    def ckd(self) -> None:
        """The GFR per 1.73 m2 and the factors of the CKD stages.

        The factors are constant per stage, steps at the stage limits of the
        GFR per 1.73 m2 (`eGFR`): a scan or fit of `f_renal_function` across
        a limit sees a discontinuity.
        """
        ckd = self.data["ckd"]
        stages = _entry(ckd, "stage_gfr")
        normal = normal_gfr(self.definition.sex)
        f = ids.F_RENAL_FUNCTION
        self.add(
            _initial(
                ids.EGFR,
                f"{normal!r} l_per_min * {f}",
                "GFR per 1.73 m2 of the renal function [l/min]",
                "The GFR per 1.73 m2 body surface area of the renal function: the "
                f"GFR per 1.73 m2 of the {self.definition.sex.value} reference "
                f"individual ({normal * 1e3:.4g} ml/min, the reference GFR scaled "
                f"to 1.73 m2 with the reference body surface area) times `{f}`; "
                "it gives the stage of chronic kidney disease, the lower limits "
                f"(ml/min/1.73 m2) of normal, mild, moderate and severe "
                f"{list(stages.values)}; {stages.notes}",
                (stages.source, "ICRP89", "DuBois1916"),
                unit=U.l_per_min,
                kind="egfr",
                organ=Organ.KIDNEY,
            )
        )
        limits = [x * 1e-3 for x in stages.values]
        for key in (
            "plasma_protein",
            "hematocrit",
            "kidney_volume",
            "kidney_flow",
            "gastric_emptying_time",
            "si_transit_time",
        ):
            self._stage(
                key, _entry(ckd, key), f"on the {key.replace('_', ' ')}", limits
            )
        for protein in self.liver_proteins(ckd["proteins"]):
            self._stage(
                protein,
                _entry(ckd["proteins"], protein),
                f"on the activity of {protein} in the liver (with "
                f"`{ids.CKD_HEPATIC}` = 1)",
                limits,
            )

    def _stage(
        self, key: str, entry: _Entry, text: str, limits: Sequence[float]
    ) -> None:
        """A factor of chronic kidney disease, constant per stage of the GFR.

        Args:
            key: the key of the factor.
            entry: the values of mild, moderate, severe and end stage CKD.
            text: the description of the factor.
            limits: the lower limits of the GFR per 1.73 m2 (l/min) of normal
                function and of mild, moderate and severe CKD.
        """
        provenance, sources = entry.provenance(_STAGES)
        terms = [_ONE, f"{ids.EGFR} >= {limits[0]!r} l_per_min"]
        for value, limit in zip(entry.values[:-1], limits[1:], strict=True):
            terms += [_number(value), f"{ids.EGFR} >= {limit!r} l_per_min"]
        terms.append(_number(entry.values[-1]))
        self.add(
            _initial(
                ids.condition_factor("ckd", key),
                f"piecewise({', '.join(terms)})",
                f"CKD factor {key}",
                f"Factor of chronic kidney disease {text}, a step per stage "
                f"(normal 1), {provenance}",
                sources,
                kind="ckd_factor",
                protein=_protein(key),
            )
        )

    # --- heart failure ---------------------------------------------------

    def heart_failure(self) -> None:
        """The factors of the NYHA classes.

        The renal blood flow is not monotonic in the class (Leithe 1984: NYHA
        IV above NYHA III), and the GFR follows it.
        """
        hf = self.data["heart_failure"]
        for key in _HF_ORGANS:
            provenance, sources = _entry(hf, key).provenance(_HF_CLASSES)
            self.add(
                _initial(
                    ids.condition_factor("hf", key),
                    _interpolation(ids.HF_SEVERITY, [1.0, *_entry(hf, key).values]),
                    f"heart failure factor {key.replace('_', ' ')}",
                    f"Factor of heart failure on the {key.replace('_', ' ')} of "
                    "the organs "
                    + ", ".join(sorted(o.name.lower() for o in _HF_ORGANS[key]))
                    + " (HF_severity 1 to 4 the NYHA classes I to IV, linear in "
                    f"between), {provenance}",
                    sources,
                    kind="heart_failure_factor",
                )
            )

    # --- the factors of all conditions ----------------------------------

    def combined(self) -> None:
        """The factors `fimp_<id>` of all conditions on the quantities."""

        def product(*factors: str | None) -> str:
            present = [f for f in factors if f is not None and f != _ONE]
            return " * ".join(present) if present else _ONE

        text = "The factor of all impairments (cirrhosis, CKD, heart failure) on"
        rest = ids.condition_factor("cp", "rest_flow")
        flows: list[str] = []
        closing: dict[str, list[str]] = {}
        for o in _flow_organs():
            hf_flow = next(
                (
                    ids.condition_factor("hf", k)
                    for k, v in _HF_ORGANS.items()
                    if o in v
                ),
                None,
            )
            ckd_flow = (
                ids.condition_factor("ckd", "kidney_flow")
                if o is Organ.KIDNEY
                else None
            )
            q = ids.organ_flow(o)
            cp_flow = self._cp_flow(o)
            self.add(
                _initial(
                    ids.impairment_factor(q),
                    product(cp_flow or rest, ckd_flow, hf_flow),
                    f"impairment factor {q}",
                    f"{text} the blood flow `{q}`; the cardiac output is the "
                    f"sum of the flows (`{ids.impairment_factor(ids.CO)}`).",
                    kind="impairment_factor",
                    organ=o,
                )
            )
            if cp_flow is None:
                closing.setdefault(product(ckd_flow, hf_flow), []).append(
                    ids.flow_fraction(o)
                )
            else:
                flows.append(
                    f"{ids.flow_fraction(o)} * ({ids.impairment_factor(q)} - {_ONE})"
                )
        # the organs which close the cardiac output under cirrhosis share the
        # factor `rest`, which enters the sum once (a compact initial
        # assignment); the sum is written as 1 + sum fQ*(m - 1), 1 exactly in
        # health
        flows += [
            f"({' + '.join(fractions)}) * ({rest} * {factor} - {_ONE})"
            if factor != _ONE
            else f"({' + '.join(fractions)}) * ({rest} - {_ONE})"
            for factor, fractions in closing.items()
        ]
        self.add(
            _initial(
                ids.impairment_factor(ids.CO),
                f"{_ONE} + " + " + ".join(flows),
                "impairment factor cardiac output",
                f"{text} the cardiac output: the sum of the blood flows of the "
                "organs with their factors as fractions of the healthy cardiac "
                "output (cirrhosis: the cardiac index; heart failure: the sum "
                "of the reduced flows, Gu 2025).",
                ("Gu2025",),
                kind="impairment_factor",
            )
        )
        quantities = {
            ids.GFR: (
                product(
                    ids.condition_factor("cp", "gfr"),
                    ids.condition_factor("hf", "kidney_flow"),
                ),
                "the GFR (heart failure: the GFR follows the renal blood flow; "
                f"chronic kidney disease is `{ids.F_RENAL_FUNCTION}`)",
            ),
            ids.HCT: (
                product(
                    ids.condition_factor("cp", "hematocrit"),
                    ids.condition_factor("ckd", "hematocrit"),
                ),
                "the hematocrit of the sex of the build; the anaemia of "
                "cirrhosis and of chronic kidney disease combine "
                "multiplicatively",
            ),
            ids.organ_mass(Organ.KIDNEY): (
                ids.condition_factor("ckd", "kidney_volume"),
                "the kidney mass (adipose, the remainder, takes the difference)",
            ),
            "albumin": (ids.condition_factor("cp", "albumin"), "the plasma albumin"),
            "agp": (
                ids.condition_factor("cp", "agp"),
                "the plasma alpha-1-acid glycoprotein",
            ),
            "plasma_protein": (
                ids.condition_factor("ckd", "plasma_protein"),
                "the binding capacity of the plasma proteins: it divides the "
                "binding factor of fu, so fu of a strongly bound drug rises by it",
            ),
            "t_stomach": (
                ids.condition_factor("ckd", "gastric_emptying_time"),
                "the gastric emptying time",
            ),
            "t_si": (
                ids.condition_factor("ckd", "si_transit_time"),
                "the transit times of the small intestine",
            ),
        }
        for key, (formula, what) in quantities.items():
            self.add(
                _initial(
                    ids.impairment_factor(key),
                    formula,
                    f"impairment factor {key}",
                    f"{text} {what}.",
                    kind="impairment_factor",
                )
            )
        cp_proteins = set(self.liver_proteins(self.data["child_pugh"]["proteins"]))
        ckd_proteins = set(self.liver_proteins(self.data["ckd"]["proteins"]))
        for protein in sorted(cp_proteins | ckd_proteins):
            hepatic = (
                f"({_ONE} + {ids.CKD_HEPATIC} * ({ids.condition_factor('ckd', protein)} - {_ONE}))"
                if protein in ckd_proteins
                else None
            )
            self.add(
                _initial(
                    ids.impairment_activity(protein, Organ.LIVER),
                    product(
                        ids.condition_factor("cp", protein)
                        if protein in cp_proteins
                        else None,
                        hepatic,
                    ),
                    f"impairment factor activity {protein} liver",
                    f"{text} the activity of {protein} in the liver, the input "
                    f"`{ids.impairment_activity(protein)}` of the liver; its rates "
                    f"run with `{ids.activity(protein)}` times it.",
                    kind="impairment_factor",
                    organ=Organ.LIVER,
                    protein=protein,
                )
            )


def impairment_objects(definition: "PBPKDefinition") -> list[Sbase]:
    """Get the parameters of the impaired organ function of the top model.

    Args:
        definition: the definition (the proteins of the liver).

    Returns:
        The severities, the factors of every condition and of all conditions
        (initial assignments) and the initial assignments of `f_shunt` and
        `f_functional`; the physiology reads the factors through
        `physiology_factors`.
    """
    impairment = _Impairment(definition)
    impairment.severities()
    impairment.child_pugh()
    impairment.ckd()
    impairment.heart_failure()
    impairment.combined()
    return impairment.objects
