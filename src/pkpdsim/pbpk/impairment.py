"""The impaired organ function: parameter vectors on the physiology.

Impairment is never a structurally different model: three severities of the
top model set factors on the quantities of the healthy physiology, the
liver and the binding to plasma proteins (the table `impairment` of the
parameter database `pkpdsim.database`, every number with its source; a
value which no publication read gives is a model choice):

- `CP_severity`: liver cirrhosis, 0 healthy, 1, 2, 3 the Child-Pugh classes
  A, B, C; a severity between the classes interpolates every factor linearly
  (scans and fits use intermediate values). It sets the functional fraction
  `f_functional` and the portosystemic shunt `f_shunt` of the liver, the
  hepatic arterial and the renal blood flow, the cardiac output (cardiac
  index; the organs without a factor of their own close the cardiac output to
  it), the GFR, albumin and alpha-1-acid glycoprotein (and so fu) and the
  activity of the enzymes and transporters of the liver: the consensus of the
  studies of the database on the abundance of the protein per amount of
  liver tissue in each class over the controls (`pkpdsim.database.integration`,
  the mRNA expression only for a protein without an abundance), every class
  with its interval and studies in the notes of its factor.
- `f_renal_function`: chronic kidney disease by the GFR relative to healthy
  (the factor of the physiology); its GFR per 1.73 m2 (`eGFR`) gives the
  stage (normal, mild, moderate, severe, end stage), which sets the binding
  capacity of the plasma proteins (and so fu). The factors are piecewise
  linear in the eGFR through anchor points at the stage midpoints (`[ckd]
  stage_anchor` of the data, a model choice: continuous covariate for scans
  and fits; stage midpoints): 1 from 90 ml/min/1.73 m2, the mild value at 75,
  the moderate value at 45, the severe value at 22.5 and the end stage value
  at 7.5 and below, linear in between (from 90 to 75 from 1 to the mild
  value), so a scan or fit of `f_renal_function` sees no step. The binding
  to albumin and alpha-1-acid glycoprotein follows the CKD profile of its
  protein (`[ckd.binding]` of the data): albumin falls in concentration and
  in the binding function of its molecule, AAG rises in concentration; the
  binding to an unknown protein follows the binding capacity of the plasma
  proteins (`plasma_protein`). The switch `CKD_hepatic` (0 off, 1 on) adds the reduced
  activity of hepatic CYP2D6 and OATP1B. The tubular capacity of the kidney
  (every metabolism and transport of the kidney) follows `f_renal_function`
  continuously, `1 - f_tubular + f_tubular*f_renal_function` with the share
  `f_tubular` 1 of the intact nephron hypothesis (`[ckd] tubular_capacity`
  of the data, a model choice; the kidney owns the factor,
  `pkpdsim.pbpk.organ`).
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
binding factor `f`, the factor of the binding protein in cirrhosis
(`fimp_albumin`, `fimp_agp`, 1 for an unknown protein) times the CKD factor of
the binding: the CKD profile of the protein (`fimp_albumin_ckd`,
`fimp_agp_ckd`, the product of its rows of `[ckd.binding]`), for an unknown
protein one over the factor of the binding capacity `fimp_plasma_protein`
(fu of a strongly bound substance rises by it). fu stays below 1. A binding
with `ckd=False` (measured unchanged in renal impairment) leaves out the CKD
factor, so `f` is the factor of the binding protein alone. The parameters of
the CKD profile of a protein are created only when a substance binds to it
with the CKD factor. Enzymes and transporters: a rate of a
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

import itertools
import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from functools import cache
from typing import TYPE_CHECKING, Any

from sbmlutils.factory import Parameter, Sbase, UnitDefinition
from sbmlutils.metadata import SBO

from pkpdsim.database import integration, studies
from pkpdsim.pbpk import ids
from pkpdsim.pbpk.adipose import (
    BINDING_FACTORS,
    CARDIAC_OUTPUT,
    FLOW_FACTORS,
    ObesityFactor,
    obesity_factor,
    obesity_factors,
    obesity_protein_factors,
)
from pkpdsim.pbpk.age import (
    AGE_REFERENCE,
    AgeCurve,
    age_factor,
    factor_curves,
    protein_factor,
    protein_factors,
)
from pkpdsim.pbpk.annotation import annotate, initial_assignment
from pkpdsim.pbpk.chemistry import Metabolism, Transport
from pkpdsim.pbpk.circadian import circadian_factor, circadian_keys, protein_rhythms
from pkpdsim.pbpk.data import cite, load
from pkpdsim.pbpk.organs import PORTAL_ORGANS, Organ
from pkpdsim.pbpk.physiology import normal_gfr
from pkpdsim.pbpk.pregnancy import (
    ORGAN_FLOW,
    PLACENTA_SHARE,
    pregnancy_factor,
)
from pkpdsim.pbpk.proteins import PROTEINS
from pkpdsim.pbpk.sex_differences import (
    ALBUMIN,
    protein_sex_factor,
    protein_sex_factors,
    sex_curves,
    sex_factor,
)
from pkpdsim.pbpk.substance import Substance
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import number
from pkpdsim.sex import Sex
from pkpdsim.species import Species

if TYPE_CHECKING:
    from pkpdsim.pbpk.definition import PBPKDefinition

__all__ = [
    "GASTRIC_EMPTYING",
    "MELD_HEALTHY",
    "SI_TRANSIT",
    "ImpairmentFactors",
    "MeldAnchor",
    "child_pugh_protein",
    "child_pugh_proteins",
    "class_factors",
    "cp_severity_of_meld",
    "impairment_objects",
    "meld_anchors",
    "meld_of_cp_severity",
    "organ_activity",
    "physiology_factors",
    "unbound_fraction",
    "unbound_fraction_notes",
]

DATA = "impairment"
"""The table of the impairment in the database."""

_CP_CLASSES = ("A", "B", "C")
"""The Child-Pugh classes, the severities 1, 2, 3."""

_CP_KINDS = ("protein", "mrna")
"""The kinds of measurement of the factors of cirrhosis on the proteins, protein first."""

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


def _binding_keys(protein: str) -> list[str]:
    """The keys of the rows of the CKD profile of the binding to a protein.

    Args:
        protein: the binding protein, `albumin`, `agp` or `unknown`.

    Returns:
        The keys of `[ckd.binding]` of the data whose protein it is, in their
        order; empty for a protein without a profile (an unknown protein).
    """
    rows = _data()["ckd"].get("binding", {})
    return [key for key, row in rows.items() if row["protein"] == protein]


def _ckd_binding(protein: str) -> str:
    """The id of the CKD factor of the binding to a protein, `fimp_albumin_ckd`."""
    return ids.impairment_factor(f"{protein}_ckd")


@cache
def _pregnancy_proteins() -> frozenset[tuple[Organ, str]]:
    """The organs and proteins with a factor of the pregnancy."""
    from pkpdsim.pbpk.pregnancy import protein_factors

    return frozenset(protein_factors())


@cache
def _obesity_proteins() -> frozenset[tuple[Organ, str]]:
    """The organs and proteins with a factor of obesity (alike in both sexes)."""
    return frozenset(obesity_protein_factors(Sex.MALE))


def organ_activity(organ: Organ, protein: str) -> str | None:
    """Get the local input of the individual on a protein in an organ.

    Args:
        organ: the organ.
        protein: the protein of a process of the organ.

    Returns:
        `fimp_<protein>`, replaced by `fimp_<organ>_<protein>` of the body,
        for a protein with data of impairment in the liver or with a factor
        of the age or the sex in the organ (`pkpdsim.pbpk.age`,
        `pkpdsim.pbpk.sex_differences`); `None` otherwise.
    """
    if (organ is Organ.LIVER and protein in _proteins()) or (
        organ,
        protein,
    ) in (
        protein_factors()
        | protein_sex_factors()
        | set(protein_rhythms())
        | _pregnancy_proteins()
        | _obesity_proteins()
    ):
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
        sum of the impaired flows), the flows, the GFR (at renal function 1,
        `GFR_base`), the hematocrit and the
        kidney mass.
    """
    co = ids.impairment_factor(ids.CO)
    factors = {
        ids.CO: co,
        ids.GFR_BASE: ids.impairment_factor(ids.GFR),
        ids.HCT: ids.impairment_factor(ids.HCT),
        ids.organ_mass(Organ.KIDNEY): ids.impairment_factor(
            ids.organ_mass(Organ.KIDNEY)
        ),
    }
    for o in _flow_organs():
        q = ids.organ_flow(o)
        factors[q] = ids.impairment_factor(q)
    return factors


def unbound_fraction_notes(s: Substance) -> str:
    """Get the notes of the fraction unbound in plasma of a substance.

    The notes say what `unbound_fraction` gives: the reference fraction scaled
    with the binding protein under impairment, scaled with the binding
    capacity of chronic kidney disease alone (an unknown protein with the CKD
    factor), the reference fraction itself where no factor applies (an
    unknown protein without the CKD factor), or no binding.

    Args:
        s: the substance.

    Returns:
        The notes of `fu_<s>`.
    """
    if s.binding is None:
        return (
            "Fraction unbound in plasma of the individual; the substance does "
            "not bind to plasma proteins."
        )
    if s.binding.protein == "unknown" and not s.binding.ckd:
        return (
            "Fraction unbound in plasma of the individual: equal to the reference "
            "fraction, no impairment factor applies (the binding protein is "
            "unknown and the binding has no CKD factor)."
        )
    if s.binding.protein == "unknown":
        return (
            "Fraction unbound in plasma of the individual: the reference fraction "
            "scaled with the binding capacity of the plasma proteins in chronic "
            "kidney disease, no factor of a binding protein applies (the binding "
            "protein is unknown; McNamara 2002, `pkpdsim.pbpk.impairment`)."
        )
    if s.binding.ckd and _binding_keys(s.binding.protein):
        return (
            "Fraction unbound in plasma of the individual: the reference fraction "
            "scaled with the binding protein under impairment, its factor of "
            "cirrhosis times its profile of chronic kidney disease "
            f"`{_ckd_binding(s.binding.protein)}` (McNamara 2002, "
            "`pkpdsim.pbpk.impairment`)."
        )
    return (
        "Fraction unbound in plasma of the individual: the reference fraction "
        "scaled with the binding protein under impairment (McNamara 2002, "
        "`pkpdsim.pbpk.impairment`)."
    )


def unbound_fraction(s: Substance) -> str:
    """Get the formula of the fraction unbound in plasma of a substance under impairment.

    Args:
        s: the substance.

    Returns:
        `fu_ref_<s>` for a substance without binding, else the McNamara 2002
        scaling `1/(1 + f*(1 - fu_ref)/fu_ref)` with the binding factor `f`:
        the factor of the binding protein in cirrhosis times the CKD profile
        of the protein (`fimp_albumin_ckd`, `fimp_agp_ckd`); an unknown
        protein has no factor of cirrhosis and its CKD factor is one over the
        factor of the binding capacity `fimp_plasma_protein`, which raises fu
        of a strongly bound substance by itself
        (`fu ~ fimp_plasma_protein*fu_ref`), fu stays below 1; fu of a weakly
        bound substance rises less (fu_ref 0.5 at a factor 2.5: 0.71, not
        1.25). A binding with `ckd=False` has the factor of the binding
        protein alone; with an unknown protein as well no factor applies and
        the result is `fu_ref_<s>` itself.
    """
    fu_ref = ids.fu_ref(s.sid)
    if s.binding is None:
        return fu_ref
    if s.binding.protein == "unknown" and not s.binding.ckd:
        return fu_ref
    protein = {
        "albumin": ids.impairment_factor("albumin"),
        "agp": ids.impairment_factor("agp"),
    }.get(s.binding.protein, _ONE)
    binding = protein
    if s.binding.ckd and _binding_keys(s.binding.protein):
        binding = f"{protein} * {_ckd_binding(s.binding.protein)}"
    elif s.binding.ckd:
        binding = f"{protein} / {ids.impairment_factor('plasma_protein')}"
    return f"{_ONE} / ({_ONE} + {binding} * ({_ONE} - {fu_ref}) / {fu_ref})"


def ckd_unbound_fractions(s: Substance) -> dict[str, float]:
    """Get the fraction unbound of a substance at the CKD stages.

    The fu of `unbound_fraction` without cirrhosis at the values of the
    stages mild, moderate, severe and end stage of the CKD factor of its
    binding (the profile of its protein, `[ckd.binding]`, or one over the
    factor of the binding capacity `plasma_protein` of an unknown protein):
    the factors are linear in the eGFR between the anchors of the stages and
    constant below the end stage, so the fu of every renal function lies
    between these values (and the reference fu of renal function 1).

    Args:
        s: the substance.

    Returns:
        The fu by stage; empty for a substance whose fu CKD does not change
        (no binding, `ckd=False`) or whose fu is a `Formula`.
    """
    if s.binding is None or not s.binding.ckd:
        return {}
    fu = number(s.binding.fu)
    if fu is None:
        return {}
    keys = _binding_keys(s.binding.protein)
    if keys:
        rows = [class_factors("ckd.binding", key).values for key in keys]
        factors = [math.prod(row[k] for row in rows) for k in range(len(_STAGES))]
    else:
        factors = [1.0 / v for v in class_factors("ckd", "plasma_protein").values]
    return {
        stage: 1.0 / (1.0 + f * (1.0 - fu) / fu)
        for stage, f in zip(_STAGES, factors, strict=True)
    }


def _number(x: float) -> str:
    """A number in a formula."""
    return f"{x!r} dimensionless"


def _l_per_min(ml_per_min: float) -> str:
    """A flow in a formula, from ml/min (the eGFR is in l/min)."""
    return f"{ml_per_min / 1000.0!r} l_per_min"


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
class ImpairmentFactors:
    """The factors of the classes of an impairment with their provenance.

    Attributes:
        values: the values of the classes in order.
        source: the source of the values.
        notes: the notes of the values.
        choice: the classes whose value is a choice of the library, not a
            value of the source.
    """

    values: tuple[float, ...]
    sources: tuple[str, ...]
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
        sources = self.sources
        if self.choice:
            text += f" ({', '.join(self.choice)}: a model choice, not a value of the source)"
            sources += ("ModelChoice",)
        return f"{text}; {self.notes}", sources


def _entry(section: dict[str, Any], key: str) -> ImpairmentFactors:
    """An entry of the data."""
    entry = section[key]
    return ImpairmentFactors(
        values=tuple(float(v) for v in entry["value"]),
        sources=(entry["source"],),
        notes=entry.get("notes", ""),
        choice=tuple(entry.get("choice", ())),
    )


_HF_CLASSES = ("I", "II", "III", "IV")


@cache
def child_pugh_proteins() -> frozenset[str]:
    """The proteins with a consensus of cirrhosis of a Child-Pugh class (`pkpdsim.database`)."""
    return frozenset(
        entity
        for kind in _CP_KINDS
        for entity in integration.entities(kind=kind, tissue="liver")
        if any(_cp_consensus(entity, c) is not None for c in _CP_CLASSES)
    )


def _cp_consensus(protein: str, cls: str) -> integration.Consensus | None:
    """The consensus of the database on a protein in a Child-Pugh class, protein before mRNA."""
    return integration.factor(
        protein,
        {"condition": "cirrhosis", "child_pugh": cls},
        tissue="liver",
        kinds=_CP_KINDS,
    )


def _consensus_factors(
    found: dict[str, integration.Consensus | None],
    classes: Sequence[str],
    what: str,
    notes: str = "",
) -> ImpairmentFactors:
    """The factors of the classes of a condition of the consensus of every class.

    A class without a consensus is the linear interpolation between its
    neighbours (healthy 1 before the first class) and the last class with a
    consensus beyond it, a model choice.

    Args:
        found: the consensus by class, `None` for a class without a study.
        classes: the classes in the order of their severity.
        what: the quantity of the ratios, for the notes.
        notes: further notes before the consensus of the classes.

    Returns:
        The factors with the studies of every class.

    Raises:
        ValueError: if no class has a consensus.
    """
    # four significant digits, far within the confidence of the consensus
    known = [(0, 1.0)] + [
        (i + 1, float(f"{f.value:.4g}"))
        for i, c in enumerate(classes)
        if (f := found[c]) is not None
    ]
    if len(known) == 1:
        raise ValueError(f"No consensus of the classes {list(classes)} for {what}.")
    values = []
    for i in range(1, len(classes) + 1):
        below = max((k for k in known if k[0] <= i), key=lambda k: k[0])
        above = [k for k in known if k[0] >= i]
        if below[0] == i or not above:
            values.append(below[1])
            continue
        hi = min(above, key=lambda k: k[0])
        value = below[1] + (hi[1] - below[1]) * (i - below[0]) / (hi[0] - below[0])
        values.append(float(f"{value:.4g}"))
    texts = []
    for c in classes:
        f = found[c]
        texts.append(f"{c} {f.text()}" if f is not None else f"{c} no study")
    choice = tuple(c for c in classes if found[c] is None)
    sources = tuple(
        dict.fromkeys(s for f in found.values() if f is not None for s in f.sources)
    )
    text = (
        "the consensus of the studies of the parameter database "
        f"(`pkpdsim.database.integration`): the ratio of {what} of the class to "
        "the controls of each study on the log scale, pooled by a random-effects "
        "meta-analysis; " + "; ".join(texts)
    )
    if notes:
        text = f"{notes}; {text}"
    if choice:
        text += (
            "; a class without a study is the linear interpolation between its "
            "neighbours (healthy 1) or the last class with a study, a model choice"
        )
    return ImpairmentFactors(
        values=tuple(values), sources=sources, notes=text, choice=choice
    )


def child_pugh_protein(protein: str) -> ImpairmentFactors:
    """The factors of a protein of the Child-Pugh classes from the consensus of the database.

    A class without a consensus is the linear interpolation between its
    neighbours (healthy 1 at severity 0) and the last class with a consensus
    beyond it, a model choice.

    Args:
        protein: the gene symbol of the protein.

    Returns:
        The factors of the classes A, B, C with the studies of every class.

    Raises:
        ValueError: if no class of the protein has a consensus.
    """
    return _consensus_factors(
        {c: _cp_consensus(protein, c) for c in _CP_CLASSES},
        _CP_CLASSES,
        f"the abundance per amount of liver tissue of {protein}",
    )


_SECTIONS: dict[str, tuple[tuple[str, ...], str, str]] = {
    "child_pugh": (_CP_CLASSES, "cirrhosis", "child_pugh"),
    "ckd": (_STAGES, "ckd", "ckd_stage"),
    "heart_failure": (_HF_CLASSES, "heart_failure", "nyha"),
}
"""The conditions of the table: their classes and the condition and the field
of the class of the groups of the studies."""


def class_factors(section: str, key: str) -> ImpairmentFactors:
    """The factors of the classes of a condition on a quantity.

    An entry of the table `impairment` states its values (`value`), or names
    the quantity of the studies of the database whose consensus they are
    (`consensus`, a key of `pkpdsim.database.studies.QUANTITIES`): the ratio of
    the class to the controls of each study, pooled over the studies, a class
    without a study interpolated as for the proteins (`child_pugh_protein`).
    The classes are the Child-Pugh classes A, B, C (`child_pugh`), the CKD
    stages mild, moderate, severe and end stage (`ckd`, `ckd.binding`) and
    the NYHA classes I to IV (`heart_failure`).

    Args:
        section: the section of the table, e.g. `child_pugh` or `ckd.binding`.
        key: the key of the factor, e.g. `albumin`.

    Returns:
        The factors of the classes with their provenance.

    Raises:
        KeyError: for a section or a key which is no factor of the table.
    """
    path = section.split(".")
    if path[0] not in _SECTIONS:
        raise KeyError(f"No section '{section}' of the impairment, {list(_SECTIONS)}.")
    table: dict[str, Any] = _data()
    for part in path:
        table = table[part]
    entry = table[key]
    quantity = entry.get("consensus")
    if quantity is None:
        return _entry(table, key)
    classes, condition, field = _SECTIONS[path[0]]
    return _consensus_factors(
        {
            c: integration.consensus(
                quantity,
                {"condition": condition, field: c},
                kind="physiology",
                tissue="body",
            )
            for c in classes
        },
        classes,
        f"the {studies.QUANTITIES[quantity].description} ({quantity})",
        entry.get("notes", ""),
    )


MELD_HEALTHY = 6.43
"""The MELD score of `CP_severity` 0: the smallest score of the formula
(`3.78 ln(bilirubin) + 11.2 ln(INR) + 9.57 ln(creatinine) + 6.43`, every
laboratory value bounded below by 1 mg/dL or 1), which subjects with a normal
bilirubin, INR and creatinine have (Yang 2020: 6.43 in every quartile of the
patients with a normal liver function)."""


@dataclass(frozen=True)
class MeldAnchor:
    """The MELD score of a severity of cirrhosis.

    Attributes:
        severity: the severity `CP_severity`, 0 healthy, 1, 2, 3 the classes.
        meld: the MELD score of the severity.
        level: the consensus of the studies on the MELD of the class, `None`
            for the healthy severity (`MELD_HEALTHY`).
    """

    severity: float
    meld: float
    level: integration.Level | None


@cache
def meld_anchors() -> tuple[MeldAnchor, ...]:
    """The MELD scores of the severities of cirrhosis, the anchors of the mapping.

    Severity 0 is `MELD_HEALTHY`; the classes A, B, C (severity 1, 2, 3) have
    the consensus of the studies of the database on the MELD of the patients
    of the class (`pkpdsim.database.integration.level`, the mean or median of
    each study pooled by a random-effects meta-analysis), to four significant
    digits.

    Returns:
        The anchors in the order of the severities, the MELD increasing.

    Raises:
        ValueError: if a class has no MELD in the database or the MELD of the
            anchors does not increase with the severity.
    """
    anchors = [MeldAnchor(0.0, MELD_HEALTHY, None)]
    for i, c in enumerate(_CP_CLASSES, start=1):
        found = integration.level("meld", {"condition": "cirrhosis", "child_pugh": c})
        if found is None:
            raise ValueError(f"No MELD of the Child-Pugh class {c} in the database.")
        anchors.append(MeldAnchor(float(i), float(f"{found.value:.4g}"), found))
    melds = [a.meld for a in anchors]
    if any(b <= a for a, b in itertools.pairwise(melds)):
        raise ValueError(f"The MELD of the severities does not increase: {melds}.")
    return tuple(anchors)


def cp_severity_of_meld(meld: float) -> float:
    """Get the severity of cirrhosis `CP_severity` of a MELD score.

    The severity is piecewise linear in the MELD through the anchors of
    `meld_anchors` (a model choice: the MELD of a class is the consensus of
    the patients of the class, between the classes the factors and the MELD
    are both linear in the severity). A MELD at or below `MELD_HEALTHY` is
    healthy (0), one at or above the MELD of class C is class C (3), since the
    library has no factors beyond the class.

    Args:
        meld: the MELD score.

    Returns:
        The severity, from 0 to 3.

    Raises:
        ValueError: for a MELD which is not finite.
    """
    if not math.isfinite(meld):
        raise ValueError(f"The MELD {meld} is not finite.")
    anchors = meld_anchors()
    if meld <= anchors[0].meld:
        return 0.0
    for lo, hi in itertools.pairwise(anchors):
        if meld <= hi.meld:
            fraction = (meld - lo.meld) / (hi.meld - lo.meld)
            return lo.severity + fraction * (hi.severity - lo.severity)
    return anchors[-1].severity


def meld_of_cp_severity(severity: float) -> float:
    """Get the MELD score of a severity of cirrhosis, the inverse of `cp_severity_of_meld`.

    Args:
        severity: the severity `CP_severity`, from 0 to 3.

    Returns:
        The MELD, piecewise linear in the severity through `meld_anchors`.

    Raises:
        ValueError: for a severity outside of 0 to 3.
    """
    anchors = meld_anchors()
    if not anchors[0].severity <= severity <= anchors[-1].severity:
        raise ValueError(
            f"The severity {severity} is outside of {anchors[0].severity} to "
            f"{anchors[-1].severity}."
        )
    for lo, hi in itertools.pairwise(anchors):
        if severity <= hi.severity:
            fraction = (severity - lo.severity) / (hi.severity - lo.severity)
            return lo.meld + fraction * (hi.meld - lo.meld)
    return anchors[-1].meld


def _proteins() -> frozenset[str]:
    """The proteins with a factor of impairment in the liver."""
    return child_pugh_proteins() | frozenset(_data()["ckd"]["proteins"])


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

    def _cp(
        self, key: str, entry: ImpairmentFactors, text: str, healthy: float = 1.0
    ) -> str:
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
                class_factors("child_pugh", key),
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
        for protein in self.liver_proteins(child_pugh_proteins()):
            self._cp(
                protein,
                child_pugh_protein(protein),
                f"Factor of liver cirrhosis on the activity of {protein} in the "
                "liver (abundance per amount of liver tissue)",
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

    def liver_proteins(self, section: Iterable[str]) -> list[str]:
        """The proteins of a set which a process of the liver uses."""
        used = {
            p.protein
            for p in self.definition.chemistry
            if isinstance(p, Metabolism | Transport)
            and p.organ is Organ.LIVER
            and p.protein is not None
        }
        return sorted(p for p in section if p in used)

    def individual_proteins(self) -> list[tuple[Organ, str]]:
        """The organs and proteins of the processes with a factor of the individual."""
        used = {
            (p.organ, p.protein)
            for p in self.definition.chemistry
            if isinstance(p, Metabolism | Transport) and p.protein is not None
        }
        return sorted(
            ((o, p) for o, p in used if organ_activity(o, p) is not None),
            key=lambda x: (x[0].value, x[1]),
        )

    def _sex_formula(self, factor: AgeCurve) -> str:
        """The factor of the sex: its curve of the age, or its value at 35 years.

        A build without the physiology over age (`PBPKDefinition.age`) has no
        `AGE`, so the factor is the constant of the adult.
        """
        if self.definition.age:
            return factor.curve.formula(ids.AGE)
        return _number(factor.curve.value(AGE_REFERENCE))

    def sex(self) -> None:
        """The factors of the sex on the albumin and on the proteins of the organs.

        A build over age has them as curves of `AGE`, a build without as the
        constants of the adult (`_sex_formula`).
        """
        sex = self.definition.sex
        albumin = sex_curves(sex).get(ALBUMIN)
        if albumin is not None:
            self.add(
                _initial(
                    sex_factor(ALBUMIN),
                    self._sex_formula(albumin),
                    "sex factor albumin",
                    f"The plasma albumin of the {sex.value} build relative to the "
                    "men of the same age (`pkpdsim.pbpk.sex_differences`): "
                    f"{albumin.notes}",
                    albumin.sources,
                    kind="sex_factor",
                )
            )
        for organ, protein in self.individual_proteins():
            factor = protein_sex_factor(sex, organ, protein)
            if factor is None:
                continue
            where = f"{protein} in the {organ.name.lower().replace('_', ' ')}"
            self.add(
                _initial(
                    sex_factor(f"{organ.value}_{protein.lower()}"),
                    self._sex_formula(factor),
                    f"sex factor {protein} {organ.name.lower()}",
                    f"The activity of {where} of the {sex.value} build relative "
                    "to the men of the same age, per amount of tissue "
                    f"(`pkpdsim.pbpk.sex_differences`): {factor.notes}",
                    factor.sources,
                    kind="sex_factor",
                    organ=organ,
                    protein=protein,
                )
            )

    def obesity(self) -> None:
        """The factors `fobes_<key>` of overweight and obesity over the body mass index.

        Each is 1 up to a body mass index of 25 kg/m2 and is the residual of
        the obese repository of Berton 2022 over what the scaling of the
        library already does (`pkpdsim.pbpk.adipose`); the factor of the
        cardiac output carries every organ flow, so the cardiac output of the
        model stays the sum of its flows.
        """
        sex = self.definition.sex
        factors = obesity_factors(sex)

        def add(
            pid: str,
            factor: ObesityFactor,
            what: str,
            organ: Organ | None = None,
            protein: str | None = None,
        ) -> None:
            self.add(
                _initial(
                    obesity_factor(pid),
                    factor.formula(),
                    f"obesity factor {pid}",
                    f"The factor of overweight and obesity on {what}, a function "
                    f"of the body mass index of the {sex.value} individual which "
                    f"is 1 up to 25 kg/m2 (`pkpdsim.pbpk.adipose`): "
                    f"{factor.notes}",
                    factor.sources,
                    kind="obesity_factor",
                    organ=organ,
                    protein=protein,
                )
            )

        add(ids.CO, factors[CARDIAC_OUTPUT], "the cardiac output")
        for key, organs in FLOW_FACTORS.items():
            factor = factors[key]
            for organ in organs:
                add(
                    ids.organ_flow(organ),
                    factor,
                    f"the blood flow fraction of the "
                    f"{organ.name.lower().replace('_', ' ')}",
                    organ=organ,
                )
        for key in BINDING_FACTORS:
            add(key, factors[key], f"the plasma {key}")
        for organ, protein in self.individual_proteins():
            factor = obesity_protein_factors(sex).get((organ, protein))
            if factor is None:
                continue
            add(
                f"{organ.value}_{protein.lower()}",
                factor,
                f"the activity of {protein} in the "
                f"{organ.name.lower().replace('_', ' ')} per amount of tissue",
                organ=organ,
                protein=protein,
            )

    def age(self) -> None:
        """The factors of the age on the binding proteins and on the proteins."""
        sex = self.definition.sex
        curves = factor_curves(sex)
        for key, name in (("albumin", "albumin"), ("agp", "alpha-1-acid glycoprotein")):
            c = curves[key]
            self.add(
                _initial(
                    age_factor(key),
                    c.curve.formula(ids.AGE),
                    f"age factor {name}",
                    f"The factor of the age on the plasma {name} of the "
                    f"{sex.value} individual, 1 at {AGE_REFERENCE:g} years "
                    f"(`pkpdsim.pbpk.age`): {c.notes}",
                    c.sources,
                    kind="age_factor",
                )
            )
        for organ, protein in self.individual_proteins():
            key = f"{organ.value}_{protein.lower()}"
            where = f"{protein} in the {organ.name.lower().replace('_', ' ')}"
            ontogeny = protein_factor(organ, protein)
            if ontogeny is not None:
                self.add(
                    _initial(
                        age_factor(key),
                        ontogeny.curve.formula(ids.AGE),
                        f"age factor {protein} {organ.name.lower()}",
                        f"The ontogeny of {where} per amount of tissue, 1 at "
                        f"{AGE_REFERENCE:g} years (`pkpdsim.pbpk.age`): "
                        f"{ontogeny.notes}",
                        ontogeny.sources,
                        kind="age_factor",
                        organ=organ,
                        protein=protein,
                    )
                )

    def binding_proteins(self) -> list[str]:
        """The proteins with a CKD profile which a substance binds to with the CKD factor."""
        return sorted(
            {
                s.binding.protein
                for s in self.definition.substances
                if s.binding is not None
                and s.binding.ckd
                and _binding_keys(s.binding.protein)
            }
        )

    # --- chronic kidney disease -----------------------------------------

    def ckd(self) -> None:
        """The GFR per 1.73 m2 and the factors of the CKD stages.

        The factors are piecewise linear in the GFR per 1.73 m2 (`eGFR`)
        through the anchors at the stage midpoints (`stage_anchor`), so a
        scan or fit of `f_renal_function` sees no step at a stage limit.
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
                (*stages.sources, "ICRP89", "DuBois1916"),
                unit=U.l_per_min,
                kind="egfr",
                organ=Organ.KIDNEY,
            )
        )
        anchors = _entry(ckd, "stage_anchor")
        if len(anchors.values) != len(_STAGES):
            raise ValueError(
                f"{DATA}.toml: {len(anchors.values)} anchors for the stages {_STAGES}."
            )
        # the eGFR (l/min) of 1 and of the values of the stages, decreasing
        points = [stages.values[0], *anchors.values]
        if any(a <= b for a, b in itertools.pairwise(points)):
            raise ValueError(f"{DATA}.toml: the anchors {points} do not decrease.")
        limits = points
        for key in (
            "plasma_protein",
            "hematocrit",
            "kidney_volume",
            "kidney_flow",
            "gastric_emptying_time",
            "si_transit_time",
        ):
            self._stage(
                key,
                class_factors("ckd", key),
                f"on the {key.replace('_', ' ')}",
                limits,
            )
        for protein in self.liver_proteins(ckd["proteins"]):
            self._stage(
                protein,
                _entry(ckd["proteins"], protein),
                f"on the activity of {protein} in the liver (with "
                f"`{ids.CKD_HEPATIC}` = 1)",
                limits,
            )
        for protein in self.binding_proteins():
            for key in _binding_keys(protein):
                self._stage(
                    key,
                    class_factors("ckd.binding", key),
                    f"on the {ckd['binding'][key]['quantity']}",
                    limits,
                )

    def _stage(
        self, key: str, entry: ImpairmentFactors, text: str, limits: Sequence[float]
    ) -> None:
        """A factor of chronic kidney disease, piecewise linear in the eGFR.

        Args:
            key: the key of the factor.
            entry: the values of mild, moderate, severe and end stage CKD.
            text: the description of the factor.
            limits: the eGFR per 1.73 m2 (ml/min), decreasing, of normal
                function (the factor 1, from there up) and of the values of
                mild, moderate, severe and end stage CKD (the factor is
                constant from the last one down).
        """
        provenance, sources = entry.provenance(_STAGES)
        values = [1.0, *entry.values]
        egfr = ids.EGFR
        terms = [_number(values[0]), f"{egfr} >= {_l_per_min(limits[0])}"]
        for i in range(len(values) - 1):
            a, b = values[i], values[i + 1]
            x0, x1 = limits[i], limits[i + 1]
            line = (
                _number(a)
                if a == b
                else f"{_number(a)} + {_number(b - a)} * ({_l_per_min(x0)} - {egfr})"
                f" / {_l_per_min(x0 - x1)}"
            )
            terms += [line, f"{egfr} >= {_l_per_min(x1)}"]
        terms.append(_number(values[-1]))
        self.add(
            _initial(
                ids.condition_factor("ckd", key),
                f"piecewise({', '.join(terms)})",
                f"CKD factor {key}",
                f"Factor of chronic kidney disease {text}, linear in the eGFR "
                f"between the anchors of the stages (normal 1), {provenance}",
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
        for key in _HF_ORGANS:
            factors = class_factors("heart_failure", key)
            provenance, sources = factors.provenance(_HF_CLASSES)
            self.add(
                _initial(
                    ids.condition_factor("hf", key),
                    _interpolation(ids.HF_SEVERITY, [1.0, *factors.values]),
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

        female_albumin = ALBUMIN in sex_curves(self.definition.sex)
        human = self.definition.age and self.definition.species is Species.HUMAN
        rhythms = circadian_keys(self.definition)
        pregnant = self.definition.pregnancy is not None

        def preg(key: str) -> str | None:
            """The factor of the pregnancy on a quantity, if the model is pregnant."""
            return pregnancy_factor(key) if pregnant else None

        def circ(key: str) -> str | None:
            """The circadian factor of a quantity, if the model has the rhythm."""
            return circadian_factor(key) if key in rhythms else None

        adipose = self.definition.adipose
        obese = adipose is not None and adipose.obesity
        obesity_keys = (
            {
                ids.CO,
                *(ids.organ_flow(o) for os in FLOW_FACTORS.values() for o in os),
                *BINDING_FACTORS,
                *(f"{o.value}_{p.lower()}" for o, p in _obesity_proteins()),
            }
            if obese
            else set()
        )

        def obes(key: str) -> str | None:
            """The factor of obesity on a quantity, if the build has the factors."""
            return obesity_factor(key) if key in obesity_keys else None

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
                    product(
                        cp_flow or rest,
                        ckd_flow,
                        hf_flow,
                        preg(ORGAN_FLOW),
                        preg(q) if o is Organ.KIDNEY else None,
                        obes(ids.CO),
                        obes(q),
                    ),
                    f"impairment factor {q}",
                    f"{text} the blood flow `{q}`; the cardiac output is the "
                    f"sum of the flows (`{ids.impairment_factor(ids.CO)}`).",
                    kind="impairment_factor",
                    organ=o,
                )
            )
            if cp_flow is None:
                closing.setdefault(
                    product(ckd_flow, hf_flow, preg(ORGAN_FLOW), obes(ids.CO), obes(q)),
                    [],
                ).append(ids.flow_fraction(o))
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
        if pregnant:
            flows.append(pregnancy_factor(PLACENTA_SHARE))
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
                    circ(ids.GFR),
                    preg(ids.GFR),
                ),
                "the GFR (heart failure: the GFR follows the renal blood flow; "
                f"chronic kidney disease is `{ids.F_RENAL_FUNCTION}`)",
            ),
            ids.HCT: (
                product(
                    ids.condition_factor("cp", "hematocrit"),
                    ids.condition_factor("ckd", "hematocrit"),
                    preg(ids.HCT),
                ),
                "the hematocrit of the sex of the build; the anaemia of "
                "cirrhosis and of chronic kidney disease combine "
                "multiplicatively",
            ),
            ids.organ_mass(Organ.KIDNEY): (
                ids.condition_factor("ckd", "kidney_volume"),
                "the kidney mass (adipose, the remainder, takes the difference)",
            ),
            "albumin": (
                product(
                    ids.condition_factor("cp", "albumin"),
                    age_factor("albumin") if human else None,
                    sex_factor(ALBUMIN) if female_albumin else None,
                    circ("albumin"),
                    preg("albumin"),
                    obes("albumin"),
                ),
                "the plasma albumin"
                + (
                    " (and the age" + (" and the sex)" if female_albumin else ")")
                    if human
                    else ""
                ),
            ),
            "agp": (
                product(
                    ids.condition_factor("cp", "agp"),
                    age_factor("agp") if human else None,
                    preg("agp"),
                    obes("agp"),
                ),
                "the plasma alpha-1-acid glycoprotein"
                + (" (and the age)" if human else ""),
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
                product(ids.condition_factor("ckd", "si_transit_time"), preg("t_si")),
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
        for protein in self.binding_proteins():
            keys = _binding_keys(protein)
            self.add(
                _initial(
                    _ckd_binding(protein),
                    " * ".join(ids.condition_factor("ckd", k) for k in keys),
                    f"impairment factor {protein} ckd",
                    "The factor of chronic kidney disease on the binding to "
                    f"{protein}, the product of its rows of the CKD profile "
                    f"({', '.join(keys)}); it multiplies the binding factor of fu "
                    f"of a substance bound to {protein} with the CKD factor, in "
                    f"place of `{ids.impairment_factor('plasma_protein')}`.",
                    kind="impairment_factor",
                )
            )
        cp_proteins = set(self.liver_proteins(child_pugh_proteins()))
        ckd_proteins = set(self.liver_proteins(self.data["ckd"]["proteins"]))
        for organ, protein in self.individual_proteins():
            liver = organ is Organ.LIVER
            hepatic = (
                f"({_ONE} + {ids.CKD_HEPATIC} * ({ids.condition_factor('ckd', protein)} - {_ONE}))"
                if liver and protein in ckd_proteins
                else None
            )
            key = f"{organ.value}_{protein.lower()}"
            self.add(
                _initial(
                    ids.impairment_activity(protein, organ),
                    product(
                        ids.condition_factor("cp", protein)
                        if liver and protein in cp_proteins
                        else None,
                        hepatic,
                        age_factor(key)
                        if human and protein_factor(organ, protein) is not None
                        else None,
                        sex_factor(key)
                        if protein_sex_factor(self.definition.sex, organ, protein)
                        is not None
                        else None,
                        circ(key),
                        preg(key)
                        if (organ, protein) in _pregnancy_proteins()
                        else None,
                        obes(key),
                    ),
                    f"individual factor activity {protein} {organ.name.lower()}",
                    f"The factor of the individual (the impairments, the age and "
                    f"the sex) on the activity of {protein} in the "
                    f"{organ.name.lower().replace('_', ' ')}, the input "
                    f"`{ids.impairment_activity(protein)}` of the organ; its rates "
                    f"run with `{ids.activity(protein)}` times it.",
                    kind="impairment_factor",
                    organ=organ,
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
    if definition.species is Species.HUMAN:
        impairment.sex()
    if definition.age and definition.species is Species.HUMAN:
        impairment.age()
    if definition.adipose is not None and definition.adipose.obesity:
        impairment.obesity()
    impairment.combined()
    return impairment.objects
