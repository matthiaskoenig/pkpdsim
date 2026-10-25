"""The annotation of every element the library creates (#63).

Every generator annotates an element when it creates it, with
`annotate(element, kind, ...)`: the kind names what the element is (an organ
volume, the plasma of an organ, a metabolism) and selects one row of the
table `KINDS`, which holds its SBO term and its ontology terms. The context
of the element adds the terms of where it is: the organ (UBERON and the
compartment terms of `pkpdsim.annotations`, the PBPKO compartment), the
blood pool, the lumen segment, the protein of a process (UniProt, NCIT) and
the substance of a species (the annotations of its specification).

The PBPKO terms are those of the pinned pymetadata (`PBPKO`): a kind gives
its PBPKO term, and an organ-specific term (`VOLUME_OF_LIVER`) where PBPKO
has one for the organ. Where PBPKO has no term, the kind takes the nearest
broader PBPKO term and names the missing one (`requested`); every missing
term is listed in `docs/pbpk/annotation.md` ("terms requested upstream",
`requested_terms`).

Every resource is written in one format, the identifiers.org
`collection/term` of `pkpdsim.annotations`: `pbpko/PBPKO:00558`,
`uberon/UBERON:0002107`, `uniprot/P08684`.

The qualifiers: `BQB.IS` for what the element is (the PBPKO term, the kind's
terms, the protein of a process, the substance of a species); the location
is `BQB.IS_PART_OF` for a compartment, `BQB.IS_PROPERTY_OF` for a parameter
and `BQB.OCCURS_IN` for a reaction; a parameter of a protein (activity,
inhibition constant) is `BQB.IS_PROPERTY_OF` the protein.
"""

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field

from pymetadata.ontologies.pbpko import PBPKO
from sbmlutils.factory import (
    Compartment,
    InitialAssignment,
    Reaction,
    Sbase,
    Species,
    UnitDefinition,
)
from sbmlutils.metadata import BQB, SBO

from pkpdsim.annotations import compartments as _compartments
from pkpdsim.pbpk.organs import Organ, Pool
from pkpdsim.pbpk.proteins import PROTEINS
from pkpdsim.pbpk.substance import Mode, Substance

__all__ = [
    "KINDS",
    "Kind",
    "Term",
    "annotate",
    "initial_assignment",
    "model_annotations",
    "pbpko",
    "requested_terms",
]

PBPKO_URL = "https://bioregistry.io/"
"""The prefix of the resources of PBPKO terms, `https://bioregistry.io/PBPKO:00558`."""

Term = tuple[BQB, str]
"""An annotation: the qualifier and the resource, `collection/term`."""


def pbpko(name: str) -> str:
    """Get the resource of a term of the pinned PBPKO by its name.

    Args:
        name: the name of the term in the pymetadata enum, `LIVER_COMPARTMENT`.

    Returns:
        The resource `https://bioregistry.io/PBPKO:00558`: identifiers.org
        does not register PBPKO, so pymetadata cannot parse a
        `pbpko/PBPKO:00558` (it logs a MIRIAM error and writes an
        identifiers.org url which does not resolve); the bioregistry url is
        written as given and resolves.

    Raises:
        AttributeError: if PBPKO has no term of the name.
    """
    term = getattr(PBPKO, name)
    return f"{PBPKO_URL}{term.curie}"


def _has(name: str) -> bool:
    """Check whether the pinned PBPKO has a term of the name."""
    return hasattr(PBPKO, name)


# --- the organs, pools and segments ----------------------------------------


@dataclass(frozen=True)
class _Place:
    """The terms of a place: an organ, a blood pool or a lumen segment.

    Attributes:
        label: the name of the place in the requested terms.
        terms: the terms of the place (`pkpdsim.annotations`).
        compartment: the PBPKO compartment of the place, `None` if PBPKO has
            none.
        tokens: the tokens of the place in the names of the organ-specific
            PBPKO terms (`VOLUME_OF_{token}`), the first one PBPKO has is
            taken; none for no such terms.
    """

    label: str
    terms: tuple[Term, ...]
    compartment: str | None
    tokens: tuple[str, ...] = ()


def _location_terms(key: str) -> tuple[Term, ...]:
    """The terms of where a compartment of `pkpdsim.annotations` is, not of what it is.

    The blood pools and the veins of `pkpdsim.annotations` are plasma
    (`BQB.IS`) which is part of a vessel (`BQB.IS_PART_OF`); the place of an
    element in a pool is the vessel, the kind of the element says whether it
    is plasma.
    """
    return tuple(t for t in _compartments[key] if t[0] is BQB.IS_PART_OF)


def _organ(
    key: str, label: str, tokens: tuple[str, ...], compartment: str | None
) -> _Place:
    """An organ, its terms from the compartments of `pkpdsim.annotations`."""
    return _Place(
        label=label,
        terms=tuple(_compartments[key]),
        compartment=compartment,
        tokens=tokens,
    )


ORGANS: dict[Organ, _Place] = {
    Organ.LUNG: _organ("lu", "lung", ("LUNG",), "LUNG_COMPARTMENT"),
    Organ.HEART: _organ("he", "heart", ("HEART",), "HEART_COMPARTMENT"),
    Organ.BRAIN: _organ("br", "brain", ("BRAIN",), "BRAIN_COMPARTMENT"),
    Organ.LIVER: _organ("li", "liver", ("LIVER",), "LIVER_COMPARTMENT"),
    Organ.KIDNEY: _organ("ki", "kidney", ("KIDNEY",), "KIDNEY_COMPARTMENT"),
    Organ.STOMACH: _organ("stomach", "stomach", ("STOMACH",), "STOMACH_COMPARTMENT"),
    Organ.SMALL_INTESTINE: _organ(
        "si",
        "small intestine",
        ("SMALL_INTESTINE", "INTESTINE"),
        "SMALL_INTESTINE_COMPARTMENT",
    ),
    Organ.LARGE_INTESTINE: _organ(
        "la",
        "large intestine",
        ("LARGE_INTESTINE", "INTESTINE"),
        "LARGE_INTESTINE_COMPARTMENT",
    ),
    Organ.SPLEEN: _organ("sp", "spleen", ("SPLEEN",), "SPLEEN_COMPARTMENT"),
    Organ.PANCREAS: _organ("pa", "pancreas", ("PANCREAS",), "PANCREAS_COMPARTMENT"),
    Organ.MUSCLE: _organ("mu", "muscle", ("MUSCLE",), "MUSCLE_COMPARTMENT"),
    Organ.SKIN: _organ("sk", "skin", ("SKIN",), "SKIN_COMPARTMENT"),
    Organ.ADIPOSE: _organ("ad", "adipose", ("FAT",), "ADIPOSE_COMPARTMENT"),
    Organ.BONE: _organ("bo", "bone", ("BONE",), None),
    Organ.GONADS: _organ("go", "gonads", ("GONAD",), "GONAD_COMPARTMENT"),
}
"""The terms of the organs; the PBPKO token of adipose tissue is `FAT`, the
intestines take the terms of the intestine where PBPKO has none of their own
(`MAXIMUM_RATE_IN_INTESTINE`); PBPKO has no compartment of the bone (only of
the bone marrow)."""

POOLS: dict[Pool, _Place] = {
    Pool.ARTERIAL: _Place(
        "arterial blood",
        ((BQB.IS_PART_OF, "uberon/UBERON:0013755"), *_location_terms("ar")),
        "ARTERIAL_BLOOD_COMPARTMENT",
    ),
    Pool.VENOUS: _Place(
        "venous blood",
        ((BQB.IS_PART_OF, "uberon/UBERON:0013756"), *_location_terms("ve")),
        "VENOUS_BLOOD_COMPARTMENT",
    ),
    Pool.PORTAL: _Place(
        "portal vein blood",
        ((BQB.IS_PART_OF, "uberon/UBERON:0002017"), *_location_terms("po")),
        None,
    ),
}
"""The terms of the blood pools; PBPKO has no compartment of the portal vein."""

SEGMENTS: dict[str, _Place] = {
    "st": _Place("stomach", tuple(_compartments["stomach"]), "STOMACH_COMPARTMENT"),
    "du": _Place("duodenum", tuple(_compartments["duodenum"]), "DUODENUM_COMPARTMENT"),
    "jej1": _Place(
        "upper jejunum", tuple(_compartments["jejunum"]), "JEJUNUM_COMPARTMENT"
    ),
    "jej2": _Place(
        "lower jejunum", tuple(_compartments["jejunum"]), "JEJUNUM_COMPARTMENT"
    ),
    "ile1": _Place("upper ileum", tuple(_compartments["ileum"]), "ILEUM_COMPARTMENT"),
    "ile2": _Place("lower ileum", tuple(_compartments["ileum"]), "ILEUM_COMPARTMENT"),
    "cae": _Place("caecum", tuple(_compartments["caecum"]), "CECUM_COMPARTMENT"),
    "asc": _Place(
        "ascending colon",
        tuple(_compartments["ascending colon"]),
        "ASCENDING_COLON_COMPARTMENT",
    ),
    "tra": _Place(
        "transverse colon",
        tuple(_compartments["transverse colon"]),
        "TRANSVERSE_COLON_COMPARTMENT",
    ),
    "des": _Place(
        "descending colon",
        tuple(_compartments["descending colon"]),
        "DESCENDING_COLON_COMPARTMENT",
    ),
    "sig": _Place(
        "sigmoid colon",
        tuple(_compartments["sigmoid colon"]),
        "SIGMOID_COLON_COMPARTMENT",
    ),
    "rec": _Place("rectum", tuple(_compartments["rectum"]), "RECTUM_COMPARTMENT"),
    "pt": _Place(
        "proximal tubule",
        tuple(_compartments["proximal tubule"]),
        "PROXIMAL_TUBULE_COMPARTMENT",
    ),
    "lh": _Place(
        "loop of Henle",
        tuple(_compartments["loop of Henle"]),
        "LOOP_OF_HENLE_COMPARTMENT",
    ),
    "dt": _Place(
        "distal tubule",
        tuple(_compartments["distal tubule"]),
        "DISTAL_TUBULE_COMPARTMENT",
    ),
    "cd": _Place(
        "collecting duct",
        tuple(_compartments["collecting duct"]),
        "COLLECTING_DUCT_COMPARTMENT",
    ),
}
"""The terms of the lumen segments of the intestine and of the tubule of the
nephron kidney, by their id token."""

SITES: dict[str, _Place] = {
    "hepatic_vein": _Place(
        "hepatic vein",
        (*_location_terms("hv"), (BQB.IS_PART_OF, "uberon/UBERON:0001143")),
        None,
    ),
    "arm_vein": _Place("arm vein", _location_terms("fov"), None),
}
"""The terms of the sampling sites of the observers which are no blood pool."""


# --- the kinds -------------------------------------------------------------


@dataclass(frozen=True)
class Kind:
    """What an element is: its SBO term and its ontology terms.

    Attributes:
        sbo: the SBO term, `None` keeps the SBO term of the element.
        pbpko: the name of the PBPKO term, the broader term where PBPKO has
            no term for the kind (`requested`), `None` for no PBPKO term.
        organ: the name of the organ-specific PBPKO term with `{}` for the
            token of the organ (`VOLUME_OF_{}`), which replaces `pbpko`
            where PBPKO has it.
        pools: the PBPKO term of the kind by blood pool, which replaces
            `pbpko`.
        terms: further terms of the kind (`BQB.IS`, other ontologies).
        requested: the PBPKO term which is missing, if `pbpko` is a broader
            term.
        per_organ: the organ-specific terms which PBPKO lacks are requested
            for every organ (the physiology), not only where PBPKO has some
            (the kinetics of the processes, `MICHAELIS_CONSTANT__LIVER`).
        is_place: the element is the place itself (the whole organ of a
            flow-limited organ), so the terms of the place are `BQB.IS`.
    """

    sbo: str | None
    pbpko: str | None
    organ: str | None = None
    pools: Mapping[Pool, str] = field(default_factory=dict)
    terms: tuple[Term, ...] = ()
    requested: str | None = None
    per_organ: bool = True
    is_place: bool = False


_PHYS = "PHYSIOLOGICAL_PARAMETER"
_BIOCHEM = "BIOCHEMICAL_PARAMETER"
_MISC = "MISCELLANEOUS_PARAMETER"
_P = SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER
_C = SBO.PHYSICAL_COMPARTMENT
_SPECIES = SBO.SIMPLE_CHEMICAL
_T = SBO.TRANSPORT_REACTION
_PLASMA = tuple(_compartments["plasma"])
_IV_ROUTE: Term = (
    BQB.IS_PROPERTY_OF,
    "ncit/C38276",  # intravenous route of administration
)
_IV_BOLUS: tuple[Term, ...] = (
    (BQB.IS_PROPERTY_OF, pbpko("INTRAVENOUS_BOLUS_EXPOSURE")),
    _IV_ROUTE,
)
_IV_INFUSION: tuple[Term, ...] = (
    (BQB.IS_PROPERTY_OF, pbpko("INTRAVENOUS_INFUSION_EXPOSURE")),
    _IV_ROUTE,
)
_BLOOD_CELLS: tuple[Term, ...] = ((BQB.IS, "cl/CL:0000232"),)  # erythrocyte
_DISEASES: dict[str, Term] = {
    "cirrhosis": (BQB.IS_PROPERTY_OF, "ncit/C2951"),  # cirrhosis
    "child_pugh": (BQB.IS, "ncit/C121007"),  # Child-Pugh clinical classification
    "ckd": (BQB.IS_PROPERTY_OF, "ncit/C80078"),  # chronic kidney disease
    "heart_failure": (BQB.IS_PROPERTY_OF, "ncit/C50577"),  # heart failure
}

KINDS: dict[str, Kind] = {
    # compartments
    "plasma": Kind(
        _C,
        "PLASMA_COMPARTMENT",
        pools={
            Pool.ARTERIAL: "ARTERIAL_PLASMA",
            Pool.VENOUS: "VENOUS_PLASMA_COMPARTMENT",
        },
        terms=_PLASMA,
    ),
    "blood_cells": Kind(_C, "RED_BLOOD_CELL_COMPARTMENT", terms=_BLOOD_CELLS),
    "capacity": Kind(
        SBO.IMPLICIT_COMPARTMENT,
        "COMPARTMENT",
        terms=((BQB.IS_VERSION_OF, pbpko("BLOOD_COMPARTMENT")),),
        requested="plasma-equivalent capacity of the blood (volume times BP)",
    ),
    "organ": Kind(_C, "COMPARTMENT", is_place=True),
    "tissue": Kind(
        _C,
        "COMPARTMENT",
        terms=((BQB.IS, "uberon/UBERON:0000479"),),  # tissue
        requested="tissue compartment (interstitial and intracellular space)",
    ),
    "interstitial": Kind(
        _C,
        "COMPARTMENT",
        terms=(
            (BQB.IS, "uberon/UBERON:0000913"),  # interstitial fluid
            (BQB.IS, "ncit/C120839"),  # interstitial fluid
        ),
        requested="interstitial compartment",
    ),
    "intracellular": Kind(
        _C,
        "COMPARTMENT",
        organ="{}_INTRACELLULAR_SPACE_COMPARTMENT",
        terms=((BQB.IS, "go/GO:0005622"),),  # intracellular anatomical structure
        requested="intracellular space compartment",
        per_organ=False,
    ),
    "urine": Kind(_C, "URINE_COMPARTMENT", terms=tuple(_compartments["urine"])),
    "feces": Kind(_C, "FECES_COMPARTMENT", terms=tuple(_compartments["feces"])),
    "bile": Kind(_C, "BILE_COMPARTMENT", terms=tuple(_compartments["bi"])),
    "lumen": Kind(_C, "GUT_LUMEN_COMPARTMENT", terms=tuple(_compartments["gu_lumen"])),
    "mucosa": Kind(
        _C,
        "COMPARTMENT",
        terms=(
            (BQB.IS, "uberon/UBERON:0001242"),  # intestinal mucosa
            (BQB.IS_PART_OF, pbpko("INTESTINAL_TISSUE_COMPARTMENT")),
        ),
        requested="intestinal mucosa compartment (enterocytes)",
    ),
    "gut_wall": Kind(
        _C,
        "COMPARTMENT",
        terms=((BQB.IS_PART_OF, pbpko("INTESTINAL_TISSUE_COMPARTMENT")),),
        requested="intestinal tissue compartment without the mucosa",
    ),
    "tubule": Kind(_C, "FILTRATE_COMPARTMENT"),
    # covariates and the body
    "body_weight": Kind(_P, "BODYWEIGHT"),
    "height": Kind(_P, "HEIGHT"),
    "age": Kind(_P, "AGE_IN_PBPK_MODEL"),
    "bmi": Kind(_P, "BODY_MASS_INDEX"),
    "bsa": Kind(_P, "SURFACE_AREA", terms=((BQB.IS, "ncit/C25157"),)),
    "lean_mass": Kind(
        _P, _PHYS, terms=((BQB.IS, "ncit/C71258"),), requested="fat-free mass"
    ),
    "scaling_factor": Kind(_P, _PHYS, requested="covariate scaling factor"),
    "organ_mass": Kind(
        _P, _PHYS, organ="{}_MASS", requested="mass of organ", per_organ=False
    ),
    "organ_volume": Kind(_P, "VOLUME_OF_COMPARTMENT", organ="VOLUME_OF_{}"),
    "blood_volume": Kind(_P, "VOLUME_OF_BLOOD"),
    "pool_volume": Kind(
        _P,
        "VOLUME_OF_COMPARTMENT",
        pools={
            Pool.ARTERIAL: "ARTERIAL_VOLUME",
            Pool.VENOUS: "VENOUS_VOLUME",
            Pool.PORTAL: "VOLUME_OF_PORTAL_VEIN",
        },
    ),
    "pool_fraction": Kind(
        _P, _PHYS, requested="fraction of the blood volume in a blood pool"
    ),
    "vascular_fraction": Kind(_P, "VOLUME_FRACTION_OF_BLOOD"),
    "interstitial_fraction": Kind(
        _P, _PHYS, requested="interstitial volume fraction of an organ"
    ),
    "tissue_volume": Kind(_P, "VOLUME_OF_COMPARTMENT"),
    "zone_fraction": Kind(
        _P, _PHYS, requested="volume fraction of a zone of the liver"
    ),
    "hematocrit": Kind(
        _P,
        "VOLUME_FRACTION_OF_RED_BLOOD_CELLS",
        terms=((BQB.IS_PROPERTY_OF, "uberon/UBERON:0000178"),),  # blood
    ),
    "cardiac_output": Kind(_P, "CARDIAC_OUTPUT_RATE"),
    "flow_fraction": Kind(
        _P,
        "FRACTION_OF_BLOOD_FLOW_TO_A_COMPARTMENT",
        organ="FRACTION_OF_BLOOD_FLOW_TO_{}",
    ),
    "organ_flow": Kind(
        _P, "BLOOD_FLOW_RATE_TO_COMPARTMENT", organ="BLOOD_FLOW_RATE_TO_{}"
    ),
    "hepatic_artery_flow": Kind(
        _P, "BLOOD_FLOW_RATE_TO_COMPARTMENT", requested="hepatic artery blood flow"
    ),
    "portal_flow": Kind(_P, "BLOOD_FLOW_RATE_TO_PORTAL_VEIN"),
    "hepatic_flow": Kind(_P, "HEPATIC_VEIN_BLOOD_FLOW"),
    "gfr": Kind(_P, "GLOMERULAR_FILTRATION_RATE"),
    "egfr": Kind(_P, "GLOMERULAR_FILTRATION_RATE", terms=((BQB.IS, "ncit/C110935"),)),
    "organ_function": Kind(_P, _PHYS, requested="relative organ function"),
    "tubular_share": Kind(
        _P,
        _PHYS,
        requested="share of the tubular capacity which follows the renal function",
    ),
    "shunt": Kind(_P, _PHYS, requested="portosystemic shunt fraction"),
    # impairment
    "child_pugh": Kind(
        _P,
        _PHYS,
        terms=(_DISEASES["child_pugh"],),
        requested="severity of hepatic impairment",
    ),
    "heart_failure": Kind(
        _P,
        _PHYS,
        terms=(_DISEASES["heart_failure"],),
        requested="severity of heart failure",
    ),
    "ckd": Kind(
        _P,
        _PHYS,
        terms=(_DISEASES["ckd"],),
        requested="stage of chronic kidney disease",
    ),
    "cirrhosis_factor": Kind(
        _P, _PHYS, terms=(_DISEASES["cirrhosis"],), requested="impairment factor"
    ),
    "ckd_factor": Kind(
        _P, _PHYS, terms=(_DISEASES["ckd"],), requested="impairment factor"
    ),
    "heart_failure_factor": Kind(
        _P, _PHYS, terms=(_DISEASES["heart_failure"],), requested="impairment factor"
    ),
    "impairment_factor": Kind(_P, _PHYS, requested="impairment factor"),
    # substances
    "molecular_weight": Kind(SBO.MOLECULAR_MASS, "MOLECULAR_WEIGHT"),
    "fu": Kind(_P, "FRACTION_UNBOUND_PLASMA"),
    "fraction_reabsorbed": Kind(_P, "RESORPTION"),
    "kp": Kind(_P, "TISSUE_PLASMA_PARTITION_COEFFICIENT"),
    "kp_organ": Kind(
        _P,
        "TISSUE_PLASMA_PARTITION_COEFFICIENT",
        organ="{}_PLASMA_PARTITION_COEFFICIENT",
    ),
    "kp_factor": Kind(_P, _BIOCHEM, requested="factor of the partition coefficients"),
    "partition": Kind(_P, "PARTITION_COEFFICIENT"),
    "tissue_composition": Kind(
        _P, _PHYS, requested="tissue composition (water, protein, lipid fractions)"
    ),
    "permeability": Kind(
        _P, _BIOCHEM, requested="permeability surface area product (exchange flow)"
    ),
    "bp": Kind(_P, "BLOOD_PLASMA_RATIO"),
    "kbc": Kind(_P, "RED_BLOOD_CELLS_PARTITION_COEFFICIENT"),
    "exchange_rate": Kind(
        _P, _BIOCHEM, requested="exchange rate constant of plasma and blood cells"
    ),
    "endothelial_exchange": Kind(
        _P, _PHYS, requested="endothelial exchange rate constant per organ volume"
    ),
    # proteins and interactions
    "activity": Kind(
        SBO.RELATIVE_ACTIVITY,
        _BIOCHEM,
        requested="relative activity of enzyme or transporter",
    ),
    "enzyme_amount": Kind(
        SBO.RELATIVE_ACTIVITY,
        _BIOCHEM,
        requested="relative amount of enzyme or transporter",
    ),
    "kdeg": Kind(SBO.DECAY_CONSTANT, "DEGRADATION_RATE_CONSTANT"),
    "ki": Kind(SBO.INHIBITORY_CONSTANT, "INHIBITION_CONSTANT"),
    "kinact": Kind(
        SBO.INACTIVATION_RATE_CONSTANT,
        _BIOCHEM,
        requested="maximal inactivation rate constant",
    ),
    "emax": Kind(_P, _BIOCHEM, requested="maximal induction"),
    "ec50": Kind(SBO.EC50, _BIOCHEM, requested="half-maximal effective concentration"),
    "placement_weight": Kind(
        _P, _MISC, requested="abundance weight of a protein or process"
    ),
    # kinetics of the processes
    "vmax": Kind(
        SBO.MAXIMAL_VELOCITY,
        "MAXIMUM_RATE",
        organ="MAXIMUM_RATE_IN_{}",
        per_organ=False,
    ),
    "vmax_basolateral": Kind(
        SBO.MAXIMAL_VELOCITY, "MAXIMUM_RATE_OF_BASOLATERAL_TRANSPORTER"
    ),
    "vmax_apical": Kind(SBO.MAXIMAL_VELOCITY, "MAXIMUM_RATE_OF_APICAL_TRANSPORTER"),
    "km": Kind(
        SBO.MICHAELIS_CONSTANT,
        _BIOCHEM,
        organ="MICHAELIS_CONSTANT__{}",
        requested="Michaelis constant",
        per_organ=False,
    ),
    "km_basolateral": Kind(SBO.MICHAELIS_CONSTANT, "MICHAELIS_CONSTANT_BASO"),
    "km_apical": Kind(SBO.MICHAELIS_CONSTANT, "MICHAELIS_CONSTANT_APICAL"),
    "hill": Kind(SBO.HILL_COEFFICIENT, "HILL_COEFFICIENT"),
    "k_metabolism": Kind(SBO.KINETIC_CONSTANT, "INTRINSIC_CLEARANCE_RATE"),
    "k_influx": Kind(SBO.KINETIC_CONSTANT, "CLEARANCE_RATE"),
    "k_efflux": Kind(SBO.KINETIC_CONSTANT, "CLEARANCE_RATE"),
    # a clearance per litre of organ (1/min) is no clearance (a volume per
    # time) nor a rate constant of an amount: the broader term (R49)
    "k_metabolism_volume": Kind(
        SBO.KINETIC_CONSTANT,
        _BIOCHEM,
        terms=((BQB.IS_VERSION_OF, pbpko("INTRINSIC_CLEARANCE_RATE")),),
        requested="intrinsic clearance per organ volume",
    ),
    "k_transport_volume": Kind(
        SBO.KINETIC_CONSTANT,
        _BIOCHEM,
        terms=((BQB.IS_VERSION_OF, pbpko("CLEARANCE_RATE")),),
        requested="transport clearance per organ volume",
    ),
    "custom_parameter": Kind(None, _MISC),
    "custom_species": Kind(None, "CONCENTRATION_OF_COMPOUND_IN_COMPARTMENT"),
    "custom_amount": Kind(None, "AMOUNT_OF_COMPOUND_IN_COMPARTMENT"),
    # routes, intestine, sinks
    "dose": Kind(_P, "PARAMETER", terms=((BQB.IS, "ncit/C25488"),), requested="dose"),
    "iv_time": Kind(
        _P,
        "PARAMETER",
        terms=_IV_BOLUS,
        requested="injection time of an intravenous dose",
    ),
    "iv_rate": Kind(
        SBO.FORWARD_RATE_CONSTANT,
        "PARAMETER",
        terms=_IV_BOLUS,
        requested="injection rate constant of an intravenous dose",
    ),
    "infusion_rate": Kind(
        _P,
        "PARAMETER",
        terms=_IV_INFUSION,
        requested="infusion rate of an intravenous infusion",
    ),
    "dissolution": Kind(_P, _MISC, requested="dissolution of a formulation"),
    "application_rate": Kind(SBO.FORWARD_RATE_CONSTANT, "ABSORPTION_RATE_CONSTANT"),
    "absorption_rate": Kind(
        SBO.FORWARD_RATE_CONSTANT, "ABSORPTION_RATE_CONSTANT_INTESTINE"
    ),
    "transit_time": Kind(_P, "INTESTINAL_TRANSIT_TIME"),
    "gastric_transit_time": Kind(_P, _PHYS, requested="gastric emptying time"),
    "fed": Kind(
        _P,
        _MISC,
        terms=((BQB.IS_VERSION_OF, "ncit/C65056"),),  # postprandial
        requested="prandial state (fed or fasted)",
    ),
    "intestine_geometry": Kind(
        _P, _PHYS, requested="intestinal mucosa and absorption geometry"
    ),
    "gut_flow": Kind(_P, "BLOOD_FLOW_RATE_TO_GUT"),
    "gut_volume": Kind(_P, "VOLUME_OF_GUT"),
    "mucosa_flow": Kind(
        _P,
        "BLOOD_FLOW_RATE_TO_COMPARTMENT",
        requested="blood flow rate to intestinal mucosa",
    ),
    "mucosa_volume": Kind(_P, "ENTEROCYTE_VOLUME"),
    "bile_transit": Kind(
        SBO.FORWARD_RATE_CONSTANT, _PHYS, requested="bile transit rate constant"
    ),
    "tubule_fraction": Kind(
        _P,
        _PHYS,
        requested="tubular flow fraction of the glomerular filtration rate",
    ),
    "fraction_absorbed": Kind(_P, "UPTAKE_FRACTION_TO_BLOOD"),
    "arm_weight": Kind(
        _P, _PHYS, requested="weight of a tissue in the arm vein sampling site"
    ),
    # observers and the outflow of the organs
    "outflow": Kind(_P, "TOTAL_VENOUS_CONCENTRATION_OF_COMPOUND_IN_COMPARTMENT"),
    "outflow_blood_cells": Kind(_P, "CONCENTRATION_OF_COMPOUND_IN_COMPARTMENT"),
    "plasma_concentration": Kind(_P, "CONCENTRATION_IN_PLASMA"),
    "free_concentration": Kind(_P, "UNBOUND_CONCENTRATION_IN_PLASMA"),
    "blood_concentration": Kind(_P, "CONCENTRATION_IN_BLOOD"),
    "unbound_concentration": Kind(
        _P, "UNBOUND_CONCENTRATION_OF_COMPUND_IN_COMPARTMENT"
    ),
    # the amount per time into the urine is a result of the model; PBPKO's
    # urinary rate is a biochemical parameter (an input), so a version of it
    "urine_rate": Kind(
        _P,
        "OUTPUT_PARAMETER",
        terms=((BQB.IS_VERSION_OF, pbpko("URINARY_RATE")),),
        requested="urinary excretion rate (amount of compound per time)",
    ),
    # an observer of a definition: what it computes is the definition's, the
    # library knows only that it is a result of the model
    "observer": Kind(_P, "OUTPUT_PARAMETER"),
    # species
    "concentration": Kind(_SPECIES, "CONCENTRATION_OF_COMPOUND_IN_COMPARTMENT"),
    "concentration_plasma": Kind(_SPECIES, "CONCENTRATION_IN_PLASMA"),
    "concentration_tissue": Kind(
        _SPECIES,
        "CONCENTRATION_OF_COMPOUND_IN_COMPARTMENT",
        organ="CONCENTRATION_IN_{}",
    ),
    "concentration_lumen": Kind(_SPECIES, "CONCENTRATION_IN_LUMEN"),
    "concentration_mucosa": Kind(_SPECIES, "CONCENTRATION_IN_INTESTINAL_TISSUE"),
    "concentration_bile": Kind(_SPECIES, "CONCENTRATION_IN_BILE"),
    "concentration_tubule": Kind(_SPECIES, "CONCENTRATION_IN_FILTRATE"),
    "amount_urine": Kind(_SPECIES, "CUMULATIVE_AMOUNT_IN_URINE"),
    "amount_feces": Kind(_SPECIES, "CUMULATIVE_AMOUNT_IN_FECES"),
    "amount_undissolved": Kind(
        _SPECIES, "AMOUNT_IN_LUMEN", requested="amount of undissolved drug"
    ),
    # reactions
    "blood_flow": Kind(
        _T,
        "DISTRIBUTION_PROCESS",
        terms=((BQB.IS, "go/GO:0008015"),),  # blood circulation
    ),
    "exchange": Kind(_T, "DISTRIBUTION_PROCESS"),
    "metabolism": Kind(
        SBO.BIOCHEMICAL_REACTION,
        "LINEAR_ELIMINATION",
        terms=((BQB.IS, "go/GO:0008152"),),  # metabolic process
    ),
    "metabolism_saturable": Kind(
        SBO.BIOCHEMICAL_REACTION,
        "NON_LINEAR_ELIMINATION",
        terms=((BQB.IS, "go/GO:0008152"),),
    ),
    "gut_metabolism": Kind(
        SBO.BIOCHEMICAL_REACTION, "GUT_METABOLISM", terms=((BQB.IS, "go/GO:0008152"),)
    ),
    "transport": Kind(
        _T,
        "DISTRIBUTION_PROCESS",
        terms=((BQB.IS, "go/GO:0055085"),),  # transmembrane transport
    ),
    "biliary_excretion": Kind(
        _T, "BILIARY_ELIMINATION", terms=((BQB.IS, "go/GO:0055085"),)
    ),
    "renal_secretion": Kind(
        _T,
        "RENAL_EXCRETION",
        terms=((BQB.IS, "go/GO:0097254"),),  # renal tubular secretion
    ),
    "reabsorption": Kind(_T, "RESORPTION", terms=((BQB.IS, "go/GO:0055085"),)),
    "filtration": Kind(
        _T, "RENAL_EXCRETION", terms=((BQB.IS, "go/GO:0003094"),)
    ),  # glomerular filtration
    "tubular_flow": Kind(
        _T,
        None,
        terms=(
            (BQB.IS_PART_OF, pbpko("RENAL_EXCRETION")),
            (BQB.IS_PART_OF, "go/GO:0007588"),  # excretion
        ),
    ),  # excretion
    "absorption": Kind(
        _T, "ABSORPTION", terms=((BQB.IS, "go/GO:0050892"),)
    ),  # intestinal absorption
    "dissolution_reaction": Kind(
        _T,
        None,
        terms=((BQB.IS_PART_OF, pbpko("ABSORPTION")),),
        requested="dissolution process",
    ),
    "transit": Kind(
        _T,
        None,
        terms=((BQB.IS_VERSION_OF, "go/GO:0022600"),),  # digestive system process
        requested="gastrointestinal transit process",
    ),
    "fecal_excretion": Kind(_T, "FECAL_EXCRETION", terms=((BQB.IS, "go/GO:0007588"),)),
    "bile_flow": Kind(
        _T, None, terms=((BQB.IS_PART_OF, pbpko("ENTEROHEPATIC_RECIRCULATION")),)
    ),
    "iv_application": Kind(
        _T,
        None,
        terms=(
            (BQB.IS_VERSION_OF, pbpko("INTRAVENOUS_BOLUS_EXPOSURE")),
            (BQB.IS_VERSION_OF, pbpko("INTRAVENOUS_INFUSION_EXPOSURE")),
        ),
    ),
    "oral_application": Kind(_T, "ORAL_EXPOSURE"),
    "custom_reaction": Kind(
        SBO.BIOCHEMICAL_OR_TRANSPORT_REACTION,
        None,
        terms=((BQB.IS_VERSION_OF, "go/GO:0008150"),),  # biological process
    ),
}
"""The kinds of the elements the library creates."""


_DEFAULT_SBO = {
    "Compartment": _C,
    "Species": _SPECIES,
    "Parameter": _P,
    "Reaction": SBO.BIOCHEMICAL_OR_TRANSPORT_REACTION,
}
"""The SBO term of an element whose kind keeps the SBO term of the element, if it has none."""


def _resolve(kind: Kind, place: _Place | None, pool: Pool | None) -> str | None:
    """The PBPKO term of a kind at a place: organ-specific, of the pool or the kind's."""
    if pool is not None and pool in kind.pools:
        return kind.pools[pool]
    if place is not None and kind.organ is not None:
        for token in place.tokens:
            specific = kind.organ.format(token)
            if _has(specific):
                return specific
    return kind.pbpko


def _location(element: Sbase) -> BQB:
    """The qualifier of the place of an element."""
    if isinstance(element, Compartment):
        return BQB.IS_PART_OF
    if isinstance(element, Reaction):
        return BQB.OCCURS_IN
    return BQB.IS_PROPERTY_OF


def initial_assignment(
    symbol: str, formula: str, unit: UnitDefinition, name: str
) -> InitialAssignment:
    """Get the initial assignment of a constant, with a name and an SBO term.

    Every quantity of the library which depends on parameters only is a
    constant with an initial assignment; the assignment carries the name of
    its target and the SBO term of a mathematical expression, so the factory
    has nothing to hint about.

    Args:
        symbol: the id of the target.
        formula: the formula.
        unit: the unit of the formula.
        name: the name of the target.

    Returns:
        The initial assignment.
    """
    return InitialAssignment(
        symbol, formula, unit, name=name, sboTerm=SBO.MATHEMATICAL_EXPRESSION
    )


def annotate[S: Sbase](
    element: S,
    kind: str,
    *,
    organ: Organ | None = None,
    pool: Pool | None = None,
    segment: str | None = None,
    site: str | None = None,
    protein: str | None = None,
    substance: Substance | None = None,
) -> S:
    """Annotate an element the library creates by its kind and its place.

    Args:
        element: the compartment, species, parameter or reaction.
        kind: the kind of the element, a key of `KINDS`.
        organ: the organ of the element.
        pool: the blood pool of the element.
        segment: the lumen segment of the element (`du`, `pt`), a key of
            `SEGMENTS`.
        site: the sampling site of an observer, a key of `SITES`.
        protein: the protein of a process or of its parameters; a protein
            which is not in the protein table adds no terms (a definition
            checks its proteins, `PBPKDefinition.validate`).
        substance: the substance of a species, whose annotations it takes.

    Returns:
        The element, with the SBO term of the kind (a kind which keeps the
        SBO term of the element gives a generic one to an element without)
        and its terms added to its annotations; a term of the place whose
        resource the element has already (as what it is) is not repeated.

    Raises:
        KeyError: if the kind, the segment or the site is unknown.
    """
    k = KINDS[kind]
    place = ORGANS[organ] if organ is not None else None
    terms: list[Term] = []
    term = _resolve(k, place, pool)
    if term is not None:
        terms.append((BQB.IS, pbpko(term)))
    terms += k.terms
    if substance is not None:
        terms += [(q, str(r)) for q, r in substance.annotations]
    where = BQB.IS if k.is_place else _location(element)
    places = [p for p in (place, POOLS.get(pool) if pool else None) if p is not None]
    if segment is not None:
        places.append(SEGMENTS[segment])
    if site is not None:
        places.append(SITES[site])
    if not isinstance(element, Species):
        present = {r for _, r in terms}
        for p in places:
            terms += [(where, r) for _, r in p.terms if r not in present]
            if isinstance(element, Compartment) and p.compartment is not None:
                terms.append((where, pbpko(p.compartment)))
    if protein is not None and protein in PROTEINS:
        entry = PROTEINS[protein]
        q = BQB.IS if isinstance(element, Reaction) else BQB.IS_PROPERTY_OF
        terms.append((q, f"uniprot/{entry.uniprot}"))
        if entry.ncit is not None:
            terms.append((q, f"ncit/{entry.ncit}"))
    if k.sbo is not None:
        element.sboTerm = k.sbo
    elif element.sboTerm is None:
        element.sboTerm = _DEFAULT_SBO[type(element).__name__]
    existing = list(element.annotations or [])
    element.annotations = existing + [
        t for t in dict.fromkeys(terms) if t not in existing
    ]
    return element


def model_annotations(substances: Iterable[Substance]) -> list[Term]:
    """Get the annotations of a PBPK model: its type, the taxon and the approach.

    A model is perfusion-limited when every substance which enters tissue
    does so flow-limited (`Mode.PERFUSION`), permeability-limited when one
    has a finite exchange.

    Args:
        substances: the substances of the model.

    Returns:
        The terms of the model.
    """
    modes = {s.distribution.mode for s in substances if s.distribution is not None}
    terms: list[Term] = [
        (BQB.IS, pbpko("PHYSIOLOGICALLY_BASED_PHARMACOKINETIC_MODEL")),
        (BQB.IS, pbpko("WHOLE_BODY_PBPK")),
        (BQB.HAS_TAXON, "taxonomy/9606"),  # human
        (BQB.HAS_PROPERTY, "mamo/MAMO_0000046"),  # ordinary differential equation model
        (BQB.HAS_PROPERTY, "ncit/C15299"),  # pharmacokinetics
    ]
    if modes == {Mode.PERFUSION}:
        terms.append((BQB.IS, pbpko("PERFUSION_LIMITED_MODEL")))
    elif modes:
        terms.append((BQB.IS, pbpko("PERMEABILITY_LIMITED_PBPK")))
    return terms


def requested_terms() -> dict[str, str]:
    """Get the PBPKO terms the library needs and the pinned PBPKO lacks.

    Returns:
        The missing terms with the broader PBPKO term the library uses
        instead: the `requested` terms of the kinds and the organ-specific
        terms of the kinds which PBPKO lacks for an organ.
    """
    missing: dict[str, str] = {}
    for k in KINDS.values():
        if k.requested is not None:
            missing[k.requested] = k.pbpko or "-"
        if k.organ is None or not k.per_organ:
            continue
        for place in ORGANS.values():
            names = [k.organ.format(token) for token in place.tokens]
            if not any(_has(name) for name in names):
                missing.setdefault(names[0], k.pbpko or "-")
    for place in [*ORGANS.values(), *POOLS.values()]:
        if place.compartment is None:
            missing[f"{place.label} compartment"] = "-"
    return missing


def _check() -> None:
    """Check that every PBPKO term the table names exists in the pinned PBPKO.

    Raises:
        AttributeError: naming the first term which does not.
    """
    names = [n for k in KINDS.values() for n in (k.pbpko, *k.pools.values()) if n]
    names += [
        p.compartment
        for p in (*ORGANS.values(), *POOLS.values(), *SEGMENTS.values())
        if p.compartment
    ]
    for name in names:
        if not _has(name):
            raise AttributeError(f"The pinned PBPKO has no term '{name}'.")


_check()
