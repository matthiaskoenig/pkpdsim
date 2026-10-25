"""The id scheme of the PBPK models.

Every id of a generated model is built here, so that the organs, the
observers and the fits agree on them. A flat id is `<ORGAN>__<local id>`
(`LI__Cpl_api`), the ids of substances are suffixed by `_<substance>`.
"""

import re
from collections.abc import Iterable

from pkpdsim.pbpk.organs import Organ, Part, Pool

__all__ = [
    "AGE",
    "BMI",
    "BSA",
    "BW",
    "CKD_HEPATIC",
    "CO",
    "CO_BASE",
    "CP_SEVERITY",
    "EGFR",
    "FED",
    "FFM",
    "F_CARDIAC_FUNCTION",
    "F_FUNCTIONAL",
    "F_MUCOSA_FLOW",
    "F_MUCOSA_VOLUME",
    "F_RENAL_FUNCTION",
    "F_SHUNT",
    "F_TUBULAR",
    "F_ZONE",
    "GFR",
    "GFR_BASE",
    "GUT",
    "GUT_FLOW",
    "GUT_INTERSTITIAL_FRACTION",
    "GUT_SEGMENTS",
    "GUT_VASCULAR_FRACTION",
    "GUT_VOLUME",
    "GUT_WALL_FLOW",
    "HCT",
    "HEIGHT",
    "HF_SEVERITY",
    "K_BILE",
    "MLEAN",
    "MUCOSA_FLOW_TOTAL",
    "MUCOSA_VOLUME_TOTAL",
    "PS_VAS",
    "QH",
    "QHA",
    "QPO",
    "SV_TIS",
    "TUBULE_SEGMENTS",
    "T_LUMEN",
    "T_STOMACH_FED",
    "VBILE",
    "VBLOOD",
    "VFECES",
    "VLUMEN",
    "VURINE",
    "WALL_FRACTION",
    "absorption",
    "absorption_reaction",
    "activity",
    "area_factor",
    "arm",
    "arm_blood",
    "arm_free",
    "arm_weight",
    "bile",
    "bile_reaction",
    "blood",
    "bp",
    "compartment",
    "condition_factor",
    "dissolution",
    "dissolution_age",
    "dissolution_clock",
    "dissolution_dose50",
    "dissolution_lag",
    "dissolution_rate",
    "dissolution_reaction",
    "dissolution_shape",
    "dissolution_t50",
    "dose_reaction",
    "ec50",
    "emax",
    "enzyme",
    "exchange_reaction",
    "feces",
    "feces_reaction",
    "filtration_reaction",
    "fkp",
    "flat",
    "flow_fraction",
    "fraction_absorbed",
    "fraction_reabsorbed",
    "free",
    "ftissue",
    "fu",
    "fu_ref",
    "hepatic_vein",
    "hepatic_vein_blood",
    "hepatic_vein_free",
    "impairment_activity",
    "impairment_factor",
    "inflow_reaction",
    "infusion_rate",
    "injection_rate",
    "injection_time",
    "interstitial_fraction",
    "iv_dose",
    "k_cell_pls",
    "k_int_pls",
    "k_organ",
    "kbc",
    "kbc_exchange",
    "kdeg",
    "ki",
    "kinact",
    "kiu",
    "kp",
    "kp_organ",
    "lumen",
    "lumen_volume",
    "mr",
    "mucosa",
    "mucosa_coefficient",
    "mucosa_flow",
    "mucosa_reaction",
    "mucosa_unbound",
    "mucosa_volume",
    "organ_flow",
    "organ_mass",
    "organ_outflow",
    "organ_outflow_bc",
    "organ_volume",
    "outflow_conc",
    "outflow_conc_bc",
    "outflow_reaction",
    "placement_profile",
    "placement_weight",
    "plasma_capacity",
    "po_dose",
    "pool_capacity",
    "pool_compartment",
    "pool_exchange_reaction",
    "pool_fraction",
    "pool_species",
    "pool_volume",
    "ps",
    "ref",
    "segment_process",
    "shunt_reaction",
    "sinusoid_reaction",
    "sol_dose",
    "solid",
    "solution",
    "species",
    "species_of",
    "submodel",
    "tissue_volume_total_ref",
    "transit_reaction",
    "transit_time",
    "tubule",
    "tubule_flow_fraction",
    "tubule_reabsorption_reaction",
    "tubule_reaction",
    "tubule_volume",
    "unbound",
    "urine",
    "urine_rate",
    "vascular_fraction",
    "wall_compartment",
    "zone",
]


BW = "BW"
"""The body weight (covariate)."""
HEIGHT = "HEIGHT"
"""The height (covariate)."""
AGE = "AGE"
"""The age (covariate)."""
BMI = "BMI"
"""The body mass index."""
BSA = "BSA"
"""The body surface area."""
FFM = "FFM"
"""The fat-free mass."""
HCT = "HCT"
"""The hematocrit."""
CO = "CO"
"""The cardiac output."""
CO_BASE = "CO_base"
"""The cardiac output without impairment, of which the organ flows are fractions."""
GFR = "GFR"
"""The glomerular filtration rate."""
GFR_BASE = "GFR_base"
"""The glomerular filtration rate at renal function 1, of which the renal function is the factor."""
VBLOOD = "Vblood"
"""The total blood volume."""
MLEAN = "Mlean"
"""The mass of the lean organ tissues without their blood."""
SV_TIS = "SV_tis"
"""The scaling factor of the lean tissues, the blood and the cardiac output."""
QHA = "Qha"
"""The hepatic arterial blood flow."""
QPO = "Qpo"
"""The portal venous blood flow."""
QH = "Qh"
"""The total liver blood flow."""
F_RENAL_FUNCTION = "f_renal_function"
"""The renal function relative to normal."""
F_CARDIAC_FUNCTION = "f_cardiac_function"
"""The cardiac function (cardiac output) relative to normal."""
F_SHUNT = "f_shunt"
"""The portosystemic shunt: the fraction of the liver inflow which bypasses the liver."""
F_FUNCTIONAL = "f_functional"
"""The functional fraction of the liver tissue."""
F_TUBULAR = "f_tubular"
"""The share of the tubular capacity of the kidney which follows the renal function."""
CP_SEVERITY = "CP_severity"
"""The severity of liver cirrhosis: 0 healthy, 1, 2, 3 the Child-Pugh classes A, B, C."""
HF_SEVERITY = "HF_severity"
"""The severity of heart failure: 0 healthy, 1 to 4 the NYHA classes I to IV."""
CKD_HEPATIC = "CKD_hepatic"
"""The switch (0 off, 1 on) of the hepatic changes of chronic kidney disease."""
EGFR = "eGFR"
"""The GFR per 1.73 m2 body surface area of the renal function, which gives the CKD stage."""
F_ZONE = "f_zone"
"""The volume fraction of a zone of a zonated organ, `1/n` of `n` zones."""
VBILE = "Vbile"
"""The volume of the bile, the compartment of `Cbile_<s>`."""
VURINE = "Vurine"
"""The compartment of the urine amounts `Aurine_<s>`."""
VLUMEN = "Vlumen"
"""The local volume of the single intestinal lumen, the compartment of `Clumen_<s>`."""
VFECES = "Vfeces"
"""The compartment of the feces amounts `Afeces_<s>`."""
K_BILE = "k_bile"
"""The rate of the transit of the bile into the intestinal lumen."""
T_LUMEN = "t_lumen"
"""The local transit time of the single intestinal lumen into the feces."""
FED = "FED"
"""The prandial state, 0 fasted (default) and 1 fed: it slows the gastric emptying."""
GUT = "GU"
"""The id of the submodel of the intestine (the small intestine, or the
segmented and single lumen intestine, which replace small and large
intestine), the prefix of the hand-written models `GU__`."""
GUT_VOLUME = "Vgu"
"""The local volume of the intestine which replaces small and large intestine, `Vsi + Vla`."""
GUT_VASCULAR_FRACTION = "f_vas_gu"
"""The local vascular fraction of the intestine of small and large intestine."""
GUT_INTERSTITIAL_FRACTION = "f_int_gu"
"""The local interstitial fraction of the intestine of small and large intestine."""
GUT_FLOW = "Qgu"
"""The local blood flow of the intestine of small and large intestine, `Qsi + Qla`."""
GUT_WALL_FLOW = "Qgu_wall"
"""The local blood flow of the intestinal wall past the mucosa, `Qgu - Qmu_total`."""
MUCOSA_FLOW_TOTAL = "Qmu_total"
"""The local blood flow of the mucosa of all segments."""
MUCOSA_VOLUME_TOTAL = "Vmu_total"
"""The local volume of the mucosa of all segments."""
WALL_FRACTION = "f_wall"
"""The local fraction of the intestinal tissue without the mucosa, `1 - Vmu_total/Vti`."""
F_MUCOSA_VOLUME = "f_mucosa_volume"
"""The local fraction of the intestinal tissue which is mucosa."""
F_MUCOSA_FLOW = "f_mucosa_flow"
"""The local fraction of the intestinal blood flow which perfuses the mucosa."""
T_STOMACH_FED = "T_st_fed"
"""The local transit time of the stomach in the fed state."""
GUT_SEGMENTS: tuple[str, ...] = (
    "st",
    "du",
    "jej1",
    "jej2",
    "ile1",
    "ile2",
    "cae",
    "asc",
    "tra",
    "des",
    "sig",
    "rec",
)
"""The lumen segments of the segmented intestine, in the order of the transit:
stomach, duodenum, upper and lower jejunum, upper and lower ileum, caecum,
colon ascendens, transversum, descendens, sigmoid and rectum."""
TUBULE_SEGMENTS: tuple[str, ...] = ("pt", "lh", "dt", "cd")
"""The tubule segments of the nephron kidney, in the order of the tubular
flow: proximal tubule, loop of Henle, distal tubule and collecting duct."""
PS_VAS = "PS_vas"
"""The local id of the exchange between plasma and interstitium of an organ."""


def impairment_factor(pid: str) -> str:
    """The id of the factor of all impairments on a quantity, `fimp_GFR`.

    The product of the factors of the conditions (Child-Pugh, chronic kidney
    disease, heart failure) on the quantity `pid` of the body.
    """
    return f"fimp_{pid}"


def condition_factor(condition: str, key: str) -> str:
    """The id of the factor of one condition, `fcp_albumin`, `fckd_hematocrit`.

    Args:
        condition: the condition, `cp` (Child-Pugh), `ckd` (chronic kidney
            disease) or `hf` (heart failure).
        key: the key of the factor in `data/impairment.toml`, a protein in
            lower case.
    """
    return f"f{condition}_{key.lower()}"


def impairment_activity(protein: str, organ: Organ | None = None) -> str:
    """The id of the factor of the impairments on the activity of a protein in an organ.

    In the organ the local input `fimp_cyp3a4`, in the body `fimp_li_cyp3a4`;
    the rates of the protein in the organ run with `f_<protein>` times it, so
    that a genotype or an interaction which sets `LI__f_cyp3a4` acts on top
    of the impairment.
    """
    if organ is None:
        return f"fimp_{protein.lower()}"
    return f"fimp_{organ.value}_{protein.lower()}"


def ref(pid: str) -> str:
    """The id of the reference value of a parameter, `FFM_ref`, `Vli_ref`."""
    return f"{pid}_ref"


def submodel(organ: Organ) -> str:
    """The id of the submodel of an organ, the organ id in capitals, `LI`.

    The small intestine is the submodel `GU` (gut) for every intestine
    variant: the segmented and the single lumen intestine replace small and
    large intestine with it; with the standard organs as intestine the small
    intestine is `GU` and the large intestine `LA`.
    """
    if organ is Organ.SMALL_INTESTINE:
        return GUT
    return organ.value.upper()


def flat(organ: Organ, local: str) -> str:
    """The id of a local id in the flattened model, `LI__<local>`."""
    return f"{submodel(organ)}__{local}"


def zone(i: int) -> str:
    """The token of zone `i` of a zonated organ, `z1` (periportal) to `z<n>`."""
    return f"z{i}"


def _zoned(prefix: str, zone: str | None, s: str | None = None) -> str:
    """An id of an organ part in a zone, `<prefix>_<zone>_<s>`, without zone `<prefix>_<s>`."""
    return "_".join(p for p in (prefix, zone, s) if p is not None)


def compartment(part: Part, zone: str | None = None) -> str:
    """The local id of a compartment of an organ, `Vpl`; of a zone `Vpl_z1`."""
    return _zoned(f"V{part.value}", zone)


def species(part: Part, s: str, zone: str | None = None) -> str:
    """The local id of a species of substance `s` in a part, `Cpl_<s>`; of a zone `Cpl_z1_<s>`."""
    return _zoned(f"C{part.value}", zone, s)


def pool_compartment(pool: Pool, part: Part) -> str:
    """The id of a compartment of a blood pool, `Vve_pl`."""
    return f"V{pool.value}_{part.value}"


def pool_species(pool: Pool, s: str, part: Part = Part.PLASMA) -> str:
    """The id of a species in a blood pool, `Cve_<s>` (plasma), `Cve_bc_<s>` (blood cells)."""
    if part is Part.PLASMA:
        return f"C{pool.value}_{s}"
    return f"C{pool.value}_{part.value}_{s}"


def organ_volume(organ: Organ) -> str:
    """The id of the volume of an organ, `Vli`."""
    return f"V{organ.value}"


def organ_flow(organ: Organ) -> str:
    """The id of the blood flow of an organ, `Qli`."""
    return f"Q{organ.value}"


def organ_mass(organ: Organ) -> str:
    """The id of the tissue mass of an organ without its blood, `Mad`."""
    return f"M{organ.value}"


def flow_fraction(organ: Organ) -> str:
    """The id of the blood flow of an organ as fraction of cardiac output, `fQli`."""
    return f"fQ{organ.value}"


def pool_volume(pool: Pool) -> str:
    """The id of the volume of a blood pool, `Vve`."""
    return f"V{pool.value}"


def pool_fraction(pool: Pool) -> str:
    """The id of the volume of a blood pool as fraction of the blood volume, `fVve`."""
    return f"fV{pool.value}"


def pool_capacity(pool: Pool, s: str) -> str:
    """The id of the compartment of a pool species of `s` with erythrocytes in equilibrium.

    The pool volume times the blood to plasma ratio, `Vve_pl_<s>`: the plasma
    species then holds the amount of the blood.
    """
    return f"V{pool.value}_{Part.PLASMA.value}_{s}"


def plasma_capacity(s: str, zone: str | None = None) -> str:
    """The local id of the compartment of `Cpl_<s>` with erythrocytes in equilibrium.

    The vascular volume of the organ times the blood to plasma ratio,
    `Vpl_<s>` (of a zone `Vpl_z1_<s>`): the plasma species then holds the
    amount of the blood.
    """
    return _zoned(f"V{Part.PLASMA.value}", zone, s)


def vascular_fraction(organ: Organ) -> str:
    """The id of the vascular fraction of the organ volume, `f_vas_li`."""
    return f"f_vas_{organ.value}"


def interstitial_fraction(organ: Organ) -> str:
    """The id of the interstitial fraction of the organ volume, `f_int_li`."""
    return f"f_int_{organ.value}"


def tissue_volume_total_ref(s: str) -> str:
    """The id of the reference tissue volume of the organs `s` enters, `Vti_total_ref_<s>`.

    The tissue volume of the reference individual summed over the organs
    whose exchange is a share of `ftissue_<s>`; the exchange of an organ is
    `ftissue_<s> * Vti_o / Vti_total_ref_<s>`, so it scales with the individual.
    """
    return f"Vti_total_ref_{s}"


def ps(s: str) -> str:
    """The local id of the exchange of `s` between plasma and tissue, `PS_<s>`."""
    return f"PS_{s}"


def k_int_pls(s: str) -> str:
    """The local id of the interstitial to plasma coefficient, `K_int_pls_<s>`."""
    return f"K_int_pls_{s}"


def k_cell_pls(s: str) -> str:
    """The local id of the intracellular to plasma coefficient, `K_cell_pls_<s>`."""
    return f"K_cell_pls_{s}"


def k_organ(s: str) -> str:
    """The local id of the organ to plasma coefficient (flow-limited), `K_or_<s>`."""
    return f"K_or_{s}"


def organ_outflow(organ: Organ, s: str) -> str:
    """The id of the plasma concentration leaving an organ, `Cout_li_<s>`."""
    return f"Cout_{organ.value}_{s}"


def organ_outflow_bc(organ: Organ, s: str) -> str:
    """The id of the blood cell concentration leaving an organ, `Cout_li_bc_<s>`.

    Only for a substance with kinetic blood cells.
    """
    return f"Cout_{organ.value}_{Part.BLOOD_CELLS.value}_{s}"


def inflow_reaction(pool: Pool, part: Part, s: str) -> str:
    """The local id of the inflow of `s` from a pool into an organ, `J_in_ar_pl_<s>`.

    `part` is the part of the blood which carries it, plasma or blood cells.
    """
    return f"J_in_{pool.value}_{part.value}_{s}"


def outflow_reaction(part: Part, s: str) -> str:
    """The local id of the outflow of `s` from an organ into its pool, `J_out_pl_<s>`."""
    return f"J_out_{part.value}_{s}"


def shunt_reaction(pool: Pool, part: Part, s: str) -> str:
    """The local id of the shunt of `s` from a pool past the liver, `J_shunt_po_pl_<s>`.

    `part` is the part of the blood which carries it, plasma or blood cells.
    """
    return f"J_shunt_{pool.value}_{part.value}_{s}"


def filtration_reaction(s: str) -> str:
    """The local id of the glomerular filtration of `s` into the urine, `J_gfr_<s>`."""
    return f"J_gfr_{s}"


def exchange_reaction(part: Part, s: str, zone: str | None = None) -> str:
    """The local id of the passive exchange of `s` into a part of an organ, `J_ti_<s>`.

    From plasma into the tissue, the interstitium or the blood cells, from
    the interstitium into the cells (`J_ic_<s>`); in a zone `J_ti_z1_<s>`.
    """
    return _zoned(f"J_{part.value}", zone, s)


def sinusoid_reaction(part: Part, s: str, zone: str) -> str:
    """The local id of the flow of `s` from a zone into the next one, `J_sin_z1_pl_<s>`.

    `part` is the part of the blood which carries it, plasma or blood cells
    (the whole organ of a flow-limited substance carries it as plasma).
    """
    return f"J_sin_{zone}_{part.value}_{s}"


def outflow_conc(s: str) -> str:
    """The local id of the concentration leaving an organ, `Cout_<s>`."""
    return f"Cout_{s}"


def outflow_conc_bc(s: str) -> str:
    """The local id of the blood cell concentration leaving an organ, `Cout_bc_<s>`."""
    return f"Cout_{Part.BLOOD_CELLS.value}_{s}"


def mr(s: str) -> str:
    """The id of the molecular weight, `Mr_<s>`."""
    return f"Mr_{s}"


def fu(s: str) -> str:
    """The id of the unbound fraction in plasma, `fu_<s>`."""
    return f"fu_{s}"


def fu_ref(s: str) -> str:
    """The id of the reference unbound fraction in plasma, `fu_ref_<s>`."""
    return f"fu_ref_{s}"


def kp(s: str) -> str:
    """The id of the tissue to plasma coefficient, `Kp_<s>`."""
    return f"Kp_{s}"


def kp_organ(organ: Organ, s: str) -> str:
    """The id of the tissue to plasma coefficient of an organ, `Kp_li_<s>`."""
    return f"Kp_{organ.value}_{s}"


def fkp(s: str) -> str:
    """The id of the factor of the coefficients of a Kp method, `fKp_<s>`."""
    return f"fKp_{s}"


def ftissue(s: str) -> str:
    """The id of the whole-body exchange factor, `ftissue_<s>`."""
    return f"ftissue_{s}"


def bp(s: str) -> str:
    """The id of the blood to plasma ratio, `BP_<s>`."""
    return f"BP_{s}"


def kbc(s: str) -> str:
    """The id of the blood cell to plasma coefficient, `Kbc_<s>`."""
    return f"Kbc_{s}"


def kbc_exchange(s: str) -> str:
    """The id of the exchange rate between plasma and blood cells, `kex_bc_<s>`."""
    return f"kex_bc_{s}"


def fraction_reabsorbed(s: str) -> str:
    """The id of the fraction of the filtrate which is reabsorbed, `f_reabsorbed_<s>`."""
    return f"f_reabsorbed_{s}"


def urine(s: str) -> str:
    """The id of the amount in urine, `Aurine_<s>`."""
    return f"Aurine_{s}"


def urine_rate(s: str) -> str:
    """The id of the urinary excretion rate, `Rurine_<s>` (mmole/min).

    The sum of the fluxes of the kidney into the urine `Aurine_<s>`, an output
    of the kidney which replaces the parameter of the top model.
    """
    return f"Rurine_{s}"


def feces(s: str) -> str:
    """The id of the amount in feces, `Afeces_<s>`."""
    return f"Afeces_{s}"


def bile(s: str) -> str:
    """The id of the concentration in bile, `Cbile_<s>`."""
    return f"Cbile_{s}"


def _segment(prefix: str, s: str | None, seg: str | None) -> str:
    """An id of a lumen segment, `<prefix>_<seg>_<s>`, without segment `<prefix>_<s>`."""
    parts = [prefix, *(p for p in (seg, s) if p is not None)]
    return "_".join(parts)


def lumen(s: str, seg: str | None = None) -> str:
    """The local id of the concentration in an intestinal lumen segment, `Clumen_du_<s>`.

    Without a segment the single lumen, `Clumen_<s>`.
    """
    return _segment("Clumen", s, seg)


def lumen_volume(seg: str | None = None) -> str:
    """The local id of the volume of a lumen segment, `Vlumen_du` (single: `Vlumen`)."""
    return _segment("Vlumen", None, seg)


def solid(s: str, seg: str) -> str:
    """The local id of the undissolved amount of `s` in a lumen segment, `Asolid_du_<s>`.

    In the stomach (and the single lumen) the dose state `PODOSE_<s>` is the
    undissolved drug.
    """
    return _segment("Asolid", s, seg)


def tubule(s: str, seg: str) -> str:
    """The local id of the concentration in the lumen of a tubule segment, `Ctub_pt_<s>`.

    The segment is always part of the id (`TUBULE_SEGMENTS`), so the id
    names its substance without ambiguity.
    """
    return f"Ctub_{seg}_{s}"


def tubule_volume(seg: str) -> str:
    """The local id of the volume of the lumen of a tubule segment, `Vtub_pt`."""
    return f"Vtub_{seg}"


def tubule_flow_fraction(seg: str) -> str:
    """The local id of the tubular flow leaving a segment as fraction of the GFR, `f_tub_pt`."""
    return f"f_tub_{seg}"


def tubule_reaction(s: str, seg: str) -> str:
    """The local id of the tubular flow of `s` out of a segment, `J_tub_pt_<s>`.

    Out of the collecting duct the urinary excretion `J_tub_cd_<s>`.
    """
    return f"J_tub_{seg}_{s}"


def tubule_reabsorption_reaction(s: str) -> str:
    """The local id of the reabsorbed fraction of `s` out of the collecting duct, `J_reab_<s>`."""
    return f"J_reab_{s}"


def mucosa(s: str, seg: str | None = None) -> str:
    """The local id of the concentration in the mucosa of a segment, `Cmu_du_<s>`."""
    return _segment("Cmu", s, seg)


def mucosa_unbound(s: str, seg: str | None = None) -> str:
    """The local id of the unbound concentration in the mucosa of a segment, `Cu_mu_du_<s>`."""
    return _segment("Cu_mu", s, seg)


def wall_compartment(part: Part) -> str:
    """The local id of a tissue compartment of the intestinal wall without the mucosa, `Vti_wall`.

    The compartment of the tissue species of a substance with mucosa species:
    the mucosa is carved out of the wall, so the capacity of the intestine
    stays that of the wall without mucosa.
    """
    return f"V{part.value}_wall"


def mucosa_volume(seg: str | None = None) -> str:
    """The local id of the volume of the mucosa of a segment, `Vmu_du`."""
    return _segment("Vmu", None, seg)


def mucosa_flow(seg: str | None = None) -> str:
    """The local id of the blood flow of the mucosa of a segment, `Qmu_du`."""
    return _segment("Qmu", None, seg)


def mucosa_coefficient(s: str) -> str:
    """The local id of the mucosa to plasma coefficient, `K_mu_<s>`."""
    return f"K_mu_{s}"


def transit_time(seg: str | None = None) -> str:
    """The local id of the transit time of a lumen segment, `T_du` (single: `t_lumen`)."""
    return T_LUMEN if seg is None else f"T_{seg}"


def area_factor(seg: str) -> str:
    """The local id of the absorption scaling of a lumen segment, `f_area_du`."""
    return f"f_area_{seg}"


def fraction_absorbed(s: str) -> str:
    """The id of the fraction of the absorbed flux entering the gut, `fabs_<s>` (single lumen)."""
    return f"fabs_{s}"


def transit_reaction(s: str, seg: str | None = None, solid: bool = False) -> str:
    """The local id of the transit of `s` out of a lumen segment, `J_tr_du_<s>`.

    The undissolved drug `J_trs_du_<s>`.
    """
    return _segment("J_trs" if solid else "J_tr", s, seg)


def dissolution_reaction(s: str, seg: str) -> str:
    """The local id of the dissolution of `s` in a lumen segment, `J_dis_du_<s>`."""
    return _segment("J_dis", s, seg)


def mucosa_reaction(kind: str, s: str, seg: str | None = None) -> str:
    """The local id of a flow of the mucosa of a segment, `J_mu_in_du_<s>`.

    `kind` is `in` (arterial inflow), `out` (outflow into the portal vein)
    or `bc` (the blood cells of the mucosal blood, which pass it).
    """
    return _segment(f"J_mu_{kind}", s, seg)


def placement_weight(key: str, seg: str | None) -> str:
    """The local id of the share of a process in a segment or a zone, `w_cyp3a4_du`.

    `key` is the protein (or profile) in lower case or `uniform` (by mucosa
    volume); in the zone of a zonated organ `w_cyp3a4_z1`.
    """
    return _segment(f"w_{key}", None, seg)


def placement_profile(key: str, zone: str) -> str:
    """The local id of the relative abundance of a profile in a zone, `p_cyp3a4_z1`."""
    return f"p_{key}_{zone}"


def segment_process(pid: str, seg: str | None) -> str:
    """The local id of a process of the chemistry in a segment or a zone, `CYP3A4_DRG_du`.

    In the mucosa of a segment of the intestine or a zone of a zonated organ
    (`CYP3A4_DRG_z1`).
    """
    return pid if seg is None else f"{pid}_{seg}"


def dissolution_rate(s: str) -> str:
    """The local id of the dissolution rate of a tablet (1/min), `kdis_<s>`.

    An assignment rule of the intestine for a Weibull or dose-dependent
    dissolution; a first order dissolution is the constant `Ka_dis_<s>`.
    """
    return f"kdis_{s}"


def dissolution_age(s: str) -> str:
    """The local id of the mean age of the undissolved dose (min), `age_dis_<s>`."""
    return f"age_dis_{s}"


def dissolution_clock(s: str) -> str:
    """The local id of the age-weighted undissolved amount (mmole*min), `W_dis_<s>`.

    The state whose ratio to the undissolved amount is the mean time since
    the dose (`age_dis_<s>`) of a Weibull dissolution.
    """
    return f"W_dis_{s}"


def iv_dose(s: str) -> str:
    """The id of the intravenous dose state (mg), `IVDOSE_<s>`."""
    return f"IVDOSE_{s}"


def injection_time(s: str) -> str:
    """The id of the half-life of an intravenous injection (s), `ti_<s>`."""
    return f"ti_{s}"


def injection_rate(s: str) -> str:
    """The id of the rate of an intravenous injection, `Kiv_<s>`."""
    return f"Kiv_{s}"


def infusion_rate(s: str) -> str:
    """The id of the rate of an intravenous infusion (mg/min), `Ri_<s>`."""
    return f"Ri_{s}"


def po_dose(s: str) -> str:
    """The id of the oral dose state of a tablet (mg), `PODOSE_<s>`."""
    return f"PODOSE_{s}"


def sol_dose(s: str) -> str:
    """The id of the oral dose state of a solution (mg), `SOLDOSE_<s>`."""
    return f"SOLDOSE_{s}"


def dissolution(s: str) -> str:
    """The id of the rate of the dissolution of a tablet, `Ka_dis_<s>`."""
    return f"Ka_dis_{s}"


def dissolution_t50(s: str) -> str:
    """The id of the time of half dissolution of a Weibull dissolution (min), `t50_dis_<s>`."""
    return f"t50_dis_{s}"


def dissolution_shape(s: str) -> str:
    """The id of the shape of a Weibull dissolution, `shape_dis_<s>`."""
    return f"shape_dis_{s}"


def dissolution_lag(s: str) -> str:
    """The id of the lag time of a Weibull dissolution (min), `lag_dis_<s>`."""
    return f"lag_dis_{s}"


def dissolution_dose50(s: str) -> str:
    """The id of the undissolved amount of half the dissolution rate (mg), `D50_dis_<s>`."""
    return f"D50_dis_{s}"


def solution(s: str) -> str:
    """The id of the rate of the application of a solution, `Ka_sol_<s>`."""
    return f"Ka_sol_{s}"


def absorption(s: str) -> str:
    """The id of the rate of the absorption from the intestinal lumen, `ka_<s>`.

    The first order absorption of the single lumen; in the segmented
    intestine the absorption rate constant of the small intestine, which the
    area factor of a segment scales.
    """
    return f"ka_{s}"


def dose_reaction(kind: str, s: str) -> str:
    """The id of the application of a dose state, `J_iv_<s>`, `J_po_<s>`, `J_sol_<s>`."""
    return f"J_{kind}_{s}"


def absorption_reaction(s: str, seg: str | None = None) -> str:
    """The local id of the absorption from a lumen segment, `J_abs_du_<s>`."""
    return _segment("J_abs", s, seg)


def feces_reaction(s: str) -> str:
    """The id of the unabsorbed fraction of the absorption into the feces, `J_feces_<s>`.

    Of the single lumen with `fabs_<s>` below 1.
    """
    return f"J_feces_{s}"


def bile_reaction(s: str) -> str:
    """The id of the transit of the bile into the intestinal lumen, `J_bile_<s>`.

    Into the feces when the intestine has no lumen (standard organs).
    """
    return f"J_bile_{s}"


def pool_exchange_reaction(pool: Pool, s: str) -> str:
    """The id of the exchange of `s` between plasma and blood cells of a pool, `J_bc_ve_<s>`."""
    return f"J_{Part.BLOOD_CELLS.value}_{pool.value}_{s}"


def arm(s: str) -> str:
    """The id of the plasma concentration at the arm-vein site, `Carm_<s>`."""
    return f"Carm_{s}"


def arm_free(s: str) -> str:
    """The id of the unbound concentration at the arm-vein site, `Carm_free_<s>`."""
    return f"Carm_free_{s}"


def arm_blood(s: str) -> str:
    """The id of the blood concentration at the arm-vein site, `Carm_blood_<s>`."""
    return f"Carm_blood_{s}"


def arm_weight(site: Organ | Pool) -> str:
    """The id of the weight of an organ outflow or a pool in the arm vein, `w_arm_sk`."""
    return f"w_arm_{site.value}"


def free(pool: Pool, s: str) -> str:
    """The id of the unbound concentration in a pool, `Cve_free_<s>`."""
    return f"C{pool.value}_free_{s}"


def blood(pool: Pool, s: str) -> str:
    """The id of the blood concentration in a pool, `Cve_blood_<s>`."""
    return f"C{pool.value}_blood_{s}"


def hepatic_vein(s: str) -> str:
    """The id of the concentration in the hepatic vein, `Chv_<s>`."""
    return f"Chv_{s}"


def hepatic_vein_free(s: str) -> str:
    """The id of the unbound concentration in the hepatic vein, `Chv_free_<s>`."""
    return f"Chv_free_{s}"


def hepatic_vein_blood(s: str) -> str:
    """The id of the blood concentration in the hepatic vein, `Chv_blood_<s>`."""
    return f"Chv_blood_{s}"


def unbound(part: Part, s: str, zone: str | None = None) -> str:
    """The local id of the unbound concentration of `s` in a tissue part, `Cu_ti_<s>`.

    In a zone of a zonated organ `Cu_ti_z1_<s>`.
    """
    return _zoned(f"Cu_{part.value}", zone, s)


def activity(protein: str) -> str:
    """The id of the base activity of a protein in an organ, `f_cyp3a4`.

    A constant of every organ in which a process uses the protein (flat
    `LI__f_cyp3a4`), relative to normal (1): the factor which genotype,
    impairment and the experiments set. Interactions act on top of it, not on
    it: mechanism-based inhibition and induction through the turnover state
    `E_<protein>` of the organ, reversible inhibition through the inhibitor
    terms in the rate of every process of the protein. A process of the
    protein runs with `f_<protein>` (times `E_<protein>`) on its Vmax or
    clearance.
    """
    return f"f_{protein.lower()}"


def enzyme(protein: str) -> str:
    """The id of the turnover state of a protein in an organ, `E_cyp3a4`.

    The amount of the protein relative to its baseline (initially 1), a state
    of the organs with a process of the protein when a substance of the
    definition inactivates (mechanism-based) or induces it.
    """
    return f"E_{protein.lower()}"


def kdeg(protein: str) -> str:
    """The id of the degradation rate of a protein in an organ, `kdeg_cyp3a4`."""
    return f"kdeg_{protein.lower()}"


def ki(protein: str, s: str) -> str:
    """The id of the inhibition constant of `s` on a protein, `Ki_cyp3a4_<s>`.

    The competitive constant of a mixed inhibition, the KI (half-maximal
    inactivation) of a mechanism-based inhibition.
    """
    return f"Ki_{protein.lower()}_{s}"


def kiu(protein: str, s: str) -> str:
    """The id of the uncompetitive constant of a mixed inhibition, `Kiu_cyp3a4_<s>`."""
    return f"Kiu_{protein.lower()}_{s}"


def kinact(protein: str, s: str) -> str:
    """The id of the maximal inactivation rate of `s` on a protein, `kinact_cyp3a4_<s>`."""
    return f"kinact_{protein.lower()}_{s}"


def emax(protein: str, s: str) -> str:
    """The id of the maximal induction of a protein by `s`, `Emax_cyp3a4_<s>`."""
    return f"Emax_{protein.lower()}_{s}"


def ec50(protein: str, s: str) -> str:
    """The id of the half-maximal inducing concentration of `s`, `EC50_cyp3a4_<s>`."""
    return f"EC50_{protein.lower()}_{s}"


_PARTS = "|".join(p.value for p in Part)
_POOLS = "|".join(p.value for p in Pool)
_SEGMENTS = "|".join(GUT_SEGMENTS)
_TUBULE = "|".join(TUBULE_SEGMENTS)


def species_of(s: str, ids: Iterable[str]) -> list[str]:
    """The flat species ids of substance `s`, by the id scheme.

    Only species of exactly this substance match: `api` does not match the
    species of `api_m2`, `drg` not those of `xdrg`.

    Args:
        s: the substance id.
        ids: the ids to select from.

    Returns:
        The ids which are an organ species (of a zone as well), a pool
        species, a blood cell pool species, a lumen, undissolved or mucosa
        species, a tubule species of `s`, in the order of `ids`.
    """
    pattern = re.compile(
        rf"(?:C(?:(?:{_PARTS})(?:_z\d+)?|(?:{_POOLS})(?:_{Part.BLOOD_CELLS.value})?)"
        rf"|(?:Clumen|Cmu|Asolid)(?:_(?:{_SEGMENTS}))?"
        rf"|Ctub_(?:{_TUBULE}))_{re.escape(s)}"
    )
    return [i for i in ids if pattern.fullmatch(i.rsplit("__", 1)[-1])]
