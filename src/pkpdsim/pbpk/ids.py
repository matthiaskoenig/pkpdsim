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
    "CO",
    "FFM",
    "F_CARDIAC_FUNCTION",
    "F_FUNCTIONAL",
    "F_RENAL_FUNCTION",
    "F_SHUNT",
    "GFR",
    "HCT",
    "HEIGHT",
    "K_BILE",
    "K_FECES",
    "MLEAN",
    "PS_VAS",
    "QH",
    "QHA",
    "QPO",
    "SV_TIS",
    "T_LUMEN",
    "VBILE",
    "VBLOOD",
    "VFECES",
    "VLUMEN",
    "VURINE",
    "absorption",
    "absorption_reaction",
    "activity",
    "arm",
    "arm_blood",
    "arm_free",
    "arm_weight",
    "bile",
    "bile_reaction",
    "blood",
    "bp",
    "compartment",
    "dissolution",
    "dose_reaction",
    "exchange_reaction",
    "feces",
    "feces_reaction",
    "filtration_reaction",
    "flat",
    "flow_fraction",
    "fraction_reabsorbed",
    "free",
    "ftissue",
    "fu",
    "fu_ref",
    "hepatic_vein",
    "hepatic_vein_blood",
    "hepatic_vein_free",
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
    "kp",
    "kp_organ",
    "lumen",
    "mr",
    "organ_flow",
    "organ_mass",
    "organ_outflow",
    "organ_outflow_bc",
    "organ_volume",
    "outflow_conc",
    "outflow_conc_bc",
    "outflow_reaction",
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
    "shunt_reaction",
    "sol_dose",
    "solution",
    "species",
    "species_of",
    "submodel",
    "tissue_volume_total_ref",
    "unbound",
    "urine",
    "vascular_fraction",
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
GFR = "GFR"
"""The glomerular filtration rate."""
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
VBILE = "Vbile"
"""The volume of the bile, the compartment of `Cbile_<s>`."""
VURINE = "Vurine"
"""The compartment of the urine amounts `Aurine_<s>`."""
VLUMEN = "Vlumen"
"""The volume of the intestinal lumen, the compartment of `Clumen_<s>`."""
VFECES = "Vfeces"
"""The compartment of the feces amounts `Afeces_<s>`."""
K_BILE = "k_bile"
"""The rate of the transit of the bile into the intestinal lumen."""
K_FECES = "k_feces"
"""The rate of the transit of the intestinal lumen into the feces, `1/t_lumen`."""
T_LUMEN = "t_lumen"
"""The transit time of the intestinal lumen."""
PS_VAS = "PS_vas"
"""The local id of the exchange between plasma and interstitium of an organ."""


def ref(pid: str) -> str:
    """The id of the reference value of a parameter, `FFM_ref`, `Vli_ref`."""
    return f"{pid}_ref"


def submodel(organ: Organ) -> str:
    """The id of the submodel of an organ, the organ id in capitals, `LI`."""
    return organ.value.upper()


def flat(organ: Organ, local: str) -> str:
    """The id of a local id in the flattened model, `LI__<local>`."""
    return f"{submodel(organ)}__{local}"


def compartment(part: Part) -> str:
    """The local id of a compartment of an organ, `Vpl`."""
    return f"V{part.value}"


def species(part: Part, s: str) -> str:
    """The local id of a species of substance `s` in a part, `Cpl_<s>`."""
    return f"C{part.value}_{s}"


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


def plasma_capacity(s: str) -> str:
    """The local id of the compartment of `Cpl_<s>` with erythrocytes in equilibrium.

    The vascular volume of the organ times the blood to plasma ratio,
    `Vpl_<s>`: the plasma species then holds the amount of the blood.
    """
    return f"V{Part.PLASMA.value}_{s}"


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


def exchange_reaction(part: Part, s: str) -> str:
    """The local id of the passive exchange of `s` into a part of an organ, `J_ti_<s>`.

    From plasma into the tissue, the interstitium or the blood cells, from
    the interstitium into the cells (`J_ic_<s>`).
    """
    return f"J_{part.value}_{s}"


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


def feces(s: str) -> str:
    """The id of the amount in feces, `Afeces_<s>`."""
    return f"Afeces_{s}"


def bile(s: str) -> str:
    """The id of the concentration in bile, `Cbile_<s>`."""
    return f"Cbile_{s}"


def lumen(s: str) -> str:
    """The id of the concentration in the intestinal lumen, `Clumen_<s>`."""
    return f"Clumen_{s}"


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


def solution(s: str) -> str:
    """The id of the rate of the application of a solution, `Ka_sol_<s>`."""
    return f"Ka_sol_{s}"


def absorption(s: str) -> str:
    """The id of the rate of the absorption from the intestinal lumen, `ka_<s>`."""
    return f"ka_{s}"


def dose_reaction(kind: str, s: str) -> str:
    """The id of the application of a dose state, `J_iv_<s>`, `J_po_<s>`, `J_sol_<s>`."""
    return f"J_{kind}_{s}"


def absorption_reaction(s: str) -> str:
    """The id of the absorption from the intestinal lumen, `J_abs_<s>`."""
    return f"J_abs_{s}"


def feces_reaction(s: str) -> str:
    """The id of the transit of the intestinal lumen into the feces, `J_feces_<s>`."""
    return f"J_feces_{s}"


def bile_reaction(s: str) -> str:
    """The id of the transit of the bile into the intestinal lumen, `J_bile_<s>`."""
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


def unbound(part: Part, s: str) -> str:
    """The local id of the unbound concentration of `s` in a tissue part, `Cu_ti_<s>`."""
    return f"Cu_{part.value}_{s}"


def activity(protein: str) -> str:
    """The id of the relative activity of a protein, `f_cyp3a4`."""
    return f"f_{protein.lower()}"


_PARTS = "|".join(p.value for p in Part)
_POOLS = "|".join(p.value for p in Pool)


def species_of(s: str, ids: Iterable[str]) -> list[str]:
    """The flat species ids of substance `s`, by the id scheme.

    Only species of exactly this substance match: `api` does not match the
    species of `api_m2`, `drg` not those of `xdrg`.

    Args:
        s: the substance id.
        ids: the ids to select from.

    Returns:
        The ids which are an organ species, a pool species or a blood cell
        pool species of `s`, in the order of `ids`.
    """
    pattern = re.compile(
        rf"C(?:(?:{_PARTS})|(?:{_POOLS})(?:_{Part.BLOOD_CELLS.value})?)_{re.escape(s)}"
    )
    return [i for i in ids if pattern.fullmatch(i.rsplit("__", 1)[-1])]
