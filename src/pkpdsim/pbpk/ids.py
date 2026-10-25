"""The id scheme of the PBPK models.

Every id of a generated model is built here, so that the organs, the
observers and the fits agree on them. A flat id is `<ORGAN>__<local id>`
(`LI__Cpl_api`), the ids of substances are suffixed by `_<substance>`.
"""

import re
from collections.abc import Iterable

from pkpdsim.pbpk.organs import Organ, Part, Pool

__all__ = [
    "activity",
    "arm",
    "arm_blood",
    "arm_free",
    "bile",
    "blood",
    "bp",
    "compartment",
    "feces",
    "flat",
    "free",
    "ftissue",
    "fu",
    "fu_ref",
    "hepatic_vein",
    "kbc",
    "kp",
    "kp_organ",
    "mr",
    "organ_flow",
    "organ_volume",
    "outflow_conc",
    "pool_compartment",
    "pool_species",
    "species",
    "species_of",
    "submodel",
    "urine",
]


def submodel(organ: Organ) -> str:
    """The id of the submodel of an organ, `LI`."""
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


def outflow_conc(s: str) -> str:
    """The local id of the concentration leaving an organ, `Cout_<s>`."""
    return f"Cout_{s}"


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


def urine(s: str) -> str:
    """The id of the amount in urine, `Aurine_<s>`."""
    return f"Aurine_{s}"


def feces(s: str) -> str:
    """The id of the amount in feces, `Afeces_<s>`."""
    return f"Afeces_{s}"


def bile(s: str) -> str:
    """The id of the concentration in bile, `Cbile_<s>`."""
    return f"Cbile_{s}"


def arm(s: str) -> str:
    """The id of the plasma concentration at the arm-vein site, `Carm_<s>`."""
    return f"Carm_{s}"


def arm_free(s: str) -> str:
    """The id of the unbound concentration at the arm-vein site, `Carm_free_<s>`."""
    return f"Carm_free_{s}"


def arm_blood(s: str) -> str:
    """The id of the blood concentration at the arm-vein site, `Carm_blood_<s>`."""
    return f"Carm_blood_{s}"


def free(pool: Pool, s: str) -> str:
    """The id of the unbound concentration in a pool, `Cve_free_<s>`."""
    return f"C{pool.value}_free_{s}"


def blood(pool: Pool, s: str) -> str:
    """The id of the blood concentration in a pool, `Cve_blood_<s>`."""
    return f"C{pool.value}_blood_{s}"


def hepatic_vein(s: str) -> str:
    """The id of the concentration in the hepatic vein, `Chv_<s>`."""
    return f"Chv_{s}"


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
