# Chemistry

The chemistry of a model is a list of processes, each defined once per organ: a metabolism, an active transport or a custom reaction. The organ (the standard organ or its variant) places the process into its parts: the standard organ into the part which holds the substrate, the zonated liver into every zone, the intestine into the mucosa of every segment. A process has no compartments; it names its organ, its substances, its kinetics and its protein.

```python
from pkpdsim.pbpk import (
    MM,
    Binding,
    Direction,
    FirstOrder,
    Hill,
    Induction,
    Inhibition,
    Kp,
    Membrane,
    Metabolism,
    Organ,
    Substance,
    Transport,
)

chemistry = [
    # liver: drg -> drg_m by CYP3A4, Michaelis-Menten on the unbound drug
    Metabolism(
        "DRG2M",
        Organ.LIVER,
        "drg",
        products={"drg_m": 1.0},
        kinetics=MM(vmax=1e-3, km=0.01),
        protein="CYP3A4",
    ),
    # liver: elimination of prf by CYP2C9, first order per litre of liver
    Metabolism(
        "PRFEL",
        Organ.LIVER,
        "prf",
        {},
        FirstOrder(k=0.2),
        protein="CYP2C9",
    ),
    # kidney: elimination of hil with Hill kinetics, no protein
    Metabolism("HILEL", Organ.KIDNEY, "hil", {}, Hill(vmax=1e-3, km=0.01, n=2.0)),
    # liver: elimination of inh by the CYP2D6 it inactivates and induces
    Metabolism(
        "INHEL", Organ.LIVER, "inh", {}, MM(vmax=1e-3, km=0.02), protein="CYP2D6"
    ),
    # liver: uptake of prm by OATP1B1 and its excretion into the bile by MRP2
    Transport(
        "PRMIN",
        Organ.LIVER,
        "prm",
        Direction.INFLUX,
        Membrane.BASOLATERAL,
        MM(vmax=1e-2, km=0.05),
        protein="SLCO1B1",
    ),
    Transport(
        "PRMBIL",
        Organ.LIVER,
        "prm",
        Direction.EFFLUX,
        Membrane.APICAL,
        FirstOrder(k=0.1),
        protein="ABCC2",
    ),
]
inhibitor = Substance(
    "inh",
    name="inhibitor",
    mr=350.0,
    binding=Binding(fu=0.3),
    distribution=Kp(1.0, ftissue=1.0),
    inhibits=(Inhibition("CYP2D6", "mechanism", ki=0.005, kinact=0.05),),
    induces=(Induction("CYP2D6", emax=2.0, ec50=0.1),),
)
print([p.sid for p in chemistry], [i.kind for i in inhibitor.inhibits])
```

```text
['DRG2M', 'PRFEL', 'HILEL', 'INHEL', 'PRMIN', 'PRMBIL'] ['mechanism']
```

The equations on this page are those which these processes generate in the toy model of `tests/docs/test_pbpk_equations.py` (with the substances `drg`, `drg_m`, `prf`, `prm`, `hil`, `inh` of [Standard organ](organ.md)), as local ids of the liver (`LI__DRG2M` in the flat model) and the kidney.

## Processes

| process | acts in | placed |
|---|---|---|
| `Metabolism(sid, organ, substrate, products, kinetics, protein, profile)` | the part which holds the tissue of the substrate: intracellular in the full form, the merged tissue, the whole organ when flow-limited, the plasma without tissue | the products by their stoichiometric fractions (sum 1) into their equivalent part, else into the plasma; no products is an elimination |
| `Transport(sid, organ, substance, direction, membrane, kinetics, protein)` | basolateral: between interstitium and cell; apical: between cell and the sink of the organ (bile, urine, intestinal lumen) | influx into the cell, efflux out of it; forces the full form of the substance in the organ |
| `CustomReaction(sid, organ, part, equation, formula, parameters, species)` | the part it names | its formula as given, the escape hatch for a mechanism without a class (saturable binding, a deep tissue pool) |

The organs of the chemistry are those of `Organ`; the chemistry of the intestine is defined on `Organ.SMALL_INTESTINE` and the intestine places it into the mucosa of its segments (see [Intestine and routes](intestine.md)). The standard kidney refuses an apical influx from the urine (the urine is cumulative); a tubular reabsorption by a transporter needs the [nephron kidney](variants.md#nephron-kidney).

## Kinetics

Every rate acts on the unbound concentration of its part and is an amount per time. The Vmax of a process is per litre of organ, so the whole organ has the capacity `Vmax*V<o>` and follows the organ of the individual. The clearance `k` of a first order process is per litre of organ as well (1/min, `FirstOrder(k)`), so the whole organ has `k*V<o>` and follows the organ of the individual as a Vmax does, as the first order constants of the hand-written models (`k*V<o>`) did; `FirstOrder(k, per_volume=False)` states the clearance of the whole organ (l/min), which does not scale with the organ volume. The activity `f_<protein>` of the protein multiplies the rate, and in the liver the impairment factor `fimp_<protein>` and the functional fraction `f_functional` multiply it as well. In the kidney every metabolism and transport runs with the tubular capacity `1 - f_tubular + f_tubular*f_renal_function`, which follows the renal function: the share `f_tubular` is 1, the intact nephron hypothesis (Bricker 1960, a choice of the library), so a tubular secretion falls with the GFR as the secretion of a fraction of the GFR of the hand-written models did; the filtration `fu*GFR` follows the GFR itself (see [Impairment](impairment.md#chronic-kidney-disease)). A custom reaction states its rate itself and is not scaled.

```equations title="Michaelis-Menten, liver"
DRG2M: Cti_drg -> Cpl_drg_m = f_cyp3a4 * fimp_cyp3a4 * f_functional * DRG2M_Vmax * Vli * (fu_drg / Kp_drg * Cti_drg) / (DRG2M_Km + fu_drg / Kp_drg * Cti_drg)
```

```equations title="Hill, kidney (no protein, the tubular capacity)"
HILEL: Cti_hil -> = (1 - f_tubular + f_tubular * f_renal_function) * HILEL_Vmax * Vki * (fu_hil / Kp_hil * Cti_hil / HILEL_Km)^HILEL_n / (1 + (fu_hil / Kp_hil * Cti_hil / HILEL_Km)^HILEL_n)
```

```equations title="transport, liver (full form of prm)"
PRMIN: Cis_prm -> Cic_prm = f_slco1b1 * fimp_slco1b1 * f_functional * PRMIN_Vmax * Vli * (fu_prm / K_int_pls_prm * Cis_prm) / (PRMIN_Km + fu_prm / K_int_pls_prm * Cis_prm)
PRMBIL: Cic_prm -> Cbile_prm = f_abcc2 * fimp_abcc2 * f_functional * PRMBIL_k * Vli * fu_prm / K_cell_pls_prm * Cic_prm
```

The parameters of a process are constants of its organ with the process id as prefix (`DRG2M_Vmax`, `DRG2M_Km`, `PRFEL_k`, `HILEL_n`), so a fit of a Vmax uses the flat id `LI__DRG2M_Vmax`, as the hand-written models do. The tubular share is the constant `f_tubular` of the kidney (`KI__f_tubular`), with its provenance in the notes; an experiment may set it to test another hypothesis.

## Proteins

A protein is named once in the protein table of the library with its UniProt and NCIT annotation; a process names it by its key. In every organ where a process uses it, the protein has the activity factor `f_<protein>` (1, a constant of the organ, `LI__f_cyp3a4` in the flat model): genotype, interactions of an experiment and other factors act on it, as `f_cyp3a4` and `f_cyp2c9` of the hand-written models. In the liver the impairment multiplies it by `fimp_<protein>` (the input `fimp_li_<protein>` of the body, see [Impairment](impairment.md)).

<!-- generated: proteins (scripts/pbpk_tables.py), do not edit -->
*Generated from `data/proteins.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the data file.*

| protein | name | UniProt | NCIT |
| --- | --- | --- | --- |
| `CYP1A2` | Cytochrome P450 1A2 | [P05177](https://www.uniprot.org/uniprotkb/P05177) | `C17673` |
| `CYP2C8` | Cytochrome P450 2C8 | [P10632](https://www.uniprot.org/uniprotkb/P10632) | `C46001` |
| `CYP2C9` | Cytochrome P450 2C9 | [P11712](https://www.uniprot.org/uniprotkb/P11712) | `C28750` |
| `CYP2C19` | Cytochrome P450 2C19 | [P33261](https://www.uniprot.org/uniprotkb/P33261) | `C25948` |
| `CYP2D6` | Cytochrome P450 2D6 | [P10635](https://www.uniprot.org/uniprotkb/P10635) | `C17287` |
| `CYP3A4` | Cytochrome P450 3A4 | [P08684](https://www.uniprot.org/uniprotkb/P08684) | `C17573` |
| `UGT1A1` | UDP-glucuronosyltransferase 1A1 | [P22309](https://www.uniprot.org/uniprotkb/P22309) | `C28288` |
| `UGT1A9` | UDP-glucuronosyltransferase 1A9 | [O60656](https://www.uniprot.org/uniprotkb/O60656) | `C105864` |
| `UGT2B7` | UDP-glucuronosyltransferase 2B7 | [P16662](https://www.uniprot.org/uniprotkb/P16662) | `C105577` |
| `SULT1A1` | Sulfotransferase 1A1 | [P50225](https://www.uniprot.org/uniprotkb/P50225) | `C107443` |
| `SLCO1B1` | Solute carrier organic anion transporter family member 1B1 | [Q9Y6L6](https://www.uniprot.org/uniprotkb/Q9Y6L6) | `C106616` |
| `SLCO1B3` | Solute carrier organic anion transporter family member 1B3 | [Q9NPD5](https://www.uniprot.org/uniprotkb/Q9NPD5) | `C106619` |
| `SLC22A1` | Solute carrier family 22 member 1 (OCT1) | [O15245](https://www.uniprot.org/uniprotkb/O15245) | `C105075` |
| `SLC22A2` | Solute carrier family 22 member 2 (OCT2) | [O15244](https://www.uniprot.org/uniprotkb/O15244) | `C107576` |
| `SLC22A6` | Solute carrier family 22 member 6 (OAT1) | [Q4U2R8](https://www.uniprot.org/uniprotkb/Q4U2R8) |  |
| `SLC22A8` | Solute carrier family 22 member 8 (OAT3) | [Q8TCC7](https://www.uniprot.org/uniprotkb/Q8TCC7) |  |
| `SLC47A1` | Multidrug and toxin extrusion protein 1 (MATE1) | [Q96FL8](https://www.uniprot.org/uniprotkb/Q96FL8) | `C113566` |
| `ABCB1` | ATP-dependent translocase ABCB1 (P-glycoprotein) | [P08183](https://www.uniprot.org/uniprotkb/P08183) | `C17746` |
| `ABCG2` | Broad substrate specificity ATP-binding cassette transporter ABCG2 | [Q9UNQ0](https://www.uniprot.org/uniprotkb/Q9UNQ0) | `C21253` |
| `ABCC2` | ATP-binding cassette sub-family C member 2 (MRP2) | [Q92887](https://www.uniprot.org/uniprotkb/Q92887) | `C103884` |

Degradation rate `kdeg_<protein>` of the turnover state (1/hr): li 0.0193, si 0.0301, other 0.0193.

??? note "Sources"

    - `Greenblatt2003`: Greenblatt DJ, von Moltke LL, Harmatz JS, Chen G, Weemhoff JL, Jen C, Kelley CJ, LeDuc BW, Zinny MA. Time course of recovery of cytochrome P450 3A function after single doses of grapefruit juice. Clin Pharmacol Ther 74(2):121-129, 2003. [doi:10.1016/S0009-9236(03)00118-8](https://doi.org/10.1016/S0009-9236(03)00118-8)
    - `ModelChoice`: Choice of the pkpdsim PBPK library, not a value of the literature; the notes of the value give the reasoning.
    - `RowlandYeo2011`: Rowland Yeo K, Walsky RL, Jamei M, Rostami-Hodjegan A, Tucker GT. Prediction of time-dependent CYP3A4 drug-drug interactions by physiologically based pharmacokinetic modelling: impact of inactivation parameters and enzyme turnover. Eur J Pharm Sci 43(3):160-173, 2011. [doi:10.1016/j.ejps.2011.04.008](https://doi.org/10.1016/j.ejps.2011.04.008)
    - `UniProt`: The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2023. Nucleic Acids Res 51(D1):D523-D531, 2023. [doi:10.1093/nar/gkac1052](https://doi.org/10.1093/nar/gkac1052)
<!-- end generated -->

## Inhibition and induction

An inhibition or an induction is declared on the inhibitor (`Substance.inhibits`, `Substance.induces`), and acts on every process of the protein in every organ, with the unbound concentration of the inhibitor in the compartment from which the process takes its substrate (the part of the substrate of a metabolism, the interstitium of a basolateral influx, the cell of an efflux): `Inhibition(protein, kind, ki, kiu, kinact)` with the kinds `competitive`, `noncompetitive`, `uncompetitive`, `mixed` and `mechanism` (mechanism-based), and `Induction(protein, emax, ec50)`. The constants are parameters of the top model, `Ki_<protein>_<s>`, `Kiu_<protein>_<s>`, `kinact_<protein>_<s>`, `Emax_<protein>_<s>`, `EC50_<protein>_<s>`. A substance inhibits a protein once and induces it once.

The rate laws of mutually exclusive inhibitors (Segel 1975), with the competitive sum `c = sum I/Ki` (competitive, mixed with its competitive constant, noncompetitive, and mechanism-based with `Ki = KI`) and the uncompetitive sum `u = sum I/Kiu` (uncompetitive, mixed, noncompetitive), are `Vmax*V*C/(Km*(1 + c) + C*(1 + u))` for Michaelis-Menten, the Hill law with `Vmax_app = Vmax/(1 + u)` and `Km_app = Km*(1 + c)/(1 + u)`, and `k*C/(1 + c)` for first order (`k*V*C/(1 + c)` per volume; the limit `C << Km`, which an uncompetitive inhibitor does not change). The competitive inhibition of CYP2C9 by `drg`:

```equations title="competitive inhibition of a first order process (per volume)"
PRFEL: Cor_prf -> = f_cyp2c9 * fimp_cyp2c9 * f_functional * PRFEL_k * Vli * fu_prf / K_or_prf * Cor_prf / (1 + fu_drg / Kp_drg * Cti_drg / Ki_cyp2c9_drg)
```

A mechanism-based inhibition or an induction gives the organ the turnover state `E_<protein>`, the amount of the protein relative to its baseline (initially 1), which multiplies the rates of the protein (Mayhew 2000, Almond 2009, Fahmi 2008); the degradation rate `kdeg_<protein>` is that of hepatic CYP3A4 (Rowland Yeo 2011) in every organ but the intestine (Greenblatt 2003), a choice of the library for the proteins without their own estimate. `inh` inactivates and induces the CYP2D6 which metabolises it:

```equations title="turnover of CYP2D6 under mechanism-based inhibition and induction"
d(E_cyp2d6)/dt = kdeg_cyp2d6 * (1 + Emax_cyp2d6_inh * (fu_inh / Kp_inh) * Cti_inh / (EC50_cyp2d6_inh + fu_inh / Kp_inh * Cti_inh)) - kdeg_cyp2d6 * E_cyp2d6 - kinact_cyp2d6_inh * (fu_inh / Kp_inh) * Cti_inh / (Ki_cyp2d6_inh + fu_inh / Kp_inh * Cti_inh) * E_cyp2d6
INHEL: Cti_inh -> = f_cyp2d6 * fimp_cyp2d6 * E_cyp2d6 * f_functional * INHEL_Vmax * Vli * (fu_inh / Kp_inh * Cti_inh) / (INHEL_Km * (1 + fu_inh / Kp_inh * Cti_inh / Ki_cyp2d6_inh) + fu_inh / Kp_inh * Cti_inh)
```

Models whose interactions are factors in the experiments (`LI__f_cyp3a4` set by a coadministration) keep them as factors; the declared inhibition is what makes interactions work across packages in a coupled physiology (`PBPKDefinition.combine`).

## How a variant places the chemistry

The chemistry is defined once per organ; the organ variant places it into its parts, scaled by volume, so that the total capacity is that of the standard organ:

- the zonated liver places every metabolism and transport into every zone with the share `w` of the zone, by volume (`f_zone = 1/n`) or by the abundance profile of the protein (`profiles`, `w = p/sum(p)` for zones of equal volume), see [Organ variants](variants.md#zonated-liver);
- the intestine places the processes of `Organ.SMALL_INTESTINE` into the mucosa of its segments by the relative expression of the protein times the mucosa volume (CYP3A4 of the data; any other protein by mucosa volume), see [Intestine and routes](intestine.md#mucosa);
- the nephron kidney places an apical transport between the cells and the lumen of the proximal tubule.

A custom reaction states its rate for the whole organ: the variants refuse it with an error which names the process, since an arbitrary formula cannot be split by volume.
