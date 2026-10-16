"""The vocabulary with which the fit mappings of the models are described.

Every fit mapping of an experiment carries metadata which says what its curve
is: the tissue, the route, the application form, the dosing, the fasting, the
health of the subjects, what was coadministered, the genotype and the phenotype
of the subjects, the site of the administration and the kind of the data. The
fits select their data by these values (see `pkpdsim.fit.metadata_filter`), so
all models describe their studies with the same words: there is one vocabulary
for all models, and every closed vocabulary is an enum of this module.
Their conventional members are `HEALTHY`, `NONE` (nothing coadministered) and
`NR` (not reported), and a member has one name and one value everywhere. A new
study which needs a new member adds it here.

A model uses `PKPDMappingMetaData` as it is. It extends it only for a facet
no other model has, in its `fitting/mapping_collections.py`, without
redefining a field or repeating a vocabulary:

```python
# not executed
@dataclass(kw_only=True)
class LocalMappingMetaData(PKPDMappingMetaData):
    diplotype: str | None = None
```

The tables of sbmlsim show the metadata by member name (`to_dict`). The values
are readable text, and they are API as much as the names: a model may look a
member up by its value, e.g. `Genotype(f"CYP2C9 {genotype}")` from the genotype
in the data of a study, so a value is not changed lightly either.
"""

from __future__ import annotations

from dataclasses import dataclass, fields
from enum import Enum, StrEnum, unique
from typing import override

from sbmlsim.fit.objects import MappingMetaData

__all__ = [
    "AdministrationSite",
    "ApplicationForm",
    "Coadministration",
    "DataType",
    "Dosing",
    "Fasting",
    "Genotype",
    "Health",
    "PKPDData",
    "PKPDMappingMetaData",
    "Phenotype",
    "Route",
    "Tissue",
]


@unique
class Tissue(StrEnum):
    """The tissue of a measurement."""

    PLASMA = "plasma"
    SERUM = "serum"
    BLOOD = "blood"
    ARTERIAL_BLOOD = "arterial blood"
    SALIVA = "saliva"
    URINE = "urine"
    FECES = "feces"
    FECES_URINE = "feces and urine"
    BILE = "bile"
    NONE = "none"  # not a sample of a tissue, e.g. the blood pressure


@unique
class Route(StrEnum):
    """The route of the administration of a substance.

    Not `pkpdutils.Route` (`ORAL`, `IV_BOLUS`), the route of a dose in the
    non-compartmental analysis; a module which needs both imports one under an
    alias.
    """

    PO = "po"
    IV = "iv"
    SC = "sc"
    IM = "im"
    SL = "sl"
    BU = "bu"  # buccal
    RE = "re"  # rectal
    INHALATION = "inhalation"
    INTESTINAL = "intestinal"
    PATCH = "patch"
    RING = "ring"
    MIXED = "mixed"  # several routes, e.g. po and iv
    NR = "not reported"


@unique
class Dosing(StrEnum):
    """The dosing of a substance."""

    SINGLE = "single"
    MULTIPLE = "multiple"
    CONSTANT_INFUSION = "infusion"
    BOTH = "single and multiple"
    NR = "not reported"


@unique
class Fasting(StrEnum):
    """The fasting state of the subjects at the administration."""

    NR = "not reported"
    FASTED = "fasted"
    FED = "fed"
    LIGHT_BREAKFAST = "light breakfast"


@unique
class ApplicationForm(StrEnum):
    """The pharmaceutical form of the administered substance."""

    TABLET = "tablet"
    EFFERVESCENT_TABLET = "effervescent tablet"
    EXTENDED_RELEASE_TABLET = "extended release tablet"
    EXTENDED_RELEASE_GRANULES = "extended release granules"
    SUSTAINED_RELEASE_TABLET = "sustained release tablet"
    PROLONGED_RELEASE_TABLET = "prolonged release tablet"
    DELAYED_RELEASE_TABLET = "delayed release tablet"
    IMMEDIATE_RELEASE = "immediate release"
    SLOW_RELEASE = "slow release"
    CONTROLLED_RELEASE = "controlled release"
    CAPSULE = "capsule"
    SOLUTION = "solution"
    BUFFERED_SOLUTION = "buffered solution"
    SUSPENSION = "suspension"
    SYRUP = "syrup"
    POWDER = "powder"
    GRANULES = "granules"
    FILM = "oral soluble film"
    NANOPARTICLES = "nanoparticles"
    SUPPOSITORY = "suppository"
    PATCH = "patch"
    RING = "ring"
    SUBCUTANEOUS = "subcutaneous"
    PROPACETAMOL = "propacetamol"  # the prodrug of acetaminophen, intravenous
    MIXED = "mixed"  # several forms, e.g. po and iv
    NR = "not reported"


@unique
class DataType(StrEnum):
    """What the values of a curve are."""

    MEAN = "mean"
    MEDIAN = "median"
    INDIVIDUAL = "individual"
    NR = "not reported"


@unique
class PKPDData(StrEnum):
    """Whether a curve is pharmacokinetic or pharmacodynamic data."""

    PK = "pk"
    PD = "pd"
    NR = "not reported"


@unique
class Health(StrEnum):
    """The health of the subjects."""

    HEALTHY = "healthy"
    NR = "not reported"
    ALCOHOLISM = "alcoholism"  # alcohol use disorder, drinking or abstinent
    BILIARY_OBSTRUCTION = "biliary obstruction"
    # adults with a cancer which has no member of its own (HCC, NSCLC, CML)
    CANCER = "cancer"
    CARDIAC_IMPAIRMENT = "cardiac impairment"
    CARDIAC_RENAL_IMPAIRMENT = "cardiac renal impairment"
    CHOLECYSTECTOMY = "cholecystectomy"  # after the removal of the gallbladder
    CIRRHOSIS = "cirrhosis"
    CIRRHOSIS_HEART_FAILURE = "cirrhosis & heart failure"
    CML = "chronic myelogenous leukemia"
    CYSTIC_FIBROSIS = "cystic fibrosis"
    DIABETES = "diabetes mellitus"
    DUODENAL_ULCER = "duodenal ulcer"
    EPILEPSY = "epilepsy"
    HCC = "hepatocellular carcinoma"
    HCC_RENAL_IMPAIRMENT = "hepatocellular carcinoma & renal impairment"
    # cardiac patients of several diseases, e.g. ischemic heart disease,
    # hypertension or mitral valve disease; not HEART_FAILURE
    HEART_DISEASE = "heart disease"
    HEART_FAILURE = "heart failure"
    HEPATIC_IMPAIRMENT = "hepatic impairment"
    HEPATITIS = "acute viral hepatitis"
    HEPATO_RENAL_IMPAIRMENT = "hepatorenal impairment"
    HOSPITALIZED = "hospitalized"
    HYPERTENSION = "hypertension"
    # regular drinking without alcohol use disorder
    MODERATE_DRINKING = "moderate alcohol consumption"
    NSCLC = "non-small cell lung cancer"
    OBESE = "obese"
    PAIN = "pain"
    PEDIATRIC_CANCER = "pediatric cancer"  # children and young adults
    PRIMARY_HYPERCHOLESTEROLEMIA = "primary hypercholesterolemia"
    RENAL_IMPAIRMENT = "renal impairment"
    RENAL_TRANSPLANT = "renal transplant"
    SURGERY = "surgical patients"
    T1DM = "type 1 diabetes mellitus"
    T2DM = "type 2 diabetes mellitus"
    T2DM_RENAL_IMPAIRMENT = "T2DM & renal impairment"


@unique
class Coadministration(StrEnum):
    """A substance administered together with the substance of the model."""

    NONE = "none"  # nothing coadministered
    ACETAMINOPHEN = "acetaminophen"
    ALFENTANIL = "alfentanil"
    ALISKIREN = "aliskiren"
    ALLOPURINOL = "allopurinol"
    AMENAMEVIR = "amenamevir"
    AMLODIPINE = "amlodipine"
    ANESTHESIA = "anesthesia"  # e.g. spinal or epidural (bupivacaine, fentanyl)
    ANTIEPILEPTICS = "antiepileptics"
    ANTIHYPERGLYCEMIC_MEDICATION = "antihyperglycemic medication"
    ANTIPYRINE = "antipyrine"
    APPLE = "apple"
    ASA = "aspirin"
    ATENOLOL = "atenolol"
    ATOGEPANT = "atogepant"
    ATORVASTATIN = "atorvastatin"
    ATROPINE = "atropine"
    BREVISCAPINE = "breviscapine"
    BUCOLOME = "bucolome"
    CABOTEGRAVIR = "cabotegravir"
    CANAGLIFLOZIN = "canagliflozin"
    CAPTOPRIL = "captopril"
    CARBAMAZEPINE = "carbamazepine"
    CARBOHYDRATE = "carbohydrate"
    CELECOXIB = "celecoxib"
    CETUXIMAB = "cetuximab"
    CHOLESTYRAMINE = "cholestyramine"
    CILAZAPRIL = "cilazapril"
    CILOSTAZOL = "cilostazol"
    CIMETIDINE = "cimetidine"
    CISAPRIDE = "cisapride"
    CISPLATIN = "cisplatin"
    CLARITHROMYCIN = "clarithromycin"
    CLOPIDOGREL = "clopidogrel"
    CLOTRIMAZOLE = "clotrimazole"
    COCA_COLA = "coca cola"  # the substance dissolved in Coca-Cola, i.e. caffeine
    # several other drugs, mostly the probe drugs of a phenotyping cocktail
    COCKTAIL = "phenotyping cocktail"
    COMEDICATION = "comedication"  # the chronic medication of patients
    CYCLOSPORINE = "cyclosporine"
    DAPAGLIFLOZIN = "dapagliflozin"
    DIGOXIN = "digoxin"
    DILTIAZEM = "diltiazem"
    EDOXABAN = "edoxaban"
    EMPAGLIFLOZIN = "empagliflozin"
    ENALAPRIL = "enalapril"
    ERYTHROMYCIN = "erythromycin"
    # the tracer of the erythromycin breath test, a probe of CYP3A4
    ERYTHROMYCIN_BREATH_TEST = "[14C N-methyl] erythromycin"
    EVOGLIPTIN = "evogliptin"
    FAMOTIDINE = "famotidine"
    FENOFIBRATE = "fenofibrate"
    FIMASARTAN = "fimasartan"
    FLIBANSERINE = "flibanserine"
    FUROSEMIDE = "furosemide"
    GABAPENTIN = "gabapentin"
    GEMFIBROZIL = "gemfibrozil"
    GEMFIBROZIL_ITRACONAZOLE = "gemfibrozil+itraconazole"
    GEMIGLIPTIN = "gemigliptin"
    GINSENG = "ginseng"
    GLIMEPIRIDE = "glimepiride"
    GRAPEFRUIT_JUICE = "grapefruit juice"
    GREEN_TEA = "green tea"
    GUAIFENESIN = "guaifenesin"
    HCTZ = "hydrochlorothiazide"
    IMATINIB = "imatinib"
    INDOMETACIN = "indometacin"
    INSULIN = "insulin"
    IRBESARTAN = "irbesartan"
    ISMN = "ismn"
    ITRACONAZOLE = "itraconazole"
    KETOCONAZOLE = "ketoconazole"
    LACIDIPINE = "lacidipine"
    LAPATINIB = "lapatinib"
    LCZ696 = "lcz696"
    LINAGLIPTIN = "linagliptin"
    LISINOPRIL = "lisinopril"
    LOBEGLITAZONE = "lobeglitazone"
    MARAVIROC = "maraviroc"
    METFORMIN = "metformin"
    METOCLOPRAMIDE = "metoclopramide"  # changes the gastric emptying
    MIDAZOLAM = "midazolam"
    MORPHINE = "morphine"
    MYCOPHENOLATE_MOFETIL = "mycophenolate mofetil"
    NELFINAVIR = "nelfinavir"
    ORAL_CONTRACEPTIVES = "oral contraceptives"
    ORANGE = "orange"
    ORLISTAT = "orlistat"
    OSIMERTINIB = "osimertinib"
    PACLITAXEL = "paclitaxel"
    PAH = "p-aminohippurate"  # a marker of the renal plasma flow
    PEMBROLIZUMAB = "pembrolizumab"
    PENTAGASTRIN = "pentagastrin"  # stimulates the gastric acid secretion
    PETHIDINE = "pethidine"
    PHENYTOIN = "phenytoin"
    PIOGLITAZONE = "pioglitazone"
    PRASUGREL = "prasugrel"
    PREDNISONE = "prednisone"
    PREGABALIN = "pregabalin"
    PROBENECID = "probenecid"
    PROTEIN = "protein"
    QUERCETIN = "quercetin"
    QUINIDINE = "quinidine"
    RAMIPRIL = "ramipril"
    RIFAMPICIN = "rifampicin"
    RITONAVIR = "ritonavir"
    ROSIGLITAZONE_METFORMIN = "rosiglitazone + metformin"
    ROSUVASTATIN = "rosuvastatin"
    SETIPIPRANT = "setipiprant"
    SILEXAN = "silexan"
    SILYMARIN = "silymarin"
    SIMVASTATIN = "simvastatin"
    SITAGLIPTIN = "sitagliptin"
    ST_JOHNS_WORT = "st john's wort extract"
    TACROLIMUS = "tacrolimus"
    TACROLIMUS_MYCOPHENOLATE_MOFETIL = "tacrolimus and mycophenolate mofetil"
    TORASEMIDE = "torasemide"
    TROGLITAZONE = "troglitazone"
    TROLEANDOMYCIN = "troleandomycin"
    UREA = "urea"
    VALSARTAN = "valsartan"
    VERAPAMIL = "verapamil"
    VINORELBINE = "vinorelbine"
    VOGLIBOSE = "voglibose"
    WARFARIN = "warfarin"


@unique
class Genotype(StrEnum):
    """The genotype of a group of subjects, the gene and its alleles."""

    NR = "not reported"
    ABCB1_2677AT_3435CT = "ABCB1 2677AT/3435CT"
    ABCB1_2677GG_3435CC = "ABCB1 2677GG/3435CC"
    ABCB1_2677GT_3435CT = "ABCB1 2677GT/3435CT"
    ABCB1_2677TT_3435TT = "ABCB1 2677TT/3435TT"
    ABCB1_3435CC = "ABCB1 3435CC"
    ABCB1_3435CT = "ABCB1 3435CT"
    ABCB1_3435TT = "ABCB1 3435TT"
    ABCB1_3435CT_OR_TT = "ABCB1 3435CT or 3435TT"
    # the haplotypes of c.1236C>T, c.2677G>T/A and c.3435C>T
    ABCB1_CGC_CGC = "ABCB1 CGC/CGC"
    ABCB1_TTT_TTT = "ABCB1 TTT/TTT"
    ABCG2_421AA = "ABCG2 c.421AA"
    ABCG2_421CA = "ABCG2 c.421CA"
    ABCG2_421CC = "ABCG2 c.421CC"
    CES1_GAIN_OF_FUNCTION = "CES1 mutation with increased activity"
    CES1_LOSS_OF_FUNCTION = "CES1 mutation with reduced activity"
    CES1_WILD_TYPE = "CES1 wild type"
    CYP2C9_1_1 = "CYP2C9 *1/*1"
    CYP2C9_1_13 = "CYP2C9 *1/*13"
    CYP2C9_1_1_AND_1_3 = "CYP2C9 *1/*1 and *1/*3"
    CYP2C9_1_2 = "CYP2C9 *1/*2"
    CYP2C9_1_3 = "CYP2C9 *1/*3"
    CYP2C9_2_2 = "CYP2C9 *2/*2"
    CYP2C9_2_3 = "CYP2C9 *2/*3"
    CYP2C9_3_3 = "CYP2C9 *3/*3"
    CYP3A4_1G_1G = "CYP3A4 *1G/*1G"
    CYP3A4_1_1 = "CYP3A4 *1/*1"
    CYP3A4_1_1G = "CYP3A4 *1/*1G"
    CYP3A5_1_1 = "CYP3A5 *1/*1"
    CYP3A5_1_3 = "CYP3A5 *1/*3"
    CYP3A5_3_3 = "CYP3A5 *3/*3"
    SLCO1B1_521CC = "SLCO1B1 c.521CC"
    SLCO1B1_521TC = "SLCO1B1 c.521TC"
    SLCO1B1_521TT = "SLCO1B1 c.521TT"
    # SLCO1B1 alone: the wild type is the reference genotype *1a/*1a (c.521TT),
    # a polymorphism any genotype with a variant allele (e.g. *1b, *5, *15, *17)
    SLCO1B1_POLYMORPHISM = "SLCO1B1 polymorphism"
    SLCO1B1_WILD_TYPE = "SLCO1B1 wild type"
    # SLCO1B1 and SLCO1B3 together: a variant is any genotype other than the
    # reference genotype of its gene (SLCO1B1 *1a/*1a)
    SLCO1B1_SLCO1B3_VARIANT = "SLCO1B1 and SLCO1B3 variant"
    SLCO1B1_SLCO1B3_WILD_TYPE = "SLCO1B1 and SLCO1B3 wild type"
    SLCO1B1_VARIANT_SLCO1B3_WILD_TYPE = "SLCO1B1 variant and SLCO1B3 wild type"
    SLCO1B1_WILD_TYPE_SLCO1B3_VARIANT = "SLCO1B1 wild type and SLCO1B3 variant"
    SLCO2B1_1457_CC = "SLCO2B1 c.1457 C/C"
    SLCO2B1_1457_CT = "SLCO2B1 c.1457 C/T"
    SLCO2B1_1457_TT = "SLCO2B1 c.1457 T/T"
    SLCO2B1_282_AA = "SLCO2B1 g.-282 A/A"
    SLCO2B1_282_GA = "SLCO2B1 g.-282 G/A"
    SLCO2B1_282_GG = "SLCO2B1 g.-282 G/G"
    SLCO2B1_601_AA = "SLCO2B1 c.601 A/A"
    SLCO2B1_601_GA = "SLCO2B1 c.601 G/A"
    SLCO2B1_601_GG = "SLCO2B1 c.601 G/G"
    SLCO2B1_935_AA = "SLCO2B1 c.935 A/A"
    SLCO2B1_935_GA = "SLCO2B1 c.935 G/A"
    SLCO2B1_935_GG = "SLCO2B1 c.935 G/G"
    UGT1A3_2_2 = "UGT1A3 *2/*2"
    UGT1A3_NONCARRIER = "UGT1A3 *1/*1 or *1/*3"
    UGT1A3_X_2 = "UGT1A3 *x/*2"


@unique
class Phenotype(StrEnum):
    """The metabolizer phenotype of the subjects for an enzyme."""

    NR = "not reported"
    CYP2C19_EM = "CYP2C19 extensive metabolizer"  # e.g. *1/*1
    CYP2C19_IM = "CYP2C19 intermediate metabolizer"  # e.g. *1/*2, *1/*3
    CYP2C19_PM = "CYP2C19 poor metabolizer"  # e.g. *2/*2, *2/*3, *3/*3
    CYP2D6_EM = "CYP2D6 extensive metabolizer"
    CYP2D6_IM = "CYP2D6 intermediate metabolizer"
    CYP2D6_PM = "CYP2D6 poor metabolizer"
    CYP2D6_UM = "CYP2D6 ultrarapid metabolizer"


@unique
class AdministrationSite(StrEnum):
    """The site of an injection or application."""

    NR = "not reported"
    ABDOMEN = "abdomen"
    THIGH = "thigh"
    UPPER_ARM = "upper arm"


def _name(value: object) -> str:
    """Get the member name of an enum value, the string of any other value."""
    return value.name if isinstance(value, Enum) else str(value)


@dataclass(kw_only=True)
class PKPDMappingMetaData(MappingMetaData):
    """The metadata of a fit mapping: the facets every model describes.

    The metadata describes the curve; what a fit does with it is decided by the
    filters of the fits. The fields are keyword-only, so a model adds fields
    without a default after the defaults of this class.

    Attributes:
        tissue: the tissue of the measurement.
        route: the route of the administration.
        application_form: the pharmaceutical form, not reported by default.
        dosing: the dosing.
        fasting: the fasting state, not reported by default.
        health: the health of the subjects.
        coadministration: the substance coadministered, none by default.
        genotype: the genotype of the subjects, not reported by default.
        phenotype: the metabolizer phenotype, not reported by default.
        administration_site: the site of an injection or application, not
            reported by default.
        data: whether the curve is pharmacokinetic or pharmacodynamic data, not
            reported by default.
        data_type: what the values of the curve are, not reported by default.
    """

    tissue: Tissue
    route: Route
    application_form: ApplicationForm = ApplicationForm.NR
    dosing: Dosing
    fasting: Fasting = Fasting.NR
    health: Health
    coadministration: Coadministration = Coadministration.NONE
    genotype: Genotype = Genotype.NR
    phenotype: Phenotype = Phenotype.NR
    administration_site: AdministrationSite = AdministrationSite.NR
    data: PKPDData = PKPDData.NR
    data_type: DataType = DataType.NR

    @override
    def to_dict(self) -> dict[str, str]:
        """Get the metadata by field: the member names of the enums.

        Any other value is its string (`"None"` for `None`). The order is that of
        the fields, the fields of this class first.
        """
        return {field.name: _name(getattr(self, field.name)) for field in fields(self)}
