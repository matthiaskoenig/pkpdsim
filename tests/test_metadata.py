"""The shared vocabulary of the metadata of the fit mappings."""

import dataclasses
import math
import re
from dataclasses import dataclass
from enum import StrEnum

import pytest
from sbmlsim.fit.objects import MappingMetaData

from pkpdsim import metadata
from pkpdsim.metadata import (
    AdministrationSite,
    ApplicationForm,
    Coadministration,
    DataType,
    Dosing,
    Fasting,
    Genotype,
    Health,
    Phenotype,
    PKPDData,
    PKPDMappingMetaData,
    Route,
    SamplingSite,
    Tissue,
)

VOCABULARIES = [
    obj
    for obj in vars(metadata).values()
    if isinstance(obj, type) and issubclass(obj, StrEnum) and obj is not StrEnum
]


@dataclass(kw_only=True)
class ModelMappingMetaData(PKPDMappingMetaData):
    diplotype: str | None = None
    activity_score: float = math.nan


def minimal() -> ModelMappingMetaData:
    return ModelMappingMetaData(
        tissue=Tissue.PLASMA,
        route=Route.PO,
        dosing=Dosing.SINGLE,
        health=Health.HEALTHY,
    )


def test_the_metadata_is_the_mapping_metadata_of_sbmlsim() -> None:
    assert issubclass(PKPDMappingMetaData, MappingMetaData)


def test_application_form_and_fasting_default_to_not_reported() -> None:
    metadata = minimal()
    assert metadata.application_form is ApplicationForm.NR
    assert metadata.fasting is Fasting.NR


def test_tissue_route_and_dosing_are_required() -> None:
    with pytest.raises(TypeError, match="dosing"):
        ModelMappingMetaData(  # ty: ignore[missing-argument]
            tissue=Tissue.PLASMA, route=Route.PO, health=Health.HEALTHY
        )


def test_to_dict_names_every_field_in_the_order_of_the_fields() -> None:
    metadata = ModelMappingMetaData(
        tissue=Tissue.URINE,
        route=Route.IV,
        application_form=ApplicationForm.SOLUTION,
        dosing=Dosing.MULTIPLE,
        fasting=Fasting.FED,
        health=Health.HEALTHY,
        coadministration=Coadministration.RIFAMPICIN,
        genotype=Genotype.NR,
    )
    assert list(metadata.to_dict().items()) == [
        ("tissue", "URINE"),
        ("route", "IV"),
        ("application_form", "SOLUTION"),
        ("dosing", "MULTIPLE"),
        ("fasting", "FED"),
        ("health", "HEALTHY"),
        ("coadministration", "RIFAMPICIN"),
        ("genotype", "NR"),
        ("phenotype", "NR"),
        ("administration_site", "NR"),
        ("data", "NR"),
        ("data_type", "NR"),
        ("sampling_site", "NOT_APPLICABLE"),
        ("diplotype", "None"),
        ("activity_score", "nan"),
    ]


def test_other_values_are_their_strings() -> None:
    metadata = ModelMappingMetaData(
        tissue=Tissue.PLASMA,
        route=Route.PO,
        dosing=Dosing.SINGLE,
        health=Health.HEALTHY,
        diplotype="*1/*4",
        activity_score=1.0,
    )
    assert metadata.to_dict()["diplotype"] == "*1/*4"
    assert metadata.to_dict()["activity_score"] == "1.0"


@pytest.mark.parametrize(
    ("enum", "members"),
    [
        (
            Tissue,
            "PLASMA SERUM BLOOD ARTERIAL_BLOOD ARTERIAL_PLASMA SALIVA URINE FECES "
            "FECES_URINE BILE NONE",
        ),
        (
            Route,
            "PO IV SC IM SL BU RE INHALATION INTESTINAL PATCH RING MIXED NR",
        ),
        (Dosing, "SINGLE MULTIPLE CONSTANT_INFUSION BOTH NR"),
        (Fasting, "NR FASTED FED LIGHT_BREAKFAST"),
        (
            ApplicationForm,
            "TABLET EFFERVESCENT_TABLET EXTENDED_RELEASE_TABLET "
            "EXTENDED_RELEASE_GRANULES SUSTAINED_RELEASE_TABLET "
            "PROLONGED_RELEASE_TABLET DELAYED_RELEASE_TABLET IMMEDIATE_RELEASE "
            "SLOW_RELEASE CONTROLLED_RELEASE CAPSULE SOLUTION BUFFERED_SOLUTION "
            "SUSPENSION SYRUP POWDER GRANULES FILM NANOPARTICLES SUPPOSITORY PATCH "
            "RING SUBCUTANEOUS PROPACETAMOL MIXED NR",
        ),
        (DataType, "MEAN MEDIAN INDIVIDUAL NR"),
        (
            SamplingSite,
            "ARM_VEIN VENOUS ARTERIAL ARTERIALISED HEPATIC_VEIN CAPILLARY "
            "NOT_APPLICABLE",
        ),
        (PKPDData, "PK PD NR"),
    ],
)
def test_the_vocabularies_hold_the_members_of_the_models(
    enum: type[StrEnum], members: str
) -> None:
    """The union of the members the models of pkdb_models use, synonyms merged.

    `SamplingSite` holds the sites of the observers of the PBPK library
    (`pkpdsim.pbpk.observers`), which the models take up when they are built
    on it; the hand-written models do not set it.
    """
    assert list(enum.__members__) == members.split()


@pytest.mark.parametrize("enum", VOCABULARIES, ids=lambda e: e.__name__)
def test_the_values_are_unique_and_text(enum: type[StrEnum]) -> None:
    """A value is one member (no alias) and readable text, not a code."""
    values = [member.value for member in enum.__members__.values()]
    assert len(values) == len(set(values))
    assert [value for value in values if "_" in value] == []


def test_the_conventional_members() -> None:
    assert Health.HEALTHY == "healthy"
    assert Health.NR == "not reported"
    assert Coadministration.NONE == "none"
    for enum in (Genotype, Phenotype, AdministrationSite, PKPDData, DataType):
        assert enum.NR == "not reported"


def test_the_members_of_the_apixaban_studies() -> None:
    """The NAFLD group of Tirona2018 and the naproxen arm of Frost2014a."""
    assert Health("non-alcoholic fatty liver disease") is Health.NAFLD
    assert Coadministration("naproxen") is Coadministration.NAPROXEN


def test_the_arterial_plasma_is_the_term_of_pk_db() -> None:
    """The tissue `arterial plasma` of the vocabulary of PK-DB."""
    assert Tissue("arterial plasma") is Tissue.ARTERIAL_PLASMA


def test_every_genotype_names_its_gene() -> None:
    """The value starts with the gene, and so does the name."""
    unnamed = [
        member
        for member in Genotype
        if member is not Genotype.NR
        and not (
            re.match(r"[A-Z][A-Z0-9]+ ", member.value)
            and member.name.startswith(member.value.split()[0])
        )
    ]
    assert unnamed == []


def test_the_fields_and_their_defaults() -> None:
    md = PKPDMappingMetaData(
        tissue=Tissue.PLASMA,
        route=Route.PO,
        dosing=Dosing.SINGLE,
        health=Health.HEALTHY,
    )
    assert md.coadministration is Coadministration.NONE
    assert md.genotype is Genotype.NR
    assert md.phenotype is Phenotype.NR
    assert md.administration_site is AdministrationSite.NR
    assert md.data is PKPDData.NR
    assert md.data_type is DataType.NR
    assert [f.name for f in dataclasses.fields(md)] == [
        "tissue",
        "route",
        "application_form",
        "dosing",
        "fasting",
        "health",
        "coadministration",
        "genotype",
        "phenotype",
        "administration_site",
        "data",
        "data_type",
        "sampling_site",
    ]


def test_health_is_required() -> None:
    with pytest.raises(TypeError, match="health"):
        PKPDMappingMetaData(  # ty: ignore[missing-argument]
            tissue=Tissue.PLASMA, route=Route.PO, dosing=Dosing.SINGLE
        )


def test_sampling_site_default_arm_vein() -> None:
    """Plasma and blood data without a stated site are arm-vein samples."""
    md = minimal()
    assert md.sampling_site is None
    assert md.resolved_sampling_site() is SamplingSite.ARM_VEIN
    assert SamplingSite("arm vein") is SamplingSite.ARM_VEIN
    urine = dataclasses.replace(
        md, tissue=Tissue.URINE, sampling_site=SamplingSite.NOT_APPLICABLE
    )
    assert urine.to_dict()["sampling_site"] == "NOT_APPLICABLE"


@pytest.mark.parametrize(
    "tissue",
    [Tissue.URINE, Tissue.FECES, Tissue.FECES_URINE, Tissue.BILE, Tissue.NONE],
)
def test_urine_without_a_site_is_not_applicable(tissue: Tissue) -> None:
    md = ModelMappingMetaData(
        tissue=tissue, route=Route.PO, dosing=Dosing.SINGLE, health=Health.HEALTHY
    )
    assert md.sampling_site is None
    assert md.resolved_sampling_site() is SamplingSite.NOT_APPLICABLE
    assert md.to_dict()["sampling_site"] == "NOT_APPLICABLE"


#: every tissue by its class: a new member of `Tissue` fails
#: `test_every_tissue_is_classified` until it is classified here and in the
#: metadata
VENOUS_BLOOD = (Tissue.PLASMA, Tissue.SERUM, Tissue.BLOOD)
ARTERIAL_BLOOD = (Tissue.ARTERIAL_BLOOD, Tissue.ARTERIAL_PLASMA)
NO_BLOOD = (
    Tissue.SALIVA,
    Tissue.URINE,
    Tissue.FECES,
    Tissue.FECES_URINE,
    Tissue.BILE,
    Tissue.NONE,
)
DEFAULT_SITE = {
    **dict.fromkeys(VENOUS_BLOOD, SamplingSite.ARM_VEIN),
    **dict.fromkeys(ARTERIAL_BLOOD, SamplingSite.ARTERIAL),
    **dict.fromkeys(NO_BLOOD, SamplingSite.NOT_APPLICABLE),
}


def test_every_tissue_is_classified() -> None:
    """Every member of `Tissue` is in exactly one class of this test."""
    classes = VENOUS_BLOOD + ARTERIAL_BLOOD + NO_BLOOD
    assert len(classes) == len(set(classes))
    assert set(classes) == set(Tissue)


@pytest.mark.parametrize("tissue", list(Tissue))
def test_the_site_without_a_statement_follows_the_tissue(tissue: Tissue) -> None:
    md = ModelMappingMetaData(
        tissue=tissue, route=Route.PO, dosing=Dosing.SINGLE, health=Health.HEALTHY
    )
    assert md.sampling_site is None
    assert md.resolved_sampling_site() is DEFAULT_SITE[tissue]
    assert md.to_dict()["sampling_site"] == DEFAULT_SITE[tissue].name


@pytest.mark.parametrize("tissue", list(Tissue))
@pytest.mark.parametrize("site", list(SamplingSite))
def test_a_stated_site_is_kept_or_refused(tissue: Tissue, site: SamplingSite) -> None:
    """A tissue takes its own site; a venous blood tissue any blood site."""
    if tissue in VENOUS_BLOOD:
        allowed = set(SamplingSite) - {SamplingSite.NOT_APPLICABLE}
    else:
        allowed = {DEFAULT_SITE[tissue]}

    def create() -> ModelMappingMetaData:
        return ModelMappingMetaData(
            tissue=tissue,
            route=Route.PO,
            dosing=Dosing.SINGLE,
            health=Health.HEALTHY,
            sampling_site=site,
        )

    if site in allowed:
        assert create().resolved_sampling_site() is site
    else:
        with pytest.raises(ValueError, match=tissue.name):
            create()


@pytest.mark.parametrize(
    "site",
    sorted(set(SamplingSite) - {SamplingSite.ARTERIAL, SamplingSite.NOT_APPLICABLE}),
)
def test_a_non_arterial_site_on_arterial_blood_is_refused(site: SamplingSite) -> None:
    with pytest.raises(ValueError, match="arterial"):
        ModelMappingMetaData(
            tissue=Tissue.ARTERIAL_PLASMA,
            route=Route.PO,
            dosing=Dosing.SINGLE,
            health=Health.HEALTHY,
            sampling_site=site,
        )


def test_replace_of_an_arterial_tissue_follows_its_site() -> None:
    arterial = dataclasses.replace(minimal(), tissue=Tissue.ARTERIAL_BLOOD)
    assert arterial.resolved_sampling_site() is SamplingSite.ARTERIAL
    assert arterial.to_dict()["sampling_site"] == "ARTERIAL"
    back = dataclasses.replace(arterial, tissue=Tissue.PLASMA)
    assert back.resolved_sampling_site() is SamplingSite.ARM_VEIN


@pytest.mark.parametrize("tissue", [Tissue.PLASMA, Tissue.SERUM, Tissue.BLOOD])
def test_a_blood_tissue_keeps_an_explicit_site(tissue: Tissue) -> None:
    md = ModelMappingMetaData(
        tissue=tissue,
        route=Route.PO,
        dosing=Dosing.SINGLE,
        health=Health.HEALTHY,
        sampling_site=SamplingSite.ARTERIAL,
    )
    assert md.sampling_site is SamplingSite.ARTERIAL


def test_a_blood_site_on_urine_is_refused() -> None:
    for site in (SamplingSite.ARM_VEIN, SamplingSite.VENOUS):
        with pytest.raises(ValueError, match="no blood sample"):
            ModelMappingMetaData(
                tissue=Tissue.URINE,
                route=Route.PO,
                dosing=Dosing.SINGLE,
                health=Health.HEALTHY,
                sampling_site=site,
            )


def test_not_applicable_on_plasma_is_refused() -> None:
    with pytest.raises(ValueError, match="blood sample"):
        ModelMappingMetaData(
            tissue=Tissue.PLASMA,
            route=Route.PO,
            dosing=Dosing.SINGLE,
            health=Health.HEALTHY,
            sampling_site=SamplingSite.NOT_APPLICABLE,
        )


def test_replace_of_the_tissue_follows_the_default_site() -> None:
    plasma = minimal()
    urine = dataclasses.replace(plasma, tissue=Tissue.URINE)
    assert urine.resolved_sampling_site() is SamplingSite.NOT_APPLICABLE
    back = dataclasses.replace(urine, tissue=Tissue.PLASMA)
    assert back.resolved_sampling_site() is SamplingSite.ARM_VEIN


def test_an_explicit_site_is_kept_by_replace_of_another_field() -> None:
    md = dataclasses.replace(minimal(), sampling_site=SamplingSite.ARTERIAL)
    other = dataclasses.replace(md, route=Route.IV)
    assert other.sampling_site is SamplingSite.ARTERIAL
    assert other.resolved_sampling_site() is SamplingSite.ARTERIAL
    with pytest.raises(ValueError, match="no blood sample"):
        dataclasses.replace(md, tissue=Tissue.URINE)


@dataclasses.dataclass(kw_only=True)
class OwnPostInit(PKPDMappingMetaData):
    extra: int = 0

    def __post_init__(self) -> None:  # no super() call
        self.extra += 1


def test_a_subclass_with_its_own_post_init_resolves() -> None:
    urine = OwnPostInit(
        tissue=Tissue.URINE, route=Route.PO, dosing=Dosing.SINGLE, health=Health.HEALTHY
    )
    assert urine.resolved_sampling_site() is SamplingSite.NOT_APPLICABLE
    assert urine.to_dict()["sampling_site"] == "NOT_APPLICABLE"
    plasma = OwnPostInit(
        tissue=Tissue.PLASMA,
        route=Route.PO,
        dosing=Dosing.SINGLE,
        health=Health.HEALTHY,
    )
    assert plasma.to_dict()["sampling_site"] == "ARM_VEIN"
