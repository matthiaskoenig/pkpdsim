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
        ("sampling_site", "ARM_VEIN"),
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
    """The union of the members the models of pkdb_models use, synonyms merged."""
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
    assert md.sampling_site is SamplingSite.ARM_VEIN
    assert SamplingSite("arm vein") is SamplingSite.ARM_VEIN
    urine = dataclasses.replace(md, sampling_site=SamplingSite.NOT_APPLICABLE)
    assert urine.to_dict()["sampling_site"] == "NOT_APPLICABLE"
