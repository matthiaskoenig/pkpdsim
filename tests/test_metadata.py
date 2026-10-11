"""The shared vocabulary of the metadata of the fit mappings."""

import math
from dataclasses import dataclass
from enum import StrEnum

import pytest
from sbmlsim.fit.objects import MappingMetaData

from pkpdsim.metadata import (
    ApplicationForm,
    DataType,
    Dosing,
    Fasting,
    PKPDData,
    PKPDMappingMetaData,
    Route,
    Tissue,
)


class Health(StrEnum):
    HEALTHY = "healthy"


class Coadministration(StrEnum):
    NONE = "none"
    RIFAMPICIN = "rifampicin"


@dataclass(kw_only=True)
class ModelMappingMetaData(PKPDMappingMetaData):
    health: Health
    coadministration: Coadministration = Coadministration.NONE
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
    )
    assert list(metadata.to_dict().items()) == [
        ("tissue", "URINE"),
        ("route", "IV"),
        ("application_form", "SOLUTION"),
        ("dosing", "MULTIPLE"),
        ("fasting", "FED"),
        ("health", "HEALTHY"),
        ("coadministration", "RIFAMPICIN"),
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
            "PLASMA SERUM BLOOD ARTERIAL_BLOOD SALIVA URINE FECES FECES_URINE BILE NONE",
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
        (DataType, "MEAN MEDIAN INDIVIDUAL"),
        (PKPDData, "PK PD"),
    ],
)
def test_the_vocabularies_hold_the_members_of_the_models(
    enum: type[StrEnum], members: str
) -> None:
    """The union of the members the models of pkdb_models use, synonyms merged."""
    assert list(enum.__members__) == members.split()


@pytest.mark.parametrize(
    "enum", [Tissue, Route, Dosing, Fasting, ApplicationForm, DataType, PKPDData]
)
def test_the_values_are_unique_and_text(enum: type[StrEnum]) -> None:
    values = [member.value for member in enum.__members__.values()]
    assert len(values) == len(set(values))
    assert all("_" not in value for value in values)
