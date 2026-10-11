"""The vocabulary with which the fit mappings of the models are described.

Every fit mapping of an experiment carries metadata which says what its curve
is: the tissue, the route, the application form, the dosing, the fasting, the
health of the subjects and what was coadministered. The fits select their data
by these values (see `pkpdsim.fit.metadata_filter`), so a model describes its
studies with the same words as the other models:

- the closed vocabularies are shared, the enums of this module;
- the open vocabularies (the health of the subjects, the coadministered
  substances, the genotypes) are the model's own `StrEnum`s, which grow with
  its studies. Their conventional members are `HEALTHY`, `NONE` (nothing
  coadministered) and `NR` (not reported), and a member has the same name and
  value in every model.

A model derives the metadata of its mappings from `PKPDMappingMetaData` and adds
its own fields:

```python
# not executed
@dataclass(kw_only=True)
class HctzMappingMetaData(PKPDMappingMetaData):
    health: Health
    coadministration: Coadministration = Coadministration.NONE
```

The tables of sbmlsim show the metadata by member name (`to_dict`); no code
reads the values, which are readable text.
"""

from __future__ import annotations

from dataclasses import dataclass, fields
from enum import Enum, StrEnum, unique
from typing import override

from sbmlsim.fit.objects import MappingMetaData

__all__ = [
    "ApplicationForm",
    "DataType",
    "Dosing",
    "Fasting",
    "PKPDData",
    "PKPDMappingMetaData",
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


@unique
class PKPDData(StrEnum):
    """Whether a curve is pharmacokinetic or pharmacodynamic data."""

    PK = "pk"
    PD = "pd"


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
    """

    tissue: Tissue
    route: Route
    application_form: ApplicationForm = ApplicationForm.NR
    dosing: Dosing
    fasting: Fasting = Fasting.NR

    @override
    def to_dict(self) -> dict[str, str]:
        """Get the metadata by field: the member names of the enums.

        Any other value is its string (`"None"` for `None`). The order is that of
        the fields, the fields of this class first.
        """
        return {field.name: _name(getattr(self, field.name)) for field in fields(self)}
