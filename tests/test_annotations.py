"""The annotations of the compartments are terms of identifiers.org.

A resource is written as `collection/term`, the collection is a namespace of
the registry of identifiers.org and the term matches the pattern of the
namespace, which is what pymetadata checks when sbmlutils writes the
annotation. The registry is in the cache of pymetadata, see `conftest.py`.
"""

import re

import pytest
from pymetadata.core.annotation import ProviderType, RDFAnnotation
from sbmlutils.metadata import BQB

from pkpdsim.annotations import compartments

RESOURCES = [
    pytest.param(qualifier, resource, id=f"{key}-{resource}")
    for key, annotations in compartments.items()
    for qualifier, resource in annotations
]

STRICT_PATTERNS = {
    # identifiers.org accepts any number of digits, the ids of UBERON have seven
    "uberon": re.compile(r"UBERON:\d{7}"),
}
"""Patterns of collections whose pattern in the registry is too permissive."""


def test_compartments() -> None:
    """Every annotation is a list of distinct qualifier and resource pairs."""
    assert {"li", "ki", "gu", "ve", "ar", "po", "plasma", "urine"} <= set(compartments)
    for key, annotations in compartments.items():
        for qualifier, resource in annotations:
            assert isinstance(qualifier, BQB), key
            assert isinstance(resource, str) and resource, key
        assert len(set(annotations)) == len(annotations), key


def check_resource(qualifier: BQB, resource: str) -> None:
    """Check that a resource is a term of a collection of identifiers.org."""
    annotation = RDFAnnotation(qualifier=qualifier, resource=resource, validate=False)
    assert annotation.provider == ProviderType.IDENTIFIERS_ORG, resource
    # the form of the module, the collection in lower case and the term as
    # identifiers.org writes it, e.g. `fma/FMA:7197` and `ncit/C12392`
    assert f"{annotation.collection}/{annotation.term}" == resource
    # the collection exists and the term matches its pattern
    assert annotation.check_miriam_term(), resource
    if pattern := STRICT_PATTERNS.get(annotation.collection or ""):
        assert pattern.fullmatch(annotation.term or ""), resource


@pytest.mark.parametrize(("qualifier", "resource"), RESOURCES)
def test_resource(qualifier: BQB, resource: str) -> None:
    """The resource is a term of a collection of identifiers.org."""
    check_resource(qualifier, resource)


@pytest.mark.parametrize(
    "resource",
    [
        "uberon/UBERON:000115",  # six digits
        "GO:0005886",  # without the collection
        "http://snomed.info/id/57041003",  # an url
        "ncit/NCIT:C12392",  # the prefix of a collection which is not in the term
        "bto/BTO:511",  # does not match the pattern
        "unknown/X:1",  # not a collection
    ],
)
def test_invalid_resource(resource: str) -> None:
    """The check of the resources rejects the errors it was written for."""
    with pytest.raises(AssertionError):
        check_resource(BQB.IS, resource)
