"""Tests of the annotation of every element the library creates (#63).

Two definitions together use every feature of the library: the organ
variants (zonated liver, nephron kidney, segmented intestine; the standard
organs with the single lumen), kinetic and equilibrium erythrocytes, a Kp
with an override, a Kp method, the flow-limited organ, active transport on
both membranes, competitive and mechanism-based inhibition, induction, the
impairment, the observers, a PD component and every route and formulation.
Every compartment, species, parameter and reaction of the comp model and of
the flat model has an SBO term and at least one ontology term.
"""

import logging
import re
from pathlib import Path

import libsbml
import pytest
from pymetadata.core.annotation import ProviderType, RDFAnnotation
from pymetadata.ontologies.pbpko import PBPKO
from sbmlutils.factory import Parameter, Species
from sbmlutils.metadata import BQB

from pkpdsim.annotations import compartments
from pkpdsim.pbpk import (
    MM,
    Binding,
    CustomReaction,
    Direction,
    DoseDependentDissolution,
    Erythrocytes,
    Excretion,
    FirstOrder,
    Hill,
    Induction,
    Inhibition,
    Kp,
    KpMethod,
    Membrane,
    Metabolism,
    Mode,
    Organ,
    Part,
    PBPKDefinition,
    Route,
    Substance,
    Transport,
    U,
    Variants,
    WeibullDissolution,
)
from pkpdsim.pbpk.annotation import (
    KINDS,
    ORGANS,
    POOLS,
    SEGMENTS,
    SITES,
    model_annotations,
    pbpko,
    requested_terms,
)
from pkpdsim.pbpk.builder import BuildResult, build

from .conftest import cardiac_pd

DOCS = Path(__file__).parents[2] / "docs" / "pbpk" / "annotation.md"
"""The documentation of the annotation, with the terms requested upstream."""

_CHEBI = [(BQB.IS, "chebi/CHEBI:15377")]
"""An annotation of a substance of its specification (water, a placeholder)."""


def _substances() -> list[Substance]:
    """A drug with kinetic erythrocytes, an inhibitor with a Kp method, a metabolite."""
    return [
        Substance(
            "drg",
            name="drug",
            mr=300.0,
            annotations=_CHEBI,
            binding=Binding(fu=0.1, protein="albumin"),
            erythrocytes=Erythrocytes("kinetic", kbc=2.0, k_exchange=1.0),
            distribution=Kp(2.0, ftissue=1.0, overrides={Organ.ADIPOSE: 0.5}),
            excretion=Excretion(urine=True, fraction_reabsorbed=0.1),
        ),
        Substance(
            "inh",
            name="inhibitor",
            mr=400.0,
            erythrocytes=Erythrocytes("equilibrium", bp=1.2),
            distribution=KpMethod(
                method="rodgers_rowland",
                logp=2.0,
                pka=(8.0,),
                compound="base",
                ftissue=1.0,
            ),
            inhibits=(
                Inhibition("CYP3A4", "mechanism", 1e-3, kinact=0.01),
                Inhibition("SLCO1B1", "competitive", 1e-3),
            ),
            induces=(Induction("CYP2C9", emax=1.0, ec50=1e-3),),
        ),
        Substance("drg_m", name="metabolite", mr=316.0),
    ]


def _variants() -> PBPKDefinition:
    """The organ variants with every process kind and every route."""
    return PBPKDefinition(
        sid="variants",
        substances=_substances(),
        chemistry=[
            Metabolism(
                "DRG2M", Organ.LIVER, "drg", {"drg_m": 1.0}, MM(1e-3, 0.01), "CYP3A4"
            ),
            Metabolism(
                "DRG2M9",
                Organ.LIVER,
                "drg",
                {"drg_m": 1.0},
                Hill(1e-3, 0.01, 2.0),
                "CYP2C9",
            ),
            Metabolism(
                "INHX", Organ.LIVER, "inh", {}, FirstOrder(0.01, per_volume=True)
            ),
            Transport(
                "DRGUP",
                Organ.LIVER,
                "drg",
                Direction.INFLUX,
                Membrane.BASOLATERAL,
                MM(1e-3, 0.01),
                "SLCO1B1",
            ),
            Transport(
                "DRGBI",
                Organ.LIVER,
                "drg",
                Direction.EFFLUX,
                Membrane.APICAL,
                FirstOrder(0.01),
                "ABCB1",
            ),
            Transport(
                "DRGSEC",
                Organ.KIDNEY,
                "drg",
                Direction.EFFLUX,
                Membrane.APICAL,
                FirstOrder(0.01),
            ),
            Metabolism(
                "DRGGUT",
                Organ.SMALL_INTESTINE,
                "drg",
                {"drg_m": 1.0},
                MM(1e-4, 0.01),
                "CYP3A4",
            ),
            Transport(
                "DRGPGP",
                Organ.SMALL_INTESTINE,
                "drg",
                Direction.EFFLUX,
                Membrane.APICAL,
                FirstOrder(0.01),
                "ABCB1",
            ),
        ],
        variants=Variants(
            organs={
                Organ.LIVER: "zonated",
                Organ.KIDNEY: "nephron",
                Organ.SMALL_INTESTINE: "segmented",
            }
        ),
        routes=[
            Route.iv("drg"),
            Route.oral("drg", WeibullDissolution(t50=30.0, shape=0.8, lag=5.0)),
            Route.oral_solution("drg"),
            Route.iv("inh"),
            Route.oral("inh", DoseDependentDissolution(k=1.0, dose50=10.0)),
        ],
        pd=[cardiac_pd()],
    )


def _standard() -> PBPKDefinition:
    """The standard organs: single lumen, flow-limited drug, custom reaction, dynamic factor."""
    flow_limited = Substance(
        "drg",
        name="drug",
        mr=300.0,
        annotations=_CHEBI,
        binding=Binding(fu=0.1, protein="albumin"),
        erythrocytes=Erythrocytes("equilibrium", kbc=2.0),
        distribution=Kp(2.0, mode=Mode.PERFUSION),
        excretion=Excretion(urine=True),
    )
    secreted = Substance(
        "sec",
        name="secreted drug",
        mr=250.0,
        distribution=Kp(
            1.0, organs=frozenset({Organ.KIDNEY}), ftissue=1.0, ps={Organ.KIDNEY: 0.2}
        ),
    )
    return PBPKDefinition(
        sid="standard",
        substances=[flow_limited, secreted, Substance("drg_m", "metabolite", 316.0)],
        chemistry=[
            Metabolism(
                "DRG2M", Organ.LIVER, "drg", {"drg_m": 1.0}, MM(1e-3, 0.01), "CYP3A4"
            ),
            Transport(
                "SECK",
                Organ.KIDNEY,
                "sec",
                Direction.EFFLUX,
                Membrane.APICAL,
                MM(1e-3, 0.01),
                "SLC22A8",
            ),
            Transport(
                "SECUP",
                Organ.KIDNEY,
                "sec",
                Direction.INFLUX,
                Membrane.BASOLATERAL,
                FirstOrder(0.1, per_volume=True),
            ),
            CustomReaction(
                sid="SECBIND",
                organ=Organ.KIDNEY,
                part=Part.INTRACELLULAR,
                equation="Cic_sec -> Cbound_sec",
                formula="kb_sec * Cic_sec",
                parameters=[Parameter("kb_sec", 0.01, U.l_per_min)],
                species=[
                    Species(
                        "Cbound_sec",
                        "Vic",
                        initialConcentration=0.0,
                        substanceUnit=U.mmole,
                        hasOnlySubstanceUnits=False,
                    )
                ],
            ),
        ],
        variants=Variants(organs={Organ.SMALL_INTESTINE: "single_lumen"}),
        routes=[Route.iv("drg"), Route.oral("drg"), Route.oral_solution("sec")],
        dynamic=("f_renal_function",),
    )


@pytest.fixture(scope="module", params=["variants", "standard"])
def built(
    request: pytest.FixtureRequest, tmp_path_factory: pytest.TempPathFactory
) -> BuildResult:
    """The build of a definition which uses a part of the features."""
    d = _variants() if request.param == "variants" else _standard()
    return build(d, tmp_path_factory.mktemp(request.param))


def _elements(model: libsbml.Model) -> list[libsbml.SBase]:
    """The compartments, species, parameters and reactions of a model."""
    return [
        *model.getListOfCompartments(),
        *model.getListOfSpecies(),
        *model.getListOfParameters(),
        *model.getListOfReactions(),
    ]


def _models(doc: libsbml.SBMLDocument) -> list[libsbml.Model]:
    """The model of a document and its comp model definitions."""
    models = [doc.getModel()]
    comp = doc.getPlugin("comp")
    if comp is not None:
        models += list(comp.getListOfModelDefinitions())
    return models


def test_every_generated_element_annotated(built: BuildResult) -> None:
    """Every element of the comp and the flat model has an SBO term and an ontology term."""
    missing: list[str] = []
    for path in (built.comp_path, built.flat_path):
        doc = libsbml.readSBMLFromFile(str(path))
        for model in _models(doc):
            for element in _elements(model):
                where = f"{path.name}:{model.getId()}:{element.getId()}"
                if not element.isSetSBOTerm():
                    missing.append(f"{where} (SBO)")
                if element.getNumCVTerms() == 0:
                    missing.append(f"{where} (ontology term)")
    assert not missing, missing[:20]


def test_model_type(built: BuildResult) -> None:
    """The model is a whole-body PBPK model of the human, perfusion- or permeability-limited."""
    doc = libsbml.readSBMLFromFile(str(built.flat_path))
    resources = set(_resources(doc.getModel()))
    assert _url(PBPKO.PHYSIOLOGICALLY_BASED_PHARMACOKINETIC_MODEL) in resources
    assert _url(PBPKO.WHOLE_BODY_PBPK) in resources
    assert "https://identifiers.org/taxonomy:9606" in resources
    # both definitions have a substance with a finite exchange
    assert _url(PBPKO.PERMEABILITY_LIMITED_PBPK) in resources
    assert _url(PBPKO.PERFUSION_LIMITED_MODEL) not in resources


def test_model_type_perfusion_limited() -> None:
    """The perfusion-limited term only when every substance is in perfusion mode."""
    perfusion = Substance("a", "a", 1.0, distribution=Kp(1.0, mode=Mode.PERFUSION))
    tissue = Substance("b", "b", 1.0, distribution=Kp(1.0, ftissue=1.0))
    vascular = Substance("c", "c", 1.0)
    perfusion_term = (BQB.IS, pbpko("PERFUSION_LIMITED_MODEL"))
    permeability_term = (BQB.IS, pbpko("PERMEABILITY_LIMITED_PBPK"))
    assert perfusion_term in model_annotations([perfusion, vascular])
    assert permeability_term not in model_annotations([perfusion, vascular])
    assert permeability_term in model_annotations([perfusion, tissue])
    assert perfusion_term not in model_annotations([perfusion, tissue])
    terms = model_annotations([vascular])
    assert perfusion_term not in terms
    assert permeability_term not in terms
    assert (BQB.IS, pbpko("PHYSIOLOGICALLY_BASED_PHARMACOKINETIC_MODEL")) in terms


def _resources(sbase: libsbml.SBase) -> list[str]:
    """The resources of the CV terms of an element."""
    result = []
    for i in range(sbase.getNumCVTerms()):
        term = sbase.getCVTerm(i)
        result += [term.getResourceURI(j) for j in range(term.getNumResources())]
    return result


_PBPKO_URL = re.compile(r"https://bioregistry\.io/PBPKO:\d{5}")


def _url(term: str) -> str:
    """The resource of a PBPKO term as the library writes it."""
    return f"https://bioregistry.io/{PBPKO[term].curie}"


def _table_resources() -> list[str]:
    """Every resource of the table of the kinds and the places."""
    resources = [
        pbpko(n) for k in KINDS.values() for n in (k.pbpko, *k.pools.values()) if n
    ]
    resources += [r for k in KINDS.values() for _, r in k.terms]
    for places in (ORGANS, POOLS, SEGMENTS, SITES):
        for place in places.values():
            resources += [r for _, r in place.terms]
            if place.compartment is not None:
                resources.append(pbpko(place.compartment))
    return sorted(set(resources))


def test_pbpko_resource_format(built: BuildResult) -> None:
    """The resource formats: identifiers.org `collection/term`, PBPKO by bioregistry.

    Every resource of a collection which identifiers.org registers is its
    `collection/term` (the format of `pkpdsim.annotations`), written as the
    identifiers.org url; identifiers.org does not register PBPKO, so a PBPKO
    term is `https://bioregistry.io/PBPKO:<5 digits>` of the pinned PBPKO,
    written as given (R48).
    """
    for resource in _table_resources():
        if "PBPKO" in resource:
            assert _PBPKO_URL.fullmatch(resource), resource
            assert resource.rsplit("/", 1)[1] in PBPKO, resource
            continue
        annotation = RDFAnnotation(qualifier=BQB.IS, resource=resource, validate=False)
        assert annotation.provider == ProviderType.IDENTIFIERS_ORG, resource
        assert f"{annotation.collection}/{annotation.term}" == resource
        assert annotation.check_miriam_term(), resource
    doc = libsbml.readSBMLFromFile(str(built.comp_path))
    pbpko_urls = set()
    for model in _models(doc):
        for element in [model, *_elements(model)]:
            for url in _resources(element):
                if "PBPKO" in url:
                    assert _PBPKO_URL.fullmatch(url), url
                    pbpko_urls.add(url)
                else:
                    assert url.startswith("https://identifiers.org/"), url
    assert pbpko_urls
    assert all(url.rsplit("/", 1)[1] in PBPKO for url in pbpko_urls)


def test_build_logs_no_annotation_errors(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """A build logs no error of the annotation (no MIRIAM namespace of pymetadata)."""
    with caplog.at_level(logging.WARNING):
        build(_variants(), tmp_path)
    errors = [
        f"{r.name}: {r.getMessage()}"
        for r in caplog.records
        if r.levelno >= logging.ERROR or "MIRIAM" in r.getMessage()
    ]
    assert not errors, errors[:5]


def test_requested_terms_documented() -> None:
    """Every PBPKO term the library lacks is listed in the documentation, requested upstream."""
    text = DOCS.read_text()
    assert "## Terms requested upstream" in text
    upstream = text.split("## Terms requested upstream", 1)[1]
    missing = [term for term in requested_terms() if term not in upstream]
    assert not missing, missing


def test_compartment_terms_resolved() -> None:
    """The compartment terms of `pkpdsim.annotations` have no open FIXME, every organ a UBERON term."""
    source = (
        Path(__file__).parents[2] / "src" / "pkpdsim" / "annotations.py"
    ).read_text()
    assert "fixme" not in source.lower()
    for key in ("fluid volume", "re", "rev", "kiv", "fov"):
        assert compartments[key], key
    for organ, place in ORGANS.items():
        assert any(r.startswith("uberon/") for _, r in place.terms), organ


def test_clearance_per_volume_is_no_clearance(built: BuildResult) -> None:
    """A first order clearance per litre of organ is not asserted as a clearance (R49).

    It is a kinetic constant and a `BIOCHEMICAL_PARAMETER` which is a version
    of the (intrinsic) clearance; the tubular share is a physiological
    parameter.
    """
    doc = libsbml.readSBMLFromFile(str(built.flat_path))
    model = doc.getModel()
    pid, version = (
        ("LI__INHX_k", "INTRINSIC_CLEARANCE_RATE")
        if model.getParameter("LI__INHX_k") is not None
        else ("KI__SECUP_k", "CLEARANCE_RATE")
    )
    k = model.getParameter(pid)
    assert k.getSBOTermID() == "SBO:0000009"  # kinetic constant
    terms = {
        (k.getCVTerm(i).getBiologicalQualifierType(), url)
        for i in range(k.getNumCVTerms())
        for url in _resources_of(k.getCVTerm(i))
    }
    assert (libsbml.BQB_IS, pbpko("BIOCHEMICAL_PARAMETER")) in terms
    assert (libsbml.BQB_IS_VERSION_OF, pbpko(version)) in terms
    assert (libsbml.BQB_IS, pbpko(version)) not in terms
    share = model.getParameter("KI__f_tubular")
    assert pbpko("PHYSIOLOGICAL_PARAMETER") in _resources(share)


def _resources_of(term: libsbml.CVTerm) -> list[str]:
    """The resources of one CV term."""
    return [term.getResourceURI(j) for j in range(term.getNumResources())]
