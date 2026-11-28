"""Tests of the annotation of every element the library creates (#63).

Two definitions together use every feature of the library: the organ
variants (zonated liver, nephron kidney, segmented intestine; the standard
organs with the single lumen), kinetic and equilibrium erythrocytes, a Kp
with an override, a Kp method, the flow-limited organ, active transport on
both membranes, first order kinetics per volume, competitive and
mechanism-based inhibition, induction, the impairment, the observers (of the
library, the urinary excretion rate and of a definition), a PD component and
every route and formulation.
Every compartment, species, parameter and reaction of the comp model and of
the flat model has an SBO term and at least one ontology term.
"""

import dataclasses
import logging
import re
from pathlib import Path

import libsbml
import pytest
from pymetadata.core.annotation import ProviderType, RDFAnnotation
from pymetadata.ontologies.pbpko import PBPKO
from sbmlutils.factory import Parameter, Reaction, Species, UncertParameter, UncertSpan
from sbmlutils.factory import Uncertainty as SbmlUncertainty
from sbmlutils.metadata import BQB, SBO

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
    Observer,
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
    annotate,
    model_annotations,
    pbpko,
    requested_terms,
)
from pkpdsim.pbpk.builder import BuildResult, build
from pkpdsim.pbpk.distribution import constant
from pkpdsim.pbpk.values import Value
from pkpdsim.sbml.domains import (
    BOOLEAN,
    DOMAIN,
    FRACTION,
    NONNEGATIVE,
    REAL,
    Domain,
    Values,
    domain_uncertainty,
)
from pkpdsim.sbml.parameters import value_parameter

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
            Metabolism("INHX", Organ.LIVER, "inh", {}, FirstOrder(0.01)),
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
                FirstOrder(0.01, per_volume=False),
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
            Route.im("drg"),
            Route.sc("drg"),
            Route.rectal("drg"),
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
                FirstOrder(0.1),
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
        observers=[
            Observer("Aurine_total", "Aurine_drg + Aurine_sec", U.mmole, "urine")
        ],
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
    written as given.
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


@pytest.mark.uncached_build
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


@pytest.mark.uncached_build
def test_an_ontology_term_of_a_substance_is_written_as_its_url(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """A term of pymetadata in the annotations of a substance is its url.

    The species of the substance carry the url which resolves the term (COB,
    which identifiers.org does not register, by bioregistry.io), not the bare
    identifier `COB_0000013`, and the build logs no annotation error.
    """
    d = _variants()
    substances = [
        dataclasses.replace(s, annotations=[*s.annotations, (BQB.IS, PBPKO.MOLECULE)])
        if s.sid == "drg"
        else s
        for s in d.substances
    ]
    with caplog.at_level(logging.WARNING):
        result = build(dataclasses.replace(d, substances=substances), tmp_path)
    errors = [
        f"{r.name}: {r.getMessage()}"
        for r in caplog.records
        if r.levelno >= logging.ERROR
    ]
    assert not errors, errors[:5]

    doc = libsbml.readSBMLFromFile(str(result.comp_path))
    urls = {
        url
        for model in _models(doc)
        for species in model.getListOfSpecies()
        for url in _resources(species)
    }
    assert "https://bioregistry.io/COB:0000013" in urls
    assert all(url.startswith("https://") for url in urls), sorted(urls)[:5]


def test_requested_terms_documented() -> None:
    """Every PBPKO term the library lacks is listed in the documentation, requested upstream."""
    text = DOCS.read_text(encoding="utf-8")
    assert "## Terms requested upstream" in text
    upstream = text.split("## Terms requested upstream", 1)[1]
    missing = [term for term in requested_terms() if term not in upstream]
    assert not missing, missing


def test_compartment_terms_resolved() -> None:
    """The compartment terms of `pkpdsim.annotations` have no open FIXME, every organ a UBERON term."""
    source = (
        Path(__file__).parents[2] / "src" / "pkpdsim" / "annotations.py"
    ).read_text(encoding="utf-8")
    assert "fixme" not in source.lower()
    for key in ("fluid volume", "re", "rev", "kiv", "fov"):
        assert compartments[key], key
    for organ, place in ORGANS.items():
        assert any(r.startswith("uberon/") for _, r in place.terms), organ


def test_clearance_per_volume_is_no_clearance(built: BuildResult) -> None:
    """A first order clearance per litre of organ is not asserted as a clearance.

    It is a kinetic constant and a `BIOCHEMICAL_PARAMETER` which is a version
    of the (intrinsic) clearance; the clearance of a whole organ
    (`per_volume=False`) is a clearance; the tubular share is a physiological
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
    whole = model.getParameter("KI__DRGSEC_k")
    if whole is not None:
        assert whole.getSBOTermID() == "SBO:0000009"
        assert pbpko("CLEARANCE_RATE") in _resources(whole)
        assert pbpko("BIOCHEMICAL_PARAMETER") not in _resources(whole)
    share = model.getParameter("KI__f_tubular")
    assert pbpko("PHYSIOLOGICAL_PARAMETER") in _resources(share)


def _resources_of(term: libsbml.CVTerm) -> list[str]:
    """The resources of one CV term."""
    return [term.getResourceURI(j) for j in range(term.getNumResources())]


def test_outputs_of_the_model_are_output_parameters(built: BuildResult) -> None:
    """The urinary excretion rate and an observer of the definition are output parameters.

    PBPKO's urinary rate is a biochemical parameter (an input): the rate is
    a version of it, not asserted as it.
    """
    doc = libsbml.readSBMLFromFile(str(built.comp_path))
    model = doc.getModel()
    rate = model.getParameter("Rurine_drg")
    assert _url_of("OUTPUT_PARAMETER") in _resources(rate)
    terms = {
        (rate.getCVTerm(i).getBiologicalQualifierType(), url)
        for i in range(rate.getNumCVTerms())
        for url in _resources_of(rate.getCVTerm(i))
    }
    assert (libsbml.BQB_IS_VERSION_OF, pbpko("URINARY_RATE")) in terms
    assert (libsbml.BQB_IS, pbpko("URINARY_RATE")) not in terms
    if model.getId() == "standard":
        observer = model.getParameter("Aurine_total")
        assert _url_of("OUTPUT_PARAMETER") in _resources(observer)
        assert observer.getSBOTermID() == "SBO:0000002"


def _url_of(name: str) -> str:
    """The resource of a PBPKO term by its name."""
    return pbpko(name)


def _uncertainty_names(element: Parameter | Reaction) -> list[str | None]:
    """The names of the uncertainties of an element."""
    return [u.name for u in element.uncertainties or []]


def test_annotate_writes_the_domain_of_the_kind() -> None:
    p = annotate(Parameter("fabs_drg", 0.94, U.dimensionless), "fraction_absorbed")
    assert _uncertainty_names(p) == [DOMAIN]


def test_annotate_writes_no_domain_for_a_kind_without_one() -> None:
    r = annotate(Reaction(sid="R", equation="A -> B"), "transport", organ=Organ.LIVER)
    assert not _uncertainty_names(r)


def test_annotate_keeps_the_domain_the_element_already_carries() -> None:
    # the precedence of a ModelParameter over the kind
    p = Parameter(
        "fm_x", 0.3, U.dimensionless, uncertainties=[domain_uncertainty(FRACTION)]
    )
    annotate(p, "organ_volume", organ=Organ.LIVER)  # the domain is NONNEGATIVE
    (span,) = (p.uncertainties or [])[0].uncertParameters or []
    assert isinstance(span, UncertSpan)
    assert _uncertainty_names(p) == [DOMAIN]
    assert (span.valueLower, span.valueUpper) == (0.0, 1.0)


def test_annotate_appends_the_domain_after_the_uncertainty_of_a_source() -> None:
    p = value_parameter(
        "Vli",
        Value(
            1.5,
            U.litre,
            "fit",
            uncertainties=[
                SbmlUncertainty(
                    uncertParameters=[
                        UncertParameter(
                            type=libsbml.DISTRIB_UNCERTTYPE_STANDARDDEVIATION,
                            value=0.1,
                            unit=U.litre,
                        )
                    ],
                    name="source fit",
                )
            ],
        ),
        U.litre,
        "liver volume",
        "SBO:0000002",
        "",
    )
    annotate(p, "organ_volume", organ=Organ.LIVER)
    assert _uncertainty_names(p) == ["source fit", DOMAIN]


#: the kinds which state no domain: the elements of a reaction hold no number, and
#: the two kinds of a constant which the `ModelParameter` states the domain of
_WITHOUT_DOMAIN = [
    "absorption",
    "bile_flow",
    "biliary_excretion",
    "blood_flow",
    "custom_parameter",
    "custom_reaction",
    "definition_parameter",
    "depot_absorption_im",
    "depot_absorption_inhaled",
    "depot_absorption_rectal",
    "depot_absorption_sc",
    "depot_application_im",
    "depot_application_inhaled",
    "depot_application_rectal",
    "depot_application_sc",
    "depot_loss_im",
    "depot_loss_inhaled",
    "depot_loss_rectal",
    "depot_loss_sc",
    "dissolution_reaction",
    "exchange",
    "fecal_excretion",
    "filtration",
    "gut_metabolism",
    "iv_application",
    "lumen_conversion",
    "metabolism",
    "metabolism_product",
    "metabolism_saturable",
    "oral_application",
    "reabsorption",
    "remodelling",
    "renal_secretion",
    "transit",
    "transport",
    "tubular_flow",
]

#: the kinds whose element is a share, a probability or a switch of 0 and 1
_FRACTIONS = [
    "arm_weight",
    "depot_fraction_im",
    "depot_fraction_inhaled",
    "depot_fraction_rectal",
    "depot_fraction_sc",
    "fed",
    "flow_fraction",
    "fraction_absorbed",
    "fraction_reabsorbed",
    "fu",
    "hematocrit",
    "inhaled_swallowed_fraction",
    "interstitial_fraction",
    "lobe_patency",
    "lobe_supply_fraction",
    "placement_weight",
    "pool_fraction",
    "product_fraction",
    "rectal_portal_fraction",
    "shunt",
    "tissue_composition",
    "tubular_share",
    "tubule_fraction",
    "vascular_fraction",
    "zone_fraction",
]

#: the kinds which take both signs
_REAL = [
    "body_composition_change",
    "body_weight_change",
    "lobe_growth",
    "observer",
    "steady_rate",
]

_REACTION_SBO = frozenset(
    {
        SBO.TRANSPORT_REACTION,
        SBO.BIOCHEMICAL_REACTION,
        SBO.DEGRADATION,
        SBO.BIOCHEMICAL_OR_TRANSPORT_REACTION,
    }
)


def test_the_kinds_which_state_no_domain() -> None:
    # a kind added to the table has to state a domain or join the list on purpose
    assert sorted(k for k, v in KINDS.items() if v.domain is None) == _WITHOUT_DOMAIN


def test_the_kinds_without_domain_are_the_reactions_and_two_constants() -> None:
    reactions = {k for k, v in KINDS.items() if v.sbo in _REACTION_SBO}
    assert set(_WITHOUT_DOMAIN) - reactions == {
        "custom_parameter",
        "definition_parameter",
    }
    assert reactions <= set(_WITHOUT_DOMAIN)


def test_the_fractions_of_the_library() -> None:
    assert sorted(k for k, v in KINDS.items() if v.domain == FRACTION) == _FRACTIONS


def test_the_kinds_which_take_both_signs() -> None:
    assert sorted(k for k, v in KINDS.items() if v.domain == REAL) == _REAL


def test_the_kinds_with_an_interval_of_their_own() -> None:
    assert KINDS["child_pugh"].domain == Domain(0.0, 3.0)
    assert KINDS["heart_failure"].domain == Domain(0.0, 4.0)
    assert KINDS["clock_time"].domain == Domain(0.0, 24.0)


def test_a_factor_is_not_bounded_by_one() -> None:
    # induction and hypertrophy are real states of the model
    for kind in ("scaling_factor", "age_factor", "pregnancy", "remnant_fraction"):
        assert KINDS[kind].domain == NONNEGATIVE


def test_a_ratio_of_concentrations_is_not_a_share() -> None:
    assert KINDS["concentration_ratio"].domain == NONNEGATIVE
    assert _uncertainty_names(constant("r_protein_int_pls")) == [DOMAIN]
    (span,) = (constant("r_protein_int_pls").uncertainties or [])[0].uncertParameters
    assert isinstance(span, UncertSpan)
    assert span.valueUpper == float("inf")
    assert KINDS["tissue_composition"].domain == FRACTION


def test_a_gross_fraction_and_an_abundance_are_not_shares() -> None:
    # `Fa` is gross (an absorbed drug which returns to the lumen counts twice) and
    # a relative abundance has a mean of 1 over the zones
    assert KINDS["absorbed_fraction"].domain == NONNEGATIVE
    assert KINDS["relative_abundance"].domain == NONNEGATIVE
    assert KINDS["placement_weight"].domain == FRACTION


def test_the_switch_of_renal_impairment_is_a_boolean() -> None:
    assert KINDS["ckd"].domain == BOOLEAN


def test_the_flags_of_the_library_are_booleans() -> None:
    assert sorted(k for k, v in KINDS.items() if v.domain == BOOLEAN) == [
        "ckd",
        "lobe_present",
    ]


def test_the_severities_stay_continuous_between_their_classes() -> None:
    # the model interpolates every factor linearly between the classes and a
    # scan or the severity of a MELD score takes values like 1.5, so an integer
    # domain would forbid what the model computes
    for kind in ("child_pugh", "heart_failure"):
        domain = KINDS[kind].domain
        assert domain is not None
        assert domain.values is Values.CONTINUOUS


def test_the_bounds_of_every_integer_kind_are_integral() -> None:
    for kind in KINDS.values():
        d = kind.domain
        if d is not None and d.values is not Values.CONTINUOUS:
            assert d.lower == int(d.lower) and d.upper == int(d.upper)


def test_the_prandial_state_and_the_patency_stay_continuous() -> None:
    # both interpolate between their two states (a light meal, a partial stenosis)
    assert KINDS["fed"].domain == FRACTION
    assert KINDS["lobe_patency"].domain == FRACTION
