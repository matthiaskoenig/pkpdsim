"""Tests of the composition and the build of a PBPK model."""

import dataclasses
import logging
import math
import os
import time
from pathlib import Path

import libsbml
import pytest
import roadrunner
from sbmlutils.factory import Creator, create_model

from pkpdsim.pbpk import (
    Direction,
    Erythrocytes,
    Excretion,
    FirstOrder,
    FirstOrderDissolution,
    Inhibition,
    Kp,
    Membrane,
    Mode,
    Organ,
    Part,
    PBPKDefinition,
    Pool,
    Route,
    StandardOrgan,
    Substance,
    Transport,
    Variants,
    WeibullDissolution,
    ids,
)
from pkpdsim.pbpk.builder import (
    BuildResult,
    build,
    built_organs,
    compose,
    organ_component,
    organ_context,
)
from pkpdsim.pbpk.physiology import Sex, reference_tissue_volumes

from .conftest import (
    definition,
    dose_amount,
    individual,
    load,
    simple_drug,
    total_amount,
)

DAY = 24 * 60.0
"""One day (min)."""


def _closed_drug(**changes: object) -> Substance:
    """The simple drug without excretion."""
    return dataclasses.replace(simple_drug(), excretion=Excretion(), **changes)


def _closed(substance: Substance, sid: str = "closed") -> PBPKDefinition:
    """A definition of one substance without elimination, given intravenously."""
    return definition(
        sid=sid,
        substances=[substance],
        chemistry=[],
        routes=[Route.iv(substance.sid)],
    )


def _build(d: PBPKDefinition, path: Path) -> "roadrunner.RoadRunner":
    """Build and load a definition."""
    return load(build(d, path))


def _drift(r: "roadrunner.RoadRunner", sids: list[str], tend: float = DAY) -> float:
    """The relative change of the amount of substances (species, sinks, doses)."""

    def amount() -> float:
        return sum(total_amount(r, s) + dose_amount(r, s) for s in sids)

    before = amount()
    r.simulate(0, tend, 50)
    return abs(amount() - before) / before


@pytest.fixture(scope="module")
def simple(tmp_path_factory: pytest.TempPathFactory) -> BuildResult:
    """The build of the simple definition of the test helpers."""
    return build(definition(), tmp_path_factory.mktemp("simple"))


def _flat(result: BuildResult) -> libsbml.SBMLDocument:
    """The document of the flat model of a build (it owns the model)."""
    return libsbml.readSBMLFromFile(str(result.flat_path))


def test_build_validates_and_flattens(simple: BuildResult) -> None:
    """Both files exist, are valid SBML and the flat model loads in roadrunner."""
    assert simple.comp_path.name == "simple.xml"
    assert simple.flat_path.name == "simple_flat.xml"
    for path in (simple.comp_path, simple.flat_path):
        doc = libsbml.readSBMLFromFile(str(path))
        doc.checkConsistency()
        errors = [
            doc.getError(i).getMessage()
            for i in range(doc.getNumErrors())
            if doc.getError(i).getSeverity() >= libsbml.LIBSBML_SEV_ERROR
        ]
        assert errors == [], path.name
    r = load(simple)
    r.simulate(0, 10, 5)


def test_comp_model_without_warnings(tmp_path: Path) -> None:
    """The comp model has no validation warnings, units included."""
    result = create_model(compose(definition()), filepath=tmp_path / "simple.xml")
    assert result.validation is not None
    assert result.validation.error_count == 0
    assert result.validation.warning_count == 0


def test_flat_ids_follow_scheme(simple: BuildResult) -> None:
    """Organ ids are prefixed by the submodel, top ids are not; the outputs are rules.

    The outputs of the organs (`Cout_<o>_<s>`, the urinary excretion rate
    `Rurine_<s>`) are the rules of the organs under the ids of the top model.
    """
    doc = _flat(simple)
    model = doc.getModel()
    species = {s.getId() for s in model.getListOfSpecies()}
    assert {"LI__Cpl_drg", "LI__Cti_drg", "GU__Cpl_drg", "Cve_drg", "Car_drg"} <= (
        species
    )
    assert not any(s.startswith("LA__") for s in species)
    assert {"Aurine_drg", "Afeces_drg", "Cbile_drg"} <= species
    lumen = {
        ids.flat(Organ.SMALL_INTESTINE, ids.lumen("drg", seg))
        for seg in ids.GUT_SEGMENTS
    }
    assert lumen <= species
    rules = {
        r.getVariable(): libsbml.formulaToL3String(r.getMath())
        for r in model.getListOfRules()
        if r.isAssignment()
    }
    organs = built_organs(definition())
    outflows = {
        ids.organ_outflow(o, s): f"{ids.submodel(o)}__Cpl_{s}"
        for o in organs
        for s in ("drg", "drg_m")
    }
    assert {k: rules[k] for k in outflows} == outflows
    observers = {
        f(s)
        for s in ("drg", "drg_m")
        for f in (
            ids.arm,
            ids.arm_free,
            ids.arm_blood,
            ids.hepatic_vein,
            ids.hepatic_vein_free,
            ids.hepatic_vein_blood,
        )
    } | {f(p, s) for s in ("drg", "drg_m") for p in Pool for f in (ids.free, ids.blood)}
    unbound = {ids.flat(o, ids.unbound(Part.TISSUE, "drg")) for o in organs}
    assert rules[ids.flat(Organ.LIVER, ids.unbound(Part.TISSUE, "drg"))] == (
        "(fu_drg / Kp_drg) * LI__Cti_drg"
    )
    assert set(rules) - set(outflows) == observers | unbound | {
        ids.injection_rate("drg"),
        ids.urine_rate("drg"),
    }
    # the urinary excretion rate of the kidney replaces the parameter of the top
    assert rules[ids.urine_rate("drg")] == "KI__J_gfr_drg"
    initial = {
        a.getSymbol(): libsbml.formulaToL3String(a.getMath())
        for a in model.getListOfInitialAssignments()
    }
    assert initial[ids.fu("drg")].replace(" ", "") == (
        "1dimensionless/(1dimensionless+(fimp_albumin/fimp_plasma_protein)*"
        "(1dimensionless-fu_ref_drg)/fu_ref_drg)"
    )
    assert initial[ids.bp("drg")].replace(" ", "") == "1dimensionless-HCT"
    assert initial["GU__Vgu"].replace(" ", "") == "Vsi+Vla"
    assert initial["LI__Vpl"].replace(" ", "") == "Vli*f_vas_li*(1dimensionless-HCT)"
    assert model.getParameter("Cout_li_drg") is not None
    for pid in ("IVDOSE_drg", "PODOSE_drg", "Ri_drg", "ti_drg", "Mr_drg", "ka_drg"):
        assert model.getParameter(pid) is not None, pid
    assert model.getParameter(ids.absorption("drg")).getValue() == 1.0
    assert model.getParameter(ids.absorption("drg_m")).getValue() == 0.0


@pytest.mark.skipif(
    "PYTEST_XDIST_WORKER" in os.environ,
    reason="CPU time inflates under the contention of the parallel workers (R53); "
    "runs serially: `pytest -n 0`, the second command of tox",
)
def test_load_time(simple: BuildResult) -> None:
    """The flat model loads fast in sbmlsim.

    Every quantity of parameters only is an initial assignment, since
    roadrunner inlines assignment rules: with the deep rules of the physiology
    the simple model took 28.8 s to compile (LLJIT O2), now about 2 s. The CPU
    time of the compile grows with the contention of the parallel test run
    (11 s observed among the workers), so the test runs only serially, in a
    run without pytest-xdist (`pytest -n 0`, the second command of tox and so
    of continuous integration); the structural guard
    `test_rules_read_a_state` runs in every run.
    """
    start = time.process_time()
    load(simple)
    assert time.process_time() - start < 10.0


def rules_without_state(doc: libsbml.SBMLDocument) -> list[str]:
    """Get the assignment rules which read no state, also through other rules.

    A state is what changes during a simulation by itself: a species, the
    variable of a rate rule, `time`, or a variable parameter or compartment
    which no rule computes (an input changed during the simulation, a dynamic
    parameter of the definition). A reaction (its id in a formula is its
    flux, the urinary excretion rate `Rurine_<s>` sums fluxes) is read as
    its kinetic law. A rule which reaches none of them, directly or through
    the assignment rules and reactions it reads, depends on constants only and
    belongs into an initial assignment: roadrunner inlines assignment rules,
    and the deep chains of the physiology as rules made the compile of the
    simple model 14 times slower (R21). The rate of an intravenous injection
    `Kiv_<s>` is a rule by design (the dose protocol changes `ti_<s>` between
    doses) and is not reported.

    Args:
        doc: the flat model.

    Returns:
        The ids of the rules which read no state.
    """
    model = doc.getModel()
    rules = {r.getVariable(): r for r in model.getListOfRules() if r.isAssignment()}
    states = {s.getId() for s in model.getListOfSpecies()}
    states |= {r.getVariable() for r in model.getListOfRules() if r.isRate()}
    states |= {
        e.getId()
        for e in [*model.getListOfParameters(), *model.getListOfCompartments()]
        if not e.getConstant() and e.getId() not in rules
    }
    states.add("time")

    def names(rule: libsbml.Rule) -> set[str]:
        found: set[str] = set()
        nodes = [rule.getMath()]
        while nodes:
            node = nodes.pop()
            if node.getType() == libsbml.AST_NAME_TIME:
                found.add("time")
            elif node.isName():
                found.add(node.getName())
            nodes.extend(node.getChild(i) for i in range(node.getNumChildren()))
        return found

    reads = {sid: names(rule) for sid, rule in rules.items()}
    reads |= {
        r.getId(): names(r.getKineticLaw())
        for r in model.getListOfReactions()
        if r.isSetKineticLaw()
    }
    stateful = set(states)
    changed = True
    while changed:
        changed = False
        for sid, symbols in reads.items():
            if sid not in stateful and symbols & stateful:
                stateful.add(sid)
                changed = True
    return sorted(
        sid
        for sid in rules
        if sid not in stateful and not sid.startswith(ids.injection_rate(""))
    )


def test_rules_read_a_state(simple: BuildResult) -> None:
    """No assignment rule of the flat model depends on constants only (R53).

    The deterministic guard of the load time: the cause of the slow compile
    (assignment rules over constants which should be initial assignments) is
    absent, whatever the contention of the test run.
    """
    assert rules_without_state(_flat(simple)) == []


def test_rules_of_constants_are_found(simple: BuildResult) -> None:
    """The guard finds the slow pattern: the physiology as assignment rules."""
    doc = _flat(simple)
    model = doc.getModel()
    for pid in (ids.BMI, ids.BSA, ids.FFM, ids.CO):
        ia = model.getInitialAssignment(pid)
        rule = model.createAssignmentRule()
        rule.setVariable(pid)
        rule.setMath(ia.getMath().deepCopy())
        model.removeInitialAssignment(pid)
        model.getParameter(pid).setConstant(False)
    assert rules_without_state(doc) == sorted([ids.BMI, ids.BSA, ids.CO, ids.FFM])


def test_build_logs_no_sbmlutils_warning(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """A build logs no warning of sbmlutils: every element is complete.

    The definition has the variants (zonated liver, nephron kidney, segmented
    intestine with a Weibull tablet and an intravenous route) and a
    mechanism-based inhibitor (a turnover rate rule), so the elements of every
    builder are created. A definition with a `dynamic` parameter gets the hint
    of sbmlutils that nothing in the model changes the variable parameter: the
    experiment changes it during the simulation, the hint is expected there.
    """
    inhibitor = Substance(
        "inh",
        name="inhibitor",
        mr=350.0,
        distribution=Kp(1.0, ftissue=1.0),
        inhibits=(Inhibition("CYP3A4", "mechanism", ki=0.005, kinact=0.05),),
    )
    d = definition(
        substances=[
            simple_drug(),
            Substance("drg_m", name="metabolite", mr=316.0),
            inhibitor,
        ],
        variants=Variants(organs={Organ.LIVER: "zonated", Organ.KIDNEY: "nephron"}),
        routes=[
            Route.iv("drg"),
            Route.oral("drg", WeibullDissolution(t50=30.0, shape=1.5)),
        ],
    )
    with caplog.at_level(logging.WARNING, logger="sbmlutils"):
        build(d, tmp_path)
    warnings = [
        r.getMessage() for r in caplog.records if r.name.startswith("sbmlutils")
    ]
    assert warnings == []


def test_covariate_change_propagates(simple: BuildResult) -> None:
    """A change of an experiment before the simulation reaches organs and pools."""
    reference = individual(simple.flat_path)
    r = individual(simple.flat_path, BW=100.0, f_cardiac_function=0.5)
    vad, qmu = ids.organ_volume(Organ.ADIPOSE), ids.organ_flow(Organ.MUSCLE)
    assert r[vad] > reference[vad] + 10.0
    assert r[qmu] < 0.6 * reference[qmu]
    for local, factor in (("Vpl", 1.1), ("Vti", 1.5), ("PS_drg", 1.5)):
        flat = ids.flat(Organ.ADIPOSE, local)
        assert r[flat] > reference[flat] * factor, flat
    assert r["AD__Vpl"] == pytest.approx(
        r[vad] * r[ids.vascular_fraction(Organ.ADIPOSE)] * (1 - r[ids.HCT]), rel=1e-12
    )
    assert r["Vve_pl"] == pytest.approx(r["Vve"] * (1 - r[ids.HCT]), rel=1e-12)
    assert r[ids.fu("drg")] == reference[ids.fu("drg")]


def test_reference_tissue_volume(simple: BuildResult) -> None:
    """`Vti_total_ref_<s>` is the tissue volume of the reference individual."""
    r = load(simple)
    tissue = sum(
        r[ids.organ_volume(o)] * (1 - r[ids.vascular_fraction(o)]) for o in Organ
    )
    assert r[ids.tissue_volume_total_ref("drg")] == pytest.approx(tissue, rel=1e-9)
    assert tissue == pytest.approx(
        sum(reference_tissue_volumes(Sex.MALE).values()), rel=1e-9
    )


_ERYTHROCYTES = {
    "none": None,
    "equilibrium": Erythrocytes("equilibrium", bp=1.3),
    "kinetic": Erythrocytes("kinetic", kbc=2.0, k_exchange=1.0),
}


@pytest.mark.parametrize(
    ("mode", "ery"),
    [
        (mode, ery)
        for mode in Mode
        for ery in _ERYTHROCYTES
        if not (mode is Mode.PERFUSION and ery == "kinetic")
    ],
)
def test_mass_balance_closed_system(tmp_path: Path, mode: Mode, ery: str) -> None:
    """Without elimination the amount of an IV dose stays in the body."""
    drug = _closed_drug(
        distribution=Kp(2.0, mode=mode, ftissue=1.0),
        erythrocytes=_ERYTHROCYTES[ery],
    )
    r = _build(_closed(drug), tmp_path)
    r[ids.iv_dose("drg")] = 10.0
    assert _drift(r, ["drg"]) < 1e-6
    assert dose_amount(r, "drg") < 1e-12
    tissue = {
        Mode.TISSUE: ["MU__Cti_drg"],
        Mode.PERMEABILITY: ["MU__Cis_drg", "MU__Cic_drg"],
        Mode.PERFUSION: ["MU__Cor_drg"],
    }[mode]
    for sid in tissue:
        assert r[f"[{sid}]"] > 0, sid


def test_mass_balance_with_elimination(tmp_path: Path) -> None:
    """Drug, metabolite, sinks and doses together keep the amount of the doses.

    The liver metabolizes the drug, the metabolite leaves the liver into the
    bile, which drains into the lumen; the IV bolus, an IV infusion and the
    tablet are given.
    """
    bile = Transport(
        sid="DRGMBILE",
        organ=Organ.LIVER,
        substance="drg_m",
        direction=Direction.EFFLUX,
        membrane=Membrane.APICAL,
        kinetics=FirstOrder(k=0.1),
    )
    d = definition(sid="elimination")
    d = dataclasses.replace(d, chemistry=[*d.chemistry, bile])
    r = _build(d, tmp_path)
    r[ids.iv_dose("drg")] = 10.0
    r[ids.po_dose("drg")] = 10.0
    r[ids.infusion_rate("drg")] = 0.01
    assert dose_amount(r, "drg") == pytest.approx(20.0 / r[ids.mr("drg")])
    r.simulate(0, DAY, 50)
    mmole = (20.0 + 0.01 * DAY) / r[ids.mr("drg")]
    total = sum(total_amount(r, s) + dose_amount(r, s) for s in ("drg", "drg_m"))
    assert total == pytest.approx(mmole, rel=1e-6)
    assert total_amount(r, "drg_m") > 0.01 * mmole
    assert r[ids.feces("drg_m")] > 0
    assert r[ids.feces("drg")] > 0


def test_oral_doses_add_up(tmp_path: Path) -> None:
    """Tablet and solution: a dose added to the dose state is applied as well.

    `PKPDSimulationExperiment.add_doses` sets `PODOSE_<s> + dose`.
    """
    drug = _closed_drug()
    d = definition(
        sid="oral",
        substances=[drug],
        chemistry=[],
        routes=[Route.oral("drg"), Route.oral_solution("drg")],
    )
    r = _build(d, tmp_path)
    mr = r[ids.mr("drg")]
    r[ids.po_dose("drg")] = 10.0
    r[ids.sol_dose("drg")] = 5.0
    r.simulate(0, 60, 10)
    assert r[ids.sol_dose("drg")] < 1e-6
    r[ids.po_dose("drg")] = r[ids.po_dose("drg")] + 10.0
    r.simulate(60, DAY, 50)
    assert total_amount(r, "drg") + dose_amount(r, "drg") == pytest.approx(
        25.0 / mr, rel=1e-6
    )
    assert r[ids.feces("drg")] > 0
    assert r["[Cve_drg]"] > 0


def test_partial_organs(tmp_path: Path) -> None:
    """A substance entering liver and muscle has tissue only there."""
    organs = frozenset({Organ.LIVER, Organ.MUSCLE})
    drug = _closed_drug(distribution=Kp(2.0, organs=organs, ftissue=1.0))
    result = build(_closed(drug, "partial"), tmp_path)
    doc = _flat(result)
    tissue = sorted(
        s.getId() for s in doc.getModel().getListOfSpecies() if "__Cti_" in s.getId()
    )
    assert tissue == ["LI__Cti_drg", "MU__Cti_drg"]
    r = load(result)
    volumes = reference_tissue_volumes(Sex.MALE)
    assert r[ids.tissue_volume_total_ref("drg")] == pytest.approx(
        volumes[Organ.LIVER] + volumes[Organ.MUSCLE], rel=1e-12
    )
    r[ids.iv_dose("drg")] = 10.0
    assert _drift(r, ["drg"]) < 1e-6


def test_two_substances_independent(tmp_path: Path) -> None:
    """A second substance does not change the curve of the first one."""
    drug = _closed_drug()
    other = dataclasses.replace(_closed_drug(), sid="oth", name="other", mr=200.0)
    alone = _build(_closed(drug, "alone"), tmp_path)
    both = _build(
        definition(
            sid="both",
            substances=[drug, other],
            chemistry=[],
            routes=[Route.iv("drg"), Route.iv("oth")],
        ),
        tmp_path,
    )
    alone[ids.iv_dose("drg")] = 10.0
    both[ids.iv_dose("drg")] = 10.0
    both[ids.iv_dose("oth")] = 3.0
    times = [0, 1, 5, 30, 120, 600, DAY]
    a = alone.simulate(times=times, selections=["time", "[Cve_drg]"])
    b = both.simulate(times=times, selections=["time", "[Cve_drg]"])
    for ca, cb in zip(a[1:, 1], b[1:, 1], strict=True):
        assert cb == pytest.approx(ca, rel=1e-7)


def _ports(md: libsbml.Model) -> dict[str, str]:
    """The ports of a model definition: id -> SBO term."""
    plugin = md.getPlugin("comp")
    return {p.getId(): p.getSBOTermID() for p in plugin.getListOfPorts()}


def test_no_hand_written_replacements(simple: BuildResult) -> None:
    """Every replacement names a port of its submodel; every input and output is replaced."""
    doc = libsbml.readSBMLFromFile(str(simple.comp_path))
    model = doc.getModel()
    comp = model.getPlugin("comp")
    definitions = doc.getPlugin("comp")
    ports: dict[str, dict[str, str]] = {}
    for submodel in comp.getListOfSubmodels():
        md = definitions.getModelDefinition(submodel.getModelRef())
        ports[submodel.getId()] = _ports(md)
    d = definition()
    assert set(ports) == {ids.submodel(o) for o in built_organs(d)}
    assert ids.submodel(Organ.LARGE_INTESTINE) not in ports

    referenced: dict[str, set[str]] = {sub: set() for sub in ports}
    for element in model.getListOfAllElements():
        plugin = element.getPlugin("comp")
        if not isinstance(plugin, libsbml.CompSBasePlugin):
            continue
        refs = [
            plugin.getReplacedElement(i) for i in range(plugin.getNumReplacedElements())
        ]
        if plugin.isSetReplacedBy():
            refs.append(plugin.getReplacedBy())
        for ref in refs:
            sub, port = ref.getSubmodelRef(), ref.getPortRef()
            assert port in ports[sub], (sub, port)
            assert port not in referenced[sub], (sub, port)
            referenced[sub].add(port)

    input_output = {"SBO:0000600", "SBO:0000601"}
    for sub, organ_ports in ports.items():
        unreplaced = {p for p, sbo in organ_ports.items() if sbo in input_output}
        assert unreplaced <= referenced[sub], sub

    for organ in built_organs(d):
        _, boundary = organ_component(d, organ).build(organ_context(d, organ))
        expected = {f"{local}_port" for local in boundary.replacements} | {
            f"{local}_port" for local in boundary.outputs
        }
        assert expected == referenced[ids.submodel(organ)], organ


def test_variants() -> None:
    """A component as variant is used, an unknown variant name raises."""
    custom = StandardOrgan(name="custom")
    d = definition(variants=Variants(organs={Organ.LIVER: custom}))
    assert organ_component(d, Organ.LIVER) is custom
    model = compose(d)
    assert "custom_li" in {md.sid for md in model.model_definitions or []}

    unknown = definition(variants=Variants(organs={Organ.KIDNEY: "glomerular"}))
    with pytest.raises(ValueError, match=r"'glomerular'.*'ki'.*'nephron'"):
        compose(unknown)


def test_formulation_rate() -> None:
    """The first order dissolution of a formulation is the rate of the tablet."""
    d = definition(routes=[Route.oral("drg", FirstOrderDissolution(k=0.5))])
    parameters = {str(p.sid): p for p in compose(d).parameters or []}
    assert parameters[ids.dissolution("drg")].value == 0.5
    assert ids.iv_dose("drg") not in parameters


def test_route_given_twice_raises() -> None:
    """Two routes of one kind for one substance raise."""
    with pytest.raises(ValueError, match="twice"):
        compose(definition(routes=[Route.iv("drg"), Route.iv("drg")]))


def test_injection_rate_follows_injection_time(simple: BuildResult) -> None:
    """`Kiv_<s> = ln(2)/ti_<s>` is a rule: a dose protocol changes `ti_<s>` between doses."""
    r = load(simple)
    kiv, ti = ids.injection_rate("drg"), ids.injection_time("drg")
    assert r[kiv] == pytest.approx(math.log(2) / r[ti] * 60, rel=1e-12)
    r.simulate(0, 10, 2)
    r[ti] = 120.0
    r.simulate(10, 20, 2)
    assert r[kiv] == pytest.approx(math.log(2) / 120.0 * 60, rel=1e-12)


def test_creators_reach_the_flat_model(tmp_path: Path) -> None:
    """The creators of the definition are the model history of the comp and the flat model."""
    creator = Creator(
        familyName="König",
        givenName="Matthias",
        email="creator@example.org",
        organization="Humboldt-University Berlin",
        orcid="0000-0003-1725-179X",
    )
    result = build(definition(creators=[creator]), tmp_path)
    for path in (result.comp_path, result.flat_path):
        doc = libsbml.readSBMLFromFile(str(path))  # the document owns the history
        history = doc.getModel().getModelHistory()
        assert history is not None, path.name
        assert history.getNumCreators() == 1, path.name
        written = history.getCreator(0)
        assert (written.getFamilyName(), written.getGivenName()) == (
            "König",
            "Matthias",
        )
        assert written.getOrganisation() == "Humboldt-University Berlin"


def test_no_creators_no_history(simple: BuildResult) -> None:
    """Without creators the flat model has no creators in its history."""
    doc = _flat(simple)
    history = doc.getModel().getModelHistory()
    assert history is None or history.getNumCreators() == 0
