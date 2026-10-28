"""Tests of the PD components: coupled to the body by their ports.

A PD component reads elements of the body through its input ports (the
element of the body replaces the input) and changes parameters of the body
through its output ports (the output replaces the parameter, which is
dynamic: everything downstream of it follows during the simulation).
"""

import dataclasses
from pathlib import Path

import libsbml
import numpy as np
import pytest
from sbmlutils.factory import (
    AssignmentRule,
    Model,
    ModelUnits,
    Parameter,
    Port,
    PortType,
    Species,
)

from pkpdsim.pbpk import (
    MODEL_UNITS,
    CustomReaction,
    Organ,
    Part,
    PDComponent,
    Pool,
    U,
    ids,
)
from pkpdsim.pbpk.builder import build, compose

from .conftest import cardiac_pd, definition, load


def _reader(sid: str = "eff") -> PDComponent:
    """A PD component which copies the unbound arm vein concentration of the drug."""
    model = Model(
        sid,
        units=U,
        model_units=MODEL_UNITS,
        objects=[
            Parameter(
                "C_drug",
                0.0,
                U.mM,
                constant=False,
                port=Port("C_drug_port", portType=PortType.INPUT_PORT),
            ),
            Parameter("E", None, U.mM, constant=False),
            AssignmentRule("E", "C_drug", U.mM),
        ],
    )
    return PDComponent(model=model, inputs={"C_drug": ids.arm_free("drg")}, outputs={})


def test_pd_reads_free_concentration(tmp_path: Path) -> None:
    """The input of a PD component is the unbound arm vein concentration of the body.

    The element of the body replaces the input port, so the component reads
    `Carm_free_drg` under its own id, at every time of the simulation.
    """
    r = load(build(definition(pd=[_reader()]), tmp_path))
    r[ids.iv_dose("drg")] = 10.0
    s = r.simulate(0, 600, 61, selections=["time", "eff__E", ids.arm_free("drg")])
    effect, free = s["eff__E"], s[ids.arm_free("drg")]
    assert free.max() > 0
    np.testing.assert_allclose(effect, free, rtol=1e-12, atol=0.0)


def test_pd_output_feeds_back(tmp_path: Path) -> None:
    """A PD output `f_cardiac` replaces `f_cardiac_function` and lowers `CO`.

    The replaced parameter is dynamic: the cardiac output and the organ
    flows (of the top model and inside the organs) are assignment rules of
    it and follow the effect of the drug during the simulation.
    """
    result = build(definition(pd=[cardiac_pd(ec50=1e-4)]), tmp_path)
    doc = libsbml.readSBMLFromFile(str(result.flat_path))
    model = doc.getModel()
    for pid in (ids.CO_BASE, ids.CO, ids.organ_flow(Organ.MUSCLE), ids.QH):
        assert model.getAssignmentRuleByVariable(pid) is not None, pid
    r = load(result)
    co_base = float(r[ids.CO])
    q_mu = float(r[ids.organ_flow(Organ.MUSCLE)])
    assert float(r[ids.F_CARDIAC_FUNCTION]) == pytest.approx(1.0)
    r[ids.iv_dose("drg")] = 10.0
    s = r.simulate(
        0,
        120,
        121,
        selections=[
            "time",
            ids.CO,
            ids.F_CARDIAC_FUNCTION,
            ids.arm_free("drg"),
            ids.organ_flow(Organ.MUSCLE),
        ],
    )
    f = s[ids.F_CARDIAC_FUNCTION]
    c = s[ids.arm_free("drg")]
    np.testing.assert_allclose(f, 1 - 0.5 * c / (1e-4 + c), rtol=1e-10)
    assert f.min() < 0.8
    # the cardiac output follows the cardiac function (no impairment)
    np.testing.assert_allclose(s[ids.CO], co_base * f, rtol=1e-10)
    np.testing.assert_allclose(s[ids.organ_flow(Organ.MUSCLE)], q_mu * f, rtol=1e-10)


@pytest.mark.parametrize(
    ("changes", "match"),
    [
        ({"inputs": {"C_unknown": ids.arm_free("drg")}}, "'C_unknown'"),
        ({"inputs": {"C_drug": "Carm_free_nothing"}}, "'Carm_free_nothing'"),
        ({"outputs": {"f_unknown": ids.F_CARDIAC_FUNCTION}}, "'f_unknown'"),
        ({"outputs": {"f_cardiac": "no_parameter"}}, "'no_parameter'"),
        ({"outputs": {"f_cardiac": ids.pool_species(Pool.VENOUS, "drg")}}, "'Cve_drg'"),
        ({"inputs": {"f_cardiac": ids.arm_free("drg")}}, "'f_cardiac'.*output"),
    ],
    ids=[
        "input port",
        "input element",
        "output port",
        "output element",
        "output species",
        "port type",
    ],
)
def test_pd_unknown_id_raises(changes: dict[str, dict[str, str]], match: str) -> None:
    """An unknown port or element of a PD component raises, naming the id."""
    pd = dataclasses.replace(cardiac_pd(), **changes)
    with pytest.raises(ValueError, match=match):
        compose(definition(pd=[pd]))


def test_pd_id_of_the_body_raises() -> None:
    """A PD component whose id is an id of the body raises, naming it."""
    pd = cardiac_pd(sid=ids.CO)
    with pytest.raises(ValueError, match=f"'{ids.CO}'"):
        compose(definition(pd=[pd]))


def test_pd_model_not_changed() -> None:
    """Composing does not change the model of the PD component (it is copied)."""
    pd = cardiac_pd()
    before = [(str(p.sid), p.metaId, p.annotations) for p in pd.model.parameters]
    compose(definition(pd=[pd]))
    after = [(str(p.sid), p.metaId, p.annotations) for p in pd.model.parameters]
    assert before == after


def test_pd_output_reaches_the_organs(tmp_path: Path) -> None:
    """The reduced cardiac output reaches the organs: the liver blood flow falls.

    A coupled PD output is not only an observer of the body: the flows which
    the organs read follow the cardiac function during the simulation, and
    so does the inflow of the drug into the liver.
    """
    r = load(build(definition(pd=[cardiac_pd()]), tmp_path))
    before = float(r[ids.QH])
    r[ids.iv_dose("drg")] = 10.0
    r.simulate(0, 10, 11)
    f = float(r[ids.F_CARDIAC_FUNCTION])
    assert f < 0.9
    assert float(r[ids.QH]) == pytest.approx(before * f, rel=1e-10)
    inflow = ids.flat(
        Organ.LIVER, ids.inflow_reaction(Pool.ARTERIAL, Part.PLASMA, "drg")
    )
    hct = float(r[ids.HCT])
    car = float(r[f"[{ids.pool_species(Pool.ARTERIAL, 'drg')}]"])
    shunt = float(r[ids.F_SHUNT])
    assert float(r[inflow]) == pytest.approx(
        (1 - shunt) * float(r[ids.QHA]) * (1 - hct) * car, rel=1e-8
    )
    assert float(r[ids.QHA]) < before


@pytest.mark.parametrize(
    "target",
    [ids.CO, ids.organ_flow(Organ.MUSCLE), ids.organ_outflow(Organ.LIVER, "drg")],
    ids=["initial assignment", "flow", "organ output"],
)
def test_pd_output_to_derived_parameter_raises(target: str) -> None:
    """An output which replaces a derived parameter raises, naming it.

    `CO` and the flows are initial assignments of the physiology,
    `Cout_li_drg` the output of the liver: the component changes the
    parameters they derive from (`f_cardiac_function`), not them.
    """
    pd = dataclasses.replace(cardiac_pd(), outputs={"f_cardiac": target})
    with pytest.raises(ValueError, match=f"'{target}'.*derived"):
        compose(definition(pd=[pd]))


def test_pd_output_to_dynamic_parameter(tmp_path: Path) -> None:
    """A parameter of `dynamic` may be the target of an output: it is dynamic either way."""
    d = definition(pd=[cardiac_pd()], dynamic=(ids.F_CARDIAC_FUNCTION,))
    model = libsbml.readSBMLFromFile(str(build(d, tmp_path).flat_path)).getModel()
    assert model.getAssignmentRuleByVariable(ids.CO) is not None
    assert model.getAssignmentRuleByVariable(ids.F_CARDIAC_FUNCTION) is not None


def test_pd_two_outputs_to_one_parameter_raises() -> None:
    """Two components cannot both replace one parameter of the body."""
    with pytest.raises(ValueError, match=f"'{ids.F_CARDIAC_FUNCTION}'.*derived"):
        compose(definition(pd=[cardiac_pd(), cardiac_pd(sid="cardio2")]))


def test_pd_model_units_raise() -> None:
    """A component in other model units (seconds) or without model units raises."""
    seconds = cardiac_pd()
    seconds.model.model_units = ModelUnits(
        time=U.second,
        extent=U.mmole,
        substance=U.mmole,
        length=U.metre,
        area=U.m2,
        volume=U.litre,
    )
    with pytest.raises(ValueError, match=r"'cardio'.*time s(econd)? \(library min\)"):
        compose(definition(pd=[seconds]))
    none = cardiac_pd()
    none.model.model_units = None
    with pytest.raises(ValueError, match="'cardio' has no model units"):
        compose(definition(pd=[none]))


def test_pd_id_of_an_organ_submodel_raises() -> None:
    """A component with the id of an organ submodel (`LI`) raises, naming it."""
    with pytest.raises(ValueError, match=r"'LI'.*submodel of an organ"):
        compose(definition(pd=[cardiac_pd(sid=ids.submodel(Organ.LIVER))]))


def _organ_reader(organ_species: str) -> PDComponent:
    """A PD component which doubles a concentration of an organ (`C_organ`)."""
    model = Model(
        "org",
        units=U,
        model_units=MODEL_UNITS,
        objects=[
            Parameter(
                "C_organ",
                0.0,
                U.mM,
                constant=False,
                port=Port("C_organ_port", portType=PortType.INPUT_PORT),
            ),
            Parameter("E", None, U.mM, constant=False),
            AssignmentRule("E", "2 dimensionless * C_organ", U.mM),
        ],
    )
    return PDComponent(model=model, inputs={"C_organ": organ_species}, outputs={})


def test_pd_reads_a_species_of_an_organ(tmp_path: Path) -> None:
    """A PD input names a species of an organ by its flat id (`MU__Cti_drg`).

    The organ exposes the species as an output, the top model has it as the
    parameter `C<o>_<part>_<s>` (`Cmu_ti_drg`), which the component reads; the
    species keeps its id in the flat model. An imaging signal of the liver
    reads the concentration in the liver cells this way.
    """
    local = ids.species(Part.TISSUE, "drg")
    species = ids.flat(Organ.MUSCLE, local)
    top = ids.organ_species(Organ.MUSCLE, local)
    assert top == "Cmu_ti_drg"
    result = build(definition(pd=[_organ_reader(species)]), tmp_path)
    model = libsbml.readSBMLFromFile(str(result.flat_path)).getModel()
    assert model.getSpecies(species) is not None
    rule = model.getAssignmentRuleByVariable(top)
    assert rule is not None
    assert libsbml.formulaToL3String(rule.getMath()) == species
    r = load(result)
    r[ids.iv_dose("drg")] = 10.0
    s = r.simulate(0, 600, 61, selections=["time", "org__E", top, f"[{species}]"])
    assert s[f"[{species}]"].max() > 0
    np.testing.assert_allclose(s[top], s[f"[{species}]"], rtol=1e-12, atol=0.0)
    np.testing.assert_allclose(s["org__E"], 2 * s[top], rtol=1e-12, atol=0.0)


def _organ_definitions(model: Model) -> dict[str, str]:
    """Get the organ model definitions of a composed model, written as SBML.

    The ports, parameters and rules of every model definition but the PD
    components, by the id of the definition; the metaids and the history,
    which differ between two builds, are left out.
    """
    result: dict[str, str] = {}
    for md in model.model_definitions or []:
        if str(md.sid) in {"eff", "org"}:
            continue
        ports = sorted(
            f"{e.sid}:{getattr(e, 'port', None)}"
            for name in ("compartments", "species", "parameters", "reactions")
            for e in getattr(md, name) or []
            if getattr(e, "port", None)
        )
        parameters = sorted(str(p.sid) for p in md.parameters or [])
        rules = sorted(f"{r.variable}={r.value}" for r in md.rules or [])
        result[str(md.sid)] = repr((ports, parameters, rules))
    return result


def test_pd_without_an_organ_input_is_unchanged(tmp_path: Path) -> None:
    """A component which reads the body only adds no output to an organ.

    The organ model definitions of a definition with a component which reads
    the top model are those of the definition without it: no organ gains a
    port, a parameter or a rule; and the flat model has no `C<o>_<part>_<s>`.
    """
    plain = _organ_definitions(compose(definition()))
    with_reader = _organ_definitions(compose(definition(pd=[_reader()])))
    assert plain
    assert with_reader == plain
    result = build(definition(pd=[_reader()]), tmp_path)
    model = libsbml.readSBMLFromFile(str(result.flat_path)).getModel()
    organ_species = {
        ids.organ_species(organ, ids.species(part, "drg"))
        for organ in Organ
        for part in Part
    }
    assert not organ_species & {p.getId() for p in model.getListOfParameters()}
    # an organ input adds exactly its output to the organ which it reads
    with_organ = _organ_definitions(
        compose(definition(pd=[_organ_reader("MU__Cti_drg")]))
    )
    changed = {sid for sid in plain if with_organ[sid] != plain[sid]}
    assert changed == {sid for sid in plain if sid.endswith("_mu")}


@pytest.mark.parametrize(
    ("organ_species", "match"),
    [
        ("MU__Cnothing_drg", r"'MU__Cnothing_drg'.*no species of the organ"),
        ("MU__Vti", r"'MU__Vti'.*no species of the organ"),
        ("XX__Cti_drg", r"reads the unknown id 'XX__Cti_drg'"),
    ],
    ids=["unknown species", "compartment", "unknown submodel"],
)
def test_pd_input_of_an_organ_raises(organ_species: str, match: str) -> None:
    """An input which names no species of a built organ raises, naming the id."""
    with pytest.raises(ValueError, match=match):
        compose(definition(pd=[_organ_reader(organ_species)]))


@pytest.mark.parametrize(
    ("local", "amount", "match"),
    [
        ("Apool_drg", True, r"'MU__Apool_drg', an amount \(hasOnlySubstanceUnits\)"),
        ("Xpool_drg", False, r"'MU__Xpool_drg', a species .* no concentration"),
    ],
    ids=["amount", "no concentration id"],
)
def test_pd_input_of_a_species_which_is_no_concentration_raises(
    local: str, amount: bool, match: str
) -> None:
    """A species of a custom reaction which is an amount or no `C...` id is named."""
    pool = Species(
        local,
        "Vti",
        initialConcentration=0.0,
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=amount,
    )
    reaction = CustomReaction(
        "POOL",
        Organ.MUSCLE,
        Part.TISSUE,
        f"Cti_drg -> {local}",
        "0.001 l_per_min * Cti_drg",
        species=[pool],
    )
    d = definition(pd=[_organ_reader(f"MU__{local}")])
    d = dataclasses.replace(d, chemistry=[*d.chemistry, reaction])
    with pytest.raises(ValueError, match=match):
        compose(d)
