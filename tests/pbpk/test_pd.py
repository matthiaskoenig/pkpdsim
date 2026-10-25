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
from sbmlutils.factory import AssignmentRule, Model, Parameter, Port, PortType

from pkpdsim.pbpk import MODEL_UNITS, Organ, Part, PDComponent, Pool, U, ids
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
