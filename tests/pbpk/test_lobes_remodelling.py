"""Tests of the remodelling of the lobes: the anatomy of the liver follows the simulation.

A lobe which loses its portal inflow atrophies and the lobes which keep theirs
hypertrophy, so the volume of a lobe is a state of the model and not a
constant of the individual. Two drivers move it: the flow driven growth law,
which predicts the remodelling from the perfusion per volume of the lobe, and
a prescribed course, which maps a measured time course onto the model.

The state of a species of the liver is its amount, so a lobe which grows
dilutes what it holds by itself. A lobe which shrinks would concentrate it
instead, which no tissue does: the tissue which is lost carries its share
away, so the remodelling emits that loss as a flux into a sink. The two
invariants are the first tests of the drivers.
"""

from pathlib import Path

import pytest

from pkpdsim.pbpk import (
    FirstOrder,
    Metabolism,
    Organ,
    PBPKDefinition,
    Route,
    Variants,
    ZonatedLiver,
    ids,
)
from pkpdsim.pbpk.builder import build
from pkpdsim.pbpk.organ import OrganComponent
from pkpdsim.pbpk.organs import Part
from pkpdsim.pbpk.remodelling import (
    FlowDriven,
    PrescribedVolume,
    lobes_dynamic,
    remodelling_sink,
)
from pkpdsim.pbpk.variants.liver_lobes import LobedLiver
from pkpdsim.sbml.values import Formula
from pkpdsim.species import Species

from .conftest import definition, individual
from .test_lobes_build import RAT_TOKENS

S = "drg"
DAY = 24 * 60.0


def _definition(liver: OrganComponent | None) -> PBPKDefinition:
    return definition(
        sid="m",
        species=Species.RAT,
        chemistry=[
            Metabolism(
                sid="DRG2M",
                organ=Organ.LIVER,
                substrate=S,
                products={"drg_m": 1.0},
                kinetics=FirstOrder(k=0.5),
                protein="CYP3A4",
            )
        ],
        routes=[Route.iv(S)],
        variants=Variants(organs={Organ.LIVER: liver})
        if liver is not None
        else Variants(),
    )


def _built(tmp_path: Path, liver: OrganComponent | None) -> Path:
    return build(_definition(liver), tmp_path).flat_path


def _rules(path: Path) -> tuple[set[str], set[str]]:
    """The ids with a rate rule and the ids with an assignment rule."""
    import libsbml

    document = libsbml.readSBMLFromFile(str(path))
    model = document.getModel()
    rates: set[str] = set()
    assignments: set[str] = set()
    for k in range(model.getNumRules()):
        rule = model.getRule(k)
        if rule.getTypeCode() == libsbml.SBML_RATE_RULE:
            rates.add(rule.getVariable())
        elif rule.getTypeCode() == libsbml.SBML_ASSIGNMENT_RULE:
            assignments.add(rule.getVariable())
    return rates, assignments


def test_a_static_anatomy_keeps_the_lobe_volume_a_constant(tmp_path: Path) -> None:
    rates, _ = _rules(_built(tmp_path, LobedLiver()))
    for token in RAT_TOKENS:
        assert ids.lobe_volume(token) not in rates


def test_the_flow_driven_driver_makes_the_lobe_volume_a_state(tmp_path: Path) -> None:
    rates, assignments = _rules(_built(tmp_path, LobedLiver(remodelling=FlowDriven())))
    for token in RAT_TOKENS:
        assert ids.lobe_volume(token) in rates
        assert ids.lobe_growth(token) in assignments
        # the share of the lobe leaves its volume behind and follows it
        assert ids.lobe_volume_fraction(token) in assignments


def _patency() -> set[str]:
    """The quantities an intervention changes, which are always dynamic."""
    return {
        name(token)
        for token in RAT_TOKENS
        for name in (
            ids.arterial_conductance,
            ids.portal_conductance,
            ids.lobe_present,
        )
    }


def test_lobes_dynamic_names_the_lobe_volumes() -> None:
    """With a driver the volumes move too, so everything which reads them follows."""
    liver = LobedLiver(remodelling=FlowDriven())
    assert set(lobes_dynamic(_definition(liver))) == _patency() | {
        ids.lobe_volume(t) for t in RAT_TOKENS
    }


def test_lobes_dynamic_names_the_patency_without_a_driver() -> None:
    """An intervention changes these whether or not the anatomy follows it.

    A ligation between two phases of a timecourse only moves the flows if the
    flows are assignment rules, which they are because these are the seeds of
    the dynamic quantities of the build; without a driver the volumes of the
    lobes stay constants of the individual.
    """
    assert set(lobes_dynamic(_definition(LobedLiver()))) == _patency()
    assert lobes_dynamic(_definition(None)) == ()


def test_a_site_volume_follows_the_lobe_volume(tmp_path: Path) -> None:
    """The compartments of a lobe leave their initial assignment and follow its volume."""
    liver = LobedLiver(zonation=ZonatedLiver(n_zones=2), remodelling=FlowDriven())
    _, assignments = _rules(_built(tmp_path, liver))
    site = ids.site("lll", ids.zone(1))
    assert site is not None
    assert ids.flat(Organ.LIVER, ids.compartment(Part.PLASMA, site)) in assignments


def test_the_flow_driven_driver_is_at_rest_at_baseline(tmp_path: Path) -> None:
    """Every lobe has the reference perfusion at baseline, so nothing grows or shrinks."""
    r = individual(_built(tmp_path, LobedLiver(remodelling=FlowDriven())))
    result = r.simulate(
        0.0, 5 * DAY, 51, selections=["time", *(ids.lobe_volume(t) for t in RAT_TOKENS)]
    )
    for k, token in enumerate(RAT_TOKENS, start=1):
        assert result[:, k][-1] == pytest.approx(result[:, k][0], rel=1e-6), token


def test_the_remnant_grows_back_to_the_reference_volume(tmp_path: Path) -> None:
    """A liver of one third of its mass is hyperperfused, so it grows back and stops there."""
    r = individual(_built(tmp_path, LobedLiver(remodelling=FlowDriven())))
    for token in RAT_TOKENS:
        r[ids.lobe_volume(token)] = r[ids.lobe_volume(token)] / 3.0
    result = r.simulate(0.0, 30 * DAY, 301, selections=["time", ids.VLI_LOBES])
    total = result[:, 1]
    assert total[0] == pytest.approx(r[ids.VLI_REF] / 3.0, rel=1e-6)
    assert total[-1] == pytest.approx(r[ids.VLI_REF], rel=1e-2)
    assert max(total) <= r[ids.VLI_REF] * (1.0 + 1e-6)


def test_a_growing_lobe_loses_nothing_to_the_sink(tmp_path: Path) -> None:
    """Growth dilutes the tissue by itself, so the remodelling flux must not fire.

    The state of a species of the liver is its amount: a lobe which grows
    spreads the same amount over the larger volume without any flux, which is
    what growing tissue does. Only atrophy carries content away, so the sink
    of a growing liver stays empty.
    """
    curve = Formula(
        f"{ids.VLI_REF} * 0.3 dimensionless * "
        f"(1 dimensionless + 1 dimensionless * time / 7200 min)",
        "l",
    )
    liver = LobedLiver(remodelling=PrescribedVolume({"lll": curve}))
    r = individual(_built(tmp_path, liver))
    r[ids.iv_dose(S)] = 10.0
    sink = remodelling_sink(S)
    result = r.simulate(
        0.0, 7200.0, 201, selections=["time", sink, ids.lobe_volume("lll")]
    )
    assert result[:, 2][-1] > 1.5 * result[:, 2][0]
    assert max(abs(result[:, 1])) == pytest.approx(0.0, abs=1e-12)


def test_an_atrophying_lobe_keeps_what_is_left_at_its_concentration(
    tmp_path: Path,
) -> None:
    """The tissue a lobe loses carries its share away, so shrinking does not concentrate."""
    curve = Formula(
        f"{ids.VLI_REF} * 0.3 dimensionless * "
        f"(1 dimensionless - 0.5 dimensionless * time / 7200 min)",
        "l",
    )
    liver = LobedLiver(remodelling=PrescribedVolume({"lll": curve}))
    r = individual(_built(tmp_path, liver))
    r[ids.iv_dose(S)] = 10.0
    cti = f"[{ids.flat(Organ.LIVER, ids.species(Part.TISSUE, S, 'lll'))}]"
    other = f"[{ids.flat(Organ.LIVER, ids.species(Part.TISSUE, S, 'rl'))}]"
    result = r.simulate(0.0, 7200.0, 201, selections=["time", cti, other])
    shrinking, intact = result[:, 1], result[:, 2]
    # the lobe halves while the intact lobe keeps its volume: without the
    # remodelling flux the shrinking lobe would hold twice the concentration
    assert shrinking[-1] == pytest.approx(intact[-1], rel=0.1)


def test_a_prescribed_volume_follows_its_curve(tmp_path: Path) -> None:
    """The lobe volume tracks the given course, with the lag of the tracking time."""
    curve = Formula(
        f"{ids.VLI_REF} * 0.3 dimensionless * "
        f"(1 dimensionless - 0.1 dimensionless * time / 7200 min)",
        "l",
    )
    liver = LobedLiver(remodelling=PrescribedVolume({"lll": curve}))
    r = individual(_built(tmp_path, liver))
    result = r.simulate(0.0, 7200.0, 101, selections=["time", ids.lobe_volume("lll")])
    volume = result[:, 1]
    assert volume[0] == pytest.approx(r[ids.VLI_REF] * 0.3, rel=1e-9)
    assert volume[-1] == pytest.approx(volume[0] * 0.9, rel=2e-3)


def test_a_prescribed_curve_of_an_unknown_lobe_is_refused(tmp_path: Path) -> None:
    """The lobes are resolved at build time, so a curve of no lobe is refused there."""
    curve = Formula("1 l", "l")
    with pytest.raises(ValueError, match="nope"):
        _built(tmp_path, LobedLiver(remodelling=PrescribedVolume({"nope": curve})))
