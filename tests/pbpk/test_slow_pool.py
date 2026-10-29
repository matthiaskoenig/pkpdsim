"""The slow (deep) pool of a tissue: `SlowPool`."""

import dataclasses
from pathlib import Path
from typing import Any

import libsbml
import pytest

from pkpdsim.pbpk import (
    Binding,
    Excretion,
    Kp,
    Mode,
    Observer,
    Organ,
    Part,
    PBPKDefinition,
    Route,
    SlowPool,
    Substance,
    Transport,
    U,
    Value,
    Variants,
    build,
    ids,
)
from pkpdsim.pbpk.chemistry import Direction, FirstOrder, Membrane

from .conftest import (
    assert_reference_sbml,
    dose_amount,
    individual,
    reference_definition,
    total_amount,
)

K = 2e-3  # the exchange of the pool per litre of tissue (1/min)
KP = 1.5  # the pool to the tissue at equilibrium


def drug(**fields: Any) -> Substance:  # noqa: ANN401
    """Get a drug in the merged tissue of every organ, filtered into the urine."""
    values: dict[str, Any] = {
        "sid": "drg",
        "name": "drug",
        "mr": 300.0,
        "binding": Binding(fu=0.5, protein="albumin"),
        "distribution": Kp(2.0, ftissue=1.0),
        "excretion": Excretion(urine=True),
    }
    values.update(fields)
    return Substance(**values)


def pool(**fields: Any) -> SlowPool:  # noqa: ANN401
    """Get the slow pool of the drug in muscle and skin."""
    values: dict[str, Any] = {
        "substance": "drg",
        "organs": frozenset({Organ.MUSCLE, Organ.SKIN}),
        "k": K,
        "kp": KP,
    }
    values.update(fields)
    return SlowPool(**values)


def with_pool(
    *pools: SlowPool,
    substance: Substance | None = None,
    **fields: Any,  # noqa: ANN401
) -> PBPKDefinition:
    """Get the definition of the drug with slow pools, an iv route."""
    values: dict[str, Any] = {
        "sid": "slow",
        "substances": [substance or drug()],
        "routes": [Route.iv("drg")],
        "slow_pools": list(pools),
    }
    values.update(fields)
    return PBPKDefinition(**values)


def reaction_formula(path: Path, rid: str) -> str:
    """Get the rate law of a reaction of an SBML file as infix."""
    doc = libsbml.readSBMLFromFile(str(path))
    reaction = doc.getModel().getReaction(rid)
    assert reaction is not None, rid
    return libsbml.formulaToL3String(reaction.getKineticLaw().getMath())


def test_every_organ_of_the_pool_has_its_species_and_exchange(tmp_path: Path) -> None:
    result = build(with_pool(pool()), tmp_path)
    doc = libsbml.readSBMLFromFile(str(result.flat_path))
    model = doc.getModel()
    for organ in (Organ.MUSCLE, Organ.SKIN):
        sid = ids.flat(organ, ids.slow_species("drg"))
        species = model.getSpecies(sid)
        assert species is not None, sid
        assert species.getCompartment() == ids.flat(organ, "Vti")
        assert not species.getHasOnlySubstanceUnits()
        rid = ids.flat(organ, ids.slow_reaction("drg"))
        assert reaction_formula(result.flat_path, rid) == (
            f"kslow_drg * {ids.flat(organ, 'Vti')} * (Kpslow_drg * "
            f"{ids.flat(organ, 'Cti_drg')} - {sid})"
        )
    assert model.getSpecies(ids.flat(Organ.LIVER, ids.slow_species("drg"))) is None
    assert model.getParameter(ids.kslow("drg")).getValue() == K
    assert model.getParameter(ids.kp_slow("drg")).getValue() == KP
    assert ids.slow_species("drg") == "Cslow_drg"
    assert ids.kslow("drg") == "kslow_drg"
    assert ids.kp_slow("drg") == "Kpslow_drg"
    assert ids.slow_reaction("drg") == "J_slow_drg"


def test_the_pool_keeps_the_mass(tmp_path: Path) -> None:
    """Dose = body + urine at every time; the pool takes up a part of it."""
    r = individual(build(with_pool(pool()), tmp_path).flat_path)
    r.integrator.absolute_tolerance = 1e-14
    r.integrator.relative_tolerance = 1e-10
    r[ids.iv_dose("drg")] = 100.0
    start = total_amount(r, "drg") + dose_amount(r, "drg")
    pools = [ids.flat(o, ids.slow_species("drg")) for o in (Organ.MUSCLE, Organ.SKIN)]
    for end in (10.0, 600.0, 6000.0):
        r.simulate(0, end, 3)
        assert total_amount(r, "drg") + dose_amount(r, "drg") == pytest.approx(
            start, rel=1e-8
        )
        assert all(r[f"[{sid}]"] > 0 for sid in pools)
        r.reset()
        r[ids.iv_dose("drg")] = 100.0


def test_the_exchange_scales_with_the_tissue_volume(tmp_path: Path) -> None:
    """`k*Vti*(kp*Cti - Cslow)`: the rate follows the tissue of the individual."""
    path = build(with_pool(pool()), tmp_path).flat_path
    cti = ids.flat(Organ.MUSCLE, ids.species(Part.TISSUE, "drg"))
    vti = ids.flat(Organ.MUSCLE, "Vti")
    rid = ids.flat(Organ.MUSCLE, ids.slow_reaction("drg"))
    rates, volumes = [], []
    for bw in (73.0, 110.0):
        r = individual(path, **{ids.BW: bw})
        r[f"[{cti}]"] = 1.0
        rates.append(r[rid])
        volumes.append(r[vti])
        assert r[rid] == pytest.approx(K * r[vti] * KP * 1.0, rel=1e-12)
    assert volumes[1] > 1.2 * volumes[0]
    assert rates[1] / rates[0] == pytest.approx(volumes[1] / volumes[0], rel=1e-12)


def test_the_pool_reaches_kp_times_the_tissue(tmp_path: Path) -> None:
    """Without elimination the pool reaches `Cslow = kp*Cti` in every organ of it."""
    d = with_pool(pool(k=0.05), substance=drug(excretion=Excretion()))
    r = individual(build(d, tmp_path).flat_path)
    r.integrator.absolute_tolerance = 1e-14
    r.integrator.relative_tolerance = 1e-10
    r[ids.iv_dose("drg")] = 10.0
    r.simulate(0, 20000, 3)
    for organ in (Organ.MUSCLE, Organ.SKIN):
        cti = r[f"[{ids.flat(organ, ids.species(Part.TISSUE, 'drg'))}]"]
        cslow = r[f"[{ids.flat(organ, ids.slow_species('drg'))}]"]
        assert cti > 0
        assert cslow == pytest.approx(KP * cti, rel=1e-6)


def test_an_endogenous_pool_starts_at_its_equilibrium(tmp_path: Path) -> None:
    """With a baseline the pool starts at `Kpslow*Kp*baseline` and stays there."""
    baseline = 0.01
    d = with_pool(pool(), substance=drug(excretion=Excretion(), baseline=baseline))
    r = individual(build(d, tmp_path).flat_path)
    sid = f"[{ids.flat(Organ.MUSCLE, ids.slow_species('drg'))}]"
    assert r[sid] == pytest.approx(KP * 2.0 * baseline, rel=1e-12)
    r.simulate(0, 1000, 3)
    assert r[sid] == pytest.approx(KP * 2.0 * baseline, rel=1e-8)


def test_the_constants_take_values_with_provenance(tmp_path: Path) -> None:
    k = Value(K, U.per_min, source="fit", notes="the fitted exchange")
    result = build(with_pool(pool(k=k)), tmp_path)
    doc = libsbml.readSBMLFromFile(str(result.flat_path))
    p = doc.getModel().getParameter(ids.kslow("drg"))
    assert p.getValue() == K
    assert "Fitted value." in p.getNotesString()
    assert "the fitted exchange" in p.getNotesString()


def test_an_observer_reads_the_pool(tmp_path: Path) -> None:
    """The pool of an organ is exposed as `C<o>_slow_<s>` (one exposure rule)."""
    flat = ids.flat(Organ.MUSCLE, ids.slow_species("drg"))
    observer = Observer("Aslow_mu_drg", f"{flat} * 1 litre", U.mmole, "pool amount")
    result = build(with_pool(pool(), observers=[observer]), tmp_path)
    r = individual(result.flat_path)
    assert ids.organ_species(Organ.MUSCLE, ids.slow_species("drg")) == "Cmu_slow_drg"
    assert "Cmu_slow_drg" in r.model.getGlobalParameterIds()


def test_the_species_of_the_substance_include_the_pool() -> None:
    assert ids.species_of("drg", ["MU__Cslow_drg", "MU__Cslow_drg_m", "Cslow_x"]) == [
        "MU__Cslow_drg"
    ]


def test_without_a_pool_the_sbml_is_unchanged(tmp_path: Path) -> None:
    assert reference_definition().slow_pools == ()
    assert_reference_sbml(build(reference_definition(), tmp_path))


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ({"substance": "xyz"}, r"slow pool of 'xyz': unknown substance"),
        ({"organs": frozenset()}, r"slow pool of 'drg': no organs"),
        ({"k": 0.0}, r"slow pool of 'drg': k must be positive"),
        ({"kp": -1.0}, r"slow pool of 'drg': kp must be positive"),
        (
            {"k": Value(K, U.l_per_min)},
            r"slow pool of 'drg': k must be in 1/min",
        ),
        (
            {"organs": frozenset({Organ.LARGE_INTESTINE})},
            r"slow pool of 'drg': organ 'la' is not built",
        ),
        (
            {"organs": frozenset({Organ.SMALL_INTESTINE})},
            r"slow pool of 'drg': organ 'si' is the variant 'segmented'",
        ),
    ],
)
def test_an_invalid_pool_is_refused(change: dict[str, Any], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        with_pool(pool(**change)).validate()


def test_a_pool_needs_the_merged_tissue() -> None:
    """The pool exchanges with `Cti`: the flow-limited and the full form have none."""
    for dist in (
        Kp(2.0, mode=Mode.PERFUSION),
        Kp(2.0, mode=Mode.PERMEABILITY, ftissue=1.0),
        Kp(2.0, ftissue=1.0, organs=frozenset({Organ.LIVER})),
    ):
        with pytest.raises(ValueError, match=r"organ 'mu' has no merged tissue"):
            with_pool(pool(), substance=drug(distribution=dist)).validate()
    transport = Transport(
        "DRGIM",
        Organ.MUSCLE,
        "drg",
        Direction.INFLUX,
        Membrane.BASOLATERAL,
        FirstOrder(k=1.0),
    )
    with pytest.raises(ValueError, match=r"organ 'mu' has no merged tissue"):
        with_pool(pool(), chemistry=[transport]).validate()


def test_a_pool_in_a_variant_is_refused() -> None:
    d = with_pool(
        pool(organs=frozenset({Organ.LIVER})),
        variants=Variants(organs={Organ.LIVER: "zonated"}),
    )
    with pytest.raises(ValueError, match=r"organ 'li' is the variant 'zonated'"):
        d.validate()


def test_one_pool_per_substance() -> None:
    with pytest.raises(ValueError, match=r"several slow pools of \['drg'\]"):
        with_pool(pool(), pool(organs=frozenset({Organ.ADIPOSE}))).validate()


def test_combine_takes_the_pools_of_the_definitions() -> None:
    first = with_pool(pool())
    other = with_pool(substance=drug(), sid="other")
    combined = PBPKDefinition.combine(first, other, sid="both")
    assert list(combined.slow_pools) == [pool()]
    changed = with_pool(pool(k=2 * K), sid="changed")
    with pytest.raises(ValueError, match="the slow pool 'drg' is defined differently"):
        PBPKDefinition.combine(first, changed, sid="both")


def test_the_organs_are_normalised() -> None:
    p = pool(organs=["mu", "skin"])
    assert p.organs == frozenset({Organ.MUSCLE, Organ.SKIN})
    with pytest.raises(ValueError, match=r"SlowPool\.organs"):
        pool(organs=["nowhere"])
    assert dataclasses.replace(p, k=1.0).k == 1.0
