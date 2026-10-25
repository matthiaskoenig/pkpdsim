"""Tests of the intestine variants and the oral routes.

The segmented intestine (default) and the single lumen replace small and
large intestine with the submodel `GU`: lumen segments with transit,
dissolution and absorption, the mucosa for a substance with chemistry of the
intestine, the bile into the duodenum, the feces.
"""

import dataclasses
import math
from pathlib import Path
from typing import TYPE_CHECKING

import libsbml
import numpy as np
import pytest

from pkpdsim.pbpk import (
    DoseDependentDissolution,
    Excretion,
    FirstOrder,
    Induction,
    Inhibition,
    Kp,
    Metabolism,
    Organ,
    Part,
    PBPKDefinition,
    Pool,
    Route,
    SingleLumen,
    Substance,
    Transport,
    Variants,
    WeibullDissolution,
    build,
    compose,
    ids,
)
from pkpdsim.pbpk.chemistry import Direction, Membrane
from pkpdsim.pbpk.impairment import GASTRIC_EMPTYING, SI_TRANSIT
from pkpdsim.pbpk.variants import lookup, names

from .conftest import definition, dose_amount, load, simple_drug, total_amount

if TYPE_CHECKING:
    import roadrunner

S, M = "drg", "drg_m"
DAY = 24 * 60.0
SMALL = ("du", "jej1", "jej2", "ile1", "ile2")
"""The segments of the small intestine."""
COLON = ("cae", "asc", "tra", "des", "sig", "rec")
"""The segments of the large intestine."""


def _gu(local: str) -> str:
    """The flat id of a local id of the intestine."""
    return ids.flat(Organ.SMALL_INTESTINE, local)


def _drug() -> Substance:
    """The simple drug without renal excretion (eliminated by the chemistry only)."""
    return dataclasses.replace(simple_drug(), excretion=Excretion())


def _metabolite() -> Substance:
    """The metabolite: no distribution, renal excretion."""
    return Substance(
        sid=M, name="metabolite", mr=316.0, excretion=Excretion(urine=True)
    )


def _gut_metabolism(k: float = 0.5) -> Metabolism:
    """A CYP3A4 metabolism of the small intestine, first order (l/min)."""
    return Metabolism(
        sid="DRG2M_GU",
        organ=Organ.SMALL_INTESTINE,
        substrate=S,
        products={M: 1.0},
        kinetics=FirstOrder(k=k),
        protein="CYP3A4",
    )


def _efflux(k: float = 0.5) -> Transport:
    """The apical P-gp efflux of the small intestine, first order (l/min)."""
    return Transport(
        sid="PGP_DRG",
        organ=Organ.SMALL_INTESTINE,
        substance=S,
        direction=Direction.EFFLUX,
        membrane=Membrane.APICAL,
        kinetics=FirstOrder(k=k),
        protein="ABCB1",
    )


def _definition(sid: str, **overrides: object) -> PBPKDefinition:
    """The drug and its metabolite with the routes of the drug, no chemistry."""
    fields: dict[str, object] = {
        "sid": sid,
        "substances": [_drug(), _metabolite()],
        "chemistry": [],
        "routes": [Route.oral(S), Route.oral_solution(S)],
    }
    fields.update(overrides)
    return definition(**fields)


def _build(tmp_path: Path, d: PBPKDefinition) -> "roadrunner.RoadRunner":
    """Build a definition and load its flat model."""
    return load(build(d, tmp_path))


def _flat_species(r: "roadrunner.RoadRunner") -> set[str]:
    """The species ids of the flat model."""
    return set(r.model.getFloatingSpeciesIds()) | set(r.model.getBoundarySpeciesIds())


def _lumen_amount(
    r: "roadrunner.RoadRunner", sid: str, segments: tuple[str, ...]
) -> float:
    """The dissolved and undissolved amount (mmole) of a substance in lumen segments."""
    species = _flat_species(r)
    total = 0.0
    for seg in segments:
        total += r[f"[{_gu(ids.lumen(sid, seg))}]"] * r[_gu(ids.lumen_volume(seg))]
        if _gu(ids.solid(sid, seg)) in species:
            total += r[_gu(ids.solid(sid, seg))]
    return total


def _fraction_absorbed(r: "roadrunner.RoadRunner", dose_mg: float) -> float:
    """The fraction of an oral solution absorbed after a day: 1 - (lumen + feces)/dose."""
    r[ids.sol_dose(S)] = dose_mg
    r.simulate(0, DAY, 50)
    left = _lumen_amount(r, S, ids.GUT_SEGMENTS) + r[ids.feces(S)]
    return 1.0 - left / (dose_mg / r[ids.mr(S)])


def test_cat_fraction_absorbed(tmp_path: Path) -> None:
    """Uniform absorption in the small intestine reproduces the CAT model.

    `Fa = 1 - (1 + ka*T_si/n)^-n` (Yu & Amidon 1999) with the n = 5 segments
    of the small intestine of equal transit time, no absorption in the
    colon (and none in the stomach) and no gut chemistry.
    """
    r = _build(tmp_path, _definition("cat", routes=[Route.oral_solution(S)]))
    t_si, ka = 200.0, 0.6
    for seg in SMALL:
        r[_gu(ids.transit_time(seg))] = t_si / 5
        assert r[_gu(ids.area_factor(seg))] == 1.0
    for seg in COLON:
        r[_gu(ids.area_factor(seg))] = 0.0
    r[ids.absorption(S)] = ka
    fa = _fraction_absorbed(r, 10.0)
    expected = 1 - (1 + ka / 60 * t_si / 5) ** -5
    assert fa == pytest.approx(expected, abs=1e-3)


def test_mass_balance_oral(tmp_path: Path) -> None:
    """Dose = body + lumen (dissolved and undissolved) + mucosa + feces + urine.

    Tablet and solution of the drug with gut wall metabolism (mucosa),
    biliary excretion of the metabolite (bile into the duodenum, transit into
    the feces) and the liver metabolism.
    """
    bile = Transport(
        sid="DRGMBILE",
        organ=Organ.LIVER,
        substance=M,
        direction=Direction.EFFLUX,
        membrane=Membrane.APICAL,
        kinetics=FirstOrder(k=0.1),
    )
    base = definition()
    d = _definition("balance", chemistry=[*base.chemistry, _gut_metabolism(), bile])
    r = _build(tmp_path, d)
    mr = r[ids.mr(S)]
    r[ids.po_dose(S)] = 10.0
    r[ids.sol_dose(S)] = 5.0
    r.simulate(0, 60, 10)
    assert r[_gu(ids.solid(S, "du"))] > 0
    assert r[f"[{_gu(ids.mucosa(S, 'du'))}]"] > 0
    r.simulate(60, 3 * DAY, 50)
    total = sum(total_amount(r, s) + dose_amount(r, s) for s in (S, M))
    assert total == pytest.approx(15.0 / mr, rel=1e-6)
    assert r[ids.feces(S)] > 0
    assert r[ids.feces(M)] > 0
    assert r[ids.urine(M)] > 0


def test_single_lumen_variant(tmp_path: Path) -> None:
    """The single lumen: absorption `ka`, of which `fabs` enters the body.

    The fraction of a solution which reaches the feces is
    `1 - fabs*ka/(ka + 1/T)` with the lumen transit time `T`.
    """
    d = _definition(
        "single",
        routes=[Route.oral_solution(S)],
        variants=Variants(organs={Organ.SMALL_INTESTINE: SingleLumen(fabs=0.5)}),
    )
    r = _build(tmp_path, d)
    species = _flat_species(r)
    assert _gu(ids.lumen(S)) in species
    assert not any(_gu(ids.lumen(S, seg)) in species for seg in ids.GUT_SEGMENTS)
    assert r[ids.fraction_absorbed(S)] == 0.5
    mr = r[ids.mr(S)]
    r[ids.sol_dose(S)] = 10.0
    r.simulate(0, 2 * DAY, 50)
    ka, t = r[ids.absorption(S)] / 60, r[_gu(ids.T_LUMEN)]
    expected = 1 - 0.5 * ka / (ka + 1 / t)
    assert r[ids.feces(S)] / (10.0 / mr) == pytest.approx(expected, rel=1e-6)
    assert total_amount(r, S) + total_amount(r, M) + dose_amount(r, S) == (
        pytest.approx(10.0 / mr, rel=1e-6)
    )
    named = compose(
        dataclasses.replace(
            d, variants=Variants(organs={Organ.SMALL_INTESTINE: "single_lumen"})
        )
    )
    assert "single_lumen_si" in {md.sid for md in named.model_definitions or []}


def test_gut_wall_metabolism_reduces_fg(tmp_path: Path) -> None:
    """Gut wall metabolism in the mucosa reduces the drug reaching the body.

    Without chemistry of the intestine the drug has no mucosa species and is
    absorbed into the plasma of the wall; with it, the mucosa metabolizes it
    on its first pass, the metabolite forms in the mucosa.
    """
    exposure = {}
    for gut in (False, True):
        d = _definition(
            f"fg{int(gut)}",
            routes=[Route.oral_solution(S)],
            chemistry=[_gut_metabolism()] if gut else [],
        )
        r = _build(tmp_path, d)
        mucosa = {_gu(ids.mucosa(S, seg)) for seg in SMALL + COLON}
        assert bool(mucosa & _flat_species(r)) is gut
        assert (mucosa <= _flat_species(r)) is gut
        r[ids.sol_dose(S)] = 10.0
        result = r.simulate(
            0, DAY, 200, selections=["time", f"[{ids.pool_species(Pool.VENOUS, S)}]"]
        )
        exposure[gut] = float(np.trapezoid(result[:, 1], result[:, 0]))
        if gut:
            assert r[f"[{_gu(ids.mucosa(M, 'du'))}]"] > 0
    assert exposure[True] < 0.9 * exposure[False]


def test_placement_by_expression(tmp_path: Path) -> None:
    """CYP3A4 is placed by expression times mucosa volume, the shares sum to 1.

    The colon has no CYP3A4 in the data: no reaction there; the duodenum has
    the highest expression per volume.
    """
    d = _definition("placement", chemistry=[_gut_metabolism()])
    result = build(d, tmp_path)
    r = load(result)
    weights = {seg: _gu(ids.placement_weight("cyp3a4", seg)) for seg in SMALL}
    total = sum(r[w] for w in weights.values())
    assert total == pytest.approx(1.0, rel=1e-12)
    doc = libsbml.readSBMLFromFile(str(result.flat_path))
    reactions = {x.getId() for x in doc.getModel().getListOfReactions()}
    assert {_gu(ids.segment_process("DRG2M_GU", seg)) for seg in SMALL} <= reactions
    assert not {_gu(ids.segment_process("DRG2M_GU", seg)) for seg in COLON} & reactions
    per_volume = {
        seg: r[weights[seg]] / r[_gu(ids.mucosa_volume(seg))] for seg in SMALL
    }
    assert max(per_volume, key=lambda k: per_volume[k]) == "du"


def test_efflux_returns_to_lumen(tmp_path: Path) -> None:
    """Apical P-gp efflux moves the drug from the mucosa back into the lumen of the segment."""
    fa = {}
    for efflux in (False, True):
        d = _definition(
            f"efflux{int(efflux)}",
            routes=[Route.oral_solution(S)],
            chemistry=[_efflux()] if efflux else [],
        )
        result = build(d, tmp_path)
        if efflux:
            doc = libsbml.readSBMLFromFile(str(result.flat_path))
            reaction = doc.getModel().getReaction(
                _gu(ids.segment_process("PGP_DRG", "jej1"))
            )
            assert [x.getSpecies() for x in reaction.getListOfReactants()] == [
                _gu(ids.mucosa(S, "jej1"))
            ]
            assert [x.getSpecies() for x in reaction.getListOfProducts()] == [
                _gu(ids.lumen(S, "jej1"))
            ]
        fa[efflux] = _fraction_absorbed(load(result), 10.0)
    assert fa[True] < fa[False]


def test_solution_skips_dissolution(tmp_path: Path) -> None:
    """A solution enters the stomach dissolved; a tablet dissolves and has undissolved drug."""
    r = _build(tmp_path, _definition("dissolution"))
    stomach = f"[{_gu(ids.lumen(S, 'st'))}]"
    r[ids.sol_dose(S)] = 10.0
    r.simulate(0, 10, 5)
    solution = r[stomach]
    assert r[_gu(ids.solid(S, "du"))] == 0.0
    r.resetAll()
    r[ids.po_dose(S)] = 10.0
    r.simulate(0, 10, 5)
    assert r[stomach] < 0.5 * solution
    assert r[_gu(ids.solid(S, "du"))] > 0


def test_fed_slows_gastric_emptying(tmp_path: Path) -> None:
    """`FED = 1` slows the gastric emptying: more of a solution stays in the stomach."""
    r = _build(tmp_path, _definition("fed", routes=[Route.oral_solution(S)]))
    stomach = f"[{_gu(ids.lumen(S, 'st'))}]"
    left = {}
    for fed in (0.0, 1.0):
        r.resetAll()
        r[ids.FED] = fed
        r[ids.sol_dose(S)] = 10.0
        r.simulate(0, 30, 5)
        left[fed] = r[stomach]
    assert left[1.0] > 1.5 * left[0.0]


def _with_transit_factors(path: Path, factor: float) -> "roadrunner.RoadRunner":
    """Load a flat model with the CKD factors of the transit times set to a value.

    The factors are initial assignments of the impairment (1 in every stage
    of the data, R40); the test replaces them by a constant.
    """
    from sbmlsim.model.model_roadrunner import load_roadrunner

    doc = libsbml.readSBMLFromFile(str(path))
    model = doc.getModel()
    for pid in (GASTRIC_EMPTYING, SI_TRANSIT):
        model.removeInitialAssignment(pid)
        model.getParameter(pid).setValue(factor)
    changed = path.with_name(f"{path.stem}_transit{factor:g}.xml")
    libsbml.writeSBMLToFile(doc, str(changed))
    return load_roadrunner(changed)


def test_longer_transit_factor_slows_transit(tmp_path: Path) -> None:
    """A transit time factor of 2 doubles the transit times: the drug moves slower.

    The factor multiplies the time `T*factor` of the stomach and the small
    intestine, so the transit rate is `V*C/(T*factor)` (R55, the precedence
    bug gave `V*C/T*factor`). After a solution the stomach keeps more drug and
    less reaches the colon.
    """
    result = build(_definition("transit", routes=[Route.oral_solution(S)]), tmp_path)
    stomach = f"[{_gu(ids.lumen(S, 'st'))}]"
    colon = f"[{_gu(ids.lumen(S, 'cae'))}]"
    left, arrived = {}, {}
    for factor in (1.0, 2.0):
        r = _with_transit_factors(result.flat_path, factor)
        r[ids.sol_dose(S)] = 10.0
        r.simulate(0, 30, 5)
        left[factor] = r[stomach]
        r.simulate(30, 240, 5)
        arrived[factor] = r[colon]
    assert left[2.0] > 1.5 * left[1.0]
    assert arrived[2.0] < arrived[1.0]


def test_add_doses_compatible(tmp_path: Path) -> None:
    """Two tablet doses `PODOSE_drg + 10.0` add up, the second on the rest of the first."""
    r = _build(tmp_path, _definition("doses", routes=[Route.oral(S)]))
    mr = r[ids.mr(S)]
    r[ids.po_dose(S)] = r[ids.po_dose(S)] + 10.0
    r.simulate(0, 15, 5)
    rest = r[ids.po_dose(S)]
    assert 0 < rest < 10.0
    r[ids.po_dose(S)] = r[ids.po_dose(S)] + 10.0
    assert r[ids.po_dose(S)] == pytest.approx(rest + 10.0)
    r.simulate(15, DAY, 50)
    total = total_amount(r, S) + total_amount(r, M) + dose_amount(r, S)
    assert total == pytest.approx(20.0 / mr, rel=1e-6)


@pytest.mark.parametrize("shape", [1.0, 2.0])
def test_weibull_dissolution(shape: float, tmp_path: Path) -> None:
    """The undissolved fraction follows `exp(-ln2*((t - lag)/t50)^shape)` after the lag.

    The time since the dose is the mean age of the undissolved drug, exact
    for one dose: nothing dissolves during the lag, half at `lag + t50`.
    """
    lag, t50 = 20.0, 40.0
    route = Route.oral(S, WeibullDissolution(t50=t50, shape=shape, lag=lag))
    r = _build(tmp_path, _definition(f"weibull{int(shape)}", routes=[route]))
    mr = r[ids.mr(S)]
    solids = [_gu(ids.solid(S, seg)) for seg in ids.GUT_SEGMENTS[1:]]

    def undissolved() -> float:
        return r[ids.po_dose(S)] / mr + sum(r[a] for a in solids)

    r[ids.po_dose(S)] = 10.0
    r.simulate(0, lag - 1, 5)
    assert undissolved() == pytest.approx(10.0 / mr, rel=1e-9)
    r.simulate(lag - 1, lag + t50, 20)
    assert undissolved() / (10.0 / mr) == pytest.approx(0.5, abs=1e-4)
    r.simulate(lag + t50, lag + 2 * t50, 20)
    expected = math.exp(-math.log(2) * 2**shape)
    assert undissolved() / (10.0 / mr) == pytest.approx(expected, abs=1e-4)


def test_dose_dependent_dissolution(tmp_path: Path) -> None:
    """The dissolution saturates with the undissolved dose: a large dose dissolves slower."""
    route = Route.oral(S, DoseDependentDissolution(k=2.0, dose50=50.0))
    r = _build(tmp_path, _definition("dosedep", routes=[route]))
    rate = _gu(ids.dissolution_rate(S))
    fraction = {}
    for dose in (1.0, 500.0):
        r.resetAll()
        r[ids.po_dose(S)] = dose
        r.simulate(0, 1e-3, 2)
        assert r[rate] == pytest.approx(2.0 / 60 * 50.0 / (50.0 + dose), rel=1e-3)
        r.simulate(1e-3, 30, 5)
        fraction[dose] = 1 - r[ids.po_dose(S)] / dose
    assert fraction[500.0] < fraction[1.0]


def test_variant_registry() -> None:
    """The intestine has the standard, segmented and single lumen variants; unknown names raise."""
    assert names(Organ.SMALL_INTESTINE) == ["segmented", "single_lumen", "standard"]
    assert names(Organ.LIVER) == ["standard", "zonated"]
    with pytest.raises(ValueError, match=r"'lobular'.*'li'"):
        lookup(Organ.LIVER, "lobular")
    with pytest.raises(ValueError, match=r"'lobular'.*'li'"):
        compose(definition(variants=Variants(organs={Organ.LIVER: "lobular"})))


def test_standard_intestine(tmp_path: Path) -> None:
    """Standard organs as intestine: no lumen, oral routes raise, the bile goes into the feces."""
    standard = Variants(organs={Organ.SMALL_INTESTINE: "standard"})
    with pytest.raises(ValueError, match="Oral routes"):
        compose(definition(variants=standard))
    with pytest.raises(ValueError, match="PGP_DRG"):
        compose(
            definition(variants=standard, chemistry=[_efflux()], routes=[Route.iv(S)])
        )
    efflux = Transport(
        sid="DRGBILE",
        organ=Organ.LIVER,
        substance=S,
        direction=Direction.EFFLUX,
        membrane=Membrane.APICAL,
        kinetics=FirstOrder(k=0.5),
    )
    d = definition(
        sid="standard_gut",
        substances=[_drug()],
        chemistry=[efflux],
        routes=[Route.iv(S)],
        variants=standard,
    )
    r = _build(tmp_path, d)
    species = _flat_species(r)
    assert ids.flat(Organ.LARGE_INTESTINE, ids.species(Part.PLASMA, S)) in species
    assert not any("Clumen" in x for x in species)
    r[ids.iv_dose(S)] = 10.0
    r.simulate(0, DAY, 20)
    assert r[ids.feces(S)] > 0
    assert total_amount(r, S) + dose_amount(r, S) == pytest.approx(
        10.0 / r[ids.mr(S)], rel=1e-7
    )


def test_large_intestine_variant_with_gut_raises() -> None:
    """The segmented intestine replaces the large intestine: a variant of it raises."""
    d = definition(variants=Variants(organs={Organ.LARGE_INTESTINE: "standard"}))
    with pytest.raises(ValueError, match="replaces"):
        compose(d)


def test_chemistry_of_large_intestine_raises() -> None:
    """A process of the large intestine is refused, the small intestine is the gut organ."""
    gut = dataclasses.replace(_gut_metabolism(), organ=Organ.LARGE_INTESTINE)
    with pytest.raises(ValueError, match="large intestine"):
        compose(_definition("large", chemistry=[gut]))


def test_gut_inhibitor_acts_in_mucosa() -> None:
    """An inhibitor of a gut protein has mucosa species and acts from the mucosa of the segment."""
    inhibitor = Substance(
        sid="inh",
        name="inhibitor",
        mr=500.0,
        inhibits=(Inhibition("CYP3A4", "mechanism", 1e-3, kinact=1e-2),),
    )
    d = _definition(
        "inhibited",
        substances=[_drug(), _metabolite(), inhibitor],
        chemistry=[_gut_metabolism()],
        routes=[Route.oral(S), Route.oral("inh")],
    )
    gu = next(
        md for md in compose(d).model_definitions or [] if md.sid == "segmented_si"
    )
    species = {str(sp.sid) for sp in gu.species or []}
    assert {ids.mucosa("inh", seg) for seg in SMALL + COLON} <= species
    reaction = next(
        x for x in gu.reactions or [] if x.sid == ids.segment_process("DRG2M_GU", "du")
    )
    assert reaction.formula is not None
    assert ids.mucosa("inh", "du") in str(reaction.formula.math)
    assert ids.enzyme("CYP3A4") in {str(p.sid) for p in gu.parameters or []}


def test_large_intestine_kp_raises() -> None:
    """A Kp of the large intestine is refused: the wall has the Kp of the small intestine."""
    drug = dataclasses.replace(
        _drug(),
        distribution=Kp(2.0, ftissue=1.0, overrides={Organ.LARGE_INTESTINE: 3.0}),
    )
    with pytest.raises(ValueError, match="large intestine"):
        compose(_definition("kp", substances=[drug, _metabolite()]))


def test_mucosa_unbound_observer(tmp_path: Path) -> None:
    """`Cu_mu_<seg>_<s> = fu/K_mu*Cmu`, an observer of every mucosa segment."""
    r = _build(tmp_path, _definition("observer", chemistry=[_gut_metabolism()]))
    r[ids.sol_dose(S)] = 10.0
    r.simulate(0, 60, 5)
    for seg in SMALL + COLON:
        c = r[f"[{_gu(ids.mucosa(S, seg))}]"]
        k = r[_gu(ids.mucosa_coefficient(S))]
        assert c > 0
        assert r[_gu(ids.mucosa_unbound(S, seg))] == pytest.approx(
            r[ids.fu(S)] / k * c, rel=1e-12
        )


def test_gut_turnover_steady_state(tmp_path: Path) -> None:
    """The turnover of CYP3A4 in the gut reaches the analytic steady state.

    A constant infusion of an inducer and mechanism-based inhibitor gives the
    same unbound concentration `I` in every mucosa segment;
    `E = (1 + Emax*I/(EC50 + I))/(1 + kinact/kdeg*I/(KI + I))`. Only the
    segments with CYP3A4 (not the colon) drive the state.
    """
    emax, ec50, ki, kinact = 2.0, 1e-3, 2e-3, 1e-3
    inhibitor = Substance(
        sid="inh",
        name="inhibitor",
        mr=400.0,
        excretion=Excretion(urine=True),
        inhibits=(Inhibition("CYP3A4", "mechanism", ki, kinact=kinact),),
        induces=(Induction("CYP3A4", emax, ec50),),
    )
    d = _definition(
        "turnover",
        substances=[_drug(), _metabolite(), inhibitor],
        chemistry=[_gut_metabolism()],
        routes=[Route.iv("inh")],
    )
    result = build(d, tmp_path)
    doc = libsbml.readSBMLFromFile(str(result.flat_path))
    rule = doc.getModel().getRateRule(_gu(ids.enzyme("CYP3A4")))
    drive = libsbml.formulaToL3String(rule.getMath())
    assert _gu(ids.mucosa("inh", "du")) in drive
    assert not any(_gu(ids.mucosa("inh", seg)) in drive for seg in COLON)
    r = load(result)
    r[ids.infusion_rate("inh")] = 0.01
    r.simulate(0, 20 * DAY, 20)
    unbound = [r[_gu(ids.mucosa_unbound("inh", seg))] for seg in SMALL]
    assert unbound == pytest.approx([unbound[0]] * len(SMALL), rel=1e-6)
    i = unbound[0]
    kdeg = r[_gu(ids.kdeg("CYP3A4"))]
    expected = (1 + emax * i / (ec50 + i)) / (1 + kinact / kdeg * i / (ki + i))
    assert r[_gu(ids.enzyme("CYP3A4"))] == pytest.approx(expected, rel=1e-6)


def test_bile_drains_for_every_substance(tmp_path: Path) -> None:
    """The bile of a substance without lumen species drains into the feces."""
    d = _definition("bile_drain", routes=[Route.oral(S)])
    doc = libsbml.readSBMLFromFile(str(build(d, tmp_path).flat_path))
    model = doc.getModel()
    for sid, into in ((S, _gu(ids.lumen(S, "du"))), (M, ids.feces(M))):
        reaction = model.getReaction(_gu(ids.bile_reaction(sid)))
        assert [x.getSpecies() for x in reaction.getListOfReactants()] == [
            ids.bile(sid)
        ]
        assert [x.getSpecies() for x in reaction.getListOfProducts()] == [into]
