"""Tests of the nephron kidney, the variant `nephron` of the kidney.

The nephron kidney filters the unbound plasma concentration into the lumen of
the proximal tubule, `fu*GFR*Cpl`; the tubular fluid flows through proximal
tubule, loop of Henle, distal tubule and collecting duct, each segment
reabsorbs water (the flow leaving a segment is the fraction `f_tub_<seg>` of
the GFR), the drug moves with the fluid and the collecting duct drains into
the urine. Without secretion and reabsorption every filtered molecule reaches
the urine: at steady state the urinary excretion rate is the filtration
rate, the renal clearance `fu*GFR` of the standard kidney, whatever the
volumes and flows of the segments. Apical transport acts between the cells of
the kidney and the lumen of the proximal tubule: an efflux secretes, an
influx reabsorbs.
"""

import dataclasses
from pathlib import Path
from typing import TYPE_CHECKING

import libsbml
import pytest

from pkpdsim.pbpk import (
    MM,
    Binding,
    CustomReaction,
    Direction,
    Erythrocytes,
    Excretion,
    FirstOrder,
    Inhibition,
    Kp,
    Membrane,
    Metabolism,
    Mode,
    NephronKidney,
    Organ,
    Part,
    PBPKDefinition,
    Pool,
    Route,
    StandardOrgan,
    Substance,
    Transport,
    Variants,
    ids,
    variants,
)
from pkpdsim.pbpk.builder import build, organ_context
from pkpdsim.pbpk.organ import OrganComponent

from .conftest import definition, dose_amount, individual, simple_drug, total_amount

if TYPE_CHECKING:
    import roadrunner

S = "drg"
"""The substance of the tests."""

DAY = 24 * 60.0
"""One day (min)."""

STEADY = 2 * DAY
"""The end of the infusions (min): the drugs of the steady state tests are
vascular (and enter only the kidney tissue), with half-lives below an hour, and
the tubule passes its fluid within minutes, so the steady state is reached far
below the tolerances."""

INFUSION = 1.0
"""The infusion rate of the steady state tests (mg/min)."""

SEGMENTS = ("pt", "lh", "dt", "cd")
"""The tubule segments in the order of the flow."""


def _kidney(local: str) -> str:
    """The flat id of a local id of the kidney."""
    return ids.flat(Organ.KIDNEY, local)


def _renal_drug(
    fraction_reabsorbed: float = 0.0, urine: bool = True, **changes: object
) -> Substance:
    """A vascular drug (fu 0.5), filtered, which reaches its steady state within hours."""
    drug = Substance(
        sid=S,
        name="drug",
        mr=300.0,
        binding=Binding(fu=0.5, protein="albumin"),
        excretion=Excretion(urine=urine, fraction_reabsorbed=fraction_reabsorbed),
    )
    return dataclasses.replace(drug, **changes)


_KIDNEY_TISSUE = Kp(
    1.0,
    organs=frozenset({Organ.KIDNEY}),
    ftissue=1.0,
    ps={Organ.KIDNEY: 0.2},
    ps_per_volume=False,
)
"""The drug enters the tissue of the kidney only, with the passive exchange
0.2 l/min of the cells instead of the share of `ftissue` (the transporters
give the kidney the full form)."""


def _secretion(k: float = 0.5, protein: str | None = None) -> Transport:
    """The apical efflux from the kidney cells into the proximal tubule (secretion)."""
    return Transport(
        sid="DRGSEC",
        organ=Organ.KIDNEY,
        substance=S,
        direction=Direction.EFFLUX,
        membrane=Membrane.APICAL,
        kinetics=FirstOrder(k=k, per_volume=False),
        protein=protein,
    )


def _reabsorption(k: float = 0.05, protein: str | None = None) -> Transport:
    """The apical influx from the proximal tubule into the kidney cells (reabsorption)."""
    return Transport(
        sid="DRGREAB",
        organ=Organ.KIDNEY,
        substance=S,
        direction=Direction.INFLUX,
        membrane=Membrane.APICAL,
        kinetics=FirstOrder(k=k, per_volume=False),
        protein=protein,
    )


def _model(
    tmp_path: Path,
    sid: str,
    kidney: str | OrganComponent,
    substances: list[Substance] | None = None,
    chemistry: list[Metabolism | Transport | CustomReaction] | None = None,
) -> Path:
    """Build the renal drug (default `_renal_drug`), given intravenously, with a kidney variant."""
    d = definition(
        sid=sid,
        substances=substances or [_renal_drug()],
        chemistry=chemistry or [],
        routes=[Route.iv(S)],
        variants=Variants(organs={Organ.KIDNEY: kidney}),
    )
    return build(d, tmp_path / sid).flat_path


def _load(path: Path, **changes: float) -> "roadrunner.RoadRunner":
    """Load a model with changes before its initialization, tight tolerances."""
    r = individual(path, **changes)
    r.integrator.absolute_tolerance = 1e-12
    r.integrator.relative_tolerance = 1e-10
    return r


def _urine_rate(r: "roadrunner.RoadRunner") -> float:
    """The rate of the urinary excretion (mmole/min), out of the collecting duct."""
    return r[_kidney(ids.tubule_reaction(S, "cd"))]


def _kidney_plasma(r: "roadrunner.RoadRunner") -> float:
    """The plasma concentration of the kidney (mM), `KI__Cpl_<s>`."""
    return r[f"[{_kidney(ids.species(Part.PLASMA, S))}]"]


def _steady(r: "roadrunner.RoadRunner") -> None:
    """Run the infusion into its steady state."""
    r[ids.infusion_rate(S)] = INFUSION
    r.simulate(0, STEADY, 2)


def test_registry_has_nephron_kidney(tmp_path: Path) -> None:
    """The kidney has the variant `nephron`; it builds by name, with the tubule species."""
    assert variants.names(Organ.KIDNEY) == ["nephron", "standard"]
    assert variants.lookup(Organ.KIDNEY, "nephron") == NephronKidney()
    r = _load(_model(tmp_path, "registry", "nephron"))
    for seg in SEGMENTS:
        assert r[f"[{_kidney(ids.tubule(S, seg))}]"] == 0.0
        assert r[_kidney(ids.tubule_volume(seg))] > 0.0
    flows = [r[_kidney(ids.tubule_flow_fraction(seg))] for seg in SEGMENTS]
    assert flows == sorted(flows, reverse=True)
    assert 0.0 < flows[-1] < flows[0] < 1.0


def _boundary_definitions() -> dict[str, PBPKDefinition]:
    """The definitions of the boundary test: every form of a filtered substance in the kidney."""
    metabolite = Substance(
        sid="drg_m", name="metabolite", mr=316.0, excretion=Excretion(urine=True)
    )
    kinetic = Erythrocytes("kinetic", kbc=2.0, k_exchange=1.0)
    uptake = Transport(
        sid="DRGUP",
        organ=Organ.KIDNEY,
        substance=S,
        direction=Direction.INFLUX,
        membrane=Membrane.BASOLATERAL,
        kinetics=MM(vmax=1e-3, km=1e-3),
        protein="SLC22A2",
    )
    renal_metabolism = Metabolism(
        sid="DRGKI",
        organ=Organ.KIDNEY,
        substrate=S,
        products={"drg_m": 1.0},
        kinetics=FirstOrder(k=0.01),
    )
    return {
        "merged": definition(),
        "reabsorbed": definition(
            substances=[
                dataclasses.replace(
                    simple_drug(), excretion=Excretion(True, fraction_reabsorbed=0.3)
                ),
                metabolite,
            ]
        ),
        "full": definition(
            substances=[
                dataclasses.replace(simple_drug(), erythrocytes=kinetic),
                metabolite,
            ],
            chemistry=[uptake, _secretion(protein="SLC47A1"), renal_metabolism],
        ),
        "perfusion": definition(
            substances=[
                dataclasses.replace(
                    simple_drug(),
                    distribution=Kp(2.0, mode=Mode.PERFUSION),
                    erythrocytes=Erythrocytes("equilibrium", bp=1.3),
                ),
                metabolite,
            ]
        ),
        "vascular": definition(substances=[_renal_drug()], chemistry=[]),
    }


@pytest.mark.parametrize(
    "form", ["merged", "reabsorbed", "full", "perfusion", "vascular"]
)
def test_boundary_same_as_standard(form: str) -> None:
    """The nephron kidney reads and gives the elements of the standard kidney.

    The same inputs replaced by the same elements of the body (organ volume,
    flow, fractions, `HCT`, `GFR`, the substance parameters with
    `f_reabsorbed_<s>`, the pools, the urine) and the same outputs, for every
    form of a filtered substance in the kidney: merged tissue, a fraction
    reabsorbed, the full form with basolateral uptake, secretion, renal
    metabolism and kinetic erythrocytes, the flow-limited kidney with
    erythrocytes in equilibrium, a vascular drug.
    """
    ctx = organ_context(_boundary_definitions()[form], Organ.KIDNEY)
    _, standard = StandardOrgan().build(ctx)
    md, nephron = NephronKidney().build(ctx)
    assert nephron.replacements == standard.replacements
    assert nephron.outputs == standard.outputs
    assert md.sid == "nephron_ki"


def test_boundary_secretion_only_reads_gfr() -> None:
    """A secreted substance which is not filtered: the tubular flow reads the GFR as well.

    The standard kidney reads the GFR only for the filtration; the tubular
    fluid of the nephron kidney is a fraction of the GFR, the only additional
    element of the body.
    """
    d = definition(
        substances=[_renal_drug(urine=False, distribution=_KIDNEY_TISSUE)],
        chemistry=[_secretion()],
    )
    ctx = organ_context(d, Organ.KIDNEY)
    _, standard = StandardOrgan().build(ctx)
    _, nephron = NephronKidney().build(ctx)
    assert dict(nephron.replacements) == {**standard.replacements, ids.GFR: ids.GFR}
    assert nephron.outputs == standard.outputs


@pytest.mark.parametrize("fraction_reabsorbed", [0.0, 0.4])
def test_no_secretion_no_reabsorption_equals_fu_gfr(
    tmp_path: Path, fraction_reabsorbed: float
) -> None:
    """Without secretion and reabsorption the renal clearance is `fu*GFR`, as in the standard kidney.

    At the steady state of an infusion every compartment of the tubule is
    stationary, so the flow out of each segment equals the flow into it: the
    urinary excretion rate is the filtration rate `fu*GFR*Cpl` (exactly,
    whatever the volumes and the water reabsorption of the segments), and it
    is the infusion. A fraction reabsorbed returns from the collecting duct
    into the plasma of the kidney, `fu*GFR*(1 - f_reabsorbed)` as the
    standard kidney. The circulation sees the same kidney: the venous
    steady state equals that of the standard kidney, and the plasma leaving
    the kidney carries the plasma flow less the net filtered amount, so the
    urine rate over the arterial plasma is the clearance in both kidneys.
    """
    drug = _renal_drug(fraction_reabsorbed)
    venous = f"[{ids.pool_species(Pool.VENOUS, S)}]"
    arterial = f"[{ids.pool_species(Pool.ARTERIAL, S)}]"
    result = {}
    for kidney in ("standard", "nephron"):
        r = _load(_model(tmp_path, kidney, kidney, substances=[drug]))
        _steady(r)
        result[kidney] = r[venous]
        clearance = r[ids.fu(S)] * r[ids.GFR] * (1.0 - fraction_reabsorbed)
        rate = r[ids.urine_rate(S)]
        assert rate / r[arterial] == pytest.approx(clearance, rel=1e-6)
    assert _urine_rate(r) == pytest.approx(INFUSION / r[ids.mr(S)], rel=1e-8)
    assert _urine_rate(r) / _kidney_plasma(r) == pytest.approx(clearance, rel=1e-3)
    assert _urine_rate(r) / _kidney_plasma(r) == pytest.approx(clearance, rel=1e-8)
    assert result["nephron"] == pytest.approx(result["standard"], rel=1e-8)
    # the drug concentrates along the tubule as the water is reabsorbed
    lumen = [r[f"[{_kidney(ids.tubule(S, seg))}]"] for seg in SEGMENTS]
    assert lumen == sorted(lumen)
    assert lumen[-1] / lumen[0] == pytest.approx(
        r[_kidney(ids.tubule_flow_fraction("pt"))]
        / r[_kidney(ids.tubule_flow_fraction("cd"))],
        rel=1e-8,
    )


def _effective_secretion(r: "roadrunner.RoadRunner", k: float) -> float:
    """The secretion clearance on the unbound plasma concentration of the kidney (l/min).

    At steady state the unbound drug passes from plasma to interstitium
    (`PS_vas`), into the cells (`PS`) and into the lumen (`k`) with one flux:
    three conductances in series, `1/(1/PS_vas + 1/PS + 1/k)`.
    """
    ps_vas, ps = r[_kidney(ids.PS_VAS)], r[_kidney(ids.ps(S))]
    return 1.0 / (1.0 / ps_vas + 1.0 / ps + 1.0 / k)


def test_secretion_increases_clearance(tmp_path: Path) -> None:
    """Secretion into the proximal tubule adds to the filtration: the analytic well-stirred value.

    At steady state the urine rate is `fu*(GFR + k_sec)*Cpl` of the kidney
    plasma, with `k_sec` the conductances of the passage plasma -
    interstitium - cell - lumen in series. The filtrate leaves the plasma of
    the kidney with the plasma flow, the secreted drug is taken out of the
    plasma flow: on the arterial inflow the kidney is well-stirred for the
    secretion, `Qpl*Car = Qpl*Cpl + CL_sec*Cpl`, so
    `CL_R = Qpl*(CL_gfr + CL_sec)/(Qpl + CL_sec)` with `CL_gfr = fu*GFR`,
    `CL_sec = fu*k_sec` and the plasma flow `Qpl = Qki*(1 - HCT)`.
    """
    k = 0.5
    drug = _renal_drug(distribution=_KIDNEY_TISSUE)
    r = _load(
        _model(tmp_path, "secretion", "nephron", [drug], chemistry=[_secretion(k)])
    )
    _steady(r)
    fu, gfr = r[ids.fu(S)], r[ids.GFR]
    secretion = fu * _effective_secretion(r, k)
    clearance = fu * gfr + secretion
    rate = _urine_rate(r)
    assert rate == pytest.approx(INFUSION / r[ids.mr(S)], rel=1e-8)
    assert rate / _kidney_plasma(r) == pytest.approx(clearance, rel=1e-8)
    assert rate / _kidney_plasma(r) > 2 * fu * gfr
    qpl = r[ids.organ_flow(Organ.KIDNEY)] * (1.0 - r[ids.HCT])
    arterial = r[f"[{ids.pool_species(Pool.ARTERIAL, S)}]"]
    assert rate / arterial == pytest.approx(
        qpl * clearance / (qpl + secretion), rel=1e-8
    )


def test_nephron_kidney_scales_secretion(tmp_path: Path) -> None:
    """The secretion of the nephron kidney follows the renal function, as in the standard kidney.

    With half the renal function (a CKD stage, which changes fu) the
    clearance of the apical efflux on its unbound substrate in the cells
    halves (the tubular capacity, the intact nephron hypothesis).
    """
    drug = _renal_drug(distribution=_KIDNEY_TISSUE)
    path = _model(tmp_path, "renal", "nephron", [drug], chemistry=[_secretion()])

    def clearance(**changes: float) -> float:
        r = _load(path, **changes)
        r[f"[{_kidney(ids.species(Part.INTRACELLULAR, S))}]"] = 1e-3
        unbound = r[ids.fu(S)] * 1e-3 / r[_kidney(ids.k_cell_pls(S))]
        return r[_kidney("DRGSEC")] / unbound

    assert clearance() == pytest.approx(0.5, rel=1e-9)
    assert clearance(f_renal_function=0.5) / clearance() == pytest.approx(0.5, rel=1e-9)


def test_tubular_reabsorption_decreases_clearance(tmp_path: Path) -> None:
    """An apical influx from the proximal tubule reabsorbs a share of the filtrate.

    At steady state the proximal tubule passes on the share
    `f_tub_pt*GFR/(f_tub_pt*GFR + k)` of the filtered drug, the influx with
    the clearance `k` on the lumen concentration takes the rest back into the
    cells; downstream every molecule reaches the urine.
    """
    k = 0.05
    drug = _renal_drug(distribution=_KIDNEY_TISSUE)
    r = _load(_model(tmp_path, "reab", "nephron", [drug], chemistry=[_reabsorption(k)]))
    _steady(r)
    fu, gfr = r[ids.fu(S)], r[ids.GFR]
    flow = r[_kidney(ids.tubule_flow_fraction("pt"))] * gfr
    rate = _urine_rate(r)
    assert rate == pytest.approx(INFUSION / r[ids.mr(S)], rel=1e-8)
    assert rate / _kidney_plasma(r) == pytest.approx(
        fu * gfr * flow / (flow + k), rel=1e-8
    )
    assert r[_kidney("DRGREAB")] > 0


def test_unplaced_process_raises(tmp_path: Path) -> None:
    """A custom reaction of the kidney is refused, naming it; the nephron is a kidney only.

    A custom reaction states its rate for a part of the standard kidney, the
    nephron kidney places metabolism and transport only (R34).
    """
    custom = CustomReaction(
        sid="DRGC",
        organ=Organ.KIDNEY,
        part=Part.TISSUE,
        equation=f"{ids.species(Part.TISSUE, S)} ->",
        formula=f"k_drgc * {ids.species(Part.TISSUE, S)}",
        parameters=[],
    )
    d = definition(chemistry=[custom])
    with pytest.raises(ValueError, match="DRGC"):
        NephronKidney().build(organ_context(d, Organ.KIDNEY))
    d = definition(variants=Variants(organs={Organ.LIVER: NephronKidney()}))
    with pytest.raises(ValueError, match="variant of the kidney"):
        build(d, tmp_path)


def test_standard_kidney_still_refuses_apical_influx() -> None:
    """The apical influx from the urine stays refused in the standard kidney (R18)."""
    d = definition(
        substances=[_renal_drug(distribution=_KIDNEY_TISSUE)],
        chemistry=[_reabsorption()],
    )
    with pytest.raises(ValueError, match=r"DRGREAB.*nephron"):
        StandardOrgan().build(organ_context(d, Organ.KIDNEY))


def test_mass_balance_nephron(tmp_path: Path) -> None:
    """The nephron kidney conserves the amounts.

    The drug with kinetic erythrocytes, filtered with a fraction reabsorbed,
    taken up into the kidney cells (OCT2), secreted (MATE1, saturable),
    reabsorbed from the lumen (OAT3 as a stand-in protein) and metabolized in
    the cells into a filtered metabolite; an inhibitor of the reabsorption,
    filtered, which acts from the lumen of the proximal tubule. The amounts of
    drug and metabolite (body, tubule, urine, bile, intestine, feces) stay the
    dose, and so does the inhibitor.
    """
    drug = dataclasses.replace(
        simple_drug(),
        erythrocytes=Erythrocytes("kinetic", kbc=2.0, k_exchange=1.0),
        excretion=Excretion(urine=True, fraction_reabsorbed=0.2),
    )
    metabolite = Substance(
        sid="drg_m", name="metabolite", mr=316.0, excretion=Excretion(urine=True)
    )
    inhibitor = Substance(
        sid="inh",
        name="inhibitor",
        mr=400.0,
        distribution=Kp(1.0, ftissue=1.0),
        excretion=Excretion(urine=True),
        inhibits=(Inhibition(protein="SLC22A8", kind="competitive", ki=1e-3),),
    )
    chemistry: list[Metabolism | Transport | CustomReaction] = [
        Transport(
            sid="DRGUP",
            organ=Organ.KIDNEY,
            substance=S,
            direction=Direction.INFLUX,
            membrane=Membrane.BASOLATERAL,
            kinetics=MM(vmax=1e-2, km=1e-3),
            protein="SLC22A2",
        ),
        dataclasses.replace(
            _secretion(), kinetics=MM(vmax=1e-2, km=1e-3), protein="SLC47A1"
        ),
        _reabsorption(0.05, protein="SLC22A8"),
        Metabolism(
            sid="DRGKI",
            organ=Organ.KIDNEY,
            substrate=S,
            products={"drg_m": 1.0},
            kinetics=FirstOrder(k=0.01),
        ),
    ]
    d = definition(
        sid="balance",
        substances=[drug, metabolite, inhibitor],
        chemistry=chemistry,
        routes=[Route.iv(S), Route.iv("inh")],
        variants=Variants(organs={Organ.KIDNEY: "nephron"}),
    )
    r = _load(build(d, tmp_path).flat_path)
    r[ids.iv_dose(S)] = 10.0
    r[ids.iv_dose("inh")] = 20.0
    r.simulate(0, DAY, 49)
    drug_amount = total_amount(r, S) + total_amount(r, "drg_m") + dose_amount(r, S)
    assert drug_amount == pytest.approx(10.0 / r[ids.mr(S)], rel=1e-7)
    inh_amount = total_amount(r, "inh") + dose_amount(r, "inh")
    assert inh_amount == pytest.approx(20.0 / r[ids.mr("inh")], rel=1e-7)
    for sid in (S, "drg_m", "inh"):
        assert r[ids.urine(sid)] > 0
    assert r[_kidney(ids.tubule_reabsorption_reaction(S))] > 0
    assert r[_kidney("DRGSEC")] > 0
    assert r[_kidney("DRGREAB")] > 0
    model = libsbml.readSBMLFromString(r.getCurrentSBML()).getModel()
    reabsorption = model.getReaction(_kidney("DRGREAB"))
    modifiers = [m.getSpecies() for m in reabsorption.getListOfModifiers()]
    assert modifiers == [_kidney(ids.tubule("inh", "pt"))]
    tubule_species = ids.species_of(S, [s.getId() for s in model.getListOfSpecies()])
    for seg in SEGMENTS:
        assert _kidney(ids.tubule(S, seg)) in tubule_species
