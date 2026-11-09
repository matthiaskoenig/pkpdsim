"""Tests of the zonated liver, the variant `zonated` of the liver.

The zonated liver is the standard liver in `n` zones in series along the
sinusoid. Its boundary is that of the standard liver; with one zone it is the
standard liver. With `n` well-stirred zones in series and linear kinetics
the steady-state availability of the liver is `(1 + x/n)^-n` with the
dimensionless clearance `x` of the liver, the well-stirred `1/(1 + x)` for
one zone: the zonated liver extracts more than the standard liver with the
same intrinsic clearance, and the difference of the extraction is about
`x/2*(1 - 1/n)` of it. The equivalence with the standard liver is therefore
shown at a low extraction (`x` about 1e-3), where the difference is below
the tolerance 1e-3, and the zones at a moderate extraction against the
closed formula of the tanks in series.
"""

import dataclasses
import itertools
import math
from pathlib import Path
from typing import TYPE_CHECKING

import libsbml
import numpy as np
import pytest
from sbmlutils.factory import Parameter as SbmlParameter

from pkpdsim.pbpk import (
    MM,
    Binding,
    CustomReaction,
    Direction,
    Erythrocytes,
    FirstOrder,
    Inhibition,
    Kp,
    Membrane,
    Metabolism,
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
    ZonatedLiver,
    distribution,
    ids,
    variants,
)
from pkpdsim.pbpk.builder import build, organ_context
from pkpdsim.pbpk.organ import OrganComponent, OrganContext
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Value
from pkpdsim.zonation import ids as zids
from pkpdsim.zonation.profiles import Linear
from pkpdsim.zonation.settings import Geometry
from pkpdsim.zonation.transport import M2_PER_S, diffusion_parameter

from .conftest import (
    definition,
    dose_amount,
    individuals,
    initialized,
    simple_drug,
    total_amount,
)

if TYPE_CHECKING:
    from collections.abc import Callable

    import roadrunner
    from sbmlutils.factory import Sbase

S = "drg"
"""The substance of the tests."""

DAY = 24 * 60.0
"""One day (min)."""

STEADY = 10 * DAY
"""The end of the infusions (min): the drug of the steady state tests (fu 0.5,
Kp 0.5, filtered) has a half-life of a few hours, so its steady state is
reached far below the tolerances."""

INFUSION = 1.0
"""The infusion rate of the steady state tests (mg/min)."""

TIMES = np.linspace(0.0, DAY, 1441)
"""The time points of the curves (min)."""

MIXED = TIMES >= 5.0
"""The time points of the comparison of the same equations, after the mixing of the bolus."""

PASSED = TIMES >= 30.0
"""The time points of the comparison of zones and the well-stirred liver.

The zones have the mean transit time of the well-stirred liver, but a
narrower distribution of the transit times, and the uptake into the liver
tissue while it fills is an extraction of the liver as well (`fu*PS/Qpl`
about 0.1 for the test drug): the curves differ (about 1e-2) as long as the
bolus mixes and the liver tissue fills, in the first 30 min."""


def _liver(local: str) -> str:
    """The flat id of a local id of the liver."""
    return ids.flat(Organ.LIVER, local)


def _zones(n: int) -> list[str]:
    """The zones of a zonated liver of `n` zones."""
    return [ids.zone(i) for i in range(1, n + 1)]


def _fast_drug(**changes: object) -> Substance:
    """A drug which reaches its steady state within hours.

    fu 0.5, Kp 0.5 in the merged tissue with the exchange `ftissue` 10 l/min,
    filtered into the urine (`fu*GFR` about 0.06 l/min).
    """
    drug = dataclasses.replace(
        simple_drug(),
        binding=Binding(fu=0.5, protein="albumin"),
        distribution=Kp(0.5, ftissue=10.0),
    )
    return dataclasses.replace(drug, **changes)


def _linear(clint: float) -> Metabolism:
    """The first order liver metabolism with the unbound intrinsic clearance `clint` (l/min)."""
    return Metabolism(
        sid="DRGX",
        organ=Organ.LIVER,
        substrate=S,
        products={},
        kinetics=FirstOrder(k=clint, per_volume=False),
        protein="CYP3A4",
    )


def _drg2m(**changes: object) -> Metabolism:
    """The Michaelis-Menten CYP3A4 metabolism `DRG2M` of the test helpers."""
    process = Metabolism(
        sid="DRG2M",
        organ=Organ.LIVER,
        substrate=S,
        products={"drg_m": 1.0},
        kinetics=MM(vmax=1e-3, km=0.01),
        protein="CYP3A4",
    )
    return dataclasses.replace(process, **changes)


def _model(
    tmp_path: Path,
    sid: str,
    liver: str | OrganComponent,
    chemistry: list[Metabolism | Transport | CustomReaction],
    substances: list[Substance] | None = None,
) -> Path:
    """Build the drug (default `_fast_drug`), given intravenously, with a liver variant."""
    d = definition(
        sid=sid,
        substances=substances or [_fast_drug()],
        chemistry=chemistry,
        routes=[Route.iv(S)],
        variants=Variants(organs={Organ.LIVER: liver}),
    )
    return build(d, tmp_path / sid).flat_path


def _load(
    model: "Path | Callable[..., roadrunner.RoadRunner]", **changes: float
) -> "roadrunner.RoadRunner":
    """Load a model with changes before its initialization, tight tolerances.

    The model is a file or the `individuals` of one, loaded once for a test.
    """
    r = initialized(model, **changes)
    r.integrator.absolute_tolerance = 1e-12
    r.integrator.relative_tolerance = 1e-10
    return r


def _curves(r: "roadrunner.RoadRunner") -> np.ndarray:
    """The venous plasma and the liver outflow concentrations (mM) after an IV bolus of 10 mg."""
    r[ids.iv_dose(S)] = 10.0
    cols = [f"[{ids.pool_species(Pool.VENOUS, S)}]", ids.organ_outflow(Organ.LIVER, S)]
    return r.simulate(times=TIMES, selections=["time", *cols])[:, 1:]


def _deviation(a: np.ndarray, b: np.ndarray, times: np.ndarray = MIXED) -> float:
    """The largest relative deviation of the curves `b` from the curves `a` at the times."""
    return float(np.max(np.abs(b[times] - a[times]) / a[times]))


def _extraction(r: "roadrunner.RoadRunner") -> float:
    """The extraction of the liver at the steady state of the infusion.

    `1 - Cout/Cin` with the plasma concentration of the inflow, the mix of
    the hepatic artery and the portal vein `(Qha*Car + Qpo*Cpo)/Qh`; the
    steady state is checked (urine and metabolism eliminate the infusion).
    """
    r[ids.infusion_rate(S)] = INFUSION
    r.simulate(0, STEADY, 2)
    metabolism = sum(
        r[i] for i in r.model.getReactionIds() if i.startswith(_liver("DRGX"))
    )
    urine = r[ids.flat(Organ.KIDNEY, ids.filtration_reaction(S))]
    assert metabolism + urine == pytest.approx(INFUSION / r[ids.mr(S)], rel=1e-9)
    car = r[f"[{ids.pool_species(Pool.ARTERIAL, S)}]"]
    cpo = r[f"[{ids.pool_species(Pool.PORTAL, S)}]"]
    cin = (r[ids.QHA] * car + r[ids.QPO] * cpo) / r[ids.QH]
    return 1.0 - r[ids.organ_outflow(Organ.LIVER, S)] / cin


def _x(r: "roadrunner.RoadRunner", clint: float) -> float:
    """The dimensionless clearance `x` of the liver of the merged tissue.

    A zone exchanges `PS/n*fu*(Cpl - Cti/Kp)` and metabolizes
    `CLint/n*fu*Cti/Kp`; at steady state its unbound tissue concentration is
    `fu*Cpl*PS/(PS + CLint)`, so it clears the plasma with
    `fu*PS*CLint/(PS + CLint)/n`. The plasma flow `Qh*(1 - HCT)` carries it,
    `x = fu*PS*CLint/((PS + CLint)*Qh*(1 - HCT))`.
    """
    ps, fu = r[_liver(ids.ps(S))], r[ids.fu(S)]
    return fu * ps * clint / ((ps + clint) * r[ids.QH] * (1.0 - r[ids.HCT]))


def test_registry_has_zonated_liver(tmp_path: Path) -> None:
    """The liver has the variant `zonated`, 10 zones placed by volume."""
    assert variants.names(Organ.LIVER) == ["standard", "zonated"]
    liver = variants.lookup(Organ.LIVER, "zonated")
    assert liver == ZonatedLiver(n_zones=10)
    path = _model(tmp_path, "registry", "zonated", [_linear(1.0)])
    r = _load(path)
    assert r[f"[{_liver(ids.species(Part.TISSUE, S, ids.zone(10)))}]"] == 0.0
    assert r[_liver(ids.F_ZONE)] == pytest.approx(1 / 10, rel=1e-15)


_UPTAKE = Transport(
    sid="DRGUP",
    organ=Organ.LIVER,
    substance=S,
    direction=Direction.INFLUX,
    membrane=Membrane.BASOLATERAL,
    kinetics=MM(vmax=1e-3, km=1e-3),
    protein="SLCO1B1",
)
"""A saturable uptake into the liver cells (full form)."""

_BILE = Transport(
    sid="DRGBILE",
    organ=Organ.LIVER,
    substance=S,
    direction=Direction.EFFLUX,
    membrane=Membrane.APICAL,
    kinetics=FirstOrder(k=0.5),
)
"""The efflux from the liver cells into the bile."""


def _boundary_definitions() -> dict[str, PBPKDefinition]:
    """The definitions of the boundary test: every form of a substance in the liver."""
    metabolite = Substance(sid="drg_m", name="metabolite", mr=316.0)
    kinetic = Erythrocytes("kinetic", kbc=2.0, k_exchange=1.0)
    return {
        "merged": definition(),
        "full": definition(
            substances=[
                dataclasses.replace(simple_drug(), erythrocytes=kinetic),
                metabolite,
            ],
            chemistry=[_drg2m(), _UPTAKE, _BILE],
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
        "mbi": definition(
            substances=[
                simple_drug(),
                metabolite,
                Substance(
                    sid="inh",
                    name="inhibitor",
                    mr=400.0,
                    distribution=Kp(1.0, ftissue=1.0),
                    inhibits=(
                        Inhibition(
                            protein="CYP3A4", kind="mechanism", ki=1e-3, kinact=0.01
                        ),
                    ),
                ),
            ]
        ),
    }


@pytest.mark.parametrize("form", ["merged", "full", "perfusion", "mbi"])
def test_boundary_same_as_standard(form: str) -> None:
    """The zonated liver reads and gives the elements of the standard liver.

    The same inputs replaced by the same elements of the body, and the same
    outputs (the plasma and the kinetic blood cells leaving the liver), for
    every form of a substance in the liver: merged tissue, the full form with
    transport into the bile and kinetic erythrocytes, the flow-limited liver
    with erythrocytes in equilibrium, a mechanism-based inhibitor.
    """
    ctx = organ_context(_boundary_definitions()[form], Organ.LIVER)
    _, standard = StandardOrgan().build(ctx)
    md, zonated = ZonatedLiver(n_zones=3, profiles={"CYP3A4": [1, 2, 3]}).build(ctx)
    assert zonated.replacements == standard.replacements
    assert zonated.outputs == standard.outputs
    assert md.sid == "zonated_li"


@pytest.mark.parametrize(
    ("form", "changes", "chemistry"),
    [
        ("merged", {}, [_drg2m()]),
        (
            "full",
            {"erythrocytes": Erythrocytes("kinetic", kbc=2.0, k_exchange=1.0)},
            [_drg2m(), _UPTAKE, _BILE],
        ),
        (
            "perfusion",
            {
                "distribution": Kp(2.0, mode=Mode.PERFUSION),
                "erythrocytes": Erythrocytes("equilibrium", bp=1.3),
            },
            [_drg2m()],
        ),
    ],
)
def test_one_zone_is_standard_liver(
    tmp_path: Path,
    form: str,
    changes: dict[str, object],
    chemistry: list[Metabolism | Transport | CustomReaction],
) -> None:
    """A zonated liver of one zone is the standard liver, for every kinetics.

    The same equations with the ids of the zone: the venous and the hepatic
    outflow curves after an IV bolus agree to the integrator (Michaelis-Menten
    metabolism, saturable uptake, bile, shunt and functional fraction).
    """
    drug = dataclasses.replace(simple_drug(), **changes)
    metabolite = Substance(sid="drg_m", name="metabolite", mr=316.0)
    curves = {}
    for name, liver in (("standard", "standard"), ("one", ZonatedLiver(n_zones=1))):
        path = _model(tmp_path, f"{form}_{name}", liver, chemistry, [drug, metabolite])
        curves[name] = _curves(_load(path, f_shunt=0.2, f_functional=0.7))
    assert _deviation(curves["standard"], curves["one"]) < 1e-8


@pytest.mark.parametrize(
    "kinetics",
    [FirstOrder(k=2e-3, per_volume=False), MM(vmax=1.2e-3, km=1.0)],
    ids=["linear", "mm"],
)
def test_uniform_zones_equal_standard_liver(
    tmp_path: Path, kinetics: FirstOrder | MM
) -> None:
    """Uniform zones give the standard liver at a low extraction.

    Zones in series are not a well-stirred liver: with the same intrinsic
    clearance they extract `1 - (1 + x/n)^-n` instead of `x/(1 + x)`, a
    relative difference of about `x/2*(1 - 1/n)` (`x` of the steady state of
    the standard liver, `E/(1 - E)`). The test drug has
    `x` about 1e-3 (asserted), so that the zones of the default variant differ
    from the standard liver by about 7e-4 in the steady-state extraction,
    below the tolerance 1e-3,
    and less in the venous and the hepatic outflow curves after an IV bolus
    (`PASSED`), where the liver is one of two eliminating organs. The Michaelis-Menten metabolism
    (`Vmax*V/Km` the linear clearance, `Km` 1 mM about 20 times the unbound
    concentration of the infusion) is nearly linear here; with a saturated
    metabolism the zones see different concentrations and differ more.
    """
    n = ZonatedLiver().n_zones
    process = dataclasses.replace(_linear(2e-3), kinetics=kinetics)
    results = {}
    for name, liver in (("standard", "standard"), ("zonated", "zonated")):
        path = _model(tmp_path, name, liver, [process])
        r = _load(path)
        x = _x(r, 2e-3 if isinstance(kinetics, FirstOrder) else 1.2e-3 * r["Vli"])
        curves = _curves(r)
        results[name] = (_extraction(_load(path)), curves, x)
    e_standard, c_standard, x = results["standard"]
    e_zonated, c_zonated, _ = results["zonated"]
    # the dimensionless clearance of the well-stirred liver at the steady state
    x_ss = e_standard / (1 - e_standard)
    assert 1e-4 < x_ss < 2e-3
    if isinstance(kinetics, FirstOrder):
        assert x_ss == pytest.approx(x, rel=1e-6)
    else:
        # saturated by `Cu/Km` (about 5 %) at the infusion
        assert x_ss == pytest.approx(x, rel=0.1)
    expected = (1 - (1 + x_ss / n) ** -n) / e_standard - 1
    assert expected == pytest.approx(x_ss / 2 * (1 - 1 / n), rel=1e-2)
    # the saturation `Cu/(Km + Cu)` (about 5 %) lowers the clearance of the
    # zones before the last one, whose concentration is higher, by about that
    # fraction of the difference
    linear = isinstance(kinetics, FirstOrder)
    assert e_zonated / e_standard - 1 == pytest.approx(
        expected, rel=1e-2 if linear else 0.1
    )
    assert abs(e_zonated / e_standard - 1) < 1e-3
    assert _deviation(c_standard, c_zonated, PASSED) < 1e-3


def test_extraction_tanks_in_series(tmp_path: Path) -> None:
    """At a moderate extraction the zones are well-stirred tanks in series.

    A zone of the merged tissue with linear metabolism clears the plasma with
    `fu*PS*CLint/(PS + CLint)/n` (`_x`), so the steady-state availability of
    `n` zones is `(1 + x/n)^-n`: the well-stirred liver `1/(1 + x)` for one
    zone, more extraction with more zones, below the parallel tube `exp(-x)`
    (`x` about 1, an extraction of 0.5; a fast exchange of the tissue,
    `ftissue` 100 l/min). The plasma concentration falls along
    the sinusoid, the last zone is the outflow.
    """
    clint = 30.0
    drug = _fast_drug(distribution=Kp(0.5, ftissue=100.0))
    path = _model(tmp_path, "standard", "standard", [_linear(clint)], [drug])
    r = _load(path)
    x = _x(r, clint)
    assert 0.5 < x < 2.0
    assert _extraction(r) == pytest.approx(x / (1 + x), rel=1e-8)
    extractions = []
    for n in (1, 2, 8):
        r = _load(
            _model(
                tmp_path, f"zones{n}", ZonatedLiver(n_zones=n), [_linear(clint)], [drug]
            )
        )
        e = _extraction(r)
        assert e == pytest.approx(1 - (1 + x / n) ** -n, rel=1e-8)
        extractions.append(e)
        plasma = [r[f"[{_liver(ids.species(Part.PLASMA, S, z))}]"] for z in _zones(n)]
        assert all(a > b for a, b in itertools.pairwise(plasma))
        assert plasma[-1] == pytest.approx(
            r[ids.organ_outflow(Organ.LIVER, S)], rel=1e-12
        )
    assert extractions == sorted(extractions)
    assert extractions[-1] < 1 - math.exp(-x)


OPEN_PERMEABILITY = 1e3
"""The permeability (m/s) of an open fenestral path.

The value of the open barrier of `tests/zonation/test_barrier.py`
(`endothelial_ps=1e3` with the porosity 1), far above the permeability of a
molecule, so that the barrier does not limit the transfer."""

TIGHT_POROSITY = 1e-6
"""The fenestral porosity of a nearly closed barrier (`tests/zonation/test_barrier.py`)."""

NO_DIFFUSION = 1e-30
"""The aqueous diffusion (m**2/s) of a substance which does not diffuse axially.

A `Substance` needs a positive diffusion coefficient
(`tests/zonation/test_transport.py::test_a_non_positive_diffusion_raises`), so
the axial diffusion is switched off with a negligible value instead of 0: the
conductance `D*A/dx` of a zone is then about 1e-21 l/min against a plasma flow
of about 1 l/min."""


def _permeable_drug(diffusion: float | None = None) -> Substance:
    """The drug in the full form (interstitium and cells), fu 0.5, Kp 2.

    The sinusoid geometry has the space of Disse and its barrier only for a
    substance with an interstitial part. `diffusion` is the aqueous diffusion
    coefficient (m**2/s) which replaces the one derived from the molar mass;
    `NO_DIFFUSION` leaves the axial diffusion out.
    """
    drug = dataclasses.replace(
        simple_drug(),
        binding=Binding(fu=0.5, protein="albumin"),
        distribution=Kp(2.0, ftissue=10.0, mode=Mode.PERMEABILITY),
    )
    if diffusion is None:
        return drug
    return dataclasses.replace(drug, diffusion=Value(diffusion, M2_PER_S))


def _barrier_changes(barrier: str, n: int) -> dict[str, float]:
    """The changes of the fenestrated barrier of the sinusoid geometry.

    `open` is the mechanism of `tests/zonation/test_barrier.py`: the fenestral
    porosity of every zone 1 and a permeability far above the one of the
    molecule. `tight` is its nearly closed barrier, the porosity 1e-6.
    """
    if barrier == "tight":
        return {_liver(zids.porosity(z)): TIGHT_POROSITY for z in _zones(n)}
    if barrier != "open":
        raise ValueError(f"The barrier is 'open' or 'tight', not '{barrier}'.")
    changes = {_liver(zids.permeability("fen", S)): OPEN_PERMEABILITY}
    return changes | {_liver(zids.porosity(z)): 1.0 for z in _zones(n)}


def _simulate(
    tmp_path: Path,
    liver: str | OrganComponent,
    diffusion: float | None = None,
    barrier: str | None = None,
) -> dict[str, np.ndarray]:
    """The curves of the drug in the full form after an IV bolus of 10 mg.

    Args:
        tmp_path: the directory of the build (one model per call).
        liver: the liver variant.
        diffusion: the aqueous diffusion of the drug (m**2/s), `None` the one
            derived from its molar mass, `NO_DIFFUSION` no axial diffusion.
        barrier: `open` or `tight` for the fenestrated barrier of the
            sinusoid geometry (`_barrier_changes`), `None` for the measured
            porosity of the table.

    Returns:
        The venous plasma and the hepatic outflow concentration (mM) by their
        selection.
    """
    sid = f"sim{sum(1 for _ in tmp_path.iterdir())}"
    metabolite = Substance(sid="drg_m", name="metabolite", mr=316.0)
    path = _model(
        tmp_path, sid, liver, [_drg2m()], [_permeable_drug(diffusion), metabolite]
    )
    zones = getattr(liver, "n_zones", 1)
    changes = _barrier_changes(barrier, zones) if barrier is not None else {}
    r = _load(path, **changes)
    r[ids.iv_dose(S)] = 10.0
    cols = [f"[{ids.pool_species(Pool.VENOUS, S)}]", ids.organ_outflow(Organ.LIVER, S)]
    result = r.simulate(times=TIMES, selections=["time", *cols])
    return {col: np.asarray(result[:, i + 1]) for i, col in enumerate(cols)}


def test_the_zone_tokens_are_those_of_the_engine() -> None:
    """The zone of the PBPK ids and of the engine are one token, which the adapter needs.

    The zonated liver numbers its zones with the engine (`pkpdsim.zonation`)
    and feeds the tokens into the id functions of the library, so the two id
    modules must agree on the token of a zone.
    """
    assert [ids.zone(k) for k in range(1, 11)] == [zids.zone(k) for k in range(1, 11)]


def _sinusoid_elements(liver: ZonatedLiver, ctx: OrganContext) -> dict[str, "Sbase"]:
    """The parameters and reactions of a liver by id."""
    md, _ = liver.build(ctx)
    elements = [*md.parameters, *md.reactions, *md.compartments, *md.species]
    return {str(e.sid): e for e in elements}


def test_the_sinusoid_geometry_adds_the_elements_of_the_barrier_and_the_diffusion() -> (
    None
):
    """The sinusoid replaces `PS_vas` by the fenestral path and adds the axial diffusion.

    The geometry of the engine adds the radius and the length of the
    sinusoid, the fenestral porosity and the lateral area of every zone, the
    fenestral permeability of the substance in the full form and the
    diffusion coefficient of every substance which diffuses, and one
    reaction per pair of zones and diffusing part; the generic endothelium of
    the standard organ (`PS_vas`, `k_vas`) is gone. Every element it adds has
    an SBO term and an ontology term, as `tests/pbpk/test_annotation.py`
    requires of every element of a model.
    """
    n = 3
    metabolite = Substance(sid="drg_m", name="metabolite", mr=316.0)
    ctx = organ_context(
        definition(substances=[_permeable_drug(), metabolite], chemistry=[_drg2m()]),
        Organ.LIVER,
    )
    stirred = _sinusoid_elements(ZonatedLiver(n_zones=n), ctx)
    sinusoid = _sinusoid_elements(
        ZonatedLiver(n_zones=n, geometry=Geometry.SINUSOID), ctx
    )
    zones = _zones(n)
    diffusions = [
        str(diffusion_parameter(s).sid) for s in (_permeable_drug(), metabolite)
    ]
    expected = {
        zids.RADIUS,
        zids.LENGTH,
        zids.permeability("fen", S),
        *diffusions,
        *(zids.porosity(z) for z in zones),
        *(zids.area("pd", z) for z in zones),
        *(
            ids.axial_reaction(part, s, z)
            for part, s in (
                (Part.PLASMA, S),
                (Part.INTERSTITIAL, S),
                (Part.PLASMA, metabolite.sid),
            )
            for z in zones[:-1]
        ),
    }
    assert set(sinusoid) - set(stirred) == expected
    assert set(stirred) - set(sinusoid) == {ids.PS_VAS, distribution.K_VAS}
    missing = [
        sid
        for sid in sorted(expected)
        if sinusoid[sid].sboTerm is None or not sinusoid[sid].annotations
    ]
    assert not missing


def test_the_sinusoid_with_diffusion_off_and_an_open_barrier_equals_well_stirred(
    tmp_path: Path,
) -> None:
    """The sinusoid without axial diffusion and with an open barrier is the well-stirred liver.

    The well-stirred zone exchanges with the interstitium across the generic
    endothelium `PS_vas = k_vas*Vli`, which is chosen large against the blood
    flows; the sinusoid exchanges through the fenestrae of its zone, the
    geometry of the liver. The two agree when the fenestral path does not
    limit the transfer either (`barrier="open"`) and the sinusoid has no
    axial diffusion (`NO_DIFFUSION`), and the test is not vacuous: a nearly
    closed barrier keeps the drug out of the cells and differs by far more
    than the tolerance.
    """
    n = 10
    stirred = _simulate(
        tmp_path, ZonatedLiver(n_zones=n, geometry=Geometry.WELL_STIRRED)
    )
    opened = _simulate(
        tmp_path,
        ZonatedLiver(n_zones=n, geometry=Geometry.SINUSOID),
        diffusion=NO_DIFFUSION,
        barrier="open",
    )
    tight = _simulate(
        tmp_path,
        ZonatedLiver(n_zones=n, geometry=Geometry.SINUSOID),
        diffusion=NO_DIFFUSION,
        barrier="tight",
    )
    key = f"[{ids.pool_species(Pool.VENOUS, S)}]"
    assert opened[key] == pytest.approx(stirred[key], rel=1e-3, abs=0)
    reference = stirred[key][MIXED]
    assert np.max(np.abs(tight[key][MIXED] - reference) / reference) > 1e-2


def test_profile_preserves_total_capacity(tmp_path: Path) -> None:
    """A profile redistributes the capacity of a process over the zones and keeps the total.

    With the same tissue concentration in every zone the zones together
    metabolize what the standard liver does, with and without the functional
    fraction; the rate of a zone follows the profile (pericentral CYP3A4,
    1 to `n`); a metabolism takes the profile of its name.
    """
    n = ZonatedLiver().n_zones
    profile = [float(i) for i in range(1, n + 1)]
    d = definition(routes=[Route.iv(S)])
    named = definition(
        routes=[Route.iv(S)],
        chemistry=[_drg2m(profile="pericentral")],
    )
    paths = {
        "standard": build(d, tmp_path / "standard").flat_path,
        "protein": build(
            dataclasses.replace(
                d,
                variants=Variants(
                    organs={Organ.LIVER: ZonatedLiver(profiles={"CYP3A4": profile})}
                ),
            ),
            tmp_path / "protein",
        ).flat_path,
        "named": build(
            dataclasses.replace(
                named,
                variants=Variants(
                    organs={
                        Organ.LIVER: ZonatedLiver(profiles={"pericentral": profile})
                    }
                ),
            ),
            tmp_path / "named",
        ).flat_path,
    }
    people = {key: individuals(path) for key, path in paths.items()}
    c = 0.005
    for f_functional in (1.0, 0.5):
        r = _load(people["standard"], f_functional=f_functional)
        r[f"[{_liver(ids.species(Part.TISSUE, S))}]"] = c
        standard = r[_liver("DRG2M")]
        for key, token in (("protein", "cyp3a4"), ("named", "pericentral")):
            r = _load(people[key], f_functional=f_functional)
            for z in _zones(n):
                r[f"[{_liver(ids.species(Part.TISSUE, S, z))}]"] = c
            rates = [r[_liver(ids.segment_process("DRG2M", z))] for z in _zones(n)]
            weights = [r[_liver(ids.placement_weight(token, z))] for z in _zones(n)]
            assert sum(weights) == pytest.approx(1.0, rel=1e-14)
            assert sum(rates) == pytest.approx(standard, rel=1e-12)
            np.testing.assert_allclose(
                np.array(rates) / rates[0], profile, rtol=1e-12, atol=0
            )


def test_a_profile_object_sets_the_zonal_levels(tmp_path: Path) -> None:
    """`profiles` takes a `Profile` of the engine, evaluated at the zone count.

    The levels of the zones are those of the shape, normalized to the mean 1
    by the engine, so the share of a zone is the level over the number of
    zones and the shares sum to 1 (the capacity of the liver is kept).
    """
    n = 4
    shape = Linear(1.0, 4.0)
    d = definition(
        routes=[Route.iv(S)],
        variants=Variants(
            organs={Organ.LIVER: ZonatedLiver(n_zones=n, profiles={"CYP3A4": shape})}
        ),
    )
    r = _load(build(d, tmp_path / "shape").flat_path)
    levels = [r[_liver(ids.placement_profile("cyp3a4", z))] for z in _zones(n)]
    np.testing.assert_allclose(levels, shape.values(n), rtol=1e-12, atol=0)
    weights = [r[_liver(ids.placement_weight("cyp3a4", z))] for z in _zones(n)]
    assert sum(weights) == pytest.approx(1.0, rel=1e-14)
    np.testing.assert_allclose(
        np.array(weights) * n, shape.values(n), rtol=1e-12, atol=0
    )


def test_first_order_per_volume_in_zones(tmp_path: Path) -> None:
    """A first order clearance per litre of liver is placed by zone volume and follows the liver.

    Without a profile a zone has the share `Vti_z/sum(Vti)` of the capacity
    `k*Vli` of the standard liver: with the same tissue concentration in every
    zone the zones together clear `k*Vli*fu/Kp*C`, what the standard liver
    clears, and the clearance scales with the liver of the individual (BW).
    """
    d = definition(
        routes=[Route.iv(S)],
        chemistry=[Metabolism("DRGX", Organ.LIVER, S, {}, FirstOrder(k=0.5))],
    )
    zonated = Variants(organs={Organ.LIVER: ZonatedLiver(n_zones=4)})
    paths = {
        "standard": build(d, tmp_path / "standard").flat_path,
        "zonated": build(
            dataclasses.replace(d, variants=zonated), tmp_path / "zonated"
        ).flat_path,
    }
    people = {key: individuals(path) for key, path in paths.items()}
    c = 0.005
    totals = {}
    for bw in (75.0, 150.0):
        r = _load(people["standard"], BW=bw)
        r[f"[{_liver(ids.species(Part.TISSUE, S))}]"] = c
        clearance = 0.5 * r["Vli"] * r[ids.fu(S)] / r[ids.kp(S)]
        assert r[_liver("DRGX")] == pytest.approx(clearance * c, rel=1e-12)
        r = _load(people["zonated"], BW=bw)
        for z in _zones(4):
            r[f"[{_liver(ids.species(Part.TISSUE, S, z))}]"] = c
        rates = np.array([r[_liver(ids.segment_process("DRGX", z))] for z in _zones(4)])
        volumes = np.array(
            [r[_liver(ids.compartment(Part.TISSUE, z))] for z in _zones(4)]
        )
        assert rates.sum() == pytest.approx(clearance * c, rel=1e-12)
        np.testing.assert_allclose(
            rates / rates.sum(), volumes / volumes.sum(), rtol=1e-12
        )
        totals[bw] = (rates.sum(), r["Vli"])
    assert totals[150.0][1] / totals[75.0][1] > 1.1
    assert totals[150.0][0] / totals[75.0][0] == pytest.approx(
        totals[150.0][1] / totals[75.0][1], rel=1e-12
    )


def _ctx_with_custom_reaction() -> OrganContext:
    """The context of the liver with a custom reaction from the tissue into the plasma."""
    custom = CustomReaction(
        sid="CUSTOM",
        organ=Organ.LIVER,
        part=Part.TISSUE,
        equation=(
            f"{ids.species(Part.TISSUE, S)} -> {ids.species(Part.PLASMA, 'drg_m')}"
        ),
        formula=(
            f"k_custom * {ids.compartment(Part.TISSUE)} * {ids.species(Part.TISSUE, S)}"
        ),
        parameters=[SbmlParameter("k_custom", 0.1, U.per_min)],
    )
    return organ_context(definition(chemistry=[_drg2m(), custom]), Organ.LIVER)


def test_a_custom_reaction_is_placed_per_zone() -> None:
    """A custom reaction is placed in every zone, reading the species of its zone.

    Its rate law is the modeller's and is not split over the zones: in a zone
    it reads the species and the compartments of that zone, so a rate which
    is extensive in the volume of the part (`k*Vti_<zone>*Cti_<zone>_<s>`)
    keeps the rate of the whole liver. The constant of the rate law is
    defined once for the liver.
    """
    model, _ = ZonatedLiver(n_zones=10).build(_ctx_with_custom_reaction())
    assert sum(1 for r in model.reactions if str(r.sid).startswith("CUSTOM_")) == 10
    by_sid = {str(r.sid): r for r in model.reactions}
    assert by_sid["CUSTOM_z1"].equation.to_string() == "Cti_z1_drg => Cpl_z1_drg_m"
    law = by_sid["CUSTOM_z10"].formula
    assert law is not None
    assert str(law.math) == "k_custom * Vti_z10 * Cti_z10_drg"
    assert by_sid["CUSTOM_z1"].compartment == "Vti_z1"
    pars = [str(q.sid) for r in model.reactions for q in r.pars or []]
    assert pars.count("k_custom") == 1


def test_a_custom_reaction_runs_in_every_zone(tmp_path: Path) -> None:
    """The zones together run a volume-extensive custom reaction as the standard liver.

    `k*Vti*Cti_<s>` in the standard liver and `k*Vti_<zone>*Cti_<zone>_<s>`
    in every zone: with the same tissue concentration in every zone the rates
    of the zones sum to the rate of the standard liver, and the metabolite
    appears in the plasma of the zone.
    """
    n = 4
    custom = CustomReaction(
        sid="CUSTOM",
        organ=Organ.LIVER,
        part=Part.TISSUE,
        equation=(
            f"{ids.species(Part.TISSUE, S)} -> {ids.species(Part.PLASMA, 'drg_m')}"
        ),
        formula=(
            f"k_custom * {ids.compartment(Part.TISSUE)} * {ids.species(Part.TISSUE, S)}"
        ),
        parameters=[SbmlParameter("k_custom", 0.1, U.per_min)],
    )
    d = definition(chemistry=[_drg2m(), custom], routes=[Route.iv(S)])
    c = 0.005
    r = _load(build(d, tmp_path / "standard").flat_path)
    r[f"[{_liver(ids.species(Part.TISSUE, S))}]"] = c
    standard = r[_liver("CUSTOM")]
    assert standard > 0
    zonated = dataclasses.replace(
        d, variants=Variants(organs={Organ.LIVER: ZonatedLiver(n_zones=n)})
    )
    r = _load(build(zonated, tmp_path / "zonated").flat_path)
    for z in _zones(n):
        r[f"[{_liver(ids.species(Part.TISSUE, S, z))}]"] = c
    rates = [r[_liver(ids.segment_process("CUSTOM", z))] for z in _zones(n)]
    assert sum(rates) == pytest.approx(standard, rel=1e-12)
    np.testing.assert_allclose(rates, [standard / n] * n, rtol=1e-12, atol=0)


def test_unplaced_process_raises() -> None:
    """A process the zonated liver cannot place is refused, naming it.

    A metabolism with a profile the liver does not give.
    """
    named = definition(chemistry=[_drg2m(profile="pericentral")])
    with pytest.raises(ValueError, match=r"DRG2M.*pericentral"):
        ZonatedLiver().build(organ_context(named, Organ.LIVER))


def test_profiles_used_once_per_protein() -> None:
    """A profile is used by a process, and a protein has one profile.

    A profile which no process uses (a misspelt protein) is refused, and so
    are the processes of one protein with different profiles.
    """
    ctx = organ_context(definition(), Organ.LIVER)
    with pytest.raises(ValueError, match="CYP3a4"):
        ZonatedLiver(n_zones=2, profiles={"CYP3a4": [1, 2]}).build(ctx)
    second = _drg2m(sid="DRG2M2", profile="pericentral")
    ctx = organ_context(definition(chemistry=[_drg2m(), second]), Organ.LIVER)
    with pytest.raises(ValueError, match="CYP3A4"):
        ZonatedLiver(n_zones=2, profiles={"pericentral": [1, 2]}).build(ctx)


@pytest.mark.parametrize(
    ("n_zones", "profiles"),
    [
        (0, {}),
        (2, {"CYP3A4": [1.0]}),
        (2, {"CYP3A4": [1.0, -1.0]}),
        (2, {"CYP3A4": [0.0, 0.0]}),
        (2, {"CYP3A4": [1.0, math.nan]}),
        (2, {"CYP 3A4": [1.0, 2.0]}),
        (2, {"CYP3A4": [1.0, 2.0], "cyp3a4": [2.0, 1.0]}),
    ],
)
def test_invalid_zonated_liver_raises(
    n_zones: int, profiles: dict[str, list[float]]
) -> None:
    """At least one zone, a finite value >= 0 per zone (not all 0), keys which are ids."""
    with pytest.raises(ValueError, match=r"zone|Profile"):
        ZonatedLiver(n_zones=n_zones, profiles=profiles)


def test_zonated_liver_only_for_the_liver(tmp_path: Path) -> None:
    """The zonated liver is no variant of another organ."""
    d = definition(variants=Variants(organs={Organ.KIDNEY: ZonatedLiver()}))
    with pytest.raises(ValueError, match="variant of the liver"):
        build(d, tmp_path)


_BALANCE: dict[str, tuple[dict[str, object], list[Metabolism | Transport], int]] = {
    "full_kinetic": (
        {"erythrocytes": Erythrocytes("kinetic", kbc=2.0, k_exchange=1.0)},
        [_drg2m(), _UPTAKE, _BILE],
        4,
    ),
    "perfusion_bp": (
        {
            "distribution": Kp(2.0, mode=Mode.PERFUSION),
            "erythrocytes": Erythrocytes("equilibrium", bp=1.3),
        },
        [_drg2m()],
        3,
    ),
    "merged_bp": (
        {"erythrocytes": Erythrocytes("equilibrium", bp=1.3)},
        [_drg2m()],
        3,
    ),
}
"""The forms of the drug in the mass balance test: changes of the drug, the
processes of the liver and the number of zones."""

_TISSUE = {
    "full_kinetic": (Part.INTERSTITIAL, Part.INTRACELLULAR),
    "perfusion_bp": (Part.ORGAN,),
    "merged_bp": (Part.TISSUE,),
}
"""The tissue parts of the drug in the liver by form."""


@pytest.mark.parametrize("form", list(_BALANCE))
def test_mass_balance_zonated(tmp_path: Path, form: str) -> None:
    """The zonated liver conserves the amounts.

    The drug with a pericentral CYP3A4 metabolism into the metabolite and a
    mechanism-based inhibitor of CYP3A4, with shunt and functional fraction,
    in the forms: the full form (saturable uptake, bile efflux) with kinetic
    erythrocytes, the flow-limited liver and the merged tissue with
    erythrocytes in equilibrium (the flows between the zones carry
    `Q*BP*C`). The amounts of drug and metabolite (body, bile, intestine,
    feces, urine) stay the dose, and so does the inhibitor. The unbound
    concentrations of every zone are observers, the one turnover state of
    CYP3A4 falls, driven by the zones weighted by the profile.
    """
    changes, chemistry, n = _BALANCE[form]
    drug = dataclasses.replace(simple_drug(), **changes)
    inhibitor = Substance(
        sid="inh",
        name="inhibitor",
        mr=400.0,
        distribution=Kp(1.0, ftissue=1.0),
        inhibits=(
            Inhibition(protein="CYP3A4", kind="mechanism", ki=1e-3, kinact=0.05),
        ),
    )
    d = definition(
        sid="balance",
        substances=[
            drug,
            Substance(sid="drg_m", name="metabolite", mr=316.0),
            inhibitor,
        ],
        chemistry=chemistry,
        routes=[Route.iv(S), Route.iv("inh")],
        variants=Variants(
            organs={
                Organ.LIVER: ZonatedLiver(
                    n_zones=n, profiles={"CYP3A4": [0.5 * (i + 1) for i in range(n)]}
                )
            }
        ),
    )
    r = _load(build(d, tmp_path).flat_path, f_shunt=0.3, f_functional=0.7)
    r[ids.iv_dose(S)] = 10.0
    r[ids.iv_dose("inh")] = 20.0
    r.simulate(0, DAY, 49)
    drug_amount = total_amount(r, S) + total_amount(r, "drg_m") + dose_amount(r, S)
    assert drug_amount == pytest.approx(10.0 / r[ids.mr(S)], rel=1e-7)
    inh_amount = total_amount(r, "inh") + dose_amount(r, "inh")
    assert inh_amount == pytest.approx(20.0 / r[ids.mr("inh")], rel=1e-7)
    assert total_amount(r, "drg_m") > 0
    if _BILE in chemistry:
        assert r[ids.feces(S)] > 0
    for upstream in _zones(n)[:-1]:
        assert r[_liver(ids.sinusoid_reaction(Part.PLASMA, S, upstream))] > 0
    assert r[_liver(ids.enzyme("CYP3A4"))] < 1.0
    model = libsbml.readSBMLFromString(r.getCurrentSBML()).getModel()
    turnover = libsbml.formulaToL3String(
        model.getRateRule(_liver(ids.enzyme("CYP3A4"))).getMath()
    )
    for z in _zones(n):
        assert _liver(ids.placement_weight("cyp3a4", z)) in turnover
        for part in _TISSUE[form]:
            assert r[_liver(ids.unbound(part, S, z))] > 0
        assert r[_liver(ids.unbound(Part.TISSUE, "inh", z))] > 0
