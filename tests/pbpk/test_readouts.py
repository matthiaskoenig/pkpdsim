"""The readouts of every key of a model figure."""

from pathlib import Path
from xml.etree import ElementTree

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.organs import Organ
from pkpdsim.pbpk.visual.content import ModelContent
from pkpdsim.pbpk.visual.interactive import (
    ReadoutGroup,
    interactive_figure,
    readouts,
)

from .test_figure import seven_substances


def _groups() -> dict[str, ReadoutGroup]:
    return readouts(ModelContent.of(seven_substances()))


def test_an_organ_reports_its_substance_in_the_compartments_of_its_form() -> None:
    groups = _groups()
    content = ModelContent.of(seven_substances())
    for organ in ("li", "ki", "mu"):
        for chip in content.organs[organ].chips:
            f = content.organs[organ].parts[chip.sid]
            got = [r.id for r in groups[f"{organ}/sub:{chip.sid}"].readouts]
            pre = ids.flat(Organ(organ), "")
            assert got[0] == f"{pre}Cpl_{chip.sid}"
            n = {"none": 1, "organ": 2, "merged": 3, "full": 5}[f]
            assert len(got) == n
    forms = {f for o in content.organs.values() for f in o.parts.values()}
    assert {"none", "full"} <= forms


def test_a_full_organ_reports_both_unbound_forms() -> None:
    content = ModelContent.of(seven_substances())
    hits = [
        (o, s)
        for o, org in content.organs.items()
        for s, f in org.parts.items()
        if f == "full" and s in {c.sid for c in org.chips}
    ]
    o, s = hits[0]
    pre = ids.flat(Organ(o), "")
    got = _groups()[f"{o}/sub:{s}"].readouts
    assert [r.id for r in got] == [
        f"{pre}{x}_{s}" for x in ("Cpl", "Cis", "Cic", "Cu_is", "Cu_ic")
    ]
    assert {r.kind for r in got} == {"concentration"}
    assert {r.unit for r in got} == {"mmole/liter"}


def test_a_plasma_only_substance_reports_its_plasma_alone() -> None:
    content = ModelContent.of(seven_substances())
    hits = [
        (o, s)
        for o, org in content.organs.items()
        for s, f in org.parts.items()
        if f == "none" and s in {c.sid for c in org.chips}
    ]
    o, s = hits[0]
    got = _groups()[f"{o}/sub:{s}"].readouts
    assert [r.id for r in got] == [ids.flat(Organ(o), f"Cpl_{s}")]


def test_the_three_ways_to_point_at_the_sample_agree() -> None:
    groups = _groups()
    want = ["Carm_api", "Carm_free_api", "Carm_blood_api"]
    for key in ("sub:api", "blood/sub:api"):
        assert [r.id for r in groups[key].readouts] == want
    content = ModelContent.of(seven_substances())
    every = [
        i
        for sub in content.substances
        for i in (
            f"Carm_{sub.sid}",
            f"Carm_free_{sub.sid}",
            f"Carm_blood_{sub.sid}",
        )
    ]
    assert [r.id for r in groups["sample"].readouts] == every
    assert [r.id for r in groups["blood"].readouts] == every


def test_the_urine_reports_the_amount_and_the_rate() -> None:
    got = _groups()["urine/sub:api"].readouts
    assert [(r.id, r.kind, r.unit) for r in got] == [
        ("Aurine_api", "amount", "mmole"),
        ("Rurine_api", "rate", "mmole/min"),
    ]


def test_the_lumen_reports_the_dissolved_amount() -> None:
    groups = _groups()
    assert [r.id for r in groups["lumen/sub:api"].readouts] == ["Alumen_api"]


def test_a_process_badge_reports_the_flux_of_its_reaction() -> None:
    content = ModelContent.of(seven_substances())
    groups = _groups()
    p = next(p for p in content.processes if p.kind == "MET")
    lines = content.organs[p.organ].lines
    i = [q for q in content.processes if q.organ == p.organ].index(p)
    key = f"{p.organ}/txt:{lines[i][0]}: {lines[i][1]}"
    got = groups[key].readouts
    assert [(r.id, r.kind) for r in got] == [(ids.flat(Organ(p.organ), p.sid), "flux")]


def test_a_line_without_a_process_has_no_readout() -> None:
    groups = _groups()
    uri = [k for k in groups if k.startswith("ki/txt:URI")]
    assert uri
    assert all(groups[k].readouts == () for k in uri)


def test_a_decoration_of_the_drawing_has_no_readout() -> None:
    assert _groups()["bladder"].readouts == ()


def _key_paths(root: ElementTree.Element) -> set[str]:
    out: set[str] = set()

    def walk(el: ElementTree.Element, path: tuple[str, ...]) -> None:
        k = el.get("data-k")
        here = (*path, k) if k else path
        if k:
            out.add("/".join(here))
        for child in el:
            walk(child, here)

    walk(root, ())
    return out


def test_every_key_path_of_the_drawn_figure_has_a_group(tmp_path: Path) -> None:
    d = seven_substances()
    html = interactive_figure(d, tmp_path / "f.html", snippet=True).read_text()
    svg = html[html.index("<svg") : html.rindex("</svg>") + 6]
    paths = _key_paths(ElementTree.fromstring(svg))
    assert paths
    missing = paths - set(_groups())
    assert not missing


def test_a_zonated_organ_reports_the_ids_without_a_zone() -> None:
    """Known limit: a zonated liver has zone ids (`LI__Cpl_z1_api`), not given here."""
    content = ModelContent.of(seven_substances())
    assert "zonated" in (content.organs["li"].variant or "")
    got = [r.id for r in _groups()["li/sub:api"].readouts]
    assert got[0] == "LI__Cpl_api"
    assert not any("_z1_" in i for i in got)


def test_the_vein_of_the_blood_line_follows_the_figure() -> None:
    content = ModelContent.of(seven_substances())
    human = readouts(content)
    rat = readouts(content, vein="tail vein")
    line = "blood/txt:OBS: {}: the observer Carm of every substance"
    assert line.format("arm vein") in human
    assert line.format("tail vein") in rat
    assert line.format("arm vein") not in rat
