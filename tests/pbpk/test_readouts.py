"""The readouts of every key of a model figure."""

import dataclasses
import json
from pathlib import Path
from xml.etree import ElementTree

from pkpdsim.pbpk import Species, build, ids
from pkpdsim.pbpk.figure import model_figure
from pkpdsim.pbpk.organs import Organ
from pkpdsim.pbpk.visual.anatomy import anatomy_figure
from pkpdsim.pbpk.visual.animal import vein_of
from pkpdsim.pbpk.visual.content import ModelContent
from pkpdsim.pbpk.visual.interactive import (
    ReadoutGroup,
    interactive_figure,
    readouts,
    write_readouts,
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
            assert got[0] == f"[{pre}Cpl_{chip.sid}]"
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
        f"[{pre}{x}_{s}]"
        if x.startswith("C") and not x.startswith("Cu_")
        else f"{pre}{x}_{s}"
        for x in ("Cpl", "Cis", "Cic", "Cu_is", "Cu_ic")
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
    assert [r.id for r in got] == [f"[{ids.flat(Organ(o), f'Cpl_{s}')}]"]


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
    assert got[0] == "[LI__Cpl_api]"
    assert not any("_z1_" in i for i in got)


def test_the_vein_of_the_blood_line_follows_the_figure() -> None:
    content = ModelContent.of(seven_substances())
    human = readouts(content)
    rat = readouts(content, vein="tail vein")
    line = "blood/txt:OBS: {}: the observer Carm of every substance"
    assert line.format("arm vein") in human
    assert line.format("tail vein") in rat
    assert line.format("arm vein") not in rat


def _key_paths_of(svg: Path) -> set[str]:
    """The paths of the data-k of every element, the ancestors joined by `/`."""
    paths: set[str] = set()

    def walk(node: ElementTree.Element, prefix: tuple[str, ...]) -> None:
        key = node.get("data-k")
        here = (*prefix, key) if key is not None else prefix
        if key is not None:
            paths.add("/".join(here))
        for child in node:
            walk(child, here)

    walk(ElementTree.parse(svg).getroot(), ())
    return paths


def test_every_key_the_figure_draws_is_a_key_of_the_sidecar(tmp_path: Path) -> None:
    definition = seven_substances()
    svg = anatomy_figure(definition, tmp_path / "m_body.svg")
    sidecar = write_readouts(ModelContent.of(definition), svg)
    assert sidecar.name == "m_body.readouts.json"
    payload = json.loads(sidecar.read_text())
    assert payload["contract"] == "figurereadouts/1"
    assert payload["model"] == definition.sid
    drawn = _key_paths_of(svg)
    assert len(drawn) > 20
    assert any("/" in k for k in drawn)
    assert drawn <= set(payload["keys"])  # a dead hotspot is a defect
    group = payload["keys"]["li"]
    assert set(group) == {"label", "readouts"}
    assert all(set(r) == {"id", "label", "unit", "kind"} for r in group["readouts"])


def test_the_sidecar_follows_the_vein_of_an_animal_figure(tmp_path: Path) -> None:
    content = ModelContent.of(seven_substances())
    sidecar = write_readouts(content, tmp_path / "m.svg", vein="tail vein")
    keys = json.loads(sidecar.read_text())["keys"]
    assert any(k.startswith("blood/txt:OBS: tail vein") for k in keys)
    assert not any("arm vein" in k for k in keys)


def test_the_archive_writes_the_sidecar_beside_the_anatomy_figure(
    tmp_path: Path,
) -> None:
    from .test_archive import _create

    _create(tmp_path)
    sidecar = tmp_path / "models" / "simple_body.readouts.json"
    assert json.loads(sidecar.read_text())["model"] == "simple"
    assert not (tmp_path / "models" / "simple.readouts.json").exists()


def test_the_sidecar_covers_the_figure_of_an_animal(tmp_path: Path) -> None:
    definition = dataclasses.replace(seven_substances(), species=Species.RAT)
    svg = anatomy_figure(definition, tmp_path / "m_body.svg")
    sidecar = write_readouts(
        ModelContent.of(definition), svg, vein=vein_of(str(definition.species))
    )
    keys = set(json.loads(sidecar.read_text())["keys"])
    drawn = _key_paths_of(svg)
    assert any("tail vein" in k for k in drawn)
    assert drawn <= keys


def test_the_schematic_draws_keys_the_readouts_do_not_return(tmp_path: Path) -> None:
    """The reason the archive writes no sidecar of the schematic.

    When this fails, `readouts` covers the keys of the schematic: write its
    sidecar in `create_archive` and test it like the anatomy figure.
    """
    definition = seven_substances()
    svg = model_figure(definition, tmp_path / "m.svg")
    missing = _key_paths_of(svg) - set(readouts(ModelContent.of(definition)))
    assert "gu" in missing
    assert "ki/urine" in missing
    assert "lumen/feces" in missing


def test_a_concentration_readout_selects_a_concentration_of_the_model(
    tmp_path: Path,
) -> None:
    """Every readout names a selection whose values carry the unit it reports.

    The plain id of a species with `hasOnlySubstanceUnits="false"` selects its
    amount in `mmole`, so a concentration readout of such a species has to be
    the bracket form; a concentration which is a parameter is the plain id.
    Checked against the compiled model rather than against the id scheme,
    because the scheme alone cannot say which ids are species.
    """
    definition = seven_substances()
    sbml = build(definition, tmp_path).flat_path.read_text(encoding="utf-8")
    root = ElementTree.fromstring(sbml)
    ns = {"s": root.tag.split("}")[0].strip("{")}
    amounts = {
        e.get("id")
        for e in root.iter(f"{{{ns['s']}}}species")
        if e.get("hasOnlySubstanceUnits") == "false"
    }
    parameters = {e.get("id") for e in root.iter(f"{{{ns['s']}}}parameter")}

    seen = 0
    for group in readouts(ModelContent.of(definition)).values():
        for readout in group.readouts:
            if readout.kind != "concentration":
                continue
            bare = readout.id[1:-1] if readout.id.startswith("[") else readout.id
            if bare not in amounts and bare not in parameters:
                continue  # an id the model does not carry, the consumer filters it
            seen += 1
            if bare in amounts:
                assert readout.id == f"[{bare}]", (
                    f"{bare} is a species whose plain id selects an amount, so the "
                    f"concentration readout has to be its bracket form"
                )
            else:
                assert readout.id == bare, (
                    f"{bare} is a parameter, so its concentration is the plain id"
                )
    assert seen > 20, f"only {seen} concentration readouts were checked"
