"""Build the 2D anatomy resources of the figures from the anatomograms of the EBI Expression Atlas.

The anatomograms (https://github.com/ebi-gene-expression-group/anatomogram, the images under CC BY 4.0,
the code under Apache 2.0) are SVG drawings of the human male and female, the mouse male and female and
the rat, every tissue a group whose id is its UBERON term. This script downloads them at a pinned
commit (checked by sha256), converts every shape into absolute polylines and writes one JSON resource
per species and sex into `src/pkpdsim/resources/anatomy/` (read by `pkpdsim.pbpk.visual.resources`):

- the outline of the body and its details (eyes, nose, mouth, hands), simplified polylines;
- the organs of the library by their key of `pkpdsim.pbpk.Organ` (skin is the outline) with their UBERON
  ids, and further anatomy (urinary bladder, gallbladder, uterus, placenta, breast, diaphragm, aorta);
- the landmarks of the human drawings along the height (top, neck, shoulder, armpit, waist, crotch,
  bottom), the waist below the kidneys, the crotch below the bladder;
- the body cavities: the thorax above a diaphragm (the lower edge of the hull of lungs and heart) and
  the abdomen above a pelvic floor (a bowl rising to the iliac crests), inside the trunk.

It needs packages which are not dependencies of pkpdsim:

```bash
uv run --with svgelements --with lxml --with shapely python scripts/visual/anatomy_assets.py
```
"""

import hashlib
import io
import json
import re
import sys
import urllib.request
from pathlib import Path

import numpy as np
import svgelements as se  # ty: ignore[unresolved-import]
from lxml import etree  # ty: ignore[unresolved-import]
from PIL import Image, ImageDraw
from shapely.geometry import LineString, Point, Polygon, box  # ty: ignore[unresolved-import]
from shapely.ops import unary_union  # ty: ignore[unresolved-import]

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "src" / "pkpdsim" / "resources" / "anatomy"
CACHE = Path.home() / ".cache" / "pkpdsim" / "anatomogram"
COMMIT = "9fcc37022cce1e2862692a5f5fbfb78572b87e67"
URL = (
    "https://raw.githubusercontent.com/ebi-gene-expression-group/anatomogram/"
    f"{COMMIT}/src/svg/{{name}}.svg"
)
#: resource key: (file of the anatomogram, species, sex, sha256 of the file)
SOURCES = {
    "human_male": (
        "homo_sapiens.male",
        "human",
        "male",
        "c19874f0ec6c6450525277d97679bfe57cb0fd280a920210d0961442d43162ab",
    ),
    "human_female": (
        "homo_sapiens.female",
        "human",
        "female",
        "895ecda0934595950746a13c16abb8bd49ae4a4f947d575878f2f95b33b18e28",
    ),
    "mouse_male": (
        "mus_musculus.male",
        "mouse",
        "male",
        "3c25d8e5d66e92d2a0d4d23764af0f0388e4f7a5244393e105ca18f93b3aaba6",
    ),
    "mouse_female": (
        "mus_musculus.female",
        "mouse",
        "female",
        "f2d277f18754994b019e723114357da8c97a105140e498796017867ddf7f163d",
    ),
    "rat_male": (
        "rattus_norvegicus",
        "rat",
        "male",
        "b70efcfaf104f100b0945a77750af3cec52b97c82f8c7331e5316413b329140e",
    ),
}
#: the organs of the library (`pkpdsim.pbpk.Organ`) by their UBERON terms and their cavity
ORGANS = {
    "br": ("brain", ["UBERON_0000955"], "head"),
    "lu": ("lung", ["UBERON_0002048"], "thorax"),
    "he": ("heart", ["UBERON_0000948"], "thorax"),
    "li": ("liver", ["UBERON_0002107"], "abdomen"),
    "ki": ("kidney", ["UBERON_0002113"], "abdomen"),
    "st": ("stomach", ["UBERON_0000945"], "abdomen"),
    "si": (
        "small intestine",
        ["UBERON_0002108", "UBERON_0002114", "UBERON_0002116", "UBERON_0002115"],
        "abdomen",
    ),
    "la": (
        "large intestine",
        [
            "UBERON_0001155",
            "UBERON_0001153",
            "UBERON_0001052",
            "UBERON_0001154",
            "UBERON_0001157",
        ],
        "abdomen",
    ),
    "sp": ("spleen", ["UBERON_0002106"], "abdomen"),
    "pa": ("pancreas", ["UBERON_0001264"], "abdomen"),
    "go": ("gonads", ["UBERON_0000473", "UBERON_0000992"], None),
    "mu": ("muscle", ["UBERON_0001134", "UBERON_0014892", "UBERON_0001630"], "limb"),
    # UBERON_00024818 (sic, eight digits) is the tibia of the human drawings
    "bo": (
        "bone",
        ["UBERON_0002371", "UBERON_00024818", "UBERON_0002481", "UBERON_0001474"],
        "limb",
    ),
    "ad": ("adipose tissue", ["UBERON_0001013", "UBERON_0001347"], None),
}
#: further anatomy of the drawings, no organs of the library
ANATOMY = {
    "bladder": ("urinary bladder", ["UBERON_0001255"], "abdomen"),
    "gallbladder": ("gallbladder", ["UBERON_0002110"], "abdomen"),
    "uterus": ("uterus", ["UBERON_0000995"], "abdomen"),
    "placenta": ("placenta", ["UBERON_0001987"], "abdomen"),
    "breast": ("breast", ["UBERON_0000310"], None),
    "diaphragm": ("diaphragm", ["UBERON_0001103"], None),
    "aorta": ("aorta", ["UBERON_0000947"], None),
}
ATTRIBUTION = (
    "Anatomograms of the Expression Atlas (EMBL-EBI) and the Human Protein Atlas, "
    "https://github.com/ebi-gene-expression-group/anatomogram, CC BY 4.0; converted to "
    "polylines, simplified, tissues grouped by the organs of pkpdsim, landmarks and cavities "
    "derived (pkpdsim, scripts/visual/anatomy_assets.py)."
)
SVG_NS = "{http://www.w3.org/2000/svg}"
RES = 0.25  # units of the drawing per pixel of the raster of the landmarks


def download(name: str, sha256: str) -> bytes:
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{COMMIT[:8]}_{name}.svg"
    if not path.exists():
        with urllib.request.urlopen(URL.format(name=name), timeout=60) as r:
            path.write_bytes(r.read())
    data = path.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    if digest != sha256:
        raise ValueError(f"{name}: sha256 {digest}, expected {sha256}")
    return data


def annotate(data: bytes) -> bytes:
    """Mark every shape with the UBERON term of its closest group; drop what svgelements cannot reify."""
    parser = etree.XMLParser(resolve_entities=False, no_network=True)
    root = etree.parse(io.BytesIO(data), parser).getroot()
    root.attrib.pop("style", None)
    for tag in ("use", "metadata"):
        for el in list(root.iter(SVG_NS + tag)):
            el.getparent().remove(el)
    for el in root.iter():
        if not isinstance(el.tag, str):
            continue
        if el.tag.replace(SVG_NS, "") not in (
            "path",
            "ellipse",
            "circle",
            "rect",
            "polygon",
            "polyline",
        ):
            continue
        organ, layer, p = None, None, el
        while p is not None:
            i = p.get("id") or ""
            if organ is None and re.match(r"^(UBERON|CL|EFO)_\d+$", i):
                organ = i
            if i.startswith("LAYER_"):
                layer = i
            p = p.getparent()
        el.set("data-organ", organ or "")
        el.set("data-layer", layer or "")
    return etree.tostring(root)


def parts(data: bytes) -> dict[str, list[np.ndarray]]:
    """Every shape as absolute polylines (transforms applied), by its UBERON term or its layer."""
    sys.setrecursionlimit(20000)
    svg = se.SVG.parse(io.BytesIO(annotate(data)), reify=True)
    out: dict[str, list[np.ndarray]] = {}
    for el in svg.elements():
        if not isinstance(el, se.Shape):
            continue
        key = el.values.get("data-organ") or el.values.get("data-layer") or ""
        p = el if isinstance(el, se.Path) else se.Path(el)
        for sub in p.as_subpaths():
            sp = se.Path(sub)
            length = sp.length(error=1e-3) if len(sp) else 0
            if not length:
                continue
            n = max(8, int(length / 0.3))
            pts = np.array(
                [(q.x, q.y) for q in (sp.point(t) for t in np.linspace(0, 1, n))]
            )
            out.setdefault(key, []).append(pts)
    return out


def polygon(q: np.ndarray) -> Polygon | None:
    p = Polygon(q).buffer(0)
    return p if p.area > 1e-4 else None


def union(shapes: dict[str, list[np.ndarray]], uids: list[str]) -> Polygon | None:
    ps = [polygon(q) for u in uids for q in shapes.get(u, []) if len(q) > 2]
    ps = [p for p in ps if p is not None]
    return unary_union(ps) if ps else None


def to_d(geom: Polygon | None, tol: float) -> str:
    if geom is None or geom.is_empty:
        return ""
    geom = geom.simplify(tol, preserve_topology=True)
    out = []
    for g in getattr(geom, "geoms", [geom]):
        if not isinstance(g, Polygon) or g.area < 1e-3:
            continue
        for ring in [g.exterior, *g.interiors]:
            c = np.asarray(ring.coords)[:-1]
            out.append("M" + " ".join(f"{x:.2f},{y:.2f}" for x, y in c) + "Z")
    return "".join(out)


def outline(shapes: dict[str, list[np.ndarray]]) -> tuple[Polygon, Polygon | None]:
    """The body (the largest closed shape of the outline layer) and the small details inside it."""
    subs = [q for q in shapes["LAYER_OUTLINE"] if len(q) > 3]
    body_pts = max(subs, key=lambda q: Polygon(q).buffer(0).area)
    body = Polygon(body_pts).buffer(0)
    details = [polygon(q) for q in subs if q is not body_pts]
    # the small shapes inside the body (eyes, nose, mouth, hands), not the icons beside it
    details = [
        d
        for d in details
        if d is not None and d.area < 6.0 and body.contains(d.centroid)
    ]
    return body, (unary_union(details) if details else None)


def runs(row: np.ndarray) -> list[tuple[float, float]]:
    r, i, n = [], 0, len(row)
    while i < n:
        if row[i]:
            j = i
            while j < n and row[j]:
                j += 1
            r.append((i * RES, j * RES))
            i = j
        else:
            i += 1
    return r


def landmarks(
    body: Polygon, shapes: dict[str, list[np.ndarray]]
) -> tuple[dict[str, float], dict[float, tuple[float, float]]]:
    """The landmarks of a human drawing and the trunk run of every raster row."""
    x0, y0, x1, y1 = body.bounds
    w, h = int((x1 + 2) / RES), int((y1 + 2) / RES)
    im = Image.new("L", (w, h), 0)
    ImageDraw.Draw(im).polygon(
        [(x / RES, y / RES) for x, y in body.exterior.coords], fill=255
    )
    mask = np.array(im) > 0
    cx = (x0 + x1) / 2
    ys = np.arange(mask.shape[0]) * RES
    top, bottom = y0, y1
    height = bottom - top
    width = np.array(
        [sum(b - a for a, b in runs(mask[i])) for i in range(mask.shape[0])]
    )
    sel = (ys > top + 0.08 * height) & (ys < top + 0.2 * height)
    neck = float(ys[sel][np.argmin(width[sel])])
    armpit = float(
        next(
            y
            for i, y in enumerate(ys)
            if y > neck + 0.03 * height and len(runs(mask[i])) >= 3
        )
    )

    def ymax(uid: str) -> float:
        return float(np.vstack(shapes[uid])[:, 1].max())

    crotch = ymax("UBERON_0001255") + 0.025 * height  # below the urinary bladder
    waist = (
        ymax("UBERON_0002113") + 0.03 * height
    )  # below the lower pole of the kidneys
    trunk = {}
    for i, y in enumerate(ys):
        if armpit <= y <= crotch:
            for a, b in runs(mask[i]):
                if a <= cx <= b:
                    trunk[round(float(y), 3)] = (a, b)
    marks = {
        "top": top,
        "neck": neck,
        "shoulder": neck + 0.35 * (armpit - neck),
        "armpit": armpit,
        "waist": waist,
        "crotch": crotch,
        "bottom": bottom,
    }
    return {k: round(float(v), 2) for k, v in marks.items()}, trunk


def cavities(
    body: Polygon,
    shapes: dict[str, list[np.ndarray]],
    species: str,
    marks: dict[str, float] | None,
    trunk: dict[float, tuple[float, float]] | None,
) -> tuple[Polygon, Polygon]:
    x0, y0, x1, y1 = body.bounds
    thor = unary_union(
        [
            g
            for g in (union(shapes, ORGANS["lu"][1]), union(shapes, ORGANS["he"][1]))
            if g is not None
        ]
    )
    if species == "human" and marks is not None and trunk is not None:
        floor = marks["crotch"] - 0.04 * (marks["crotch"] - marks["waist"])
        rows = [box(a, y, b, y + 0.25) for y, (a, b) in trunk.items() if y < floor]
        trunk_rows = unary_union(rows).buffer(0.3).buffer(-0.3)
        cage = thor.convex_hull.buffer(2.2)
        upper = cage.intersection(body).intersection(
            box(x0, y0, x1, marks["armpit"] + 0.25)
        )
        cav = unary_union([trunk_rows, upper]).intersection(body.buffer(-0.6))
        cx = (x0 + x1) / 2
        hw = 0.5 * (thor.bounds[2] - thor.bounds[0]) + 2
        xs = np.linspace(cx - 3 * hw, cx + 3 * hw, 121)
        bowl = Polygon(
            [(x, floor - 7.0 * min(1.0, ((x - cx) / hw) ** 2)) for x in xs]
            + [(xs[-1], y1 + 5), (xs[0], y1 + 5)]
        )
        cav = (
            cav.difference(bowl)
            .intersection(box(x0, y0, x1, floor))
            .buffer(1.0)
            .buffer(-1.0)
        )
        hull = thor.convex_hull.buffer(1.4)
        hx0, _, hx1, _ = hull.bounds
        dx = np.linspace(hx0 + 0.05, hx1 - 0.05, 80)
        dy = [hull.intersection(box(x - 0.01, y0, x + 0.01, y1)).bounds[3] for x in dx]
        line = [(x0 - 5, dy[0]), *zip(dx, dy, strict=False), (x1 + 5, dy[-1])]
        above = Polygon([*line, (x1 + 5, y0 - 5), (x0 - 5, y0 - 5)])
        thorax = cav.intersection(above)
    else:
        trunk_organs = [
            union(shapes, ORGANS[k][1])
            for k in ("lu", "he", "li", "ki", "st", "si", "la", "sp", "pa")
        ]
        block = unary_union([g for g in trunk_organs if g is not None])
        cav = block.convex_hull.buffer(1.2).intersection(body.buffer(-0.4))
        thorax = cav.intersection(thor.convex_hull.buffer(0.8))
    abdomen = cav.difference(thorax).buffer(0.2).buffer(-0.2)
    return thorax, abdomen


def build(key: str) -> dict:
    name, species, sex, sha256 = SOURCES[key]
    shapes = parts(download(name, sha256))
    body, details = outline(shapes)
    x0, y0, x1, y1 = body.bounds
    res: dict = {
        "id": key,
        "species": species,
        "sex": sex,
        "source": {
            "file": f"{name}.svg",
            "commit": COMMIT,
            "license": "CC BY 4.0",
            "attribution": ATTRIBUTION,
        },
        "viewbox": [round(x0, 2), round(y0, 2), round(x1 - x0, 2), round(y1 - y0, 2)],
        "midline": round((x0 + x1) / 2, 2),
    }
    marks, trunk = landmarks(body, shapes) if species == "human" else (None, None)
    if marks is not None:
        res["landmarks"] = marks
    thorax, abdomen = cavities(body, shapes, species, marks, trunk)
    res["cavities"] = {"thorax": to_d(thorax, 0.03), "abdomen": to_d(abdomen, 0.03)}
    res["outline"] = to_d(body, 0.03)
    res["details"] = to_d(details, 0.02)
    for table, field in ((ORGANS, "organs"), (ANATOMY, "anatomy")):
        res[field] = {}
        for k, (label, uids, cavity) in table.items():
            g = union(shapes, uids)
            if g is not None:
                res[field][k] = {
                    "name": label,
                    "uberon": [u for u in uids if u in shapes],
                    "cavity": cavity,
                    "d": to_d(g, 0.04),
                }
    return res


# --- the pig: an own silhouette (dorsal view) with the organs of the human drawings ------------
#: the legs of the pig: x, y, side and the angle of the leg [degrees]
PIG_LEGS = ((15, 66, -1, 25), (65, 66, 1, -25), (16, 120, -1, -20), (64, 120, 1, 20))
#: the organs of the trunk of the pig, fitted from the human drawing of the sex
PIG_TRUNK = ("lu", "he", "li", "st", "sp", "pa", "si", "la", "ki")
PIG_LAYOUT = {"ki": (1.4, 6.5), "he": (0.0, 4.0)}
PIG_ATTRIBUTION = (
    "Silhouette of the pig drawn by pkpdsim (scripts/visual/anatomy_assets.py, CC BY 4.0); the organs "
    "after the human anatomograms of the Expression Atlas (EMBL-EBI) and the Human Protein Atlas, "
    "https://github.com/ebi-gene-expression-group/anatomogram, CC BY 4.0, fitted into its trunk."
)


def ellipse(cx: float, cy: float, rx: float, ry: float, rot: float = 0.0) -> Polygon:
    from shapely.affinity import rotate, scale  # ty: ignore[unresolved-import]

    e = scale(Point(cx, cy).buffer(1.0, 64), rx, ry)
    return rotate(e, rot, origin=(cx, cy)) if rot else e


def pig_silhouette() -> Polygon:
    """A pig in dorsal view, the snout up: an oval trunk, a head tapering to the snout disc, floppy
    ears, the legs at the sides and a curly tail."""
    trunk = ellipse(40, 92, 26, 45)
    head = Polygon(
        [(24, 58), (26, 42), (33, 22), (36, 17), (44, 17), (47, 22), (54, 42), (56, 58)]
    )
    snout = ellipse(40, 11, 7.5, 5.2)
    ears = [
        Polygon([(31, 40), (14, 30), (16, 22), (33, 31)]).buffer(2.2),
        Polygon([(49, 40), (66, 30), (64, 22), (47, 31)]).buffer(2.2),
    ]
    legs = [ellipse(x, y, 5.5, 9.5, rot=r) for x, y, _, r in PIG_LEGS]
    body = unary_union([trunk, head, snout, *ears, *legs]).buffer(2.5).buffer(-2.5)
    t = np.linspace(0, 2.6 * np.pi, 160)
    r = 1.2 + 2.6 * t / (2.6 * np.pi)
    tail = LineString(np.c_[40 + r * np.sin(t), 140 + r * np.cos(t) + t * 0.55]).buffer(
        0.9
    )
    return unary_union([body, LineString([(40, 135), (40, 141)]).buffer(1.3), tail])


def _polylines(d: str) -> list[np.ndarray]:
    return [
        np.array(
            [[float(v) for v in p.split(",")] for p in sub.strip().rstrip("Z").split()]
        )
        for sub in d.split("M")[1:]
    ]


def _d(polys: list[np.ndarray]) -> str:
    return "".join(
        "M" + " ".join(f"{x:.2f},{y:.2f}" for x, y in q) + "Z" for q in polys
    )


def pig(human: dict, sex: str) -> dict:
    """The resource of the pig of a sex: the own silhouette, the organs of the human drawing."""
    body = pig_silhouette()
    if not isinstance(body, Polygon):
        body = max(body.geoms, key=lambda g: g.area)
    details = [
        ellipse(37.6, 10.5, 1.0, 1.6),
        ellipse(42.4, 10.5, 1.0, 1.6),
    ]  # the nostrils
    details += [ellipse(34.5, 27, 1.2, 1.0), ellipse(45.5, 27, 1.2, 1.0)]  # the eyes
    for x, y, _, r in PIG_LEGS:  # the split hooves
        a = np.radians(-r)
        tip = np.array([x, y]) + 8.5 * np.array([np.sin(a), np.cos(a)])
        details += [
            ellipse(tip[0] - 1.1, tip[1], 1.0, 1.6, rot=r),
            ellipse(tip[0] + 1.1, tip[1], 1.0, 1.6, rot=r),
        ]
    cx = human["midline"]
    pieces = {}
    for k in (
        *PIG_TRUNK,
        "gallbladder",
        "bladder",
        *(("uterus", "go") if sex == "female" else ()),
    ):
        entry = human["organs"].get(k) or human["anatomy"].get(k)
        if entry is None:
            continue
        polys = _polylines(entry["d"])
        if k in PIG_LAYOUT:
            dx, dy = PIG_LAYOUT[k]
            polys = [
                q + np.array([dx * np.sign(q[:, 0].mean() - cx), dy]) for q in polys
            ]
        pieces[k] = polys
    block = np.vstack([q for v in pieces.values() for q in v])
    lo, hi = block.min(0), block.max(0)
    k_block = min(42 / (hi[0] - lo[0]), 80 / (hi[1] - lo[1]))
    centre, target = (lo + hi) / 2, np.array([40.0, 92.0])
    moved = {k: [(q - centre) * k_block + target for q in v] for k, v in pieces.items()}
    brain = _polylines(human["organs"]["br"]["d"])
    moved["br"] = [
        (q - np.vstack(brain).mean(0)) * 0.75 + np.array([40.0, 40.0]) for q in brain
    ]
    if sex == "male":
        testes = _polylines(human["organs"]["go"]["d"])
        moved["go"] = [
            (q - np.vstack(testes).mean(0)) * 1.6 + np.array([40.0, 131.0])
            for q in testes
        ]
    organs, anatomy = {}, {}
    for k, polys in moved.items():
        src = human["organs"].get(k) or human["anatomy"].get(k)
        entry = {**src, "d": _d(polys)}
        (organs if k in human["organs"] else anatomy)[k] = entry
    x0, y0, x1, y1 = body.bounds
    return {
        "id": f"pig_{sex}",
        "species": "pig",
        "sex": sex,
        "source": {
            "file": "",
            "commit": COMMIT,
            "license": "CC BY 4.0",
            "attribution": PIG_ATTRIBUTION,
        },
        "viewbox": [round(x0, 2), round(y0, 2), round(x1 - x0, 2), round(y1 - y0, 2)],
        "midline": 40.0,
        # the points of the routes: the snout (oral gavage) and the ear vein (dose and sample)
        "landmarks": {"snout_x": 40.0, "snout_y": 1.5, "vein_x": 64.0, "vein_y": 27.0},
        "cavities": {"thorax": "", "abdomen": ""},
        "outline": to_d(body, 0.03),
        "details": to_d(unary_union(details), 0.02),
        "organs": organs,
        "anatomy": anatomy,
    }


NOTICE = """# Anatomy resources

The files of this folder are derived from the anatomograms of the Expression Atlas (EMBL-EBI), which the Human Protein Atlas uses as well: https://github.com/ebi-gene-expression-group/anatomogram (commit {commit}). The anatomograms are distributed under the Creative Commons Attribution 4.0 International license (CC BY 4.0, https://creativecommons.org/licenses/by/4.0/), the code of that repository under the Apache License 2.0.

Changes: the drawings were converted into simplified polylines, the tissues grouped by the organs of pkpdsim, the landmarks and the body cavities derived, by `scripts/visual/anatomy_assets.py` of pkpdsim. The pig (`pig_male.json`, `pig_female.json`) is a silhouette drawn by pkpdsim with the organs of the human drawings fitted into its trunk. The resources are distributed under CC BY 4.0 as well; a figure made from them names the anatomograms in its legend.
"""


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    built = {}
    for key in SOURCES:
        built[key] = build(key)
    for sex in ("male", "female"):
        built[f"pig_{sex}"] = pig(built[f"human_{sex}"], sex)
    for key, res in built.items():
        text = json.dumps(res, separators=(",", ":"), sort_keys=True)
        (OUT / f"{key}.json").write_text(text + "\n", encoding="utf-8")
        print(
            key, len(text) // 1024, "kB", sorted(res["organs"]), sorted(res["anatomy"])
        )
    (OUT / "NOTICE.md").write_text(NOTICE.format(commit=COMMIT), encoding="utf-8")


if __name__ == "__main__":
    sys.exit(main())
