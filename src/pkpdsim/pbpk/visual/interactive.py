"""The interactive pages of the figures: hover cards on the static SVG.

One self-contained HTML file per figure (no external script or style), for the
docs and GitHub Pages: hovering an element with a `data-k` key shows its card,
an organ (its volumes in the reference man and woman, its UBERON ids, the form
of every substance, the processes), a substance (its mechanisms), a route, the
blood sample, the lumen or the urine; a key `txt:<text>` shows its own text
(a process line). The credit (pkpdsim, the anatomy and its license) sits in
the legend below the figure, never in a card.
"""

import json
from collections.abc import Mapping
from pathlib import Path
from xml.sax.saxutils import escape

from pkpdsim import __version__
from pkpdsim.pbpk.definition import PBPKDefinition
from pkpdsim.pbpk.figure import schematic_svg
from pkpdsim.pbpk.visual.anatomy import TITLES
from pkpdsim.pbpk.visual.animal import figure_of, reference_volumes
from pkpdsim.pbpk.visual.body import reference_body
from pkpdsim.pbpk.visual.content import ModelContent
from pkpdsim.pbpk.visual.resources import anatomy

__all__ = [
    "CREDIT",
    "LIBRARY",
    "REPOSITORY",
    "cards",
    "html_page",
    "interactive_figure",
]

REPOSITORY = "https://github.com/matthiaskoenig/pkpdsim"
LIBRARY = "https://matthiaskoenig.github.io/pkpdsim/pbpk/"
CREDIT = (
    f'Created with <a href="{REPOSITORY}">pkpdsim</a> {__version__} '
    f'(<a href="{LIBRARY}">PBPK library</a>) from the model definition; anatomy after the anatomograms '
    'of the <a href="https://www.ebi.ac.uk/gxa">Expression Atlas</a> and the '
    '<a href="https://www.proteinatlas.org">Human Protein Atlas</a> '
    '(<a href="https://github.com/ebi-gene-expression-group/anatomogram">anatomogram</a>, '
    '<a href="https://creativecommons.org/licenses/by/4.0/">CC BY 4.0</a>).'
)
"""The credit of a figure, for its legend (HTML)."""

_FORMS = {
    "none": "plasma only",
    "merged": "tissue (merged)",
    "full": "interstitium and cells",
    "organ": "one organ species",
}


def cards(
    content: ModelContent, volumes: Mapping[str, Mapping[str, float]]
) -> dict[str, list[str]]:
    """The hover cards of the figures by key, each a list of HTML lines.

    Args:
        content: the content of the definition.
        volumes: the organ volumes [l] (or masses [kg]) of the individuals
            the figure shows, by their name (`reference man`).

    Returns:
        The cards of the organs, substances (`sub:<sid>`), routes
        (`route:oral`, `route:iv`, `iv`), the sample (`sample`, `blood`), the
        lumen and the urine.
    """
    shapes = anatomy("human", "male")
    out: dict[str, list[str]] = {}
    for key, org in content.organs.items():
        rows = [f"<b>{escape(TITLES.get(key, key))}</b> <code>{key}</code>"]
        if org.variant:
            rows[0] += f" <i>({escape(org.variant)})</i>"
        sizes = [f"{v[key]:.3g} ({escape(n)})" for n, v in volumes.items() if key in v]
        if sizes:
            rows.append("volume [l]: " + ", ".join(sizes))
        shape = shapes.organs.get(key)
        if shape is not None:
            rows.append(
                "UBERON: " + ", ".join(u.replace("_", ":") for u in shape.uberon)
            )
        rows.append(
            "substances: "
            + "; ".join(f"{escape(s)}: {_FORMS[f]}" for s, f in org.parts.items())
        )
        rows += [
            f'<span class="b">{escape(b)}</span> {escape(t)}' for b, t in org.lines
        ]
        out[key] = rows
    for s in content.substances:
        out[f"sub:{s.sid}"] = [
            f"<b>{escape(s.sid)}</b> {escape(s.name)}",
            escape(s.details),
        ]
    out["route:oral"] = ["<b>oral dose</b>"] + [
        f"{k.replace('_', ' ')}: {escape(s)}"
        for k, s in content.routes
        if k.startswith("oral")
    ]
    out["route:iv"] = ["<b>intravenous dose</b>"] + [
        f"{escape(s)} into the venous blood" for k, s in content.routes if k == "iv"
    ]
    out["iv"] = out["route:iv"]
    out["sample"] = [
        "<b>blood sample</b>",
        "arm vein: the observer Carm of every substance, the default target of the fits",
    ]
    out["blood"] = out["sample"]
    if content.lumen_lines:
        out["lumen"] = ["<b>intestinal lumen</b>"] + [
            f'<span class="b">{escape(b)}</span> {escape(t)}'
            for b, t in content.lumen_lines
        ]
    out["urine"] = ["<b>urine</b>", "excreted: " + ", ".join(content.urine)]
    out["feces"] = [
        "<b>feces</b>",
        "from the last segment: " + ", ".join(content.luminal),
    ]
    for k in ("bladder", "gallbladder", "uterus"):
        out[k] = [f"<b>{k}</b>", "anatomy of the drawing, no organ of the model"]
    return out


def html_page(
    svg: str, cards: dict[str, list[str]], title: str, width: int, legend: str
) -> str:
    """A self-contained page of a figure with hover cards and its legend.

    Args:
        svg: the figure, whose elements carry `data-k` keys.
        cards: the cards by key.
        title: the title of the page.
        width: the width of the figure [px], the maximal width of the page.
        legend: the legend below the figure (HTML), the credit appended.

    The page posts its height to the page which embeds it in a frame
    (`{"pkpdFigure": true, "height": ...}`), so the frame can fit it
    (`docs/javascripts/figures.js`).

    Returns:
        The HTML, deterministic.
    """
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{escape(title)}</title>
<style>
body {{ margin:0; background:#fff; color:#1f2a2e; font:14px/1.45 system-ui, sans-serif; }}
.wrap {{ max-width:{width}px; margin:0 auto; padding:8px; }}
.wrap svg {{ width:100%; height:auto; display:block; }}
[data-k] {{ cursor:help; }}
path[data-k]:hover {{ filter:brightness(1.12) drop-shadow(0 0 2px #0006); }}
.legend {{ max-width:{width}px; margin:4px auto 16px; padding:0 8px; font-size:12.5px; color:#555; }}
.legend a {{ color:#00796b; }}
#card {{ position:fixed; pointer-events:none; max-width:380px; background:#fff; border:1px solid #d5dde0; border-radius:8px;
  box-shadow:0 6px 24px #0002; padding:8px 10px; font-size:12.5px; display:none; z-index:10; }}
#card div {{ margin:2px 0; }} #card code {{ font-size:11.5px; color:#666; }}
#card .b {{ display:inline-block; font:700 10px monospace; color:#fff; background:#6d6d6d; border-radius:3px; padding:0 4px; }}
</style></head><body><div class="wrap">{svg}</div>
<p class="legend">{legend}<br>{CREDIT}</p><div id="card" role="tooltip"></div>
<script>
const CARDS = {json.dumps(cards, sort_keys=True)};
const card = document.getElementById('card');
const esc = s => s.replace(/[&<>]/g, c => ({{'&':'&amp;','<':'&lt;','>':'&gt;'}})[c]);
document.querySelectorAll('.wrap [data-k]').forEach(el => {{
  el.addEventListener('mousemove', ev => {{
    const k = el.dataset.k;
    const rows = k.startsWith('txt:') ? [esc(k.slice(4))] : CARDS[k];
    if (!rows) return;
    ev.stopPropagation();
    card.innerHTML = rows.map(r => `<div>${{r}}</div>`).join('');
    card.style.display = 'block';
    card.style.left = Math.min(ev.clientX + 14, window.innerWidth - card.offsetWidth - 8) + 'px';
    card.style.top = Math.min(ev.clientY + 14, window.innerHeight - card.offsetHeight - 8) + 'px';
  }});
  el.addEventListener('mouseleave', () => card.style.display = 'none');
}});
// the height of the page for the frame which embeds it (the docs size their frames by it)
// (of its content: the scroll height of the document is at least the height of the frame)
const height = () => parent !== window && parent.postMessage(
  {{pkpdFigure: true, height: document.body.getBoundingClientRect().height + 16}}, '*');
addEventListener('load', height);
addEventListener('resize', height);
</script></body></html>
"""


def interactive_figure(
    definition: PBPKDefinition, path: Path, *, figure: str = "anatomy"
) -> Path:
    """Write the interactive page of a figure of a definition.

    Args:
        definition: the definition.
        path: the HTML file to write.
        figure: `anatomy` (the anatomical figure of both builds) or
            `schematic` (the figure of `pkpdsim.pbpk.figure`).

    Returns:
        The path of the page.

    Raises:
        ValueError: for another figure.
    """
    content = ModelContent.of(definition)
    species = str(definition.species)
    if species == "human":
        volumes = {
            "reference man": reference_body("male").quantities.volumes,
            "reference woman": reference_body("female").quantities.volumes,
        }
    else:
        volumes = {f"reference {species}": reference_volumes(species)}
    if figure == "anatomy":
        fig = figure_of(definition, content)
        svg = fig.svg(tags=True)
        width = fig.size[0]
        legend = "Hover the organs, labels, substances and symbols for details."
    elif figure == "schematic":
        svg = schematic_svg(definition, (definition.sex,))
        width = int(float(svg.split('width="', 1)[1].split('"', 1)[0]))
        legend = "Hover the organ rows, substance chips and process lines for details."
    else:
        raise ValueError(
            f"No interactive figure '{figure}'; the figures are 'anatomy' and 'schematic'."
        )
    page = html_page(
        svg, cards(content, volumes), f"{content.name} ({figure})", width, legend
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(page, encoding="utf-8")
    return path
