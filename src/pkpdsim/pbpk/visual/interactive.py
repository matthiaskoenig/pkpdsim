"""The interactive figures: hover cards on the static SVG.

One self-contained HTML file per figure (no external script or style), for the
docs and GitHub Pages: hovering an element with a `data-k` key shows its card,
an organ (its volumes in the reference man and woman, its UBERON ids, the form
of every substance, the processes), a substance (its mechanisms), a route, the
blood sample, the lumen or the urine; a key `txt:<text>` shows its own text
(a process line). The credit (pkpdsim, the anatomy or the organ icons and the
license of the anatomograms, `CREDIT`, `ICON_CREDIT`) sits in
the legend below the figure, never in a card.

The file is a whole page (`html_page`), which a docs page embeds in a frame,
or a snippet (`html_snippet`, `interactive_figure(..., snippet=True)`), one
`<div>` which a page of markdown or HTML embeds as it is: everything of the
snippet, its style and the ids of its SVG, is scoped to that `<div>`, so two
snippets of one model stand on one page without touching it or each other.
The factory of a model writes the snippets beside its SVG figures
(`pkpdsim.pbpk.create_archive`) and the release of a model repository shows
them on its index page.
"""

import json
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal
from xml.sax.saxutils import escape

from pkpdsim import __version__
from pkpdsim.pbpk import ids
from pkpdsim.pbpk.definition import PBPKDefinition
from pkpdsim.pbpk.figure import schematic_svg
from pkpdsim.pbpk.organs import Organ, Part
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.visual.anatomy import (
    OBSERVER_TEXT,
    TITLES,
    line_key,
    substance_key,
)
from pkpdsim.pbpk.visual.animal import figure_of, reference_volumes
from pkpdsim.pbpk.visual.body import reference_body
from pkpdsim.pbpk.visual.content import ModelContent
from pkpdsim.pbpk.visual.resources import anatomy
from pkpdsim.sbml.values import unit_text
from pkpdsim.sex import Sex

__all__ = [
    "CREDIT",
    "ICON_CREDIT",
    "LIBRARY",
    "REPOSITORY",
    "Readout",
    "ReadoutGroup",
    "cards",
    "figure_credit",
    "html_page",
    "html_snippet",
    "interactive_figure",
    "readouts",
    "write_readouts",
]

REPOSITORY = "https://github.com/matthiaskoenig/pkpdsim"
LIBRARY = "https://matthiaskoenig.github.io/pkpdsim/pbpk/"


def figure_credit(
    drawn: str, *, library: str = "PBPK library", url: str = LIBRARY
) -> str:
    """Get the credit of a figure whose `drawn` part follows the anatomograms (HTML).

    Args:
        drawn: what of the figure follows the anatomograms (`anatomy`, `organ
            icons`).
        library: the name of the library of pkpdsim which drew the figure.
        url: the docs page of that library.

    Returns:
        The credit: pkpdsim with its version, the library with its page, and
        the attribution of the anatomograms with their license.
    """
    return (
        f'Created with <a href="{REPOSITORY}">pkpdsim</a> {__version__} '
        f'(<a href="{url}">{escape(library)}</a>) from the model definition; {drawn} after '
        'the anatomograms of the <a href="https://www.ebi.ac.uk/gxa">Expression Atlas</a> '
        'and the <a href="https://www.proteinatlas.org">Human Protein Atlas</a> '
        '(<a href="https://github.com/ebi-gene-expression-group/anatomogram">anatomogram</a>, '
        '<a href="https://creativecommons.org/licenses/by/4.0/">CC BY 4.0</a>).'
    )


CREDIT = figure_credit("anatomy")
"""The credit of an anatomical figure, for its legend (HTML)."""

ICON_CREDIT = figure_credit("organ icons")
"""The credit of a schematic, for its legend (HTML): its organ icons are the
shapes of the anatomograms (`pkpdsim.pbpk.visual.icons`), so the license asks
for the attribution, but the figure draws no anatomy."""

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


@dataclass(frozen=True)
class Readout:
    """One timecourse of the model which a place of the figure reports.

    Attributes:
        id: the id of the selection in the compiled model (`LI__Cpl_api`), as
            the id scheme gives it; it is not checked against a compiled model,
            the consumer filters the ids the model does not contain.
        label: what it is, in words (`apixaban, liver plasma`).
        unit: the unit the model uses for it (`mmole/liter`).
        kind: `concentration`, `amount`, `rate` or `flux`.
    """

    id: str
    label: str
    unit: str
    kind: Literal["concentration", "amount", "rate", "flux"]


@dataclass(frozen=True)
class ReadoutGroup:
    """What a key of the figure reports: its label and its readouts.

    Attributes:
        label: the name of the place, for the heading of its panel.
        readouts: the timecourses of the place, empty for a decoration of the
            drawing or a place of the model which carries no timecourse.
    """

    label: str
    readouts: tuple[Readout, ...]


_CONC = unit_text(U.mM)
_AMOUNT = unit_text(U.mmole)
_RATE = unit_text(U.mmole_per_min)


def _compartments(organ: Organ, sid: str, name: str, form: str) -> tuple[Readout, ...]:
    """The readouts of a substance in an organ, by the tissue form of the substance.

    Known limit: the ids are those of an organ without zones. A zonated organ
    (`OrganContent.variant` `zonated, N zones`) has the zone in its ids
    (`ids.species(..., zone=)`, `LI__Cpl_z1_api`), which are not given here.
    """
    title = TITLES.get(organ.value, organ.value)

    def conc(part: Part, text: str, *, unbound: bool = False) -> Readout:
        local = ids.unbound(part, sid) if unbound else ids.species(part, sid)
        return Readout(
            ids.flat(organ, local), f"{name}, {title} {text}", _CONC, "concentration"
        )

    out = [conc(Part.PLASMA, "plasma")]
    if form == "organ":
        out.append(conc(Part.ORGAN, "organ"))
    elif form == "merged":
        out += [
            conc(Part.TISSUE, "tissue"),
            conc(Part.TISSUE, "tissue unbound", unbound=True),
        ]
    elif form == "full":
        out += [
            conc(Part.INTERSTITIAL, "interstitium"),
            conc(Part.INTRACELLULAR, "cells"),
            conc(Part.INTERSTITIAL, "interstitium unbound", unbound=True),
            conc(Part.INTRACELLULAR, "cells unbound", unbound=True),
        ]
    return tuple(out)


def _sample(sid: str, name: str) -> tuple[Readout, ...]:
    """The readouts of a substance at the sampling site, the arm vein."""
    return (
        Readout(ids.arm(sid), f"{name}, plasma (sample)", _CONC, "concentration"),
        Readout(ids.arm_free(sid), f"{name}, unbound (sample)", _CONC, "concentration"),
        Readout(ids.arm_blood(sid), f"{name}, blood (sample)", _CONC, "concentration"),
    )


def _urine(sid: str, name: str) -> tuple[Readout, ...]:
    """The readouts of a substance in the urine, the amount and the rate."""
    return (
        Readout(ids.urine(sid), f"{name}, urine amount", _AMOUNT, "amount"),
        Readout(ids.urine_rate(sid), f"{name}, urine rate", _RATE, "rate"),
    )


def _feces(sid: str, name: str) -> tuple[Readout, ...]:
    """The readouts of a substance in the feces, the amount and the rate."""
    return (
        Readout(ids.feces(sid), f"{name}, feces amount", _AMOUNT, "amount"),
        Readout(ids.feces_rate(sid), f"{name}, feces rate", _RATE, "rate"),
    )


def _lumen(sid: str, name: str) -> tuple[Readout, ...]:
    """The readout of the dissolved amount of a substance in the lumen."""
    return (Readout(ids.lumen_amount(sid), f"{name}, lumen amount", _AMOUNT, "amount"),)


def readouts(
    content: ModelContent, *, vein: str = "arm vein"
) -> dict[str, ReadoutGroup]:
    """The timecourses which every key of the figure reports.

    The keys are the paths of the `data-k` keys of the figure from the outside
    in, joined by `/`, as the nesting of its SVG gives them: an organ (`li`),
    a substance of an organ (`li/sub:gli`), a process line of an organ
    (`li/txt:MET: ...`), the blood sample (`blood`, `blood/sub:gli`,
    `blood/txt:OBS: ...`), the legend chip of a substance (`sub:gli`), the
    sample icon (`sample`), the lumen, the urine and the feces with their
    substances, the routes and the decorations of the drawing (`bladder`).
    The three ways to point at a substance where it is sampled (`sub:<sid>`,
    `sample`, `blood/sub:<sid>`) give the same readouts, the arm vein
    observers.

    The compartments of a substance in an organ follow its tissue form
    (`OrganContent.parts`): the plasma, then for `organ` the one species, for
    `merged` the tissue and its unbound form, for `full` the interstitium and
    the cells with both unbound forms. A process line reports the flux of its
    reaction (`flat(organ, process.sid)`); a line which carries no process
    (`URI`, `INH`, `SLW`), a route and a decoration have no readout. An organ
    reports the compartments of every substance it shows. Ids are built with
    `pkpdsim.pbpk.ids`, the keys with the helpers of the figure writer
    (`substance_key`, `line_key`). The ids are what the id scheme gives, not
    checked against a compiled model: the consumer filters the ids the model
    does not contain (a zonated organ, whose ids carry a zone, is not
    supported, see `_compartments`).

    Args:
        content: the content of the definition.
        vein: the name of the sampling site in the line of the blood sample,
            `arm vein` of the human figure, that of the animal figures their
            own (`VEINS`).

    Returns:
        The group of every key, by the key path.
    """
    names = {s.sid: s.name or s.sid for s in content.substances}
    out: dict[str, ReadoutGroup] = {}

    def add(key: str, label: str, items: tuple[Readout, ...] = ()) -> None:
        out[key] = ReadoutGroup(label, items)

    def chips(
        base: str,
        label: str,
        sids: Sequence[str],
        per: Mapping[str, tuple[Readout, ...]],
    ) -> None:
        """Add the group of a place and its substance chips."""
        add(base, label, tuple(r for s in sids for r in per[s]))
        for s in sids:
            add(f"{base}/{substance_key(s)}", f"{names[s]}, {label}", per[s])

    for key, org in content.organs.items():
        organ = Organ(key)
        title = TITLES.get(key, key)
        shown = [c.sid for c in org.chips]
        per = {
            s: _compartments(organ, s, names[s], org.parts.get(s, "none"))
            for s in shown
        }
        chips(key, title, shown, per)
        flux = [p for p in content.processes if p.organ == key]
        for p, (badge, text) in zip(flux, org.lines, strict=False):
            add(
                f"{key}/{line_key(badge, text)}",
                f"{p.sid}, {title}",
                (
                    Readout(
                        ids.flat(organ, p.sid),
                        f"{p.sid}, {title} flux",
                        _RATE,
                        "flux",
                    ),
                ),
            )
        for badge, text in org.lines[len(flux) :]:
            add(f"{key}/{line_key(badge, text)}", title)

    sampled = {s.sid: _sample(s.sid, names[s.sid]) for s in content.substances}
    chips("blood", "blood sample", list(sampled), sampled)
    out["blood/" + line_key("OBS", OBSERVER_TEXT.replace("arm vein", vein, 1))] = out[
        "blood"
    ]
    add("sample", "blood sample", out["blood"].readouts)
    for s, items in sampled.items():
        add(substance_key(s), names[s], items)

    lumen = {s: _lumen(s, names[s]) for s in content.luminal}
    chips("lumen", "intestinal lumen", content.luminal, lumen)
    for badge, text in content.lumen_lines:
        add(f"lumen/{line_key(badge, text)}", "intestinal lumen")
    chips(
        "urine", "urine", content.urine, {s: _urine(s, names[s]) for s in content.urine}
    )
    chips(
        "feces",
        "feces",
        content.luminal,
        {s: _feces(s, names[s]) for s in content.luminal},
    )

    for key in ("iv", "route:iv"):
        add(key, "intravenous dose")
    for kind, sid in content.routes:
        if kind == "iv":
            add(
                "iv/" + line_key("IV", f"{sid} into the venous blood"),
                "intravenous dose",
            )
    add("route:oral", "oral dose")
    for key in ("bladder", "gallbladder", "uterus"):
        add(key, key)
    return out


CONTRACT = "figurereadouts/1"
"""The version of the payload of the sidecar `write_readouts` writes."""


def write_readouts(content: ModelContent, svg: Path, *, vein: str = "arm vein") -> Path:
    """Write the readouts of the keys of a figure beside it, for other tools.

    The sidecar `<svg stem>.readouts.json` is `{"contract": "figurereadouts/1",
    "model": <sid>, "keys": {<path>: {"label": str, "readouts": [{"id",
    "label", "unit", "kind"}]}}}`, the paths and groups of `readouts`, so a
    reader of the SVG needs no import of pkpdsim: it joins the `data-k` of the
    ancestors of a clicked element into a path and looks it up. Every path the
    figure draws is a key (a decoration has no readouts). The ids are those of
    the id scheme and are not checked against a compiled model, the reader
    filters the ones its model lacks.

    Args:
        content: the content of the model the figure shows.
        svg: the SVG file of the figure.
        vein: the sampling site named in the figure (`readouts`), the vein of
            the animal of an animal figure.

    Returns:
        The path of the sidecar.
    """
    payload = {
        "contract": CONTRACT,
        "model": content.sid,
        "keys": {
            key: {
                "label": group.label,
                "readouts": [
                    {"id": r.id, "label": r.label, "unit": r.unit, "kind": r.kind}
                    for r in group.readouts
                ],
            }
            for key, group in readouts(content, vein=vein).items()
        },
    }
    path = svg.with_suffix(".readouts.json")
    path.write_text(json.dumps(payload, indent=1) + "\n", encoding="utf-8")
    return path


def _script_json(value: object) -> str:
    """Serialize a value as the JSON of a `<script>` block of a page.

    An HTML parser ends a script element at the first `</script`, wherever it
    stands, also inside a JSON string, and a JavaScript string ends at a line
    separator (U+2028, U+2029) which JSON writes raw: the characters which can
    end the element or the string are written as their JSON escapes, so the
    block holds the same value and no text of a card can leave it.

    Args:
        value: what the block assigns, e.g. the cards.

    Returns:
        The JSON, with `<`, `>`, `&`, U+2028 and U+2029 escaped.
    """
    return (
        json.dumps(value, sort_keys=True)
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
        .replace("&", "\\u0026")
        .replace("\u2028", "\\u2028")
        .replace("\u2029", "\\u2029")
    )


#: the script of the hover cards, the elements with a `data-k` key of
#: `{query}` and the element `card` which shows the rows of a key
_CARDS_JS = """\
const esc = s => s.replace(/[&<>]/g, c => ({{'&':'&amp;','<':'&lt;','>':'&gt;'}})[c]);
{query}.forEach(el => {{
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
}});"""


def html_page(
    svg: str,
    cards: dict[str, list[str]],
    title: str,
    width: int,
    legend: str,
    credit: str = CREDIT,
) -> str:
    """A self-contained page of a figure with hover cards and its legend.

    Args:
        svg: the figure, whose elements carry `data-k` keys.
        cards: the cards by key.
        title: the title of the page.
        width: the width of the figure [px], the maximal width of the page.
        legend: the legend below the figure (HTML), the credit appended.
        credit: the credit (HTML), `CREDIT` of an anatomical figure,
            `ICON_CREDIT` of a schematic.

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
<p class="legend">{legend}<br>{credit}</p><div id="card" role="tooltip"></div>
<script>
const CARDS = {_script_json(cards)};
const card = document.getElementById('card');
{_CARDS_JS.format(query="document.querySelectorAll('.wrap [data-k]')")}
// the height of the page for the frame which embeds it (the docs size their frames by it)
// (of its content: the scroll height of the document is at least the height of the frame)
const height = () => parent !== window && parent.postMessage(
  {{pkpdFigure: true, height: document.body.getBoundingClientRect().height + 16}}, '*');
addEventListener('load', height);
addEventListener('resize', height);
</script></body></html>
"""


_STYLE = re.compile(r"(<style[^>]*>)(.*?)(</style>)", re.S)
_ID = re.compile(r'\bid="([^"]+)"')


def _scoped_css(css: str, scope: str) -> str:
    """Prefix every selector of a style sheet with the id of its container.

    The style sheets of the figures are generated and simple: class and type
    selectors, no at-rule, no comment and no nesting, so a rule is scoped by
    its selectors alone.

    Args:
        css: the style sheet.
        scope: the id of the element the rules are scoped to.

    Returns:
        The rules, each selector prefixed, on one line.

    Raises:
        ValueError: if the style sheet has an at-rule or a comment, which the
            scope cannot be written into.
    """
    if "@" in css or "/*" in css:
        raise ValueError(
            "The style sheet of the figure has an at-rule or a comment, which "
            f"the snippet cannot scope to '{scope}'."
        )
    rules = []
    for rule in css.split("}"):
        selectors, brace, body = rule.partition("{")
        if not brace:
            continue
        scoped = ",".join(f"#{scope} {s.strip()}" for s in selectors.split(","))
        rules.append(f"{scoped}{{{body.strip()}}}")
    return "".join(rules)


def _scoped_svg(svg: str, scope: str) -> str:
    """Scope the ids and the style sheet of an SVG to its container.

    An SVG in a page of HTML shares the ids and the style sheet of the page,
    so a snippet prefixes every id of the figure with the id of its container
    (`id`, `url(#...)`, `href`, `aria-labelledby`) and scopes every rule of its
    style sheet (`_scoped_css`).

    Args:
        svg: the figure.
        scope: the id of the container.

    Returns:
        The SVG, scoped.
    """
    for sid in sorted(set(_ID.findall(svg)), key=len, reverse=True):
        for text in (f'id="{sid}"', f"url(#{sid})", f'href="#{sid}"'):
            svg = svg.replace(text, text.replace(sid, f"{scope}-{sid}"))
        svg = svg.replace(
            f'aria-labelledby="{sid}"', f'aria-labelledby="{scope}-{sid}"'
        )
    return _STYLE.sub(
        lambda m: m[1] + _scoped_css(m[2], scope) + m[3],
        svg,
    )


def html_snippet(
    svg: str,
    cards: dict[str, list[str]],
    scope: str,
    width: int,
    legend: str,
    credit: str = CREDIT,
) -> str:
    """A snippet of a figure with hover cards and its legend, for a page.

    One `<div>` with the figure, its legend, its cards and its script, which a
    page of markdown or HTML embeds as it is, e.g. the index of a model
    repository. Everything of the snippet is scoped to the `<div>`: the style
    sheets of the snippet and of the SVG are written under its id and the ids
    of the SVG carry it (`_scoped_svg`), so a page shows two snippets without
    their styles meeting, and the script is a closure over the `<div>`. The
    snippet brings its own light background, so it reads in a dark theme as
    well, and it has no blank line, which would end it in markdown.

    Args:
        svg: the figure, whose elements carry `data-k` keys.
        cards: the cards by key.
        scope: the id of the `<div>`, unique on a page, e.g.
            `pkf-<model>-anatomy`.
        width: the width of the figure [px], the maximal width of the snippet.
        legend: the legend below the figure (HTML), the credit appended.
        credit: the credit (HTML), `CREDIT` of an anatomical figure,
            `ICON_CREDIT` of a schematic.

    Returns:
        The HTML of the snippet, deterministic.

    Raises:
        ValueError: if the style sheet of the SVG cannot be scoped
            (`_scoped_css`) or the snippet has a blank line.
    """
    style = (
        f"#{scope}{{display:block;max-width:{width}px;margin:1.2em auto;padding:8px;"
        "background:#fff;color:#1f2a2e;border-radius:8px;"
        "font:14px/1.45 system-ui, sans-serif}"
        # the figure alone: a nested SVG (an icon) keeps the size of its attributes
        f"#{scope} > svg{{width:100%;height:auto;display:block}}"
        f"#{scope} [data-k]{{cursor:help}}"
        f"#{scope} path[data-k]:hover{{filter:brightness(1.12) drop-shadow(0 0 2px #0006)}}"
        f"#{scope} .pkf-legend{{margin:4px 0 0;font-size:12.5px;color:#555}}"
        f"#{scope} .pkf-legend a{{color:#00796b}}"
        f"#{scope} .pkf-card{{position:fixed;pointer-events:none;max-width:380px;"
        "background:#fff;color:#1f2a2e;border:1px solid #d5dde0;border-radius:8px;"
        "box-shadow:0 6px 24px #0002;padding:8px 10px;font-size:12.5px;"
        "display:none;z-index:10}"
        f"#{scope} .pkf-card div{{margin:2px 0}}"
        f"#{scope} .pkf-card code{{font-size:11.5px;color:#666}}"
        f"#{scope} .pkf-card .b{{display:inline-block;font:700 10px monospace;"
        "color:#fff;background:#6d6d6d;border-radius:3px;padding:0 4px}"
    )
    script = _CARDS_JS.format(query="root.querySelectorAll('svg [data-k]')")
    snippet = f"""<div id="{scope}" class="pkf-figure">
<style>{style}</style>
{_scoped_svg(svg.strip(), scope)}
<p class="pkf-legend">{legend}<br>{credit}</p><div class="pkf-card" role="tooltip"></div>
<script>
(() => {{
const root = document.getElementById('{scope}');
const CARDS = {_script_json(cards)};
const card = root.querySelector('.pkf-card');
{script}
}})();
</script>
</div>
"""
    if any(not line.strip() for line in snippet.splitlines()):
        raise ValueError(
            f"The snippet '{scope}' has a blank line, which ends it in markdown."
        )
    return snippet


def interactive_figure(
    definition: PBPKDefinition,
    path: Path,
    *,
    figure: str = "anatomy",
    snippet: bool = False,
    sexes: Sequence[Sex] | None = None,
) -> Path:
    """Write the interactive page or snippet of a figure of a definition.

    Args:
        definition: the definition.
        path: the HTML file to write.
        figure: `anatomy` (the anatomical figure of both builds) or
            `schematic` (the figure of `pkpdsim.pbpk.figure`).
        snippet: write the snippet of the figure (`html_snippet`), the `<div>`
            which a page embeds, under the id `pkf-<sid>-<figure>`, instead of
            a whole page (`html_page`).
        sexes: the sexes which are built (`create_archive(..., sexes=...)`),
            which the figure names as the static figures do; `None` is the sex
            of the definition only.

    Returns:
        The path of the file.

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
        svg = fig.svg(tags=True, sexes=sexes)
        width = fig.size[0]
        legend = "Hover the organs, labels, substances and symbols for details."
        credit = CREDIT
    elif figure == "schematic":
        svg = schematic_svg(definition, sexes or (definition.sex,))
        width = int(float(svg.split('width="', 1)[1].split('"', 1)[0]))
        legend = "Hover the organ rows, substance chips and process lines for details."
        credit = ICON_CREDIT
    else:
        raise ValueError(
            f"No interactive figure '{figure}'; the figures are 'anatomy' and 'schematic'."
        )
    rows = cards(content, volumes)
    if snippet:
        text = html_snippet(
            svg, rows, f"pkf-{definition.sid}-{figure}", width, legend, credit
        )
    else:
        text = html_page(svg, rows, f"{content.name} ({figure})", width, legend, credit)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path
