"""The publication figure of a model: the anatomical figure and the schematic on one A4 page.

Panel A is the anatomical figure of both builds in compact form (the substances
and the submodels, without the process text), panel B the schematic in compact
form (the substance dots and the key in its top band), both at the width of the
page; only the panels have to fit on the page. The legend sits below them, a
caption generated from the content which explains every symbol and names the
credit, also written as Markdown for the caption of a manuscript.
"""

import re
from pathlib import Path
from xml.sax.saxutils import escape

from pkpdsim import __version__
from pkpdsim.pbpk.definition import PBPKDefinition
from pkpdsim.pbpk.figure import schematic_svg
from pkpdsim.pbpk.visual.animal import figure_of
from pkpdsim.pbpk.visual.body import reference_body
from pkpdsim.pbpk.visual.content import ModelContent
from pkpdsim.pbpk.visual.interactive import LIBRARY, REPOSITORY

__all__ = ["PAGE", "legend", "publication_figure"]

PAGE = (210.0, 297.0, 10.0)
"""The page: width, height and margin [mm]."""

_FONT = "'DejaVu Sans', 'Noto Sans', sans-serif"
_TYPE = 2.3  # the type of the legend [mm]


def _inner(svg: str) -> tuple[float, float, str]:
    """The width, height and content of an SVG document."""
    head = re.match(r'<svg[^>]*width="([\d.]+)"[^>]*height="([\d.]+)"', svg)
    if head is None:
        raise ValueError("An SVG with a width and a height is expected.")
    width, height = float(head.group(1)), float(head.group(2))
    body = re.sub(r"^.*?<svg[^>]*>", "", svg, count=1, flags=re.S)
    return width, height, body[: body.rindex("</svg>")]


def legend(content: ModelContent) -> list[tuple[str, str]]:
    """The legend of the publication figure: paragraphs of (kind, text), kind `b` bold, `i` the credit."""
    man, woman = reference_body("male").quantities, reference_body("female").quantities
    substances = "; ".join(f"{s.sid}, {s.name}" for s in content.substances)
    routes = content.route_text[0].upper() + content.route_text[1:]
    pd = " ".join(
        f"PD component {sid}: reads {', '.join(reads) or 'nothing'}; sets {', '.join(sets) or 'nothing'}."
        for sid, reads, sets in content.pd
    )
    return [
        ("b", f"The PBPK model {content.name}."),
        (
            "",
            f"(A) Anatomical overview of the male and the female build on the reference man ({man.height:.0f} cm, "
            f"{man.bw:.0f} kg) and woman ({woman.height:.0f} cm, {woman.bw:.0f} kg). The organs keep their "
            "anatomical forms and are sized by their volumes in the model; subcutaneous and visceral adipose "
            "tissue, muscle, bones and arteries are layers of the body. On the body: portal vein (purple), bile "
            "(green), urine (yellow), feces (brown, dashed), oral route (orange, dashed); circled symbols mark the "
            "oral dose (tablet), the intravenous dose (syringe) and the blood sample at the arm vein (drop). "
            "Labels: the substances of an organ (filled: enters the tissue; dashed: plasma only); organs with "
            "processes show their submodel with a row per substance and the compartments plasma, interstitium, "
            "cells and the apical side (bile, urine, lumen), arrows for exchange, metabolism (enzyme), transport "
            "(transporter) and filtration.",
        ),
        (
            "",
            "(B) Model structure (the substances and the key of chips and badges at the top): venous and "
            "arterial blood and the portal vein with the circulating substances, the lung, the systemic organs, "
            "the liver, the portal organs and the intestinal lumen in its segments, urine and feces as "
            "compartments; per organ the substances and the processes: MET metabolism, TRA transport, URI urinary "
            "excretion, ORL oral dose, ABS absorption, BIL biliary excretion, FEC fecal excretion.",
        ),
        ("", f"Substances: {substances}. {routes}. {pd}".strip()),
        (
            "i",
            f"Created with pkpdsim {__version__} (PBPK library) from the model definition; anatomy after the "
            "anatomograms of the Expression Atlas and the Human Protein Atlas (CC BY 4.0).",
        ),
    ]


def _wrap(text: str, n: int) -> list[str]:
    out, cur = [], ""
    for w in text.split():
        if cur and len(cur) + 1 + len(w) > n:
            out.append(cur)
            cur = w
        else:
            cur = f"{cur} {w}" if cur else w
    return [*out, cur]


def _linked(text: str) -> str:
    """The credit with its links (SVG)."""
    style = 'fill="#00796b" text-decoration="underline"'
    t = escape(text)
    t = t.replace("pkpdsim", f'<a href="{REPOSITORY}" {style}>pkpdsim</a>', 1)
    return t.replace(
        "(PBPK library)", f'(<a href="{LIBRARY}" {style}>PBPK library</a>)', 1
    )


def publication_figure(
    definition: PBPKDefinition, path: Path, *, pdf: bool = False
) -> Path:
    """Write the publication figure of a definition.

    Writes `<path>` (the panels on the A4 page and the legend below it),
    `<stem>_panels.svg` (the A4 page of the panels), `<stem>_legend.md` (the
    legend as Markdown) and, with `pdf`, `<stem>_panels.pdf`.

    Args:
        definition: the definition.
        path: the SVG file of the figure with its legend.
        pdf: also write the PDF of the panels (needs the extra `pdf`).

    Returns:
        The path of the figure.

    Raises:
        ImportError: for `pdf` without cairosvg (`pip install pkpdsim[pdf]`).
    """
    content = ModelContent.of(definition)
    a = figure_of(definition, content, compact=True).svg(title=False, substances=False)
    b = schematic_svg(definition, (definition.sex,), compact=True)
    wa, ha, ba = _inner(a)
    wb, hb, bb = _inner(b)
    page_w, page_h, margin = PAGE
    width = page_w - 2 * margin
    ya = margin
    h_a = ha * width / wa
    yb = ya + h_a + 4
    h_b = page_h - margin - yb
    w_b = wb * h_b / hb
    if w_b > width:
        w_b, h_b = width, hb * width / wb
    xb = margin + (width - w_b) / 2
    panels = (
        f'<rect width="{page_w}" height="{page_h}" fill="#fff"/>'
        f'<text x="3" y="{ya + 3:.1f}" font-size="5.5" font-weight="700">A</text>'
        f'<svg x="{margin}" y="{ya:.1f}" width="{width:.1f}" height="{h_a:.1f}" viewBox="0 0 {wa:.0f} {ha:.0f}">{ba}</svg>'
        f'<text x="3" y="{yb + 3:.1f}" font-size="5.5" font-weight="700">B</text>'
        f'<svg x="{xb:.1f}" y="{yb:.1f}" width="{w_b:.1f}" height="{h_b:.1f}" viewBox="0 0 {wb:.0f} {hb:.0f}">{bb}</svg>'
    )
    stem = path.with_suffix("")
    page = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{page_w}mm" height="{page_h}mm" '
        f'viewBox="0 0 {page_w} {page_h}" font-family="{_FONT}">{panels}</svg>'
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    Path(f"{stem}_panels.svg").write_text(page + "\n", encoding="utf-8")
    paragraphs = legend(content)
    n = int(width / (_TYPE * 0.52))
    lines = []
    y = page_h + 6
    for kind, text in paragraphs:
        for line in _wrap(text, n):
            attr = (
                ' font-weight="700"'
                if kind == "b"
                else (' font-style="italic" fill="#555"' if kind == "i" else "")
            )
            body = _linked(line) if kind == "i" else escape(line)
            lines.append(
                f'<text x="{margin}" y="{y:.1f}" font-size="{_TYPE}"{attr}>{body}</text>'
            )
            y += _TYPE * 1.3
        y += 1.0
    total = y + 6
    path.write_text(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{page_w}mm" height="{total:.0f}mm" '
        f'viewBox="0 0 {page_w} {total:.0f}" font-family="{_FONT}"><rect width="{page_w}" height="{total:.0f}" fill="#fff"/>'
        f'{panels}<path d="M0,{page_h} h{page_w}" stroke="#bbb" stroke-width="0.3" stroke-dasharray="2,1.5"/>'
        f"{''.join(lines)}</svg>\n",
        encoding="utf-8",
    )
    md = []
    for kind, text in paragraphs:
        if kind == "b":
            md.append(f"**{text}**")
        elif kind == "i":
            linked = text.replace("pkpdsim", f"[pkpdsim]({REPOSITORY})", 1)
            md.append(
                "*"
                + linked.replace("(PBPK library)", f"([PBPK library]({LIBRARY}))", 1)
                + "*"
            )
        else:
            md.append(text)
    Path(f"{stem}_legend.md").write_text("\n\n".join(md) + "\n", encoding="utf-8")
    if pdf:
        try:
            import cairosvg
        except ImportError as err:
            raise ImportError(
                "The PDF of the publication figure needs cairosvg: pip install pkpdsim[pdf]"
            ) from err
        cairosvg.svg2pdf(bytestring=page.encode(), write_to=f"{stem}_panels.pdf")
    return path
