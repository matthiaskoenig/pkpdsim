"""Write the generated figures of the docs.

The figures of the page `docs/pbpk/figures.md`, of the glimepiride-like model of
the tutorial: the anatomical figure (A) and the schematic (B) as interactive pages,
the publication figure (the panels on an A4 page with the legend below it) and a
strip of individuals of the library; the anatomical figures of the species of
`docs/pbpk/species.md`, interactive as well.

Every page of `docs/pbpk/` whose python fences create a model shows the content of
that model: after each fence which creates a `PBPKDefinition` the script runs the
fences of the page up to it (in a temporary directory, as the snippet tests do),
takes the definition it created and writes its figures A and B as interactive pages
(`img/figures/models/<page>-<fence>-{anatomy,schematic}.html`) and the frames which
show them into a generated region after the fence. The docs show the interactive
figures only. Run it after a change of the figures, the anatomy resources or a
model of a page:

```bash
uv run python scripts/visual/docs_figures.py
```

`tests/docs/test_docs_figures.py` fails while a file is out of date.
"""

import ast
import contextlib
import dataclasses
import re
import sys
import tempfile
from pathlib import Path

from pkpdsim.pbpk import (
    MM,
    Binding,
    Excretion,
    FirstOrder,
    Kp,
    Metabolism,
    Organ,
    PBPKDefinition,
    Route,
    Species,
    Substance,
)
from pkpdsim.pbpk.visual.body import Individual, body_strip
from pkpdsim.pbpk.visual.interactive import interactive_figure
from pkpdsim.pbpk.visual.publication import publication_figure
from pkpdsim.sex import Sex

ROOT = Path(__file__).resolve().parents[2]
PBPK = ROOT / "docs" / "pbpk"
OUT = PBPK / "img" / "figures"
MODELS = "models"
"""The folder of the figures of the models of the pages, below `OUT`."""

BLOCK = re.compile(r"^```python\n(.*?)^```\n", re.S | re.M)
SKIP_MARKER = "# not executed"
OUTPUT = re.compile(r"\n```text\n.*?^```\n", re.S | re.M)
REGION = re.compile(
    r"\n<!-- generated: figures-[\w-]+ \(scripts/visual/docs_figures\.py\), "
    r"do not edit -->\n.*?^<!-- end generated -->\n",
    re.S | re.M,
)
_SOURCES = re.compile(rf'src="img/figures/{MODELS}/([\w.-]+)"')
NAMES = ("definition", "DEFINITION")
"""The names of the definition of a fence, preferred to another new definition."""


def glimepiride() -> PBPKDefinition:
    """The glimepiride-like definition of the tutorial, given orally and intravenously."""
    gli = Substance(
        "gli",
        name="glimepiride",
        mr=490.6,
        binding=Binding(fu=0.005, protein="albumin"),
        distribution=Kp(0.25, ftissue=50.0, ps={Organ.LIVER: 42.9}),
    )
    m1 = Substance(
        "gli_m1",
        name="hydroxy glimepiride (M1)",
        mr=506.6,
        excretion=Excretion(urine=True),
    )
    m2 = Substance(
        "gli_m2",
        name="carboxy glimepiride (M2)",
        mr=520.6,
        excretion=Excretion(urine=True),
    )
    chemistry = [
        Metabolism(
            "GLI2M1",
            Organ.LIVER,
            "gli",
            products={"gli_m1": 1.0},
            kinetics=MM(vmax=0.04, km=0.01),
            protein="CYP2C9",
        ),
        Metabolism(
            "M12M2",
            Organ.LIVER,
            "gli_m1",
            products={"gli_m2": 1.0},
            kinetics=FirstOrder(k=0.2),
        ),
    ]
    return PBPKDefinition(
        sid="glimepiride_body",
        name="Glimepiride body model",
        substances=[gli, m1, m2],
        chemistry=chemistry,
        routes=[Route.oral("gli"), Route.iv("gli")],
    )


def individuals() -> list[Individual]:
    """The individuals of the strip: the references, children, old age, obesity, cirrhosis."""
    return [
        Individual(Sex.MALE, label="reference man"),
        Individual(Sex.FEMALE, label="reference woman"),
        Individual(Sex.MALE, age=1.0, label="1 year"),
        Individual(Sex.FEMALE, age=5.0, label="5 years"),
        Individual(Sex.MALE, age=85.0, label="85 years"),
        Individual(Sex.MALE, bw=40 * 1.76**2, label="BMI 40"),
        Individual(Sex.MALE, liver=0.6655, label="Child-Pugh C"),
    ]


def _fences(text: str) -> list[tuple[int, str, int]]:
    """The python fences of a page: number (from 1), source, end of the fence."""
    return [
        (number, m.group(1), m.end())
        for number, m in enumerate(BLOCK.finditer(text), start=1)
    ]


def model_fences(page: Path) -> dict[int, PBPKDefinition]:
    """The definitions which the python fences of a page create, by fence.

    The fences run in order in one namespace in a temporary directory, as the
    snippet tests run them, up to the last fence which creates a definition; a
    fence creates the definition named `definition` (or `DEFINITION`) it binds,
    else the last new one.
    """
    fences = [
        (n, code)
        for n, code, _ in _fences(page.read_text(encoding="utf-8"))
        if not code.lstrip().startswith(SKIP_MARKER)
    ]
    last = max((n for n, code in fences if "PBPKDefinition(" in code), default=0)
    found: dict[int, PBPKDefinition] = {}
    namespace: dict[str, object] = {"__name__": "__main__"}
    seen: set[int] = set()
    with tempfile.TemporaryDirectory() as tmp, contextlib.chdir(tmp):
        for number, code in fences:
            if number > last:
                break
            if number == last:
                code = _up_to_definition(code)
            exec(compile(code, f"<{page.name} block {number}>", "exec"), namespace)
            new = {
                name: value
                for name, value in namespace.items()
                if isinstance(value, PBPKDefinition) and id(value) not in seen
            }
            seen.update(id(v) for v in new.values())
            if "PBPKDefinition(" in code and new:
                named = [new[k] for k in NAMES if k in new]
                found[number] = named[0] if named else list(new.values())[-1]
    return found


def _up_to_definition(code: str) -> str:
    """The statements of a fence up to the last one which creates a definition.

    The rest of the last fence (a build, a simulation) does not change the
    definition and is not run.
    """
    tree = ast.parse(code)
    lines = code.splitlines(keepends=True)
    ends = [
        node.end_lineno
        for node in tree.body
        if "PBPKDefinition(" in "".join(lines[node.lineno - 1 : node.end_lineno])
    ]
    return "".join(lines[: ends[-1]]) if ends else code


def model_pages() -> list[Path]:
    """The pages of the PBPK library with a python fence which creates a model."""
    return sorted(
        page
        for page in PBPK.glob("*.md")
        if any(
            "PBPKDefinition(" in code and not code.lstrip().startswith(SKIP_MARKER)
            for _, code, _ in _fences(page.read_text(encoding="utf-8"))
        )
    )


def _frame(src: str, title: str, html: str) -> str:
    """The frame of an interactive figure, its first size from the figure."""
    m = re.search(r'<svg[^>]* width="([\d.]+)" height="([\d.]+)"', html)
    ratio = f"{float(m.group(1)):.0f} / {float(m.group(2)) + 70:.0f}" if m else "4 / 3"
    return (
        f'<iframe class="pkpd-figure" src="{src}" title="{title}" '
        f'style="aspect-ratio:{ratio}" loading="lazy"></iframe>'
    )


def region(page: Path, number: int, definition: PBPKDefinition, out: Path) -> str:
    """Write the figures A and B of a definition of a page; get the region of the page."""
    stem = f"{page.stem}-{number}"
    frames = []
    for figure, title in (("anatomy", "anatomical figure"), ("schematic", "schematic")):
        path = interactive_figure(
            definition, out / MODELS / f"{stem}-{figure}.html", figure=figure
        )
        frames.append(
            _frame(
                f"img/figures/{MODELS}/{path.name}",
                f"The {title} of the model {definition.sid}",
                path.read_text(encoding="utf-8"),
            )
        )
    return (
        f"\n<!-- generated: figures-{stem} (scripts/visual/docs_figures.py), "
        "do not edit -->\n"
        f"The content of the model `{definition.sid}`: the anatomical figure (A) and "
        "the schematic (B), hover them for details ([figures](figures.md)).\n\n"
        + "\n\n".join(frames)
        + "\n<!-- end generated -->\n"
    )


def page_text(page: Path, out: Path) -> str:
    """The text of a page with the figures of its models after their fences."""
    text = REGION.sub("", page.read_text(encoding="utf-8"))
    definitions = model_fences(page)
    ends = {n: end for n, _, end in _fences(text)}
    for number in sorted(definitions, reverse=True):
        end = ends[number]
        # after the output of the fence, which the page pastes below it
        output = OUTPUT.match(text, end)
        if output:
            end = output.end()
        text = text[:end] + region(page, number, definitions[number], out) + text[end:]
    return text


def write(out: Path) -> list[Path]:
    """Write every figure of the docs into a folder; get the files.

    The pages of the models are not changed; `write_pages` writes them.
    """
    out.mkdir(parents=True, exist_ok=True)
    d = glimepiride()
    files = [
        interactive_figure(d, out / "anatomy.html"),
        interactive_figure(d, out / "schematic.html", figure="schematic"),
        publication_figure(d, out / "publication.svg"),
    ]
    # the legend of a manuscript is not a page of the docs, which shows it inline
    (out / "publication_legend.md").unlink()
    files.append(out / "publication_panels.svg")
    for species in (Species.MOUSE, Species.RAT, Species.PIG):
        files.append(
            interactive_figure(
                dataclasses.replace(d, species=species), out / f"{species.value}.html"
            )
        )
    strip = out / "individuals.svg"
    strip.write_text(body_strip(individuals()) + "\n", encoding="utf-8")
    return [*files, strip]


def write_pages(out: Path) -> dict[Path, str]:
    """Write the figures of the models of the pages into a folder; get the texts of the pages.

    A figure of a model which no page shows any more is removed.
    """
    folder = out / MODELS
    before = set(folder.glob("*.html")) if folder.exists() else set()
    texts = {page: page_text(page, out) for page in model_pages()}
    shown = {folder / src for text in texts.values() for src in _SOURCES.findall(text)}
    for stale in before - shown:
        stale.unlink()
    return texts


def main() -> None:
    for path in write(OUT):
        print(path.relative_to(ROOT))
    for page, text in write_pages(OUT).items():
        page.write_text(text, encoding="utf-8")
        print(page.relative_to(ROOT))


if __name__ == "__main__":
    sys.exit(main())
