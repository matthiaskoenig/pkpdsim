"""Write the generated figures of the docs page `docs/pbpk/figures.md`.

The figures of the glimepiride-like model of the tutorial (`docs/pbpk/tutorial.md`):
the anatomical figure and the schematic, each static (SVG) and interactive (HTML),
the publication figure (the panels on an A4 page with the legend below it) and a
strip of individuals of the library. Run it after a change of the figures or the
anatomy resources:

```bash
uv run python scripts/visual/docs_figures.py
```

`tests/docs/test_docs_figures.py` fails while a file is out of date.
"""

import sys
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
    Substance,
    model_figure,
)
from pkpdsim.pbpk.visual.anatomy import anatomy_figure
from pkpdsim.pbpk.visual.body import Individual, body_strip
from pkpdsim.pbpk.visual.interactive import interactive_figure
from pkpdsim.pbpk.visual.publication import publication_figure
from pkpdsim.sex import Sex

OUT = Path(__file__).resolve().parents[2] / "docs" / "pbpk" / "img" / "figures"


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


def write(out: Path) -> list[Path]:
    """Write every figure of the page into a folder; get the files."""
    out.mkdir(parents=True, exist_ok=True)
    d = glimepiride()
    files = [
        anatomy_figure(d, out / "anatomy.svg"),
        interactive_figure(d, out / "anatomy.html"),
        model_figure(d, out / "schematic.svg"),
        interactive_figure(d, out / "schematic.html", figure="schematic"),
        publication_figure(d, out / "publication.svg"),
    ]
    # the legend of a manuscript is not a page of the docs, which shows it inline
    (out / "publication_legend.md").unlink()
    files.append(out / "publication_panels.svg")
    strip = out / "individuals.svg"
    strip.write_text(body_strip(individuals()) + "\n", encoding="utf-8")
    return [*files, strip]


def main() -> None:
    for path in write(OUT):
        print(path.relative_to(OUT.parents[3]))


if __name__ == "__main__":
    sys.exit(main())
