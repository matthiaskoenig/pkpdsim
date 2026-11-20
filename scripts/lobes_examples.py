"""Generate the reference data and the browser models of the liver lobe examples.

Run ``uv run python scripts/lobes_examples.py`` after changing the lobes of
the liver, their remodelling or their interventions. The browser runs the same
compiled SBML through the same SUNDIALS IDA of sbmlsim-core as this script; it
holds no second, JavaScript model of a lobed liver. The precomputed courses of
this dataset are what the page shows before a reader moves a control, and the
test of the page compares them against a fresh run.

The runtime of the browser (`engine.wasm`, `mathlib.wasm`, `sbmlsim-core.js`)
is the one of the cardiovascular examples, `docs/pbpk/examples/solver`: both
exports come from the same installed release and the test asserts that they
agree, so the bytes of the engine are committed once.
"""

import argparse
import hashlib
import json
from importlib.metadata import version
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

import numpy as np
import sbmlsim_core
from sbmlsim.model import AbstractModel
from sbmlsim.model.model_sbml import IntegratorSettings
from sbmlsim.simulation import Simulation
from sbmlsim.simulator import Simulator

from pkpdsim.pbpk import (
    Binding,
    Excretion,
    FirstOrder,
    Kp,
    Metabolism,
    PBPKDefinition,
    PressureFlow,
    Route,
    Sex,
    Substance,
    Variants,
    build,
    ids,
)
from pkpdsim.pbpk.data import load
from pkpdsim.pbpk.interventions import (
    Intervention,
    PortalLigation,
    Resection,
    scenarios,
)
from pkpdsim.pbpk.lobes import HABR, lobe_table
from pkpdsim.pbpk.organs import Organ
from pkpdsim.pbpk.physiology import physiology_table
from pkpdsim.pbpk.remodelling import FlowDriven
from pkpdsim.pbpk.variants.liver_lobes import LobedLiver
from pkpdsim.species import Species

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "docs/pbpk/examples"
OUTPUT = EXAMPLES / "lobes-data.json"
SOLVERS = EXAMPLES / "solver_lobes"
RUNTIME = "solver"
"""The directory of the shared sbmlsim-core runtime, relative to the examples."""

S = "drg"
DAY = 24 * 60.0
END = 5 * DAY
POINTS = 121
DOSE = 10.0
SPECIES = (Species.RAT, Species.MOUSE, Species.HUMAN)
INTEGRATOR: IntegratorSettings = {
    "absolute_tolerance": 1e-10,
    "relative_tolerance": 1e-8,
}


def sexes(species: Species) -> tuple[Sex, ...]:
    """The sexes the parameter database has a reference physiology of.

    The human has both, the rat and the mouse have the male alone
    (`physiology_rat.toml`, `physiology_mouse.toml`), so the page offers the
    sex of an individual where the data of the library carries it and says so
    where it does not.
    """
    table = load(physiology_table(species))
    return tuple(sex for sex in Sex if sex.value in table)


def definition(species: Species, sex: Sex = Sex.MALE) -> PBPKDefinition:
    """The lobed liver of a species on the pressure network, with a toy drug.

    One drug cleared by the liver, so the exposure of a lobe says what the
    intervention did to the clearance, with the remodelling of the anatomy
    and the hepatic arterial buffer response both on: the explorer is about
    what a ligation and a resection do over the days which follow.
    """
    return PBPKDefinition(
        sid=f"lobes_{species.value}",
        species=species,
        sex=sex,
        substances=[
            Substance(
                sid=S,
                name="drug",
                mr=300.0,
                binding=Binding(fu=0.1, protein="albumin"),
                distribution=Kp(2.0, ftissue=1.0),
                excretion=Excretion(urine=True),
            ),
            Substance(sid="drg_m", name="metabolite", mr=316.0),
        ],
        chemistry=[
            Metabolism(
                sid="DRG2M",
                organ=Organ.LIVER,
                substrate=S,
                products={"drg_m": 1.0},
                kinetics=FirstOrder(k=0.5),
                protein="CYP3A4",
            )
        ],
        routes=[Route.iv(S)],
        variants=Variants(
            organs={
                Organ.LIVER: LobedLiver(remodelling=FlowDriven(), habr=HABR()),
            },
            circulation=PressureFlow(),
        ),
    )


def columns(species: Species) -> tuple[str, ...]:
    """The ids the page shows: the anatomy, the flows and the exposure per lobe."""
    tokens = [lobe.token for lobe in lobe_table(species)]
    return (
        "time",
        ids.VLI_LOBES,
        ids.F_REMNANT,
        "CV_P_po",
        "CV_P_li",
        *ids.HEPATIC_FLOWS,
        f"[{ids.pool_species(ids.Pool.VENOUS, S)}]",
        *(ids.lobe_volume(t) for t in tokens),
        *(ids.lobe_volume_fraction(t) for t in tokens),
        *(ids.lobe_perfusion(t) for t in tokens),
        *(ids.lobe_pressure(t) for t in tokens),
        *(ids.buffer_factor(t) for t in tokens),
        *(ids.lobe_flow(flow, t) for flow in ids.HEPATIC_FLOWS for t in tokens),
        *(ids.lobe_outflow(Organ.LIVER, t, S) for t in tokens),
    )


def cases(species: Species) -> dict[str, Intervention | None]:
    """The scenarios of the page: the intact liver and the named procedures.

    The procedures of the species come with it
    (`pkpdsim.pbpk.interventions.scenarios`); the ligation and the resection
    of the single largest lobe are added, because they are the smallest
    intervention which shows what the anatomy does with a time course.
    """
    lobes = lobe_table(species)
    largest = max(lobes, key=lambda lobe: lobe.mass_fraction.value).token
    found: dict[str, Intervention | None] = {"intact": None}
    found |= dict(scenarios(species))
    found[f"pvl_{largest}"] = PortalLigation([largest], name=f"pvl_{largest}")
    found[f"resection_{largest}"] = Resection([largest], name=f"resection_{largest}")
    return found


def _simulate(
    path: Path, selections: tuple[str, ...], changes: dict[str, float]
) -> dict[str, list[float]]:
    """The course of the five days after an intervention applied at the start."""
    result = Simulator(engine="core", n_workers=1, **INTEGRATOR).simulate(
        AbstractModel(source=path, selections=list(selections)),
        Simulation(
            end=END,
            steps=POINTS - 1,
            preinit_changes={ids.iv_dose(S): DOSE, **changes},
        ),
    )
    # seven significant digits are far more than a plot or a comparison of
    # the page needs and they keep the dataset small
    return {
        sid: [float(f"{float(v):.7g}") for v in np.asarray(result[sid])]
        for sid in selections
    }


def _export(path: Path, output: Path, species: Species) -> str:
    """Compile the SBML for the browser and write it beside the shared runtime.

    Only the model is written: the host script and the runtime of sbmlsim-core
    are the ones of the cardiovascular export, which the page loads from
    `../solver`.

    Returns:
        The SHA-256 of the SBML the model was compiled from.
    """
    output.mkdir(parents=True, exist_ok=True)
    sbml = path.read_text(encoding="utf-8")
    compiled = sbmlsim_core.compile_sbml(sbml)
    (output / "model.wasm").write_bytes(compiled.wasm)
    (output / "model.json").write_text(
        json.dumps(compiled.description, separators=(",", ":"), allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    fingerprint = hashlib.sha256(sbml.encode("utf-8")).hexdigest()
    config = {
        "generated_by": "scripts/lobes_examples.py",
        "engine": "SUNDIALS IDA",
        "species": species.value,
        "sbmlsim_core_version": version("sbmlsim-core"),
        "sbml_sha256": fingerprint,
        "runtime": f"../{RUNTIME}",
        "columns": list(columns(species)),
        "start": 0.0,
        "end": END,
        "points": POINTS,
        "dose": DOSE,
        "settings": {
            **INTEGRATOR,
            "absolute_tolerance": {"scalar": INTEGRATOR["absolute_tolerance"]},
            "maximum_num_steps": 1000000,
            "initial_time_step": 0,
            "maximum_time_step": 0,
            "stiff": True,
            "linear_solver": "auto",
            "jacobian": "auto",
        },
    }
    (output / "config.json").write_text(
        json.dumps(config, separators=(",", ":"), allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return fingerprint


def generate(*, sex: Sex = Sex.MALE, browser: Path | None = None) -> dict[str, Any]:
    """Build the lobed liver of every species and generate the dataset.

    Args:
        sex: the sex of the individual every model of the dataset is built for.
        browser: the directory of the compiled browser models, `None` to
            generate the dataset only.

    Returns:
        The dataset of the page.
    """
    data: dict[str, Any] = {
        "schema_version": 1,
        "provenance": {
            "generator": "scripts/lobes_examples.py",
            "engine": "core",
            "solver": "SUNDIALS IDA",
            "sex": sex.value,
            "sbmlsim_version": version("sbmlsim"),
            "sbmlsim_core_version": version("sbmlsim-core"),
            "end_time": END,
            "points": POINTS,
            "time_unit": "min",
            "dose": DOSE,
            "dose_target": ids.iv_dose(S),
            "substance": S,
            "notice": (
                "Illustrative lobed livers of a toy drug, not a reproduction "
                "of a measured hepatectomy or portal vein ligation study."
            ),
        },
        "species": {},
    }
    with TemporaryDirectory(prefix="pkpdsim-lobes-examples-") as directory:
        folder = Path(directory)
        for species in (name for name in SPECIES if sex in sexes(name)):
            result = build(definition(species, sex), folder / species.value)
            selections = columns(species)
            fingerprint = (
                _export(result.flat_path, browser / species.value, species)
                if browser is not None
                else hashlib.sha256(
                    result.flat_path.read_text(encoding="utf-8").encode("utf-8")
                ).hexdigest()
            )
            lobes = lobe_table(species)
            data["species"][species.value] = {
                "sbml_sha256": fingerprint,
                "sexes": [one.value for one in sexes(species)],
                "lobes": [
                    {
                        "token": lobe.token,
                        "name": lobe.name,
                        "mass_fraction": lobe.mass_fraction.value,
                        "source": lobe.mass_fraction.source,
                    }
                    for lobe in lobes
                ],
                "columns": list(selections),
                "cases": {
                    name: {
                        "label": name.replace("_", " "),
                        "changes": (
                            {}
                            if intervention is None
                            else intervention.changes([lobe.token for lobe in lobes])
                        ),
                        "kind": (
                            "intact"
                            if intervention is None
                            else type(intervention).__name__
                        ),
                        "lobes": (
                            [] if intervention is None else list(intervention.lobes)
                        ),
                        "series": _simulate(
                            result.flat_path,
                            selections,
                            {}
                            if intervention is None
                            else intervention.changes([lobe.token for lobe in lobes]),
                        ),
                    }
                    for name, intervention in cases(species).items()
                },
            }
    return data


def anatomy_of(species: Species, sex: Sex = Sex.MALE) -> dict[str, Any]:
    """The drawn animal of a species and a sex, its box and where its liver sits.

    The figure of the library (`pkpdsim.pbpk.visual`), which is the anatomogram
    of the EBI Expression Atlas for the mouse and the rat and the body of
    pkpdsim for the human: the explorer draws the same animal the physiology
    figures of the documentation draw, so a reader sees which individual the
    lobes belong to. The drawings are each scaled to their own box, so the
    species are not on one scale. The body of the human is drawn per sex; the
    anatomogram of the rat and of the mouse is one drawing for both, which is
    what the library holds, so there only the numbers of the run differ.
    """
    from pkpdsim.pbpk.visual.animal import AnimalBody
    from pkpdsim.pbpk.visual.body import reference_body

    if species is Species.HUMAN:
        body = reference_body(sex.value)
        drawing = body.svg(arteries=False, tags=True)
        points = np.asarray(body.outline)
        liver = body.centroid("li")
    else:
        animal = AnimalBody(species.value, sex.value)
        drawing = animal.svg(tags=True)
        points = np.asarray(animal.a.outline)
        liver = animal.centroid("li")
    low, high = points.min(0), points.max(0)
    return {
        "svg": drawing,
        "box": [
            round(float(v), 3)
            for v in (low[0], low[1], high[0] - low[0], high[1] - low[1])
        ],
        "liver": [round(float(v), 3) for v in liver],
    }


def dataset_path(output: Path, sex: Sex) -> Path:
    """Where the dataset of a sex lives: the male beside the page, the female under `female/`.

    The convention of the cardiovascular examples, which the page shares its
    runtime with: the explorer loads the dataset of the sex a reader asks for
    and not both, so a page costs what it costs today.
    """
    return output if sex is Sex.MALE else output.parent / "female" / output.name


def anatomy(output: Path = OUTPUT) -> dict[str, Any]:
    """The drawings of every species and sex, with the licence they travel under."""
    return {
        "generated_by": "scripts/lobes_examples.py",
        "attribution": "Anatomy adapted from EMBL-EBI Expression Atlas anatomograms (CC BY 4.0), using pkpdsim physiology figures.",
        "source": "https://github.com/ebi-gene-expression-group/anatomogram",
        "license": "https://creativecommons.org/licenses/by/4.0/",
        "bodies": {
            species.value: {
                sex.value: anatomy_of(species, sex) for sex in sexes(species)
            }
            for species in SPECIES
        },
    }


def _write_json(path: Path, content: dict[str, Any]) -> None:
    """Write one generated file of the page."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(content, separators=(",", ":"), ensure_ascii=False, allow_nan=False)
        + "\n",
        encoding="utf-8",
        newline="\n",
    )


def write(output: Path = OUTPUT, *, browser: Path | None = SOLVERS) -> None:
    """Write the datasets of both sexes, the drawings and the browser models.

    Args:
        output: the path of the dataset of the male individual.
        browser: the directory of the compiled browser models.
    """
    for sex in Sex:
        models = (
            None
            if browser is None
            else (
                browser if sex is Sex.MALE else browser.parent / "female" / browser.name
            )
        )
        _write_json(dataset_path(output, sex), generate(sex=sex, browser=models))
    _write_json(output.parent / "lobes-anatomy.json", anatomy(output))


def main() -> None:
    """Write the dataset and the browser models of the liver lobe examples."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="write the dataset only, not the compiled browser models",
    )
    args = parser.parse_args()
    write(args.output, browser=None if args.no_browser else SOLVERS)


if __name__ == "__main__":
    main()
