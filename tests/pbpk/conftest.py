"""Shared helpers of the PBPK tests.

`simple_drug` and `definition` are plain functions (not fixtures), so a test
varies them with arguments: `definition(routes=[Route.iv("drg")])`.
"""

import difflib
import gzip
import re
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pytest
from sbmlutils.factory import AssignmentRule, Model, Parameter, Port, PortType
from sbmlutils.metadata import BQB, SBO

from pkpdsim.pbpk import (
    MM,
    MODEL_UNITS,
    Binding,
    Erythrocytes,
    Excretion,
    Kp,
    Metabolism,
    Mode,
    Organ,
    PBPKDefinition,
    PDComponent,
    Route,
    Substance,
    U,
    ids,
)
from pkpdsim.pbpk.builder import BuildResult

if TYPE_CHECKING:
    import roadrunner


def simple_drug() -> Substance:
    """Get the simple drug: Mr 300, fu 0.1 (albumin), Kp 2, merged tissue, urine."""
    return Substance(
        sid="drg",
        name="drug",
        mr=300.0,
        binding=Binding(fu=0.1, protein="albumin"),
        distribution=Kp(2.0, ftissue=1.0),
        excretion=Excretion(urine=True),
    )


def definition(**overrides: Any) -> PBPKDefinition:  # noqa: ANN401
    """Get the definition of the simple drug with a liver metabolism `DRG2M`.

    Args:
        **overrides: fields of `PBPKDefinition` which replace the defaults.

    Returns:
        The definition.
    """
    fields: dict[str, Any] = {
        "sid": "simple",
        "substances": [
            simple_drug(),
            Substance(sid="drg_m", name="metabolite", mr=316.0),
        ],
        "chemistry": [
            Metabolism(
                sid="DRG2M",
                organ=Organ.LIVER,
                substrate="drg",
                products={"drg_m": 1.0},
                kinetics=MM(vmax=1e-3, km=0.01),
                protein="CYP3A4",
            )
        ],
        "routes": [Route.iv("drg"), Route.oral("drg")],
    }
    fields.update(overrides)
    return PBPKDefinition(**fields)


def cardiac_pd(
    sid: str = "cardio", emax: float = 0.5, ec50: float = 1e-4, reads: str = "drg"
) -> PDComponent:
    """Get a PD component which lowers the cardiac function with the drug.

    `f_cardiac = 1 - Emax*C/(EC50 + C)` of the unbound arm vein concentration
    `C` of the substance `reads`; the output replaces `f_cardiac_function`.
    Its elements are annotated, as a package would annotate its PD model.

    Args:
        sid: the id of the component.
        emax: the maximal reduction of the cardiac function.
        ec50: the concentration of half the reduction (mM).
        reads: the substance whose unbound arm vein concentration it reads.

    Returns:
        The component.
    """
    pd = [(BQB.IS, "ncit/C15720")]  # pharmacodynamics
    model = Model(
        sid,
        units=U,
        model_units=MODEL_UNITS,
        objects=[
            Parameter(
                "C_drug",
                0.0,
                U.mM,
                constant=False,
                name="unbound concentration",
                sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                annotations=pd,
                port=Port("C_drug_port", portType=PortType.INPUT_PORT),
            ),
            Parameter(
                "Emax",
                emax,
                U.dimensionless,
                sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                annotations=pd,
            ),
            Parameter("EC50", ec50, U.mM, sboTerm=SBO.EC50, annotations=pd),
            Parameter(
                "f_cardiac",
                None,
                U.dimensionless,
                constant=False,
                name="cardiac function",
                sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                annotations=pd,
                port=Port("f_cardiac_port", portType=PortType.OUTPUT_PORT),
            ),
            AssignmentRule(
                "f_cardiac",
                "1 dimensionless - Emax * C_drug / (EC50 + C_drug)",
                U.dimensionless,
            ),
        ],
    )
    return PDComponent(
        model=model,
        inputs={"C_drug": ids.arm_free(reads)},
        outputs={"f_cardiac": ids.F_CARDIAC_FUNCTION},
    )


def load(result: BuildResult) -> "roadrunner.RoadRunner":
    """Load the flat model of a build as sbmlsim does, with tight integrator tolerances.

    `load_roadrunner` of sbmlsim is the production path: LLJIT, always compiled.
    """
    from sbmlsim.model.model_roadrunner import load_roadrunner

    r = load_roadrunner(result.flat_path)
    r.integrator.absolute_tolerance = 1e-12
    r.integrator.relative_tolerance = 1e-10
    return r


def individual(path: Path, **changes: float) -> "roadrunner.RoadRunner":
    """Load a model with changes before its initialization, as an experiment does.

    The changes are the pre-initialization assignments of sbmlsim
    (`RoadrunnerSBMLModel.initialize`), which evaluates the initial
    assignments which read a changed parameter again.

    Args:
        path: the SBML file.
        **changes: the values of parameters by id, e.g. `BW=100.0`.

    Returns:
        The roadrunner instance of the initialized model.
    """
    from sbmlsim.model.model_roadrunner import RoadrunnerSBMLModel
    from sbmlsim.model.symbols import TargetKind
    from sbmlsim.simulator.plan import Assignment

    model = RoadrunnerSBMLModel(source=path)
    model.initialize(
        [Assignment(k, TargetKind.PARAMETER, value=v) for k, v in changes.items()]
    )
    return model.r_loaded


def total_amount(r: "roadrunner.RoadRunner", sid: str) -> float:
    """Get the total amount (mmole) of a substance: all species and the sinks.

    The amount of a species is decided from the SBML: a species with
    `hasOnlySubstanceUnits` is an amount, any other a concentration which is
    multiplied with the volume of its compartment. The species of the id
    scheme (`ids.species_of`: organs, pools, intestinal lumen, undissolved
    drug, mucosa) and the sinks: the urine and feces amounts and the bile;
    the dose states are not counted (`dose_amount`).
    """
    import libsbml

    doc = libsbml.readSBMLFromString(r.getCurrentSBML())
    model = doc.getModel()
    by_id = {sp.getId(): sp for sp in model.getListOfSpecies()}
    sinks = {ids.urine(sid), ids.feces(sid), ids.bile(sid)}
    selected = ids.species_of(sid, by_id) + [
        i for i in by_id if i.rsplit("__", 1)[-1] in sinks
    ]
    total = 0.0
    for i in selected:
        sp = by_id[i]
        if sp.getHasOnlySubstanceUnits():
            total += float(r[i])
        else:
            total += float(r[f"[{i}]"]) * float(r[sp.getCompartment()])
    return total


def dose_amount(r: "roadrunner.RoadRunner", sid: str) -> float:
    """Get the amount (mmole) of a substance in its dose states, not yet applied."""
    doses = (ids.iv_dose(sid), ids.po_dose(sid), ids.sol_dose(sid))
    mg = sum(float(r[d]) for d in doses if d in r.model.getGlobalParameterIds())
    return mg / float(r[ids.mr(sid)])


REFERENCE = Path(__file__).parent / "data"
"""The comp and flat SBML of `reference_definition()` built at pkpdsim f008d56, gzipped."""

_DATE = re.compile(r"<dcterms:(created|modified)\b|<dcterms:W3CDTF>")
"""A line of the model history with a date, which changes with every build."""


def reference_definition() -> PBPKDefinition:
    """The definition of the stored reference: every option of the library unset.

    `end` in the full form with kinetic erythrocytes (`Kbc`, `K_int_pls`,
    `K_cell_pls`), `mrg` in the merged tissue with a Kp override and
    erythrocytes in equilibrium (`Kp_<s>`, `Kp_ad_<s>`), with an IV and an
    oral route; no baseline, no parameters of the definition.
    """
    return PBPKDefinition(
        sid="baseline_reference",
        substances=[
            Substance(
                "end",
                name="endogenous substance",
                mr=150.0,
                binding=Binding(fu=0.5, protein="albumin"),
                erythrocytes=Erythrocytes("kinetic", kbc=2.0, k_exchange=1.0),
                distribution=Kp(2.0, mode=Mode.PERMEABILITY, ftissue=1.0),
            ),
            Substance(
                "mrg",
                name="merged substance",
                mr=200.0,
                erythrocytes=Erythrocytes("equilibrium", kbc=1.5),
                distribution=Kp(2.0, overrides={Organ.ADIPOSE: 5.0}, ftissue=1.0),
            ),
        ],
        routes=[Route.iv("mrg"), Route.oral("mrg")],
    )


def dateless(text: str) -> list[str]:
    """The lines of an SBML text without the dates of the model history."""
    return [line for line in text.splitlines() if not _DATE.search(line)]


def assert_reference_sbml(result: BuildResult) -> None:
    """Assert that a build is the stored reference of f008d56, line by line.

    The order of the elements counts. The references are
    `tests/pbpk/data/<file>.gz` (`gzip -n`), built from
    `reference_definition()` at f008d56; a change of the library which
    changes the SBML of every model on purpose rebuilds them and says so.
    """
    for path in (result.comp_path, result.flat_path):
        reference = gzip.decompress((REFERENCE / f"{path.name}.gz").read_bytes())
        expected = dateless(reference.decode("utf-8"))
        actual = dateless(path.read_text(encoding="utf-8"))
        if actual != expected:
            diff = list(
                difflib.unified_diff(
                    expected, actual, "f008d56", "now", n=0, lineterm=""
                )
            )
            pytest.fail(f"{path.name} differs:\n" + "\n".join(diff[:40]))
