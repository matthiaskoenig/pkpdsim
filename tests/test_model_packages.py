"""The model objects of the packages of pkdb_models are created as they are written.

The keywords of three model packages of pkdb_models (`develop` of 2026-11-06,
their `__init__.py`), with the paths of a temporary package: `glimepiride`
(a female build, adopted fits, an `AbsoluteTolerance`, extra actions with their
arguments, known exceptions, the release), `hctz` (two adopted fits, the
credibility declaration) and the `endoglin` of `misc` (no fits, an extra action
of `all`). Only the imports differ, which is what the first change of the
packages after the move into sbmlsim is.
"""

from pathlib import Path

import pytest
from sbmlsim.model.tolerances import AbsoluteTolerance
from sbmlsim.plot import FigureSettings as SbmlsimFigureSettings
from sbmlsim.project import CheckSettings, FigureOutputs, SimulatorSettings
from sbmlsim.release.settings import Author

from pkpdsim.model import ExtraAction, FigureSettings, PKPDModel, PKPDReleaseSettings
from pkpdsim.sex import Sex


@pytest.fixture
def package(tmp_path: Path) -> Path:
    """Get the directory of a model package."""
    return tmp_path / "package"


def test_glimepiride(package: Path) -> None:
    model_base_path = package / "models" / "output" / "models"
    model_path = model_base_path / "glimepiride_body_flat.xml"
    female_model_path = model_base_path / "glimepiride_body_female_flat.xml"
    data_paths = [package / "data", package / "data_extra"]
    absolute_tolerance = AbsoluteTolerance(
        amount=1e-14, concentration=1e-14, other=1e-10
    )
    model = PKPDModel(
        name="glimepiride",
        package="pkdb_models.models.glimepiride",
        title="Glimepiride PBPK model",
        path=package,
        model_path=model_path,
        female_model_path=female_model_path,
        data_paths=tuple(data_paths),
        archive_path=model_base_path.parent / "glimepiride_model.omex",
        factory=".models.factory:create_glimepiride_models",
        experiments=".simulations:EXPERIMENTS",
        fit_definitions=".fitting.fitting:FIT_DEFINITIONS",
        # the contexts of use of the model and their evidence
        credibility=".credibility:CREDIBILITY",
        # the adopted fits, whose values the factory takes (`models/fits.py`)
        fits=("PK_20261031_024416__3c43d",),
        simulator=SimulatorSettings(
            absolute_tolerance=absolute_tolerance, relative_tolerance=1e-10
        ),
        figures=FigureSettings(
            formats=("svg",),
            fig_dpi=300,
            legend_fontsize=10,
            axes_labelsize=16,
            xtick_labelsize=12,
            ytick_labelsize=12,
        ),
        extra_actions=(
            ExtraAction(
                name="list_cyp2c9",
                help="print the CYP2C9 analyses which can be run",
                function=".run_cyp2c9_analyses:list_cyp2c9",
            ),
            ExtraAction(
                name="analyze_cyp2c9",
                help="run the CYP2C9 analyses given by --analyses",
                function=".run_cyp2c9_analyses:analyze_cyp2c9",
                in_all=True,
            ),
        ),
        extra_arguments=".run_cyp2c9_analyses:add_arguments",
        fit_identifiability=True,
        checks=CheckSettings(
            not_registered={
                "GlimepirideCYP2C9Scan": (
                    "population scan which `run_cyp2c9_scan` runs on its own, the "
                    "CYP2C9 analysis `pk` reads its tables"
                )
            },
        ),
        # the release of the model repository
        release=PKPDReleaseSettings(
            title="Physiologically based pharmacokinetic (PBPK) model of glimepiride",
            repository="matthiaskoenig/glimepiride-model",
            doi="10.5281/zenodo.14628700",
            authors=(
                Author(
                    "König, Matthias",
                    "0000-0003-1725-179X",
                    "Humboldt-Universität zu Berlin",
                ),
            ),
            keywords=("model", "glimepiride", "PBPK", "SBML"),
            docker_image="matthiaskoenig/glimepiride",
            funding="Matthias König was supported by the BMFTR.",
            description="The glimepiride PBPK model.",
            first_year=2024,
        ),
    )
    assert model.root == model.path == package
    assert model.model_paths == {Sex.MALE: model_path, Sex.FEMALE: female_model_path}
    assert model.model_path == model_path
    assert model.model_path_of(Sex.FEMALE) == female_model_path
    assert model.data_paths == tuple(data_paths)
    assert model.fits == ".fitting.fitting:FIT_DEFINITIONS"
    assert model.adopted_ids == ("PK_20261031_024416__3c43d",)
    assert model.adopted_fits() == {"PK": "PK_20261031_024416__3c43d"}
    assert model.factory_dir == model_base_path
    assert model.fits_dir == package / "models" / "output" / "fits"
    assert model.simulator.absolute_tolerance == absolute_tolerance
    assert model.figures == FigureOutputs(formats=("svg",))
    assert model.figure_settings == SbmlsimFigureSettings(
        fig_dpi=300,
        legend_fontsize=10,
        axes_labelsize=16,
        xtick_labelsize=12,
        ytick_labelsize=12,
    )
    assert [a.name for a in model.extra_actions] == ["list_cyp2c9", "analyze_cyp2c9"]
    assert model.fit_identifiability is True
    assert model.release is not None
    assert model.release.repository == "matthiaskoenig/glimepiride-model"
    assert set(model.checks.not_registered) == {"GlimepirideCYP2C9Scan"}


def test_hctz(package: Path) -> None:
    model_base_path = package / "models" / "output" / "models"
    data_paths = [package / "data"]
    model = PKPDModel(
        name="hctz",
        package="pkdb_models.models.hctz",
        title="HCTZ PBPK/PD model",
        path=package,
        model_path=model_base_path / "hctz_body_flat.xml",
        data_paths=tuple(data_paths),
        archive_path=model_base_path.parent / "hctz_model.omex",
        factory=".models.factory:create_hctz_models",
        experiments=".simulations:EXPERIMENTS",
        fit_definitions=".fitting.fitting:FIT_DEFINITIONS",
        # the contexts of use of the model and their evidence
        credibility=".credibility:CREDIBILITY",
        simulator=SimulatorSettings(absolute_tolerance=1e-10, relative_tolerance=1e-10),
        figures=FigureSettings(formats=("svg",), fig_dpi=300, legend_fontsize=10),
        fit_identifiability=True,
        # the adopted fits in `models/output/fits`, whose values and uncertainty
        # the factory takes (`models/fits.py`)
        fits=("PK_20261028_003811__eb362", "PD_20261031_014547__bed5d"),
        # the studies left out of the experiments and of the adopted fits
        checks=CheckSettings(
            not_registered={
                "Dussol2005": "left out of the experiments and the fits",
                "Ernst2006": "left out of the experiments and the fits",
                "Januszewicz1959": "left out of the experiments and the fits",
            },
        ),
        release=PKPDReleaseSettings(
            title="PBPK/PD model of hydrochlorothiazide",
            repository="matthiaskoenig/hctz-model",
            doi="10.5281/zenodo.15427418",
            authors=(Author("Schwaiger, Amanda", "0000-0002-8993-8569", "Charité"),),
            keywords=("model", "hctz"),
            docker_image="matthiaskoenig/hctz",
            funding="None.",
            description="The hctz PBPK/PD model.",
            first_year=2024,
        ),
    )
    assert model.female_model_path is None
    assert model.credibility == ".credibility:CREDIBILITY"
    assert model.adopted_fits() == {
        "PK": "PK_20261028_003811__eb362",
        "PD": "PD_20261031_014547__bed5d",
    }
    assert model.archive_path == package / "models" / "output" / "hctz_model.omex"
    assert model.simulator == SimulatorSettings(
        absolute_tolerance=1e-10, relative_tolerance=1e-10
    )
    assert model.figure_settings == SbmlsimFigureSettings(
        fig_dpi=300, legend_fontsize=10
    )
    assert len(model.checks.not_registered) == 3


def test_endoglin(package: Path) -> None:
    model_base_path = package / "models" / "results" / "models"
    model = PKPDModel(
        name="endoglin",
        package="pkdb_models.misc.endoglin",
        title="Endoglin TGF-beta signalling model",
        path=package,
        model_path=model_base_path / "endoglin.xml",
        data_paths=(package / "data",),
        archive_path=model_base_path / "endoglin_model.omex",
        factory=".models.factory:create_endoglin_models",
        experiments=".simulations:EXPERIMENTS",
        extra_actions=(
            ExtraAction(
                name="data",
                help="plot the data of the studies by readout into results/data",
                function=".data.analysis:plot_data",
                in_all=True,
            ),
        ),
        simulator=SimulatorSettings(absolute_tolerance=1e-8, relative_tolerance=1e-8),
        figures=FigureSettings(formats=("svg",), legend_fontsize=8),
    )
    assert model.fits == {}
    assert model.load_fit_definitions() == {}
    assert model.adopted_ids == ()
    assert model.adopted_fits() == {}
    assert model.release is None
    assert model.credibility is None
    assert model.extra_actions[0].in_all is True
    assert model.figure_settings == SbmlsimFigureSettings(legend_fontsize=8)
