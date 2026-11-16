"""The default web app of a model: pages of its doses, its diseases and its scans."""

from functools import cached_property
from pathlib import Path
from types import SimpleNamespace
from typing import TYPE_CHECKING, Any, ClassVar, cast

import pytest
from sbmlsim.project import Project
from sbmlutils.factory import (
    Compartment,
    Model,
    ModelUnits,
    Parameter,
    RateRule,
    Species,
    UnitDefinition,
    Units,
    create_model,
)

from pkpdsim import webapp

if TYPE_CHECKING:
    from pkpdsim.model import PKPDModel
from pkpdsim.webapp import INTERVALS, PROTOCOL_DAYS, _Studied, default_app


class U(Units):
    """The units of the dosed toy model."""

    mmole = UnitDefinition("mmole")
    min = UnitDefinition("min")
    m2 = UnitDefinition("m2", "meter^2")
    mg = UnitDefinition("mg")
    kg = UnitDefinition("kg")
    mM = UnitDefinition("mM", "mmole/liter")
    per_min = UnitDefinition("per_min", "1/min")
    mg_per_min = UnitDefinition("mg_per_min", "mg/min")
    mmole_per_min = UnitDefinition("mmole_per_min", "mmole/min")


def _dosed_model() -> Model:
    """A model of the doses of the library: an oral and an intravenous dose of `x`."""
    return Model(
        sid="dosed",
        name="dosed toy model",
        units=U,
        model_units=ModelUnits(
            time=U.min,
            extent=U.mmole,
            substance=U.mmole,
            length=U.meter,
            area=U.m2,
            volume=U.liter,
        ),
        compartments=[
            Compartment("Vve", 1.0, name="venous blood", unit=U.liter),
            Compartment("Vurine", 1.0, name="urine", unit=U.liter),
        ],
        species=[
            Species(
                "Cve_x",
                initialConcentration=0.0,
                compartment="Vve",
                substanceUnit=U.mmole,
                hasOnlySubstanceUnits=False,
                name="xanthol (venous blood)",
            ),
            Species(
                "Aurine_x",
                initialAmount=0.0,
                compartment="Vurine",
                substanceUnit=U.mmole,
                hasOnlySubstanceUnits=True,
                name="xanthol (urine)",
            ),
        ],
        parameters=[
            Parameter(
                "PODOSE_x", 0.0, unit=U.mg, constant=False, name="oral dose x [mg]"
            ),
            Parameter(
                "IVDOSE_x",
                0.0,
                unit=U.mg,
                constant=False,
                name="intravenous dose x [mg]",
            ),
            Parameter("ka", 0.05, unit=U.per_min),
            Parameter("kr", 0.01, unit=U.per_min),
            Parameter("f_renal_function", 1.0, unit=U.dimensionless),
            Parameter("BW", 75.0, unit=U.kg, name="body weight"),
            Parameter("Mr", 100.0, unit=U.dimensionless),
        ],
        rate_rules=[
            RateRule("PODOSE_x", "-ka * PODOSE_x", unit=U.mg_per_min),
            RateRule("IVDOSE_x", "-ka * 10 * IVDOSE_x", unit=U.mg_per_min),
            RateRule(
                "Cve_x",
                "(ka * PODOSE_x + ka * 10 * IVDOSE_x) / Mr / Vve - f_renal_function * kr * Cve_x",
                unit=U.mmole_per_min,
            ),
            RateRule(
                "Aurine_x", "f_renal_function * kr * Cve_x * Vve", unit=U.mmole_per_min
            ),
        ],
    )


@pytest.fixture(scope="module")
def sbml(tmp_path_factory: pytest.TempPathFactory) -> Path:
    directory = tmp_path_factory.mktemp("dosed")
    return create_model(
        model=_dosed_model(),
        filepath=directory / "dosed.xml",
        sbml_level=3,
        sbml_version=2,
    ).sbml_path


class _Base:
    """The base of the experiments of the dosed model: its classes of the renal function."""

    renal_map: ClassVar[dict[str, float]] = {
        "Normal renal function": 1.0,
        "Mild renal impairment": 0.69,
        "Severe renal impairment": 0.19,
    }
    cirrhosis_map: ClassVar[dict[str, float]] = {"Control": 0.0, "Mild cirrhosis": 0.4}


def _model(sbml: Path, variants: dict[str, Path]) -> PKPDModel:
    """A model of the dosed toy, with what `default_app` reads of a `PKPDModel`."""
    return cast(
        "PKPDModel",
        SimpleNamespace(
            name="dosed",
            title="Dosed toy model",
            variants=variants,
            release=None,
            variant_path=lambda variant=None: sbml,
            default_panels=lambda path: {},
        ),
    )


@pytest.fixture
def studied(monkeypatch: pytest.MonkeyPatch) -> _Studied:
    studied = _Studied(
        readouts=("[Cve_x]", "Aurine_x", "[Cve_y]"),
        doses={"PODOSE_x": (25.0, 50.0, 12.5), "IVDOSE_x": (10.0,)},
        ends=(24.0, 48.0, 24.0),
        health=frozenset({"healthy", "renal impairment"}),
        scans=(
            {
                "parameter": "BW",
                "default": 75.0,
                "range": [50, 75, 100],
                "units": "kg",
                "label": "body weight [kg]",
            },
            {
                "parameter": "f_renal_function",
                "default": 1.0,
                "range": [0.1, 1.0],
                "units": "dimensionless",
                "label": "renal",
            },
        ),
        base=_Base,
    )
    monkeypatch.setattr(webapp, "_studied", lambda model, ids: studied)
    return studied


def test_the_pages_of_the_doses_the_diseases_and_the_scans(
    sbml: Path, studied: _Studied
) -> None:
    panels, app = default_app(_model(sbml, {"male": sbml, "female": sbml}), sbml)
    assert app is not None
    assert [page.id for page in app.pages] == ["dose", "protocol", "renal", "scan_BW"]
    assert [page.title for page in app.pages] == [
        "Dose",
        "Dosing protocol",
        "Renal impairment",
        "body weight",
    ]
    assert all(page.size == "small" for page in app.pages)
    # a panel per readout the model has, `[Cve_y]` it has not
    assert [panels[i].curves[0].id for i in app.pages[0].panels] == [
        "[Cve_x]",
        "Aurine_x",
    ]
    plasma = panels["dose_0"]
    assert plasma.title == "xanthol (venous blood)"
    assert plasma.curves[0].unit == "uM"
    assert panels["dose_1"].curves[0].unit == "umol"
    # the studied routes, the most studied dosed by default, its studied doses as stops
    oral, iv = plasma.knobs
    assert (oral.target, oral.label, oral.unit) == (
        "PODOSE_x",
        "oral dose xanthol",
        "mg",
    )
    assert (oral.min, oral.max, oral.default) == (0.0, 50.0, 25.0)
    assert [c.value for c in oral.choices] == [12.5, 25.0, 50.0]
    assert (iv.target, iv.default) == ("IVDOSE_x", 0.0)
    # a single dose to the median end of the studies, 24 h in minutes
    assert plasma.simulation.end == 24 * 60.0
    assert plasma.simulation.x_unit == "hr"
    # the variants of the model, which a reader switches between
    assert [v.id for v in app.variants] == ["male", "female"]


def test_the_dosing_protocol_repeats_every_dose(sbml: Path, studied: _Studied) -> None:
    panels, app = default_app(_model(sbml, {"model": sbml}), sbml)
    assert app is not None and app.variants == ()
    protocol = panels["protocol_0"]
    interval, days = protocol.knobs[-2:]
    assert [(c.label, c.value) for c in interval.choices] == list(INTERVALS)
    assert interval.default == 24.0 and interval.target is None
    assert (days.min, days.max, days.default) == (1.0, PROTOCOL_DAYS, PROTOCOL_DAYS)
    assert protocol.simulation.end == PROTOCOL_DAYS * 24 * 60.0
    assert [
        (d.target, d.dose, d.interval, d.duration) for d in protocol.simulation.dosing
    ] == [
        ("PODOSE_x", "PODOSE_x", "interval", "days"),
        ("IVDOSE_x", "IVDOSE_x", "interval", "days"),
    ]


def test_a_disease_page_has_the_classes_of_the_model(
    sbml: Path, studied: _Studied
) -> None:
    panels, app = default_app(_model(sbml, {"model": sbml}), sbml)
    assert app is not None
    renal = panels["renal_0"].knobs[0]
    assert renal.kind == "choice" and renal.target == "f_renal_function"
    assert [(c.label, c.value) for c in renal.choices] == list(_Base.renal_map.items())
    # the scan of the renal function is the page of the disease, no page of a scan
    assert "scan_f_renal_function" not in {page.id for page in app.pages}
    body = panels["scan_BW_0"].knobs[0]
    assert (body.target, body.min, body.max, body.unit, body.label) == (
        "BW",
        50.0,
        100.0,
        "kg",
        "body weight",
    )


def test_the_service_serves_the_default_app(sbml: Path, studied: _Studied) -> None:
    from sbmlsim.webapp.page import build_app

    panels, app = default_app(_model(sbml, {"model": sbml}), sbml)
    project = Project(
        name="dosed", root=sbml.parent, variants={"model": sbml}, panels=panels, app=app
    )
    bundle = build_app(project)
    answers = bundle.content["answers"][""]
    assert [page["id"] for page in answers["app"]["pages"]] == [
        "dose",
        "protocol",
        "renal",
        "scan_BW",
    ]
    assert (
        answers["panels"]["protocol_0"]["run"]["dosing"][0]["interval_factor"] == 60.0
    )


def test_a_model_without_a_dose_shows_its_readouts(
    sbml: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    studied = _Studied(("[Cve_x]",), {}, (), frozenset(), (), None)
    monkeypatch.setattr(webapp, "_studied", lambda model, ids: studied)
    panels, app = default_app(_model(sbml, {"model": sbml}), sbml)
    assert app is not None
    assert [page.id for page in app.pages] == ["model"]
    assert panels["model_0"].knobs == ()


def test_a_scan_map_which_is_a_property_is_no_scan_of_the_class() -> None:
    class Experiment:
        @cached_property
        def scan_map(self) -> dict[str, Any]:
            return {"bw": {"parameter": "BW"}}

    class Model:
        name = "m"

        def experiment_classes(self) -> list[type]:
            return [Experiment]

    # the experiments cannot be read (no project), so only the scans are collected
    studied = webapp._studied(Model(), set())  # ty: ignore[invalid-argument-type]
    assert studied.scans == ()


def test_the_default_app_of_a_model_of_two_builds() -> None:
    """The toy per sex: its experiments compare A with data and dose nothing."""
    from pkpdsim_sexes import MODEL

    panels, app = MODEL.default_app(MODEL.variant_path(None))
    assert app is not None
    assert [page.id for page in app.pages] == ["model"]
    assert [v.id for v in app.variants] == list(MODEL.variants)
    readouts = [panels[i].curves[0].id for i in app.pages[0].panels]
    assert readouts and all(r in ("A", "[A]") for r in readouts)
