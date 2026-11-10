"""The provenance of every number of the zonation table in the model it builds.

Every other table-backed parameter of the PBPK and the perfusion library
carries the citation of its source in its notes
(`pkpdsim.sbml.parameters.value_parameter`), and a released zonated liver has
to do the same: a credibility or a release review has to be able to tell
stereology from a model choice without opening the code. These tests walk the
keys of `tables/zonation.toml` and the model which is built from them, so a
value added to the table without provenance fails here and not at a review.
"""

import dataclasses
import html
import re

import pytest
from sbmlutils.factory import Model, Parameter

from pkpdsim.database.tables import cite, load, sources
from pkpdsim.pbpk.substance import Substance
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import CHOICE, Value
from pkpdsim.species import Species
from pkpdsim.zonation.builder import ZonationBuilder
from pkpdsim.zonation.flow import ImposedFlow, PressureFlow
from pkpdsim.zonation.geometry import (
    Layer,
    bile_flow,
    fractions,
    geometry,
    table_value,
)
from pkpdsim.zonation.settings import BileTransport, Geometry, ZonationSettings

from .conftest import liver_flow, mass_g

TRACER = Substance("tr", "tracer", Value(180.0, U.g_per_mole))
"""A substance with a molar mass, so the fenestral path of the barrier is built."""

UNASSIGNED = frozenset({"kupffer", "stellate"})
"""The keys of the table which deliberately reach no parameter.

The Kupffer cells (2.1 %) and the stellate cells (1.4 %) of the stereology get
no compartment of the engine, by decision: their volume stays unassigned so
that the represented volume does not inflate (`docs/zonation.md`). They are in
the table to say what the unassigned 3.5 % of the parenchyma is.
"""

SECTIONS = {Species.RAT: "liver_rat", Species.HUMAN: "liver_human"}
"""The section of the table of every species the engine has a liver of."""


def _keys(species: Species) -> list[str]:
    """Every key of the section of the table of a species, in the order of the file."""
    return list(load("zonation")[SECTIONS[species]])


def _model(species: Species = Species.RAT, **kwargs: object) -> Model:
    """A zonated liver of a species with every path of the engine built.

    The pressure-driven flow, so that the leaky tube writes the sinusoid
    geometry and the fenestrae; the bile, so that the bile flow is written;
    the haematocrit of the table, so that the lumen is split.
    """
    kwargs.setdefault("reference_flow", liver_flow(species))
    kwargs.setdefault("flow", PressureFlow())
    settings = ZonationSettings(
        n_zones=3,
        geometry=Geometry.SINUSOID,
        bile=BileTransport.ADVECTIVE,
        species=species,
        liver_mass=mass_g(1000.0),
        **kwargs,  # ty: ignore[invalid-argument-type]
    )
    return ZonationBuilder(settings, [TRACER], None).build()


def _parameters(model: Model) -> list[Parameter]:
    """The parameters of a model."""
    return [p for p in model.parameters or [] if isinstance(p, Parameter)]


def _text(notes: object) -> str:
    """The plain text of the notes of an element.

    `value_parameter` writes the notes as markdown and sbmlutils converts
    them to the XHTML of the SBML element, so a note of a table is compared
    as text with its tags removed, its entities unescaped and its whitespace
    collapsed.
    """
    plain = re.sub(r"<[^>]+>", " ", str(notes))
    return " ".join(html.unescape(plain).split())


def _row_note(species: Species, key: str) -> str:
    """The note of a row of the table, whitespace collapsed.

    The bile flow is the one row the geometry reads with a reader of its own
    (`bile_flow`, which converts ml/min), every other row with `table_value`.
    """
    value = (
        bile_flow(species)
        if key == "bile_flow"
        else table_value(SECTIONS[species], key)
    )
    return " ".join(str(value.notes).split())


@pytest.mark.parametrize("species", list(SECTIONS))
def test_every_value_of_the_table_reaches_the_model_with_its_source(
    species: Species,
) -> None:
    """Every key of the zonation table produces a parameter which names its source.

    The test walks the keys of the table and not a list of parameter ids, so
    a value added to the table without provenance fails: its note has to
    arrive in the notes of a parameter which carries a `Source:` line, either
    unchanged (the geometry, the fractions, the bile flow) or inside the note
    of a value derived from it (the shares of the lumen, the porosity of a
    zone, the count of the fenestrae), and the source has to be a key of
    `sources.toml`, which `cite()` resolves.
    """
    model = _model(species)
    sourced = [
        _text(p.notes) for p in _parameters(model) if "Source:" in _text(p.notes)
    ]
    assert sourced, "no parameter of the zonated liver carries a source"
    citations = [_text(cite(key)) for key in sources()]
    for notes in sourced:
        assert any(c in notes for c in citations), notes
    keys = _keys(species)
    assert len(keys) > len(UNASSIGNED)
    for key in keys:
        if key in UNASSIGNED:
            continue
        note = _row_note(species, key)
        assert note, f"the row '{key}' of the table has no note"
        assert any(note in text for text in sourced), key


@pytest.mark.parametrize("species", list(SECTIONS))
def test_no_parameter_of_a_table_value_has_empty_notes(species: Species) -> None:
    """The parameters the table feeds carry notes, not an empty element.

    The named parameters of the geometry, the stereology and the flow, by the
    ids of the engine: an unsourced number among them is the defect this
    group of tests exists for, and a missing id is a parameter the engine
    stopped writing.
    """
    model = _model(species)
    by_id = {str(p.sid): p for p in _parameters(model)}
    for pid in (
        "rho_li",
        "f_par",
        "f_plasma",
        "f_rbc",
        "f_endothelium",
        "f_disse",
        "f_cell",
        "f_bile",
        "r_sin",
        "L_sin",
        "t_endo",
        "t_disse",
        "r_fen",
        "N_fen",
        "eps_fen_li",
        "eps_fen_z1",
        "Hct_li",
        "Q_bile_ref",
    ):
        assert pid in by_id, pid
        assert "Source:" in str(by_id[pid].notes), pid


def test_the_imposed_flow_carries_the_provenance_of_the_reference_flow() -> None:
    """`Q_li` of an imposed flow names the physiology its reference comes from.

    The flow of the organ is the reference flow of the caller scaled with the
    liver mass, and the reference flow of the PBPK library is an aggregate of
    the physiology table with a source of its own
    (`pkpdsim.pbpk.physiology.reference_flow`), so the parameter of the flow
    carries the citation of that physiology and not an empty note.
    """
    imposed = _model(flow=ImposedFlow())
    flow = next(p for p in _parameters(imposed) if str(p.sid) == "Q_li")
    notes = _text(flow.notes)
    assert "Source:" in notes
    assert "the reference flow of the whole organ scaled with the mass" in notes
    assert "total liver blood flow of the reference male rat" in notes


def test_a_flow_the_caller_states_is_the_model_choice_it_is() -> None:
    """A flow the caller imposes says that it is a choice and does not look sourced."""
    imposed = _model(flow=ImposedFlow(1.0))
    flow = next(p for p in _parameters(imposed) if str(p.sid) == "Q_li")
    notes = _text(flow.notes)
    assert _text(cite(CHOICE)) in notes
    assert "caller of the engine imposed" in notes


@pytest.mark.parametrize("species", list(SECTIONS))
def test_the_derived_shares_of_the_lumen_name_both_rows(species: Species) -> None:
    """The plasma and the erythrocytes rest on the lumen and the haematocrit.

    Neither row alone backs the split, so the two shares are derived values
    which carry both notes; where the two rows have different sources the
    aggregate is the model choice it is and names both keys.
    """
    section = SECTIONS[species]
    lumen, hct = table_value(section, "lumen"), table_value(section, "haematocrit")
    fr = fractions(species)
    for layer in (Layer.PLASMA, Layer.RBC):
        notes = str(fr[layer].notes)
        assert notes.startswith("derived: ")
        assert str(lumen.notes) in notes
        assert str(hct.notes) in notes
        if lumen.source == hct.source:
            assert fr[layer].source == lumen.source
        else:
            assert fr[layer].source == CHOICE
            assert str(lumen.source) in notes
            assert str(hct.source) in notes


def test_a_haematocrit_the_caller_states_is_the_model_choice_it_is() -> None:
    """A haematocrit a caller passes never borrows the source of the table."""
    fr = fractions(Species.RAT, haematocrit=0.3)
    plasma = fr[Layer.PLASMA]
    assert plasma.source == CHOICE
    assert "supplied by the caller" in str(plasma.notes)


def test_a_number_a_caller_puts_in_the_geometry_is_the_model_choice_it_is() -> None:
    """A plain number in `SinusoidGeometry` is wrapped and does not look sourced.

    A caller may replace a field of the geometry with a number (a perfused
    organ its density, a test its parenchymal fraction); the number is the
    choice of that caller and gets the source `ModelChoice` with a note which
    says so, never the source of the row it replaces.
    """
    geo = dataclasses.replace(geometry(Species.RAT), density=1.2, f_par=0.9)
    assert geo.density.source == CHOICE
    assert geo.f_par.source == CHOICE
    assert "supplied by the caller" in str(geo.density.notes)
    assert geo.density.value == pytest.approx(1.2)
    assert geo.f_par.value == pytest.approx(0.9)
    assert geo.r_sin.source == geometry(Species.RAT).r_sin.source


def test_the_human_fallback_reaches_the_notes_of_the_model() -> None:
    """A human build says in its own notes that its stereology is the rat's.

    The design requires the human entries to fall back to the rat fractions
    "with that stated in the table and in the model notes", so a reader of
    the flat SBML of a human zonated liver sees `fallback: rat` on every
    fraction and not only in the table.
    """
    model = _model(Species.HUMAN)
    by_id = {str(p.sid): _text(p.notes) for p in _parameters(model)}
    for pid in (
        "f_plasma",
        "f_rbc",
        "f_endothelium",
        "f_disse",
        "f_cell",
        "f_bile",
        "f_par",
        "r_sin",
        "L_sin",
        "t_endo",
        "t_disse",
        "r_fen",
    ):
        assert "fallback: rat" in by_id[pid], pid
    assert "fallback: rat" not in by_id["rho_li"], "the human density is no fallback"
    for key in ("density", "haematocrit"):
        assert "fallback: rat" not in _row_note(Species.HUMAN, key), key
    # `Hct_li` is derived from the fraction of the lumen, which is a fallback,
    # so its note carries that of the lumen as well as the human haematocrit


def test_the_rat_fractions_claim_no_fallback() -> None:
    """The rat is the source of the fractions, so no fraction of it is a fallback."""
    model = _model(Species.RAT)
    by_id = {str(p.sid): _text(p.notes) for p in _parameters(model)}
    for pid in ("f_plasma", "f_rbc", "f_disse", "f_cell", "f_bile", "rho_li", "r_sin"):
        assert "fallback: rat" not in by_id[pid], pid
