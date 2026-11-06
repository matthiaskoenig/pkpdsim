"""The physiologies of the laboratory animals: rat, mouse and pig."""

from pathlib import Path

import libsbml
import pytest

from pkpdsim.pbpk import (
    NephronKidney,
    Organ,
    Pool,
    Sex,
    Species,
    Variants,
    build,
    ids,
    load,
)
from pkpdsim.pbpk.annotation import TAXA
from pkpdsim.pbpk.organs import PORTAL_ORGANS
from pkpdsim.pbpk.physiology import (
    _density,
    _references,
    physiology_objects,
    physiology_table,
    reference_tissue_volumes,
)
from pkpdsim.pbpk.variants.segments import intestine_table, segments

from .conftest import definition, individuals, total_amount

ANIMALS = (Species.RAT, Species.MOUSE, Species.PIG)
"""The laboratory animals of the library."""


def _v(table: dict, key: str) -> float:
    """The value of an entry of a table."""
    return float(table[key]["value"])


def _l_per_min(entry: dict) -> float:
    """A flow of a table in l/min."""
    factor = {"l/min": 1.0, "ml/min": 1e-3}[entry["unit"]]
    return float(entry["value"]) * factor


def _sexes(species: Species) -> list[Sex]:
    """The sexes with a reference physiology of a species."""
    data = load(physiology_table(species))
    return [sex for sex in Sex if sex.value in data]


@pytest.fixture(scope="module")
def builds(tmp_path_factory: pytest.TempPathFactory) -> dict[Species, Path]:
    """The simple definition built for every laboratory animal (male)."""
    return {
        species: build(
            definition(sid=f"simple_{species.value}", species=species),
            tmp_path_factory.mktemp(species.value),
        ).flat_path
        for species in ANIMALS
    }


@pytest.mark.parametrize("species", ANIMALS)
def test_every_reference_animal_is_consistent(species: Species) -> None:
    """Every sex of a species has adipose, a venous pool and flows below the cardiac output."""
    for sex in _sexes(species):
        refs = _references(sex, load(physiology_table(species)), species)
        volumes = reference_tissue_volumes(sex, species=species)
        assert volumes[Organ.ADIPOSE] > 0
        assert refs[ids.pool_fraction(Pool.VENOUS)].value > 0
        fractions = [
            refs[ids.flow_fraction(o)].value for o in Organ if o is not Organ.LUNG
        ]
        assert sum(fractions) == pytest.approx(1.0, abs=1e-12)
        assert refs[ids.ref(ids.GFR)].value < refs[ids.ref(ids.CO)].value
    assert Sex.MALE in _sexes(species)


@pytest.mark.parametrize("species", ANIMALS)
def test_the_reference_animal_is_the_default_individual(
    builds: dict[Species, Path], species: Species
) -> None:
    """The build has the body weight, cardiac output, GFR, blood and hematocrit of its table."""
    p = load(physiology_table(species))["male"]
    r = individuals(builds[species])()
    assert r[ids.BW] == _v(p, "bodyweight")
    assert r[ids.CO] == pytest.approx(_l_per_min(p["cardiac_output"]), rel=1e-9)
    assert r[ids.GFR] == pytest.approx(_l_per_min(p["gfr"]), rel=1e-9)
    assert r[ids.VBLOOD] == pytest.approx(_v(p["blood"], "volume"), rel=1e-9)
    assert r[ids.HCT] == _v(p, "hematocrit")
    # tissues and blood are the body weight, organ blood and pools the blood
    vblood = r[ids.VBLOOD]
    organ_blood = {o: _v(p["organs"][o], "blood") * vblood for o in Organ}
    tissue = sum(
        (r[ids.organ_volume(o)] - organ_blood[o]) * _density(p["organs"][o]).value
        for o in Organ
    )
    assert tissue + vblood * _v(p["blood"], "density") == pytest.approx(
        r[ids.BW], rel=1e-9
    )
    pools = sum(r[ids.pool_volume(pool)] for pool in Pool)
    assert sum(organ_blood.values()) + pools == pytest.approx(vblood, rel=1e-9)
    # the organs in parallel receive the cardiac output
    total = sum(r[ids.organ_flow(o)] for o in Organ if o is not Organ.LUNG)
    assert total == pytest.approx(r[ids.CO], rel=1e-9)
    portal = sum(r[ids.organ_flow(o)] for o in PORTAL_ORGANS)
    assert r[ids.QPO] == pytest.approx(portal, rel=1e-12)


@pytest.mark.parametrize("species", ANIMALS)
def test_an_animal_scales_with_its_body_weight(
    builds: dict[Species, Path], species: Species
) -> None:
    """Organs and blood in proportion to the body weight, cardiac output and GFR allometric."""
    at = individuals(builds[species])
    reference = at()
    bw = reference[ids.BW]
    values = {
        k: reference[k]
        for k in (ids.CO, ids.GFR, ids.VBLOOD, ids.organ_volume(Organ.LIVER))
    }
    b_co, b_gfr = reference[ids.CO_EXPONENT], reference[ids.GFR_EXPONENT]
    heavier = at(BW=2 * bw)
    assert heavier[ids.CO] / values[ids.CO] == pytest.approx(2**b_co, rel=1e-12)
    assert heavier[ids.GFR] / values[ids.GFR] == pytest.approx(2**b_gfr, rel=1e-12)
    assert heavier[ids.VBLOOD] / values[ids.VBLOOD] == pytest.approx(2.0, rel=1e-12)
    liver = ids.organ_volume(Organ.LIVER)
    assert heavier[liver] / values[liver] == pytest.approx(2.0, rel=1e-12)


@pytest.mark.parametrize("species", ANIMALS)
def test_an_animal_has_no_human_covariates(
    builds: dict[Species, Path], species: Species
) -> None:
    """No height, no age, no factors of the age or the sex; the taxon of the species."""
    doc = libsbml.readSBMLFromFile(str(builds[species]))
    model = doc.getModel()
    parameters = {p.getId() for p in model.getListOfParameters()}
    assert ids.BW in parameters
    for pid in (ids.HEIGHT, ids.AGE, ids.BSA, ids.FFM):
        assert pid not in parameters
    assert not [p for p in parameters if p.startswith(("fage_", "fsex_"))]
    cv = model.getCVTerms()
    resources = {
        cv.get(i).getResourceURI(j)
        for i in range(cv.getSize())
        for j in range(cv.get(i).getNumResources())
    }
    taxon = TAXA[species].split("/")[1]
    assert any(r.endswith(taxon) for r in resources)
    assert not any(r.endswith("/9606") for r in resources)


@pytest.mark.parametrize("species", ANIMALS)
def test_an_animal_has_its_intestine(
    builds: dict[Species, Path], species: Species
) -> None:
    """The transit times and volumes of the lumen are those of the intestine of the species."""
    r = individuals(builds[species])()
    table = intestine_table(species)
    for seg in segments(table):
        assert seg.table == table
    stomach = load(table)["stomach"]
    assert r[ids.T_STOMACH] == pytest.approx(_v(stomach, "transit_time"), rel=1e-12)


@pytest.mark.parametrize("species", ANIMALS)
def test_an_intravenous_dose_is_conserved(
    builds: dict[Species, Path], species: Species
) -> None:
    """The drug of an intravenous dose stays in the body and its sinks over a day."""
    mg = _v(load(physiology_table(species))["male"], "bodyweight")  # 1 mg/kg
    r = individuals(builds[species])(**{ids.iv_dose("drg"): mg})
    r.simulate(0, 24 * 60, 50)
    dose = mg / r[ids.mr("drg")]
    total = total_amount(r, "drg") + total_amount(r, "drg_m")
    assert total == pytest.approx(dose, rel=1e-6)


def test_the_options_of_the_human_are_refused() -> None:
    """A human population, a sex without data, the human nephron: errors naming the species."""
    with pytest.raises(ValueError, match="human"):
        physiology_objects(
            Sex.MALE, population="white_american_nhanes", species=Species.RAT
        )
    missing = [(s, sex) for s in ANIMALS for sex in Sex if sex not in _sexes(s)]
    for species, sex in missing:
        with pytest.raises(
            ValueError, match=f"no reference physiology of the sex '{sex.value}'"
        ):
            physiology_objects(sex, species=species)
    nephron = definition(
        species=Species.RAT, variants=Variants(organs={Organ.KIDNEY: NephronKidney()})
    )
    with pytest.raises(ValueError, match="nephron"):
        build(nephron, Path("unused"))


def test_the_species_is_a_choice_of_the_definition() -> None:
    """The species is normalised from its value; combined definitions must agree on it."""
    assert definition(species="rat").species is Species.RAT
    assert definition().species is Species.HUMAN
