"""The extracellular distribution: Kp at a scale of the interstitial shares."""

import dataclasses
import inspect
import itertools
import re
import warnings
from pathlib import Path

import libsbml
import pytest

from pkpdsim.pbpk import (
    Binding,
    Formula,
    KCellFloorWarning,
    Kp,
    Mode,
    ModelParameter,
    Organ,
    PBPKDefinition,
    Sex,
    Substance,
    U,
    Value,
    build,
    distribution,
    ids,
    interstitial_share,
    reference_individual,
)
from pkpdsim.pbpk.builder import built_organs
from pkpdsim.pbpk.data import load
from pkpdsim.pbpk.impairment import ckd_unbound_fractions
from pkpdsim.pbpk.physiology import reference_fractions

from .conftest import (
    assert_reference_sbml,
    definition,
    individual,
    reference_definition,
    simple_drug,
)

FU = 0.1  # the fu of the simple drug


def with_kp(kp: Kp) -> PBPKDefinition:
    """Get the definition of the simple drug with another distribution."""
    d = definition()
    drug = dataclasses.replace(simple_drug(), distribution=kp)
    return dataclasses.replace(d, substances=[drug, *d.substances[1:]])


def parameters(path: Path) -> dict[str, tuple[float, str]]:
    """Get the value and the notes of the parameters of an SBML file by id."""
    doc = libsbml.readSBMLFromFile(str(path))
    return {
        p.getId(): (p.getValue(), p.getNotesString())
        for p in doc.getModel().getListOfParameters()
    }


FEMALE_INTERSTITIUM = 1.25
"""The factor of the female interstitial fractions of `distinct_shares`."""


@pytest.fixture
def distinct_shares(monkeypatch: pytest.MonkeyPatch) -> None:
    """Give the female reference individual 1.25 times the interstitial fractions of the male.

    The shares of the data differ between the sexes only by the rounding of
    the tabulated `f_int` (`Fint` of Schmitt 2008 is one value per tissue), so
    a test of the sex of the build reads fractions which differ by a factor:
    it tests which sex is read, not the rounding.
    """
    original = distribution.reference_fractions

    def fractions(organ: Organ) -> dict[Sex, tuple[float, float]]:
        f_vas, f_int = original(organ)[Sex.MALE]
        return {
            Sex.MALE: (f_vas, f_int),
            Sex.FEMALE: (f_vas, FEMALE_INTERSTITIUM * f_int),
        }

    monkeypatch.setattr(distribution, "reference_fractions", fractions)


@pytest.mark.usefixtures("distinct_shares")
def test_the_share_is_that_of_the_sex_of_the_build() -> None:
    """The share reads the fractions of the sex it is given."""
    k_int = distribution.k_int_pls_value(FU)
    for organ in Organ:
        fractions = distribution.reference_fractions(organ)
        for sex, (f_vas, f_int) in fractions.items():
            share = interstitial_share(organ, FU, sex)
            assert share == f_int * k_int / (1.0 - f_vas)
            assert 0 < share < 1
        assert interstitial_share(organ, FU, Sex.FEMALE) == pytest.approx(
            FEMALE_INTERSTITIUM * interstitial_share(organ, FU, Sex.MALE), rel=1e-14
        )


def test_the_sex_of_the_share_has_no_default() -> None:
    """Every caller names the sex, so a female build never reads the male share (R2)."""
    sex = inspect.signature(interstitial_share).parameters["sex"]
    assert sex.default is inspect.Parameter.empty
    with pytest.raises(TypeError, match="sex"):
        interstitial_share(Organ.MUSCLE, FU)  # ty: ignore[missing-argument]


@pytest.mark.parametrize("organ", list(Organ))
def test_the_shares_of_the_sexes_differ_by_the_rounding_of_f_int(organ: Organ) -> None:
    """The share is `Fint*K_int`, `Fint` of Schmitt 2008 per tissue, not per sex.

    `f_int = Fint*(1 - f_vas)` of each sex is tabulated to 4 significant
    digits, so `f_int/(1 - f_vas)` gives back `Fint` of the notes of the data
    to that rounding in both sexes; the shares of the sexes differ by it only.
    """
    data = load("physiology_human")
    k_int = distribution.k_int_pls_value(FU)
    for sex in Sex:
        notes = data[sex.value]["organs"][organ]["f_int"]["notes"]
        found = re.search(r"Fint = ([0-9.]+)", notes)
        assert found is not None
        f_int_tissue = float(found.group(1))
        f_vas, f_int = reference_fractions(organ)[sex]
        assert f_int / (1.0 - f_vas) == pytest.approx(f_int_tissue, rel=5e-4)
        assert interstitial_share(organ, FU, sex) == pytest.approx(
            f_int_tissue * k_int, rel=5e-4
        )


@pytest.mark.usefixtures("distinct_shares")
def test_the_female_build_has_the_female_share(tmp_path: Path) -> None:
    """Each build gives every organ the share of its own sex (D18)."""
    kp = Kp(interstitial=1.0, ftissue=1.0)
    values: dict[Sex, dict[Organ, float]] = {}
    for sex in Sex:
        d = dataclasses.replace(with_kp(kp), sex=sex)
        r = individual(build(d, tmp_path / sex.value).flat_path)
        values[sex] = {
            organ: r[ids.kp_organ(organ, "drg")] for organ in built_organs(d)
        }
        for organ, value in values[sex].items():
            assert value == pytest.approx(interstitial_share(organ, FU, sex), rel=1e-15)
    for organ, value in values[Sex.FEMALE].items():
        assert value == pytest.approx(
            FEMALE_INTERSTITIUM * values[Sex.MALE][organ], rel=1e-12
        )


def test_a_value_with_the_shares_is_refused() -> None:
    """With `interstitial` no organ reads `Kp_<s>`: a value is refused."""
    kp = Kp(2.0, interstitial=1.0, ftissue=1.0)
    with pytest.raises(ValueError, match="drg: Kp value and interstitial exclude"):
        with_kp(kp).validate()


def test_every_organ_has_its_share_times_the_scale(tmp_path: Path) -> None:
    """The organs without an override have `fKp_<s> * share`, no `Kp_<s>`."""
    kp = Kp(interstitial=1.0, overrides={Organ.LIVER: 2.0}, ftissue=1.0)
    result = build(with_kp(kp), tmp_path)
    flat = parameters(result.flat_path)
    assert ids.kp("drg") not in flat
    assert flat[ids.fkp("drg")][0] == 1.0
    r = individual(result.flat_path)
    organs = built_organs(with_kp(kp))
    assert Organ.LARGE_INTESTINE not in organs  # the segmented intestine
    assert ids.kp_organ(Organ.LARGE_INTESTINE, "drg") not in flat
    assert r[ids.kp_organ(Organ.LIVER, "drg")] == 2.0
    for organ in organs:
        if organ is not Organ.LIVER:
            pid = ids.kp_organ(organ, "drg")
            share = interstitial_share(organ, FU, Sex.MALE)
            assert r[pid] == pytest.approx(share, rel=1e-15)


def test_the_scale_is_a_parameter_a_fit_changes(tmp_path: Path) -> None:
    scale = Value(1.5, U.dimensionless, source="fit", notes="the fitted scale")
    kp = Kp(interstitial=scale, overrides={Organ.LIVER: 2.0}, ftissue=1.0)
    result = build(with_kp(kp), tmp_path)
    value, notes = parameters(result.flat_path)[ids.fkp("drg")]
    assert value == 1.5
    assert "Fitted value." in notes
    skin = ids.kp_organ(Organ.SKIN, "drg")
    share = interstitial_share(Organ.SKIN, FU, Sex.MALE)
    assert individual(result.flat_path)[skin] == pytest.approx(1.5 * share)
    changed = individual(result.flat_path, **{ids.fkp("drg"): 3.0})
    assert changed[skin] == pytest.approx(3.0 * share)


def test_the_scale_may_be_a_formula_of_a_definition_parameter(tmp_path: Path) -> None:
    phi = ModelParameter("phi_drg", 0.5, U.dimensionless)
    kp = Kp(
        interstitial=Formula("2 dimensionless * phi_drg", U.dimensionless),
        overrides={Organ.LIVER: 2.0},
        ftissue=1.0,
    )
    result = build(dataclasses.replace(with_kp(kp), parameters=[phi]), tmp_path)
    r = individual(result.flat_path)
    share = interstitial_share(Organ.MUSCLE, FU, Sex.MALE)
    assert r[ids.kp_organ(Organ.MUSCLE, "drg")] == pytest.approx(share)


def test_a_kp_without_a_value_needs_every_organ() -> None:
    """Without a value every entered organ needs an override (or the shares)."""
    entered = frozenset({Organ.MUSCLE, Organ.ADIPOSE})
    complete = Kp(overrides={Organ.MUSCLE: 2.0, Organ.ADIPOSE: 3.0}, organs=entered)
    with_kp(dataclasses.replace(complete, ftissue=1.0)).validate()
    partial = Kp(overrides={Organ.MUSCLE: 2.0}, organs=entered, ftissue=1.0)
    with pytest.raises(ValueError, match=r"drg: Kp without a value.*\['ad'\]"):
        with_kp(partial).validate()
    with pytest.raises(ValueError, match="drg: Kp without a value"):
        with_kp(Kp(ftissue=1.0)).validate()


def test_a_kp_of_overrides_only_builds_no_kp_of_the_substance(tmp_path: Path) -> None:
    organs = [o for o in built_organs(definition()) if o is not Organ.LIVER]
    kp = Kp(
        overrides=dict.fromkeys(organs, 2.0),
        organs=frozenset(organs),
        ftissue=1.0,
    )
    flat = parameters(build(with_kp(kp), tmp_path).flat_path)
    assert ids.kp("drg") not in flat
    assert flat[ids.kp_organ(Organ.SKIN, "drg")][0] == 2.0


def test_the_shares_need_a_number_for_the_fu() -> None:
    drug = simple_drug()
    assert drug.binding is not None
    binding = dataclasses.replace(
        drug.binding, fu=Formula("0.1 dimensionless", U.dimensionless)
    )
    d = definition()
    kp = Kp(interstitial=1.0, ftissue=1.0)
    other = dataclasses.replace(drug, binding=binding, distribution=kp)
    with pytest.raises(ValueError, match="interstitial needs the fu as a number"):
        dataclasses.replace(d, substances=[other, *d.substances[1:]]).validate()


def test_a_scale_must_be_positive() -> None:
    with pytest.raises(ValueError, match="Kp interstitial must be positive"):
        with_kp(Kp(interstitial=0.0, ftissue=1.0)).validate()


def test_the_full_form_at_the_share_is_refused() -> None:
    """A full-form organ needs a Kp above its share: the scale 1 is refused."""
    kp = Kp(interstitial=1.0, mode=Mode.PERMEABILITY, ftissue=1.0)
    drug = dataclasses.replace(simple_drug(), distribution=kp)
    with pytest.raises(ValueError, match=r"drg.*'sk'.*interstitial share"):
        distribution.check_kp(drug, Organ.SKIN, Sex.MALE)
    above = dataclasses.replace(kp, interstitial=1.01)
    distribution.check_kp(
        dataclasses.replace(drug, distribution=above), Organ.SKIN, Sex.MALE
    )


@pytest.mark.usefixtures("distinct_shares")
def test_the_full_form_is_checked_for_the_sex_of_the_build() -> None:
    """A Kp between the shares of both sexes is valid for the sex with the smaller one."""
    organ = max(
        Organ,
        key=lambda o: abs(
            interstitial_share(o, FU, Sex.MALE) - interstitial_share(o, FU, Sex.FEMALE)
        ),
    )
    low, high = sorted(Sex, key=lambda sex: interstitial_share(organ, FU, sex))
    between = (
        interstitial_share(organ, FU, low) + interstitial_share(organ, FU, high)
    ) / 2
    kp = Kp(between, mode=Mode.PERMEABILITY, ftissue=1.0)
    drug = dataclasses.replace(simple_drug(), distribution=kp)
    distribution.check_kp(drug, organ, low)
    with pytest.raises(ValueError, match=rf"{high.value} reference individual"):
        distribution.check_kp(drug, organ, high)


# --- the floor of the cell partition in chronic kidney disease (ruling R26) -------


#: the fu of a metabolite of torasemide (M1), bound to albumin with the CKD profile
FU_M1 = 0.14


def _floor_drug(
    kp: float, *, ckd: bool = True, accept: dict[Organ, str] | None = None
) -> Substance:
    """A drug of the full form in skin with a Kp just above its share there."""
    return dataclasses.replace(
        simple_drug(),
        binding=Binding(fu=FU_M1, protein="albumin", ckd=ckd),
        distribution=Kp(
            10.0,
            mode=Mode.PERMEABILITY,
            ftissue=1.0,
            overrides={Organ.SKIN: kp},
            ckd_floor=accept or {},
        ),
    )


def _floor_definition(drug: Substance, sex: Sex = Sex.MALE) -> PBPKDefinition:
    d = definition(sid=f"floor_{sex}", sex=sex)
    return dataclasses.replace(d, substances=[drug, *d.substances[1:]])


def _floor_warnings(record: list[warnings.WarningMessage]) -> list[str]:
    return [str(w.message) for w in record if issubclass(w.category, KCellFloorWarning)]


def test_the_fu_of_the_ckd_stages() -> None:
    """The albumin profile raises fu at every stage, `ckd=False` and no binding none."""
    fu = ckd_unbound_fractions(_floor_drug(1.0))
    assert list(fu) == ["mild", "moderate", "severe", "end stage"]
    # end stage: the albumin 0.905 times its binding function 0.669
    factor = 0.905 * 0.669
    assert fu["end stage"] == pytest.approx(
        1 / (1 + factor * (1 - FU_M1) / FU_M1), rel=1e-12
    )
    assert fu["end stage"] == pytest.approx(0.212, abs=1e-3)
    assert all(a <= b for a, b in itertools.pairwise(fu.values()))
    assert ckd_unbound_fractions(_floor_drug(1.0, ckd=False)) == {}
    assert ckd_unbound_fractions(dataclasses.replace(simple_drug(), binding=None)) == {}


def test_a_kp_below_the_share_of_ckd_is_named(tmp_path: Path) -> None:
    """The build warns: organ, substance, stages, the build's sex and the floor.

    The Kp of skin is 1 % above its share at the reference fu, so the build
    accepts it, and below the share at the fu of the late CKD stages, where
    the rule floors `K_cell_pls` at `K_CELL_FLOOR`.
    """
    for sex in Sex:
        share = interstitial_share(Organ.SKIN, FU_M1, sex)
        drug = _floor_drug(1.01 * share)
        stages = distribution.kp_floor_stages(drug, Organ.SKIN, sex)
        assert [stage for stage, _, _ in stages][-1] == "end stage"
        assert "mild" not in [stage for stage, _, _ in stages]
        with warnings.catch_warnings(record=True) as record:
            warnings.simplefilter("always")
            build(_floor_definition(drug, sex), tmp_path)
        (message,) = _floor_warnings(record)
        assert message.startswith("drg: the Kp ")
        assert "organ 'sk'" in message
        assert f"{sex.value} build" in message
        assert "end stage" in message
        assert f"floored at {distribution.K_CELL_FLOOR}" in message


def test_a_kp_with_room_for_ckd_is_not_named(tmp_path: Path) -> None:
    """A Kp above the share at the end stage, and the binding without CKD: no warning."""
    share = interstitial_share(Organ.SKIN, FU_M1, Sex.MALE)
    for drug in (_floor_drug(2.0 * share), _floor_drug(1.01 * share, ckd=False)):
        assert distribution.kp_floor_stages(drug, Organ.SKIN, Sex.MALE) == []
        with warnings.catch_warnings(record=True) as record:
            warnings.simplefilter("always")
            build(_floor_definition(drug), tmp_path)
        assert _floor_warnings(record) == []


REASON = "the binding of the metabolite in renal impairment is not measured"


def test_an_accepted_floor_is_silent_and_in_the_notes(tmp_path: Path) -> None:
    """`ckd_floor` names the organ with its reason: no warning, the notes say why.

    The acceptance is per organ: the same drug with the floor accepted in
    another organ still warns of skin.
    """
    share = interstitial_share(Organ.SKIN, FU_M1, Sex.MALE)
    accepted = _floor_drug(1.01 * share, accept={Organ.SKIN: REASON})
    # the organs are normalised as those of `ps`
    named = Kp(1.0, ckd_floor={"skin": REASON})  # ty: ignore[invalid-argument-type]
    assert dict(named.ckd_floor) == {Organ.SKIN: REASON}
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        result = build(_floor_definition(accepted), tmp_path / "accepted")
    assert _floor_warnings(record) == []
    notes = parameters(result.flat_path)[f"SK__{ids.k_cell_pls('drg')}"][1]
    assert "The model accepts the floor" in notes
    assert REASON in notes
    assert "end stage (fu 0.2119" in notes
    other = _floor_drug(1.01 * share, accept={Organ.MUSCLE: REASON})
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        result = build(_floor_definition(other), tmp_path / "other")
    (message,) = _floor_warnings(record)
    assert "organ 'sk'" in message
    assert "`ckd_floor`" in message
    sk = parameters(result.flat_path)[f"SK__{ids.k_cell_pls('drg')}"][1]
    assert "accepts the floor" not in sk


def test_an_accepted_floor_needs_a_reason() -> None:
    share = interstitial_share(Organ.SKIN, FU_M1, Sex.MALE)
    with pytest.raises(ValueError, match=r"ckd_floor 'sk' needs the reason"):
        _floor_definition(
            _floor_drug(1.01 * share, accept={Organ.SKIN: " "})
        ).validate()
    with pytest.raises(ValueError, match=r"ckd_floor for organs not entered"):
        d = _floor_definition(_floor_drug(1.01 * share))
        drug = d.substances[0]
        assert isinstance(drug.distribution, Kp)
        entered = dataclasses.replace(
            drug.distribution,
            organs=frozenset({Organ.SKIN}),
            ckd_floor={Organ.LIVER: REASON},
        )
        dataclasses.replace(
            d,
            substances=[
                dataclasses.replace(drug, distribution=entered),
                *d.substances[1:],
            ],
        ).validate()


def test_without_the_option_the_sbml_is_unchanged(tmp_path: Path) -> None:
    """A Kp with a value and no `interstitial` is the SBML of f008d56."""
    d = reference_definition()
    assert all(
        isinstance(s.distribution, Kp) and s.distribution.interstitial is None
        for s in d.substances
    )
    assert_reference_sbml(build(d, tmp_path))


def test_the_reference_individual_is_the_default_individual(tmp_path: Path) -> None:
    """Body weight, height and organ volumes are those of the built model."""
    ref = reference_individual(Sex.MALE)
    assert (ref.bw, ref.height) == (73.0, 176.0)
    r = individual(build(definition(), tmp_path).flat_path)
    assert r[ids.BW] == ref.bw
    assert r[ids.HEIGHT] == ref.height
    for organ in Organ:
        assert r[ids.organ_volume(organ)] == pytest.approx(ref.volumes[organ])
        assert ref.tissue_volumes[organ] < ref.volumes[organ]
    with pytest.raises(TypeError):
        ref.volumes[Organ.LIVER] = 1.0  # ty: ignore[invalid-assignment]


def test_the_reference_individual_has_the_healthy_gfr(tmp_path: Path) -> None:
    """The GFR at renal function 1 is `healthy_gfr` of the build of each sex.

    A model reads it instead of a literal (the GFR of the fluid model of hctz,
    0.125 l/min): `GFR_base` of the default individual, in l/min.
    """
    from pkpdsim.pbpk.experiment import healthy_gfr

    for sex in Sex:
        ref = reference_individual(sex)
        path = build(definition(sid=f"gfr_{sex}", sex=sex), tmp_path).flat_path
        healthy = healthy_gfr(path).m_as("l/min")
        assert ref.gfr == pytest.approx(healthy, rel=1e-12), sex
        assert individual(path)[ids.GFR] == pytest.approx(ref.gfr, rel=1e-12)
    assert reference_individual(Sex.MALE).gfr == pytest.approx(0.125, rel=1e-12)
    assert reference_individual(Sex.FEMALE).gfr < reference_individual(Sex.MALE).gfr


def test_the_reference_individual_has_the_cells_and_the_lumen(tmp_path: Path) -> None:
    """The intracellular volumes and the lumen volumes are those of the built model.

    A model converts a value per litre of the cells of an organ (the
    hepatocytes of gdeobdtpa) or per litre of the lumen (the absorption of
    creatine) with them instead of copies of the numbers.
    """
    ref = reference_individual(Sex.MALE)
    r = individual(build(definition(), tmp_path).flat_path)
    for organ in Organ:
        cells = r[ids.organ_volume(organ)] * (
            1.0 - r[ids.vascular_fraction(organ)] - r[ids.interstitial_fraction(organ)]
        )
        assert ref.intracellular_volumes[organ] == pytest.approx(cells, rel=1e-12)
    assert ref.intracellular_volumes[Organ.LIVER] == pytest.approx(1.3859, abs=1e-4)
    assert tuple(ref.lumen_volumes) == ids.GUT_SEGMENTS
    for seg, volume in ref.lumen_volumes.items():
        assert r[f"GU__{ids.lumen_volume(seg)}"] == pytest.approx(volume, rel=1e-12)
    small = sum(
        ref.lumen_volumes[seg] for seg in ("du", "jej1", "jej2", "ile1", "ile2")
    )
    assert small == pytest.approx(0.10499, abs=1e-5)
    female = reference_individual(Sex.FEMALE)
    assert female.lumen_volumes == ref.lumen_volumes
    assert (
        female.intracellular_volumes[Organ.LIVER]
        != ref.intracellular_volumes[Organ.LIVER]
    )
    with pytest.raises(TypeError):
        ref.lumen_volumes["du"] = 1.0  # ty: ignore[invalid-assignment]


def test_a_substance_keeps_its_fields() -> None:
    """`Substance` with the option is the substance as given."""
    kp = Kp(interstitial=1.0, ftissue=1.0)
    s = Substance("drg", name="drug", mr=300.0, distribution=kp)
    assert s.distribution is kp
    assert kp.value is None
