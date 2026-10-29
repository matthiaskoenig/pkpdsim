"""The extracellular distribution: Kp at a scale of the interstitial shares."""

import dataclasses
from pathlib import Path

import libsbml
import pytest

from pkpdsim.pbpk import (
    Formula,
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


def test_the_share_is_that_of_the_sex_of_the_build() -> None:
    k_int = distribution.k_int_pls_value(FU)
    differs = False
    for organ in Organ:
        fractions = reference_fractions(organ)
        for sex, (f_vas, f_int) in fractions.items():
            share = interstitial_share(organ, FU, sex)
            assert share == f_int * k_int / (1.0 - f_vas)
            assert 0 < share < 1
        assert interstitial_share(organ, FU) == interstitial_share(organ, FU, Sex.MALE)
        differs |= interstitial_share(organ, FU, Sex.FEMALE) != interstitial_share(
            organ, FU, Sex.MALE
        )
    assert differs


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
    assert values[Sex.MALE] != values[Sex.FEMALE]


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


def test_a_substance_keeps_its_fields() -> None:
    """`Substance` with the option is the substance as given."""
    kp = Kp(interstitial=1.0, ftissue=1.0)
    s = Substance("drg", name="drug", mr=300.0, distribution=kp)
    assert s.distribution is kp
    assert kp.value is None
