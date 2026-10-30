# Experiments and their metadata

An experiment of a model is a published study as code: it loads the digitized
data of the study, encodes its dosing protocol as simulations and pairs the data
with the observables of the model (the fit mappings). pkpdsim holds what the
experiments of the models share: the base class, the impairment presets, the
vocabulary which describes a mapping and the filters which select the mappings
of a fit by that vocabulary.

## The base of the experiments

The base of the experiments of a model derives from `PKPDSimulationExperiment`
and binds it to the model object with the class attribute `pkpd_model`; the
experiments of the model derive from that base:

```python
# not executed
from typing import ClassVar, NamedTuple

from pkpdsim.experiment import PKPDSimulationExperiment
from sbmlsim import Q
from sbmlsim.units import Quantity

from pkdb_models.models.hctz import MODEL


class MolecularWeights(NamedTuple):
    hctz: Quantity


class HctzSimulationExperiment(PKPDSimulationExperiment):
    pkpd_model = MODEL
    selections = ("time", "[Cve_hctz]", "Aurine_hctz")
    Mr: ClassVar[MolecularWeights] = MolecularWeights(hctz=Q(297.74, "g/mole"))
    labels = {"[Cve_hctz]": "Hydrochlorothiazide"}
    units = {"[Cve_hctz]": "µmole/l"}
```

The class attributes are the conventions of the experiments:

| attribute | meaning |
| --- | --- |
| `pkpd_model` | the model object, required: an experiment of a base without it raises a `TypeError` |
| `selections` | the selections added to every task, e.g. `("time", "[Cve_hctz]")` |
| `Mr` | the molecular weights of the substances for the conversion of the data, a `NamedTuple` of the model; a model of the PBPK library reads them from its definition, `molecular_weights(DEFINITION)` of `pkpdsim.pbpk` (`Mapping[str, Quantity]` by substance, `self.Mr["tor"]`) |
| `labels`, `units` | the labels and units of the observables by id |
| `sex` | the sex of the subjects, whose build of the model the simulations run, `Sex.MALE` by default, see [Builds per sex](#builds-per-sex) |
| `simulation_sexes` | the sex of single simulations which differs from `sex`, by the key of the simulation |

A base which annotates a convention keeps the type of `PKPDSimulationExperiment`,
e.g. `labels: ClassVar[Mapping[str, str]]`, also for the maps of the impaired
organ functions and the labels of the pharmacokinetic parameters: a narrower
type, `ClassVar[dict[str, str]]`, is an invalid override of a mutable attribute
(ty, `invalid-mutable-override`). Only `Mr` has the type of the model, its
`NamedTuple`, which `PKPDSimulationExperiment` does not declare for this reason.

The base provides what is the same in every model: `models` (the build of the
model object which the experiment simulates, without changes), `tasks` (a task
per simulation, on the build of its sex),
`default_changes` (none, a fitted parameter is part of the model and must not be
set by a simulation), `data` (the selections for every task),
`load_dataframe(fig_id, study=None, *, na_values=None)` and
`load_dataframes_by_substance(fig_id, study=None)` (the digitized data of the
figure or table `Fig1` of the study, by default that of the experiment, read
from the `data_path` of the experiment, a path or paths, with the unit of the
geometric means in `gmean_unit`, `add_gmean_unit`; by substance without the
rows which have none), `figures_by_sid(*figures)`, `across()` (the dimension
of the draws of the fitted parameters of a run over them, `None` without,
which a figure of a scan leaves out of its dimensions) and
`add_function(plot, task_id, y, label, color, linestyle="-")` (a function of
the results of a task over the time as a curve of the width of the curves of
the selections, the band of the draws in a run over them; `linestyle=None`
leaves the line style to sbmlsim). A model overrides what it does
differently, e.g. `default_changes` for a model whose simulations set a
parameter.

## Builds per sex

A model with a female build (`PKPDModel.female_model_path`, see
[Models and their commands](models.md#builds-per-sex)) runs a study of women on
it: the experiment sets the sex of its subjects, the class attribute `sex`, and
every simulation runs the build of that sex. A study of a group of men and a
group of women declares the sex of the simulations of the other sex,
`simulation_sexes`, by the key of the simulation:

```python
# not executed
from pkpdsim.sex import Sex


class Smith2001(HctzSimulationExperiment):
    """Hydrochlorothiazide in healthy women."""

    sex = Sex.FEMALE


class Jones2005(HctzSimulationExperiment):
    """Hydrochlorothiazide in a group of men and a group of women."""

    simulation_sexes: ClassVar[Mapping[str, Sex]] = {"po_women": Sex.FEMALE}
```

The build of `sex` is the model `model` of the experiment, as in a model with
one build, so `add_doses` and the units of `Q` need no model key; the build of
another sex of `simulation_sexes` is `model_<sex>` (`model_id(sex)`,
`model_female`), and `sex_of(key)` gives the sex of a simulation. Only the
builds the experiment simulates are loaded. The models of an experiment are
loaded before its simulations are created, so the sexes are declared and not
computed in `simulations`; a key of `simulation_sexes` which is no simulation
of the experiment, and a sex the model has no build of, fail the experiment
with a `ValueError` which names it. A group whose majority is women runs the
female build, a mixed group without a majority of women the male one, which the
study module states. The experiments of a model without a female build are
unchanged.

The experiments of the PBPK library take the healthy GFR, the renal function of
a measured clearance and the renal classes of the build of their sex, see
[Impairment](pbpk/impairment.md#classes-in-experiments). These are class
methods of the experiment, so in a study of both sexes they follow `sex` of the
experiment, not `simulation_sexes`: the simulation of the group of the other sex
passes its sex, else it gets the renal function on the GFR of the wrong build:

```python
# not executed
class Jones2005(HctzSimulationExperiment):
    simulation_sexes: ClassVar[Mapping[str, Sex]] = {"po_women": Sex.FEMALE}

    def simulations(self) -> dict[str, Simulation]:
        men = self.renal_function(45.0)  # the male build, `sex`
        women = self.renal_function(38.0, sex=Sex.FEMALE)
        classes = self.renal_map_of(Sex.FEMALE)  # the renal classes of women
        ...
```

## Keys

The keys of an experiment (its datasets, simulations, tasks, figures and fit mappings) are SIds, which sbmlsim checks when the experiment initializes: letters, digits and `_`, not starting with a digit, and unique across the kinds. Keys are often built from the texts of the data (a dose, a label of a legend), so `pkpdsim.experiment.to_sid(text)` maps a text to an SId. A decimal point between two digits becomes `p`, every other character an SId does not allow becomes `_`, and one character becomes one character. The function maps characters, so a part of a key gives the same part of the SId of the key (`f"task_po_{to_sid(dose)}"`). Whether a whole key is valid (not empty, not starting with a digit) is not checked by `to_sid`, because the text of a part of a key may start with a digit; sbmlsim checks every key at initialize and names all invalid ones:

```python
from pkpdsim.experiment import to_sid

assert to_sid("po 2.5 mg") == "po_2p5_mg"
assert to_sid("0.5mg") == "0p5mg"
assert to_sid("µg/ml") == "_g_ml"
assert to_sid("Cve_hctz") == "Cve_hctz"
assert to_sid("v1.2.3") == "v1p2p3"
assert f"task_po_{to_sid('2.5 mg')}" == "task_po_2p5_mg"
```

Two texts of an experiment can give one SId (`po 2.5` and `po-2.5`), and the second key would silently replace the first entry of a dict. The method `to_sid` of the base gives the same SId and raises a `ValueError` naming both texts, so use it for the keys of an experiment:

```python
# not executed
def simulations(self) -> dict[str, Simulation]:
    return {
        self.to_sid(f"po {dose}"): Simulation(...)  # "po 2.5" and "po-2.5" are an error
        for dose in (2.5, 5.0)
    }
```

## Repeated doses

A dose is given by its dose state, e.g. `PODOSE_hctz`, from which a rate rule of the model dissolves, absorbs or infuses it. The first dose is usually a change before the initialization (`preinit_changes`). A later dose must add to the dose state: a change of the state to the dose, `Change(times, {"PODOSE_hctz": Q(12.5, "mg")})`, discards what is left of the earlier doses, which loses drug at every dose of a slowly dissolving, absorbed or infused form ([pkdb_models#105](https://github.com/matthiaskoenig/pkdb_models/issues/105)). `add_doses` of the base returns the change which adds the dose, a formula of the dose state in its unit in the model:

```python
# not executed
def simulations(self) -> dict[str, Simulation]:
    times = [24 * k for k in range(1, 7)]  # [hr], the doses after the first
    return {
        "hctz12p5_daily": Simulation(
            end=7 * 24,
            preinit_changes={"PODOSE_hctz": Q(12.5, "mg")},
            changes=[self.add_doses(times, {"PODOSE_hctz": Q(12.5, "mg")})],
        )
    }
```

The change is `Change(times, {"PODOSE_hctz": "PODOSE_hctz + 12.5"})` for a model whose dose state is in mg. The models are loaded before `simulations()`, so the unit of the dose state is known there.

## Units of the data

The data of a study can use units which sbmlsim does not define. An international unit (IU) is specific to a substance (the activity of an enzyme, the effect of a coagulation factor), so it has no value without it. `data_units` of an experiment maps the unit as the data write it to the definition of a unit of its own, which the base defines in the registry of the experiment before the data are read and which `load_dataframe` and `load_dataframes_by_substance` rename in the unit columns of the data (`unit` and every column ending in `_unit`). The name carries the model and the substance (`IU_ramipril_renin`, not `IU`): the registry of sbmlsim exists once per process (pytest runs all models in one) and ignores a redefinition, so two models which defined `IU` with different values would silently get the first one. `define_unit` therefore raises for a name which is defined with another value:

```python
from pkpdsim.experiment import define_unit, rename_units
from sbmlsim.units import UnitRegistry
import pandas as pd

ureg = UnitRegistry()
name = define_unit(ureg, "IU_ramipril_renin = 0.6 * microgram")
assert name == "IU_ramipril_renin"
assert define_unit(ureg, "IU_ramipril_renin = 0.6 * microgram") == name

try:
    define_unit(ureg, "IU_ramipril_renin = 1.0 * microgram")
except ValueError as err:
    print(err)

df = pd.DataFrame(
    {"unit": ["IU/ml", "µIU/ml", "IU_per_ml"], "mean_unit": ["IU", "2IU", None]}
)
print(rename_units(df, {"IU": name}))
```

`rename_units` replaces a name where it is a unit, alone or after a metric prefix (`IU/ml` and `µIU/ml`), and leaves it where it is part of another name (`IU_per_ml`, `2IU`). In the experiment:

```python
# not executed
from types import MappingProxyType


class RamiprilSimulationExperiment(PKPDSimulationExperiment):
    pkpd_model = MODEL
    # a read-only mapping, which ruff (RUF012) asks of a class attribute
    data_units = MappingProxyType({"IU": "IU_ramipril_renin = 0.6 * microgram"})
```

`initialize` defines the units before it reads the data, and a loader which is called outside of `initialize` defines them on its first call. `data_unit_names()` gives the names of the definitions (`{"IU": "IU_ramipril_renin"}`, read-only) and defines the units as well, for a table which the model reads itself, e.g. the data of the figure of a scan: `rename_units(df, self.data_unit_names())`. A definition which pint cannot evaluate (an unknown unit in the expression, a syntax error) or which `define_unit` refuses is a `ValueError` naming the unit, which fails the initialization of that experiment only: a run reports it with its traceback, runs the other experiments and raises at the end, and `pkpdsim.testing.initialize_experiments` names it with every other failure.

## Impairment presets

The classes of impaired organ function with the factors of their parameters and
the colors of the figures are the presets of `pkpdsim.impairment`; the base has
the renal function (`RENAL_FUNCTION`, `RENAL_COLORS`), the cirrhosis
(`CIRRHOSIS`, `CIRRHOSIS_COLORS`) and the cardiac function (`CARDIAC`,
`CARDIAC_COLORS`) as its maps (`renal_map`, `renal_colors`, `cirrhosis_map`,
`cirrhosis_colors`, `cardiac_map`, `cardiac_colors`), and the labels and units of
the pharmacokinetic parameters (`PK_LABELS`, `PK_UNITS` of `pkpdsim.pk`) as
`pk_labels` and `pk_units`. A model which classifies otherwise chooses another
preset in its base:

```python
from pkpdsim.impairment import (
    RENAL_COLORS,
    RENAL_COLORS_END_STAGE,
    RENAL_FUNCTION,
    RENAL_FUNCTION_90,
)

# the presets with the end stage class have their own colors
assert set(RENAL_FUNCTION_90) - set(RENAL_COLORS) == {"End stage renal disease"}
assert set(RENAL_FUNCTION_90) == set(RENAL_COLORS_END_STAGE)

# a model without end stage data uses the classes of a preset without that class
renal_map = {
    name: factor
    for name, factor in RENAL_FUNCTION_90.items()
    if name != "End stage renal disease"
}
assert list(renal_map) == list(RENAL_FUNCTION)

# a preset plus a class of the model
renal_map = {**RENAL_FUNCTION, "Hemodialysis": 0.05}
```

```python
# not executed
class HctzSimulationExperiment(PKPDSimulationExperiment):
    pkpd_model = MODEL
    renal_map = renal_map
    renal_colors = {**RENAL_COLORS, "Hemodialysis": "#444444"}
```

## The vocabulary of the metadata

Every fit mapping carries metadata which says what its curve is. The fits select
their data by these values, so all models describe their studies with the same
words: there is one vocabulary for all models. `pkpdsim.metadata` holds it as
`StrEnum`s: `Tissue`, `Route`, `Dosing`, `Fasting`, `ApplicationForm`, `Health`,
`Coadministration`, `Genotype`, `Phenotype`, `AdministrationSite`, `PKPDData`,
`DataType` and `SamplingSite`. Their conventional members are `HEALTHY`, `NONE` (nothing
coadministered) and `NR` (not reported), and a member has one name and one value
in every model. A study which needs a new member adds it to `pkpdsim.metadata`,
not to the model.

A model uses `PKPDMappingMetaData`, which has the fields `tissue`, `route`,
`application_form` (default `NR`), `dosing`, `fasting` (default `NR`), `health`,
`coadministration` (default `NONE`), `genotype`, `phenotype`, `administration_site`,
`data` and `data_type` (default `NR`) and `sampling_site` (stated or `None`; read the
site with `resolved_sampling_site()`. Without a stated site it follows the tissue:
plasma, serum and blood are typically drawn from an arm vein, so they are
`ARM_VEIN`; arterial blood and plasma are `ARTERIAL`; saliva, urine, feces, feces
and urine, bile and no tissue have no blood sample and are `NOT_APPLICABLE`; a blood
site on such a tissue, `NOT_APPLICABLE` on a blood tissue and a site other than
`ARTERIAL` on arterial blood or plasma are a `ValueError`). All fields
are keyword-only, so a field without a default may follow the defaults:

```python
from pkpdsim.metadata import (
    Coadministration,
    Dosing,
    Health,
    PKPDMappingMetaData,
    Route,
    Tissue,
)

metadata = PKPDMappingMetaData(
    tissue=Tissue.PLASMA,
    route=Route.PO,
    dosing=Dosing.SINGLE,
    health=Health.HEALTHY,
)
print(metadata.to_dict())
```

A model extends the class only for a facet no other model has, in its
`fitting/mapping_collections.py`, without redefining a field or repeating a
vocabulary:

```python
from dataclasses import dataclass


@dataclass(kw_only=True)
class LocalMappingMetaData(PKPDMappingMetaData):
    diplotype: str | None = None
```

`to_dict` is what the tables of sbmlsim show: the member names of the enums
(`PO`, not `po`), the string of any other value, in the order of the fields with
those of `PKPDMappingMetaData` first. The values are readable text, and they are
API as much as the names: a model may look a member up by its value, e.g.
`Genotype(f"CYP2C9 {genotype}")` from the genotype in the data of a study, so a
value is not changed lightly either.

## Filters

A `fitting/mapping_collections.py` selects the mappings of a fit by what a study
is, not by listing names. `metadata_filter(metadata_type=PKPDMappingMetaData, /, *,
excluded=None, **allowed)` builds a filter on the fields of the metadata (the class only for a
local extension): a mapping passes when
every allowed field has one of its values and no excluded field has one of its
values. `observable_filter(observables)` selects by the index of the data of the
observable (`observable_id`, i.e. `Data.index` of the `y` of the mapping: the
selection of the model without its brackets and without the task, `Cve_hctz` for
the selection `[Cve_hctz]`; a filter of `[Cve_hctz]` selects nothing). Both are callables of sbmlsim's
`MappingFilter` and have a readable `repr`. `check_sampling_sites(fit_mappings, sites)` checks the `sampling_site` of every mapping against a table of the site of each observable id and returns one problem per mismatch or unknown observable. An unknown field, a class which is no
dataclass or is an instance, a filter without clauses, an empty set of allowed
values (it selects nothing) or of excluded values (it excludes nothing) and values
which are no collection raise when the module of the filters is imported. A single
value is an error, not a set of its characters: write `route={Route.IV}`, not
`route=Route.IV`:

```python
from types import SimpleNamespace

from pkpdsim.fit import metadata_filter

filter_control = metadata_filter(
    route={Route.PO, Route.IV},
    coadministration={Coadministration.NONE},
    excluded={"dosing": {Dosing.MULTIPLE}},
)
print(filter_control)

mapping = SimpleNamespace(metadata=metadata)  # stands for a FitMapping
assert filter_control("Fig1", mapping)

try:
    metadata_filter(rout={Route.PO})
except TypeError as err:
    print(err)

try:
    metadata_filter(route=Route.PO)
except TypeError as err:
    print(err)
```

`select_mappings(model, experiment_classes, filters=(), outliers=(),
validation=())` selects the data of a fit with the paths of the model object
(`path` and `data_paths`) and returns the fit mapping collections by experiment
id; it prints the counts of the mappings per experiment and kind. The filters
select the training data, `outliers` are the keys of mappings which are not
usable and `validation` are keys or a filter of the validation data:

```python
# not executed
from pkpdsim.fit import metadata_filter, observable_filter, select_mappings

from pkdb_models.models.hctz import MODEL
from pkdb_models.models.hctz.experiments.studies import Anderson1961, Barbhaiya1982

PK = observable_filter({"Cve_hctz", "Aurine_hctz"})
CONTROL = metadata_filter(
    route={Route.PO, Route.IV},
    coadministration={Coadministration.NONE},
)

collections = select_mappings(
    MODEL,
    [Anderson1961, Barbhaiya1982],
    filters=[PK, CONTROL],
    outliers=["Barbhaiya1982_Fig2"],
)
```

## Exceptions of the tests

The tests of `pkpdsim.testing` read the known exceptions of a model from
`PKPDModel.checks`, a `CheckSettings`. Every exception carries its reason, and an
exception which no longer applies (a class which is registered now, a module
which is gone) fails the test, so the list cannot go stale:

| field | meaning |
| --- | --- |
| `not_registered` | experiment classes which are not in the group `all` of the registry on purpose, class name to reason |
| `not_importable` | modules below the package which do not import in the test environment (an optional dependency), name relative to the package to reason |
| `one_at_a_time` | the reason why the registered experiments run one at a time (the memory a run of all of them needs), `None` to run them together; an empty reason is an error |

```python
from pkpdsim.model import CheckSettings

checks = CheckSettings(
    not_registered={"Anderson1961Scan": "a scan which takes hours"},
    not_importable={"experiments.zonation.plot_pyvista": "needs pyvista"},
    one_at_a_time="a run of all experiments needs 20 GB",
)
print(checks)
```

```python
# not executed
MODEL = PKPDModel(
    name="hctz",
    ...,
    checks=checks,
)
```

## Renames in a verification

A migration of the metadata renames fields and members, and `pkpdsim-verify
compare` would report every mapping as changed. `compare --renames FILE` applies
the expected renames to the old dump first, so that only the unexpected changes
remain:

```bash
pkpdsim-verify compare base.pkl new.pkl --renames renames.json
```

The file is a json object with the parts `metadata` and `strings` (and, across
versions of sbmlsim, `ids`, `sids`, `figure_fields` and `expected`, see
[Verification across versions of sbmlsim](models.md#verification-across-versions-of-sbmlsim));
an unknown part or a part of the wrong type is an error (exit code 2):

```json
{
  "metadata": {
    "fields": {"tissue_type": "tissue"},
    "members": {"route": {"ORAL": "PO", "INTRAVENOUS": "IV"}},
    "added": {"fasting": "NR"}
  },
  "strings": {"Anderson 1961": "Anderson1961"}
}
```

- `metadata.fields`: the old name of a field to its new name,
- `metadata.members`: per (new) field, the old member name to the new one,
- `metadata.added`: a field the old metadata did not have to its member name, set
  where the old mapping has no such field,
- `strings`: an old text to the new one, in every key and string of the old dump
  (legend keys, labels); the arrays are untouched. A rename which gives two keys of
  a part one key is an error.

The renamed mappings are sorted like a dump, so the comparison does not depend on
the order the renames change.
