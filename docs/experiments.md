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
| `Mr` | the molecular weights of the substances for the conversion of the data, a `NamedTuple` of the model |
| `labels`, `units` | the labels and units of the observables by id |

The base provides what is the same in every model: `models` (the SBML model of
the model object, without changes), `tasks` (a task per simulation),
`default_changes` (none, a fitted parameter is part of the model and must not be
set by a simulation), `data` (the selections for every task),
`load_dataframe(fig_id, study=None, *, na_values=None)` and
`load_dataframes_by_substance(fig_id, study=None)` (the digitized data of the
figure or table `Fig1` of the study, by default that of the experiment, read
from the `data_path` of the experiment, a path or paths) and `figures_by_sid(*figures)`. A model
overrides what it does differently, e.g. `default_changes` for a model whose
simulations set a parameter.

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
their data by these values, so the models describe their studies with the same
words. `pkpdsim.metadata` holds the closed vocabularies as `StrEnum`s, shared by
all models: `Tissue`, `Route`, `Dosing`, `Fasting`, `ApplicationForm`, `DataType`
and `PKPDData`. The open vocabularies (the health of the subjects, the
coadministered substances, the genotypes) are the model's own `StrEnum`s, which
grow with its studies. Their conventional members are `HEALTHY`, `NONE` (nothing
coadministered) and `NR` (not reported), with the same name and value in every
model.

A model derives the metadata of its mappings from `PKPDMappingMetaData`, which
has the fields `tissue`, `route`, `application_form` (default `NR`), `dosing` and
`fasting` (default `NR`), and adds its own fields. All fields are keyword-only,
so the model's fields without a default may follow the defaults:

```python
from dataclasses import dataclass
from enum import StrEnum, unique

from pkpdsim.metadata import Dosing, PKPDMappingMetaData, Route, Tissue


@unique
class Health(StrEnum):
    HEALTHY = "healthy"
    RENAL_IMPAIRMENT = "renal impairment"
    NR = "not reported"


@unique
class Coadministration(StrEnum):
    NONE = "none"
    PROBENECID = "probenecid"


@dataclass(kw_only=True)
class HctzMappingMetaData(PKPDMappingMetaData):
    health: Health
    coadministration: Coadministration = Coadministration.NONE


metadata = HctzMappingMetaData(
    tissue=Tissue.PLASMA,
    route=Route.PO,
    dosing=Dosing.SINGLE,
    health=Health.HEALTHY,
)
print(metadata.to_dict())
```

`to_dict` is what the tables of sbmlsim show: the member names of the enums
(`PO`, not `po`), the string of any other value, in the order of the fields with
those of `PKPDMappingMetaData` first. No code reads the values, which are readable
text.

## Filters

A `fitting/mapping_collections.py` selects the mappings of a fit by what a study
is, not by listing names. `metadata_filter(metadata_type, *, excluded=None,
**allowed)` builds a filter on the fields of the metadata: a mapping passes when
every allowed field has one of its values and no excluded field has one of its
values. `observable_filter(observables)` selects by the id of the observable
(`observable_id`, the id without the task). Both are callables of sbmlsim's
`MappingFilter` and have a readable `repr`. An unknown field, a class which is no
dataclass or is an instance, a filter without clauses, an empty set of allowed
values (it selects nothing) or of excluded values (it excludes nothing) and values
which are no collection raise when the module of the filters is imported. A single
value is an error, not a set of its characters: write `route={Route.IV}`, not
`route=Route.IV`:

```python
from types import SimpleNamespace

from pkpdsim.fit import metadata_filter

filter_control = metadata_filter(
    HctzMappingMetaData,
    route={Route.PO, Route.IV},
    coadministration={Coadministration.NONE},
    excluded={"dosing": {Dosing.MULTIPLE}},
)
print(filter_control)

mapping = SimpleNamespace(metadata=metadata)  # stands for a FitMapping
assert filter_control("Fig1", mapping)

try:
    metadata_filter(HctzMappingMetaData, rout={Route.PO})
except TypeError as err:
    print(err)

try:
    metadata_filter(HctzMappingMetaData, route=Route.PO)
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
    HctzMappingMetaData,
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
| `one_at_a_time` | the reason why the registered experiments run one at a time (the memory a run of all of them needs), `None` to run them together |

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

The file is a json object with the parts `metadata` and `strings`; an unknown part
or a part of the wrong type is an error (exit code 2):

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
  (legend keys, labels); the arrays are untouched.

The renamed mappings are sorted like a dump, so the comparison does not depend on
the order the renames change.
