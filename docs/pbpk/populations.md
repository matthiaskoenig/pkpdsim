# Virtual populations

A virtual individual is a combination of factors, each of which the library knows: the population (the distribution of the body weight and the height over age, [Populations](physiology.md#populations)), the sex (the build of the model, [Sex](physiology.md#sex)), the age (the [physiology over age](physiology.md#physiology-over-age)) and the impaired organ function (the severities of the [impairment](impairment.md)). `virtual_population` draws them:

- the sex by its share (`female`, 0 for men only, 1 for women only);
- the age uniformly within a range of years or an age class of ICH E11 and ICH E7 (`ages=(20, 50)`, `ages="old_old"`, [Age classes](classes.md#age-classes));
- the height and the logarithm of the body weight from the bivariate normal distribution of the age and the sex of the population (the moments interpolated linearly between the ages of the data, the correlation of height and body weight included where the data give it), truncated at three standard deviations of each;
- a severity as a fixed value or uniformly within a range (`conditions={"CP_severity": (1, 3), "f_renal_function": 0.5}`; `CP_severity` 0 to 3, `HF_severity` 0 to 4, `f_renal_function` 0 to 1.5).

The draws are reproducible by their seed. The changes of an individual (`AGE`, `BW`, `HEIGHT` and the severities) are initial values: set before the initialization, the physiology of the build of its sex follows, the organs, the flows, the GFR, the binding proteins and the activity of the proteins of its age.

```python
import numpy as np

from pkpdsim.pbpk import Sex, virtual_population

vp = virtual_population(
    "white_american_nhanes",
    1000,
    ages="old_old",
    female=0.5,
    conditions={"CP_severity": (0.0, 1.0)},
    seed=1234,
)
women = vp.of_sex(Sex.FEMALE)
print(len(vp), len(women))
print(round(float(np.mean([i.age for i in vp.individuals])), 1))
print(
    round(float(np.median([i.bw for i in women])), 1),
    round(float(np.mean([i.height for i in women])), 1),
)
print(sorted(vp.values(Sex.FEMALE)))
```

```text
1000 498
80.1
67.6 156.8
['AGE', 'BW', 'CP_severity', 'HEIGHT']
```

## Simulating a population

A study of both sexes simulates each sex on the build of its sex: a scan over the individuals of the sex, whose dimension `population_dimension` gives (the covariates in the units of the model), the simulation of the other sex declared in `simulation_sexes` of the experiment.

```python
# not executed
from sbmlsim import Q
from sbmlsim.simulation import Scan, Simulation

from pkpdsim.pbpk import Sex, virtual_population


class ElderlyPopulation(HctzSimulationExperiment):
    """A dose of 25 mg in 200 elderly women and men."""

    simulation_sexes = {"po_women": Sex.FEMALE}
    population = virtual_population("white_american_nhanes", 200, ages=(65, 85))

    def simulations(self) -> dict[str, Scan]:
        dose = Simulation(end=48 * 60, preinit_changes={"PODOSE_hctz": Q(25, "mg")})
        return {
            "po_men": Scan(
                dose, dimensions=[self.population_dimension(self.population, Sex.MALE)]
            ),
            "po_women": Scan(
                dose,
                dimensions=[self.population_dimension(self.population, Sex.FEMALE)],
            ),
        }
```

A condition with a homeostasis of a PD component starts each individual from its drug-free steady state, `self.equilibrated(scan)`. A scan over a population is checked as every simulation: an individual whose filtered plasma water exceeds the plasma flow of its kidney is refused when the experiment is initialized (`kidney_outflow`).

## The populations

| population | ages | correlation of height and body weight |
| --- | --- | --- |
| `white_american_nhanes`, `black_american_nhanes`, `mexican_american_nhanes` | birth to 100 years (the row of 80 years and over holds beyond 80) | yes, of every age |
| `east_asian_tanaka` | 20 to 49 years | none given, drawn independently |
| `european_icrp` | reference individuals only, no distributions | - |

The children of a population have the body size of their age, the physiology of the reference children of ICRP 89 scaled to it; a child of the size of its age is the reference child. The ages beyond the data of a population are refused.

## Limitations

The covariates of a population are the body weight, the height and the age; the further variability of the physiology at a fixed body size (organ masses, blood flows, protein abundances, the sampling of PK-Sim and Simcyp) is not drawn, every individual has the physiology of its covariates. The age is drawn uniformly, not with the age structure of a population. The severities are independent of each other and of the age (a population of cirrhosis drawn with its own age structure gives its ages as the range). The covariate distributions are data of the library and never SBML `distrib` uncertainties: `run_simulation --uncertainty` draws only the uncertainty of the fits.
