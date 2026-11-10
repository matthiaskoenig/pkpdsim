# Zonated liver

`pkpdsim.zonation` is the engine of the zonated liver models: a one-dimensional sinusoid of `n` zones from periportal (`z1`) to pericentral (`z<n>`), in which a substance travels down the sinusoid by flow and diffusion and every zone has a cell model of its own. The design behind it is that of pkpdsim#168 (one zonation builder for all zonated liver models); this page is the reader's view of what the engine builds, what is set, and what to take care of.

The engine writes plain SBML and imports neither the PBPK library nor the perfusion library: a caller (a zonated liver of a body model, a perfused organ) hands over its own boundary and cell model. What it does import of pkpdsim is the SBML writing kit (`pkpdsim.sbml`), the parameter database (`pkpdsim.database.tables`), the species (`pkpdsim.species`) and the zone vocabulary (`pkpdsim.zones`) - and itself. The reference flow is the clearest case of the division: the engine reads no physiology, so the caller states the flow of the reference organ, and the snippet of the Peclet number below takes it from the PBPK library (`pkpdsim.pbpk.physiology.reference_flow`) as a caller does.

To build one, see [the `ZonatedLiver` variant](pbpk/variants.md) of the PBPK library (`Geometry.SINUSOID` for the six layers), [machine perfusion](perfusion.md) for a liver in a circuit, and [the build](pbpk/build.md) for what a build writes; the snippet under [One zonated liver](#one-zonated-liver) builds a small one with the engine alone.

## One zonated liver

The engine is one class: `ZonationBuilder(settings, substances, cell_model)` and its `build()`, which returns the SBML model of the chain. A liver of ten zones with gadoxetate, a pericentral porosity gradient and the reference flow of the rat:

```python
from pkpdsim.pbpk.physiology import reference_flow
from pkpdsim.pbpk.substance import Substance
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Value
from pkpdsim.species import Species
from pkpdsim.zonation.builder import ZonationBuilder
from pkpdsim.zonation.flow import ImposedFlow
from pkpdsim.zonation.profiles import Linear
from pkpdsim.zonation.settings import Geometry, ZonationSettings

settings = ZonationSettings(
    n_zones=10,
    geometry=Geometry.SINUSOID,
    flow=ImposedFlow(),  # no value: the reference flow, scaled with the mass
    species=Species.RAT,
    liver_mass=Value(0.01, U.kg),
    reference_flow=reference_flow(Species.RAT),  # the caller states it
    profiles={"porosity": Linear(1.0, 2.0)},  # twice as porous pericentral
)
gadoxetate = Substance("ged", "gadoxetate", Value(682.0, U.g_per_mole))
model = ZonationBuilder(settings, [gadoxetate], None).build()

# six layers per zone, the chain between them, and the fenestral path
compartments = {str(c.sid) for c in model.compartments}
assert {"V_plasma_z1", "V_disse_z1", "V_cell_z10", "V_bile_z10"} <= compartments
reactions = {str(r.sid) for r in model.reactions}
assert {"FLOW_ged_plasma_z1_plasma_z2", "FEN_ged_z1"} <= reactions

# every number of the table reaches the model with the citation of its source
notes = {str(p.sid): str(p.notes) for p in model.parameters}
assert "Source:" in notes["rho_li"]
```

`cell_model=None` is a chain without cells, which is what a tracer needs; a real model hands over a cell model and the ports of [the cell model contract](#the-cell-model-contract). Inside a body model this is done for the caller by the [`ZonatedLiver` variant](pbpk/variants.md), inside a circuit by [machine perfusion](perfusion.md).

## The layers

Every zone has up to six layers, each a compartment `V<layer>_z<k>` which carries every substance:

| Layer | Carries substance | Fraction of the parenchyma |
|---|---|---|
| erythrocytes | only if the haematocrit is above 0 | `Hct` x 10.6 % |
| plasma (the sinusoidal lumen) | yes | `(1 - Hct)` x 10.6 % |
| endothelial cell | yes | 2.8 % |
| space of Disse | yes | 4.9 % |
| hepatocyte (the cell model) | yes | 78.0 % |
| bile canaliculus | yes | 0.4 % |

The fractions are those of the rat liver, Blouin, Bolender and Weibel 1977 (J Cell Biol 72:441), and are entered in the parameter table `tables/zonation.toml` (see [Parameter database](database.md)), the lumen 10.6 % split by the haematocrit. Blouin reports them as fractions of the parenchyma, which is what the engine computes with (`f_par * mass / density`, `f_par` 1.0 until a model enters a parenchymal fraction). The Kupffer cells (2.1 %) and the stellate cells (1.4 %) are in the table and get no compartment, by decision: their volume stays unassigned, so the represented volume does not inflate. The human liver carries the rat fractions with the note `fallback: rat` until a sourced human stereology is entered; its density (1.0 kg/l, the convention of `physiology_human.toml` for an organ to which ICRP 89 gives no specific gravity, a model choice) and its haematocrit (0.43, the male reference of ICRP 89) are human values. The rat values are the library's: the density 1.05 kg/l of `physiology_rat.toml` and the systemic haematocrit 0.46 of He 2017, an upper bound for the lumen, whose haematocrit is lower because of the plasma skimming the engine does not model. There is one haematocrit per species and no sex, so a female human build splits the lumen with the male 0.43 unless the caller passes its own.

A volume is the fraction of the layer times the parenchyma of the liver mass over the density, shared evenly over the zones, so a change of the mass before the simulation sets all volumes. The geometry (`SinusoidGeometry`, per species from the table) gives only areas and lengths: the sinusoid radius `r_sin` 4.4 um and length `L_sin` 300 um, the endothelium of 0.165 um, the space of Disse of 1.0 um (model choices, taken over from the models of pkdb_models), and the number of sinusoids follows from the lumen volume. The axial exchange uses the prism identity, the conductance `D * A / dx` with `A / dx = n**2 * V_k / L_sin**2`, so no radius enters it.

The `WELL_STIRRED` geometry is the same engine with a vascular compartment (lumen, endothelium and space of Disse in one), a cell compartment and the bile per zone, without axial diffusion and without the barrier. Switching is one setting, the experiments do not change.

## The settings

`ZonationSettings` holds what is set; the defaults are the sinusoid of the rat at 10 zones.

| Setting | Meaning |
|---|---|
| `n_zones` | the number of zones, at least 1 (default 10) |
| `geometry` | `Geometry.SINUSOID` (six layers) or `Geometry.WELL_STIRRED` |
| `flow` | `ImposedFlow` (the flow is given) or `PressureFlow` (the flow follows the pressure drop) |
| `bile` | `BileTransport.ADVECTIVE` (default), `DIFFUSIVE`, or `None` for no bile |
| `haematocrit` | the haematocrit of the lumen; `None` takes the species, 0 leaves out the erythrocytes |
| `profiles` | the zonal profile of a process by its name |
| `species`, `liver_mass`, `reference_flow` | the species of the table, the mass of the liver, the blood flow of the reference organ |
| `geometry_constants`, `fractions` | a geometry or fractions replacing the table (a perfused organ states its own) |
| `boundary` | who supplies the ends of the chain: the engine (`OPEN`) or the caller (`EXTERNAL`) |

Both flow modes emit the same flow across the faces of the zones, so nothing downstream depends on which one was chosen. A `PressureFlow` is the leaky tube of the spt model: the pressure falls from `p_pp` (default 1333 Pa, 10 mmHg) to `p_pc` (266 Pa, 2 mmHg), the plasma viscosity is 1.8e-3 Pa s, and the resistance factor of the blood is calibrated per species so that the periportal inflow at the default pressures is the reference flow of the settings. The chain evaluates its flows at the faces of its zones (`Q_li_in` at `x = 0`, `Q_li_<zone>` at the pericentral face of a zone), so the flow it receives is that calibrated inflow at every zone count, which is the reference flow to 1e-12 in both modes. The mean flow over the zones, `Q_li`, is lower than the inflow (0.8315 of it at 10 zones, 0.8298 at 40): it is the mean of the two faces of a zone, averaged over the zones, and sets the velocity and the Peclet number.

**Scanning the flow.** The flow chain is assignment rules, so a change after the initialization (what a simulator applies) reaches all of it: change `Q_li` and the plasma and erythrocyte flow (`Q_plasma_li`, `Q_rbc_li`) follow, and a pressure flow follows a change of `P_pp` or `P_pc` in the flow of every face (`Q_sin_<site>`, `Q_li_<site>`, the carriers of every face, the site `in` or a zone) and in `Q_li`, their mean over the zones. One id is enough for each of these. Two things are separate by design. The bile flow `Q_bile` follows its own reference `Q_bile_ref` and the liver mass and not the blood flow (the bile is not a fraction of the blood flow), so a scan of the blood flow leaves it, and a scan of the bile flow changes `Q_bile_ref`. And under a pressure flow `Q_li` is an output, the mean over the zones: scan the pressures, not `Q_li`. The geometry (the number of sinusoids, the viscosities, the resistances `W_sin`, `w_sin`, `lambda_sin`, the fenestra count, the oncotic pressure) stays an initial assignment of the build: it derives from the mass and the geometry the model was built with, not from a knob a run turns, so a change of the mass is a change before the initialization (as `individual` of the tests does) and a change of the geometry is a new build.

The engine reads no physiology of its own. A default flow needs the `reference_flow` of the caller; for a species of the PBPK library that is `pkpdsim.pbpk.physiology.reference_flow(species)`.

The diffusion coefficients are derived from the molar mass `Mr` and need nothing else for a new substance. The aqueous coefficient at 37 C is `D = D_REF * (Mr / MR_REF) ** -0.4`, anchored on glucose (`MR_REF` 180.16 g/mole, `D_REF` 0.90e-9 m2/s) and checked against sucrose, a Gd chelate, albumin and insulin (`tests/zonation/test_transport.py`, within a factor of 1.5, insulin within 2); the cytosolic coefficient is 0.27 of it (Seksek, Biwersi and Verkman 1997). A `Substance.diffusion` overrides the derived value, and every derived value is written into the model as a parameter whose notes name the correlation, so a fit can replace it. Gadoxetate (682 Da), for example, gets about 0.53e-9 m2/s.

The fenestrae are a path from plasma to the space of Disse in parallel to the endothelial cell. Its permeability is `P_fen = porosity * (D_w / t_endo) * F(lambda)` with the Renkin sieving `F(lambda) = (1 - lambda)**2 * (1 - 2.104 lambda + 2.09 lambda**3 - 0.95 lambda**5)` (Renkin 1954) and `lambda` the ratio of the radius of the molecule to that of the fenestra. The measured gradient is the porosity of Braet and Wisse 2002 (Comp Hepatol 1:1): 6 % periportal rising to 8 % pericentral, a fenestral diameter of 110.7 nm periportal to 104.8 nm pericentral and 9 to 13 fenestrae per um2. The table holds the porosities and, as model choices of the spt model, the radius of 53.5 nm and the density of 10 per um2. Transcellular permeabilities (`P_endo_in`, `P_endo_out`, `P_hep`) are 0, the path is off, until a model sets them.

## The profiles

A profile is a gradient as a shape of the normalized position `x`, 0 periportal and 1 pericentral, so a fit estimates the two parameters of the shape and not one level per zone, and the same declaration gives the same gradient at any number of zones. The shapes are `Constant`, `Linear(pp, pc)`, `Exponential`, `Sigmoid(midpoint, steepness)` (a Hill form, optionally reversed), a list of weights for measured zonation and `Random(seed, cv)`. Every shape is normalized to the mean 1 over the volumes of the zones, so the total amount of a protein does not change with its distribution:

```python
from pkpdsim.zonation.profiles import Linear, Sigmoid, zone_areas, zone_positions

print(zone_positions(4))
levels = Linear(1.0, 4.0).values(4)
print(levels)
assert abs(sum(levels) / 4 - 1.0) < 1e-12
assert abs(sum(Sigmoid(0.5, 6.0).values(10)) / 10 - 1.0) < 1e-12
assert abs(sum(zone_areas(10)) - 1.0) < 1e-12
```

The position of a zone is its centre, `(k - 0.5) / n`. `lobular_weight` is the area of the idealized hexagonal lobule at a position and `zone_areas(n)` the exact area share of every zone, for a comparison with literature values which are weighted by the area.

## The figures

Five visualizations are functions over a simulation result, usable by every zonated model. The result is a `pandas.DataFrame` with a `time` column and a column per selection named by its plain id (`result_frame` converts a roadrunner result); the number of zones and the layers are read from the ids present. A figure is a `matplotlib.figure.Figure` owned by the caller, `matplotlib.pyplot` is not used, and every figure is vector: the heatmap and the colorbar are not rasterized, so an SVG holds no embedded image, which a viewer that drops raster data would show as an empty panel. A position axis names its ends (0 periportal, 1 pericentral), a time carries its unit (`time_unit`, default min) and a concentration its unit (`unit`, default mM) and an optional display `name`. A line axis starts at zero when no value is negative; a colour scale spans the range of the values it draws (the zero is explicit on the ticks of a colorbar, and a scale from zero would show a small gradient as one shade). `space_time` and `lobule_map` take `vmin` and `vmax` to clip a scale which a spike dominates, e.g. the initial condition in one zone; a caption can state the pair. `data_axes(fig)` gives the panels without the colorbars.

1. `zone_profiles`: a row per observable, the zones on the x axis, a line per time, data overlaid; a column per condition with the bands of a fit, below.
2. `space_time`: the zones against the time, the colour the concentration, a heatmap per layer and the outflow curve below.
3. `layer_profile`: one substance in every layer which carries it, at one time, on one scale.
4. `lobule_map`: the idealized lobule hexagon, below; `draw_lobule` draws one into the axes of a caller and `colorbar` is the colorbar every figure here goes through.
5. `tube_mesh` (`pkpdsim.zonation.render`): the points and cells of the sinusoid tube, below.

### The grid of conditions

`zone_profiles` draws more than one condition: in place of a frame it takes a sequence of `ZoneColumn`, one column of the grid each, and draws `len(profiles)` rows by `len(columns)` columns with the x axis shared and the y axis shared per row. The label of a column is the title of its panel in the top row, the colour of a line is the colour of its column where the column has one, and the line style stays the style of the time, so a reader tells the times apart inside a column and the columns apart by colour. A column carries a `Band` per profile and time, the median and, over the draws of an adopted fit, the ends of the band, which are filled in the colour of the line. A row draws times of its own where it has them (`ZoneProfile.times`) instead of the times of the figure, with its own legend right of the last column, and a row of the data names the column it belongs to in the column `column` (a row without one is drawn in every column).

```python
import numpy as np

from pkpdsim.zonation.figures import Band, ZoneColumn, ZoneProfile, zone_profiles

zones = 4
drug = ZoneProfile("Drug", "mM", tuple(f"drg_cell_z{k}" for k in range(1, zones + 1)))
dead = ZoneProfile("Dead cells", "%", tuple(f"alive_z{k}" for k in range(1, zones + 1)))
profiles = [drug, dead]
times = (60.0, 360.0)
columns = []
for label, color, level in (("1 ml/min", "#1f77b4", 1.0), ("3 ml/min", "#d62728", 1.6)):
    values = {}
    for i, profile in enumerate(profiles):
        for time in times:
            median = level * (i + 1) * time / 60.0 * np.linspace(1.0, 0.4, zones)
            values[(profile.label, time)] = Band(median, 0.8 * median, 1.2 * median)
    columns.append(ZoneColumn(label, values, color))
figure = zone_profiles(columns, profiles, times=times)
figure.savefig("zones.svg")
```

The bands are plain arrays, so the figures need no result object of a simulator; `pkpdsim.zonation.results` is the bridge which reads them out of a `sbmlsim.result.scan.ScanResult`, with the draws of an adopted fit as the band. `zone_values` reads the zones at a time, matched and not interpolated: the nearest point of the result, and of several points at the same distance the last, because a change of the simulation (a medium change, a second dose) gives two points at its time and the state after the change is the one to draw. `zone_band` is the same over the draws of an `UncertaintyRun` (the median and the quantiles of the run), `series_band` a time course over the draws, `zone_column` the whole column of a condition (every profile at its own times, read in `source_unit or unit` and transformed into `unit`, so a living fraction becomes the percentage of dead cells) and `band_text` the text a title names the band with:

```python
# not executed: the result and the draws come from a run of the model
from pkpdsim.zonation.figures import zone_profiles
from pkpdsim.zonation.results import band_text, zone_column

columns = [
    zone_column(
        results[group.key],
        profiles,
        times=TIMES,
        run=run,
        time_unit="hr",
        label=group.text,
        color=group.color,
    )
    for group in groups
]
figure = zone_profiles(columns, profiles, times=TIMES, data=measured, time_unit="hr")
figure.suptitle("The zones of the liver" + band_text(run))
```

### The lobule hexagon

`lobule_map(result, quantity, *, time, zones)` draws the value of every zone on the hexagonal lobule: the central vein is the centre, the six portal triads are the corners, and a point has the normalized distance `d_pp / (d_pp + d_pc)` of its nearest triad and the central vein. The rim is zone 1 (periportal), the centre is zone `zones` (pericentral). The bands of the zones are the intervals `[(k - 1) / n, k / n]` of that distance, bounded by its isolines (Apollonius circles of a triad and the centre, cut off by the edges of the hexagon), so the area of a band is exactly `zone_areas(n)`: the periportal and pericentral zones are small and the zones in between, where the isoline reaches the edges, are large. The area of a zone is not the weight `lobular_weight` at its centre: at 10 zones that midpoint rule overstates the zone of the peak by 5.6 %, and the map does not use it. The map is meant for patterns such as necrosis, which are results of the one-dimensional model and not a two-dimensional model. It draws the hepatocytes by default (`layer=Layer.CELL`) and a vector patch per band and wedge.

```python
import numpy as np
import pandas as pd

from pkpdsim.zonation import observers
from pkpdsim.zonation.figures import lobule_map
from pkpdsim.zonation.geometry import Layer

time = np.linspace(0.0, 120.0, 25)
columns = {"time": time}
for k, sid in enumerate(observers.profile_ids("drg", Layer.CELL, 10), start=1):
    columns[sid] = k * time / 60.0
figure = lobule_map(pd.DataFrame(columns), "drg", time=60, zones=10)
figure.savefig("lobule.svg")
```

`draw_lobule(ax, values, mappable=...)` is the renderer behind the map, for a model which composes a grid of lobules of its own (its conditions, its measured data and its legend stay in the model): it draws one compound patch per zone (the twelve wedges of its band as one path, so a grid of 53 panels is 530 artists and not 6360) into the axes its caller made and sets them up for the hexagon, the limits, the equal aspect, no ticks and no spines. A value which is not finite is a zone without data and is drawn in `missing`, a grey (`MISSING`) which the caller may change; the colour is taken here and not from the `bad` of the colormap, whose default is transparent, because a white hole in the hexagon reads as a defect of the geometry. `thirds=True` adds the isolines of 1/3 and 2/3 as dotted rings, the thirds a histology counts its zones in; `anatomy=True` the outline of the hexagon and the three dots of the portal triad of every corner (portal vein, hepatic artery, bile duct); `labels=True` names the vessels, `CV` at the centre and `PT` below the lower right corner, for the first panel of a grid; both are written over a white stroke, so they read on a dark map and on a light one. `points` is the sampling of the isolines (the default of `lobule_bands`), which a grid of small panels sets lower. The centre of this idealization is the central vein itself, so no disc of a vein lumen covers the pericentral bands. `fill_band(ax, x, band, color, alpha=...)` fills a `Band` between its ends the way the profiles of the grid do, for a model which draws a profile of its own beside a map. `colorbar(figure, mappable, ax, label=...)` is the one colorbar of the package: it un-rasterizes the solids of the bar, which matplotlib draws as an embedded image, and marks its axes `COLORBAR`; it passes its keyword arguments on to `Figure.colorbar`, so a grid takes one bar over its panels (`ax=[...]`, `location`, `shrink`, `aspect`).

`lobule_isoline(position)` (`pkpdsim.zonation.render`) is one contour of the lobular position over the same hexagon, the twelve wedges of `lobular_radius` joined into a closed ring around the centre: the isoline of 0 is the hexagon itself (the radius 1 towards a triad, the apothem `cos(pi / 6)` towards the middle of an edge), the isoline of 1 is the centre, the radius towards a triad is `1 - position`, and the area inside the ring is `1 - lobular_area(position)` of the area of the hexagon. A band of `lobule_bands` is two such isolines, so the inner arc of a band is point for point the outer arc of the band inside it and the two leave no seam.

Both take the sampling of their isolines as `points`, the directions over the arc of an Apollonius circle; the part which runs along the edge of the hexagon is straight and needs none. Where the arc ends is `lobular_edge(x)` of the profiles, the direction the engine integrates the area of a zone up to, so the polygons and `zone_areas` cut the lobule in the same place. The default is what a drawing needs: the area of a zone is then within 8.6e-4 of `zone_areas` and the SVG of a panel of ten zones is about 113 kB, where the accuracy a measurement wants (`points=256`, within 3e-6) is 1.4 MB. The areas of the engine come from `lobular_area` and never from these polygons.

```python
import numpy as np
from matplotlib import colormaps
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.figure import Figure

from pkpdsim.zonation.figures import colorbar, draw_lobule

values = np.linspace(0.2, 1.0, 10)
values[3] = np.nan  # a zone without data, drawn in the grey of `MISSING`
scale = ScalarMappable(Normalize(0.2, 1.0), colormaps["viridis"])
figure = Figure(figsize=(6.4, 3.2), layout="constrained")
left, right = figure.subplots(1, 2)
draw_lobule(left, values, mappable=scale, thirds=True, anatomy=True, labels=True)
draw_lobule(right, values[::-1], mappable=scale, anatomy=True)
colorbar(figure, scale, [left, right], label="Drug [mM]", location="right", shrink=0.8)
figure.savefig("lobules.svg")
```

### The tube

`tube_mesh(geo, n_zones)` returns the points (m) and cells of the tube of the geometry, the sinusoid radius and length with one segment per zone from the portal end (x = 0) to the central end, as plain numpy arrays; `cells[k]` lists the points of the two rings of the zone `k + 1`. The engine does not depend on a renderer. A script which draws the tube (colour on the real radii, frames or a video) installs pyvista for its run, as the videos of gdeobdtpa do. VTK has no cell type of a prism with `sides` corners, so the script builds a polyhedron per segment from the two rings (the caps) and the side faces which join ring point `j` to `j + 1` modulo `sides`:

```python
# not executed: pyvista is not a dependency of pkpdsim
import numpy as np
import pyvista as pv

from pkpdsim.species import Species
from pkpdsim.zonation.geometry import geometry
from pkpdsim.zonation.render import tube_mesh

sides = 24
points, cells = tube_mesh(geometry(Species.RAT), 10, sides=sides)
vtk_cells = []
for cell in cells:
    a, b = cell[:sides], cell[sides:]
    faces = [a[::-1], b] + [
        [a[j], a[(j + 1) % sides], b[(j + 1) % sides], b[j]] for j in range(sides)
    ]
    body = [len(faces)] + [i for f in faces for i in (len(f), *f)]
    vtk_cells += [len(body), *body]
types = np.full(len(cells), pv.CellType.POLYHEDRON)
grid = pv.UnstructuredGrid(np.array(vtk_cells), types, points)
grid.cell_data["concentration"] = np.linspace(0.0, 1.0, len(cells))  # a value per zone
```

## The cell model contract

A zone is a `comp` submodel of a standalone cell model, which stays loadable on its own for the couplings to meshes. The engine places it once per zone and replaces its ports:

| Port | Replaced by |
|---|---|
| `Vext` | the space of Disse of the zone |
| `Vli` | the hepatocyte compartment of the zone |
| `Vbi` (optional) | the canaliculus of the zone |
| `[x_ext]` | the Disse species of the substance |
| a bound port (`bind`) | the species of the substance in the layer the caller names (cell, bile, ...) |
| `position` | `(k - 0.5) / n_zones` |
| `temperature` | the temperature of the organ or the circuit |

**A cell which holds the substance.** A cell model which only reads a substance needs `[x_ext]`; a cell model which holds it in a layer of the chain (the Gd-EOB-DTPA cell, whose `ged` in the cell and `ged_bi` in the bile are ports) binds that port with `compose(bind={"ged": ("ged", Layer.CELL), "ged_bi": ("ged", Layer.BILE)})`, or the cell model and the chain would hold two copies of the same substance. The port is replaced by the species of the chain in the layer of the zone, so the id which survives the flattening is the chain's (`ged_cell_z1`). `bind` is not `expose`: `expose` makes a state of the cell addressable by the liver and is wired with `ReplacedBy` onto an element which is no port, while a port cannot carry a `ReplacedBy`. The engine checks `bind` before it composes: the port has to be a port of the cell model (`cell_model=`) the engine does not wire itself, the substance a solute of the chain and the layer one the geometry has (the well-stirred geometry has one vascular site and no plasma, endothelium or Disse layer, the bile layer needs a bile).

**Reaching a cell state from the liver.** The math of the liver (an observer of the cells, a mean over the zones) cannot name `M01__glc_cell`, the id which a flattening gives an element of a submodel: it exists in no model until the flattening has happened, and the composed model would be invalid SBML. The caller exposes the states of the cell it observes (`compose(expose=, cell_model=)`): `expose` maps the id in the cell model to a function of the zone which gives the name of the element in the liver. The engine adds that element per zone, a proxy, of the class of the cell model's element with its unit, compartment, flags and initial value (read from the cell model on every build, never from a table), and wires it onto the element of the submodel with `comp:ReplacedBy`. The flattening merges the two and **the id of the proxy survives**: the flat id of an exposed element is the name the caller gave (`atp_z1_cell`), and every element which is not exposed is `M{k:02}__<id>`, as the couplings address it. A quantity which is a function of something the liver already owns (the temperature factor of its temperature) is computed by the liver, not exposed. An element which is a `Port` cannot be exposed (a port and a `ReplacedBy` exclude each other).

The submodel of zone `k` is named `M{k:02}`. The engine exposes the layer species of a zone to the reactions of the cell model, so a transport may span plasma, endothelial cell, Disse, hepatocyte and canaliculus, and custom kinetics are accepted. The bile runs retrograde, against the blood: from `z<n>` to `z1` and then into a duct (`<s>_bile_duct`, an amount which gives nothing back), in both modes. The default `BileTransport.ADVECTIVE` moves the canalicular fluid with the bile flow of the species (a flow of the whole liver which follows the liver mass, the table value of the rat liver, 0.01562 ml/min, the human canalicular 0.27 ml/min), the substance leaves a zone with `Q_bile * c`. `BileTransport.DIFFUSIVE` replaces the flow by diffusion along the canaliculus towards the duct, an absorbing boundary half a zone beyond `z1`, and needs no bile flow; `bile=None` omits the bile. The direction is the same in both, what a caller sets is the mechanism (`tests/zonation/test_bile.py`).

The ids follow one convention: zones are `z<k>`, a profile factor is `p_<process>_z<k>`, a compartment `V<layer>_z<k>`, a reaction `<name>_z<k>`; the observers (`pkpdsim.zonation.observers`) give the profile of every layer, the outflow concentration, the amount per layer and the mass balance residual `balance_<s>`.

## One zone

`n_zones=1` is the unzonated organ: one zone, the position 0.5, and every profile is a constant 1 whatever its shape. There is no axial reaction, the chain has one tank (a Peclet number of 2), and `lobule_map` draws one colour. The same code path builds it, so it is the check that a zonated model reduces to its well-mixed form (`tests/zonation/test_hexagon.py` pins the Peclet number, the profiles and the colour).

## Velocity

The velocity of the blood in the sinusoid is the reference flow of the caller over the lumen of the table geometry: 108 um/s for the rat and 47 um/s for the human (default settings, the reference flow of the PBPK library). Both are outputs of the engine and not measurements. The intravital literature gives roughly 100 to 500 um/s for the rat and 100 to 300 um/s for the human, so the rat value is at the low end of its range and the human value, 47 um/s, lies below its range of 100 to 300 um/s, by more than a factor of two against its lower end. The engine states this gap and leaves it. The cause is the pairing of the reference flow with the geometry of the table (the sinusoid count follows from the lumen volume, the radius and the length are model choices); the geometry is not adjusted to hit a velocity, and the question is open (the comment of `ZonationBuilder.velocity`).

## Numerical dispersion

The convection is a first order upwind scheme, so a chain of `n` stirred zones disperses like `n` tanks in series. Two Peclet numbers are to be told apart:

- the Peclet number of the discretization alone, `v * L / (v * dx / 2) = 2 n`, with no physical diffusion: 20 at the default 10 zones. This is the value the notes of a built model record ("a Peclet number of 20 without physical diffusion") and what the design's rule of thumb of about `2 n` says;
- the Peclet number with the physical diffusion of a solute, `v * L / (d_eff + v * dx / 2)`, which `effective_peclet(settings, d_eff)` returns: 12.87 at 10 zones for `d_eff` of glucose (0.9e-9 m2/s), lower than 20 because the diffusion adds to the dispersion.

This is not small. At 10 zones, for the rat with the velocity above, the numerical dispersion `v * dx / 2` is 1.63e-9 m2/s, the same order as the aqueous diffusion of a small solute (0.9e-9 m2/s for glucose, the reference of the derivation above). The zone count therefore shares the setting of the axial gradient with the physics, and a finer chain is what separates them: at 40 zones the numerical dispersion is 4.06e-10 m2/s, a Peclet number of 25 with glucose, 80 without physical diffusion. The build warns when the dispersion dominates a diffusion by more than `DISPERSION_FACTOR` (ten times, `pkpdsim.zonation.builder`): the two being of the same order at the defaults, a warning at the bare dispersion fires for every physiological build and says nothing, while a coefficient ten times below it is decoration and the zone count alone sets the gradient.

```python
from pkpdsim.pbpk.physiology import reference_flow
from pkpdsim.species import Species
from pkpdsim.zonation import effective_peclet
from pkpdsim.zonation.settings import ZonationSettings

settings = ZonationSettings(n_zones=10, reference_flow=reference_flow(Species.RAT))
assert effective_peclet(settings, 0.0) == 20.0
assert round(effective_peclet(settings, 0.9e-9), 2) == 12.87

fine = ZonationSettings(n_zones=40, reference_flow=reference_flow(Species.RAT))
assert effective_peclet(fine, 0.0) == 80.0
assert round(effective_peclet(fine, 0.9e-9), 1) == 24.9
```

## Water in the pressure mode

In the pressure flow mode the lumen loses plasma through the fenestrae where its flow falls and takes it back where the flow rises, and the engine accounts for all of it. The chain evaluates its flows at the faces of its zones, so the water of a zone is the difference of its two faces and the first zone filters like any other. The plasma flow across a face is `Q_li_<site> - Hct * Q_li_in`, the filtrate which the space of Disse carries across it is `Q_fil_<zone> = Q_plasma(0) - Q_plasma(x)`, and the erythrocytes keep a constant flow, no erythrocyte passing a fenestra: the three add up to the periportal inflow at every face, measured to the last bit (exactly 0.0) at 1, 10, 40 and 160 zones. The filtrate flows along the space of Disse towards the pericentral end (`FLOW_<s>_disse_z<k>_disse_z<k+1>`) and carries solute with it (`SD_<s>_<zone>`, the solvent drag across the endothelium, sieved by the reflection coefficient of the fenestrae); at the default pressures it peaks in the middle of the sinusoid, where 46.8 % of the plasma which entered the lumen is in the space of Disse, at 10 zones as at 40 (an even zone count has a face at the middle, where the peak is).

What leaves the model is the net filtration of the whole sinusoid, the hepatic lymph (`FLOW_<s>_disse_z<n>_disse_out`, at the flow `Q_fil_z<n>` and the concentration of the last zone), which the observers count in `cum_out_<s>` so that `balance_<s>` stays exact. At the default oncotic pressure, the mean of the two pressures, the sinusoid filters and reabsorbs the same amount and that net is zero: it cancels to a few 1e-16 of the plasma inflow in floating point, which the floor of `FILTRATION_FLOOR` (a billionth of the inflow) reads as none, and the lymph reaction is then not emitted at all. A lower oncotic pressure is the choice which makes the sinusoid filter, and the lymph then takes solute out of the model.

One thing the mode does not do: the volumes are rigid. Every compartment is set once by its initial assignment from the liver mass, so a zone whose lumen has lost a third of its flow does not shrink and the space of Disse which carries the filtrate does not swell - the water is a flow field over fixed volumes. What the chain receives is the whole calibrated inflow at every zone count, the flow at the portal face: while the engine read the flow at the centre of the first zone instead, a pressure-driven build received 0.950 of the reference flow at ten zones (0.747 at one, 0.987 at forty) and the shortfall went nowhere, it was never delivered.

## Limits

The limits the sections above state, with the two which belong nowhere else: the [velocity](#velocity) of the human sinusoid lies below its literature range, the [numerical dispersion](#numerical-dispersion) of the upwind chain is of the order of the aqueous diffusion of a small solute at the default zone count, the [water of the pressure mode](#water-in-the-pressure-mode) is a flow field over rigid volumes, the Kupffer and the stellate cells have no compartment, the transcellular permeabilities are 0 until a model sets them, the human fractions are the rat ones (`fallback: rat`), and the lobule map is a plot of a one-dimensional model and no two-dimensional model.

- **The transporter and enzyme profiles of Gd-EOB-DTPA are constant.** No zonation of OATP1B1, OATP1B3 or MRP2 for gadoxetate was found in the literature, so the gadoxetate model gives those processes a constant profile: its zones differ in their geometry, their porosity and their flow, and not in their transporters. The profiles are the one place to change it - a sourced gradient is entered as the profile of the process and nothing else moves - and until one is entered a zonal difference of the uptake or the biliary excretion of gadoxetate is not something the model has.
- **The flattening and the ODE export of a large zonated model are slow.** The composition writes one `comp` submodel per zone, and the cost of flattening it and of exporting the ODEs grows with the zone count times the size of the cell model; the predecessor spt model is where that was first felt. The merge of the duplicated unit definitions relieves the largest part of it (`pkpdsim.sbml.flatten`: the flat model of a liver of ten zones carried 429 unit definitions against the 39 the model has, one copy per submodel, and keeps 39 now), and what is left is the number of elements itself. A chain of a few dozen zones with a large cell model is a build to run once and keep, not one to put in a loop.

## References

- Blouin, Bolender, Weibel 1977, J Cell Biol 72:441: volume fractions of the rat liver parenchyma.
- Braet, Wisse 2002, Comp Hepatol 1:1: fenestrae, their diameter, frequency and porosity from periportal to pericentral.
- Renkin 1954, J Gen Physiol 38:225: hindered diffusion through pores, the sieving factor.
- Seksek, Biwersi, Verkman 1997, J Cell Biol 138:131: the diffusion in the cytosol.
- Schwen et al. 2015, PLoS ONE 10:e0133653: representative sinusoids, the layers of the sinusoid.
- Berndt et al. 2018, PLoS Comput Biol 14:e1006005: zonated carbohydrate metabolism, the thickness of the space of Disse and the fenestration, pressure-driven convection.
- Goresky 1963 and the multiple-indicator dilution literature: flow-limited and barrier-limited distribution.
