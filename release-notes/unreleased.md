# Release notes for pkpdsim (unreleased)

The body figure of a model app is built by pkpdsim itself: the drawing, the readouts sidecar and the build of a model are all written here, so the panels of the figure are read back from the build rather than kept as a hand-maintained list in every model package which shows one.

## Features

- **The body figure of a model app is built by pkpdsim** ([#222](https://github.com/matthiaskoenig/pkpdsim/issues/222)). `pkpdsim.webapp.body_figure(svg, sbml, *, simulation, title, credit, substances, knobs, readouts)` turns the drawing of a model and its readouts sidecar into the panels of an app: it restricts the ids of the sidecar to the entities of the flat build (read with libsbml, no compile), puts the knobs of the app on every panel it builds, and labels the curves as the panels of an app label theirs, the compartment alone where a panel names one substance and the short name of the substance beside it where it draws several. It returns the panels by their key, the figure and the keys whose readouts no build has, which stay decorations of the drawing. The three steps were a copy in every model package which showed its figure, which is what a model package may not share (a model folder has to be exportable on its own), and the restriction belongs here in any case: pkpdsim writes both the figure and the model, so it can read the ids back instead of asking a package to keep a hand-maintained list of what its model lacks. The page is `docs/webapp.md`.

## Fixes

- **The lumen key of a model figure reports the drug of every segment** ([#222](https://github.com/matthiaskoenig/pkpdsim/issues/222)). The summed dissolved amount `Alumen_<s>` is a readout of the library which a model requests with `PBPKDefinition.readouts`, so a model which does not request it has no such id and the lumen key of its figure proposed the one id it had: a consumer which checks the sidecar against its build dropped it and clicking the intestine of an oral model reported nothing at all. The key now names the species of every lumen compartment of the segmented intestine (`GU__Clumen_du_<s>`, whose plain id is its amount in `mmole`, labelled in the names of the segments, the transit of the drug through the gut) or the one species of a single lumen, and keeps the summed amount where the model declares it. `ModelContent` carries the `lumen_segments` of the intestine and the `readouts` the model requests, and the new tests check the ids of the lumen against a compiled model rather than against the id scheme, which is what let the gap through. The urine and the feces rates stay proposals of the sidecar, dropped where a model does not declare them, and their amounts are unaffected; `docs/pbpk/figures.md` says which ids a model has to ask for.

## Dependencies

- `sbmlsim>=0.34.0`, whose `panels_from_figure` takes the knobs of the app.
