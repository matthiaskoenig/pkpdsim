# Anatomy resources

The files of this folder are derived from the anatomograms of the Expression Atlas (EMBL-EBI), which the Human Protein Atlas uses as well: https://github.com/ebi-gene-expression-group/anatomogram (commit 9fcc37022cce1e2862692a5f5fbfb78572b87e67). The anatomograms are distributed under the Creative Commons Attribution 4.0 International license (CC BY 4.0, https://creativecommons.org/licenses/by/4.0/), the code of that repository under the Apache License 2.0.

Changes: the drawings were converted into simplified polylines, the tissues grouped by the organs of pkpdsim, the landmarks and the body cavities derived, by `scripts/visual/anatomy_assets.py` of pkpdsim. The pig (`pig_male.json`, `pig_female.json`) is a silhouette drawn by pkpdsim with the organs of the human drawings fitted into its trunk. The resources are distributed under CC BY 4.0 as well; a figure made from them names the anatomograms in its legend.
