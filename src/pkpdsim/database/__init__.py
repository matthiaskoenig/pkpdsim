"""The parameter database of pkpdsim.

The numbers the PBPK library and the models are built from, apart from the
code which uses them. The database has three parts.

- `sources.toml`: the publications every number is taken from, by key, and
  `ModelChoice` for a number the library chooses where the literature gives
  none (`pkpdsim.database.tables.sources`, `cite`).
- `tables/`: the tables of the library, one TOML file per topic (the
  physiology of the reference individuals, the populations, the
  distribution, the intestine, the kidney, the routes, the proteins, the
  impairment). A number is `{ value = 1.8, unit = "kg", source = "ICRP89" }`
  and may state its uncertainty (`sd`, `se`, `cv`, `n`, `ci`), read by
  `pkpdsim.database.tables.entry`.
- `studies/`: the data of the studies as reported, one folder per source
  (`study.toml` with the groups and the measurements, one TSV file of
  records per measurement), read by `pkpdsim.database.studies`. A quantity
  which several studies report is integrated into a consensus with its
  uncertainty by `pkpdsim.database.integration`: the ratio of a group to the
  control group of its study on the log scale, pooled by a random-effects
  meta-analysis.

The database depends on no other module of pkpdsim; the library reads it and
holds no numbers of its own.
"""
