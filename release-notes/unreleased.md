# Release notes for pkpdsim (unreleased)

Two follow-ups of the domain contract of 0.32.0 ([#226](https://github.com/matthiaskoenig/pkpdsim/issues/226)): the two forms of `Domain.clip` disagreed about the value set of a domain, and one of the two public constants named `SEVERITIES` is renamed, which the hash seed of a run had decided the order of.

## Changes

- **`pkpdsim.sensitivity.SEVERITIES` is `SEVERITY_IDS`, and a tuple** ([#231](https://github.com/matthiaskoenig/pkpdsim/issues/231), [#232](https://github.com/matthiaskoenig/pkpdsim/issues/232)). The constant holds the ids of the severities and shared its name with `pkpdsim.pbpk.virtual.SEVERITIES`, which maps each condition to the range a virtual population samples it from, so one of the two had to change and the newer one did. It is a `tuple[str, ...]` where it was a `frozenset[str]`, so the order the ids are reported in is the order they are written in and no longer depends on the hash seed of the run.

## Fixes

- **The interval form of `Domain.clip` respects the value set of a domain** ([#230](https://github.com/matthiaskoenig/pkpdsim/issues/230)). `BOOLEAN.clip(0.85, 1.15)` returned `(0.85, 1)`, a pair whose lower bound the domain itself rejects, and `integers(0, 4).clip(0.3, 3.7)` returned its argument unchanged, while the single value form of the same method has always snapped onto the values of the domain. Each finite bound of an interval is now widened onto a value the domain holds, the lower one down and the upper one up, and the widening is logged at WARNING, because a fractional span of a switch or of a count is a modelling error and not a rounding one. An interval whose bounds are reversed is still returned as it is, so the inverted pair sbmlsim builds for a negative value is not repaired into one a sampler would accept.
