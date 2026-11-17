# Unreleased changes

## Fixes

- Protect model input data and unfamiliar nonempty directories from CLI output replacement. Concurrent writers use a per-result lock, interrupted publication restores the previous result, and failed simulation retries retain the earlier report.
- Keep the model clock and delayed events continuous across segmented timecourses, and reset initial event delays before repeating a simulation.
- Align uncertainty draws by time, preserve discontinuities, skip failed draws, and return an empty band when every draw fails.
- Treat a study's collection date as unknown when a candidate data copy lacks provenance, so temporal validation cannot claim that evidence was collected after fitting.

## Performance

- Reuse checksum-verified, validated PBPK model pairs when the composed SBML, software versions, and validation settings match. Unchanged outputs retain their modification times. Set `PKPDSIM_BUILD_CACHE` to a shared directory or to `off`; models with external references always revalidate.
- Read kidney validation values from the initialized model instead of integrating a one-minute simulation for each condition. Keep recent values in a bounded cache.
