"""The release of a model repository with pkpdsim, the command `pkpd_release`.

The release is sbmlsim's (`sbmlsim.release`); pkpdsim adds its release hook
(`pkpdsim.release.hook`), which a `PKPDModel` names in its release settings,
and the command `pkpd_release` (`pkpdsim.release.cli`), which builds a release
with the kind `release` of sbmlsim and publishes it.
"""
