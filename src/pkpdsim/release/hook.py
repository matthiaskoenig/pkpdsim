"""What pkpdsim adds to the release of a model repository, the release hook.

The release of a model package is sbmlsim's (`sbmlsim.release`, the kind
`release` of `python -m sbmlsim.run` and `python -m sbmlsim.release
publish`); `pkpd_release` runs both (`pkpdsim.release.cli`). A `PKPDModel`
names `HOOK` in its release settings (`ReleaseSettings.hook`, see
`pkpdsim.model`), which adds the content of the pharmacokinetic models:

- the evidence of the extensions of pkpdsim on the page of the credibility
  (`pkpdsim.credibility`): the ratio pairs and the temporal validation of a
  context of use after its acceptance table, the inputs of its decisive
  mappings at the end of its section and the study data of the adopted fits
  with the records of their collection from pkdb_data after the sections of
  sbmlsim, the templates of `pkpdsim/resources/release`; the rule of the ratio
  pairs is the row `(ratio pairs)` of the acceptance table and the rules of
  the other extensions are left out of it, they have tables of their own, as
  pkpdsim 0.19 showed them;
- the study of a file of the data, `<substance>/<study>` for a file of the
  folder `data/<substance>/<study>` of the package, so the release notes name
  the studies as pkpdsim 0.19 named them.

The framework of the release, `pkpdsim`, and its core libraries are fields of
the release settings, which `PKPDModel` sets. The module imports the release
of sbmlsim only for its types.
"""

from __future__ import annotations

from collections.abc import Sequence
from importlib.resources import files
from importlib.resources.abc import Traversable
from pathlib import Path, PurePosixPath
from typing import TYPE_CHECKING, Any, Final, cast, override

from sbmlsim.release.hooks import ReleaseHook

if TYPE_CHECKING:
    from sbmlsim.credibility.assessment import CredibilityAssessment

__all__ = ["DATA_FOLDER", "HOOK", "TEMPLATES", "PKReleaseHook"]

#: the templates of the release which pkpdsim adds to the ones of sbmlsim
TEMPLATES: Traversable = files("pkpdsim").joinpath("resources", "release")

#: the folder of the data of a model package, `data/<substance>/<study>`
DATA_FOLDER: Final = "data"

#: the evidence extension of the ratio pairs and its name in the acceptance table
_RATIOS: Final = "ratios"
_RATIO_PAIRS: Final = "ratio pairs"

#: the sections of a context of use which the page shows, no ratio pairs an
#: empty list and no other section `None`
_SECTIONS: Final = ("ratios", "temporal", "inputs")


class PKReleaseHook(ReleaseHook):
    """The release hook of pkpdsim, see the module."""

    @override
    def templates(self) -> Sequence[Path | Traversable]:
        """Get the templates of pkpdsim, the parts of the page of the credibility.

        Returns:
            `TEMPLATES`, with `docs/credibility-acceptance.md.j2`,
            `docs/credibility-context.md.j2` and
            `docs/credibility-sections.md.j2`.
        """
        return (TEMPLATES,)

    @override
    def credibility(
        self, assessment: CredibilityAssessment, data: dict[str, Any]
    ) -> dict[str, Any]:
        """Add the evidence of the extensions of pkpdsim to the page of the credibility.

        Args:
            assessment: the assessment of the contexts of use of the release.
            data: the data of the page as sbmlsim shows it.

        Returns:
            `data` with the study data of the adopted fits (`studies`), and for
            every context of use the ratio pairs (`ratios`), the temporal
            validation (`temporal`) and the inputs (`inputs`) of its sections,
            and its acceptance rows of the core followed by the ones of the
            ratio pairs, named `ratio pairs`.
        """
        from sbmlsim.credibility.extensions import ACCEPTANCE_EXTENSION

        studies = assessment.extensions.get("studies") or {}
        data["studies"] = studies.get("studies", [])
        contexts = cast("list[dict[str, Any]]", data["contexts"])
        for entry, result in zip(contexts, assessment.contexts, strict=True):
            rows = cast("list[dict[str, Any]]", entry["acceptance"])
            core = [row for row in rows if ACCEPTANCE_EXTENSION not in row]
            ratios = [
                {**row, ACCEPTANCE_EXTENSION: _RATIO_PAIRS}
                for row in rows
                if row.get(ACCEPTANCE_EXTENSION) == _RATIOS
            ]
            entry["acceptance"] = core + ratios
            for name in _SECTIONS:
                section = result.section(name) or {}
                entry[name] = section.get(name, [] if name == _RATIOS else None)
        return data

    @override
    def data_study(self, data_dir: PurePosixPath, path: PurePosixPath) -> str | None:
        """Get the study of a file of the data, `<substance>/<study>`.

        Args:
            data_dir: the data folder of the file relative to the root of the
                package, e.g. `data/hctz`.
            path: the path of the file relative to `data_dir`.

        Returns:
            `<substance>/<study>` for a file below `data/<substance>/<study>`
            of the package, e.g. `hctz/Beermann1976`; `None` for any other
            file, which belongs to no study.
        """
        parts = (data_dir / path).parts
        if len(parts) < 4 or parts[0] != DATA_FOLDER:
            return None
        return f"{parts[1]}/{parts[2]}"


#: the release hook of pkpdsim, which the release settings of a model name
HOOK: Final = PKReleaseHook()
