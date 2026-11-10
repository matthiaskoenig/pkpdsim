"""The annotations every library of pkpdsim writes.

The initial assignment of a constant (`initial_assignment`, the SBO term of a
mathematical expression and the name of its target) and the NCBI taxon of the
species of a build (`TAXA`, the `BQB.HAS_TAXON` of a model).

The annotation of an element by its kind is the PBPK library's
(`pkpdsim.pbpk.annotation.annotate`, which resolves the PBPKO term of a kind
at its place through the vocabulary of the library): no other library of
pkpdsim annotates by kind, so the table and its resolver stay there.
"""

from sbmlutils.factory import InitialAssignment
from sbmlutils.factory.units import UnitType
from sbmlutils.metadata import SBO

from pkpdsim.species import Species

__all__ = ["TAXA", "initial_assignment"]


def initial_assignment(
    symbol: str, formula: str, unit: UnitType, name: str
) -> InitialAssignment:
    """Get the initial assignment of a constant, with a name and an SBO term.

    Every quantity of the library which depends on parameters only is a
    constant with an initial assignment; the assignment carries the name of
    its target and the SBO term of a mathematical expression, so the factory
    has nothing to hint about.

    Args:
        symbol: the id of the target.
        formula: the formula.
        unit: the unit of the formula.
        name: the name of the target.

    Returns:
        The initial assignment.
    """
    return InitialAssignment(
        symbol, formula, unit, name=name, sboTerm=SBO.MATHEMATICAL_EXPRESSION
    )


TAXA: dict[Species, str] = {
    Species.HUMAN: "taxonomy/9606",  # Homo sapiens
    Species.RAT: "taxonomy/10116",  # Rattus norvegicus
    Species.MOUSE: "taxonomy/10090",  # Mus musculus
    Species.PIG: "taxonomy/9823",  # Sus scrofa
}
"""The NCBI taxon of every species of a build."""
