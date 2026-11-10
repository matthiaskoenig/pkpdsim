"""The parameters and the observers every library of pkpdsim writes.

A constant of a model is the parameter of a number with its provenance
(`value_parameter`, whose notes carry the citation of the source and whose
`distrib` uncertainty carries the uncertainty of a fitted value,
`value_notes`); a constant which is tied to other parameters of the model is
a `FormulaParameter`, whose value is the initial assignment of its formula
(`formula_assignment`). An observer is a parameter which an assignment rule
sets (`observer`, `RULE`): read-only, so a model with observers has the same
states as a model without them.
"""

from sbmlutils.factory import (
    AssignmentRule,
    InitialAssignment,
    Parameter,
    UnitDefinition,
)
from sbmlutils.metadata import SBO

from pkpdsim.database.tables import cite
from pkpdsim.sbml.annotation import initial_assignment
from pkpdsim.sbml.values import Formula, Q, Value, as_value

__all__ = [
    "RULE",
    "FormulaParameter",
    "formula_assignment",
    "observer",
    "value_notes",
    "value_parameter",
]

RULE = SBO.MATHEMATICAL_EXPRESSION


class FormulaParameter(Parameter):
    """A constant parameter whose value is a `Formula` over ids of the top model.

    The parameter has no value; the build gives it the initial assignment of
    the formula (`formula_assignment`): the top model in its objects, an organ
    in its assignments, with an input port for every id of the top model the
    formula reads (`OrganBuilder.model_definition`).

    Attributes:
        formula: the formula.
    """

    def __init__(
        self, pid: str, formula: Formula, name: str, sbo: str, notes: str
    ) -> None:
        """Create the parameter of a formula."""
        super().__init__(
            pid,
            None,
            formula.unit,
            constant=True,
            name=name,
            sboTerm=sbo,
            notes=notes,
        )
        self.formula = formula


def formula_assignment(p: FormulaParameter) -> InitialAssignment:
    """Get the initial assignment of the parameter of a formula."""
    return initial_assignment(
        str(p.sid), p.formula.expression, p.formula.unit, str(p.name)
    )


def value_parameter(
    pid: str,
    q: Q,
    unit: UnitDefinition,
    name: str,
    sbo: str,
    text: str = "",
) -> Parameter:
    """Get a constant parameter of a number with its provenance.

    Args:
        pid: the id of the parameter.
        q: the number, plain (in `unit`) or a `Value` with unit and source,
            or a `Formula` (a `FormulaParameter`, the notes name the formula).
        unit: the unit of a plain number.
        name: the name of the parameter.
        sbo: the SBO term.
        text: the description of the parameter.

    Returns:
        The parameter, the citation of the source in the notes and the
        uncertainty of a fitted value.
    """
    v = as_value(q, unit)
    lines = [text] if text else []
    if isinstance(v, Formula):
        if v.notes:
            lines.append(v.notes)
        lines.append(
            f"Tied to other parameters of the model: `{pid} = {v.expression}`."
        )
        return FormulaParameter(pid, v, name, sbo, "\n\n".join(lines))
    return Parameter(
        pid,
        v.value,
        v.unit,
        constant=True,
        name=name,
        sboTerm=sbo,
        notes=value_notes(v, text),
        uncertainties=v.uncertainties,
    )


def value_notes(v: Value, text: str = "") -> str | None:
    """Get the notes of a number with its provenance.

    Args:
        v: the value.
        text: the description of the element which carries the value.

    Returns:
        The description, the notes of the value and the citation of its
        source (`Fitted value.` for a fit), `None` without any.
    """
    lines = [text] if text else []
    if v.notes:
        lines.append(v.notes)
    if v.source == "fit":
        lines.append("Fitted value.")
    elif v.source is not None:
        lines.append(f"Source:\n\n- {cite(v.source)}")
    return "\n\n".join(lines) if lines else None


def observer(
    sid: str, formula: str, unit: UnitDefinition, name: str, sbo: str
) -> tuple[Parameter, AssignmentRule]:
    """Get the parameter of an observer and its assignment rule."""
    return (
        Parameter(sid, None, unit, constant=False, name=name, sboTerm=sbo),
        AssignmentRule(sid, formula, unit, name=name, sboTerm=RULE),
    )
