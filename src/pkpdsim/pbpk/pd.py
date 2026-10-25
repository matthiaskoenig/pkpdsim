"""Pharmacodynamic and physiology components, coupled to the body by their ports.

A component is an sbmlutils model of the package (a PD model, a model of
body weight or of the electrolytes) with declared ports. It reads quantities
of the body through its input ports and changes quantities of the body
through its output ports; the builder adds it as a comp `ModelDefinition`
with a `Submodel` of the id of the model and generates the replacements, so
the package writes no `ReplacedElement`:

- an input: the element of the body (`Carm_free_api`, a species, an
  observer) replaces the element of the input port of the component, which
  reads it under its own id;
- an output: the element of the output port of the component replaces the
  parameter of the body (`f_cardiac_function`), which every quantity of the
  body reads; the parameter is dynamic (`PBPKDefinition.dynamic`), so the
  quantities which depend on it are assignment rules and follow it during
  the simulation.

The flat model has the elements of the component prefixed by its id
(`<sid>__E`), the replaced parameter of the body keeps its id.
"""

from collections.abc import Mapping
from dataclasses import dataclass

from sbmlutils.factory import PORT_SUFFIX, Model, Port, PortType

__all__ = ["PDComponent", "ports"]


@dataclass(frozen=True, eq=False)
class PDComponent:
    """A model with declared ports, connected to the body by name.

    Two components are equal when their models have the same id and they
    have the same inputs and outputs (the sbmlutils `Model` compares by
    identity), so the component which two packages construct each for
    itself is one component of a coupled model (`PBPKDefinition.combine`).

    Attributes:
        model: the model of the component; its id is the id of its submodel
            and the prefix of its elements in the flat model.
        inputs: the elements of the body which the component reads, by the
            id of the element of its input port: `{"C_drug": "Carm_free_api"}`
            reads the unbound arm vein concentration as `C_drug`.
        outputs: the parameters of the body which the component changes, by
            the id of the element of its output port: `{"f_cardiac":
            "f_cardiac_function"}` replaces the cardiac function of the body
            by `f_cardiac`.
    """

    model: Model
    inputs: Mapping[str, str]
    outputs: Mapping[str, str]

    def _key(
        self,
    ) -> tuple[str, frozenset[tuple[str, str]], frozenset[tuple[str, str]]]:
        """The model id, the inputs and the outputs, by which components compare."""
        return (
            str(self.model.sid),
            frozenset(self.inputs.items()),
            frozenset(self.outputs.items()),
        )

    def __eq__(self, other: object) -> bool:
        """Compare by the model id, the inputs and the outputs."""
        if not isinstance(other, PDComponent):
            return NotImplemented
        return self._key() == other._key()

    def __hash__(self) -> int:
        """Hash the model id, the inputs and the outputs."""
        return hash(self._key())


def ports(model: Model) -> dict[str, tuple[str, PortType | None]]:
    """Get the ports of a model by the id of their element.

    Args:
        model: the model.

    Returns:
        The id and the type of the port of every element with a port
        (`port=True` gives `<id>_port`, a `Port` its id and type) and of
        every port of the model which names its element (`idRef`).
    """
    result: dict[str, tuple[str, PortType | None]] = {}
    for name in ("compartments", "species", "parameters", "reactions"):
        for element in getattr(model, name) or []:
            port = getattr(element, "port", None)
            if port is True:
                result[str(element.sid)] = (f"{element.sid}{PORT_SUFFIX}", None)
            elif isinstance(port, Port):
                result[str(element.sid)] = (str(port.sid), port.portType)
    for port in model.ports or []:
        if port.idRef is not None:
            result[str(port.idRef)] = (str(port.sid), port.portType)
    return result
