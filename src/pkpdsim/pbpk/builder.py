"""The composition of a PBPK model from its definition, and its build.

The comp model has a top model and one submodel per organ:

- the top model holds the physiology (covariates, volumes, flows), the blood
  pools (arterial, venous, portal; plasma `C<pool>_<s>` and kinetic blood
  cells `C<pool>_bc_<s>`), the parameters of the substances (`Mr`, `fu`, `Kp`,
  `BP`, ...), the sinks (urine, feces, bile), the intestinal lumen and the
  routes;
- every organ is a `ModelDefinition` with a `Submodel` (`LI`, `KI`, ...), the
  standard organ or a variant; it reads the elements of the body through its
  ports, and its outflow concentrations replace the parameters
  `Cout_<o>_<s>` of the top model.

Every port is declared by the organ (its `Boundary`) and every replacement is
generated here from it; the flat model has the ids of the scheme
(`LI__Cpl_<s>` in the organs, the ids of the top model unprefixed).

Until the segmented intestine the lumen is a single compartment `Vlumen`
(the fluid of the small intestine) in the top model: the oral doses enter it,
a substance is absorbed from it with first order `ka_<s>` into the portal
vein (the gut wall is passed without distribution) and leaves it with the
small intestinal transit into the feces; the bile drains into it.
"""

from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from sbmlutils.comp import flatten_sbml
from sbmlutils.factory import (
    PORT_SUFFIX,
    AssignmentRule,
    Compartment,
    Model,
    ModelDefinition,
    Parameter,
    Port,
    Reaction,
    ReplacedBy,
    ReplacedElement,
    Sbase,
    Species,
    Submodel,
    UnitDefinition,
    create_model,
)
from sbmlutils.metadata import SBO
from sbmlutils.validation import validate_doc

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.data import cite, value
from pkpdsim.pbpk.definition import PBPKDefinition
from pkpdsim.pbpk.organ import Boundary, OrganComponent, OrganContext, StandardOrgan
from pkpdsim.pbpk.organs import Organ, Part, Pool
from pkpdsim.pbpk.parts import parts_for
from pkpdsim.pbpk.physiology import physiology_objects, reference_tissue_volumes
from pkpdsim.pbpk.placement import kp_id, value_parameter
from pkpdsim.pbpk.routes import route_objects
from pkpdsim.pbpk.substance import Kp, KpMethod, Substance
from pkpdsim.pbpk.units import MODEL_UNITS, U
from pkpdsim.pbpk.values import Value

__all__ = ["BuildResult", "build", "compose", "organ_component", "organ_context"]

_ONE = "1 dimensionless"
_SBO_PARAMETER = SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER
_VARIANTS: dict[str, OrganComponent] = {"standard": StandardOrgan()}
"""The organ variants by name; the registry of the variants extends it."""


@dataclass(frozen=True)
class BuildResult:
    """The files of a built model.

    Attributes:
        comp_path: the comp model, the organs as submodels.
        flat_path: the flattened model, which is simulated.
    """

    comp_path: Path
    flat_path: Path


def organ_component(definition: PBPKDefinition, organ: Organ) -> OrganComponent:
    """Get the component of an organ: its variant, the standard organ by default.

    Args:
        definition: the definition.
        organ: the organ.

    Returns:
        The component.

    Raises:
        ValueError: if the variant of the organ is a name of no variant.
    """
    variant = definition.variants.organs.get(organ, "standard")
    if not isinstance(variant, str):
        return variant
    try:
        return _VARIANTS[variant]
    except KeyError:
        raise ValueError(
            f"Unknown variant '{variant}' of organ '{organ.value}', the variants "
            f"are {sorted(_VARIANTS)}."
        ) from None


def organ_context(definition: PBPKDefinition, organ: Organ) -> OrganContext:
    """Get the context of an organ: every substance flows through every organ.

    Args:
        definition: the definition.
        organ: the organ.

    Returns:
        The substances with their parts in the organ, the chemistry and the
        inhibitions and inductions of the substances.
    """
    chemistry = tuple(definition.chemistry)
    substances = tuple(definition.substances)
    return OrganContext(
        organ=organ,
        substances=substances,
        parts={s.sid: parts_for(s, organ, chemistry) for s in substances},
        chemistry=chemistry,
        inhibitions=tuple((s.sid, i) for s in substances for i in s.inhibits),
        inductions=tuple((s.sid, i) for s in substances for i in s.induces),
    )


def compose(definition: PBPKDefinition) -> Model:
    """Compose the comp model of a definition.

    Args:
        definition: the definition.

    Returns:
        The comp model: the top model with one submodel per organ.

    Raises:
        ValueError: if the definition is invalid, an organ variant is unknown,
            an organ reads an element the body does not have, or a route is
            given twice.
        NotImplementedError: for a Kp method or a PD component, which later
            versions of the library build.
    """
    definition.validate()
    if definition.pd:
        raise NotImplementedError(
            f"{definition.sid}: PD components are not built yet by the library."
        )
    top = _Top(definition)
    for s in definition.substances:
        top.substance(s)
    for pool in Pool:
        top.pool(pool)
    top.sinks()
    top.routes()

    model_definitions: list[ModelDefinition] = []
    submodels: list[Submodel] = []
    replaced: list[ReplacedElement] = []
    for organ in Organ:
        md, boundary = organ_component(definition, organ).build(
            organ_context(definition, organ)
        )
        if any(md.sid == other.sid for other in model_definitions):
            raise ValueError(f"Two organs build the model definition '{md.sid}'.")
        _scope_port_metaids(md)
        model_definitions.append(md)
        submodel = ids.submodel(organ)
        if submodel in top.known:
            raise ValueError(
                f"The submodel id '{submodel}' of organ '{organ.value}' is an id "
                "of the top model."
            )
        submodels.append(Submodel(submodel, modelRef=md.sid, name=organ.name.lower()))
        replaced += top.connect(organ, submodel, boundary)

    return Model(
        sid=definition.sid,
        name=definition.name or definition.sid,
        notes=definition.notes,
        units=U,
        model_units=MODEL_UNITS,
        objects=top.objects,
        model_definitions=model_definitions,
        submodels=submodels,
        replaced_elements=replaced,
    )


def build(definition: PBPKDefinition, output_dir: Path) -> BuildResult:
    """Build a definition: the comp model and the flat model.

    Args:
        definition: the definition.
        output_dir: the folder of the files `<sid>.xml` and `<sid>_flat.xml`.

    Returns:
        The paths of the comp and the flat model.

    Raises:
        ValueError: if the definition cannot be composed, or a model has SBML
            errors.
    """
    model = compose(definition)
    comp_path = output_dir / f"{definition.sid}.xml"
    flat_path = output_dir / f"{definition.sid}_flat.xml"
    create_model(model, filepath=comp_path, validate=True, raise_on_error=True)
    doc = flatten_sbml(comp_path, flat_path)
    validation = validate_doc(doc, title=flat_path.name)
    if validation.error_count:
        raise ValueError(
            f"The flat model '{flat_path}' has {validation.error_count} SBML errors."
        )
    return BuildResult(comp_path=comp_path, flat_path=flat_path)


class _Top:
    """The objects of the top model, collected while the definition is composed."""

    def __init__(self, definition: PBPKDefinition) -> None:
        self.definition = definition
        self.objects: list[Sbase] = list(physiology_objects(definition.sex))
        self.known: set[str] = {
            str(o.sid) for o in self.objects if isinstance(o, Parameter)
        }

    def add(self, objects: Iterable[Sbase]) -> None:
        """Add objects, recording the ids of compartments, species and parameters."""
        for o in objects:
            if isinstance(o, Parameter | Compartment | Species):
                self.known.add(str(o.sid))
            self.objects.append(o)

    def _rule(
        self,
        pid: str,
        formula: str,
        unit: UnitDefinition,
        name: str,
        text: str,
    ) -> None:
        """A parameter with its assignment rule."""
        self.add(
            [
                Parameter(
                    pid,
                    None,
                    unit,
                    constant=False,
                    name=name,
                    sboTerm=_SBO_PARAMETER,
                    notes=f"{text}\n\n`{pid} = {formula}`",
                ),
                AssignmentRule(pid, formula, unit),
            ]
        )

    def _compartment(
        self,
        cid: str,
        formula: str,
        name: str,
        text: str,
        sbo: str = SBO.PHYSICAL_COMPARTMENT,
    ) -> None:
        """A compartment with its volume as assignment rule."""
        self.add(
            [
                Compartment(
                    cid,
                    None,
                    U.litre,
                    constant=False,
                    name=name,
                    sboTerm=sbo,
                    notes=f"{text}\n\n`{cid} = {formula}`",
                ),
                AssignmentRule(cid, formula, U.litre),
            ]
        )

    def _species(
        self, sid: str, compartment: str, name: str, amount: bool = False
    ) -> None:
        """A species of the top model, initially zero."""
        self.add(
            [
                Species(
                    sid,
                    compartment,
                    initialAmount=0.0 if amount else None,
                    initialConcentration=None if amount else 0.0,
                    substanceUnit=U.mmole,
                    hasOnlySubstanceUnits=amount,
                    name=name,
                    sboTerm=SBO.SIMPLE_CHEMICAL,
                )
            ]
        )

    def _reaction(self, sid: str, equation: str, rate: str, name: str) -> None:
        """A transport of the top model, the rate in mmole/min."""
        self.add(
            [
                Reaction(
                    sid=sid,
                    equation=equation,
                    formula=(rate, U.mmole_per_min),
                    name=name,
                    sboTerm=SBO.TRANSPORT_REACTION,
                )
            ]
        )

    # --- substances ------------------------------------------------------

    def substance(self, s: Substance) -> None:
        """The parameters of a substance: Mr, fu, the distribution and the blood.

        Raises:
            NotImplementedError: for a Kp method.
        """
        sid = s.sid
        self.add(
            [
                value_parameter(
                    ids.mr(sid),
                    s.mr,
                    U.g_per_mole,
                    f"molecular weight {sid} [g/mole]",
                    SBO.MOLECULAR_MASS,
                ),
                value_parameter(
                    ids.fu_ref(sid),
                    s.binding.fu if s.binding is not None else 1.0,
                    U.dimensionless,
                    f"reference fraction unbound {sid}",
                    _SBO_PARAMETER,
                    "Fraction unbound in plasma of the reference individual"
                    + (
                        f" (binding protein {s.binding.protein})."
                        if s.binding is not None
                        else "; the substance does not bind to plasma proteins."
                    ),
                ),
            ]
        )
        self._rule(
            ids.fu(sid),
            ids.fu_ref(sid),
            U.dimensionless,
            f"fraction unbound {sid}",
            "Fraction unbound in plasma of the individual.",
        )
        self._distribution(s)
        self._erythrocytes(s)

    def _distribution(self, s: Substance) -> None:
        """The tissue to plasma coefficients, ftissue and the reference tissue volume."""
        dist = s.distribution
        if dist is None:
            return
        sid = s.sid
        if isinstance(dist, KpMethod):
            raise NotImplementedError(
                f"{sid}: Kp methods are not built yet by the library, give a Kp."
            )
        assert isinstance(dist, Kp)
        objects = [
            value_parameter(
                ids.kp(sid),
                dist.value,
                U.dimensionless,
                f"tissue to plasma coefficient {sid}",
                _SBO_PARAMETER,
                "Whole-tissue (interstitial and intracellular) to plasma "
                "coefficient of the organs without an own value.",
            )
        ]
        for organ, q in dist.overrides.items():
            objects.append(
                value_parameter(
                    ids.kp_organ(organ, sid),
                    q,
                    U.dimensionless,
                    f"tissue to plasma coefficient {sid} {organ.name.lower()}",
                    _SBO_PARAMETER,
                    "Whole-tissue to plasma coefficient of the organ.",
                )
            )
        if dist.ftissue is not None:
            objects.append(
                value_parameter(
                    ids.ftissue(sid),
                    dist.ftissue,
                    U.l_per_min,
                    f"whole-body exchange {sid} [l/min]",
                    _SBO_PARAMETER,
                    "Exchange between plasma and tissue of the whole body of the "
                    "reference individual; an organ has the share of its tissue "
                    "volume. It multiplies `fu*(Cpl - Cti/Kp)`, so its fitted "
                    "values differ from those of the hand-written models.",
                )
            )
        shared = self._ftissue_organs(s)
        if shared:
            volumes = reference_tissue_volumes(self.definition.sex)
            total = sum(volumes[o] for o in shared)
            objects.append(
                Parameter(
                    ids.tissue_volume_total_ref(sid),
                    total,
                    U.litre,
                    constant=True,
                    name=f"reference tissue volume of the organs {sid} enters [l]",
                    sboTerm=_SBO_PARAMETER,
                    notes="Tissue volume `V<o>*(1 - f_vas_<o>)` of the "
                    f"{self.definition.sex.value} reference individual summed "
                    "over the organs whose exchange is a share of "
                    f"`{ids.ftissue(sid)}`: "
                    + ", ".join(o.name.lower() for o in shared)
                    + ". A constant, so that the exchange scales with the "
                    "tissue volume of the individual.",
                )
            )
        self.add(objects)

    def _ftissue_organs(self, s: Substance) -> list[Organ]:
        """The organs whose exchange of a substance is a share of `ftissue_<s>`.

        The merged and the full form of an organ the substance enters, unless
        the distribution gives the exchange of the organ (`ps`).
        """
        dist = s.distribution
        if dist is None:
            return []
        return [
            organ
            for organ in Organ
            if parts_for(s, organ, self.definition.chemistry).tissue
            in {"merged", "full"}
            and kp_id(s, organ) is not None
            and organ not in dist.ps
        ]

    def _erythrocytes(self, s: Substance) -> None:
        """The blood to plasma ratio and the parameters of the blood cells."""
        sid, ery = s.sid, s.erythrocytes
        bp, name = ids.bp(sid), f"blood to plasma ratio {sid}"
        if ery is None:
            self._rule(
                bp,
                f"{_ONE} - {ids.HCT}",
                U.dimensionless,
                name,
                "Blood to plasma ratio of a substance which does not enter the "
                "erythrocytes; an observer value, the flows carry plasma only.",
            )
            return
        if ery.bp is not None:
            self.add(
                [
                    value_parameter(
                        bp,
                        ery.bp,
                        U.dimensionless,
                        name,
                        _SBO_PARAMETER,
                        "Blood to plasma ratio, erythrocytes in equilibrium.",
                    )
                ]
            )
            return
        assert ery.kbc is not None
        self.add(
            [
                value_parameter(
                    ids.kbc(sid),
                    ery.kbc,
                    U.dimensionless,
                    f"blood cell to plasma coefficient {sid}",
                    _SBO_PARAMETER,
                )
            ]
        )
        if ery.k_exchange is not None:
            self.add(
                [
                    value_parameter(
                        ids.kbc_exchange(sid),
                        ery.k_exchange,
                        U.per_min,
                        f"exchange plasma and blood cells {sid} [1/min]",
                        _SBO_PARAMETER,
                        "Rate of the exchange between plasma and blood cells per "
                        "blood cell volume.",
                    )
                ]
            )
        self._rule(
            bp,
            f"{_ONE} - {ids.HCT} + {ids.HCT} * {ids.kbc(sid)}",
            U.dimensionless,
            name,
            "Blood to plasma ratio of the blood cell to plasma coefficient"
            + (
                " (at equilibrium of the kinetic blood cells)."
                if ery.mode == "kinetic"
                else "."
            ),
        )

    # --- pools -----------------------------------------------------------

    def pool(self, pool: Pool) -> None:
        """The compartments and species of a blood pool."""
        label = pool.name.lower()
        volume = ids.pool_volume(pool)
        plasma = ids.pool_compartment(pool, Part.PLASMA)
        cells = ids.pool_compartment(pool, Part.BLOOD_CELLS)
        self._compartment(
            plasma,
            f"{volume} * ({_ONE} - {ids.HCT})",
            f"{label} plasma",
            "Plasma of the blood pool.",
        )
        self._compartment(
            cells,
            f"{volume} * {ids.HCT}",
            f"{label} blood cells",
            "Blood cells of the blood pool.",
        )
        for s in self.definition.substances:
            ery = s.erythrocytes
            compartment = plasma
            if ery is not None and ery.mode == "equilibrium":
                compartment = ids.pool_capacity(pool, s.sid)
                self._compartment(
                    compartment,
                    f"{volume} * {ids.bp(s.sid)}",
                    f"plasma-equivalent capacity V*BP of {s.sid}, {label} pool",
                    "Not a physical volume: the plasma-equivalent capacity of the "
                    "blood of a substance with erythrocytes in equilibrium, the "
                    "pool volume times the blood to plasma ratio, so that the "
                    "plasma concentration in it gives the amount in the blood.",
                    SBO.IMPLICIT_COMPARTMENT,
                )
            cpl = ids.pool_species(pool, s.sid)
            self._species(cpl, compartment, f"{s.name} ({label} plasma)")
            if ery is None or ery.mode != "kinetic":
                continue
            cbc = ids.pool_species(pool, s.sid, Part.BLOOD_CELLS)
            self._species(cbc, cells, f"{s.name} ({label} blood cells)")
            self._reaction(
                ids.pool_exchange_reaction(pool, s.sid),
                f"{cpl} <-> {cbc}",
                f"{ids.kbc_exchange(s.sid)} * {cells} * "
                f"({cpl} - {cbc} / {ids.kbc(s.sid)})",
                f"exchange {s.sid} plasma and blood cells, {label} pool",
            )

    # --- sinks and lumen -------------------------------------------------

    def sinks(self) -> None:
        """The urine, the feces, the bile and the single intestinal lumen."""
        nominal = "Nominal volume of a compartment of amounts, which enters no rate."
        transit = value("intestine_human", "lumen", "transit_time")
        self.add(
            [
                Compartment(
                    ids.VURINE,
                    1.0,
                    U.litre,
                    constant=True,
                    name="urine",
                    sboTerm=SBO.PHYSICAL_COMPARTMENT,
                    notes=nominal,
                ),
                Compartment(
                    ids.VFECES,
                    1.0,
                    U.litre,
                    constant=True,
                    name="feces",
                    sboTerm=SBO.PHYSICAL_COMPARTMENT,
                    notes=nominal,
                ),
                _value_compartment(
                    ids.VBILE,
                    value("intestine_human", "bile", "volume"),
                    "bile",
                    "Bile, which drains into the intestinal lumen.",
                ),
                _value_compartment(
                    ids.VLUMEN,
                    value("intestine_human", "lumen", "volume"),
                    "intestinal lumen",
                    "The single intestinal lumen: the fluid of the small "
                    "intestine, into which the oral doses and the bile enter.",
                ),
                value_parameter(
                    ids.K_BILE,
                    value("intestine_human", "bile", "k_bile"),
                    U.per_min,
                    "transit bile to intestinal lumen [1/min]",
                    _SBO_PARAMETER,
                ),
                value_parameter(
                    ids.K_FECES,
                    1.0 / transit.value,
                    U.per_min,
                    "transit intestinal lumen to feces [1/min]",
                    _SBO_PARAMETER,
                    f"The reciprocal of the small intestinal transit time "
                    f"{transit.value} min ({transit.notes}).\n\nSource:\n\n- "
                    f"{cite(str(transit.source))}",
                ),
            ]
        )
        ka = value("intestine_human", "lumen", "ka")
        for s in self.definition.substances:
            sid = s.sid
            self._species(ids.urine(sid), ids.VURINE, f"{s.name} (urine)", amount=True)
            self._species(ids.feces(sid), ids.VFECES, f"{s.name} (feces)", amount=True)
            self._species(ids.bile(sid), ids.VBILE, f"{s.name} (bile)")
            self._species(ids.lumen(sid), ids.VLUMEN, f"{s.name} (intestinal lumen)")
            self.add(
                [
                    value_parameter(
                        ids.absorption(sid),
                        ka,
                        U.per_hr,
                        f"absorption {sid} [1/hr]",
                        _SBO_PARAMETER,
                        "First order absorption from the intestinal lumen into "
                        "the portal vein; a model fits it.",
                    )
                ]
            )
            self._reaction(
                ids.bile_reaction(sid),
                f"{ids.bile(sid)} -> {ids.lumen(sid)}",
                f"{ids.K_BILE} * {ids.VBILE} * {ids.bile(sid)}",
                f"transit {sid} bile to intestinal lumen",
            )
            self._reaction(
                ids.absorption_reaction(sid),
                f"{ids.lumen(sid)} -> {ids.pool_species(Pool.PORTAL, sid)}",
                f"{ids.absorption(sid)} / 60 min_per_hr * {ids.VLUMEN} * "
                f"{ids.lumen(sid)}",
                f"absorption {sid} intestinal lumen to portal vein",
            )
            self._reaction(
                ids.feces_reaction(sid),
                f"{ids.lumen(sid)} -> {ids.feces(sid)}",
                f"{ids.K_FECES} * {ids.VLUMEN} * {ids.lumen(sid)}",
                f"transit {sid} intestinal lumen to feces",
            )

    # --- routes ----------------------------------------------------------

    def routes(self) -> None:
        """The dose states of the routes.

        Raises:
            ValueError: if a route of a substance is given twice.
        """
        by_id = {s.sid: s for s in self.definition.substances}
        seen: set[tuple[str, str]] = set()
        for route in self.definition.routes:
            key = (route.kind, route.substance)
            if key in seen:
                raise ValueError(
                    f"Route '{route.kind}' of substance '{route.substance}' is "
                    "given twice."
                )
            seen.add(key)
            self.add(route_objects(route, by_id[route.substance]))

    # --- organs ----------------------------------------------------------

    def connect(
        self, organ: Organ, submodel: str, boundary: Boundary
    ) -> list[ReplacedElement]:
        """The replacements of the boundary of an organ.

        Every element of the body which the organ reads replaces the port of
        the organ; every output of the organ replaces a parameter
        `Cout_<o>_<s>` of the top model, which is created here.

        Raises:
            ValueError: if the organ reads an element the body does not have.
        """
        missing = sorted(
            t for t in boundary.replacements.values() if t not in self.known
        )
        if missing:
            raise ValueError(
                f"Organ '{organ.value}' reads elements which the body does not "
                f"have: {missing}."
            )
        for local, top in boundary.outputs.items():
            self.add(
                [
                    Parameter(
                        top,
                        0.0,
                        U.mM,
                        constant=False,
                        name=f"outflow plasma concentration {top} [mM]",
                        sboTerm=_SBO_PARAMETER,
                        notes=f"Plasma concentration leaving the organ "
                        f"{organ.name.lower()}, replaced by the output of the "
                        "organ; the value is a placeholder.",
                        replacedBy=ReplacedBy(
                            elementRef=top,
                            submodelRef=submodel,
                            portRef=f"{local}{PORT_SUFFIX}",
                        ),
                    )
                ]
            )
        return [
            ReplacedElement(
                elementRef=top,
                submodelRef=submodel,
                portRef=f"{local}{PORT_SUFFIX}",
            )
            for local, top in boundary.replacements.items()
        ]


def _scope_port_metaids(md: ModelDefinition) -> None:
    """Give the ports of a model definition metaids unique in the document.

    The port of an element is `<id>_port`, and sbmlutils gives a port the
    metaid of its id; a metaid is unique in the whole document, so the ports
    of two organs with the same local ids (`Vpl_port`) would collide. The
    metaid is prefixed with the id of the model definition, the port ids, to
    which the replacements refer, stay.
    """
    for name in ("compartments", "species", "parameters", "reactions"):
        for element in getattr(md, name) or []:
            port = getattr(element, "port", None)
            if port is None or port is False:
                continue
            if port is True:
                port = Port(
                    f"{element.sid}{PORT_SUFFIX}",
                    name=f"Port of {element.sid}",
                    sboTerm=SBO.PORT,
                )
                element.port = port
            port.metaId = f"{md.sid}__{port.sid}"


def _value_compartment(cid: str, v: Value, name: str, text: str) -> Compartment:
    """A constant compartment of a value of the data, its source in the notes."""
    return Compartment(
        cid,
        v.value,
        U.litre,
        constant=True,
        name=name,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        notes=f"{text} {v.notes}\n\nSource:\n\n- {cite(str(v.source))}",
    )
