from __future__ import annotations
from collections import Counter

import importlib
import importlib.util
import os
import sys
from typing import Any, Callable
from enum import StrEnum
import yaml
import cmath
from symbolica import Expression as SBE, AtomType  # type: ignore
from gammaloop.misc.common import pjoin, DATA_PATH, GammaLoopError, logger, GL_WARNINGS_ISSUED, GammaLoopWarning
import gammaloop.misc.utils as utils
from gammaloop.base_objects.param_card import ParamCard
pjoin = os.path.join


class ParameterType(StrEnum):
    COMPLEX = 'complex'
    REAL = 'real'


class ParameterNature(StrEnum):
    INTERNAL = 'internal'
    EXTERNAL = 'external'


class SerializableVertexRule(object):
    def __init__(self, name: str, particles: list[str], color_structures: list[str], lorentz_structures: list[str], couplings: list[list[str | None]]):
        self.name: str = name
        self.particles: list[str] = particles
        self.color_structures: list[str] = color_structures
        self.lorentz_structures: list[str] = lorentz_structures
        self.couplings: list[list[str | None]] = couplings

    @staticmethod
    def from_vertex_rule(vertex_rule: VertexRule) -> SerializableVertexRule:
        return SerializableVertexRule(
            vertex_rule.name,
            [particle.name for particle in vertex_rule.particles],
            [utils.expression_to_string_safe(c)
             for c in vertex_rule.color_structures],
            [lorentz.name for lorentz in vertex_rule.lorentz_structures],
            [
                [
                    None if (i_col, j_lor) not in vertex_rule.couplings else vertex_rule.couplings[(
                        i_col, j_lor)].name
                    for j_lor in range(len(vertex_rule.lorentz_structures))
                ] for i_col in range(len(vertex_rule.color_structures))
            ]
        )

    @staticmethod
    def from_dict(dict_repr: dict[str, Any]) -> SerializableVertexRule:
        return SerializableVertexRule(
            dict_repr['name'],
            dict_repr['particles'],
            dict_repr['color_structures'],
            dict_repr['lorentz_structures'],
            dict_repr['couplings']
        )


class VertexRule(object):
    def __init__(self, name: str, particles: list[Particle], color_structures: list[str], lorentz_structures: list[LorentzStructure], couplings: dict[tuple[int, int], Coupling]):
        self.name: str = name
        self.particles: list[Particle] = particles
        self.color_structures: list[SBE] = [utils.parse_python_expression_safe(
            color_structure) for color_structure in color_structures]
        self.lorentz_structures: list[LorentzStructure] = lorentz_structures
        self.couplings: dict[tuple[int, int], Coupling] = couplings

    @staticmethod
    def from_ufo_object(model: Model, ufo_object: Any) -> VertexRule:
        return VertexRule(
            ufo_object.name,
            [model.get_particle(particle.name)
             for particle in ufo_object.particles],
            ufo_object.color,
            [model.get_lorentz_structure(lorentz.name)
             for lorentz in ufo_object.lorentz],
            {(i, j): model.get_coupling(coupling.name)
             for (i, j), coupling in ufo_object.couplings.items()}
        )

    @staticmethod
    def from_serializable_vertex_rule(model: Model, serializable_vertex_rule: SerializableVertexRule) -> VertexRule:
        return VertexRule(
            serializable_vertex_rule.name,
            [model.get_particle(particle_name)
             for particle_name in serializable_vertex_rule.particles],
            serializable_vertex_rule.color_structures,
            [model.get_lorentz_structure(
                lorentz_name) for lorentz_name in serializable_vertex_rule.lorentz_structures],
            dict(((i, j), model.get_coupling(serializable_vertex_rule.couplings[i][j]))  # type: ignore
                 for i in range(len(serializable_vertex_rule.color_structures))
                 for j in range(len(serializable_vertex_rule.lorentz_structures)) if serializable_vertex_rule.couplings[i][j] is not None)
        )

    def to_serializable_vertex_rule(self) -> SerializableVertexRule:
        return SerializableVertexRule.from_vertex_rule(self)


class Propagator(object):
    def __init__(self, name: str, particle: Particle, numerator: str, denominator: str):
        self.name: str = name
        self.particle: Particle = particle
        self.numerator: SBE = utils.parse_python_expression_safe(numerator)
        self.denominator: SBE = utils.parse_python_expression_safe(denominator)

    @staticmethod
    def from_particle(particle: Particle, gauge: str) -> Propagator:

        name = f'{particle.name}_prop{gauge}'
        if particle.spin == -1:  # ghost
            numerator = '1'
            if particle.is_massive():
                denominator = f'P(1)**2-{particle.mass.name}**2'
            else:
                denominator = 'P(1)**2'
        elif particle.spin == 1:  # scalar 2s+1=1
            if particle.is_massive():
                numerator = '1'
                denominator = f'P(1)**2-{particle.mass.name}**2'
            else:
                numerator = '1'
                denominator = 'P(1)**2'
        elif particle.spin == 2:  # spinor 2s+1=2
            if particle.is_massive():
                numerator = f'PSlash(1,2)+{particle.mass.name}*Identity(1,2)'
                denominator = f'P(1)**2-{particle.mass.name}**2'
            else:
                numerator = 'PSlash(1,2)'
                denominator = 'P(1)**2'
        elif particle.spin == 3:  # vector 2s+1=3
            if particle.is_massive():
                numerator = f'-Metric(1,2)+P(1)*P(2)/{particle.mass.name}**2'
                denominator = f'P(1)**2-{particle.mass.name}**2'
            else:
                denominator = 'P(1)**2'
                if gauge == 'Feynman':
                    numerator = '-Metric(1,2)'
                else:
                    raise GammaLoopError(
                        f'Gauge {gauge} not implemented for vector particles')
        else:
            raise GammaLoopError(
                f'Particle spin {particle.spin} not implemented')

        return Propagator(
            name,
            particle,
            numerator,
            denominator
        )

    @staticmethod
    def from_serializable_propagator(model: Model, serializable_propagator: SerializablePropagator, ) -> Propagator:
        return Propagator(
            serializable_propagator.name,
            model.get_particle(serializable_propagator.particle),
            serializable_propagator.numerator,
            serializable_propagator.denominator
        )

    def to_serializable_propagator(self) -> SerializablePropagator:
        return SerializablePropagator.from_propagator(self)


class SerializablePropagator(object):
    def __init__(self, name: str, particle: str, numerator: str, denominator: str):
        self.name: str = name
        self.particle: str = particle
        self.numerator: str = numerator
        self.denominator: str = denominator

    @staticmethod
    def from_propagator(propagator: Propagator) -> SerializablePropagator:
        return SerializablePropagator(propagator.name, propagator.particle.name, utils.expression_to_string_safe(propagator.numerator), utils.expression_to_string_safe(propagator.denominator))

    @staticmethod
    def from_dict(dict_repr: dict[str, Any]) -> SerializablePropagator:
        return SerializablePropagator(
            dict_repr['name'],
            dict_repr['particle'],
            dict_repr['numerator'],
            dict_repr['denominator']
        )


class SerializableCoupling(object):
    def __init__(self, name: str, expression: str, orders: list[tuple[str, int]], value: tuple[float, float] | None):
        self.name: str = name
        self.expression: str = expression
        self.value: tuple[float, float] | None = value
        self.orders: list[tuple[str, int]] = orders

    @staticmethod
    def from_coupling(coupling: Coupling) -> SerializableCoupling:
        return SerializableCoupling(
            coupling.name,
            utils.expression_to_string_safe(coupling.expression),
            list(coupling.orders.items()),
            None if coupling.value is None else (
                coupling.value.real, coupling.value.imag)
        )

    @staticmethod
    def from_dict(dict_repr: dict[str, Any]) -> SerializableCoupling:
        return SerializableCoupling(
            dict_repr['name'],
            dict_repr['expression'],
            dict_repr['orders'],
            dict_repr['value']
        )


class Coupling(object):
    def __init__(self, name: str, expression: str, orders: dict[str, int], value: complex | None = None):
        self.name: str = name
        self.expression: SBE = utils.parse_python_expression_safe(expression)
        self.value: complex | None = value
        self.orders: dict[str, int] = orders

    @staticmethod
    def from_ufo_object(ufo_object: Any) -> Coupling:
        expression = ufo_object.value
        if isinstance(ufo_object.value, dict):
            if list(ufo_object.value.keys()) == [0,]:
                expression = ufo_object.value[0]
            else:
                if GammaLoopWarning.DroppingEpsilonTerms not in GL_WARNINGS_ISSUED:
                    GL_WARNINGS_ISSUED.add(
                        GammaLoopWarning.DroppingEpsilonTerms)
                    logger.warning(
                        "Only finite terms of order ε^0 are retained in model expressions")
                expression = ufo_object.value.get(0, '0')
        return Coupling(ufo_object.name, expression, ufo_object.order)

    @staticmethod
    def from_serializable_coupling(serializable_coupling: SerializableCoupling) -> Coupling:
        return Coupling(
            serializable_coupling.name,
            serializable_coupling.expression,
            dict(serializable_coupling.orders),
            None if serializable_coupling.value is None else complex(
                *serializable_coupling.value)
        )

    def to_serializable_coupling(self) -> SerializableCoupling:
        return SerializableCoupling.from_coupling(self)


class SerializableLorentzStructure(object):
    def __init__(self, name: str, spins: list[int], structure: str):
        self.name: str = name
        self.spins: list[int] = spins
        self.structure: str = structure

    @staticmethod
    def from_lorentz_structure(lorentz_structure: LorentzStructure) -> SerializableLorentzStructure:
        return SerializableLorentzStructure(
            lorentz_structure.name,
            lorentz_structure.spins,
            utils.expression_to_string_safe(lorentz_structure.structure)
        )

    @staticmethod
    def from_dict(dict_repr: dict[str, Any]) -> SerializableLorentzStructure:
        return SerializableLorentzStructure(
            dict_repr['name'],
            dict_repr['spins'],
            dict_repr['structure']
        )


class LorentzStructure(object):
    def __init__(self, name: str, spins: list[int], structure: str):
        self.name: str = name
        self.spins: list[int] = spins
        self.structure: SBE = utils.parse_python_expression_safe(structure)

    @staticmethod
    def from_ufo_object(ufo_object: Any) -> LorentzStructure:
        return LorentzStructure(ufo_object.name, ufo_object.spins, ufo_object.structure)

    @staticmethod
    def from_serializable_lorentz_structure(serializable_lorentz_structure: SerializableLorentzStructure) -> LorentzStructure:
        return LorentzStructure(
            serializable_lorentz_structure.name,
            serializable_lorentz_structure.spins,
            serializable_lorentz_structure.structure
        )

    def to_serializable_lorentz_structure(self) -> SerializableLorentzStructure:
        return SerializableLorentzStructure.from_lorentz_structure(self)


class SerializableParameter(object):
    def __init__(self, lhablock: str | None, lhacode: tuple[int, ...] | None, name: str, nature: str, parameter_type: str, value: tuple[float, float] | None, expression: str | None):
        self.lhablock: str | None = lhablock
        self.lhacode: tuple[int, ...] | None = lhacode
        self.name: str = name
        self.nature: str = nature
        self.parameter_type: str = parameter_type
        self.value: tuple[float, float] | None = value
        self.expression: str | None = expression

    @staticmethod
    def from_parameter(parameter: Parameter) -> SerializableParameter:
        return SerializableParameter(
            parameter.lhablock, parameter.lhacode, parameter.name,
            str(parameter.nature), str(
                parameter.parameter_type),
            None if parameter.value is None else (
                parameter.value.real, parameter.value.imag),
            utils.expression_to_string(parameter.expression)
        )

    @staticmethod
    def from_dict(dict_repr: dict[str, Any]) -> SerializableParameter:
        return SerializableParameter(
            dict_repr['lhablock'], dict_repr['lhacode'], dict_repr['name'],
            dict_repr['nature'], dict_repr['parameter_type'], dict_repr['value'],
            dict_repr['expression']
        )


class Parameter(object):

    # Useful for assigning typed defaults
    @staticmethod
    def default() -> Parameter:
        return Parameter(lhablock='DUMMY', lhacode=(0,), name='DUMMY', nature=ParameterNature.INTERNAL, parameter_type=ParameterType.REAL, value=0., expression='0')

    def __init__(self, lhablock: str | None, lhacode: tuple[int, ...] | None, name: str, nature: ParameterNature, parameter_type: ParameterType, value: complex | None, expression: str | None):
        self.lhablock: str | None = lhablock
        self.lhacode: tuple[int, ...] | None = lhacode
        self.name: str = name
        self.nature: ParameterNature = nature
        self.parameter_type: ParameterType = parameter_type
        self.value: complex | None = value
        self.expression: SBE | None = utils.parse_python_expression(expression)

    @staticmethod
    def from_ufo_object(ufo_object: Any) -> Parameter:
        # This ugly hack is to fix an annoying persistent typo in many NLO models
        if ufo_object.nature == 'interal':
            ufo_object.nature = 'internal'

        if ufo_object.nature == 'external':
            param_value = complex(ufo_object.value)
            param_expression: str | None = None
        elif ufo_object.nature == 'internal':
            expression: Any = ufo_object.value
            if isinstance(ufo_object.value, dict):
                if list(ufo_object.value.keys()) == [0,]:
                    expression = ufo_object.value[0]
                else:
                    if GammaLoopWarning.DroppingEpsilonTerms not in GL_WARNINGS_ISSUED:
                        GL_WARNINGS_ISSUED.add(
                            GammaLoopWarning.DroppingEpsilonTerms)
                        logger.warning(
                            "Only finite terms of order ε^0 are retained in model expressions")
                    expression = ufo_object.value.get(0, '0')
            param_expression: str | None = expression
            try:
                param_value = complex(expression)
                param_expression = None
            except ValueError:
                param_value = None
        else:
            raise GammaLoopError(
                f"Invalid parameter nature '{ufo_object.nature}'")

        return Parameter(
            ufo_object.lhablock if hasattr(ufo_object, 'lhablock') else None,
            (None if ufo_object.lhacode is None else tuple(ufo_object.lhacode)
             ) if hasattr(ufo_object, 'lhacode') else None,
            ufo_object.name,
            ParameterNature(ufo_object.nature), ParameterType(ufo_object.type),
            param_value, param_expression
        )

    @staticmethod
    def from_serializable_parameter(serializable_parameter: SerializableParameter) -> Parameter:
        return Parameter(
            serializable_parameter.lhablock, serializable_parameter.lhacode, serializable_parameter.name,
            ParameterNature(serializable_parameter.nature), ParameterType(
                serializable_parameter.parameter_type),
            None if serializable_parameter.value is None else complex(
                *serializable_parameter.value),
            serializable_parameter.expression
        )

    def to_serializable_parameter(self) -> SerializableParameter:
        return SerializableParameter.from_parameter(self)

    def __str__(self) -> str:
        if self.nature == ParameterNature.INTERNAL:
            return f"Internal parameter({self.name}, {self.expression})"
        else:
            return f"External parameter({self.name}, {self.value})"


class SerializableParticle(object):
    def __init__(self, pdg_code: int, name: str, antiname: str, spin: int, color: int, mass: str, width: str, texname: str, antitexname: str, charge: float, ghost_number: int, lepton_number: int, y_charge: int):
        self.pdg_code: int = pdg_code
        self.name: str = name
        self.antiname: str = antiname
        self.spin: int = spin
        self.color: int = color
        self.mass: str = mass
        self.width: str = width
        self.texname: str = texname
        self.antitexname: str = antitexname
        self.charge: float = charge
        self.ghost_number: int = ghost_number
        self.lepton_number: int = lepton_number
        self.y_charge: int = y_charge

    @classmethod
    def from_particle(cls, particle: Particle) -> SerializableParticle:
        return SerializableParticle(
            particle.pdg_code, particle.name, particle.antiname, particle.spin, particle.color,
            particle.mass.name,
            particle.width.name,
            particle.texname, particle.antitexname,
            particle.charge,
            particle.ghost_number,
            particle.lepton_number,
            particle.y_charge
        )

    @classmethod
    def from_dict(cls, dict_repr: dict[str, Any]) -> SerializableParticle:
        return SerializableParticle(
            dict_repr['pdg_code'], dict_repr['name'], dict_repr['antiname'], dict_repr['spin'], dict_repr['color'],
            dict_repr['mass'],
            dict_repr['width'],
            dict_repr['texname'], dict_repr['antitexname'],
            dict_repr['charge'],
            dict_repr['ghost_number'],
            dict_repr['lepton_number'],
            dict_repr['y_charge']
        )


class Particle(object):
    def __init__(self, pdg_code: int, name: str, antiname: str, spin: int, color: int, mass: Parameter, width: Parameter, texname: str, antitexname: str, charge: float, ghost_number: int, lepton_number: int, y_charge: int):
        self.pdg_code: int = pdg_code
        self.name: str = name
        self.antiname: str = antiname
        self.spin: int = spin
        self.color: int = color
        self.mass: Parameter = mass
        self.width: Parameter = width
        self.texname: str = texname
        self.antitexname: str = antitexname
        self.charge: float = charge
        self.ghost_number: int = ghost_number
        self.lepton_number: int = lepton_number
        self.y_charge: int = y_charge

    @staticmethod
    def default() -> Particle:
        return Particle(0, '', '', 1, 1, Parameter.default(), Parameter.default(), '', '', 0, 0, 0, 0)

    def is_ghost(self) -> bool:
        return self.ghost_number != 0

    def is_massive(self) -> bool:
        return self.mass.value is None or abs(self.mass.value) > 0.

    @staticmethod
    def sanitize_texname(texname: str) -> str:
        # Some UFO model unfortunately include ~ and _ in texnames, which are not allowed in LaTeX
        return texname.replace('~', 'x').replace('_', '')

    @staticmethod
    def from_ufo_object(model: Model, ufo_object: Any) -> Particle:

        return Particle(
            ufo_object.pdg_code, ufo_object.name, ufo_object.antiname, ufo_object.spin, ufo_object.color,
            model.get_parameter(ufo_object.mass.name),
            model.get_parameter(ufo_object.width.name),
            Particle.sanitize_texname(ufo_object.texname),
            Particle.sanitize_texname(ufo_object.antitexname),
            ufo_object.charge,
            ufo_object.GhostNumber,
            ufo_object.LeptonNumber,
            ufo_object.Y if hasattr(ufo_object, 'Y') else 0
        )

    @staticmethod
    def from_serializable_particle(model: Model, serializable_particle: SerializableParticle) -> Particle:
        return Particle(
            serializable_particle.pdg_code, serializable_particle.name, serializable_particle.antiname, serializable_particle.spin, serializable_particle.color,
            model.get_parameter(serializable_particle.mass),
            model.get_parameter(serializable_particle.width),
            serializable_particle.texname, serializable_particle.antitexname,
            serializable_particle.charge,
            serializable_particle.ghost_number,
            serializable_particle.lepton_number,
            serializable_particle.y_charge
        )

    def to_serializable_particle(self) -> SerializableParticle:
        return SerializableParticle.from_particle(self)

    def get_pdg_code(self) -> int:
        return self.pdg_code

    def get_anti_pdg_code(self, model: Model) -> int:
        return self.get_anti_particle(model).get_pdg_code()

    def get_anti_particle(self, model: Model) -> Particle:
        return model.get_particle(self.antiname)


class Order(object):
    def __init__(self, expansion_order: int, hierarchy: int, name: str):
        self.expansion_order: int = expansion_order
        self.hierarchy: int = hierarchy
        self.name: str = name

    @staticmethod
    def from_ufo_object(ufo_object: Any) -> Order:
        return Order(ufo_object.expansion_order, ufo_object.hierarchy, ufo_object.name)

    @staticmethod
    def from_dict(dict_repr: Any) -> Order:
        return Order(dict_repr['expansion_order'], dict_repr['hierarchy'], dict_repr['name'])


class SerializableModel(object):

    def __init__(self, name: str):
        self.name: str = name
        self.restriction: str | None = None
        self.orders: list[Order] = []
        self.parameters: list[SerializableParameter] = []
        self.particles: list[SerializableParticle] = []
        self.propagators: list[SerializablePropagator] = []
        self.lorentz_structures: list[SerializableLorentzStructure] = []
        self.couplings: list[SerializableCoupling] = []
        self.vertex_rules: list[SerializableVertexRule] = []

    @staticmethod
    def from_model(model: Model) -> SerializableModel:
        serializable_model = SerializableModel(model.name)
        serializable_model.restriction = model.restriction
        serializable_model.orders = model.orders
        serializable_model.parameters = [
            parameter.to_serializable_parameter() for parameter in model.parameters]
        serializable_model.particles = [
            particle.to_serializable_particle() for particle in model.particles]
        serializable_model.propagators = [
            propagator.to_serializable_propagator() for propagator in model.propagators]
        serializable_model.lorentz_structures = [lorentz_structure.to_serializable_lorentz_structure(
        ) for lorentz_structure in model.lorentz_structures]
        serializable_model.couplings = [
            coupling.to_serializable_coupling() for coupling in model.couplings]
        serializable_model.vertex_rules = [
            vertex_rule.to_serializable_vertex_rule() for vertex_rule in model.vertex_rules]
        return serializable_model

    def to_yaml(self) -> str:
        return utils.verbose_yaml_dump({
            'name': self.name,
            'restriction': self.restriction,
            'orders': [order.__dict__ for order in self.orders],
            'parameters': [parameter.__dict__ for parameter in self.parameters],
            'particles': [particle.__dict__ for particle in self.particles],
            'propagators': [propagator.__dict__ for propagator in self.propagators],
            'lorentz_structures': [lorentz_structure.__dict__ for lorentz_structure in self.lorentz_structures],
            'couplings': [coupling.__dict__ for coupling in self.couplings],
            'vertex_rules': [vertex_rule.__dict__ for vertex_rule in self.vertex_rules]
        })

    @staticmethod
    def from_yaml(yaml_str: str) -> SerializableModel:
        yaml_model = yaml.load(yaml_str, Loader=yaml.FullLoader)
        serializable_model = SerializableModel(yaml_model['name'])
        serializable_model.restriction = yaml_model['restriction']
        serializable_model.orders = [Order.from_dict(
            order) for order in yaml_model['orders']]
        serializable_model.parameters = [SerializableParameter.from_dict(
            parameter) for parameter in yaml_model['parameters']]
        serializable_model.particles = [SerializableParticle.from_dict(
            particle) for particle in yaml_model['particles']]
        serializable_model.propagators = [SerializablePropagator.from_dict(
            propagator) for propagator in yaml_model['propagators']]
        serializable_model.lorentz_structures = [SerializableLorentzStructure.from_dict(
            lorentz_structure) for lorentz_structure in yaml_model['lorentz_structures']]
        serializable_model.couplings = [SerializableCoupling.from_dict(
            coupling) for coupling in yaml_model['couplings']]
        serializable_model.vertex_rules = [SerializableVertexRule.from_dict(
            vertex_rule) for vertex_rule in yaml_model['vertex_rules']]
        return serializable_model


class Model(object):

    # Because parameters are now sorted accorded to their dependency,
    # it should not be necessary to go deeper than 1 level of recursion
    MAX_ALLOWED_RECURSION_IN_PARAMETER_EVALUATION = 1

    def __init__(self, name: str):
        self.name: str = name
        self.restriction: str | None = None
        self.orders: list[Order] = []
        self.parameters: list[Parameter] = []
        self.particles: list[Particle] = []
        self.propagators: list[Propagator] = []
        self.lorentz_structures: list[LorentzStructure] = []
        self.couplings: list[Coupling] = []
        self.vertex_rules: list[VertexRule] = []

        self.name_to_position: dict[str, dict[str | int, int]] = {}

    def get_full_name(self) -> str:
        return f"{self.name}-{'full' if self.restriction is None else self.restriction}"

    # We must delay the building of model functions to make sure Symbolica has been registered at this point.
    @staticmethod
    def model_function_reglog(xs: list[complex]) -> complex:
        arg: list[float] = [xs[0].real, xs[0].imag]
        if abs(arg[0]) == 0.:
            arg[0] = 0.
        if abs(arg[1]) == 0.:
            arg[1] = 0.
        if arg == [0., 0.]:
            return 0.
        else:
            return cmath.log(complex(arg[0], arg[1]))

    @staticmethod
    def model_function_reglogp(xs: list[complex]) -> complex:
        arg: list[float] = [xs[0].real, xs[0].imag]
        if abs(arg[0]) == 0.:
            arg[0] = 0.
        if abs(arg[1]) == 0.:
            arg[1] = 0.
        if arg == [0., 0.]:
            return 0.
        else:
            if arg[0] < 0. and arg[1] < 0.:
                return cmath.log(complex(arg[0], arg[1])) + complex(0., 2*cmath.pi)
            else:
                return cmath.log(complex(arg[0], arg[1]))

    @staticmethod
    def model_function_reglogm(xs: list[complex]) -> complex:
        arg: list[float] = [xs[0].real, xs[0].imag]
        if abs(arg[0]) == 0.:
            arg[0] = 0.
        if abs(arg[1]) == 0.:
            arg[1] = 0.
        if arg == [0., 0.]:
            return 0.
        else:
            if arg[0] < 0. and arg[1] > 0.:
                return cmath.log(complex(arg[0], arg[1])) - complex(0., 2*cmath.pi)
            else:
                return cmath.log(complex(arg[0], arg[1]))

    @classmethod
    def set_model_functions(cls):
        cls.model_functions: dict[Callable[[SBE], SBE], Callable[[list[complex]], complex]] = {  # type: ignore
            SBE.symbol('Theta'): lambda xs: 0. if xs[0].real < 0. else 1.,
            SBE.symbol('tan'): lambda xs: cmath.tan(xs[0]),
            SBE.symbol('acos'): lambda xs: cmath.acos(xs[0]),
            SBE.symbol('asin'): lambda xs: cmath.asin(xs[0]),
            SBE.symbol('atan'): lambda xs: cmath.atan(xs[0]),
            SBE.symbol('complexconjugate'): lambda xs: xs[0].conjugate(),
            SBE.symbol('complex'): lambda xs: xs[0]+complex(0, 1)*xs[1],
            # These functions are typically only used in NLO models
            SBE.symbol('cond'): lambda xs: xs[1] if abs(xs[0]) == 0. else xs[2],
            SBE.symbol('reglog'): lambda xs: Model.model_function_reglog(xs),
            SBE.symbol('reglogp'): lambda xs: Model.model_function_reglogp(xs),
            SBE.symbol('reglogm'): lambda xs: Model.model_function_reglogm(xs),
        }

    @classmethod
    def get_model_functions(cls) -> dict[Callable[[SBE], SBE], Callable[[list[complex]], complex]]:
        if not hasattr(cls, 'model_functions'):
            cls.set_model_functions()
        return dict(cls.model_functions)  # type: ignore

    @classmethod
    def set_model_variables(cls):
        cls.model_variables: dict[SBE, complex] = {  # type: ignore
            SBE.symbol('I'): complex(0, 1),
            SBE.symbol('pi'): cmath.pi
        }

    @classmethod
    def get_model_variables(cls) -> dict[SBE, complex]:
        if not hasattr(cls, 'model_variables'):
            cls.set_model_variables()
        return dict(cls.model_variables)  # type: ignore

    def is_empty(self) -> bool:
        return self.name == 'NotLoaded' or len(self.particles) == 0

    def get_internal_parameters(self) -> list[Parameter]:
        return [parameter for parameter in self.parameters if parameter.nature == ParameterNature.INTERNAL]

    def get_external_parameters(self) -> list[Parameter]:
        return [parameter for parameter in self.parameters if parameter.nature == ParameterNature.EXTERNAL]

    def get_coupling_orders(self) -> list[str]:
        return list(set(sum([list(coupling.orders.keys()) for coupling in self.couplings], [])))

    def sort_parameters(self):
        """ Sort parameters to linearize their evaluations, meaning all dependent parameters of parameter i appear before i """

        sorted_parameters: list[Parameter] = self.get_external_parameters()

        remaining_parameters: list[Parameter] = self.get_internal_parameters()
        for p in list(remaining_parameters):
            if p.expression is None:
                if p.value is None:
                    raise GammaLoopError(
                        "Parameter {} has no expression and no value".format(p.name))
                else:
                    sorted_parameters.append(
                        remaining_parameters.pop(remaining_parameters.index(p)))

        param_variables: dict[str, dict[str, Any]] = {
            p.name: {
                'var': SBE.symbol(p.name),
                'dependent_params': p.expression.get_all_symbols(False)
            } for p in remaining_parameters if p.expression is not None
        }


        while len(remaining_parameters) > 0:
            n_added: int = 0
            for p in list(remaining_parameters):
                if not any(
                        any(dep_param == param_variables[p.name]['var']
                            for p in remaining_parameters)
                        for dep_param in param_variables[p.name]['dependent_params']):
                    n_added += 1
                    sorted_parameters.append(
                        remaining_parameters.pop(remaining_parameters.index(p)))
            if n_added == 0:
                raise GammaLoopError(
                    "Circular dependency in parameters of model {}", self.name)

        self.parameters = sorted_parameters

    @staticmethod
    def from_ufo_model(ufo_model_path: str) -> Model:

        model_path = os.path.abspath(
            pjoin(DATA_PATH, 'models', ufo_model_path))
        sys.path.insert(0, os.path.dirname(model_path))
        sys.path.insert(0, model_path)
        model_name = os.path.basename(model_path)
        ufo_model = importlib.import_module(model_name)
        del sys.path[0]
        del sys.path[0]

        # Note that this approach does not work because of relative imports in the UFO model __init__.py file
        # spec = importlib.util.spec_from_file_location(os.path.basename(model_path), pjoin(model_path,'__init__.py'))
        # ufo_model = importlib.util.module_from_spec(spec)
        # spec.loader.exec_module(ufo_model)

        model = Model(model_name)

        # Load coupling orders
        model.orders = [Order.from_ufo_object(
            order) for order in ufo_model.all_orders]
        model.name_to_position['orders'] = {
            order.name: i for i, order in enumerate(model.orders)}
        # Load parameters
        model.parameters = [Parameter.from_ufo_object(
            param) for param in ufo_model.all_parameters]
        model.name_to_position['parameters'] = {
            param.name: i for i, param in enumerate(model.parameters)}
        if hasattr(ufo_model, 'all_CTparameters'):
            ct_parameters = [Parameter.from_ufo_object(
                param) for param in ufo_model.all_CTparameters]
            model.name_to_position['parameters'].update({
                param.name: len(model.parameters)+i for i, param in enumerate(ct_parameters)})
            model.parameters.extend(ct_parameters)

        # Load particles
        model.particles = [Particle.from_ufo_object(
            model, particle) for particle in ufo_model.all_particles]
        model.name_to_position['particles'] = {
            particle.name: i for i, particle in enumerate(model.particles)}
        model.name_to_position['particles_from_PDG'] = {
            particle.pdg_code: i for i, particle in enumerate(model.particles)}
        model.set_case_insensitive_particle_dictionary()

        # Load propagators
        model.propagators = [Propagator.from_particle(
            particle, 'Feynman') for particle in model.particles]
        model.name_to_position['propagators'] = {
            propagator.name: i for i, propagator in enumerate(model.propagators)}
        # Load Lorentz structures
        model.lorentz_structures = [LorentzStructure.from_ufo_object(
            lorentz) for lorentz in ufo_model.all_lorentz]
        model.name_to_position['lorentz_structures'] = {
            lorentz.name: i for i, lorentz in enumerate(model.lorentz_structures)}
        # Load Couplings
        model.couplings = [Coupling.from_ufo_object(
            coupling) for coupling in ufo_model.all_couplings]
        model.name_to_position['couplings'] = {
            coupling.name: i for i, coupling in enumerate(model.couplings)}
        # Load Vertices (ignore counterterms)
        model.vertex_rules = [VertexRule.from_ufo_object(
            model, vertex) for vertex in ufo_model.all_vertices if vertex.__class__.__name__ != 'CTVertex']
        model.name_to_position['vertex_rules'] = {
            vertex_rule.name: i for i, vertex_rule in enumerate(model.vertex_rules)}

        model.sort_parameters()

        return model

    def sync_name_to_position_dict(self):
        self.name_to_position.clear()
        self.name_to_position['orders'] = {
            order.name: i for i, order in enumerate(self.orders)}
        self.name_to_position['parameters'] = {
            param.name: i for i, param in enumerate(self.parameters)}
        self.name_to_position['particles'] = {
            particle.name: i for i, particle in enumerate(self.particles)}
        self.name_to_position['particles_from_PDG'] = {
            particle.pdg_code: i for i, particle in enumerate(self.particles)}
        self.set_case_insensitive_particle_dictionary()

        self.name_to_position['propagators'] = {
            propagator.name: i for i, propagator in enumerate(self.propagators)}
        self.name_to_position['lorentz_structures'] = {
            lorentz.name: i for i, lorentz in enumerate(self.lorentz_structures)}
        self.name_to_position['couplings'] = {
            coupling.name: i for i, coupling in enumerate(self.couplings)}
        self.name_to_position['vertex_rules'] = {
            vertex_rule.name: i for i, vertex_rule in enumerate(self.vertex_rules)}

    def set_case_insensitive_particle_dictionary(self):
        case_sensitivity_check: dict[str, str] = {}
        for k in self.name_to_position['particles'].keys():
            if str(k).lower() in case_sensitivity_check:
                raise GammaLoopError(
                    f"Particle {k} is case-insensitive equivalent to {case_sensitivity_check[str(k).lower()]}")
            case_sensitivity_check[str(k).lower()] = str(k)

        self.name_to_position['particles_case_insensitive'] = {
            particle.name.lower(): i for i, particle in enumerate(self.particles)}

    @staticmethod
    def from_serializable_model(serializable_model: SerializableModel) -> Model:
        model = Model(serializable_model.name)
        model.restriction = serializable_model.restriction
        model.orders = serializable_model.orders

        model.name_to_position['orders'] = {
            order.name: i for i, order in enumerate(model.orders)}
        model.parameters = [Parameter.from_serializable_parameter(
            param) for param in serializable_model.parameters]
        model.name_to_position['parameters'] = {
            param.name: i for i, param in enumerate(model.parameters)}
        model.particles = [Particle.from_serializable_particle(
            model, particle) for particle in serializable_model.particles]
        model.name_to_position['particles'] = {
            particle.name: i for i, particle in enumerate(model.particles)}
        model.name_to_position['particles_from_PDG'] = {
            particle.pdg_code: i for i, particle in enumerate(model.particles)}
        model.set_case_insensitive_particle_dictionary()

        model.propagators = [Propagator.from_serializable_propagator(
            model, propagator) for propagator in serializable_model.propagators]
        model.lorentz_structures = [LorentzStructure.from_serializable_lorentz_structure(
            lorentz) for lorentz in serializable_model.lorentz_structures]
        model.name_to_position['lorentz_structures'] = {
            lorentz.name: i for i, lorentz in enumerate(model.lorentz_structures)}
        model.couplings = [Coupling.from_serializable_coupling(
            coupling) for coupling in serializable_model.couplings]
        model.name_to_position['couplings'] = {
            coupling.name: i for i, coupling in enumerate(model.couplings)}
        model.vertex_rules = [VertexRule.from_serializable_vertex_rule(
            model, vertex_rule) for vertex_rule in serializable_model.vertex_rules]
        model.name_to_position['vertex_rules'] = {
            vertex_rule.name: i for i, vertex_rule in enumerate(model.vertex_rules)}

        return model

    def to_serializable_model(self) -> SerializableModel:
        return SerializableModel.from_model(self)

    @staticmethod
    def from_yaml(yaml_str: str) -> Model:
        return Model.from_serializable_model(SerializableModel.from_yaml(yaml_str))

    def to_yaml(self) -> str:
        return self.to_serializable_model().to_yaml()

    def get_vertices_from_particles(self, particles: list[Particle]) -> list[VertexRule]:
        return [vertex_rule for vertex_rule in self.vertex_rules if Counter(vertex_rule.particles) == Counter(particles)]

    def get_order(self, order_name: str) -> Order:
        return self.orders[self.name_to_position['orders'][order_name]]

    def get_parameter(self, parameter_name: str) -> Parameter:
        return self.parameters[self.name_to_position['parameters'][parameter_name]]

    def get_parameter_from_lha_specification(self, lhablock: str, lhacode: tuple[int, ...]) -> Parameter | None:
        for parameter in self.parameters:
            if parameter.lhablock is None:
                continue
            if parameter.lhablock.lower() == lhablock.lower() and parameter.lhacode == lhacode:
                return parameter
        return None

    def get_particle(self, particle_name: str) -> Particle:
        return self.particles[self.name_to_position['particles_case_insensitive'][particle_name.lower()]]

    def get_particle_from_pdg(self, pdg: int) -> Particle:
        return self.particles[self.name_to_position['particles_from_PDG'][pdg]]

    def get_propagator(self, propagator_name: str) -> Propagator:
        return self.propagators[self.name_to_position['propagators'][propagator_name]]

    def get_lorentz_structure(self, lorentz_name: str) -> LorentzStructure:
        return self.lorentz_structures[self.name_to_position['lorentz_structures'][lorentz_name]]

    def get_coupling(self, coupling_name: str) -> Coupling:
        return self.couplings[self.name_to_position['couplings'][coupling_name]]

    def get_vertex_rule(self, vertex_name: str) -> VertexRule:
        return self.vertex_rules[self.name_to_position['vertex_rules'][vertex_name]]

    def remove_zero_couplings(self) -> tuple[list[Coupling], list[VertexRule]]:

        removed_couplings: list[Coupling] = []
        remaining_couplings: list[Coupling] = []
        for coupling in self.couplings:
            if coupling.value == 0.:
                removed_couplings.append(coupling)
            else:
                remaining_couplings.append(coupling)
        self.couplings = remaining_couplings

        removed_vertex_rules: list[VertexRule] = []
        remaining_vertex_rules: list[VertexRule] = []
        for vertex_rule in self.vertex_rules:
            vertex_rule.couplings = {
                key: coupling
                for key, coupling in vertex_rule.couplings.items() if coupling.value != 0.
            }
            if len(vertex_rule.couplings) == 0:
                removed_vertex_rules.append(vertex_rule)
            else:
                remaining_vertex_rules.append(vertex_rule)
        self.vertex_rules = remaining_vertex_rules
        self.sync_name_to_position_dict()

        return (removed_couplings, removed_vertex_rules)

    def get_parameter_values(self) -> dict[SBE, complex]:
        parameter_map: dict[SBE, complex] = {}
        for parameter in self.parameters:
            if parameter.value is None:
                raise GammaLoopError(
                    f"The value of parameter {parameter.name} has not been set yet")
            parameter_map[SBE.symbol(parameter.name)] = parameter.value
        return parameter_map

    def update_coupling_values(self) -> None:
        parameter_values = self.get_model_variables()
        parameter_values.update(self.get_parameter_values())
        for coupling in self.couplings:
            coupling.value = utils.evaluate_symbolica_expression_safe(
                coupling.expression, parameter_values, self.get_model_functions())

    def update_internal_parameters(self, input_card: InputParamCard) -> None:
        evaluation_variables: dict[SBE, complex] = self.get_model_variables()
        for parameter in self.get_external_parameters():
            parameter.value = input_card[parameter.name.lower()]
            evaluation_variables[SBE.symbol(parameter.name)] = parameter.value

        # Collect constant internal parameters such as ZERO
        for parameter in self.get_internal_parameters():
            if parameter.expression is None:
                if parameter.value is None:
                    raise GammaLoopError(
                        "Internal parameter '{parameter.name}' has no value nor expression.")
                evaluation_variables[SBE.symbol(
                    parameter.name)] = parameter.value

        # Now update all other dependent variables
        # Note that this is done in a loop because some parameters may depend on other parameters
        evaluation_round: int = 0
        while True:
            evaluation_round += 1
            found_new_evaluation = False
            found_unevaluated = False
            for parameter in self.get_internal_parameters():
                parameter_SBE = SBE.symbol(parameter.name)
                if parameter_SBE not in evaluation_variables:
                    if parameter.expression is None:
                        raise GammaLoopError(
                            f"Internal parameter '{parameter.name}' has no value nor expression.")
                    # print([utils.expression_to_string_safe(f(SBE.symbol('x')))
                    #       for f in self.get_model_functions()])
                    # print(parameter.expression)
                    if evaluation_round == Model.MAX_ALLOWED_RECURSION_IN_PARAMETER_EVALUATION:
                        eval_result: complex | None = utils.evaluate_symbolica_expression_safe(
                            parameter.expression, evaluation_variables, self.get_model_functions())
                    else:
                        eval_result = utils.evaluate_symbolica_expression(
                            parameter.expression, evaluation_variables, self.get_model_functions())
                    if eval_result is not None:
                        parameter.value = eval_result
                        evaluation_variables[parameter_SBE] = parameter.value
                        found_new_evaluation = True
                    else:
                        found_unevaluated = True

            if not found_unevaluated:
                break
            if not found_new_evaluation:
                logger.critical(
                    "Could not evaluate all internal parameters. Reprocessing evaluation now to show latest error:")
                for parameter in self.get_internal_parameters():
                    if parameter.expression is not None:
                        eval_result = utils.evaluate_symbolica_expression_safe(
                            parameter.expression, evaluation_variables, self.get_model_functions())
                raise GammaLoopError(
                    "Could not evaluate all internal parameters.")

            if evaluation_round >= Model.MAX_ALLOWED_RECURSION_IN_PARAMETER_EVALUATION:
                raise GammaLoopError(
                    "Maximum number of allowed recursions in parameter evaluation reached.")

    def apply_input_param_card(self, input_card: InputParamCard, simplify: bool = False, update: bool = True) -> None:

        external_parameters = self.get_external_parameters()
        zero_parameters: list[str] = []
        for param in external_parameters:
            if param.name.lower() not in input_card:
                if param.value is not None:
                    logger.debug(
                        "Parameter '%s' not set in input parameter card. It will be added now with its default value of %s", param.name, param.value)
                    input_card[param.name.lower()] = param.value
                else:
                    raise GammaLoopError(
                        f"Parameter '{param.name}' is external with no default value, but not set in input parameter card.")
            else:
                param.value = input_card[param.name.lower()]

            if param.value == 0. and simplify:
                zero_parameters.append(param.name)
                param.nature = ParameterNature.INTERNAL
                param.expression = utils.parse_python_expression_safe('ZERO')
                del input_card[param.name.lower()]

        if len(zero_parameters) > 0:
            logger.debug("The following %d external parameters were forced to zero by the restriction card:\n%s", len(
                zero_parameters), ', '.join(zero_parameters))

        if update:
            self.update_internal_parameters(input_card)

            self.update_coupling_values()

            if simplify:
                removed_couplings, removed_vertices = self.remove_zero_couplings()
                if len(removed_couplings) > 0:
                    logger.debug("A total of %d couplings have been removed due restriction card specification.", len(
                        removed_couplings))
                if len(removed_vertices) > 0:
                    logger.debug("A total of %d vertices have been removed due restriction card specification:\n%s",
                                 len(removed_vertices),
                                 ', '.join(f"{v_r.name}->({'|'.join(p.name for p in v_r.particles)})" for v_r in removed_vertices))


class InputParamCard(dict[str, complex]):

    @staticmethod
    def from_param_card(param_card: ParamCard, model: Model | None = None) -> InputParamCard:
        input_param_card = InputParamCard()
        all_parameters, _restrictions = param_card.analyze_param_card(model)
        for param_name, locations in all_parameters.items():
            if len(locations) != 1:
                raise GammaLoopError(
                    f"Composite parameter defined across multiple LHA entries not supported. Parameter name: {param_name} and locations: {locations}")
            block_name, block_lhaid = locations[0]
            param_value = param_card.get_value(
                block_name, block_lhaid, 0.0)
            if not isinstance(param_value, float):
                raise GammaLoopError(
                    f"Parameter value must be float, this is not the case for parameter '{param_name}' with value '{param_value}'.")
            input_param_card[param_name] = complex(param_value, 0.)
        return input_param_card

    @staticmethod
    def from_model(model: Model) -> InputParamCard:
        input_param_card = InputParamCard()
        for param in model.get_external_parameters():
            if param.value is not None:
                input_param_card[param.name.lower()] = param.value
            else:
                raise GammaLoopError(
                    f"Model '{model.name}' cannot be loaded without any restriction card since at least one parameter, '{param.name}' lacks a default value.")
        return input_param_card

    def from_yaml(self, yaml_str: str) -> InputParamCard:
        yaml_dict = yaml.load(yaml_str, Loader=yaml.FullLoader)
        for param_name, param_value in yaml_dict.items():
            self[param_name] = param_value
        return self

    def from_yaml_file(self, yaml_path: str) -> InputParamCard:
        with open(yaml_path, 'r') as f:
            return self.from_yaml(f.read())

    def to_yaml(self) -> str:
        return utils.verbose_yaml_dump(self)
