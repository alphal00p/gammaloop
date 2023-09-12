from __future__ import annotations

import importlib
import importlib.util
import os
import sys
import yaml
from pprint import pprint
from gammaloop.misc.common import *
import symbolica as sb
import gammaloop.misc.utils as utils
from symbolica import Expression as SBE
pjoin = os.path.join

class SerializableVertexRule(object):
    def __init__(self, name: str, particles: list[str], color_structures: list[str], lorentz_structures: list[str], couplings: list[list[str]]):
        self.name: str = name
        self.particles: list[str] = particles
        self.color_structures : list[str] = color_structures
        self.lorentz_structures : list[str] = lorentz_structures
        self.couplings: list[list[str]] = couplings
    
    @staticmethod
    def from_vertex_rule(vertex_rule: VertexRule) -> SerializableVertexRule:
        return SerializableVertexRule(
            vertex_rule.name, 
            [particle.name for particle in vertex_rule.particles], 
            [utils.expression_to_string(c) for c in vertex_rule.color_structures],
            [lorentz.name for lorentz in vertex_rule.lorentz_structures],
            [
                [
                    None if (i_col, j_lor) not in vertex_rule.couplings else vertex_rule.couplings[(i_col, j_lor)].name 
                      for j_lor in range(len(vertex_rule.lorentz_structures))
                ] for i_col in range(len(vertex_rule.color_structures))
            ]
        )
    
    @staticmethod
    def from_dict(dict) -> SerializableVertexRule:
        return SerializableVertexRule(
            dict['name'], 
            dict['particles'], 
            dict['color_structures'],
            dict['lorentz_structures'],
            dict['couplings']
        )

class VertexRule(object):
    def __init__(self, name: str, particles: list[Particle], color_structures: list[str], lorentz_structures: list[LorentzStructure], couplings: dict[(int,int),Coupling]):
        self.name: str = name
        self.particles: list[Particle] = particles
        self.color_structures : list[SBE] = [utils.parse_python_expression(color_structure) for color_structure  in color_structures]
        self.lorentz_structures : list[LorentzStructure] = lorentz_structures
        self.couplings: dict[(int,int),Coupling] = couplings

    @staticmethod
    def from_UFO_object(model: Model, ufo_object) -> VertexRule:
        return VertexRule(
            ufo_object.name, 
            [model.get_particle(particle.name) for particle in ufo_object.particles], 
            ufo_object.color,
            [model.get_lorentz_structure(lorentz.name) for lorentz in ufo_object.lorentz],
            {(i,j): model.get_coupling(coupling.name) for (i,j), coupling in ufo_object.couplings.items()}
        )

    @staticmethod
    def from_serializable_vertex_rule(model: Model, serializable_vertex_rule: SerializableVertexRule) -> VertexRule:
        return VertexRule(
            serializable_vertex_rule.name, 
            [model.get_particle(particle_name) for particle_name in serializable_vertex_rule.particles], 
            serializable_vertex_rule.color_structures,
            [model.get_lorentz_structure(lorentz_name) for lorentz_name in serializable_vertex_rule.lorentz_structures],
            dict( ((i, j), model.get_coupling(serializable_vertex_rule.couplings[i][j])) 
                for i in range(len(serializable_vertex_rule.color_structures)) 
                for j in range(len(serializable_vertex_rule.lorentz_structures)) if serializable_vertex_rule.couplings[i][j] is not None )
        )
    
    def to_serializable_vertex_rule(self) -> SerializableVertexRule:
        return SerializableVertexRule.from_vertex_rule(self)

class SerializableCoupling(object):
    def __init__(self, name: str, expression: str, orders: list[(str,int)]):
        self.name: str = name
        self.expression: str = expression
        self.orders : list[(str,int)] = orders

    @staticmethod
    def from_coupling(coupling: Coupling) -> SerializableCoupling:
        return SerializableCoupling(
            coupling.name, 
            utils.expression_to_string(coupling.expression), 
            list(coupling.orders.items())
        )
    
    @staticmethod
    def from_dict(dict) -> SerializableCoupling:
        return SerializableCoupling(
            dict['name'], 
            dict['expression'], 
            dict['orders']
        )

class Coupling(object):
    def __init__(self, name: str, expression: str, orders: dict[str,int]):
        self.name: str = name
        self.expression: SBE = utils.parse_python_expression(expression)
        self.orders : dict[str,int] = orders

    @staticmethod
    def from_UFO_object(ufo_object) -> Coupling:
        return Coupling(ufo_object.name, ufo_object.value, ufo_object.order)

    @staticmethod
    def from_serializable_coupling(serializable_coupling: SerializableCoupling) -> Coupling:
        return Coupling(
            serializable_coupling.name, 
            serializable_coupling.expression, 
            dict(serializable_coupling.orders)
        )

    def to_serializable_coupling(self) -> SerializableCoupling:
        return SerializableCoupling.from_coupling(self)

class SerializableLorentzStructure(object):
    def __init__(self, name: str, spins: list[int], structure: str):
        self.name: str = name
        self.spins: list[int] = []
        self.structure: str = structure
    
    @staticmethod
    def from_lorentz_structure(lorentz_structure: LorentzStructure) -> SerializableLorentzStructure:
        return SerializableLorentzStructure(
            lorentz_structure.name, 
            lorentz_structure.spins, 
            utils.expression_to_string(lorentz_structure.structure)
        )

    @staticmethod
    def from_dict(dict) -> SerializableLorentzStructure:
        return SerializableLorentzStructure(
            dict['name'], 
            dict['spins'], 
            dict['structure']
        )

class LorentzStructure(object):
    def __init__(self, name: str, spins: list[int], structure: str):
        self.name: str = name
        self.spins: list[int] = []
        self.structure: SBE = utils.parse_python_expression(structure)

    @staticmethod
    def from_UFO_object(ufo_object) -> LorentzStructure:
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
    def __init__(self, lhablock: str | None, lhacode: int | None, name: str, nature: str, parameter_type: str, value: float | None, expression: str | None):
        self.lhablock: str | None = lhablock
        self.lhacode: list[int] | None = lhacode
        self.name: str = name
        self.nature: str = nature
        self.parameter_type: str = parameter_type
        self.value: float | None = value
        self.expression: str | None = expression

    @staticmethod
    def from_parameter(parameter: Parameter) -> SerializableParameter:
        return SerializableParameter(
            parameter.lhablock, parameter.lhacode, parameter.name, 
            parameter.nature, parameter.parameter_type, parameter.value, 
            utils.expression_to_string(parameter.expression)
        )

    @staticmethod
    def from_dict(dict) -> SerializableParameter:
        return SerializableParameter(
            dict['lhablock'], dict['lhacode'], dict['name'], 
            dict['nature'], dict['parameter_type'], dict['value'], 
            dict['expression']
        )

class Parameter(object):
    def __init__(self, lhablock: str | None, lhacode: int | None, name: str, nature: str, parameter_type: str, value: float | None, expression: str | None):
        self.lhablock: str | None = lhablock
        self.lhacode: list[int] | None = lhacode
        self.name: str = name
        self.nature: str = nature
        self.parameter_type: str = parameter_type
        self.value: float | None = value
        self.expression: SBE | None = utils.parse_python_expression(expression) 

    @staticmethod
    def from_UFO_object(ufo_object) -> Parameter:
        if ufo_object.nature == 'external':
            param_value = float(ufo_object.value)
            param_expression : str | None = None
            match param_value:
                case 0.: param_expression = '0'
                case 1.: param_expression = '1'
                case _: param_expression = None
        elif ufo_object.nature == 'internal':
            param_expression = ufo_object.value
            try:
                param_value = float(ufo_object.value)
            except:
                param_value = None
        else:
            raise GammaLoopError(f"Invalid parameter nature '{ufo_object.nature}'")

        return Parameter(
            ufo_object.lhablock, ufo_object.lhacode, ufo_object.name, ufo_object.nature, ufo_object.type, param_value, param_expression
        )
    
    @staticmethod
    def from_serializable_parameter(serializable_parameter: SerializableParameter) -> Parameter:
        return Parameter(
            serializable_parameter.lhablock, serializable_parameter.lhacode, serializable_parameter.name, 
            serializable_parameter.nature, serializable_parameter.parameter_type, serializable_parameter.value, serializable_parameter.expression
        )
    
    def to_serializable_parameter(self) -> SerializableParameter:
        return SerializableParameter.from_parameter(self)

    def __str__(self) -> str:
        if self.nature == 'internal':
            return f"Internal parameter({self.name}, {self.expression})"
        else:
            return f"External parameter({self.name}, {self.value})"

class SerializableParticle(object):
    def __init__(self, pdg_code: int, name: str, antiname: str, spin: int, color: int, mass: str, width: str, texname: str, antitexname: str, charge: int, ghost_number: int, lepton_number: int, Y: int):
        self.pdg_code: int = pdg_code
        self.name: str = name
        self.antiname: str = antiname
        self.spin: int = 3
        self.color: int = 1
        self.mass: str = mass
        self.width: str = width
        self.texname: str = texname
        self.antitexname: str = antitexname
        self.charge: float = charge
        self.ghost_number: int = ghost_number
        self.lepton_number: int = lepton_number
        self.y: int = Y

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
            particle.y
        )

    @classmethod
    def from_dict(cls, dict) -> SerializableParticle:
        return SerializableParticle(
            dict['pdg_code'], dict['name'], dict['antiname'], dict['spin'], dict['color'], 
            dict['mass'], 
            dict['width'],
            dict['texname'], dict['antitexname'], 
            dict['charge'], 
            dict['ghost_number'], 
            dict['lepton_number'], 
            dict['y']
        )

class Particle(object):
    def __init__(self, pdg_code: int, name: str, antiname: str, spin: int, color: int, mass: Parameter, width: Parameter, texname: str, antitexname: str, charge: int, ghost_number: int, lepton_number: int, Y: int):
        self.pdg_code: int = pdg_code
        self.name: str = name
        self.antiname: str = antiname
        self.spin: int = 3
        self.color: int = 1
        self.mass: Parameter = mass
        self.width: Parameter = width
        self.texname: str = texname
        self.antitexname: str = antitexname
        self.charge: float = charge
        self.ghost_number: int = ghost_number
        self.lepton_number: int = lepton_number
        self.y: int = Y

    @staticmethod
    def from_UFO_object(model: Model, ufo_object) -> Particle:

        return Particle(
            ufo_object.pdg_code, ufo_object.name, ufo_object.antiname, ufo_object.spin, ufo_object.color, 
            model.get_parameter(ufo_object.mass.name),
            model.get_parameter(ufo_object.width.name), 
            ufo_object.texname, ufo_object.antitexname, 
            ufo_object.charge, 
            ufo_object.GhostNumber, 
            ufo_object.LeptonNumber, 
            ufo_object.Y
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
            serializable_particle.y
        )
    
    def to_serializable_particle(self) -> SerializableParticle:
        return SerializableParticle.from_particle(self)

class Order(object):
    def __init__(self, expansion_order: int, hierarchy: int, name: str):
        self.expansion_order: int = expansion_order
        self.hierarchy: int = hierarchy
        self.name: str = name

    @staticmethod
    def from_UFO_object(ufo_object) -> Order:
        return Order(ufo_object.expansion_order, ufo_object.hierarchy, ufo_object.name)

    @staticmethod
    def from_dict(dict) -> Order:
        return Order(dict['expansion_order'], dict['hierarchy'], dict['name'])

class SerializableModel(object):

    def __init__(self, name: str):
        self.name: str = name
        self.orders: list[Order] = []
        self.parameters: list[SerializableParameter] = []
        self.particles: list[SerializableParticle] = []
        self.lorentz_structures: list[SerializableLorentzStructure] = []
        self.couplings: list[SerializableCoupling] = []
        self.vertex_rules: list[SerializableVertexRule] = []
    
    @staticmethod
    def from_model(model: Model) -> SerializableModel:
        serializable_model = SerializableModel(model.name)
        serializable_model.orders = model.orders
        serializable_model.parameters = [ parameter.to_serializable_parameter() for parameter in model.parameters ]
        serializable_model.particles = [ particle.to_serializable_particle() for particle in model.particles ]
        serializable_model.lorentz_structures = [ lorentz_structure.to_serializable_lorentz_structure() for lorentz_structure in model.lorentz_structures ]
        serializable_model.couplings = [ coupling.to_serializable_coupling() for coupling in model.couplings ]
        serializable_model.vertex_rules = [ vertex_rule.to_serializable_vertex_rule() for vertex_rule in model.vertex_rules ]
        return serializable_model

    def to_yaml(self) -> str:
        return utils.verbose_yaml_dump({
            'name': self.name,
            'orders': [ order.__dict__ for order in self.orders ],
            'parameters': [ parameter.__dict__ for parameter in self.parameters ],
            'particles': [ particle.__dict__ for particle in self.particles ],
            'lorentz_structures': [ lorentz_structure.__dict__ for lorentz_structure in self.lorentz_structures ],
            'couplings': [ coupling.__dict__ for coupling in self.couplings ],
            'vertex_rules': [ vertex_rule.__dict__ for vertex_rule in self.vertex_rules ]
        })

    @staticmethod
    def from_yaml(yaml_str) -> SerializableModel:
        yaml_model = yaml.load(yaml_str, Loader=yaml.FullLoader)
        serializable_model = SerializableModel(yaml_model['name'])
        serializable_model.orders = [ Order.from_dict(order) for order in yaml_model['orders'] ]
        serializable_model.parameters = [ SerializableParameter.from_dict(parameter) for parameter in yaml_model['parameters'] ]
        serializable_model.particles = [ SerializableParticle.from_dict(particle) for particle in yaml_model['particles'] ]
        serializable_model.lorentz_structures = [ SerializableLorentzStructure.from_dict(lorentz_structure) for lorentz_structure in yaml_model['lorentz_structures'] ]
        serializable_model.couplings = [ SerializableCoupling.from_dict(coupling) for coupling in yaml_model['couplings'] ]
        serializable_model.vertex_rules = [ SerializableVertexRule.from_dict(vertex_rule) for vertex_rule in yaml_model['vertex_rules'] ]
        return serializable_model


class Model(object):

    def __init__(self, name: str):
        self.name: str = name
        self.orders: list[Order] = []
        self.parameters: list[Parameter] = []
        self.particles: list[Particle] = []
        self.lorentz_structures: list[LorentzStructure] = []
        self.couplings: list[Coupling] = []
        self.vertex_rules: list[VertexRule] = []

        self.name_to_position: dict[str,dict[str | int,int]] = {}

    @staticmethod
    def from_UFOModel(ufo_model_path: str) -> Model:

        model_path = os.path.abspath(pjoin(DATA_PATH,'models',ufo_model_path))
        sys.path.insert(0,os.path.dirname(model_path))
        model_name = os.path.basename(model_path)
        ufo_model = importlib.import_module(model_name)
        del sys.path[0]
        
        # Note that this approach does not work because of relative imports in the UFO model __init__.py file
        #spec = importlib.util.spec_from_file_location(os.path.basename(model_path), pjoin(model_path,'__init__.py'))
        #ufo_model = importlib.util.module_from_spec(spec)
        #spec.loader.exec_module(ufo_model)

        model = Model(model_name)

        # Load coupling orders
        model.orders = [ Order.from_UFO_object(order) for order in ufo_model.all_orders ]
        model.name_to_position['orders'] = { order.name: i for i, order in enumerate(model.orders) }
        # Load parameters
        model.parameters = [ Parameter.from_UFO_object(param) for param in ufo_model.all_parameters ]
        model.name_to_position['parameters'] = { param.name: i for i, param in enumerate(model.parameters) }
        # Load particles
        model.particles = [ Particle.from_UFO_object(model, particle) for particle in ufo_model.all_particles ]
        model.name_to_position['particles'] = { particle.name: i for i, particle in enumerate(model.particles) }
        model.name_to_position['particles_from_PDG'] = { particle.pdg_code: i for i, particle in enumerate(model.particles) }
        # Load Lorentz structures
        model.lorentz_structures = [ LorentzStructure.from_UFO_object(lorentz) for lorentz in ufo_model.all_lorentz ]
        model.name_to_position['lorentz_structures'] = { lorentz.name: i for i, lorentz in enumerate(model.lorentz_structures) }
        # Load Couplings
        model.couplings = [ Coupling.from_UFO_object(coupling) for coupling in ufo_model.all_couplings ]
        model.name_to_position['couplings'] = { coupling.name: i for i, coupling in enumerate(model.couplings) }
        # Load Vertices (ignore counterterms)
        model.vertex_rules = [ VertexRule.from_UFO_object(model, vertex) for vertex in ufo_model.all_vertices if vertex.__class__.__name__ != 'CTVertex' ]
        model.name_to_position['vertex_rules'] = { vertex_rule.name: i for i, vertex_rule in enumerate(model.vertex_rules) }
        
        return model

    @staticmethod
    def from_serializable_model(serializable_model: SerializableModel) -> Model:
        model = Model(serializable_model.name)
        model.orders = serializable_model.orders
        print(model.orders)
        model.name_to_position['orders'] = { order.name: i for i, order in enumerate(model.orders) }
        model.parameters = [ Parameter.from_serializable_parameter(param) for param in serializable_model.parameters ]
        model.name_to_position['parameters'] = { param.name: i for i, param in enumerate(model.parameters) }
        model.particles = [ Particle.from_serializable_particle(model, particle) for particle in serializable_model.particles ]
        model.name_to_position['particles'] = { particle.name: i for i, particle in enumerate(model.particles) }
        model.name_to_position['particles_from_PDG'] = { particle.pdg_code: i for i, particle in enumerate(model.particles) }
        model.lorentz_structures = [ LorentzStructure.from_serializable_lorentz_structure(lorentz) for lorentz in serializable_model.lorentz_structures ]
        model.name_to_position['lorentz_structures'] = { lorentz.name: i for i, lorentz in enumerate(model.lorentz_structures) }
        model.couplings = [ Coupling.from_serializable_coupling(coupling) for coupling in serializable_model.couplings ]
        model.name_to_position['couplings'] = { coupling.name: i for i, coupling in enumerate(model.couplings) }
        model.vertex_rules = [ VertexRule.from_serializable_vertex_rule(model, vertex_rule) for vertex_rule in serializable_model.vertex_rules ]
        model.name_to_position['vertex_rules'] = { vertex_rule.name: i for i, vertex_rule in enumerate(model.vertex_rules) }
        return model

    def to_serializable_model(self) -> SerializableModel:
        return SerializableModel.from_model(self)

    @staticmethod
    def from_yaml(yaml_str: str) -> Model:
        return Model.from_serializable_model(SerializableModel.from_yaml(yaml_str))
    
    def to_yaml(self) -> str:
        return self.to_serializable_model().to_yaml()

    def get_order(self, order_name: str) -> Order:
        return self.orders[self.name_to_position['orders'][order_name]]
    
    def get_parameter(self, parameter_name: str) -> Parameter:
        return self.parameters[self.name_to_position['parameters'][parameter_name]]

    def get_particle(self, particle_name: str) -> Particle:
        return self.particles[self.name_to_position['particles'][particle_name]]

    def get_particle_from_pdg(self, pdg: int) -> Particle:
        return self.particles[self.name_to_position['particles_from_PDG'][pdg]]

    def get_lorentz_structure(self, lorentz_name: str) -> LorentzStructure:
        return self.lorentz_structures[self.name_to_position['lorentz_structures'][lorentz_name]]

    def get_coupling(self, coupling_name: str) -> Coupling:
        return self.couplings[self.name_to_position['couplings'][coupling_name]]

    def get_vertex_rule(self, vertex_name: str) -> VertexRule:
        return self.vertex_rules[self.name_to_position['vertex_rules'][vertex_name]]