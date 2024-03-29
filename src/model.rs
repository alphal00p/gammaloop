use crate::utils;
use ahash::RandomState;
use color_eyre::{Help, Report};
use eyre::{eyre, Context};
use num::Complex;
use serde::{Deserialize, Serialize};
use serde_yaml::Error;
use smartstring::{LazyCompact, SmartString};
use std::sync::Arc;
use std::{collections::HashMap, fs::File};
use symbolica::representations::Atom;

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub enum ParameterNature {
    #[default]
    #[serde(rename = "external")]
    External,
    #[serde(rename = "internal")]
    Internal,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub enum ParameterType {
    #[default]
    #[serde(rename = "real")]
    Real,
    #[serde(rename = "complex")]
    Imaginary,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerializableVertexRule {
    pub name: SmartString<LazyCompact>,
    pub particles: Vec<SmartString<LazyCompact>>,
    pub color_structures: Vec<SmartString<LazyCompact>>,
    pub lorentz_structures: Vec<SmartString<LazyCompact>>,
    pub couplings: Vec<Vec<Option<SmartString<LazyCompact>>>>,
}

impl SerializableVertexRule {
    pub fn from_vertex_rule(vertex_rule: &VertexRule) -> SerializableVertexRule {
        SerializableVertexRule {
            name: vertex_rule.name.clone(),
            particles: vertex_rule
                .particles
                .iter()
                .map(|particle| particle.name.clone())
                .collect(),
            color_structures: vertex_rule
                .color_structures
                .iter()
                .map(utils::to_str_expression)
                .collect(),
            lorentz_structures: vertex_rule
                .lorentz_structures
                .iter()
                .map(|lorentz_structure| lorentz_structure.name.clone())
                .collect(),
            couplings: vertex_rule
                .couplings
                .iter()
                .map(|couplings| {
                    couplings
                        .iter()
                        .map(|coupling| coupling.as_ref().map(|cpl| cpl.name.clone()))
                        .collect()
                })
                .collect(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct VertexRule {
    pub name: SmartString<LazyCompact>,
    pub particles: Vec<Arc<Particle>>,
    pub color_structures: Vec<Atom>,
    pub lorentz_structures: Vec<Arc<LorentzStructure>>,
    pub couplings: Vec<Vec<Option<Arc<Coupling>>>>,
}

impl VertexRule {
    pub fn from_serializable_vertex_rule(
        model: &Model,
        vertex_rule: &SerializableVertexRule,
    ) -> VertexRule {
        VertexRule {
            name: vertex_rule.name.clone(),
            particles: vertex_rule
                .particles
                .iter()
                .map(|particle_name| model.get_particle(particle_name).clone())
                .collect(),
            color_structures: vertex_rule
                .color_structures
                .iter()
                .map(|color_structure_name| {
                    utils::parse_python_expression(color_structure_name.as_str())
                })
                .collect(),
            lorentz_structures: vertex_rule
                .lorentz_structures
                .iter()
                .map(|lorentz_structure_name| {
                    model.get_lorentz_structure(lorentz_structure_name).clone()
                })
                .collect(),
            couplings: vertex_rule
                .couplings
                .iter()
                .map(|coupling_names| {
                    coupling_names
                        .iter()
                        .map(|coupling_name| {
                            coupling_name
                                .as_ref()
                                .map(|cpl_name| model.get_coupling(cpl_name))
                        })
                        .collect()
                })
                .collect(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerializableCoupling {
    name: SmartString<LazyCompact>,
    expression: SmartString<LazyCompact>,
    #[serde(with = "vectorize")]
    orders: HashMap<SmartString<LazyCompact>, usize, RandomState>,
    value: Option<(f64, f64)>,
}

impl SerializableCoupling {
    pub fn from_coupling(coupling: &Coupling) -> SerializableCoupling {
        SerializableCoupling {
            name: coupling.name.clone(),
            expression: utils::to_str_expression(&coupling.expression),
            orders: coupling.orders.clone(),
            value: coupling.value.map(|value| (value.re, value.im)),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Coupling {
    pub name: SmartString<LazyCompact>,
    pub expression: Atom,
    pub orders: HashMap<SmartString<LazyCompact>, usize, RandomState>,
    pub value: Option<Complex<f64>>,
}

impl Coupling {
    pub fn from_serializable_coupling(coupling: &SerializableCoupling) -> Coupling {
        Coupling {
            name: coupling.name.clone(),
            expression: utils::parse_python_expression(coupling.expression.as_str()),
            orders: coupling.orders.clone(),
            value: coupling.value.map(|value| Complex::new(value.0, value.1)),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerializableParticle {
    pdg_code: isize,
    name: SmartString<LazyCompact>,
    antiname: SmartString<LazyCompact>,
    spin: isize,
    color: isize,
    mass: SmartString<LazyCompact>,
    width: SmartString<LazyCompact>,
    texname: SmartString<LazyCompact>,
    antitexname: SmartString<LazyCompact>,
    charge: f64,
    ghost_number: isize,
    lepton_number: isize,
    y_charge: isize,
}

impl SerializableParticle {
    pub fn from_particle(particle: &Particle) -> SerializableParticle {
        SerializableParticle {
            pdg_code: particle.pdg_code,
            name: particle.name.clone(),
            antiname: particle.antiname.clone(),
            spin: particle.spin,
            color: particle.color,
            mass: particle.mass.name.clone(),
            width: particle.width.name.clone(),
            texname: particle.texname.clone(),
            antitexname: particle.antitexname.clone(),
            charge: particle.charge,
            ghost_number: particle.ghost_number,
            lepton_number: particle.lepton_number,
            y_charge: particle.y_charge,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Particle {
    pub pdg_code: isize,
    pub name: SmartString<LazyCompact>,
    pub antiname: SmartString<LazyCompact>,
    pub spin: isize,
    pub color: isize,
    pub mass: Arc<Parameter>,
    pub width: Arc<Parameter>,
    pub texname: SmartString<LazyCompact>,
    pub antitexname: SmartString<LazyCompact>,
    pub charge: f64,
    pub ghost_number: isize,
    pub lepton_number: isize,
    pub y_charge: isize,
}

impl Particle {
    pub fn from_serializable_particle(model: &Model, particle: &SerializableParticle) -> Particle {
        Particle {
            pdg_code: particle.pdg_code,
            name: particle.name.clone(),
            antiname: particle.antiname.clone(),
            spin: particle.spin,
            color: particle.color,
            mass: model.get_parameter(&particle.mass),
            width: model.get_parameter(&particle.width),
            texname: particle.texname.clone(),
            antitexname: particle.antitexname.clone(),
            charge: particle.charge,
            ghost_number: particle.ghost_number,
            lepton_number: particle.lepton_number,
            y_charge: particle.y_charge,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerializableLorentzStructure {
    name: SmartString<LazyCompact>,
    spins: Vec<isize>,
    structure: SmartString<LazyCompact>,
}

impl SerializableLorentzStructure {
    pub fn from_lorentz_structure(ls: &LorentzStructure) -> SerializableLorentzStructure {
        SerializableLorentzStructure {
            name: ls.name.clone(),
            spins: ls.spins.clone(),
            structure: utils::to_str_expression(&ls.structure),
        }
    }
}

#[derive(Debug, Clone)]
pub struct LorentzStructure {
    pub name: SmartString<LazyCompact>,
    pub spins: Vec<isize>,
    pub structure: Atom,
}

impl LorentzStructure {
    pub fn from_serializable_lorentz_structure(
        ls: &SerializableLorentzStructure,
    ) -> LorentzStructure {
        LorentzStructure {
            name: ls.name.clone(),
            spins: ls.spins.clone(),
            structure: utils::parse_python_expression(ls.structure.as_str()),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerializableParameter {
    name: SmartString<LazyCompact>,
    lhablock: Option<SmartString<LazyCompact>>,
    lhacode: Option<Vec<usize>>,
    nature: ParameterNature,
    parameter_type: ParameterType,
    value: Option<(f64, f64)>,
    expression: Option<SmartString<LazyCompact>>,
}

impl SerializableParameter {
    pub fn from_parameter(param: &Parameter) -> SerializableParameter {
        SerializableParameter {
            name: param.name.clone(),
            lhablock: param.lhablock.clone(),
            lhacode: param.lhacode.clone(),
            nature: param.nature.clone(),
            parameter_type: param.parameter_type.clone(),
            value: param.value.map(|value| (value.re, value.im)),
            expression: param.expression.as_ref().map(utils::to_str_expression),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Parameter {
    pub name: SmartString<LazyCompact>,
    pub lhablock: Option<SmartString<LazyCompact>>,
    pub lhacode: Option<Vec<usize>>,
    pub nature: ParameterNature,
    pub parameter_type: ParameterType,
    pub value: Option<Complex<f64>>,
    pub expression: Option<Atom>,
}

impl Parameter {
    pub fn from_serializable_parameter(param: &SerializableParameter) -> Parameter {
        Parameter {
            name: param.name.clone(),
            lhablock: param.lhablock.clone(),
            lhacode: param.lhacode.clone(),
            nature: param.nature.clone(),
            parameter_type: param.parameter_type.clone(),
            value: param.value.map(|value| Complex::new(value.0, value.1)),
            expression: param
                .expression
                .as_ref()
                .map(|expr| utils::parse_python_expression(expr.as_str())),
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Order {
    pub name: SmartString<LazyCompact>,
    pub expansion_order: isize,
    pub hierarchy: isize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerializableModel {
    pub name: SmartString<LazyCompact>,
    pub restriction: Option<SmartString<LazyCompact>>,
    orders: Vec<Order>,
    parameters: Vec<SerializableParameter>,
    particles: Vec<SerializableParticle>,
    lorentz_structures: Vec<SerializableLorentzStructure>,
    couplings: Vec<SerializableCoupling>,
    vertex_rules: Vec<SerializableVertexRule>,
}

impl SerializableModel {
    pub fn from_file(file_path: String) -> Result<SerializableModel, Report> {
        let f = File::open(file_path.clone())
            .wrap_err_with(|| format!("Could not open model yaml file {}", file_path))
            .suggestion("Does the path exist?")?;
        serde_yaml::from_reader(f)
            .map_err(|e| eyre!(format!("Error parsing model yaml: {}", e)))
            .suggestion("Is it a correct yaml file")
    }

    pub fn from_yaml_str(yaml_str: String) -> Result<SerializableModel, Report> {
        serde_yaml::from_str(yaml_str.as_str())
            .map_err(|e| eyre!(format!("Error parsing model yaml: {}", e)))
            .suggestion("Is it a correct yaml file")
    }

    pub fn from_model(model: &Model) -> SerializableModel {
        SerializableModel {
            name: model.name.clone(),
            restriction: model.restriction.clone(),
            orders: model
                .orders
                .iter()
                .map(|order| order.as_ref().clone())
                .collect(),
            parameters: model
                .parameters
                .iter()
                .map(|parameter| SerializableParameter::from_parameter(parameter.as_ref()))
                .collect(),
            particles: model
                .particles
                .iter()
                .map(|particle| SerializableParticle::from_particle(particle.as_ref()))
                .collect(),
            lorentz_structures: model
                .lorentz_structures
                .iter()
                .map(|lorentz_structure| {
                    SerializableLorentzStructure::from_lorentz_structure(lorentz_structure.as_ref())
                })
                .collect(),
            couplings: model
                .couplings
                .iter()
                .map(|coupling| SerializableCoupling::from_coupling(coupling.as_ref()))
                .collect(),
            vertex_rules: model
                .vertex_rules
                .iter()
                .map(|vertex_rule| SerializableVertexRule::from_vertex_rule(vertex_rule.as_ref()))
                .collect(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Model {
    pub name: SmartString<LazyCompact>,
    pub restriction: Option<SmartString<LazyCompact>>,
    pub orders: Vec<Arc<Order>>,
    pub parameters: Vec<Arc<Parameter>>,
    pub particles: Vec<Arc<Particle>>,
    pub lorentz_structures: Vec<Arc<LorentzStructure>>,
    pub couplings: Vec<Arc<Coupling>>,
    pub vertex_rules: Vec<Arc<VertexRule>>,
    pub order_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState>,
    pub parameter_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState>,
    pub lorentz_structure_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState>,
    pub particle_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState>,
    pub particle_pdg_to_position: HashMap<isize, usize, RandomState>,
    pub coupling_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState>,
    pub vertex_rule_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState>,
}

impl Default for Model {
    fn default() -> Self {
        Model {
            name: SmartString::<LazyCompact>::from("ModelNotLoaded"),
            restriction: None,
            orders: vec![],
            parameters: vec![],
            particles: vec![],
            lorentz_structures: vec![],
            couplings: vec![],
            vertex_rules: vec![],
            order_name_to_position:
                HashMap::<SmartString<LazyCompact>, usize, RandomState>::default(),
            parameter_name_to_position:
                HashMap::<SmartString<LazyCompact>, usize, RandomState>::default(),
            lorentz_structure_name_to_position: HashMap::<
                SmartString<LazyCompact>,
                usize,
                RandomState,
            >::default(),
            particle_name_to_position:
                HashMap::<SmartString<LazyCompact>, usize, RandomState>::default(),
            particle_pdg_to_position: HashMap::<isize, usize, RandomState>::default(),
            coupling_name_to_position:
                HashMap::<SmartString<LazyCompact>, usize, RandomState>::default(),
            vertex_rule_name_to_position:
                HashMap::<SmartString<LazyCompact>, usize, RandomState>::default(),
        }
    }
}
impl Model {
    pub fn is_empty(&self) -> bool {
        self.name == "ModelNotLoaded" || self.particles.is_empty()
    }

    pub fn from_serializable_model(serializable_model: SerializableModel) -> Model {
        let mut model: Model = Model::default();
        model.name = serializable_model.name;
        model.restriction = serializable_model.restriction;

        // Extract coupling orders
        model.orders = serializable_model
            .orders
            .iter()
            .enumerate()
            .map(|(i_order, serializable_order)| {
                let order = Arc::new(Order {
                    name: serializable_order.name.clone(),
                    expansion_order: serializable_order.expansion_order,
                    hierarchy: serializable_order.hierarchy,
                });
                model
                    .order_name_to_position
                    .insert(order.name.clone(), i_order);
                order
            })
            .collect();

        // Extract parameters
        model.parameters = serializable_model
            .parameters
            .iter()
            .enumerate()
            .map(|(i_param, serializable_param)| {
                let parameter =
                    Arc::new(Parameter::from_serializable_parameter(serializable_param));
                model
                    .parameter_name_to_position
                    .insert(parameter.name.clone(), i_param);
                parameter
            })
            .collect();

        // Extract particles
        model.particles = serializable_model
            .particles
            .iter()
            .enumerate()
            .map(|(i_part, serializable_particle)| {
                let particle = Arc::new(Particle::from_serializable_particle(
                    &model,
                    serializable_particle,
                ));
                model
                    .particle_name_to_position
                    .insert(particle.name.clone(), i_part);
                model
                    .particle_pdg_to_position
                    .insert(particle.pdg_code, i_part);
                particle
            })
            .collect();

        // Extract Lorentz structures
        model.lorentz_structures = serializable_model
            .lorentz_structures
            .iter()
            .enumerate()
            .map(|(i_lor, serializable_lorentz_structure)| {
                let lorentz_structure =
                    Arc::new(LorentzStructure::from_serializable_lorentz_structure(
                        serializable_lorentz_structure,
                    ));
                model
                    .lorentz_structure_name_to_position
                    .insert(lorentz_structure.name.clone(), i_lor);
                lorentz_structure
            })
            .collect();

        // Extract couplings
        model.couplings = serializable_model
            .couplings
            .iter()
            .enumerate()
            .map(|(i_coupl, serializable_coupling)| {
                let coupling =
                    Arc::new(Coupling::from_serializable_coupling(serializable_coupling));
                model
                    .coupling_name_to_position
                    .insert(coupling.name.clone(), i_coupl);
                coupling
            })
            .collect();

        // Extract vertex rules
        model.vertex_rules = serializable_model
            .vertex_rules
            .iter()
            .enumerate()
            .map(|(i_vr, serializable_vertex_rule)| {
                let vertex_rule = Arc::new(VertexRule::from_serializable_vertex_rule(
                    &model,
                    serializable_vertex_rule,
                ));
                model
                    .vertex_rule_name_to_position
                    .insert(vertex_rule.name.clone(), i_vr);
                vertex_rule
            })
            .collect();

        model
    }

    pub fn to_serializable(&self) -> SerializableModel {
        SerializableModel::from_model(self)
    }

    pub fn to_yaml(&self) -> Result<String, Error> {
        serde_yaml::to_string(&self.to_serializable())
    }

    pub fn from_file(file_path: String) -> Result<Model, Report> {
        SerializableModel::from_file(file_path).map(Model::from_serializable_model)
    }

    pub fn from_yaml_str(yaml_str: String) -> Result<Model, Report> {
        SerializableModel::from_yaml_str(yaml_str).map(Model::from_serializable_model)
    }

    #[inline]
    pub fn get_particle(&self, name: &SmartString<LazyCompact>) -> Arc<Particle> {
        if let Some(position) = self.particle_name_to_position.get(name) {
            self.particles[*position].clone()
        } else {
            panic!("Particle '{}' not found in model '{}'.", name, self.name);
        }
    }
    #[inline]
    pub fn get_particle_from_pdg(&self, pdg: isize) -> Arc<Particle> {
        if let Some(position) = self.particle_pdg_to_position.get(&pdg) {
            self.particles[*position].clone()
        } else {
            panic!(
                "Particle with PDG {} not found in model '{}'.",
                pdg, self.name
            );
        }
    }
    #[inline]
    pub fn get_parameter(&self, name: &SmartString<LazyCompact>) -> Arc<Parameter> {
        if let Some(position) = self.parameter_name_to_position.get(name) {
            self.parameters[*position].clone()
        } else {
            panic!("Parameter '{}' not found in model '{}'.", name, self.name);
        }
    }
    #[inline]
    pub fn get_order(&self, name: &SmartString<LazyCompact>) -> Arc<Order> {
        if let Some(position) = self.order_name_to_position.get(name) {
            self.orders[*position].clone()
        } else {
            panic!(
                "Coupling order '{}' not found in model '{}'.",
                name, self.name
            );
        }
    }
    #[inline]
    pub fn get_lorentz_structure(&self, name: &SmartString<LazyCompact>) -> Arc<LorentzStructure> {
        if let Some(position) = self.lorentz_structure_name_to_position.get(name) {
            self.lorentz_structures[*position].clone()
        } else {
            panic!(
                "Lorentz structure '{}' not found in model '{}'.",
                name, self.name
            );
        }
    }
    #[inline]
    pub fn get_coupling(&self, name: &SmartString<LazyCompact>) -> Arc<Coupling> {
        if let Some(position) = self.coupling_name_to_position.get(name) {
            self.couplings[*position].clone()
        } else {
            panic!("Coupling '{}' not found in model '{}'.", name, self.name);
        }
    }
    #[inline]
    pub fn get_vertex_rule(&self, name: &SmartString<LazyCompact>) -> Arc<VertexRule> {
        if let Some(position) = self.vertex_rule_name_to_position.get(name) {
            self.vertex_rules[*position].clone()
        } else {
            panic!("Vertex rule '{}' not found in model '{}'.", name, self.name);
        }
    }
}
