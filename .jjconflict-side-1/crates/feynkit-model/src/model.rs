use std::{
    collections::{BTreeMap, BTreeSet, HashMap},
    fs,
    path::Path,
    sync::OnceLock,
};

use serde::{Deserialize, Deserializer, Serialize, Serializer, de};
use symbolica::{
    atom::{Atom, AtomCore},
    id::Replacement,
    parser::ParseSettings,
    printer::PrintOptions,
    symbol,
};

use crate::{
    ComplexValue, EntityKind, EvaluatedValues, EvaluationRequest, ModelError, ModelEvaluator,
    ModelExpression, ModelValidationError, ParameterCard, RecomputeError,
};

const fn default_true() -> bool {
    true
}

fn parse_expression(
    kind: EntityKind,
    name: &str,
    field: &'static str,
    expression: &str,
) -> Result<Atom, ModelError> {
    let normalized = expression
        .replace("**", "^")
        .replace("cmath.sqrt", "sqrt")
        .replace("cmath.pi", "pi")
        .replace("math.sqrt", "sqrt")
        .replace("math.pi", "pi");
    Atom::parse(&normalized, "UFO", ParseSettings::default()).map_err(|message| {
        ModelError::SymbolicParse {
            kind,
            name: name.to_owned(),
            field,
            expression: expression.to_owned(),
            message,
        }
    })
}

/// Export an expression in the UFO interchange syntax.
///
/// Runtime atoms deliberately live in Symbolica's `UFO` namespace, but the
/// UFO JSON and evaluator boundaries use the original, namespaceless symbol
/// spelling. Keeping that conversion in one place also gives every exported
/// expression a deterministic ordering.
fn export_expression(expression: &Atom) -> String {
    expression
        .printer(PrintOptions::file_no_namespace())
        .to_string()
}

/// The UFO convention stores twice the spin plus one (scalar = 1, spinor = 2,
/// vector = 3, and tensor = 5).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub(crate) struct ParticleDefinition {
    pub pdg_code: i64,
    pub name: String,
    pub antiname: String,
    pub spin: i64,
    pub color: i64,
    pub mass: String,
    pub width: String,
    pub texname: String,
    pub antitexname: String,
    pub charge: f64,
    pub ghost_number: i64,
    pub lepton_number: i64,
    pub y_charge: i64,
    #[serde(default = "default_true")]
    pub propagating: bool,
    #[serde(
        default,
        rename = "goldstoneboson",
        alias = "goldstone",
        alias = "GoldstoneBoson"
    )]
    pub goldstone: bool,
    #[serde(default)]
    pub propagator: Option<String>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub(crate) struct PropagatorDefinition {
    pub name: String,
    pub particle: String,
    pub numerator: String,
    pub denominator: String,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub(crate) struct LorentzStructureDefinition {
    pub name: String,
    pub spins: Vec<i64>,
    pub structure: String,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub(crate) struct CouplingDefinition {
    pub name: String,
    pub expression: String,
    #[serde(with = "vectorize")]
    pub orders: BTreeMap<String, usize>,
    #[serde(default, with = "crate::card::optional_complex")]
    pub value: Option<ComplexValue>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum ParameterNature {
    #[default]
    External,
    Internal,
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum ParameterType {
    #[default]
    Real,
    Complex,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub(crate) struct ParameterDefinition {
    pub name: String,
    pub lhablock: Option<String>,
    pub lhacode: Option<Vec<usize>>,
    pub nature: ParameterNature,
    pub parameter_type: ParameterType,
    #[serde(default, with = "crate::card::optional_complex")]
    pub value: Option<ComplexValue>,
    pub expression: Option<String>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Order {
    pub name: String,
    pub expansion_order: i64,
    pub hierarchy: i64,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub(crate) struct VertexRuleDefinition {
    pub name: String,
    pub particles: Vec<String>,
    pub color_structures: Vec<String>,
    pub lorentz_structures: Vec<String>,
    pub couplings: Vec<Vec<Option<String>>>,
}

/// A callable function declared by a UFO model.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub(crate) struct ModelFunctionDefinition {
    pub name: String,
    #[serde(default)]
    pub arguments: Vec<String>,
    pub expression: Option<String>,
}

/// Metadata for a UFO form factor.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub(crate) struct ModelFormFactorDefinition {
    pub name: String,
    #[serde(rename = "type")]
    pub type_name: Option<String>,
    pub value: Option<String>,
}

/// Serializable data from which a validated [`Model`] is constructed.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub(crate) struct ModelDefinition {
    pub name: String,
    pub restriction: Option<String>,
    pub orders: Vec<Order>,
    pub parameters: Vec<ParameterDefinition>,
    pub particles: Vec<ParticleDefinition>,
    pub propagators: Vec<PropagatorDefinition>,
    pub lorentz_structures: Vec<LorentzStructureDefinition>,
    pub couplings: Vec<CouplingDefinition>,
    pub vertex_rules: Vec<VertexRuleDefinition>,
    #[serde(default)]
    pub functions: Vec<ModelFunctionDefinition>,
    #[serde(default)]
    pub form_factors: Vec<ModelFormFactorDefinition>,
}

macro_rules! entity_id {
    ($name:ident) => {
        #[derive(
            Clone,
            Copy,
            Debug,
            PartialEq,
            Eq,
            PartialOrd,
            Ord,
            Hash,
            Serialize,
            Deserialize,
            bincode_trait_derive::Encode,
            bincode_trait_derive::Decode,
        )]
        pub struct $name(usize);

        impl $name {
            pub const fn from_index(index: usize) -> Self {
                Self(index)
            }

            pub const fn index(self) -> usize {
                self.0
            }
        }
    };
}

entity_id!(OrderId);
entity_id!(ParameterId);
entity_id!(ParticleId);
entity_id!(PropagatorId);
entity_id!(LorentzStructureId);
entity_id!(CouplingId);
entity_id!(VertexRuleId);
entity_id!(ModelFunctionId);
entity_id!(ModelFormFactorId);

/// Content identity of a validated model.
#[derive(
    Clone,
    Copy,
    Debug,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    Serialize,
    Deserialize,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
)]
pub struct ModelFingerprint([u8; 32]);

impl ModelFingerprint {
    pub const fn as_bytes(&self) -> &[u8; 32] {
        &self.0
    }
}

impl std::fmt::Display for ModelFingerprint {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for byte in self.0 {
            write!(formatter, "{byte:02x}")?;
        }
        Ok(())
    }
}

/// A particle whose links into the model are stable typed identifiers.
#[derive(Clone, Debug, PartialEq)]
pub struct Particle {
    pub pdg_code: i64,
    pub name: String,
    pub antiparticle: ParticleId,
    pub spin: i64,
    pub color: i64,
    pub mass: ParameterId,
    pub width: ParameterId,
    pub texname: String,
    pub antitexname: String,
    pub charge: f64,
    pub ghost_number: i64,
    pub lepton_number: i64,
    pub y_charge: i64,
    pub propagating: bool,
    pub goldstone: bool,
    pub propagator: Option<PropagatorId>,
}

impl Particle {
    pub fn is_antiparticle(&self) -> bool {
        self.pdg_code < 0
    }

    pub fn is_fermion(&self) -> bool {
        self.spin % 2 == 0
    }

    pub fn is_vector(&self) -> bool {
        self.spin == 3
    }

    pub fn is_tensor(&self) -> bool {
        self.spin == 5
    }

    pub fn is_scalar(&self) -> bool {
        self.spin == 1
    }

    pub fn is_spinor(&self) -> bool {
        self.spin == 2
    }

    pub fn is_ghost(&self) -> bool {
        self.ghost_number != 0
    }

    pub fn is_goldstone(&self) -> bool {
        self.goldstone
    }

    pub fn is_propagating(&self) -> bool {
        self.propagating
    }

    pub fn is_qcd_charged(&self) -> bool {
        self.color != 1
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Propagator {
    pub name: String,
    pub particle: ParticleId,
    pub numerator: Atom,
    pub denominator: Atom,
}

#[derive(Clone, Debug, PartialEq)]
pub struct LorentzStructure {
    pub name: String,
    pub spins: Vec<i64>,
    pub structure: Atom,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Coupling {
    pub name: String,
    pub expression: Atom,
    pub orders: BTreeMap<String, usize>,
    pub value: Option<ComplexValue>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Parameter {
    pub name: String,
    pub lhablock: Option<String>,
    pub lhacode: Option<Vec<usize>>,
    pub nature: ParameterNature,
    pub parameter_type: ParameterType,
    pub value: Option<ComplexValue>,
    pub expression: Option<Atom>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct VertexRule {
    pub name: String,
    pub particles: Vec<ParticleId>,
    pub color_structures: Vec<Atom>,
    pub lorentz_structures: Vec<LorentzStructureId>,
    pub couplings: Vec<Vec<Option<CouplingId>>>,
}

impl VertexRule {
    /// Return the largest power of every coupling order represented by an
    /// entry in this vertex's color/Lorentz coupling matrix.
    pub fn coupling_orders(&self, model: &Model) -> BTreeMap<String, usize> {
        let mut orders: BTreeMap<String, usize> = BTreeMap::new();
        for coupling in self.couplings.iter().flatten().flatten() {
            for (name, power) in &model.coupling_by_id(*coupling).unwrap().orders {
                orders
                    .entry(name.clone())
                    .and_modify(|current| *current = (*current).max(*power))
                    .or_insert(*power);
            }
        }
        orders
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct ModelFunction {
    pub name: String,
    pub arguments: Vec<String>,
    pub expression: Option<Atom>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ModelFormFactor {
    pub name: String,
    pub type_name: Option<String>,
    pub value: Option<Atom>,
}

#[derive(Clone, Debug, Default)]
struct ModelIndexes {
    orders: HashMap<String, OrderId>,
    parameters: HashMap<String, ParameterId>,
    particles: HashMap<String, ParticleId>,
    particles_by_pdg: HashMap<i64, ParticleId>,
    propagators: HashMap<String, PropagatorId>,
    lorentz_structures: HashMap<String, LorentzStructureId>,
    couplings: HashMap<String, CouplingId>,
    vertex_rules: HashMap<String, VertexRuleId>,
    functions: HashMap<String, ModelFunctionId>,
    form_factors: HashMap<String, ModelFormFactorId>,
}

/// A validated physics model with fallible, indexed lookups.
#[derive(Clone, Debug)]
pub struct Model {
    name: String,
    restriction: Option<String>,
    orders: Vec<Order>,
    parameters: Vec<Parameter>,
    particles: Vec<Particle>,
    propagators: Vec<Propagator>,
    lorentz_structures: Vec<LorentzStructure>,
    couplings: Vec<Coupling>,
    vertex_rules: Vec<VertexRule>,
    functions: Vec<ModelFunction>,
    form_factors: Vec<ModelFormFactor>,
    indexes: ModelIndexes,
    coupling_replacement_rules: OnceLock<Vec<Replacement>>,
}

impl Model {
    /// Construct an empty, validated model for application state that has not
    /// loaded a physics model yet.
    pub fn empty(name: impl Into<String>) -> Self {
        Self::new(ModelDefinition {
            name: name.into(),
            restriction: None,
            orders: Vec::new(),
            parameters: Vec::new(),
            particles: Vec::new(),
            propagators: Vec::new(),
            lorentz_structures: Vec::new(),
            couplings: Vec::new(),
            vertex_rules: Vec::new(),
            functions: Vec::new(),
            form_factors: Vec::new(),
        })
        .expect("an empty model definition is valid")
    }

    fn new(definition: ModelDefinition) -> Result<Self, ModelError> {
        let indexes = ModelIndexes::build(&definition)?;
        let parameters = definition
            .parameters
            .iter()
            .map(|parameter| {
                Ok(Parameter {
                    name: parameter.name.clone(),
                    lhablock: parameter.lhablock.clone(),
                    lhacode: parameter.lhacode.clone(),
                    nature: parameter.nature.clone(),
                    parameter_type: parameter.parameter_type.clone(),
                    value: parameter.value,
                    expression: parameter
                        .expression
                        .as_deref()
                        .map(|expression| {
                            parse_expression(
                                EntityKind::Parameter,
                                &parameter.name,
                                "expression",
                                expression,
                            )
                        })
                        .transpose()?,
                })
            })
            .collect::<Result<Vec<_>, ModelError>>()?;
        let particles = definition
            .particles
            .iter()
            .map(|particle| Particle {
                pdg_code: particle.pdg_code,
                name: particle.name.clone(),
                antiparticle: indexes.particles[&particle.antiname],
                spin: particle.spin,
                color: particle.color,
                mass: indexes.parameters[&particle.mass],
                width: indexes.parameters[&particle.width],
                texname: particle.texname.clone(),
                antitexname: particle.antitexname.clone(),
                charge: particle.charge,
                ghost_number: particle.ghost_number,
                lepton_number: particle.lepton_number,
                y_charge: particle.y_charge,
                propagating: particle.propagating,
                goldstone: particle.goldstone,
                propagator: particle
                    .propagator
                    .as_ref()
                    .map(|name| indexes.propagators[name]),
            })
            .collect();
        let propagators = definition
            .propagators
            .iter()
            .map(|propagator| {
                Ok(Propagator {
                    name: propagator.name.clone(),
                    particle: indexes.particles[&propagator.particle],
                    numerator: parse_expression(
                        EntityKind::Propagator,
                        &propagator.name,
                        "numerator",
                        &propagator.numerator,
                    )?,
                    denominator: parse_expression(
                        EntityKind::Propagator,
                        &propagator.name,
                        "denominator",
                        &propagator.denominator,
                    )?,
                })
            })
            .collect::<Result<Vec<_>, ModelError>>()?;
        let lorentz_structures = definition
            .lorentz_structures
            .iter()
            .map(|lorentz| {
                Ok(LorentzStructure {
                    name: lorentz.name.clone(),
                    spins: lorentz.spins.clone(),
                    structure: parse_expression(
                        EntityKind::LorentzStructure,
                        &lorentz.name,
                        "structure",
                        &lorentz.structure,
                    )?,
                })
            })
            .collect::<Result<Vec<_>, ModelError>>()?;
        let couplings = definition
            .couplings
            .iter()
            .map(|coupling| {
                Ok(Coupling {
                    name: coupling.name.clone(),
                    expression: parse_expression(
                        EntityKind::Coupling,
                        &coupling.name,
                        "expression",
                        &coupling.expression,
                    )?,
                    orders: coupling.orders.clone(),
                    value: coupling.value,
                })
            })
            .collect::<Result<Vec<_>, ModelError>>()?;
        let vertex_rules = definition
            .vertex_rules
            .iter()
            .map(|vertex| {
                Ok(VertexRule {
                    name: vertex.name.clone(),
                    particles: vertex
                        .particles
                        .iter()
                        .map(|name| indexes.particles[name])
                        .collect(),
                    color_structures: vertex
                        .color_structures
                        .iter()
                        .map(|expression| {
                            parse_expression(
                                EntityKind::VertexRule,
                                &vertex.name,
                                "color_structures",
                                expression,
                            )
                        })
                        .collect::<Result<Vec<_>, _>>()?,
                    lorentz_structures: vertex
                        .lorentz_structures
                        .iter()
                        .map(|name| indexes.lorentz_structures[name])
                        .collect(),
                    couplings: vertex
                        .couplings
                        .iter()
                        .map(|row| {
                            row.iter()
                                .map(|name| name.as_ref().map(|name| indexes.couplings[name]))
                                .collect()
                        })
                        .collect(),
                })
            })
            .collect::<Result<Vec<_>, ModelError>>()?;
        let functions = definition
            .functions
            .iter()
            .map(|function| {
                Ok(ModelFunction {
                    name: function.name.clone(),
                    arguments: function.arguments.clone(),
                    expression: function
                        .expression
                        .as_deref()
                        .map(|expression| {
                            parse_expression(
                                EntityKind::ModelFunction,
                                &function.name,
                                "expression",
                                expression,
                            )
                        })
                        .transpose()?,
                })
            })
            .collect::<Result<Vec<_>, ModelError>>()?;
        let form_factors = definition
            .form_factors
            .iter()
            .map(|form_factor| {
                Ok(ModelFormFactor {
                    name: form_factor.name.clone(),
                    type_name: form_factor.type_name.clone(),
                    value: form_factor
                        .value
                        .as_deref()
                        .map(|value| {
                            parse_expression(
                                EntityKind::FormFactor,
                                &form_factor.name,
                                "value",
                                value,
                            )
                        })
                        .transpose()?,
                })
            })
            .collect::<Result<Vec<_>, ModelError>>()?;
        Ok(Self {
            name: definition.name,
            restriction: definition.restriction,
            orders: definition.orders,
            parameters,
            particles,
            propagators,
            lorentz_structures,
            couplings,
            vertex_rules,
            functions,
            form_factors,
            indexes,
            coupling_replacement_rules: OnceLock::new(),
        })
    }

    pub fn from_json(json: &str) -> Result<Self, ModelError> {
        Self::new(serde_json::from_str(json)?)
    }

    pub fn from_path(path: impl AsRef<Path>) -> Result<Self, ModelError> {
        let path = path.as_ref();
        let json = fs::read_to_string(path).map_err(|source| ModelError::Read {
            path: path.to_path_buf(),
            source,
        })?;
        Self::from_json(&json)
    }

    pub fn to_json(&self) -> Result<String, ModelError> {
        Ok(serde_json::to_string(self)?)
    }

    pub fn to_json_pretty(&self) -> Result<String, ModelError> {
        Ok(serde_json::to_string_pretty(self)?)
    }

    pub fn write_json(&self, path: impl AsRef<Path>) -> Result<(), ModelError> {
        let path = path.as_ref();
        fs::write(path, self.to_json_pretty()?).map_err(|source| ModelError::Write {
            path: path.to_path_buf(),
            source,
        })
    }

    fn definition(&self) -> ModelDefinition {
        ModelDefinition {
            name: self.name.clone(),
            restriction: self.restriction.clone(),
            orders: self.orders.clone(),
            parameters: self
                .parameters
                .iter()
                .map(|parameter| ParameterDefinition {
                    name: parameter.name.clone(),
                    lhablock: parameter.lhablock.clone(),
                    lhacode: parameter.lhacode.clone(),
                    nature: parameter.nature.clone(),
                    parameter_type: parameter.parameter_type.clone(),
                    value: parameter.value,
                    expression: parameter.expression.as_ref().map(export_expression),
                })
                .collect(),
            particles: self
                .particles
                .iter()
                .map(|particle| ParticleDefinition {
                    pdg_code: particle.pdg_code,
                    name: particle.name.clone(),
                    antiname: self.particles[particle.antiparticle.index()].name.clone(),
                    spin: particle.spin,
                    color: particle.color,
                    mass: self.parameters[particle.mass.index()].name.clone(),
                    width: self.parameters[particle.width.index()].name.clone(),
                    texname: particle.texname.clone(),
                    antitexname: particle.antitexname.clone(),
                    charge: particle.charge,
                    ghost_number: particle.ghost_number,
                    lepton_number: particle.lepton_number,
                    y_charge: particle.y_charge,
                    propagating: particle.propagating,
                    goldstone: particle.goldstone,
                    propagator: particle
                        .propagator
                        .map(|id| self.propagators[id.index()].name.clone()),
                })
                .collect(),
            propagators: self
                .propagators
                .iter()
                .map(|propagator| PropagatorDefinition {
                    name: propagator.name.clone(),
                    particle: self.particles[propagator.particle.index()].name.clone(),
                    numerator: export_expression(&propagator.numerator),
                    denominator: export_expression(&propagator.denominator),
                })
                .collect(),
            lorentz_structures: self
                .lorentz_structures
                .iter()
                .map(|lorentz| LorentzStructureDefinition {
                    name: lorentz.name.clone(),
                    spins: lorentz.spins.clone(),
                    structure: export_expression(&lorentz.structure),
                })
                .collect(),
            couplings: self
                .couplings
                .iter()
                .map(|coupling| CouplingDefinition {
                    name: coupling.name.clone(),
                    expression: export_expression(&coupling.expression),
                    orders: coupling.orders.clone(),
                    value: coupling.value,
                })
                .collect(),
            vertex_rules: self
                .vertex_rules
                .iter()
                .map(|vertex| VertexRuleDefinition {
                    name: vertex.name.clone(),
                    particles: vertex
                        .particles
                        .iter()
                        .map(|id| self.particles[id.index()].name.clone())
                        .collect(),
                    color_structures: vertex
                        .color_structures
                        .iter()
                        .map(export_expression)
                        .collect(),
                    lorentz_structures: vertex
                        .lorentz_structures
                        .iter()
                        .map(|id| self.lorentz_structures[id.index()].name.clone())
                        .collect(),
                    couplings: vertex
                        .couplings
                        .iter()
                        .map(|row| {
                            row.iter()
                                .map(|id| id.map(|id| self.couplings[id.index()].name.clone()))
                                .collect()
                        })
                        .collect(),
                })
                .collect(),
            functions: self
                .functions
                .iter()
                .map(|function| ModelFunctionDefinition {
                    name: function.name.clone(),
                    arguments: function.arguments.clone(),
                    expression: function.expression.as_ref().map(export_expression),
                })
                .collect(),
            form_factors: self
                .form_factors
                .iter()
                .map(|form_factor| ModelFormFactorDefinition {
                    name: form_factor.name.clone(),
                    type_name: form_factor.type_name.clone(),
                    value: form_factor.value.as_ref().map(export_expression),
                })
                .collect(),
        }
    }

    #[cfg(test)]
    pub(crate) fn into_definition(self) -> ModelDefinition {
        self.definition()
    }

    pub fn fingerprint(&self) -> ModelFingerprint {
        let bytes = serde_json::to_vec(&self.definition())
            .expect("canonical model definitions are always JSON serializable");
        ModelFingerprint(*blake3::hash(&bytes).as_bytes())
    }

    /// Remove parameters and interactions that evaluate to zero, rebuilding
    /// all typed identifiers atomically after compaction.
    pub fn simplify_zero_values(&mut self) -> Result<(), ModelError> {
        let mut definition = self.definition();
        for parameter in &mut definition.parameters {
            if parameter.nature == ParameterNature::External
                && parameter
                    .value
                    .is_some_and(|value| value.re == 0.0 && value.im == 0.0)
            {
                parameter.nature = ParameterNature::Internal;
                parameter.expression = Some("0".to_owned());
            }
        }
        let zero_couplings = self
            .couplings()
            .iter()
            .filter(|coupling| {
                coupling
                    .value
                    .is_some_and(|value| value.re == 0.0 && value.im == 0.0)
            })
            .map(|coupling| coupling.name.clone())
            .collect::<BTreeSet<_>>();
        for vertex in &mut definition.vertex_rules {
            for coupling in vertex.couplings.iter_mut().flatten() {
                if coupling
                    .as_ref()
                    .is_some_and(|name| zero_couplings.contains(name))
                {
                    *coupling = None;
                }
            }
        }
        definition
            .vertex_rules
            .retain(|vertex| vertex.couplings.iter().flatten().any(Option::is_some));
        definition
            .couplings
            .retain(|coupling| !zero_couplings.contains(&coupling.name));
        *self = Self::new(definition)?;
        Ok(())
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn restriction(&self) -> Option<&str> {
        self.restriction.as_deref()
    }

    /// Return an independently validated model with updated restriction metadata.
    ///
    /// The original model is left unchanged, and the returned model has its own
    /// canonical fingerprint. This keeps restriction changes from bypassing the
    /// model's validation and indexing boundary.
    pub fn with_restriction(&self, restriction: Option<String>) -> Result<Self, ModelError> {
        let mut definition = self.definition();
        definition.restriction = restriction;
        Self::new(definition)
    }

    pub fn orders(&self) -> &[Order] {
        &self.orders
    }

    pub fn parameters(&self) -> &[Parameter] {
        &self.parameters
    }

    pub fn particles(&self) -> &[Particle] {
        &self.particles
    }

    pub fn propagators(&self) -> &[Propagator] {
        &self.propagators
    }

    pub fn lorentz_structures(&self) -> &[LorentzStructure] {
        &self.lorentz_structures
    }

    pub fn couplings(&self) -> &[Coupling] {
        &self.couplings
    }

    pub fn vertex_rules(&self) -> &[VertexRule] {
        &self.vertex_rules
    }

    pub fn functions(&self) -> &[ModelFunction] {
        &self.functions
    }

    pub fn form_factors(&self) -> &[ModelFormFactor] {
        &self.form_factors
    }

    pub fn particle_id(&self, name: &str) -> Result<ParticleId, ModelError> {
        self.indexes
            .particles
            .get(name)
            .copied()
            .ok_or_else(|| self.not_found(EntityKind::Particle, name))
    }

    pub fn particle_id_by_pdg(&self, pdg: i64) -> Result<ParticleId, ModelError> {
        self.indexes
            .particles_by_pdg
            .get(&pdg)
            .copied()
            .ok_or_else(|| self.not_found(EntityKind::Particle, &pdg.to_string()))
    }

    pub fn particle(&self, name: &str) -> Result<&Particle, ModelError> {
        self.particle_by_id(self.particle_id(name)?)
    }

    pub fn particle_by_pdg(&self, pdg: i64) -> Result<&Particle, ModelError> {
        self.particle_by_id(self.particle_id_by_pdg(pdg)?)
    }

    pub fn particle_by_id(&self, id: ParticleId) -> Result<&Particle, ModelError> {
        self.particles
            .get(id.0)
            .ok_or_else(|| self.not_found(EntityKind::Particle, &id.0.to_string()))
    }

    pub fn particle_id_at(&self, index: usize) -> Result<ParticleId, ModelError> {
        self.particle_by_id(ParticleId(index))?;
        Ok(ParticleId(index))
    }

    pub fn antiparticle(&self, particle: &Particle) -> Result<&Particle, ModelError> {
        self.particle_by_id(particle.antiparticle)
    }

    pub fn particle_is_self_conjugate(&self, particle: ParticleId) -> bool {
        self.particle_by_id(particle)
            .is_ok_and(|record| record.antiparticle == particle)
    }

    pub fn particle_mass(&self, particle: ParticleId) -> Result<&Parameter, ModelError> {
        let particle = self.particle_by_id(particle)?;
        self.parameter_by_id(particle.mass)
    }

    pub fn particle_width(&self, particle: ParticleId) -> Result<&Parameter, ModelError> {
        let particle = self.particle_by_id(particle)?;
        self.parameter_by_id(particle.width)
    }

    pub fn particle_is_massless(&self, particle: ParticleId) -> bool {
        self.particle_mass(particle).is_ok_and(|mass| {
            mass.name == "ZERO"
                || mass
                    .value
                    .is_some_and(|value| value.re == 0.0 && value.im == 0.0)
                || mass
                    .expression
                    .as_ref()
                    .is_some_and(|expression| expression.is_zero())
        })
    }

    pub fn parameter_id(&self, name: &str) -> Result<ParameterId, ModelError> {
        self.lookup_id(EntityKind::Parameter, name, &self.indexes.parameters)
    }

    pub fn parameter(&self, name: &str) -> Result<&Parameter, ModelError> {
        self.parameter_by_id(self.parameter_id(name)?)
    }

    pub fn parameter_by_id(&self, id: ParameterId) -> Result<&Parameter, ModelError> {
        self.parameters
            .get(id.index())
            .ok_or_else(|| self.not_found(EntityKind::Parameter, &id.index().to_string()))
    }

    pub fn coupling_id(&self, name: &str) -> Result<CouplingId, ModelError> {
        self.lookup_id(EntityKind::Coupling, name, &self.indexes.couplings)
    }

    pub fn coupling(&self, name: &str) -> Result<&Coupling, ModelError> {
        self.coupling_by_id(self.coupling_id(name)?)
    }

    pub fn coupling_by_id(&self, id: CouplingId) -> Result<&Coupling, ModelError> {
        self.couplings
            .get(id.index())
            .ok_or_else(|| self.not_found(EntityKind::Coupling, &id.index().to_string()))
    }

    /// Replace named UFO vertex coefficients by their analytic model expressions.
    ///
    /// Generated diagrams deliberately retain symbols such as `UFO::GC_11` so
    /// numerical consumers can use the model's precomputed coupling values.
    /// This non-mutating expansion is the explicit boundary for symbolic
    /// inspection and algebra in terms of the underlying model parameters.
    pub fn expand_couplings(&self, expression: &Atom) -> Atom {
        expression.replace_multiple(self.coupling_replacement_rules())
    }

    fn coupling_replacement_rules(&self) -> &[Replacement] {
        self.coupling_replacement_rules.get_or_init(|| {
            self.couplings
                .iter()
                .map(|coupling| {
                    Replacement::new(
                        Atom::var(symbol!(&format!("UFO::{}", coupling.name))).to_pattern(),
                        coupling.expression.clone(),
                    )
                })
                .collect()
        })
    }

    pub fn order_id(&self, name: &str) -> Result<OrderId, ModelError> {
        self.lookup_id(EntityKind::Order, name, &self.indexes.orders)
    }

    pub fn order(&self, name: &str) -> Result<&Order, ModelError> {
        self.order_by_id(self.order_id(name)?)
    }

    pub fn order_by_id(&self, id: OrderId) -> Result<&Order, ModelError> {
        self.orders
            .get(id.index())
            .ok_or_else(|| self.not_found(EntityKind::Order, &id.index().to_string()))
    }

    pub fn propagator_id(&self, name: &str) -> Result<PropagatorId, ModelError> {
        self.lookup_id(EntityKind::Propagator, name, &self.indexes.propagators)
    }

    pub fn propagator(&self, name: &str) -> Result<&Propagator, ModelError> {
        self.propagator_by_id(self.propagator_id(name)?)
    }

    pub fn propagator_by_id(&self, id: PropagatorId) -> Result<&Propagator, ModelError> {
        self.propagators
            .get(id.index())
            .ok_or_else(|| self.not_found(EntityKind::Propagator, &id.index().to_string()))
    }

    pub fn propagator_id_at(&self, index: usize) -> Result<PropagatorId, ModelError> {
        self.propagator_by_id(PropagatorId(index))?;
        Ok(PropagatorId(index))
    }

    pub fn lorentz_structure_id(&self, name: &str) -> Result<LorentzStructureId, ModelError> {
        self.lookup_id(
            EntityKind::LorentzStructure,
            name,
            &self.indexes.lorentz_structures,
        )
    }

    pub fn lorentz_structure(&self, name: &str) -> Result<&LorentzStructure, ModelError> {
        self.lorentz_structure_by_id(self.lorentz_structure_id(name)?)
    }

    pub fn lorentz_structure_by_id(
        &self,
        id: LorentzStructureId,
    ) -> Result<&LorentzStructure, ModelError> {
        self.lorentz_structures
            .get(id.index())
            .ok_or_else(|| self.not_found(EntityKind::LorentzStructure, &id.index().to_string()))
    }

    pub fn vertex_rule_id(&self, name: &str) -> Result<VertexRuleId, ModelError> {
        self.indexes
            .vertex_rules
            .get(name)
            .copied()
            .ok_or_else(|| self.not_found(EntityKind::VertexRule, name))
    }

    pub fn vertex_rule(&self, name: &str) -> Result<&VertexRule, ModelError> {
        self.vertex_rule_by_id(self.vertex_rule_id(name)?)
    }

    pub fn vertex_rule_by_id(&self, id: VertexRuleId) -> Result<&VertexRule, ModelError> {
        self.vertex_rules
            .get(id.0)
            .ok_or_else(|| self.not_found(EntityKind::VertexRule, &id.0.to_string()))
    }

    pub fn vertex_rule_id_at(&self, index: usize) -> Result<VertexRuleId, ModelError> {
        self.vertex_rule_by_id(VertexRuleId(index))?;
        Ok(VertexRuleId(index))
    }

    pub fn function(&self, name: &str) -> Result<&ModelFunction, ModelError> {
        self.function_by_id(self.function_id(name)?)
    }

    pub fn function_id(&self, name: &str) -> Result<ModelFunctionId, ModelError> {
        self.lookup_id(EntityKind::ModelFunction, name, &self.indexes.functions)
    }

    pub fn function_by_id(&self, id: ModelFunctionId) -> Result<&ModelFunction, ModelError> {
        self.functions
            .get(id.index())
            .ok_or_else(|| self.not_found(EntityKind::ModelFunction, &id.index().to_string()))
    }

    pub fn form_factor(&self, name: &str) -> Result<&ModelFormFactor, ModelError> {
        self.form_factor_by_id(self.form_factor_id(name)?)
    }

    pub fn form_factor_id(&self, name: &str) -> Result<ModelFormFactorId, ModelError> {
        self.lookup_id(EntityKind::FormFactor, name, &self.indexes.form_factors)
    }

    pub fn form_factor_by_id(&self, id: ModelFormFactorId) -> Result<&ModelFormFactor, ModelError> {
        self.form_factors
            .get(id.index())
            .ok_or_else(|| self.not_found(EntityKind::FormFactor, &id.index().to_string()))
    }

    pub fn default_parameter_card(&self) -> Result<ParameterCard, ModelError> {
        self.parameters
            .iter()
            .filter(|parameter| parameter.nature == ParameterNature::External)
            .map(|parameter| {
                parameter
                    .value
                    .map(|value| (parameter.name.clone(), value))
                    .ok_or_else(|| ModelError::MissingCardValue {
                        name: parameter.name.clone(),
                    })
            })
            .collect()
    }

    /// Apply a parameter card atomically and invalidate dependent values not
    /// explicitly supplied by the card.
    ///
    /// Use [`Self::apply_parameter_card_with`] to recompute those values in the
    /// same atomic operation.
    pub fn apply_parameter_card(&mut self, card: &ParameterCard) -> Result<(), ModelError> {
        let updates =
            card.iter()
                .map(|(name, value)| {
                    let index =
                        self.indexes.parameters.get(name).copied().ok_or_else(|| {
                            ModelError::UnknownCardParameter { name: name.clone() }
                        })?;
                    Ok((index, *value))
                })
                .collect::<Result<Vec<_>, ModelError>>()?;

        if updates.is_empty() {
            return Ok(());
        }
        for parameter in &mut self.parameters {
            if parameter.nature == ParameterNature::Internal && parameter.expression.is_some() {
                parameter.value = None;
            }
        }
        for coupling in &mut self.couplings {
            coupling.value = None;
        }
        for (index, value) in updates {
            self.parameters[index.index()].value = Some(value);
        }
        Ok(())
    }

    /// Recompute all expression-backed parameters and couplings atomically.
    pub fn recompute_with<E>(&mut self, evaluator: &mut E) -> Result<(), RecomputeError<E::Error>>
    where
        E: ModelEvaluator,
    {
        self.recompute_with_overrides(evaluator, &BTreeSet::new())
    }

    fn recompute_with_overrides<E>(
        &mut self,
        evaluator: &mut E,
        overrides: &BTreeSet<String>,
    ) -> Result<(), RecomputeError<E::Error>>
    where
        E: ModelEvaluator,
    {
        let request = self.evaluation_request(overrides);
        let evaluated = evaluator
            .evaluate(request)
            .map_err(RecomputeError::Evaluator)?;
        self.validate_evaluated_values(&evaluated, overrides)?;

        for parameter in &mut self.parameters {
            if parameter.nature == ParameterNature::Internal
                && parameter.expression.is_some()
                && !overrides.contains(&parameter.name)
            {
                parameter.value = evaluated.internal_parameters.get(&parameter.name).copied();
            }
        }
        for coupling in &mut self.couplings {
            coupling.value = evaluated.couplings.get(&coupling.name).copied();
        }
        Ok(())
    }

    /// Apply parameter values and recompute all dependents as one atomic model
    /// update.
    pub fn apply_parameter_card_with<E>(
        &mut self,
        card: &ParameterCard,
        evaluator: &mut E,
    ) -> Result<(), RecomputeError<E::Error>>
    where
        E: ModelEvaluator,
    {
        let mut updated = self.clone();
        let overrides = card
            .keys()
            .filter(|name| {
                updated.indexes.parameters.get(*name).is_some_and(|index| {
                    let parameter = &updated.parameters[index.index()];
                    parameter.nature == ParameterNature::Internal && parameter.expression.is_some()
                })
            })
            .cloned()
            .collect();
        updated.apply_parameter_card(card)?;
        updated.recompute_with_overrides(evaluator, &overrides)?;
        *self = updated;
        Ok(())
    }

    fn evaluation_request(&self, overrides: &BTreeSet<String>) -> EvaluationRequest {
        let known_parameters = self
            .parameters
            .iter()
            .filter(|parameter| {
                parameter.nature == ParameterNature::External
                    || parameter.expression.is_none()
                    || overrides.contains(&parameter.name)
            })
            .filter_map(|parameter| parameter.value.map(|value| (parameter.name.clone(), value)))
            .collect();
        let internal_parameters = self
            .parameters
            .iter()
            .filter(|parameter| {
                parameter.nature == ParameterNature::Internal
                    && !overrides.contains(&parameter.name)
            })
            .filter_map(|parameter| {
                parameter
                    .expression
                    .as_ref()
                    .map(|expression| ModelExpression {
                        name: parameter.name.clone(),
                        expression: export_expression(expression),
                    })
            })
            .collect();
        let couplings = self
            .couplings
            .iter()
            .map(|coupling| ModelExpression {
                name: coupling.name.clone(),
                expression: export_expression(&coupling.expression),
            })
            .collect();
        EvaluationRequest {
            known_parameters,
            internal_parameters,
            couplings,
            functions: self
                .functions
                .iter()
                .map(|function| crate::EvaluationFunction {
                    name: function.name.clone(),
                    arguments: function.arguments.clone(),
                    expression: function.expression.as_ref().map(export_expression),
                })
                .collect(),
            form_factors: self
                .form_factors
                .iter()
                .map(|form_factor| crate::EvaluationFormFactor {
                    name: form_factor.name.clone(),
                    type_name: form_factor.type_name.clone(),
                    value: form_factor.value.as_ref().map(export_expression),
                })
                .collect(),
        }
    }

    fn validate_evaluated_values<E>(
        &self,
        evaluated: &EvaluatedValues,
        overrides: &BTreeSet<String>,
    ) -> Result<(), RecomputeError<E>>
    where
        E: std::error::Error + 'static,
    {
        for parameter in self.parameters.iter().filter(|parameter| {
            parameter.nature == ParameterNature::Internal
                && parameter.expression.is_some()
                && !overrides.contains(&parameter.name)
        }) {
            if !evaluated.internal_parameters.contains_key(&parameter.name) {
                return Err(RecomputeError::MissingValue {
                    kind: EntityKind::Parameter,
                    name: parameter.name.clone(),
                });
            }
        }
        for name in evaluated.internal_parameters.keys() {
            let Some(index) = self.indexes.parameters.get(name) else {
                return Err(RecomputeError::UnexpectedValue {
                    kind: EntityKind::Parameter,
                    name: name.clone(),
                });
            };
            let parameter = &self.parameters[index.index()];
            if parameter.nature != ParameterNature::Internal
                || parameter.expression.is_none()
                || overrides.contains(name)
            {
                return Err(RecomputeError::UnexpectedValue {
                    kind: EntityKind::Parameter,
                    name: name.clone(),
                });
            }
        }
        for coupling in &self.couplings {
            if !evaluated.couplings.contains_key(&coupling.name) {
                return Err(RecomputeError::MissingValue {
                    kind: EntityKind::Coupling,
                    name: coupling.name.clone(),
                });
            }
        }
        for name in evaluated.couplings.keys() {
            if !self.indexes.couplings.contains_key(name) {
                return Err(RecomputeError::UnexpectedValue {
                    kind: EntityKind::Coupling,
                    name: name.clone(),
                });
            }
        }
        Ok(())
    }

    fn lookup_id<I: Copy>(
        &self,
        kind: EntityKind,
        name: &str,
        index: &HashMap<String, I>,
    ) -> Result<I, ModelError> {
        index
            .get(name)
            .copied()
            .ok_or_else(|| self.not_found(kind, name))
    }

    fn not_found(&self, kind: EntityKind, key: &str) -> ModelError {
        ModelError::NotFound {
            model: self.name.clone(),
            kind,
            key: key.to_owned(),
        }
    }
}

impl Serialize for Model {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        self.definition().serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for Model {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let definition = ModelDefinition::deserialize(deserializer)?;
        Self::new(definition).map_err(de::Error::custom)
    }
}

impl ModelIndexes {
    fn build(definition: &ModelDefinition) -> Result<Self, ModelValidationError> {
        if definition.name.trim().is_empty() {
            return Err(ModelValidationError::EmptyModelName);
        }

        let orders = typed_index(
            EntityKind::Order,
            &definition.orders,
            |value| &value.name,
            OrderId,
        )?;
        let parameters = typed_index(
            EntityKind::Parameter,
            &definition.parameters,
            |value| &value.name,
            ParameterId,
        )?;
        let particles = typed_index(
            EntityKind::Particle,
            &definition.particles,
            |value| &value.name,
            ParticleId,
        )?;
        let propagators = typed_index(
            EntityKind::Propagator,
            &definition.propagators,
            |value| &value.name,
            PropagatorId,
        )?;
        let lorentz_structures = typed_index(
            EntityKind::LorentzStructure,
            &definition.lorentz_structures,
            |value| &value.name,
            LorentzStructureId,
        )?;
        let couplings = typed_index(
            EntityKind::Coupling,
            &definition.couplings,
            |value| &value.name,
            CouplingId,
        )?;
        let vertex_rules = typed_index(
            EntityKind::VertexRule,
            &definition.vertex_rules,
            |value| &value.name,
            VertexRuleId,
        )?;
        let functions = typed_index(
            EntityKind::ModelFunction,
            &definition.functions,
            |value| &value.name,
            ModelFunctionId,
        )?;
        let form_factors = typed_index(
            EntityKind::FormFactor,
            &definition.form_factors,
            |value| &value.name,
            ModelFormFactorId,
        )?;

        let mut particles_by_pdg = HashMap::new();
        for (index, particle) in definition.particles.iter().enumerate() {
            if particles_by_pdg
                .insert(particle.pdg_code, ParticleId(index))
                .is_some()
            {
                return Err(ModelValidationError::DuplicatePdg {
                    pdg: particle.pdg_code,
                });
            }
        }

        validate_parameters(&definition.parameters)?;
        validate_couplings(&definition.couplings, &orders)?;
        validate_particles(
            &definition.particles,
            &parameters,
            &particles,
            &definition.propagators,
            &propagators,
        )?;
        validate_propagators(&definition.propagators, &particles)?;
        validate_vertices(
            &definition.vertex_rules,
            &definition.particles,
            &particles,
            &definition.lorentz_structures,
            &lorentz_structures,
            &couplings,
        )?;

        Ok(Self {
            orders,
            parameters,
            particles,
            particles_by_pdg,
            propagators,
            lorentz_structures,
            couplings,
            vertex_rules,
            functions,
            form_factors,
        })
    }
}

fn typed_index<T, I: Copy>(
    kind: EntityKind,
    values: &[T],
    name: impl Fn(&T) -> &str,
    id: impl Fn(usize) -> I,
) -> Result<HashMap<String, I>, ModelValidationError> {
    let mut result = HashMap::with_capacity(values.len());
    for (index, value) in values.iter().enumerate() {
        let name = name(value);
        if result.insert(name.to_owned(), id(index)).is_some() {
            return Err(ModelValidationError::DuplicateName {
                kind,
                name: name.to_owned(),
            });
        }
    }
    Ok(result)
}

fn validate_parameters(parameters: &[ParameterDefinition]) -> Result<(), ModelValidationError> {
    for parameter in parameters {
        match parameter.nature {
            ParameterNature::External if parameter.value.is_none() => {
                return Err(ModelValidationError::ExternalParameterWithoutValue {
                    name: parameter.name.clone(),
                });
            }
            ParameterNature::Internal
                if parameter.value.is_none() && parameter.expression.is_none() =>
            {
                return Err(ModelValidationError::UnresolvedInternalParameter {
                    name: parameter.name.clone(),
                });
            }
            ParameterNature::External | ParameterNature::Internal => {}
        }
    }
    Ok(())
}

fn validate_couplings(
    couplings: &[CouplingDefinition],
    orders: &HashMap<String, OrderId>,
) -> Result<(), ModelValidationError> {
    for coupling in couplings {
        for order in coupling.orders.keys() {
            reference(
                EntityKind::Coupling,
                &coupling.name,
                "orders",
                EntityKind::Order,
                order,
                orders,
            )?;
        }
    }
    Ok(())
}

fn validate_particles(
    particles: &[ParticleDefinition],
    parameters: &HashMap<String, ParameterId>,
    particle_names: &HashMap<String, ParticleId>,
    propagator_definitions: &[PropagatorDefinition],
    propagators: &HashMap<String, PropagatorId>,
) -> Result<(), ModelValidationError> {
    for particle in particles {
        reference(
            EntityKind::Particle,
            &particle.name,
            "mass",
            EntityKind::Parameter,
            &particle.mass,
            parameters,
        )?;
        reference(
            EntityKind::Particle,
            &particle.name,
            "width",
            EntityKind::Parameter,
            &particle.width,
            parameters,
        )?;
        let antiparticle_id = *reference(
            EntityKind::Particle,
            &particle.name,
            "antiname",
            EntityKind::Particle,
            &particle.antiname,
            particle_names,
        )?;
        let antiparticle = &particles[antiparticle_id.0];
        if antiparticle.antiname != particle.name {
            return Err(ModelValidationError::NonReciprocalAntiparticle {
                particle: particle.name.clone(),
                antiparticle: antiparticle.name.clone(),
                actual_antiparticle: antiparticle.antiname.clone(),
            });
        }
        if let Some(propagator) = &particle.propagator {
            let propagator_index = *reference(
                EntityKind::Particle,
                &particle.name,
                "propagator",
                EntityKind::Propagator,
                propagator,
                propagators,
            )?;
            let propagator_definition = &propagator_definitions[propagator_index.index()];
            if propagator_definition.particle != particle.name {
                return Err(ModelValidationError::PropagatorParticleMismatch {
                    particle: particle.name.clone(),
                    propagator: propagator.clone(),
                    propagator_particle: propagator_definition.particle.clone(),
                });
            }
        }
    }
    Ok(())
}

fn validate_propagators(
    propagators: &[PropagatorDefinition],
    particles: &HashMap<String, ParticleId>,
) -> Result<(), ModelValidationError> {
    for propagator in propagators {
        reference(
            EntityKind::Propagator,
            &propagator.name,
            "particle",
            EntityKind::Particle,
            &propagator.particle,
            particles,
        )?;
    }
    Ok(())
}

fn validate_vertices(
    vertices: &[VertexRuleDefinition],
    particle_definitions: &[ParticleDefinition],
    particles: &HashMap<String, ParticleId>,
    lorentz_definitions: &[LorentzStructureDefinition],
    lorentz_structures: &HashMap<String, LorentzStructureId>,
    couplings: &HashMap<String, CouplingId>,
) -> Result<(), ModelValidationError> {
    for vertex in vertices {
        let mut particle_spins = Vec::with_capacity(vertex.particles.len());
        for particle in &vertex.particles {
            let particle_id = *reference(
                EntityKind::VertexRule,
                &vertex.name,
                "particles",
                EntityKind::Particle,
                particle,
                particles,
            )?;
            particle_spins.push(particle_definitions[particle_id.0].spin);
        }
        for lorentz_structure in &vertex.lorentz_structures {
            let lorentz_index = *reference(
                EntityKind::VertexRule,
                &vertex.name,
                "lorentz_structures",
                EntityKind::LorentzStructure,
                lorentz_structure,
                lorentz_structures,
            )?;
            let lorentz_spins = &lorentz_definitions[lorentz_index.index()].spins;
            if lorentz_spins != &particle_spins {
                return Err(ModelValidationError::LorentzSpinMismatch {
                    vertex: vertex.name.clone(),
                    lorentz_structure: lorentz_structure.clone(),
                    particle_spins,
                    lorentz_spins: lorentz_spins.clone(),
                });
            }
        }
        if vertex.couplings.len() != vertex.color_structures.len() {
            return Err(ModelValidationError::CouplingRowCount {
                name: vertex.name.clone(),
                expected: vertex.color_structures.len(),
                actual: vertex.couplings.len(),
            });
        }
        for (row, entries) in vertex.couplings.iter().enumerate() {
            if entries.len() != vertex.lorentz_structures.len() {
                return Err(ModelValidationError::CouplingColumnCount {
                    name: vertex.name.clone(),
                    row,
                    expected: vertex.lorentz_structures.len(),
                    actual: entries.len(),
                });
            }
            for coupling in entries.iter().flatten() {
                reference(
                    EntityKind::VertexRule,
                    &vertex.name,
                    "couplings",
                    EntityKind::Coupling,
                    coupling,
                    couplings,
                )?;
            }
        }
    }
    Ok(())
}

fn reference<'a, V>(
    kind: EntityKind,
    name: &str,
    field: &'static str,
    target_kind: EntityKind,
    target: &str,
    index: &'a HashMap<String, V>,
) -> Result<&'a V, ModelValidationError> {
    index
        .get(target)
        .ok_or_else(|| ModelValidationError::UnknownReference {
            kind,
            name: name.to_owned(),
            field,
            target_kind,
            target: target.to_owned(),
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn model_json(coupling_name: &str) -> String {
        format!(
            r#"{{
                "name": "scalar",
                "restriction": null,
                "orders": [{{"name":"QED","expansion_order":99,"hierarchy":1}}],
                "parameters": [{{
                    "name":"ZERO","lhablock":null,"lhacode":null,"nature":"internal",
                    "parameter_type":"real","value":[0.0,0.0],"expression":null
                }},{{
                    "name":"mass","lhablock":"MASS","lhacode":[1],"nature":"external",
                    "parameter_type":"real","value":[1.0,0.0],"expression":null
                }},{{
                    "name":"double_mass","lhablock":null,"lhacode":null,"nature":"internal",
                    "parameter_type":"real","value":[2.0,0.0],"expression":"2*mass"
                }}],
                "particles": [{{
                    "pdg_code":1,"name":"s","antiname":"s","spin":1,"color":1,
                    "mass":"mass","width":"ZERO","texname":"s","antitexname":"s",
                    "charge":0.0,"ghost_number":0,"lepton_number":0,"y_charge":0,
                    "propagating":false,"goldstoneboson":true,"propagator":"s_prop"
                }}],
                "propagators": [{{
                    "name":"s_prop","particle":"s","numerator":"1","denominator":"P^2-mass^2"
                }}],
                "lorentz_structures": [{{"name":"L1","spins":[1,1,1],"structure":"1"}}],
                "couplings": [{{"name":"GC1","expression":"double_mass","orders":[["QED",1]],"value":[2.0,0.0]}}],
                "vertex_rules": [{{
                    "name":"V1","particles":["s","s","s"],"color_structures":["1"],
                    "lorentz_structures":["L1"],"couplings":[["{coupling_name}"]]
                }}],
                "functions":[{{"name":"twice","arguments":["x"],"expression":"2*x"}}],
                "form_factors":[{{"name":"FF1","type":"complex","value":"twice(P(1)^2)"}}]
            }}"#
        )
    }

    #[test]
    fn parses_ufo_json_and_builds_typed_indexes() {
        let model = Model::from_json(&model_json("GC1")).unwrap();

        let particle_id = model.particle_id_by_pdg(1).unwrap();
        assert_eq!(model.particle_by_id(particle_id).unwrap().name, "s");
        assert_eq!(
            model
                .antiparticle(model.particle("s").unwrap())
                .unwrap()
                .name,
            "s"
        );
        assert!(model.particle("s").unwrap().is_scalar());
        assert!(!model.particle("s").unwrap().is_propagating());
        assert!(model.particle("s").unwrap().is_goldstone());
        assert_eq!(
            model.particle("s").unwrap().propagator.map(|id| model
                .propagator_by_id(id)
                .unwrap()
                .name
                .as_str()),
            Some("s_prop")
        );
        assert_eq!(model.coupling("GC1").unwrap().orders["QED"], 1);
        assert_eq!(model.function("twice").unwrap().arguments, ["x"]);
        assert_eq!(
            model.form_factor("FF1").unwrap().type_name.as_deref(),
            Some("complex")
        );

        assert_eq!(
            model.order_by_id(model.order_id("QED").unwrap()).unwrap(),
            model.order("QED").unwrap()
        );
        assert_eq!(
            model
                .parameter_by_id(model.parameter_id("mass").unwrap())
                .unwrap(),
            model.parameter("mass").unwrap()
        );
        assert_eq!(
            model
                .propagator_by_id(model.propagator_id("s_prop").unwrap())
                .unwrap(),
            model.propagator("s_prop").unwrap()
        );
        assert_eq!(
            model
                .lorentz_structure_by_id(model.lorentz_structure_id("L1").unwrap())
                .unwrap(),
            model.lorentz_structure("L1").unwrap()
        );
        assert_eq!(
            model
                .coupling_by_id(model.coupling_id("GC1").unwrap())
                .unwrap(),
            model.coupling("GC1").unwrap()
        );
        assert_eq!(
            model
                .vertex_rule_by_id(model.vertex_rule_id("V1").unwrap())
                .unwrap(),
            model.vertex_rule("V1").unwrap()
        );
        assert_eq!(
            model
                .function_by_id(model.function_id("twice").unwrap())
                .unwrap(),
            model.function("twice").unwrap()
        );
        assert_eq!(
            model
                .form_factor_by_id(model.form_factor_id("FF1").unwrap())
                .unwrap(),
            model.form_factor("FF1").unwrap()
        );
    }

    #[test]
    fn expands_named_couplings_without_mutating_the_source_expression() {
        let model = Model::from_json(&model_json("GC1")).unwrap();
        let named_coupling = Atom::var(symbol!("UFO::GC1"));
        let unknown_coupling = Atom::var(symbol!("UFO::GC_UNKNOWN"));
        let source = named_coupling.clone() * 2 + unknown_coupling.clone();

        let expanded = model.expand_couplings(&source);
        let expected = model.coupling("GC1").unwrap().expression.clone() * 2 + unknown_coupling;

        assert_eq!(expanded, expected);
        assert_eq!(
            source,
            named_coupling * 2 + Atom::var(symbol!("UFO::GC_UNKNOWN"))
        );
        assert_eq!(model.expand_couplings(&expanded), expanded);
    }

    #[test]
    fn restriction_updates_are_immutable_and_revalidated() {
        let model = Model::from_json(&model_json("GC1")).unwrap();
        let original_fingerprint = model.fingerprint();

        let restricted = model.with_restriction(Some("massless".to_owned())).unwrap();

        assert_eq!(model.restriction(), None);
        assert_eq!(model.fingerprint(), original_fingerprint);
        assert_eq!(restricted.restriction(), Some("massless"));
        assert_ne!(restricted.fingerprint(), original_fingerprint);
        assert_eq!(
            restricted.particle_id("s").unwrap(),
            model.particle_id("s").unwrap()
        );
    }

    #[test]
    fn rejects_unknown_vertex_couplings_with_context() {
        let error = Model::from_json(&model_json("missing")).unwrap_err();

        assert!(matches!(
            error,
            ModelError::Validation(ModelValidationError::UnknownReference {
                kind: EntityKind::VertexRule,
                target_kind: EntityKind::Coupling,
                ref name,
                ref target,
                ..
            }) if name == "V1" && target == "missing"
        ));
    }

    #[test]
    fn rejects_nonreciprocal_antiparticle_mappings() {
        let mut definition: ModelDefinition = serde_json::from_str(&model_json("GC1")).unwrap();
        let mut antiparticle = definition.particles[0].clone();
        antiparticle.pdg_code = -1;
        antiparticle.name = "s~".to_owned();
        antiparticle.antiname = "s~".to_owned();
        antiparticle.propagator = None;
        definition.particles[0].antiname = antiparticle.name.clone();
        definition.particles.push(antiparticle);

        assert!(matches!(
            Model::new(definition),
            Err(ModelError::Validation(
                ModelValidationError::NonReciprocalAntiparticle {
                    ref particle,
                    ref antiparticle,
                    ref actual_antiparticle,
                }
            )) if particle == "s" && antiparticle == "s~" && actual_antiparticle == "s~"
        ));
    }

    #[test]
    fn rejects_propagators_owned_by_another_particle() {
        let mut definition: ModelDefinition = serde_json::from_str(&model_json("GC1")).unwrap();
        let mut owner = definition.particles[0].clone();
        owner.pdg_code = 2;
        owner.name = "other".to_owned();
        owner.antiname = owner.name.clone();
        owner.propagator = None;
        definition.particles.push(owner);
        definition.propagators[0].particle = "other".to_owned();

        assert!(matches!(
            Model::new(definition),
            Err(ModelError::Validation(
                ModelValidationError::PropagatorParticleMismatch {
                    ref particle,
                    ref propagator,
                    ref propagator_particle,
                }
            )) if particle == "s" && propagator == "s_prop" && propagator_particle == "other"
        ));
    }

    #[test]
    fn rejects_lorentz_spins_that_disagree_with_vertex_particles() {
        let mut definition: ModelDefinition = serde_json::from_str(&model_json("GC1")).unwrap();
        definition.lorentz_structures[0].spins = vec![2, 2, 2];

        assert!(matches!(
            Model::new(definition),
            Err(ModelError::Validation(
                ModelValidationError::LorentzSpinMismatch {
                    ref vertex,
                    ref lorentz_structure,
                    ref particle_spins,
                    ref lorentz_spins,
                }
            )) if vertex == "V1"
                && lorentz_structure == "L1"
                && particle_spins.as_slice() == [1, 1, 1]
                && lorentz_spins.as_slice() == [2, 2, 2]
        ));
    }

    #[test]
    fn parameter_card_round_trips_and_applies_atomically() {
        let mut model = Model::from_json(&model_json("GC1")).unwrap();
        let card = ParameterCard::from_json(r#"{"mass":[2.5,0.25]}"#).unwrap();
        assert_eq!(card.to_json().unwrap(), r#"{"mass":[2.5,0.25]}"#);

        model.apply_parameter_card(&card).unwrap();
        assert_eq!(
            model.parameter("mass").unwrap().value,
            Some(ComplexValue::new(2.5, 0.25))
        );
        assert_eq!(model.parameter("double_mass").unwrap().value, None);
        assert_eq!(model.coupling("GC1").unwrap().value, None);

        let normalized = ParameterCard::from_iter([
            ("mass".to_owned(), ComplexValue::new(3.0, 0.0)),
            ("double_mass".to_owned(), ComplexValue::new(6.0, 0.0)),
        ]);
        model.apply_parameter_card(&normalized).unwrap();
        assert_eq!(
            model.parameter("mass").unwrap().value,
            Some(ComplexValue::new(3.0, 0.0))
        );
        assert_eq!(
            model.parameter("double_mass").unwrap().value,
            Some(ComplexValue::new(6.0, 0.0))
        );
        assert_eq!(model.coupling("GC1").unwrap().value, None);

        let invalid = ParameterCard::from_iter([
            ("mass".to_owned(), ComplexValue::new(7.0, 0.0)),
            ("unknown".to_owned(), ComplexValue::new(0.0, 0.0)),
        ]);
        let before_invalid = model.to_json().unwrap();
        assert!(model.apply_parameter_card(&invalid).is_err());
        assert_eq!(model.to_json().unwrap(), before_invalid);
    }

    #[test]
    fn rejects_duplicate_function_and_form_factor_names() {
        let model = Model::from_json(&model_json("GC1")).unwrap();
        let mut definition = model.clone().into_definition();
        definition.functions.push(definition.functions[0].clone());
        assert!(matches!(
            Model::new(definition),
            Err(ModelError::Validation(
                ModelValidationError::DuplicateName {
                    kind: EntityKind::ModelFunction,
                    ref name,
                }
            )) if name == "twice"
        ));

        let mut definition = model.into_definition();
        definition
            .form_factors
            .push(definition.form_factors[0].clone());
        assert!(matches!(
            Model::new(definition),
            Err(ModelError::Validation(
                ModelValidationError::DuplicateName {
                    kind: EntityKind::FormFactor,
                    ref name,
                }
            )) if name == "FF1"
        ));
    }

    #[test]
    fn recomputes_card_updates_through_an_owned_evaluator_request() {
        let mut model = Model::from_json(&model_json("GC1")).unwrap();
        let card = ParameterCard::from_iter([("mass".to_owned(), ComplexValue::new(2.5, 0.25))]);
        let mut evaluator = |request: EvaluationRequest| {
            assert_eq!(
                request.known_parameters["mass"],
                ComplexValue::new(2.5, 0.25)
            );
            assert_eq!(
                request.known_parameters["ZERO"],
                ComplexValue::new(0.0, 0.0)
            );
            assert_eq!(
                request.internal_parameters,
                [ModelExpression {
                    name: "double_mass".to_owned(),
                    expression: "2*mass".to_owned(),
                }]
            );
            assert_eq!(
                request.couplings,
                [ModelExpression {
                    name: "GC1".to_owned(),
                    expression: "double_mass".to_owned(),
                }]
            );
            assert_eq!(
                request.functions,
                [crate::EvaluationFunction {
                    name: "twice".to_owned(),
                    arguments: vec!["x".to_owned()],
                    expression: Some("2*x".to_owned()),
                }]
            );
            assert_eq!(
                request.form_factors,
                [crate::EvaluationFormFactor {
                    name: "FF1".to_owned(),
                    type_name: Some("complex".to_owned()),
                    value: Some("twice(P(1)^2)".to_owned()),
                }]
            );
            Ok::<_, std::convert::Infallible>(EvaluatedValues {
                internal_parameters: BTreeMap::from([(
                    "double_mass".to_owned(),
                    ComplexValue::new(5.0, 0.5),
                )]),
                couplings: BTreeMap::from([("GC1".to_owned(), ComplexValue::new(5.0, 0.5))]),
            })
        };

        model
            .apply_parameter_card_with(&card, &mut evaluator)
            .unwrap();

        assert_eq!(
            model.parameter("mass").unwrap().value,
            Some(ComplexValue::new(2.5, 0.25))
        );
        assert_eq!(
            model.parameter("double_mass").unwrap().value,
            Some(ComplexValue::new(5.0, 0.5))
        );
        assert_eq!(
            model.coupling("GC1").unwrap().value,
            Some(ComplexValue::new(5.0, 0.5))
        );
    }

    #[test]
    fn normalized_internal_card_values_override_recomputation() {
        let mut model = Model::from_json(&model_json("GC1")).unwrap();
        let card = ParameterCard::from_iter([
            ("mass".to_owned(), ComplexValue::new(3.0, 0.0)),
            ("double_mass".to_owned(), ComplexValue::new(99.0, 0.0)),
        ]);
        let mut evaluator = |request: EvaluationRequest| {
            assert_eq!(
                request.known_parameters["double_mass"],
                ComplexValue::new(99.0, 0.0)
            );
            assert!(request.internal_parameters.is_empty());
            Ok::<_, std::convert::Infallible>(EvaluatedValues {
                internal_parameters: BTreeMap::new(),
                couplings: BTreeMap::from([("GC1".to_owned(), ComplexValue::new(99.0, 0.0))]),
            })
        };

        model
            .apply_parameter_card_with(&card, &mut evaluator)
            .unwrap();

        assert_eq!(
            model.parameter("double_mass").unwrap().value,
            Some(ComplexValue::new(99.0, 0.0))
        );
        assert_eq!(
            model.coupling("GC1").unwrap().value,
            Some(ComplexValue::new(99.0, 0.0))
        );
    }

    #[test]
    fn failed_recomputation_leaves_the_model_unchanged() {
        let mut model = Model::from_json(&model_json("GC1")).unwrap();
        let original = model.to_json().unwrap();
        let card = ParameterCard::from_iter([("mass".to_owned(), ComplexValue::new(7.0, 0.0))]);
        let mut evaluator = |_request: EvaluationRequest| {
            Err::<EvaluatedValues, _>(std::io::Error::other("evaluation failed"))
        };

        assert!(matches!(
            model.apply_parameter_card_with(&card, &mut evaluator),
            Err(RecomputeError::Evaluator(_))
        ));
        assert_eq!(model.to_json().unwrap(), original);

        let mut incomplete = |_request: EvaluationRequest| {
            Ok::<_, std::convert::Infallible>(EvaluatedValues {
                internal_parameters: BTreeMap::from([(
                    "double_mass".to_owned(),
                    ComplexValue::new(14.0, 0.0),
                )]),
                couplings: BTreeMap::new(),
            })
        };
        assert!(matches!(
            model.apply_parameter_card_with(&card, &mut incomplete),
            Err(RecomputeError::MissingValue {
                kind: EntityKind::Coupling,
                ref name,
            }) if name == "GC1"
        ));
        assert_eq!(model.to_json().unwrap(), original);
    }
}
