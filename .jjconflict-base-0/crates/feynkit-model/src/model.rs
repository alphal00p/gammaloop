use std::{
    collections::{BTreeMap, BTreeSet, HashMap},
    fs,
    path::Path,
};

use serde::{Deserialize, Deserializer, Serialize, Serializer, de};

use crate::{
    ComplexValue, EntityKind, EvaluatedValues, EvaluationRequest, ModelError, ModelEvaluator,
    ModelExpression, ModelValidationError, ParameterCard, RecomputeError,
};

const fn default_true() -> bool {
    true
}

/// The UFO convention stores twice the spin plus one (scalar = 1, spinor = 2,
/// vector = 3, and tensor = 5).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Particle {
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

impl Particle {
    pub fn is_antiparticle(&self) -> bool {
        self.pdg_code < 0
    }

    pub fn is_self_antiparticle(&self) -> bool {
        self.name == self.antiname
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

    pub fn is_massless(&self) -> bool {
        self.mass == "ZERO"
    }

    pub fn is_qcd_charged(&self) -> bool {
        self.color != 1
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Propagator {
    pub name: String,
    pub particle: String,
    pub numerator: String,
    pub denominator: String,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct LorentzStructure {
    pub name: String,
    pub spins: Vec<i64>,
    pub structure: String,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Coupling {
    pub name: String,
    pub expression: String,
    #[serde(with = "vectorize")]
    pub orders: BTreeMap<String, usize>,
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
pub struct Parameter {
    pub name: String,
    pub lhablock: Option<String>,
    pub lhacode: Option<Vec<usize>>,
    pub nature: ParameterNature,
    pub parameter_type: ParameterType,
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
pub struct VertexRule {
    pub name: String,
    pub particles: Vec<String>,
    pub color_structures: Vec<String>,
    pub lorentz_structures: Vec<String>,
    pub couplings: Vec<Vec<Option<String>>>,
}

impl VertexRule {
    /// Return the largest power of every coupling order represented by an
    /// entry in this vertex's color/Lorentz coupling matrix.
    pub fn coupling_orders(&self, model: &Model) -> Result<BTreeMap<String, usize>, ModelError> {
        let mut orders: BTreeMap<String, usize> = BTreeMap::new();
        for coupling_name in self.couplings.iter().flatten().flatten() {
            for (name, power) in &model.coupling(coupling_name)?.orders {
                orders
                    .entry(name.clone())
                    .and_modify(|current| *current = (*current).max(*power))
                    .or_insert(*power);
            }
        }
        Ok(orders)
    }
}

/// A callable function declared by a UFO model.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ModelFunction {
    pub name: String,
    #[serde(default)]
    pub arguments: Vec<String>,
    pub expression: Option<String>,
}

/// Metadata for a UFO form factor.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ModelFormFactor {
    pub name: String,
    #[serde(rename = "type")]
    pub type_name: Option<String>,
    pub value: Option<String>,
}

/// Serializable data from which a validated [`Model`] is constructed.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ModelDefinition {
    pub name: String,
    pub restriction: Option<String>,
    pub orders: Vec<Order>,
    pub parameters: Vec<Parameter>,
    pub particles: Vec<Particle>,
    pub propagators: Vec<Propagator>,
    pub lorentz_structures: Vec<LorentzStructure>,
    pub couplings: Vec<Coupling>,
    pub vertex_rules: Vec<VertexRule>,
    #[serde(default)]
    pub functions: Vec<ModelFunction>,
    #[serde(default)]
    pub form_factors: Vec<ModelFormFactor>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct ParticleId(usize);

impl ParticleId {
    pub fn index(self) -> usize {
        self.0
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct VertexRuleId(usize);

impl VertexRuleId {
    pub fn index(self) -> usize {
        self.0
    }
}

#[derive(Clone, Debug, Default)]
struct ModelIndexes {
    orders: HashMap<String, usize>,
    parameters: HashMap<String, usize>,
    particles: HashMap<String, ParticleId>,
    particles_by_pdg: HashMap<i64, ParticleId>,
    propagators: HashMap<String, usize>,
    lorentz_structures: HashMap<String, usize>,
    couplings: HashMap<String, usize>,
    vertex_rules: HashMap<String, VertexRuleId>,
    functions: HashMap<String, usize>,
    form_factors: HashMap<String, usize>,
}

/// A validated physics model with fallible, indexed lookups.
#[derive(Clone, Debug)]
pub struct Model {
    definition: ModelDefinition,
    indexes: ModelIndexes,
}

impl Model {
    pub fn new(definition: ModelDefinition) -> Result<Self, ModelError> {
        let indexes = ModelIndexes::build(&definition)?;
        Ok(Self {
            definition,
            indexes,
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
        Ok(serde_json::to_string(&self.definition)?)
    }

    pub fn to_json_pretty(&self) -> Result<String, ModelError> {
        Ok(serde_json::to_string_pretty(&self.definition)?)
    }

    pub fn write_json(&self, path: impl AsRef<Path>) -> Result<(), ModelError> {
        let path = path.as_ref();
        fs::write(path, self.to_json_pretty()?).map_err(|source| ModelError::Write {
            path: path.to_path_buf(),
            source,
        })
    }

    pub fn definition(&self) -> &ModelDefinition {
        &self.definition
    }

    pub fn into_definition(self) -> ModelDefinition {
        self.definition
    }

    pub fn name(&self) -> &str {
        &self.definition.name
    }

    pub fn restriction(&self) -> Option<&str> {
        self.definition.restriction.as_deref()
    }

    pub fn orders(&self) -> &[Order] {
        &self.definition.orders
    }

    pub fn parameters(&self) -> &[Parameter] {
        &self.definition.parameters
    }

    pub fn particles(&self) -> &[Particle] {
        &self.definition.particles
    }

    pub fn propagators(&self) -> &[Propagator] {
        &self.definition.propagators
    }

    pub fn lorentz_structures(&self) -> &[LorentzStructure] {
        &self.definition.lorentz_structures
    }

    pub fn couplings(&self) -> &[Coupling] {
        &self.definition.couplings
    }

    pub fn vertex_rules(&self) -> &[VertexRule] {
        &self.definition.vertex_rules
    }

    pub fn functions(&self) -> &[ModelFunction] {
        &self.definition.functions
    }

    pub fn form_factors(&self) -> &[ModelFormFactor] {
        &self.definition.form_factors
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
        self.definition
            .particles
            .get(id.0)
            .ok_or_else(|| self.not_found(EntityKind::Particle, &id.0.to_string()))
    }

    pub fn antiparticle(&self, particle: &Particle) -> Result<&Particle, ModelError> {
        self.particle(&particle.antiname)
    }

    pub fn parameter(&self, name: &str) -> Result<&Parameter, ModelError> {
        self.lookup(EntityKind::Parameter, name, &self.indexes.parameters)
            .map(|index| &self.definition.parameters[index])
    }

    pub fn coupling(&self, name: &str) -> Result<&Coupling, ModelError> {
        self.lookup(EntityKind::Coupling, name, &self.indexes.couplings)
            .map(|index| &self.definition.couplings[index])
    }

    pub fn order(&self, name: &str) -> Result<&Order, ModelError> {
        self.lookup(EntityKind::Order, name, &self.indexes.orders)
            .map(|index| &self.definition.orders[index])
    }

    pub fn propagator(&self, name: &str) -> Result<&Propagator, ModelError> {
        self.lookup(EntityKind::Propagator, name, &self.indexes.propagators)
            .map(|index| &self.definition.propagators[index])
    }

    pub fn lorentz_structure(&self, name: &str) -> Result<&LorentzStructure, ModelError> {
        self.lookup(
            EntityKind::LorentzStructure,
            name,
            &self.indexes.lorentz_structures,
        )
        .map(|index| &self.definition.lorentz_structures[index])
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
        self.definition
            .vertex_rules
            .get(id.0)
            .ok_or_else(|| self.not_found(EntityKind::VertexRule, &id.0.to_string()))
    }

    pub fn function(&self, name: &str) -> Result<&ModelFunction, ModelError> {
        self.lookup(EntityKind::ModelFunction, name, &self.indexes.functions)
            .map(|index| &self.definition.functions[index])
    }

    pub fn form_factor(&self, name: &str) -> Result<&ModelFormFactor, ModelError> {
        self.lookup(EntityKind::FormFactor, name, &self.indexes.form_factors)
            .map(|index| &self.definition.form_factors[index])
    }

    pub fn default_parameter_card(&self) -> Result<ParameterCard, ModelError> {
        self.definition
            .parameters
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
        for parameter in &mut self.definition.parameters {
            if parameter.nature == ParameterNature::Internal && parameter.expression.is_some() {
                parameter.value = None;
            }
        }
        for coupling in &mut self.definition.couplings {
            coupling.value = None;
        }
        for (index, value) in updates {
            self.definition.parameters[index].value = Some(value);
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

        for parameter in &mut self.definition.parameters {
            if parameter.nature == ParameterNature::Internal
                && parameter.expression.is_some()
                && !overrides.contains(&parameter.name)
            {
                parameter.value = evaluated.internal_parameters.get(&parameter.name).copied();
            }
        }
        for coupling in &mut self.definition.couplings {
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
                    let parameter = &updated.definition.parameters[*index];
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
            .definition
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
            .definition
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
                        expression: expression.clone(),
                    })
            })
            .collect();
        let couplings = self
            .definition
            .couplings
            .iter()
            .map(|coupling| ModelExpression {
                name: coupling.name.clone(),
                expression: coupling.expression.clone(),
            })
            .collect();
        EvaluationRequest {
            known_parameters,
            internal_parameters,
            couplings,
            functions: self.definition.functions.clone(),
            form_factors: self.definition.form_factors.clone(),
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
        for parameter in self.definition.parameters.iter().filter(|parameter| {
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
            let parameter = &self.definition.parameters[*index];
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
        for coupling in &self.definition.couplings {
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

    fn lookup(
        &self,
        kind: EntityKind,
        name: &str,
        index: &HashMap<String, usize>,
    ) -> Result<usize, ModelError> {
        index
            .get(name)
            .copied()
            .ok_or_else(|| self.not_found(kind, name))
    }

    fn not_found(&self, kind: EntityKind, key: &str) -> ModelError {
        ModelError::NotFound {
            model: self.definition.name.clone(),
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
        self.definition.serialize(serializer)
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

impl TryFrom<ModelDefinition> for Model {
    type Error = ModelError;

    fn try_from(definition: ModelDefinition) -> Result<Self, Self::Error> {
        Self::new(definition)
    }
}

impl From<Model> for ModelDefinition {
    fn from(model: Model) -> Self {
        model.definition
    }
}

impl ModelIndexes {
    fn build(definition: &ModelDefinition) -> Result<Self, ModelValidationError> {
        if definition.name.trim().is_empty() {
            return Err(ModelValidationError::EmptyModelName);
        }

        let orders = named_index(EntityKind::Order, &definition.orders, |value| &value.name)?;
        let parameters = named_index(EntityKind::Parameter, &definition.parameters, |value| {
            &value.name
        })?;
        let particles = typed_index(
            EntityKind::Particle,
            &definition.particles,
            |value| &value.name,
            ParticleId,
        )?;
        let propagators = named_index(EntityKind::Propagator, &definition.propagators, |value| {
            &value.name
        })?;
        let lorentz_structures = named_index(
            EntityKind::LorentzStructure,
            &definition.lorentz_structures,
            |value| &value.name,
        )?;
        let couplings = named_index(EntityKind::Coupling, &definition.couplings, |value| {
            &value.name
        })?;
        let vertex_rules = typed_index(
            EntityKind::VertexRule,
            &definition.vertex_rules,
            |value| &value.name,
            VertexRuleId,
        )?;
        let functions = named_index(EntityKind::ModelFunction, &definition.functions, |value| {
            &value.name
        })?;
        let form_factors =
            named_index(EntityKind::FormFactor, &definition.form_factors, |value| {
                &value.name
            })?;

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

fn named_index<T>(
    kind: EntityKind,
    values: &[T],
    name: impl Fn(&T) -> &str,
) -> Result<HashMap<String, usize>, ModelValidationError> {
    let mut result = HashMap::with_capacity(values.len());
    for (index, value) in values.iter().enumerate() {
        let name = name(value);
        if result.insert(name.to_owned(), index).is_some() {
            return Err(ModelValidationError::DuplicateName {
                kind,
                name: name.to_owned(),
            });
        }
    }
    Ok(result)
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

fn validate_parameters(parameters: &[Parameter]) -> Result<(), ModelValidationError> {
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
    couplings: &[Coupling],
    orders: &HashMap<String, usize>,
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
    particles: &[Particle],
    parameters: &HashMap<String, usize>,
    particle_names: &HashMap<String, ParticleId>,
    propagator_definitions: &[Propagator],
    propagators: &HashMap<String, usize>,
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
            let propagator_definition = &propagator_definitions[propagator_index];
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
    propagators: &[Propagator],
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
    vertices: &[VertexRule],
    particle_definitions: &[Particle],
    particles: &HashMap<String, ParticleId>,
    lorentz_definitions: &[LorentzStructure],
    lorentz_structures: &HashMap<String, usize>,
    couplings: &HashMap<String, usize>,
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
            let lorentz_spins = &lorentz_definitions[lorentz_index].spins;
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
            model.particle("s").unwrap().propagator.as_deref(),
            Some("s_prop")
        );
        assert_eq!(model.coupling("GC1").unwrap().orders["QED"], 1);
        assert_eq!(model.function("twice").unwrap().arguments, ["x"]);
        assert_eq!(
            model.form_factor("FF1").unwrap().type_name.as_deref(),
            Some("complex")
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
            ("unknown".to_owned(), ComplexValue::ZERO),
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
            assert_eq!(request.known_parameters["ZERO"], ComplexValue::ZERO);
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
                [ModelFunction {
                    name: "twice".to_owned(),
                    arguments: vec!["x".to_owned()],
                    expression: Some("2*x".to_owned()),
                }]
            );
            assert_eq!(
                request.form_factors,
                [ModelFormFactor {
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
