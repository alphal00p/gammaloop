use std::{
    collections::{BTreeMap, BTreeSet, HashMap},
    fmt::{Display, Formatter},
    ops::{Deref, DerefMut},
    path::Path,
};

use color_eyre::Report;
use eyre::eyre;
use idenso::{
    dirac::AGS,
    representations::{Bispinor, ColorAdjoint, ColorFundamental, ColorSextet},
};
use itertools::Itertools;
use linnet::half_edge::involution::{EdgeIndex, Flow};
use rand::{Rng, SeedableRng, rngs::SmallRng};
use serde::{Deserialize, Serialize, de::DeserializeOwned};
use spenso::algebra::complex::Complex;
use spenso::network::library::symbolic::ETS;
use spenso::structure::{
    IndexLess, OrderedStructure, PermutedStructure,
    representation::{Euclidean, LibraryRep, Minkowski, RepName},
    slot::{DummyAind, IsAbstractSlot, Slot},
};
use spenso::tensors::{
    data::{DataTensor, DenseTensor, SetTensorData, SparseTensor},
    parametric::ParamTensor,
};
use symbolica::domains::rational::Fraction;
use symbolica::prelude::*;
use symbolica_utils::Replaces;

use crate::{
    momentum::Helicity,
    numerator::aind::Aind,
    settings::global::VectorPolarizationSumGauge,
    utils::{F, GS, W_, serde_utils::SmartSerde, symbolica_ext::DOD},
};

pub use feynkit_model::{
    Coupling, CouplingId, LorentzStructure, LorentzStructureId, Model, ModelFingerprint,
    ModelFormFactor, ModelFormFactorId, ModelFunction, ModelFunctionId, Order, OrderId, Parameter,
    ParameterId, ParameterNature, ParameterType, Particle, ParticleId, Propagator, PropagatorId,
    VertexRule, VertexRuleId,
};

/// Parameter-card representation used at GammaLoop's precision boundary.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SerializableInputParamCard<T> {
    pub data: BTreeMap<String, (T, T)>,
}

impl<T> SerializableInputParamCard<T>
where
    T: From<f64> + Clone + Serialize + DeserializeOwned,
{
    pub(crate) fn from_file(file_path: impl AsRef<Path>) -> Result<Self, Report> {
        let data = SmartSerde::from_file(file_path, "model")?;
        Ok(Self { data })
    }

    pub fn from_str(source: String, format: &str) -> Result<Self, Report> {
        let data = SmartSerde::from_str(source, format, "model_parameters")?;
        Ok(Self { data })
    }

    pub fn from_input_param_card(card: &InputParamCard<T>) -> Self {
        Self {
            data: card
                .data
                .iter()
                .map(|(name, value)| {
                    (
                        name.namespaceless_string().to_owned(),
                        (value.re.clone(), value.im.clone()),
                    )
                })
                .collect(),
        }
    }

    pub fn to_file<P: AsRef<Path>>(&self, path: P, overwrite: bool) -> Result<(), Report> {
        SmartSerde::to_file(&self.data, path, overwrite)
    }
}

#[derive(Debug, Clone)]
pub struct InputParamCard<T>
where
    T: From<f64> + Clone + Serialize + DeserializeOwned,
{
    data: HashMap<UFOSymbol, Complex<T>>,
}

impl<T> Default for InputParamCard<T>
where
    T: From<f64> + Clone + Serialize + DeserializeOwned,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T> InputParamCard<T>
where
    T: From<f64> + Clone + Serialize + DeserializeOwned,
{
    pub fn new() -> Self {
        Self {
            data: HashMap::new(),
        }
    }

    pub fn insert(&mut self, symbol: UFOSymbol, value: Complex<T>) -> Option<Complex<T>> {
        self.data.insert(symbol, value)
    }

    pub fn get(&self, symbol: &UFOSymbol) -> Option<&Complex<T>> {
        self.data.get(symbol)
    }

    pub fn get_mut(&mut self, symbol: &UFOSymbol) -> Option<&mut Complex<T>> {
        self.data.get_mut(symbol)
    }

    pub fn remove(&mut self, symbol: &UFOSymbol) -> Option<Complex<T>> {
        self.data.remove(symbol)
    }

    pub fn iter(&self) -> impl Iterator<Item = (&UFOSymbol, &Complex<T>)> {
        self.data.iter()
    }

    pub fn to_serializable(&self) -> SerializableInputParamCard<T> {
        SerializableInputParamCard::from_input_param_card(self)
    }

    pub fn from_serializable(card: &SerializableInputParamCard<T>) -> Self {
        Self {
            data: card
                .data
                .iter()
                .map(|(name, (re, im))| {
                    (UFOSymbol::from(name), Complex::new(re.clone(), im.clone()))
                })
                .collect(),
        }
    }

    pub fn from_file(file_path: impl AsRef<Path>) -> Result<Self, Report> {
        Ok(Self::from_serializable(
            &SerializableInputParamCard::from_file(file_path)?,
        ))
    }

    pub fn from_str(source: String, format: &str) -> Result<Self, Report> {
        Ok(Self::from_serializable(
            &SerializableInputParamCard::from_str(source, format)?,
        ))
    }

    pub fn to_file<P: AsRef<Path>>(&self, path: P, overwrite: bool) -> Result<(), Report> {
        self.to_serializable().to_file(path, overwrite)
    }
}

impl<T> Deref for InputParamCard<T>
where
    T: From<f64> + Clone + Serialize + DeserializeOwned,
{
    type Target = HashMap<UFOSymbol, Complex<T>>;

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl<T> DerefMut for InputParamCard<T>
where
    T: From<f64> + Clone + Serialize + DeserializeOwned,
{
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data
    }
}

impl From<feynkit_model::ParameterCard> for InputParamCard<F<f64>> {
    fn from(card: feynkit_model::ParameterCard) -> Self {
        Self {
            data: card
                .into_iter()
                .map(|(name, value)| {
                    (
                        UFOSymbol::from(name),
                        Complex::new(F(value.re), F(value.im)),
                    )
                })
                .collect(),
        }
    }
}

impl From<&InputParamCard<F<f64>>> for feynkit_model::ParameterCard {
    fn from(card: &InputParamCard<F<f64>>) -> Self {
        card.iter()
            .map(|(name, value)| {
                (
                    name.namespaceless_string().to_owned(),
                    feynkit_model::ComplexValue::new(value.re.0, value.im.0),
                )
            })
            .collect()
    }
}

/// A UFO symbol interned in Symbolica's `UFO` namespace.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct UFOSymbol(pub Symbol);

impl UFOSymbol {
    pub fn zero() -> Self {
        Self(symbol!("UFO::ZERO"))
    }

    pub fn is_zero(self) -> bool {
        self.0 == Self::zero().0
    }

    pub fn namespaceless_string(&self) -> &str {
        self.0.get_stripped_name()
    }
}

impl Display for UFOSymbol {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        Display::fmt(&self.0, formatter)
    }
}

impl From<UFOSymbol> for Atom {
    fn from(value: UFOSymbol) -> Self {
        Atom::var(value.0)
    }
}

impl<T: AsRef<str>> From<T> for UFOSymbol {
    fn from(name: T) -> Self {
        if name.as_ref() == "ZERO" {
            Self::zero()
        } else {
            Self(symbol!(&format!("UFO::{}", name.as_ref())))
        }
    }
}

/// GammaLoop operations over the canonical FeynKit model.
pub trait ModelGammaLoopExt {
    fn apply_param_card(&mut self, card: &InputParamCard<F<f64>>) -> Result<(), Report>;
    fn recompute_dependents(&mut self) -> Result<(), Report>;
    fn simplify(&mut self, card: &mut InputParamCard<F<f64>>) -> Result<(), Report>;
    fn default_param_card(&self) -> InputParamCard<F<f64>>;
    fn get_particle(&self, name: impl AsRef<str>) -> ParticleId;
    fn try_get_particle(&self, name: impl AsRef<str>) -> Result<ParticleId, Report>;
    fn get_particle_from_pdg(&self, pdg: isize) -> ParticleId;
    fn try_get_particle_from_pdg(&self, pdg: isize) -> Result<ParticleId, Report>;
    fn get_propagator(&self, name: impl AsRef<str>) -> PropagatorId;
    fn get_propagator_for_particle(&self, particle: ParticleId) -> PropagatorId;
    fn get_parameter_opt(&self, name: impl AsRef<str>) -> Option<&Parameter>;
    fn get_parameter(&self, name: impl AsRef<str>) -> &Parameter;
    fn get_coupling(&self, name: impl AsRef<str>) -> &Coupling;
    fn get_vertex_rule(&self, name: impl AsRef<str>) -> VertexRuleId;
    fn vertex_rules_for_particles(&self, particles: &[ParticleId]) -> Vec<VertexRuleId>;
    fn unresolved_particles(&self) -> HashMap<String, BTreeSet<ParticleId>>;
    fn contains_symbol(&self, symbol: &UFOSymbol) -> bool;
    fn get_symbol_value(&self, symbol: UFOSymbol) -> Option<Complex<F<f64>>>;
    fn generate_params(&self) -> Vec<Atom>;
    fn apply_coupling_replacement_rules(&self, atom: &Atom) -> Atom;
    fn apply_parameter_replacement_rules(&self, atom: &Atom) -> Atom;
    fn get_description(
        &self,
        show_particles: bool,
        show_parameters: bool,
        show_vertices: bool,
        show_couplings: bool,
    ) -> String;
    fn generate_edge_style_template(
        &self,
        template_path: impl AsRef<std::path::Path>,
    ) -> Result<(), std::io::Error>;
    fn is_empty(&self) -> bool;
}

impl ModelGammaLoopExt for Model {
    fn apply_param_card(&mut self, card: &InputParamCard<F<f64>>) -> Result<(), Report> {
        self.apply_parameter_card(&card.into())?;
        self.recompute_dependents()
    }

    fn recompute_dependents(&mut self) -> Result<(), Report> {
        let replacements = self
            .couplings()
            .iter()
            .map(|coupling| {
                let symbol = UFOSymbol::from(&coupling.name).0;
                Replacement::new(Atom::var(symbol).to_pattern(), function!(symbol))
            })
            .chain(self.parameters().iter().map(|parameter| {
                let symbol = UFOSymbol::from(&parameter.name).0;
                Replacement::new(Atom::var(symbol).to_pattern(), function!(symbol))
            }))
            .collect_vec();
        let mut function_map = FunctionMap::new();
        let mut expressions = Vec::new();
        for coupling in self.couplings() {
            let symbol = UFOSymbol::from(&coupling.name).0;
            expressions.push(function!(symbol));
            function_map
                .add_function::<Symbol, Symbol>(
                    symbol,
                    Vec::new(),
                    coupling.expression.replace_multiple(&replacements),
                )
                .map_err(|error| eyre!("{error}"))?;
        }

        let mut parameters = Vec::new();
        let mut parameter_values = Vec::new();
        let mut dependent_parameters = Vec::new();
        for (index, parameter) in self.parameters().iter().enumerate() {
            let symbol = UFOSymbol::from(&parameter.name).0;
            let key = function!(symbol);
            match parameter.nature {
                ParameterNature::External => {
                    parameters.push(key);
                    let value = parameter.value.ok_or_else(|| {
                        eyre!("External parameter '{}' has no value", parameter.name)
                    })?;
                    parameter_values.push(Complex::new(F(value.re), F(value.im)));
                }
                ParameterNature::Internal if parameter.name == "ZERO" => {}
                ParameterNature::Internal => {
                    if let Some(body) = &parameter.expression {
                        expressions.push(key.clone());
                        dependent_parameters.push(ParameterId::from_index(index));
                        function_map
                            .add_function::<Symbol, Symbol>(
                                symbol,
                                Vec::new(),
                                body.replace_multiple(&replacements),
                            )
                            .map_err(|error| eyre!("{error}"))?;
                    } else {
                        let value = parameter.value.ok_or_else(|| {
                            eyre!(
                                "Internal parameter '{}' has neither expression nor value",
                                parameter.name
                            )
                        })?;
                        let value = symbolica::domains::float::Complex::new(
                            Fraction::<IntegerRing>::try_from(value.re).unwrap(),
                            Fraction::<IntegerRing>::try_from(value.im).unwrap(),
                        );
                        function_map
                            .add_aliases([(key, Atom::num(value))])
                            .map_err(|error| eyre!("{error}"))?;
                    }
                }
            }
        }
        function_map
            .add_aliases([(
                Atom::var(Symbol::PI),
                Atom::num(symbolica::domains::float::Complex::new(
                    Rational::try_from(std::f64::consts::PI).unwrap(),
                    Rational::zero(),
                )),
            )])
            .map_err(|error| eyre!("{error}"))?;

        let evaluator = AtomView::to_eval_tree_multiple(&expressions, &function_map, &parameters)?
            .linearize(
                &OptimizationSettings::new()
                    .cpe_iterations(Some(1))
                    .verbose(false),
            );
        let mut evaluator =
            evaluator.map_coeff(&|value| Complex::new(F(value.re.to_f64()), F(value.im.to_f64())));
        let mut values = vec![Complex::new(F(0.0), F(0.0)); expressions.len()];
        evaluator.evaluate(&parameter_values, &mut values);
        let coupling_count = self.couplings().len();
        let evaluated = feynkit_model::EvaluatedValues {
            couplings: self
                .couplings()
                .iter()
                .zip(&values[..coupling_count])
                .map(|(coupling, value)| {
                    (
                        coupling.name.clone(),
                        feynkit_model::ComplexValue::new(value.re.0, value.im.0),
                    )
                })
                .collect(),
            internal_parameters: dependent_parameters
                .into_iter()
                .zip(&values[coupling_count..])
                .map(|(parameter, value)| {
                    (
                        self.parameter_by_id(parameter).unwrap().name.clone(),
                        feynkit_model::ComplexValue::new(value.re.0, value.im.0),
                    )
                })
                .collect(),
        };
        self.recompute_with(&mut |_| Ok::<_, std::convert::Infallible>(evaluated.clone()))
            .map_err(|error| eyre!(error))?;
        Ok(())
    }

    fn simplify(&mut self, card: &mut InputParamCard<F<f64>>) -> Result<(), Report> {
        self.apply_param_card(card)?;
        card.data
            .retain(|_, value| value.re != F(0.0) || value.im != F(0.0));

        self.simplify_zero_values()?;
        Ok(())
    }

    fn default_param_card(&self) -> InputParamCard<F<f64>> {
        self.default_parameter_card()
            .expect("validated models have complete external parameter cards")
            .into()
    }

    fn get_particle(&self, name: impl AsRef<str>) -> ParticleId {
        self.try_get_particle(name).unwrap()
    }

    fn try_get_particle(&self, name: impl AsRef<str>) -> Result<ParticleId, Report> {
        self.particle_id(name.as_ref()).map_err(Into::into)
    }

    fn get_particle_from_pdg(&self, pdg: isize) -> ParticleId {
        self.try_get_particle_from_pdg(pdg).unwrap()
    }

    fn try_get_particle_from_pdg(&self, pdg: isize) -> Result<ParticleId, Report> {
        let pdg = i64::try_from(pdg).map_err(|_| eyre!("PDG code {pdg} is out of range"))?;
        self.particle_id_by_pdg(pdg).map_err(Into::into)
    }

    fn get_propagator(&self, name: impl AsRef<str>) -> PropagatorId {
        self.propagator_id(name.as_ref()).unwrap()
    }

    fn get_propagator_for_particle(&self, particle: ParticleId) -> PropagatorId {
        if let Some(propagator) = particle.resolve(self).propagator {
            return propagator;
        }
        self.propagators()
            .iter()
            .position(|propagator| propagator.particle == particle)
            .and_then(|index| self.propagator_id(&self.propagators()[index].name).ok())
            .unwrap_or_else(|| {
                panic!(
                    "Propagator for particle '{}' not found in model '{}'",
                    particle.resolve(self).name,
                    self.name()
                )
            })
    }

    fn get_parameter_opt(&self, name: impl AsRef<str>) -> Option<&Parameter> {
        self.parameter(name.as_ref()).ok()
    }

    fn get_parameter(&self, name: impl AsRef<str>) -> &Parameter {
        self.parameter(name.as_ref()).unwrap()
    }

    fn get_coupling(&self, name: impl AsRef<str>) -> &Coupling {
        self.coupling(name.as_ref()).unwrap()
    }

    fn get_vertex_rule(&self, name: impl AsRef<str>) -> VertexRuleId {
        self.vertex_rule_id(name.as_ref()).unwrap()
    }

    fn vertex_rules_for_particles(&self, particles: &[ParticleId]) -> Vec<VertexRuleId> {
        let mut requested = particles.to_vec();
        requested.sort_unstable();

        self.vertex_rules()
            .iter()
            .filter(|rule| {
                let mut rule_particles = rule.particles.clone();
                rule_particles.sort_unstable();
                rule_particles == requested
            })
            .map(|rule| self.vertex_rule_id(&rule.name).unwrap())
            .collect()
    }

    fn unresolved_particles(&self) -> HashMap<String, BTreeSet<ParticleId>> {
        let mut unresolved = HashMap::<String, BTreeSet<ParticleId>>::new();
        for vertex_rule in self.vertex_rules() {
            let massless = vertex_rule
                .particles
                .iter()
                .copied()
                .filter(|particle| particle.resolve(self).is_massless(self))
                .collect::<BTreeSet<_>>();
            for order in vertex_rule.coupling_orders(self).into_keys() {
                unresolved.entry(order).or_default().extend(&massless);
            }
        }
        unresolved
    }

    fn contains_symbol(&self, symbol: &UFOSymbol) -> bool {
        self.parameter(symbol.namespaceless_string()).is_ok()
            || self.coupling(symbol.namespaceless_string()).is_ok()
    }

    fn get_symbol_value(&self, symbol: UFOSymbol) -> Option<Complex<F<f64>>> {
        self.coupling(symbol.namespaceless_string())
            .ok()
            .and_then(|coupling| coupling.value)
            .or_else(|| {
                self.parameter(symbol.namespaceless_string())
                    .ok()
                    .and_then(|parameter| parameter.value)
            })
            .map(|value| Complex::new(F(value.re), F(value.im)))
    }

    fn generate_params(&self) -> Vec<Atom> {
        self.couplings()
            .iter()
            .filter(|coupling| coupling.value.is_some())
            .map(|coupling| Atom::from(UFOSymbol::from(&coupling.name)))
            .chain(
                self.parameters()
                    .iter()
                    .filter(|parameter| parameter.value.is_some())
                    .map(|parameter| Atom::from(UFOSymbol::from(&parameter.name))),
            )
            .collect()
    }

    fn apply_coupling_replacement_rules(&self, atom: &Atom) -> Atom {
        let replacements = self
            .couplings()
            .iter()
            .map(|coupling| {
                Replacement::new(
                    Atom::from(UFOSymbol::from(&coupling.name)).to_pattern(),
                    coupling.expression.clone(),
                )
            })
            .collect::<Vec<_>>();
        atom.replace_multiple(&replacements)
    }

    fn apply_parameter_replacement_rules(&self, atom: &Atom) -> Atom {
        let replacements = self
            .parameters()
            .iter()
            .filter_map(|parameter| {
                parameter.expression.as_ref().map(|expression| {
                    Replacement::new(
                        Atom::from(UFOSymbol::from(&parameter.name)).to_pattern(),
                        expression.clone(),
                    )
                })
            })
            .collect::<Vec<_>>();
        atom.replace_multiple(&replacements)
    }

    fn get_description(
        &self,
        show_particles: bool,
        show_parameters: bool,
        show_vertices: bool,
        show_couplings: bool,
    ) -> String {
        use std::fmt::Write as _;

        let restriction = self.restriction().unwrap_or("None");
        let orders = if self.orders().is_empty() {
            "None".to_owned()
        } else {
            self.orders()
                .iter()
                .map(|order| {
                    format!(
                        "{} (expansion order: {}, hierarchy: {})",
                        order.name, order.expansion_order, order.hierarchy
                    )
                })
                .join(", ")
        };
        let mut output = format!(
            "\nModel name                    : {}\nRestriction                   : {}\nCoupling orders               : {}\n\n{} particles",
            self.name(),
            restriction,
            orders,
            self.particles().len()
        );

        if show_particles {
            for particle in self.particles() {
                let mass = self
                    .parameter_by_id(particle.mass)
                    .expect("validated particle mass references resolve");
                let width = self
                    .parameter_by_id(particle.width)
                    .expect("validated particle width references resolve");
                let mass_value = mass
                    .value
                    .map_or_else(|| "None".to_owned(), |value| format!("{:.6}", value.re));
                let width_value = width
                    .value
                    .map_or_else(|| "None".to_owned(), |value| format!("{:.6}", value.re));
                let _ = write!(
                    output,
                    "\n  {:<16} PDG {:+8}  mass {}={}  width {}={}",
                    particle.name,
                    particle.pdg_code,
                    mass.name,
                    mass_value,
                    width.name,
                    width_value
                );
            }
        } else {
            output.push_str(" [hidden]");
        }

        let _ = write!(output, "\n\n{} parameters", self.parameters().len());
        if show_parameters {
            let mut parameters = self.parameters().iter().collect_vec();
            parameters.sort_by_key(|parameter| {
                (
                    !matches!(parameter.nature, ParameterNature::External),
                    parameter.name.as_str(),
                )
            });
            for parameter in parameters {
                let value = parameter.value.map_or_else(String::new, |value| {
                    if value.im == 0.0 {
                        format!("{:.6}", value.re)
                    } else {
                        format!("{:.6}{:+.6}i", value.re, value.im)
                    }
                });
                let expression = parameter
                    .expression
                    .as_ref()
                    .map(ToString::to_string)
                    .unwrap_or_default();
                let _ = write!(
                    output,
                    "\n  {:<20} {:?} {:?}  {}{}{}",
                    parameter.name,
                    parameter.nature,
                    parameter.parameter_type,
                    value,
                    if !value.is_empty() && !expression.is_empty() {
                        " = "
                    } else {
                        ""
                    },
                    expression
                );
            }
        } else {
            output.push_str(" [hidden]");
        }

        let _ = write!(output, "\n\n{} vertices", self.vertex_rules().len());
        if show_vertices {
            for rule in self.vertex_rules() {
                let particles = rule
                    .particles
                    .iter()
                    .map(|particle| {
                        self.particle_by_id(*particle)
                            .expect("validated vertex particle references resolve")
                            .name
                            .as_str()
                    })
                    .join(", ");
                let _ = write!(output, "\n  {:<20} {}", rule.name, particles);
            }
        } else {
            output.push_str(" [hidden]");
        }

        let _ = write!(output, "\n\n{} couplings", self.couplings().len());
        if show_couplings {
            for coupling in self.couplings() {
                let orders = coupling
                    .orders
                    .iter()
                    .map(|(name, power)| format!("{name}={power}"))
                    .join(" ");
                let _ = write!(
                    output,
                    "\n  {:<20} {:<16} {}",
                    coupling.name, coupling.expression, orders
                );
            }
        } else {
            output.push_str(" [hidden]");
        }
        output.push('\n');
        output
    }

    fn generate_edge_style_template(
        &self,
        template_path: impl AsRef<std::path::Path>,
    ) -> Result<(), std::io::Error> {
        use std::fs;

        if let Some(parent) = template_path.as_ref().parent() {
            fs::create_dir_all(parent)?;
        }

        let mut content = String::from(
            r#"#import "crates/linnest/typst/src/physics-edge-style.typ": mi, massive, massless, dashed, dotted, stroke-style, source-stroke, sink-stroke, fermion-flow, wave, coil, zigzag, default-edge, style

// Auto-generated particle styles from the canonical model. The reusable
// drawing callbacks live in physics-edge-style.typ.
#let map = (
"#,
        );
        for particle in self.particles() {
            content.push_str(&format!(
                "  \"{}\": {},\n",
                particle.name,
                particle.generate_edge_typst_dict(self)
            ));
        }
        content.push_str(
            r#")

#let source-style(edge, typst-fields: "plain", ..options) = {
  let callbacks = style(map: map, typst-fields: typst-fields, ..options.named())
  (callbacks.source-style)(edge)
}

#let sink-style(edge, typst-fields: "plain", ..options) = {
  let callbacks = style(map: map, typst-fields: typst-fields, ..options.named())
  (callbacks.sink-style)(edge)
}

#let edge-label(edge, typst-fields: "plain", ..options) = {
  let callbacks = style(map: map, typst-fields: typst-fields, ..options.named())
  (callbacks.edge-label)(edge)
}
"#,
        );
        fs::write(template_path, content)
    }

    fn is_empty(&self) -> bool {
        self.name() == "ModelNotLoaded" || self.particles().is_empty()
    }
}

/// GammaLoop runtime operations on a stable particle reference.
pub trait ParticleIdGammaLoopExt {
    fn resolve(self, model: &Model) -> &Particle;
    fn antiparticle(self, model: &Model) -> ParticleId;
    fn is_self_antiparticle(&self, model: &Model) -> bool;
    fn is_massless(&self, model: &Model) -> bool;
}

impl ParticleIdGammaLoopExt for ParticleId {
    fn resolve(self, model: &Model) -> &Particle {
        model
            .particle_by_id(self)
            .expect("particle IDs stored with a model must resolve in that model")
    }

    fn antiparticle(self, model: &Model) -> ParticleId {
        self.resolve(model).antiparticle
    }

    fn is_self_antiparticle(&self, model: &Model) -> bool {
        model.particle_is_self_conjugate(*self)
    }

    fn is_massless(&self, model: &Model) -> bool {
        model.particle_is_massless(*self)
    }
}

/// GammaLoop runtime operations on a stable vertex-rule reference.
pub trait VertexRuleIdGammaLoopExt {
    fn resolve(self, model: &Model) -> &VertexRule;
}

impl VertexRuleIdGammaLoopExt for VertexRuleId {
    fn resolve(self, model: &Model) -> &VertexRule {
        model
            .vertex_rule_by_id(self)
            .expect("vertex-rule IDs stored with a model must resolve in that model")
    }
}

/// GammaLoop runtime operations on a stable propagator reference.
pub trait PropagatorIdGammaLoopExt {
    fn resolve(self, model: &Model) -> &Propagator;
}

impl PropagatorIdGammaLoopExt for PropagatorId {
    fn resolve(self, model: &Model) -> &Propagator {
        model
            .propagator_by_id(self)
            .expect("propagator IDs stored with a model must resolve in that model")
    }
}

/// Symbolica and tensor operations derived from a canonical particle record.
pub trait ParticleGammaLoopExt {
    fn symbolic_mass(&self, model: &Model) -> Atom;
    fn is_massless(&self, model: &Model) -> bool;
    fn resolved_mass_value(&self, model: &Model) -> Result<Complex<F<f64>>, Report>;
    fn has_zero_resolved_mass(&self, model: &Model) -> Result<bool, Report>;
    fn random_helicity(&self, seed: u64) -> Helicity;
    fn spin_reps(&self) -> IndexLess<LibraryRep, Aind>;
    fn color_reps(&self, flow: Flow) -> IndexLess;
    fn polarization_sum(
        &self,
        model: &Model,
        edge: EdgeIndex,
        average: bool,
        gauge: VectorPolarizationSumGauge,
    ) -> Result<Option<Replacement>, Report>;
    fn generate_edge_typst_dict(&self, model: &Model) -> String;
}

impl ParticleGammaLoopExt for Particle {
    fn symbolic_mass(&self, model: &Model) -> Atom {
        Atom::from(UFOSymbol::from(
            &model.parameter_by_id(self.mass).unwrap().name,
        ))
    }

    fn is_massless(&self, model: &Model) -> bool {
        let mass = model.parameter_by_id(self.mass).unwrap();
        mass.name == "ZERO"
            || mass
                .value
                .is_some_and(|value| value.re == 0.0 && value.im == 0.0)
            || mass
                .expression
                .as_ref()
                .is_some_and(|expression| expression.is_zero())
    }

    fn resolved_mass_value(&self, model: &Model) -> Result<Complex<F<f64>>, Report> {
        model
            .parameter_by_id(self.mass)
            .ok()
            .and_then(|parameter| parameter.value)
            .map(|mass| Complex::new(F(mass.re), F(mass.im)))
            .ok_or_else(|| {
                eyre!(
                    "Particle '{}' (PDG {}) requires a resolved value for mass parameter '{}' in model '{}'.",
                    self.name,
                    self.pdg_code,
                    model.parameter_by_id(self.mass).unwrap().name,
                    model.name()
                )
            })
    }

    fn has_zero_resolved_mass(&self, model: &Model) -> Result<bool, Report> {
        let mass = self.resolved_mass_value(model)?;
        Ok(mass.re == F(0.0) && mass.im == F(0.0))
    }

    fn random_helicity(&self, seed: u64) -> Helicity {
        let mut rng = SmallRng::seed_from_u64(seed);
        if self.is_spinor() {
            if rng.random_bool(0.5) {
                Helicity::PLUS
            } else {
                Helicity::MINUS
            }
        } else if self.is_vector() {
            Helicity::try_from(rng.random_range(1..=1)).unwrap()
        } else {
            Helicity::ZERO
        }
    }

    fn spin_reps(&self) -> IndexLess<LibraryRep, Aind> {
        PermutedStructure::<IndexLess<LibraryRep, Aind>>::from_iter(match self.spin {
            -1..=1 => vec![],
            spin if spin > 0 && spin % 2 == 0 => {
                vec![Bispinor {}.new_rep(4).cast(); spin.div_euclid(2) as usize]
            }
            spin if spin > 0 => {
                vec![Minkowski {}.new_rep(4).cast(); spin.div_euclid(2) as usize]
            }
            _ => vec![],
        })
        .structure
    }

    fn color_reps(&self, flow: Flow) -> IndexLess {
        let reps = match (flow, self.color) {
            (Flow::Source, 3) | (Flow::Sink, -3) => {
                vec![ColorFundamental {}.new_rep(3).cast()]
            }
            (Flow::Source, -3) | (Flow::Sink, 3) => {
                vec![ColorFundamental {}.dual().new_rep(3).cast()]
            }
            (Flow::Source, 6) | (Flow::Sink, -6) => {
                vec![ColorSextet {}.new_rep(6).cast()]
            }
            (Flow::Source, -6) | (Flow::Sink, 6) => {
                vec![ColorSextet {}.dual().new_rep(6).cast()]
            }
            (_, 8) => vec![ColorAdjoint {}.new_rep(8).cast()],
            _ => vec![],
        };
        PermutedStructure::<IndexLess>::from_iter(reps).structure
    }

    fn polarization_sum(
        &self,
        model: &Model,
        edge: EdgeIndex,
        average: bool,
        gauge: VectorPolarizationSumGauge,
    ) -> Result<Option<Replacement>, Report> {
        let average_factor = if average {
            match self.spin {
                1 => Atom::one(),
                2 => Atom::num(1) / Atom::num(2),
                3 => {
                    let states = if self.is_massless(model) { 2 } else { 3 };
                    Atom::num(1) / Atom::num(states)
                }
                4 => Atom::num(1) / Atom::num(4),
                5 => {
                    let states = if self.is_massless(model) { 4 } else { 5 };
                    Atom::num(1) / Atom::num(states)
                }
                spin => {
                    return Err(eyre!(
                        "Polarization averaging for particle '{}' (PDG {}, spin {}) is not supported yet.",
                        self.name,
                        self.pdg_code,
                        spin
                    ));
                }
            }
        } else {
            Atom::one()
        };

        Ok(match self.spin {
            1 => None,
            2 => {
                let mu: Slot<Minkowski, Aind> = Minkowski {}.new_rep(4).slot(Aind::new_dummy());
                let mass_sign = if self.is_antiparticle() { -1 } else { 1 };
                let rhs = GS.emr_mom(edge, mu.to_atom())
                    * function!(AGS.gamma, W_.a_, W_.b_, mu.to_atom())
                    + Atom::num(mass_sign)
                        * self.symbolic_mass(model)
                        * function!(ETS.metric, W_.a_, W_.b_);
                Some(if !self.is_antiparticle() {
                    (function!(GS.u, edge.0, W_.a_) * function!(GS.ubar, edge.0, W_.b_))
                        .replace_with(rhs * average_factor)
                } else {
                    (function!(GS.v, edge.0, W_.a_) * function!(GS.vbar, edge.0, W_.b_))
                        .replace_with(rhs * average_factor)
                })
            }
            3 => {
                let minus_metric = -function!(ETS.metric, W_.a_, W_.b_);
                let rhs = match gauge {
                    VectorPolarizationSumGauge::Feynman => minus_metric,
                    VectorPolarizationSumGauge::LightLikeAxial if !self.is_massless(model) => {
                        minus_metric
                            + GS.emr_mom(edge, W_.a_) * GS.emr_mom(edge, W_.b_)
                                / self.symbolic_mass(model).pow(2)
                    }
                    VectorPolarizationSumGauge::LightLikeAxial => {
                        let temporal_component = GS.emr_mom(edge, GS.cind(0));
                        let n_a = temporal_component.clone() * GS.energy_delta(W_.a_)
                            - GS.emr_vec_index(edge, W_.a_);
                        let n_b = temporal_component.clone() * GS.energy_delta(W_.b_)
                            - GS.emr_vec_index(edge, W_.b_);
                        let q_dot_n = temporal_component.pow(2)
                            + Euclidean {}
                                .new_rep(4)
                                .inner_product(GS.emr_vec(edge), GS.emr_vec(edge));
                        minus_metric
                            + (GS.emr_mom(edge, W_.a_) * n_b + n_a * GS.emr_mom(edge, W_.b_))
                                / q_dot_n
                    }
                };
                Some(
                    (function!(GS.epsilon, edge.0, W_.a_)
                        * function!(GS.epsilonbar, edge.0, W_.b_))
                    .replace_with(rhs * average_factor),
                )
            }
            spin => {
                return Err(eyre!(
                    "Polarization sum replacement for particle '{}' (PDG {}, spin {}) is not implemented yet.",
                    self.name,
                    self.pdg_code,
                    spin
                ));
            }
        })
    }

    fn generate_edge_typst_dict(&self, model: &Model) -> String {
        let label = format!("mi(`{}`)", self.texname);
        let thickness = if self.is_massless(model) {
            "massless"
        } else {
            "massive"
        };
        let color = if self.charge.abs() > 0.0 {
            "blue"
        } else {
            "black"
        };
        let source = format!("source-stroke(c: {color}, thickness: {thickness})");
        let sink = format!("sink-stroke(c: {color}, thickness: {thickness})");
        let (source, sink) = if self.is_ghost() {
            (
                format!("source-stroke(c: {color}, thickness: {thickness}, dash: dotted)"),
                format!("sink-stroke(c: {color}, thickness: {thickness}, dash: dotted)"),
            )
        } else if self.is_vector() && self.charge == 0.0 && self.color == 1 {
            (format!("{source} + wave"), format!("{sink} + wave"))
        } else if self.is_vector() && self.charge == 0.0 && self.color == 8 {
            (format!("{source} + coil"), format!("{sink} + coil"))
        } else if self.is_vector() {
            (format!("{source} + zigzag"), format!("{sink} + zigzag"))
        } else if self.is_scalar() {
            (
                format!("source-stroke(c: {color}, thickness: {thickness}, dash: dashed)"),
                format!("sink-stroke(c: {color}, thickness: {thickness}, dash: dashed)"),
            )
        } else {
            (source, sink)
        };
        let flow_marker = if self.is_fermion() && !self.is_ghost() {
            " + fermion-flow"
        } else {
            ""
        };
        format!("(source:{source}, sink:{sink}, label:{label}){flow_marker}")
    }
}

/// Parsed tensor data derived from a canonical UFO vertex rule.
pub trait VertexRuleGammaLoopExt {
    fn tensors(
        &self,
        model: &Model,
        i: Aind,
        j: Aind,
    ) -> [ParamTensor<OrderedStructure<Euclidean, Aind>>; 3];
    fn degree_of_divergence(&self, model: &Model) -> i32;
}

impl VertexRuleGammaLoopExt for VertexRule {
    fn tensors(
        &self,
        model: &Model,
        i: Aind,
        j: Aind,
    ) -> [ParamTensor<OrderedStructure<Euclidean, Aind>>; 3] {
        let spin_structure = self
            .lorentz_structures
            .iter()
            .map(|id| {
                model
                    .lorentz_structure_by_id(*id)
                    .unwrap()
                    .structure
                    .clone()
            })
            .collect_vec();
        let color_structure = self.color_structures.clone();

        let color_slot = Euclidean {}.new_slot(color_structure.len(), i);
        let spin_slot = Euclidean {}.new_slot(spin_structure.len(), j);
        let color_structure = ParamTensor::composite(DataTensor::Dense(
            DenseTensor::from_data(
                color_structure,
                PermutedStructure::from_iter([color_slot]).structure,
            )
            .unwrap(),
        ));
        let spin_structure = ParamTensor::composite(DataTensor::Dense(
            DenseTensor::from_data(
                spin_structure,
                PermutedStructure::from_iter([spin_slot]).structure,
            )
            .unwrap(),
        ));
        let mut couplings = ParamTensor::composite(DataTensor::Sparse(SparseTensor::empty(
            PermutedStructure::from_iter([color_slot, spin_slot]).structure,
            Atom::Zero,
        )));
        for (color, row) in self.couplings.iter().enumerate() {
            for (spin, coupling) in row.iter().enumerate() {
                if let Some(coupling) = coupling {
                    couplings
                        .set(
                            &[color, spin],
                            Atom::from(UFOSymbol::from(
                                &model.coupling_by_id(*coupling).unwrap().name,
                            )),
                        )
                        .unwrap();
                }
            }
        }
        [color_structure, couplings, spin_structure]
    }

    fn degree_of_divergence(&self, model: &Model) -> i32 {
        self.lorentz_structures
            .iter()
            .map(|id| {
                model
                    .lorentz_structure_by_id(*id)
                    .unwrap()
                    .structure
                    .all_dod()
            })
            .max()
            .unwrap_or_default()
    }
}

/// Parsed expressions derived from a canonical UFO propagator.
pub trait PropagatorGammaLoopExt {
    fn numerator_atom(&self) -> Atom;
    fn denominator_atom(&self) -> Atom;
    fn degree_of_divergence(&self) -> i32;
}

impl PropagatorGammaLoopExt for Propagator {
    fn numerator_atom(&self) -> Atom {
        self.numerator.clone()
    }

    fn denominator_atom(&self) -> Atom {
        self.denominator.clone()
    }

    fn degree_of_divergence(&self) -> i32 {
        (self.numerator_atom() / self.denominator_atom()).all_dod()
    }
}
