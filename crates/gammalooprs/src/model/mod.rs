use crate::HasModel;
use crate::momentum::Helicity;
use crate::numerator::aind::Aind;
use crate::utils::serde_utils::SmartSerde;
use crate::utils::symbolica_ext::{CallSymbol, Replaces};
use crate::utils::{self, F, W_};
use ahash::{AHashMap, HashSet, RandomState};
use bincode::{Decode, Encode};
use color_eyre::Report;
use color_eyre::owo_colors::OwoColorize;
use colored::Colorize;
use eyre::eyre;
use itertools::Itertools;
// use linnet::half_edge::drawing::Decoration;
use linnet::half_edge::involution::{EdgeIndex, Flow};
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use serde::de::DeserializeOwned;
use spenso::shadowing::symbolica_utils::SpensoPrintSettings;
use spenso::structure::{IndexLess, PermutedStructure};
use tabled::settings::Modify;
use tabled::{
    builder::Builder,
    settings::{Span, object::Cell},
    settings::{Style, style::VerticalLine},
};

// use log::{info, trace};
use idenso::representations::{Bispinor, ColorAdjoint, ColorFundamental, ColorSextet};
use serde::{Deserialize, Serialize};
use smartstring::{LazyCompact, SmartString};
use spenso::algebra::complex::Complex;
use spenso::network::library::symbolic::ETS;
use spenso::structure::OrderedStructure;
use spenso::structure::representation::Euclidean;
use spenso::structure::representation::{LibraryRep, Minkowski};
use spenso::structure::{
    representation::BaseRepName, representation::Lorentz, representation::RepName,
    slot::IsAbstractSlot, slot::Slot,
};
use spenso::tensors::data::{DataTensor, DenseTensor, SetTensorData, SparseTensor};
use spenso::tensors::parametric::ParamTensor;
use std::collections::BTreeMap;
use std::fmt::{Display, Formatter};
use std::fs;
use symbolica::domains::float::Real;
use symbolica::domains::integer::IntegerRing;
use symbolica::domains::rational::{Fraction, Rational};

type ExternalFunctionMap = HashMap<String, Box<dyn ExternalFunction<Complex<F<f64>>>>, RandomState>;
use symbolica::evaluate::{ExternalFunction, FunctionMap, OptimizationSettings};
use symbolica::id::Replacement;
use tracing::info;

use color_eyre::Result;
use std::collections::HashMap;
use std::ops::Deref;
use std::path::Path;
use std::sync::Arc;
use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};

use crate::utils::GS;

use symbolica::printer::{AtomPrinter, PrintOptions};
use symbolica::{function, get_symbol, parse, symbol};

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SerializableInputParamCard<T> {
    pub data: BTreeMap<String, (T, T)>,
}

impl<T> SerializableInputParamCard<T>
where
    T: From<f64> + Clone + Serialize + DeserializeOwned,
{
    pub(crate) fn from_file(file_path: impl AsRef<Path>) -> Result<Self, Report> {
        let hashmap_card: BTreeMap<String, (T, T)> = SmartSerde::from_file(file_path, "model")?;
        Ok(SerializableInputParamCard { data: hashmap_card })
    }

    pub fn from_str(s: String, format: &str) -> Result<Self, Report> {
        let hashmap_card: BTreeMap<String, (T, T)> =
            SmartSerde::from_str(s, format, "model_parameters")?;
        Ok(SerializableInputParamCard { data: hashmap_card })
    }

    pub fn from_input_param_card(card: &InputParamCard<T>) -> Self {
        let serializeable_card: BTreeMap<String, (T, T)> = card
            .data
            .iter()
            .map(|(k, v)| {
                (
                    k.namespaceless_string().to_string(),
                    (v.re.clone(), v.im.clone()),
                )
            })
            .collect();
        SerializableInputParamCard {
            data: serializeable_card,
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
    /// Create empty card
    pub fn new() -> Self {
        Self {
            data: HashMap::new(),
        }
    }

    /// Insert a parameter
    pub fn insert(&mut self, sym: UFOSymbol, value: Complex<T>) -> Option<Complex<T>> {
        self.data.insert(sym, value)
    }

    /// Get immutable reference
    pub fn get(&self, sym: &UFOSymbol) -> Option<&Complex<T>> {
        self.data.get(sym)
    }

    /// Get mutable reference
    pub fn get_mut(&mut self, sym: &UFOSymbol) -> Option<&mut Complex<T>> {
        self.data.get_mut(sym)
    }

    /// Remove a parameter
    pub fn remove(&mut self, sym: &UFOSymbol) -> Option<Complex<T>> {
        self.data.remove(sym)
    }

    /// Iterate over all parameters
    pub fn iter(&self) -> impl Iterator<Item = (&UFOSymbol, &Complex<T>)> {
        self.data.iter()
    }

    pub fn to_serializable(&self) -> SerializableInputParamCard<T> {
        SerializableInputParamCard::from_input_param_card(self)
    }

    pub fn from_serializable(
        serializable_input_param_card: &SerializableInputParamCard<T>,
    ) -> Self {
        let data: HashMap<UFOSymbol, Complex<T>> = serializable_input_param_card
            .data
            .iter()
            .map(|(k, v)| {
                (
                    UFOSymbol::from(k.as_str()),
                    Complex::new(v.0.clone(), v.1.clone()),
                )
            })
            .collect();
        InputParamCard { data }
    }

    pub fn from_file(file_path: impl AsRef<Path>) -> Result<Self, Report> {
        let serializable_input_param_card = SerializableInputParamCard::from_file(file_path)?;
        Ok(Self::from_serializable(&serializable_input_param_card))
    }
    pub fn from_str(s: String, format: &str) -> Result<Self, Report> {
        let serializable_input_param_card = SerializableInputParamCard::from_str(s, format)?;
        Ok(Self::from_serializable(&serializable_input_param_card))
    }

    pub fn to_file<P: AsRef<Path>>(&self, path: P, overwrite: bool) -> Result<(), Report> {
        let serializable_card = self.to_serializable();
        serializable_card.to_file(path, overwrite)?;
        Ok(())
    }
}

impl InputParamCard<F<f64>> {
    pub fn apply_to_model(&self, model: &mut Model) -> Result<(), Report> {
        for (param, value) in &self.data {
            if let Some(model_param) = model.get_parameter_mut_opt(param.namespaceless_string()) {
                model_param.value = Some(*value);
            } else {
                return Err(eyre!(
                    "Parameter {} not found in model when applying input parameter card",
                    param
                ));
            }
        }
        model.recompute_dependents()
    }
}

impl InputParamCard<F<f64>> {
    pub fn default_from_model(model: &Model) -> Self {
        let mut card = InputParamCard::new();
        for param in model.parameters.values() {
            if param.nature == ParameterNature::External
                && !param.name.is_zero()
                && let Some(value) = param.value
            {
                card.insert(param.name, value);
            }
        }
        card
    }
}

/// Let it behave like a map if desired
impl<T> std::ops::Deref for InputParamCard<T>
where
    T: From<f64> + Clone + Serialize + DeserializeOwned,
{
    type Target = HashMap<UFOSymbol, Complex<T>>;
    fn deref(&self) -> &Self::Target {
        &self.data
    }
}
impl<T> std::ops::DerefMut for InputParamCard<T>
where
    T: From<f64> + Clone + Serialize + DeserializeOwned,
{
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct UFOSymbol(pub Symbol);

impl Display for UFOSymbol {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl From<UFOSymbol> for Atom {
    fn from(value: UFOSymbol) -> Self {
        Atom::var(value.0)
    }
}

impl UFOSymbol {
    pub fn zero() -> Self {
        UFOSymbol(GS.ufozero)
    }

    pub fn is_zero(&self) -> bool {
        self.0 == Self::zero().0
    }

    pub fn namespaceless_string(&self) -> &str {
        self.0.get_stripped_name()
    }
}

#[test]
fn zerosym() {
    println!("{}", Atom::num(1) * Atom::from(UFOSymbol::zero()));
    println!("{}", Atom::from(UFOSymbol::zero()));
    // assert_eq!(UFOSymbol::zero().namespaceless_string(), "ZERO");
}

impl<T> From<T> for UFOSymbol
where
    T: AsRef<str>,
{
    fn from(s: T) -> Self {
        let is_zero = s.as_ref() == "ZERO";
        if is_zero {
            UFOSymbol::zero()
        } else {
            let name = format!("UFO::{}", s.as_ref());
            if let Some(a) = get_symbol!(&name) {
                UFOSymbol(a)
            } else {
                UFOSymbol(symbol!(
                    &name,
                    print = |a, opt| {
                        let AtomView::Var(a) = a else {
                            return None;
                        };
                        match opt.custom_print_mode {
                            Some(("spenso", i)) => {
                                let SpensoPrintSettings { .. } = SpensoPrintSettings::from(i);
                                if SpensoPrintSettings::from(i).is_typst() {
                                    Some(format!("\"{}\"", a.get_symbol().get_stripped_name()))
                                } else {
                                    None
                                }
                            }
                            _ => None,
                        }
                    }
                ))
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct ArcPropagator(pub Arc<Propagator>);
impl Deref for ArcPropagator {
    type Target = Propagator;
    fn deref(&self) -> &Self::Target {
        self.0.deref()
    }
}

impl Encode for ArcPropagator {
    fn encode<E: bincode::enc::Encoder>(
        &self,
        encoder: &mut E,
    ) -> std::result::Result<(), bincode::error::EncodeError> {
        Encode::encode(&self.0.name.to_string(), encoder)?;
        Ok(())
    }
}

impl<T: HasModel> Decode<T> for ArcPropagator {
    fn decode<D: bincode::de::Decoder<Context = T>>(
        decoder: &mut D,
    ) -> std::result::Result<Self, bincode::error::DecodeError> {
        let name: String = Decode::decode(decoder)?;
        let context = decoder.context();
        let model = context.get_model();
        let prop = model.get_propagator(&name);
        Ok(prop)
    }
}

#[derive(Debug, Clone, Eq, Ord, PartialEq, PartialOrd, Hash)]
pub struct ArcParticle(pub Arc<Particle>);

impl Deref for ArcParticle {
    type Target = Particle;
    fn deref(&self) -> &Self::Target {
        self.0.deref()
    }
}

#[derive(Debug, Clone, Eq, Ord, PartialEq, PartialOrd, Hash)]
pub struct ArcVertexRule(pub Arc<VertexRule>);
impl Deref for ArcVertexRule {
    type Target = VertexRule;
    fn deref(&self) -> &Self::Target {
        self.0.deref()
    }
}

impl Encode for ArcParticle {
    fn encode<E: bincode::enc::Encoder>(
        &self,
        encoder: &mut E,
    ) -> std::result::Result<(), bincode::error::EncodeError> {
        Encode::encode(&self.0.pdg_code, encoder)?;
        Ok(())
    }
}

impl<T: HasModel> Decode<T> for ArcParticle {
    fn decode<D: bincode::de::Decoder<Context = T>>(
        decoder: &mut D,
    ) -> std::result::Result<Self, bincode::error::DecodeError> {
        let pdg_code: isize = Decode::decode(decoder)?;
        let context = decoder.context();
        let model = context.get_model();
        let particle = model.get_particle_from_pdg(pdg_code);
        Ok(particle)
    }
}

impl Encode for ArcVertexRule {
    fn encode<E: bincode::enc::Encoder>(
        &self,
        encoder: &mut E,
    ) -> std::result::Result<(), bincode::error::EncodeError> {
        Encode::encode(&self.0.name.to_string(), encoder)?;
        Ok(())
    }
}

impl<T: HasModel> Decode<T> for ArcVertexRule {
    fn decode<D: bincode::de::Decoder<Context = T>>(
        decoder: &mut D,
    ) -> std::result::Result<Self, bincode::error::DecodeError> {
        let vertex_rule_name: String = Decode::decode(decoder)?;
        let context = decoder.context();
        let model = context.get_model();
        let vertex_rule = model.get_vertex_rule(vertex_rule_name);
        Ok(vertex_rule)
    }
}

#[allow(unused)]
pub(crate) fn normalise_complex(atom: &Atom) -> Atom {
    let re = parse!("re_");
    let im = parse!("im_");

    let comp_id = symbol!("complex");

    let complexfn = function!(comp_id, re, im).to_pattern();

    let i = Atom::i();
    let complexpanded = &re + i * &im;

    atom.replace(&complexfn).with(complexpanded.to_pattern())
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq, Eq)]
pub enum ParameterNature {
    #[default]
    #[serde(rename = "external")]
    External,
    #[serde(rename = "internal")]
    Internal,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq, Eq)]
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
    pub(crate) fn from_vertex_rule(vertex_rule: &VertexRule) -> SerializableVertexRule {
        SerializableVertexRule {
            name: vertex_rule.name.clone(),
            particles: vertex_rule
                .particles
                .iter()
                .map(|particle| particle.0.name.clone())
                .collect(),
            color_structures: vertex_rule
                .color_structures
                .iter()
                .map(|a| a.to_canonical_string())
                .map(SmartString::from)
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
                        .map(|coupling| {
                            coupling
                                .as_ref()
                                .map(|cpl| cpl.namespaceless_string().into())
                        })
                        .collect()
                })
                .collect(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ColorStructure {
    pub color_structure: Vec<Atom>,
}

impl ColorStructure {
    pub(crate) fn iter(&'_ self) -> std::slice::Iter<'_, Atom> {
        self.color_structure.iter()
    }
}

impl FromIterator<Atom> for ColorStructure {
    fn from_iter<T: IntoIterator<Item = Atom>>(iter: T) -> Self {
        ColorStructure {
            color_structure: iter.into_iter().collect(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct VertexRule {
    pub name: SmartString<LazyCompact>,
    pub particles: Vec<ArcParticle>,
    pub color_structures: ColorStructure,
    pub lorentz_structures: Vec<Arc<LorentzStructure>>,
    pub couplings: Vec<Vec<Option<CouplingName>>>,
}

impl Eq for VertexRule {}

impl PartialEq for VertexRule {
    fn eq(&self, other: &Self) -> bool {
        self.name == other.name
    }
}

impl PartialOrd for VertexRule {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for VertexRule {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.name.cmp(&other.name)
    }
}

impl std::hash::Hash for VertexRule {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.name.hash(state);
    }
}

impl VertexRule {
    pub(crate) fn tensors(
        &self,
        i: Aind,
        j: Aind,
    ) -> [ParamTensor<OrderedStructure<Euclidean, Aind>>; 3] {
        let spin_structure = self
            .lorentz_structures
            .iter()
            .map(|ls| ls.structure.clone())
            .collect_vec();

        let color_structure: Vec<Atom> = self.color_structures.color_structure.clone();

        let i = Euclidean {}.new_slot(color_structure.len(), i);
        let j = Euclidean {}.new_slot(spin_structure.len(), j);

        let color_structure: ParamTensor<OrderedStructure<Euclidean, Aind>> =
            ParamTensor::composite(DataTensor::Dense(
                DenseTensor::from_data(
                    color_structure,
                    PermutedStructure::from_iter([i]).structure,
                )
                .unwrap(),
            ));

        let spin_structure: ParamTensor<OrderedStructure<Euclidean, Aind>> =
            ParamTensor::composite(DataTensor::Dense(
                DenseTensor::from_data(spin_structure, PermutedStructure::from_iter([j]).structure)
                    .unwrap(),
            ));

        let mut couplings: ParamTensor<OrderedStructure<Euclidean, Aind>> =
            ParamTensor::composite(DataTensor::Sparse(SparseTensor::empty(
                PermutedStructure::from_iter([i, j]).structure,
                Atom::Zero,
            )));

        for (i, row) in self.couplings.iter().enumerate() {
            for (j, col) in row.iter().enumerate() {
                if let Some(atom) = col {
                    couplings.set(&[i, j], atom.0.into()).unwrap();
                }
            }
        }

        [color_structure, couplings, spin_structure]
    }

    // #[allow(clippy::complexity)]
    // pub fn get_coupling_orders(&self) -> Vec<Vec<Option<CouplingName>>> {
    //     self.couplings
    //         .iter()
    //         .map(|row| {
    //             row.iter()
    //                 .map(|co| co.clone().map(|c| c.clone()))
    //                 .collect::<Vec<_>>()
    //         })
    //         .collect::<Vec<_>>()
    // }
}

impl VertexRule {
    pub(crate) fn coupling_orders(
        &self,
        model: &Model,
    ) -> AHashMap<SmartString<LazyCompact>, usize> {
        let mut node_coupling_orders = AHashMap::default();
        self.couplings.iter().for_each(|cs| {
            cs.iter().for_each(|c_opt| {
                if let Some(c) = c_opt {
                    model.couplings[c]
                        .orders
                        .iter()
                        .for_each(|(coupling_order, &weight)| {
                            let w = node_coupling_orders
                                .entry(coupling_order.clone())
                                .or_insert(weight);
                            if *w < weight {
                                *w = weight;
                            }
                        });
                }
            })
        });
        node_coupling_orders
    }

    pub(crate) fn from_serializable_vertex_rule(
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
                                .map(|cpl_name| CouplingName(UFOSymbol::from(cpl_name)))
                        })
                        .collect()
                })
                .collect(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerializablePropagator {
    pub name: SmartString<LazyCompact>,
    pub particle: SmartString<LazyCompact>,
    pub numerator: SmartString<LazyCompact>,
    pub denominator: SmartString<LazyCompact>,
}

impl SerializablePropagator {
    pub(crate) fn from_propagator(propagator: &Propagator) -> SerializablePropagator {
        SerializablePropagator {
            name: propagator.name.clone(),
            particle: propagator.particle.0.name.clone(),
            numerator: propagator.numerator.to_canonical_string().into(),
            denominator: propagator.denominator.to_canonical_string().into(),
        }
    }
}

#[derive(Debug, Clone, Encode)]
pub struct Propagator {
    #[bincode(with_serde)]
    pub name: SmartString<LazyCompact>,
    pub particle: ArcParticle,
    pub numerator: Atom,
    pub denominator: Atom,
}

impl Propagator {
    pub(crate) fn from_serializable_propagator(
        model: &Model,
        propagator: &SerializablePropagator,
    ) -> Propagator {
        Propagator {
            name: propagator.name.clone(),
            particle: model.get_particle(&propagator.particle).clone(),
            numerator: utils::parse_python_expression(propagator.numerator.as_str()),
            denominator: utils::parse_python_expression(propagator.denominator.as_str()),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerializableCoupling {
    name: SmartString<LazyCompact>,
    expression: SmartString<LazyCompact>,
    #[serde(with = "vectorize")]
    orders: BTreeMap<SmartString<LazyCompact>, usize>,
    value: Option<(f64, f64)>,
}

impl SerializableCoupling {
    pub(crate) fn from_coupling(coupling: &Coupling) -> SerializableCoupling {
        SerializableCoupling {
            name: coupling.name.namespaceless_string().into(),
            expression: coupling.expression.to_canonical_string().into(),
            orders: coupling.orders.clone(),
            value: coupling.value.map(|value| (value.re, value.im)),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CouplingName(pub UFOSymbol);

impl Deref for CouplingName {
    type Target = UFOSymbol;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[derive(Debug, Clone)]
pub struct Coupling {
    pub name: UFOSymbol,
    pub expression: Atom,
    pub orders: BTreeMap<SmartString<LazyCompact>, usize>,
    pub value: Option<Complex<f64>>,
}

impl Coupling {
    pub(crate) fn from_serializable_coupling(coupling: &SerializableCoupling) -> Coupling {
        Coupling {
            name: (&coupling.name).into(),
            expression: utils::parse_python_expression(coupling.expression.as_str()),
            orders: coupling.orders.clone(),
            value: coupling.value.map(|value| Complex::new(value.0, value.1)),
        }
    }

    pub(crate) fn rep_rule(&self) -> [Atom; 2] {
        let lhs = self.name.into();
        //let rhs = normalise_complex(&self.expression);
        let rhs = self.expression.clone();

        [lhs, rhs]
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
    pub(crate) fn from_particle(particle: &Particle) -> SerializableParticle {
        SerializableParticle {
            pdg_code: particle.pdg_code,
            name: particle.name.clone(),
            antiname: particle.antiname.clone(),
            spin: particle.spin,
            color: particle.color,
            mass: particle.mass.namespaceless_string().into(),
            width: particle.width.namespaceless_string().into(),
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
    pub mass: ParameterName,
    pub width: ParameterName,
    pub texname: SmartString<LazyCompact>,
    pub antitexname: SmartString<LazyCompact>,
    pub charge: f64,
    pub ghost_number: isize,
    pub lepton_number: isize,
    pub y_charge: isize,
}

impl Particle {
    pub(crate) fn spin_sum(&self, eid: EdgeIndex) -> Option<Replacement> {
        match self.spin {
            2 => {
                if !self.is_antiparticle() {
                    Some(
                        (function!(GS.u, eid.0, W_.a_) * function!(GS.ubar, eid.0, W_.b_))
                            .replace_with(function!(ETS.metric, W_.a_, W_.b_)),
                    )
                } else {
                    Some(
                        (function!(GS.v, eid.0, W_.a_) * function!(GS.vbar, eid.0, W_.b_))
                            .replace_with(function!(ETS.metric, W_.a_, W_.b_)),
                    )
                }
            }
            3 => Some(
                (function!(GS.epsilon, eid.0, W_.a_) * function!(GS.epsilonbar, eid.0, W_.b_))
                    .replace_with(function!(ETS.metric, W_.a_, W_.b_)),
            ),
            _ => None,
        }
    }
    pub(crate) fn random_helicity(&self, seed: u64) -> Helicity {
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

    pub(crate) fn is_massless(&self) -> bool {
        !self.is_massive()
    }
    // pub fn decoration(&self) -> Decoration {
    //     match self.spin {
    //         0 => Decoration::Dashed,
    //         1 => Decoration::None,
    //         2 => Decoration::Arrow,
    //         3 => {
    //             if self.pdg_code.abs() == 9 || self.pdg_code.abs() == 21 {
    //                 Decoration::Coil
    //             } else {
    //                 Decoration::Wave
    //             }
    //         }
    //         _ => Decoration::None,
    //     }
    // }

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

    pub fn symbolic_mass(&self) -> Atom {
        self.mass.0.into()
    }
}

impl Ord for Particle {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.pdg_code.cmp(&other.pdg_code)
    }
}

impl PartialOrd for Particle {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl std::hash::Hash for Particle {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.pdg_code.hash(state);
    }
}
impl Eq for Particle {}
impl PartialEq for Particle {
    fn eq(&self, other: &Self) -> bool {
        if self.pdg_code == other.pdg_code {
            // The checks below are too slow
            // if self.name != other.name {
            //     panic!(
            //         "Particle with same pdg code but different names: {} and {}",
            //         self.name, other.name
            //     );
            // }
            // if self.spin != other.spin {
            //     panic!(
            //         "Particle with same pdg code but different spins: {} and {}",
            //         self.spin, other.spin
            //     );
            // }
            // if self.color != other.color {
            //     panic!(
            //         "Particle with same pdg code but different colors: {} and {}",
            //         self.color, other.color
            //     );
            // }
            // if self.mass != other.mass {
            //     panic!(
            //         "Particle with same pdg code but different masses: {} and {}",
            //         self.mass, other.mass
            //     );
            // }
            // if self.width != other.width {
            //     panic!(
            //         "Particle with same pdg code but different widths: {} and {}",
            //         self.width, other.width
            //     );
            // }
            // if self.texname != other.texname {
            //     panic!(
            //         "Particle with same pdg code but different texnames: {} and {}",
            //         self.texname, other.texname
            //     );
            // }
            // if self.antitexname != other.antitexname {
            //     panic!(
            //         "Particle with same pdg code but different antitexnames: {} and {}",
            //         self.antitexname, other.antitexname
            //     );
            // }
            // if self.charge != other.charge {
            //     panic!(
            //         "Particle with same pdg code but different charges: {} and {}",
            //         self.charge, other.charge
            //     );
            // }
            // if self.ghost_number != other.ghost_number {
            //     panic!(
            //         "Particle with same pdg code but different ghost_numbers: {} and {}",
            //         self.ghost_number, other.ghost_number
            //     );
            // }
            // if self.lepton_number != other.lepton_number {
            //     panic!(
            //         "Particle with same pdg code but different lepton_numbers: {} and {}",
            //         self.lepton_number, other.lepton_number
            //     );
            // }
            // if self.y_charge != other.y_charge {
            //     panic!(
            //         "Particle with same pdg code but different y_charges: {} and {}",
            //         self.y_charge, other.y_charge
            //     );
            // }
            true
        } else {
            false
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InOutIndex {
    incoming: Slot<LibraryRep>,
    outgoing: Slot<LibraryRep>,
}

#[derive(Debug, Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = symbolica::state::HasStateMap)]
pub struct EdgeSlots<LorRep: RepName> {
    pub lorentz: Vec<Slot<LorRep>>,
    spin: Vec<Slot<Bispinor>>,
    pub color: Vec<Slot<LibraryRep>>,
}

impl From<EdgeSlots<Minkowski>> for OrderedStructure {
    fn from(value: EdgeSlots<Minkowski>) -> Self {
        PermutedStructure::<OrderedStructure>::from(
            value
                .lorentz
                .into_iter()
                .map(|x| x.to_lib())
                .chain(value.spin.into_iter().map(|x| x.to_lib()))
                .chain(value.color)
                .collect_vec(),
        )
        .structure
    }
}

impl<LorRep: BaseRepName> Display for EdgeSlots<LorRep> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "Lorentz: ")?;
        for l in &self.lorentz {
            write!(f, "{} ", l)?;
        }
        write!(f, "Spin: ")?;
        for s in &self.spin {
            write!(f, "{} ", s)?;
        }
        write!(f, "Color: ")?;
        for c in &self.color {
            write!(f, "{} ", c)?;
        }
        Ok(())
    }
}

impl From<EdgeSlots<Lorentz>> for OrderedStructure {
    fn from(value: EdgeSlots<Lorentz>) -> Self {
        PermutedStructure::<OrderedStructure>::from(
            value
                .lorentz
                .into_iter()
                .map(|x| x.to_lib())
                .chain(value.spin.into_iter().map(|a| a.to_lib()))
                .chain(value.color)
                .collect_vec(),
        )
        .structure
    }
}

impl Particle {
    pub fn is_antiparticle(&self) -> bool {
        self.pdg_code < 0
    }

    pub(crate) fn get_anti_particle(&self, model: &Model) -> ArcParticle {
        model.get_particle(&self.antiname)
    }

    pub fn is_self_antiparticle(&self) -> bool {
        self.name == self.antiname
    }

    pub(crate) fn spin_reps(&self) -> IndexLess<LibraryRep, Aind> {
        PermutedStructure::<IndexLess<LibraryRep, Aind>>::from_iter(match self.spin {
            -1..=1 => vec![],
            a => {
                if a > 0 {
                    if a % 2 == 0 {
                        vec![Bispinor {}.new_rep(4).cast(); a.div_euclid(2) as usize]
                    } else {
                        vec![Minkowski {}.new_rep(4).cast(); a.div_euclid(2) as usize]
                    }
                } else {
                    vec![]
                }
            }
        })
        .structure
    }

    pub(crate) fn is_massive(&self) -> bool {
        !self.mass.is_zero()
    }

    /// Generate edge styles for visualization based on particle properties
    pub fn generate_edge_typst_dict(&self) -> String {
        let label = format!("mi(`{}`)", self.texname);
        // Determine line thickness based on mass
        let thickness = if self.is_massive() {
            "massive"
        } else {
            "massless"
        };

        // Determine color based on charge
        let color = if self.charge.abs() > 0.0 {
            "blue"
        } else {
            "black"
        };

        // Generate base styles
        let base_source = format!("source_stroke(c: {}, thickness: {})", color, thickness);
        let base_sink = format!("sink_stroke(c: {}, thickness: {})", color, thickness);

        let (source, sink) = if self.is_ghost() {
            (
                format!("source_stroke(c: {color}, thickness: {thickness},dash: dotted)",),
                format!("sink_stroke(c: {color}, thickness: {thickness},dash: dotted)"),
            )
        } else if self.is_fermion() {
            (base_source, base_sink)
        } else if self.is_vector() {
            // Vector bosons: differentiate based on charge and color properties
            if self.charge == 0.0 && self.color == 1 {
                // Neutral color singlet (photon): wavy line
                (
                    format!("{} + wave", base_source),
                    format!("{} + wave", base_sink),
                )
            } else if self.charge == 0.0 && self.color == 8 {
                // Neutral color octet (gluon): coiled line
                (
                    format!("{} + coil", base_source),
                    format!("{} + coil", base_sink),
                )
            } else {
                // Charged vector bosons (W+/W-/Z with mass): zigzag line
                (
                    format!("{} + zigzag", base_source),
                    format!("{} + zigzag", base_sink),
                )
            }
        } else if self.is_scalar() {
            // Scalar particles: dashed lines
            (
                format!("source_stroke(c: {color}, thickness: {thickness},dash: dashed)",),
                format!("sink_stroke(c: {color}, thickness: {thickness},dash: dashed)"),
            )
        } else {
            // Default: solid line
            (base_source, base_sink)
        };

        format!("(source:{}, sink:{}, label:{})", source, sink, label)
    }

    pub(crate) fn color_reps(&self, flow: Flow) -> IndexLess {
        let reps = match flow {
            Flow::Source => match self.color {
                3 => vec![ColorFundamental {}.new_rep(3).cast()],

                -3 => vec![ColorFundamental {}.dual().new_rep(3).cast()],
                6 => vec![ColorSextet {}.new_rep(6).cast()],
                -6 => vec![ColorSextet {}.dual().new_rep(6).cast()],
                8 => vec![ColorAdjoint {}.new_rep(8).cast()],
                _ => vec![],
            },
            Flow::Sink => match self.color {
                -3 => vec![ColorFundamental {}.new_rep(3).cast()],
                3 => vec![ColorFundamental {}.dual().new_rep(3).cast()],
                -6 => vec![ColorSextet {}.new_rep(6).cast()],
                6 => vec![ColorSextet {}.dual().new_rep(6).cast()],
                8 => vec![ColorAdjoint {}.new_rep(8).cast()],
                _ => vec![],
            },
        };
        PermutedStructure::<IndexLess>::from_iter(reps).structure
    }

    pub(crate) fn from_serializable_particle(particle: &SerializableParticle) -> Particle {
        Particle {
            pdg_code: particle.pdg_code,
            name: particle.name.clone(),
            antiname: particle.antiname.clone(),
            spin: particle.spin,
            color: particle.color,
            mass: ParameterName((&particle.mass).into()),
            width: ParameterName((&particle.width).into()),
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
    pub(crate) fn from_lorentz_structure(ls: &LorentzStructure) -> SerializableLorentzStructure {
        SerializableLorentzStructure {
            name: ls.name.clone(),
            spins: ls.spins.clone(),
            structure: ls.structure.to_canonical_string().into(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LorentzStructure {
    pub name: SmartString<LazyCompact>,
    pub spins: Vec<isize>,
    pub structure: Atom,
}

impl LorentzStructure {
    pub(crate) fn from_serializable_lorentz_structure(
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
    value: Option<(F<f64>, F<f64>)>,
    expression: Option<SmartString<LazyCompact>>,
}

impl SerializableParameter {
    pub(crate) fn from_parameter(param: &Parameter) -> SerializableParameter {
        SerializableParameter {
            name: param.name.namespaceless_string().into(),
            lhablock: param.lhablock.clone(),
            lhacode: param.lhacode.clone(),
            nature: param.nature.clone(),
            parameter_type: param.parameter_type.clone(),
            value: param.value.map(|value| (value.re, value.im)),
            expression: param
                .expression
                .as_ref()
                .map(|a| a.to_canonical_string())
                .map(SmartString::from),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct ParameterName(pub UFOSymbol);

impl Deref for ParameterName {
    type Target = UFOSymbol;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[derive(Debug, Clone)]
pub struct Parameter {
    pub name: UFOSymbol,
    pub lhablock: Option<SmartString<LazyCompact>>,
    pub lhacode: Option<Vec<usize>>,
    pub nature: ParameterNature,
    pub parameter_type: ParameterType,
    pub value: Option<Complex<F<f64>>>,
    pub expression: Option<Atom>,
}

impl std::fmt::Display for Parameter {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.name)
    }
}

fn parameter_display_sort_key(parameter: &Parameter) -> (u8, String) {
    let nature_rank = match parameter.nature {
        ParameterNature::External => 0,
        ParameterNature::Internal => 1,
    };
    (nature_rank, parameter.name.to_string())
}

impl PartialEq for Parameter {
    fn eq(&self, other: &Self) -> bool {
        self.name == other.name
            && self.nature == other.nature
            && self.parameter_type == other.parameter_type
            // && self.value == other.value
            && self.expression == other.expression
            && self.lhablock == other.lhablock
            && self.lhacode == other.lhacode
    }
}

impl Eq for Parameter {}

impl Parameter {
    pub(crate) fn from_serializable_parameter(param: &SerializableParameter) -> Parameter {
        Parameter {
            name: (&param.name).into(),
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

    pub(crate) fn rep_rule(&self) -> Option<[Atom; 2]> {
        let lhs = self.name.into();
        let rhs = self.expression.clone();

        //Some([lhs, normalise_complex(&rhs?)])
        Some([lhs, rhs?])
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
    propagators: Vec<SerializablePropagator>,
    lorentz_structures: Vec<SerializableLorentzStructure>,
    couplings: Vec<SerializableCoupling>,
    vertex_rules: Vec<SerializableVertexRule>,
}

impl SerializableModel {
    pub(crate) fn from_file(file_path: impl AsRef<Path>) -> Result<SerializableModel, Report> {
        SmartSerde::from_file(file_path, "model")
    }
    pub(crate) fn from_str(s: String, format: &str) -> Result<SerializableModel, Report> {
        SmartSerde::from_str(s, format, "model")
    }

    pub(crate) fn from_model(model: &Model) -> SerializableModel {
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
                .values()
                .map(SerializableParameter::from_parameter)
                .collect(),
            particles: model
                .particles
                .iter()
                .map(|particle| SerializableParticle::from_particle(particle.0.as_ref()))
                .collect(),
            propagators: model
                .propagators
                .iter()
                .map(|propagator| SerializablePropagator::from_propagator(propagator.as_ref()))
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
                .values()
                .map(SerializableCoupling::from_coupling)
                .collect(),
            vertex_rules: model
                .vertex_rules
                .iter()
                .map(|vertex_rule| SerializableVertexRule::from_vertex_rule(vertex_rule.0.as_ref()))
                .collect(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Model {
    pub name: SmartString<LazyCompact>,
    pub restriction: Option<SmartString<LazyCompact>>,
    pub orders: Vec<Arc<Order>>,
    pub parameters: BTreeMap<ParameterName, Parameter>,
    pub particles: Vec<ArcParticle>,
    pub propagators: Vec<Arc<Propagator>>,
    pub lorentz_structures: Vec<Arc<LorentzStructure>>,
    pub couplings: BTreeMap<CouplingName, Coupling>,
    pub vertex_rules: Vec<ArcVertexRule>,
    pub unresolved_particles: HashMap<SmartString<LazyCompact>, HashSet<ArcParticle>>,
    pub particle_set_to_vertex_rules_map: HashMap<Vec<ArcParticle>, Vec<ArcVertexRule>>,
    pub order_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState>,
    pub lorentz_structure_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState>,
    pub particle_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState>,
    pub particle_pdg_to_position: HashMap<isize, usize, RandomState>,
    pub propagator_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState>,
    pub vertex_rule_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState>,
    pub particle_name_to_propagator_position: HashMap<SmartString<LazyCompact>, usize, RandomState>,
}

impl Default for Model {
    fn default() -> Self {
        Model {
            name: SmartString::<LazyCompact>::from("ModelNotLoaded"),
            restriction: None,
            orders: vec![],
            parameters: BTreeMap::new(),
            particles: vec![],
            propagators: vec![],
            lorentz_structures: vec![],
            couplings: BTreeMap::new(),
            vertex_rules: vec![],
            order_name_to_position:
                HashMap::<SmartString<LazyCompact>, usize, RandomState>::default(),
            lorentz_structure_name_to_position: HashMap::<
                SmartString<LazyCompact>,
                usize,
                RandomState,
            >::default(),
            unresolved_particles: HashMap::new(),
            particle_set_to_vertex_rules_map: HashMap::new(),
            particle_name_to_position:
                HashMap::<SmartString<LazyCompact>, usize, RandomState>::default(),
            particle_pdg_to_position: HashMap::<isize, usize, RandomState>::default(),
            propagator_name_to_position:
                HashMap::<SmartString<LazyCompact>, usize, RandomState>::default(),
            vertex_rule_name_to_position:
                HashMap::<SmartString<LazyCompact>, usize, RandomState>::default(),
            particle_name_to_propagator_position: HashMap::<
                SmartString<LazyCompact>,
                usize,
                RandomState,
            >::default(),
        }
    }
}
impl Model {
    pub fn apply_param_card(
        &mut self,
        input_param_card: &InputParamCard<F<f64>>,
    ) -> Result<(), Report> {
        input_param_card.apply_to_model(self)?;
        Ok(())
    }

    pub fn get_description(
        &self,
        show_particles: bool,
        show_parameters: bool,
        show_vertices: bool,
        show_couplings: bool,
    ) -> String {
        let name = self.name.clone().green().to_string();
        let restriction = match &self.restriction {
            Some(r) => r.clone().blue().to_string(),
            None => "None".into(),
        };
        let coupling_orders_value = if self.orders.is_empty() {
            "None".into()
        } else {
            self.orders
                .iter()
                .map(|order| {
                    format!(
                        "{} (expansion order: {}, hierarchy: {})",
                        order.name.green(),
                        order.expansion_order,
                        order.hierarchy
                    )
                })
                .collect::<Vec<_>>()
                .join(", ")
        };

        let particle_list = if !show_particles {
            "[ hiddem ]".blue().to_string()
        } else {
            let mut particle_table = Builder::new();

            particle_table.push_record([
                "Name".green().to_string(),
                "PDG code".normal().to_string(),
                "Mass name".blue().to_string(),
                "Mass value".normal().to_string(),
                "Width name".blue().to_string(),
                "Width value".normal().to_string(),
            ]);
            for p in &self.particles {
                let mass_value = if let Some(value) = self.get_parameter(p.mass.0.to_string()).value
                {
                    format!("{:.6}", value.re)
                } else {
                    "None".into()
                };
                let width_value =
                    if let Some(value) = self.get_parameter(p.width.0.to_string()).value {
                        format!("{:.6}", value.re)
                    } else {
                        "None".into()
                    };

                particle_table.push_record(&[
                    p.name.to_string().green().to_string(),
                    format!("{:+}", p.pdg_code).normal().to_string(),
                    p.mass.0.to_string().blue().to_string(),
                    mass_value.normal().to_string(),
                    p.width.0.to_string().blue().to_string(),
                    width_value.normal().to_string(),
                ]);
            }
            let mut particle_table_built = particle_table.build();

            format!(
                "\n{}\n",
                particle_table_built.with(Style::rounded()).to_string()
            )
        };

        let parameter_list = if !show_parameters {
            "[ hidden ]".blue().to_string()
        } else {
            let mut parameter_table = Builder::new();
            parameter_table.push_record([
                "Name".green().to_string(),
                "Nature".normal().to_string(),
                "Type".normal().to_string(),
                "Value".normal().to_string(),
                "Expression".to_string(),
            ]);
            let mut parameters = self.parameters.values().collect_vec();
            parameters.sort_by_key(|parameter| parameter_display_sort_key(parameter));

            for param in parameters {
                parameter_table.push_record(&[
                    param.name.to_string().green().to_string(),
                    if param.nature == ParameterNature::External {
                        format!("{:?}", param.nature).green().to_string()
                    } else {
                        format!("{:?}", param.nature).yellow().to_string()
                    },
                    if param.parameter_type == ParameterType::Real {
                        format!("{:?}", param.parameter_type).green().to_string()
                    } else {
                        format!("{:?}", param.parameter_type).yellow().to_string()
                    },
                    if let Some(value) = param.value {
                        format!("{:.6}", value.re)
                    } else {
                        "".into()
                    },
                    if let Some(expr) = &param.expression {
                        expr.to_string()
                    } else {
                        "".into()
                    },
                ]);
            }
            format!(
                "\n{}\n",
                parameter_table.build().with(Style::rounded()).to_string()
            )
        };

        let vertex_list = if !show_vertices {
            "[ hidden ]".blue().to_string()
        } else {
            let mut vertex_table = Builder::new();

            let max_n_particles = self
                .vertex_rules
                .iter()
                .map(|vr| vr.particles.len())
                .max()
                .unwrap_or(3);
            let mut header = vec!["Name".green().to_string(), "Particles".blue().to_string()];
            for _ in 0..max_n_particles - 1 {
                header.push("".to_string());
            }
            vertex_table.push_record(header);
            // for vr in &self.vertex_rules {
            //     vertex_table.push_record(&[
            //         vr.name.to_string().green().to_string(),
            //         vr.particles
            //             .iter()
            //             .map(|p| format!("{:<6}", p.name.to_string()).blue().to_string())
            //             .collect::<Vec<_>>()
            //             .join(", "),
            //     ]);
            // }
            for vr in &self.vertex_rules {
                let mut record = vec![vr.name.to_string().green().to_string()];
                for p in &vr.particles {
                    record.push(p.name.to_string().blue().to_string());
                }
                // Pad the record with empty strings if necessary
                while record.len() < max_n_particles + 1 {
                    record.push("".to_string());
                }
                vertex_table.push_record(&record);
            }
            let mut vertex_table_built = vertex_table.build();

            vertex_table_built
                .with(
                    Style::rounded()
                        .remove_vertical()
                        // keep only the split after column 0 (between col 0 and 1)
                        .verticals([(1, VerticalLine::inherit(Style::rounded()))]),
                )
                .with(Modify::new(Cell::new(0, 1)).with(Span::column(max_n_particles as isize)));
            format!("\n{}\n", vertex_table_built.to_string())
        };

        let coupling_list = if !show_couplings {
            "[ hidden ]".blue().to_string()
        } else {
            let mut coupling_table = Builder::new();
            coupling_table.push_record([
                "Name".green().to_string(),
                "Orders".blue().to_string(),
                "Expression".normal().to_string(),
            ]);
            for c in self.couplings.values() {
                coupling_table.push_record(&[
                    c.name.to_string().green().to_string(),
                    c.orders
                        .iter()
                        .map(|(k, v)| format!("{}={}", k.blue(), v))
                        .collect::<Vec<_>>()
                        .join(" "),
                    c.expression.to_string(),
                ]);
            }
            format!(
                "\n{}\n",
                coupling_table.build().with(Style::rounded()).to_string(),
            )
        };

        #[rustfmt::skip]
        return format!("
{model_name_label:<30}: {name}
{restriction_label:<30}: {restriction}
{coupling_orders_label:<30}: {coupling_orders_value}

{n_particles} particles : {particle_list}
{n_parameters} parameters :{parameter_list}
{n_vertices} vertices : {vertex_list}
{n_couplings} couplings : {coupling_list}
",
model_name_label = "Model name",
restriction_label = "Restriction",
coupling_orders_label = "Coupling orders",
n_particles = format!("{}", self.particles.len()).green(),
n_parameters = format!("{}", self.parameters.len()).green(),
n_vertices = format!("{}", self.vertex_rules.len()).green(),
n_couplings = format!("{}", self.couplings.len()).green(),
);
    }

    /// Generate edge-style.typ template file with styles for all particles in the model
    pub fn generate_edge_style_template(
        &self,
        template_path: impl AsRef<std::path::Path>,
    ) -> Result<(), std::io::Error> {
        use std::fs;

        // Create the directory if it doesn't exist
        if let Some(parent) = template_path.as_ref().parent() {
            fs::create_dir_all(parent)?;
        }

        let mut edge_style_content = String::new();
        edge_style_content.push_str(r#"#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, cetz,edge,hide
#import "@preview/mitex:0.2.6": *

#let massive = 1mm
#let massless = 0.5mm
#let source_stroke(c:black, thickness:0.5mm,dash:none) = (stroke:(paint:c,thickness:thickness)+dash)
#let sink_stroke(c:black, thickness:0.5mm,dash:none) = (stroke:source_stroke(c:c.lighten(50%), thickness:thickness,dash:dash).stroke)
#let wave = (decorations:cetz.decorations.wave.with(amplitude: 4pt,segment-length:0.2))
#let double = (extrude:(-0.5mm, 0.5mm))
#let arrow = (marks:((inherit:"solid",rev:false,pos:1.1,scale:50%),))
#let antiarrow = (marks:((inherit:"solid",rev:true,pos:1.1,scale:50%),))
#let arrowmap = orientation => if orientation == "Default"{
  arrow
} else if orientation == "Reversed"{
  antiarrow
} else{
  (:)
}
#let coil = (decorations:cetz.decorations.coil.with(amplitude: 4pt,segment-length:0.2))
#let zigzag = (decorations:cetz.decorations.zigzag.with(amplitude: 4pt,segment-length:0.2))
#let dashed = (dash: (0.1em, 0.5em))
#let dotted = (dash: (0.01em, 0.3em))
// Auto-generated particle styles from model (computed in Rust)
#let map = (
"#);

        // Generate styles for all particles in the model
        for particle in self.particles.iter() {
            edge_style_content.push_str(&format!(
                r#"  "{}": {},
"#,
                particle.name,
                particle.generate_edge_typst_dict()
            ));
        }

        edge_style_content.push_str(")\n");

        fs::write(&template_path, edge_style_content)?;
        info!(
            "Generated dynamic edge styles for {} particles",
            self.particles.len()
        );

        Ok(())
    }

    pub fn simplify(
        &mut self,
        model_parameters: &mut InputParamCard<F<f64>>,
    ) -> Result<(), Report> {
        self.apply_param_card(model_parameters)?;
        self.recompute_dependents()?;

        // Remove zero parameters from the input card
        model_parameters.data = model_parameters
            .data
            .iter()
            .filter(|(_, v)| v.re != F::<f64>::from_f64(0.0) || v.im != F::<f64>::from_f64(0.0))
            .map(|(k, v)| (*k, *v))
            .collect::<HashMap<UFOSymbol, Complex<F<f64>>>>();

        // Set all external parameters with value 0 to constant internal parameters with expression zero.
        let mut removed_parameters = vec![];
        for (_p_name, param) in self.parameters.iter_mut() {
            if param.nature == ParameterNature::External
                && let Some(value) = param.value
                && value == Complex::new(F(0.0), F(0.0))
            {
                param.value = Some(Complex::new(F(0.0), F(0.0)));
                param.expression = Some(parse!("UFO::ZERO"));
                param.nature = ParameterNature::Internal;
                removed_parameters.push(param.name);
            }
        }

        if !removed_parameters.is_empty() {
            info!(
                "The following {} external parameters were forced to zero by the restriction card:\n{}",
                format!("{}", removed_parameters.len()).green(),
                removed_parameters
                    .iter()
                    .map(|p| p.to_string())
                    .collect::<Vec<_>>()
                    .join(", ")
                    .to_string()
                    .blue(),
            );
        }

        // Remove all vertices with zero couplings
        let mut retained_vertex_rules = vec![];
        let mut removed_vertex_rules = vec![];
        for v in self.vertex_rules.iter() {
            let mut new_vr = v.0.as_ref().clone();
            for row in new_vr.couplings.iter_mut() {
                *row = row
                    .iter()
                    .map(|c_opt| {
                        if let Some(c) = c_opt {
                            if let Some(cpl) = self.couplings.get(c) {
                                if let Some(value) = cpl.value {
                                    if value == Complex::new(0.0, 0.0) {
                                        None
                                    } else {
                                        Some(*c)
                                    }
                                } else {
                                    Some(*c)
                                }
                            } else {
                                Some(*c)
                            }
                        } else {
                            None
                        }
                    })
                    .collect::<Vec<_>>();
            }
            if new_vr
                .couplings
                .iter()
                .all(|row| row.iter().all(|c| c.is_none()))
            {
                removed_vertex_rules.push(new_vr);
            } else {
                retained_vertex_rules.push(new_vr);
            }
        }

        self.vertex_rules = retained_vertex_rules
            .into_iter()
            .map(|vr| ArcVertexRule(Arc::new(vr)))
            .collect::<Vec<_>>();

        if !removed_vertex_rules.is_empty() {
            info!(
                "The following {} vertex rules were removed by the restriction card:\n{}",
                format!("{}", removed_vertex_rules.len()).green(),
                removed_vertex_rules
                    .iter()
                    .map(|v| format!(
                        "{} -> ({})",
                        v.name,
                        v.particles
                            .iter()
                            .map(|p| p.name.to_string())
                            .collect::<Vec<_>>()
                            .join(", ")
                    ))
                    .collect::<Vec<_>>()
                    .join(" | ")
                    .to_string()
                    .blue(),
            );
        }

        let removed_couplings = self
            .couplings
            .iter()
            .filter(|(_, cpl)| cpl.value == Some(Complex::new(0.0, 0.0)))
            .map(|(name, _)| *name)
            .collect::<Vec<_>>();

        // Now remove all couplings that are zero
        self.couplings
            .retain(|_, cpl| cpl.value != Some(Complex::new(0.0, 0.0)));

        if !removed_couplings.is_empty() {
            info!(
                "The following {} couplings were removed by the restriction card:\n{}",
                format!("{}", removed_couplings.len()).green(),
                removed_couplings
                    .iter()
                    .map(|c| c.to_string())
                    .collect::<Vec<_>>()
                    .join(", ")
                    .to_string()
                    .blue(),
            );
        }

        self.update_name_dictionaries();

        Ok(())
    }

    pub fn default_param_card(&self) -> InputParamCard<F<f64>> {
        InputParamCard::default_from_model(self)
    }

    pub fn contains_symbol(&self, symbol: &UFOSymbol) -> bool {
        self.couplings.contains_key(&CouplingName(*symbol))
            || self.parameters.contains_key(&ParameterName(*symbol))
    }
    pub fn get_symbol_value(&self, symbol: UFOSymbol) -> Option<Complex<F<f64>>> {
        if let Some(cpl) = self.couplings.get(&CouplingName(symbol)) {
            return cpl.value.map(|a| a.map(F));
        }
        if let Some(param) = self.parameters.get(&ParameterName(symbol))
            && let Some(value) = param.value
        {
            return Some(value);
        }
        None
    }

    fn parameters_to_empty_fns(&self) -> Vec<Replacement> {
        let mut reps = vec![];
        for n in self.couplings.keys() {
            reps.push(Replacement::new(
                Atom::from(n.0).to_pattern(),
                function!(n.0.0),
            ))
        }
        for n in self.parameters.keys() {
            reps.push(Replacement::new(
                Atom::from(n.0).to_pattern(),
                function!(n.0.0),
            ))
        }

        reps
    }

    pub fn recompute_dependents(&mut self) -> Result<()> {
        let mut fn_map = FunctionMap::new();
        let reps = self.parameters_to_empty_fns();

        let mut expr = vec![];
        let mut new_values_len = 0;

        for (n, c) in &self.couplings {
            let key = n.0.0;
            expr.push(function!(key));

            fn_map
                .add_function::<Symbol>(
                    key,
                    c.name.namespaceless_string().into(),
                    vec![],
                    c.expression.replace_multiple(&reps),
                )
                .map_err(|e| eyre!(" {}", e))?;
            new_values_len += 1;
        }

        let mut params = vec![];
        let mut param_values = vec![];

        for (n, p) in &self.parameters {
            let key = function!(n.0.0);
            match p.nature {
                ParameterNature::External => {
                    params.push(key);
                    if let Some(value) = p.value {
                        param_values.push(value);
                    } else {
                        return Err(eyre!("External parameter {} has no value", p.name));
                    }
                }
                ParameterNature::Internal => {
                    new_values_len += 1;
                    expr.push(key.clone());
                    if let Some(body) = p.expression.clone() {
                        fn_map
                            .add_function::<Symbol>(
                                n.0.0,
                                p.name.namespaceless_string().into(),
                                vec![],
                                body.replace_multiple(&reps),
                            )
                            .map_err(|e| eyre!(" {}", e))?;
                    } else {
                        let value = p
                            .value
                            .ok_or(eyre!("internal param {} has no expression or value", key))?;
                        let value_rat = (
                            Fraction::<IntegerRing>::try_from(value.re.0).unwrap(),
                            Fraction::<IntegerRing>::try_from(value.im.0).unwrap(),
                        )
                            .into();
                        fn_map.add_constant(key, value_rat);
                    }
                }
            }
        }

        fn_map
            .add_external_function(
                UFOSymbol::from("complexconjugate").0,
                "complexconjugate".into(),
            )
            .unwrap();

        // fn_map.add_external_function(name, rename)

        fn_map.add_constant(
            Atom::var(Symbol::PI),
            (Rational::try_from(0.0.pi()).unwrap(), Rational::zero()).into(),
        );

        let evaluator = AtomView::to_eval_tree_multiple(&expr, &fn_map, &params)
            .unwrap()
            .linearize(&OptimizationSettings {
                cpe_iterations: Some(1),
                verbose: false,
                ..OptimizationSettings::default()
            });
        let mut ext: ExternalFunctionMap = HashMap::default();
        ext.insert(
            "complexconjugate".to_string(),
            Box::new(|a| {
                if a.len() > 1 {
                    panic!("complex_conjugate takes one argument, got {}", a.len());
                };
                a[0].conj()
            }),
        );
        let mut evaluator = evaluator
            .map_coeff(&|f| Complex::new(F(f.re.to_f64()), F(f.im.to_f64())))
            .with_external_functions(ext)
            .unwrap();

        let mut new_values = vec![Complex::new(F(0.0), F(0.0)); new_values_len];
        evaluator.evaluate(&param_values, &mut new_values);

        for (i, c) in self.couplings.values_mut().enumerate() {
            c.value = Some(new_values[i].map(|f| f.0));
        }
        let mut i = self.couplings.len();

        for c in self.parameters.values_mut() {
            if c.nature == ParameterNature::External {
                continue;
            }
            c.value = Some(new_values[i]);
            i += 1;
        }

        Ok(())
    }

    pub(crate) fn generate_params(&self) -> Vec<Atom> {
        let mut params = vec![];

        for cpl in self.couplings.values().filter(|c| c.value.is_some()) {
            if cpl.value.is_some() {
                params.push(cpl.name.into());
            }
        }
        for param in self.parameters.values().filter(|p| p.value.is_some()) {
            if param.value.is_some() {
                let name = param.name.into();
                params.push(name);
            }
        }

        params
    }

    pub fn is_empty(&self) -> bool {
        self.name == "ModelNotLoaded" || self.particles.is_empty()
    }

    pub fn apply_coupling_replacement_rules(&self, a: &Atom) -> Atom {
        let mut reps = vec![];
        for cpl in self.couplings.values() {
            let [a, b] = cpl.rep_rule();
            reps.push(Replacement::new(a.to_pattern(), b));
        }

        a.replace_multiple(&reps)
    }
    pub fn apply_parameter_replacement_rules(&self, a: &Atom) -> Atom {
        let mut reps = vec![];
        for p in self.parameters.values() {
            let Some([a, b]) = p.rep_rule() else {
                continue;
            };
            reps.push(Replacement::new(a.to_pattern(), b));
        }

        a.replace_multiple(&reps)
    }

    pub fn export_coupling_replacement_rules(
        &self,
        export_root: &str,
        print_ops: PrintOptions,
    ) -> Result<(), Report> {
        let path = Path::new(export_root).join("sources").join("model");

        if !path.exists() {
            fs::create_dir_all(&path)?;
        }
        let mut reps = Vec::new();

        for cpl in self.couplings.values() {
            reps.push(
                cpl.rep_rule()
                    .map(|a| format!("{}", AtomPrinter::new_with_options(a.as_view(), print_ops))),
            );
        }

        for para in self.parameters.values() {
            if let Some(rule) = para.rep_rule() {
                reps.push(
                    rule.map(|a| {
                        format!("{}", AtomPrinter::new_with_options(a.as_view(), print_ops))
                    }),
                );
            }
        }

        fs::write(
            path.join("model_replacements.json"),
            serde_json::to_string_pretty(&reps)?,
        )?;

        Ok(())
    }

    fn generate_particle_set_to_vertex_rules_map(&mut self) {
        let mut map = HashMap::new();

        for vertex in self.vertex_rules.iter() {
            let mut vertex_particles = vertex.0.particles.clone();
            vertex_particles.sort();
            map.entry(vertex_particles)
                .and_modify(|l: &mut Vec<ArcVertexRule>| l.push(vertex.clone()))
                .or_insert(vec![vertex.clone()]);
        }
        self.particle_set_to_vertex_rules_map = map;
    }

    fn generate_unresolved_particles(&mut self) {
        let mut map = HashMap::new();

        for v in &self.vertex_rules {
            let mut set = HashSet::default();
            for p in &v.0.particles {
                if p.0.is_massless() {
                    set.insert(p.clone());
                }
            }
            for (k, _) in v.0.coupling_orders(self) {
                let current_set = map.entry(k).or_insert(HashSet::<ArcParticle>::default());

                set.iter().for_each(|d| {
                    current_set.insert(d.clone());
                });
            }
        }

        self.unresolved_particles = map;
    }

    pub(crate) fn update_name_dictionaries(&mut self) {
        self.order_name_to_position = self
            .orders
            .iter()
            .enumerate()
            .map(|(i, o)| (o.name.clone(), i))
            .collect();

        self.lorentz_structure_name_to_position = self
            .lorentz_structures
            .iter()
            .enumerate()
            .map(|(i, ls)| (ls.name.clone(), i))
            .collect();

        self.particle_name_to_position = self
            .particles
            .iter()
            .enumerate()
            .map(|(i, p)| (p.0.name.clone(), i))
            .collect();

        self.particle_pdg_to_position = self
            .particles
            .iter()
            .enumerate()
            .map(|(i, p)| (p.0.pdg_code, i))
            .collect();

        self.propagator_name_to_position = self
            .propagators
            .iter()
            .enumerate()
            .map(|(i, pr)| (pr.name.clone(), i))
            .collect();

        self.vertex_rule_name_to_position = self
            .vertex_rules
            .iter()
            .enumerate()
            .map(|(i, vr)| (vr.0.name.clone(), i))
            .collect();
    }

    pub(crate) fn from_serializable_model(serializable_model: SerializableModel) -> Model {
        //initialize the UFO and ETS symbols

        // let _ = *UFO;
        let _ = *ETS;

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
            .map(|serializable_param| {
                let parameter = Parameter::from_serializable_parameter(serializable_param);

                (ParameterName(parameter.name), parameter)
            })
            .collect();

        // Extract particles
        model.particles = serializable_model
            .particles
            .iter()
            .enumerate()
            .map(|(i_part, serializable_particle)| {
                let particle =
                    Arc::new(Particle::from_serializable_particle(serializable_particle));
                model
                    .particle_name_to_position
                    .insert(particle.name.clone(), i_part);
                model
                    .particle_pdg_to_position
                    .insert(particle.pdg_code, i_part);
                ArcParticle(particle)
            })
            .collect();

        // Extract propagators

        model.propagators = serializable_model
            .propagators
            .iter()
            .enumerate()
            .map(|(i_prop, serializable_propagator)| {
                let propagator = Arc::new(Propagator::from_serializable_propagator(
                    &model,
                    serializable_propagator,
                ));
                model
                    .propagator_name_to_position
                    .insert(propagator.name.clone(), i_prop);
                propagator
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
            .map(|serializable_coupling| {
                let coupling = Coupling::from_serializable_coupling(serializable_coupling);

                (CouplingName(coupling.name), coupling)
            })
            .collect();

        // Extract vertex rules
        model.vertex_rules = serializable_model
            .vertex_rules
            .iter()
            .enumerate()
            .map(|(i_vr, serializable_vertex_rule)| {
                let vertex_rule = ArcVertexRule(Arc::new(
                    VertexRule::from_serializable_vertex_rule(&model, serializable_vertex_rule),
                ));
                model
                    .vertex_rule_name_to_position
                    .insert(vertex_rule.0.name.clone(), i_vr);
                vertex_rule
            })
            .collect();

        // Set propagator mapping
        model.particle_name_to_propagator_position =
            serializable_model.propagators.iter().enumerate().fold(
                HashMap::<SmartString<LazyCompact>, usize, RandomState>::default(),
                |mut map, (i_prop, serializable_propagator)| {
                    map.insert(serializable_propagator.particle.clone(), i_prop);
                    map
                },
            );

        model.generate_unresolved_particles();
        model.generate_particle_set_to_vertex_rules_map();

        model
    }

    pub fn to_serializable(&self) -> SerializableModel {
        SerializableModel::from_model(self)
    }

    pub fn from_file(file_path: impl AsRef<Path>) -> Result<Model, Report> {
        let mut model =
            SerializableModel::from_file(file_path).map(Model::from_serializable_model)?;

        model.recompute_dependents()?;
        Ok(model)
    }

    pub fn from_str(s: String, format: &str) -> Result<Model, Report> {
        let mut model =
            SerializableModel::from_str(s, format).map(Model::from_serializable_model)?;

        model.recompute_dependents()?;
        Ok(model)
    }

    #[inline]
    pub(crate) fn get_propagator_for_particle<S: AsRef<str>>(&self, name: S) -> Arc<Propagator> {
        if let Some(position) = self.particle_name_to_propagator_position.get(name.as_ref()) {
            self.propagators[*position].clone()
        } else {
            panic!(
                "Propagator for particle '{}' not found in model '{}'. Valid entries are:\n{}",
                name.as_ref(),
                self.name,
                self.particle_name_to_propagator_position.keys().join(", ")
            );
        }
    }

    #[inline]
    pub fn get_particle<S: AsRef<str>>(&self, name: S) -> ArcParticle {
        if let Some(position) = self.particle_name_to_position.get(name.as_ref()) {
            self.particles[*position].clone()
        } else {
            panic!(
                "Particle '{}' not found in model '{}'. Valid entries are:\n{}",
                name.as_ref(),
                self.name,
                self.particle_name_to_position.keys().join(", ")
            );
        }
    }

    #[inline]
    pub fn try_get_particle<S: AsRef<str>>(&self, name: S) -> Result<ArcParticle> {
        if let Some(position) = self.particle_name_to_position.get(name.as_ref()) {
            Ok(self.particles[*position].clone())
        } else {
            Err(eyre!(
                "Particle '{}' not found in model '{}'. Valid entries are:\n{}",
                name.as_ref(),
                self.name,
                self.particle_name_to_position.keys().join(", ")
            ))
        }
    }
    #[inline]
    pub(crate) fn get_particle_from_pdg(&self, pdg: isize) -> ArcParticle {
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
    pub(crate) fn try_get_particle_from_pdg(&self, pdg: isize) -> Result<ArcParticle> {
        if let Some(position) = self.particle_pdg_to_position.get(&pdg) {
            Ok(self.particles[*position].clone())
        } else {
            Err(eyre!(
                "Particle with PDG {} not found in model '{}'.",
                pdg,
                self.name
            ))
        }
    }

    #[inline]
    pub fn get_propagator<S: AsRef<str>>(&self, name: S) -> ArcPropagator {
        if let Some(position) = self.propagator_name_to_position.get(name.as_ref()) {
            ArcPropagator(self.propagators[*position].clone())
        } else {
            panic!(
                "Propagator '{}' not found in model '{}'.",
                name.as_ref(),
                self.name
            );
        }
    }

    #[inline]
    pub fn get_parameter_opt<S: AsRef<str>>(&self, name: S) -> Option<&Parameter> {
        if let Some(position) = self
            .parameters
            .get(&ParameterName(UFOSymbol::from(name.as_ref())))
        {
            Some(position)
        } else {
            None
        }
    }

    #[inline]
    pub fn get_parameter<S: AsRef<str>>(&self, name: S) -> &Parameter {
        if let Some(position) = self
            .parameters
            .get(&ParameterName(UFOSymbol::from(name.as_ref())))
        {
            position
        } else {
            panic!(
                "Parameter '{}' not found in model '{}'.",
                name.as_ref(),
                self.name
            );
        }
    }

    #[inline]
    pub fn get_parameter_mut_opt<S: AsRef<str>>(&mut self, name: S) -> Option<&mut Parameter> {
        if let Some(position) = self
            .parameters
            .get_mut(&ParameterName(UFOSymbol::from(name.as_ref())))
        {
            Some(position)
        } else {
            None
        }
    }

    #[inline]
    pub fn get_parameter_mut<S: AsRef<str>>(&mut self, name: S) -> Result<&mut Parameter> {
        if let Some(position) = self
            .parameters
            .get_mut(&ParameterName(UFOSymbol::from(name.as_ref())))
        {
            Ok(position)
        } else {
            Err(eyre!(
                "Parameter '{}' not found in model '{}'.",
                name.as_ref(),
                self.name
            ))
        }
    }

    #[inline]
    pub fn get_order<S: AsRef<str>>(&self, name: S) -> Arc<Order> {
        if let Some(position) = self.order_name_to_position.get(name.as_ref()) {
            self.orders[*position].clone()
        } else {
            panic!(
                "Coupling order '{}' not found in model '{}'.",
                name.as_ref(),
                self.name
            );
        }
    }
    #[inline]
    pub fn get_lorentz_structure<S: AsRef<str>>(&self, name: S) -> Arc<LorentzStructure> {
        if let Some(position) = self.lorentz_structure_name_to_position.get(name.as_ref()) {
            self.lorentz_structures[*position].clone()
        } else {
            panic!(
                "Lorentz structure '{}' not found in model '{}'.",
                name.as_ref(),
                self.name
            );
        }
    }
    #[inline]
    pub fn get_coupling<S: AsRef<str>>(&self, name: S) -> &Coupling {
        if let Some(coupling) = self
            .couplings
            .get(&CouplingName(UFOSymbol::from(name.as_ref())))
        {
            coupling
        } else {
            panic!(
                "Coupling '{}' not found in model '{}'.",
                name.as_ref(),
                self.name
            );
        }
    }
    #[inline]
    pub fn get_vertex_rule<S: AsRef<str>>(&self, name: S) -> ArcVertexRule {
        if let Some(position) = self.vertex_rule_name_to_position.get(name.as_ref()) {
            self.vertex_rules[*position].clone()
        } else {
            panic!(
                "Vertex rule '{}' not found in model '{}'.",
                name.as_ref(),
                self.name
            );
        }
    }
}

#[cfg(test)]
mod tests {

    use crate::{
        model::{
            ArcPropagator, ArcVertexRule, Parameter, ParameterName, ParameterNature, ParameterType,
        },
        momentum::{Helicity, ThreeMomentum},
        utils::{F, load_generic_model},
    };

    use super::{ArcParticle, Model, UFOSymbol, parameter_display_sort_key};

    #[test]
    fn test_encode_decode_arc_particle() {
        let model = load_generic_model("scalars");
        let particle = model.get_particle("scalar_0");
        let particle_encoded =
            bincode::encode_to_vec(&particle, bincode::config::standard()).unwrap();

        let particle_decoded: ArcParticle = bincode::decode_from_slice_with_context(
            &particle_encoded,
            bincode::config::standard(),
            model,
        )
        .unwrap()
        .0;

        assert_eq!(particle, particle_decoded);
    }

    #[test]
    fn test_encode_decode_vertex_rule() {
        let model = load_generic_model("sm");
        let vertex = model.get_vertex_rule("V_141");
        let vertex_encoded = bincode::encode_to_vec(&vertex, bincode::config::standard()).unwrap();

        let vertex_decoded: ArcVertexRule = bincode::decode_from_slice_with_context(
            &vertex_encoded,
            bincode::config::standard(),
            model,
        )
        .unwrap()
        .0;
        assert_eq!(vertex, vertex_decoded);
    }

    #[test]
    fn test_encode_decode_arc_propagator() {
        let model = load_generic_model(
            "sm
        ",
        );
        let propagator = model.get_propagator("t_propFeynman");
        let propagator_encoded =
            bincode::encode_to_vec(&propagator, bincode::config::standard()).unwrap();

        let propagator_decoded: ArcPropagator = bincode::decode_from_slice_with_context(
            &propagator_encoded,
            bincode::config::standard(),
            model,
        )
        .unwrap()
        .0;

        assert_eq!(propagator.name, propagator_decoded.name);
    }

    #[test]
    fn test_pol_limit_vector() {
        let vals = [0., 3., 4.].into_iter().map(F).collect::<Vec<_>>();
        let mom = ThreeMomentum::new(vals[0], vals[1], vals[2]).into_on_shell_four_momentum(None); //Some(F(8.66025)));

        let hel = Helicity::MINUS;
        println!("mom:{mom}");
        println!("hel{hel}");
        println!("{}", mom.eps_pol(hel));
        println!("{}", mom.eps_pol(hel).bar());
        // let vals = [0.01, 0., 1.].into_iter().map(F).collect::<Vec<_>>();
        // let mom = ThreeMomentum::new(vals[0], vals[1], vals[2]).into_on_shell_four_momentum(None);

        let hel = Helicity::PLUS;
        println!("mom:{mom}");
        println!("hel{hel}");
        println!("{}", mom.eps_pol(hel));
        println!("{}", mom.eps_pol(hel).bar());
        let hel = Helicity::MINUS;

        println!("hel{hel}");
        println!("{}", mom.eps_pol(hel));
        println!("{}", mom.eps_pol(hel).bar());

        let vals = [0., 0., -1.].into_iter().map(F).collect::<Vec<_>>();
        let mom = ThreeMomentum::new(vals[0], vals[1], vals[2]).into_on_shell_four_momentum(None);

        let hel = Helicity::PLUS;

        println!("mom:{mom}");
        println!("hel{hel}");

        println!("{}", mom.eps_pol(hel));
        println!("{}", mom.eps_pol(hel).bar());

        let hel = Helicity::MINUS;

        println!("hel{hel}");
        println!("{}", mom.eps_pol(hel));
        println!("{}", mom.eps_pol(hel).bar());
    }

    #[test]
    fn display_orders_external_parameters_first() {
        let mut model = Model::default();
        let external = Parameter {
            name: UFOSymbol::from("zeta"),
            lhablock: None,
            lhacode: None,
            nature: ParameterNature::External,
            parameter_type: ParameterType::Real,
            value: None,
            expression: None,
        };
        let internal = Parameter {
            name: UFOSymbol::from("alpha"),
            lhablock: None,
            lhacode: None,
            nature: ParameterNature::Internal,
            parameter_type: ParameterType::Real,
            value: None,
            expression: None,
        };
        model
            .parameters
            .insert(ParameterName(external.name.clone()), external.clone());
        model
            .parameters
            .insert(ParameterName(internal.name.clone()), internal.clone());

        assert!(parameter_display_sort_key(&external) < parameter_display_sort_key(&internal));

        let description = model.get_description(false, true, false, false);
        let external_index = description.find("zeta").unwrap();
        let internal_index = description.find("alpha").unwrap();
        assert!(external_index < internal_index);
    }
}
