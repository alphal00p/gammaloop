#![allow(dead_code)]

use aind::Aind;
use idenso::color::ColorSimplifier;
use idenso::gamma::GammaSimplifier;
use idenso::representations::Bispinor;
use linnet::half_edge::involution::EdgeIndex;
use log::warn;
use schemars::JsonSchema;

use spenso::network::library::{DummyLibrary, TensorLibraryData};
use spenso::network::parsing::{ParseSettings, ShadowedStructure};
use spenso::network::store::NetworkStore;
use spenso::network::{ContractScalars, Sequential, SingleSmallestDegree, SmallestDegree, Steps};
use spenso::shadowing::symbolica_utils::SerializableAtom;
use spenso::shadowing::symbolica_utils::SerializableSymbol;

use spenso::tensors::data::DataTensor;
use spenso::tensors::data::GetTensorData;
use spenso::tensors::data::StorageTensor;
use spenso::tensors::parametric::MixedTensor;
use spenso::tensors::parametric::atomcore::TensorAtomMaps;

use spenso::tensors::parametric::ParamTensor;
use spenso::tensors::parametric::TensorSet;
use std::fmt::Debug;
use std::ops::Deref;
use std::sync::{Arc, Mutex};
use symbolica_ext::AtomCoreExt;
use thiserror::Error;
use tracing::{debug, instrument};
// use crate::feyngen::dis::{DisEdge, DisVertex};

use crate::momentum::{PolDef, PolType};
use crate::utils::{FUN_LIB, GS, TENSORLIB, W_};
use crate::{
    model::Model,
    utils::{F, serde_utils::IsDefault},
};

use crate::{GammaLoopContextContainer, disable};
use ahash::AHashMap;
use bincode::{Decode, Encode};
use color_eyre::{Report, Result};
use eyre::eyre;
// use gxhash::GxBuildHasher;
use itertools::Itertools;

use serde::de::DeserializeOwned;
use serde::ser::SerializeStruct;
use serde::{Deserialize, Serialize};

use spenso::contraction::Contract;

use spenso::network::library::symbolic::{ETS, ExplicitKey};

use spenso::structure::concrete_index::{ExpandedIndex, FlatIndex};

use spenso::structure::representation::{LibraryRep, Minkowski};
use spenso::structure::{HasStructure, ScalarTensor, SmartShadowStructure};

use spenso::{
    shadowing::Shadowable,
    structure::{
        NamedStructure, TensorStructure,
        representation::{Lorentz, RepName},
    },
};

use symbolica::domains::rational::Rational;
use symbolica::poly::PolyVariable;
use symbolica::printer::PrintOptions;
use symbolica::state::Workspace;

use crate::numerator::ufo::UFO;
use symbolica::atom::{AtomCore, AtomOrView, AtomView, Symbol};
use symbolica::evaluate::{CompileOptions, InlineASM};
use symbolica::{
    atom::{Atom, FunctionBuilder},
    function, parse, symbol,
};

pub mod symbolica_ext;

use symbolica::{evaluate::FunctionMap, id::Replacement};
pub mod aind;
pub mod ufo;
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq)]
/// Settings for the numerator
pub struct NumeratorSettings {
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub eval_settings: NumeratorEvaluatorOptions,
    /// Parse mode for the numerator, once all processing is done. `Polynomial` turns it into a polynomial in the energies, while `Direct` keeps it as is
    pub parse_mode: NumeratorParseMode,
    /// If set, dump the expression the expression at each step in this format
    pub dump_expression: Option<ExpressionFormat>,
    /// If set, instead of deriving the numerator from feynman rules, use this as the numerator
    /// Will be parsed to a symbolica expression
    // pub global_numerator: Option<String>,
    /// If set, multiply the numerator by this prefactor
    // pub global_prefactor: GlobalPrefactor,
    /// Type of Gamma algebra to use, either symbolic (replacement rules) or concrete (replace by value using spenso)
    pub gamma_algebra: GammaAlgebraMode,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, Encode, Decode, PartialEq)]
pub enum ExpressionFormat {
    Mathematica,
    #[default]
    Symbolica,
}

impl From<ExpressionFormat> for PrintOptions {
    fn from(value: ExpressionFormat) -> Self {
        match value {
            ExpressionFormat::Symbolica => PrintOptions::file(),
            ExpressionFormat::Mathematica => PrintOptions::mathematica(),
        }
    }
}

impl Default for NumeratorSettings {
    fn default() -> Self {
        NumeratorSettings {
            eval_settings: Default::default(),
            // global_numerator: None,
            // global_prefactor: GlobalPrefactor::default(),
            dump_expression: None,
            gamma_algebra: GammaAlgebraMode::Symbolic,
            parse_mode: NumeratorParseMode::Polynomial,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq)]
pub enum GammaAlgebraMode {
    Symbolic,
    Concrete,
}

pub type AtomStructure = SmartShadowStructure<SerializableSymbol, Vec<SerializableAtom>>;

pub struct RepeatingIterator<T> {
    elements: Vec<T>,
    positions: std::vec::IntoIter<usize>,
}

pub enum RepeatingIteratorTensorOrScalar<T: HasStructure> {
    Tensors(RepeatingIterator<T>),
    Scalars(RepeatingIterator<T::Scalar>),
}

impl<T> RepeatingIterator<T> {
    pub(crate) fn new(positions: Vec<usize>, elements: Vec<T>) -> Self {
        RepeatingIterator {
            elements,
            positions: positions.into_iter(),
        }
    }

    pub(crate) fn new_not_repeating(elements: Vec<T>) -> Self {
        let positions: Vec<usize> = (0..elements.len()).collect();
        RepeatingIterator {
            elements,
            positions: positions.into_iter(),
        }
    }
}

// #[test]
// fn rep_iter(){
//     let mut r= RepeatingIterator::new_not_repeating(vec![1,2,3,4,5]);
//     while let Some(s) =r.next()  {
//         println!("{}",s);
//     }
// }

impl<T: HasStructure> From<(TensorSet<T>, Vec<usize>)> for RepeatingIteratorTensorOrScalar<T> {
    fn from(value: (TensorSet<T>, Vec<usize>)) -> Self {
        match value.0 {
            TensorSet::Tensors(t) => {
                RepeatingIteratorTensorOrScalar::Tensors(RepeatingIterator::new(value.1, t))
            }
            TensorSet::Scalars(s) => {
                RepeatingIteratorTensorOrScalar::Scalars(RepeatingIterator::new(value.1, s))
            }
        }
    }
}

impl<T: HasStructure> From<TensorSet<T>> for RepeatingIteratorTensorOrScalar<T> {
    fn from(value: TensorSet<T>) -> Self {
        match value {
            TensorSet::Tensors(t) => {
                RepeatingIteratorTensorOrScalar::Tensors(RepeatingIterator::new_not_repeating(t))
            }
            TensorSet::Scalars(s) => {
                RepeatingIteratorTensorOrScalar::Scalars(RepeatingIterator::new_not_repeating(s))
            }
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
pub struct Numerator<State> {
    pub state: State,
}

impl<S: NumeratorState> Numerator<S> {
    pub(crate) fn export(&self) -> String {
        self.state.export()
    }

    pub(crate) fn forget_type(self) -> Numerator<PythonState> {
        Numerator {
            state: self.state.forget_type(),
        }
    }

    pub(crate) fn update_model(&mut self, model: &Model) -> Result<()> {
        self.state.update_model(model)
    }

    fn add_consts_to_fn_map(fn_map: &mut FunctionMap) {
        fn_map.add_constant(
            parse!("Nc"),
            symbolica::domains::float::Complex::from(Rational::from(3)),
        );

        fn_map.add_constant(
            parse!("TR"),
            symbolica::domains::float::Complex::from(Rational::from((1, 2))),
        );

        fn_map.add_constant(
            parse!("pi"),
            symbolica::domains::float::Complex::from(
                Rational::try_from(std::f64::consts::PI).unwrap(),
            ),
        );
    }
}

impl<S: GetSingleAtom> Numerator<S> {
    pub(crate) fn get_single_atom(&self) -> Result<Atom, NumeratorStateError> {
        self.state.get_single_atom()
    }
}

pub trait TypedNumeratorState:
    NumeratorState + TryFrom<PythonState, Error: std::error::Error + Send + Sync + 'static>
{
    fn apply<F, S: TypedNumeratorState>(
        state: &mut Numerator<PythonState>,
        f: F,
    ) -> Result<(), NumeratorStateError>
    where
        F: FnMut(Numerator<Self>) -> Numerator<S>;
}

impl Numerator<PythonState> {
    pub(crate) fn try_from<S: TypedNumeratorState>(self) -> Result<Numerator<S>, Report> {
        Ok(Numerator {
            state: self.state.try_into()?,
        })
    }

    pub(crate) fn apply<F, S: TypedNumeratorState, T: TypedNumeratorState>(
        &mut self,
        f: F,
    ) -> Result<(), NumeratorStateError>
    where
        F: FnMut(Numerator<S>) -> Numerator<T>,
    {
        S::apply(self, f)
    }
}
pub trait NumeratorState:
    Clone + Debug + Encode + for<'a> Decode<GammaLoopContextContainer<'a>>
{
    fn export(&self) -> String;

    fn forget_type(self) -> PythonState;

    fn update_model(&mut self, model: &Model) -> Result<()>;
    // fn try_from(state: PythonState) -> Result<Self>;
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
pub struct UnInit;

impl Default for UnInit {
    fn default() -> Self {
        let _ = *UFO;
        let _ = *ETS;
        UnInit
    }
}

impl TryFrom<PythonState> for UnInit {
    type Error = NumeratorStateError;

    fn try_from(value: PythonState) -> std::result::Result<Self, Self::Error> {
        match value {
            PythonState::UnInit(s) => {
                if let Some(s) = s {
                    Ok(s)
                } else {
                    Err(NumeratorStateError::NoneVariant)
                }
            }
            _ => Err(NumeratorStateError::NotUnit),
        }
    }
}

impl NumeratorState for UnInit {
    fn export(&self) -> String {
        "Uninitialized".to_string()
    }

    fn forget_type(self) -> PythonState {
        PythonState::UnInit(Some(self))
    }

    fn update_model(&mut self, _model: &Model) -> Result<()> {
        Err(eyre!("Uninitialized, nothing to update"))
    }
}

impl TypedNumeratorState for UnInit {
    fn apply<F, S: TypedNumeratorState>(
        num: &mut Numerator<PythonState>,
        mut f: F,
    ) -> Result<(), NumeratorStateError>
    where
        F: FnMut(Numerator<Self>) -> Numerator<S>,
    {
        if let PythonState::UnInit(s) = &mut num.state {
            if let Some(s) = s.take() {
                *num = f(Numerator { state: s }).forget_type();
                return Ok(());
            } else {
                return Err(NumeratorStateError::NoneVariant);
            }
        }
        Err(NumeratorStateError::NotUnit)
    }
}

#[allow(dead_code)]
#[derive(JsonSchema)]
struct _GlobalPrefactorAny {
    projector: serde_json::Value,
    num: serde_json::Value,
}

#[derive(Debug, Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode, PartialEq)]
#[trait_decode(trait = crate::GammaLoopContext)]
pub struct GlobalPrefactor {
    pub projector: Atom,
    pub num: Atom,
}

impl JsonSchema for GlobalPrefactor {
    fn json_schema(generated: &mut schemars::SchemaGenerator) -> schemars::Schema {
        generated.subschema_for::<_GlobalPrefactorAny>()
    }
    fn schema_name() -> std::borrow::Cow<'static, str> {
        "GlobalPrefactor".into()
    }
}

impl GlobalPrefactor {
    pub fn polarizations(&self) -> Vec<(PolDef, Atom)> {
        let mut pols = Vec::new();
        let full_prefactor = &self.projector * &self.num;
        let pat = function!(GS.epsilon, W_.e_, W_.i_).to_pattern();
        for m in full_prefactor.pattern_match(&pat, None, None) {
            let Some(e) = m.get(&W_.e_) else {
                continue;
            };
            let Ok(e) = i64::try_from(e) else {
                continue;
            };

            pols.push((
                PolDef {
                    pol_type: PolType::Epsilon,
                    eid: EdgeIndex(e as usize),
                    hel: None,
                },
                pat.replace_wildcards(&m),
            ));
        }

        let pat = function!(GS.epsilonbar, W_.e_, W_.i_).to_pattern();
        for m in full_prefactor.pattern_match(&pat, None, None) {
            let Some(e) = m.get(&W_.e_) else {
                continue;
            };
            let Ok(e) = i64::try_from(e) else {
                continue;
            };

            pols.push((
                PolDef {
                    pol_type: PolType::EpsilonBar,
                    eid: EdgeIndex(e as usize),
                    hel: None,
                },
                pat.replace_wildcards(&m),
            ));
        }

        let pat = function!(GS.u, W_.e_, W_.i_).to_pattern();
        for m in full_prefactor.pattern_match(&pat, None, None) {
            let Some(e) = m.get(&W_.e_) else {
                continue;
            };
            let Ok(e) = i64::try_from(e) else {
                continue;
            };

            pols.push((
                PolDef {
                    pol_type: PolType::U,
                    eid: EdgeIndex(e as usize),
                    hel: None,
                },
                pat.replace_wildcards(&m),
            ));
        }

        let pat = function!(GS.v, W_.e_, W_.i_).to_pattern();
        for m in full_prefactor.pattern_match(&pat, None, None) {
            let Some(e) = m.get(&W_.e_) else {
                continue;
            };
            let Ok(e) = i64::try_from(e) else {
                continue;
            };

            pols.push((
                PolDef {
                    pol_type: PolType::V,
                    eid: EdgeIndex(e as usize),
                    hel: None,
                },
                pat.replace_wildcards(&m),
            ));
        }

        let pat = function!(GS.ubar, W_.e_, W_.i_).to_pattern();
        for m in full_prefactor.pattern_match(&pat, None, None) {
            let Some(e) = m.get(&W_.e_) else {
                continue;
            };
            let Ok(e) = i64::try_from(e) else {
                continue;
            };

            pols.push((
                PolDef {
                    pol_type: PolType::UBar,
                    eid: EdgeIndex(e as usize),
                    hel: None,
                },
                pat.replace_wildcards(&m),
            ));
        }

        let pat = function!(GS.vbar, W_.e_, W_.i_).to_pattern();
        for m in full_prefactor.pattern_match(&pat, None, None) {
            let Some(e) = m.get(&W_.e_) else {
                continue;
            };
            let Ok(e) = i64::try_from(e) else {
                continue;
            };

            pols.push((
                PolDef {
                    pol_type: PolType::VBar,
                    eid: EdgeIndex(e as usize),
                    hel: None,
                },
                pat.replace_wildcards(&m),
            ));
        }

        pols.sort_by(|a, b| a.0.cmp(&b.0));

        pols
    }
}

impl Default for GlobalPrefactor {
    fn default() -> Self {
        GlobalPrefactor {
            projector: Atom::num(1),
            num: Atom::num(1),
        }
    }
}

impl Serialize for GlobalPrefactor {
    fn serialize<S>(&self, serializer: S) -> std::result::Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let mut state = serializer.serialize_struct("GlobalPrefactor", 2)?;
        state.serialize_field("projector", &self.projector.to_string())?;
        state.serialize_field("num", &self.num.to_string())?;

        state.end()
    }
}

impl<'de> Deserialize<'de> for GlobalPrefactor {
    fn deserialize<D>(deserializer: D) -> std::result::Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        #[derive(Deserialize)]
        struct GlobalPrefactorHelper {
            num: String,
            projector: String,
        }
        let helper = GlobalPrefactorHelper::deserialize(deserializer)?;
        Ok(GlobalPrefactor {
            projector: parse!(&helper.projector),
            num: parse!(&helper.num),
        })
    }
}

pub trait NumeratorNode {}

pub mod graph;
pub mod uninit;
#[allow(clippy::default_constructed_unit_structs)]
impl Default for Numerator<UnInit> {
    fn default() -> Self {
        Numerator {
            state: UnInit::default(),
        }
    }
}

#[derive(Debug, Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = symbolica::state::HasStateMap)]
pub struct SymbolicExpression<State> {
    // pub colorless: ParamTensor<OrderedStructure<Euclidean, Aind>>,
    // pub color: ParamTensor<OrderedStructure<Euclidean, Aind>>,
    pub expr: Atom,
    pub state: State,
}

pub trait GetSingleAtom {
    fn get_single_atom(&self) -> Result<Atom, NumeratorStateError>;
}

pub trait UnexpandedNumerator: NumeratorState + GetSingleAtom {
    // fn expr(&self) -> Result<SerializableAtom, NumeratorStateError>;

    fn map_color(self, f: impl Fn(Atom) -> Atom) -> Self;

    fn map_color_mut(&mut self, f: impl FnMut(&mut Atom));

    fn map_colorless(self, f: impl Fn(Atom) -> Atom) -> Self;
}

impl<E: ExpressionState> GetSingleAtom for SymbolicExpression<E> {
    fn get_single_atom(&self) -> Result<Atom, NumeratorStateError> {
        Ok(self.expr.clone())
    }
}

impl<E: ExpressionState> SymbolicExpression<E> {
    // // #[allow(dead_code)]
    // fn map_color(self, f: impl Fn(Atom) -> Atom) -> Self {
    //     SymbolicExpression {
    //         colorless: self.colorless,
    //         color: self.color.map_data_self(f),
    //         state: self.state,
    //     }
    // }

    // // #[allow(dead_code)]
    // fn map_color_mut(&mut self, f: impl FnMut(&mut Atom)) {
    //     self.color.map_data_mut(f);
    // }
    // // #[allow(dead_code)]
    // fn map_colorless(self, f: impl Fn(Atom) -> Atom) -> Self {
    //     SymbolicExpression {
    //         colorless: self.colorless.map_data_self(f),
    //         color: self.color,
    //         state: self.state,
    //     }
    // }
}

impl<E: ExpressionState> SymbolicExpression<E> {
    pub(crate) fn new(expression: SerializableAtom) -> Self {
        E::new(expression)
    }
    pub(crate) fn new_color(expression: SerializableAtom) -> Self {
        E::new_color(expression)
    }
}

pub trait ExpressionState:
    Serialize
    + Clone
    + DeserializeOwned
    + Debug
    + Encode
    + for<'a> Decode<GammaLoopContextContainer<'a>>
    + Default
{
    fn forget_type(data: SymbolicExpression<Self>) -> PythonState;

    fn new(expression: SerializableAtom) -> SymbolicExpression<Self> {
        SymbolicExpression {
            // colorless: ParamTensor::composite(DataTensor::new_scalar(expression.0)),
            // color: ParamTensor::composite(DataTensor::new_scalar(Atom::num(1))),
            expr: expression.0,
            state: Self::default(),
        }
    }

    fn new_color(expression: SerializableAtom) -> SymbolicExpression<Self> {
        SymbolicExpression {
            // color: ParamTensor::composite(DataTensor::new_scalar(expression.0)),
            // colorless: ParamTensor::composite(DataTensor::new_scalar(Atom::num(1))),
            expr: expression.0,
            state: Self::default(),
        }
    }

    fn get_expression(
        num: &mut PythonState,
    ) -> Result<SymbolicExpression<Self>, NumeratorStateError>;
}

impl<E: ExpressionState> TypedNumeratorState for SymbolicExpression<E> {
    fn apply<F, S: TypedNumeratorState>(
        num: &mut Numerator<PythonState>,
        mut f: F,
    ) -> Result<(), NumeratorStateError>
    where
        F: FnMut(Numerator<Self>) -> Numerator<S>,
    {
        let s = Self::try_from(&mut num.state)?;
        *num = f(Numerator { state: s }).forget_type();
        Ok(())
    }
}

impl<E: ExpressionState> TryFrom<&mut PythonState> for SymbolicExpression<E> {
    type Error = NumeratorStateError;

    fn try_from(value: &mut PythonState) -> std::result::Result<Self, Self::Error> {
        let a = E::get_expression(value)?;
        Ok(a)
    }
}

impl<E: ExpressionState> TryFrom<PythonState> for SymbolicExpression<E> {
    type Error = NumeratorStateError;

    fn try_from(mut value: PythonState) -> std::result::Result<Self, Self::Error> {
        let a = E::get_expression(&mut value)?;
        Ok(a)
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Encode, Decode, Default)]
pub struct Local {}
pub type AppliedFeynmanRule = SymbolicExpression<Local>;

impl ExpressionState for Local {
    fn forget_type(data: SymbolicExpression<Self>) -> PythonState {
        PythonState::AppliedFeynmanRule(Some(data))
    }

    fn get_expression(
        num: &mut PythonState,
    ) -> Result<SymbolicExpression<Self>, NumeratorStateError> {
        if let PythonState::AppliedFeynmanRule(s) = num {
            if let Some(s) = s.take() {
                Ok(s)
            } else {
                Err(NumeratorStateError::NoneVariant)
            }
        } else {
            Err(NumeratorStateError::NotAppliedFeynmanRule)
        }
    }
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, Encode, Decode, Default)]
pub struct NonLocal {}
pub type Global = SymbolicExpression<NonLocal>;

impl ExpressionState for NonLocal {
    fn forget_type(data: SymbolicExpression<Self>) -> PythonState {
        PythonState::Global(Some(data))
    }

    fn get_expression(
        num: &mut PythonState,
    ) -> Result<SymbolicExpression<Self>, NumeratorStateError> {
        if let PythonState::Global(s) = num {
            if let Some(s) = s.take() {
                Ok(s)
            } else {
                Err(NumeratorStateError::NoneVariant)
            }
        } else {
            Err(NumeratorStateError::NotGlobal)
        }
    }
}

impl Numerator<Global> {
    #[instrument(skip(self))]
    pub(crate) fn color_simplify(self) -> Numerator<ColorSimplified> {
        // debug!("Color simplifying global numerator");
        // let mut fully_simplified = true;

        let state = ColorSimplified {
            expr: self.state.expr.simplify_color(),
            state: Color::Fully,
        };
        // debug!("Color simplified numerator:{}", state.expr);
        Numerator { state }
    }
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, Encode, Decode, Default)]
pub enum Color {
    #[default]
    Fully,
    Partially,
}
pub type ColorSimplified = SymbolicExpression<Color>;

impl ExpressionState for Color {
    fn forget_type(data: SymbolicExpression<Self>) -> PythonState {
        PythonState::ColorSimplified(Some(data))
    }

    fn get_expression(
        num: &mut PythonState,
    ) -> Result<SymbolicExpression<Self>, NumeratorStateError> {
        if let PythonState::ColorSimplified(s) = num {
            if let Some(s) = s.take() {
                Ok(s)
            } else {
                Err(NumeratorStateError::NoneVariant)
            }
        } else {
            Err(NumeratorStateError::NotColorSimplified)
        }
    }
}

// #[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, Default)]
// pub struct Projected {}
// pub type ColorProjected = SymbolicExpression<Projected>;

// impl ExpressionState for Projected {
//     fn forget_type(data: SymbolicExpression<Self>) -> PythonState {
//         PythonState::ColorProjected(Some(data))
//     }

//     fn get_expression(
//         num: &mut PythonState,
//     ) -> Result<SymbolicExpression<Self>, NumeratorStateError> {
//         if let PythonState::ColorProjected(s) = num {
//             if let Some(s) = s.take() {
//                 Ok(s)
//             } else {
//                 Err(NumeratorStateError::NoneVariant)
//             }
//         } else {
//             Err(NumeratorStateError::NotColorProjected)
//         }
//     }
// }

#[derive(Debug, Copy, Clone, Serialize, Deserialize, Encode, Decode, Default)]
pub struct Gamma {}
pub type GammaSimplified = SymbolicExpression<Gamma>;

impl ExpressionState for Gamma {
    fn forget_type(data: SymbolicExpression<Self>) -> PythonState {
        PythonState::GammaSimplified(Some(data))
    }

    fn get_expression(
        num: &mut PythonState,
    ) -> Result<SymbolicExpression<Self>, NumeratorStateError> {
        if let PythonState::GammaSimplified(s) = num {
            if let Some(s) = s.take() {
                Ok(s)
            } else {
                Err(NumeratorStateError::NoneVariant)
            }
        } else {
            Err(NumeratorStateError::NotGammaSymplified)
        }
    }
}

impl<State: ExpressionState> NumeratorState for SymbolicExpression<State> {
    fn export(&self) -> String {
        self.get_single_atom().unwrap().to_canonical_string()
    }

    fn forget_type(self) -> PythonState {
        State::forget_type(self)
    }

    fn update_model(&mut self, _model: &Model) -> Result<()> {
        Err(eyre!("Only an expression, nothing to update"))
    }
}

impl<T: Copy + Default> Numerator<SymbolicExpression<T>> {
    pub(crate) fn canonize_lorentz(&self) -> Result<Self, String> {
        let pats: Vec<LibraryRep> = vec![Minkowski {}.into(), Bispinor {}.into()];

        let mut indices_map = AHashMap::new();

        for p in &pats {
            for a in self.state.expr.pattern_match(
                &p.to_symbolic([W_.x_, W_.y_]).to_pattern(),
                None,
                None,
            ) {
                indices_map.insert(
                    p.to_symbolic([a[&W_.x_].clone(), a[&W_.y_].clone()]),
                    p.to_symbolic([a[&W_.x_].clone()]),
                );
            }
        }

        let sorted = indices_map.into_iter().sorted().collect::<Vec<_>>();

        let expr = self.state.expr.canonize_tensors(sorted)?.canonical_form;

        Ok(Self {
            state: SymbolicExpression {
                expr,
                state: T::default(),
            },
        })
    }

    // pub(crate) fn canonize_color(&self) -> Result<Self, String> {

    //     let pats: Vec<_> = vec![ColorAdjoint{}];
    //     let dualizablepats: Vec<_> = vec![
    //         ColorFundamental::selfless_symbol(),
    //         ColorSextet::selfless_symbol(),
    //     ];

    //     let mut color = self.state.expr.replace(&function!(symbol!(DOWNIND), GS.x__).to_pattern())
    //                 .with(Atom::var(GS.x__).to_pattern())
    //     ;

    //     let mut indices_map = AHashMap::new();

    //     color.iter_flat().for_each(|(_, v)| {
    //         for p in pats.iter().chain(&dualizablepats) {
    //             for a in
    //                 v.0.pattern_match(&function!(*p, GS.x_, GS.y_).to_pattern(), None, None)
    //             {
    //                 indices_map.insert(
    //                     function!(*p, a[&GS.x_], a[&GS.y_]),
    //                     function!(*p, a[&GS.x_]),
    //                 );
    //             }
    //         }
    //     });

    //     let sorted = indices_map.into_iter().sorted().collect::<Vec<_>>();
    //     // println!(
    //     //     "indices sorted [{}]",
    //     //     sorted
    //     //         .iter()
    //     //         .map(|(a, b)| format!(
    //     //             "(Atom::parse(\"{}\").unwrap(),Atom::parse(\"{}\").unwrap())",
    //     //             a, b
    //     //         ))
    //     //         .collect::<Vec<_>>()
    //     //         .join(", ")
    //     // );

    //     color = color.map_data_ref_result(|a| a.0.canonize_tensors(&sorted).map(|a| a.into()))?;

    //     let colorless = self.state.colorless.clone();
    //     Ok(Self {
    //         state: SymbolicExpression {
    //             colorless,
    //             color,
    //             state: T::default(),
    //         },
    //     }));
    //     todo!()
    // }
}

impl Numerator<AppliedFeynmanRule> {
    pub(crate) fn to_d_dim<'a>(mut self, dim: impl Into<AtomOrView<'a>>) -> Self {
        self.state.expr = self.state.expr.map_mink_dim(dim);
        self
    }

    pub(crate) fn color_simplify(self) -> Numerator<ColorSimplified> {
        debug!("Color simplifying global numerator");
        // let mut fully_simplified = true;

        let state = ColorSimplified {
            expr: self.state.expr.simplify_color(),
            state: Color::Fully,
        };
        debug!("Color simplified numerator:{}", state.expr);
        Numerator { state }
    }
}

// #[derive(Debug, Error)]
// pub enum ColorError {
//     #[error("Not fully simplified: {0}")]
//     NotFully(SerializableAtom),
// }

impl Numerator<ColorSimplified> {
    pub(crate) fn gamma_simplify(self) -> Numerator<GammaSimplified> {
        debug!("Gamma simplifying color symplified numerator");

        Numerator {
            state: GammaSimplified {
                expr: self.state.expr.simplify_gamma(),
                state: Default::default(),
            },
        }
    }
}

pub type Gloopoly =
    symbolica::poly::polynomial::MultivariatePolynomial<symbolica::domains::atom::AtomField, u8>;

#[derive(Debug, Clone)]
// #[trait_decode(trait = symbolica::state::HasStateMap)]
pub struct PolySplit {
    pub colorless: DataTensor<Gloopoly, ShadowedStructure<Aind>>,
    pub var_map: Arc<Vec<PolyVariable>>,
    pub energies: Vec<usize>,
    pub color: ParamTensor<ShadowedStructure<Aind>>,
    pub colorsimplified: Color,
}

impl Encode for PolySplit {
    fn encode<E: bincode::enc::Encoder>(
        &self,
        _encoder: &mut E,
    ) -> std::result::Result<(), bincode::error::EncodeError> {
        todo!()
    }
}

impl PolySplit {
    pub(crate) fn from_color_out(
        _color_simplified: Numerator<ColorSimplified>,
    ) -> DataTensor<Atom> {
        disable!(

        let colorless_parsed = color_simplified
            .state
            .colorless
            .map_data(|a| {
                let mut net =
                    TensorNetwork::<MixedTensor<f64, AtomStructure>, SerializableAtom>::try_from(
                        a.as_view(),
                    )
                    .unwrap();
                net.contract().unwrap();
                net.to_fully_parametric()
                    .result_tensor_smart()
                    .unwrap()
                    .tensor
                    .map_structure(OrderedStructure::from)
            })
            .flatten(&Atom::num(0))
            .unwrap();

        colorless_parsed
        );
        todo!()
    }

    fn to_shadowed_poly_impl(
        poly: &Gloopoly,
        workspace: &Workspace,
        reps: Arc<Mutex<Vec<Atom>>>,
    ) -> Atom {
        if poly.is_zero() {
            return Atom::num(0);
        }

        let mut add = Atom::num(0);
        let coef = symbol!("coef");
        let shift = reps.as_ref().lock().unwrap().len();

        let mut mul_h;
        let mut var_h = workspace.new_atom();
        let mut num_h = workspace.new_atom();
        let mut pow_h = workspace.new_atom();

        for (i, monomial) in poly.into_iter().enumerate() {
            mul_h = Atom::num(1);
            for (var_id, &pow) in poly.variables.iter().zip(monomial.exponents) {
                if pow > 0 {
                    match var_id {
                        PolyVariable::Symbol(v) => {
                            var_h.to_var(*v);
                        }
                        PolyVariable::Temporary(_) => {
                            unreachable!("Temporary variable in expression")
                        }
                        PolyVariable::Function(_, a) | PolyVariable::Power(a) => {
                            var_h.set_from_view(&a.as_view());
                        }
                    }

                    if pow > 0 {
                        num_h.to_num((pow as i64).into());
                        pow_h.to_pow(var_h.as_view(), num_h.as_view());
                        mul_h = mul_h * pow_h.as_view();
                    } else {
                        mul_h = mul_h * var_h.as_view();
                    }
                }
            }

            reps.lock()
                .as_mut()
                .unwrap()
                .push(monomial.coefficient.clone());

            mul_h = mul_h * function!(coef, Atom::num((i + shift) as i64));
            add = add + mul_h.as_view();
        }

        add
    }
    fn shadow_poly(poly: Gloopoly, reps: Arc<Mutex<Vec<Atom>>>) -> Atom {
        Workspace::get_local().with(|ws| Self::to_shadowed_poly_impl(&poly, ws, reps))
    }

    pub(crate) fn optimize(self) -> PolyContracted {
        let reps = Arc::new(Mutex::new(Vec::new()));

        let colorless = self
            .colorless
            .map_data(|a| Self::shadow_poly(a, reps.clone()));

        let out = match self.colorsimplified {
            Color::Fully => ParamTensor::composite(colorless)
                .contract(&self.color)
                .unwrap(),
            Color::Partially => {
                warn!(
                    "Not fully color-simplified, taking the colorless part associated to {}",
                    self.color.get_owned_linear(FlatIndex::from(0)).unwrap()
                );
                ParamTensor::new_scalar(colorless.get_owned_linear(FlatIndex::from(0)).unwrap())
            }
        };

        let reps = Arc::try_unwrap(reps).unwrap().into_inner().unwrap();

        PolyContracted {
            tensor: out,
            coef_map: reps,
        }
    }
}

impl Numerator<PolySplit> {
    pub(crate) fn contract(self) -> Result<Numerator<PolyContracted>> {
        match self.validate_squared_energies_impl() {
            Err(_) => {
                debug!("Trying to contract polynomial");
                let state = self.state.optimize();
                debug!("PolyContracted: {}", state.tensor);
                Ok(Numerator { state })
            }
            Ok(r) => Err(eyre!("has higher powers here: {}", r)),
        }
    }

    fn validate_squared_energies_impl(
        &self,
    ) -> Result<DataTensor<ExpandedIndex, ShadowedStructure<Aind>>, ()> {
        self.state.colorless.map_data_ref_result(|p| {
            let mut square: Result<Vec<usize>, ()> = Err(());
            for (i, &e) in p.exponents.iter().enumerate() {
                if e > 1 {
                    if let Ok(sq) = &mut square {
                        sq.push(i);
                    } else {
                        square = Ok(vec![i]);
                    }
                }
            }
            square.map(ExpandedIndex::from)
        })
    }

    pub(crate) fn validate_squared_energies(&self) -> bool {
        self.validate_squared_energies_impl().is_err()
    }
}

#[derive(Debug, Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = symbolica::state::HasStateMap)]
pub struct PolyContracted {
    pub tensor: ParamTensor<ShadowedStructure<Aind>>,
    pub coef_map: Vec<Atom>,
}

impl Numerator<PolyContracted> {
    pub(crate) fn to_contracted(self) -> Numerator<Contracted> {
        let coefs: Vec<_> = (0..self.state.coef_map.len())
            .map(|i| function!(GS.coeff, Atom::num(i as i64)).to_pattern())
            .collect();

        let coefs_reps: Vec<_> = self.state.coef_map.iter().map(|a| a.to_pattern()).collect();

        let reps: Vec<_> = coefs
            .into_iter()
            .zip(coefs_reps)
            .map(|(p, rhs)| Replacement::new(p, rhs))
            .collect();

        Numerator {
            state: Contracted {
                tensor: self.state.tensor.replace_multiple(&reps),
            },
        }
    }

    fn generate_fn_map(&self) -> FunctionMap {
        let mut fn_map = FunctionMap::new();

        for (v, k) in self.state.coef_map.clone().iter().enumerate() {
            fn_map
                .add_tagged_function(
                    GS.coeff,
                    vec![Atom::num(v as i64)],
                    format!("coef{v}"),
                    vec![],
                    k.clone(),
                )
                .unwrap();
        }

        Numerator::<Contracted>::add_consts_to_fn_map(&mut fn_map);
        fn_map
    }
}

impl NumeratorState for PolyContracted {
    fn export(&self) -> String {
        self.tensor.to_string()
    }
    fn forget_type(self) -> PythonState {
        PythonState::PolyContracted(Some(self))
    }

    fn update_model(&mut self, _model: &Model) -> Result<()> {
        Err(eyre!(
            "Only applied feynman rule, simplified color, gamma, parsed into network and contracted, nothing to update"
        ))
    }
}

impl PolyContracted {}

impl GammaSimplified {
    pub(crate) fn parse(self) -> Network {
        let lib = DummyLibrary::<(), _>::new();
        let net = StandardTensorNet::try_from_view(
            self.get_single_atom().unwrap().as_view(),
            &lib,
            &ParseSettings::default(),
        )
        .unwrap();

        // println!("net scalar{}", net.scalar.as_ref().unwrap());
        Network { net }
    }

    // pub(crate) fn parse_only_colorless(self) -> Network {
    //     let lib = DummyLibrary::<(), _>::new();
    //     let net = StandardTensorNet::try_from_view(
    //         self.colorless
    //             .clone()
    //             .scalar()
    //             .ok_or(NumeratorStateError::Any(eyre!("not a scalar")))
    //             .unwrap()
    //             .as_view(),
    //         &lib,
    //     )
    //     .unwrap();

    //     // println!("net scalar{}", net.scalar.as_ref().unwrap());
    //     Network { net }
    // }
}

impl Numerator<GammaSimplified> {
    pub(crate) fn parse(self) -> Result<Numerator<Network>> {
        // debug!("GammaSymplified numerator: {}", self.export());
        debug!("Parsing gamma simplified numerator into tensor network");
        Ok(Numerator {
            state: Network {
                net: ParsingNet::try_from_view(
                    self.state.expr.as_view(),
                    TENSORLIB.read().unwrap().deref(),
                    &ParseSettings::default(),
                )?,
            },
        })
    }

    // pub(crate) fn parse_only_colorless(self) -> Numerator<Network> {
    //     debug!("Parsing only colorless gamma simplified numerator into tensor network");
    //     Numerator {
    //         state: self.state.parse_only_colorless(),
    //     }
    // }
}

pub type ParsingNet = spenso::network::Network<
    NetworkStore<MixedTensor<F<f64>, ShadowedStructure<Aind>>, Atom>,
    ExplicitKey<Aind>,
    Symbol,
    Aind,
>;

pub type ParamParsingNet = spenso::network::Network<
    NetworkStore<ParamTensor<ShadowedStructure<Aind>>, Atom>,
    ExplicitKey<Aind>,
    Symbol,
    Aind,
>;

pub type IntParsingNet = spenso::network::Network<
    NetworkStore<MixedTensor<i64, ShadowedStructure<Aind>>, Atom>,
    ExplicitKey<Aind>,
    Symbol,
    Aind,
>;

#[derive(Debug, Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = symbolica::state::HasStateMap)]
pub struct Network {
    // #[bincode(with_serde
    pub net: ParsingNet,
}

impl TryFrom<PythonState> for Network {
    type Error = NumeratorStateError;

    fn try_from(value: PythonState) -> std::result::Result<Self, Self::Error> {
        match value {
            PythonState::Network(s) => {
                if let Some(s) = s {
                    Ok(s)
                } else {
                    Err(NumeratorStateError::NoneVariant)
                }
            }
            _ => Err(NumeratorStateError::NotNetwork),
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct ContractionSettings {
    n_steps: Option<usize>,
    mode: ExecutionMode,
}

impl From<()> for ContractionSettings {
    fn from(_: ()) -> Self {
        Self {
            n_steps: None,
            mode: ExecutionMode::All,
        }
    }
}

impl Default for ContractionSettings {
    fn default() -> Self {
        Self {
            n_steps: None,
            mode: ExecutionMode::All,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub enum ExecutionMode {
    Single,
    Scalar,
    All,
}

pub type StandardTensorNet = spenso::network::Network<
    NetworkStore<MixedTensor<F<f64>, ShadowedStructure<Aind>>, Atom>,
    ExplicitKey<Aind>,
    Symbol,
    Aind,
>;

impl Network {
    pub(crate) fn parse_impl(expr: AtomView) -> Self {
        let lib = DummyLibrary::<(), _>::new();
        let net = StandardTensorNet::try_from_view(expr, &lib, &ParseSettings::default()).unwrap();

        // println!("net scalar{}", net.scalar.as_ref().unwrap());
        Network { net }
    }

    pub(crate) fn contract(&mut self, settings: impl Into<ContractionSettings>) -> Result<()> {
        let lib = TENSORLIB.read().unwrap();
        let fnlib = FUN_LIB.deref();
        let settings = settings.into();
        debug!(
            "Contracting network:{} with settings: {:#?}",
            self.net.dot(),
            settings
        );
        if let Some(n) = settings.n_steps {
            for _ in 0..n {
                match settings.mode {
                    ExecutionMode::All => self
                        .net
                        .execute::<Steps<1>, SmallestDegree, _, _, _>(lib.deref(), fnlib)?,
                    ExecutionMode::Scalar => self
                        .net
                        .execute::<Steps<1>, ContractScalars, _, _, _>(lib.deref(), fnlib)?,
                    ExecutionMode::Single => self
                        .net
                        .execute::<Steps<1>, SingleSmallestDegree<false>, _, _, _>(
                            lib.deref(),
                            fnlib,
                        )?,
                }
            }
        } else {
            match settings.mode {
                ExecutionMode::All => {
                    self.net
                        .execute::<Sequential, SmallestDegree, _, _, _>(lib.deref(), fnlib)?;
                }
                ExecutionMode::Scalar => {
                    self.net
                        .execute::<Sequential, ContractScalars, _, _, _>(lib.deref(), fnlib)?;
                }
                ExecutionMode::Single => {
                    self.net
                        .execute::<Sequential, SingleSmallestDegree<false>, _, _, _>(
                            lib.deref(),
                            fnlib,
                        )?;
                }
            }
        }
        Ok(())
    }
}

impl NumeratorState for Network {
    fn export(&self) -> String {
        " self.expression.to_string()".to_string()
    }

    fn forget_type(self) -> PythonState {
        PythonState::Network(Some(self))
    }

    fn update_model(&mut self, _model: &Model) -> Result<()> {
        Err(eyre!(
            "Only applied feynman rule, simplified color, gamma and parsed into network, nothing to update"
        ))
    }
}

impl TypedNumeratorState for Network {
    fn apply<F, S: TypedNumeratorState>(
        num: &mut Numerator<PythonState>,
        mut f: F,
    ) -> Result<(), NumeratorStateError>
    where
        F: FnMut(Numerator<Self>) -> Numerator<S>,
    {
        if let PythonState::Network(s) = &mut num.state {
            if let Some(s) = s.take() {
                *num = f(Numerator { state: s }).forget_type();
                return Ok(());
            } else {
                return Err(NumeratorStateError::NoneVariant);
            }
        }
        Err(NumeratorStateError::NotNetwork)
    }
}

#[derive(Debug, Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = symbolica::state::HasStateMap)]
pub struct Contracted {
    pub tensor: ParamTensor<ShadowedStructure<Aind>>,
}

impl TryFrom<PythonState> for Contracted {
    type Error = NumeratorStateError;

    fn try_from(value: PythonState) -> std::result::Result<Self, Self::Error> {
        match value {
            PythonState::Contracted(s) => {
                if let Some(s) = s {
                    Ok(s)
                } else {
                    Err(NumeratorStateError::NoneVariant)
                }
            }
            _ => Err(NumeratorStateError::NotContracted),
        }
    }
}

impl Contracted {
    pub(crate) fn one() -> Self {
        Contracted {
            tensor: ParamTensor::new_scalar(Atom::one()),
        }
    }

    pub(crate) fn zero() -> Self {
        Contracted {
            tensor: ParamTensor::new_scalar(Atom::Zero),
        }
    }

    pub(crate) fn generate_kinematic_params_impl(
        n_edges: usize,
        pol_data: Vec<(String, i64, usize)>,
    ) -> Vec<Atom> {
        fn atoms_for_pol(name: String, num: i64, size: usize) -> Vec<Atom> {
            let mut data = vec![];
            for index in 0..size {
                let e = FunctionBuilder::new(symbol!(&name));
                data.push(
                    e.add_arg(Atom::num(num))
                        .add_arg(parse!(&format!("cind({})", index)))
                        .finish(),
                );
            }
            data
        }
        let mut params: Vec<Atom> = vec![];

        let mut pols = Vec::new();

        for i in 0..n_edges {
            let named_structure: NamedStructure<String> =
                NamedStructure::from_iter([Lorentz {}.new_slot(4, i)], "Q".into(), Some(i))
                    .structure;
            params.extend(
                named_structure
                    .to_shell()
                    .expanded_shadow()
                    .unwrap()
                    .data
                    .clone(),
            );
        }

        for (name, num, size) in pol_data {
            pols.extend(atoms_for_pol(name, num, size));
        }

        params.extend(pols);
        params.push(Atom::i());

        params
    }
}

impl NumeratorState for Contracted {
    fn export(&self) -> String {
        self.tensor.to_string()
    }
    fn forget_type(self) -> PythonState {
        PythonState::Contracted(Some(self))
    }

    fn update_model(&mut self, _model: &Model) -> Result<()> {
        Err(eyre!(
            "Only applied feynman rule, simplified color, gamma, parsed into network and contracted, nothing to update"
        ))
    }
}

impl TypedNumeratorState for Contracted {
    fn apply<F, S: TypedNumeratorState>(
        num: &mut Numerator<PythonState>,
        mut f: F,
    ) -> Result<(), NumeratorStateError>
    where
        F: FnMut(Numerator<Self>) -> Numerator<S>,
    {
        if let PythonState::Contracted(s) = &mut num.state {
            if let Some(s) = s.take() {
                *num = f(Numerator { state: s }).forget_type();
                return Ok(());
            } else {
                return Err(NumeratorStateError::NoneVariant);
            }
        }
        Err(NumeratorStateError::NotContracted)
    }
}

impl Numerator<Contracted> {
    fn generate_fn_map(&self) -> FunctionMap {
        let mut map = FunctionMap::new();
        Numerator::<Contracted>::add_consts_to_fn_map(&mut map);
        map
    }
    // #[allow(clippy::too_many_arguments)]
    // pub(crate) fn generate_evaluators_from_params(
    //     self,
    //     n_edges: usize,
    //     name: &str,
    //     model_params_start: usize,
    //     params: &[Atom],
    //     double_param_values: Vec<Complex<F<f64>>>,
    //     quad_param_values: Vec<Complex<F<f128>>>,
    //     extra_info: &ExtraInfo,
    //     export_settings: &GenerationSettings,
    // ) -> Numerator<Evaluators> {
    //     let o = &export_settings.numerator_settings.eval_settings;
    //     let owned_fn_map = self.generate_fn_map();

    //     let inline_asm = export_settings.compile.inline_asm();
    //     let compile_options = export_settings
    //         .compile
    //         .to_symbolica_compile_options();

    //     let eval_settings = NumeratorEvaluatorSettings {
    //         options: o,
    //         inline_asm,
    //         compile_options,
    //     };

    //     let fn_map: FunctionMap = owned_fn_map;
    //     let single = self.state.evaluator(
    //         extra_info.path.clone(),
    //         name,
    //         &eval_settings,
    //         params,
    //         &fn_map,
    //     );

    //     match o {
    //         NumeratorEvaluatorOptions::Joint(_) => Numerator {
    //             state: Evaluators {
    //                 orientated: Some(single.orientated_joint_impl(
    //                     n_edges,
    //                     name,
    //                     params,
    //                     extra_info,
    //                     &eval_settings,
    //                     &fn_map,
    //                 )),
    //                 single,
    //                 choice: SingleOrCombined::Combined,
    //                 orientations: extra_info.orientations.clone(),
    //                 quad_param_values,
    //                 double_param_values,
    //                 model_params_start,
    //                 emr_len: n_edges,
    //                 fn_map,
    //             },
    //         },
    //         NumeratorEvaluatorOptions::Iterative(IterativeOptions {
    //             iterations,
    //             n_cores,
    //             verbose,
    //             ..
    //         }) => Numerator {
    //             state: Evaluators {
    //                 orientated: Some(single.orientated_iterative_impl(
    //                     n_edges,
    //                     name,
    //                     params,
    //                     extra_info,
    //                     &eval_settings,
    //                     *iterations,
    //                     *n_cores,
    //                     *verbose,
    //                     &fn_map,
    //                 )),
    //                 single,
    //                 choice: SingleOrCombined::Combined,
    //                 orientations: extra_info.orientations.clone(),
    //                 quad_param_values,
    //                 double_param_values,
    //                 model_params_start,
    //                 emr_len: n_edges,
    //                 fn_map,
    //             },
    //         },
    //         _ => Numerator {
    //             state: Evaluators {
    //                 orientated: None,
    //                 single,
    //                 choice: SingleOrCombined::Single,
    //                 orientations: extra_info.orientations.clone(),
    //                 quad_param_values,
    //                 double_param_values,
    //                 model_params_start,
    //                 emr_len: n_edges,
    //                 fn_map,
    //             },
    //         },
    //     }
    // }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq)]
pub struct IterativeOptions {
    pub eval_options: EvaluatorOptions,
    pub iterations: usize,
    pub n_cores: usize,
    pub verbose: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq)]
pub enum NumeratorParseMode {
    Polynomial,
    Direct,
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq)]
#[serde(tag = "type")]
pub enum NumeratorEvaluatorOptions {
    #[serde(rename = "Single")]
    Single(EvaluatorOptions),
    #[serde(rename = "Joint")]
    Joint(EvaluatorOptions),
    #[serde(rename = "Iterative")]
    Iterative(IterativeOptions),
}

#[derive(Clone)]
pub struct NumeratorEvaluatorSettings<'a> {
    pub options: &'a NumeratorEvaluatorOptions,
    pub inline_asm: InlineASM,
    pub compile_options: CompileOptions,
}

impl Default for NumeratorEvaluatorOptions {
    fn default() -> Self {
        NumeratorEvaluatorOptions::Single(EvaluatorOptions::default())
    }
}

impl NumeratorEvaluatorOptions {
    pub(crate) fn compile_options(&self) -> NumeratorCompileOptions {
        match self {
            NumeratorEvaluatorOptions::Single(options) => options.compile_options,
            NumeratorEvaluatorOptions::Joint(options) => options.compile_options,
            NumeratorEvaluatorOptions::Iterative(options) => options.eval_options.compile_options,
        }
    }

    pub(crate) fn cpe_rounds(&self) -> Option<usize> {
        match self {
            NumeratorEvaluatorOptions::Single(options) => options.cpe_rounds,
            NumeratorEvaluatorOptions::Joint(options) => options.cpe_rounds,
            NumeratorEvaluatorOptions::Iterative(options) => options.eval_options.cpe_rounds,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Copy, PartialEq, Eq, Hash, Encode, Decode)]
pub struct EvaluatorOptions {
    pub cpe_rounds: Option<usize>,
    pub compile_options: NumeratorCompileOptions,
}

impl Default for EvaluatorOptions {
    fn default() -> Self {
        EvaluatorOptions {
            cpe_rounds: Some(1),
            compile_options: NumeratorCompileOptions::Compiled,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Copy, PartialEq, Eq, Hash, Encode, Decode)]
#[serde(tag = "subtype")]
pub enum NumeratorCompileOptions {
    #[serde(rename = "Compiled")]
    Compiled,
    #[serde(rename = "NotCompiled")]
    NotCompiled,
}

impl NumeratorCompileOptions {
    pub(crate) fn compile(&self) -> bool {
        matches!(self, NumeratorCompileOptions::Compiled)
    }
}

#[derive(Error, Debug)]
pub enum NumeratorStateError {
    #[error("Not UnInit")]
    NotUnit,
    #[error("Not AppliedFeynmanRule")]
    NotAppliedFeynmanRule,
    #[error("Not ColorProjected")]
    NotColorProjected,
    #[error("Not Global")]
    NotGlobal,
    #[error("Not ColorSimplified")]
    NotColorSimplified,
    #[error("Not GammaSymplified")]
    NotGammaSymplified,
    #[error("Not Network")]
    NotNetwork,
    #[error("Not Contracted")]
    NotContracted,
    #[error("Not Evaluators")]
    NotEvaluators,
    #[error("None variant")]
    NoneVariant,
    #[error("Expanded")]
    Expanded,
    #[error("Any")]
    Any(#[from] eyre::Report),
}

#[derive(Debug, Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = symbolica::state::HasStateMap)]
#[allow(clippy::large_enum_variant)]
pub enum PythonState {
    UnInit(Option<UnInit>),
    Global(Option<Global>),
    AppliedFeynmanRule(Option<AppliedFeynmanRule>),
    ColorSimplified(Option<ColorSimplified>),
    // ColorProjected(Option<ColorProjected>),
    GammaSimplified(Option<GammaSimplified>),
    Network(Option<Network>),
    Contracted(Option<Contracted>),
    PolyContracted(Option<PolyContracted>),
}

impl Default for PythonState {
    fn default() -> Self {
        PythonState::UnInit(Some(UnInit))
    }
}

impl NumeratorState for PythonState {
    fn export(&self) -> String {
        match self {
            PythonState::UnInit(state) => {
                if let Some(s) = state {
                    s.export()
                } else {
                    "None".into()
                }
            }
            PythonState::Global(state) => {
                if let Some(s) = state {
                    s.export()
                } else {
                    "None".into()
                }
            }
            PythonState::AppliedFeynmanRule(state) => {
                if let Some(s) = state {
                    s.export()
                } else {
                    "None".into()
                }
            }
            PythonState::ColorSimplified(state) => {
                if let Some(s) = state {
                    s.export()
                } else {
                    "None".into()
                }
            }

            // PythonState::ColorProjected(state) => {
            //     if let Some(s) = state {
            //         s.export()
            //     } else {
            //         "None".into()
            //     }
            // }
            PythonState::GammaSimplified(state) => {
                if let Some(s) = state {
                    s.export()
                } else {
                    "None".into()
                }
            }
            PythonState::Network(state) => {
                if let Some(s) = state {
                    s.export()
                } else {
                    "None".into()
                }
            }
            PythonState::Contracted(state) => {
                if let Some(s) = state {
                    s.export()
                } else {
                    "None".into()
                }
            }
            PythonState::PolyContracted(state) => {
                if let Some(s) = state {
                    s.export()
                } else {
                    "None".into()
                }
            }
        }
    }

    fn forget_type(self) -> PythonState {
        self
    }

    fn update_model(&mut self, model: &Model) -> Result<()> {
        match self {
            PythonState::Global(state) => {
                if let Some(s) = state {
                    s.update_model(model)
                } else {
                    Err(NumeratorStateError::NoneVariant.into())
                }
            }
            PythonState::AppliedFeynmanRule(state) => {
                if let Some(s) = state {
                    s.update_model(model)
                } else {
                    Err(NumeratorStateError::NoneVariant.into())
                }
            }
            PythonState::ColorSimplified(state) => {
                if let Some(s) = state {
                    s.update_model(model)
                } else {
                    Err(NumeratorStateError::NoneVariant.into())
                }
            }
            PythonState::GammaSimplified(state) => {
                if let Some(s) = state {
                    s.update_model(model)
                } else {
                    Err(NumeratorStateError::NoneVariant.into())
                }
            }
            PythonState::Network(state) => {
                if let Some(s) = state {
                    s.update_model(model)
                } else {
                    Err(NumeratorStateError::NoneVariant.into())
                }
            }
            PythonState::Contracted(state) => {
                if let Some(s) = state {
                    s.update_model(model)
                } else {
                    Err(NumeratorStateError::NoneVariant.into())
                }
            }
            PythonState::PolyContracted(state) => {
                if let Some(s) = state {
                    s.update_model(model)
                } else {
                    Err(NumeratorStateError::NoneVariant.into())
                }
            }

            _ => Err(eyre!("No model to update")),
        }
    }
}

impl GetSingleAtom for PythonState {
    fn get_single_atom(&self) -> Result<Atom, NumeratorStateError> {
        match self {
            PythonState::Global(state) => {
                if let Some(s) = state {
                    s.get_single_atom()
                } else {
                    Err(NumeratorStateError::NoneVariant)
                }
            }
            PythonState::AppliedFeynmanRule(state) => {
                if let Some(s) = state {
                    s.get_single_atom()
                } else {
                    Err(NumeratorStateError::NoneVariant)
                }
            }
            PythonState::ColorSimplified(state) => {
                if let Some(s) = state {
                    s.get_single_atom()
                } else {
                    Err(NumeratorStateError::NoneVariant)
                }
            }
            PythonState::GammaSimplified(state) => {
                if let Some(s) = state {
                    s.get_single_atom()
                } else {
                    Err(NumeratorStateError::NoneVariant)
                }
            }
            _ => Err(NumeratorStateError::Expanded),
        }
    }
}

impl PythonState {}
// #[cfg(test)]
// mod tests;
