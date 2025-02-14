use crate::debug_info::DEBUG_LOGGER;
use crate::feyngen::FeynGenError;
use log::warn;
use std::fmt::Debug;
use std::fs;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};
use std::time::Instant;
// use crate::feyngen::dis::{DisEdge, DisVertex};
use crate::graph::{BareGraph, VertexInfo};
use crate::model::normalise_complex;
use crate::momentum::Polarization;
use crate::utils::{f128, GS};
use crate::{
    graph::EdgeType,
    model::Model,
    momentum::FourMomentum,
    utils::{FloatLike, F},
};
use crate::{ExportSettings, Settings};
use ahash::{AHashMap, HashSet};
use bincode::{Decode, Encode};
use color_eyre::{Report, Result};
use eyre::eyre;
use gat_lending_iterator::LendingIterator;
// use gxhash::GxBuildHasher;
use indexmap::IndexSet;
use itertools::Itertools;

use log::{debug, trace};

use ref_ops::RefNeg;
use serde::de::DeserializeOwned;
use serde::ser::SerializeStruct;
use serde::{Deserialize, Serialize};
use spenso::arithmetic::ScalarMul;
use spenso::contraction::Contract;
use spenso::data::{DataTensor, GetTensorData, StorageTensor};

use spenso::iterators::IteratableTensor;
use spenso::parametric::atomcore::{PatternReplacement, TensorAtomMaps, TensorAtomOps};
use spenso::parametric::{
    EvalTensor, EvalTensorSet, LinearizedEvalTensorSet, MixedTensor, ParamTensorSet,
    SerializableCompiledEvaluator, TensorSet,
};
use spenso::shadowing::ETS;
use spenso::structure::abstract_index::DOWNIND;
use spenso::structure::concrete_index::{ExpandedIndex, FlatIndex};

use spenso::structure::representation::{
    BaseRepName, Bispinor, ColorAdjoint, ColorFundamental, ColorSextet, Minkowski,
};
use spenso::structure::{HasStructure, ScalarTensor, SmartShadowStructure, VecStructure};
use spenso::symbolica_utils::SerializableAtom;
use spenso::symbolica_utils::SerializableSymbol;
use spenso::upgrading_arithmetic::FallibleSub;
use spenso::{
    complex::Complex,
    network::TensorNetwork,
    parametric::ParamTensor,
    shadowing::Shadowable,
    structure::{
        representation::{Lorentz, PhysReps, RepName},
        NamedStructure, TensorStructure,
    },
};
use symbolica::domains::finite_field::PrimeIteratorU64;
use symbolica::domains::rational::Rational;
use symbolica::poly::Variable;
use symbolica::printer::{AtomPrinter, PrintOptions};
use symbolica::state::Workspace;

use crate::numerator::ufo::UFO;
use symbolica::atom::{AtomCore, AtomView, Symbol};
use symbolica::evaluate::{CompileOptions, ExpressionEvaluator, InlineASM};
use symbolica::id::{Context, Match, MatchSettings};

use symbolica::{
    atom::{Atom, FunctionBuilder},
    function, symb,
};
use symbolica::{
    domains::float::NumericalFloatLike,
    evaluate::FunctionMap,
    id::{Pattern, Replacement},
};

pub mod ufo;
#[derive(Debug, Clone, Serialize, Deserialize)]
/// Settings for the numerator
pub struct NumeratorSettings {
    pub eval_settings: NumeratorEvaluatorOptions,
    /// Parse mode for the numerator, once all processing is done. `Polynomial` turns it into a polynomial in the energies, while `Direct` keeps it as is
    pub parse_mode: NumeratorParseMode,
    /// If set, dump the expression the expression at each step in this format
    pub dump_expression: Option<ExpressionFormat>,
    /// If set, instead of deriving the numerator from feynman rules, use this as the numerator
    /// Will be parsed to a symbolica expression
    pub global_numerator: Option<String>,
    /// If set, multiply the numerator by this prefactor
    pub global_prefactor: GlobalPrefactor,
    /// Type of Gamma algebra to use, either symbolic (replacement rules) or concrete (replace by value using spenso)
    pub gamma_algebra: GammaAlgebraMode,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default)]
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
            global_numerator: None,
            global_prefactor: GlobalPrefactor::default(),
            dump_expression: None,
            gamma_algebra: GammaAlgebraMode::Symbolic,
            parse_mode: NumeratorParseMode::Polynomial,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum GammaAlgebraMode {
    Symbolic,
    Concrete,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExtraInfo {
    pub path: PathBuf,
    pub orientations: Vec<Vec<bool>>,
}

pub type AtomStructure = SmartShadowStructure<SerializableSymbol, Vec<SerializableAtom>>;

pub trait Evaluate<T: FloatLike> {
    fn evaluate_all_orientations(
        &mut self,
        emr: &[FourMomentum<F<T>>],
        polarizations: &[Polarization<Complex<F<T>>>],
        tag: Option<Uuid>,
        setting: &Settings,
    ) -> Result<RepeatingIteratorTensorOrScalar<DataTensor<Complex<F<T>>>>>;

    fn evaluate_single(
        &mut self,
        emr: &[FourMomentum<F<T>>],
        polarizations: &[Polarization<Complex<F<T>>>],
        tag: Option<Uuid>,
        setting: &Settings,
    ) -> DataTensor<Complex<F<T>>>;
}

impl<T: FloatLike> Evaluate<T> for Numerator<Evaluators> {
    fn evaluate_all_orientations(
        &mut self,
        emr: &[FourMomentum<F<T>>],
        polarizations: &[Polarization<Complex<F<T>>>],
        _tag: Option<Uuid>,
        settings: &Settings,
    ) -> Result<RepeatingIteratorTensorOrScalar<DataTensor<Complex<F<T>>>>> {
        <T as NumeratorEvaluateFloat>::update_params(self, emr, polarizations, settings);

        if !settings.general.load_compiled_numerator {
            self.disable_compiled();
        }
        if !settings.general.joint_numerator_eval {
            self.disable_combined();
        }
        <T as NumeratorEvaluateFloat>::evaluate_all_orientations(self)
    }

    fn evaluate_single(
        &mut self,
        emr: &[FourMomentum<F<T>>],
        polarizations: &[Polarization<Complex<F<T>>>],
        _tag: Option<Uuid>,
        setting: &Settings,
    ) -> DataTensor<Complex<F<T>>> {
        if !setting.general.load_compiled_numerator {
            self.disable_compiled();
        }

        if !setting.general.joint_numerator_eval {
            self.disable_combined();
        }
        <T as NumeratorEvaluateFloat>::update_params(self, emr, polarizations, setting);
        <T as NumeratorEvaluateFloat>::evaluate_single(self)
    }
}

pub trait NumeratorEvaluateFloat<T: FloatLike = Self> {
    fn evaluate_all_orientations(
        num: &mut Numerator<Evaluators>,
    ) -> Result<RepeatingIteratorTensorOrScalar<DataTensor<Complex<F<T>>>>>;

    fn evaluate_single(num: &mut Numerator<Evaluators>) -> DataTensor<Complex<F<T>>>;

    fn update_params(
        num: &mut Numerator<Evaluators>,
        emr: &[FourMomentum<F<T>>],
        polarizations: &[Polarization<Complex<F<T>>>],
        settings: &Settings,
    );
}

pub struct RepeatingIterator<T> {
    elements: Vec<T>,
    positions: std::vec::IntoIter<usize>,
}

pub enum RepeatingIteratorTensorOrScalar<T: HasStructure> {
    Tensors(RepeatingIterator<T>),
    Scalars(RepeatingIterator<T::Scalar>),
}

impl<T> RepeatingIterator<T> {
    pub fn new(positions: Vec<usize>, elements: Vec<T>) -> Self {
        RepeatingIterator {
            elements,
            positions: positions.into_iter(),
        }
    }

    pub fn new_not_repeating(elements: Vec<T>) -> Self {
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

impl<T> LendingIterator for RepeatingIterator<T> {
    type Item<'a>
        = &'a T
    where
        Self: 'a;

    fn next(&mut self) -> Option<Self::Item<'_>> {
        Some(&self.elements[self.positions.next()?])
    }
}

impl NumeratorEvaluateFloat for f64 {
    fn update_params(
        num: &mut Numerator<Evaluators>,
        emr: &[FourMomentum<F<Self>>],
        polarizations: &[Polarization<Complex<F<Self>>>],
        settings: &Settings,
    ) {
        let params = &mut num.state.double_param_values;
        emr.iter()
            .flat_map(|p| p.into_iter().cloned().map(Complex::new_re))
            .enumerate()
            .for_each(|(i, c)| params[i] = c);

        assert_eq!(emr.len(), num.state.emr_len);

        let mut i = 4 * emr.len();

        for p in polarizations {
            if !p.is_scalar() {
                for &pi in p {
                    assert!(
                        i < num.state.model_params_start - 1,
                        "polarization index out of bounds"
                    );
                    params[i] = pi;
                    i += 1;
                }
            }
        }

        if settings.general.debug > 0 {
            DEBUG_LOGGER.write("params", &params);
        }
    }

    fn evaluate_all_orientations(
        num: &mut Numerator<Evaluators>,
    ) -> Result<RepeatingIteratorTensorOrScalar<DataTensor<Complex<F<Self>>>>> {
        let params = &mut num.state.double_param_values;
        if num.state.single.param_len != params.len() {
            return Err(eyre!(
                "params length mismatch {} not equal to {}",
                num.state.single.param_len,
                params.len()
            ));
        }

        if let (Some(orientation_evaluator), SingleOrCombined::Combined) =
            (&mut num.state.orientated, &num.state.choice)
        {
            let pos = orientation_evaluator.positions.clone();
            let res = match &mut orientation_evaluator.compiled {
                CompiledEvaluator {
                    state: CompiledState::Enabled,
                    evaluator: Some(s),
                } => {
                    trace!("using compiled evaluator");
                    s.evaluate(params)
                }
                _ => orientation_evaluator.eval_double.evaluate(params),
            };

            Ok((res, pos).into())
        } else {
            trace!("no orientation evaluator using single evaluator");
            let base_params = params.clone();
            let mut tensors = Vec::new();
            let orientations = &num.state.orientations;
            for o in orientations {
                for (i, &sign) in o.iter().enumerate() {
                    if !sign {
                        params[4 * i] = -base_params[4 * i];
                    } else {
                        params[4 * i] = base_params[4 * i];
                    }
                }

                let t = if let Some(c) = &mut num.state.single.compiled.evaluator {
                    c.evaluate(params)
                } else {
                    num.state.single.eval_double.evaluate(params)
                };
                tensors.push(t);
            }
            let tensors = TensorSet::from_iter(tensors);
            Ok(tensors.into())
        }
    }

    fn evaluate_single(num: &mut Numerator<Evaluators>) -> DataTensor<Complex<F<Self>>> {
        let params = &num.state.double_param_values;
        if num.state.single.param_len != params.len() {
            panic!("params length mismatch");
        }
        if let Some(c) = &mut num.state.single.compiled.evaluator {
            c.evaluate(params)
        } else {
            num.state.single.eval_double.evaluate(params)
        }
    }
}

impl NumeratorEvaluateFloat for f128 {
    fn update_params(
        num: &mut Numerator<Evaluators>,
        emr: &[FourMomentum<F<Self>>],
        polarizations: &[Polarization<Complex<F<Self>>>],
        settings: &Settings,
    ) {
        let params = &mut num.state.quad_param_values;
        emr.iter()
            .flat_map(|p| p.into_iter().cloned().map(Complex::new_re))
            .enumerate()
            .for_each(|(i, c)| params[i] = c);

        let mut i = 4 * emr.len();

        for p in polarizations {
            if !p.is_scalar() {
                for pi in p {
                    params[i] = pi.clone();
                    i += 1;
                }
            }
        }

        if settings.general.debug > 0 {
            DEBUG_LOGGER.write("params", &params);
        }
    }

    fn evaluate_all_orientations(
        num: &mut Numerator<Evaluators>,
    ) -> Result<RepeatingIteratorTensorOrScalar<DataTensor<Complex<F<Self>>>>> {
        let params = &mut num.state.quad_param_values;
        if num.state.single.param_len != params.len() {
            return Err(eyre!("params length mismatch"));
        }
        if let (Some(orientation_evaluator), SingleOrCombined::Combined) =
            (&mut num.state.orientated, &num.state.choice)
        {
            let pos = orientation_evaluator.positions.clone();
            let res = orientation_evaluator.eval_quad.evaluate(params);
            Ok((res, pos).into())
        } else {
            let base_params = params.clone();
            let mut tensors = Vec::new();
            let orientations = &num.state.orientations;
            for o in orientations {
                for (i, &sign) in o.iter().enumerate() {
                    if !sign {
                        params[4 * i] = -(base_params[4 * i].clone());
                    } else {
                        params[4 * i] = base_params[4 * i].clone();
                    }
                }
                let t = num.state.single.eval_quad.evaluate(params);
                tensors.push(t);
            }
            let tensors = TensorSet::from_iter(tensors);
            Ok(tensors.into())
        }
    }

    fn evaluate_single(num: &mut Numerator<Evaluators>) -> DataTensor<Complex<F<Self>>> {
        let params = &num.state.quad_param_values;
        if num.state.single.param_len != params.len() {
            panic!("params length mismatch");
        }
        num.state.single.eval_quad.evaluate(params)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
pub struct Numerator<State> {
    pub state: State,
}

impl<S: NumeratorState> Numerator<S> {
    pub fn export(&self) -> String {
        self.state.export()
    }

    pub fn forget_type(self) -> Numerator<PythonState> {
        Numerator {
            state: self.state.forget_type(),
        }
    }

    pub fn update_model(&mut self, model: &Model) -> Result<()> {
        self.state.update_model(model)
    }

    fn add_consts_to_fn_map(fn_map: &mut FunctionMap) {
        fn_map.add_constant(Atom::parse("Nc").unwrap(), Rational::from_unchecked(3, 1));

        fn_map.add_constant(Atom::parse("TR").unwrap(), Rational::from_unchecked(1, 2));

        fn_map.add_constant(Atom::parse("pi").unwrap(), std::f64::consts::PI.into());
    }
}

impl<S: GetSingleAtom> Numerator<S> {
    pub fn get_single_atom(&self) -> Result<SerializableAtom, NumeratorStateError> {
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
    pub fn try_from<S: TypedNumeratorState>(self) -> Result<Numerator<S>, Report> {
        Ok(Numerator {
            state: self.state.try_into()?,
        })
    }

    pub fn apply<F, S: TypedNumeratorState, T: TypedNumeratorState>(
        &mut self,
        f: F,
    ) -> Result<(), NumeratorStateError>
    where
        F: FnMut(Numerator<S>) -> Numerator<T>,
    {
        S::apply(self, f)
    }
}
pub trait NumeratorState: Clone + Debug + Encode + Decode<StateMap> {
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

#[derive(Debug, Clone)]
pub struct GlobalPrefactor {
    pub color: Atom,
    pub colorless: Atom,
}

impl Default for GlobalPrefactor {
    fn default() -> Self {
        GlobalPrefactor {
            color: Atom::new_num(1),
            colorless: Atom::new_num(1),
        }
    }
}

impl Serialize for GlobalPrefactor {
    fn serialize<S>(&self, serializer: S) -> std::result::Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let mut state = serializer.serialize_struct("GlobalPrefactor", 2)?;
        state.serialize_field("color", &self.color.to_string())?;
        state.serialize_field("colorless", &self.colorless.to_string())?;
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
            color: String,
            colorless: String,
        }
        let helper = GlobalPrefactorHelper::deserialize(deserializer)?;
        Ok(GlobalPrefactor {
            color: Atom::parse(&helper.color).unwrap(),
            colorless: Atom::parse(&helper.colorless).unwrap(),
        })
    }
}

impl Numerator<UnInit> {
    pub fn from_graph(
        self,
        graph: &BareGraph,
        prefactor: &GlobalPrefactor,
    ) -> Numerator<AppliedFeynmanRule> {
        debug!("Applying feynman rules");
        let state = AppliedFeynmanRule::from_graph(graph, prefactor);
        debug!(
            "Applied feynman rules:\n\tcolor:{}\n\tcolorless:{}",
            state.color, state.colorless
        );
        Numerator { state }
    }

    pub fn from_global(
        self,
        global: Atom,
        // _graph: &BareGraph,
        prefactor: &GlobalPrefactor,
    ) -> Numerator<Global> {
        debug!("Setting global numerator");
        let state = {
            let mut global = global;
            global = global * &prefactor.color * &prefactor.colorless;
            Global::new(global.into())
        };
        debug!(
            "Global numerator:\n\tcolor:{}\n\tcolorless:{}",
            state.color, state.colorless
        );
        Numerator { state }
    }

    pub fn from_global_color(
        self,
        global: Atom,
        // _graph: &BareGraph,
        prefactor: &GlobalPrefactor,
    ) -> Numerator<Global> {
        debug!("Setting global numerator");
        let state = {
            let mut global = global;
            global = global * &prefactor.color * &prefactor.colorless;
            Global::new_color(global.into())
        };
        debug!(
            "Global numerator:\n\tcolor:{}\n\tcolorless:{}",
            state.color, state.colorless
        );
        Numerator { state }
    }
}

#[allow(clippy::default_constructed_unit_structs)]
impl Default for Numerator<UnInit> {
    fn default() -> Self {
        Numerator {
            state: UnInit::default(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
pub struct SymbolicExpression<State> {
    #[bincode(with_serde)]
    pub colorless: DataTensor<SerializableAtom>,
    #[bincode(with_serde)]
    pub color: DataTensor<SerializableAtom>,
    pub state: State,
}

pub trait GetSingleAtom {
    fn get_single_atom(&self) -> Result<SerializableAtom, NumeratorStateError>;
}

pub trait UnexpandedNumerator: NumeratorState + GetSingleAtom {
    // fn expr(&self) -> Result<SerializableAtom, NumeratorStateError>;

    fn map_color(self, f: impl Fn(SerializableAtom) -> SerializableAtom) -> Self;

    fn map_color_mut(&mut self, f: impl FnMut(&mut SerializableAtom));

    fn map_colorless(self, f: impl Fn(SerializableAtom) -> SerializableAtom) -> Self;
}

impl<E: ExpressionState> GetSingleAtom for SymbolicExpression<E> {
    fn get_single_atom(&self) -> Result<SerializableAtom, NumeratorStateError> {
        self.colorless
            .contract(&self.color)
            .map_err(|n| NumeratorStateError::Any(n.into()))?
            .scalar()
            .ok_or(NumeratorStateError::Any(eyre!("not a scalar")))
    }
}

impl<E: ExpressionState> SymbolicExpression<E> {
    #[allow(dead_code)]
    fn map_color(self, f: impl Fn(SerializableAtom) -> SerializableAtom) -> Self {
        SymbolicExpression {
            colorless: self.colorless,
            color: self.color.map_data(f),
            state: self.state,
        }
    }

    #[allow(dead_code)]
    fn map_color_mut(&mut self, f: impl FnMut(&mut SerializableAtom)) {
        self.color.map_data_mut(f);
    }
    #[allow(dead_code)]
    fn map_colorless(self, f: impl Fn(SerializableAtom) -> SerializableAtom) -> Self {
        SymbolicExpression {
            colorless: self.colorless.map_data(f),
            color: self.color,
            state: self.state,
        }
    }
}

impl<E: ExpressionState> SymbolicExpression<E> {
    pub fn new(expression: SerializableAtom) -> Self {
        E::new(expression)
    }
    pub fn new_color(expression: SerializableAtom) -> Self {
        E::new_color(expression)
    }
}

pub trait ExpressionState:
    Serialize + Clone + DeserializeOwned + Debug + Encode + Decode<StateMap> + Default
{
    fn forget_type(data: SymbolicExpression<Self>) -> PythonState;

    fn new(expression: SerializableAtom) -> SymbolicExpression<Self> {
        SymbolicExpression {
            colorless: DataTensor::new_scalar(expression),
            color: DataTensor::new_scalar(Atom::new_num(1).into()),
            state: Self::default(),
        }
    }

    fn new_color(expression: SerializableAtom) -> SymbolicExpression<Self> {
        SymbolicExpression {
            color: DataTensor::new_scalar(expression),
            colorless: DataTensor::new_scalar(Atom::new_num(1).into()),
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
    pub fn write(
        &self,
        extra_info: &ExtraInfo,
        bare_graph: &BareGraph,
        numerator_format: ExpressionFormat,
    ) -> std::io::Result<()> {
        debug!("Writing global numerator");
        let file_path = extra_info.path.join("expressions");

        let printer_ops = numerator_format.into();

        let out = format!(
            "{}",
            AtomPrinter::new_with_options(self.get_single_atom().unwrap().0.as_view(), printer_ops)
        );

        fs::write(
            file_path.join(format!("global_{}_exp.json", bare_graph.name)),
            serde_json::to_string_pretty(&out).unwrap(),
        )?;

        Ok(())
    }

    pub fn color_simplify(mut self) -> Numerator<ColorSimplified> {
        debug!("Color simplifying global numerator");
        let mut fully_simplified = true;

        let colorless =
            self.state
                .colorless
                .map_data_ref_mut(|a| match Self::color_simplify_global_impl(a) {
                    Ok(expression) => expression,
                    Err(ColorError::NotFully(expression)) => {
                        fully_simplified = false;
                        expression
                    }
                });
        let color =
            self.state
                .color
                .map_data_ref_mut(|a| match ColorSimplified::color_simplify_impl(a) {
                    Ok(expression) => expression,
                    Err(ColorError::NotFully(expression)) => {
                        fully_simplified = false;
                        expression
                    }
                });
        let state = ColorSimplified {
            colorless,
            color,
            state: if fully_simplified {
                Color::Fully
            } else {
                Color::Partially
            },
        };
        debug!(
            "Color simplified numerator: color:{}\n colorless:{}",
            state.color, state.colorless
        );
        Numerator { state }
    }

    fn color_simplify_global_impl(
        expression: &SerializableAtom,
    ) -> Result<SerializableAtom, ColorError> {
        let mut expression = expression.clone();
        ColorSimplified::isolate_color(&mut expression);
        let color_pat = function!(GS.color_wrap, GS.x_).to_pattern();

        let mut fully = true;

        let color_simplified = {
            let pat = expression
                .0
                .pattern_match(&color_pat, None, None)
                .next()
                .unwrap();

            match ColorSimplified::color_simplify_impl(&pat[&GS.x_].clone().into()) {
                Ok(s) => s,
                Err(ColorError::NotFully(s)) => {
                    fully = false;
                    s
                }
            }
            .0
            .factor()
        };

        let mut expression =
            expression
                .0
                .replace_all(&color_pat, Atom::new_num(1).to_pattern(), None, None);
        expression = expression * color_simplified;

        if fully {
            Ok(expression.into())
        } else {
            Err(ColorError::NotFully(expression.into()))
        }
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
        self.get_single_atom().unwrap().0.to_canonical_string()
    }

    fn forget_type(self) -> PythonState {
        State::forget_type(self)
    }

    fn update_model(&mut self, _model: &Model) -> Result<()> {
        Err(eyre!("Only an expression, nothing to update"))
    }
}

impl AppliedFeynmanRule {
    pub fn simplify_ids(&mut self) {
        let replacements: Vec<Replacement> = vec![
            Replacement::new(
                Pattern::parse("f_(i_,aind(loru(a__)))*id(aind(lord(a__),loru(b__)))").unwrap(),
                Pattern::parse("f_(i_,aind(loru(b__)))").unwrap(),
            ),
            Replacement::new(
                Pattern::parse("f_(i_,aind(lord(a__)))*id(aind(loru(a__),lord(b__)))").unwrap(),
                Pattern::parse("f_(i_,aind(lord(b__)))").unwrap(),
            ),
            Replacement::new(
                Pattern::parse("f_(i_,aind(loru(a__)))*id(aind(loru(b__),lord(a__)))").unwrap(),
                Pattern::parse("f_(i_,aind(loru(b__)))").unwrap(),
            ),
            Replacement::new(
                Pattern::parse("f_(i_,aind(lord(a__)))*id(aind(lord(b__),loru(a__)))").unwrap(),
                Pattern::parse("f_(i_,aind(lord(b__)))").unwrap(),
            ),
        ];

        let settings = MatchSettings {
            rhs_cache_size: 0,
            ..Default::default()
        };

        let reps: Vec<Replacement> = replacements
            .into_iter()
            .map(|r| r.with_settings(settings.clone()))
            .collect();

        self.colorless
            .map_data_mut(|a| a.replace_all_multiple_repeat_mut(&reps));
    }
    pub fn from_graph(graph: &BareGraph, prefactor: &GlobalPrefactor) -> Self {
        debug!("Generating numerator for graph: {}", graph.name);
        // debug!("momentum: {}", graph.dot_lmb());

        let vatoms: Vec<_> = graph
            .vertices
            .iter()
            .flat_map(|v| v.colorless_vertex_rule(graph))
            .collect();

        let mut eatoms: Vec<_> = vec![];
        let i = Atom::new_var(Atom::I);
        for (j, e) in graph.edges.iter().enumerate() {
            let [n, c] = e.color_separated_numerator(graph, j);
            let n = if matches!(e.edge_type, EdgeType::Virtual) {
                &n * &i
            } else {
                n
            };
            eatoms.push([n, c]);
            // shift += s;
            // graph.shifts.0 += shift;
        }
        let mut colorless_builder = DataTensor::new_scalar(Atom::new_num(1));

        let mut colorful_builder = DataTensor::new_scalar(Atom::new_num(1));

        for [colorless, color] in &vatoms {
            // println!("colorless vertex: {}", colorless);
            // println!("colorful vertex: {}", color);
            colorless_builder = colorless_builder.contract(colorless).unwrap();
            colorful_builder = colorful_builder.contract(color).unwrap();
            // println!("vertex: {v}");
            // builder = builder * v;
        }

        for [n, c] in &eatoms {
            // println!("colorless edge {n}");
            // println!("colorfull edge {c}");
            colorless_builder = colorless_builder.scalar_mul(n).unwrap();
            colorful_builder = colorful_builder.scalar_mul(c).unwrap();
        }

        colorless_builder = colorless_builder.scalar_mul(&prefactor.colorless).unwrap();
        colorful_builder = colorful_builder.scalar_mul(&prefactor.color).unwrap();

        let mut num = AppliedFeynmanRule {
            colorless: colorless_builder.map_data(|a| normalise_complex(&a).into()),
            color: colorful_builder.map_data(|a| normalise_complex(&a).into()),
            state: Default::default(),
        };
        num.simplify_ids();
        num
    }
}

impl<T: Copy + Default> Numerator<SymbolicExpression<T>> {
    pub fn canonize_lorentz(&self) -> Result<Self, String> {
        let pats: Vec<_> = vec![Minkowski::selfless_symbol(), Bispinor::selfless_symbol()];

        let mut indices_map = AHashMap::new();

        self.state.colorless.iter_flat().for_each(|(_, v)| {
            for p in &pats {
                for a in
                    v.0.pattern_match(&function!(*p, GS.x_, GS.y_).to_pattern(), None, None)
                {
                    indices_map.insert(
                        function!(*p, a[&GS.x_], a[&GS.y_]),
                        function!(*p, a[&GS.x_]),
                    );
                }
            }
        });

        let sorted = indices_map.into_iter().sorted().collect::<Vec<_>>();

        let color = self.state.color.clone(); //map_data_ref_result(|a| a.0.canonize_tensors(&indices, None).map(|a| a.into()))?;

        let colorless = self
            .state
            .colorless
            .map_data_ref_result(|a| a.0.canonize_tensors(&sorted).map(|a| a.into()))?;

        Ok(Self {
            state: SymbolicExpression {
                colorless,
                color,
                state: T::default(),
            },
        })
    }

    pub fn canonize_color(&self) -> Result<Self, String> {
        let pats: Vec<_> = vec![ColorAdjoint::selfless_symbol()];
        let dualizablepats: Vec<_> = vec![
            ColorFundamental::selfless_symbol(),
            ColorSextet::selfless_symbol(),
        ];

        let mut color = self.state.color.map_data_ref(|a| {
            SerializableAtom::from(a.0.replace_all(
                &function!(symb!(DOWNIND), GS.x__).to_pattern(),
                Atom::new_var(GS.x__).to_pattern(),
                None,
                None,
            ))
        });

        let mut indices_map = AHashMap::new();

        color.iter_flat().for_each(|(_, v)| {
            for p in pats.iter().chain(&dualizablepats) {
                for a in
                    v.0.pattern_match(&function!(*p, GS.x_, GS.y_).to_pattern(), None, None)
                {
                    indices_map.insert(
                        function!(*p, a[&GS.x_], a[&GS.y_]),
                        function!(*p, a[&GS.x_]),
                    );
                }
            }
        });

        let sorted = indices_map.into_iter().sorted().collect::<Vec<_>>();
        // println!(
        //     "indices sorted [{}]",
        //     sorted
        //         .iter()
        //         .map(|(a, b)| format!(
        //             "(Atom::parse(\"{}\").unwrap(),Atom::parse(\"{}\").unwrap())",
        //             a, b
        //         ))
        //         .collect::<Vec<_>>()
        //         .join(", ")
        // );

        color = color.map_data_ref_result(|a| a.0.canonize_tensors(&sorted).map(|a| a.into()))?;

        let colorless = self.state.colorless.clone();
        Ok(Self {
            state: SymbolicExpression {
                colorless,
                color,
                state: T::default(),
            },
        })
    }
}

impl Numerator<AppliedFeynmanRule> {
    pub fn write(
        &self,
        extra_info: &ExtraInfo,
        bare_graph: &BareGraph,
        numerator_format: ExpressionFormat,
    ) -> std::io::Result<()> {
        debug!("Writing local numerator");
        let file_path = extra_info.path.join("expressions");

        let printer_ops = numerator_format.into();

        let out = format!(
            "{}",
            AtomPrinter::new_with_options(self.get_single_atom().unwrap().0.as_view(), printer_ops)
        );

        fs::write(
            file_path.join(format!("local_{}_exp.json", bare_graph.name)),
            serde_json::to_string_pretty(&out).unwrap(),
        )?;

        Ok(())
    }

    pub fn color_simplify(self) -> Numerator<ColorSimplified> {
        debug!("Color simplifying local numerator");

        let state = ColorSimplified::color_simplify(self.state);

        debug!(
            "Color simplified numerator: color:{}\n colorless:{}",
            state.color, state.colorless
        );

        Numerator { state }
    }
}

#[derive(Debug, Error)]
pub enum ColorError {
    #[error("Not fully simplified: {0}")]
    NotFully(SerializableAtom),
}

impl ColorSimplified {
    fn isolate_color(expression: &mut SerializableAtom) {
        let color_fn = FunctionBuilder::new(Symbol::new("color"))
            .add_arg(Atom::new_num(1))
            .finish();
        expression.0 = &expression.0 * color_fn;
        let replacements = vec![
            (
                Pattern::parse("f_(x___,cof(3,i_),z___)*color(a___)").unwrap(),
                Pattern::parse("color(a___*f_(x___,cof(3,i_),z___))").unwrap(),
            ),
            (
                Pattern::parse("f_(x___,dind(cof(3,i_)),z___)*color(a___)").unwrap(),
                Pattern::parse("color(a___*f_(x___,dind(cof(3,i_)),z___))").unwrap(),
            ),
            (
                Pattern::parse("f_(x___,coad(i__),z___)*color(a___)").unwrap(),
                Pattern::parse("color(a___*f_(x___,coad(i__),z___))").unwrap(),
            ),
        ];
        let settings = MatchSettings {
            rhs_cache_size: 0,
            ..Default::default()
        };

        let reps: Vec<Replacement> = replacements
            .into_iter()
            .map(|(lhs, rhs)| Replacement::new(lhs, rhs).with_settings(settings.clone()))
            .collect();

        expression.replace_all_multiple_repeat_mut(&reps);

        // let mut color = Atom::new();

        // if let AtomView::Mul(mul) = expression.0.as_view() {
        //     for a in mul {
        //         if let AtomView::Fun(f) = a {
        //             if f.get_symbol() == Symbol::new("color") {
        //                 color.set_from_view(&f.iter().next().unwrap());
        //             }
        //         }
        //     }
        // }

        // expression.0.replace_all(
        //     &Pattern::Fn(Symbol::new("color"), vec![]),
        //     &Pattern::Literal(Atom::new_num(1)).into(),
        //     None,
        //     None,
        // );

        // color.into()
    }

    pub fn color_simplify_impl(
        expression: &SerializableAtom,
    ) -> Result<SerializableAtom, ColorError> {
        let f_ = symb!("f_");
        let cof = ColorFundamental::rep(3);
        let coaf = ColorFundamental::rep(3).dual();
        let coad = ColorAdjoint::rep(8);
        let a___ = Atom::new_var(symb!("a___"));
        let a_ = Atom::new_var(symb!("a_"));
        let b_ = Atom::new_var(symb!("b_"));
        let c___ = Atom::new_var(symb!("c___"));
        let c_ = Atom::new_var(symb!("c_"));
        let d_ = Atom::new_var(symb!("d_"));
        let e_ = Atom::new_var(symb!("e_"));
        let i_ = Atom::new_var(symb!("i_"));
        let i = symb!("i");
        let j = symb!("j");
        let k = symb!("k");

        let tr = Atom::new_var(symb!("TR"));
        let nc = Atom::new_var(symb!("Nc"));
        let reps = vec![
            (
                function!(f_, a___, cof.pattern(&b_), c___)
                    * function!(ETS.id, coaf.pattern(&b_), cof.pattern(&c_)),
                function!(f_, a___, cof.pattern(&c_), c___),
            ),
            // (
            //     function!(f_, a___, cof.pattern(&c_), c___)
            //         * function!(ETS.id, cof.pattern(&b_), coaf.pattern(&c_)),
            //     function!(f_, a___, cof.pattern(&b_), c___),
            // ),
            (
                function!(f_, a___, coaf.pattern(&b_), c___)
                    * function!(ETS.id, cof.pattern(&b_), coaf.pattern(&c_)),
                function!(f_, a___, coaf.pattern(&c_), c___),
            ),
            // (
            //     function!(f_, a___, coaf.pattern(&c_), c___)
            //         * function!(ETS.id, coaf.pattern(&b_), cof.pattern(&c_)),
            //     function!(f_, a___, coaf.pattern(&b_), c___),
            // ),
            (
                function!(f_, a___, coad.pattern(&b_), c___)
                    * function!(ETS.id, coad.pattern(&b_), coad.pattern(&a_)),
                function!(f_, a___, coad.pattern(&a_), c___),
            ),
            (
                function!(f_, a___, coad.pattern(&a_), c___)
                    * function!(ETS.id, coad.pattern(&b_), coad.pattern(&a_)),
                function!(f_, a___, coad.pattern(&b_), c___),
            ),
            (
                function!(ETS.id, coaf.pattern(&a_), cof.pattern(&a_)),
                nc.clone(),
            ),
            (
                function!(ETS.id, cof.pattern(&a_), coaf.pattern(&a_)),
                nc.clone(),
            ),
            (
                function!(ETS.id, coad.pattern(&a_), coad.pattern(&a_)),
                (&nc * &nc) - 1,
            ),
            (
                function!(UFO.t, a_, cof.pattern(&b_), coaf.pattern(&b_)),
                Atom::new_num(0),
            ),
            (
                function!(UFO.t, a_, cof.pattern(&c_), coaf.pattern(&e_))
                    * function!(UFO.t, b_, cof.pattern(&e_), coaf.pattern(&c_)),
                &tr * function!(ETS.id, a_, b_),
            ),
            (
                function!(UFO.t, a_, cof.pattern(&c_), coaf.pattern(&e_)).pow(Atom::new_num(2)),
                &tr * function!(ETS.id, a_, a_),
            ),
            (
                function!(UFO.t, &e_, a_, b_) * function!(UFO.t, &e_, c_, d_),
                &tr * (function!(ETS.id, a_, d_) * function!(ETS.id, c_, b_)
                    - (function!(ETS.id, a_, b_) * function!(ETS.id, c_, d_) / &nc)),
            ),
            (
                function!(UFO.t, i_, a_, coaf.pattern(&b_))
                    * function!(UFO.t, e_, cof.pattern(&b_), coaf.pattern(&c_))
                    * function!(UFO.t, i_, cof.pattern(&c_), d_),
                -(&tr / &nc) * function!(UFO.t, e_, a_, d_),
            ),
            (
                function!(
                    UFO.f,
                    coad.pattern(&a_),
                    coad.pattern(&b_),
                    coad.pattern(&c_)
                )
                .pow(Atom::new_num(2)),
                &nc * (&nc * &nc - 1),
            ),
        ];

        let frep = [Replacement::new(
            function!(
                UFO.f,
                coad.pattern(&a_),
                coad.pattern(&b_),
                coad.pattern(&c_)
            )
            .to_pattern(),
            (((function!(
                UFO.t,
                coad.pattern(&a_),
                cof.pattern(&function!(i, a_, b_, c_)),
                coaf.pattern(&function!(j, a_, b_, c_))
            ) * function!(
                UFO.t,
                coad.pattern(&b_),
                cof.pattern(&function!(j, a_, b_, c_)),
                coaf.pattern(&function!(k, a_, b_, c_))
            ) * function!(
                UFO.t,
                coad.pattern(&c_),
                cof.pattern(&function!(k, a_, b_, c_)),
                coaf.pattern(&function!(i, a_, b_, c_))
            ) - function!(
                UFO.t,
                coad.pattern(&a_),
                cof.pattern(&function!(i, a_, b_, c_)),
                coaf.pattern(&function!(j, a_, b_, c_))
            ) * function!(
                UFO.t,
                coad.pattern(&c_),
                cof.pattern(&function!(j, a_, b_, c_)),
                coaf.pattern(&function!(k, a_, b_, c_))
            ) * function!(
                UFO.t,
                coad.pattern(&b_),
                cof.pattern(&function!(k, a_, b_, c_)),
                coaf.pattern(&function!(i, a_, b_, c_))
            )) / &tr)
                * Atom::new_var(Atom::I).ref_neg())
            .to_pattern(),
        )];

        let settings = MatchSettings {
            rhs_cache_size: 0,
            ..Default::default()
        };
        let replacements: Vec<Replacement> = reps
            .into_iter()
            .map(|(a, b)| {
                Replacement::new(a.to_pattern(), b.to_pattern()).with_settings(settings.clone())
            })
            .collect();

        let mut atom = Atom::new_num(0);
        // for r in &replacements {
        //     println!("{r}")
        // }

        let mut expression = expression.clone();
        let mut first = true;
        while first
            || expression
                .0
                .replace_all_multiple_into(&replacements, &mut atom)
        {
            if !first {
                std::mem::swap(&mut expression.0, &mut atom)
            };
            first = false;
            expression.0 = expression.0.replace_all_multiple(&frep);
            expression.0 = expression.0.expand();
        }

        let pats: Vec<_> = vec![ColorAdjoint::selfless_symbol()];
        let dualizablepats: Vec<_> = vec![
            ColorFundamental::selfless_symbol(),
            ColorSextet::selfless_symbol(),
        ];

        let mut fully_simplified = true;
        for p in pats.iter().chain(&dualizablepats) {
            if expression
                .0
                .pattern_match(&function!(*p, GS.x_, GS.y_).to_pattern(), None, None)
                .next()
                .is_some()
            {
                fully_simplified = false;
            }
        }

        if fully_simplified {
            Ok(expression)
        } else {
            Err(ColorError::NotFully(expression))
        }
    }

    pub fn color_simplify<T: ExpressionState>(mut expr: SymbolicExpression<T>) -> ColorSimplified {
        let colorless = expr.colorless;
        let mut fully_simplified = true;

        expr.color.map_data_mut(|a| {
            *a = match Self::color_simplify_impl(a) {
                Ok(expression) => expression,
                Err(ColorError::NotFully(expression)) => {
                    fully_simplified = false;
                    expression
                }
            }
        });

        ColorSimplified {
            colorless,
            color: expr.color,
            state: if fully_simplified {
                Color::Fully
            } else {
                Color::Partially
            },
        }
    }

    /// Parses the color simplified numerator into a network,
    /// If the color hasn't been fully simplified,
    pub fn parse(self) -> Network {
        let a = match self.state {
            Color::Fully => self.get_single_atom().unwrap().0,
            Color::Partially => {
                warn!(
                    "Not fully simplified, taking the colorless part associated to {}",
                    self.color.get_owned_linear(FlatIndex::from(0)).unwrap()
                );
                self.colorless
                    .get_owned_linear(FlatIndex::from(0))
                    .unwrap()
                    .0
            }
        };

        let net = TensorNetwork::<MixedTensor<f64, AtomStructure>, SerializableAtom>::try_from(
            a.as_view(),
        )
        .unwrap()
        .to_fully_parametric()
        .cast();

        // println!("net scalar{}", net.scalar.as_ref().unwrap());
        Network { net }
    }
}

pub type Gloopoly =
    symbolica::poly::polynomial::MultivariatePolynomial<symbolica::domains::atom::AtomField, u8>;
impl Numerator<ColorSimplified> {
    pub fn validate_against_branches(&self, seed: usize) -> bool {
        let gamma = self.clone().gamma_simplify().parse();
        let nogamma = self.clone().parse();
        let mut iter = PrimeIteratorU64::new(seed as u64);
        let reps = nogamma.random_concretize_reps(Some(&mut iter), false);

        let reps_view = reps
            .iter()
            .map(|(a, b)| (a.as_view(), b.as_view()))
            .collect::<Vec<_>>();
        let gammat = gamma
            .apply_reps(reps_view.clone())
            .contract::<Rational>(ContractionSettings::Normal)
            .unwrap()
            .state
            .tensor;

        let nogammat = nogamma
            .apply_reps(reps_view)
            .contract::<Rational>(ContractionSettings::Normal)
            .unwrap()
            .state
            .tensor;

        // println!("{gammat}\n{nogammat}");

        gammat
            .tensor
            .sub_fallible(&nogammat.tensor)
            .unwrap()
            .iter_flat()
            .all(|(_, a)| {
                let a = a.expand();
                // println!("{a}");
                a.is_zero()
            })
    }

    pub fn write(
        &self,
        extra_info: &ExtraInfo,
        bare_graph: &BareGraph,
        numerator_format: ExpressionFormat,
    ) -> std::io::Result<()> {
        debug!("Writing color simplified numerator");
        let file_path = extra_info.path.join("expressions");

        let printer_ops = numerator_format.into();

        let out = format!(
            "{}",
            AtomPrinter::new_with_options(self.get_single_atom().unwrap().0.as_view(), printer_ops)
        );

        fs::write(
            file_path.join(format!("color_simplified_{}_exp.json", bare_graph.name)),
            serde_json::to_string_pretty(&out).unwrap(),
        )?;

        Ok(())
    }

    pub fn gamma_simplify(self) -> Numerator<GammaSimplified> {
        debug!("Gamma simplifying color symplified numerator");

        let gamma_simplified = self
            .state
            .colorless
            .map_data(GammaSimplified::gamma_symplify_impl);

        Numerator {
            state: GammaSimplified {
                colorless: gamma_simplified,
                color: self.state.color,
                state: Default::default(),
            },
        }
    }

    pub fn parse_poly(self, bare_graph: &BareGraph) -> Numerator<PolySplit> {
        debug!("Parsing color simplified numerator into polynomial");
        let state = PolySplit::from_color_simplified(self, bare_graph);
        Numerator { state }
    }

    pub fn parse(self) -> Numerator<Network> {
        debug!("Parsing color simplified numerator into network");

        // debug!("colorless: {}", self.state.colorless);
        // debug!("color: {}", self.state.color);
        let state = self.state.parse();

        // debug!("");
        Numerator { state }
    }
}

#[derive(Debug, Clone)]
pub struct PolySplit {
    pub colorless: DataTensor<Gloopoly>,
    pub var_map: Arc<Vec<Variable>>,
    pub energies: Vec<usize>,
    pub color: ParamTensor,
    pub colorsimplified: Color,
}

impl PolySplit {
    fn generate_variables(bare_graph: &BareGraph) -> Vec<Variable> {
        let mut vars: Vec<Variable> = vec![];

        for (i, e) in bare_graph.edges.iter().enumerate() {
            if let EdgeType::Virtual = e.edge_type {
                let named_structure: NamedStructure<String> = NamedStructure::from_iter(
                    [PhysReps::new_slot(Lorentz {}.into(), 4, i)],
                    "Q".into(),
                    Some(i),
                );

                vars.push(
                    named_structure
                        .to_shell()
                        .expanded_shadow()
                        .unwrap()
                        .get_owned_linear(0.into())
                        .unwrap()
                        .into(),
                );
            }
        }

        vars
    }

    pub fn from_color_out(color_simplified: Numerator<ColorSimplified>) -> DataTensor<Atom> {
        let colorless_parsed = color_simplified
            .state
            .colorless
            .map_data(|a| {
                let mut net =
                    TensorNetwork::<MixedTensor<f64, AtomStructure>, SerializableAtom>::try_from(
                        a.0.as_view(),
                    )
                    .unwrap();
                net.contract().unwrap();
                net.to_fully_parametric()
                    .result_tensor_smart()
                    .unwrap()
                    .tensor
                    .map_structure(VecStructure::from)
            })
            .flatten(&Atom::new_num(0))
            .unwrap();

        colorless_parsed
    }

    fn from_color_simplified(
        color_simplified: Numerator<ColorSimplified>,
        bare_graph: &BareGraph,
    ) -> Self {
        let var_map = Arc::new(Self::generate_variables(bare_graph));
        let energies: Vec<_> = (0..var_map.len()).collect();

        let colorless_parsed = color_simplified
            .state
            .colorless
            .map_data(|a| {
                let mut net =
                    TensorNetwork::<MixedTensor<f64, AtomStructure>, SerializableAtom>::try_from(
                        a.0.as_view(),
                    )
                    .unwrap();
                net.contract().unwrap();
                net.to_fully_parametric()
                    .result_tensor_smart()
                    .unwrap()
                    .tensor
                    .map_structure(VecStructure::from)
            })
            .flatten(&Atom::new_num(0))
            .unwrap()
            .map_data(|a| a.as_view().to_polynomial_in_vars::<u8>(&var_map));

        PolySplit {
            colorless: colorless_parsed,
            var_map,
            energies,
            color: ParamTensor::composite(color_simplified.state.color.map_data(|a| a.0)),
            colorsimplified: color_simplified.state.state,
        }
    }

    fn to_shadowed_poly_impl(
        poly: &Gloopoly,
        workspace: &Workspace,
        reps: Arc<Mutex<Vec<Atom>>>,
    ) -> Atom {
        if poly.is_zero() {
            return Atom::new_num(0);
        }

        let mut add = Atom::new_num(0);
        let coef = symb!("coef");
        let shift = reps.as_ref().lock().unwrap().len();

        let mut mul_h;
        let mut var_h = workspace.new_atom();
        let mut num_h = workspace.new_atom();
        let mut pow_h = workspace.new_atom();

        for (i, monomial) in poly.into_iter().enumerate() {
            mul_h = Atom::new_num(1);
            for (var_id, &pow) in poly.variables.iter().zip(monomial.exponents) {
                if pow > 0 {
                    match var_id {
                        Variable::Symbol(v) => {
                            var_h.to_var(*v);
                        }
                        Variable::Temporary(_) => {
                            unreachable!("Temporary variable in expression")
                        }
                        Variable::Function(_, a) | Variable::Other(a) => {
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

            mul_h = mul_h * function!(coef, Atom::new_num((i + shift) as i64));
            add = add + mul_h.as_view();
        }

        add
    }
    fn shadow_poly(poly: Gloopoly, reps: Arc<Mutex<Vec<Atom>>>) -> Atom {
        Workspace::get_local().with(|ws| Self::to_shadowed_poly_impl(&poly, ws, reps))
    }

    pub fn optimize(self) -> PolyContracted {
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
    pub fn contract(self) -> Result<Numerator<PolyContracted>> {
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

    fn validate_squared_energies_impl(&self) -> Result<DataTensor<ExpandedIndex>, ()> {
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

    pub fn validate_squared_energies(&self) -> bool {
        self.validate_squared_energies_impl().is_err()
    }
}

#[derive(Debug, Clone, Encode, Decode)]
#[bincode(decode_context = "StateMap")]
pub struct PolyContracted {
    #[bincode(with_serde)]
    pub tensor: ParamTensor,
    pub coef_map: Vec<Atom>,
}

impl Numerator<PolyContracted> {
    pub fn write(
        &self,
        extra_info: &ExtraInfo,
        bare_graph: &BareGraph,
        numerator_format: ExpressionFormat,
    ) -> std::io::Result<()> {
        debug!("Writing poly contracted numerator");
        let file_path = extra_info.path.join("expressions");

        let printer_ops = numerator_format.into();

        let out = (
            self.state
                .coef_map
                .iter()
                .map(|a| {
                    format!(
                        "{}",
                        AtomPrinter::new_with_options(a.as_view(), printer_ops)
                    )
                })
                .collect_vec(),
            format!(
                "{}",
                AtomPrinter::new_with_options(
                    self.state.tensor.iter_flat().next().unwrap().1,
                    printer_ops
                )
            ),
        );

        fs::write(
            file_path.join(format!("poly_contracted_{}_exp.json", bare_graph.name)),
            serde_json::to_string_pretty(&out).unwrap(),
        )?;

        Ok(())
    }

    pub fn to_contracted(self) -> Numerator<Contracted> {
        let coefs: Vec<_> = (0..self.state.coef_map.len())
            .map(|i| function!(GS.coeff, Atom::new_num(i as i64)).to_pattern())
            .collect();

        let coefs_reps: Vec<_> = self.state.coef_map.iter().map(|a| a.to_pattern()).collect();

        let reps: Vec<_> = coefs
            .into_iter()
            .zip(coefs_reps)
            .map(|(p, rhs)| Replacement::new(p, rhs))
            .collect();

        Numerator {
            state: Contracted {
                tensor: self.state.tensor.replace_all_multiple(&reps),
            },
        }
    }

    fn generate_fn_map(&self) -> FunctionMap {
        let mut fn_map = FunctionMap::new();

        for (v, k) in self.state.coef_map.clone().iter().enumerate() {
            fn_map
                .add_tagged_function(
                    GS.coeff,
                    vec![Atom::new_num(v as i64)],
                    format!("coef{v}"),
                    vec![],
                    k.clone(),
                )
                .unwrap();
        }

        Numerator::<Contracted>::add_consts_to_fn_map(&mut fn_map);
        fn_map
    }
    pub fn generate_evaluators(
        self,
        model: &Model,
        graph: &BareGraph,
        extra_info: &ExtraInfo,
        export_settings: &ExportSettings,
    ) -> Numerator<Evaluators> {
        let s = &export_settings.numerator_settings.eval_settings;
        debug!("generating evaluators for contracted numerator");
        let (params, double_param_values, quad_param_values, model_params_start) =
            Contracted::generate_params(graph, model);

        let emr_len = graph.edges.len();

        let fn_map = self.generate_fn_map();

        // let fn_map: FunctionMap = (&owned_fn_map).into();
        let settings = NumeratorEvaluatorSettings {
            options: s,
            inline_asm: export_settings.gammaloop_compile_options.inline_asm(),
            compile_options: export_settings
                .gammaloop_compile_options
                .to_symbolica_compile_options(),
        };

        let single = self.state.evaluator(
            extra_info.path.clone(),
            &graph.name,
            &settings,
            &params,
            &fn_map,
        );

        match s {
            NumeratorEvaluatorOptions::Joint(_) => Numerator {
                state: Evaluators {
                    orientated: Some(
                        single.orientated_joint(graph, &params, extra_info, &settings, &fn_map),
                    ),
                    single,
                    choice: SingleOrCombined::Combined,
                    orientations: extra_info.orientations.clone(),
                    quad_param_values,
                    double_param_values,
                    model_params_start,
                    emr_len,
                    fn_map,
                },
            },
            NumeratorEvaluatorOptions::Iterative(IterativeOptions {
                iterations,
                n_cores,
                verbose,
                ..
            }) => Numerator {
                state: Evaluators {
                    orientated: Some(single.orientated_iterative(
                        graph,
                        &params,
                        extra_info,
                        &settings,
                        *iterations,
                        *n_cores,
                        *verbose,
                        &fn_map,
                    )),
                    single,
                    choice: SingleOrCombined::Combined,
                    orientations: extra_info.orientations.clone(),
                    quad_param_values,
                    double_param_values,
                    model_params_start,
                    emr_len,
                    fn_map,
                },
            },
            _ => Numerator {
                state: Evaluators {
                    orientated: None,
                    single,
                    choice: SingleOrCombined::Single,
                    orientations: extra_info.orientations.clone(),
                    quad_param_values,
                    double_param_values,
                    model_params_start,
                    emr_len,
                    fn_map,
                },
            },
        }

        // else {
        //     let first = self.state.tensor.iter_flat().next().unwrap().1.to_owned();

        //     Err(self)
        // }
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
        Err(eyre!("Only applied feynman rule, simplified color, gamma, parsed into network and contracted, nothing to update"))
    }
}

impl PolyContracted {
    pub fn evaluator(
        self,
        path: PathBuf,
        name: &str,
        eval_options: &NumeratorEvaluatorSettings,
        params: &[Atom],
        fn_map: &FunctionMap,
    ) -> EvaluatorSingle {
        debug!("Generating single evaluator");
        //

        debug!("Generate eval tree set with {} params", params.len());

        let mut eval_tree = self.tensor.to_evaluation_tree(fn_map, params).unwrap();
        debug!("Horner scheme");

        eval_tree.horner_scheme();
        debug!("Common subexpression elimination");
        eval_tree.common_subexpression_elimination();
        debug!("Linearize double");
        let cpe_rounds = eval_options.options.cpe_rounds();
        let eval_double = eval_tree
            .map_coeff::<Complex<F<f64>>, _>(&|r| Complex {
                re: F(r.into()),
                im: F(0.),
            })
            .linearize(cpe_rounds);
        debug!("Linearize quad");

        let eval_quad = eval_tree
            .map_coeff::<Complex<F<f128>>, _>(&|r| Complex {
                re: F(r.into()),
                im: F(f128::new_zero()),
            })
            .linearize(cpe_rounds);

        let eval = eval_tree
            .map_coeff::<F<f64>, _>(&|r| r.into())
            .linearize(cpe_rounds);
        let compiled = if eval_options.options.compile_options().compile() {
            debug!("Compiling iterative evaluator");
            let path = path.join("compiled");
            // let res = std::fs::create_dir_all(&path);
            match std::fs::create_dir(&path) {
                Ok(_) => {}
                Err(e) => match e.kind() {
                    std::io::ErrorKind::AlreadyExists => {}
                    _ => {
                        panic!("Error creating directory: {} at path {}", e, path.display())
                    }
                },
            }
            let mut filename = path.clone();
            filename.push(format!("{}_numerator_single.cpp", name));
            debug!("Compiling  single evaluator to {}", filename.display());
            let filename = filename.to_string_lossy();

            let function_name = format!("{}_numerator_single", name);

            let library_name = path.join(format!("{}_numerator_single.so", name));
            let library_name = library_name.to_string_lossy();

            CompiledEvaluator::new(
                eval.export_cpp(&filename, &function_name, true, eval_options.inline_asm)
                    .unwrap()
                    .compile(&library_name, eval_options.compile_options.clone())
                    .unwrap()
                    .load()
                    .unwrap(),
            )
        } else {
            CompiledEvaluator::default()
        };

        EvaluatorSingle {
            tensor: self.tensor,
            eval_double,
            eval_quad,
            compiled,
            param_len: params.len(),
        }
    }
}

impl GammaSimplified {
    pub fn gamma_symplify_impl(mut expr: SerializableAtom) -> SerializableAtom {
        // let mink = Minkowski::rep(4);
        fn mink(wildcard: Symbol) -> Atom {
            Minkowski::rep(4).pattern(Atom::new_var(wildcard))
        }
        expr.0 = expr.0.expand();
        let pats = [
            Replacement::new(
                Pattern::parse("id(a_,b_)*t_(d___,b_,c___)").unwrap(),
                Pattern::parse("t_(d___,a_,c___)").unwrap(),
            ),
            Replacement::new(
                Pattern::parse("Metric(mink(a_),mink(b_))*t_(d___,mink(b_),c___)").unwrap(),
                Pattern::parse("t_(d___,mink(a_),c___)").unwrap(),
            ),
        ];

        expr = expr.replace_all_multiple_repeat(&pats);
        let pats = vec![
            (
                Pattern::parse("ProjP(a_,b_)").unwrap(),
                Pattern::parse("1/2*id(a_,b_)-1/2*gamma5(a_,b_)").unwrap(),
            ),
            (
                Pattern::parse("ProjM(a_,b_)").unwrap(),
                Pattern::parse("1/2*id(a_,b_)+1/2*gamma5(a_,b_)").unwrap(),
            ),
            (
                Pattern::parse("id(a_,b_)*f_(d___,b_,e___)").unwrap(),
                Pattern::parse("f_(c___,a_,e___)").unwrap(),
            ),
            // (
            //     Pattern::parse("id(aind(a_,b_))*f_(c___,aind(d___,a_,e___))").unwrap(),
            //     Pattern::parse("f_(c___,aind(d___,b_,e___))")
            //         .unwrap()
            //         ,
            // ),
            (
                Pattern::parse("γ(a_,b_,c_)*γ(d_,c_,e_)").unwrap(),
                Pattern::parse("gamma_chain(a_,d_,b_,e_)").unwrap(),
            ),
            (
                Pattern::parse("gamma_chain(a__,b_,c_)*gamma_chain(d__,c_,e_)").unwrap(),
                Pattern::parse("gamma_chain(a__,d__,b_,e_)").unwrap(),
            ),
            (
                Pattern::parse("γ(a_,b_,c_)*gamma_chain(d__,c_,e_)").unwrap(),
                Pattern::parse("gamma_chain(a_,d__,b_,e_)").unwrap(),
            ),
            (
                Pattern::parse("gamma_chain(a__,b_,c_)*γ(d_,c_,e_)").unwrap(),
                Pattern::parse("gamma_chain(a__,d_,b_,e_)").unwrap(),
            ),
            (
                Pattern::parse("gamma_chain(a__,b_,b_)").unwrap(),
                Pattern::parse("gamma_trace(a__)").unwrap(),
            ),
        ];
        let reps: Vec<Replacement> = pats
            .into_iter()
            .map(|(lhs, rhs)| Replacement::new(lhs, rhs))
            .collect();
        expr.0 = expr.0.expand();
        expr.replace_all_multiple_repeat_mut(&reps);
        expr.0 = expr.0.expand();
        expr.replace_all_multiple_repeat_mut(&reps);

        let pat = Pattern::parse("gamma_trace(a__)").unwrap();

        let mut it = expr.0.pattern_match(&pat, None, None);

        let mut max_nargs = 0;

        while let Some(p) = it.next_detailed() {
            for (_, v) in p.match_stack {
                match v {
                    Match::Single(_) => {
                        if max_nargs < 1 {
                            max_nargs = 1;
                        }
                    }
                    Match::Multiple(_, v) => {
                        if max_nargs < v.len() {
                            max_nargs = v.len();
                        }
                    }
                    _ => panic!("should be a single match"),
                }
            }
        }

        expr.0 = expr.0.expand();
        // let pats: Vec<_> = [
        //     (
        //         function!(ETS.id, GS.a_, GS.b_) * function!(GS.f_, GS.d___, GS.b_, GS.c___),
        //         function!(GS.f_, GS.d___, GS.a_, GS.c___),
        //     ),
        //     (
        //         function!(ETS.metric, mink(GS.a_), mink(GS.b_))
        //             * function!(GS.f_, GS.d___, mink(GS.b_), GS.c___),
        //         function!(GS.f_, GS.d___, mink(GS.a_), GS.c___),
        //     ),
        // ]
        // .iter()
        // .map(|(a, b)| Replacement::new(a.to_pattern(), b.to_pattern()))
        // .collect();

        // expr = expr.replace_all_multiple_repeat(&pats);

        let gamma_chain = symb!("gamma_chain");
        let gamma_trace = symb!("gamma_trace");
        let reps: Vec<_> = [
            (
                function!(ETS.id, GS.a_, GS.b_) * function!(GS.f_, GS.d___, GS.b_, GS.c___),
                function!(GS.f_, GS.d___, GS.a_, GS.c___),
            ),
            (
                function!(ETS.metric, mink(GS.a_), mink(GS.b_))
                    * function!(GS.f_, GS.d___, mink(GS.b_), GS.c___),
                function!(GS.f_, GS.d___, mink(GS.a_), GS.c___),
            ),
            (
                function!(UFO.projp, GS.a_, GS.b_),
                (function!(ETS.id, GS.a_, GS.b_) - function!(ETS.gamma5, GS.a_, GS.b_)) / 2,
            ),
            (
                function!(UFO.projm, GS.a_, GS.b_),
                (function!(ETS.id, GS.a_, GS.b_) + function!(ETS.gamma5, GS.a_, GS.b_)) / 2,
            ),
            (
                function!(ETS.gamma, GS.a_, GS.b_, GS.c_)
                    * function!(ETS.gamma, GS.d_, GS.c_, GS.e_),
                function!(gamma_chain, GS.a_, GS.d_, GS.b_, GS.e_),
            ),
            (function!(ETS.gamma, GS.a_, GS.b_, GS.b_), Atom::Zero),
            (
                function!(gamma_chain, GS.a__, GS.a_, GS.b_)
                    * function!(gamma_chain, GS.b__, GS.b_, GS.c_),
                function!(gamma_chain, GS.a__, GS.b__, GS.a_, GS.c_),
            ),
            (
                function!(gamma_chain, GS.a__, GS.a_, GS.b_)
                    * function!(ETS.gamma, GS.y_, GS.b_, GS.c_),
                function!(gamma_chain, GS.a__, GS.y_, GS.a_, GS.c_),
            ),
            (
                function!(ETS.gamma, GS.a_, GS.a_, GS.b_)
                    * function!(gamma_chain, GS.y__, GS.b_, GS.c_),
                function!(gamma_chain, GS.a_, GS.y__, GS.a_, GS.c_),
            ),
        ]
        .iter()
        .map(|(a, b)| Replacement::new(a.to_pattern(), b.to_pattern()))
        .collect();

        expr.0 = expr.0.expand();
        expr.replace_all_multiple_repeat_mut(&reps);
        expr.0 = expr.0.expand();
        expr.replace_all_multiple_repeat_mut(&reps);

        let _pat = function!(gamma_chain, GS.a_, GS.a___, GS.b_, GS.a_).to_pattern();
        // let patodd = (-2 * function!(gamma_chain, GS.a___, GS.b_)).to_pattern();
        // let pateven = (2 * (function!(gamma_chain, GS.b_, GS.a___))).to_pattern();

        // let rhs = PatternOrMap::Map(Box::new(move |m| m.));
        //
        fn gamma_chain_perm(arg: AtomView, _context: &Context, out: &mut Atom) -> bool {
            let gamma_chain = symb!("gamma_chain");
            let mut found = false;
            if let AtomView::Fun(f) = arg {
                if f.get_symbol() == gamma_chain {
                    found = true;
                    let args = f.iter().collect::<Vec<_>>();
                    let len = args.len();
                    if len <= 3 {
                        return false;
                    }

                    if args[len - 3] == args[0] {
                        if len == 4 {
                            *out = function!(ETS.id, args[len - 2], args[len - 1]) * 4;
                            return true;
                        } else if len % 2 == 0 {
                            let mut gcn = FunctionBuilder::new(gamma_chain);
                            let mut gcnp = FunctionBuilder::new(gamma_chain);
                            gcn = gcn.add_arg(args[len - 4]);
                            gcn = gcn.add_args(&args[1..(len - 4)]);
                            for a in &args[1..(len - 4)] {
                                gcnp = gcnp.add_arg(*a);
                            }
                            gcnp = gcnp.add_arg(args[len - 4]);
                            gcnp = gcnp.add_args(&args[(len - 2)..len]);
                            gcn = gcn.add_args(&args[(len - 2)..len]);
                            *out = (gcn.finish() + gcnp.finish()) * 2;
                        } else {
                            let mut gcn = FunctionBuilder::new(gamma_chain);
                            gcn = gcn.add_args(&args[1..(len - 4)]);
                            gcn = gcn.add_args(&args[(len - 2)..len]);
                            *out = gcn.finish() * -2;
                        }

                        // println!("{}->{}", arg, out);
                    } else {
                        return false;
                    }
                }
            }
            found
        }

        loop {
            let new = expr.0.replace_map(&gamma_chain_perm);
            if new == expr.0 {
                break;
            } else {
                expr.0 = new;
            }
        }

        expr.replace_all_repeat_mut(
            &(function!(gamma_chain, GS.a__, GS.x_, GS.x_).to_pattern()),
            function!(gamma_trace, GS.a__).to_pattern(),
            None,
            None,
        );

        // //Chisholm identity:
        // expr.replace_all_repeat_mut(
        //     &(function!(ETS.gamma, GS.a_, GS.x_, GS.y_) * function!(gamma_trace, GS.a_, GS.a__)).to_pattern(),
        //     (function!(gamma_chain, GS.a__)).to_pattern(),
        //     None,
        //     None,
        // );
        //
        fn gamma_tracer(arg: AtomView, _context: &Context, out: &mut Atom) -> bool {
            let gamma_trace = symb!("gamma_trace");

            let mut found = false;
            if let AtomView::Fun(f) = arg {
                if f.get_symbol() == gamma_trace {
                    found = true;
                    let mut sum = Atom::Zero;

                    if f.get_nargs() == 1 {
                        *out = Atom::Zero;
                    }
                    let args = f.iter().collect::<Vec<_>>();

                    for i in 1..args.len() {
                        let sign = if i % 2 == 0 { -1 } else { 1 };

                        let mut gcn = FunctionBuilder::new(gamma_trace);
                        #[allow(clippy::needless_range_loop)]
                        for j in 1..args.len() {
                            if i != j {
                                gcn = gcn.add_arg(args[j]);
                            }
                        }

                        let metric = if args[0] == args[i] {
                            Atom::new_num(4)
                            // Atom::new_var(GS.dim)
                        } else {
                            function!(ETS.metric, args[0], args[i])
                        };
                        if args.len() == 2 {
                            sum = sum + metric * sign * 4;
                        } else {
                            sum = sum + metric * gcn.finish() * sign;
                        }
                    }
                    *out = sum;

                    // println!("{}->{}", arg, out);
                }
            }

            found
        }

        loop {
            let new = expr.0.replace_map(&gamma_tracer);
            if new == expr.0 {
                break;
            } else {
                expr.0 = new;
            }
        }

        let reps: Vec<_> = [
            (
                function!(ETS.metric, mink(GS.a_), mink(GS.b_))
                    * function!(GS.f_, GS.d___, mink(GS.b_), GS.c___),
                function!(GS.f_, GS.d___, mink(GS.a_), GS.c___),
            ),
            (
                function!(ETS.metric, mink(GS.a_), mink(GS.b_)).pow(Atom::new_num(2)),
                Atom::new_num(4),
            ),
            (
                function!(ETS.gamma, GS.a__).pow(Atom::new_num(2)),
                Atom::new_num(16),
            ),
        ]
        .iter()
        .map(|(a, b)| Replacement::new(a.to_pattern(), b.to_pattern()))
        .collect();

        expr.replace_all_multiple_repeat_mut(&reps);
        expr.0 = expr.0.expand();
        expr.replace_all_multiple_repeat_mut(&reps);
        expr
    }

    // pub fn gamma_symplify(color: ColorProjected) -> GammaSimplified {
    //     Self::gamma_symplify_impl(color.colorless)
    // }

    pub fn parse(self) -> Network {
        let net = TensorNetwork::<MixedTensor<f64, AtomStructure>, SerializableAtom>::try_from(
            self.get_single_atom().unwrap().0.as_view(),
        )
        .unwrap()
        .to_fully_parametric()
        .cast();

        // println!("net scalar{}", net.scalar.as_ref().unwrap());
        Network { net }
    }

    pub fn parse_only_colorless(self) -> Network {
        let net = TensorNetwork::<MixedTensor<f64, AtomStructure>, SerializableAtom>::try_from(
            self.colorless
                .clone()
                .scalar()
                .ok_or(NumeratorStateError::Any(eyre!("not a scalar")))
                .unwrap()
                .0
                .as_view(),
        )
        .unwrap()
        .to_fully_parametric()
        .cast();

        // println!("net scalar{}", net.scalar.as_ref().unwrap());
        Network { net }
    }
}

impl Numerator<GammaSimplified> {
    pub fn write(
        &self,
        extra_info: &ExtraInfo,
        bare_graph: &BareGraph,
        numerator_format: ExpressionFormat,
    ) -> std::io::Result<()> {
        debug!("Writing gamma simplified numerator");
        let file_path = extra_info.path.join("expressions");

        let printer_ops = numerator_format.into();

        let out = format!(
            "{}",
            AtomPrinter::new_with_options(self.get_single_atom().unwrap().0.as_view(), printer_ops)
        );

        fs::write(
            file_path.join(format!("gamma_simplified_{}_exp.json", bare_graph.name)),
            serde_json::to_string_pretty(&out).unwrap(),
        )?;

        Ok(())
    }
    pub fn parse(self) -> Numerator<Network> {
        // debug!("GammaSymplified numerator: {}", self.export());
        debug!("Parsing gamma simplified numerator into tensor network");
        Numerator {
            state: self.state.parse(),
        }
    }

    pub fn parse_only_colorless(self) -> Numerator<Network> {
        debug!("Parsing only colorless gamma simplified numerator into tensor network");
        Numerator {
            state: self.state.parse_only_colorless(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
pub struct Network {
    #[bincode(with_serde)]
    pub net: TensorNetwork<ParamTensor<AtomStructure>, SerializableAtom>,
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
pub enum ContractionSettings<R = Rational> {
    Levelled((usize, FunctionMap<R>)),
    Normal,
}

impl Network {
    pub fn parse_impl(expr: AtomView) -> Self {
        let net =
            TensorNetwork::<MixedTensor<f64, AtomStructure>, SerializableAtom>::try_from(expr)
                .unwrap()
                .to_fully_parametric()
                .cast();

        // println!("net scalar{}", net.scalar.as_ref().unwrap());
        Network { net }
    }

    pub fn contract<R>(mut self, settings: ContractionSettings<R>) -> Result<Contracted> {
        match settings {
            ContractionSettings::Levelled((_depth, _fn_map)) => {
                // let levels: Levels<_, _> = self.net.into();
                // levels.contract(depth, fn_map);
                unimplemented!("cannot because of attached dummy lifetime...")
            }
            ContractionSettings::Normal => {
                self.net.contract().unwrap();
                let tensor = self
                    .net
                    .result_tensor_smart()?
                    .map_structure(VecStructure::from);
                // debug!("contracted tensor: {}", tensor);
                Ok(Contracted { tensor })
            }
        }
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

impl<T: Copy + Default> Numerator<SymbolicExpression<T>> {
    pub fn apply_reps(&self, rep_atoms: Vec<(AtomView, AtomView)>) -> Self {
        // println!(
        //     "REPLACEMENTS:\n{}",
        //     rep_atoms
        //         .iter()
        //         .map(|(a, b)| format!("{} -> {}", a, b))
        //         .collect::<Vec<_>>()
        //         .join("\n")
        // );
        // let reps: Vec<Replacement> = rep_atoms
        //     .iter()
        //     .map(|(l, r)| Replacement::new(l.to_pattern(), r.to_pattern()))
        //     .collect::<Vec<_>>();
        let expr = SymbolicExpression::<T> {
            colorless: self.state.colorless.map_data_ref_self(|a| {
                let mut b: Atom = a.0.clone();
                for (src, trgt) in rep_atoms.iter() {
                    b = b.replace_all(&src.to_pattern(), trgt.to_pattern(), None, None);
                }
                //let b = a.0.replace_all_multiple(&reps).into();
                // println!("AFTER: {}", b);
                // Expand each inner level
                b = b.replace_map(&|term, ctx, out| {
                    if ctx.function_level == 0
                        && ctx.parent_type == Some(symbolica::atom::AtomType::Mul)
                    {
                        *out = term.expand();
                        true
                    } else {
                        false
                    }
                });
                b.into()
            }),
            color: self.state.color.map_data_ref_self(|a| {
                let mut b: Atom = a.0.clone();
                // println!("BEFORE: {}", a.0);
                for (src, trgt) in rep_atoms.iter() {
                    b = b.replace_all(&src.to_pattern(), trgt.to_pattern(), None, None);
                }
                b = b.replace_map(&|term, ctx, out| {
                    if ctx.function_level == 0
                        && ctx.parent_type == Some(symbolica::atom::AtomType::Mul)
                    {
                        *out = term.expand();
                        true
                    } else {
                        false
                    }
                });
                //let b = a.0.replace_all_multiple(&reps).into();
                // println!("AFTER: {}", b);
                b.into()
            }),
            state: self.state.state,
        };
        Self { state: { expr } }
    }
}

impl Numerator<Network> {
    pub fn evaluate_with_replacements(
        &self,
        replacements: Vec<(AtomView, AtomView)>,
        fully_numerical_substitutions: bool,
    ) -> Result<Atom, FeynGenError> {
        if !fully_numerical_substitutions {
            let t = self
                .apply_reps(replacements)
                .contract::<Rational>(ContractionSettings::Normal)
                .map_err(|e| FeynGenError::NumeratorEvaluationError(e.to_string()))?;
            if let Some(s) = t.state.tensor.scalar() {
                Ok(s.expand())
            } else {
                Err(FeynGenError::NumeratorEvaluationError(
                    "Could not simplify to a scalar.".into(),
                ))
            }
        } else {
            let g = self.state.net.graph.map_nodes_ref(|(_, d)| {
                d.map_data_ref_result::<_, FeynGenError>(|a| {
                    let mut b = a.clone();
                    for (src, trgt) in replacements.iter() {
                        b = b.replace_all(&src.to_pattern(), trgt.to_pattern(), None, None);
                    }
                    b = b.expand();
                    let mut re = Rational::zero();
                    let mut im = Rational::zero();
                    for (var, coeff) in b.coefficient_list::<u8,_>(&[Atom::new_var(Atom::I)]).iter() {
                        let c = coeff.try_into().map_err(|e| {
                            FeynGenError::NumeratorEvaluationError(format!(
                                "Could not convert tensor coefficient to integer: error: {}, expresssion: {}",
                                e, coeff
                            ))
                        })?;
                        if *var == Atom::new_var(Atom::I) {
                            re = c;
                        } else if *var == Atom::new_num(1) {
                            im = c;
                        } else {
                            return Err(FeynGenError::NumeratorEvaluationError(format!(
                                "Could not convert the following tensor coefficient to a complex integer: {}",
                                b
                            )));
                        }
                    }
                    Ok(Complex::<Rational>::new(re, im))
                })
                .unwrap()
            });

            let mut net = TensorNetwork {
                graph: g,
                scalar: self.state.net.scalar.clone(),
            };

            net.contract().unwrap();

            let (result_tensor, scalar) = net
                .result()
                .map_err(|e| FeynGenError::NumeratorEvaluationError(e.to_string()))?;
            if let Some(s) = result_tensor.clone().scalar() {
                let factor = scalar.unwrap_or(Atom::new_num(1).into()).0;
                let res = (Atom::new_num(s.re) + Atom::new_num(s.im) * Atom::I) * factor;
                Ok(res.expand())
            } else {
                Err(FeynGenError::NumeratorEvaluationError(
                    "Could not simplify to a scalar.".into(),
                ))
            }
        }
    }

    pub fn apply_reps(&self, rep_atoms: Vec<(AtomView, AtomView)>) -> Self {
        let net = TensorNetwork {
            graph: self.state.net.graph.map_nodes_ref(|(_, d)| {
                d.map_data_ref_self(|a| {
                    let mut b = a.clone();
                    for (src, trgt) in rep_atoms.iter() {
                        b = b.replace_all(&src.to_pattern(), trgt.to_pattern(), None, None);
                    }
                    b
                    //let b = a.replace_all_multiple(&reps);
                })
            }),
            scalar: self.state.net.scalar.clone(),
        };

        Self {
            state: Network { net },
        }
    }
    // pub fn random_concretize_reps(&mut self, seed: usize) -> Vec<Replacement> {
    //     let mut prime = PrimeIteratorU64::new(1)
    //         .skip(seed)
    //         .map(|u| Atom::new_num(symbolica::domains::integer::Integer::new(u as i64)));

    //     let mut reps = vec![];

    //     let pat = function!(
    //         GS.f_,
    //         Atom::new_var(GS.y___),
    //         function!(symb!("cind"), Atom::new_var(GS.x_))
    //     )
    //     .to_pattern();

    //     for (n, d) in self.state.net.graph.nodes.iter() {
    //         for (_, a) in d.tensor.iter_flat() {
    //             for m in a.pattern_match(&pat, None, None) {
    //                 if let Atom::Var(f) = m[&GS.f_].to_atom() {
    //                     let mat = function!(
    //                         f.get_symbol(),
    //                         m[&GS.y___].to_atom(),
    //                         function!(symb!("cind"), m[&GS.x_].to_atom())
    //                     );
    //                     // println!("{mat}");

    //                     reps.push(Replacement::new(
    //                         mat.to_pattern(),
    //                         prime.next().unwrap().to_pattern(),
    //                     ));
    //                 } else {
    //                     println!("{}", m[&GS.f_].to_atom());
    //                 }
    //             }
    //         }
    //     }
    //     reps
    // }

    pub fn random_concretize_reps(
        &self,
        sample_iterator: Option<&mut PrimeIteratorU64>,
        fully_numerical_substitution: bool,
    ) -> Vec<(Atom, Atom)> {
        let prime_iterator = if let Some(iterator) = sample_iterator {
            iterator
        } else {
            &mut PrimeIteratorU64::new(1)
        };

        let mut prime = prime_iterator
            .map(|u| Atom::new_num(symbolica::domains::integer::Integer::new(u as i64)));
        let mut reps = vec![];

        if !fully_numerical_substitution {
            let variable = function!(
                GS.f_,
                Atom::new_var(GS.y_),
                function!(symb!("cind"), Atom::new_var(GS.x_))
            );
            let pat = variable.to_pattern();

            for (_n, d) in self.state.net.graph.nodes.iter() {
                for (_, a) in d.tensor.iter_flat() {
                    for m in a.pattern_match(&pat, None, None) {
                        reps.push(pat.replace_wildcards(&m));
                    }
                }
            }
        } else {
            for (_n, d) in self.state.net.graph.nodes.iter() {
                for (_, a) in d.tensor.iter_flat() {
                    let all_symbols = a.get_all_symbols(true);
                    for s in all_symbols {
                        if s == Atom::I {
                            continue;
                        }
                        let mut found_it = false;
                        let pat = function!(s, symb!("x__")).to_pattern();
                        for m in a.pattern_match(&pat, None, None) {
                            reps.push(pat.replace_wildcards(&m));
                            found_it = true;
                        }
                        if !found_it {
                            reps.push(Atom::new_var(s));
                        }
                    }
                }
            }
        }

        reps = reps
            .iter()
            .collect::<HashSet<_>>()
            .iter()
            .map(|&a| a.to_owned())
            .collect::<Vec<_>>();
        reps.sort();
        reps.iter()
            .map(|a: &Atom| (a.clone(), prime.next().unwrap()))
            .collect()
    }

    pub fn contract<R>(self, settings: ContractionSettings<R>) -> Result<Numerator<Contracted>> {
        // debug!(
        //     "contracting network {}",
        //     self.state.net.rich_graph().dot_nodes()
        // );
        let contracted = self.state.contract(settings)?;
        Ok(Numerator { state: contracted })
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
pub struct Contracted {
    #[bincode(with_serde)]
    pub tensor: ParamTensor,
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
    pub fn evaluator(
        self,
        path: PathBuf,
        name: &str,
        eval_options: &NumeratorEvaluatorSettings,
        params: &[Atom],
        fn_map: &FunctionMap,
    ) -> EvaluatorSingle {
        debug!("generating single evaluator");
        // println!("{}", self.tensor);

        debug!("Generate eval tree set with {} params", params.len());

        let mut eval_tree = self.tensor.to_evaluation_tree(fn_map, params).unwrap();
        debug!("Horner scheme");

        eval_tree.horner_scheme();
        debug!("Common subexpression elimination");
        eval_tree.common_subexpression_elimination();
        debug!("Linearize double");
        let cpe_rounds = eval_options.options.cpe_rounds();
        let eval_double = eval_tree
            .map_coeff::<Complex<F<f64>>, _>(&|r| Complex {
                re: F(r.into()),
                im: F(0.),
            })
            .linearize(cpe_rounds);
        debug!("Linearize quad");

        let eval_quad = eval_tree
            .map_coeff::<Complex<F<f128>>, _>(&|r| Complex {
                re: F(r.into()),
                im: F(f128::new_zero()),
            })
            .linearize(cpe_rounds);

        let eval = eval_tree
            .map_coeff::<F<f64>, _>(&|r| r.into())
            .linearize(cpe_rounds);
        let compiled = if eval_options.options.compile_options().compile() {
            debug!("compiling iterative evaluator");
            let path = path.join("compiled");
            // let res = std::fs::create_dir_all(&path);
            match std::fs::create_dir(&path) {
                Ok(_) => {}
                Err(e) => match e.kind() {
                    std::io::ErrorKind::AlreadyExists => {}
                    _ => {
                        panic!("Error creating directory: {} at path {}", e, path.display())
                    }
                },
            }
            let mut filename = path.clone();
            filename.push(format!("{}_numerator_single.cpp", name));
            debug!("Compiling  single evaluator to {}", filename.display());
            let filename = filename.to_string_lossy();

            let function_name = format!("{}_numerator_single", name);

            let library_name = path.join(format!("{}_numerator_single.so", name));
            let library_name = library_name.to_string_lossy();

            CompiledEvaluator::new(
                eval.export_cpp(&filename, &function_name, true, eval_options.inline_asm)
                    .unwrap()
                    .compile(&library_name, eval_options.compile_options.clone())
                    .unwrap()
                    .load()
                    .unwrap(),
            )
        } else {
            CompiledEvaluator::default()
        };

        EvaluatorSingle {
            tensor: self.tensor,
            eval_double,
            eval_quad,
            compiled,
            param_len: params.len(),
        }
    }

    pub fn generate_kinematic_params_impl(
        n_edges: usize,
        pol_data: Vec<(String, i64, usize)>,
    ) -> Vec<Atom> {
        fn atoms_for_pol(name: String, num: i64, size: usize) -> Vec<Atom> {
            let mut data = vec![];
            for index in 0..size {
                let e = FunctionBuilder::new(Symbol::new(&name));
                data.push(
                    e.add_arg(Atom::new_num(num))
                        .add_arg(Atom::parse(&format!("cind({})", index)).unwrap())
                        .finish(),
                );
            }
            data
        }
        let mut params: Vec<Atom> = vec![];

        let mut pols = Vec::new();

        for i in 0..n_edges {
            let named_structure: NamedStructure<String> = NamedStructure::from_iter(
                [PhysReps::new_slot(Lorentz {}.into(), 4, i)],
                "Q".into(),
                Some(i),
            );
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
        params.push(Atom::new_var(Atom::I));

        params
    }
    #[allow(clippy::type_complexity)]
    pub fn generate_params(
        graph: &BareGraph,
        model: &Model,
    ) -> (
        Vec<Atom>,
        Vec<Complex<F<f64>>>,
        Vec<Complex<F<f128>>>,
        usize,
    ) {
        let mut params: Vec<Atom> = vec![];
        let mut param_values: Vec<Complex<F<f64>>> = vec![];
        let mut pols = Vec::new();

        for (i, _) in graph.edges.iter().enumerate() {
            let named_structure: NamedStructure<String> = NamedStructure::from_iter(
                [PhysReps::new_slot(Lorentz {}.into(), 4, i)],
                "Q".into(),
                Some(i),
            );
            params.extend(
                named_structure
                    .to_shell()
                    .expanded_shadow()
                    .unwrap()
                    .data
                    .clone(),
            );
            param_values.extend([Complex::new_zero(); 4]);
        }

        for i in graph
            .vertices
            .iter()
            .filter(|v| matches!(v.vertex_info, VertexInfo::ExternalVertexInfo(_)))
        {
            let pos = graph.get_vertex_position(&i.name).unwrap();
            let vertex_slot = &graph.vertex_slots[pos];
            if let VertexInfo::ExternalVertexInfo(e) = &i.vertex_info {
                pols.extend(e.get_concrete_polarization_atom(pos, vertex_slot));
            }
        }

        param_values.extend(vec![Complex::new_zero(); pols.len()]);
        params.extend(pols);

        param_values.push(Complex::new_i());
        params.push(Atom::new_var(Atom::I));

        let model_params_start = params.len();
        param_values.extend(model.generate_values());
        params.extend(model.generate_params());

        let quad_values = param_values.iter().map(|c| c.map(|f| f.higher())).collect();

        // for p in params.iter() {
        //     println!("Param: {}", p);
        // }

        (params, param_values, quad_values, model_params_start)
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
        Err(eyre!("Only applied feynman rule, simplified color, gamma, parsed into network and contracted, nothing to update"))
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
    pub fn write(
        &self,
        extra_info: &ExtraInfo,
        bare_graph: &BareGraph,
        numerator_format: ExpressionFormat,
    ) -> std::io::Result<()> {
        debug!("Writing contracted numerator");
        let file_path = extra_info.path.join("expressions");

        let printer_ops = numerator_format.into();

        let out = format!(
            "{}",
            AtomPrinter::new_with_options(
                self.state.tensor.iter_flat().next().unwrap().1,
                printer_ops
            )
        );

        fs::write(
            file_path.join(format!("contracted_{}_exp.json", bare_graph.name)),
            serde_json::to_string_pretty(&out).unwrap(),
        )?;

        Ok(())
    }
    fn generate_fn_map(&self) -> FunctionMap {
        let mut map = FunctionMap::new();
        Numerator::<Contracted>::add_consts_to_fn_map(&mut map);
        map
    }
    #[allow(clippy::too_many_arguments)]
    pub fn generate_evaluators_from_params(
        self,
        n_edges: usize,
        name: &str,
        model_params_start: usize,
        params: &[Atom],
        double_param_values: Vec<Complex<F<f64>>>,
        quad_param_values: Vec<Complex<F<f128>>>,
        extra_info: &ExtraInfo,
        export_settings: &ExportSettings,
    ) -> Numerator<Evaluators> {
        let o = &export_settings.numerator_settings.eval_settings;
        let owned_fn_map = self.generate_fn_map();

        let inline_asm = export_settings.gammaloop_compile_options.inline_asm();
        let compile_options = export_settings
            .gammaloop_compile_options
            .to_symbolica_compile_options();

        let eval_settings = NumeratorEvaluatorSettings {
            options: o,
            inline_asm,
            compile_options,
        };

        let fn_map: FunctionMap = owned_fn_map;
        let single = self.state.evaluator(
            extra_info.path.clone(),
            name,
            &eval_settings,
            params,
            &fn_map,
        );

        match o {
            NumeratorEvaluatorOptions::Joint(_) => Numerator {
                state: Evaluators {
                    orientated: Some(single.orientated_joint_impl(
                        n_edges,
                        name,
                        params,
                        extra_info,
                        &eval_settings,
                        &fn_map,
                    )),
                    single,
                    choice: SingleOrCombined::Combined,
                    orientations: extra_info.orientations.clone(),
                    quad_param_values,
                    double_param_values,
                    model_params_start,
                    emr_len: n_edges,
                    fn_map,
                },
            },
            NumeratorEvaluatorOptions::Iterative(IterativeOptions {
                iterations,
                n_cores,
                verbose,
                ..
            }) => Numerator {
                state: Evaluators {
                    orientated: Some(single.orientated_iterative_impl(
                        n_edges,
                        name,
                        params,
                        extra_info,
                        &eval_settings,
                        *iterations,
                        *n_cores,
                        *verbose,
                        &fn_map,
                    )),
                    single,
                    choice: SingleOrCombined::Combined,
                    orientations: extra_info.orientations.clone(),
                    quad_param_values,
                    double_param_values,
                    model_params_start,
                    emr_len: n_edges,
                    fn_map,
                },
            },
            _ => Numerator {
                state: Evaluators {
                    orientated: None,
                    single,
                    choice: SingleOrCombined::Single,
                    orientations: extra_info.orientations.clone(),
                    quad_param_values,
                    double_param_values,
                    model_params_start,
                    emr_len: n_edges,
                    fn_map,
                },
            },
        }
    }

    pub fn generate_evaluators(
        self,
        model: &Model,
        graph: &BareGraph,
        extra_info: &ExtraInfo,
        export_settings: &ExportSettings,
    ) -> Numerator<Evaluators> {
        let o = &export_settings.numerator_settings.eval_settings;
        debug!("generating evaluators for contracted numerator");
        let (params, double_param_values, quad_param_values, model_params_start) =
            Contracted::generate_params(graph, model);

        trace!("params length:{}", params.len());

        let emr_len = graph.edges.len();

        let fn_map = self.generate_fn_map();

        let settings = NumeratorEvaluatorSettings {
            options: o,
            inline_asm: export_settings.gammaloop_compile_options.inline_asm(),
            compile_options: export_settings
                .gammaloop_compile_options
                .to_symbolica_compile_options(),
        };

        let single = self.state.evaluator(
            extra_info.path.clone(),
            &graph.name,
            &settings,
            &params,
            &fn_map,
        );

        match o {
            NumeratorEvaluatorOptions::Joint(_) => Numerator {
                state: Evaluators {
                    orientated: Some(
                        single.orientated_joint(graph, &params, extra_info, &settings, &fn_map),
                    ),
                    single,
                    choice: SingleOrCombined::Combined,
                    orientations: extra_info.orientations.clone(),
                    quad_param_values,
                    double_param_values,
                    model_params_start,
                    emr_len,
                    fn_map,
                },
            },
            NumeratorEvaluatorOptions::Iterative(IterativeOptions {
                iterations,
                n_cores,
                verbose,
                ..
            }) => Numerator {
                state: Evaluators {
                    orientated: Some(single.orientated_iterative(
                        graph,
                        &params,
                        extra_info,
                        &settings,
                        *iterations,
                        *n_cores,
                        *verbose,
                        &fn_map,
                    )),
                    single,
                    choice: SingleOrCombined::Combined,
                    orientations: extra_info.orientations.clone(),
                    double_param_values,
                    quad_param_values,
                    model_params_start,
                    emr_len,
                    fn_map,
                },
            },
            _ => Numerator {
                state: Evaluators {
                    orientated: None,
                    single,
                    choice: SingleOrCombined::Single,
                    orientations: extra_info.orientations.clone(),
                    double_param_values,
                    quad_param_values,
                    model_params_start,
                    emr_len,
                    fn_map,
                },
            },
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
pub struct IterativeOptions {
    pub eval_options: EvaluatorOptions,
    pub iterations: usize,
    pub n_cores: usize,
    pub verbose: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
pub enum NumeratorParseMode {
    Polynomial,
    Direct,
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
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
    pub fn compile_options(&self) -> NumeratorCompileOptions {
        match self {
            NumeratorEvaluatorOptions::Single(options) => options.compile_options,
            NumeratorEvaluatorOptions::Joint(options) => options.compile_options,
            NumeratorEvaluatorOptions::Iterative(options) => options.eval_options.compile_options,
        }
    }

    pub fn cpe_rounds(&self) -> Option<usize> {
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
    pub fn compile(&self) -> bool {
        matches!(self, NumeratorCompileOptions::Compiled)
    }
}

#[derive(Clone, Serialize, Deserialize)]
pub struct EvaluatorOrientations {
    pub eval_double: LinearizedEvalTensorSet<Complex<F<f64>>, VecStructure>,
    pub eval_quad: LinearizedEvalTensorSet<Complex<F<f128>>, VecStructure>,
    pub compiled: CompiledEvaluator<EvalTensorSet<SerializableCompiledEvaluator, VecStructure>>,
    pub positions: Vec<usize>,
}

impl Debug for EvaluatorOrientations {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("EvaluatorOrientations")
            .field(
                "eval_double",
                &"EvalTensorSet<ExpressionEvaluator<Complex<F<f64>>>, AtomStructure>",
            )
            .field(
                "eval_quad",
                &"EvalTensorSet<ExpressionEvaluator<Complex<F<f128>>>, AtomStructure>",
            )
            .field("compiled", &self.compiled.is_compiled())
            .finish()
    }
}

#[derive(Clone, Serialize, Deserialize)]
pub struct EvaluatorSingle {
    tensor: ParamTensor,
    eval_double: EvalTensor<ExpressionEvaluator<Complex<F<f64>>>, VecStructure>,
    eval_quad: EvalTensor<ExpressionEvaluator<Complex<F<f128>>>, VecStructure>,
    compiled: CompiledEvaluator<EvalTensor<SerializableCompiledEvaluator, VecStructure>>,
    param_len: usize, //to validate length of params at run time
}

impl EvaluatorSingle {
    #[allow(clippy::too_many_arguments)]
    pub fn orientated_iterative_impl(
        &self,
        n_edges: usize,
        name: &str,
        params: &[Atom],
        extra_info: &ExtraInfo,
        eval_options: &NumeratorEvaluatorSettings,
        iterations: usize,
        n_cores: usize,
        verbose: bool,
        fn_map: &FunctionMap,
    ) -> EvaluatorOrientations {
        let mut seen = 0;

        let mut index_map = IndexSet::new();
        let mut positions = Vec::new();

        let len = extra_info.orientations.len();

        let reps = (0..n_edges)
            .map(|i| {
                (
                    Pattern::parse(&format!("Q({},cind(0))", i)).unwrap(),
                    Pattern::parse(&format!("-Q({},cind(0))", i)).unwrap(),
                )
            })
            .collect_vec();

        let settings = MatchSettings {
            rhs_cache_size: 0,
            ..Default::default()
        };

        for (ni, o) in extra_info.orientations.iter().enumerate() {
            let time = Instant::now();
            let elapsed_parse = time.elapsed();

            let reps = reps
                .iter()
                .enumerate()
                .filter_map(|(i, (lhs, rhs))| {
                    if o[i] {
                        None
                    } else {
                        Some(
                            Replacement::new(lhs.clone(), rhs.clone())
                                .with_settings(settings.clone()),
                        )
                    }
                })
                .collect_vec();

            let time = Instant::now();
            let orientation_replaced_net = self.tensor.replace_all_multiple(&reps);
            let elapsed = time.elapsed();

            let time = Instant::now();
            let (entry, is_new) = index_map.insert_full(orientation_replaced_net);
            let hash_time = time.elapsed();
            if !is_new {
                seen += 1;
                debug!("Recycled orientation");
            }
            debug!(
                "Contracted an orientation {:.1}%, parse:{:?},reps:{:?},hash:{:?}",
                100. * (ni as f64) / (len as f64),
                elapsed_parse,
                elapsed,
                hash_time,
            );
            positions.push(entry);
        }

        debug!(
            "Recycled orientations: {:.1}%",
            100. * (seen as f64) / (len as f64)
        );

        let new_index_map = index_map.split_off(1);
        let init = vec![index_map.pop().unwrap()];

        let set = ParamTensorSet::new(init);

        debug!("{} tensors in set", set.tensors.len());

        let cpe_rounds = eval_options.options.cpe_rounds();

        debug!("Generate eval tree set with {} params", params.len());

        let mut eval_tree = set.eval_tree(fn_map, params).unwrap();

        debug!("{} tensors in eval_tree", eval_tree.len());
        debug!("Horner scheme");

        eval_tree.horner_scheme();
        debug!("Common subexpression elimination");
        eval_tree.common_subexpression_elimination();
        debug!("Linearize double");

        let mut linearized = eval_tree.linearize(cpe_rounds);

        debug!("{} linearized tensors", linearized.len());

        for (i, t) in new_index_map.iter().enumerate() {
            let eval_tree = t.to_evaluation_tree(fn_map, params).unwrap();
            debug!("Push optimizing :{}", i + 1);
            linearized.push_optimize(eval_tree, cpe_rounds, iterations, n_cores, verbose);
        }

        let eval_double = linearized
            .clone()
            .map_coeff::<Complex<F<f64>>, _>(&|r| Complex {
                re: F(r.into()),
                im: F(0.),
            });
        debug!("Linearize quad");

        let eval_quad = linearized
            .clone()
            .map_coeff::<Complex<F<f128>>, _>(&|r| Complex {
                re: F(r.into()),
                im: F(f128::new_zero()),
            });

        let eval = linearized.clone().map_coeff::<F<f64>, _>(&|r| r.into());

        let compiled = if eval_options.options.compile_options().compile() {
            debug!("compiling iterative evaluator");
            let path = extra_info.path.join("compiled");
            // let res = std::fs::create_dir_all(&path);
            match std::fs::create_dir(&path) {
                Ok(_) => {}
                Err(e) => match e.kind() {
                    std::io::ErrorKind::AlreadyExists => {}
                    _ => {
                        panic!("Error creating directory: {}", e)
                    }
                },
            }

            let mut filename = path.clone();

            filename.push(format!("{}_numerator_iterative.cpp", name));
            let filename = filename.to_string_lossy();

            let function_name = format!("{}_numerator_iterative", name);

            let library_name = path.join(format!("{}_numerator_iterative.so", name));
            let library_name = library_name.to_string_lossy();
            CompiledEvaluator::new(
                eval.export_cpp(&filename, &function_name, true, eval_options.inline_asm)
                    .unwrap()
                    .compile(&library_name, eval_options.compile_options.clone())
                    .unwrap()
                    .load()
                    .unwrap(),
            )
        } else {
            CompiledEvaluator::default()
        };

        EvaluatorOrientations {
            positions,
            eval_double,
            eval_quad,
            compiled,
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub fn orientated_iterative(
        &self,
        graph: &BareGraph,
        params: &[Atom],
        extra_info: &ExtraInfo,
        export_settings: &NumeratorEvaluatorSettings,
        iterations: usize,
        n_cores: usize,
        verbose: bool,
        fn_map: &FunctionMap,
    ) -> EvaluatorOrientations {
        debug!("generate iterative evaluator");
        self.orientated_iterative_impl(
            graph.edges.len(),
            &graph.name,
            params,
            extra_info,
            export_settings,
            iterations,
            n_cores,
            verbose,
            fn_map,
        )
    }

    pub fn orientated_joint_impl(
        &self,
        n_edges: usize,
        name: &str,
        params: &[Atom],
        extra_info: &ExtraInfo,
        eval_options: &NumeratorEvaluatorSettings,
        fn_map: &FunctionMap,
    ) -> EvaluatorOrientations {
        let mut seen = 0;

        let mut index_map = IndexSet::new();
        let mut positions = Vec::new();

        let len = extra_info.orientations.len();

        let reps = (0..n_edges)
            .map(|i| {
                (
                    Pattern::parse(&format!("Q({},cind(0))", i)).unwrap(),
                    Pattern::parse(&format!("-Q({},cind(0))", i)).unwrap(),
                )
            })
            .collect_vec();

        let settings = MatchSettings {
            rhs_cache_size: 0,
            ..Default::default()
        };

        for (ni, o) in extra_info.orientations.iter().enumerate() {
            let time = Instant::now();
            let elapsed_parse = time.elapsed();

            let reps = reps
                .iter()
                .enumerate()
                .filter_map(|(i, (lhs, rhs))| {
                    if o[i] {
                        None
                    } else {
                        Some(
                            Replacement::new(lhs.clone(), rhs.clone())
                                .with_settings(settings.clone()),
                        )
                    }
                })
                .collect_vec();

            let time = Instant::now();
            let orientation_replaced_net = self.tensor.replace_all_multiple(&reps);
            let elapsed = time.elapsed();

            let time = Instant::now();
            let (entry, is_new) = index_map.insert_full(orientation_replaced_net);
            let hash_time = time.elapsed();
            if !is_new {
                seen += 1;
            }
            debug!(
                "Contracted an orientation {:.1}%, parse:{:?},reps:{:?},hash:{:?}",
                100. * (ni as f64) / (len as f64),
                elapsed_parse,
                elapsed,
                hash_time,
            );
            positions.push(entry);
        }

        debug!(
            "Recycled orientations: {:.1}%",
            100. * (seen as f64) / (len as f64)
        );

        let set = ParamTensorSet::new(index_map.into_iter().collect_vec());

        let cpe_rounds = eval_options.options.cpe_rounds();

        debug!("Generate eval tree set with {} params", params.len());

        let mut eval_tree = set.eval_tree(fn_map, params).unwrap();
        debug!("Horner scheme");

        eval_tree.horner_scheme();
        debug!("Common subexpression elimination");
        eval_tree.common_subexpression_elimination();
        debug!("Linearize double");
        let eval_double = eval_tree
            .map_coeff::<Complex<F<f64>>, _>(&|r| Complex {
                re: F(r.into()),
                im: F(0.),
            })
            .linearize(cpe_rounds);
        debug!("Linearize quad");

        let eval_quad = eval_tree
            .map_coeff::<Complex<F<f128>>, _>(&|r| Complex {
                re: F(r.into()),
                im: F(f128::new_zero()),
            })
            .linearize(cpe_rounds);

        let eval = eval_tree
            .map_coeff::<F<f64>, _>(&|r| r.into())
            .linearize(cpe_rounds);

        let compiled = if eval_options.options.compile_options().compile() {
            debug!("Compiling joint evaluator");
            let path = extra_info.path.join("compiled");
            // let res = std::fs::create_dir_all(&path);
            match std::fs::create_dir(&path) {
                Ok(_) => {}
                Err(e) => match e.kind() {
                    std::io::ErrorKind::AlreadyExists => {}
                    _ => {
                        panic!("Error creating directory: {}", e)
                    }
                },
            }

            let mut filename = path.clone();

            filename.push(format!("{}_numerator_joint.cpp", name));
            let filename = filename.to_string_lossy();

            let function_name = format!("{}_numerator_joint", name);

            let library_name = path.join(format!("{}_numerator_joint.so", name));
            let library_name = library_name.to_string_lossy();

            CompiledEvaluator::new(
                eval.export_cpp(&filename, &function_name, true, eval_options.inline_asm)
                    .unwrap()
                    .compile(&library_name, eval_options.compile_options.clone())
                    .unwrap()
                    .load()
                    .unwrap(),
            )
        } else {
            CompiledEvaluator::default()
        };

        EvaluatorOrientations {
            positions,
            eval_double,
            eval_quad,
            compiled,
        }
    }

    pub fn orientated_joint(
        &self,
        graph: &BareGraph,
        params: &[Atom],
        extra_info: &ExtraInfo,
        eval_options: &NumeratorEvaluatorSettings,
        fn_map: &FunctionMap,
    ) -> EvaluatorOrientations {
        debug!("Generate joint evaluator");
        self.orientated_joint_impl(
            graph.edges.len(),
            &graph.name,
            params,
            extra_info,
            eval_options,
            fn_map,
        )
    }
}

impl Debug for EvaluatorSingle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("EvaluatorSingle")
            .field(
                "eval_double",
                &"EvalTensor<ExpressionEvaluator<Complex<F<f64>>>, AtomStructure>",
            )
            .field(
                "eval_quad",
                &"EvalTensor<ExpressionEvaluator<Complex<F<f128>>>, AtomStructure>",
            )
            .field("compiled", &self.compiled.is_compiled())
            .finish()
    }
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct CompiledEvaluator<E> {
    pub state: CompiledState,
    pub evaluator: Option<E>,
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub enum CompiledState {
    Enabled,
    Disabled,
}

impl<E> Default for CompiledEvaluator<E> {
    fn default() -> Self {
        CompiledEvaluator {
            state: CompiledState::Disabled,
            evaluator: None,
        }
    }
}

impl<E> CompiledEvaluator<E> {
    pub fn new(evaluator: E) -> Self {
        CompiledEvaluator {
            state: CompiledState::Enabled,
            evaluator: Some(evaluator),
        }
    }
    pub fn disable(&mut self) -> Result<()> {
        if self.evaluator.is_some() {
            self.state = CompiledState::Disabled;
            Ok(())
        } else {
            Err(eyre!("Cannot disable evaluator that is not compiled"))
        }
    }

    pub fn enable(&mut self) -> Result<()> {
        if self.evaluator.is_some() {
            self.state = CompiledState::Enabled;
            Ok(())
        } else {
            Err(eyre!("Cannot enable evaluator that is not compiled"))
        }
    }

    pub fn is_compiled(&self) -> bool {
        self.evaluator.is_some()
    }

    pub fn is_enabled(&self) -> bool {
        matches!(self.state, CompiledState::Enabled)
    }
}
use symbolica::state::StateMap;

#[derive(Clone, Debug, Encode, Decode)]
#[bincode(decode_context = "StateMap")]
pub struct Evaluators {
    #[bincode(with_serde)]
    orientated: Option<EvaluatorOrientations>,
    #[bincode(with_serde)]
    pub single: EvaluatorSingle,
    choice: SingleOrCombined,
    orientations: Vec<Vec<bool>>,
    #[bincode(with_serde)]
    pub double_param_values: Vec<Complex<F<f64>>>,
    #[bincode(with_serde)]
    quad_param_values: Vec<Complex<F<f128>>>,
    model_params_start: usize,
    emr_len: usize,
    fn_map: FunctionMap,
}

#[derive(Clone, Serialize, Deserialize, Debug, Encode, Decode)]
pub enum SingleOrCombined {
    Single,
    Combined,
}

impl TryFrom<PythonState> for Evaluators {
    type Error = NumeratorStateError;

    fn try_from(value: PythonState) -> std::result::Result<Self, Self::Error> {
        match value {
            PythonState::Evaluators(s) => {
                if let Some(s) = s {
                    Ok(s)
                } else {
                    Err(NumeratorStateError::NoneVariant)
                }
            }
            _ => Err(NumeratorStateError::NotEvaluators),
        }
    }
}

impl NumeratorState for Evaluators {
    fn export(&self) -> String {
        "evaluators".to_string()
    }

    fn forget_type(self) -> PythonState {
        PythonState::Evaluators(Some(self))
    }

    fn update_model(&mut self, model: &Model) -> Result<()> {
        for (i, &v) in model.generate_values().iter().enumerate() {
            self.double_param_values[self.model_params_start + i] = v.map(|f| f);
            self.quad_param_values[self.model_params_start + i] = v.map(|f| f.higher());
        }
        Ok(())
    }
}

impl TypedNumeratorState for Evaluators {
    fn apply<F, S: TypedNumeratorState>(
        num: &mut Numerator<PythonState>,
        mut f: F,
    ) -> Result<(), NumeratorStateError>
    where
        F: FnMut(Numerator<Self>) -> Numerator<S>,
    {
        if let PythonState::Evaluators(s) = &mut num.state {
            if let Some(s) = s.take() {
                *num = f(Numerator { state: s }).forget_type();
                return Ok(());
            } else {
                return Err(NumeratorStateError::NoneVariant);
            }
        }
        Err(NumeratorStateError::NotEvaluators)
    }
}

impl Numerator<Evaluators> {
    pub fn disable_compiled(&mut self) {
        if let Some(orientated) = &mut self.state.orientated {
            orientated.compiled.disable().unwrap();
        }
        self.state.single.compiled.disable().unwrap();
    }

    pub fn disable_combined(&mut self) {
        if self.state.orientated.is_some() {
            self.state.choice = SingleOrCombined::Single;
        }
    }

    pub fn enable_compiled(&mut self) {
        if let Some(orientated) = &mut self.state.orientated {
            orientated.compiled.enable().unwrap();
        }
        self.state.single.compiled.enable().unwrap();
    }

    pub fn enable_combined(
        &mut self,
        generate: Option<(&Model, &BareGraph, &ExtraInfo, &NumeratorEvaluatorSettings)>,
    ) {
        if self.state.orientated.is_some() {
            self.state.choice = SingleOrCombined::Combined;
        } else if let Some((model, graph, extra_info, export_settings)) = generate {
            let (params, _, _, _) = Contracted::generate_params(graph, model);
            let orientated = self.state.single.orientated_joint(
                graph,
                &params,
                extra_info,
                export_settings,
                &self.state.fn_map,
            );
            self.state.orientated = Some(orientated);
            self.state.choice = SingleOrCombined::Combined;
        } else {
            panic!("Cannot enable combined without generating the evaluators")
        }
    }
}

use thiserror::Error;
use uuid::Uuid;

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

#[derive(Debug, Clone, Encode, Decode)]
#[bincode(decode_context = "StateMap")]
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
    Evaluators(Option<Evaluators>),
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
            PythonState::Evaluators(state) => {
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

            PythonState::Evaluators(state) => {
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
    fn get_single_atom(&self) -> Result<SerializableAtom, NumeratorStateError> {
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
#[cfg(test)]
mod tests;
