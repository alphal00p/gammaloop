use std::fmt::Debug;
use std::path::PathBuf;
use std::time::Instant;

use crate::debug_info::DEBUG_LOGGER;
use crate::graph::BareGraph;
use crate::momentum::Polarization;
use crate::utils::f128;
use crate::{
    graph::EdgeType,
    model::Model,
    momentum::FourMomentum,
    utils::{FloatLike, F},
};
use crate::{ExportSettings, Settings};
use bincode::{Decode, Encode};
use color_eyre::{Report, Result};
use eyre::eyre;
use gat_lending_iterator::LendingIterator;
// use gxhash::GxBuildHasher;
use indexmap::IndexSet;
use itertools::Itertools;

use log::{debug, trace};

use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};
use spenso::data::DataTensor;

use spenso::network::Levels;
use spenso::parametric::{
    EvalTensor, EvalTensorSet, LinearizedEvalTensorSet, ParamTensorSet, SerializableAtom,
    SerializableCompiledEvaluator, TensorSet,
};
use spenso::structure::{HasStructure, SerializableSymbol, SmartShadowStructure};
use spenso::{
    complex::Complex,
    network::TensorNetwork,
    parametric::{ParamTensor, PatternReplacement},
    structure::{Lorentz, NamedStructure, PhysReps, RepName, Shadowable, TensorStructure},
};

use symbolica::atom::AtomView;
use symbolica::evaluate::ExpressionEvaluator;
use symbolica::id::{Condition, Match, MatchSettings};
use symbolica::{
    atom::{Atom, FunctionBuilder},
    state::State,
};
use symbolica::{
    domains::{float::NumericalFloatLike, rational::Rational},
    evaluate::FunctionMap,
    id::{Pattern, Replacement},
};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NumeratorSettings {
    pub eval_settings: NumeratorEvaluatorOptions,
    pub global_numerator: Option<String>,
    pub gamma_algebra: GammaAlgebraMode,
}

impl Default for NumeratorSettings {
    fn default() -> Self {
        NumeratorSettings {
            eval_settings: Default::default(),
            global_numerator: None,
            gamma_algebra: GammaAlgebraMode::Symbolic,
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
    ) -> Result<RepeatingIteratorTensorOrScalar<DataTensor<Complex<F<T>>, AtomStructure>>>;

    fn evaluate_single(
        &mut self,
        emr: &[FourMomentum<F<T>>],
        polarizations: &[Polarization<Complex<F<T>>>],
        tag: Option<Uuid>,
        setting: &Settings,
    ) -> DataTensor<Complex<F<T>>, AtomStructure>;
}

impl<T: FloatLike> Evaluate<T> for Numerator<Evaluators> {
    fn evaluate_all_orientations(
        &mut self,
        emr: &[FourMomentum<F<T>>],
        polarizations: &[Polarization<Complex<F<T>>>],
        _tag: Option<Uuid>,
        settings: &Settings,
    ) -> Result<RepeatingIteratorTensorOrScalar<DataTensor<Complex<F<T>>, AtomStructure>>> {
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
    ) -> DataTensor<Complex<F<T>>, AtomStructure> {
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
    ) -> Result<RepeatingIteratorTensorOrScalar<DataTensor<Complex<F<T>>, AtomStructure>>>;

    fn evaluate_single(num: &mut Numerator<Evaluators>)
        -> DataTensor<Complex<F<T>>, AtomStructure>;

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
    type Item<'a> = &'a T where Self:'a ;

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
    ) -> Result<RepeatingIteratorTensorOrScalar<DataTensor<Complex<F<Self>>, AtomStructure>>> {
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

    fn evaluate_single(
        num: &mut Numerator<Evaluators>,
    ) -> DataTensor<Complex<F<Self>>, AtomStructure> {
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
    ) -> Result<RepeatingIteratorTensorOrScalar<DataTensor<Complex<F<Self>>, AtomStructure>>> {
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

    fn evaluate_single(
        num: &mut Numerator<Evaluators>,
    ) -> DataTensor<Complex<F<Self>>, AtomStructure> {
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
        fn_map.add_constant(Atom::parse("Nc").unwrap(), 3.into());

        fn_map.add_constant(Atom::parse("TR").unwrap(), Rational::from((1, 2)));

        fn_map.add_constant(
            Atom::parse("pi").unwrap(),
            Rational::from(std::f64::consts::PI),
        );
    }
}

impl<S: UnexpandedNumerator> Numerator<S> {
    pub fn expr(&self) -> Result<&SerializableAtom, NumeratorStateError> {
        self.state.expr()
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
pub trait NumeratorState: Serialize + Clone + DeserializeOwned + Debug + Encode + Decode {
    fn export(&self) -> String;

    fn forget_type(self) -> PythonState;

    fn update_model(&mut self, model: &Model) -> Result<()>;
    // fn try_from(state: PythonState) -> Result<Self>;
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
pub struct UnInit;

impl Default for UnInit {
    fn default() -> Self {
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

impl Numerator<UnInit> {
    pub fn from_graph(self, graph: &mut BareGraph) -> Numerator<AppliedFeynmanRule> {
        debug!("applying feynman rules ");
        Numerator {
            state: AppliedFeynmanRule::from_graph(graph),
        }
    }

    pub fn from_global(self, global: Atom, _graph: &BareGraph) -> Numerator<Global> {
        debug!("setting global numerator");
        Numerator {
            state: Global::new(global.into()),
        }
    }
}

impl Default for Numerator<UnInit> {
    fn default() -> Self {
        Numerator { state: UnInit }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
pub struct SingleExpression<State> {
    #[bincode(with_serde)]
    pub expression: SerializableAtom,
    pub state: State,
}

pub trait UnexpandedNumerator: NumeratorState {
    fn expr(&self) -> Result<&SerializableAtom, NumeratorStateError>;
}

impl<E: ExpressionState> UnexpandedNumerator for SingleExpression<E> {
    fn expr(&self) -> Result<&SerializableAtom, NumeratorStateError> {
        Ok(&self.expression)
    }
}

impl<E: ExpressionState> SingleExpression<E> {
    pub fn new(expression: SerializableAtom) -> Self {
        E::new(expression)
    }
}

pub trait ExpressionState:
    Serialize + Clone + DeserializeOwned + Debug + Encode + Decode + Default
{
    fn forget_type(self, expression: SerializableAtom) -> PythonState;

    fn new(expression: SerializableAtom) -> SingleExpression<Self> {
        SingleExpression {
            expression,
            state: Self::default(),
        }
    }

    fn get_expression(num: &mut PythonState)
        -> Result<SingleExpression<Self>, NumeratorStateError>;
}

impl<E: ExpressionState> TypedNumeratorState for SingleExpression<E> {
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

impl<E: ExpressionState> TryFrom<&mut PythonState> for SingleExpression<E> {
    type Error = NumeratorStateError;

    fn try_from(value: &mut PythonState) -> std::result::Result<Self, Self::Error> {
        let a = E::get_expression(value)?;
        Ok(a)
    }
}

impl<E: ExpressionState> TryFrom<PythonState> for SingleExpression<E> {
    type Error = NumeratorStateError;

    fn try_from(mut value: PythonState) -> std::result::Result<Self, Self::Error> {
        let a = E::get_expression(&mut value)?;
        Ok(a)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, Default)]
pub struct Local {}
pub type AppliedFeynmanRule = SingleExpression<Local>;

impl ExpressionState for Local {
    fn forget_type(self, expression: SerializableAtom) -> PythonState {
        PythonState::AppliedFeynmanRule(Some(AppliedFeynmanRule {
            expression,
            state: self,
        }))
    }

    fn new(expression: SerializableAtom) -> SingleExpression<Self> {
        SingleExpression {
            expression,
            state: Local {},
        }
    }

    fn get_expression(
        num: &mut PythonState,
    ) -> Result<SingleExpression<Self>, NumeratorStateError> {
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

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, Default)]
pub struct NonLocal {}
pub type Global = SingleExpression<NonLocal>;

impl ExpressionState for NonLocal {
    fn forget_type(self, expression: SerializableAtom) -> PythonState {
        PythonState::Global(Some(Global {
            expression,
            state: self,
        }))
    }

    fn get_expression(
        num: &mut PythonState,
    ) -> Result<SingleExpression<Self>, NumeratorStateError> {
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
    pub fn color_symplify(self) -> Numerator<ColorSymplified> {
        debug!("color symplifying global numerator");
        Numerator {
            state: ColorSymplified::color_symplify(self.state),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, Default)]
pub struct Color {}
pub type ColorSymplified = SingleExpression<Color>;

impl ExpressionState for Color {
    fn forget_type(self, expression: SerializableAtom) -> PythonState {
        PythonState::ColorSymplified(Some(ColorSymplified {
            expression,
            state: self,
        }))
    }

    fn get_expression(
        num: &mut PythonState,
    ) -> Result<SingleExpression<Self>, NumeratorStateError> {
        if let PythonState::ColorSymplified(s) = num {
            if let Some(s) = s.take() {
                Ok(s)
            } else {
                Err(NumeratorStateError::NoneVariant)
            }
        } else {
            Err(NumeratorStateError::NotColorSymplified)
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, Default)]
pub struct Projected {}
pub type ColorProjected = SingleExpression<Projected>;

impl ExpressionState for Projected {
    fn forget_type(self, expression: SerializableAtom) -> PythonState {
        PythonState::ColorProjected(Some(ColorProjected {
            expression,
            state: self,
        }))
    }

    fn get_expression(
        num: &mut PythonState,
    ) -> Result<SingleExpression<Self>, NumeratorStateError> {
        if let PythonState::ColorProjected(s) = num {
            if let Some(s) = s.take() {
                Ok(s)
            } else {
                Err(NumeratorStateError::NoneVariant)
            }
        } else {
            Err(NumeratorStateError::NotColorProjected)
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, Default)]
pub struct Gamma {}
pub type GammaSymplified = SingleExpression<Gamma>;

impl ExpressionState for Gamma {
    fn forget_type(self, expression: SerializableAtom) -> PythonState {
        PythonState::GammaSymplified(Some(GammaSymplified {
            expression,
            state: self,
        }))
    }

    fn get_expression(
        num: &mut PythonState,
    ) -> Result<SingleExpression<Self>, NumeratorStateError> {
        if let PythonState::GammaSymplified(s) = num {
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

impl<State: ExpressionState> NumeratorState for SingleExpression<State> {
    fn export(&self) -> String {
        self.expression.to_string()
    }

    fn forget_type(self) -> PythonState {
        self.state.forget_type(self.expression)
    }

    fn update_model(&mut self, _model: &Model) -> Result<()> {
        Err(eyre!("Only an expression, nothing to update"))
    }
}

impl AppliedFeynmanRule {
    pub fn from_graph(graph: &mut BareGraph) -> Self {
        debug!("generating numerator for graph: {}", graph.name);
        debug!("momentum: {}", graph.dot_lmb());

        let vatoms: Vec<Atom> = graph
            .vertices
            .iter()
            .flat_map(|v| v.contracted_vertex_rule(graph))
            .collect();

        let mut eatoms: Vec<Atom> = vec![];
        let mut shift = 0;
        let i = Atom::new_var(State::I);
        for e in &graph.edges {
            let (n, s) = e.numerator(graph);
            let n = if matches!(e.edge_type, EdgeType::Virtual) {
                &n * &i
            } else {
                n
            };
            eatoms.push(n);
            shift += s;
            graph.shifts.0 += shift;
        }
        let mut builder = Atom::new_num(1);

        for v in &vatoms {
            // println!("vertex: {v}");
            builder = builder * v;
        }

        for e in &eatoms {
            builder = builder * e;
        }

        let a = Atom::new_var(State::get_symbol("a_"));
        let b = Atom::new_var(State::get_symbol("b_"));

        let complex = FunctionBuilder::new(State::get_symbol("complex"))
            .add_arg(&a)
            .add_arg(&b)
            .finish();

        builder = complex.into_pattern().replace_all(
            builder.as_view(),
            &(&a + &b * &i).into_pattern().into(),
            None,
            None,
        );

        AppliedFeynmanRule::new(builder.into())
    }
}

impl Numerator<AppliedFeynmanRule> {
    pub fn color_symplify(self) -> Numerator<ColorSymplified> {
        debug!("Applied feynman rules: {}", self.export());
        debug!("color symplifying local numerator");

        Numerator {
            state: ColorSymplified::color_symplify(self.state),
        }
    }
}

impl ColorSymplified {
    fn isolate_color(expression: &mut SerializableAtom) {
        let color_fn = FunctionBuilder::new(State::get_symbol("color"))
            .add_arg(&Atom::new_num(1))
            .finish();
        expression.0 = &expression.0 * color_fn;
        let replacements = vec![
            (
                Pattern::parse("f_(x___,aind(y___,cof(i__),z___))*color(a___)").unwrap(),
                Pattern::parse("color(a___*f_(x___,aind(y___,cof(i__),z___)))")
                    .unwrap()
                    .into(),
            ),
            (
                Pattern::parse("f_(x___,aind(y___,coaf(i__),z___))*color(a___)").unwrap(),
                Pattern::parse("color(a___*f_(x___,aind(y___,coaf(i__),z___)))")
                    .unwrap()
                    .into(),
            ),
            (
                Pattern::parse("f_(x___,aind(y___,coad(i__),z___))*color(a___)").unwrap(),
                Pattern::parse("color(a___*f_(x___,aind(y___,coad(i__),z___)))")
                    .unwrap()
                    .into(),
            ),
        ];
        let reps: Vec<Replacement> = replacements
            .iter()
            .map(|(lhs, rhs)| Replacement::new(lhs, rhs))
            .collect();

        expression.replace_repeat_multiple(&reps);
    }

    pub fn color_symplify_impl(mut expression: SerializableAtom) -> SerializableAtom {
        let (mut coefs, rem) = expression.0.coefficient_list(State::get_symbol("color"));

        let replacements = vec![
            (Pattern::parse("color(a___)").unwrap(),Pattern::parse("a___").unwrap().into()),
            (Pattern::parse("f_(x___,aind(y___,cof(i__),z___))*id(aind(coaf(i__),cof(j__)))").unwrap(),Pattern::parse("f_(x___,aind(y___,cof(j__),z___))").unwrap().into()),
            (Pattern::parse("f_(x___,aind(y___,cof(j__),z___))*id(aind(cof(i__),coaf(j__)))").unwrap(),Pattern::parse("f_(x___,aind(y___,cof(i__),z___))").unwrap().into()),
            (Pattern::parse("f_(x___,aind(y___,coaf(i__),z___))*id(aind(cof(j__),coaf(i__)))").unwrap(),Pattern::parse("f_(x___,aind(y___,coaf(j__),z___))").unwrap().into()),
            (Pattern::parse("f_(x___,aind(y___,coaf(i__),z___))*id(aind(coaf(j__),cof(i__)))").unwrap(),Pattern::parse("f_(x___,aind(y___,coaf(j__),z___))").unwrap().into()),
            (Pattern::parse("f_(x___,aind(y___,coad(i__),z___))*id(aind(coad(j__),coad(i__)))").unwrap(),Pattern::parse("f_(x___,aind(y___,coad(j__),z___))").unwrap().into()),
            (
                Pattern::parse("id(aind(coaf(3,a_),cof(3,a_)))").unwrap(),
                Pattern::parse("Nc").unwrap().into(),
            ),
            (
                Pattern::parse("id(aind(cof(3,a_),coaf(3,a_)))").unwrap(),
                Pattern::parse("Nc").unwrap().into(),
            ),
            (
                Pattern::parse("id(aind(coad(8,a_),coad(8,a_)))").unwrap(),
                Pattern::parse("Nc*Nc -1").unwrap().into(),
            ),
            (
                Pattern::parse("T(aind(coad(8,b_),cof(3,a_),coaf(3,a_)))").unwrap(),
                Pattern::parse("0").unwrap().into(),
            ),
            (
                Pattern::parse("T(aind(coad(8,c_),cof(3,a_),coaf(3,b_)))T(aind(coad(8,d_),cof(3,b_),coaf(3,a_)))").unwrap(),
                Pattern::parse("TR* id(aind(coad(8,c_),coad(8,d_)))").unwrap().into(),
            ),
            (
                Pattern::parse("T(aind(coad(8,g_),cof(3,a_),coaf(3,b_)))*T(aind(coad(8,g_),cof(3,c_),coaf(3,d_)))").unwrap(),
                Pattern::parse("TR* (id(aind(cof(3,a_),coaf(3,d_)))* id(aind(cof(3,c_),coaf(3,b_)))-1/Nc id(aind(cof(3,a_),coaf(3,b_)))* id(aind(cof(3,c_),coaf(3,d_))))").unwrap().into(),
            )
        ];

        let reps: Vec<Replacement> = replacements
            .iter()
            .map(|(lhs, rhs)| Replacement::new(lhs, rhs))
            .collect();

        let mut atom = Atom::new_num(0);
        for (key, coef) in coefs.iter_mut() {
            SerializableAtom::replace_repeat_multiple_atom_expand(key, &reps);

            atom = atom + coef.factor() * key.factor();
            // println!("coef {i}:{}\n", coef.factor());
        }
        atom = atom + rem;
        expression.0 = atom;
        expression
    }

    pub fn color_symplify<T: UnexpandedNumerator>(expr: T) -> ColorSymplified {
        let mut expr = expr.expr().unwrap().clone();
        Self::isolate_color(&mut expr);
        Self::new(Self::color_symplify_impl(expr))
    }

    pub fn parse(self) -> Network {
        let net = TensorNetwork::try_from(self.expression.0.as_view())
            .unwrap()
            .to_fully_parametric()
            .cast();

        // println!("net scalar{}", net.scalar.as_ref().unwrap());
        Network { net }
    }
}

impl Numerator<ColorSymplified> {
    pub fn gamma_symplify(self) -> Numerator<GammaSymplified> {
        debug!("ColorSymplified numerator: {}", self.export());
        debug!("gamma symplifying color symplified numerator");
        Numerator {
            state: GammaSymplified::gamma_symplify_impl(self.state.expression),
        }
    }

    pub fn color_project(self) -> Numerator<ColorProjected> {
        debug!("ColorSymplified numerator: {}", self.export());
        debug!("projecting color symplified numerator");
        Numerator {
            state: ColorProjected::color_project(self.state),
        }
    }

    pub fn parse(self) -> Numerator<Network> {
        debug!("ColorSymplified numerator: {}", self.export());
        Numerator {
            state: self.state.parse(),
        }
    }
}

impl ColorProjected {
    pub fn color_project(mut color_simple: ColorSymplified) -> ColorProjected {
        let reps = vec![
            (
                Pattern::parse("f_(a___,aind(b___,cof(c__),d___))").unwrap(),
                Pattern::parse("1").unwrap().into(),
            ),
            (
                Pattern::parse("f_(a___,aind(b___,cos(c__),d___))").unwrap(),
                Pattern::parse("1").unwrap().into(),
            ),
            (
                Pattern::parse("f_(a___,aind(b___,coaf(c__),d___))").unwrap(),
                Pattern::parse("1").unwrap().into(),
            ),
        ];

        let reps = reps
            .iter()
            .map(|(lhs, rhs)| Replacement::new(lhs, rhs))
            .collect_vec();
        color_simple.expression.replace_repeat_multiple(&reps);

        ColorProjected::new(color_simple.expression)
    }
}

impl Numerator<ColorProjected> {
    pub fn parse(self) -> Numerator<Network> {
        // debug!("ColorProjected numerator: {}", self.export());
        Numerator {
            state: Network::parse_impl(self.state.expression.0.as_view()),
        }
    }

    pub fn gamma_symplify(self) -> Numerator<GammaSymplified> {
        // debug!("ColorSymplified numerator: {}", self.export());
        debug!("gamma symplifying color symplified numerator");
        Numerator {
            state: GammaSymplified::gamma_symplify(self.state),
        }
    }
}

impl GammaSymplified {
    pub fn gamma_symplify_impl(mut expr: SerializableAtom) -> Self {
        expr.0 = expr.0.expand();
        let pats = [(
            Pattern::parse("id(aind(a_,b_))*t_(aind(d___,b_,c___))").unwrap(),
            Pattern::parse("t_(aind(d___,a_,c___))").unwrap().into(),
        )];

        let reps: Vec<Replacement> = pats
            .iter()
            .map(|(lhs, rhs)| Replacement::new(lhs, rhs))
            .collect();
        expr.replace_repeat_multiple(&reps);
        let pats = vec![
            (
                Pattern::parse("ProjP(aind(a_,b_))").unwrap(),
                Pattern::parse("1/2*id(aind(a_,b_))-1/2*gamma5(aind(a_,b_))")
                    .unwrap()
                    .into(),
            ),
            (
                Pattern::parse("ProjM(aind(a_,b_))").unwrap(),
                Pattern::parse("1/2*id(aind(a_,b_))+1/2*gamma5(aind(a_,b_))")
                    .unwrap()
                    .into(),
            ),
            (
                Pattern::parse("id(aind(a_,b_))*f_(c___,aind(d___,b_,e___))").unwrap(),
                Pattern::parse("f_(c___,aind(d___,a_,e___))")
                    .unwrap()
                    .into(),
            ),
            (
                Pattern::parse("id(aind(a_,b_))*f_(c___,aind(d___,a_,e___))").unwrap(),
                Pattern::parse("f_(c___,aind(d___,b_,e___))")
                    .unwrap()
                    .into(),
            ),
            (
                Pattern::parse("γ(aind(a_,b_,c_))*γ(aind(d_,c_,e_))").unwrap(),
                Pattern::parse("gamma_chain(aind(a_,d_,b_,e_))")
                    .unwrap()
                    .into(),
            ),
            (
                Pattern::parse("gamma_chain(aind(a__,b_,c_))*gamma_chain(aind(d__,c_,e_))")
                    .unwrap(),
                Pattern::parse("gamma_chain(aind(a__,d__,b_,e_))")
                    .unwrap()
                    .into(),
            ),
            (
                Pattern::parse("γ(aind(a_,b_,c_))*gamma_chain(aind(d__,c_,e_))").unwrap(),
                Pattern::parse("gamma_chain(aind(a_,d__,b_,e_))")
                    .unwrap()
                    .into(),
            ),
            (
                Pattern::parse("gamma_chain(aind(a__,b_,c_))*γ(aind(d_,c_,e_))").unwrap(),
                Pattern::parse("gamma_chain(aind(a__,d_,b_,e_))")
                    .unwrap()
                    .into(),
            ),
            (
                Pattern::parse("gamma_chain(aind(a__,b_,b_))").unwrap(),
                Pattern::parse("gamma_trace(aind(a__))").unwrap().into(),
            ),
        ];
        let reps: Vec<Replacement> = pats
            .iter()
            .map(|(lhs, rhs)| Replacement::new(lhs, rhs))
            .collect();
        expr.0 = expr.0.expand();
        expr.replace_repeat_multiple(&reps);
        expr.0 = expr.0.expand();
        expr.replace_repeat_multiple(&reps);

        let pat = Pattern::parse("gamma_trace(a__)").unwrap();

        let set = MatchSettings::default();
        let cond = Condition::default();

        let mut it = pat.pattern_match(expr.0.as_view(), &cond, &set);

        let mut max_nargs = 0;
        while let Some(a) = it.next() {
            for (_, v) in a.match_stack {
                match v {
                    Match::Single(s) => {
                        match s {
                            AtomView::Fun(f) => {
                                let a = f.get_nargs();
                                if a > max_nargs {
                                    max_nargs = a;
                                }
                            }
                            _ => {
                                panic!("should be a function")
                            }
                        }
                        // print!("{}", s)
                    }
                    _ => panic!("should be a single match"),
                }
                // println!();
            }
        }

        let mut reps = vec![];
        for n in 1..=max_nargs {
            if n % 2 == 0 {
                let mut sum = Atom::new_num(0);

                // sum((-1)**(k+1) * d(p_[0], p_[k]) * f(*p_[1:k], *p_[k+1:l])
                for j in 1..n {
                    let gamma_chain_builder =
                        FunctionBuilder::new(State::get_symbol("gamma_trace"));

                    let mut gamma_chain_builder_slots =
                        FunctionBuilder::new(State::get_symbol("aind"));

                    let metric_builder = FunctionBuilder::new(State::get_symbol("g"));

                    let metric_builder_slots = FunctionBuilder::new(State::get_symbol("aind"));

                    for k in 1..j {
                        let mu = Atom::parse(&format!("a{}_", k)).unwrap();
                        gamma_chain_builder_slots = gamma_chain_builder_slots.add_arg(&mu);
                    }

                    for k in (j + 1)..n {
                        let mu = Atom::parse(&format!("a{}_", k)).unwrap();
                        gamma_chain_builder_slots = gamma_chain_builder_slots.add_arg(&mu);
                    }

                    let metric = metric_builder
                        .add_arg(
                            &metric_builder_slots
                                .add_args(&[
                                    &Atom::parse(&format!("a{}_", 0)).unwrap(),
                                    &Atom::parse(&format!("a{}_", j)).unwrap(),
                                ])
                                .finish(),
                        )
                        .finish();

                    let gamma = &gamma_chain_builder
                        .add_arg(&gamma_chain_builder_slots.finish())
                        .finish()
                        * &metric;

                    if j % 2 == 0 {
                        sum = &sum - &gamma;
                    } else {
                        sum = &sum + &gamma;
                    }
                }

                let gamma_chain_builder = FunctionBuilder::new(State::get_symbol("gamma_trace"));
                let mut gamma_chain_builder_slots = FunctionBuilder::new(State::get_symbol("aind"));
                for k in 0..n {
                    let mu = Atom::parse(&format!("a{}_", k)).unwrap();
                    gamma_chain_builder_slots = gamma_chain_builder_slots.add_arg(&mu);
                }
                let a = gamma_chain_builder
                    .add_arg(&gamma_chain_builder_slots.finish())
                    .finish();

                reps.push((a.into_pattern(), sum.into_pattern().into()));
            } else {
                let gamma_chain_builder = FunctionBuilder::new(State::get_symbol("gamma_trace"));
                let mut gamma_chain_builder_slots = FunctionBuilder::new(State::get_symbol("aind"));
                for k in 0..n {
                    let mu = Atom::parse(&format!("a{}_", k)).unwrap();
                    gamma_chain_builder_slots = gamma_chain_builder_slots.add_arg(&mu);
                }
                let a = gamma_chain_builder
                    .add_arg(&gamma_chain_builder_slots.finish())
                    .finish();
                // println!("{}", a);
                reps.push((a.into_pattern(), Atom::new_num(0).into_pattern().into()));
            }
        }

        reps.push((
            Pattern::parse("gamma_trace(aind())").unwrap(),
            Pattern::parse("4").unwrap().into(),
        ));

        // Dd
        reps.push((
            Pattern::parse("f_(i_,aind(loru(a__)))*g(aind(lord(a__),lord(b__)))").unwrap(),
            Pattern::parse("f_(i_,aind(lord(b__)))").unwrap().into(),
        ));
        // Du
        reps.push((
            Pattern::parse("f_(i_,aind(loru(a__)))*g(aind(lord(a__),loru(b__)))").unwrap(),
            Pattern::parse("f_(i_,aind(loru(b__)))").unwrap().into(),
        ));
        // Uu
        reps.push((
            Pattern::parse("f_(i_,aind(lord(a__)))*g(aind(loru(a__),loru(b__)))").unwrap(),
            Pattern::parse("f_(i_,aind(loru(b__)))").unwrap().into(),
        ));
        // Ud
        reps.push((
            Pattern::parse("f_(i_,aind(lord(a__)))*g(aind(loru(a__),lord(b__)))").unwrap(),
            Pattern::parse("f_(i_,aind(lord(b__)))").unwrap().into(),
        ));

        // dD
        reps.push((
            Pattern::parse("f_(i_,aind(loru(a__)))*g(aind(lord(b__),lord(a__)))").unwrap(),
            Pattern::parse("f_(i_,aind(lord(b__)))").unwrap().into(),
        ));
        // uD
        reps.push((
            Pattern::parse("f_(i_,aind(loru(a__)))*g(aind(loru(b__),lord(a__)))").unwrap(),
            Pattern::parse("f_(i_,aind(loru(b__)))").unwrap().into(),
        ));
        // uU
        reps.push((
            Pattern::parse("f_(i_,aind(lord(a__)))*g(aind(loru(b__),loru(a__)))").unwrap(),
            Pattern::parse("f_(i_,aind(loru(b__)))").unwrap().into(),
        ));
        // dU
        reps.push((
            Pattern::parse("f_(i_,aind(lord(a__)))*g(aind(lord(b__),loru(a__)))").unwrap(),
            Pattern::parse("f_(i_,aind(lord(b__)))").unwrap().into(),
        ));

        let reps = reps
            .iter()
            .map(|(lhs, rhs)| Replacement::new(lhs, rhs))
            .collect_vec();
        expr.replace_repeat_multiple(&reps);
        expr.0 = expr.0.expand();
        expr.replace_repeat_multiple(&reps);

        GammaSymplified::new(expr)
    }

    pub fn gamma_symplify(color: ColorProjected) -> GammaSymplified {
        Self::gamma_symplify_impl(color.expression)
    }

    pub fn parse(self) -> Network {
        let net = TensorNetwork::try_from(self.expression.0.as_view())
            .unwrap()
            .to_fully_parametric()
            .cast();

        // println!("net scalar{}", net.scalar.as_ref().unwrap());
        Network { net }
    }
}

impl Numerator<GammaSymplified> {
    pub fn parse(self) -> Numerator<Network> {
        // debug!("GammaSymplified numerator: {}", self.export());
        debug!("parsing numerator into tensor network");
        Numerator {
            state: self.state.parse(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
pub struct Network {
    #[bincode(with_serde)]
    net: TensorNetwork<ParamTensor<AtomStructure>, SerializableAtom>,
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
pub enum ContractionSettings<'a, 'b, R> {
    Levelled((usize, &'a mut FunctionMap<'b, R>)),
    Normal,
}

impl Network {
    pub fn parse_impl(expr: AtomView) -> Self {
        let net = TensorNetwork::try_from(expr)
            .unwrap()
            .to_fully_parametric()
            .cast();

        // println!("net scalar{}", net.scalar.as_ref().unwrap());
        Network { net }
    }

    pub fn contract<R>(mut self, settings: ContractionSettings<R>) -> Result<Contracted> {
        match settings {
            ContractionSettings::Levelled((_depth, _fn_map)) => {
                let _levels: Levels<_, _> = self.net.into();
                // levels.contract(depth, fn_map);
                unimplemented!("cannot because of attached dummy lifetime...")
            }
            ContractionSettings::Normal => {
                self.net.contract();
                let tensor = self.net.result_tensor_smart()?;
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

impl Numerator<Network> {
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
    pub tensor: ParamTensor<AtomStructure>,
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
        export_settings: &ExportSettings,
        params: &[Atom],
    ) -> EvaluatorSingle {
        debug!("generating single evaluator");
        let mut fn_map: FunctionMap = FunctionMap::new();

        Numerator::<Contracted>::add_consts_to_fn_map(&mut fn_map);

        debug!("Generate eval tree set with {} params", params.len());

        let mut eval_tree = self.tensor.eval_tree(&fn_map, params).unwrap();
        debug!("Horner scheme");

        eval_tree.horner_scheme();
        debug!("Common subexpression elimination");
        eval_tree.common_subexpression_elimination();
        debug!("Linearize double");
        let cpe_rounds = export_settings
            .numerator_settings
            .eval_settings
            .cpe_rounds();
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
        let compiled = if export_settings
            .numerator_settings
            .eval_settings
            .compile_options()
            .compile()
        {
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
            let inline_asm = export_settings.gammaloop_compile_options.inline_asm();

            let compile_options = export_settings
                .gammaloop_compile_options
                .to_symbolica_compile_options();

            CompiledEvaluator::new(
                eval.export_cpp(&filename, &function_name, true, inline_asm)
                    .unwrap()
                    .compile(&library_name, compile_options)
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
                let e = FunctionBuilder::new(State::get_symbol(&name));
                data.push(
                    e.add_arg(&Atom::new_num(num))
                        .add_arg(&Atom::parse(&format!("cind({})", index)).unwrap())
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
        params.push(Atom::new_var(State::I));

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

        for (i, ext) in graph.edges.iter().enumerate() {
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
            match ext.edge_type {
                EdgeType::Incoming => pols.extend(
                    ext.particle
                        .incoming_polarization_atom_concrete(&ext.in_slot(graph), i),
                ),

                EdgeType::Outgoing => pols.extend(
                    ext.particle
                        .outgoing_polarization_atom_concrete(&ext.in_slot(graph), i),
                ),
                _ => {}
            }
        }

        param_values.extend(vec![Complex::new_zero(); pols.len()]);
        params.extend(pols);

        param_values.push(Complex::new_i());
        params.push(Atom::new_var(State::I));

        let model_params_start = params.len();
        param_values.extend(model.generate_values());
        params.extend(model.generate_params());

        let quad_values = param_values.iter().map(|c| c.map(|f| f.higher())).collect();

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
        let single = self
            .state
            .evaluator(extra_info.path.clone(), name, export_settings, params);

        match export_settings.numerator_settings.eval_settings {
            NumeratorEvaluatorOptions::Joint(_) => Numerator {
                state: Evaluators {
                    orientated: Some(single.orientated_joint_impl(
                        n_edges,
                        name,
                        params,
                        extra_info,
                        export_settings,
                    )),
                    single,
                    choice: SingleOrCombined::Combined,
                    orientations: extra_info.orientations.clone(),
                    quad_param_values,
                    double_param_values,
                    model_params_start,
                    emr_len: n_edges,
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
                        export_settings,
                        iterations,
                        n_cores,
                        verbose,
                    )),
                    single,
                    choice: SingleOrCombined::Combined,
                    orientations: extra_info.orientations.clone(),
                    quad_param_values,
                    double_param_values,
                    model_params_start,
                    emr_len: n_edges,
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
        debug!("generating evaluators for contracted numerator");
        let (params, double_param_values, quad_param_values, model_params_start) =
            Contracted::generate_params(graph, model);

        trace!("params length:{}", params.len());
        for (i, p) in params.iter().enumerate() {
            println!("\t{i} {p}");
        }

        let emr_len = graph.edges.len();

        let single = self.state.evaluator(
            extra_info.path.clone(),
            &graph.name,
            export_settings,
            &params,
        );

        match export_settings.numerator_settings.eval_settings {
            NumeratorEvaluatorOptions::Joint(_) => Numerator {
                state: Evaluators {
                    orientated: Some(single.orientated_joint(
                        graph,
                        &params,
                        extra_info,
                        export_settings,
                    )),
                    single,
                    choice: SingleOrCombined::Combined,
                    orientations: extra_info.orientations.clone(),
                    quad_param_values,
                    double_param_values,
                    model_params_start,
                    emr_len,
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
                        export_settings,
                        iterations,
                        n_cores,
                        verbose,
                    )),
                    single,
                    choice: SingleOrCombined::Combined,
                    orientations: extra_info.orientations.clone(),
                    quad_param_values,
                    double_param_values,
                    model_params_start,
                    emr_len,
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
#[serde(tag = "type")]
pub enum NumeratorEvaluatorOptions {
    #[serde(rename = "Single")]
    Single(EvaluatorOptions),
    #[serde(rename = "Joint")]
    Joint(EvaluatorOptions),
    #[serde(rename = "Iterative")]
    Iterative(IterativeOptions),
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
    pub eval_double: LinearizedEvalTensorSet<Complex<F<f64>>, AtomStructure>,
    pub eval_quad: LinearizedEvalTensorSet<Complex<F<f128>>, AtomStructure>,
    pub compiled: CompiledEvaluator<EvalTensorSet<SerializableCompiledEvaluator, AtomStructure>>,
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
    tensor: ParamTensor<AtomStructure>,
    eval_double: EvalTensor<ExpressionEvaluator<Complex<F<f64>>>, AtomStructure>,
    eval_quad: EvalTensor<ExpressionEvaluator<Complex<F<f128>>>, AtomStructure>,
    compiled: CompiledEvaluator<EvalTensor<SerializableCompiledEvaluator, AtomStructure>>,
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
        export_settings: &ExportSettings,
        iterations: usize,
        n_cores: usize,
        verbose: bool,
    ) -> EvaluatorOrientations {
        let mut fn_map: FunctionMap = FunctionMap::new();

        Numerator::<Contracted>::add_consts_to_fn_map(&mut fn_map);

        let mut seen = 0;

        let mut index_map = IndexSet::new();
        let mut positions = Vec::new();

        let len = extra_info.orientations.len();

        let reps = (0..n_edges)
            .map(|i| {
                (
                    Pattern::parse(&format!("Q({},cind(0))", i)).unwrap(),
                    Pattern::parse(&format!("-Q({},cind(0))", i))
                        .unwrap()
                        .into(),
                )
            })
            .collect_vec();

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
                        Some(Replacement::new(lhs, rhs))
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

        let cpe_rounds = export_settings
            .numerator_settings
            .eval_settings
            .cpe_rounds();

        debug!("Generate eval tree set with {} params", params.len());

        let mut eval_tree = set.eval_tree(&fn_map, params).unwrap();

        debug!("{} tensors in eval_tree", eval_tree.len());
        debug!("Horner scheme");

        eval_tree.horner_scheme();
        debug!("Common subexpression elimination");
        eval_tree.common_subexpression_elimination();
        debug!("Linearize double");

        let mut linearized = eval_tree.linearize(cpe_rounds);

        debug!("{} linearized tensors", linearized.len());

        for (i, t) in new_index_map.iter().enumerate() {
            let eval_tree = t.eval_tree(&fn_map, params).unwrap();
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

        let compiled = if export_settings
            .numerator_settings
            .eval_settings
            .compile_options()
            .compile()
        {
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
            let inline_asm = export_settings.gammaloop_compile_options.inline_asm();

            let compile_options = export_settings
                .gammaloop_compile_options
                .to_symbolica_compile_options();
            CompiledEvaluator::new(
                eval.export_cpp(&filename, &function_name, true, inline_asm)
                    .unwrap()
                    .compile(&library_name, compile_options)
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
        export_settings: &ExportSettings,
        iterations: usize,
        n_cores: usize,
        verbose: bool,
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
        )
    }

    pub fn orientated_joint_impl(
        &self,
        n_edges: usize,
        name: &str,
        params: &[Atom],
        extra_info: &ExtraInfo,
        export_settings: &ExportSettings,
    ) -> EvaluatorOrientations {
        let mut fn_map: FunctionMap = FunctionMap::new();

        Numerator::<Contracted>::add_consts_to_fn_map(&mut fn_map);

        let mut seen = 0;

        let mut index_map = IndexSet::new();
        let mut positions = Vec::new();

        let len = extra_info.orientations.len();

        let reps = (0..n_edges)
            .map(|i| {
                (
                    Pattern::parse(&format!("Q({},cind(0))", i)).unwrap(),
                    Pattern::parse(&format!("-Q({},cind(0))", i))
                        .unwrap()
                        .into(),
                )
            })
            .collect_vec();

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
                        Some(Replacement::new(lhs, rhs))
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

        let cpe_rounds = export_settings
            .numerator_settings
            .eval_settings
            .cpe_rounds();

        debug!("Generate eval tree set with {} params", params.len());

        let mut eval_tree = set.eval_tree(&fn_map, params).unwrap();
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

        let compiled = if export_settings
            .numerator_settings
            .eval_settings
            .compile_options()
            .compile()
        {
            debug!("compiling joint evaluator");
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
            let inline_asm = export_settings.gammaloop_compile_options.inline_asm();

            let compile_options = export_settings
                .gammaloop_compile_options
                .to_symbolica_compile_options();
            CompiledEvaluator::new(
                eval.export_cpp(&filename, &function_name, true, inline_asm)
                    .unwrap()
                    .compile(&library_name, compile_options)
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
        export_settings: &ExportSettings,
    ) -> EvaluatorOrientations {
        debug!("generate joint evaluator");
        self.orientated_joint_impl(
            graph.edges.len(),
            &graph.name,
            params,
            extra_info,
            export_settings,
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

#[derive(Clone, Serialize, Deserialize, Debug, Encode, Decode)]
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
        generate: Option<(&Model, &BareGraph, &ExtraInfo, &ExportSettings)>,
    ) {
        if self.state.orientated.is_some() {
            self.state.choice = SingleOrCombined::Combined;
        } else if let Some((model, graph, extra_info, export_settings)) = generate {
            let (params, _, _, _) = Contracted::generate_params(graph, model);
            let orientated =
                self.state
                    .single
                    .orientated_joint(graph, &params, extra_info, export_settings);
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
    #[error("Not ColorSymplified")]
    NotColorSymplified,
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

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
#[allow(clippy::large_enum_variant)]
pub enum PythonState {
    UnInit(Option<UnInit>),
    Global(Option<Global>),
    AppliedFeynmanRule(Option<AppliedFeynmanRule>),
    ColorSymplified(Option<ColorSymplified>),
    ColorProjected(Option<ColorProjected>),
    GammaSymplified(Option<GammaSymplified>),
    Network(Option<Network>),
    Contracted(Option<Contracted>),
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
            PythonState::ColorSymplified(state) => {
                if let Some(s) = state {
                    s.export()
                } else {
                    "None".into()
                }
            }

            PythonState::ColorProjected(state) => {
                if let Some(s) = state {
                    s.export()
                } else {
                    "None".into()
                }
            }
            PythonState::GammaSymplified(state) => {
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
            PythonState::ColorSymplified(state) => {
                if let Some(s) = state {
                    s.update_model(model)
                } else {
                    Err(NumeratorStateError::NoneVariant.into())
                }
            }
            PythonState::GammaSymplified(state) => {
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

impl UnexpandedNumerator for PythonState {
    fn expr(&self) -> Result<&SerializableAtom, NumeratorStateError> {
        match self {
            PythonState::Global(state) => {
                if let Some(s) = state {
                    s.expr()
                } else {
                    Err(NumeratorStateError::NoneVariant)
                }
            }
            PythonState::AppliedFeynmanRule(state) => {
                if let Some(s) = state {
                    s.expr()
                } else {
                    Err(NumeratorStateError::NoneVariant)
                }
            }
            PythonState::ColorSymplified(state) => {
                if let Some(s) = state {
                    s.expr()
                } else {
                    Err(NumeratorStateError::NoneVariant)
                }
            }
            PythonState::GammaSymplified(state) => {
                if let Some(s) = state {
                    s.expr()
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
