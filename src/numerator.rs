use std::fmt::Debug;
use std::path::PathBuf;
use std::time::Instant;

use crate::graph::BareGraph;
use crate::momentum::Polarization;
use crate::utils::f128;
use crate::ExportSettings;
use crate::{
    graph::EdgeType,
    model::Model,
    momentum::FourMomentum,
    utils::{FloatLike, F},
};
use bincode::{Decode, Encode};
use color_eyre::{Report, Result};
use eyre::eyre;
use gat_lending_iterator::LendingIterator;
// use gxhash::GxBuildHasher;
use indexmap::IndexSet;
use itertools::Itertools;

use log::{debug, info};
use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};
use spenso::data::DataTensor;

use spenso::network::Levels;
use spenso::parametric::{
    EvalTensor, EvalTensorSet, ParamTensorSet, SerializableAtom, SerializableCompiledEvaluator,
    TensorSet,
};
use spenso::structure::{HasStructure, SerializableSymbol, SmartShadowStructure};
use spenso::{
    complex::Complex,
    network::TensorNetwork,
    parametric::{ParamTensor, PatternReplacement},
    structure::{Lorentz, NamedStructure, PhysReps, RepName, Shadowable, TensorStructure},
};
use symbolica::evaluate::{CompileOptions, ExpressionEvaluator, InlineASM};
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
    ) -> Result<RepeatingIteratorTensorOrScalar<DataTensor<Complex<F<T>>, AtomStructure>>>;

    fn evaluate_single(
        &mut self,
        emr: &[FourMomentum<F<T>>],
        polarizations: &[Polarization<Complex<F<T>>>],
    ) -> DataTensor<Complex<F<T>>, AtomStructure>;
}

impl<T: FloatLike> Evaluate<T> for Num<Evaluators> {
    fn evaluate_all_orientations(
        &mut self,
        emr: &[FourMomentum<F<T>>],
        polarizations: &[Polarization<Complex<F<T>>>],
    ) -> Result<RepeatingIteratorTensorOrScalar<DataTensor<Complex<F<T>>, AtomStructure>>> {
        <T as NumeratorEvaluateFloat>::update_params(self, emr, polarizations);
        <T as NumeratorEvaluateFloat>::evaluate_all_orientations(self)
    }

    fn evaluate_single(
        &mut self,
        emr: &[FourMomentum<F<T>>],
        polarizations: &[Polarization<Complex<F<T>>>],
    ) -> DataTensor<Complex<F<T>>, AtomStructure> {
        <T as NumeratorEvaluateFloat>::update_params(self, emr, polarizations);
        <T as NumeratorEvaluateFloat>::evaluate_single(self)
    }
}

pub trait NumeratorEvaluateFloat<T: FloatLike = Self> {
    fn evaluate_all_orientations(
        num: &mut Num<Evaluators>,
    ) -> Result<RepeatingIteratorTensorOrScalar<DataTensor<Complex<F<T>>, AtomStructure>>>;

    fn evaluate_single(num: &mut Num<Evaluators>) -> DataTensor<Complex<F<T>>, AtomStructure>;

    fn update_params(
        num: &mut Num<Evaluators>,
        emr: &[FourMomentum<F<T>>],
        polarizations: &[Polarization<Complex<F<T>>>],
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

impl<T> LendingIterator for RepeatingIterator<T> {
    type Item<'a> = &'a T where Self:'a ;

    fn next(&mut self) -> Option<Self::Item<'_>> {
        Some(&self.elements[self.positions.next()?])
    }
}

impl NumeratorEvaluateFloat for f64 {
    fn update_params(
        num: &mut Num<Evaluators>,
        emr: &[FourMomentum<F<Self>>],
        polarizations: &[Polarization<Complex<F<Self>>>],
    ) {
        let params = &mut num.state.double_param_values;
        emr.iter()
            .flat_map(|p| p.into_iter().cloned().map(Complex::new_re))
            .enumerate()
            .for_each(|(i, c)| params[i] = c);

        let mut i = 4 * emr.len();

        for p in polarizations {
            if !p.is_scalar() {
                for &pi in p {
                    params[i] = pi;
                    i += 1;
                }
            }
        }
    }

    fn evaluate_all_orientations(
        num: &mut Num<Evaluators>,
    ) -> Result<RepeatingIteratorTensorOrScalar<DataTensor<Complex<F<Self>>, AtomStructure>>> {
        let params = &num.state.double_param_values;
        if num.state.single.param_len != params.len() {
            return Err(eyre!(
                "params length mismatch {} not equal to {}",
                num.state.single.param_len,
                params.len()
            ));
        }
        if let Some(orientation_evaluator) = &mut num.state.orientated {
            let pos = orientation_evaluator.positions.clone();
            let res = if let Some(c) = &mut orientation_evaluator.compiled.evaluator {
                c.evaluate(params)
            } else {
                orientation_evaluator.eval_double.evaluate(params)
            };

            Ok((res, pos).into())
        } else {
            let mut oriented_params = params.to_vec();
            let base_params = params.clone();
            let mut tensors = Vec::new();
            let orientations = num.state.orientations.clone();
            for o in orientations {
                for (i, &sign) in o.iter().enumerate() {
                    if sign {
                        oriented_params[i] = -oriented_params[i];
                    }
                }
                tensors.push(<Self as NumeratorEvaluateFloat>::evaluate_single(num));
                oriented_params = base_params.clone();
            }
            Ok(RepeatingIteratorTensorOrScalar::Tensors(
                RepeatingIterator::new_not_repeating(tensors),
            ))
        }
    }

    fn evaluate_single(num: &mut Num<Evaluators>) -> DataTensor<Complex<F<Self>>, AtomStructure> {
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
        num: &mut Num<Evaluators>,
        emr: &[FourMomentum<F<Self>>],
        polarizations: &[Polarization<Complex<F<Self>>>],
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
    }

    fn evaluate_all_orientations(
        num: &mut Num<Evaluators>,
    ) -> Result<RepeatingIteratorTensorOrScalar<DataTensor<Complex<F<Self>>, AtomStructure>>> {
        let params = &num.state.quad_param_values;
        if num.state.single.param_len != params.len() {
            return Err(eyre!("params length mismatch"));
        }
        if let Some(orientation_evaluator) = &mut num.state.orientated {
            let pos = orientation_evaluator.positions.clone();
            let res = orientation_evaluator.eval_quad.evaluate(params);
            Ok((res, pos).into())
        } else {
            let mut oriented_params = params.to_vec();
            let base_params = params.clone();
            let mut tensors = Vec::new();
            let orientations = num.state.orientations.clone();
            for o in orientations {
                for (i, &sign) in o.iter().enumerate() {
                    if sign {
                        oriented_params[i] *= F(f128::from_f64(-1.0));
                    }
                }
                tensors.push(<Self as NumeratorEvaluateFloat>::evaluate_single(num));
                oriented_params = base_params.clone();
            }
            Ok(RepeatingIteratorTensorOrScalar::Tensors(
                RepeatingIterator::new_not_repeating(tensors),
            ))
        }
    }

    fn evaluate_single(num: &mut Num<Evaluators>) -> DataTensor<Complex<F<Self>>, AtomStructure> {
        let params = &num.state.quad_param_values;
        if num.state.single.param_len != params.len() {
            panic!("params length mismatch");
        }
        num.state.single.eval_quad.evaluate(params)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
pub struct Num<State> {
    pub state: State,
}

impl<S: NumeratorState> Num<S> {
    pub fn export(&self) -> String {
        self.state.export()
    }

    pub fn forget_type(self) -> Num<PythonState> {
        Num {
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

impl<S: UnexpandedNumerator> Num<S> {
    pub fn expr(&self) -> Result<SerializableAtom, NumeratorStateError> {
        self.state.expr()
    }
}

pub trait TypedNumeratorState:
    NumeratorState + TryFrom<PythonState, Error: std::error::Error + Send + Sync + 'static>
{
    fn apply<F, S: TypedNumeratorState>(
        state: &mut Num<PythonState>,
        f: F,
    ) -> Result<(), NumeratorStateError>
    where
        F: FnMut(Num<Self>) -> Num<S>;
}

impl Num<PythonState> {
    pub fn try_from<S: TypedNumeratorState>(self) -> Result<Num<S>, Report> {
        Ok(Num {
            state: self.state.try_into()?,
        })
    }

    pub fn apply<F, S: TypedNumeratorState, T: TypedNumeratorState>(
        &mut self,
        f: F,
    ) -> Result<(), NumeratorStateError>
    where
        F: FnMut(Num<S>) -> Num<T>,
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

pub trait UnexpandedNumerator: NumeratorState {
    fn expr(&self) -> Result<SerializableAtom, NumeratorStateError>;
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

    fn update_model(&mut self, model: &Model) -> Result<()> {
        Err(eyre!("Uninitialized, nothing to update"))
    }
}

impl TypedNumeratorState for UnInit {
    fn apply<F, S: TypedNumeratorState>(
        num: &mut Num<PythonState>,
        mut f: F,
    ) -> Result<(), NumeratorStateError>
    where
        F: FnMut(Num<Self>) -> Num<S>,
    {
        if let PythonState::UnInit(s) = &mut num.state {
            if let Some(s) = s.take() {
                *num = f(Num { state: s }).forget_type();
                return Ok(());
            } else {
                return Err(NumeratorStateError::NoneVariant);
            }
        }
        Err(NumeratorStateError::NotUnit)
    }
}

impl Num<UnInit> {
    pub fn from_graph(self, graph: &mut BareGraph) -> Num<AppliedFeynmanRule> {
        Num {
            state: AppliedFeynmanRule::from_graph(graph),
        }
    }
}

impl Default for Num<UnInit> {
    fn default() -> Self {
        Num { state: UnInit }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
pub struct AppliedFeynmanRule {
    #[bincode(with_serde)]
    expression: SerializableAtom,
}

impl TryFrom<PythonState> for AppliedFeynmanRule {
    type Error = NumeratorStateError;

    fn try_from(value: PythonState) -> std::result::Result<Self, Self::Error> {
        match value {
            PythonState::AppliedFeynmanRule(s) => {
                if let Some(s) = s {
                    Ok(s)
                } else {
                    Err(NumeratorStateError::NoneVariant)
                }
            }
            _ => Err(NumeratorStateError::NotAppliedFeynmanRule),
        }
    }
}

impl AppliedFeynmanRule {
    pub fn from_graph(graph: &mut BareGraph) -> Self {
        debug!("generating numerator for graph: {}", graph.name);

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
            &(&a + &b * &i).into_pattern(),
            None,
            None,
        );

        AppliedFeynmanRule {
            expression: builder.into(),
        }
    }

    fn isolate_color(&mut self) {
        let color_fn = FunctionBuilder::new(State::get_symbol("color"))
            .add_arg(&Atom::new_num(1))
            .finish();
        self.expression.0 = &self.expression.0 * color_fn;
        let replacements = vec![
            (
                Pattern::parse("f_(x___,aind(y___,cof(i__),z___))*color(a___)").unwrap(),
                Pattern::parse("color(a___*f_(x___,aind(y___,cof(i__),z___)))").unwrap(),
            ),
            (
                Pattern::parse("f_(x___,aind(y___,coaf(i__),z___))*color(a___)").unwrap(),
                Pattern::parse("color(a___*f_(x___,aind(y___,coaf(i__),z___)))").unwrap(),
            ),
            (
                Pattern::parse("f_(x___,aind(y___,coad(i__),z___))*color(a___)").unwrap(),
                Pattern::parse("color(a___*f_(x___,aind(y___,coad(i__),z___)))").unwrap(),
            ),
        ];
        let reps: Vec<Replacement> = replacements
            .iter()
            .map(|(lhs, rhs)| Replacement::new(lhs, rhs))
            .collect();

        self.expression.replace_repeat_multiple(&reps);
    }

    pub fn color_symplify(mut self) -> ColorSymplified {
        self.isolate_color();
        let (mut coefs, rem) = self
            .expression
            .0
            .coefficient_list(State::get_symbol("color"));

        let replacements = vec![
            (Pattern::parse("color(a___)").unwrap(),Pattern::parse("a___").unwrap()),
            (Pattern::parse("f_(x___,aind(y___,cof(i__),z___))*id(aind(coaf(i__),cof(j__)))").unwrap(),Pattern::parse("f_(x___,aind(y___,cof(j__),z___))").unwrap()),
            (Pattern::parse("f_(x___,aind(y___,cof(j__),z___))*id(aind(cof(i__),coaf(j__)))").unwrap(),Pattern::parse("f_(x___,aind(y___,cof(i__),z___))").unwrap()),
            (Pattern::parse("f_(x___,aind(y___,coaf(i__),z___))*id(aind(cof(j__),coaf(i__)))").unwrap(),Pattern::parse("f_(x___,aind(y___,coaf(j__),z___))").unwrap()),
            (Pattern::parse("f_(x___,aind(y___,coaf(i__),z___))*id(aind(coaf(j__),cof(i__)))").unwrap(),Pattern::parse("f_(x___,aind(y___,coaf(j__),z___))").unwrap()),
            (Pattern::parse("f_(x___,aind(y___,coad(i__),z___))*id(aind(coad(j__),coad(i__)))").unwrap(),Pattern::parse("f_(x___,aind(y___,coad(j__),z___))").unwrap()),
            (
                Pattern::parse("id(aind(coaf(3,a_),cof(3,a_)))").unwrap(),
                Pattern::parse("Nc").unwrap(),
            ),
            (
                Pattern::parse("id(aind(cof(3,a_),coaf(3,a_)))").unwrap(),
                Pattern::parse("Nc").unwrap(),
            ),
            (
                Pattern::parse("id(aind(coad(8,a_),coad(8,a_)))").unwrap(),
                Pattern::parse("Nc*Nc -1").unwrap(),
            ),
            (
                Pattern::parse("T(aind(coad(8,b_),cof(3,a_),coaf(3,a_)))").unwrap(),
                Pattern::parse("0").unwrap(),
            ),
            (
                Pattern::parse("T(aind(coad(8,c_),cof(3,a_),coaf(3,b_)))T(aind(coad(8,d_),cof(3,b_),coaf(3,a_)))").unwrap(),
                Pattern::parse("TR* id(aind(coad(8,c_),coad(8,d_)))").unwrap(),
            ),
            (
                Pattern::parse("T(aind(coad(8,g_),cof(3,a_),coaf(3,b_)))*T(aind(coad(8,g_),cof(3,c_),coaf(3,d_)))").unwrap(),
                Pattern::parse("TR* (id(aind(cof(3,a_),coaf(3,d_)))* id(aind(cof(3,c_),coaf(3,b_)))-1/Nc id(aind(cof(3,a_),coaf(3,b_)))* id(aind(cof(3,c_),coaf(3,d_))))").unwrap(),
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
        self.expression.0 = atom;
        ColorSymplified {
            expression: self.expression,
        }
    }
}

impl UnexpandedNumerator for AppliedFeynmanRule {
    fn expr(&self) -> Result<SerializableAtom, NumeratorStateError> {
        Ok(self.expression.clone())
    }
}

impl NumeratorState for AppliedFeynmanRule {
    fn export(&self) -> String {
        self.expression.to_string()
    }

    fn forget_type(self) -> PythonState {
        PythonState::AppliedFeynmanRule(Some(self))
    }

    fn update_model(&mut self, model: &Model) -> Result<()> {
        Err(eyre!("Only applied feynman rule, nothing to update"))
    }
}

impl TypedNumeratorState for AppliedFeynmanRule {
    fn apply<F, S: TypedNumeratorState>(
        num: &mut Num<PythonState>,
        mut f: F,
    ) -> Result<(), NumeratorStateError>
    where
        F: FnMut(Num<Self>) -> Num<S>,
    {
        if let PythonState::AppliedFeynmanRule(s) = &mut num.state {
            if let Some(s) = s.take() {
                *num = f(Num { state: s }).forget_type();
                return Ok(());
            } else {
                return Err(NumeratorStateError::NoneVariant);
            }
        }
        Err(NumeratorStateError::NotAppliedFeynmanRule)
    }
}

impl Num<AppliedFeynmanRule> {
    pub fn color_symplify(self) -> Num<ColorSymplified> {
        Num {
            state: self.state.color_symplify(),
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
pub struct ColorSymplified {
    #[bincode(with_serde)]
    pub expression: SerializableAtom,
}

impl TryFrom<PythonState> for ColorSymplified {
    type Error = NumeratorStateError;

    fn try_from(value: PythonState) -> std::result::Result<Self, Self::Error> {
        match value {
            PythonState::ColorSymplified(s) => {
                if let Some(s) = s {
                    Ok(s)
                } else {
                    Err(NumeratorStateError::NoneVariant)
                }
            }
            _ => Err(NumeratorStateError::NotColorSymplified),
        }
    }
}

impl ColorSymplified {
    pub fn gamma_symplify(self) -> GammaSymplified {
        let replacements = vec![(
            Pattern::parse("color(a___)").unwrap(),
            Pattern::parse("a___").unwrap(),
        )];

        let reps: Vec<Replacement> = replacements
            .iter()
            .map(|(lhs, rhs)| Replacement::new(lhs, rhs))
            .collect();

        GammaSymplified {
            expression: self.expression,
        }
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

impl NumeratorState for ColorSymplified {
    fn export(&self) -> String {
        self.expression.0.to_string()
    }

    fn forget_type(self) -> PythonState {
        PythonState::ColorSymplified(Some(self))
    }

    fn update_model(&mut self, model: &Model) -> Result<()> {
        Err(eyre!(
            "Only applied feynman rule,and simplified color, nothing to update"
        ))
    }
}

impl TypedNumeratorState for ColorSymplified {
    fn apply<F, S: TypedNumeratorState>(
        num: &mut Num<PythonState>,
        mut f: F,
    ) -> Result<(), NumeratorStateError>
    where
        F: FnMut(Num<Self>) -> Num<S>,
    {
        if let PythonState::ColorSymplified(s) = &mut num.state {
            if let Some(s) = s.take() {
                *num = f(Num { state: s }).forget_type();
                return Ok(());
            } else {
                return Err(NumeratorStateError::NoneVariant);
            }
        }
        Err(NumeratorStateError::NotColorSymplified)
    }
}

impl UnexpandedNumerator for ColorSymplified {
    fn expr(&self) -> Result<SerializableAtom, NumeratorStateError> {
        Ok(self.expression.clone())
    }
}

impl Num<ColorSymplified> {
    pub fn gamma_symplify(self) -> Num<GammaSymplified> {
        Num {
            state: self.state.gamma_symplify(),
        }
    }

    pub fn parse(self) -> Num<Network> {
        Num {
            state: self.state.parse(),
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
pub struct GammaSymplified {
    #[bincode(with_serde)]
    expression: SerializableAtom,
}

impl TryFrom<PythonState> for GammaSymplified {
    type Error = NumeratorStateError;

    fn try_from(value: PythonState) -> std::result::Result<Self, Self::Error> {
        match value {
            PythonState::GammaSymplified(s) => {
                if let Some(s) = s {
                    Ok(s)
                } else {
                    Err(NumeratorStateError::NoneVariant)
                }
            }
            _ => Err(NumeratorStateError::NotGammaSymplified),
        }
    }
}

impl NumeratorState for GammaSymplified {
    fn export(&self) -> String {
        self.expression.0.to_string()
    }

    fn forget_type(self) -> PythonState {
        PythonState::GammaSymplified(Some(self))
    }

    fn update_model(&mut self, model: &Model) -> Result<()> {
        Err(eyre!(
            "Only applied feynman rule, simplified color and gamma, nothing to update"
        ))
    }
}

impl GammaSymplified {
    pub fn parse(self) -> Network {
        let net = TensorNetwork::try_from(self.expression.0.as_view())
            .unwrap()
            .to_fully_parametric()
            .cast();

        // println!("net scalar{}", net.scalar.as_ref().unwrap());
        Network { net }
    }
}

impl UnexpandedNumerator for GammaSymplified {
    fn expr(&self) -> Result<SerializableAtom, NumeratorStateError> {
        Ok(self.expression.clone())
    }
}

impl TypedNumeratorState for GammaSymplified {
    fn apply<F, S: TypedNumeratorState>(
        num: &mut Num<PythonState>,
        mut f: F,
    ) -> Result<(), NumeratorStateError>
    where
        F: FnMut(Num<Self>) -> Num<S>,
    {
        if let PythonState::GammaSymplified(s) = &mut num.state {
            if let Some(s) = s.take() {
                *num = f(Num { state: s }).forget_type();
                return Ok(());
            } else {
                return Err(NumeratorStateError::NoneVariant);
            }
        }
        Err(NumeratorStateError::NotGammaSymplified)
    }
}

impl Num<GammaSymplified> {
    pub fn parse(self) -> Num<Network> {
        Num {
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
                // println!("contracted tensor: {}", tensor);
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

    fn update_model(&mut self, model: &Model) -> Result<()> {
        Err(eyre!(
            "Only applied feynman rule, simplified color, gamma and parsed into network, nothing to update"
        ))
    }
}

impl TypedNumeratorState for Network {
    fn apply<F, S: TypedNumeratorState>(
        num: &mut Num<PythonState>,
        mut f: F,
    ) -> Result<(), NumeratorStateError>
    where
        F: FnMut(Num<Self>) -> Num<S>,
    {
        if let PythonState::Network(s) = &mut num.state {
            if let Some(s) = s.take() {
                *num = f(Num { state: s }).forget_type();
                return Ok(());
            } else {
                return Err(NumeratorStateError::NoneVariant);
            }
        }
        Err(NumeratorStateError::NotNetwork)
    }
}

impl Num<Network> {
    pub fn contract<R>(self, settings: ContractionSettings<R>) -> Result<Num<Contracted>> {
        let contracted = self.state.contract(settings)?;
        Ok(Num { state: contracted })
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
pub struct Contracted {
    #[bincode(with_serde)]
    tensor: ParamTensor<AtomStructure>,
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
        export_settings: &ExportSettings,
        params: &[Atom],
    ) -> EvaluatorSingle {
        let mut fn_map: FunctionMap = FunctionMap::new();

        Num::<Contracted>::add_consts_to_fn_map(&mut fn_map);

        debug!("Generate eval tree set with {} params", params.len());

        let mut eval_tree = self.tensor.eval_tree(&fn_map, params).unwrap();
        debug!("Horner scheme");

        eval_tree.horner_scheme();
        debug!("Common subexpression elimination");
        eval_tree.common_subexpression_elimination();
        debug!("Linearize double");
        let cpe_rounds = export_settings.numerator_settings.cpe_rounds();
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
            .compile_options()
            .compile()
        {
            let mut filename = path.clone();
            filename.push("numerator_single.cpp");
            let filename = filename.to_string_lossy();

            let function_name = "numerator_single";

            let library_name = path.join("numerator_single.so");
            let library_name = library_name.to_string_lossy();
            let inline_asm = export_settings.gammaloop_compile_options.inline_asm();

            let compile_options = export_settings
                .gammaloop_compile_options
                .to_symbolica_compile_options();

            CompiledEvaluator::new(
                eval.export_cpp(&filename, function_name, true, inline_asm)
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

    fn update_model(&mut self, model: &Model) -> Result<()> {
        Err(eyre!("Only applied feynman rule, simplified color, gamma, parsed into network and contracted, nothing to update"))
    }
}

impl TypedNumeratorState for Contracted {
    fn apply<F, S: TypedNumeratorState>(
        num: &mut Num<PythonState>,
        mut f: F,
    ) -> Result<(), NumeratorStateError>
    where
        F: FnMut(Num<Self>) -> Num<S>,
    {
        if let PythonState::Contracted(s) = &mut num.state {
            if let Some(s) = s.take() {
                *num = f(Num { state: s }).forget_type();
                return Ok(());
            } else {
                return Err(NumeratorStateError::NoneVariant);
            }
        }
        Err(NumeratorStateError::NotContracted)
    }
}

impl Num<Contracted> {
    pub fn generate_evaluators(
        self,
        model: &Model,
        graph: &BareGraph,
        extra_info: &ExtraInfo,
        export_settings: &ExportSettings,
    ) -> Num<Evaluators> {
        let (params, double_param_values, quad_param_values, model_params_start) =
            Contracted::generate_params(graph, model);

        let single = self
            .state
            .evaluator(extra_info.path.clone(), export_settings, &params);

        match export_settings.numerator_settings {
            NumeratorEvaluatorOptions::Combined(_) => Num {
                state: Evaluators {
                    orientated: Some(single.orientated(
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
                },
            },
            _ => Num {
                state: Evaluators {
                    orientated: None,
                    single,
                    choice: SingleOrCombined::Single,
                    orientations: extra_info.orientations.clone(),
                    quad_param_values,
                    double_param_values,
                    model_params_start,
                },
            },
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
#[serde(tag = "type")]
pub enum NumeratorEvaluatorOptions {
    #[serde(rename = "Single")]
    Single(EvaluatorOptions),
    #[serde(rename = "Combined")]
    Combined(EvaluatorOptions),
}

impl Default for NumeratorEvaluatorOptions {
    fn default() -> Self {
        NumeratorEvaluatorOptions::Combined(EvaluatorOptions::default())
    }
}

impl NumeratorEvaluatorOptions {
    pub fn compile_options(&self) -> NumeratorCompileOptions {
        match self {
            NumeratorEvaluatorOptions::Single(options) => options.compile_options,
            NumeratorEvaluatorOptions::Combined(options) => options.compile_options,
        }
    }

    pub fn cpe_rounds(&self) -> Option<usize> {
        match self {
            NumeratorEvaluatorOptions::Single(options) => options.cpe_rounds,
            NumeratorEvaluatorOptions::Combined(options) => options.cpe_rounds,
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
    pub eval_double: EvalTensorSet<ExpressionEvaluator<Complex<F<f64>>>, AtomStructure>,
    pub eval_quad: EvalTensorSet<ExpressionEvaluator<Complex<F<f128>>>, AtomStructure>,
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
    pub fn orientated(
        &self,
        graph: &BareGraph,
        params: &[Atom],
        extra_info: &ExtraInfo,
        export_settings: &ExportSettings,
    ) -> EvaluatorOrientations {
        let mut fn_map: FunctionMap = FunctionMap::new();

        Num::<Contracted>::add_consts_to_fn_map(&mut fn_map);

        let mut seen = 0;

        let mut index_map = IndexSet::new();
        let mut positions = Vec::new();

        let len = extra_info.orientations.len();

        let reps = graph
            .edges
            .iter()
            .enumerate()
            .map(|(i, _)| {
                (
                    Pattern::parse(&format!("Q({},cind(1))", i)).unwrap(),
                    Pattern::parse(&format!("-Q({},cind(1))", i)).unwrap(),
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

        let cpe_rounds = export_settings.numerator_settings.cpe_rounds();

        debug!("Generate eval tree set with {} params", params.len());

        let mut eval_tree = set.eval_tree(&fn_map, &params).unwrap();
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
            .compile_options()
            .compile()
        {
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

            let mut filename = extra_info.path.clone();
            filename.push("numerator.cpp");
            let filename = filename.to_string_lossy();

            let function_name = "numerator";

            let library_name = extra_info.path.join("numerator.so");
            let library_name = library_name.to_string_lossy();
            let inline_asm = export_settings.gammaloop_compile_options.inline_asm();

            let compile_options = export_settings
                .gammaloop_compile_options
                .to_symbolica_compile_options();
            CompiledEvaluator::new(
                eval.export_cpp(&filename, function_name, true, inline_asm)
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
    state: CompiledState,
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
            self.double_param_values[self.model_params_start + i] = v.map(|f| f.into());
            self.quad_param_values[self.model_params_start + i] = v.map(|f| f.higher().into());
        }
        Ok(())
    }
}

impl TypedNumeratorState for Evaluators {
    fn apply<F, S: TypedNumeratorState>(
        num: &mut Num<PythonState>,
        mut f: F,
    ) -> Result<(), NumeratorStateError>
    where
        F: FnMut(Num<Self>) -> Num<S>,
    {
        if let PythonState::Evaluators(s) = &mut num.state {
            if let Some(s) = s.take() {
                *num = f(Num { state: s }).forget_type();
                return Ok(());
            } else {
                return Err(NumeratorStateError::NoneVariant);
            }
        }
        Err(NumeratorStateError::NotEvaluators)
    }
}

impl Num<Evaluators> {
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
                    .orientated(graph, &params, extra_info, export_settings);
            self.state.orientated = Some(orientated);
            self.state.choice = SingleOrCombined::Combined;
        } else {
            panic!("Cannot enable combined without generating the evaluators")
        }
    }
}

use thiserror::Error;

#[derive(Error, Debug)]
pub enum NumeratorStateError {
    #[error("Not UnInit")]
    NotUnit,
    #[error("Not AppliedFeynmanRule")]
    NotAppliedFeynmanRule,
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
    AppliedFeynmanRule(Option<AppliedFeynmanRule>),
    ColorSymplified(Option<ColorSymplified>),
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
    fn expr(&self) -> Result<SerializableAtom, NumeratorStateError> {
        match self {
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
