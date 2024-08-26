use std::fmt::Debug;
use std::path::PathBuf;
use std::time::Instant;

use crate::graph::{BareGraph, Edge};
use crate::momentum::Polarization;
use crate::utils::f128;
use crate::ExportSettings;
use crate::{
    graph::{EdgeType, LoopMomentumBasis},
    model::Model,
    momentum::FourMomentum,
    utils::{FloatLike, F},
};
use ahash::AHashMap;
use bincode::{Decode, Encode};
use color_eyre::Report;
use eyre::{eyre, Result};
use gat_lending_iterator::LendingIterator;
// use gxhash::GxBuildHasher;
use indexmap::IndexSet;
use itertools::Itertools;

use log::{debug, info};
use serde::de::DeserializeOwned;
use serde::{ser::SerializeStruct, Deserialize, Serialize};
use spenso::data::DataTensor;

use spenso::network::Levels;
use spenso::parametric::{
    CompiledEvalTensorSet, EvalTensor, EvalTensorSet, EvalTreeTensorSet, ParamTensorSet,
    SerializableAtom, SerializableCompiledEvaluator, TensorSet,
};
use spenso::structure::{HasStructure, SerializableSymbol, SmartShadowStructure};
use spenso::{
    complex::Complex,
    network::TensorNetwork,
    parametric::{ParamTensor, PatternReplacement},
    structure::{Lorentz, NamedStructure, PhysReps, RepName, Shadowable, TensorStructure},
    symbolic::SymbolicTensor,
};
use symbolica::evaluate::{CompileOptions, ExpressionEvaluator, InlineASM};
use symbolica::{
    atom::{Atom, FunctionBuilder},
    printer::{AtomPrinter, PrintOptions},
    state::State,
};
use symbolica::{
    domains::{float::NumericalFloatLike, rational::Rational},
    evaluate::FunctionMap,
    id::{Pattern, Replacement},
};

pub fn apply_replacements(
    graph: &BareGraph,
    model: &Model,
    lmb: &LoopMomentumBasis,
    mut atom: Atom,
) -> Atom {
    atom = model.substitute_model_params(&atom);

    for edge in &graph.edges {
        atom = edge.substitute_lmb(atom, graph, lmb);
    }
    atom
}

#[allow(clippy::large_enum_variant)]
pub enum CompiledNumerator {
    Compiled(CompiledEvalTensorSet<AtomStructure>),
    UnInit,
}

#[allow(clippy::large_enum_variant, clippy::type_complexity)]
pub enum EagerNumerator<T: FloatLike> {
    Eager(EvalTensorSet<ExpressionEvaluator<Complex<F<T>>>, AtomStructure>),
    UnInit,
}

pub enum EvalNumerator {
    Eval(EvalTreeTensorSet<Rational, AtomStructure>),
    UnInit,
}

impl<T: FloatLike> Debug for EagerNumerator<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            EagerNumerator::Eager(_e) => write!(f, "eager"),
            EagerNumerator::UnInit => write!(f, "uninit"),
        }
    }
}

impl Debug for CompiledNumerator {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CompiledNumerator::Compiled(e) => e.fmt(f),
            CompiledNumerator::UnInit => write!(f, "uninit"),
        }
    }
}

impl Debug for EvalNumerator {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            EvalNumerator::Eval(_e) => write!(f, "eval"),
            EvalNumerator::UnInit => write!(f, "uninit"),
        }
    }
}

impl Clone for CompiledNumerator {
    fn clone(&self) -> Self {
        CompiledNumerator::UnInit
    }
}
impl Clone for EvalNumerator {
    fn clone(&self) -> Self {
        EvalNumerator::UnInit
    }
}

impl<T: FloatLike> Clone for EagerNumerator<T> {
    fn clone(&self) -> Self {
        EagerNumerator::<T>::UnInit
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExtraInfo {
    pub path: PathBuf,
    pub orientations: Vec<Vec<bool>>,
}

impl From<&Edge> for NumeratorEdge {
    fn from(value: &Edge) -> Self {
        NumeratorEdge {
            edge_type: value.edge_type,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NumeratorEdge {
    pub edge_type: EdgeType,
}

#[derive(Debug, Clone)]
#[allow(clippy::type_complexity)]
pub struct Numerator {
    pub expression: Atom,
    pub network: Option<TensorNetwork<ParamTensor<AtomStructure>, SerializableAtom>>,
    pub extra_info: ExtraInfo,
    pub const_map: AHashMap<Atom, Complex<F<f64>>>,
    pub base_eval: EvalNumerator,
    pub positions: Option<Vec<usize>>,
    pub eval_double: EagerNumerator<f64>,
    pub eval_quad: EagerNumerator<f128>,
    pub compiled: CompiledNumerator,
}

pub struct NumeratorEvaluator {
    pub compiled: CompiledNumerator,
    pub eval_double: EagerNumerator<f64>,
    pub eval_quad: EagerNumerator<f128>,
}

impl Serialize for Numerator {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let expression =
            AtomPrinter::new_with_options(self.expression.as_view(), PrintOptions::file())
                .to_string();

        let const_map: AHashMap<String, Complex<F<f64>>> = self
            .const_map
            .iter()
            .map(|(k, &v)| {
                (
                    AtomPrinter::new_with_options(k.as_view(), PrintOptions::file()).to_string(),
                    v,
                )
            })
            .collect();

        let mut state = serializer.serialize_struct("Numerator", 3)?;
        state.serialize_field("expression", &expression)?;
        state.serialize_field("const_map", &const_map)?;
        state.serialize_field("extra_info", &self.extra_info)?;
        state.end()
    }
}

impl<'de> Deserialize<'de> for Numerator {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        #[derive(Deserialize)]
        struct NumeratorData {
            expression: String,
            extra_info: ExtraInfo,
            const_map: AHashMap<String, Complex<F<f64>>>,
        }

        let data = NumeratorData::deserialize(deserializer)?;

        let expression = Atom::parse(&data.expression).map_err(serde::de::Error::custom)?;

        let const_map: AHashMap<Atom, Complex<F<f64>>> = data
            .const_map
            .into_iter()
            .map(|(k, v)| (Atom::parse(&k).unwrap(), v))
            .collect();

        let sym_tensor: SymbolicTensor = expression
            .clone()
            .try_into()
            .map_err(serde::de::Error::custom)?;
        let network: TensorNetwork<ParamTensor<AtomStructure>, SerializableAtom> = sym_tensor
            .to_network()
            .map_err(serde::de::Error::custom)?
            .to_fully_parametric()
            .cast();

        Ok(Numerator {
            expression,
            network: Some(network),
            const_map,
            positions: None,
            extra_info: data.extra_info,
            base_eval: EvalNumerator::UnInit,
            eval_double: EagerNumerator::<f64>::UnInit,
            eval_quad: EagerNumerator::<f128>::UnInit,
            compiled: CompiledNumerator::UnInit,
        })
    }
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
        <T as NumeratorEvaluateFloat>::evaluate_all_orientations(
            self,
            &<T as NumeratorEvaluateFloat>::params(emr, polarizations),
        )
    }

    fn evaluate_single(
        &mut self,
        emr: &[FourMomentum<F<T>>],
        polarizations: &[Polarization<Complex<F<T>>>],
    ) -> DataTensor<Complex<F<T>>, AtomStructure> {
        <T as NumeratorEvaluateFloat>::evaluate_single(
            self,
            &<T as NumeratorEvaluateFloat>::params(emr, polarizations),
        )
    }
}

pub trait NumeratorEvaluateFloat<T: FloatLike = Self> {
    fn evaluate_all_orientations(
        num: &mut Num<Evaluators>,
        params: &[Complex<F<T>>],
    ) -> Result<RepeatingIteratorTensorOrScalar<DataTensor<Complex<F<T>>, AtomStructure>>>;

    fn evaluate_single(
        num: &mut Num<Evaluators>,
        params: &[Complex<F<T>>],
    ) -> DataTensor<Complex<F<T>>, AtomStructure>;

    fn params(
        emr: &[FourMomentum<F<T>>],
        polarizations: &[Polarization<Complex<F<T>>>],
    ) -> Vec<Complex<F<T>>>;
}

// pub trait Compile<T: FloatLike> {
//     fn compile(&mut self, emr: &[FourMomentum<F<T>>]) -> Result<Complex<F<T>>>;
// }

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

// impl<T: Clone> Iterator for RepeatingIterator<T> {
//     type Item = T;
//     fn next(&mut self) -> Option<Self::Item> {
//         Some(self.elements[self.positions.next()?].clone())
//     }
// }

impl<T> LendingIterator for RepeatingIterator<T> {
    type Item<'a> = &'a T where Self:'a ;

    fn next(&mut self) -> Option<Self::Item<'_>> {
        Some(&self.elements[self.positions.next()?])
    }
}

impl NumeratorEvaluateFloat for f64 {
    fn params(
        emr: &[FourMomentum<F<Self>>],
        polarizations: &[Polarization<Complex<F<Self>>>],
    ) -> Vec<Complex<F<Self>>> {
        let mut params: Vec<Complex<F<Self>>> = emr
            .iter()
            .flat_map(|p| p.into_iter().cloned().map(Complex::new_re))
            .collect_vec();

        for p in polarizations {
            for &pi in p {
                params.push(pi)
            }
        }

        let zero = f64::new_zero();
        params.push(Complex {
            im: F(zero.one()),
            re: F(zero),
        });

        params
    }

    fn evaluate_all_orientations(
        num: &mut Num<Evaluators>,
        params: &[Complex<F<Self>>],
    ) -> Result<RepeatingIteratorTensorOrScalar<DataTensor<Complex<F<Self>>, AtomStructure>>> {
        if let Some(orientation_evaluator) = &mut num.state.orientated {
            let pos = orientation_evaluator.positions.clone();
            let res = if let Some(c) = &mut orientation_evaluator.compiled.evaluator {
                c.evaluate(params)
            } else {
                orientation_evaluator.eval_double.evaluate(params)
            };
            Ok((res, pos).into())
        } else {
            Err(eyre!("No multi-orientation evaluator for numerator"))
        }
    }

    fn evaluate_single(
        num: &mut Num<Evaluators>,
        params: &[Complex<F<Self>>],
    ) -> DataTensor<Complex<F<Self>>, AtomStructure> {
        if let Some(c) = &mut num.state.single.compiled.evaluator {
            c.evaluate(params)
        } else {
            num.state.single.eval_double.evaluate(params)
        }
    }
}

impl NumeratorEvaluateFloat for f128 {
    fn params(
        emr: &[FourMomentum<F<Self>>],
        polarizations: &[Polarization<Complex<F<Self>>>],
    ) -> Vec<Complex<F<Self>>> {
        let mut params: Vec<Complex<F<Self>>> = emr
            .iter()
            .flat_map(|p| p.into_iter().cloned().map(Complex::new_re))
            .collect_vec();

        for p in polarizations {
            for pi in p {
                params.push(pi.clone())
            }
        }

        let zero = f128::new_zero();
        params.push(Complex {
            im: F(zero.one()),
            re: F(zero),
        });

        params
    }

    fn evaluate_all_orientations(
        num: &mut Num<Evaluators>,
        params: &[Complex<F<Self>>],
    ) -> Result<RepeatingIteratorTensorOrScalar<DataTensor<Complex<F<Self>>, AtomStructure>>> {
        if let Some(orientation_evaluator) = &mut num.state.orientated {
            let pos = orientation_evaluator.positions.clone();
            let res = orientation_evaluator.eval_quad.evaluate(params);
            Ok((res, pos).into())
        } else {
            Err(eyre!("No multi-orientation evaluator for numerator"))
        }
    }

    fn evaluate_single(
        num: &mut Num<Evaluators>,
        params: &[Complex<F<Self>>],
    ) -> DataTensor<Complex<F<Self>>, AtomStructure> {
        num.state.single.eval_quad.evaluate(params)
    }
}

impl Numerator {
    pub fn substitute_model_params(&mut self, model: &Model) {
        self.expression = model.substitute_model_params(&self.expression);
    }

    pub fn build_const_fn_map_and_split(&mut self, fn_map: &mut FunctionMap, model: &Model) {
        let mut split_reps = vec![];
        info!("splitting");
        split_reps.extend(model.valued_coupling_re_im_split(fn_map));
        split_reps.extend(model.valued_parameter_re_im_split(fn_map));

        let reps = split_reps
            .iter()
            .map(|(lhs, rhs)| Replacement::new(lhs, rhs))
            .collect_vec();

        fn_map.add_constant(Atom::parse("Nc").unwrap(), 3.into());

        fn_map.add_constant(Atom::parse("TR").unwrap(), Rational::from((1, 2)));

        fn_map.add_constant(
            Atom::parse("pi").unwrap(),
            Rational::from(std::f64::consts::PI),
        );

        if let Some(net) = &mut self.network {
            net.replace_all_multiple_repeat_mut(&reps);
        }
    }

    #[allow(dead_code)]
    fn replace_repeat(&mut self, lhs: Pattern, rhs: Pattern) {
        let atom = self.expression.replace_all(&lhs, &rhs, None, None);
        if atom != self.expression {
            self.expression = atom;
            self.replace_repeat(lhs, rhs);
        }
    }

    fn replace_repeat_multiple(&mut self, reps: &[Replacement<'_>]) {
        let atom = self.expression.replace_all_multiple(reps);
        // info!("expanded rep");
        if atom != self.expression {
            // info!("applied replacement");
            self.expression = atom;
            self.replace_repeat_multiple(reps);
        }
    }

    #[allow(dead_code)]
    fn replace_repeat_multiple_atom(expr: &mut Atom, reps: &[Replacement<'_>]) {
        let atom = expr.replace_all_multiple(reps);
        if atom != *expr {
            *expr = atom;
            Self::replace_repeat_multiple_atom(expr, reps)
        }
    }

    fn replace_repeat_multiple_atom_expand(expr: &mut Atom, reps: &[Replacement<'_>]) {
        let a = expr.expand();
        let atom = a.replace_all_multiple(reps);
        if atom != *expr {
            *expr = atom;
            Self::replace_repeat_multiple_atom_expand(expr, reps)
        }
    }

    pub fn isolate_color(&mut self) {
        let color_fn = FunctionBuilder::new(State::get_symbol("color"))
            .add_arg(&Atom::new_num(1))
            .finish();
        self.expression = &self.expression * color_fn;
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

        self.replace_repeat_multiple(&reps)
    }

    pub fn process_color_simple(&mut self) {
        self.isolate_color();
        let (mut coefs, rem) = self.expression.coefficient_list(State::get_symbol("color"));

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
            Self::replace_repeat_multiple_atom_expand(key, &reps);

            atom = atom + coef.factor() * key.factor();
            // println!("coef {i}:{}\n", coef.factor());
        }
        atom = atom + rem;
        self.expression = atom;

        // println!("rem:{}", rem);

        // self.replace_repeat_multiple(&reps);
    }

    #[allow(dead_code)]
    fn process_color(&mut self) {
        let idlhs = Pattern::parse("T(aind(y___,a_,z___))*id(aind(a_,b_))").unwrap();

        let idrhs = Pattern::parse("T(aind(y___,b_,z___))").unwrap();

        self.replace_repeat(idlhs, idrhs);

        // // T(a,...,i,j)T(b,...,j,k) = T(a,...,b,...,i,k)
        // let contractfundlhs = Pattern::parse("T(aind(y___,i_,j_))*T(aind(z___,j_,k_))").unwrap();

        // let contractfundrhs = Pattern::parse("T(aind(y___,z___,i_,k_))").unwrap();

        // self.replace_repeat(contractfundlhs, contractfundrhs);

        // //T(a,x,b,i,j)T(c,x,d,k,l) = 1/2(T(a,d,i,l)T(c,b,k,j)
        // //-1/Nc T(a,b,i,j)T(c,d,k,l))
        // let contractadjlhs =
        //     Pattern::parse("T(aind(a___,x__,b___,i_,j_))T(aind(c___,x__,d___,k_,l_))").unwrap();

        // let contractadjrhs =
        //     Pattern::parse("1/2(T(aind(a___,d___,i_,l_))T(aind(c___,b___,k_,j_))-1/Nc T(aind(a___,b___,i_,j_))T(aind(c___,d___,k_,l_)))").unwrap();

        // self.replace_repeat(contractadjlhs, contractadjrhs);

        // //T(a,b,c,...,i,i) = Tr(a,b,c,...)
        // let tracefundlhs = Pattern::parse("T(aind(a___,i_,i_)").unwrap();

        // let tracefundrhs = Pattern::parse("Tr(a___)").unwrap();

        // self.replace_repeat(tracefundlhs, tracefundrhs);

        // //T(a,x,b,x,c,i,j) = 1/2(T(a,c,i,j)Tr(b)-1/Nc T(a,b,c,i,j))

        // let traceadjlhs = Pattern::parse("T(a_,x_,b_,x_,c_,i_,j_)").unwrap();

        // let traceadjrhs =
        //     Pattern::parse("1/2(T(a_,c_,i_,j_)Tr(b_)-1/Nc T(aind(a_,b_,c_,i_,j_)))").unwrap();

        // self.replace_repeat(traceadjlhs, traceadjrhs);
    }

    // pub fn generate(graph: &mut BareGraph, orientations: Vec<Vec<bool>>) -> Self {
    //     debug!("generating numerator for graph: {}", graph.name);

    //     let vatoms: Vec<Atom> = graph
    //         .vertices
    //         .iter()
    //         .flat_map(|v| v.contracted_vertex_rule(graph))
    //         .collect();

    //     let mut eatoms: Vec<Atom> = vec![];
    //     let mut shift = 0;
    //     for e in &graph.edges {
    //         let (n, i) = e.numerator(graph);
    //         eatoms.push(n);
    //         shift += i;
    //         graph.shifts.0 += shift;
    //     }
    //     let mut builder = Atom::new_num(1);

    //     for v in &vatoms {
    //         builder = builder * v;
    //     }

    //     for e in &eatoms {
    //         builder = builder * e;
    //     }

    //     let i = Atom::new_var(State::I);
    //     let a = Atom::new_var(State::get_symbol("a_"));
    //     let b = Atom::new_var(State::get_symbol("b_"));

    //     let complex = FunctionBuilder::new(State::get_symbol("complex"))
    //         .add_arg(&a)
    //         .add_arg(&b)
    //         .finish();

    //     builder = complex.into_pattern().replace_all(
    //         builder.as_view(),
    //         &(&a + &b * &i).into_pattern(),
    //         None,
    //         None,
    //     );

    //     let expression = builder.clone();

    //     Numerator {
    //         expression,
    //         network: None,
    //         positions: None,
    //         const_map: AHashMap::new(),
    //         base_eval: EvalNumerator::UnInit,
    //         eval_double: EagerNumerator::<f64>::UnInit,
    //         eval_quad: EagerNumerator::<f128>::UnInit,
    //         compiled: CompiledNumerator::UnInit,
    //         extra_info: ExtraInfo {
    //             orientations,
    //             edges: graph.edges.iter().map(NumeratorEdge::from).collect(),
    //             graph_name: graph.name.clone().into(),
    //         },
    //     }
    // }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
pub struct Num<State> {
    pub state: State,
}

// impl<'de, State: NumeratorState> Deserialize<'de> for Num<State> {
//     fn deserialize<D>(deserializer: D) -> std::result::Result<Self, D::Error>
//     where
//         D: serde::Deserializer<'de>,
//     {
//         let state = State::deserialize(deserializer)?;
//         Ok(Num { state })
//     }
// }

impl<S: NumeratorState> Num<S> {
    pub fn export(&self) -> String {
        self.state.export()
    }

    pub fn forget_type(self) -> Num<PythonState> {
        Num {
            state: self.state.forget_type(),
        }
    }

    fn build_const_fn_map_with_split_reps(
        fn_map: &mut FunctionMap,
        model: &Model,
    ) -> Vec<(Pattern, Pattern)> {
        let mut split_reps = vec![];
        info!("splitting");
        split_reps.extend(model.valued_coupling_re_im_split(fn_map));
        split_reps.extend(model.valued_parameter_re_im_split(fn_map));

        fn_map.add_constant(Atom::parse("Nc").unwrap(), 3.into());

        fn_map.add_constant(Atom::parse("TR").unwrap(), Rational::from((1, 2)));

        fn_map.add_constant(
            Atom::parse("pi").unwrap(),
            Rational::from(std::f64::consts::PI),
        );

        split_reps
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
        for e in &graph.edges {
            let (n, i) = e.numerator(graph);
            eatoms.push(n);
            shift += i;
            graph.shifts.0 += shift;
        }
        let mut builder = Atom::new_num(1);

        for v in &vatoms {
            builder = builder * v;
        }

        for e in &eatoms {
            builder = builder * e;
        }

        let i = Atom::new_var(State::I);
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
    expression: SerializableAtom,
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
        GammaSymplified {
            expression: self.expression,
        }
    }

    pub fn parse(self) -> Network {
        Network {
            net: TensorNetwork::try_from(self.expression.0.as_view())
                .unwrap()
                .to_fully_parametric()
                .cast(),
        }
    }
}

impl NumeratorState for ColorSymplified {
    fn export(&self) -> String {
        self.expression.0.to_string()
    }

    fn forget_type(self) -> PythonState {
        PythonState::ColorSymplified(Some(self))
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
}

impl GammaSymplified {
    pub fn parse(self) -> Network {
        Network {
            net: TensorNetwork::try_from(self.expression.0.as_view())
                .unwrap()
                .to_fully_parametric()
                .cast(),
        }
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
                Ok(Contracted {
                    tensor: self.net.result_tensor_smart()?,
                })
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
        mut self,
        model: &Model,
        graph: &BareGraph,
        path: PathBuf,
        export_settings: &ExportSettings,
    ) -> EvaluatorSingle {
        let mut fn_map: FunctionMap = FunctionMap::new();

        let replacements =
            Num::<Contracted>::build_const_fn_map_with_split_reps(&mut fn_map, model);

        let reps = replacements
            .iter()
            .map(|(lhs, rhs)| Replacement::new(lhs, rhs))
            .collect_vec();

        self.tensor.replace_all_multiple_repeat_mut(&reps);

        let params = Self::generate_params(graph);

        debug!("Generate eval tree set with {} params", params.len());

        let mut eval_tree = self.tensor.eval_tree(&fn_map, &params).unwrap();
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

            println!("filename: {}", filename);

            let function_name = "numerator_single";

            let library_name = "libneval_single.so";
            let inline_asm = InlineASM::X64;
            CompiledEvaluator::new(
                eval.export_cpp(&filename, function_name, true, inline_asm)
                    .unwrap()
                    .compile(library_name, CompileOptions::default())
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
        }
    }

    pub fn generate_params(graph: &BareGraph) -> Vec<Atom> {
        let mut params: Vec<Atom> = vec![];
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

        params.extend(pols);
        params.push(Atom::new_var(State::I));
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
        let single = self
            .state
            .evaluator(model, graph, extra_info.path.clone(), export_settings);

        match export_settings.numerator_settings {
            NumeratorEvaluatorOptions::Combined(_) => Num {
                state: Evaluators {
                    orientated: Some(single.orientated(model, graph, extra_info, export_settings)),
                    single,
                },
            },
            _ => Num {
                state: Evaluators {
                    orientated: None,
                    single,
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
}

impl EvaluatorSingle {
    pub fn orientated(
        &self,
        model: &Model,
        graph: &BareGraph,
        extra_info: &ExtraInfo,
        export_settings: &ExportSettings,
    ) -> EvaluatorOrientations {
        let mut fn_map: FunctionMap = FunctionMap::new();

        let _ = Num::<Contracted>::build_const_fn_map_with_split_reps(&mut fn_map, model);

        let params = Contracted::generate_params(graph);
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
            let inline_asm = InlineASM::X64;
            CompiledEvaluator::new(
                eval.export_cpp(&filename, function_name, true, inline_asm)
                    .unwrap()
                    .compile(&library_name, CompileOptions::default())
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
    single: EvaluatorSingle,
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
