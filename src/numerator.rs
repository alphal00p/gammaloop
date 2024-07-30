use std::fmt::{Debug, Pointer};

use crate::graph::{BareGraph, Edge};
use crate::utils::f128;
use crate::{
    graph::{EdgeType, LoopMomentumBasis},
    model::Model,
    momentum::FourMomentum,
    utils::{FloatLike, F},
};
use ahash::AHashMap;
use eyre::{eyre, Result};
use itertools::Itertools;
use log::debug;
use serde::{ser::SerializeStruct, Deserialize, Serialize};
use spenso::{
    complex::Complex,
    network::TensorNetwork,
    parametric::{CompiledEvalTensor, EvalTreeTensor, ParamTensor, PatternReplacement},
    structure::{Lorentz, NamedStructure, PhysReps, RepName, Shadowable, TensorStructure},
    symbolic::SymbolicTensor,
};
use symbolica::{
    atom::AtomView,
    domains::{
        float::{Complex as SymComplex, NumericalFloatLike},
        rational::Rational,
    },
    evaluate::{CompiledEvaluator, EvalTree, FunctionMap},
    id::{Pattern, Replacement},
};
use symbolica::{
    atom::{Atom, FunctionBuilder, Symbol},
    printer::{AtomPrinter, PrintOptions},
    state::State,
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
pub enum EvaluatorNumerator {
    Compiled(TensorNetwork<CompiledEvalTensor<AtomStructure>, CompiledEvaluator>),
    Eager(TensorNetwork<EvalTreeTensor<Rational, AtomStructure>, EvalTree<Rational>>),
    UnInit,
}

impl Debug for EvaluatorNumerator {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            EvaluatorNumerator::Compiled(c) => c.fmt(f),
            EvaluatorNumerator::Eager(e) => e.fmt(f),
            EvaluatorNumerator::UnInit => write!(f, "uninit"),
        }
    }
}

impl Clone for EvaluatorNumerator {
    fn clone(&self) -> Self {
        EvaluatorNumerator::UnInit
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExtraInfo {
    edges: Vec<NumeratorEdge>,
    graph_name: String,
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
    pub network: Option<TensorNetwork<ParamTensor<NamedStructure<Symbol, Vec<Atom>>>, Atom>>,
    pub extra_info: ExtraInfo,
    pub const_map: AHashMap<Atom, Complex<F<f64>>>,
    pub eval: EvaluatorNumerator,
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
        let network: TensorNetwork<
            ParamTensor<NamedStructure<symbolica::atom::Symbol, Vec<Atom>>>,
            Atom,
        > = sym_tensor
            .to_network()
            .map_err(serde::de::Error::custom)?
            .to_fully_parametric();

        Ok(Numerator {
            expression,
            network: Some(network),
            const_map,
            extra_info: data.extra_info,
            eval: EvaluatorNumerator::UnInit,
        })
    }
}
pub type AtomStructure = NamedStructure<Symbol, Vec<Atom>>;

pub trait Evaluate<T: FloatLike> {
    fn evaluate(&mut self, emr: &[FourMomentum<F<T>>]) -> Result<Complex<F<T>>>;
}

// pub trait Compile<T: FloatLike> {
//     fn compile(&mut self, emr: &[FourMomentum<F<T>>]) -> Result<Complex<F<T>>>;
// }

impl Evaluate<f64> for Numerator {
    fn evaluate(&mut self, emr: &[FourMomentum<F<f64>>]) -> Result<Complex<F<f64>>> {
        let mut params: Vec<Complex<F<f64>>> = emr
            .iter()
            .flat_map(|&p| p.into_iter().map(Complex::new_re))
            .collect_vec();
        params.push(Complex {
            re: F(1.),
            im: F(0.),
        });
        match &mut self.eval {
            EvaluatorNumerator::Compiled(c) => Ok(c.evaluate(&params).result()?),
            EvaluatorNumerator::Eager(e) => Ok(e
                .map_coeff::<Complex<F<f64>>, _>(&|r| Complex {
                    re: F(r.into()),
                    im: F(r.zero().into()),
                })
                .linearize()
                .evaluate(&params)
                .result()?),
            EvaluatorNumerator::UnInit => Err(eyre!("Uninitialized numerator")),
        }
    }
}

impl Evaluate<f128> for Numerator {
    fn evaluate(&mut self, emr: &[FourMomentum<F<f128>>]) -> Result<Complex<F<f128>>> {
        let params: Vec<Complex<F<f128>>> = emr
            .iter()
            .flat_map(|p| p.into_iter().cloned().map(Complex::new_re))
            .collect_vec();
        match &mut self.eval {
            EvaluatorNumerator::Eager(e) => Ok(e
                .map_coeff::<Complex<F<f128>>, _>(&|r| Complex {
                    re: F(r.into()),
                    im: F(f128::new_zero()),
                })
                .linearize()
                .evaluate(&params)
                .result()?),
            EvaluatorNumerator::UnInit => Err(eyre!("Uninitialized numerator")),
            _ => Err(eyre!("Compiled eval not supported for f128")),
        }
    }
}

impl Numerator {
    pub fn substitute_model_params(&mut self, model: &Model) {
        self.expression = model.substitute_model_params(&self.expression);
    }

    pub fn build_const_fn_map_and_split(&mut self, fn_map: &mut FunctionMap, model: &Model) {
        // let mut replacements = vec![];
        // replacements.extend(model.dependent_coupling_replacements());
        // replacements.extend(model.internal_parameter_replacements());

        // let reps = replacements
        //     .iter()
        //     .map(|(lhs, rhs)| Replacement::new(lhs, rhs))
        //     .collect_vec();

        // if let Some(net) = &mut self.network {
        //     net.replace_all_multiple_repeat_mut(&reps);
        // }

        let mut split_reps = vec![];
        split_reps.extend(model.valued_coupling_re_im_split(fn_map));
        split_reps.extend(model.valued_parameter_re_im_split(fn_map));

        let reps = split_reps
            .iter()
            .map(|(lhs, rhs)| Replacement::new(lhs, rhs))
            .collect_vec();

        if let Some(net) = &mut self.network {
            net.replace_all_multiple_repeat_mut(&reps);
        }
    }

    pub fn compile<T: FloatLike + Default>(&mut self, model: &Model, graph: &BareGraph) {
        self.generate_evaluator(model);
        if let EvaluatorNumerator::Eager(eval) = &mut self.eval {
            self.eval = EvaluatorNumerator::Compiled(
                eval.map_coeff::<F<T>, _>(&|r| r.into())
                    .linearize()
                    .compile_asm(
                        &(self.extra_info.graph_name.clone() + "numerator_asm"),
                        &(self.extra_info.graph_name.clone() + "libneval_asm"),
                    ),
            );
        }
    }

    pub fn generate_evaluator(&mut self, model: &Model) {
        let mut fn_map: FunctionMap = FunctionMap::new();
        self.build_const_fn_map_and_split(&mut fn_map, model);

        let mut params: Vec<Atom> = vec![];

        for (i, _) in self.extra_info.edges.iter().enumerate() {
            let named_structure: NamedStructure<&str> = NamedStructure::from_iter(
                [PhysReps::new_slot(Lorentz {}.into(), 4, i)],
                "Q",
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

        params.push(Atom::new_var(State::I));

        if let Some(net) = self.network.as_mut() {
            net.contract();
            let mut eval_tree = net.eval_tree(|a| a.clone(), &fn_map, &params).unwrap();
            eval_tree.horner_scheme();
            eval_tree.common_subexpression_elimination(1);
            self.eval = EvaluatorNumerator::Eager(eval_tree);
        }
    }

    pub fn evaluate<T: FloatLike>(
        &self,
        emr: &[FourMomentum<F<T>>],
        graph: &BareGraph,
    ) -> Complex<F<T>> {
        let emr_params = graph
            .edges
            .iter()
            .enumerate()
            .map(|(i, _)| {
                let named_structure: NamedStructure<&str> = NamedStructure::from_iter(
                    [PhysReps::new_slot(Lorentz {}.into(), 4, i)],
                    "Q",
                    Some(i),
                );
                debug!("Q{}", i);
                named_structure.to_shell().expanded_shadow().unwrap()
            })
            .collect_vec();

        let mut constmap: AHashMap<AtomView<'_>, SymComplex<F<T>>> = AHashMap::new();

        let mut pols = Vec::new();
        for (i, ext) in graph
            .edges
            .iter()
            .enumerate()
            .filter(|(_, e)| !matches!(e.edge_type, EdgeType::Virtual))
        {
            match ext.edge_type {
                EdgeType::Incoming => {
                    for (a, p) in ext
                        .particle
                        .incoming_polarization_match(i, &emr[i])
                        .into_iter()
                    {
                        pols.push((a, p));
                    }
                }

                EdgeType::Outgoing => {
                    for (a, p) in ext
                        .particle
                        .outgoing_polarization_match(i, &emr[i])
                        .into_iter()
                    {
                        pols.push((a, p));
                    }
                }
                _ => {}
            }
        }

        for (a, p) in &pols {
            constmap.insert(a.as_view(), p.clone().into());
        }

        for (key, value) in self.const_map.iter() {
            let v = Complex::new(F::<T>::from_ff64(value.re), F::<T>::from_ff64(value.im));
            constmap.insert(key.as_view(), v.into());
        }

        for (i, pe) in emr.iter().enumerate() {
            constmap.insert(
                emr_params[i][1.into()].as_view(),
                Complex::new_re(pe.spatial.px.clone()).into(),
            );

            constmap.insert(
                emr_params[i][2.into()].as_view(),
                Complex::new_re(pe.spatial.py.clone()).into(),
            );

            constmap.insert(
                emr_params[i][3.into()].as_view(),
                Complex::new_re(pe.spatial.pz.clone()).into(),
            );

            constmap.insert(
                emr_params[i][0.into()].as_view(),
                Complex::new_re(pe.temporal.value.clone()).into(),
            );
        }

        // let mut file = File::create("serializable_map.json").unwrap();

        // debug!("Serializable map: {}", serialized);
        // file.write_all(serialized.as_bytes()).unwrap();

        // println!(
        //     "Numerator evaluated: {}",
        //     numerator_network.scalar.as_ref().unwrap()
        // );

        let serializable_map: AHashMap<String, Complex<F<T>>> = constmap
            .iter()
            .map(|(k, v)| (k.to_string(), v.clone().into()))
            .collect();

        let serialized = serde_json::to_string(&serializable_map).unwrap();
        debug!("Serializable map: {:?}", serialized);

        let mut evaluated = self
            .network
            .as_ref()
            .unwrap()
            .evaluate::<SymComplex<F<T>>, _>(|c| c.into(), &constmap);

        evaluated.contract();

        evaluated.result().unwrap().into()
    }

    pub fn process(&mut self, model: &Model, graph: &BareGraph) {
        let mut const_map = AHashMap::new();
        model.append_parameter_map(&mut const_map);

        let i_float = Complex::new_i();
        let i = Atom::new_var(State::I);
        const_map.insert(i, i_float);
        self.const_map.extend(const_map);

        self.substitute_model_params(model);
        self.process_color_simple();
        self.fill_network();
        self.generate_evaluator(model);
        self.compile::<f64>(model, graph);
    }

    // pub fn generate_fn_map(&self mut){
    //     let mut replacements = vec![];
    //     let mut fn_map: FunctionMap<Rational> = FunctionMap::new();

    //     for (k, v) in const_atom_map.iter() {
    //     let name_re = Atom::new_var(State::get_symbol(k.to_string() + "_re"));
    //     let name_im = Atom::new_var(State::get_symbol(k.to_string() + "_im"));
    //     let i = Atom::new_var(State::I);
    //     let pat = &name_re + i * &name_im;
    //     replacements.push((Atom::new_var(*k).into_pattern(), pat.into_pattern()));

    //     fn_map.add_constant(name_re.into(), Rational::from(v.re));
    //     fn_map.add_constant(name_im.into(), Rational::from(v.im));
    // }
    // }

    pub fn fill_network(&mut self) {
        let sym_tensor: SymbolicTensor = self.expression.clone().try_into().unwrap();

        let network = sym_tensor.to_network().unwrap().to_fully_parametric();
        self.network = Some(network);
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

    pub fn generate(graph: &mut BareGraph) -> Self {
        debug!("generating numerator for graph: {}", graph.name);

        let vatoms: Vec<Atom> = graph
            .vertices
            .iter()
            .enumerate()
            .flat_map(|(i, v)| v.contracted_vertex_rule(&graph))
            .collect();
        // graph
        //     .edges
        //     .iter()
        //     .filter(|e| e.edge_type == EdgeType::Virtual)
        //     .map(|e| e.particle);
        let mut eatoms: Vec<Atom> = vec![];
        let mut shift = 0;
        for e in &graph.edges {
            let (n, i) = e.numerator(&graph);
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

        debug!("numerator: {}", builder);

        let expression = builder.clone();

        Numerator {
            expression,
            network: None,
            const_map: AHashMap::new(),
            eval: EvaluatorNumerator::UnInit,
            extra_info: ExtraInfo {
                edges: graph.edges.iter().map(NumeratorEdge::from).collect(),
                graph_name: graph.name.clone().into(),
            },
        }
    }
}

#[cfg(test)]
mod tests;
