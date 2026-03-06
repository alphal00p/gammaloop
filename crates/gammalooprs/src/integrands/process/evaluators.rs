use std::{mem::transmute, ops::Neg, path::Path};

use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;
use eyre::eyre;
use linnet::half_edge::{
    involution::{EdgeVec, Orientation},
    subgraph::{SubSetIter, SubSetLike, subset::SubSet},
    typed_vec::IndexLike,
};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use spenso::algebra::{
    algebraic_traits::RefOne,
    complex::{Complex, symbolica_traits::CompiledComplexEvaluatorSpenso},
};
use symbolica::{
    atom::{Atom, AtomCore, FunctionBuilder, Indeterminate, Symbol},
    domains::{
        dual::HyperDual,
        float::Complex as SymComplex,
        integer::IntegerRing,
        rational::{Fraction, Rational},
    },
    evaluate::{
        CompileOptions, CompiledComplexEvaluator, Dualizer, ExportSettings, ExpressionEvaluator,
        FunctionMap, OptimizationSettings,
    },
};
use tracing::{debug, instrument};
use typed_index_collections::TiVec;

use crate::{
    GammaLoopContext,
    cff::expression::GraphOrientation,
    graph::Graph,
    integrands::process::{
        amplitude::load::set_override_if,
        param_builder::{FnMapEntry, LUParams},
    },
    momentum::Helicity,
    momentum::sample::MomentumSample,
    processes::EvaluatorSettings,
    settings::{GlobalSettings, RuntimeSettings},
    utils::{
        ArbPrec, F, FloatLike, GS, Length, W_, f128,
        symbolica_ext::{CallSymbol, LogPrint},
    },
};

use super::{
    ParamBuilder,
    param_builder::{ThresholdParams, UpdateAndGetParams},
};

#[derive(Clone, Copy)]
pub enum SingleOrAllOrientations<'a, OID> {
    Single {
        orientation: &'a EdgeVec<Orientation>,
        id: OID,
    },
    All {
        all: &'a TiVec<OID, EdgeVec<Orientation>>,
        filter: &'a SubSet<OID>,
    },
}
impl<'a, OID: IndexLike> SingleOrAllOrientations<'a, OID> {
    pub fn is_all(&self) -> bool {
        let SingleOrAllOrientations::All { filter, .. } = self else {
            return false;
        };
        (**filter).is_full()
    }
    pub fn iter(&self) -> SingleOrAllOrientationsIterator<'_, OID> {
        match self {
            SingleOrAllOrientations::All { all, filter } => SingleOrAllOrientationsIterator::All {
                all: *all,
                filter: (*filter).included_iter(),
            },
            SingleOrAllOrientations::Single { orientation, id } => {
                SingleOrAllOrientationsIterator::Single {
                    orientation: Some(orientation),
                    id: *id,
                }
            }
        }
    }
}

impl<OID> Length for SingleOrAllOrientations<'_, OID> {
    fn len(&self) -> usize {
        match self {
            SingleOrAllOrientations::Single { .. } => 1,
            SingleOrAllOrientations::All { all, .. } => all.len(),
        }
    }
}

#[derive(Clone)]
pub enum SingleOrAllOrientationsIterator<'a, OID> {
    Single {
        orientation: Option<&'a EdgeVec<Orientation>>,
        id: OID,
    },
    All {
        all: &'a TiVec<OID, EdgeVec<Orientation>>,
        filter: SubSetIter<'a, OID>,
    },
}

impl<'a, OID: IndexLike> Iterator for SingleOrAllOrientationsIterator<'a, OID>
where
    usize: From<OID>,
{
    type Item = (OID, &'a EdgeVec<Orientation>);

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            SingleOrAllOrientationsIterator::Single { orientation, id } => {
                orientation.take().map(|a| (*id, a))
            }
            SingleOrAllOrientationsIterator::All { all, filter } => {
                let a = filter.next()?;
                Some((a, &all[a]))
            }
        }
    }
}

#[derive(Clone, Debug, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
pub struct CompiledComplexEvaluatorGL(CompiledComplexEvaluator);

impl CompiledComplexEvaluatorGL {
    pub fn evaluate(&mut self, args: &[Complex<F<f64>>], out: &mut [Complex<F<f64>>]) {
        unsafe {
            self.0.evaluate(
                transmute::<&[Complex<F<f64>>], &[SymComplex<f64>]>(args),
                transmute::<&mut [Complex<F<f64>>], &mut [SymComplex<f64>]>(out),
            );
        }
    }
}
#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[derive(Debug, Clone, Deserialize, Serialize, Encode, Decode, PartialEq, JsonSchema)]
pub enum EvaluatorMethod {
    SingleParametric,
    Iterative,
    SummedFunctionMap,
    Summed,
}

impl Default for EvaluatorMethod {
    fn default() -> Self {
        EvaluatorMethod::SingleParametric
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct EvaluatorStack {
    pub single_parametric: GenericEvaluator,
    pub iterative: Option<(GenericEvaluator, usize)>,
    // pub iterative_function_map: Option<GenericEvaluator>,
    pub summed_function_map: Option<GenericEvaluator>,
    pub summed: Option<GenericEvaluator>,
}

impl EvaluatorStack {
    #[instrument(
        skip_all,
          fields(
           indicatif.pb_show = true, indicatif.pb_msg = "Generating Single Parametric Evaluator",
          )
      )]
    fn new_single_parametric<A: AtomCore>(
        parametric_atom: &[A],
        param_builder: &ParamBuilder,
        settings: &EvaluatorSettings,
    ) -> Result<GenericEvaluator> {
        let opt_settings = settings.optimization_settings();

        GenericEvaluator::new_from_builder(
            parametric_atom
                .iter()
                .map(|atom| GS.collect_orientation_if(atom.as_atom_view(), false)),
            &param_builder,
            None,
            opt_settings.clone(),
            settings.store_atom,
        )
    }
    #[instrument(
        skip_all,
          fields(
           indicatif.pb_show = true, indicatif.pb_msg = "Generating Iterative Evaluator",
          )
      )]
    fn new_iterative<A: AtomCore>(
        parametric_atom: &[A],
        param_builder: &ParamBuilder,
        orientations: &[EdgeVec<Orientation>],
        settings: &EvaluatorSettings,
    ) -> Result<(GenericEvaluator, usize)> {
        // let  n_orientations=;

        Ok((
            GenericEvaluator::new_from_builder(
                parametric_atom
                    .iter()
                    .map(|atom| {
                        orientations.iter().map(|a| {
                            let selected = a.select(atom.as_atom_view());
                            debug!(selected_expr = %selected.log_print(), "Iterative");
                            selected
                        })
                    })
                    .flatten(),
                &param_builder,
                None,
                settings.optimization_settings(),
                settings.store_atom,
            )?,
            orientations.len(),
        ))
    }

    #[instrument(
        skip_all,
          fields(
           indicatif.pb_show = true, indicatif.pb_msg = "Generating Summed Function Map Evaluator",
          )
      )]
    fn new_summed_function_map<A: AtomCore>(
        atoms: &[A],
        param_builder: &ParamBuilder,
        orientations: &[EdgeVec<Orientation>],
        settings: &EvaluatorSettings,
    ) -> Result<GenericEvaluator> {
        let params: Vec<Atom> = (&param_builder.pairs)
            .into_iter()
            .flat_map(|p| p.params.clone())
            .collect();
        let mut fn_map = param_builder.fn_map.clone();

        //I(sign(1), sign(2), sign(3),...) -> I(σ1, σ2, σ3,...)

        let entries: Vec<FnMapEntry> = atoms
            .iter()
            .enumerate()
            .map(|(i, a)| {
                let mut args = vec![];
                let mut lhs = FunctionBuilder::new(GS.integrand);
                lhs = lhs.add_arg(i);
                for (e, _) in &orientations[0] {
                    lhs = lhs.add_arg(GS.sign(e));
                    args.push(Indeterminate::try_from(GS.sign(e)).unwrap());
                }
                let param_integrand = GS.collect_orientation_if(a.as_atom_view(), false);
                fn_map
                    .add_tagged_function(
                        GS.integrand,
                        vec![Atom::num(i)],
                        "integrand".into(),
                        args.clone(),
                        param_integrand.clone(),
                    )
                    .map_err(|a| eyre!(a))?;
                Ok(FnMapEntry {
                    lhs: lhs.finish(),
                    rhs: param_integrand.clone(),
                    tags: vec![Atom::num(i)],
                    args: args.into_iter().map(|a| a.into()).collect(),
                })
            })
            .collect::<Result<_>>()?;

        // fn_map.add_conditional(name)

        let sum = (0..entries.len()).map(|i| {
            orientations
                .iter()
                .map(|a| {
                    GS.collect_orientation_if(a.orientation_thetas() * GS.integrand(i, a), true)
                    // GS.integrand(a)
                })
                .fold(Atom::Zero, |acc, n| acc + n)
                .replace(
                    Symbol::IF.f([GS.override_if, W_.b_, W_.c_])
                        + Symbol::IF.f([GS.override_if, W_.d_, W_.e_]),
                )
                .repeat()
                .with(Symbol::IF.f([
                    Atom::var(GS.override_if),
                    Atom::var(W_.b_) + Atom::var(W_.d_),
                    Atom::var(W_.c_) + Atom::var(W_.e_),
                ]))
        });

        GenericEvaluator::new_from_raw_params(
            sum,
            &params,
            &fn_map,
            entries,
            settings.optimization_settings(),
            None,
            settings.store_atom,
        )
    }

    #[instrument(
        skip_all,
          fields(
           indicatif.pb_show = true, indicatif.pb_msg = "Generating Summed Evaluator",
          )
      )]
    fn new_summed<A: AtomCore>(
        atoms: &[A],
        param_builder: &ParamBuilder,
        orientations: &[EdgeVec<Orientation>],
        settings: &EvaluatorSettings,
    ) -> Result<GenericEvaluator> {
        let sum = atoms.iter().map(|atom| {
            orientations
                .iter()
                .map(|a| {
                    let selected = GS.collect_orientation_if(
                        a.orientation_thetas() * a.select(atom.as_atom_view()),
                        true,
                    );
                    debug!(selected_expr = %selected.log_print(), "Iterative");
                    selected
                })
                .fold(Atom::Zero, |acc, n| acc + n)
                .replace(
                    Symbol::IF.f([GS.override_if, W_.b_, W_.c_])
                        + Symbol::IF.f([GS.override_if, W_.d_, W_.e_]),
                )
                .repeat()
                .with(Symbol::IF.f([
                    Atom::var(GS.override_if),
                    Atom::var(W_.b_) + Atom::var(W_.d_),
                    Atom::var(W_.c_) + Atom::var(W_.e_),
                ]))
        });

        GenericEvaluator::new_from_builder(
            sum,
            &param_builder,
            None,
            settings.optimization_settings(),
            settings.store_atom,
        )
    }
    #[instrument(
           skip_all,
           fields(indicatif.pb_show = true, indicatif.pb_msg = "Building Evaluator Stack"),
           err
       )]
    pub fn new<A: AtomCore>(
        atoms: &[A],
        param_builder: &ParamBuilder,
        orientations: &[EdgeVec<Orientation>],
        settings: &EvaluatorSettings,
    ) -> Result<Self> {
        let iterative = if settings.iterative_orientation_optimization {
            Some(Self::new_iterative(
                atoms,
                param_builder,
                orientations,
                settings,
            )?)
        } else {
            None
        };

        let summed_function_map = if settings.summed_function_map {
            Some(Self::new_summed_function_map(
                atoms,
                param_builder,
                orientations,
                settings,
            )?)
        } else {
            None
        };

        let summed = if settings.summed {
            Some(Self::new_summed(
                atoms,
                param_builder,
                orientations,
                settings,
            )?)
        } else {
            None
        };

        let single_parametric = Self::new_single_parametric(atoms, param_builder, settings)?;

        Ok(EvaluatorStack {
            single_parametric,
            iterative,
            summed_function_map,
            summed,
        })
    }

    fn evaluate_parametric<'a, T: FloatLike, OID: IndexLike>(
        &'a mut self,
        mut input: InputParams<'a, T>,
        orientations: SingleOrAllOrientations<'a, OID>,
    ) -> Vec<Complex<F<T>>>
    where
        usize: From<OID>,
    {
        let mut result: Option<Vec<Complex<F<T>>>> = None;
        for (_, e) in orientations.iter() {
            input.set_orientation_values(e);
            let output = <T as GenericEvaluatorFloat>::get_evaluator(&mut self.single_parametric)(
                input.as_slice(),
            );
            if let Some(result) = &mut result {
                for (r, v) in result.iter_mut().zip(output) {
                    *r += v;
                }
            } else {
                result = Some(output)
            }
        }
        result.unwrap()
    }

    fn evaluate_iterative<'a, T: FloatLike, OID: IndexLike>(
        &'a mut self,
        mut input: InputParams<'a, T>,
        orientations: SingleOrAllOrientations<'a, OID>,
    ) -> Result<Vec<Complex<F<T>>>>
    where
        usize: From<OID>,
    {
        let Some((iterative, len)) = &mut self.iterative else {
            return Err(eyre!(
                "Iterative evaluator not available. Regenerate with iterative set to true."
            ));
        };

        let mut result = vec![];
        if orientations.is_all() {
            let mut push = 0;
            let mut val = None;
            for r in <T as GenericEvaluatorFloat>::get_evaluator(iterative)(input.as_slice()) {
                push += 1;
                if let Some(mut value) = val.take() {
                    value += r;
                    if push == *len {
                        push = 0;
                        result.push(value);
                    } else {
                        val = Some(value)
                    }
                } else {
                    val = Some(r);
                }
            }
            Ok(result)
        } else {
            for (_, e) in orientations.iter() {
                input.set_orientation_values(e);

                let mut val = None;
                for (i, r) in
                    <T as GenericEvaluatorFloat>::get_evaluator(iterative)(input.as_slice())
                        .into_iter()
                        .enumerate()
                {
                    if let Some(mut value) = val.take() {
                        value += r;
                        if i % *len == 0 {
                            let pos = i / (*len);
                            if result.len() >= pos {
                                result.push(value);
                            } else {
                                result[pos] += value;
                            }
                            val = None;
                        } else {
                            val = Some(value)
                        }
                    } else {
                        val = Some(r);
                    }
                }
            }
            Ok(result)
        }
    }

    fn evaluate_summed_fnmap<'a, T: FloatLike, OID: IndexLike>(
        &'a mut self,
        mut input: InputParams<'a, T>,
        orientations: SingleOrAllOrientations<'a, OID>,
    ) -> Result<Vec<Complex<F<T>>>>
    where
        usize: From<OID>,
    {
        let Some(summed_function_map) = &mut self.summed_function_map else {
            return Err(eyre!(
                "Summed function map evaluator not available. Regenerate with summed_function_map set to true."
            ));
        };

        if orientations.is_all() {
            input.override_if(true);
            // if let Some(exprs) = &summed_function_map.exprs {
            //     for e in exprs {
            //         debug!(expr=%e.log_print(),"Summed evaluator");
            //     }
            // }

            Ok(<T as GenericEvaluatorFloat>::get_evaluator(
                summed_function_map,
            )(input.as_slice()))
        } else {
            let mut result: Option<Vec<Complex<F<T>>>> = None;
            for (_, e) in orientations.iter() {
                input.set_orientation_values(e);
                let output = <T as GenericEvaluatorFloat>::get_evaluator(summed_function_map)(
                    input.as_slice(),
                );
                if let Some(result) = &mut result {
                    for (r, v) in result.iter_mut().zip(output) {
                        *r += v;
                    }
                } else {
                    result = Some(output)
                }
            }
            Ok(result.unwrap())
        }
    }

    fn evaluate_summed<'a, T: FloatLike, OID: IndexLike>(
        &'a mut self,
        mut input: InputParams<'a, T>,
        orientations: SingleOrAllOrientations<'a, OID>,
    ) -> Result<Vec<Complex<F<T>>>>
    where
        usize: From<OID>,
    {
        let Some(summed) = &mut self.summed else {
            return Err(eyre!(
                "Summed evaluator not available. Regenerate with summed set to true."
            ));
        };

        if orientations.is_all() {
            input.override_if(true);

            Ok(<T as GenericEvaluatorFloat>::get_evaluator(summed)(
                input.as_slice(),
            ))
        } else {
            let mut result: Option<Vec<Complex<F<T>>>> = None;
            for (_, e) in orientations.iter() {
                input.set_orientation_values(e);
                let output = <T as GenericEvaluatorFloat>::get_evaluator(summed)(input.as_slice());
                if let Some(result) = &mut result {
                    for (r, v) in result.iter_mut().zip(output) {
                        *r += v;
                    }
                } else {
                    result = Some(output)
                }
            }
            Ok(result.unwrap())
        }
    }
    #[instrument(
        name = "evaluate",
        level = "info",
        skip(self, input, orientations, settings),
        fields(
            num_orientations = orientations.len(),
            method = ?settings.general.evaluator_method,
        )
    )]
    pub fn evaluate<'a, T: FloatLike, OID: IndexLike>(
        &'a mut self,
        input: InputParams<'a, T>,
        orientations: SingleOrAllOrientations<'a, OID>,
        settings: &RuntimeSettings,
    ) -> Result<Vec<Complex<F<T>>>>
    where
        usize: From<OID>,
    {
        match settings.general.evaluator_method {
            EvaluatorMethod::SingleParametric => Ok(self.evaluate_parametric(input, orientations)),
            EvaluatorMethod::Iterative => self.evaluate_iterative(input, orientations),
            EvaluatorMethod::SummedFunctionMap => self.evaluate_summed_fnmap(input, orientations),
            EvaluatorMethod::Summed => self.evaluate_summed(input, orientations),
        }
    }

    #[instrument(
          name = "compile",
          level = "info",
          skip(self, path,name, settings),
          fields(
              name = %name.as_ref(),
              path = %path.as_ref().display(),
          )
      )]
    pub fn compile(
        &mut self,
        name: impl AsRef<str>,
        path: impl AsRef<Path>,
        settings: &GlobalSettings,
    ) -> Result<()> {
        let name = name.as_ref();
        self.single_parametric.compile(
            path.as_ref().join(name).with_extension("cpp"),
            name,
            path.as_ref().join(name).with_extension("so"),
            settings.generation.compile.export_settings(),
        );

        if let Some((iterative, _)) = &mut self.iterative {
            iterative.compile(
                path.as_ref()
                    .join(format!("{}_iterative", name))
                    .with_extension("cpp"),
                &format!("{}_iterative", name),
                path.as_ref()
                    .join(format!("{}_iterative", name))
                    .with_extension("so"),
                settings.generation.compile.export_settings(),
            );
        }

        if let Some(summed_function_map) = &mut self.summed_function_map {
            summed_function_map.compile(
                path.as_ref()
                    .join(format!("{}_summed_function_map", name))
                    .with_extension("cpp"),
                &format!("{}_summed_function_map", name),
                path.as_ref()
                    .join(format!("{}_summed_function_map", name))
                    .with_extension("so"),
                settings.generation.compile.export_settings(),
            );
        }

        if let Some(summed) = &mut self.summed {
            summed.compile(
                path.as_ref()
                    .join(format!("{}_summed", name))
                    .with_extension("cpp"),
                &format!("{}_summed", name),
                path.as_ref()
                    .join(format!("{}_summed", name))
                    .with_extension("so"),
                settings.generation.compile.export_settings(),
            );
        }
        Ok(())
    }
}

#[derive(Clone, Encode, Decode, Debug)]
#[trait_decode(trait = GammaLoopContext)]
pub struct GenericEvaluator {
    pub exprs: Option<Vec<Atom>>,
    pub fn_map_entries: Vec<FnMapEntry>,
    pub exprs_len: usize,
    pub rational: Option<ExpressionEvaluator<symbolica::domains::float::Complex<Rational>>>,
    pub f64_compiled: Option<CompiledComplexEvaluatorSpenso>,
    pub f64_eager: ExpressionEvaluator<Complex<F<f64>>>,
    pub f128: ExpressionEvaluator<Complex<F<f128>>>,
    pub dual_shape: Option<Vec<Vec<usize>>>,
    pub arb: ExpressionEvaluator<Complex<F<ArbPrec>>>,
}

impl GenericEvaluator {
    pub(crate) fn compute_out_size(&self) -> usize {
        let number_type_size = if let Some(dual_shape) = &self.dual_shape {
            dual_shape.iter().map(|vec| vec.len()).sum()
        } else {
            1
        };

        number_type_size * self.exprs_len
    }

    pub(crate) fn compile(
        &mut self,
        cpp_path: impl AsRef<Path>,
        function_name: impl AsRef<str>,
        lib_path: impl AsRef<Path>,
        settings: ExportSettings,
    ) {
        let compile = self
            .f64_eager
            .export_cpp::<Complex<f64>>(cpp_path.as_ref(), function_name.as_ref(), settings)
            .unwrap()
            .compile(lib_path.as_ref(), CompileOptions::default())
            .unwrap()
            .load()
            .unwrap();

        self.f64_compiled = Some(compile);
    }

    pub(crate) fn new_from_builder<I: IntoIterator<Item = Atom>>(
        atoms: I,
        builder: &ParamBuilder<f64>,
        dual_shape: Option<Vec<Vec<usize>>>,
        optimization_settings: OptimizationSettings,
        store_atom: bool,
    ) -> Result<Self> {
        let params: Vec<Atom> = (&builder.pairs)
            .into_iter()
            .flat_map(|p| p.params.clone())
            .collect();

        Self::new_from_raw_params(
            atoms,
            &params,
            &builder.fn_map,
            builder.reps.clone(),
            optimization_settings,
            dual_shape,
            store_atom,
        )
    }

    pub(crate) fn new_from_raw_params<I: IntoIterator<Item = Atom>>(
        atoms: I,
        params: &[Atom],
        fn_map: &FunctionMap,
        fn_map_entries: Vec<FnMapEntry>,
        optimization_settings: OptimizationSettings,
        dual_shape: Option<Vec<Vec<usize>>>,
        store_atom: bool,
    ) -> Result<Self> {
        let exprs: Vec<Atom> = atoms.into_iter().collect();

        let mut tree: Option<ExpressionEvaluator<SymComplex<Fraction<IntegerRing>>>> = None;
        for n in exprs.iter() {
            let eval = n
                .evaluator(fn_map, params, optimization_settings.clone())
                .map_err(|e| eyre!("Failed to create evaluator for atom: {:120}\n: {}", n, e))?;

            tree = Some(if let Some(mut tree) = tree {
                tree.merge(eval, optimization_settings.cpe_iterations)
                    .map_err(|e| eyre!("Failed to merge evaluators: {}", e))?;
                tree
            } else {
                eval
            });
        }

        let mut tree = tree.ok_or_else(|| eyre!("No expressions to evaluate"))?;

        if let Some(dual_shape) = &dual_shape {
            let dual = HyperDual::<SymComplex<Rational>>::new(dual_shape.clone());
            let dualizer = Dualizer::new(dual, vec![]);
            tree = tree
                .vectorize(&dualizer, ahash::HashMap::default())
                .unwrap();
        }

        let rational = tree.clone();
        let f64_eager = tree
            .clone()
            .map_coeff(&|r| Complex::new(F::from(&r.re), F::from(&r.im)));

        let f128 = tree
            .clone()
            .map_coeff(&|r| Complex::new(F::from(&r.re), F::from(&r.im)));
        let arb: ExpressionEvaluator<Complex<F<ArbPrec>>> =
            tree.map_coeff(&|r| Complex::new(F::from(&r.re), F::from(&r.im)));

        let evaluator = GenericEvaluator {
            exprs_len: exprs.len(),
            fn_map_entries,
            exprs: if store_atom { Some(exprs) } else { None },
            rational: Some(rational),
            f64_compiled: None,
            f64_eager,
            f128,
            dual_shape,
            arb,
        };

        Ok(evaluator)
    }
}

pub enum SliceMut<'a, T: FloatLike> {
    Borrowed(&'a mut [Complex<F<T>>]),
    Owned(Vec<Complex<F<T>>>),
}

pub struct InputParams<'a, T: FloatLike> {
    pub values: SliceMut<'a, T>,
    pub orientations_start: usize,
    pub override_pos: usize,
    pub multiplicative_offset: usize,
}

impl<'a, T: FloatLike> InputParams<'a, T> {
    pub(crate) fn set_orientation_values_impl<A: Clone + Neg<Output = A>, O: GraphOrientation>(
        values: &mut [A],
        one: A,
        zero: A,
        mult_offset: usize,
        start: usize,
        orientation: &O,
    ) {
        let minusone = -(one.clone());
        let mut o_start = start * mult_offset;

        for (_eid, i) in orientation.orientation() {
            // debug!("Setting orientation input for edge {}: {:?}", eid, i);
            match i {
                Orientation::Default => {
                    values[o_start] = one.clone();
                    o_start += mult_offset;
                }
                Orientation::Reversed => {
                    values[o_start] = minusone.clone();
                    o_start += mult_offset;
                }
                Orientation::Undirected => {
                    values[o_start] = zero.clone();
                    o_start += mult_offset;
                }
            }
        }
    }

    pub(crate) fn set_orientation_values<O: GraphOrientation>(&mut self, orientation: &O) {
        let zero: Complex<F<T>> = Complex::new_re(F(T::from_f64(0.)));
        let one = zero.ref_one();
        let mult_offset = self.multiplicative_offset;
        let start = self.orientations_start;
        Self::set_orientation_values_impl(
            self.as_mut(),
            one,
            zero,
            mult_offset,
            start,
            orientation,
        );
    }

    pub(crate) fn override_if(&mut self, over_ride: bool) {
        let zero: Complex<F<T>> = Complex::new_re(F(T::from_f64(0.)));
        let one = zero.ref_one();
        let multiplicative_offset = self.multiplicative_offset;
        let start = self.override_pos;
        set_override_if(
            self.as_mut(),
            one,
            zero,
            over_ride,
            start,
            multiplicative_offset,
        );
    }
    pub fn as_mut(&mut self) -> &mut [Complex<F<T>>] {
        match &mut self.values {
            SliceMut::Borrowed(s) => s,
            SliceMut::Owned(v) => v,
        }
    }

    pub fn as_slice(&self) -> &[Complex<F<T>>] {
        match &self.values {
            SliceMut::Borrowed(s) => s,
            SliceMut::Owned(v) => v,
        }
    }
}

pub trait GenericEvaluatorFloat<T: FloatLike = Self> {
    fn get_evaluator_single(
        generic_evaluator: &mut GenericEvaluator,
    ) -> impl FnMut(&[Complex<F<T>>]) -> Complex<F<T>>;

    #[allow(clippy::type_complexity)]
    fn get_evaluator(
        generic_evaluator: &mut GenericEvaluator,
    ) -> impl FnMut(&[Complex<F<T>>]) -> Vec<Complex<F<T>>>;

    #[allow(clippy::too_many_arguments)]
    fn get_parameters<'a>(
        param_builder: &'a mut ParamBuilder,
        cache: (bool, bool),
        graph: &'a Graph,
        sample: &'a MomentumSample<T>,
        helicities: &[Helicity],
        additional_params: &[F<T>],
        left_threshold_params: Option<&ThresholdParams<T>>,
        right_threshold_params: Option<&ThresholdParams<T>>,
        lu_params: Option<&LUParams<T>>,
    ) -> InputParams<'a, T>;
}

impl GenericEvaluatorFloat for f64 {
    #[inline(always)]
    fn get_evaluator_single(
        generic_evaluator: &mut GenericEvaluator,
    ) -> impl FnMut(&[Complex<F<f64>>]) -> Complex<F<f64>> {
        #[inline(always)]
        |params: &[Complex<F<f64>>]| {
            if let Some(compiled) = &mut generic_evaluator.f64_compiled {
                // info!("USING COMPILED F64 SINGLE");
                let mut out = [Complex::default()];

                unsafe {
                    compiled.evaluate(
                        transmute::<&[Complex<F<f64>>], &[Complex<f64>]>(params),
                        transmute::<&mut [Complex<F<f64>>], &mut [Complex<f64>]>(&mut out),
                    );
                }
                out[0]
            } else {
                // info!("USING EAGER F64 SINGLE");
                generic_evaluator.f64_eager.evaluate_single(params)
            }
        }
    }

    fn get_evaluator(
        generic_evaluator: &mut GenericEvaluator,
    ) -> impl FnMut(&[Complex<F<Self>>]) -> Vec<Complex<F<Self>>> {
        |params: &[Complex<F<f64>>]| {
            let mut out = vec![Complex::default(); generic_evaluator.compute_out_size()];
            if let Some(compiled) = &mut generic_evaluator.f64_compiled {
                // info!("USING COMPILED COMPLEX SINGLE");
                //
                unsafe {
                    compiled.evaluate(
                        transmute::<&[Complex<F<f64>>], &[Complex<f64>]>(params),
                        transmute::<&mut [Complex<F<f64>>], &mut [Complex<f64>]>(&mut out),
                    );
                }
                out
            } else {
                // info!("USING EAGER COMPLEX SINGLE");
                generic_evaluator.f64_eager.evaluate(params, &mut out);
                out
            }
        }
    }

    fn get_parameters<'a>(
        param_builder: &'a mut ParamBuilder,
        cache: (bool, bool),
        graph: &'a Graph,
        sample: &'a MomentumSample<Self>,
        helicities: &[Helicity],
        additional_params: &[F<f64>],
        left_threshold_params: Option<&ThresholdParams<f64>>,
        right_threshold_params: Option<&ThresholdParams<f64>>,
        lu_params: Option<&LUParams<f64>>,
    ) -> InputParams<'a, f64> {
        param_builder.update_emr_and_get_params(
            cache,
            sample,
            graph,
            helicities,
            additional_params,
            left_threshold_params,
            right_threshold_params,
            lu_params,
        )
    }

    // fn get_debug_evaluator(
    //     generic_evaluator: &GenericEvaluatorDebug,
    // ) -> impl Fn(&[Complex<F<Self>>]) -> Complex<F<Self>> {
    //     #[inline(always)]
    //     |params: &[Complex<F<f64>>]| {
    //         // generic_evaluator
    //         //     .builder
    //         //     .borrow_mut()
    //         //     .fill_in_values(Vec::from_iter(params.iter().cloned()));

    //         // let a = generic_evaluator
    //         //     .builder
    //         //     .borrow()
    //         //     .replace(&generic_evaluator.expr);

    //         // debug!("Replaced atom:{:+>}", a);
    //         // generic_evaluator
    //         //     .expr
    //         //     .evaluate(
    //         //         |c| Complex::new_re(F::<f64>::from(c)),
    //         //         const_map,
    //         //         function_map,
    //         //     )
    //         //     .unwrap()

    //         // generic_evaluator.expr.evaluate(coeff_map, const_map, function_map)

    //         if let Some(compiled) = &generic_evaluator.f64_compiled {
    //             let mut out = [Complex::default()];
    //             compiled.borrow_mut().evaluate(params, &mut out);
    //             out[0]
    //         } else {
    //             generic_evaluator
    //                 .f64_eager
    //                 .borrow_mut()
    //                 .evaluate_single(params)
    //         }
    //     }
    // }
}

impl GenericEvaluatorFloat for f128 {
    #[inline(always)]
    fn get_evaluator_single(
        generic_evaluator: &mut GenericEvaluator,
    ) -> impl FnMut(&[Complex<F<f128>>]) -> Complex<F<f128>> {
        // info!("USING COMPLEX F128 SINGLE");
        #[inline(always)]
        |params: &[Complex<F<f128>>]| generic_evaluator.f128.evaluate_single(params)
    }

    fn get_evaluator(
        generic_evaluator: &mut GenericEvaluator,
    ) -> impl FnMut(&[Complex<F<f128>>]) -> Vec<Complex<F<f128>>> {
        |params: &[Complex<F<f128>>]| {
            // info!("USING COMPLEX F128 MULTIPLE");
            let mut out = vec![Complex::default(); generic_evaluator.compute_out_size()];
            generic_evaluator.f128.evaluate(params, &mut out);
            out
        }
    }

    fn get_parameters<'a>(
        param_builder: &'a mut ParamBuilder,
        cache: (bool, bool),
        graph: &'a Graph,
        sample: &'a MomentumSample<Self>,
        helicities: &[Helicity],
        additional_params: &[F<f128>],
        left_threshold_params: Option<&ThresholdParams<f128>>,
        right_threshold_params: Option<&ThresholdParams<f128>>,
        lu_params: Option<&LUParams<f128>>,
    ) -> InputParams<'a, Self> {
        param_builder.update_emr_and_get_params(
            cache,
            sample,
            graph,
            helicities,
            additional_params,
            left_threshold_params,
            right_threshold_params,
            lu_params,
        )
    }
}

impl GenericEvaluatorFloat for ArbPrec {
    #[inline(always)]
    fn get_evaluator_single(
        generic_evaluator: &mut GenericEvaluator,
    ) -> impl FnMut(&[Complex<F<ArbPrec>>]) -> Complex<F<ArbPrec>> {
        #[inline(always)]
        |params: &[Complex<F<ArbPrec>>]| generic_evaluator.arb.evaluate_single(params)
    }

    fn get_evaluator(
        generic_evaluator: &mut GenericEvaluator,
    ) -> impl FnMut(&[Complex<F<ArbPrec>>]) -> Vec<Complex<F<ArbPrec>>> {
        |params: &[Complex<F<ArbPrec>>]| {
            let mut out = vec![Complex::default(); generic_evaluator.compute_out_size()];
            generic_evaluator.arb.evaluate(params, &mut out);
            out
        }
    }

    fn get_parameters<'a>(
        param_builder: &'a mut ParamBuilder,
        cache: (bool, bool),
        graph: &'a Graph,
        sample: &'a MomentumSample<Self>,
        helicities: &[Helicity],
        additional_params: &[F<ArbPrec>],
        left_threshold_params: Option<&ThresholdParams<ArbPrec>>,
        right_threshold_params: Option<&ThresholdParams<ArbPrec>>,
        lu_params: Option<&LUParams<ArbPrec>>,
    ) -> InputParams<'a, Self> {
        param_builder.update_emr_and_get_params(
            cache,
            sample,
            graph,
            helicities,
            additional_params,
            left_threshold_params,
            right_threshold_params,
            lu_params,
        )
    }
}

#[cfg(test)]
mod tests {
    use symbolica::atom::Symbol;

    use crate::utils::{W_, symbolica_ext::CallSymbol};

    use super::*;

    #[test]
    fn test_function_map_summed() {
        fn integer_to_orientation(i: isize) -> Orientation {
            match i {
                1 => Orientation::Default,
                -1 => Orientation::Reversed,
                0 => Orientation::Undirected,
                _ => panic!("Invalid orientation index"),
            }
        }

        let a = [
            EdgeVec::from_iter([1, -1, 1, 1].into_iter().map(integer_to_orientation)),
            EdgeVec::from_iter([1, -1, -1, 1].into_iter().map(integer_to_orientation)),
            EdgeVec::from_iter([-1, -1, 1, 1].into_iter().map(integer_to_orientation)),
        ]
        .iter()
        .map(|a| {
            GS.collect_orientation_if(a.orientation_thetas() * GS.integrand(1, a), true)

            // GS.integrand(a)
        })
        .fold(Atom::Zero, |acc, n| acc + n)
        .replace(
            Symbol::IF.f([GS.override_if, W_.b_, W_.c_])
                + Symbol::IF.f([GS.override_if, W_.d_, W_.e_]),
        )
        .repeat()
        .with(Symbol::IF.f([
            Atom::var(GS.override_if),
            Atom::var(W_.b_) + Atom::var(W_.d_),
            Atom::var(W_.c_) + Atom::var(W_.e_),
        ]));

        println!("{:100}", a);
    }
}
