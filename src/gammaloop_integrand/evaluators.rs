use std::{borrow::Cow, mem::transmute, path::Path};

use bincode_trait_derive::{Decode, Encode};
use bitvec::{order::Lsb0, slice::IterOnes, vec::BitVec};
use linnet::half_edge::involution::{EdgeVec, Orientation};
use spenso::algebra::complex::{symbolica_traits::CompiledComplexEvaluatorSpenso, Complex};
use symbolica::{
    atom::{Atom, AtomCore},
    domains::{float::Complex as SymComplex, rational::Rational},
    evaluate::{
        CompileOptions, CompiledComplexEvaluator, ExportSettings, ExpressionEvaluator, FunctionMap,
        OptimizationSettings,
    },
};
use typed_index_collections::TiVec;

use crate::{
    gammaloop_integrand::param_builder::LUParams,
    graph::Graph,
    momentum::Helicity,
    momentum_sample::MomentumSample,
    utils::{f128, FloatLike, F},
    GammaLoopContext,
};

use super::{
    param_builder::{ThresholdParams, UpdateAndGetParams},
    ParamBuilder,
};

#[derive(Clone, Copy)]
pub enum SingleOrAllOrientations<'a, OID> {
    Single {
        orientation: &'a EdgeVec<Orientation>,
        id: OID,
    },
    All {
        all: &'a TiVec<OID, EdgeVec<Orientation>>,
        filter: &'a BitVec,
    },
}
impl<'a, OID: Copy> SingleOrAllOrientations<'a, OID> {
    pub fn iter(&self) -> SingleOrAllOrientationsIterator<'_, OID> {
        match self {
            SingleOrAllOrientations::All { all, filter } => SingleOrAllOrientationsIterator::All {
                all: *all,
                filter: filter.iter_ones(),
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

#[derive(Clone, Copy)]
pub enum SingleOrAllOrientationsIterator<'a, OID> {
    Single {
        orientation: Option<&'a EdgeVec<Orientation>>,
        id: OID,
    },
    All {
        all: &'a TiVec<OID, EdgeVec<Orientation>>,
        filter: IterOnes<'a, usize, Lsb0>,
    },
}

impl<'a, OID: From<usize> + Copy> Iterator for SingleOrAllOrientationsIterator<'a, OID>
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
                let a = OID::from(filter.next()?);
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

#[derive(Clone, Encode, Decode, Debug)]
#[trait_decode(trait = GammaLoopContext)]
pub struct GenericEvaluator {
    pub exprs: Vec<Atom>,
    pub rational: Option<ExpressionEvaluator<symbolica::domains::float::Complex<Rational>>>,
    pub f64_compiled: Option<CompiledComplexEvaluatorSpenso>,
    pub f64_eager: ExpressionEvaluator<Complex<F<f64>>>,
    pub f128: ExpressionEvaluator<Complex<F<f128>>>,
}

impl GenericEvaluator {
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
        optimization_settings: OptimizationSettings,
    ) -> Option<Self> {
        let params: Vec<Atom> = (&builder.pairs)
            .into_iter()
            .flat_map(|p| p.params.clone())
            .collect();

        Self::new_from_raw_params(atoms, &params, &builder.fn_map, optimization_settings)
    }

    pub(crate) fn new_from_raw_params<I: IntoIterator<Item = Atom>>(
        atoms: I,
        params: &[Atom],
        fn_map: &FunctionMap,
        optimization_settings: OptimizationSettings,
    ) -> Option<Self> {
        let exprs: Vec<Atom> = atoms.into_iter().collect();
        let tree = exprs
            .iter()
            .map(|n| {
                n.evaluator(fn_map, params, optimization_settings.clone())
                    .unwrap()
            })
            .reduce(|mut acc, n| {
                acc.merge(n, optimization_settings.cpe_iterations).unwrap();
                acc
            })?;

        let rational = tree.clone();
        let f64_eager = tree
            .clone()
            .map_coeff(&|r| Complex::new(F::from(&r.re), F::from(&r.im)));
        let f128 = tree.map_coeff(&|r| Complex::new(F::from(&r.re), F::from(&r.im)));

        let evaluator = GenericEvaluator {
            exprs,
            rational: Some(rational),
            f64_compiled: None,
            f64_eager,
            f128,
        };

        Some(evaluator)
    }
}

pub trait GenericEvaluatorFloat<T: FloatLike = Self> {
    fn get_evaluator_single(
        generic_evaluator: &mut GenericEvaluator,
    ) -> impl FnMut(&[Complex<F<T>>]) -> Complex<F<T>>;

    fn get_evaluator(
        generic_evaluator: &mut GenericEvaluator,
    ) -> impl FnMut(&[Complex<F<T>>]) -> Vec<Complex<F<T>>>;

    fn get_parameters<'a>(
        param_builder: &'a mut ParamBuilder,
        cache: bool,
        graph: &'a Graph,
        sample: &'a MomentumSample<T>,
        helicities: &[Helicity],
        threshold_params: Option<&ThresholdParams<T>>,
        lu_params: Option<&LUParams<T>>,
    ) -> Cow<'a, Vec<Complex<F<T>>>>;
}

impl GenericEvaluatorFloat for f64 {
    #[inline(always)]
    fn get_evaluator_single(
        generic_evaluator: &mut GenericEvaluator,
    ) -> impl FnMut(&[Complex<F<f64>>]) -> Complex<F<f64>> {
        #[inline(always)]
        |params: &[Complex<F<f64>>]| {
            if let Some(compiled) = &mut generic_evaluator.f64_compiled {
                // status_info!("USING COMPILED F64 SINGLE");
                let mut out = [Complex::default()];

                unsafe {
                    compiled.evaluate(
                        transmute::<&[Complex<F<f64>>], &[Complex<f64>]>(params),
                        transmute::<&mut [Complex<F<f64>>], &mut [Complex<f64>]>(&mut out),
                    );
                }
                out[0]
            } else {
                // status_info!("USING EAGER F64 SINGLE");
                generic_evaluator.f64_eager.evaluate_single(params)
            }
        }
    }

    fn get_evaluator(
        generic_evaluator: &mut GenericEvaluator,
    ) -> impl FnMut(&[Complex<F<Self>>]) -> Vec<Complex<F<Self>>> {
        |params: &[Complex<F<f64>>]| {
            let mut out = vec![Complex::default(); generic_evaluator.exprs.len()];
            if let Some(compiled) = &mut generic_evaluator.f64_compiled {
                // status_info!("USING COMPILED COMPLEX SINGLE");
                //
                unsafe {
                    compiled.evaluate(
                        transmute::<&[Complex<F<f64>>], &[Complex<f64>]>(params),
                        transmute::<&mut [Complex<F<f64>>], &mut [Complex<f64>]>(&mut out),
                    );
                }
                out
            } else {
                // status_info!("USING EAGER COMPLEX SINGLE");
                generic_evaluator.f64_eager.evaluate(params, &mut out);
                out
            }
        }
    }

    fn get_parameters<'a>(
        param_builder: &'a mut ParamBuilder,
        cache: bool,
        graph: &'a Graph,
        sample: &'a MomentumSample<Self>,
        helicities: &[Helicity],
        threshold_params: Option<&ThresholdParams<f64>>,
        lu_params: Option<&LUParams<f64>>,
    ) -> Cow<'a, Vec<Complex<F<Self>>>> {
        let params = param_builder.update_emr_and_get_params(
            cache,
            sample,
            graph,
            helicities,
            threshold_params,
            lu_params,
        );

        params
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
        // status_info!("USING COMPLEX F128 SINGLE");
        #[inline(always)]
        |params: &[Complex<F<f128>>]| generic_evaluator.f128.evaluate_single(params)
    }

    fn get_evaluator(
        generic_evaluator: &mut GenericEvaluator,
    ) -> impl FnMut(&[Complex<F<f128>>]) -> Vec<Complex<F<f128>>> {
        |params: &[Complex<F<f128>>]| {
            // status_info!("USING COMPLEX F128 MULTIPLE");
            let mut out = vec![Complex::default(); generic_evaluator.exprs.len()];
            generic_evaluator.f128.evaluate(params, &mut out);
            out
        }
    }

    fn get_parameters<'a>(
        param_builder: &'a mut ParamBuilder,
        cache: bool,
        graph: &'a Graph,
        sample: &'a MomentumSample<Self>,
        helicities: &[Helicity],
        threshold_params: Option<&ThresholdParams<f128>>,
        lu_params: Option<&LUParams<f128>>,
    ) -> Cow<'a, Vec<Complex<F<Self>>>> {
        param_builder.update_emr_and_get_params(
            cache,
            sample,
            graph,
            helicities,
            threshold_params,
            lu_params,
        )
    }
}
