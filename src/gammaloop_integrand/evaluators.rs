use std::{borrow::Cow, cell::RefCell, path::Path};

use bincode_trait_derive::{Decode, Encode};
use bitvec::{order::Lsb0, slice::IterOnes, vec::BitVec};
use linnet::half_edge::involution::{EdgeVec, Orientation};
use spenso::{
    algebra::complex::Complex,
    tensors::parametric::{SerializableCompiledCode, SerializableCompiledEvaluator},
};
use symbolica::{
    atom::{Atom, AtomCore},
    domains::rational::Rational,
    evaluate::{CompileOptions, ExpressionEvaluator, InlineASM, OptimizationSettings},
};
use typed_index_collections::TiVec;

use crate::{
    graph::Graph,
    momentum::Helicity,
    momentum_sample::MomentumSample,
    status_info,
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

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct GenericEvaluator {
    pub exprs: Vec<Atom>,
    pub rational:
        RefCell<Option<ExpressionEvaluator<symbolica::domains::float::Complex<Rational>>>>,
    pub f64_compiled: Option<RefCell<SerializableCompiledEvaluator>>,
    pub f64_eager: RefCell<ExpressionEvaluator<Complex<F<f64>>>>,
    pub f128: RefCell<ExpressionEvaluator<Complex<F<f128>>>>,
}

pub trait GenericEvaluate<T> {
    fn evaluate(&self, params: &[T], out: &mut [T]);
    fn evaluate_single(&self, params: &[T]) -> T;
}

// impl<T: Real + Default> GenericEvaluate<T> for GenericEvaluator
// where
//     T: for<'a> From<&'a symbolica::domains::float::Complex<Rational>>,
//     symbolica::domains::float::Complex<Rational>: for<'a> From<&'a T>,
// {
//     fn evaluate(&self, params: &[T], out: &mut [T]) {
//         let rational = self.rational.take().unwrap();
//         let mut t_eval = rational.map_coeff(&|t| T::from(t));
//         t_eval.evaluate(params, out);
//         self.rational.replace(Some(t_eval.map_coeff(&|t| t.into())));
//     }

//     fn evaluate_single(&self, params: &[T]) -> T {
//         let rational = self.rational.take().unwrap();
//         let mut t_eval = rational.map_coeff(&|t| T::from(t));
//         let out = t_eval.evaluate_single(params);
//         self.rational.replace(Some(t_eval.map_coeff(&|t| t.into())));
//         out
//     }
// }

impl GenericEvaluate<Complex<F<f64>>> for GenericEvaluator {
    fn evaluate(&self, params: &[Complex<F<f64>>], out: &mut [Complex<F<f64>>]) {
        if let Some(f64_compiled) = &self.f64_compiled {
            // status_info!("USING COMPLEX COMPILED MULTIPLE");
            f64_compiled.borrow_mut().evaluate(params, out);
        } else {
            // status_info!("USING COMPLEX COMPILED SINGLE");
            self.f64_eager.borrow_mut().evaluate(params, out);
        }
    }

    fn evaluate_single(&self, params: &[Complex<F<f64>>]) -> Complex<F<f64>> {
        if let Some(f64_compiled) = &self.f64_compiled {
            // status_info!("USING COMPLEX COMPILED SINGLE");
            let mut out = [Complex::default()];
            f64_compiled.borrow_mut().evaluate(params, &mut out);
            out[0]
        } else {
            // status_info!("USING COMPLEX EAGER SINGLE");
            self.f64_eager.borrow_mut().evaluate_single(params)
        }
    }
}

impl GenericEvaluate<Complex<F<f128>>> for GenericEvaluator {
    fn evaluate(&self, params: &[Complex<F<f128>>], out: &mut [Complex<F<f128>>]) {
        // status_info!("USING COMPLEX F128 MULTIPLE");
        self.f128.borrow_mut().evaluate(params, out);
    }

    fn evaluate_single(&self, params: &[Complex<F<f128>>]) -> Complex<F<f128>> {
        // status_info!("USING COMPLEX F128 SINGLE");
        self.f128.borrow_mut().evaluate_single(params)
    }
}

impl GenericEvaluator {
    pub(crate) fn compile(
        &mut self,
        cpp_path: impl AsRef<Path>,
        function_name: impl AsRef<str>,
        lib_path: impl AsRef<Path>,
        inline_asm: InlineASM,
    ) {
        let compile = self
            .f64_eager
            .borrow()
            .export_cpp(
                &cpp_path.as_ref().to_string_lossy(),
                function_name.as_ref(),
                true,
                inline_asm,
            )
            .unwrap()
            .compile(
                &lib_path.as_ref().to_string_lossy(),
                CompileOptions::default(),
            )
            .unwrap()
            .load()
            .unwrap();

        self.f64_compiled = Some(RefCell::new(SerializableCompiledEvaluator {
            evaluator: compile,
            compiled_code: SerializableCompiledCode {
                library_filename: lib_path.as_ref().to_path_buf(),
                function_name: function_name.as_ref().to_string(),
            },
        }));
    }

    pub(crate) fn new_from_builder<I: IntoIterator<Item = Atom>>(
        atoms: I,
        builder: &ParamBuilder<f64>,
        optimization_settings: OptimizationSettings,
    ) -> Option<Self> {
        let params: Vec<Atom> = builder.into_iter().flat_map(|p| p.params.clone()).collect();

        let exprs: Vec<Atom> = atoms.into_iter().collect();
        let tree = exprs
            .iter()
            .map(|n| {
                n.evaluator(&builder.fn_map, &params, optimization_settings.clone())
                    .unwrap()
            })
            .reduce(|mut acc, n| {
                acc.merge(n, optimization_settings.cpe_iterations);
                acc
            })?;

        let rational = tree.clone();
        let f64_eager = tree
            .clone()
            .map_coeff(&|r| Complex::new(F::from(&r.re), F::from(&r.im)));
        let f128 = tree.map_coeff(&|r| Complex::new(F::from(&r.re), F::from(&r.im)));

        let evaluator = GenericEvaluator {
            exprs,
            rational: RefCell::new(Some(rational)),
            f64_compiled: None,
            f64_eager: RefCell::new(f64_eager),
            f128: RefCell::new(f128),
        };

        Some(evaluator)
    }
}

pub trait GenericEvaluatorFloat<T: FloatLike = Self> {
    fn get_evaluator_single(
        generic_evaluator: &GenericEvaluator,
    ) -> impl Fn(&[Complex<F<T>>]) -> Complex<F<T>>;

    fn get_evaluator(
        generic_evaluator: &GenericEvaluator,
    ) -> impl Fn(&[Complex<F<T>>]) -> Vec<Complex<F<T>>>;

    fn get_parameters<'a>(
        param_builder: &'a mut ParamBuilder,
        graph: &'a Graph,
        sample: &'a MomentumSample<T>,
        helicities: &[Helicity],
        threshold_params: Option<&ThresholdParams<T>>,
    ) -> Cow<'a, Vec<Complex<F<T>>>>;
}

impl GenericEvaluatorFloat for f64 {
    #[inline(always)]
    fn get_evaluator_single(
        generic_evaluator: &GenericEvaluator,
    ) -> impl Fn(&[Complex<F<f64>>]) -> Complex<F<f64>> {
        #[inline(always)]
        |params: &[Complex<F<f64>>]| {
            if let Some(compiled) = &generic_evaluator.f64_compiled {
                // status_info!("USING COMPILED F64 SINGLE");
                let mut out = [Complex::default()];
                compiled.borrow_mut().evaluate(params, &mut out);
                out[0]
            } else {
                // status_info!("USING EAGER F64 SINGLE");
                generic_evaluator
                    .f64_eager
                    .borrow_mut()
                    .evaluate_single(params)
            }
        }
    }

    fn get_evaluator(
        generic_evaluator: &GenericEvaluator,
    ) -> impl Fn(&[Complex<F<Self>>]) -> Vec<Complex<F<Self>>> {
        |params: &[Complex<F<f64>>]| {
            let mut out = vec![Complex::default(); generic_evaluator.exprs.len()];
            if let Some(compiled) = &generic_evaluator.f64_compiled {
                // status_info!("USING COMPILED COMPLEX SINGLE");
                compiled.borrow_mut().evaluate(params, &mut out);
                out
            } else {
                // status_info!("USING EAGER COMPLEX SINGLE");
                generic_evaluator
                    .f64_eager
                    .borrow_mut()
                    .evaluate(params, &mut out);
                out
            }
        }
    }

    fn get_parameters<'a>(
        param_builder: &'a mut ParamBuilder,
        graph: &'a Graph,
        sample: &'a MomentumSample<Self>,
        helicities: &[Helicity],
        threshold_params: Option<&ThresholdParams<f64>>,
    ) -> Cow<'a, Vec<Complex<F<Self>>>> {
        let params =
            param_builder.update_emr_and_get_params(sample, graph, helicities, threshold_params);

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
        generic_evaluator: &GenericEvaluator,
    ) -> impl Fn(&[Complex<F<f128>>]) -> Complex<F<f128>> {
        // status_info!("USING COMPLEX F128 SINGLE");
        #[inline(always)]
        |params: &[Complex<F<f128>>]| generic_evaluator.f128.borrow_mut().evaluate_single(params)
    }

    fn get_evaluator(
        generic_evaluator: &GenericEvaluator,
    ) -> impl Fn(&[Complex<F<f128>>]) -> Vec<Complex<F<f128>>> {
        |params: &[Complex<F<f128>>]| {
            // status_info!("USING COMPLEX F128 MULTIPLE");
            let mut out = vec![Complex::default(); generic_evaluator.exprs.len()];
            generic_evaluator
                .f128
                .borrow_mut()
                .evaluate(params, &mut out);
            out
        }
    }

    fn get_parameters<'a>(
        param_builder: &'a mut ParamBuilder,
        graph: &'a Graph,
        sample: &'a MomentumSample<Self>,
        helicities: &[Helicity],
        threshold_params: Option<&ThresholdParams<f128>>,
    ) -> Cow<'a, Vec<Complex<F<Self>>>> {
        param_builder.update_emr_and_get_params(sample, graph, helicities, threshold_params)
    }
}
