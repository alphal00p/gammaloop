use std::{path::Path, sync::Arc, time::Instant};

use bincode_trait_derive::{Decode, Encode};
use color_eyre::eyre::{Result, ensure, eyre};
use linnet::half_edge::{
    involution::{EdgeVec, Orientation},
    typed_vec::IndexLike,
};
use spenso::algebra::complex::Complex;
use symbolica::atom::Atom;

use crate::{
    GammaLoopContext,
    cff::expression::OrientationID,
    graph::Graph,
    integrands::{
        evaluation::EvaluationMetaData,
        process::{
            ParamBuilder,
            evaluators::{
                EvaluatorStack, GenericEvaluator, SingleOrAllOrientations, evaluate_evaluator,
            },
            param_builder::FnMapEntry,
        },
    },
    momentum::sample::{LoopMomenta, MomentumSample},
    processes::{EvaluatorBuildTimings, EvaluatorSettings},
    settings::{RuntimeSettings, global::FrozenCompilationMode},
    utils::{F, FloatLike, h_dual, hyperdual_utils::DualOrNot},
    uv::uv_graph::UVE,
};

use super::{FermiSurfaceProduct, FermiSurfaceSector};

/// A smooth compiled coefficient together with its localized Fermi support.
/// Shell parameters use the same parameter layout and native precision as the
/// coefficient; no model expression is evaluated through an f64 intermediate.
#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub(crate) struct FermiSurfaceEvaluator {
    product: FermiSurfaceProduct,
    shell_parameters: GenericEvaluator,
    coefficient: EvaluatorStack,
}

impl FermiSurfaceEvaluator {
    pub(crate) fn new(
        graph: &Graph,
        sector: FermiSurfaceSector,
        param_builder: &ParamBuilder,
        numerator_definitions: &[Arc<FnMapEntry>],
        orientation_catalog: Option<(&[EdgeVec<Orientation>], &[OrientationID])>,
        settings: &EvaluatorSettings,
    ) -> Result<(Self, EvaluatorBuildTimings)> {
        ensure!(
            sector.orientations.len() == sector.product.factors().len()
                && sector
                    .orientations
                    .iter()
                    .all(|sign| matches!(sign, -1 | 1)),
            "Fermi evaluation requires one fixed energy orientation per surface"
        );
        let (coefficient, mut timings) = EvaluatorStack::from_integrand_with_timings(
            &sector.coefficient,
            param_builder,
            numerator_definitions,
            orientation_catalog,
            Some(sector.product.shape.clone()),
            settings,
        )?;
        let started = Instant::now();
        let shell_parameters = GenericEvaluator::new_from_builder(
            sector
                .product
                .factors()
                .iter()
                .zip(&sector.orientations)
                .flat_map(|(factor, orientation)| {
                    let edge = &graph[factor.edge_id];
                    [
                        edge.mass_atom(),
                        Atom::num(*orientation)
                            * edge
                                .chemical_potential_atom()
                                .expect("validated Fermi edge has a chemical potential"),
                    ]
                }),
            param_builder,
            None,
            settings.optimization_settings(),
            settings,
        )?;
        timings.symbolica_time += started.elapsed();
        Ok((
            Self {
                product: sector.product,
                shell_parameters,
                coefficient,
            },
            timings,
        ))
    }

    pub(crate) fn dual_len(&self) -> usize {
        self.product.shape.len()
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn evaluate<T: FloatLike, OID: IndexLike>(
        &mut self,
        sample: &MomentumSample<T>,
        graph: &Graph,
        settings: &RuntimeSettings,
        param_builder: &mut ParamBuilder,
        orientations: SingleOrAllOrientations<'_, OID>,
        metadata: &mut EvaluationMetaData,
    ) -> Result<Complex<F<T>>>
    where
        usize: From<OID>,
    {
        ensure!(
            sample.sample.dual_loop_moms.is_none(),
            "Fermi evaluation requires an unlocalized sampling point"
        );
        let helicities = settings.kinematics.externals.get_helicities();
        let additional_params = settings.additional_params();
        let cache = (settings.general.enable_cache, settings.general.debug_cache);
        let input = T::get_parameters(
            param_builder,
            cache,
            graph,
            sample,
            helicities,
            &additional_params,
            None,
            None,
            None,
        );
        let parameters = evaluate_evaluator(&mut self.shell_parameters, input.as_slice(), metadata);
        let mut masses = Vec::with_capacity(self.product.factors().len());
        let mut chemical_potentials = Vec::with_capacity(self.product.factors().len());
        for pair in parameters.as_chunks::<2>().0 {
            let (DualOrNot::NonDual(mass), DualOrNot::NonDual(nu)) = (&pair[0], &pair[1]) else {
                return Err(eyre!("Fermi shell parameters must be scalar"));
            };
            ensure!(
                mass.im == mass.im.zero() && nu.im == nu.im.zero(),
                "Fermi evaluation requires real masses and chemical potentials"
            );
            masses.push(mass.re.clone());
            chemical_potentials.push(nu.re.clone());
        }
        let adapted = sample.lmb_transform(&graph.loop_momentum_basis, self.product.lmb());
        self.product.localize(
            adapted.loop_moms(),
            &masses,
            &chemical_potentials,
            |scale| Ok(h_dual(scale, None, None, &settings.lu_h_function)),
            |momenta| {
                let mut projected = adapted.clone();
                projected.sample.loop_moms = LoopMomenta::from_iter(
                    momenta
                        .iter()
                        .map(|momentum| momentum.map_ref(&|value| value.values[0].clone())),
                );
                projected.sample.dual_loop_moms = Some(momenta.clone());
                let projected =
                    projected.lmb_transform(self.product.lmb(), &graph.loop_momentum_basis);
                let input = T::get_parameters(
                    param_builder,
                    cache,
                    graph,
                    &projected,
                    helicities,
                    &additional_params,
                    None,
                    None,
                    None,
                );
                match self
                    .coefficient
                    .evaluate(input, orientations, settings, metadata)?
                    .pop()
                {
                    Some(DualOrNot::Dual(value)) => Ok(value),
                    _ => Err(eyre!(
                        "Fermi coefficient evaluator must return a derivative jet"
                    )),
                }
            },
        )
    }

    pub(crate) fn compile(
        &mut self,
        name: impl AsRef<str>,
        path: impl AsRef<Path>,
        frozen_mode: &FrozenCompilationMode,
    ) -> Result<()> {
        self.coefficient.compile(&name, &path, frozen_mode)?;
        let name = format!("{}_shell_parameters", name.as_ref());
        self.shell_parameters.compile_external(
            path.as_ref().join(&name).with_extension("cpp"),
            &name,
            path.as_ref().join(&name).with_extension("so"),
            frozen_mode,
        )
    }

    pub(crate) fn for_each_generic_evaluator_mut(
        &mut self,
        mut f: impl FnMut(&mut GenericEvaluator) -> Result<()>,
    ) -> Result<()> {
        f(&mut self.shell_parameters)?;
        self.coefficient.for_each_generic_evaluator_mut(f)
    }

    pub(crate) fn generic_evaluator_count(&self) -> usize {
        1 + self.coefficient.generic_evaluator_count()
    }
}

#[cfg(test)]
mod tests {
    use linnet::half_edge::{involution::EdgeIndex, subgraph::subset::SubSet};
    use symbolica::atom::AtomCore;
    use three_dimensional_reps::ThermalDistributionFactor;
    use typed_index_collections::TiVec;

    use super::*;
    use crate::{
        dot,
        graph::parse::IntoGraph,
        initialisation::test_initialise,
        integrands::process::sampling::maps::SamplingEvaluationError,
        momentum::{ThreeMomentum, sample::BareMomentumSample},
        settings::{
            global::{CompilationOptimizationLevel, CompilationOptionsSnapshot},
            runtime::HFunction,
        },
        utils::{GS, h},
    };

    #[test]
    fn compiled_fermi_sector_uses_fixed_orientation_and_projected_energy_jets() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph fermi_runtime_cycle {
            node [num=1]
            edge [num=1 particle="b"]
            A -> B [id=0 lmb_id=0]
            B -> A [id=1]
        })?;
        let sample = MomentumSample {
            sample: BareMomentumSample {
                loop_moms: LoopMomenta(vec![ThreeMomentum::new(F(2.0), F(0.0), F(0.0))]),
                dual_loop_moms: None,
                loop_mom_cache_id: 0,
                loop_mom_base_cache_id: 0,
                external_moms: Vec::new().into(),
                external_mom_cache_id: 0,
                external_mom_base_cache_id: 0,
                jacobian: F(1.0),
                orientation: None,
                parameterization_branch: None,
            },
        };
        let catalogue = TiVec::<OrientationID, EdgeVec<Orientation>>::new();
        let filter = SubSet::full(0);
        let orientations = SingleOrAllOrientations::All {
            all: &catalogue,
            filter: &filter,
        };
        let mut runtime = RuntimeSettings::default();
        runtime.lu_h_function.function = HFunction::Exponential;
        runtime.lu_h_function.sigma = 1.0;
        for orientation in [-1, 1] {
            let mut builder = graph.param_builder.clone();
            let sector = FermiSurfaceSector {
                product: FermiSurfaceProduct::new(
                    &graph,
                    &[ThermalDistributionFactor {
                        edge_id: EdgeIndex(0),
                        sign: -1,
                        derivative_order: 2,
                    }],
                )?,
                orientations: vec![orientation],
                coefficient: GS.ose(EdgeIndex(1)).pow(2),
            };
            let (mut evaluator, _) = FermiSurfaceEvaluator::new(
                &graph,
                sector,
                &builder,
                &[],
                None,
                &EvaluatorSettings::default(),
            )?;
            builder.initialize_duals(evaluator.dual_len());
            for compile in [false, true] {
                if compile {
                    evaluator.for_each_generic_evaluator_mut(|evaluator| {
                        evaluator.activate_symjit(&CompilationOptionsSnapshot {
                            optimization_level: CompilationOptimizationLevel::O0,
                            ..Default::default()
                        })
                    })?;
                }
                for mu in [-5.0, 5.0] {
                    for (atom, value) in [
                        (graph[EdgeIndex(0)].mass_atom(), 3.0),
                        (graph[EdgeIndex(0)].chemical_potential_atom().unwrap(), mu),
                    ] {
                        let position = builder
                            .pairs
                            .model_parameters
                            .params
                            .iter()
                            .position(|parameter| *parameter == atom)
                            .unwrap()
                            + builder.pairs.model_parameters.value_range.start;
                        for (index, values) in builder.values.iter_mut().enumerate() {
                            values[position * (index + 1)] = Complex::new_re(F(value));
                        }
                    }
                    let actual = evaluator.evaluate(
                        &sample,
                        &graph,
                        &runtime,
                        &mut builder,
                        orientations,
                        &mut EvaluationMetaData::new_empty(),
                    )?;
                    // -d_E [h(sqrt(E²-9)/2) E³(E²-9)/16] at E=5.
                    let expected = if orientation as f64 * mu > 0.0 {
                        1275.0 / 8.0 * h(&F(2.0), None, None, &runtime.lu_h_function).0
                    } else {
                        0.0
                    };
                    assert!(
                        (actual.re.0 - expected).abs() < 2e-13,
                        "orientation={orientation}, mu={mu}, compiled={compile}: {actual:?} != {expected}"
                    );
                    assert_eq!(actual.im, F(0.0));
                }
            }
        }
        Ok(())
    }

    #[test]
    fn fermi_kernel_preserves_callback_errors_and_marks_numerical_failure() -> Result<()> {
        let graph = super::super::routing::test_graph()?;
        let product = FermiSurfaceProduct::new(
            &graph,
            &[ThermalDistributionFactor {
                edge_id: EdgeIndex(1),
                sign: 1,
                derivative_order: 1,
            }],
        )?;
        for (momentum, numerical) in [(0.0, true), (1.0, false)] {
            let momenta = LoopMomenta(vec![ThreeMomentum::new(F(momentum), F(0.0), F(0.0)); 2]);
            let error = product
                .localize(
                    &momenta,
                    &[F(0.0)],
                    &[F(1.0)],
                    |_| Err(eyre!("original user callback error")),
                    |_| unreachable!(),
                )
                .unwrap_err();
            assert_eq!(
                error.downcast_ref::<SamplingEvaluationError>().is_some(),
                numerical
            );
            if !numerical {
                assert_eq!(error.to_string(), "original user callback error");
            }
        }
        Ok(())
    }
}
