use std::{collections::BTreeMap, path::Path, slice};

use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;
use eyre::{Context, eyre};
use linnet::half_edge::involution::{EdgeVec, Orientation};
use spenso::algebra::{algebraic_traits::IsZero, complex::Complex};
use tracing::{debug, instrument, warn};
use typed_index_collections::{TiVec, ti_vec};

use crate::{
    GammaLoopContext,
    cff::{
        CutCFFIndex,
        esurface::{
            Esurface, EsurfaceCollection, EsurfaceID, ExistingEsurfaceId, ExistingEsurfaces,
            ExistingThresholds, GroupEsurfaceId, RaisedEsurfaceData, RaisedEsurfaceId,
            esurface_value_is_strictly_inside,
        },
        expression::OrientationID,
    },
    graph::{FeynmanGraph, Graph, GraphGroupPosition, LmbIndex, LoopMomentumBasis},
    integrands::{
        evaluation::EvaluationMetaData,
        process::{
            GenericEvaluator, ParamBuilder, ThresholdParams,
            evaluators::{
                EvaluatorStack, SingleOrAllOrientations, evaluate_evaluator,
                evaluate_evaluator_single,
            },
            threshold_multiplier::{
                ThresholdMultiplierEvaluatorCollection, ThresholdMultiplierInputWorkspace,
            },
        },
    },
    model::Model,
    momentum::{
        Energy, FourMomentum, Rotation, ThreeMomentum,
        sample::{
            BareMomentumSample, ExternalFourMomenta, LoopIndex, LoopMomenta, MomentumSample,
            SubspaceData,
        },
    },
    processes::{
        EvaluatorBuildTimings, ResolvedThresholdCountertermVariant, SingleThresholdPieces,
        ThresholdCountertermComponentKind, ThresholdCountertermMetadataRegistry,
        ThresholdCountertermSide, ThresholdCountertermVariantId,
    },
    settings::{GlobalSettings, RuntimeSettings},
    subtraction::{
        evaluate_integrated_ct_normalisation, evaluate_uv_damper,
        overlap::{OverlapGroup, OverlapStructure},
        overlap_subspace::{
            OverlapInput as SubspaceOverlapInput, OverlapStructure as SubspaceOverlapStructure,
            find_maximal_overlap as find_maximal_subspace_overlap,
        },
    },
    utils::{
        F, FloatLike,
        hyperdual_utils::{
            DualOrNot, extract_t_derivatives, extract_t_derivatives_complex, new_constant,
            shape_from_cut_cff_index, simple_n_deriv_shape,
        },
        newton_solver::{NewtonIterationResult, RadialRootDiagnostics, RadialRootIdentity},
    },
    uv::Integrands,
};
use symbolica::domains::dual::{DualNumberStructure, HyperDual};

const MAX_ITERATIONS: usize = 40;

fn multiply_dual_or_not_complex<T: FloatLike>(
    lhs: DualOrNot<Complex<F<T>>>,
    rhs: &DualOrNot<Complex<F<T>>>,
) -> DualOrNot<Complex<F<T>>> {
    match (lhs, rhs) {
        (DualOrNot::NonDual(lhs), DualOrNot::NonDual(rhs)) => DualOrNot::NonDual(lhs * rhs),
        (DualOrNot::Dual(lhs), DualOrNot::Dual(rhs)) => DualOrNot::Dual(lhs * rhs.clone()),
        (DualOrNot::Dual(lhs), DualOrNot::NonDual(rhs)) => {
            let rhs_dual = new_constant(&lhs, rhs);
            DualOrNot::Dual(lhs * rhs_dual)
        }
        (DualOrNot::NonDual(lhs), DualOrNot::Dual(rhs)) => {
            let lhs_dual = new_constant(rhs, &lhs);
            DualOrNot::Dual(lhs_dual * rhs.clone())
        }
    }
}

fn combine_completed_threshold_pieces<T: FloatLike>(
    local: Complex<F<T>>,
    integrated: Complex<F<T>>,
    local_multiplier: &F<T>,
    integrated_multiplier: &F<T>,
) -> Complex<F<T>> {
    let mut result = Complex::new_re(local.re.zero());
    if !local_multiplier.is_zero() {
        result += local * Complex::new_re(local_multiplier.clone());
    }
    if !integrated_multiplier.is_zero() {
        result += integrated * Complex::new_re(integrated_multiplier.clone());
    }
    result
}

struct EvaluatedAmplitudeThreshold<T: FloatLike> {
    weighted: Complex<F<T>>,
    bare_pieces: Option<SingleThresholdPieces<Complex<F<T>>>>,
}

impl<T: FloatLike> EvaluatedAmplitudeThreshold<T> {
    fn new(zero: F<T>, record_components: bool) -> Self {
        Self {
            weighted: Complex::new_re(zero.clone()),
            bare_pieces: record_components.then(|| SingleThresholdPieces {
                local: Complex::new_re(zero.clone()),
                integrated: Complex::new_re(zero),
            }),
        }
    }

    fn add_completed_pieces(
        &mut self,
        local: Complex<F<T>>,
        integrated: Complex<F<T>>,
        local_multiplier: &F<T>,
        integrated_multiplier: &F<T>,
    ) {
        if let Some(bare_pieces) = &mut self.bare_pieces {
            if !local_multiplier.is_zero() {
                bare_pieces.local += local.clone();
            }
            if !integrated_multiplier.is_zero() {
                bare_pieces.integrated += integrated.clone();
            }
        }
        self.weighted += combine_completed_threshold_pieces(
            local,
            integrated,
            local_multiplier,
            integrated_multiplier,
        );
    }

    fn finish(mut self, local_multiplier: &F<T>, integrated_multiplier: &F<T>) -> Self {
        if let Some(bare_pieces) = &self.bare_pieces {
            self.weighted = combine_completed_threshold_pieces(
                bare_pieces.local.clone(),
                bare_pieces.integrated.clone(),
                local_multiplier,
                integrated_multiplier,
            );
        }
        self
    }
}

fn append_amplitude_component_evaluations<T: FloatLike>(
    components: &mut Option<Vec<AmplitudeCountertermComponentEvaluation<T>>>,
    evaluation: &EvaluatedAmplitudeThreshold<T>,
    variant_id: ThresholdCountertermVariantId,
    esurface_id: RaisedEsurfaceId,
    overlap_group: usize,
    local_multiplier: &F<T>,
    integrated_multiplier: &F<T>,
) {
    let Some(components) = components else {
        return;
    };
    let bare_pieces = evaluation
        .bare_pieces
        .as_ref()
        .expect("recorded amplitude components require completed bare pieces");
    for (kind, multiplier, bare_piece) in [
        (
            ThresholdCountertermComponentKind::Local,
            local_multiplier,
            &bare_pieces.local,
        ),
        (
            ThresholdCountertermComponentKind::Integrated,
            integrated_multiplier,
            &bare_pieces.integrated,
        ),
    ] {
        let evaluation_skipped = multiplier.is_zero();
        let (bare, weighted) = if evaluation_skipped {
            (None, Complex::new_re(multiplier.zero()))
        } else {
            (
                Some(bare_piece.clone()),
                bare_piece.clone() * Complex::new_re(multiplier.clone()),
            )
        };
        components.push(AmplitudeCountertermComponentEvaluation {
            variant_id,
            kind,
            esurface_id,
            overlap_group,
            multiplier_value: multiplier.clone(),
            bare,
            weighted,
            evaluation_skipped,
        });
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct AmplitudeCountertermData {
    /// True retains the historical grouped raised-surface runtime without generalized allocation.
    pub legacy_equivalent: bool,
    pub overlap: OverlapStructure,
    pub evaluators: TiVec<RaisedEsurfaceId, AmplitudeCountertermEvaluator>,
    pub helper_evaluators: Vec<GenericEvaluator>,
    // `generated_mask` tracks whether a threshold slot actually has symbolic content.
    // This is independent of any orientation filtering and reflects the underlying
    // threshold-generation logic itself.
    pub generated_mask: TiVec<RaisedEsurfaceId, bool>,
    // `active_mask` is an additional generation-time gate coming from the selected
    // orientation subset. A slot can be generated in principle but inactive for the
    // current generated evaluator set, in which case we keep its index but compile it
    // to a zero/dummy evaluator and hard-fail if runtime ever reaches it.
    pub active_mask: TiVec<RaisedEsurfaceId, bool>,
    /// Generalized data indexed exactly like the graph's resolved variant registry.
    pub variant_evaluators: TiVec<ThresholdCountertermVariantId, AmplitudeCountertermEvaluator>,
    pub variant_helper_evaluators: TiVec<ThresholdCountertermVariantId, Vec<GenericEvaluator>>,
    pub variant_generated_mask: TiVec<ThresholdCountertermVariantId, bool>,
    pub variant_active_mask: TiVec<ThresholdCountertermVariantId, bool>,
    pub variant_raised_esurfaces: TiVec<ThresholdCountertermVariantId, RaisedEsurfaceId>,
    pub variant_subspaces: TiVec<ThresholdCountertermVariantId, SubspaceData>,
    /// Stable graph-local metadata registry for diagnostics, events, and archives.
    pub variant_metadata: TiVec<ThresholdCountertermVariantId, ResolvedThresholdCountertermVariant>,
    /// Eager-only opaque multiplier evaluators. `None` is the identity/no-allocation path.
    pub threshold_multipliers: Option<ThresholdMultiplierEvaluatorCollection>,
    pub metadata_registry: Option<ThresholdCountertermMetadataRegistry>,
    pub lmbs: TiVec<LmbIndex, LoopMomentumBasis>,
    pub raised_data: RaisedEsurfaceData,
    pub esurface_map: TiVec<GroupEsurfaceId, TiVec<GraphGroupPosition, Option<RaisedEsurfaceId>>>,
    pub local_esurface_exists: TiVec<GroupEsurfaceId, bool>,
    pub own_group_position: GraphGroupPosition,
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct AmplitudeCountertermAtom {
    pub parametric: Integrands,
}

impl AmplitudeCountertermAtom {
    pub(crate) fn is_generated(&self) -> bool {
        self.parametric.iter().next().is_some()
    }

    pub(crate) fn zero_like(&self) -> Self {
        Self {
            parametric: self.parametric.map(|_| symbolica::atom::Atom::Zero),
        }
    }

    #[instrument(skip_all)]
    pub(crate) fn to_evaluator_with_timings(
        &self,
        param_builder: &ParamBuilder,
        orientations: &TiVec<OrientationID, EdgeVec<Orientation>>,
        global_settings: &GlobalSettings,
    ) -> (AmplitudeCountertermEvaluator, EvaluatorBuildTimings) {
        let _progress_guard =
            crate::processes::enter_detailed_progress_span("Building Threshold CT Evaluator");
        let mut evaluator_stacks = BTreeMap::new();
        let mut timings = EvaluatorBuildTimings::default();

        for (index, integrand) in self.parametric.iter() {
            let dual_shape = shape_from_cut_cff_index(index);

            let (evaluator_stack, evaluator_timings) = EvaluatorStack::new_with_timings(
                slice::from_ref(integrand),
                param_builder,
                orientations.as_slice().as_ref(),
                dual_shape,
                &global_settings.generation.evaluator,
            )
            .unwrap();
            timings += evaluator_timings;
            evaluator_stacks.insert(*index, evaluator_stack);
        }

        (AmplitudeCountertermEvaluator { evaluator_stacks }, timings)
    }

    pub(crate) fn new() -> Self {
        Self {
            parametric: std::iter::empty().collect(),
        }
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct AmplitudeCountertermEvaluator {
    pub evaluator_stacks: BTreeMap<CutCFFIndex, EvaluatorStack>,
}

impl AmplitudeCountertermEvaluator {
    pub(crate) fn generic_evaluator_count(&self) -> usize {
        self.evaluator_stacks
            .values()
            .map(EvaluatorStack::generic_evaluator_count)
            .sum()
    }
}

#[derive(Debug, Clone)]
pub struct AmplitudeLocalCountertermEvaluation<T: FloatLike> {
    pub variant_id: Option<ThresholdCountertermVariantId>,
    pub esurface_id: RaisedEsurfaceId,
    pub overlap_group: usize,
    pub value: Complex<F<T>>,
}

#[derive(Debug, Clone)]
pub struct AmplitudeCountertermComponentEvaluation<T: FloatLike> {
    pub variant_id: ThresholdCountertermVariantId,
    pub kind: ThresholdCountertermComponentKind,
    pub esurface_id: RaisedEsurfaceId,
    pub overlap_group: usize,
    pub multiplier_value: F<T>,
    pub bare: Option<Complex<F<T>>>,
    pub weighted: Complex<F<T>>,
    pub evaluation_skipped: bool,
}

#[derive(Debug, Clone)]
pub struct AmplitudeCountertermEvaluation<T: FloatLike> {
    pub total: Complex<F<T>>,
    pub local_counterterms: Vec<AmplitudeLocalCountertermEvaluation<T>>,
    pub components: Option<Vec<AmplitudeCountertermComponentEvaluation<T>>>,
}

impl AmplitudeCountertermData {
    fn radial_root_identity(
        graph_name: &str,
        overlap_group: usize,
        raised_esurface_id: RaisedEsurfaceId,
        rotation: &Rotation,
    ) -> RadialRootIdentity {
        RadialRootIdentity::new(format!(
            "amplitude graph '{graph_name}' overlap group {overlap_group} raised E-surface {} probe rotation {}",
            raised_esurface_id.0, rotation.method,
        ))
    }

    fn variant_radial_root_identity(
        graph_name: &str,
        overlap_group: usize,
        variant_id: ThresholdCountertermVariantId,
        raised_esurface_id: RaisedEsurfaceId,
        rotation: &Rotation,
    ) -> RadialRootIdentity {
        RadialRootIdentity::new(format!(
            "amplitude graph '{graph_name}' overlap group {overlap_group} threshold variant {} raised E-surface {} probe rotation {}",
            variant_id.0, raised_esurface_id.0, rotation.method,
        ))
    }

    pub(crate) fn generic_evaluator_count(&self) -> usize {
        let stack_count = self
            .evaluators
            .iter()
            .map(AmplitudeCountertermEvaluator::generic_evaluator_count)
            .chain(
                self.variant_evaluators
                    .iter()
                    .map(AmplitudeCountertermEvaluator::generic_evaluator_count),
            )
            .sum::<usize>();
        let overlap_count = self
            .overlap
            .overlap_groups
            .iter()
            .map(|group| {
                group
                    .prefactor_evaluator
                    .as_ref()
                    .map(|evaluators| evaluators.len())
                    .unwrap_or(0)
            })
            .sum::<usize>();
        stack_count
            + overlap_count
            + self.helper_evaluators.len()
            + self
                .variant_helper_evaluators
                .iter()
                .map(Vec::len)
                .sum::<usize>()
    }

    pub fn new_empty(own_group_position: GraphGroupPosition) -> Self {
        Self {
            legacy_equivalent: true,
            overlap: OverlapStructure::new_empty(),
            evaluators: TiVec::new(),
            helper_evaluators: vec![],
            generated_mask: TiVec::new(),
            active_mask: TiVec::new(),
            variant_evaluators: TiVec::new(),
            variant_helper_evaluators: TiVec::new(),
            variant_generated_mask: TiVec::new(),
            variant_active_mask: TiVec::new(),
            variant_raised_esurfaces: TiVec::new(),
            variant_subspaces: TiVec::new(),
            variant_metadata: TiVec::new(),
            threshold_multipliers: None,
            metadata_registry: None,
            lmbs: TiVec::new(),
            raised_data: RaisedEsurfaceData {
                raised_groups: TiVec::new(),
                pass_two_evaluator: None,
            },
            esurface_map: TiVec::new(),
            local_esurface_exists: TiVec::new(),
            own_group_position,
        }
    }

    fn ensure_active_raised_esurface(&self, raised_esurface_id: RaisedEsurfaceId) -> Result<()> {
        if self.active_mask[raised_esurface_id] {
            return Ok(());
        }

        Err(color_eyre::eyre::eyre!(
            "Amplitude threshold evaluator {} was reached at runtime even though generation marked it inactive for the selected orientation subset",
            raised_esurface_id.0
        ))
    }

    fn ensure_active_variant(&self, variant_id: ThresholdCountertermVariantId) -> Result<()> {
        if self.variant_active_mask[variant_id] {
            return Ok(());
        }

        Err(eyre!(
            "Amplitude threshold variant {} was reached at runtime even though generation marked it inactive for the selected orientation subset",
            variant_id.0,
        ))
    }

    fn metadata_variant_for_raised_esurface(
        &self,
        raised_esurface_id: RaisedEsurfaceId,
    ) -> Result<ThresholdCountertermVariantId> {
        let registry = self.metadata_registry.as_ref().ok_or_else(|| {
            eyre!("Amplitude threshold metadata was requested without a graph-local registry")
        })?;
        let raised_group = &self.raised_data.raised_groups[raised_esurface_id];
        let mut matching = registry.variants.iter().filter(|variant| {
            variant.side == ThresholdCountertermSide::Amplitude
                && variant.threshold_esurface_ids.len() == raised_group.esurface_ids.len()
                && variant
                    .threshold_esurface_ids
                    .iter()
                    .zip(&raised_group.esurface_ids)
                    .all(|(metadata_id, esurface_id)| *metadata_id == esurface_id.0)
        });
        let variant = matching.next().ok_or_else(|| {
            eyre!(
                "Amplitude threshold registry has no variant for raised E-surface {}",
                raised_esurface_id.0,
            )
        })?;
        if matching.next().is_some() {
            return Err(eyre!(
                "Amplitude threshold registry has multiple variants for legacy raised E-surface {}",
                raised_esurface_id.0,
            ));
        }
        Ok(ThresholdCountertermVariantId(variant.variant_id))
    }

    #[allow(clippy::too_many_arguments)]
    fn evaluate_variant_multiplier<T: FloatLike>(
        threshold_multipliers: &mut Option<ThresholdMultiplierEvaluatorCollection>,
        multiplier_workspace: Option<&mut ThresholdMultiplierInputWorkspace<T>>,
        variant_id: ThresholdCountertermVariantId,
        graph: &Graph,
        edge_masses: &EdgeVec<F<T>>,
        model_values: &[Complex<F<f64>>],
        additional_values: &[F<T>],
        base_sample: &MomentumSample<T>,
        root_sample: &MomentumSample<T>,
        evaluation_metadata: &mut EvaluationMetaData,
        record_primary_timing: bool,
    ) -> Result<(F<T>, F<T>)> {
        let Some(collection) = threshold_multipliers.as_mut() else {
            return Ok((base_sample.one(), base_sample.one()));
        };
        let reference = collection
            .left_variants()
            .iter()
            .find(|reference| reference.variant_id == variant_id)
            .copied()
            .ok_or_else(|| {
                eyre!(
                    "Amplitude threshold variant {} has no multiplier-registry entry",
                    variant_id.0,
                )
            })?;
        let Some(evaluator_id) = reference.evaluator_id else {
            return Ok((base_sample.one(), base_sample.one()));
        };

        let workspace = multiplier_workspace.ok_or_else(|| {
            eyre!(
                "Amplitude threshold variant {} requires a missing multiplier workspace",
                variant_id.0,
            )
        })?;
        let local_values = {
            let layout = collection.layout();
            workspace.bind(
                layout,
                &graph.loop_momentum_basis,
                edge_masses,
                model_values,
                additional_values,
                base_sample,
                root_sample,
            )
            .with_context(|| {
                format!(
                    "Failed to bind local threshold multiplier for amplitude graph '{}' variant {}",
                    graph.name, variant_id.0,
                )
            })?
        };
        let local = collection.evaluators_mut()[evaluator_id.0].evaluate(
            local_values,
            evaluation_metadata,
            record_primary_timing,
        )
        .with_context(|| {
            format!(
                "Failed to evaluate local threshold multiplier for amplitude graph '{}' variant {}",
                graph.name, variant_id.0,
            )
        })?;

        let integrated_values = {
            let layout = collection.layout();
            workspace.bind(
                layout,
                &graph.loop_momentum_basis,
                edge_masses,
                model_values,
                additional_values,
                root_sample,
                root_sample,
            )
            .with_context(|| {
                format!(
                    "Failed to bind integrated threshold multiplier for amplitude graph '{}' variant {}",
                    graph.name, variant_id.0,
                )
            })?
        };
        let integrated = collection.evaluators_mut()[evaluator_id.0].evaluate(
            integrated_values,
            evaluation_metadata,
            record_primary_timing,
        )
        .with_context(|| {
            format!(
                "Failed to evaluate integrated threshold multiplier for amplitude graph '{}' variant {}",
                graph.name, variant_id.0,
            )
        })?;

        Ok((local, integrated))
    }

    pub fn compile(
        &mut self,
        path: impl AsRef<Path>,
        _override_existing: bool,
        frozen_mode: &crate::settings::global::FrozenCompilationMode,
    ) -> Result<()> {
        for (i, e) in self.evaluators.iter_mut_enumerated() {
            for (cff_index, evaluator_stack) in e.evaluator_stacks.iter_mut() {
                let order_index = cff_index.left_threshold_order.unwrap() - 1;
                evaluator_stack.compile(
                    format!("esurface_{}_order_{}", i.0, order_index + 1),
                    path.as_ref(),
                    frozen_mode,
                )?;
            }
        }

        for (variant_id, evaluator) in self.variant_evaluators.iter_mut_enumerated() {
            for (cff_index, evaluator_stack) in &mut evaluator.evaluator_stacks {
                let order_index = cff_index.left_threshold_order.unwrap() - 1;
                evaluator_stack.compile(
                    format!(
                        "threshold_variant_{}_order_{}",
                        variant_id.0,
                        order_index + 1,
                    ),
                    path.as_ref(),
                    frozen_mode,
                )?;
            }
        }

        for (group_index, group) in self.overlap.overlap_groups.iter_mut().enumerate() {
            if let Some(prefactor_evaluators) = group.prefactor_evaluator.as_mut() {
                for (order_index, prefactor_evaluator) in
                    prefactor_evaluators.iter_mut().enumerate()
                {
                    prefactor_evaluator.borrow_mut().compile_external(
                        path.as_ref()
                            .join(format!(
                                "overlap_prefactor_{group_index}_order_{}",
                                order_index + 1
                            ))
                            .with_extension("cpp"),
                        format!("overlap_prefactor_{group_index}_order_{}", order_index + 1),
                        path.as_ref()
                            .join(format!(
                                "overlap_prefactor_{group_index}_order_{}",
                                order_index + 1
                            ))
                            .with_extension("so"),
                        frozen_mode,
                    )?;
                }
            }
        }

        for (order_index, evaluator) in self.helper_evaluators.iter_mut().enumerate() {
            evaluator.compile_external(
                path.as_ref()
                    .join(format!("threshold_helper_{order_index}"))
                    .with_extension("cpp"),
                format!("threshold_helper_{order_index}"),
                path.as_ref()
                    .join(format!("threshold_helper_{order_index}"))
                    .with_extension("so"),
                frozen_mode,
            )?;
        }
        for (variant_id, helpers) in self.variant_helper_evaluators.iter_mut_enumerated() {
            for (order_index, evaluator) in helpers.iter_mut().enumerate() {
                evaluator.compile_external(
                    path.as_ref()
                        .join(format!(
                            "threshold_variant_{}_helper_{order_index}",
                            variant_id.0,
                        ))
                        .with_extension("cpp"),
                    format!("threshold_variant_{}_helper_{order_index}", variant_id.0,),
                    path.as_ref()
                        .join(format!(
                            "threshold_variant_{}_helper_{order_index}",
                            variant_id.0,
                        ))
                        .with_extension("so"),
                    frozen_mode,
                )?;
            }
        }
        Ok(())
    }

    pub(crate) fn for_each_generic_evaluator_mut(
        &mut self,
        mut f: impl FnMut(&mut crate::integrands::process::GenericEvaluator) -> Result<()>,
    ) -> Result<()> {
        for evaluator in self.evaluators.iter_mut() {
            for evaluator_stack in evaluator.evaluator_stacks.values_mut() {
                evaluator_stack.for_each_generic_evaluator_mut(&mut f)?;
            }
        }
        for evaluator in &mut self.variant_evaluators {
            for evaluator_stack in evaluator.evaluator_stacks.values_mut() {
                evaluator_stack.for_each_generic_evaluator_mut(&mut f)?;
            }
        }

        for group in &mut self.overlap.overlap_groups {
            if let Some(prefactor_evaluators) = group.prefactor_evaluator.as_mut() {
                for prefactor_evaluator in prefactor_evaluators.iter_mut() {
                    f(prefactor_evaluator.get_mut())?;
                }
            }
        }

        for evaluator in &mut self.helper_evaluators {
            f(evaluator)?;
        }
        for helpers in &mut self.variant_helper_evaluators {
            for evaluator in helpers {
                f(evaluator)?;
            }
        }
        if let Some(multipliers) = &mut self.threshold_multipliers {
            multipliers.for_each_generic_evaluator_mut(&mut f)?;
        }

        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub fn evaluate<T: FloatLike>(
        &mut self,
        momentum_sample: &MomentumSample<T>,
        graph: &Graph,
        model: &Model,
        esurfaces: &EsurfaceCollection,
        rotation: &Rotation,
        settings: &RuntimeSettings,
        param_builder: &mut ParamBuilder<f64>,
        orientation: SingleOrAllOrientations<'_, OrientationID>,
        evaluation_metadata: &mut EvaluationMetaData,
        record_primary_timing: bool,
        record_components: bool,
    ) -> Result<AmplitudeCountertermEvaluation<T>> {
        if record_components && self.metadata_registry.is_none() {
            return Err(eyre!(
                "Amplitude graph '{}' requested detailed threshold components without metadata",
                graph.name,
            ));
        }
        if !self.legacy_equivalent {
            return self.evaluate_variants(
                momentum_sample,
                graph,
                model,
                esurfaces,
                rotation,
                settings,
                param_builder,
                orientation,
                evaluation_metadata,
                record_primary_timing,
                record_components,
            );
        }

        debug!("start evaluate threshold counterterm");
        let existing_esurfaces = self
            .overlap
            .existing_esurfaces
            .iter()
            .map(|e| e.0)
            .collect::<Vec<_>>();
        debug!("subtracting esurfaces: {:?}", existing_esurfaces);
        debug!("overlap structure\n: {}", self.overlap);

        let counter_term_builder = CounterTermBuilder::new(
            graph,
            model,
            rotation,
            settings,
            esurfaces,
            momentum_sample,
            &self.overlap,
            &self.raised_data,
            self.own_group_position,
            &self.esurface_map,
        );

        let mut result = Complex::new_re(momentum_sample.zero());
        let mut local_counterterms = Vec::new();
        let mut components = record_components.then(Vec::new);
        let mut radial_root_failed = false;
        let split_helpers = self.metadata_registry.is_some();

        for (overlap_group, group) in self.overlap.overlap_groups.iter().enumerate() {
            let overlap_builder = counter_term_builder.new_overlap_builder(group);

            for existing_esurface_id in group.existing_esurfaces.iter() {
                let group_esurface_id = self.overlap.existing_esurfaces[*existing_esurface_id];
                if self
                    .local_esurface_exists
                    .get(group_esurface_id)
                    .is_some_and(|exists| !*exists)
                {
                    continue;
                }

                let Some(esurface_builder) =
                    overlap_builder.new_esurface_builder(*existing_esurface_id)
                else {
                    continue;
                };

                let raised_esurface_id = esurface_builder.raised_esurface_id;
                self.ensure_active_raised_esurface(raised_esurface_id)?;
                let radial_root_identity = Self::radial_root_identity(
                    &graph.name,
                    overlap_group,
                    raised_esurface_id,
                    rotation,
                );
                let Some(rstar_solution) = esurface_builder.solve_rstar(
                    &radial_root_identity,
                    &mut evaluation_metadata.radial_root_diagnostics,
                ) else {
                    evaluation_metadata.record_threshold_counterterm_error(format!(
                        "amplitude graph '{}' overlap group {} raised E-surface {} failed center or radial-root validation in probe rotation {}",
                        graph.name,
                        overlap_group,
                        raised_esurface_id.0,
                        rotation.method,
                    ));
                    radial_root_failed = true;
                    continue;
                };
                if radial_root_failed {
                    continue;
                }
                let single_evaluation = rstar_solution.rstar_samples().evaluate(
                    param_builder,
                    orientation,
                    evaluation_metadata,
                    record_primary_timing,
                    &mut self.evaluators[raised_esurface_id],
                    &mut self.helper_evaluators,
                    split_helpers,
                    record_components,
                )?;
                let single_result = single_evaluation.weighted.clone();

                if record_components {
                    let variant_id = self
                        .metadata_variant_for_raised_esurface(raised_esurface_id)
                        .with_context(|| {
                            format!(
                                "while recording amplitude graph '{}' overlap group {}",
                                graph.name, overlap_group,
                            )
                        })?;
                    let one = momentum_sample.one();
                    append_amplitude_component_evaluations(
                        &mut components,
                        &single_evaluation,
                        variant_id,
                        raised_esurface_id,
                        overlap_group,
                        &one,
                        &one,
                    );
                }

                if !single_result.is_zero() {
                    //    debug!(
                    //        "Param Builder for {}:\n{}",
                    //        existing_esurface_id, param_builder
                    //    );
                    debug!(
                        "Counterterm for esurface {}: {:+16e}",
                        existing_esurface_id, single_result
                    );
                }

                local_counterterms.push(AmplitudeLocalCountertermEvaluation {
                    variant_id: None,
                    esurface_id: raised_esurface_id,
                    overlap_group,
                    value: single_result.clone(),
                });
                result += single_result;
            }
        }
        if radial_root_failed {
            return Ok(AmplitudeCountertermEvaluation {
                total: Complex::new_re(F::from_f64(f64::NAN)),
                local_counterterms,
                components: None,
            });
        }
        Ok(AmplitudeCountertermEvaluation {
            total: result,
            local_counterterms,
            components,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn evaluate_variants<T: FloatLike>(
        &mut self,
        momentum_sample: &MomentumSample<T>,
        graph: &Graph,
        model: &Model,
        esurfaces: &EsurfaceCollection,
        rotation: &Rotation,
        settings: &RuntimeSettings,
        param_builder: &mut ParamBuilder<f64>,
        orientation: SingleOrAllOrientations<'_, OrientationID>,
        evaluation_metadata: &mut EvaluationMetaData,
        record_primary_timing: bool,
        record_components: bool,
    ) -> Result<AmplitudeCountertermEvaluation<T>> {
        if self.variant_evaluators.len() != self.variant_subspaces.len()
            || self.variant_evaluators.len() != self.variant_raised_esurfaces.len()
            || self.variant_evaluators.len() != self.variant_helper_evaluators.len()
            || self.variant_evaluators.len() != self.variant_generated_mask.len()
            || self.variant_evaluators.len() != self.variant_active_mask.len()
            || self.variant_evaluators.len() != self.variant_metadata.len()
        {
            return Err(eyre!(
                "Amplitude graph '{}' has inconsistent threshold-variant runtime dimensions",
                graph.name,
            ));
        }
        if let Some(multipliers) = &self.threshold_multipliers {
            multipliers.validate(self.variant_evaluators.len(), 0)?;
            if multipliers
                .left_variants()
                .iter()
                .enumerate()
                .any(|(variant_id, reference)| reference.variant_id.0 != variant_id)
            {
                return Err(eyre!(
                    "Amplitude graph '{}' threshold multiplier references are not ordered by contiguous variant ID",
                    graph.name,
                ));
            }
        }

        let candidate_variant_ids = self
            .variant_evaluators
            .keys()
            .filter(|&variant_id| {
                self.variant_generated_mask[variant_id] && self.variant_active_mask[variant_id]
            })
            .collect::<Vec<_>>();
        if candidate_variant_ids.is_empty() {
            return Ok(AmplitudeCountertermEvaluation {
                total: Complex::new_re(momentum_sample.zero()),
                local_counterterms: Vec::new(),
                components: record_components.then(Vec::new),
            });
        }
        for &variant_id in &candidate_variant_ids {
            self.ensure_active_variant(variant_id)?;
        }

        // The union selects only a common coordinate frame for the SOCP. Each projected
        // E-surface below retains its own exact subspace and sampled fixed complement.
        let common_subspace = SubspaceData::union_in_common_parent(
            candidate_variant_ids
                .iter()
                .map(|&variant_id| &self.variant_subspaces[variant_id]),
            graph,
            &self.lmbs,
        )
        .with_context(|| {
            format!(
                "Amplitude graph '{}' threshold variants cannot share an overlap-center frame",
                graph.name,
            )
        })?;
        let common_lmb = common_subspace.get_lmb(&self.lmbs);
        let sample_in_common_lmb =
            momentum_sample.lmb_transform(&graph.loop_momentum_basis, common_lmb);
        let real_mass_vector = graph.get_real_mass_vector(model);
        let e_cm = F::from_f64(settings.kinematics.e_cm);
        let existence_tolerance = F::from_f64(settings.subtraction.esurface_existence_threshold);

        let mut thresholds = EsurfaceCollection::new();
        let mut threshold_subspaces = Vec::new();
        let mut synthetic_variants = TiVec::<EsurfaceID, ThresholdCountertermVariantId>::new();
        for variant_id in candidate_variant_ids {
            let raised_esurface_id = self.variant_raised_esurfaces[variant_id];
            let esurface_id = self.raised_data.raised_groups[raised_esurface_id].esurface_ids[0];
            let esurface = &esurfaces[esurface_id];
            let subspace = &self.variant_subspaces[variant_id];
            let existence = esurface.classify_existence_subspace(
                sample_in_common_lmb.loop_moms(),
                sample_in_common_lmb.external_moms(),
                subspace,
                &self.lmbs,
                graph,
                &real_mass_vector,
                &[],
                &e_cm,
                &existence_tolerance,
            );
            if !existence.is_existing() {
                continue;
            }
            thresholds.push(esurface.clone());
            threshold_subspaces.push(subspace.clone());
            synthetic_variants.push(variant_id);
        }
        if thresholds.is_empty() {
            return Ok(AmplitudeCountertermEvaluation {
                total: Complex::new_re(momentum_sample.zero()),
                local_counterterms: Vec::new(),
                components: record_components.then(Vec::new),
            });
        }

        let existing_thresholds: ExistingThresholds = thresholds.keys().collect();
        let sample_f64 = momentum_sample_to_f64(&sample_in_common_lmb);
        let overlap_input = SubspaceOverlapInput {
            graph,
            settings,
            subspace: &common_subspace,
            threshold_subspaces: Some(&threshold_subspaces),
            lmbs: &self.lmbs,
            thresholds: &thresholds,
            edge_masses: graph.get_real_mass_vector::<f64>(model),
        };
        let overlap = find_maximal_subspace_overlap(
            &overlap_input,
            &existing_thresholds,
            sample_f64.loop_moms(),
            sample_f64.external_moms(),
            rotation,
        )
        .with_context(|| {
            format!(
                "Failed to find graph-local projected threshold overlap for amplitude graph '{}'",
                graph.name,
            )
        })?;

        let max_prefactor_power = synthetic_variants
            .iter()
            .map(|&variant_id| {
                self.raised_data.raised_groups[self.variant_raised_esurfaces[variant_id]]
                    .max_occurence
                    + 1
            })
            .max()
            .unwrap_or(1);
        let mut total = Complex::new_re(momentum_sample.zero());
        let mut local_counterterms = Vec::new();
        let mut components = record_components.then(Vec::new);
        let mut radial_root_failed = false;
        let mut multiplier_workspace = self.threshold_multipliers.as_ref().map(|multipliers| {
            ThresholdMultiplierInputWorkspace::new(multipliers.layout(), momentum_sample.zero())
        });
        let additional_params = settings.additional_params::<T>();

        for (overlap_group_index, overlap_group) in overlap.overlap_groups.iter().enumerate() {
            for &existing_esurface_id in &overlap_group.existing_esurfaces {
                let synthetic_esurface_id = overlap.existing_esurfaces[existing_esurface_id];
                let variant_id = synthetic_variants[synthetic_esurface_id];
                let raised_esurface_id = self.variant_raised_esurfaces[variant_id];
                let subspace = &self.variant_subspaces[variant_id];
                let esurface = &thresholds[synthetic_esurface_id];

                let mut center = overlap_group.center.cast::<T>();
                for loop_index in (0..center.0.len()).map(LoopIndex::from) {
                    if !subspace.contains_loop_index(loop_index) {
                        center[loop_index] = ThreeMomentum::new(
                            momentum_sample.zero(),
                            momentum_sample.zero(),
                            momentum_sample.zero(),
                        );
                    }
                }
                let shifted_loop_momenta = sample_in_common_lmb.loop_moms() - &center;
                let radius = shifted_loop_momenta
                    .hyper_radius_squared(subspace.as_subspace_simple())
                    .sqrt();
                let unit_shifted_momenta =
                    shifted_loop_momenta.rescale(&radius.inv(), subspace.as_subspace_simple());
                let (raw_radius_guess, _) = esurface.get_radius_guess_subspace(
                    &unit_shifted_momenta,
                    sample_in_common_lmb.external_moms(),
                    subspace,
                    &self.lmbs,
                    graph,
                    &real_mass_vector,
                );
                let function = |r: &_| {
                    esurface.compute_self_and_r_derivative_subspace(
                        r,
                        &unit_shifted_momenta,
                        &center,
                        sample_in_common_lmb.external_moms(),
                        &real_mass_vector,
                        subspace,
                        &self.lmbs,
                        graph,
                    )
                };
                let zero = raw_radius_guess.zero();
                let mut radius_guess = raw_radius_guess.clone();
                if radius_guess.is_nan() || radius_guess.is_infinite() || radius_guess <= zero {
                    radius_guess = e_cm.clone();
                }
                let identity = Self::variant_radial_root_identity(
                    &graph.name,
                    overlap_group_index,
                    variant_id,
                    raised_esurface_id,
                    rotation,
                );
                let tolerance = F::from_f64(settings.subtraction.radial_root_residual_tolerance);
                let solution = match evaluation_metadata.radial_root_diagnostics.solve(
                    &identity,
                    &zero,
                    &radius_guess,
                    function,
                    &tolerance,
                    MAX_ITERATIONS,
                    64,
                    &e_cm,
                ) {
                    Ok(solution) => solution,
                    Err(error) => {
                        evaluation_metadata.record_threshold_counterterm_error(format!(
                            "amplitude graph '{}' threshold variant {} overlap group {} failed radial-root validation: {error}",
                            graph.name, variant_id.0, overlap_group_index,
                        ));
                        radial_root_failed = true;
                        continue;
                    }
                };
                if radial_root_failed {
                    continue;
                }

                let root_loop_momenta = &unit_shifted_momenta
                    .rescale(&solution.solution, subspace.as_subspace_simple())
                    + &center;
                let mut root_in_common_lmb = sample_in_common_lmb.clone();
                root_in_common_lmb.sample.loop_moms = root_loop_momenta;
                let root_in_generation_lmb =
                    root_in_common_lmb.lmb_transform(common_lmb, &graph.loop_momentum_basis);
                let (local_multiplier, integrated_multiplier) = Self::evaluate_variant_multiplier(
                    &mut self.threshold_multipliers,
                    multiplier_workspace.as_mut(),
                    variant_id,
                    graph,
                    &real_mass_vector,
                    param_builder.model_values(),
                    &additional_params,
                    momentum_sample,
                    &root_in_generation_lmb,
                    evaluation_metadata,
                    record_primary_timing,
                )?;
                let both_zero = local_multiplier.is_zero() && integrated_multiplier.is_zero();
                let evaluation = if both_zero {
                    EvaluatedAmplitudeThreshold::new(momentum_sample.zero(), record_components)
                        .finish(&local_multiplier, &integrated_multiplier)
                } else {
                    evaluate_generalized_rstar(
                        graph,
                        settings,
                        param_builder,
                        orientation,
                        evaluation_metadata,
                        record_primary_timing,
                        &real_mass_vector,
                        &self.lmbs,
                        common_lmb,
                        &root_in_common_lmb,
                        &center,
                        &unit_shifted_momenta,
                        &radius,
                        &solution,
                        esurface,
                        synthetic_esurface_id,
                        raised_esurface_id,
                        subspace,
                        &thresholds,
                        &threshold_subspaces,
                        &overlap,
                        overlap_group_index,
                        max_prefactor_power,
                        &mut self.variant_evaluators[variant_id],
                        &mut self.variant_helper_evaluators[variant_id],
                        &local_multiplier,
                        &integrated_multiplier,
                        record_components,
                    )?
                };
                let value = evaluation.weighted.clone();
                append_amplitude_component_evaluations(
                    &mut components,
                    &evaluation,
                    variant_id,
                    raised_esurface_id,
                    overlap_group_index,
                    &local_multiplier,
                    &integrated_multiplier,
                );

                local_counterterms.push(AmplitudeLocalCountertermEvaluation {
                    variant_id: Some(variant_id),
                    esurface_id: raised_esurface_id,
                    overlap_group: overlap_group_index,
                    value: value.clone(),
                });
                total += value;
            }
        }

        if radial_root_failed {
            return Ok(AmplitudeCountertermEvaluation {
                total: Complex::new_re(F::from_f64(f64::NAN)),
                local_counterterms,
                components: None,
            });
        }
        Ok(AmplitudeCountertermEvaluation {
            total,
            local_counterterms,
            components,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn kinematics_for_variant_approach<T: FloatLike>(
        &self,
        momentum_sample: &MomentumSample<T>,
        graph: &Graph,
        model: &Model,
        esurfaces: &EsurfaceCollection,
        rotation: &Rotation,
        settings: &RuntimeSettings,
    ) -> Result<OverlapStructureWithKinematics<T>> {
        if self.variant_evaluators.len() != self.variant_subspaces.len()
            || self.variant_evaluators.len() != self.variant_raised_esurfaces.len()
            || self.variant_evaluators.len() != self.variant_generated_mask.len()
            || self.variant_evaluators.len() != self.variant_active_mask.len()
        {
            return Err(eyre!(
                "Amplitude graph '{}' has inconsistent threshold-variant approach dimensions",
                graph.name,
            ));
        }

        let candidate_variant_ids = self
            .variant_evaluators
            .keys()
            .filter(|&variant_id| {
                self.variant_generated_mask[variant_id] && self.variant_active_mask[variant_id]
            })
            .collect::<Vec<_>>();
        if candidate_variant_ids.is_empty() {
            return Ok(OverlapStructureWithKinematics {
                existing_esurfaces: ExistingEsurfaces::new(),
                variant_ids: Some(TiVec::new()),
                overlap_groups_with_kinematics: Vec::new(),
            });
        }
        for &variant_id in &candidate_variant_ids {
            self.ensure_active_variant(variant_id)?;
        }

        let common_subspace = SubspaceData::union_in_common_parent(
            candidate_variant_ids
                .iter()
                .map(|&variant_id| &self.variant_subspaces[variant_id]),
            graph,
            &self.lmbs,
        )
        .with_context(|| {
            format!(
                "Amplitude graph '{}' threshold variants cannot share an approach-center frame",
                graph.name,
            )
        })?;
        let common_lmb = common_subspace.get_lmb(&self.lmbs);
        let sample_in_common_lmb =
            momentum_sample.lmb_transform(&graph.loop_momentum_basis, common_lmb);
        let real_mass_vector = graph.get_real_mass_vector(model);
        let e_cm = F::from_f64(settings.kinematics.e_cm);
        let existence_tolerance = F::from_f64(settings.subtraction.esurface_existence_threshold);

        let mut thresholds = EsurfaceCollection::new();
        let mut threshold_subspaces = Vec::new();
        let mut synthetic_variants = TiVec::<EsurfaceID, ThresholdCountertermVariantId>::new();
        let mut existing_esurfaces = ExistingEsurfaces::new();
        for variant_id in candidate_variant_ids {
            let raised_esurface_id = self.variant_raised_esurfaces[variant_id];
            let esurface_id = self.raised_data.raised_groups[raised_esurface_id].esurface_ids[0];
            let esurface = &esurfaces[esurface_id];
            let subspace = &self.variant_subspaces[variant_id];
            let existence = esurface.classify_existence_subspace(
                sample_in_common_lmb.loop_moms(),
                sample_in_common_lmb.external_moms(),
                subspace,
                &self.lmbs,
                graph,
                &real_mass_vector,
                &[],
                &e_cm,
                &existence_tolerance,
            );
            if !existence.is_existing() {
                continue;
            }
            let group_esurface_id = self
                .esurface_map
                .iter_enumerated()
                .find_map(|(group_esurface_id, graph_map)| {
                    (graph_map[self.own_group_position] == Some(raised_esurface_id))
                        .then_some(group_esurface_id)
                })
                .ok_or_else(|| {
                    eyre!(
                        "Amplitude graph '{}' threshold variant {} raised E-surface {} has no group E-surface mapping",
                        graph.name,
                        variant_id.0,
                        raised_esurface_id.0,
                    )
                })?;
            thresholds.push(esurface.clone());
            threshold_subspaces.push(subspace.clone());
            synthetic_variants.push(variant_id);
            existing_esurfaces.push(group_esurface_id);
        }

        let variant_ids = synthetic_variants.iter().copied().collect::<TiVec<_, _>>();
        if thresholds.is_empty() {
            return Ok(OverlapStructureWithKinematics {
                existing_esurfaces,
                variant_ids: Some(variant_ids),
                overlap_groups_with_kinematics: Vec::new(),
            });
        }

        let existing_thresholds: ExistingThresholds = thresholds.keys().collect();
        let sample_f64 = momentum_sample_to_f64(&sample_in_common_lmb);
        let overlap_input = SubspaceOverlapInput {
            graph,
            settings,
            subspace: &common_subspace,
            threshold_subspaces: Some(&threshold_subspaces),
            lmbs: &self.lmbs,
            thresholds: &thresholds,
            edge_masses: graph.get_real_mass_vector::<f64>(model),
        };
        let overlap = find_maximal_subspace_overlap(
            &overlap_input,
            &existing_thresholds,
            sample_f64.loop_moms(),
            sample_f64.external_moms(),
            rotation,
        )
        .with_context(|| {
            format!(
                "Failed to find graph-local projected threshold approach overlap for amplitude graph '{}'",
                graph.name,
            )
        })?;

        let mut radial_root_diagnostics = RadialRootDiagnostics::default();
        let mut overlap_groups_with_kinematics = Vec::with_capacity(overlap.overlap_groups.len());
        for (overlap_group_index, overlap_group) in overlap.overlap_groups.iter().enumerate() {
            let mut group_center_in_common_lmb = sample_f64.clone();
            for loop_index in common_subspace.iter_lmb_indices() {
                group_center_in_common_lmb.sample.loop_moms[loop_index] =
                    overlap_group.center[loop_index];
            }
            let group_center = group_center_in_common_lmb
                .lmb_transform(common_lmb, &graph.loop_momentum_basis)
                .loop_moms()
                .clone();
            let mut loop_momenta_at_esurface = TiVec::new();
            let mut approach_centers_at_esurface = TiVec::new();

            for &existing_esurface_id in &overlap_group.existing_esurfaces {
                let synthetic_esurface_id = overlap.existing_esurfaces[existing_esurface_id];
                let variant_id = synthetic_variants[synthetic_esurface_id];
                let raised_esurface_id = self.variant_raised_esurfaces[variant_id];
                let subspace = &self.variant_subspaces[variant_id];
                let esurface = &thresholds[synthetic_esurface_id];

                let mut projected_center = overlap_group.center.cast::<T>();
                for loop_index in (0..projected_center.0.len()).map(LoopIndex::from) {
                    if !subspace.contains_loop_index(loop_index) {
                        projected_center[loop_index] = ThreeMomentum::new(
                            momentum_sample.zero(),
                            momentum_sample.zero(),
                            momentum_sample.zero(),
                        );
                    }
                }
                let shifted_loop_momenta = sample_in_common_lmb.loop_moms() - &projected_center;
                let radius = shifted_loop_momenta
                    .hyper_radius_squared(subspace.as_subspace_simple())
                    .sqrt();
                let unit_shifted_momenta =
                    shifted_loop_momenta.rescale(&radius.inv(), subspace.as_subspace_simple());
                let (raw_radius_guess, _) = esurface.get_radius_guess_subspace(
                    &unit_shifted_momenta,
                    sample_in_common_lmb.external_moms(),
                    subspace,
                    &self.lmbs,
                    graph,
                    &real_mass_vector,
                );
                let function = |r: &_| {
                    esurface.compute_self_and_r_derivative_subspace(
                        r,
                        &unit_shifted_momenta,
                        &projected_center,
                        sample_in_common_lmb.external_moms(),
                        &real_mass_vector,
                        subspace,
                        &self.lmbs,
                        graph,
                    )
                };
                let zero = raw_radius_guess.zero();
                let mut radius_guess = raw_radius_guess.clone();
                if radius_guess.is_nan() || radius_guess.is_infinite() || radius_guess <= zero {
                    radius_guess = e_cm.clone();
                }
                let identity = Self::variant_radial_root_identity(
                    &graph.name,
                    overlap_group_index,
                    variant_id,
                    raised_esurface_id,
                    rotation,
                );
                let tolerance = F::from_f64(settings.subtraction.radial_root_residual_tolerance);
                let solution = radial_root_diagnostics
                    .solve(
                        &identity,
                        &zero,
                        &radius_guess,
                        function,
                        &tolerance,
                        MAX_ITERATIONS,
                        64,
                        &e_cm,
                    )
                    .map_err(|error| {
                        eyre!(
                            "Could not construct threshold-counterterm approach kinematics for amplitude graph '{}' variant {} raised E-surface {} in overlap group {}: {error}",
                            graph.name,
                            variant_id.0,
                            raised_esurface_id.0,
                            overlap_group_index,
                        )
                    })?;

                let root_loop_momenta = &unit_shifted_momenta
                    .rescale(&solution.solution, subspace.as_subspace_simple())
                    + &projected_center;
                let mut root_in_common_lmb = sample_in_common_lmb.clone();
                root_in_common_lmb.sample.loop_moms = root_loop_momenta;
                let root_in_generation_lmb =
                    root_in_common_lmb.lmb_transform(common_lmb, &graph.loop_momentum_basis);

                let mut full_center_in_common_lmb = sample_in_common_lmb.clone();
                for loop_index in subspace.iter_lmb_indices() {
                    full_center_in_common_lmb.sample.loop_moms[loop_index] =
                        projected_center[loop_index].clone();
                }
                let full_center_in_generation_lmb = full_center_in_common_lmb
                    .lmb_transform(common_lmb, &graph.loop_momentum_basis)
                    .loop_moms()
                    .clone();

                loop_momenta_at_esurface.push(Some(root_in_generation_lmb));
                approach_centers_at_esurface.push(Some(full_center_in_generation_lmb));
            }

            overlap_groups_with_kinematics.push(OverlapGroupWithKinematics {
                overlap_group: OverlapGroup {
                    existing_esurfaces: overlap_group.existing_esurfaces.clone(),
                    complement: overlap_group.complement.clone(),
                    center: group_center,
                    prefactor_evaluator: None,
                },
                loop_momenta_at_esurface,
                approach_centers_at_esurface: Some(approach_centers_at_esurface),
            });
        }

        Ok(OverlapStructureWithKinematics {
            existing_esurfaces,
            variant_ids: Some(variant_ids),
            overlap_groups_with_kinematics,
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub fn kinematics_for_approach<T: FloatLike>(
        &mut self,
        momentum_sample: &MomentumSample<T>,
        graph: &Graph,
        model: &Model,
        esurfaces: &EsurfaceCollection,
        rotation: &Rotation,
        settings: &RuntimeSettings,
    ) -> Result<OverlapStructureWithKinematics<T>> {
        if !self.legacy_equivalent {
            return self.kinematics_for_variant_approach(
                momentum_sample,
                graph,
                model,
                esurfaces,
                rotation,
                settings,
            );
        }

        let counter_term_builder = CounterTermBuilder::new(
            graph,
            model,
            rotation,
            settings,
            esurfaces,
            momentum_sample,
            &self.overlap,
            &self.raised_data,
            self.own_group_position,
            &self.esurface_map,
        );

        let mut overlap_sturcture_with_kinematics = OverlapStructureWithKinematics {
            existing_esurfaces: self.overlap.existing_esurfaces.clone(),
            variant_ids: None,
            overlap_groups_with_kinematics: Vec::new(),
        };
        let mut radial_root_diagnostics = RadialRootDiagnostics::default();

        for (overlap_group, group) in self.overlap.overlap_groups.iter().enumerate() {
            let overlap_builder = counter_term_builder.new_overlap_builder(group);

            let mut loop_momenta_at_esurface: TiVec<ExistingEsurfaceId, Option<MomentumSample<T>>> =
                ti_vec![];
            for existing_esurface_id in group.existing_esurfaces.iter() {
                let single_result = overlap_builder
                    .new_esurface_builder(*existing_esurface_id)
                    .map(|esurface_builder| -> Result<_> {
                        self.ensure_active_raised_esurface(esurface_builder.raised_esurface_id)?;
                        let raised_esurface_id = esurface_builder.raised_esurface_id;
                        let radial_root_identity = Self::radial_root_identity(
                            &graph.name,
                            overlap_group,
                            raised_esurface_id,
                            rotation,
                        );
                        let rstar_sample = esurface_builder
                            .solve_rstar(
                                &radial_root_identity,
                                &mut radial_root_diagnostics,
                            )
                            .ok_or_else(|| {
                                eyre!(
                                    "Could not construct threshold-counterterm kinematics for raised E-surface {} in probe rotation {}",
                                    raised_esurface_id.0,
                                    rotation.method,
                                )
                            })?
                            .rstar_samples();
                        Result::Ok(rstar_sample.rstar_sample)
                    })
                    .transpose()?;

                loop_momenta_at_esurface.push(single_result);
            }

            overlap_sturcture_with_kinematics
                .overlap_groups_with_kinematics
                .push(OverlapGroupWithKinematics {
                    overlap_group: group.clone(),
                    loop_momenta_at_esurface,
                    approach_centers_at_esurface: None,
                });
        }

        Ok(overlap_sturcture_with_kinematics)
    }
}

pub struct OverlapGroupWithKinematics<T: FloatLike> {
    pub overlap_group: OverlapGroup,
    pub loop_momenta_at_esurface: TiVec<ExistingEsurfaceId, Option<MomentumSample<T>>>,
    /// Full fixed-complement centers aligned with `loop_momenta_at_esurface`.
    /// `None` retains the homogeneous legacy group center exactly.
    pub approach_centers_at_esurface: Option<TiVec<ExistingEsurfaceId, Option<LoopMomenta<F<T>>>>>,
}

pub struct OverlapStructureWithKinematics<T: FloatLike> {
    pub existing_esurfaces: ExistingEsurfaces,
    /// Stable semantic variant IDs aligned with `existing_esurfaces` on the generalized path.
    /// `None` identifies the untouched legacy group-E-surface representation.
    pub variant_ids: Option<TiVec<ExistingEsurfaceId, ThresholdCountertermVariantId>>,
    pub overlap_groups_with_kinematics: Vec<OverlapGroupWithKinematics<T>>,
}

fn momentum_sample_to_f64<T: FloatLike>(sample: &MomentumSample<T>) -> MomentumSample<f64> {
    MomentumSample {
        sample: BareMomentumSample {
            loop_moms: sample
                .loop_moms()
                .iter()
                .map(|momentum| momentum.map_ref(&F::into_ff64))
                .collect(),
            dual_loop_moms: None,
            loop_mom_cache_id: sample.sample.loop_mom_cache_id,
            loop_mom_base_cache_id: sample.sample.loop_mom_base_cache_id,
            external_moms: sample
                .external_moms()
                .iter()
                .map(|momentum| momentum.map_ref(&F::into_ff64))
                .collect(),
            external_mom_cache_id: sample.sample.external_mom_cache_id,
            external_mom_base_cache_id: sample.sample.external_mom_base_cache_id,
            jacobian: sample.sample.jacobian.into_ff64(),
            orientation: sample.sample.orientation,
            parameterization_branch: sample.sample.parameterization_branch,
        },
    }
}

fn dualize_external_momenta<T: FloatLike>(
    shape: &HyperDual<F<T>>,
    external_moms: &ExternalFourMomenta<F<T>>,
) -> ExternalFourMomenta<HyperDual<F<T>>> {
    external_moms
        .iter()
        .map(|momentum| momentum.map_ref(&|component| new_constant(shape, component)))
        .collect()
}

fn real_dual_to_complex<T: FloatLike>(dual: HyperDual<F<T>>) -> HyperDual<Complex<F<T>>> {
    let shape = dual
        .get_shape()
        .iter()
        .map(|entry| entry.to_vec())
        .collect();
    let values = dual.values.into_iter().map(Complex::new_re).collect();
    HyperDual::from_values(shape, values)
}

fn repeated_product<T>(value: T, power: usize, one: T) -> T
where
    T: Clone + std::ops::Mul<T, Output = T>,
{
    (0..power).fold(one, |product, _| product * value.clone())
}

#[allow(clippy::too_many_arguments)]
fn evaluate_generalized_multichanneling_prefactor<T: FloatLike>(
    momentum_sample: &MomentumSample<T>,
    real_mass_vector: &EdgeVec<F<T>>,
    all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
    thresholds: &EsurfaceCollection,
    threshold_subspaces: &[SubspaceData],
    overlap: &SubspaceOverlapStructure,
    overlap_group_index: usize,
    power: usize,
) -> DualOrNot<Complex<F<T>>> {
    if overlap.overlap_groups.len() < 2 {
        return DualOrNot::NonDual(Complex::new_re(momentum_sample.one()));
    }

    if let Some(dual_loop_momenta) = &momentum_sample.sample.dual_loop_moms {
        let shape = &dual_loop_momenta[LoopIndex(0)].px;
        let dual_external_momenta =
            dualize_external_momenta(shape, momentum_sample.external_moms());
        let group_products = overlap
            .overlap_groups
            .iter()
            .map(|group| {
                group.complement.iter().fold(
                    new_constant(shape, &momentum_sample.one()),
                    |product, &existing_esurface_id| {
                        let esurface_id = overlap.existing_esurfaces[existing_esurface_id];
                        let subspace = &threshold_subspaces[esurface_id.0];
                        let value = thresholds[esurface_id].compute_from_dual_momenta(
                            subspace.get_lmb(all_lmbs),
                            real_mass_vector,
                            dual_loop_momenta,
                            &dual_external_momenta,
                        );
                        product
                            * repeated_product(
                                value,
                                power,
                                new_constant(shape, &momentum_sample.one()),
                            )
                    },
                )
            })
            .collect::<Vec<_>>();
        let denominator = group_products
            .iter()
            .skip(1)
            .fold(group_products[0].clone(), |sum, product| sum + product);
        return DualOrNot::Dual(real_dual_to_complex(
            group_products[overlap_group_index].clone() / denominator,
        ));
    }

    let group_products = overlap
        .overlap_groups
        .iter()
        .map(|group| {
            group
                .complement
                .iter()
                .map(|&existing_esurface_id| {
                    let esurface_id = overlap.existing_esurfaces[existing_esurface_id];
                    let subspace = &threshold_subspaces[esurface_id.0];
                    let value = thresholds[esurface_id].compute_from_momenta(
                        subspace.get_lmb(all_lmbs),
                        real_mass_vector,
                        momentum_sample.loop_moms(),
                        momentum_sample.external_moms(),
                    );
                    repeated_product(value, power, momentum_sample.one())
                })
                .fold(momentum_sample.one(), |product, value| product * value)
        })
        .collect::<Vec<_>>();
    let denominator = group_products
        .iter()
        .skip(1)
        .fold(group_products[0].clone(), |sum, value| sum + value);
    DualOrNot::NonDual(Complex::new_re(
        group_products[overlap_group_index].clone() / denominator,
    ))
}

#[allow(clippy::too_many_arguments)]
fn evaluate_generalized_rstar<T: FloatLike>(
    graph: &Graph,
    settings: &RuntimeSettings,
    param_builder: &mut ParamBuilder<f64>,
    orientations: SingleOrAllOrientations<'_, OrientationID>,
    evaluation_metadata: &mut EvaluationMetaData,
    record_primary_timing: bool,
    real_mass_vector: &EdgeVec<F<T>>,
    all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
    common_lmb: &LoopMomentumBasis,
    root_sample: &MomentumSample<T>,
    center: &LoopMomenta<F<T>>,
    unit_shifted_momenta: &LoopMomenta<F<T>>,
    radius: &F<T>,
    solution: &NewtonIterationResult<T>,
    esurface: &Esurface,
    esurface_id: EsurfaceID,
    raised_esurface_id: RaisedEsurfaceId,
    subspace: &SubspaceData,
    thresholds: &EsurfaceCollection,
    threshold_subspaces: &[SubspaceData],
    overlap: &SubspaceOverlapStructure,
    overlap_group_index: usize,
    prefactor_power: usize,
    ct_evaluator: &mut AmplitudeCountertermEvaluator,
    helper_evaluators: &mut [GenericEvaluator],
    local_multiplier: &F<T>,
    integrated_multiplier: &F<T>,
    record_components: bool,
) -> Result<EvaluatedAmplitudeThreshold<T>> {
    let radius_star = &solution.solution;
    let e_cm = F::from_f64(settings.kinematics.e_cm);
    let uv_settings = &settings.subtraction.local_ct_settings.uv_localisation;
    let integrated_settings = &settings.subtraction.integrated_ct_settings;
    let uv_damp_plus = evaluate_uv_damper(radius, radius_star, &e_cm, uv_settings);
    let uv_damp_minus = evaluate_uv_damper(&-radius, radius_star, &e_cm, uv_settings);
    let h_function =
        evaluate_integrated_ct_normalisation(radius, radius_star, &e_cm, integrated_settings);
    let mut evaluation = EvaluatedAmplitudeThreshold::new(root_sample.zero(), record_components);

    for (cut_cff_index, evaluator_stack) in &mut ct_evaluator.evaluator_stacks {
        let order_index = cut_cff_index.left_threshold_order.unwrap() - 1;
        let (sample_for_order_in_common_lmb, threshold_params) = if order_index == 0 {
            (
                root_sample.clone(),
                ThresholdParams {
                    radius: DualOrNot::NonDual(radius.clone()),
                    radius_star: DualOrNot::NonDual(radius_star.clone()),
                    esurface_derivative: DualOrNot::NonDual(
                        solution.derivative_at_solution.clone(),
                    ),
                    uv_damp_plus: DualOrNot::NonDual(uv_damp_plus.clone()),
                    uv_damp_minus: DualOrNot::NonDual(uv_damp_minus.clone()),
                    h_function: DualOrNot::NonDual(h_function.clone()),
                },
            )
        } else {
            let dual_shape = HyperDual::<F<T>>::new(simple_n_deriv_shape(order_index));
            let dual_radius_star = dual_shape.variable(0, radius_star.clone());
            let dual_center = center
                .iter()
                .map(|momentum| momentum.map_ref(&|value| new_constant(&dual_radius_star, value)))
                .collect::<LoopMomenta<_>>();
            let dual_loop_momenta = unit_shifted_momenta
                .rescale_with_hyper_dual(&dual_radius_star, subspace.as_subspace_simple())
                .iter()
                .zip(dual_center.iter())
                .map(|(momentum, center)| momentum.clone() + center.clone())
                .collect::<LoopMomenta<_>>();
            let mut sample = root_sample.clone();
            sample.sample.dual_loop_moms = Some(dual_loop_momenta);
            (
                sample,
                ThresholdParams {
                    radius: DualOrNot::Dual(new_constant(&dual_radius_star, radius)),
                    radius_star: DualOrNot::Dual(dual_radius_star.clone()),
                    esurface_derivative: DualOrNot::Dual(new_constant(
                        &dual_radius_star,
                        &solution.derivative_at_solution,
                    )),
                    uv_damp_plus: DualOrNot::Dual(new_constant(&dual_radius_star, &uv_damp_plus)),
                    uv_damp_minus: DualOrNot::Dual(new_constant(&dual_radius_star, &uv_damp_minus)),
                    h_function: DualOrNot::Dual(new_constant(&dual_radius_star, &h_function)),
                },
            )
        };

        let esurface_derivatives = if order_index == 0 {
            DualOrNot::NonDual(solution.derivative_at_solution.clone())
        } else {
            let dual_shape = HyperDual::<F<T>>::new(simple_n_deriv_shape(order_index + 1));
            let dual_radius_star = dual_shape.variable(0, radius_star.clone());
            let dual_center = center
                .iter()
                .map(|momentum| momentum.map_ref(&|value| new_constant(&dual_radius_star, value)))
                .collect::<LoopMomenta<_>>();
            let dual_loop_momenta = unit_shifted_momenta
                .rescale_with_hyper_dual(&dual_radius_star, subspace.as_subspace_simple())
                .iter()
                .zip(dual_center.iter())
                .map(|(momentum, center)| momentum.clone() + center.clone())
                .collect::<LoopMomenta<_>>();
            let dual_externals =
                dualize_external_momenta(&dual_radius_star, root_sample.external_moms());
            DualOrNot::Dual(esurface.compute_from_dual_momenta(
                common_lmb,
                real_mass_vector,
                &dual_loop_momenta,
                &dual_externals,
            ))
        };

        let sample_for_order =
            sample_for_order_in_common_lmb.lmb_transform(common_lmb, &graph.loop_momentum_basis);
        let params = T::get_parameters(
            param_builder,
            (false, false),
            graph,
            &sample_for_order,
            settings.kinematics.externals.get_helicities(),
            &settings.additional_params(),
            Some(&threshold_params),
            None,
            None,
        );
        let pass_one = evaluator_stack
            .evaluate(
                params,
                orientations,
                settings,
                evaluation_metadata,
                record_primary_timing,
            )
            .expect("Amplitude variant counterterm evaluator stack failed")
            .pop()
            .unwrap();
        let prefactor = evaluate_generalized_multichanneling_prefactor(
            &sample_for_order_in_common_lmb,
            real_mass_vector,
            all_lmbs,
            thresholds,
            threshold_subspaces,
            overlap,
            overlap_group_index,
            prefactor_power,
        );
        let pass_one = multiply_dual_or_not_complex(pass_one, &prefactor);

        let mut helper_params = Vec::new();
        match pass_one {
            DualOrNot::Dual(dual) => {
                helper_params.extend_from_slice(&extract_t_derivatives_complex(dual));
            }
            DualOrNot::NonDual(value) => helper_params.push(value),
        }
        match esurface_derivatives {
            DualOrNot::Dual(dual) => {
                helper_params.extend(
                    extract_t_derivatives(dual)[1..]
                        .iter()
                        .cloned()
                        .map(Complex::new_re),
                );
            }
            DualOrNot::NonDual(value) => helper_params.push(Complex::new_re(value)),
        }
        helper_params.push(Complex::new_re(radius.clone()));
        helper_params.push(Complex::new_re(radius_star.clone()));
        helper_params.push(Complex::new_re(uv_damp_plus.clone()));
        helper_params.push(Complex::new_re(uv_damp_minus.clone()));
        helper_params.push(Complex::new_re(h_function.clone()));

        let mut pieces = evaluate_evaluator(
            &mut helper_evaluators[order_index],
            &helper_params,
            evaluation_metadata,
            record_primary_timing,
        );
        assert_eq!(
            pieces.len(),
            2,
            "amplitude threshold helper must return local and integrated pieces",
        );
        let integrated = pieces.pop().unwrap().unwrap_real();
        let local = pieces.pop().unwrap().unwrap_real();
        evaluation.add_completed_pieces(local, integrated, local_multiplier, integrated_multiplier);
    }

    crate::debug_tags!(#integration, #subtraction, #threshold, #inspect;
        stage = "amplitude_threshold_variant_total",
        graph = %graph.name,
        esurface_id = esurface_id.0,
        raised_esurface_id = raised_esurface_id.0,
        overlap_group = overlap_group_index,
        result = %format!("{:+16e}", evaluation.weighted),
        "amplitude threshold variant total"
    );
    Ok(evaluation.finish(local_multiplier, integrated_multiplier))
}

struct CounterTermBuilder<'a, T: FloatLike> {
    overlap_structure: &'a OverlapStructure,
    raised_data: &'a RaisedEsurfaceData,
    own_group_position: GraphGroupPosition,
    esurface_map: &'a TiVec<GroupEsurfaceId, TiVec<GraphGroupPosition, Option<RaisedEsurfaceId>>>,
    real_mass_vector: EdgeVec<F<T>>,
    e_cm: F<T>,
    graph: &'a Graph,
    rotation_for_overlap: &'a Rotation,
    settings: &'a RuntimeSettings,
    esurface_collection: &'a EsurfaceCollection,
    sample: &'a MomentumSample<T>,
}

impl<'a, T: FloatLike> CounterTermBuilder<'a, T> {
    #[allow(clippy::too_many_arguments)]
    fn new(
        graph: &'a Graph,
        model: &'a Model,
        rotation_for_overlap: &'a Rotation,
        settings: &'a RuntimeSettings,
        esurface_collection: &'a EsurfaceCollection,
        sample: &'a MomentumSample<T>,
        overlap_structure: &'a OverlapStructure,
        raised_data: &'a RaisedEsurfaceData,
        own_group_position: GraphGroupPosition,
        esurface_map: &'a TiVec<
            GroupEsurfaceId,
            TiVec<GraphGroupPosition, Option<RaisedEsurfaceId>>,
        >,
    ) -> Self {
        let real_mass_vector = graph.get_real_mass_vector(model);
        let e_cm = F::from_f64(settings.kinematics.e_cm);

        Self {
            real_mass_vector,
            e_cm,
            graph,
            rotation_for_overlap,
            settings,
            esurface_collection,
            overlap_structure,
            raised_data,
            sample,
            own_group_position,
            esurface_map,
        }
    }

    fn new_overlap_builder(&'a self, overlap_group: &'a OverlapGroup) -> OverlapBuilder<'a, T> {
        let center = &overlap_group.center;

        // Overlap construction stores amplitude centers in the identity frame. Rotate the center
        // exactly once here, alongside the sample rotation used by this stability probe.
        let (unrotated_center, rotated_center) = (
            center.cast(),
            center.rotate(self.rotation_for_overlap).cast(),
        );

        let shifted_loop_momenta = self.sample.loop_moms() - &rotated_center;
        let radius = shifted_loop_momenta.hyper_radius_squared(None).sqrt();
        let unit_shifted_momenta = shifted_loop_momenta.rescale(&radius.inv(), None);

        OverlapBuilder {
            counterterm_builder: self,
            overlap_group,
            rotated_center,
            _unrotated_center: unrotated_center,
            unit_shifted_momenta,
            radius,
        }
    }
}

struct OverlapBuilder<'a, T: FloatLike> {
    counterterm_builder: &'a CounterTermBuilder<'a, T>,
    overlap_group: &'a OverlapGroup,
    /// The center after its single identity-to-probe-frame rotation.
    rotated_center: LoopMomenta<F<T>>,
    /// The stored identity-frame center, retained for frame diagnostics.
    _unrotated_center: LoopMomenta<F<T>>,
    unit_shifted_momenta: LoopMomenta<F<T>>,
    radius: F<T>,
}

impl<'a, T: FloatLike> OverlapBuilder<'a, T> {
    fn new_esurface_builder(
        &'a self,
        existing_esurface_id: ExistingEsurfaceId,
    ) -> Option<EsurfaceCTBuilder<'a, T>> {
        let group_esurface_id = self
            .counterterm_builder
            .overlap_structure
            .existing_esurfaces[existing_esurface_id];

        let raised_esurface_id = self.counterterm_builder.esurface_map[group_esurface_id]
            [self.counterterm_builder.own_group_position];

        raised_esurface_id.map(|raised_esurface_id| {
            let esurface_id = self.counterterm_builder.raised_data.raised_groups
                [raised_esurface_id]
                .esurface_ids[0];

            EsurfaceCTBuilder {
                overlap_builder: self,
                _existing_esurface_id: existing_esurface_id,
                group_esurface_id,
                esurface: &self.counterterm_builder.esurface_collection[esurface_id],
                esurface_id,
                raised_esurface_id,
            }
        })
    }
}

struct EsurfaceCTBuilder<'a, T: FloatLike> {
    overlap_builder: &'a OverlapBuilder<'a, T>,
    _existing_esurface_id: ExistingEsurfaceId,
    group_esurface_id: GroupEsurfaceId,
    esurface: &'a Esurface,
    esurface_id: EsurfaceID,
    raised_esurface_id: RaisedEsurfaceId,
}

impl<'a, T: FloatLike> EsurfaceCTBuilder<'a, T> {
    fn solve_rstar(
        self,
        radial_root_identity: &RadialRootIdentity,
        radial_root_diagnostics: &mut RadialRootDiagnostics,
    ) -> Option<RstarSolution<'a, T>> {
        let mut all_center_values_valid = true;
        let center_surface_values = self
            .overlap_builder
            .overlap_group
            .existing_esurfaces
            .iter()
            .map(|&existing_esurface_id| {
                let group_esurface_id = self
                    .overlap_builder
                    .counterterm_builder
                    .overlap_structure
                    .existing_esurfaces[existing_esurface_id];

                match self.overlap_builder.counterterm_builder.esurface_map[group_esurface_id]
                    [self.overlap_builder.counterterm_builder.own_group_position]
                {
                    Some(raised_esurface_id) => {
                        let esurface_id = self
                            .overlap_builder
                            .counterterm_builder
                            .raised_data
                            .raised_groups[raised_esurface_id]
                            .esurface_ids[0];
                        let esurface =
                            &self.overlap_builder.counterterm_builder.esurface_collection
                                [esurface_id];
                        let value = esurface.compute_from_momenta(
                            &self
                                .overlap_builder
                                .counterterm_builder
                                .graph
                                .loop_momentum_basis,
                            &self.overlap_builder.counterterm_builder.real_mass_vector,
                            &self.overlap_builder.rotated_center,
                            self.overlap_builder
                                .counterterm_builder
                                .sample
                                .external_moms(),
                        );
                        let is_valid = esurface_value_is_strictly_inside(
                            &value,
                            &self.overlap_builder.counterterm_builder.e_cm,
                        );
                        all_center_values_valid &= is_valid;
                        format!(
                            "existing={} group={} raised={} local={} edges={:?} value={:+16e} inside={}",
                            usize::from(existing_esurface_id),
                            group_esurface_id.0,
                            raised_esurface_id.0,
                            esurface_id.0,
                            esurface.energies,
                            value,
                            is_valid
                        )
                    }
                    None => {
                        format!(
                            "existing={} group={} absent_in_current_graph",
                            usize::from(existing_esurface_id),
                            group_esurface_id.0
                        )
                    }
                }
            })
            .collect::<Vec<_>>();
        crate::debug_tags!(#integration, #subtraction, #threshold, #inspect, #center;
            stage = "amplitude_threshold_center_values",
            graph = %self.overlap_builder.counterterm_builder.graph.name,
            selected_existing_esurface_id = usize::from(self._existing_esurface_id),
            selected_group_esurface_id = self.group_esurface_id.0,
            selected_raised_esurface_id = self.raised_esurface_id.0,
            selected_esurface_id = self.esurface_id.0,
            rotation_id = %self.overlap_builder.counterterm_builder.rotation_for_overlap.method,
            center_provenance = "identity_frame_rotated_once",
            overlap_group_size = self.overlap_builder.overlap_group.existing_esurfaces.len(),
            radius = %format!("{:+16e}", self.overlap_builder.radius),
            file.rotated_center = %format!("{}", self.overlap_builder.rotated_center),
            file.surface_values = %center_surface_values.join("\n"),
            "amplitude threshold center values"
        );

        if !all_center_values_valid {
            warn!(
                graph = %self.overlap_builder.counterterm_builder.graph.name,
                selected_esurface_id = self.esurface_id.0,
                rotation_id = %self.overlap_builder.counterterm_builder.rotation_for_overlap.method,
                center_provenance = "identity_frame_rotated_once",
                center = %self.overlap_builder.rotated_center,
                surface_values = %center_surface_values.join("; "),
                "refusing to evaluate an amplitude threshold counterterm with an invalid probe-frame overlap center"
            );
            return None;
        }

        let (raw_radius_guess, _) = self.esurface.get_radius_guess(
            &self.overlap_builder.unit_shifted_momenta,
            self.overlap_builder
                .counterterm_builder
                .sample
                .external_moms(),
            &self
                .overlap_builder
                .counterterm_builder
                .graph
                .loop_momentum_basis,
        );

        let function = |r: &_| {
            self.esurface.compute_self_and_r_derivative(
                r,
                &self.overlap_builder.unit_shifted_momenta,
                &self.overlap_builder.rotated_center,
                self.overlap_builder
                    .counterterm_builder
                    .sample
                    .external_moms(),
                &self.overlap_builder.counterterm_builder.real_mass_vector,
                &self
                    .overlap_builder
                    .counterterm_builder
                    .graph
                    .loop_momentum_basis,
            )
        };

        let zero = raw_radius_guess.zero();
        let mut radius_guess = raw_radius_guess.clone();
        if radius_guess.is_nan() || radius_guess.is_infinite() || radius_guess <= zero {
            radius_guess = self.overlap_builder.counterterm_builder.e_cm.clone();
        }
        let tolerance = F::from_f64(
            self.overlap_builder
                .counterterm_builder
                .settings
                .subtraction
                .radial_root_residual_tolerance,
        );
        let solution = match radial_root_diagnostics.solve(
            radial_root_identity,
            &zero,
            &radius_guess,
            function,
            &tolerance,
            MAX_ITERATIONS,
            64,
            &self.overlap_builder.counterterm_builder.e_cm,
        ) {
            Ok(solution) => solution,
            Err(error) => {
                warn!(
                    graph = %self.overlap_builder.counterterm_builder.graph.name,
                    esurface_id = self.esurface_id.0,
                    rotation_id = %self.overlap_builder.counterterm_builder.rotation_for_overlap.method,
                    center_provenance = "identity_frame_rotated_once",
                    raw_radius_guess = %raw_radius_guess,
                    radius_guess = %radius_guess,
                    error = %error,
                    "refusing to evaluate an amplitude threshold counterterm with an invalid radial solution"
                );
                return None;
            }
        };
        crate::debug_tags!(#integration, #subtraction, #threshold, #inspect;
            stage = "amplitude_threshold_rstar_solution",
            graph = %self.overlap_builder.counterterm_builder.graph.name,
            existing_esurface_id = %self._existing_esurface_id,
            group_esurface_id = self.group_esurface_id.0,
            raised_esurface_id = self.raised_esurface_id.0,
            esurface_id = self.esurface_id.0,
            rotation_id = %self.overlap_builder.counterterm_builder.rotation_for_overlap.method,
            center_provenance = "identity_frame_rotated_once",
            radius_guess = %format!("{:+16e}", radius_guess),
            radius_star = %format!("{:+16e}", solution.solution),
            derivative = %format!("{:+16e}", solution.derivative_at_solution),
            error = %format!("{:+16e}", solution.error_of_function),
            iterations = solution.num_iterations_used,
            nonfinite = false,
            "amplitude threshold rstar solution"
        );

        Some(RstarSolution {
            esurface_ct_builder: self,
            solution,
        })
    }
}

struct RstarSolution<'a, T: FloatLike> {
    esurface_ct_builder: EsurfaceCTBuilder<'a, T>,
    solution: NewtonIterationResult<T>,
}

impl<'a, T: FloatLike> RstarSolution<'a, T> {
    fn rstar_samples(self) -> RstarSample<'a, T> {
        let rstar_loop_momenta = &self
            .esurface_ct_builder
            .overlap_builder
            .unit_shifted_momenta
            .rescale(&self.solution.solution, None)
            + &self.esurface_ct_builder.overlap_builder.rotated_center;

        let mut rstar_sample = self
            .esurface_ct_builder
            .overlap_builder
            .counterterm_builder
            .sample
            .clone();

        rstar_sample.sample.loop_moms = rstar_loop_momenta;

        RstarSample {
            rstar_solution: self,
            rstar_sample,
        }
    }
}

struct RstarSample<'a, T: FloatLike> {
    rstar_solution: RstarSolution<'a, T>,
    rstar_sample: MomentumSample<T>,
}

impl<'a, T: FloatLike> RstarSample<'a, T> {
    #[allow(clippy::too_many_arguments)]
    fn evaluate<'b, 'c: 'b>(
        self,
        param_builder: &mut ParamBuilder<f64>,
        orientations: SingleOrAllOrientations<'a, OrientationID>,
        evaluation_metadata: &mut EvaluationMetaData,
        record_primary_timing: bool,
        ct_evaluator: &mut AmplitudeCountertermEvaluator,
        helper_evaluators: &mut [GenericEvaluator],
        split_helpers: bool,
        record_components: bool,
    ) -> Result<EvaluatedAmplitudeThreshold<T>> {
        assert!(
            split_helpers || !record_components,
            "recorded legacy amplitude components require split helpers",
        );
        let esurface_ct_builder = &self.rstar_solution.esurface_ct_builder;
        let ct_builder = esurface_ct_builder.overlap_builder.counterterm_builder;

        let esurface_id = esurface_ct_builder.esurface_id;

        let model_params = param_builder
            .model_values()
            .iter()
            .map(|c| Complex::new(F::from_ff64(c.re), F::from_ff64(c.im)))
            .collect::<Vec<_>>();

        let radius = self
            .rstar_solution
            .esurface_ct_builder
            .overlap_builder
            .radius
            .clone();

        let radius_star = self.rstar_solution.solution.solution.clone();
        let e_cm = &ct_builder.e_cm;
        let settings = &ct_builder
            .settings
            .subtraction
            .local_ct_settings
            .uv_localisation;

        let integrated_settings = &ct_builder.settings.subtraction.integrated_ct_settings;

        let uv_damp_plus = evaluate_uv_damper(&radius, &radius_star, e_cm, settings);
        let uv_damp_minus = evaluate_uv_damper(&-&radius, &radius_star, e_cm, settings);

        debug!("uv_damp_plus: {:?}", uv_damp_plus);
        debug!("uv_damp_minus: {:?}", uv_damp_minus);
        debug!("radius: {:?}", radius);
        debug!("radius_star: {:?}", radius_star);

        let h_function =
            evaluate_integrated_ct_normalisation(&radius, &radius_star, e_cm, integrated_settings);
        crate::debug_tags!(#integration, #subtraction, #threshold, #inspect;
            stage = "amplitude_threshold_prefactors",
            graph = %ct_builder.graph.name,
            esurface_id = esurface_id.0,
            raised_esurface_id = esurface_ct_builder.raised_esurface_id.0,
            radius = %format!("{:+16e}", radius),
            radius_star = %format!("{:+16e}", radius_star),
            derivative = %format!(
                "{:+16e}",
                self.rstar_solution.solution.derivative_at_solution
            ),
            uv_damp_plus = %format!("{:+16e}", uv_damp_plus),
            uv_damp_minus = %format!("{:+16e}", uv_damp_minus),
            h_function = %format!("{:+16e}", h_function),
            "amplitude threshold prefactors"
        );

        let debug_diagnostics_enabled = tracing::event_enabled!(tracing::Level::DEBUG);
        if debug_diagnostics_enabled {
            let coincidence_tolerance = F::from_f64(1.0e-8) * e_cm;
            for (candidate_esurface_id, candidate_esurface) in
                ct_builder.esurface_collection.iter_enumerated()
            {
                let value = candidate_esurface.compute_from_momenta(
                    &ct_builder.graph.loop_momentum_basis,
                    &ct_builder.real_mass_vector,
                    self.rstar_sample.loop_moms(),
                    self.rstar_sample.external_moms(),
                );
                if value.abs() < coincidence_tolerance {
                    let candidate_raised_esurface_id = ct_builder
                        .raised_data
                        .raised_groups
                        .iter_enumerated()
                        .find_map(|(raised_esurface_id, raised_group)| {
                            raised_group
                                .esurface_ids
                                .contains(&candidate_esurface_id)
                                .then_some(raised_esurface_id.0)
                        });
                    let candidate_atom = candidate_esurface.to_atom(&[]).to_string();
                    crate::debug_tags!(#integration, #subtraction, #threshold, #inspect, #esurface;
                        stage = "amplitude_threshold_rstar_coincident_esurface",
                        graph = %ct_builder.graph.name,
                        selected_esurface_id = esurface_id.0,
                        selected_raised_esurface_id = esurface_ct_builder.raised_esurface_id.0,
                        candidate_esurface_id = candidate_esurface_id.0,
                        candidate_raised_esurface_id = ?candidate_raised_esurface_id,
                        selected = candidate_esurface_id == esurface_id,
                        value = %format!("{:+16e}", value),
                        tolerance = %format!("{:+16e}", coincidence_tolerance),
                        file.atom = %candidate_atom,
                        "threshold rstar coincident esurface graph={} selected={} raised={} candidate={} candidate_raised={:?} value={} tol={} atom={}",
                        ct_builder.graph.name,
                        esurface_id.0,
                        esurface_ct_builder.raised_esurface_id.0,
                        candidate_esurface_id.0,
                        candidate_raised_esurface_id,
                        format!("{:+16e}", value),
                        format!("{:+16e}", coincidence_tolerance),
                        candidate_atom
                    );
                }
            }
        }

        let mut evaluation =
            EvaluatedAmplitudeThreshold::new(self.rstar_sample.zero(), record_components);

        for (cut_cff_index, evaluator_stack) in ct_evaluator.evaluator_stacks.iter_mut() {
            let order_index = cut_cff_index.left_threshold_order.unwrap() - 1;
            let (sample_for_order, threshold_params) = if order_index == 0 {
                debug!(
                    "rescaled loop momenta at rstar:\n{}",
                    self.rstar_sample.loop_moms()
                );
                (
                    self.rstar_sample.clone(),
                    ThresholdParams {
                        radius: DualOrNot::NonDual(radius.clone()),
                        radius_star: DualOrNot::NonDual(radius_star.clone()),
                        esurface_derivative: DualOrNot::NonDual(
                            self.rstar_solution.solution.derivative_at_solution.clone(),
                        ),
                        uv_damp_plus: DualOrNot::NonDual(uv_damp_plus.clone()),
                        uv_damp_minus: DualOrNot::NonDual(uv_damp_minus.clone()),
                        h_function: DualOrNot::NonDual(h_function.clone()),
                    },
                )
            } else {
                let dual_shape = HyperDual::<F<T>>::new(simple_n_deriv_shape(order_index));
                let dual_radius_star = dual_shape.variable(0, radius_star.clone());

                let dualized_center = self
                    .rstar_solution
                    .esurface_ct_builder
                    .overlap_builder
                    .rotated_center
                    .iter()
                    .map(|momentum| {
                        momentum.map_ref(&|value| new_constant(&dual_radius_star, value))
                    })
                    .collect::<LoopMomenta<_>>();
                let dual_loop_momenta = self
                    .rstar_solution
                    .esurface_ct_builder
                    .overlap_builder
                    .unit_shifted_momenta
                    .rescale_with_hyper_dual(&dual_radius_star, None)
                    .iter()
                    .zip(dualized_center.iter())
                    .map(|(momentum, center)| momentum.clone() + center.clone())
                    .collect::<LoopMomenta<_>>();
                debug!("rescaled loop momenta at rstar:\n{}", dual_loop_momenta);
                let mut sample_with_duals = self.rstar_sample.clone();
                sample_with_duals.sample.dual_loop_moms = Some(dual_loop_momenta);

                (
                    sample_with_duals,
                    ThresholdParams {
                        radius: DualOrNot::Dual(new_constant(&dual_radius_star, &radius)),
                        radius_star: DualOrNot::Dual(dual_radius_star.clone()),
                        esurface_derivative: DualOrNot::Dual(new_constant(
                            &dual_radius_star,
                            &self.rstar_solution.solution.derivative_at_solution.clone(),
                        )),
                        uv_damp_plus: DualOrNot::Dual(new_constant(
                            &dual_radius_star,
                            &uv_damp_plus,
                        )),
                        uv_damp_minus: DualOrNot::Dual(new_constant(
                            &dual_radius_star,
                            &uv_damp_minus,
                        )),
                        h_function: DualOrNot::Dual(new_constant(&dual_radius_star, &h_function)),
                    },
                )
            };

            let esurface_derivatives = if order_index == 0 {
                DualOrNot::NonDual(self.rstar_solution.solution.derivative_at_solution.clone())
            } else {
                let dual_shape_for_esurface =
                    HyperDual::<F<T>>::new(simple_n_deriv_shape(order_index + 1));
                let dual_radius_star_for_esurface =
                    dual_shape_for_esurface.variable(0, radius_star.clone());
                let dualized_center_for_esurface = self
                    .rstar_solution
                    .esurface_ct_builder
                    .overlap_builder
                    .rotated_center
                    .iter()
                    .map(|momentum| {
                        momentum
                            .map_ref(&|value| new_constant(&dual_radius_star_for_esurface, value))
                    })
                    .collect::<LoopMomenta<_>>();
                let dual_loop_momenta_for_esurface = self
                    .rstar_solution
                    .esurface_ct_builder
                    .overlap_builder
                    .unit_shifted_momenta
                    .rescale_with_hyper_dual(&dual_radius_star_for_esurface, None)
                    .iter()
                    .zip(dualized_center_for_esurface.iter())
                    .map(|(momentum, center)| momentum.clone() + center.clone())
                    .collect::<LoopMomenta<_>>();

                let dualized_externals = self
                    .rstar_sample
                    .external_moms()
                    .iter()
                    .map(|momentum| FourMomentum {
                        temporal: Energy {
                            value: new_constant(
                                &dual_radius_star_for_esurface,
                                &momentum.temporal.value,
                            ),
                        },
                        spatial: momentum
                            .spatial
                            .map_ref(&|value| new_constant(&dual_radius_star_for_esurface, value)),
                    })
                    .collect();

                DualOrNot::Dual(
                    self.rstar_solution
                        .esurface_ct_builder
                        .esurface
                        .compute_from_dual_momenta(
                            &ct_builder.graph.loop_momentum_basis,
                            &ct_builder.real_mass_vector,
                            &dual_loop_momenta_for_esurface,
                            &dualized_externals,
                        ),
                )
            };

            let params = T::get_parameters(
                param_builder,
                (false, false),
                ct_builder.graph,
                &sample_for_order,
                ct_builder.settings.kinematics.externals.get_helicities(),
                &ct_builder.settings.additional_params(),
                Some(&threshold_params),
                None,
                None,
            );
            let params_slice = params.as_slice();
            let params_nonfinite_count = params_slice
                .iter()
                .filter(|value| {
                    value.re.is_nan()
                        || value.re.is_infinite()
                        || value.im.is_nan()
                        || value.im.is_infinite()
                })
                .count();
            let first_params_nonfinite_index = params_slice.iter().position(|value| {
                value.re.is_nan()
                    || value.re.is_infinite()
                    || value.im.is_nan()
                    || value.im.is_infinite()
            });
            crate::debug_tags!(#integration, #subtraction, #threshold, #inspect;
                stage = "amplitude_threshold_pass_one_params",
                graph = %ct_builder.graph.name,
                esurface_id = esurface_id.0,
                raised_esurface_id = esurface_ct_builder.raised_esurface_id.0,
                order = order_index + 1,
                params_len = params_slice.len(),
                params_nonfinite_count,
                first_params_nonfinite_index = ?first_params_nonfinite_index,
                "amplitude threshold pass one params"
            );

            let pass_one_result = evaluator_stack
                .evaluate(
                    params,
                    orientations,
                    ct_builder.settings,
                    evaluation_metadata,
                    record_primary_timing,
                )
                .expect("Amplitude counterterm evaluator stack failed")
                .pop()
                .unwrap();

            let prefactor = self.evaluate_multichanneling_prefactor(
                &sample_for_order,
                &model_params,
                evaluation_metadata,
                record_primary_timing,
                order_index,
            );

            let raw_pass_one_is_nonfinite = match &pass_one_result {
                DualOrNot::Dual(dual_result) => dual_result.values.iter().any(|value| {
                    value.re.is_nan()
                        || value.re.is_infinite()
                        || value.im.is_nan()
                        || value.im.is_infinite()
                }),
                DualOrNot::NonDual(non_dual_result) => {
                    non_dual_result.re.is_nan()
                        || non_dual_result.re.is_infinite()
                        || non_dual_result.im.is_nan()
                        || non_dual_result.im.is_infinite()
                }
            };
            let prefactor_is_nonfinite = match &prefactor {
                DualOrNot::Dual(dual_result) => dual_result.values.iter().any(|value| {
                    value.re.is_nan()
                        || value.re.is_infinite()
                        || value.im.is_nan()
                        || value.im.is_infinite()
                }),
                DualOrNot::NonDual(non_dual_result) => {
                    non_dual_result.re.is_nan()
                        || non_dual_result.re.is_infinite()
                        || non_dual_result.im.is_nan()
                        || non_dual_result.im.is_infinite()
                }
            };
            crate::debug_tags!(#integration, #subtraction, #threshold, #inspect;
                stage = "amplitude_threshold_pass_one_raw",
                graph = %ct_builder.graph.name,
                esurface_id = esurface_id.0,
                raised_esurface_id = esurface_ct_builder.raised_esurface_id.0,
                order = order_index + 1,
                raw_result = %format!("{pass_one_result}"),
                prefactor = %format!("{prefactor}"),
                raw_nonfinite = raw_pass_one_is_nonfinite,
                prefactor_nonfinite = prefactor_is_nonfinite,
                "amplitude threshold pass one raw"
            );

            let pass_one_result = multiply_dual_or_not_complex(pass_one_result, &prefactor);
            let pass_one_is_nonfinite = match &pass_one_result {
                DualOrNot::Dual(dual_result) => dual_result.values.iter().any(|value| {
                    value.re.is_nan()
                        || value.re.is_infinite()
                        || value.im.is_nan()
                        || value.im.is_infinite()
                }),
                DualOrNot::NonDual(non_dual_result) => {
                    non_dual_result.re.is_nan()
                        || non_dual_result.re.is_infinite()
                        || non_dual_result.im.is_nan()
                        || non_dual_result.im.is_infinite()
                }
            };
            crate::debug_tags!(#integration, #subtraction, #threshold, #inspect;
                stage = "amplitude_threshold_pass_one",
                graph = %ct_builder.graph.name,
                esurface_id = esurface_id.0,
                raised_esurface_id = esurface_ct_builder.raised_esurface_id.0,
                order = order_index + 1,
                prefactor = %format!("{prefactor}"),
                result = %format!("{pass_one_result}"),
                nonfinite = pass_one_is_nonfinite,
                "amplitude threshold pass one"
            );

            debug!(
                "Pass one result for esurface {} order {}: {}",
                esurface_id.0,
                order_index + 1,
                pass_one_result
            );

            let mut params_for_pass_two = vec![];
            match pass_one_result {
                DualOrNot::Dual(dual_result) => {
                    params_for_pass_two
                        .extend_from_slice(&extract_t_derivatives_complex(dual_result));
                }
                DualOrNot::NonDual(non_dual_result) => {
                    params_for_pass_two.push(non_dual_result);
                }
            }

            match esurface_derivatives {
                DualOrNot::Dual(dual_e_surface) => {
                    extract_t_derivatives(dual_e_surface)[1..]
                        .iter()
                        .for_each(|value| params_for_pass_two.push(Complex::new_re(value.clone())));
                }
                DualOrNot::NonDual(non_dual_e_surface) => {
                    params_for_pass_two.push(Complex::new_re(non_dual_e_surface));
                }
            }

            params_for_pass_two.push(Complex::new_re(radius.clone()));
            params_for_pass_two.push(Complex::new_re(radius_star.clone()));
            params_for_pass_two.push(Complex::new_re(uv_damp_plus.clone()));
            params_for_pass_two.push(Complex::new_re(uv_damp_minus.clone()));
            params_for_pass_two.push(Complex::new_re(h_function.clone()));

            let pass_two_result = if split_helpers {
                let mut pieces = evaluate_evaluator(
                    &mut helper_evaluators[order_index],
                    &params_for_pass_two,
                    evaluation_metadata,
                    record_primary_timing,
                );
                assert_eq!(
                    pieces.len(),
                    3,
                    "recording amplitude threshold helper must return legacy, local, and integrated pieces",
                );
                let integrated = pieces.pop().unwrap().unwrap_real();
                let local = pieces.pop().unwrap().unwrap_real();
                let legacy = pieces.pop().unwrap().unwrap_real();
                if record_components {
                    let one = self.rstar_sample.one();
                    let combined = combine_completed_threshold_pieces(
                        local.clone(),
                        integrated.clone(),
                        &one,
                        &one,
                    );
                    evaluation.add_completed_pieces(local, integrated, &one, &one);
                    combined
                } else {
                    evaluation.weighted += legacy.clone();
                    legacy
                }
            } else {
                let combined = evaluate_evaluator_single(
                    &mut helper_evaluators[order_index],
                    &params_for_pass_two,
                    evaluation_metadata,
                    record_primary_timing,
                );
                evaluation.weighted += combined.clone();
                combined
            };

            debug!(
                "Pass two result for esurface {} order {}: {}",
                esurface_id.0,
                order_index + 1,
                pass_two_result
            );
            let pass_two_is_nonfinite = pass_two_result.re.is_nan()
                || pass_two_result.re.is_infinite()
                || pass_two_result.im.is_nan()
                || pass_two_result.im.is_infinite();
            crate::debug_tags!(#integration, #subtraction, #threshold, #inspect;
                stage = "amplitude_threshold_pass_two",
                graph = %ct_builder.graph.name,
                esurface_id = esurface_id.0,
                raised_esurface_id = esurface_ct_builder.raised_esurface_id.0,
                order = order_index + 1,
                result = %format!("{:+16e}", pass_two_result),
                nonfinite = pass_two_is_nonfinite,
                "amplitude threshold pass two"
            );
        }

        debug!(
            ct_eval = format!("{:+16e}", evaluation.weighted),
            "esurface {}", esurface_id.0
        );
        let total_ct_is_nonfinite = evaluation.weighted.re.is_nan()
            || evaluation.weighted.re.is_infinite()
            || evaluation.weighted.im.is_nan()
            || evaluation.weighted.im.is_infinite();
        crate::debug_tags!(#integration, #subtraction, #threshold, #inspect;
            stage = "amplitude_threshold_total",
            graph = %ct_builder.graph.name,
            esurface_id = esurface_id.0,
            raised_esurface_id = esurface_ct_builder.raised_esurface_id.0,
            result = %format!("{:+16e}", evaluation.weighted),
            nonfinite = total_ct_is_nonfinite,
            "amplitude threshold total"
        );

        let one = self.rstar_sample.one();
        Ok(evaluation.finish(&one, &one))
    }

    fn evaluate_multichanneling_prefactor(
        &self,
        momentum_sample: &MomentumSample<T>,
        model_params: &[Complex<F<T>>],
        evaluation_metadata: &mut EvaluationMetaData,
        record_primary_timing: bool,
        order_index: usize,
    ) -> DualOrNot<Complex<F<T>>> {
        let overlap_builder = self.rstar_solution.esurface_ct_builder.overlap_builder;
        let overlap = overlap_builder.counterterm_builder.overlap_structure;

        if overlap.overlap_groups.len() < 2 {
            return DualOrNot::NonDual(Complex::new_re(momentum_sample.one()));
        }

        let multiplicative_offset = momentum_sample
            .sample
            .dual_loop_moms
            .as_ref()
            .map(|dual_loop_moms| dual_loop_moms.first().unwrap().px.values.len())
            .unwrap_or(1);
        let zero = Complex::new_re(momentum_sample.zero());
        let mut params = if let Some(dual_loop_moms) = &momentum_sample.sample.dual_loop_moms {
            dual_loop_moms
                .iter()
                .flat_map(|mom| {
                    [
                        mom.px.values.clone(),
                        mom.py.values.clone(),
                        mom.pz.values.clone(),
                    ]
                    .into_iter()
                    .flatten()
                    .map(Complex::new_re)
                    .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        } else {
            momentum_sample
                .loop_moms()
                .iter()
                .flat_map(|momentum| {
                    [
                        momentum.px.clone(),
                        momentum.py.clone(),
                        momentum.pz.clone(),
                    ]
                    .into_iter()
                    .map(Complex::new_re)
                    .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        };

        params.extend(momentum_sample.external_moms().iter().flat_map(|momentum| {
            [
                momentum.temporal.value.clone(),
                momentum.spatial.px.clone(),
                momentum.spatial.py.clone(),
                momentum.spatial.pz.clone(),
            ]
            .into_iter()
            .flat_map(|value| {
                std::iter::once(Complex::new_re(value))
                    .chain((1..multiplicative_offset).map(|_| zero.clone()))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>()
        }));
        params.extend(model_params.iter().cloned().flat_map(|value| {
            std::iter::once(value)
                .chain((1..multiplicative_offset).map(|_| zero.clone()))
                .collect::<Vec<_>>()
        }));

        let evaluator = overlap_builder
            .overlap_group
            .prefactor_evaluator
            .as_ref()
            .unwrap()
            .get(order_index)
            .expect("missing overlap prefactor evaluator for amplitude threshold order");

        evaluate_evaluator(
            &mut evaluator.borrow_mut(),
            &params,
            evaluation_metadata,
            record_primary_timing,
        )
        .pop()
        .expect("overlap prefactor evaluator should return exactly one value")
    }
}

#[cfg(test)]
mod tests {
    use super::{
        AmplitudeCountertermAtom, AmplitudeCountertermData, EvaluatedAmplitudeThreshold,
        append_amplitude_component_evaluations, combine_completed_threshold_pieces,
    };
    use crate::{
        cff::{CutCFFIndex, esurface::RaisedEsurfaceId},
        dot,
        graph::{Graph, GraphGroupPosition, parse::from_dot::IntoGraph},
        initialisation::test_initialise,
        integrands::{
            evaluation::EvaluationMetaData,
            process::threshold_multiplier::{
                ThresholdMultiplierEvaluatorCollection, ThresholdMultiplierInput,
                ThresholdMultiplierInputWorkspace, ThresholdMultiplierLayout,
                ThresholdMultiplierPoint,
            },
        },
        momentum::{
            FourMomentum, ThreeMomentum,
            sample::{BareMomentumSample, ExternalFourMomenta, LoopMomenta, MomentumSample},
        },
        processes::{EvaluatorSettings, ThresholdCountertermVariantId},
        utils::F,
        uv::Integrands,
    };
    use linnet::half_edge::involution::EdgeVec;
    use spenso::algebra::complex::Complex;
    use symbolica::{atom::Atom, symbol};
    use typed_index_collections::ti_vec;

    #[test]
    fn empty_amplitude_counterterm_atom_is_not_generated() {
        let atom = AmplitudeCountertermAtom {
            parametric: std::iter::empty().collect(),
        };

        assert!(!atom.is_generated());
    }

    #[test]
    fn non_empty_amplitude_counterterm_atom_is_generated() {
        let atom = AmplitudeCountertermAtom {
            parametric: Integrands::from_iter([(
                CutCFFIndex::new_all_none(),
                Atom::var(symbol!("x")),
            )]),
        };

        assert!(atom.is_generated());
    }

    #[test]
    fn zero_like_amplitude_counterterm_atom_is_zero() {
        let atom = AmplitudeCountertermAtom {
            parametric: Integrands::from_iter([(
                CutCFFIndex::new_all_none(),
                Atom::var(symbol!("x")),
            )]),
        };

        let zeroed = atom.zero_like();

        assert_eq!(
            zeroed.parametric,
            Integrands::from_iter([(CutCFFIndex::new_all_none(), Atom::Zero)])
        );
    }

    #[test]
    fn inactive_amplitude_esurface_guard_reports_runtime_access() {
        let mut data = AmplitudeCountertermData::new_empty(GraphGroupPosition(0));
        data.active_mask = ti_vec![false];

        let error = data
            .ensure_active_raised_esurface(RaisedEsurfaceId(0))
            .unwrap_err();
        assert!(error.to_string().contains("generation marked it inactive"));
    }

    #[test]
    fn generalized_multiplier_uses_component_contexts_after_helper_evaluation() {
        test_initialise().unwrap();
        let graph: Graph = dot!(digraph amplitude_multiplier_contexts {
            ext [style=invis]
            v1 [num=1]
            v2 [num=1]
            edge [num=1 mass=0]
            ext -> v1 [id=0]
            v1 -> v2 [id=1]
            v1 -> v2 [id=2]
            v2 -> ext [id=3]
        })
        .unwrap();
        let edge = graph.loop_momentum_basis.loop_edges.raw[0];
        let layout = ThresholdMultiplierLayout::new(
            Vec::new(),
            Vec::new(),
            graph.loop_momentum_basis.ext_edges.len(),
            (0..graph.n_edges()).collect(),
            Vec::new(),
        )
        .unwrap();
        let expression = layout
            .parse_expression(&format!(
                "Q3(effective, {}, cind(1)) + 10 * Q3(star, {}, cind(1))",
                edge.0, edge.0,
            ))
            .unwrap();
        let mut collection = Some(
            ThresholdMultiplierEvaluatorCollection::build(
                layout,
                vec![(ThresholdCountertermVariantId(0), Some(expression))],
                Vec::new(),
                &EvaluatorSettings::default(),
            )
            .unwrap()
            .unwrap(),
        );

        let make_sample = |loop_scale: f64| MomentumSample {
            sample: BareMomentumSample {
                loop_moms: (0..graph.loop_momentum_basis.loop_edges.len())
                    .map(|index| ThreeMomentum::new(F(loop_scale + index as f64), F(0.25), F(-0.5)))
                    .collect::<LoopMomenta<_>>(),
                dual_loop_moms: None,
                loop_mom_cache_id: 0,
                loop_mom_base_cache_id: 0,
                external_moms: (0..graph.loop_momentum_basis.ext_edges.len())
                    .map(|position| {
                        FourMomentum::from_args(
                            F(5.0 + position as f64),
                            F(0.1 * position as f64),
                            F(0.0),
                            F(0.0),
                        )
                    })
                    .collect::<ExternalFourMomenta<_>>(),
                external_mom_cache_id: 0,
                external_mom_base_cache_id: 0,
                jacobian: F(1.0),
                orientation: None,
                parameterization_branch: None,
            },
        };
        let base = make_sample(1.0);
        let root = make_sample(2.0);
        let masses = (0..graph.n_edges()).map(|_| F(0.0)).collect::<EdgeVec<_>>();
        let mut expected_workspace =
            ThresholdMultiplierInputWorkspace::new(collection.as_ref().unwrap().layout(), F(0.0));
        let expected_for = |values: &[Complex<F<f64>>]| {
            collection
                .as_ref()
                .unwrap()
                .layout()
                .inputs()
                .iter()
                .zip(values)
                .filter_map(|(input, value)| match input {
                    ThresholdMultiplierInput::EdgeMomentum {
                        point,
                        edge: input_edge,
                        component: 1,
                    } if *input_edge == edge.0 => Some(match point {
                        ThresholdMultiplierPoint::Effective => value.re,
                        ThresholdMultiplierPoint::Star => value.re * F(10.0),
                    }),
                    _ => None,
                })
                .fold(F(0.0), |sum, value| sum + value)
        };
        let expected_local = expected_for(
            expected_workspace
                .bind(
                    collection.as_ref().unwrap().layout(),
                    &graph.loop_momentum_basis,
                    &masses,
                    &[],
                    &[],
                    &base,
                    &root,
                )
                .unwrap()
                .as_slice(),
        );
        let expected_integrated = expected_for(
            expected_workspace
                .bind(
                    collection.as_ref().unwrap().layout(),
                    &graph.loop_momentum_basis,
                    &masses,
                    &[],
                    &[],
                    &root,
                    &root,
                )
                .unwrap()
                .as_slice(),
        );
        let mut workspace =
            ThresholdMultiplierInputWorkspace::new(collection.as_ref().unwrap().layout(), F(0.0));
        let (local, integrated) = AmplitudeCountertermData::evaluate_variant_multiplier(
            &mut collection,
            Some(&mut workspace),
            ThresholdCountertermVariantId(0),
            &graph,
            &masses,
            &[],
            &[],
            &base,
            &root,
            &mut EvaluationMetaData::new_empty(),
            false,
        )
        .unwrap();
        assert_eq!(local, expected_local);
        assert_eq!(integrated, expected_integrated);
        assert_ne!(local, integrated);

        let completed = combine_completed_threshold_pieces(
            Complex::new_re(F(2.0)),
            Complex::new_re(F(3.0)),
            &F(5.0),
            &F(7.0),
        );
        assert_eq!(completed, Complex::new_re(F(31.0)));

        let mut data = AmplitudeCountertermData::new_empty(GraphGroupPosition(0));
        data.threshold_multipliers = collection;
        let mut visited = 0;
        data.for_each_generic_evaluator_mut(|_| {
            visited += 1;
            Ok(())
        })
        .unwrap();
        assert_eq!(visited, 1);
        assert_eq!(data.generic_evaluator_count(), 0);
    }

    #[test]
    fn exact_zero_multiplier_skips_a_nan_completed_piece() {
        let completed = combine_completed_threshold_pieces(
            Complex::new_re(F(f64::NAN)),
            Complex::new_re(F(3.0)),
            &F(0.0),
            &F(7.0),
        );

        assert_eq!(completed, Complex::new_re(F(21.0)));
    }

    #[test]
    fn partial_zero_records_a_skipped_component_without_its_nan_bare_value() {
        let local_multiplier = F(0.0);
        let integrated_multiplier = F(7.0);
        let mut evaluation = EvaluatedAmplitudeThreshold::new(F(0.0), true);
        evaluation.add_completed_pieces(
            Complex::new_re(F(f64::NAN)),
            Complex::new_re(F(3.0)),
            &local_multiplier,
            &integrated_multiplier,
        );
        let evaluation = evaluation.finish(&local_multiplier, &integrated_multiplier);
        let mut components = Some(Vec::new());
        append_amplitude_component_evaluations(
            &mut components,
            &evaluation,
            ThresholdCountertermVariantId(0),
            RaisedEsurfaceId(2),
            3,
            &local_multiplier,
            &integrated_multiplier,
        );
        let components = components.unwrap();

        assert_eq!(evaluation.weighted, Complex::new_re(F(21.0)));
        assert_eq!(components.len(), 2);
        assert!(components[0].evaluation_skipped);
        assert!(components[0].bare.is_none());
        assert_eq!(components[0].weighted, Complex::new_re(F(0.0)));
        assert!(!components[1].evaluation_skipped);
        assert_eq!(components[1].bare, Some(Complex::new_re(F(3.0))));
        assert_eq!(components[1].weighted, Complex::new_re(F(21.0)));
    }
}
