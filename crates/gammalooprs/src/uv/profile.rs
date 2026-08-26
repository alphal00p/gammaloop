//! UV Profile Analysis Module
//!
//! This module provides functionality for analyzing ultraviolet behavior of loop integrands
//! by evaluating them at different momentum scalings and computing degrees of divergence.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Display;
use std::path::Path;
use std::sync::{Arc, Mutex};

use crate::DependentMomentaConstructor;
use crate::cff::expression::OrientationData;
use crate::cff::orientations::GraphOrientation;
use crate::graph::parse::string_utils::ToOrderedSimple;
use crate::graph::{Graph, LmbIndex, LoopMomentumBasis};
use crate::integrands::evaluation::EvaluationResult;
use crate::integrands::process::{
    OrientationProfileMode, ProcessIntegrand, evaluate_profile_momentum_point,
    orientation_labels_for_graph,
};
use crate::model::Model;
use crate::momentum::ThreeMomentum;
use crate::momentum::sample::{LoopIndex, LoopMomenta, MomentumSample};
use crate::processes::{Amplitude, AmplitudeGraph, CrossSection};
use crate::settings::RuntimeSettings;
use crate::utils::F;
use crate::uv::UltravioletGraph;
use clap::ValueEnum;
use color_eyre::{Result, eyre::Context};
use colored::Colorize;
use eyre::eyre;
use itertools::Itertools;
use linnet::half_edge::PowersetIterator;
use linnet::half_edge::involution::{EdgeIndex, EdgeVec, HedgePair, Orientation, SignOrZero};
use linnet::half_edge::subgraph::subset::SubSet;
use linnet::half_edge::subgraph::{Inclusion, ModifySubSet, SuBitGraph, SubSetLike, SubSetOps};
use rand::Rng;
use rayon::prelude::*;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use spenso::algebra::complex::Complex;
use symbolica::domains::atom::AtomField;
use symbolica::numerical_integration::MonteCarloRng;
use symbolica::poly::series::Series;
use symbolica::symbol;
use tabled::{
    Table, Tabled,
    builder::Builder,
    settings::{Modify, Span, Style},
};
use tracing::{debug, info_span, instrument};
use tracing_indicatif::{span_ext::IndicatifSpanExt, style::ProgressStyle};
use typed_index_collections::TiVec;

type LoopMomentumSample = TiVec<LoopIndex, ThreeMomentum<F<f64>>>;
type ProfileLmbLimits = Vec<(LoopMomentumBasis, Vec<(SubSet<LoopIndex>, i32)>)>;
const UV_PROFILE_RETRY_MAX_DOD: f64 = -0.9;

#[derive(
    Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, JsonSchema, ValueEnum,
)]
#[serde(rename_all = "kebab-case")]
pub enum UVLimitSelection {
    #[default]
    OnlyDivergent,
    All,
}

impl UVLimitSelection {
    fn includes(self, bare_dod: i32) -> bool {
        matches!(self, Self::All) || bare_dod >= 0
    }

    fn cycle_union(
        graph: &Graph,
        lmb: &LoopMomentumBasis,
        subset: &SubSet<LoopIndex>,
    ) -> SuBitGraph {
        let mut subgraph: SuBitGraph = graph.empty_subgraph();
        for (pair, edge_id, _) in graph.underlying.iter_edges() {
            if matches!(pair, HedgePair::Paired { .. })
                && subset.included_iter().any(|loop_index| {
                    lmb.edge_signatures[edge_id].internal[loop_index] != SignOrZero::Zero
                })
            {
                subgraph.add(pair);
            }
        }
        subgraph
    }

    fn all_lmb_limits(
        graph: &Graph,
        lmb: &LoopMomentumBasis,
        compatible_cuts: Option<&[SuBitGraph]>,
    ) -> Vec<(SubSet<LoopIndex>, i32)> {
        let mut loops = PowersetIterator::<LoopIndex>::new(lmb.loop_edges.len() as u8);
        loops.next();
        loops
            .filter_map(|subset| {
                let cycle = Self::cycle_union(graph, lmb, &subset);
                let is_cut_compatible = compatible_cuts
                    .is_none_or(|cuts| cuts.iter().any(|cut| !cycle.intersects(cut)));
                is_cut_compatible.then(|| (subset, graph.compute_dod(&cycle)))
            })
            .collect()
    }

    fn lmb_limits(
        self,
        graph: &Graph,
        all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
        compatible_cuts: Option<&[SuBitGraph]>,
    ) -> Result<ProfileLmbLimits> {
        if matches!(self, Self::All) {
            return Ok(all_lmbs
                .iter()
                .cloned()
                .map(|lmb| {
                    let limits = Self::all_lmb_limits(graph, &lmb, compatible_cuts);
                    (lmb, limits)
                })
                .collect());
        }

        let mut examined = BTreeSet::new();
        let mut selected: ProfileLmbLimits = Vec::new();
        for lmb in std::iter::once(&graph.loop_momentum_basis).chain(all_lmbs.iter()) {
            let mut subsets = PowersetIterator::<LoopIndex>::new(lmb.loop_edges.len() as u8);
            subsets.next();
            let limits = subsets
                .filter_map(|subset| {
                    let cycle = Self::cycle_union(graph, lmb, &subset);
                    if !examined.insert(cycle.string_label()) {
                        return None;
                    }
                    let is_cut_compatible = compatible_cuts
                        .is_none_or(|cuts| cuts.iter().any(|cut| !cycle.intersects(cut)));
                    if !is_cut_compatible {
                        return None;
                    }
                    let dod = graph.compute_dod(&cycle);
                    self.includes(dod).then_some((subset, dod))
                })
                .collect_vec();
            if !limits.is_empty() {
                selected.push((lmb.clone(), limits));
            }
        }
        Ok(selected)
    }
}

pub struct ProfileSettings {
    pub n_points: usize,
    pub min_scale_exponent: f64,
    pub max_scale_exponent: f64,
    pub seed: u64,
    pub use_f128: bool,
    pub analyse_analytically: bool,
    pub orientation_mode: OrientationProfileMode,
    pub fixed_uv_ray: Option<UVProfileFixedRay>,
    pub graph_id: Option<usize>,
    pub cutkosky_cut: Option<Vec<EdgeIndex>>,
    pub selected_limits: UVLimitSelection,
}

impl Default for ProfileSettings {
    fn default() -> Self {
        ProfileSettings {
            n_points: 15,
            min_scale_exponent: 3.0,
            max_scale_exponent: 6.0,
            seed: 42,
            analyse_analytically: false,
            use_f128: false,
            orientation_mode: OrientationProfileMode::Summed,
            fixed_uv_ray: None,
            graph_id: None,
            cutkosky_cut: None,
            selected_limits: UVLimitSelection::OnlyDivergent,
        }
    }
}

#[derive(Debug, Clone)]
pub struct UVProfileFixedRay {
    directions: Vec<[f64; 3]>,
    norms: Vec<f64>,
}

impl UVProfileFixedRay {
    pub fn from_flat_components(directions: &[f64], norms: &[f64]) -> Result<Self> {
        if directions.is_empty() {
            return Err(eyre!(
                "Fixed UV ray directions cannot be empty when the option is used."
            ));
        }
        if !directions.len().is_multiple_of(3) {
            return Err(eyre!(
                "Fixed UV ray directions must contain a multiple of 3 components, got {}.",
                directions.len()
            ));
        }
        if norms.is_empty() {
            return Err(eyre!(
                "Fixed UV ray norms cannot be empty when directions are supplied."
            ));
        }

        let directions = directions
            .chunks_exact(3)
            .map(|direction| {
                let norm = (direction[0] * direction[0]
                    + direction[1] * direction[1]
                    + direction[2] * direction[2])
                    .sqrt();
                if norm == 0.0 {
                    return Err(eyre!(
                        "Fixed UV ray directions cannot contain a zero vector."
                    ));
                }
                Ok([
                    direction[0] / norm,
                    direction[1] / norm,
                    direction[2] / norm,
                ])
            })
            .collect::<Result<Vec<_>>>()?;

        if norms.iter().any(|norm| !norm.is_finite() || *norm <= 0.0) {
            return Err(eyre!(
                "Fixed UV ray norms must be finite and strictly positive."
            ));
        }

        Ok(Self {
            directions,
            norms: norms.to_vec(),
        })
    }

    fn sample(&self, loop_count: usize) -> Result<LoopMomentumSample> {
        if self.directions.len() != 1 && self.directions.len() != loop_count {
            return Err(eyre!(
                "Fixed UV ray directions contain {} ray(s), but this LMB has {} loop variable(s). Use either one direction or one direction per loop variable.",
                self.directions.len(),
                loop_count
            ));
        }
        if self.norms.len() != 1 && self.norms.len() != loop_count {
            return Err(eyre!(
                "Fixed UV ray norms contain {} value(s), but this LMB has {} loop variable(s). Use either one norm or one norm per loop variable.",
                self.norms.len(),
                loop_count
            ));
        }

        (0..loop_count)
            .map(|i| {
                let direction = self.directions[if self.directions.len() == 1 { 0 } else { i }];
                let norm = self.norms[if self.norms.len() == 1 { 0 } else { i }];
                Ok(ThreeMomentum {
                    px: F(norm * direction[0]),
                    py: F(norm * direction[1]),
                    pz: F(norm * direction[2]),
                })
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use std::sync::OnceLock;

    use linnet::half_edge::involution::EdgeIndex;
    use linnet::half_edge::subgraph::{
        Inclusion, ModifySubSet, SuBitGraph, SubSetLike, subset::SubSet,
    };
    use typed_index_collections::TiVec;

    use crate::dot;
    use crate::graph::parse::from_dot::IntoGraph;
    use crate::graph::{Graph, LMBext, LmbIndex, LoopMomentumBasis};
    use crate::initialisation::test_initialise;
    use crate::momentum::ThreeMomentum;
    use crate::momentum::sample::{
        BareMomentumSample, ExternalThreeMomenta, LoopIndex, LoopMomenta, MomentumSample,
    };
    use crate::utils::F;

    use super::{
        SubsetOrientationInput, UVLimitSelection, UVProfile, UVProfileFailure, UVProfileFixedRay,
        UVProfilePassFail, UVSamplingResult,
    };

    static THETA_GRAPH: OnceLock<(Graph, TiVec<LmbIndex, LoopMomentumBasis>)> = OnceLock::new();

    fn theta_graph() -> (Graph, TiVec<LmbIndex, LoopMomentumBasis>) {
        THETA_GRAPH
            .get_or_init(|| {
                test_initialise().unwrap();
                let mut graph: Graph = dot!(digraph theta_profile {
                    edge [num=1 mass=0]
                    node [num=1]
                    A -> B [id=0]
                    A -> B [id=1]
                    A -> B [id=2]
                })
                .unwrap();
                let lmbs = graph.generate_loop_momentum_bases();
                graph.loop_momentum_basis = lmbs[LmbIndex::from(0)].clone();
                (graph, lmbs)
            })
            .clone()
    }

    #[test]
    fn selected_limits_default_to_bare_uv_divergences() {
        let divergent = UVLimitSelection::default();
        assert!(divergent.includes(1));
        assert!(divergent.includes(0));
        assert!(!divergent.includes(-1));
        assert!(UVLimitSelection::All.includes(-1));
    }

    #[test]
    fn cycle_union_uses_lmb_signatures_for_a_self_loop() {
        let _ = theta_graph();
        let graph: Graph = dot!(digraph tadpole_profile {
            edge [num=1 mass=0]
            node [num=1]
            A -> A [id=0]
        })
        .unwrap();
        assert_eq!(graph.loop_momentum_basis.loop_edges.len(), 1);

        let mut subset = SubSet::empty(1);
        subset.add(LoopIndex::from(0));
        let cycle = UVLimitSelection::cycle_union(&graph, &graph.loop_momentum_basis, &subset);
        assert!(cycle.includes(&graph[&EdgeIndex::from(0)].1));
    }

    #[test]
    fn only_divergent_profiles_all_theta_cycles_once() {
        let (graph, lmbs) = theta_graph();
        let selected = UVLimitSelection::OnlyDivergent
            .lmb_limits(&graph, &lmbs, None)
            .unwrap();
        let represented = selected
            .iter()
            .flat_map(|(lmb, limits)| {
                limits.iter().map(|(subset, _)| {
                    UVLimitSelection::cycle_union(&graph, lmb, subset).string_label()
                })
            })
            .collect::<std::collections::BTreeSet<_>>();
        let selected_limit_count: usize = selected.iter().map(|(_, limits)| limits.len()).sum();

        assert_eq!(selected_limit_count, 4);
        assert_eq!(represented.len(), selected_limit_count);
        assert!(
            selected.len() > 1,
            "the third one-loop theta cycle requires another generated physical basis"
        );
        assert!(
            selected
                .iter()
                .all(|(selected_lmb, _)| lmbs.iter().any(|lmb| lmb == selected_lmb)),
            "divergent-only selection must not synthesize loop-momentum bases"
        );
        assert_eq!(
            selected
                .iter()
                .flat_map(|(_, limits)| limits)
                .filter(|(subset, _)| subset.n_included() == 1)
                .count(),
            3,
            "all three logarithmically divergent one-loop theta cycles must be selected"
        );

        let generation_only: TiVec<LmbIndex, LoopMomentumBasis> =
            vec![graph.loop_momentum_basis.clone()].into();
        let generation_selected = UVLimitSelection::OnlyDivergent
            .lmb_limits(&graph, &generation_only, None)
            .unwrap();
        assert_eq!(
            generation_selected
                .iter()
                .flat_map(|(_, limits)| limits)
                .filter(|(subset, _)| subset.n_included() == 1)
                .count(),
            2,
            "an unrepresentable theta cycle must be skipped when only the generation basis is physical"
        );

        let exhaustive = UVLimitSelection::All
            .lmb_limits(&graph, &lmbs, None)
            .unwrap();
        assert_eq!(exhaustive.len(), lmbs.len());
        assert!(exhaustive.iter().all(
            |(_, limits)| limits.len() == (1 << graph.loop_momentum_basis.loop_edges.len()) - 1
        ));

        let mut cut_union: SuBitGraph = graph.empty_subgraph();
        cut_union.add(graph[&EdgeIndex::from(0)].1);
        let cut_compatible_exhaustive = UVLimitSelection::All
            .lmb_limits(&graph, &lmbs, Some(std::slice::from_ref(&cut_union)))
            .unwrap();
        let cut_compatible_count: usize = cut_compatible_exhaustive
            .iter()
            .map(|(_, limits)| limits.len())
            .sum();
        assert!(cut_compatible_count > 0 && cut_compatible_count < exhaustive.len() * 3);
        assert!(
            cut_compatible_exhaustive
                .iter()
                .all(|(lmb, limits)| limits.iter().all(|(subset, _)| {
                    !UVLimitSelection::cycle_union(&graph, lmb, subset).intersects(&cut_union)
                }))
        );
    }

    #[test]
    fn candidate_lmb_sample_is_transformed_to_the_generation_lmb() {
        let (graph, lmbs) = theta_graph();
        let candidate = lmbs
            .iter()
            .find(|lmb| *lmb != &graph.loop_momentum_basis)
            .unwrap();
        let loop_moms: LoopMomenta<F<f64>> = [
            ThreeMomentum::new(F(0.25), F(-0.5), F(0.75)),
            ThreeMomentum::new(F(-0.125), F(0.375), F(0.625)),
        ]
        .into_iter()
        .collect();
        let sample = MomentumSample {
            sample: BareMomentumSample {
                loop_moms: loop_moms.clone(),
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
        let subset = SubSet::empty(candidate.loop_edges.len());
        let candidate_sample = candidate
            .loop_edges
            .iter()
            .map(|_| ThreeMomentum::default())
            .collect();
        let input = SubsetOrientationInput {
            graph_id: 0,
            subset: &subset,
            lmb: candidate,
            generation_lmb: &graph.loop_momentum_basis,
            sample: &candidate_sample,
            orientation: None,
        };
        let transformed: LoopMomenta<F<f64>> =
            input.generation_loop_momenta(&sample).into_iter().collect();
        let externals: ExternalThreeMomenta<F<f64>> = Vec::new().into();

        for (_, edge_id, _) in graph.underlying.iter_edges() {
            assert_eq!(
                candidate.edge_signatures[edge_id].compute_momentum(&loop_moms, &externals),
                graph.loop_momentum_basis.edge_signatures[edge_id]
                    .compute_momentum(&transformed, &externals)
            );
        }
    }

    #[test]
    fn analysis_preserves_selected_graph_and_cut_identity() {
        let profile = UVProfile {
            per_graph: vec![UVSamplingResult {
                graph_index: 2,
                graph_name: "GL2".to_string(),
                cutkosky_cut: Some(vec![EdgeIndex::from(5), EdgeIndex::from(2)]),
                per_lmb: Vec::new(),
            }],
            scales: Vec::new(),
            allow_vanishing_missing_fits: true,
        };

        let analysis = profile.analyse();
        let graph = &analysis.graphs[0];
        assert_eq!(graph.graph_index, 2);
        assert_eq!(graph.graph_name, "GL2");
        assert_eq!(
            graph.cutkosky_cut,
            Some(vec![EdgeIndex::from(5), EdgeIndex::from(2)])
        );
    }

    #[test]
    fn pass_fail_display_lists_only_failures_with_context() {
        let report = UVProfilePassFail {
            max_dod: -0.9,
            total: 3,
            failed: 1,
            failures: vec![UVProfileFailure {
                graph_index: 2,
                graph_name: "GL2".to_string(),
                cutkosky_cut: Some(vec![EdgeIndex::from(5), EdgeIndex::from(2)]),
                lmb_index: 0,
                fixed: vec![EdgeIndex::from(1)],
                free: vec![EdgeIndex::from(3)],
                orientation_label: None,
                reason: "estimated DOD exceeds threshold".to_string(),
            }],
        };

        let rendered = report.to_string();
        assert!(rendered.contains("UV limit tests:"));
        assert!(rendered.contains("GL2"));
        assert!(rendered.contains("e5,e2"));
        assert!(rendered.contains("estimated DOD exceeds threshold"));
    }

    #[test]
    fn fixed_uv_ray_repeats_single_direction_and_norm() {
        let ray = UVProfileFixedRay::from_flat_components(&[0.0, 0.0, 2.0], &[3.0]).unwrap();
        let sample = ray.sample(2).unwrap();

        assert_eq!(sample.len(), 2);
        for momentum in &sample {
            assert_eq!(momentum.px.0, 0.0);
            assert_eq!(momentum.py.0, 0.0);
            assert_eq!(momentum.pz.0, 3.0);
        }
    }

    #[test]
    fn fixed_uv_ray_rejects_invalid_input_shapes() {
        assert!(UVProfileFixedRay::from_flat_components(&[1.0, 0.0], &[1.0]).is_err());
        assert!(UVProfileFixedRay::from_flat_components(&[0.0, 0.0, 0.0], &[1.0]).is_err());
        assert!(UVProfileFixedRay::from_flat_components(&[1.0, 0.0, 0.0], &[0.0]).is_err());
    }
}

pub fn logspace(start: f64, stop: f64, num: usize, base: f64) -> Vec<f64> {
    let log_start = start;
    let log_stop = stop;
    let step = (log_stop - log_start) / (num - 1) as f64;

    (0..num)
        .map(|i| {
            let exponent = log_start + step * i as f64;
            base.powf(exponent)
        })
        .collect()
}

fn lmb_seed(base_seed: u64, graph_id: usize, lmb_index: usize) -> u64 {
    base_seed
        .wrapping_add((graph_id as u64).wrapping_mul(0x9E3779B97F4A7C15))
        .wrapping_add(lmb_index as u64)
}

struct UVProfileRunner<'a> {
    integrand: &'a Arc<Mutex<ProcessIntegrand>>,
    scales: &'a [f64],
    model: &'a Model,
    settings: &'a RuntimeSettings,
    profile_settings: &'a ProfileSettings,
    base_seed: u64,
    selected_cut_id: Option<usize>,
}

pub trait UVProfileable {
    fn profile(
        &mut self,
        model: &Model,
        // settings: &RuntimeSettings,
        profile_settings: &ProfileSettings,
    ) -> Result<UVProfile>;
}

impl UVProfileable for Amplitude {
    #[instrument(skip_all)]
    fn profile(
        &mut self,
        model: &Model,
        // settings: &RuntimeSettings,
        profile_settings: &ProfileSettings,
    ) -> Result<UVProfile> {
        if profile_settings.cutkosky_cut.is_some() {
            return Err(eyre!(
                "Cutkosky-cut selection is only supported for cross sections"
            ));
        }

        let scales = logspace(
            profile_settings.min_scale_exponent,
            profile_settings.max_scale_exponent,
            profile_settings.n_points,
            10.0,
        );

        let settings = self
            .integrand
            .as_ref()
            .ok_or(eyre!("Integrand Not built yet"))?
            .get_settings()
            .clone();
        let graph_inputs = self
            .graphs
            .iter()
            .enumerate()
            .filter(|(graph_id, _)| profile_settings.graph_id.is_none_or(|id| id == *graph_id))
            .collect::<Vec<_>>();
        if let Some(graph_id) = profile_settings.graph_id
            && graph_inputs.is_empty()
        {
            return Err(eyre!(
                "Graph id {} is out of range for amplitude '{}'; it has {} graphs",
                graph_id,
                self.name,
                self.graphs.len()
            ));
        }

        let base_seed = profile_settings.seed;
        let integrand = Arc::new(Mutex::new(self.integrand.take().unwrap()));

        let profile_span = info_span!("Profiling graphs", indicatif.pb_show = true);
        profile_span.pb_set_style(&ProgressStyle::with_template(
            "{wide_bar} {pos}/{len} {msg}",
        )?);

        profile_span.pb_set_length(graph_inputs.len() as u64);
        profile_span.pb_set_message("Profiling graphs");
        profile_span.pb_set_finish_message("all graphs profiled");
        let _profile_span_enter = profile_span.enter();

        let runner = UVProfileRunner {
            integrand: &integrand,
            scales: &scales,
            model,
            settings: &settings,
            profile_settings,
            base_seed,
            selected_cut_id: None,
        };

        let per_graph = graph_inputs
            .par_iter()
            .map(|(graph_id, graph)| {
                let res = runner.sample_graph(*graph_id, graph)?;
                profile_span.pb_inc(1);
                Ok(res)
            })
            .collect::<Result<Vec<_>>>()?;

        drop(_profile_span_enter);
        drop(profile_span);

        let integrand = Arc::try_unwrap(integrand)
            .map_err(|_| color_eyre::eyre::eyre!("integrand still shared"))?
            .into_inner()
            .expect("integrand mutex poisoned");
        self.integrand = Some(integrand);

        Ok(UVProfile {
            per_graph,
            scales,
            allow_vanishing_missing_fits: false,
        })
        // results.push((inspect_res, analytic_res));
    }
}

impl UVProfileable for CrossSection {
    #[instrument(skip_all)]
    fn profile(&mut self, model: &Model, profile_settings: &ProfileSettings) -> Result<UVProfile> {
        let scales = logspace(
            profile_settings.min_scale_exponent,
            profile_settings.max_scale_exponent,
            profile_settings.n_points,
            10.0,
        );

        let integrand = self
            .integrand
            .as_ref()
            .ok_or(eyre!("Integrand Not built yet"))?;
        let settings = integrand.get_settings().clone();
        let (graph_inputs, selected_cut_id) = match integrand {
            ProcessIntegrand::CrossSection(cross_section) => {
                if let Some(graph_id) = profile_settings.graph_id
                    && graph_id >= cross_section.data.graph_terms.len()
                {
                    return Err(eyre!(
                        "Graph id {} is out of range for cross section '{}'; it has {} graphs",
                        graph_id,
                        self.name,
                        cross_section.data.graph_terms.len()
                    ));
                }
                if profile_settings.cutkosky_cut.is_some() && profile_settings.graph_id.is_none() {
                    return Err(eyre!("A Cutkosky-cut selection requires selecting a graph"));
                }

                let selected_cut_id = profile_settings
                    .cutkosky_cut
                    .as_ref()
                    .map(|requested_edges| {
                        let graph_id = profile_settings
                            .graph_id
                            .expect("Cutkosky-cut selection requires a graph");
                        let graph_term = &cross_section.data.graph_terms[graph_id];
                        let requested_edges = requested_edges.iter().copied().sorted().collect_vec();
                        let matching_cuts = graph_term
                            .cuts
                            .iter_enumerated()
                            .filter_map(|(cut_id, cut)| {
                                let cut_edges = graph_term
                                    .graph
                                    .underlying
                                    .iter_edges_of(&cut.cut)
                                    .map(|(_, edge_id, _)| edge_id)
                                    .sorted()
                                    .collect_vec();
                                (cut_edges == requested_edges).then_some(cut_id)
                            })
                            .collect_vec();
                        let cut_id = match matching_cuts.as_slice() {
                            [cut_id] => *cut_id,
                            [] => {
                                let available = graph_term
                                    .cuts
                                    .iter_enumerated()
                                    .map(|(cut_id, cut)| {
                                        let edges = graph_term
                                            .graph
                                            .underlying
                                            .iter_edges_of(&cut.cut)
                                            .map(|(_, edge_id, _)| edge_id.to_string())
                                            .sorted()
                                            .join(",");
                                        format!("{}:[{}]", cut_id.0, edges)
                                    })
                                    .join(", ");
                                return Err(eyre!(
                                    "No Cutkosky cut with edges [{}] exists for graph '{}'. Available cuts: {}",
                                    requested_edges.iter().map(ToString::to_string).join(","),
                                    graph_term.graph.name,
                                    available
                                ));
                            }
                            _ => {
                                return Err(eyre!(
                                    "Cutkosky-cut edges [{}] ambiguously match cuts {} for graph '{}'",
                                    requested_edges.iter().map(ToString::to_string).join(","),
                                    matching_cuts.iter().map(|cut_id| cut_id.0).join(","),
                                    graph_term.graph.name
                                ));
                            }
                        };
                        let cut_group = graph_term
                            .cut_group_data
                            .cut_groups
                            .iter()
                            .find(|cut_group| cut_group.cuts.contains(&cut_id))
                            .ok_or_else(|| {
                                eyre!(
                                    "Cut {} is not assigned to a cut group for graph '{}'",
                                    cut_id.0,
                                    graph_term.graph.name
                                )
                            })?;
                        // LU evaluates one contribution per cut group and tags its generated
                        // event with the group's first (representative) physical cut.
                        cut_group
                            .cuts
                            .first()
                            .map(|cut_id| cut_id.0)
                            .ok_or_else(|| {
                                eyre!(
                                    "The selected cut group is empty for graph '{}'",
                                    graph_term.graph.name
                                )
                            })
                    })
                    .transpose()?;

                let graph_inputs = cross_section
                    .data
                    .graph_terms
                    .iter()
                    .enumerate()
                    .filter(|(graph_id, _)| {
                        profile_settings.graph_id.is_none_or(|id| id == *graph_id)
                    })
                    .map(|(graph_id, graph_term)| {
                        let cuts = graph_term
                            .cut_group_data
                            .cut_groups
                            .iter()
                            .filter(|group| {
                                selected_cut_id.is_none_or(|cut_id| {
                                    group.cuts.first().is_some_and(|id| id.0 == cut_id)
                                })
                            })
                            .map(|group| {
                                group
                                    .cuts
                                    .iter()
                                    .map(|cut_id| graph_term.cuts[*cut_id].cut.as_subgraph())
                                    .fold(
                                        graph_term.graph.empty_subgraph(),
                                        |union: SuBitGraph, cut| union.union(&cut),
                                    )
                            })
                            .collect_vec();
                        (
                            graph_id,
                            graph_term.graph.clone(),
                            graph_term.lmbs.clone(),
                            cuts,
                        )
                    })
                    .collect::<Vec<_>>();
                (graph_inputs, selected_cut_id)
            }
            ProcessIntegrand::Amplitude(_) => {
                unreachable!("cross-section UV profiling expects cross-section integrands")
            }
        };
        let base_seed = profile_settings.seed;
        let integrand = Arc::new(Mutex::new(self.integrand.take().unwrap()));

        let profile_span = info_span!("Profiling graphs", indicatif.pb_show = true);
        profile_span.pb_set_style(&ProgressStyle::with_template(
            "{wide_bar} {pos}/{len} {msg}",
        )?);
        profile_span.pb_set_length(graph_inputs.len() as u64);
        profile_span.pb_set_message("Profiling graphs");
        profile_span.pb_set_finish_message("all graphs profiled");
        let _profile_span_enter = profile_span.enter();

        let runner = UVProfileRunner {
            integrand: &integrand,
            scales: &scales,
            model,
            settings: &settings,
            profile_settings,
            base_seed,
            selected_cut_id,
        };

        let per_graph = graph_inputs
            .par_iter()
            .map(|(graph_id, graph, lmbs, cuts)| {
                let res = runner.sample_cross_section_graph(*graph_id, graph, lmbs, cuts)?;
                profile_span.pb_inc(1);
                Ok(res)
            })
            .collect::<Result<Vec<_>>>()?;

        drop(_profile_span_enter);
        drop(profile_span);

        let integrand = Arc::try_unwrap(integrand)
            .map_err(|_| color_eyre::eyre::eyre!("integrand still shared"))?
            .into_inner()
            .expect("integrand mutex poisoned");
        self.integrand = Some(integrand);

        Ok(UVProfile {
            per_graph,
            scales,
            allow_vanishing_missing_fits: true,
        })
    }
}

pub struct UVProfile {
    pub per_graph: Vec<UVSamplingResult>,
    pub scales: Vec<f64>,
    pub allow_vanishing_missing_fits: bool,
}

impl UVProfile {
    pub fn analyse(&self) -> UVProfileAnalysis {
        let graphs = self
            .per_graph
            .iter()
            .map(|graph| {
                let graph_index = graph.graph_index;
                let lmbs: Vec<UVProfileLmbAnalysis> = graph
                    .per_lmb
                    .iter()
                    .enumerate()
                    .map(|(lmb_index, lmb)| {
                        let lmb_label = lmb_label(&lmb.lmb);
                        let subsets: Vec<UVProfileSubsetAnalysis> = lmb
                            .per_subsets
                            .iter()
                            .enumerate()
                            .map(|(subset_index, (subset, subset_result))| {
                                let mut not_included: SubSet<LoopIndex> =
                                    SubSet::full(subset.size());
                                not_included.subtract_with(subset);
                                let free: Vec<EdgeIndex> = subset
                                    .included_iter()
                                    .map(|loop_index| lmb.lmb.loop_edges[loop_index])
                                    .collect();

                                let fixed: Vec<EdgeIndex> = not_included
                                    .included_iter()
                                    .map(|loop_index| lmb.lmb.loop_edges[loop_index])
                                    .collect();
                                let analysis = subset_result.analyse(&self.scales);
                                let per_orientation_inspect_entries =
                                    analysis.per_orientation_inspect_entries();
                                let analytic_entries =
                                    analysis.analytic.as_ref().and_then(|analytic| {
                                        let entries = analytic
                                            .per_orientation
                                            .iter()
                                            .map(|(orientation, orientation_analysis)| {
                                                let (orientation_edges, orientation_signs) =
                                                    orientation_signs(orientation);
                                                UVProfileAnalyticEntry {
                                                    graph_index,
                                                    lmb_index,
                                                    subset_index,
                                                    fixed: fixed.clone(),
                                                    free: free.clone(),
                                                    orientation_edges,
                                                    orientation_signs,
                                                    is_constant: orientation_analysis.is_constant,
                                                    leading_coef: orientation_analysis
                                                        .leading_coef
                                                        .to_string(),
                                                }
                                            })
                                            .collect::<Vec<_>>();
                                        if entries.is_empty() {
                                            None
                                        } else {
                                            Some(entries)
                                        }
                                    });
                                UVProfileSubsetAnalysis {
                                    subset_index,
                                    fixed,
                                    free,
                                    initial_dod: subset_result.initial_dod,
                                    analysis,
                                    per_orientation_inspect_entries,
                                    analytic_entries,
                                }
                            })
                            .collect();

                        UVProfileLmbAnalysis {
                            lmb_index,
                            lmb_label,
                            subsets,
                        }
                    })
                    .collect();

                UVProfileGraphAnalysis {
                    graph_index,
                    graph_name: graph.graph_name.clone(),
                    cutkosky_cut: graph.cutkosky_cut.clone(),
                    lmbs,
                }
            })
            .collect();

        UVProfileAnalysis {
            scales: self.scales.clone(),
            graphs,
            allow_vanishing_missing_fits: self.allow_vanishing_missing_fits,
        }
    }

    pub fn write_profile_data<P: AsRef<Path>>(
        &self,
        _settings: &ProfileSettings,
        out_dir: P,
    ) -> Result<()> {
        self.analyse().write_profile_data(out_dir)
    }

    pub fn write_typst_bundle<P: AsRef<Path>>(
        &self,
        settings: &ProfileSettings,
        out_dir: P,
    ) -> Result<()> {
        self.write_profile_data(settings, out_dir)
    }

    pub fn pass_fail(&self, max_dod: f64, _settings: &ProfileSettings) -> UVProfilePassFail {
        self.analyse().pass_fail(max_dod)
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct UVProfileAnalysis {
    pub scales: Vec<f64>,
    pub graphs: Vec<UVProfileGraphAnalysis>,
    #[serde(skip_serializing)]
    pub allow_vanishing_missing_fits: bool,
}

#[derive(Debug, Clone, Serialize)]
pub struct UVProfileGraphAnalysis {
    pub graph_index: usize,
    pub graph_name: String,
    pub cutkosky_cut: Option<Vec<EdgeIndex>>,
    pub lmbs: Vec<UVProfileLmbAnalysis>,
}

#[derive(Debug, Clone, Serialize)]
pub struct UVProfileLmbAnalysis {
    pub lmb_index: usize,
    pub lmb_label: String,
    pub subsets: Vec<UVProfileSubsetAnalysis>,
}

#[derive(Debug, Clone, Serialize)]
pub struct UVProfileSubsetAnalysis {
    pub subset_index: usize,
    pub fixed: Vec<EdgeIndex>,
    pub free: Vec<EdgeIndex>,
    pub initial_dod: i32,
    // pub subset_label: String,
    pub analysis: Analysis,
    pub per_orientation_inspect_entries: Option<Vec<UVProfileOrientationInspectEntry>>,
    pub analytic_entries: Option<Vec<UVProfileAnalyticEntry>>,
}

impl UVProfileSubsetAnalysis {
    pub fn estimated_dod(&self) -> Option<i64> {
        self.analysis
            .inspect_level
            .as_ref()
            .map(|analysis| analysis.estimated_dod)
    }

    pub fn bare_dod_matches_estimate(&self) -> bool {
        self.estimated_dod() == Some(i64::from(self.initial_dod))
    }
}

// #[derive(Debug, Clone, Serialize)]
// pub struct UVProfilePoint {
//     norm:f64,
// }

#[derive(Debug, Clone, Serialize)]
pub struct UVProfileSummary {
    pub subset_count: usize,
    pub slope_min: Option<f64>,
    pub slope_max: Option<f64>,
    pub dod_min: Option<i64>,
    pub dod_max: Option<i64>,
    pub stable_min: Option<f64>,
    pub stable_max: Option<f64>,
}

#[derive(Debug, Clone, Serialize)]
pub struct UVProfilePassFail {
    pub max_dod: f64,
    pub total: usize,
    pub failed: usize,
    pub failures: Vec<UVProfileFailure>,
}

#[derive(Debug, Clone, Serialize)]
pub struct UVProfileFailure {
    pub graph_index: usize,
    pub graph_name: String,
    pub cutkosky_cut: Option<Vec<EdgeIndex>>,
    pub lmb_index: usize,
    pub fixed: Vec<EdgeIndex>,
    pub free: Vec<EdgeIndex>,
    pub orientation_label: Option<String>,
    pub reason: String,
}

impl Display for UVProfilePassFail {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let status = if self.failed == 0 {
            "PASS".green().bold()
        } else {
            "FAIL".red().bold()
        };
        writeln!(
            f,
            "UV limit tests: {} ({}/{})",
            status,
            self.total.saturating_sub(self.failed),
            self.total
        )?;

        if self.failures.is_empty() {
            return Ok(());
        }

        let mut table = Builder::new();
        table.push_record([
            "status",
            "graph",
            "cut edges",
            "LMB",
            "fixed",
            "→ ∞",
            "orientation",
            "reason",
        ]);
        for failure in &self.failures {
            table.push_record([
                "FAIL".red().bold().to_string(),
                format!("{}:{}", failure.graph_index, failure.graph_name),
                failure
                    .cutkosky_cut
                    .as_ref()
                    .map(|edges| edges.iter().map(ToString::to_string).join(","))
                    .unwrap_or_else(|| "all".to_string()),
                failure.lmb_index.to_string(),
                format!(
                    "{{{}}}",
                    failure.fixed.iter().map(ToString::to_string).join(",")
                ),
                format!(
                    "{{{}}}",
                    failure.free.iter().map(ToString::to_string).join(",")
                ),
                failure
                    .orientation_label
                    .clone()
                    .unwrap_or_else(|| "sum".to_string()),
                failure.reason.red().to_string(),
            ]);
        }

        write!(f, "{}", table.build().with(Style::rounded()))
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct UVProfileAnalyticEntry {
    pub graph_index: usize,
    pub lmb_index: usize,
    pub subset_index: usize,
    pub fixed: Vec<EdgeIndex>,
    pub free: Vec<EdgeIndex>,
    pub orientation_edges: Vec<String>,
    pub orientation_signs: Vec<String>,
    pub is_constant: bool,
    pub leading_coef: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct UVProfileOrientationInspectEntry {
    pub orientation_label: String,
    pub analysis: Option<InspectAnalysis>,
}

#[derive(Tabled)]
struct UVProfileSubsetRow {
    #[tabled(rename = "fixed")]
    fixed: String,
    #[tabled(rename = "→ ∞")]
    free: String,
    #[tabled(rename = "slope")]
    slope: String,
    #[tabled(rename = "r2")]
    r_squared: String,
    #[tabled(rename = "DOD")]
    estimated_dod: String,
    #[tabled(rename = "bare DOD")]
    initial_dod: String,
    #[tabled(rename = "inspect")]
    inspect: String,
}

#[derive(Tabled)]
struct UVProfileOrientationSubsetRow {
    #[tabled(rename = "fixed")]
    fixed: String,
    #[tabled(rename = "→ ∞")]
    free: String,
    #[tabled(rename = "orientation")]
    orientation_label: String,
    #[tabled(rename = "slope")]
    slope: String,
    #[tabled(rename = "r2")]
    r_squared: String,
    #[tabled(rename = "DOD")]
    estimated_dod: String,
    #[tabled(rename = "bare DOD")]
    initial_dod: String,
    #[tabled(rename = "inspect")]
    inspect: String,
}

#[derive(Debug, Clone, Serialize)]
struct FitResult {
    slope: f64,
    points: Vec<f64>,
    points_detail: Vec<EvaluationResult>,
    intercept: f64,
    r_squared: f64,
}

impl UVProfileAnalysis {
    pub fn tables_per_graph(&self, max_dod: f64) -> Vec<Table> {
        self.graphs
            .iter()
            .map(|graph| {
                let rows = graph
                    .lmbs
                    .iter()
                    .flat_map(|lmb| {
                        lmb.subsets.iter().map(|subset| {
                            let (slope, r_squared, estimated_dod) =
                                match &subset.analysis.inspect_level {
                                    Some(analysis) => {
                                        let r2_text = format!("{:.3}", analysis.result.r_squared);
                                        let r2_text = if analysis.result.r_squared >= 0.99 {
                                            r2_text.green()
                                        } else {
                                            r2_text.red()
                                        }
                                        .to_string();

                                        let dod = analysis.estimated_dod;
                                        let dod_text = dod.to_string();
                                        let dod_text = if (dod as f64) <= max_dod {
                                            dod_text.green()
                                        } else {
                                            dod_text.red()
                                        }
                                        .to_string();

                                        (format!("{:.6}", analysis.result.slope), r2_text, dod_text)
                                    }
                                    None => ("-".to_string(), "-".to_string(), "-".to_string()),
                                };

                            UVProfileSubsetRow {
                                fixed: format!(
                                    "{{{}}}",
                                    subset.fixed.iter().map(ToString::to_string).join(",")
                                ),
                                free: format!(
                                    "{{{}}}",
                                    subset.free.iter().map(ToString::to_string).join(",")
                                ),
                                slope,
                                r_squared,
                                estimated_dod,
                                initial_dod: if subset.initial_dod >= 0 {
                                    subset.initial_dod.to_string().red().to_string()
                                } else {
                                    subset.initial_dod.to_string().green().to_string()
                                },
                                inspect: inspect_retry_label(
                                    subset.analysis.inspect_level.as_ref(),
                                ),
                            }
                        })
                    })
                    .collect::<Vec<_>>();

                let mut table = Table::new(rows);
                table.with(Style::rounded());
                table
            })
            .collect()
    }

    pub fn analytic_tables_per_graph(&self) -> Vec<Option<Table>> {
        self.graphs
            .iter()
            .map(|graph| {
                let mut groups: Vec<Vec<&UVProfileAnalyticEntry>> = Vec::new();
                let mut orientation_headers: Option<Vec<String>> = None;

                for lmb in &graph.lmbs {
                    for subset in &lmb.subsets {
                        if let Some(entries) = &subset.analytic_entries {
                            if orientation_headers.is_none() {
                                orientation_headers =
                                    entries.first().map(|entry| entry.orientation_edges.clone());
                            }
                            groups.push(entries.iter().collect());
                        }
                    }
                }

                if groups.is_empty() {
                    return None;
                }

                let mut builder = Builder::new();
                let mut header = vec!["fixed".to_string(), "→ ∞".to_string()];
                header.extend(orientation_headers.unwrap_or_default());
                header.extend(["const".to_string(), "leading coef".to_string()]);
                builder.push_record(header);

                let mut span_ops: Vec<(usize, usize, usize)> = Vec::new();
                let mut row_index = 1;

                for group in groups {
                    let span_len = group.len();
                    let start_row = row_index;

                    for (i, entry) in group.into_iter().enumerate() {
                        let fixed = if i == 0 {
                            format!(
                                "{{{}}}",
                                entry.fixed.iter().map(ToString::to_string).join(",")
                            )
                        } else {
                            String::new()
                        };
                        let free = if i == 0 {
                            format!(
                                "{{{}}}",
                                entry.free.iter().map(ToString::to_string).join(",")
                            )
                        } else {
                            String::new()
                        };
                        let mut row = vec![fixed, free];
                        row.extend(entry.orientation_signs.iter().cloned());
                        row.extend([entry.is_constant.to_string(), "".to_string()]);
                        builder.push_record(row);
                        row_index += 1;
                    }

                    if span_len > 1 {
                        span_ops.push((start_row, 0, span_len));
                        span_ops.push((start_row, 1, span_len));
                    }
                }

                let mut table = builder.build();
                for (row, col, span_len) in span_ops {
                    table.with(Modify::new((row, col)).with(Span::row(span_len as isize)));
                }
                table.with(Style::rounded());
                Some(table)
            })
            .collect()
    }

    pub fn per_orientation_tables_per_graph(&self, max_dod: f64) -> Vec<Option<Table>> {
        self.graphs
            .iter()
            .map(|graph| {
                let rows = graph
                    .lmbs
                    .iter()
                    .flat_map(|lmb| {
                        lmb.subsets.iter().flat_map(|subset| {
                            subset
                                .per_orientation_inspect_entries
                                .iter()
                                .flatten()
                                .map(|entry| {
                                    let (slope, r_squared, estimated_dod) = match &entry.analysis {
                                        Some(analysis) => {
                                            let r2_text =
                                                format!("{:.3}", analysis.result.r_squared);
                                            let r2_text = if analysis.result.r_squared >= 0.99 {
                                                r2_text.green()
                                            } else {
                                                r2_text.red()
                                            }
                                            .to_string();

                                            let dod = analysis.estimated_dod;
                                            let dod_text = dod.to_string();
                                            let dod_text = if (dod as f64) <= max_dod {
                                                dod_text.green()
                                            } else {
                                                dod_text.red()
                                            }
                                            .to_string();

                                            (
                                                format!("{:.6}", analysis.result.slope),
                                                r2_text,
                                                dod_text,
                                            )
                                        }
                                        None => ("-".to_string(), "-".to_string(), "-".to_string()),
                                    };

                                    UVProfileOrientationSubsetRow {
                                        fixed: format!(
                                            "{{{}}}",
                                            subset.fixed.iter().map(ToString::to_string).join(",")
                                        ),
                                        free: format!(
                                            "{{{}}}",
                                            subset.free.iter().map(ToString::to_string).join(",")
                                        ),
                                        orientation_label: entry.orientation_label.clone(),
                                        slope,
                                        r_squared,
                                        estimated_dod,
                                        initial_dod: if subset.initial_dod >= 0 {
                                            subset.initial_dod.to_string().red().to_string()
                                        } else {
                                            subset.initial_dod.to_string().green().to_string()
                                        },
                                        inspect: inspect_retry_label(entry.analysis.as_ref()),
                                    }
                                })
                                .collect::<Vec<_>>()
                        })
                    })
                    .collect::<Vec<_>>();

                if rows.is_empty() {
                    None
                } else {
                    let mut table = Table::new(rows);
                    table.with(Style::rounded());
                    Some(table)
                }
            })
            .collect()
    }

    pub fn write_profile_data<P: AsRef<Path>>(&self, out_dir: P) -> Result<()> {
        let out_dir = out_dir.as_ref();
        std::fs::create_dir_all(out_dir).context("failed to create UV profile output directory")?;

        let json_path = out_dir.join("uv_profile.json");
        let json = serde_json::to_string_pretty(self)
            .context("failed to serialize UV profile analysis to JSON")?;
        std::fs::write(&json_path, json).context("failed to write UV profile JSON output")?;

        Ok(())
    }

    pub fn pass_fail(&self, max_dod: f64) -> UVProfilePassFail {
        let mut failures = Vec::new();
        let mut total = 0;

        for graph in &self.graphs {
            for lmb in &graph.lmbs {
                for subset in &lmb.subsets {
                    let per_orientation_reason_count = subset
                        .analysis
                        .per_orientation_inspect
                        .as_ref()
                        .map(|entries| {
                            entries
                                .iter()
                                .filter(|entry| {
                                    inspect_failure_reason(
                                        entry.analysis.as_ref(),
                                        entry.inspect_fit_status,
                                        max_dod,
                                        self.allow_vanishing_missing_fits,
                                    )
                                    .is_some()
                                })
                                .count()
                        })
                        .unwrap_or(0);
                    let all_orientations_pass = subset.analysis.per_orientation_inspect.is_some()
                        && per_orientation_reason_count == 0;

                    total += 1;
                    let reason = inspect_failure_reason(
                        subset.analysis.inspect_level.as_ref(),
                        subset.analysis.inspect_fit_status,
                        max_dod,
                        self.allow_vanishing_missing_fits,
                    );

                    if let Some(reason) = reason.filter(|_| !all_orientations_pass) {
                        failures.push(UVProfileFailure {
                            graph_index: graph.graph_index,
                            graph_name: graph.graph_name.clone(),
                            cutkosky_cut: graph.cutkosky_cut.clone(),
                            lmb_index: lmb.lmb_index,
                            fixed: subset.fixed.clone(),
                            free: subset.free.clone(),
                            orientation_label: None,
                            reason: reason.to_string(),
                        });
                    }

                    for entry in subset.analysis.per_orientation_inspect.iter().flatten() {
                        total += 1;
                        let reason = inspect_failure_reason(
                            entry.analysis.as_ref(),
                            entry.inspect_fit_status,
                            max_dod,
                            self.allow_vanishing_missing_fits,
                        );

                        if let Some(reason) = reason {
                            failures.push(UVProfileFailure {
                                graph_index: graph.graph_index,
                                graph_name: graph.graph_name.clone(),
                                cutkosky_cut: graph.cutkosky_cut.clone(),
                                lmb_index: lmb.lmb_index,
                                fixed: subset.fixed.clone(),
                                free: subset.free.clone(),
                                orientation_label: Some(entry.orientation_label.clone()),
                                reason: reason.to_string(),
                            });
                        }
                    }
                }
            }
        }

        UVProfilePassFail {
            max_dod,
            total,
            failed: failures.len(),
            failures,
        }
    }
}

fn inspect_failure_reason(
    analysis: Option<&InspectAnalysis>,
    fit_status: InspectFitStatus,
    max_dod: f64,
    allow_vanishing_missing_fits: bool,
) -> Option<&'static str> {
    match analysis {
        None if allow_vanishing_missing_fits && fit_status.missing_fit_is_vanishing() => None,
        None => Some("missing_fit"),
        Some(analysis) if analysis.result.slope > max_dod || analysis.result.slope.is_nan() => {
            Some("dod_exceeds_threshold")
        }
        _ => None,
    }
}

fn lmb_label(lmb: &LoopMomentumBasis) -> String {
    let edges: Vec<String> = lmb.loop_edges.iter().map(|edge| edge.to_string()).collect();
    format!("loop_edges=[{}]", edges.join(","))
}

fn orientation_signs(orientation: &OrientationData) -> (Vec<String>, Vec<String>) {
    let edges = orientation
        .orientation
        .iter()
        .map(|(edge, _)| edge.to_string())
        .collect::<Vec<_>>();
    let signs = orientation
        .orientation
        .iter()
        .map(|(_, sign)| SignOrZero::from(*sign).to_string())
        .collect::<Vec<_>>();
    (edges, signs)
}

pub struct UVSamplingResult {
    pub graph_index: usize,
    pub graph_name: String,
    pub cutkosky_cut: Option<Vec<EdgeIndex>>,
    pub per_lmb: Vec<LMBResult>,
}

impl<'a> UVProfileRunner<'a> {
    fn sample_graph(&self, graph_id: usize, g: &AmplitudeGraph) -> Result<UVSamplingResult> {
        let deferred_integrands = g.derived_data.deferred_integrands.as_ref();
        let profiles_per_orientation = self
            .profile_settings
            .orientation_mode
            .profiles_per_orientation();
        if profiles_per_orientation {
            let integrand = self.integrand.lock().expect("integrand mutex poisoned");
            match &*integrand {
                ProcessIntegrand::Amplitude(amplitude)
                    if amplitude.data.explicit_orientation_sum_only =>
                {
                    return Err(eyre!(
                        "an explicit orientation sum cannot be profiled by production orientation"
                    ));
                }
                ProcessIntegrand::Amplitude(_) => {}
                ProcessIntegrand::CrossSection(_) => {
                    unreachable!("UV profiling expects amplitudes")
                }
            }
        }
        let analytic_integrands = if !self.profile_settings.analyse_analytically {
            Vec::new()
        } else if let Some(deferred_integrands) = deferred_integrands {
            let materialized = deferred_integrands.materialize();
            let mut roots = materialized.iter();
            let (index, integrand) = roots
                .next()
                .ok_or_else(|| eyre!("deferred amplitude integrand has no root residue"))?;
            if *index != crate::cff::CutCFFIndex::new_all_none() || roots.next().is_some() {
                return Err(eyre!(
                    "deferred amplitude integrand must contain exactly one root residue"
                ));
            }
            let mut summed =
                OrientationData::new(EdgeVec::from_iter(std::iter::empty::<Orientation>()));
            summed.label = Some("explicit source-local sum".to_string());
            vec![(summed, integrand.clone())]
        } else {
            g.derived_data
                .cff_expression
                .as_ref()
                .unwrap()
                .expression
                .orientations
                .iter()
                .map(|orientation| {
                    (
                        orientation.data.clone(),
                        orientation
                            .data
                            .orientation
                            .select(&g.derived_data.all_mighty_integrand),
                    )
                })
                .collect()
        };
        let orientation_labels = if profiles_per_orientation {
            let integrand = self.integrand.lock().expect("integrand mutex poisoned");
            match &*integrand {
                ProcessIntegrand::Amplitude(amplitude) => {
                    Some(orientation_labels_for_graph(amplitude, graph_id)?)
                }
                ProcessIntegrand::CrossSection(_) => {
                    unreachable!("UV profiling expects amplitudes")
                }
            }
        } else {
            None
        };
        let all_lmbs = g
            .derived_data
            .lmbs
            .as_ref()
            .ok_or_else(|| eyre!("Loop momentum bases have not been generated"))?;
        let lmbs = self
            .profile_settings
            .selected_limits
            .lmb_limits(&g.graph, all_lmbs, None)?;

        let lmb_span = info_span!(
            "Profiling loop momentum bases",
            indicatif.pb_show = true,
            graph_id = graph_id
        );
        lmb_span.pb_set_style(
            &ProgressStyle::with_template("{wide_bar} {pos}/{len} {msg}")
                .expect("invalid progress bar template"),
        );
        lmb_span.pb_set_length(lmbs.len() as u64);
        lmb_span.pb_set_message("Profiling loop momentum bases");
        lmb_span.pb_set_finish_message("all loop momentum bases profiled");
        let _lmb_span_enter = lmb_span.enter();

        let per_lmb = lmbs
            .par_iter()
            .enumerate()
            .map(|(lmb_index, (lmb, subsets))| {
                let mut result = self.sample_lmb(
                    graph_id,
                    &g.graph,
                    lmb_index,
                    lmb,
                    subsets,
                    orientation_labels.as_deref(),
                )?;

                if self.profile_settings.analyse_analytically {
                    let orientation_limits: Vec<(
                        SubSet<LoopIndex>,
                        OrientationData,
                        Series<AtomField>,
                    )> = analytic_integrands
                        .par_iter()
                        .map(|(orientation, integrand)| {
                            g.graph
                                .all_limits(
                                    &g.graph.full_filter(),
                                    integrand,
                                    symbol!("lambd"),
                                    lmb,
                                )
                                .into_iter()
                                .map(|(limit, value)| (limit, orientation.clone(), value))
                                .collect::<Vec<_>>()
                        })
                        .reduce(Vec::new, |mut limits, mut more_limits| {
                            limits.append(&mut more_limits);
                            limits
                        });

                    for (limit, orientation, value) in orientation_limits {
                        let Some(subset) = result.per_subsets.get_mut(&limit) else {
                            continue;
                        };
                        let analytic = subset.analytic.get_or_insert_with(|| AnalyticResult {
                            per_orientations: BTreeMap::new(),
                        });
                        analytic.per_orientations.insert(orientation, value);
                    }
                }
                lmb_span.pb_inc(1);
                Ok(result)
            })
            .collect::<Result<Vec<_>>>()?;

        drop(_lmb_span_enter);
        drop(lmb_span);

        Ok(UVSamplingResult {
            graph_index: graph_id,
            graph_name: g.graph.name.clone(),
            cutkosky_cut: None,
            per_lmb,
        })
    }

    fn sample_cross_section_graph(
        &self,
        graph_id: usize,
        graph: &Graph,
        all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
        cuts: &[SuBitGraph],
    ) -> Result<UVSamplingResult> {
        if self.profile_settings.analyse_analytically {
            return Err(eyre!(
                "Analytic UV profiling is not implemented for cross sections"
            ));
        }

        let orientation_labels = if self
            .profile_settings
            .orientation_mode
            .profiles_per_orientation()
        {
            let integrand = self.integrand.lock().expect("integrand mutex poisoned");
            match &*integrand {
                ProcessIntegrand::CrossSection(cross_section) => {
                    if cross_section.data.explicit_orientation_sum_only {
                        return Err(eyre!(
                            "an explicit orientation sum cannot be profiled by production orientation"
                        ));
                    }
                    Some(orientation_labels_for_graph(cross_section, graph_id)?)
                }
                ProcessIntegrand::Amplitude(_) => {
                    unreachable!("cross-section UV profiling expects cross sections")
                }
            }
        } else {
            None
        };
        let lmbs = self
            .profile_settings
            .selected_limits
            .lmb_limits(graph, all_lmbs, Some(cuts))?;

        let lmb_span = info_span!(
            "Profiling loop momentum bases",
            indicatif.pb_show = true,
            graph_id = graph_id
        );
        lmb_span.pb_set_style(
            &ProgressStyle::with_template("{wide_bar} {pos}/{len} {msg}")
                .expect("invalid progress bar template"),
        );
        lmb_span.pb_set_length(lmbs.len() as u64);
        lmb_span.pb_set_message("Profiling loop momentum bases");
        lmb_span.pb_set_finish_message("all loop momentum bases profiled");
        let _lmb_span_enter = lmb_span.enter();

        let per_lmb = lmbs
            .par_iter()
            .enumerate()
            .map(|(lmb_index, (lmb, subsets))| {
                let result = self.sample_lmb(
                    graph_id,
                    graph,
                    lmb_index,
                    lmb,
                    subsets,
                    orientation_labels.as_deref(),
                )?;
                lmb_span.pb_inc(1);
                Ok(result)
            })
            .collect::<Result<Vec<_>>>()?;

        drop(_lmb_span_enter);
        drop(lmb_span);

        Ok(UVSamplingResult {
            graph_index: graph_id,
            graph_name: graph.name.clone(),
            cutkosky_cut: self.profile_settings.cutkosky_cut.clone(),
            per_lmb,
        })
    }
}

pub struct LMBResult {
    pub(crate) lmb: LoopMomentumBasis,
    pub(crate) per_subsets: BTreeMap<SubSet<LoopIndex>, SubSetResult>,
}

impl<'a> UVProfileRunner<'a> {
    fn sample_lmb(
        &self,
        graph_id: usize,
        graph: &Graph,
        lmb_index: usize,
        lmb: &LoopMomentumBasis,
        subsets: &[(SubSet<LoopIndex>, i32)],
        orientation_labels: Option<&[String]>,
    ) -> Result<LMBResult> {
        let sample: LoopMomentumSample =
            if let Some(fixed_uv_ray) = &self.profile_settings.fixed_uv_ray {
                fixed_uv_ray.sample(lmb.loop_edges.len())?
            } else {
                let mut rng = MonteCarloRng::new(lmb_seed(self.base_seed, graph_id, lmb_index), 0);
                lmb.loop_edges
                    .iter()
                    .map(|_| ThreeMomentum {
                        px: F(rng.random_range(
                            -self.settings.kinematics.e_cm..self.settings.kinematics.e_cm,
                        )),
                        py: F(rng.random_range(
                            -self.settings.kinematics.e_cm..self.settings.kinematics.e_cm,
                        )),
                        pz: F(rng.random_range(
                            -self.settings.kinematics.e_cm..self.settings.kinematics.e_cm,
                        )),
                    })
                    .collect()
            };

        let subset_span = info_span!(
            "Profiling subsets",
            indicatif.pb_show = true,
            graph_id = graph_id,
            lmb_index = lmb_index
        );
        subset_span.pb_set_style(
            &ProgressStyle::with_template("{wide_bar} {pos}/{len} {msg}")
                .expect("invalid progress bar template"),
        );
        subset_span.pb_set_length(subsets.len() as u64);
        subset_span.pb_set_message("Profiling subsets");
        let _subset_span_enter = subset_span.enter();

        let per_subsets_vec: Vec<(SubSet<LoopIndex>, SubSetResult)> = subsets
            .par_iter()
            .map_init(
                || {
                    let mut integrand = self
                        .integrand
                        .lock()
                        .expect("integrand mutex poisoned")
                        .clone();
                    // A selected LU cut is projected from the production event weights.
                    if self.selected_cut_id.is_some() {
                        integrand.get_mut_settings().general.generate_events = true;
                    }
                    integrand
                },
                |integrand, (subset, initial_dod)| {
                    let res = self.sample_subset(
                        integrand,
                        SubsetSampleInput {
                            graph_id,
                            subset,
                            initial_dod: *initial_dod,
                            lmb,
                            generation_lmb: &graph.loop_momentum_basis,
                            sample: &sample,
                            orientation_labels,
                        },
                    )?;
                    subset_span.pb_inc(1);
                    Ok((subset.clone(), res))
                },
            )
            .collect::<Result<Vec<_>>>()?;
        let per_subsets = per_subsets_vec.into_iter().collect();

        drop(_subset_span_enter);
        drop(subset_span);

        Ok(LMBResult {
            lmb: lmb.clone(),
            per_subsets,
        })
    }
}

pub struct SubSetResult {
    pub(crate) initial_dod: i32,
    pub(crate) inspect: InspectSamples,
    pub(crate) analytic: Option<AnalyticResult>,
}

#[derive(Debug, Clone, Default)]
pub(crate) struct InspectSamples {
    pub(crate) summed: Vec<InspectResult>,
    pub(crate) summed_used_arb_prec_retry: bool,
    pub(crate) per_orientation: Vec<OrientationInspectSamples>,
}

#[derive(Debug, Clone)]
pub(crate) struct OrientationInspectSamples {
    pub(crate) label: String,
    pub(crate) inspect: Vec<InspectResult>,
    pub(crate) used_arb_prec_retry: bool,
}

#[derive(Debug, Clone)]
struct InspectRun {
    inspect: Vec<InspectResult>,
    used_arb_prec_retry: bool,
}

#[derive(Clone, Copy)]
struct SubsetOrientationInput<'a> {
    graph_id: usize,
    subset: &'a SubSet<LoopIndex>,
    lmb: &'a LoopMomentumBasis,
    generation_lmb: &'a LoopMomentumBasis,
    sample: &'a LoopMomentumSample,
    orientation: Option<usize>,
}

impl SubsetOrientationInput<'_> {
    fn generation_loop_momenta(self, sample: &MomentumSample<f64>) -> Vec<ThreeMomentum<F<f64>>> {
        sample
            .lmb_transform(self.lmb, self.generation_lmb)
            .loop_moms()
            .iter()
            .cloned()
            .collect()
    }
}

struct SubsetSampleInput<'a> {
    graph_id: usize,
    subset: &'a SubSet<LoopIndex>,
    initial_dod: i32,
    lmb: &'a LoopMomentumBasis,
    generation_lmb: &'a LoopMomentumBasis,
    sample: &'a LoopMomentumSample,
    orientation_labels: Option<&'a [String]>,
}

impl<'a> UVProfileRunner<'a> {
    fn sample_subset(
        &self,
        integrand: &mut ProcessIntegrand,
        input: SubsetSampleInput<'_>,
    ) -> Result<SubSetResult> {
        let inspect = if let Some(orientation_labels) = input.orientation_labels {
            let per_orientation = orientation_labels
                .iter()
                .enumerate()
                .map(|(orientation_id, label)| {
                    let inspect = self.sample_subset_orientation(
                        integrand,
                        SubsetOrientationInput {
                            graph_id: input.graph_id,
                            subset: input.subset,
                            lmb: input.lmb,
                            generation_lmb: input.generation_lmb,
                            sample: input.sample,
                            orientation: Some(orientation_id),
                        },
                    )?;
                    Ok(OrientationInspectSamples {
                        label: label.clone(),
                        inspect: inspect.inspect,
                        used_arb_prec_retry: inspect.used_arb_prec_retry,
                    })
                })
                .collect::<Result<Vec<_>>>()?;
            let summed = sum_orientation_inspect_samples(&per_orientation);
            InspectSamples {
                summed,
                summed_used_arb_prec_retry: per_orientation
                    .iter()
                    .any(|orientation| orientation.used_arb_prec_retry),
                per_orientation,
            }
        } else {
            let inspect = self.sample_subset_orientation(
                integrand,
                SubsetOrientationInput {
                    graph_id: input.graph_id,
                    subset: input.subset,
                    lmb: input.lmb,
                    generation_lmb: input.generation_lmb,
                    sample: input.sample,
                    orientation: None,
                },
            )?;
            InspectSamples {
                summed: inspect.inspect,
                summed_used_arb_prec_retry: inspect.used_arb_prec_retry,
                per_orientation: Vec::new(),
            }
        };
        let analytic = None;

        Ok(SubSetResult {
            inspect,
            initial_dod: input.initial_dod,
            analytic,
        })
    }

    fn sample_subset_orientation(
        &self,
        integrand: &mut ProcessIntegrand,
        input: SubsetOrientationInput<'_>,
    ) -> Result<InspectRun> {
        let inspect = self.sample_subset_orientation_with_precision(
            integrand,
            input,
            self.profile_settings.use_f128,
        )?;
        if !self.profile_settings.use_f128
            && inspect_results_need_arbprec_retry(&inspect, self.scales)
        {
            return Ok(InspectRun {
                inspect: self.sample_subset_orientation_with_precision(integrand, input, true)?,
                used_arb_prec_retry: true,
            });
        }
        Ok(InspectRun {
            inspect,
            used_arb_prec_retry: false,
        })
    }

    fn sample_subset_orientation_with_precision(
        &self,
        integrand: &mut ProcessIntegrand,
        input: SubsetOrientationInput<'_>,
        use_arb_prec: bool,
    ) -> Result<Vec<InspectResult>> {
        let n_included = input.subset.n_included() as i32;
        self.scales
            .iter()
            .map(|s| {
                let prefactor = s.powi(3 * n_included);
                let mut scaled_sample = input.sample.clone();
                for l in input.subset.included_iter() {
                    scaled_sample[l] = scaled_sample[l].map_ref(&|a| a * F(*s));
                }
                let dependent_momenta_constructor = match &*integrand {
                    ProcessIntegrand::Amplitude(amplitude) => {
                        DependentMomentaConstructor::Amplitude(&amplitude.data.external_signature)
                    }
                    ProcessIntegrand::CrossSection(_) => DependentMomentaConstructor::CrossSection,
                };
                let sample_in_lmb = MomentumSample::new(
                    scaled_sample.iter().cloned().collect::<LoopMomenta<_>>(),
                    0,
                    &self.settings.kinematics.externals,
                    0,
                    F(1.0),
                    dependent_momenta_constructor,
                    input.orientation,
                )?;
                let loop_momenta = input.generation_loop_momenta(&sample_in_lmb);

                let inspect_res_eval = evaluate_momentum_space_point(
                    integrand,
                    self.model,
                    loop_momenta,
                    input.graph_id,
                    input.orientation,
                    use_arb_prec,
                    self.selected_cut_id,
                )?;

                Ok(InspectResult {
                    result: inspect_res_eval,
                    prefactor,
                })
            })
            .collect()
    }
}

impl SubSetResult {
    pub fn analyse_inspect(&self, scales: &[f64]) -> Option<InspectAnalysis> {
        analyse_inspect_results(
            &self.inspect.summed,
            scales,
            self.inspect.summed_used_arb_prec_retry,
        )
    }

    pub fn analyse(&self, scales: &[f64]) -> Analysis {
        Analysis {
            inspect_level: self.analyse_inspect(scales),
            inspect_fit_status: InspectFitStatus::from_results(&self.inspect.summed),
            per_orientation_inspect: (!self.inspect.per_orientation.is_empty()).then(|| {
                self.inspect
                    .per_orientation
                    .iter()
                    .map(|orientation| OrientationInspectAnalysis {
                        orientation_label: orientation.label.clone(),
                        analysis: analyse_inspect_results(
                            &orientation.inspect,
                            scales,
                            orientation.used_arb_prec_retry,
                        ),
                        inspect_fit_status: InspectFitStatus::from_results(&orientation.inspect),
                    })
                    .collect()
            }),
            analytic: self.analyse_analytic(),
        }
    }

    pub fn analyse_analytic(&self) -> Option<AnalyticAnalysis> {
        //             .derived_data
        //             .cff_expression
        //             .as_ref()
        //             .unwrap()
        //             .orientations
        //             .iter()
        //             .enumerate()
        //         {
        //             let expansion = symbol!("lambd");

        //             for (ls, res) in &analytic_res[i_lmb][i] {
        //                 // print!("res:{res}");
        //                 let l = res.coefficient_list::<i8>(&[Atom::var(expansion)]);

        //                 println!(
        //                     "In the limit of {:?} going to infinity for orientation \n{}:",
        //                     ls.included_iter()
        //                         .map(|l| lmb.loop_edges[l].to_string())
        //                         .collecqt::<Vec<_>>(),
        //                     o.data
        //                 );
        //                 if l.is_empty() {
        //                     println!("\tFull cancellation to order 1");
        //                 }
        //                 for (t, a) in l {
        //                     println!("\t{}: {}", t, a);
        //                 }
        //             }
        //         }

        Some(AnalyticAnalysis {
            per_orientation: self
                .analytic
                .as_ref()?
                .per_orientations
                .par_iter()
                .map(|(k, v)| {
                    (
                        k.clone(),
                        OrientationAnalyticAnalysis {
                            is_constant: v.is_constant(),
                            leading_coef: v.lcoeff().to_ordered_simple(),
                        },
                    )
                })
                .collect(),
        })
    }
}

fn analyse_inspect_results(
    inspect: &[InspectResult],
    scales: &[f64],
    used_arb_prec_retry: bool,
) -> Option<InspectAnalysis> {
    let result = log_log_slope(inspect, scales)?;
    let dod = result.slope.round() as i64;
    Some(InspectAnalysis {
        result,
        estimated_dod: dod,
        used_arb_prec_retry,
    })
}

fn inspect_results_need_arbprec_retry(inspect: &[InspectResult], scales: &[f64]) -> bool {
    match analyse_inspect_results(inspect, scales, false) {
        None => true,
        Some(analysis) => {
            analysis.result.slope.is_nan() || analysis.result.slope > UV_PROFILE_RETRY_MAX_DOD
        }
    }
}

fn log_log_slope(inspect: &[InspectResult], scales: &[f64]) -> Option<FitResult> {
    let mut valid_samples = Vec::new();
    let mut sum_x = 0.0;
    let mut sum_y = 0.0;
    let mut sum_xy = 0.0;
    let mut sum_x2 = 0.0;
    let mut points = vec![];
    let mut points_detail = vec![];

    for (x, s) in inspect.iter().zip(scales) {
        let norm = x.magnitude();
        if norm <= 0.0 {
            debug!("{s}:\t{}", x.result.evaluation_metadata);
            continue;
        }
        if !norm.is_finite() {
            continue;
        }
        points_detail.push(x.result.clone());
        points.push(norm);
        let y = (norm).log10();
        let x = s.log10();
        if !y.is_finite() || !x.is_finite() {
            continue;
        }
        valid_samples.push((x, y));
        sum_x += x;
        sum_y += y;
        sum_xy += x * y;
        sum_x2 += x * x;
    }

    if valid_samples.len() < 2 {
        return None;
    }

    let n = valid_samples.len() as f64;
    let denominator = n * sum_x2 - sum_x * sum_x;
    if denominator.abs() < 1e-15 {
        return None;
    }

    let slope = (n * sum_xy - sum_x * sum_y) / denominator;
    let intercept = (sum_y - slope * sum_x) / n;

    let y_mean = sum_y / n;
    let mut ss_tot = 0.0;
    let mut ss_res = 0.0;
    for (x, y) in &valid_samples {
        let y_pred = intercept + slope * x;
        ss_tot += (y - y_mean).powi(2);
        ss_res += (y - y_pred).powi(2);
    }

    let r_squared = if ss_tot > 1e-15 {
        1.0 - ss_res / ss_tot
    } else {
        0.0
    };

    Some(FitResult {
        points,
        points_detail,
        slope,
        intercept,
        r_squared,
    })
}

fn sum_orientation_inspect_samples(
    per_orientation: &[OrientationInspectSamples],
) -> Vec<InspectResult> {
    let Some((first, rest)) = per_orientation.split_first() else {
        return Vec::new();
    };

    (0..first.inspect.len())
        .map(|point_index| {
            let mut summed = first.inspect[point_index].clone();
            for orientation in rest {
                summed.result.integrand_result +=
                    orientation.inspect[point_index].result.integrand_result;
            }
            summed
        })
        .collect()
}

fn inspect_retry_label(analysis: Option<&InspectAnalysis>) -> String {
    match analysis {
        Some(analysis) if analysis.used_arb_prec_retry => "arb retry".yellow().to_string(),
        Some(_) => String::new(),
        None => "-".to_string(),
    }
}

fn evaluate_momentum_space_point(
    integrand: &mut ProcessIntegrand,
    model: &Model,
    loop_momenta: Vec<ThreeMomentum<F<f64>>>,
    graph_id: usize,
    orientation: Option<usize>,
    use_arb_prec: bool,
    selected_cut_id: Option<usize>,
) -> Result<EvaluationResult> {
    let mut result = match integrand {
        ProcessIntegrand::Amplitude(amplitude) => evaluate_profile_momentum_point(
            amplitude,
            model,
            graph_id,
            orientation,
            loop_momenta,
            use_arb_prec,
        ),
        ProcessIntegrand::CrossSection(cross_section) => evaluate_profile_momentum_point(
            cross_section,
            model,
            graph_id,
            orientation,
            loop_momenta,
            use_arb_prec,
        ),
    }?;

    if let Some(selected_cut_id) = selected_cut_id {
        result.integrand_result = result
            .event_groups
            .iter()
            .flat_map(|event_group| event_group.iter())
            .filter(|event| event.cut_info.cut_id == selected_cut_id)
            .fold(Complex::new_re(F(0.0)), |sum, event| sum + event.weight);
    }

    Ok(result)
}

#[derive(Debug, Clone)]
pub struct InspectResult {
    pub(crate) result: EvaluationResult,
    pub(crate) prefactor: f64,
}

impl InspectResult {
    fn magnitude(&self) -> f64 {
        self.result.integrand_result.norm_squared().sqrt().0 * self.prefactor
    }
}

#[derive(Debug, Clone, Copy, Default)]
struct InspectFitStatus {
    finite_samples: usize,
    positive_finite_samples: usize,
}

impl InspectFitStatus {
    fn from_results(inspect: &[InspectResult]) -> Self {
        inspect.iter().fold(Self::default(), |mut status, result| {
            let norm = result.magnitude();
            if norm.is_finite() {
                status.finite_samples += 1;
                if norm > 0.0 {
                    status.positive_finite_samples += 1;
                }
            }
            status
        })
    }

    fn missing_fit_is_vanishing(self) -> bool {
        self.finite_samples > self.positive_finite_samples && self.positive_finite_samples < 2
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct Analysis {
    ///Is None if the fit hasn't worked
    inspect_level: Option<InspectAnalysis>,
    #[serde(skip_serializing)]
    inspect_fit_status: InspectFitStatus,
    #[serde(skip_serializing)]
    per_orientation_inspect: Option<Vec<OrientationInspectAnalysis>>,
    ///Is None if the analytic analysis is disabled
    #[serde(skip_serializing)]
    analytic: Option<AnalyticAnalysis>,
}

impl Analysis {
    fn per_orientation_inspect_entries(&self) -> Option<Vec<UVProfileOrientationInspectEntry>> {
        self.per_orientation_inspect.as_ref().map(|entries| {
            entries
                .iter()
                .map(|entry| UVProfileOrientationInspectEntry {
                    orientation_label: entry.orientation_label.clone(),
                    analysis: entry.analysis.clone(),
                })
                .collect()
        })
    }
}

#[derive(Debug, Clone)]
struct OrientationInspectAnalysis {
    orientation_label: String,
    analysis: Option<InspectAnalysis>,
    inspect_fit_status: InspectFitStatus,
}

#[derive(Debug, Clone, Serialize)]
pub struct AnalyticAnalysis {
    per_orientation: BTreeMap<OrientationData, OrientationAnalyticAnalysis>,
}

#[derive(Debug, Clone, Serialize)]
pub struct OrientationAnalyticAnalysis {
    is_constant: bool,
    leading_coef: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct InspectAnalysis {
    result: FitResult,
    estimated_dod: i64,
    used_arb_prec_retry: bool,
}

pub struct AnalyticResult {
    pub(crate) per_orientations: BTreeMap<OrientationData, Series<AtomField>>,
}
