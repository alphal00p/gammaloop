//! UV Profile Analysis Module
//!
//! This module provides functionality for analyzing ultraviolet behavior of loop integrands
//! by evaluating them at different momentum scalings and computing degrees of divergence.

use std::collections::{BTreeMap, HashSet};
use std::fmt::Display;
use std::path::Path;
use std::sync::{Arc, Mutex};

use crate::DependentMomentaConstructor;
use crate::cff::expression::{OrientationData, OrientationID};
use crate::cff::orientations::GraphOrientation;
use crate::graph::parse::string_utils::ToOrderedSimple;
use crate::graph::{FeynmanGraph, Graph, LmbIndex, LoopMomentumBasis};
use crate::integrands::evaluation::EvaluationResult;
use crate::integrands::process::{
    OrientationProfileMode, ProcessIntegrand, evaluate_profile_momentum_point,
    orientation_labels_for_graph,
};
use crate::model::Model;
use crate::momentum::ThreeMomentum;
use crate::momentum::sample::{LoopIndex, LoopMomenta, MomentumSample};
use crate::processes::{Amplitude, AmplitudeGraph, CrossSection, CutId};
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
use symbolica::atom::Atom;
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
const UV_PROFILE_RETRY_MIN_R_SQUARED: f64 = 0.9;
const UV_PROFILE_ASYMPTOTIC_MIN_R_SQUARED: f64 = 0.99;
const UV_PROFILE_MIN_POINTS: usize = 5;
const UV_PROFILE_THREAD_STACK_SIZE_BYTES: usize = 32 * 1024 * 1024;

#[derive(
    Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, JsonSchema, ValueEnum,
)]
#[serde(rename_all = "kebab-case")]
pub enum UVLimitSelection {
    /// One representative for every expected divergent cycle union: the generation LMB when
    /// suitable, otherwise the first suitable basis in the deterministically sorted full list.
    #[default]
    OnlyDivergent,
    /// Every nonempty loop subset in every generated physical LMB.
    All,
}

impl UVLimitSelection {
    #[cfg(test)]
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
        profile_domain: &SuBitGraph,
        compatible_cuts: Option<&[SuBitGraph]>,
    ) -> Vec<(SubSet<LoopIndex>, i32)> {
        let mut loops = PowersetIterator::<LoopIndex>::new(lmb.loop_edges.len() as u8);
        loops.next();
        loops
            .filter_map(|subset| {
                let cycle = Self::cycle_union(graph, lmb, &subset);
                let is_in_profile_domain = cycle.subtract(profile_domain).is_empty();
                let is_cut_compatible = compatible_cuts
                    .is_none_or(|cuts| cuts.iter().any(|cut| !cycle.intersects(cut)));
                (is_in_profile_domain && is_cut_compatible)
                    .then(|| (subset, graph.compute_dod(&cycle)))
            })
            .collect()
    }

    fn lmb_limits(
        self,
        graph: &Graph,
        all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
        profile_domain: &SuBitGraph,
        compatible_cuts: Option<&[SuBitGraph]>,
    ) -> Result<ProfileLmbLimits> {
        if matches!(self, Self::All) {
            return Ok(all_lmbs
                .iter()
                .cloned()
                .map(|lmb| {
                    let limits = Self::all_lmb_limits(graph, &lmb, profile_domain, compatible_cuts);
                    (lmb, limits)
                })
                .collect());
        }

        let mut sorted_lmbs = all_lmbs.iter().collect_vec();
        sorted_lmbs.sort_by_key(|lmb| {
            (
                lmb.loop_edges.iter().map(|edge| edge.0).collect_vec(),
                lmb.ext_edges.iter().map(|edge| edge.0).collect_vec(),
                lmb.tree.string_label(),
            )
        });
        let mut seen_lmbs = HashSet::new();
        let candidate_lmbs = std::iter::once(&graph.loop_momentum_basis)
            .chain(sorted_lmbs)
            .filter(|lmb| seen_lmbs.insert(*lmb))
            .collect_vec();
        let profile_domains = compatible_cuts
            .map(|cuts| {
                cuts.iter()
                    .map(|cut| profile_domain.subtract(cut))
                    .collect_vec()
            })
            .unwrap_or_else(|| vec![profile_domain.clone()]);
        let mut targets = profile_domains
            .into_iter()
            .filter(|domain| !domain.is_empty())
            .flat_map(|domain| graph.all_cycle_unions(&domain))
            .filter_map(|cycle| {
                let cycle = cycle.filter;
                (!cycle.is_empty()).then(|| {
                    let dod = graph.compute_dod(&cycle);
                    (cycle, dod)
                })
            })
            .filter(|(_, dod)| *dod >= 0)
            .collect_vec();
        targets.sort_by_key(|(cycle, _)| cycle.string_label());
        targets.dedup_by(|(left, _), (right, _)| left == right);
        // Most graphs are fully represented by the generation LMB. Build support maps lazily so
        // stored spanning-tree bases are touched only when an earlier basis cannot represent a
        // divergent cycle union.
        let mut candidate_supports = vec![None; candidate_lmbs.len()];

        let mut selected: ProfileLmbLimits = Vec::new();
        for (target, dod) in targets {
            let mut representation = None;
            for (candidate_index, &lmb) in candidate_lmbs.iter().enumerate() {
                let supports = candidate_supports[candidate_index].get_or_insert_with(|| {
                    (0..lmb.loop_edges.len())
                        .map(|loop_index| {
                            let loop_index = LoopIndex(loop_index);
                            let mut support: SuBitGraph = graph.empty_subgraph();
                            for (pair, edge_id, _) in graph.underlying.iter_edges() {
                                if matches!(pair, HedgePair::Paired { .. })
                                    && lmb.edge_signatures[edge_id].internal[loop_index]
                                        != SignOrZero::Zero
                                {
                                    support.add(pair);
                                }
                            }
                            support
                        })
                        .collect_vec()
                });
                let mut subset = SubSet::empty(lmb.loop_edges.len());
                let mut represented_cycle: SuBitGraph = graph.empty_subgraph();
                for (loop_index, support) in supports.iter().enumerate() {
                    let loop_index = LoopIndex(loop_index);
                    if !support.is_empty() && support.subtract(&target).is_empty() {
                        subset.add(loop_index);
                        represented_cycle.union_with(support);
                    }
                }

                if represented_cycle == target {
                    representation = Some((lmb, subset));
                    break;
                }
            }

            let Some((lmb, subset)) = representation else {
                return Err(eyre!(
                    "Divergent UV cycle {} (DOD {dod}) in graph '{}' is not represented by the generation LMB or any basis in the complete generated LMB list",
                    target.string_label(),
                    graph.name,
                ));
            };
            if let Some((_, limits)) = selected
                .iter_mut()
                .find(|(selected_lmb, _)| selected_lmb == lmb)
            {
                limits.push((subset, dod));
            } else {
                selected.push((lmb.clone(), vec![(subset, dod)]));
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

    use itertools::Itertools;
    use linnet::half_edge::involution::{EdgeIndex, EdgeVec, Orientation};
    use linnet::half_edge::subgraph::{
        Inclusion, ModifySubSet, SuBitGraph, SubSetLike, subset::SubSet,
    };
    use typed_index_collections::TiVec;

    use crate::dot;
    use crate::graph::parse::from_dot::IntoGraph;
    use crate::graph::{FeynmanGraph, Graph, LMBext, LmbIndex, LoopMomentumBasis};
    use crate::initialisation::test_initialise;
    use crate::integrands::evaluation::EvaluationResult;
    use crate::momentum::ThreeMomentum;
    use crate::momentum::sample::{
        BareMomentumSample, ExternalThreeMomenta, LoopIndex, LoopMomenta, MomentumSample,
    };
    use crate::observables::events::{AdditionalWeightKey, Event};
    use crate::utils::F;
    use crate::uv::UltravioletGraph;
    use spenso::algebra::complex::Complex;
    use symbolica::atom::{Atom, AtomCore};
    use three_dimensional_reps::symbols::sign;

    use super::{
        Analysis, FitResult, InspectAnalysis, InspectFitStatus, InspectResult,
        OrientationInspectAnalysis, SubsetOrientationInput, UV_PROFILE_ASYMPTOTIC_MIN_R_SQUARED,
        UV_PROFILE_RETRY_MIN_R_SQUARED, UVLimitSelection, UVProfile, UVProfileAnalysis,
        UVProfileFailure, UVProfileFixedRay, UVProfileGraphAnalysis, UVProfileLmbAnalysis,
        UVProfilePassFail, UVProfileSubsetAnalysis, UVSamplingResult,
        analytic_integrand_for_orientation, inspect_results_need_arbprec_retry, log_log_slope,
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
    fn analytic_integrand_selects_physical_orientation_and_exact_residue_map_key() {
        let edge = EdgeIndex(0);
        let mut orientation =
            super::OrientationData::new(EdgeVec::from_iter([Orientation::Default]));
        orientation.label = Some("+|N0".to_string());
        let integrand =
            sign(edge) * (super::OrientationID(4).atom() * 2 + super::OrientationID(9).atom() * 3);

        let (four_data, four) =
            analytic_integrand_for_orientation(super::OrientationID(4), &orientation, &integrand);
        let (nine_data, nine) =
            analytic_integrand_for_orientation(super::OrientationID(9), &orientation, &integrand);

        assert_eq!(four_data.label.as_deref(), Some("+|N0|sigma(4)"));
        assert_eq!(nine_data.label.as_deref(), Some("+|N0|sigma(9)"));
        assert_eq!(four.expand(), Atom::num(2));
        assert_eq!(nine.expand(), Atom::num(3));
        assert!(!four.contains_symbol(super::OrientationID::symbol()));
        assert!(!nine.contains_symbol(super::OrientationID::symbol()));
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
    fn only_divergent_targets_stay_inside_the_production_uv_domain() {
        let _ = theta_graph();
        let mut graph: Graph = dot!(digraph theta_profile_domain {
            edge [num=1 mass=0]
            node [num=1]
            A -> B [id=0]
            A -> B [id=1]
            A -> B [id=2 is_dummy=true]
        })
        .unwrap();
        let lmbs = graph.generate_loop_momentum_bases();
        graph.loop_momentum_basis = lmbs[LmbIndex::from(0)].clone();
        let selected = UVLimitSelection::OnlyDivergent
            .lmb_limits(&graph, &lmbs, &graph.no_dummy(), None)
            .unwrap();
        assert!(!selected.is_empty());
        assert!(
            selected
                .iter()
                .all(|(lmb, limits)| limits.iter().all(|(subset, _)| {
                    !UVLimitSelection::cycle_union(&graph, lmb, subset)
                        .includes(&graph[&EdgeIndex::from(2)].1)
                }))
        );

        let mut incoming_domain = graph.full_filter();
        incoming_domain.sub(graph[&EdgeIndex::from(0)].1);
        let selected = UVLimitSelection::OnlyDivergent
            .lmb_limits(&graph, &lmbs, &incoming_domain, None)
            .unwrap();
        assert!(!selected.is_empty());
        assert!(
            selected
                .iter()
                .all(|(lmb, limits)| limits.iter().all(|(subset, _)| {
                    !UVLimitSelection::cycle_union(&graph, lmb, subset)
                        .includes(&graph[&EdgeIndex::from(0)].1)
                }))
        );
    }

    #[test]
    fn only_divergent_profiles_all_theta_cycles_once() {
        let (graph, lmbs) = theta_graph();
        let selected = UVLimitSelection::OnlyDivergent
            .lmb_limits(&graph, &lmbs, &graph.full_filter(), None)
            .unwrap();
        let represented = selected
            .iter()
            .flat_map(|(lmb, limits)| {
                limits.iter().map(|(subset, _)| {
                    UVLimitSelection::cycle_union(&graph, lmb, subset).string_label()
                })
            })
            .collect::<std::collections::BTreeSet<_>>();
        let expected = graph
            .all_cycle_unions(&graph.full_filter())
            .into_iter()
            .filter_map(|cycle| {
                let cycle = cycle.filter;
                (!cycle.is_empty() && graph.compute_dod(&cycle) >= 0).then(|| cycle.string_label())
            })
            .collect::<std::collections::BTreeSet<_>>();
        let selected_limit_count: usize = selected.iter().map(|(_, limits)| limits.len()).sum();

        assert_eq!(represented, expected);
        assert_eq!(selected_limit_count, 4);
        assert_eq!(represented.len(), selected_limit_count);
        assert!(
            selected.len() > 1,
            "the third one-loop theta cycle requires another generated physical basis"
        );
        let generation_limits = selected
            .iter()
            .find(|(lmb, _)| lmb == &graph.loop_momentum_basis)
            .map(|(_, limits)| limits)
            .expect("the generation LMB must represent every suitable divergent theta limit");
        assert_eq!(
            generation_limits.len(),
            3,
            "the generation LMB must be preferred for its two one-loop cycles and their union"
        );
        let (fallback_lmb, fallback_limits) = selected
            .iter()
            .find(|(lmb, _)| lmb != &graph.loop_momentum_basis)
            .expect("the third one-loop theta cycle needs a fallback LMB");
        assert_eq!(fallback_limits.len(), 1);
        let fallback_target =
            UVLimitSelection::cycle_union(&graph, fallback_lmb, &fallback_limits[0].0);
        let first_suitable_fallback = lmbs
            .iter()
            .sorted_by_key(|lmb| {
                (
                    lmb.loop_edges.iter().map(|edge| edge.0).collect_vec(),
                    lmb.ext_edges.iter().map(|edge| edge.0).collect_vec(),
                    lmb.tree.string_label(),
                )
            })
            .find(|lmb| {
                UVLimitSelection::all_lmb_limits(&graph, lmb, &graph.full_filter(), None)
                    .iter()
                    .any(|(subset, _)| {
                        UVLimitSelection::cycle_union(&graph, lmb, subset) == fallback_target
                    })
            })
            .expect("a generated physical LMB must represent the fallback theta cycle");
        assert_eq!(
            fallback_lmb, first_suitable_fallback,
            "the first suitable deterministically sorted LMB must represent a fallback limit"
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

        let reversed_lmbs: TiVec<LmbIndex, LoopMomentumBasis> =
            lmbs.iter().cloned().rev().collect::<Vec<_>>().into();
        let reversed_selected = UVLimitSelection::OnlyDivergent
            .lmb_limits(&graph, &reversed_lmbs, &graph.full_filter(), None)
            .unwrap();
        assert_eq!(
            reversed_selected, selected,
            "fallback LMB selection must be independent of the generated-list order"
        );

        let exhaustive = UVLimitSelection::All
            .lmb_limits(&graph, &lmbs, &graph.full_filter(), None)
            .unwrap();
        assert_eq!(exhaustive.len(), lmbs.len());
        assert!(exhaustive.iter().all(
            |(_, limits)| limits.len() == (1 << graph.loop_momentum_basis.loop_edges.len()) - 1
        ));

        let mut first_physical_cut: SuBitGraph = graph.empty_subgraph();
        first_physical_cut.add(graph[&EdgeIndex::from(0)].1);
        let mut second_physical_cut: SuBitGraph = graph.empty_subgraph();
        second_physical_cut.add(graph[&EdgeIndex::from(1)].1);
        let physical_cuts = [first_physical_cut.clone(), second_physical_cut.clone()];
        let cut_compatible_exhaustive = UVLimitSelection::All
            .lmb_limits(&graph, &lmbs, &graph.full_filter(), Some(&physical_cuts))
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
                    let cycle = UVLimitSelection::cycle_union(&graph, lmb, subset);
                    physical_cuts.iter().any(|cut| !cycle.intersects(cut))
                }))
        );
        assert!(
            cut_compatible_exhaustive.iter().any(|(lmb, limits)| {
                limits.iter().any(|(subset, _)| {
                    let cycle = UVLimitSelection::cycle_union(&graph, lmb, subset);
                    cycle.intersects(&first_physical_cut) && !cycle.intersects(&second_physical_cut)
                })
            }),
            "a limit compatible with one exact physical cut must not be rejected because it intersects another cut in the same residue group"
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
            compatible_event_cut_ids: None,
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
    fn pass_fail_display_lists_multiple_failures_with_context() {
        let report = UVProfilePassFail {
            max_dod: -0.9,
            total: 4,
            failed: 2,
            failures: vec![
                UVProfileFailure {
                    graph_index: 2,
                    graph_name: "GL2".to_string(),
                    cutkosky_cut: Some(vec![EdgeIndex::from(5), EdgeIndex::from(2)]),
                    lmb_index: 0,
                    fixed: vec![EdgeIndex::from(1)],
                    free: vec![EdgeIndex::from(3)],
                    orientation_label: None,
                    reason: "estimated DOD exceeds threshold".to_string(),
                },
                UVProfileFailure {
                    graph_index: 0,
                    graph_name: "GL0".to_string(),
                    cutkosky_cut: None,
                    lmb_index: 1,
                    fixed: vec![EdgeIndex::from(4)],
                    free: vec![EdgeIndex::from(6), EdgeIndex::from(7)],
                    orientation_label: Some("+--+".to_string()),
                    reason: "missing fit".to_string(),
                },
            ],
        };

        let rendered = report.to_string();
        assert!(rendered.contains("UV limit tests:"));
        assert!(rendered.contains("2/4"));
        assert!(rendered.contains("GL2"));
        assert!(rendered.contains("e5,e2"));
        assert!(rendered.contains("estimated DOD exceeds threshold"));
        assert!(rendered.contains("GL0"));
        assert!(rendered.contains("+--+"));
        assert!(rendered.contains("missing fit"));
        assert_eq!(rendered.matches("FAIL").count(), 3);
    }

    #[test]
    fn passing_orientations_do_not_hide_a_failing_summed_limit() {
        let fit = |slope, estimated_dod| InspectAnalysis {
            result: FitResult {
                slope,
                points: Vec::new(),
                scales: Vec::new(),
                fit_start: 0,
                intercept: 0.0,
                r_squared: 1.0,
                full_range_slope: slope,
                full_range_intercept: 0.0,
                full_range_r_squared: 1.0,
            },
            estimated_dod,
            used_arb_prec_retry: false,
        };
        let analysis = UVProfileAnalysis {
            scales: Vec::new(),
            graphs: vec![UVProfileGraphAnalysis {
                graph_index: 0,
                graph_name: "GL0".to_string(),
                cutkosky_cut: None,
                lmbs: vec![UVProfileLmbAnalysis {
                    lmb_index: 0,
                    lmb_label: "generation LMB".to_string(),
                    subsets: vec![UVProfileSubsetAnalysis {
                        subset_index: 0,
                        fixed: Vec::new(),
                        free: vec![EdgeIndex::from(1)],
                        initial_dod: 0,
                        analysis: Analysis {
                            inspect_level: Some(fit(0.0, 0)),
                            inspect_fit_status: InspectFitStatus::default(),
                            per_orientation_inspect: Some(vec![OrientationInspectAnalysis {
                                orientation_label: "+-".to_string(),
                                analysis: Some(fit(-1.0, -1)),
                                inspect_fit_status: InspectFitStatus::default(),
                            }]),
                            analytic: None,
                        },
                        per_orientation_inspect_entries: None,
                        analytic_entries: None,
                    }],
                }],
            }],
            allow_vanishing_missing_fits: true,
        };

        let report = analysis.pass_fail(-0.9);
        assert_eq!(report.total, 2);
        assert_eq!(report.failed, 1);
        assert_eq!(report.failures.len(), 1);
        assert_eq!(report.failures[0].orientation_label, None);
        assert_eq!(report.failures[0].reason, "dod_exceeds_threshold");
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

    #[test]
    fn uv_fit_retries_low_quality_double_results_before_using_an_asymptotic_tail() {
        let scales = (0..17)
            .map(|i| 10.0_f64.powf(4.0 + 0.25 * i as f64))
            .collect_vec();
        let double_with_interior_outlier = scales
            .iter()
            .enumerate()
            .map(|(i, scale)| {
                let mut result = EvaluationResult::zero();
                let magnitude = if i == 7 {
                    scale.recip() * 1.0e10
                } else {
                    scale.recip()
                };
                result.integrand_result = Complex::new_re(F(magnitude));
                InspectResult {
                    result,
                    prefactor: 1.0,
                }
            })
            .collect_vec();

        let fit = log_log_slope(&double_with_interior_outlier, &scales).unwrap();
        assert!((fit.slope + 1.0).abs() < 1.0e-12);
        assert!(fit.full_range_r_squared < 0.9);
        assert!(inspect_results_need_arbprec_retry(
            &double_with_interior_outlier,
            &scales
        ));

        // A fit can clear the coarse full-range retry gate while still lacking a trustworthy
        // asymptotic suffix. Retry that case as well instead of reporting the Double fit as an
        // unstable physical limit.
        let double_with_distributed_noise = scales
            .iter()
            .enumerate()
            .map(|(i, scale)| {
                let mut result = EvaluationResult::zero();
                let log_noise = if i % 2 == 0 { 0.2 } else { -0.2 };
                result.integrand_result =
                    Complex::new_re(F(scale.recip() * 10.0_f64.powf(log_noise)));
                InspectResult {
                    result,
                    prefactor: 1.0,
                }
            })
            .collect_vec();
        let fit = log_log_slope(&double_with_distributed_noise, &scales).unwrap();
        assert!(fit.full_range_r_squared > UV_PROFILE_RETRY_MIN_R_SQUARED);
        assert!(fit.r_squared < UV_PROFILE_ASYMPTOTIC_MIN_R_SQUARED);
        assert!(inspect_results_need_arbprec_retry(
            &double_with_distributed_noise,
            &scales
        ));

        // Arb evaluation of an actual banana ray has one pre-asymptotic point followed by a
        // clean asymptotic suffix. These logarithms are the retained regression data from that
        // ray; no precision retry is required once the full-precision samples are available.
        let arb_log_magnitudes = [
            -25.361635842,
            -24.003219677,
            -24.055611693,
            -24.224420479,
            -24.434654742,
            -24.663795182,
            -24.902490667,
            -25.146260599,
            -25.392796041,
            -25.640859851,
            -25.889774837,
            -26.139165860,
            -26.388823797,
            -26.638631544,
            -26.888523551,
            -27.138462785,
            -27.388428677,
        ];
        let arb_transient = arb_log_magnitudes
            .into_iter()
            .map(|log_magnitude| {
                let mut result = EvaluationResult::zero();
                result.integrand_result = Complex::new_re(F(10.0_f64.powf(log_magnitude)));
                InspectResult {
                    result,
                    prefactor: 1.0,
                }
            })
            .collect_vec();
        let fit = log_log_slope(&arb_transient, &scales).unwrap();
        assert!((fit.slope + 0.945490).abs() < 1.0e-5);
        assert!(fit.r_squared > 0.99);
        assert!(fit.full_range_r_squared < 0.9);
    }

    #[test]
    fn uv_fit_json_excludes_raw_lu_evaluation_results() {
        let scales = [1.0e3, 1.0e4, 1.0e5, 1.0e6, 1.0e7];
        let inspect = scales
            .iter()
            .map(|scale| {
                let mut result = EvaluationResult::zero();
                result.integrand_result = Complex::new_re(F(1.0 / scale));
                let mut event = Event::default();
                event.additional_weights.weights.insert(
                    AdditionalWeightKey::ThresholdCounterterm { subset_index: 0 },
                    Complex::new_re(F(1.0)),
                );
                result.event_groups.push_singleton(event);
                InspectResult {
                    result,
                    prefactor: 1.0,
                }
            })
            .collect_vec();

        let analysis = super::analyse_inspect_results(&inspect, &scales, false).unwrap();
        let json = serde_json::to_value(analysis)
            .expect("UV fit JSON must not include raw LU evaluations with non-string map keys");
        assert_eq!(
            json["result"]["points"].as_array().unwrap().len(),
            scales.len()
        );
        assert_eq!(json["result"]["scales"], serde_json::json!(scales));
        assert_eq!(json["result"]["fit_start"], 0);
        assert!(json["result"].get("full_range_slope").is_some());
        assert!(json["result"].get("full_range_intercept").is_some());
        assert!(json["result"].get("full_range_r_squared").is_some());
        assert!(json["result"].get("points_detail").is_none());
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
        if profile_settings.n_points < UV_PROFILE_MIN_POINTS {
            return Err(eyre!(
                "UV profiling needs at least {UV_PROFILE_MIN_POINTS} scale points for a reliable asymptotic fit"
            ));
        }
        if !profile_settings.min_scale_exponent.is_finite()
            || !profile_settings.max_scale_exponent.is_finite()
            || profile_settings.min_scale_exponent >= profile_settings.max_scale_exponent
        {
            return Err(eyre!(
                "UV profiling needs finite scale exponents with min < max, got {} and {}",
                profile_settings.min_scale_exponent,
                profile_settings.max_scale_exponent,
            ));
        }
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
        // Profiling owns a private evaluator clone, so an evaluation error cannot leave the
        // process without its production integrand.
        let integrand = Arc::new(Mutex::new(self.integrand.as_ref().unwrap().clone()));

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
        };

        // Generated symbolic evaluators need the same worker-stack headroom while profiling as
        // they do in the production generation pool.
        let profile_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(rayon::current_num_threads())
            .stack_size(UV_PROFILE_THREAD_STACK_SIZE_BYTES)
            .build()?;
        let per_graph = profile_pool.install(|| {
            graph_inputs
                .par_iter()
                .map(|(graph_id, graph)| {
                    let res = runner.sample_graph(*graph_id, graph)?;
                    profile_span.pb_inc(1);
                    Ok(res)
                })
                .collect::<Result<Vec<_>>>()
        })?;

        drop(_profile_span_enter);
        drop(profile_span);

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
        if profile_settings.n_points < UV_PROFILE_MIN_POINTS {
            return Err(eyre!(
                "UV profiling needs at least {UV_PROFILE_MIN_POINTS} scale points for a reliable asymptotic fit"
            ));
        }
        if !profile_settings.min_scale_exponent.is_finite()
            || !profile_settings.max_scale_exponent.is_finite()
            || profile_settings.min_scale_exponent >= profile_settings.max_scale_exponent
        {
            return Err(eyre!(
                "UV profiling needs finite scale exponents with min < max, got {} and {}",
                profile_settings.min_scale_exponent,
                profile_settings.max_scale_exponent,
            ));
        }
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
        let graph_inputs = match integrand {
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

                let selected_cut = profile_settings
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
                        Ok(cut_id)
                    })
                    .transpose()?;

                cross_section
                    .data
                    .graph_terms
                    .iter()
                    .enumerate()
                    .filter(|(graph_id, _)| {
                        profile_settings.graph_id.is_none_or(|id| id == *graph_id)
                    })
                    .map(|(graph_id, graph_term)| -> Result<_> {
                        // LU evaluates one contribution per cut group and tags its generated
                        // event with the group's first (representative) physical cut. Retain the
                        // exact physical support separately because UV-cycle compatibility is not
                        // a property of the union of a group's raised-edge aliases.
                        let cuts = graph_term
                            .cuts
                            .iter_enumerated()
                            .filter(|(cut_id, _)| {
                                selected_cut.is_none_or(|selected| selected == *cut_id)
                            })
                            .map(|(cut_id, cut)| {
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
                                let representative_cut_id =
                                    cut_group.cuts.first().copied().ok_or_else(|| {
                                        eyre!(
                                            "The cut group containing cut {} is empty for graph '{}'",
                                            cut_id.0,
                                            graph_term.graph.name
                                        )
                                    })?;
                                Ok((cut.cut.as_subgraph(), representative_cut_id))
                            })
                            .collect::<Result<Vec<_>>>()?;
                        Ok((
                            graph_id,
                            graph_term.graph.clone(),
                            graph_term.lmbs.clone(),
                            cuts,
                        ))
                    })
                    .collect::<Result<Vec<_>>>()?
            }
            ProcessIntegrand::Amplitude(_) => {
                unreachable!("cross-section UV profiling expects cross-section integrands")
            }
        };
        let base_seed = profile_settings.seed;
        // Profiling owns a private evaluator clone, so an evaluation error cannot leave the
        // process without its production integrand.
        let integrand = Arc::new(Mutex::new(self.integrand.as_ref().unwrap().clone()));

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
        };

        // Generated symbolic evaluators need the same worker-stack headroom while profiling as
        // they do in the production generation pool.
        let profile_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(rayon::current_num_threads())
            .stack_size(UV_PROFILE_THREAD_STACK_SIZE_BYTES)
            .build()?;
        let per_graph = profile_pool.install(|| {
            graph_inputs
                .par_iter()
                .map(|(graph_id, graph, lmbs, cuts)| {
                    let res = runner.sample_cross_section_graph(*graph_id, graph, lmbs, cuts)?;
                    profile_span.pb_inc(1);
                    Ok(res)
                })
                .collect::<Result<Vec<_>>>()
        })?;

        drop(_profile_span_enter);
        drop(profile_span);

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
                                                    orientation_label: orientation
                                                        .label
                                                        .clone()
                                                        .unwrap_or_else(|| orientation.to_string()),
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
    pub orientation_label: String,
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
    scales: Vec<f64>,
    fit_start: usize,
    intercept: f64,
    r_squared: f64,
    full_range_slope: f64,
    full_range_intercept: f64,
    full_range_r_squared: f64,
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
                let mut header = vec![
                    "fixed".to_string(),
                    "→ ∞".to_string(),
                    "orientation key".to_string(),
                ];
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
                        let mut row = vec![fixed, free, entry.orientation_label.clone()];
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
                    total += 1;
                    let reason = inspect_failure_reason(
                        subset.analysis.inspect_level.as_ref(),
                        subset.analysis.inspect_fit_status,
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
        Some(analysis)
            if !analysis.result.r_squared.is_finite()
                || analysis.result.r_squared < UV_PROFILE_ASYMPTOTIC_MIN_R_SQUARED =>
        {
            Some("unstable_fit")
        }
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

fn analytic_integrand_for_orientation(
    orientation_id: OrientationID,
    orientation: &OrientationData,
    integrand: &Atom,
) -> (OrientationData, Atom) {
    let mut orientation = orientation.clone();
    orientation.label = Some(match orientation.label {
        Some(label) => format!("{label}|sigma({})", orientation_id.0),
        None => format!("sigma({})", orientation_id.0),
    });
    // The full residue-map key is the primary selector. Drop inactive keyed
    // branches before resolving physical-direction metadata, which can be
    // singular outside the branch that owns it.
    let selected = orientation
        .orientation
        .select(orientation_id.select(integrand));
    (orientation, selected)
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
                .iter_enumerated()
                .map(|(orientation_id, orientation)| {
                    analytic_integrand_for_orientation(
                        orientation_id,
                        &orientation.data,
                        &g.derived_data.all_mighty_integrand,
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
        let lmbs = self.profile_settings.selected_limits.lmb_limits(
            &g.graph,
            all_lmbs,
            &g.graph.no_dummy(),
            None,
        )?;

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
                    None,
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
        cuts: &[(SuBitGraph, CutId)],
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
        let mut profile_domain = graph.full_filter();
        profile_domain.subtract_with(&graph.initial_state_cut.left);
        let cut_supports = cuts.iter().map(|(cut, _)| cut.clone()).collect_vec();
        let lmbs = self.profile_settings.selected_limits.lmb_limits(
            graph,
            all_lmbs,
            &profile_domain,
            Some(&cut_supports),
        )?;

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
                    Some(cuts),
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
    #[allow(clippy::too_many_arguments)]
    fn sample_lmb(
        &self,
        graph_id: usize,
        graph: &Graph,
        lmb_index: usize,
        lmb: &LoopMomentumBasis,
        subsets: &[(SubSet<LoopIndex>, i32)],
        orientation_labels: Option<&[String]>,
        compatible_cuts: Option<&[(SuBitGraph, CutId)]>,
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
                    // LU cut compatibility is projected from the production event weights.
                    if compatible_cuts.is_some() {
                        integrand.get_mut_settings().general.generate_events = true;
                    }
                    integrand
                },
                |integrand, (subset, initial_dod)| {
                    let compatible_event_cut_ids = compatible_cuts.map(|cuts| {
                        let cycle = UVLimitSelection::cycle_union(graph, lmb, subset);
                        cuts
                            .iter()
                            .filter(|(cut, _)| !cycle.intersects(cut))
                            .map(|(_, representative_cut_id)| *representative_cut_id)
                            .sorted_by_key(|cut_id| cut_id.0)
                            .dedup()
                            .collect_vec()
                    });
                    if compatible_event_cut_ids
                        .as_ref()
                        .is_some_and(Vec::is_empty)
                    {
                        return Err(eyre!(
                            "UV subset {subset:?} in LMB '{}' of graph '{}' has no compatible Cutkosky cut after cut-compatible limit selection",
                            lmb.tree.string_label(),
                            graph.name,
                        ));
                    }
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
                            compatible_event_cut_ids: compatible_event_cut_ids.as_deref(),
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
    compatible_event_cut_ids: Option<&'a [CutId]>,
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
    compatible_event_cut_ids: Option<&'a [CutId]>,
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
                            compatible_event_cut_ids: input.compatible_event_cut_ids,
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
                    compatible_event_cut_ids: input.compatible_event_cut_ids,
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
                    input.compatible_event_cut_ids,
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
            analysis.result.slope.is_nan()
                || analysis.result.slope > UV_PROFILE_RETRY_MAX_DOD
                || !analysis.result.r_squared.is_finite()
                || analysis.result.r_squared < UV_PROFILE_ASYMPTOTIC_MIN_R_SQUARED
                || analysis.result.full_range_slope.is_nan()
                || analysis.result.full_range_slope > UV_PROFILE_RETRY_MAX_DOD
                || !analysis.result.full_range_r_squared.is_finite()
                || analysis.result.full_range_r_squared < UV_PROFILE_RETRY_MIN_R_SQUARED
        }
    }
}

fn log_log_slope(inspect: &[InspectResult], scales: &[f64]) -> Option<FitResult> {
    let mut valid_samples = Vec::new();
    let mut points = vec![];
    let mut valid_scales = vec![];

    for (x, s) in inspect.iter().zip(scales) {
        let norm = x.magnitude();
        if norm <= 0.0 {
            debug!("{s}:\t{}", x.result.evaluation_metadata);
            continue;
        }
        if !norm.is_finite() {
            continue;
        }
        let y = (norm).log10();
        let x = s.log10();
        if !y.is_finite() || !x.is_finite() {
            continue;
        }
        points.push(norm);
        valid_scales.push(*s);
        valid_samples.push((x, y));
    }

    if valid_samples.len() < 2 {
        return None;
    }

    let fit = |samples: &[(f64, f64)]| {
        let n = samples.len() as f64;
        let (sum_x, sum_y, sum_xy, sum_x2) = samples.iter().fold(
            (0.0, 0.0, 0.0, 0.0),
            |(sum_x, sum_y, sum_xy, sum_x2), (x, y)| {
                (sum_x + x, sum_y + y, sum_xy + x * y, sum_x2 + x * x)
            },
        );
        let denominator = n * sum_x2 - sum_x * sum_x;
        if denominator.abs() < 1e-15 {
            return None;
        }
        let slope = (n * sum_xy - sum_x * sum_y) / denominator;
        let intercept = (sum_y - slope * sum_x) / n;
        let y_mean = sum_y / n;
        let (ss_tot, ss_res) = samples.iter().fold((0.0, 0.0), |(ss_tot, ss_res), (x, y)| {
            let residual = y - (intercept + slope * x);
            (ss_tot + (y - y_mean).powi(2), ss_res + residual.powi(2))
        });
        let r_squared = if ss_tot > 1e-15 {
            1.0 - ss_res / ss_tot
        } else {
            0.0
        };
        Some((slope, intercept, r_squared))
    };

    let full_range_fit = fit(&valid_samples)?;
    // UV power counting is asymptotic. Retain the longest high-scale suffix which is already
    // linear, dropping only a pre-asymptotic prefix. At least five points are kept whenever the
    // sampling range provides them.
    let minimum_suffix_len = valid_samples
        .len()
        .min(5.max(valid_samples.len().div_ceil(2)));
    let (fit_start, (slope, intercept, r_squared)) = (0..=valid_samples.len() - minimum_suffix_len)
        .find_map(|start| {
            let candidate = fit(&valid_samples[start..])?;
            (candidate.2 >= UV_PROFILE_ASYMPTOTIC_MIN_R_SQUARED).then_some((start, candidate))
        })
        .unwrap_or((0, full_range_fit));
    if slope > UV_PROFILE_RETRY_MAX_DOD || r_squared < UV_PROFILE_ASYMPTOTIC_MIN_R_SQUARED {
        debug!(
            slope,
            intercept,
            r_squared,
            full_range_r_squared = full_range_fit.2,
            fit_start,
            ?valid_samples,
            "unstable UV profile log-log fit"
        );
    }

    Some(FitResult {
        points,
        scales: valid_scales,
        fit_start,
        slope,
        intercept,
        r_squared,
        full_range_slope: full_range_fit.0,
        full_range_intercept: full_range_fit.1,
        full_range_r_squared: full_range_fit.2,
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
    compatible_event_cut_ids: Option<&[CutId]>,
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

    if let Some(compatible_event_cut_ids) = compatible_event_cut_ids {
        result.integrand_result = result
            .event_groups
            .iter()
            .flat_map(|event_group| event_group.iter())
            .filter(|event| compatible_event_cut_ids.contains(&CutId(event.cut_info.cut_id)))
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
