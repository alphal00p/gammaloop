use std::collections::BTreeMap;

use crate::{
    cff::{CutCFFIndex, expression::GraphOrientation},
    debug_tags,
    graph::{Graph, cuts::CutSet},
    momentum::Sign,
    numerator::symbolica_ext::AtomCoreExt,
    settings::global::OrientationPattern,
    utils::{GS, W_},
    uv::{
        IntegrandExpr, UltravioletGraph,
        approx::{ApproxOp, ForestNodeLike},
    },
};
use color_eyre::Result;
use eyre::eyre;
use gammaloop_tracing_filter::{LogMessage, debug_instrument};
use idenso::{color::ColorSimplifier, shorthands::metric::MetricSimplifier};
use linnet::half_edge::{
    involution::{EdgeIndex, EdgeVec, Orientation},
    subgraph::{SuBitGraph, SubSetOps},
};
use symbolica::{
    atom::{Atom, AtomCore},
    function, symbol,
};

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct RepresentativeScore {
    undirected_count: usize,
    default_count: usize,
    leading_defaults: Vec<bool>,
}

pub(crate) struct FinalIntegrand<'a> {
    valid_orientations: &'a [EdgeVec<Orientation>],
    orientation_pattern: &'a OrientationPattern,
    uv_marker: Option<&'a Atom>,
    tag_local_terms: bool,
}

pub(crate) struct LocalizedIntegratedCt {
    pub active: IntegrandExpr,
    pub frozen_integrands: BTreeMap<CutCFFIndex, Atom>,
}

impl TryFrom<LocalizedIntegratedCt> for IntegrandExpr {
    type Error = eyre::Report;

    fn try_from(value: LocalizedIntegratedCt) -> Result<Self, Self::Error> {
        let integrands = value
            .active
            .integrands
            .into_iter()
            .zip(value.frozen_integrands)
            .map(|((active_index, active), (frozen_index, frozen))| {
                if active_index != frozen_index {
                    return Err(eyre!(
                        "Mismatched indices for localized integrated CT active and frozen factors: {:?} vs {:?}",
                        active_index,
                        frozen_index
                    ));
                }
                Ok((active_index, active * frozen))
            })
            .collect::<Result<BTreeMap<_, _>>>()?;

        Ok(IntegrandExpr { integrands })
    }
}

impl<'a> FinalIntegrand<'a> {
    pub(crate) fn new(
        valid_orientations: &'a [EdgeVec<Orientation>],
        orientation_pattern: &'a OrientationPattern,
        uv_marker: Option<&'a Atom>,
        tag_local_terms: bool,
    ) -> Self {
        Self {
            valid_orientations,
            orientation_pattern,
            uv_marker,
            tag_local_terms,
        }
    }

    fn representative_score(
        orientation: &EdgeVec<Orientation>,
        internal_edges: &[EdgeIndex],
    ) -> RepresentativeScore {
        let mut undirected_count = 0;
        let mut default_count = 0;
        let mut leading_defaults = Vec::with_capacity(internal_edges.len());

        for edge in internal_edges {
            let is_default = match orientation[*edge] {
                Orientation::Undirected => {
                    undirected_count += 1;
                    false
                }
                Orientation::Default => {
                    default_count += 1;
                    true
                }
                Orientation::Reversed => false,
            };
            leading_defaults.push(is_default);
        }

        RepresentativeScore {
            undirected_count,
            default_count,
            leading_defaults,
        }
    }

    fn orientations_match_outside_integrated_subgraph(
        reduced_orientation: &EdgeVec<Orientation>,
        global_orientation: &EdgeVec<Orientation>,
    ) -> bool {
        reduced_orientation.iter().all(|(edge, reduced)| {
            matches!(reduced, Orientation::Undirected) || *reduced == global_orientation[edge]
        })
    }

    /// Select the deterministic full-graph representative used to host the integrated UV term.
    ///
    /// The reduced orientation already fixes all edges that remain explicit after contracting the
    /// integrated subgraph. Any contracted edge appears as `Undirected` and is therefore ignored when
    /// matching compatible full orientations.
    ///
    /// Among the compatible candidates, the representative is chosen by maximizing:
    /// 1. the number of `Undirected` internal edges
    /// 2. then the number of `Default` internal edges
    /// 3. then the lexicographically maximal `is_default` pattern in ascending internal-edge order
    ///
    /// The first criterion is currently future-proof only: the canonical full-graph acyclic basis is
    /// built from `Default` / `Reversed` assignments on paired internal edges, so `Undirected`
    /// internal entries are not expected there unless the global orientation-generation step changes.
    fn select_representative_orientation<'b>(
        reduced_orientation: &EdgeVec<Orientation>,
        valid_global_orientations: &'b [EdgeVec<Orientation>],
        internal_edges: &[EdgeIndex],
    ) -> Result<&'b EdgeVec<Orientation>> {
        valid_global_orientations
            .iter()
            .filter(|candidate| {
                Self::orientations_match_outside_integrated_subgraph(reduced_orientation, candidate)
            })
            .max_by_key(|candidate| Self::representative_score(candidate, internal_edges))
            .ok_or_else(|| {
                eyre!(
                    "no valid global orientation matches reduced orientation {}",
                    GS.orientation_delta(reduced_orientation)
                )
            })
    }

    /// Build only the selector factors for the internal edges of the integrated subgraph.
    ///
    /// This intentionally does **not** use `orientation_thetas()` on the whole representative
    /// orientation. The reduced-graph CFF term already carries selectors for the edges that remain
    /// explicit after contraction, so re-applying them here would duplicate selectors on external
    /// edges.
    fn internal_orientation_selector(
        representative: &EdgeVec<Orientation>,
        internal_edges: &[EdgeIndex],
    ) -> Atom {
        let mut selector = Atom::num(1);

        for edge in internal_edges {
            match representative[*edge] {
                Orientation::Default => selector *= GS.sign_theta(GS.sign(*edge)),
                Orientation::Reversed => selector *= GS.sign_theta(-GS.sign(*edge)),
                Orientation::Undirected => {}
            }
        }

        selector
    }

    /// Localize a reduced-graph orientation term onto compatible full orientations.
    ///
    /// The reduced term already contains the external-edge selector structure through
    /// `reduced_orientation.orientation_thetas()`. We only add selectors for the integrated
    /// subgraph's internal edges.
    fn localized_orientation_term(
        reduced_expression: &Atom,
        reduced_orientation: &EdgeVec<Orientation>,
        valid_global_orientations: &[EdgeVec<Orientation>],
        internal_edges: &[EdgeIndex],
        average_compatible_orientations: bool,
    ) -> Result<Atom> {
        let compatible: Vec<_> = valid_global_orientations
            .iter()
            .filter(|candidate| {
                Self::orientations_match_outside_integrated_subgraph(reduced_orientation, candidate)
            })
            .collect();
        let compatible_count = compatible.len();
        if compatible_count == 0 {
            return Err(eyre!(
                "no valid global orientation matches reduced orientation {}",
                GS.orientation_delta(reduced_orientation)
            ));
        }

        let representatives: Vec<_> = if average_compatible_orientations {
            compatible
        } else {
            vec![Self::select_representative_orientation(
                reduced_orientation,
                valid_global_orientations,
                internal_edges,
            )?]
        };

        let reduced_selector = reduced_orientation.orientation_thetas();
        let mut localized = Atom::Zero;
        let mut internal_selector_sum = Atom::Zero;
        let mut representative_delta_sum = Atom::Zero;
        for representative in representatives {
            let internal_selector =
                Self::internal_orientation_selector(representative, internal_edges);
            representative_delta_sum += GS.orientation_delta(representative);
            internal_selector_sum += internal_selector.clone();
            localized += reduced_expression.clone() * &reduced_selector * internal_selector;
        }
        if average_compatible_orientations {
            localized /= Atom::num(compatible_count as i64);
            internal_selector_sum /= Atom::num(compatible_count as i64);
            representative_delta_sum /= Atom::num(compatible_count as i64);
        }
        let active_orientation_count = valid_global_orientations
            .iter()
            .filter(|orientation| !orientation.select(localized.as_atom_view()).is_zero())
            .count();

        debug_tags!(#generation, #profile, #uv, #integrated, #local, #orientation, #trace;
            stage = "localized_integrated_ct_orientation_term",
            internal_edges = ?internal_edges,
            compatible_count,
            active_orientation_count,
            average_compatible_orientations,
            log.reduced_delta = GS.orientation_delta(reduced_orientation),
            log.representative_delta = representative_delta_sum,
            log.reduced_selector = reduced_selector,
            log.internal_selector = internal_selector_sum,
            log.localized = localized,
            "Localized integrated UV CT orientation term"
        );

        Ok(localized)
    }

    fn internal_paired_edges_of_subgraph(graph: &Graph, subgraph: &SuBitGraph) -> Vec<EdgeIndex> {
        let mut edge_ids: Vec<_> = graph
            .as_ref()
            .iter_edges_of(subgraph)
            .filter_map(|(pair, edge_id, _)| pair.is_paired().then_some(edge_id))
            .collect();
        edge_ids.sort_by_key(|edge| usize::from(*edge));
        edge_ids
    }

    fn finite_integrated_ct(&self, integrated_4d: &ApproxOp) -> Atom {
        let mut t4 = match integrated_4d {
            ApproxOp::Root => BTreeMap::from([(CutCFFIndex::new_all_none(), Atom::Zero)]),
            _ => integrated_4d
                .expr()
                .unwrap_or(BTreeMap::from([(CutCFFIndex::new_all_none(), Atom::Zero)])),
        };
        if t4.len() != 1 {
            panic!("Should only have one t_arg for the 4d approximation");
        }

        let t4 = t4.remove(&CutCFFIndex::new_all_none()).unwrap();
        let finite = t4
            .series(GS.dim_epsilon, Atom::Zero, 0)
            .unwrap()
            .coefficient((0, 1).into());
        debug_tags!(#generation, #profile, #uv, #integrated, #local, #trace;
            stage = "finite_integrated_ct",
            log.t4 = t4,
            log.finite = finite,
            "Extracted finite integrated UV CT"
        );

        finite
            .map_mink_dim(4)
            .replace(function!(symbol!("vakint::g"), W_.a__))
            .with(function!(symbol!("spenso::g"), W_.a__))
    }

    fn localize_integrated_ct<S: ForestNodeLike>(
        &self,
        graph: &mut Graph,
        integrated_node: &S,
        finite_ct: &Atom,
        cuts: &CutSet,
    ) -> Result<LocalizedIntegratedCt> {
        debug_tags!(#uv; log.expr = finite_ct, "Localizing integrated UV CT");
        if finite_ct.is_zero() {
            let indices = cuts.residue_selector.generate_allowed_keys();
            return Ok(LocalizedIntegratedCt {
                active: IntegrandExpr {
                    integrands: indices.iter().map(|index| (*index, Atom::Zero)).collect(),
                },
                frozen_integrands: indices
                    .into_iter()
                    .map(|index| (index, Atom::one()))
                    .collect(),
            });
        }

        let to_contract = integrated_node.subgraph();
        let integrated_loop_count = graph.n_loops(to_contract);
        // Keep finite addbacks for nested multi-loop entries. Integrated expressions carry
        // their forest-composition signs where needed, so dropping the localized finite
        // representative here removes the Tint(T(...)) terms.
        let finite_ct = finite_ct.clone();
        debug_tags!(#generation, #profile, #uv, #integrated, #local, #summary;
            stage = "localize_integrated_ct_forest_overlap",
            integrated_loop_count,
            "Applied integrated UV forest-overlap addback rule"
        );
        let cff = graph.cff(
            &to_contract
                .union(&graph.tree_edges)
                .subtract(&graph.initial_state_cut),
            cuts,
            self.orientation_pattern,
        )?;

        let fourddenoms = GS.wrap_tree_denoms(
            graph.denominator(&graph.tree_edges.subtract(&graph.initial_state_cut), |_| -1),
        );

        let internal_edges = Self::internal_paired_edges_of_subgraph(graph, to_contract);
        let average_compatible_orientations = true;
        let localizing_integrand = GS.localizing_integrand(integrated_node.lmb());
        debug_tags!(#generation, #profile, #uv, #integrated, #local, #summary;
            stage = "localize_integrated_ct_factors",
            finite_ct_byte_size = finite_ct.as_view().get_byte_size(),
            fourddenoms_byte_size = fourddenoms.as_view().get_byte_size(),
            localizing_integrand_byte_size = localizing_integrand.as_view().get_byte_size(),
            internal_edges = ?internal_edges,
            average_compatible_orientations,
            "Integrated UV CT localization size checkpoint"
        );

        let mut active_integrands = BTreeMap::new();
        let mut frozen_integrands = BTreeMap::new();
        for (index, term) in cff.terms {
            let mut localized = Atom::Zero;
            for (cff_term, orientation) in term.expression.into_iter().zip(term.orientations) {
                localized += Self::localized_orientation_term(
                    &(cff_term * &fourddenoms),
                    &orientation,
                    self.valid_orientations,
                    &internal_edges,
                    average_compatible_orientations,
                )?;
            }

            let localized_cff_byte_size = localized.as_view().get_byte_size();
            let active_ct = localized * finite_ct.clone();
            let localized_ct = &active_ct * &localizing_integrand;
            debug_tags!(#generation, #profile, #uv, #integrated, #local, #term, #summary;
                stage = "localize_integrated_ct_term",
                cut_index = ?index,
                localized_cff_byte_size,
                active_ct_byte_size = active_ct.as_view().get_byte_size(),
                localized_ct_byte_size = localized_ct.as_view().get_byte_size(),
                "Integrated UV CT localization size checkpoint"
            );

            active_integrands.insert(index, active_ct);
            frozen_integrands.insert(index, localizing_integrand.clone());
        }

        Ok(LocalizedIntegratedCt {
            active: IntegrandExpr {
                integrands: active_integrands,
            },
            frozen_integrands,
        })
    }

    #[debug_instrument(
        graph = %graph.log_display(),
        current = %integrated_node.log_display(),
    )]
    pub(crate) fn localized_integrated_ct<S: ForestNodeLike>(
        &self,
        graph: &mut Graph,
        integrated_node: &S,
        integrated_4d: &ApproxOp,
        cuts: &CutSet,
    ) -> Result<LocalizedIntegratedCt> {
        let finite = self.finite_integrated_ct(integrated_4d);

        // let n_loops = graph.n_loops(integrated_node.subgraph());

        // let series = finite
        //     .series(GS.dim_epsilon, Atom::Zero, n_loops as i64 + 1)
        //     .unwrap();
        // let series_atom = series.to_atom();

        // debug_tags!(#uv, #integrated, #inspect, #series;
        //     log.series = series_atom,
        //     "dim epsilon Series "
        // );

        // let mut pole_stripped = Atom::Zero;

        // for (power, p) in series.terms() {
        //     if power >= 0 {
        //         pole_stripped += p * Atom::var(GS.dim_epsilon).pow(power);
        //     }
        // }

        debug_tags!(
            #uv;
            log.finite = finite,
            "Computing localized integrated UV CT",
        );

        self.localize_integrated_ct(graph, integrated_node, &finite, cuts)
    }

    #[debug_instrument(
        graph = %graph.log_display(),
        current = %current.log_display(),
        term_count = local_terms.len(),
    )]
    pub(crate) fn build<S: ForestNodeLike>(
        &self,
        graph: &mut Graph,
        current: &S,
        local_terms: &BTreeMap<CutCFFIndex, Atom>,
        local_sign: Sign,
        integrated_4d: &ApproxOp,
        cutset: &CutSet,
    ) -> Result<BTreeMap<CutCFFIndex, Atom>> {
        let global_num = graph.global_atom();
        debug_tags!(#generation, #profile, #uv, #graph, #summary;
            global_num = %global_num.log_display(),
            "Computed global numerator"
        );

        let integrated_t: IntegrandExpr = self
            .localized_integrated_ct(graph, current, integrated_4d, cutset)?
            .try_into()?;

        let reduced = graph
            .full_filter()
            .subtract(current.subgraph())
            .subtract(&graph.initial_state_cut);

        let mut integrands = BTreeMap::new();

        for (term_index, ((local_index, local), (integrated_index, integ))) in local_terms
            .iter()
            .zip(integrated_t.integrands.iter())
            .enumerate()
        {
            let mut term_started = std::time::Instant::now();
            let mut cff;
            if let Some(marker) = self.uv_marker {
                let local_term = if self.tag_local_terms {
                    local * function!(GS.uv_local, marker.clone())
                } else {
                    local.clone()
                };
                let integrated_term = integ * function!(GS.uv_integrated, marker.clone());
                cff = local_sign * (&local_term - &integrated_term);
                debug_tags!(#generation, #profile, #uv, #graph, #term, #trace;
                    stage = "final_integrand_cff_inputs",
                    cut_index = ?local_index,
                    local_sign = ?local_sign,
                    tag_local_terms = self.tag_local_terms,
                    log.local = local_term,
                    log.integrated = integrated_term,
                    log.cff_before_replacements = cff,
                    "Built final UV CFF inputs"
                );
            } else {
                cff = local_sign * (local - integ);
                debug_tags!(#generation, #profile, #uv, #graph, #term, #trace;
                    stage = "final_integrand_cff_inputs",
                    cut_index = ?local_index,
                    local_sign = ?local_sign,
                    tag_local_terms = self.tag_local_terms,
                    log.local = local,
                    log.integrated = integ,
                    log.cff_before_replacements = cff,
                    "Built final UV CFF inputs"
                );
            }

            if local_index != integrated_index {
                return Err(eyre!(
                    "Mismatched indices for local and integrated approximations: {:?} vs {:?}",
                    local_index,
                    integrated_index
                ));
            }

            for (p, eid, _) in graph.as_ref().iter_edges_of(&reduced) {
                let eid = usize::from(eid) as i64;
                if p.is_paired() {
                    cff = cff
                        .replace(function!(GS.energy, eid))
                        .with(function!(GS.ose, eid));
                }
            }

            cff = cff
                .replace(function!(GS.ose, W_.mass_, W_.prop_))
                .with(W_.prop_);
            debug_tags!(#generation, #profile, #uv, #graph, #term, #summary;
                stage = "final_integrand_cff_done",
                cut_index = ?local_index,
                log.expression = cff,
                elapsed_ms = term_started.elapsed().as_secs_f64() * 1000.0,
                "Built final UV CFF"
            );
            term_started += term_started.elapsed();

            let mut resnum = graph
                .numerator(&reduced, current.subgraph())
                .get_single_atom()
                .unwrap();
            debug_tags!(#generation, #profile, #uv, #graph, #term, #summary;
                stage = "final_integrand_reduced_numerator_done",
                log.expression = resnum,
                cut_index = ?local_index,
                elapsed_ms = term_started.elapsed().as_secs_f64() * 1000.0,
                "Built reduced numerator for final UV integrand"
            );
            term_started += term_started.elapsed();

            let bridgeless_reduced = reduced.subtract(&graph.tree_edges);

            let mut reps = Vec::new();
            for (p, eid, _) in graph.as_ref().iter_edges_of(&bridgeless_reduced) {
                if p.is_paired() {
                    reps.push(GS.add_parametric_sign(eid));
                }
            }

            resnum = resnum.replace_multiple(&reps);
            resnum *= cff * &global_num;

            let color_simplify_input = resnum.replace(GS.dim).with(4);

            debug_tags!(#generation,#algebra,#color,#before, #uv, #inspect, #dump;
                log.expression = color_simplify_input,
                term_index = term_index,
                "Final UV term before color simplification"
            );

            term_started += term_started.elapsed();
            resnum = color_simplify_input
                .collect_factors()
                .simplify_metrics()
                .simplify_color();
            debug_tags!(#generation, #algebra,#color, #uv, #graph, #term, #summary;
                stage = "final_integrand_color_simplify_done",
                cut_index = ?local_index,
                log.expression = resnum,
                elapsed_ms = term_started.elapsed().as_secs_f64() * 1000.0,
                "Simplified final UV term color"
            );

            term_started += term_started.elapsed();
            resnum = resnum.expand_dots()?;
            debug_tags!(#generation, #profile, #uv, #graph, #term, #summary;
                stage = "final_integrand_expand_dots_done",
                cut_index = ?local_index,
                log.expression = resnum,
                elapsed_ms = term_started.elapsed().as_secs_f64() * 1000.0,
                "Expanded final UV term dots"
            );

            resnum = resnum.replace_multiple(&reps);

            debug_tags!(#generation, #profile, #uv, #graph, #term, #summary;
                stage = "final_integrand_replace_multiple_done",
                log.expression = resnum,
                cut_index = ?local_index,
                "Applied final UV parametric replacements"
            );
            integrands.insert(
                *local_index,
                resnum.replace(GS.m_uv_expansion).with(GS.m_uv_vacuum),
            );
        }

        Ok(integrands)
    }
}

#[cfg(test)]
mod tests {
    use super::FinalIntegrand;
    use crate::utils::GS;
    use linnet::half_edge::involution::{EdgeIndex, EdgeVec, Orientation};
    use symbolica::function;

    fn orientation(value: i8) -> Orientation {
        match value {
            1 => Orientation::Default,
            -1 => Orientation::Reversed,
            0 => Orientation::Undirected,
            _ => panic!("invalid orientation encoding"),
        }
    }

    fn edgevec(values: impl IntoIterator<Item = i8>) -> EdgeVec<Orientation> {
        EdgeVec::from_iter(values.into_iter().map(orientation))
    }

    fn edges(values: impl IntoIterator<Item = usize>) -> Vec<EdgeIndex> {
        values.into_iter().map(EdgeIndex).collect()
    }

    #[test]
    fn picks_first_valid_internal_representatives_per_external_class() {
        let valid = vec![
            edgevec([1, -1, 1, -1, 1]),
            edgevec([1, 1, -1, -1, 1]),
            edgevec([1, 1, -1, -1, -1]),
            edgevec([-1, 1, 1, 1, -1]),
        ];
        let internal = edges([2, 3, 4]);

        let reduced_pm = edgevec([1, -1, 0, 0, 0]);
        let reduced_pp = edgevec([1, 1, 0, 0, 0]);
        let reduced_mp = edgevec([-1, 1, 0, 0, 0]);

        assert_eq!(
            FinalIntegrand::select_representative_orientation(&reduced_pm, &valid, &internal)
                .unwrap(),
            &valid[0]
        );
        assert_eq!(
            FinalIntegrand::select_representative_orientation(&reduced_pp, &valid, &internal)
                .unwrap(),
            &valid[1]
        );
        assert_eq!(
            FinalIntegrand::select_representative_orientation(&reduced_mp, &valid, &internal)
                .unwrap(),
            &valid[3]
        );
    }

    #[test]
    fn prefers_more_undirected_internal_edges_first() {
        let valid = vec![edgevec([1, 0, 1]), edgevec([1, 1, 1]), edgevec([1, -1, 1])];
        let reduced = edgevec([1, 0, 1]);
        let internal = edges([1]);

        assert_eq!(
            FinalIntegrand::select_representative_orientation(&reduced, &valid, &internal).unwrap(),
            &valid[0]
        );
    }

    #[test]
    fn prefers_more_default_internal_edges_when_undirected_counts_tie() {
        let valid = vec![
            edgevec([1, 1, -1]),
            edgevec([1, -1, -1]),
            edgevec([1, -1, 1]),
        ];
        let reduced = edgevec([1, 0, 0]);
        let internal = edges([1, 2]);

        assert_eq!(
            FinalIntegrand::select_representative_orientation(&reduced, &valid, &internal).unwrap(),
            &valid[0]
        );
    }

    #[test]
    fn breaks_default_ties_by_lexicographically_maximal_leading_defaults() {
        let valid = vec![
            edgevec([1, 1, -1, 1]),
            edgevec([1, -1, 1, 1]),
            edgevec([1, 1, 1, -1]),
        ];
        let reduced = edgevec([1, 0, 0, 0]);
        let internal = edges([1, 2, 3]);

        assert_eq!(
            FinalIntegrand::select_representative_orientation(&reduced, &valid, &internal).unwrap(),
            &valid[2]
        );
    }

    #[test]
    fn errors_for_missing_external_orientation_class() {
        let valid = vec![edgevec([1, 1, 1]), edgevec([1, -1, 1])];
        let reduced = edgevec([-1, 0, 1]);
        let internal = edges([1]);

        assert!(
            FinalIntegrand::select_representative_orientation(&reduced, &valid, &internal).is_err()
        );
    }

    #[test]
    fn selector_only_depends_on_internal_edges() {
        let representative = edgevec([1, -1, 1, -1]);
        let selector =
            FinalIntegrand::internal_orientation_selector(&representative, &edges([1, 3]));

        let expected =
            GS.sign_theta(-GS.sign(EdgeIndex(1))) * GS.sign_theta(-GS.sign(EdgeIndex(3)));
        assert_eq!(selector, expected);
    }

    #[test]
    fn orientation_term_keeps_external_selectors_and_adds_internal_ones() {
        let reduced_expression = function!(GS.ose, 0);
        let reduced_orientation = edgevec([1, 0, -1]);
        let valid = vec![edgevec([1, 1, -1]), edgevec([1, -1, -1])];
        let localized = FinalIntegrand::localized_orientation_term(
            &reduced_expression,
            &reduced_orientation,
            &valid,
            &edges([1]),
            false,
        )
        .unwrap();

        let expected = reduced_expression
            * GS.sign_theta(GS.sign(EdgeIndex(0)))
            * GS.sign_theta(-GS.sign(EdgeIndex(2)))
            * GS.sign_theta(GS.sign(EdgeIndex(1)));
        assert_eq!(localized, expected);
    }
}
