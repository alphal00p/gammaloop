use color_eyre::Result;
use linnet::half_edge::subgraph::{SubSetLike, SubSetOps};
use symbolica::atom::Atom;
use three_dimensional_reps::CffGenerationContext;

use crate::{
    debug_tags,
    graph::Graph,
    utils::GS,
    uv::{UltravioletGraph, approx::ForestNodeLike},
};

use super::{FrozenActiveCt, Localizer, OrientationIntegrands};

impl Localizer<'_> {
    pub(crate) fn localize<S: ForestNodeLike>(
        &self,
        expr: &Atom,
        graph: &mut Graph,
        integrated_node: &S,
    ) -> Result<FrozenActiveCt> {
        if expr.is_zero() {
            let indices = self.cutset.residue_selector.generate_allowed_keys();
            return Ok(FrozenActiveCt {
                active: OrientationIntegrands::from_ids_and_indices(
                    self.orientation.orientation_ids(),
                    &indices,
                ),
                frozen_integrands: indices
                    .into_iter()
                    .map(|index| (index, Atom::one()))
                    .collect(),
                active_lmb: None,
            });
        }

        let to_contract = integrated_node.subgraph();
        let integrated_loop_count = graph.n_loops(to_contract);
        // Keep finite addbacks for nested multi-loop entries. Integrated expressions carry
        // their forest-composition signs where needed, so dropping the localized finite
        // representative here removes the Tint(T(...)) terms.
        let finite_ct = expr.clone();
        debug_tags!(#generation, #profile, #uv, #integrated, #local, #summary;
            stage = "localize_integrated_ct_forest_overlap",
            integrated_loop_count,
            "Applied integrated UV forest-overlap addback rule"
        );
        let fourddenoms = GS.wrap_tree_denoms(
            graph.denominator(&graph.tree_edges.subtract(&graph.initial_state_cut), |_| -1),
        );

        // CFF capacity belongs to the complete factorized numerator of this
        // reduced term. The integrated finite CT replaces the numerator of
        // the contracted spinney; the remaining graph and global factors are
        // still grown later by their existing owners.
        let reduced = graph
            .full_filter()
            .subtract(to_contract)
            .subtract(&graph.initial_state_cut);
        let outside_numerator = graph
            .numerator(&reduced, to_contract)
            .get_single_atom()
            .expect("Graph numerator should be available")
            * graph.global_atom();
        let analysis_numerator = expr * &outside_numerator;

        let localizing_integrand = GS.localizing_integrand(integrated_node.lmb());
        let active = self
            .projected_cff(
                graph,
                to_contract,
                &analysis_numerator,
                CffGenerationContext::Standalone,
            )?
            .fallible_map(|orientation_id, source_edge_energy_map, localized| {
                let localized = localized * &fourddenoms;
                let localized_cff_byte_size = localized.as_view().get_byte_size();
                let mapped_finite_ct = self.orientation.map_numerator(
                    graph,
                    orientation_id,
                    source_edge_energy_map,
                    &finite_ct,
                )?;
                let active_ct = localized * mapped_finite_ct;
                let localized_ct = &active_ct * &localizing_integrand;
                debug_tags!(#generation, #profile, #uv, #integrated, #local, #term, #summary;
                    stage = "localize_integrated_ct_term",
                    integrated_node = %integrated_node.log_display(),
                    contracted = %to_contract.string_label(),
                    reduced = %reduced.string_label(),
                    residue_map_key = orientation_id.0,
                    localized_cff_byte_size,
                    active_ct_byte_size = active_ct.as_view().get_byte_size(),
                    localized_ct_byte_size = localized_ct.as_view().get_byte_size(),
                    "Integrated UV CT localization size checkpoint"
                );
                Ok(active_ct)
            })?;
        let frozen_integrands = self
            .cutset
            .residue_selector
            .generate_allowed_keys()
            .into_iter()
            .map(|index| (index, localizing_integrand.clone()))
            .collect();

        Ok(FrozenActiveCt {
            active,
            frozen_integrands,
            active_lmb: None,
        })
    }
}
