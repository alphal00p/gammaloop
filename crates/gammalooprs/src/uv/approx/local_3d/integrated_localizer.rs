use color_eyre::Result;
use idenso::tensor::SymbolicTensor;
use linnet::half_edge::subgraph::{SubSetLike, SubSetOps};
use spenso::structure::partial::PartialStructure;
use symbolica::atom::Atom;
use three_dimensional_reps::CffGenerationContext;

use crate::{
    debug_tags,
    graph::Graph,
    utils::GS,
    uv::{
        UltravioletGraph,
        approx::{ForestNodeLike, direct_3d::DirectResidueBranches},
    },
};

use super::{FrozenActiveCt, Localizer, OrientationIntegrandBranch, OrientationIntegrands};

impl Localizer<'_> {
    pub(crate) fn localize<S: ForestNodeLike>(
        &self,
        expr: &SymbolicTensor<PartialStructure>,
        graph: &mut Graph,
        integrated_node: &S,
    ) -> Result<FrozenActiveCt> {
        if expr.is_zero() {
            let indices = self.cutset.residue_selector.generate_allowed_keys();
            return Ok(FrozenActiveCt {
                active: OrientationIntegrands::from_ids_and_indices(
                    self.orientation.orientation_ids()?,
                    &indices,
                ),
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
        let finite_ct = expr;
        debug_tags!(#generation, #profile, #uv, #integrated, #local, #summary;
            stage = "localize_integrated_ct_forest_overlap",
            integrated_loop_count,
            "Applied integrated UV forest-overlap addback rule"
        );
        let fourddenoms = GS.wrap_tree_denoms(graph.denominator(
            &graph.tree_edges.subtract(&graph.initial_state_cut),
            &graph.model,
            |_| -1,
        ));

        // CFF capacity belongs to the complete factorized numerator of this
        // reduced term. The integrated finite CT replaces the numerator of
        // the contracted spinney; the remaining graph and global factors are
        // still grown later by their existing owners.
        let reduced = graph
            .full_filter()
            .subtract(to_contract)
            .subtract(&graph.initial_state_cut);
        let outside_numerator = graph.numerator(&reduced, to_contract)?.state.expr;
        let global = graph.global_numerator()?;
        let analysis_numerator = expr.multiply(&outside_numerator)?.multiply(&global)?;

        let localizing_integrand = GS.localizing_integrand(integrated_node.lmb());
        // Integration replaces the contracted spinney by an independent
        // coefficient. Its cograph residue maps need not extend to a complete
        // production map, but each still belongs to the finite addback once.
        let mut active = self.with_independent_source_sum().projected_cff(
            graph,
            to_contract,
            [analysis_numerator.expression()],
            CffGenerationContext::Standalone,
        )?;
        active = active.map(|localized| localized * &fourddenoms);
        // The finite effective vertex has the same tensor ownership as every
        // other mapped numerator. Keep its shared body separate from scalar
        // CFF weights, so addition to a local numerator family cannot mix an
        // explicit matrix with an unresolved scalar family call. The existing
        // residue-key owner maps it exactly once using the source chart.
        let mapped = DirectResidueBranches::from_transient(&active)?.multiply_key_mapped(
            self.orientation,
            graph,
            finite_ct,
            DirectResidueBranches::numerator_scope(),
        )?;
        let active = OrientationIntegrands(
            mapped
                .iter_keys()
                .map(|(key, integrands)| {
                    debug_tags!(#generation, #profile, #uv, #integrated, #local, #term, #summary;
                        stage = "localize_integrated_ct_term",
                        integrated_node = %integrated_node.log_display(),
                        contracted = %to_contract.string_label(),
                        reduced = %reduced.string_label(),
                        residue_map_key = key.selector_host.0,
                        scalar_weight_bytes = integrands.iter().map(|(_, atom)| atom.as_view().get_byte_size()).sum::<usize>(),
                        shared_numerator_bytes = integrands.numerators().iter().map(|entry| entry.rhs.expression().as_view().get_byte_size()).sum::<usize>(),
                        frozen_kernel_bytes = localizing_integrand.as_view().get_byte_size(),
                        "Factorized integrated UV CT localization size checkpoint"
                    );
                    OrientationIntegrandBranch {
                        selector_id: key.selector_host,
                        source_edge_energy_map: key.source_edge_energy_map().map(<[_]>::to_vec),
                        integrands: integrands.clone(),
                    }
                })
                .collect(),
        );
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
        })
    }
}
