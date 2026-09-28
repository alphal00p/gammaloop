use color_eyre::Result;
use idenso::IndexTooling;
use linnet::half_edge::subgraph::{SubSetLike, SubSetOps};
use symbolica::atom::Atom;
use three_dimensional_reps::CffGenerationContext;

use crate::{
    debug_tags,
    graph::Graph,
    numerator::{aind::Aind, symbolica_ext::NumeratorAtomExt},
    utils::GS,
    uv::{
        UltravioletGraph,
        approx::{ForestNodeLike, direct_3d::DirectResidueBranches},
    },
};

use super::{FrozenActiveCt, Localizer, OrientationIntegrandBranch, OrientationIntegrands};

impl Localizer<'_> {
    /// A supplied final cograph numerator is attached together with the finite
    /// coefficient as one generic local-4D product. Without it, direct CFF3D
    /// replay retains its Taylor-visible, already mapped finite coefficient.
    pub(crate) fn localize<S: ForestNodeLike>(
        &self,
        expr: &Atom,
        graph: &mut Graph,
        integrated_node: &S,
        cograph_numerator: Option<&Atom>,
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
        let fourddenoms = GS.wrap_tree_denoms(
            graph.denominator(&graph.tree_edges.subtract(&graph.initial_state_cut), |_| -1),
        );

        // Capacity belongs to the complete factorized numerator. The finite CT
        // replaces the contracted spinney. Final local-4D callers supply the
        // exact remaining graph/global product; direct CFF3D uses that product
        // only for analysis here and attaches it later during forest replay.
        let reduced = graph
            .full_filter()
            .subtract(to_contract)
            .subtract(&graph.initial_state_cut);
        let analysis_numerator = if let Some(outside_numerator) = cograph_numerator {
            // This is the exact final cograph/global product supplied by its
            // owner. Preserve each factor's private tensor contractions while
            // retaining the free physical interface between the two factors.
            let factors = [expr, outside_numerator]
                .into_iter()
                .map(|body| {
                    let interface = body.list_dangling::<Aind>()?.into_iter().collect();
                    Ok(body.freshen_private_indices(&interface, || {
                        DirectResidueBranches::numerator_scope().1
                    }))
                })
                .collect::<Result<Vec<_>>>()?;
            Atom::mul_many(factors)
        } else {
            let outside_numerator = graph
                .numerator(&reduced, to_contract)
                .get_single_atom()
                .expect("Graph numerator should be available")
                * graph.global_atom();
            expr * outside_numerator
        };

        let localizing_integrand = GS.localizing_integrand(integrated_node.lmb());
        // Integration replaces the contracted spinney by an independent
        // coefficient. Its cograph residue maps need not extend to a complete
        // production map, but each still belongs to the finite addback once.
        let mut active = self.with_independent_source_sum().projected_cff(
            graph,
            to_contract,
            [&analysis_numerator],
            CffGenerationContext::Standalone,
        )?;
        active = active.map(|localized| localized * &fourddenoms);
        if cograph_numerator.is_some() {
            let retained = DirectResidueBranches::from_transient(&active)?.multiply_key_mapped(
                self.orientation,
                graph,
                &analysis_numerator,
                DirectResidueBranches::numerator_scope(),
                true,
            )?;
            active = OrientationIntegrands(
                retained
                    .iter_keys()
                    .map(|(key, integrands)| OrientationIntegrandBranch {
                        selector_id: key.selector_host,
                        source_edge_energy_map: key.source_edge_energy_map().map(<[_]>::to_vec),
                        integrands: integrands.clone(),
                    })
                    .collect(),
            );
        } else {
            for branch in &mut active.0 {
                // A source map owns the finite coefficient for every cut order;
                // attaching it must not repeat the same numerator substitution.
                let mapped_finite_ct = self.orientation.map_numerator(
                    graph,
                    branch.selector_id,
                    branch.source_edge_energy_map.as_deref(),
                    finite_ct,
                )?;
                branch.integrands = branch.integrands.map(|localized| {
                    let localized_cff_byte_size = localized.as_view().get_byte_size();
                    let active_ct = localized * &mapped_finite_ct;
                    debug_tags!(#generation, #profile, #uv, #integrated, #local, #term, #summary;
                        stage = "localize_integrated_ct_term",
                        integrated_node = %integrated_node.log_display(),
                        contracted = %to_contract.string_label(),
                        reduced = %reduced.string_label(),
                        residue_map_key = branch.selector_id.0,
                        localized_cff_byte_size,
                        active_ct_byte_size = active_ct.as_view().get_byte_size(),
                        localized_ct_byte_size = (&active_ct * &localizing_integrand).as_view().get_byte_size(),
                        "Integrated UV CT localization size checkpoint"
                    );
                    active_ct
                });
            }
        }
        // A localized LTD cograph can have higher residue orders than its
        // parent's cut selector. The smooth localization factor multiplies
        // every order of this actual source, with its numerator maps intact.
        let mut localized = FrozenActiveCt::from(active);
        localized.frozen_integrands = localized
            .frozen_integrands
            .map(|_| localizing_integrand.clone());
        Ok(localized)
    }
}
