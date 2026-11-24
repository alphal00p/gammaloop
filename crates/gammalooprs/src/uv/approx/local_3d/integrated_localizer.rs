use std::collections::BTreeSet;

use color_eyre::Result;
use itertools::Itertools;
use linnet::half_edge::involution::EdgeIndex;
use linnet::half_edge::subgraph::{SubSetLike, SubSetOps};
use symbolica::atom::Atom;
use thiserror::Error;
use three_dimensional_reps::{CffGenerationContext, GenerationError, UnlocalizedEnergyContacts};

use crate::{
    debug_tags,
    graph::Graph,
    utils::GS,
    uv::{UltravioletGraph, approx::ForestNodeLike},
};

use super::{FrozenActiveCt, Localizer, OrientationIntegrands};

/// A failure to localize an integrated counterterm on its cograph.
#[derive(Debug, Error)]
pub(crate) enum IntegratedLocalizationError {
    /// The counterterm's numerator raises the cograph's loop-energy degree
    /// above its propagators', so contracting the last of them leaves a
    /// polynomial energy contact and no denominator. It has no residue and a
    /// strict 3D integrand cannot represent it; vacuum subtraction removes it
    /// exactly. A contact beside surviving denominators, which vacuum
    /// subtraction keeps, is not this error.
    #[error(
        "the integrated counterterm of {{{}}} cannot be localized on its cograph {{{}}} \
         without vacuum subtraction: the cograph's numerator exceeds the loop-energy degree \
         of its propagators, which leaves a polynomial energy contact of power {power} that \
         only vacuum subtraction removes; set `global.generation.medium.vacuum_subtraction = \
         true`",
        .node.iter().join(","),
        .cograph.iter().join(",")
    )]
    EnergyContactRequiresVacuumSubtraction {
        node: BTreeSet<EdgeIndex>,
        cograph: BTreeSet<EdgeIndex>,
        power: usize,
    },
}

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
        // The finite counterterm can raise the cograph's loop-energy degree
        // above its propagators', as the three-gluon vertices of a gluon loop
        // between two integrated bubbles do. Contracting the last of them then
        // leaves a quotient polynomial in that energy, which has no residue:
        // beside no denominator its weight is the bare medium unit, which the
        // vacuum subtraction of the localized weight maps to zero, so
        // generation may omit it there. Without vacuum subtraction no 3D
        // integrand represents it; beside surviving denominators vacuum
        // subtraction keeps it, so the generator's error stays.
        let vacuum_subtraction = self.orientation.vacuum_subtraction;
        let contacts = if vacuum_subtraction {
            UnlocalizedEnergyContacts::OmitUnderVacuumSubtraction
        } else {
            UnlocalizedEnergyContacts::Reject
        };
        // Integration replaces the contracted spinney by an independent
        // coefficient. Its cograph residue maps need not extend to a complete
        // production map, but each still belongs to the finite addback once.
        let mut active = self
            .with_independent_source_sum()
            .with_unlocalized_energy_contacts(contacts)
            .projected_cff(
                graph,
                to_contract,
                [&analysis_numerator],
                CffGenerationContext::Standalone,
            )
            .map_err(|error| {
                let removable_contact_power = match error.downcast_ref::<GenerationError>() {
                    Some(GenerationError::UnlocalizedEnergyContact {
                        power, remainder, ..
                    }) if remainder.is_removed_by_vacuum_subtraction() => Some(*power),
                    _ => None,
                };
                match removable_contact_power {
                    Some(power) if !vacuum_subtraction => error.wrap_err(
                        IntegratedLocalizationError::EnergyContactRequiresVacuumSubtraction {
                            node: graph.paired_edges(to_contract).into_iter().collect(),
                            cograph: graph.paired_edges(&reduced).into_iter().collect(),
                            power,
                        },
                    ),
                    _ => error,
                }
            })?;
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
                let localized = localized * &fourddenoms;
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
        if self.orientation.cff_options()?.representation
            == three_dimensional_reps::RepresentationMode::Ltd
        {
            // A localized LTD source can have higher residue orders than the
            // parent. The smooth factor multiplies every actual source order.
            let mut localized = FrozenActiveCt::from(active);
            localized.frozen_integrands = localized
                .frozen_integrands
                .map(|_| localizing_integrand.clone());
            Ok(localized)
        } else {
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
}
