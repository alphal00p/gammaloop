#[cfg(test)]
use crate::cff::CutCFFIndex;
use crate::{
    debug_tags,
    graph::Graph,
    utils::{GS, W_},
    uv::{
        Integrands, UVgenerationSettings, UltravioletGraph,
        approx::{
            ForestNodeLike,
            direct_3d::{Direct3dCts, DirectResidueBranches},
            integrated::IntegratedCts,
            local_3d::{Localizer, OrientationIntegrands},
            projected_4d::Projected4dCts,
        },
        marker::UvMarker,
    },
};
use color_eyre::Result;
use gammaloop_tracing_filter::{LogMessage, debug_instrument};
use idenso::{
    color::{ColorSimplifier, ColorSimplifySettings},
    shorthands::metric::MetricSimplifier,
};
use linnet::half_edge::subgraph::{SuBitGraph, SubSetLike, SubSetOps};
use symbolica::{
    atom::{Atom, AtomCore},
    function,
    id::Replacement,
};
use three_dimensional_reps::CffGenerationContext;
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct FinalIntegrands(Integrands);

impl FinalIntegrands {
    /// Iterate over finalized semantic expressions for diagnostics. Deferred
    /// evaluator calls are materialized before leaving this stage.
    #[cfg(test)]
    pub(crate) fn iter(&self) -> impl Iterator<Item = (CutCFFIndex, Atom)> {
        self.0.materialize().compact.into_iter()
    }

    pub(crate) fn map(&self, f: impl FnMut(&Atom) -> Atom) -> Self {
        Self(self.0.map(f))
    }

    pub(crate) fn zip_add(self, other: Self) -> Result<Self> {
        Ok(Self(self.0.zip_add(other.0)?))
    }

    pub(crate) fn materialize(&self) -> Integrands {
        self.0.materialize()
    }

    pub(crate) fn into_integrands(self) -> Integrands {
        self.0
    }
}

pub(crate) struct FinalIntegrandBuilder<'a> {
    localizer: Localizer<'a>,
    marker: UvMarker,
}

impl<'a> FinalIntegrandBuilder<'a> {
    pub(crate) fn new(localizer: Localizer<'a>, settings: &UVgenerationSettings) -> Self {
        Self {
            localizer,
            marker: UvMarker::new(settings),
        }
    }

    #[debug_instrument(
        graph = %graph.log_display(),
        current = %current.log_display(),
    )]
    pub(crate) fn build_direct<S: ForestNodeLike>(
        &self,
        graph: &mut Graph,
        current: &S,
        local_terms: &Direct3dCts,
        integrated: &IntegratedCts,
    ) -> Result<FinalIntegrands> {
        let reduced = graph
            .full_filter()
            .subtract(current.subgraph())
            .subtract(&graph.initial_state_cut);
        let full_graph = graph.full_filter();

        let global_num = graph.global_atom();
        debug_tags!(#generation, #profile, #uv, #graph, #summary;
            global_num = %global_num.log_display(),
            "Computed global numerator"
        );

        let resnum = graph
            .numerator(&reduced, current.subgraph())
            .get_single_atom()
            .expect("graph numerator should be available")
            * global_num;
        let localized_integrated = self
            .localizer
            .localize(
                &integrated.physical_finite_counterterm_atom(),
                graph,
                current,
            )?
            .combine()?;
        let localized_integrated = DirectResidueBranches::from_transient(&localized_integrated)?
            .map(|atom| self.marker.prefix(&full_graph, current.subgraph(), atom));
        let localized_local = local_terms
            .branches()?
            .map(|atom| self.marker.prefix(&full_graph, current.subgraph(), atom));
        let final_branches = localized_integrated
            .zip_add(&localized_local)?
            .multiply_key_mapped(self.localizer.orientation, graph, &resnum)?;

        // `DirectResidueBranches` is the sparse, factorization-preserving
        // representation of sum_k sigma(k) I_k while the Taylor forest is
        // built. Materialize that opaque residue-map-key selector only at the
        // evaluator boundary, after every branch-owned numerator map has been
        // applied. An explicit sum simply replaces every sigma by one; physical
        // edge directions are separate sign metadata.
        let selected = final_branches.materialize(
            self.localizer.uses_exact_maps()
                && !self.localizer.orientation.explicit_orientation_sum_only,
        )?;
        let coarse_energy_replacements = if self.localizer.uses_exact_maps() {
            None
        } else {
            let bridgeless_reduced = reduced.subtract(&graph.tree_edges);
            let mut replacements = Vec::new();
            for (pair, edge_id, _) in graph.as_ref().iter_edges_of(&bridgeless_reduced) {
                if pair.is_paired() {
                    replacements.push(GS.add_parametric_sign(edge_id));
                }
            }
            Some(replacements)
        };
        Self::simplify_final(
            graph,
            &reduced,
            selected,
            coarse_energy_replacements.as_deref(),
        )
    }

    #[debug_instrument(
        graph = %graph.log_display(),
        current = %current.log_display(),
    )]
    pub(crate) fn build_projected<S: ForestNodeLike>(
        &self,
        graph: &mut Graph,
        current: &S,
        local_terms: &Projected4dCts,
        integrated: &IntegratedCts,
    ) -> Result<FinalIntegrands> {
        if current.subgraph().is_empty() {
            return Err(eyre::eyre!(
                "the empty forest must own the complete production CFF, not projected local-4D coefficients"
            ));
        }

        let reduced = graph
            .full_filter()
            .subtract(current.subgraph())
            .subtract(&graph.initial_state_cut);
        let full_graph = graph.full_filter();

        let global_num = graph.global_atom();
        debug_tags!(#generation, #profile, #uv, #graph, #summary;
            global_num = %global_num.log_display(),
            "Computed global numerator"
        );
        let resnum = graph
            .numerator(&reduced, current.subgraph())
            .get_single_atom()
            .expect("graph numerator should be available")
            * global_num;

        // Only the projected local-4D route reaches this assembly boundary.
        // Its child Taylor coefficient deliberately omits the untouched
        // cograph; construct that outer CFF exactly once here.
        let active_sectors = local_terms.sectors();
        // Restore the production-tree denominators after the outer CFF
        // contracts its loop-energy dependence.  The exact DDx GL0 UV ray
        // verifies that this is the full production tree, including the
        // carrier shared with the factorized self-energy coefficient.
        let fourddenoms = GS.wrap_tree_denoms(
            graph.denominator(&graph.tree_edges.subtract(&graph.initial_state_cut), |_| -1),
        );
        let mut localized_local: Option<OrientationIntegrands> = None;
        for (_, sector) in active_sectors {
            let localized = if sector
                .active
                .iter_orientations()
                .all(|(_, _, integrands)| integrands.iter().all(|(_, atom)| atom.is_zero()))
            {
                // A disabled integrated prefix can deliberately retain a
                // typed zero sector for later forest replay. Preserve all
                // selector and residue keys, but do not ask the outer CFF
                // to resolve an orientation for a coefficient that is
                // identically zero in every hosted branch.
                sector.combine()?
            } else {
                // The child Taylor coefficient remains factorized while its
                // crown-energy dependence contributes to the outer CFF bound.
                // Cut orders and selector/source-map branches are evaluated
                // independently. Keep them independent in the capacity
                // oracle too, so opposite leading powers cannot cancel
                // before generalized-CFF generation.
                let analysis_numerator = sector.active.factorized_capacity_envelope() * &resnum;
                let outer = self
                    .localizer
                    .projected_cff(
                        graph,
                        current.subgraph(),
                        &analysis_numerator,
                        CffGenerationContext::FactorizedContour,
                    )?
                    .map(|atom| atom * &fourddenoms);
                let active = outer.zip_mul_mapped_factor(
                    &sector.active,
                    |orientation_id, source_edge_energy_map, active| {
                        let factorized_numerator = active * &resnum;
                        self.localizer.map_numerator(
                            graph,
                            orientation_id,
                            source_edge_energy_map,
                            &factorized_numerator,
                        )
                    },
                )?;
                active.zip_mul_unmapped(&sector.frozen_integrands)?
            };
            localized_local = Some(match localized_local {
                Some(sum) => sum.zip_add(&localized)?,
                None => localized,
            });
        }
        let localized_local = localized_local
            .ok_or_else(|| eyre::eyre!("factorized local term has no active UV sectors"))?
            .map(|atom| self.marker.prefix(&full_graph, current.subgraph(), atom));
        let localized_integrated = self
            .localizer
            .localize(
                &integrated.physical_finite_counterterm_atom(),
                graph,
                current,
            )?
            .combine()?
            .multiply_mapped(|orientation_id, source_edge_energy_map| {
                self.localizer
                    .map_numerator(graph, orientation_id, source_edge_energy_map, &resnum)
            })?
            .map(|atom| self.marker.prefix(&full_graph, current.subgraph(), atom));
        let final_branches = localized_integrated.zip_add(&localized_local)?;

        // Both legitimate projected maps have now consumed every still-unmapped
        // numerator factor. Production hosts are mapping metadata only in this
        // lane: sum them explicitly without ever materializing a selector or
        // traversing another numerator map.
        let mut hosted = final_branches.iter_orientations();
        let (_, _, first) = hosted.next().ok_or_else(|| {
            eyre::eyre!("final 3D UV integrand contains no production energy maps")
        })?;
        let mut selector_free = first.clone();
        for (_, _, integrands) in hosted {
            selector_free = selector_free.zip_add(integrands.clone())?;
        }
        Self::simplify_final(graph, &reduced, selector_free, None)
    }

    /// Normalize an already mapped and selector-assembled final integrand. This
    /// tail is deliberately blind to residue maps and cannot map a numerator.
    fn simplify_final(
        graph: &Graph,
        reduced: &SuBitGraph,
        integrands: Integrands,
        coarse_energy_replacements: Option<&[Replacement]>,
    ) -> Result<FinalIntegrands> {
        let simplified = integrands.fallible_map(|atom| {
            let mut atom = atom.clone();
            for (pair, edge_id, _) in graph.as_ref().iter_edges_of(reduced) {
                if pair.is_paired() {
                    let edge_id = usize::from(edge_id) as i64;
                    atom = atom
                        .replace(function!(GS.energy, edge_id))
                        .with(function!(GS.ose, edge_id));
                }
            }
            atom = atom
                .replace(function!(GS.ose, W_.mass_, W_.prop_))
                .with(W_.prop_);
            atom = atom
                .replace(GS.dim)
                .with(4)
                .collect_factors()
                .simplify_metrics()
                .simplify_color_with(
                    ColorSimplifySettings::default().with_cof_dimension_invariants(),
                )
                .expand_dots()?;

            // Coarse export projectors still introduce parametric signs at
            // this legacy boundary. Exact production branches have already
            // mapped every owned numerator fragment and must not fall back
            // to a second coarse replacement here.
            if let Some(replacements) = coarse_energy_replacements {
                atom = atom.replace_multiple(replacements);
            }
            Ok(atom
                .replace(GS.m_uv_expansion)
                .with(GS.m_uv_vacuum)
                .replace(GS.dim_epsilon)
                .with(0))
        })?;
        Ok(FinalIntegrands(simplified))
    }
}
