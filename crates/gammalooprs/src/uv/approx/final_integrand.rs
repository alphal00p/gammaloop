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
            local_3d::Localizer,
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
use linnet::half_edge::subgraph::{Inclusion, SuBitGraph, SubSetLike, SubSetOps};
use symbolica::{
    atom::{Atom, AtomCore},
    function,
    id::Replacement,
};
use three_dimensional_reps::CffGenerationContext;
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct FinalIntegrands(Integrands);

impl FinalIntegrands {
    /// Iterate over finalized semantic expressions for diagnostics.
    #[cfg(test)]
    pub(crate) fn iter(&self) -> impl Iterator<Item = (CutCFFIndex, Atom)> {
        self.0.iter().map(|(index, atom)| (*index, atom.clone()))
    }

    pub(crate) fn map(&self, f: impl FnMut(&Atom) -> Atom) -> Self {
        Self(self.0.map(f))
    }

    pub(crate) fn zip_add(self, other: Self) -> Result<Self> {
        Ok(Self(self.0.zip_add(other.0)?))
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
        let selected = final_branches
            .materialize(!self.localizer.orientation.explicit_orientation_sum_only)?;
        Self::simplify_final(graph, &reduced, selected)
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
        let localizer = self.localizer.with_independent_source_sum();
        // Only the projected local-4D route reaches this assembly boundary.
        // Its child Taylor coefficient deliberately omits the untouched
        // cograph; choose its outer CFF per independent sector here, converting
        // only the selected raw routing proposal into physical surfaces.
        let active_sectors = local_terms.sectors();
        if active_sectors.is_empty() {
            return Err(eyre::eyre!(
                "factorized local term has no active UV sectors"
            ));
        }
        // Restore the production-tree denominators after the outer CFF
        // contracts its loop-energy dependence.  The exact DDx GL0 UV ray
        // verifies that this is the full production tree, including the
        // carrier shared with the factorized self-energy coefficient.
        let fourddenoms = GS.wrap_tree_denoms(
            graph.denominator(&graph.tree_edges.subtract(&graph.initial_state_cut), |_| -1),
        );
        let allowed_zero: Integrands = localizer
            .cutset
            .residue_selector
            .generate_allowed_keys()
            .into_iter()
            .map(|index| (index, Atom::Zero))
            .collect();
        let active_edges = graph
            .iter_edges_of(&reduced)
            .filter_map(|(pair, edge, data)| {
                (pair.is_paired()
                    && !data.data.is_dummy
                    && !graph.tree_edges.includes(&graph[&edge].1))
                .then_some(edge)
            })
            .collect::<Vec<_>>();
        let mut selector_free: Option<Integrands> = None;
        for sector in active_sectors {
            if sector.coefficient.is_zero() {
                // A disabled integrated prefix can deliberately retain a
                // typed zero sector for later forest replay. Preserve all
                // allowed cut orders without asking the outer CFF to resolve
                // a map for an identically zero coefficient.
                selector_free.get_or_insert_with(|| allowed_zero.clone());
                continue;
            }
            // Choose the soft Taylor routing once with the untouched cograph
            // numerator present. The child remains factorized, and both the
            // capacity oracle and every outer map consume this same expression.
            // Independent sectors are never summed before deriving their ranks.
            let (numerator, localized) = localizer.projected_cff_from_soft_momentum_proposals(
                graph,
                current.subgraph(),
                &(&sector.coefficient * &resnum),
                active_edges.iter().copied(),
                CffGenerationContext::EmbeddedCffFactor,
            )?;
            let localized = localized
                .map(|atom| atom * &fourddenoms)
                .multiply_mapped(|orientation_id, source_edge_energy_map| {
                    localizer.map_numerator(
                        graph,
                        orientation_id,
                        source_edge_energy_map,
                        &numerator,
                    )
                })?
                .map(|atom| {
                    self.marker.prefix(
                        &full_graph,
                        current.subgraph(),
                        &(atom * &sector.frozen_factor),
                    )
                });
            // Child contours have already been summed. Only actual outer
            // branches carry maps and cut orders; consume each map once before
            // adding its value, checking the complete allowed cut-key shape.
            for (_, _, integrands) in localized.iter_orientations() {
                let sum = selector_free.take().unwrap_or_else(|| allowed_zero.clone());
                selector_free = Some(sum.zip_add(integrands.clone())?);
            }
        }
        // The integrated addback is shared with the direct route. Its localizer
        // hosts the independent cograph sum in both routes, retaining every
        // source map even when no compatible physical-prefix host survives.
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
        // Both legitimate projected maps have now consumed every still-unmapped
        // numerator factor. Production hosts are mapping metadata only in this
        // lane: sum them explicitly without ever materializing a selector or
        // traversing another numerator map.
        for (_, _, integrands) in localized_integrated.iter_orientations() {
            let sum = selector_free.take().unwrap_or_else(|| allowed_zero.clone());
            selector_free = Some(sum.zip_add(integrands.clone())?);
        }
        let selector_free = selector_free.ok_or_else(|| {
            eyre::eyre!("final 3D UV integrand contains no production energy maps")
        })?;
        Self::simplify_final(graph, &reduced, selector_free)
    }

    /// Normalize an already mapped and selector-assembled final integrand. This
    /// tail is deliberately blind to residue maps and cannot map a numerator.
    fn simplify_final(
        graph: &Graph,
        reduced: &SuBitGraph,
        integrands: Integrands,
    ) -> Result<FinalIntegrands> {
        let energy_replacements = graph
            .as_ref()
            .iter_edges_of(reduced)
            .filter(|(pair, _, _)| pair.is_paired())
            .map(|(_, edge_id, _)| {
                let edge_id = usize::from(edge_id) as i64;
                Replacement::new(function!(GS.energy, edge_id), function!(GS.ose, edge_id))
            })
            .collect::<Vec<_>>();
        let simplified = integrands.fallible_map(|atom| {
            let mut atom = atom
                .replace_multiple(&energy_replacements)
                .replace(function!(GS.ose, W_.mass_, W_.prop_))
                .with(W_.prop_);
            // Preserve the sum of CFF denominators after residue mapping, just
            // as the Taylor stage preserves its separate propagator topologies.
            atom = atom
                .replace(GS.dim)
                .with(4)
                .simplify_metrics()
                .simplify_color_with(
                    ColorSimplifySettings::default().with_cof_dimension_invariants(),
                )
                .expand_dots()?;

            // Exact production branches have already mapped every owned
            // numerator fragment. The former coarse export sign replacement
            // has no role here and must not remap these factors a second time.
            Ok(atom
                .replace(GS.m_uv_expansion)
                .with(GS.m_uv_vacuum)
                .replace(GS.dim_epsilon)
                .with(0))
        })?;
        Ok(FinalIntegrands(simplified))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        cff::esurface::RaisedEsurfaceGroup,
        dot,
        graph::{cuts::CutSet, parse::IntoGraph},
        initialisation::test_initialise,
        settings::global::OrientationPattern,
        uv::{
            ApproximationType, Spinney,
            approx::{
                OrientationProjection, Rooted, UVCtx,
                local_4d::{Full4dCts, Local4dCts, uv_limit},
                projected_4d::{Projected4dApproximation, Projected4dSector},
            },
            hedge_poset::OwnedForestNode,
        },
    };
    use linnet::half_edge::subgraph::InternalSubGraph;

    #[test]
    fn projected_zero_sectors_preserve_cut_orders_without_energy_maps() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph projected_zero {
            edge [num=1 mass=1];
            node [num=1];
            a -> b [id=0 lmb_id=0];
            a -> b [id=1];
        })?;
        let current = OwnedForestNode {
            spinney: Spinney::new(
                InternalSubGraph::cleaned_filter_optimist(graph.full_filter(), graph.as_ref()),
                &graph,
                &graph.loop_momentum_basis,
            )
            .expect("the bubble has a compatible UV spinney"),
            topo_order: 1,
        };
        let mut cutset = CutSet::empty(graph.n_hedges());
        // This is a key-shape diagnostic: a typed zero needs neither actual
        // threshold surfaces nor a production map to preserve two cut orders.
        cutset.residue_selector.left_th_cut = Some(RaisedEsurfaceGroup {
            esurface_ids: Vec::new(),
            max_occurence: 2,
        });
        let production = Default::default();
        let pattern = OrientationPattern::default();
        let options = graph.denominator_only_cff_3d_expression_options();
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &pattern, true),
        );
        let settings = UVgenerationSettings {
            local_uv_cts_from_expanded_4d_integrands: true,
            ..Default::default()
        };
        let builder = FinalIntegrandBuilder::new(localizer, &settings);
        let local = Projected4dCts::new(vec![Projected4dSector {
            coefficient: Atom::Zero,
            frozen_factor: Atom::var(GS.numerator_sampling_scale),
        }]);
        let finalized =
            builder.build_projected(&mut graph, &current, &local, &IntegratedCts::root())?;
        assert_eq!(
            finalized.into_integrands(),
            cutset
                .residue_selector
                .generate_allowed_keys()
                .into_iter()
                .map(|index| (index, Atom::Zero))
                .collect(),
        );
        // A disabled integrated pole prefix is a legitimate zero recursion
        // input. Its next Taylor operation prunes the zero FourDSector, so
        // exercise the actual producer before the final assembly boundary.
        let zero_prefix = Full4dCts::recursion_input(
            &Local4dCts::root(),
            &IntegratedCts::root(),
            ApproximationType::PolePart,
            false,
            current.lmb(),
        )?;
        let given = OwnedForestNode {
            spinney: Spinney::empty(&graph),
            topo_order: 0,
        };
        let zero_local = uv_limit(
            &zero_prefix,
            &UVCtx::new(&graph, &settings),
            &current,
            &given,
            &current,
            &given,
        )?;
        assert!(zero_local.atom().is_zero());
        assert!(zero_local.active_sectors().is_empty());
        let produced = Projected4dApproximation::new(localizer, &mut graph, &settings)
            .project_local_4d(&zero_local)?;
        assert_eq!(produced.sectors().len(), 1);
        assert!(produced.sectors()[0].coefficient.is_zero());
        assert_eq!(produced.sectors()[0].frozen_factor, Atom::one());
        assert_eq!(
            builder
                .build_projected(&mut graph, &current, &produced, &IntegratedCts::root())?
                .into_integrands(),
            cutset
                .residue_selector
                .generate_allowed_keys()
                .into_iter()
                .map(|index| (index, Atom::Zero))
                .collect(),
            "a pruned local zero must preserve every allowed cut order without energy maps",
        );
        let missing_maps = Projected4dApproximation::new(localizer, &mut graph, &settings)
            .project_local_4d(&Local4dCts::root())
            .expect_err("a nonzero local source still requires production maps");
        assert!(
            missing_maps
                .to_string()
                .contains("no production energy maps")
        );
        let error = builder
            .build_projected(
                &mut graph,
                &current,
                &Projected4dCts::new(Vec::new()),
                &IntegratedCts::root(),
            )
            .expect_err("absent sectors are different from a deliberate typed zero");
        assert!(error.to_string().contains("no active UV sectors"));
        Ok(())
    }
}
