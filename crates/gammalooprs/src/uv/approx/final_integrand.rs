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
            integrated::IntegratedCts,
            local_3d::{Local3DCts, Localizer},
            local_4d::Full4dCts,
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
use linnet::half_edge::subgraph::{SubSetLike, SubSetOps};
use symbolica::{
    atom::{Atom, AtomCore},
    function,
};
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
    project_from_4d: bool,
}

pub(crate) struct LocalizedIntegratedCt {
    pub active: Integrands,
    pub frozen_integrands: Integrands,
}

impl TryFrom<LocalizedIntegratedCt> for Integrands {
    type Error = eyre::Report;

    fn try_from(value: LocalizedIntegratedCt) -> Result<Self, Self::Error> {
        value.active.zip_mul(&value.frozen_integrands)
    }
}

impl<'a> FinalIntegrandBuilder<'a> {
    pub(crate) fn new(localizer: Localizer<'a>, settings: &UVgenerationSettings) -> Self {
        Self {
            localizer,
            marker: UvMarker::new(settings),
            project_from_4d: settings.local_uv_cts_from_expanded_4d_integrands,
        }
    }

    #[debug_instrument(
        graph = %graph.log_display(),
        current = %current.log_display(),
    )]
    pub(crate) fn build_3d<S: ForestNodeLike>(
        &self,
        graph: &mut Graph,
        current: &S,
        local_terms: &Local3DCts,
        integrated: &IntegratedCts,
    ) -> Result<FinalIntegrands> {
        let reduced = graph
            .full_filter()
            .subtract(current.subgraph())
            .subtract(&graph.initial_state_cut);
        let full_graph = graph.full_filter();

        // The empty forest is the original factorized 3D integrand. The
        // expanded-4D setting governs only proper UV approximations, so the
        // root must retain the ordinary production-orientation assembly.
        if self.project_from_4d && !current.subgraph().is_empty() {
            // Keep finite addbacks for nested multi-loop entries. Integrated
            // coefficients carry their forest-composition signs; projecting the
            // complete typed coefficient preserves the Tint(T(...)) terms. The
            // normalized spatial factor restores the integrated loop variables
            // without sending it through numerator energy-map substitution.
            let integrated_source = Full4dCts::from_frozen_coefficient(
                &integrated.physical_finite_counterterm_atom(),
                graph,
                &reduced,
                current.lmb(),
            );
            let projected_integrated =
                self.localizer
                    .project_4d(&integrated_source, graph, current.subgraph())?;
            let localized_integrated = projected_integrated
                .projected_integrands()?
                .map(|atom| self.marker.prefix(&full_graph, current.subgraph(), atom));
            let local_terms = local_terms
                .projected_integrands()?
                .map(|atom| self.marker.prefix(&full_graph, current.subgraph(), atom));
            let final_int = localized_integrated.zip_add(local_terms)?;

            let result = final_int.fallible_map(|a| {
                let mut a = a.clone();

                for (p, eid, _) in graph.as_ref().iter_edges_of(&reduced) {
                    let eid = usize::from(eid) as i64;
                    if p.is_paired() {
                        a = a
                            .replace(function!(GS.energy, eid))
                            .with(function!(GS.ose, eid));
                    }
                }

                a = a
                    .replace(function!(GS.ose, W_.mass_, W_.prop_))
                    .with(W_.prop_);

                let color_simplify_input = a.replace(GS.dim).with(4);

                a = color_simplify_input
                    .collect_factors()
                    .simplify_metrics()
                    .simplify_color_with(
                        ColorSimplifySettings::default().with_cof_dimension_invariants(),
                    );

                a = a.expand_dots()?;

                Ok(a.replace(GS.m_uv_expansion)
                    .with(GS.m_uv_vacuum)
                    .replace(GS.dim_epsilon)
                    .with(0))
            })?;

            return Ok(FinalIntegrands(result));
        }

        let global_num = graph.global_atom();
        debug_tags!(#generation, #profile, #uv, #graph, #summary;
            global_num = %global_num.log_display(),
            "Computed global numerator"
        );

        let localized_integrated = self
            .localizer
            .localize(
                &integrated.physical_finite_counterterm_atom(),
                graph,
                current,
            )?
            .combine()?
            .map(|atom| self.marker.prefix(&full_graph, current.subgraph(), atom));
        let local_terms = local_terms
            .integrands()
            .map(|atom| self.marker.prefix(&full_graph, current.subgraph(), atom));
        let final_int = localized_integrated.zip_add(&local_terms)?;
        let resnum = graph
            .numerator(&reduced, current.subgraph())
            .get_single_atom()
            .unwrap()
            * &global_num;
        let bridgeless_reduced = reduced.subtract(&graph.tree_edges);

        let coarse_energy_replacements = if self.localizer.uses_exact_maps() {
            None
        } else {
            let mut replacements = Vec::new();
            for (pair, edge_id, _) in graph.as_ref().iter_edges_of(&bridgeless_reduced) {
                if pair.is_paired() {
                    replacements.push(GS.add_parametric_sign(edge_id));
                }
            }
            Some(replacements)
        };

        let mut result: Option<Integrands> = None;
        for (orientation_id, source_edge_energy_map, integrands) in final_int.iter_orientations() {
            let mapped_resnum = self.localizer.map_numerator(
                graph,
                orientation_id,
                source_edge_energy_map,
                &resnum,
            )?;
            let mapped = integrands.fallible_map(|a| {
                let mut a = a.clone();

                for (p, eid, _) in graph.as_ref().iter_edges_of(&reduced) {
                    let eid = usize::from(eid) as i64;
                    if p.is_paired() {
                        a = a
                            .replace(function!(GS.energy, eid))
                            .with(function!(GS.ose, eid));
                    }
                }

                a = a
                    .replace(function!(GS.ose, W_.mass_, W_.prop_))
                    .with(W_.prop_);

                a *= &mapped_resnum;

                let color_simplify_input = a.replace(GS.dim).with(4);

                a = color_simplify_input
                    .collect_factors()
                    .simplify_metrics()
                    .simplify_color_with(
                        ColorSimplifySettings::default().with_cof_dimension_invariants(),
                    );

                a = a.expand_dots()?;

                // Coarse export projectors still introduce parametric signs at
                // this legacy boundary. Exact production branches have already
                // mapped every owned numerator fragment and must not fall back
                // to a second coarse replacement here.
                if let Some(replacements) = &coarse_energy_replacements {
                    a = a.replace_multiple(replacements);
                }

                Ok(a.replace(GS.m_uv_expansion)
                    .with(GS.m_uv_vacuum)
                    .replace(GS.dim_epsilon)
                    .with(0))
            })?;
            result = Some(match result {
                Some(sum) => sum.zip_add(mapped)?,
                None => mapped,
            });
        }

        Ok(FinalIntegrands(result.ok_or_else(|| {
            eyre::eyre!("final 3D UV integrand contains no production energy maps")
        })?))
    }
}
