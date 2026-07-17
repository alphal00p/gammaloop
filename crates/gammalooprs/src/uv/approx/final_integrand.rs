use crate::{
    cff::CutCFFIndex,
    debug_tags,
    graph::Graph,
    utils::{GS, W_},
    uv::{
        Integrands, UVgenerationSettings, UltravioletGraph,
        approx::{
            ForestNodeLike,
            integrated::IntegratedCts,
            local_3d::{Local3DCts, Localizer},
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
use linnet::half_edge::subgraph::SubSetOps;
use symbolica::{
    atom::{Atom, AtomCore},
    function,
};
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct FinalIntegrands(Integrands);

impl FinalIntegrands {
    pub(crate) fn iter(&self) -> impl Iterator<Item = (&CutCFFIndex, &Atom)> {
        self.0.iter()
    }
}

pub(crate) struct FinalIntegrandBuilder<'a> {
    localizer: Localizer<'a>,
    marker: UvMarker,
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
        let global_num = graph.global_atom();
        debug_tags!(#generation, #profile, #uv, #graph, #summary;
            global_num = %global_num.log_display(),
            "Computed global numerator"
        );

        let localized_integrated: Integrands = self
            .localizer
            .localize(
                &integrated.physical_finite_counterterm_atom(),
                graph,
                current,
            )?
            .combine()?;

        let reduced = graph
            .full_filter()
            .subtract(current.subgraph())
            .subtract(&graph.initial_state_cut);

        let full_graph = graph.full_filter();
        let localized_integrated = localized_integrated
            .map(|atom| self.marker.prefix(&full_graph, current.subgraph(), atom));
        let local_terms = local_terms
            .integrands()
            .map(|atom| self.marker.prefix(&full_graph, current.subgraph(), atom));
        let final_int = localized_integrated.zip_add(&local_terms)?;
        let mut resnum = graph
            .numerator(&reduced, current.subgraph())
            .get_single_atom()
            .unwrap();
        let bridgeless_reduced = reduced.subtract(&graph.tree_edges);

        let mut reps = Vec::new();
        for (p, eid, _) in graph.as_ref().iter_edges_of(&bridgeless_reduced) {
            if p.is_paired() {
                reps.push(GS.add_parametric_sign(eid));
            }
        }

        resnum = resnum.replace_multiple(&reps) * &global_num;
        Ok(FinalIntegrands(final_int.fallible_map(|a| {
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

            a *= &resnum;

            let color_simplify_input = a.replace(GS.dim).with(4);

            a = color_simplify_input
                .collect_factors()
                .simplify_metrics()
                .simplify_color_with(
                    ColorSimplifySettings::default().with_cof_dimension_invariants(),
                );

            a = a.expand_dots()?;

            a = a.replace_multiple(&reps);

            Ok(a.replace(GS.m_uv_expansion)
                .with(GS.m_uv_vacuum)
                .replace(GS.dim_epsilon)
                .with(0))
        })?))
    }
}
