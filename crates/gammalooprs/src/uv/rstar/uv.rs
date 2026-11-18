//! Lower-loop UV maps evaluated by the existing massive R operation.

use std::sync::Arc;

use eyre::{Result, WrapErr, ensure, eyre};
use linnet::half_edge::{
    involution::EdgeIndex,
    subgraph::{Inclusion, InternalSubGraph, SuBitGraph, SubSetOps},
};
use symbolica::atom::{Atom, AtomCore};
use vakint::Vakint;

use crate::{
    graph::{Autogen, Graph, LMBext},
    model::Model,
    utils::GS,
    uv::{
        ApproximationType, RenormalizationPrescriptionSettings, UVgenerationSettings,
        UltravioletGraph, approx::CutStructure, counterterm_cache::CountertermCache,
        hedge_poset::Wood,
    },
};

use super::{CountertermCall, GraphRegion};

/// Adapter from a scoped R* counterterm argument to the established UV forest.
///
/// Contractions retain original edge and half-edge labels, so local insertions
/// keep their momentum and tensor attachment identities. A decorated graph's
/// complete numerator participates in its cache identity: support alone cannot
/// distinguish different insertion polynomials.
pub struct MassiveUvEvaluator<'a> {
    graph: &'a Graph,
    model: &'a Model,
    vakint: &'a Vakint,
    settings: UVgenerationSettings,
    reference_mass: Option<(EdgeIndex, Atom)>,
    counterterm_cache: Option<Arc<CountertermCache>>,
}

impl<'a> MassiveUvEvaluator<'a> {
    pub fn new(
        graph: &'a Graph,
        model: &'a Model,
        vakint: &'a Vakint,
        settings: &UVgenerationSettings,
    ) -> Self {
        let settings = UVgenerationSettings {
            generate_integrated: true,
            subtract_uv: true,
            add_marker: false,
            renormalization_prescription: RenormalizationPrescriptionSettings {
                log_divergent: ApproximationType::PolePart,
                massive_power_divergent: ApproximationType::PolePart,
                massless_power_divergent: ApproximationType::PolePart,
                overrides: Vec::new(),
            },
            ..settings.clone()
        };
        Self {
            graph,
            model,
            vakint,
            settings,
            reference_mass: None,
            counterterm_cache: None,
        }
    }

    /// Share completed decorated UV maps between parent graphs. Insertions
    /// participate in graph identity before lookup, including their epsilon
    /// dependence and open tensor structure.
    pub fn with_counterterm_cache(mut self, cache: Arc<CountertermCache>) -> Self {
        self.counterterm_cache = Some(cache);
        self
    }

    /// Use the parent's reference-family mass on this edge when it survives.
    ///
    /// `mass` is the mass itself. The regular denominator constructor squares
    /// it. Inserted polynomials may retain this parameter after contraction;
    /// those occurrences are never replaced by a new marked propagator.
    pub fn with_reference_mass(mut self, edge: EdgeIndex, mass: Atom) -> Self {
        self.reference_mass = Some((edge, mass));
        self
    }

    /// Evaluate the complete UV counterterm of the decorated argument.
    ///
    /// `insertions` follows `call.argument.insertions`. For an insertion-monomial
    /// scope, call this separately for each polynomial basis element and keep
    /// its Laurent coefficient and hard loops outside this call. Within the
    /// resulting UV forest, child poles stay inside their parent projection.
    /// No graph-wide prefactor or external projector is included.
    pub fn counterterm(&self, call: &CountertermCall, insertions: &[Atom]) -> Result<Atom> {
        ensure!(
            insertions.len() == call.argument.insertions.len(),
            "each contracted R* support needs exactly one insertion polynomial"
        );
        if call.argument.loops == 0 || insertions.iter().any(AtomCore::is_zero) {
            return Ok(Atom::Zero);
        }

        let graph = self.decorated_graph(&call.argument, insertions)?;
        let support = &call.argument.remaining;
        if let Some(cache) = &self.counterterm_cache {
            return Ok(cache
                .counterterm(&graph, support, self.model, self.vakint, &self.settings)?
                .replace(GS.m_uv_expansion)
                .with(GS.m_uv_vacuum));
        }
        if graph.compute_dod(support) < 0 {
            return Ok(Atom::Zero);
        }
        let catalogue = graph.classified_spinneys(
            support,
            self.model,
            &self.settings,
            &graph.loop_momentum_basis,
        );
        ensure!(
            catalogue.iter().any(|spinney| spinney.filter() == support),
            "the decorated UV argument has no overall spinney with a compatible momentum basis"
        );
        let mut forest = Wood::from_spinneys(
            catalogue,
            &graph,
            CutStructure::empty(&graph),
            &self.settings.vakint,
        )
        .unfold();
        forest.integrate(&graph, self.model, self.vakint, &self.settings)?;
        // Complete the regular R operation's auxiliary-mass identification,
        // just as its terminal renormalization sum does. Parent reference
        // masses are distinct parameters and remain in the local polynomial.
        Ok(forest
            .counterterm_for_support(support)?
            .replace(GS.m_uv_expansion)
            .with(GS.m_uv_vacuum))
    }

    fn decorated_graph(&self, region: &GraphRegion, insertions: &[Atom]) -> Result<Graph> {
        let mut graph = self.graph.clone();
        let mut contracted: SuBitGraph = graph.empty_subgraph();
        for slot in &region.insertions {
            ensure!(
                region.support.includes(&slot.support),
                "an R* insertion lies outside its counterterm support"
            );
            contracted.union_with(&slot.support);
        }
        ensure!(
            region.support.subtract(&contracted) == region.remaining,
            "R* remaining edges do not match its contracted supports"
        );
        let contracted = InternalSubGraph::try_new(contracted, &graph.underlying)
            .ok_or_else(|| eyre!("an R* contracted support contains an unpaired half-edge"))?;
        // Route before identifying nodes, using the original boundary rather
        // than treating the retained, consumed edges as new external carriers.
        graph.loop_momentum_basis = super::routing::remainder_basis(
            &self.graph,
            &region.support,
            &contracted,
            region.boundary.clone(),
            self.reference_mass.as_ref().map(|(edge, _)| *edge),
        )
        .wrap_err("could not route the decorated UV quotient")?;

        for (slot, polynomial) in region.insertions.iter().zip(insertions) {
            let nodes = graph
                .iter_nodes_of(&slot.support)
                .map(|(node, _, _)| node)
                .collect::<Vec<_>>();
            let first = nodes
                .first()
                .ok_or_else(|| eyre!("an R* insertion support contains no vertices"))?;
            let mut vertex = graph.underlying[*first].clone();
            vertex.num = Autogen::explicit(polynomial.clone());
            vertex.vertex_rule = None;
            // Spinney enumeration first uses this structural upper bound,
            // then computes the actual degree from the explicit polynomial.
            vertex.dod = Autogen::explicit(self.graph.local_dod(&slot.support).max(0));
            graph.underlying.identify_nodes(&nodes, vertex);
        }
        // Historical vertices still carry their old numerators until removed.
        // Compacting vertices preserves every edge and half-edge identity.
        graph.underlying.forget_identification_history();
        let consumed_edges = graph
            .iter_edges_of(&contracted)
            .map(|(_, edge, _)| edge)
            .collect::<Vec<_>>();
        for edge in consumed_edges {
            graph.underlying[edge].is_dummy = true;
        }
        if let Some((edge, mass)) = &self.reference_mass
            && region.contains_edge(&graph, *edge)
        {
            let data = &mut graph.underlying[*edge];
            data.particle = data.particle.clone().override_mass(Some(mass.clone()));
        }
        ensure!(
            graph.n_loops(&region.remaining) == region.loops,
            "decorated UV quotient has an inconsistent remaining loop count"
        );
        Ok(graph)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        finalized_runtime_dot,
        graph::parse::IntoFinalizedRuntimeGraph,
        initialisation::test_initialise,
        utils::load_generic_model,
        uv::{
            CTRenormalizationRule,
            rstar::{CountertermScope, InsertionSlot, RStarPlan},
            uv_graph::UVE,
        },
    };
    use linnet::half_edge::subgraph::ModifySubSet;
    use spenso::{
        network::tags::SPENSO_TAG,
        structure::representation::{Minkowski, RepName},
    };
    use symbolica::{function, parse, symbol};

    #[test]
    fn convergent_decorated_quotient_caches_zero_without_aliasing_divergent_insertion() -> Result<()>
    {
        test_initialise()?;
        // Contracting the bubble leaves a one-loop triangle. Its constant
        // insertion is convergent; a quadratic insertion makes it logarithmic.
        let mut graph: Graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=1];
                node [num=1];
                edge [particle="scalar_1" num=1];
                A -> B [id=0 source="{ufo_order:0}" sink="{ufo_order:0}"];
                A -> B [id=1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
                B -> C [id=2 source="{ufo_order:2}" sink="{ufo_order:0}"];
                C -> D [id=3 source="{ufo_order:1}" sink="{ufo_order:0}"];
                D -> A [id=4 lmb_id=1 source="{ufo_order:1}" sink="{ufo_order:2}"];
            }, "scalars"
        )?;
        let first = graph
            .iter_nodes()
            .find_map(|(node, _, data)| (data.name.value == "A").then_some(node))
            .unwrap();
        let hard = function!(GS.emr_mom, 0, Minkowski {}.new_rep(GS.dim).to_symbolic([]));
        graph.underlying[first].num = Autogen::explicit(function!(SPENSO_TAG.dot, &hard, &hard));
        graph.underlying[first].dod = Autogen::explicit(2);
        let support = graph.full_filter();
        let mut contracted = graph.empty_subgraph::<SuBitGraph>();
        for edge in [0, 1] {
            contracted.add(graph[&EdgeIndex::from(edge)].1);
        }
        let call = CountertermCall {
            argument: GraphRegion {
                remaining: support.subtract(&contracted),
                insertions: vec![InsertionSlot {
                    boundary: graph.dummy_less_full_crown(&contracted),
                    support: contracted,
                }],
                boundary: graph.dummy_less_full_crown(&support),
                support,
                loops: 1,
            },
            taylor_order: 0,
            scope: CountertermScope::InsertionMonomial,
        };
        let model = load_generic_model("scalars");
        let settings = super::super::tests::settings();
        let vakint = Vakint::new()?;
        let cache = Arc::new(CountertermCache::default());
        let ordinary = MassiveUvEvaluator::new(&graph, &model, &vakint, &settings);
        let cached = MassiveUvEvaluator::new(&graph, &model, &vakint, &settings)
            .with_counterterm_cache(cache.clone());
        let constant = Atom::one();
        let quotient = ordinary.decorated_graph(&call.argument, std::slice::from_ref(&constant))?;
        assert_eq!(quotient.compute_dod(&call.argument.remaining), -2);
        assert!(
            ordinary
                .counterterm(&call, std::slice::from_ref(&constant))?
                .is_zero()
        );
        for _ in 0..2 {
            assert!(
                cached
                    .counterterm(&call, std::slice::from_ref(&constant))?
                    .is_zero()
            );
        }
        assert_eq!(cache.stats().misses, 1);
        assert_eq!(cache.stats().stores, 1);
        assert_eq!(cache.stats().hits, 1);

        let momentum = function!(GS.emr_mom, 2, Minkowski {}.new_rep(GS.dim).to_symbolic([]));
        let quadratic = function!(SPENSO_TAG.dot, &momentum, &momentum);
        let expected = ordinary.counterterm(&call, std::slice::from_ref(&quadratic))?;
        assert!(!expected.is_zero());
        let actual = cached.counterterm(&call, std::slice::from_ref(&quadratic))?;
        super::super::tests::assert_same_counterterm(&actual, &expected)?;
        let repeated = cached.counterterm(&call, std::slice::from_ref(&quadratic))?;
        super::super::tests::assert_same_counterterm(&repeated, &expected)?;
        assert_eq!(cache.stats().misses, 2);
        assert_eq!(cache.stats().stores, 2);
        assert_eq!(cache.stats().hits, 2);

        // A convergent overall support is absent from the catalogue. Its own
        // non-MS override must still be rejected before a zero can be stored.
        let mut unsupported = settings;
        unsupported
            .renormalization_prescription
            .overrides
            .push(CTRenormalizationRule::new(
                quotient.ct_identifier(&call.argument.remaining, &model),
                ApproximationType::MUV,
            ));
        let invalid_cache = CountertermCache::default();
        let error = invalid_cache
            .counterterm(
                &quotient,
                &call.argument.remaining,
                &model,
                &vakint,
                &unsupported,
            )
            .unwrap_err();
        assert!(error.to_string().contains("pole-part subtraction"));
        assert_eq!(invalid_cache.stats().stores, 0);
        Ok(())
    }

    #[test]
    fn decorated_quotient_preserves_labels_and_replaces_all_contracted_vertices() -> Result<()> {
        test_initialise()?;
        let graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=1];
                edge [particle="scalar_1" num=1];
                A [num=2];
                B [num=3];
                C [num=5];
                A -> B [id=0 source="{ufo_order:0}" sink="{ufo_order:0}"];
                A -> B [id=1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
                A -> C [id=2 source="{ufo_order:2}" sink="{ufo_order:0}"];
                C -> B [id=3 lmb_id=1 source="{ufo_order:1}" sink="{ufo_order:2}"];
            }, "scalars"
        )?;
        let model = load_generic_model("scalars");
        let settings = UVgenerationSettings::default();
        let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(0), &settings)?;
        let call = &plan.ir[0].forest[0].steps[0];
        let vakint = Vakint::new()?;
        let evaluator = MassiveUvEvaluator::new(&graph, &model, &vakint, &settings)
            .with_reference_mass(EdgeIndex::from(2), Atom::var(symbol!("mu")));
        let quotient = evaluator.decorated_graph(&call.argument, &[parse!("P+mu^2")])?;
        assert_eq!(quotient.n_hedges(), graph.n_hedges());
        assert_eq!(quotient.n_edges(), graph.n_edges());
        assert_eq!(quotient.n_nodes(), graph.n_nodes() - 1);
        assert_eq!(quotient.n_loops(&call.argument.remaining), 1);
        assert_eq!(quotient.loop_momentum_basis.loop_edges.len(), 1);
        let numerator = quotient
            .numerator(&call.argument.remaining, &quotient.empty_subgraph())
            .get_single_atom()
            .unwrap();
        assert_eq!(numerator.expand(), parse!("5*(P+mu^2)").expand());
        assert!(quotient.underlying[EdgeIndex::from(0)].is_dummy);
        assert!(quotient.underlying[EdgeIndex::from(1)].is_dummy);
        assert_eq!(
            quotient.underlying[EdgeIndex::from(2)].mass_atom(&model),
            parse!("mu")
        );
        Ok(())
    }
}
