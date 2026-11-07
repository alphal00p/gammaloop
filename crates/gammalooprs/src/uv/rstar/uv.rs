//! Lower-loop UV maps evaluated by the existing massive R operation.

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
        UltravioletGraph, approx::CutStructure, hedge_poset::Wood,
    },
};

use super::{CountertermCall, GraphRegion};

/// Adapter from a scoped R* counterterm argument to the established UV forest.
///
/// Contractions retain original edge and half-edge labels, so local insertions
/// keep their momentum and tensor attachment identities. Each decorated graph
/// receives its own forest: support alone cannot cache different polynomials.
pub struct MassiveUvEvaluator<'a> {
    graph: &'a Graph,
    model: &'a Model,
    vakint: &'a Vakint,
    settings: UVgenerationSettings,
    reference_mass: Option<(EdgeIndex, Atom)>,
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
        }
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
        graph.loop_momentum_basis = self
            .graph
            .shrunken_sub_lmb(&region.support, &contracted, region.boundary.clone())
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
        uv::{rstar::RStarPlan, uv_graph::UVE},
    };
    use symbolica::{parse, symbol};

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
