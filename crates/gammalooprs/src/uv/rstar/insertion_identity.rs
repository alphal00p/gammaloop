//! Dressed bare skeletons and transport between their insertion occurrences.

use std::{collections::BTreeSet, sync::LazyLock};

use eyre::{Result, ensure};
use linnet::half_edge::{
    NodeIndex,
    subgraph::{SuBitGraph, SubSetOps},
};
use symbolica::{
    atom::{Atom, Symbol},
    function, symbol,
};

use crate::{
    graph::{FeynmanGraph, Graph},
    model::Model,
    uv::counterterm_identity::{CountertermIdentity, CountertermKey},
};

static INSERTION_HOLE: LazyLock<Symbol> =
    LazyLock::new(|| symbol!("rstar_insertion_skeleton::hole"));

/// Equality includes the dressed skeleton, global tensors and hole loop orders.
/// Insertion coefficients and the separate `Graph::overall_factor` are excluded;
/// all scalar and tensor parts of `global_prefactor` remain part of the key.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct InsertionSkeletonKey(CountertermKey);

/// One chosen isomorphism into a common dressed skeleton.
///
/// Hole vertices have unit numerators. Their coefficients remain separate, so
/// they can be summed without changing the skeleton's equivalence certificate.
/// A single isomorphism transports an occurrence; automorphisms introduce no
/// multiplicity or symmetry factor.
pub struct InsertionSkeletonIdentity {
    key: InsertionSkeletonKey,
    identity: CountertermIdentity,
}

impl InsertionSkeletonIdentity {
    pub fn new(graph: &Graph, model: &Model, holes: &[(NodeIndex, usize)]) -> Result<Self> {
        let mut decorated = graph.clone();
        let vertices = graph
            .underlying
            .iter_nodes()
            .map(|(vertex, _, _)| vertex)
            .collect::<BTreeSet<_>>();
        let mut seen = BTreeSet::new();
        let hole_symbol = *INSERTION_HOLE;
        for &(vertex, loops) in holes {
            ensure!(vertices.contains(&vertex), "insertion hole has no vertex");
            ensure!(seen.insert(vertex), "duplicate insertion hole vertex");
            ensure!(loops > 0, "insertion hole must contain a loop");
            ensure!(
                decorated.underlying[vertex].num.value.is_one(),
                "insertion skeleton still contains a hole coefficient"
            );
            decorated.underlying[vertex].num.value = function!(hole_symbol, loops as i64);
        }
        let support = decorated
            .no_dummy()
            .subtract(&decorated.external_filter::<SuBitGraph>());
        let identity = CountertermIdentity::for_bare_skeleton(&decorated, &support, model)?;
        Ok(Self {
            key: InsertionSkeletonKey(identity.key.clone()),
            identity,
        })
    }

    pub fn key(&self) -> &InsertionSkeletonKey {
        &self.key
    }

    /// Transport one complete correlated numerator into this representative.
    ///
    /// All retained internal and external momenta are mapped with orientation
    /// signs. Open indices keep their typed skeleton ports, while contracted
    /// indices receive fresh names. Supplying the complete product preserves
    /// correlations between several insertion sites; this operation does not
    /// form products of separately summed families.
    pub fn transport_from(&self, source: &Self, expression: &Atom) -> Result<Atom> {
        ensure!(self.key == source.key, "insertion skeletons do not match");
        let canonical = source.identity.canonicalize_delta(expression)?;
        self.identity.restore_delta(&canonical)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        finalized_runtime_dot,
        graph::{Autogen, parse::IntoFinalizedRuntimeGraph},
        initialisation::test_initialise,
        numerator::aind::Aind,
        utils::{GS, load_generic_model},
    };
    use linnet::half_edge::involution::{EdgeIndex, HedgePair};
    use spenso::{
        network::library::symbolic::ETS,
        structure::representation::{Minkowski, RepName},
    };
    use symbolica::atom::AtomCore;

    fn bubble(model: &Model, permuted: bool) -> Result<(Graph, NodeIndex, Aind)> {
        let mut graph: Graph = if permuted {
            finalized_runtime_dot!(
                digraph G {
                    graph [projector=1]; node [num=1]; edge [particle="scalar_1" num=1];
                    e [style=invis]; B; A;
                    A:0 -> e [id=0 source="{ufo_order:0}"];
                    e -> B:1 [id=1 sink="{ufo_order:0}"];
                    B -> A [id=3 source="{ufo_order:1}" sink="{ufo_order:1}"];
                    A -> B [id=2 lmb_id=0 source="{ufo_order:2}" sink="{ufo_order:2}"];
                }, model
            )?
        } else {
            finalized_runtime_dot!(
                digraph G {
                    graph [projector=1]; node [num=1]; edge [particle="scalar_1" num=1];
                    e [style=invis]; A; B;
                    e -> A:0 [id=2 sink="{ufo_order:0}"];
                    B:1 -> e [id=3 source="{ufo_order:0}"];
                    A -> B [id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
                    A -> B [id=1 lmb_id=0 source="{ufo_order:2}" sink="{ufo_order:2}"];
                }, model
            )?
        };
        let edges = graph
            .iter_edges()
            .map(|(_, edge, _)| edge)
            .collect::<Vec<_>>();
        for edge in edges {
            let data = &mut graph.underlying[edge];
            data.particle = data.particle.clone().override_mass(Some(Atom::Zero));
        }
        graph.underlying[EdgeIndex::from(if permuted { 2 } else { 1 })]
            .extra_data
            .vakint_edge_power = Some(2);
        let external_vertex = |edge: usize| {
            graph
                .iter_edges()
                .find_map(|(pair, id, _)| {
                    if usize::from(id) == edge
                        && let HedgePair::Unpaired { hedge, .. } = pair
                    {
                        Some((graph.underlying.node_id(hedge), hedge))
                    } else {
                        None
                    }
                })
                .unwrap()
        };
        let (hole, _) = external_vertex(if permuted { 0 } else { 2 });
        let (outside, hedge) = external_vertex(if permuted { 1 } else { 3 });
        let port = Aind::Hedge(hedge.0 as u16, 1);
        let mink = Minkowski {}.new_rep(GS.dim);
        graph.underlying[outside].num = Autogen::explicit(function!(
            ETS.metric,
            mink.to_symbolic([Atom::from(port)]),
            mink.to_symbolic([Atom::from(Aind::Dummy(if permuted { 900 } else { 400 }))])
        ));
        Ok((graph, hole, port))
    }

    #[test]
    fn skeleton_transports_internal_momentum_reversal_and_renamed_ports() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("scalars");
        let (first, first_hole, first_port) = bubble(&model, false)?;
        let (second, second_hole, second_port) = bubble(&model, true)?;
        let first = InsertionSkeletonIdentity::new(&first, &model, &[(first_hole, 1)])?;
        let second = InsertionSkeletonIdentity::new(&second, &model, &[(second_hole, 1)])?;
        assert_eq!(first.key(), second.key());
        let mink = Minkowski {}.new_rep(GS.dim);
        let numerator = function!(GS.emr_mom, 0, mink.to_symbolic([Atom::from(first_port)]));
        assert_eq!(
            second.transport_from(&first, &numerator)?,
            -function!(GS.emr_mom, 3, mink.to_symbolic([Atom::from(second_port)]))
        );
        let loop_coordinate = function!(
            feynkit_graph::symbols::loop_momentum(),
            0,
            mink.to_symbolic([Atom::from(first_port)])
        );
        assert_eq!(
            second.transport_from(&first, &loop_coordinate)?,
            function!(GS.emr_mom, 2, mink.to_symbolic([Atom::from(second_port)]))
        );
        let index = mink.to_symbolic([Atom::from(Aind::Dummy(99001))]);
        let contracted = function!(GS.emr_mom, 0, index).pow(2);
        assert_ne!(
            second.transport_from(&first, &contracted)?,
            second.transport_from(&first, &contracted)?,
            "separate transported products must not share private contractions"
        );
        Ok(())
    }

    #[test]
    fn skeleton_distinguishes_hole_order_powers_and_outside_rules() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("scalars");
        let (graph, hole, _) = bubble(&model, false)?;
        let key = InsertionSkeletonIdentity::new(&graph, &model, &[(hole, 1)])?;
        let mut differently_weighted = graph.clone();
        differently_weighted.overall_factor *= Atom::num(-1);
        assert_eq!(
            key.key(),
            InsertionSkeletonIdentity::new(&differently_weighted, &model, &[(hole, 1)])?.key()
        );
        assert_ne!(
            key.key(),
            InsertionSkeletonIdentity::new(&graph, &model, &[(hole, 2)])?.key()
        );
        let mut changed = graph.clone();
        changed.underlying[EdgeIndex::from(0)]
            .extra_data
            .vakint_edge_power = Some(3);
        assert_ne!(
            key.key(),
            InsertionSkeletonIdentity::new(&changed, &model, &[(hole, 1)])?.key()
        );
        let mut changed = graph.clone();
        changed.underlying[EdgeIndex::from(0)].num.value *= Atom::one() + Atom::var(GS.dim_epsilon);
        assert_ne!(
            key.key(),
            InsertionSkeletonIdentity::new(&changed, &model, &[(hole, 1)])?.key()
        );
        let mut changed = graph.clone();
        changed.global_prefactor.projector *= Atom::num(2);
        assert_ne!(
            key.key(),
            InsertionSkeletonIdentity::new(&changed, &model, &[(hole, 1)])?.key()
        );
        assert!(InsertionSkeletonIdentity::new(&graph, &model, &[(hole, 0)]).is_err());
        assert!(InsertionSkeletonIdentity::new(&graph, &model, &[(hole, 1), (hole, 1)]).is_err());
        changed.underlying[hole].num.value = Atom::num(2);
        assert!(InsertionSkeletonIdentity::new(&changed, &model, &[(hole, 1)]).is_err());
        Ok(())
    }

    #[test]
    fn two_holes_transport_typed_ports_exposed_by_global_tensors() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("scalars");
        let mink = Minkowski {}.new_rep(GS.dim);
        let color = idenso::representations::ColorAdjoint {}.new_rep(8);
        let prepare = |permuted: bool| -> Result<_> {
            let (mut graph, first_hole, lorentz) = bubble(&model, permuted)?;
            let Aind::Hedge(owner, _) = lorentz else {
                unreachable!()
            };
            let second_hole = graph
                .underlying
                .node_id(linnet::half_edge::involution::Hedge(usize::from(owner)));
            graph.underlying[second_hole].num.value = Atom::one();
            let color_port = Aind::Hedge(owner, 7);
            let color_partner = Aind::Dummy(if permuted { 901 } else { 401 });
            graph.global_prefactor.num = function!(
                ETS.metric,
                mink.to_symbolic([Atom::from(lorentz)]),
                mink.to_symbolic([Atom::from(Aind::Dummy(if permuted { 900 } else { 400 }))])
            );
            graph.global_prefactor.projector = function!(
                ETS.metric,
                color.to_symbolic([Atom::from(color_port)]),
                color.to_symbolic([Atom::from(color_partner)])
            );
            let identity = InsertionSkeletonIdentity::new(
                &graph,
                &model,
                &[(first_hole, 1), (second_hole, 1)],
            )?;
            Ok((identity, lorentz, color_port, color_partner))
        };
        let (first, first_lorentz, first_color, first_partner) = prepare(false)?;
        let (second, second_lorentz, second_color, second_partner) = prepare(true)?;
        assert_eq!(first.key(), second.key());
        let value = |edge, lorentz, color_port, color_partner| {
            function!(GS.emr_mom, edge, mink.to_symbolic([Atom::from(lorentz)]))
                * function!(
                    ETS.metric,
                    color.to_symbolic([Atom::from(color_port)]),
                    color.to_symbolic([Atom::from(color_partner)])
                )
        };
        assert_eq!(
            second.transport_from(&first, &value(0, first_lorentz, first_color, first_partner))?,
            -value(3, second_lorentz, second_color, second_partner)
        );
        Ok(())
    }

    #[test]
    fn ordinary_counterterm_identity_keeps_boundary_only_contract() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("scalars");
        let (graph, _, port) = bubble(&model, false)?;
        let support = graph
            .no_dummy()
            .subtract(&graph.external_filter::<SuBitGraph>());
        let ordinary = CountertermIdentity::new(&graph, &support, &model)?;
        let mut changed = graph.clone();
        changed.global_prefactor.projector *= Atom::num(2);
        let other = CountertermIdentity::new(&changed, &support, &model)?;
        assert_eq!(ordinary.key, other.key);
        let slot = Minkowski {}.new_rep(GS.dim).to_symbolic([Atom::from(port)]);
        let internal = function!(GS.emr_mom, 0, &slot);
        let error = ordinary
            .canonicalize_delta(&internal)
            .expect_err("ordinary counterterms must reject surviving loop momenta");
        assert!(
            error.to_string().contains("non-boundary momentum"),
            "a valid vector must reach the momentum-boundary check: {error:#}"
        );
        let external = function!(GS.emr_mom, 2, &slot);
        assert_eq!(
            ordinary.restore_delta(&ordinary.canonicalize_delta(&external)?)?,
            external
        );
        Ok(())
    }
}
