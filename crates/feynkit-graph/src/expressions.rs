//! Symbolic graph traversal shared with runtime graph consumers.
//!
//! A selected paired edge contributes its local numerator and the numerators of
//! both incident vertices. Remaining vertices contribute only when their whole
//! crown is selected, allowing dummy half-edges to be absent. An excluded region
//! removes its edges and all vertices incident to that region.

use linnet::half_edge::{
    HedgeGraph, NodeIndex,
    involution::{EdgeIndex, HedgePair},
    subgraph::{ModifySubSet, SubGraphLike, SubSetLike, SubSetOps, subset::SubSet},
};
use spenso::structure::{
    dimension::Dimension,
    representation::{Minkowski, RepName},
};
use symbolica::atom::{Atom, AtomCore, Symbol};

/// Application data access is supplied by callers; the selection rules are shared.
pub trait GraphExpressions<E, V> {
    fn numerator_of<S: SubGraphLike + SubSetOps>(
        &self,
        subgraph: &S,
        without: &S,
        vertex_numerator: impl Fn(&V) -> Atom,
        edge_numerator: impl Fn(&E) -> Atom,
        is_dummy: impl Fn(&E) -> bool,
    ) -> Atom;

    /// Only edges paired inside the selection contribute a propagator. Powers
    /// may be negative or zero, and the factor callback is not run for zero powers.
    fn denominator_of<S: SubGraphLike, Error>(
        &self,
        subgraph: &S,
        propagator: impl Fn(EdgeIndex, &E) -> Result<Atom, Error>,
        edge_power: impl Fn(EdgeIndex, &E) -> isize,
    ) -> Result<Atom, Error>;
}

impl<E, V, H> GraphExpressions<E, V> for HedgeGraph<E, V, H> {
    fn numerator_of<S: SubGraphLike + SubSetOps>(
        &self,
        subgraph: &S,
        without: &S,
        vertex_numerator: impl Fn(&V) -> Atom,
        edge_numerator: impl Fn(&E) -> Atom,
        is_dummy: impl Fn(&E) -> bool,
    ) -> Atom {
        let mut numerator = Atom::one();
        let mut seen: SubSet<NodeIndex> = SubSet::empty(self.n_nodes());
        for (node, _, _) in self.iter_nodes_of(without) {
            seen.add(node);
        }
        let selected = subgraph.subtract(without);
        for (pair, _, edge) in self.iter_edges_of(&selected) {
            if let HedgePair::Paired { source, sink } = pair {
                for node in [self.node_id(source), self.node_id(sink)] {
                    if !seen[node] {
                        seen.add(node);
                        numerator *= vertex_numerator(&self[node]);
                    }
                }
                numerator *= edge_numerator(edge.data);
            }
        }
        let unseen = !seen;
        // From all the nodes not yet covered by paired edges, include those
        // included in the subgraph, ignoring dummies.
        for node in unseen.included_iter() {
            if self
                .iter_crown(node)
                .all(|hedge| subgraph.includes(&hedge) || is_dummy(&self[self[&hedge]]))
            {
                numerator *= vertex_numerator(&self[node]);
            }
        }
        numerator
    }

    fn denominator_of<S: SubGraphLike, Error>(
        &self,
        subgraph: &S,
        propagator: impl Fn(EdgeIndex, &E) -> Result<Atom, Error>,
        edge_power: impl Fn(EdgeIndex, &E) -> isize,
    ) -> Result<Atom, Error> {
        let mut denominator = Atom::one();
        for (pair, edge_id, edge) in self.iter_edges_of(subgraph) {
            if pair.is_paired() {
                let power = edge_power(edge_id, edge.data);
                if power != 0 {
                    denominator *= propagator(edge_id, edge.data)?.pow(power as i64);
                }
            }
        }
        Ok(denominator)
    }
}

/// Symbol families used by the canonical tagged quadratic propagator.
#[derive(Clone, Copy)]
pub struct PropagatorSymbols {
    pub momentum: Symbol,
    pub denominator: Symbol,
}

impl PropagatorSymbols {
    /// Keep the identity, uncontracted momentum and mass alongside the explicit
    /// quadratic expression so UV derivatives can act on the fourth argument.
    pub fn denominator(
        &self,
        edge: EdgeIndex,
        mass_squared: &Atom,
        dimension: impl Into<Dimension>,
    ) -> Atom {
        let momentum = self.momentum.call(edge.0);
        let quadratic = Minkowski {}
            .new_rep(dimension)
            .inner_product(&momentum, &momentum)
            - mass_squared;
        self.denominator
            .call_args([Atom::num(edge.0), momentum, mass_squared.clone(), quadratic])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use linnet::half_edge::{
        builder::HedgeGraphBuilder,
        involution::{Flow, Hedge, Orientation},
        subgraph::SuBitGraph,
    };

    fn graph() -> HedgeGraph<(i32, bool), i32> {
        let mut builder = HedgeGraphBuilder::new();
        let a = builder.add_node(2);
        let b = builder.add_node(3);
        let c = builder.add_node(5);
        builder.add_edge(a, b, (7, false), Orientation::Default);
        builder.add_edge(b, c, (11, false), Orientation::Default);
        builder.add_external_edge(a, (13, false), Orientation::Default, Flow::Sink);
        builder.add_external_edge(c, (17, true), Orientation::Default, Flow::Source);
        builder.build()
    }

    fn numerator(
        graph: &HedgeGraph<(i32, bool), i32>,
        selected: &SuBitGraph,
        without: &SuBitGraph,
    ) -> Atom {
        graph.numerator_of(
            selected,
            without,
            |n| Atom::num(*n),
            |e| Atom::num(e.0),
            |e| e.1,
        )
    }

    #[test]
    fn numerator_preserves_paired_boundary_and_exclusion_rules() {
        let graph = graph();
        let empty = graph.empty_subgraph();
        let full = graph.full_filter();
        assert_eq!(
            numerator(&graph, &full, &empty),
            Atom::num(2 * 3 * 5 * 7 * 11)
        );
        assert_eq!(numerator(&graph, &empty, &empty), Atom::one());
        assert_eq!(numerator(&graph, &full, &full), Atom::one());

        let mut first: SuBitGraph = graph.empty_subgraph();
        first.add(graph[&EdgeIndex(0)].1);
        assert_eq!(numerator(&graph, &first, &empty), Atom::num(2 * 3 * 7));
        assert_eq!(numerator(&graph, &full, &first), Atom::num(5 * 11));

        // Both selected boundary half-edges together cover the middle vertex.
        let mut middle: SuBitGraph = graph.empty_subgraph();
        middle.add(Hedge(1));
        middle.add(Hedge(2));
        assert_eq!(numerator(&graph, &middle, &empty), Atom::num(3));

        // A missing dummy does not prevent the remaining crown from contributing.
        let mut last: SuBitGraph = graph.empty_subgraph();
        last.add(Hedge(3));
        assert_eq!(numerator(&graph, &last, &empty), Atom::num(5));
    }

    #[test]
    fn denominator_preserves_signed_powers_and_ignores_boundaries() {
        let graph = graph();
        let factor = |_: EdgeIndex, edge: &(i32, bool)| Ok::<_, ()>(Atom::num(edge.0));
        let denominator = graph
            .denominator_of(&graph.full_filter(), factor, |edge, _| {
                if edge == EdgeIndex(0) { -2 } else { 3 }
            })
            .unwrap();
        assert_eq!(denominator, Atom::num(11).pow(3) / Atom::num(7).pow(2));
        let mut split: SuBitGraph = graph.empty_subgraph();
        split.add(Hedge(0));
        assert_eq!(
            graph.denominator_of(&split, factor, |_, _| 1).unwrap(),
            Atom::one()
        );
        assert_eq!(
            graph.denominator_of(&graph.full_filter(), |_, _| Err::<Atom, _>(()), |_, _| 0),
            Ok(Atom::one())
        );
    }
}

impl crate::FeynmanDiagram {
    /// Construct the selected local numerator using the same boundary and
    /// exclusion rules as GammaLoop. Diagram-wide factors remain separate.
    pub fn numerator_of<S: SubGraphLike + SubSetOps>(&self, subgraph: &S, without: &S) -> Atom {
        self.graph.numerator_of(
            subgraph,
            without,
            |vertex| vertex.numerator.clone(),
            |edge| edge.numerator.clone(),
            |edge| edge.is_dummy,
        )
    }

    pub fn denominator_of<S: SubGraphLike>(
        &self,
        subgraph: &S,
        edge_powers: &std::collections::BTreeMap<crate::EdgeId, isize>,
    ) -> Result<Atom, crate::DiagramError> {
        self.denominator_of_in_dimension(subgraph, edge_powers, crate::symbols::dimension().into())
    }

    /// Symbolic dimension and signed propagator powers follow the runtime UV
    /// representation, retaining its four-argument `denom` wrapper.
    pub fn denominator_of_in_dimension<S: SubGraphLike>(
        &self,
        subgraph: &S,
        edge_powers: &std::collections::BTreeMap<crate::EdgeId, isize>,
        dimension: Dimension,
    ) -> Result<Atom, crate::DiagramError> {
        let symbols = PropagatorSymbols {
            momentum: crate::symbols::momentum(),
            denominator: crate::symbols::denominator(),
        };
        self.graph.denominator_of(
            subgraph,
            |edge_id, edge| {
                let particle = self.model.particle_by_id(edge.particle)?;
                let mass = particle
                    .symbolic_mass(&self.model)
                    .replace(symbolica::symbol!("UFO::ZERO"))
                    .with(Atom::Zero);
                Ok(symbols.denominator(edge_id, &mass.pow(2), dimension))
            },
            |edge_id, _| {
                edge_powers
                    .get(&crate::EdgeId(edge_id.0))
                    .copied()
                    .unwrap_or(1)
            },
        )
    }

    /// Local superficial degree for a half-edge selection, including its
    /// boundary vertices and excluding physical external momentum carriers.
    pub fn superficial_degree_of_divergence_of(
        &self,
        subgraph: &linnet::half_edge::subgraph::SuBitGraph,
        dimension: i32,
    ) -> Result<i32, crate::DiagramError> {
        use crate::DOD;
        let loops = self.momentum_basis_of(subgraph)?.loop_edges.len();
        let mut degree = i64::from(dimension) * loops as i64;
        for (_, _, vertex) in self.graph.iter_nodes_of(subgraph) {
            degree += i64::from(vertex.numerator.all_dod(crate::symbols::momentum())?);
        }
        for (pair, edge_id, edge) in self.graph.iter_edges_of(subgraph) {
            if pair.is_paired() && edge.data.external.is_none() && !edge.data.is_dummy {
                degree += i64::from(
                    edge.data
                        .numerator
                        .edge_dod(crate::symbols::momentum(), edge_id.0)?,
                ) - 2;
            }
        }
        i32::try_from(degree).map_err(|_| {
            crate::DiagramError::UvPowerCounting("degree exceeds the integer range".into())
        })
    }
}
