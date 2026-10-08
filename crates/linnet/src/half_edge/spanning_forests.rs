//! Spanning forest enumeration by edge contraction.

use super::{
    involution::{Flow, Hedge, Involution},
    nodestore::Checkpoint,
    subgraph::{Inclusion, ModifySubSet, SuBitGraph, SubGraphLike, SubSetLike, SubSetOps},
    tree::SimpleTraversalTree,
    HedgeGraph, NodeIndex,
};

impl<E, V, H, N: Checkpoint> HedgeGraph<E, V, H, N> {
    /// Gives all subgraphs corresponding to all the spanning trees of the graph.
    ///
    /// See Winter, Pawel, “An Algorithm for the Enumeration of Spanning Trees,”
    /// *BIT Numerical Mathematics* 26, no. 1 (1986): 44–62,
    /// <https://doi.org/10.1007/BF01939361>.
    pub fn all_spanning_forests_of<S: SubGraphLike>(&self, subgraph: &S) -> Vec<S::Base>
    where
        S::Base: SubSetLike<Base = S::Base> + Clone + ModifySubSet<Hedge>,
    {
        if subgraph.is_empty() {
            return vec![];
        }
        let mut visited_edges: SuBitGraph = self.empty_subgraph();

        let mut nodes = vec![];

        // Iterate over all edges in the subgraph
        for hedge_index in subgraph.included_iter() {
            if !visited_edges.includes(&hedge_index) {
                // Perform DFS to find all reachable edges from this edge
                let root_node = self.node_id(hedge_index);
                let tree =
                    SimpleTraversalTree::depth_first_traverse(self, subgraph, &root_node, None)
                        .unwrap();

                let reachable_edges = tree.covers(subgraph);

                visited_edges.union_with(&reachable_edges);
                let mut nodes_of_tree = tree.node_order();
                nodes_of_tree.remove(0);
                nodes.extend(nodes_of_tree);
            }
        }

        let mut forests = SpanningForests {
            involution: self.as_ref(),
            topology: self.node_store.topology().clone(),
            available: subgraph.included().clone(),
            removed: vec![],
            crowns: vec![],
            parallel: vec![],
            levels: vec![],
            digits: vec![],
            forests: vec![],
        };
        // An outer checkpoint keeps the journal allocated across branches.
        forests.topology.checkpoint();
        forests.contract(&nodes);
        forests.forests
    }
}

/// Winter's enumeration on a private copy of the topology: each branch
/// identifies the endpoints of an edge in place and rolls the identification
/// back afterwards, and the available half-edges are restored from a log.
struct SpanningForests<'a, N, S> {
    involution: &'a Involution,
    topology: N,
    /// The half-edges that may still join a forest; only edges with both halves
    /// available count.
    available: S,
    /// The half-edges taken out of `available`, restored when their level ends.
    removed: Vec<Hedge>,
    /// The crowns being branched over, one run per level.
    crowns: Vec<Hedge>,
    /// The source half-edges of the edges contracted at each level, ascending
    /// within a level; `levels` holds where each level's run starts.
    parallel: Vec<Hedge>,
    levels: Vec<usize>,
    /// The parallel edge chosen at each level while emitting forests.
    digits: Vec<usize>,
    forests: Vec<S>,
}

impl<N: Checkpoint, S: SubSetLike + ModifySubSet<Hedge>> SpanningForests<'_, N, S> {
    /// Contracts the last of `nodes` with each neighbour in turn. The forests of
    /// the contracted graph, completed by one of the edges joining the two
    /// nodes, are forests of the graph; those edges then leave the graph.
    fn contract(&mut self, nodes: &[NodeIndex]) {
        let Some((&last, rest)) = nodes.split_last() else {
            self.emit();
            return;
        };
        let node = self.topology.representative(last);
        let (crown, removed) = (self.crowns.len(), self.removed.len());
        self.crowns.extend(self.topology.crown(node));
        for at in crown..self.crowns.len() {
            let Some(neighbor) = self.neighbor(self.crowns[at]) else {
                continue;
            };
            if neighbor == node {
                // Self-loops join no forest.
                continue;
            }
            let level = self.parallel.len();
            self.levels.push(level);
            for later in at..self.crowns.len() {
                let hedge = self.crowns[later];
                if self.neighbor(hedge) == Some(neighbor) {
                    let twin = self.involution.inv(hedge);
                    for half in [hedge, twin] {
                        self.available.sub(half);
                        self.removed.push(half);
                    }
                    let source = match self.involution.flow(hedge) {
                        Flow::Source => hedge,
                        Flow::Sink => twin,
                    };
                    self.parallel.push(source);
                }
            }
            self.parallel[level..].sort_unstable();

            self.topology.checkpoint();
            self.topology.identify(&[node, neighbor]);
            self.contract(rest);
            self.topology.rollback();
            self.parallel.truncate(level);
            self.levels.pop();
        }
        for hedge in self.removed.drain(removed..) {
            self.available.add(hedge);
        }
        self.crowns.truncate(crown);
    }

    /// The node across the edge of `hedge` if both its halves are available.
    fn neighbor(&self, hedge: Hedge) -> Option<NodeIndex> {
        let twin = self.involution.inv(hedge);
        (twin != hedge && self.available.includes(&hedge) && self.available.includes(&twin))
            .then(|| self.topology.node_of(twin))
    }

    /// Adds a forest for every choice of one edge per level, the shallowest
    /// level varying fastest.
    fn emit(&mut self) {
        self.digits.clear();
        self.digits.resize(self.levels.len(), 0);
        loop {
            let mut forest = S::empty(self.available.size());
            for (&start, &digit) in self.levels.iter().zip(&self.digits) {
                let source = self.parallel[start + digit];
                forest.add(source);
                forest.add(self.involution.inv(source));
            }
            self.forests.push(forest);

            let mut level = 0;
            loop {
                let Some(digit) = self.digits.get_mut(level) else {
                    return;
                };
                let end = self
                    .levels
                    .get(level + 1)
                    .copied()
                    .unwrap_or(self.parallel.len());
                *digit += 1;
                if self.levels[level] + *digit < end {
                    break;
                }
                *digit = 0;
                level += 1;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use rand::{rngs::SmallRng, seq::SliceRandom, Rng, SeedableRng};

    use crate::half_edge::{
        builder::HedgeGraphBuilder,
        involution::{Flow, Hedge},
        remap::{HedgeRemap, NodeRemap},
        subgraph::{Inclusion, ModifySubSet, SuBitGraph, SubSetLike},
        HedgeGraph,
    };

    type Graph = HedgeGraph<(), (), ()>;

    /// Two nodes joined by a double edge hang off a triangle with a double
    /// edge, so the parallel edges of three contraction levels multiply, the
    /// first level's choice varying fastest. A self-loop and dangling
    /// half-edges never join a forest.
    #[test]
    fn forests_keep_their_order() {
        let mut builder = HedgeGraphBuilder::new();
        let nodes: Vec<_> = (0..4).map(|_| builder.add_node(())).collect();
        for (a, b) in [(0, 1), (0, 1), (1, 2), (2, 0), (1, 1), (2, 3), (2, 3)] {
            builder.add_edge(nodes[a], nodes[b], (), true);
        }
        builder.add_external_edge(nodes[0], (), true, Flow::Sink);
        builder.add_external_edge(nodes[3], (), true, Flow::Source);
        let graph: Graph = builder.build();
        assert_eq!(
            ordered_forests(&graph),
            [
                [0, 1, 4, 5, 10, 11],
                [2, 3, 4, 5, 10, 11],
                [0, 1, 4, 5, 12, 13],
                [2, 3, 4, 5, 12, 13],
                [0, 1, 6, 7, 10, 11],
                [2, 3, 6, 7, 10, 11],
                [0, 1, 6, 7, 12, 13],
                [2, 3, 6, 7, 12, 13],
                [4, 5, 6, 7, 10, 11],
                [4, 5, 6, 7, 12, 13],
            ]
        );
    }

    /// Within a level the parallel edges follow their source half-edges, not
    /// the crown being contracted: each double edge has opposite directions
    /// and is renumbered so that the order of its half-edges at either end
    /// disagrees with the order of its sources.
    #[test]
    fn parallel_edges_follow_their_sources() {
        let mut builder = HedgeGraphBuilder::new();
        let nodes: Vec<_> = (0..4).map(|_| builder.add_node(())).collect();
        for (a, b) in [(0, 1), (1, 0), (1, 2), (2, 0), (1, 1), (2, 3), (3, 2)] {
            builder.add_edge(nodes[a], nodes[b], (), true);
        }
        builder.add_external_edge(nodes[0], (), true, Flow::Sink);
        let mut graph: Graph = builder.build();
        // Half-edge `i` moves to `new[i]`: the sources of each double edge keep
        // their order and its sinks move around them.
        let new = vec![1, 3, 2, 0, 4, 5, 6, 7, 8, 9, 11, 13, 12, 10, 14];
        graph.permute_hedges(&HedgeRemap::bijection(new)).unwrap();
        assert_eq!(
            ordered_forests(&graph),
            [
                [1, 3, 4, 5, 11, 13],
                [0, 2, 4, 5, 11, 13],
                [1, 3, 4, 5, 10, 12],
                [0, 2, 4, 5, 10, 12],
                [1, 3, 6, 7, 11, 13],
                [0, 2, 6, 7, 11, 13],
                [1, 3, 6, 7, 10, 12],
                [0, 2, 6, 7, 10, 12],
                [4, 5, 6, 7, 11, 13],
                [4, 5, 6, 7, 10, 12],
            ]
        );
    }

    /// The half-edges of each spanning forest of `graph`, in enumeration order.
    fn ordered_forests(graph: &Graph) -> Vec<Vec<usize>> {
        let forests = graph.all_spanning_forests_of(&graph.full_filter());
        let hedges = |forest: &SuBitGraph| forest.included_iter().map(|hedge| hedge.0).collect();
        forests.iter().map(hedges).collect()
    }

    /// The spanning forests of random multigraphs with self-loops, dangling
    /// half-edges, pending identifications and partial subgraphs are the
    /// acyclic sets of full rank among the edges with both halves selected,
    /// each listed once.
    #[test]
    fn random_forests_match_brute_force() {
        let mut random = SmallRng::seed_from_u64(7);
        for _ in 0..500 {
            let mut builder = HedgeGraphBuilder::new();
            let nodes: Vec<_> = (0..1 + random.gen_range(0..6))
                .map(|_| builder.add_node(()))
                .collect();
            for _ in 0..random.gen_range(0..10) {
                let [a, b] = [(); 2].map(|_| nodes[random.gen_range(0..nodes.len())]);
                builder.add_edge(a, b, (), random.gen_bool(0.5));
            }
            for _ in 0..random.gen_range(0..3) {
                let node = nodes[random.gen_range(0..nodes.len())];
                builder.add_external_edge(node, (), true, Flow::Sink);
            }
            let mut graph: Graph = builder.build();
            // Shuffled numberings put sinks before sources and nodes out of
            // insertion order.
            let mut hedges: Vec<u32> = (0..graph.n_hedges() as u32).collect();
            hedges.shuffle(&mut random);
            graph
                .permute_hedges(&HedgeRemap::bijection(hedges))
                .unwrap();
            let mut nodes: Vec<u32> = (0..graph.node_bound() as u32).collect();
            nodes.shuffle(&mut random);
            graph.permute_nodes(&NodeRemap::bijection(nodes)).unwrap();
            if graph.n_hedges() > 0 && random.gen_bool(0.3) {
                let hedge = Hedge(random.gen_range(0..graph.n_hedges()));
                let ends = [hedge, graph.inv(hedge)].map(|hedge| graph.node_id(hedge));
                graph.identify_nodes(&ends, ());
            }
            let mut subgraph: SuBitGraph = graph.empty_subgraph();
            for hedge in (0..graph.n_hedges()).map(Hedge) {
                if random.gen_bool(0.8) {
                    subgraph.add(hedge);
                }
            }

            let mut forests = graph.all_spanning_forests_of(&subgraph);
            let n_forests = forests.len();
            forests.sort();
            forests.dedup();
            assert_eq!(forests.len(), n_forests, "forests repeat");
            assert_eq!(forests, brute_force(&graph, &subgraph));
        }
    }

    /// Every subset of the selected edges that is acyclic and of full rank,
    /// sorted; none for an empty subgraph.
    fn brute_force(graph: &Graph, subgraph: &SuBitGraph) -> Vec<SuBitGraph> {
        if subgraph.is_empty() {
            return vec![];
        }
        let edges: Vec<Hedge> = (0..graph.n_hedges())
            .map(Hedge)
            .filter(|&hedge| {
                let twin = graph.inv(hedge);
                hedge < twin && subgraph.includes(&hedge) && subgraph.includes(&twin)
            })
            .collect();
        // The number of union-find joins, and whether some edge closed a cycle.
        let joins = |chosen: &mut dyn Iterator<Item = Hedge>| {
            let mut parent: Vec<usize> = (0..graph.node_bound()).collect();
            let root = |parent: &[usize], mut node: usize| {
                while parent[node] != node {
                    node = parent[node];
                }
                node
            };
            let (mut joins, mut cyclic) = (0, false);
            for hedge in chosen {
                let [a, b] = [hedge, graph.inv(hedge)].map(|hedge| graph.node_id(hedge).0);
                let (a, b) = (root(&parent, a), root(&parent, b));
                if a == b {
                    cyclic = true;
                } else {
                    parent[a] = b;
                    joins += 1;
                }
            }
            (joins, cyclic)
        };
        let (rank, _) = joins(&mut edges.iter().copied());
        let mut forests: Vec<SuBitGraph> = (0..1usize << edges.len())
            .map(|chosen| {
                let mut forest: SuBitGraph = graph.empty_subgraph();
                for (_, &hedge) in
                    (edges.iter().enumerate()).filter(|(at, _)| chosen >> at & 1 == 1)
                {
                    forest.add(hedge);
                    forest.add(graph.inv(hedge));
                }
                forest
            })
            .filter(|forest| {
                let mut chosen = edges.iter().copied().filter(|hedge| forest.includes(hedge));
                joins(&mut chosen) == (rank, false)
            })
            .collect();
        forests.sort();
        forests
    }
}
