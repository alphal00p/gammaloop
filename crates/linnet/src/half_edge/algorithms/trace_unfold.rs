use std::cmp::Ordering;
use std::collections::VecDeque;
use std::fmt::{self, Display};
use std::hash::{Hash, Hasher};
use std::marker::PhantomData;
use std::ops::Deref;
use std::slice;

use indexmap::set::MutableValues;
use indexmap::IndexSet;

use crate::half_edge::builder::HedgeGraphBuilder;
use crate::half_edge::involution::{EdgeIndex, Flow};
use crate::half_edge::nodestore::{NodeStorageOps, NodeStorageVec};
use crate::half_edge::subgraph::SubSetLike;
use crate::half_edge::{HedgeGraph, NoData, NodeIndex, NodeVec};

/// An owned trace-unfolded DAG together with the source-graph node each unfolded node came from.
///
/// The wrapper keeps the structural graph and the `source_nodes` mapping together across
/// structure-preserving transformations such as [`UnfoldedTraceGraph::map`]. It also provides
/// traversal helpers that are meaningful only for trace-unfolded graphs.
pub struct UnfoldedTraceGraph<G, V, N: NodeStorageOps<NodeData = V> = NodeStorageVec<V>> {
    graph: HedgeGraph<EdgeIndex, V, NoData, N>,
    source_nodes: NodeVec<NodeIndex>,
    source_ty: PhantomData<fn() -> G>,
}

impl<G, V, N: NodeStorageOps<NodeData = V>> UnfoldedTraceGraph<G, V, N> {
    /// Creates a wrapper from an already-built unfolded graph and its source-node mapping.
    pub fn new(
        graph: HedgeGraph<EdgeIndex, V, NoData, N>,
        source_nodes: NodeVec<NodeIndex>,
    ) -> Self {
        Self {
            graph,
            source_nodes,
            source_ty: PhantomData,
        }
    }

    /// Maps node data while preserving the unfolded graph structure and source-node mapping.
    pub fn map<V2>(
        self,
        f: impl FnMut(&crate::half_edge::involution::Involution, NodeIndex, V) -> V2,
    ) -> UnfoldedTraceGraph<G, V2, N::OpStorage<V2>> {
        UnfoldedTraceGraph {
            graph: self.graph.map(f, |_, _, _, _, e| e, |_, h| h),
            source_nodes: self.source_nodes,
            source_ty: PhantomData,
        }
    }

    /// Returns the source-graph node from which an unfolded node originated.
    pub fn source_node(&self, node: NodeIndex) -> NodeIndex {
        self.source_nodes[node]
    }
}

impl<G, V, N: NodeStorageOps<NodeData = V>> AsRef<HedgeGraph<EdgeIndex, V, NoData, N>>
    for UnfoldedTraceGraph<G, V, N>
{
    fn as_ref(&self) -> &HedgeGraph<EdgeIndex, V, NoData, N> {
        &self.graph
    }
}

impl<G, V, N: NodeStorageOps<NodeData = V>> Deref for UnfoldedTraceGraph<G, V, N> {
    type Target = HedgeGraph<EdgeIndex, V, NoData, N>;

    fn deref(&self) -> &Self::Target {
        &self.graph
    }
}

/// Ops must be owned, hashable, and totally ordered for canonicalization.
pub trait Op: Clone + Eq + Ord {}
impl<T: Clone + Eq + Ord> Op for T {}

pub trait Independence<O> {
    fn independent(&self, a: &O, b: &O) -> bool;
}

pub trait Key<K> {
    fn key(&self, e: EdgeIndex) -> K;
}

/// Unfolds a graph of operations into a graph of trace keys.
///
/// The trait separates three concerns:
///
/// - [`TraceUnfold::key`] maps one traversed edge to its semantic operation label.
/// - [`Independence`] decides when two labels commute and may live in the same Foata level.
/// - [`TraceUnfold::realizations_for_op`] optionally expands one traversed edge into several
///   equivalent operation sequences before canonicalization.
///
/// In the simplest case, each edge has exactly one realization and the unfolding merges paths
/// that differ only by commuting independent operations.
///
/// # Examples
///
/// A minimal unfold where the two ways of reaching `AB` collapse to the same trace key:
///
/// ```
/// use linnet::half_edge::{
///     HedgeGraph, NoData, algorithms::trace_unfold::{HiddenData, Independence, TraceUnfold},
///     builder::HedgeGraphBuilder, involution::EdgeIndex, nodestore::DefaultNodeStore,
/// };
///
/// struct ExampleGraph {
///     graph: HedgeGraph<String, &'static str, NoData, DefaultNodeStore<&'static str>>,
/// }
///
/// impl Independence<HiddenData<String, EdgeIndex>> for ExampleGraph {
///     fn independent(
///         &self,
///         a: &HiddenData<String, EdgeIndex>,
///         b: &HiddenData<String, EdgeIndex>,
///     ) -> bool {
///         !a.order.chars().any(|c| b.order.contains(c))
///     }
/// }
///
/// impl TraceUnfold<String> for ExampleGraph {
///     type EdgeData = String;
///     type HedgeData = NoData;
///     type NodeData = &'static str;
///     type NodeStorage = DefaultNodeStore<&'static str>;
///
///     fn graph(
///         &self,
///     ) -> &HedgeGraph<Self::EdgeData, Self::NodeData, Self::HedgeData, Self::NodeStorage> {
///         &self.graph
///     }
///
///     fn key(&self, e: EdgeIndex) -> String {
///         self.graph[e].clone()
///     }
/// }
///
/// let mut builder = HedgeGraphBuilder::<String, &'static str>::new();
/// let empty = builder.add_node("empty");
/// let a = builder.add_node("A");
/// let b = builder.add_node("B");
/// let ab = builder.add_node("AB");
///
/// builder.add_edge(empty, a, "A".to_string(), true);
/// builder.add_edge(empty, b, "B".to_string(), true);
/// builder.add_edge(a, ab, "B".to_string(), true);
/// builder.add_edge(b, ab, "A".to_string(), true);
///
/// let graph = ExampleGraph {
///     graph: builder.build::<DefaultNodeStore<&'static str>>(),
/// };
/// let unfolded = graph.trace_unfold::<DefaultNodeStore<usize>>(empty);
/// let labels: Vec<_> = unfolded.iter_nodes().map(|(_, _, key)| key.to_string()).collect();
///
/// assert_eq!(unfolded.n_nodes(), 4);
/// assert!(labels.iter().any(|label| label == "∅"));
/// assert!(labels.iter().any(|label| label == "{A}"));
/// assert!(labels.iter().any(|label| label == "{B}"));
/// assert!(labels.iter().any(|label| label == "{A,B}"));
/// ```
pub trait TraceUnfold<Key>: Independence<HiddenData<Key, EdgeIndex>> + Sized
where
    Key: Eq + Hash + Clone + Ord,
{
    type EdgeData: Ord + Hash;
    type NodeData;
    type NodeStorage: NodeStorageOps<NodeData = Self::NodeData>;
    type HedgeData;

    fn graph(
        &self,
    ) -> &HedgeGraph<Self::EdgeData, Self::NodeData, Self::HedgeData, Self::NodeStorage>;

    fn key(&self, e: EdgeIndex) -> Key;

    /// Returns all operation sequences represented by one traversed edge.
    ///
    /// The default is a single realization containing the edge's own label unchanged. Override
    /// this when one graph transition may also be interpreted as a factorized sequence of more
    /// primitive operations.
    ///
    /// # Examples
    ///
    /// Here one edge labeled `AB` is unfolded both atomically and as the commuting pair
    /// `A, B`, producing two distinct trace nodes over the same underlying graph node.
    ///
    /// ```
    /// use linnet::half_edge::{
    ///     HedgeGraph, NoData, algorithms::trace_unfold::{HiddenData, Independence, TraceUnfold},
    ///     builder::HedgeGraphBuilder, involution::EdgeIndex, nodestore::DefaultNodeStore,
    /// };
    ///
    /// struct FactoringGraph {
    ///     graph: HedgeGraph<String, &'static str, NoData, DefaultNodeStore<&'static str>>,
    /// }
    ///
    /// impl Independence<HiddenData<String, EdgeIndex>> for FactoringGraph {
    ///     fn independent(
    ///         &self,
    ///         a: &HiddenData<String, EdgeIndex>,
    ///         b: &HiddenData<String, EdgeIndex>,
    ///     ) -> bool {
    ///         !a.order.chars().any(|c| b.order.contains(c))
    ///     }
    /// }
    ///
    /// impl TraceUnfold<String> for FactoringGraph {
    ///     type EdgeData = String;
    ///     type HedgeData = NoData;
    ///     type NodeData = &'static str;
    ///     type NodeStorage = DefaultNodeStore<&'static str>;
    ///
    ///     fn graph(
    ///         &self,
    ///     ) -> &HedgeGraph<Self::EdgeData, Self::NodeData, Self::HedgeData, Self::NodeStorage> {
    ///         &self.graph
    ///     }
    ///
    ///     fn key(&self, e: EdgeIndex) -> String {
    ///         self.graph[e].clone()
    ///     }
    ///
    ///     fn realizations_for_op(
    ///         &self,
    ///         op: HiddenData<String, EdgeIndex>,
    ///     ) -> Vec<Vec<HiddenData<String, EdgeIndex>>> {
    ///         if op.order != "AB" {
    ///             return vec![vec![op]];
    ///         }
    ///
    ///         vec![
    ///             vec![op.clone()],
    ///             vec![
    ///                 HiddenData {
    ///                     order: "A".to_string(),
    ///                     data: op.data,
    ///                 },
    ///                 HiddenData {
    ///                     order: "B".to_string(),
    ///                     data: op.data,
    ///                 },
    ///             ],
    ///         ]
    ///     }
    /// }
    ///
    /// let mut builder = HedgeGraphBuilder::<String, &'static str>::new();
    /// let empty = builder.add_node("empty");
    /// let ab = builder.add_node("AB");
    /// builder.add_edge(empty, ab, "AB".to_string(), true);
    ///
    /// let graph = FactoringGraph {
    ///     graph: builder.build::<DefaultNodeStore<&'static str>>(),
    /// };
    /// let unfolded = graph.trace_unfold::<DefaultNodeStore<usize>>(empty);
    /// let labels: Vec<_> = unfolded.iter_nodes().map(|(_, _, key)| key.to_string()).collect();
    ///
    /// assert_eq!(unfolded.n_nodes(), 3);
    /// assert!(labels.iter().any(|label| label == "∅"));
    /// assert!(labels.iter().any(|label| label == "{AB}"));
    /// assert!(labels.iter().any(|label| label == "{A,B}"));
    /// ```
    fn realizations_for_op(
        &self,
        op: HiddenData<Key, EdgeIndex>,
    ) -> Vec<Vec<HiddenData<Key, EdgeIndex>>> {
        vec![vec![op]]
    }

    #[allow(clippy::type_complexity)]
    fn trace_unfold_of<M: NodeStorageOps<NodeData = usize>, S: SubSetLike>(
        &self,
        subgraph: &S,
        start: NodeIndex,
    ) -> UnfoldedTraceGraph<Self, TraceKey<Key, EdgeIndex>, M::OpStorage<TraceKey<Key, EdgeIndex>>>
    {
        // `traces` is the canonical dedup table. Each unfolded node is identified by the pair
        // `(source_graph_node, canonical_trace_key)`, so paths that differ only by commuting
        // independent operations collapse to the same unfolded node.
        let root = (start, TraceKey::empty());
        let mut q = VecDeque::new();
        let mut traces: IndexSet<(NodeIndex, TraceKey<Key, EdgeIndex>)> = IndexSet::new();
        let mut builder: HedgeGraphBuilder<EdgeIndex, usize, NoData> = HedgeGraphBuilder::new();
        let (ind, _) = traces.insert_full(root);
        let nid = builder.add_node(ind);
        q.push_back((nid, start, TraceKey::empty()));
        let g = self.graph();

        while let Some((bnid, nid, key)) = q.pop_front() {
            for hedge in g.iter_crown_in(subgraph, nid) {
                if g.flow(hedge) == Flow::Source {
                    if let Some(to_node) = g.involved_node_id(hedge) {
                        let op = HiddenData {
                            order: self.key(g[&hedge]),
                            data: g[&hedge],
                        };
                        for realization in self.realizations_for_op(op) {
                            let new_key = realization
                                .into_iter()
                                .fold(key.clone(), |acc, factor| acc.push(self, factor));
                            let (ind, is_new) = traces.insert_full((to_node, new_key.clone()));
                            if is_new {
                                let bbnid = builder.add_node(ind);
                                debug_assert_eq!(bbnid.0, ind);
                                q.push_back((bbnid, to_node, new_key.clone()))
                            }

                            builder.add_edge(bnid, NodeIndex(ind), g[&hedge], true);
                        }
                    }
                }
            }
        }

        let source_nodes = traces.iter().map(|(source_node, _)| *source_node).collect();
        let unfolded = UnfoldedTraceGraph::new(
            builder.build::<M>().map(
                |_, _, v| {
                    let mut trace: TraceKey<Key, EdgeIndex> = TraceKey::empty();
                    let v = &mut traces.get_index_mut2(v).unwrap().1;
                    std::mem::swap(v, &mut trace);
                    trace
                },
                |_, _, _, _, a| a,
                |_, h| h,
            ),
            source_nodes,
        );
        unfolded
            .validate_invariant()
            .expect("trace-unfolded graph violates non-union unique-parent invariant");
        unfolded
    }

    #[allow(clippy::type_complexity)]
    fn trace_unfold<M: NodeStorageOps<NodeData = usize>>(
        &self,
        start: NodeIndex,
    ) -> UnfoldedTraceGraph<Self, TraceKey<Key, EdgeIndex>, M::OpStorage<TraceKey<Key, EdgeIndex>>>
    {
        self.trace_unfold_of::<M, _>(&self.graph().full_filter(), start)
    }
}

impl<G, V, N: NodeStorageOps<NodeData = V>> UnfoldedTraceGraph<G, V, N> {
    /// Returns all unfolded parents of `node` together with the source-graph edge realized by
    /// that branch.
    fn incoming_parents(&self, node: NodeIndex) -> Vec<(NodeIndex, EdgeIndex)> {
        self.graph
            .iter_crown(node)
            .filter(|hedge| self.graph.flow(*hedge) == Flow::Sink)
            .filter_map(|hedge| {
                let parent = self.graph.involved_node_id(hedge)?;
                let unfolded_edge = self.graph[&hedge];
                Some((parent, self.graph[unfolded_edge]))
            })
            .collect()
    }

    fn all_ops<'b, K>(
        &'b self,
        node: NodeIndex,
    ) -> impl Iterator<Item = &'b HiddenData<K, EdgeIndex>> + 'b
    where
        K: Eq + Hash + Clone + Ord + 'b,
        V: AsRef<TraceKey<K, EdgeIndex>>,
    {
        self.graph[node]
            .as_ref()
            .iter_levels_top_down()
            .flat_map(|level| level.iter_leaf_ops())
    }

    fn node_contains_op<K>(&self, node: NodeIndex, target: &HiddenData<K, EdgeIndex>) -> bool
    where
        K: Eq + Hash + Clone + Ord,
        V: AsRef<TraceKey<K, EdgeIndex>>,
    {
        self.all_ops(node).any(|op| op == target)
    }

    fn introduced_ops<K>(
        &self,
        parent: NodeIndex,
        child: NodeIndex,
    ) -> Vec<HiddenData<K, EdgeIndex>>
    where
        K: Eq + Hash + Clone + Ord,
        V: AsRef<TraceKey<K, EdgeIndex>>,
    {
        let mut parent_ops: Vec<_> = self.all_ops(parent).cloned().collect();
        let mut introduced = Vec::new();
        for op in self.all_ops(child) {
            // `child` is the canonicalized extension of `parent`, so the multiset difference
            // between their operations is exactly the set of newly introduced operations on that
            // step.
            if let Some(index) = parent_ops.iter().position(|parent_op| parent_op == op) {
                parent_ops.swap_remove(index);
            } else {
                introduced.push(op.clone());
            }
        }
        introduced
    }

    fn closest_common_ancestor(&self, nodes: &[NodeIndex]) -> Option<NodeIndex> {
        let first = *nodes.first()?;
        let first_chain: Vec<_> = std::iter::once(first)
            .chain(self.path_to_root(first).map(|(parent, _)| parent))
            .collect();

        // `path_to_root` follows the unique-parent spine to the empty trace, so intersecting
        // those spines is enough to recover the closest shared prefix of the given nodes.
        first_chain.into_iter().find(|candidate| {
            nodes.iter().skip(1).all(|node| {
                std::iter::once(*node)
                    .chain(self.path_to_root(*node).map(|(parent, _)| parent))
                    .any(|ancestor| ancestor == *candidate)
            })
        })
    }

    fn strip_independent_prefix<K>(
        &self,
        start: NodeIndex,
        target: &HiddenData<K, EdgeIndex>,
        indep: &impl Independence<HiddenData<K, EdgeIndex>>,
    ) -> NodeIndex
    where
        K: Eq + Hash + Clone + Ord,
        V: AsRef<TraceKey<K, EdgeIndex>>,
    {
        let mut current = start;

        loop {
            if self.is_disjoint_union(current) {
                let leaf_ops: Vec<_> = self.leaf_ops(current).cloned().collect();
                // If every operation introduced by this union is independent of `target`, the
                // relevant reusable prefix is shared by all incoming branches, so we can replace
                // the union by their closest common ancestor.
                if leaf_ops.iter().all(|op| indep.independent(op, target)) {
                    let parents: Vec<_> = self
                        .incoming_parents(current)
                        .into_iter()
                        .map(|(parent, _)| parent)
                        .collect();
                    if let Some(prefix) = self.closest_common_ancestor(&parents) {
                        current = prefix;
                        continue;
                    }
                }
                return current;
            }

            let Some((parent, _)) = self.unique_parent(current) else {
                return current;
            };
            let introduced = self.introduced_ops(parent, current);
            // As long as the whole step from `parent -> current` is independent of `target`, that
            // step can be factored out of the reusable prefix and we continue with `parent`.
            if introduced.iter().all(|op| indep.independent(op, target)) {
                current = parent;
            } else {
                return current;
            }
        }
    }

    pub fn leaf_ops<'b, K>(
        &'b self,
        node: NodeIndex,
    ) -> impl DoubleEndedIterator<Item = &'b HiddenData<K, EdgeIndex>> + 'b
    where
        K: Eq + Hash + Clone + Ord + 'b,
        V: AsRef<TraceKey<K, EdgeIndex>>,
    {
        self.graph[node].as_ref().iter_leaf_ops()
    }

    /// Returns whether `node` represents a disjoint-union step in the unfolded graph.
    pub fn is_disjoint_union(&self, node: NodeIndex) -> bool {
        self.incoming_parents(node).len() > 1
    }

    /// Returns the unique parent of `node` when the node is not a disjoint union.
    ///
    /// The accompanying [`EdgeIndex`] is the original source-graph edge carried by that branch.
    pub fn unique_parent(&self, node: NodeIndex) -> Option<(NodeIndex, EdgeIndex)> {
        let mut incoming = self.incoming_parents(node).into_iter();
        let parent = incoming.next()?;
        if incoming.next().is_some() {
            return None;
        }
        Some(parent)
    }

    /// Walks the unique-parent chain from `node` to the root.
    ///
    /// The iterator stops before a disjoint union, because such a node no longer has a single
    /// canonical parent chain.
    pub fn path_to_root(
        &self,
        node: NodeIndex,
    ) -> impl Iterator<Item = (NodeIndex, EdgeIndex)> + '_ {
        let mut current = Some(node);
        std::iter::from_fn(move || {
            let node = current?;
            let parent = self.unique_parent(node)?;
            current = Some(parent.0);
            Some(parent)
        })
    }

    /// Returns the closest reusable prefix node for the contribution of `target` at `node`.
    ///
    /// The traversal is specific to trace-unfolded graphs:
    ///
    /// - it follows the unique branch that still contains `target`
    /// - across a disjoint union, it continues with the branch containing `target`
    /// - once `target` disappears, it strips off any prefix whose newly introduced operations are
    ///   all independent of `target`
    ///
    /// This is the graph-level query used by clients that cache prefix computations per unfolded
    /// node and need the closest reusable prefix for one leaf operation.
    pub fn op_dependency_frontier<K>(
        &self,
        node: NodeIndex,
        target: &HiddenData<K, EdgeIndex>,
        indep: &impl Independence<HiddenData<K, EdgeIndex>>,
    ) -> Option<NodeIndex>
    where
        K: Eq + Hash + Clone + Ord,
        V: AsRef<TraceKey<K, EdgeIndex>>,
    {
        let mut current = node;
        loop {
            if self.is_disjoint_union(current) {
                if let Some((parent, _)) = self
                    .incoming_parents(current)
                    .into_iter()
                    .find(|(parent, _)| self.node_contains_op(*parent, target))
                {
                    // The target operation is still present on exactly one incoming branch, so
                    // that branch remains the relevant prefix chain.
                    current = parent;
                    continue;
                }
                // The target no longer appears on any incoming branch, so the reusable prefix is
                // the closest common ancestor of all incoming branches rather than any particular
                // parent.
                return self.closest_common_ancestor(
                    &self
                        .incoming_parents(current)
                        .into_iter()
                        .map(|(parent, _)| parent)
                        .collect::<Vec<_>>(),
                );
            }

            let Some((parent, _)) = self.unique_parent(current) else {
                return Some(current);
            };
            if self.node_contains_op(parent, target) {
                // We have not yet reached the step at which `target` was introduced.
                current = parent;
                continue;
            }

            return Some(self.strip_independent_prefix(parent, target, indep));
        }
    }

    /// Checks the structural invariant used by the traversal helpers above.
    ///
    /// Distinct incoming branches of the same disjoint union must carry distinct source-graph
    /// edges; otherwise helpers such as [`unique_parent`](Self::unique_parent) and
    /// [`op_dependency_frontier`](Self::op_dependency_frontier) become ambiguous.
    pub fn validate_invariant<K>(&self) -> Result<(), String>
    where
        K: Eq + Hash + Clone + Ord,
        V: AsRef<TraceKey<K, EdgeIndex>>,
    {
        for (node, _, _) in self.graph.iter_nodes() {
            let incoming = self.incoming_parents(node);
            for (i, (_, edge)) in incoming.iter().enumerate() {
                if incoming
                    .iter()
                    .skip(i + 1)
                    .any(|(_, other_edge)| other_edge == edge)
                {
                    return Err(format!(
                        "node {node} has duplicate incoming branches for edge {edge}"
                    ));
                }
            }
        }

        Ok(())
    }
}

impl<E: Ord + Hash, H, V, N: NodeStorageOps<NodeData = V>, K: Eq + Hash + Clone + Ord>
    TraceUnfold<K> for HedgeGraph<E, V, H, N>
where
    Self: Independence<HiddenData<K, EdgeIndex>>,
    Self: Key<K>,
{
    type EdgeData = E;
    type NodeData = V;
    type NodeStorage = N;
    type HedgeData = H;

    fn graph(
        &self,
    ) -> &HedgeGraph<Self::EdgeData, Self::NodeData, Self::HedgeData, Self::NodeStorage> {
        self
    }

    fn key(&self, e: EdgeIndex) -> K {
        Key::key(self, e)
    }
}

// impl<E, H, V> Independence<HedgePair> for HedgeGraph<E, V, H>
// where
//     Self: Independence<E> + Independence<V> + Independence<H>,
// {
//     fn independent(&self, a: &HedgePair, b: &HedgePair) -> bool {
//         match (a, b) {
//             (
//                 HedgePair::Paired { source, sink },
//                 HedgePair::Paired {
//                     source: so,
//                     sink: sk,
//                 },
//             ) => {

//                 let a = self[source];
//                 let v = self.node_id(*source);
//                 self.independent(&self[source], &self[so])
//                     && self.independent(&self[sink], &self[sk]) && self.independent(&self[self[]], b)
//             }
//             _ => {}
//         }

//         false
//     }
// }

#[derive(Clone, Debug)]
pub struct HiddenData<O, D> {
    /// the semantic label used for identity + ordering
    pub order: O,
    /// extra context used only for independence tests
    pub data: D,
}

// Equality/hash only by `order` so merging works.
impl<O: PartialEq, D> PartialEq for HiddenData<O, D> {
    fn eq(&self, other: &Self) -> bool {
        self.order == other.order
    }
}
impl<O: Eq, D> Eq for HiddenData<O, D> {}

impl<O: Hash, D> Hash for HiddenData<O, D> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.order.hash(state);
    }
}

// Sorting only by `order` so canonicalization works.
impl<O: Ord, D> Ord for HiddenData<O, D> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.order.cmp(&other.order)
    }
}
impl<O: PartialOrd, D> PartialOrd for HiddenData<O, D> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.order.partial_cmp(&other.order)
    }
}

/// A borrowed view over a trace key's Foata levels.
///
/// # Examples
///
/// ```
/// use linnet::half_edge::algorithms::trace_unfold::{HiddenData, TraceKey};
///
/// let trace = TraceKey::from_levels(vec![
///     vec![HiddenData { order: "A", data: () }],
///     vec![HiddenData { order: "B", data: () }],
/// ]);
/// let levels: Vec<_> = trace
///     .view()
///     .iter_levels_top_down()
///     .map(|level| level.to_string())
///     .collect();
///
/// assert_eq!(levels, vec!["{A}", "{B}"]);
/// ```
#[derive(Clone, Copy, Debug)]
pub struct TraceKeyView<'a, O, D> {
    levels: &'a [Vec<HiddenData<O, D>>],
}

impl<'a, O, D> TraceKeyView<'a, O, D> {
    /// Returns whether the viewed trace has no levels.
    pub fn is_empty(&self) -> bool {
        self.levels.is_empty()
    }

    /// Returns the total number of operations across all viewed levels.
    pub fn op_count(&self) -> usize {
        self.levels.iter().map(Vec::len).sum()
    }

    /// Iterates over the Foata levels in their stored order.
    ///
    /// Each item is itself a [`TraceKeyView`] over exactly one level.
    pub fn iter_levels_top_down(
        &self,
    ) -> impl DoubleEndedIterator<Item = TraceKeyView<'a, O, D>> + ExactSizeIterator + 'a {
        self.levels.iter().map(|level| TraceKeyView {
            levels: slice::from_ref(level),
        })
    }

    /// Iterates over the Foata levels in reverse order.
    ///
    /// Each item is itself a [`TraceKeyView`] over exactly one level.
    pub fn iter_levels_bottom_up(
        &self,
    ) -> impl DoubleEndedIterator<Item = TraceKeyView<'a, O, D>> + ExactSizeIterator + 'a {
        self.levels.iter().rev().map(|level| TraceKeyView {
            levels: slice::from_ref(level),
        })
    }

    /// Iterates over the operations in the leaf level of the viewed trace.
    pub fn iter_leaf_ops(&self) -> impl DoubleEndedIterator<Item = &'a HiddenData<O, D>> + 'a {
        self.levels
            .last()
            .into_iter()
            .flat_map(|level| level.iter())
    }

    /// Splits the viewed trace into its strict prefix and its leaf level.
    pub fn split_last_level(&self) -> Option<(TraceKeyView<'a, O, D>, TraceKeyView<'a, O, D>)> {
        let (leaf, prefix) = self.levels.split_last()?;
        Some((
            TraceKeyView { levels: prefix },
            TraceKeyView {
                levels: slice::from_ref(leaf),
            },
        ))
    }

    /// Writes the viewed trace using explicit Foata levels and a custom label mapping.
    pub fn write_foata_like<W: fmt::Write>(
        &self,
        f: &mut W,
        mut map: impl FnMut(&O) -> String,
    ) -> fmt::Result
    where
        O: Op,
    {
        if self.is_empty() {
            return write!(f, "∅");
        }

        for (i, level) in self.iter_levels_top_down().enumerate() {
            if i > 0 {
                write!(f, " · ")?;
            }
            write!(f, "{{")?;
            for (j, op) in level.iter_leaf_ops().enumerate() {
                if j > 0 {
                    write!(f, ",")?;
                }
                write!(f, "{}", map(&op.order))?;
            }
            write!(f, "}}")?;
        }
        Ok(())
    }
}

impl<'a, O: Clone, D: Clone> TraceKeyView<'a, O, D> {
    /// Clones the viewed levels into an owned [`TraceKey`].
    ///
    /// # Examples
    ///
    /// ```
    /// use linnet::half_edge::algorithms::trace_unfold::{HiddenData, TraceKey};
    ///
    /// let trace = TraceKey::from_levels(vec![
    ///     vec![HiddenData { order: "A", data: () }],
    ///     vec![HiddenData { order: "B", data: () }],
    /// ]);
    ///
    /// let cloned = trace.view().to_owned();
    ///
    /// assert_eq!(cloned, trace);
    /// ```
    pub fn to_owned(&self) -> TraceKey<O, D> {
        TraceKey::from_levels(self.levels.to_vec())
    }
}

impl<'a, O: Clone, D: Clone> From<TraceKeyView<'a, O, D>> for TraceKey<O, D> {
    fn from(value: TraceKeyView<'a, O, D>) -> Self {
        value.to_owned()
    }
}

impl<O, D> Display for TraceKeyView<'_, O, D>
where
    O: Op + Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.write_foata_like(f, |o| o.to_string())
    }
}

/// A canonical Foata-style trace normal form.
///
/// `TraceKey` stores a trace as a list of levels, where each inner vector contains operations
/// that are pairwise independent and therefore commute. The levels are ordered, so operations in
/// later levels could not be moved further left without breaking the independence relation used
/// when the key was constructed.
///
/// `TraceKey` is the owning form. Use [`TraceKey::view`] when you only need borrowed access to
/// the levels, and [`TraceKey::from_levels`] when you already have an explicit canonical level
/// structure and want to build an owned key directly.
///
/// The intended split is:
///
/// - [`TraceKey`] owns a canonical trace.
/// - [`TraceKeyView`] borrows the same trace without exposing the private storage.
///
/// # Examples
///
/// ```
/// use linnet::half_edge::algorithms::trace_unfold::{HiddenData, Independence, TraceKey};
///
/// struct CommutesOnlyAandC;
///
/// impl Independence<HiddenData<&'static str, ()>> for CommutesOnlyAandC {
///     fn independent(
///         &self,
///         a: &HiddenData<&'static str, ()>,
///         b: &HiddenData<&'static str, ()>,
///     ) -> bool {
///         matches!(
///             (a.order, b.order),
///             ("A", "C") | ("C", "A")
///         )
///     }
/// }
///
/// let trace = TraceKey::empty()
///     .push(&CommutesOnlyAandC, HiddenData { order: "A", data: () })
///     .push(&CommutesOnlyAandC, HiddenData { order: "B", data: () })
///     .push(&CommutesOnlyAandC, HiddenData { order: "C", data: () });
///
/// assert_eq!(trace.to_string(), "{A,C} · {B}");
/// assert_eq!(trace.view().to_owned(), trace);
/// ```
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct TraceKey<O, D> {
    levels: Vec<Vec<HiddenData<O, D>>>,
}

impl<O, D> AsRef<TraceKey<O, D>> for TraceKey<O, D> {
    fn as_ref(&self) -> &TraceKey<O, D> {
        self
    }
}

impl<O, D> TraceKey<O, D> {
    /// Returns the empty trace.
    ///
    /// # Examples
    ///
    /// ```
    /// use linnet::half_edge::algorithms::trace_unfold::TraceKey;
    ///
    /// let empty = TraceKey::<&'static str, ()>::empty();
    /// assert!(empty.is_empty());
    /// assert_eq!(empty.to_string(), "∅");
    /// ```
    pub fn empty() -> Self {
        Self { levels: Vec::new() }
    }

    /// Builds a trace key from explicit Foata levels.
    ///
    /// This does not validate that the levels satisfy any independence relation; it assumes the
    /// caller is providing an already-canonical arrangement.
    ///
    /// # Examples
    ///
    /// ```
    /// use linnet::half_edge::algorithms::trace_unfold::{HiddenData, TraceKey};
    ///
    /// let trace = TraceKey::from_levels(vec![
    ///     vec![HiddenData { order: "A", data: () }],
    ///     vec![
    ///         HiddenData { order: "B", data: () },
    ///         HiddenData { order: "C", data: () },
    ///     ],
    /// ]);
    ///
    /// assert_eq!(trace.to_string(), "{A} · {B,C}");
    /// ```
    pub fn from_levels(levels: Vec<Vec<HiddenData<O, D>>>) -> Self {
        Self { levels }
    }

    /// Borrows the trace as a lightweight view over its levels.
    ///
    /// # Examples
    ///
    /// ```
    /// use linnet::half_edge::algorithms::trace_unfold::{HiddenData, TraceKey};
    ///
    /// let trace = TraceKey::from_levels(vec![
    ///     vec![HiddenData { order: "A", data: () }],
    ///     vec![HiddenData { order: "B", data: () }],
    /// ]);
    ///
    /// assert!(!trace.view().is_empty());
    /// let levels: Vec<_> = trace
    ///     .view()
    ///     .iter_levels_bottom_up()
    ///     .map(|level| level.to_string())
    ///     .collect();
    /// assert_eq!(levels, vec!["{B}", "{A}"]);
    /// ```
    pub fn view(&self) -> TraceKeyView<'_, O, D> {
        TraceKeyView {
            levels: &self.levels,
        }
    }

    /// Returns whether the trace has no levels.
    ///
    /// # Examples
    ///
    /// ```
    /// use linnet::half_edge::algorithms::trace_unfold::{HiddenData, TraceKey};
    ///
    /// let empty = TraceKey::<&'static str, ()>::empty();
    /// let non_empty = TraceKey::from_levels(vec![vec![HiddenData {
    ///     order: "A",
    ///     data: (),
    /// }]]);
    ///
    /// assert!(empty.is_empty());
    /// assert!(!non_empty.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.levels.is_empty()
    }

    /// Returns the total number of operations across all levels.
    pub fn op_count(&self) -> usize {
        self.view().op_count()
    }

    /// Iterates over the Foata levels in their stored order.
    ///
    /// This is the same top-down order used by [`Display`] and [`TraceKey::write_foata_like`]:
    /// earlier levels are yielded first, later levels last.
    ///
    /// # Examples
    ///
    /// ```
    /// use linnet::half_edge::algorithms::trace_unfold::{HiddenData, TraceKey};
    ///
    /// let trace = TraceKey::from_levels(vec![
    ///         vec![HiddenData { order: "A", data: () }],
    ///         vec![
    ///             HiddenData { order: "B", data: () },
    ///             HiddenData { order: "C", data: () },
    ///         ],
    ///     ]);
    ///
    /// let levels: Vec<_> = trace
    ///     .iter_levels_top_down()
    ///     .map(|level| level.to_string())
    ///     .collect();
    ///
    /// assert_eq!(levels, vec!["{A}", "{B,C}"]);
    /// ```
    pub fn iter_levels_top_down(
        &self,
    ) -> impl DoubleEndedIterator<Item = TraceKeyView<'_, O, D>> + ExactSizeIterator + '_ {
        self.view().iter_levels_top_down()
    }

    /// Iterates over the Foata levels in reverse order.
    ///
    /// This yields the last level first and works entirely by reference.
    ///
    /// # Examples
    ///
    /// ```
    /// use linnet::half_edge::algorithms::trace_unfold::{HiddenData, TraceKey};
    ///
    /// let trace = TraceKey::from_levels(vec![
    ///         vec![HiddenData { order: "A", data: () }],
    ///         vec![
    ///             HiddenData { order: "B", data: () },
    ///             HiddenData { order: "C", data: () },
    ///         ],
    ///     ]);
    ///
    /// let levels: Vec<_> = trace
    ///     .iter_levels_bottom_up()
    ///     .map(|level| level.to_string())
    ///     .collect();
    ///
    /// assert_eq!(levels, vec!["{B,C}", "{A}"]);
    /// ```
    pub fn iter_levels_bottom_up(
        &self,
    ) -> impl DoubleEndedIterator<Item = TraceKeyView<'_, O, D>> + ExactSizeIterator + '_ {
        self.view().iter_levels_bottom_up()
    }

    /// Iterates over the operations in the leaf level of the trace.
    ///
    /// # Examples
    ///
    /// ```
    /// use linnet::half_edge::algorithms::trace_unfold::{HiddenData, TraceKey};
    ///
    /// let trace = TraceKey::from_levels(vec![
    ///     vec![HiddenData { order: "A", data: () }],
    ///     vec![
    ///         HiddenData { order: "B", data: () },
    ///         HiddenData { order: "C", data: () },
    ///     ],
    /// ]);
    ///
    /// let leaf: Vec<_> = trace.iter_leaf_ops().map(|op| op.order).collect();
    ///
    /// assert_eq!(leaf, vec!["B", "C"]);
    /// ```
    pub fn iter_leaf_ops(&self) -> impl DoubleEndedIterator<Item = &HiddenData<O, D>> + '_ {
        self.view().iter_leaf_ops()
    }

    /// Splits the trace into its strict prefix and its leaf level.
    pub fn split_last_level(&self) -> Option<(TraceKeyView<'_, O, D>, TraceKeyView<'_, O, D>)> {
        self.view().split_last_level()
    }

    /// Inserts one operation into the trace, placing it as late as possible.
    ///
    /// The new operation is appended to the latest level whose members are all independent of it.
    /// If no existing level is compatible, a new level is created at the end.
    ///
    /// # Examples
    ///
    /// ```
    /// use linnet::half_edge::algorithms::trace_unfold::{HiddenData, Independence, TraceKey};
    ///
    /// struct CommutesOnlyAandC;
    ///
    /// impl Independence<HiddenData<&'static str, ()>> for CommutesOnlyAandC {
    ///     fn independent(
    ///         &self,
    ///         a: &HiddenData<&'static str, ()>,
    ///         b: &HiddenData<&'static str, ()>,
    ///     ) -> bool {
    ///         matches!(
    ///             (a.order, b.order),
    ///             ("A", "C") | ("C", "A")
    ///         )
    ///     }
    /// }
    ///
    /// let trace = TraceKey::empty()
    ///     .push(&CommutesOnlyAandC, HiddenData { order: "A", data: () })
    ///     .push(&CommutesOnlyAandC, HiddenData { order: "B", data: () })
    ///     .push(&CommutesOnlyAandC, HiddenData { order: "C", data: () });
    ///
    /// assert_eq!(trace.to_string(), "{A,C} · {B}");
    /// ```
    pub fn push<I>(&self, indep: &I, op: HiddenData<O, D>) -> Self
    where
        I: Independence<HiddenData<O, D>>,
        O: Clone + Eq + Ord,
        D: Clone,
    {
        let mut levels = self.levels.clone();

        for i in (0..levels.len()).rev() {
            let ok = levels[i].iter().all(|b| indep.independent(&op, b));
            if ok {
                levels[i].push(op);
                levels[i].sort(); // uses O only
                return Self { levels };
            }
        }

        levels.push(vec![op]);
        Self { levels }
    }

    /// Writes the trace using explicit Foata levels and a custom label mapping.
    ///
    /// # Examples
    ///
    /// ```
    /// use linnet::half_edge::algorithms::trace_unfold::{HiddenData, TraceKey};
    ///
    /// let trace = TraceKey::from_levels(vec![
    ///         vec![HiddenData { order: 1_u8, data: () }],
    ///         vec![
    ///             HiddenData { order: 2_u8, data: () },
    ///             HiddenData { order: 3_u8, data: () },
    ///         ],
    ///     ]);
    ///
    /// let mut rendered = String::new();
    /// trace
    ///     .write_foata_like(&mut rendered, |order| format!("op{order}"))
    ///     .unwrap();
    ///
    /// assert_eq!(rendered, "{op1} · {op2,op3}");
    /// ```
    pub fn write_foata_like<W: fmt::Write>(
        &self,
        f: &mut W,
        mut map: impl FnMut(&O) -> String,
    ) -> fmt::Result
    where
        O: Op,
    {
        self.view().write_foata_like(f, &mut map)
    }
}

/// Display as Foata-like levels: {A,C} · {B} · {D,E}; empty = ∅
impl<O, D> Display for TraceKey<O, D>
where
    O: Op + Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.write_foata_like(f, |o| o.to_string())
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use crate::{
        dot,
        half_edge::{
            builder::HedgeGraphBuilder, involution::HedgePair, nodestore::DefaultNodeStore, NoData,
        },
        parser::{DotEdgeData, DotGraph, DotHedgeData, DotVertexData},
    };

    use super::*;
    use ahash::HashSet;

    fn has_overlap(a: &str, b: &str) -> bool {
        println!("Checking overlap between '{a}' and '{b}'");
        let set_a: BTreeSet<char> = a.chars().collect();
        b.chars().any(|c| set_a.contains(&c))
    }

    fn remove_chars(source: &str, remove: &str) -> String {
        let to_remove: HashSet<char> = remove.chars().collect();
        source.chars().filter(|c| !to_remove.contains(c)).collect()
    }

    struct FrontierGraph {
        graph: HedgeGraph<&'static str, &'static str, NoData, DefaultNodeStore<&'static str>>,
        extra_independent_pairs: &'static [(&'static str, &'static str)],
    }

    impl Key<&'static str> for FrontierGraph {
        fn key(&self, e: EdgeIndex) -> &'static str {
            self.graph[e]
        }
    }

    impl Independence<HiddenData<&'static str, EdgeIndex>> for FrontierGraph {
        fn independent(
            &self,
            a: &HiddenData<&'static str, EdgeIndex>,
            b: &HiddenData<&'static str, EdgeIndex>,
        ) -> bool {
            self.extra_independent_pairs.contains(&(a.order, b.order))
        }
    }

    impl TraceUnfold<&'static str> for FrontierGraph {
        type EdgeData = &'static str;
        type NodeData = &'static str;
        type NodeStorage = DefaultNodeStore<&'static str>;
        type HedgeData = NoData;

        fn graph(
            &self,
        ) -> &HedgeGraph<Self::EdgeData, Self::NodeData, Self::HedgeData, Self::NodeStorage>
        {
            &self.graph
        }

        fn key(&self, e: EdgeIndex) -> &'static str {
            self.graph[e]
        }
    }

    impl Key<String> for HedgeGraph<DotEdgeData, DotVertexData, DotHedgeData> {
        fn key(&self, e: EdgeIndex) -> String {
            let HedgePair::Paired { source, sink } = self[&e].1 else {
                return String::new();
            };

            let source_name = self[self.node_id(source)].name.clone().unwrap_or_default();
            let sink_name = self[self.node_id(sink)].name.clone().unwrap_or_default();

            // Remove "empty" from the name for key purposes
            remove_chars(&sink_name, &source_name)
        }
    }

    impl Independence<HiddenData<String, EdgeIndex>>
        for HedgeGraph<DotEdgeData, DotVertexData, DotHedgeData>
    {
        fn independent(
            &self,
            a: &HiddenData<String, EdgeIndex>,
            b: &HiddenData<String, EdgeIndex>,
        ) -> bool {
            let HedgePair::Paired {
                source: sourcea, ..
            } = self[&a.data].1
            else {
                return true;
            };

            let HedgePair::Paired {
                source: sourceb, ..
            } = self[&b.data].1
            else {
                return true;
            };

            let sourcan = self[self.node_id(sourcea)].name.clone().unwrap_or_default();
            let sourcbn = self[self.node_id(sourceb)].name.clone().unwrap_or_default();

            sourcan.as_str() == "empty"
                || sourcbn.as_str() == "empty"
                || !has_overlap(&sourcan, &sourcbn)
        }
    }

    // #[test]
    // fn overlappping_trace() {
    //     let indep = OverlapIndep;
    //     let t_ac = TraceKey::empty()
    //         .push(&indep, "A".to_string())
    //         .push(&indep, "B".to_string());
    //     let t_ca = TraceKey::empty()
    //         .push(&indep, "B".to_string())
    //         .push(&indep, "A".to_string());

    //     // Same trace key => same display
    //     assert_snapshot!(t_ac.to_string(), @"{A,B}");
    //     assert_snapshot!(t_ca.to_string(), @"{A,B}");
    // }
    #[test]
    fn test_tbt() {
        let graph: DotGraph = dot!(digraph{
            empty [id=0 name=""]
            empty -> AB [id=0 ];
            empty -> CD [id=1];
            AB -> ABCD [id=2]
            AB -> ABCEF [id=3 ]
            CD -> BCDEF [id=4 ]
            CD -> ABCD [id=5]
            ABCD -> ABCDEF [id=6]
            BCDEF -> ABCDEF [id=7]
            ABCEF -> ABCDEF [id=8]
        })
        .unwrap();

        let g = graph.graph.transitive_closure().unwrap();
        let mut output = String::new();
        g.dot_impl_fmt(
            &mut output,
            &g.full_filter(),
            "start=2;\n",
            &|_| None,
            &|_| None,
            &|v| Some(format!("label=\"{}\"", v.name.clone().unwrap_or_default())),
        )
        .unwrap();
        insta::assert_snapshot!(output, @r#"
        digraph {
          node	 [shape=circle,height=0.1,label=""];
          overlap = "scale";
          layout = "neato";
          start=2;
          1	 [label="AB"];
          2	 [label="ABCD"];
          4	 [label="ABCEF"];
          3	 [label="ABCDEF"];
          5	 [label="BCDEF"];
          6	 [label="CD"];
          0	 [label="empty"];
          0:14:s	-> 1:15:s	 [id=0  color="red:blue;0.5"];
          0:16:s	-> 6:17:s	 [id=1  color="red:blue;0.5"];
          1:0:s	-> 2:1:s	 [id=2  color="red:blue;0.5"];
          1:2:s	-> 4:3:s	 [id=3  color="red:blue;0.5"];
          6:12:s	-> 5:13:s	 [id=4  color="red:blue;0.5"];
          6:10:s	-> 2:11:s	 [id=5  color="red:blue;0.5"];
          2:4:s	-> 3:5:s	 [id=6  color="red:blue;0.5"];
          5:8:s	-> 3:9:s	 [id=7  color="red:blue;0.5"];
          4:6:s	-> 3:7:s	 [id=8  color="red:blue;0.5"];
          0:18:s	-> 4:19:s	 [id=9  color="red:blue;0.5"];
          0:20:s	-> 2:21:s	 [id=10  color="red:blue;0.5"];
          0:22:s	-> 5:23:s	 [id=11  color="red:blue;0.5"];
          0:24:s	-> 3:25:s	 [id=12  color="red:blue;0.5"];
          1:26:s	-> 3:27:s	 [id=13  color="red:blue;0.5"];
          6:28:s	-> 3:29:s	 [id=14  color="red:blue;0.5"];
        }
        "#);

        let g = g.trace_unfold::<DefaultNodeStore<usize>>(NodeIndex(0));
        let mut output = String::new();
        g.dot_impl_fmt(
            &mut output,
            &g.full_filter(),
            "start=2;\n",
            &|_| None,
            &|a| Some(format!("label=\"{a}\"")),
            &|v| Some(format!("label=\"{v}\"")),
        )
        .unwrap();

        insta::assert_snapshot!(output, @r#"
        digraph {
          node	 [shape=circle,height=0.1,label=""];
          overlap = "scale";
          layout = "neato";
          start=2;
          0	 [label="∅"];
          1	 [label="{AB}"];
          2	 [label="{CD}"];
          3	 [label="{ABCEF}"];
          4	 [label="{ABCD}"];
          5	 [label="{BCDEF}"];
          6	 [label="{ABCDEF}"];
          7	 [label="{AB,CD}"];
          8	 [label="{AB,CEF}"];
          9	 [label="{AB,CDEF}"];
          10	 [label="{BEF,CD}"];
          11	 [label="{ABEF,CD}"];
          12	 [label="{ABCEF,D}"];
          13	 [label="{ABCD,EF}"];
          14	 [label="{A,BCDEF}"];
          15	 [label="{AB,CD} · {EF}"];
          16	 [label="{AB,CEF} · {D}"];
          17	 [label="{BEF,CD} · {A}"];
          0:0:s	-> 1:1:s	 [id=0  color="red:blue;0.5" label="e0"];
          0:2:s	-> 2:3:s	 [id=1  color="red:blue;0.5" label="e1"];
          0:4:s	-> 3:5:s	 [id=2  color="red:blue;0.5" label="e9"];
          0:6:s	-> 4:7:s	 [id=3  color="red:blue;0.5" label="e10"];
          0:8:s	-> 5:9:s	 [id=4  color="red:blue;0.5" label="e11"];
          0:10:s	-> 6:11:s	 [id=5  color="red:blue;0.5" label="e12"];
          1:12:s	-> 7:13:s	 [id=6  color="red:blue;0.5" label="e2"];
          1:14:s	-> 8:15:s	 [id=7  color="red:blue;0.5" label="e3"];
          1:16:s	-> 9:17:s	 [id=8  color="red:blue;0.5" label="e13"];
          2:18:s	-> 7:19:s	 [id=9  color="red:blue;0.5" label="e5"];
          2:20:s	-> 10:21:s	 [id=10  color="red:blue;0.5" label="e4"];
          2:22:s	-> 11:23:s	 [id=11  color="red:blue;0.5" label="e14"];
          3:24:s	-> 12:25:s	 [id=12  color="red:blue;0.5" label="e8"];
          4:26:s	-> 13:27:s	 [id=13  color="red:blue;0.5" label="e6"];
          5:28:s	-> 14:29:s	 [id=14  color="red:blue;0.5" label="e7"];
          7:30:s	-> 15:31:s	 [id=15  color="red:blue;0.5" label="e6"];
          8:32:s	-> 16:33:s	 [id=16  color="red:blue;0.5" label="e8"];
          10:34:s	-> 17:35:s	 [id=17  color="red:blue;0.5" label="e7"];
        }
        "#);
    }

    #[test]
    fn dependency_frontier_walkthrough() {
        let mut builder = HedgeGraphBuilder::<&'static str, &'static str>::new();
        let empty = builder.add_node("empty");
        let a = builder.add_node("A");
        let b = builder.add_node("B");
        let ab = builder.add_node("AB");
        let bd = builder.add_node("BD");
        let abd = builder.add_node("ABD");

        builder.add_edge(empty, a, "A", true);
        builder.add_edge(empty, b, "B", true);
        builder.add_edge(a, ab, "B", true);
        builder.add_edge(b, ab, "A", true);
        builder.add_edge(b, bd, "D", true);
        builder.add_edge(ab, abd, "D", true);
        builder.add_edge(bd, abd, "A", true);

        let graph = FrontierGraph {
            graph: builder.build::<DefaultNodeStore<&'static str>>(),
            extra_independent_pairs: &[("A", "B"), ("B", "A"), ("A", "D"), ("D", "A")],
        };
        println!("{}", graph.graph.dot_display(&graph.graph.full_filter()));
        let unfolded = graph.trace_unfold::<DefaultNodeStore<usize>>(empty);
        println!(
            "{}",
            unfolded.graph.dot_label(&unfolded.graph.full_filter())
        );

        let mut lines = unfolded
            .iter_nodes()
            .map(|(node, _, key)| {
                let mut entry = key.to_string();
                let frontiers = key
                    .view()
                    .iter_leaf_ops()
                    .map(|op| {
                        let frontier = unfolded
                            .op_dependency_frontier(node, op, &graph)
                            .expect("frontier must exist");
                        format!("{} -> {}", op.order, unfolded[frontier])
                    })
                    .collect::<Vec<_>>();
                // frontiers.sort();
                if !frontiers.is_empty() {
                    entry.push_str(" | ");
                    entry.push_str(&frontiers.join(", "));
                }
                entry
            })
            .collect::<Vec<_>>();
        lines.sort();

        insta::assert_snapshot!(lines.join("\n"), @r#"
        {A,B} | A -> ∅, B -> ∅
        {A,B} · {D} | D -> {A,B}
        {A} | A -> ∅
        {B} | B -> ∅
        {B} · {A,D} | A -> ∅, D -> {B}
        {B} · {D} | D -> {B}
        ∅
        "#);
    }

    #[test]
    fn multi_parent_node_need_not_have_matching_leaf_count() {
        let mut builder = HedgeGraphBuilder::<&'static str, &'static str>::new();
        let empty = builder.add_node("empty");
        let a = builder.add_node("A");
        let b = builder.add_node("B");
        let ab = builder.add_node("AB");
        let bd = builder.add_node("BD");
        let abd = builder.add_node("ABD");

        builder.add_edge(empty, a, "A", true);
        builder.add_edge(empty, b, "B", true);
        builder.add_edge(a, ab, "B", true);
        builder.add_edge(b, ab, "A", true);
        builder.add_edge(b, bd, "D", true);
        builder.add_edge(ab, abd, "D", true);
        builder.add_edge(bd, abd, "A", true);

        let graph = FrontierGraph {
            graph: builder.build::<DefaultNodeStore<&'static str>>(),
            extra_independent_pairs: &[("A", "B"), ("B", "A")],
        };
        println!("{}", graph.graph.dot_label(&graph.graph.full_filter()));
        let unfolded = graph.trace_unfold::<DefaultNodeStore<usize>>(empty);
        println!(
            "{}",
            unfolded.graph.dot_label(&unfolded.graph.full_filter())
        );
        let node = unfolded
            .iter_nodes()
            .find_map(|(node, _, key)| (key.to_string() == "{A,B} · {D}").then_some(node))
            .expect("expected canonical node {A,B} · {D}");

        let parent_labels = unfolded
            .incoming_parents(node)
            .into_iter()
            .map(|(parent, _)| unfolded[parent].to_string())
            .collect::<BTreeSet<_>>();
        let leaf_labels = unfolded
            .leaf_ops(node)
            .map(|op| op.order)
            .collect::<Vec<_>>();

        insta::assert_snapshot!(
            format!(
                "node={}\nparents={:?}\nleaf_ops={:?}",
                unfolded[node], parent_labels, leaf_labels
            ),
            @r#"
        node={A,B} · {D}
        parents={"{A,B}", "{B} · {D}"}
        leaf_ops=["D"]
        "#
        );

        assert_eq!(unfolded.incoming_parents(node).len(), 2);
        assert_eq!(unfolded.leaf_ops(node).count(), 1);
    }

    // #[test]
    // fn snapshot_empty() {
    //     let t = TraceKey::<MyOp>::empty();
    //     assert_snapshot!(t.to_string(), @"∅");
    // }

    // #[test]
    // fn snapshot_commutation_merges_ac_and_ca() {
    //     let indep = IndepACOnly;
    //     let t_ac = TraceKey::empty()
    //         .push(&indep, MyOp::A)
    //         .push(&indep, MyOp::C);
    //     let t_ca = TraceKey::empty()
    //         .push(&indep, MyOp::C)
    //         .push(&indep, MyOp::A);

    //     // Same trace key => same display
    //     assert_snapshot!(t_ac.to_string(), @"{A,C}");
    //     assert_snapshot!(t_ca.to_string(), @"{A,C}");
    // }

    // #[test]
    // fn snapshot_dependent_orders_differ() {
    //     let indep = IndepACOnly;

    //     let t_ab = TraceKey::empty()
    //         .push(&indep, MyOp::A)
    //         .push(&indep, MyOp::B);
    //     let t_ba = TraceKey::empty()
    //         .push(&indep, MyOp::B)
    //         .push(&indep, MyOp::A);

    //     assert_snapshot!(t_ab.to_string(), @"{A} · {B}");
    //     assert_snapshot!(t_ba.to_string(), @"{B} · {A}");
    // }

    // #[test]
    // fn snapshot_ac_then_b_creates_new_level() {
    //     let indep = IndepACOnly;

    //     let t = TraceKey::empty()
    //         .push(&indep, MyOp::A)
    //         .push(&indep, MyOp::C)
    //         .push(&indep, MyOp::B);

    //     assert_snapshot!(t.to_string(), @"{A,C} · {B}");
    // }

    // #[test]
    // fn snapshot_place_as_late_as_possible() {
    //     let indep = IndepACOnly;

    //     // Start with {B} · {C}; then push A.
    //     // A can't join {B} but can join {C} => {B} · {A,C}
    //     let t = TraceKey::empty()
    //         .push(&indep, MyOp::B)
    //         .push(&indep, MyOp::C)
    //         .push(&indep, MyOp::A);

    //     assert_snapshot!(t.to_string(), @"{B} · {A,C}");
    // }
}
