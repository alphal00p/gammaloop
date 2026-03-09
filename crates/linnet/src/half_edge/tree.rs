use std::{cmp::Ordering, collections::VecDeque};

use itertools::Itertools;

use crate::{
    half_edge::subgraph::{SuBitGraph, SubGraphLike},
    tree::{
        parent_pointer::ParentPointerStore, Forest, ForestNodeStore, ForestNodeStoreBfs,
        ForestNodeStoreDown, ForestNodeStorePreorder, RootId,
    },
};

use super::{
    involution::{Hedge, Involution},
    subgraph::{Cycle, Inclusion, InternalSubGraph, ModifySubSet, SubSetLike, SubSetOps},
    HedgeGraph, HedgeGraphError, NodeIndex, NodeStorageOps,
};

// pub struct HedgeTree<V, P: ForestNodeStore<NodeData = ()>, E> {
//     graph: HedgeGraph<E, V, Forest<V, P>>,
//     tree_subgraph: InternalSubGraph,
//     covers: SuBitGraph,
// }

// impl<V, P: ForestNodeStore<NodeData = ()>, E> HedgeTree<V, P, E> {
//     pub fn ancestor_nodes(&self,) {
//         self.
//     }
// }

// #[derive(Debug, Clone)]
// pub struct TraversalTreeRef<
//     'a,
//     E,
//     V,
//     N: NodeStorage<NodeData = V>,
//     P: ForestNodeStore<NodeData = ()>,
// > {
//     graph: &'a HedgeGraph<E, V, N>,
//     simple: SimpleTraversalTree<P>,
//     tree_subgraph: InternalSubGraph,
// }

#[derive(Debug, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
/// Represents the role of a graph node within a [`SimpleTraversalTree`].
///
/// When a traversal (like DFS or BFS) is performed on a graph, each node
/// of the graph can be classified based on its position and discovery time
/// within that traversal. This enum is typically used as the root data in the
/// `Forest` underlying a `SimpleTraversalTree`.
pub enum TTRoot {
    /// The node is a child of another node in the traversal tree.
    /// The `usize` value might represent discovery order or depth, depending on context.
    Child(usize),
    /// The node is the root of this specific traversal tree.
    Root,
    /// The node is not part of this traversal tree (e.g., it was not visited).
    None,
}

impl TTRoot {
    pub fn includes(&self) -> bool {
        match self {
            TTRoot::Child(_) => true,
            TTRoot::Root => true,
            TTRoot::None => false,
        }
    }

    pub fn is_root(&self) -> bool {
        match self {
            TTRoot::Child(_) => false,
            TTRoot::Root => true,
            TTRoot::None => false,
        }
    }

    pub fn is_child(&self) -> bool {
        match self {
            TTRoot::Child(_) => true,
            TTRoot::Root => false,
            TTRoot::None => false,
        }
    }
}

#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
// #[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))] // Manual Bincode implemented
/// Represents a traversal tree (e.g., DFS or BFS tree) derived from a [`HedgeGraph`].
///
/// This structure uses a generic [`Forest`] to store the tree topology. Each "root" in the
/// `Forest` corresponds to a node in the original `HedgeGraph`, and its data is a [`TTRoot`]
/// enum indicating its role (Root, Child, or None) in this specific traversal. The nodes
/// within each tree of the `Forest` typically represent the half-edges that form the
/// traversal path.
///
/// # Type Parameters
///
/// - `P`: The type of [`ForestNodeStore`] used by the underlying `Forest`.
///   It stores `()` as `NodeData` because the tree nodes (representing half-edges)
///   don't need additional data beyond their structural role in the traversal tree.
///   Defaults to [`ParentPointerStore<()>`].
pub struct SimpleTraversalTree<P: ForestNodeStore<NodeData = ()> = ParentPointerStore<()>> {
    /// The underlying forest structure representing the traversal tree(s).
    /// Each root in this forest is a graph node, and its data is its [`TTRoot`] status.
    /// Nodes within each tree of this forest represent half-edges.
    forest: Forest<TTRoot, P>,
    /// A bitmask representing the set of half-edges that constitute the actual
    /// edges of this traversal tree.
    // #[cfg_attr(feature = "bincode", bincode(with_serde))] // Handled by manual impl
    pub tree_subgraph: SuBitGraph,
}
#[cfg(feature = "bincode")]
impl<P: ForestNodeStore<NodeData = ()>> ::bincode::Encode for SimpleTraversalTree<P>
where
    P: ::bincode::Encode,
{
    fn encode<__E: ::bincode::enc::Encoder>(
        &self,
        encoder: &mut __E,
    ) -> core::result::Result<(), ::bincode::error::EncodeError> {
        ::bincode::Encode::encode(&self.forest, encoder)?;
        ::bincode::Encode::encode(&::bincode::serde::Compat(&self.tree_subgraph), encoder)?;
        core::result::Result::Ok(())
    }
}

#[cfg(feature = "bincode")]
impl<P: ForestNodeStore<NodeData = ()>, __Context> ::bincode::Decode<__Context>
    for SimpleTraversalTree<P>
where
    P: ::bincode::Decode<__Context>,
{
    fn decode<__D: ::bincode::de::Decoder<Context = __Context>>(
        decoder: &mut __D,
    ) -> core::result::Result<Self, ::bincode::error::DecodeError> {
        core::result::Result::Ok(Self {
            forest: ::bincode::Decode::decode(decoder)?,
            tree_subgraph: (<::bincode::serde::Compat<_> as ::bincode::Decode<__Context>>::decode(
                decoder,
            )?)
            .0,
        })
    }
}

#[cfg(feature = "bincode")]
impl<'__de, P: ForestNodeStore<NodeData = ()>, __Context> ::bincode::BorrowDecode<'__de, __Context>
    for SimpleTraversalTree<P>
where
    P: ::bincode::de::BorrowDecode<'__de, __Context>,
{
    fn borrow_decode<__D: ::bincode::de::BorrowDecoder<'__de, Context = __Context>>(
        decoder: &mut __D,
    ) -> core::result::Result<Self, ::bincode::error::DecodeError> {
        core::result::Result::Ok(Self {
            forest: ::bincode::BorrowDecode::<'_, __Context>::borrow_decode(decoder)?,
            tree_subgraph: (<::bincode::serde::BorrowCompat<_> as ::bincode::BorrowDecode<
                '_,
                __Context,
            >>::borrow_decode(decoder)?)
            .0,
        })
    }
}
impl<P: ForestNodeStore<NodeData = ()>> SimpleTraversalTree<P> {
    pub fn cast<P2: ForestNodeStore<NodeData = ()>>(self) -> SimpleTraversalTree<P2>
    where
        Forest<TTRoot, P2>: From<Forest<TTRoot, P>>,
    {
        SimpleTraversalTree {
            forest: self.forest.into(),
            tree_subgraph: self.tree_subgraph,
        }
    }
}

impl<P: ForestNodeStore<NodeData = ()>> SimpleTraversalTree<P> {
    pub fn node_order(&self) -> Vec<NodeIndex> {
        self.forest
            .iter_roots()
            .enumerate()
            .filter_map(|(i, (a, _))| match a {
                TTRoot::Root => Some((0, i)),
                TTRoot::Child(a) => Some((*a, i)),
                _ => None,
            })
            .sorted_by(|a, b| a.0.cmp(&b.0))
            .map(|a| NodeIndex(a.1))
            .collect()
    }
    pub fn covers<S: SubSetLike>(&self, subgraph: &S) -> SuBitGraph {
        // println!("calculating covers..");
        let mut covers = SuBitGraph::empty(self.tree_subgraph.size());

        // self.tree_subgraph.covers(graph)

        for i in 0..self.tree_subgraph.size() {
            if self.node_data(self.node_id(Hedge(i))).includes() && subgraph.includes(&Hedge(i)) {
                covers.add(Hedge(i));
            }
        }

        covers
    }

    // fn root_hedge(&self, node: NodeIndex) -> Hedge {
    //     self.forest[&RootId(node.0)].into()
    // }

    pub fn node_data(&self, node: NodeIndex) -> &TTRoot {
        &self.forest[RootId(node.0)]
    }

    pub fn internal<I: AsRef<Involution>>(&self, hedge: Hedge, _inv: I) -> bool {
        self.tree_subgraph.includes(&hedge)
    }

    pub fn tree_subgraph<I: AsRef<Involution>>(&self, inv: I) -> InternalSubGraph {
        let mut tree = SuBitGraph::empty(self.tree_subgraph.size());
        let inv = inv.as_ref();

        for (r, h) in self.forest.iter_roots() {
            let h: Hedge = (*h).into();
            if r.is_child() {
                tree.add(h);
                tree.add(inv.inv(h))
            }
        }

        unsafe { InternalSubGraph::new_unchecked(tree) }
    }

    pub fn get_cycle<I: AsRef<Involution>>(&self, cut: Hedge, inv: &I) -> Option<Cycle> {
        if self.internal(cut, inv) {
            return None;
        }
        let mut cycle = self.path_to_root(cut, inv);
        cycle.sym_diff_with(&self.path_to_root(inv.as_ref().inv(cut), inv));
        let mut cycle = Cycle::new_unchecked(cycle);
        cycle.loop_count = Some(1);

        Some(cycle)
    }

    pub fn node_id(&self, hedge: Hedge) -> NodeIndex {
        NodeIndex::from(self.forest.root(hedge.into()))
    }

    /// get the Swapping of the current node
    pub fn node_parent<I: AsRef<Involution>>(&self, from: NodeIndex, inv: I) -> Option<NodeIndex> {
        let root = RootId::from(from);

        if self.forest[root].is_child() {
            // if the current node is in the forest
            // get the involved hedge connected to the root pointing hedge of the current node
            let involved = inv.as_ref().inv(self.forest[&root].into());
            Some(self.node_id(involved))
        } else {
            None
        }
    }

    fn path_to_root<I: AsRef<Involution>>(&self, start: Hedge, inv: I) -> SuBitGraph {
        let mut path = SuBitGraph::empty(self.tree_subgraph.size());

        self.ancestor_iter_hedge(start, inv.as_ref())
            .for_each(|a| path.add(a));
        path
    }

    /// among the hedges of the same node get the one pointing to the root node.
    /// if the from hedge is that pointing to the root node, go to the involved hedge.
    pub fn hedge_parent<I: AsRef<Involution>>(&self, from: Hedge, inv: I) -> Option<Hedge> {
        let root = self.forest.root(from.into()); //Get "NodeId/RootId"

        match self.forest[root] {
            TTRoot::Child(_) => {
                let roothedge = self.forest[&root].into(); //Get "chosen" root among node hairs

                if from == roothedge {
                    //if it is the same as the input, go to the involved hedge
                    Some(inv.as_ref().inv(from))
                } else {
                    // else go to the "chosen" hedge
                    Some(roothedge)
                }
            }
            TTRoot::None => {
                // println!("None");
                None
            }
            TTRoot::Root => {
                // println!("Root");
                None
            } // if it is attached to the root node, it has no parent
        }
    }
}

/// An iterator that traverses upwards from a given starting half-edge to the root
/// of its traversal tree, yielding each [`Hedge`] along the path.
///
/// # Type Parameters
///
/// - `'a`: The lifetime of the borrowed [`SimpleTraversalTree`] and [`Involution`].
/// - `P`: The type of [`ForestNodeStore`] used by the `SimpleTraversalTree`.
pub struct TraversalTreeAncestorHedgeIterator<'a, P: ForestNodeStore<NodeData = ()>> {
    /// A reference to the traversal tree being traversed.
    tt: &'a SimpleTraversalTree<P>,
    /// A reference to the graph's involution, needed to find opposite half-edges.
    inv: &'a Involution,
    /// The current [`Hedge`] in the traversal, or `None` if iteration is finished.
    current: Option<Hedge>,
}

impl<P: ForestNodeStore<NodeData = ()>> Iterator for TraversalTreeAncestorHedgeIterator<'_, P> {
    type Item = Hedge;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(c) = self.current {
            self.current = self.tt.hedge_parent(c, self.inv);
            Some(c)
        } else {
            None
        }
    }
}

#[derive(Clone)]
/// An iterator that performs a pre-order traversal of the nodes in a [`SimpleTraversalTree`].
///
/// Starting from a given graph node, it visits the node itself and then recursively
/// visits its children in the traversal tree.
///
/// # Type Parameters
///
/// - `'a`: The lifetime of the borrowed [`SimpleTraversalTree`] and [`Involution`].
/// - `S`: The type of [`ForestNodeStore`] used by the `SimpleTraversalTree`, which
///   must also implement [`ForestNodeStoreDown`] to allow child iteration.
pub struct PreorderTraversalIter<'a, S: ForestNodeStore<NodeData = ()>>
where
    S: ForestNodeStoreBfs,
{
    /// A reference to the traversal tree being traversed.
    store: &'a SimpleTraversalTree<S>,
    /// A reference to the graph's involution, used for navigating child relationships.
    inv: &'a Involution,
    /// A stack to keep track of nodes to visit during the pre-order traversal.
    stack: Vec<NodeIndex>,
}

impl<'a, S: ForestNodeStoreBfs + ForestNodeStore<NodeData = ()>> PreorderTraversalIter<'a, S> {
    /// Create a new pre-order iterator starting at `start`.
    pub fn new(
        store: &'a SimpleTraversalTree<S>,
        inv: &'a impl AsRef<Involution>,
        start: NodeIndex,
    ) -> Self {
        PreorderTraversalIter {
            inv: inv.as_ref(),
            store,
            stack: vec![start],
        }
    }
}

impl<S: ForestNodeStoreBfs + ForestNodeStore<NodeData = ()>> Iterator
    for PreorderTraversalIter<'_, S>
{
    type Item = NodeIndex;
    fn next(&mut self) -> Option<Self::Item> {
        // Pop the next node from the stack
        let node = self.stack.pop()?;

        // Push children onto the stack in reverse order so the first child is processed next
        // (This assumes iter_children returns them in the desired forward order)
        for child in self
            .store
            .iter_children(node, &self.inv)
            .collect::<Vec<_>>()
            .into_iter()
            .rev()
        {
            self.stack.push(child);
        }

        Some(node)
    }
}

impl<P: ForestNodeStore<NodeData = ()>> SimpleTraversalTree<P> {
    pub fn debug_draw<FR>(
        &self,
        format_root: FR, // Closure to format root data (&R) -> String
    ) -> String
    where
        FR: FnMut(&TTRoot) -> String,

        P: ForestNodeStoreDown,
    {
        self.forest.debug_draw(format_root, |_| "".into())
    }
    pub fn iter_children<'a, I: AsRef<Involution>>(
        &'a self,
        node_id: NodeIndex,
        inv: &'a I,
    ) -> impl Iterator<Item = NodeIndex> + 'a
    where
        P: ForestNodeStoreBfs,
    {
        let root_node = self.forest[&RootId::from(node_id)];

        let root_node_data = self.forest[RootId::from(node_id)];

        // println!("{:?}", root_node);
        self.forest.iter_bfs(root_node).filter_map(move |a| {
            if root_node == a && !root_node_data.is_root() {
                return None;
            }
            let h = Hedge::from(a);
            // println!("Hedge:{h}");
            if self.tree_subgraph.includes(&h) {
                let invh = inv.as_ref().inv(a.into());
                // println!("Included");
                if invh != h {
                    // println!("Hi");
                    if self.tree_subgraph.includes(&invh) {
                        let child = self.node_id(invh);
                        // println!("Node:{child}");

                        let root = RootId::from(child);

                        if self.forest[root].is_child() {
                            Some(child)
                        } else {
                            None
                        }
                    } else {
                        None
                    }
                } else {
                    None
                }
            } else {
                None
            }
        })
    }

    // pub fn iter_descendents<'a, I: AsRef<Involution>>(
    //     &'a self,
    //     node_id: NodeIndex,
    //     inv: &'a I,
    // ) -> impl Iterator<Item = NodeIndex> + 'a
    // where
    //     P: ForestNodeStoreDown,
    // {
    //     self.iter_children(node_id, inv)
    //         .map(|n| self.iter_descendents(n, inv))
    //         .flatten()
    // }
    /// Returns a pre-order iterator over the hedges *within* the traversal tree.
    /// Requires the underlying Forest store `P` to support pre-order traversal.
    pub fn iter_preorder_tree_nodes<'a>(
        &'a self,
        inv: &'a impl AsRef<Involution>,
        start: NodeIndex,
    ) -> PreorderTraversalIter<'a, P>
    where
        P: ForestNodeStorePreorder<NodeData = ()> + ForestNodeStoreBfs,
    {
        PreorderTraversalIter::new(self, inv, start)
    }

    pub fn ancestor_iter_hedge<'a>(
        &'a self,
        start: Hedge,
        inv: &'a Involution,
    ) -> impl Iterator<Item = Hedge> + 'a {
        TraversalTreeAncestorHedgeIterator {
            tt: self,
            inv,
            current: Some(start),
        }
    }

    pub fn dot<E, V, H, N: NodeStorageOps<NodeData = V>>(
        &self,
        graph: &HedgeGraph<E, V, H, N>,
    ) -> String {
        let mut structure = graph.just_structure();
        structure.align_superficial_to_tree(self);

        structure.base_dot()
    }

    pub fn tree_order<I: AsRef<Involution>>(
        &self,
        left: NodeIndex,
        right: NodeIndex,
        inv: I,
    ) -> Option<Ordering> {
        if left == right {
            return Some(Ordering::Equal);
        }

        let left_to_root_passes_through_right = self
            .ancestor_iter_node(left, inv.as_ref())
            .any(|a| a == right);
        if left_to_root_passes_through_right {
            return Some(Ordering::Less);
        }

        let right_to_root_passes_through_left = self
            .ancestor_iter_node(right, inv.as_ref())
            .any(|a| a == left);
        if right_to_root_passes_through_left {
            return Some(Ordering::Greater);
        }

        None
    }

    /// Iterate over all half-edges in the tree.
    ///
    /// Each iteration yields a three-tuple:
    ///
    /// - the half-edge
    /// - whether it is part of a child, or root node, or outside of the tree
    /// - the root-pointing half-edge if it is actually part of the tree
    ///
    pub fn iter_hedges(&self) -> impl Iterator<Item = (Hedge, TTRoot, Option<Hedge>)> + '_ {
        self.forest.iter_node_ids().map(|a| {
            let root = self.forest.root(a); //Get "NodeId/RootId"
            let hedge: Hedge = a.into();
            let ttype = self.forest[root];
            let root_among_hedges = if self.tree_subgraph.includes(&hedge) {
                Some(self.forest[&root].into())
            } else {
                None
            };

            (hedge, ttype, root_among_hedges)
        })
    }
}

pub struct TraversalTreeAncestorNodeIterator<'a, P: ForestNodeStore<NodeData = ()>> {
    tt: &'a SimpleTraversalTree<P>,
    inv: &'a Involution,
    current: Option<NodeIndex>,
}

impl<P: ForestNodeStore<NodeData = ()>> Iterator for TraversalTreeAncestorNodeIterator<'_, P> {
    type Item = NodeIndex;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(c) = self.current {
            self.current = self.tt.node_parent(c, self.inv);
            Some(c)
        } else {
            None
        }
    }
}

impl<P: ForestNodeStore<NodeData = ()>> SimpleTraversalTree<P> {
    pub fn ancestor_iter_node<'a>(
        &'a self,
        from: NodeIndex,
        inv: &'a Involution,
    ) -> impl Iterator<Item = NodeIndex> + 'a {
        TraversalTreeAncestorNodeIterator {
            tt: self,
            inv,
            current: Some(from),
        }
    }
}

impl SimpleTraversalTree {
    pub fn empty<E, V, H, N: NodeStorageOps<NodeData = V>>(graph: &HedgeGraph<E, V, H, N>) -> Self {
        let forest = graph.node_store.to_forest(|_| TTRoot::None);
        SimpleTraversalTree {
            forest,
            tree_subgraph: graph.empty_subgraph(),
        }
    }

    // pub fn root(&self,Hedge) -> Hedge {
    //     self.forest.root()
    // }

    pub fn depth_first_traverse<S: SubGraphLike, E, V, H, N: NodeStorageOps<NodeData = V>>(
        graph: &HedgeGraph<E, V, H, N>,
        subgraph: &S,
        root_node: &NodeIndex,
        include_hedge: Option<Hedge>,
    ) -> Result<Self, HedgeGraphError> {
        let mut seen = subgraph.hairs(graph.iter_crown(*root_node));

        if seen.is_empty() {
            // if the root node is not in the subgraph
            return Err(HedgeGraphError::RootNodeNotInSubgraph(*root_node));
        }

        let mut stack = seen.included_iter().collect::<Vec<_>>();
        if let Some(r) = include_hedge {
            if graph.inv(r) == r {
                return Err(HedgeGraphError::ExternalHedgeIncluded(r)); //cannot include an external hedge in the traversal
            }

            let pos = stack.iter().find_position(|a| **a == r).map(|a| a.0);

            let last = stack.len() - 1;
            if let Some(pos) = pos {
                stack.swap(pos, last);
            } else {
                return Err(HedgeGraphError::NotInNode(r, *root_node));
            }
        }

        let mut init = Self::empty(graph);
        init.forest[RootId(root_node.0)] = TTRoot::Root;

        let mut iter_order = 0;

        while let Some(hedge) = stack.pop() {
            // if the hedge is not external get the neighbors of the paired hedge
            if let Some(cn) = graph.involved_node_crown(hedge) {
                let connected = graph.inv(hedge);

                if !seen.includes(&connected) && subgraph.includes(&connected) {
                    // if this new hedge hasn't been seen before, it means the node it belongs to
                    // is a new node in the traversal
                    init.tree_subgraph.add(connected);
                    init.tree_subgraph.add(hedge);
                    let node_id = init.forest.change_to_root(connected.into());
                    iter_order += 1;
                    init.forest[node_id] = TTRoot::Child(iter_order);
                } else {
                    continue;
                }

                // mark the new node as seen
                // seen.union_with(&cn.hairs);

                for i in cn {
                    if subgraph.includes(&i) {
                        seen.add(i);
                        if !seen.includes(&graph.inv(i)) {
                            stack.push(i);
                        }
                    }
                }
            }
        }

        // init.tree_subgraph = tree_subgraph;
        Ok(init)
    }

    pub fn breadth_first_traverse<S: SubGraphLike, E, V, H, N: NodeStorageOps<NodeData = V>>(
        graph: &HedgeGraph<E, V, H, N>,
        subgraph: &S,
        root_node: &NodeIndex,
        include_hedge: Option<Hedge>,
    ) -> Result<Self, HedgeGraphError> {
        let mut seen = subgraph.hairs(graph.iter_crown(*root_node));

        if seen.is_empty() {
            // if the root node is not in the subgraph
            return Err(HedgeGraphError::InvalidNode(*root_node));
        }

        let mut queue = seen.included_iter().collect::<VecDeque<_>>();
        if let Some(r) = include_hedge {
            if graph.inv(r) == r {
                return Err(HedgeGraphError::InvalidHedge(r)); //cannot include an external hedge in the traversal
            }
            let pos = queue.iter().find_position(|a| **a == r).map(|a| a.0);
            if let Some(pos) = pos {
                queue.swap(pos, 0);
            } else {
                return Err(HedgeGraphError::InvalidHedge(r));
            }
        }

        let mut init = Self::empty(graph);

        let mut iter_order = 0;
        init.forest[RootId(root_node.0)] = TTRoot::Root;
        while let Some(hedge) = queue.pop_front() {
            // if the hedge is not external get the neighbors of the paired hedge
            if let Some(cn) = graph.connected_neighbors(subgraph, hedge) {
                let connected = graph.inv(hedge);

                if !seen.includes(&connected) && subgraph.includes(&connected) {
                    // if this new hedge hasn't been seen before, it means the node it belongs to
                    //  a new node in the traversal
                    //
                    init.tree_subgraph.add(connected);
                    init.tree_subgraph.add(hedge);
                    let node_id = init.forest.change_to_root(connected.into());
                    iter_order += 1;
                    init.forest[node_id] = TTRoot::Child(iter_order);
                } else {
                    continue;
                }
                // mark the new node as seen
                seen.union_with(&cn);

                // for all hedges in this new node, they have a parent, the initial hedge
                for i in cn.included_iter() {
                    // if they lead to a new node, they are potential branches, add them to the queue
                    if !seen.includes(&graph.inv(i)) && subgraph.includes(&i) {
                        queue.push_back(i);
                    }
                }
            }
        }
        Ok(init)
    }
}

#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
// #[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
/// Represents a traversal tree that owns its graph representation, typically simplified.
///
/// This structure might be used for scenarios where a traversal tree needs to be
/// stored or manipulated independently of the original, more complex graph from
/// which it was derived. The graph it owns might store simplified data (e.g., `()` for
/// edge data, `bool` for node data indicating inclusion).
///
/// # Type Parameters
///
/// - `P`: The type of [`ForestNodeStore`] used by the underlying `Forest` within the
///   owned `HedgeGraph`. `P` must also implement `Clone` and [`ForestNodeStorePreorder`].
pub struct OwnedTraversalTree<P: ForestNodeStore<NodeData = ()> + Clone + ForestNodeStorePreorder> {
    /// The graph representation owned by this traversal tree.
    /// Edge data is `()`, and node data is `bool` (likely indicating tree membership or similar).
    /// The node storage of this graph is a `Forest<bool, P>`.
    pub graph: HedgeGraph<(), bool, (), Forest<bool, P>>,
    /// An [`InternalSubGraph`] representing the edges that form this traversal tree
    /// within the context of its owned `graph`.
    pub tree_subgraph: InternalSubGraph,
    /// A bitmask indicating the set of half-edges covered by this traversal tree,
    /// potentially in the context of a larger original graph.
    // #[cfg_attr(feature = "bincode", bincode(with_serde))]
    pub covers: SuBitGraph,
}

// impl TraversalTree {
//     pub fn children(&self, hedge: Hedge) -> SuBitGraph {
//         let mut children = <SuBitGraph as SubGraph>::empty(self.inv.inv.len());

//         for (i, m) in self.parents.iter().enumerate() {
//             if let Parent::Hedge { hedge_to_root, .. } = m {
//                 if *hedge_to_root == hedge {
//                     children.set(i, true);
//                 }
//             }
//         }

//         children.set(hedge.0, false);

//         children
//     }

//     pub fn leaf_edges(&self) -> SuBitGraph {
//         let mut leaves = <SuBitGraph as SubGraph>::empty(self.inv.inv.len());
//         for hedge in self.covers().included_iter() {
//             let is_not_parent = !self.parent_iter().any(|(_, p)| {
//                 if let Parent::Hedge { hedge_to_root, .. } = p {
//                     *hedge_to_root == hedge
//                 } else {
//                     false
//                 }
//             });
//             if is_not_parent {
//                 leaves.set(hedge.0, true);
//             }
//         }
//         leaves
//     }

//     pub fn leaf_nodes<V, N: NodeStorageOps<NodeData = V>, E>(
//         &self,
//         graph: &HedgeGraph<E, V, N>,
//     ) -> Vec<NodeIndex> {
//         let mut leaves = IndexSet::new();

//         for hedge in self.covers().included_iter() {
//             if let Parent::Hedge { hedge_to_root, .. } = self.parent(hedge) {
//                 if *hedge_to_root == hedge {
//                     let mut sect = self
//                         .tree
//                         .filter
//                         .intersection(&graph.node_hairs(hedge).hairs);

//                     sect.set(hedge.0, false);

//                     if sect.count_ones() == 0 {
//                         leaves.insert(graph.node_id(hedge));
//                     }
//                 }
//             }
//         }

//         leaves.into_iter().collect()
//     }

//     pub fn child_nodes<V, N: NodeStorageOps<NodeData = V>, E>(
//         &self,
//         parent: NodeIndex,
//         graph: &HedgeGraph<E, V, N>,
//     ) -> Vec<NodeIndex> {
//         let mut children = IndexSet::new();

//         for h in graph.hairs_from_id(parent).hairs.included_iter() {
//             if let Parent::Hedge { hedge_to_root, .. } = self.parent(h) {
//                 if *hedge_to_root != h {
//                     if let Some(c) = graph.involved_node_id(h) {
//                         children.insert(c);
//                     }
//                 }
//             }
//         }

//         children.into_iter().collect()
//     }
// }
#[cfg(test)]
pub mod tests {
    use crate::{dot, half_edge::involution::Orientation, parser::DotGraph};

    use super::*;

    #[test]
    fn double_pentagon_tree() {
        let mut graph: DotGraph = dot!(
            digraph {
        node [shape=circle,height=0.1,label=""];  overlap="scale"; layout="neato";
        00 -> 07[ dir=none,label="a"];
        00 -> 12[ dir=forward,label="d"];
        01 -> 00[ dir=forward,label="d"];
        05 -> 02[ dir=forward,label="d"];
        02 -> 01[ dir=forward,label="d"];
        09 -> 04[ dir=forward,label="d"];
        02 -> 06[ dir=none,label="a"];
        01 -> 03[ dir=none,label="a"];
        03 -> 13[ dir=forward,label="d"];
        04 -> 03[ dir=forward,label="d"];
        04 -> 05[ dir=none,label="g"];
        06 -> 07[ dir=forward,label="e-"];
        07 -> 11[ dir=forward,label="e-"];
        08 -> 06[ dir=forward,label="e-"];
        10 -> 05[ dir=forward,label="d"];
        }
        )
        .unwrap();

        // println!("{}", graph.dot_display(&graph.full_filter()));

        let tree = SimpleTraversalTree::depth_first_traverse(
            &graph,
            &graph.full_filter(),
            &NodeIndex(0),
            None,
        )
        .unwrap();
        // println!("{}", tree.dot(&graph));

        graph.align_underlying_to_tree(&tree);
        graph.align_superficial_to_tree(&tree);
        assert!(graph
            .iter_edges()
            .all(|(_, _, d)| !matches!(d.orientation, Orientation::Reversed)));
        // println!("{}", tree.dot(&graph))
        // Test implementation
    }
}
