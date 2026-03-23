//! # Half-Edge Graph Representation
//!
//! This module provides the core data structures and algorithms for representing
//! and manipulating graphs using the **half-edge** (or doubly connected edge list - DCEL)
//! data structure. This representation is particularly useful for algorithms that
//! require efficient traversal of graph topology, such as iterating around faces,
//! finding incident edges to a node, or quickly accessing an edge's opposite pair.
//!
//! ## Key Concepts and Components:
//!
//! ### 1. `HedgeGraph<E, V, S>`
//! This is the central struct representing a graph.
//!   - `E`: Generic type for data associated with edges.
//!   - `V`: Generic type for data associated with nodes (vertices).
//!   - `S`: A type implementing `NodeStorage`, which defines how node data and their connectivity to half-edges are stored.
//!
//! `HedgeGraph` stores half-edges in an `edge_store` (typically `SmartHedgeVec`)
//! and nodes in a `node_store`.
//!
//! ### 2. `Hedge` and `NodeIndex`
//!   - `Hedge`: An identifier (often an index) for a single directed half-edge.
//!     Each undirected edge in the graph is typically represented by two oppositely
//!     directed half-edges.
//!   - `NodeIndex`: A simple wrapper around `usize` used to identify nodes in the graph.
//!
//! ### 3. Involutions (`involution` module)
//! A crucial concept in this implementation is the **involution** on half-edges.
//! An involution is a function that, when applied twice, returns the original input.
//! In the context of a half-edge graph, the primary involution is finding the
//! **opposite half-edge**.
//!   - `Involution`: The struct managing this mapping.
//!   - `HedgePair`: Represents how hedges are paired (e.g., forming a full edge,
//!     being unpaired/dangling, or split).
//!   - `Flow`: Indicates the underlying directionality of a half-edge (Source or Sink)
//!     relative to its edge.
//!   - `Orientation`: Indicates the superficial orientation of an edge (Default, Reversed,
//!     Undirected).
//!
//! The `inv(hedge: Hedge) -> Hedge` method on `HedgeGraph` is the primary way to
//! access the opposite half-edge.
//!
//! ### 4. Node Storage (`nodestore` module)
//! This module defines how nodes are stored and how they relate to the half-edges
//! that originate from them.
//!   - `NodeStorage` trait: An abstraction for different node storage strategies.
//!   - `NodeStorageOps` trait: Provides operations on node storage.
//!   - Each node typically stores a collection of its outgoing (or incident) half-edges.
//!
//! ### 5. Subgraphs (`subgraph` module)
//! This module provides extensive support for defining, manipulating, and querying
//! subgraphs. A subgraph is essentially a subset of the half-edges of a parent graph.
//!   - `SubGraph` trait: Defines the basic interface for a subgraph.
//!   - `SubGraphOps` trait: Provides operations applicable to subgraphs.
//!   - Various concrete subgraph types exist, such as:
//!     - `InternalSubGraph`: Represents a subgraph composed of edges internal to a region.
//!     - `HedgeNode`: Represents a "node" in a higher-level graph structure, possibly
//!       formed by contracting a subgraph. It contains an `internal_graph` and `hairs`
//!       (external connections).
//!     - `Cycle`: Represents a cycle in the graph.
//!     - `OrientedCut`: Represents a cut separating the graph.
//!
//! Operations include creating subgraphs from filters, extracting subgraphs,
//! performing boolean operations (union, intersection, difference), and analyzing
//! subgraph properties.
//!
//! ### 6. Graph Operations
//! The `HedgeGraph` provides a rich API for:
//!   - **Construction**: Building graphs using `HedgeGraphBuilder`, adding nodes and edges
//!     (`add_dangling_edge`, `add_pair`).
//!   - **Modification**: Deleting hedges (`delete_hedges`), splitting edges (`split_edge`),
//!     joining graphs (`join`), sewing dangling edges (`sew`), and identifying/contracting
//!     nodes (`identify_nodes`).
//!   - **Traversal & Access**: Getting an edge's opposite (`inv`), finding the node a
//!     hedge points to (`node_id`), iterating over neighbors (`iter_crown`), and
//!     accessing edge/node data.
//!   - **Analysis**: Checking connectivity (`is_connected`), counting components
//!     (`count_connected_components`), finding cycle bases (`cycle_basis`), spanning trees
//!     (`all_spanning_forests_of`), and calculating cyclotomatic numbers.
//!   - **Serialization/Visualization**: Generating DOT format strings for graph visualization.
//!
//! This module forms the backbone of the `linnet` library, enabling complex graph
//! manipulations with a focus on subgraphs and topological features.

use core::panic;
use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::fmt::Display;
use std::hash::Hash;
use std::num::TryFromIntError;
use std::ops::{Index, IndexMut, RangeBounds};

use ahash::{AHashMap, AHashSet};

use bitvec::prelude::*;
use bitvec::slice::IterOnes;
use builder::{HedgeData, HedgeGraphBuilder};
use hedgevec::{Accessors, SmartEdgeVec};
use indexmap::IndexSet;
use involution::{
    EdgeData, EdgeIndex, EdgeVec, Flow, Hedge, HedgePair, HedgeVec, Involution, InvolutionError,
    InvolutiveMapping, Orientation,
};

use itertools::Itertools;
use nodestore::{DefaultNodeStore, NodeStorage, NodeStorageOps};
use rand::{rngs::SmallRng, Rng, SeedableRng};
use swap::Swap;

define_indexed_vec!(
    /// A type-safe wrapper around a `usize` to represent the index of a node
    /// in a graph.
    ///
    /// This helps prevent accidental misuse of raw `usize` values where a node
    /// identifier is expected.
    pub struct NodeIndex; // new‑type around usize

    /// Vector whose items are `T`, indexable only by `EntityId`.
    pub struct NodeVec;
);

impl NodeIndex {
    pub fn add_data<H>(self, data: H) -> HedgeData<H> {
        HedgeData {
            data,
            node: self,
            is_in_subgraph: false,
        }
    }
}

impl std::fmt::Display for NodeIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// Iterator over the powerset of a bitvec, of size n < 64.
///
/// Generates all possible subsets of a set of up to 63 elements.
/// The subsets are represented as `SuBitGraph`s.
///
/// **Note:** The maximum number of elements is limited due to the `usize`
/// representation of the current subset index. On a 64-bit system, this
/// could theoretically support up to 64 elements, but practically, iterating
/// 2^64 items is not feasible. The comment "n < 64" suggests a practical limit.
pub struct PowersetIterator<ID = Hedge> {
    /// The total number of subsets, equal to 2^(number of elements).
    size: usize,
    /// The current subset index, ranging from 0 to `size - 1`.
    /// Each bit in `current` corresponds to an element's presence in the subset.
    current: usize,
    len: u8,
    id: std::marker::PhantomData<ID>,
}

impl<ID> PowersetIterator<ID> {
    pub fn new(n_elements: u8) -> Self {
        PowersetIterator {
            size: 1 << n_elements,
            current: 0,
            len: n_elements,
            id: std::marker::PhantomData,
        }
    }
}

impl<ID> Iterator for PowersetIterator<ID> {
    type Item = SubSet<ID>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current < self.size {
            let out = SubSet::from_usize(self.current, self.len as usize);
            self.current += 1;
            Some(out)
        } else {
            None
        }
    }
}
pub mod involution;

#[derive(Clone, Debug, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
/// Represents attributes for an edge when generating Graphviz DOT output.
///
/// This struct allows specifying common DOT attributes like label, color,
/// and other custom attributes for an edge.
pub struct GVEdgeAttrs {
    /// The label to be displayed on the edge.
    pub label: Option<String>,
    /// The color of the edge. Can be a color name (e.g., "red") or other
    /// Graphviz color specifications.
    pub color: Option<String>,
    /// A field for any other custom DOT attributes for the edge,
    /// represented as a single string (e.g., "style=dashed, arrowhead=vee").
    pub other: Option<String>,
}

impl std::fmt::Display for GVEdgeAttrs {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let out = [
            ("label=", self.label.as_ref()),
            (
                "color=",
                self.color.as_ref().map(|str| format!("\"{str}\"")).as_ref(),
            ),
            ("", self.other.as_ref()),
        ]
        .iter()
        .filter_map(|(prefix, x)| x.map(|s| format!("{prefix}{s}")))
        .join(" ")
        .to_string();
        if out.is_empty() {
            write!(f, "")
        } else {
            write!(f, " {out}")
        }
    }
}
pub mod builder;
pub mod nodestore;
pub mod subgraph;
pub mod swap;

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(
    feature = "bincode",
    derive(bincode_trait_derive::Encode, bincode_trait_derive::Decode)
)]
/// The main graph data structure, representing a graph using the half-edge
/// (or doubly connected edge list - DCEL) principle.
///
/// A `HedgeGraph` stores nodes and edges (represented as pairs of half-edges).
/// It is generic over the data associated with edges (`E`), nodes/vertices (`V`),
/// and the specific storage strategy used for nodes (`S`).
///
/// The half-edge representation allows for efficient traversal of graph topology,
/// such as iterating around faces (in planar embeddings), finding incident edges
/// to a node, or quickly accessing an edge's opposite half-edge.
///
/// # Type Parameters
///
/// - `E`: The type of data associated with each edge (or pair of half-edges).
/// - `V`: The type of data associated with each node (vertex).
/// - `S`: The node storage strategy, implementing the [`NodeStorage`] trait.
///   This determines how node data and their connectivity to half-edges
///   are stored. Defaults to [`DefaultNodeStore<V>`] (feature-selected; `nodestore-vec` uses
///   [`NodeStorageVec<V>`]).
pub struct HedgeGraph<E, V, H = NoData, S: NodeStorage<NodeData = V> = DefaultNodeStore<V>> {
    hedge_data: HedgeVec<H>,
    /// Internal storage for all half-edges, their data, and their topological
    /// relationships (e.g., opposite half-edge, next half-edge around a node).
    /// This is typically a [`SmartHedgeVec<E>`].
    edge_store: SmartEdgeVec<E>,
    /// Storage for all nodes in the graph, including their data (`V`) and
    /// information about the half-edges incident to them.
    /// The specific implementation is determined by the `S` type parameter.
    pub node_store: S,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(
    feature = "bincode",
    derive(bincode_trait_derive::Encode, bincode_trait_derive::Decode)
)]
pub struct NoData {}
impl Display for NoData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "")
    }
}

impl<E, V, H, S: NodeStorage<NodeData = V>> AsRef<Involution> for HedgeGraph<E, V, H, S> {
    fn as_ref(&self) -> &Involution {
        self.edge_store.as_ref()
    }
}

impl<E, V: Default, H: Default, N: NodeStorageOps<NodeData = V>> HedgeGraph<E, V, H, N> {
    /// Creates a random graph with a specified number of nodes and edges.
    ///
    /// The graph generated will have `()` as edge data and default values for node data `V`.
    /// The structure of the graph (connectivity) is determined randomly based on the provided seed.
    ///
    /// # Parameters
    /// - `nodes`: The desired number of nodes in the graph.
    /// - `edges`: The desired number of edges (pairs of half-edges) in the graph.
    /// - `seed`: A seed for the random number generator to ensure reproducibility.
    ///
    /// # Returns
    /// A new `HedgeGraph<(), V, N>` with a random structure.
    ///
    /// # Constraints
    /// - `N::Neighbors` must implement `BaseSubgraph` and `SubGraphOps`.
    /// - `V` must implement `Default`.
    ///
    /// # Panics
    /// Panics if `sources.len()` (derived from `edges`) is 0 when `externals` (unpaired half-edges)
    /// is not empty, or if `sinks.len()` is 0 when `externals` is not empty. This can occur if
    /// `edges` is too small relative to the internal logic of pairing, leading to attempts to access
    /// empty vectors `sources` or `sinks`. Also panics if node merging logic (to reduce node count
    /// to `nodes`) attempts to operate on empty or insufficiently small `sources` or `sinks` lists.
    pub fn random(nodes: usize, edges: usize, seed: u64) -> HedgeGraph<(), V, H, N>
    where
        N::Neighbors: BaseSubgraph + SubSetOps + ModifySubSet<Hedge>,
        <<N as NodeStorage>::Neighbors as SubSetLike>::Base: SubSetOps,
    {
        let inv: Involution<()> = Involution::<()>::random(edges, seed);

        let mut rng = SmallRng::seed_from_u64(seed);

        let mut externals = Vec::new();
        let mut sources = Vec::new();
        let mut sinks = Vec::new();

        for (i, e) in inv.iter() {
            let n_h: Hedge = inv.len();
            let mut nodeid = N::Neighbors::empty(n_h.0);
            nodeid.add(i);
            match e {
                InvolutiveMapping::Identity { .. } => externals.push(nodeid),
                InvolutiveMapping::Source { .. } => sources.push(nodeid),
                InvolutiveMapping::Sink { .. } => sinks.push(nodeid),
            }
        }

        assert_eq!(sources.len(), sinks.len());

        while !externals.is_empty() {
            if rng.gen_bool(0.5) {
                let source_i = rng.gen_range(0..sources.len());

                sources[source_i].union_with(&externals.pop().unwrap());
            } else {
                let sink_i = rng.gen_range(0..sinks.len());

                sinks[sink_i].union_with(&externals.pop().unwrap());
            }
        }

        let mut lengthone = false;

        while sources.len() + sinks.len() > nodes {
            if rng.gen_bool(0.5) {
                if sources.len() <= 1 {
                    if lengthone {
                        break;
                    }
                    lengthone = true;
                    continue;
                }

                let idx1 = rng.gen_range(0..sources.len());
                let idx2 = rng.gen_range(0..sources.len() - 1);

                let n_i = sources.swap_remove(idx1);
                sources[idx2].union_with(&n_i);
            } else {
                if sinks.len() <= 1 {
                    if lengthone {
                        break;
                    }
                    lengthone = true;
                    continue;
                }
                let idx1 = rng.gen_range(0..sinks.len());
                let idx2 = rng.gen_range(0..sinks.len() - 1);
                let n_i = sinks.swap_remove(idx1);
                sinks[idx2].union_with(&n_i);
            }
        }

        let mut hedge_data = HedgeVec::new();
        let n_h: Hedge = inv.len();
        for _ in 0..n_h.0 {
            hedge_data.push(H::default());
        }

        HedgeGraph {
            hedge_data,
            node_store: N::random(&sources, &sinks),
            edge_store: SmartEdgeVec::new(inv),
        }
    }
}

impl<E, V, H, N: NodeStorageOps<NodeData = V>> HedgeGraph<E, V, H, N> {
    pub fn check(&self) -> Result<(), HedgeGraphError> {
        // Check involution
        for (h, _) in self.iter_hedges() {
            if h != self.inv(self.inv(h)) {
                Err(InvolutionError::NonInvolutive(h))?;
            }
        }

        // Check smart edge vec
        self.edge_store.check_hedge_pairs()?;
        Ok(())
    }

    /// Deletes all half-edges specified in the `subgraph` from the graph.
    ///
    /// This operation modifies both the `edge_store` and `node_store` to remove
    /// the specified half-edges and update connectivity information.
    ///
    /// # Parameters
    /// - `subgraph`: A subgraph specifying the set of half-edges to delete.
    ///   `S` must implement `SubGraph` with `Base = N::Base`.
    pub fn delete_hedges<S: SubSetLike<Base = N::Base>>(&mut self, subgraph: &S) {
        let mut left = Hedge(0);
        let mut extracted: Hedge = self.len();
        while left < extracted {
            if !subgraph.includes(&left) {
                //left is in the right place
                left.0 += 1;
            } else {
                //left needs to be swapped
                extracted.0 -= 1;
                if !subgraph.includes(&extracted) {
                    // println!("{extracted}<=>{left}");
                    //only with an extracted that is in the wrong spot
                    self.hedge_data.swap(left, extracted);
                    left.0 += 1;
                }
            }
        }
        let _ = self.hedge_data.split_off(left);

        self.edge_store.delete(subgraph);
        self.node_store.delete(subgraph);
    }

    /// Creates a new `HedgeGraph` instance representing a "concretized" view of a given `subgraph`.
    ///
    /// This method effectively extracts the specified `subgraph` into a new, independent graph.
    /// The new graph contains copies of the nodes and edges (and their data via references)
    /// that are part of the `subgraph`.
    ///
    /// Node and edge data in the new graph are references (`&'a V`, `&'a E`) to the data
    /// in the original graph.
    ///
    /// # Parameters
    /// - `subgraph`: A reference to a subgraph `S` within the current graph.
    ///
    /// # Returns
    /// A new `HedgeGraph` containing only the elements of the `subgraph`.
    /// The node storage type of the new graph is `N::OpStorage<&'a V>`.
    pub fn concretize<'a, S: SubSetLike>(
        &'a self,
        subgraph: &'a S,
    ) -> HedgeGraph<&'a E, &'a V, &'a H, N::OpStorage<&'a V>> {
        let mut builder = HedgeGraphBuilder::new();

        let mut node_map = AHashMap::new();

        for (n, _, d) in self.iter_nodes_of(subgraph) {
            node_map.insert(n, builder.add_node(d));
        }

        for (pair, _, d) in self.iter_edges_of(subgraph) {
            match pair {
                HedgePair::Paired { source, sink } => {
                    let src = node_map[&self.node_id(source)].add_data(&self[source]);
                    let dst = node_map[&self.node_id(sink)].add_data(&self[sink]);

                    builder.add_edge(src, dst, d.data, d.orientation);
                }
                HedgePair::Unpaired { hedge, flow } => {
                    let src = node_map[&self.node_id(hedge)].add_data(&self[hedge]);

                    builder.add_external_edge(src, d.data, d.orientation, flow);
                }
                HedgePair::Split {
                    source,
                    sink,
                    split,
                } => match split {
                    Flow::Sink => {
                        let src = node_map[&self.node_id(sink)].add_data(&self[sink]);
                        builder.add_external_edge(src, d.data, d.orientation, split);
                    }
                    Flow::Source => {
                        let src = node_map[&self.node_id(source)].add_data(&self[source]);
                        builder.add_external_edge(src, d.data, d.orientation, split);
                    }
                },
            }
        }

        builder.build()
    }

    /// Extracts a subgraph and transforms its data, potentially splitting edges and nodes
    /// that are on the boundary of the subgraph.
    ///
    /// This is a more advanced operation than `concretize`. It allows for:
    /// - Transforming edge data for edges fully within the subgraph (`internal_data` closure).
    /// - Transforming edge data for edges that are split by the subgraph boundary (`split_edge_fn` closure).
    /// - Transforming node data for nodes that are part of split edges (`split_node` closure).
    /// - Transforming node data for nodes fully within the subgraph and taking ownership (`owned_node` closure).
    ///
    /// The result is a new `HedgeGraph` with transformed edge data `O` and node data `V2`.
    ///
    /// # Parameters
    /// - `subgraph`: The subgraph `S` to extract.
    /// - `split_edge_fn`: A closure that transforms `EdgeData<&E>` (for edges on the boundary)
    ///   to `EdgeData<O>`.
    /// - `internal_data`: A closure that transforms `EdgeData<E>` (for edges fully inside)
    ///   to `EdgeData<O>`.
    /// - `split_node`: A closure that transforms node data `&V` (for nodes on the boundary)
    ///   to `V2`.
    /// - `owned_node`: A closure that transforms node data `V` (for nodes fully inside, taking ownership)
    ///   to `V2`.
    ///
    /// # Returns
    /// A new `HedgeGraph<O, V2, N::OpStorage<V2>>` representing the extracted and transformed subgraph.
    pub fn extract<O, V2, S: SubSetLike<Base = N::Base>>(
        &mut self,
        subgraph: &S,
        split_edge_fn: impl FnMut(EdgeData<&E>) -> EdgeData<O>,
        internal_data: impl FnMut(EdgeData<E>) -> EdgeData<O>,
        split_node: impl FnMut(&V) -> V2,
        owned_node: impl FnMut(V) -> V2,
    ) -> HedgeGraph<O, V2, H, N::OpStorage<V2>> {
        let new_edge_store = self
            .edge_store
            .extract(subgraph, split_edge_fn, internal_data);

        let mut left = Hedge(0);
        let mut extracted: Hedge = self.len();
        while left < extracted {
            if !subgraph.includes(&left) {
                //left is in the right place
                left.0 += 1;
            } else {
                //left needs to be swapped
                extracted.0 -= 1;
                if !subgraph.includes(&extracted) {
                    // println!("{extracted}<=>{left}");
                    //only with an extracted that is in the wrong spot
                    self.hedge_data.swap(left, extracted);
                    left.0 += 1;
                }
            }
        }
        let new_hedge_data = self.hedge_data.split_off(left);

        let new_node_store = self.node_store.extract(subgraph, split_node, owned_node);
        HedgeGraph {
            hedge_data: new_hedge_data,
            node_store: new_node_store,
            edge_store: new_edge_store,
        }
    }

    pub fn extract_nodes<O>(
        &mut self,
        nodes: impl IntoIterator<Item = NodeIndex>,
        split_edge_fn: impl FnMut(EdgeData<&E>) -> EdgeData<O>,
        internal_data: impl FnMut(EdgeData<E>) -> EdgeData<O>,
    ) -> HedgeGraph<O, V, H, N>
    where
        N: NodeStorageOps<Base = SuBitGraph>,
    {
        let (extracted, nodes) = self.node_store.extract_nodes(nodes);

        let left = self.hedge_data.partition(|h| !extracted.includes(h));

        let ne_hedge = self.hedge_data.split_off(left);

        let new_edge_store = self
            .edge_store
            .extract(&extracted, split_edge_fn, internal_data);
        HedgeGraph {
            hedge_data: ne_hedge,
            node_store: nodes,
            edge_store: new_edge_store,
        }
    }
    /// Gives the involved hedge.
    /// Returns the opposite (or "twin") half-edge of the given `hedge`.
    ///
    /// - If `hedge` is part of a paired edge, this returns its sibling half-edge.
    /// - If `hedge` is an identity (unpaired) half-edge, it returns itself.
    ///
    /// # Parameters
    /// - `hedge`: The half-edge for which to find the opposite.
    ///
    /// # Returns
    /// The opposite `Hedge`.
    pub fn inv(&self, hedge: Hedge) -> Hedge {
        self.edge_store.inv(hedge)
    }

    /// Splits a paired edge (identified by one of its half-edges, `hedge`) into two
    /// new dangling (identity) half-edges.
    ///
    /// The original edge data is typically associated with one of the new dangling edges,
    /// and the provided `data` is associated with the other. The underlying `Flow` of the
    /// new half-edges corresponds to their roles in the original paired edge.
    ///
    /// # Parameters
    /// - `hedge`: One of the half-edges of the paired edge to be split.
    /// - `data`: The `EdgeData<E>` to be associated with the new dangling half-edge
    ///   created from the `hedge` side of the split. The original edge's data will be
    ///   associated with the other new dangling half-edge (created from `inv(hedge)`).
    ///
    /// # Returns
    /// - `Ok(())` if the edge was successfully split.
    /// - `Err(HedgeGraphError::InvolutionError(InvolutionError::NotPaired))` if `hedge` is an identity edge.
    pub fn split_edge(&mut self, hedge: Hedge, data: EdgeData<E>) -> Result<(), HedgeGraphError> {
        Ok(self.edge_store.split_edge(hedge, data)?)
    }

    /// Consumes two graphs (`self` and `other`) and joins them into a new graph.
    ///
    /// Dangling half-edges from both graphs are matched using `matching_fn`. If a match
    /// is found, these two dangling edges are "sewn" together to form a new paired edge
    /// in the resulting graph. The data for this new edge is determined by `merge_fn`.
    ///
    /// Node stores are extended, and the edge store is created by joining the two
    /// original edge stores.
    ///
    /// # Parameters
    /// - `other`: The other graph to join with `self`.
    /// - `matching_fn`: A closure `(Flow, EdgeData<&E>, Flow, EdgeData<&E>) -> bool`
    ///   that determines if two dangling half-edges (one from `self`, one from `other`)
    ///   should be matched and sewn together. It takes the flow and edge data of both.
    /// - `merge_fn`: A closure `(Flow, EdgeData<E>, Flow, EdgeData<E>) -> (Flow, EdgeData<E>)`
    ///   that takes the flow and owned edge data of two matched dangling half-edges and
    ///   returns the flow and edge data for the new combined edge.
    ///
    /// # Returns
    /// - `Ok(Self)`: The new graph resulting from the join.
    /// - `Err(HedgeGraphError)`: If an error occurs during the join process (e.g.,
    ///   from `self.edge_store.join` or if node validation fails).
    pub fn join(
        mut self,
        other: Self,
        matching_fn: impl Fn(Flow, EdgeData<&E>, Flow, EdgeData<&E>) -> bool,
        merge_fn: impl Fn(Flow, EdgeData<E>, Flow, EdgeData<E>) -> (Flow, EdgeData<E>),
    ) -> Result<Self, HedgeGraphError> {
        self.hedge_data.extend(other.hedge_data);
        let mut g = HedgeGraph {
            hedge_data: self.hedge_data,
            node_store: self.node_store.extend(other.node_store),
            edge_store: self
                .edge_store
                .join(other.edge_store, matching_fn, merge_fn)?,
        };
        g.node_store.check_and_set_nodes()?;

        Ok(g)
    }

    /// Joins another graph (`other`) into `self` by consuming `other`.
    ///
    /// This is similar to `join`, but modifies `self` in place instead of returning a new graph.
    /// Dangling edges are matched and merged according to `matching_fn` and `merge_fn`.
    ///
    /// # Parameters
    /// - `other`: The graph to be consumed and joined into `self`.
    /// - `matching_fn`: A closure to determine if two dangling half-edges should match.
    /// - `merge_fn`: A closure to determine the data for the newly formed paired edge.
    ///
    /// # Returns
    /// - `Ok(())` if the join was successful.
    /// - `Err(HedgeGraphError)` if an error occurs (e.g. from `self.edge_store.join_mut` or
    ///   if node validation fails).
    pub fn join_mut(
        &mut self,
        other: Self,
        matching_fn: impl Fn(Flow, EdgeData<&E>, Flow, EdgeData<&E>) -> bool,
        merge_fn: impl Fn(Flow, EdgeData<E>, Flow, EdgeData<E>) -> (Flow, EdgeData<E>),
    ) -> Result<(), HedgeGraphError> {
        // self.node_store.check_and_set_nodes()?;

        // other.node_store.check_and_set_nodes()?;
        self.node_store.extend_mut(other.node_store);
        self.edge_store
            .join_mut(other.edge_store, matching_fn, merge_fn)?;

        self.hedge_data.extend(other.hedge_data);

        self.node_store.check_and_set_nodes()?;
        Ok(())
    }

    /// Sews dangling edges internal to the graph, matching edges with the given function and merging them with the given function.
    /// "Sews" together pairs of dangling (identity) half-edges within the graph `self`.
    ///
    /// This operation attempts to find pairs of dangling half-edges that satisfy `matching_fn`
    /// and converts them into internal, paired edges. The data for the new paired edge
    /// is determined by `merge_fn`.
    ///
    /// # Parameters
    /// - `matching_fn`: A closure `(Flow, EdgeData<&E>, Flow, EdgeData<&E>) -> bool`
    ///   that determines if two dangling half-edges within `self` should be matched.
    /// - `merge_fn`: A closure `(Flow, EdgeData<E>, Flow, EdgeData<E>) -> (Flow, EdgeData<E>)`
    ///   that provides the data for the new paired edge formed from two matched dangling edges.
    ///
    /// # Returns
    /// - `Ok(())` if the sewing operation completes.
    /// - `Err(HedgeGraphError::InvolutionError)` if an error occurs within the underlying
    ///   `edge_store.sew` operation (e.g. `InvolutionError::NotIdentity` if an attempt is made
    ///   to sew an edge that is not an identity edge).
    pub fn sew(
        &mut self,
        matching_fn: impl Fn(Flow, EdgeData<&E>, Flow, EdgeData<&E>) -> bool,
        merge_fn: impl Fn(Flow, EdgeData<E>, Flow, EdgeData<E>) -> (Flow, EdgeData<E>),
    ) -> Result<(), HedgeGraphError> {
        self.edge_store.sew(matching_fn, merge_fn)
    }

    /// Adds a new dangling (identity) half-edge to the graph, incident to the specified `source` node.
    ///
    /// This consumes the graph and returns a new graph with the added edge.
    ///
    /// # Parameters
    /// - `source`: The `NodeIndex` to which the new dangling edge will be attached.
    /// - `data`: The data `E` for the new edge.
    /// - `flow`: The [`Flow`] (Source or Sink) of the new dangling half-edge.
    /// - `orientation`: The [`Orientation`] of the new edge.
    ///
    /// # Returns
    /// - `Ok((Hedge, Self))`: A tuple containing the `Hedge` identifier of the new
    ///   dangling half-edge and the modified graph.
    /// - `Err(HedgeGraphError::NoNode)`: If the `source` node index is invalid for the node store.
    /// - Other `HedgeGraphError` variants may arise from internal node store operations.
    pub fn add_dangling_edge(
        mut self,
        source: impl Into<HedgeData<H>>,
        data: E,
        flow: Flow,
        orientation: impl Into<Orientation>,
    ) -> Result<(Hedge, Self), HedgeGraphError> {
        let source = source.into();

        self.hedge_data.push(source.data);
        let (edge_store, hedge) = self.edge_store.add_dangling_edge(data, flow, orientation);
        let mut g = HedgeGraph {
            hedge_data: self.hedge_data,
            edge_store,
            node_store: self.node_store.add_dangling_edge(source.node)?,
        };

        g.node_store.check_and_set_nodes()?;

        Ok((hedge, g))
    }

    /// Adds a new paired edge between the `source` and `sink` nodes.
    ///
    /// This consumes the graph and returns a new graph with the added edge.
    /// The new edge consists of two half-edges, one originating from `source` and
    /// the other from `sink`, pointing towards each other.
    ///
    /// # Parameters
    /// - `source`: The `NodeIndex` of the source node for the new edge.
    /// - `sink`: The `NodeIndex` of the sink node for the new edge.
    /// - `data`: The data `E` for the new edge.
    /// - `orientation`: The [`Orientation`] of the new edge.
    ///
    /// # Returns
    /// - `Ok((Hedge, Hedge, Self))`: A tuple containing the `Hedge` identifiers for the
    ///   source-incident half-edge, the sink-incident half-edge, and the modified graph.
    /// - `Err(HedgeGraphError::NoNode)`: If `source` or `sink` node indices are invalid.
    /// - Other `HedgeGraphError` variants may arise from internal node store operations.
    pub fn add_pair(
        mut self,
        source: impl Into<HedgeData<H>>,
        sink: impl Into<HedgeData<H>>,
        data: E,
        orientation: impl Into<Orientation>,
    ) -> Result<(Hedge, Hedge, Self), HedgeGraphError> {
        let source = source.into();
        let sink = sink.into();
        self.hedge_data.push(source.data);
        self.hedge_data.push(sink.data);
        let (edge_store, sourceh, sinkh) = self.edge_store.add_paired(data, orientation);
        let mut g = HedgeGraph {
            hedge_data: self.hedge_data,
            edge_store,
            node_store: self
                .node_store
                .add_dangling_edge(source.node)?
                .add_dangling_edge(sink.node)?,
        };

        g.node_store.check_and_set_nodes()?;

        Ok((sourceh, sinkh, g))
    }

    /// Checks if the specified `subgraph` is connected.
    ///
    /// A subgraph is connected if there is a path between any two half-edges (or nodes
    /// they are incident to) within that subgraph, using only edges also within the subgraph.
    ///
    /// # Parameters
    /// - `subgraph`: The subgraph `S` to check for connectivity.
    ///
    /// # Returns
    /// `true` if the subgraph is connected, `false` otherwise.
    /// Returns `true` for an empty subgraph.
    ///
    /// # Panics
    /// Panics if the traversal (used internally) fails, which can happen if the
    /// starting node for traversal is not part of the `subgraph`.
    pub fn is_connected<S: SubGraphLike>(&self, subgraph: &S) -> bool {
        let n_edges = subgraph.nedges(self);
        if let Some(start) = subgraph.included_iter().next() {
            SimpleTraversalTree::depth_first_traverse(self, subgraph, &self.node_id(start), None)
                .unwrap()
                .covers(subgraph)
                .nedges(self)
                == n_edges
        } else {
            true
        }
    }

    /// Modifies a [`HedgeNode`] subgraph by removing "branches" or "tendrils".
    ///
    /// A branch is typically a path of edges within the `subgraph` that ultimately
    /// connects to the rest of the `subgraph` at only one point (one node with degree > 1
    /// within the branch, relative to other branch edges). This operation iteratively
    /// removes edges that are part of such terminal paths until no more such branches exist.
    ///
    /// It first removes purely external edges from the `subgraph`'s internal part,
    /// then iteratively prunes edges that form degree-1 connections within the subgraph's context.
    /// Finally, it fixes the `hairs` of the `HedgeNode` to be consistent.
    ///
    /// # Parameters
    /// - `subgraph`: A mutable reference to the `HedgeNode` to be pruned.
    pub fn cut_branches(&self, subgraph: &mut HedgeNode) {
        let nodes = AHashSet::<NodeIndex>::from_iter(
            subgraph
                .internal_graph
                .included_iter()
                .map(|i| self.node_id(i)),
        );
        self.remove_externals(subgraph);

        let mut has_branch = true;
        while has_branch {
            has_branch = false;

            for n in &nodes {
                let int = self.iter_crown_in(subgraph, *n).collect::<Vec<_>>();
                let first = int.first();
                let next = int.get(1);

                if let Some(first) = first {
                    if next.is_none() {
                        subgraph.internal_graph.filter.sub(*first);
                        subgraph.internal_graph.filter.sub(self.inv(*first));
                        has_branch = true;
                    }
                }
            }
        }

        self.nesting_node_fix(subgraph);
    }

    // fn _set_hedge_data(&mut self, hedge: Hedge, nodeid: NodeIndex) {
    //     self.node_store.set_hedge_data(hedge, nodeid);
    // }
}

// Subgraphs
impl<E, V, H, N: NodeStorageOps<NodeData = V>> HedgeGraph<E, V, H, N> {
    /// Calculates the "internal crown" of a given `subgraph`.
    ///
    /// The internal crown consists of all half-edges within the `subgraph` that are
    /// either unpaired (dangling/external) or part of a "split" edge (i.e., their
    /// opposite half-edge is *not* in the `subgraph`). These are effectively the
    /// boundary half-edges of the `subgraph` from its own perspective.
    ///
    /// # Parameters
    /// - `subgraph`: The subgraph `S` for which to find the internal crown.
    ///
    /// # Returns
    /// A new subgraph of type `S::Base` containing the internal crown half-edges.
    /// `S::Base` must implement `ModifySubgraph<HedgePair>`.
    pub fn internal_crown<S: SubSetLike>(&self, subgraph: &S) -> S::Base
    where
        S::Base: ModifySubSet<HedgePair>,
    {
        let mut crown = S::Base::empty(self.n_hedges());

        for (p, _, _) in self.iter_edges_of(subgraph) {
            if !p.is_paired() {
                crown.add(p)
            }
        }

        crown
    }

    /// Calculates the "full crown" of a given `subgraph`.
    ///
    /// The full crown consists of all half-edges that are incident to any node
    /// touched by the `subgraph`, provided that these half-edges are either
    /// unpaired (identity) or their opposite half-edge is not included in the `subgraph`.
    ///
    /// This is different from `internal_crown` as it considers all incident edges to
    /// nodes in the subgraph's footprint, not just edges within the subgraph itself.
    /// It might include edges not present in the initial `subgraph`.
    ///
    /// # Parameters
    /// - `subgraph`: The subgraph `S` for which to find the full crown.
    ///
    /// # Returns
    /// A new subgraph of type `S::Base` containing the full crown half-edges.
    /// `S::Base` must implement `ModifySubgraph<Hedge>`.
    pub fn full_crown<S: SubSetLike>(&self, subgraph: &S) -> S::Base
    where
        S::Base: ModifySubSet<Hedge>,
    {
        let mut crown = S::Base::empty(self.n_hedges());

        for (_, n, _) in self.iter_nodes_of(subgraph) {
            for h in n {
                let invh = self.inv(h);
                if h == invh || !subgraph.includes(&invh) {
                    crown.add(h);
                }
            }
        }

        crown
    }

    /// Add all half-edges that are incident to any node
    /// touched by the `subgraph`, provided that these half-edges are either
    /// unpaired (identity) or their opposite half-edge is not included in the `subgraph`.
    pub fn add_crown<S>(&self, subgraph: &mut S)
    where
        S: ModifySubSet<Hedge> + SubSetLike,
    {
        let nodes: Vec<_> = self.iter_nodes_of(subgraph).map(|(n, _, _)| n).collect();
        for n in nodes {
            for h in self.iter_crown(n) {
                let invh = self.inv(h);
                if h == invh || !subgraph.includes(&invh) {
                    subgraph.add(h);
                }
            }
        }
    }
    /// Creates a `SuBitGraph` representing a subgraph containing all external (identity/dangling)
    /// half-edges in the entire graph.
    ///
    /// # Returns
    /// A `SuBitGraph` where bits corresponding to external half-edges are set to true.
    pub fn external_filter<S: ModifySubSet<Hedge> + SubSetLike>(&self) -> S {
        let mut filter: S = self.empty_subgraph();

        for (i, _, _) in self.iter_edges() {
            if i.is_unpaired() {
                filter.add(i.any_hedge());
            }
        }

        filter
    }

    /// Creates a `SuBitGraph` representing a subgraph containing all half-edges in the graph.
    ///
    /// # Returns
    /// A `SuBitGraph` of length `self.n_hedges()` with all bits set to true.
    pub fn full_filter(&self) -> SuBitGraph {
        SuBitGraph::full(self.n_hedges())
    }

    /// Returns a [`FullOrEmpty`] subgraph representing the entire graph (all hedges included).
    pub fn full(&self) -> FullOrEmpty {
        FullOrEmpty::full(self.n_hedges())
    }

    /// Returns a [`FullOrEmpty`] subgraph representing an empty graph (no hedges included).
    pub fn empty(&self) -> FullOrEmpty {
        FullOrEmpty::empty(self.n_hedges())
    }

    /// Creates an [`InternalSubGraph`] from a `SuBitGraph` filter, ensuring it has no "hairs".
    ///
    /// This uses a "pessimistic" approach: an edge is included only if both its
    /// half-edges are set in the input `filter`. Dangling edges are removed.
    ///
    /// # Parameters
    /// - `filter`: A `SuBitGraph` representing the desired set of half-edges.
    ///
    /// # Returns
    /// A new `InternalSubGraph`.
    pub fn clean_subgraph(&self, filter: SuBitGraph) -> InternalSubGraph {
        InternalSubGraph::cleaned_filter_pessimist(filter, self)
    }

    /// Returns a [`HedgeNode`] that represents the entire graph.
    /// The `internal_graph` of this `HedgeNode` will include all internal edges,
    /// and its `hairs` will include all external (dangling) edges of the graph.
    pub fn full_node(&self) -> HedgeNode {
        self.nesting_node_from_subgraph(self.full_graph())
    }

    /// Returns an [`InternalSubGraph`] that includes all fully internal edges of the graph.
    /// External (dangling) edges are excluded.
    pub fn full_graph(&self) -> InternalSubGraph {
        InternalSubGraph::cleaned_filter_optimist(self.full_filter(), self)
    }

    /// Creates an empty subgraph of a specific type `S`.
    ///
    /// # Type Parameters
    /// - `S`: The type of subgraph to create, must implement `SubGraph`.
    ///
    /// # Returns
    /// A new, empty subgraph of type `S`, sized for this graph.
    pub fn empty_subgraph<S: SubSetLike>(&self) -> S {
        S::empty(self.n_hedges())
    }

    /// Creates a subgraph of type `S` by filtering edges based on their data.
    ///
    /// # Type Parameters
    /// - `S`: The type of subgraph to create, must implement `BaseSubgraph`.
    ///
    /// # Parameters
    /// - `filter`: A closure that takes edge data `&E` and returns `true` if the
    ///   edge should be included in the subgraph.
    ///
    /// # Returns
    /// A new subgraph of type `S`.
    pub fn from_filter<S: BaseSubgraph>(&self, filter: impl FnMut(&E) -> bool) -> S
    where
        S::Base: SubSetOps<Hedge>,
    {
        S::from_filter(self, filter)
    }

    /// Creates a [`HedgeNode`] from a given [`InternalSubGraph`].
    ///
    /// The `internal_graph` of the new `HedgeNode` is the one provided.
    /// The `hairs` of the `HedgeNode` are calculated as all half-edges incident to the
    /// `internal_graph` that are not part of the `internal_graph` itself.
    ///
    /// # Parameters
    /// - `internal_graph`: The `InternalSubGraph` to form the core of the `HedgeNode`.
    ///
    /// # Panics
    /// Panics if the provided `internal_graph` is not valid for this graph (e.g.,
    /// if it refers to hedges outside the graph's bounds or is not truly internal).
    pub fn nesting_node_from_subgraph(&self, internal_graph: InternalSubGraph) -> HedgeNode {
        let mut hairs: SuBitGraph = self.empty_subgraph();

        if !internal_graph.valid::<E, V, H, N>(self) {
            panic!("Invalid subgraph")
        }

        for i in internal_graph.included_iter() {
            hairs.union_with_iter(self.neighbors(i));
        }

        let mut nh = !hairs;
        nh.union_with(&internal_graph.filter);
        HedgeNode {
            hairs: !nh,
            internal_graph,
        }
    }

    pub fn remove_internal_hedges(&self, subgraph: &SuBitGraph) -> SuBitGraph {
        let mut hairs = subgraph.clone();
        for i in subgraph.included_iter() {
            if subgraph.includes(&self.inv(i)) {
                hairs.sub(i);
                hairs.sub(self.inv(i));
            }
        }
        hairs
    }

    pub(crate) fn split_hairs_and_internal_hedges(
        &self,
        mut subgraph: SuBitGraph,
    ) -> (SuBitGraph, InternalSubGraph) {
        let mut internal: InternalSubGraph = self.empty_subgraph();
        for i in subgraph.included_iter() {
            let invh = self.inv(i);
            if subgraph.includes(&invh) {
                internal.filter.add(i);
                internal.filter.add(invh);
            }
        }
        for i in internal.filter.included_iter() {
            subgraph.sub(i);
        }
        (subgraph, internal)
    }

    /// Adjusts the `hairs` of a `HedgeNode` to be consistent with its `internal_graph`.
    ///
    /// Ensures that `hairs` only contains true external connections relative to the
    /// `internal_graph`. It recalculates hairs based on neighbors of the internal graph
    /// and removes any that are already part of the internal graph.
    ///
    /// # Parameters
    /// - `node`: A mutable reference to the `HedgeNode` to fix.
    fn nesting_node_fix(&self, node: &mut HedgeNode) {
        let mut externalhedges: SuBitGraph = self.empty_subgraph();

        for i in node.internal_graph.filter.included_iter() {
            externalhedges.union_with_iter(self.neighbors(i));
        }

        let mut ne = !externalhedges;
        ne.union_with(&node.internal_graph.filter);

        node.hairs = !ne;
    }

    fn remove_externals(&self, subgraph: &mut HedgeNode) {
        let externals: SuBitGraph = self.external_filter();

        subgraph.internal_graph.filter.subtract_with(&externals);
    }
}

// Counts
impl<E, V, H, N: NodeStorageOps<NodeData = V>> HedgeGraph<E, V, H, N> {
    /// Counts the number of full internal edges within the given `subgraph`.
    ///
    /// An edge is considered internal if both its half-edges are included in the `subgraph`.
    /// This method avoids double-counting by only counting an edge once.
    ///
    /// # Parameters
    /// - `subgraph`: The subgraph `S` in which to count internal edges.
    ///
    /// # Returns
    /// The number of full internal edges.
    pub fn count_internal_edges<S: SubSetLike>(&self, subgraph: &S) -> usize {
        let mut internal_edge_count = 0;
        // Iterate over all half-edges in the subgraph
        for hedge_index in subgraph.included_iter() {
            let inv_hedge_index = self.inv(hedge_index);

            // Check if the involuted half-edge is also in the subgraph
            if subgraph.includes(&inv_hedge_index) {
                // To avoid double-counting, only count when hedge_index < inv_hedge_index
                if hedge_index < inv_hedge_index {
                    internal_edge_count += 1;
                }
            }
        }
        internal_edge_count
    }

    /// Returns the total number of half-edges in the graph.
    pub fn n_hedges(&self) -> usize {
        let n_h: Hedge = self.len();
        n_h.0
    }

    pub fn n_edges(&self) -> usize {
        let n_e: EdgeIndex = self.len();
        n_e.0
    }

    /// Returns the total number of nodes in the graph.
    pub fn n_nodes(&self) -> usize {
        let n_n: NodeIndex = self.len();
        n_n.0
    }

    /// Returns the number of external (dangling/identity) half-edges in the graph.
    pub fn n_externals(&self) -> usize {
        self.edge_store.n_dangling()
    }

    /// Returns the number of internal (paired) half-edges in the graph.
    /// Note that this counts half-edges, so a single full internal edge contributes 2 to this count.
    pub fn n_internals(&self) -> usize {
        self.edge_store.n_paired()
    }

    // pub fn n_base_nodes(&self) -> usize {
    //     self.nodes.iter().filter(|n| n.is_node()).count()
    // }

    /// Counts the number of distinct nodes that are incident to at least one half-edge
    /// in the given `subgraph`.
    ///
    /// # Parameters
    /// - `subgraph`: The subgraph `S` to consider.
    ///
    /// # Returns
    /// The number of unique nodes touched by the `subgraph`.
    pub fn number_of_nodes_in_subgraph<S: SubSetLike>(&self, subgraph: &S) -> usize {
        self.iter_nodes_of(subgraph).count()
    }

    /// Calculates the degree of each node within the context of a given `InternalSubGraph`.
    ///
    /// The degree of a node in this context is the number of half-edges from the `subgraph`
    /// that are incident to that node.
    ///
    /// # Parameters
    /// - `subgraph`: The `InternalSubGraph` to calculate node degrees from.
    ///
    /// # Returns
    /// An `AHashMap` mapping each `NodeIndex` (for nodes involved in the `subgraph`)
    /// to its degree within that `subgraph`.
    pub fn node_degrees_in_subgraph(
        &self,
        subgraph: &InternalSubGraph,
    ) -> AHashMap<NodeIndex, usize> {
        let mut degrees = AHashMap::new();

        for (_, node, _) in self.iter_nodes_of(subgraph) {
            let node_pos = self.id_from_crown(node).unwrap();

            // Count the number of edges in the subgraph incident to this node
            let incident_edges = SuBitGraph::from_hedge_iter(
                self.iter_crown_in(subgraph, node_pos),
                subgraph.size(),
            );
            let degree = incident_edges.n_included();

            degrees.insert(node_pos, degree);
        }

        degrees
    }
}

pub trait EdgeAccessors<Index> {
    fn orientation(&self, index: Index) -> Orientation;

    fn set_orientation(&mut self, index: Index, orientation: Orientation);
}

impl<E, V, H, N: NodeStorageOps<NodeData = V>> EdgeAccessors<Hedge> for HedgeGraph<E, V, H, N> {
    fn orientation(&self, index: Hedge) -> Orientation {
        self.edge_store.orientation(index)
    }

    fn set_orientation(&mut self, index: Hedge, orientation: Orientation) {
        self.edge_store.set_orientation(index, orientation);
    }
}

impl<E, V, H, N: NodeStorageOps<NodeData = V>> EdgeAccessors<HedgePair> for HedgeGraph<E, V, H, N> {
    fn orientation(&self, index: HedgePair) -> Orientation {
        self.edge_store.orientation(index)
    }

    fn set_orientation(&mut self, index: HedgePair, orientation: Orientation) {
        self.edge_store.set_orientation(index, orientation);
    }
}

impl<E, V, H, N: NodeStorageOps<NodeData = V>> EdgeAccessors<EdgeIndex> for HedgeGraph<E, V, H, N> {
    fn orientation(&self, index: EdgeIndex) -> Orientation {
        self.edge_store.orientation(index)
    }

    fn set_orientation(&mut self, index: EdgeIndex, orientation: Orientation) {
        self.edge_store.set_orientation(index, orientation);
    }
}

// Accessors
impl<E, V, H, N: NodeStorageOps<NodeData = V>> HedgeGraph<E, V, H, N> {
    /// including pos
    pub fn owned_neighbors<S: SubGraphLike>(&self, subgraph: &S, pos: Hedge) -> SuBitGraph {
        subgraph.hairs(self.neighbors(pos))
    }

    pub fn connected_neighbors<S: SubGraphLike>(
        &self,
        subgraph: &S,
        pos: Hedge,
    ) -> Option<SuBitGraph> {
        Some(subgraph.hairs(self.involved_node_crown(pos)?))
    }
    pub fn get_edge_data(&self, edge: Hedge) -> &E {
        &self[self[&edge]]
    }

    pub fn hedge_pair(&self, hedge: Hedge) -> HedgePair {
        self.edge_store[&self.edge_store[&hedge]].1
    }

    pub fn get_edge_data_full(&self, hedge: Hedge) -> EdgeData<&E> {
        let orientation = self.edge_store.orientation(hedge);
        EdgeData::new(&self[self[&hedge]], orientation)
    }

    /// Gives the underlying orientation of this half-edge.
    pub fn flow(&self, hedge: Hedge) -> Flow {
        self.edge_store.flow(hedge)
    }

    pub fn superficial_hedge_orientation(&self, hedge: Hedge) -> Option<Flow> {
        self.edge_store.superficial_hedge_orientation(hedge)
    }

    pub fn underlying_hedge_orientation(&self, hedge: Hedge) -> Flow {
        self.edge_store.underlying_hedge_orientation(hedge)
    }

    pub fn neighbors(&self, hedge: Hedge) -> N::NeighborsIter<'_> {
        self.iter_crown(self.node_id(hedge))
    }

    pub fn iter_crown(&self, id: NodeIndex) -> N::NeighborsIter<'_> {
        self.node_store.get_neighbor_iterator(id)
    }

    pub fn iter_crown_in<'a, S: SubSetLike>(
        &'a self,
        subgraph: &'a S,
        id: NodeIndex,
    ) -> impl Iterator<Item = Hedge> + 'a {
        self.iter_crown(id).filter(|i| subgraph.includes(i))
    }

    pub fn id_from_crown<'a>(&'a self, mut neighbors: N::NeighborsIter<'a>) -> Option<NodeIndex> {
        let e = neighbors.next()?;
        Some(self.node_id(e))
    }

    pub fn involved_node_crown(&self, hedge: Hedge) -> Option<N::NeighborsIter<'_>> {
        self.involved_node_id(hedge).map(|id| self.iter_crown(id))
    }

    pub fn involved_node_id(&self, hedge: Hedge) -> Option<NodeIndex> {
        let invh = self.inv(hedge);
        if invh == hedge {
            return None;
        }
        Some(self.node_id(invh))
    }

    pub fn node_id(&self, hedge: Hedge) -> NodeIndex {
        self.node_store.node_id_ref(hedge)
    }

    pub fn is_self_loop(&self, hedge: Hedge) -> bool {
        !self.is_dangling(hedge) && self.node_id(hedge) == self.node_id(self.inv(hedge))
    }

    pub fn is_dangling(&self, hedge: Hedge) -> bool {
        self.inv(hedge) == hedge
    }

    /// Collect all nodes in the subgraph (all nodes that the hedges are connected to)
    pub fn nodes<S: SubSetLike>(&self, subgraph: &S) -> Vec<NodeIndex> {
        let mut nodes = IndexSet::new();
        for i in subgraph.included_iter() {
            let node = self.node_id(i);
            nodes.insert(node);
        }

        nodes.into_iter().collect()
    }

    pub fn set_flow(&mut self, hedge: Hedge, flow: Flow) {
        self.edge_store.set_flow(hedge, flow);
    }

    ///Permutes nodes not pointing to any root anymore to end of nodestore and then extract it
    pub fn forget_identification_history(&mut self) -> NodeVec<V> {
        self.node_store.forget_identification_history()
    }

    ///Retains the NodeIndex ordering and just appends a new node.
    pub fn identify_nodes(&mut self, nodes: &[NodeIndex], node_data_merge: V) -> NodeIndex {
        self.node_store.identify_nodes(nodes, node_data_merge)
    }

    ///Identifies all nodes in this subgraph and gives value node_data_merge to identified node.
    ///Deletes all edges in the subgraph.
    ///This invalidates both hedge indices and node indices
    pub fn contract_subgraph<S: SubSetLike<Base = N::Base>>(
        &mut self,
        subgraph: &S,
        node_data_merge: V,
    ) {
        let nodes: Vec<_> = self.iter_nodes_of(subgraph).map(|(a, _, _)| a).collect();

        self.identify_nodes(&nodes, node_data_merge);
        self.forget_identification_history();
        self.delete_hedges(subgraph);
    }

    ///Retains the NodeIndex ordering and just appends a new node.
    pub fn identify_nodes_without_self_edges<S>(
        &mut self,
        nodes: &[NodeIndex],
        node_data_merge: V,
    ) -> (NodeIndex, S)
    where
        S: ModifySubSet<Hedge> + SubSetLike,
    {
        let mut self_edges: S = self.empty_subgraph();
        for n in nodes {
            for h in self.iter_crown(*n) {
                if self.is_self_loop(h) {
                    self_edges.add(h);
                }
            }
        }
        let n = self.node_store.identify_nodes(nodes, node_data_merge);

        for h in self.iter_crown(n) {
            if self.is_self_loop(h) {
                if self_edges.includes(&h) {
                    self_edges.sub(h);
                } else {
                    self_edges.add(h);
                }
            }
        }

        // self_edges

        (n, self_edges)
    }

    /// Collect all edges in the subgraph
    /// (This is without double counting, i.e. if two half-edges are part of the same edge, only one `EdgeIndex` will be collected)
    pub fn edges<S: SubSetLike>(&self, subgraph: &S) -> Vec<EdgeIndex> {
        self.iter_edges_of(subgraph).map(|(_, i, _)| i).collect()
    }
}

pub enum NodeKind<Data> {
    Internal(Data),
    External { data: Data, flow: Flow },
}

pub enum DanglingMatcher<Data> {
    Actual { hedge: Hedge, data: Data },
    Internal { data: Data },
    Saturator { hedge: Hedge },
}

impl<Data> DanglingMatcher<Data> {
    fn new(pair: HedgePair, data: Data) -> Self {
        match pair {
            HedgePair::Unpaired { hedge, .. } => DanglingMatcher::Actual { hedge, data },
            HedgePair::Paired { .. } => DanglingMatcher::Internal { data },
            _ => panic!("Split"),
        }
    }
    pub fn matches(&self, other: &Self) -> bool {
        match (self, other) {
            (
                DanglingMatcher::Actual { hedge: h1, .. },
                DanglingMatcher::Saturator { hedge: h2 },
            ) => h1 == h2,
            _ => false,
        }
    }

    pub fn unwrap(self) -> Data {
        match self {
            DanglingMatcher::Actual { data, .. } => data,
            DanglingMatcher::Internal { data } => data,
            DanglingMatcher::Saturator { .. } => panic!("Cannot unwrap a saturator"),
        }
    }
}

// Mapping
impl<E, V, H, N: NodeStorageOps<NodeData = V>> HedgeGraph<E, V, H, N> {
    pub fn saturate_dangling<V2>(
        self,
        node_map: impl FnMut(&Involution, NodeIndex, V) -> V2,
        mut dangling_map: impl FnMut(&Involution, &N, Hedge, Flow, EdgeData<&E>) -> (V2, H),
    ) -> HedgeGraph<E, V2, H, <<N as NodeStorageOps>::OpStorage<V2> as NodeStorageOps>::OpStorage<V2>>
    {
        let ext: SuBitGraph = self.external_filter();

        let mut saturator = HedgeGraphBuilder::new();
        for i in ext.included_iter() {
            let flow = self.flow(i);
            let d = &self[[&i]];
            let orientation = self.orientation(i);
            let (new_node_data, h) = dangling_map(
                self.edge_store.as_ref(),
                &self.node_store,
                i,
                flow,
                EdgeData::new(d, orientation),
            );

            let n = saturator.add_node(new_node_data).add_data(h);

            saturator.add_external_edge(
                n,
                DanglingMatcher::Saturator { hedge: i },
                orientation,
                -flow,
            );
        }

        let saturator = saturator.build();
        let new_graph = self.map(
            node_map,
            |_, _, p, _, e| e.map(|d| DanglingMatcher::new(p, d)),
            |_, h| h,
        );
        new_graph
            .join(
                saturator,
                |_, dl, _, dr| dl.data.matches(dr.data),
                |fl, dl, _, _| (fl, dl),
            )
            .unwrap()
            .map(|_, _, v| v, |_, _, _, _, e| e.map(|d| d.unwrap()), |_, h| h)
    }

    pub fn map_data_ref<'a, E2, V2, H2>(
        &'a self,
        node_map: impl FnMut(&'a Self, N::NeighborsIter<'a>, &'a V) -> V2,
        edge_map: impl FnMut(&'a Self, EdgeIndex, HedgePair, EdgeData<&'a E>) -> EdgeData<E2>,
        mut hedge_map: impl FnMut(Hedge, &'a H) -> H2,
    ) -> HedgeGraph<E2, V2, H2, N::OpStorage<V2>> {
        HedgeGraph {
            hedge_data: self
                .hedge_data
                .iter()
                .map(|(i, a)| hedge_map(i, a))
                .collect(),
            node_store: self.node_store.map_data_ref_graph(self, node_map),
            edge_store: self.edge_store.map_data_ref(self, edge_map),
        }
    }

    pub fn to_ref(&self) -> HedgeGraph<&E, &V, &H, N::OpStorage<&V>> {
        self.map_data_ref(|_, _, v| v, |_, _, _, e| e, |_, h| h)
    }

    pub fn map_data_ref_mut<'a, E2, V2, H2>(
        &'a mut self,
        node_map: impl FnMut(N::NeighborsIter<'a>, &'a mut V) -> V2,
        edge_map: impl FnMut(EdgeIndex, HedgePair, EdgeData<&'a mut E>) -> EdgeData<E2>,
        hedge_map: impl FnMut((Hedge, &'a H)) -> H2,
    ) -> HedgeGraph<E2, V2, H2, N::OpStorage<V2>> {
        HedgeGraph {
            hedge_data: self.hedge_data.iter().map(hedge_map).collect(),
            node_store: self.node_store.map_data_ref_mut_graph(node_map),
            edge_store: self.edge_store.map_data_ref_mut(edge_map),
        }
    }

    pub fn map_data_ref_result<'a, E2, V2, H2, Er>(
        &'a self,
        node_map: impl FnMut(&'a Self, N::NeighborsIter<'a>, &'a V) -> Result<V2, Er>,
        edge_map: impl FnMut(
            &'a Self,
            EdgeIndex,
            HedgePair,
            EdgeData<&'a E>,
        ) -> Result<EdgeData<E2>, Er>,
        hedge_map: impl FnMut((Hedge, &'a H)) -> Result<H2, Er>,
    ) -> Result<HedgeGraph<E2, V2, H2, N::OpStorage<V2>>, Er> {
        let hedge_data = self
            .hedge_data
            .iter()
            .map(hedge_map)
            .collect::<Result<HedgeVec<_>, _>>()?;
        Ok(HedgeGraph {
            hedge_data,
            node_store: self.node_store.map_data_ref_graph_result(self, node_map)?,
            edge_store: self.edge_store.map_data_ref_result(self, edge_map)?,
        })
    }

    pub fn just_structure(&self) -> HedgeGraph<(), (), (), N::OpStorage<()>> {
        self.map_data_ref(|_, _, _| (), |_, _, _, d| d.map(|_| ()), |_, _| ())
    }

    pub fn map_nodes_ref<'a, V2>(
        &'a self,
        f: impl FnMut(&'a Self, N::NeighborsIter<'a>, &'a V) -> V2,
    ) -> HedgeGraph<&'a E, V2, &'a H, N::OpStorage<V2>> {
        HedgeGraph {
            hedge_data: self.hedge_data.iter().collect(),
            node_store: self.node_store.map_data_ref_graph(self, f),
            edge_store: self.edge_store.map_data_ref(self, &|_, _, _, e| e),
        }
    }
    pub fn map<E2, V2, H2>(
        self,
        f: impl FnMut(&Involution, NodeIndex, V) -> V2,
        g: impl FnMut(&Involution, &N, HedgePair, EdgeIndex, EdgeData<E>) -> EdgeData<E2>,
        mut h: impl FnMut(Hedge, H) -> H2,
    ) -> HedgeGraph<E2, V2, H2, N::OpStorage<V2>> {
        let edge_store = self.edge_store.map_data(&self.node_store, g);
        HedgeGraph {
            hedge_data: self
                .hedge_data
                .into_iter()
                .map(|(heg, i)| h(heg, i))
                .collect(),
            node_store: self.node_store.map_data_graph(edge_store.as_ref(), f),
            edge_store,
        }
    }
    pub fn map_result<E2, V2, H2, Err>(
        self,
        f: impl FnMut(&Involution, NodeIndex, V) -> Result<V2, Err>,
        g: impl FnMut(&Involution, &N, HedgePair, EdgeIndex, EdgeData<E>) -> Result<EdgeData<E2>, Err>,
        mut h: impl FnMut(Hedge, H) -> Result<H2, Err>,
    ) -> Result<HedgeGraph<E2, V2, H2, N::OpStorage<V2>>, Err> {
        let edge_store = self.edge_store.map_data_result(&self.node_store, g)?;
        Ok(HedgeGraph {
            hedge_data: self
                .hedge_data
                .into_iter()
                .map(|(heg, i)| h(heg, i))
                .collect::<Result<HedgeVec<_>, Err>>()?,
            node_store: self
                .node_store
                .map_data_graph_result(edge_store.as_ref(), f)?,
            edge_store,
        })
    }
    pub fn new_smart_hedgevec<T>(
        &self,
        f: &impl Fn(HedgePair, EdgeData<&E>) -> EdgeData<T>,
    ) -> SmartEdgeVec<T> {
        self.edge_store.map(f)
    }

    pub fn new_edgevec<T>(&self, f: impl FnMut(&E, EdgeIndex, &HedgePair) -> T) -> EdgeVec<T> {
        self.edge_store.new_edgevec(f)
    }

    pub fn new_nodevec<'a, T>(
        &'a self,
        f: impl FnMut(NodeIndex, N::NeighborsIter<'a>, &'a V) -> T,
    ) -> NodeVec<T> {
        self.node_store.new_nodevec(f)
    }

    pub fn new_hedgevec<T>(&self, mut f: impl FnMut(Hedge, &H) -> T) -> HedgeVec<T> {
        self.hedge_data.iter().map(|(i, h)| f(i, h)).collect()
    }

    pub fn new_edgevec_from_iter<T, I: IntoIterator<Item = T>>(
        &self,
        iter: I,
    ) -> Result<EdgeVec<T>, HedgeGraphError> {
        self.edge_store.new_edgevec_from_iter(iter)
    }
}

// Cuts
impl<E, V, H, N: NodeStorageOps<NodeData = V>> HedgeGraph<E, V, H, N> {
    fn non_cut_edges_impl(
        &self,
        connected_components: usize,
        cyclotomatic_number: usize,
        start: Hedge,
        current: &mut SuBitGraph,
        set: &mut AHashSet<SuBitGraph>,
    ) {
        if current.n_included() > 2 * cyclotomatic_number {
            return;
        }

        let complement = current.complement(self);

        if current.n_included() > 0
            && self.count_connected_components(&complement) == connected_components
            && complement.covers::<_, _, _, _, SuBitGraph>(self) == self.full_filter()
        {
            // println!("//inserted with {con_comp}");
            set.insert(current.clone());
        }

        for i in (start.0..self.n_hedges()).map(Hedge) {
            let j = self.inv(i);
            if i > j {
                current.add(i);
                current.add(j);
                self.non_cut_edges_impl(
                    connected_components,
                    cyclotomatic_number,
                    Hedge(i.0 + 1),
                    current,
                    set,
                );
                current.sub(i);
                current.sub(j);
            }
        }
    }

    /// all sets of full edges that do not disconnect the graph/ increase its connected components
    pub fn non_cut_edges(&self) -> AHashSet<SuBitGraph> {
        let connected_components = self.count_connected_components(&self.full_filter());

        let cyclotomatic_number = self.cyclotomatic_number(&self.full_node().internal_graph);

        let mut current = self.empty_subgraph::<SuBitGraph>();
        let mut set = AHashSet::new();

        self.non_cut_edges_impl(
            connected_components,
            cyclotomatic_number,
            Hedge(0),
            &mut current,
            &mut set,
        );

        set
    }

    // pub fn backtracking_cut_set(
    //     &self,
    //     source: HedgeNode,
    //     target: HedgeNode,
    // ) -> Vec<InternalSubGraph> {
    // }

    pub fn non_bridges(&self) -> SuBitGraph {
        let (c, _) = self.cycle_basis();
        let mut cycle_cover: SuBitGraph = self.empty_subgraph();
        for cycle in c {
            cycle_cover.union_with(&cycle.filter);
        }

        cycle_cover
    }

    pub fn non_bridges_of<S: SubSetLike<Base = SuBitGraph> + SubGraphLike>(
        &self,
        subgraph: &S,
    ) -> SuBitGraph {
        let (c, _) = self.cycle_basis_of(subgraph);
        let mut cycle_cover: SuBitGraph = self.empty_subgraph();
        for cycle in c {
            cycle_cover.union_with(&cycle.filter);
        }

        cycle_cover
    }

    pub fn bridges_of<S: SubSetLike<Base = SuBitGraph> + SubGraphLike>(
        &self,
        subgraph: &S,
    ) -> SuBitGraph {
        subgraph.included().subtract(&self.non_bridges_of(subgraph))
    }

    pub fn bridges(&self) -> SuBitGraph {
        self.non_bridges().complement(self)
    }

    pub fn combine_to_single_hedgenode(&self, source: &[NodeIndex]) -> HedgeNode {
        let s: SuBitGraph =
            source
                .iter()
                .map(|a| self.iter_crown(*a))
                .fold(self.empty_subgraph(), |mut acc, e| {
                    acc.union_with_iter(e);
                    acc
                });

        let (hairs, internal_graph) = self.split_hairs_and_internal_hedges(s);

        HedgeNode {
            internal_graph,
            hairs,
        }
    }

    pub fn all_cuts_from_ids(
        &self,
        source: &[NodeIndex],
        target: &[NodeIndex],
    ) -> Vec<(SuBitGraph, OrientedCut, SuBitGraph)>
    where
        N: NodeStorageOps,
    {
        let source = self.combine_to_single_hedgenode(source);
        let target = self.combine_to_single_hedgenode(target);
        self.all_cuts(source, target)
    }

    pub fn tadpoles(&self, externals: &[NodeIndex]) -> Vec<SuBitGraph> {
        let mut identified: HedgeGraph<(), (), (), N::OpStorage<()>> = self.just_structure();

        let n = identified.identify_nodes(externals, ());
        let hairs = identified.iter_crown(n).next().unwrap();

        let non_bridges = identified.non_bridges();
        identified
            .connected_components(&non_bridges)
            .into_iter()
            .filter_map(|mut a| {
                if !a.includes(&hairs) {
                    let full: SuBitGraph = a.covers(self);

                    for i in full.included_iter() {
                        a.add(i);
                        a.add(self.inv(i));
                    }
                    Some(a)
                } else {
                    None
                }
            })
            .collect()
    }

    pub fn all_bonds<R: RangeBounds<usize>>(&self, size: &R) -> Vec<SuBitGraph> {
        self.all_bonds_of(&self.full_filter(), size)
    }
    pub fn a_bond(&self, cond: &impl Fn(&SuBitGraph) -> bool) -> Option<SuBitGraph> {
        self.a_bond_of(&self.full_filter(), cond)
    }

    fn all_bonds_of_impl<R: RangeBounds<usize>>(
        &self,
        subgraph: &SuBitGraph,
        size: &R,
    ) -> Vec<SuBitGraph> {
        let Some((growth, cut)) = self
            .iter_nodes_of(subgraph)
            .map(|(_, crown, _)| {
                let mut subcrown: SuBitGraph = self.empty_subgraph();
                for i in crown.filter(|i| subgraph.includes(i)) {
                    subcrown.add(i);
                }
                let cut: SuBitGraph = self.internal_crown(&subcrown);
                (subcrown, cut)
            })
            .next()
        else {
            return vec![];
        };

        let mut regions = BTreeSet::new();

        self.all_bonds_backtrack_of(&mut regions, subgraph, cut, growth, size)
    }
    fn a_bond_of_impl(
        &self,
        subgraph: &SuBitGraph,
        cond: &impl Fn(&SuBitGraph) -> bool,
    ) -> Option<SuBitGraph> {
        let (growth, cut) = self
            .iter_nodes_of(subgraph)
            .map(|(_, crown, _)| {
                let mut subcrown: SuBitGraph = self.empty_subgraph();
                for i in crown.filter(|i| subgraph.includes(i)) {
                    subcrown.add(i);
                }
                let cut: SuBitGraph = self.internal_crown(&subcrown);
                (subcrown, cut)
            })
            .next()?;

        let mut regions = BTreeSet::new();

        self.a_bond_backtrack_of(&mut regions, subgraph, cut, growth, cond)
    }

    fn a_bond_backtrack_of(
        &self,
        regions: &mut BTreeSet<SuBitGraph>,
        subgraph: &SuBitGraph,
        cut: SuBitGraph,
        growth: SuBitGraph,
        cond: &impl Fn(&SuBitGraph) -> bool,
    ) -> Option<SuBitGraph> {
        if regions.contains(&growth) {
            return None;
        } else {
            regions.insert(growth.clone());
        }

        if cond(&cut) && self.is_connected(&subgraph.subtract(&growth)) {
            return Some(cut);
        }

        // println!(
        //     "//cut: \n{}, //growth: \n{}",
        //     self.dot(&cut),
        //     self.dot(&growth)
        // );
        let mut seen = BTreeSet::new();
        for h in cut.included_iter() {
            let invh = self.inv(h);
            if invh == h || growth.includes(&invh) {
                continue;
            }
            let new_node = self.node_id(invh);

            if !seen.insert(new_node) {
                continue;
            }

            let mut new_growth = growth.clone();
            for c in self.iter_crown(new_node) {
                new_growth.add(c);
            }
            let new_cut = self.internal_crown(&new_growth);
            let Some(cut) = self.a_bond_backtrack_of(regions, subgraph, new_cut, new_growth, cond)
            else {
                continue;
            };

            return Some(cut);
        }

        None
    }

    fn all_bonds_backtrack_of<R: RangeBounds<usize>>(
        &self,
        regions: &mut BTreeSet<SuBitGraph>,
        subgraph: &SuBitGraph,
        cut: SuBitGraph,
        growth: SuBitGraph,
        size: &R,
    ) -> Vec<SuBitGraph> {
        let mut cuts = vec![];
        if regions.contains(&growth) {
            return cuts;
        } else {
            regions.insert(growth.clone());
        }

        // println!(
        //     "//cut: \n{}, //growth: \n{}",
        //     self.dot(&cut),
        //     self.dot(&growth)
        // );
        let mut seen = BTreeSet::new();
        for h in cut.included_iter() {
            let invh = self.inv(h);
            if invh == h || growth.includes(&invh) {
                continue;
            }
            let new_node = self.node_id(invh);

            if !seen.insert(new_node) {
                continue;
            }

            let mut new_growth = growth.clone();
            for c in self.iter_crown(new_node) {
                new_growth.add(c);
            }
            let new_cut = self.internal_crown(&new_growth);
            cuts.extend(self.all_bonds_backtrack_of(regions, subgraph, new_cut, new_growth, size));
        }

        // println!("Hi");
        if size.contains(&cut.n_included()) && self.is_connected(&subgraph.subtract(&growth)) {
            cuts.push(cut);
        }
        cuts
    }

    pub fn all_bonds_of<S: SubGraphLike<Base = SuBitGraph>, R: RangeBounds<usize>>(
        &self,
        subgraph: &S,
        size: &R,
    ) -> Vec<SuBitGraph> {
        let mut cuts = vec![];
        for c in self.connected_components(subgraph) {
            cuts.extend(self.all_bonds_of_impl(c.included(), size))
        }
        cuts
    }
    pub fn a_bond_of<S: SubGraphLike<Base = SuBitGraph>>(
        &self,
        subgraph: &S,
        cond: &impl Fn(&SuBitGraph) -> bool,
    ) -> Option<SuBitGraph> {
        for c in self.connected_components(subgraph) {
            let Some(cut) = self.a_bond_of_impl(c.included(), cond) else {
                continue;
            };
            return Some(cut);
        }
        None
    }

    pub fn all_cuts(
        &self,
        source: HedgeNode,
        target: HedgeNode,
    ) -> Vec<(SuBitGraph, OrientedCut, SuBitGraph)>
    where
        N: NodeStorageOps,
    {
        // println!("//Source\n{}", self.dot(&source.hairs));
        // println!("//Target\n{}", self.dot(&target.hairs));

        let full_source = source.internal_and_hairs();
        let full_target = target.internal_and_hairs();
        let s_connectivity = self.count_connected_components(&full_source);

        let t_connectivity = self.count_connected_components(&full_target);

        let augmented: HedgeGraph<(), (), (), N::OpStorage<()>> = self.just_structure();
        let s_nodes = self
            .iter_nodes_of(&source)
            .map(|a| self.id_from_crown(a.1).unwrap())
            .collect::<Vec<_>>();
        let t_nodes = self
            .iter_nodes_of(&target)
            .map(|a| self.id_from_crown(a.1).unwrap())
            .collect::<Vec<_>>();

        let t_node = t_nodes[0];
        let s_node = s_nodes[0];

        let augmented = s_nodes.iter().fold(augmented, |aug, n| {
            let (_, _, augmented) = aug.add_pair(*n, t_node, (), false).unwrap();
            augmented
        });

        let augmented = t_nodes.iter().fold(augmented, |aug, n| {
            let (_, _, augmented) = aug.add_pair(s_node, *n, (), false).unwrap();
            augmented
        });

        let mut non_bridges = augmented.non_bridges();

        for _ in &s_nodes {
            non_bridges.pop();
            non_bridges.pop();
        }
        for _ in &t_nodes {
            non_bridges.pop();
            non_bridges.pop();
        }

        // println!("//non_bridges:\n{}", self.dot(&non_bridges));

        non_bridges.union_with(&source.hairs);
        non_bridges.union_with(&target.hairs);

        // println!("//non_bridges:\n{}", self.dot(&non_bridges));

        let mut regions = AHashSet::new();
        self.all_s_t_cuts_impl(
            &non_bridges,
            s_connectivity,
            source,
            &target,
            t_connectivity,
            &mut regions,
        );

        let mut cuts = vec![];
        let bridges = non_bridges.complement(self);

        for mut r in regions.drain() {
            // let disconnected = r.complement(self);

            // let mut s_side_covers = if s.hairs.intersects(&r) {
            //     s.hairs.clone()
            // } else {
            //     SimpleTraversalTree::depth_first_traverse(self, &disconnected, &source, None)
            //         .unwrap()
            //         .covers()
            // };
            let cut = OrientedCut::from_underlying_coerce(r.hairs.clone(), self).unwrap();
            r.add_all_hairs(self);
            let mut s_side_covers = r.internal_and_hairs();
            for i in bridges.included_iter() {
                if s_side_covers.includes(&self.inv(i)) {
                    s_side_covers.add(i);
                }
            }

            let t_side_covers = s_side_covers.complement(self);

            // let internal = InternalSubGraph::cleaned_filter_pessimist(t_side_covers, self);
            // let mut t_side = self.nesting_node_from_subgraph(internal);
            // t_side.hairs.union_with(&t.hairs);
            //

            cuts.push((s_side_covers, cut, t_side_covers));
        }

        cuts
    }

    pub fn all_s_t_cuts_impl<S: SubSetLike<Base = SuBitGraph>>(
        &self,
        subgraph: &S,
        s_connectivity: usize,
        s: HedgeNode, // will grow
        t: &HedgeNode,
        t_connectivity: usize,
        regions: &mut AHashSet<HedgeNode>,
    ) {
        // println!("regions size:{}", regions.len());
        //

        let hairy = s.internal_graph.filter.union(&s.hairs);
        let mut complement = hairy.complement(self);
        complement.intersect_with(subgraph.included());
        if !complement.includes(&t.hairs) {
            return;
        }

        let t_count = self.count_connected_components(&complement);
        let s_count = self.count_connected_components(&hairy);

        if t_count <= t_connectivity && s_count <= s_connectivity && regions.get(&s).is_none() {
            for h in s.hairs.included_iter() {
                let invh = self.inv(h);

                if invh != h && !t.hairs.includes(&invh) && subgraph.includes(&invh) {
                    let mut new_node = s.clone();

                    new_node.hairs.union_with_iter(self.neighbors(invh));
                    new_node.hairs.intersect_with(subgraph.included());

                    new_node.fix(self);
                    self.all_s_t_cuts_impl(
                        subgraph,
                        s_connectivity,
                        new_node,
                        t,
                        t_connectivity,
                        regions,
                    );
                }
            }
            regions.insert(s);
        }
    }
}

// Cycles
impl<E, V, H, N: NodeStorageOps<NodeData = V>> HedgeGraph<E, V, H, N> {
    ///Gives all subgraphs corresponding to all the spanning trees of the graph
    ///Winter, Pawel. “An Algorithm for the Enumeration of Spanning Trees.” BIT Numerical Mathematics 26, no. 1 (March 1, 1986): 44–62. https://doi.org/10.1007/BF01939361.
    pub fn all_spanning_forests_of<S: SubGraphLike>(&self, subgraph: &S) -> Vec<S::Base>
    where
        for<'a> N::OpStorage<&'a V>: Clone,
        S::Base: SubSetLike<Base = S::Base>
            + SubSetOps
            + Clone
            + ModifySubSet<HedgePair>
            + ModifySubSet<Hedge>,
    {
        let ref_self = self.to_ref();

        let exts = self.internal_crown(subgraph);

        if subgraph.is_empty() {
            return vec![];
        }
        let mut visited_edges: SuBitGraph = self.empty_subgraph();

        let mut nodes = vec![];

        // Iterate over all edges in the subgraph
        for hedge_index in subgraph.included_iter() {
            if !visited_edges.includes(&hedge_index) {
                // Perform DFS to find all reachable edges from this edge

                //
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

        let mut included = subgraph.included().clone();
        for h in exts.included_iter() {
            included.sub(h);
        }

        ref_self.contract_edge_for_spanning(nodes, included)
    }

    fn contract_edge_for_spanning<S>(self, mut nodes: Vec<NodeIndex>, mut subgraph: S) -> Vec<S>
    where
        V: Clone,
        E: Clone,
        H: Clone,
        N: Clone,
        S: SubSetOps + Clone + SubSetLike<Base = S> + ModifySubSet<HedgePair> + ModifySubSet<Hedge>,
    {
        let mut trees = vec![];
        if let Some(node) = nodes.pop() {
            for h in self.iter_crown(node) {
                //This is an internal edge
                if subgraph.includes(&h) && subgraph.includes(&self.inv(h)) {
                    //we found an neighboring vertex in the graph
                    if let Some(new_node) = self.involved_node_id(h) {
                        // This is not a dangling edge
                        // if new_node != node {

                        let node_data_merge = self[node].clone();
                        let mut new_self = self.clone();

                        // We now identify node and new_node. All spanning trees of this contracted graph will be spanning trees of the full graph, when adding one of the parallel edges
                        let (mapped_node, parallel): (_, S::Base) = new_self
                            .identify_nodes_without_self_edges(&[node, new_node], node_data_merge);

                        let mapped_nodes = nodes
                            .iter()
                            .map(|a| {
                                if *a == new_node || *a == node {
                                    mapped_node
                                } else {
                                    *a
                                }
                            })
                            .collect();

                        for v in new_self
                            .contract_edge_for_spanning(mapped_nodes, subgraph.subtract(&parallel))
                        // recurse with a the node identified graph and the contracted edge subgraph
                        {
                            // v is the set of spanning trees of the contracted graph
                            for (p, _, _) in self.iter_edges_of(&parallel.intersection(&subgraph)) {
                                let mut with_p = v.clone();
                                with_p.add(p);
                                trees.push(with_p);
                            }
                        }

                        // Remove all these edges from consideration.
                        subgraph.subtract_with(&parallel);
                    }
                }
            }
        } else {
            trees.push(self.empty_subgraph())
        }

        trees
    }
}

// Cycles
impl<E, V, H, N: NodeStorageOps<NodeData = V>> HedgeGraph<E, V, H, N> {
    pub fn cyclotomatic_number<S: SubGraphLike>(&self, subgraph: &S) -> usize {
        let n_hedges = self.count_internal_edges(subgraph);
        // println!("n_hedges: {}", n_hedges);
        let n_nodes = self.number_of_nodes_in_subgraph(subgraph);
        // println!("n_nodes: {}", n_nodes);
        let n_components = self.count_connected_components(subgraph);

        n_hedges + n_components - n_nodes
    }

    pub fn cycle_basis(&self) -> (Vec<Cycle>, SuBitGraph) {
        self.paton_cycle_basis(&self.full_filter()).unwrap()
    }

    pub fn cycle_basis_of<S: SubSetLike<Base = SuBitGraph> + SubGraphLike>(
        &self,
        subgraph: &S,
    ) -> (Vec<Cycle>, SuBitGraph) {
        self.paton_cycle_basis(subgraph.included()).unwrap()
    }

    pub fn order_basis(&self, basis: &[HedgeNode]) -> Vec<Vec<InternalSubGraph>> {
        let mut seen = vec![basis[0].internal_graph.clone()];
        let mut partitions = vec![seen.clone()];

        for cycle in basis.iter() {
            if seen
                .iter()
                .any(|p| !p.empty_intersection(&cycle.internal_graph))
            {
                partitions
                    .last_mut()
                    .unwrap()
                    .push(cycle.internal_graph.clone());
            } else {
                for p in partitions.last().unwrap() {
                    seen.push(p.clone());
                }
                partitions.push(vec![cycle.internal_graph.clone()]);
            }
        }

        partitions
    }

    pub fn all_cycles(&self) -> Vec<Cycle> {
        Cycle::all_sum_powerset_filter_map(&self.cycle_basis().0, &|mut c| {
            if c.is_circuit(self) {
                c.loop_count = Some(1);
                Some(c)
            } else {
                None
            }
        })
        .unwrap()
        .into_iter()
        .collect()
    }

    pub fn all_cycle_sym_diffs(&self) -> Result<Vec<InternalSubGraph>, TryFromIntError> {
        Cycle::all_sum_powerset_filter_map(&self.cycle_basis().0, &Some)
            .map(|a| a.into_iter().map(|c| c.internal_graph(self)).collect())
    }

    pub fn all_cycle_unions(&self) -> AHashSet<InternalSubGraph> {
        InternalSubGraph::all_unions_iterative(&self.all_cycle_sym_diffs().unwrap())
    }
    pub fn paton_count_loops(
        &self,
        subgraph: &InternalSubGraph,
        start: &NodeIndex,
    ) -> Result<usize, HedgeGraphError> {
        let tree = SimpleTraversalTree::depth_first_traverse(self, subgraph, start, None)?;

        let cuts = subgraph.subtract(&tree.tree_subgraph(self));
        Ok(self.edge_store.n_internals(&cuts))
    }

    fn paton_cycle_basis(
        &self,
        subgraph: &SuBitGraph,
    ) -> Result<(Vec<Cycle>, SuBitGraph), HedgeGraphError> {
        let mut visited_edges: SuBitGraph = self.empty_subgraph();

        let mut cycle_basis = vec![];
        let mut forest: SuBitGraph = self.empty_subgraph();

        // Iterate over all edges in the subgraph
        for hedge_index in subgraph.included_iter() {
            if visited_edges.includes(&hedge_index) {
                continue;
            }

            let root_node = self.node_id(hedge_index);
            let tree = SimpleTraversalTree::depth_first_traverse(self, subgraph, &root_node, None)?;

            let reachable_edges = tree.covers(subgraph);
            let cuts = subgraph.included().subtract(&tree.tree_subgraph);

            visited_edges.union_with(&reachable_edges);
            for c in cuts.included_iter() {
                if c > self.inv(c) && subgraph.includes(&self.inv(c)) {
                    cycle_basis.push(tree.get_cycle(c, self).unwrap());
                }
            }
            forest.union_with(&tree.tree_subgraph);
        }

        Ok((cycle_basis, forest))
    }

    pub fn all_spinneys_with_basis(&self, basis: &[&InternalSubGraph]) -> AHashSet<HedgeNode> {
        let mut spinneys = AHashSet::new();
        let mut base_cycle: InternalSubGraph = self.empty_subgraph();

        for cycle in basis {
            base_cycle.sym_diff_with(cycle);
        }

        spinneys.insert(self.nesting_node_from_subgraph(base_cycle.clone()));

        if basis.len() == 1 {
            return spinneys;
        }

        for i in 0..basis.len() {
            for s in self.all_spinneys_with_basis(
                &basis
                    .iter()
                    .enumerate()
                    .filter_map(|(j, s)| if j != i { Some(*s) } else { None })
                    .collect_vec(),
            ) {
                spinneys
                    .insert(self.nesting_node_from_subgraph(s.internal_graph.union(&base_cycle)));
                spinneys.insert(s);
            }
        }

        spinneys
    }

    pub fn all_spinneys_rec(&self, spinneys: &mut AHashSet<HedgeNode>, cycle_sums: Vec<HedgeNode>) {
        let _len = spinneys.len();

        // let mut pset:Su = PowersetIterator::new(cycle_sums.len() as u8);

        // pset.next(); //Skip empty set

        for (ci, cj) in cycle_sums.iter().tuple_combinations() {
            let _union = ci.internal_graph.union(&cj.internal_graph);

            // spinneys.insert(union);
        }
    }

    pub fn all_spinneys(
        &self,
    ) -> AHashMap<InternalSubGraph, Vec<(InternalSubGraph, Option<InternalSubGraph>)>> {
        let basis_cycles = self.cycle_basis().0;

        let mut all_combinations = PowersetIterator::<usize>::new(basis_cycles.len() as u8);
        all_combinations.next(); //Skip empty set

        let mut spinneys: AHashMap<
            InternalSubGraph,
            Vec<(InternalSubGraph, Option<InternalSubGraph>)>,
        > = AHashMap::new();

        let mut cycles: Vec<InternalSubGraph> = Vec::new();
        for p in all_combinations {
            let mut base_cycle: InternalSubGraph = self.empty_subgraph();

            for i in p.included_iter() {
                base_cycle.sym_diff_with(&basis_cycles[i].clone().internal_graph(self));
            }

            cycles.push(base_cycle);
        }

        for (ci, cj) in cycles.iter().tuple_combinations() {
            let union = ci.union(cj);

            if let Some(v) = spinneys.get_mut(&union) {
                v.push((ci.clone(), Some(cj.clone())));
            } else {
                spinneys.insert(union, vec![(ci.clone(), Some(cj.clone()))]);
            }
        }

        for c in cycles {
            spinneys.insert(c.clone(), vec![(c.clone(), None)]);
        }
        spinneys
    }

    pub fn all_spinneys_alt(&self) -> AHashSet<InternalSubGraph> {
        let mut spinneys = AHashSet::new();
        let cycles = self.all_cycles();

        let mut pset = PowersetIterator::<usize>::new(cycles.len() as u8);
        pset.next(); //Skip empty set

        for p in pset {
            let mut union: InternalSubGraph = self.empty_subgraph();

            for i in p.included_iter() {
                union.union_with(&cycles[i].clone().internal_graph(self));
            }

            spinneys.insert(union);
        }

        for c in cycles {
            spinneys.insert(c.internal_graph(self));
        }
        spinneys
    }
}

// Traversal Trees
impl<E, V, H, N: NodeStorageOps<NodeData = V>> HedgeGraph<E, V, H, N> {
    pub fn count_connected_components<S: SubGraphLike>(&self, subgraph: &S) -> usize {
        self.connected_components(subgraph).len()
    }

    pub fn connected_components<S: SubGraphLike>(&self, subgraph: &S) -> Vec<SuBitGraph> {
        let mut visited_edges: SuBitGraph = self.empty_subgraph();

        let mut components = vec![];

        // Iterate over all edges in the subgraph
        for hedge_index in subgraph.included_iter() {
            if !visited_edges.includes(&hedge_index) {
                // Perform DFS to find all reachable edges from this edge

                //
                let root_node = self.node_id(hedge_index);
                let reachable_edges =
                    SimpleTraversalTree::depth_first_traverse(self, subgraph, &root_node, None)
                        .unwrap()
                        .covers(subgraph);

                visited_edges.union_with(&reachable_edges);

                components.push(reachable_edges);
            }
        }
        components
    }
    pub fn align_underlying_to_superficial(&mut self) {
        self.edge_store.align_underlying_to_superficial();
    }

    /// aligns the underlying orientation of the graph to the tree, such that all tree edges are oriented towards the root, and all others point towards the leaves
    ///
    /// Relies on the tree being tremaux (i.e. the tree-order is total)
    /// This is the case for depth-first traversal.
    pub fn align_underlying_to_tree<P: ForestNodeStore<NodeData = ()>>(
        &mut self,
        tree: &SimpleTraversalTree<P>,
    ) {
        for (h, tt, i) in tree.iter_hedges() {
            // println!("hedge: {}, tt: {:?}, i: {:?}\n", h, tt, i);
            match tt {
                tree::TTRoot::Root => {
                    if i.is_some() {
                        self.edge_store.set_flow(h, Flow::Sink);
                    } else {
                        self.edge_store.set_flow(h, Flow::Source);
                    }
                }
                tree::TTRoot::Child(_) => {
                    if let Some(root_pointer) = i {
                        if root_pointer == h {
                            self.edge_store.set_flow(h, Flow::Source);
                        } else {
                            self.edge_store.set_flow(h, Flow::Sink);
                        }
                    } else {
                        let current_node_id = tree.node_id(h);
                        let involved_node_id = tree.node_id(self.inv(h));

                        let order =
                            tree.tree_order(current_node_id, involved_node_id, &self.edge_store);
                        match order {
                            Some(Ordering::Equal) => {
                                // self.edge_store.set_flow(h);
                            }
                            Some(Ordering::Less) => {
                                //the path to the root from the current node, passes through the involved node
                                self.edge_store.set_flow(h, Flow::Sink);
                            }
                            Some(Ordering::Greater) => {
                                self.edge_store.set_flow(h, Flow::Source);
                            }
                            None => {}
                        }
                    }
                }
                tree::TTRoot::None => {}
            }
        }
    }

    /// aligns the superficial orientation of the graph to the tree,
    ///
    /// such that all tree edges are oriented towards the root, and all others are unoriented.
    pub fn align_superficial_to_tree<P: ForestNodeStore<NodeData = ()>>(
        &mut self,
        tree: &SimpleTraversalTree<P>,
    ) {
        for (h, tt, i) in tree.iter_hedges() {
            match tt {
                tree::TTRoot::Root => {
                    if i.is_some() {
                        let flow = self.edge_store.flow(h);
                        match flow {
                            Flow::Source => {
                                self.edge_store.set_orientation(h, Orientation::Reversed);
                            }
                            Flow::Sink => {
                                self.edge_store.set_orientation(h, Orientation::Default);
                            }
                        }
                    } else {
                        self.edge_store.set_orientation(h, Orientation::Undirected);
                    }
                }
                tree::TTRoot::Child(_) => {
                    if let Some(root_pointer) = i {
                        let flow = self.edge_store.flow(h);
                        if root_pointer == h {
                            match flow {
                                Flow::Source => {
                                    self.edge_store.set_orientation(h, Orientation::Default);
                                }
                                Flow::Sink => {
                                    self.edge_store.set_orientation(h, Orientation::Reversed);
                                }
                            }
                        } else {
                            match flow {
                                Flow::Source => {
                                    self.edge_store.set_orientation(h, Orientation::Reversed);
                                }
                                Flow::Sink => {
                                    self.edge_store.set_orientation(h, Orientation::Default);
                                }
                            }

                            // self.edge_store.involution.set_as_sink(h);
                        }
                    } else {
                        self.edge_store.set_orientation(h, Orientation::Undirected);
                    }
                }
                tree::TTRoot::None => {}
            }
        }
    }

    // pub fn align_to_tree_superficial(&mut self, tree: &TraversalTree) {
    //     for (i, p) in tree.parent_iter() {
    //         match self.edge_store.involution.hedge_data_mut(i) {
    //             InvolutiveMapping::Identity { data, .. } => match p {
    //                 Parent::Root => {}
    //                 Parent::Hedge { hedge_to_root, .. } => {
    //                     if *hedge_to_root == i {
    //                         data.orientation = Orientation::Default;
    //                     } else {
    //                         data.orientation = Orientation::Reversed;
    //                     }
    //                 }
    //                 Parent::Unset => {}
    //             },
    //             InvolutiveMapping::Source { data, .. } => match p {
    //                 Parent::Root => {}
    //                 Parent::Hedge { hedge_to_root, .. } => {
    //                     if *hedge_to_root == i {
    //                         data.orientation = Orientation::Default;
    //                     } else {
    //                         data.orientation = Orientation::Reversed;
    //                     }
    //                 }
    //                 Parent::Unset => {}
    //             },
    //             _ => {}
    //         }
    //     }
    // }
}

// Iterators
impl<E, V, H, N: NodeStorageOps<NodeData = V>> HedgeGraph<E, V, H, N> {
    pub fn iter_hedges(&self) -> impl Iterator<Item = (Hedge, &H)> {
        self.hedge_data.iter()
    }

    ///Iterate over all nodes, returns an iterator that yields
    pub fn iter_nodes(&self) -> impl Iterator<Item = (NodeIndex, N::NeighborsIter<'_>, &V)> {
        self.node_store.iter_nodes()
    }

    pub fn iter_nodes_mut(
        &mut self,
    ) -> impl Iterator<Item = (NodeIndex, N::NeighborsIter<'_>, &mut V)> {
        self.node_store.iter_nodes_mut()
    }

    pub fn iter_node_ids(&self) -> impl Iterator<Item = NodeIndex> + '_ {
        self.node_store.iter_node_id()
    }

    pub fn iter_edge_ids_of<'a, S: SubSetLike>(
        &'a self,
        subgraph: &'a S,
    ) -> EdgeIter<'a, E, V, H, S, N, S::BaseIter<'a>> {
        EdgeIter::new(self, subgraph)
    }

    pub fn iter_edges_of<'a, S: SubSetLike>(
        &'a self,
        subgraph: &'a S,
    ) -> impl Iterator<Item = (HedgePair, EdgeIndex, EdgeData<&'a E>)> + 'a {
        self.edge_store.iter_edges_of(subgraph)
    }

    pub fn iter_edges(&self) -> impl Iterator<Item = (HedgePair, EdgeIndex, EdgeData<&E>)> {
        self.edge_store.iter_edges()
    }

    pub fn iter_nodes_of<'a, S: SubSetLike>(
        &'a self,
        subgraph: &'a S,
    ) -> impl Iterator<Item = (NodeIndex, N::NeighborsIter<'a>, &'a V)>
    where
        N::NeighborsIter<'a>: Clone,
    {
        NodeIterator {
            graph: self,
            edges: subgraph.included_iter(),
            seen: SubSet::empty(self.n_nodes()),
        }
    }
}

// Display
impl<E, V, H, N: NodeStorageOps<NodeData = V>> HedgeGraph<E, V, H, N> {
    pub fn dot_impl_fmt<S: SubGraphLike, Str1: AsRef<str>>(
        &self,
        writer: &mut impl std::fmt::Write,
        subgraph: &S,
        graph_info: Str1,
        hedge_attr: &impl Fn(&H) -> Option<String>,
        edge_attr: &impl Fn(&E) -> Option<String>,
        node_attr: &impl Fn(&V) -> Option<String>,
    ) -> Result<(), std::fmt::Error> {
        subgraph.dot_fmt(writer, self, graph_info, hedge_attr, edge_attr, node_attr)
    }
    pub fn dot_impl_io<S: SubGraphLike, Str1: AsRef<str>>(
        &self,
        writer: &mut impl std::io::Write,
        subgraph: &S,
        graph_info: Str1,
        hedge_attr: &impl Fn(&H) -> Option<String>,
        edge_attr: &impl Fn(&E) -> Option<String>,
        node_attr: &impl Fn(&V) -> Option<String>,
    ) -> Result<(), std::io::Error> {
        subgraph.dot_io(writer, self, graph_info, hedge_attr, edge_attr, node_attr)
    }

    pub fn dot_impl<S: SubGraphLike, Str1: AsRef<str>>(
        &self,
        subgraph: &S,
        graph_info: Str1,
        hedge_attr: &impl Fn(&H) -> Option<String>,
        edge_attr: &impl Fn(&E) -> Option<String>,
        node_attr: &impl Fn(&V) -> Option<String>,
    ) -> String {
        let mut output = String::new();
        subgraph
            .dot_fmt(
                &mut output,
                self,
                graph_info,
                hedge_attr,
                edge_attr,
                node_attr,
            )
            .unwrap();
        output
    }

    pub fn dot<S: SubGraphLike>(&self, node_as_graph: &S) -> String {
        let mut output = String::new();
        self.dot_impl_fmt(
            &mut output,
            node_as_graph,
            "start=2;\n",
            &|_| None,
            &|_| None,
            &|_| None,
        )
        .unwrap();
        output
    }

    pub fn dot_display<S: SubGraphLike>(&self, node_as_graph: &S) -> String
    where
        E: Display,
        V: Display,
        H: Display,
    {
        let mut output = String::new();
        self.dot_impl_fmt(
            &mut output,
            node_as_graph,
            "start=2;\n",
            &|h| Some(format!("{h}")),
            &|a| Some(format!("{a}")),
            &|v| Some(format!("{v}")),
        )
        .unwrap();
        output
    }

    pub fn dot_label<S: SubGraphLike>(&self, node_as_graph: &S) -> String
    where
        E: Display,
        V: Display,
    {
        let mut output = String::new();
        self.dot_impl_fmt(
            &mut output,
            node_as_graph,
            "start=2;\n",
            &|_| None,
            &|a| Some(format!("label=\"{a}\"")),
            &|v| Some(format!("label=\"{v}\"")),
        )
        .unwrap();
        output
    }

    pub fn base_dot(&self) -> String {
        self.dot(&self.full_filter())
    }
}

impl<E, V, H, N: NodeStorage<NodeData = V>> Index<Hedge> for HedgeGraph<E, V, H, N> {
    type Output = H;
    fn index(&self, index: Hedge) -> &Self::Output {
        &self.hedge_data[index]
    }
}

impl<E, V, H, N: NodeStorage<NodeData = V>> IndexMut<Hedge> for HedgeGraph<E, V, H, N> {
    fn index_mut(&mut self, index: Hedge) -> &mut Self::Output {
        &mut self.hedge_data[index]
    }
}

impl<E, V, H, N: NodeStorage<NodeData = V>> Index<&Hedge> for HedgeGraph<E, V, H, N> {
    type Output = EdgeIndex;
    fn index(&self, index: &Hedge) -> &Self::Output {
        &self.edge_store[index]
    }
}

impl<E, V, H, N: NodeStorage<NodeData = V>> Index<[&Hedge; 1]> for HedgeGraph<E, V, H, N> {
    type Output = E;
    fn index(&self, index: [&Hedge; 1]) -> &Self::Output {
        &self[self[index[0]]]
    }
}

impl<E, V, H, N: NodeStorage<NodeData = V>> IndexMut<[&Hedge; 1]> for HedgeGraph<E, V, H, N> {
    // type Output = E;
    fn index_mut(&mut self, index: [&Hedge; 1]) -> &mut Self::Output {
        let edgeid = self[index[0]];
        &mut self[edgeid]
    }
}

impl<E, V, H, N: NodeStorageOps<NodeData = V>> Index<NodeIndex> for HedgeGraph<E, V, H, N> {
    type Output = V;
    fn index(&self, index: NodeIndex) -> &Self::Output {
        self.node_store.get_node_data(index)
    }
}
impl<E, V, H, N: NodeStorageOps<NodeData = V>> IndexMut<NodeIndex> for HedgeGraph<E, V, H, N> {
    fn index_mut(&mut self, index: NodeIndex) -> &mut Self::Output {
        self.node_store.get_node_data_mut(index)
    }
}
// impl<E, V, N: NodeStorageOps<NodeData = V>> Index<&NodeIndex> for HedgeGraph<E, V, N> {
//     type Output = HedgeNode;
//     fn index(&self, index: &NodeIndex) -> &Self::Output {
//         self.node_store.get_neighbor_iterator(*index)
//     }
// // }
// impl<E, V, N: NodeStorageOps<NodeData = V>> Index<&HedgeNode> for HedgeGraph<E, V, N> {
//     type Output = V;
//     fn index(&self, index: &HedgeNode) -> &Self::Output {
//         let id = self.id_from_hairs(index).unwrap();
//         &self[id]
//     }
// }

// impl<E, V, N: NodeStorageOps<NodeData = V>> IndexMut<&HedgeNode> for HedgeGraph<E, V, N> {
//     fn index_mut(&mut self, index: &HedgeNode) -> &mut Self::Output {
//         let id = self.id_from_hairs(index).unwrap();
//         &mut self[id]
//     }
// }

impl<E, V, H, N: NodeStorage<NodeData = V>> Index<EdgeIndex> for HedgeGraph<E, V, H, N> {
    type Output = E;
    fn index(&self, index: EdgeIndex) -> &Self::Output {
        &self.edge_store[index]
    }
}

impl<E, V, H, N: NodeStorage<NodeData = V>> Index<&EdgeIndex> for HedgeGraph<E, V, H, N> {
    type Output = (E, HedgePair);
    fn index(&self, index: &EdgeIndex) -> &Self::Output {
        &self.edge_store[index]
    }
}

impl<E, V, H, N: NodeStorage<NodeData = V>> IndexMut<EdgeIndex> for HedgeGraph<E, V, H, N> {
    fn index_mut(&mut self, index: EdgeIndex) -> &mut Self::Output {
        &mut self.edge_store[index]
    }
}

impl<E, V, H, N: NodeStorage<NodeData = V>> IndexMut<&EdgeIndex> for HedgeGraph<E, V, H, N> {
    fn index_mut(&mut self, index: &EdgeIndex) -> &mut Self::Output {
        &mut self.edge_store[index]
    }
}

pub struct NodeIterator<'a, E, V, H, N: NodeStorage<NodeData = V>, I = IterOnes<'a, usize, Lsb0>> {
    graph: &'a HedgeGraph<E, V, H, N>,
    edges: I,
    seen: SubSet<NodeIndex>,
}

impl<'a, E, V, H, I: Iterator<Item = Hedge>, N: NodeStorageOps<NodeData = V>> Iterator
    for NodeIterator<'a, E, V, H, N, I>
where
    N::NeighborsIter<'a>: Clone,
{
    type Item = (NodeIndex, N::NeighborsIter<'a>, &'a V);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(next) = self.edges.next() {
            let node = self.graph.neighbors(next);
            let node_pos = self.graph.id_from_crown(node.clone()).unwrap();

            if self.seen[node_pos] {
                self.next()
            } else {
                self.seen.add(node_pos);
                Some((node_pos, node, &self.graph[node_pos]))
            }
        } else {
            None
        }
    }
}

#[cfg(feature = "symbolica")]
pub mod symbolica_interop;

use subgraph::{
    BaseSubgraph, Cycle, FullOrEmpty, HedgeNode, Inclusion, InternalSubGraph, ModifySubSet,
    OrientedCut, PairwiseSubSetOps, SubSetLike, SubSetOps,
};

use thiserror::Error;
use tree::SimpleTraversalTree;

use crate::define_indexed_vec;
use crate::half_edge::subgraph::subset::SubSet;
use crate::half_edge::subgraph::{SuBitGraph, SubGraphLike, SubGraphOps};
use crate::tree::ForestNodeStore;

#[derive(Error, Debug)]
pub enum HedgeError {
    #[error("Invalid start node")]
    InvalidStart,
}

pub struct EdgeIter<'a, E, V, H, S, N: NodeStorage<NodeData = V>, I: Iterator<Item = Hedge> + 'a> {
    graph: &'a HedgeGraph<E, V, H, N>,
    included_iter: I,
    subgraph: &'a S,
}
impl<'a, E, V, H, S, N: NodeStorage<NodeData = V>> EdgeIter<'a, E, V, H, S, N, S::BaseIter<'a>>
where
    S: SubSetLike,
{
    pub fn new(graph: &'a HedgeGraph<E, V, H, N>, subgraph: &'a S) -> Self {
        EdgeIter {
            graph,
            subgraph,
            included_iter: subgraph.included_iter(),
        }
    }
}

impl<'a, E, V, H, S, N: NodeStorage<NodeData = V>> Iterator
    for EdgeIter<'a, E, V, H, S, N, S::BaseIter<'a>>
where
    S: SubSetLike,
{
    type Item = (HedgePair, EdgeData<&'a E>);

    fn next(&mut self) -> Option<Self::Item> {
        let i = self.included_iter.next()?;
        let orientation = self.graph.edge_store.orientation(i);
        let data = &self.graph[self.graph[&i]];
        if let Some(e) =
            HedgePair::from_source_with_subgraph(i, &self.graph.edge_store, self.subgraph)
        {
            Some((e, EdgeData::new(data, orientation)))
        } else {
            self.next()
        }
    }
}

#[derive(Debug, Error)]
pub enum HedgeGraphError {
    #[error("Node ({0}) that Hedge {1} points to does not contain it")]
    NodeDoesNotContainHedge(NodeIndex, Hedge),
    #[error("Nodes do not partition: {0}")]
    NodesDoNotPartition(String),
    #[error("Invalid node")]
    NoNode,
    #[error("Invalid hedge {0}")]
    InvalidHedge(Hedge),
    #[error("External hedge as included: {0}")]
    ExternalHedgeIncluded(Hedge),

    #[error("Included hedge {0} is not in node {1}")]
    NotInNode(Hedge, NodeIndex),

    #[error("Traversal Root node not in subgraph {0}")]
    RootNodeNotInSubgraph(NodeIndex),
    #[error("Invalid node {0}")]
    InvalidNode(NodeIndex),
    #[error("Invalid edge")]
    NoEdge,
    #[error("Dangling Half edge present")]
    HasIdentityHedge,
    #[error("SymbolicaError: {0}")]
    SymbolicaError(&'static str),
    #[error("InvolutionError: {0}")]
    InvolutionError(#[from] InvolutionError),
    #[error("Data length mismatch")]
    DataLengthMismatch,
    #[error("Invalid hedge pair {0:?} not equal to {1:?} for edge {2}")]
    InvalidHedgePair(HedgePair, HedgePair, EdgeIndex),
    // #[error("From file error: {0}")]
    // FromFileError(#[from] GraphFromFileError),
    // #[error("Parse error: {0}")]
    // ParseError(#[from] PestError),
}

pub mod algorithms;
pub mod hedgevec;
pub mod tree;
pub mod typed_vec;

#[cfg(feature = "drawing")]
pub mod layout;
#[cfg(test)]
mod test_graphs;
#[cfg(test)]
mod tests;
