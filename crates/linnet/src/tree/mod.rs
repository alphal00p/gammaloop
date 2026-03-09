//! Defines core data structures and traits for representing forests (collections of trees).
//!
//! This module provides multiple ways to store tree node relationships, each with different
//! trade-offs in terms of memory usage and traversal performance:
//!
//! *   [`ParentPointerStore`]: Each node stores only its data and a pointer to its parent.
//!     Memory efficient, fast upward traversal, slow downward traversal (finding children).
//! *   [`ParentChildStore`]: Uses a first-child, next-sibling representation with cyclic
//!     sibling links. Balances memory and allows efficient child/sibling iteration.
//! *   [`ChildVecStore`]: Each node stores its parent, data, and an explicit `Vec` of its children.
//!     Simple child iteration, potentially higher memory usage for nodes with many children.
//!
//! The core components are:
//! *   [`Forest<R, N>`]: The main struct representing a collection of trees. It stores root data (`R`)
//!     and uses a specific node storage strategy (`N`) that implements [`ForestNodeStore`].
//! *   [`ForestNodeStore`] trait: Defines the basic interface for node storage (adding nodes,
//!     accessing data/parent, mapping data).
//! *   [`ForestNodeStoreDown`] trait: Extends `ForestNodeStore` with methods for downward traversal
//!     (iterating children and leaves).
//! *   [`TreeNodeId`], [`RootId`], [`ParentId`]: Typed identifiers for nodes, roots, and parent links.
//!
//! Conversions between the different store types are provided via `From` implementations.

use std::{
    fmt::{Display, Write},
    ops::{Index, IndexMut, Range},
};

use child_pointer::ParentChildStore;
use child_vec::ChildVecStore;
use parent_pointer::{PPNode, ParentId, ParentPointerStore};
use thiserror::Error;

use crate::half_edge::{
    involution::Hedge,
    subgraph::{Inclusion, SuBitGraph, SubSetLike, SubSetOps},
    NodeIndex,
};

pub mod child_pointer;
pub mod child_vec;
pub mod iterato;
pub mod parent_pointer;
/// A type-safe identifier for a node within a `Forest`.
/// Wraps a `usize` index into the underlying node storage vector.
#[derive(Clone, Debug, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
pub struct TreeNodeId(pub usize);

impl Display for TreeNodeId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Node ({})", self.0)
    }
}

impl From<usize> for TreeNodeId {
    fn from(i: usize) -> Self {
        TreeNodeId(i)
    }
}

impl From<TreeNodeId> for usize {
    fn from(i: TreeNodeId) -> Self {
        i.0
    }
}

impl From<Hedge> for TreeNodeId {
    fn from(h: Hedge) -> Self {
        h.0.into()
    }
}

impl From<TreeNodeId> for Hedge {
    fn from(h: TreeNodeId) -> Self {
        Hedge(h.0)
    }
}

/// Internal data associated with the root of a tree in the `Forest`.
#[derive(Clone, Debug, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
pub struct RootData<R> {
    pub(crate) data: R,
    pub(crate) root_id: TreeNodeId,
}

/// A type-safe identifier for a tree within a `Forest`.
/// Wraps a `usize` index into the `Forest`'s roots vector.
#[derive(Clone, Debug, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
pub struct RootId(pub(crate) usize);

impl Display for RootId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Root ({})", self.0)
    }
}
impl From<NodeIndex> for RootId {
    fn from(n: NodeIndex) -> Self {
        n.0.into()
    }
}

impl From<usize> for RootId {
    fn from(i: usize) -> Self {
        RootId(i)
    }
}

/// Represents a forest (a collection of disjoint trees).
///
/// `R` is the type of data associated with each root (tree).
/// `N` is the storage implementation for the nodes, which must implement [`ForestNodeStore`].
/// Typically `N` will be one of [`ParentPointerStore<V>`], [`ParentChildStore<V>`], or [`ChildVecStore<V>`],
/// where `V` is the type of data associated with each node.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
pub struct Forest<R, N> {
    /// The underlying storage for all nodes in the forest.
    pub(crate) nodes: N,
    /// Metadata for each root/tree in the forest.
    pub(crate) roots: Vec<RootData<R>>,
}

impl<R, N> Index<&RootId> for Forest<R, N> {
    type Output = TreeNodeId;
    fn index(&self, index: &RootId) -> &Self::Output {
        &self.roots[index.0].root_id
    }
}

impl<R, N> Index<RootId> for Forest<R, N> {
    type Output = R;
    fn index(&self, index: RootId) -> &Self::Output {
        &self.roots[index.0].data
    }
}

impl<R, N> IndexMut<RootId> for Forest<R, N> {
    fn index_mut(&mut self, index: RootId) -> &mut Self::Output {
        &mut self.roots[index.0].data
    }
}

impl<R, N: ForestNodeStore> Index<&TreeNodeId> for Forest<R, N> {
    type Output = ParentId;
    fn index(&self, index: &TreeNodeId) -> &Self::Output {
        &self.nodes[index]
    }
}

impl<R, N: ForestNodeStore> Index<TreeNodeId> for Forest<R, N> {
    type Output = Option<N::NodeData>;
    fn index(&self, index: TreeNodeId) -> &Self::Output {
        &self.nodes[index]
    }
}

impl<R, N: ForestNodeStore> IndexMut<TreeNodeId> for Forest<R, N> {
    fn index_mut(&mut self, index: TreeNodeId) -> &mut Self::Output {
        &mut self.nodes[index]
    }
}

impl<R, N: Default> Default for Forest<R, N> {
    fn default() -> Self {
        Self::new()
    }
}

// --- Forest Methods (Trait-Bounded) ---

/// Methods available when the node store supports downward traversal.
impl<R, N: ForestNodeStoreDown> Forest<R, N> {
    /// Returns an iterator over the direct children of the node `start`.
    pub fn iter_children(&self, start: TreeNodeId) -> impl Iterator<Item = TreeNodeId> + '_ {
        self.nodes.iter_children(start)
    }
}

/// Methods available when the node store supports iterating leaves of a specific root.
/// (Note: This is automatically implemented if ForestNodeStoreDown is implemented).
impl<R, N: ForestNodeStoreRootLeaves> Forest<R, N> {
    /// Returns an iterator over all leaf nodes within the tree identified by `root_id`.
    pub fn iter_root_leaves(&self, root_id: RootId) -> impl Iterator<Item = TreeNodeId> + '_ {
        self.nodes.iter_root_leaves(root_id)
    }
}

/// Methods available when the node store supports pre-order traversal.
impl<R, N: ForestNodeStorePreorder> Forest<R, N> {
    /// Returns a pre-order DFS iterator starting from the given node.
    pub fn iter_preorder(&self, start: TreeNodeId) -> impl Iterator<Item = TreeNodeId> + '_ {
        self.nodes.iter_preorder(start)
    }
}

/// Methods available when the node store supports BFS traversal.
impl<R, N: ForestNodeStoreBfs> Forest<R, N> {
    /// Returns a BFS iterator starting at the given node.
    pub fn iter_bfs(&self, start: TreeNodeId) -> impl Iterator<Item = TreeNodeId> + '_ {
        self.nodes.iter_bfs(start)
    }
}

/// Methods available when the node store supports ancestor traversal.
/// (Note: This is automatically implemented for any ForestNodeStore).
impl<R, N: ForestNodeStoreAncestors> Forest<R, N> {
    /// Returns an iterator that traverses upwards from `start_node` towards its root.
    pub fn iter_ancestors(&self, start_node: TreeNodeId) -> impl Iterator<Item = TreeNodeId> + '_ {
        self.nodes.iter_ancestors(start_node)
    }
}

pub struct NodeIdIter {
    pub(crate) range: Range<usize>,
    pub(crate) current: TreeNodeId,
}

impl Iterator for NodeIdIter {
    type Item = TreeNodeId;
    fn next(&mut self) -> Option<Self::Item> {
        let out = if self.range.contains(&self.current.0) {
            Some(self.current)
        } else {
            None
        };
        self.current.0 += 1;
        out
    }
}

/// General Forest methods available for any `ForestNodeStore`.
impl<R, N: ForestNodeStore> Forest<R, N> {
    pub fn iter_roots(&self) -> impl Iterator<Item = (&R, &TreeNodeId)> {
        self.roots.iter().map(|r| (&r.data, &r.root_id))
    }
    // pub fn iter_root_ids(&self) -> impl Iterator<Item = RootId> { /* ... */
    // }

    pub fn iter_root_ids(&self) -> impl Iterator<Item = RootId> {
        (0..self.roots.len()).map(RootId)
    }

    pub fn iter_nodes(&self) -> impl Iterator<Item = (TreeNodeId, Option<&N::NodeData>)> {
        self.nodes.iter_nodes()
    }
    pub fn iter_node_ids(&self) -> NodeIdIter {
        self.nodes.iter_node_id()
    }

    pub fn root(&self, nodeid: TreeNodeId) -> RootId {
        self.nodes.root(nodeid)
    }

    pub fn map_roots<R2>(self, mut transform: impl FnMut(R) -> R2) -> Forest<R2, N> {
        Forest {
            roots: self
                .roots
                .into_iter()
                .map(|a| RootData {
                    data: transform(a.data),
                    root_id: a.root_id,
                })
                .collect(),
            nodes: self.nodes,
        }
    }

    pub fn map_roots_ref<R2>(&self, mut transform: impl FnMut(&R) -> R2) -> Forest<R2, N>
    where
        N: Clone,
    {
        Forest {
            roots: self
                .roots
                .iter()
                .map(|a| RootData {
                    data: transform(&a.data),
                    root_id: a.root_id,
                })
                .collect(),
            nodes: self.nodes.clone(),
        }
    }

    pub fn map_nodes<F, T, U>(self, transform: F) -> Forest<R, N::Store<U>>
    where
        N: ForestNodeStore<NodeData = T>,
        N::Store<U>: ForestNodeStore<NodeData = U>,
        F: FnMut(T) -> U,
    {
        Forest {
            nodes: self.nodes.map(transform),
            roots: self.roots,
        }
    }

    pub fn map_nodes_ref<F, T, U>(&self, transform: F) -> Forest<R, N::Store<U>>
    where
        N: ForestNodeStore<NodeData = T>,
        N::Store<U>: ForestNodeStore<NodeData = U>,
        F: FnMut(&T) -> U,
        R: Clone,
    {
        Forest {
            nodes: self.nodes.map_ref(transform),
            roots: self.roots.clone(),
        }
    }
}

/// Methods for constructing/modifying the forest.
impl<R, N> Forest<R, N> {
    /// Adds a new root node to the forest, creating a new tree.
    ///
    /// Takes the data for the new node (`node_data`) and the data for the new tree (`tree_data`).
    /// Returns the `TreeNodeId` of the new root node and the `RootId` of the new tree.
    /// Requires `N` to implement `ForestNodeStore<NodeData = T>`.
    pub fn add_root<T>(&mut self, node_data: T, tree_data: R) -> (TreeNodeId, RootId)
    where
        N: ForestNodeStore<NodeData = T>, // Constraint on the node store type
    {
        let root_id = RootId(self.roots.len());
        let root_node_id = self.nodes.add_root(node_data, root_id);
        self.roots.push(RootData {
            data: tree_data,
            root_id: root_node_id,
        });
        (root_node_id, root_id)
    }

    /// Adds a new child node to the specified parent node.
    ///
    /// Takes the `TreeNodeId` of the parent and the data for the new child node (`node_data`).
    /// Returns the `TreeNodeId` of the newly added child node.
    /// Requires `N` to implement `ForestNodeStore<NodeData = T>`.
    pub fn add_child<T>(&mut self, parent_id: TreeNodeId, node_data: T) -> TreeNodeId
    where
        N: ForestNodeStore<NodeData = T>, // Constraint on the node store type
    {
        self.nodes.add_child(node_data, parent_id)
    }
}

impl<R, N: Default> Forest<R, N> {
    pub fn new() -> Self {
        Forest {
            nodes: N::default(),
            roots: Vec::new(),
        }
    }
}

/// Errors that can occur during forest operations.
#[derive(Debug, Clone, Error)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
pub enum ForestError {
    #[error("Length mismatch")]
    LengthMismatch,
    #[error("Parent pointers are cyclic")]
    CyclicPP,
    #[error("Child pointers are cyclic")]
    CyclicCP,
    #[error("Self loop Parent pointers at: {0}")]
    SelfLoopPP(TreeNodeId),
    #[error("Parent pointer of child of {0} is not correct")]
    WrongPP(TreeNodeId),
    #[error("Does not partition: Sets overlap")]
    DoesNotPartion,
    #[error("Does not partition: Sets do not cover the domain")]
    IncompletePartition,
    #[error("Invalid TreeNodeId: {0:?}")]
    InvalidNodeId(TreeNodeId),
    #[error("Invalid RootId: {0:?}")]
    InvalidRootId(RootId),
    #[error("Mismatch roots len {expected}!={obtained}")]
    MismatchRootsLen { expected: usize, obtained: usize },
    #[error("Root node pointer wrong for root {0}")]
    WrongRootPointer(RootId),
}

impl<U, P: ForestNodeStore> Forest<U, P> {
    pub fn root_node(&self, nodeid: TreeNodeId) -> TreeNodeId {
        self.nodes.root_node(nodeid)
    }
    pub fn validate_structure(&self) -> Result<(), ForestError> {
        let roots = self.nodes.validate()?;

        for (rid, nid) in roots {
            if self.roots[rid.0].root_id != nid {
                return Err(ForestError::WrongRootPointer(rid));
            }
        }

        Ok(())
    }

    pub fn debug_draw<FR, FN>(
        &self,
        mut format_root: FR, // Closure to format root data (&R) -> String
        mut format_node: FN, // Closure to format node data (&N::NodeData) -> String
    ) -> String
    where
        FR: FnMut(&U) -> String,
        FN: FnMut(&Option<P::NodeData>) -> String,
        P: ForestNodeStoreDown,
    {
        let mut output = String::new();
        let _ = writeln!(output, "Forest ({} roots):", self.roots.len());

        // Recursive helper function - structure remains the same, handles sub-trees
        fn draw_subtree_recursive<W: Write, N: ForestNodeStoreDown, FmtNode>(
            f: &mut W,                 // Output writer
            nodes: &N,                 // Node store access
            node_id: TreeNodeId,       // Current node to draw
            prefix: &str,              // Current line prefix (e.g., "  │   ")
            is_last_child: bool,       // Is this the last child of its parent?
            format_node: &mut FmtNode, // Mutable reference to node formatter closure
        ) -> Result<(), std::fmt::Error>
        where
            FmtNode: FnMut(&Option<N::NodeData>) -> String, // Closure signature
        {
            // Determine the connector shape based on whether it's the last child
            let connector = if is_last_child {
                "└── "
            } else {
                "├── "
            };
            // Use Index trait which should be implemented for the store
            let node_data_ref = &nodes[node_id]; // Assuming Index<TreeNodeId> exists for N
            let node_data_str = format_node(node_data_ref);

            // Write the line for the current node
            writeln!(
                f,
                "{prefix}{connector}{node_id}:{node_data_str}", // TreeNodeId needs Debug
            )?;

            // Prepare the prefix for the children of this node
            let child_prefix = format!("{}{}", prefix, if is_last_child { "    " } else { "│   " });

            let children: Vec<_> = nodes.iter_children(node_id).collect();
            let num_children = children.len();
            for (i, child_id) in children.into_iter().enumerate() {
                // Recurse: pass the new prefix and whether this child is the last one
                draw_subtree_recursive(
                    f,
                    nodes,
                    child_id,
                    &child_prefix,
                    i == num_children - 1,
                    format_node,
                )?;
            }
            Ok(())
        } // End of helper fn definition

        let num_roots = self.roots.len();
        for (root_idx, root_meta) in self.roots.iter().enumerate() {
            let root_id = RootId(root_idx);
            let start_node_id = root_meta.root_id;

            // Format data using closures
            let root_data_str = format_root(&root_meta.data);
            // Use Index trait for Forest to access node data via store
            let root_node_data = &self[start_node_id];
            let root_node_data_str = format_node(root_node_data);

            // 1. Print the Root line
            let _ = writeln!(
                output,
                "{root_id}:{root_data_str}", // RootId needs Debug
            );

            // 2. Print the connector line
            let _ = writeln!(output, "  │"); // Fixed indent for the connector

            // 3. Define the prefix for the root node line itself and its direct children
            let node_line_prefix = "  "; // Prefix for the root TreeNodeId line

            // 4. Print the root TreeNodeId line
            let _ = writeln!(
                output,
                "{node_line_prefix}{start_node_id}:{root_node_data_str}", // TreeNodeId needs Debug
            );

            // 5. Draw children, starting the recursion with the correct prefix
            let children: Vec<_> = self.nodes.iter_children(start_node_id).collect();
            let num_children = children.len();

            // Prefix for the *children's connectors* (└── or ├──)
            // It's based on the root node's prefix
            // let children_connector_prefix = format!("{}{}", node_line_prefix, "    "); // Indent matching the root node + 4 spaces for connector alignment

            for (i, child_id) in children.into_iter().enumerate() {
                // Pass the prefix that aligns the connectors (└──, ├──)
                // under the root node line.
                let _ = draw_subtree_recursive(
                    &mut output,
                    &self.nodes,
                    child_id,
                    node_line_prefix, // Prefix aligns the connector itself
                    i == num_children - 1,
                    &mut format_node,
                );
            }

            // 6. Add separation between trees (only if not the last root)
            if root_idx < num_roots - 1 {
                let _ = writeln!(output); // Just a blank line for separation
            }
        }

        output
    }
    pub fn swap_roots(&mut self, a: RootId, b: RootId) {
        if a == b {
            return;
        }

        let roota = self.roots[a.0].root_id;
        if roota == self.nodes.root_node(roota) {
            self.nodes.set_root(roota, b);
        }
        let rootb = self.roots[b.0].root_id;
        if rootb == self.nodes.root_node(rootb) {
            self.nodes.set_root(rootb, a);
        }

        self.roots.swap(a.0, b.0);
    }

    pub fn from_bitvec_partition<I: IntoIterator<Item = (U, SuBitGraph)>>(
        bitvec_part: I,
    ) -> Result<Self, ForestError> {
        Forest::<U, ParentPointerStore<P::NodeData>>::from_bitvec_partition_impl(bitvec_part).map(
            |a| Forest {
                nodes: P::from_pp(a.nodes),
                roots: a.roots,
            },
        )
    }
}

/// Construction from BitVec partition (remains the same)
impl<U, V> Forest<U, ParentPointerStore<V>> {
    /// Creates a forest from a partition represented by a vector of (RootData, BitVec).
    /// Each BitVec represents the set of node indices belonging to that tree/root.
    /// Nodes within a set are linked arbitrarily (first encountered becomes root, others point to it).
    /// Node data is set to `()`.
    ///
    /// Returns `Err(ForestError)` if the BitVecs have different lengths or overlap.
    fn from_bitvec_partition_impl<I: IntoIterator<Item = (U, SuBitGraph)>>(
        bitvec_part: I,
    ) -> Result<Self, ForestError> {
        let mut nodes = vec![];
        let mut roots = vec![];
        let mut cover: Option<SuBitGraph> = None;

        for (d, set) in bitvec_part {
            let len = set.size();
            if let Some(c) = &mut cover {
                if c.size() != len {
                    return Err(ForestError::LengthMismatch);
                }
                if c.intersects(&set) {
                    return Err(ForestError::DoesNotPartion);
                }
                c.union_with(&set);
            } else {
                cover = Some(SuBitGraph::empty(len));

                nodes.resize_with(len, || None);
                // nodes = vec![None; len];
            }
            let mut first = None;
            for i in set.included_iter() {
                if let Some(root) = first {
                    nodes[i.0] = Some(PPNode::dataless_child(root))
                } else {
                    first = Some(i.into());
                    nodes[i.0] = Some(PPNode::dataless_root(RootId(roots.len())));
                }
            }
            roots.push(RootData {
                root_id: first.unwrap(),
                data: d,
            });
        }
        Ok(Forest {
            nodes: nodes
                .into_iter()
                .collect::<Option<Vec<_>>>()
                .unwrap()
                .into_iter()
                .collect(),
            roots,
        })
    }
}

// --- Core Traits ---

/// The core trait defining the interface for node storage within a `Forest`.
/// Requires Index and IndexMut support for accessing node data and ParentId.
pub trait ForestNodeStore:
    for<'a> Index<&'a TreeNodeId, Output = ParentId>      // `store[&node_id]` -> ParentId
    + Index<TreeNodeId, Output = Option<Self::NodeData>>          // `store[node_id]` -> NodeData
    + IndexMut<TreeNodeId, Output = Option<Self::NodeData>>   // `store[node_id] = ...`
    + Sized // Needed for some default implementations returning Self iterators
{
    type NodeData;
    type Store<T>: ForestNodeStore<NodeData = T>+From<ParentPointerStore<T>>+Into<ParentPointerStore<T>>;

    fn parent(&self,child:TreeNodeId)->ParentId{
        self[&child]
    }





    fn to_store(self)->Self::Store<Self::NodeData>;
    fn from_store(store:Self::Store<Self::NodeData>)->Self;
    fn debug_draw(&self,node_display:impl FnMut(Option<&Self::NodeData>)->Option<String>)->String;
    fn validate(&self)->Result<Vec<(RootId,TreeNodeId)>,ForestError>;

    fn panicing_validate(&self){
        self.validate()
            .unwrap_or_else(|a| panic!("Err:{}\nforest:\n{}", a, self.debug_draw(|_| None)));
    }

    fn n_nodes(&self)->usize;
    fn swap(&mut self,a:TreeNodeId,b:TreeNodeId);
    fn set_root(&mut self,a:TreeNodeId,root:RootId);

    fn from_pp(pp:ParentPointerStore<Self::NodeData>)->Self{
        Self::from_store(Self::Store::<Self::NodeData>::from(pp))
    }



    fn split_off(&mut self,at:TreeNodeId)->Self;


    /// Makes x's root, the root of y's root.
    fn make_root_of(&mut self, x: TreeNodeId, y: TreeNodeId) -> TreeNodeId{
        let rootx = self.root_node(x);
        let rooty = self.root_node(y);
        self.reparent(rootx,rooty);
        rootx
    }

    fn reparent(&mut self,parent:TreeNodeId,child:TreeNodeId);

    /// Finds the `RootId` by traversing upwards. Default implementation provided.
    fn root(&self, nodeid: TreeNodeId) -> RootId { let mut current = nodeid;
            loop {
                match self[&current] {
                    ParentId::PointingRoot(p)=>current=p,
                    ParentId::Root(root_id) => return root_id,
                    ParentId::Node(parent_node_id) => current = parent_node_id,
                }
            } }

    fn root_node(&self,nodeid:TreeNodeId)->TreeNodeId{
        let mut current = nodeid;
        loop {
            match self[&current] {
                ParentId::PointingRoot(r)=> current=r,
                ParentId::Root(_) => return current,
                ParentId::Node(parent_node_id) => current = parent_node_id,
            }
        }
    }

    fn iter_nodes(&self) -> impl Iterator<Item = (TreeNodeId, Option<&Self::NodeData>)>;
    fn iter_node_id(&self) -> NodeIdIter {
        NodeIdIter {
            range: 0..self.n_nodes(),
            current: TreeNodeId(0),
        }
    }
    /// Adds a new root node.
    fn add_root(&mut self, data: Self::NodeData, root_id: RootId) -> TreeNodeId;

    /// Adds a new child node as the *last* child of the parent.
    fn add_child(&mut self, data: Self::NodeData, parent: TreeNodeId) -> TreeNodeId;


    fn set_node_data(&mut self,data:Self::NodeData,node_id:TreeNodeId)->Option<Self::NodeData>;
    /// Adds a new child node as the *last* child of the parent.
    fn add_dataless_child(&mut self, parent: TreeNodeId) -> TreeNodeId;

    fn map<F, U>(self, transform: F) -> Self::Store<U> where F: FnMut(Self::NodeData) -> U;
    fn map_ref<F, U>(&self, transform: F) -> Self::Store<U> where F: FnMut(&Self::NodeData) -> U;
    fn forgetful_map<U>(&self)->Self::Store<U>;

    fn extend(&mut self,other:Self,shift_roots_by:RootId);
}

/// Trait extension for `ForestNodeStore` providing downward traversal capabilities.
pub trait ForestNodeStoreDown: ForestNodeStore {
    /// Returns an iterator over all leaf nodes in the entire store.
    fn iter_leaves(&self) -> impl Iterator<Item = TreeNodeId> + '_;

    /// Returns an iterator over the direct children of the given `node_id`.
    fn iter_children(&self, node_id: TreeNodeId) -> impl Iterator<Item = TreeNodeId> + '_;
}

// --- Traversal Traits ---

/// Trait for stores supporting ancestor iteration (upward traversal).
pub trait ForestNodeStoreAncestors: ForestNodeStore {
    /// Returns an iterator from `start_node` up to its root (inclusive).
    fn iter_ancestors(&self, start_node: TreeNodeId) -> impl Iterator<Item = TreeNodeId> + '_;
}

/// Trait for stores supporting pre-order traversal.
pub trait ForestNodeStorePreorder: ForestNodeStore {
    type Iterator<'a>: Iterator<Item = TreeNodeId> + Clone
    where
        Self: 'a;
    /// Returns a pre-order DFS iterator starting from `start`.
    fn iter_preorder(&self, start: TreeNodeId) -> Self::Iterator<'_>;
}

/// Trait for stores supporting breadth-first traversal.
pub trait ForestNodeStoreBfs: ForestNodeStore {
    /// Returns a BFS iterator starting from `start`.
    fn iter_bfs(&self, start: TreeNodeId) -> impl Iterator<Item = TreeNodeId> + '_;
}

/// Trait for stores supporting iteration over leaves belonging to a specific root.
pub trait ForestNodeStoreRootLeaves: ForestNodeStoreDown {
    /// Returns an iterator over the leaf nodes belonging only to the tree with `root_id`.
    fn iter_root_leaves(&self, root_id: RootId) -> impl Iterator<Item = TreeNodeId> + '_;
}

// --- Default Trait Implementations ---

/// Default implementation for ancestor iteration for any `ForestNodeStore`.
impl<N: ForestNodeStore> ForestNodeStoreAncestors for N {
    fn iter_ancestors(&self, start_node: TreeNodeId) -> impl Iterator<Item = TreeNodeId> + '_ {
        // Use the helper struct from the iterators module
        iterato::AncestorsIter::new(self, start_node)
    }
}

/// Default implementation for root leaves iteration for any `ForestNodeStoreDown`.
impl<N: ForestNodeStoreDown> ForestNodeStoreRootLeaves for N {
    fn iter_root_leaves(&self, root_id: RootId) -> impl Iterator<Item = TreeNodeId> + '_ {
        self.iter_leaves()
            .filter(move |&leaf_id| self.root(leaf_id) == root_id)
    }
}

impl<P, V> Forest<V, P> {
    pub fn cast<PP>(self) -> Forest<V, PP>
    where
        PP: From<P>,
    {
        Forest {
            nodes: self.nodes.into(),
            roots: self.roots,
        }
    }
}

impl<R, V> From<Forest<R, ParentPointerStore<V>>> for Forest<R, ParentChildStore<V>> {
    fn from(forest: Forest<R, ParentPointerStore<V>>) -> Self {
        Forest {
            nodes: forest.nodes.into(),
            roots: forest.roots,
        }
    }
}

impl<R, V> From<Forest<R, ParentChildStore<V>>> for Forest<R, ParentPointerStore<V>> {
    fn from(forest: Forest<R, ParentChildStore<V>>) -> Self {
        Forest {
            nodes: forest.nodes.into(),
            roots: forest.roots,
        }
    }
}

impl<R, V> From<Forest<R, ChildVecStore<V>>> for Forest<R, ParentChildStore<V>> {
    fn from(forest: Forest<R, ChildVecStore<V>>) -> Self {
        Forest {
            nodes: forest.nodes.into(),
            roots: forest.roots,
        }
    }
}

impl<R, V> From<Forest<R, ParentChildStore<V>>> for Forest<R, ChildVecStore<V>> {
    fn from(forest: Forest<R, ParentChildStore<V>>) -> Self {
        Forest {
            nodes: forest.nodes.into(),
            roots: forest.roots,
        }
    }
}

impl<R, V> From<Forest<R, ChildVecStore<V>>> for Forest<R, ParentPointerStore<V>> {
    fn from(forest: Forest<R, ChildVecStore<V>>) -> Self {
        Forest {
            nodes: forest.nodes.into(),
            roots: forest.roots,
        }
    }
}

impl<R, V> From<Forest<R, ParentPointerStore<V>>> for Forest<R, ChildVecStore<V>> {
    fn from(forest: Forest<R, ParentPointerStore<V>>) -> Self {
        Forest {
            nodes: forest.nodes.into(),
            roots: forest.roots,
        }
    }
}
