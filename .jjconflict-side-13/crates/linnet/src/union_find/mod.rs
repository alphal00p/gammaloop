//! # Union-Find (Disjoint Set Union - DSU) Data Structure
//!
//! This module implements a Union-Find data structure, also known as Disjoint
//! Set Union (DSU). It is used to keep track of a collection of disjoint sets
//! and supports two primary operations:
//!
//! 1.  **Find**: Determine which set a particular element belongs to (i.e., find the
//!     representative or root of the set).
//! 2.  **Union**: Merge two sets into a single set.
//!
//! This implementation includes common optimizations:
//! - **Path Compression**: During a `find` operation, nodes along the path to the
//!   root are made to point directly to the root. This flattens the tree structure,
//!   speeding up future operations.
//! - **Union by Rank**: When merging two sets, the root of the tree with a smaller
//!   "rank" (typically related to tree depth or size) is made a child of the root
//!   with a larger rank. This helps keep the trees relatively shallow.
//!
//! ## Core Components:
//!
//! - **`UnionFind<U>`**: The main struct.
//!   - `nodes: Vec<Cell<UFNode>>`: Stores the parent-pointer forest. Each node is
//!     wrapped in a `Cell` to allow for path compression (interior mutability)
//!     during `find` operations.
//!   - `set_data: Vec<SetData<U>>`: Stores data of type `U` associated with the
//!     representative (root) of each disjoint set.
//!
//! - **`UFNode`**: An enum representing a node within the Union-Find forest:
//!   - `Root { set_data_idx: SetIndex, rank: usize }`: Indicates that this node is
//!     the representative of its set. It holds an index (`set_data_idx`) to its
//!     associated data in `set_data` and its `rank`.
//!   - `Child(ParentPointer)`: Indicates that this node is a child of another node,
//!     specified by `ParentPointer`.
//!
//! - **`SetData<U>`**: A struct holding the actual data `Option<U>` for a set and
//!   the `ParentPointer` to the root node that currently owns this data slot.
//!
//! - **`ParentPointer`**: A type alias (currently for `Hedge`) used as an identifier
//!   for nodes within the Union-Find structure.
//! - **`SetIndex`**: A type alias for `usize`, used to index into the `set_data` vector.
//!
//! ## Key Operations:
//!
//! - `UnionFind::new(associated_data: Vec<U>)`: Creates a new Union-Find structure
//!   where each initial element is in its own set, with the provided associated data.
//! - `find(node: ParentPointer) -> ParentPointer`: Finds the representative of the
//!   set containing `node`. Implements path compression.
//! - `union(x: ParentPointer, y: ParentPointer, merge_fn: F) -> ParentPointer`: Merges
//!   the sets containing `x` and `y`. Uses union-by-rank. The `merge_fn` is used
//!   to combine the data associated with the two sets.
//! - `find_data(node: ParentPointer) -> &U`: Returns a reference to the data
//!   associated with the set containing `node`.
//! - `map_set_data_of(node: ParentPointer, f: F)`: Allows mutable access to the data
//!   of `node`'s set.
//! - `extract(...)`: A more advanced operation to partition the Union-Find structure
//!   itself based on a predicate, potentially creating a new `UnionFind` instance.
//! - `from_bitvec_partition(partitions: Vec<(U, BitVec)>)`: Constructs a `UnionFind`
//!   from a set of partitions defined by bitvectors.
//!
//! This data structure is commonly used in algorithms for graph connectivity,
//! Kruskal's algorithm for minimum spanning trees, and other problems involving
//! partitioning or equivalence classes.

use std::{
    cell::Cell,
    ops::{Index, IndexMut},
};

use thiserror::Error;

use crate::half_edge::{
    involution::Hedge,
    subgraph::{subset::SubSet, Inclusion, SubSetLike, SubSetOps},
    typed_vec::IndexLike,
};

/// A newtype for a node (index into `self.nodes`).
// #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub type ParentPointer = Hedge;

/// A newtype for the set–data index (index into `UnionFind::set_data`).
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
pub struct SetIndex(pub usize);

impl From<usize> for SetIndex {
    fn from(x: usize) -> Self {
        SetIndex(x)
    }
}

/// A node can be:
/// - `Root { set_data_idx, rank }`: this node is a root, with `rank` for union–by–rank,
///   and it owns the data at `set_data_idx`.
/// - `Child(parent)`: a non–root pointing to another node's index.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
pub enum UFNode {
    Root { set_data_idx: SetIndex, rank: usize },
    Child(ParentPointer),
}

#[derive(Debug, Clone, Error)]
pub enum UnionFindError {
    #[error("The set of bitvecs does not partion the elements")]
    DoesNotPartion,
    #[error("The set of bitvecs does not have the same length")]
    LengthMismatch,
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
/// A UnionFind structure that:
/// - Stores a parallel `elements: Vec<T>` (one per node).
/// - Maintains a parent–pointer forest (`Vec<Cell<UFNode>>`).
/// - Stores associated data (`U`) for each root in `set_data: Vec<Option<U>>`,
///   using swap–removal when merging.
/// - `data_to_node` is the inverse of `set_data_idx`, telling us which node currently owns
///   each slot in `set_data`.
pub struct UnionFind<U> {
    /// Each node is a `Cell<UFNode>` for in-place mutation during path compression.
    pub nodes: Vec<Cell<UFNode>>,

    /// For each root, there's exactly one `Some(U)` slot here.
    /// Non–roots may have been swapped out to maintain compactness.
    pub(crate) set_data: Vec<SetData<U>>,
    // forest:Forest<>
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
pub struct SetData<U> {
    pub(crate) root_pointer: ParentPointer,
    pub(crate) data: Option<U>,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
enum SetSplit {
    Left,
    FullyExtracted,
    Split,
    Unset,
}

pub fn left<E>(l: E, _: E) -> E {
    l
}

pub fn right<E>(_: E, r: E) -> E {
    r
}

impl<U> UnionFind<U> {
    pub fn n_elements(&self) -> usize {
        self.nodes.len()
    }

    pub fn swap_set_data(&mut self, a: SetIndex, b: SetIndex) {
        self.set_data.swap(a.0, b.0);
        for n in &mut self.nodes {
            if let UFNode::Root { set_data_idx, .. } = n.get_mut() {
                if *set_data_idx == a {
                    *set_data_idx = a;
                } else if *set_data_idx == b {
                    *set_data_idx = b;
                }
            }
        }
    }

    // if true,then extract
    #[allow(clippy::manual_map)]
    pub fn extract<O>(
        &mut self,
        mut part: impl FnMut(ParentPointer) -> bool,
        mut split_set_data: impl FnMut(&U) -> O,
        mut extract_set_data: impl FnMut(U) -> O,
    ) -> UnionFind<O> {
        let mut left = Hedge(0);
        let mut extracted = Hedge(self.nodes.len());
        let mut set_split = vec![SetSplit::Unset; self.set_data.len()];

        while left < extracted {
            if !part(left) {
                //left is in the right place
                left.0 += 1;
                let root = self.find_data_index(left);
                match set_split[root.0] {
                    SetSplit::Left => {}
                    SetSplit::Unset => set_split[root.0] = SetSplit::Left,
                    SetSplit::FullyExtracted => set_split[root.0] = SetSplit::Split,
                    SetSplit::Split => {}
                }
            } else {
                //left needs to be swapped
                extracted.0 -= 1;
                let root = self.find_data_index(extracted);
                if !part(extracted) {
                    //only with an extracted that is in the wrong spot
                    match set_split[root.0] {
                        SetSplit::Left => {}
                        SetSplit::Unset => set_split[root.0] = SetSplit::Left,
                        SetSplit::FullyExtracted => set_split[root.0] = SetSplit::Split,
                        SetSplit::Split => {}
                    }
                    self.swap(left, extracted);
                    left.0 += 1;
                } else {
                    match set_split[root.0] {
                        SetSplit::Left => set_split[root.0] = SetSplit::Split,
                        SetSplit::Unset => set_split[root.0] = SetSplit::FullyExtracted,
                        SetSplit::FullyExtracted => {}
                        SetSplit::Split => {}
                    }
                }
            }
        }

        let mut left_nodes = SetIndex(0);
        let mut extracted_nodes = SetIndex(self.n_sets());
        while left_nodes < extracted_nodes {
            if let SetSplit::Left = set_split[left_nodes.0] {
                //left is in the right place
                left_nodes.0 += 1;
            } else {
                //left needs to be swapped
                extracted_nodes.0 -= 1;
                if let SetSplit::Left = set_split[extracted_nodes.0] {
                    //only with an extracted that is in the wrong spot
                    self.swap_set_data(left_nodes, extracted_nodes);
                    left_nodes.0 += 1;
                }
            }
        }
        let mut overlapping_nodes = left_nodes;
        let mut non_overlapping_extracted = SetIndex(self.n_sets());

        while overlapping_nodes < non_overlapping_extracted {
            if let SetSplit::Split = set_split[overlapping_nodes.0] {
                //overlapping is in the right place, is a split node
                overlapping_nodes.0 += 1;
            } else {
                //overlapping needs to be swapped
                non_overlapping_extracted.0 -= 1;
                if let SetSplit::Split = set_split[non_overlapping_extracted.0] {
                    //only with an extracted that is in the wrong spot
                    self.swap_set_data(overlapping_nodes, non_overlapping_extracted);
                    overlapping_nodes.0 += 1;
                }
            }
        }

        let mut extracted_nodes = self.nodes.split_off(overlapping_nodes.0);
        let extracted_data: Vec<_> = self
            .set_data
            .split_off(overlapping_nodes.0)
            .into_iter()
            .map(|a| SetData {
                data: if let Some(data) = a.data {
                    Some(extract_set_data(data))
                } else {
                    None
                },
                root_pointer: a.root_pointer,
            })
            .collect();

        let mut overlapping_data = vec![];

        for i in (left_nodes.0)..(overlapping_nodes.0) {
            let data = Some(split_set_data(&self[SetIndex(i)]));
            let root_pointer = self[&SetIndex(i)].root_pointer;

            overlapping_data.push(SetData { root_pointer, data })
        }

        overlapping_data.extend(extracted_data);

        for n in &mut extracted_nodes {
            if let UFNode::Root { set_data_idx, .. } = n.get_mut() {
                set_data_idx.0 -= left_nodes.0
            }
        }

        UnionFind {
            nodes: extracted_nodes,
            set_data: overlapping_data,
        }
    }

    pub fn swap(&mut self, a: ParentPointer, b: ParentPointer) {
        match (self.is_child(a), self.is_child(b)) {
            (true, true) => {}
            (true, false) => {
                for n in &mut self.nodes {
                    if let UFNode::Child(pp) = n.get_mut() {
                        if *pp == a {
                            *pp = b;
                        }
                    }
                }
            }
            (false, true) => {
                for n in &mut self.nodes {
                    if let UFNode::Child(pp) = n.get_mut() {
                        if *pp == b {
                            *pp = a;
                        }
                    }
                }
            }
            (false, false) => {
                for n in &mut self.nodes {
                    if let UFNode::Child(pp) = n.get_mut() {
                        if *pp == a {
                            *pp = b;
                        } else if *pp == b {
                            *pp = a;
                        }
                    }
                }
            }
        }
        self.nodes.swap(a.0, b.0);
    }

    pub fn n_sets(&self) -> usize {
        self.set_data.len()
    }

    pub fn extend(&mut self, mut other: Self) {
        let shift_nodes = self.set_data.len();
        let shift_set = self.nodes.len();
        other.nodes.iter_mut().for_each(|a| match a.get_mut() {
            UFNode::Child(a) => a.0 += shift_set,
            UFNode::Root { set_data_idx, .. } => set_data_idx.0 += shift_nodes,
        });
        other
            .set_data
            .iter_mut()
            .for_each(|a| a.root_pointer.0 += shift_set);
        self.set_data.extend(other.set_data);
    }

    pub fn iter_set_data(&self) -> impl Iterator<Item = (SetIndex, &U)> {
        self.set_data
            .iter()
            .enumerate()
            .map(|(i, d)| (SetIndex(i), d.data.as_ref().unwrap()))
    }

    pub fn iter_set_data_mut(&mut self) -> impl Iterator<Item = (SetIndex, &mut U)> {
        self.set_data
            .iter_mut()
            .enumerate()
            .map(|(i, d)| (SetIndex(i), d.data.as_mut().unwrap()))
    }

    pub fn drain_set_data(self) -> impl Iterator<Item = (SetIndex, U)> {
        self.set_data
            .into_iter()
            .enumerate()
            .map(|(i, d)| (SetIndex(i), d.data.unwrap()))
    }

    pub fn from_partition<ID: IndexLike>(
        bitvec_part: Vec<(U, SubSet<ID>)>,
    ) -> Result<Self, UnionFindError> {
        let mut nodes = vec![];
        let mut set_data = vec![];
        let mut cover: Option<SubSet<ID>> = None;

        for (d, set) in bitvec_part {
            let len = set.size();
            if let Some(c) = &mut cover {
                if c.size() != len {
                    return Err(UnionFindError::LengthMismatch);
                }
                if c.intersects(&set) {
                    return Err(UnionFindError::DoesNotPartion);
                }
                c.union_with(&set);
            } else {
                cover = Some(SubSet::empty(len));
                nodes = vec![None; len];
            }
            let mut first = None;
            for i in set.included_iter() {
                if let Some(root) = first {
                    nodes[i.into()] = Some(Cell::new(UFNode::Child(root)))
                } else {
                    first = Some(Hedge(i.into()));
                    nodes[i.into()] = Some(Cell::new(UFNode::Root {
                        set_data_idx: SetIndex(set_data.len()),
                        rank: set.n_included(),
                    }))
                }
            }
            set_data.push(SetData {
                root_pointer: first.unwrap(),
                data: Some(d),
            });
        }
        Ok(UnionFind {
            nodes: nodes.into_iter().collect::<Option<_>>().unwrap(),
            set_data,
        })
    }

    /// Builds a union-find where each node is its own set, with `rank=0` and `SetIndex(i)` owning
    /// the `i`th slot in `set_data`.
    pub fn new(associated: Vec<U>) -> Self {
        // let n = elements.len();
        // assert_eq!(n, associated.len());
        let nodes = (0..associated.len())
            .map(|i| {
                Cell::new(UFNode::Root {
                    set_data_idx: SetIndex(i),
                    rank: 0,
                })
            })
            .collect();

        let set_data = associated
            .into_iter()
            .enumerate()
            .map(|(i, d)| SetData {
                root_pointer: Hedge(i),
                data: Some(d),
            })
            .collect();

        Self {
            // elements,
            nodes,
            set_data,
        }
    }

    /// **Find** the representative (root) of the set containing `x`, path compressing along the way.
    pub fn find(&self, x: ParentPointer) -> ParentPointer {
        match self[&x].get() {
            UFNode::Root { .. } => x,
            UFNode::Child(parent) => {
                let root = self.find(parent);
                // path compression
                self[&x].set(UFNode::Child(root));
                root
            }
        }
    }

    pub fn is_child(&self, x: ParentPointer) -> bool {
        matches!(self[&x].get(), UFNode::Child(_))
        // if let UFNode::Child(_) = self[&x].get() {
        //     true
        // } else {
        //     false
        // }
    }

    /// Returns the `SetIndex` for the set containing `x`.
    pub fn find_data_index(&self, x: ParentPointer) -> SetIndex {
        let root = self.find(x);
        match self[&root].get() {
            UFNode::Root { set_data_idx, .. } => set_data_idx,
            UFNode::Child(_) => unreachable!("find always returns a root"),
        }
    }

    /// Returns a shared reference to the data for `x`'s set, unwrapping the `Option`.
    /// Panics if no data is present (which shouldn't happen for a valid root).
    pub fn find_data(&self, x: ParentPointer) -> &U {
        &self[self.find_data_index(x)]
    }

    /// **Union** the sets containing `x` and `y`, merging their data with `merge(U, U) -> U`.
    ///
    /// - Union–by–rank
    /// - Merged data is placed in the winner’s slot.
    /// - Loser’s slot is swap–removed from `set_data`.
    /// - Returns the new root.
    pub fn union<F>(&mut self, x: ParentPointer, y: ParentPointer, merge: F) -> ParentPointer
    where
        F: FnOnce(U, U) -> U,
    {
        let rx = self.find(x);
        let ry = self.find(y);
        if rx == ry {
            return rx;
        }

        // Extract rank + data index from each root
        let (rank_x, data_x) = match self[&rx].get() {
            UFNode::Root { rank, set_data_idx } => (rank, set_data_idx),
            _ => unreachable!(),
        };
        let (rank_y, data_y) = match self[&ry].get() {
            UFNode::Root { rank, set_data_idx } => (rank, set_data_idx),
            _ => unreachable!(),
        };

        let (winner, loser, winner_data_idx, loser_data_idx, same_rank) = match rank_x.cmp(&rank_y)
        {
            std::cmp::Ordering::Less => (ry, rx, data_y, data_x, false),
            std::cmp::Ordering::Greater => (rx, ry, data_x, data_y, false),
            std::cmp::Ordering::Equal => (rx, ry, data_x, data_y, true),
        };

        if same_rank {
            if let UFNode::Root { set_data_idx, rank } = self[&winner].get() {
                self[&winner].set(UFNode::Root {
                    set_data_idx,
                    rank: rank + 1,
                });
            }
        }

        // Make loser point to winner
        self[&loser].set(UFNode::Child(winner));

        // Merge their data
        // We can now fetch the `Option<U>` for each side using `uf[&SetIndex]`.
        let winner_opt = self[&winner_data_idx].data.take();
        let loser_opt = self[&loser_data_idx].data.take();

        // Take out the two `U`s, merge them, store back in winner's slot
        let merged = merge(
            winner_opt.expect("winner has no data?"),
            loser_opt.expect("loser has no data?"),
        );
        self[&winner_data_idx].data = Some(merged);

        // Swap–remove from the losing slot
        let last_idx = self.set_data.len() - 1;
        if loser_data_idx.0 != last_idx {
            // Swap data
            self.set_data.swap(loser_data_idx.0, last_idx);

            // Fix the node that got swapped in
            let swapped_node = self.set_data[loser_data_idx.0].root_pointer;
            if let UFNode::Root { set_data_idx, rank } = self[&swapped_node].get() {
                // set_data_idx might have been the old `last_idx`.
                // Now it must be `loser_data_idx`.
                if set_data_idx.0 == last_idx {
                    self[&swapped_node].set(UFNode::Root {
                        set_data_idx: loser_data_idx,
                        rank,
                    });
                }
            }

            // If the winner's data was at last_idx, fix it
            if winner_data_idx.0 == last_idx {
                if let UFNode::Root { set_data_idx, rank } = self[&winner].get() {
                    if set_data_idx.0 == last_idx {
                        // point it to the new location
                        self[&winner].set(UFNode::Root {
                            set_data_idx: loser_data_idx,
                            rank,
                        });
                    }
                }
            }
        }
        self.set_data.pop();

        winner
    }

    /// Allows mutating the set–data of `x` in place, unwrapping the `Option`.
    pub fn map_set_data_of<F>(&mut self, x: ParentPointer, f: F)
    where
        F: FnOnce(&mut U),
    {
        let idx = self.find_data_index(x);
        let data_ref = &mut self[idx]; // &mut U
        f(data_ref);
    }

    pub fn map_set_data<V, F>(self, mut f: F) -> UnionFind<V>
    where
        F: FnMut(SetIndex, U) -> V,
    {
        UnionFind {
            nodes: self.nodes,
            set_data: self
                .set_data
                .into_iter()
                .enumerate()
                .map(|(i, a)| SetData {
                    data: Some(f(SetIndex(i), a.data.unwrap())),
                    root_pointer: a.root_pointer,
                })
                .collect(),
        }
    }

    pub fn map_set_data_ref<'a, F, V>(&'a self, mut f: F) -> UnionFind<V>
    where
        F: FnMut(&'a U) -> V,
    {
        UnionFind {
            nodes: self.nodes.clone(),
            set_data: self
                .set_data
                .iter()
                .map(|d| SetData {
                    root_pointer: d.root_pointer,
                    data: d.data.as_ref().map(&mut f),
                })
                .collect(),
        }
    }

    pub fn map_set_data_ref_mut<'a, F, V>(&'a mut self, mut f: F) -> UnionFind<V>
    where
        F: FnMut(&'a mut U) -> V,
    {
        UnionFind {
            nodes: self.nodes.clone(),
            set_data: self
                .set_data
                .iter_mut()
                .map(|d| SetData {
                    root_pointer: d.root_pointer,
                    data: d.data.as_mut().map(&mut f),
                })
                .collect(),
        }
    }

    pub fn map_set_data_ref_result<'a, F, V, Er>(&'a self, mut f: F) -> Result<UnionFind<V>, Er>
    where
        F: FnMut(&'a U) -> Result<V, Er>,
    {
        let r: Result<Vec<_>, Er> = self
            .set_data
            .iter()
            .map(|d| {
                let data = d.data.as_ref().map(&mut f).transpose();
                match data {
                    Ok(data) => Ok(SetData {
                        root_pointer: d.root_pointer,
                        data,
                    }),
                    Err(err) => Err(err),
                }
            })
            .collect();
        Ok(UnionFind {
            nodes: self.nodes.clone(),
            set_data: r?,
        })
    }

    /// Takes ownership of the old data for `x`'s set, applies a function, and replaces it.
    pub fn replace_set_data_of<F>(&mut self, x: ParentPointer, f: F)
    where
        F: FnOnce(U) -> U,
    {
        let idx = self.find_data_index(x);
        let old_data = self[&idx].data.take().expect("no data to replace");
        self[&idx].data.replace(f(old_data));
    }

    pub fn add_child(&mut self, set_id: SetIndex) -> ParentPointer {
        let root = self[&set_id].root_pointer;
        let h = Hedge(self.nodes.len());
        self.nodes.push(Cell::new(UFNode::Child(root)));
        h
    }
}

// -------------------------------------------------------------------
// Index impls
// -------------------------------------------------------------------

/// 1) `impl Index<SetIndex>` => returns `&U` (unwrapped from `Option<U>`).
impl<U> Index<SetIndex> for UnionFind<U> {
    type Output = U;
    fn index(&self, idx: SetIndex) -> &Self::Output {
        self.set_data[idx.0]
            .data
            .as_ref()
            .expect("no data in that slot!")
    }
}

/// 1b) `impl IndexMut<SetIndex>` => returns `&mut U`.
impl<U> IndexMut<SetIndex> for UnionFind<U> {
    fn index_mut(&mut self, idx: SetIndex) -> &mut Self::Output {
        self.set_data[idx.0]
            .data
            .as_mut()
            .expect("no data in that slot!")
    }
}

/// 2) `impl Index<&SetIndex>` => returns `&Option<U>`.
///    This lets you see if it’s Some/None, or call methods like `.take()`.
impl<U> Index<&SetIndex> for UnionFind<U> {
    type Output = SetData<U>;
    fn index(&self, idx: &SetIndex) -> &Self::Output {
        &self.set_data[idx.0]
    }
}

/// 2b) `impl IndexMut<&SetIndex>` => returns `&mut Option<U>`.
impl<U> IndexMut<&SetIndex> for UnionFind<U> {
    fn index_mut(&mut self, idx: &SetIndex) -> &mut Self::Output {
        &mut self.set_data[idx.0]
    }
}

// /// `impl Index<ParentPointer>` => returns `&T`.
// /// This is the direct element in `elements`.
// impl<T, U> Index<ParentPointer> for UnionFind<T, U> {
//     type Output = T;
//     fn index(&self, idx: ParentPointer) -> &Self::Output {
//         &self.elements[idx.0]
//     }
// }

// /// `impl IndexMut<ParentPointer>` => returns `&mut T`.
// impl<T, U> IndexMut<ParentPointer> for UnionFind<T, U> {
//     fn index_mut(&mut self, idx: ParentPointer) -> &mut Self::Output {
//         &mut self.elements[idx.0]
//     }
// }

/// `impl Index<&ParentPointer>` => returns `&Cell<UFNode>`,
/// allowing `self[&x].get()` or `self[&x].set(...)`.
impl<U> Index<&ParentPointer> for UnionFind<U> {
    type Output = Cell<UFNode>;
    fn index(&self, idx: &ParentPointer) -> &Self::Output {
        &self.nodes[idx.0]
    }
}

// pub mod bitvec_find;
#[cfg(test)]
pub mod test;
