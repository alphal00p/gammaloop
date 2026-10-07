//! Implements a tree structure using a first-child, next-sibling representation.
//! Each node stores its parent, data, a pointer to its first child, and pointers
//! to its left and right siblings using a *cyclic* doubly linked list for siblings.

use std::{collections::VecDeque, fmt::Write};

use crate::half_edge::subgraph::{subset::SubSet, Inclusion, ModifySubSet, SubSetLike};

use super::{
    child_vec::ChildVecStore,
    parent_pointer::{PPNode, ParentId, ParentPointerStore},
    ForestError, ForestNodeStore, ForestNodeStoreAncestors, ForestNodeStoreBfs,
    ForestNodeStoreDown, ForestNodeStorePreorder, RootId, TreeNodeId,
};

/// Represents a node within a `ParentChildStore`.
///
/// Contains the parent pointer and data (`PPNode`), an optional pointer to the
/// *first* child, and pointers to the left and right siblings. The sibling
/// pointers form a cyclic doubly linked list among all children of the same parent.
/// If a node is the only child, its `neighbor_left` and `neighbor_right` point to itself.
#[derive(Clone, Debug, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub struct PCNode<V> {
    /// Parent pointer and node data.
    pub(crate) parent_pointer: PPNode<V>,
    /// Pointer to the first child of this node, if any.
    pub(crate) child: Option<TreeNodeId>, // To get the other children, follow the sibling links of this child
    /// Pointer to the previous (left) sibling in the cyclic list. Points to self if only child.
    pub(crate) neighbor_left: TreeNodeId, // Will form a cyclic linked list, if equal to right, then it is the only child
    /// Pointer to the next (right) sibling in the cyclic list. Points to self if only child.
    pub(crate) neighbor_right: TreeNodeId, // previous sibling, if any. Will form a cyclic linked list, if equal to left, then it is the only child
}

impl<V> PCNode<V> {
    pub fn debug_display(
        &self,
        writer: &mut impl std::fmt::Write,
        draw_data: &mut impl FnMut(Option<&V>) -> Option<String>,
    ) {
        write!(writer, "{}<- ->{}", self.neighbor_left, self.neighbor_right).unwrap();
        if let Some(child) = self.child {
            write!(writer, "child: {child}").unwrap();
        }
        self.parent_pointer.debug_display(writer, draw_data);
    }

    pub fn shift(&mut self, root: RootId, nodes: TreeNodeId) {
        self.neighbor_left.0 += nodes.0;
        self.neighbor_right.0 += nodes.0;
        self.parent_pointer.shift(root, nodes);
        if let Some(child) = &mut self.child {
            child.0 += nodes.0;
        }
    }

    pub fn forget<U>(&self) -> PCNode<U> {
        PCNode {
            parent_pointer: self.parent_pointer.forget(),
            child: self.child,
            neighbor_left: self.neighbor_left,
            neighbor_right: self.neighbor_right,
        }
    }

    pub fn map<F, U>(self, transform: F) -> PCNode<U>
    where
        F: FnMut(V) -> U,
    {
        PCNode {
            parent_pointer: self.parent_pointer.map(transform),
            child: self.child,
            neighbor_left: self.neighbor_left,
            neighbor_right: self.neighbor_right,
        }
    }

    pub fn map_ref<F, U>(&self, transform: F) -> PCNode<U>
    where
        F: FnMut(&V) -> U,
    {
        PCNode {
            parent_pointer: self.parent_pointer.map_ref(transform),
            child: self.child,
            neighbor_right: self.neighbor_right,
            neighbor_left: self.neighbor_left,
        }
    }
}

/// A forest data structure using a first-child, next-sibling representation with cyclic sibling links.
///
/// Each node stores its parent, data, a pointer to its first child, and pointers to its left and
/// right siblings. Siblings under a parent form a cyclic doubly linked list.
///
/// This representation allows for efficient iteration over children and siblings.
/// It implements `ForestNodeStore` and `ForestNodeStoreDown`.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub struct ParentChildStore<V> {
    /// The flat list of nodes. The index in the vector corresponds to the `TreeNodeId`.
    pub(crate) nodes: Vec<PCNode<V>>,
}

impl<V> std::ops::Index<&TreeNodeId> for ParentChildStore<V> {
    type Output = ParentId;
    fn index(&self, index: &TreeNodeId) -> &Self::Output {
        &self.nodes[index.0].parent_pointer.parent
    }
}

impl<V> std::ops::Index<TreeNodeId> for ParentChildStore<V> {
    type Output = Option<V>;
    fn index(&self, index: TreeNodeId) -> &Self::Output {
        &self.nodes[index.0].parent_pointer.data
    }
}

impl<V> std::ops::IndexMut<TreeNodeId> for ParentChildStore<V> {
    fn index_mut(&mut self, index: TreeNodeId) -> &mut Self::Output {
        &mut self.nodes[index.0].parent_pointer.data
    }
}

impl<V> FromIterator<PCNode<V>> for ParentChildStore<V> {
    fn from_iter<I: IntoIterator<Item = PCNode<V>>>(iter: I) -> Self {
        ParentChildStore {
            nodes: iter.into_iter().collect(),
        }
    }
}

/// An iterator over the siblings of a node (including the node itself).
/// Uses the cyclic `neighbor_right` pointers.
pub enum NeighborIter<'a, V> {
    /// Case where the node has no children/siblings to iterate over.
    None,
    Some {
        initial: TreeNodeId,
        current: Option<TreeNodeId>,
        first: bool,
        store: &'a ParentChildStore<V>,
    },
}

impl<V> Iterator for NeighborIter<'_, V> {
    type Item = TreeNodeId;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            NeighborIter::None => None,
            NeighborIter::Some {
                initial,
                current,
                first,
                store,
            } => {
                let next = *current;
                if let Some(c) = next {
                    if c == *initial && !*first {
                        return None;
                    }
                    *first = false;
                    *current = store.right_neighbor_cyclic(c);
                }
                next
            }
        }
    }
}

impl<V> ParentChildStore<V> {
    /// Returns the left sibling of `node_id`.
    /// Returns `None` if `node_id` is the only child of its parent (based on cyclic link).
    pub fn left_neighbor_cyclic(&self, node_id: TreeNodeId) -> Option<TreeNodeId> {
        if self.nodes[node_id.0].neighbor_left == node_id {
            None
        } else {
            Some(self.nodes[node_id.0].neighbor_left)
        }
    }

    /// Returns the left sibling of `node_id`, stopping at the "first" child.
    ///
    /// This differs from `left_neighbor_cyclic` by returning `None` if moving left
    /// would wrap around from the *parent's first recorded child* back to the last child.
    /// Useful if a non-cyclic view of siblings is needed, although the internal
    /// representation *is* cyclic.
    /// Returns `None` if the node is a root or the only child.
    pub fn left_neighbor_ending(&self, node_id: TreeNodeId) -> Option<TreeNodeId> {
        let left = self.nodes[node_id.0].neighbor_left;
        if left == node_id {
            None
        } else {
            let first_neighbor = match self[&node_id] {
                ParentId::Node(p) => p,
                ParentId::Root(_) => return None,
                ParentId::PointingRoot(p) => p,
            };
            if left == first_neighbor {
                None
            } else {
                Some(left)
            }
        }
    }

    /// Returns the first child of `node_id`, if one exists.
    pub fn first_child(&self, node_id: TreeNodeId) -> Option<TreeNodeId> {
        if node_id.is_empty() {
            return None;
        }
        self.nodes[node_id.0].child
    }

    /// Returns the right sibling of `node_id`.
    /// Returns `None` if `node_id` is the only child of its parent (based on cyclic link).
    pub fn right_neighbor_cyclic(&self, node_id: TreeNodeId) -> Option<TreeNodeId> {
        if self.nodes[node_id.0].neighbor_right == node_id {
            None
        } else {
            Some(self.nodes[node_id.0].neighbor_right)
        }
    }

    /// Returns the right sibling of `node_id`, stopping before wrapping around the cycle.
    ///
    /// This differs from `right_neighbor_cyclic` by returning `None` if moving right
    /// would wrap around from the *last* child back to the *parent's first recorded child*.
    /// Useful if a non-cyclic view of siblings is needed.
    /// Returns `None` if the node is a root or the only child.
    pub fn right_neighbor_ending(&self, node_id: TreeNodeId) -> Option<TreeNodeId> {
        let right = self.nodes[node_id.0].neighbor_right;
        if right == node_id {
            None
        } else {
            let first_neighbor = match self[&node_id] {
                ParentId::Node(p) => p,
                ParentId::Root(_) => return None,
                ParentId::PointingRoot(p) => p,
            };
            if right == first_neighbor {
                None
            } else {
                Some(right)
            }
        }
    }
}

impl<V> ForestNodeStoreDown for ParentChildStore<V> {
    fn iter_leaves(&self) -> impl Iterator<Item = TreeNodeId> {
        self.nodes.iter().enumerate().filter_map(|(i, a)| {
            if a.child.is_none() {
                Some(TreeNodeId(i))
            } else {
                None
            }
        })
    }

    fn iter_children(&self, node_id: TreeNodeId) -> impl Iterator<Item = TreeNodeId> {
        if let Some(c) = self.first_child(node_id) {
            NeighborIter::Some {
                initial: c,
                current: Some(c),
                store: self,
                first: true,
            }
        } else {
            NeighborIter::None
        }
    }
}

impl<V> ForestNodeStore for ParentChildStore<V> {
    type NodeData = V;
    type Store<T> = ParentChildStore<T>;

    fn from_parents(nodes: impl IntoIterator<Item = PPNode<V>>) -> Self {
        let mut nodes: Vec<_> = nodes
            .into_iter()
            .enumerate()
            .map(|(index, parent_pointer)| {
                let id = TreeNodeId(index);
                PCNode {
                    parent_pointer,
                    child: None,
                    neighbor_left: id,
                    neighbor_right: id,
                }
            })
            .collect();

        // Parent-pointer conversion orders siblings by their original node ID.
        // Link that same cyclic list directly, without per-parent child vectors
        // or moving the payload through an intermediate forest store.
        for index in 0..nodes.len() {
            if let ParentId::Node(parent) = nodes[index].parent_pointer.parent {
                let id = TreeNodeId(index);
                if let Some(first) = nodes[parent.0].child {
                    let last = nodes[first.0].neighbor_left;
                    nodes[last.0].neighbor_right = id;
                    nodes[index].neighbor_left = last;
                    nodes[index].neighbor_right = first;
                    nodes[first.0].neighbor_left = id;
                } else {
                    nodes[parent.0].child = Some(id);
                }
            }
        }
        Self { nodes }
    }

    fn validate(&self) -> Result<Vec<(RootId, TreeNodeId)>, super::ForestError> {
        let mut roots = vec![];
        let mut seen: SubSet<TreeNodeId> = SubSet::full(self.n_nodes());

        while let Some(next) = seen.included_iter().next() {
            let current = next;
            seen.sub(next);
            let mut children: SubSet<TreeNodeId> = SubSet::empty(self.n_nodes());
            for p in self.iter_ancestors(current) {
                seen.sub(p);
                if children.includes(&p) {
                    return Err(ForestError::CyclicPP);
                }
                children.add(p);
            }
            for c in self.iter_children(current) {
                if let ParentId::Node(n) = self[&c] {
                    if n != current {
                        return Err(ForestError::WrongPP(current));
                    }
                }
            }
            let root_id = self.root_node(current);
            let root = self.root(root_id);

            roots.push((root, root_id));
        }

        Ok(roots)
    }

    fn debug_draw(
        &self,
        mut node_display: impl FnMut(Option<&Self::NodeData>) -> Option<String>,
    ) -> String {
        // fn debug_draw(&self, node_display: impl FnMut()) -> String {
        let mut output = String::new();

        // Recursive helper function - structure remains the same, handles sub-trees
        fn draw_subtree_recursive<W: Write, V>(
            f: &mut W,                                                  // Output writer
            nodes: &ParentChildStore<V>,                                // Node store access
            node_id: TreeNodeId,                                        // Current node to draw
            prefix: &str,        // Current line prefix (e.g., "  │   ")
            is_last_child: bool, // Is this the last child of its parent?
            format_node: &mut impl FnMut(Option<&V>) -> Option<String>, // Mutable reference to node formatter closure
            seen: &mut SubSet<TreeNodeId>,
        ) -> Result<(), std::fmt::Error> {
            // Determine the connector shape based on whether it's the last child
            let connector = if is_last_child {
                "└── "
            } else {
                "├── "
            };

            write!(f, "{prefix}{connector}{node_id}:")?;
            nodes.nodes[node_id.0].debug_display(f, format_node);

            seen.sub(node_id);
            // Prepare the prefix for the children of this node
            let child_prefix = format!("{}{}", prefix, if is_last_child { "    " } else { "│   " });

            let children: Vec<_> = nodes.iter_children(node_id).collect();
            let num_children = children.len();
            for (i, child_id) in children.into_iter().enumerate() {
                seen.add(child_id);
                // Recurse: pass the new prefix and whether this child is the last one
                draw_subtree_recursive(
                    f,
                    nodes,
                    child_id,
                    &child_prefix,
                    i == num_children - 1,
                    format_node,
                    seen,
                )?;
            }
            Ok(())
        } // End of helper fn definition

        let mut seen: SubSet<TreeNodeId> = SubSet::full(self.n_nodes());

        while let Some(next) = seen.included_iter().next() {
            let root_id = self.root_node(next);
            // println!("{root_id}");
            seen.sub(next);

            seen.sub(root_id);

            let node_line_prefix = "  ";

            write!(output, "{node_line_prefix}{root_id}:").unwrap();
            self.nodes[root_id.0].debug_display(&mut output, &mut node_display);

            let children: Vec<_> = self
                .iter_children(root_id)
                .inspect(|a| {
                    seen.sub(*a);
                })
                .collect();
            let num_children = children.len();

            for (i, child_id) in children.into_iter().enumerate() {
                // Pass the prefix that aligns the connectors (└──, ├──)
                // under the root node line.
                seen.sub(child_id);
                let _ = draw_subtree_recursive(
                    &mut output,
                    self,
                    child_id,
                    node_line_prefix, // Prefix aligns the connector itself
                    i == num_children - 1,
                    &mut node_display,
                    &mut seen,
                );
            }

            let _ = writeln!(output);
        }

        output
    }
    fn set_root(&mut self, a: TreeNodeId, root: RootId) {
        let rootx = self.root_node(a);
        self.nodes[rootx.0].parent_pointer.parent = ParentId::Root(root);
    }

    fn forgetful_map<U>(&self) -> Self::Store<U> {
        ParentChildStore {
            nodes: self.nodes.iter().map(PCNode::forget).collect(),
        }
    }

    fn n_nodes(&self) -> usize {
        self.nodes.len()
    }

    fn split_off(&mut self, at: TreeNodeId) -> Self {
        let tmp = ParentChildStore {
            nodes: std::mem::take(&mut self.nodes),
        };
        let mut as_vec: ChildVecStore<V> = tmp.into();
        let extracted = as_vec.split_off(at);
        *self = ParentChildStore::from(as_vec);
        ParentChildStore::from(extracted)
    }

    fn extend(&mut self, other: Self, shift_roots_by: RootId) {
        let nodeshift = TreeNodeId(self.nodes.len());
        self.nodes.extend(other.nodes.into_iter().map(|mut r| {
            r.shift(shift_roots_by, nodeshift);
            r
        }));
    }

    fn set_node_data(&mut self, data: Self::NodeData, node_id: TreeNodeId) -> Option<V> {
        let old = self.nodes[node_id.0].parent_pointer.data.take();
        self.nodes[node_id.0].parent_pointer.data = Some(data);
        old
    }

    fn reparent(&mut self, parent: TreeNodeId, child: TreeNodeId) {
        let parent_node = &mut self.nodes[parent.0];
        if let Some(first_child_id) = parent_node.child {
            // Parent already has children. Insert new child *after* the current last child.
            let last_child_id = self.nodes[first_child_id.0].neighbor_left; // Get current last child (left of first)

            self.nodes[last_child_id.0].neighbor_right = child;
            self.nodes[first_child_id.0].neighbor_left = child;
            self.nodes[child.0].neighbor_left = last_child_id;
            self.nodes[child.0].neighbor_right = first_child_id;
        } else {
            parent_node.child = Some(child);

            self.nodes[child.0].neighbor_left = child;
            self.nodes[child.0].neighbor_right = child;
        }

        self.nodes[child.0].parent_pointer.parent = ParentId::Node(parent);
    }

    fn swap(&mut self, a: TreeNodeId, b: TreeNodeId) {
        if a == b {
            return;
        }

        // Apply the same pointer renaming as ChildVecStore's swap without converting
        // the whole store. Only these nodes can contain references to a or b.
        let mut affected = vec![a, b];
        for id in [a, b] {
            affected.extend(self.iter_children(id));
            let node = &self.nodes[id.0];
            affected.extend([node.neighbor_left, node.neighbor_right]);
            if let ParentId::Node(parent) = node.parent_pointer.parent {
                affected.push(parent);
            }
        }
        affected.sort_unstable();
        affected.dedup();

        let rename = |id: &mut TreeNodeId| {
            if *id == a {
                *id = b;
            } else if *id == b {
                *id = a;
            }
        };
        for id in affected {
            let node = &mut self.nodes[id.0];
            if let ParentId::Node(parent) = &mut node.parent_pointer.parent {
                rename(parent);
            }
            if let Some(child) = &mut node.child {
                rename(child);
            }
            rename(&mut node.neighbor_left);
            rename(&mut node.neighbor_right);
        }
        self.nodes.swap(a.0, b.0);
    }

    fn to_store(self) -> Self::Store<Self::NodeData> {
        self
    }

    fn from_store(store: Self::Store<Self::NodeData>) -> Self {
        store
    }

    fn iter_nodes(&self) -> impl Iterator<Item = (TreeNodeId, Option<&Self::NodeData>)> {
        self.nodes
            .iter()
            .enumerate()
            .map(|(i, n)| (TreeNodeId(i), n.parent_pointer.data.as_ref()))
    }

    fn add_root(&mut self, data: Self::NodeData, root_id: super::RootId) -> TreeNodeId {
        let node_id = TreeNodeId(self.nodes.len());
        self.nodes.push(PCNode {
            parent_pointer: PPNode::root(data, root_id),
            child: None,
            neighbor_left: node_id,
            neighbor_right: node_id,
        });
        node_id
    }

    fn add_dataless_root(&mut self, root_id: super::RootId) -> TreeNodeId {
        let node_id = TreeNodeId(self.nodes.len());
        self.nodes.push(PCNode {
            parent_pointer: PPNode::dataless_root(root_id),
            child: None,
            neighbor_left: node_id,
            neighbor_right: node_id,
        });
        node_id
    }

    fn add_dataless_child(&mut self, parent: TreeNodeId) -> TreeNodeId {
        let child = TreeNodeId(self.nodes.len());
        debug_assert!(parent.0 < self.nodes.len(), "Parent ID out of bounds");

        self.nodes.push(PCNode {
            parent_pointer: PPNode::dataless_child(parent),
            child: None,
            neighbor_left: child,
            neighbor_right: child,
        });

        self.reparent(parent, child);

        child
    }

    fn add_child(&mut self, data: Self::NodeData, parent: TreeNodeId) -> TreeNodeId {
        let child = TreeNodeId(self.nodes.len());
        debug_assert!(parent.0 < self.nodes.len(), "Parent ID out of bounds");

        self.nodes.push(PCNode {
            parent_pointer: PPNode::child(data, parent),
            child: None, // New node initially has no children
            neighbor_left: child,
            neighbor_right: child,
        });
        self.reparent(parent, child);

        child
    }

    fn map<F, U>(self, mut transform: F) -> Self::Store<U>
    where
        F: FnMut(Self::NodeData) -> U,
    {
        self.nodes
            .into_iter()
            .map(|node| node.map(&mut transform))
            .collect()
    }

    fn map_ref<F, U>(&self, mut transform: F) -> Self::Store<U>
    where
        F: FnMut(&Self::NodeData) -> U,
    {
        self.nodes
            .iter()
            .map(|node| node.map_ref(&mut transform))
            .collect()
    }
}

// Preorder Iterator specific to ParentChildStore internal structure
mod pc_preorder {
    use super::*;
    pub struct PreorderIter<'a, V> {
        store: &'a ParentChildStore<V>,
        current: Option<TreeNodeId>,
        start: TreeNodeId,
        first_sibling: TreeNodeId,
    }

    impl<V> Clone for PreorderIter<'_, V> {
        fn clone(&self) -> Self {
            *self
        }
    }

    impl<V> Copy for PreorderIter<'_, V> {}

    impl<'a, V> PreorderIter<'a, V> {
        pub fn new(store: &'a ParentChildStore<V>, start: TreeNodeId) -> Self {
            PreorderIter {
                store,
                current: (!start.is_empty()).then_some(start),
                start,
                first_sibling: start,
            }
        }

        // Cache the first sibling so advancing along a crown does not repeatedly
        // load the parent and its child pointer.
        fn next_node_logic(&mut self, node: TreeNodeId) -> Option<TreeNodeId> {
            if let Some(child) = self.store.nodes[node.0].child {
                self.first_sibling = child;
                return Some(child);
            }
            let mut current = node;
            while current != self.start {
                let entry = &self.store.nodes[current.0];
                if entry.neighbor_right != self.first_sibling {
                    return Some(entry.neighbor_right);
                }
                current = match entry.parent_pointer.parent {
                    ParentId::Node(parent) | ParentId::PointingRoot(parent) => parent,
                    ParentId::Root(_) => return None,
                };
                if current != self.start {
                    self.first_sibling = match self.store[&current] {
                        ParentId::Node(parent) | ParentId::PointingRoot(parent) => {
                            self.store.nodes[parent.0].child.unwrap()
                        }
                        ParentId::Root(_) => return None,
                    };
                }
            }
            None
        }
    }

    impl<V> Iterator for PreorderIter<'_, V> {
        type Item = TreeNodeId;
        fn next(&mut self) -> Option<Self::Item> {
            let node_to_return = self.current?;
            self.current = self.next_node_logic(node_to_return);
            Some(node_to_return)
        }
    }
}

// Bfs Iterator specific to ParentChildStore
mod pc_bfs {
    use super::*;
    #[derive(Clone)]
    pub struct BfsTreeIter<'a, V> {
        store: &'a ParentChildStore<V>,
        queue: VecDeque<TreeNodeId>,
    }
    impl<'a, V> BfsTreeIter<'a, V> {
        pub fn new(store: &'a ParentChildStore<V>, start: TreeNodeId) -> Self {
            let mut queue = VecDeque::new();
            if !start.is_empty() {
                queue.push_back(start);
            }
            BfsTreeIter { store, queue }
        }
    }
    impl<V> Iterator for BfsTreeIter<'_, V> {
        type Item = TreeNodeId;
        fn next(&mut self) -> Option<Self::Item> {
            let node = self.queue.pop_front()?;
            // Use the store's iter_children which correctly uses NeighborIter
            for child in self.store.iter_children(node) {
                self.queue.push_back(child);
            }
            Some(node)
        }
    }
}

impl<V> ForestNodeStorePreorder for ParentChildStore<V> {
    type Iterator<'a>
        = pc_preorder::PreorderIter<'a, V>
    where
        Self: 'a;
    fn iter_preorder(&self, start: TreeNodeId) -> Self::Iterator<'_> {
        pc_preorder::PreorderIter::new(self, start)
    }

    fn flatten_branches(&mut self, nodes: impl IntoIterator<Item = TreeNodeId>) {
        for node in nodes {
            self.flatten_branch(node);
        }
    }

    #[inline]
    fn flatten_branch(&mut self, node: TreeNodeId) -> TreeNodeId {
        let Some((root, branch)) = self.branch(node) else {
            return node;
        };
        if self.nodes[branch.0].child.is_none() {
            return root;
        }
        let stop = self.nodes[branch.0].neighbor_right;
        let mut current = branch;
        loop {
            // Splice each child ring immediately after its parent, retaining
            // preorder without a stack or an intermediate incidence vector.
            if let Some(child) = self.nodes[current.0].child.take() {
                let next = self.nodes[current.0].neighbor_right;
                let last = self.nodes[child.0].neighbor_left;
                self.nodes[current.0].neighbor_right = child;
                self.nodes[child.0].neighbor_left = current;
                self.nodes[last.0].neighbor_right = next;
                self.nodes[next.0].neighbor_left = last;
            }
            self.nodes[current.0].parent_pointer.parent = ParentId::Node(root);
            current = self.nodes[current.0].neighbor_right;
            if current == stop {
                break;
            }
        }
        root
    }
}

impl<V> ForestNodeStoreBfs for ParentChildStore<V> {
    fn iter_bfs(&self, start: TreeNodeId) -> impl Iterator<Item = TreeNodeId> + '_ {
        pc_bfs::BfsTreeIter::new(self, start)
    }
}

impl<V> From<ParentChildStore<V>> for ParentPointerStore<V> {
    fn from(child_store: ParentChildStore<V>) -> Self {
        child_store
            .nodes
            .into_iter()
            .map(|pc_node| pc_node.parent_pointer)
            .collect()
    }
}

impl<V> From<ParentPointerStore<V>> for ParentChildStore<V> {
    fn from(value: ParentPointerStore<V>) -> Self {
        Self::from_parents(value.nodes)
    }
}

#[cfg(test)]
mod test {
    use crate::{
        half_edge::{
            involution::Hedge,
            nodestore::NodeStorageOps,
            subgraph::{ModifySubSet, SuBitGraph, SubSetLike},
        },
        tree::{
            child_vec::ChildVecStore, Forest, ForestNodeStore, ForestNodeStoreDown, RootId,
            TreeNodeId,
        },
    };
    use std::collections::BTreeSet;

    use super::ParentChildStore;

    #[test]
    fn parent_pointer_conversion_preserves_all_node_fields() {
        use crate::tree::parent_pointer::{PPNode, ParentId, ParentPointerStore};

        // Include forward parents, cycles, self-parents, dataless nodes and
        // pointing roots. Conversion preserves these raw fields; validation of
        // a tree's shape remains the existing separate boundary.
        for count in 0..=5usize {
            let choices = count + 4;
            for mut code in 0..choices.pow(count as u32) {
                let source: ParentPointerStore<_> = (0..count)
                    .map(|index| {
                        let choice = code % choices;
                        code /= choices;
                        let parent = match choice.checked_sub(count) {
                            None => ParentId::Node(TreeNodeId(choice)),
                            Some(0) => ParentId::Root(RootId(0)),
                            Some(1) => ParentId::Root(RootId(usize::MAX)),
                            Some(2) => ParentId::PointingRoot(TreeNodeId(0)),
                            _ => ParentId::PointingRoot(TreeNodeId(usize::MAX)),
                        };
                        PPNode {
                            parent,
                            data: (index % 2 == 0).then_some(index),
                        }
                    })
                    .collect();
                let expected: ParentChildStore<_> = ChildVecStore::from(source.clone()).into();
                let actual = ParentChildStore::from(source);
                assert_eq!(actual.nodes, expected.nodes);
            }
        }
    }

    #[test]
    fn parent_pointer_conversion_preserves_nonclone_payloads_and_invalid_parent_drop() {
        use crate::tree::parent_pointer::{PPNode, ParentId, ParentPointerStore};
        use std::{cell::RefCell, panic::AssertUnwindSafe, rc::Rc};

        struct Payload(usize, Rc<RefCell<Vec<usize>>>);
        impl Drop for Payload {
            fn drop(&mut self) {
                self.1.borrow_mut().push(self.0);
            }
        }
        for parent in [0, 1, 2, 3, usize::MAX] {
            let mut results = Vec::new();
            for direct in [false, true] {
                let drops = Rc::new(RefCell::new(Vec::new()));
                let source: ParentPointerStore<_> = (0..3)
                    .map(|index| PPNode {
                        parent: if index == 1 {
                            ParentId::Node(TreeNodeId(parent))
                        } else {
                            ParentId::Root(RootId(index))
                        },
                        data: Some(Payload(index, drops.clone())),
                    })
                    .collect();
                let outcome = std::panic::catch_unwind(AssertUnwindSafe(|| {
                    let converted: ParentChildStore<_> = if direct {
                        ParentChildStore::from(source)
                    } else {
                        ChildVecStore::from(source).into()
                    };
                    assert!(drops.borrow().is_empty());
                    drop(converted);
                }));
                results.push((outcome.is_ok(), drops.borrow().clone()));
            }
            assert_eq!(results[0], results[1], "parent {parent}");
            assert_eq!(results[0].0, parent < 3);
            assert_eq!(results[0].1, [0, 1, 2]);
        }
    }

    #[test]
    fn subtree_traversal_and_branch_flattening_preserve_preorder() {
        use crate::tree::{parent_pointer::ParentId, ForestNodeStorePreorder};

        for forest_code in 0..720 {
            let mut code = forest_code;
            let mut original = ParentChildStore { nodes: Vec::new() };
            for id in 0..6 {
                let parent = code % (id + 1);
                code /= id + 1;
                if parent == id {
                    original.add_root(id, RootId(id));
                } else {
                    original.add_child(id, TreeNodeId(parent));
                }
            }
            // Exercise both ascending and shuffled node IDs.
            for shuffled in [false, true] {
                if shuffled {
                    original.swap(TreeNodeId(0), TreeNodeId(5));
                    original.swap(TreeNodeId(1), TreeNodeId(3));
                }
                let oracle = ChildVecStore::from(original.clone());
                for start in (0..6).map(TreeNodeId).chain([TreeNodeId::EMPTY]) {
                    let expected: Vec<_> = oracle.iter_preorder(start).collect();
                    let mut iter = original.iter_preorder(start);
                    for offset in 0..=expected.len() {
                        assert_eq!(iter.collect::<Vec<_>>(), expected[offset..]);
                        assert_eq!(iter.next(), expected.get(offset).copied());
                    }
                    assert_eq!(iter.next(), None);
                    if start.is_empty() {
                        continue;
                    }

                    // Flatten the branch containing `start`, then all branches at once.
                    for nodes in [vec![start], (0..6).map(TreeNodeId).collect()] {
                        let mut actual = original.clone();
                        let mut flattened_oracle = oracle.clone();
                        actual.flatten_branches(nodes.clone());
                        flattened_oracle.flatten_branches(nodes.clone());
                        assert_eq!(ChildVecStore::from(actual.clone()), flattened_oracle);
                        // Single-branch calls return the root. The parent/child
                        // store flattens just that branch; the child-vector store
                        // flattens the whole tree once the branch is not flat.
                        let (mut pc, mut cv, mut cv_expected) =
                            (original.clone(), oracle.clone(), oracle.clone());
                        for &node in &nodes {
                            let root = original.root_node(node);
                            assert_eq!(pc.flatten_branch(node), root);
                            assert_eq!(cv.flatten_branch(node), root);
                            if cv_expected.branch(node).is_some_and(|(_, branch)| {
                                cv_expected.iter_children(branch).next().is_some()
                            }) {
                                let tree = original.iter_node_id();
                                cv_expected.flatten_branches(
                                    tree.filter(|&n| original.root_node(n) == root),
                                );
                            }
                        }
                        assert_eq!(pc, actual);
                        assert_eq!(cv, cv_expected);
                        actual.validate().unwrap();
                        let branches: Vec<_> =
                            nodes.iter().filter_map(|&n| original.branch(n)).collect();
                        for (id, data) in original.iter_nodes() {
                            assert_eq!(actual[id].as_ref(), data);
                            let root = original.root_node(id);
                            assert_eq!(
                                actual.iter_preorder(root).collect::<Vec<_>>(),
                                oracle.iter_preorder(root).collect::<Vec<_>>()
                            );
                            if original.branch(id).is_some_and(|b| branches.contains(&b)) {
                                assert_eq!(actual.parent(id), ParentId::Node(root));
                                assert_eq!(actual.iter_children(id).count(), 0);
                            } else {
                                assert_eq!(actual.parent(id), original.parent(id));
                            }
                        }
                        let once = actual.clone();
                        actual.flatten_branches(nodes);
                        assert_eq!(actual, once);
                    }
                }
            }
        }
    }

    #[test]
    fn swaps_match_child_vec_for_all_small_forests() {
        // Each node is either a new root or a child of an earlier node. Enumerate
        // all such six-node forests, then also change their index order serially.
        for forest_code in 0..720 {
            let mut code = forest_code;
            let mut original = ParentChildStore { nodes: Vec::new() };
            let mut roots = 0;
            for id in 0..6 {
                let parent = code % (id + 1);
                code /= id + 1;
                if parent == id {
                    original.add_root(id, RootId(roots));
                    roots += 1;
                } else {
                    original.add_child(id, TreeNodeId(parent));
                }
            }
            let mut serial = original.clone();
            let mut serial_oracle = ChildVecStore::from(original.clone());
            for a in 0..6 {
                for b in 0..6 {
                    let a = TreeNodeId(a);
                    let b = TreeNodeId(b);
                    let mut actual = original.clone();
                    let mut oracle = ChildVecStore::from(original.clone());
                    actual.swap(a, b);
                    oracle.swap(a, b);
                    assert_eq!(
                        actual
                            .iter_nodes()
                            .map(|(node, _)| node)
                            .collect::<BTreeSet<_>>(),
                        oracle
                            .iter_nodes()
                            .map(|(node, _)| node)
                            .collect::<BTreeSet<_>>(),
                    );
                    for (node, data) in actual.iter_nodes() {
                        assert_eq!(data, oracle[node].as_ref());
                        assert_eq!(actual.parent(node), oracle.parent(node));
                        assert_eq!(
                            actual.iter_children(node).collect::<BTreeSet<_>>(),
                            oracle.iter_children(node).collect::<BTreeSet<_>>(),
                        );
                    }
                    actual.validate().unwrap();

                    serial.swap(a, b);
                    serial_oracle.swap(a, b);
                    assert_eq!(
                        serial
                            .iter_nodes()
                            .map(|(node, _)| node)
                            .collect::<BTreeSet<_>>(),
                        serial_oracle
                            .iter_nodes()
                            .map(|(node, _)| node)
                            .collect::<BTreeSet<_>>(),
                    );
                    for (node, data) in serial.iter_nodes() {
                        assert_eq!(data, serial_oracle[node].as_ref());
                        assert_eq!(serial.parent(node), serial_oracle.parent(node));
                        assert_eq!(
                            serial.iter_children(node).collect::<BTreeSet<_>>(),
                            serial_oracle.iter_children(node).collect::<BTreeSet<_>>(),
                        );
                    }
                    serial.validate().unwrap();
                }
            }
        }
    }

    #[test]
    fn deletion_matches_child_vec_for_every_subset() {
        let mut forest: Forest<usize, ParentChildStore<usize>> = Forest::new();
        let (root, _) = forest.add_root(0, 0);
        let first = forest.add_child(root, 1);
        forest.add_child(root, 2);
        forest.add_child(first, 3);
        forest.add_child(first, 4);
        let (other, _) = forest.add_root(5, 1);
        forest.add_child(other, 6);
        forest.add_child(other, 7);
        forest.add_root(8, 2);

        // Include empty/full deletion, deleted roots, split trees and isolated roots.
        for mask in 0..512 {
            let mut deleted = SuBitGraph::empty(9);
            for id in 0..9 {
                if mask & (1 << id) != 0 {
                    deleted.add(Hedge(id));
                }
            }
            let mut actual = forest.clone();
            let mut oracle = forest.clone().cast::<ChildVecStore<usize>>();
            actual.delete(&deleted);
            oracle.delete(&deleted);
            let actual_values = actual
                .iter_nodes()
                .map(|(node, data)| {
                    (
                        data.copied(),
                        actual[actual.root(node)],
                        actual
                            .iter_children(node)
                            .map(|child| actual[child])
                            .collect::<BTreeSet<_>>(),
                    )
                })
                .collect::<BTreeSet<_>>();
            let expected_values = oracle
                .iter_nodes()
                .map(|(node, data)| {
                    (
                        data.copied(),
                        oracle[oracle.root(node)],
                        oracle
                            .iter_children(node)
                            .map(|child| oracle[child])
                            .collect::<BTreeSet<_>>(),
                    )
                })
                .collect::<BTreeSet<_>>();
            assert_eq!(actual_values, expected_values);
            actual.validate_structure().unwrap();
        }
    }
}
