//! Implements a tree structure where each node only stores a pointer to its parent.

use std::{
    fmt::Display,
    ops::{Index, IndexMut},
};

use crate::half_edge::subgraph::{subset::SubSet, ModifySubSet, SubSetLike};

use super::{Forest, ForestNodeStore, ForestNodeStoreAncestors, RootId, TreeNodeId};

/// Represents a node within a `ParentPointerStore`.
///
/// Contains the actual data (`V`) and a [ParentId] which points
/// either to the parent [TreeNodeId] or identifies this node as a root via `RootId`.
#[derive(Clone, Debug, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
pub struct PPNode<V> {
    /// Pointer to the parent node or the root ID if this is a root node.
    pub(crate) parent: ParentId,
    /// The data associated with this node.
    pub(crate) data: Option<V>,
}

impl<V> PPNode<V> {
    pub fn child(data: V, parent: TreeNodeId) -> Self {
        PPNode {
            parent: ParentId::Node(parent),
            data: Some(data),
        }
    }

    pub fn debug_display(
        &self,
        writer: &mut impl std::fmt::Write,
        draw_data: &mut impl FnMut(Option<&V>) -> Option<String>,
    ) {
        write!(writer, "{}", self.parent).unwrap();
        if let Some(root_data) = draw_data(self.data.as_ref()) {
            writeln!(writer, "{root_data}").unwrap();
        } else {
            writeln!(writer).unwrap();
        }
    }

    pub fn forget<U>(&self) -> PPNode<U> {
        PPNode {
            parent: self.parent,
            data: None,
        }
    }

    pub fn shift(&mut self, root: RootId, nodes: TreeNodeId) {
        match &mut self.parent {
            ParentId::Node(n) => {
                n.0 += nodes.0;
            }
            ParentId::PointingRoot(n) => {
                n.0 += nodes.0;
            }
            ParentId::Root(r) => {
                r.0 += root.0;
            }
        }
    }

    pub fn dataless_child(parent: TreeNodeId) -> Self {
        PPNode {
            parent: ParentId::Node(parent),
            data: None,
        }
    }

    pub fn root(data: V, root_id: RootId) -> Self {
        PPNode {
            parent: ParentId::Root(root_id),
            data: Some(data),
        }
    }

    pub fn dataless_root(root_id: RootId) -> Self {
        PPNode {
            parent: ParentId::Root(root_id),
            data: None,
        }
    }

    pub fn map<F, U>(self, transform: F) -> PPNode<U>
    where
        F: FnMut(V) -> U,
    {
        PPNode {
            parent: self.parent,
            data: self.data.map(transform),
        }
    }

    pub fn map_ref<F, U>(&self, transform: F) -> PPNode<U>
    where
        F: FnMut(&V) -> U,
    {
        PPNode {
            parent: self.parent,
            data: self.data.as_ref().map(transform),
        }
    }
}

/// Identifies the parent of a node.
///
/// A node is either a child of another `Node` or it's the `Root` of a tree,
/// identified by a `RootId`.
#[derive(Clone, Debug, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
pub enum ParentId {
    /// This node is a root node belonging to the tree identified by `RootId`.
    Root(RootId),

    /// This node is a child of the node identified by `TreeNodeId`.
    Node(TreeNodeId),
    /// for use when splitting
    PointingRoot(TreeNodeId),
}

impl Display for ParentId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParentId::Node(n) => write!(f, "parent node: {n}"),
            ParentId::Root(r) => write!(f, "root id:{r}"),
            ParentId::PointingRoot(n) => write!(f, "root node:{n}"),
        }
    }
}

impl ParentId {
    pub fn is_root(&self) -> bool {
        match self {
            ParentId::Root(_) | ParentId::PointingRoot(_) => true,
            ParentId::Node(_) => false,
        }
    }

    pub fn is_node(&self) -> bool {
        !self.is_root()
    }
}

/// A forest data structure where each node only stores its data and a pointer to its parent.
///
/// This representation is memory-efficient, especially for sparse trees, and allows for
/// very fast traversal *upwards* towards the root. However, traversing *downwards*
/// (finding children) requires iterating through all nodes to find those pointing to a
/// specific parent, which can be slow.
///
/// It implements the `ForestNodeStore` trait.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
pub struct ParentPointerStore<V> {
    /// The flat list of nodes. The index in the vector corresponds to the `TreeNodeId`.
    pub(crate) nodes: Vec<PPNode<V>>,
}

impl<V, R> Forest<R, ParentPointerStore<V>> {
    /// Changes the root of the tree containing `new_root`.
    ///
    /// All parent pointers on the path from the `new_root` to the original root
    /// are reversed. The `new_root` becomes the root node associated with the
    /// original tree's `RootId`. Updates the `Forest`'s root tracking accordingly.
    ///
    /// Returns the `RootId` of the tree that was modified.
    pub fn change_to_root(&mut self, new_root: TreeNodeId) -> RootId {
        let root_id = self.nodes.change_root(new_root);
        self.roots[root_id.0].root_id = new_root;
        root_id
    }
}

impl<V> ParentPointerStore<V> {
    /// Re–roots the tree at the given node (making it a root).
    /// Along the chain from `new_root` to the old root, the parent pointers are reversed.
    /// The `new_root` becomes associated with the original `RootId`.
    ///
    /// Returns the `RootId` of the affected tree.
    ///
    /// **Note:** This modifies the store directly. If using within a `Forest`,
    /// prefer `Forest::change_to_root` which also updates the `Forest`'s root list.
    pub fn change_root(&mut self, new_root: TreeNodeId) -> RootId {
        let mut current = new_root;
        let root_id = self.root(new_root);
        let mut prev = None;

        loop {
            let orig_parent = self[&current];
            self.nodes[current.0].parent = match prev {
                None => ParentId::Root(root_id),
                Some(p) => ParentId::Node(p),
            };
            prev = Some(current);

            match orig_parent {
                ParentId::Node(p) => current = p,

                ParentId::PointingRoot(p) => current = p,
                ParentId::Root(r) => return r,
            }
        }
    }
}

impl<V> FromIterator<PPNode<V>> for ParentPointerStore<V> {
    fn from_iter<I: IntoIterator<Item = PPNode<V>>>(iter: I) -> Self {
        ParentPointerStore {
            nodes: iter.into_iter().collect(),
        }
    }
}

impl<V> Index<&TreeNodeId> for ParentPointerStore<V> {
    type Output = ParentId;
    fn index(&self, index: &TreeNodeId) -> &Self::Output {
        &self.nodes[index.0].parent
    }
}

impl<V> Index<TreeNodeId> for ParentPointerStore<V> {
    type Output = Option<V>;
    fn index(&self, index: TreeNodeId) -> &Self::Output {
        &self.nodes[index.0].data
    }
}

impl<V> IndexMut<TreeNodeId> for ParentPointerStore<V> {
    fn index_mut(&mut self, index: TreeNodeId) -> &mut Self::Output {
        &mut self.nodes[index.0].data
    }
}

impl<V> ForestNodeStore for ParentPointerStore<V> {
    type NodeData = V;
    type Store<T> = ParentPointerStore<T>;

    fn set_root(&mut self, a: TreeNodeId, root: RootId) {
        let rootx = self.root_node(a);
        self.nodes[rootx.0].parent = ParentId::Root(root);
    }

    fn validate(&self) -> Result<Vec<(RootId, TreeNodeId)>, super::ForestError> {
        let mut roots = vec![];
        let mut seen: SubSet<usize> = !SubSet::empty(self.n_nodes());

        while let Some(next) = seen.included_iter().next() {
            let current = TreeNodeId(next);
            seen.sub(next);
            self.iter_ancestors(current).for_each(|a| {
                seen.sub(a.0);
            });
            let root_id = self.root_node(current);
            let root = self.root(root_id);

            roots.push((root, root_id));
        }

        Ok(roots)
    }

    fn debug_draw(
        &self,
        _node_display: impl FnMut(Option<&Self::NodeData>) -> Option<String>,
    ) -> String {
        "".to_string()
    }
    fn n_nodes(&self) -> usize {
        self.nodes.len()
    }

    fn set_node_data(&mut self, data: Self::NodeData, node_id: TreeNodeId) -> Option<V> {
        let old = self.nodes[node_id.0].data.take();
        self.nodes[node_id.0].data = Some(data);
        old
    }

    fn forgetful_map<U>(&self) -> Self::Store<U> {
        ParentPointerStore {
            nodes: self
                .nodes
                .iter()
                .map(|a| PPNode {
                    parent: a.parent,
                    data: None,
                })
                .collect(),
        }
    }

    // with path compression for split nodes
    fn split_off(&mut self, at: TreeNodeId) -> Self {
        for mut current in (0..self.nodes.len()).map(TreeNodeId) {
            while let ParentId::Node(p) = self[&current] {
                let mut newp = p;
                while (newp <= at) ^ (current > at) {
                    // the pointer crosses the splitting boundary, we need to find a new pointer.

                    match self[&newp] {
                        ParentId::Node(next) => {
                            newp = next;
                        }
                        ParentId::Root(_) => {
                            // split_root = true;
                            break;
                        }
                        ParentId::PointingRoot(next) => {
                            newp = next;
                        }
                    }
                }

                if let ParentId::Root(r) = self[&newp] {
                    // set current to now be the root, and original root to point here.

                    self.nodes[newp.0].parent = ParentId::PointingRoot(current);
                    self.nodes[current.0].parent = ParentId::Root(r);
                } else {
                    self.reparent(newp, current);
                }

                current = p
            }
        }

        for n in (0..self.nodes.len()).map(TreeNodeId) {
            if let ParentId::Node(n) = &mut self.nodes[n.0].parent {
                if *n > at {
                    n.0 -= at.0;
                }
            }

            if let ParentId::PointingRoot(rp) = self[&n] {
                let root = self[&rp];
                self.nodes[n.0].parent = root;
            }
        }
        Self {
            nodes: self.nodes.split_off(at.0),
        }
    }

    fn extend(&mut self, other: Self, shift_roots_by: RootId) {
        let nodeshift = TreeNodeId(self.nodes.len());
        self.nodes.extend(other.nodes.into_iter().map(|mut r| {
            r.shift(shift_roots_by, nodeshift);
            r
        }));
    }

    fn reparent(&mut self, parent: TreeNodeId, child: TreeNodeId) {
        self.nodes[child.0].parent = ParentId::Node(parent);
    }

    fn swap(&mut self, a: TreeNodeId, b: TreeNodeId) {
        println!("Swap PP");
        for n in &mut self.nodes {
            if let ParentId::Node(pp) = &mut n.parent {
                if *pp == a {
                    *pp = b;
                } else if *pp == b {
                    *pp = a;
                }
            }
        }

        self.nodes.swap(a.0, b.0);
    }
    fn from_store(store: Self::Store<Self::NodeData>) -> Self {
        store
    }

    fn to_store(self) -> Self::Store<Self::NodeData> {
        self
    }

    fn add_root(&mut self, data: Self::NodeData, root_id: RootId) -> TreeNodeId {
        let node_id = TreeNodeId(self.nodes.len());
        self.nodes.push(PPNode::root(data, root_id));
        node_id
    }

    fn iter_nodes(&self) -> impl Iterator<Item = (TreeNodeId, Option<&Self::NodeData>)> {
        self.nodes
            .iter()
            .enumerate()
            .map(|(i, n)| (TreeNodeId(i), n.data.as_ref()))
    }

    fn add_child(&mut self, data: Self::NodeData, parent: TreeNodeId) -> TreeNodeId {
        let node_id = TreeNodeId(self.nodes.len());
        self.nodes.push(PPNode::child(data, parent));
        node_id
    }

    fn add_dataless_child(&mut self, parent: TreeNodeId) -> TreeNodeId {
        let node_id = TreeNodeId(self.nodes.len());
        self.nodes.push(PPNode::dataless_child(parent));
        node_id
    }

    fn map<F, U>(self, mut transform: F) -> Self::Store<U>
    where
        F: FnMut(Self::NodeData) -> U,
    {
        self.nodes
            .into_iter()
            .map(|n| n.map(&mut transform))
            .collect()
    }

    fn map_ref<F, U>(&self, mut transform: F) -> Self::Store<U>
    where
        F: FnMut(&Self::NodeData) -> U,
    {
        self.nodes
            .iter()
            .map(|n| n.map_ref(&mut transform))
            .collect()
    }
}

#[cfg(test)]
mod test {
    use crate::tree::Forest;

    use super::ParentPointerStore;

    #[test]
    fn reroot() {
        let mut tree: Forest<i8, ParentPointerStore<i8>> = Forest::new();

        let (a, _) = tree.add_root(1, 1);
        let b = tree.add_child(a, 2);
        let c = tree.add_child(b, 3);
        let d = tree.add_child(c, 4);

        assert_eq!(tree[&tree.root(d)], a);
        tree.change_to_root(d);
        assert_eq!(tree[&tree.root(a)], d);
    }
}
