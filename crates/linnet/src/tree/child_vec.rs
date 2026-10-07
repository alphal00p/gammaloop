use std::fmt::Write;

use itertools::Itertools;

use crate::half_edge::subgraph::{subset::SubSet, Inclusion, ModifySubSet, SubSetLike};

use super::{
    child_pointer::{PCNode, ParentChildStore},
    iterato::{BfsIter, PreorderIter},
    parent_pointer::{PPNode, ParentId, ParentPointerStore},
    ForestError, ForestNodeStore, ForestNodeStoreAncestors, ForestNodeStoreBfs,
    ForestNodeStoreDown, ForestNodeStorePreorder, RootId, TreeNodeId,
};

/// A node in the ChildVecStore. It contains a PPNode plus an ordered vector of children.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub struct CVNode<V> {
    pub parent_pointer: PPNode<V>,
    pub children: Vec<TreeNodeId>,
}
impl<V> CVNode<V> {
    pub fn shift(&mut self, root: RootId, nodes: TreeNodeId) {
        self.parent_pointer.shift(root, nodes);
        for c in &mut self.children {
            c.0 += nodes.0;
        }
    }

    pub fn debug_display(
        &self,
        writer: &mut impl std::fmt::Write,
        draw_data: &mut impl FnMut(Option<&V>) -> Option<String>,
    ) {
        write!(writer, "children:").unwrap();
        for child in &self.children {
            write!(writer, "{child},").unwrap();
        }
        self.parent_pointer.debug_display(writer, draw_data);
    }

    pub fn forget<U>(&self) -> CVNode<U> {
        CVNode {
            parent_pointer: self.parent_pointer.forget(),
            children: self.children.clone(),
        }
    }
}

/// The ChildVecStore itself.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub struct ChildVecStore<V> {
    pub nodes: Vec<CVNode<V>>,
}

//
// Indexing implementations, similar to ParentChildStore/ParentPointerStore.
//
impl<V> std::ops::Index<&TreeNodeId> for ChildVecStore<V> {
    type Output = ParentId;
    fn index(&self, index: &TreeNodeId) -> &Self::Output {
        &self.nodes[index.0].parent_pointer.parent
    }
}

impl<V> std::ops::Index<TreeNodeId> for ChildVecStore<V> {
    type Output = Option<V>;
    fn index(&self, index: TreeNodeId) -> &Self::Output {
        &self.nodes[index.0].parent_pointer.data
    }
}

impl<V> std::ops::IndexMut<TreeNodeId> for ChildVecStore<V> {
    fn index_mut(&mut self, index: TreeNodeId) -> &mut Self::Output {
        &mut self.nodes[index.0].parent_pointer.data
    }
}

impl<V> FromIterator<CVNode<V>> for ChildVecStore<V> {
    fn from_iter<I: IntoIterator<Item = CVNode<V>>>(iter: I) -> Self {
        ChildVecStore {
            nodes: iter.into_iter().collect(),
        }
    }
}

//
// Implement the ForestNodeStore trait for ChildVecStore:
//
impl<V> ForestNodeStore for ChildVecStore<V> {
    type NodeData = V;
    type Store<T> = ChildVecStore<T>;
    fn debug_draw(
        &self,
        mut node_display: impl FnMut(Option<&Self::NodeData>) -> Option<String>,
    ) -> String {
        // fn debug_draw(&self, node_display: impl FnMut()) -> String {
        let mut output = String::new();
        writeln!(output, "Number of nodes:{}", self.n_nodes()).unwrap();

        // Recursive helper function - structure remains the same, handles sub-trees
        fn draw_subtree_recursive<W: Write, V>(
            f: &mut W,                                                  // Output writer
            nodes: &ChildVecStore<V>,                                   // Node store access
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
                seen.sub(child_id);
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

        let mut seen: SubSet<TreeNodeId> = !SubSet::empty(self.n_nodes());

        while let Some(next) = seen.included_iter().next() {
            let root_id = self.root_node(next);
            // println!("{}", root_id);

            seen.sub(root_id);

            seen.sub(next);

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

    fn validate(&self) -> Result<Vec<(RootId, TreeNodeId)>, super::ForestError> {
        let mut roots = vec![];
        let mut seen: SubSet<TreeNodeId> = SubSet::full(self.n_nodes());

        while let Some(next) = seen.included_iter().next() {
            let current = next;
            seen.sub(next);

            if ParentId::Node(current) == self.parent(current) {
                return Err(ForestError::SelfLoopPP(current));
            }
            if ParentId::PointingRoot(current) == self.parent(current) {
                return Err(ForestError::SelfLoopPP(current));
            }
            let mut parents: SubSet<TreeNodeId> = SubSet::empty(self.n_nodes());
            for p in self.iter_ancestors(current) {
                if parents.includes(&p) {
                    return Err(ForestError::CyclicPP);
                }
                parents.add(p);
            }
            let mut children: SubSet<TreeNodeId> = SubSet::empty(self.n_nodes());

            for c in self.iter_preorder(current) {
                if children.includes(&c) {
                    return Err(ForestError::CyclicCP);
                }
                children.add(c);
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

    fn set_root(&mut self, a: TreeNodeId, root: RootId) {
        let rootx = self.root_node(a);
        self.nodes[rootx.0].parent_pointer.parent = ParentId::Root(root);
    }

    fn forgetful_map<U>(&self) -> Self::Store<U> {
        ChildVecStore {
            nodes: self.nodes.iter().map(CVNode::forget).collect(),
        }
    }

    fn n_nodes(&self) -> usize {
        self.nodes.len()
    }

    fn set_node_data(&mut self, data: Self::NodeData, node_id: TreeNodeId) -> Option<V> {
        let old = self.nodes[node_id.0].parent_pointer.data.take();
        self.nodes[node_id.0].parent_pointer.data = Some(data);
        old
    }

    fn split_off(&mut self, at: TreeNodeId) -> Self {
        // Splitting only rewrites parent pointers; children follow from them.
        let mut parents = ParentPointerStore::from(ChildVecStore {
            nodes: std::mem::take(&mut self.nodes),
        });
        let split = parents.split_off(at);
        *self = parents.into();
        split.into()
    }

    fn extend(&mut self, other: Self, shift_roots_by: RootId) {
        let nodeshift = TreeNodeId(self.nodes.len());
        self.nodes.extend(other.nodes.into_iter().map(|mut r| {
            r.shift(shift_roots_by, nodeshift);
            r
        }));
    }
    fn reparent(&mut self, parent: TreeNodeId, child: TreeNodeId) {
        // First remove child pointer of old parent:
        let old_parent = self.parent(child);
        if let ParentId::Node(old_parent) = old_parent {
            if old_parent == parent {
                return;
            }

            let pos_in_children = self
                .iter_children(old_parent)
                .find_position(|a| *a == child);
            if let Some((pos, _)) = pos_in_children {
                self.nodes[old_parent.0].children.swap_remove(pos);
            }
        }
        self.nodes[child.0].parent_pointer.parent = ParentId::Node(parent);
        self.nodes[parent.0].children.push(child);
    }

    fn swap(&mut self, a: TreeNodeId, b: TreeNodeId) {
        // println!("Swapping {a} <-> {b}");
        if a == b {
            return;
        }

        // println!("{}", self.debug_draw(|_| None));
        // println!("swapping {a} <-> {b}");
        // Before modifying the parent pointers of the children, we save them.
        let parent_a = self.parent(a);
        let parent_b = self.parent(b);

        // Now modify the parent pointers of the children.
        {
            let (bef, after) = self.nodes.split_at_mut(a.0);
            let (anode, after) = after.split_first_mut().unwrap();

            for c in &anode.children {
                if *c <= a {
                    if let Some(child) = bef.get_mut(c.0) {
                        // println!("{}", b);
                        if let ParentId::Node(pp) = &mut child.parent_pointer.parent {
                            if *pp == a {
                                *pp = b;
                            }
                        }
                    }
                } else if let Some(child) = after.get_mut(c.0 - (a.0 + 1)) {
                    if let ParentId::Node(pp) = &mut child.parent_pointer.parent {
                        if *pp == a {
                            *pp = b;
                        }
                    }
                }
            }
        }
        {
            let (bef, after) = self.nodes.split_at_mut(b.0);
            let (bnode, after) = after.split_first_mut().unwrap();

            for c in &bnode.children {
                if *c <= b {
                    if let Some(child) = bef.get_mut(c.0) {
                        if let ParentId::Node(pp) = &mut child.parent_pointer.parent {
                            if *pp == b {
                                *pp = a;
                            }
                        }
                    }
                } else if let Some(child) = after.get_mut(c.0 - b.0 - 1) {
                    if let ParentId::Node(pp) = &mut child.parent_pointer.parent {
                        if *pp == b {
                            *pp = a;
                        }
                    }
                }
            }
        }

        // println!("SS");

        //Now we modify the child pointers of the parents. Here a swap is dangerous as we could have that parent(a) = parent(b)
        if parent_b == parent_a {
            let dummya = TreeNodeId(self.n_nodes());
            let dummyb = TreeNodeId(self.n_nodes() + 1);
            if let ParentId::Node(parent_a) = parent_a {
                for c in &mut self.nodes[parent_a.0].children {
                    if *c == a {
                        *c = dummya;
                    } else if *c == b {
                        *c = dummyb
                    }
                }
                for c in &mut self.nodes[parent_a.0].children {
                    if *c == dummya {
                        *c = b;
                    } else if *c == dummyb {
                        *c = a;
                    }
                }
            }
        } else {
            // it is safe to swap directly
            if let ParentId::Node(parent_b) = parent_b {
                for c in &mut self.nodes[parent_b.0].children {
                    if *c == b {
                        *c = a;
                        break;
                    }
                }
            }

            if let ParentId::Node(parent_a) = parent_a {
                for c in &mut self.nodes[parent_a.0].children {
                    if *c == a {
                        *c = b;
                        break;
                    }
                }
            }
        }

        // println!("SSSS");
        self.nodes.swap(a.0, b.0);
        #[cfg(test)]
        self.panicing_validate();
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
            .map(|(i, node)| (TreeNodeId(i), node.parent_pointer.data.as_ref()))
    }

    fn add_root(&mut self, data: Self::NodeData, root_id: RootId) -> TreeNodeId {
        let node_id = TreeNodeId(self.nodes.len());
        self.nodes.push(CVNode {
            parent_pointer: PPNode::root(data, root_id),
            children: Vec::new(),
        });
        node_id
    }

    fn add_dataless_root(&mut self, root_id: RootId) -> TreeNodeId {
        let node_id = TreeNodeId(self.nodes.len());
        self.nodes.push(CVNode {
            parent_pointer: PPNode::dataless_root(root_id),
            children: Vec::new(),
        });
        node_id
    }

    fn add_child(&mut self, data: Self::NodeData, parent: TreeNodeId) -> TreeNodeId {
        let node_id = TreeNodeId(self.nodes.len());
        self.nodes.push(CVNode {
            parent_pointer: PPNode::child(data, parent),
            children: Vec::new(),
        });
        // Add this child to its parent's list.
        self.nodes[parent.0].children.push(node_id);
        node_id
    }

    fn add_dataless_child(&mut self, parent: TreeNodeId) -> TreeNodeId {
        let node_id = TreeNodeId(self.nodes.len());
        self.nodes.push(CVNode {
            parent_pointer: PPNode::dataless_child(parent),
            children: Vec::new(),
        });
        // Add this child to its parent's list.
        self.nodes[parent.0].children.push(node_id);
        node_id
    }

    fn map<F, U>(self, mut transform: F) -> Self::Store<U>
    where
        F: FnMut(Self::NodeData) -> U,
    {
        let nodes = self
            .nodes
            .into_iter()
            .map(|node| CVNode {
                parent_pointer: node.parent_pointer.map(&mut transform),
                children: node.children, // children remain the same id values.
            })
            .collect();
        ChildVecStore { nodes }
    }

    fn map_ref<F, U>(&self, mut transform: F) -> Self::Store<U>
    where
        F: FnMut(&Self::NodeData) -> U,
    {
        let nodes = self
            .nodes
            .iter()
            .map(|node| CVNode {
                parent_pointer: node.parent_pointer.map_ref(&mut transform),
                children: node.children.clone(),
            })
            .collect();
        ChildVecStore { nodes }
    }
}

impl<V> ForestNodeStoreDown for ChildVecStore<V> {
    fn iter_leaves(&self) -> impl Iterator<Item = TreeNodeId> {
        self.nodes.iter().enumerate().filter_map(|(i, node)| {
            if node.children.is_empty() {
                Some(TreeNodeId(i))
            } else {
                None
            }
        })
    }

    fn iter_children(&self, node_id: TreeNodeId) -> impl Iterator<Item = TreeNodeId> {
        // Clone the children vector (or borrow as iterator) – they are stored in order.
        if node_id.is_empty() {
            Vec::new()
        } else {
            self.nodes[node_id.0].children.clone()
        }
        .into_iter()
    }
}

impl<V> ForestNodeStorePreorder for ChildVecStore<V> {
    type Iterator<'a>
        = PreorderIter<'a, Self>
    where
        Self: 'a;
    fn iter_preorder(&self, start: TreeNodeId) -> Self::Iterator<'_> {
        // Use the generic PreorderIter helper
        PreorderIter::new(self, start)
    }

    fn flatten_branches(&mut self, nodes: impl IntoIterator<Item = TreeNodeId>) {
        // Lift branches directly below a root first. Nodes nested in them are
        // then one hop from the root, so locating the remaining branches does
        // not walk deep hierarchies.
        let (mut direct, mut nested) = (Vec::new(), Vec::new());
        for node in nodes {
            match self.parent(node) {
                ParentId::Node(root) if !matches!(self.parent(root), ParentId::Node(_)) => {
                    direct.push((root, node))
                }
                ParentId::Node(_) => nested.push(node),
                _ => {}
            }
        }
        if direct.is_empty() && nested.is_empty() {
            return;
        }
        let mut selected = SubSet::empty(self.n_nodes());
        self.lift_branches(direct, &mut selected);
        let nested = nested.into_iter().filter_map(|node| self.branch(node));
        self.lift_branches(nested.collect(), &mut selected);
    }
}

impl<V> ChildVecStore<V> {
    /// Move the descendants of each `(root, branch)` directly after the branch
    /// among the root's children, in preorder. `selected` must be empty and is
    /// left empty.
    fn lift_branches(
        &mut self,
        mut branches: Vec<(TreeNodeId, TreeNodeId)>,
        selected: &mut SubSet<TreeNodeId>,
    ) {
        branches.retain(|&(_, branch)| {
            let lift = !self.nodes[branch.0].children.is_empty() && !selected.includes(&branch);
            if lift {
                selected.add(branch);
            }
            lift
        });
        branches.sort_unstable_by_key(|&(root, _)| root);
        let mut stack = Vec::new();
        for group in branches.chunk_by(|a, b| a.0 == b.0) {
            let root = group[0].0;
            let mut children = std::mem::take(&mut self.nodes[root.0].children);
            // Reparenting appends, so branches attached since the last
            // flattening form a suffix of the root's children. Rebuild only
            // that suffix, moving each branch's descendants after it in preorder.
            let mut pending = group.len();
            let mut start = children.len();
            while pending > 0 {
                start -= 1;
                pending -= usize::from(selected.includes(&children[start]));
            }
            for child in children.split_off(start) {
                children.push(child);
                if !selected.includes(&child) {
                    continue;
                }
                selected.sub(child);
                self.append_descendants(root, child, &mut children, &mut stack);
            }
            self.nodes[root.0].children = children;
        }
    }

    /// Push the descendants of `branch` onto `out` in preorder, making each a
    /// direct child of `root`. `stack` must be empty and is left empty.
    fn append_descendants(
        &mut self,
        root: TreeNodeId,
        branch: TreeNodeId,
        out: &mut Vec<TreeNodeId>,
        stack: &mut Vec<(TreeNodeId, usize)>,
    ) {
        // Index the child vectors instead of taking them, so their
        // allocations are kept for reuse rather than freed.
        stack.push((branch, 0));
        while let Some((node, next)) = stack.last_mut() {
            let Some(&descendant) = self.nodes[node.0].children.get(*next) else {
                self.nodes[node.0].children.clear();
                stack.pop();
                continue;
            };
            *next += 1;
            self.nodes[descendant.0].parent_pointer.parent = ParentId::Node(root);
            out.push(descendant);
            stack.push((descendant, 0));
        }
    }
}
impl<V> ForestNodeStoreBfs for ChildVecStore<V> {
    fn iter_bfs(&self, start: TreeNodeId) -> impl Iterator<Item = TreeNodeId> + '_ {
        // Use the generic BfsIter helper
        BfsIter::new(self, start)
    }
}
//
// Conversions
//

// 1. Conversion from ParentPointerStore to ChildVecStore.
//    In a ParentPointerStore, the nodes only have a PPNode (with parent pointer and data).
//    To build the children vector, we iterate over all nodes and add each node as a child
//    of its parent if applicable.
impl<V> From<ParentPointerStore<V>> for ChildVecStore<V> {
    fn from(pp_store: ParentPointerStore<V>) -> Self {
        let n = pp_store.nodes.len();
        let mut some_pp = pp_store.map(|a| Some(a));
        let mut nodes = Vec::with_capacity(n);
        for pp in some_pp.nodes.iter_mut() {
            let data = pp.data.take().flatten();
            let parent_pointer = PPNode {
                data,
                parent: pp.parent,
            };
            nodes.push(CVNode {
                parent_pointer,
                children: vec![],
            });
        }

        for (i, pp) in some_pp.nodes.iter().enumerate() {
            if let ParentId::Node(n) = pp.parent {
                nodes[n.0].children.push(TreeNodeId(i));
            }
        }

        ChildVecStore { nodes }
    }
}

// 2. Conversion from ParentChildStore to ChildVecStore.
//    The ParentChildStore stores child pointers via cyclic neighbor links. We can recover
//    an ordered children vector by using its provided methods.
impl<V> From<ParentChildStore<V>> for ChildVecStore<V> {
    fn from(pc_store: ParentChildStore<V>) -> Self {
        let n = pc_store.nodes.len();
        let mut some_pc = pc_store.map(|a| a);
        let mut nodes = Vec::with_capacity(n);
        for nid in (0..n).map(TreeNodeId) {
            let data = some_pc[nid].take();
            let parent_pointer = PPNode {
                data,
                parent: some_pc.nodes[nid.0].parent_pointer.parent,
            };

            let children = some_pc.iter_children(nid).collect();

            nodes.push(CVNode {
                parent_pointer,
                children,
            });
        }

        ChildVecStore { nodes }
    }
}

// 3. Conversion from ChildVecStore to ParentPointerStore.
//    This conversion simply drops the children vector and keeps only the PPNode.
impl<V> From<ChildVecStore<V>> for ParentPointerStore<V> {
    fn from(store: ChildVecStore<V>) -> Self {
        store
            .nodes
            .into_iter()
            .map(|cv_node| cv_node.parent_pointer)
            .collect()
    }
}

// 4. Conversion from ChildVecStore to ParentChildStore.
//    Here we build a PCNode from each CVNode. For each node, if its children vec is non-empty,
//    we set its 'child' field to the first child, and, for each child, we compute neighbor links
//    (cyclically) using the ordering in the vector.
impl<V> From<ChildVecStore<V>> for ParentChildStore<V> {
    fn from(store: ChildVecStore<V>) -> Self {
        let n = store.nodes.len();
        let mut some_store = store.map(|a| Some(a));
        let mut nodes = Vec::with_capacity(n);

        for nid in (0..n).map(TreeNodeId) {
            let data = some_store.nodes[nid.0].parent_pointer.data.take().flatten();
            let parent_pointer = PPNode {
                data,
                parent: some_store.nodes[nid.0].parent_pointer.parent,
            };

            let child = some_store.nodes[nid.0].children.iter().cloned().next();

            nodes.push(PCNode {
                parent_pointer,
                child,
                neighbor_left: nid,
                neighbor_right: nid,
            });
        }

        for node in some_store.nodes.iter() {
            let len = node.children.len();
            for (j, &child) in node.children.iter().enumerate() {
                let left = node.children[(j + len - 1) % len];
                let right = node.children[(j + 1) % len];
                nodes[child.0].neighbor_left = left;
                nodes[child.0].neighbor_right = right;
            }
        }

        ParentChildStore { nodes }
    }
}

#[cfg(test)]
mod test {
    // use super::*;

    #[test]
    fn swap() {}
}
