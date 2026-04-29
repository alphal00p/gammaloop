use std::collections::HashSet;

use bincode_trait_derive::{Decode, Encode};
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomCore},
    parse,
};
use typed_index_collections::TiVec;

#[derive(
    Debug, Copy, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, Eq, PartialOrd, Ord, Hash,
)]
pub struct NodeId(pub usize);

impl From<usize> for NodeId {
    fn from(value: usize) -> Self {
        Self(value)
    }
}

impl From<NodeId> for usize {
    fn from(value: NodeId) -> Self {
        value.0
    }
}

impl NodeId {
    pub const fn root() -> Self {
        Self(0)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, Eq, PartialEq)]
pub struct TreeNode<T> {
    pub data: T,
    pub node_id: NodeId,
    pub children: Vec<NodeId>,
    pub parent: Option<NodeId>,
}

fn determine_shifted_id(removed_ids_sorted: &[NodeId], original_id: NodeId) -> NodeId {
    let shift = removed_ids_sorted
        .iter()
        .filter(|&&removed_id| removed_id < original_id)
        .count();
    NodeId(original_id.0 - shift)
}

impl<T> TreeNode<T> {
    fn update_node_ids(&mut self, removed_ids: &[NodeId]) {
        self.node_id = determine_shifted_id(removed_ids, self.node_id);

        self.children
            .retain(|child_id| !removed_ids.contains(child_id));
        self.children.iter_mut().for_each(|child_id| {
            *child_id = determine_shifted_id(removed_ids, *child_id);
        });

        if let Some(parent_id) = self.parent {
            if removed_ids.contains(&parent_id) {
                self.parent = None;
            } else {
                self.parent = Some(determine_shifted_id(removed_ids, parent_id));
            }
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, Eq, PartialEq)]
pub struct Tree<T> {
    nodes: TiVec<NodeId, TreeNode<T>>,
}

impl<T> Tree<T> {
    pub fn from_root(data: T) -> Self {
        let node_id = NodeId(0);
        let root_node = TreeNode {
            data,
            node_id,
            children: Vec::new(),
            parent: None,
        };
        Self {
            nodes: vec![root_node].into(),
        }
    }

    pub fn insert_node(&mut self, parent_id: NodeId, data: T) {
        let node_id = NodeId(self.nodes.len());
        let new_node = TreeNode {
            data,
            node_id,
            children: Vec::new(),
            parent: Some(parent_id),
        };
        self.nodes.push(new_node);
        self.nodes[parent_id].children.push(node_id);
    }

    pub fn get_node(&self, node_id: NodeId) -> &TreeNode<T> {
        &self.nodes[node_id]
    }

    pub fn get_data_node_mut(&mut self, node_id: NodeId) -> &mut T {
        &mut self.nodes[node_id].data
    }

    pub fn apply_mut_closure(&mut self, node_id: NodeId, closure: impl Fn(&mut T)) {
        closure(&mut self.nodes[node_id].data);
    }

    pub fn get_bottom_layer(&self) -> Vec<NodeId> {
        self.nodes
            .iter()
            .filter(|node| node.children.is_empty())
            .map(|node| node.node_id)
            .collect()
    }

    pub fn map<G>(self, f: impl Fn(T) -> G) -> Tree<G> {
        Tree {
            nodes: self
                .nodes
                .into_iter()
                .map(|node| TreeNode {
                    data: f(node.data),
                    node_id: node.node_id,
                    children: node.children,
                    parent: node.parent,
                })
                .collect(),
        }
    }

    pub fn map_mut(&mut self, f: impl Fn(&mut T)) {
        for node in &mut self.nodes {
            f(&mut node.data);
        }
    }

    pub fn get_num_nodes(&self) -> usize {
        self.nodes.len()
    }

    pub fn iter_nodes(&self) -> impl Iterator<Item = &TreeNode<T>> {
        self.nodes.iter()
    }

    fn path_to_root(&self, mut node_id: NodeId) -> Vec<NodeId> {
        let mut path = Vec::new();
        loop {
            path.push(node_id);
            if let Some(parent_id) = self.nodes[node_id].parent {
                node_id = parent_id;
            } else {
                break;
            }
        }
        path.reverse();
        path
    }

    fn obtain_subtree_node_ids_impl(&self, node_id: NodeId, collected_ids: &mut Vec<NodeId>) {
        for &child_id in &self.nodes[node_id].children {
            collected_ids.push(child_id);
            self.obtain_subtree_node_ids_impl(child_id, collected_ids);
        }
    }

    fn obtain_subtree_node_ids(&self, node_id: NodeId) -> Vec<NodeId> {
        let mut collected_ids = vec![node_id];
        self.obtain_subtree_node_ids_impl(node_id, &mut collected_ids);
        collected_ids
    }

    pub fn remove_node(&mut self, node_id: NodeId) {
        let mut subtree = self.obtain_subtree_node_ids(node_id);
        subtree.sort();

        for id in subtree.iter().rev() {
            self.nodes.remove(*id);
        }

        for node in &mut self.nodes {
            node.update_node_ids(&subtree);
        }
    }

    pub fn filter_mut(&mut self, predicate: impl Fn(&T) -> bool) {
        let mut nodes_to_remove = HashSet::new();

        for node in &self.nodes {
            if !predicate(&node.data) {
                nodes_to_remove.extend(self.obtain_subtree_node_ids(node.node_id));
            }
        }

        let nodes_to_remove: Vec<NodeId> = nodes_to_remove.into_iter().sorted().collect();

        for id in nodes_to_remove.iter().rev() {
            self.nodes.remove(*id);
        }

        for node in &mut self.nodes {
            node.update_node_ids(&nodes_to_remove);
        }
    }

    pub fn keep_branches_with_value_count_mut(&mut self, value: &T, n: usize)
    where
        T: Eq,
    {
        if self.nodes.is_empty() {
            return;
        }

        let leaves = self.get_bottom_layer();
        let mut nodes_to_keep = HashSet::new();
        let mut has_match = false;

        for leaf in leaves {
            let path = self.path_to_root(leaf);
            let count = path
                .iter()
                .filter(|&&node_id| self.nodes[node_id].data == *value)
                .count();

            if count == n {
                has_match = true;
                nodes_to_keep.extend(path);
            }
        }

        if !has_match {
            self.nodes.clear();
            return;
        }

        let nodes_to_remove: Vec<NodeId> = self
            .nodes
            .iter()
            .filter(|node| !nodes_to_keep.contains(&node.node_id))
            .map(|node| node.node_id)
            .sorted()
            .collect();

        for id in nodes_to_remove.iter().rev() {
            self.nodes.remove(*id);
        }

        for node in &mut self.nodes {
            node.update_node_ids(&nodes_to_remove);
        }
    }

    pub fn max_value_count_on_branch(&self, value: &T) -> usize
    where
        T: Eq,
    {
        if self.nodes.is_empty() {
            return 0;
        }

        self.get_bottom_layer()
            .into_iter()
            .map(|leaf| {
                self.path_to_root(leaf)
                    .into_iter()
                    .filter(|&node_id| self.nodes[node_id].data == *value)
                    .count()
            })
            .max()
            .unwrap_or(0)
    }
}

impl<T> Tree<T>
where
    Atom: From<T>,
    T: Copy,
{
    fn to_atom_inv_impl(&self, cur_node: NodeId) -> Atom {
        let node = &self.nodes[cur_node];
        let inv_data_esurface = (Atom::num(1) / Atom::from(node.data))
            .replace(parse!("η_inf^-1"))
            .with(Atom::num(0));

        let child_sum = node
            .children
            .iter()
            .map(|&child| self.to_atom_inv_impl(child))
            .reduce(|acc, x| acc + x)
            .unwrap_or(Atom::num(1));

        inv_data_esurface * child_sum
    }

    pub fn to_atom_inv(&self) -> Atom {
        if self.nodes.is_empty() {
            return Atom::num(0);
        }
        self.to_atom_inv_impl(NodeId::root())
    }
}

#[cfg(test)]
mod tests {
    use super::{NodeId, Tree};

    #[test]
    fn test_remove_node() {
        let mut tree = Tree::from_root(0);
        tree.insert_node(NodeId(0), 1);
        tree.insert_node(NodeId(0), 2);
        tree.insert_node(NodeId(2), 3);
        tree.insert_node(NodeId(2), 4);
        tree.insert_node(NodeId(4), 5);
        tree.insert_node(NodeId(4), 6);
        tree.insert_node(NodeId(4), 7);
        tree.insert_node(NodeId(7), 8);
        tree.insert_node(NodeId(7), 9);

        tree.remove_node(NodeId(4));

        let mut expected_tree = Tree::from_root(0);
        expected_tree.insert_node(NodeId(0), 1);
        expected_tree.insert_node(NodeId(0), 2);
        expected_tree.insert_node(NodeId(2), 3);

        assert_eq!(tree, expected_tree);
    }

    #[test]
    fn test_keep_branches_with_value_count_mut() {
        let mut tree = Tree::from_root(0);
        tree.insert_node(NodeId(0), 1);
        tree.insert_node(NodeId(1), 2);
        tree.insert_node(NodeId(2), 1);
        tree.insert_node(NodeId(0), 1);
        tree.insert_node(NodeId(4), 3);
        tree.insert_node(NodeId(5), 4);

        tree.keep_branches_with_value_count_mut(&1, 1);

        let mut expected_tree = Tree::from_root(0);
        expected_tree.insert_node(NodeId(0), 1);
        expected_tree.insert_node(NodeId(1), 3);
        expected_tree.insert_node(NodeId(2), 4);

        assert_eq!(tree, expected_tree);
    }
}
