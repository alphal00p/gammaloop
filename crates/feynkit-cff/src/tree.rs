use std::collections::{BTreeMap, BTreeSet};

use serde::{Deserialize, Deserializer, Serialize, de};
use symbolica::{
    atom::{Atom, AtomCore},
    parse,
};
use thiserror::Error;

/// Stable index of a node within an expression tree.
#[derive(
    Clone,
    Copy,
    Debug,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    Serialize,
    Deserialize,
    bincode::Encode,
    bincode::Decode,
)]
pub struct NodeId(pub usize);

impl NodeId {
    pub const ROOT: Self = Self(0);

    pub const fn root() -> Self {
        Self::ROOT
    }

    pub const fn index(self) -> usize {
        self.0
    }
}

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

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, bincode::Encode, bincode::Decode)]
pub struct TreeNode<T> {
    pub data: T,
    pub id: NodeId,
    pub children: Vec<NodeId>,
    pub parent: Option<NodeId>,
}

/// Rooted tree encoding a sum of products of inverse surfaces.
///
/// Every root-to-leaf path is one product of inverse surface factors, while
/// sibling branches are summed.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, bincode::Encode, bincode::Decode)]
pub struct ExpressionTree<T> {
    nodes: Vec<TreeNode<T>>,
}

#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum ExpressionTreeError {
    #[error("an expression tree must contain a root node")]
    Empty,
    #[error("tree node at position {position} carries id {id:?}")]
    NodeIdMismatch { position: usize, id: NodeId },
    #[error("the root node must not have a parent")]
    RootHasParent,
    #[error("non-root tree node {0:?} has no parent")]
    MissingParent(NodeId),
    #[error("tree node {node:?} refers to unknown parent {parent:?}")]
    UnknownParent { node: NodeId, parent: NodeId },
    #[error("tree node {node:?} refers to unknown child {child:?}")]
    UnknownChild { node: NodeId, child: NodeId },
    #[error("tree node {child:?} occurs under more than one parent")]
    DuplicateChild { child: NodeId },
    #[error("tree node {child:?} does not point back to parent {parent:?}")]
    ParentMismatch { parent: NodeId, child: NodeId },
    #[error("tree node {0:?} is not reachable from the root")]
    Unreachable(NodeId),
}

impl<T> ExpressionTree<T> {
    pub fn from_root(data: T) -> Self {
        Self {
            nodes: vec![TreeNode {
                data,
                id: NodeId::ROOT,
                children: Vec::new(),
                parent: None,
            }],
        }
    }

    pub fn nodes(&self) -> &[TreeNode<T>] {
        &self.nodes
    }

    pub fn node(&self, id: NodeId) -> Option<&TreeNode<T>> {
        self.nodes.get(id.index())
    }

    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }

    pub fn validate(&self) -> Result<(), ExpressionTreeError> {
        if self.nodes.is_empty() {
            return Err(ExpressionTreeError::Empty);
        }
        for (position, node) in self.nodes.iter().enumerate() {
            if node.id.index() != position {
                return Err(ExpressionTreeError::NodeIdMismatch {
                    position,
                    id: node.id,
                });
            }
        }
        if self.nodes[NodeId::ROOT.index()].parent.is_some() {
            return Err(ExpressionTreeError::RootHasParent);
        }

        let mut owners = vec![None; self.nodes.len()];
        for node in &self.nodes {
            for child in &node.children {
                let Some(child_node) = self.nodes.get(child.index()) else {
                    return Err(ExpressionTreeError::UnknownChild {
                        node: node.id,
                        child: *child,
                    });
                };
                if owners[child.index()].replace(node.id).is_some() {
                    return Err(ExpressionTreeError::DuplicateChild { child: *child });
                }
                if child_node.parent != Some(node.id) {
                    return Err(ExpressionTreeError::ParentMismatch {
                        parent: node.id,
                        child: *child,
                    });
                }
            }
        }
        for node in self.nodes.iter().skip(1) {
            let parent = node
                .parent
                .ok_or(ExpressionTreeError::MissingParent(node.id))?;
            if parent.index() >= self.nodes.len() {
                return Err(ExpressionTreeError::UnknownParent {
                    node: node.id,
                    parent,
                });
            }
            if owners[node.id.index()] != Some(parent) {
                return Err(ExpressionTreeError::ParentMismatch {
                    parent,
                    child: node.id,
                });
            }
        }

        let mut seen = vec![false; self.nodes.len()];
        let mut pending = vec![NodeId::ROOT];
        while let Some(node) = pending.pop() {
            if seen[node.index()] {
                return Err(ExpressionTreeError::DuplicateChild { child: node });
            }
            seen[node.index()] = true;
            pending.extend(self.nodes[node.index()].children.iter().copied());
        }
        if let Some(position) = seen.iter().position(|seen| !seen) {
            return Err(ExpressionTreeError::Unreachable(NodeId(position)));
        }
        Ok(())
    }

    pub fn leaves(&self) -> impl Iterator<Item = &TreeNode<T>> {
        self.nodes.iter().filter(|node| node.children.is_empty())
    }

    pub fn leaf_ids(&self) -> Vec<NodeId> {
        self.leaves().map(|node| node.id).collect()
    }

    pub fn insert(&mut self, parent: NodeId, data: T) -> Option<NodeId> {
        if parent.index() >= self.nodes.len() {
            return None;
        }
        let id = NodeId(self.nodes.len());
        self.nodes.push(TreeNode {
            data,
            id,
            children: Vec::new(),
            parent: Some(parent),
        });
        self.nodes[parent.index()].children.push(id);
        Some(id)
    }

    pub fn insert_node(&mut self, parent: NodeId, data: T) {
        self.insert(parent, data)
            .expect("cannot insert below a missing expression-tree node");
    }

    pub fn get_node(&self, id: NodeId) -> &TreeNode<T> {
        self.node(id)
            .expect("expression tree does not contain the requested node")
    }

    pub fn get_bottom_layer(&self) -> Vec<NodeId> {
        self.leaf_ids()
    }

    pub fn get_num_nodes(&self) -> usize {
        self.node_count()
    }

    pub fn iter_nodes(&self) -> impl Iterator<Item = &TreeNode<T>> {
        self.nodes.iter()
    }

    pub fn apply_mut_closure(&mut self, id: NodeId, f: impl FnOnce(&mut T)) {
        let data = self
            .data_mut(id)
            .expect("expression tree does not contain the requested node");
        f(data);
    }

    pub(crate) fn data_mut(&mut self, id: NodeId) -> Option<&mut T> {
        self.nodes.get_mut(id.index()).map(|node| &mut node.data)
    }

    pub(crate) fn try_map<U, E>(
        self,
        mut f: impl FnMut(T) -> Result<U, E>,
    ) -> Result<ExpressionTree<U>, E> {
        Ok(ExpressionTree {
            nodes: self
                .nodes
                .into_iter()
                .map(|node| {
                    Ok(TreeNode {
                        data: f(node.data)?,
                        id: node.id,
                        children: node.children,
                        parent: node.parent,
                    })
                })
                .collect::<Result<_, E>>()?,
        })
    }

    pub fn map<U>(self, mut f: impl FnMut(T) -> U) -> ExpressionTree<U> {
        ExpressionTree {
            nodes: self
                .nodes
                .into_iter()
                .map(|node| TreeNode {
                    data: f(node.data),
                    id: node.id,
                    children: node.children,
                    parent: node.parent,
                })
                .collect(),
        }
    }

    pub fn map_mut(&mut self, mut f: impl FnMut(&mut T)) {
        for node in &mut self.nodes {
            f(&mut node.data);
        }
    }

    pub fn root_to_leaf_paths(&self) -> Vec<Vec<&T>> {
        self.leaves()
            .map(|leaf| {
                let mut path = Vec::new();
                let mut current = Some(leaf.id);
                while let Some(id) = current {
                    let node = &self.nodes[id.index()];
                    path.push(&node.data);
                    current = node.parent;
                }
                path.reverse();
                path
            })
            .collect()
    }

    /// Keep exactly the root-to-leaf branches containing `count` occurrences of `value`.
    ///
    /// Shared prefixes are retained once. If no branch matches, the tree becomes empty.
    pub fn retain_branches_with_value_count(&mut self, value: &T, count: usize)
    where
        T: Eq,
    {
        if self.nodes.is_empty() {
            return;
        }

        let mut retained = BTreeSet::new();
        for leaf in self.leaf_ids() {
            let path = self.path_to_root(leaf);
            if path
                .iter()
                .filter(|id| self.nodes[id.index()].data == *value)
                .count()
                == count
            {
                retained.extend(path);
            }
        }

        if retained.is_empty() {
            self.nodes.clear();
            return;
        }

        let remap = retained
            .iter()
            .enumerate()
            .map(|(new, old)| (*old, NodeId(new)))
            .collect::<BTreeMap<_, _>>();
        let old_nodes = std::mem::take(&mut self.nodes);
        self.nodes = old_nodes
            .into_iter()
            .filter(|node| retained.contains(&node.id))
            .map(|node| TreeNode {
                data: node.data,
                id: remap[&node.id],
                children: node
                    .children
                    .into_iter()
                    .filter_map(|child| remap.get(&child).copied())
                    .collect(),
                parent: node.parent.map(|parent| remap[&parent]),
            })
            .collect();
    }

    pub fn keep_branches_with_value_count_mut(&mut self, value: &T, count: usize)
    where
        T: Eq,
    {
        self.retain_branches_with_value_count(value, count);
    }

    /// Largest number of occurrences of `value` on any root-to-leaf branch.
    pub fn max_value_count_on_branch(&self, value: &T) -> usize
    where
        T: Eq,
    {
        self.leaf_ids()
            .into_iter()
            .map(|leaf| {
                self.path_to_root(leaf)
                    .into_iter()
                    .filter(|id| self.nodes[id.index()].data == *value)
                    .count()
            })
            .max()
            .unwrap_or(0)
    }

    fn path_to_root(&self, mut node: NodeId) -> Vec<NodeId> {
        let mut path = Vec::new();
        loop {
            path.push(node);
            match self.nodes[node.index()].parent {
                Some(parent) => node = parent,
                None => break,
            }
        }
        path.reverse();
        path
    }
}

impl<T> ExpressionTree<T>
where
    Atom: From<T>,
    T: Copy,
{
    /// Lower the tree to a sum of products of inverse Symbolica surfaces.
    pub fn to_atom_inverse(&self) -> Atom {
        fn lower<T>(tree: &ExpressionTree<T>, current: NodeId) -> Atom
        where
            Atom: From<T>,
            T: Copy,
        {
            let node = tree.get_node(current);
            let inverse_surface = (Atom::num(1) / Atom::from(node.data))
                .replace(parse!("η_inf^-1"))
                .with(Atom::num(0));
            let children = node
                .children
                .iter()
                .map(|child| lower(tree, *child))
                .reduce(|sum, term| sum + term)
                .unwrap_or_else(|| Atom::num(1));
            inverse_surface * children
        }

        if self.is_empty() {
            Atom::Zero
        } else {
            lower(self, NodeId::ROOT)
        }
    }
}

#[derive(Deserialize)]
struct ExpressionTreeData<T> {
    nodes: Vec<TreeNode<T>>,
}

impl<'de, T> Deserialize<'de> for ExpressionTree<T>
where
    T: Deserialize<'de>,
{
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let data = ExpressionTreeData::deserialize(deserializer)?;
        let tree = Self { nodes: data.nodes };
        tree.validate().map_err(de::Error::custom)?;
        Ok(tree)
    }
}

#[cfg(test)]
mod tests {
    use serde_json::json;

    use super::*;

    #[test]
    fn deserialization_revalidates_tree_links() {
        let invalid = json!({
            "nodes": [{
                "data": 1,
                "id": 0,
                "children": [1],
                "parent": null
            }]
        });
        assert!(serde_json::from_value::<ExpressionTree<u8>>(invalid).is_err());
    }

    #[test]
    fn retaining_branches_preserves_shared_prefixes_and_reindexes_nodes() {
        let mut tree = ExpressionTree::from_root(9);
        let prefix = tree.insert(NodeId::ROOT, 8).unwrap();
        let matching = tree.insert(prefix, 1).unwrap();
        tree.insert(matching, 5).unwrap();
        let rejected = tree.insert(prefix, 2).unwrap();
        tree.insert(rejected, 6).unwrap();

        tree.retain_branches_with_value_count(&1, 1);

        assert_eq!(tree.root_to_leaf_paths(), vec![vec![&9, &8, &1, &5]]);
        tree.validate().unwrap();
    }

    #[test]
    fn retaining_an_absent_count_can_empty_a_tree() {
        let mut tree = ExpressionTree::from_root(1);
        tree.insert(NodeId::ROOT, 2).unwrap();

        tree.retain_branches_with_value_count(&1, 2);

        assert!(tree.is_empty());
    }

    #[test]
    fn maximum_value_count_is_computed_per_branch() {
        let mut tree = ExpressionTree::from_root(1);
        let repeated = tree.insert(NodeId::ROOT, 1).unwrap();
        tree.insert(repeated, 1).unwrap();
        tree.insert(NodeId::ROOT, 2).unwrap();

        assert_eq!(tree.max_value_count_on_branch(&1), 3);
    }

    #[test]
    fn symbolic_lowering_sums_sibling_branches() {
        let mut tree = ExpressionTree::from_root(1_i64);
        tree.insert(NodeId::ROOT, 2).unwrap();
        tree.insert(NodeId::ROOT, 3).unwrap();

        assert_eq!(tree.to_atom_inverse(), Atom::num(5) / Atom::num(6));
    }
}
