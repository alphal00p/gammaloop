use serde::{Deserialize, Deserializer, Serialize, de};
use thiserror::Error;

/// Stable index of a node within an expression tree.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct NodeId(pub usize);

impl NodeId {
    pub const ROOT: Self = Self(0);

    pub const fn index(self) -> usize {
        self.0
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
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
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
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
    pub(crate) fn from_root(data: T) -> Self {
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

    pub(crate) fn leaf_ids(&self) -> Vec<NodeId> {
        self.leaves().map(|node| node.id).collect()
    }

    pub(crate) fn insert(&mut self, parent: NodeId, data: T) -> Option<NodeId> {
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
}
