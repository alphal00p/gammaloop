use std::ops::{Index, IndexMut};

use serde::{Deserialize, Serialize};

/// data structure for a tree

#[derive(Debug, Copy, Clone, Serialize, Deserialize)]
pub struct NodeId(usize);

impl NodeId {
    pub const fn root() -> Self {
        NodeId(0)
    }
}

impl From<usize> for NodeId {
    fn from(value: usize) -> Self {
        NodeId(value)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TreeNode<T> {
    pub data: T,
    pub node_id: NodeId,
    pub children: Vec<NodeId>,
    pub parent: Option<NodeId>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Tree<T> {
    nodes: Vec<TreeNode<T>>,
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
        Tree {
            nodes: vec![root_node],
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
        self.nodes[parent_id.0].children.push(node_id);
    }

    pub fn get_node(&self, node_id: NodeId) -> &TreeNode<T> {
        &self.nodes[node_id.0]
    }

    pub fn get_data_node_mut(&mut self, node_id: NodeId) -> &mut T {
        &mut self.nodes[node_id.0].data
    }

    pub fn apply_mut_closure(&mut self, node_id: NodeId, closure: impl Fn(&mut T)) {
        closure(&mut self.nodes[node_id.0].data);
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

    pub fn get_num_nodes(&self) -> usize {
        self.nodes.len()
    }
}

pub struct NodeCache<T> {
    cache: Vec<T>,
}

impl<T: Clone> NodeCache<T> {
    pub fn new_default(num_nodes: usize, default: T) -> Self {
        NodeCache {
            cache: vec![default; num_nodes],
        }
    }

    pub fn insert(&mut self, node_id: NodeId, data: T) {
        self.cache[node_id.0] = data;
    }
}

impl<T> Index<NodeId> for NodeCache<T> {
    type Output = T;

    fn index(&self, node_id: NodeId) -> &T {
        &self.cache[node_id.0]
    }
}

impl<T> IndexMut<NodeId> for NodeCache<T> {
    fn index_mut(&mut self, node_id: NodeId) -> &mut T {
        &mut self.cache[node_id.0]
    }
}
