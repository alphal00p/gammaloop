#[derive(Clone, Copy, Debug)]
struct NodeIndex(usize);

#[derive(Clone, Copy, Debug)]
struct LeafIndex(usize);

#[derive(Clone)]
struct ForestNode {
    parent: Option<NodeIndex>,
    content: NodeContent,
}

#[derive(Clone)]
enum NodeContent {
    Internal(Vec<NodeIndex>),
    Leaf(LeafIndex),
}

struct Forest<T> {
    nodes: Vec<ForestNode>,
    leaf_data: Vec<T>, // Store leaf data separately for sequential access
}

impl<T> Forest<T> {
    fn new() -> Self {
        Forest {
            nodes: Vec::new(),
            leaf_data: Vec::new(), // Initialize empty leaf data storage
        }
    }

    // Add a new leaf node to the forest
    fn add_leaf(&mut self, data: T) -> LeafIndex {
        let data_index = self.leaf_data.len();
        self.leaf_data.push(data);
        let node_index = self.nodes.len();
        self.nodes.push(ForestNode {
            parent: None,
            content: NodeContent::Leaf(LeafIndex(data_index)),
        });
        LeafIndex(data_index)
    }

    // Sequential access to all leaf data
    pub fn all_leaves(&self) -> &Vec<T> {
        &self.leaf_data
    }

    fn collect_leaf_indices(&self, node_index: NodeIndex) -> Vec<LeafIndex> {
        let mut leaf_indices = Vec::new();
        self.dfs_collect_leaf_indices(node_index, &mut leaf_indices);
        leaf_indices
    }

    // Helper for recursively collecting leaf indices.
    fn dfs_collect_leaf_indices(&self, node_index: NodeIndex, leaf_indices: &mut Vec<LeafIndex>) {
        match &self.nodes[node_index.0].content {
            NodeContent::Leaf(leaf_index) => leaf_indices.push(*leaf_index),
            NodeContent::Internal(children) => {
                for &child_index in children {
                    self.dfs_collect_leaf_indices(child_index, leaf_indices);
                }
            }
        }
    }

    // Access all leaves' data starting from a specific node.
    pub fn get_leaves_from_node(&self, node_index: NodeIndex) -> Vec<&T> {
        let leaf_indices = self.collect_leaf_indices(node_index);
        leaf_indices
            .iter()
            .map(|&leaf_idx| &self.leaf_data[leaf_idx.0])
            .collect()
    }

    // Existing implementation...

    // Union of a vector of leaf indices, creating a new internal node.
    // Returns the NodeIndex of the newly created internal node.
    pub fn union_of_leaves(&mut self, leaf_indices: Vec<LeafIndex>) -> NodeIndex {
        // Convert LeafIndex to NodeIndex for internal representation.
        let child_indices: Vec<NodeIndex> = leaf_indices
            .into_iter()
            .map(|leaf| NodeIndex(leaf.0))
            .collect();

        // Create a new internal node with these children.
        self.union_internal(child_indices)
    }

    // Generalized union function for n-ary unions, works for any node indices.
    pub fn union(&mut self, node_indices: Vec<NodeIndex>) -> NodeIndex {
        self.union_internal(node_indices)
    }

    // Internal function to handle the creation of a new internal node.
    fn union_internal(&mut self, node_indices: Vec<NodeIndex>) -> NodeIndex {
        let new_node = ForestNode {
            parent: None,
            content: NodeContent::Internal(node_indices.clone()),
        };
        let new_node_index = NodeIndex(self.nodes.len());
        self.nodes.push(new_node);

        // Update parent references for each child node.
        for node_index in node_indices {
            if let Some(node) = self.nodes.get_mut(node_index.0) {
                node.parent = Some(new_node_index);
            }
        }

        new_node_index
    }

    // Other methods as previously defined...
}
