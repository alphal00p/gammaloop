//! Helper iterator implementations used by different ForestNodeStore types.

use super::{ForestNodeStore, ForestNodeStoreDown, ParentId, TreeNodeId};
use std::collections::VecDeque;

// --- Ancestors Iterator ---

/// An iterator that traverses upwards from a starting node to its root.
/// Used by the default implementation of `ForestNodeStoreAncestors`.
#[derive(Clone)]
pub struct AncestorsIter<'a, S: ForestNodeStore> {
    store: &'a S,
    /// The next node ID to yield. `None` when the root has been yielded.
    current: Option<TreeNodeId>,
}

impl<'a, S: ForestNodeStore> AncestorsIter<'a, S> {
    /// Creates a new ancestor iterator.
    pub fn new(store: &'a S, start_node: TreeNodeId) -> Self {
        AncestorsIter {
            store,
            current: Some(start_node),
        }
    }
}

impl<S: ForestNodeStore> Iterator for AncestorsIter<'_, S> {
    type Item = TreeNodeId;

    fn next(&mut self) -> Option<Self::Item> {
        let node_to_return = self.current?;
        match self.store[&node_to_return] {
            ParentId::Root(_) => self.current = None,
            ParentId::Node(parent_id) => self.current = Some(parent_id),
            ParentId::PointingRoot(parent_id) => self.current = Some(parent_id),
        }
        Some(node_to_return)
    }
}

// --- BFS Iterator Helper (Example for ChildVecStore / ParentChildStore) ---
// These could live here or within their respective store modules.
// Kept generic here for demonstration.

/// A Breadth-First Search (BFS) iterator state.
#[derive(Clone)]
pub struct BfsIter<'a, S: ForestNodeStore>
where
    S: ForestNodeStoreDown, // BFS requires iterating children
{
    store: &'a S,
    queue: VecDeque<TreeNodeId>,
}

impl<'a, S: ForestNodeStoreDown> BfsIter<'a, S> {
    /// Create a new BFS iterator starting at `start`.
    pub fn new(store: &'a S, start: TreeNodeId) -> Self {
        let mut queue = VecDeque::new();
        queue.push_back(start);
        BfsIter { store, queue }
    }
}

impl<S: ForestNodeStoreDown> Iterator for BfsIter<'_, S> {
    type Item = TreeNodeId;
    fn next(&mut self) -> Option<Self::Item> {
        let node = self.queue.pop_front()?;
        for child in self.store.iter_children(node) {
            self.queue.push_back(child);
        }
        Some(node)
    }
}

// --- Preorder Iterator Helper (Example for ChildVecStore / ParentChildStore) ---

/// A pre–order DFS iterator state.
// #[derive(Clone)]
pub struct PreorderIter<'a, S: ForestNodeStore>
where
    S: ForestNodeStoreDown, // Preorder needs children access
{
    store: &'a S,
    /// Stack for DFS traversal. Stores nodes to visit.
    stack: Vec<TreeNodeId>,
    // Alternatively, could use the ParentChildStore's original logic without explicit stack
    // if implemented specifically for that store.
}

impl<S: ForestNodeStoreDown> Clone for PreorderIter<'_, S> {
    fn clone(&self) -> Self {
        Self {
            store: self.store,
            stack: self.stack.clone(),
        }
    }
}

impl<'a, S: ForestNodeStoreDown> PreorderIter<'a, S> {
    /// Create a new pre-order iterator starting at `start`.
    pub fn new(store: &'a S, start: TreeNodeId) -> Self {
        PreorderIter {
            store,
            stack: vec![start],
        }
    }
}

impl<S: ForestNodeStoreDown> Iterator for PreorderIter<'_, S> {
    type Item = TreeNodeId;
    fn next(&mut self) -> Option<Self::Item> {
        // Pop the next node from the stack
        let node = self.stack.pop()?;

        // Push children onto the stack in reverse order so the first child is processed next
        // (This assumes iter_children returns them in the desired forward order)
        for child in self
            .store
            .iter_children(node)
            .collect::<Vec<_>>()
            .into_iter()
            .rev()
        {
            self.stack.push(child);
        }

        Some(node)
    }
}

// --- Neighbor Iterator (Specific to ParentChildStore) ---
// Belongs in child_pointer.rs
