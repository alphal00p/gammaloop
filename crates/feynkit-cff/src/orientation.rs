use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::EdgeId;

/// Orientation relative to the direction stored in the input graph.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum EdgeOrientation {
    #[default]
    Default,
    Reversed,
    /// The edge is absent after contraction and carries no orientation.
    Undirected,
}

/// A complete, edge-id-addressable orientation of the input graph.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct GraphOrientation {
    edges: BTreeMap<EdgeId, EdgeOrientation>,
}

impl GraphOrientation {
    pub fn new(edges: BTreeMap<EdgeId, EdgeOrientation>) -> Self {
        Self { edges }
    }

    pub fn get(&self, edge: EdgeId) -> Option<EdgeOrientation> {
        self.edges.get(&edge).copied()
    }

    pub fn iter(&self) -> impl Iterator<Item = (EdgeId, EdgeOrientation)> + '_ {
        self.edges
            .iter()
            .map(|(edge, orientation)| (*edge, *orientation))
    }

    /// Compare orientations while treating undirected edges in `reference` as wildcards.
    pub fn is_compatible_with(&self, reference: &Self) -> bool {
        reference.iter().all(|(edge, expected)| {
            expected == EdgeOrientation::Undirected || self.get(edge) == Some(expected)
        })
    }
}

/// Stable index of an orientation within a generated expression.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct OrientationId(pub usize);

impl OrientationId {
    pub const fn index(self) -> usize {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn undirected_reference_edges_are_wildcards() {
        let edge = EdgeId::new(3);
        let orientation =
            GraphOrientation::new(BTreeMap::from([(edge, EdgeOrientation::Reversed)]));
        let reference =
            GraphOrientation::new(BTreeMap::from([(edge, EdgeOrientation::Undirected)]));

        assert!(orientation.is_compatible_with(&reference));
        assert!(!reference.is_compatible_with(&orientation));
    }
}
