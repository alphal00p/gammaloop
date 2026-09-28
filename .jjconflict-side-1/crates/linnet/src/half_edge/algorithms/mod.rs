//! # Graph Algorithms for Half-Edge Data Structure
//!
//! This module provides various graph algorithms specifically designed to work with
//! the half-edge representation of graphs. All algorithms respect the half-edge
//! invariants and handle edge pairs correctly.
//!
//! ## Available Algorithms
//!
//! ### Topological Operations
//! - [`topological_order`]: Topological sorting using Kahn's algorithm
//!
//! ### Transitive Operations
//! - [`transitive_ops`]: Transitive closure and transitive reduction algorithms
//!
//! ## Half-Edge Considerations
//!
//! All algorithms in this module are designed to work correctly with the half-edge
//! representation where each directed edge consists of a pair of half-edges. When
//! algorithms modify the graph structure (adding or removing edges), they properly
//! handle both half-edges of each pair to maintain graph integrity.

use super::{
    involution::{Flow, Hedge},
    nodestore::NodeStorageOps,
    HedgeGraph,
};

/// Selects which notion of edge direction a directed graph algorithm follows.
///
/// The underlying basis follows the source/sink roles stored by the half-edge
/// involution. The superficial basis also applies each edge's orientation; an
/// undirected superficial orientation therefore has no directed flow.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum DirectionBasis {
    /// Follow the source/sink roles stored by the half-edge involution.
    #[default]
    Underlying,
    /// Follow the displayed orientation, ignoring superficially undirected edges.
    Superficial,
}

impl DirectionBasis {
    pub(super) fn flow<E, V, H, N: NodeStorageOps<NodeData = V>>(
        self,
        graph: &HedgeGraph<E, V, H, N>,
        hedge: Hedge,
    ) -> Option<Flow> {
        match self {
            Self::Underlying => Some(graph.flow(hedge)),
            Self::Superficial => graph.superficial_hedge_orientation(hedge),
        }
    }
}

pub mod topological_order;
pub mod trace_unfold;
pub mod transitive_ops;
