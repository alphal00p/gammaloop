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

pub mod topological_order;
pub mod trace_unfold;
pub mod transitive_ops;
