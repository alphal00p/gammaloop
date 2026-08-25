//! Standalone Cross-Free Family (CFF) combinatorics.
//!
//! The crate owns only topology-level CFF generation. It deliberately has no
//! knowledge of GammaLoop settings, numerical evaluators, UV subtraction, or
//! runtime threshold classification. Callers provide a [`CffGraph`] and receive
//! typed orientation and denominator trees that can be translated to their own
//! symbolic representation.
//!
//! ```
//! use feynkit_cff::{
//!     CffEdge, CffGenerator, CffGraph, EdgeFlow, EdgeId, VertexId,
//! };
//!
//! let graph = CffGraph::new(
//!     2,
//!     [
//!         CffEdge::external(EdgeId::new(0), VertexId::new(0), EdgeFlow::Incoming),
//!         CffEdge::external(EdgeId::new(1), VertexId::new(1), EdgeFlow::Outgoing),
//!         CffEdge::internal(EdgeId::new(2), VertexId::new(0), VertexId::new(1)),
//!     ],
//! )?;
//! let result = CffGenerator::default().generate(&graph)?;
//! assert_eq!(result.report.acyclic_orientations, 2);
//! # Ok::<_, feynkit_cff::CffError>(())
//! ```

#![forbid(unsafe_code)]

mod error;
mod expression;
mod generation;
mod graph;
mod ids;
mod orientation;
mod surface;
mod tree;

pub use error::CffError;
pub use expression::{CffExpression, OrientationExpression};
pub use generation::{CffGenerator, CffOptions, CffReport, CffResult, ShiftRewrite};
pub use graph::{CffEdge, CffGraph, EdgeFlow, EdgeKind};
pub use ids::{EdgeId, VertexId};
pub use orientation::{EdgeOrientation, GraphOrientation, OrientationId};
pub use surface::{
    EnergySurface, EnergySurfaceId, ExternalShift, HSurface, HSurfaceId, Surface, SurfaceCache,
    SurfaceId, VertexSet,
};
pub use tree::{ExpressionTree, ExpressionTreeError, NodeId, TreeNode};
