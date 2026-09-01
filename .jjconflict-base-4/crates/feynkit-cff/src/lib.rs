//! Standalone Cross-Free Family (CFF) combinatorics.
//!
//! The crate owns the canonical CFF graph, generation algorithms, surface
//! arena, orientation and denominator trees, residue construction, and
//! Symbolica lowering. It deliberately has no knowledge of GammaLoop settings,
//! numerical evaluators, or runtime threshold classification. Callers may
//! construct a [`CffGraph`] explicitly or invoke the extension trait directly
//! on a finalized FeynKit diagram.
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
pub use expression::{CffExpression, OrientationData, OrientationExpression};
pub use generation::{
    CffGeneration, CffGenerator, CffOptions, CffReport, CffResult, FeynmanDiagramCffExt,
    HedgeEdgeRole, HedgeGraphCffExt, ShiftRewrite,
};
pub use graph::{CffEdge, CffGraph, EdgeFlow, EdgeKind};
pub use ids::{EdgeId, VertexId};
pub use orientation::{EdgeOrientation, GraphOrientation, GraphOrientationExt, OrientationId};
pub use surface::{
    EnergySurface, EnergySurfaceId, ExternalShift, HSurface, HSurfaceId, RaisedEnergySurfaceData,
    RaisedEnergySurfaceGroup, RaisedEnergySurfaceId, Surface, SurfaceCache, SurfaceId,
    SurfaceIdMap, VertexSet,
};
pub use tree::{ExpressionTree, ExpressionTreeError, NodeId, TreeNode};
