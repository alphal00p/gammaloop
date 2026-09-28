//! Standalone Cross-Free Family (CFF) combinatorics.
//!
//! The crate exposes topology-only and generalized CFF generation. Both use
//! the recursion owned by `three-dimensional-reps`; FeynKit provides graph
//! adapters, surface views, residue construction, and Symbolica lowering. It deliberately has no knowledge of GammaLoop settings,
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

mod algebra;
mod error;
mod expression;
mod generation;
mod graph;
mod ids;
mod orientation;
mod surface;
pub mod symbols;
mod tree;

pub use algebra::{CutPropagator, SurfacePole};
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

/// Generalized CFF generation, including repeated poles and bounded numerator
/// energy dependence, shared with GammaLoop's production runtime.
pub use three_dimensional_reps as generalized;
pub use three_dimensional_reps::{
    Generate3DExpressionOptions, GeneratedThreeDExpression, ParsedGraph, ThreeDGraphSource,
    generate_3d_expression,
};
