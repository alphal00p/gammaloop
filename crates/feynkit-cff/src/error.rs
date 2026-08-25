use crate::{EdgeId, VertexId};
use thiserror::Error;

/// Errors reported while validating a topology or generating a CFF expression.
#[derive(Clone, Debug, Error, PartialEq, Eq)]
pub enum CffError {
    #[error("a CFF graph must contain at least one vertex")]
    EmptyGraph,
    #[error("CFF generation supports at most {maximum} vertices, received {actual}")]
    TooManyVertices { actual: usize, maximum: usize },
    #[error("edge {edge} refers to vertex {vertex}, but the graph has {vertex_count} vertices")]
    UnknownVertex {
        edge: EdgeId,
        vertex: VertexId,
        vertex_count: usize,
    },
    #[error("edge id {0} occurs more than once")]
    DuplicateEdge(EdgeId),
    #[error("the internal topology is disconnected")]
    DisconnectedGraph,
    #[error("edge {0} is not an internal edge")]
    NotInternalEdge(EdgeId),
    #[error("edge {0} does not exist")]
    UnknownEdge(EdgeId),
    #[error("diagram edge {0} connects two external vertices")]
    ExternalToExternalEdge(EdgeId),
    #[error("diagram edge {0} has an endpoint excluded as external without external metadata")]
    MissingExternalMetadata(EdgeId),
    #[error("contracted edge {0} cannot also have a fixed orientation")]
    OrientedContractedEdge(EdgeId),
    #[error("fixed orientation for edge {0} must be default or reversed")]
    UndirectedFixedEdge(EdgeId),
    #[error(
        "generating {candidate_count} orientations exceeds the configured maximum of {maximum}"
    )]
    TooManyOrientations {
        candidate_count: usize,
        maximum: usize,
    },
    #[error("the number of free internal edges cannot be represented on this platform")]
    OrientationCountOverflow,
    #[error("acyclic orientation {orientation} has no source or sink with a connected complement")]
    NoSourceOrSink { orientation: usize },
    #[error("acyclic orientation {orientation} has no valid contraction branch")]
    NoContractionBranch { orientation: usize },
    #[error("CFF branches ended at different depths for orientation {orientation}")]
    UnevenBranchDepth { orientation: usize },
    #[error("the external-shift coefficient for edge {edge} exceeds the i64 range")]
    ExternalShiftOverflow { edge: EdgeId },
    #[error("internal CFF invariant failed: {0}")]
    Invariant(String),
}
