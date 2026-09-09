use thiserror::Error;

#[derive(Debug, Error)]
pub enum OneLoopError {
    /// The input graph is not one loop.
    #[error("unsupported loop order: graph has {found} loops, only one-loop is supported")]
    UnsupportedLoopOrder { found: usize },

    /// A `gammalooprs::Graph` could not be turned into an `IntegralFamily`.
    #[error("failed to extract integral family from graph: {reason}")]
    ExtractionFailed { reason: String },

    /// Wraps an underlying Symbolica error surfaced during extraction.
    #[error("symbolica error: {0}")]
    Symbolica(String),
}
