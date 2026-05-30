//! The crate's single error type. Every fallible `oneloop` operation returns
//! `Result<_, OneLoopError>`

use thiserror::Error;

#[derive(Debug, Error)]
pub enum OneLoopError {
    /// The input graph is not one loop.
    #[error("unsupported loop order: graph has {found} loops, only one-loop is supported")]
    UnsupportedLoopOrder { found: usize },

    /// A `gammalooprs::Graph` could not be turned into an `IntegralFamily`
    #[error("failed to extract integral family from graph: {reason}")]
    ExtractionFailed { reason: String },

    /// The IBP system did not close onto the registered masters
    #[error("IBP reduction did not close onto masters: {0}")]
    SolverFailed(String),

    /// A required master configuration has no analytic closed form yet
    #[error("no analytic closed form registered for master: {which}")]
    MasterNotInLibrary { which: String },

    /// Wraps an underlying Symbolica error.
    #[error("symbolica error: {0}")]
    Symbolica(String),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unsupported_loop_order_message_names_the_count() {
        let e = OneLoopError::UnsupportedLoopOrder { found: 2 };
        let msg = e.to_string();
        assert!(msg.contains("2"));
        assert!(msg.to_lowercase().contains("loop"));
    }

    #[test]
    fn master_not_in_library_names_the_config() {
        let e = OneLoopError::MasterNotInLibrary {
            which: "C0(massive triangle)".to_string(),
        };
        assert!(e.to_string().contains("C0(massive triangle)"));
    }
}
