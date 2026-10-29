//! Shared exact represented-sample statistics owner.
//!
//! Re-exported here so existing histogram consumers retain the same API and
//! numerical types; independent consumers can depend on `gammaloop-statistics`
//! without importing GammaLoop or Symbolica.

pub use gammaloop_statistics::{
    ExactF64Moments, ExactMomentError, MAX_CANONICAL_OBJECT_BYTES, MAX_EXACT_OBJECTS,
    MAX_JSON_FRAME_BYTES, MAX_SAMPLE_COUNT, MomentOrder, ReportPolicy, certified_round_f64,
    certified_sqrt_f64,
};
