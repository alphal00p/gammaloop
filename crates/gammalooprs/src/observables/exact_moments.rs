//! Shared exact represented-sample statistics owner.
//!
//! Re-exported here so existing histogram consumers retain the same API and
//! numerical types; independent consumers can depend on `gammaloop-statistics`
//! without importing GammaLoop or Symbolica.

pub use gammaloop_statistics::*;
