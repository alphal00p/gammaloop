//! Shared display notation for Symbolica symbols, below tensor and graph renderers.
//!
//! Structured labels retain the existing portable math-display encoding. The
//! process-local notation owner is independent of symbol serialization; a later
//! opt-in portable attachment will export notation explicitly.
//!
//! ```
//! use symbolica_notation::IndexDisplay;
//! let label = IndexDisplay::symbol("mu")?;
//! assert_eq!(label.to_mathml().as_deref(), Some("<mi>μ</mi>"));
//! # Ok::<(), symbolica_notation::IndexDisplayError>(())
//! ```
mod index_display;
mod mathml;

pub use index_display::{IndexDisplay, IndexDisplayError};
pub use mathml::MathMl;
