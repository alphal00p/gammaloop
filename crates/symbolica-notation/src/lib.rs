//! Shared display notation for Symbolica symbols, below tensor and graph renderers.
//!
//! Structured labels retain the existing portable math-display encoding. The
//! process-local notation owner is independent of symbol serialization; a later
//! opt-in portable attachment will export notation explicitly. [`Notation::register`]
//! stores immutable shared notation in this process without changing symbol tags
//! or user data; [`Notation::lookup`] retrieves it. Renderers do not consult the
//! registry yet.
//!
//! ```
//! use symbolica_notation::{IndexDisplay, IndexDisplayError};
//! let label = IndexDisplay::symbol("mu")?;
//! assert_eq!(label.to_mathml().as_deref(), Some("<mi>μ</mi>"));
//! # Ok::<(), IndexDisplayError>(())
//! ```
mod index_display;
mod mathml;
mod registry;

pub use index_display::{IndexDisplay, IndexDisplayError};
pub use mathml::MathMl;
pub use registry::{DisplayText, DisplayTextError, Notation, RegistryError, TextOverrides};
