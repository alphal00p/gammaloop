use symbolica::atom::{Atom, AtomView};
use thiserror::Error;

use crate::{Vakint, VakintError, VakintExpression, VakintSettings};

/// Options for the native RustRed tensor-reduction backend.
///
/// This type is intentionally opaque so that backend controls can be added
/// without changing Vakint's existing settings layout.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RustRedOptions {
    _private: (),
}

impl RustRedOptions {
    /// Creates options for the native backend.
    pub const fn new() -> Self {
        Self { _private: () }
    }
}

/// Backend used for Lorentz tensor reduction.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum TensorReductionMode {
    /// Preserve Vakint's existing FORM-backed tensor reduction.
    #[default]
    Form,
    /// Use RustRed's native Rust and Symbolica implementation.
    RustRed(RustRedOptions),
}

/// Failure reported by the additive tensor-backend interface.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum TensorReductionError {
    /// The existing Vakint input or FORM path failed.
    #[error(transparent)]
    Vakint(#[from] VakintError),
    /// The selected RustRed dependency does not yet expose the native tensor service.
    #[error(
        "RustRed tensor reduction is unavailable because the current RustRed core dependency does not yet expose its native tensor-reduction service"
    )]
    RustRedUnavailable,
}

/// A tensor-reduction operation configured independently of [`VakintSettings`].
///
/// Construct this through [`Vakint::tensor_reducer`]. The default mode is
/// [`TensorReductionMode::Form`], preserving the historical behavior.
#[derive(Debug)]
pub struct TensorReducer<'a> {
    pub(crate) vakint: &'a Vakint,
    pub(crate) settings: &'a VakintSettings,
    pub(crate) mode: TensorReductionMode,
}

impl TensorReducer<'_> {
    /// Returns the backend selected for this operation.
    pub const fn selected_mode(&self) -> &TensorReductionMode {
        &self.mode
    }

    /// Selects a backend without changing Vakint's persistent settings.
    pub const fn mode(mut self, mode: TensorReductionMode) -> Self {
        self.mode = mode;
        self
    }

    /// Reduces every term in the input expression through the selected backend.
    pub fn reduce(self, input: AtomView) -> Result<Atom, TensorReductionError> {
        let mut expression = VakintExpression::try_from(input)?;
        expression.tensor_reduce_in_mode(self.vakint, self.settings, self.mode)?;
        Ok(expression.into())
    }
}

pub(crate) fn reduce_with_rustred(
    _expression: &mut VakintExpression,
    _vakint: &Vakint,
    _settings: &VakintSettings,
    _options: RustRedOptions,
) -> Result<(), TensorReductionError> {
    // RustRed does not expose its native tensor service yet. Keep this boundary
    // explicit: the future implementation belongs in RustRed core and this
    // adapter must never contain a second CAS algorithm or fall back to FORM.
    Err(TensorReductionError::RustRedUnavailable)
}
