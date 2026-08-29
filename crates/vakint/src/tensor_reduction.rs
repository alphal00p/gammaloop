use ::rustred::algebra::CoefficientContextError;
use ::rustred::family::IntegralFamilyError;
use ::rustred::family::isp::IspCompletionError;
use ::rustred::family::presentation::FamilyPresentationError;
use ::rustred::tensor::TensorError;
use symbolica::atom::{Atom, AtomView};
use thiserror::Error;

use crate::{Vakint, VakintError, VakintExpression, VakintSettings};

mod rustred;

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
    /// The topology matcher could not identify this term.
    #[error("RustRed tensor term {term} has no topology match: {integral}")]
    RustRedUnrecognizedIntegral { term: usize, integral: String },
    /// A matched family has no physical vacuum basis for this native slice.
    #[error(
        "RustRed tensor term {term} has {loop_count} loops and {propagator_count} physical propagators; the native vacuum bridge needs at least one loop and one propagator"
    )]
    RustRedUnsupportedFamily {
        term: usize,
        loop_count: usize,
        propagator_count: usize,
    },
    /// Unknown topologies do not yet carry enough authenticated conventions.
    #[error("RustRed tensor term {term} matched an unknown topology")]
    RustRedUnknownTopology { term: usize },
    /// A canonical propagator could not be recovered after matching.
    #[error("RustRed tensor term {term} is missing canonical propagator {propagator}")]
    RustRedMissingPropagator { term: usize, propagator: usize },
    /// An explicit routing is admitted only when Vakint's matcher exports a
    /// complete replayable loop-basis witness.
    #[error(
        "RustRed tensor term {term} uses an explicit {loop_count}-loop propagator routing that cannot be replayed from a complete matcher basis witness; use a registered short topology instead"
    )]
    RustRedExplicitRoutingUnsupported { term: usize, loop_count: usize },
    /// A physical propagator momentum is not a nonzero integer-linear
    /// combination of the canonical loop basis.
    #[error(
        "RustRed tensor term {term} has unsupported momentum {momentum} in propagator {propagator}; expected a nonzero integer-linear combination of canonical loop momenta"
    )]
    RustRedUnsupportedMomentum {
        term: usize,
        propagator: usize,
        momentum: String,
    },
    /// Adapter-private momentum labels are reserved and cannot be supplied by
    /// a caller, even though Symbolica function heads are globally spellable.
    #[error(
        "RustRed tensor term {term} contains adapter-reserved momentum label head {head}; reserved labels cannot appear in input"
    )]
    RustRedReservedMomentumLabel { term: usize, head: String },
    /// The vacuum lane accepts one common nonzero exact scalar mass squared.
    #[error(
        "RustRed tensor term {term} has unsupported mass squared {mass} in propagator {propagator}; expected the same nonzero exact number or symbol on every physical propagator"
    )]
    RustRedUnsupportedMass {
        term: usize,
        propagator: usize,
        mass: String,
    },
    /// Propagator powers must be exact machine-sized integers at this boundary.
    #[error(
        "RustRed tensor term {term} has unsupported power {power} in propagator {propagator}; expected an exact integer"
    )]
    RustRedUnsupportedPower {
        term: usize,
        propagator: usize,
        power: String,
    },
    /// The configured epsilon spelling could not be parsed exactly.
    #[error("RustRed could not parse Vakint's epsilon symbol {symbol:?}: {detail}")]
    RustRedInvalidDimension { symbol: String, detail: String },
    /// A required projector guard vanished under Vakint's dimension map.
    #[error("RustRed tensor term {term} has singular specialized dimension {dimension}")]
    RustRedSingularDimension { term: usize, dimension: String },
    /// RustRed's exact coefficient domain could not be constructed.
    #[error("RustRed coefficient-domain construction failed: {source}")]
    RustRedCoefficientContext {
        #[source]
        source: CoefficientContextError,
    },
    /// RustRed rejected the affine integral-family bridge.
    #[error("RustRed tensor term {term} rejected its affine family: {source}")]
    RustRedFamily {
        term: usize,
        #[source]
        source: IntegralFamilyError,
    },
    /// RustRed could not complete the independent physical rows with ISPs.
    #[error("RustRed tensor term {term} could not complete its vacuum family: {source}")]
    RustRedIspCompletion {
        term: usize,
        #[source]
        source: IspCompletionError,
    },
    /// RustRed rejected the physical presentation supplied by Vakint.
    #[error("RustRed tensor term {term} rejected its family presentation: {source}")]
    RustRedPresentation {
        term: usize,
        #[source]
        source: FamilyPresentationError,
    },
    /// RustRed's native projector rejected this numerator.
    #[error("RustRed tensor projection failed for term {term}: {source}")]
    RustRedTensor {
        term: usize,
        #[source]
        source: TensorError,
    },
    /// The native result crossed a boundary not representable by this adapter.
    #[error("RustRed tensor term {term} produced unsupported output: {detail}")]
    RustRedUnsupportedOutput { term: usize, detail: String },
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
    expression: &mut VakintExpression,
    vakint: &Vakint,
    settings: &VakintSettings,
    options: RustRedOptions,
) -> Result<(), TensorReductionError> {
    rustred::reduce(expression, vakint, settings, options)
}
