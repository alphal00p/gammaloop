//! Snapshot-bound inspection of registered RustRed constants and masters.

use std::{collections::BTreeMap, fmt, str::FromStr, sync::Arc};

use feynkit_graph::FeynmanDiagram;
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    domains::float::Float,
};
use thiserror::Error;

use crate::{
    VakintError,
    precision::{
        AvailablePrecision, CoefficientPrecision, NumericalCoefficient, PrecisionSource,
        PrecisionWarning,
    },
    rustred_evaluation::masters::MasterCatalog,
    utils::vakint_macros::vk_parse,
};

pub(crate) mod constants;
mod graph;
mod registry;

#[cfg(test)]
mod acceptance_tests;
#[cfg(test)]
mod parameter_tests;

/// A lookup failure, distinct from missing numerical precision or epsilon orders.
#[derive(Debug, Error)]
pub enum InspectionError {
    #[error("unknown RustRed constant {symbol}")]
    UnknownConstant { symbol: String },
    #[error("unknown RustRed master {symbol}")]
    UnknownMaster { symbol: String },
    #[error("ambiguous RustRed master {symbol}; qualified matches: {candidates:?}")]
    AmbiguousMaster {
        symbol: String,
        candidates: Vec<String>,
    },
    #[error("invalid {kind} symbol {symbol}")]
    InvalidSymbol { kind: &'static str, symbol: String },
    #[error("master {master} has no usable native graph: {detail}")]
    UnavailableGraph { master: String, detail: String },
    #[error("registered constant {symbol} does not have a real numerical value")]
    NonRealConstant { symbol: String },
    #[error(
        "parameter {symbol} would override a registered constant or master source; use an ordinary parameter symbol"
    )]
    ReservedParameter { symbol: String },
}

/// Parameters supply independent inputs; registered constants retain their definitions.
/// Run this before exact replacements as well as before numerical normalization.
pub(crate) fn validate_parameter(
    identity: AtomView<'_>,
    catalog: Option<&MasterCatalog>,
) -> Result<(), VakintError> {
    crate::Vakint::initialize_vakint_symbols();
    let reserved_source = crate::precision::legacy::is_source(identity)
        || crate::precision::snapshot_reference(identity);
    let native_period = matches!(identity, AtomView::Fun(function)
        if function.get_symbol() == symbolica::transcendental::hpl_at_sixth_root()
            || function.get_symbol() == symbolica::transcendental::gpl_at_one());
    let registered_period = match catalog {
        Some(catalog) => catalog.period_binding(identity)?.is_some(),
        None => false,
    };
    if reserved_source
        || native_period
        || constants::builtin_definition(identity)?.is_some()
        || registered_period
    {
        return Err(InspectionError::ReservedParameter {
            symbol: identity.to_canonical_string(),
        }
        .into());
    }
    Ok(())
}

macro_rules! symbolic_identity {
    ($name:ident, $kind:literal) => {
        #[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(Atom);

        impl $name {
            /// Canonical Symbolica identity, including any family or function arguments.
            pub fn as_atom(&self) -> &Atom {
                &self.0
            }
        }

        impl TryFrom<Atom> for $name {
            type Error = VakintError;

            fn try_from(value: Atom) -> Result<Self, Self::Error> {
                if matches!(value.as_view(), AtomView::Var(_) | AtomView::Fun(_)) {
                    Ok(Self(value))
                } else {
                    Err(InspectionError::InvalidSymbol {
                        kind: $kind,
                        symbol: value.to_canonical_string(),
                    }
                    .into())
                }
            }
        }

        impl FromStr for $name {
            type Err = VakintError;

            fn from_str(value: &str) -> Result<Self, Self::Err> {
                crate::Vakint::initialize_vakint_symbols();
                let atom = vk_parse!(value).map_err(|_| InspectionError::InvalidSymbol {
                    kind: $kind,
                    symbol: value.to_owned(),
                })?;
                Self::try_from(atom)
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str(&self.0.to_canonical_string())
            }
        }
    };
}

symbolic_identity!(ConstantSymbol, "constant");
symbolic_identity!(MasterSymbol, "master");

/// Where the inspected value came from. Versions alone do not identify mutable checkouts.
#[derive(Clone, Debug)]
pub struct DataProvenance {
    pub(crate) version: String,
    pub(crate) source_sha256: String,
    pub(crate) master_id: Option<String>,
    pub(crate) epsilon_power: Option<i64>,
}

impl DataProvenance {
    pub fn version(&self) -> &str {
        &self.version
    }
    pub fn source_sha256(&self) -> &str {
        &self.source_sha256
    }
    pub fn master_id(&self) -> Option<&str> {
        self.master_id.as_deref()
    }
    pub fn epsilon_power(&self) -> Option<i64> {
        self.epsilon_power
    }
}

/// An exact formula may depend on finite periods; availability records that distinction.
#[derive(Clone, Debug)]
pub enum MasterSymbolicExpression {
    Expression(Atom),
    NumericalOnly(u32),
}

#[derive(Clone, Debug)]
pub struct VakintConstant {
    symbol: ConstantSymbol,
    requested_decimal_precision: u32,
    available_precision: AvailablePrecision,
    numerical: NumericalCoefficient,
    warnings: Vec<PrecisionWarning>,
    provenance: DataProvenance,
    // Keep finite period bindings alive when the engine activates another release.
    _snapshot: Arc<MasterCatalog>,
}

impl VakintConstant {
    pub fn requested_decimal_precision(&self) -> u32 {
        self.requested_decimal_precision
    }
    pub fn symbol(&self) -> &ConstantSymbol {
        &self.symbol
    }
    pub fn available_precision(&self) -> AvailablePrecision {
        self.available_precision
    }
    pub fn to_atom(&self) -> Result<Atom, VakintError> {
        self.numerical.to_atom()
    }
    pub fn numerical_value(&self) -> &Float {
        &self.numerical.value().re
    }
    pub fn precision(&self) -> &CoefficientPrecision {
        self.numerical.precision()
    }
    pub fn warnings(&self) -> &[PrecisionWarning] {
        &self.warnings
    }
    pub fn provenance(&self) -> &DataProvenance {
        &self.provenance
    }
}

#[derive(Clone, Debug)]
pub struct MasterCoefficient {
    symbol: Atom,
    requested_decimal_precision: u32,
    epsilon_power: i64,
    available_precision: AvailablePrecision,
    numerical: NumericalCoefficient,
    warnings: Vec<PrecisionWarning>,
    provenance: DataProvenance,
    symbolic_expression: MasterSymbolicExpression,
    _snapshot: Arc<MasterCatalog>,
}

impl MasterCoefficient {
    pub fn requested_decimal_precision(&self) -> u32 {
        self.requested_decimal_precision
    }
    pub fn symbol(&self) -> &Atom {
        &self.symbol
    }
    pub fn epsilon_power(&self) -> i64 {
        self.epsilon_power
    }
    pub fn available_precision(&self) -> AvailablePrecision {
        self.available_precision
    }
    pub fn to_atom(&self) -> Result<Atom, VakintError> {
        self.numerical.to_atom()
    }
    pub fn numerical_value(&self) -> &Float {
        &self.numerical.value().re
    }
    pub fn precision(&self) -> &CoefficientPrecision {
        self.numerical.precision()
    }
    pub fn warnings(&self) -> &[PrecisionWarning] {
        &self.warnings
    }
    pub fn provenance(&self) -> &DataProvenance {
        &self.provenance
    }
    pub fn symbolic_expression(&self) -> &MasterSymbolicExpression {
        &self.symbolic_expression
    }
}

#[derive(Clone, Debug)]
pub struct VakintMaster {
    symbol: MasterSymbol,
    requested_decimal_precision: u32,
    master_id: String,
    family_sha256: String,
    loop_count: usize,
    normalization: String,
    provenance: DataProvenance,
    graph: FeynmanDiagram,
    coefficients: BTreeMap<i64, MasterCoefficient>,
    _snapshot: Arc<MasterCatalog>,
}

impl VakintMaster {
    pub fn requested_decimal_precision(&self) -> u32 {
        self.requested_decimal_precision
    }
    pub fn symbol(&self) -> &MasterSymbol {
        &self.symbol
    }
    pub fn master_id(&self) -> &str {
        &self.master_id
    }
    pub fn family_sha256(&self) -> &str {
        &self.family_sha256
    }
    pub fn loop_count(&self) -> usize {
        self.loop_count
    }
    pub fn normalization(&self) -> &str {
        &self.normalization
    }
    pub fn provenance(&self) -> &DataProvenance {
        &self.provenance
    }
    pub fn graph(&self) -> &FeynmanDiagram {
        &self.graph
    }
    pub fn coefficients(&self) -> &BTreeMap<i64, MasterCoefficient> {
        &self.coefficients
    }
    pub fn is_computed(&self) -> bool {
        !self.coefficients.is_empty()
    }
}

/// Decimal source evidence stays separate from arithmetic working precision.
#[derive(Clone, Debug)]
pub(crate) struct FiniteSource {
    pub value: String,
    pub significant_digits: u32,
    pub absolute_error: Option<String>,
    pub error_is_bounded: bool,
    pub identity: String,
    pub master_id: Option<String>,
    pub epsilon_power: Option<i64>,
    pub data_version: String,
    pub source_sha256: String,
}

#[derive(Clone, Debug)]
pub(crate) struct PreparedExpression {
    pub expression: Atom,
    pub sources: BTreeMap<Atom, FiniteSource>,
    pub origins: Vec<PrecisionSource>,
}

impl PreparedExpression {
    pub(crate) fn available_precision(&self) -> AvailablePrecision {
        self.sources
            .values()
            .map(|source| source.significant_digits)
            .min()
            .map_or(
                AvailablePrecision::ArbitraryPrecision,
                AvailablePrecision::FixedPrecision,
            )
    }
}
