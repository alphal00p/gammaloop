//! Model-driven Feynman-diagram generation without GammaLoop application state.

#![forbid(unsafe_code)]

use spenso::{network::tags::SPENSO_TAG, symbolica_init::in_symbolica_initializer};
use symbolica::{atom::Symbol, initialize};

mod generation;
mod grouping;
mod options;
mod process;

// Symbolica does not permit adding tags after a bare symbol with the same name
// has been registered. Claim the public momentum head during global state
// initialization so parsing and generation always agree on its tensor type.
initialize!(
    || {
        in_symbolica_initializer(|| {
            let _ = SPENSO_TAG
                .force_in_initializer()
                .rank_one_tensor_symbol("FeynKit::Momentum");
        });
    },
    "spenso"
);

pub(crate) fn momentum_symbol() -> Symbol {
    SPENSO_TAG.rank_one_tensor_symbol("FeynKit::Momentum")
}

pub use generation::{
    DiagramGroup, GenerationError, GenerationReport, GenerationResult, Generator, GroupMember,
};
pub use grouping::GroupingError;
pub use options::{
    CancellationToken, FilterScope, GenerationControl, GenerationFilter, GenerationFilterKind,
    GenerationOptions, GenerationProgress, GraphGroupingOptions, NumeratorGrouping,
    SelfEnergyFilterOptions, SewnFilterOptions, SnailFilterOptions, TadpoleFilterOptions,
};
pub use process::{GenerationType, ParticleSelector, Process, ProcessError};
