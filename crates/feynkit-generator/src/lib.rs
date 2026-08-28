//! Model-driven Feynman-diagram generation without GammaLoop application state.

#![forbid(unsafe_code)]

use symbolica::{atom::Symbol, initialize};

mod generation;
mod grouping;
mod options;
mod process;

// Symbolica does not permit adding tags after a bare symbol with the same name
// has been registered. Claim the public momentum head during global state
// initialization so parsing and generation always agree on its tensor type.
// Use Spenso's canonical tag names directly: FeynKit must interoperate with
// whichever Spenso instance the host embeds instead of linking its own copy.
initialize!(|| {
    let _ = momentum_symbol();
});

pub(crate) fn momentum_symbol() -> Symbol {
    symbolica::symbol!(
        "FeynKit::Momentum",
        tags = ["spenso::tensor", "spenso::rank1"]
    )
}

pub use generation::{
    DiagramGroup, GenerationError, GenerationReport, GenerationResult, Generator, GroupMember,
    UnresolvedCutContent, unresolved_cut_content,
};
pub use grouping::GroupingError;
pub use options::{
    CancellationToken, FilterScope, GenerationControl, GenerationFilter, GenerationFilterKind,
    GenerationOptions, GenerationProgress, GraphGroupingOptions, NumeratorGrouping,
    SelfEnergyFilterOptions, SewnFilterOptions, SnailFilterOptions, TadpoleFilterOptions,
};
pub use process::{
    GenerationType, ParticleSelector, Process, ProcessError, SelectorError, VertexSelector,
};
