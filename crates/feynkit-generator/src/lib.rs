//! Model-driven Feynman-diagram generation without GammaLoop application state.

#![forbid(unsafe_code)]

mod generation;
mod grouping;
mod options;
mod process;

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
