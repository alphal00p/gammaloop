//! Model-driven Feynman-diagram generation without GammaLoop application state.

#![forbid(unsafe_code)]

// Unit tests share the case table of the generation oracle tests, which names this crate.
#[cfg(test)]
extern crate self as feynkit_generator;

mod generation;
mod grouping;
mod options;
mod process;

use feynkit_graph::momentum_symbol;

pub use generation::{
    DiagramGroup, ExternalLabel, ExternalNode, FeynmanDiagramRulesExt, GenerationDiagramView,
    GenerationDressingView, GenerationError, GenerationReport, GenerationResult,
    GenerationSkeletonView, GroupMember, SourceIdentity, UnresolvedCutContent,
    unresolved_cut_content,
};
pub use grouping::GroupingError;
pub use options::{
    CancellationToken, FilterScope, GenerationControl, GenerationFilter, GenerationFilterKind,
    GenerationOptions, GenerationOutput, GenerationProgress, GraphGroupingOptions,
    NumeratorGrouping, SelfEnergyFilterOptions, SewnFilterOptions, SnailFilterOptions,
    TadpoleFilterOptions,
};
pub use process::{
    GenerationType, ModelProcessExt, ParticleSelector, Process, ProcessError, SelectorError,
    VertexSelector,
};
