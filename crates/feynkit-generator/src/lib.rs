//! Model-driven Feynman-diagram generation without GammaLoop application state.

#![forbid(unsafe_code)]

// Unit tests share the case table of the generation oracle tests, which names this crate.
#[cfg(test)]
extern crate self as feynkit_generator;

mod generation;
mod grouping;
mod options;
mod parser;
mod process;

use feynkit_graph::momentum_symbol;

pub use generation::{
    BuiltDiagramSet, Counts, DiagramGroup, DiagramSet, ExternalLabel, ExternalNode,
    FeynmanDiagramRulesExt, GenerationDiagramView, GenerationDressingView, GenerationError,
    GenerationReport, GenerationSkeletonView, GroupMember, LoopCounts, PreparedDiagram,
    PreparedEdge, PreparedEndpoint, SourceIdentity, UnresolvedCutContent, unresolved_cut_content,
};
pub use grouping::GroupingError;
pub use options::{
    Bridges, CancellationToken, Classes, CutAmplitude, Cuts, GenerationControl, GenerationProgress,
    GraphGroupingOptions, Grouping, Loops, Naming, NumeratorOptions, Orders, Range, Request,
    Restrict, Run, Selection, Stage, Symmetrize, TopologyRules, Veto,
};
pub use parser::RequestParseError;
pub use process::{
    GenerationType, ModelProcessExt, ParticleSelector, Process, ProcessError, SelectorError,
    VertexSelector,
};
