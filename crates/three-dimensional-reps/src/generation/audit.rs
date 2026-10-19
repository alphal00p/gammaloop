//! Invocation evidence for tests that run in their own process.
//!
//! Recording is process-wide so Rayon workers cannot escape the audit. There is
//! deliberately no reset or mutable enable flag: an isolated test snapshots all
//! calls in its process, while ordinary builds do not compile this module.

use std::{sync::Mutex, thread::ThreadId};

use super::{CffGenerationContext, Generate3DExpressionOptions, ParsedGraph};

#[derive(Clone, Debug)]
pub struct LtdInvocation {
    pub loops: usize,
    pub denominator_edges: usize,
    pub repeated_poles: bool,
    pub context: CffGenerationContext,
    pub thread: ThreadId,
}

#[derive(Clone, Debug, Default)]
pub struct GenerationInvocations {
    pub cff_builders: usize,
    pub ltd: Vec<LtdInvocation>,
    pub ltd_residue_selections: usize,
}

static INVOCATIONS: Mutex<GenerationInvocations> = Mutex::new(GenerationInvocations {
    cff_builders: 0,
    ltd: Vec::new(),
    ltd_residue_selections: 0,
});

impl GenerationInvocations {
    pub fn snapshot() -> Self {
        INVOCATIONS.lock().unwrap().clone()
    }
}

pub(crate) fn record_cff() {
    INVOCATIONS.lock().unwrap().cff_builders += 1;
}

pub(super) fn record_ltd(
    graph: &ParsedGraph,
    options: &Generate3DExpressionOptions,
    repeated_poles: bool,
) {
    INVOCATIONS.lock().unwrap().ltd.push(LtdInvocation {
        loops: graph.loop_names.len(),
        denominator_edges: graph.denominator_internal_edge_ids().len(),
        repeated_poles,
        context: options.cff_generation_context,
        thread: std::thread::current().id(),
    });
}

pub(crate) fn record_ltd_residue_selection() {
    INVOCATIONS.lock().unwrap().ltd_residue_selections += 1;
}
