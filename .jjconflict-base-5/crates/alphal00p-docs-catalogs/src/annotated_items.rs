//! Source-backed proc-macro descriptors for the supported Rust surfaces.
//!
//! Each marker is consumed by the attribute macro. At compile time the macro
//! parses the named product source item, captures its real signature, members,
//! and existing Rustdoc, and emits a static descriptor plus a dependency edge.
//! Product crates therefore remain independently publishable: only this
//! internal adapter depends on the documentation crates.

use alphal00p_docs_schema::DocItem;

#[alphal00p_docs::trait_item(
    id = "GammaLoopContext",
    title = "GammaLoopContext",
    summary = "Combines the state-map and model capabilities required by GammaLoop computations.",
    format = "rust-markdown",
    source = "crates/gammalooprs/src/lib.rs",
    source_id = "gammalooprs::GammaLoopContext"
)]
fn gammalooprs_gamma_loop_context() {}

#[alphal00p_docs::trait_item(
    id = "HasModel",
    title = "HasModel",
    summary = "Provides shared access to the physics model owned by a computation context.",
    format = "rust-markdown",
    source = "crates/gammalooprs/src/lib.rs",
    source_id = "gammalooprs::HasModel"
)]
fn gammalooprs_has_model() {}

#[alphal00p_docs::ty(
    id = "GammaLoopContextContainer",
    title = "Concrete GammaLoop context",
    summary = "Pairs a Symbolica state map with the physics model required by GammaLoop computations.",
    format = "rust-markdown",
    source = "crates/gammalooprs/src/lib.rs",
    source_id = "gammalooprs::GammaLoopContextContainer"
)]
fn gammalooprs_gamma_loop_context_container() {}

#[alphal00p_docs::trait_item(
    id = "HasIntegrand",
    title = "Integrand contract",
    summary = "Defines grid construction, sample evaluation, dimensions, and result merging for an integrand.",
    format = "rust-markdown",
    source = "crates/gammalooprs/src/integrands/mod.rs",
    source_id = "gammalooprs::integrands::HasIntegrand"
)]
fn gammalooprs_has_integrand() {}

#[alphal00p_docs::ty(
    id = "Integrand",
    title = "Runtime integrand",
    summary = "Dispatches evaluation and observable handling across built-in and generated process integrands.",
    format = "rust-markdown",
    source = "crates/gammalooprs/src/integrands/mod.rs",
    source_id = "gammalooprs::integrands::Integrand"
)]
fn gammalooprs_integrand() {}

#[alphal00p_docs::ty(
    id = "GlobalSettings",
    title = "Global settings",
    summary = "Configures logging, process generation, worker counts, and other execution-wide GammaLoop defaults.",
    kind = "setting",
    format = "rust-markdown",
    source = "crates/gammalooprs/src/settings/mod.rs",
    source_id = "gammalooprs::settings::GlobalSettings"
)]
fn gammalooprs_global_settings() {}

#[alphal00p_docs::ty(
    id = "RuntimeSettings",
    title = "Runtime settings",
    summary = "Collects kinematics, integration, subtraction, observables, and evaluation settings for one integrand.",
    kind = "setting",
    format = "rust-markdown",
    source = "crates/gammalooprs/src/settings/mod.rs",
    source_id = "gammalooprs::settings::RuntimeSettings"
)]
fn gammalooprs_runtime_settings() {}

#[alphal00p_docs::ty(
    id = "SingleSampleEvaluationResult",
    title = "Single-sample result",
    summary = "Returns one f64 evaluation together with the observable snapshot produced for that sample.",
    format = "rust-markdown",
    source = "crates/gammalooprs/src/integrands/evaluation.rs",
    source_id = "gammalooprs::integrands::evaluation::SingleSampleEvaluationResult"
)]
fn gammalooprs_single_sample_evaluation_result() {}

#[alphal00p_docs::ty(
    id = "BatchSampleEvaluationResult",
    title = "Batch evaluation result",
    summary = "Returns f64 sample evaluations and the observable snapshot accumulated across the complete batch.",
    format = "rust-markdown",
    source = "crates/gammalooprs/src/integrands/evaluation.rs",
    source_id = "gammalooprs::integrands::evaluation::BatchSampleEvaluationResult"
)]
fn gammalooprs_batch_sample_evaluation_result() {}

#[alphal00p_docs::ty(
    id = "PreciseSingleSampleEvaluationResult",
    title = "Precise single-sample result",
    summary = "Retains the active numerical precision for one evaluation and its observable snapshot.",
    format = "rust-markdown",
    source = "crates/gammalooprs/src/integrands/evaluation.rs",
    source_id = "gammalooprs::integrands::evaluation::PreciseSingleSampleEvaluationResult"
)]
fn gammalooprs_precise_single_sample_evaluation_result() {}

#[alphal00p_docs::ty(
    id = "PreciseBatchSampleEvaluationResult",
    title = "Precise batch result",
    summary = "Retains the active numerical precision for every sample in a batch and its observable snapshot.",
    format = "rust-markdown",
    source = "crates/gammalooprs/src/integrands/evaluation.rs",
    source_id = "gammalooprs::integrands::evaluation::PreciseBatchSampleEvaluationResult"
)]
fn gammalooprs_precise_batch_sample_evaluation_result() {}

#[alphal00p_docs::ty(
    id = "HistogramSnapshot",
    title = "Mergeable histogram snapshot",
    summary = "Stores raw per-bin statistics, underflow, overflow, and metadata in a serializable snapshot.",
    format = "rust-markdown",
    source = "crates/gammalooprs/src/observables/mod.rs",
    source_id = "gammalooprs::observables::HistogramSnapshot"
)]
fn gammalooprs_histogram_snapshot() {}

#[alphal00p_docs::func(
    id = "HistogramSnapshot::merge",
    title = "Merge histogram snapshots",
    summary = "Combines compatible raw-stat snapshots without losing reconstructible accumulator state.",
    format = "rust-markdown",
    source = "crates/gammalooprs/src/observables/mod.rs",
    source_id = "gammalooprs::observables::HistogramSnapshot::merge"
)]
fn gammalooprs_histogram_snapshot_merge() {}

#[alphal00p_docs::func(
    id = "HistogramSnapshot::rebin",
    title = "Rebin a histogram snapshot",
    summary = "Combines contiguous continuous bins while preserving raw statistics and rejecting invalid factors.",
    format = "rust-markdown",
    source = "crates/gammalooprs/src/observables/mod.rs",
    source_id = "gammalooprs::observables::HistogramSnapshot::rebin"
)]
fn gammalooprs_histogram_snapshot_rebin() {}

#[alphal00p_docs::ty(
    id = "HistogramAccumulatorState",
    title = "Histogram accumulator",
    summary = "Accumulates statistically independent samples before producing mergeable histogram snapshots.",
    format = "rust-markdown",
    source = "crates/gammalooprs/src/observables/mod.rs",
    source_id = "gammalooprs::observables::HistogramAccumulatorState"
)]
fn gammalooprs_histogram_accumulator_state() {}

#[alphal00p_docs::func(
    id = "HistogramAccumulatorState::snapshot",
    title = "Snapshot a histogram accumulator",
    summary = "Captures the accumulator's current raw statistics in a serializable histogram snapshot.",
    format = "rust-markdown",
    source = "crates/gammalooprs/src/observables/mod.rs",
    source_id = "gammalooprs::observables::HistogramAccumulatorState::snapshot"
)]
fn gammalooprs_histogram_accumulator_state_snapshot() {}

#[alphal00p_docs::func(
    id = "set_interrupt_handler",
    title = "Install interrupt handling",
    summary = "Installs GammaLoop's process interrupt handler before a long-running integration.",
    format = "rust-markdown",
    source = "crates/gammalooprs/src/lib.rs",
    source_id = "gammalooprs::set_interrupt_handler"
)]
fn gammalooprs_set_interrupt_handler() {}

#[alphal00p_docs::func(
    id = "request_interrupt",
    title = "Request interruption",
    summary = "Signals cooperative cancellation to long-running GammaLoop evaluation and integration code.",
    format = "rust-markdown",
    source = "crates/gammalooprs/src/lib.rs",
    source_id = "gammalooprs::request_interrupt"
)]
fn gammalooprs_request_interrupt() {}

#[alphal00p_docs::func(
    id = "is_interrupt_requested",
    title = "Inspect interruption state",
    summary = "Reports whether cooperative cancellation has been requested for the current process.",
    format = "rust-markdown",
    source = "crates/gammalooprs/src/lib.rs",
    source_id = "gammalooprs::is_interrupt_requested"
)]
fn gammalooprs_is_interrupt_requested() {}

#[alphal00p_docs::func(
    id = "clear_interrupt_request",
    title = "Clear interruption state",
    summary = "Clears a previous cooperative-cancellation request before beginning new work.",
    format = "rust-markdown",
    source = "crates/gammalooprs/src/lib.rs",
    source_id = "gammalooprs::clear_interrupt_request"
)]
fn gammalooprs_clear_interrupt_request() {}

#[alphal00p_docs::ty(
    id = "StateLoadOption",
    title = "State loading options",
    summary = "Selects the state folder, boot commands, settings overrides, logging, and read-only mode.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/lib.rs",
    source_id = "gammaloop_api::StateLoadOption"
)]
fn gammaloop_api_state_load_option() {}

#[alphal00p_docs::ty(
    id = "LoadedState",
    title = "Loaded GammaLoop state",
    summary = "Owns the physics state, run history, settings, and session state returned by the loading facade.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/lib.rs",
    source_id = "gammaloop_api::LoadedState"
)]
fn gammaloop_api_loaded_state() {}

#[alphal00p_docs::func(
    id = "StateLoadOption::load",
    title = "Load a GammaLoop state",
    summary = "Loads the selected state and returns the supported session boundary.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/lib.rs",
    source_id = "gammaloop_api::StateLoadOption::load"
)]
fn gammaloop_api_state_load_option_load() {}

#[alphal00p_docs::func(
    id = "LoadedState::cli_session",
    title = "Create a CLI session",
    summary = "Borrows a loaded state as a command session while keeping state ownership explicit.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/lib.rs",
    source_id = "gammaloop_api::LoadedState::cli_session"
)]
fn gammaloop_api_loaded_state_cli_session() {}

#[alphal00p_docs::ty(
    id = "CliSession",
    title = "Command session",
    summary = "Borrows a loaded state's coordinated state, history, settings, and transient session data.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/session.rs",
    source_id = "gammaloop_api::session::CliSession"
)]
fn gammaloop_api_cli_session() {}

#[alphal00p_docs::func(
    id = "CliSession::execute_command",
    title = "Execute a parsed command",
    summary = "Runs one command through the same stateful execution path used by the CLI and Python facade.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/session.rs",
    source_id = "gammaloop_api::session::CliSession::execute_command"
)]
fn gammaloop_api_cli_session_execute_command() {}

#[alphal00p_docs::ty(
    id = "CommandHistory",
    title = "Parsed command record",
    summary = "Stores a typed command with its optional original spelling for replay and serialization.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/state.rs",
    source_id = "gammaloop_api::state::CommandHistory"
)]
fn gammaloop_api_command_history() {}

#[alphal00p_docs::func(
    id = "CommandHistory::from_raw_string",
    title = "Parse command text",
    summary = "Parses shell-style GammaLoop command text into the typed command representation used by a session.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/state.rs",
    source_id = "gammaloop_api::state::CommandHistory::from_raw_string"
)]
fn gammaloop_api_command_history_from_raw_string() {}

#[alphal00p_docs::ty(
    id = "RunHistory",
    title = "Replayable run history",
    summary = "Combines frozen settings, named command blocks, and ordered commands into a reusable run card.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/state.rs",
    source_id = "gammaloop_api::state::RunHistory"
)]
fn gammaloop_api_run_history() {}

#[alphal00p_docs::func(
    id = "RunHistory::load",
    title = "Load a run card",
    summary = "Loads and validates a TOML, YAML, or JSON run history before it is replayed.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/state.rs",
    source_id = "gammaloop_api::state::RunHistory::load"
)]
fn gammaloop_api_run_history_load() {}

#[alphal00p_docs::ty(
    id = "State",
    title = "GammaLoop physics state",
    summary = "Owns the model, model parameters, processes, generated integrands, and generation summaries.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/state.rs",
    source_id = "gammaloop_api::state::State"
)]
fn gammaloop_api_state() {}

#[alphal00p_docs::func(
    id = "State::get_integrand_info",
    title = "Inspect a generated integrand",
    summary = "Returns structured backend, graph-group, orientation, loop-basis, cut, and threshold metadata.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/state.rs",
    source_id = "gammaloop_api::state::State::get_integrand_info"
)]
fn gammaloop_api_state_get_integrand_info() {}

#[alphal00p_docs::ty(
    id = "IntegrandInfo",
    title = "Structured integrand information",
    summary = "Describes one generated integrand without exposing its internal evaluator representation.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/integrand_info.rs",
    source_id = "gammaloop_api::integrand_info::IntegrandInfo"
)]
fn gammaloop_api_integrand_info() {}

#[alphal00p_docs::ty(
    id = "ImportModel",
    title = "Model import request",
    summary = "Selects a built-in JSON model, a JSON model path, or a UFO directory and controls simplification before replacing the active model.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/commands/import/model.rs",
    source_id = "gammaloop_api::commands::import::model::ImportModel"
)]
fn gammaloop_api_import_model() {}

#[alphal00p_docs::func(
    id = "ImportModel::run",
    title = "Import a model",
    summary = "Parses the model specification, loads its model and parameter card, optionally simplifies it, and installs it in the active state.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/commands/import/model.rs",
    source_id = "gammaloop_api::commands::import::model::ImportModel::run"
)]
fn gammaloop_api_import_model_run() {}

#[alphal00p_docs::ty(
    id = "Generate",
    title = "Process-generation request",
    summary = "Selects cross-section, amplitude, or existing-process generation and whether generated compilation sources are retained.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/commands/generate.rs",
    source_id = "gammaloop_api::commands::generate::Generate"
)]
fn gammaloop_api_generate() {}

#[alphal00p_docs::func(
    id = "Generate::run",
    title = "Generate a process",
    summary = "Generates or reuses graphs, builds the requested integrand, performs configured compilation, and records generation summaries in state.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/commands/generate.rs",
    source_id = "gammaloop_api::commands::generate::Generate::run"
)]
fn gammaloop_api_generate_run() {}

#[alphal00p_docs::ty(
    id = "Integrate",
    title = "Integration request",
    summary = "Selects generated integrands and configures their workspace, targets, resume behavior, sampling correlation, batching, and status renderer.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/commands/integrate.rs",
    source_id = "gammaloop_api::commands::integrate::Integrate"
)]
fn gammaloop_api_integrate() {}

#[alphal00p_docs::func(
    id = "Integrate::run",
    title = "Run an integration",
    summary = "Resolves the selected integrands and targets, prepares or resumes their workspace, executes integration, and returns results with observable snapshots.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/commands/integrate.rs",
    source_id = "gammaloop_api::commands::integrate::Integrate::run"
)]
fn gammaloop_api_integrate_run() {}

#[alphal00p_docs::ty(
    id = "IntegrationOutput",
    title = "Integration output",
    summary = "Returns the final integration result, per-slot observable snapshots, and the workspace path that owns resumable integration state.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/commands/integrate.rs",
    source_id = "gammaloop_api::commands::integrate::IntegrationOutput"
)]
fn gammaloop_api_integration_output() {}

#[alphal00p_docs::ty(
    id = "EvaluateSamples",
    title = "f64 sample-evaluation request",
    summary = "Selects an integrand, coordinates, weights, discrete dimensions, and optional momentum-space channels for a batch.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/commands/evaluate_samples.rs",
    source_id = "gammaloop_api::commands::evaluate_samples::EvaluateSamples"
)]
fn gammaloop_api_evaluate_samples_request() {}

#[alphal00p_docs::func(
    id = "evaluate_sample",
    title = "Evaluate one f64 sample",
    summary = "Runs a one-row request and returns the sample plus its observable snapshot through the f64 output contract.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/commands/evaluate_samples.rs",
    source_id = "gammaloop_api::commands::evaluate_samples::evaluate_sample"
)]
fn gammaloop_api_evaluate_sample() {}

#[alphal00p_docs::func(
    id = "evaluate_samples",
    title = "Evaluate an f64 batch",
    summary = "Runs a batch request and returns per-sample values plus one batch-level observable snapshot.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/commands/evaluate_samples.rs",
    source_id = "gammaloop_api::commands::evaluate_samples::evaluate_samples"
)]
fn gammaloop_api_evaluate_samples() {}

#[alphal00p_docs::ty(
    id = "EvaluateSamplesPrecise",
    title = "Precision-preserving evaluation request",
    summary = "Accepts f64 integration or momentum-space coordinates while preserving the selected precision for internal computation and returned values.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/commands/evaluate_samples.rs",
    source_id = "gammaloop_api::commands::evaluate_samples::EvaluateSamplesPrecise"
)]
fn gammaloop_api_evaluate_samples_precise_request() {}

#[alphal00p_docs::func(
    id = "evaluate_sample_precise",
    title = "Evaluate one precision-preserving sample",
    summary = "Returns a single result in the active f64, f128, or arbitrary-precision representation.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/commands/evaluate_samples.rs",
    source_id = "gammaloop_api::commands::evaluate_samples::evaluate_sample_precise"
)]
fn gammaloop_api_evaluate_sample_precise() {}

#[alphal00p_docs::func(
    id = "evaluate_samples_precise",
    title = "Evaluate a precision-preserving batch",
    summary = "Returns every sample in the active precision together with one batch-level observable snapshot.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/commands/evaluate_samples.rs",
    source_id = "gammaloop_api::commands::evaluate_samples::evaluate_samples_precise"
)]
fn gammaloop_api_evaluate_samples_precise() {}

#[alphal00p_docs::ty(
    id = "StateFolderKind",
    title = "State-folder classification",
    summary = "Distinguishes missing, scratch, unmanifested, saved, and invalid state directories.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/state.rs",
    source_id = "gammaloop_api::state::StateFolderKind"
)]
fn gammaloop_api_state_folder_kind() {}

#[alphal00p_docs::func(
    id = "classify_state_folder",
    title = "Classify a state folder",
    summary = "Checks manifest and required-file evidence without treating empty or log-only folders as saved state.",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/state.rs",
    source_id = "gammaloop_api::state::classify_state_folder"
)]
fn gammaloop_api_classify_state_folder() {}

#[alphal00p_docs::ty(
    id = "CLISettings",
    title = "Command-line settings",
    summary = "Owns the generated command-line, global, session, and state-loading settings.",
    kind = "setting",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/lib.rs",
    source_id = "gammaloop_api::CLISettings"
)]
fn gammaloop_api_cli_settings() {}

#[alphal00p_docs::ty(
    id = "gammaloop",
    name = "gammaloop",
    title = "GammaLoop command",
    summary = "Runs the Clap-derived command tree for state creation, inspection, generation, and evaluation.",
    kind = "command",
    format = "rust-markdown",
    source = "crates/gammaloop-api/src/lib.rs",
    source_id = "gammaloop_api::CLISettings"
)]
fn gammaloop_api_command() {}

#[alphal00p_docs::ty(
    id = "HedgeGraph",
    title = "Half-edge graph",
    summary = "Stores graph topology as paired half-edges with independent edge, vertex, and half-edge payloads.",
    format = "rust-markdown",
    source = "crates/linnet/src/half_edge.rs",
    source_id = "linnet::half_edge::HedgeGraph"
)]
fn linnet_hedge_graph() {}

#[alphal00p_docs::ty(
    id = "HedgeGraphBuilder",
    title = "Graph builder",
    summary = "Builds a valid half-edge graph while assigning typed node, edge, and hedge identifiers.",
    format = "rust-markdown",
    source = "crates/linnet/src/half_edge/builder.rs",
    source_id = "linnet::half_edge::builder::HedgeGraphBuilder"
)]
fn linnet_hedge_graph_builder() {}

#[alphal00p_docs::ty(
    id = "SuBitGraph",
    title = "Bit-vector subgraph",
    summary = "Represents a subgraph as an efficient half-edge inclusion filter.",
    format = "rust-markdown",
    source = "crates/linnet/src/half_edge/subgraph/subset.rs",
    source_id = "linnet::half_edge::subgraph::subset::SuBitGraph"
)]
fn linnet_su_bit_graph() {}

#[alphal00p_docs::ty(
    id = "SimpleTraversalTree",
    title = "Traversal tree",
    summary = "Records parent, root, and traversal information produced by breadth-first or depth-first graph walks.",
    format = "rust-markdown",
    source = "crates/linnet/src/half_edge/tree.rs",
    source_id = "linnet::half_edge::tree::SimpleTraversalTree"
)]
fn linnet_simple_traversal_tree() {}

#[alphal00p_docs::ty(
    id = "DotGraph",
    title = "DOT-backed graph",
    summary = "Combines a half-edge graph with DOT attributes and parser/serializer support.",
    format = "rust-markdown",
    source = "crates/linnet/src/parser/mod.rs",
    source_id = "linnet::parser::DotGraph"
)]
fn linnet_dot_graph() {}

#[alphal00p_docs::ty(
    id = "HedgePair",
    title = "Half-edge pairing",
    summary = "Represents an unpaired external half-edge, a paired full edge, or an edge split across a subgraph boundary.",
    format = "rust-markdown",
    source = "crates/linnet/src/half_edge/involution.rs",
    source_id = "linnet::half_edge::involution::HedgePair"
)]
fn linnet_hedge_pair() {}

#[alphal00p_docs::ty(
    id = "Flow",
    title = "Intrinsic half-edge flow",
    summary = "Marks a half-edge as the intrinsic source or sink within its full edge.",
    format = "rust-markdown",
    source = "crates/linnet/src/half_edge/involution.rs",
    source_id = "linnet::half_edge::involution::Flow"
)]
fn linnet_flow() {}

#[alphal00p_docs::ty(
    id = "Orientation",
    title = "Edge orientation",
    summary = "Records default, reversed, or undirected presentation independently of an edge's intrinsic half-edge flow.",
    format = "rust-markdown",
    source = "crates/linnet/src/half_edge/involution.rs",
    source_id = "linnet::half_edge::involution::Orientation"
)]
fn linnet_orientation() {}

#[alphal00p_docs::trait_item(
    id = "SubSetLike",
    title = "Half-edge subset interface",
    summary = "Defines the core inclusion queries, iteration, joining, and labeling contract shared by half-edge subset representations.",
    format = "rust-markdown",
    source = "crates/linnet/src/half_edge/subgraph.rs",
    source_id = "linnet::half_edge::subgraph::SubSetLike"
)]
fn linnet_subset_like() {}

#[alphal00p_docs::trait_item(
    id = "SubSetOps",
    title = "Half-edge subset operations",
    summary = "Provides in-place and value-returning intersection, union, symmetric-difference, and subtraction operations.",
    format = "rust-markdown",
    source = "crates/linnet/src/half_edge/subgraph.rs",
    source_id = "linnet::half_edge::subgraph::SubSetOps"
)]
fn linnet_subset_ops() {}

#[alphal00p_docs::trait_item(
    id = "Inclusion",
    title = "Subset inclusion queries",
    summary = "Abstracts full containment and non-empty intersection checks between a subgraph and another value.",
    format = "rust-markdown",
    source = "crates/linnet/src/half_edge/subgraph.rs",
    source_id = "linnet::half_edge::subgraph::Inclusion"
)]
fn linnet_inclusion() {}

#[alphal00p_docs::trait_item(
    id = "BaseSubgraph",
    title = "Mutable base subgraph",
    summary = "Combines graph-aware subset behavior and mutation with construction from edge filters or half-edge iterators.",
    format = "rust-markdown",
    source = "crates/linnet/src/half_edge/subgraph.rs",
    source_id = "linnet::half_edge::subgraph::BaseSubgraph"
)]
fn linnet_base_subgraph() {}

#[alphal00p_docs::trait_item(
    id = "SubGraphLike",
    title = "Graph-aware subgraph interface",
    summary = "Extends half-edge subsets with graph coverage, edge counts, external hairs, and DOT rendering.",
    format = "rust-markdown",
    source = "crates/linnet/src/half_edge/subgraph.rs",
    source_id = "linnet::half_edge::subgraph::SubGraphLike"
)]
fn linnet_subgraph_like() {}

#[alphal00p_docs::trait_item(
    id = "NodeStorage",
    title = "Node-storage interface",
    summary = "Defines a graph's node data, transformed storage family, neighbor representation, and neighbor iterator.",
    format = "rust-markdown",
    source = "crates/linnet/src/half_edge/nodestore.rs",
    source_id = "linnet::half_edge::nodestore::NodeStorage"
)]
fn linnet_node_storage() {}

#[alphal00p_docs::trait_item(
    id = "NodeStorageOps",
    title = "Mutable node-storage operations",
    summary = "Adds construction, extraction, node identification, validation, traversal, and data mapping to a node-storage backend.",
    format = "rust-markdown",
    source = "crates/linnet/src/half_edge/nodestore.rs",
    source_id = "linnet::half_edge::nodestore::NodeStorageOps"
)]
fn linnet_node_storage_ops() {}

#[alphal00p_docs::trait_item(
    id = "TensorStructure",
    title = "Tensor structure",
    summary = "Describes ordered slots, dimensions, representations, and contraction compatibility independently of data.",
    format = "rust-markdown",
    source = "crates/spenso/src/structure.rs",
    source_id = "spenso::structure::TensorStructure"
)]
fn spenso_tensor_structure() {}

#[alphal00p_docs::ty(
    id = "DenseTensor",
    title = "Dense tensor",
    summary = "Stores every tensor component in structure-defined flat order.",
    format = "rust-markdown",
    source = "crates/spenso/src/tensors/data/dense.rs",
    source_id = "spenso::tensors::data::dense::DenseTensor"
)]
fn spenso_dense_tensor() {}

#[alphal00p_docs::ty(
    id = "SparseTensor",
    title = "Sparse tensor",
    summary = "Stores only explicit tensor components keyed by flat index.",
    format = "rust-markdown",
    source = "crates/spenso/src/tensors/data/sparse.rs",
    source_id = "spenso::tensors::data::sparse::SparseTensor"
)]
fn spenso_sparse_tensor() {}

#[alphal00p_docs::trait_item(
    id = "Contract",
    title = "Pairwise contraction",
    summary = "Defines pairwise tensor contraction with an explicit settings type and output.",
    format = "rust-markdown",
    source = "crates/spenso/src/contraction.rs",
    source_id = "spenso::contraction::Contract"
)]
fn spenso_contract() {}

#[alphal00p_docs::ty(
    id = "Network",
    title = "Tensor network",
    summary = "Binds symbolic tensor expressions to libraries and an executable contraction graph.",
    format = "rust-markdown",
    source = "crates/spenso/src/network/mod.rs",
    source_id = "spenso::network::Network"
)]
fn spenso_network() {}

#[alphal00p_docs::trait_item(
    id = "ContractionStrategy",
    title = "Network contraction strategy",
    summary = "Defines how a ready network operation is contracted and whether product subgraphs can be rewritten incrementally.",
    format = "rust-markdown",
    source = "crates/spenso/src/network/contract.rs",
    source_id = "spenso::network::contract::ContractionStrategy"
)]
fn spenso_contraction_strategy() {}

#[alphal00p_docs::trait_item(
    id = "ExecutionStrategy",
    title = "Network execution strategy",
    summary = "Defines how an executor repeatedly applies a selected contraction strategy until the strategy's stopping condition is reached.",
    format = "rust-markdown",
    source = "crates/spenso/src/network/mod.rs",
    source_id = "spenso::network::ExecutionStrategy"
)]
fn spenso_execution_strategy() {}

#[alphal00p_docs::ty(
    id = "ExecutionMode",
    title = "Network execution mode",
    summary = "Enumerates sequential, parallel, borrowed-reference, and extracting execution modes.",
    format = "rust-markdown",
    source = "crates/spenso/src/network/mod.rs",
    source_id = "spenso::network::ExecutionMode"
)]
fn spenso_execution_mode() {}

#[alphal00p_docs::ty(
    id = "ExecutionResult",
    title = "Network execution result",
    summary = "Represents a completed network as the multiplicative identity, the additive identity, or a concrete value.",
    format = "rust-markdown",
    source = "crates/spenso/src/network/mod.rs",
    source_id = "spenso::network::ExecutionResult"
)]
fn spenso_execution_result() {}

#[alphal00p_docs::ty(
    id = "TensorNetworkError",
    title = "Tensor-network error",
    summary = "Classifies invalid topology, contraction, library, structure, I/O, and result-shape failures during network processing.",
    format = "rust-markdown",
    source = "crates/spenso/src/network/mod.rs",
    source_id = "spenso::network::TensorNetworkError"
)]
fn spenso_tensor_network_error() {}

#[alphal00p_docs::ty(
    id = "ParseSettings",
    title = "Network parsing settings",
    summary = "Controls scalar precontraction, depth limits, shorthand expansion and inference, scalar-leaf handling, and tensor-head recognition.",
    format = "rust-markdown",
    source = "crates/spenso/src/network/parsing/mod.rs",
    source_id = "spenso::network::parsing::ParseSettings"
)]
fn spenso_parse_settings() {}

#[alphal00p_docs::trait_item(
    id = "NetworkParse",
    title = "Symbolic network parsing",
    summary = "Converts Symbolica atoms and atom views into parametric tensor networks under explicit parsing settings.",
    format = "rust-markdown",
    source = "crates/spenso/src/network/parsing/mod.rs",
    source_id = "spenso::network::parsing::NetworkParse"
)]
fn spenso_network_parse() {}

#[alphal00p_docs::ty(
    id = "SymbolicParallelism",
    title = "Symbolic parallelism policy",
    summary = "Chooses automatic license-aware, serial, or forced-parallel Rayon behavior for operations over Symbolica atoms.",
    format = "rust-markdown",
    source = "crates/spenso/src/symbolic_parallelism.rs",
    source_id = "spenso::symbolic_parallelism::SymbolicParallelism"
)]
fn spenso_symbolic_parallelism() {}

#[alphal00p_docs::ty(
    id = "ParamTensor",
    title = "Parametric tensor",
    summary = "Stores tensor components as parameter-dependent symbolic coefficients.",
    format = "rust-markdown",
    source = "crates/spenso/src/tensors/parametric.rs",
    source_id = "spenso::tensors::parametric::ParamTensor"
)]
fn spenso_param_tensor() {}

#[alphal00p_docs::macro_item(
    id = "SimpleRepresentation",
    title = "SimpleRepresentation derive",
    summary = "Derives representation naming and duality boilerplate from the representation helper attribute.",
    format = "rust-markdown",
    source = "crates/spenso-macros/src/lib.rs",
    source_id = "spenso_macros::SimpleRepresentation"
)]
fn spenso_macros_simple_representation() {}

#[alphal00p_docs::func(
    id = "hep_lib",
    title = "HEP tensor library",
    summary = "Constructs the standard high-energy-physics tensor library for a chosen scalar data type.",
    format = "rust-markdown",
    source = "crates/spenso-hep-lib/src/lib.rs",
    source_id = "spenso_hep_lib::hep_lib"
)]
fn spenso_hep_lib_hep_lib() {}

#[alphal00p_docs::func(
    id = "gamma_data_dirac",
    title = "Dirac gamma data",
    summary = "Creates sparse Dirac-basis gamma-matrix data for the supplied tensor structure.",
    format = "rust-markdown",
    source = "crates/spenso-hep-lib/src/lib.rs",
    source_id = "spenso_hep_lib::gamma_data_dirac"
)]
fn spenso_hep_lib_gamma_data_dirac() {}

#[alphal00p_docs::func(
    id = "su3_generator_data",
    title = "SU(3) generator data",
    summary = "Creates the standard sparse SU(3) fundamental-generator tensor.",
    format = "rust-markdown",
    source = "crates/spenso-hep-lib/src/lib.rs",
    source_id = "spenso_hep_lib::su3_generator_data"
)]
fn spenso_hep_lib_su3_generator_data() {}

#[alphal00p_docs::trait_item(
    id = "IndexTooling",
    title = "Index tooling",
    summary = "Adds canonicalization, wrapping, conjugation, and dangling-index inspection to Symbolica atoms.",
    format = "rust-markdown",
    source = "crates/idenso/src/lib.rs",
    source_id = "idenso::IndexTooling"
)]
fn idenso_index_tooling() {}

#[alphal00p_docs::ty(
    id = "CookSettings",
    title = "Index cooking settings",
    summary = "Controls reversible, flattened, and filtered encodings of tensor indices.",
    format = "rust-markdown",
    source = "crates/idenso/src/cook.rs",
    source_id = "idenso::cook::CookSettings"
)]
fn idenso_cook_settings() {}

#[alphal00p_docs::trait_item(
    id = "Cookable",
    title = "Cookable expressions",
    summary = "Applies or reverses index-cooking transformations on symbolic expressions.",
    format = "rust-markdown",
    source = "crates/idenso/src/cook.rs",
    source_id = "idenso::cook::Cookable"
)]
fn idenso_cookable() {}

#[alphal00p_docs::trait_item(
    id = "SelectiveExpand",
    title = "Selective expansion",
    summary = "Expands only tensor families selected by their representation.",
    format = "rust-markdown",
    source = "crates/idenso/src/selective_expand.rs",
    source_id = "idenso::selective_expand::SelectiveExpand"
)]
fn idenso_selective_expand() {}

#[alphal00p_docs::macro_item(
    id = "bis",
    name = "bis!",
    title = "Bispinor representation macro",
    summary = "Builds stripped or indexed bispinor representation atoms in Spenso syntax.",
    format = "rust-markdown",
    source = "crates/idenso/src/lib.rs",
    source_id = "idenso::bis"
)]
fn idenso_bis() {}

#[alphal00p_docs::macro_item(
    id = "cof",
    name = "cof!",
    title = "Color-fundamental macro",
    summary = "Builds stripped or indexed color-fundamental representation atoms.",
    format = "rust-markdown",
    source = "crates/idenso/src/lib.rs",
    source_id = "idenso::cof"
)]
fn idenso_cof() {}

#[alphal00p_docs::macro_item(
    id = "coad",
    name = "coad!",
    title = "Color-adjoint representation macro",
    summary = "Builds stripped or indexed color-adjoint representation atoms.",
    format = "rust-markdown",
    source = "crates/idenso/src/lib.rs",
    source_id = "idenso::coad"
)]
fn idenso_coad() {}

#[alphal00p_docs::macro_item(
    id = "gamma",
    name = "gamma!",
    title = "Dirac gamma macro",
    summary = "Builds a Dirac gamma chain factor or an explicitly indexed gamma tensor.",
    format = "rust-markdown",
    source = "crates/idenso/src/dirac/macros.rs",
    source_id = "idenso::gamma"
)]
fn idenso_gamma() {}

#[alphal00p_docs::macro_item(
    id = "gamma0",
    name = "gamma0!",
    title = "Gamma-zero macro",
    summary = "Builds a gamma-zero chain factor or an explicitly indexed gamma-zero tensor.",
    format = "rust-markdown",
    source = "crates/idenso/src/dirac/macros.rs",
    source_id = "idenso::gamma0"
)]
fn idenso_gamma0() {}

#[alphal00p_docs::macro_item(
    id = "gamma5",
    name = "gamma5!",
    title = "Gamma-five macro",
    summary = "Builds a gamma-five chain factor or an explicitly indexed gamma-five tensor.",
    format = "rust-markdown",
    source = "crates/idenso/src/dirac/macros.rs",
    source_id = "idenso::gamma5"
)]
fn idenso_gamma5() {}

#[alphal00p_docs::macro_item(
    id = "epsilon",
    name = "epsilon!",
    title = "Levi-Civita tensor macro",
    summary = "Builds an arbitrary-rank antisymmetric epsilon tensor from atom-like arguments.",
    format = "rust-markdown",
    source = "crates/idenso/src/epsilon.rs",
    source_id = "idenso::epsilon"
)]
fn idenso_epsilon() {}

#[alphal00p_docs::macro_item(
    id = "color_t",
    name = "color_t!",
    title = "Color-generator macro",
    summary = "Builds a color-generator chain factor or an explicitly indexed generator tensor.",
    format = "rust-markdown",
    source = "crates/idenso/src/color/macros.rs",
    source_id = "idenso::color_t"
)]
fn idenso_color_t() {}

#[alphal00p_docs::macro_item(
    id = "color_f",
    name = "color_f!",
    title = "Color structure-constant macro",
    summary = "Builds an explicitly indexed adjoint color structure-constant tensor.",
    format = "rust-markdown",
    source = "crates/idenso/src/color/macros.rs",
    source_id = "idenso::color_f"
)]
fn idenso_color_f() {}

#[alphal00p_docs::func(
    id = "representations::initialize",
    title = "Initialize representations",
    summary = "Initializes Idenso's standard representations and tensor symbols.",
    format = "rust-markdown",
    source = "crates/idenso/src/representations.rs",
    source_id = "idenso::representations::initialize"
)]
fn idenso_initialize() {}

#[alphal00p_docs::macro_item(
    id = "vakint_parse",
    name = "vakint_parse!",
    title = "Vakint expression parser",
    summary = "Parses a Symbolica expression using Vakint's namespace by default or a caller-supplied namespace.",
    format = "rust-markdown",
    source = "crates/vakint/src/utils.rs",
    source_id = "vakint::vakint_parse"
)]
fn vakint_parse_macro() {}

#[alphal00p_docs::macro_item(
    id = "vakint_symbol",
    name = "vakint_symbol!",
    title = "Vakint symbol parser",
    summary = "Resolves a symbol name in Vakint's namespace while preserving an explicit namespace.",
    format = "rust-markdown",
    source = "crates/vakint/src/utils.rs",
    source_id = "vakint::vakint_symbol"
)]
fn vakint_symbol_macro() {}

#[alphal00p_docs::ty(
    id = "Vakint",
    title = "Vakint engine",
    summary = "Owns supported topology matching, canonicalization, reduction, and evaluation operations.",
    format = "rust-markdown",
    source = "crates/vakint/src/lib.rs",
    source_id = "vakint::Vakint"
)]
fn vakint_engine() {}

#[alphal00p_docs::ty(
    id = "VakintSettings",
    title = "Vakint settings",
    summary = "Controls symbols, precision, normalization, external tools, and backend evaluation order.",
    kind = "setting",
    format = "rust-markdown",
    source = "crates/vakint/src/lib.rs",
    source_id = "vakint::VakintSettings"
)]
fn vakint_settings() {}

#[alphal00p_docs::ty(
    id = "VakintExpression",
    title = "Vakint expression",
    summary = "Stores a sum of matched vacuum-integral terms with separate numerators.",
    format = "rust-markdown",
    source = "crates/vakint/src/lib.rs",
    source_id = "vakint::VakintExpression"
)]
fn vakint_expression() {}

#[alphal00p_docs::ty(
    id = "EvaluationOrder",
    title = "Evaluation order",
    summary = "Defines the ordered analytic and numerical backends available to an evaluation.",
    format = "rust-markdown",
    source = "crates/vakint/src/lib.rs",
    source_id = "vakint::EvaluationOrder"
)]
fn vakint_evaluation_order() {}

#[alphal00p_docs::func(
    id = "VakintExpression::canonicalize",
    title = "Canonicalize an expression",
    summary = "Matches and rewrites each integral into Vakint's canonical topology representation.",
    format = "rust-markdown",
    source = "crates/vakint/src/lib.rs",
    source_id = "vakint::VakintExpression::canonicalize"
)]
fn vakint_expression_canonicalize() {}

#[alphal00p_docs::func(
    id = "VakintExpression::tensor_reduce",
    title = "Tensor-reduce an expression",
    summary = "Reduces tensor numerators to scalar-integral combinations.",
    format = "rust-markdown",
    source = "crates/vakint/src/lib.rs",
    source_id = "vakint::VakintExpression::tensor_reduce"
)]
fn vakint_expression_tensor_reduce() {}

#[alphal00p_docs::func(
    id = "VakintExpression::evaluate_integral",
    title = "Evaluate an expression",
    summary = "Evaluates canonical scalar integrals using the configured ordered backend policy.",
    format = "rust-markdown",
    source = "crates/vakint/src/lib.rs",
    source_id = "vakint::VakintExpression::evaluate_integral"
)]
fn vakint_expression_evaluate_integral() {}

#[alphal00p_docs::ty(
    id = "NumericalEvaluationResult",
    title = "Numerical Laurent result",
    summary = "Stores numerical Laurent coefficients by power of the dimensional regulator.",
    format = "rust-markdown",
    source = "crates/vakint/src/lib.rs",
    source_id = "vakint::NumericalEvaluationResult"
)]
fn vakint_numerical_evaluation_result() {}

pub(super) fn for_component(component: &str) -> Option<Vec<DocItem>> {
    Some(match component {
        "gammalooprs" => vec![
            __alphal00p_docs_trait_gammalooprs_gamma_loop_context(),
            __alphal00p_docs_trait_gammalooprs_has_model(),
            __alphal00p_docs_ty_gammalooprs_gamma_loop_context_container(),
            __alphal00p_docs_trait_gammalooprs_has_integrand(),
            __alphal00p_docs_ty_gammalooprs_integrand(),
            __alphal00p_docs_ty_gammalooprs_global_settings(),
            __alphal00p_docs_ty_gammalooprs_runtime_settings(),
            __alphal00p_docs_ty_gammalooprs_single_sample_evaluation_result(),
            __alphal00p_docs_ty_gammalooprs_batch_sample_evaluation_result(),
            __alphal00p_docs_ty_gammalooprs_precise_single_sample_evaluation_result(),
            __alphal00p_docs_ty_gammalooprs_precise_batch_sample_evaluation_result(),
            __alphal00p_docs_ty_gammalooprs_histogram_snapshot(),
            __alphal00p_docs_func_gammalooprs_histogram_snapshot_merge(),
            __alphal00p_docs_func_gammalooprs_histogram_snapshot_rebin(),
            __alphal00p_docs_ty_gammalooprs_histogram_accumulator_state(),
            __alphal00p_docs_func_gammalooprs_histogram_accumulator_state_snapshot(),
            __alphal00p_docs_func_gammalooprs_set_interrupt_handler(),
            __alphal00p_docs_func_gammalooprs_request_interrupt(),
            __alphal00p_docs_func_gammalooprs_is_interrupt_requested(),
            __alphal00p_docs_func_gammalooprs_clear_interrupt_request(),
        ],
        "gammaloop-api" => vec![
            __alphal00p_docs_ty_gammaloop_api_state_load_option(),
            __alphal00p_docs_ty_gammaloop_api_loaded_state(),
            __alphal00p_docs_func_gammaloop_api_state_load_option_load(),
            __alphal00p_docs_func_gammaloop_api_loaded_state_cli_session(),
            __alphal00p_docs_ty_gammaloop_api_cli_session(),
            __alphal00p_docs_func_gammaloop_api_cli_session_execute_command(),
            __alphal00p_docs_ty_gammaloop_api_command_history(),
            __alphal00p_docs_func_gammaloop_api_command_history_from_raw_string(),
            __alphal00p_docs_ty_gammaloop_api_run_history(),
            __alphal00p_docs_func_gammaloop_api_run_history_load(),
            __alphal00p_docs_ty_gammaloop_api_state(),
            __alphal00p_docs_func_gammaloop_api_state_get_integrand_info(),
            __alphal00p_docs_ty_gammaloop_api_integrand_info(),
            __alphal00p_docs_ty_gammaloop_api_import_model(),
            __alphal00p_docs_func_gammaloop_api_import_model_run(),
            __alphal00p_docs_ty_gammaloop_api_generate(),
            __alphal00p_docs_func_gammaloop_api_generate_run(),
            __alphal00p_docs_ty_gammaloop_api_integrate(),
            __alphal00p_docs_func_gammaloop_api_integrate_run(),
            __alphal00p_docs_ty_gammaloop_api_integration_output(),
            __alphal00p_docs_ty_gammaloop_api_evaluate_samples_request(),
            __alphal00p_docs_func_gammaloop_api_evaluate_sample(),
            __alphal00p_docs_func_gammaloop_api_evaluate_samples(),
            __alphal00p_docs_ty_gammaloop_api_evaluate_samples_precise_request(),
            __alphal00p_docs_func_gammaloop_api_evaluate_sample_precise(),
            __alphal00p_docs_func_gammaloop_api_evaluate_samples_precise(),
            __alphal00p_docs_ty_gammaloop_api_state_folder_kind(),
            __alphal00p_docs_func_gammaloop_api_classify_state_folder(),
            __alphal00p_docs_ty_gammaloop_api_cli_settings(),
            __alphal00p_docs_ty_gammaloop_api_command(),
        ],
        "linnet" => vec![
            __alphal00p_docs_ty_linnet_hedge_graph(),
            __alphal00p_docs_ty_linnet_hedge_graph_builder(),
            __alphal00p_docs_ty_linnet_su_bit_graph(),
            __alphal00p_docs_ty_linnet_simple_traversal_tree(),
            __alphal00p_docs_ty_linnet_dot_graph(),
            __alphal00p_docs_ty_linnet_hedge_pair(),
            __alphal00p_docs_ty_linnet_flow(),
            __alphal00p_docs_ty_linnet_orientation(),
            __alphal00p_docs_trait_linnet_subset_like(),
            __alphal00p_docs_trait_linnet_subset_ops(),
            __alphal00p_docs_trait_linnet_inclusion(),
            __alphal00p_docs_trait_linnet_base_subgraph(),
            __alphal00p_docs_trait_linnet_subgraph_like(),
            __alphal00p_docs_trait_linnet_node_storage(),
            __alphal00p_docs_trait_linnet_node_storage_ops(),
        ],
        "spenso" => vec![
            __alphal00p_docs_trait_spenso_tensor_structure(),
            __alphal00p_docs_ty_spenso_dense_tensor(),
            __alphal00p_docs_ty_spenso_sparse_tensor(),
            __alphal00p_docs_trait_spenso_contract(),
            __alphal00p_docs_ty_spenso_network(),
            __alphal00p_docs_trait_spenso_contraction_strategy(),
            __alphal00p_docs_trait_spenso_execution_strategy(),
            __alphal00p_docs_ty_spenso_execution_mode(),
            __alphal00p_docs_ty_spenso_execution_result(),
            __alphal00p_docs_ty_spenso_tensor_network_error(),
            __alphal00p_docs_ty_spenso_parse_settings(),
            __alphal00p_docs_trait_spenso_network_parse(),
            __alphal00p_docs_ty_spenso_symbolic_parallelism(),
            __alphal00p_docs_ty_spenso_param_tensor(),
        ],
        "spenso-macros" => {
            vec![__alphal00p_docs_macro_spenso_macros_simple_representation()]
        }
        "spenso-hep-lib" => vec![
            __alphal00p_docs_func_spenso_hep_lib_hep_lib(),
            __alphal00p_docs_func_spenso_hep_lib_gamma_data_dirac(),
            __alphal00p_docs_func_spenso_hep_lib_su3_generator_data(),
        ],
        "idenso" => vec![
            __alphal00p_docs_trait_idenso_index_tooling(),
            __alphal00p_docs_ty_idenso_cook_settings(),
            __alphal00p_docs_trait_idenso_cookable(),
            __alphal00p_docs_trait_idenso_selective_expand(),
            __alphal00p_docs_macro_idenso_bis(),
            __alphal00p_docs_macro_idenso_cof(),
            __alphal00p_docs_macro_idenso_coad(),
            __alphal00p_docs_macro_idenso_gamma(),
            __alphal00p_docs_macro_idenso_gamma0(),
            __alphal00p_docs_macro_idenso_gamma5(),
            __alphal00p_docs_macro_idenso_epsilon(),
            __alphal00p_docs_macro_idenso_color_t(),
            __alphal00p_docs_macro_idenso_color_f(),
            __alphal00p_docs_func_idenso_initialize(),
        ],
        "vakint" => vec![
            __alphal00p_docs_macro_vakint_parse_macro(),
            __alphal00p_docs_macro_vakint_symbol_macro(),
            __alphal00p_docs_ty_vakint_engine(),
            __alphal00p_docs_ty_vakint_settings(),
            __alphal00p_docs_ty_vakint_expression(),
            __alphal00p_docs_ty_vakint_evaluation_order(),
            __alphal00p_docs_func_vakint_expression_canonicalize(),
            __alphal00p_docs_func_vakint_expression_tensor_reduce(),
            __alphal00p_docs_func_vakint_expression_evaluate_integral(),
            __alphal00p_docs_ty_vakint_numerical_evaluation_result(),
        ],
        _ => return None,
    })
}
