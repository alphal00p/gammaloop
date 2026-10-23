use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use alphal00p_persistence::OutputPath;

// use bincode::{Decode, Encode};
use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;
use eyre::eyre;
use rayon::ThreadPool;
use schemars::JsonSchema;
use std::result::Result as StdResult;
use symbolica::evaluate::OptimizationSettings;
use thiserror::Error;

use crate::{
    GammaLoopContext, GammaLoopContextContainer,
    settings::{GlobalSettings, runtime::LockedRuntimeSettings},
    utils::serde_utils::{IsDefault, is_false, is_true, is_usize, show_defaults_helper},
    uv::export::UVForestExportSettings,
};
use serde::{Deserialize, Serialize};

use crate::model::Model;
use feynkit_generator::GenerationType;

#[cfg(test)]
mod output_tests;

/// Directory containing process payloads below a saved generation or export root.
///
/// ```
/// use gammalooprs::processes::PROCESSES_DIRECTORY;
/// use std::path::Path;
/// assert_eq!(
///     Path::new("export").join(PROCESSES_DIRECTORY),
///     Path::new("export/processes"),
/// );
/// ```
pub const PROCESSES_DIRECTORY: &str = "processes";

pub(crate) const PROCESS_DEFINITION_FILE: &str = "def.bin";
pub(crate) const AMPLITUDE_PAYLOAD_FILE: &str = "amp.bin";
pub(crate) const CROSS_SECTION_PAYLOAD_FILE: &str = "cs.bin";
/// Directory containing a generated integrand's binary payload and settings.
/// ```
/// use gammalooprs::processes::INTEGRAND_DIRECTORY;
/// use std::path::Path;
/// assert_eq!(Path::new("child").join(INTEGRAND_DIRECTORY), Path::new("child/integrand"));
/// ```
pub const INTEGRAND_DIRECTORY: &str = "integrand";
pub(crate) const INTEGRAND_BINARY_FILE: &str = "integrand/integrand.bin";
pub(crate) const INTEGRAND_SETTINGS_FILE: &str = "integrand/settings.toml";

mod snapshot_error;
pub use snapshot_error::{ProcessPayloadError, ProcessPayloadFailure, ProcessSnapshotError};

mod generation_report;
pub use generation_report::{
    EvaluatorBuildTimings, GeneratedGraphKey, GeneratedGraphReport, GraphGenerationStats,
    NamedGraphGenerationReport, merge_generated_graph_reports,
};
mod generation_progress;
pub use generation_progress::{
    GenerationProcessKind, GenerationProgressMode, GenerationProgressModeGuard,
    GenerationProgressObserver, GenerationProgressObserverGuard, GenerationProgressPhase,
    begin_phase, cut_finished, detailed_progress_enabled, detailed_progress_message,
    enter_detailed_progress_span, enter_progress_context,
};
mod selection;
pub use selection::{
    CycleSignature, GraphGroupSelectionMode, GraphGroupSelectionPlan, GraphGroupSelectionReport,
    GraphGroupSelectionSpec, GraphSelectionSignatureInventory, ParticleSignature,
    RaisedCutSignatureInventory, RaisedPropagatorScope, RaisedPropagatorSignature,
    SelectionPolarity, VertexSignature,
};
pub(crate) use selection::{GraphCutSelectionSubject, GraphSelectionSubject};

#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[derive(
    Debug, Clone, Copy, Default, Serialize, Deserialize, Encode, Decode, PartialEq, Eq, JsonSchema,
)]
#[serde(rename_all = "snake_case")]
pub enum TensorNetworkContractionOrder {
    SparseAtomAware,
    AtomAware,
    ResultRankOnly,
    EntryAware,
    #[default]
    IntermediateCost,
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[derive(Debug, Clone, Copy, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
pub enum ExecutionMode {
    Sequential,
    Parallel,
    SequentialRef,
    SequentialExtract,
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[derive(Debug, Clone, Copy, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
pub enum ContractionMode {
    SmallestDegree,
    MinResultRank,
}
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[derive(Debug, Clone, Copy, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
pub struct EvaluatorSettings {
    /// Perform symbolic tensor and scalar algebra before lowering to an evaluator.
    #[serde(default, skip_serializing_if = "is_false")]
    pub do_algebra: bool,
    /// Re-optimize orientation expressions iteratively while building the evaluator.
    #[serde(
        default = "evaluator_default_iterative_orientation_optimization",
        skip_serializing_if = "is_true"
    )]
    pub iterative_orientation_optimization: bool,
    /// Build one evaluator for the sum of all contributions instead of separate evaluators.
    #[serde(default, skip_serializing_if = "is_false")]
    pub summed: bool,
    /// Sum all orientations, retaining each selected orientation as a separate
    /// Symbolica function-map alias instead of merging them into one expression.
    #[serde(default, skip_serializing_if = "is_false")]
    pub summed_function_map: bool,
    /// Compile the lowered evaluator using the configured generation compilation backend.
    #[serde(default, skip_serializing_if = "is_false")]
    pub compile: bool,
    /// Retain the intermediate Symbolica atom alongside the executable evaluator.
    #[serde(default, skip_serializing_if = "is_false")]
    pub store_atom: bool,
    /// Apply function-map replacements before numerical lowering.
    #[serde(default, skip_serializing_if = "is_false")]
    pub do_fn_map_replacements: bool,
    /// Inline generated numerator functions instead of retaining shared evaluator bodies.
    #[serde(
        default = "evaluator_default_inline_numerator_functions",
        skip_serializing_if = "is_true"
    )]
    pub inline_numerator_functions: bool,
    /// Translate compatible symbolic expressions directly into Spenso networks.
    #[serde(
        default = "evaluator_default_direct_translation",
        skip_serializing_if = "is_true"
    )]
    pub direct_translation: bool,
    /// Number of Horner-scheme optimization passes applied to scalar expressions.
    #[serde(
        default = "evaluator_default_horner_iterations",
        skip_serializing_if = "is_usize::<1>"
    )]
    pub horner_iterations: usize,
    /// Worker threads used while constructing one evaluator.
    #[serde(
        default = "evaluator_default_n_cores",
        skip_serializing_if = "is_usize::<1>"
    )]
    pub n_cores: usize,
    /// Optional number of common-pair-elimination optimization passes.
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub cpe_iterations: Option<usize>,
    /// Optimization stage at which evaluator construction aborts for diagnostic isolation.
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub abort_level: usize,

    /// Maximum variable count for attempting a Horner-scheme optimization.
    #[serde(
        default = "evaluator_default_max_horner_scheme_variables",
        skip_serializing_if = "is_usize::<500>"
    )]
    pub max_horner_scheme_variables: usize,

    /// Maximum cached candidate pairs during common-pair elimination.
    #[serde(
        default = "evaluator_default_max_common_pair_cache_entries",
        skip_serializing_if = "is_usize::<1000000>"
    )]
    pub max_common_pair_cache_entries: usize,

    /// Maximum expression-tree distance considered during common-pair elimination.
    #[serde(
        default = "evaluator_default_max_common_pair_distance",
        skip_serializing_if = "is_usize::<1000>"
    )]
    pub max_common_pair_distance: usize,
    /// Heuristic used to order tensor-network contractions.
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub tensor_network_contraction_order: TensorNetworkContractionOrder,
    /// Emit detailed evaluator-construction diagnostics.
    #[serde(default, skip_serializing_if = "is_false")]
    pub verbose: bool,

    /// Spenso execution and contraction modes used by the lowered tensor network.
    #[serde(
        default = "evaluator_default_spenso_execution_mode",
        skip_serializing_if = "is_default_spenso_execution_mode"
    )]
    pub spenso_execution_mode: (ExecutionMode, ContractionMode),
}

const fn evaluator_default_iterative_orientation_optimization() -> bool {
    true
}

const fn evaluator_default_direct_translation() -> bool {
    true
}

const fn evaluator_default_inline_numerator_functions() -> bool {
    true
}

const fn evaluator_default_horner_iterations() -> usize {
    1
}

const fn evaluator_default_n_cores() -> usize {
    1
}

const fn evaluator_default_max_horner_scheme_variables() -> usize {
    500
}

const fn evaluator_default_max_common_pair_cache_entries() -> usize {
    1_000_000
}

const fn evaluator_default_max_common_pair_distance() -> usize {
    1000
}

const fn evaluator_default_spenso_execution_mode() -> (ExecutionMode, ContractionMode) {
    (ExecutionMode::Sequential, ContractionMode::MinResultRank)
}

fn is_default_spenso_execution_mode(mode: &(ExecutionMode, ContractionMode)) -> bool {
    show_defaults_helper(mode == &evaluator_default_spenso_execution_mode())
}

impl Default for EvaluatorSettings {
    fn default() -> Self {
        Self {
            iterative_orientation_optimization:
                evaluator_default_iterative_orientation_optimization(),
            summed: false,
            do_algebra: false,
            summed_function_map: false,
            direct_translation: evaluator_default_direct_translation(),
            compile: false,
            do_fn_map_replacements: false,
            inline_numerator_functions: evaluator_default_inline_numerator_functions(),
            store_atom: false,
            horner_iterations: evaluator_default_horner_iterations(),
            n_cores: evaluator_default_n_cores(),
            cpe_iterations: None,
            abort_level: 0,
            max_horner_scheme_variables: evaluator_default_max_horner_scheme_variables(),
            max_common_pair_cache_entries: evaluator_default_max_common_pair_cache_entries(),
            max_common_pair_distance: evaluator_default_max_common_pair_distance(),
            tensor_network_contraction_order: TensorNetworkContractionOrder::default(),
            verbose: false,
            spenso_execution_mode: evaluator_default_spenso_execution_mode(),
        }
    }
}

impl EvaluatorSettings {
    pub fn optimization_settings(&self) -> OptimizationSettings {
        OptimizationSettings::new()
            .horner_iterations(self.horner_iterations)
            .cores(self.n_cores)
            .cpe_iterations(self.cpe_iterations)
            .abort_check(Some(Box::new(
                crate::is_interrupt_requested as fn() -> bool,
            )))
            .abort_level(self.abort_level)
            .max_horner_scheme_variables(self.max_horner_scheme_variables)
            .max_common_pair_cache_entries(self.max_common_pair_cache_entries)
            .max_common_pair_distance(self.max_common_pair_distance)
            .verbose(self.verbose)
            .direct_translation(self.direct_translation)
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct ProcessList {
    // pub amplitude_process: Vec<Process<Amplitude>>
    // pub cross_section_process: Vec<Process<CrossSection>>
    pub processes: Vec<Process>,
}

/// A process-list operation requested an index outside its current membership.
///
/// See [`ProcessList::remove_process`] for an example.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Error)]
#[error("Process ID {index} is out of range for a list of {len} processes")]
pub struct ProcessIndexError {
    /// Requested zero-based process index; valid indices are strictly less than `len`.
    pub index: usize,
    /// Process count when the operation rejected the index.
    pub len: usize,
}

/// Restricts which saved processes and integrands are loaded into a process list.
///
/// An unqualified integrand selector applies to every process containing an
/// integrand with that name, while a qualified selector applies only to its
/// named process. Process selectors take precedence and load every integrand
/// in the selected process.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ProcessLoadSelection {
    pub process_names: BTreeSet<String>,
    pub integrand_selectors: Vec<ProcessLoadIntegrandSelector>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ProcessLoadIntegrandSelector {
    Name(String),
    Qualified {
        process_name: String,
        integrand_name: String,
    },
}

impl ProcessLoadSelection {
    pub fn is_unrestricted(&self) -> bool {
        self.process_names.is_empty() && self.integrand_selectors.is_empty()
    }

    pub fn selected_integrand_names<'a, I>(
        &self,
        process_name: &str,
        available_names: I,
    ) -> BTreeSet<String>
    where
        I: IntoIterator<Item = &'a str>,
    {
        let available_names = available_names.into_iter();
        if self.process_names.contains(process_name) {
            return available_names.map(str::to_string).collect();
        }

        available_names
            .filter(|integrand_name| {
                self.integrand_selectors
                    .iter()
                    .any(|selector| match selector {
                        ProcessLoadIntegrandSelector::Name(name) => name == *integrand_name,
                        ProcessLoadIntegrandSelector::Qualified {
                            process_name: selected_process,
                            integrand_name: selected_integrand,
                        } => {
                            selected_process == process_name
                                && selected_integrand == *integrand_name
                        }
                    })
            })
            .map(str::to_string)
            .collect()
    }
}

#[cfg(test)]
mod process_load_selection_tests {
    use super::{ProcessLoadIntegrandSelector, ProcessLoadSelection};
    use std::collections::BTreeSet;

    #[test]
    fn process_and_integrand_selectors_are_unioned() {
        let selection = ProcessLoadSelection {
            process_names: BTreeSet::from(["all".to_string()]),
            integrand_selectors: vec![
                ProcessLoadIntegrandSelector::Name("shared".to_string()),
                ProcessLoadIntegrandSelector::Qualified {
                    process_name: "one".to_string(),
                    integrand_name: "specific".to_string(),
                },
            ],
        };

        assert_eq!(
            selection.selected_integrand_names("all", ["a", "b"]),
            BTreeSet::from(["a".to_string(), "b".to_string()])
        );
        assert_eq!(
            selection.selected_integrand_names("one", ["shared", "specific", "other"]),
            BTreeSet::from(["shared".to_string(), "specific".to_string()])
        );
        assert_eq!(
            selection.selected_integrand_names("two", ["shared", "specific"]),
            BTreeSet::from(["shared".to_string()])
        );
    }
}

/// Exact, ordered membership of a saved process list; position is the process ID.
/// Produced by [`ProcessList::write_snapshot`] and consumed by [`ProcessList::load`].
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ProcessSnapshotIndex {
    pub processes: Vec<ProcessSnapshotEntry>,
}

/// One process in [`ProcessSnapshotIndex`]. See [`ProcessList::write_snapshot`].
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ProcessSnapshotEntry {
    pub kind: GenerationType,
    pub folder_name: String,
    pub children: Vec<IntegrandSnapshotEntry>,
    pub settings_history: bool,
}

/// One child and its required generated-payload status in a process snapshot.
/// See [`ProcessList::write_snapshot`].
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct IntegrandSnapshotEntry {
    pub name: String,
    pub payload: IntegrandPayload,
}

/// Whether a snapshot child requires a generated integrand payload.
///
/// Metadata-only children contain graphs but no executable integrand. Generated
/// children require both their binary integrand and its pinned runtime settings.
///
/// ```
/// use gammalooprs::processes::{IntegrandPayload, IntegrandSnapshotEntry};
/// let child = IntegrandSnapshotEntry {
///     name: "LO".into(),
///     payload: IntegrandPayload::Generated,
/// };
/// assert_eq!(child.payload, IntegrandPayload::Generated);
/// ```
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Hash)]
#[serde(rename_all = "snake_case")]
pub enum IntegrandPayload {
    /// Only the child's graph metadata is present.
    MetadataOnly,
    /// Binary integrand data and pinned runtime settings are required.
    Generated,
}

impl ProcessSnapshotEntry {
    /// Storage category for this process. See [`ProcessList::write_snapshot`].
    pub fn kind_directory(&self) -> &'static str {
        ProcessCollection::directory_for_kind(self.kind)
    }

    /// Required paths relative to a complete state generation.
    /// See [`ProcessList::write_snapshot`].
    pub fn required_files(&self) -> Vec<PathBuf> {
        let root = Path::new(PROCESSES_DIRECTORY)
            .join(self.kind_directory())
            .join(&self.folder_name);
        let mut files = vec![root.join(PROCESS_DEFINITION_FILE)];
        if self.settings_history {
            files.push(root.join("settings_history.toml"));
        }
        for child in &self.children {
            let child_root = root.join(&child.name);
            files.push(child_root.join(match self.kind {
                GenerationType::Amplitude => AMPLITUDE_PAYLOAD_FILE,
                GenerationType::CrossSection => CROSS_SECTION_PAYLOAD_FILE,
            }));
            match child.payload {
                IntegrandPayload::Generated => {
                    files.push(child_root.join(INTEGRAND_BINARY_FILE));
                    files.push(child_root.join(INTEGRAND_SETTINGS_FILE));
                }
                IntegrandPayload::MetadataOnly => {}
            }
        }
        files
    }
}

impl ProcessSnapshotIndex {
    /// Validate unique process names and strictly ordered child names.
    ///
    /// # Errors
    /// Returns a typed membership error or the original path validation error.
    ///
    /// ```
    /// use alphal00p_persistence::OutputPathError;
    /// use feynkit_generator::GenerationType;
    /// use gammalooprs::processes::{
    ///     ProcessSnapshotEntry, ProcessSnapshotError, ProcessSnapshotIndex,
    /// };
    /// let index = ProcessSnapshotIndex {
    ///     processes: vec![ProcessSnapshotEntry {
    ///         kind: GenerationType::Amplitude,
    ///         folder_name: "../outside".into(),
    ///         children: vec![],
    ///         settings_history: false,
    ///     }],
    /// };
    /// assert!(matches!(
    ///     index.validate(),
    ///     Err(ProcessSnapshotError::Output(OutputPathError::InvalidComponent(_)))
    /// ));
    /// ```
    pub fn validate(&self) -> StdResult<(), ProcessSnapshotError> {
        let mut names = BTreeSet::new();
        for entry in &self.processes {
            OutputPath::validate_component(&entry.folder_name)?;
            if !names.insert(&entry.folder_name) {
                return Err(ProcessSnapshotError::DuplicateProcess {
                    process: entry.folder_name.clone(),
                });
            }
            let mut previous: Option<&str> = None;
            for child in &entry.children {
                OutputPath::validate_component(&child.name)?;
                if let Some(previous) = previous.filter(|name| *name >= child.name.as_str()) {
                    return Err(ProcessSnapshotError::InvalidChildOrder {
                        process: entry.folder_name.clone(),
                        previous: previous.to_owned(),
                        child: child.name.clone(),
                    });
                }
                previous = Some(&child.name);
            }
        }
        Ok(())
    }

    fn validate_selection(
        &self,
        selection: &ProcessLoadSelection,
    ) -> StdResult<(), ProcessSnapshotError> {
        let available: BTreeMap<_, BTreeSet<_>> = self
            .processes
            .iter()
            .map(|entry| {
                (
                    entry.folder_name.as_str(),
                    entry
                        .children
                        .iter()
                        .map(|child| child.name.as_str())
                        .collect(),
                )
            })
            .collect();
        for name in &selection.process_names {
            if !available.contains_key(name.as_str()) {
                return Err(ProcessSnapshotError::UnknownProcess {
                    process: name.clone(),
                    available: available.keys().map(|name| (*name).to_owned()).collect(),
                });
            }
        }
        for selector in &selection.integrand_selectors {
            let found = match selector {
                ProcessLoadIntegrandSelector::Name(name) => available
                    .values()
                    .any(|names| names.contains(name.as_str())),
                ProcessLoadIntegrandSelector::Qualified {
                    process_name,
                    integrand_name,
                } => available
                    .get(process_name.as_str())
                    .is_some_and(|names| names.contains(integrand_name.as_str())),
            };
            if !found {
                return Err(ProcessSnapshotError::UnknownIntegrand {
                    selector: selector.clone(),
                });
            }
        }
        Ok(())
    }
}

/// Controls which graph details and algebraic transformations are included in DOT exports.
#[derive(Serialize, Deserialize, Debug, Clone)]
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(from_py_object, get_all, set_all)
)]
pub struct DotExportSettings {
    // pub root_folder: PathBuf,
    /// Write all diagrams of an integrand to one file during filesystem export.
    pub combine_diagrams: bool,
    /// Request UV counterterm graph output; currently retained for compatibility.
    pub with_uv: bool,
    /// Add the complete symbolic numerator as the graph-level `full_num` field.
    pub output_full_numerator: bool,
    /// Split cross-section graphs by distinct initial-state assignments.
    pub split_xs_by_initial_states: bool,
    /// Simplify gamma-matrix algebra in the full numerator before formatting it.
    pub do_gamma_algebra: bool,
    /// Simplify color algebra in the full numerator before formatting it.
    pub do_color_algebra: bool,
    #[serde(default, skip_serializing_if = "is_false")]
    /// Include vertex and edge fields whose values GammaLoop generated automatically.
    pub include_autogenerated_fields: bool,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(from_py_object, get_all, set_all)
)]
pub struct StandaloneExportSettings {
    #[serde(default)]
    pub mode: StandaloneExportMode,
    #[serde(default)]
    pub format: StandaloneDataFormat,
    #[serde(default)]
    pub precision: StandaloneNumericTarget,
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[derive(
    Serialize, Deserialize, Debug, Clone, Copy, PartialEq, Eq, Default, Encode, Decode, JsonSchema,
)]
#[serde(rename_all = "snake_case")]
pub enum StandaloneExportMode {
    #[default]
    Rust,
    Python,
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[derive(
    Serialize, Deserialize, Debug, Clone, Copy, PartialEq, Eq, Default, Encode, Decode, JsonSchema,
)]
#[serde(rename_all = "snake_case")]
pub enum StandaloneDataFormat {
    #[default]
    Binary,
    Json,
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[derive(
    Serialize, Deserialize, Debug, Clone, Copy, PartialEq, Eq, Default, Encode, Decode, JsonSchema,
)]
#[serde(rename_all = "snake_case")]
pub enum StandaloneNumericTarget {
    #[default]
    Double,
    Quad,
    Arb,
}

#[cfg(feature = "python_api")]
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(feature = "python_api", pyo3::pymethods)]
impl DotExportSettings {
    /// Construct DOT export settings with GammaLoop's defaults.
    #[new]
    fn new() -> Self {
        DotExportSettings::default()
    }
}

impl Default for DotExportSettings {
    fn default() -> Self {
        Self {
            // root_folder: PathBuf::from("."),
            split_xs_by_initial_states: false,
            output_full_numerator: false,
            do_gamma_algebra: false,
            do_color_algebra: true,
            combine_diagrams: false,
            with_uv: false,
            include_autogenerated_fields: false,
        }
    }
}

impl Default for StandaloneExportSettings {
    fn default() -> Self {
        Self {
            format: StandaloneDataFormat::Binary,
            mode: StandaloneExportMode::Rust,
            precision: StandaloneNumericTarget::Double,
        }
    }
}

impl Default for ProcessList {
    fn default() -> Self {
        Self::new()
    }
}

// Python Api
impl ProcessList {
    /// Generates a new empty process list
    pub fn new() -> Self {
        ProcessList { processes: vec![] }
    }

    pub fn warm_up_all_processes(&mut self, model: &Model) -> Result<()> {
        for process in &mut self.processes.iter_mut() {
            process.warm_up(model)?;
        }

        Ok(())
    }

    /// Load only indexed membership, ignoring stray directories.
    /// See [`Self::write_snapshot`] for the matching serializer.
    ///
    /// # Errors
    /// Distinguishes invalid membership and paths, unreadable or malformed
    /// definitions/settings, and nested payload errors with their original sources.
    pub fn load(
        root: &OutputPath,
        index: &ProcessSnapshotIndex,
        context: GammaLoopContextContainer<'_>,
    ) -> StdResult<Self, ProcessSnapshotError> {
        Self::load_with_selection(root, index, context, None)
    }

    /// Load the indexed union of selected processes and children, with dense IDs.
    /// See [`Self::write_snapshot`] for the matching serializer.
    ///
    /// # Errors
    /// Returns [`ProcessSnapshotError::UnknownProcess`] or
    /// [`ProcessSnapshotError::UnknownIntegrand`] for unknown selectors, in
    /// addition to the failures documented on [`Self::load`].
    pub fn load_with_selection(
        root: &OutputPath,
        index: &ProcessSnapshotIndex,
        context: GammaLoopContextContainer<'_>,
        selection: Option<&ProcessLoadSelection>,
    ) -> StdResult<Self, ProcessSnapshotError> {
        index.validate()?;
        let selection = selection.filter(|selection| !selection.is_unrestricted());
        if let Some(selection) = selection {
            index.validate_selection(selection)?;
        }
        let mut list = Self::new();
        for (process_id, entry) in index.processes.iter().enumerate() {
            let selected = selection.map(|selection| {
                selection.selected_integrand_names(
                    &entry.folder_name,
                    entry.children.iter().map(|child| child.name.as_str()),
                )
            });
            if selected.as_ref().is_some_and(BTreeSet::is_empty)
                && !selection
                    .is_some_and(|selection| selection.process_names.contains(&entry.folder_name))
            {
                continue;
            }
            let process =
                Process::load_snapshot(root, entry, process_id, context, selected.as_ref())?;
            list.processes.push(process);
        }
        if selection.is_some() {
            list.reindex_processes();
        }
        Ok(list)
    }

    /// Serialize into a fresh owned directory, never over an existing payload.
    ///
    /// This stages payloads only; the application publishes the complete session.
    /// The returned index is the only source of load membership.
    ///
    /// # Errors
    /// Validates membership before writing. Refuses existing payload directories;
    /// path and nested serializer errors retain their original sources.
    ///
    /// ```no_run
    /// # use gammalooprs::processes::ProcessList;
    /// # use std::error::Error;
    /// # use alphal00p_persistence::OutputPath;
    /// # fn main() -> Result<(), Box<dyn Error>> {
    /// let stage = OutputPath::new("fresh-stage".as_ref())?;
    /// stage.create_directory()?;
    /// let mut list = ProcessList::new();
    /// let index = list.write_snapshot(&stage)?;
    /// assert!(index.processes.is_empty());
    /// # Ok(()) }
    /// ```
    pub fn write_snapshot(
        &mut self,
        root: &OutputPath,
    ) -> StdResult<ProcessSnapshotIndex, ProcessSnapshotError> {
        let index = self.snapshot_index()?;
        let path = root.join(PROCESSES_DIRECTORY)?;
        path.create_directory()?;
        for kind in [GenerationType::Amplitude, GenerationType::CrossSection] {
            path.join(ProcessCollection::directory_for_kind(kind))?
                .create_directory()?;
        }
        for process in &mut self.processes {
            process.write_snapshot(&path)?;
        }
        Ok(index)
    }

    fn snapshot_index(&self) -> StdResult<ProcessSnapshotIndex, ProcessSnapshotError> {
        let mut entries = Vec::new();
        for (id, process) in self.processes.iter().enumerate() {
            if process.definition.process_id != id {
                return Err(ProcessSnapshotError::InvalidProcessId {
                    process: process.definition.folder_name.clone(),
                    expected: id,
                    actual: process.definition.process_id,
                });
            }
            entries.push(process.snapshot_entry()?);
        }
        let index = ProcessSnapshotIndex { processes: entries };
        index.validate()?;
        Ok(index)
    }

    /// Compile selected integrands, creating their artifact directories as needed.
    ///
    /// # Errors
    /// Returns path validation, directory creation, or evaluator compilation failures.
    ///
    /// ```no_run
    /// # use gammalooprs::processes::ProcessList;
    /// # use rayon::ThreadPoolBuilder;
    /// # use alphal00p_persistence::OutputPath;
    /// # fn main() -> eyre::Result<()> {
    /// let pool = ThreadPoolBuilder::new().build()?;
    /// let output = OutputPath::new("compiled".as_ref())?;
    /// ProcessList::new().compile(&output, None, None, &pool)?;
    /// # Ok(()) }
    /// ```
    pub fn compile(
        &mut self,
        folder: &OutputPath,
        process_id: Option<usize>,
        integrand_name: Option<String>,
        thread_pool: &ThreadPool,
    ) -> Result<Vec<GeneratedGraphReport>> {
        let path = folder.join(PROCESSES_DIRECTORY)?;
        path.ensure_directory()?;

        let mut reports = Vec::new();
        for p in self.processes.iter_mut() {
            if let Some(id) = process_id
                && p.definition.process_id != id
            {
                continue;
            }
            reports.extend(p.compile(&path, integrand_name.clone(), thread_pool)?);
        }

        Ok(reports)
    }

    pub fn activate_loaded_integrand_backends(
        &mut self,
        allow_symjit_fallback: bool,
    ) -> Result<()> {
        for process in &mut self.processes {
            process.activate_loaded_integrand_backends(allow_symjit_fallback)?;
        }
        Ok(())
    }

    pub fn get_integrand(
        &self,
        process_id: usize,
        integrand_name: impl AsRef<str>,
    ) -> Result<crate::processes::process::ResolvedIntegrandRef<'_>> {
        let process = &self.processes[process_id];
        process.get_integrand(integrand_name)
    }

    pub fn get_integrand_mut(
        &mut self,
        process_id: usize,
        integrand_name: impl AsRef<str>,
    ) -> Result<&mut crate::integrands::process::ProcessIntegrand> {
        let process = &mut self.processes[process_id];
        process.get_integrand_mut(integrand_name)
    }

    /// Export graph DOT files below the checked root's `processes` directory.
    ///
    /// Construct `path` at the user-selected export root so descendant validation
    /// remains anchored there while process and integrand directories are created.
    ///
    /// # Errors
    /// Returns an error for invalid names, symlinks or incompatible file types below
    /// the root, graph serialization failures, or filesystem failures.
    ///
    /// ```no_run
    /// use alphal00p_persistence::OutputPath;
    /// use eyre::Result;
    /// use gammalooprs::processes::{DotExportSettings, ProcessList};
    /// # fn export(processes: &ProcessList) -> Result<()> {
    /// let output = OutputPath::new("exports".as_ref())?;
    /// processes.export_dot(&output, &DotExportSettings::default())?;
    /// # Ok(()) }
    /// ```
    pub fn export_dot(&self, path: &OutputPath, settings: &DotExportSettings) -> Result<()> {
        let path = path.join(PROCESSES_DIRECTORY)?;
        path.ensure_directory()?;

        for p in self.processes.iter() {
            p.export_dot(&path, settings)?;
        }
        Ok(())
    }

    /// Export selected graphs' UV forests below the checked export root.
    ///
    /// The process and generated integrand must exist. Graph IDs refer to that
    /// integrand; structure-only output uses `computed: false`.
    ///
    /// # Errors
    /// Returns an error for an invalid process, integrand or graph selection, an
    /// ungenerated integrand, forest evaluation failures, rejected output paths,
    /// or filesystem failures.
    ///
    /// ```no_run
    /// use alphal00p_persistence::OutputPath;
    /// use eyre::Result;
    /// use gammalooprs::model::Model;
    /// use gammalooprs::processes::ProcessList;
    /// use gammalooprs::uv::export::UVForestExportSettings;
    /// # fn export(processes: &ProcessList, model: &Model) -> Result<()> {
    /// let output = OutputPath::new("exports".as_ref())?;
    /// let settings = UVForestExportSettings { computed: false };
    /// processes.export_uv_forests(model, &output, 0, "LO", &[0], &settings)?;
    /// # Ok(()) }
    /// ```
    pub fn export_uv_forests(
        &self,
        model: &Model,
        path: &OutputPath,
        process_id: usize,
        integrand_name: &str,
        graph_ids: &[usize],
        settings: &UVForestExportSettings,
    ) -> Result<()> {
        let path = path.join(PROCESSES_DIRECTORY)?;
        path.ensure_directory()?;
        let process = self.processes.get(process_id).ok_or_else(|| {
            eyre!(
                "Process id {} is out of range; there are {} processes",
                process_id,
                self.processes.len()
            )
        })?;
        process.export_uv_forests(model, &path, integrand_name, graph_ids, settings)
    }

    /// Export generated integrands and loaders below the checked export root.
    ///
    /// Output is grouped by process kind, process name and integrand name under
    /// `processes`. Settings select the archive format, precision and loader mode.
    ///
    /// # Errors
    /// Returns an error for ungenerated integrands, unsupported export settings,
    /// serialization failures, rejected output paths, or filesystem failures.
    ///
    /// ```no_run
    /// use alphal00p_persistence::OutputPath;
    /// use eyre::Result;
    /// use gammalooprs::processes::{ProcessList, StandaloneExportSettings};
    /// # fn export(processes: &ProcessList, settings: &StandaloneExportSettings) -> Result<()> {
    /// let output = OutputPath::new("exports".as_ref())?;
    /// processes.export_standalone(&output, settings)?;
    /// # Ok(()) }
    /// ```
    pub fn export_standalone(
        &self,
        path: &OutputPath,
        settings: &StandaloneExportSettings,
    ) -> Result<()> {
        let path = path.join(PROCESSES_DIRECTORY)?;
        path.ensure_directory()?;

        for p in self.processes.iter() {
            p.export_standalone(&path, settings)?;
        }
        Ok(())
    }

    pub fn add_process(&mut self, mut process: Process) {
        process.definition.process_id = self.processes.len();
        self.processes.push(process);
    }

    /// Remove a process and restore dense positional IDs.
    ///
    /// # Errors
    /// An out-of-range index leaves the list unchanged.
    ///
    /// ```
    /// # use gammalooprs::processes::ProcessList;
    /// let mut processes = ProcessList::new();
    /// assert!(processes.remove_process(0).is_err());
    /// ```
    pub fn remove_process(&mut self, index: usize) -> StdResult<Process, ProcessIndexError> {
        if index >= self.processes.len() {
            return Err(ProcessIndexError {
                index,
                len: self.processes.len(),
            });
        }
        let removed = self.processes.remove(index);
        self.reindex_processes();
        Ok(removed)
    }

    fn reindex_processes(&mut self) {
        for (id, process) in self.processes.iter_mut().enumerate() {
            process.definition.process_id = id;
        }
    }

    ///preprocesses the process list according to the settings
    pub fn preprocess(
        &mut self,
        model: &Model,
        settings: &GlobalSettings,
        locked_runtime_settings: &LockedRuntimeSettings,
        thread_pool: &ThreadPool,
    ) -> Result<Vec<GeneratedGraphReport>> {
        // Reject unsupported settings before any process changes its graphs or history.
        for process in &self.processes {
            settings
                .generation
                .validate_for_process(process.definition.generation_type)?;
        }
        let mut reports = Vec::new();
        for process in self.processes.iter_mut() {
            reports.extend(process.preprocess(
                model,
                settings,
                locked_runtime_settings,
                thread_pool,
            )?);
        }
        Ok(reports)
    }

    pub fn generate_integrands(
        &mut self,
        model: &Model,
        global_settings: &GlobalSettings,
        runtime_default: LockedRuntimeSettings,
        thread_pool: &ThreadPool,
    ) -> Result<Vec<GeneratedGraphReport>> {
        let mut reports = Vec::new();
        for process in &mut self.processes {
            reports.extend(process.generate_integrands(
                model,
                global_settings,
                runtime_default,
                thread_pool,
            )?);
        }
        Ok(reports)
    }

    pub fn find_process(&self, process_id: Option<usize>) -> Result<usize> {
        if let Some(id) = process_id {
            if id >= self.processes.len() {
                return Err(color_eyre::eyre::eyre!(
                    "Invalid process id {}. Number of processes: {}",
                    id,
                    self.processes.len()
                ));
            }
            Ok(id)
        } else if self.processes.is_empty() {
            Err(color_eyre::eyre::eyre!("No processes generated yet."))
        } else if self.processes.len() > 1 {
            Err(color_eyre::eyre::eyre!(
                "There are {} processes available. Please specify a process.",
                self.processes.len()
            ))
        } else {
            Ok(0)
        }
    }

    pub fn find_integrand(
        &self,
        process_id: Option<usize>,
        integrand_name: Option<&String>,
    ) -> Result<(usize, String)> {
        let p_id = self.find_process(process_id)?;
        let integrand_name = self.processes[p_id]
            .collection
            .find_integrand(integrand_name.cloned())?;
        Ok((p_id, integrand_name))
    }
}

pub mod process;
pub use process::*;
pub mod amplitude;
pub use amplitude::*;
pub mod cross_section;
pub use cross_section::*;
pub mod threshold_counterterms;
pub use threshold_counterterms::*;

#[cfg(test)]
mod no_directive_equivalence_tests;
#[cfg(test)]
mod raised_cross_section_tests;

#[cfg(test)]
mod tests {

    use symbolica::state::State;

    use crate::{
        GammaLoopContextContainer,
        graph::Graph,
        initialisation::test_initialise,
        settings::{
            RuntimeSettings,
            global::{
                CompilationMode, CompilationOptimizationLevel, GammaloopCompileOptions,
                GenerationSettings, TropicalSubgraphTableSettings,
            },
            runtime::LockedRuntimeSettings,
        },
        utils::load_generic_model,
    };

    use super::{AmplitudeGraph, EvaluatorSettings};

    #[test]
    fn numerator_inlining_defaults_and_explicit_opt_out_round_trip() {
        use crate::utils::serde_utils::ShowDefaultsGuard;

        assert!(EvaluatorSettings::default().inline_numerator_functions);
        for input in ["", "n_cores = 2"] {
            let settings: EvaluatorSettings = toml::from_str(input).unwrap();
            assert!(settings.inline_numerator_functions);
        }

        {
            let _guard = ShowDefaultsGuard::new(false);
            let settings = EvaluatorSettings::default();
            let serialized = toml::to_string(&settings).unwrap();
            assert!(!serialized.contains("inline_numerator_functions"));
            assert_eq!(
                toml::from_str::<EvaluatorSettings>(&serialized).unwrap(),
                settings
            );

            let settings: EvaluatorSettings =
                toml::from_str("inline_numerator_functions = false").unwrap();
            let serialized = toml::to_string(&settings).unwrap();
            assert!(serialized.contains("inline_numerator_functions = false"));
            assert_eq!(
                toml::from_str::<EvaluatorSettings>(&serialized).unwrap(),
                settings
            );
        }

        let _guard = ShowDefaultsGuard::new(true);
        let serialized = toml::to_string(&EvaluatorSettings::default()).unwrap();
        assert!(serialized.contains("inline_numerator_functions = true"));
    }

    #[test]
    fn test_encode_decode_amplitude_graph() {
        test_initialise().unwrap();
        // Use the model routing produced by graph import for the persistence roundtrip.
        let model = load_generic_model("sm");

        let graph: Graph = crate::graph::Graph::import_dot(
            stringify!(
                digraph G{
                    e1      [style=invis]
                    e2      [style=invis]
                    e1 -> n1  [particle=H]
                    n4 -> e2    [particle=H]
                    n1 -> n2    [particle=H]
                    n1 -> n3    [particle=H]
                    n2 -> n3    [particle=t]
                    n3 -> n4    [particle=t]
                    n4 -> n2    [particle=t]
                }
            ),
            &crate::utils::load_generic_model("sm"),
        )
        .map(|mut graphs| {
            assert_eq!(graphs.len(), 1);
            graphs.remove(0)
        })
        .unwrap();
        let mut amplitude: AmplitudeGraph = AmplitudeGraph::new(graph.clone());

        amplitude
            .preprocess(
                &model,
                &GenerationSettings {
                    compile: GammaloopCompileOptions {
                        compilation_mode: CompilationMode::Cpp,
                        jit_direct_translation: true,
                        fast_math: false,
                        optimization_level: CompilationOptimizationLevel::O0,
                        unsafe_math: false,
                        compiler: crate::settings::global::default_external_compiler().to_owned(),
                        custom: Vec::new(),
                    },
                    tropical_subgraph_table: TropicalSubgraphTableSettings {
                        panic_on_fail: false,
                        target_omega: 1.0,
                        ..Default::default()
                    },
                    ..Default::default()
                },
                &LockedRuntimeSettings::from(&RuntimeSettings::default()),
            )
            .unwrap();

        let mut symbols = Vec::new();
        State::export(&mut symbols).unwrap();
        let state_map = State::import(&mut std::io::Cursor::new(symbols), None).unwrap();

        let context = GammaLoopContextContainer {
            model: &model,
            state_map: &state_map,
        };

        let encoded_amplitude =
            bincode::encode_to_vec(&amplitude, bincode::config::standard()).unwrap();

        let decoded: AmplitudeGraph = bincode::decode_from_slice_with_context(
            &encoded_amplitude,
            bincode::config::standard(),
            context,
        )
        .expect("amplitude decode failed")
        .0;

        assert_eq!(decoded.graph.name, amplitude.graph.name);
        assert_eq!(decoded.graph.overall_factor, amplitude.graph.overall_factor);
        assert_eq!(
            decoded.graph.loop_momentum_basis,
            amplitude.graph.loop_momentum_basis
        );
        assert_eq!(
            decoded.graph.param_builder.reps,
            amplitude.graph.param_builder.reps
        );
        assert_eq!(
            decoded.graph.param_builder.values,
            amplitude.graph.param_builder.values
        );
        assert_eq!(
            decoded.derived_data.all_mighty_integrand,
            amplitude.derived_data.all_mighty_integrand
        );
        assert_eq!(
            decoded.derived_data.all_mighty_numerators,
            amplitude.derived_data.all_mighty_numerators
        );
    }
}
