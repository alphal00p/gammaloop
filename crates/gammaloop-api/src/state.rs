use std::{
    collections::{BTreeMap, BTreeSet, HashSet},
    fs::{self},
    io::{self, IsTerminal},
    ops::ControlFlow,
    path::{Path, PathBuf},
    str::FromStr,
    sync::{
        atomic::{AtomicBool, AtomicU64, Ordering},
        Arc, Mutex,
    },
    thread,
    time::{Duration, Instant},
};

use clap::Args;
use color_eyre::{Result, Section};
use colored::Colorize;
use eyre::{eyre, Context};
use gammalooprs::{
    processes::{Amplitude, CrossSection},
    utils::serde_utils::IsDefault,
};
use linnet::half_edge::subgraph::SubGraphLike;
use schemars::{schema_for, JsonSchema, Schema};
use serde::{Deserialize, Serialize};
use spenso::algebra::complex::Complex;
use symbolica::numerical_integration::Sample;
use sysinfo::{get_current_pid, ProcessRefreshKind, ProcessesToUpdate, System};
use toml::Value as TomlValue;
use tracing::{debug, info, info_span, Span};
use tracing_indicatif::span_ext::IndicatifSpanExt;

use gammalooprs::{
    clear_interrupt_request,
    feyngen::GenerationType,
    graph::Graph,
    initialisation::initialise,
    integrands::{process::ProcessIntegrand, HasIntegrand},
    is_interrupt_requested,
    model::{InputParamCard, Model, SerializableInputParamCard, UFOSymbol},
    processes::{
        begin_phase, merge_generated_graph_reports, DotExportSettings, GeneratedGraphReport,
        GenerationProcessKind, GenerationProgressMode, GenerationProgressModeGuard,
        GenerationProgressObserver, GenerationProgressObserverGuard, GenerationProgressPhase,
        GraphGenerationStats, GraphGroupSelectionMode, GraphGroupSelectionPlan,
        GraphGroupSelectionReport, GraphGroupSelectionSpec, NamedGraphGenerationReport, Process,
        ProcessCollection, ProcessDefinition, ProcessList,
    },
    settings::{
        global::GenerationSettings, runtime::LockedRuntimeSettings, GlobalSettings, RuntimeSettings,
    },
    utils::{
        serde_utils::{get_schema_folder, SmartSerde},
        tracing::{init_bench_tracing, init_test_tracing},
        F,
    },
    GammaLoopContextContainer,
};

use crate::{
    command_parser::{normalize_clap_args, split_command_line},
    commands::{save::SaveState, Commands},
    integrand_info::{collect_integrand_info, IntegrandInfo},
    model_parameters::{external_model_parameter_type, validate_model_parameter_type},
    render_smart_toml,
    tracing::{set_file_log_filter, set_log_style, set_stderr_log_filter},
    CLISettings,
};

#[derive(Debug, Clone, Copy, Default)]
pub struct GenerationResourceSummary {
    pub peak_ram_bytes: u64,
    pub generation_cores: usize,
}

#[derive(Debug, Clone, Default)]
pub struct GenerationReports {
    pub reports: Vec<GeneratedGraphReport>,
    pub resources: GenerationResourceSummary,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct IntegrandGenerationSummaryKey {
    pub process_id: usize,
    pub integrand_name: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IntegrandGenerationSummary {
    pub peak_ram_bytes: u64,
    pub reports: Vec<GeneratedGraphReport>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SelectedGraphGroups {
    pub source_process_id: usize,
    pub source_process_name: String,
    pub source_integrand_name: String,
    pub process_id: usize,
    pub process_name: String,
    pub integrand_name: String,
    pub report: GraphGroupSelectionReport,
    pub copied_to_output: bool,
    pub replaced_existing_target: bool,
    pub removed_target_artifacts: bool,
    pub discarded_generated_integrand: bool,
    pub removed_generated_artifacts: bool,
    pub removed_generation_summary: bool,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct GraphGroupSelectionTarget {
    pub output_process_name: Option<String>,
    pub output_integrand_name: Option<String>,
    pub clear_existing_processes: bool,
}

impl GraphGroupSelectionTarget {
    pub fn in_place() -> Self {
        Self::default()
    }

    pub fn copy(
        output_process_name: Option<String>,
        output_integrand_name: Option<String>,
        clear_existing_processes: bool,
    ) -> Self {
        Self {
            output_process_name,
            output_integrand_name,
            clear_existing_processes,
        }
    }

    fn is_copy_mode(&self) -> bool {
        self.output_process_name.is_some() || self.output_integrand_name.is_some()
    }
}

#[derive(Debug, Clone, Copy)]
pub struct GraphGroupSelectionContext<'a> {
    pub generation_settings: &'a GenerationSettings,
    pub state_folder: &'a Path,
    pub read_only_state: bool,
}

impl<'a> GraphGroupSelectionContext<'a> {
    pub fn new(
        generation_settings: &'a GenerationSettings,
        state_folder: &'a Path,
        read_only_state: bool,
    ) -> Self {
        Self {
            generation_settings,
            state_folder,
            read_only_state,
        }
    }
}

struct IntegrandCopyInsertion {
    process_id: usize,
    process_name: String,
    integrand_name: String,
    replaced_existing_target: bool,
    removed_target_artifacts: bool,
}

struct SelectionSource {
    process_id: usize,
    process_name: String,
    integrand_name: String,
}

enum IntegrandCopyPayload {
    Amplitude(Amplitude),
    CrossSection(CrossSection),
}

impl IntegrandCopyPayload {
    fn rename(&mut self, new_name: &str) {
        match self {
            Self::Amplitude(amplitude) => {
                amplitude.name = new_name.to_string();
                rename_process_integrand(amplitude.integrand.as_mut(), new_name);
            }
            Self::CrossSection(cross_section) => {
                cross_section.name = new_name.to_string();
                rename_process_integrand(cross_section.integrand.as_mut(), new_name);
            }
        }
    }

    fn apply_graph_group_selection(&mut self, plan: &GraphGroupSelectionPlan) -> Result<()> {
        match self {
            Self::Amplitude(amplitude) => amplitude.apply_graph_group_selection(plan),
            Self::CrossSection(cross_section) => cross_section.apply_graph_group_selection(plan),
        }
    }

    fn is_compatible_with(&self, collection: &ProcessCollection) -> bool {
        matches!(
            (self, collection),
            (Self::Amplitude(_), ProcessCollection::Amplitudes(_))
                | (Self::CrossSection(_), ProcessCollection::CrossSections(_))
        )
    }

    fn kind_name(&self) -> &'static str {
        match self {
            Self::Amplitude(_) => "amplitudes",
            Self::CrossSection(_) => "cross sections",
        }
    }

    fn into_collection(self) -> ProcessCollection {
        match self {
            Self::Amplitude(amplitude) => {
                let mut collection = ProcessCollection::Amplitudes(BTreeMap::new());
                collection.add_amplitude(amplitude);
                collection
            }
            Self::CrossSection(cross_section) => {
                let mut collection = ProcessCollection::CrossSections(BTreeMap::new());
                collection.add_cross_section(cross_section);
                collection
            }
        }
    }

    fn insert_into_process(self, process: &mut Process, process_name: &str) -> Result<()> {
        match (self, &mut process.collection) {
            (Self::Amplitude(amplitude), ProcessCollection::Amplitudes(amplitudes)) => {
                amplitudes.insert(amplitude.name.clone(), amplitude);
                Ok(())
            }
            (
                Self::CrossSection(cross_section),
                ProcessCollection::CrossSections(cross_sections),
            ) => {
                cross_sections.insert(cross_section.name.clone(), cross_section);
                Ok(())
            }
            (payload, _) => Err(eyre!(
                "Destination process '{}' exists but does not contain {}",
                process_name,
                payload.kind_name()
            )),
        }
    }
}

struct GenerationMonitor {
    peak_ram_bytes: Arc<AtomicU64>,
    current_ram_bytes: Arc<AtomicU64>,
    stop_requested: Arc<AtomicBool>,
    handle: Option<thread::JoinHandle<()>>,
}

struct AggregateGenerationProgressReporter {
    progress_span: Span,
    current_ram_bytes: Arc<AtomicU64>,
    peak_ram_bytes: Arc<AtomicU64>,
    state: Mutex<AggregateGenerationProgressState>,
}

#[derive(Default)]
struct AggregateGenerationProgressState {
    phase: Option<GenerationProgressPhase>,
    kind: Option<GenerationProcessKind>,
    process: String,
    integrand: String,
    total_graphs: usize,
    done_graphs: usize,
    total_cuts: Option<usize>,
    done_cuts: usize,
    discovered_st_cuts: usize,
    discovered_valid_cuts: usize,
    active_graphs: BTreeSet<String>,
    last_graph: Option<String>,
    stats: GraphGenerationStats,
}

impl AggregateGenerationProgressReporter {
    fn new(
        current_ram_bytes: Arc<AtomicU64>,
        peak_ram_bytes: Arc<AtomicU64>,
        eta_warmup_steps: u64,
    ) -> Arc<Self> {
        let span = info_span!(
            "Integrand generation",
            indicatif.pb_show = true,
            indicatif.pb_msg = "Starting integrand generation"
        );
        span.pb_set_style(
            &gammalooprs::utils::long_running_progress_style_with_eta_warmup(eta_warmup_steps),
        );
        span.pb_start();
        span.pb_set_length(0);
        span.pb_set_position(0);
        span.pb_tick();
        Self::new_with_span(current_ram_bytes, peak_ram_bytes, span)
    }

    #[cfg(test)]
    fn new_hidden(current_ram_bytes: Arc<AtomicU64>, peak_ram_bytes: Arc<AtomicU64>) -> Arc<Self> {
        Self::new_with_span(current_ram_bytes, peak_ram_bytes, Span::none())
    }

    fn new_with_span(
        current_ram_bytes: Arc<AtomicU64>,
        peak_ram_bytes: Arc<AtomicU64>,
        progress_span: Span,
    ) -> Arc<Self> {
        Arc::new(Self {
            progress_span,
            current_ram_bytes,
            peak_ram_bytes,
            state: Mutex::new(AggregateGenerationProgressState::default()),
        })
    }

    fn fixed_field(value: &str, width: usize) -> String {
        let mut chars = value.chars().collect::<Vec<_>>();
        if chars.len() > width {
            chars.truncate(width.saturating_sub(1));
            chars.push('~');
        }
        format!("{:<width$}", chars.into_iter().collect::<String>())
    }

    fn progress_memory(bytes: u64) -> String {
        const MIB: f64 = 1024.0 * 1024.0;
        const GIB: f64 = MIB * 1024.0;
        let value = bytes as f64;
        if value >= GIB {
            format!("{:>6.2} GiB", value / GIB)
        } else {
            format!("{:>6.0} MiB", value / MIB)
        }
    }

    fn progress_percent(percent: f64) -> String {
        if percent < 0.01 {
            "0%".to_string()
        } else if percent < 0.1 {
            format!("{percent:.3}%")
        } else if percent < 1.0 {
            format!("{percent:.2}%")
        } else if percent < 10.0 {
            format!("{percent:.1}%")
        } else {
            format!("{percent:.0}%")
        }
    }

    fn progress_time_share(duration: Duration, total: Duration) -> String {
        if total.is_zero() {
            "--".to_string()
        } else {
            Self::progress_percent(duration.as_secs_f64() * 100.0 / total.as_secs_f64())
        }
    }

    fn progress_units(state: &AggregateGenerationProgressState) -> (u64, u64) {
        match (state.phase, state.total_cuts) {
            (Some(GenerationProgressPhase::GraphGeneration), Some(total_cuts)) => {
                let total = total_cuts.saturating_add(state.total_graphs) as u64;
                let completed = state
                    .done_cuts
                    .saturating_add(state.done_graphs)
                    .min(total_cuts.saturating_add(state.total_graphs))
                    as u64;
                (completed, total)
            }
            _ => (state.done_graphs as u64, state.total_graphs as u64),
        }
    }

    fn graph_progress_counts(state: &AggregateGenerationProgressState) -> (usize, usize, usize) {
        let active = state
            .active_graphs
            .len()
            .min(state.total_graphs.saturating_sub(state.done_graphs));
        (state.done_graphs, active, state.total_graphs)
    }

    fn refresh(&self, state: &AggregateGenerationProgressState) {
        let (completed_units, total_units) = Self::progress_units(state);
        self.progress_span.pb_set_length(total_units);
        self.progress_span.pb_set_position(completed_units);

        let current_ram = Self::progress_memory(self.current_ram_bytes.load(Ordering::Relaxed));
        let peak_ram = Self::progress_memory(self.peak_ram_bytes.load(Ordering::Relaxed));
        let kind = match state.kind {
            Some(GenerationProcessKind::Amplitude) => "AMP",
            Some(GenerationProcessKind::CrossSection) => "XS",
            None => "GEN",
        };
        let phase = match state.phase {
            Some(GenerationProgressPhase::GraphPreprocessing) => "prep",
            Some(GenerationProgressPhase::GraphGeneration) => "graphs",
            Some(GenerationProgressPhase::Backend) => "backend",
            None => "gen",
        };
        let last_graph = state.last_graph.as_deref().unwrap_or("-");
        let cut_context = match (state.phase, state.kind, state.total_cuts) {
            (
                Some(GenerationProgressPhase::GraphPreprocessing),
                Some(GenerationProcessKind::Amplitude),
                _,
            ) => None,
            (Some(GenerationProgressPhase::GraphPreprocessing), _, _) => Some(format!(
                "{:>4}/{:<4}",
                state.discovered_valid_cuts, state.discovered_st_cuts
            )),
            (_, _, Some(total_cuts)) => Some(format!("{:>4}/{:<4}", state.done_cuts, total_cuts)),
            _ => Some(format!("{:>4}/{:<4}", "-", "-")),
        };
        let phase = Self::fixed_field(phase, 7).bold().blue();
        let kind = Self::fixed_field(kind, 3).cyan();
        let identifier = if state.process.is_empty() {
            state.integrand.clone()
        } else {
            format!("{}@{}", state.integrand, state.process)
        };
        let identifier = Self::fixed_field(&identifier, 24).yellow();
        let last_graph = Self::fixed_field(last_graph, 8).green();
        let (done_graphs, active_graphs, total_graphs) = Self::graph_progress_counts(state);
        let graph_ratio = format!(
            "{} {} / {}",
            format!("{done_graphs:>3}").green(),
            format!("({active_graphs:>3})").dimmed(),
            format!("{total_graphs:<3}").green(),
        );
        let ram = format!("{current_ram}/{peak_ram}").yellow();
        let phase_detail = if state.phase == Some(GenerationProgressPhase::GraphPreprocessing) {
            let steps = match state.kind {
                Some(GenerationProcessKind::Amplitude) => {
                    "CFFs + LMBs + integrands + threshold CTs + tropical samplers"
                }
                _ => "cuts + CFFs + LMBs + integrands + threshold CTs",
            };
            format!("{} {}", "step".bold().blue(), steps.cyan())
        } else {
            let total_time = state.stats.total_time;
            let expression_share =
                Self::progress_time_share(state.stats.expression_build_time(), total_time);
            let spenso_share =
                Self::progress_time_share(state.stats.evaluator_spenso_time, total_time);
            let symbolica_share =
                Self::progress_time_share(state.stats.evaluator_symbolica_time, total_time);
            let compile_share =
                Self::progress_time_share(state.stats.evaluator_compile_time, total_time);
            format!(
                "{} {} {} / {} {} / {} {} / {} {}",
                "time".bold().blue(),
                "expr".bold().blue(),
                expression_share.magenta(),
                "spenso".bold().blue(),
                spenso_share.cyan(),
                "eval".bold().blue(),
                symbolica_share.green(),
                "compile".bold().blue(),
                compile_share.yellow(),
            )
        };
        let mut sections = vec![
            format!("{} {} {}", phase, kind, identifier),
            format!("{} {}", "g".bold().blue(), graph_ratio),
            format!("{} {}", "last".bold().blue(), last_graph),
        ];
        if let Some(cut_context) = cut_context {
            sections.push(format!("{} {}", "cut".bold().blue(), cut_context.green()));
        }
        sections.push(format!("{} {}", "ram".bold().blue(), ram));
        sections.push(phase_detail);
        self.progress_span.pb_set_message(&sections.join(" | "));
        self.progress_span.pb_tick();
    }
}

impl GenerationProgressObserver for AggregateGenerationProgressReporter {
    fn begin_phase(
        &self,
        phase: GenerationProgressPhase,
        kind: GenerationProcessKind,
        process: &str,
        integrand: &str,
        total_graphs: usize,
        total_cuts: Option<usize>,
    ) {
        let mut state = self
            .state
            .lock()
            .expect("aggregate generation progress state mutex is poisoned");
        state.phase = Some(phase);
        state.kind = Some(kind);
        state.process.clear();
        state.process.push_str(process);
        state.integrand.clear();
        state.integrand.push_str(integrand);
        state.total_graphs = total_graphs;
        state.done_graphs = 0;
        state.total_cuts = total_cuts;
        state.done_cuts = 0;
        state.discovered_st_cuts = 0;
        state.discovered_valid_cuts = 0;
        state.active_graphs.clear();
        state.last_graph = None;
        state.stats = GraphGenerationStats::default();
        self.progress_span.pb_reset_elapsed();
        self.refresh(&state);
    }

    fn graph_started(
        &self,
        _kind: GenerationProcessKind,
        _integrand: &str,
        graph: &str,
        _cut_count: Option<usize>,
    ) {
        let mut state = self
            .state
            .lock()
            .expect("aggregate generation progress state mutex is poisoned");
        state.active_graphs.insert(graph.to_string());
        state.last_graph = Some(graph.to_string());
        self.refresh(&state);
    }

    fn graph_finished(
        &self,
        _kind: GenerationProcessKind,
        _integrand: &str,
        graph: &str,
        stats: &GraphGenerationStats,
        completed_cuts: Option<usize>,
    ) {
        let mut state = self
            .state
            .lock()
            .expect("aggregate generation progress state mutex is poisoned");
        state.done_graphs = state.done_graphs.saturating_add(1).min(state.total_graphs);
        state.done_cuts += completed_cuts.unwrap_or(0);
        state.active_graphs.remove(graph);
        state.last_graph = Some(graph.to_string());
        state.stats.merge_in_place(stats);
        self.refresh(&state);
    }

    fn cuts_discovered(
        &self,
        _integrand: &str,
        graph: &str,
        st_cut_count: usize,
        valid_cut_count: usize,
    ) {
        let mut state = self
            .state
            .lock()
            .expect("aggregate generation progress state mutex is poisoned");
        state.discovered_st_cuts += st_cut_count;
        state.discovered_valid_cuts += valid_cut_count;
        state.last_graph = Some(graph.to_string());
        self.refresh(&state);
    }

    fn cut_finished(&self, _integrand: &str, graph: &str, cut_count: usize) {
        let mut state = self
            .state
            .lock()
            .expect("aggregate generation progress state mutex is poisoned");
        state.done_cuts = state.done_cuts.saturating_add(cut_count);
        if let Some(total_cuts) = state.total_cuts {
            state.done_cuts = state.done_cuts.min(total_cuts);
        }
        state.last_graph = Some(graph.to_string());
        self.refresh(&state);
    }

    fn backend_started(&self, kind: GenerationProcessKind, integrand: &str, graph_count: usize) {
        let mut state = self
            .state
            .lock()
            .expect("aggregate generation progress state mutex is poisoned");
        state.phase = Some(GenerationProgressPhase::Backend);
        state.kind = Some(kind);
        state.integrand.clear();
        state.integrand.push_str(integrand);
        state.total_graphs = graph_count;
        state.done_graphs = 0;
        state.total_cuts = None;
        state.done_cuts = 0;
        state.active_graphs.clear();
        state.last_graph = None;
        self.progress_span.pb_reset_elapsed();
        self.refresh(&state);
    }

    fn backend_finished(&self, _kind: GenerationProcessKind, _integrand: &str, elapsed: Duration) {
        let mut state = self
            .state
            .lock()
            .expect("aggregate generation progress state mutex is poisoned");
        state.done_graphs = state.total_graphs;
        state.stats.total_time += elapsed;
        state.stats.evaluator_compile_time += elapsed;
        self.refresh(&state);
    }
}

impl GenerationMonitor {
    const POLL_INTERVAL: Duration = Duration::from_millis(100);

    fn start() -> Result<Self> {
        let pid = get_current_pid().map_err(|err| eyre!("Failed to resolve current pid: {err}"))?;
        let peak_ram_bytes = Arc::new(AtomicU64::new(0));
        let current_ram_bytes = Arc::new(AtomicU64::new(0));
        let stop_requested = Arc::new(AtomicBool::new(false));
        let peak_ram_bytes_for_thread = Arc::clone(&peak_ram_bytes);
        let current_ram_bytes_for_thread = Arc::clone(&current_ram_bytes);
        let stop_requested_for_thread = Arc::clone(&stop_requested);

        let handle = thread::Builder::new()
            .name("generation-ram-monitor".to_string())
            .spawn(move || {
                let mut system = System::new();
                loop {
                    if stop_requested_for_thread.load(Ordering::Relaxed) || is_interrupt_requested()
                    {
                        break;
                    }

                    system.refresh_processes_specifics(
                        ProcessesToUpdate::Some(&[pid]),
                        true,
                        ProcessRefreshKind::nothing().with_memory(),
                    );
                    if let Some(process) = system.process(pid) {
                        let memory = process.memory();
                        current_ram_bytes_for_thread.store(memory, Ordering::Relaxed);
                        peak_ram_bytes_for_thread.fetch_max(memory, Ordering::Relaxed);
                    }
                    thread::sleep(Self::POLL_INTERVAL);
                }
            })
            .map_err(|err| eyre!("Failed to spawn generation RAM monitor: {err}"))?;

        Ok(Self {
            peak_ram_bytes,
            current_ram_bytes,
            stop_requested,
            handle: Some(handle),
        })
    }

    fn current_ram_bytes(&self) -> Arc<AtomicU64> {
        Arc::clone(&self.current_ram_bytes)
    }

    fn peak_ram_bytes(&self) -> Arc<AtomicU64> {
        Arc::clone(&self.peak_ram_bytes)
    }

    fn finish(&mut self) -> u64 {
        self.stop_requested.store(true, Ordering::Relaxed);
        if let Some(handle) = self.handle.take() {
            let _ = handle.join();
        }
        self.peak_ram_bytes.load(Ordering::Relaxed)
    }
}

impl Drop for GenerationMonitor {
    fn drop(&mut self) {
        let _ = self.finish();
    }
}

#[derive(Debug, Clone, PartialEq, Eq, JsonSchema)]
pub enum ProcessRef {
    Id(usize),
    Name(String),
    Unqualified(String),
}

impl Serialize for ProcessRef {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        match self {
            ProcessRef::Id(id) => serializer.serialize_u64(*id as u64),
            ProcessRef::Name(name) => serializer.serialize_str(&format!("name:{name}")),
            ProcessRef::Unqualified(value) => serializer.serialize_str(value),
        }
    }
}

impl<'de> Deserialize<'de> for ProcessRef {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        use serde::de::{self, Visitor};
        use std::fmt;

        struct ProcessRefVisitor;

        impl<'de> Visitor<'de> for ProcessRefVisitor {
            type Value = ProcessRef;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("a process reference string or numeric id")
            }

            fn visit_u64<E>(self, value: u64) -> Result<Self::Value, E>
            where
                E: de::Error,
            {
                Ok(ProcessRef::Id(value as usize))
            }

            fn visit_str<E>(self, value: &str) -> Result<Self::Value, E>
            where
                E: de::Error,
            {
                ProcessRef::from_str(value).map_err(E::custom)
            }

            fn visit_string<E>(self, value: String) -> Result<Self::Value, E>
            where
                E: de::Error,
            {
                self.visit_str(&value)
            }
        }

        deserializer.deserialize_any(ProcessRefVisitor)
    }
}

#[cfg(feature = "python_api")]
impl<'a, 'py> pyo3::FromPyObject<'a, 'py> for ProcessRef {
    type Error = pyo3::PyErr;

    fn extract(obj: pyo3::Borrowed<'a, 'py, pyo3::types::PyAny>) -> pyo3::PyResult<Self> {
        if let Ok(process_id) = <usize as pyo3::FromPyObject<'a, 'py>>::extract(obj) {
            return Ok(ProcessRef::Id(process_id));
        }

        let selector = <String as pyo3::FromPyObject<'a, 'py>>::extract(obj).map_err(|_| {
            pyo3::exceptions::PyTypeError::new_err(
                "process selectors must be either an integer process id or a string selector",
            )
        })?;
        ProcessRef::from_str(&selector).map_err(pyo3::exceptions::PyValueError::new_err)
    }
}

impl FromStr for ProcessRef {
    type Err = String;

    fn from_str(value: &str) -> std::result::Result<Self, Self::Err> {
        if let Some(rest) = value.strip_prefix('#') {
            let id = rest
                .parse::<usize>()
                .map_err(|_| format!("Invalid process id in '{value}'"))?;
            return Ok(ProcessRef::Id(id));
        }
        if let Some(rest) = value.strip_prefix("name:") {
            if rest.is_empty() {
                return Err("Process name cannot be empty".to_string());
            }
            return Ok(ProcessRef::Name(rest.to_string()));
        }
        if value.is_empty() {
            return Err("Process reference cannot be empty".to_string());
        }
        Ok(ProcessRef::Unqualified(value.to_string()))
    }
}

#[test]
fn try_complicated() {
    "GGHHH3loop_no_iterative_optimization_3L"
        .parse::<ProcessRef>()
        .unwrap();
    // ProcessRef::
}

impl ProcessRef {
    pub fn resolve(&self, process_list: &ProcessList) -> Result<usize> {
        let process_count = process_list.processes.len();
        match self {
            ProcessRef::Id(id) => {
                if *id >= process_count {
                    return Err(eyre!(
                        "Process ID {} invalid, only {} processes available",
                        id,
                        process_count
                    ));
                }
                Ok(*id)
            }
            ProcessRef::Name(name) => process_list
                .processes
                .iter()
                .position(|p| p.definition.folder_name == *name)
                .ok_or_else(|| {
                    eyre!(
                        "No process named '{}'. Use 'display processes' to list available processes.",
                        name
                    )
                }),
            ProcessRef::Unqualified(value) => {
                let name_match = process_list
                    .processes
                    .iter()
                    .position(|p| p.definition.folder_name == *value);
                if let Ok(id) = value.parse::<usize>() {
                    let id_valid = id < process_count;
                    match (id_valid, name_match) {
                        (true, Some(_)) => Err(eyre!(
                            "Ambiguous process reference '{}'. Use '#{}' or 'name:{}' to disambiguate.",
                            value,
                            id,
                            value
                        )),
                        (true, None) => Ok(id),
                        (false, Some(index)) => Ok(index),
                        (false, None) => Err(eyre!(
                            "No process named '{}'. Use 'display processes' to list available processes.",
                            value
                        )),
                    }
                } else if let Some(index) = name_match {
                    Ok(index)
                } else {
                    Err(eyre!(
                        "No process named '{}'. Use 'display processes' to list available processes.",
                        value
                    ))
                }
            }
        }
    }
}

pub trait ProcessListExt {
    fn find_integrand_ref(
        &self,
        process: Option<&ProcessRef>,
        integrand_name: Option<&String>,
    ) -> Result<(usize, String)>;
    fn get_amplitude_mut_ref(
        &mut self,
        process: Option<&ProcessRef>,
        integrand_name: Option<&String>,
    ) -> Result<&mut Amplitude>;
    fn get_cross_section_mut_ref(
        &mut self,
        process: Option<&ProcessRef>,
        integrand_name: Option<&String>,
    ) -> Result<&mut CrossSection>;
}

impl ProcessListExt for ProcessList {
    fn find_integrand_ref(
        &self,
        process: Option<&ProcessRef>,
        integrand_name: Option<&String>,
    ) -> Result<(usize, String)> {
        let process_id = match process {
            Some(process_ref) => process_ref.resolve(self)?,
            None => self.find_process(None)?,
        };
        let integrand_name = self.processes[process_id]
            .collection
            .find_integrand(integrand_name.cloned())
            .with_note(|| format!("in process id {process_id}"))?;
        Ok((process_id, integrand_name))
    }

    fn get_amplitude_mut_ref(
        &mut self,
        process: Option<&ProcessRef>,
        integrand_name: Option<&String>,
    ) -> Result<&mut Amplitude> {
        let (process_id, integrand_name) = self.find_integrand_ref(process, integrand_name)?;
        let process = &mut self.processes[process_id];
        match &mut process.collection {
            ProcessCollection::Amplitudes(amplitudes) => {
                amplitudes.get_mut(&integrand_name).ok_or_else(|| {
                    eyre!(
                        "No amplitude named '{}' in process '{}'",
                        integrand_name,
                        process.definition.folder_name
                    )
                })
            }
            ProcessCollection::CrossSections(_) => Err(eyre!(
                "Process '{}' does not contain amplitudes",
                process.definition.folder_name
            )),
        }
    }

    fn get_cross_section_mut_ref(
        &mut self,
        process: Option<&ProcessRef>,
        integrand_name: Option<&String>,
    ) -> Result<&mut CrossSection> {
        let (process_id, integrand_name) = self.find_integrand_ref(process, integrand_name)?;
        let process = &mut self.processes[process_id];
        match &mut process.collection {
            ProcessCollection::CrossSections(crosssections) => {
                crosssections.get_mut(&integrand_name).ok_or_else(|| {
                    eyre!(
                        "No cross section named '{}' in process '{}'",
                        integrand_name,
                        process.definition.folder_name
                    )
                })
            }
            ProcessCollection::Amplitudes(_) => Err(eyre!(
                "Process '{}' does not contain crosssections",
                process.definition.folder_name
            )),
        }
    }
}

pub trait SyncSettings {
    fn sync_settings(&self) -> Result<()>;
}

impl SyncSettings for CLISettings {
    fn sync_settings(&self) -> Result<()> {
        // println!("Syncing settings {}", self.global.logfile_directive);
        set_file_log_filter(&self.global.logfile_directive)?;
        set_stderr_log_filter(&self.global.display_directive)?;
        set_log_style(self.global.log_style.to_runtime());
        Ok(())
    }
}

// Static flag to control serialization behavior
static SERIALIZE_COMMANDS_AS_STRINGS: AtomicBool = AtomicBool::new(false);

fn is_command_blocks_empty(command_blocks: &[CommandsBlock]) -> bool {
    command_blocks.is_empty()
}

fn should_persist_command(command: &Commands) -> bool {
    !matches!(
        command,
        Commands::Quit(_) | Commands::StartCommandsBlock(_) | Commands::FinishCommandsBlock
    )
}

/// Set whether CommandHistory should serialize as strings when the raw_string is available
pub fn set_serialize_commands_as_strings(value: bool) {
    SERIALIZE_COMMANDS_AS_STRINGS.store(value, std::sync::atomic::Ordering::Relaxed);
}

/// Get the current setting for CommandHistory serialization behavior
pub fn get_serialize_commands_as_strings() -> bool {
    SERIALIZE_COMMANDS_AS_STRINGS.load(std::sync::atomic::Ordering::Relaxed)
}

pub struct SerializeCommandsAsStringsGuard {
    previous: bool,
}

impl SerializeCommandsAsStringsGuard {
    pub fn new(value: bool) -> Self {
        let previous = get_serialize_commands_as_strings();
        set_serialize_commands_as_strings(value);
        Self { previous }
    }
}

impl Drop for SerializeCommandsAsStringsGuard {
    fn drop(&mut self) {
        set_serialize_commands_as_strings(self.previous);
    }
}

/// Represents a command with optional raw string representation
///
/// This struct stores both the parsed command and optionally the original
/// string that was used to create it. This allows for preserving the exact
/// user input while still having access to the structured command data.
#[derive(Debug, Clone, JsonSchema, PartialEq)]
pub struct CommandHistory {
    /// The parsed command
    pub command: Commands,
    /// The original string representation of the command, if available
    pub raw_string: Option<String>,
}

impl Serialize for CommandHistory {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        if get_serialize_commands_as_strings() {
            if let Some(ref raw_string) = self.raw_string {
                raw_string.serialize(serializer)
            } else {
                self.command.serialize(serializer)
            }
        } else {
            self.command.serialize(serializer)
        }
    }
}

impl<'de> Deserialize<'de> for CommandHistory {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        use serde::de::{self, Visitor};
        use std::fmt;

        struct CommandHistoryVisitor;

        impl<'de> Visitor<'de> for CommandHistoryVisitor {
            type Value = CommandHistory;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("a string or a Commands structure")
            }

            fn visit_str<E>(self, value: &str) -> Result<Self::Value, E>
            where
                E: de::Error,
            {
                CommandHistory::from_raw_string(value).map_err(|err| {
                    E::custom(format!(
                        "Failed to parse command string '{}': {}",
                        value, err
                    ))
                })
            }

            fn visit_string<E>(self, value: String) -> Result<Self::Value, E>
            where
                E: de::Error,
            {
                self.visit_str(&value)
            }

            fn visit_seq<A>(self, mut seq: A) -> Result<Self::Value, A::Error>
            where
                A: de::SeqAccess<'de>,
            {
                // Handle TOML array format for enums like [Quit]
                let command =
                    Commands::deserialize(de::value::SeqAccessDeserializer::new(&mut seq))?;
                Ok(CommandHistory {
                    command,
                    raw_string: None,
                })
            }

            fn visit_map<A>(self, map: A) -> Result<Self::Value, A::Error>
            where
                A: de::MapAccess<'de>,
            {
                // Handle map format for enums
                let command = Commands::deserialize(de::value::MapAccessDeserializer::new(map))?;
                Ok(CommandHistory {
                    command,
                    raw_string: None,
                })
            }
        }

        deserializer.deserialize_any(CommandHistoryVisitor)
    }
}

impl CommandHistory {
    /// Create a new CommandHistory with just a command (no raw string)
    pub fn new(command: Commands) -> Self {
        Self {
            command,
            raw_string: None,
        }
    }

    /// Create a new CommandHistory with both command and raw string
    pub fn new_with_raw(command: Commands, raw_string: String) -> Self {
        Self {
            command,
            raw_string: Some(raw_string),
        }
    }

    /// Create a CommandHistory from a command (alias for new)
    pub fn from_command(command: Commands) -> Self {
        Self::new(command)
    }

    /// Parse a raw string into a CommandHistory
    ///
    /// This function attempts to parse the raw string using clap, and if successful,
    /// creates a CommandHistory with both the parsed command and the original string.
    pub fn from_raw_string(raw_string: &str) -> Result<Self, clap::Error> {
        use crate::Repl;
        use clap::error::ErrorKind;
        use clap::Parser;

        let args = split_command_line(raw_string)
            .map(normalize_clap_args)
            .map_err(|_| {
                clap::Error::raw(
                    ErrorKind::InvalidValue,
                    "Could not parse command: unmatched quotes or trailing escape",
                )
            })?;
        let cli = Repl::try_parse_from(
            std::iter::once("gammaloop").chain(args.iter().map(String::as_str)),
        )?;

        Ok(Self::new_with_raw(cli.command, raw_string.into()))
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, JsonSchema, PartialEq)]
#[serde(default, deny_unknown_fields)]
pub struct RunHistory {
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub default_runtime_settings: RuntimeSettings,

    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub cli_settings: CLISettings,

    #[serde(skip_serializing_if = "is_command_blocks_empty")]
    pub command_blocks: Vec<CommandsBlock>,
    // #[serde(with = "serde_yaml::with::singleton_map_recursive")]
    // #[schemars(with = "Vec<CommandHistory>")]
    pub commands: Vec<CommandHistory>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, JsonSchema, PartialEq)]
#[serde(default, deny_unknown_fields)]
pub struct CommandsBlock {
    pub name: String,
    pub commands: Vec<CommandHistory>,
}

#[derive(Debug, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
struct RawRunHistoryToml {
    default_runtime_settings: RuntimeSettings,
    cli_settings: CLISettings,
    command_blocks: Vec<RawCommandsBlockToml>,
    commands: Vec<TomlValue>,
}

#[derive(Debug, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
struct RawCommandsBlockToml {
    name: String,
    commands: Vec<TomlValue>,
}

impl CommandsBlock {
    pub fn semantically_eq(&self, other: &Self) -> bool {
        self.name == other.name
            && self.commands.len() == other.commands.len()
            && self
                .commands
                .iter()
                .zip(other.commands.iter())
                .all(|(left, right)| left.command == right.command)
    }
}

impl SmartSerde for RunHistory {
    fn has_schema_path(&self, online: bool) -> Option<Result<PathBuf>> {
        Some(get_schema_folder(online).map(|f| f.join("runhistory.json")))
    }
}

impl RunHistory {
    pub fn freeze_boot_settings_from(&mut self, boot_run_history: &RunHistory) {
        self.cli_settings.global = boot_run_history.cli_settings.global.clone();
        self.default_runtime_settings = boot_run_history.default_runtime_settings.clone();
    }

    pub fn frozen_boot_settings_match(&self, boot_run_history: &RunHistory) -> bool {
        self.cli_settings.global == boot_run_history.cli_settings.global
            && self.default_runtime_settings == boot_run_history.default_runtime_settings
    }

    /// Add a command to the run history
    pub fn push(&mut self, command: Commands) {
        self.push_with_raw(command, None);
    }

    /// Add a command with optional raw string to the run history
    ///
    /// If raw_string is provided, it will be stored alongside the command
    /// for potential later serialization as a string.
    pub fn push_with_raw(&mut self, command: Commands, raw_string: Option<String>) {
        if should_persist_command(&command) {
            self.commands.push(CommandHistory {
                command,
                raw_string,
            });
        }
    }

    pub fn schema() -> Schema {
        schema_for!(RunHistory)
    }

    pub fn validate(&self) -> Result<()> {
        let mut seen_names = HashSet::with_capacity(self.command_blocks.len());
        for block in &self.command_blocks {
            if block.name.trim().is_empty() {
                return Err(eyre!(
                    "Run card `command_blocks` contains a block with an empty name"
                ));
            }
            if !seen_names.insert(block.name.clone()) {
                return Err(eyre!(
                    "Run card `command_blocks` contains duplicate block name '{}'",
                    block.name
                ));
            }
        }
        Ok(())
    }

    pub fn command_block(&self, name: &str) -> Option<&CommandsBlock> {
        self.command_blocks.iter().find(|block| block.name == name)
    }

    pub fn select_command_blocks(
        &self,
        selected_block_names: &[String],
    ) -> Result<Vec<CommandsBlock>> {
        let mut selected = Vec::with_capacity(selected_block_names.len());
        for name in selected_block_names {
            let block = self.command_block(name).ok_or_else(|| {
                eyre!(
                    "Unknown command block '{}'. Available command blocks: {}",
                    name,
                    self.command_blocks
                        .iter()
                        .map(|block| block.name.as_str())
                        .collect::<Vec<_>>()
                        .join(", ")
                )
            })?;
            selected.push(block.clone());
        }
        Ok(selected)
    }

    pub fn conflicting_command_block_names(&self, command_blocks: &[CommandsBlock]) -> Vec<String> {
        command_blocks
            .iter()
            .filter_map(|new_block| match self.command_block(&new_block.name) {
                Some(existing_block) if !existing_block.semantically_eq(new_block) => {
                    Some(new_block.name.clone())
                }
                _ => None,
            })
            .collect()
    }

    pub fn merge_command_blocks_with_overwrite(
        &mut self,
        command_blocks: &[CommandsBlock],
        overwrite_conflicts: bool,
    ) -> Result<()> {
        for new_block in command_blocks {
            match self
                .command_blocks
                .iter()
                .position(|existing| existing.name == new_block.name)
            {
                Some(index) if self.command_blocks[index].semantically_eq(new_block) => {}
                Some(index) if overwrite_conflicts => {
                    self.command_blocks[index] = new_block.clone();
                }
                Some(_) => {
                    return Err(eyre!(
                        "Run card command block '{}' redefines an existing block with different commands",
                        new_block.name
                    ));
                }
                None => self.command_blocks.push(new_block.clone()),
            }
        }
        self.validate()
    }

    pub fn merge_command_blocks(&mut self, command_blocks: &[CommandsBlock]) -> Result<()> {
        self.merge_command_blocks_with_overwrite(command_blocks, false)
    }

    pub fn apply_session_settings(
        &self,
        global_settings: &mut CLISettings,
        default_runtime_settings: &mut RuntimeSettings,
    ) -> Result<()> {
        if self.cli_settings.global != GlobalSettings::default() {
            global_settings.global = self.cli_settings.global.clone();
            global_settings.sync_settings()?;
        }
        if self.default_runtime_settings != RuntimeSettings::default() {
            *default_runtime_settings = self.default_runtime_settings.clone();
        }
        Ok(())
    }

    pub fn run(
        &mut self,
        state: &mut State,
        global_settings: &mut CLISettings,
        default_runtime_settings: &mut RuntimeSettings,
    ) -> Result<ControlFlow<SaveState>> {
        let mut session_state = crate::session::CliSessionState::default();
        let mut session = crate::session::CliSession::new(
            state,
            self,
            global_settings,
            default_runtime_settings,
            &mut session_state,
        );
        session.replay_run_history()
    }

    pub(crate) fn filtered_for_save(&self) -> Self {
        let mut filtered = self.clone();
        filtered
            .commands
            .retain(|command_history| should_persist_command(&command_history.command));
        filtered
    }

    pub(crate) fn to_toml_string(&self, serialize_commands_as_strings: bool) -> Result<String> {
        let _serialize_commands_guard =
            SerializeCommandsAsStringsGuard::new(serialize_commands_as_strings);
        render_smart_toml(&self.filtered_for_save())
    }

    pub fn load(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        debug!("Loaded run history from file {}", path.display());

        let runhistory = match path.extension().and_then(|ext| ext.to_str()) {
            Some("toml") => Self::load_toml(path)?,
            _ => Self::from_file(path, "run history")?,
        };
        runhistory.validate()?;
        Ok(runhistory)
    }

    pub fn save_toml(
        &self,
        root_folder: &Path,
        override_state_file: bool,
        _strict: bool,
    ) -> Result<()> {
        self.filtered_for_save()
            .to_file(root_folder.join("run.toml"), override_state_file)?;

        //Self::schema().to_file(root_folder.join("run_schema.json"))?;
        Ok(())
    }

    pub fn save_yaml(
        &self,
        root_folder: &Path,
        override_state_file: bool,
        _strict: bool,
    ) -> Result<()> {
        self.filtered_for_save()
            .to_file(root_folder.join("run.yaml"), override_state_file)?;
        //Self::schema().to_file(root_folder.join("run_schema.json"))?;
        Ok(())
    }

    fn load_toml(path: &Path) -> Result<Self> {
        let raw = fs::read_to_string(path)
            .with_context(|| format!("Could not open run history file {}", path.display()))?;
        let raw_run_history: RawRunHistoryToml =
            toml::from_str(&raw).wrap_err("Error parsing run history toml")?;

        let commands = raw_run_history
            .commands
            .into_iter()
            .enumerate()
            .map(|(index, value)| {
                parse_toml_command_history(value, &format!("top-level command #{}", index + 1))
            })
            .collect::<Result<Vec<_>>>()?;

        let command_blocks = raw_run_history
            .command_blocks
            .into_iter()
            .map(|block| {
                let RawCommandsBlockToml { name, commands } = block;
                let commands = commands
                    .into_iter()
                    .enumerate()
                    .map(|(index, value)| {
                        parse_toml_command_history(
                            value,
                            &format!("command block '{}' command #{}", name, index + 1),
                        )
                    })
                    .collect::<Result<Vec<_>>>()?;
                Ok(CommandsBlock { name, commands })
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Self {
            default_runtime_settings: raw_run_history.default_runtime_settings,
            cli_settings: raw_run_history.cli_settings,
            command_blocks,
            commands,
        })
    }
}

fn parse_toml_command_history(value: TomlValue, context: &str) -> Result<CommandHistory> {
    match value {
        TomlValue::String(raw) => CommandHistory::from_raw_string(&raw)
            .map_err(|err| eyre!("Failed to parse {} '{}': {}", context, raw, err)),
        other => other
            .try_into::<CommandHistory>()
            .map_err(|err| eyre!("Failed to parse {}: {}", context, err)),
    }
}

#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(from_py_object, unsendable, name = "GammaLoopState")
)]
#[derive(Clone)]
pub struct State {
    pub model: Model,
    pub model_parameters: InputParamCard<F<f64>>,
    pub process_list: ProcessList,
    pub generation_summaries: BTreeMap<IntegrandGenerationSummaryKey, IntegrandGenerationSummary>,
}

const STATE_MANIFEST_FILE: &str = "state_manifest.toml";
const INTEGRAND_GENERATION_SUMMARY_FILE: &str = "generation_summary.json";
const CURRENT_STATE_MANIFEST_VERSION: u32 = 1;
const GENERATION_THREAD_STACK_SIZE_BYTES: usize = 32 * 1024 * 1024;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(default, deny_unknown_fields)]
struct StateManifest {
    version: u32,
}

impl Default for StateManifest {
    fn default() -> Self {
        Self {
            version: CURRENT_STATE_MANIFEST_VERSION,
        }
    }
}

fn ensure_supported_state_manifest_version(manifest: &StateManifest) -> Result<()> {
    if manifest.version > CURRENT_STATE_MANIFEST_VERSION {
        return Err(eyre!(
            "State version {} is newer than this binary supports (max {}). Please upgrade gammaloop.",
            manifest.version,
            CURRENT_STATE_MANIFEST_VERSION
        ));
    }

    Ok(())
}

fn run_state_migration_checks(manifest: &StateManifest, save_path: &Path) -> Result<()> {
    ensure_supported_state_manifest_version(manifest)?;

    match manifest.version {
        1 => {
            if !save_path.join("symbolica_state.bin").exists() {
                return Err(eyre!(
                    "Saved state at '{}' is missing required file symbolica_state.bin",
                    save_path.display()
                ));
            }
            if !save_path.join("processes").exists() {
                return Err(eyre!(
                    "Saved state at '{}' is missing required folder processes/",
                    save_path.display()
                ));
            }
            if !save_path.join("model.json").exists() {
                return Err(eyre!(
                    "Saved state at '{}' is missing required file model.json",
                    save_path.display()
                ));
            }
            Ok(())
        }
        _ => Err(eyre!(
            "State version {} is not supported by this binary.",
            manifest.version
        )),
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StateFolderKind {
    Missing,
    Scratch,
    Unmanifested,
    Saved,
    Invalid(String),
}

fn is_scratch_state_entry(entry: &fs::DirEntry) -> bool {
    entry.file_type().map(|ft| ft.is_dir()).unwrap_or(false)
        && entry.file_name().to_string_lossy() == "logs"
}

pub fn classify_state_folder(save_path: &Path) -> Result<StateFolderKind> {
    if !save_path.exists() {
        return Ok(StateFolderKind::Missing);
    }
    if !save_path.is_dir() {
        return Ok(StateFolderKind::Invalid(format!(
            "'{}' exists but is not a directory",
            save_path.display()
        )));
    }

    let manifest_path = save_path.join(STATE_MANIFEST_FILE);
    if manifest_path.exists() {
        let manifest = load_state_manifest(save_path)?;
        return Ok(match run_state_migration_checks(&manifest, save_path) {
            Ok(()) => StateFolderKind::Saved,
            Err(err) => StateFolderKind::Invalid(err.to_string()),
        });
    }

    let mut entries = fs::read_dir(save_path)
        .with_context(|| format!("Trying to read state folder '{}'", save_path.display()))?;
    if entries.by_ref().all(|entry| {
        entry
            .map(|entry| is_scratch_state_entry(&entry))
            .unwrap_or(false)
    }) {
        return Ok(StateFolderKind::Scratch);
    }

    Ok(StateFolderKind::Unmanifested)
}

fn load_state_manifest(save_path: &Path) -> Result<StateManifest> {
    let manifest_path = save_path.join(STATE_MANIFEST_FILE);
    let raw_manifest = fs::read_to_string(&manifest_path).with_context(|| {
        format!(
            "Trying to read state manifest file {}",
            manifest_path.display()
        )
    })?;
    let manifest = toml::from_str::<StateManifest>(&raw_manifest).with_context(|| {
        format!(
            "Trying to parse state manifest file {}",
            manifest_path.display()
        )
    })?;
    ensure_supported_state_manifest_version(&manifest)?;
    Ok(manifest)
}

fn save_state_manifest(save_path: &Path) -> Result<()> {
    let manifest = StateManifest::default();
    let raw_manifest =
        toml::to_string_pretty(&manifest).context("Trying to serialize state manifest to TOML")?;
    fs::write(save_path.join(STATE_MANIFEST_FILE), raw_manifest).with_context(|| {
        format!(
            "Trying to write state manifest file {}",
            save_path.join(STATE_MANIFEST_FILE).display()
        )
    })?;
    Ok(())
}

fn integrand_generation_summary_path(
    root_folder: &Path,
    process: &Process,
    integrand_name: &str,
) -> PathBuf {
    let process_kind_folder = match &process.collection {
        ProcessCollection::Amplitudes(_) => "amplitudes",
        ProcessCollection::CrossSections(_) => "cross_sections",
    };

    root_folder
        .join("processes")
        .join(process_kind_folder)
        .join(&process.definition.folder_name)
        .join(integrand_name)
        .join(INTEGRAND_GENERATION_SUMMARY_FILE)
}

fn process_kind_folder(root_folder: &Path, process: &Process) -> PathBuf {
    let process_kind_folder = match &process.collection {
        ProcessCollection::Amplitudes(_) => "amplitudes",
        ProcessCollection::CrossSections(_) => "cross_sections",
    };

    root_folder.join("processes").join(process_kind_folder)
}

fn process_artifact_folder(root_folder: &Path, process: &Process) -> PathBuf {
    process_kind_folder(root_folder, process).join(&process.definition.folder_name)
}

fn integrand_artifact_folder(
    root_folder: &Path,
    process: &Process,
    integrand_name: &str,
) -> PathBuf {
    process_artifact_folder(root_folder, process).join(integrand_name)
}

fn generated_integrand_artifact_path(
    root_folder: &Path,
    process: &Process,
    integrand_name: &str,
) -> PathBuf {
    integrand_artifact_folder(root_folder, process, integrand_name).join("integrand")
}

fn ensure_existing_path_under_root(path: &Path, root: &Path, description: &str) -> Result<PathBuf> {
    let canonical_root = root
        .canonicalize()
        .with_context(|| format!("Trying to resolve artifact cleanup root {}", root.display()))?;
    let canonical_path = path.canonicalize().with_context(|| {
        format!(
            "Trying to resolve {} path {} before removal",
            description,
            path.display()
        )
    })?;
    if !canonical_path.starts_with(&canonical_root) {
        return Err(eyre!(
            "Refusing to remove {} path {} because it resolves outside {}",
            description,
            canonical_path.display(),
            canonical_root.display()
        ));
    }
    Ok(canonical_path)
}

fn remove_file_if_exists_under_root(path: &Path, root: &Path, description: &str) -> Result<bool> {
    if !path.try_exists().with_context(|| {
        format!(
            "Trying to check whether {} path {} exists",
            description,
            path.display()
        )
    })? {
        return Ok(false);
    }
    let canonical_path = ensure_existing_path_under_root(path, root, description)?;
    fs::remove_file(path).with_context(|| {
        format!(
            "Trying to remove {} path {}",
            description,
            canonical_path.display()
        )
    })?;
    Ok(true)
}

fn remove_dir_if_exists_under_root(path: &Path, root: &Path, description: &str) -> Result<bool> {
    if !path.try_exists().with_context(|| {
        format!(
            "Trying to check whether {} path {} exists",
            description,
            path.display()
        )
    })? {
        return Ok(false);
    }
    let canonical_path = ensure_existing_path_under_root(path, root, description)?;
    fs::remove_dir_all(path).with_context(|| {
        format!(
            "Trying to remove {} path {}",
            description,
            canonical_path.display()
        )
    })?;
    Ok(true)
}

fn remove_saved_process_artifacts(root_folder: &Path, process: &Process) -> Result<bool> {
    let root = process_kind_folder(root_folder, process);
    let path = process_artifact_folder(root_folder, process);
    remove_dir_if_exists_under_root(&path, &root, "saved process artifacts")
}

fn remove_saved_integrand_artifacts(
    root_folder: &Path,
    process: &Process,
    integrand_name: &str,
) -> Result<bool> {
    let root = process_artifact_folder(root_folder, process);
    let path = integrand_artifact_folder(root_folder, process, integrand_name);
    remove_dir_if_exists_under_root(&path, &root, "saved integrand artifacts")
}

fn validate_output_name(value: &str, flag_name: &str) -> Result<()> {
    if value.trim().is_empty() {
        return Err(eyre!("{flag_name} must not be empty"));
    }
    Ok(())
}

fn rename_process_integrand(integrand: Option<&mut ProcessIntegrand>, new_name: &str) {
    let Some(integrand) = integrand else {
        return;
    };
    let _ = integrand.get_mut_settings();
    match integrand {
        ProcessIntegrand::Amplitude(amplitude) => {
            amplitude.data.name = new_name.to_string();
        }
        ProcessIntegrand::CrossSection(cross_section) => {
            cross_section.data.name = new_name.to_string();
        }
    }
}

fn load_integrand_generation_summaries(
    root_folder: &Path,
    process_list: &ProcessList,
) -> Result<BTreeMap<IntegrandGenerationSummaryKey, IntegrandGenerationSummary>> {
    let mut summaries = BTreeMap::new();

    for (process_id, process) in process_list.processes.iter().enumerate() {
        for integrand_name in process.collection.get_integrand_names() {
            let summary_path =
                integrand_generation_summary_path(root_folder, process, integrand_name);
            if !summary_path.exists() {
                continue;
            }

            let raw_summary = fs::read_to_string(&summary_path).with_context(|| {
                format!(
                    "Trying to read integrand generation summary {}",
                    summary_path.display()
                )
            })?;
            let summary = serde_json::from_str(&raw_summary).with_context(|| {
                format!(
                    "Trying to parse integrand generation summary {}",
                    summary_path.display()
                )
            })?;
            summaries.insert(
                IntegrandGenerationSummaryKey {
                    process_id,
                    integrand_name: integrand_name.to_string(),
                },
                summary,
            );
        }
    }

    Ok(summaries)
}

impl State {
    fn overridable_model_parameter_names(&self) -> Vec<String> {
        let mut names = self
            .model_parameters
            .keys()
            .map(|symbol| symbol.to_string())
            .collect::<Vec<_>>();
        names.sort();
        names
    }

    fn resolve_overridable_model_parameter_type(
        &self,
        parameter_name: &str,
    ) -> Result<gammalooprs::model::ParameterType> {
        let possibilities = self.overridable_model_parameter_names();
        let parameter_type = external_model_parameter_type(&self.model, parameter_name)
            .ok_or_else(|| eyre!("No model parameter named '{parameter_name}'"))
            .with_note(|| {
                format!(
                    "Possible model parameters are: {}",
                    possibilities.join(", ")
                )
            })?;

        let symbol = UFOSymbol::from(parameter_name);
        if !self.model_parameters.contains_key(&symbol) {
            return Err(eyre!(
                "Model parameter '{parameter_name}' cannot be overridden because it is not present in the shared top-level model_parameters.json"
            ))
            .with_note(|| format!("Possible model parameters are: {}", possibilities.join(", ")));
        }

        Ok(parameter_type)
    }

    pub fn find_generated_integrand_ref_by_name(
        &self,
        integrand_name: &str,
    ) -> Result<(usize, String)> {
        let matches = self
            .process_list
            .processes
            .iter()
            .enumerate()
            .filter_map(|(process_id, process)| {
                let found = match &process.collection {
                    ProcessCollection::Amplitudes(amplitudes) => amplitudes
                        .get(integrand_name)
                        .filter(|amplitude| amplitude.integrand.is_some())
                        .map(|_| (process_id, process.definition.folder_name.clone())),
                    ProcessCollection::CrossSections(cross_sections) => cross_sections
                        .get(integrand_name)
                        .filter(|cross_section| cross_section.integrand.is_some())
                        .map(|_| (process_id, process.definition.folder_name.clone())),
                };
                found.map(|(process_id, process_name)| {
                    (process_id, process_name, integrand_name.to_string())
                })
            })
            .collect::<Vec<_>>();

        match matches.as_slice() {
            [(process_id, _, canonical_name)] => Ok((*process_id, canonical_name.clone())),
            [] => Err(eyre!(
                "No generated integrand named '{integrand_name}' was found. Per-integrand model parameters are only supported for generated integrands."
            )),
            _ => {
                let names = matches
                    .iter()
                    .map(|(process_id, process_name, _)| format!("#{process_id} ({process_name})"))
                    .collect::<Vec<_>>();
                Err(eyre!(
                    "Integrand name '{integrand_name}' is ambiguous across generated integrands"
                ))
                .with_note(|| format!("Matching processes: {}", names.join(", ")))
            }
        }
    }

    pub fn import_model(&mut self, path: impl AsRef<Path>) -> Result<()> {
        self.model = Model::from_file(path)?;
        Ok(())
    }

    fn matching_process_ids(&self, process: Option<&ProcessRef>) -> Result<Vec<usize>> {
        let Some(process) = process else {
            return Ok((0..self.process_list.processes.len()).collect());
        };

        match process {
            ProcessRef::Id(id) => Ok((*id < self.process_list.processes.len())
                .then_some(*id)
                .into_iter()
                .collect()),
            ProcessRef::Name(name) => Ok(self
                .process_list
                .processes
                .iter()
                .position(|p| p.definition.folder_name == *name)
                .into_iter()
                .collect()),
            ProcessRef::Unqualified(value) => {
                let name_match = self
                    .process_list
                    .processes
                    .iter()
                    .position(|p| p.definition.folder_name == *value);
                if let Ok(id) = value.parse::<usize>() {
                    let id_valid = id < self.process_list.processes.len();
                    match (id_valid, name_match) {
                        (true, Some(_)) => Err(eyre!(
                            "Ambiguous process reference '{}'. Use '#{}' or 'name:{}' to disambiguate.",
                            value,
                            id,
                            value
                        )),
                        (true, None) => Ok(vec![id]),
                        (false, Some(index)) => Ok(vec![index]),
                        (false, None) => Ok(vec![]),
                    }
                } else {
                    Ok(name_match.into_iter().collect())
                }
            }
        }
    }

    pub fn remove_selected_integrands(
        &mut self,
        process: Option<&ProcessRef>,
        integrand_name: Option<&str>,
    ) -> Result<Vec<RemovedIntegrand>> {
        let process_ids = self.matching_process_ids(process)?;
        let mut removed = Vec::new();

        for process_id in process_ids.into_iter().rev() {
            let process_name = self.process_list.processes[process_id]
                .definition
                .folder_name
                .clone();
            let integrand_names = {
                let collection = &self.process_list.processes[process_id].collection;
                match integrand_name {
                    Some(name) => collection
                        .get_integrand_names()
                        .into_iter()
                        .find(|candidate| *candidate == name)
                        .map(|name| vec![name.to_string()])
                        .unwrap_or_default(),
                    None => collection
                        .get_integrand_names()
                        .into_iter()
                        .map(str::to_string)
                        .collect(),
                }
            };

            if integrand_names.is_empty() {
                continue;
            }

            {
                let process_entry = &mut self.process_list.processes[process_id];
                for integrand_name in &integrand_names {
                    process_entry.collection.remove_integrand(integrand_name)?;
                }
            }
            for integrand_name in &integrand_names {
                self.generation_summaries
                    .remove(&IntegrandGenerationSummaryKey {
                        process_id,
                        integrand_name: integrand_name.clone(),
                    });
            }

            let removed_empty_process = self.process_list.processes[process_id]
                .collection
                .get_integrand_names()
                .is_empty();
            if removed_empty_process {
                self.process_list.processes.remove(process_id);
                self.remove_generation_summaries_for_process(process_id);
            }

            for integrand_name in integrand_names {
                removed.push(RemovedIntegrand {
                    process_id,
                    process_name: process_name.clone(),
                    integrand_name,
                    removed_empty_process,
                });
            }
        }

        removed.reverse();
        Ok(removed)
    }

    pub fn remove_process(&mut self, process: Option<&ProcessRef>) -> Result<RemovedProcess> {
        let process_id = self.resolve_process_ref(process)?;
        let removed = self.process_list.processes.remove(process_id);
        self.remove_generation_summaries_for_process(process_id);
        Ok(RemovedProcess {
            process_id,
            process_name: removed.definition.folder_name,
        })
    }

    pub fn remove_integrand(
        &mut self,
        process: &ProcessRef,
        integrand_name: &str,
    ) -> Result<RemovedIntegrand> {
        let process_id = process.resolve(&self.process_list)?;
        let (process_name, canonical_integrand_name, removed_empty_process) = {
            let process_entry = &mut self.process_list.processes[process_id];
            let process_name = process_entry.definition.folder_name.clone();
            let canonical_integrand_name = process_entry
                .collection
                .find_integrand(Some(integrand_name.to_string()))?;
            process_entry
                .collection
                .remove_integrand(&canonical_integrand_name)?;
            let removed_empty_process = process_entry.collection.get_integrand_names().is_empty();
            (
                process_name,
                canonical_integrand_name,
                removed_empty_process,
            )
        };
        self.generation_summaries
            .remove(&IntegrandGenerationSummaryKey {
                process_id,
                integrand_name: canonical_integrand_name.clone(),
            });

        if removed_empty_process {
            self.process_list.processes.remove(process_id);
            self.remove_generation_summaries_for_process(process_id);
        }

        Ok(RemovedIntegrand {
            process_id,
            process_name,
            integrand_name: canonical_integrand_name,
            removed_empty_process,
        })
    }

    pub fn resolve_process_ref(&self, process: Option<&ProcessRef>) -> Result<usize> {
        match process {
            Some(process_ref) => process_ref.resolve(&self.process_list),
            None => self.process_list.find_process(None),
        }
    }

    pub fn find_integrand_ref(
        &self,
        process: Option<&ProcessRef>,
        integrand_name: Option<&String>,
    ) -> Result<(usize, String)> {
        let process_id = self.resolve_process_ref(process)?;
        let integrand_name = self.process_list.processes[process_id]
            .collection
            .find_integrand(integrand_name.cloned())
            .with_note(|| format!("in process id {process_id}"))?;
        Ok((process_id, integrand_name))
    }

    pub fn get_integrand_info(
        &self,
        process: Option<&ProcessRef>,
        integrand_name: Option<&String>,
    ) -> Result<IntegrandInfo> {
        let (process_id, integrand_name) = self.find_integrand_ref(process, integrand_name)?;
        collect_integrand_info(self, process_id, &integrand_name)
    }

    pub fn duplicate_integrand(
        &mut self,
        process: Option<&ProcessRef>,
        integrand_name: Option<&String>,
        output_process_name: &str,
        output_integrand_name: &str,
    ) -> Result<()> {
        validate_output_name(output_process_name, "--output_process_name")?;
        validate_output_name(output_integrand_name, "--output_integrand_name")?;

        let (source_process_id, source_integrand_name) =
            self.find_integrand_ref(process, integrand_name)?;
        let mut payload =
            self.cloned_integrand_payload(source_process_id, &source_integrand_name)?;
        payload.rename(output_integrand_name);
        self.insert_integrand_copy(
            source_process_id,
            output_process_name,
            output_integrand_name,
            payload,
            false,
            false,
            None,
            false,
        )?;
        Ok(())
    }

    pub fn select_integrand_graph_groups(
        &mut self,
        process: Option<&ProcessRef>,
        integrand_name: Option<&String>,
        selection: &GraphGroupSelectionSpec,
        target: &GraphGroupSelectionTarget,
        context: GraphGroupSelectionContext<'_>,
    ) -> Result<SelectedGraphGroups> {
        if target.clear_existing_processes && !target.is_copy_mode() {
            return Err(eyre!(
                "--clear-existing-processes requires --output_process or --output_integrand for `select`"
            ));
        }

        let (source_process_id, source_integrand_name) =
            self.find_integrand_ref(process, integrand_name)?;
        let (source_process_name, plan, discarded_generated_integrand) = {
            let process_entry = &self.process_list.processes[source_process_id];
            let process_name = process_entry.definition.folder_name.clone();
            match &process_entry.collection {
                ProcessCollection::Amplitudes(amplitudes) => {
                    if selection.has_raised_cut_rules() {
                        return Err(eyre!(
                            "Raised-cut signature selection can only be used with cross-section integrands."
                        ));
                    }
                    if selection.mode() == GraphGroupSelectionMode::CrossSectionAmplitudeGraphs {
                        return Err(eyre!(
                            "`select --amplitude-graphs` can only be used with cross-section integrands."
                        ));
                    }
                    let amplitude = amplitudes.get(&source_integrand_name).ok_or_else(|| {
                        eyre!(
                            "No amplitude named '{}' in process '{}'",
                            source_integrand_name,
                            process_name
                        )
                    })?;
                    let plan = amplitude.plan_graph_group_selection(selection)?;
                    amplitude.validate_graph_group_selection_plan(&plan)?;
                    (process_name, plan, amplitude.integrand.is_some())
                }
                ProcessCollection::CrossSections(cross_sections) => {
                    let cross_section =
                        cross_sections.get(&source_integrand_name).ok_or_else(|| {
                            eyre!(
                                "No cross section named '{}' in process '{}'",
                                source_integrand_name,
                                process_name
                            )
                        })?;
                    let plan = cross_section.plan_graph_group_selection_with_context(
                        selection,
                        &self.model,
                        &process_entry.definition,
                        context.generation_settings,
                    )?;
                    cross_section.validate_graph_group_selection_plan(&plan)?;
                    (process_name, plan, cross_section.integrand.is_some())
                }
            }
        };
        let report = plan.report().clone();
        let source = SelectionSource {
            process_id: source_process_id,
            process_name: source_process_name,
            integrand_name: source_integrand_name,
        };

        if target.is_copy_mode() {
            return self.select_integrand_graph_groups_to_output(
                &source,
                &plan,
                report,
                target,
                context.state_folder,
                context.read_only_state,
            );
        }

        let mut removed_generated_artifacts = false;
        let mut removed_generation_summary = false;
        if discarded_generated_integrand {
            if context.read_only_state {
                return Err(eyre!(
                    "Cannot select graph groups for generated integrand '{}' in process '{}' because this session was started with --read-only-state. Restart without --read-only-state or select before generation.",
                    source.integrand_name,
                    source.process_name
                ));
            }

            let process_entry = &self.process_list.processes[source.process_id];
            let integrand_artifact_root = integrand_artifact_folder(
                context.state_folder,
                process_entry,
                &source.integrand_name,
            );
            let generated_artifact_path = generated_integrand_artifact_path(
                context.state_folder,
                process_entry,
                &source.integrand_name,
            );
            removed_generated_artifacts = remove_dir_if_exists_under_root(
                &generated_artifact_path,
                &integrand_artifact_root,
                "generated integrand artifact",
            )?;
            let generation_summary_path = integrand_generation_summary_path(
                context.state_folder,
                process_entry,
                &source.integrand_name,
            );
            removed_generation_summary = remove_file_if_exists_under_root(
                &generation_summary_path,
                &integrand_artifact_root,
                "integrand generation summary",
            )?;
        }

        {
            let process_entry = &mut self.process_list.processes[source.process_id];
            match &mut process_entry.collection {
                ProcessCollection::Amplitudes(amplitudes) => {
                    let amplitude =
                        amplitudes.get_mut(&source.integrand_name).ok_or_else(|| {
                            eyre!(
                                "No amplitude named '{}' in process '{}'",
                                source.integrand_name,
                                source.process_name
                            )
                        })?;
                    amplitude.apply_graph_group_selection(&plan)?;
                }
                ProcessCollection::CrossSections(cross_sections) => {
                    let cross_section =
                        cross_sections
                            .get_mut(&source.integrand_name)
                            .ok_or_else(|| {
                                eyre!(
                                    "No cross section named '{}' in process '{}'",
                                    source.integrand_name,
                                    source.process_name
                                )
                            })?;
                    cross_section.apply_graph_group_selection(&plan)?;
                }
            }
        }

        if discarded_generated_integrand {
            self.generation_summaries
                .remove(&IntegrandGenerationSummaryKey {
                    process_id: source.process_id,
                    integrand_name: source.integrand_name.clone(),
                });
        }

        Ok(SelectedGraphGroups {
            source_process_id: source.process_id,
            source_process_name: source.process_name.clone(),
            source_integrand_name: source.integrand_name.clone(),
            process_id: source.process_id,
            process_name: source.process_name,
            integrand_name: source.integrand_name,
            report,
            copied_to_output: false,
            replaced_existing_target: false,
            removed_target_artifacts: false,
            discarded_generated_integrand,
            removed_generated_artifacts,
            removed_generation_summary,
        })
    }

    fn select_integrand_graph_groups_to_output(
        &mut self,
        source: &SelectionSource,
        plan: &GraphGroupSelectionPlan,
        report: GraphGroupSelectionReport,
        target: &GraphGroupSelectionTarget,
        state_folder: &Path,
        read_only_state: bool,
    ) -> Result<SelectedGraphGroups> {
        let output_process_name = target
            .output_process_name
            .as_deref()
            .unwrap_or(&source.process_name);
        let output_integrand_name = target
            .output_integrand_name
            .as_deref()
            .unwrap_or(&source.integrand_name);

        validate_output_name(output_process_name, "--output_process")?;
        validate_output_name(output_integrand_name, "--output_integrand")?;

        if output_process_name == source.process_name
            && output_integrand_name == source.integrand_name
        {
            return Err(eyre!(
                "Copy-mode select target '{} / {}' is the selected source integrand. Omit --output_process/--output_integrand for in-place selection or choose a different output target.",
                output_process_name,
                output_integrand_name
            ));
        }

        let target_process_id = self
            .process_list
            .processes
            .iter()
            .position(|process| process.definition.folder_name == output_process_name);
        let replace_existing_process = target.output_process_name.is_some()
            && target.clear_existing_processes
            && target_process_id.is_some_and(|process_id| process_id != source.process_id);

        let mut payload =
            self.cloned_integrand_payload(source.process_id, &source.integrand_name)?;
        payload.rename(output_integrand_name);
        payload.apply_graph_group_selection(plan)?;

        let insertion = self.insert_integrand_copy(
            source.process_id,
            output_process_name,
            output_integrand_name,
            payload,
            target.clear_existing_processes,
            replace_existing_process,
            Some(state_folder),
            read_only_state,
        )?;

        Ok(SelectedGraphGroups {
            source_process_id: source.process_id,
            source_process_name: source.process_name.clone(),
            source_integrand_name: source.integrand_name.clone(),
            process_id: insertion.process_id,
            process_name: insertion.process_name,
            integrand_name: insertion.integrand_name,
            report,
            copied_to_output: true,
            replaced_existing_target: insertion.replaced_existing_target,
            removed_target_artifacts: insertion.removed_target_artifacts,
            discarded_generated_integrand: false,
            removed_generated_artifacts: false,
            removed_generation_summary: false,
        })
    }

    fn cloned_integrand_payload(
        &self,
        source_process_id: usize,
        source_integrand_name: &str,
    ) -> Result<IntegrandCopyPayload> {
        let source_process = &self.process_list.processes[source_process_id];
        match &source_process.collection {
            ProcessCollection::Amplitudes(amplitudes) => amplitudes
                .get(source_integrand_name)
                .cloned()
                .map(IntegrandCopyPayload::Amplitude)
                .ok_or_else(|| eyre!("Missing source amplitude '{}'", source_integrand_name)),
            ProcessCollection::CrossSections(cross_sections) => cross_sections
                .get(source_integrand_name)
                .cloned()
                .map(IntegrandCopyPayload::CrossSection)
                .ok_or_else(|| eyre!("Missing source cross section '{}'", source_integrand_name)),
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn insert_integrand_copy(
        &mut self,
        source_process_id: usize,
        output_process_name: &str,
        output_integrand_name: &str,
        payload: IntegrandCopyPayload,
        clear_existing_processes: bool,
        replace_existing_process: bool,
        state_folder: Option<&Path>,
        read_only_state: bool,
    ) -> Result<IntegrandCopyInsertion> {
        if let Some(destination_process_id) = self
            .process_list
            .processes
            .iter()
            .position(|process| process.definition.folder_name == output_process_name)
        {
            if replace_existing_process {
                if read_only_state {
                    return Err(eyre!(
                        "Cannot overwrite output process '{}' because this session was started with --read-only-state.",
                        output_process_name
                    ));
                }
                let removed_target_artifacts = if let Some(state_folder) = state_folder {
                    remove_saved_process_artifacts(
                        state_folder,
                        &self.process_list.processes[destination_process_id],
                    )?
                } else {
                    false
                };
                self.remove_generation_summaries_for_process_without_shifting(
                    destination_process_id,
                );
                let process = self.single_integrand_process(
                    source_process_id,
                    destination_process_id,
                    output_process_name,
                    payload,
                );
                self.process_list.processes[destination_process_id] = process;
                return Ok(IntegrandCopyInsertion {
                    process_id: destination_process_id,
                    process_name: output_process_name.to_string(),
                    integrand_name: output_integrand_name.to_string(),
                    replaced_existing_target: true,
                    removed_target_artifacts,
                });
            }

            let target_exists = self.process_list.processes[destination_process_id]
                .collection
                .get_integrand_names()
                .contains(&output_integrand_name);
            if target_exists {
                if !clear_existing_processes {
                    return Err(eyre!(
                        "An integrand '{}' already exists in process '{}'",
                        output_integrand_name,
                        output_process_name
                    ));
                }
                if read_only_state {
                    return Err(eyre!(
                        "Cannot overwrite output integrand '{}' in process '{}' because this session was started with --read-only-state.",
                        output_integrand_name,
                        output_process_name
                    ));
                }
            }

            if !payload
                .is_compatible_with(&self.process_list.processes[destination_process_id].collection)
            {
                return Err(eyre!(
                    "Destination process '{}' exists but does not contain {}",
                    output_process_name,
                    payload.kind_name()
                ));
            }

            let removed_target_artifacts = if target_exists {
                if let Some(state_folder) = state_folder {
                    remove_saved_integrand_artifacts(
                        state_folder,
                        &self.process_list.processes[destination_process_id],
                        output_integrand_name,
                    )?
                } else {
                    false
                }
            } else {
                false
            };

            if target_exists {
                self.generation_summaries
                    .remove(&IntegrandGenerationSummaryKey {
                        process_id: destination_process_id,
                        integrand_name: output_integrand_name.to_string(),
                    });
            }

            payload.insert_into_process(
                &mut self.process_list.processes[destination_process_id],
                output_process_name,
            )?;
            Ok(IntegrandCopyInsertion {
                process_id: destination_process_id,
                process_name: output_process_name.to_string(),
                integrand_name: output_integrand_name.to_string(),
                replaced_existing_target: target_exists,
                removed_target_artifacts,
            })
        } else {
            let process_id = self.process_list.processes.len();
            let process = self.single_integrand_process(
                source_process_id,
                process_id,
                output_process_name,
                payload,
            );
            self.process_list.add_process(process);
            Ok(IntegrandCopyInsertion {
                process_id,
                process_name: output_process_name.to_string(),
                integrand_name: output_integrand_name.to_string(),
                replaced_existing_target: false,
                removed_target_artifacts: false,
            })
        }
    }

    fn single_integrand_process(
        &self,
        source_process_id: usize,
        process_id: usize,
        process_name: &str,
        payload: IntegrandCopyPayload,
    ) -> Process {
        let source_process = &self.process_list.processes[source_process_id];
        let mut definition = source_process.definition.clone();
        definition.folder_name = process_name.to_string();
        definition.process_id = process_id;
        Process {
            definition,
            settings_history: source_process.settings_history.clone(),
            collection: payload.into_collection(),
        }
    }

    pub fn resolve_effective_model_parameter_card_for_settings(
        &self,
        settings: &RuntimeSettings,
    ) -> Result<InputParamCard<F<f64>>> {
        let mut card = self.model_parameters.clone();

        for (parameter_name, value) in &settings.model.external_parameters {
            let parameter_type = self.resolve_overridable_model_parameter_type(parameter_name)?;
            let value = Complex::new(value.0, value.1);
            validate_model_parameter_type(parameter_name, parameter_type, &value)?;

            let parameter = card
                .get_mut(&UFOSymbol::from(parameter_name.as_str()))
                .ok_or_else(|| {
                    eyre!(
                        "Model parameter '{parameter_name}' is missing from the shared top-level model_parameters.json"
                    )
                })?;
            *parameter = value;
        }

        Ok(card)
    }

    pub fn resolve_serializable_model_parameter_card_for_settings(
        &self,
        settings: &RuntimeSettings,
    ) -> Result<SerializableInputParamCard<F<f64>>> {
        Ok(self
            .resolve_effective_model_parameter_card_for_settings(settings)?
            .to_serializable())
    }

    pub fn resolve_effective_model_parameter_card_for_integrand(
        &self,
        process_id: usize,
        integrand_name: &str,
    ) -> Result<SerializableInputParamCard<F<f64>>> {
        let resolved = self
            .process_list
            .get_integrand(process_id, integrand_name)?;
        match resolved.get_settings() {
            Some(settings) => self.resolve_serializable_model_parameter_card_for_settings(settings),
            None => Ok(self.model_parameters.to_serializable()),
        }
    }

    pub fn resolve_model_for_settings(&self, settings: &RuntimeSettings) -> Result<Model> {
        let mut model = self.model.clone();
        self.resolve_effective_model_parameter_card_for_settings(settings)?
            .apply_to_model(&mut model)?;
        Ok(model)
    }

    pub fn resolve_model_for_integrand(
        &self,
        process_id: usize,
        integrand_name: &str,
    ) -> Result<Model> {
        let resolved = self
            .process_list
            .get_integrand(process_id, integrand_name)?;
        match resolved.get_settings() {
            Some(settings) => self.resolve_model_for_settings(settings),
            None => Ok(self.model.clone()),
        }
    }

    pub fn generate_integrands(
        &mut self,
        global_settings: &GlobalSettings,
        runtime_default: LockedRuntimeSettings,
    ) -> Result<GenerationReports> {
        self.run_generation_with_monitor(global_settings, move |state, generation_pool| {
            let mut reports = state.process_list.preprocess(
                &state.model,
                global_settings,
                &runtime_default,
                generation_pool,
            )?;
            merge_generated_graph_reports(
                &mut reports,
                state.process_list.generate_integrands(
                    &state.model,
                    global_settings,
                    runtime_default,
                    generation_pool,
                )?,
            );
            Ok(reports)
        })
    }

    fn attach_process_id_to_named_reports(
        process_id: usize,
        reports: Vec<NamedGraphGenerationReport>,
    ) -> Vec<GeneratedGraphReport> {
        reports
            .into_iter()
            .map(|report| GeneratedGraphReport {
                process_id,
                integrand_name: report.integrand_name,
                graph_name: report.graph_name,
                stats: report.stats,
            })
            .collect()
    }

    pub fn generate_integrand(
        &mut self,
        global_settings: &GlobalSettings,
        runtime_default: LockedRuntimeSettings,
        process_id: usize,
        integrand_name: Option<String>,
    ) -> Result<GenerationReports> {
        self.run_generation_with_monitor(global_settings, move |state, generation_pool| {
            let p = &mut state.process_list.processes[process_id];
            let process_name = p.definition.folder_name.clone();
            if let Some(name) = &integrand_name {
                let mut reports = Vec::new();
                match &mut p.collection {
                    ProcessCollection::Amplitudes(a) => {
                        if let Some(a) = a.get_mut(name) {
                            begin_phase(
                                GenerationProgressPhase::GraphPreprocessing,
                                GenerationProcessKind::Amplitude,
                                &process_name,
                                &a.name,
                                a.graphs.len(),
                                None,
                            );
                            merge_generated_graph_reports(
                                &mut reports,
                                Self::attach_process_id_to_named_reports(
                                    process_id,
                                    a.preprocess(
                                        &state.model,
                                        &global_settings.generation,
                                        &runtime_default,
                                        generation_pool,
                                    )?,
                                ),
                            );
                            merge_generated_graph_reports(
                                &mut reports,
                                Self::attach_process_id_to_named_reports(
                                    process_id,
                                    a.build_integrand(
                                        &state.model,
                                        &process_name,
                                        global_settings,
                                        runtime_default,
                                        generation_pool,
                                    )?,
                                ),
                            );
                        } else {
                            return Err(eyre!(
                                "No amplitude named '{}' in process id {}",
                                name,
                                process_id
                            ));
                        }
                    }
                    ProcessCollection::CrossSections(cs) => {
                        if let Some(cs) = cs.get_mut(name) {
                            merge_generated_graph_reports(
                                &mut reports,
                                Self::attach_process_id_to_named_reports(
                                    process_id,
                                    cs.preprocess(
                                        &state.model,
                                        &p.definition,
                                        &global_settings.generation,
                                        runtime_default,
                                        generation_pool,
                                    )?,
                                ),
                            );
                            merge_generated_graph_reports(
                                &mut reports,
                                Self::attach_process_id_to_named_reports(
                                    process_id,
                                    cs.build_integrand(
                                        &state.model,
                                        &process_name,
                                        global_settings,
                                        runtime_default,
                                        generation_pool,
                                    )?,
                                ),
                            );
                        } else {
                            return Err(eyre!(
                                "No cross section named '{}' in process id {}",
                                name,
                                process_id
                            ));
                        }
                    }
                }
                Ok(reports)
            } else {
                let mut reports = p.preprocess(
                    &state.model,
                    global_settings,
                    &runtime_default,
                    generation_pool,
                )?;
                merge_generated_graph_reports(
                    &mut reports,
                    p.generate_integrands(
                        &state.model,
                        global_settings,
                        runtime_default,
                        generation_pool,
                    )?,
                );
                Ok(reports)
            }
        })
    }

    fn run_generation_with_monitor<F>(
        &mut self,
        global_settings: &GlobalSettings,
        generation: F,
    ) -> Result<GenerationReports>
    where
        F: FnOnce(&mut Self, &rayon::ThreadPool) -> Result<Vec<GeneratedGraphReport>>,
    {
        let generation_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(global_settings.n_cores.generate)
            .stack_size(GENERATION_THREAD_STACK_SIZE_BYTES)
            .build()?;
        let generation_cores = generation_pool.current_num_threads();
        clear_interrupt_request();
        let mut monitor = GenerationMonitor::start()?;
        let stderr_is_terminal = io::stderr().is_terminal();
        let progress_mode = if generation_cores == 1 && stderr_is_terminal {
            GenerationProgressMode::Detailed
        } else {
            GenerationProgressMode::Aggregate
        };
        let _progress_mode_guard = GenerationProgressModeGuard::set(progress_mode);
        let _progress_observer_guard = if generation_cores > 1 && stderr_is_terminal {
            let reporter = AggregateGenerationProgressReporter::new(
                monitor.current_ram_bytes(),
                monitor.peak_ram_bytes(),
                generation_cores as u64,
            );
            Some(GenerationProgressObserverGuard::set(reporter))
        } else {
            None
        };
        let generation_result = generation(self, &generation_pool);
        let peak_ram_bytes = monitor.finish();
        clear_interrupt_request();

        generation_result.map(|reports| GenerationReports {
            reports,
            resources: GenerationResourceSummary {
                peak_ram_bytes,
                generation_cores,
            },
        })
    }

    pub fn compile_integrands(
        &mut self,
        folder: impl AsRef<Path>,
        override_existing: bool,
        global_settings: &GlobalSettings,
        process_id: Option<usize>,
        integrand_name: Option<String>,
    ) -> Result<Vec<GeneratedGraphReport>> {
        let compile_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(global_settings.n_cores.compile)
            .build()?;
        self.process_list.compile(
            folder,
            override_existing,
            process_id,
            integrand_name,
            &compile_pool,
        )
    }

    pub fn generation_summary(
        &self,
        process_id: usize,
        integrand_name: &str,
    ) -> Option<&IntegrandGenerationSummary> {
        self.generation_summaries
            .get(&IntegrandGenerationSummaryKey {
                process_id,
                integrand_name: integrand_name.to_string(),
            })
    }

    pub fn record_generation_summary(
        &mut self,
        reports: &[GeneratedGraphReport],
        resources: GenerationResourceSummary,
    ) {
        let mut reports_by_integrand: BTreeMap<
            IntegrandGenerationSummaryKey,
            Vec<GeneratedGraphReport>,
        > = BTreeMap::new();

        for report in reports {
            let Some(process) = self.process_list.processes.get(report.process_id) else {
                continue;
            };
            if !process
                .collection
                .get_integrand_names()
                .contains(&report.integrand_name.as_str())
            {
                continue;
            }

            reports_by_integrand
                .entry(IntegrandGenerationSummaryKey {
                    process_id: report.process_id,
                    integrand_name: report.integrand_name.clone(),
                })
                .or_default()
                .push(report.clone());
        }

        for (key, reports) in reports_by_integrand {
            self.generation_summaries.insert(
                key,
                IntegrandGenerationSummary {
                    peak_ram_bytes: resources.peak_ram_bytes,
                    reports,
                },
            );
        }
    }

    fn remove_generation_summaries_for_process(&mut self, process_id: usize) {
        let mut updated = BTreeMap::new();
        for (mut key, mut summary) in std::mem::take(&mut self.generation_summaries) {
            if key.process_id == process_id {
                continue;
            }
            if key.process_id > process_id {
                key.process_id -= 1;
                for report in &mut summary.reports {
                    if report.process_id > process_id {
                        report.process_id -= 1;
                    }
                }
            }
            updated.insert(key, summary);
        }
        self.generation_summaries = updated;
    }

    fn remove_generation_summaries_for_process_without_shifting(&mut self, process_id: usize) {
        self.generation_summaries
            .retain(|key, _| key.process_id != process_id);
    }

    pub fn export_dots(
        &mut self,
        path: impl AsRef<Path>,
        settings: &DotExportSettings,
    ) -> Result<()> {
        self.process_list.export_dot(path, settings)?;
        Ok(())
    }

    pub fn import_graphs(
        &mut self,
        graphs: Vec<Graph>,
        process_name: Option<String>,
        process_id: Option<usize>,
        integrand_name: Option<String>,
        overwrite: bool,
        append: bool,
    ) -> Result<()> {
        let generation_type = if graphs.iter().all(|g| g.initial_state_cut.nedges(g) == 0) {
            GenerationType::Amplitude
        } else if graphs.iter().all(|g| g.initial_state_cut.nedges(g) > 0) {
            GenerationType::CrossSection
        } else {
            return Err(eyre!(
                "Mix of amplitude and cross section graphs in the same file is not supported"
            ));
        };

        let integrand_base_name = integrand_name.clone().unwrap_or("default".to_string());
        let process = if let Some(proc_id) = process_id {
            if proc_id >= self.process_list.processes.len() {
                return Err(eyre!(
                    "Process ID {} invalid, only {} processes available",
                    proc_id,
                    self.process_list.processes.len()
                ));
            }
            Some(&mut self.process_list.processes[proc_id])
        } else {
            let p_name = match process_name {
                Some(n) => n,
                None => {
                    return Err(eyre!(
                        "Either process ID or process name must be provided when importing graphs"
                    ));
                }
            };
            if let Some(existing_proc) = self
                .process_list
                .processes
                .iter_mut()
                .find(|p| p.definition.folder_name == p_name)
            {
                Some(existing_proc)
            } else {
                let process_defintion =
                    ProcessDefinition::from_graph_list(&graphs, generation_type, &self.model)?;
                let process = Process::from_graph_list(
                    p_name,
                    integrand_base_name.clone(),
                    // TODO: avoid clone here
                    graphs.clone(),
                    generation_type,
                    Some(process_defintion),
                    None,
                    &self.model,
                )?;

                self.process_list.add_process(process);
                None
            }
        };
        if let Some(p) = process {
            let existing_names = p.get_integrand_names();
            let integrand_name = if existing_names.contains(&integrand_base_name.as_str()) {
                if append {
                    let mut integrand_i = 0;
                    while existing_names
                        .iter()
                        .any(|ce| *ce == format!("{}_{}", integrand_base_name, integrand_i))
                    {
                        integrand_i += 1;
                    }
                    format!("{}_{}", integrand_base_name, integrand_i)
                } else if overwrite {
                    p.collection.remove_integrand(&integrand_base_name)?;
                    integrand_base_name.clone()
                } else {
                    return Err(eyre!(
                        "Integrand name '{}' already exists in process '{}', use either --overwrite or --append flag when loading graphs",
                        integrand_base_name,
                        p.definition.folder_name
                    ));
                }
            } else {
                integrand_base_name.clone()
            };

            match generation_type {
                GenerationType::Amplitude => p
                    .collection
                    .add_amplitude(Amplitude::from_graph_list(integrand_name.clone(), graphs)?),
                GenerationType::CrossSection => {
                    p.collection
                        .add_cross_section(CrossSection::from_graph_list(
                            integrand_name.clone(),
                            graphs,
                            &self.model,
                        )?)
                }
            }
        }

        Ok(())
    }

    pub fn bench(
        &mut self,
        samples: usize,
        process_id: usize,
        integrand_name: String,
        _n_cores: usize,
    ) -> Result<()> {
        let integrand = self
            .process_list
            .get_integrand_mut(process_id, integrand_name)?;
        let name = integrand.name();

        info!(
            "\nBenchmarking runtime of integrand '{}' over {} samples...\n",
            name.green(),
            samples.to_string().blue()
        );

        let now = Instant::now();
        for _ in 0..samples {
            let _ = integrand.evaluate_sample(
                &Sample::Continuous(
                    F(1.),
                    (0..integrand.get_n_dim())
                        .map(|_| F(rand::random::<f64>()))
                        .collect(),
                ),
                &self.model,
                F(1.),
                1,
                false,
                Complex::new_zero(),
            );
        }
        let total_time = now.elapsed().as_secs_f64();
        info!(
            "\n> Total time: {} s for {} samples, {} ms per sample\n",
            format!("{:.1}", total_time).blue(),
            format!("{}", samples).blue(),
            format!("{:.5}", total_time * 1000. / (samples as f64)).green(),
        );

        Ok(())
    }

    pub fn new(log_dir: impl AsRef<Path>, log_file_name: Option<String>) -> Self {
        super::tracing::init_tracing(log_dir.as_ref().join("logs"), log_file_name);
        let _ = initialise();

        Self {
            model: Model::default(),
            process_list: ProcessList::default(),
            model_parameters: InputParamCard::default(),
            generation_summaries: BTreeMap::new(),
        }
    }

    pub fn new_test() -> Self {
        init_test_tracing();

        Self {
            model: Model::default(),
            process_list: ProcessList::default(),
            model_parameters: InputParamCard::default(),
            generation_summaries: BTreeMap::new(),
        }
    }

    pub fn new_bench() -> Self {
        init_bench_tracing();

        Self {
            model: Model::default(),
            process_list: ProcessList::default(),
            model_parameters: InputParamCard::default(),
            generation_summaries: BTreeMap::new(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RemovedProcess {
    pub process_id: usize,
    pub process_name: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RemovedIntegrand {
    pub process_id: usize,
    pub process_name: String,
    pub integrand_name: String,
    pub removed_empty_process: bool,
}

impl State {
    pub fn load(
        save_path: PathBuf,
        model_path: Option<PathBuf>,
        trace_logs_filename: Option<String>,
    ) -> Result<Self> {
        // let root_folder = root_folder.join("gammaloop_state");
        let manifest = load_state_manifest(&save_path)?;
        run_state_migration_checks(&manifest, &save_path)?;
        // Install GammaLoop's subscriber before importing Symbolica state. Symbolica warnings
        // initialize its fallback subscriber on first use, which would otherwise claim the global
        // tracing dispatch and escape ANSI styling in all subsequent GammaLoop output.
        let mut loaded_state = State::new(&save_path, trace_logs_filename);
        debug!("Loading state manifest version {}", manifest.version);

        let mut model = if let Some(model_path) = &model_path {
            info!("Loading model from {}", model_path.display());
            Model::from_file(model_path)?
        } else {
            let model_dir = save_path.join("model.json");
            info!(
                "Loading model from default location: {}",
                model_dir.display()
            );
            Model::from_file(model_dir)?
        };

        debug!("Loaded model: {}", model.name);

        let input_param_card = if save_path.join("model_parameters.json").exists() {
            let a = InputParamCard::from_file(save_path.join("model_parameters.json"))?;

            let _ = model.apply_param_card(&a);
            a
        } else {
            InputParamCard::default_from_model(&model)
        };

        let symbolica_state = symbolica::state::State::import(
            &mut fs::File::open(save_path.join("symbolica_state.bin"))
                .context("Trying to open symbolica state binary")?,
            None,
        )?;

        let context: GammaLoopContextContainer<'_> = GammaLoopContextContainer {
            state_map: &symbolica_state,
            model: &model,
        };

        let process_list =
            ProcessList::load(&save_path, context).context("Trying to load processList")?;

        loaded_state.process_list = process_list;
        loaded_state.model = model;
        loaded_state.model_parameters = input_param_card;
        loaded_state.generation_summaries =
            load_integrand_generation_summaries(&save_path, &loaded_state.process_list)?;
        Ok(loaded_state)
    }

    pub fn activate_loaded_integrand_backends(
        &mut self,
        allow_symjit_fallback: bool,
    ) -> Result<()> {
        self.process_list
            .activate_loaded_integrand_backends(allow_symjit_fallback)
    }

    pub fn compile(
        &mut self,
        root_folder: &Path,
        override_compiled: bool,
        settings: &GlobalSettings,
    ) -> Result<()> {
        fs::create_dir_all(root_folder)?;

        let compile_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(settings.n_cores.compile)
            .build()?;
        self.process_list
            .compile(root_folder, override_compiled, None, None, &compile_pool)?;
        Ok(())
    }

    fn save_generation_summaries(&self, root_folder: &Path) -> Result<()> {
        for (process_id, process) in self.process_list.processes.iter().enumerate() {
            for integrand_name in process.collection.get_integrand_names() {
                let summary_path =
                    integrand_generation_summary_path(root_folder, process, integrand_name);
                let key = IntegrandGenerationSummaryKey {
                    process_id,
                    integrand_name: integrand_name.to_string(),
                };
                if let Some(summary) = self.generation_summaries.get(&key) {
                    if let Some(parent) = summary_path.parent() {
                        fs::create_dir_all(parent)?;
                    }
                    let raw_summary = serde_json::to_string_pretty(summary).with_context(|| {
                        format!(
                            "Trying to serialize integrand generation summary {}",
                            summary_path.display()
                        )
                    })?;
                    fs::write(&summary_path, raw_summary).with_context(|| {
                        format!(
                            "Trying to write integrand generation summary {}",
                            summary_path.display()
                        )
                    })?;
                } else if summary_path.exists() {
                    fs::remove_file(&summary_path).with_context(|| {
                        format!(
                            "Trying to remove stale integrand generation summary {}",
                            summary_path.display()
                        )
                    })?;
                }
            }
        }

        Ok(())
    }

    pub fn save(
        &mut self,
        root_folder: &Path,
        override_state_file: bool,
        strict: bool,
    ) -> Result<()> {
        // let root_folder = root_folder.join("gammaloop_state");

        // check if the export root exists, if not create it, if it does return error
        let mut selected_root_folder = PathBuf::from(root_folder);
        let mut user_input = String::new();
        if !root_folder.exists() {
            fs::create_dir_all(root_folder)?;
        } else {
            if strict {
                return Err(eyre!(
                    "Export root already exists, please choose a different path or remove the existing directory",
                ));
            }

            if !override_state_file {
                while selected_root_folder.exists() {
                    println!(
                        "Gammaloop export root {} already exists. Specify 'o' for overwriting, 'n' for not saving, or '<NEW_PATH>' to specify where to save current state to:",
                        selected_root_folder.display()
                    );
                    user_input.clear();
                    io::stdin()
                        .read_line(&mut user_input)
                        .expect("Could not read user-specified gammaloop state export destination");
                    //user_input = user_input.trim().into();
                    match user_input.trim() {
                        "o" => break,
                        "n" => {
                            return Ok(());
                        }
                        new_path => {
                            selected_root_folder = PathBuf::from(new_path);
                            continue;
                        }
                    }
                }
            }
        }

        fs::create_dir_all(&selected_root_folder)?;

        let mut state_file =
        // info!("Hi");
            fs::File::create(selected_root_folder.join("symbolica_state.bin"))?;

        symbolica::state::State::export(&mut state_file)?;
        self.process_list
            .save(&selected_root_folder, override_state_file)?;
        self.save_generation_summaries(&selected_root_folder)?;

        // let binary = bincode::encode_to_vec(&self.integrands, bincode::config::standard())?;
        // fs::write(root_folder.join("process_list.bin"), binary)?;?
        self.model
            .to_serializable()
            .to_file(selected_root_folder.join("model.json"), override_state_file)?;
        self.model_parameters.to_file(
            selected_root_folder.join("model_parameters.json"),
            override_state_file,
        )?;
        save_state_manifest(&selected_root_folder)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::fs;

    use gammalooprs::{
        graph::Graph,
        initialisation::test_initialise,
        integrands::process::ActiveF64Backend,
        model::InputParamCard,
        momentum::{Dep, ExternalMomenta, Helicity},
        processes::{
            process::ProcessCollection, RaisedPropagatorScope, RaisedPropagatorSignature,
            SelectionPolarity,
        },
        settings::global::{CompilationMode, FrozenCompilationMode},
        settings::{
            runtime::kinematic::{improvement::PhaseSpaceImprovementSettings, Externals},
            KinematicsSettings, RuntimeSettings,
        },
        utils::{load_generic_model, serde_utils::SHOWDEFAULTS},
    };
    use tempfile::tempdir;

    use crate::commands::{
        display::Display,
        save::SaveState,
        set::{ProcessSetArgs, Set, SetArgs},
    };

    use super::*;

    fn build_generated_scalar_bubble_state_with_external_backend() -> State {
        test_initialise().expect("test initialisation should succeed");
        let mut state = State::new_test();
        state.model = load_generic_model("scalars");
        state.model_parameters = InputParamCard::default_from_model(&state.model);

        let graph_path =
            crate::test_workspace_root().join("tests/resources/graphs/scalar_bubble.dot");
        let graphs = Graph::from_path(&graph_path, &state.model)
            .expect("scalar bubble graph fixture should load");

        state
            .import_graphs(
                graphs,
                Some("scalar_bubble".to_string()),
                None,
                Some("default".to_string()),
                false,
                false,
            )
            .expect("graph import should succeed");

        let mut cli_settings = CLISettings::default();
        cli_settings.global.generation.evaluator.compile = true;
        cli_settings.global.generation.compile.compilation_mode = CompilationMode::Assembly;

        let runtime_defaults = RuntimeSettings::default();
        state
            .generate_integrands(&cli_settings.global, (&runtime_defaults).into())
            .expect("integrand generation should succeed");

        state
    }

    #[test]
    fn aggregate_generation_progress_tracks_counts_and_timings() {
        let current_ram_bytes = Arc::new(AtomicU64::new(256 * 1024 * 1024));
        let peak_ram_bytes = Arc::new(AtomicU64::new(512 * 1024 * 1024));
        let reporter =
            AggregateGenerationProgressReporter::new_hidden(current_ram_bytes, peak_ram_bytes);

        reporter.begin_phase(
            GenerationProgressPhase::GraphPreprocessing,
            GenerationProcessKind::CrossSection,
            "proc",
            "itg",
            2,
            None,
        );
        reporter.graph_started(GenerationProcessKind::CrossSection, "itg", "GL01", None);
        reporter.cuts_discovered("itg", "GL01", 4, 2);
        {
            let state = reporter
                .state
                .lock()
                .expect("aggregate generation progress state mutex is poisoned");
            assert_eq!(state.done_graphs, 0);
            assert_eq!(state.discovered_st_cuts, 4);
            assert_eq!(state.discovered_valid_cuts, 2);
            assert_eq!(
                AggregateGenerationProgressReporter::progress_units(&state),
                (0, 2)
            );
            assert_eq!(
                AggregateGenerationProgressReporter::graph_progress_counts(&state),
                (0, 1, 2)
            );
        }
        reporter.graph_finished(
            GenerationProcessKind::CrossSection,
            "itg",
            "GL01",
            &GraphGenerationStats {
                total_time: Duration::from_secs(3),
                ..GraphGenerationStats::default()
            },
            None,
        );
        {
            let state = reporter
                .state
                .lock()
                .expect("aggregate generation progress state mutex is poisoned");
            assert_eq!(state.done_graphs, 1);
            assert_eq!(state.stats.total_time, Duration::from_secs(3));
            assert_eq!(
                AggregateGenerationProgressReporter::graph_progress_counts(&state),
                (1, 0, 2)
            );
        }

        reporter.begin_phase(
            GenerationProgressPhase::GraphGeneration,
            GenerationProcessKind::CrossSection,
            "proc",
            "itg",
            2,
            Some(3),
        );
        reporter.graph_started(GenerationProcessKind::CrossSection, "itg", "GL01", Some(1));
        reporter.graph_started(GenerationProcessKind::CrossSection, "itg", "GL02", Some(2));
        {
            let state = reporter
                .state
                .lock()
                .expect("aggregate generation progress state mutex is poisoned");
            assert_eq!(
                AggregateGenerationProgressReporter::graph_progress_counts(&state),
                (0, 2, 2)
            );
        }
        reporter.graph_finished(
            GenerationProcessKind::CrossSection,
            "itg",
            "GL01",
            &GraphGenerationStats {
                evaluator_count: 2,
                total_time: Duration::from_secs(4),
                evaluator_spenso_time: Duration::from_secs(1),
                evaluator_symbolica_time: Duration::from_secs(1),
                evaluator_compile_time: Duration::ZERO,
            },
            None,
        );
        reporter.cut_finished("itg", "GL01", 1);
        reporter.graph_finished(
            GenerationProcessKind::CrossSection,
            "itg",
            "GL02",
            &GraphGenerationStats {
                evaluator_count: 3,
                total_time: Duration::from_secs(6),
                evaluator_spenso_time: Duration::from_secs(2),
                evaluator_symbolica_time: Duration::ZERO,
                evaluator_compile_time: Duration::ZERO,
            },
            None,
        );
        reporter.cut_finished("itg", "GL02", 2);

        {
            let state = reporter
                .state
                .lock()
                .expect("aggregate generation progress state mutex is poisoned");
            assert_eq!(state.done_graphs, 2);
            assert_eq!(state.done_cuts, 3);
            assert!(state.active_graphs.is_empty());
            assert_eq!(state.last_graph.as_deref(), Some("GL02"));
            assert_eq!(state.stats.evaluator_count, 5);
            assert_eq!(state.stats.total_time, Duration::from_secs(10));
            assert_eq!(state.stats.expression_build_time(), Duration::from_secs(6));
            assert_eq!(state.stats.evaluator_spenso_time, Duration::from_secs(3));
            assert_eq!(state.stats.evaluator_symbolica_time, Duration::from_secs(1));
            assert_eq!(
                AggregateGenerationProgressReporter::progress_units(&state),
                (5, 5)
            );
            assert_eq!(
                AggregateGenerationProgressReporter::graph_progress_counts(&state),
                (2, 0, 2)
            );
        }

        reporter.backend_started(GenerationProcessKind::CrossSection, "itg", 2);
        reporter.backend_finished(
            GenerationProcessKind::CrossSection,
            "itg",
            Duration::from_secs(2),
        );

        let state = reporter
            .state
            .lock()
            .expect("aggregate generation progress state mutex is poisoned");
        assert_eq!(state.phase, Some(GenerationProgressPhase::Backend));
        assert_eq!(state.done_graphs, 2);
        assert_eq!(state.stats.total_time, Duration::from_secs(12));
        assert_eq!(state.stats.evaluator_compile_time, Duration::from_secs(2));
    }

    #[test]
    fn aggregate_generation_progress_tracks_amplitude_preprocessing() {
        let reporter = AggregateGenerationProgressReporter::new_hidden(
            Arc::new(AtomicU64::new(0)),
            Arc::new(AtomicU64::new(0)),
        );

        reporter.begin_phase(
            GenerationProgressPhase::GraphPreprocessing,
            GenerationProcessKind::Amplitude,
            "proc",
            "amp",
            3,
            None,
        );
        reporter.graph_started(GenerationProcessKind::Amplitude, "amp", "GL01", None);
        reporter.graph_started(GenerationProcessKind::Amplitude, "amp", "GL02", None);

        {
            let state = reporter
                .state
                .lock()
                .expect("aggregate generation progress state mutex is poisoned");
            assert_eq!(state.kind, Some(GenerationProcessKind::Amplitude));
            assert_eq!(state.process, "proc");
            assert_eq!(state.integrand, "amp");
            assert_eq!(
                AggregateGenerationProgressReporter::progress_units(&state),
                (0, 3)
            );
            assert_eq!(
                AggregateGenerationProgressReporter::graph_progress_counts(&state),
                (0, 2, 3)
            );
        }

        reporter.graph_finished(
            GenerationProcessKind::Amplitude,
            "amp",
            "GL01",
            &GraphGenerationStats::default(),
            None,
        );

        let state = reporter
            .state
            .lock()
            .expect("aggregate generation progress state mutex is poisoned");
        assert_eq!(
            AggregateGenerationProgressReporter::progress_units(&state),
            (1, 3)
        );
        assert_eq!(
            AggregateGenerationProgressReporter::graph_progress_counts(&state),
            (1, 1, 3)
        );
        assert_eq!(state.last_graph.as_deref(), Some("GL01"));
    }

    #[test]
    fn aggregate_generation_progress_formats_memory_and_time_shares() {
        assert_eq!(
            AggregateGenerationProgressReporter::progress_memory(512 * 1024 * 1024),
            "   512 MiB"
        );
        assert_eq!(
            AggregateGenerationProgressReporter::progress_memory(5 * 1024 * 1024 * 1024),
            "  5.00 GiB"
        );

        assert_eq!(
            AggregateGenerationProgressReporter::progress_percent(0.0),
            "0%"
        );
        assert_eq!(
            AggregateGenerationProgressReporter::progress_percent(0.009),
            "0%"
        );
        assert_eq!(
            AggregateGenerationProgressReporter::progress_percent(0.024),
            "0.024%"
        );
        assert_eq!(
            AggregateGenerationProgressReporter::progress_percent(0.24),
            "0.24%"
        );
        assert_eq!(
            AggregateGenerationProgressReporter::progress_percent(2.4),
            "2.4%"
        );
        assert_eq!(
            AggregateGenerationProgressReporter::progress_percent(54.0),
            "54%"
        );

        assert_eq!(
            AggregateGenerationProgressReporter::progress_time_share(
                Duration::ZERO,
                Duration::ZERO,
            ),
            "--"
        );
        assert_eq!(
            AggregateGenerationProgressReporter::progress_time_share(
                Duration::ZERO,
                Duration::from_secs(4),
            ),
            "0%"
        );
        assert_eq!(
            AggregateGenerationProgressReporter::progress_time_share(
                Duration::from_secs(1),
                Duration::from_secs(4),
            ),
            "25%"
        );
    }

    fn build_scalar_bubble_diagram_state() -> State {
        test_initialise().expect("test initialisation should succeed");
        let mut state = State::new_test();
        state.model = load_generic_model("scalars");
        state.model_parameters = InputParamCard::default_from_model(&state.model);

        let graph_path =
            crate::test_workspace_root().join("tests/resources/graphs/scalar_bubble.dot");
        let graphs = Graph::from_path(&graph_path, &state.model)
            .expect("scalar bubble graph fixture should load");
        state
            .import_graphs(
                graphs,
                Some("scalar_bubble".to_string()),
                None,
                Some("default".to_string()),
                false,
                false,
            )
            .expect("graph import should succeed");
        state
    }

    fn scalar_bubble_master_graph_name(state: &State, integrand_name: &str) -> String {
        let process = &state.process_list.processes[0];
        match &process.collection {
            ProcessCollection::Amplitudes(amplitudes) => amplitudes
                .get(integrand_name)
                .expect("amplitude should exist")
                .graphs[0]
                .graph
                .name
                .clone(),
            ProcessCollection::CrossSections(cross_sections) => cross_sections
                .get(integrand_name)
                .expect("cross section should exist")
                .supergraphs[0]
                .graph
                .name
                .clone(),
        }
    }

    #[test]
    fn select_copy_mode_creates_pregenerated_target_without_mutating_source() {
        let _guard = crate::LOG_TEST_MUTEX
            .lock()
            .unwrap_or_else(|err| err.into_inner());
        let mut state = build_scalar_bubble_diagram_state();
        let temp = tempdir().unwrap();
        let graph_name = scalar_bubble_master_graph_name(&state, "default");
        let selection = GraphGroupSelectionSpec::new().with_master_graph_names(vec![graph_name]);
        let target = GraphGroupSelectionTarget::copy(None, Some("selected".to_string()), false);

        let selected = state
            .select_integrand_graph_groups(
                Some(&ProcessRef::Name("scalar_bubble".to_string())),
                Some(&"default".to_string()),
                &selection,
                &target,
                GraphGroupSelectionContext::new(&GenerationSettings::default(), temp.path(), false),
            )
            .unwrap();

        assert!(selected.copied_to_output);
        assert_eq!(selected.process_id, 0);
        assert_eq!(selected.integrand_name, "selected");
        let process = &state.process_list.processes[0];
        assert!(process
            .collection
            .get_integrand_names()
            .contains(&"default"));
        assert!(process
            .collection
            .get_integrand_names()
            .contains(&"selected"));
        match &process.collection {
            ProcessCollection::Amplitudes(amplitudes) => {
                assert!(amplitudes["default"].integrand.is_none());
                assert!(amplitudes["selected"].integrand.is_none());
                assert_eq!(amplitudes["default"].graphs.len(), 1);
                assert_eq!(amplitudes["selected"].graphs.len(), 1);
            }
            ProcessCollection::CrossSections(cross_sections) => {
                assert!(cross_sections["default"].integrand.is_none());
                assert!(cross_sections["selected"].integrand.is_none());
                assert_eq!(cross_sections["default"].supergraphs.len(), 1);
                assert_eq!(cross_sections["selected"].supergraphs.len(), 1);
            }
        }
    }

    #[test]
    fn select_copy_mode_rejects_and_clears_existing_target_integrand() {
        let _guard = crate::LOG_TEST_MUTEX
            .lock()
            .unwrap_or_else(|err| err.into_inner());
        let mut state = build_scalar_bubble_diagram_state();
        let temp = tempdir().unwrap();
        let graph_name = scalar_bubble_master_graph_name(&state, "default");
        let selection = GraphGroupSelectionSpec::new().with_master_graph_names(vec![graph_name]);
        let target = GraphGroupSelectionTarget::copy(None, Some("selected".to_string()), false);

        state
            .select_integrand_graph_groups(
                Some(&ProcessRef::Name("scalar_bubble".to_string())),
                Some(&"default".to_string()),
                &selection,
                &target,
                GraphGroupSelectionContext::new(&GenerationSettings::default(), temp.path(), false),
            )
            .unwrap();
        let err = state
            .select_integrand_graph_groups(
                Some(&ProcessRef::Name("scalar_bubble".to_string())),
                Some(&"default".to_string()),
                &selection,
                &target,
                GraphGroupSelectionContext::new(&GenerationSettings::default(), temp.path(), false),
            )
            .unwrap_err();
        assert!(format!("{err}").contains("already exists"));

        let process = &state.process_list.processes[0];
        let kind_folder = match &process.collection {
            ProcessCollection::Amplitudes(_) => "amplitudes",
            ProcessCollection::CrossSections(_) => "cross_sections",
        };
        let stale_integrand_folder = temp
            .path()
            .join("processes")
            .join(kind_folder)
            .join("scalar_bubble")
            .join("selected");
        fs::create_dir_all(stale_integrand_folder.join("integrand")).unwrap();
        state.generation_summaries.insert(
            IntegrandGenerationSummaryKey {
                process_id: 0,
                integrand_name: "selected".to_string(),
            },
            IntegrandGenerationSummary {
                peak_ram_bytes: 1,
                reports: Vec::new(),
            },
        );

        let clear_target =
            GraphGroupSelectionTarget::copy(None, Some("selected".to_string()), true);
        let selected = state
            .select_integrand_graph_groups(
                Some(&ProcessRef::Name("scalar_bubble".to_string())),
                Some(&"default".to_string()),
                &selection,
                &clear_target,
                GraphGroupSelectionContext::new(&GenerationSettings::default(), temp.path(), false),
            )
            .unwrap();
        assert!(selected.replaced_existing_target);
        assert!(selected.removed_target_artifacts);
        assert!(!stale_integrand_folder.exists());
        assert!(!state
            .generation_summaries
            .contains_key(&IntegrandGenerationSummaryKey {
                process_id: 0,
                integrand_name: "selected".to_string(),
            }));

        let err = state
            .select_integrand_graph_groups(
                Some(&ProcessRef::Name("scalar_bubble".to_string())),
                Some(&"default".to_string()),
                &selection,
                &clear_target,
                GraphGroupSelectionContext::new(&GenerationSettings::default(), temp.path(), true),
            )
            .unwrap_err();
        assert!(format!("{err}").contains("--read-only-state"));
    }

    #[test]
    fn select_amplitude_graphs_mode_rejects_amplitude_integrands() {
        let _guard = crate::LOG_TEST_MUTEX
            .lock()
            .unwrap_or_else(|err| err.into_inner());
        let mut state = build_scalar_bubble_diagram_state();
        let temp = tempdir().unwrap();
        let graph_name = scalar_bubble_master_graph_name(&state, "default");
        let selection = GraphGroupSelectionSpec::new()
            .with_mode(GraphGroupSelectionMode::CrossSectionAmplitudeGraphs)
            .with_master_graph_names(vec![graph_name]);
        let target = GraphGroupSelectionTarget::in_place();

        let err = state
            .select_integrand_graph_groups(
                Some(&ProcessRef::Name("scalar_bubble".to_string())),
                Some(&"default".to_string()),
                &selection,
                &target,
                GraphGroupSelectionContext::new(&GenerationSettings::default(), temp.path(), false),
            )
            .unwrap_err();

        assert!(format!("{err}").contains("--amplitude-graphs"));
    }

    #[test]
    fn select_raised_cut_filters_reject_amplitude_integrands() {
        let _guard = crate::LOG_TEST_MUTEX
            .lock()
            .unwrap_or_else(|err| err.into_inner());
        let mut state = build_scalar_bubble_diagram_state();
        let temp = tempdir().unwrap();
        let selection = GraphGroupSelectionSpec::new().with_raised_cut_signatures(
            SelectionPolarity::With,
            RaisedPropagatorScope::All,
            vec![RaisedPropagatorSignature::from_str("[]").unwrap()],
        );
        let target = GraphGroupSelectionTarget::in_place();

        let err = state
            .select_integrand_graph_groups(
                Some(&ProcessRef::Name("scalar_bubble".to_string())),
                Some(&"default".to_string()),
                &selection,
                &target,
                GraphGroupSelectionContext::new(&GenerationSettings::default(), temp.path(), false),
            )
            .unwrap_err();

        assert!(format!("{err}").contains("cross-section integrands"));
    }

    #[test]
    fn test_run_history() {
        use crate::state::RunHistory;
        //SHOWDEFAULTS.store(true, std::sync::atomic::Ordering::Relaxed);
        let mut run_history: RunHistory = Default::default();
        let kinematics_settings = KinematicsSettings {
            e_cm: 100.0,
            externals: Externals::Constant {
                momenta: vec![
                    ExternalMomenta::Independent([F(1.), F(2.), F(3.), F(4.)]),
                    ExternalMomenta::Dependent(Dep::Dep),
                ],
                improvement_settings: PhaseSpaceImprovementSettings::default(),
                helicities: vec![Helicity::PLUS, Helicity::MINUS],
                f_64_cache: None,
                f_128_cache: None,
            },
        };

        run_history.push(Commands::Set(Set::Global {
            input: SetArgs::Stored,
        }));

        run_history.default_runtime_settings.kinematics = kinematics_settings;
        set_serialize_commands_as_strings(true);
        let toml = toml::to_string_pretty(&run_history).unwrap();
        println!("{}", toml);
        let deserialized: RunHistory = toml::from_str(&toml).unwrap();
        assert_eq!(run_history, deserialized);
        SHOWDEFAULTS.store(false, std::sync::atomic::Ordering::Relaxed);

        run_history.to_file("test_path.toml", true).unwrap();
        let deserialized_from_file = RunHistory::from_file("test_path.toml", " ").unwrap();
        assert_eq!(run_history, deserialized_from_file);
    }

    #[test]
    fn run_history_applies_default_runtime_before_commands() {
        let mut run_history = RunHistory {
            default_runtime_settings: toml::from_str(
                r#"
[sampling]
graphs = "monte_carlo"
orientations = "summed"
lmb_multichanneling = false
lmb_channels = "summed"
coordinate_system = "tropical"
mapping = "linear"
b = 1.0
"#,
            )
            .unwrap(),
            ..Default::default()
        };

        let mut state = State::new_test();
        let mut cli_settings = CLISettings::default();
        let mut default_runtime_settings = RuntimeSettings::default();

        let _ = run_history
            .run(&mut state, &mut cli_settings, &mut default_runtime_settings)
            .unwrap();

        assert_eq!(
            default_runtime_settings.sampling,
            run_history.default_runtime_settings.sampling
        );
    }

    #[test]
    fn test_command_history_serialization() {
        use super::{set_serialize_commands_as_strings, CommandHistory};
        use crate::commands::Commands;

        // Test basic construction
        let cmd_history = CommandHistory::new(Commands::Quit(SaveState::default()));
        assert_eq!(cmd_history.raw_string, None);

        // Test with raw string
        let cmd_history_with_raw =
            CommandHistory::new_with_raw(Commands::Quit(SaveState::default()), "quit".to_string());
        assert_eq!(cmd_history_with_raw.raw_string, Some("quit".to_string()));

        // Test serialization as Commands (default behavior)
        set_serialize_commands_as_strings(false);

        let json = serde_json::to_string(&cmd_history).unwrap();
        let deserialized_json: CommandHistory = serde_json::from_str(&json).unwrap();
        let toml = toml::to_string(&cmd_history).unwrap();
        let deserialized_toml: CommandHistory = toml::from_str(&toml).unwrap();
        assert_eq!(cmd_history, deserialized_json);
        assert_eq!(cmd_history, deserialized_toml);

        // Test serialization as string
        set_serialize_commands_as_strings(true);

        let json_string = serde_json::to_string_pretty(&cmd_history_with_raw).unwrap();
        assert!(json_string.contains("quit"));

        // Reset flag
        set_serialize_commands_as_strings(false);
    }

    #[test]
    fn test_command_history_toml_and_json_formats() {
        use super::{set_serialize_commands_as_strings, CommandHistory};
        use crate::commands::Commands;

        // Test different command types
        let quit_cmd = CommandHistory::new(Commands::Quit(SaveState::default()));
        let quit_with_raw =
            CommandHistory::new_with_raw(Commands::Quit(SaveState::default()), "quit".to_string());

        // Test JSON serialization/deserialization
        {
            // Test Commands format in JSON
            set_serialize_commands_as_strings(false);
            let json = serde_json::to_string_pretty(&quit_cmd).unwrap();
            let deserialized: CommandHistory = serde_json::from_str(&json).unwrap();
            assert_eq!(quit_cmd, deserialized);

            // Test string format in JSON
            set_serialize_commands_as_strings(true);
            let json_string = serde_json::to_string_pretty(&quit_with_raw).unwrap();
            assert!(json_string.contains("quit"));
            let deserialized_string: CommandHistory = serde_json::from_str(&json_string).unwrap();
            assert_eq!(quit_with_raw, deserialized_string);
        }

        // Test TOML serialization/deserialization with wrapper struct
        {
            #[derive(serde::Serialize, serde::Deserialize)]
            struct CommandWrapper {
                command: CommandHistory,
            }

            // Test Commands format in TOML
            set_serialize_commands_as_strings(false);
            let wrapper = CommandWrapper {
                command: quit_cmd.clone(),
            };
            let toml = toml::to_string_pretty(&wrapper).unwrap();
            let deserialized_wrapper: CommandWrapper = toml::from_str(&toml).unwrap();
            assert_eq!(quit_cmd, deserialized_wrapper.command);

            // Test string format in TOML
            set_serialize_commands_as_strings(true);
            let wrapper_string = CommandWrapper {
                command: quit_with_raw.clone(),
            };
            let toml_string = toml::to_string_pretty(&wrapper_string).unwrap();
            assert!(toml_string.contains("quit"));
            let deserialized_string_wrapper: CommandWrapper = toml::from_str(&toml_string).unwrap();
            assert_eq!(quit_with_raw, deserialized_string_wrapper.command);
        }

        // Test cross-format compatibility: serialize in one format, deserialize in another
        {
            set_serialize_commands_as_strings(false);

            // Serialize as JSON, deserialize the Commands directly from JSON Value
            let json = serde_json::to_string(&quit_cmd).unwrap();
            let json_value: serde_json::Value = serde_json::from_str(&json).unwrap();
            let from_json: CommandHistory = serde_json::from_value(json_value).unwrap();
            assert_eq!(quit_cmd, from_json);
        }
        // Reset flag
        set_serialize_commands_as_strings(false);
    }

    #[test]
    fn run_history_push_with_raw_skips_quit_and_definition_commands() {
        use super::RunHistory;
        use crate::commands::{run::Run, StartCommandsBlock};

        let mut run_history = RunHistory::default();
        run_history.push_with_raw(
            Commands::Run(Run {
                block_names: vec!["block_a".to_string()],
                commands: None,
            }),
            Some("run block_a".to_string()),
        );
        run_history.push_with_raw(
            Commands::Quit(SaveState::default()),
            Some("quit".to_string()),
        );
        run_history.push_with_raw(
            Commands::StartCommandsBlock(StartCommandsBlock {
                name: "block_a".to_string(),
            }),
            Some("start_commands_block block_a".to_string()),
        );
        run_history.push_with_raw(
            Commands::FinishCommandsBlock,
            Some("finish_commands_block".to_string()),
        );

        assert_eq!(run_history.commands.len(), 1);
        assert_eq!(
            run_history.commands[0].raw_string.as_deref(),
            Some("run block_a")
        );
    }

    #[test]
    fn run_history_filtered_for_save_preserves_command_blocks() {
        let mut run_history = RunHistory {
            command_blocks: vec![
                CommandsBlock {
                    name: "generate".to_string(),
                    commands: vec![CommandHistory::new_with_raw(
                        Commands::Display(Display::Processes),
                        "display processes".to_string(),
                    )],
                },
                CommandsBlock {
                    name: "integrate".to_string(),
                    commands: vec![CommandHistory::new_with_raw(
                        CommandHistory::from_raw_string("quit -o").unwrap().command,
                        "quit -o".to_string(),
                    )],
                },
            ],
            ..Default::default()
        };
        run_history.push_with_raw(
            CommandHistory::from_raw_string("display processes")
                .unwrap()
                .command,
            Some("display processes".to_string()),
        );

        let filtered = run_history.filtered_for_save();
        set_serialize_commands_as_strings(true);
        let toml = toml::to_string_pretty(&filtered).unwrap();
        set_serialize_commands_as_strings(false);

        assert_eq!(filtered.command_blocks.len(), 2);
        assert!(toml.contains("commands = ["));
        assert!(toml.contains("[[command_blocks]]"));
        assert!(toml.contains("quit -o"));
    }

    #[test]
    fn command_history_parses_multiline_set_string() {
        let raw = "set process -p epem_a_tth -i LO string '[integrator]\nn_start = 1000\n'";
        let cmd = CommandHistory::from_raw_string(raw).unwrap();
        assert_eq!(cmd.raw_string.as_deref(), Some(raw));

        match cmd.command {
            Commands::Set(Set::Process { input, .. }) => match input {
                ProcessSetArgs::String { string } => {
                    assert_eq!(string, "[integrator]\nn_start = 1000\n");
                }
                other => panic!("Expected string set input, got {other:?}"),
            },
            other => panic!("Expected set process command, got {other:?}"),
        }
    }

    #[test]
    fn command_history_parses_hash_process_refs() {
        let cmd = CommandHistory::from_raw_string("display integrand -p #12").unwrap();
        match cmd.command {
            Commands::Display(Display::Integrands {
                process,
                integrand_name,
                graphs,
                categories,
                hide_non_existing_thresholds,
            }) => {
                assert_eq!(process, Some(ProcessRef::Id(12)));
                assert_eq!(integrand_name, None);
                assert!(graphs.is_empty());
                assert!(categories.is_empty());
                assert!(!hide_non_existing_thresholds);
            }
            other => panic!("Expected display integrand command, got {other:?}"),
        }
    }

    #[test]
    fn command_history_requires_explicit_generate_mode() {
        assert!(CommandHistory::from_raw_string("generate e+ e- > d d~").is_err());
    }

    #[test]
    fn run_history_parses_triple_quoted_set_kv_command() {
        let toml = r#"
commands = [
    """set default-runtime kv kinematics.externals='{"type":"constant","data":{"momenta":[[1.0,2.0,3.0,4.0],[5.0,6.0,-7.0,-8.0]],"helicities":[1,1]}}'""",
]

[default_runtime_settings.general]
integral_unit = "picobarn"
"#;

        let run_history: RunHistory = toml::from_str(toml).unwrap();
        assert_eq!(run_history.commands.len(), 1);

        let expected_cmd = r#"set default-runtime kv kinematics.externals='{"type":"constant","data":{"momenta":[[1.0,2.0,3.0,4.0],[5.0,6.0,-7.0,-8.0]],"helicities":[1,1]}}'"#;
        let command_history = &run_history.commands[0];
        assert_eq!(command_history.raw_string.as_deref(), Some(expected_cmd));

        match &command_history.command {
            Commands::Set(Set::DefaultRuntime {
                input: SetArgs::Kv { pairs },
            }) => {
                assert_eq!(pairs.len(), 1);
                assert_eq!(pairs[0].key, "kinematics.externals");
                assert_eq!(
                    pairs[0].value,
                    r#"{"type":"constant","data":{"momenta":[[1.0,2.0,3.0,4.0],[5.0,6.0,-7.0,-8.0]],"helicities":[1,1]}}"#
                );
            }
            other => panic!("Expected set default-runtime kv command, got {other:?}"),
        }

        assert_eq!(
            run_history.default_runtime_settings.general.integral_unit,
            gammalooprs::settings::runtime::IntegralUnit::Picobarn
        );
    }

    #[test]
    fn run_history_load_preserves_command_blocks() {
        let temp = tempdir().unwrap();
        let run_path = temp.path().join("run.toml");
        fs::write(
            &run_path,
            r#"
[[command_blocks]]
name = "zeta"
commands = ["quit -n"]

[[command_blocks]]
name = "alpha"
commands = ["quit -o"]
"#,
        )
        .unwrap();

        let run_history = RunHistory::load(&run_path).unwrap();
        assert!(run_history.commands.is_empty());
        assert_eq!(run_history.command_blocks.len(), 2);
    }

    #[test]
    fn run_history_selects_named_command_blocks_in_order() {
        let temp = tempdir().unwrap();
        let run_path = temp.path().join("run.toml");
        fs::write(
            &run_path,
            r#"
[[command_blocks]]
name = "first"
commands = ["quit -n"]

[[command_blocks]]
name = "second"
commands = ["quit -o"]
"#,
        )
        .unwrap();

        let requested = vec!["second".to_string(), "first".to_string()];
        let run_history = RunHistory::load(&run_path).unwrap();
        let selected = run_history
            .select_command_blocks(requested.as_slice())
            .unwrap();
        assert_eq!(selected.len(), 2);
        assert_eq!(
            selected[0].commands[0].raw_string.as_deref(),
            Some("quit -o")
        );
        assert_eq!(
            selected[1].commands[0].raw_string.as_deref(),
            Some("quit -n")
        );
    }

    #[test]
    fn run_history_selection_rejects_unknown_command_block() {
        let temp = tempdir().unwrap();
        let run_path = temp.path().join("run.toml");
        fs::write(
            &run_path,
            r#"
[[command_blocks]]
name = "first"
commands = ["quit -n"]
"#,
        )
        .unwrap();

        let requested = vec!["missing".to_string()];
        let run_history = RunHistory::load(&run_path).unwrap();
        let err = run_history
            .select_command_blocks(requested.as_slice())
            .unwrap_err();
        let message = format!("{err}");
        assert!(message.contains("Unknown command block 'missing'"));
        assert!(message.contains("first"));
    }

    #[test]
    fn run_history_selection_rejects_missing_command_block_when_only_legacy_commands_exist() {
        let temp = tempdir().unwrap();
        let run_path = temp.path().join("run.toml");
        fs::write(
            &run_path,
            r#"
commands = ["quit -o"]
"#,
        )
        .unwrap();

        let requested = vec!["first".to_string()];
        let run_history = RunHistory::load(&run_path).unwrap();
        let err = run_history
            .select_command_blocks(requested.as_slice())
            .unwrap_err();
        assert!(format!("{err}").contains("Unknown command block"));
    }

    #[test]
    fn run_history_load_accepts_commands_and_command_blocks_together() {
        let temp = tempdir().unwrap();
        let run_path = temp.path().join("run.toml");
        fs::write(
            &run_path,
            r#"
commands = ["quit -o"]

[[command_blocks]]
name = "first"
commands = ["quit -n"]
"#,
        )
        .unwrap();

        let run_history = RunHistory::load(&run_path).unwrap();
        assert_eq!(run_history.commands.len(), 1);
        assert_eq!(run_history.command_blocks.len(), 1);
        assert_eq!(
            run_history.commands[0].raw_string.as_deref(),
            Some("quit -o")
        );
    }

    #[test]
    fn run_history_load_rejects_duplicate_command_block_names() {
        let temp = tempdir().unwrap();
        let run_path = temp.path().join("run.toml");
        fs::write(
            &run_path,
            r#"
[[command_blocks]]
name = "first"
commands = ["quit -o"]

[[command_blocks]]
name = "first"
commands = ["quit -n"]
"#,
        )
        .unwrap();

        let err = RunHistory::load(&run_path).unwrap_err();
        assert!(format!("{err}").contains("duplicate block name"));
    }

    #[test]
    fn state_manifest_roundtrip_current_version() {
        let temp = tempdir().unwrap();
        save_state_manifest(temp.path()).unwrap();

        let manifest = load_state_manifest(temp.path()).unwrap();
        assert_eq!(manifest.version, CURRENT_STATE_MANIFEST_VERSION);
    }

    #[test]
    fn state_manifest_rejects_future_versions() {
        let temp = tempdir().unwrap();
        let future_manifest = StateManifest {
            version: CURRENT_STATE_MANIFEST_VERSION + 1,
        };
        fs::write(
            temp.path().join(STATE_MANIFEST_FILE),
            toml::to_string_pretty(&future_manifest).unwrap(),
        )
        .unwrap();

        let err = load_state_manifest(temp.path()).unwrap_err();
        assert!(format!("{err}").contains("newer than this binary supports"));
    }

    #[test]
    fn state_folder_classifies_saved_layout() {
        let temp = tempdir().unwrap();
        save_state_manifest(temp.path()).unwrap();
        fs::write(temp.path().join("model.json"), "{}").unwrap();
        fs::write(temp.path().join("symbolica_state.bin"), []).unwrap();
        fs::create_dir_all(temp.path().join("processes")).unwrap();

        assert_eq!(
            classify_state_folder(temp.path()).unwrap(),
            StateFolderKind::Saved
        );
    }

    #[test]
    fn state_folder_classifies_logs_only_folder_as_scratch() {
        let temp = tempdir().unwrap();
        fs::create_dir_all(temp.path().join("logs")).unwrap();
        fs::write(temp.path().join("logs").join("gammalog.jsonl"), "").unwrap();

        assert_eq!(
            classify_state_folder(temp.path()).unwrap(),
            StateFolderKind::Scratch
        );
    }

    #[test]
    fn state_folder_classifies_non_manifest_contents_as_unmanifested() {
        let temp = tempdir().unwrap();
        fs::create_dir_all(temp.path().join("processes").join("amplitudes")).unwrap();

        assert_eq!(
            classify_state_folder(temp.path()).unwrap(),
            StateFolderKind::Unmanifested
        );
    }

    #[test]
    fn activate_loaded_integrand_backends_falls_back_to_eager_when_external_artifacts_are_missing()
    {
        let mut state = build_generated_scalar_bubble_state_with_external_backend();

        state.activate_loaded_integrand_backends(false).unwrap();

        for process in &state.process_list.processes {
            match &process.collection {
                ProcessCollection::Amplitudes(amplitudes) => {
                    for amplitude in amplitudes.values() {
                        let integrand = amplitude.integrand.as_ref().unwrap();
                        if matches!(
                            integrand.frozen_compilation(),
                            FrozenCompilationMode::Cpp(_) | FrozenCompilationMode::Assembly(_)
                        ) {
                            assert_eq!(integrand.active_f64_backend(), ActiveF64Backend::Eager);
                        }
                    }
                }
                ProcessCollection::CrossSections(cross_sections) => {
                    for cross_section in cross_sections.values() {
                        let integrand = cross_section.integrand.as_ref().unwrap();
                        if matches!(
                            integrand.frozen_compilation(),
                            FrozenCompilationMode::Cpp(_) | FrozenCompilationMode::Assembly(_)
                        ) {
                            assert_eq!(integrand.active_f64_backend(), ActiveF64Backend::Eager);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn resolve_effective_model_parameter_card_overlays_runtime_model_settings() {
        let mut state = State::new_test();
        state.model = load_generic_model("scalars");
        state.model_parameters = InputParamCard::default_from_model(&state.model);

        let mut settings = RuntimeSettings::default();
        settings
            .model
            .external_parameters
            .insert("mass_scalar_2".to_string(), (F(7.5), F(0.0)));

        let resolved = state
            .resolve_effective_model_parameter_card_for_settings(&settings)
            .unwrap();

        assert_eq!(
            resolved[&UFOSymbol::from("mass_scalar_2")],
            Complex::new(F(7.5), F(0.0))
        );
        assert_eq!(
            resolved[&UFOSymbol::from("mass_scalar_1")],
            state.model_parameters[&UFOSymbol::from("mass_scalar_1")]
        );
    }

    #[test]
    fn resolve_effective_model_parameter_card_rejects_non_overridable_parameters() {
        let mut state = State::new_test();
        state.model = load_generic_model("scalars");
        state.model_parameters = InputParamCard::default_from_model(&state.model);
        state
            .model_parameters
            .remove(&UFOSymbol::from("mass_scalar_2"));

        let mut settings = RuntimeSettings::default();
        settings
            .model
            .external_parameters
            .insert("mass_scalar_2".to_string(), (F(7.5), F(0.0)));

        let err = state
            .resolve_effective_model_parameter_card_for_settings(&settings)
            .unwrap_err();

        assert!(err
            .to_string()
            .contains("cannot be overridden because it is not present"));
    }
}

#[derive(Args, Debug, Clone)]
pub struct ExistingArgs {
    pub process_id: u32,
    pub name: Option<String>,
}
