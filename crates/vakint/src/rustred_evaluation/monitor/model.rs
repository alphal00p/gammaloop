//! The typed monitor model: one value read by the terminal renderer, the
//! supervisor viewer and the published JSON (`vakint.monitor.v2`). Labels and
//! styling live in the renderer; nothing here is matched by display text.

use rustred::solver::CandidateOrderedPhaseKind;
use serde::{Deserialize, Serialize};

use super::{ApiStatus, MonitorStatus, ScopeKind};

/// Version of [`MonitorModel`] inside the published `monitor` object.
pub(super) const MODEL_SCHEMA: &str = "vakint.monitor.v2";

/// Everything the monitor shows, in display-independent units (seconds,
/// counts, bytes). Optional parts are absent when nothing was observed.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub(crate) struct MonitorModel {
    pub scope: ScopeKind,
    pub status: MonitorStatus,
    pub api_status: ApiStatus,
    /// Wall time since the monitor started.
    pub elapsed_seconds: f64,
    pub run: Option<RunConfig>,
    pub stage: StageModel,
    /// Cumulative wall time per stage kind, running stage included.
    pub stages: Vec<StageTotal>,
    pub batch: Option<Counter>,
    pub cohort: Option<Counter>,
    pub ordered: Option<OrderedModel>,
    pub dfs: Option<DfsModel>,
    pub memory: MemoryModel,
    pub details: Details,
    /// Most recent notable events, oldest first.
    pub events: Vec<String>,
}

/// Settings of the evaluation being observed, as supplied by the caller.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(default)]
pub(crate) struct RunConfig {
    pub n_cores: usize,
    pub arithmetic: String,
    pub step_cache: String,
}

/// Broad steps of an evaluation, in pipeline order.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum StageKind {
    #[default]
    Setup,
    DataLoad,
    Ingress,
    Lowering,
    Native,
    PostProcessing,
    Output,
    Other,
}

impl StageKind {
    pub(super) const ALL: [Self; 8] = [
        Self::Setup,
        Self::DataLoad,
        Self::Ingress,
        Self::Lowering,
        Self::Native,
        Self::PostProcessing,
        Self::Output,
        Self::Other,
    ];

    #[cfg(feature = "progress-tui")]
    pub(super) fn label(self) -> &'static str {
        match self {
            Self::Setup => "setup",
            Self::DataLoad => "data load",
            Self::Ingress => "ingress",
            Self::Lowering => "scalar lowering",
            Self::Native => "native reduction",
            Self::PostProcessing => "post-processing",
            Self::Output => "output",
            Self::Other => "other",
        }
    }

    #[cfg(feature = "progress-tui")]
    pub(super) fn short_label(self) -> &'static str {
        match self {
            Self::Setup => "setup",
            Self::DataLoad => "load",
            Self::Ingress => "ingress",
            Self::Lowering => "lowering",
            Self::Native => "native",
            Self::PostProcessing => "post",
            Self::Output => "output",
            Self::Other => "other",
        }
    }
}

/// Stage identifier of every native activity (envelopes without a caller
/// phase). Activity changes inside the native core do not restart its clock.
pub(super) const NATIVE_STAGE: &str = "native";

/// The stage kind of a caller phase identifier. Identifiers are the static
/// `snake_case` names passed to the observer, never rendered text.
pub(super) fn classify_stage(id: &str) -> StageKind {
    match id {
        "" | "request_start" => StageKind::Setup,
        "data_loading"
        | "external_data_loading"
        | "loading"
        | "prepared_input_loading"
        | "prepared_input_selection"
        | "prepared_input_authentication_and_report"
        | "saved_application_preparation"
        | "routing_verification" => StageKind::DataLoad,
        "input_preparation"
        | "canonicalization"
        | "topology_matching"
        | "tensor_prepass"
        | "diagram_dressing"
        | "numerator_replacement"
        | "scalar_input"
        | "scalar_input_preparation"
        | "factored_input"
        | "matched_input_batching"
        | "integral_reduction"
        | "ordered_batch_preparation"
        | "coverage_mass_batching" => StageKind::Ingress,
        "scalar_lowering" | "lowering_guard_materialization" => StageKind::Lowering,
        NATIVE_STAGE | "entry_admission" | "routed_reduction" | "saved_rule_application" => {
            StageKind::Native
        }
        "durable_report_write"
        | "final_durable_report"
        | "report_serialization"
        | "guard_serialization"
        | "running_cache_close"
        | "incomplete_record_cleanup" => StageKind::Output,
        id if id.starts_with("result_") || id.ends_with("_serialization") => StageKind::Output,
        "terminal_substitution"
        | "terminal_geometry"
        | "unprojected_terminal_materialization"
        | "numerical_materialization"
        | "normalized_retained_output"
        | "lower_factor_normalization_unobserved"
        | "evaluation_output"
        | "evaluation_statistics"
        | "witness_cache_statistics"
        | "warm_exact_parity" => StageKind::PostProcessing,
        "retained_api_evaluation" => StageKind::Other,
        id if id.starts_with("retained_") => StageKind::PostProcessing,
        _ => StageKind::Other,
    }
}

/// `snake_case` identifiers and `CamelCase` native activities as words.
pub(super) fn humanize(id: &str) -> String {
    let mut words = String::with_capacity(id.len() + 4);
    for (index, character) in id.chars().enumerate() {
        if character == '_' {
            words.push(' ');
        } else if character.is_uppercase() {
            if index > 0 && !words.ends_with(' ') {
                words.push(' ');
            }
            words.extend(character.to_lowercase());
        } else {
            words.push(character);
        }
    }
    words
}

/// The current stage: a caller phase or the native core.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub(crate) struct StageModel {
    pub kind: StageKind,
    /// The caller phase identifier, or [`NATIVE_STAGE`].
    pub id: String,
    /// The native activity of the last envelope (native stage only).
    pub activity: Option<String>,
    pub seconds: f64,
    /// A caller-owned completed/total counter of this stage, when it has one.
    pub work: Option<Counter>,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub(crate) struct StageTotal {
    pub kind: StageKind,
    pub seconds: f64,
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(default)]
pub(crate) struct Counter {
    pub completed: usize,
    pub total: usize,
    pub label: String,
}

/// One Ordered phase kind as shown to people: the five-loop routed schedule
/// has Route and owner (Apply) phases; generic adapters only reduce.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum PhaseKind {
    #[default]
    Route,
    Apply,
    Reduce,
    Done,
}

impl PhaseKind {
    pub(super) fn of(
        kind: CandidateOrderedPhaseKind,
        owner: Option<usize>,
        completed: usize,
        total: usize,
    ) -> Self {
        match kind {
            _ if completed >= total => Self::Done,
            CandidateOrderedPhaseKind::Routing => Self::Route,
            CandidateOrderedPhaseKind::Reduction if owner.is_some() => Self::Apply,
            CandidateOrderedPhaseKind::Reduction => Self::Reduce,
        }
    }

    #[cfg(feature = "progress-tui")]
    pub(super) fn label(self) -> &'static str {
        match self {
            Self::Route => "Route",
            Self::Apply => "Apply",
            Self::Reduce => "Reduce",
            Self::Done => "done",
        }
    }
}

/// The phase running now (or the last one, once the schedule is done).
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub(crate) struct PhaseModel {
    /// Zero-based index of the running phase (= completed phases).
    pub index: usize,
    pub total: usize,
    pub kind: PhaseKind,
    pub support: Option<usize>,
    pub owner: Option<usize>,
    /// Time in this phase; absent until the phase start was observed.
    pub seconds: Option<f64>,
    /// Keys processed in this phase so far.
    pub keys: u64,
    pub step_cache_hits: u64,
    pub step_cache_misses: u64,
}

/// A finished phase of the current request.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub(crate) struct PhaseRecord {
    pub index: usize,
    pub kind: PhaseKind,
    pub support: Option<usize>,
    pub owner: Option<usize>,
    pub seconds: f64,
    pub keys: u64,
    pub step_cache_hits: u64,
    pub step_cache_misses: u64,
    /// Frontier keys and coefficient monomials when the phase ended.
    pub frontier_keys: u64,
    pub frontier_monomials: u64,
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(default)]
pub(crate) struct Census {
    pub keys: u64,
    pub peak_keys: u64,
    /// Numerator plus denominator monomials of the stored coefficients.
    pub monomials: u64,
    pub peak_monomials: u64,
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(default)]
pub(crate) struct Grouping {
    /// Nonzero checked contributions before frontier/output grouping.
    pub emitted: u64,
    /// Contributions merged into an existing frontier key.
    pub merged: u64,
    /// Frontier keys removed by exact cancellation.
    pub cancelled: u64,
    pub output_merged: u64,
    pub output_cancelled: u64,
}

/// Ordered engine census of the current native request.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub(crate) struct OrderedModel {
    /// Native request number within this monitor (phases restart per request).
    pub request: u64,
    pub phase: PhaseModel,
    pub input_keys: u64,
    pub processed_keys: u64,
    /// Keys per second over the recent window (about 30 s).
    pub recent_rate: Option<f64>,
    /// Keys per second over the whole request.
    pub request_rate: Option<f64>,
    pub frontier: Census,
    pub output: Census,
    pub grouping: Grouping,
    pub step_cache_hits: u64,
    pub step_cache_misses: u64,
    /// The most recent finished phases with work, oldest first.
    pub phases: Vec<PhaseRecord>,
    pub census_age_seconds: f64,
}

/// Depth-first engine counters (RustRed's DFS strategy; not used by Vakint's
/// own evaluations, which always reduce with the Ordered scheduler).
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub(crate) struct DfsModel {
    pub roots_completed: Option<u64>,
    pub roots_total: Option<u64>,
    pub entered: u64,
    pub completed: Option<u64>,
    pub attempts: Option<u64>,
    pub active: u64,
    pub pending_frames: u64,
    pub pending_unique: Option<u64>,
    pub pending_exact: bool,
    pub persistent_hits: Option<u64>,
    pub persistent_misses: Option<u64>,
    pub ram_hits: Option<u64>,
    pub cached_entries: Option<u64>,
    pub since_completion_seconds: Option<f64>,
    pub estimates: DfsEstimates,
    pub collection: Option<Collection>,
    pub path: Vec<String>,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub(crate) struct Collection {
    pub completed_roots: Option<u64>,
    pub total_roots: Option<u64>,
    pub contribution: Option<usize>,
    pub contributions: Option<usize>,
    pub age_seconds: Option<f64>,
}

/// Rate-based DFS estimates. They extrapolate observed throughput; future
/// discovery and caller post-processing are excluded.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub(crate) struct DfsEstimates {
    pub state: EstimateState,
    pub window_seconds: Option<f64>,
    pub samples: usize,
    pub sample_age_seconds: Option<f64>,
    pub completion_rate: Option<f64>,
    pub root_rate: Option<f64>,
    /// Net frontier growth (entered minus completed identities per second).
    pub growth_rate: Option<f64>,
    pub frontier_eta_seconds: Option<f64>,
    pub root_eta_seconds: Option<f64>,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum EstimateState {
    /// Not enough samples yet.
    #[default]
    Warmup,
    Live,
    /// The newest rate sample is older than its own window.
    Stale,
    /// No native work is running (caller stages are not estimated).
    Inactive,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub(crate) struct MemoryModel {
    /// Process values are only shown for a verified live worker identity.
    pub verified: bool,
    pub rss_bytes: Option<u64>,
    pub peak_rss_bytes: Option<u64>,
    pub cpu_percent: Option<f64>,
    pub available_bytes: Option<u64>,
    pub total_bytes: Option<u64>,
    pub sample_age_seconds: Option<f64>,
    pub guard: Option<GuardSummary>,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub(crate) struct GuardSummary {
    pub status: String,
    pub threshold_bytes: Option<u64>,
    pub covered: bool,
}

/// Secondary information for the details view.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub(crate) struct Details {
    pub event_age_seconds: Option<f64>,
    pub native_observation_age_seconds: Option<f64>,
    pub family: Option<String>,
    pub current_key: Option<String>,
    /// Canonical numerator rank R, positive power sum A and D = A - R.
    pub current_degrees: Option<[String; 3]>,
    pub caller_expression: Option<ExpressionSize>,
    pub worker_pid: Option<u32>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(default)]
pub(crate) struct ExpressionSize {
    pub bytes: usize,
    pub summands: Option<usize>,
    pub stage: String,
}

/// A family fingerprint can be several kilobytes; show its head and length.
pub(super) fn short_fingerprint(fingerprint: &str) -> Option<String> {
    const HEAD: usize = 40;
    if fingerprint.is_empty() {
        return None;
    }
    let mut characters = fingerprint.chars();
    let head = characters.by_ref().take(HEAD).collect::<String>();
    Some(if characters.next().is_none() {
        head
    } else {
        format!("{head}… ({} chars)", fingerprint.chars().count())
    })
}

impl MonitorModel {
    /// The model shown while the display is frozen: native census panels
    /// from `frozen`; status, clocks, stage, batch counters and memory live.
    #[cfg(any(feature = "progress-tui", test))]
    pub(super) fn with_frozen_native(mut self, frozen: &Self) -> Self {
        // Collection notices have their own clock and stay current.
        let collection = self.dfs.as_mut().and_then(|dfs| dfs.collection.take());
        self.ordered.clone_from(&frozen.ordered);
        self.dfs.clone_from(&frozen.dfs);
        if let Some(dfs) = &mut self.dfs {
            dfs.collection = collection;
        }
        self
    }

    /// Advance every running clock by `seconds` (time since the model was
    /// built, e.g. the publication age seen by a supervisor).
    #[cfg(any(feature = "progress-tui", test))]
    pub(super) fn advance(&mut self, seconds: f64) {
        if !(seconds.is_finite() && seconds > 0.0) {
            return;
        }
        let running = self.status == MonitorStatus::Running;
        if running {
            self.elapsed_seconds += seconds;
            self.stage.seconds += seconds;
            if let Some(total) = self
                .stages
                .iter_mut()
                .find(|total| total.kind == self.stage.kind)
            {
                total.seconds += seconds;
            }
            if let Some(phase) = self
                .ordered
                .as_mut()
                .and_then(|ordered| ordered.phase.seconds.as_mut())
            {
                *phase += seconds;
            }
        }
        if let Some(ordered) = &mut self.ordered {
            ordered.census_age_seconds += seconds;
        }
        for age in [
            &mut self.details.event_age_seconds,
            &mut self.details.native_observation_age_seconds,
        ]
        .into_iter()
        .flatten()
        {
            *age += seconds;
        }
        if let Some(dfs) = &mut self.dfs {
            for age in [
                &mut dfs.since_completion_seconds,
                &mut dfs.estimates.sample_age_seconds,
            ]
            .into_iter()
            .flatten()
            {
                *age += seconds;
            }
            if let Some(age) = dfs
                .collection
                .as_mut()
                .and_then(|collection| collection.age_seconds.as_mut())
            {
                *age += seconds;
            }
        }
    }
}

/// Compact human-readable formatting shared by the renderers.
#[cfg(any(feature = "progress-tui", test))]
pub(super) mod format {
    /// `9999`, `92.1k`, `26.2M`, `1.3G`.
    pub(in super::super) fn count(value: u64) -> String {
        let value_f = value as f64;
        match value {
            0..=9_999 => value.to_string(),
            10_000..=999_949 => format!("{:.1}k", value_f / 1e3),
            999_950..=999_949_999 => format!("{:.1}M", value_f / 1e6),
            _ => format!("{:.1}G", value_f / 1e9),
        }
    }

    /// `36.2s`, `5:50`, `4:14:24`.
    pub(in super::super) fn duration(seconds: f64) -> String {
        if !seconds.is_finite() || seconds < 0.0 {
            return "—".into();
        }
        if seconds < 60.0 {
            return format!("{seconds:.1}s");
        }
        let total = seconds as u64;
        let (hours, minutes, rest) = (total / 3600, (total / 60) % 60, total % 60);
        if hours == 0 {
            format!("{minutes}:{rest:02}")
        } else {
            format!("{hours}:{minutes:02}:{rest:02}")
        }
    }

    /// Per-second rates: `1.2k/s`, `263/s`, `9.7/s`.
    pub(in super::super) fn rate(value: f64) -> String {
        if !value.is_finite() || value < 0.0 {
            "—".into()
        } else if value >= 1000.0 {
            format!("{}/s", count(value.round() as u64))
        } else if value >= 100.0 {
            format!("{value:.0}/s")
        } else if value >= 0.1 {
            format!("{value:.1}/s")
        } else {
            format!("{value:.2}/s")
        }
    }

    /// IEC units, e.g. `9.9 GiB`.
    pub(in super::super) fn bytes(value: u64) -> String {
        let units = ["B", "KiB", "MiB", "GiB", "TiB", "PiB", "EiB"];
        let exponent = (value.max(1).ilog2() / 10) as usize;
        if exponent == 0 {
            format!("{value} B")
        } else {
            let scale = (1_u64 << (10 * exponent)) as f64;
            format!("{:.1} {}", value as f64 / scale, units[exponent])
        }
    }

    pub(in super::super) fn percent(part: u64, whole: u64) -> Option<String> {
        (whole > 0).then(|| format!("{:.1}%", 100.0 * part as f64 / whole as f64))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stage_identifiers_classify_without_display_text() {
        for (id, kind) in [
            ("request_start", StageKind::Setup),
            ("data_loading", StageKind::DataLoad),
            ("prepared_input_loading", StageKind::DataLoad),
            ("canonicalization", StageKind::Ingress),
            ("ordered_batch_preparation", StageKind::Ingress),
            ("scalar_lowering", StageKind::Lowering),
            (NATIVE_STAGE, StageKind::Native),
            ("saved_rule_application", StageKind::Native),
            ("terminal_geometry", StageKind::PostProcessing),
            (
                "retained_catalogue_exact_combination",
                StageKind::PostProcessing,
            ),
            ("result_expression_serialization", StageKind::Output),
            ("guard_serialization", StageKind::Output),
            ("retained_api_evaluation", StageKind::Other),
            ("synthetic_opaque_call", StageKind::Other),
        ] {
            assert_eq!(classify_stage(id), kind, "{id}");
        }
        assert_eq!(humanize("prepared_input_loading"), "prepared input loading");
        assert_eq!(
            humanize("CoefficientCombination"),
            "coefficient combination"
        );
    }

    #[test]
    fn compact_numbers_durations_and_rates() {
        use format::*;
        assert_eq!(count(9_999), "9999");
        assert_eq!(count(92_100), "92.1k");
        assert_eq!(count(999_960), "1.0M");
        assert_eq!(count(26_171_827), "26.2M");
        assert_eq!(count(3_000_000_000), "3.0G");
        assert_eq!(duration(36.24), "36.2s");
        assert_eq!(duration(350.0), "5:50");
        assert_eq!(duration(15_264.0), "4:14:24");
        assert_eq!(rate(9.7), "9.7/s");
        assert_eq!(rate(263.4), "263/s");
        assert_eq!(rate(1234.0), "1234/s");
        assert_eq!(rate(12_345.0), "12.3k/s");
        assert_eq!(bytes(1024), "1.0 KiB");
        assert_eq!(percent(1, 4).as_deref(), Some("25.0%"));
        assert_eq!(percent(1, 0), None);
    }

    #[test]
    fn long_fingerprints_are_shortened() {
        assert_eq!(short_fingerprint(""), None);
        assert_eq!(short_fingerprint("abc").as_deref(), Some("abc"));
        let long = "x".repeat(5104);
        let short = short_fingerprint(&long).unwrap();
        assert!(short.ends_with("(5104 chars)") && short.len() < 64);
    }

    #[test]
    fn advancing_moves_running_clocks_only() {
        let mut model = MonitorModel {
            status: MonitorStatus::Running,
            elapsed_seconds: 10.0,
            stage: StageModel {
                kind: StageKind::Native,
                seconds: 4.0,
                ..Default::default()
            },
            stages: vec![StageTotal {
                kind: StageKind::Native,
                seconds: 4.0,
            }],
            ordered: Some(OrderedModel {
                phase: PhaseModel {
                    seconds: Some(2.0),
                    ..Default::default()
                },
                census_age_seconds: 1.0,
                ..Default::default()
            }),
            ..Default::default()
        };
        model.advance(3.0);
        assert_eq!(model.elapsed_seconds, 13.0);
        assert_eq!(model.stage.seconds, 7.0);
        assert_eq!(model.stages[0].seconds, 7.0);
        let ordered = model.ordered.as_ref().unwrap();
        assert_eq!(ordered.phase.seconds, Some(5.0));
        assert_eq!(ordered.census_age_seconds, 4.0);
        model.status = MonitorStatus::Succeeded;
        model.advance(3.0);
        assert_eq!(model.elapsed_seconds, 13.0);
        assert_eq!(model.ordered.unwrap().census_age_seconds, 7.0);
    }
}
