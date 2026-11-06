//! In-memory latest model of one monitored call, for in-process consumers.

use std::sync::{
    Arc, Mutex, TryLockError,
    atomic::{AtomicU64, Ordering},
};

/// The latest `vakint.monitor.v2` model of a monitored call, kept in memory.
///
/// Set [`RustRedMonitorOptions::progress_sink`](super::RustRedMonitorOptions::progress_sink)
/// to receive the same model as the `model` object of a progress snapshot
/// file, at the same cadence (the monitor interval, between 0.1 s and 1 s,
/// plus request start, finish and admission), without a file or terminal.
/// After the call returns the sink holds the call's final model. Clones share
/// the model. The monitor thread is the only writer; readers never block it.
///
/// ```
/// use vakint::{RustRedMonitorOptions, RustRedProgressSink};
///
/// let sink = RustRedProgressSink::default();
/// let options = RustRedMonitorOptions {
///     progress_sink: Some(sink.clone()),
///     ..Default::default()
/// };
/// assert!(options.enabled());
/// assert!(sink.latest_since(0).is_none());
/// ```
#[derive(Clone, Debug, Default)]
pub struct RustRedProgressSink(Arc<SinkState>);

#[derive(Debug, Default)]
struct SinkState {
    /// Number of stored models; read without the lock.
    sequence: AtomicU64,
    latest: Mutex<Option<serde_json::Value>>,
}

impl RustRedProgressSink {
    /// Replace the latest model.
    pub(crate) fn store(&self, model: serde_json::Value) {
        let mut latest = self
            .0
            .latest
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        *latest = Some(model);
        // Under the lock, so a reader sees a sequence with its own model.
        self.0.sequence.fetch_add(1, Ordering::Release);
    }

    /// Number of models stored so far.
    pub fn sequence(&self) -> u64 {
        self.0.sequence.load(Ordering::Acquire)
    }

    /// The latest model with its sequence number, if newer than `seen`.
    ///
    /// Returns `None` while the monitor is storing a model; ask again later.
    pub fn latest_since(&self, seen: u64) -> Option<(u64, serde_json::Value)> {
        if self.sequence() <= seen {
            return None;
        }
        let latest = match self.0.latest.try_lock() {
            Ok(latest) => latest,
            Err(TryLockError::Poisoned(error)) => error.into_inner(),
            Err(TryLockError::WouldBlock) => return None,
        };
        let sequence = self.sequence();
        latest.clone().map(|model| (sequence, model))
    }
}

#[cfg(feature = "symbolica_community_module")]
pub(crate) use metrics::ProgressMetrics;

/// The simplified progress view of the Python binding.
#[cfg(feature = "symbolica_community_module")]
mod metrics {
    use super::super::{
        MonitorStatus,
        model::{self, MonitorModel, format},
    };

    /// The most important progress metrics of one model, flattened.
    #[derive(Clone, Debug, Default, PartialEq)]
    pub(crate) struct ProgressMetrics {
        /// `running`, `succeeded`, `failed` or `partial`; the caller may
        /// override it (e.g. `interrupted`).
        pub status: String,
        pub elapsed_seconds: f64,
        /// The stage kind in `snake_case` (`setup`, `data_load`, `native`, ...).
        pub stage: String,
        /// The stage identifier (a caller phase, or `native`).
        pub stage_id: String,
        /// The native activity, else the caller stage, in words.
        pub activity: String,
        pub batch: Option<(usize, usize, String)>,
        /// Ordered schedule position: running phase index and phase total.
        pub phase: Option<(usize, usize)>,
        pub phase_kind: Option<String>,
        pub processed_keys: Option<u64>,
        pub input_keys: Option<u64>,
        pub keys_per_second: Option<f64>,
        pub frontier_keys: Option<u64>,
        pub peak_frontier_keys: Option<u64>,
        pub step_cache_hits: Option<u64>,
        pub step_cache_misses: Option<u64>,
        pub rss_bytes: Option<u64>,
        pub peak_rss_bytes: Option<u64>,
    }

    fn snake_case<T: serde::Serialize>(value: T) -> String {
        serde_json::to_value(value)
            .ok()
            .and_then(|value| value.as_str().map(str::to_owned))
            .unwrap_or_default()
    }

    impl ProgressMetrics {
        /// Read a published model (the `model` object of `vakint.monitor.v2`).
        pub(crate) fn from_model(model: &serde_json::Value) -> Option<Self> {
            let model: MonitorModel = serde_json::from_value(model.clone()).ok()?;
            let status = match model.status {
                MonitorStatus::Running => "running",
                MonitorStatus::Succeeded => "succeeded",
                MonitorStatus::Failed => "failed",
                MonitorStatus::Partial => "partial",
            };
            let activity = model
                .stage
                .activity
                .clone()
                .filter(|activity| !activity.is_empty())
                .unwrap_or_else(|| match model.stage.id.as_str() {
                    "" | "request_start" => "starting".into(),
                    id => model::humanize(id),
                });
            let ordered = model.ordered.as_ref();
            Some(Self {
                status: status.into(),
                elapsed_seconds: model.elapsed_seconds,
                stage: snake_case(model.stage.kind),
                stage_id: model.stage.id.clone(),
                activity,
                batch: model
                    .batch
                    .as_ref()
                    .map(|batch| (batch.completed, batch.total, batch.label.clone())),
                phase: ordered.map(|ordered| (ordered.phase.index, ordered.phase.total)),
                phase_kind: ordered.map(|ordered| snake_case(ordered.phase.kind)),
                processed_keys: ordered.map(|ordered| ordered.processed_keys),
                input_keys: ordered.map(|ordered| ordered.input_keys),
                keys_per_second: ordered
                    .and_then(|ordered| ordered.recent_rate.or(ordered.request_rate))
                    .map(|rate| rate.keys as f64 / rate.seconds),
                frontier_keys: ordered.map(|ordered| ordered.frontier.keys),
                peak_frontier_keys: ordered.map(|ordered| ordered.frontier.peak_keys),
                step_cache_hits: ordered.map(|ordered| ordered.step_cache_hits),
                step_cache_misses: ordered.map(|ordered| ordered.step_cache_misses),
                rss_bytes: model.memory.rss_bytes,
                peak_rss_bytes: model.memory.peak_rss_bytes,
            })
        }

        /// A caller-only view (no native monitor) of stage `stage_id` (a
        /// caller phase identifier, e.g. `external_data_loading`).
        pub(crate) fn caller(stage_id: &str, activity: &str, elapsed_seconds: f64) -> Self {
            Self {
                status: "running".into(),
                elapsed_seconds,
                stage: snake_case(model::classify_stage(stage_id)),
                stage_id: stage_id.into(),
                activity: activity.into(),
                ..Default::default()
            }
        }

        /// The model of [`Self::caller`], in the published schema.
        pub(crate) fn caller_model(&self) -> serde_json::Value {
            let mut model = MonitorModel {
                elapsed_seconds: self.elapsed_seconds,
                ..Default::default()
            };
            model.stage.kind = model::classify_stage(&self.stage_id);
            model.stage.id = self.stage_id.clone();
            model.stage.activity = Some(self.activity.clone());
            model.stage.seconds = self.elapsed_seconds;
            serde_json::to_value(model).unwrap_or_default()
        }

        /// One line, absent parts omitted, e.g. `Apply 38/83 · 6.3k keys/s ·
        /// frontier 12.3k (peak 174.9k) · cache 0/361.4k hits · RSS 12.8 GiB
        /// (peak 15.0 GiB) · 2:34`.
        pub(crate) fn summary(&self) -> String {
            let mut parts = Vec::new();
            match (self.phase, self.phase_kind.as_deref()) {
                (Some((_, total)), Some("done")) => {
                    parts.push(format!("Ordered schedule done ({total} phases)"));
                }
                (Some((index, total)), Some(kind)) => {
                    let label = match kind {
                        "route" => "Route",
                        "apply" => "Apply",
                        _ => "Reduce",
                    };
                    parts.push(format!("{label} {}/{total}", index + 1));
                }
                _ => parts.push(self.activity.clone()),
            }
            if let Some((completed, total, label)) = &self.batch
                && *total > 1
            {
                let label = if label.is_empty() { "batch" } else { label };
                parts.push(format!("{label} {completed}/{total}"));
            }
            if let Some(rate) = self.keys_per_second {
                parts.push(format::rate(rate).replacen("/s", " keys/s", 1));
            }
            if let Some(frontier) = self.frontier_keys {
                let mut text = format!("frontier {}", format::count(frontier));
                if let Some(peak) = self.peak_frontier_keys.filter(|peak| *peak > frontier) {
                    text.push_str(&format!(" (peak {})", format::count(peak)));
                }
                parts.push(text);
            }
            if let (Some(hits), Some(misses)) = (self.step_cache_hits, self.step_cache_misses)
                && hits + misses > 0
            {
                parts.push(format!(
                    "cache {}/{} hits",
                    format::count(hits),
                    format::count(hits + misses)
                ));
            }
            if let Some(rss) = self.rss_bytes {
                let mut text = format!("RSS {}", format::bytes(rss));
                if let Some(peak) = self.peak_rss_bytes.filter(|peak| *peak > rss) {
                    text.push_str(&format!(" (peak {})", format::bytes(peak)));
                }
                parts.push(text);
            }
            parts.push(format::duration(self.elapsed_seconds));
            parts.join(" · ")
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn metrics_flatten_the_ordered_census_and_summarize_it() {
            let mut model = MonitorModel {
                status: MonitorStatus::Running,
                elapsed_seconds: 154.2,
                ..Default::default()
            };
            model.stage.kind = model::StageKind::Native;
            model.stage.id = "native".into();
            model.stage.activity = Some("ordered reduction".into());
            let mut ordered = model::OrderedModel {
                input_keys: 79,
                processed_keys: 361_435,
                step_cache_misses: 361_435,
                recent_rate: Some(model::Rate {
                    keys: 63_120,
                    seconds: 10.0,
                }),
                ..Default::default()
            };
            ordered.phase.index = 37;
            ordered.phase.total = 83;
            ordered.phase.kind = model::PhaseKind::Apply;
            ordered.frontier.keys = 12_345;
            ordered.frontier.peak_keys = 174_856;
            model.ordered = Some(ordered);
            model.memory.rss_bytes = Some(12_800 << 20);
            model.memory.peak_rss_bytes = Some(15_000 << 20);
            let metrics = ProgressMetrics::from_model(&serde_json::to_value(&model).unwrap())
                .expect("a published model");
            assert_eq!(metrics.status, "running");
            assert_eq!(metrics.stage, "native");
            assert_eq!(metrics.activity, "ordered reduction");
            assert_eq!(metrics.phase, Some((37, 83)));
            assert_eq!(metrics.phase_kind.as_deref(), Some("apply"));
            assert_eq!(metrics.keys_per_second, Some(6312.0));
            assert_eq!(
                metrics.summary(),
                "Apply 38/83 · 6312 keys/s · frontier 12.3k (peak 174.9k) · \
                 cache 0/361.4k hits · RSS 12.5 GiB (peak 14.6 GiB) · 2:34"
            );
            // Caller stages without a census: the stage in words and the clock.
            let caller =
                ProgressMetrics::caller("external_data_loading", "loading 1L-4L data", 3.31);
            assert_eq!(caller.stage, "data_load");
            assert_eq!(caller.summary(), "loading 1L-4L data · 3.3s");
            let model = caller.caller_model();
            assert_eq!(model["stage"]["kind"], "data_load");
            assert_eq!(model["status"], "Running");
            let back = ProgressMetrics::from_model(&model).unwrap();
            assert_eq!(back.stage, "data_load");
            assert_eq!(back.activity, "loading 1L-4L data");
            assert!(ProgressMetrics::from_model(&serde_json::json!([1])).is_none());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn readers_see_each_stored_model_once_with_its_sequence() {
        let sink = RustRedProgressSink::default();
        assert_eq!(sink.sequence(), 0);
        assert!(sink.latest_since(0).is_none());
        sink.store(serde_json::json!({"elapsed_seconds": 1.0}));
        let (seen, model) = sink.latest_since(0).unwrap();
        assert_eq!((seen, model["elapsed_seconds"].as_f64()), (1, Some(1.0)));
        assert!(sink.latest_since(seen).is_none());
        let clone = sink.clone();
        clone.store(serde_json::json!({"elapsed_seconds": 2.0}));
        assert_eq!(sink.latest_since(seen).unwrap().0, 2);
        // A reader never blocks a held lock.
        let guard = sink.0.latest.lock().unwrap();
        assert!(clone.latest_since(0).is_none());
        drop(guard);
        assert!(clone.latest_since(0).is_some());
    }
}
