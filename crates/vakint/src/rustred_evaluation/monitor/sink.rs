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
