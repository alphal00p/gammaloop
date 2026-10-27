//! Shared lossless capture through the existing benchmark logger, without another observer.
use super::{digest, validate_hit_diagnostic, write_json};
use eyre::{Result, ensure};
use serde_json::{Value, json};
use std::{
    collections::BTreeSet,
    path::Path,
    sync::atomic::{AtomicBool, Ordering},
};

#[derive(Default)]
pub(super) struct CacheHitDiagnostic {
    pub(super) input_roots: BTreeSet<Vec<i64>>,
    pub(super) input_roots_seen: bool,
    pub(super) input_roots_error: Option<String>,
    pub(super) input_dispatch_count: u64,
    pub(super) hit_events: Vec<Value>,
    pub(super) hit_error: Option<String>,
}

impl CacheHitDiagnostic {
    pub(super) const fn new() -> Self {
        Self {
            input_roots: BTreeSet::new(),
            input_roots_seen: false,
            input_roots_error: None,
            input_dispatch_count: 0,
            hit_events: Vec::new(),
            hit_error: None,
        }
    }
    pub(super) fn capture_input_roots(&mut self, message: &str) -> Result<()> {
        let message: Value = serde_json::from_str(message)?;
        ensure!(
            message["phase"] == "input_roots",
            "unknown input-root telemetry phase"
        );
        let roots: Vec<Vec<i64>> = serde_json::from_value(message["roots"].clone())?;
        ensure!(
            roots.iter().all(|root| root.len() == 15),
            "invalid input-root telemetry arity"
        );
        self.input_roots_seen = true;
        self.input_dispatch_count = self
            .input_dispatch_count
            .checked_add(1)
            .ok_or_else(|| eyre::eyre!("input dispatch counter overflow"))?;
        self.input_roots.extend(roots);
        Ok(())
    }

    pub(super) fn capture_hit_event(&mut self, message: &str) -> Result<()> {
        let event: Value = serde_json::from_str(message)?;
        ensure!(
            event["schema"] == "vakint.cache-hit-events.v1",
            "unknown cache-hit telemetry schema"
        );
        self.hit_events.try_reserve(1)?;
        self.hit_events.push(event);
        Ok(())
    }

    pub(super) fn hit_diagnostic(&self, success: bool) -> Value {
        json!({
            "schema":"vakint.cache-hit-diagnostic.v1",
            "interval_ns":0,"input_dispatch_count":self.input_dispatch_count,
            "events":self.hit_events,
            "error":self.hit_error.as_ref().or(self.input_roots_error.as_ref()),
            "success":success,
            "diagnostic_only":true,
            "timing_scope":"Lossless observation overhead; not an unobserved performance sample",
        })
    }

    pub(super) fn save_hit_diagnostic(&self, directory: &Path, success: bool) -> Value {
        let document = self.hit_diagnostic(success);
        let validation = validate_hit_diagnostic(&document);
        let mut receipt = json!({
            "path":"cache-hit-events.json","sha256":null,
            "diagnostic_only":true,"interval_ns":0,
            "stream_complete":validation.is_ok(),
            "classification_available":false,
            "classification_requires":"Same-cohort recursive mapper certificate and offline hit join",
            "error":validation.as_ref().err().map(ToString::to_string),
        });
        let path = directory.join("cache-hit-events.json");
        // Evidence is published only after the full buffered stream has been
        // written, synced and atomically renamed. Failure never invents counts
        // or replaces the mathematical outcome with a telemetry error.
        match write_json(&path, &document).and_then(|()| digest(&path)) {
            Ok(hash) => receipt["sha256"] = hash.into(),
            Err(error) => {
                receipt["stream_complete"] = false.into();
                receipt["error"] =
                    format!("cache-hit evidence could not be committed: {error}").into();
            }
        }
        receipt
    }
}

static CLASSIFY_CACHE_HITS: AtomicBool = AtomicBool::new(false);
pub(super) struct ClassificationLogging;

impl ClassificationLogging {
    pub(super) fn enabled() -> bool {
        CLASSIFY_CACHE_HITS.load(Ordering::Relaxed)
    }
    pub(super) fn enter(enabled: bool) -> Self {
        CLASSIFY_CACHE_HITS.store(enabled, Ordering::Relaxed);
        Self
    }
}

impl Drop for ClassificationLogging {
    fn drop(&mut self) {
        CLASSIFY_CACHE_HITS.store(false, Ordering::Relaxed);
    }
}
