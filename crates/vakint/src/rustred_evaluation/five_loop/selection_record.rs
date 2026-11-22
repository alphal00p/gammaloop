//! Opt-in record of the five-loop witness selection (diagnostics only).
//!
//! `VAKINT_WITNESS_LOG=<file>` appends one JSON line per matched input
//! (`event: "input"`, with the integral and the numerator), per full weighted
//! lowering of a witness (`event: "full_key"`: the progress counter of the
//! selection loop) and per selection (`event: "selection"`: how many witnesses
//! the exact screen decided, how many were fully lowered, the selected witness
//! and its exact support). Unset, nothing is formatted or written, and the
//! selection is unchanged either way.

use std::io::Write;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Mutex, OnceLock};

use sha2::{Digest, Sha256};

use crate::ReplacementRules;

static SINK: OnceLock<Option<Mutex<std::fs::File>>> = OnceLock::new();
static SEQUENCE: AtomicU64 = AtomicU64::new(0);

fn sink() -> Option<&'static Mutex<std::fs::File>> {
    SINK.get_or_init(|| {
        let path = std::env::var_os("VAKINT_WITNESS_LOG")?;
        std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(path)
            .ok()
            .map(Mutex::new)
    })
    .as_ref()
}

/// Whether the record is on (`VAKINT_WITNESS_LOG` names a writable file).
pub(super) fn enabled() -> bool {
    sink().is_some()
}

/// Process-wide sequence number of a selection, in start order.
pub(super) fn next_sequence() -> u64 {
    SEQUENCE.fetch_add(1, Ordering::Relaxed)
}

/// Append one record (an object) with the process id and the wall clock.
pub(super) fn write(mut record: serde_json::Value) {
    let Some(sink) = sink() else {
        return;
    };
    if let Some(object) = record.as_object_mut() {
        object.insert("pid".into(), std::process::id().into());
        let now = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0.0, |elapsed| elapsed.as_secs_f64());
        object.insert("unix_seconds".into(), now.into());
    }
    let line = format!("{record}\n");
    if let Ok(mut file) = sink.lock() {
        let _ = file.write_all(line.as_bytes());
        let _ = file.flush();
    }
}

/// First 16 hex digits of the sha256 of `text`.
pub(super) fn digest(text: &str) -> String {
    Sha256::digest(text.as_bytes())
        .iter()
        .take(8)
        .map(|byte| format!("{byte:02x}"))
        .collect()
}

/// Order-independent text of a witness: its sorted numerator substitutions.
pub(super) fn witness_text(witness: &ReplacementRules) -> String {
    let mut rules = witness
        .numerator_substitutions
        .iter()
        .map(|(source, target)| format!("{source}->{target}"))
        .collect::<Vec<_>>();
    rules.sort();
    rules.join(";")
}
