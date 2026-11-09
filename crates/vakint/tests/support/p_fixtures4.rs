//! Shared helpers for the four-loop Forcer fixtures of the M3 plan
//! (`tests/fixtures/forcer/p4l_*.json`, M3-V7).
//!
//! Include it after the general fixture helpers:
//!
//! ```text
//! #[path = "support/p_fixtures.rs"]
//! mod p_fixtures;
//! #[path = "support/p_fixtures4.rs"]
//! mod p_fixtures4;
//! ```
//!
//! A file is regenerated with FORM 5 by `VAKINT_FORCER_FIXTURES=write`
//! (`p_fixtures::writer`), with `VAKINT_FORCER_FIXTURES_WORKERS` TFORM
//! workers per FORM process (default 1; the M3 runs use 4, `tform -w4`).

#![allow(dead_code)]

use vakint::p_integral::forcer_fixtures::{
    ANCHOR4_PROVENANCE, ForcerCase, ForcerFixture, anchor_cases_4l, write_or_read_with,
};

use super::p_fixtures::{fixture_path, writer};

/// TFORM workers per FORM process of a regeneration.
pub fn workers() -> usize {
    std::env::var("VAKINT_FORCER_FIXTURES_WORKERS")
        .ok()
        .and_then(|workers| workers.parse().ok())
        .unwrap_or(1)
}

fn written_or_read(file: &str, provenance: &str, definitions: Vec<ForcerCase>) -> ForcerFixture {
    let writer = writer();
    write_or_read_with(
        &fixture_path(file),
        provenance,
        definitions,
        writer
            .as_ref()
            .map(|(form, jobs)| (form.as_str(), *jobs, workers())),
    )
    .unwrap_or_else(|error| panic!("{file}: {error}"))
}

/// `p4l_anchors.json` (regenerated first in writer mode), checked against
/// its definition (`anchor_cases_4l`).
pub fn anchor4_fixture() -> ForcerFixture {
    written_or_read("p4l_anchors.json", ANCHOR4_PROVENANCE, anchor_cases_4l())
}
