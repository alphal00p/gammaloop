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
//! `VAKINT_FORCER_FIXTURES=write` regenerates a file with FORM 5
//! (`p_fixtures::writer`), with `VAKINT_FORCER_FIXTURES_WORKERS` TFORM
//! workers per FORM process (default 1; the M3 runs use 4, `tform -w4`), and
//! then its depths; `VAKINT_FORCER_FIXTURES=depths` recomputes only the
//! depths (and the 11-term values) from the frozen outputs, without FORM, on
//! `VAKINT_FORCER_FIXTURES_WORKERS` threads.

#![allow(dead_code)]

use symbolica::atom::Atom;
use vakint::p_integral::forcer_fixtures::{
    ANCHOR4_PROVENANCE, BATTERY4_PROVENANCE, ForcerCase, ForcerFixture, GRID4_PROVENANCE,
    anchor_cases_4l, grid_cases_4l, write_or_read_with,
};

use super::p_fixtures::{benchmark_inputs_text, fixture_path, writer};

/// A stored value (Symbolica's canonical text) as an atom.
pub fn parse(text: &str) -> Atom {
    vakint::vakint_parse!(text).unwrap_or_else(|error| panic!("{text}: {error}"))
}

/// TFORM workers per FORM process of a regeneration, and threads of a depth
/// pass.
pub fn workers() -> usize {
    std::env::var("VAKINT_FORCER_FIXTURES_WORKERS")
        .ok()
        .and_then(|workers| workers.parse().ok())
        .unwrap_or(1)
}

/// Whether the depths are (re)computed: `VAKINT_FORCER_FIXTURES=write` or
/// `depths`.
fn depths_requested() -> bool {
    std::env::var("VAKINT_FORCER_FIXTURES").is_ok_and(|mode| mode == "write" || mode == "depths")
}

/// `fixture` with its depths (and, with `eleven`, its 11-term values)
/// recomputed and saved to `file` when [`depths_requested`].
fn with_depths(file: &str, mut fixture: ForcerFixture, eleven: bool) -> ForcerFixture {
    if depths_requested() {
        let started = std::time::Instant::now();
        fixture
            .write_depths(eleven, workers())
            .unwrap_or_else(|error| panic!("{file}: {error}"));
        std::fs::write(fixture_path(file), fixture.to_json()).unwrap();
        eprintln!("{file}: depths in {:.1} s", started.elapsed().as_secs_f64());
    }
    fixture
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

/// `p4l_grid.json` (regenerated first in writer mode), checked against its
/// definition (`grid_cases_4l`).
pub fn grid4_fixture() -> ForcerFixture {
    let definitions = grid_cases_4l(&benchmark_inputs_text()).unwrap();
    let fixture = written_or_read("p4l_grid.json", GRID4_PROVENANCE, definitions);
    with_depths("p4l_grid.json", fixture, false)
}

/// `p4l_battery.json`, whose cases are defined by the file itself
/// (regenerated first in writer mode, outputs and depths).
pub fn battery4_fixture() -> ForcerFixture {
    let file = "p4l_battery.json";
    let read = || {
        ForcerFixture::parse(&std::fs::read_to_string(fixture_path(file)).unwrap())
            .unwrap_or_else(|error| panic!("{file}: {error}"))
    };
    let mut fixture = read();
    if let Some((form, jobs)) = writer() {
        fixture = fixture.unwritten_like(BATTERY4_PROVENANCE);
        fixture
            .write_outputs_with(&form, jobs, workers())
            .unwrap_or_else(|error| panic!("{file}: {error}"));
        std::fs::write(fixture_path(file), fixture.to_json()).unwrap();
    }
    with_depths(file, fixture, true);
    let fixture = read();
    assert_eq!(
        fixture.provenance.as_deref(),
        Some(BATTERY4_PROVENANCE),
        "{file}: provenance"
    );
    fixture
}
