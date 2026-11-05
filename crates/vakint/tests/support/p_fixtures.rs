//! Shared helpers for the frozen Forcer fixtures (`tests/fixtures/forcer/`)
//! and the benchmark inputs (`tests/fixtures/benchmark_p_integrals.json`).
//!
//! Include it with `#[path = "support/p_fixtures.rs"] mod p_fixtures;`.
//! The fixture types live in `vakint::p_integral::forcer_fixtures`; a file is
//! regenerated with FORM 5 by `VAKINT_FORCER_FIXTURES=write` (FORM from
//! `VAKINT_TEST_FORM5`, else `FORM_PATH`, else `form`; parallel FORM processes
//! from `VAKINT_FORCER_FIXTURES_JOBS`, default 1).

#![allow(dead_code)]

use std::path::PathBuf;

use symbolica::atom::Atom;
use vakint::p_integral::forcer_fixtures::{
    ANCHOR_PROVENANCE, BenchmarkInputs, ForcerCase, ForcerFixture, GRID_PROVENANCE, anchor_cases,
    grid_cases, write_or_read,
};
use vakint::{RustRedMasterMode, VakintDependency, VakintSettings};

/// `tests/fixtures/forcer/<file>`.
pub fn fixture_path(file: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/forcer")
        .join(file)
}

/// The text of `tests/fixtures/benchmark_p_integrals.json`.
pub fn benchmark_inputs_text() -> String {
    std::fs::read_to_string(
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/benchmark_p_integrals.json"),
    )
    .unwrap()
}

/// The parsed benchmark inputs.
pub fn benchmark_inputs() -> BenchmarkInputs {
    BenchmarkInputs::parse(&benchmark_inputs_text()).unwrap()
}

/// FORM 5 for live Forcer runs, if there is one.
pub fn form5() -> Option<String> {
    std::env::var("VAKINT_TEST_FORM5").ok().or_else(|| {
        let form = std::env::var("FORM_PATH").unwrap_or_else(|_| "form".into());
        VakintDependency::Forcer
            .version(&form, "python3")
            .is_ok()
            .then_some(form)
    })
}

/// `Some((form, jobs))` when `VAKINT_FORCER_FIXTURES=write` asks for a
/// regeneration; it then needs FORM 5.
pub fn writer() -> Option<(String, usize)> {
    if std::env::var("VAKINT_FORCER_FIXTURES").is_ok_and(|mode| mode == "write") {
        let form = form5().expect("regenerating Forcer fixtures needs FORM 5");
        let jobs = std::env::var("VAKINT_FORCER_FIXTURES_JOBS")
            .ok()
            .and_then(|jobs| jobs.parse().ok())
            .unwrap_or(1);
        Some((form, jobs))
    } else {
        None
    }
}

fn written_or_read(file: &str, provenance: &str, definitions: Vec<ForcerCase>) -> ForcerFixture {
    let writer = writer();
    write_or_read(
        &fixture_path(file),
        provenance,
        definitions,
        writer.as_ref().map(|(form, jobs)| (form.as_str(), *jobs)),
    )
    .unwrap_or_else(|error| panic!("{file}: {error}"))
}

/// `p{loops}l_grid.json` (regenerated first in writer mode), checked against
/// its definition.
pub fn grid_fixture(loops: usize) -> ForcerFixture {
    let definitions = grid_cases(loops, &benchmark_inputs_text()).unwrap();
    written_or_read(&format!("p{loops}l_grid.json"), GRID_PROVENANCE, definitions)
}

/// `p_anchors.json` (regenerated first in writer mode), checked against its
/// definition.
pub fn anchor_fixture() -> ForcerFixture {
    written_or_read("p_anchors.json", ANCHOR_PROVENANCE, anchor_cases())
}

/// One of the M1 batteries `p{1,2,3,4}l_battery.json`.
pub fn battery(file: &str) -> ForcerFixture {
    ForcerFixture::parse(&std::fs::read_to_string(fixture_path(file)).unwrap()).unwrap()
}

/// The case `name` of `fixture`.
pub fn case<'a>(fixture: &'a ForcerFixture, name: &str) -> &'a ForcerCase {
    fixture
        .cases
        .iter()
        .find(|case| case.name == name)
        .unwrap_or_else(|| panic!("no fixture case {name}"))
}

/// The frozen `(scale_power, reduced)` values of a case in vakint's `ε`.
pub fn frozen_values(case: &ForcerCase) -> Vec<(i64, Atom)> {
    case.frozen_values(&VakintSettings::default().epsilon_symbol)
        .unwrap_or_else(|error| panic!("{}: {error}", case.name))
}

/// The frozen value of a case finalized like a Forcer evaluation.
pub fn finalized(case: &ForcerCase, settings: &VakintSettings, mode: RustRedMasterMode) -> Atom {
    case.finalized(settings, mode)
        .unwrap_or_else(|error| panic!("{} ({mode}): {error}", case.name))
}
