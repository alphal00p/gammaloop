//! Tier C of the Forcer battery (M2-V5; M2 plan, sections 8 and 9): the full
//! frozen grids and the benchmark cells, RustRed against Forcer through
//! vakint's public API (`support/p_backends.rs`: `Retain` exactly in the tag
//! basis, `Symbolic` exactly, `Numerical` at 40 digits, through `ε^(L+2)`),
//! and the exact-agreement file of the benchmark.
//!
//! Full profile only (`just test_vakint_full`): `.config/nextest.toml`
//! excludes this binary from `test_vakint` by name.
//!
//! - Grids: every case of `tests/fixtures/forcer/p{1,2,3}l_grid.json` (17,
//!   39, 417) against its frozen Forcer output, in parts of at most
//!   [`GRID_PART`] cases (one test each).
//! - Benchmark cells: every 1L-3L cell of
//!   `tests/fixtures/benchmark_p_integrals.json`, its frozen input (what the
//!   benchmark times) evaluated by RustRed against the sum of the frozen
//!   Forcer outputs of its tops (`bench_<cell>_<family>` in the grids).
//! - Physics cells: the 2L and 3L inputs of
//!   `tests/fixtures/benchmark_p_physics.json` against a live Forcer (FORM 5).
//! - With `VAKINT_P_EXACT_OUT=<path>` (never set by default),
//!   [`the_exact_agreement_file_is_written_on_request`] computes every 1L-3L
//!   cell and both physics cells and writes `<path>` (schema
//!   `vakint.p-exact.v1`, read by `examples/benchmark_p_integrals` and
//!   `docs/benchmark/benchmark.typ`); it refuses a dirty tree and writes
//!   nothing unless every cell was computed. M2-T10c runs it with
//!   `VAKINT_P_EXACT_OUT=crates/vakint/docs/benchmark/data/p_exact.json` at the
//!   source of the measurement, e.g.
//!   `cargo test --profile ci-optim -p vakint --features vakint/feynkit-ingress
//!   --test p_integral_rustred_full -- --exact the_exact_agreement_file_is_written_on_request`
//!   (one process, no per-test time limit).
//!
//! STAGE 1 (until M2-V4 ships the built-in 1L-3L library) checks the one-loop
//! parts against RustRed's one-loop reference library; the two- and
//! three-loop parts are ignored and the file writer refuses to run.

#![cfg(feature = "feynkit-ingress")]

#[path = "support/p_integrals.rs"]
mod p_integrals;
#[path = "support/p_rustred.rs"]
mod p_rustred;
#[path = "support/p_backends.rs"]
mod p_backends;
#[path = "support/p_fixtures.rs"]
mod p_fixtures;

use std::path::{Path, PathBuf};
use std::process::Command;

use p_backends::*;
use p_fixtures::{benchmark_inputs, grid_fixture};
use p_integrals::{forcer_form, parse};
use serde_json::{Value, json};
use symbolica::atom::{Atom, AtomCore};
use vakint::p_integral::forcer_fixtures::{BenchmarkCell, ForcerCase};
use vakint::utils::simplify_real;
use vakint::{LoopNormalizationFactor, RustRedMasterMode};

/// At most this many grid cases per test (the full profile's time limit is
/// per test).
const GRID_PART: usize = 40;

/// The Forcer side of a frozen case in `mode`: the frozen output finalized
/// like a Forcer evaluation; in `Numerical` mode the exact (`Symbolic`) value,
/// which [`values_agree`] evaluates at 40 digits.
fn frozen_value(case: &ForcerCase, mode: RustRedMasterMode, terms: i64) -> Result<Atom, String> {
    let mode = match mode {
        RustRedMasterMode::Numerical => RustRedMasterMode::Symbolic,
        mode => mode,
    };
    let settings = rustred_p_settings(mode, LoopNormalizationFactor::GScheme, terms);
    case.finalized(&settings, mode)
        .map(|value| simplify_real(value.as_view()))
        .map_err(|error| format!("Forcer fixture {}: {error}", case.name))
}

/// The sum of [`frozen_value`] over `cases`.
fn frozen_sum(cases: &[&ForcerCase], mode: RustRedMasterMode, terms: i64) -> Result<Atom, String> {
    cases.iter().try_fold(Atom::Zero, |sum, case| {
        Ok(sum + frozen_value(case, mode, terms)?)
    })
}

/// Part `part` of the `loops`-loop grid: every case agrees with its frozen
/// Forcer output; the parts cover the grid exactly once.
fn assert_grid_part(loops: usize, part: usize, parts: usize) {
    let fixture = grid_fixture(loops);
    assert_eq!(
        fixture.cases.len().div_ceil(GRID_PART),
        parts,
        "{loops}-loop grid parts"
    );
    let terms = battery_terms(loops);
    for case in fixture.cases.iter().skip(part * GRID_PART).take(GRID_PART) {
        let input = case.expression().unwrap();
        agreement(&input, terms, |mode| frozen_value(case, mode, terms))
            .unwrap_or_else(|difference| panic!("p{loops}l_grid.json: {}: {difference}", case.name));
    }
}

#[test]
fn the_one_loop_grid_agrees_with_forcer() {
    assert_grid_part(1, 0, 1);
}

#[test]
#[ignore = "M2-V5 STAGE 2: needs the built-in 1L-3L p library (M2-V4)"]
fn the_two_loop_grid_agrees_with_forcer() {
    assert_grid_part(2, 0, 1);
}

/// The 417 three-loop cases in 11 parts.
macro_rules! three_loop_grid_parts {
    ($($name:ident: $part:expr),* $(,)?) => {$(
        #[test]
        #[ignore = "M2-V5 STAGE 2: needs the built-in 1L-3L p library (M2-V4)"]
        fn $name() {
            assert_grid_part(3, $part, 11);
        }
    )*};
}

three_loop_grid_parts!(
    the_three_loop_grid_part_01_agrees_with_forcer: 0,
    the_three_loop_grid_part_02_agrees_with_forcer: 1,
    the_three_loop_grid_part_03_agrees_with_forcer: 2,
    the_three_loop_grid_part_04_agrees_with_forcer: 3,
    the_three_loop_grid_part_05_agrees_with_forcer: 4,
    the_three_loop_grid_part_06_agrees_with_forcer: 5,
    the_three_loop_grid_part_07_agrees_with_forcer: 6,
    the_three_loop_grid_part_08_agrees_with_forcer: 7,
    the_three_loop_grid_part_09_agrees_with_forcer: 8,
    the_three_loop_grid_part_10_agrees_with_forcer: 9,
    the_three_loop_grid_part_11_agrees_with_forcer: 10,
);

/// One cell of the exact-agreement file: a benchmark cell or a physics cell.
#[derive(Clone, Debug, PartialEq)]
struct ExactCell {
    loops: usize,
    /// `(dots, rank)` of a benchmark cell, `None` for the physics cell.
    grid: Option<(usize, usize)>,
    agree: bool,
    terms: i64,
}

impl ExactCell {
    fn to_json(&self) -> Value {
        match self.grid {
            Some((dots, rank)) => json!({
                "L": self.loops, "dots": dots, "rank": rank, "agree": self.agree, "terms": self.terms
            }),
            None => json!({
                "L": self.loops, "physics": true, "agree": self.agree, "terms": self.terms
            }),
        }
    }
}

/// The agreement of a benchmark cell: its frozen input with RustRed against
/// the frozen Forcer outputs of its tops. A disagreement is reported on
/// stderr and recorded as `agree: false`.
fn benchmark_cell(cell: &BenchmarkCell, grid: &[ForcerCase]) -> ExactCell {
    let tops = cell
        .tops
        .iter()
        .map(|top| {
            let name = format!("bench_{}_{}", cell.name, top.family);
            let case = grid
                .iter()
                .find(|case| case.name == name)
                .unwrap_or_else(|| panic!("no grid case {name}"));
            assert_eq!(
                (&case.powers, &case.numerator),
                (&top.powers, &cell.numerator),
                "{name}"
            );
            case
        })
        .collect::<Vec<_>>();
    let terms = battery_terms(cell.loops);
    let input = parse(&cell.input);
    let agree = agreement(&input, terms, |mode| frozen_sum(&tops, mode, terms))
        .inspect_err(|difference| eprintln!("benchmark cell {}: {difference}", cell.name))
        .is_ok();
    ExactCell {
        loops: cell.loops,
        grid: Some((cell.dots, cell.rank)),
        agree,
        terms,
    }
}

/// The benchmark cells of `loops` loops, in table order.
fn benchmark_cells(loops: usize) -> Vec<ExactCell> {
    let grid = grid_fixture(loops).cases;
    let cells = benchmark_inputs()
        .cells
        .into_iter()
        .filter(|cell| cell.loops == loops)
        .collect::<Vec<_>>();
    assert_eq!(cells.len(), 15, "{loops}-loop benchmark cells");
    cells.iter().map(|cell| benchmark_cell(cell, &grid)).collect()
}

fn assert_all_agree(cells: &[ExactCell]) {
    let failed = cells.iter().filter(|cell| !cell.agree).collect::<Vec<_>>();
    assert!(failed.is_empty(), "disagreeing cells (details above): {failed:?}");
}

#[test]
fn the_one_loop_benchmark_cells_agree_with_forcer() {
    assert_all_agree(&benchmark_cells(1));
}

#[test]
#[ignore = "M2-V5 STAGE 2: needs the built-in 1L-3L p library (M2-V4)"]
fn the_two_loop_benchmark_cells_agree_with_forcer() {
    assert_all_agree(&benchmark_cells(2));
}

#[test]
#[ignore = "M2-V5 STAGE 2: needs the built-in 1L-3L p library (M2-V4)"]
fn the_three_loop_benchmark_cells_agree_with_forcer() {
    assert_all_agree(&benchmark_cells(3));
}

/// The physics inputs `tests/fixtures/benchmark_p_physics.json`
/// (schema `vakint.p-physics-inputs.v1`) as `(loops, input)`.
fn physics_inputs() -> Vec<(usize, String)> {
    let text = std::fs::read_to_string(
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/benchmark_p_physics.json"),
    )
    .unwrap();
    let document: Value = serde_json::from_str(&text).unwrap();
    assert_eq!(document["schema"], "vakint.p-physics-inputs.v1");
    document["cells"]
        .as_array()
        .unwrap()
        .iter()
        .map(|cell| {
            (
                cell["loops"].as_u64().unwrap() as usize,
                cell["input"].as_str().unwrap().to_owned(),
            )
        })
        .collect()
}

/// The agreement of the physics cell of `loops` loops with a live Forcer at
/// `form`.
fn physics_cell(form: &str, loops: usize) -> ExactCell {
    let input = physics_inputs()
        .into_iter()
        .find(|(l, _)| *l == loops)
        .unwrap_or_else(|| panic!("no {loops}-loop physics input"))
        .1;
    let terms = battery_terms(loops);
    let input = parse(&input);
    let agree = agreement(&input, terms, |mode| live_forcer_value(form, &input, mode, terms))
        .inspect_err(|difference| eprintln!("physics cell L{loops}: {difference}"))
        .is_ok();
    ExactCell {
        loops,
        grid: None,
        agree,
        terms,
    }
}

/// The 2L and 3L physics cells (the 4L one is Forcer-only in benchmark v1).
const PHYSICS_LOOPS: [usize; 2] = [2, 3];

#[test]
#[ignore = "M2-V5 STAGE 2: needs the built-in 1L-3L p library (M2-V4)"]
fn the_physics_cells_agree_with_forcer() {
    let Some(form) = forcer_form() else {
        eprintln!("FORM 5 is unavailable: skipped the physics cells");
        return;
    };
    for loops in PHYSICS_LOOPS {
        assert_all_agree(&[physics_cell(&form, loops)]);
    }
}

/// The sources an exact-agreement file marks, as the benchmark harness
/// determines them: `gammaloop_sha` is the last commit outside
/// `crates/vakint/docs/benchmark`, `rustred_sha` the rustred commit pinned in
/// `Cargo.lock`.
struct Sources {
    gammaloop_sha: String,
    rustred_sha: String,
}

fn git(root: &Path, arguments: &[&str]) -> String {
    let output = Command::new("git")
        .arg("-C")
        .arg(root)
        .args(arguments)
        .output()
        .unwrap();
    assert!(output.status.success(), "git {arguments:?} failed");
    String::from_utf8(output.stdout).unwrap().trim().to_owned()
}

/// The sources of this checkout; a tree with changes outside
/// `crates/vakint/docs/benchmark` is refused (its result belongs to no
/// commit).
fn sources() -> Sources {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let paths = [".", ":(exclude)crates/vakint/docs/benchmark"];
    let status = git(&root, &[&["status", "--porcelain", "--"], &paths[..]].concat());
    assert!(
        status.is_empty(),
        "the exact-agreement file marks a commit: commit or stash first\n{status}"
    );
    let gammaloop_sha = git(&root, &[&["log", "-1", "--format=%H", "--"], &paths[..]].concat());
    let lock = std::fs::read_to_string(root.join("Cargo.lock")).unwrap();
    let rustred_sha = lock
        .lines()
        .find_map(|line| {
            let source = line.strip_prefix("source = \"git+")?;
            source.contains("/rustred?").then_some(source)?;
            Some(source.rsplit_once('#')?.1.trim_end_matches('"').to_owned())
        })
        .expect("Cargo.lock pins no rustred commit");
    for sha in [&gammaloop_sha, &rustred_sha] {
        assert!(
            sha.len() == 40 && sha.chars().all(|c| c.is_ascii_hexdigit()),
            "{sha}"
        );
    }
    Sources {
        gammaloop_sha,
        rustred_sha,
    }
}

/// The exact-agreement document (schema `vakint.p-exact.v1`, no `stub` key).
fn exact_document(sources: &Sources, cells: &[ExactCell]) -> Value {
    json!({
        "schema": "vakint.p-exact.v1",
        "gammaloop_sha": sources.gammaloop_sha,
        "rustred_sha": sources.rustred_sha,
        "cells": cells.iter().map(ExactCell::to_json).collect::<Vec<_>>(),
    })
}

/// Write `document` to `path` (pretty JSON, final newline) through a
/// temporary file next to it.
fn write_document(path: &Path, document: &Value) {
    let temporary = path.with_extension("json.partial");
    std::fs::write(
        &temporary,
        serde_json::to_string_pretty(document).unwrap() + "\n",
    )
    .unwrap();
    std::fs::rename(&temporary, path).unwrap();
}

/// Every 1L-3L benchmark cell (45) and both physics cells, computed and
/// written to `VAKINT_P_EXACT_OUT` when it is set; nothing is written unless
/// the engine covers 1L-3L, FORM 5 is there for the physics cells and the
/// tree is clean. Disagreeing cells are written as `agree: false`, and the
/// test then fails. Without the variable, only the document format is
/// checked.
#[test]
fn the_exact_agreement_file_is_written_on_request() {
    let Some(path) = std::env::var_os("VAKINT_P_EXACT_OUT").map(PathBuf::from) else {
        let fake = Sources {
            gammaloop_sha: "a".repeat(40),
            rustred_sha: "b".repeat(40),
        };
        let cells = [
            ExactCell { loops: 1, grid: Some((2, 8)), agree: true, terms: 5 },
            ExactCell { loops: 3, grid: None, agree: false, terms: 9 },
        ];
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("p_exact.json");
        write_document(&path, &exact_document(&fake, &cells));
        let written: Value = serde_json::from_str(&std::fs::read_to_string(&path).unwrap()).unwrap();
        assert_eq!(
            written,
            json!({
                "schema": "vakint.p-exact.v1",
                "gammaloop_sha": "a".repeat(40),
                "rustred_sha": "b".repeat(40),
                "cells": [
                    {"L": 1, "dots": 2, "rank": 8, "agree": true, "terms": 5},
                    {"L": 3, "physics": true, "agree": false, "terms": 9},
                ],
            })
        );
        assert!(written.get("stub").is_none());
        return;
    };
    assert_eq!(
        P_BATTERY_LOOPS,
        1..=3,
        "the exact-agreement file needs the 1L-3L library (M2-V5 STAGE 2)"
    );
    let sources = sources();
    let form = forcer_form().expect("the physics cells of the exact-agreement file need FORM 5");
    let mut cells = Vec::new();
    for loops in 1..=3 {
        cells.extend(benchmark_cells(loops));
    }
    for loops in PHYSICS_LOOPS {
        cells.push(physics_cell(&form, loops));
    }
    assert_eq!(cells.len(), 45 + 2);
    write_document(&path, &exact_document(&sources, &cells));
    eprintln!("wrote {} cells to {}", cells.len(), path.display());
    assert_all_agree(&cells);
}
