//! Tier C of the Forcer battery (M2-V5, extended to four loops by M3-V9c;
//! M2 plan, sections 8 and 9; M3 plan, gate G7): the full frozen grids and
//! the benchmark cells, RustRed against Forcer through vakint's public API
//! (`support/p_backends.rs`: `Retain` exactly in the tag basis, `Symbolic`
//! exactly, `Numerical` at 40 digits), and the exact-agreement file of the
//! benchmark.
//!
//! Full profile only (`just test_vakint_full`): `.config/nextest.toml`
//! excludes this binary from `test_vakint` by name.
//!
//! - Grids: every case of `tests/fixtures/forcer/p{1,2,3}l_grid.json` (17,
//!   39, 417) against its frozen Forcer output through `ε^(L+2)`, in parts of
//!   at most [`GRID_PART`] cases (one test each).
//! - Benchmark cells: every 1L-3L cell of
//!   `tests/fixtures/benchmark_p_integrals.json`, its frozen input (what the
//!   benchmark times) evaluated by RustRed against the sum of the frozen
//!   Forcer outputs of its tops (`bench_<cell>_<family>` in the grids).
//! - Four loops (M3-V9c), on M3's scope 1-4 package
//!   (`VAKINT_TEST_P_LIBRARY_1_4`, injected; each four-loop test reports a
//!   skip on stderr without it), on [`four_loop_cores`] workers under the
//!   default Ordered arithmetic: every case of `p4l_grid.json` (220) against
//!   its frozen Forcer output through its frozen depth (the deepest at which
//!   Forcer's value finalizes), at most `ε^6` ([`case_terms`]), in 20 parts
//!   ([`GRID4_PART`]: the 11 tops of each benchmark cell, then the sampled
//!   and dotted sectors; seconds and peak RSS per case on stderr; optional
//!   per-case limits for long measured runs, [`CASE_SECONDS_VARIABLE`],
//!   [`CASE_RSS_GB_VARIABLE`], [`SKIP_CASES_VARIABLE`]); the 15 four-loop
//!   benchmark cells are the sums of their grid tops (canonical forms, FORM-
//!   and library-free); RustRed's values are invariant under
//!   `q -> -q` and a unimodular change of loop momenta, and through the DOT
//!   round trip, for one case per family top.
//! - Physics cells: the 2L, 3L and 4L inputs of
//!   `tests/fixtures/benchmark_p_physics.json` against a live Forcer (FORM 5;
//!   skipped without it unless `VAKINT_REQUIRE_FORM5=1`); the 4L one on the
//!   1-4 package, through the deepest depth (at most `ε^6`) at which Forcer's
//!   value finalizes.
//! - With `VAKINT_P_EXACT_OUT=<path>` (never set by default),
//!   [`the_exact_agreement_file_is_written_on_request`] computes every 1L-3L
//!   cell, the 2L and 3L physics cells (default arithmetic), and the 4L cells
//!   and the 4L physics cell under the `exact` arithmetic within the caps of
//!   M3 G7 ([`EXACT4_CELL_CAP`] per cell, [`EXACT4_TOTAL_CAP`] in total; a
//!   cell beyond them is written as not exact-checked, with its cost), and
//!   writes `<path>` (schema `vakint.p-exact.v1`, read by
//!   `examples/benchmark_p_integrals` and `docs/benchmark/benchmark.typ`). A
//!   relative `<path>` is relative to the repository root (test processes
//!   run in `crates/vakint`), and its directory must exist. It refuses a
//!   dirty tree, needs a FORM 5 that runs Forcer and the 1-4 package, and
//!   writes nothing unless every cell was computed or skipped by a cap (an
//!   error aborts; a disagreement is written as `agree: false`, and the test
//!   then fails). It runs at the source of the measurement, from the
//!   repository root, e.g.
//!   `VAKINT_P_EXACT_OUT=crates/vakint/docs/benchmark/data/p_exact.json
//!   cargo test --profile ci-optim -p vakint --features vakint/feynkit-ingress
//!   --test p_integral_rustred_full -- --exact the_exact_agreement_file_is_written_on_request`
//!   (one process, no per-test time limit).
//!
//! The 1L-3L RustRed side is `Vakint::new()` with its built-in 1L-3L p
//! library (M2-V4). The four-loop tests take far longer than the full
//! profile's per-test limit at the top ranks: run them as their own nextest
//! run without that limit (`-E 'binary(p_integral_rustred_full) and
//! test(/four_loop/)'`).

#![cfg(feature = "feynkit-ingress")]

#[path = "support/p_backends.rs"]
mod p_backends;
#[path = "support/p_fixtures.rs"]
mod p_fixtures;
#[path = "support/p_fixtures4.rs"]
mod p_fixtures4;
#[path = "support/p_integrals.rs"]
mod p_integrals;
#[path = "support/p_library4.rs"]
mod p_library4;
#[path = "support/p_rustred.rs"]
mod p_rustred;
#[path = "support/p_tierc4.rs"]
mod p_tierc4;

use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{Duration, Instant};

use p_backends::*;
use p_fixtures::{benchmark_inputs, grid_fixture};
use p_fixtures4::grid4_fixture;
use p_integrals::{canonical, forcer_form, parse};
use p_library4::{PLibraryRoute, require_p_library_1_4};
use p_tierc4::*;
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use symbolica::atom::{Atom, AtomCore};
use vakint::p_integral::forcer_fixtures::{BenchmarkCell, ForcerCase};
use vakint::utils::simplify_real;
use vakint::{
    DEFAULT_ORDERED_WEIGHT_ARITHMETIC, LoopNormalizationFactor, RustRedMasterMode, Vakint,
};

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
        agreement(&input, terms, |mode| frozen_value(case, mode, terms)).unwrap_or_else(
            |difference| panic!("p{loops}l_grid.json: {}: {difference}", case.name),
        );
    }
}

#[test]
fn the_one_loop_grid_agrees_with_forcer() {
    assert_grid_part(1, 0, 1);
}

#[test]
fn the_two_loop_grid_agrees_with_forcer() {
    assert_grid_part(2, 0, 1);
}

/// The 417 three-loop cases in 11 parts.
macro_rules! three_loop_grid_parts {
    ($($name:ident: $part:expr),* $(,)?) => {$(
        #[test]
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

/// Grid cases per four-loop part: the 11 tops of one benchmark cell, so that
/// a part's engine replays the steps its tops share.
const GRID4_PART: usize = 11;

/// The two benchmark tops whose Forcer value is an exact 0 (`p4l_grid.json`,
/// M3-V7: the J11 top with `dot(k(1),p(1))*dot(k(2),p(1))`, zero and one
/// dot); RustRed's value is an exact 0 too.
const ZERO_GRID4_CASES: [&str; 2] = ["bench_L4_d0_r2_P4L_J11", "bench_L4_d1_r2_P4L_J11"];

/// The Laurent depth a frozen four-loop case is compared through: its frozen
/// depth (the deepest at which Forcer's value finalizes), at most the depth
/// of the battery, [`battery_terms`] (`ε^6`, M3 G7's "at least `ε^(L+2)`").
/// Deeper, RustRed's `Symbolic` series grows about fivefold per order (M3-V9c
/// notes: a 4L sector at 12 terms 20 s, at 16 terms over 11 min); Tier A
/// owns the deepest orders.
fn case_terms(case: &ForcerCase) -> i64 {
    let frozen = (case.depth.as_ref())
        .unwrap_or_else(|| panic!("{}: no frozen depth", case.name))
        .terms;
    frozen.min(battery_terms(4))
}

/// Optional limits of one four-loop grid case, for long measured runs
/// (unset by default: no limit): RustRed's wall seconds over its three modes
/// (`rustred_resources.max_elapsed`) and this process's RSS in GB
/// (`max_process_rss_bytes`, checked by RustRed between its steps).
const CASE_SECONDS_VARIABLE: &str = "VAKINT_TEST_P4L_CASE_SECONDS";
const CASE_RSS_GB_VARIABLE: &str = "VAKINT_TEST_P4L_CASE_RSS_GB";
/// Comma-separated four-loop grid cases not to run (e.g. ones whose cost
/// another case over the limits bounds from below); unset by default.
const SKIP_CASES_VARIABLE: &str = "VAKINT_TEST_P4L_SKIP_CASES";

fn case_limit(variable: &str) -> Option<f64> {
    let text = std::env::var(variable).ok()?;
    let limit: f64 = (text.parse().ok())
        .filter(|limit: &f64| *limit > 0.0)
        .unwrap_or_else(|| panic!("{variable}={text}: not a positive number"));
    Some(limit)
}

/// Part `part` of the four-loop grid on a fresh engine on the 1-4 package:
/// every case agrees with its frozen Forcer output through [`case_terms`]
/// under the default arithmetic (an exact 0 for [`ZERO_GRID4_CASES`]); a
/// benchmark part holds exactly the tops of the cell `cell`. The parts cover
/// the grid exactly once. Seconds and the process's peak RSS per case on
/// stderr. A case over the limits of [`CASE_SECONDS_VARIABLE`] and
/// [`CASE_RSS_GB_VARIABLE`], or listed in [`SKIP_CASES_VARIABLE`], is not
/// checked: reported on stderr, it fails the part once the other cases ran.
fn assert_grid4_part(part: usize, cell: Option<&str>) {
    let what = format!("part {} of the four-loop grid", part + 1);
    let Some(engine) = p_library4_engine(PLibraryRoute::Injected, &what) else {
        return;
    };
    let fixture = grid4_fixture();
    assert_eq!(fixture.cases.len(), 220, "p4l_grid.json");
    assert_eq!(fixture.cases.len().div_ceil(GRID4_PART), 20, "parts");
    let cases = &fixture.cases[part * GRID4_PART..(part + 1) * GRID4_PART];
    match cell {
        Some(cell) => assert!(
            (cases.iter()).all(|case| case.name.starts_with(&format!("bench_{cell}_"))),
            "{what} is not cell {cell}"
        ),
        None => assert!(cases.iter().all(|case| !case.name.starts_with("bench_"))),
    }
    let cores = four_loop_cores();
    let seconds_limit = case_limit(CASE_SECONDS_VARIABLE).map(Duration::from_secs_f64);
    let rss_limit = case_limit(CASE_RSS_GB_VARIABLE).map(|gb| (gb * (1u64 << 30) as f64) as u64);
    let limit_text = |variable: &str| std::env::var(variable).unwrap_or_else(|_| "none".into());
    let limits = format!(
        "{} s, {} GB",
        limit_text(CASE_SECONDS_VARIABLE),
        limit_text(CASE_RSS_GB_VARIABLE)
    );
    let skipped = std::env::var(SKIP_CASES_VARIABLE).unwrap_or_default();
    let skipped: Vec<&str> = skipped.split(',').map(str::trim).collect();
    let mut unchecked = Vec::new();
    for case in cases {
        if skipped.contains(&case.name.as_str()) {
            eprintln!(
                "[tier C 4L] {}: NOT CHECKED (listed in {SKIP_CASES_VARIABLE})",
                case.name
            );
            unchecked.push(format!("{} (listed)", case.name));
            continue;
        }
        let terms = case_terms(case);
        let input = case.expression().unwrap();
        let rustred = |mode| {
            let mut settings = settings4(mode, terms, None, cores);
            settings.rustred_resources.max_process_rss_bytes = rss_limit;
            settings
        };
        reset_peak_rss();
        let started = Instant::now();
        if ZERO_GRID4_CASES.contains(&case.name.as_str()) {
            for mode in MODES {
                assert!(
                    frozen_value(case, mode, terms).unwrap().is_zero(),
                    "{}",
                    case.name
                );
            }
            assert_rustred_zero(&engine, &case.name, &input, rustred);
        } else {
            let agreement = match rustred_values(&engine, &input, rustred, seconds_limit) {
                (Ok(values), _) => {
                    compare_with_forcer(&values, |mode| frozen_value(case, mode, terms))
                }
                (Err(RustRedSide::Failed(error)), _) => Err(Disagreement::Failed(error)),
                (Err(RustRedSide::OverBudget), seconds) => {
                    eprintln!(
                        "[tier C 4L] {} ({terms} terms, {cores} cores): NOT CHECKED, over the \
                         case limits ({limits}) after {seconds:.1} s, peak RSS {}",
                        case.name,
                        peak_rss_text()
                    );
                    unchecked.push(format!("{} (over the limits, {seconds:.1} s)", case.name));
                    continue;
                }
            };
            agreement.unwrap_or_else(|difference| {
                panic!("p4l_grid.json: {} ({terms} terms): {difference}", case.name)
            });
        }
        eprintln!(
            "[tier C 4L] {} ({terms} terms, {cores} cores): {:.1} s, peak RSS {}",
            case.name,
            started.elapsed().as_secs_f64(),
            peak_rss_text()
        );
    }
    assert!(
        unchecked.is_empty(),
        "{what}: {} of {} cases not checked: {}",
        unchecked.len(),
        cases.len(),
        unchecked.join(", ")
    );
}

/// The four-loop grid in 20 parts: the tops of each benchmark cell (table
/// order), then the 55 sampled and dotted sectors in 5 parts.
macro_rules! four_loop_grid_parts {
    ($($name:ident: $part:expr, $cell:expr),* $(,)?) => {$(
        #[test]
        fn $name() {
            assert_grid4_part($part, $cell);
        }
    )*};
}

four_loop_grid_parts!(
    the_four_loop_grid_l4_d0_r0_tops_agree_with_forcer: 0, Some("L4_d0_r0"),
    the_four_loop_grid_l4_d0_r2_tops_agree_with_forcer: 1, Some("L4_d0_r2"),
    the_four_loop_grid_l4_d0_r4_tops_agree_with_forcer: 2, Some("L4_d0_r4"),
    the_four_loop_grid_l4_d0_r6_tops_agree_with_forcer: 3, Some("L4_d0_r6"),
    the_four_loop_grid_l4_d0_r8_tops_agree_with_forcer: 4, Some("L4_d0_r8"),
    the_four_loop_grid_l4_d1_r0_tops_agree_with_forcer: 5, Some("L4_d1_r0"),
    the_four_loop_grid_l4_d1_r2_tops_agree_with_forcer: 6, Some("L4_d1_r2"),
    the_four_loop_grid_l4_d1_r4_tops_agree_with_forcer: 7, Some("L4_d1_r4"),
    the_four_loop_grid_l4_d1_r6_tops_agree_with_forcer: 8, Some("L4_d1_r6"),
    the_four_loop_grid_l4_d1_r8_tops_agree_with_forcer: 9, Some("L4_d1_r8"),
    the_four_loop_grid_l4_d2_r0_tops_agree_with_forcer: 10, Some("L4_d2_r0"),
    the_four_loop_grid_l4_d2_r2_tops_agree_with_forcer: 11, Some("L4_d2_r2"),
    the_four_loop_grid_l4_d2_r4_tops_agree_with_forcer: 12, Some("L4_d2_r4"),
    the_four_loop_grid_l4_d2_r6_tops_agree_with_forcer: 13, Some("L4_d2_r6"),
    the_four_loop_grid_l4_d2_r8_tops_agree_with_forcer: 14, Some("L4_d2_r8"),
    the_four_loop_grid_sectors_part_1_agree_with_forcer: 15, None,
    the_four_loop_grid_sectors_part_2_agree_with_forcer: 16, None,
    the_four_loop_grid_sectors_part_3_agree_with_forcer: 17, None,
    the_four_loop_grid_sectors_part_4_agree_with_forcer: 18, None,
    the_four_loop_grid_sectors_part_5_agree_with_forcer: 19, None,
);

/// RustRed's values on the 1-4 package are invariant under `q -> -q` and a
/// unimodular change of loop momenta (with line reversals), one case per
/// four-loop family top: the top times `1 + dot(k(1),p(1))` (the odd part
/// alone is 0 on the LA4, J4 and J11 tops, for Forcer too).
#[test]
fn four_loop_values_are_routing_invariant() {
    let what = "the four-loop routing invariance";
    let Some(engine) = p_library4_engine(PLibraryRoute::Injected, what) else {
        return;
    };
    for family in FOUR_LOOP_FAMILIES {
        let input = parse("1+dot(k(1),p(1))") * long_form_top(family, &[1; 11]);
        assert_routing_invariant(&engine, &input, 0..1);
    }
}

/// The DOT round trip at four loops, one case per family top (a dot on line
/// 1): the DOT export read back is the same integral with RustRed's same
/// values on the 1-4 package.
#[test]
fn four_loop_dot_graphs_keep_the_integral_and_its_value() {
    let what = "the four-loop DOT round trip";
    let Some(engine) = p_library4_engine(PLibraryRoute::Injected, what) else {
        return;
    };
    let mut powers = [1; 11];
    powers[0] = 2;
    for family in FOUR_LOOP_FAMILIES {
        assert_dot_round_trip_value(&engine, &long_form_top(family, &powers));
    }
}

/// One cell of the exact-agreement file: a benchmark cell or a physics cell,
/// with the Ordered arithmetic of its RustRed side and, at four loops,
/// RustRed's seconds (the cost the caps of M3 G7 bound).
#[derive(Clone, Debug, PartialEq)]
struct ExactCell {
    loops: usize,
    /// `(dots, rank)` of a benchmark cell, `None` for the physics cell.
    grid: Option<(usize, usize)>,
    agree: bool,
    terms: i64,
    arithmetic: &'static str,
    rustred_seconds: Option<f64>,
}

/// The keys `L` and `dots`, `rank` (benchmark cell) or `physics` of a cell.
fn cell_key(loops: usize, grid: Option<(usize, usize)>) -> serde_json::Map<String, Value> {
    let key = match grid {
        Some((dots, rank)) => json!({"L": loops, "dots": dots, "rank": rank}),
        None => json!({"L": loops, "physics": true}),
    };
    key.as_object().unwrap().clone()
}

impl ExactCell {
    fn to_json(&self) -> Value {
        let mut row = cell_key(self.loops, self.grid);
        row.insert("agree".into(), json!(self.agree));
        row.insert("terms".into(), json!(self.terms));
        row.insert("arithmetic".into(), json!(self.arithmetic));
        if let Some(seconds) = self.rustred_seconds {
            row.insert(
                "rustred_seconds".into(),
                json!((seconds * 10.0).round() / 10.0),
            );
        }
        Value::Object(row)
    }
}

/// A four-loop cell that is not exact-checked: its exact run exceeded, or
/// would exceed, a cap of M3 G7; `cost` says by how much (`> 600 s`, ...).
#[derive(Clone, Debug, PartialEq)]
struct NotExactChecked {
    loops: usize,
    grid: Option<(usize, usize)>,
    cost: String,
}

impl NotExactChecked {
    fn to_json(&self) -> Value {
        let mut row = cell_key(self.loops, self.grid);
        let note = format!("not exact-checked (cost {})", self.cost);
        row.insert("note".into(), json!(note));
        Value::Object(row)
    }
}

/// Whether a cell agrees: a disagreement of computed values is reported on
/// stderr and recorded as `false`; a cell that was not computed (an error of
/// either side, or a Forcer value 0) aborts, so that no file marks it.
fn cell_agrees(cell: &str, agreement: Result<(), Disagreement>) -> bool {
    match agreement {
        Ok(()) => true,
        Err(Disagreement::Differs(difference)) => {
            eprintln!("{cell}: {difference}");
            false
        }
        Err(Disagreement::Failed(error)) => panic!("{cell} was not computed: {error}"),
    }
}

/// The grid cases of the tops of a benchmark cell (`bench_<cell>_<family>`),
/// each with the cell's numerator and its top's powers.
fn cell_tops<'a>(cell: &BenchmarkCell, grid: &'a [ForcerCase]) -> Vec<&'a ForcerCase> {
    (cell.tops.iter())
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
        .collect()
}

/// The agreement of a benchmark cell: its frozen input with RustRed against
/// the frozen Forcer outputs of its tops ([`cell_agrees`]).
fn benchmark_cell(cell: &BenchmarkCell, grid: &[ForcerCase]) -> ExactCell {
    let tops = cell_tops(cell, grid);
    let terms = battery_terms(cell.loops);
    let input = parse(&cell.input);
    let agree = cell_agrees(
        &format!("benchmark cell {}", cell.name),
        agreement(&input, terms, |mode| frozen_sum(&tops, mode, terms)),
    );
    ExactCell {
        loops: cell.loops,
        grid: Some((cell.dots, cell.rank)),
        agree,
        terms,
        arithmetic: DEFAULT_ORDERED_WEIGHT_ARITHMETIC,
        rustred_seconds: None,
    }
}

/// The benchmark cells of `loops` loops (15), in table order.
fn benchmark_cells_of(loops: usize) -> Vec<BenchmarkCell> {
    let cells = (benchmark_inputs().cells.into_iter())
        .filter(|cell| cell.loops == loops)
        .collect::<Vec<_>>();
    assert_eq!(cells.len(), 15, "{loops}-loop benchmark cells");
    cells
}

/// The agreements of the benchmark cells of `loops` (1 to 3) loops, in table
/// order.
fn benchmark_cells(loops: usize) -> Vec<ExactCell> {
    let grid = grid_fixture(loops).cases;
    (benchmark_cells_of(loops).iter())
        .map(|cell| benchmark_cell(cell, &grid))
        .collect()
}

fn assert_all_agree(cells: &[ExactCell]) {
    let failed = cells.iter().filter(|cell| !cell.agree).collect::<Vec<_>>();
    assert!(
        failed.is_empty(),
        "disagreeing cells (details above): {failed:?}"
    );
}

#[test]
fn the_one_loop_benchmark_cells_agree_with_forcer() {
    assert_all_agree(&benchmark_cells(1));
}

#[test]
fn the_two_loop_benchmark_cells_agree_with_forcer() {
    assert_all_agree(&benchmark_cells(2));
}

#[test]
fn the_three_loop_benchmark_cells_agree_with_forcer() {
    assert_all_agree(&benchmark_cells(3));
}

/// The four-loop benchmark cells (their frozen inputs, what the benchmark
/// times) are the sums of their tops in `p4l_grid.json`: equal canonical
/// forms, so the four-loop grid parts check every cell under the default
/// arithmetic. FORM- and library-free.
#[test]
fn the_four_loop_benchmark_cells_are_the_sums_of_their_grid_tops() {
    let grid = grid4_fixture().cases;
    for cell in benchmark_cells_of(4) {
        let tops = cell_tops(&cell, &grid);
        let sum = (tops.iter()).fold(Atom::Zero, |sum, top| sum + top.expression().unwrap());
        let difference = (canonical(&parse(&cell.input)) - canonical(&sum)).expand();
        assert!(difference.is_zero(), "{}: {difference}", cell.name);
    }
}

/// The input of the physics cell of `loops` loops
/// (`tests/fixtures/benchmark_p_physics.json`, schema
/// `vakint.p-physics-inputs.v1`).
fn physics_input(loops: usize) -> Atom {
    let text = std::fs::read_to_string(
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/benchmark_p_physics.json"),
    )
    .unwrap();
    let document: Value = serde_json::from_str(&text).unwrap();
    assert_eq!(document["schema"], "vakint.p-physics-inputs.v1");
    let cell = (document["cells"].as_array().unwrap().iter())
        .find(|cell| cell["loops"].as_u64() == Some(loops as u64))
        .unwrap_or_else(|| panic!("no {loops}-loop physics input"));
    parse(cell["input"].as_str().unwrap())
}

/// The agreement of the physics cell of `loops` (2 or 3) loops with a live
/// Forcer at `form` ([`cell_agrees`]).
fn physics_cell(form: &str, loops: usize) -> ExactCell {
    let terms = battery_terms(loops);
    let input = physics_input(loops);
    let agree = cell_agrees(
        &format!("physics cell L{loops}"),
        agreement(&input, terms, |mode| {
            live_forcer_value(form, &input, mode, terms)
        }),
    );
    ExactCell {
        loops,
        grid: None,
        agree,
        terms,
        arithmetic: DEFAULT_ORDERED_WEIGHT_ARITHMETIC,
        rustred_seconds: None,
    }
}

/// The physics cells: 2L and 3L on the built-in library, 4L on the 1-4
/// package (M3 G7; no frozen Forcer value, so a live Forcer).
const PHYSICS_LOOPS: [usize; 3] = [2, 3, 4];

#[test]
fn the_physics_cells_agree_with_forcer() {
    let Some(form) = live_form("the physics cells") else {
        return;
    };
    for loops in PHYSICS_LOOPS
        .into_iter()
        .filter(|loops| P_BATTERY_LOOPS.contains(loops))
    {
        assert_all_agree(&[physics_cell(&form, loops)]);
    }
}

/// The 4L physics cell on the 1-4 package under the default arithmetic
/// against a live Forcer, through the deepest depth (at most `ε^6`) at which
/// Forcer's value finalizes ([`live_terms`]).
#[test]
fn the_four_loop_physics_cell_agrees_with_forcer() {
    let what = "the four-loop physics cell";
    let Some(engine) = p_library4_engine(PLibraryRoute::Injected, what) else {
        return;
    };
    let Some(form) = live_form(what) else {
        return;
    };
    let input = physics_input(4);
    let terms = live_terms(&form, &input, 4);
    let cores = four_loop_cores();
    let seconds = agreement4(
        &engine,
        &input,
        |mode| settings4(mode, terms, None, cores),
        |mode| live_forcer_value(&form, &input, mode, terms),
    )
    .unwrap_or_else(|difference| panic!("physics cell L4 ({terms} terms): {difference}"));
    eprintln!("[tier C 4L] physics_L4 ({terms} terms, {cores} cores): RustRed {seconds:.1} s");
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

/// The root of the gammaloop repository (test processes run in
/// `crates/vakint`).
fn repository_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

/// `VAKINT_P_EXACT_OUT`, a relative path taken relative to
/// [`repository_root`]; its directory must exist (checked before any cell
/// is computed).
fn exact_output(path: PathBuf) -> PathBuf {
    let path = if path.is_relative() {
        repository_root().join(path)
    } else {
        path
    };
    let directory = path.parent().filter(|parent| parent.is_dir());
    assert!(
        directory.is_some() && path.file_name().is_some(),
        "VAKINT_P_EXACT_OUT={}: no such directory",
        path.display()
    );
    path
}

/// A relative `VAKINT_P_EXACT_OUT` whose directory does not exist below the
/// repository root is refused before any computation.
#[test]
#[should_panic(expected = "no such directory")]
fn an_exact_output_without_its_directory_is_refused() {
    exact_output(PathBuf::from("crates/vakint/crates/vakint/p_exact.json"));
}

/// The sources of this checkout; a tree with changes outside
/// `crates/vakint/docs/benchmark` is refused (its result belongs to no
/// commit).
fn sources() -> Sources {
    let root = repository_root();
    let paths = [".", ":(exclude)crates/vakint/docs/benchmark"];
    let status = git(
        &root,
        &[&["status", "--porcelain", "--"], &paths[..]].concat(),
    );
    assert!(
        status.is_empty(),
        "the exact-agreement file marks a commit: commit or stash first\n{status}"
    );
    let gammaloop_sha = git(
        &root,
        &[&["log", "-1", "--format=%H", "--"], &paths[..]].concat(),
    );
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

/// The four-loop part of an exact-agreement file: the 1-4 package it was
/// computed on (`library.json`: its identity and sha256), the cells that
/// are not exact-checked, and the exact seconds spent.
struct FourLoopExact {
    identity: String,
    library_sha256: String,
    not_exact_checked: Vec<NotExactChecked>,
    rustred_seconds: f64,
}

/// The exact-agreement document (schema `vakint.p-exact.v1`, no `stub` key).
fn exact_document(sources: &Sources, cells: &[ExactCell], four: &FourLoopExact) -> Value {
    json!({
        "schema": "vakint.p-exact.v1",
        "gammaloop_sha": sources.gammaloop_sha,
        "rustred_sha": sources.rustred_sha,
        "p_library_1_4": {"identity": four.identity, "library_sha256": four.library_sha256},
        "exact_4l": {
            "arithmetic": EXACT4_ARITHMETIC,
            "cell_cap_seconds": EXACT4_CELL_CAP.as_secs(),
            "total_cap_seconds": EXACT4_TOTAL_CAP.as_secs(),
            "rustred_seconds": four.rustred_seconds.round(),
        },
        "cells": cells.iter().map(ExactCell::to_json).collect::<Vec<_>>(),
        "not_exact_checked": (four.not_exact_checked.iter())
            .map(NotExactChecked::to_json)
            .collect::<Vec<_>>(),
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

/// The cap of one four-loop cell's exact run (M3 G7): its three RustRed
/// evaluations under [`EXACT4_ARITHMETIC`].
const EXACT4_CELL_CAP: Duration = Duration::from_secs(600);
/// The cap of every four-loop exact run together (M3 G7).
const EXACT4_TOTAL_CAP: Duration = Duration::from_secs(7200);
/// The Ordered arithmetic of the four-loop cells of the exact-agreement file.
const EXACT4_ARITHMETIC: &str = "exact";

/// The Forcer side of a cell, per mode.
type Forcer<'a> = dyn Fn(RustRedMasterMode) -> Result<Atom, String> + 'a;

/// The four-loop cells of the exact-agreement file on `engine` (the 1-4
/// package), under [`EXACT4_ARITHMETIC`] on [`four_loop_cores`] workers,
/// through their deepest common depth (the least [`case_terms`] of the tops;
/// a live probe for the physics cell): the 15 benchmark cells in increasing
/// dots and rank, then the physics cell. Each RustRed side runs within
/// [`EXACT4_CELL_CAP`] and what is left of [`EXACT4_TOTAL_CAP`]; a cell that
/// runs out is not exact-checked, and neither is a later benchmark cell with
/// at least its dots and rank (its run costs more). Returns the checked cells
/// and the others with their cost, and the exact seconds spent.
fn four_loop_exact_cells(
    engine: &Vakint,
    form: &str,
) -> (Vec<ExactCell>, Vec<NotExactChecked>, f64) {
    let grid = grid4_fixture().cases;
    let mut cells = benchmark_cells_of(4);
    cells.sort_by_key(|cell| (cell.dots, cell.rank));
    let cores = four_loop_cores();
    let (mut checked, mut skipped, mut spent) = (Vec::new(), Vec::new(), 0.0);
    // The cells that ran out of a whole EXACT4_CELL_CAP: `(dots, rank, name)`.
    let mut over = Vec::<(usize, usize, String)>::new();
    let mut run = |grid, name: &str, input: &Atom, terms, forcer: &Forcer| {
        let left = EXACT4_TOTAL_CAP.saturating_sub(Duration::from_secs_f64(spent));
        let budget = EXACT4_CELL_CAP.min(left);
        if budget.is_zero() {
            let cost = format!(
                "unknown: the {} s total cap is spent",
                EXACT4_TOTAL_CAP.as_secs()
            );
            return Err((cost, false));
        }
        let exact = |mode| settings4(mode, terms, Some(EXACT4_ARITHMETIC), cores);
        let outcome = match capped_agreement(engine, input, exact, forcer, budget) {
            Capped::Done(agreement, seconds) => {
                spent += seconds;
                Ok(ExactCell {
                    loops: 4,
                    grid,
                    agree: cell_agrees(name, agreement),
                    terms,
                    arithmetic: EXACT4_ARITHMETIC,
                    rustred_seconds: Some(seconds),
                })
            }
            Capped::OverBudget(seconds) => {
                spent += seconds;
                let cost = format!("> {} s", budget.as_secs());
                Err((cost, budget == EXACT4_CELL_CAP))
            }
        };
        eprintln!("[exact 4L] {name} ({terms} terms): {outcome:?}; {spent:.0} s spent");
        outcome
    };
    for cell in &cells {
        let key = Some((cell.dots, cell.rank));
        let lower = (over.iter())
            .find(|(dots, rank, _)| *dots <= cell.dots && *rank <= cell.rank)
            .map(|(_, _, name)| name.clone());
        let outcome = match lower {
            Some(lower) => Err((
                format!(
                    "> {} s, like {lower}, which has no more dots and no higher rank",
                    EXACT4_CELL_CAP.as_secs()
                ),
                false,
            )),
            None => {
                let tops = cell_tops(cell, &grid);
                let terms = tops.iter().map(|top| case_terms(top)).min().unwrap();
                let forcer = |mode| frozen_sum(&tops, mode, terms);
                let name = format!("benchmark cell {}", cell.name);
                run(key, &name, &parse(&cell.input), terms, &forcer)
            }
        };
        match outcome {
            Ok(exact) => checked.push(exact),
            Err((cost, whole_cap)) => {
                if whole_cap {
                    over.push((cell.dots, cell.rank, cell.name.clone()));
                }
                skipped.push(NotExactChecked {
                    loops: 4,
                    grid: key,
                    cost,
                });
            }
        }
    }
    let input = physics_input(4);
    let terms = live_terms(form, &input, 4);
    let forcer = |mode| live_forcer_value(form, &input, mode, terms);
    match run(None, "physics cell L4", &input, terms, &forcer) {
        Ok(exact) => checked.push(exact),
        Err((cost, _)) => skipped.push(NotExactChecked {
            loops: 4,
            grid: None,
            cost,
        }),
    }
    (checked, skipped, spent)
}

/// The identity and the sha256 of the `library.json` of the package in
/// `directory`.
fn library_mark(directory: &Path) -> (String, String) {
    let bytes = std::fs::read(directory.join("library.json")).unwrap();
    let description: Value = serde_json::from_slice(&bytes).unwrap();
    let identity = description["identity"].as_str().unwrap().to_owned();
    (identity, format!("{:x}", Sha256::digest(&bytes)))
}

/// Every benchmark cell (45 at 1L-3L, 15 at 4L) and the physics cells (2L,
/// 3L, 4L), computed and written to `VAKINT_P_EXACT_OUT` ([`exact_output`])
/// when it is set: 1L-3L on the built-in library under the default
/// arithmetic, 4L on the 1-4 package under the caps of M3 G7
/// ([`four_loop_exact_cells`]). Nothing is written unless the default engine
/// covers 1L-3L, the 1-4 package is given, FORM 5 runs Forcer for the
/// physics cells, the tree is clean and every cell was computed or capped.
/// Disagreeing cells are written as `agree: false`, and the test then fails.
/// Without the variable, only the document format and the output path rule
/// are checked.
#[test]
fn the_exact_agreement_file_is_written_on_request() {
    let Some(path) = std::env::var_os("VAKINT_P_EXACT_OUT").map(PathBuf::from) else {
        assert_eq!(
            exact_output(PathBuf::from("crates/vakint/Cargo.toml")),
            PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../crates/vakint/Cargo.toml")
        );
        assert_eq!(
            exact_output(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("p_exact.json")),
            PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("p_exact.json")
        );
        let fake = Sources {
            gammaloop_sha: "a".repeat(40),
            rustred_sha: "b".repeat(40),
        };
        let cell = |loops, grid, agree, terms, arithmetic, rustred_seconds| ExactCell {
            loops,
            grid,
            agree,
            terms,
            arithmetic,
            rustred_seconds,
        };
        let cells = [
            cell(1, Some((2, 8)), true, 5, "modular_polynomial", None),
            cell(3, None, false, 9, "modular_polynomial", None),
            cell(4, Some((0, 2)), true, 10, "exact", Some(12.34)),
        ];
        let four = FourLoopExact {
            identity: "p-library:scope=1-4".into(),
            library_sha256: "c".repeat(64),
            not_exact_checked: vec![
                NotExactChecked {
                    loops: 4,
                    grid: Some((1, 4)),
                    cost: "> 600 s".into(),
                },
                NotExactChecked {
                    loops: 4,
                    grid: None,
                    cost: "> 551 s".into(),
                },
            ],
            rustred_seconds: 1163.4,
        };
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("p_exact.json");
        write_document(&path, &exact_document(&fake, &cells, &four));
        let written: Value =
            serde_json::from_str(&std::fs::read_to_string(&path).unwrap()).unwrap();
        assert_eq!(
            written,
            json!({
                "schema": "vakint.p-exact.v1",
                "gammaloop_sha": "a".repeat(40),
                "rustred_sha": "b".repeat(40),
                "p_library_1_4": {"identity": "p-library:scope=1-4", "library_sha256": "c".repeat(64)},
                "exact_4l": {
                    "arithmetic": "exact",
                    "cell_cap_seconds": 600,
                    "total_cap_seconds": 7200,
                    "rustred_seconds": 1163.0,
                },
                "cells": [
                    {"L": 1, "dots": 2, "rank": 8, "agree": true, "terms": 5,
                     "arithmetic": "modular_polynomial"},
                    {"L": 3, "physics": true, "agree": false, "terms": 9,
                     "arithmetic": "modular_polynomial"},
                    {"L": 4, "dots": 0, "rank": 2, "agree": true, "terms": 10,
                     "arithmetic": "exact", "rustred_seconds": 12.3},
                ],
                "not_exact_checked": [
                    {"L": 4, "dots": 1, "rank": 4, "note": "not exact-checked (cost > 600 s)"},
                    {"L": 4, "physics": true, "note": "not exact-checked (cost > 551 s)"},
                ],
            })
        );
        assert!(written.get("stub").is_none());
        return;
    };
    assert_eq!(
        P_BATTERY_LOOPS,
        1..=3,
        "the exact-agreement file needs the 1L-3L p library"
    );
    let library = require_p_library_1_4("the four-loop cells of the exact-agreement file");
    let path = exact_output(path);
    let sources = sources();
    let form = forcer_form().expect("the physics cells of the exact-agreement file need FORM 5");
    assert_forcer_runs(&form);
    let mut cells = Vec::new();
    for loops in 1..=3 {
        cells.extend(benchmark_cells(loops));
    }
    for loops in PHYSICS_LOOPS
        .into_iter()
        .filter(|loops| P_BATTERY_LOOPS.contains(loops))
    {
        cells.push(physics_cell(&form, loops));
    }
    let engine = p_library4_engine(PLibraryRoute::Injected, "").unwrap();
    let (four_loop, not_exact_checked, rustred_seconds) = four_loop_exact_cells(&engine, &form);
    cells.extend(four_loop);
    assert_eq!(
        cells.len() + not_exact_checked.len(),
        45 + 15 + PHYSICS_LOOPS.len()
    );
    let (identity, library_sha256) = library_mark(&library);
    let four = FourLoopExact {
        identity,
        library_sha256,
        not_exact_checked,
        rustred_seconds,
    };
    write_document(&path, &exact_document(&sources, &cells, &four));
    eprintln!(
        "wrote {} cells and {} not exact-checked to {} ({rustred_seconds:.0} s exact at 4L)",
        cells.len(),
        four.not_exact_checked.len(),
        path.display()
    );
    assert_all_agree(&cells);
}
