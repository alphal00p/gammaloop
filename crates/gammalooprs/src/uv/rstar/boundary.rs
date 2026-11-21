//! The remainder-boundary producer (`GAMMALOOP_REMAINDER_BOUNDARY=1|check`; on by
//! default, `0`/`off` restores the expanded path).
//!
//! Today one root-scope R* remainder reaches vakint as the expanded
//! `Σ_{i+j≤deg} P_i D_j`: the contracted numerator times the Taylor
//! coefficients of its denominator product, expanded into monomials and
//! grouped per topology. On the new path (beta_computation DESIGN.md v2) the
//! same remainder is handed over factored, as a vakint [`RemainderBoundary`]:
//!
//! * the denominator product alone is Taylor expanded ([`denominator_orders`]);
//!   its monomials give the power patterns, their vacuum frames
//!   (`TermFrame`) and their cut (`cut_frame`), with one composite routing;
//! * the numerator `P = C(σ(N))` is contracted without the denominators and
//!   split by one pass of a recursive grammar into blocks `S_g` of grade
//!   `g = (T, s, b)` (k·p structure, power of p², power of the reference mass
//!   squared), which are free of `p` and of the reference mass ([`Grammar`]);
//! * per pattern and grade, the small coefficient `R_{π,g}` carries the
//!   structure `t_g`, the mass compensation and the denominator orders,
//!   projected onto the vacuum tensors and multiplied by the scalar root
//!   normalisation.
//!
//! Nothing here expands a product `P_i × D_j` or a sum of loop dots. Every
//! block, coefficient and routing image is checked against the vakint
//! contract before the hand-over ([`verify_vakint_scalar`]).

use std::{
    collections::{BTreeMap, BTreeSet, HashMap},
    sync::LazyLock,
    time::Instant,
};

use eyre::{Result, WrapErr, eyre};
use linnet::half_edge::{involution::EdgeIndex, subgraph::SuBitGraph};
use spenso::{
    network::{library::symbolic::ETS, tags::SPENSO_TAG},
    structure::{
        representation::{Minkowski, RepName},
        slot::{SlotMatch, SlotMatcher},
    },
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol},
    function,
    id::Replacement,
    symbol,
};
use vakint::{
    TensorReductionMethod, Vakint, VakintSettings,
    p_integral::{BoundaryOptions, BoundaryPattern, MassiveClosure, RemainderBoundary},
    symbols::S,
};

use crate::{
    utils::GS,
    uv::approx::{
        integrated::{TermFrame, TermRewrites, term_frame},
        streamed,
    },
};

use super::{
    TaylorExpansion,
    expression::RStarIntegrand,
    integrate::{VacuumCut, cut_frame, tensor_reduce_by_factors},
};

// ---------------------------------------------------------------------------
// Switches
// ---------------------------------------------------------------------------

/// `GAMMALOOP_REMAINDER_BOUNDARY`: the expanded path (`Off`: an explicit `0`,
/// `off` or empty value), the boundary for every admitted remainder (`On`:
/// `1`, `on`, or unset, the production setting since r5), or both compared
/// (`Check`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BoundaryMode {
    Off,
    On,
    Check,
}

impl BoundaryMode {
    pub fn as_str(&self) -> &'static str {
        match self {
            BoundaryMode::Off => "off",
            BoundaryMode::On => "on",
            BoundaryMode::Check => "check",
        }
    }
}

fn parse_mode(value: Option<&str>) -> std::result::Result<BoundaryMode, String> {
    match value {
        Some("" | "0" | "off") => Ok(BoundaryMode::Off),
        None | Some("1" | "on") => Ok(BoundaryMode::On),
        Some("check") => Ok(BoundaryMode::Check),
        Some(other) => Err(format!(
            "GAMMALOOP_REMAINDER_BOUNDARY={other:?}: expected unset, 0, off, 1, on or check"
        )),
    }
}

static MODE: LazyLock<std::result::Result<BoundaryMode, String>> =
    LazyLock::new(|| match std::env::var("GAMMALOOP_REMAINDER_BOUNDARY") {
        Ok(value) => parse_mode(Some(&value)),
        Err(std::env::VarError::NotPresent) => parse_mode(None),
        Err(std::env::VarError::NotUnicode(value)) => Err(format!(
            "GAMMALOOP_REMAINDER_BOUNDARY={value:?} is not unicode"
        )),
    });

static STRICT: LazyLock<bool> =
    LazyLock::new(|| std::env::var("GAMMALOOP_REMAINDER_BOUNDARY_STRICT").as_deref() == Ok("1"));

#[cfg(test)]
thread_local! {
    static TEST_MODE: std::cell::Cell<Option<(BoundaryMode, bool)>> = const { std::cell::Cell::new(None) };
}

/// The mode of `GAMMALOOP_REMAINDER_BOUNDARY`, read once; a value other than
/// unset, empty, `0`, `off`, `1`, `on` or `check` is an error at first use.
pub fn remainder_boundary_mode() -> Result<BoundaryMode> {
    #[cfg(test)]
    if let Some((mode, _)) = TEST_MODE.with(|slot| slot.get()) {
        return Ok(mode);
    }
    MODE.clone().map_err(|error| eyre!(error))
}

/// `GAMMALOOP_REMAINDER_BOUNDARY_STRICT=1`: every decline after admission is
/// an error, and vakint pattern declines are errors too.
pub(crate) fn boundary_strict() -> bool {
    #[cfg(test)]
    if let Some((_, strict)) = TEST_MODE.with(|slot| slot.get()) {
        return strict;
    }
    *STRICT
}

/// Tests: run `f` with the mode (and strictness) set on this thread, never
/// through the environment.
#[cfg(test)]
pub(crate) fn with_boundary_mode<T>(mode: BoundaryMode, strict: bool, f: impl FnOnce() -> T) -> T {
    struct Restore(Option<(BoundaryMode, bool)>);
    impl Drop for Restore {
        fn drop(&mut self) {
            let previous = self.0;
            TEST_MODE.with(|slot| slot.set(previous));
        }
    }
    let _restore = Restore(TEST_MODE.with(|slot| slot.replace(Some((mode, strict)))));
    f()
}

/// `GAMMALOOP_BOUNDARY_FALLBACK_MAX_SUMMANDS` (default 1,000,000 summands of
/// `P`; 0 = no limit): above it a decline is an error in every mode and vakint
/// runs strict (DESIGN 10c).
pub(crate) fn fallback_max_summands() -> usize {
    static MAX: LazyLock<usize> = LazyLock::new(|| {
        match std::env::var("GAMMALOOP_BOUNDARY_FALLBACK_MAX_SUMMANDS")
            .ok()
            .and_then(|value| value.parse::<usize>().ok())
        {
            Some(0) => usize::MAX,
            Some(max) => max,
            None => 1_000_000,
        }
    });
    *MAX
}

/// The vakint environment gates that `lower_direct_blocks` does not apply
/// itself (DESIGN A10, vakint `evaluate_direct`): `VAKINT_LOWERING` unset or
/// `direct`, the fold not off, no `VAKINT_P_FOLD_CHECK`. Read once.
pub(crate) fn vakint_gates() -> std::result::Result<(), &'static str> {
    static GATES: LazyLock<std::result::Result<(), &'static str>> = LazyLock::new(|| {
        match std::env::var("VAKINT_LOWERING").as_deref() {
            Err(_) | Ok("direct") => {}
            Ok(_) => return Err("lowering"),
        }
        if matches!(
            std::env::var("VAKINT_P_FOLD").as_deref(),
            Ok("off" | "0" | "false")
        ) {
            return Err("fold_off");
        }
        if std::env::var_os("VAKINT_P_FOLD_CHECK").is_some() {
            return Err("fold_check");
        }
        Ok(())
    });
    *GATES
}

/// The settings part of the A10 gate: RustRed present and before Forcer, as
/// the direct reader requires (`evaluate_direct`).
pub(crate) fn settings_gate(settings: &VakintSettings) -> std::result::Result<(), &'static str> {
    let first = settings.evaluation_order.0.iter().find(|method| {
        matches!(
            method,
            vakint::EvaluationMethod::RustRed(_) | vakint::EvaluationMethod::Forcer(_)
        )
    });
    match first {
        Some(vakint::EvaluationMethod::RustRed(_)) => Ok(()),
        Some(_) => Err("forcer_before_rustred"),
        None => Err("no_rustred"),
    }
}

/// `GAMMALOOP_BOUNDARY_FALLBACK_MAX_BYTES` (default 256 MiB of `P`; 0 = no
/// limit). idenso returns the contracted closed scalar factorized, often as a
/// single top-level product, so the summand count of `P` is 1 and says
/// nothing about its size; its byte size does.
pub(crate) fn fallback_max_bytes() -> usize {
    static MAX: LazyLock<usize> = LazyLock::new(|| {
        match std::env::var("GAMMALOOP_BOUNDARY_FALLBACK_MAX_BYTES")
            .ok()
            .and_then(|value| value.parse::<usize>().ok())
        {
            Some(0) => usize::MAX,
            Some(max) => max,
            None => 256 << 20,
        }
    });
    *MAX
}

/// Remainders with at least this many loops are always large (DESIGN
/// Addendum C, extended after the contraction by the step 1-3 review): today's
/// path would contract or expand the whole remainder (1237: about 40 GB).
pub(crate) const LARGE_LOOPS: usize = 6;

/// Why a remainder is large (DESIGN 10c, Addendum C), or `None`: at least
/// [`LARGE_LOOPS`] loops; any remainder of a census dump; after the
/// contraction, `P` above the byte or the summand limit. One predicate for
/// both decisions: a large decline never falls back to today's path, and a
/// large remainder runs vakint strict (no `ordinary()` expansion of a pattern).
pub(super) fn large_reason(
    stats: &BoundaryStats,
    loops: usize,
    census_dump: bool,
) -> Option<String> {
    if loops >= LARGE_LOOPS {
        return Some(format!("{loops} loops (large from {LARGE_LOOPS})"));
    }
    if census_dump {
        return Some("a census dump (PINT6_CLOSURE_DUMP)".to_owned());
    }
    if stats.contracted {
        if stats.p_bytes > fallback_max_bytes() {
            return Some(format!(
                "P of {} bytes, limit {}",
                stats.p_bytes,
                fallback_max_bytes()
            ));
        }
        if stats.p_summands > fallback_max_summands() {
            return Some(format!(
                "P of {} summands, limit {}",
                stats.p_summands,
                fallback_max_summands()
            ));
        }
    }
    None
}

/// [`large_reason`] for the running process (census dump from the
/// environment).
pub(super) fn is_large(stats: &BoundaryStats, loops: usize) -> bool {
    large_reason(stats, loops, census_dump_active()).is_some()
}

#[cfg(test)]
thread_local! {
    static TEST_CENSUS_DUMP: std::cell::Cell<Option<bool>> = const { std::cell::Cell::new(None) };
}

/// Whether a census dump (`PINT6_CLOSURE_DUMP`) is running.
pub(super) fn census_dump_active() -> bool {
    #[cfg(test)]
    if let Some(active) = TEST_CENSUS_DUMP.with(|slot| slot.get()) {
        return active;
    }
    super::pint6_dump::settings().is_some()
}

/// Tests: run `f` as if a census dump were (or were not) running on this
/// thread. Only the remainder-boundary decisions read it; nothing is written.
#[cfg(test)]
pub(crate) fn with_census_dump<T>(active: bool, f: impl FnOnce() -> T) -> T {
    struct Restore(Option<bool>);
    impl Drop for Restore {
        fn drop(&mut self) {
            let previous = self.0;
            TEST_CENSUS_DUMP.with(|slot| slot.set(previous));
        }
    }
    let _restore = Restore(TEST_CENSUS_DUMP.with(|slot| slot.replace(Some(active))));
    f()
}

/// The remainder-boundary switches, checked at start-up: the mode parses.
/// Mode on (and check) with a census dump write the boundary record (DESIGN
/// 6b); the refusal of 6e49bae5e5 is lifted with the dump writer.
pub fn validate_remainder_boundary() -> Result<BoundaryMode> {
    child_boundary_mode()?;
    remainder_boundary_mode()
}

/// `GAMMALOOP_CHILD_BOUNDARY=off|dirac|closed` (default `closed`: production
/// since r6, results/child-boundary/REPORT.md): which closed
/// proper UV children (support roots evaluated as their own R* root, whose
/// momentum-free root closure closes every root-scope remainder) take the
/// remainder-boundary path for their admitted ordinary remainders, as a
/// closed diagram root does.
///
/// - empty, `0` or `off`: none, the path before the production default
///   (every child keeps `not_admitted=no_closed_root`);
/// - `dirac`: children closed by the Dirac closure
///   (`GAMMALOOP_UV_CHILD_CLOSURE=dirac`) of Taylor order 0, the vertex
///   `c gamma^mu` (the degree-1 self-energy keeps the open index `nu` of
///   `Tr[gamma_nu X]` and is not admitted: `child_open_index`);
/// - `closed` (also unset): also children closed by the isotropic closure
///   (`GAMMALOOP_UV_CHILD_CLOSURE=isotropic,dirac`): a degree-0 four-slot
///   child through its three scalar metric channels, whose channel markers
///   are momentum-free scalar spectators of ONE boundary per remainder.
///
/// Acts only with `GAMMALOOP_REMAINDER_BOUNDARY=1|check`. Any other value is
/// an error at the first use (and at start-up).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum ChildBoundary {
    Off,
    Dirac,
    Closed,
}

impl ChildBoundary {
    pub fn as_str(&self) -> &'static str {
        match self {
            ChildBoundary::Off => "off",
            ChildBoundary::Dirac => "dirac",
            ChildBoundary::Closed => "closed",
        }
    }
}

fn parse_child_boundary(value: Option<&str>) -> std::result::Result<ChildBoundary, String> {
    match value.map(str::trim) {
        Some("" | "0" | "off") => Ok(ChildBoundary::Off),
        Some("dirac") => Ok(ChildBoundary::Dirac),
        None | Some("closed") => Ok(ChildBoundary::Closed),
        Some(other) => Err(format!(
            "GAMMALOOP_CHILD_BOUNDARY={other:?}: expected unset, 0, off, dirac or closed"
        )),
    }
}

static CHILD_BOUNDARY: LazyLock<std::result::Result<ChildBoundary, String>> =
    LazyLock::new(|| match std::env::var("GAMMALOOP_CHILD_BOUNDARY") {
        Ok(value) => parse_child_boundary(Some(&value)),
        Err(std::env::VarError::NotPresent) => parse_child_boundary(None),
        Err(std::env::VarError::NotUnicode(value)) => {
            Err(format!("GAMMALOOP_CHILD_BOUNDARY={value:?} is not unicode"))
        }
    });

#[cfg(test)]
thread_local! {
    static TEST_CHILD_BOUNDARY: std::cell::Cell<Option<ChildBoundary>> =
        const { std::cell::Cell::new(None) };
}

/// The mode of `GAMMALOOP_CHILD_BOUNDARY`, read once.
pub fn child_boundary_mode() -> Result<ChildBoundary> {
    #[cfg(test)]
    if let Some(mode) = TEST_CHILD_BOUNDARY.with(|slot| slot.get()) {
        return Ok(mode);
    }
    CHILD_BOUNDARY.clone().map_err(|error| eyre!(error))
}

/// Tests: run `f` with the child-boundary mode set on this thread.
#[cfg(test)]
#[allow(dead_code)]
pub(crate) fn with_child_boundary<T>(mode: ChildBoundary, f: impl FnOnce() -> T) -> T {
    let previous = TEST_CHILD_BOUNDARY.with(|slot| slot.replace(Some(mode)));
    let result = f();
    TEST_CHILD_BOUNDARY.with(|slot| slot.set(previous));
    result
}

#[cfg(test)]
#[test]
fn child_boundary_values_parse() {
    assert_eq!(parse_child_boundary(None), Ok(ChildBoundary::Closed));
    for off in ["", "0", "off", " off "] {
        assert_eq!(parse_child_boundary(Some(off)), Ok(ChildBoundary::Off));
    }
    assert_eq!(
        parse_child_boundary(Some("dirac")),
        Ok(ChildBoundary::Dirac)
    );
    assert_eq!(
        parse_child_boundary(Some("closed")),
        Ok(ChildBoundary::Closed)
    );
    assert!(parse_child_boundary(Some("1")).is_err());
    assert!(ChildBoundary::Dirac < ChildBoundary::Closed);
}

/// `VAKINT_P_BOUNDARY_BLOCKS`: `reuse` (default) or `off`.
pub(crate) fn vakint_boundary_blocks() -> &'static str {
    match std::env::var("VAKINT_P_BOUNDARY_BLOCKS").as_deref() {
        Ok("off" | "0") => "off",
        _ => "reuse",
    }
}

/// The non-identity checkpoint field `{"mode": .., "strict": ..}` (DESIGN 0).
pub fn remainder_boundary_provenance() -> serde_json::Value {
    serde_json::json!({
        "mode": match remainder_boundary_mode() {
            Ok(mode) => mode.as_str().to_owned(),
            Err(error) => format!("error({error})"),
        },
        "strict": boundary_strict(),
    })
}

/// The start-up line of a binary that may take the new path (DESIGN 0).
pub fn startup_line() -> String {
    let mut line = format!(
        "remainder_boundary={} strict={} vakint_gates={} vakint_boundary_blocks={} fallback_max_summands={} fallback_max_bytes={} large_loops={LARGE_LOOPS}",
        match remainder_boundary_mode() {
            Ok(mode) => mode.as_str().to_owned(),
            Err(error) => format!("error({error})"),
        },
        u8::from(boundary_strict()),
        match vakint_gates() {
            Ok(()) => "ok",
            Err(reason) => reason,
        },
        vakint_boundary_blocks(),
        match fallback_max_summands() {
            usize::MAX => "unlimited".to_owned(),
            max => max.to_string(),
        },
        match fallback_max_bytes() {
            usize::MAX => "unlimited".to_owned(),
            max => max.to_string(),
        }
    );
    // Only when on: the line is unchanged with the switch off.
    match child_boundary_mode() {
        Ok(ChildBoundary::Off) => {}
        Ok(mode) => line.push_str(&format!(" child_boundary={}", mode.as_str())),
        Err(error) => line.push_str(&format!(" child_boundary=error({error})")),
    }
    line
}

// ---------------------------------------------------------------------------
// The scalar root normalisation (beta)
// ---------------------------------------------------------------------------

/// The closed metric-root normalisation as data: for every closed,
/// colour-free, index-free scalar `X` in dots, the root normalisation callback
/// gives `factor · X.replace_multiple(&external_routing)` (by value).
#[derive(Clone, Debug)]
pub struct ClosedRootNormalization {
    /// The physical external edge `x0` (beta's `p`).
    pub external: EdgeIndex,
    /// `σ`: `Q(e, args) -> ±Q(x0, args)` for both external edges.
    pub external_routing: Vec<Replacement>,
    /// `c`, with the scalar dimension `dim` not yet substituted.
    pub factor: Atom,
}

// ---------------------------------------------------------------------------
// Declines and counts
// ---------------------------------------------------------------------------

/// A remainder admitted to the new path that the build cannot represent.
#[derive(Clone, Debug)]
pub(super) struct Decline {
    pub stage: &'static str,
    pub reason: String,
}

impl Decline {
    pub(super) fn new(stage: &'static str, reason: impl Into<String>) -> Self {
        Self {
            stage,
            reason: reason.into(),
        }
    }
}

/// The counts of one built remainder (`BOUNDARY_STATS`, DESIGN 8).
#[derive(Clone, Debug, Default)]
pub(super) struct BoundaryStats {
    pub loops: usize,
    pub degree: usize,
    pub p_summands: usize,
    pub p_bytes: usize,
    pub p_truncated: usize,
    pub certified: bool,
    pub notation_retries: usize,
    pub grades: usize,
    pub graded_entries: usize,
    pub d_orders: Vec<usize>,
    pub raw_patterns: usize,
    pub patterns: usize,
    pub zero_patterns: usize,
    pub classes: usize,
    pub blocks: usize,
    pub block_max: (usize, usize),
    pub blocks_total: (usize, usize),
    pub r_terms: usize,
    pub r_summands_total: usize,
    pub r_summands_max: usize,
    pub images: usize,
    pub routing_terms: usize,
    pub block_reads: usize,
    pub block_hits: usize,
    /// Whether the numerator was contracted (the certificate is meaningful).
    pub contracted: bool,
    /// Stage name and seconds, in order.
    pub seconds: Vec<(&'static str, f64)>,
    /// RSS (GiB) after the contraction, the split, the build and the evaluation.
    pub rss_gb: [f64; 4],
}

/// `BOUNDARY_STATS`: one stderr line per built remainder (DESIGN 8).
pub(super) fn log_stats(graph: &str, id: usize, stats: &BoundaryStats) {
    eprintln!(
        "BOUNDARY_STATS graph={graph} id={id} loops={} deg={} P={}/{} P_truncated={} certified={} notation_retries={} grades={} graded_entries={} D_j={} raw_patterns={} patterns={} zero_patterns={} classes={} blocks={} block_max={}/{} blocks_total={}/{} R={}/{}/{} images={} routing_terms={} block_reads={} block_hits={} boundary_blocks={} seconds={} rss_gb={}",
        stats.loops,
        stats.degree,
        stats.p_summands,
        stats.p_bytes,
        stats.p_truncated,
        stats.certified,
        stats.notation_retries,
        stats.grades,
        stats.graded_entries,
        stats
            .d_orders
            .iter()
            .map(|n| n.to_string())
            .collect::<Vec<_>>()
            .join(","),
        stats.raw_patterns,
        stats.patterns,
        stats.zero_patterns,
        stats.classes,
        stats.blocks,
        stats.block_max.0,
        stats.block_max.1,
        stats.blocks_total.0,
        stats.blocks_total.1,
        stats.r_terms,
        stats.r_summands_total,
        stats.r_summands_max,
        stats.images,
        stats.routing_terms,
        stats.block_reads,
        stats.block_hits,
        vakint_boundary_blocks(),
        stats
            .seconds
            .iter()
            .map(|(stage, seconds)| format!("{stage}:{seconds:.3}"))
            .collect::<Vec<_>>()
            .join(","),
        stats
            .rss_gb
            .iter()
            .map(|gb| format!("{gb:.2}"))
            .collect::<Vec<_>>()
            .join(","),
    );
}

/// The remainder-boundary counts of one evaluation (`BOUNDARY_SUMMARY`,
/// DESIGN 8).
#[derive(Clone, Debug, Default)]
pub(super) struct BoundarySummary {
    /// Ordinary remainders seen with the mode on or check.
    pub remainders: usize,
    /// Built boundaries (evaluated or checked).
    pub boundary: usize,
    /// Built boundaries that vanish (scaleless or zero).
    pub zero: usize,
    /// Declines after admission (today's path ran).
    pub fallback: usize,
    /// Admitted remainders whose processing ended in an error (large or
    /// strict declines, check mismatches, hard errors of the build, the
    /// evaluation or today's path).
    pub errors: usize,
    pub not_admitted: BTreeMap<String, usize>,
    pub fallback_reasons: BTreeMap<String, usize>,
    /// Contracted numerators with and without idenso's closed-scalar
    /// certificate (DESIGN Addendum A(2)).
    pub certified: usize,
    pub uncertified: usize,
    pub seconds_boundary: f64,
    /// Every today's-path outer integration (ordinary, IR, fallbacks, checks).
    pub seconds_today: f64,
    /// vakint::timing over the same intervals as `seconds_boundary` and
    /// `seconds_today`: reduction (with finalisation) and lowering seconds.
    pub reduce_boundary: f64,
    pub lower_boundary: f64,
    pub reduce_today: f64,
    pub lower_today: f64,
    /// vakint::timing reduction entries (key-plan reductions) over the
    /// same intervals (`CHILD_BOUNDARY` lines only).
    pub reduce_calls_boundary: u64,
    pub reduce_calls_today: u64,
    /// Check mode: remainders compared, and those with a skipped input pair.
    pub checks: usize,
    pub checks_skipped: usize,
    /// `GAMMALOOP_UNCONTRACTED_NUMERATOR=on`: remainders vakint's factorised
    /// entry took (today's path did not run), their wall seconds and
    /// vakint::timing counters (U7a review M3; printed only when nonzero).
    pub uncontracted: usize,
    pub seconds_uncontracted: f64,
    pub reduce_uncontracted: f64,
    pub lower_uncontracted: f64,
}

impl BoundarySummary {
    /// Add the vakint::timing counters of an uncontracted record's entry.
    pub(super) fn add_uncontracted(&mut self, delta: vakint::timing::Delta) {
        self.reduce_uncontracted += delta.reduce_seconds;
        self.lower_uncontracted += delta.lower_seconds;
    }

    /// Add the vakint::timing counters of a boundary-path interval.
    pub(super) fn add_boundary(&mut self, delta: vakint::timing::Delta) {
        self.reduce_boundary += delta.reduce_seconds;
        self.lower_boundary += delta.lower_seconds;
        self.reduce_calls_boundary += delta.reduce_calls;
    }

    /// Add the vakint::timing counters of a today's-path interval.
    pub(super) fn add_today(&mut self, delta: vakint::timing::Delta) {
        self.reduce_today += delta.reduce_seconds;
        self.lower_today += delta.lower_seconds;
        self.reduce_calls_today += delta.reduce_calls;
    }

    pub(super) fn declined(&mut self, decline: &Decline) {
        self.fallback += 1;
        *self
            .fallback_reasons
            .entry(format!("{}:{}", decline.stage, decline.reason))
            .or_default() += 1;
    }

    pub(super) fn contracted(&mut self, stats: &BoundaryStats) {
        if stats.contracted {
            if stats.certified {
                self.certified += 1;
            } else {
                self.uncertified += 1;
            }
        }
    }

    /// `outcome`: `ok`, or `error` when the evaluation ends in an error (the
    /// summary is printed before the error propagates).
    /// `eval_seconds`, `eval`: the wall seconds and the vakint::timing counters
    /// of the whole evaluation so far (its UV children, IR parts and outer
    /// integrations; `eval_s`, `reduce_eval_s`, `lower_eval_s`).
    pub(super) fn log(
        &self,
        graph: &str,
        root: &str,
        outcome: &str,
        eval_seconds: f64,
        eval: vakint::timing::Delta,
    ) {
        let list = |map: &BTreeMap<String, usize>| {
            if map.is_empty() {
                "none".to_owned()
            } else {
                map.iter()
                    .map(|(reason, count)| format!("{reason}:{count}"))
                    .collect::<Vec<_>>()
                    .join(",")
            }
        };
        eprintln!(
            "BOUNDARY_SUMMARY graph={graph} root={root} remainders={} boundary={} zero={} fallback={} errors={} not_admitted={} fallback_reasons={} certified={} uncertified={} checks={} checks_skipped={} seconds_boundary={:.3} seconds_today={:.3} reduce_boundary_s={:.3} lower_boundary_s={:.3} reduce_today_s={:.3} lower_today_s={:.3} eval_s={eval_seconds:.3} reduce_eval_s={:.3} lower_eval_s={:.3}{} outcome={outcome}",
            self.remainders,
            self.boundary,
            self.zero,
            self.fallback,
            self.errors,
            list(&self.not_admitted),
            list(&self.fallback_reasons),
            self.certified,
            self.uncertified,
            self.checks,
            self.checks_skipped,
            self.seconds_boundary,
            self.seconds_today,
            self.reduce_boundary,
            self.lower_boundary,
            self.reduce_today,
            self.lower_today,
            eval.reduce_seconds,
            eval.lower_seconds,
            if self.uncontracted > 0 {
                format!(
                    " uncontracted={} seconds_uncontracted={:.3} reduce_uncontracted_s={:.3} lower_uncontracted_s={:.3}",
                    self.uncontracted,
                    self.seconds_uncontracted,
                    self.reduce_uncontracted,
                    self.lower_uncontracted
                )
            } else {
                String::new()
            },
        );
        #[cfg(test)]
        LOGGED_SUMMARIES.with(|logged| {
            logged
                .borrow_mut()
                .push((root.to_owned(), outcome.to_owned(), self.clone()))
        });
    }
}

#[cfg(test)]
thread_local! {
    static LOGGED_SUMMARIES: std::cell::RefCell<Vec<(String, String, BoundarySummary)>> =
        const { std::cell::RefCell::new(Vec::new()) };
}

/// Tests: the summaries logged on this thread since the last call, as
/// `(root, outcome, summary)`.
#[cfg(test)]
pub(super) fn take_logged_summaries() -> Vec<(String, String, BoundarySummary)> {
    LOGGED_SUMMARIES.with(|logged| std::mem::take(&mut *logged.borrow_mut()))
}

/// The process resident set size in GiB (`/proc/self/statm`).
pub(super) fn rss_gb() -> f64 {
    std::fs::read_to_string("/proc/self/statm")
        .ok()
        .and_then(|statm| statm.split_whitespace().nth(1)?.parse::<u64>().ok())
        .map_or(0.0, |pages| (pages * 4096) as f64 / (1u64 << 30) as f64)
}

/// The stage clock of one build: every stage end is logged as
/// `BOUNDARY_STAGE graph= id= stage= seconds= rss_gb=`.
pub(super) struct StageClock<'a> {
    graph: &'a str,
    id: usize,
    last: Instant,
}

impl<'a> StageClock<'a> {
    pub(super) fn new(graph: &'a str, id: usize) -> Self {
        Self {
            graph,
            id,
            last: Instant::now(),
        }
    }

    /// End `stage`: record and log its seconds and the RSS.
    pub(super) fn lap(&mut self, stats: &mut BoundaryStats, stage: &'static str) -> f64 {
        let now = Instant::now();
        let seconds = now.duration_since(self.last).as_secs_f64();
        self.last = now;
        let rss = rss_gb();
        eprintln!(
            "BOUNDARY_STAGE graph={} id={} stage={stage} seconds={seconds:.3} rss_gb={rss:.2}",
            self.graph, self.id
        );
        stats.seconds.push((stage, seconds));
        seconds
    }
}

// ---------------------------------------------------------------------------
// Grades and the recursive grammar (DESIGN 2a)
// ---------------------------------------------------------------------------

/// `g = (T, s, b)`: the sorted basis positions of the `k·p` factors, the power
/// of `p²` and the power of the reference mass squared.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Hash, Default)]
pub(super) struct Grade {
    pub structure: Vec<usize>,
    pub p_squared: usize,
    pub mass_power: usize,
}

impl Grade {
    /// The λ-degree `a(g) = |T| + 2 s`.
    pub(super) fn degree(&self) -> usize {
        self.structure.len() + 2 * self.p_squared
    }

    fn is_zero(&self) -> bool {
        self.structure.is_empty() && self.p_squared == 0 && self.mass_power == 0
    }

    fn times(&self, other: &Grade) -> Grade {
        let mut structure = self.structure.clone();
        structure.extend_from_slice(&other.structure);
        structure.sort_unstable();
        Grade {
            structure,
            p_squared: self.p_squared + other.p_squared,
            mass_power: self.mass_power + other.mass_power,
        }
    }

    /// `t_{T,s} = Π_{j∈T} dot(k(j+1), p) · dot(p, p)^s` in vakint symbols.
    pub(super) fn structure_atom(&self, p: &Atom) -> Atom {
        let mut factors = self
            .structure
            .iter()
            .map(|j| function!(S.dot, function!(S.k, *j as i64 + 1), p))
            .collect::<Vec<_>>();
        if self.p_squared > 0 {
            factors.push(function!(S.dot, p, p).pow(self.p_squared as i64));
        }
        Atom::mul_many(factors.iter())
    }
}

/// A small map grade -> coefficient (p-free, reference-mass-free, vakint
/// spelling).
pub(super) type GradedMap = BTreeMap<Grade, Atom>;

/// Why a view has no graded map: a decline, an index slot that the
/// notation retry may still remove, or a hard error (an exhausted resource
/// budget), which is never a decline.
#[derive(Debug)]
enum Fail {
    Decline(Decline),
    Indexed(String),
    Error(eyre::Report),
}

/// Grammar steps (sum terms visited, convolution products) between two
/// resource checkpoints (DESIGN 9).
const CHECKPOINT_STEPS: usize = 4096;

/// A resource checkpoint of the boundary build, as a hard error.
fn checkpoint(stage: &'static str) -> Result<()> {
    symbolica::utils::resources::checkpoint()
        .map_err(|error| eyre!("{error}"))
        .wrap_err_with(|| format!("while building an R* remainder boundary ({stage})"))
}

fn decline<T>(reason: impl Into<String>) -> std::result::Result<T, Fail> {
    Err(Fail::Decline(Decline::new("collect", reason)))
}

/// A vector of a `dot` argument: the basis loop momentum at position `j`, or
/// the external `x0`.
#[derive(Clone, Copy, Debug)]
enum Vector {
    Loop(usize),
    External,
}

fn head_name(symbol: Symbol) -> String {
    symbol
        .get_name()
        .rsplit("::")
        .next()
        .unwrap_or("")
        .to_owned()
}

/// The λ symbols of the outer and hard Taylor operations.
pub(super) fn lambda_symbols() -> [Symbol; 2] {
    [
        symbol!("rstar_outer_taylor"; Scalar),
        symbol!("rstar_hard_taylor"; Scalar),
    ]
}

/// The symbols that must never reach vakint in a block, coefficient or
/// routing image (DESIGN 2d), with the reference mass squared.
fn forbidden_symbols(reference_mass_squared: Symbol) -> Vec<Symbol> {
    let minkowski = Minkowski {}
        .new_rep(GS.dim)
        .to_symbolic([])
        .get_symbol()
        .expect("the Minkowski representation is a function");
    let ags = &idenso::dirac::AGS;
    let cs = &idenso::color::CS;
    let mut symbols = vec![
        GS.emr_mom,
        GS.loop_mom,
        SPENSO_TAG.dot,
        minkowski,
        ETS.metric,
        ags.gamma,
        ags.gammaconj,
        ags.gammaadj,
        ags.gamma0,
        ags.projp,
        ags.projm,
        ags.gamma5,
        ags.charge_conjugation,
        ags.sigma,
        SPENSO_TAG.trace,
        SPENSO_TAG.chain,
        cs.cas,
        cs.idx,
        cs.gram,
        GS.den,
        GS.dim,
        reference_mass_squared,
    ];
    symbols.extend(lambda_symbols());
    symbols
}

/// The grammar of one remainder: its basis positions, the external `x0`, the
/// truncation degree and the caches of the split.
pub(super) struct Grammar {
    position: HashMap<usize, usize>,
    x0: usize,
    degree: usize,
    reference_mass_squared: Symbol,
    lambdas: [Symbol; 2],
    forbidden: Vec<Symbol>,
    minkowski: Symbol,
    regulated_dimension: Atom,
    p: Atom,
    cache: HashMap<Atom, GradedMap>,
    retry_images: HashMap<Atom, Atom>,
    matcher: SlotMatcher,
    /// Grammar steps since the start of the split ([`CHECKPOINT_STEPS`]).
    work: usize,
}

/// Non-sum factors up to this size are cached by atom (as `map_factors_cached`).
const CACHE_BYTES: usize = 4096;

impl Grammar {
    /// `basis`: the R* basis loop edges `b_j` in `j` order; `x0`: the external
    /// edge; `degree`: the outer Taylor degree.
    pub(super) fn new(
        basis: &[usize],
        x0: usize,
        degree: usize,
        reference_mass_squared: Symbol,
    ) -> Self {
        Self {
            position: basis
                .iter()
                .enumerate()
                .map(|(j, edge)| (*edge, j))
                .collect(),
            x0,
            degree,
            reference_mass_squared,
            lambdas: lambda_symbols(),
            forbidden: forbidden_symbols(reference_mass_squared),
            minkowski: Minkowski {}
                .new_rep(GS.dim)
                .to_symbolic([])
                .get_symbol()
                .expect("the Minkowski representation is a function"),
            regulated_dimension: Atom::num(4) - Atom::num(2) * Atom::var(GS.dim_epsilon),
            p: function!(S.p, x0 as i64),
            cache: HashMap::new(),
            retry_images: HashMap::new(),
            matcher: SlotMatcher::default(),
            work: 0,
        }
    }

    /// One grammar step; every [`CHECKPOINT_STEPS`] steps a resource
    /// checkpoint, whose failure is a hard error.
    fn tick(&mut self) -> std::result::Result<(), Fail> {
        self.work += 1;
        if self.work.is_multiple_of(CHECKPOINT_STEPS) {
            checkpoint("collect").map_err(Fail::Error)?;
        }
        Ok(())
    }

    /// The vacuum label `k(j+1)` of basis position `j`.
    fn loop_vector(j: usize) -> Atom {
        function!(S.k, j as i64 + 1)
    }

    /// The external vector `p(x0)` in vakint symbols.
    pub(super) fn external(&self) -> &Atom {
        &self.p
    }

    /// `Q(e)` or `Q(e, mink(d))` (an index-free payload): the vector.
    fn vector(&self, view: AtomView<'_>) -> std::result::Result<Vector, Fail> {
        let AtomView::Fun(function) = view else {
            return decline(format!("dot_argument:{}", kind(view)));
        };
        if function.get_symbol() != GS.emr_mom {
            return decline(format!("dot_argument:{}", head_name(function.get_symbol())));
        }
        let nargs = function.get_nargs();
        if nargs == 0 || nargs > 2 {
            return decline("dot_argument:Q");
        }
        if nargs == 2 {
            match function.get(1) {
                AtomView::Fun(slot)
                    if slot.get_symbol() == self.minkowski && slot.get_nargs() == 1 => {}
                _ => return decline("dot_argument:indexed_Q"),
            }
        }
        let Ok(edge) = usize::try_from(function.get(0)) else {
            return decline("dot_argument:Q");
        };
        if let Some(j) = self.position.get(&edge) {
            Ok(Vector::Loop(*j))
        } else if edge == self.x0 {
            Ok(Vector::External)
        } else {
            decline(format!("dot_argument:Q({edge})"))
        }
    }

    /// A `dot` argument as `Σ α_e v_e` with numbers `α_e`.
    fn linear_vectors(&self, view: AtomView<'_>) -> std::result::Result<Vec<(Atom, Vector)>, Fail> {
        let terms = match view {
            AtomView::Add(sum) => sum.iter().collect::<Vec<_>>(),
            other => vec![other],
        };
        let mut vectors = Vec::with_capacity(terms.len());
        for term in terms {
            let (coefficient, vector) = match term {
                AtomView::Mul(product) => {
                    let factors = product.iter().collect::<Vec<_>>();
                    match factors.as_slice() {
                        [n @ AtomView::Num(_), v] | [v, n @ AtomView::Num(_)] => (n.to_owned(), *v),
                        _ => return decline("dot_argument:nonlinear"),
                    }
                }
                other => (Atom::num(1), other),
            };
            vectors.push((coefficient, self.vector(vector)?));
        }
        Ok(vectors)
    }

    /// The bilinear expansion of a spenso `dot(u, v)` into leaves. Graded:
    /// `k·p` and `p²` become grades; ungraded (`map_scalar`): vakint dots.
    fn dot(
        &self,
        left: AtomView<'_>,
        right: AtomView<'_>,
        graded: bool,
        truncated: &mut usize,
    ) -> std::result::Result<GradedMap, Fail> {
        let left = self.linear_vectors(left)?;
        let right = self.linear_vectors(right)?;
        let mut accumulated = BTreeMap::<Grade, Vec<Atom>>::new();
        for (a, u) in &left {
            for (b, v) in &right {
                let coefficient = a * b;
                let (grade, leaf) = match (u, v) {
                    (Vector::Loop(j), Vector::Loop(l)) => (
                        Grade::default(),
                        function!(S.dot, Self::loop_vector(*j), Self::loop_vector(*l)),
                    ),
                    (Vector::Loop(j), Vector::External) | (Vector::External, Vector::Loop(j)) => {
                        if graded {
                            (
                                Grade {
                                    structure: vec![*j],
                                    ..Default::default()
                                },
                                Atom::num(1),
                            )
                        } else {
                            (
                                Grade::default(),
                                function!(S.dot, Self::loop_vector(*j), &self.p),
                            )
                        }
                    }
                    (Vector::External, Vector::External) => {
                        if graded {
                            (
                                Grade {
                                    p_squared: 1,
                                    ..Default::default()
                                },
                                Atom::num(1),
                            )
                        } else {
                            (Grade::default(), function!(S.dot, &self.p, &self.p))
                        }
                    }
                };
                if grade.degree() > self.degree {
                    *truncated += 1;
                    continue;
                }
                accumulated
                    .entry(grade)
                    .or_default()
                    .push(coefficient * leaf);
            }
        }
        Ok(collect_map(accumulated))
    }

    /// Whether a function (outside a `dot`) carries an explicit index slot or
    /// an opaque representation, or is an indexed head (metric, gamma, trace,
    /// chain): the notation retry applies.
    fn indexed(&mut self, view: AtomView<'_>) -> Option<String> {
        let indexed_heads = [
            ETS.metric,
            idenso::dirac::AGS.gamma,
            idenso::dirac::AGS.gamma5,
            idenso::dirac::AGS.gamma0,
            SPENSO_TAG.trace,
            SPENSO_TAG.chain,
        ];
        let mut found = None;
        let matcher = &mut self.matcher;
        view.visitor(&mut |part| {
            if found.is_some() {
                return false;
            }
            let AtomView::Fun(function) = part else {
                return true;
            };
            let symbol = function.get_symbol();
            if symbol == SPENSO_TAG.dot {
                return false;
            }
            if indexed_heads.contains(&symbol) {
                found = Some(head_name(symbol));
                return false;
            }
            match matcher.classify(part) {
                SlotMatch::Explicit(_) | SlotMatch::Opaque => {
                    found = Some(head_name(symbol));
                    false
                }
                SlotMatch::Other => true,
            }
        });
        found
    }

    /// The graded map of one view (DESIGN 2a, the recursive grammar).
    fn grade(
        &mut self,
        view: AtomView<'_>,
        stats: &mut BoundaryStats,
    ) -> std::result::Result<GradedMap, Fail> {
        match view {
            AtomView::Num(_) => Ok(single(Grade::default(), view.to_owned())),
            AtomView::Var(variable) => {
                let symbol = variable.get_symbol();
                if symbol == self.reference_mass_squared {
                    Ok(single(
                        Grade {
                            mass_power: 1,
                            ..Default::default()
                        },
                        Atom::num(1),
                    ))
                } else if symbol == GS.dim {
                    Ok(single(Grade::default(), self.regulated_dimension.clone()))
                } else if self.lambdas.contains(&symbol) {
                    decline("lambda")
                } else if self.forbidden.contains(&symbol) {
                    decline(format!("symbol:{}", head_name(symbol)))
                } else {
                    Ok(single(Grade::default(), view.to_owned()))
                }
            }
            AtomView::Add(sum) => {
                let mut accumulated = BTreeMap::<Grade, Vec<Atom>>::new();
                for term in sum.iter() {
                    self.tick()?;
                    for (grade, coefficient) in self.grade(term, stats)? {
                        accumulated.entry(grade).or_default().push(coefficient);
                    }
                }
                Ok(collect_map(accumulated))
            }
            AtomView::Mul(product) => {
                let mut factors = Vec::with_capacity(product.get_nargs());
                for factor in product.iter() {
                    factors.push(self.factor(factor, stats)?);
                }
                self.product(factors, stats)
            }
            AtomView::Pow(_) | AtomView::Fun(_) => self.factor(view, stats),
        }
    }

    /// The graded map of a non-sum factor, cached by atom up to 4 KiB.
    fn factor(
        &mut self,
        view: AtomView<'_>,
        stats: &mut BoundaryStats,
    ) -> std::result::Result<GradedMap, Fail> {
        match view {
            AtomView::Pow(_) | AtomView::Fun(_) if view.get_byte_size() <= CACHE_BYTES => {
                let key = view.to_owned();
                if let Some(map) = self.cache.get(&key) {
                    return Ok(map.clone());
                }
                let map = self.uncached_factor(view, stats)?;
                self.cache.insert(key, map.clone());
                Ok(map)
            }
            AtomView::Pow(_) | AtomView::Fun(_) => self.uncached_factor(view, stats),
            other => self.grade(other, stats),
        }
    }

    fn uncached_factor(
        &mut self,
        view: AtomView<'_>,
        stats: &mut BoundaryStats,
    ) -> std::result::Result<GradedMap, Fail> {
        match view {
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let base_map = self.grade(base, stats)?;
                let integer = i64::try_from(exponent).ok();
                if base_map.is_empty() {
                    // A vanishing (or wholly truncated) base.
                    return match integer {
                        Some(n) if n >= 1 => Ok(GradedMap::new()),
                        _ => decline("graded_power"),
                    };
                }
                if base_map.keys().all(Grade::is_zero) {
                    let exponent = self.scalar(exponent)?;
                    let (grade, coefficient) = base_map.into_iter().next().unwrap();
                    return Ok(single(grade, coefficient.pow(&exponent)));
                }
                match integer {
                    Some(n) if n >= 1 => self.power(&base_map, n as u64, stats),
                    Some(_) => decline("graded_power"),
                    None if base_map.keys().any(|grade| grade.mass_power > 0) => {
                        decline("fractional_reference_mass_power")
                    }
                    None => decline("graded_power"),
                }
            }
            AtomView::Fun(function) => {
                let symbol = function.get_symbol();
                if symbol == SPENSO_TAG.dot {
                    if function.get_nargs() != 2 {
                        return decline("dot_argument:arity");
                    }
                    return self.dot(
                        function.get(0),
                        function.get(1),
                        true,
                        &mut stats.p_truncated,
                    );
                }
                if let Some(head) = self.indexed(view) {
                    return Err(Fail::Indexed(head));
                }
                if self.forbidden.contains(&symbol) {
                    return decline(format!("function:{}", head_name(symbol)));
                }
                // No argument may hold a momentum, the reference mass, a λ or
                // any other forbidden head (the scalar `dim` is mapped).
                for argument in function.iter() {
                    if self
                        .forbidden
                        .iter()
                        .any(|f| *f != GS.dim && argument.contains_symbol(*f))
                    {
                        return decline(format!("function:{}", head_name(symbol)));
                    }
                }
                let mut builder = FunctionBuilder::new(symbol);
                for argument in function.iter() {
                    builder = builder.add_arg(self.scalar(argument)?);
                }
                Ok(single(Grade::default(), builder.finish()))
            }
            other => self.grade(other, stats),
        }
    }

    /// A view that must be of grade zero only (exponents, function
    /// arguments): its mapped atom.
    fn scalar(&mut self, view: AtomView<'_>) -> std::result::Result<Atom, Fail> {
        let mut stats = BoundaryStats::default();
        let map = self.grade(view, &mut stats)?;
        if map.keys().any(|grade| !grade.is_zero()) {
            return decline("graded_scalar");
        }
        Ok(map.into_values().next().unwrap_or(Atom::Zero))
    }

    /// The convolution of the factors' maps, truncating `a(g) > degree`.
    fn product(
        &mut self,
        factors: Vec<GradedMap>,
        stats: &mut BoundaryStats,
    ) -> std::result::Result<GradedMap, Fail> {
        // The common case: every factor has one grade.
        if factors.iter().all(|map| map.len() == 1) {
            let mut grade = Grade::default();
            let mut coefficients = Vec::with_capacity(factors.len());
            for map in factors {
                let (g, c) = map.into_iter().next().unwrap();
                grade = grade.times(&g);
                coefficients.push(c);
            }
            if grade.degree() > self.degree {
                stats.p_truncated += 1;
                return Ok(GradedMap::new());
            }
            let coefficient = Atom::mul_many(coefficients.iter());
            return Ok(if coefficient.is_zero() {
                GradedMap::new()
            } else {
                single(grade, coefficient)
            });
        }
        let mut accumulated = single(Grade::default(), Atom::num(1));
        for factor in factors {
            accumulated = self.convolve(&accumulated, &factor, stats)?;
            if accumulated.is_empty() {
                break;
            }
        }
        Ok(accumulated)
    }

    fn convolve(
        &mut self,
        a: &GradedMap,
        b: &GradedMap,
        stats: &mut BoundaryStats,
    ) -> std::result::Result<GradedMap, Fail> {
        let mut accumulated = BTreeMap::<Grade, Vec<Atom>>::new();
        for (ga, ca) in a {
            for (gb, cb) in b {
                self.tick()?;
                let grade = ga.times(gb);
                if grade.degree() > self.degree {
                    stats.p_truncated += 1;
                    continue;
                }
                accumulated.entry(grade).or_default().push(ca * cb);
            }
        }
        Ok(collect_map(accumulated))
    }

    /// `base^n` by binary powering, truncated.
    fn power(
        &mut self,
        base: &GradedMap,
        mut n: u64,
        stats: &mut BoundaryStats,
    ) -> std::result::Result<GradedMap, Fail> {
        let mut result = single(Grade::default(), Atom::num(1));
        let mut square = base.clone();
        while n > 0 {
            if n & 1 == 1 {
                result = self.convolve(&result, &square, stats)?;
            }
            n >>= 1;
            if n > 0 {
                square = self.convolve(&square, &square, stats)?;
            }
        }
        Ok(result)
    }

    /// The graded map of one summand of `P`, with the notation retry: the
    /// summand's indexed factors through the R* remainder contraction
    /// ([`retry_contraction`]: cook, infer, `simplify_algebra`, `to_dots`,
    /// uncook; once per distinct product of indexed factors), then graded
    /// again. `Err` is a hard error (a failing retry contraction or an
    /// exhausted resource budget).
    pub(super) fn grade_summand(
        &mut self,
        summand: AtomView<'_>,
        stats: &mut BoundaryStats,
    ) -> Result<std::result::Result<GradedMap, Decline>> {
        let head = match self.grade(summand, stats) {
            Ok(map) => return Ok(Ok(map)),
            Err(Fail::Decline(decline)) => return Ok(Err(decline)),
            Err(Fail::Error(error)) => return Err(error),
            Err(Fail::Indexed(head)) => head,
        };
        stats.notation_retries += 1;
        let notated = match streamed::map_indexed_factors(
            summand,
            &mut self.retry_images,
            &retry_contraction,
        ) {
            Ok(Some(notated)) => notated,
            // The retry would introduce dummy indices: a decline.
            Ok(None) => return Ok(Err(Decline::new("collect", format!("indexed:{head}")))),
            // A failing contraction (cook, infer, simplify_algebra) or an
            // exhausted budget is a hard error, as the remainder contraction
            // itself: never a silent fallback, and its cause is kept.
            Err(error) => {
                return Err(error.wrap_err(format!(
                    "while retrying the contraction of an indexed factor ({head}) of an R* remainder boundary"
                )));
            }
        };
        match self.grade(notated.as_view(), stats) {
            Ok(map) => Ok(Ok(map)),
            Err(Fail::Decline(decline)) => Ok(Err(decline)),
            Err(Fail::Error(error)) => Err(error),
            Err(Fail::Indexed(head)) => Ok(Err(Decline::new("collect", format!("indexed:{head}")))),
        }
    }

    /// The leaf rules without grading (DESIGN 2b): `p(x0)` is allowed, dots
    /// of linear combinations are expanded bilinearly, `dim -> 4-2ε`. For the
    /// denominator coefficients and the normalisation.
    pub(super) fn map_scalar(
        &mut self,
        view: AtomView<'_>,
    ) -> Result<std::result::Result<Atom, Decline>> {
        match self.map_scalar_inner(view) {
            Ok(atom) => Ok(Ok(atom)),
            Err(Fail::Decline(decline)) => Ok(Err(Decline::new("project", decline.reason))),
            Err(Fail::Indexed(head)) => Ok(Err(Decline::new("project", format!("indexed:{head}")))),
            Err(Fail::Error(error)) => Err(error),
        }
    }

    fn map_scalar_inner(&mut self, view: AtomView<'_>) -> std::result::Result<Atom, Fail> {
        match view {
            AtomView::Num(_) => Ok(view.to_owned()),
            AtomView::Var(variable) => {
                let symbol = variable.get_symbol();
                if symbol == GS.dim {
                    Ok(self.regulated_dimension.clone())
                } else if self.forbidden.contains(&symbol) {
                    decline(format!("symbol:{}", head_name(symbol)))
                } else {
                    Ok(view.to_owned())
                }
            }
            AtomView::Add(sum) => {
                let mut terms = Vec::with_capacity(sum.get_nargs());
                for term in sum.iter() {
                    terms.push(self.map_scalar_inner(term)?);
                }
                Ok(Atom::add_many(terms.iter()))
            }
            AtomView::Mul(product) => {
                let mut factors = Vec::with_capacity(product.get_nargs());
                for factor in product.iter() {
                    factors.push(self.map_scalar_inner(factor)?);
                }
                Ok(Atom::mul_many(factors.iter()))
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                Ok(self
                    .map_scalar_inner(base)?
                    .pow(&self.map_scalar_inner(exponent)?))
            }
            AtomView::Fun(function) => {
                let symbol = function.get_symbol();
                if symbol == SPENSO_TAG.dot && function.get_nargs() == 2 {
                    // Ungraded leaves all have grade zero: nothing is truncated.
                    let mut truncated = 0;
                    let map = self.dot(function.get(0), function.get(1), false, &mut truncated)?;
                    return Ok(map.into_values().next().unwrap_or(Atom::Zero));
                }
                if self.forbidden.contains(&symbol)
                    || function.iter().any(|argument| {
                        argument.contains_symbol(GS.emr_mom)
                            || argument.contains_symbol(GS.loop_mom)
                    })
                {
                    return decline(format!("function:{}", head_name(symbol)));
                }
                let mut builder = FunctionBuilder::new(symbol);
                for argument in function.iter() {
                    builder = builder.add_arg(self.map_scalar_inner(argument)?);
                }
                Ok(builder.finish())
            }
        }
    }
}

/// The notation retry's map of a product of indexed factors: the R*
/// remainder contraction itself (cook, infer, `simplify_algebra`, `to_dots`,
/// uncook; [`RStarIntegrand::simplify_view_certified`]). Deviation from
/// DESIGN 2a, which names the dot-notation pass alone (cook, infer,
/// `to_dots`): that pass leaves a closed vector pair `Q(a,μ) Q(b,μ)`
/// unchanged, so the design's own case (c) would decline; both maps keep the
/// value, the contraction also closes such pairs.
pub(super) fn retry_contraction(numerator: AtomView<'_>) -> Result<Atom> {
    Ok(RStarIntegrand::simplify_view_certified(numerator)?.0)
}

fn kind(view: AtomView<'_>) -> &'static str {
    match view {
        AtomView::Num(_) => "number",
        AtomView::Var(_) => "symbol",
        AtomView::Add(_) => "sum",
        AtomView::Mul(_) => "product",
        AtomView::Pow(_) => "power",
        AtomView::Fun(_) => "function",
    }
}

fn single(grade: Grade, coefficient: Atom) -> GradedMap {
    let mut map = GradedMap::new();
    if !coefficient.is_zero() {
        map.insert(grade, coefficient);
    }
    map
}

fn collect_map(accumulated: BTreeMap<Grade, Vec<Atom>>) -> GradedMap {
    accumulated
        .into_iter()
        .filter_map(|(grade, mut coefficients)| {
            let coefficient = if coefficients.len() == 1 {
                coefficients.pop().unwrap()
            } else {
                Atom::add_many(coefficients.iter())
            };
            (!coefficient.is_zero()).then_some((grade, coefficient))
        })
        .collect()
}

/// The graded split of `P` (consumed): every summand's pieces pushed into the
/// block of their grade, then one sum per block. Blocks that vanish are
/// dropped. `Err(decline)` when a summand has no graded map.
pub(super) fn split_into_blocks(
    p: Atom,
    grammar: &mut Grammar,
    stats: &mut BoundaryStats,
) -> Result<std::result::Result<BTreeMap<Grade, Atom>, Decline>> {
    let mut vectors = BTreeMap::<Grade, Vec<Atom>>::new();
    {
        let mut emit = |index: usize, summand: AtomView<'_>| -> Result<Option<Decline>> {
            if index % CHECKPOINT_STEPS == CHECKPOINT_STEPS - 1 {
                checkpoint("collect")?;
            }
            let map = match grammar.grade_summand(summand, stats)? {
                Ok(map) => map,
                Err(decline) => return Ok(Some(decline)),
            };
            for (grade, coefficient) in map {
                stats.graded_entries += 1;
                vectors.entry(grade).or_default().push(coefficient);
            }
            Ok(None)
        };
        match p.as_view() {
            AtomView::Add(sum) => {
                for (index, summand) in sum.iter().enumerate() {
                    if let Some(decline) = emit(index, summand)? {
                        return Ok(Err(decline));
                    }
                }
            }
            single => {
                if let Some(decline) = emit(0, single)? {
                    return Ok(Err(decline));
                }
            }
        }
    }
    drop(p);
    let mut blocks = BTreeMap::new();
    while let Some((grade, coefficients)) = vectors.pop_first() {
        checkpoint("collect: block sums")?;
        let block = Atom::add_many(coefficients.iter());
        drop(coefficients);
        if !block.is_zero() {
            blocks.insert(grade, block);
        }
    }
    stats.grades = blocks.len();
    Ok(Ok(blocks))
}

// ---------------------------------------------------------------------------
// Postconditions (DESIGN 2d)
// ---------------------------------------------------------------------------

/// Where a scalar goes in the boundary.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Role {
    /// A block `S_g`.
    Block,
    /// A coefficient `R_{π,g}`.
    Coefficient,
    /// A routing image.
    Routing,
}

impl Role {
    fn as_str(&self) -> &'static str {
        match self {
            Role::Block => "block",
            Role::Coefficient => "coefficient",
            Role::Routing => "routing",
        }
    }
}

/// The vakint contract of one block, coefficient or routing image: nothing
/// forbidden anywhere; `k(n)` only as a direct `dot` argument with
/// `1 <= n <= loops` (blocks, coefficients) or as a linear term with
/// `1 <= n < loops` (routing); `p(n)`: none in blocks, only `p(x0)` inside a
/// `dot` in coefficients, only the closing head `p(head)` in routing images.
pub(super) fn verify_vakint_scalar(
    view: AtomView<'_>,
    role: Role,
    loops: usize,
    x0: usize,
    head: usize,
    reference_mass_squared: Symbol,
) -> std::result::Result<(), Decline> {
    let forbidden = forbidden_symbols(reference_mass_squared);
    let fail =
        |what: String| Decline::new("verify", format!("postcondition:{}:{what}", role.as_str()));
    let vector_index = |view: AtomView<'_>, symbol: Symbol| -> Option<i64> {
        match view {
            AtomView::Fun(function)
                if function.get_symbol() == symbol && function.get_nargs() == 1 =>
            {
                i64::try_from(function.get(0)).ok()
            }
            _ => None,
        }
    };
    if role == Role::Routing {
        // A linear combination of cut loop momenta and the closing head.
        let terms = match view {
            AtomView::Add(sum) => sum.iter().collect::<Vec<_>>(),
            other => vec![other],
        };
        for term in terms {
            let vector = match term {
                AtomView::Mul(product) => {
                    let factors = product.iter().collect::<Vec<_>>();
                    match factors.as_slice() {
                        [AtomView::Num(_), v] | [v, AtomView::Num(_)] => *v,
                        _ => return Err(fail("nonlinear".into())),
                    }
                }
                other => other,
            };
            if let Some(n) = vector_index(vector, S.k) {
                if n < 1 || n as usize >= loops {
                    return Err(fail(format!("k({n})")));
                }
            } else if let Some(n) = vector_index(vector, S.p) {
                if n as usize != head {
                    return Err(fail(format!("p({n})")));
                }
            } else {
                return Err(fail(kind(vector).into()));
            }
        }
        return Ok(());
    }
    fn walk(
        view: AtomView<'_>,
        inside_dot: bool,
        check: &mut dyn FnMut(AtomView<'_>, bool) -> std::result::Result<bool, Decline>,
    ) -> std::result::Result<(), Decline> {
        if !check(view, inside_dot)? {
            return Ok(());
        }
        match view {
            AtomView::Num(_) | AtomView::Var(_) => Ok(()),
            AtomView::Add(sum) => sum.iter().try_for_each(|v| walk(v, false, check)),
            AtomView::Mul(product) => product.iter().try_for_each(|v| walk(v, false, check)),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                walk(base, false, check)?;
                walk(exponent, false, check)
            }
            AtomView::Fun(function) => {
                let dot = function.get_symbol() == S.dot;
                function.iter().try_for_each(|v| walk(v, dot, check))
            }
        }
    }
    let mut check = |view: AtomView<'_>, inside_dot: bool| -> std::result::Result<bool, Decline> {
        match view {
            AtomView::Var(variable) => {
                let symbol = variable.get_symbol();
                if forbidden.contains(&symbol) {
                    return Err(fail(head_name(symbol)));
                }
                Ok(false)
            }
            AtomView::Fun(function) => {
                let symbol = function.get_symbol();
                if forbidden.contains(&symbol) {
                    return Err(fail(head_name(symbol)));
                }
                if symbol == S.dot && function.get_nargs() != 2 {
                    return Err(fail("dot_arity".into()));
                }
                if symbol == S.k {
                    let Some(n) = vector_index(view, S.k) else {
                        return Err(fail("k".into()));
                    };
                    if !inside_dot || n < 1 || n as usize > loops {
                        return Err(fail(format!("k({n})")));
                    }
                    return Ok(false);
                }
                if symbol == S.p {
                    let Some(n) = vector_index(view, S.p) else {
                        return Err(fail("p".into()));
                    };
                    if role == Role::Block || !inside_dot || n as usize != x0 {
                        return Err(fail(format!("p({n})")));
                    }
                    return Ok(false);
                }
                Ok(true)
            }
            _ => Ok(true),
        }
    };
    walk(view, false, &mut check)
}

// ---------------------------------------------------------------------------
// Denominator orders (DESIGN 2b)
// ---------------------------------------------------------------------------

/// Whether a factor is a `den` or a power of one.
fn is_denominator(factor: AtomView<'_>) -> bool {
    match factor {
        AtomView::Fun(function) => function.get_symbol() == GS.den,
        AtomView::Pow(power) => matches!(
            power.get_base_exp().0,
            AtomView::Fun(function) if function.get_symbol() == GS.den
        ),
        _ => false,
    }
}

/// The Taylor orders `D_0 .. D_deg` of the scaled and compensated
/// denominator product, each split into its denominator monomials `μ`:
/// `μ -> [d_{0,μ} .. d_{deg,μ}]` (coefficients unmapped), and the summand
/// count of every order.
#[allow(clippy::too_many_arguments)]
pub(super) fn denominator_orders(
    denominators: &[Atom],
    root_edges: &[EdgeIndex],
    degree: usize,
    reference_mass_squared: Symbol,
    auxiliary_mass_squared: &Atom,
) -> Result<std::result::Result<(BTreeMap<Atom, Vec<Atom>>, Vec<usize>), Decline>> {
    denominator_orders_with(
        denominators,
        root_edges,
        degree,
        reference_mass_squared,
        auxiliary_mass_squared,
        true,
    )
}

/// [`denominator_orders`], with (`compensate`) or without the mass
/// compensation (the reference mass is then the auxiliary one, as at `λ=0`).
pub(super) fn denominator_orders_with(
    denominators: &[Atom],
    root_edges: &[EdgeIndex],
    degree: usize,
    reference_mass_squared: Symbol,
    auxiliary_mass_squared: &Atom,
    compensate: bool,
) -> Result<std::result::Result<(BTreeMap<Atom, Vec<Atom>>, Vec<usize>), Decline>> {
    let parameter = symbol!("rstar_outer_taylor"; Scalar);
    let taylor = TaylorExpansion::new(parameter, degree);
    let product = Atom::mul_many(denominators.iter());
    let scaled = RStarIntegrand::scale_edges(&product, root_edges, parameter);
    let compensated = if compensate {
        taylor.compensate_mass(&scaled, reference_mass_squared, auxiliary_mass_squared)
    } else {
        // At fixed reference mass (`GAMMALOOP_RSTAR_TOP_MASS_TERMS=drop`).
        scaled
            .replace(reference_mass_squared)
            .with(auxiliary_mass_squared.clone())
    };
    let orders = match taylor.orders(&compensated) {
        Ok(orders) => orders,
        Err(error) => return Ok(Err(Decline::new("den", format!("series: {error}")))),
    };
    let mut by_monomial = BTreeMap::<Atom, Vec<Vec<Atom>>>::new();
    let mut sizes = Vec::with_capacity(orders.len());
    for (j, order) in orders.iter().enumerate() {
        let expanded = order.expand();
        let summands = match expanded.as_view() {
            AtomView::Add(sum) => sum.iter().collect::<Vec<_>>(),
            AtomView::Num(_) if expanded.is_zero() => Vec::new(),
            single => vec![single],
        };
        sizes.push(summands.len());
        for summand in summands {
            let factors = match summand {
                AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
                single => vec![single],
            };
            let (dens, others): (Vec<_>, Vec<_>) = factors
                .into_iter()
                .partition(|factor| is_denominator(*factor));
            if dens.is_empty() {
                return Ok(Err(Decline::new(
                    "den",
                    "a denominator order term without denominators",
                )));
            }
            if others.iter().any(|factor| factor.contains_symbol(GS.den)) {
                return Ok(Err(Decline::new(
                    "den",
                    "a denominator inside a coefficient",
                )));
            }
            let monomial = Atom::mul_many(dens.iter().copied());
            let coefficient = Atom::mul_many(others.iter().copied());
            by_monomial
                .entry(monomial)
                .or_insert_with(|| vec![Vec::new(); degree + 1])[j]
                .push(coefficient);
        }
    }
    let by_monomial = by_monomial
        .into_iter()
        .map(|(monomial, orders)| {
            (
                monomial,
                orders
                    .into_iter()
                    .map(|coefficients| Atom::add_many(coefficients.iter()))
                    .collect(),
            )
        })
        .collect();
    Ok(Ok((by_monomial, sizes)))
}

// ---------------------------------------------------------------------------
// The build
// ---------------------------------------------------------------------------

/// Everything the build needs from the evaluation, besides the remainder.
pub(super) struct BoundaryRequest<'b> {
    pub graph_name: &'b str,
    pub marked_edge: EdgeIndex,
    /// The integration id this remainder receives (`integration_id + 1`).
    pub integration_id: usize,
    pub root_edges: &'b [EdgeIndex],
    pub degree: usize,
    pub reference_mass_squared: Symbol,
    pub auxiliary_mass_squared: &'b Atom,
    pub vakint: &'b Vakint,
    /// `VakintSettings` with the remainder's Laurent depth.
    pub settings: &'b VakintSettings,
    /// The vakint `additional_normalization` (one power per vacuum loop).
    pub additional_normalization: &'b str,
    pub epsilon: Symbol,
    pub closed_root: &'b ClosedRootNormalization,
    /// `GAMMALOOP_RSTAR_TOP_MASS_TERMS` of this (top-level) root.
    pub mass_terms: super::TopMassTerms,
}

/// One kept or zero pattern, for check mode and records (the key, the zero
/// patterns and the grades are for the census record of DESIGN 6b).
#[allow(dead_code)]
#[derive(Clone, Debug)]
pub(super) struct PatternMeta {
    /// The position of the pattern among all patterns (kept and zero), in key
    /// order: today's closure order.
    pub order: usize,
    /// The cut topology `topo(prop(..))` with this pattern's powers.
    pub topology: Atom,
    /// The grouping key (today's topology key, `topo(prop(..))`).
    pub key: Atom,
    /// `I_π`: today's `VakintTerm.integral` (rebuilt, routed, vakint symbols).
    pub integral: Atom,
    pub closure_power: i64,
    pub unprojected_nonzero: bool,
}

/// A built boundary with its metadata.
#[allow(dead_code)]
pub(super) struct BuiltBoundary {
    pub boundary: RemainderBoundary,
    /// Per kept pattern, in the order of `boundary.patterns`.
    pub patterns: Vec<PatternMeta>,
    pub zero_patterns: Vec<PatternMeta>,
    /// Per block index, its grade.
    pub grades: Vec<Grade>,
    /// The vacuum loop count `L`.
    pub loops: usize,
    /// The R* loop basis edge ids `b_j`, in order (`k(j+1) = Q(b_j)`).
    pub basis: Vec<usize>,
    /// `c̃` with the vakint normalisation of every vacuum loop (folded into
    /// every coefficient; informational in the record).
    pub normalisation: Atom,
    /// `p(x0)` in vakint symbols (the `k·p` structures of the grades).
    pub external: Atom,
    /// The outer Taylor degree.
    pub degree: usize,
    /// The integration id of this remainder.
    pub integration_id: usize,
}

/// The outcome of a build.
pub(super) enum Built {
    Boundary(BuiltBoundary, BoundaryStats),
    /// Every pattern is scaleless, or the remainder vanishes: the value is 0.
    Zero(BoundaryStats),
    Declined(Decline, BoundaryStats),
}

/// One pattern under construction.
struct Pattern {
    key: Atom,
    /// `(μ, rest)` of the raw patterns merged into this one.
    raw: Vec<(Atom, Atom)>,
    frame: TermFrame,
    cut: VacuumCut,
}

/// Build the boundary of one admitted remainder, cheap stages first (DESIGN
/// 6a): den, frames, cut, routing, contract, collect, project, assemble,
/// verify. Hard errors (a failing contraction, a failing projection) are
/// errors; anything the boundary cannot represent is a decline.
pub(super) fn build_boundary(
    remainder: &RStarIntegrand,
    request: &BoundaryRequest<'_>,
) -> Result<Built> {
    let loops = remainder.region.loops;
    let mut stats = BoundaryStats {
        loops,
        degree: request.degree,
        ..Default::default()
    };
    let mut clock = StageClock::new(request.graph_name, request.integration_id);
    macro_rules! declined {
        ($decline:expr) => {
            return Ok(Built::Declined($decline, stats))
        };
    }

    // den: the denominator product alone.
    if remainder.routed_is_zero() {
        clock.lap(&mut stats, "den");
        return Ok(Built::Zero(stats));
    }
    let (numerator, denominators) = match remainder.split_routed() {
        Ok(split) => split,
        Err(reason) => declined!(Decline::new("den", reason)),
    };
    let (by_monomial, sizes) = match denominator_orders_with(
        &denominators,
        request.root_edges,
        request.degree,
        request.reference_mass_squared,
        request.auxiliary_mass_squared,
        request.mass_terms.compensates(),
    )? {
        Ok(orders) => orders,
        Err(decline) => declined!(decline),
    };
    // GAMMALOOP_TOP_MASS_SPLIT: the den orders of the other mode as well (the
    // split of the projected coefficients by origin, `TOP_MASS_BOUNDARY`).
    let other_orders = if super::top_mass::split_enabled() {
        match denominator_orders_with(
            &denominators,
            request.root_edges,
            request.degree,
            request.reference_mass_squared,
            request.auxiliary_mass_squared,
            !request.mass_terms.compensates(),
        )? {
            Ok((orders, _)) => Some(orders),
            Err(decline) => declined!(decline),
        }
    } else {
        None
    };
    stats.d_orders = sizes;
    stats.raw_patterns = by_monomial.len();
    clock.lap(&mut stats, "den");

    // frames: today's topology key per monomial, merged; the vacuum frame.
    let converter = remainder.converter_graph();
    let rewrites = TermRewrites::new(
        &converter.underlying,
        &remainder.region.remaining,
        &converter.empty_subgraph::<SuBitGraph>(),
        false,
    );
    let mut grouped = BTreeMap::<Atom, Vec<(Atom, Atom)>>::new();
    for monomial in by_monomial.keys() {
        let Some((key, rest)) = rewrites.denominator_image(monomial) else {
            declined!(Decline::new("frames", "no single topology"));
        };
        grouped
            .entry(key)
            .or_default()
            .push((monomial.clone(), rest));
    }
    let mut frames = Vec::with_capacity(grouped.len());
    for (key, raw) in grouped {
        let frame = match term_frame(&key) {
            Ok(frame) => frame,
            Err(error) => declined!(Decline::new("frames", format!("momentum solve: {error:#}"))),
        };
        frames.push((key, raw, frame));
    }
    clock.lap(&mut stats, "frames");

    // cut: once per pattern, at its massive line.
    let x0 = usize::from(request.closed_root.external);
    let probe = function!(S.p, x0 as i64);
    let mut patterns = Vec::with_capacity(frames.len());
    for (key, raw, frame) in frames {
        let cut = match cut_frame(&frame.integral, &probe) {
            Ok(Some(cut)) => cut,
            Ok(None) => declined!(Decline::new("cut", "no massive line")),
            Err(error) => declined!(Decline::new("cut", format!("{error:#}"))),
        };
        patterns.push(Pattern {
            key,
            raw,
            frame,
            cut,
        });
    }
    let scaleless = patterns
        .iter()
        .filter(|pattern| pattern.cut.scaleless)
        .count();
    if scaleless == patterns.len() {
        clock.lap(&mut stats, "cut");
        return Ok(Built::Zero(stats));
    }
    if scaleless != 0 {
        declined!(Decline::new("cut", "scaleless patterns beside others"));
    }
    if patterns
        .iter()
        .any(|pattern| !matches!(pattern.cut.integrand.as_view(), AtomView::Fun(_)))
    {
        declined!(Decline::new("cut", "a non-topo cut"));
    }
    clock.lap(&mut stats, "cut");

    // routing: the composite vacuum-label routing, the same for every pattern.
    let basis = remainder
        .loop_basis_edges()
        .iter()
        .map(|edge| usize::from(*edge))
        .collect::<Vec<_>>();
    if basis.len() != loops {
        declined!(Decline::new(
            "routing",
            "the loop basis does not match the loop count"
        ));
    }
    let mut routing: Option<Vec<(Atom, Atom)>> = None;
    for pattern in &patterns {
        let mut images = Vec::with_capacity(loops);
        for (j, edge) in basis.iter().enumerate() {
            let Some(image) = pattern.frame.loop_image(*edge) else {
                declined!(Decline::new(
                    "routing",
                    format!("no loop image of Q({edge})")
                ));
            };
            images.push((
                Grammar::loop_vector(j),
                image.replace_multiple(&pattern.cut.routing).expand(),
            ));
        }
        match &routing {
            None => routing = Some(images),
            Some(first) if *first == images => {}
            Some(_) => declined!(Decline::new("routing", "routing_differs")),
        }
    }
    let routing = routing.expect("at least one pattern");
    let first = &patterns[0].cut.closure;
    if patterns.iter().any(|pattern| {
        pattern.cut.closure.momentum != first.momentum
            || pattern.cut.closure.mass_squared != first.mass_squared
    }) {
        declined!(Decline::new(
            "routing",
            "the closing momentum or mass differs"
        ));
    }
    let closing = MassiveClosure {
        momentum: first.momentum.clone(),
        mass_squared: first.mass_squared.clone(),
        power: first.power,
    };
    let head = match closing.momentum.as_view() {
        AtomView::Fun(function) if function.get_symbol() == S.p && function.get_nargs() == 1 => {
            usize::try_from(function.get(0)).ok()
        }
        _ => None,
    };
    let Some(head) = head else {
        declined!(Decline::new("routing", "the closing momentum is not p(h)"));
    };
    stats.routing_terms = routing
        .iter()
        .map(|(_, image)| match image.as_view() {
            AtomView::Add(sum) => sum.get_nargs(),
            _ => 1,
        })
        .sum();
    clock.lap(&mut stats, "routing");

    // contract: P = C(σ(N)), without the denominators.
    let (p, certified) = remainder
        .numerator_contracted(numerator, &request.closed_root.external_routing)
        .wrap_err("while contracting the numerator of an R* remainder boundary")?;
    stats.certified = certified;
    stats.contracted = true;
    stats.p_summands = vakint::census::summands(p.as_view());
    stats.p_bytes = p.as_view().get_byte_size();
    remainder.assert_no_combinations()?;
    clock.lap(&mut stats, "contract");
    stats.rss_gb[0] = rss_gb();

    // collect: the graded split of P into blocks.
    let mut grammar = Grammar::new(&basis, x0, request.degree, request.reference_mass_squared);
    let blocks = match split_into_blocks(p, &mut grammar, &mut stats)? {
        Ok(blocks) => blocks,
        Err(decline) => declined!(decline),
    };
    clock.lap(&mut stats, "collect");
    stats.rss_gb[1] = rss_gb();
    if blocks.is_empty() {
        return Ok(Built::Zero(stats));
    }

    // project: R_{π,g} = c̃ · P[ρ_{π,g}], expanded.
    let normalization = match grammar.map_scalar(
        super::linearize::linearize(
            request
                .closed_root
                .factor
                .replace_multiple(&request.closed_root.external_routing),
        )
        .as_view(),
    )? {
        Ok(factor) => factor,
        Err(decline) => declined!(decline),
    };
    let additional = symbolica::parse!(request.additional_normalization);
    let mut tensor_settings = request.settings.clone();
    tensor_settings.tensor_reduction_method = TensorReductionMethod::FeynKit;
    tensor_settings.use_dot_product_notation = true;
    let mut images = HashMap::<(Atom, Atom), Atom>::new();
    let auxiliary = request.auxiliary_mass_squared;
    let mut mapped_orders = HashMap::<Atom, Vec<Atom>>::new();
    for (monomial, orders) in &by_monomial {
        let mut mapped = Vec::with_capacity(orders.len());
        for order in orders {
            let routed = super::linearize::linearize(
                order.replace_multiple(&request.closed_root.external_routing),
            );
            match grammar.map_scalar(routed.as_view())? {
                Ok(order) => mapped.push(order),
                Err(decline) => declined!(decline),
            }
        }
        mapped_orders.insert(monomial.clone(), mapped);
    }
    // GAMMALOOP_TOP_MASS_SPLIT: the other mode's den orders, mapped alike.
    let mut mapped_other = HashMap::<Atom, Vec<Atom>>::new();
    for (monomial, orders) in other_orders.iter().flatten() {
        let mut mapped = Vec::with_capacity(orders.len());
        for order in orders {
            let routed = order.replace_multiple(&request.closed_root.external_routing);
            match grammar.map_scalar(routed.as_view())? {
                Ok(order) => mapped.push(order),
                Err(decline) => declined!(decline),
            }
        }
        mapped_other.insert(monomial.clone(), mapped);
    }
    // Per origin: projected terms, their summands and projection seconds.
    let mut origin_stats = BTreeMap::<&'static str, (usize, usize, f64)>::new();
    let header_normalisation = &normalization * additional.pow(patterns[0].frame.nloops as i64);
    let block_grades = blocks.keys().cloned().collect::<Vec<_>>();
    let mut kept = Vec::new();
    let mut zero_patterns = Vec::new();
    let mass_terms = request.mass_terms;
    for (order, pattern) in patterns.iter().enumerate() {
        // d_{j,π} = Σ_{μ∈π} rest_μ · d_{j,μ}.
        let mut orders = vec![Vec::<Atom>::new(); request.degree + 1];
        let mut other = vec![Vec::<Atom>::new(); request.degree + 1];
        for (monomial, rest) in &pattern.raw {
            let rest = match grammar.map_scalar(rest.as_view())? {
                Ok(rest) => rest,
                Err(decline) => declined!(decline),
            };
            for (j, order) in mapped_orders[monomial].iter().enumerate() {
                if !order.is_zero() {
                    orders[j].push(&rest * order);
                }
            }
            for (j, order) in mapped_other.get(monomial).into_iter().flatten().enumerate() {
                if !order.is_zero() {
                    other[j].push(&rest * order);
                }
            }
        }
        let orders = orders
            .iter()
            .map(|terms| Atom::add_many(terms.iter()))
            .collect::<Vec<_>>();
        let other = other
            .iter()
            .map(|terms| Atom::add_many(terms.iter()))
            .collect::<Vec<_>>();
        let normalization = &normalization * additional.pow(pattern.frame.nloops as i64);
        let mut terms = Vec::new();
        let mut unprojected_nonzero = false;
        for (index, grade) in block_grades.iter().enumerate() {
            // ρ = t_g (M²)^b Σ_c C(b,c) (-1)^c Σ_{j ≤ deg-a-2c, j ≡ a (2)} d_j.
            let a = grade.degree();
            let b = grade.mass_power;
            if other_orders.is_some() {
                // The split by origin (GAMMALOOP_TOP_MASS_SPLIT), from the
                // den orders with (`comp`) and without (`fixed`) compensation.
                let (comp, fixed) = if mass_terms.compensates() {
                    (&orders, &other)
                } else {
                    (&other, &orders)
                };
                let parity = |j: usize| j % 2 == a % 2 && a + j <= request.degree;
                let mut origins = vec![
                    ("l0", Vec::new()),
                    ("l2_mom", Vec::new()),
                    ("comp_den", Vec::new()),
                    ("comp_ins", Vec::new()),
                ];
                for j in (0..=request.degree).filter(|j| parity(*j)) {
                    origins[if a + j == 0 { 0 } else { 1 }]
                        .1
                        .push(fixed[j].clone());
                    origins[2].1.push(&comp[j] - &fixed[j]);
                }
                for c in 1..=b {
                    let Some(limit) = request.degree.checked_sub(a + 2 * c) else {
                        break;
                    };
                    let binomial = binomial(b, c) * if c % 2 == 0 { 1 } else { -1 };
                    for j in (0..=limit).filter(|j| j % 2 == a % 2) {
                        origins[3].1.push(Atom::num(binomial) * &comp[j]);
                    }
                }
                for (origin, parts) in origins {
                    let started = Instant::now();
                    let rho = (grade.structure_atom(grammar.external())
                        * auxiliary.pow(b as i64)
                        * Atom::add_many(parts.iter()))
                    .expand();
                    let entry = origin_stats.entry(origin).or_default();
                    if !rho.is_zero() {
                        let vectors =
                            vakint::VakintTerm::identify_vectors_in_numerator(rho.as_view())
                                .map_err(|error| eyre!("{error}"))?;
                        let term = vakint::VakintTerm {
                            integral: patterns[0].frame.integral.clone(),
                            numerator: rho,
                            vectors,
                        };
                        let projected = tensor_reduce_by_factors(
                            request.vakint,
                            &tensor_settings,
                            &term,
                            &mut images,
                        )?
                        .unwrap_or_else(|| term.numerator.clone());
                        let coefficient = (&normalization * projected).expand();
                        if !coefficient.is_zero() {
                            entry.0 += 1;
                            entry.1 += vakint::census::summands(coefficient.as_view());
                        }
                    }
                    entry.2 += started.elapsed().as_secs_f64();
                }
            }
            let mut sum = Vec::new();
            for c in 0..=b {
                // GAMMALOOP_RSTAR_TOP_MASS_TERMS=drop: no compensation (c = 0).
                if c > 0 && !mass_terms.compensates() {
                    break;
                }
                let Some(limit) = request.degree.checked_sub(a + 2 * c) else {
                    break;
                };
                let binomial = binomial(b, c) * if c % 2 == 0 { 1 } else { -1 };
                for (j, order) in orders.iter().enumerate().take(limit + 1) {
                    if j % 2 == a % 2 && !order.is_zero() && mass_terms.keeps_order(a + j + 2 * c) {
                        sum.push(Atom::num(binomial) * order);
                    }
                }
            }
            let sum = Atom::add_many(sum.iter());
            if sum.is_zero() {
                continue;
            }
            unprojected_nonzero = true;
            let rho =
                (grade.structure_atom(grammar.external()) * auxiliary.pow(b as i64) * sum).expand();
            if rho.is_zero() {
                continue;
            }
            let vectors = vakint::VakintTerm::identify_vectors_in_numerator(rho.as_view())
                .map_err(|error| eyre!("{error}"))?;
            let term = vakint::VakintTerm {
                integral: patterns[0].frame.integral.clone(),
                numerator: rho,
                vectors,
            };
            let projected =
                tensor_reduce_by_factors(request.vakint, &tensor_settings, &term, &mut images)?
                    .unwrap_or_else(|| term.numerator.clone());
            let coefficient = (&normalization * projected).expand();
            if coefficient.is_zero() {
                continue;
            }
            let summands = vakint::census::summands(coefficient.as_view());
            stats.r_terms += 1;
            stats.r_summands_total += summands;
            stats.r_summands_max = stats.r_summands_max.max(summands);
            terms.push((index, coefficient));
        }
        let meta = PatternMeta {
            order,
            topology: pattern.cut.integrand.clone(),
            key: pattern.key.clone(),
            integral: pattern.frame.integral.clone(),
            closure_power: pattern.cut.closure.power,
            unprojected_nonzero,
        };
        if terms.is_empty() {
            zero_patterns.push(meta);
        } else {
            kept.push((
                meta,
                BoundaryPattern {
                    topology: pattern.cut.integrand.clone(),
                    closure_power: pattern.cut.closure.power,
                    terms,
                },
            ));
        }
        checkpoint("project")?;
    }
    stats.images = images.len();
    stats.patterns = kept.len();
    stats.zero_patterns = zero_patterns.len();
    stats.classes = kept
        .iter()
        .map(|(meta, _)| meta.closure_power)
        .collect::<BTreeSet<_>>()
        .len();
    clock.lap(&mut stats, "project");
    for (origin, (terms, summands, seconds)) in &origin_stats {
        eprintln!(
            "TOP_MASS_BOUNDARY graph={} id={} loops={loops} deg={} mode={} origin={origin} terms={terms} summands={summands} project_s={seconds:.3}",
            request.graph_name,
            request.integration_id,
            request.degree,
            mass_terms.as_str()
        );
    }
    if kept.is_empty() {
        return Ok(Built::Zero(stats));
    }

    // assemble: the blocks in grade order (every block is used by some
    // pattern or kept for the record), the patterns in key order.
    let blocks = blocks.into_values().collect::<Vec<_>>();
    stats.blocks = blocks.len();
    for block in &blocks {
        let size = (
            vakint::census::summands(block.as_view()),
            block.as_view().get_byte_size(),
        );
        stats.block_max = (stats.block_max.0.max(size.0), stats.block_max.1.max(size.1));
        stats.blocks_total = (stats.blocks_total.0 + size.0, stats.blocks_total.1 + size.1);
    }
    let label = format!(
        "graph={} marked_edge={} integration_id={} loops={loops}",
        request.graph_name,
        usize::from(request.marked_edge),
        request.integration_id
    );
    let (metas, boundary_patterns): (Vec<_>, Vec<_>) = kept.into_iter().unzip();
    let boundary = match RemainderBoundary::new(
        label,
        closing.momentum.clone(),
        closing.mass_squared.clone(),
        routing,
        blocks,
        boundary_patterns,
    ) {
        Ok(boundary) => boundary,
        Err(error) => declined!(Decline::new("assemble", format!("{error}"))),
    };
    clock.lap(&mut stats, "assemble");

    // verify: the vakint contract of every block, coefficient and image.
    for block in &boundary.blocks {
        if let Err(decline) = verify_vakint_scalar(
            block.as_view(),
            Role::Block,
            loops,
            x0,
            head,
            request.reference_mass_squared,
        ) {
            declined!(decline);
        }
    }
    for pattern in &boundary.patterns {
        for (_, coefficient) in &pattern.terms {
            if let Err(decline) = verify_vakint_scalar(
                coefficient.as_view(),
                Role::Coefficient,
                loops,
                x0,
                head,
                request.reference_mass_squared,
            ) {
                declined!(decline);
            }
        }
        // The vakint contract of the cut and the closing line.
        let mut heads = BTreeSet::new();
        pattern.topology.visitor(&mut |part| {
            if let AtomView::Fun(function) = part
                && function.get_symbol() == S.p
            {
                heads.insert(part.to_owned());
                return false;
            }
            true
        });
        if heads.len() != 1 || !heads.contains(&closing.momentum) {
            declined!(Decline::new("verify", "vakint_contract:cut_head"));
        }
        if pattern.closure_power < 1 {
            declined!(Decline::new("verify", "vakint_contract:closure_power"));
        }
    }
    if closing.mass_squared.is_zero()
        || closing.mass_squared.contains_symbol(S.k)
        || closing.mass_squared.contains_symbol(S.p)
        || closing.mass_squared.contains_symbol(request.epsilon)
    {
        declined!(Decline::new("verify", "vakint_contract:closing_mass"));
    }
    for (_, image) in &boundary.routing {
        if let Err(decline) = verify_vakint_scalar(
            image.as_view(),
            Role::Routing,
            loops,
            x0,
            head,
            request.reference_mass_squared,
        ) {
            declined!(decline);
        }
    }
    clock.lap(&mut stats, "verify");
    stats.rss_gb[2] = rss_gb();
    Ok(Built::Boundary(
        BuiltBoundary {
            boundary,
            patterns: metas,
            zero_patterns,
            grades: block_grades,
            loops,
            basis: basis.clone(),
            normalisation: header_normalisation,
            external: grammar.external().clone(),
            degree: request.degree,
            integration_id: request.integration_id,
        },
        stats,
    ))
}

fn binomial(n: usize, k: usize) -> i64 {
    (0..k).fold(1i64, |value, i| value * (n - i) as i64 / (i as i64 + 1))
}

/// The wall seconds of the phases of [`evaluate_boundary`] (`BOUNDARY_STATS`
/// `lower:`, `reduce:`, `ordinary:`, `expand:`, which split `evaluate:`): the
/// key-plan lowering, the reduction and finalisation of the key plans, the
/// ordinary path of the declined patterns, the final expansion.
#[derive(Clone, Copy, Debug, Default)]
pub(super) struct EvaluateSplit {
    pub lower: f64,
    pub reduce: f64,
    pub ordinary: f64,
    pub expand: f64,
}

impl EvaluateSplit {
    pub(super) fn push(&self, stats: &mut BoundaryStats) {
        stats.seconds.push(("lower", self.lower));
        stats.seconds.push(("reduce", self.reduce));
        stats.seconds.push(("ordinary", self.ordinary));
        stats.seconds.push(("expand", self.expand));
    }
}

/// [`RemainderBoundary::evaluate`] with the block read counts: the value
/// (expanded), `(reads, hits)` and the phase seconds.
pub(super) fn evaluate_boundary(
    boundary: &RemainderBoundary,
    vakint: &Vakint,
    settings: &VakintSettings,
    options: BoundaryOptions,
) -> Result<(Atom, usize, usize, EvaluateSplit)> {
    let mut split = EvaluateSplit::default();
    // As `RemainderBoundary::evaluate` (the settings' RustRed resource
    // budget around the whole lowering and reduction), keeping the counts.
    symbolica::utils::resources::ResourceBudget::new(settings.rustred_resources).run(|| {
        symbolica::utils::resources::checkpoint().map_err(|error| eyre!("{error}"))?;
        let clock = Instant::now();
        let (classes, declined, (reads, hits)) = boundary
            .lower_counted(vakint, settings, options)
            .map_err(|error| eyre!("{error}"))
            .wrap_err("while lowering an R* remainder boundary")?;
        split.lower = clock.elapsed().as_secs_f64();
        let clock = Instant::now();
        let mut result = Atom::Zero;
        for class in classes {
            result += class
                .evaluate(vakint, settings)
                .map_err(|error| eyre!("{error}"))
                .wrap_err("while reducing an R* remainder boundary")?;
        }
        split.reduce = clock.elapsed().as_secs_f64();
        let clock = Instant::now();
        for index in declined {
            let pattern = &boundary.patterns[index];
            result += MassiveClosure::evaluate_batch(
                vakint,
                settings,
                &[(
                    boundary.closure(pattern),
                    boundary.direct_input(pattern).ordinary(),
                )],
            )
            .map_err(|error| eyre!("{error}"))?
            .pop()
            .unwrap_or_default();
        }
        split.ordinary = clock.elapsed().as_secs_f64();
        super::integrate::log_unknown_masters("boundary", vakint, settings, &result);
        symbolica::utils::resources::checkpoint().map_err(|error| eyre!("{error}"))?;
        let clock = Instant::now();
        let result = result.expand();
        split.expand = clock.elapsed().as_secs_f64();
        Ok((result, reads, hits, split))
    })
}

// ---------------------------------------------------------------------------
// Check mode (DESIGN 7)
// ---------------------------------------------------------------------------

/// `GAMMALOOP_BOUNDARY_CHECK_PAIR_MAX_BYTES` (default 16 MiB; 0 = unlimited):
/// the largest numerator, measured before any expansion, compared by value.
pub(super) fn check_pair_max_bytes() -> usize {
    static MAX: LazyLock<usize> = LazyLock::new(|| {
        match std::env::var("GAMMALOOP_BOUNDARY_CHECK_PAIR_MAX_BYTES")
            .ok()
            .and_then(|value| value.parse::<usize>().ok())
        {
            Some(0) => usize::MAX,
            Some(max) => max,
            None => 16 << 20,
        }
    });
    *MAX
}

/// `GAMMALOOP_BOUNDARY_CHECK_DUMP=<dir>`: both sides of every failing pair.
fn check_dump_directory() -> Option<&'static std::path::Path> {
    static DIRECTORY: LazyLock<Option<std::path::PathBuf>> = LazyLock::new(|| {
        std::env::var_os("GAMMALOOP_BOUNDARY_CHECK_DUMP")
            .filter(|value| !value.is_empty())
            .map(Into::into)
    });
    DIRECTORY.as_deref()
}

thread_local! {
    /// The closed-root normalisation of the remainder being checked and the
    /// outcome of the callback comparison (inert when unset).
    static CALLBACK_CHECK: std::cell::RefCell<Option<(ClosedRootNormalization, Option<&'static str>)>> =
        const { std::cell::RefCell::new(None) };
}

/// Run `f` (today's outer integration of one remainder) comparing the root
/// callback with `factor · σ` on the real Taylor polynomial; returns the
/// outcome (`identical`, `equal`, `MISMATCH`, `skipped`, or `none` when the
/// callback did not run).
pub(super) fn with_callback_check<T>(
    normalization: &ClosedRootNormalization,
    f: impl FnOnce() -> T,
) -> (T, &'static str) {
    // Restored on unwinding too: a check left set would compare every later
    // root projection of the thread against another graph's σ and factor.
    let scope = super::ScopedSlot::enter(&CALLBACK_CHECK, (normalization.clone(), None));
    let result = f();
    let outcome = scope.exit();
    (
        result,
        outcome.and_then(|(_, outcome)| outcome).unwrap_or("none"),
    )
}

/// Whether a callback comparison is pending on this thread.
pub(super) fn callback_check_active() -> bool {
    CALLBACK_CHECK.with(|slot| slot.borrow().is_some())
}

/// The callback comparison: `projected` (the callback's output after
/// `dimensionally_regularized`) against `c · σ(expanded)` (also regularized).
pub(super) fn check_callback(expanded: &Atom, projected: &Atom) {
    CALLBACK_CHECK.with(|slot| {
        let mut slot = slot.borrow_mut();
        let Some((normalization, outcome)) = slot.as_mut() else {
            return;
        };
        let result = if expanded.as_view().get_byte_size() > check_pair_max_bytes() {
            "skipped"
        } else {
            let expected = crate::uv::approx::integrated::Integrated::dimensionally_regularized(
                &(expanded.replace_multiple(&normalization.external_routing)
                    * &normalization.factor),
            );
            if expected == *projected {
                "identical"
            } else if (&expected - projected).expand().is_zero() {
                "equal"
            } else {
                "MISMATCH"
            }
        };
        // Several calls (none expected): the worst outcome wins.
        let rank = |value: &str| match value {
            "MISMATCH" => 3,
            "skipped" => 2,
            "equal" => 1,
            _ => 0,
        };
        if outcome.is_none_or(|previous| rank(result) > rank(previous)) {
            *outcome = Some(result);
        }
    });
}

/// The comparison of a built boundary with today's captured closures.
pub(super) struct CheckReport {
    pub old_closures: usize,
    pub pairs: usize,
    /// `identical`, `equal`, `MISMATCH` or `skipped:<n>`.
    pub input: String,
    pub value: &'static str,
    /// One `BOUNDARY_CHECK_PAIR` line per failing pair.
    pub failing: Vec<String>,
    /// Pattern order -> sequence of the v1 record of the same closure.
    pub ref_closures: BTreeMap<usize, usize>,
}

impl CheckReport {
    pub(super) fn mismatch(&self, callback: &str) -> bool {
        self.input == "MISMATCH" || self.value == "MISMATCH" || callback == "MISMATCH"
    }
}

/// Compare two atoms by value: `identical` (atom equality), `equal` (the
/// expanded difference vanishes, or its `together()` below 64 KiB), else
/// `MISMATCH`.
fn compare_atoms(new: &Atom, old: &Atom) -> &'static str {
    if new == old {
        return "identical";
    }
    let difference = (new - old).expand();
    if difference.is_zero() {
        return "equal";
    }
    if difference.as_view().get_byte_size() < 64 << 10 && difference.together().expand().is_zero() {
        return "equal";
    }
    "MISMATCH"
}

/// Group the boundary's patterns and today's captured closures by the vacuum
/// integral `I_π` (today's `VakintTerm.integral`) and compare per group:
/// the closing line, the cut topology and the routed, expanded inputs (after
/// relabelling today's closing head to the boundary's when they differ).
/// `new_value`/`old_value`: the two values before `finish_tensor` (`None`:
/// not evaluated).
pub(super) fn compare_pairs(
    built: Option<&BuiltBoundary>,
    captured: &[super::integrate::CapturedClosure],
    new_value: Option<&Atom>,
    old_value: Option<&Atom>,
    max_bytes: usize,
    label: &str,
) -> CheckReport {
    #[derive(Default)]
    struct Group {
        kept: Vec<usize>,
        captured: Vec<usize>,
        /// Pattern orders of the kept patterns and of the zero patterns that
        /// today's path still writes (unprojected non-zero).
        orders: Vec<(usize, bool)>,
    }
    let mut groups = BTreeMap::<Atom, Group>::new();
    if let Some(built) = built {
        for (index, meta) in built.patterns.iter().enumerate() {
            let group = groups.entry(meta.integral.clone()).or_default();
            group.kept.push(index);
            group.orders.push((meta.order, true));
        }
        for meta in built
            .zero_patterns
            .iter()
            .filter(|meta| meta.unprojected_nonzero)
        {
            groups
                .entry(meta.integral.clone())
                .or_default()
                .orders
                .push((meta.order, false));
        }
    }
    for (index, closure) in captured.iter().enumerate() {
        groups
            .entry(closure.integral.clone())
            .or_default()
            .captured
            .push(index);
    }
    let mut failing = Vec::new();
    let mut skipped = 0;
    let mut worst = "identical";
    let fail = |message: String, failing: &mut Vec<String>| {
        failing.push(format!("BOUNDARY_CHECK_PAIR {label} {message}"));
    };
    for (pair, (integral, group)) in groups.iter().enumerate() {
        let boundary = built.map(|built| &built.boundary);
        let new_head = boundary.map(|boundary| boundary.momentum.clone());
        // Today's closures, the head relabelled to the boundary's.
        let today = group
            .captured
            .iter()
            .map(|index| {
                let closure = &captured[*index];
                let relabel = match &new_head {
                    Some(head) if *head != closure.closure.momentum => {
                        vec![Replacement::new(
                            closure.closure.momentum.to_pattern(),
                            head.clone(),
                        )]
                    }
                    _ => Vec::new(),
                };
                let input = vakint::p_integral::DirectClosureInput {
                    topology: closure.input.topology.replace_multiple(&relabel),
                    routing: closure
                        .input
                        .routing
                        .iter()
                        .map(|(key, image)| (key.clone(), image.replace_multiple(&relabel)))
                        .collect(),
                    numerator: closure.input.numerator.clone(),
                };
                (closure, input)
            })
            .collect::<Vec<_>>();
        if let (Some(boundary), Some(first)) = (boundary, group.kept.first()) {
            let pattern = &boundary.patterns[*first];
            for (closure, input) in &today {
                if closure.closure.mass_squared != boundary.mass_squared
                    || closure.closure.power != pattern.closure_power
                {
                    fail(
                        format!("pair={pair} closure differs: today {:?}", closure.closure),
                        &mut failing,
                    );
                    worst = "MISMATCH";
                }
                if input.topology != pattern.topology {
                    fail(
                        format!(
                            "pair={pair} topology differs: today {} boundary {}",
                            input.topology, pattern.topology
                        ),
                        &mut failing,
                    );
                    worst = "MISMATCH";
                }
            }
            for index in &group.kept[1..] {
                if boundary.patterns[*index].topology != pattern.topology {
                    fail(
                        format!("pair={pair} patterns of one integral with two topologies"),
                        &mut failing,
                    );
                    worst = "MISMATCH";
                }
            }
        }
        // Inputs, by value, below the size cap (measured before expansion).
        let new_bytes = boundary.map_or(0, |boundary| {
            group
                .kept
                .iter()
                .flat_map(|index| boundary.patterns[*index].terms.iter())
                .map(|(block, coefficient)| {
                    boundary.blocks[*block].as_view().get_byte_size()
                        + coefficient.as_view().get_byte_size()
                })
                .sum::<usize>()
        });
        let old_bytes = today
            .iter()
            .map(|(_, input)| input.numerator.as_view().get_byte_size())
            .sum::<usize>();
        if new_bytes > max_bytes || old_bytes > max_bytes {
            skipped += 1;
            continue;
        }
        let new = boundary.map_or(Atom::Zero, |boundary| {
            Atom::add_many(
                group
                    .kept
                    .iter()
                    .map(|index| boundary.direct_input(&boundary.patterns[*index]).ordinary())
                    .collect::<Vec<_>>()
                    .iter(),
            )
        });
        let old = Atom::add_many(
            today
                .iter()
                .map(|(_, input)| input.ordinary())
                .collect::<Vec<_>>()
                .iter(),
        );
        let outcome = compare_atoms(&new, &old);
        match outcome {
            "MISMATCH" => {
                fail(
                    format!(
                        "pair={pair} input MISMATCH integral={integral} boundary_patterns={} today_closures={}",
                        group.kept.len(),
                        group.captured.len()
                    ),
                    &mut failing,
                );
                if let Some(directory) = check_dump_directory() {
                    let stem = format!(
                        "{}/check-{}-{}-{pair}",
                        directory.display(),
                        std::process::id(),
                        label.replace([' ', '='], "_")
                    );
                    let _ =
                        std::fs::write(format!("{stem}.boundary.txt"), new.to_canonical_string());
                    let _ = std::fs::write(format!("{stem}.today.txt"), old.to_canonical_string());
                }
                worst = "MISMATCH";
            }
            "equal" if worst == "identical" => worst = "equal",
            _ => {}
        }
    }
    // `ref_closure` (DESIGN 6b): the k-th pattern of a group (kept and
    // unprojected-non-zero zero patterns, key order) with the k-th captured
    // closure; failing equal counts, the kept patterns alone.
    let mut ref_closures = BTreeMap::new();
    for group in groups.values() {
        let mut orders = group.orders.clone();
        orders.sort_unstable();
        let sequences = group
            .captured
            .iter()
            .map(|index| captured[*index].sequence)
            .collect::<Vec<_>>();
        let paired = if orders.len() == sequences.len() {
            orders.iter().map(|(order, _)| *order).collect::<Vec<_>>()
        } else {
            orders
                .iter()
                .filter(|(_, kept)| *kept)
                .map(|(order, _)| *order)
                .collect()
        };
        if paired.len() == sequences.len() {
            for (order, sequence) in paired.into_iter().zip(sequences) {
                if let Some(sequence) = sequence {
                    ref_closures.insert(order, sequence);
                }
            }
        }
    }
    let input = if worst == "MISMATCH" {
        "MISMATCH".to_owned()
    } else if skipped > 0 {
        format!("skipped:{skipped}")
    } else {
        worst.to_owned()
    };
    let value = match (new_value, old_value) {
        (Some(new), Some(old)) => {
            let outcome = compare_atoms(new, old);
            if outcome == "MISMATCH" {
                fail(
                    format!("value MISMATCH boundary={new} today={old}"),
                    &mut failing,
                );
            }
            outcome
        }
        _ => "skipped",
    };
    CheckReport {
        old_closures: captured.len(),
        pairs: groups
            .values()
            .filter(|group| !group.kept.is_empty() || !group.captured.is_empty())
            .count(),
        input,
        value,
        failing,
        ref_closures,
    }
}

/// The census record `remainder_boundary` of one built remainder (DESIGN 8).
pub(super) fn record_stats(graph: &str, id: usize, stats: &BoundaryStats) {
    if !vakint::census::enabled() {
        return;
    }
    vakint::census::record("remainder_boundary", stats_json(graph, id, stats));
}

/// The counts of one built remainder as JSON (the census record, and the
/// `counts` line of the boundary record).
pub(super) fn stats_json(graph: &str, id: usize, stats: &BoundaryStats) -> serde_json::Value {
    serde_json::json!({
        "graph": graph,
        "integration_id": id,
        "loops": stats.loops,
        "degree": stats.degree,
        "p_summands": stats.p_summands,
        "p_bytes": stats.p_bytes,
        "p_truncated": stats.p_truncated,
        "certified": stats.certified,
        "notation_retries": stats.notation_retries,
        "grades": stats.grades,
        "graded_entries": stats.graded_entries,
        "d_orders": stats.d_orders,
        "raw_patterns": stats.raw_patterns,
        "patterns": stats.patterns,
        "zero_patterns": stats.zero_patterns,
        "classes": stats.classes,
        "blocks": stats.blocks,
        "block_max": [stats.block_max.0, stats.block_max.1],
        "blocks_total": [stats.blocks_total.0, stats.blocks_total.1],
        "r": [stats.r_terms, stats.r_summands_total, stats.r_summands_max],
        "images": stats.images,
        "routing_terms": stats.routing_terms,
        "block_reads": stats.block_reads,
        "block_hits": stats.block_hits,
        "boundary_blocks": vakint_boundary_blocks(),
        "seconds": stats.seconds.iter().map(|(stage, seconds)| serde_json::json!([stage, seconds])).collect::<Vec<_>>(),
        "rss_gb": stats.rss_gb,
    })
}
