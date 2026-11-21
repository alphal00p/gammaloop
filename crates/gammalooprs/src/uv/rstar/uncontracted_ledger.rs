//! The per-diagram accounting of `GAMMALOOP_UNCONTRACTED_NUMERATOR` (U7a
//! delta: Astra's U7a verdict findings 1 and 4 and its timing note, helper
//! VH's round 2 R2-I1, R2-M2 and R2-M3; the beta project's
//! results/consolidation/U7-PLAN.md requirements 16 and 17).
//!
//! # The ledger of a diagram
//!
//! A diagram's records are counted in ONE ledger ([`DiagramLedger`]), keyed
//! by the graph name and the diagram key: the beta opens one per diagram job
//! with the diagram's integral index (the `NNNNNN` of `integral-NNNNNN.json`);
//! a root evaluation outside such a ledger opens its own
//! (`diagram=root#<n>`, a process-wide ordinal), so that two root
//! evaluations of one graph name never share or split a line (R2-M2). The
//! ledger is a thread-local stack of frames: an evaluation, the R* of its UV
//! children and their Δ-cache fills run on one thread (a plan computes its
//! children in a serial loop, a fill runs on the claiming thread), so every
//! record of the diagram lands in its frame, whichever of its root
//! evaluations (insertion surrogates included) made it. The frames are
//! restored on unwinding (a guard), so a panicking diagram never leaves a
//! frame to the next diagram of its worker. A record outside every frame of
//! its thread prints an `UNCONTRACTED_ORPHAN` line (never silently
//! uncounted).
//!
//! Every EXPECTED record (a remainder of an enabled kind at the hook) is
//! exactly one disposition:
//!
//! * `compared` (check: our value compared with today's), `taken` (on:
//!   vakint's value, or a zero Δ, replaced today's);
//! * `declined:<code>` (a designed refusal; today's path took it);
//! * `error:<code>` (not designed: the diagram fails);
//! * `unchecked` (reached the hook but never compared: check mode whose
//!   today's path failed after our attempt, or a count the ledger misses);
//! * `hit:<the fill's disposition>`: a Δ-cache hit serves a child Δ that
//!   another evaluation filled; the hit counts the fill's records with the
//!   fill's dispositions (so `expected` does not depend on the cache order),
//!   and prints one `UNCONTRACTED_CACHE_HIT` line. A hit whose fill left no
//!   note (a fill made with the switch off) is `hit:unchecked`.
//!
//! A K3 term that builds no record by design (`zero_loops`: Δ(γ) inserted,
//! no vakint call today; `order_zero`: K8) is not expected; the lines list
//! it as `no_record=<code>:<n>`.
//!
//! The close prints ONE `UNCONTRACTED_DISPOSITION` line per diagram and
//! enabled kind (also when nothing was expected, so a missing line is a
//! missing diagram) and one `UNCONTRACTED_SUMMARY` line per diagram.
//!
//! # Check mode coverage (U7a.1)
//!
//! Check mode is strict (the only behaviour: R2-I1): per diagram and kind,
//! compared + declined + the hits whose fill compared or declined must be
//! `expected`; and a diagram-kind with records but none
//! compared (here or in a hit's fill) fails, unless every code of its
//! declines is named for that graph and kind by the run's exception list,
//! `GAMMALOOP_UNCONTRACTED_EXCEPTIONS=<file>` (lines
//! `graph<TAB>kind<TAB>code<TAB>reason`, `#` comments; the reason is
//! required). Every comparison counted as done must have agreed
//! (`equal == compared`, hits included). Mode on reports the same
//! dispositions; an error outcome fails the diagram there too, and the
//! all-declined rule applies to its taken records (a setting-wide decline
//! never makes mode on fall back entirely with outcome ok). Over a run
//! ([`RunTally`]), an enabled kind with no compared (check) / taken (on)
//! record fails, unless a run-level exception row names it.
//!
//! # Resume evidence
//!
//! The close returns the dispositions as JSON (`u7a-disposition.v2`) for the
//! beta's integral record; [`resumed`] reads them back when the beta reuses
//! a saved record: matching evidence (the same switch identity, graph and
//! diagram, counts that pass the verdict again) re-prints its lines with
//! `source=checkpoint`; anything else makes the beta recompute the diagram.
//! The code that wrote the evidence is bound by the beta's checkpoint
//! identity (its pipeline fingerprint hashes every path dependency's
//! sources), so a rebuilt binary never reuses another binary's records.
//!
//! # Occurrence keys and `UNCONTRACTED_TIMING`
//!
//! Every hook occurrence of a record has a key, the same in an off and an on
//! run of one input: `occurrence=<digest>#<ordinal>`, the digest (FNV-1a
//! 64) of the graph, the diagram key, the kind, the evaluation path
//! (`root`, then `/child{<edges>}` per R* child fill), the integration id and
//! the forest-term label (`rem{<remaining edges>}ct{<child supports>}`; K3:
//! the correction, term, attempt and hard part), the ordinal counting equal
//! digests in the diagram. The record lines carry it with the record text's
//! digest (`record=`). `GAMMALOOP_UNCONTRACTED_TIMING=1` (any mode, off
//! included; nothing is computed without it in mode off) prints one
//! `UNCONTRACTED_TIMING` line per hook occurrence with the hook's wall time
//! and its route, for timing claims from uniquely paired off / on
//! occurrences only.

use std::{
    any::Any,
    cell::RefCell,
    collections::BTreeMap,
    sync::{
        Arc, LazyLock,
        atomic::{AtomicUsize, Ordering},
    },
    time::Instant,
};

use eyre::{Result, eyre};
use linnet::half_edge::subgraph::SuBitGraph;
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use symbolica::atom::Atom;

use crate::{
    graph::Graph,
    model::Model,
    uv::counterterm_cache::{CacheLookup, CountertermCache},
};

use super::{
    plan::GraphRegion,
    uncontracted::{self, Attempt, Check, Kind, Mode, Outcome, Switch},
};

// ---------------------------------------------------------------------------
// The switches of the ledger: strict (a no-op), the exceptions, the timing
// ---------------------------------------------------------------------------

/// The exception list of a check run.
pub(crate) const EXCEPTIONS_VARIABLE: &str = "GAMMALOOP_UNCONTRACTED_EXCEPTIONS";

/// The per-occurrence timing lines.
pub(crate) const TIMING_VARIABLE: &str = "GAMMALOOP_UNCONTRACTED_TIMING";

/// `GAMMALOOP_UNCONTRACTED_STRICT` (R2-I1): on and check always fail a
/// diagram on an error outcome, so `1` / `on` (and unset, empty) are
/// accepted as a no-op and `0` / `off` are refused with the reason.
pub(crate) fn parse_strict(value: Option<&str>) -> std::result::Result<(), String> {
    match value.map(str::trim) {
        None | Some("" | "1" | "on") => Ok(()),
        Some("0" | "off") => Err(format!(
            "{}={} is refused: the non-strict mode is removed (helper VH round 2 R2-I1); with {} on or check an error outcome always fails its diagram. Unset {} (1 is accepted as a no-op)",
            uncontracted::STRICT_VARIABLE,
            value.unwrap_or_default().trim(),
            uncontracted::VARIABLE,
            uncontracted::STRICT_VARIABLE
        )),
        Some(other) => Err(format!(
            "{} must be unset or 1 (a no-op; strict is the only behaviour), not `{other}`",
            uncontracted::STRICT_VARIABLE
        )),
    }
}

/// One named exception: the graph's `kind` may have no compared record
/// when every one of its records is declined by `code`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Exception {
    pub graph: String,
    pub kind: Kind,
    pub code: String,
    pub reason: String,
}

/// The exception list of a check run and the digest of its text.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct Exceptions {
    pub entries: Vec<Exception>,
    /// FNV-1a 64 of the file's bytes (the run's identity), `None` without
    /// a list.
    pub digest: Option<String>,
}

impl Exceptions {
    /// Whether `code` is named for `graph`'s `kind`.
    fn names(&self, graph: &str, kind: Kind, code: &str) -> bool {
        self.entries
            .iter()
            .any(|entry| entry.graph == graph && entry.kind == kind && entry.code == code)
    }
}

/// The code of a run-level exception row (`*<TAB><kind><TAB>none_compared
/// <TAB><reason>`): the run may end with no compared (check) or taken (on)
/// record of `kind` ([`RunTally`]; pre-review of the U7a delta, finding 6).
pub(crate) const NONE_COMPARED: &str = "none_compared";

/// The exceptions of a file's text (`graph<TAB>kind<TAB>code<TAB>reason`;
/// blank lines and `#` comments skipped; every field required). The code
/// must be a designed refusal's ([`uncontracted::is_designed_code`]): an
/// exception excuses only declined records, so a code that is an error
/// (e.g. `vakint:all_massive` since requirement 13's interim rule; the
/// single-propagator shape is now the emitter's `one_loop_tree_cut`; the
/// IR hard child route's `ir_hard_route:<head>` and the guard's
/// `ir_hard_guard:<panic|error>`) would name nothing and is refused
/// instead. The one run-level row is graph `*` with the code
/// [`NONE_COMPARED`] (and that code only with graph `*`).
pub(crate) fn parse_exceptions(text: &str) -> std::result::Result<Exceptions, String> {
    let mut entries = Vec::new();
    for (n, line) in text.lines().enumerate() {
        let line = line.trim_end_matches('\r');
        if line.trim().is_empty() || line.trim_start().starts_with('#') {
            continue;
        }
        let fields = line.split('\t').map(str::trim).collect::<Vec<_>>();
        let bad = |why: &str| format!("{EXCEPTIONS_VARIABLE} line {}: {why}: `{line}`", n + 1);
        let [graph, kind, code, reason] = fields[..] else {
            return Err(bad("expected graph<TAB>kind<TAB>code<TAB>reason"));
        };
        let kind = Kind::ALL
            .into_iter()
            .find(|candidate| candidate.as_str() == kind)
            .ok_or_else(|| bad("unknown kind"))?;
        if graph.is_empty() || code.is_empty() || reason.is_empty() {
            return Err(bad("empty graph, code or reason"));
        }
        if (graph == "*") != (code == NONE_COMPARED) {
            return Err(bad(&format!(
                "graph `*` is the run-level row and takes exactly the code {NONE_COMPARED}"
            )));
        }
        if code != NONE_COMPARED && !uncontracted::is_designed_code(code) {
            return Err(bad(
                "the code is not a designed refusal (only declined records can be excepted; an error always fails its diagram)",
            ));
        }
        entries.push(Exception {
            graph: graph.to_owned(),
            kind,
            code: code.to_owned(),
            reason: reason.to_owned(),
        });
    }
    Ok(Exceptions {
        entries,
        digest: Some(format!("{:016x}", fnv64([text.as_bytes()]))),
    })
}

/// `GAMMALOOP_UNCONTRACTED_TIMING`: unset, empty, `0`, `off` (no lines) or
/// `1`, `on`.
pub(crate) fn parse_timing(value: Option<&str>) -> std::result::Result<bool, String> {
    match value.map(str::trim) {
        None | Some("" | "0" | "off") => Ok(false),
        Some("1" | "on") => Ok(true),
        Some(other) => Err(format!("{TIMING_VARIABLE} must be 0 or 1, not `{other}`")),
    }
}

/// The routes that produce the records of a kind (review of the U7a delta,
/// zero-comparison lens, finding 2): `child` records exist only when UV
/// children take the R* route (`GAMMALOOP_UV_CHILDREN` not `massive`),
/// `ir-hard` records only on the IR hard child route (`GAMMALOOP_IR_HARD`
/// `child` or `child-check`). `ir_quotient` is the decorated IR quotient's
/// route (`GAMMALOOP_IR_QUOTIENT`; `Err`: a value that variable refuses),
/// which the records of every kind can see ([`initialise`]).
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Routes {
    pub uv_children_rstar: bool,
    pub ir_hard_child: bool,
    pub ir_quotient: std::result::Result<super::top_mass::IrQuotient, String>,
}

impl Routes {
    /// The production routes (both on, the quotient by R*).
    pub(crate) const PRODUCTION: Routes = Routes {
        uv_children_rstar: true,
        ir_hard_child: true,
        ir_quotient: Ok(super::top_mass::IrQuotient::Rstar),
    };
}

/// The value of `GAMMALOOP_IR_QUOTIENT` that selects `route`.
fn ir_quotient_name(route: super::top_mass::IrQuotient) -> &'static str {
    use super::top_mass::IrQuotient;
    match route {
        IrQuotient::Off => "off",
        IrQuotient::Massive => "massive",
        IrQuotient::Rstar => "rstar",
        IrQuotient::Auto => "auto",
        IrQuotient::Check => "check",
    }
}

/// What the switch's initialisation checks besides the switch itself
/// ([`uncontracted::validate`]): the timing switch (every mode, off
/// included: it keys an off run's occurrences), and with the switch on or
/// check `GAMMALOOP_UNCONTRACTED_STRICT`, an enabled kind whose route is
/// off (refused: the kind would expect nothing, a check that compares
/// nothing), the PINT6 dump (refused in on and check: every record would
/// decline `pint6_dump`; since the pre-review of the U7a delta, finding 7,
/// mode on too) and the exception list (on and check: mode on applies the
/// all-declined rule too; read here), and `GAMMALOOP_IR_QUOTIENT` other
/// than `rstar` (refused in on and check: helper VH's C7 settings
/// exception E2). With the switch off nothing else is
/// read, as before the switch existed. `pint6` is whether
/// `PINT6_CLOSURE_DUMP` is set, `routes` the routes in force; the values
/// are the environment's.
pub(crate) fn initialise(
    switch: Switch,
    strict: Option<&str>,
    timing: Option<&str>,
    routes: impl FnOnce() -> Routes,
    pint6: impl FnOnce() -> bool,
    exceptions: Option<&str>,
) -> std::result::Result<Exceptions, String> {
    parse_timing(timing)?;
    if switch.mode == Mode::Off {
        return Ok(Exceptions::default());
    }
    parse_strict(strict)?;
    if switch.enabled_kinds().is_empty() {
        return Err(format!(
            "{}: `{}` enables no kind this build emits (inert: {}); the switch would build no record",
            uncontracted::VARIABLE,
            switch.mode.as_str(),
            switch
                .inert_kinds()
                .iter()
                .map(|kind| kind.as_str())
                .collect::<Vec<_>>()
                .join(",")
        ));
    }
    let routes = routes();
    for (kind, active, variable) in [
        (
            Kind::Child,
            routes.uv_children_rstar,
            "GAMMALOOP_UV_CHILDREN=massive",
        ),
        (
            Kind::IrHard,
            routes.ir_hard_child,
            "GAMMALOOP_IR_HARD (not child or child-check)",
        ),
    ] {
        if switch.enabled(kind) && !active {
            return Err(format!(
                "{}={}:{} enables `{}`, whose route is off ({variable}): the kind would expect no record and its check would compare nothing; enable the route or drop the kind from the list",
                uncontracted::VARIABLE,
                switch.mode.as_str(),
                switch
                    .enabled_kinds()
                    .iter()
                    .map(|kind| kind.as_str())
                    .collect::<Vec<_>>()
                    .join(","),
                kind.as_str()
            ));
        }
    }
    // helper VH's C7, settings exception E2 (c7-normalisation.md sections
    // 7-8): only the R* quotient route. The massive route's Z̃_F(G/γ)
    // (`massive`, `off`, `check`, which uses the massive value, and `auto`,
    // massive for every quotient of at most four loops) carries non-local
    // log(m_uv_vacuum) pieces that cancel only in the complete sum
    // (`ir_hard_child`'s `equal_mod_log_mUV`), and m_uv_vacuum is not a
    // mass symbol of the records (no `insertion_mass_symbol`). Z̃ enters
    // the R̄ of every UV child with an IR correction (the IR rearrangement
    // of the child's R*), whose Δ the records of EVERY kind carry as `ct:`
    // insertions (root: the root's children; child: the child's; ir-hard:
    // the hard part's), and it multiplies the K3 values after the call. So
    // refused for any kind list, at start-up: a per-record decline would be
    // setting-wide (pre-review of the U7a delta, finding 7). The start-up
    // refusal is not the only defence: under `rstar` a quotient the R*
    // route declines (`IR_QUOTIENT ... declined=`) still falls back to the
    // massive Z̃ (`massive_declined`), so a Δ with m_uv_vacuum (polynomial
    // or in a log) is the record's ERROR `insertion_nonlocal`
    // (`uncontracted::nonlocal_form`).
    match &routes.ir_quotient {
        Ok(super::top_mass::IrQuotient::Rstar) => {}
        Ok(route) => {
            return Err(format!(
                "{}={}:{} with GAMMALOOP_IR_QUOTIENT={} is refused: the massive quotient route's Z̃_F carries non-local log(m_uv_vacuum) pieces that cancel only in the complete sum, and Z̃ enters the Δ of every UV child with an IR correction, which the records of every kind carry as insertions (C7 settings exception E2); unset GAMMALOOP_IR_QUOTIENT or set it to rstar (production)",
                uncontracted::VARIABLE,
                switch.mode.as_str(),
                switch
                    .enabled_kinds()
                    .iter()
                    .map(|kind| kind.as_str())
                    .collect::<Vec<_>>()
                    .join(","),
                ir_quotient_name(*route)
            ));
        }
        Err(error) => {
            return Err(format!("{}: {error}", uncontracted::VARIABLE));
        }
    }
    if pint6() {
        // mode on too since the pre-review of the U7a delta (finding 7): a
        // setting-wide decline would make it fall back entirely
        return Err(format!(
            "{}={} with PINT6_CLOSURE_DUMP set is refused: the PINT6 dump keeps today's records, so every record would decline (pint6_dump) and the run would {} nothing; unset one of them",
            uncontracted::VARIABLE,
            switch.mode.as_str(),
            if switch.mode == Mode::Check {
                "compare"
            } else {
                "take"
            }
        ));
    }
    match exceptions.map(str::trim).filter(|path| !path.is_empty()) {
        Some(path) => {
            let text = std::fs::read_to_string(path)
                .map_err(|error| format!("{EXCEPTIONS_VARIABLE}={path}: {error}"))?;
            parse_exceptions(&text)
        }
        None => Ok(Exceptions::default()),
    }
}

/// `GAMMALOOP_UNCONTRACTED_TIMING=1`: the per-occurrence timing lines.
pub(crate) fn timing_enabled() -> bool {
    // validated at the switch's initialisation (a bad value refuses the run
    // at start-up); a value that is not unicode is off
    static TIMING: LazyLock<bool> = LazyLock::new(|| {
        parse_timing(std::env::var(TIMING_VARIABLE).ok().as_deref())
            .unwrap_or_else(|error| panic!("{error}"))
    });
    TEST_TIMING
        .with(std::cell::Cell::get)
        .unwrap_or_else(|| *TIMING)
}

thread_local! {
    static TEST_TIMING: std::cell::Cell<Option<bool>> = const { std::cell::Cell::new(None) };
    static TEST_EXCEPTIONS: RefCell<Option<Exceptions>> = const { RefCell::new(None) };
}

/// The exception list in force (the environment's, or a test's).
fn exceptions() -> Exceptions {
    TEST_EXCEPTIONS
        .with(|slot| slot.borrow().clone())
        .unwrap_or_else(uncontracted::exceptions)
}

/// Tests: run `f` with the timing lines `on` on this thread.
#[cfg(test)]
pub(crate) fn with_timing<T>(on: bool, f: impl FnOnce() -> T) -> T {
    struct Restore(Option<bool>);
    impl Drop for Restore {
        fn drop(&mut self) {
            let previous = self.0;
            TEST_TIMING.with(|slot| slot.set(previous));
        }
    }
    let _restore = Restore(TEST_TIMING.with(|slot| slot.replace(Some(on))));
    f()
}

/// Tests: run `f` with the exception list `exceptions` on this thread.
#[cfg(test)]
pub(crate) fn with_exceptions<T>(exceptions: Exceptions, f: impl FnOnce() -> T) -> T {
    struct Restore(Option<Exceptions>);
    impl Drop for Restore {
        fn drop(&mut self) {
            let previous = self.0.take();
            TEST_EXCEPTIONS.with(|slot| *slot.borrow_mut() = previous);
        }
    }
    let _restore = Restore(TEST_EXCEPTIONS.with(|slot| slot.replace(Some(exceptions))));
    f()
}

/// The identity of the switch for a checkpoint (`None` off): the mode, the
/// enabled kinds, strictness (always, R2-I1) and the exception list's
/// digest (on and check); e.g. `check:root,child,ir-hard:strict`. The
/// beta appends it to its computation fingerprint (U7a.1 (i)); the code
/// (gammaloop, vakint) is bound there by the beta's source provenance,
/// which hashes the sources of every path dependency (pre-review of the
/// U7a delta, finding 1).
pub(crate) fn identity_of(switch: Switch, exceptions: &Exceptions) -> Option<String> {
    if switch.mode == Mode::Off {
        return None;
    }
    let kinds = switch
        .enabled_kinds()
        .iter()
        .map(|kind| kind.as_str())
        .collect::<Vec<_>>()
        .join(",");
    let mut identity = format!("{}:{kinds}:strict", switch.mode.as_str());
    if let Some(witness) = witness() {
        identity.push_str(&format!(":witness={witness}"));
    }
    if let Some(digest) = &exceptions.digest {
        identity.push_str(&format!(":exceptions={digest}"));
    }
    Some(identity)
}

/// The identity of the switch in force ([`identity_of`]).
pub(crate) fn current_identity() -> Option<String> {
    identity_of(uncontracted::switch(), &exceptions())
}

/// FNV-1a 64 of the concatenated parts (a field separator between them).
pub(crate) fn fnv64<'a>(parts: impl IntoIterator<Item = &'a [u8]>) -> u64 {
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for (n, part) in parts.into_iter().enumerate() {
        if n > 0 {
            hash ^= 0x1f;
            hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
        }
        for byte in part {
            hash ^= u64::from(*byte);
            hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    hash
}

// ---------------------------------------------------------------------------
// Counts and frames
// ---------------------------------------------------------------------------

/// The records of one kind in a frame: the dispositions (exclusive; module
/// doc) and what the attempts gave.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub(crate) struct Counts {
    /// Records that reached the hook (a hit: its fill's records).
    pub expected: usize,
    /// Records built and handed to vakint, or refused before (this frame).
    pub attempted: usize,
    /// Attempts vakint evaluated (outcome uncontracted) and zero Δs, in
    /// either mode (check mode compares them).
    pub uncontracted: usize,
    pub zero: usize,
    /// Comparisons that agreed (check mode: must be `compared`, the
    /// verdict's rule; pre-review of the U7a delta, finding 5).
    pub equal: usize,
    /// The agreeing comparisons the Δ-cache hits stand for (the fills'
    /// `equal`): must be the `hit:compared` records.
    #[serde(default)]
    pub hit_equal: usize,
    /// Disposition -> records (`compared`, `taken`, `declined:<code>`,
    /// `error:<code>`, `unchecked`, `hit:<...>`).
    pub classes: BTreeMap<String, usize>,
    /// K3 terms of the kind that build no record by design, by code
    /// (`zero_loops`: Δ(γ) inserted, no vakint call today; `order_zero`: K8):
    /// not expected and never compared, today's path takes them.
    #[serde(default)]
    pub no_record: BTreeMap<String, usize>,
    /// Remainders of this kind that reached a hook while the kind is not
    /// enabled (inert: today's path, never expected): listed on the
    /// summary line as `inert=<kind>:<n>`, so a remainder classified into an
    /// inert kind is visible (review of the U7a delta, zero-comparison lens,
    /// finding 3).
    #[serde(default)]
    pub inert: usize,
}

impl Counts {
    fn add(&mut self, class: impl Into<String>, n: usize) {
        if n > 0 {
            *self.classes.entry(class.into()).or_default() += n;
        }
    }

    fn merge(&mut self, other: &Counts) {
        self.expected += other.expected;
        self.attempted += other.attempted;
        self.uncontracted += other.uncontracted;
        self.zero += other.zero;
        self.equal += other.equal;
        self.hit_equal += other.hit_equal;
        for (class, n) in &other.classes {
            self.add(class.clone(), *n);
        }
        for (code, n) in &other.no_record {
            *self.no_record.entry(code.clone()).or_default() += n;
        }
        self.inert += other.inert;
    }

    fn class(&self, name: &str) -> usize {
        self.classes.get(name).copied().unwrap_or(0)
    }

    /// The classes `<prefix><rest>` as `rest -> n`.
    fn with_prefix(&self, prefix: &str) -> BTreeMap<String, usize> {
        self.classes
            .iter()
            .filter_map(|(class, n)| Some((class.strip_prefix(prefix)?.to_owned(), *n)))
            .collect()
    }

    fn accounted(&self) -> usize {
        self.classes.values().sum()
    }

    /// `unchecked` records and any record no disposition accounts for.
    fn unchecked(&self) -> usize {
        self.class("unchecked") + self.expected.saturating_sub(self.accounted())
    }

    fn compared(&self) -> usize {
        self.class("compared")
    }
}

/// A frame of this thread's ledger stack.
enum Scope {
    /// A diagram (closed by its [`DiagramLedger`]), with the ordinals of its
    /// occurrence digests.
    Diagram {
        graph: String,
        diagram: String,
        ordinals: BTreeMap<u64, usize>,
        /// The special scope (`surrogate`, `quotient`, `ir-hard`) this
        /// thread was in when the beta opened the ledger: a scope left set
        /// by an earlier diagram (Astra U7a finding 3), whose records would
        /// all be classified into it. The close fails the diagram.
        leaked: Option<Kind>,
    },
    /// A Δ-cache fill: its counts become the fill's note and join the frame
    /// below.
    Fill,
    /// One IR correction on the child route ([`Section`]): its K3 records
    /// join the frame below when it ends, reclassified when the correction
    /// declines.
    Section,
}

struct Frame {
    id: usize,
    scope: Scope,
    counts: BTreeMap<Kind, Counts>,
}

thread_local! {
    static FRAMES: RefCell<Vec<Frame>> = const { RefCell::new(Vec::new()) };
    /// The evaluation path (module doc), one entry per R* child fill.
    static PATH: RefCell<Vec<String>> = const { RefCell::new(Vec::new()) };
    /// The occurrences of the hooks in progress (innermost last).
    static HOOKS: RefCell<Vec<Current>> = const { RefCell::new(Vec::new()) };
}

static FRAME_IDS: AtomicUsize = AtomicUsize::new(1);

/// Whether the ledger runs: the switch on or check, or the timing lines.
pub(crate) fn active() -> bool {
    uncontracted::switch().mode != Mode::Off || timing_enabled()
}

fn push(scope: Scope) -> usize {
    let id = FRAME_IDS.fetch_add(1, Ordering::Relaxed);
    FRAMES.with(|frames| {
        frames.borrow_mut().push(Frame {
            id,
            scope,
            counts: BTreeMap::new(),
        })
    });
    id
}

/// Remove the frame `id` (and any frame left above it); never panics (it
/// runs on unwinding).
fn pop(id: usize) -> Option<Frame> {
    FRAMES
        .try_with(|frames| {
            let mut frames = frames.try_borrow_mut().ok()?;
            let at = frames.iter().rposition(|frame| frame.id == id)?;
            let mut popped = frames.split_off(at);
            Some(popped.swap_remove(0))
        })
        .ok()
        .flatten()
}

/// `f` on the innermost frame (nothing without one: a record outside any
/// diagram or fill, e.g. a test calling the emitter directly).
fn with_top(f: impl FnOnce(&mut Frame)) {
    let _ = FRAMES.try_with(|frames| {
        if let Ok(mut frames) = frames.try_borrow_mut()
            && let Some(frame) = frames.last_mut()
        {
            f(frame);
        }
    });
}

fn with_counts(kind: Kind, f: impl FnOnce(&mut Counts)) {
    let mut counted = false;
    with_top(|frame| {
        counted = true;
        f(frame.counts.entry(kind).or_default())
    });
    if !counted && uncontracted::switch().mode != Mode::Off {
        // a record outside every diagram ledger and fill of this thread (a
        // test calling the emitter directly; in a run, a thread the ledger
        // does not reach): loud, never silently uncounted
        eprintln!(
            "UNCONTRACTED_ORPHAN kind={} thread={}",
            kind.as_str(),
            std::thread::current().name().unwrap_or("?")
        );
    }
}

/// The innermost diagram frame's graph and key (`-` outside one).
fn diagram_key() -> (String, String) {
    FRAMES
        .try_with(|frames| {
            frames.try_borrow().ok().and_then(|frames| {
                frames.iter().rev().find_map(|frame| match &frame.scope {
                    Scope::Diagram { graph, diagram, .. } => Some((graph.clone(), diagram.clone())),
                    Scope::Fill | Scope::Section => None,
                })
            })
        })
        .ok()
        .flatten()
        .unwrap_or_else(|| ("-".to_owned(), "-".to_owned()))
}

fn in_diagram() -> bool {
    FRAMES.with(|frames| {
        frames
            .borrow()
            .iter()
            .any(|frame| matches!(frame.scope, Scope::Diagram { .. }))
    })
}

/// The next ordinal of `digest` in the innermost diagram frame.
fn next_ordinal(digest: u64) -> usize {
    FRAMES.with(|frames| {
        let mut frames = frames.borrow_mut();
        match frames
            .iter_mut()
            .rev()
            .find_map(|frame| match &mut frame.scope {
                Scope::Diagram { ordinals, .. } => Some(ordinals),
                Scope::Fill | Scope::Section => None,
            }) {
            Some(ordinals) => {
                let ordinal = ordinals.entry(digest).or_default();
                *ordinal += 1;
                *ordinal
            }
            None => 0,
        }
    })
}

// ---------------------------------------------------------------------------
// The records
// ---------------------------------------------------------------------------

/// A remainder of the enabled `kind` reached the hook.
pub(super) fn expect(kind: Kind) {
    with_counts(kind, |counts| counts.expected += 1);
}

/// A remainder of `kind` reached a hook while the kind is not enabled
/// (switch on or check): counted as inert, today's path takes it.
pub(super) fn inert(kind: Kind) {
    if uncontracted::switch().mode != Mode::Off {
        with_counts(kind, |counts| counts.inert += 1);
    }
}

/// The disposition of one reported attempt ([`uncontracted::report`]).
pub(super) fn record(attempt: &Attempt, check: Option<&Check>, mode: Mode) {
    with_counts(attempt.kind, |counts| {
        counts.attempted += 1;
        let class = match &attempt.outcome {
            Outcome::Value(_) => {
                counts.uncontracted += 1;
                None
            }
            Outcome::Zero => {
                counts.zero += 1;
                None
            }
            Outcome::Declined(decline) => Some(format!("declined:{}", decline.code)),
            Outcome::Error(decline) => Some(format!("error:{}", decline.code)),
        };
        let class = class.unwrap_or_else(|| match (mode, check) {
            (Mode::Check, Some(_)) => "compared".to_owned(),
            // check mode without a comparison: today's path failed after
            // our attempt
            (Mode::Check, None) => "unchecked".to_owned(),
            _ => "taken".to_owned(),
        });
        counts.equal += usize::from(check.is_some_and(|check| check.equal));
        counts.add(class, 1);
    });
}

/// Tests: one attempt of `kind` that vakint evaluated, compared (check) or
/// taken (on).
#[cfg(test)]
pub(crate) fn record_value(kind: Kind, mode: Mode) {
    with_counts(kind, |counts| {
        counts.attempted += 1;
        counts.uncontracted += 1;
        if mode == Mode::Check {
            counts.equal += 1;
            counts.add("compared", 1);
        } else {
            counts.add("taken", 1);
        }
    });
}

/// Tests: one check-mode comparison of `kind` that agreed (`equal`) or not,
/// as [`record`] counts it (a disagreeing one is otherwise a diagnostic
/// failure of the hook: the verdict must not rely on that).
#[cfg(test)]
pub(crate) fn record_compared(kind: Kind, equal: bool) {
    with_counts(kind, |counts| {
        counts.attempted += 1;
        counts.uncontracted += 1;
        counts.equal += usize::from(equal);
        counts.add("compared", 1);
    });
}

/// Tests: one check-mode attempt of `kind` whose today's path failed
/// (`unchecked`, as the hook reports it).
#[cfg(test)]
pub(crate) fn record_unchecked(kind: Kind) {
    with_counts(kind, |counts| {
        counts.attempted += 1;
        counts.uncontracted += 1;
        counts.add("unchecked", 1);
    });
}

/// vakint's unknown-master refusal of an attempt of `kind` (a hard error).
pub(super) fn record_unknown_masters(kind: Kind) {
    with_counts(kind, |counts| {
        counts.attempted += 1;
        counts.add("error:unknown_masters", 1);
    });
}

/// An error of the enabled `kind` that is not one record's outcome (check
/// mode: the IR hard child route on our K3 values declines where today's
/// pass, [`super::evaluate`]): `error:<code>`, the diagram fails.
pub(super) fn record_error(kind: Kind, code: &str) {
    with_counts(kind, |counts| counts.add(format!("error:{code}"), 1));
}

/// A K3 term of the enabled `kind` that builds no record by design (not
/// expected: listed as `no_record=`).
pub(super) fn record_no_record(kind: Kind, code: &str) {
    with_counts(kind, |counts| {
        *counts.no_record.entry(code.to_owned()).or_default() += 1
    });
}

/// What a Δ-cache fill noted: the records it evaluated, by kind, and where.
#[derive(Clone, Debug)]
pub(crate) struct FillNote {
    pub counts: BTreeMap<Kind, Counts>,
    pub graph: String,
    pub diagram: String,
    pub path: String,
}

/// A Δ-cache fill in progress (with the switch on or check).
struct FillFrame {
    id: usize,
    done: bool,
}

impl FillFrame {
    fn push() -> Self {
        Self {
            id: push(Scope::Fill),
            done: false,
        }
    }

    /// The fill's counts join the frame below and become its note.
    fn finish(mut self) -> FillNote {
        self.done = true;
        let counts = Self::merge_down(self.id);
        let (graph, diagram) = diagram_key();
        FillNote {
            counts,
            graph,
            diagram,
            path: current_path(),
        }
    }

    fn merge_down(id: usize) -> BTreeMap<Kind, Counts> {
        let counts = pop(id).map(|frame| frame.counts).unwrap_or_default();
        for (kind, counts) in &counts {
            let counts = counts.clone();
            with_counts(*kind, |below| below.merge(&counts));
        }
        counts
    }
}

impl Drop for FillFrame {
    fn drop(&mut self) {
        if !self.done {
            Self::merge_down(self.id);
        }
    }
}

/// A Δ-cache hit on `support` of `graph`: the fill's records join this
/// frame as `hit:<their disposition>` (`hit:unchecked` for one Δ whose fill
/// left no note, under every enabled cacheable kind), and one
/// `UNCONTRACTED_CACHE_HIT` line says what the hit stands for.
fn hit(note: Option<&FillNote>, graph: &Graph, support: &SuBitGraph) {
    let switch = uncontracted::switch();
    let mut stands = BTreeMap::<Kind, Counts>::new();
    match note {
        Some(note) => {
            for (kind, fill) in &note.counts {
                let counts = stands.entry(*kind).or_default();
                counts.expected += fill.expected;
                counts.hit_equal += fill.equal + fill.hit_equal;
                for (class, n) in &fill.classes {
                    let class = if class.starts_with("hit:") {
                        class.clone()
                    } else {
                        format!("hit:{class}")
                    };
                    counts.add(class, *n);
                }
            }
        }
        None => {
            for kind in switch.enabled_kinds() {
                if kind != Kind::Root {
                    let counts = stands.entry(kind).or_default();
                    counts.expected += 1;
                    counts.add("hit:unchecked", 1);
                }
            }
        }
    }
    let records = stands.values().map(|counts| counts.expected).sum::<usize>();
    for (kind, counts) in &stands {
        with_counts(*kind, |below| below.merge(counts));
    }
    if records == 0 {
        return;
    }
    let (here_graph, diagram) = diagram_key();
    let fill = stands
        .iter()
        .flat_map(|(kind, counts)| {
            counts.classes.iter().map(move |(class, n)| {
                format!(
                    "{}/{}:{n}",
                    kind.as_str(),
                    class.strip_prefix("hit:").unwrap_or(class)
                )
            })
        })
        .collect::<Vec<_>>()
        .join(",");
    print(format!(
        "UNCONTRACTED_CACHE_HIT graph={here_graph} diagram={diagram} path={} support={{{}}} child_graph={} records={records} fill={fill} fill_graph={} fill_diagram={} fill_path={}",
        current_path(),
        super::plan::join_edge_ids(&super::plan::uv_child_edge_ids(graph, support)),
        graph.name,
        note.map_or("none", |note| note.graph.as_str()),
        note.map_or("none", |note| note.diagram.as_str()),
        note.map_or("none", |note| note.path.as_str()),
    ));
}

/// `cache.counterterm_with` for an R* child Δ (`RStarPlan::rstar_child_counterterm`):
/// with the switch on or check, the fill is a ledger frame and leaves its
/// note with the entry, and a hit counts the fill's records
/// ([`hit`]); the evaluation path names the child inside the fill. Off
/// (without the timing lines): exactly `cache.counterterm_with`.
pub(super) fn through_cache(
    cache: &CountertermCache,
    graph: &Graph,
    support: &SuBitGraph,
    model: &Model,
    settings_key: String,
    compute: impl FnOnce() -> Result<Atom>,
) -> Result<Atom> {
    if uncontracted::switch().mode == Mode::Off {
        if !timing_enabled() {
            return cache.counterterm_with(graph, support, model, settings_key, compute);
        }
        return cache.counterterm_with(graph, support, model, settings_key, || {
            let _path = ChildPath::enter(graph, support);
            compute()
        });
    }
    let (value, lookup) = cache.counterterm_noted(graph, support, model, settings_key, || {
        let _path = ChildPath::enter(graph, support);
        let fill = FillFrame::push();
        let value = compute();
        let note = fill.finish();
        value.map(|value| (value, Arc::new(note) as Arc<dyn Any + Send + Sync>))
    })?;
    if let CacheLookup::Hit(note) = &lookup {
        hit(
            note.as_deref()
                .and_then(|note| note.downcast_ref::<FillNote>()),
            graph,
            support,
        );
    }
    Ok(value)
}

/// How an IR correction on the child route ended ([`Section::close`]).
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) enum SectionEnd {
    /// The correction's value is the child route's: the K3 records keep
    /// their dispositions.
    Kept,
    /// The route declined the correction by its own reason (`moments:`,
    /// `dirac_*`, `mass_power`, `depth:`, ...) after its K3 records ran;
    /// today's closed path computes the correction. Mode on with taken
    /// records: an error (the decline follows the new values; pre-review of
    /// the U7a delta, finding 4); check mode (the route had today's values):
    /// the comparisons stand.
    RouteDeclined(String),
    /// `GAMMALOOP_IR_HARD_CHILD_GUARD=1`: the scoped budget was exhausted.
    GuardBudget,
    /// `GAMMALOOP_IR_HARD_CHILD_GUARD=1`: the guard caught a panic or an
    /// ordinary error (`panic:..`, `error:..`): always an error with the
    /// switch on or check, inside a record or not (pre-review of the U7a
    /// delta, finding 3).
    GuardFailure(String),
}

/// The code head of a decline reason: up to its first `:` (`depth:3` ->
/// `depth`), whitespace-free.
fn reason_head(reason: &str) -> String {
    reason
        .split(':')
        .next()
        .unwrap_or_default()
        .chars()
        .filter(|c| !c.is_whitespace())
        .collect()
}

/// One IR correction on the child route, a frame of the diagram ledger
/// (review of the U7a delta: refusal-integrity findings 2 and 4,
/// values-and-off finding 1; pre-review findings 3 and 4). Its K3 records
/// join the frame below when it closes; when the correction declines, they
/// are reclassified first, so the lines never say `taken` for a value
/// today's closed path replaced, and an interrupted record is never
/// silently dropped:
///
/// * a route decline ([`SectionEnd::RouteDeclined`]): in mode on each
///   `taken` record becomes the ERROR `error:ir_hard_route:<reason head>`
///   and the close returns the failure (the route declined after it took
///   the new values: a failure of the new path never falls back silently);
///   check mode keeps `compared` (the route had today's values: the
///   decline is today's; the route on our values is checked separately);
/// * an exhausted guard budget ([`SectionEnd::GuardBudget`], a designed
///   decline of the whole correction): `taken` (on) and the interrupted
///   records become `declined:ir_hard_guard:budget`;
/// * a panic or an error the guard caught ([`SectionEnd::GuardFailure`]):
///   every record of the section (taken, compared, interrupted) becomes
///   `error:ir_hard_guard:<panic|error>` and the close returns the failure,
///   whether the failure happened inside a record or after them.
///
/// An interrupted record is one that reached the hook with no disposition,
/// or one counted `unchecked` (check mode: its today's path failed inside
/// the guarded correction).
///
/// Dropped without a close (an error return or an unwind), the counts join
/// the frame below unchanged.
pub(super) struct Section {
    id: usize,
    done: bool,
}

impl Section {
    /// A section with the switch on or check (`None` off).
    pub(super) fn open() -> Option<Self> {
        (uncontracted::switch().mode != Mode::Off).then(|| Self {
            id: push(Scope::Section),
            done: false,
        })
    }

    /// Close the section by `end`; `Err` (the correction's failure, already
    /// in the counts as an error) when the end is not a designed one.
    pub(super) fn close(mut self, end: SectionEnd) -> Result<()> {
        self.done = true;
        let mut counts = pop(self.id).map(|frame| frame.counts).unwrap_or_default();
        let mode = uncontracted::switch().mode;
        let mut errors = Vec::new();
        for (kind, counts) in counts.iter_mut() {
            if let Some(error) = reclassify(counts, mode, &end) {
                errors.push(format!("kind={}: {error}", kind.as_str()));
            }
        }
        for (kind, counts) in &counts {
            with_counts(*kind, |below| below.merge(counts));
        }
        match (&end, errors.is_empty()) {
            (SectionEnd::GuardFailure(reason), _) => Err(eyre!(
                "UNCONTRACTED: GAMMALOOP_IR_HARD_CHILD_GUARD caught {reason} in an IR correction with {}={}: a failure, never a decline ({})",
                uncontracted::VARIABLE,
                mode.as_str(),
                if errors.is_empty() {
                    "no record of the correction".to_owned()
                } else {
                    errors.join("; ")
                }
            )),
            (_, false) => Err(eyre!(
                "UNCONTRACTED: the IR hard child route declined after it took the uncontracted K3 values ({}): mode on would fall back to the closed path on a failure of the new path; {}",
                match &end {
                    SectionEnd::RouteDeclined(reason) => reason.as_str(),
                    _ => "?",
                },
                errors.join("; ")
            )),
            (_, true) => Ok(()),
        }
    }
}

/// The reclassification of one kind's counts of a [`Section`] by its end;
/// the error it makes, if any. The interrupted records of a guard's end are
/// those without a disposition and those `unchecked` (check mode: today's
/// path of the record failed inside the guarded correction, e.g. on the
/// exhausted budget).
fn reclassify(counts: &mut Counts, mode: Mode, end: &SectionEnd) -> Option<String> {
    if matches!(end, SectionEnd::GuardBudget | SectionEnd::GuardFailure(_)) {
        // an `unchecked` record is one without a disposition from here on
        counts.classes.remove("unchecked");
    }
    let interrupted = counts.expected.saturating_sub(counts.accounted());
    match end {
        SectionEnd::Kept => None,
        SectionEnd::RouteDeclined(reason) => {
            let taken = (mode == Mode::On)
                .then(|| counts.classes.remove("taken"))
                .flatten()?;
            let class = format!("error:ir_hard_route:{}", reason_head(reason));
            counts.add(class.clone(), taken);
            Some(format!("{taken} taken records became {class}"))
        }
        SectionEnd::GuardBudget => {
            let class = "declined:ir_hard_guard:budget";
            if mode == Mode::On
                && let Some(taken) = counts.classes.remove("taken")
            {
                counts.add(class, taken);
            }
            counts.add(class, interrupted);
            None
        }
        SectionEnd::GuardFailure(reason) => {
            let class = format!("error:ir_hard_guard:{}", reason_head(reason));
            let mut records = interrupted;
            for done in ["taken", "compared"] {
                records += counts.classes.remove(done).unwrap_or(0);
            }
            if counts.compared() == 0 {
                // the comparisons went with their records
                counts.equal = 0;
            }
            counts.add(class.clone(), records);
            (records > 0).then(|| format!("{records} records became {class}"))
        }
    }
}

impl Drop for Section {
    fn drop(&mut self) {
        if !self.done
            && let Some(frame) = pop(self.id)
        {
            for (kind, counts) in &frame.counts {
                with_counts(*kind, |below| below.merge(counts));
            }
        }
    }
}

/// The evaluation path inside an R* child fill.
struct ChildPath;

impl ChildPath {
    fn enter(graph: &Graph, support: &SuBitGraph) -> Self {
        let segment = format!(
            "{}/child{{{}}}",
            current_path(),
            super::plan::join_edge_ids(&super::plan::uv_child_edge_ids(graph, support))
        );
        PATH.with(|path| path.borrow_mut().push(segment));
        ChildPath
    }
}

impl Drop for ChildPath {
    fn drop(&mut self) {
        let _ = PATH.try_with(|path| {
            if let Ok(mut path) = path.try_borrow_mut() {
                path.pop();
            }
        });
    }
}

/// The evaluation path of this thread (`root` outside every child fill).
pub(super) fn current_path() -> String {
    PATH.try_with(|path| path.try_borrow().ok().and_then(|path| path.last().cloned()))
        .ok()
        .flatten()
        .unwrap_or_else(|| "root".to_owned())
}

// ---------------------------------------------------------------------------
// The diagram ledger
// ---------------------------------------------------------------------------

/// One diagram's ledger (module doc): opened by the beta per diagram job, or
/// by a root evaluation outside one. [`Self::close`] prints its lines and
/// applies check mode's coverage rule; dropped without a close (an error
/// return or an unwind) it prints them with `outcome=error`.
pub struct DiagramLedger {
    id: usize,
    closed: bool,
}

impl DiagramLedger {
    /// A ledger for the diagram `diagram` (the beta: its integral index) of
    /// `graph`; `None` with the switch off (and no timing lines).
    ///
    /// The beta opens it between diagrams, outside every special scope: a
    /// thread still in one (a scope an earlier diagram left set, Astra U7a
    /// finding 3) is recorded and fails the diagram at the close, whatever
    /// its records say (a classification check independent of the hooks).
    pub fn open(graph: &str, diagram: &str) -> Option<Self> {
        Self::open_with(graph, diagram, true)
    }

    fn open_with(graph: &str, diagram: &str, between_diagrams: bool) -> Option<Self> {
        active().then(|| {
            let leaked = Some(uncontracted::classify(false))
                .filter(|kind| between_diagrams && *kind != Kind::Root);
            Self {
                id: push(Scope::Diagram {
                    graph: graph.to_owned(),
                    diagram: diagram.to_owned(),
                    ordinals: BTreeMap::new(),
                    leaked,
                }),
                closed: false,
            }
        })
    }

    /// The ledger of a root evaluation of `graph` that runs outside every
    /// diagram ledger (`diagram=root#<n>`); `None` inside one.
    pub(super) fn implicit(graph: &str) -> Option<Self> {
        static ROOTS: AtomicUsize = AtomicUsize::new(1);
        if !active() || in_diagram() {
            return None;
        }
        Self::open_with(
            graph,
            &format!("root#{}", ROOTS.fetch_add(1, Ordering::Relaxed)),
            false,
        )
    }

    /// Close the ledger of a diagram whose evaluation succeeded (`ok`) or
    /// failed: one `UNCONTRACTED_DISPOSITION` line per enabled kind, one
    /// `UNCONTRACTED_SUMMARY` line; in check mode an `ok` diagram whose
    /// coverage fails is an error. The dispositions as JSON (`None` with
    /// the switch off).
    pub fn close(mut self, ok: bool) -> Result<Option<Value>> {
        self.closed = true;
        let Some(frame) = pop(self.id) else {
            return Err(eyre!(
                "UNCONTRACTED: the diagram ledger is not on this thread's stack"
            ));
        };
        finish(frame, ok)
    }
}

impl Drop for DiagramLedger {
    fn drop(&mut self) {
        if !self.closed
            && let Some(frame) = pop(self.id)
        {
            let _ = finish(frame, false);
        }
    }
}

/// The verdict on one diagram-kind.
struct Verdict {
    outcome: &'static str,
    excepted: BTreeMap<String, usize>,
    failure: Option<String>,
}

fn list(map: &BTreeMap<String, usize>) -> String {
    if map.is_empty() {
        "none".to_owned()
    } else {
        map.iter()
            .map(|(code, n)| format!("{code}:{n}"))
            .collect::<Vec<_>>()
            .join(",")
    }
}

/// The records of `counts` that are evidence of the new path in `mode`:
/// compared (check) or taken (on), here and in the Δ-cache hits' fills.
fn evidence_of(counts: &Counts, mode: Mode) -> usize {
    let class = if mode == Mode::Check {
        "compared"
    } else {
        "taken"
    };
    counts.class(class) + counts.class(&format!("hit:{class}"))
}

fn verdict(
    counts: &Counts,
    mode: Mode,
    ok: bool,
    graph: &str,
    kind: Kind,
    exceptions: &Exceptions,
) -> Verdict {
    let mut verdict = Verdict {
        outcome: "ok",
        excepted: BTreeMap::new(),
        failure: None,
    };
    let errors = counts.with_prefix("error:");
    let hits = counts.with_prefix("hit:");
    let unchecked = counts.unchecked();
    if !ok || !errors.is_empty() || hits.keys().any(|class| class.starts_with("error:")) {
        // the diagram fails with its own error
        verdict.outcome = "error";
        return verdict;
    }
    if mode == Mode::Check {
        // every comparison counted as done agreed (a disagreeing one is a
        // diagnostic failure of the hook today; the verdict does not rely on
        // it: pre-review of the U7a delta, finding 5)
        let hit_compared = hits.get("compared").copied().unwrap_or(0);
        if counts.equal != counts.compared() || counts.hit_equal != hit_compared {
            verdict.outcome = "error";
            verdict.failure = Some(format!(
                "comparisons that did not agree: compared {} equal {}, cache hits compared {hit_compared} equal {}",
                counts.compared(),
                counts.equal,
                counts.hit_equal
            ));
            return verdict;
        }
    } else if unchecked > 0 || hits.contains_key("unchecked") {
        verdict.outcome = "unchecked";
        return verdict;
    }
    let declined = counts.with_prefix("declined:");
    if mode == Mode::Check {
        let covered_hits = hits
            .iter()
            .filter(|(class, _)| *class == "compared" || class.starts_with("declined:"))
            .map(|(_, n)| *n)
            .sum::<usize>();
        let covered = counts.compared() + declined.values().sum::<usize>() + covered_hits;
        if covered != counts.expected || unchecked > 0 {
            verdict.outcome = "uncovered";
            verdict.failure = Some(format!(
                "expected {} records, accounted compared {} + declined {} + cache hits {} = {covered} (unchecked {unchecked}, cache hits by the fill's disposition {})",
                counts.expected,
                counts.compared(),
                declined.values().sum::<usize>(),
                covered_hits,
                list(&hits)
            ));
            return verdict;
        }
    }
    if counts.expected > 0 && evidence_of(counts, mode) == 0 {
        // every record declined: accepted only when every code is named
        // (check, and since the pre-review of the U7a delta, finding 7, on:
        // a setting-wide decline never makes mode on fall back entirely
        // with outcome ok)
        let mut codes = BTreeMap::<String, usize>::new();
        for (code, n) in &declined {
            *codes.entry(code.clone()).or_default() += n;
        }
        for (class, n) in &hits {
            if let Some(code) = class.strip_prefix("declined:") {
                *codes.entry(code.to_owned()).or_default() += n;
            }
        }
        if !codes.is_empty() && codes.keys().all(|code| exceptions.names(graph, kind, code)) {
            verdict.outcome = "excepted";
            verdict.excepted = codes;
        } else {
            verdict.outcome = "uncovered";
            verdict.failure = Some(format!(
                "all {} records declined, none {} (by code: {}); a named exception ({EXCEPTIONS_VARIABLE}: graph, kind, code, reason) is required",
                counts.expected,
                if mode == Mode::Check {
                    "compared"
                } else {
                    "taken"
                },
                list(&codes)
            ));
        }
    }
    verdict
}

/// The `UNCONTRACTED_DISPOSITION` line of one diagram-kind.
fn disposition_line(
    graph: &str,
    diagram: &str,
    mode: Mode,
    kind: Kind,
    counts: &Counts,
    verdict: &Verdict,
) -> String {
    let hits = counts.with_prefix("hit:");
    format!(
        "UNCONTRACTED_DISPOSITION graph={graph} diagram={diagram} mode={} kind={} expected={} attempted={} uncontracted={} zero={} compared={} equal={} declined={} no_record={} cache_hit={} cache_hit_fill={} unchecked={} error={} excepted={} outcome={}",
        mode.as_str(),
        kind.as_str(),
        counts.expected,
        counts.attempted,
        counts.uncontracted,
        counts.zero,
        counts.compared(),
        counts.equal,
        list(&counts.with_prefix("declined:")),
        list(&counts.no_record),
        hits.values().sum::<usize>(),
        list(&hits),
        counts.unchecked(),
        list(&counts.with_prefix("error:")),
        list(&verdict.excepted),
        verdict.outcome
    )
}

/// The `UNCONTRACTED_SUMMARY` line of a diagram (its kinds summed).
fn summary_line(
    graph: &str,
    diagram: &str,
    mode: Mode,
    kinds: &[Kind],
    total: &Counts,
    outcome: &str,
) -> String {
    format!(
        "UNCONTRACTED_SUMMARY graph={graph} diagram={diagram} mode={} kinds={} attempted={} uncontracted={} zero={} compared={} equal={} declined={} error={} expected={} no_record={} cache_hit={} unchecked={} outcome={outcome}",
        mode.as_str(),
        kinds
            .iter()
            .map(|kind| kind.as_str())
            .collect::<Vec<_>>()
            .join(","),
        total.attempted,
        total.uncontracted,
        total.zero,
        total.compared(),
        total.equal,
        list(&total.with_prefix("declined:")),
        list(&total.with_prefix("error:")),
        total.expected,
        list(&total.no_record),
        total.with_prefix("hit:").values().sum::<usize>(),
        total.unchecked(),
    )
}

/// The diagram outcome of its kinds' outcomes (the worst; `excepted` is ok).
fn worst<'o>(outcomes: impl IntoIterator<Item = &'o str>) -> &'static str {
    let rank = |outcome: &str| match outcome {
        "error" => 3,
        "uncovered" => 2,
        "unchecked" => 1,
        _ => 0,
    };
    match outcomes.into_iter().map(rank).max().unwrap_or(0) {
        3 => "error",
        2 => "uncovered",
        1 => "unchecked",
        _ => "ok",
    }
}

fn print(line: String) {
    eprintln!("{line}");
    #[cfg(test)]
    {
        if line.starts_with("UNCONTRACTED_SUMMARY ") {
            SUMMARIES.with(|lines| lines.borrow_mut().push(line.clone()));
        }
        LINES.with(|lines| lines.borrow_mut().push(line));
    }
}

#[cfg(test)]
thread_local! {
    static LINES: RefCell<Vec<String>> = const { RefCell::new(Vec::new()) };
    static SUMMARIES: RefCell<Vec<String>> = const { RefCell::new(Vec::new()) };
}

/// Tests: the `UNCONTRACTED_SUMMARY` lines printed on this thread since the
/// last call (kept apart from [`take_lines`]).
#[cfg(test)]
pub(crate) fn take_summaries() -> Vec<String> {
    SUMMARIES.with(|lines| std::mem::take(&mut *lines.borrow_mut()))
}

/// Tests: the ledger's lines (dispositions, summaries, resumed lines)
/// printed on this thread since the last call.
#[cfg(test)]
pub(crate) fn take_lines() -> Vec<String> {
    LINES.with(|lines| std::mem::take(&mut *lines.borrow_mut()))
}

/// The lines, the coverage rule and the evidence of a closed diagram frame.
fn finish(frame: Frame, ok: bool) -> Result<Option<Value>> {
    let switch = uncontracted::switch();
    if switch.mode == Mode::Off {
        return Ok(None);
    }
    let Scope::Diagram {
        graph,
        diagram,
        leaked,
        ..
    } = frame.scope
    else {
        return Ok(None);
    };
    let exceptions = exceptions();
    let kinds = switch.enabled_kinds();
    let mut total = Counts::default();
    let mut outcomes = Vec::new();
    let mut failures = Vec::new();
    let mut per_kind = serde_json::Map::new();
    for kind in &kinds {
        let counts = frame.counts.get(kind).cloned().unwrap_or_default();
        let verdict = verdict(&counts, switch.mode, ok, &graph, *kind, &exceptions);
        print(disposition_line(
            &graph,
            &diagram,
            switch.mode,
            *kind,
            &counts,
            &verdict,
        ));
        if let Some(failure) = &verdict.failure {
            failures.push(format!("kind={}: {failure}", kind.as_str()));
        } else if verdict.outcome == "error" {
            // an error outcome fails the diagram here too, whether or not
            // the caller propagated the record's error (review of the U7a
            // delta, zero-comparison lens, finding 4)
            failures.push(format!(
                "kind={}: error outcome ({})",
                kind.as_str(),
                list(&counts.with_prefix("error:"))
            ));
        } else if verdict.outcome == "unchecked" {
            // mode on: a record without a disposition (R2-I1: loud)
            failures.push(format!(
                "kind={}: {} records without a disposition (unchecked)",
                kind.as_str(),
                counts.unchecked()
            ));
        }
        let mut value = serde_json::to_value(&counts)?;
        value["outcome"] = json!(verdict.outcome);
        value["excepted"] = json!(verdict.excepted);
        if !verdict.excepted.is_empty() {
            // the reasons the run's list gives for the codes it accepted
            value["excepted_reasons"] = json!(
                exceptions
                    .entries
                    .iter()
                    .filter(|entry| entry.graph == graph
                        && entry.kind == *kind
                        && verdict.excepted.contains_key(&entry.code))
                    .map(|entry| format!("{}: {}", entry.code, entry.reason))
                    .collect::<Vec<_>>()
            );
        }
        per_kind.insert(kind.as_str().to_owned(), value);
        outcomes.push(verdict.outcome);
        total.merge(&counts);
    }
    if let Some(kind) = leaked {
        failures.push(format!(
            "the diagram ledger was opened inside the {} scope (a scope left set on this thread), so its records are classified into it",
            kind.as_str()
        ));
        outcomes.push("error");
    }
    let outcome = worst(outcomes);
    let inert = frame
        .counts
        .iter()
        .filter(|(kind, counts)| !switch.enabled(**kind) && counts.inert > 0)
        .map(|(kind, counts)| (kind.as_str().to_owned(), counts.inert))
        .collect::<BTreeMap<_, _>>();
    print(format!(
        "{} inert={}{}",
        summary_line(&graph, &diagram, switch.mode, &kinds, &total, outcome),
        list(&inert),
        leaked.map_or_else(String::new, |kind| format!(
            " leaked_scope={}",
            kind.as_str()
        ))
    ));
    let evidence = json!({
        "schema": EVIDENCE_SCHEMA,
        "identity": current_identity(),
        "graph": graph,
        "diagram": diagram,
        "outcome": outcome,
        "kinds": per_kind,
        "inert": inert,
    });
    if ok && !failures.is_empty() {
        return Err(eyre!(
            "UNCONTRACTED_DISPOSITION graph={graph} diagram={diagram} mode={} outcome={outcome}: the diagram fails its ledger: {}",
            switch.mode.as_str(),
            failures.join("; ")
        ));
    }
    Ok(Some(evidence))
}

/// The schema of the dispositions in a beta integral record.
/// v2 (pre-review of the U7a delta): `hit_equal` per kind, the verdict's
/// agreement rule; v1 evidence never matches.
pub(crate) const EVIDENCE_SCHEMA: &str = "u7a-disposition.v2";

/// What a resumed (reused) diagram record says about the switch in force.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Resumed {
    /// The switch is off: nothing to say.
    Off,
    /// The record carries matching dispositions (re-printed with
    /// `source=checkpoint`).
    Checked,
    /// No matching evidence (`reason`): the beta must recompute the diagram
    /// (one `UNCONTRACTED_RESUME ... outcome=recompute` line; pre-review of
    /// the U7a delta, finding 2: a reused record is never reported without
    /// matching evidence, and a run never passes on one).
    Unchecked(String),
}

/// A diagram record reused by the beta (resume or import) under the switch
/// in force: its saved dispositions (`evidence`, the record's
/// `uncontracted` field) are matching evidence when they have this
/// schema, this switch identity, this graph and diagram key, outcome ok and
/// a line per enabled kind whose saved counts pass the verdict again (the
/// coverage rule and the exception list in force, recomputed from the
/// counts, not read from the saved outcome: review of the U7a delta,
/// zero-comparison lens, finding 5); then the saved lines are printed
/// again with `source=checkpoint`.
/// Otherwise one `UNCONTRACTED_RESUME ... outcome=recompute reason=<why>`
/// line, and the caller recomputes the diagram (which prints its own
/// lines): a reused graph is never reported ok without matching evidence
/// (U7a.1 (ii)), and a check run never ends on reused records it did not
/// check (pre-review of the U7a delta, finding 2).
///
/// What matching evidence proves (pre-review finding 8): the saved counts
/// are self-consistent and pass the verdict in force; they are bound to the
/// switch identity, the graph and diagram key and, through the beta's
/// checkpoint identity (pipeline fingerprint), to the code that wrote them,
/// but not to the record's values: an edit that keeps the counts consistent
/// is not detected here.
pub fn resumed(graph: &str, diagram: &str, evidence: Option<&Value>) -> Resumed {
    let switch = uncontracted::switch();
    if switch.mode == Mode::Off {
        return Resumed::Off;
    }
    let kinds = switch.enabled_kinds();
    let identity = current_identity();
    let reason = match evidence {
        None | Some(Value::Null) => Some("no_evidence".to_owned()),
        Some(evidence) if evidence["schema"] != EVIDENCE_SCHEMA => {
            Some(format!("schema:{}", evidence["schema"]))
        }
        Some(evidence) if evidence["identity"] != json!(identity) => Some(format!(
            "identity:{}",
            evidence["identity"].as_str().unwrap_or("none")
        )),
        Some(evidence) if evidence["outcome"] != "ok" => Some(format!(
            "outcome:{}",
            evidence["outcome"].as_str().unwrap_or("none")
        )),
        Some(evidence) if evidence["graph"] != graph || evidence["diagram"] != diagram => {
            Some(format!(
                "diagram:{}/{}",
                evidence["graph"].as_str().unwrap_or("none"),
                evidence["diagram"].as_str().unwrap_or("none")
            ))
        }
        Some(evidence)
            if kinds
                .iter()
                .any(|kind| !evidence["kinds"][kind.as_str()].is_object()) =>
        {
            Some("kinds".to_owned())
        }
        Some(evidence) => {
            let exceptions = exceptions();
            kinds.iter().find_map(|kind| {
                let saved = &evidence["kinds"][kind.as_str()];
                let Ok(counts) = serde_json::from_value::<Counts>(saved.clone()) else {
                    return Some(format!("counts:{}", kind.as_str()));
                };
                let verdict = verdict(&counts, switch.mode, true, graph, *kind, &exceptions);
                (!matches!(verdict.outcome, "ok" | "excepted")
                    || saved["outcome"] != verdict.outcome)
                    .then(|| format!("verdict:{}:{}", kind.as_str(), verdict.outcome))
            })
        }
    };
    let source = "source=checkpoint";
    match (reason, evidence) {
        (None, Some(evidence)) => {
            let mut total = Counts::default();
            for kind in &kinds {
                let saved = &evidence["kinds"][kind.as_str()];
                let counts = serde_json::from_value::<Counts>(saved.clone()).unwrap_or_default();
                let verdict = Verdict {
                    outcome: if saved["outcome"] == "excepted" {
                        "excepted"
                    } else {
                        "ok"
                    },
                    excepted: serde_json::from_value(saved["excepted"].clone()).unwrap_or_default(),
                    failure: None,
                };
                print(format!(
                    "{} {source}",
                    disposition_line(graph, diagram, switch.mode, *kind, &counts, &verdict)
                ));
                total.merge(&counts);
            }
            print(format!(
                "{} {source}",
                summary_line(graph, diagram, switch.mode, &kinds, &total, "ok")
            ));
            Resumed::Checked
        }
        (reason, _) => {
            let reason = reason.unwrap_or_else(|| "no_evidence".to_owned());
            print(format!(
                "UNCONTRACTED_RESUME graph={graph} diagram={diagram} mode={} kinds={} outcome=recompute reason={reason}",
                switch.mode.as_str(),
                kinds
                    .iter()
                    .map(|kind| kind.as_str())
                    .collect::<Vec<_>>()
                    .join(",")
            ));
            Resumed::Unchecked(reason)
        }
    }
}

// ---------------------------------------------------------------------------
// The run: per-kind totals and the run-level coverage rule
// ---------------------------------------------------------------------------

/// The totals of one enabled kind over a run.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct KindTotals {
    expected: usize,
    compared: usize,
    hit_compared: usize,
    taken: usize,
    hit_taken: usize,
    declined: usize,
    errors: usize,
}

/// The dispositions of a run's diagram records, computed and reused (the
/// beta adds every record it aggregates), for the `UNCONTRACTED_RUN` line
/// (pre-review of the U7a delta, finding 6): per enabled kind the expected,
/// compared / taken, hit and declined totals; and the run-level rules:
///
/// * an enabled kind with no compared (check) or taken (on) record, direct
///   or through a Δ-cache hit, over the whole run fails the run, unless the
///   exception list has the run-level row `*<TAB><kind><TAB>none_compared`
///   (a run whose diagrams legitimately have no record of the kind, e.g.
///   a subset without IR corrections of order ω >= 1);
/// * with every kind this build emits enabled (`root,child,ir-hard`), a
///   diagram whose counterterm is nonzero with no compared / taken record in
///   any kind and no excepted kind fails the run (with a partial kind list a
///   diagram legitimately has no record of the listed kinds: the kind-level
///   rule covers it);
/// * a record without evidence (or evidence of another identity) fails it.
pub struct RunTally {
    mode: Mode,
    kinds: Vec<Kind>,
    totals: BTreeMap<Kind, KindTotals>,
    diagrams: usize,
    failures: Vec<String>,
}

impl RunTally {
    /// A tally for the switch in force (`None` off).
    pub fn new() -> Option<Self> {
        let switch = uncontracted::switch();
        (switch.mode != Mode::Off).then(|| Self {
            mode: switch.mode,
            kinds: switch.enabled_kinds(),
            totals: BTreeMap::new(),
            diagrams: 0,
            failures: Vec::new(),
        })
    }

    /// One diagram record: its `uncontracted` evidence and whether its
    /// counterterm is nonzero.
    pub fn add(&mut self, label: &str, evidence: Option<&Value>, nonzero: bool) {
        self.diagrams += 1;
        let identity = current_identity();
        let Some(evidence) = evidence.filter(|evidence| {
            evidence["schema"] == EVIDENCE_SCHEMA && evidence["identity"] == json!(identity)
        }) else {
            self.failures.push(format!(
                "{label}: a record without the run's evidence ({})",
                evidence.map_or("none".to_owned(), |evidence| format!(
                    "schema {} identity {}",
                    evidence["schema"], evidence["identity"]
                ))
            ));
            return;
        };
        let mut evidence_any = 0;
        let mut excepted = false;
        for kind in &self.kinds {
            let saved = &evidence["kinds"][kind.as_str()];
            excepted |= saved["outcome"] == "excepted";
            let counts = serde_json::from_value::<Counts>(saved.clone()).unwrap_or_default();
            let totals = self.totals.entry(*kind).or_default();
            totals.expected += counts.expected;
            totals.compared += counts.class("compared");
            totals.hit_compared += counts.class("hit:compared");
            totals.taken += counts.class("taken");
            totals.hit_taken += counts.class("hit:taken");
            totals.declined += counts.with_prefix("declined:").values().sum::<usize>()
                + counts.with_prefix("hit:declined:").values().sum::<usize>();
            totals.errors += counts.with_prefix("error:").values().sum::<usize>()
                + counts.with_prefix("hit:error:").values().sum::<usize>();
            evidence_any += evidence_of(&counts, self.mode);
        }
        let every_kind = Kind::ALL
            .into_iter()
            .filter(|kind| kind.implemented())
            .all(|kind| self.kinds.contains(&kind));
        if every_kind && nonzero && evidence_any == 0 && !excepted {
            self.failures.push(format!(
                "{label}: the counterterm is nonzero, and no record of any kind was {}",
                self.verb()
            ));
        }
    }

    fn verb(&self) -> &'static str {
        if self.mode == Mode::Check {
            "compared"
        } else {
            "taken"
        }
    }

    /// The ` kinds=<kind>:expected=..,...;...` field of the run line.
    pub fn fields(&self) -> String {
        format!(
            " records={} kinds={}",
            self.diagrams,
            self.kinds
                .iter()
                .map(|kind| {
                    let totals = self.totals.get(kind).cloned().unwrap_or_default();
                    format!(
                        "{}:expected={},compared={},hit_compared={},taken={},hit_taken={},declined={},error={}",
                        kind.as_str(),
                        totals.expected,
                        totals.compared,
                        totals.hit_compared,
                        totals.taken,
                        totals.hit_taken,
                        totals.declined,
                        totals.errors
                    )
                })
                .collect::<Vec<_>>()
                .join(";")
        )
    }

    /// The run-level failures (empty: the run passes these rules).
    pub fn failures(&self) -> Vec<String> {
        let exceptions = exceptions();
        let mut failures = self.failures.clone();
        for kind in &self.kinds {
            let totals = self.totals.get(kind).cloned().unwrap_or_default();
            let evidence = if self.mode == Mode::Check {
                totals.compared + totals.hit_compared
            } else {
                totals.taken + totals.hit_taken
            };
            if evidence == 0 && !exceptions.names("*", *kind, NONE_COMPARED) {
                failures.push(format!(
                    "kind {}: no record {} over the run ({} expected, {} declined); a run-level exception (*<TAB>{}<TAB>{NONE_COMPARED}<TAB>reason in {EXCEPTIONS_VARIABLE}) is required",
                    kind.as_str(),
                    self.verb(),
                    totals.expected,
                    totals.declined,
                    kind.as_str()
                ));
            }
        }
        failures
    }
}

// ---------------------------------------------------------------------------
// Occurrences and UNCONTRACTED_TIMING
// ---------------------------------------------------------------------------

/// The key of the hook occurrence in progress (read by the record lines).
#[derive(Clone, Debug)]
pub(super) struct Current {
    pub diagram: String,
    pub path: String,
    pub occurrence: String,
}

/// One hook occurrence: its key, and its `UNCONTRACTED_TIMING` line when
/// the timing lines are on.
pub(super) struct Hook {
    started: Instant,
    graph: String,
    kind: Kind,
    id: usize,
    term: String,
    current: Current,
    done: bool,
}

/// The forest-term label of a K1 / K2 remainder:
/// `rem{<remaining edges>}ct{<child supports, ';' separated>}`.
pub(super) fn term_label(graph: &Graph, region: &GraphRegion) -> String {
    let edges = |support: &SuBitGraph| {
        super::plan::join_edge_ids(&super::plan::uv_child_edge_ids(graph, support))
    };
    format!(
        "rem{{{}}}ct{{{}}}",
        edges(&region.remaining),
        region
            .insertions
            .iter()
            .map(|slot| edges(&slot.support))
            .collect::<Vec<_>>()
            .join(";")
    )
}

/// A hook occurrence of `kind` (integration id `id`) on `graph`; `None`
/// when neither a record (`enabled`) nor a timing line needs its key.
pub(super) fn hook(
    graph: &Graph,
    kind: Kind,
    id: usize,
    enabled: bool,
    term: impl FnOnce() -> String,
) -> Option<Hook> {
    if !enabled && !timing_enabled() {
        return None;
    }
    let started = Instant::now();
    let term = term();
    let path = current_path();
    let (_, diagram) = diagram_key();
    let id_text = id.to_string();
    let digest = fnv64([
        graph.name.as_bytes(),
        diagram.as_bytes(),
        kind.as_str().as_bytes(),
        path.as_bytes(),
        id_text.as_bytes(),
        term.as_bytes(),
    ]);
    let current = Current {
        diagram,
        path,
        occurrence: format!("{digest:016x}#{}", next_ordinal(digest)),
    };
    HOOKS.with(|hooks| hooks.borrow_mut().push(current.clone()));
    Some(Hook {
        started,
        graph: graph.name.clone(),
        kind,
        id,
        term,
        current,
        done: false,
    })
}

impl Hook {
    /// The occurrence ends by `route` (`today`, `uncontracted`, `zero`,
    /// `check`, `declined:<code>`, `error:<code>`, `failed`): its timing
    /// line, when on.
    pub(super) fn finish(mut self, route: &str) {
        self.done = true;
        Self::pop();
        if timing_enabled() {
            let line = format!(
                "UNCONTRACTED_TIMING graph={} diagram={} kind={} id={} path={} term={} occurrence={} mode={} route={route} seconds={:.6}",
                self.graph,
                self.current.diagram,
                self.kind.as_str(),
                self.id,
                self.current.path,
                self.term,
                self.current.occurrence,
                uncontracted::switch().mode.as_str(),
                self.started.elapsed().as_secs_f64()
            );
            eprintln!("{line}");
            #[cfg(test)]
            TIMING_LINES.with(|lines| lines.borrow_mut().push(line));
        }
    }

    fn pop() {
        let _ = HOOKS.try_with(|hooks| {
            if let Ok(mut hooks) = hooks.try_borrow_mut() {
                hooks.pop();
            }
        });
    }
}

impl Drop for Hook {
    fn drop(&mut self) {
        if !self.done {
            Self::pop();
        }
    }
}

#[cfg(test)]
thread_local! {
    static TIMING_LINES: RefCell<Vec<String>> = const { RefCell::new(Vec::new()) };
}

/// Tests: the `UNCONTRACTED_TIMING` lines printed on this thread since the
/// last call.
#[cfg(test)]
pub(crate) fn take_timing_lines() -> Vec<String> {
    TIMING_LINES.with(|lines| std::mem::take(&mut *lines.borrow_mut()))
}

/// The key of the hook occurrence in progress on this thread.
pub(super) fn current() -> Option<Current> {
    HOOKS
        .try_with(|hooks| {
            hooks
                .try_borrow()
                .ok()
                .and_then(|hooks| hooks.last().cloned())
        })
        .ok()
        .flatten()
}

/// The ` diagram= path= occurrence= record=` fields of a record line.
pub(super) fn occurrence_fields(text: Option<&str>) -> String {
    let current = current();
    format!(
        " diagram={} path={} occurrence={} record={}",
        current
            .as_ref()
            .map_or("-", |current| current.diagram.as_str()),
        current
            .as_ref()
            .map_or("-", |current| current.path.as_str()),
        current
            .as_ref()
            .map_or("-", |current| current.occurrence.as_str()),
        text.map_or_else(
            || "-".to_owned(),
            |text| format!("{:016x}", fnv64([text.as_bytes()]))
        )
    )
}

// ---------------------------------------------------------------------------
// Tests: an injected outcome at the hook
// ---------------------------------------------------------------------------

#[cfg(test)]
thread_local! {
    static FAULT: RefCell<Option<(Kind, String, bool)>> = const { RefCell::new(None) };
}

/// Tests: run `f` with every attempt of `kind` refused by `code` at the
/// hook (`designed`: a designed decline, else an error outcome), as if the
/// emitter or vakint's entry had refused it.
#[cfg(test)]
pub(crate) fn with_fault<T>(kind: Kind, code: &str, designed: bool, f: impl FnOnce() -> T) -> T {
    struct Restore(Option<(Kind, String, bool)>);
    impl Drop for Restore {
        fn drop(&mut self) {
            let previous = self.0.take();
            FAULT.with(|slot| *slot.borrow_mut() = previous);
        }
    }
    let _restore =
        Restore(FAULT.with(|slot| slot.replace(Some((kind, code.to_owned(), designed)))));
    f()
}

/// Tests: the attempt with the injected outcome of [`with_fault`], if any.
#[cfg(test)]
pub(super) fn inject(mut attempt: Attempt) -> Attempt {
    FAULT.with(|slot| {
        if let Some((kind, code, designed)) = slot.borrow().as_ref()
            && *kind == attempt.kind
        {
            let decline = uncontracted::Decline {
                code: code.clone(),
                detail: "injected by the test".to_owned(),
                designed: *designed,
            };
            attempt.text = None;
            attempt.outcome = if *designed {
                Outcome::Declined(decline)
            } else {
                Outcome::Error(decline)
            };
        }
    });
    attempt
}

// ---------------------------------------------------------------------------
// Tests: a changed K3 value, a panic after the records (pre-review of the
// U7a delta, findings 3 and 4)
// ---------------------------------------------------------------------------

#[cfg(test)]
thread_local! {
    static K3_MUTATION: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    static OUTER_PANIC: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// The witness variable of a build with the `uncontracted-witness` feature.
pub(crate) const WITNESS_VARIABLE: &str = "GAMMALOOP_UNCONTRACTED_WITNESS";

/// A witness build's `GAMMALOOP_UNCONTRACTED_WITNESS=k3_beyond_order`: the
/// change of [`with_k3_mutation`] on every K3 value of the run (never in a
/// build without the feature; part of the switch identity).
pub(crate) fn witness() -> Option<&'static str> {
    #[cfg(feature = "uncontracted-witness")]
    {
        static WITNESS: LazyLock<bool> = LazyLock::new(|| {
            std::env::var(WITNESS_VARIABLE).is_ok_and(|value| value.trim() == "k3_beyond_order")
        });
        if *WITNESS {
            return Some("k3_beyond_order");
        }
    }
    None
}

/// Whether the K3 values are changed ([`mutate_k3`]).
#[cfg(any(test, feature = "uncontracted-witness"))]
fn k3_mutation_on() -> bool {
    #[cfg(test)]
    if K3_MUTATION.with(std::cell::Cell::get) {
        return true;
    }
    witness().is_some()
}

/// Tests: run `f` with every K3 value vakint returns changed ABOVE its valid
/// order: `+ ε^(order+1) · u7a_mutation(Q(999))`, a term whose hard momentum
/// is not projected, so the IR hard child route's moment split refuses the
/// value (`moments:`) while the comparison through `ε^order` agrees.
#[cfg(test)]
pub(crate) fn with_k3_mutation<T>(f: impl FnOnce() -> T) -> T {
    struct Restore(bool);
    impl Drop for Restore {
        fn drop(&mut self) {
            let previous = self.0;
            K3_MUTATION.with(|slot| slot.set(previous));
        }
    }
    let _restore = Restore(K3_MUTATION.with(|slot| slot.replace(true)));
    f()
}

/// Tests and witness builds: the attempt with the change of
/// [`with_k3_mutation`] (a K3 value only).
#[cfg(any(test, feature = "uncontracted-witness"))]
pub(super) fn mutate_k3(mut attempt: Attempt) -> Attempt {
    use spenso::structure::representation::{Minkowski, RepName};
    use symbolica::{atom::AtomCore, function, symbol};
    if k3_mutation_on()
        && attempt.kind == Kind::IrHard
        && let Outcome::Value(value) = &attempt.outcome
    {
        let vector = function!(
            crate::utils::GS.emr_mom,
            999,
            Minkowski {}.new_rep(crate::utils::GS.dim).to_symbolic([])
        );
        let change = Atom::var(crate::utils::GS.dim_epsilon).pow(Atom::num(attempt.order + 1))
            * function!(symbol!("u7a_mutation"), vector);
        attempt.outcome = Outcome::Value(value + &change);
    }
    attempt
}

/// Tests: run `f` with a panic at the start of the IR hard child route's
/// outer step (after every K3 record of the correction ran).
#[cfg(test)]
pub(crate) fn with_outer_panic<T>(f: impl FnOnce() -> T) -> T {
    struct Restore(bool);
    impl Drop for Restore {
        fn drop(&mut self) {
            let previous = self.0;
            OUTER_PANIC.with(|slot| slot.set(previous));
        }
    }
    let _restore = Restore(OUTER_PANIC.with(|slot| slot.replace(true)));
    f()
}

/// Tests: the panic of [`with_outer_panic`].
#[cfg(test)]
pub(super) fn outer_panic() {
    if OUTER_PANIC.with(std::cell::Cell::get) {
        panic!("injected by the test: a panic in the IR hard child route's outer step");
    }
}
