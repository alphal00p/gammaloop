//! `GAMMALOOP_UNCONTRACTED_NUMERATOR=off|on|check[:<kinds>]` (U7a; the design
//! is results/consolidation/U7a-DESIGN.md of the beta project): an R*
//! remainder is handed to vakint's factorised entry FACTORISED and
//! UNCONTRACTED, and vakint contracts it (idenso) knowing the target topology.
//!
//! One record per ordinary remainder ([`super::evaluate`]'s
//! `integrate_ordinary`, the only hook): a `pint6.factorized-numerator.v3`
//! record of kind closure (WH-RECORD-FORMAT v1.4) built from what the call
//! site holds BEFORE any contraction: the graph's own edge and vertex
//! factors, the routed momentum vectors of the remainder's loop basis, the
//! closing line, the contracted children as `ct:` vertices carrying their Δ as
//! an insertion, the ports and our closure as the projector, the prefactor in
//! convention B and the outer Taylor step unapplied (form A). The record
//! makes a text round trip through vakint's frozen writer and reader, is
//! evaluated by `evaluate_factorized_closure` with our settings, and its value
//! goes through the same `RStarIntegrand::finish_tensor` as today's value,
//! under the same integration id.
//!
//! * Off (unset, `""`, `0`, `off`; the default): no record is built,
//!   nothing is allocated and nothing is printed; the hook calls today's body
//!   unchanged. Only the switch is read (a thread-local and a `LazyLock`), at
//!   the hook, the quotient scope, the Δ-cache key, the closed children's
//!   [`closure_data`] and the root's summary.
//! * On: the record's value replaces today's for that remainder. A DESIGNED
//!   refusal (the closed lists [`DESIGNED_LOCAL`] and [`DESIGNED_VAKINT`])
//!   falls back to today's path for that record only, with one
//!   `UNCONTRACTED ... outcome=declined:<code>` line. Any other refusal is
//!   `outcome=error:<code>` (a bug of the emitter or of vakint's entry) and
//!   fails the diagram: strict is the only behaviour of on and check (helper
//!   VH round 2 R2-I1; `GAMMALOOP_UNCONTRACTED_STRICT=1` is accepted as a
//!   no-op, `0` is refused when the switch initialises). Unknown master
//!   coefficients are a hard error (the format: an ALARM fails the diagram).
//! * Check: our attempt first (it changes no state), then today's path; the
//!   parts through the record's order (the poles for K1 / K2, `ε^(extra-1)`
//!   for K3) are compared algebraically (`UNCONTRACTED_CHECK`, which also
//!   records the record's depth, that order and the order of the first
//!   difference: U7-PLAN requirement 15); a difference at or below that
//!   order fails the diagram, and so does an error outcome (a designed
//!   refusal is the only way a record goes uncompared); today's value is
//!   returned. Check-mode times are NOT a timing comparison (our attempt runs
//!   first and warms vakint's caches): use paired off / on runs.
//!
//! Every diagram prints one `UNCONTRACTED_DISPOSITION` line per enabled kind
//! and one `UNCONTRACTED_SUMMARY` line (expected, attempted, uncontracted,
//! compared, equal, declined, no-record, Δ-cache-hit, unchecked and error
//! counts by code) with the switch on or check, and check mode fails a
//! diagram its comparisons do not cover: the ledger,
//! [`super::uncontracted_ledger`] (U7a delta, U7a.1 and U7a.4), which also
//! keys every record occurrence for paired off / on timings. The switch's
//! initialisation ([`validate`]) refuses `GAMMALOOP_UNCONTRACTED_STRICT=0`,
//! check mode with `PINT6_CLOSURE_DUMP`, a kind list without a kind this
//! build emits and a bad exception list.
//!
//! Kinds (`root` = K1 root remainders, `child` = K2 UV children, `ir-hard` =
//! K3 IR hard parts, and the inert `quotient`, `surrogate`, `insertion`):
//! the record's kind is the innermost special scope, else root / child. In
//! this build `root`, `child` and `ir-hard` build records (K1: one external
//! folded by σ, the root closure as projector; K2: every crown edge an
//! independent external, the closed child's closure as projector or none for
//! an open child; K3: below). In mode on the Δ-cache key carries
//! `#uncontracted=<enabled cacheable kinds>` ([`cache_key_suffix`]).
//!
//! K3 (U7b): one record per remainder term of an IR hard part `H(γ)` on the
//! child route (`GAMMALOOP_IR_HARD=child`; `evaluate.rs` `ir_hard_child`, at
//! its call of `integrate_ir_hard_term`): the hard part as its own support
//! root, every crown leg (`region.boundary`: soft lines and the scope's root
//! legs) an independent external whose `soft` line says `role=root` when it
//! is a root leg of the scope, the closing line at the FIXED reference mass
//! (`compensate_mass=0`, `orders=all`; the compensation is the caller's),
//! degree `ω`, depth = remainder loops + the quotient's extra loops (raised
//! on a depth retry), the UV children inside γ as `ct:` insertions. Its
//! ports are the open crown slots: Lorentz ones with no projector (gluonic
//! `h_s`), or the two bispinor slots of the open quark line closed by the IR
//! Dirac closure `gamma(c, r, mink(D, ir_dirac_closure_nu))` as the projector
//! (`GAMMALOOP_IR_DIRAC_CLOSURE`, slash kind; `nu` stays open, `1/Tr[1]` and
//! the reassembly stay the caller's, `ir_dirac::reassemble`). The value
//! replaces the return of `integrate_ir_hard_term`; check mode compares
//! through `ε^(depth - loops - 1)`, the orders the quotient product needs.
//! No record is built for a zero-loop term (Δ(γ) inserted, no vakint call
//! today) nor for a correction of order `ω = 0` (the `ir_term_today` path,
//! K8, out of scope; under the production
//! `GAMMALOOP_IR_HARD_CHILD_SELECT=omega1` it never reaches the child route).
//!
//! `GAMMALOOP_UNCONTRACTED_DUMP=<dir>` writes every record
//! (`<dir>/<kind>/<graph>-id<id>-<pid>-<seq>.record.txt`, a drop-in
//! `input.txt`) and its value (`.value.txt`, schema `u7a-value.v0`), at most
//! `GAMMALOOP_UNCONTRACTED_DUMP_MAX` records (default 20000) and 1 GiB.

use std::{
    cell::{Cell, RefCell},
    collections::{BTreeMap, BTreeSet},
    sync::{
        LazyLock,
        atomic::{AtomicUsize, Ordering},
    },
    time::Instant,
};

use eyre::{Result, eyre};
use linnet::half_edge::{
    involution::{EdgeIndex, Flow, HedgePair},
    subgraph::{SuBitGraph, SubSetOps},
};
use spenso::{
    network::tags::SPENSO_TAG,
    structure::representation::{Minkowski, RepName},
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    function,
};
use vakint::{
    Vakint, VakintSettings,
    p_integral::factorized::{
        ClosureRecordValue, ClosureSpec, Endpoint, FactorizedEdge, FactorizedError, FactorizedKind,
        FactorizedNumerator, FactorizedVertex, Insertion, LineKind, Prefactor, Projector, SoftLeg,
        TaylorOrders, TaylorSpec, evaluate_factorized_closure,
    },
};

use crate::{
    graph::Graph,
    model::{Model, ParticleIdGammaLoopExt, VertexRuleIdGammaLoopExt},
    numerator::symbolica_ext::NumeratorAtomExt,
    utils::{GS, W_},
    uv::{UVgenerationSettings, UltravioletGraph},
};

use super::{
    RStarPlan, ScopedSlot, UvPreparationTerm, boundary::ClosedRootNormalization,
    expression::RStarIntegrand, top_mass::TopMassTerms,
};

/// U7c: the K5 insertion-surrogate records (kind `surrogate`).
#[path = "uncontracted_k5.rs"]
pub(super) mod k5;

// ---------------------------------------------------------------------------
// The switch
// ---------------------------------------------------------------------------

/// The environment variable of the switch.
pub(crate) const VARIABLE: &str = "GAMMALOOP_UNCONTRACTED_NUMERATOR";

/// What the switch does with an enabled record.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Mode {
    /// No record is built (byte-identical to today).
    Off,
    /// vakint's value is used; a refusal falls back to today's path.
    On,
    /// Both are computed and compared; today's value is used.
    Check,
}

impl Mode {
    pub(crate) fn as_str(self) -> &'static str {
        match self {
            Self::Off => "off",
            Self::On => "on",
            Self::Check => "check",
        }
    }
}

/// The record kinds of the switch, in their canonical order.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) enum Kind {
    /// K1: an ordinary remainder of a diagram root (`a`).
    Root,
    /// K2: an ordinary remainder of a proper UV child evaluated by R* (`b`).
    Child,
    /// K3: a remainder term of an IR hard part `H(γ)` on the child route
    /// (`c`; U7b).
    IrHard,
    /// Inside the R* of a decorated IR quotient (`d`; inert here).
    Quotient,
    /// K5 (U7c): an insertion surrogate's records, its direct integral at
    /// <= 4 loops and the remainders of its own R* above (`e`; [`k5`]).
    Surrogate,
    /// Insertion normalisations (`f`; inert here).
    Insertion,
}

impl Kind {
    pub(crate) const ALL: [Kind; 6] = [
        Kind::Root,
        Kind::Child,
        Kind::IrHard,
        Kind::Quotient,
        Kind::Surrogate,
        Kind::Insertion,
    ];

    pub(crate) fn as_str(self) -> &'static str {
        match self {
            Self::Root => "root",
            Self::Child => "child",
            Self::IrHard => "ir-hard",
            Self::Quotient => "quotient",
            Self::Surrogate => "surrogate",
            Self::Insertion => "insertion",
        }
    }

    fn letter(self) -> &'static str {
        match self {
            Self::Root => "a",
            Self::Child => "b",
            Self::IrHard => "c",
            Self::Quotient => "d",
            Self::Surrogate => "e",
            Self::Insertion => "f",
        }
    }

    fn parse(text: &str) -> Option<Self> {
        Self::ALL
            .into_iter()
            .find(|kind| kind.as_str() == text || kind.letter() == text)
    }

    /// Whether this build emits records of this kind (the others are
    /// accepted and inert).
    pub(crate) fn implemented(self) -> bool {
        matches!(
            self,
            Self::Root | Self::Child | Self::IrHard | Self::Surrogate
        )
    }

    /// Whether `on` / `check` without a kind list enables this kind: the
    /// implemented kinds except `surrogate` (U7c), which must be listed
    /// (U8 decides the defaults).
    pub(crate) fn in_default(self) -> bool {
        self.implemented() && self != Self::Surrogate
    }

    fn bit(self) -> u8 {
        1 << (self as u8)
    }
}

/// `GAMMALOOP_UNCONTRACTED_NUMERATOR`, parsed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Switch {
    pub(crate) mode: Mode,
    /// The listed kinds (bits of [`Kind::bit`]); 0 when off.
    listed: u8,
}

impl Switch {
    pub(crate) const OFF: Switch = Switch {
        mode: Mode::Off,
        listed: 0,
    };

    /// A switch with `mode` and the listed `kinds` (tests).
    #[cfg(test)]
    pub(crate) fn new(mode: Mode, kinds: &[Kind]) -> Self {
        if mode == Mode::Off {
            return Self::OFF;
        }
        Self {
            mode,
            listed: kinds.iter().fold(0, |bits, kind| bits | kind.bit()),
        }
    }

    /// The switch of a value of the variable (`None`: unset).
    pub(crate) fn parse(value: Option<&str>) -> std::result::Result<Self, String> {
        let bad = |value: &str| {
            format!(
                "{VARIABLE} must be off, on or check[:<kinds>] (kinds: root|a, child|b, ir-hard|c, quotient|d, surrogate|e, insertion|f, comma separated), not `{value}`"
            )
        };
        let Some(value) = value.map(str::trim) else {
            return Ok(Self::OFF);
        };
        let (head, kinds) = match value.split_once(':') {
            Some((head, kinds)) => (head.trim(), Some(kinds)),
            None => (value, None),
        };
        let mode = match head {
            "" | "0" | "off" if kinds.is_none() => return Ok(Self::OFF),
            "on" | "1" => Mode::On,
            "check" => Mode::Check,
            _ => return Err(bad(value)),
        };
        let listed = match kinds {
            None => Kind::ALL
                .into_iter()
                .filter(|kind| kind.in_default())
                .fold(0, |bits, kind| bits | kind.bit()),
            Some(list) => {
                let mut bits = 0;
                for entry in list.split(',') {
                    let Some(kind) = Kind::parse(entry.trim()) else {
                        return Err(bad(value));
                    };
                    bits |= kind.bit();
                }
                bits
            }
        };
        Ok(Self { mode, listed })
    }

    /// Whether a record of `kind` is built under this switch.
    pub(crate) fn enabled(self, kind: Kind) -> bool {
        self.mode != Mode::Off && self.listed & kind.bit() != 0 && kind.implemented()
    }

    /// The enabled kinds, in canonical order.
    pub(crate) fn enabled_kinds(self) -> Vec<Kind> {
        Kind::ALL
            .into_iter()
            .filter(|kind| self.enabled(*kind))
            .collect()
    }

    /// The listed kinds this build does not emit (accepted and inert).
    pub(crate) fn inert_kinds(self) -> Vec<Kind> {
        Kind::ALL
            .into_iter()
            .filter(|kind| self.listed & kind.bit() != 0 && !kind.implemented())
            .collect()
    }

    /// `#uncontracted=<enabled cacheable kinds>` in mode on with at least one
    /// enabled kind whose values can reach the Δ cache (every kind but
    /// `root`: a diagram root's value never does, while an IR hard part's
    /// does, through the Δ of the child it belongs to); `None` otherwise,
    /// check mode included (every value is then today's), so the Δ-cache key
    /// is today's (U7a review M2).
    pub(crate) fn cache_key_suffix(self) -> Option<String> {
        if self.mode != Mode::On {
            return None;
        }
        let kinds = self
            .enabled_kinds()
            .into_iter()
            .filter(|kind| *kind != Kind::Root)
            .collect::<Vec<_>>();
        (!kinds.is_empty()).then(|| format!("#uncontracted={}", join_kinds(&kinds)))
    }

    /// The `UNCONTRACTED_SWITCH` line (printed once, at the first read).
    pub(crate) fn line(self) -> String {
        let inert = self.inert_kinds();
        format!(
            "UNCONTRACTED_SWITCH mode={} kinds={} inert={}",
            self.mode.as_str(),
            join_kinds(&self.enabled_kinds()),
            if inert.is_empty() {
                "none".to_owned()
            } else {
                join_kinds(&inert)
            }
        )
    }
}

fn join_kinds(kinds: &[Kind]) -> String {
    kinds
        .iter()
        .map(|kind| kind.as_str())
        .collect::<Vec<_>>()
        .join(",")
}

/// The switch of the environment and the check's exception list.
struct Environment {
    switch: Switch,
    exceptions: super::uncontracted_ledger::Exceptions,
}

/// The switch of the environment, read and checked once
/// ([`super::uncontracted_ledger::initialise`]: strict, the PINT6 dump, the
/// kinds, the exception list): an error refuses the run at [`validate`] (the
/// beta's start-up) and panics at the first use otherwise, as
/// `GAMMALOOP_INSERTION_CLOSURE` does.
static ENVIRONMENT: LazyLock<std::result::Result<Environment, String>> = LazyLock::new(|| {
    let read = |name: &str| match std::env::var(name) {
        Ok(value) => Ok(Some(value)),
        Err(std::env::VarError::NotPresent) => Ok(None),
        Err(std::env::VarError::NotUnicode(value)) => {
            Err(format!("{name}={value:?} is not unicode"))
        }
    };
    let switch = Switch::parse(read(VARIABLE)?.as_deref())?;
    let exceptions = super::uncontracted_ledger::initialise(
        switch,
        read(STRICT_VARIABLE)?.as_deref(),
        read(super::uncontracted_ledger::TIMING_VARIABLE)?.as_deref(),
        || super::uncontracted_ledger::Routes {
            uv_children_rstar: matches!(
                super::plan::UvChildrenRoute::from_environment(),
                super::plan::UvChildrenRoute::Rstar { .. }
            ),
            ir_hard_child: super::closed::child_mode() != super::closed::ChildMode::Off,
        },
        || super::pint6_dump::settings().is_some(),
        read(super::uncontracted_ledger::EXCEPTIONS_VARIABLE)?.as_deref(),
    )?;
    if switch.mode != Mode::Off {
        eprintln!(
            "{} strict=1 exceptions={}{}",
            switch.line(),
            match &exceptions.digest {
                Some(digest) => format!("{}:{digest}", exceptions.entries.len()),
                None => "none".to_owned(),
            },
            // a witness build only (feature `uncontracted-witness`)
            super::uncontracted_ledger::witness()
                .map_or_else(String::new, |witness| format!(" witness={witness}"))
        );
    }
    Ok(Environment { switch, exceptions })
});

fn environment() -> Switch {
    match &*ENVIRONMENT {
        Ok(environment) => environment.switch,
        Err(error) => panic!("{error}"),
    }
}

/// The switch's initialisation: `Err` for a value the run refuses (call at
/// start-up for an early error; the beta does).
pub(crate) fn validate() -> Result<()> {
    match &*ENVIRONMENT {
        Ok(_) => Ok(()),
        Err(error) => Err(eyre!("{error}")),
    }
}

/// The exception list of the environment (empty unless check mode reads
/// one).
pub(super) fn exceptions() -> super::uncontracted_ledger::Exceptions {
    match &*ENVIRONMENT {
        Ok(environment) => environment.exceptions.clone(),
        Err(error) => panic!("{error}"),
    }
}

thread_local! {
    static TEST_SWITCH: Cell<Option<Switch>> = const { Cell::new(None) };
}

/// The switch in force on this thread (the environment's, or a test's).
pub(crate) fn switch() -> Switch {
    TEST_SWITCH.with(Cell::get).unwrap_or_else(environment)
}

/// Tests: run `f` with `switch` on this thread, never through the
/// environment (restored on unwinding too).
#[cfg(test)]
pub(crate) fn with_switch<T>(switch: Switch, f: impl FnOnce() -> T) -> T {
    struct Restore(Option<Switch>);
    impl Drop for Restore {
        fn drop(&mut self) {
            let previous = self.0;
            TEST_SWITCH.with(|slot| slot.set(previous));
        }
    }
    let _restore = Restore(TEST_SWITCH.with(|slot| slot.replace(Some(switch))));
    f()
}

/// The Δ-cache key suffix of the switch in force (`None` when off), and,
/// on or check, `#uncontracted-inert=<kind>` for an R* child Δ looked up
/// inside a scope whose records are a kind the switch does not enable (an
/// insertion surrogate's R*, a decorated quotient's R*, an IR hard part
/// with `ir-hard` not listed) while a cacheable kind (`child`, `ir-hard`) is
/// enabled: such a fill builds no record, so it must never serve the same
/// child of an enabled kind, which would then be neither compared nor
/// counted (review of the U7a delta, refusal-integrity lens, finding 1).
/// Off (and with no cacheable kind enabled): nothing, today's key.
pub(crate) fn cache_key_suffix() -> Option<String> {
    let switch = switch();
    let mut suffix = switch.cache_key_suffix();
    // only when a fill can hold records at all (a cacheable kind enabled):
    // without one no lookup counts a hit, and the key stays today's
    if switch.mode != Mode::Off
        && (switch.enabled(Kind::Child)
            || switch.enabled(Kind::IrHard)
            || switch.enabled(Kind::Surrogate))
    {
        let kind = classify(true);
        if !switch.enabled(kind) {
            suffix
                .get_or_insert_with(String::new)
                .push_str(&format!("#uncontracted-inert={}", kind.as_str()));
        }
    }
    suffix
}

// ---------------------------------------------------------------------------
// Scopes and kinds
// ---------------------------------------------------------------------------

/// A special scope that makes the records inside it their own kind.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Scope {
    /// The R* of a decorated IR quotient (`quotient_rstar`).
    Quotient,
}

thread_local! {
    static SCOPE: RefCell<Option<Scope>> = const { RefCell::new(None) };
}

/// Run `f` inside `scope` (restored on unwinding too). With the switch off it
/// just calls `f`, unless `GAMMALOOP_UNCONTRACTED_TIMING=1` keys the off
/// run's occurrences, which then need the kinds an on run gives them.
pub(super) fn scoped<T>(scope: Scope, f: impl FnOnce() -> T) -> T {
    if switch().mode == Mode::Off && !super::uncontracted_ledger::timing_enabled() {
        return f();
    }
    let slot = ScopedSlot::enter(&SCOPE, scope);
    let out = f();
    slot.exit();
    out
}

fn in_scope(scope: Scope) -> bool {
    SCOPE.with(|slot| *slot.borrow() == Some(scope))
}

/// The kind of a K3 record (the IR hard part hook of `ir_hard_child`):
/// `ir-hard`, unless the hard part belongs to the R* of a decorated IR
/// quotient or of an insertion surrogate, whose kinds take every record
/// inside them (U7c / U7d enable those, with their own vertices); so
/// `ir-hard` alone never reaches a graph that carries a decorated or
/// surrogate vertex.
pub(super) fn classify_ir_hard() -> Kind {
    if in_scope(Scope::Quotient) {
        Kind::Quotient
    } else if super::expression::dressed_scope() {
        Kind::Surrogate
    } else {
        Kind::IrHard
    }
}

/// The kind of a record at the hook: the innermost special scope, else root
/// (a diagram root, `root_boundary` none) or child (a proper UV child).
pub(super) fn classify(support_root: bool) -> Kind {
    if super::expression::ir_hard_fixed_mass() {
        Kind::IrHard
    } else if in_scope(Scope::Quotient) {
        Kind::Quotient
    } else if super::expression::dressed_scope() {
        Kind::Surrogate
    } else if support_root {
        Kind::Child
    } else {
        Kind::Root
    }
}

// ---------------------------------------------------------------------------
// The call site
// ---------------------------------------------------------------------------

/// The closure of a closed K2 child, as data (today it exists only inside
/// the opaque normalisation callback and, when admitted, in `closed_root`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum ChildClosure {
    Isotropic,
    DiracVertex,
    DiracSelfEnergy,
}

impl ChildClosure {
    pub(super) fn as_str(self) -> &'static str {
        match self {
            Self::Isotropic => "isotropic",
            Self::DiracVertex => "dirac_vertex",
            Self::DiracSelfEnergy => "dirac_self_energy",
        }
    }
}

/// A closed child's closure: its kind, its normalisation and the slot it
/// leaves open (the Dirac self-energy's fresh `nu`).
#[derive(Clone, Debug)]
pub(super) struct ClosureData {
    pub kind: ChildClosure,
    pub normalization: Atom,
    /// The open slot (diagnostics; vakint derives `open` from the record).
    #[allow(dead_code)]
    pub open: Option<Atom>,
}

/// A closed child's closure as data for its records (read only at the hook,
/// on the evaluation's own thread); `None` with the switch off, so nothing is
/// cloned on today's path.
pub(super) fn closure_data(
    kind: ChildClosure,
    normalization: &Atom,
    open: Option<&Atom>,
) -> Option<ClosureData> {
    (switch().mode != Mode::Off).then(|| ClosureData {
        kind,
        normalization: normalization.clone(),
        open: open.cloned(),
    })
}

/// A borrowed view of one evaluation at the hook ([`super::evaluate`] builds
/// it): this module never touches the evaluation's own fields.
pub(super) struct Site<'a> {
    pub graph: &'a Graph,
    pub model: &'a Model,
    pub plan: &'a RStarPlan,
    pub term: &'a UvPreparationTerm,
    pub children: &'a BTreeMap<SuBitGraph, Atom>,
    pub root_closure: Option<&'a Atom>,
    pub closed_root: Option<&'a ClosedRootNormalization>,
    pub closure_data: Option<&'a ClosureData>,
    /// Whether the evaluation has a root projector (callback).
    pub root_projector: bool,
    /// The carriers of the outer Taylor step (the K2 record's crown).
    pub root_edges: &'a [EdgeIndex],
    pub reference_mass_squared: Symbol,
    pub auxiliary_mass_squared: &'a Atom,
    pub settings: &'a VakintSettings,
    pub vakint: &'a Vakint,
    pub uv_settings: &'a UVgenerationSettings,
    pub top_mass: Option<TopMassTerms>,
    /// The settings' ε (diagnostics).
    #[allow(dead_code)]
    pub epsilon: Symbol,
    /// Which route today's path takes (`boundary`, `whole`,
    /// `declined:<admission reason>`; `ir_hard` for K3), for the log.
    pub today: String,
    /// A K3 record's IR hard part (`None` at the ordinary hook).
    pub ir_hard: Option<IrHardSite<'a>>,
}

/// What the IR hard part hook (`ir_hard_child`) adds to a [`Site`] for one
/// remainder term of `H(γ)`.
#[derive(Clone, Copy, Debug)]
pub(super) struct IrHardSite<'a> {
    /// The correction, its term and the depth attempt (the label).
    pub correction: usize,
    pub term: usize,
    pub attempt: usize,
    /// The hard part's edge ids (the label).
    pub hard: &'a str,
    /// `ω`, the Taylor degree.
    pub degree: usize,
    /// The Laurent terms of the returned value: the remainder's loops plus
    /// the quotient's extra loops (`integrate_ir_hard_term`).
    pub depth: usize,
    /// The IR Dirac closure multiplied into the remainder (the projector).
    pub closure: Option<&'a super::expression::DiracBoundary>,
}

impl Site<'_> {
    /// The record's depth: K3 its own, else the remainder's loops.
    fn depth(&self) -> usize {
        self.ir_hard
            .map_or(self.term.remainder.loops, |ir_hard| ir_hard.depth)
    }
}

/// The Laurent order through which a value of `depth` terms over `loops`
/// loops is complete (`ε^(depth - loops - 1)`; `ε^-1` for K1 / K2).
pub(super) fn valid_order(depth: usize, loops: usize) -> i64 {
    depth as i64 - loops as i64 - 1
}

/// A record that is not built or not evaluated: today's path takes it when
/// the refusal is DESIGNED ([`DESIGNED_LOCAL`], [`DESIGNED_VAKINT`]); any
/// other is an error (U7a review B1).
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct Decline {
    pub code: String,
    pub detail: String,
    /// A designed refusal (the closed lists), not a bug.
    pub designed: bool,
}

/// The emitter's designed refusals (U7a review B1), a CLOSED list: a record
/// the format or this build does not state. Every other local code
/// (`routing`, `crown` (a K2 / K3 crown that is not the remainder's soft
/// edges, or a graph whose edge ids reach [`PAIRED_LEG_OFFSET`]; helper VH
/// round 2, R2-M4 / R2-M5), `missing_child`, `product_mismatch`, `text`,
/// `taylor`, `closure_data_missing`, `normalisation_mismatch`,
/// `finish_tensor`, `vakint:Empty`, and K3's `ir_hard_site_missing`,
/// `ir_dirac_closure`, ...) is a bug of the emitter or of the site:
/// `outcome=error:<code>`. So are the record builder's states `sigma` (σ
/// does not fold a root leg onto ±p(1)), `root_crown_paired` (a K1 crown
/// edge with both ends in the root) and `root_normalisation_momentum` (a
/// closed-root factor that is not a scalar times a power of p²): they never
/// occur on a correct site (helper VH round 2, R2-I3).
///
/// `one_loop_tree_cut` (U7-PLAN 14(d), narrowed as Astra U7a asks: "keep the
/// rebase test narrow"): ONLY a remainder whose one internal propagator is
/// the closing line ([`single_propagator`]; after the Taylor step a massive
/// tadpole, which vakint's closure entry refuses as "every line is
/// massive"; B_self_energy_L4_0049 id 8), of any kind. Every other one-loop
/// remainder is a record, although today's route names each of them
/// `one_loop_tree_cut` (it integrates them whole).
///
/// K3 (U7b) adds NO code: a K3 record can hit only rows already here,
/// `insertion_mass_symbol` (a UV child's Δ with a mass symbol at a `ct:`
/// node where independent crown legs meet, as K2),
/// `insertions_share_a_vertex`, `pint6_dump` and `one_loop_tree_cut` (a
/// one-loop term whose only line is the closing one); its zero-loop terms
/// and its order-0 corrections are not attempted at all ([`NoRecord`]); the
/// root-only codes (`open_root`, `root_normalisation_opaque`) and
/// `top_mass_mode` cannot occur (its crown is independent, its mass is fixed
/// whatever the top-mass mode).
pub(super) const DESIGNED_LOCAL: &[&str] = &[
    "pint6_dump",
    "zero_loops",
    "kind_not_built",
    "child_not_built",
    "open_root",
    "root_normalisation_opaque",
    "additional_normalization",
    "top_mass_mode",
    "insertion_mass_symbol",
    "insertions_share_a_vertex",
    "no_closing_line",
    "one_loop_tree_cut",
    // U7c (K5): a K3-shaped IR hard term inside an insertion surrogate's R*
    // (no verified fixture; Astra #10 keeps such terms out of `ir-hard`):
    // declined before its record is built, today's path takes it
    "surrogate_ir_hard",
];

/// vakint's designed refusals, matched by the error variant AND the prefix of
/// its detail AND a phrase the detail must contain (`(variant, detail prefix,
/// phrase, code)`, phrase `""`: none; the code is `vakint:<code>`): the rows
/// of WH-RECORD-FORMAT v1.4 section 8 a K1 / K2 record can reach. The
/// open-spinor row needs its phrase `is not a Lorentz slot` (helper VH round
/// 2, R2-I4): the prefix `open slot ` alone also starts vakint's vacuum-scope
/// refusal "open slot <port>: only Lorentz slots ... stay open".
///
/// Interim rule (U7-PLAN requirement 13, until the rebase onto helper VH's
/// frozen REFUSALS list): a `Closure` refusal is designed ONLY with a
/// section-8 prefix (`Taylor degree `, `the projector contains Q(`, the
/// open-spinor port `open slot `); any other `Closure` detail is an error,
/// among them the entry's wrapped internal failures and "every line is
/// massive: not a closure" (vakint closure.rs `prepare`: a remainder whose
/// only line is the closing one; NOT a section-8 row; the emitter declines
/// that shape itself, `one_loop_tree_cut`, so vakint refusing it means the
/// emitter's shape test missed one). Every other vakint error (the reader's
/// `Text*`, `IndexUncontracted`, `MomentumNotConserved`, ...) is an error:
/// `outcome=error:vakint:<Variant>`.
///
/// K3 (U7b) adds NO row. Of these, a K3 record can legitimately hit
/// `taylor_degree` (`ω > 2`; r6 has `ω <= 2`), `open_spinor` (a hard part
/// whose open quark line the IR Dirac closure does not close:
/// `GAMMALOOP_IR_DIRAC_CLOSURE=off`, or a declined shape other than r6's
/// `no_spinor`), `insertion_mass_symbol`, `insertion_slots`, `gamma5` and
/// `levi_civita` (a one-loop term whose only remaining line is the closing
/// one is the emitter's `one_loop_tree_cut`); `projector_momentum` and
/// `compact_slash` cannot occur (the IR Dirac closure
/// `gamma(c, r, mink(D, nu))` is momentum free, and a local check makes a
/// momentum in it an error, `ir_dirac_closure`).
pub(super) const DESIGNED_VAKINT: &[(&str, &str, &str, &str)] = &[
    ("Closure", "Taylor degree ", "", "taylor_degree"),
    (
        "Closure",
        "the projector contains Q(",
        "",
        "projector_momentum",
    ),
    ("Ports", "a compact slash gamma(", "", "compact_slash"),
    (
        "Closure",
        "open slot ",
        "is not a Lorentz slot",
        "open_spinor",
    ),
    (
        "Insertion",
        "contains a mass symbol",
        "",
        "insertion_mass_symbol",
    ),
    ("Insertion", "the factor has ", "", "insertion_slots"),
    ("Gamma5", "", "", "gamma5"),
    ("LeviCivita", "", "", "levi_civita"),
];

/// Whether `code` is the code of a designed refusal as an outcome names it
/// (`declined:<code>`): a [`DESIGNED_LOCAL`] code or `vakint:<code>` of a
/// [`DESIGNED_VAKINT`] row. Only these can be named by an exception of a
/// check run ([`super::uncontracted_ledger::parse_exceptions`]): an error is
/// never excepted.
pub(super) fn is_designed_code(code: &str) -> bool {
    DESIGNED_LOCAL.contains(&code)
        || is_designed_section_code(code)
        || code
            .strip_prefix("vakint:")
            .is_some_and(|code| DESIGNED_VAKINT.iter().any(|row| row.3 == code))
}

/// The designed code a K3 record gets when its IR correction declines as a
/// whole after the record ran ([`super::uncontracted_ledger::Section`]):
/// only `ir_hard_guard:budget` (the exhausted budget of
/// `GAMMALOOP_IR_HARD_CHILD_GUARD=1`, a resource decline; today's closed
/// path computes the correction). Never an attempt's outcome.
///
/// Pre-review of the U7a delta (findings 3 and 4): the route's own declines
/// are NOT designed codes of a record. In mode on a route decline after its
/// K3 records were taken is caused by (or at least follows) the new values
/// and is the error `ir_hard_route:<reason head>`; a panic or an error the
/// guard caught is the error `ir_hard_guard:<panic|error>`, inside a record
/// or not. Neither can be named by an exception.
pub(super) fn is_designed_section_code(code: &str) -> bool {
    code == "ir_hard_guard:budget"
}

impl Decline {
    fn new(code: impl Into<String>, detail: impl Into<String>) -> Self {
        let code = code.into();
        Self {
            designed: DESIGNED_LOCAL.contains(&code.as_str()),
            code,
            detail: detail.into(),
        }
    }

    /// A refusal of vakint's frozen reader or entry: designed when its
    /// variant, detail prefix and phrase are a row of [`DESIGNED_VAKINT`],
    /// else an error by variant.
    pub(super) fn vakint(error: &FactorizedError) -> Self {
        let (name, detail) = variant_detail(error);
        let row = DESIGNED_VAKINT.iter().find(|(variant, prefix, phrase, _)| {
            *variant == name
                && detail
                    .is_some_and(|detail| detail.starts_with(prefix) && detail.contains(phrase))
        });
        match row {
            Some((_, _, _, code)) => Self {
                code: format!("vakint:{code}"),
                detail: error.to_string(),
                designed: true,
            },
            None => Self {
                code: format!("vakint:{name}"),
                detail: error.to_string(),
                designed: false,
            },
        }
    }
}

/// The variant name of a vakint factorised error (`Closure`, `Ports`, ...)
/// and its `detail` field ("" for the variants matched by name only).
fn variant_detail(error: &FactorizedError) -> (String, Option<&str>) {
    match error {
        FactorizedError::Closure { detail } => ("Closure".into(), Some(detail)),
        FactorizedError::Ports { detail } => ("Ports".into(), Some(detail)),
        FactorizedError::Insertion { detail, .. } => ("Insertion".into(), Some(detail)),
        FactorizedError::Gamma5 { .. } => ("Gamma5".into(), Some("")),
        FactorizedError::LeviCivita { .. } => ("LeviCivita".into(), Some("")),
        other => (
            format!("{other:?}")
                .split([' ', '{', '('])
                .next()
                .unwrap_or("Unknown")
                .to_owned(),
            None,
        ),
    }
}

/// Whether vakint's error text is the unknown-master refusal (strict
/// unknown-master policy, `UnknownMasterCoefficients`) or its ALARM
/// (`MissingMasterOrders`): a hard error of the record and the diagram, never
/// a fallback (WH-RECORD-FORMAT v1.4 section 7). By text until vakint gives a
/// typed predicate.
pub(super) fn is_unknown_master_text(text: &str) -> bool {
    text.contains("expected master coefficient(s) that are not known yet")
        || text.contains("needs master Laurent coefficients beyond the available expansion depth")
}

/// The hard error of vakint's unknown-master refusal of a record, a
/// diagnostic failure: no fallback may decline it (the IR hard child guard
/// `GAMMALOOP_IR_HARD_CHILD_GUARD=1`, the Δ cache; review of the U7a delta,
/// refusal-integrity finding 2).
pub(super) fn unknown_masters_failure(
    kind: Kind,
    graph: &str,
    id: usize,
    text: &str,
) -> eyre::Report {
    super::ir_dirac::diagnostic_failure(eyre!(
        "UNCONTRACTED kind={} graph={graph} id={id}: vakint refuses the record's unknown master coefficients: {text}",
        kind.as_str()
    ))
}

/// Whether vakint's error text is an exhausted scoped resource budget (the
/// RSS budget of symbolica's `ResourceBudget`, as the closed path's guard
/// detects it by text: `resource budget exhausted`, `exceeds RAM budget`).
pub(super) fn is_budget_text(text: &str) -> bool {
    text.contains("resource budget exhausted") || text.contains("exceeds RAM budget")
}

/// What one attempt gave.
#[derive(Clone, Debug)]
pub(super) enum Outcome {
    /// vakint's value, after `finish_tensor` (today's notation).
    Value(Atom),
    /// A child counterterm Δ is literally zero: the value is 0, no call.
    Zero,
    /// A designed refusal: today's path takes the record.
    Declined(Decline),
    /// Any other refusal (a bug on either side): fails the diagram in check
    /// mode and in mode on (strict is the only behaviour, R2-I1).
    Error(Decline),
}

impl Outcome {
    fn refused(decline: Decline) -> Self {
        if decline.designed {
            Self::Declined(decline)
        } else {
            Self::Error(decline)
        }
    }
}

/// One record attempt at the hook, with everything its log line, its dump
/// and the check need.
#[derive(Clone, Debug)]
pub(super) struct Attempt {
    pub kind: Kind,
    pub graph: String,
    pub id: usize,
    pub loops: usize,
    /// The record's depth and the order `ε^order` through which its value
    /// is complete (and compared in check mode).
    pub depth: usize,
    pub order: i64,
    pub closure: &'static str,
    /// The record's text (what vakint read), when built.
    pub text: Option<String>,
    pub outcome: Outcome,
    /// The log fields after `loops=` of a record vakint evaluated.
    fields: String,
    /// Our build (routing, record, round trip), vakint's entry and our
    /// `finish_tensor`, in seconds (0 when not reached).
    pub build_seconds: f64,
    pub evaluate_seconds: f64,
    pub finish_seconds: f64,
    /// vakint::timing over the entry.
    pub delta: Option<vakint::timing::Delta>,
    /// Insertions rewritten by the Δ momentum-conservation step.
    pub conserved: usize,
}

impl Attempt {
    /// The value replacing today's: vakint's, or 0 for a zero insertion.
    pub(super) fn value(&self) -> Option<Atom> {
        match &self.outcome {
            Outcome::Value(value) => Some(value.clone()),
            Outcome::Zero => Some(Atom::Zero),
            Outcome::Declined(_) | Outcome::Error(_) => None,
        }
    }

    /// Whether the attempt is an error (not a designed refusal).
    pub(super) fn is_error(&self) -> bool {
        matches!(self.outcome, Outcome::Error(_))
    }

    pub(super) fn outcome_name(&self) -> String {
        match &self.outcome {
            Outcome::Value(_) => "uncontracted".into(),
            Outcome::Zero => "zero".into(),
            Outcome::Declined(decline) => format!("declined:{}", decline.code),
            Outcome::Error(decline) => format!("error:{}", decline.code),
        }
    }
}

/// Our attempt at one remainder of `kind` (integration id `id`). The routed
/// remainder is built by `remainder` (the untruncated `new_with_closure`
/// remainder today's path routes; K3: the one `ir_hard_child` holds) only
/// when a record is built. `Err`: a hard error (unknown masters), which fails
/// the diagram.
pub(super) fn attempt<R: std::borrow::Borrow<RStarIntegrand>>(
    site: &Site<'_>,
    kind: Kind,
    id: usize,
    mode: Mode,
    remainder: impl FnOnce() -> Result<R>,
) -> Result<Attempt> {
    let started = Instant::now();
    let region = &site.term.remainder;
    let depth = site.depth();
    let mut attempt = Attempt {
        kind,
        graph: site.graph.name.clone(),
        id,
        loops: region.loops,
        depth,
        order: valid_order(depth, region.loops),
        closure: closure_name(site, kind),
        text: None,
        outcome: Outcome::Error(Decline::new("unbuilt", "")),
        fields: String::new(),
        build_seconds: 0.0,
        evaluate_seconds: 0.0,
        finish_seconds: 0.0,
        delta: None,
        conserved: 0,
    };
    let built = if super::pint6_dump::settings().is_some() {
        Err(Decline::new(
            "pint6_dump",
            "PINT6 dumps keep today's records",
        ))
    } else if !kind.implemented() {
        Err(Decline::new("kind_not_built", kind.as_str()))
    } else if region.loops == 0 {
        Err(Decline::new("zero_loops", "a zero-loop remainder"))
    } else {
        match remainder() {
            Ok(remainder) => build(site, remainder.borrow(), kind, id, mode),
            Err(error) => Err(Decline::new("routing", format!("{error:#}"))),
        }
    };
    let built = match built {
        Ok(Some(built)) => built,
        Ok(None) => {
            attempt.outcome = Outcome::Zero;
            return Ok(attempt);
        }
        Err(decline) => {
            attempt.outcome = Outcome::refused(decline);
            attempt.build_seconds = started.elapsed().as_secs_f64();
            return Ok(attempt);
        }
    };
    let build_seconds = started.elapsed().as_secs_f64();
    attempt.build_seconds = build_seconds;
    attempt.text = Some(built.text.clone());
    attempt.conserved = built.conserved;
    let mut settings = site.settings.clone();
    settings.number_of_terms_in_epsilon_expansion = i64::try_from(depth)?;
    let prefactor = built.record.prefactor.total();
    let counted = vakint::timing::Snapshot::now();
    let evaluate_started = Instant::now();
    let result = evaluate_factorized_closure(
        site.vakint,
        &settings,
        std::slice::from_ref(&built.record),
        std::slice::from_ref(&prefactor),
    );
    let delta = counted.elapsed();
    let evaluate_seconds = evaluate_started.elapsed().as_secs_f64();
    attempt.evaluate_seconds = evaluate_seconds;
    attempt.delta = Some(delta);
    let value = match result {
        Ok(value) => value,
        Err(error) => {
            let text = error.to_string();
            if is_budget_text(&text) {
                // an exhausted scoped resource budget inside vakint (only
                // under GAMMALOOP_IR_HARD_CHILD_GUARD=1's budget): not a
                // refusal of the record; a plain error carrying the budget
                // phrase, which the guard declines as `budget` (the record
                // then is `declined:ir_hard_guard:budget`) and which fails
                // the diagram anywhere else
                return Err(eyre!(
                    "resource budget exhausted in the uncontracted record kind={} graph={} id={id}: {text}",
                    kind.as_str(),
                    attempt.graph
                ));
            }
            if is_unknown_master_text(&text) {
                tally_unknown_masters(kind);
                eprintln!(
                    "UNCONTRACTED kind={} graph={} id={id} outcome=error:unknown_masters detail=\"{}\"",
                    kind.as_str(),
                    attempt.graph,
                    escape(&text)
                );
                return Err(unknown_masters_failure(kind, &attempt.graph, id, &text));
            }
            attempt.outcome = Outcome::refused(Decline::vakint(&error));
            return Ok(attempt);
        }
    };
    let Some(record) = value.records.into_iter().next() else {
        attempt.outcome = Outcome::refused(Decline::new("vakint:Empty", "no record value"));
        return Ok(attempt);
    };
    let finish_started = Instant::now();
    let finished = match RStarIntegrand::finish_tensor(record.value.clone(), id) {
        Ok(finished) => finished,
        Err(error) => {
            attempt.outcome = Outcome::refused(Decline::new("finish_tensor", format!("{error:#}")));
            return Ok(attempt);
        }
    };
    let finish_seconds = finish_started.elapsed().as_secs_f64();
    attempt.finish_seconds = finish_seconds;
    let wall_seconds = started.elapsed().as_secs_f64();
    attempt.fields = timing_fields(
        site,
        &built,
        &record,
        &delta,
        build_seconds,
        finish_seconds,
        wall_seconds,
    );
    if vakint::timing::record_log_enabled() {
        eprintln!(
            "{}",
            vakint::timing::record_line(
                kind.as_str(),
                &built.record.label,
                region.loops,
                Some(i64::try_from(depth)?),
                Some(record.patterns),
                &delta,
                evaluate_seconds,
                built.text.len(),
            )
        );
    }
    attempt.outcome = Outcome::Value(finished);
    Ok(attempt)
}

/// The closure of the record's label and log line, from the actual closure:
/// a closed child's kind; the root's `metric` when its closure is one metric
/// `g(a, b)`, else `closed`; `none` without a closure.
fn closure_name(site: &Site<'_>, kind: Kind) -> &'static str {
    let metric = |closure: &Atom| {
        matches!(closure.as_view(), AtomView::Fun(f)
            if f.get_nargs() == 2 && f.get_symbol().get_name().ends_with("::g"))
    };
    if kind == Kind::IrHard {
        return match site.ir_hard.and_then(|ir_hard| ir_hard.closure) {
            Some(boundary) => match boundary.kind {
                super::expression::DiracKind::SelfEnergy => "ir_dirac_slash",
                super::expression::DiracKind::Vertex => "ir_dirac_vertex",
            },
            None => "none",
        };
    }
    match (kind, site.closure_data, site.root_closure) {
        (Kind::Child, Some(data), _) => data.kind.as_str(),
        (_, _, Some(closure)) if metric(closure) => "metric",
        (_, _, Some(_)) => "closed",
        _ => "none",
    }
}

#[allow(clippy::too_many_arguments)]
fn timing_fields(
    site: &Site<'_>,
    built: &Built,
    record: &ClosureRecordValue,
    delta: &vakint::timing::Delta,
    build_seconds: f64,
    finish_seconds: f64,
    wall_seconds: f64,
) -> String {
    let rss = vakint::timing::Stage::ALL
        .iter()
        .filter_map(|stage| delta.stage(*stage).peak_rss_bytes)
        .max()
        .map_or_else(
            || "-".to_owned(),
            |bytes| format!("{:.2}", bytes as f64 / (1u64 << 30) as f64),
        );
    format!(
        " depth={} degree={} externals={} insertions={} conserved={} today={} patterns={} blocks={} key_plans={} lowered_entries={} build_s={build_seconds:.3} finish_s={finish_seconds:.3}{}{} wall_s={wall_seconds:.3} rss_gib={rss} record_bytes={}",
        site.depth(),
        built.degree,
        built.externals,
        built.insertions,
        built.conserved,
        site.today,
        record.patterns,
        record.blocks,
        record.key_plans,
        record.lowered_entries,
        delta.fields(),
        delta.stage_fields(),
        built.text.len()
    )
}

/// A log detail: quotes and line breaks escaped, at most 200 characters.
fn escape(text: &str) -> String {
    let escaped = text
        .replace('\\', "\\\\")
        .replace('"', "\\\"")
        .replace('\n', "\\n")
        .replace('\t', " ");
    match escaped.char_indices().nth(200) {
        Some((cut, _)) => format!("{}...", &escaped[..cut]),
        None => escaped,
    }
}

// ---------------------------------------------------------------------------
// Check mode
// ---------------------------------------------------------------------------

/// The comparison of today's value with ours (check mode).
#[derive(Clone, Debug)]
pub(super) struct Check {
    /// Structurally equal.
    pub identical: bool,
    /// The parts through ε^order agree (expanded difference, `together`,
    /// the tensor canonical form): the criterion, since only those reach a
    /// consumer (the poles for K1 / K2).
    pub equal: bool,
    /// The whole values agree (informational: ε^0 pieces are incomplete).
    pub full: bool,
    /// The lowest Laurent order at which the values differ (U7-PLAN
    /// requirement 15, the built-in depth check): recorded for every
    /// `identical=false` comparison; at or below the compared order it fails
    /// the check.
    pub first_difference: FirstDifference,
    pub contracted: Atom,
    pub contracted_seconds: f64,
    pub uncontracted_seconds: f64,
    /// Why the comparison itself failed, if it did.
    pub error: Option<String>,
}

/// The lowest Laurent order at which today's value and ours differ.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum FirstDifference {
    /// The values agree in full (or are identical).
    None,
    /// The lowest order `ε^k` whose parts through it differ.
    At(i64),
    /// They differ only above `ε^k`, the last order scanned.
    Above(i64),
}

impl FirstDifference {
    /// The `first_difference=` value of the check line: `none`, `eps^k` or
    /// `above:eps^k`.
    pub(crate) fn as_text(self) -> String {
        match self {
            Self::None => "none".to_owned(),
            Self::At(order) => format!("eps^{order}"),
            Self::Above(order) => format!("above:eps^{order}"),
        }
    }
}

/// The orders scanned above the compared one for the first difference of two
/// values that agree through it (a record's value carries at most a few more).
const FIRST_DIFFERENCE_SCAN: i64 = 8;

/// The Laurent part through ε^order (the pole part, as `Evaluation::poles`,
/// for `order = -1`).
pub(super) fn through(value: &Atom, epsilon: Symbol, order: i64) -> Result<Atom> {
    Ok(value.series(epsilon, 0, order)?.to_atom().expand())
}

/// The lowest order at which `today` and `ours` (not equal in full) differ:
/// above `order` when they agree through it (`equal`), else from the lowest
/// order of their expanded difference (below which every part agrees term by
/// term), each order by [`through`] and [`agree`].
fn first_difference(
    today: &Atom,
    ours: &Atom,
    epsilon: Symbol,
    order: i64,
    equal: bool,
) -> Result<FirstDifference> {
    let start = if equal {
        order + 1
    } else {
        let difference = (today - ours).expand();
        let series = difference.series(
            epsilon,
            0,
            symbolica::poly::series::SeriesDepth::relative(1),
        )?;
        if series.is_zero() {
            order + 1
        } else {
            let exponent = series.get_trailing_exponent();
            exponent
                .is_integer()
                .then(|| exponent.numerator().to_i64())
                .flatten()
                .map_or(order - FIRST_DIFFERENCE_SCAN, |trailing| {
                    trailing.min(order + 1)
                })
        }
    };
    let end = order + FIRST_DIFFERENCE_SCAN;
    for k in start..=end {
        if !agree(&through(today, epsilon, k)?, &through(ours, epsilon, k)?)? {
            return Ok(FirstDifference::At(k));
        }
    }
    Ok(FirstDifference::Above(end))
}

/// Whether two values agree after expansion, `together` and the tensor
/// canonical form ([`super::insertions::insertion_coefficients_agree`]).
pub(super) fn agree(a: &Atom, b: &Atom) -> Result<bool> {
    if a == b {
        return Ok(true);
    }
    Ok(super::insertions::insertion_coefficients_agree(a, b)?.1)
}

/// Compare today's value `today` with ours `ours` (both after
/// `finish_tensor` under the same id) through ε^order: the poles (`-1`) for
/// K1 / K2, whose consumers take only those; K3 through `ε^(extra - 1)`,
/// the orders its product with the quotient counterterms reaches
/// ([`Attempt::order`]).
pub(super) fn compare(
    today: &Atom,
    ours: &Atom,
    epsilon: Symbol,
    order: i64,
    contracted_seconds: f64,
    uncontracted_seconds: f64,
) -> Check {
    let identical = today == ours;
    let mut error = None;
    let equal = identical
        || match through(today, epsilon, order)
            .and_then(|a| through(ours, epsilon, order).map(|b| (a, b)))
        {
            Ok((a, b)) => agree(&a, &b).unwrap_or_else(|e| {
                error = Some(format!("{e:#}"));
                false
            }),
            Err(e) => {
                error = Some(format!("{e:#}"));
                false
            }
        };
    let full = identical || agree(today, ours).unwrap_or(false);
    // requirement 15: the order of the first difference of every
    // identical=false comparison; at or below `order` the check fails
    let first = if full {
        FirstDifference::None
    } else {
        first_difference(today, ours, epsilon, order, equal).unwrap_or_else(|e| {
            error.get_or_insert_with(|| format!("first difference: {e:#}"));
            FirstDifference::Above(order)
        })
    };
    let equal = equal && !matches!(first, FirstDifference::At(k) if k <= order);
    Check {
        identical,
        equal,
        full,
        first_difference: first,
        contracted: today.clone(),
        contracted_seconds,
        uncontracted_seconds,
        error,
    }
}

/// The error of a failed check (`equal=false`): it fails the diagram.
pub(super) fn check_failure(attempt: &Attempt, check: &Check, epsilon: Symbol) -> eyre::Report {
    let ours = attempt.value().unwrap_or_default();
    let show = |value: &Atom| {
        through(value, epsilon, attempt.order)
            .map(|p| p.to_canonical_string())
            .unwrap_or_else(|e| format!("error:{e:#}"))
    };
    eyre!(
        "UNCONTRACTED_CHECK: kind={} graph={} id={}: the uncontracted value differs from the contracted one through eps^{} (depth {}, first difference {}){}\ncontracted: {}\nuncontracted: {}",
        attempt.kind.as_str(),
        attempt.graph,
        attempt.id,
        attempt.order,
        attempt.depth,
        check.first_difference.as_text(),
        check
            .error
            .as_ref()
            .map(|e| format!(" (comparison error: {e})"))
            .unwrap_or_default(),
        show(&check.contracted),
        show(&ours)
    )
}

/// The error of an attempt whose refusal is not designed
/// ([`Outcome::Error`]): it fails the diagram, in check mode and in mode on
/// (R2-I1).
pub(super) fn error_failure(attempt: &Attempt, mode: Mode) -> eyre::Report {
    let (code, detail) = match &attempt.outcome {
        Outcome::Error(decline) => (decline.code.as_str(), decline.detail.as_str()),
        _ => ("none", ""),
    };
    eyre!(
        "UNCONTRACTED kind={} graph={} id={} mode={}: outcome=error:{code}: the record is refused for a reason that is not a designed refusal (a bug of the emitter or of vakint's entry): {detail}",
        attempt.kind.as_str(),
        attempt.graph,
        attempt.id,
        mode.as_str()
    )
}

// ---------------------------------------------------------------------------
// GAMMALOOP_UNCONTRACTED_STRICT (a no-op) and the per-diagram lines
// ---------------------------------------------------------------------------

/// The environment variable of the former strict switch: strict is the only
/// behaviour of on and check (helper VH round 2 R2-I1); `1` is accepted as a
/// no-op and `0` is refused when the switch initialises
/// ([`super::uncontracted_ledger::parse_strict`]).
pub(crate) const STRICT_VARIABLE: &str = "GAMMALOOP_UNCONTRACTED_STRICT";

/// Tests: run `f` "strict" (the only behaviour; `strict = false` no longer
/// exists, R2-I1).
#[cfg(test)]
pub(crate) fn with_strict<T>(strict: bool, f: impl FnOnce() -> T) -> T {
    assert!(
        strict,
        "the non-strict mode is removed (R2-I1): an error outcome always fails the diagram"
    );
    f()
}

/// vakint's unknown-master refusal of an attempt of `kind`, counted in its
/// diagram's lines.
pub(super) fn tally_unknown_masters(kind: Kind) {
    super::uncontracted_ledger::record_unknown_masters(kind);
}

/// Tests: the `UNCONTRACTED_SUMMARY` lines printed on this thread since the
/// last call (one per diagram ledger, [`super::uncontracted_ledger`]).
#[cfg(test)]
pub(crate) fn take_summaries() -> Vec<String> {
    super::uncontracted_ledger::take_summaries()
}

/// The log lines, the dump and the test capture of one attempt (and its
/// check in check mode).
pub(super) fn report(attempt: &Attempt, check: Option<&Check>) {
    // the occurrence key (U7a delta, the timing note) and the record digest
    let key = super::uncontracted_ledger::occurrence_fields(attempt.text.as_deref());
    let head = format!(
        "UNCONTRACTED kind={} graph={} id={} outcome={} closure={} loops={}",
        attempt.kind.as_str(),
        attempt.graph,
        attempt.id,
        attempt.outcome_name(),
        attempt.closure,
        attempt.loops
    );
    match &attempt.outcome {
        Outcome::Value(_) => eprintln!("{head}{}{key}", attempt.fields),
        Outcome::Zero => eprintln!("{head}{key}"),
        Outcome::Declined(decline) | Outcome::Error(decline) => {
            eprintln!(
                "{head} build_s={:.3} evaluate_s={:.3} detail=\"{}\"{key}",
                attempt.build_seconds,
                attempt.evaluate_seconds,
                escape(&decline.detail)
            )
        }
    }
    super::uncontracted_ledger::record(attempt, check, switch().mode);
    if let Some(check) = check {
        eprintln!(
            "UNCONTRACTED_CHECK kind={} graph={} id={} identical={} equal={} full={} contracted_s={:.3} uncontracted_s={:.3} depth={} valid_through=eps^{} first_difference={}{key}",
            attempt.kind.as_str(),
            attempt.graph,
            attempt.id,
            check.identical,
            check.equal,
            check.full,
            check.contracted_seconds,
            check.uncontracted_seconds,
            attempt.depth,
            attempt.order,
            check.first_difference.as_text()
        );
    }
    dump(attempt, check);
    #[cfg(test)]
    capture(attempt, check);
}

// ---------------------------------------------------------------------------
// GAMMALOOP_UNCONTRACTED_DUMP
// ---------------------------------------------------------------------------

fn dump_directory() -> Option<&'static str> {
    static DIRECTORY: LazyLock<Option<String>> = LazyLock::new(|| {
        std::env::var("GAMMALOOP_UNCONTRACTED_DUMP")
            .ok()
            .filter(|value| !value.is_empty())
    });
    DIRECTORY.as_deref()
}

static SEQUENCE: AtomicUsize = AtomicUsize::new(0);

/// The `u7a-value.v0` text of an attempt.
pub(super) fn value_text(attempt: &Attempt, check: Option<&Check>) -> String {
    let mut lines = vec![
        "schema\tu7a-value.v0".to_owned(),
        format!("kind\t{}", attempt.kind.as_str()),
        format!("graph\t{}", attempt.graph),
        format!("id\t{}", attempt.id),
        format!("depth\t{}", attempt.depth),
        format!("outcome\t{}", attempt.outcome_name()),
        format!("valid_through\teps^{}", attempt.order),
    ];
    if let Some(check) = check {
        lines.push(format!("value\t{}", check.contracted.to_canonical_string()));
        lines.push(format!(
            "contracted\t{}",
            check.contracted.to_canonical_string()
        ));
    } else if let Some(value) = attempt.value() {
        lines.push(format!("value\t{}", value.to_canonical_string()));
    }
    if let Some(value) = attempt.value() {
        lines.push(format!("uncontracted\t{}", value.to_canonical_string()));
    }
    if let Some(check) = check {
        lines.push(format!("identical\t{}", check.identical));
        lines.push(format!("equal\t{}", check.equal));
        lines.push(format!(
            "first_difference\t{}",
            check.first_difference.as_text()
        ));
    }
    if let Outcome::Declined(decline) | Outcome::Error(decline) = &attempt.outcome {
        lines.push(format!("detail\t{}", escape(&decline.detail)));
    }
    lines.push("end".to_owned());
    lines.join("\n") + "\n"
}

/// `GAMMALOOP_UNCONTRACTED_DUMP_MAX=<records>` (default 20000): at most that
/// many attempts are dumped per process, and at most 1 GiB of record text in
/// all (a five-loop record is MBs); the first refused dump prints one
/// `UNCONTRACTED_DUMP cap reached` line (U7a review M5).
fn dump_cap() -> usize {
    static CAP: LazyLock<usize> = LazyLock::new(|| {
        std::env::var("GAMMALOOP_UNCONTRACTED_DUMP_MAX")
            .ok()
            .map(|value| {
                value.trim().parse::<usize>().unwrap_or_else(|_| {
                    panic!("GAMMALOOP_UNCONTRACTED_DUMP_MAX must be a count, not `{value}`")
                })
            })
            .unwrap_or(20_000)
    });
    *CAP
}

const DUMP_MAX_BYTES: usize = 1 << 30;

static DUMPED_BYTES: AtomicUsize = AtomicUsize::new(0);

fn dump(attempt: &Attempt, check: Option<&Check>) {
    let Some(directory) = dump_directory() else {
        return;
    };
    let n = SEQUENCE.fetch_add(1, Ordering::Relaxed);
    let bytes = attempt.text.as_ref().map_or(0, String::len);
    if n >= dump_cap() || DUMPED_BYTES.fetch_add(bytes, Ordering::Relaxed) + bytes > DUMP_MAX_BYTES
    {
        static WARNED: std::sync::Once = std::sync::Once::new();
        WARNED.call_once(|| {
            eprintln!(
                "UNCONTRACTED_DUMP cap reached (records={} max={} max_bytes={DUMP_MAX_BYTES}): later attempts are not dumped",
                n,
                dump_cap()
            )
        });
        return;
    }
    let folder = format!("{directory}/{}", attempt.kind.as_str());
    let _ = std::fs::create_dir_all(&folder);
    let name: String = attempt
        .graph
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '-' || c == '_' {
                c
            } else {
                '_'
            }
        })
        .collect();
    let stem = format!(
        "{folder}/{name}-id{}-{}-{n:06}",
        attempt.id,
        std::process::id()
    );
    if let Some(text) = &attempt.text {
        let _ = std::fs::write(format!("{stem}.record.txt"), text);
    }
    let _ = std::fs::write(format!("{stem}.value.txt"), value_text(attempt, check));
}

/// Tests: one reported attempt.
#[cfg(test)]
#[derive(Clone, Debug)]
pub(crate) struct Captured {
    pub kind: Kind,
    pub graph: String,
    pub id: usize,
    pub outcome: String,
    pub text: Option<String>,
    pub value: Option<Atom>,
    pub contracted: Option<Atom>,
    pub identical: Option<bool>,
    pub equal: Option<bool>,
    pub first_difference: Option<FirstDifference>,
    pub conserved: usize,
    /// The remainder loops, the depth and the order through which the value
    /// is complete.
    pub loops: usize,
    pub depth: usize,
    pub order: i64,
}

#[cfg(test)]
thread_local! {
    static CAPTURE: RefCell<Option<Vec<Captured>>> = const { RefCell::new(None) };
}

#[cfg(test)]
fn capture(attempt: &Attempt, check: Option<&Check>) {
    CAPTURE.with(|slot| {
        if let Some(list) = slot.borrow_mut().as_mut() {
            list.push(Captured {
                kind: attempt.kind,
                graph: attempt.graph.clone(),
                id: attempt.id,
                outcome: attempt.outcome_name(),
                text: attempt.text.clone(),
                value: attempt.value(),
                contracted: check.map(|check| check.contracted.clone()),
                identical: check.map(|check| check.identical),
                equal: check.map(|check| check.equal),
                first_difference: check.map(|check| check.first_difference),
                conserved: attempt.conserved,
                loops: attempt.loops,
                depth: attempt.depth,
                order: attempt.order,
            });
        }
    });
}

/// Tests: run `f` and return every attempt reported on this thread inside it.
#[cfg(test)]
pub(crate) fn with_capture<T>(f: impl FnOnce() -> T) -> (T, Vec<Captured>) {
    let previous = CAPTURE.with(|slot| slot.replace(Some(Vec::new())));
    let out = f();
    let captured = CAPTURE
        .with(|slot| slot.replace(previous))
        .unwrap_or_default();
    (out, captured)
}

/// Why the K3 hook builds no record for a remainder term of an IR hard part
/// whose kind is enabled (today's path takes it, no line is printed): a
/// zero-loop term (Δ(γ) inserted, no vakint call today, OQ-25) or a
/// correction of order `ω = 0` (K8, out of scope).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum NoRecord {
    ZeroLoops,
    OrderZero,
}

#[cfg(test)]
thread_local! {
    static NO_RECORDS: RefCell<Vec<(String, NoRecord)>> = const { RefCell::new(Vec::new()) };
}

/// The K3 hook gives no record for a term of `graph` (tests: noted on this
/// thread; nothing otherwise).
pub(super) fn no_record(graph: &str, why: NoRecord) {
    #[cfg(test)]
    NO_RECORDS.with(|list| list.borrow_mut().push((graph.to_owned(), why)));
    #[cfg(not(test))]
    let _ = (graph, why);
}

/// Tests: the terms the K3 hook gave no record on this thread since the last
/// call.
#[cfg(test)]
pub(crate) fn take_no_records() -> Vec<(String, NoRecord)> {
    NO_RECORDS.with(|list| std::mem::take(&mut *list.borrow_mut()))
}

/// Tests: a change of one record's emitted edges and of the Δ at one of its
/// `ct:` nodes, applied by [`build`] before the Δ momentum-conservation step
/// (the contracted source and today's path untouched): `(kind, node, Δ,
/// edges)`. E.g. a signed mutation of an emitted loop momentum, or a zero by
/// momentum conservation added to a Δ (Astra U7a finding 2).
#[cfg(test)]
pub(crate) type InsertionHook = std::rc::Rc<dyn Fn(Kind, usize, &mut Atom, &mut [FactorizedEdge])>;

#[cfg(test)]
thread_local! {
    static INSERTION_HOOK: RefCell<Option<InsertionHook>> = const { RefCell::new(None) };
}

/// Tests: run `f` with `hook` applied to every insertion of every record
/// built on this thread (restored on unwinding too).
#[cfg(test)]
pub(crate) fn with_insertion_hook<T>(hook: InsertionHook, f: impl FnOnce() -> T) -> T {
    struct Restore(Option<InsertionHook>);
    impl Drop for Restore {
        fn drop(&mut self) {
            let previous = self.0.take();
            INSERTION_HOOK.with(|slot| *slot.borrow_mut() = previous);
        }
    }
    let _restore = Restore(INSERTION_HOOK.with(|slot| slot.borrow_mut().replace(hook)));
    f()
}

#[cfg(test)]
fn insertion_hook(kind: Kind, node: usize, delta: &mut Atom, edges: &mut [FactorizedEdge]) {
    let hook = INSERTION_HOOK.with(|slot| slot.borrow().clone());
    if let Some(hook) = hook {
        hook(kind, node, delta, edges);
    }
}

// ---------------------------------------------------------------------------
// The record (K1)
// ---------------------------------------------------------------------------

/// A record built and read back through vakint's frozen reader.
pub(super) struct Built {
    pub record: FactorizedNumerator,
    pub text: String,
    pub degree: u32,
    pub externals: usize,
    pub insertions: usize,
    /// Insertions rewritten by the momentum-conservation step.
    pub conserved: usize,
}

/// One edge of the diagram, for the line kinds and the closed cycles.
struct EdgeInfo {
    line: LineKind,
    particle: String,
    source: Option<usize>,
    target: Option<usize>,
}

fn particle_line(model: &Model, particle: Option<crate::model::ParticleId>) -> (String, LineKind) {
    let Some(particle) = particle else {
        return ("none".into(), LineKind::Scalar);
    };
    let particle = particle.resolve(model);
    let line = if particle.is_ghost() {
        LineKind::Ghost
    } else if particle.spin == 2 {
        LineKind::Fermion
    } else if particle.spin == 3 {
        LineKind::Vector
    } else {
        LineKind::Scalar
    };
    (particle.name.clone(), line)
}

/// Every edge of the graph: line kind, particle and end vertices (`None` on
/// an external side).
fn edge_table(graph: &Graph, model: &Model) -> BTreeMap<usize, EdgeInfo> {
    let mut table = BTreeMap::new();
    for (pair, edge, data) in graph.iter_edges() {
        let (source, target) = match pair {
            HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => (
                Some(usize::from(graph.node_id(source))),
                Some(usize::from(graph.node_id(sink))),
            ),
            HedgePair::Unpaired { hedge, flow } => {
                let node = Some(usize::from(graph.node_id(hedge)));
                match flow {
                    Flow::Sink => (None, node),
                    Flow::Source => (node, None),
                }
            }
        };
        let (particle, line) = particle_line(model, data.data.particle());
        table.insert(
            usize::from(edge),
            EdgeInfo {
                line,
                particle,
                source,
                target,
            },
        );
    }
    table
}

/// The closed `kind` cycles of `edges` (k2-child-ins mkfixture rule): a
/// component of the `kind` edges of `edges` is closed when every vertex has
/// exactly two `kind` half-edges in it and no other `kind` edge of the graph
/// (an external one included) touches it.
fn closed_cycles(
    table: &BTreeMap<usize, EdgeInfo>,
    edges: &BTreeSet<usize>,
    kind: LineKind,
) -> u32 {
    fn find(parent: &mut BTreeMap<usize, usize>, x: usize) -> usize {
        let mut x = x;
        while let Some(&p) = parent.get(&x) {
            if p == x {
                break;
            }
            x = p;
        }
        parent.entry(x).or_insert(x);
        x
    }
    let mut degree = BTreeMap::<usize, u32>::new();
    let mut parent = BTreeMap::<usize, usize>::new();
    for edge in edges {
        let Some(info) = table.get(edge) else {
            continue;
        };
        if info.line != kind {
            continue;
        }
        let (Some(a), Some(b)) = (info.source, info.target) else {
            continue;
        };
        *degree.entry(a).or_default() += 1;
        *degree.entry(b).or_default() += 1;
        let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
        parent.insert(ra, rb);
    }
    let mut outside = BTreeSet::new();
    for (edge, info) in table {
        if info.line == kind && !edges.contains(edge) {
            outside.extend(info.source);
            outside.extend(info.target);
        }
    }
    let mut components = BTreeMap::<usize, Vec<usize>>::new();
    for vertex in degree.keys().copied().collect::<Vec<_>>() {
        let root = find(&mut parent, vertex);
        components.entry(root).or_default().push(vertex);
    }
    components
        .values()
        .filter(|vertices| {
            vertices
                .iter()
                .all(|v| degree[v] == 2 && !outside.contains(v))
        })
        .count() as u32
}

/// The integer coefficients of a sum of bare `Q(e)` (a routed momentum).
fn q_coefficients(atom: &Atom) -> Option<BTreeMap<usize, i64>> {
    let mut out = BTreeMap::<usize, i64>::new();
    if atom.is_zero() {
        return Some(out);
    }
    let mut term = |view: AtomView<'_>| -> Option<()> {
        let (coefficient, momentum) = match view {
            AtomView::Fun(_) => (1, view),
            AtomView::Mul(product) => {
                let factors = product.iter().collect::<Vec<_>>();
                if factors.len() != 2 {
                    return None;
                }
                let (number, momentum) = if matches!(factors[0], AtomView::Num(_)) {
                    (factors[0], factors[1])
                } else {
                    (factors[1], factors[0])
                };
                (i64::try_from(number).ok()?, momentum)
            }
            _ => return None,
        };
        let AtomView::Fun(function) = momentum else {
            return None;
        };
        if function.get_symbol() != GS.emr_mom || function.get_nargs() != 1 {
            return None;
        }
        let edge = usize::try_from(function.iter().next()?).ok()?;
        *out.entry(edge).or_default() += coefficient;
        Some(())
    };
    match atom.as_view() {
        AtomView::Add(sum) => {
            for view in sum.iter() {
                term(view)?;
            }
        }
        view => term(view)?,
    }
    out.retain(|_, coefficient| *coefficient != 0);
    Some(out)
}

/// `mink(dim)`, the representation of a momentum's slot in a dot.
fn minkowski() -> Atom {
    Minkowski {}.new_rep(GS.dim).to_symbolic([])
}

/// The c-split of the closed-root factor: `c = s · dot(Q(x0),Q(x0))^k`, with
/// `k` the exponent of that dot in `c`, and `s` free of momenta.
pub(super) fn split_root_factor(
    factor: &Atom,
    x0: usize,
) -> std::result::Result<(i64, Atom), Decline> {
    let p2 = function!(
        SPENSO_TAG.dot,
        function!(GS.emr_mom, x0, minkowski()),
        function!(GS.emr_mom, x0, minkowski())
    );
    let factors = match factor.as_view() {
        AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
        view => vec![view],
    };
    let mut k = 0i64;
    for view in factors {
        if view == p2.as_view() {
            k += 1;
        } else if let AtomView::Pow(power) = view {
            let (base, exponent) = power.get_base_exp();
            if base == p2.as_view() {
                k += i64::try_from(exponent).map_err(|_| {
                    Decline::new(
                        "root_normalisation_momentum",
                        format!("a non-integer power of p^2 in {factor}"),
                    )
                })?;
            }
        }
    }
    let scalars = factor * p2.pow(Atom::num(-k));
    if scalars.contains_symbol(GS.emr_mom) {
        return Err(Decline::new(
            "root_normalisation_momentum",
            format!("the closed-root factor {factor} is not a scalar times a power of p^2"),
        ));
    }
    Ok((k, scalars))
}

/// The `mink(.., hedge(h, k))` / `bis(.., hedge(h, k))` slots of `atoms`
/// whose half-edge `h` is in `hedges`, sorted by canonical string.
fn slots_of(atoms: &[&Atom], hedges: &BTreeSet<usize>) -> Vec<Atom> {
    let mink = symbolica::symbol!("spenso::mink");
    let bis = symbolica::symbol!("spenso::bis");
    let mut found = BTreeMap::new();
    for atom in atoms {
        atom.as_view().visitor(&mut |view| {
            if let AtomView::Fun(function) = view
                && (function.get_symbol() == mink || function.get_symbol() == bis)
                && function.get_nargs() == 2
                && let Some(AtomView::Fun(index)) = function.iter().nth(1)
                && index.get_symbol() == GS.hedgeaind
                && let Some(h) = index.iter().next()
                && let Ok(h) = usize::try_from(h)
                && hedges.contains(&h)
            {
                found.insert(view.to_canonical_string(), view.to_owned());
                return false;
            }
            true
        });
    }
    found.into_values().collect()
}

/// Whether `atom` contains `Q(edge, ..)`.
fn mentions_momentum(atom: &Atom, edge: usize) -> bool {
    let mut found = false;
    atom.as_view().visitor(&mut |view| {
        if let AtomView::Fun(function) = view
            && function.get_symbol() == GS.emr_mom
            && let Some(first) = function.iter().next()
            && usize::try_from(first).ok() == Some(edge)
        {
            found = true;
        }
        !found
    });
    found
}

/// The index- and momentum-free factors of one term (the converters'
/// `scalar_part`).
fn scalar_part(term: AtomView<'_>) -> Atom {
    const INDEXED: [&str; 8] = [
        "mink(", "bis(", "::Q(", "::g(", "::gamma(", "coad(", "cof(", "dind(",
    ];
    let factors = match term {
        AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
        view => vec![view],
    };
    let mut scalar = Atom::num(1);
    for factor in factors {
        let text = factor.to_canonical_string();
        if !INDEXED.iter().any(|key| text.contains(key)) {
            scalar *= factor.to_owned();
        }
    }
    scalar
}

/// Whether the canonical text of `atom` is a rational number `-?n(/d)?`.
fn is_rational(atom: &Atom) -> bool {
    let text = atom.to_canonical_string();
    let text = text.strip_prefix('-').unwrap_or(&text);
    let mut parts = text.splitn(2, '/');
    let digits = |part: Option<&str>| {
        part.is_some_and(|p| !p.is_empty() && p.bytes().all(|b| b.is_ascii_digit()))
    };
    let numerator = parts.next();
    match parts.next() {
        None => digits(numerator),
        Some(denominator) => digits(numerator) && digits(Some(denominator)),
    }
}

/// The insertion split `Δ = coefficient × factor` (port of the converters'
/// `split_insertion`): the scalar part of the first expanded term, when every
/// term's scalar part is a rational multiple of it; else `1 × Δ`.
pub(super) fn split_insertion(delta: &Atom) -> (Atom, Atom) {
    let expanded = delta.expand();
    if expanded.is_zero() {
        return (Atom::num(1), delta.clone());
    }
    let terms = match expanded.as_view() {
        AtomView::Add(sum) => sum.iter().collect::<Vec<_>>(),
        view => vec![view],
    };
    let coefficient = scalar_part(terms[0]);
    let mut factor = Atom::Zero;
    for term in &terms {
        let ratio = scalar_part(*term) / &coefficient;
        if !is_rational(&ratio) {
            return (Atom::num(1), delta.clone());
        }
        factor += term.to_owned() / &coefficient;
    }
    (coefficient, factor.expand())
}

/// What the Δ momentum-conservation step did at one `ct:` node
/// ([`conserve`]).
#[derive(Clone, Debug, PartialEq)]
pub(super) enum Conservation {
    /// The Δ with the eliminated momentum rewritten, expanded.
    Rewritten(Atom),
    /// The node's record momenta do not conserve, but only in independent
    /// external components (crown legs a K2 / K3 record keeps independent,
    /// `externals=independent`): no relation holds in the record; the Δ is
    /// unchanged.
    Independent,
    /// No momentum of the node occurs in the Δ: unchanged.
    Absent,
}

/// The Δ momentum-conservation step (U7a-DESIGN section 10): at a `ct:` node
/// whose record momenta conserve, eliminate the largest edge momentum of
/// the node that occurs in `delta` through `Σ σ_e Q(e) = 0` and expand.
///
/// `conserved`: the leading momentum components (loop components first)
/// that the record must conserve at every vertex, as vakint's reader checks
/// them (`validate_closure`: every component of a record with a Taylor step
/// and conserved externals, K1; the loop components of one with independent
/// externals, K2 / K3). A nonzero residual there is a MISROUTED record (the
/// reader would refuse it as `MomentumNotConserved`): `Err`, the routing
/// error `outcome=error:routing` (Astra U7a finding 2), never the designed
/// `insertion_mass_symbol` the caller gives a Δ whose mass symbol stays. A
/// residual in independent external components only is
/// [`Conservation::Independent`].
pub(super) fn conserve(
    delta: &Atom,
    node: usize,
    edges: &[FactorizedEdge],
    conserved: usize,
) -> std::result::Result<Conservation, Decline> {
    let mut net = BTreeMap::<usize, i64>::new();
    for edge in edges {
        let mut sign = 0;
        if edge.target == Endpoint::Vertex(node) {
            sign += 1;
        }
        if edge.source == Endpoint::Vertex(node) {
            sign -= 1;
        }
        if sign != 0 {
            net.insert(edge.id, sign);
        }
    }
    let Some(length) = edges.first().map(|edge| edge.momentum.len()) else {
        return Ok(Conservation::Absent);
    };
    let mut residual = vec![0i64; length];
    for edge in edges {
        if let Some(sign) = net.get(&edge.id) {
            for (slot, c) in residual.iter_mut().zip(&edge.momentum) {
                *slot += sign * c;
            }
        }
    }
    if residual.iter().take(conserved).any(|c| *c != 0) {
        return Err(Decline::new(
            "routing",
            format!(
                "the record momenta at the ct: node {node} do not conserve: residual {residual:?} is nonzero in the first {conserved} components (loop components first), which the record must conserve (vakint's reader: MomentumNotConserved); a misrouted record, not a mass refusal"
            ),
        ));
    }
    if residual.iter().any(|c| *c != 0) {
        return Ok(Conservation::Independent);
    }
    let Some(star) = net
        .keys()
        .rev()
        .find(|edge| mentions_momentum(delta, **edge))
        .copied()
    else {
        return Ok(Conservation::Absent);
    };
    let sign = net[&star];
    let mut image = Atom::Zero;
    for (edge, s) in &net {
        if *edge != star {
            image += Atom::num(-sign * s) * function!(GS.emr_mom, *edge, W_.x___);
        }
    }
    Ok(Conservation::Rewritten(
        delta
            .replace(function!(GS.emr_mom, star, W_.x___))
            .with(image)
            .expand(),
    ))
}

/// `Err` (the routing error) when the record momenta do not conserve at
/// some vertex in their first `conserved` components (loop components
/// first), as vakint's reader checks them ([`conserve`] checks one `ct:`
/// node).
pub(super) fn routing_conserved(
    edges: &[FactorizedEdge],
    conserved: usize,
) -> std::result::Result<(), Decline> {
    let mut residuals = BTreeMap::<usize, Vec<i64>>::new();
    for edge in edges {
        for (end, sign) in [(&edge.target, 1_i64), (&edge.source, -1_i64)] {
            if let Endpoint::Vertex(vertex) = end {
                let residual = residuals
                    .entry(*vertex)
                    .or_insert_with(|| vec![0; edge.momentum.len()]);
                for (slot, c) in residual.iter_mut().zip(&edge.momentum) {
                    *slot += sign * c;
                }
            }
        }
    }
    match residuals
        .iter()
        .find(|(_, residual)| residual.iter().take(conserved).any(|c| *c != 0))
    {
        Some((vertex, residual)) => Err(Decline::new(
            "routing",
            format!(
                "the record momenta at vertex {vertex} do not conserve: residual {residual:?} is nonzero in the first {conserved} components (loop components first), which the record must conserve (vakint's reader: MomentumNotConserved); a misrouted record, not a mass refusal"
            ),
        )),
        None => Ok(()),
    }
}

/// Whether a record's only internal propagator is its closing line
/// `closing`: a one-loop remainder whose loop is the closing line alone
/// (after the Taylor step a massive tadpole, which vakint's closure entry
/// refuses: "every line is massive: not a closure"; U7-PLAN 14(d),
/// B_self_energy_L4_0049 id 8). The emitter declines exactly this shape,
/// `one_loop_tree_cut` (today's route integrates it whole). Every other
/// one-loop remainder (two or more propagators; today's route names every
/// one-loop remainder `one_loop_tree_cut`) is a record: Astra U7a, "keep the
/// rebase test narrow".
pub(super) fn single_propagator(edges: &[FactorizedEdge], closing: usize) -> bool {
    let mut internal = edges.iter().filter(|edge| !edge.is_external());
    matches!((internal.next(), internal.next()), (Some(edge), None) if edge.id == closing)
}

/// The mass symbols of a record whose occurrence in an insertion vakint
/// refuses ("contains a mass symbol"), as vakint's `mass_symbol_in`
/// (factorized.rs): every symbol of an edge's mass, and R*'s reference and
/// auxiliary mass symbols by name (U7a review I5).
fn insertion_masses(edges: &[FactorizedEdge]) -> BTreeSet<String> {
    let mut masses = BTreeSet::new();
    for edge in edges {
        edge.mass_squared.as_view().visitor(&mut |part| {
            if let AtomView::Var(v) = part {
                masses.insert(v.get_symbol().get_name().to_string());
            }
            true
        });
    }
    masses
}

/// Whether `atom` contains a mass symbol of [`insertion_masses`] (or one of
/// R*'s reference / auxiliary mass symbols by name).
fn mentions_mass(atom: &Atom, masses: &BTreeSet<String>) -> bool {
    let mut found = false;
    atom.as_view().visitor(&mut |part| {
        if let AtomView::Var(v) = part {
            let symbol = v.get_symbol();
            let name = symbol.get_name();
            if masses.contains(name)
                || name.ends_with("::rstar_reference_mass_squared")
                || name.ends_with("::rstar_mass_squared")
            {
                found = true;
            }
        }
        !found
    });
    found
}

/// How a record states its externals and its closure (U7a-DESIGN 4.1).
enum Frame<'s> {
    /// K1: one external `p(1) = Q(x0)`; σ folds every root leg onto ±p(1);
    /// the root closure is the projector with the c-split.
    Root {
        closure: &'s Atom,
        closed_root: &'s ClosedRootNormalization,
    },
    /// K2: every crown edge `x_j` is its own independent `p(j)` (sorted);
    /// the projector is the closed child's closure (`p2_power=0`, its
    /// normalisation in `prefactor_scalars`), or none for an open child.
    Child {
        externals: Vec<usize>,
        projector: Option<(Atom, Atom)>,
    },
    /// K3: every crown leg of the hard part `x_j` its own independent `p(j)`
    /// (sorted), `role=root` for the scope's root legs (`roots`); the IR
    /// Dirac closure as projector (`p2_power=0`, its `1/Tr[1]` the caller's)
    /// or none; degree `ω`, fixed mass, the site's depth.
    IrHard {
        externals: Vec<usize>,
        roots: BTreeSet<usize>,
        projector: Option<Atom>,
        degree: usize,
    },
}

/// The K3 frame of one remainder term of an IR hard part: its crown legs
/// (the remainder's soft edges, `region.boundary`) and the IR Dirac closure.
fn ir_hard_frame<'s>(
    site: &Site<'s>,
    remainder: &RStarIntegrand,
) -> std::result::Result<Frame<'s>, Decline> {
    let Some(ir_hard) = site.ir_hard else {
        return Err(Decline::new(
            "ir_hard_site_missing",
            "a K3 record outside the IR hard part hook",
        ));
    };
    let soft = remainder
        .soft_edges()
        .iter()
        .map(|edge| usize::from(*edge))
        .collect::<BTreeSet<_>>();
    let boundary = site
        .graph
        .iter_edges_of(&site.term.remainder.boundary)
        .map(|(_, edge, _)| usize::from(edge))
        .collect::<BTreeSet<_>>();
    if soft != boundary {
        return Err(Decline::new(
            "crown",
            format!("the soft edges {soft:?} are not the hard part's crown {boundary:?}"),
        ));
    }
    if soft.is_empty() {
        return Err(Decline::new("crown", "a hard part without crown legs"));
    }
    let projector = match ir_hard.closure {
        None => None,
        Some(boundary) => {
            if boundary.closure.contains_symbol(GS.emr_mom) {
                return Err(Decline::new(
                    "ir_dirac_closure",
                    format!(
                        "the IR Dirac closure {} carries a momentum",
                        boundary.closure
                    ),
                ));
            }
            Some(boundary.closure.clone())
        }
    };
    Ok(Frame::IrHard {
        externals: soft.into_iter().collect(),
        roots: site
            .root_edges
            .iter()
            .map(|edge| usize::from(*edge))
            .collect(),
        projector,
        degree: ir_hard.degree,
    })
}

/// The K2 frame of a child remainder: its crown and its closure.
fn child_frame<'s>(
    site: &Site<'s>,
    remainder: &RStarIntegrand,
) -> std::result::Result<Frame<'s>, Decline> {
    if site.top_mass.is_some_and(|mode| mode != TopMassTerms::Keep) {
        return Err(Decline::new(
            "top_mass_mode",
            format!("top-mass terms {:?}", site.top_mass),
        ));
    }
    let externals = site
        .root_edges
        .iter()
        .map(|edge| usize::from(*edge))
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    let soft = remainder
        .soft_edges()
        .iter()
        .map(|edge| usize::from(*edge))
        .collect::<BTreeSet<_>>();
    if soft != externals.iter().copied().collect::<BTreeSet<_>>() {
        return Err(Decline::new(
            "crown",
            format!(
                "the crown {externals:?} is not the remainder's soft edges {:?}",
                soft
            ),
        ));
    }
    if externals.is_empty() {
        return Err(Decline::new("crown", "a child without crown edges"));
    }
    let projector = match (site.root_closure, site.root_projector, site.closure_data) {
        (None, false, _) => None,
        (Some(closure), true, Some(data)) => {
            if let Some(closed_root) = site.closed_root
                && closed_root.factor != data.normalization
                && !(&closed_root.factor - &data.normalization)
                    .expand()
                    .is_zero()
            {
                return Err(Decline::new(
                    "normalisation_mismatch",
                    format!(
                        "closed child factor {} against the closure normalisation {}",
                        closed_root.factor, data.normalization
                    ),
                ));
            }
            Some((closure.clone(), data.normalization.clone()))
        }
        (closure, projector, _) => {
            return Err(Decline::new(
                "closure_data_missing",
                format!(
                    "a child closure {} with a projector {projector} and no closure data",
                    closure.map_or_else(|| "none".to_owned(), |c| c.to_canonical_string())
                ),
            ));
        }
    };
    Ok(Frame::Child {
        externals,
        projector,
    })
}

/// The record id offset of the INCOMING leg of a crown edge with both ends
/// in the scope (K2 / K3, `HedgePair::Paired`): the outgoing leg keeps the
/// edge's id `id`, the incoming one is `PAIRED_LEG_OFFSET + id` (the
/// VERIFIED wh-fixtures/k2-child states legs 4 and 100004). Unique only
/// while every edge id of the graph stays below the offset
/// ([`paired_leg_id`]; helper VH round 2, R2-M4).
pub(super) const PAIRED_LEG_OFFSET: usize = 100_000;

/// The record id of the incoming leg of the paired crown edge `id` in a
/// graph whose largest edge id is `last_edge`; the `crown` error when an
/// edge id reaches [`PAIRED_LEG_OFFSET`] (the ids could collide).
pub(super) fn paired_leg_id(id: usize, last_edge: usize) -> std::result::Result<usize, Decline> {
    if last_edge >= PAIRED_LEG_OFFSET {
        return Err(Decline::new(
            "crown",
            format!(
                "the graph's edge {last_edge} reaches the paired-leg offset {PAIRED_LEG_OFFSET}: the incoming leg of crown edge {id} could collide with an edge id"
            ),
        ));
    }
    Ok(PAIRED_LEG_OFFSET + id)
}

/// The record of one K1 root or K2 child remainder (`Ok(None)`: a zero
/// insertion).
fn build(
    site: &Site<'_>,
    remainder: &RStarIntegrand,
    kind: Kind,
    id: usize,
    mode: Mode,
) -> std::result::Result<Option<Built>, Decline> {
    let graph = site.graph;
    let region = &site.term.remainder;
    if site.uv_settings.vakint.additional_normalization.trim() != "1" {
        return Err(Decline::new(
            "additional_normalization",
            site.uv_settings.vakint.additional_normalization.clone(),
        ));
    }
    let frame = match kind {
        Kind::Root => {
            let (Some(closure), true) = (site.root_closure, site.root_projector) else {
                return Err(Decline::new("open_root", "no root closure or projector"));
            };
            let Some(closed_root) = site.closed_root else {
                return Err(Decline::new(
                    "root_normalisation_opaque",
                    "a root projector without the closed-root normalisation as data",
                ));
            };
            if site.top_mass != Some(TopMassTerms::Keep) {
                return Err(Decline::new(
                    "top_mass_mode",
                    format!("top-mass terms {:?}", site.top_mass),
                ));
            }
            Frame::Root {
                closure,
                closed_root,
            }
        }
        Kind::Child => child_frame(site, remainder)?,
        Kind::IrHard => ir_hard_frame(site, remainder)?,
        Kind::Surrogate => k5::frame(site, remainder)?,
        other => return Err(Decline::new("kind_not_built", other.as_str())),
    };
    let externals = match &frame {
        Frame::Root { closed_root, .. } => vec![usize::from(closed_root.external)],
        Frame::Child { externals, .. } | Frame::IrHard { externals, .. } => externals.clone(),
    };
    let basis = remainder.basis();
    let loop_edges = basis
        .loop_edges
        .iter()
        .map(|edge| usize::from(*edge))
        .collect::<Vec<_>>();
    if loop_edges.len() != region.loops {
        return Err(Decline::new(
            "routing",
            format!(
                "{} loop edges for {} remainder loops",
                loop_edges.len(),
                region.loops
            ),
        ));
    }
    let loops = loop_edges.len();
    let width = loops + externals.len();
    let loop_index = loop_edges
        .iter()
        .enumerate()
        .map(|(i, edge)| (*edge, i))
        .collect::<BTreeMap<_, _>>();
    let marked = site.plan.marked_edge;
    // The external slot of a crown momentum Q(e): K1 σ (every leg ->
    // ±Q(x0) = ±p(1)); K2 the crown edge's own p(j).
    let slot_of = |edge: usize| -> std::result::Result<(usize, i64), Decline> {
        match &frame {
            Frame::Root { closed_root, .. } => {
                let x0 = function!(GS.emr_mom, usize::from(closed_root.external));
                let image =
                    function!(GS.emr_mom, edge).replace_multiple(&closed_root.external_routing);
                if image == x0 {
                    Ok((loops, 1))
                } else if image == -x0.clone() {
                    Ok((loops, -1))
                } else {
                    Err(Decline::new(
                        "sigma",
                        format!(
                            "sigma maps Q({edge}) to {image}, not ±Q({})",
                            usize::from(closed_root.external)
                        ),
                    ))
                }
            }
            Frame::Child { externals, .. } | Frame::IrHard { externals, .. } => externals
                .iter()
                .position(|x| *x == edge)
                .map(|j| (loops + j, 1))
                .ok_or_else(|| {
                    Decline::new(
                        "routing",
                        format!("Q({edge}) is not a crown momentum {externals:?}"),
                    )
                }),
        }
    };
    let independent = matches!(frame, Frame::Child { .. } | Frame::IrHard { .. });
    // the momentum components every vertex of the record conserves, as
    // vakint's reader checks them (`validate_closure`): all of them with
    // conserved externals (K1), the loop components with independent ones
    let conserved_components = if independent { loops } else { width };
    let vector = |edge: EdgeIndex| -> std::result::Result<Vec<i64>, Decline> {
        let mut out = vec![0i64; width];
        let loop_part = basis.loop_atom::<Atom>(edge, GS.emr_mom, &[], true);
        let ext_part = basis.ext_atom::<Atom>(edge, GS.emr_mom, &[], true);
        let parse = |part: &Atom| {
            q_coefficients(part).ok_or_else(|| {
                Decline::new(
                    "routing",
                    format!("edge {}: momentum {part}", usize::from(edge)),
                )
            })
        };
        for (e, c) in parse(&loop_part)? {
            let Some(i) = loop_index.get(&e) else {
                return Err(Decline::new(
                    "routing",
                    format!("edge {}: Q({e}) is not a loop edge", usize::from(edge)),
                ));
            };
            out[*i] += c;
        }
        for (e, c) in parse(&ext_part)? {
            let (slot, sign) = slot_of(e)?;
            out[slot] += c * sign;
        }
        Ok(out)
    };
    let table = edge_table(graph, site.model);
    let last_edge = table.keys().next_back().copied().unwrap_or(0);
    // contracted supports -> ct: nodes (mkfixture-k1 rule)
    let mut contracted: SuBitGraph = graph.empty_subgraph();
    for slot in &region.insertions {
        contracted.union_with(&slot.support);
    }
    let mut paired = Vec::new();
    for (pair, edge, data) in graph.iter_edges_of(&region.remaining) {
        let HedgePair::Paired { source, sink } = pair else {
            continue;
        };
        paired.push((edge, source, sink, data.data));
    }
    let mut node_index = BTreeMap::new();
    for (_, source, sink, _) in &paired {
        for hedge in [*source, *sink] {
            let node = graph.node_id(hedge);
            node_index.insert(usize::from(node), node);
        }
    }
    let endpoints = node_index.keys().copied().collect::<BTreeSet<_>>();
    let mut merge = BTreeMap::<usize, usize>::new();
    let mut nodes = Vec::new();
    for (index, slot) in region.insertions.iter().enumerate() {
        let vertices = graph
            .iter_nodes_of(&slot.support)
            .map(|(node, _, _)| usize::from(node))
            .collect::<BTreeSet<_>>();
        let Some(node) = vertices
            .iter()
            .find(|v| endpoints.contains(v))
            .or_else(|| vertices.first())
            .copied()
        else {
            return Err(Decline::new("missing_child", "a support without vertices"));
        };
        for vertex in &vertices {
            if merge.insert(*vertex, node).is_some() {
                return Err(Decline::new(
                    "insertions_share_a_vertex",
                    format!("vertex {vertex}"),
                ));
            }
        }
        nodes.push((node, index));
    }
    let image = |vertex: usize| *merge.get(&vertex).unwrap_or(&vertex);
    // legs
    let mut edges = Vec::new();
    let mut ends = Vec::<(usize, usize)>::new();
    let mut legs = BTreeSet::new();
    for (pair, edge, _) in graph.iter_edges_of(&region.boundary) {
        let id = usize::from(edge);
        let info = &table[&id];
        let (slot, sign) = slot_of(id)?;
        let mut momentum = vec![0i64; width];
        momentum[slot] = sign;
        // (record id, the crown half-edge, whether the leg enters the scope)
        let pieces = match pair {
            HedgePair::Unpaired { hedge, flow } => vec![(id, hedge, flow == Flow::Sink)],
            HedgePair::Split {
                source,
                sink,
                split,
            } if independent => match split {
                Flow::Sink => vec![(id, sink, true)],
                Flow::Source => vec![(id, source, false)],
            },
            HedgePair::Paired { source, sink } if independent => {
                vec![
                    (id, source, false),
                    (paired_leg_id(id, last_edge)?, sink, true),
                ]
            }
            _ => {
                return Err(Decline::new(
                    "root_crown_paired",
                    format!("crown edge {id} has both ends in the root"),
                ));
            }
        };
        for (record_id, hedge, inward) in pieces {
            let h = usize::from(hedge);
            let raw_node = graph.node_id(hedge);
            let raw = usize::from(raw_node);
            node_index.insert(raw, raw_node);
            let v = image(raw);
            let (source, target, source_hedge, target_hedge) = if inward {
                (Endpoint::External, Endpoint::Vertex(v), None, Some(h))
            } else {
                (Endpoint::Vertex(v), Endpoint::External, Some(h), None)
            };
            edges.push(FactorizedEdge {
                id: record_id,
                source,
                target,
                source_hedge,
                target_hedge,
                particle: info.particle.clone(),
                line: info.line,
                power: 1,
                mass_squared: Atom::Zero,
                momentum: momentum.clone(),
                factor: Atom::num(1),
            });
            ends.push((v, h));
            legs.insert(h);
        }
    }
    // internal edges
    let reference = Atom::var(site.reference_mass_squared);
    let mut closing = false;
    for (edge, source, sink, data) in &paired {
        let id = usize::from(*edge);
        let (s, t) = (
            usize::from(graph.node_id(*source)),
            usize::from(graph.node_id(*sink)),
        );
        let (a, b) = (image(s), image(t));
        let (hs, ht) = (usize::from(*source), usize::from(*sink));
        let info = &table[&id];
        let mass_squared = if *edge == marked {
            closing = true;
            reference.clone()
        } else {
            Atom::Zero
        };
        edges.push(FactorizedEdge {
            id,
            source: Endpoint::Vertex(a),
            target: Endpoint::Vertex(b),
            source_hedge: Some(hs),
            target_hedge: Some(ht),
            particle: info.particle.clone(),
            line: info.line,
            power: 1,
            mass_squared,
            momentum: vector(*edge)?,
            factor: data.num.value.map_mink_dim(GS.dim),
        });
        ends.push((a, hs));
        ends.push((b, ht));
    }
    if !closing {
        return Err(Decline::new(
            "no_closing_line",
            format!(
                "the marked edge {} is not a remaining propagator",
                usize::from(marked)
            ),
        ));
    }
    edges.sort_by_key(|edge| edge.id);
    // vertices
    let mut vertices = Vec::new();
    for (raw, node) in &node_index {
        if merge.contains_key(raw) {
            continue;
        }
        let vertex = &graph[*node];
        let rule = vertex.vertex_rule.map_or_else(
            || vertex.name.value.clone(),
            |rule| rule.resolve(site.model).name.clone(),
        );
        vertices.push(FactorizedVertex {
            id: *raw,
            rule,
            factor: vertex.get_num().map_mink_dim(GS.dim),
        });
    }
    // insertions
    let masses = insertion_masses(&edges);
    let mut insertions = Vec::new();
    let mut conserved = 0;
    for (node, index) in &nodes {
        let slot = &region.insertions[*index];
        let Some(delta) = site.children.get(&slot.support) else {
            return Err(Decline::new(
                "missing_child",
                super::plan::join_edge_ids(&super::plan::uv_child_edge_ids(graph, &slot.support)),
            ));
        };
        if delta.is_zero() {
            return Ok(None);
        }
        let mut delta = delta.clone();
        #[cfg(test)]
        insertion_hook(kind, *node, &mut delta, &mut edges);
        if mentions_mass(&delta, &masses) {
            // a misrouted node is an error (`?`), never the mass refusal
            match conserve(&delta, *node, &edges, conserved_components)? {
                Conservation::Rewritten(rewritten) if !mentions_mass(&rewritten, &masses) => {
                    if rewritten.is_zero() {
                        return Ok(None);
                    }
                    delta = rewritten;
                    conserved += 1;
                }
                _ => {
                    // a designed decline only for a record whose momenta
                    // conserve at EVERY vertex: a misrouted edge elsewhere is
                    // the routing error vakint's reader would raise (review
                    // of the U7a delta, refusal-integrity finding 3)
                    routing_conserved(&edges, conserved_components)?;
                    return Err(Decline::new(
                        "insertion_mass_symbol",
                        format!(
                            "the counterterm of {} contains a mass symbol",
                            super::plan::join_edge_ids(&super::plan::uv_child_edge_ids(
                                graph,
                                &slot.support
                            ))
                        ),
                    ));
                }
            }
        }
        let support_edges = super::plan::uv_child_edge_ids(graph, &slot.support);
        let support_set = support_edges.iter().copied().collect::<BTreeSet<_>>();
        let mut hedges = ends
            .iter()
            .filter(|(v, _)| v == node)
            .map(|(_, h)| *h)
            .collect::<Vec<_>>();
        hedges.sort_unstable();
        let (coefficient, factor) = split_insertion(&delta);
        insertions.push(Insertion {
            vertex: *node,
            hedges,
            cycles: (
                closed_cycles(&table, &support_set, LineKind::Fermion),
                closed_cycles(&table, &support_set, LineKind::Ghost),
            ),
            coefficient,
            factor,
        });
        vertices.push(FactorizedVertex {
            id: *node,
            rule: format!(
                "ct:edges-{}",
                support_edges
                    .iter()
                    .map(usize::to_string)
                    .collect::<Vec<_>>()
                    .join("-")
            ),
            factor: Atom::num(1),
        });
    }
    vertices.sort_by_key(|vertex| vertex.id);
    insertions.sort_by_key(|insertion| insertion.vertex);
    // PRODUCT (check mode): the record's factors are the site's numerator
    if mode == Mode::Check {
        let product = edges
            .iter()
            .map(|edge| &edge.factor)
            .chain(vertices.iter().map(|vertex| &vertex.factor))
            .fold(Atom::num(1), |product, factor| product * factor);
        let source = graph
            .numerator(&region.remaining, &contracted)
            .to_d_dim(GS.dim)
            .get_single_atom()
            .map_err(|error| Decline::new("product_mismatch", format!("{error}")))?;
        if product != source && !(&product - &source).expand().is_zero() {
            return Err(Decline::new(
                "product_mismatch",
                "the record's factors are not the site's numerator",
            ));
        }
    }
    // K5 (U7c): coefficient vertices -> `ins:`, bubble lines, cycles
    let k5_cycles = if kind == Kind::Surrogate {
        Some(k5::finish(
            graph,
            &table,
            region,
            &merge,
            &nodes,
            &mut edges,
            &mut vertices,
            &mut insertions,
        )?)
    } else {
        None
    };
    // ports, projector, prefactor
    let factor_atoms = edges
        .iter()
        .map(|edge| &edge.factor)
        .chain(vertices.iter().map(|vertex| &vertex.factor))
        .chain(insertions.iter().map(|insertion| &insertion.factor))
        .collect::<Vec<_>>();
    let ports = slots_of(&factor_atoms, &legs);
    let (projector, scalars) = match &frame {
        Frame::Root {
            closure,
            closed_root,
        } => {
            let (p2_power, scalars) =
                split_root_factor(&closed_root.factor, usize::from(closed_root.external))?;
            (
                Some(Projector {
                    tensor: (*closure).clone(),
                    p2_power,
                }),
                scalars,
            )
        }
        Frame::Child {
            projector: Some((tensor, normalization)),
            ..
        } => (
            Some(Projector {
                tensor: tensor.clone(),
                p2_power: 0,
            }),
            normalization.clone(),
        ),
        Frame::Child {
            projector: None, ..
        } => (None, Atom::num(1)),
        // the IR Dirac closure's 1/Tr[1] is the caller's (ir_dirac::reassemble)
        Frame::IrHard { projector, .. } => (
            projector.as_ref().map(|tensor| Projector {
                tensor: tensor.clone(),
                p2_power: 0,
            }),
            Atom::num(1),
        ),
    };
    let internal = table
        .iter()
        .filter(|(_, info)| info.source.is_some() && info.target.is_some())
        .map(|(edge, _)| *edge)
        .collect::<BTreeSet<_>>();
    let original_cycles = (
        closed_cycles(&table, &internal, LineKind::Fermion),
        closed_cycles(&table, &internal, LineKind::Ghost),
    );
    let original_cycles = k5_cycles.unwrap_or(original_cycles);
    let degree = match &frame {
        Frame::IrHard { degree, .. } => u32::try_from(*degree)
            .map_err(|_| Decline::new("taylor", format!("hard Taylor order {degree}")))?,
        _ => u32::try_from(site.plan.taylor_order).map_err(|_| {
            Decline::new("taylor", format!("Taylor order {}", site.plan.taylor_order))
        })?,
    };
    let contracted_label = if nodes.is_empty() {
        "none".to_owned()
    } else {
        nodes
            .iter()
            .map(|(_, index)| {
                super::plan::join_edge_ids(&super::plan::uv_child_edge_ids(
                    graph,
                    &region.insertions[*index].support,
                ))
            })
            .collect::<Vec<_>>()
            .join(";")
    };
    let scope = match site.plan.root_boundary {
        Some(_) => "child",
        None => "root",
    };
    // the evaluation path: a child record names the forest of its parents
    // (helper VH round 2 R2-M3)
    let path = super::uncontracted_ledger::current_path();
    let label = match site.ir_hard {
        Some(ir_hard) if kind == Kind::IrHard => format!(
            "u7b kind=ir-hard graph={} scope={scope}:{} corr={} term={} attempt={} hard={} remaining={} contracted={contracted_label} marked_edge={} omega={} depth={} closure={} path={path} id={id}",
            graph.name,
            super::plan::join_edge_ids(&super::plan::uv_child_edge_ids(graph, &site.plan.support)),
            ir_hard.correction,
            ir_hard.term,
            ir_hard.attempt,
            ir_hard.hard,
            super::plan::join_edge_ids(&super::plan::uv_child_edge_ids(graph, &region.remaining)),
            usize::from(marked),
            ir_hard.degree,
            ir_hard.depth,
            closure_name(site, kind)
        ),
        _ => format!(
            "u7a kind={} graph={} scope={}:{} term={} remaining={} contracted={contracted_label} marked_edge={} closure={} path={path} id={id}",
            kind.as_str(),
            graph.name,
            kind.as_str(),
            super::plan::join_edge_ids(&super::plan::uv_child_edge_ids(graph, &region.support)),
            if region.insertions.is_empty() {
                "bare".to_owned()
            } else {
                format!("remainder_loops={}", region.loops)
            },
            super::plan::join_edge_ids(&super::plan::uv_child_edge_ids(graph, &region.remaining)),
            usize::from(marked),
            closure_name(site, kind)
        ),
    };
    // K2: every crown leg is a root leg of the child; K3: a root leg of the
    // scope (`role=root`) or a soft line of it (`role=soft`)
    let softs = match &frame {
        Frame::Root { .. } => Vec::new(),
        Frame::Child { .. } => externals
            .iter()
            .enumerate()
            .map(|(j, edge)| SoftLeg {
                external: j + 1,
                edge: *edge,
                root: true,
            })
            .collect(),
        Frame::IrHard { roots, .. } => externals
            .iter()
            .enumerate()
            .map(|(j, edge)| SoftLeg {
                external: j + 1,
                edge: *edge,
                root: roots.contains(edge),
            })
            .collect(),
    };
    // K1 / K2: the reference mass compensated (ChildMass::Reference); K3: an
    // IR hard part is integrated at the fixed reference mass, every order
    // kept, the root's compensation applied by the caller afterwards
    let compensate_mass = !matches!(frame, Frame::IrHard { .. });
    let n_insertions = insertions.len();
    let n_externals = externals.len();
    let record = FactorizedNumerator {
        label: label.replace(['\t', '\n', '\r'], " "),
        kind: FactorizedKind::Closure(ClosureSpec {
            edge: usize::from(marked),
            mass_squared: reference.clone(),
            auxiliary_mass_squared: site.auxiliary_mass_squared.clone(),
            power: 1,
            independent_externals: independent,
            taylor: Some(TaylorSpec {
                degree,
                scaled_edges: externals.clone(),
                compensate_mass,
                orders: TaylorOrders::All,
            }),
            softs,
            original_cycles: Some(original_cycles),
            depth: Some(site.depth()),
            weights: Vec::new(),
        }),
        dimension: GS.dim,
        loop_edges: loop_edges.iter().map(|edge| Some(*edge)).collect(),
        external_edges: externals.iter().map(|edge| Some(*edge)).collect(),
        edges,
        vertices,
        ports,
        projector,
        prefactor: Prefactor {
            overall: Atom::num(1),
            loop_sign: 1,
            nf_power: 0,
            scalars,
            colour: Atom::num(1),
        },
        colour_branch: None,
        insertions,
        record: None,
    };
    let text = record
        .to_text()
        .map_err(|error| Decline::new("text", error.to_string()))?;
    let record = FactorizedNumerator::from_text(&text).map_err(|error| Decline::vakint(&error))?;
    // last, so that every error of the record above stays loud
    if single_propagator(&record.edges, usize::from(marked)) {
        return Err(Decline::new(
            "one_loop_tree_cut",
            format!(
                "the remainder {}: its only propagator is the closing line {} (after the Taylor step a massive tadpole; today's route integrates it whole)",
                super::plan::join_edge_ids(&super::plan::uv_child_edge_ids(
                    graph,
                    &region.remaining
                )),
                usize::from(marked)
            ),
        ));
    }
    Ok(Some(Built {
        record,
        text,
        degree,
        externals: n_externals,
        insertions: n_insertions,
        conserved,
    }))
}
