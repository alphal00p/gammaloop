//! Explicit hard Taylor coefficients with the existing UV tensor conventions.

use std::{
    borrow::Cow,
    collections::BTreeSet,
    sync::{LazyLock, OnceLock},
    time::Instant,
};

use eyre::{Result, WrapErr, ensure, eyre};
use feynkit_tensor::TensorReducer;
use idenso::{
    CookMode, CookSettings,
    color::ColorSimplifySettings,
    dirac::GammaSimplifySettings,
    tensor::{AlgebraSettings, ContractSettings, SymbolicTensor},
};
use linnet::half_edge::{
    involution::{EdgeIndex, HedgePair},
    subgraph::{InternalSubGraph, SuBitGraph, SubSetOps},
};
use spenso::{
    network::{library::symbolic::ETS, tags::SPENSO_TAG},
    structure::{
        partial::PartialStructureExt,
        representation::{Minkowski, RepName},
        slot::{IsAbstractSlot, SlotMatch, SlotMatcher},
    },
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    function,
    id::Replacement,
    symbol,
};
use vakint::{Vakint, VakintSettings, stage_timing::StageSpan, symbols::S};

use crate::{
    debug_tags,
    graph::{Autogen, Graph, LMBext, LoopMomentumBasis},
    model::Model,
    utils::{GS, W_},
    uv::{
        UVgenerationSettings, UltravioletGraph,
        approx::{
            counts,
            integrated::{Integrated, to_vakint_integrand},
            streamed,
        },
        uv_graph::UVE,
    },
};

use super::{GraphRegion, TaylorExpansion, integrate_single_mass, project::project_hard_vacuum};

pub(super) type RootProjector<'a> = dyn Fn(&Atom) -> Result<Atom> + 'a;

// The shared transport knows explicit edge momenta and these FeynKit basis
// coordinates. Other live momentum conventions must keep their own occurrence.
static UNSUPPORTED_PHYSICAL_MOMENTA: LazyLock<Vec<Symbol>> = LazyLock::new(|| {
    let supported = [
        GS.emr_mom,
        feynkit_graph::symbols::loop_momentum(),
        feynkit_graph::symbols::external_momentum(),
    ];
    [GS.loop_mom, GS.external_mom, S.k, S.p]
        .into_iter()
        .filter(|head| !supported.contains(head))
        .collect()
});

static MEMORY_RELEASE: std::sync::OnceLock<fn()> = std::sync::OnceLock::new();

/// Set by the host binary to return freed allocator memory to the kernel at
/// once; used before the hard Taylor fallback (`simplify_hard_taylor`).
pub fn set_memory_release(release: fn()) {
    let _ = MEMORY_RELEASE.set(release);
}

/// Exclusive wall-clock stages within one outer integral. Backend time may
/// contain its own finer timing events; those must not be added to this total.
pub(super) struct IntegrationTimings {
    started: Instant,
    last: Instant,
    routing: f64,
    taylor: f64,
    materialization: f64,
    tensor_projection: f64,
    algebra: f64,
    conversion: f64,
    backend: f64,
    finish: f64,
}

impl IntegrationTimings {
    pub(super) fn new() -> Self {
        let started = Instant::now();
        Self {
            started,
            last: started,
            routing: 0.0,
            taylor: 0.0,
            materialization: 0.0,
            tensor_projection: 0.0,
            algebra: 0.0,
            conversion: 0.0,
            backend: 0.0,
            finish: 0.0,
        }
    }

    fn lap(&mut self) -> f64 {
        let now = Instant::now();
        let elapsed = now.duration_since(self.last).as_secs_f64();
        self.last = now;
        elapsed
    }

    /// The timings of one outer integration on the remainder-boundary path
    /// (kind `single_mass`): taylor = den + collect, conversion = frames +
    /// cut + routing + contract + project + assemble + verify, backend =
    /// evaluate, finish = finish.
    pub(super) fn boundary(
        started: Instant,
        taylor: f64,
        conversion: f64,
        backend: f64,
        finish: f64,
    ) -> Self {
        Self {
            started,
            last: Instant::now(),
            routing: 0.0,
            taylor,
            materialization: 0.0,
            tensor_projection: 0.0,
            algebra: 0.0,
            conversion,
            backend,
            finish,
        }
    }

    pub(super) fn log(&self, graph: &Graph, integration_id: usize, loops: usize, kind: &str) {
        tracing::info!(
            uv = true,
            rstar = true,
            timing = true,
            stage = "rstar_outer_integration_detail",
            graph = %graph.name,
            integration_id,
            loops,
            kind,
            total_seconds = self.started.elapsed().as_secs_f64(),
            routing_seconds = self.routing,
            taylor_seconds = self.taylor,
            materialization_seconds = self.materialization,
            tensor_projection_seconds = self.tensor_projection,
            algebra_seconds = self.algebra,
            conversion_seconds = self.conversion,
            backend_seconds = self.backend,
            finish_seconds = self.finish,
            "Outer integration stage timings"
        );
        if vakint::census::enabled() {
            vakint::census::record(
                "outer_integration",
                serde_json::json!({
                    "kind": kind,
                    "loops": loops,
                    "total_seconds": self.started.elapsed().as_secs_f64(),
                    "routing_seconds": self.routing,
                    "taylor_seconds": self.taylor,
                    "tensor_projection_seconds": self.tensor_projection,
                    "algebra_seconds": self.algebra,
                    "conversion_seconds": self.conversion,
                    "backend_seconds": self.backend,
                    "finish_seconds": self.finish,
                }),
            );
        }
    }
}

/// `GAMMALOOP_EDGE_CONTRACTION=<list>` (`hard`, `remainder`, `all`): contract
/// those R* expressions with edge momenta independent
/// ([`RStarIntegrand::simplify_edge_basis`]).
/// The edge-contraction checks skip inputs above
/// `GAMMALOOP_EDGE_CHECK_MAX_BYTES` (those the loop basis cannot finish).
fn edge_check_admits(expression: &Atom) -> bool {
    static MAX: LazyLock<usize> = LazyLock::new(|| {
        std::env::var("GAMMALOOP_EDGE_CHECK_MAX_BYTES")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(usize::MAX)
    });
    let admitted = expression.as_view().get_byte_size() <= *MAX;
    if !admitted {
        eprintln!(
            "STREAM_CHECK stage=edge_skipped bytes={}",
            expression.as_view().get_byte_size()
        );
    }
    admitted
}

/// `GAMMALOOP_EDGE_DEBUG_DUMP=<dir>`: write the input of a failing step.
fn edge_debug_dump(stage: &str, expression: &Atom) {
    if let Ok(dir) = std::env::var("GAMMALOOP_EDGE_DEBUG_DUMP") {
        let _ = std::fs::write(
            format!("{dir}/{stage}-{}.txt", std::process::id()),
            format!("{expression}"),
        );
    }
}

thread_local! {
    static DRESSED_SCOPE: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Run `f` as the R* evaluation of a dressed skeleton (an insertion surrogate):
/// `GAMMALOOP_EDGE_CONTRACTION=auto` keeps edge momenta in its ordinary
/// remainders too (measured: fewer terms there, more in insertion-free graphs).
pub(crate) fn in_dressed_scope<T>(f: impl FnOnce() -> T) -> T {
    let previous = DRESSED_SCOPE.with(|scope| scope.replace(true));
    let result = f();
    DRESSED_SCOPE.with(|scope| scope.set(previous));
    result
}

thread_local! {
    static UV_CHILD_SCOPE: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Run `f` as the R* evaluation of a proper UV child (its own complete-graph
/// R*, nested children included) whose Taylor order is `taylor_order`. See
/// [`uv_child_edge_remainder`]. A nested child sets its own scope and the
/// enclosing one is restored afterwards.
pub(crate) fn in_uv_child_scope<T>(taylor_order: i32, f: impl FnOnce() -> T) -> T {
    // Restored on unwinding too: beta catches a panicking diagram and the worker
    // thread goes on, so a scope left set would treat later root remainders as
    // children.
    struct Restore(bool);
    impl Drop for Restore {
        fn drop(&mut self) {
            let previous = self.0;
            UV_CHILD_SCOPE.with(|scope| scope.set(previous));
        }
    }
    let edge_basis = child_edge_basis(*UV_CHILD_EDGE_REMAINDER, taylor_order);
    let _restore = Restore(UV_CHILD_SCOPE.with(|scope| scope.replace(edge_basis)));
    f()
}

/// Whether a child of `taylor_order` contracts its ordinary remainders in the
/// edge basis under `mode`.
fn child_edge_basis(mode: ChildEdgeRemainder, taylor_order: i32) -> bool {
    match mode {
        ChildEdgeRemainder::Off => false,
        ChildEdgeRemainder::DegreeZero => taylor_order <= 0,
        ChildEdgeRemainder::All => true,
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ChildEdgeRemainder {
    Off,
    /// `1`: only logarithmically divergent (Taylor order 0) children. Their
    /// open-slot remainders collapse in the edge basis; children of degree >= 1
    /// need the boundary-momentum dependence in the outer Taylor step, where
    /// the edge-kept form is slower (FIX1-VALIDATION.md section 5).
    DegreeZero,
    /// `all`: every child, as before.
    All,
}

static UV_CHILD_EDGE_REMAINDER: LazyLock<ChildEdgeRemainder> =
    LazyLock::new(|| {
        match std::env::var("GAMMALOOP_UV_CHILD_EDGE_REMAINDER")
            .as_deref()
            .map(str::trim)
        {
            Ok("1") => ChildEdgeRemainder::DegreeZero,
            Ok("all") => ChildEdgeRemainder::All,
            _ => ChildEdgeRemainder::Off,
        }
    });

#[cfg(test)]
#[test]
fn child_edge_basis_follows_the_mode_and_the_taylor_order() {
    use ChildEdgeRemainder::*;
    for order in [-1, 0, 1, 2] {
        assert!(!child_edge_basis(Off, order));
        assert!(child_edge_basis(All, order));
        assert_eq!(child_edge_basis(DegreeZero, order), order <= 0);
    }
}

#[cfg(test)]
#[test]
fn uv_child_scope_is_restored_when_the_child_panics() {
    // The mode is read once per process; check against whatever it is here.
    let expected = |order| child_edge_basis(*UV_CHILD_EDGE_REMAINDER, order);
    for enclosing in [false, true] {
        for order in [0, 1] {
            UV_CHILD_SCOPE.with(|scope| scope.set(enclosing));
            let seen = std::cell::Cell::new(None);
            let unwound = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                in_uv_child_scope(order, || {
                    seen.set(Some(UV_CHILD_SCOPE.with(|scope| scope.get())));
                    panic!("child integration failed");
                })
            }));
            assert!(unwound.is_err());
            assert_eq!(seen.get(), Some(expected(order)));
            // Restored to the enclosing value, which differs from the child's in
            // at least one of the four cases whatever the mode is.
            assert_eq!(UV_CHILD_SCOPE.with(|scope| scope.get()), enclosing);
            // Normal return: the child sees its own value, the enclosing one comes back.
            let inside = in_uv_child_scope(order, || UV_CHILD_SCOPE.with(|scope| scope.get()));
            assert_eq!(inside, expected(order));
            assert_eq!(UV_CHILD_SCOPE.with(|scope| scope.get()), enclosing);
            // A nested child of the other degree sets and restores its own value.
            let nested = in_uv_child_scope(order, || {
                let inner =
                    in_uv_child_scope(1 - order, || UV_CHILD_SCOPE.with(|scope| scope.get()));
                (inner, UV_CHILD_SCOPE.with(|scope| scope.get()))
            });
            assert_eq!(nested, (expected(1 - order), expected(order)));
        }
    }
    UV_CHILD_SCOPE.with(|scope| scope.set(false));
}

/// `GAMMALOOP_UV_CHILD_EDGE_REMAINDER=1|all` (opt-in): inside an R* child of
/// Taylor order 0 (`1`) or inside every R* child (`all`), the ordinary
/// remainders are contracted in the edge basis, exactly as
/// `GAMMALOOP_EDGE_CONTRACTION=remainder` does everywhere. A child keeps its
/// boundary slots open, and the edge-basis contraction of those open-index
/// remainders is much cheaper; the closed root keeps the loop basis, where
/// the edge basis is slower. Any other value or unset: unchanged.
fn uv_child_edge_remainder() -> bool {
    UV_CHILD_SCOPE.with(|scope| scope.get())
}

/// `GAMMALOOP_UV_STREAM=off|on|auto` (default off): how the edge-basis
/// remainder of a degree-0 UV child (`GAMMALOOP_UV_CHILD_EDGE_REMAINDER=1`) is
/// integrated.
///
/// - `off`: contract the whole remainder, substitute the edge combinations
///   back, integrate the result (today's path);
/// - `on`: stream it ([`RStarIntegrand::stream_contracted`]): split the input
///   product into pieces, then contract, substitute and integrate one piece at
///   a time and add up only the integrated results;
/// - `auto`: stream only when the input product distributes into at least
///   `GAMMALOOP_UV_STREAM_MIN_COMBOS` (default 729 = 3^6) term combinations
///   (the product of the term counts of its sums); otherwise the whole path.
///
/// `GAMMALOOP_UV_STREAM_PIECES` (default 81) is the target number of pieces.
/// Contracted pieces are added up into batches of at least
/// `GAMMALOOP_UV_STREAM_BATCH_BYTES` (default 256 MiB) of emitted expression,
/// and each batch is integrated as one: a remainder that emits less than one
/// batch is integrated once, as on the whole path; a large one in bounded
/// batches. (A first-piece extrapolation is not robust: a first piece can
/// contract to zero while the remainder emits hundreds of MB.)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum UvStream {
    Off,
    On,
    Auto,
}

pub(super) fn uv_stream() -> UvStream {
    static MODE: LazyLock<UvStream> = LazyLock::new(|| {
        match std::env::var("GAMMALOOP_UV_STREAM")
            .as_deref()
            .map(str::trim)
        {
            Ok("on") | Ok("1") => UvStream::On,
            Ok("auto") => UvStream::Auto,
            _ => UvStream::Off,
        }
    });
    *MODE
}

fn uv_stream_pieces() -> usize {
    static PIECES: LazyLock<usize> = LazyLock::new(|| {
        std::env::var("GAMMALOOP_UV_STREAM_PIECES")
            .ok()
            .and_then(|value| value.parse().ok())
            .filter(|&pieces: &usize| pieces >= 2)
            .unwrap_or(81)
    });
    *PIECES
}

fn uv_stream_batch_bytes() -> usize {
    static BYTES: LazyLock<usize> = LazyLock::new(|| {
        std::env::var("GAMMALOOP_UV_STREAM_BATCH_BYTES")
            .ok()
            .and_then(|value| value.parse().ok())
            .unwrap_or(256 << 20)
    });
    *BYTES
}

fn uv_stream_min_combos() -> usize {
    static COMBOS: LazyLock<usize> = LazyLock::new(|| {
        std::env::var("GAMMALOOP_UV_STREAM_MIN_COMBOS")
            .ok()
            .and_then(|value| value.parse().ok())
            .unwrap_or(729)
    });
    *COMBOS
}

/// Whether hard components keep edge momenta.
fn keeps_hard() -> bool {
    edge_contraction("keep") || edge_contraction("keep_hard") || edge_contraction("auto")
}

/// Whether ordinary remainders keep edge momenta.
fn keeps_remainder() -> bool {
    edge_contraction("keep")
        || edge_contraction("keep_remainder")
        || (edge_contraction("auto") && DRESSED_SCOPE.with(|scope| scope.get()))
}

/// `GAMMALOOP_EDGE_PROJECT_LOG=1`: one stderr line before and after each hard projection.
fn edge_project_log() -> bool {
    static LOG: LazyLock<bool> =
        LazyLock::new(|| std::env::var_os("GAMMALOOP_EDGE_PROJECT_LOG").is_some());
    *LOG
}

/// Largest projected hard term (bytes) rewritten in the loop basis by `back`
/// (`GAMMALOOP_EDGE_BACK_MAX_BYTES`, default 64 MB).
fn edge_back_max_bytes() -> usize {
    static MAX: LazyLock<usize> = LazyLock::new(|| {
        std::env::var("GAMMALOOP_EDGE_BACK_MAX_BYTES")
            .ok()
            .and_then(|value| value.parse().ok())
            .unwrap_or(64_000_000)
    });
    *MAX
}

/// Largest contracted hard Taylor term (bytes) that `pre` rewrites in the loop
/// basis before the projection (`GAMMALOOP_EDGE_PRE_MAX_BYTES`, default 16 MB).
fn edge_pre_max_bytes() -> usize {
    static MAX: LazyLock<usize> = LazyLock::new(|| {
        std::env::var("GAMMALOOP_EDGE_PRE_MAX_BYTES")
            .ok()
            .and_then(|value| value.parse().ok())
            .unwrap_or(16_000_000)
    });
    *MAX
}

/// RSS growth (bytes) allowed to a loop-basis hard Taylor contraction under
/// `fallback` (`GAMMALOOP_EDGE_FALLBACK_GIB`, default 8).
fn edge_fallback_bytes() -> u64 {
    static BYTES: LazyLock<u64> = LazyLock::new(|| {
        let gib: f64 = std::env::var("GAMMALOOP_EDGE_FALLBACK_GIB")
            .ok()
            .and_then(|value| value.parse().ok())
            .unwrap_or(8.0);
        (gib * (1u64 << 30) as f64) as u64
    });
    *BYTES
}

/// Resident set size of this process (Linux `/proc/self/statm`).
fn process_rss_bytes() -> Option<u64> {
    let statm = std::fs::read_to_string("/proc/self/statm").ok()?;
    let pages: u64 = statm.split_whitespace().nth(1)?.parse().ok()?;
    Some(pages * 4096)
}

/// First id of a kept edge-momentum combination `Q(EDGE_OFFSET + n, ...)`.
const EDGE_OFFSET: i64 = 1_000_000;

/// With `flatten` (`IDENSO_EDGE_OUTPUT=flatten`, read by the callers): the
/// expanded form of `view` when it is a nested linear combination of momenta
/// (numeric multiples and sums of vectors `vector` accepts, not all of them at
/// the top level). Routing leaves
/// vertex momenta such as `Q(2,μ) - (Q(2,μ) + Q(3,μ) - Q(6,μ))`; renamed as
/// they are, only the inner sum becomes one vector and the contraction
/// distributes the outer one (nf0 32's hottest UV-child remainder: 4.03 GB
/// instead of 9.0 MB). Expanded, the whole combination becomes one vector.
/// Without `flatten`, `None`.
fn flat_combination(
    view: AtomView<'_>,
    flatten: bool,
    vector: &dyn Fn(AtomView<'_>) -> bool,
) -> Option<Atom> {
    fn linear(view: AtomView<'_>, vector: &dyn Fn(AtomView<'_>) -> bool) -> bool {
        match view {
            AtomView::Add(sum) => sum.iter().all(|term| linear(term, vector)),
            // Exactly one non-numeric factor, checked before descending into it.
            AtomView::Mul(product) => {
                let mut factors = product.iter().filter(|f| !matches!(f, AtomView::Num(_)));
                match (factors.next(), factors.next()) {
                    (Some(factor), None) => linear(factor, vector),
                    _ => false,
                }
            }
            other => vector(other),
        }
    }
    if !flatten {
        return None;
    }
    let AtomView::Add(sum) = view else {
        return None;
    };
    (sum.iter().any(|term| !vector(term)) && linear(view, vector)).then(|| view.expand())
}

/// First id of a renamed momentum combination `Q(EDGE_BASIS_OFFSET + n, ...)`
/// of [`RStarIntegrand::simplify_edge_basis`].
const EDGE_BASIS_OFFSET: i64 = 2_000_000;

/// The rename of [`RStarIntegrand::simplify_edge_basis`]: every linear
/// combination of momenta in one index slot, `Σ c_e Q(e, slot)` with numeric
/// `c_e`, becomes one vector `Q(EDGE_BASIS_OFFSET + n, slot)`, `n` its
/// position in `combinations`. With `flatten`, a nested combination is
/// expanded first ([`flat_combination`]).
fn rename_edge_basis(
    expression: &Atom,
    flatten: bool,
    combinations: &mut std::collections::BTreeMap<Vec<(i64, Atom)>, i64>,
) -> Atom {
    let head = GS.emr_mom;
    let term_of = |term: AtomView<'_>| -> Option<(i64, Atom, Vec<Atom>)> {
        let (coefficient, vector) = match term {
            AtomView::Mul(product) => {
                let factors = product.iter().collect::<Vec<_>>();
                match factors.as_slice() {
                    [n @ AtomView::Num(_), v] | [v, n @ AtomView::Num(_)] => (n.to_owned(), *v),
                    _ => return None,
                }
            }
            other => (Atom::num(1), other),
        };
        let AtomView::Fun(function) = vector else {
            return None;
        };
        if function.get_symbol() != head || function.get_nargs() < 1 {
            return None;
        }
        let edge = i64::try_from(function.get(0)).ok()?;
        if edge >= EDGE_BASIS_OFFSET {
            return None;
        }
        Some((
            edge,
            coefficient,
            function.iter().skip(1).map(|a| a.to_owned()).collect(),
        ))
    };
    expression.replace_map(|view, _, out| {
        let flat = flat_combination(view, flatten, &|leaf| term_of(leaf).is_some());
        let AtomView::Add(sum) = flat.as_ref().map_or(view, Atom::as_view) else {
            if let Some(flat) = flat {
                **out = flat;
            }
            return;
        };
        let mut tail = None;
        let mut terms = Vec::new();
        for term in sum.iter() {
            let Some((edge, coefficient, this_tail)) = term_of(term) else {
                return;
            };
            match &tail {
                None => tail = Some(this_tail),
                Some(t) if *t == this_tail => {}
                _ => return,
            }
            terms.push((edge, coefficient));
        }
        terms.sort_by(|a, b| a.0.cmp(&b.0));
        let next = EDGE_BASIS_OFFSET + combinations.len() as i64;
        let id = *combinations.entry(terms).or_insert(next);
        **out = symbolica::atom::FunctionBuilder::new(head)
            .add_arg(Atom::num(id))
            .add_args(&tail.unwrap_or_default())
            .finish();
    })
}

/// The rename of [`RStarIntegrand::contract_keeping_edges`]: the part of a
/// linear combination of momenta in one index slot that lies on `loop_edges`
/// (at least two of them) becomes one vector `Q(EDGE_OFFSET + n, slot)`, `n`
/// its position in `combinations`; other momenta stay explicit and
/// denominators are left as they are. With `flatten`, a nested combination is
/// expanded first ([`flat_combination`]).
fn rename_kept_edges(
    expression: &Atom,
    loop_edges: &BTreeSet<i64>,
    flatten: bool,
    combinations: &mut std::collections::BTreeMap<Vec<(i64, Atom)>, i64>,
) -> Atom {
    let head = GS.emr_mom;
    let term_of = |term: AtomView<'_>| -> Option<(i64, Atom, Vec<Atom>)> {
        let (coefficient, vector) = match term {
            AtomView::Mul(product) => {
                let factors = product.iter().collect::<Vec<_>>();
                match factors.as_slice() {
                    [n @ AtomView::Num(_), v] | [v, n @ AtomView::Num(_)] => (n.to_owned(), *v),
                    _ => return None,
                }
            }
            other => (Atom::num(1), other),
        };
        let AtomView::Fun(function) = vector else {
            return None;
        };
        if function.get_symbol() != head || function.get_nargs() < 1 {
            return None;
        }
        let edge = i64::try_from(function.get(0)).ok()?;
        Some((
            edge,
            coefficient,
            function.iter().skip(1).map(|a| a.to_owned()).collect(),
        ))
    };
    expression.replace_map(|view, _, out| {
        // Denominators keep their loop-basis momenta: their topology and
        // Taylor arguments are read from them.
        if matches!(view, AtomView::Fun(f) if f.get_symbol() == GS.den) {
            **out = view.to_owned();
            return;
        }
        let flat = flat_combination(view, flatten, &|leaf| term_of(leaf).is_some());
        let AtomView::Add(sum) = flat.as_ref().map_or(view, Atom::as_view) else {
            if let Some(flat) = flat {
                **out = flat;
            }
            return;
        };
        let mut tail = None;
        let mut kept = Vec::new();
        let mut others = Vec::new();
        for term in sum.iter() {
            let Some((edge, coefficient, this_tail)) = term_of(term) else {
                return;
            };
            match &tail {
                None => tail = Some(this_tail),
                Some(t) if *t == this_tail => {}
                _ => return,
            }
            if loop_edges.contains(&edge) {
                kept.push((edge, coefficient));
            } else {
                others.push(term);
            }
        }
        if kept.len() < 2 || tail.as_ref().is_none_or(|t| t.is_empty()) {
            return;
        }
        kept.sort_by(|a, b| a.0.cmp(&b.0));
        let next = EDGE_OFFSET + combinations.len() as i64;
        let id = *combinations.entry(kept).or_insert(next);
        let vector = symbolica::atom::FunctionBuilder::new(head)
            .add_arg(Atom::num(id))
            .add_args(&tail.unwrap_or_default())
            .finish();
        **out = if others.is_empty() {
            vector
        } else {
            Atom::add_many(
                others
                    .iter()
                    .copied()
                    .chain(std::iter::once(vector.as_view())),
            )
        };
    })
}

/// `GAMMALOOP_EDGE_CONTRACTION`: comma-separated modes. Unset means
/// `fallback,back`: hard components are contracted in the loop basis as
/// before and only a contraction that exhausts its memory budget is redone
/// and projected with edge momenta kept (a projected result of at most
/// `GAMMALOOP_EDGE_BACK_MAX_BYTES` rewritten in the loop basis); `off` (or any
/// value naming no mode) restores the loop-basis contraction everywhere.
/// `keep_hard` keeps edge momenta in every hard component.
fn edge_contraction(name: &str) -> bool {
    static LIST: LazyLock<Vec<String>> = LazyLock::new(|| {
        std::env::var("GAMMALOOP_EDGE_CONTRACTION")
            .map(|value| value.split(',').map(|v| v.trim().to_owned()).collect())
            .unwrap_or_else(|_| vec!["fallback".to_owned(), "back".to_owned()])
    });
    LIST.iter().any(|entry| entry == name || entry == "all")
}

/// `GAMMALOOP_HARD_TAYLOR_CHECK=1` also runs the contracted hard Taylor path on
/// every hard component and reports `HARD_TAYLOR_CHECK ... equal=...` on stderr.
fn hard_taylor_check() -> bool {
    static CHECK: LazyLock<bool> =
        LazyLock::new(|| std::env::var_os("GAMMALOOP_HARD_TAYLOR_CHECK").is_some());
    *CHECK
}

/// An opt-in skip of a whole-expression tensor simplification whose input is a
/// closed scalar ([`closed_dot_scalar`]); see [`closed_skip_mode`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ClosedSkipMode {
    /// Unset (or any value but `1`/`check`): the pass always runs, unchanged.
    Off,
    /// `1`: a closed scalar input is used as it is; anything else is simplified.
    Skip,
    /// `check`: the pass runs; a closed input is compared with its result
    /// and the outcome reported on stderr (`CLOSED_SKIP ... identical= equal=`).
    Check,
}

/// The [`ClosedSkipMode`] of the environment variable `name` (`1`, `check`).
pub fn closed_skip_mode(name: &str) -> ClosedSkipMode {
    match std::env::var(name).as_deref() {
        Ok("1") => ClosedSkipMode::Skip,
        Ok("check") => ClosedSkipMode::Check,
        _ => ClosedSkipMode::Off,
    }
}

/// `GAMMALOOP_CLOSED_SKIP_SUMMAND_SIMPLIFY` (`1` or `check`): the skip of the
/// whole-expression [`RStarIntegrand::simplify_summands`] pass after the outer
/// Taylor step of a root-projected integral whose remainder idenso certified a
/// closed scalar when it contracted it (the root-closed path: root metric
/// contracted in the root scope, remainder already contracted to dots). The
/// certificate is a cached producer fact
/// ([`SymbolicTensor::is_certified_closed_scalar`], read before the uncook in
/// [`RStarIntegrand::simplify_view_certified`]); the Taylor step introduces no
/// index (`den` is differentiated in its scalar fourth argument), so it
/// carries over to the Taylor polynomial. `check` runs the pass anyway and
/// compares, also against the syntactic [`closed_dot_scalar`] test. Default off.
fn closed_summand_skip() -> ClosedSkipMode {
    static MODE: LazyLock<ClosedSkipMode> =
        LazyLock::new(|| closed_skip_mode("GAMMALOOP_CLOSED_SKIP_SUMMAND_SIMPLIFY"));
    *MODE
}

/// Whether `expression` is a closed scalar in dot notation, with nothing left
/// for the tensor algebra (no open port, no dummy pair, no gamma, metric,
/// colour, Levi-Civita, trace or chain): every representation occurs as the
/// index-free payload of a rank-one tensor inside a `dot` (`Q(e, mink(d))`),
/// and every other function (`den`, an index-free `Q(e)`, ...) has closed
/// arguments. A conservative syntactic test of one pass over the expression;
/// any other form, including scalar-wrapper markers, counts as not closed.
pub fn closed_dot_scalar(expression: AtomView<'_>) -> bool {
    fn marker(symbol: Symbol) -> bool {
        symbol == SPENSO_TAG.chain
            || symbol == SPENSO_TAG.trace
            || symbol == SPENSO_TAG.pure_scalar
            || symbol == SPENSO_TAG.scalar
            || symbol == SPENSO_TAG.dot
    }
    fn representation(symbol: Symbol) -> bool {
        symbol.has_tag(&SPENSO_TAG.representation)
    }
    // An operand of a dot: a rank-one tensor (or a linear combination of
    // them) whose only representation is an index-free payload.
    fn operand(view: AtomView<'_>) -> bool {
        match view {
            AtomView::Num(_) => true,
            AtomView::Var(variable) => !representation(variable.get_symbol()),
            AtomView::Add(sum) => sum.iter().all(operand),
            AtomView::Mul(product) => product.iter().all(operand),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                operand(base) && operand(exponent)
            }
            AtomView::Fun(function) => {
                let symbol = function.get_symbol();
                if representation(symbol) {
                    function.get_nargs() <= 1 && function.iter().all(closed)
                } else {
                    !marker(symbol) && function.iter().all(operand)
                }
            }
        }
    }
    fn closed(view: AtomView<'_>) -> bool {
        match view {
            AtomView::Num(_) => true,
            AtomView::Var(variable) => !representation(variable.get_symbol()),
            AtomView::Add(sum) => sum.iter().all(closed),
            AtomView::Mul(product) => product.iter().all(closed),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                closed(base) && closed(exponent)
            }
            AtomView::Fun(function) => {
                let symbol = function.get_symbol();
                if symbol == SPENSO_TAG.dot {
                    function.get_nargs() == 2 && function.iter().all(operand)
                } else {
                    !representation(symbol) && !marker(symbol) && function.iter().all(closed)
                }
            }
        }
    }
    closed(expression)
}

/// The original tensor slots and loop carriers survive all graph contractions.
/// Only propagator endpoints are identified in the converter's graph.
pub(super) struct RStarIntegrand {
    pub region: GraphRegion,
    /// The routed product of the remaining numerator, insertions and
    /// denominators, before any index contraction. A hard component is only
    /// ever Taylor projected, and its truncation keeps a tiny fraction of the
    /// contracted terms: it is contracted after the Taylor step.
    routed: Atom,
    /// [`Self::routed`] contracted to dots, built on first use, and whether
    /// idenso certified that contraction a closed scalar (a cached producer
    /// fact, see [`Self::simplify_view_certified`]; only the loop-basis
    /// contraction reports it, the edge-keeping contractions report `false`).
    contracted: OnceLock<(Atom, bool)>,
    /// Edge momenta kept as independent vectors (`GAMMALOOP_EDGE_CONTRACTION=keep`):
    /// each distinct loop-momentum combination `Σ c_e Q(e)` -> its id.
    combinations: std::sync::Mutex<std::collections::BTreeMap<Vec<(i64, Atom)>, i64>>,
    graph: Graph,
    basis: LoopMomentumBasis,
    soft_edges: Vec<EdgeIndex>,
}

impl RStarIntegrand {
    pub fn new(
        graph: &Graph,
        model: &Model,
        region: &GraphRegion,
        insertions: &[Atom],
        marked: EdgeIndex,
        reference_mass_squared: Atom,
    ) -> Result<Self> {
        Self::new_with_closure(
            graph,
            model,
            region,
            insertions,
            marked,
            reference_mass_squared,
            None,
        )
    }

    /// The metric that closes the two open external Lorentz slots of the
    /// graph numerator over `support`, spelled as in every R* remainder (the
    /// d-dimensional continuation of the vertex slots).
    ///
    /// It is momentum-free, so it commutes with routing, the outer Taylor
    /// operation and the reference mass. Applied in a root-scope remainder,
    /// together with its insertions, it contracts the root indices before the
    /// remainder algebra and inside the root `K` only, exactly where the
    /// metric root projector contracted them; child and quotient
    /// counterterms keep their open indices.
    pub fn root_closure(graph: &Graph, support: &SuBitGraph) -> Result<Atom> {
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let numerator = graph
            .numerator(support, &graph.empty_subgraph())
            .to_d_dim(GS.dim)
            .get_single_atom()?;
        let tensor = SymbolicTensor::infer(cooking.try_cook(numerator.as_view())?)?;
        let slots = tensor
            .structure()
            .logical_slots()
            .into_iter()
            .map(|slot| match slot.aind {
                spenso::structure::partial::PartialIndex::Explicit(index) => {
                    Ok(slot
                        .rep()
                        .slot::<spenso::structure::abstract_index::AbstractIndex, _>(index))
                }
                _ => Err(eyre!(
                    "the root closure needs explicitly named external slots"
                )),
            })
            .collect::<Result<Vec<_>>>()?;
        ensure!(
            slots.len() == 2 && slots[0].rep() == slots[1].rep(),
            "the root closure needs exactly two open external slots of one representation, found {}",
            slots.len()
        );
        Ok(cooking.uncook(ETS.metric(slots[0].to_atom(), slots[1].to_atom()).as_view()))
    }

    /// The three metric pairings `g12 g34`, `g13 g24`, `g14 g23` of the four
    /// open slots of the numerator over `support`, when they are exactly four
    /// Minkowski slots of one representation; `None` otherwise. They span the
    /// momentum-independent rank-four tensors in `d` dimensions, the tensor
    /// structure of a logarithmically divergent four-point counterterm.
    pub fn isotropic_boundary(graph: &Graph, support: &SuBitGraph) -> Result<Option<[Atom; 3]>> {
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let numerator = graph
            .numerator(support, &graph.empty_subgraph())
            .to_d_dim(GS.dim)
            .get_single_atom()?;
        // The metric pairings span only the parity-even rank-four tensors; a
        // gamma5 or Levi-Civita symbol could leave an epsilon part they miss.
        let parity_odd = [
            idenso::dirac::AGS.gamma5,
            idenso::dirac::AGS.projp,
            idenso::dirac::AGS.projm,
            *idenso::epsilon::EPSILON_SYMBOL,
        ];
        if parity_odd
            .iter()
            .any(|symbol| numerator.contains_symbol(*symbol))
        {
            return Ok(None);
        }
        let tensor = SymbolicTensor::infer(cooking.try_cook(numerator.as_view())?)?;
        let mut slots = Vec::new();
        for slot in tensor.structure().logical_slots() {
            let spenso::structure::partial::PartialIndex::Explicit(index) = slot.aind else {
                return Ok(None);
            };
            slots.push(
                slot.rep()
                    .slot::<spenso::structure::abstract_index::AbstractIndex, _>(index),
            );
        }
        let minkowski: spenso::structure::representation::LibraryRep = Minkowski {}.into();
        if slots.len() != 4
            || slots
                .iter()
                .any(|slot| slot.rep() != slots[0].rep() || slot.rep().rep != minkowski)
        {
            return Ok(None);
        }
        let metric = |a: usize, b: usize| ETS.metric(slots[a].to_atom(), slots[b].to_atom());
        Ok(Some(
            [
                metric(0, 1) * metric(2, 3),
                metric(0, 2) * metric(1, 3),
                metric(0, 3) * metric(1, 2),
            ]
            .map(|pairing| cooking.uncook(pairing.as_view())),
        ))
    }

    /// As [`Self::new`], with an optional root closure multiplied into the
    /// remainder together with its insertions (see [`Self::root_closure`]).
    pub fn new_with_closure(
        graph: &Graph,
        model: &Model,
        region: &GraphRegion,
        insertions: &[Atom],
        marked: EdgeIndex,
        reference_mass_squared: Atom,
        closure: Option<&Atom>,
    ) -> Result<Self> {
        let _stage = StageSpan::new("gamma.remainder.construct");
        ensure!(
            insertions.len() == region.insertions.len(),
            "each R* contraction requires an insertion polynomial"
        );
        let mut contracted: SuBitGraph = graph.empty_subgraph();
        for slot in &region.insertions {
            contracted.union_with(&slot.support);
        }
        ensure!(
            region.support.subtract(&contracted) == region.remaining,
            "R* remainder does not match its contracted supports"
        );
        let numerator_stage = StageSpan::new("gamma.remainder.numerator");
        let mut expression = graph
            .numerator(&region.remaining, &contracted)
            .to_d_dim(GS.dim)
            .get_single_atom()?;
        drop(numerator_stage);
        let routing_stage = StageSpan::new("gamma.remainder.routing_and_denominators");
        for insertion in insertions {
            expression *= insertion;
        }
        if let Some(closure) = closure {
            expression *= closure;
        }
        // Diagnostics: the remainder before routing (edge momenta intact),
        // numerator times insertions; its routed form is `hard_routed`.
        streamed::dump("hard_source", expression.as_view());
        let contracted = InternalSubGraph::try_new(contracted, &graph.underlying)
            .ok_or_else(|| eyre!("an R* insertion contains an unpaired half-edge"))?;
        let basis = if region.loops == 0 {
            graph.empty_lmb()
        } else {
            graph
                .shrunken_sub_lmb(&region.support, &contracted, region.boundary.clone())
                .wrap_err("could not route an R* prepared component")?
        };
        let soft_edges = graph
            .iter_edges_of(&region.boundary)
            .map(|(_, edge, _)| edge)
            .collect();
        let experiment_input =
            (super::edge_experiment::enabled() && region.loops != 0).then(|| expression.clone());
        if region.loops != 0 {
            // Boundary Q labels are the independent source variables of the
            // local polynomial. A subgraph basis supplies signatures only for
            // its internal edges; routing its boundary would erase them.
            expression = expression.replace_multiple(graph.normal_emr_replacement(
                &region.remaining,
                &basis,
                &[W_.x___],
                |_| true,
            ));
        }
        let minkowski = Minkowski {}.new_rep(GS.dim).to_symbolic([]);
        for (pair, edge, data) in graph.iter_edges_of(&region.remaining) {
            if !pair.is_paired() {
                continue;
            }
            ensure!(
                data.data.mass_atom(model).is_zero(),
                "R* reference preparation currently requires massless physical propagators"
            );
            let hard = basis.loop_atom::<Atom>(edge, GS.emr_mom, &[], true);
            ensure!(
                !hard.is_zero(),
                "an R* remaining propagator has no hard loop momentum"
            );
            let momentum =
                basis.loop_atom(edge, GS.emr_mom, std::slice::from_ref(&minkowski), true)
                    + basis.ext_atom(edge, GS.emr_mom, std::slice::from_ref(&minkowski), true);
            let mass_squared = if edge == marked {
                reference_mass_squared.clone()
            } else {
                Atom::Zero
            };
            // `den` differentiates only its fourth argument. The third is
            // converter metadata and is the mass itself, squared by Vakint's
            // existing graph bridge after the Taylor operation.
            expression /= function!(
                GS.den,
                usize::from(edge),
                hard,
                mass_squared.sqrt(),
                function!(SPENSO_TAG.dot, &momentum, &momentum) - &mass_squared
            );
        }
        let mut converter_graph = graph.clone();
        for slot in &region.insertions {
            let nodes = converter_graph
                .iter_nodes_of(&slot.support)
                .map(|(node, _, _)| node)
                .collect::<Vec<_>>();
            let first = nodes
                .first()
                .ok_or_else(|| eyre!("an R* contracted support contains no vertices"))?;
            let vertex = converter_graph.underlying[*first].clone();
            converter_graph.underlying.identify_nodes(&nodes, vertex);
        }
        converter_graph.underlying.forget_identification_history();
        drop(routing_stage);
        if let Some(input) = experiment_input {
            let routing = graph
                .normal_emr_replacement(&region.remaining, &basis, &[W_.x___], |_| true)
                .into_iter()
                .collect::<Vec<_>>();
            super::edge_experiment::compare(
                graph,
                &converter_graph,
                &region.remaining,
                &input,
                &routing,
                region.loops,
                &|view| Self::simplify_view(view),
            );
        }
        Ok(Self {
            region: region.clone(),
            routed: expression,
            contracted: OnceLock::new(),
            combinations: Default::default(),
            graph: converter_graph,
            basis,
            soft_edges,
        })
    }

    pub fn soft_edges(&self) -> &[EdgeIndex] {
        &self.soft_edges
    }

    /// Whether the routed remainder vanishes (nothing is integrated).
    pub(super) fn routed_is_zero(&self) -> bool {
        self.routed.is_zero()
    }

    /// The routed numerator `N` (numerator × insertions × closure) and the
    /// denominator factors `den(e, ..)^-1`, split off the top-level factors
    /// of the routed product on demand (remainder boundary, DESIGN 1): the
    /// product is exactly that by construction. Exactly one denominator per
    /// paired remaining edge, none inside `N`; `Err(reason)` otherwise.
    pub(super) fn split_routed(&self) -> std::result::Result<(Atom, Vec<Atom>), String> {
        let factors = match self.routed.as_view() {
            AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
            single => vec![single],
        };
        let paired = self
            .graph
            .iter_edges_of(&self.region.remaining)
            .filter(|(pair, _, _)| pair.is_paired())
            .map(|(_, edge, _)| usize::from(edge))
            .collect::<BTreeSet<_>>();
        let mut seen = BTreeSet::new();
        let mut denominators = Vec::new();
        let mut numerator = Vec::new();
        for factor in factors {
            if let AtomView::Pow(power) = factor
                && let (AtomView::Fun(base), exponent) = power.get_base_exp()
                && base.get_symbol() == GS.den
            {
                if i64::try_from(exponent).ok() != Some(-1) {
                    return Err("den_split: a denominator power other than -1".into());
                }
                let Ok(edge) = usize::try_from(base.get(0)) else {
                    return Err("den_split: a denominator without an edge id".into());
                };
                if !paired.contains(&edge) || !seen.insert(edge) {
                    return Err(format!("den_split: edge {edge}"));
                }
                denominators.push(factor.to_owned());
                continue;
            }
            if factor.contains_symbol(GS.den) {
                return Err("den_split: a denominator inside the numerator".into());
            }
            numerator.push(factor);
        }
        if seen != paired {
            return Err("den_split: not one denominator per paired edge".into());
        }
        Ok((Atom::mul_many(numerator.iter().copied()), denominators))
    }

    /// `P = C(σ(N))`: the routed numerator `N` of [`Self::split_routed`] with
    /// the external identification `σ` applied, contracted to dots in the
    /// loop basis as [`Self::expression`] contracts the whole remainder, and
    /// idenso's closed-scalar certificate of this contraction (DESIGN
    /// Addendum A(1): read on the contraction's own result, never by a scan).
    pub(super) fn numerator_contracted(
        &self,
        numerator: Atom,
        external_routing: &[Replacement],
    ) -> Result<(Atom, bool)> {
        let started = Instant::now();
        let routed = numerator.replace_multiple(external_routing);
        drop(numerator);
        let contracted = {
            let _algebra = StageSpan::new("gamma.remainder.algebra");
            Self::simplify_view_certified(routed.as_view())
                .wrap_err("while contracting the routed R* numerator")?
        };
        if vakint::census::enabled() {
            vakint::census::record(
                "remainder_contraction",
                serde_json::json!({
                    "loops": self.region.loops,
                    "seconds": started.elapsed().as_secs_f64(),
                    "routed_bytes": routed.as_view().get_byte_size(),
                    "contracted_bytes": contracted.0.as_view().get_byte_size(),
                    "certified_closed_scalar": contracted.1,
                    "boundary": true,
                }),
            );
        }
        Ok(contracted)
    }

    /// Whether an ordinary remainder is contracted in the loop basis, without
    /// kept edge momenta (remainder-boundary admission A7).
    pub(super) fn plain_remainder_contraction() -> bool {
        !keeps_remainder() && !edge_contraction("remainder")
    }

    /// No kept edge-momentum combination was introduced.
    pub(super) fn assert_no_combinations(&self) -> Result<()> {
        ensure!(
            self.combinations.lock().unwrap().is_empty(),
            "a remainder boundary contraction kept edge-momentum combinations"
        );
        Ok(())
    }

    /// The loop edges `b_j` of the R* basis, in `j` order.
    pub(super) fn loop_basis_edges(&self) -> Vec<EdgeIndex> {
        self.basis.loop_edges.iter().copied().collect()
    }

    /// The converter graph (insertion supports identified to one vertex).
    pub(super) fn converter_graph(&self) -> &Graph {
        &self.graph
    }

    /// The routed expression contracted to dots, as the outer and physical
    /// paths consume it. It is contracted once, on first use.
    pub fn expression(&self) -> Result<&Atom> {
        Ok(&self.contracted_certified()?.0)
    }

    /// [`Self::expression`] and whether idenso certified its contraction a
    /// closed scalar: read from the facts cached on the contracted tensor,
    /// never by a scan (see [`Self::simplify_view_certified`]). Only the
    /// loop-basis contraction reports a certificate; the edge-keeping
    /// contractions rewrite their result afterwards and report `false`.
    pub(super) fn contracted_certified(&self) -> Result<&(Atom, bool)> {
        if let Some(contracted) = self.contracted.get() {
            return Ok(contracted);
        }
        let started = Instant::now();
        let contracted = {
            let _algebra = StageSpan::new("gamma.remainder.algebra");
            if keeps_remainder() {
                (
                    self.contract_keeping_edges(&self.routed)
                        .wrap_err("while preparing the routed R* tensor expression")?,
                    false,
                )
            } else if edge_contraction("remainder") || uv_child_edge_remainder() {
                (
                    Self::simplify_edge_basis(&self.routed)
                        .wrap_err("while preparing the routed R* tensor expression")?,
                    false,
                )
            } else {
                Self::simplify_view_certified(self.routed.as_view())
                    .wrap_err("while preparing the routed R* tensor expression")?
            }
        };
        if vakint::census::enabled() {
            vakint::census::record(
                "remainder_contraction",
                serde_json::json!({
                    "loops": self.region.loops,
                    "seconds": started.elapsed().as_secs_f64(),
                    "routed_bytes": self.routed.as_view().get_byte_size(),
                    "contracted_bytes": contracted.0.as_view().get_byte_size(),
                    "certified_closed_scalar": contracted.1,
                }),
            );
        }
        Ok(self.contracted.get_or_init(|| contracted))
    }

    /// Whether `expression` is this remainder's own contracted expression (the
    /// same object, as the ordinary outer integration passes it) and idenso
    /// certified that contraction a closed scalar. Any other expression, e.g.
    /// an IR correction's quotient product, has no certificate.
    fn certified_contraction(&self, expression: &Atom) -> bool {
        self.contracted
            .get()
            .is_some_and(|(contracted, certified)| {
                *certified && std::ptr::eq(contracted, expression)
            })
    }

    /// Mutable access to the contracted expression, for tests decorating it.
    /// A decorated expression keeps no certificate.
    #[cfg(test)]
    pub fn expression_mut(&mut self) -> Result<&mut Atom> {
        self.expression()?;
        let (expression, certified) = self.contracted.get_mut().expect("contracted above");
        *certified = false;
        Ok(expression)
    }

    /// Express a two-point amplitude in its first boundary momentum, keeping
    /// the edge orientation convention used by the graph's external routing.
    pub fn physical_external_routing(
        graph: &Graph,
        boundary: &SuBitGraph,
    ) -> Result<(EdgeIndex, Vec<Replacement>)> {
        let externals = graph
            .iter_edges_of(boundary)
            .map(|(pair, edge, _)| match pair {
                HedgePair::Unpaired { flow, .. } => Ok((edge, flow)),
                _ => Err(eyre!(
                    "physical pole reconstruction needs true external edges"
                )),
            })
            .collect::<Result<Vec<_>>>()?;
        ensure!(
            externals.len() == 2,
            "physical pole reconstruction currently requires a two-point graph"
        );
        let (external, flow) = externals[0];
        let routing = externals
            .iter()
            .map(|(edge, other_flow)| {
                let sign = if *edge == external || *other_flow != flow {
                    1
                } else {
                    -1
                };
                Replacement::new(
                    function!(GS.emr_mom, usize::from(*edge), W_.x___).to_pattern(),
                    Atom::num(sign) * function!(GS.emr_mom, usize::from(external), W_.x___),
                )
            })
            .collect::<Vec<_>>();
        Ok((external, routing))
    }

    /// Integrate a prepared propagator graph at its original external momentum.
    ///
    /// The vacuum converter intentionally forgets external flows and merges
    /// equal-mass two-edge bonds. Here the explicit quotient routing supplies
    /// the complete propagator momenta directly, preserving the physical scale.
    /// Tensor reduction averages transverse to that external momentum. All
    /// insertion coefficients enter Vakint before its Laurent expansion, so a
    /// child pole automatically promotes the required scalar-integral depth.
    pub fn integrate_physical(
        &self,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        vakint_settings: &VakintSettings,
        unique_id: usize,
    ) -> Result<Atom> {
        self.integrate_physical_numerator(
            self.physical_numerator()?,
            vakint,
            uv_settings,
            vakint_settings,
            unique_id,
        )
    }

    pub(super) fn physical_numerator(&self) -> Result<Atom> {
        // The physical path transports its numerator through its own momentum
        // maps; it keeps the loop-basis contraction.
        let contracted;
        let expression = if !self.combinations.lock().unwrap().is_empty() || keeps_remainder() {
            contracted = Self::simplify(&self.routed)
                .wrap_err("while preparing the routed R* tensor expression")?;
            &contracted
        } else {
            self.expression()?
        };
        Ok(expression
            .replace(function!(GS.den, W_.x___))
            .with(Atom::one()))
    }

    /// Certify a closed numerator in one canonical denominator frame.
    ///
    /// Only the compact particle graph enters graph isomorphism. Its chosen
    /// map transports the complete numerator, which is compared as an exact
    /// Atom after typed dummy normalization. Remaining automorphisms can miss
    /// matches conservatively. Open or unsupported tensors form singletons;
    /// they never fall back to canonizing a large expression-incidence graph.
    pub(super) fn physical_identity(
        &self,
        model: &Model,
        numerator: &Atom,
    ) -> Result<Option<(crate::uv::counterterm_identity::CountertermKey, Atom)>> {
        if UNSUPPORTED_PHYSICAL_MOMENTA
            .iter()
            .any(|&momentum| numerator.contains_symbol(momentum))
        {
            return Ok(None);
        }
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let tensor = SymbolicTensor::infer(cooking.try_cook(numerator.as_view())?)?;
        if !tensor.structure().logical_slots().is_empty() {
            return Ok(None);
        }
        let mut shell = self.graph.clone();
        let edges = shell
            .iter_edges()
            .map(|(pair, edge, _)| (pair.is_paired(), edge))
            .collect::<Vec<_>>();
        for (internal, edge) in edges {
            let retained = self.region.contains_edge(&shell, edge);
            let data = &mut shell.underlying[edge];
            data.num = Autogen::explicit(Atom::one());
            data.dod = Autogen::explicit(0);
            if internal && !retained {
                data.is_dummy = true;
            }
        }
        let vertices = shell
            .underlying
            .iter_nodes()
            .map(|(vertex, _, _)| vertex)
            .collect::<Vec<_>>();
        for vertex in vertices {
            let data = &mut shell.underlying[vertex];
            data.num = Autogen::explicit(Atom::one());
            data.dod = Autogen::explicit(0);
            data.vertex_rule = None;
        }
        shell.global_prefactor.num = Atom::one();
        shell.global_prefactor.projector = Atom::one();
        shell.loop_momentum_basis = self.basis.clone();
        let identity =
            crate::uv::counterterm_identity::CountertermIdentity::for_fixed_boundary_skeleton(
                &shell,
                &self.region.remaining,
                model,
            )?;
        let mut canonical = identity.canonicalize_integrand(numerator)?.expand();
        let mut matcher = SlotMatcher::default();
        let mut indexed = false;
        canonical.visitor(&mut |part| {
            if matches!(matcher.classify(part), SlotMatch::Explicit(_)) {
                indexed = true;
                false
            } else {
                !indexed
            }
        });
        if indexed {
            let tensor = SymbolicTensor::infer(cooking.try_cook(canonical.as_view())?)?;
            // The shared typed canonicalizer retains sum scopes. Symbolica
            // canonizes each top-level summand separately and treats subtrees
            // without dummy indices as scalar leaves. Invalid tensor/interface
            // inputs remain errors rather than silently becoming singletons.
            let normalized = tensor.canonize()?;
            canonical = cooking.uncook(normalized.expression().as_view());
        }
        Ok(Some((identity.key, canonical)))
    }

    pub(super) fn integrate_physical_numerator(
        &self,
        mut numerator: Atom,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        vakint_settings: &VakintSettings,
        unique_id: usize,
    ) -> Result<Atom> {
        let _stage = StageSpan::new("gamma.physical_integral");
        let routing_stage = StageSpan::new("gamma.physical_integral.routing");
        let mut timing = IntegrationTimings::new();
        if self.region.loops == 0 {
            drop(routing_stage);
            timing.routing = timing.lap();
            let result = Self::finish_tensor(numerator, unique_id)?;
            timing.finish = timing.lap();
            timing.log(&self.graph, unique_id, 0, "physical");
            return Ok(result);
        }
        let (external, external_routing) =
            Self::physical_external_routing(&self.graph, &self.region.boundary)?;
        numerator = numerator.replace_multiple(&external_routing);
        let minkowski = Minkowski {}.new_rep(GS.dim).to_symbolic([]);
        let mut reducer = TensorReducer::new(Atom::var(GS.dim)).with_external_vector(function!(
            GS.emr_mom,
            usize::from(external),
            &minkowski
        ))?;
        for edge in &self.basis.loop_edges {
            reducer = reducer.with_integrated_vector(function!(
                GS.emr_mom,
                usize::from(*edge),
                &minkowski
            ))?;
        }
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        timing.routing = timing.lap();
        drop(routing_stage);
        let materialization_stage =
            StageSpan::new("gamma.physical_integral.tensor_materialization");
        // Compact slash momenta inside gamma chains still carry tensor rank.
        // Expose those vectors before projection, using the same materializer
        // as the vacuum bridge. Scalar products and their inverses stay compact.
        let explicit = Integrated::materialize_tensor(&numerator)?;
        let source = reducer.admit(cooking.try_cook(explicit.as_view())?)?;
        timing.materialization = timing.lap();
        drop(materialization_stage);
        let projection_stage = StageSpan::new("gamma.physical_integral.tensor_reduction");
        let reduced = reducer
            .reduce(&source)
            .wrap_err("could not reduce a physical propagator tensor")?
            .into_expression()?;
        timing.tensor_projection = timing.lap();
        drop(projection_stage);
        let algebra_stage = StageSpan::new("gamma.physical_integral.algebra");
        numerator = Self::finish_tensor(cooking.uncook(reduced.expression().as_view()), unique_id)?;
        // Tensor restoration can spell scalar products as g(p,k). Normalize
        // them before mapping to Vakint's dot notation; represented metrics
        // and inverse external invariants retain their tensor/scalar roles.
        let scalar_products =
            SymbolicTensor::infer_indexed(cooking.try_cook(numerator.as_view())?)?.to_dots()?;
        numerator = cooking.uncook(scalar_products.expression().as_view());
        timing.algebra = timing.lap();
        drop(algebra_stage);
        let conversion_stage = StageSpan::new("gamma.physical_integral.vakint_conversion");

        let mut routing = self
            .basis
            .loop_edges
            .iter()
            .enumerate()
            .map(|(index, edge)| {
                Replacement::new(
                    function!(GS.emr_mom, usize::from(*edge), W_.x___).to_pattern(),
                    function!(S.k, index + 1, W_.x___),
                )
            })
            .collect::<Vec<_>>();
        routing.push(Replacement::new(
            function!(GS.emr_mom, usize::from(external), W_.x___).to_pattern(),
            function!(S.p, usize::from(external), W_.x___),
        ));
        numerator = numerator
            .replace_multiple(&routing)
            .replace(function!(
                SPENSO_TAG.dot,
                function!(W_.a_, W_.a___, &minkowski),
                function!(W_.b_, W_.b___, &minkowski)
            ))
            .with(function!(
                S.dot,
                function!(W_.a_, W_.a___),
                function!(W_.b_, W_.b___)
            ))
            .replace(function!(
                ETS.metric,
                Minkowski {}.to_symbolic([W_.a__]),
                Minkowski {}.to_symbolic([W_.b__])
            ))
            .with(function!(
                S.metric,
                Minkowski {}.to_symbolic([W_.a__]),
                Minkowski {}.to_symbolic([W_.b__])
            ));
        let mut propagators = Atom::one();
        for (index, (pair, edge, _)) in self.graph.iter_edges_of(&self.region.remaining).enumerate()
        {
            let HedgePair::Paired { source, sink } = pair else {
                continue;
            };
            let momentum = (self.basis.loop_atom::<Atom>(edge, GS.emr_mom, &[], true)
                + self.basis.ext_atom::<Atom>(edge, GS.emr_mom, &[], true))
            .replace_multiple(&external_routing)
            .replace_multiple(&routing);
            propagators *= function!(
                S.prop,
                index + 1,
                function!(
                    S.edge,
                    usize::from(self.graph.node_id(source)),
                    usize::from(self.graph.node_id(sink))
                ),
                momentum,
                0,
                self.graph[edge].extra_data.vakint_edge_power.unwrap_or(1)
            );
        }
        numerator *=
            symbolica::parse!(&uv_settings.vakint.additional_normalization).pow(self.region.loops);
        let input = numerator * function!(S.topo, propagators);
        debug_tags!(#uv, #integrated, #rstar, #trace;
            stage = "rstar_physical_remainder",
            log.input = input,
            "Integrating a UV preparation at its physical external momentum"
        );
        let mut settings = vakint_settings.clone();
        settings.use_dot_product_notation = true;
        timing.conversion = timing.lap();
        drop(conversion_stage);
        let backend_stage = StageSpan::new("gamma.physical_integral.backend");
        let integrated = vakint
            .evaluate_integral(&settings, input.as_view())
            .wrap_err("while integrating a physical UV spinney remainder")?;
        timing.backend = timing.lap();
        drop(backend_stage);
        let result = Self::finish_tensor(integrated, unique_id)?;
        timing.finish = timing.lap();
        timing.log(&self.graph, unique_id, self.region.loops, "physical");
        Ok(result)
    }

    /// `GAMMALOOP_IR_HARD=closed`: the `(coefficient, insertion)` moments of
    /// this hard component at Taylor order `degree` through polarized metric
    /// closures ([`super::closed`]), or a decline reason (today's path then).
    pub(super) fn closed_moments(
        &self,
        degree: usize,
        stats: &mut super::closed::ClosedStats,
    ) -> Result<std::result::Result<Vec<(Atom, Atom)>, String>> {
        if !self.combination_ids().is_empty() {
            return Ok(Err("kept_edges".into()));
        }
        // Fail-closed under a scoped budget (GAMMALOOP_IR_HARD_BUDGET_GIB over
        // the RSS now): an exhausted budget, an error or a panic declines the
        // term to today's path.
        let limit = super::closed::budget_extra_bytes()
            .and_then(|extra| process_rss_bytes().map(|rss| rss + extra));
        Ok(super::closed::guarded(
            limit,
            MEMORY_RELEASE.get().copied(),
            || {
                super::closed::closed_moments(
                    &self.routed,
                    &self.soft_edges,
                    degree,
                    // One root scope: idenso certifies a closed scalar result
                    // (frontier certificate), read before the uncook.
                    &|expression| Self::simplify_view_certified(expression.as_view()),
                    stats,
                )
            },
        ))
    }

    /// Put the bare hard term and its local counterterms in the same covariant
    /// basis before a quotient's pole projection. Hard scalar integrals remain
    /// unevaluated, including the dimension-dependent projector coefficients.
    pub fn project_hard(&self, expression: &Atom) -> Result<Atom> {
        let ids = self.combination_ids();
        let started = Instant::now();
        let input_bytes = expression.as_view().get_byte_size();
        if edge_project_log() {
            eprintln!(
                "EDGE_PROJECT start loops={} kept={} input_bytes={input_bytes}",
                self.region.loops,
                ids.len()
            );
        }
        if ids.is_empty() {
            let projected = project_hard_vacuum(expression, &self.basis.loop_edges.raw)?;
            if edge_project_log() {
                eprintln!(
                    "EDGE_PROJECT done loops={} kept=0 input_bytes={input_bytes} output_bytes={} seconds={:.3}",
                    self.region.loops,
                    projected.as_view().get_byte_size(),
                    started.elapsed().as_secs_f64()
                );
            }
            return Ok(projected);
        }
        // Kept edge-momentum combinations are integrated vectors too.
        let hard = self
            .basis
            .loop_edges
            .raw
            .iter()
            .copied()
            .chain(ids.iter().copied().map(EdgeIndex::from))
            .collect::<Vec<_>>();
        let projected = project_hard_vacuum(expression, &hard)
            .wrap_err("while projecting with kept edge momenta")?;
        if streamed::check("edge_contraction") && edge_check_admits(expression) {
            let rules = self.combination_rules();
            let loop_input = expression.replace_multiple(&rules);
            let reference = project_hard_vacuum(&loop_input, &self.basis.loop_edges.raw)?;
            let substituted = projected.replace_multiple(&rules).expand();
            let equal = (&substituted - &reference.expand()).expand().is_zero();
            eprintln!(
                "STREAM_CHECK stage=edge_projection bytes={} equal={equal}",
                expression.as_view().get_byte_size()
            );
        }
        let output_bytes = projected.as_view().get_byte_size();
        // The projected hard term is a scalar polynomial: when it is small,
        // rewrite it in the loop basis so that its integration sees today's
        // input (fewer distinct scalar products); a large one stays in the edge
        // basis, where its expansion stays bounded.
        let back = edge_contraction("back") && output_bytes <= edge_back_max_bytes();
        if edge_project_log() {
            eprintln!(
                "EDGE_PROJECT done loops={} kept={} input_bytes={input_bytes} output_bytes={output_bytes} seconds={:.3} back={back}",
                self.region.loops,
                ids.len(),
                started.elapsed().as_secs_f64()
            );
        }
        if back {
            return Ok(projected
                .replace_multiple(&self.combination_rules())
                .expand());
        }
        Ok(projected)
    }

    /// [`Self::simplify`] with every linear combination of this component's
    /// loop-edge momenta in an index slot kept as one independent vector
    /// `Q(OFFSET + n, slot)` (registered in `self.combinations`); momenta of
    /// other edges (root, boundary, soft) stay explicit, so Taylor expansions
    /// and projections in them act as before. The combinations are substituted
    /// at the routing step of `to_vakint_integrand`.
    fn contract_keeping_edges(&self, expression: &Atom) -> Result<Atom> {
        let loop_edges = self
            .basis
            .loop_edges
            .iter()
            .map(|edge| usize::from(*edge) as i64)
            .collect::<BTreeSet<_>>();
        let mut combinations = self.combinations.lock().unwrap();
        let renamed = rename_kept_edges(
            expression,
            &loop_edges,
            idenso::edge_output::enabled("flatten"),
            &mut combinations,
        );
        drop(combinations);
        Self::simplify_view(renamed.as_view()).map_err(|error| {
            if let Ok(dir) = std::env::var("GAMMALOOP_EDGE_DEBUG_DUMP") {
                let _ = std::fs::write(
                    format!("{dir}/keep-failure-{}.txt", std::process::id()),
                    format!("input {expression}\nrenamed {renamed}\n"),
                );
            }
            error.wrap_err("while contracting with kept edge momenta")
        })
    }

    /// `Q(OFFSET + n, x___) -> Σ c_e Q(e, x___)` for every kept combination.
    pub(super) fn combination_rules(&self) -> Vec<Replacement> {
        self.combinations
            .lock()
            .unwrap()
            .iter()
            .map(|(terms, id)| {
                let pattern = function!(GS.emr_mom, Atom::num(*id), W_.x___);
                let image = Atom::add_many(
                    terms
                        .iter()
                        .map(|(edge, coefficient)| {
                            coefficient * function!(GS.emr_mom, Atom::num(*edge), W_.x___)
                        })
                        .collect::<Vec<_>>()
                        .iter(),
                );
                Replacement::new(pattern.to_pattern(), image.to_pattern())
            })
            .collect()
    }

    fn combination_ids(&self) -> Vec<usize> {
        self.combinations
            .lock()
            .unwrap()
            .values()
            .map(|id| *id as usize)
            .collect()
    }

    /// Hard Taylor expansion leaves the reference mass fixed. These momenta
    /// are external to this component, even when they are loops of its parent.
    ///
    /// The scaling and the truncation act on the routed product before any
    /// index is contracted, so only the kept Taylor terms are contracted
    /// (as `approx::local_4d` does for local counterterms). Contraction is
    /// linear and the parameter is a scalar, so this is the Taylor
    /// polynomial of the contracted expression, [`Self::taylor_contracted`].
    pub fn taylor(&self, degree: usize) -> Result<Atom> {
        let _stage = StageSpan::new("gamma.remainder.hard_taylor");
        let parameter = symbol!("rstar_hard_taylor"; Scalar);
        streamed::dump("hard_routed", self.routed.as_view());
        let scaled = Self::scale_edges(&self.routed, &self.soft_edges, parameter);
        let expanded = TaylorExpansion::new(parameter, degree).evaluate(&scaled)?;
        drop(scaled);
        let started = Instant::now();
        let truncated = {
            let _algebra = StageSpan::new("gamma.remainder.hard_taylor_algebra");
            streamed::report_size("hard_taylor", expanded.as_view());
            streamed::dump("hard_taylor", expanded.as_view());
            if keeps_hard() {
                let kept = self
                    .contract_keeping_edges(&expanded)
                    .wrap_err("while simplifying the R* hard Taylor polynomial")?;
                if streamed::check("edge_contraction") && edge_check_admits(&expanded) {
                    let reference = Self::simplify(&expanded)
                        .wrap_err("while simplifying the R* hard Taylor polynomial")?;
                    let substituted = kept.replace_multiple(&self.combination_rules()).expand();
                    let equal = (&substituted - &reference.expand()).expand().is_zero();
                    eprintln!(
                        "STREAM_CHECK stage=edge_keep bytes={} equal={equal}",
                        expanded.as_view().get_byte_size()
                    );
                }
                let kept_bytes = kept.as_view().get_byte_size();
                if edge_contraction("pre") && kept_bytes <= edge_pre_max_bytes() {
                    // Small contracted hard term: rewrite it in the loop basis
                    // and project it there (the loop-basis projection is the
                    // cheaper one for few distinct scalar products).
                    let rules = self.combination_rules();
                    self.combinations.lock().unwrap().clear();
                    if edge_project_log() {
                        eprintln!(
                            "EDGE_PRE loops={} kept_bytes={kept_bytes} back=true",
                            self.region.loops
                        );
                    }
                    kept.replace_multiple(&rules).expand()
                } else {
                    if edge_project_log() {
                        eprintln!(
                            "EDGE_PRE loops={} kept_bytes={kept_bytes} back=false",
                            self.region.loops
                        );
                    }
                    kept
                }
            } else if edge_contraction("hard") {
                let edge = Self::simplify_edge_basis(&expanded)
                    .wrap_err("while simplifying the R* hard Taylor polynomial")?;
                if streamed::check("edge_contraction") {
                    let reference = Self::simplify(&expanded)
                        .wrap_err("while simplifying the R* hard Taylor polynomial")?;
                    let identical = edge == reference;
                    let equal = identical || (&edge - &reference).expand().is_zero();
                    eprintln!(
                        "STREAM_CHECK stage=edge_contraction bytes={} identical={identical} equal={equal}",
                        expanded.as_view().get_byte_size()
                    );
                }
                edge
            } else if edge_contraction("fallback") {
                self.simplify_hard_taylor_or_keep(&expanded)
                    .wrap_err("while simplifying the R* hard Taylor polynomial")?
            } else {
                Self::simplify_hard_taylor(&expanded)
                    .wrap_err("while simplifying the R* hard Taylor polynomial")?
            }
        };
        if vakint::census::enabled() {
            vakint::census::record(
                "hard_taylor_contraction",
                serde_json::json!({
                    "loops": self.region.loops,
                    "seconds": started.elapsed().as_secs_f64(),
                    "taylor_bytes": expanded.as_view().get_byte_size(),
                    "contracted_bytes": truncated.as_view().get_byte_size(),
                }),
            );
        }
        if cfg!(test) || hard_taylor_check() {
            let reference = self.taylor_contracted(degree)?;
            let identical = truncated == reference;
            let equal = identical || (&truncated - &reference).expand().is_zero();
            if cfg!(test) {
                assert!(
                    equal,
                    "hard Taylor before contraction differs from the contracted path\nnew: {truncated}\nold: {reference}"
                );
            } else {
                eprintln!(
                    "HARD_TAYLOR_CHECK loops={} degree={degree} identical={identical} equal={equal}",
                    self.region.loops
                );
            }
        }
        Ok(truncated)
    }

    /// `GAMMALOOP_EDGE_CONTRACTION=fallback`: today's loop-basis contraction
    /// under a scoped memory budget (process RSS may grow by
    /// `GAMMALOOP_EDGE_FALLBACK_GIB`, default 8, over its value at the start).
    /// When that budget or the process budget is exhausted (the 27105
    /// signature: a frontier that grows without bound), this component is
    /// contracted with its edge-momentum combinations kept as independent
    /// vectors instead, and [`Self::project_hard`] projects it in that basis.
    /// Components that finish within the budget are unchanged.
    fn simplify_hard_taylor_or_keep(&self, expanded: &Atom) -> Result<Atom> {
        use symbolica::utils::resources::{ResourceBudget, ResourceLimits};
        let limit = process_rss_bytes().map(|rss| rss + edge_fallback_bytes());
        let budget = ResourceBudget::new(ResourceLimits {
            max_process_rss_bytes: limit,
            ..ResourceLimits::default()
        });
        let started = Instant::now();
        let attempt = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            budget.run(|| Self::simplify(expanded))
        }));
        let payload = match attempt {
            Ok(result) => return result,
            Err(payload) => payload,
        };
        let message = payload
            .downcast_ref::<String>()
            .cloned()
            .or_else(|| {
                payload
                    .downcast_ref::<&str>()
                    .map(|text| (*text).to_owned())
            })
            .unwrap_or_default();
        if !message.starts_with("resource budget exhausted") {
            std::panic::resume_unwind(payload);
        }
        drop(payload);
        eprintln!(
            "EDGE_FALLBACK loops={} taylor_bytes={} after={:.1}s reason={message}",
            self.region.loops,
            expanded.as_view().get_byte_size(),
            started.elapsed().as_secs_f64()
        );
        if let Some(release) = MEMORY_RELEASE.get() {
            release();
            std::thread::sleep(std::time::Duration::from_millis(1500));
        }
        let started = Instant::now();
        let kept = self.contract_keeping_edges(expanded)?;
        eprintln!(
            "EDGE_FALLBACK done kept={} contracted_bytes={} seconds={:.1}",
            self.combination_ids().len(),
            kept.as_view().get_byte_size(),
            started.elapsed().as_secs_f64()
        );
        Ok(kept)
    }

    /// [`Self::simplify`] of a hard Taylor polynomial. The factorized
    /// contraction of a deeply nested product of many index-sharing sums
    /// (27105: a 41 KB polynomial) can grow without bound; under an active
    /// resource budget (`BETA_PROCESS_RSS_LIMIT_GIB`) it aborts, and the
    /// polynomial is then expanded and simplified in chunks of summands:
    /// the same value, not the same factored form, so this path runs only
    /// where the whole simplification cannot finish
    /// (`GAMMALOOP_ALGEBRA_LEGACY=hard_taylor_fallback` lets the abort fail
    /// the diagram instead).
    fn simplify_hard_taylor(expanded: &Atom) -> Result<Atom> {
        let attempt =
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| Self::simplify(expanded)));
        match attempt {
            Ok(result) => result,
            Err(payload) => {
                let message = payload
                    .downcast_ref::<String>()
                    .cloned()
                    .or_else(|| {
                        payload
                            .downcast_ref::<&str>()
                            .map(|text| (*text).to_owned())
                    })
                    .unwrap_or_default();
                if !message.starts_with("resource budget exhausted")
                    || streamed::legacy("hard_taylor_fallback")
                {
                    std::panic::resume_unwind(payload);
                }
                eprintln!(
                    "HARD_TAYLOR_FALLBACK bytes={} reason={message}",
                    expanded.as_view().get_byte_size()
                );
                // Return the aborted frontier's pages to the kernel before the
                // budget is checked again (the allocator keeps freed memory).
                if let Some(release) = MEMORY_RELEASE.get() {
                    release();
                    std::thread::sleep(std::time::Duration::from_millis(1500));
                }
                let started = std::time::Instant::now();
                let expanded = expanded.expand();
                let simplified = Self::simplify_summands(expanded.as_view())?;
                eprintln!(
                    "HARD_TAYLOR_FALLBACK done terms_bytes={} result_bytes={} seconds={:.1}",
                    expanded.as_view().get_byte_size(),
                    simplified.as_view().get_byte_size(),
                    started.elapsed().as_secs_f64()
                );
                Ok(simplified)
            }
        }
    }

    /// The previous hard Taylor path: contract the whole routed expression,
    /// then scale and truncate it. Kept as the reference for [`Self::taylor`].
    fn taylor_contracted(&self, degree: usize) -> Result<Atom> {
        let parameter = symbol!("rstar_hard_taylor"; Scalar);
        let scaled = Self::scale_edges(self.expression()?, &self.soft_edges, parameter);
        let expanded = TaylorExpansion::new(parameter, degree).evaluate(&scaled)?;
        Self::simplify(&expanded).wrap_err("while simplifying the R* hard Taylor polynomial")
    }

    #[allow(clippy::too_many_arguments)]
    pub fn integrate(
        &self,
        expression: &Atom,
        root_graph: &Graph,
        root_edges: &[EdgeIndex],
        degree: usize,
        reference_mass_squared: Symbol,
        auxiliary_mass_squared: &Atom,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        vakint_settings: &VakintSettings,
        unique_id: usize,
        root_projector: Option<&RootProjector<'_>>,
    ) -> Result<Atom> {
        let mut timing = IntegrationTimings::new();
        let polynomial = self.outer_taylor(
            expression,
            root_edges,
            degree,
            reference_mass_squared,
            auxiliary_mass_squared,
            &mut timing,
        )?;
        // GAMMALOOP_CLOSED_SKIP_SUMMAND_SIMPLIFY (default off): on the closed
        // root path idenso contracted the remainder (with the root closure
        // metric) to a closed scalar before the Taylor step and certified it
        // so. The Taylor step only routes, scales and truncates it (dot is
        // linear, so routing keeps it in dot form; den is differentiated in
        // its scalar fourth argument): it introduces no index, and this pass
        // has nothing to contract. The skip is taken only for a root-projected
        // integral of the remainder's own certified contraction; an IR
        // correction (a quotient product, closed after its K) has no
        // certificate and is simplified as before. The certificate is a cached
        // fact, never a scan; `check` adds the syntactic closed_dot_scalar test
        // as a cross-check.
        let skip = closed_summand_skip();
        let certified = skip != ClosedSkipMode::Off
            && root_projector.is_some()
            && self.certified_contraction(expression);
        let closed = skip == ClosedSkipMode::Check
            && root_projector.is_some()
            && closed_dot_scalar(polynomial.as_view());
        if skip != ClosedSkipMode::Off {
            eprintln!(
                "CLOSED_SKIP site=outer_taylor id={unique_id} loops={} bytes={} root_projected={} certified={certified} closed={} action={}",
                self.region.loops,
                polynomial.as_view().get_byte_size(),
                root_projector.is_some(),
                if skip == ClosedSkipMode::Check {
                    if closed { "true" } else { "false" }
                } else {
                    "unchecked"
                },
                if certified && skip == ClosedSkipMode::Skip {
                    "skipped"
                } else {
                    "simplified"
                },
            );
        }
        let expanded = if certified && skip == ClosedSkipMode::Skip {
            polynomial
        } else {
            let expanded = Self::simplify_summands(polynomial.as_view())
                .inspect_err(|_| edge_debug_dump("outer_simplify", &polynomial))
                .wrap_err("while simplifying the compensated R* outer Taylor polynomial")?;
            if certified || closed {
                let identical = expanded == polynomial;
                let equal = identical || (&expanded - &polynomial).expand().is_zero();
                eprintln!(
                    "CLOSED_SKIP site=outer_taylor id={unique_id} check certified={certified} closed={closed} identical={identical} equal={equal}"
                );
            } else if skip == ClosedSkipMode::Check && root_projector.is_some() {
                // Diagnostics (with GAMMALOOP_ALGEBRA_DUMP): a root-projected
                // polynomial that failed the closed-scalar test, before and
                // after this pass, and whether the pass changed it.
                streamed::dump("outer_taylor_open", polynomial.as_view());
                streamed::dump("outer_taylor_open_simplified", expanded.as_view());
                eprintln!(
                    "CLOSED_SKIP site=outer_taylor id={unique_id} open unchanged={} closed_after={}",
                    expanded == polynomial,
                    closed_dot_scalar(expanded.as_view())
                );
            }
            drop(polynomial);
            expanded
        };
        timing.algebra = timing.lap();
        self.integrate_polynomial(
            expanded,
            root_graph,
            vakint,
            uv_settings,
            vakint_settings,
            unique_id,
            root_projector,
            timing,
        )
    }

    /// The outer Taylor polynomial of `expression` as [`Self::integrate`]
    /// takes it: routed in this component's basis, the root edges scaled by
    /// the Taylor parameter, the reference mass compensated, truncated at total
    /// degree `degree` and evaluated at one. Laps the routing and Taylor
    /// timings.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn outer_taylor(
        &self,
        expression: &Atom,
        root_edges: &[EdgeIndex],
        degree: usize,
        reference_mass_squared: Symbol,
        auxiliary_mass_squared: &Atom,
        timing: &mut IntegrationTimings,
    ) -> Result<Atom> {
        // The root's boundary carriers: the graph's external edges, or the full
        // crown of a proper support evaluated as its own R* root.
        // Quotient UV maps return original Q labels. Route newly introduced
        // component edges, preserving this component's independent hard loops.
        let routed = if self.region.loops == 0 {
            Cow::Borrowed(expression)
        } else {
            let routing = self.graph.normal_emr_replacement(
                &self.region.remaining,
                &self.basis,
                &[W_.x___],
                HedgePair::is_paired,
            );
            Cow::Owned(streamed::map_term_local(
                "the outer-Taylor routing",
                expression.as_view(),
                &|piece| piece.replace_multiple(&routing),
            )?)
        };
        let allowed_edges = self
            .basis
            .loop_edges
            .iter()
            .chain(root_edges)
            .map(|edge| usize::from(*edge))
            .chain(self.combination_ids())
            .collect::<BTreeSet<_>>();
        let mut unresolved = BTreeSet::new();
        routed.visitor(&mut |part| {
            if let AtomView::Fun(momentum) = part
                && momentum.get_symbol() == GS.emr_mom
                && let Ok(edge) = usize::try_from(momentum.get(0))
                && !allowed_edges.contains(&edge)
            {
                unresolved.insert(edge);
            }
            true
        });
        ensure!(
            unresolved.is_empty(),
            "an R* quotient left unintegrated boundary loop carriers: {unresolved:?}"
        );
        timing.routing = timing.lap();
        let parameter = symbol!("rstar_outer_taylor"; Scalar);
        let taylor = TaylorExpansion::new(parameter, degree);
        // Each intermediate is dropped as soon as its successor exists: on the
        // largest diagrams each of them is hundreds of megabytes.
        let scaled = streamed::map_term_local(
            "the outer-Taylor edge scaling",
            routed.as_view(),
            &|piece| Self::scale_edges_view(piece, root_edges, parameter),
        )?;
        drop(routed);
        streamed::guard_checkpoint("the outer-Taylor mass compensation")?;
        let compensated =
            taylor.compensate_mass(&scaled, reference_mass_squared, auxiliary_mass_squared);
        drop(scaled);
        streamed::guard_checkpoint("the outer-Taylor expansion")?;
        // Instrumentation (opt-in, no value depends on it): the P_i x D_j
        // shape of this product, reported by the next to_vakint call.
        if self.region.loops != 0 && counts::enabled() {
            counts::set_taylor_shape(counts::taylor_shape(
                compensated.as_view(),
                parameter,
                degree,
            ));
        }
        let polynomial = taylor.evaluate(&compensated)?;
        drop(compensated);
        timing.taylor = timing.lap();
        Ok(polynomial)
    }

    /// The rest of [`Self::integrate`] after the outer Taylor step and its
    /// summand simplification: root projection, conversion and integration.
    #[allow(clippy::too_many_arguments)]
    fn integrate_polynomial(
        &self,
        expanded: Atom,
        root_graph: &Graph,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        vakint_settings: &VakintSettings,
        unique_id: usize,
        root_projector: Option<&RootProjector<'_>>,
        mut timing: IntegrationTimings,
    ) -> Result<Atom> {
        // The root is the only boundary which may be projected here. Child
        // and quotient pole maps are already complete, and no later pole map
        // can multiply finite contraction terms by another insertion pole.
        // Project after Taylor expansion so physical external momenta in a
        // projector are spectators, rather than differentiated source labels.
        let expanded = if let Some(projector) = root_projector {
            let projected = projector(&expanded)
                .inspect_err(|_| edge_debug_dump("root_projector", &expanded))
                .wrap_err("while projecting the complete R* outer Taylor polynomial")?;
            if super::boundary::callback_check_active() {
                // Remainder-boundary check mode: the callback against c·σ.
                let regulated = Integrated::dimensionally_regularized(&projected);
                super::boundary::check_callback(&expanded, &regulated);
                drop(expanded);
                regulated
            } else {
                drop(expanded);
                Integrated::dimensionally_regularized(&projected)
            }
        } else {
            expanded
        };
        timing.tensor_projection = timing.lap();
        let integrated = if self.region.loops == 0 {
            expanded
        } else {
            // By value: to_vakint drops it once it has been materialized.
            let rules = self.combination_rules();
            if !rules.is_empty() {
                edge_debug_dump("to_vakint_input", &expanded);
            }
            let input = crate::uv::approx::integrated::with_momentum_combinations(
                &rules,
                || -> Result<_> {
                    to_vakint_integrand(
                        expanded,
                        &self.graph.underlying,
                        &self.region.remaining,
                        &self.graph.empty_subgraph::<SuBitGraph>(),
                        &uv_settings.vakint,
                        false,
                    )
                    .wrap_err("while converting with kept edge momenta")
                },
            )?;
            timing.conversion = timing.lap();
            let integrated = integrate_single_mass(vakint, vakint_settings, &input)
                .inspect_err(|_| edge_debug_dump("single_mass_error", &Atom::num(1)))?;
            timing.backend = timing.lap();
            integrated
        };
        let result = Self::finish_tensor(integrated, unique_id)?;
        timing.finish = timing.lap();
        timing.log(root_graph, unique_id, self.region.loops, "single_mass");
        Ok(result)
    }

    /// Use exactly the ordinary R operation's tensor and dimension restoration.
    pub fn finish_tensor(expression: Atom, unique_id: usize) -> Result<Atom> {
        let _stage = StageSpan::new("gamma.tensor_restoration");
        // Restore compact external scalar products before opening Vakint's
        // remaining dot shorthand. Both physical quotients and normalized
        // root projectors can contain reciprocals of these invariants;
        // expanding them would create inverse indexed tensors. Mixed loop
        // products and already indexed arguments retain the strict converter.
        let minkowski = Minkowski {}.new_rep(GS.dim).to_symbolic([]);
        let expression = expression
            .replace(function!(
                S.dot,
                function!(S.p, W_.i_),
                function!(S.p, W_.j_)
            ))
            .with(function!(
                SPENSO_TAG.dot,
                function!(S.p, W_.i_, &minkowski),
                function!(S.p, W_.j_, &minkowski)
            ));
        let indexed = Vakint::convert_from_dot_notation(expression.as_view());
        let restored = Integrated::restore_numerator(indexed, unique_id);
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let reduced = SymbolicTensor::infer_indexed(cooking.try_cook(restored.as_view())?)?
            .contract(ContractSettings {
                rank_one: false,
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            })?
            .to_dots()?;
        let regulated =
            Integrated::dimensionally_regularized(&cooking.uncook(reduced.expression().as_view()));
        // Keep collected trace/chain powers for their domain algebra. Opening
        // a powered trace first would give every copy the same spinor dummies,
        // turning two independent closed traces into repeated gamma factors.
        let collected = SymbolicTensor::infer_indexed(cooking.try_cook(regulated.as_view())?)?
            .contract(ContractSettings::default())?
            .simplify_algebra(&AlgebraSettings {
                gamma: Some(GammaSimplifySettings::canonical()),
                epsilon: true,
                color: Some(ColorSimplifySettings::default()),
                ..Default::default()
            })?;
        // Clifford contraction can create new scalar dimensions. Expand
        // those before the caller takes poles or inserts this tensor in an
        // enclosing graph; open representation dimensions remain metadata.
        let reduced = cooking.uncook(collected.expression().as_view());
        Integrated::warn_repeated_lorentz_indices(
            &reduced,
            "R*",
            Some(collected.reduction_status()),
        );
        Ok(Integrated::dimensionally_regularized(&reduced))
    }

    pub(super) fn scale_edges(expression: &Atom, edges: &[EdgeIndex], parameter: Symbol) -> Atom {
        Self::scale_edges_view(expression.as_view(), edges, parameter)
    }

    /// [`Self::scale_edges`] of a view (a term-local map, see
    /// [`streamed::map_term_local`]).
    pub(super) fn scale_edges_view(
        expression: AtomView<'_>,
        edges: &[EdgeIndex],
        parameter: Symbol,
    ) -> Atom {
        let replacements = edges
            .iter()
            .map(|edge| {
                let momentum = function!(GS.emr_mom, usize::from(*edge), W_.x___);
                Replacement::new(momentum.to_pattern(), Atom::var(parameter) * momentum)
            })
            .collect::<Vec<_>>();
        expression.replace_multiple(&replacements)
    }

    /// The edge-momentum combinations of `expression` (sums `Σ c_e Q(e, slot)`
    /// in an index slot) renamed to independent vectors, and the rules that
    /// substitute them back; `None` when there is no combination.
    fn edge_rename(expression: &Atom) -> Option<(Atom, Vec<Replacement>)> {
        let head = GS.emr_mom;
        let mut combinations = std::collections::BTreeMap::new();
        let renamed = rename_edge_basis(
            expression,
            idenso::edge_output::enabled("flatten"),
            &mut combinations,
        );
        if combinations.is_empty() {
            return None;
        }
        let rules = combinations
            .iter()
            .map(|(terms, id)| {
                let pattern = function!(head, Atom::num(*id), W_.x___);
                let image = Atom::add_many(
                    terms
                        .iter()
                        .map(|(edge, coefficient)| {
                            coefficient * function!(head, Atom::num(*edge), W_.x___)
                        })
                        .collect::<Vec<_>>()
                        .iter(),
                );
                Replacement::new(pattern.to_pattern(), image.to_pattern())
            })
            .collect::<Vec<_>>();
        Some((renamed, rules))
    }

    /// Substitute the edge combinations back ([`Self::edge_rename`]): a
    /// term-local replacement, applied in checkpointed pieces.
    fn edge_substitute(expression: AtomView<'_>, rules: &[Replacement]) -> Result<Atom> {
        streamed::map_term_local("the edge-combination substitution", expression, &|piece| {
            piece.replace_multiple(rules)
        })
    }

    /// [`Self::simplify`] with every momentum linear combination in an index
    /// slot (an edge momentum after routing, `Σ c_e Q(e, slot)`) replaced by
    /// one independent vector before the contraction, and substituted back
    /// afterwards, factor by factor and without expanding. Contraction is
    /// linear in each vector, so this is the same tensor expression; the
    /// contraction no longer distributes the sums of the loop-momentum basis
    /// (on 27105's hard Taylor input: 47 s and 2.7 GB instead of > 90 GiB).
    pub(super) fn simplify_edge_basis(expression: &Atom) -> Result<Atom> {
        let Some((renamed, rules)) = Self::edge_rename(expression) else {
            return Self::simplify_view(expression.as_view());
        };
        let simplified = Self::simplify_view(renamed.as_view())?;
        drop(renamed);
        Self::edge_substitute(simplified.as_view(), &rules)
    }

    /// `GAMMALOOP_UV_STREAM=on|auto` ([`uv_stream`]): this remainder's
    /// edge-basis contraction as a stream of batches of pieces, each piece
    /// contracted and substituted, each batch handed to `sink`.
    ///
    /// The routed remainder of a degree-0 UV child is a small product (a few
    /// KB) of Feynman-rule factors with several momentum sums, whose
    /// contraction distributes into GBs (nf0 32: 10 KB → 4 GB, 40 GB peak,
    /// then > 64 GB in the substitution). The sums with the most terms are split
    /// into their terms until there are at least `GAMMALOOP_UV_STREAM_PIECES`
    /// pieces (the cartesian product of the chosen terms times the other
    /// factors). Contraction, substitution and every later step of the
    /// integration are linear, and a normalized sum is canonical, so the sum of
    /// the integrated pieces is the integrated whole.
    ///
    /// Returns false, having emitted nothing, when the stream does not apply:
    /// switch off, not a degree-0 UV child in the edge basis, nothing to split,
    /// the first piece introduces a dummy index, or (`auto`) fewer than
    /// `GAMMALOOP_UV_STREAM_MIN_COMBOS` term combinations. A later piece that
    /// introduces a dummy index is an error (never observed).
    pub(super) fn stream_contracted(
        &self,
        sink: &mut dyn FnMut(Atom) -> Result<()>,
    ) -> Result<bool> {
        let mode = uv_stream();
        if mode == UvStream::Off
            || !uv_child_edge_remainder()
            || keeps_remainder()
            || self.contracted.get().is_some()
        {
            return Ok(false);
        }
        let _algebra = StageSpan::new("gamma.remainder.algebra");
        let renamed = Self::edge_rename(&self.routed);
        let input = match &renamed {
            Some((renamed, _)) => renamed.as_view(),
            None => self.routed.as_view(),
        };
        let rules = renamed
            .as_ref()
            .map(|(_, rules)| rules.as_slice())
            .unwrap_or(&[]);
        let AtomView::Mul(product) = input else {
            return Ok(false);
        };
        let factors = product.iter().collect::<Vec<_>>();
        let mut sums = factors
            .iter()
            .enumerate()
            .filter_map(|(i, f)| match f {
                AtomView::Add(sum) if sum.get_nargs() > 1 => Some((i, sum.get_nargs())),
                _ => None,
            })
            .collect::<Vec<_>>();
        let combos = sums
            .iter()
            .fold(1usize, |product, (_, terms)| product.saturating_mul(*terms));
        if mode == UvStream::Auto && combos < uv_stream_min_combos() {
            return Ok(false);
        }
        sums.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        let target = uv_stream_pieces();
        let mut chosen = Vec::new();
        let mut count = 1usize;
        for (index, terms) in sums {
            if count >= target {
                break;
            }
            chosen.push(index);
            count *= terms;
        }
        if chosen.is_empty() {
            return Ok(false);
        }
        let choices = chosen
            .iter()
            .map(|&i| match factors[i] {
                AtomView::Add(sum) => sum.iter().collect::<Vec<_>>(),
                _ => unreachable!("chosen factors are sums"),
            })
            .collect::<Vec<_>>();
        let rest = Atom::mul_many(
            factors
                .iter()
                .enumerate()
                .filter(|(i, _)| !chosen.contains(i))
                .map(|(_, f)| *f),
        );
        let started = Instant::now();
        let contract_piece = |digits: &[usize]| -> Result<Option<Atom>> {
            let mut piece = rest.clone();
            for (terms, &digit) in choices.iter().zip(digits) {
                piece = piece * terms[digit].to_owned();
            }
            let contracted = Self::simplify_view(piece.as_view())?;
            if streamed::introduces_dummies(piece.as_view(), contracted.as_view()) {
                return Ok(None);
            }
            drop(piece);
            Ok(Some(if rules.is_empty() {
                contracted
            } else {
                Self::edge_substitute(contracted.as_view(), rules)?
            }))
        };
        let batch_limit = uv_stream_batch_bytes();
        let mut batch = streamed::SumAccumulator::default();
        let (mut batch_bytes, mut emitted, mut pieces, mut batches) =
            (0usize, 0usize, 0usize, 0usize);
        let mut digits = vec![0usize; choices.len()];
        'pieces: loop {
            let Some(mapped) = contract_piece(&digits)? else {
                if pieces == 0 {
                    eprintln!("UV_STREAM declined=dummy loops={}", self.region.loops);
                    return Ok(false);
                }
                return Err(eyre!(
                    "a UV stream piece introduced a dummy index after {pieces} pieces"
                ));
            };
            let bytes = mapped.as_view().get_byte_size();
            emitted += bytes;
            batch_bytes += bytes;
            pieces += 1;
            batch.push(mapped);
            if batch_bytes >= batch_limit {
                sink(std::mem::take(&mut batch).finish())?;
                batch_bytes = 0;
                batches += 1;
            }
            let mut k = 0;
            loop {
                if k == digits.len() {
                    break 'pieces;
                }
                digits[k] += 1;
                if digits[k] < choices[k].len() {
                    break;
                }
                digits[k] = 0;
                k += 1;
            }
        }
        if batch_bytes > 0 || batches == 0 {
            sink(batch.finish())?;
            batches += 1;
        }
        eprintln!(
            "UV_STREAM pieces={pieces} batches={batches} combos={combos} in_bytes={} emitted_bytes={emitted} seconds={:.1} loops={}",
            input.get_byte_size(),
            started.elapsed().as_secs_f64(),
            self.region.loops
        );
        Ok(true)
    }

    pub(super) fn simplify(expression: &Atom) -> Result<Atom> {
        debug_tags!(#uv, #integrated, #rstar, #trace;
            stage = "rstar_tensor_input",
            log.expression = expression,
            "Preparing an explicit R* expression for native tensor algebra"
        );
        Self::simplify_view(expression.as_view())
    }

    /// [`Self::simplify`] summand by summand: the simplification acts on
    /// every summand alone, and a normalized sum is canonical. A large sum
    /// is simplified in chunks of summands and added up once, unless a chunk
    /// introduces a dummy index, which makes the whole sum be simplified.
    pub(super) fn simplify_summands(expression: AtomView<'_>) -> Result<Atom> {
        let limit = streamed::piece_limit();
        if expression.get_byte_size() < limit || !matches!(expression, AtomView::Add(_)) {
            return Self::simplify_view(expression);
        }
        let mut accumulated = streamed::SumAccumulator::default();
        let mut introduced = false;
        let outcome = streamed::for_each_summand_chunk(expression, limit, &mut |chunk| {
            let simplified = Self::simplify_view(chunk)?;
            if streamed::introduces_dummies(chunk, simplified.as_view()) {
                introduced = true;
                return Err(eyre!("a simplification introduced a dummy index"));
            }
            accumulated.push(simplified);
            Ok(())
        });
        let simplified = match outcome {
            Ok(()) => accumulated.finish(),
            Err(_) if introduced => return Self::simplify_view(expression),
            Err(error) => return Err(error),
        };
        if streamed::check_streamed(expression) {
            let whole = Self::simplify_view(expression)?;
            streamed::report_check("rstar_simplify", expression, &simplified, &whole);
        }
        Ok(simplified)
    }

    fn simplify_view(expression: AtomView<'_>) -> Result<Atom> {
        Ok(Self::simplify_view_certified(expression)?.0)
    }

    /// [`Self::simplify_view`], also returning whether idenso certified its
    /// result a closed scalar: no open slot, no contraction or identity
    /// candidate left. The flag is the producer fact cached on the tensor
    /// that `simplify_algebra(..).to_dots()` returns, read before the uncook
    /// (which leaves a plain Atom without facts); it is never computed by a
    /// scan, so an uncertified result reports `false`.
    pub(super) fn simplify_view_certified(expression: AtomView<'_>) -> Result<(Atom, bool)> {
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        // As in the ordinary R simplifier, retain open color tensors here.
        // Reduce color identities after the quotient has closed their indices.
        // Leave trace/chain powers intact for the same reason as finish_tensor:
        // their internal spinor sums must remain independent between copies.
        // The certificate accessors are on the indexed type (DESIGN Addendum
        // A(3)): no conversion is needed to read them.
        let simplified = SymbolicTensor::infer_indexed(cooking.try_cook(expression)?)?
            .simplify_algebra(&AlgebraSettings {
                gamma: Some(GammaSimplifySettings::canonical()),
                epsilon: true,
                // `IDENSO_EDGE_OUTPUT=dots`: the frontier emits the dot notation itself.
                contract: if idenso::edge_output::enabled("dots") {
                    idenso::tensor::AlgebraContraction::Dots
                } else {
                    Default::default()
                },
                ..Default::default()
            })?
            .to_dots()?;
        let certified = simplified.is_certified_closed_scalar();
        Ok((cooking.uncook(simplified.expression().as_view()), certified))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        finalized_runtime_dot, graph::parse::IntoFinalizedRuntimeGraph,
        initialisation::test_initialise, utils::load_generic_model, uv::rstar::RStarPlan,
    };
    use vakint::symbols::S;

    #[test]
    fn simplify_powered_dirac_trace_preserves_open_color_tensor() -> Result<()> {
        use crate::numerator::aind::Aind;
        use idenso::{
            IndexTooling,
            representations::{Bispinor, ColorAdjoint, ColorAntiFundamental, ColorFundamental},
        };
        use spenso::structure::slot::{DummyAind, ParseableAind};

        test_initialise()?;
        let mink = Minkowski {}.new_rep(GS.dim);
        let mu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let nu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let trace = spenso::trace!(
            Bispinor {}.new_rep(4).to_symbolic([]),
            idenso::gamma!(&mu),
            idenso::gamma!(&nu)
        );
        let spectator = idenso::color_t!(
            ColorAdjoint {}
                .new_rep(8)
                .to_symbolic([Aind::Hedge(70, 1).to_atom()]),
            ColorFundamental {}
                .new_rep(3)
                .to_symbolic([Aind::Hedge(71, 1).to_atom()]),
            ColorAntiFundamental {}
                .new_rep(3)
                .to_symbolic([Aind::Hedge(72, 1).to_atom()])
        );
        let source = trace.pow(2) * &spectator;
        let reduced = RStarIntegrand::simplify(&source)?;
        assert!(!reduced.has_repeated_lorentz_indices());
        // The early Taylor boundary retains symbolic D and the open color
        // generator; only its two independent Dirac traces are evaluated.
        let expected = Atom::num(16) * Atom::var(GS.dim) * spectator;
        super::super::tests::assert_same_counterterm(&reduced, &expected)?;
        Ok(())
    }

    #[test]
    fn finish_tensor_two_traces_contract_shared_lorentz_indices_before_poles() -> Result<()> {
        use crate::numerator::aind::Aind;
        use idenso::{IndexTooling, representations::Bispinor};
        use spenso::structure::slot::{DummyAind, ParseableAind};

        test_initialise()?;
        let epsilon = Atom::var(GS.dim_epsilon);
        let mink = Minkowski {}.new_rep(GS.dim);
        let spin = Bispinor {}.new_rep(4);
        let ports: [_; 4] =
            std::array::from_fn(|i| spin.to_symbolic([Aind::Hedge((i + 3) as u16, 1).to_atom()]));
        let mu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let nu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let gamma =
            |a: &Atom, b: &Atom, index: &Atom| function!(idenso::dirac::AGS.gamma, a, b, index);
        let source = gamma(&ports[0], &ports[1], &mu)
            * gamma(&ports[1], &ports[0], &nu)
            * gamma(&ports[2], &ports[3], &mu)
            * gamma(&ports[3], &ports[2], &nu)
            / epsilon.pow(2);
        let reduced = RStarIntegrand::finish_tensor(source, 9954)?;
        assert!(!reduced.has_repeated_lorentz_indices());
        let poles = reduced.series(GS.dim_epsilon, 0, -1)?.to_atom();
        let expected = Atom::num(64) / epsilon.pow(2) - Atom::num(32) / epsilon;
        super::super::tests::assert_same_counterterm(&poles, &expected)?;
        Ok(())
    }

    #[test]
    fn finish_tensor_reduces_clifford_words_before_caller_poles() -> Result<()> {
        use crate::numerator::aind::Aind;
        use idenso::representations::Bispinor;
        use spenso::structure::slot::{DummyAind, ParseableAind};

        test_initialise()?;
        let epsilon = Atom::var(GS.dim_epsilon);
        let mink = Minkowski {}.new_rep(GS.dim);
        let spin = Bispinor {}.new_rep(4);
        let a = spin.to_symbolic([Aind::Hedge(3, 1).to_atom()]);
        let b = spin.to_symbolic([Aind::Hedge(4, 1).to_atom()]);
        let mu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let momentum = function!(GS.emr_mom, 7, mink.to_symbolic([]));
        let gamma = |argument: &Atom| {
            function!(
                idenso::dirac::AGS.gamma,
                SPENSO_TAG.chain_in,
                SPENSO_TAG.chain_out,
                argument
            )
        };
        let word = spenso::chain!(&a, &b, gamma(&mu), gamma(&momentum), gamma(&mu));
        let source = word / epsilon.pow(2);
        let actual = RStarIntegrand::finish_tensor(source, 9951)?
            .series(GS.dim_epsilon, 0, -1)?
            .to_atom();
        let expected = (Atom::num(-2) + Atom::num(2) * &epsilon)
            * function!(idenso::dirac::AGS.gamma, &a, &b, &momentum)
            / epsilon.pow(2);
        super::super::tests::assert_same_counterterm(&actual, &expected)?;
        Ok(())
    }

    #[test]
    fn finish_tensor_preserves_external_scalar_denominators() -> Result<()> {
        use crate::uv::rstar::tests::assert_same_counterterm;

        test_initialise()?;
        let minkowski = Minkowski {}.new_rep(GS.dim);
        let mu = minkowski.to_symbolic([Atom::var(symbol!("rstar_scalar_exit::mu"))]);
        let nu = minkowski.to_symbolic([Atom::var(symbol!("rstar_scalar_exit::nu"))]);
        let spectator = function!(ETS.metric, &mu, &nu);
        let source = |a, b| function!(S.dot, function!(S.p, a), function!(S.p, b));
        let restored = |a, b| {
            function!(
                SPENSO_TAG.dot,
                function!(GS.emr_mom, a, minkowski.to_symbolic([])),
                function!(GS.emr_mom, b, minkowski.to_symbolic([]))
            )
        };
        // Cover a diagonal denominator, a mixed external invariant and a
        // reciprocal polynomial, with a retained open tensor boundary.
        for (input, expected) in [
            (source(0, 0).pow(-1), restored(0, 0).pow(-1)),
            (source(0, 1).pow(-2), restored(0, 1).pow(-2)),
            (
                (source(0, 0) + source(0, 1)).pow(-2),
                (restored(0, 0) + restored(0, 1)).pow(-2),
            ),
        ] {
            let actual = RStarIntegrand::finish_tensor(input * &spectator, 9911)?;
            assert_same_counterterm(&actual, &(expected * &spectator))?;
        }
        // Scalar restoration must not conceal unsupported inverse tensors
        // or an inverse product carrying unintegrated loop momentum.
        for invalid in [
            function!(S.p, 0, 77).pow(-1),
            function!(S.dot, function!(S.k, 1), function!(S.p, 0)).pow(-1),
        ] {
            assert!(RStarIntegrand::finish_tensor(invalid, 9912).is_err());
        }
        Ok(())
    }

    #[test]
    fn quadratic_gluon_hard_preparation_has_no_uv_pole_at_fixed_mass() -> Result<()> {
        use crate::uv::rstar::tests::{gluon_self_energy, parallel_pair, settings};

        test_initialise()?;
        let model = load_generic_model("sm");
        let graph =
            gluon_self_energy(&model, 2, 4, "V_36", |graph| parallel_pair(graph).is_some())?;
        let bubble = parallel_pair(&graph).unwrap();
        let settings = settings();
        let vakint = Vakint::new()?;
        let marked = bubble[0];
        let plan = RStarPlan::new(&graph, &model, marked, &settings)?;
        let correction = plan
            .ir
            .iter()
            .find(|correction| correction.taylor_order == 2 && graph.n_loops(&correction.hard) == 1)
            .ok_or_else(|| eyre!("fixture has no quadratic one-loop hard component"))?;
        let bare = correction
            .preparation
            .iter()
            .find(|term| term.counterterms.is_empty())
            .ok_or_else(|| eyre!("hard preparation is missing its bare term"))?;
        let mass_squared = symbol!("rstar_fixed_hard_mass_squared"; Scalar);
        let children = plan.evaluate_uv_children(
            &graph,
            &model,
            &Atom::var(mass_squared).sqrt(),
            &vakint,
            &settings,
        )?;
        let hard = RStarIntegrand::new(
            &graph,
            &model,
            &bare.remainder,
            &[],
            marked,
            Atom::var(mass_squared),
        )?;
        let polynomial = hard.taylor(correction.taylor_order as usize)?;
        let input = to_vakint_integrand(
            &polynomial,
            &hard.graph.underlying,
            &hard.region.remaining,
            &hard.graph.empty_subgraph::<SuBitGraph>(),
            &settings.vakint,
            false,
        )?;
        let mut integration_settings = settings.vakint.true_settings();
        integration_settings.number_of_terms_in_epsilon_expansion = 4;
        let integrated = RStarIntegrand::finish_tensor(
            integrate_single_mass(&vakint, &integration_settings, &input)?,
            9851,
        )?;
        let counterterm = children
            .get(&correction.hard)
            .ok_or_else(|| eyre!("hard component has no completed massive UV counterterm"))?;
        let prepared = RStarIntegrand::finish_tensor(&integrated + counterterm, 9852)?;
        let poles = prepared.series(GS.dim_epsilon, 0, -1)?.to_atom();
        debug_tags!(#uv, #integrated, #rstar, #trace;
            stage = "rstar_fixed_mass_hard_preparation",
            log.bare = integrated,
            log.counterterm = counterterm,
            log.poles = poles,
            "Checking hard preparation before any quotient UV map or outer Taylor operation"
        );
        assert!(
            poles.together().is_zero(),
            "fixed-mass hard preparation retained UV poles: {poles}"
        );
        Ok(())
    }

    /// Hard Taylor before contraction against the contracted reference, on
    /// every hard preparation term of gluonic two- and three-loop graphs,
    /// bare and with completed counterterm insertions, at every order up to
    /// the one the plan uses.
    #[test]
    fn hard_taylor_before_contraction_matches_the_contracted_path() -> Result<()> {
        use crate::uv::rstar::tests::{gluon_self_energy, parallel_pair, settings};

        test_initialise()?;
        let model = load_generic_model("sm");
        let settings = settings();
        let vakint = Vakint::new()?;
        let mass_squared = symbol!("rstar_hard_check_mass_squared"; Scalar);
        let graphs = [
            gluon_self_energy(&model, 2, 4, "V_36", |graph| parallel_pair(graph).is_some())?,
            gluon_self_energy(&model, 3, 6, "V_36", |graph| parallel_pair(graph).is_some())?,
        ];
        let mut compared = 0;
        let mut with_insertions = 0;
        for graph in &graphs {
            let marked = parallel_pair(graph).unwrap()[0];
            let plan = RStarPlan::new(graph, &model, marked, &settings)?;
            let children = plan.evaluate_uv_children(
                graph,
                &model,
                &Atom::var(mass_squared).sqrt(),
                &vakint,
                &settings,
            )?;
            for correction in &plan.ir {
                for term in &correction.preparation {
                    if term.remainder.is_scaleless(graph, plan.marked_edge) {
                        continue;
                    }
                    let Some(insertions) = term
                        .counterterms
                        .iter()
                        .map(|child| children.get(&child.argument.support).cloned())
                        .collect::<Option<Vec<_>>>()
                    else {
                        continue;
                    };
                    let hard = RStarIntegrand::new(
                        graph,
                        &model,
                        &term.remainder,
                        &insertions,
                        plan.marked_edge,
                        Atom::var(mass_squared),
                    )?;
                    for degree in 0..=correction.taylor_order as usize {
                        let new = hard.taylor(degree)?;
                        let old = hard.taylor_contracted(degree)?;
                        assert!(
                            (&new - &old).expand().is_zero(),
                            "hard Taylor order {degree} differs:\nnew: {new}\nold: {old}"
                        );
                        compared += 1;
                    }
                    with_insertions += usize::from(!insertions.is_empty());
                }
            }
        }
        assert!(
            compared > 0 && with_insertions > 0,
            "{compared} {with_insertions}"
        );
        Ok(())
    }

    #[test]
    fn hard_taylor_preserves_the_single_mass_squared_in_the_existing_bridge() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("scalars");
        let mut graph: Graph = finalized_runtime_dot!(
            digraph bubble {
                graph [projector=1];
                node [num=1];
                edge [particle="scalar_1" num=1];
                e [style=invis];
                e -> A:0 [id=2 sink="{ufo_order:0}"];
                B:1 -> e [id=3 source="{ufo_order:0}"];
                A -> B [id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
                A -> B [id=1 lmb_id=0 source="{ufo_order:2}" sink="{ufo_order:2}"];
            }, &model
        )?;
        let edges = graph
            .iter_edges()
            .map(|(_, edge, _)| edge)
            .collect::<Vec<_>>();
        for edge in edges {
            graph.underlying[edge].particle =
                graph[edge].particle.clone().override_mass(Some(Atom::Zero));
        }
        let settings = UVgenerationSettings::default();
        let marked = EdgeIndex::from(0);
        let plan = RStarPlan::new(&graph, &model, marked, &settings)?;
        let mass_squared = symbol!("rstar_bridge_mass_squared");
        let input = RStarIntegrand::new(
            &graph,
            &model,
            &plan.ordinary_uv[0].remainder,
            &[],
            marked,
            Atom::var(mass_squared),
        )?;
        assert!(!input.soft_edges().is_empty());
        // A positive Taylor order differentiates shifted Lorentz products.
        // Its parameter must be scalar so it factors out of those products
        // before Symbolica computes the series.
        let quadratic = input.taylor(2)?;
        assert!(input.soft_edges().iter().any(|edge| {
            quadratic
                .replace(function!(GS.emr_mom, usize::from(*edge), W_.x___))
                .match_iter()
                .next()
                .is_some()
        }));
        let hard = input.taylor(0)?;
        for edge in input.soft_edges() {
            assert!(
                hard.replace(function!(GS.emr_mom, usize::from(*edge), W_.x___))
                    .match_iter()
                    .next()
                    .is_none()
            );
        }
        let converted = to_vakint_integrand(
            &hard,
            &input.graph.underlying,
            &input.region.remaining,
            &input.graph.empty_subgraph::<SuBitGraph>(),
            &settings.vakint,
            false,
        )?;
        assert_eq!(converted.0.len(), 1);
        let mut masses = BTreeSet::new();
        converted.0[0].integral.visitor(&mut |part| {
            if let AtomView::Fun(prop) = part
                && prop.get_symbol() == S.prop
            {
                masses.insert(prop.get(3).to_owned());
                return false;
            }
            true
        });
        assert_eq!(
            masses,
            BTreeSet::from([Atom::Zero, Atom::var(mass_squared)])
        );
        // Boundary sources are retained individually, including the edge that
        // is dependent in the internal momentum routing. They must never be
        // replaced by the zero signatures of an internal-subgraph basis.
        for edge in input.soft_edges() {
            let vector = function!(
                GS.emr_mom,
                usize::from(*edge),
                Minkowski {}.new_rep(GS.dim).to_symbolic([])
            );
            let factor = function!(SPENSO_TAG.dot, &vector, &vector);
            let mut decorated = graph.clone();
            decorated.underlying[marked].num = crate::graph::Autogen::explicit(factor.clone());
            let preserved = RStarIntegrand::new(
                &decorated,
                &model,
                &input.region,
                &[],
                marked,
                Atom::var(mass_squared),
            )?;
            assert!(
                (preserved.expression()? - input.expression()? * &factor)
                    .together()
                    .is_zero()
            );
        }
        Ok(())
    }

    #[test]
    fn tensor_finish_restores_dots_and_regulates_only_scalar_dimension() -> Result<()> {
        test_initialise()?;
        let raw = function!(S.dot, function!(S.p, 7), function!(S.p, 7)) * function!(S.g, 1, 1);
        let actual = RStarIntegrand::simplify(&RStarIntegrand::finish_tensor(raw, 543)?)?;
        let vector = function!(GS.emr_mom, 7, Minkowski {}.new_rep(GS.dim).to_symbolic([]));
        let expected = (Atom::num(4) - Atom::num(2) * Atom::var(GS.dim_epsilon))
            * function!(SPENSO_TAG.dot, &vector, &vector);
        assert!(
            (actual - expected).expand().is_zero(),
            "Vakint spectators must retain the original native Lorentz representation"
        );
        Ok(())
    }

    #[test]
    fn physical_remainder_preserves_scale_and_promotes_insertion_poles() -> Result<()> {
        use crate::numerator::aind::Aind;
        use spenso::structure::slot::ParseableAind;

        test_initialise()?;
        let model = load_generic_model("scalars");
        let mut graph: Graph = finalized_runtime_dot!(
            digraph bubble {
                graph [projector=1];
                node [num=1];
                edge [particle="scalar_1" num=1];
                e [style=invis];
                e -> A:0 [id=2 sink="{ufo_order:0}"];
                B:1 -> e [id=3 source="{ufo_order:0}"];
                A -> B [id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
                A -> B [id=1 lmb_id=0 source="{ufo_order:2}" sink="{ufo_order:2}"];
            }, &model
        )?;
        let edges = graph
            .iter_edges()
            .map(|(_, edge, _)| edge)
            .collect::<Vec<_>>();
        for edge in edges {
            graph.underlying[edge].particle =
                graph[edge].particle.clone().override_mass(Some(Atom::Zero));
        }
        let settings = super::super::tests::settings();
        let marked = EdgeIndex::from(0);
        let plan = RStarPlan::new(&graph, &model, marked, &settings)?;
        let mut input = RStarIntegrand::new(
            &graph,
            &model,
            &plan.ordinary_uv[0].remainder,
            &[],
            marked,
            Atom::Zero,
        )?;
        let epsilon = Atom::var(GS.dim_epsilon);
        *input.expression_mut()? /= epsilon.pow(2);
        let mut integration = settings.vakint.true_settings();
        integration.integral_normalization_factor = vakint::LoopNormalizationFactor::GScheme;
        integration.number_of_terms_in_epsilon_expansion = 1;
        let vakint = Vakint::new()?;
        let actual = input.integrate_physical(&vakint, &settings, &integration, 5431)?;
        let actual = super::super::tests::canonical_tensor(&actual)?;
        let external = function!(
            GS.emr_mom,
            usize::from(input.soft_edges[0]),
            Minkowski {}.new_rep(GS.dim).to_symbolic([])
        );
        let logarithm = function!(
            Symbol::LOG,
            -function!(SPENSO_TAG.dot, &external, &external) / Atom::var(GS.mu_r_sq)
        );
        let expected = epsilon.pow(-3) - &logarithm / epsilon.pow(2)
            + logarithm.pow(2) / (Atom::num(2) * &epsilon);
        assert!(
            (&actual - &expected).expand().is_zero(),
            "physical remainders must retain scale logarithms and positive orders promoted by poles: actual={actual}, expected={expected}"
        );
        // A scalar loop-external product can be restored as a compact metric
        // by finish_tensor. It must reach Forcer as a dot, while an open
        // Lorentz metric and a reciprocal external invariant stay spectators.
        let left = Minkowski {}
            .new_rep(GS.dim)
            .to_symbolic([Aind::Hedge(87, 0).to_atom()]);
        let right = Minkowski {}
            .new_rep(GS.dim)
            .to_symbolic([Aind::Hedge(88, 0).to_atom()]);
        let metric = function!(ETS.metric, &left, &right);
        let loop_momentum = function!(
            GS.emr_mom,
            usize::from(input.basis.loop_edges.raw[0]),
            Minkowski {}.new_rep(GS.dim).to_symbolic([])
        );
        let original = input.expression()?.clone();
        *input.expression_mut()? *= function!(SPENSO_TAG.dot, &loop_momentum, &external)
            / function!(SPENSO_TAG.dot, &external, &external)
            * &metric;
        let actual_scalar = input.integrate_physical(&vakint, &settings, &integration, 5433)?;
        let actual_scalar = super::super::tests::canonical_tensor(&actual_scalar)?;
        let expected_scalar = &expected * metric / Atom::num(2);
        assert!(
            (&actual_scalar - &expected_scalar).expand().is_zero(),
            "physical scalar products and their tensor spectators must survive ingress: actual={actual_scalar}, expected={expected_scalar}"
        );
        *input.expression_mut()? = original;
        // A vacuum projector would wrongly discard this odd loop tensor.
        // With the physical external momentum, symmetry of the two lines
        // instead gives integral k^mu = p^mu/2 times the scalar bubble.
        let slot = Minkowski {}
            .new_rep(GS.dim)
            .to_symbolic([Aind::Hedge(86, 0).to_atom()]);
        let loop_vector = function!(
            GS.emr_mom,
            usize::from(input.basis.loop_edges.raw[0]),
            &slot
        );
        *input.expression_mut()? *= loop_vector;
        let actual_tensor = input.integrate_physical(&vakint, &settings, &integration, 5432)?;
        let actual_tensor = super::super::tests::canonical_tensor(&actual_tensor)?;
        let expected_tensor = expected
            * function!(GS.emr_mom, usize::from(input.soft_edges[0]), &slot)
            / Atom::num(2);
        assert!(
            (&actual_tensor - &expected_tensor).expand().is_zero(),
            "the physical tensor projector must retain the external direction: actual={actual_tensor}, expected={expected_tensor}"
        );
        Ok(())
    }
}

/// The edge-basis renames ([`rename_edge_basis`], [`rename_kept_edges`]) with
/// and without `flatten` (`IDENSO_EDGE_OUTPUT=flatten`).
#[cfg(test)]
mod edge_rename_tests {
    use std::collections::{BTreeMap, BTreeSet};

    use symbolica::{
        atom::{Atom, AtomCore, AtomView},
        function, symbol,
    };

    use super::{
        EDGE_BASIS_OFFSET, EDGE_OFFSET, flat_combination, rename_edge_basis, rename_kept_edges,
    };
    use crate::utils::GS;

    type Combinations = BTreeMap<Vec<(i64, Atom)>, i64>;

    /// `Q(edge, slot)`.
    fn q(edge: i64, slot: &str) -> Atom {
        function!(GS.emr_mom, Atom::num(edge), Atom::var(symbol!(slot)))
    }

    fn mu() -> &'static str {
        "edge_rename_test::mu"
    }

    fn nu() -> &'static str {
        "edge_rename_test::nu"
    }

    fn renamed(expression: &Atom, flatten: bool) -> (Atom, Combinations) {
        let mut combinations = Combinations::new();
        let renamed = rename_edge_basis(expression, flatten, &mut combinations);
        (renamed, combinations)
    }

    fn combination(terms: &[(i64, i64)], id: i64) -> Combinations {
        let key = terms
            .iter()
            .map(|(edge, c)| (*edge, Atom::num(*c)))
            .collect();
        BTreeMap::from([(key, id)])
    }

    /// `x - (y)` with `y` a sum stays nested (no automatic distribution).
    fn is_nested(expression: &Atom) -> bool {
        let AtomView::Add(sum) = expression.as_view() else {
            return false;
        };
        sum.iter().any(|term| matches!(term, AtomView::Mul(_)))
    }

    #[test]
    fn flatten_renames_a_nested_combination_as_one_vector() {
        // Q(2,μ) - (Q(2,μ) + Q(3,μ) - Q(6,μ)) = Q(6,μ) - Q(3,μ).
        let nested = q(2, mu()) - (q(2, mu()) + q(3, mu()) - q(6, mu()));
        assert!(is_nested(&nested), "{nested}");
        let (flat, flat_combinations) = renamed(&nested, true);
        assert_eq!(flat, q(EDGE_BASIS_OFFSET, mu()));
        assert_eq!(
            flat_combinations,
            combination(&[(3, -1), (6, 1)], EDGE_BASIS_OFFSET)
        );
        // Unset: only the inner sum is renamed, the outer one stays.
        let (inner, inner_combinations) = renamed(&nested, false);
        assert_eq!(inner, q(2, mu()) - q(EDGE_BASIS_OFFSET, mu()));
        assert_eq!(
            inner_combinations,
            combination(&[(2, 1), (3, 1), (6, -1)], EDGE_BASIS_OFFSET)
        );
        // Both are the same vector once the combinations are put back.
        for (renamed, combinations) in [(flat, flat_combinations), (inner, inner_combinations)] {
            let (terms, id) = combinations.into_iter().next().unwrap();
            let image = Atom::add_many(
                terms
                    .iter()
                    .map(|(edge, c)| c * q(*edge, mu()))
                    .collect::<Vec<_>>()
                    .iter(),
            );
            let back = renamed.replace(q(id, mu())).with(image);
            assert!((&back - &nested).expand().is_zero(), "{back}");
        }
    }

    #[test]
    fn flatten_cancellation_leaves_one_term_or_zero() {
        // Q(2,μ) - (Q(2,μ) - 2 Q(3,μ)) = 2 Q(3,μ): one term, no combination.
        let single = q(2, mu()) - (q(2, mu()) - Atom::num(2) * q(3, mu()));
        assert!(is_nested(&single), "{single}");
        let (flat, combinations) = renamed(&single, true);
        assert_eq!(flat, Atom::num(2) * q(3, mu()));
        assert!(combinations.is_empty());
        let (inner, combinations) = renamed(&single, false);
        assert_eq!(inner, q(2, mu()) - q(EDGE_BASIS_OFFSET, mu()));
        assert_eq!(
            combinations,
            combination(&[(2, 1), (3, -2)], EDGE_BASIS_OFFSET)
        );

        // Q(2,μ) + Q(3,μ) - (Q(2,μ) + Q(3,μ)) = 0.
        let zero = q(2, mu()) + q(3, mu()) - (q(2, mu()) + q(3, mu()));
        assert!(is_nested(&zero), "{zero}");
        let (flat, combinations) = renamed(&zero, true);
        assert!(flat.is_zero(), "{flat}");
        assert!(combinations.is_empty());
        // Inside a product the cancelled slot makes the product zero.
        let product = function!(symbol!("edge_rename_test::f"), q(4, nu())) * &zero;
        let (flat, combinations) = renamed(&product, true);
        assert!(flat.is_zero(), "{flat}");
        assert!(combinations.is_empty());
    }

    #[test]
    fn flatten_with_mixed_tails_falls_back_to_the_inner_rename() {
        // Q(2,ν) - (Q(3,μ) + Q(4,μ)): the expanded sum has two tails, so the
        // outer sum is not renamed and the inner one is, as without flatten.
        let mixed = q(2, nu()) - (q(3, mu()) + q(4, mu()));
        assert!(is_nested(&mixed), "{mixed}");
        let with = renamed(&mixed, true);
        let without = renamed(&mixed, false);
        assert_eq!(with.0, q(2, nu()) - q(EDGE_BASIS_OFFSET, mu()));
        assert_eq!(with.1, combination(&[(3, 1), (4, 1)], EDGE_BASIS_OFFSET));
        assert_eq!(with, without);
    }

    #[test]
    fn flatten_rejects_a_non_numeric_coefficient() {
        let vector =
            |leaf: AtomView<'_>| matches!(leaf, AtomView::Fun(f) if f.get_symbol() == GS.emr_mom);
        let x = Atom::var(symbol!("edge_rename_test::x"));
        // Q(2,μ) + x (Q(2,μ) + Q(3,μ)) is not a numeric combination of vectors.
        let symbolic = q(2, mu()) + &x * (q(2, mu()) + q(3, mu()));
        assert!(is_nested(&symbolic), "{symbolic}");
        assert_eq!(flat_combination(symbolic.as_view(), true, &vector), None);
        let with = renamed(&symbolic, true);
        assert_eq!(with.0, q(2, mu()) + &x * q(EDGE_BASIS_OFFSET, mu()));
        assert_eq!(with.1, combination(&[(2, 1), (3, 1)], EDGE_BASIS_OFFSET));
        assert_eq!(with, renamed(&symbolic, false));
        // The numeric multiple is flattened; without the flag nothing is.
        let numeric = q(2, mu()) + Atom::num(3) * (q(2, mu()) + q(3, mu()));
        assert_eq!(
            flat_combination(numeric.as_view(), true, &vector),
            Some(Atom::num(4) * q(2, mu()) + Atom::num(3) * q(3, mu()))
        );
        assert_eq!(flat_combination(numeric.as_view(), false, &vector), None);
        // A sum of vectors only is not nested: nothing to flatten.
        let plain = q(2, mu()) + q(3, mu());
        assert_eq!(flat_combination(plain.as_view(), true, &vector), None);
    }

    #[test]
    fn kept_edges_rename_flattens_numerators_and_leaves_denominators() {
        let loop_edges = BTreeSet::from([2, 3, 6]);
        let nested = |slot| q(2, slot) - (q(2, slot) + q(3, slot) - q(6, slot));
        // A denominator holding a renamable (and a nested) combination.
        let den = function!(
            GS.den,
            Atom::num(7),
            q(3, mu()) + q(6, mu()),
            nested(mu()),
            Atom::num(0)
        );
        let expression = &den * nested(nu()) * (q(2, nu()) + q(3, nu()) + q(9, nu()));
        for flatten in [true, false] {
            let mut combinations = Combinations::new();
            let renamed = rename_kept_edges(&expression, &loop_edges, flatten, &mut combinations);
            let key = |terms: &[(i64, i64)]| -> Vec<(i64, Atom)> {
                terms
                    .iter()
                    .map(|(edge, c)| (*edge, Atom::num(*c)))
                    .collect()
            };
            // With flatten -Q(3,ν) + Q(6,ν) is one kept vector; without, only
            // the inner Q(2) + Q(3) - Q(6) is, and the outer Q(2) stays.
            let first_key = if flatten {
                key(&[(3, -1), (6, 1)])
            } else {
                key(&[(2, 1), (3, 1), (6, -1)])
            };
            let second_key = key(&[(2, 1), (3, 1)]);
            assert_eq!(
                combinations.keys().cloned().collect::<BTreeSet<_>>(),
                BTreeSet::from([first_key.clone(), second_key.clone()]),
                "flatten={flatten}"
            );
            let mut ids = combinations.values().copied().collect::<Vec<_>>();
            ids.sort();
            assert_eq!(ids, vec![EDGE_OFFSET, EDGE_OFFSET + 1]);
            let first = q(combinations[&first_key], nu());
            let first = if flatten { first } else { q(2, nu()) - first };
            // Q(9) is not a loop edge of the component: it stays explicit.
            let second = q(9, nu()) + q(combinations[&second_key], nu());
            let expected = &den * first * second;
            assert_eq!(renamed, expected, "flatten={flatten}");
        }
    }
}
