//! Diagnostics (`GAMMALOOP_EDGE_EXPERIMENT=1` with `VAKINT_LOWERING_CENSUS`):
//! the Feynman numerator of an R* remainder taken through the numerator
//! algebra (index contraction, Dirac algebra, dot notation, expansion) twice
//! on the same input: in the loop-momentum basis, as `RStarIntegrand::new`
//! routes it today, and with every edge momentum `Q(e)` kept as an
//! independent vector, the routing applied to the scalar result afterwards.
//! Records time, terms, bytes, distinct scalar products, peak memory and the
//! graph class of every edge scalar product. Changes nothing.

use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use linnet::half_edge::involution::HedgePair;
use linnet::half_edge::subgraph::SuBitGraph;
use spenso::network::tags::SPENSO_TAG;
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::id::Replacement;

use crate::graph::Graph;
use crate::utils::GS;

pub(super) fn enabled() -> bool {
    static ON: std::sync::LazyLock<bool> =
        std::sync::LazyLock::new(|| std::env::var("GAMMALOOP_EDGE_EXPERIMENT").is_ok());
    *ON
}

fn status_kib(field: &str) -> u64 {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|status| {
            status.lines().find_map(|line| {
                line.strip_prefix(field)
                    .and_then(|rest| rest.trim().trim_end_matches("kB").trim().parse().ok())
            })
        })
        .unwrap_or(0)
}

fn reset_peak() -> u64 {
    let _ = std::fs::write("/proc/self/clear_refs", "5");
    status_kib("VmRSS:")
}

/// Distinct scalar products (spenso or vakint `dot`) of an expression.
fn dots(expression: AtomView<'_>) -> BTreeSet<Atom> {
    let mut found = BTreeSet::new();
    expression.visitor(&mut |part| {
        if let AtomView::Fun(function) = part
            && (function.get_symbol() == SPENSO_TAG.dot || function.get_symbol() == vakint::symbols::S.dot)
        {
            found.insert(part.to_owned());
            return false;
        }
        true
    });
    found
}

fn summands(view: AtomView<'_>) -> usize {
    match view {
        AtomView::Add(sum) => sum.get_nargs(),
        _ if view.is_zero() => 0,
        _ => 1,
    }
}

/// The edge of `Q(e, ...)`.
fn edge_of(vector: AtomView<'_>) -> Option<usize> {
    let AtomView::Fun(function) = vector else {
        return None;
    };
    if function.get_symbol() != GS.emr_mom {
        return None;
    }
    usize::try_from(function.get(0)).ok()
}

struct Stage {
    seconds: f64,
    terms: usize,
    bytes: usize,
    dots: usize,
    peak_kib: u64,
}

fn run(input: &Atom, simplify: &dyn Fn(AtomView<'_>) -> eyre::Result<Atom>) -> eyre::Result<(Atom, Stage)> {
    let base = reset_peak();
    let started = Instant::now();
    let simplified = simplify(input.as_view())?;
    let expanded = simplified.expand();
    let seconds = started.elapsed().as_secs_f64();
    let peak_kib = status_kib("VmHWM:").saturating_sub(base);
    let stage = Stage {
        seconds,
        terms: summands(expanded.as_view()),
        bytes: expanded.as_view().get_byte_size(),
        dots: dots(expanded.as_view()).len(),
        peak_kib,
    };
    Ok((expanded, stage))
}

fn json(stage: &Stage) -> serde_json::Value {
    serde_json::json!({
        "seconds": stage.seconds,
        "terms": stage.terms,
        "bytes": stage.bytes,
        "dots": stage.dots,
        "peak_kib": stage.peak_kib,
    })
}

/// Run the comparison on `unrouted` (numerator times insertions, before
/// `RStarIntegrand::new` routes it) and record it.
pub(super) fn compare(
    graph: &Graph,
    converter: &Graph,
    remaining: &SuBitGraph,
    unrouted: &Atom,
    routing: &[Replacement],
    loops: usize,
    simplify: &dyn Fn(AtomView<'_>) -> eyre::Result<Atom>,
) {
    let outcome = (|| -> eyre::Result<serde_json::Value> {
        let routing_started = Instant::now();
        let routed = unrouted.replace_multiple(routing);
        let routing_seconds = routing_started.elapsed().as_secs_f64();
        let (loop_basis, loop_stage) = run(&routed, simplify)?;
        drop(routed);
        let (edge_basis, edge_stage) = run(unrouted, simplify)?;
        // The deferred basis change: route the scalar edge polynomial.
        let base = reset_peak();
        let started = Instant::now();
        let rerouted = edge_basis.replace_multiple(routing).expand();
        let reroute_seconds = started.elapsed().as_secs_f64();
        let reroute_peak = status_kib("VmHWM:").saturating_sub(base);
        let identical = rerouted == loop_basis;
        let equal = identical || (&rerouted - &loop_basis).expand().is_zero();

        // Classes of the edge scalar products of the edge-basis result.
        let underlying = &converter.underlying;
        let mut ends = BTreeMap::<usize, (usize, usize)>::new();
        for (pair, edge, _) in underlying.iter_edges_of(remaining) {
            if let HedgePair::Paired { source, sink } = pair {
                ends.insert(
                    usize::from(edge),
                    (usize::from(underlying.node_id(source)), usize::from(underlying.node_id(sink))),
                );
            }
        }
        let mut degree = BTreeMap::<usize, usize>::new();
        for (pair, _, _) in underlying.iter_edges() {
            match pair {
                HedgePair::Paired { source, sink } => {
                    *degree.entry(usize::from(underlying.node_id(source))).or_default() += 1;
                    *degree.entry(usize::from(underlying.node_id(sink))).or_default() += 1;
                }
                HedgePair::Unpaired { hedge, .. } => {
                    *degree.entry(usize::from(underlying.node_id(hedge))).or_default() += 1;
                }
                _ => {}
            }
        }
        let mut classes = BTreeMap::<&str, usize>::new();
        for dot in dots(edge_basis.as_view()) {
            let AtomView::Fun(function) = dot.as_view() else { continue };
            let class = if function.get_nargs() != 2 {
                "other"
            } else {
                match (edge_of(function.get(0)), edge_of(function.get(1))) {
                    (Some(a), Some(b)) => match (ends.get(&a), ends.get(&b)) {
                        (Some(_), Some(_)) if a == b => "same internal edge (inverse propagator)",
                        (Some(&(a0, a1)), Some(&(b0, b1))) => {
                            let shared = [a0, a1].into_iter().find(|n| *n == b0 || *n == b1);
                            match shared {
                                Some(node) if degree.get(&node) == Some(&3) => {
                                    "adjacent at a 3-point vertex (three inverse propagators)"
                                }
                                Some(_) => "adjacent at a higher vertex",
                                None => "non-adjacent internal edges (needs ISPs)",
                            }
                        }
                        _ if a == b => "same boundary/external edge",
                        _ => "involves a boundary/external edge",
                    },
                    _ => "non-edge vector",
                }
            };
            *classes.entry(class).or_default() += 1;
        }
        let _ = graph;
        Ok(serde_json::json!({
            "loops": loops,
            "input_bytes": unrouted.as_view().get_byte_size(),
            "input_terms": summands(unrouted.as_view()),
            "routing_seconds": routing_seconds,
            "loop_basis": json(&loop_stage),
            "edge_basis": json(&edge_stage),
            "reroute_seconds": reroute_seconds,
            "reroute_peak_kib": reroute_peak,
            "rerouted_terms": summands(rerouted.as_view()),
            "identical": identical,
            "equal": equal,
            "edge_dot_classes": classes,
        }))
    })();
    match outcome {
        Ok(record) => vakint::census::record("edge_experiment", record),
        Err(error) => vakint::census::record(
            "edge_experiment",
            serde_json::json!({ "error": error.to_string() }),
        ),
    }
}
