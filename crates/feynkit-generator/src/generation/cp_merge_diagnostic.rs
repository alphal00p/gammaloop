//! Bounded, opt-in scalar CP preparation observation through the ordinary owners.
//! Compiled only in library tests; no production tracing or serialization API.

use super::{ColoredTopology, DiagramInput, Endpoint, GeneratedGraph};
use serde_json::{Value, json};
use std::cell::RefCell;
use std::thread_local;

thread_local! {
    static EVENTS: RefCell<Option<Vec<Value>>> = const { RefCell::new(None) };
}

pub(crate) fn record(event: impl FnOnce() -> Value) {
    EVENTS.with_borrow_mut(|events| {
        if let Some(events) = events {
            events.push(event());
        }
    });
}

pub(super) fn start() {
    EVENTS.with_borrow_mut(|events| *events = Some(Vec::new()));
}

pub(super) fn finish() -> Vec<Value> {
    EVENTS.with_borrow_mut(|events| events.take().unwrap())
}

pub(super) fn graph(graph: &GeneratedGraph) -> Value {
    let endpoint = |endpoint| match endpoint {
        Endpoint::Vertex(vertex) => json!({"vertex": vertex}),
        Endpoint::External(external) => json!({"external": {
            "index": external.index, "state": external.state,
            "particle": external.particle.index(),
        }}),
    };
    json!({
        "rules": graph.vertices().map(|rule| rule.index()).collect::<Vec<_>>(),
        "edges": graph.edges().map(|edge| json!({
            "source": endpoint(edge.source), "sink": endpoint(edge.sink),
            "particle": edge.particle.index(), "directed": edge.directed,
        })).collect::<Vec<_>>(),
    })
}

pub(super) fn input(input: &DiagramInput) -> Value {
    let topology = |topology: &ColoredTopology| {
        json!({
            "graph": graph(&topology.graph), "aut": topology.symmetry,
            "cuts": topology.cut_partitions,
            "thresholds": topology.topology_threshold_candidates,
        })
    };
    json!({
        "representative": topology(&input.representative),
        "comparison": topology(&input.comparison),
        "signs": {
            "internal_loop_negative": input.signs.internal_loop_negative,
            "include_external_ordering": input.signs.include_external_ordering,
            "external_ordering_negative": input.signs.external_ordering_negative,
            "include_antifermion_spin_sum": input.signs.include_antifermion_spin_sum,
            "antifermion_spin_sum_negative": input.signs.antifermion_spin_sum_negative,
        },
        "reversed_edges": input.reversed_edges, "multiplicity": input.multiplicity,
    })
}
