//! Native topology costs; setup and graph cloning are outside mutation timings.
use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use linnet::{
    half_edge::{
        builder::HedgeGraphBuilder,
        involution::Hedge,
        nodestore::{NodeForest, NodeStorageOps},
        swap::Swap,
        HedgeGraph, NoData, NodeIndex,
    },
    permutation::Permutation,
    tree::{child_pointer::ParentChildStore, child_vec::ChildVecStore, Forest},
};

type Graph<S> = HedgeGraph<usize, usize, NoData, S>;

fn build<S: NodeStorageOps<NodeData = usize>>(nodes: usize, matching: bool) -> Graph<S> {
    let mut builder = HedgeGraphBuilder::new();
    for node in 0..nodes {
        builder.add_node(node);
    }
    for node in 0..nodes {
        builder.add_edge(NodeIndex(node), NodeIndex((node + 1) % nodes), 1, false);
    }
    if matching {
        for node in 0..nodes / 2 {
            builder.add_edge(NodeIndex(node), NodeIndex(node + nodes / 2), 1, false);
        }
    }
    builder.build()
}

fn dfs<S: NodeStorageOps<NodeData = usize>>(graph: &Graph<S>) -> usize {
    let mut seen = vec![false; graph.n_nodes()];
    let mut stack = vec![NodeIndex(0)];
    seen[0] = true;
    let mut count = 0;
    while let Some(node) = stack.pop() {
        count += 1;
        for hedge in graph.iter_crown(node) {
            let next = graph.node_id(graph.inv(hedge));
            if !seen[next.0] {
                seen[next.0] = true;
                stack.push(next);
            }
        }
    }
    count
}

fn append<S: NodeStorageOps<NodeData = usize>>(mut graph: Graph<S>) -> Graph<S> {
    let n = graph.n_nodes();
    for i in 0..32 {
        graph = graph
            .add_pair(NodeIndex(i % n), NodeIndex((i + n / 2) % n), 1, false)
            .unwrap()
            .2;
    }
    graph
}

/// A square grid, with one half-edge per node to follow it through merges.
fn grid<S: NodeStorageOps<NodeData = usize>>(side: usize) -> (Graph<S>, Vec<Hedge>) {
    let mut builder = HedgeGraphBuilder::new();
    for node in 0..side * side {
        builder.add_node(node);
    }
    for node in 0..side * side {
        if node % side + 1 < side {
            builder.add_edge(NodeIndex(node), NodeIndex(node + 1), 1, false);
        }
        if node + side < side * side {
            builder.add_edge(NodeIndex(node), NodeIndex(node + side), 1, false);
        }
    }
    let graph: Graph<S> = builder.build();
    let representatives = (0..side * side)
        .map(|node| graph.iter_crown(NodeIndex(node)).next().unwrap())
        .collect();
    (graph, representatives)
}

/// Merge every node into the first one, discarding history after each batch,
/// as batched tensor-network contraction does.
fn hub_merge<S: NodeStorageOps<NodeData = usize>>(
    graph: &mut Graph<S>,
    representatives: &[Hedge],
    batch: usize,
) {
    for (merged, &hedge) in representatives[1..].iter().enumerate() {
        let hub = graph.node_id(representatives[0]);
        graph.identify_nodes(&[hub, graph.node_id(hedge)], merged);
        if (merged + 1) % batch == 0 {
            graph.forget_identification_history();
        }
    }
    graph.forget_identification_history();
}

/// A fixed pseudo-random permutation of `n` half-edges.
fn shuffle(n: usize) -> Permutation {
    let mut map: Vec<usize> = (0..n).collect();
    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
    for i in (1..n).rev() {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        map.swap(i, (state % (i as u64 + 1)) as usize);
    }
    Permutation::from_map(map)
}

fn bench_store<S: NodeStorageOps<NodeData = usize> + Clone>(c: &mut Criterion, name: &str) {
    let mut group = c.benchmark_group(name);
    for n in [12, 1000, 10000] {
        let graph: Graph<S> = build(n, true);
        assert_eq!(dfs(&graph), n);
        group.bench_with_input(BenchmarkId::new("build", n), &n, |b, &n| {
            b.iter(|| build::<S>(black_box(n), true));
        });
        group.bench_function(BenchmarkId::new("dfs", n), |b| {
            b.iter(|| dfs(black_box(&graph)));
        });
        group.bench_function(BenchmarkId::new("adjacency", n), |b| {
            b.iter(|| {
                let graph = black_box(&graph);
                graph
                    .iter_nodes()
                    .flat_map(|(_, crown, _)| crown)
                    .map(|hedge| graph.node_id(graph.inv(hedge)).0 + 1)
                    .sum::<usize>()
            });
        });
        for warmed in [false, true] {
            group.bench_function(
                BenchmarkId::new(
                    if warmed {
                        "append_spare_capacity"
                    } else {
                        "append_after_clone"
                    },
                    n,
                ),
                |b| {
                    b.iter_batched(
                        || {
                            let copy = graph.clone();
                            if warmed {
                                append(copy)
                            } else {
                                copy
                            }
                        },
                        append,
                        BatchSize::SmallInput,
                    );
                },
            );
        }
    }
    for depth in [2, 8, 32] {
        let mut history: Graph<S> = build(128 * depth, false);
        for group in 0..128 {
            let mut previous = NodeIndex(group);
            for layer in 1..depth {
                previous =
                    history.identify_nodes(&[NodeIndex(group + 128 * layer), previous], group);
            }
        }
        group.bench_function(BenchmarkId::new("discard_history", depth), |b| {
            b.iter_batched(
                || history.clone(),
                |mut graph| {
                    black_box(graph.forget_identification_history());
                    graph
                },
                BatchSize::SmallInput,
            );
        });
        history.forget_identification_history();
        history.check().unwrap();
        assert_eq!(history.n_nodes(), 128);
        assert_eq!(dfs(&history), 128);
        group.bench_function(BenchmarkId::new("dfs_after_identification", depth), |b| {
            b.iter(|| dfs(black_box(&history)));
        });
    }
    // 960 and 9,800 half-edges, forgetting after every 64 merges or only once.
    for side in [16, 50] {
        let (graph, representatives) = grid::<S>(side);
        let mut merged = graph.clone();
        hub_merge(&mut merged, &representatives, 64);
        merged.check().unwrap();
        assert_eq!(merged.n_nodes(), 1);
        for (label, batch) in [("hub_merge_batch64", 64), ("hub_merge_once", usize::MAX)] {
            group.bench_function(BenchmarkId::new(label, graph.n_hedges()), |b| {
                b.iter_batched(
                    || graph.clone(),
                    |mut graph| {
                        hub_merge(&mut graph, &representatives, batch);
                        graph
                    },
                    BatchSize::SmallInput,
                );
            });
        }
        // Half-edge swaps, which require identification history to be forgotten.
        let permutation = shuffle(graph.n_hedges());
        let mut permuted = graph.clone();
        Swap::<Hedge>::permute(&mut permuted, &permutation);
        permuted.check().unwrap();
        group.bench_function(BenchmarkId::new("permute", graph.n_hedges()), |b| {
            b.iter_batched(
                || graph.clone(),
                |mut graph| {
                    Swap::<Hedge>::permute(&mut graph, black_box(&permutation));
                    graph
                },
                BatchSize::SmallInput,
            );
        });
    }
    group.finish();
}

fn nodestores(c: &mut Criterion) {
    bench_store::<Forest<usize, ParentChildStore<()>>>(c, "parent_child");
    bench_store::<Forest<usize, ChildVecStore<()>>>(c, "child_vec");
    bench_store::<NodeForest<usize>>(c, "node_forest");
}

criterion_group!(benches, nodestores);
criterion_main!(benches);
