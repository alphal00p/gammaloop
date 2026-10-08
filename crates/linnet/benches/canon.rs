//! Native canonization against Graphica on whole generated Feynman-graph sets (externals are
//! degree-1 nodes on both sides, Linnet's input converted once outside the timing).
use std::hint::black_box;

use criterion::{criterion_group, criterion_main, Criterion};
use linnet::half_edge::HedgeGraph;
use symbolica::graph::{GenerationSettings, Graph, HalfEdge};

type Plain = Graph<usize, u8>;

fn generate(
    externals: &[(usize, HalfEdge<u8>)],
    vertices: &[Vec<HalfEdge<u8>>],
    loops: usize,
) -> Vec<Plain> {
    let min_degree = vertices.iter().map(Vec::len).min().unwrap();
    let max_vertices = (2 * loops + externals.len() - 2) / (min_degree - 2);
    let settings = GenerationSettings::new()
        .max_loops(loops)
        .max_vertices(externals.len() + max_vertices)
        .allow_self_loops(true)
        .allow_zero_flow_edges(true);
    Plain::generate(externals, vertices, settings)
        .unwrap()
        .into_keys()
        .filter(|g| g.num_loops() == loops)
        .collect()
}

fn canon(c: &mut Criterion) {
    let phi = HalfEdge::undirected(0);
    let (photon, electron) = (HalfEdge::undirected(1), HalfEdge::incoming(2));
    let sets = [
        ("phi3 vacuum L6", generate(&[], &[vec![phi; 3]], 6)),
        (
            "phi34 2-point L3",
            generate(&[(1, phi), (2, phi)], &[vec![phi; 3], vec![phi; 4]], 3),
        ),
        (
            "QED photon self-energy L3",
            generate(
                &[(1, photon), (2, photon)],
                &[vec![electron, electron.flip(), photon]],
                3,
            ),
        ),
    ];
    for (name, graphs) in sets {
        let linnet: Vec<HedgeGraph<u8, usize, ()>> =
            graphs.iter().cloned().map(HedgeGraph::from).collect();
        let mut group = c.benchmark_group(name);
        group.bench_function("linnet", |b| {
            b.iter(|| {
                for g in &linnet {
                    black_box(g.canonize());
                }
            })
        });
        group.bench_function("graphica", |b| {
            b.iter(|| {
                for g in &graphs {
                    black_box(g.canonize());
                }
            })
        });
        group.finish();
    }
}

criterion_group!(benches, canon);
criterion_main!(benches);
