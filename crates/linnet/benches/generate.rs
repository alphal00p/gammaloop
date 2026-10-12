//! Native generation against Graphica generate on whole Feynman-graph sets. Both produce one graph
//! per isomorphism class with its symmetry factor; the native side consumes each compact diagram.
//! Graphica returns all loop orders up to the bound, so its set is filtered to the requested order
//! like the native one.
use std::{hint::black_box, time::Duration};

use criterion::{criterion_group, criterion_main, Criterion};
use linnet::half_edge::generate::{Bridges, Legs, Plan, Query, Run, Structure, Stub, Theory};
use symbolica::graph::{GenerationSettings, Graph, HalfEdge};

struct Case {
    name: &'static str,
    /// legs as Graphica sees them: the stub at the external node
    externals: Vec<HalfEdge<u8>>,
    signatures: Vec<Vec<HalfEdge<u8>>>,
    loops: usize,
    one_pi: bool,
    zero_flow_edges: bool,
}

impl Case {
    fn max_vertices(&self) -> usize {
        let min_degree = self.signatures.iter().map(Vec::len).min().unwrap();
        (2 * self.loops + self.externals.len() - 2) / (min_degree - 2)
    }

    fn stub(h: &HalfEdge<u8>) -> Stub<u8> {
        match h.direction {
            None => Stub::undirected(h.data),
            Some(true) => Stub::outgoing(h.data),
            Some(false) => Stub::incoming(h.data),
        }
    }

    fn native(&self) -> usize {
        let legs = self
            .externals
            .iter()
            .enumerate()
            .fold(Legs::new(), |legs, (i, h)| {
                legs.leg(i + 1, Self::stub(&h.flip()))
            });
        let theory = self
            .signatures
            .iter()
            .fold(Theory::builder(), |theory, signature| {
                theory.vertex(signature.iter().map(Self::stub), [])
            })
            .build()
            .unwrap();
        let query = Query::new(&legs)
            .loops(self.loops..=self.loops)
            .structure(Structure {
                self_loops: true,
                zero_flow_edges: self.zero_flow_edges,
                bridges: if self.one_pi {
                    Bridges::OnePI
                } else {
                    Bridges::Any
                },
                ..Structure::default()
            });
        let plan = Plan::new(&theory, &query).unwrap();
        let run = Run::default();
        plan.diagrams(&run).unwrap().map(black_box).count()
    }

    fn graphica(&self) -> usize {
        let externals: Vec<_> = self
            .externals
            .iter()
            .enumerate()
            .map(|(i, &h)| (i + 1, h))
            .collect();
        let mut settings = GenerationSettings::new()
            .max_loops(self.loops)
            .max_vertices(self.externals.len() + self.max_vertices())
            .allow_self_loops(true)
            .allow_zero_flow_edges(self.zero_flow_edges);
        if self.one_pi {
            settings = settings.max_bridges(0);
        }
        Graph::<usize, u8>::generate(&externals, &self.signatures, settings)
            .unwrap()
            .into_iter()
            .filter(|(g, _)| g.num_loops() == self.loops)
            .map(black_box)
            .count()
    }
}

fn generate(c: &mut Criterion) {
    let phi = HalfEdge::undirected(0);
    let (photon, electron) = (HalfEdge::undirected(1), HalfEdge::incoming(2));
    let (gluon, quark, ghost) = (
        HalfEdge::undirected(3),
        HalfEdge::incoming(4),
        HalfEdge::incoming(7),
    );
    let cases = [
        Case {
            name: "phi3 vacuum L5",
            externals: vec![],
            signatures: vec![vec![phi; 3]],
            loops: 5,
            one_pi: false,
            zero_flow_edges: true,
        },
        Case {
            name: "phi34 2-point L3",
            externals: vec![phi; 2],
            signatures: vec![vec![phi; 3], vec![phi; 4]],
            loops: 3,
            one_pi: false,
            zero_flow_edges: true,
        },
        Case {
            name: "QED photon self-energy L3",
            externals: vec![photon; 2],
            signatures: vec![vec![electron, electron.flip(), photon]],
            loops: 3,
            one_pi: true,
            zero_flow_edges: false,
        },
        Case {
            name: "gg -> gg L1",
            externals: vec![gluon; 4],
            signatures: vec![
                vec![quark.flip(), quark, gluon],
                vec![gluon; 3],
                vec![gluon; 4],
                vec![ghost.flip(), ghost, gluon],
            ],
            loops: 1,
            one_pi: false,
            zero_flow_edges: false,
        },
    ];
    for case in &cases {
        assert_eq!(case.native(), case.graphica(), "{}", case.name);
        let mut group = c.benchmark_group(case.name);
        group.bench_function("linnet", |b| b.iter(|| case.native()));
        group.bench_function("graphica", |b| b.iter(|| case.graphica()));
        group.finish();
    }
}

/// Native trees from one signature per degree: the partial stars have factorially many
/// automorphisms and most partial graphs that close a cycle cannot complete.
fn trees(c: &mut Criterion) {
    let undirected = Stub::undirected(0u8);
    let theory = (1..=15)
        .fold(Theory::builder(), |theory, degree| {
            theory.vertex(vec![undirected.clone(); degree], [])
        })
        .build()
        .unwrap();
    let query = Query::new(&Legs::<u8, u8>::new())
        .loops(0..=0)
        .structure(Structure {
            max_vertices: Some(14),
            zero_flow_edges: true,
            ..Structure::default()
        });
    let plan = Plan::new(&theory, &query).unwrap();
    let run = Run::default();
    let count = || plan.diagrams(&run).unwrap().map(black_box).count();
    // OEIS A000055: the free trees with 2 to 14 vertices
    assert_eq!(count(), 5446);
    let mut group = c.benchmark_group("trees with degrees 1 to 15 up to 14 vertices");
    group
        .sample_size(10)
        .measurement_time(Duration::from_secs(8));
    group.bench_function("linnet", |b| b.iter(count));
    group.finish();
}

criterion_group!(benches, generate, trees);
criterion_main!(benches);
