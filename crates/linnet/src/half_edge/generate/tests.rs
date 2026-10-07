use std::{
    collections::HashSet,
    sync::{
        atomic::{AtomicUsize, Ordering},
        Arc,
    },
};

use super::{Diagram, GenerationError, GraphGenerator, Stub};
use crate::half_edge::{involution::HedgePair, HedgeGraph};

/// Field codes of the test models.
const PHI: u8 = 0;
const PHOTON: u8 = 1;
const ELECTRON: u8 = 2;
const GLUON: u8 = 3;
const QUARK: u8 = 4;
const GHOST: u8 = 7;

/// Exact sum of `1 / |Aut|`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Ratio(u128, u128);

impl Ratio {
    fn sum<'a, F: 'a, L: 'a>(diagrams: impl IntoIterator<Item = &'a Diagram<F, L>>) -> Ratio {
        fn gcd(a: u128, b: u128) -> u128 {
            if b == 0 {
                a
            } else {
                gcd(b, a % b)
            }
        }
        diagrams.into_iter().fold(Ratio(0, 1), |Ratio(p, q), d| {
            let s = d.automorphism_count.expect("small symmetry factor");
            let (p, q) = (p * s + q, q * s);
            let g = gcd(p, q);
            Ratio(p / g, q / g)
        })
    }
}

fn run<F, L>(generator: GraphGenerator<F, L>) -> Vec<Diagram<F, L>>
where
    F: Clone + Ord + Send + Sync,
    L: Clone + Ord + Send + Sync,
{
    generator.iter().unwrap().collect()
}

/// Scalar theory with the given vertex degrees and labelled legs.
fn phi(degrees: &[usize], labels: &[u8]) -> GraphGenerator<u8, u8> {
    let phi = Stub::undirected(PHI);
    GraphGenerator::new(
        labels.iter().map(|&l| (l, phi.clone())),
        degrees.iter().map(|&d| vec![phi.clone(); d]),
    )
    .allow_self_loops(true)
    .allow_zero_flow_edges(true)
}

fn qed(externals: Vec<(u8, Stub<u8>)>) -> GraphGenerator<u8, u8> {
    let electron = Stub::incoming(ELECTRON);
    GraphGenerator::new(
        externals,
        [vec![
            electron.clone(),
            electron.flipped(),
            Stub::undirected(PHOTON),
        ]],
    )
    .allow_self_loops(true)
}

/// QCD with gluons, one quark flavour and ghosts.
fn qcd(externals: Vec<(u8, Stub<u8>)>) -> GraphGenerator<u8, u8> {
    let (g, q, c) = (
        Stub::undirected(GLUON),
        Stub::incoming(QUARK),
        Stub::incoming(GHOST),
    );
    GraphGenerator::new(
        externals,
        [
            vec![q.clone().flipped(), q, g.clone()],
            vec![g.clone(); 3],
            vec![g.clone(); 4],
            vec![c.clone().flipped(), c, g],
        ],
    )
    .allow_self_loops(true)
}

fn gluons(labels: &[u8]) -> Vec<(u8, Stub<u8>)> {
    labels
        .iter()
        .map(|&l| (l, Stub::undirected(GLUON)))
        .collect()
}

/// The certificate and count stored with each graph are those of canonizing it, and
/// certificates are pairwise distinct.
fn check_canonical(diagrams: &[Diagram<u8, u8>]) {
    let mut seen = HashSet::new();
    for d in diagrams {
        let canon = d.graph.canonize_by(|&s| s, |&f| f, |&h| h);
        assert_eq!(canon.automorphism_count(), d.automorphism_count);
        assert_eq!(*canon.certificate(), d.certificate);
        assert!(seen.insert(d.certificate.clone()), "duplicate graph");
    }
}

fn loop_count<F, L>(g: &HedgeGraph<F, usize, Option<L>>) -> usize {
    let paired = g
        .iter_edges()
        .filter(|(pair, ..)| matches!(pair, HedgePair::Paired { .. }))
        .count();
    paired + 1 - g.n_nodes()
}

#[test]
fn phi3_vacuum_matches_oeis_and_zero_dimensional_sums() {
    // OEIS A005967 and the zero-dimensional phi^3 free energy
    for (loops, count, sum) in [
        (2, 2, Ratio(5, 24)),
        (3, 5, Ratio(5, 16)),
        (4, 17, Ratio(1105, 1152)),
        (5, 71, Ratio(565, 128)),
    ] {
        let diagrams = run(phi(&[3], &[]).loops(loops..=loops));
        assert_eq!(diagrams.len(), count, "phi3 vacuum L{loops}");
        assert_eq!(Ratio::sum(&diagrams), sum, "phi3 vacuum L{loops}");
        check_canonical(&diagrams);
    }
}

#[test]
fn phi4_vacuum_and_phi34_two_point_match_zero_dimensional_sums() {
    for (loops, count, sum) in [
        (2, 1, Ratio(1, 8)),
        (3, 2, Ratio(1, 12)),
        (4, 4, Ratio(11, 96)),
        (5, 10, Ratio(17, 72)),
    ] {
        let diagrams = run(phi(&[4], &[]).loops(loops..=loops));
        assert_eq!(diagrams.len(), count, "phi4 vacuum L{loops}");
        assert_eq!(Ratio::sum(&diagrams), sum, "phi4 vacuum L{loops}");
    }
    for (loops, count, sum) in [
        (1, 3, Ratio(3, 2)),
        (2, 29, Ratio(25, 3)),
        (3, 351, Ratio(1741, 24)),
    ] {
        let diagrams = run(phi(&[3, 4], &[1, 2]).loops(loops..=loops));
        assert_eq!(diagrams.len(), count, "phi34 2-point L{loops}");
        assert_eq!(Ratio::sum(&diagrams), sum, "phi34 2-point L{loops}");
        check_canonical(&diagrams);
    }
}

/// Counts and sums verified against Graphica (and the QGRAF-validated FeynGraph tests).
#[test]
fn typed_processes_match_graphica_counts() {
    let photons = || vec![(1, Stub::undirected(PHOTON)), (2, Stub::undirected(PHOTON))];
    for (loops, count) in [(1, 1), (2, 3), (3, 20)] {
        let diagrams = run(qed(photons()).loops(loops..=loops).max_bridges(0));
        assert_eq!(diagrams.len(), count, "QED photon self-energy L{loops}");
        assert!(diagrams.iter().all(|d| d.automorphism_count == Some(1)));
        check_canonical(&diagrams);
    }
    for (loops, count, sum) in [(0, 4, Ratio(4, 1)), (1, 148, Ratio(237, 2))] {
        let diagrams = run(qcd(gluons(&[1, 2, 3, 4])).loops(loops..=loops));
        assert_eq!(diagrams.len(), count, "gg -> gg L{loops}");
        assert_eq!(Ratio::sum(&diagrams), sum, "gg -> gg L{loops}");
        check_canonical(&diagrams);
    }
    // three quark flavours
    let (g, flavours) = (Stub::undirected(GLUON), [4, 5, 6]);
    let mut signatures = vec![vec![g.clone(); 3], vec![g.clone(); 4]];
    signatures.extend(flavours.map(|f| {
        let q = Stub::incoming(f);
        vec![q.clone().flipped(), q, g.clone()]
    }));
    for (loops, count, sum) in [(1, 5, Ratio(4, 1)), (2, 30, Ratio(275, 12))] {
        let generator = GraphGenerator::new(gluons(&[1, 2]), signatures.clone())
            .allow_self_loops(true)
            .loops(loops..=loops)
            .max_bridges(0);
        let diagrams = run(generator);
        assert_eq!(diagrams.len(), count, "QCD gluon self-energy L{loops}");
        assert_eq!(Ratio::sum(&diagrams), sum, "QCD gluon self-energy L{loops}");
    }
}

#[test]
fn directed_external_legs_keep_their_direction() {
    let electron = Stub::incoming(ELECTRON);
    // e- e+ -> photon photon: one fermion line enters the diagram and one leaves it
    let externals = vec![
        (1, electron.clone()),
        (2, electron.clone().flipped()),
        (3, Stub::undirected(PHOTON)),
        (4, Stub::undirected(PHOTON)),
    ];
    let diagrams = run(qed(externals.clone()).loops(1..=1));
    assert_eq!(diagrams.len(), 18);
    check_canonical(&diagrams);
    // electron self-energy, one-particle irreducible
    let diagrams = run(qed(externals[..2].to_vec()).loops(3..=3).max_bridges(0));
    assert_eq!(diagrams.len(), 20);
    // two incoming electrons cannot attach to QED vertices alone
    let both_in = vec![(1, electron.clone()), (2, electron)];
    assert!(run(qed(both_in).loops(0..=2).max_vertices(6)).is_empty());
}

#[test]
fn legs_with_equal_labels_are_interchangeable() {
    let diagrams = run(phi(&[3, 4], &[1, 1]).loops(3..=3));
    assert_eq!(diagrams.len(), 240);
    assert_eq!(Ratio::sum(&diagrams), Ratio(1741, 48));
    check_canonical(&diagrams);
    let diagrams = run(qcd(gluons(&[1, 1, 1, 1])).loops(1..=1));
    assert_eq!(diagrams.len(), 22);
    assert_eq!(Ratio::sum(&diagrams), Ratio(79, 16));
    let diagrams = run(qcd(gluons(&[1, 1, 2, 2])).loops(1..=1));
    assert_eq!(diagrams.len(), 61);
    check_canonical(&diagrams);
}

#[test]
fn self_loop_bridge_and_vertex_limits() {
    let diagrams = run(phi(&[3], &[1, 2]).allow_self_loops(false).loops(4..=4));
    assert_eq!(
        (diagrams.len(), Ratio::sum(&diagrams)),
        (143, Ratio(1307, 32))
    );
    let diagrams = run(phi(&[3, 4], &[1, 2]).max_bridges(1).loops(3..=3));
    assert_eq!(
        (diagrams.len(), Ratio::sum(&diagrams)),
        (183, Ratio(201, 4))
    );
    let diagrams = run(phi(&[3, 4], &[1, 2])
        .max_bridges(1)
        .allow_zero_flow_edges(false)
        .loops(3..=3));
    assert_eq!(
        (diagrams.len(), Ratio::sum(&diagrams)),
        (116, Ratio(103, 3))
    );
    // the four-gluon contact graph has one vertex
    assert_eq!(
        run(qcd(gluons(&[1, 2, 3, 4])).loops(0..=0).max_vertices(1)).len(),
        1
    );
    assert!(run(qcd(gluons(&[1, 2, 3, 4])).loops(0..=0).max_vertices(0)).is_empty());
    let diagrams = run(phi(&[3, 4], &[1, 2]).max_vertices(4).loops(3..=3));
    assert_eq!(
        (diagrams.len(), Ratio::sum(&diagrams)),
        (115, Ratio(339, 16))
    );
    assert!(diagrams.iter().all(|d| d.graph.n_nodes() <= 4));
}

#[test]
fn loop_ranges_and_low_degree_vertices() {
    let diagrams = run(phi(&[3], &[]).loops(..=5));
    assert_eq!(diagrams.len(), 2 + 5 + 17 + 71);
    assert!(diagrams.windows(2).all(|w| w[0].loops <= w[1].loops));
    assert!(diagrams.iter().all(|d| loop_count(&d.graph) == d.loops));
    // a vertex bound alone bounds the loop order: eight cubic vertices reach five loops
    let bounded = run(phi(&[3], &[]).max_vertices(8));
    assert_eq!(bounded.len(), diagrams.len());
    // mass insertions (degree two) need a vertex bound
    assert_eq!(
        phi(&[3, 2], &[1, 2]).loops(2..=2).iter().err(),
        Some(GenerationError::Unbounded)
    );
    let diagrams = run(phi(&[3, 2], &[1, 2]).loops(2..=2).max_vertices(5));
    assert_eq!((diagrams.len(), Ratio::sum(&diagrams)), (68, Ratio(25, 1)));
    let diagrams = run(phi(&[3, 2], &[]).loops(2..=2).max_vertices(5));
    assert_eq!((diagrams.len(), Ratio::sum(&diagrams)), (20, Ratio(25, 6)));
    check_canonical(&diagrams);
}

#[test]
fn trees_of_every_degree() {
    // one signature per degree: the free trees, OEIS A000055 by number of vertices
    let degrees: Vec<usize> = (1..=15).collect();
    let diagrams = run(phi(&degrees, &[]).loops(0..=0).max_vertices(11));
    let mut counts = [0; 12];
    for d in &diagrams {
        counts[d.graph.n_nodes()] += 1;
    }
    assert_eq!(counts, [0, 0, 1, 1, 2, 3, 6, 11, 23, 47, 106, 235]);
    check_canonical(&diagrams);
    // the partial stars on the way to the star with 13 leaves have up to 12! automorphisms, which
    // the search must not list
    let diagrams = run(phi(&[1, 13], &[]).loops(0..=0).max_vertices(14));
    let mut symmetries: Vec<_> = diagrams.iter().map(|d| d.automorphism_count).collect();
    symmetries.sort();
    assert_eq!(symmetries, [Some(2), Some((1..=13).product())]);
}

#[test]
fn identical_signatures_are_distinct_vertex_rules() {
    // two rules with the stubs of phi^3: every vertex colouring up to symmetry, so the sum of
    // 1/|Aut| is 2^V times the sum over phi^3 graphs
    let phi = Stub::undirected(PHI);
    let generator = GraphGenerator::<u8, u8>::new([], [vec![phi.clone(); 3], vec![phi; 3]])
        .allow_self_loops(true)
        .allow_zero_flow_edges(true)
        .loops(3..=3);
    let diagrams = run(generator);
    assert_eq!(Ratio::sum(&diagrams), Ratio(5, 1));
    check_canonical(&diagrams);
    let rules: HashSet<usize> = diagrams
        .iter()
        .flat_map(|d| d.graph.iter_nodes().map(|(_, _, &r)| r))
        .collect();
    assert_eq!(rules, HashSet::from([0, 1]));
}

#[test]
fn filter_keeps_exactly_the_accepted_graphs() {
    let no_loops = |g: &HedgeGraph<u8, usize, Option<u8>>| {
        g.iter_edges().all(|(pair, ..)| match pair {
            HedgePair::Paired { source, sink } => g.node_id(source) != g.node_id(sink),
            _ => true,
        })
    };
    let filtered = run(phi(&[3], &[1, 2]).loops(4..=4).filter(no_loops));
    let disallowed = run(phi(&[3], &[1, 2]).allow_self_loops(false).loops(4..=4));
    let certificates = |d: &[Diagram<u8, u8>]| {
        d.iter()
            .map(|d| d.certificate.clone())
            .collect::<HashSet<_>>()
    };
    assert_eq!(certificates(&filtered), certificates(&disallowed));
}

#[test]
fn abort_check_and_progress_stop_the_generation() {
    let polls = Arc::new(AtomicUsize::new(0));
    let counter = polls.clone();
    let generator = phi(&[3, 4], &[1, 2])
        .loops(3..=3)
        .abort_check(move || counter.fetch_add(1, Ordering::Relaxed) >= 20);
    let mut diagrams = generator.iter().unwrap();
    let count = diagrams.by_ref().count();
    assert!(diagrams.aborted());
    assert!(count < 351);

    let generator = phi(&[3, 4], &[1, 2])
        .loops(3..=3)
        .progress(|progress| progress.diagrams < 3);
    let mut diagrams = generator.iter().unwrap();
    assert_eq!(diagrams.by_ref().count(), 3);
    assert!(diagrams.aborted());
}

#[test]
fn generation_streams() {
    // each graph is searched for when it is requested: phi^3 + phi^4 has one graph per
    // topology, so the first three of the 84,749 two-point graphs at five loops cost three
    // topologies
    let topologies = Arc::new(AtomicUsize::new(0));
    let seen = topologies.clone();
    let generator = phi(&[3, 4], &[1, 2])
        .loops(5..=5)
        .progress(move |progress| {
            seen.store(progress.topologies, Ordering::Relaxed);
            true
        });
    let first: Vec<_> = generator.iter().unwrap().take(3).collect();
    assert_eq!(first.len(), 3);
    assert_eq!(topologies.load(Ordering::Relaxed), 3);
    check_canonical(&first);

    // progress sees every topology, so it stops a run whose graphs are all filtered out
    let seen = topologies.clone();
    let generator = phi(&[3, 4], &[1, 2])
        .loops(3..=3)
        .filter(|_| false)
        .progress(move |progress| {
            seen.store(progress.topologies, Ordering::Relaxed);
            progress.topologies < 5
        });
    let mut diagrams = generator.iter().unwrap();
    assert!(diagrams.next().is_none());
    assert!(diagrams.aborted());
    assert_eq!(topologies.load(Ordering::Relaxed), 5);
}

#[test]
fn without_signatures_nothing_is_generated() {
    let legs = [(1, Stub::undirected(PHI)), (2, Stub::undirected(PHI))];
    let none = Vec::<Vec<Stub<u8>>>::new();
    assert!(run(GraphGenerator::new(legs.clone(), none.clone()).loops(0..=2)).is_empty());
    assert!(run(GraphGenerator::<u8, u8>::new([], none.clone()).loops(0..=2)).is_empty());
    assert!(run(GraphGenerator::new(legs, none)).is_empty());
}

#[test]
fn one_leg_graphs_without_bridges_are_kept() {
    // the leg is not a bridge: without zero-flow edges, or without bridges, phi^3 keeps the
    // tadpole (leg and self-loop on one vertex) at one loop and the triangle with a doubled side
    // at two loops; Graphica returns neither
    for generator in [
        phi(&[3], &[1]).allow_zero_flow_edges(false),
        phi(&[3], &[1]).max_bridges(0),
    ] {
        let diagrams = run(generator.loops(1..=2));
        let summary: Vec<_> = diagrams
            .iter()
            .map(|d| (d.loops, d.graph.n_nodes(), d.automorphism_count))
            .collect();
        assert_eq!(summary, [(1, 1, Some(2)), (2, 3, Some(4))]);
    }
    // zero-flow bridges add the two graphs with tadpoles on bridges
    assert_eq!(run(phi(&[3], &[1]).loops(2..=2)).len(), 3);
}

#[test]
fn output_order_does_not_depend_on_threads() {
    let generators = || {
        [
            phi(&[3], &[]).loops(5..=5),
            phi(&[3, 4], &[1, 2]).loops(3..=3),
            qcd(gluons(&[1, 2, 3, 4])).loops(1..=1),
        ]
    };
    for (sequential, parallel) in generators().into_iter().zip(generators()) {
        let order =
            |d: Vec<Diagram<u8, u8>>| d.into_iter().map(|d| d.certificate).collect::<Vec<_>>();
        let expected = order(run(sequential));
        let parallel = parallel.parallel(true);
        for threads in [1, 2, 4] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            assert_eq!(pool.install(|| order(run_ref(&parallel))), expected);
        }
    }
}

fn run_ref<F, L>(generator: &GraphGenerator<F, L>) -> Vec<Diagram<F, L>>
where
    F: Clone + Ord + Send + Sync,
    L: Clone + Ord + Send + Sync,
{
    generator.iter().unwrap().collect()
}

#[test]
fn invalid_settings_are_errors() {
    let phi = Stub::undirected(PHI);
    let unbounded = GraphGenerator::<u8, u8>::new([], [vec![phi.clone(); 3]]);
    assert_eq!(unbounded.iter().err(), Some(GenerationError::Unbounded));
    let wide = GraphGenerator::<u8, u8>::new([], [vec![phi.clone(); 16]]).loops(..=2);
    assert_eq!(
        wide.iter().err(),
        Some(GenerationError::UnsupportedDegree(16))
    );
    let large = GraphGenerator::<u8, u8>::new([], [vec![phi; 3]]).loops(40..=40);
    assert_eq!(large.iter().err(), Some(GenerationError::TooManyVertices));
}

#[cfg(feature = "symbolica")]
mod graphica {
    use std::collections::HashMap;

    use symbolica::graph::{GenerationSettings, Graph, HalfEdge};

    use super::*;
    use crate::half_edge::{
        builder::HedgeGraphBuilder,
        canon::Certificate,
        involution::{Flow, Orientation},
    };

    type Plain = Graph<u8, u8>;
    type Key = Certificate<usize, u8, Option<u8>>;

    /// A case for both generators: externals as Graphica sees them (stub at the external node,
    /// node data the label) and the vertex signatures.
    struct Case {
        externals: Vec<(u8, HalfEdge<u8>)>,
        signatures: Vec<Vec<HalfEdge<u8>>>,
        loops: usize,
        max_bridges: Option<usize>,
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

        fn native(&self) -> HashMap<Key, u128> {
            // Graphica's external half-edge is the stub at the external node
            let externals = self
                .externals
                .iter()
                .map(|(l, h)| (*l, Self::stub(&h.flip())));
            let signatures = self
                .signatures
                .iter()
                .map(|s| s.iter().map(Self::stub).collect());
            let mut generator = GraphGenerator::new(externals, signatures)
                .loops(self.loops..=self.loops)
                .allow_self_loops(true)
                .allow_zero_flow_edges(self.zero_flow_edges);
            if let Some(b) = self.max_bridges {
                generator = generator.max_bridges(b);
            }
            generator
                .iter()
                .unwrap()
                .map(|d| (d.certificate, d.automorphism_count.unwrap()))
                .collect()
        }

        fn graphica(&self) -> HashMap<Key, u128> {
            let mut settings = GenerationSettings::new()
                .max_loops(self.loops)
                .max_vertices(self.externals.len() + self.max_vertices())
                .allow_self_loops(true)
                .allow_zero_flow_edges(self.zero_flow_edges);
            if let Some(b) = self.max_bridges {
                settings = settings.max_bridges(b);
            }
            Plain::generate(&self.externals, &self.signatures, settings)
                .unwrap()
                .into_iter()
                .filter(|(g, _)| g.num_loops() == self.loops)
                .map(|(g, aut)| {
                    let (certificate, _) =
                        self.to_linnet(&g).certificate_by(|&s| s, |&f| f, |&h| h);
                    (certificate, aut.to_string().parse().unwrap())
                })
                .collect()
        }

        /// Graphica's output in the generator's convention: external nodes become dangling
        /// hedges keyed by their label, vertices carry the index of their signature.
        fn to_linnet(&self, g: &Plain) -> HedgeGraph<u8, usize, Option<u8>> {
            let mut stubs = vec![vec![]; g.nodes().len()];
            for e in g.edges() {
                let (a, b) = e.vertices;
                let (at_a, at_b) = if e.directed {
                    (HalfEdge::outgoing(e.data), HalfEdge::incoming(e.data))
                } else {
                    (HalfEdge::undirected(e.data), HalfEdge::undirected(e.data))
                };
                stubs[a].push(at_a);
                stubs[b].push(at_b);
            }
            let mut builder = HedgeGraphBuilder::new();
            let nodes: Vec<_> = (0..g.nodes().len())
                .map(|v| {
                    (g.nodes()[v].data == 0).then(|| {
                        stubs[v].sort();
                        let rule = self
                            .signatures
                            .iter()
                            .position(|s| {
                                let mut s = s.clone();
                                s.sort();
                                s == stubs[v]
                            })
                            .expect("a vertex matches a signature");
                        builder.add_node(rule)
                    })
                })
                .collect();
            for e in g.edges() {
                let (x, y) = e.vertices;
                let o = if e.directed {
                    Orientation::Default
                } else {
                    Orientation::Undirected
                };
                match (nodes[x], nodes[y]) {
                    (Some(a), Some(b)) => {
                        builder.add_edge(a.add_data(None), b.add_data(None), e.data, o)
                    }
                    (Some(a), None) => builder.add_external_edge(
                        a.add_data(Some(g.nodes()[y].data)),
                        e.data,
                        o,
                        Flow::Source,
                    ),
                    (None, Some(b)) => builder.add_external_edge(
                        b.add_data(Some(g.nodes()[x].data)),
                        e.data,
                        o,
                        Flow::Sink,
                    ),
                    (None, None) => unreachable!("legs join vertices"),
                }
            }
            builder.build()
        }
    }

    #[test]
    fn sets_and_symmetry_factors_match_graphica() {
        let phi = HalfEdge::undirected(PHI);
        let photon = HalfEdge::undirected(PHOTON);
        let gluon = HalfEdge::undirected(GLUON);
        let electron = HalfEdge::incoming(ELECTRON);
        let qed = vec![vec![electron, electron.flip(), photon]];
        let qcd = vec![
            vec![HalfEdge::outgoing(QUARK), HalfEdge::incoming(QUARK), gluon],
            vec![gluon; 3],
            vec![gluon; 4],
            vec![HalfEdge::outgoing(GHOST), HalfEdge::incoming(GHOST), gluon],
        ];
        let labelled = |legs: &[HalfEdge<u8>]| {
            legs.iter()
                .enumerate()
                .map(|(i, &h)| (i as u8 + 1, h))
                .collect::<Vec<_>>()
        };
        let cases = [
            Case {
                externals: vec![],
                signatures: vec![vec![phi; 3]],
                loops: 4,
                max_bridges: None,
                zero_flow_edges: true,
            },
            Case {
                externals: labelled(&[phi, phi]),
                signatures: vec![vec![phi; 3], vec![phi; 4]],
                loops: 3,
                max_bridges: None,
                zero_flow_edges: true,
            },
            Case {
                externals: vec![(1, phi), (1, phi)],
                signatures: vec![vec![phi; 3], vec![phi; 4]],
                loops: 2,
                max_bridges: Some(1),
                zero_flow_edges: false,
            },
            Case {
                externals: labelled(&[photon, photon]),
                signatures: qed.clone(),
                loops: 3,
                max_bridges: Some(0),
                zero_flow_edges: false,
            },
            Case {
                externals: labelled(&[electron.flip(), electron, photon, photon]),
                signatures: qed,
                loops: 1,
                max_bridges: None,
                zero_flow_edges: true,
            },
            Case {
                externals: labelled(&[gluon, gluon, gluon, gluon]),
                signatures: qcd.clone(),
                loops: 1,
                max_bridges: None,
                zero_flow_edges: false,
            },
            Case {
                externals: labelled(&[HalfEdge::outgoing(GHOST), HalfEdge::incoming(GHOST)]),
                signatures: qcd,
                loops: 2,
                max_bridges: Some(0),
                zero_flow_edges: false,
            },
        ];
        for (i, case) in cases.iter().enumerate() {
            let (native, graphica) = (case.native(), case.graphica());
            assert!(!native.is_empty(), "case {i} is empty");
            assert_eq!(native, graphica, "case {i}");
        }
    }

    /// The one documented difference: Graphica counts the leg of a one-leg process as a
    /// zero-flow bridge.
    #[test]
    fn one_leg_processes_differ_only_without_zero_flow_edges_or_bridges() {
        let phi = HalfEdge::undirected(PHI);
        let case = |loops, max_bridges, zero_flow_edges| Case {
            externals: vec![(1, phi)],
            signatures: vec![vec![phi; 3]],
            loops,
            max_bridges,
            zero_flow_edges,
        };
        for loops in 1..=3 {
            let all = case(loops, None, true);
            assert_eq!(all.native(), all.graphica(), "L{loops}");
            for restricted in [case(loops, None, false), case(loops, Some(0), true)] {
                assert!(restricted.graphica().is_empty(), "L{loops}");
                let native = restricted.native();
                assert!(!native.is_empty(), "L{loops}");
                assert!(native.keys().all(|k| all.native().contains_key(k)));
            }
        }
    }
}
