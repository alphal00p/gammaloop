use std::{
    collections::{BTreeMap, HashSet},
    sync::{
        atomic::{AtomicUsize, Ordering},
        Arc,
    },
};

use proptest::prelude::*;

use super::{
    Compiled, Countable, Diagram, EdgeView, GenerationError, GraphGenerator, GraphView, Stage,
    Stub, VertexSet,
};
use crate::half_edge::{canon::CanonicalForm, involution::HedgePair, HedgeGraph};

/// Field codes of the test models.
const PHI: u8 = 0;
const PHOTON: u8 = 1;
const ELECTRON: u8 = 2;
const GLUON: u8 = 3;
const QUARK: u8 = 4;
const GHOST: u8 = 7;
/// A second scalar.
const CHI: u8 = 8;

/// Exact sum of `1 / |Aut|`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Ratio(u128, u128);

impl Ratio {
    fn new(p: u128, q: u128) -> Ratio {
        fn gcd(a: u128, b: u128) -> u128 {
            if b == 0 {
                a
            } else {
                gcd(b, a % b)
            }
        }
        let g = gcd(p, q);
        Ratio(p / g, q / g)
    }

    fn sum<'a, F: 'a, L: 'a>(diagrams: impl IntoIterator<Item = &'a Diagram<F, L>>) -> Ratio {
        diagrams.into_iter().fold(Ratio(0, 1), |Ratio(p, q), d| {
            let s = d.automorphism_count.expect("small symmetry factor");
            Ratio::new(p * s + q, q * s)
        })
    }

    fn divided(self, k: u128) -> Ratio {
        Ratio::new(self.0, self.1 * k)
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

/// QCD (signatures 0 to 3: qqg, ggg, gggg, ccg) and, with `qed`, the photon couplings of the quark
/// and the electron (4: qqa, 5: eea).
fn qcd_qed(externals: Vec<(u8, Stub<u8>)>, qed: bool) -> GraphGenerator<u8, u8> {
    let (g, q, c, a, e) = (
        Stub::undirected(GLUON),
        Stub::incoming(QUARK),
        Stub::incoming(GHOST),
        Stub::undirected(PHOTON),
        Stub::incoming(ELECTRON),
    );
    let mut signatures = vec![
        vec![q.clone().flipped(), q.clone(), g.clone()],
        vec![g.clone(); 3],
        vec![g.clone(); 4],
        vec![c.clone().flipped(), c, g],
    ];
    if qed {
        signatures.push(vec![q.clone().flipped(), q, a.clone()]);
        signatures.push(vec![e.clone().flipped(), e, a]);
    }
    GraphGenerator::new(externals, signatures).allow_self_loops(true)
}

/// QCD with gluons, one quark flavour and ghosts.
fn qcd(externals: Vec<(u8, Stub<u8>)>) -> GraphGenerator<u8, u8> {
    qcd_qed(externals, false)
}

fn gluons(labels: &[u8]) -> Vec<(u8, Stub<u8>)> {
    labels
        .iter()
        .map(|&l| (l, Stub::undirected(GLUON)))
        .collect()
}

/// The canonical form and count stored with each graph are those of canonizing it, and
/// canonical forms are pairwise distinct.
fn check_canonical(diagrams: &[Diagram<u8, u8>]) {
    let mut seen = HashSet::new();
    for d in diagrams {
        let canon = d.graph.canonize_by(|&s| s, |&f| f, |&h| h);
        assert_eq!(canon.automorphism_count(), d.automorphism_count);
        assert_eq!(*canon.canonical_form(), d.canonical_form);
        assert!(seen.insert(d.canonical_form.clone()), "duplicate graph");
    }
}

fn loop_count<F, L>(g: &HedgeGraph<F, usize, Option<L>>) -> usize {
    let paired = g
        .iter_edges()
        .filter(|(pair, ..)| matches!(pair, HedgePair::Paired { .. }))
        .count();
    paired + 1 - g.n_live()
}

/// The |Aut| multiset: (automorphism count, number of graphs).
type Histogram = &'static [(u128, usize)];

/// `diagrams` are `count` graphs with the oracle's sum of `1/|Aut|` and |Aut| multiset, and each
/// carries its canonizer's canonical form and automorphism count.
fn assert_oracle(
    case: &str,
    diagrams: &[Diagram<u8, u8>],
    count: usize,
    sum: Ratio,
    histogram: &[(u128, usize)],
) {
    assert_eq!(diagrams.len(), count, "{case}");
    assert_eq!(Ratio::sum(diagrams), sum, "{case}");
    let mut found = BTreeMap::<u128, usize>::new();
    for d in diagrams {
        *found.entry(d.automorphism_count.unwrap()).or_default() += 1;
    }
    assert_eq!(found.into_iter().collect::<Vec<_>>(), histogram, "{case}");
    check_canonical(diagrams);
}

// The oracles below are FeynGraph and Graphica, which agree on every count and |Aut| multiset
// (Graphica alone for the phi^3 vacua at six and seven loops, where FeynGraph did not finish),
// together with OEIS A005967 for the phi^3 vacuum counts and zero-dimensional field theory for
// the scalar sums.

#[test]
fn phi3_vacuum_matches_oeis_and_zero_dimensional_sums() {
    let cases: [(usize, usize, Ratio, Histogram); 6] = [
        (2, 2, Ratio(5, 24), &[(8, 1), (12, 1)]),
        (3, 5, Ratio(5, 16), &[(8, 1), (16, 2), (24, 1), (48, 1)]),
        (
            4,
            17,
            Ratio(1105, 1152),
            &[
                (8, 3),
                (12, 1),
                (16, 5),
                (32, 4),
                (48, 2),
                (72, 1),
                (128, 1),
            ],
        ),
        (
            5,
            71,
            Ratio(565, 128),
            &[
                (4, 4),
                (8, 10),
                (12, 1),
                (16, 21),
                (32, 17),
                (48, 1),
                (64, 11),
                (96, 2),
                (128, 2),
                (256, 2),
            ],
        ),
        (
            6,
            388,
            Ratio(82825, 3072),
            &[
                (2, 6),
                (4, 25),
                (6, 2),
                (8, 60),
                (12, 1),
                (16, 88),
                (20, 2),
                (32, 89),
                (48, 5),
                (64, 67),
                (120, 1),
                (128, 27),
                (256, 6),
                (320, 2),
                (384, 2),
                (512, 4),
                (3072, 1),
            ],
        ),
        (
            7,
            2592,
            Ratio(19675, 96),
            &[
                (1, 5),
                (2, 70),
                (4, 219),
                (8, 408),
                (12, 4),
                (16, 529),
                (18, 1),
                (24, 3),
                (32, 524),
                (36, 1),
                (48, 6),
                (64, 425),
                (96, 6),
                (128, 246),
                (256, 103),
                (384, 3),
                (512, 22),
                (768, 2),
                (1024, 12),
                (2048, 2),
                (3072, 1),
            ],
        ),
    ];
    for (loops, count, sum, histogram) in cases {
        let diagrams = run(phi(&[3], &[]).loops(loops..=loops));
        assert_oracle(
            &format!("phi3 vacuum L{loops}"),
            &diagrams,
            count,
            sum,
            histogram,
        );
    }
}

#[test]
fn phi4_vacuum_and_phi34_two_point_match_zero_dimensional_sums() {
    let vacuum: [(usize, usize, Ratio, Histogram); 6] = [
        (2, 1, Ratio(1, 8), &[(8, 1)]),
        (3, 2, Ratio(1, 12), &[(16, 1), (48, 1)]),
        (4, 4, Ratio(11, 96), &[(24, 1), (32, 1), (48, 2)]),
        (
            5,
            10,
            Ratio(17, 72),
            &[(16, 1), (32, 3), (48, 2), (64, 1), (128, 2), (144, 1)],
        ),
        (
            6,
            28,
            Ratio(619, 960),
            &[
                (16, 4),
                (32, 6),
                (48, 3),
                (64, 4),
                (96, 4),
                (120, 1),
                (128, 2),
                (144, 1),
                (288, 1),
                (320, 2),
            ],
        ),
        (
            7,
            97,
            Ratio(709, 324),
            &[
                (8, 4),
                (16, 7),
                (24, 2),
                (32, 13),
                (48, 6),
                (64, 23),
                (96, 13),
                (128, 9),
                (144, 1),
                (192, 6),
                (256, 3),
                (288, 3),
                (384, 1),
                (576, 2),
                (768, 3),
                (1296, 1),
            ],
        ),
    ];
    for (loops, count, sum, histogram) in vacuum {
        let diagrams = run(phi(&[4], &[]).loops(loops..=loops));
        assert_oracle(
            &format!("phi4 vacuum L{loops}"),
            &diagrams,
            count,
            sum,
            histogram,
        );
    }
    let two_point: [(usize, usize, Ratio, Histogram); 5] = [
        (1, 3, Ratio(3, 2), &[(2, 3)]),
        (2, 29, Ratio(25, 3), &[(2, 6), (4, 19), (6, 2), (8, 2)]),
        (
            3,
            351,
            Ratio(1741, 24),
            &[
                (1, 2),
                (2, 36),
                (4, 125),
                (6, 3),
                (8, 141),
                (12, 22),
                (16, 21),
                (48, 1),
            ],
        ),
        (
            4,
            5076,
            Ratio(80299, 96),
            &[
                (1, 40),
                (2, 370),
                (4, 1180),
                (6, 13),
                (8, 1719),
                (12, 130),
                (16, 1157),
                (24, 216),
                (32, 211),
                (36, 5),
                (48, 13),
                (64, 5),
                (72, 2),
                (96, 13),
                (128, 2),
            ],
        ),
        (
            5,
            84749,
            Ratio(6869123, 576),
            &[
                (1, 727),
                (2, 5060),
                (4, 15053),
                (6, 145),
                (8, 23639),
                (12, 1092),
                (16, 21069),
                (24, 2619),
                (32, 10408),
                (36, 27),
                (48, 2146),
                (64, 2057),
                (72, 115),
                (96, 296),
                (128, 108),
                (144, 22),
                (192, 130),
                (256, 22),
                (288, 7),
                (384, 7),
            ],
        ),
    ];
    for (loops, count, sum, histogram) in two_point {
        let diagrams = run(phi(&[3, 4], &[1, 2]).loops(loops..=loops));
        assert_oracle(
            &format!("phi34 2-point L{loops}"),
            &diagrams,
            count,
            sum,
            histogram,
        );
    }
}

/// QED with one charged fermion, QCD with one quark flavour and ghosts (gg -> gg without
/// zero-flow edges), and QCD with three quark flavours and no ghosts (the gluon self-energy at
/// three loops is the QGRAF-validated count of FeynGraph's tests).
#[test]
fn typed_processes_match_graphica_counts() {
    let photons = || vec![(1, Stub::undirected(PHOTON)), (2, Stub::undirected(PHOTON))];
    for (loops, count) in [(1, 1), (2, 3), (3, 20), (4, 189)] {
        let diagrams = run(qed(photons()).loops(loops..=loops).max_bridges(0));
        let case = format!("QED photon self-energy L{loops}");
        assert_oracle(
            &case,
            &diagrams,
            count,
            Ratio(count as u128, 1),
            &[(1, count)],
        );
    }
    let gg: [(usize, usize, Ratio, Histogram); 3] = [
        (0, 4, Ratio(4, 1), &[(1, 4)]),
        (1, 148, Ratio(237, 2), &[(1, 89), (2, 59)]),
        (
            2,
            4276,
            Ratio(36209, 12),
            &[(1, 2036), (2, 1692), (4, 529), (6, 19)],
        ),
    ];
    for (loops, count, sum, histogram) in gg {
        let diagrams = run(qcd(gluons(&[1, 2, 3, 4])).loops(loops..=loops));
        assert_oracle(
            &format!("gg -> gg L{loops}"),
            &diagrams,
            count,
            sum,
            histogram,
        );
    }
    // three quark flavours
    let (g, flavours) = (Stub::undirected(GLUON), [4, 5, 6]);
    let mut signatures = vec![vec![g.clone(); 3], vec![g.clone(); 4]];
    signatures.extend(flavours.map(|f| {
        let q = Stub::incoming(f);
        vec![q.clone().flipped(), q, g.clone()]
    }));
    let self_energy: [(usize, usize, Ratio, Histogram); 3] = [
        (1, 5, Ratio(4, 1), &[(1, 3), (2, 2)]),
        (2, 30, Ratio(275, 12), &[(1, 18), (2, 8), (4, 3), (6, 1)]),
        (
            3,
            479,
            Ratio(1445, 4),
            &[(1, 272), (2, 155), (4, 42), (6, 1), (8, 8), (12, 1)],
        ),
    ];
    for (loops, count, sum, histogram) in self_energy {
        let generator = GraphGenerator::new(gluons(&[1, 2]), signatures.clone())
            .allow_self_loops(true)
            .loops(loops..=loops)
            .max_bridges(0);
        let case = format!("QCD gluon self-energy L{loops}");
        assert_oracle(&case, &run(generator), count, sum, histogram);
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

/// Distinct legs with a common stub are labelled on shaped topologies: phi^3 has `(2n - 5)!!`
/// trees with `n` distinct legs, phi^3 + phi^4 the counts verified against Graphica, and every
/// tree is rigid.
#[test]
fn distinct_legs_are_labelled_on_shaped_topologies() {
    for (n, phi3, phi34) in [(4, 3, 4), (5, 15, 25), (6, 105, 220), (7, 945, 2485)] {
        let labels: Vec<u8> = (1..=n).collect();
        for (degrees, count) in [(&[3][..], phi3), (&[3, 4][..], phi34)] {
            let diagrams = run(phi(degrees, &labels).loops(0..=0));
            assert_eq!(diagrams.len(), count, "{degrees:?} trees with {n} legs");
            assert!(diagrams.iter().all(|d| d.automorphism_count == Some(1)));
            if n <= 6 {
                check_canonical(&diagrams);
            }
        }
    }
    // e- e+ -> n photons: n! trees, one per order of the photons along the electron line
    let electron = Stub::incoming(ELECTRON);
    for n in 1..=5u8 {
        let mut externals = vec![(1, electron.clone()), (2, electron.clone().flipped())];
        externals.extend((3..3 + n).map(|l| (l, Stub::undirected(PHOTON))));
        let diagrams = run(qed(externals).loops(0..=0));
        assert_eq!(diagrams.len(), (1..=n as usize).product::<usize>());
        check_canonical(&diagrams);
    }
    // the pendants completing a tree hang on bridges: five legs make two internal bridges
    let five = || phi(&[3], &[1, 2, 3, 4, 5]).loops(0..=0);
    assert!(run(five().max_bridges(1)).is_empty());
    assert_eq!(run(five().max_bridges(2)).len(), 15);
}

/// Equal labels divide the sum of `1/|Aut|` by the factorials of their multiplicities (orbit
/// counting), whether the other legs are labelled on shaped topologies or not.
#[test]
fn equal_labels_divide_the_sums_of_distinct_ones() {
    let sum = |generator: GraphGenerator<u8, u8>| Ratio::sum(&run(generator));
    for loops in 0..=1 {
        let gg = |labels: &[u8]| sum(qcd(gluons(labels)).loops(loops..=loops));
        let distinct = gg(&[1, 2, 3, 4, 5]);
        assert_eq!(gg(&[1, 1, 2, 3, 4]), distinct.divided(2), "L{loops}");
        assert_eq!(gg(&[1, 1, 1, 2, 3]), distinct.divided(6), "L{loops}");
    }
    for (loops, n) in [(1, 5u8), (2, 4)] {
        let phi3 = |labels: &[u8]| {
            let generator = phi(&[3], labels).allow_zero_flow_edges(false);
            sum(generator.loops(loops..=loops))
        };
        let distinct: Vec<u8> = (1..=n).collect();
        let factorial = (1..=n as u128).product();
        assert_eq!(
            phi3(&vec![1; n as usize]),
            phi3(&distinct).divided(factorial)
        );
    }
}

#[test]
fn many_distinct_legs_stream() {
    // forty distinct legs: phi^3 has one graph per labelled topology, and the first three trees
    // cost three of them
    let topologies = Arc::new(AtomicUsize::new(0));
    let seen = topologies.clone();
    let labels: Vec<u8> = (1..=40).collect();
    let generator = phi(&[3], &labels).loops(0..=0).progress(move |progress| {
        seen.store(progress.topologies, Ordering::Relaxed);
        true
    });
    let first: Vec<_> = generator.iter().unwrap().take(3).collect();
    assert_eq!(topologies.load(Ordering::Relaxed), 3);
    assert!(first
        .iter()
        .all(|d| d.graph.n_live() == 38 && d.automorphism_count == Some(1)));
    check_canonical(&first);
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
    assert!(diagrams.iter().all(|d| d.graph.n_live() <= 4));
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
        counts[d.graph.n_live()] += 1;
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
fn models_with_more_labels_than_a_mask_holds() {
    // 65 charged fermions and the photon give 131 edge labels, more than the label masks hold, so
    // the dressing scans the labels; up to two loops every one-particle-irreducible photon
    // self-energy has one fermion loop, so each graph of one flavour comes in 65 flavours
    let photon = Stub::undirected(0u8);
    let self_energies = |flavours: u8| {
        let signatures = (1..=flavours)
            .map(|f| vec![Stub::incoming(f), Stub::outgoing(f), photon.clone()])
            .collect::<Vec<_>>();
        let legs = [(1u8, photon.clone()), (2, photon.clone())];
        run(GraphGenerator::new(legs, signatures)
            .loops(1..=2)
            .max_bridges(0))
    };
    let one = self_energies(1);
    let many = self_energies(65);
    assert_eq!(one.len(), 4);
    assert_eq!(many.len(), 65 * one.len());
    let Ratio(p, q) = Ratio::sum(&one);
    assert_eq!(Ratio::sum(&many), Ratio::new(65 * p, q));
    check_canonical(&many);
}

/// The graphs of `generator` and the number of partial field assignments its search shows the
/// filter, less the number of edges of the graphs: zero when it shows each graph's prefixes once
/// and nothing else.
fn dead_field_assignments(generator: GraphGenerator<u8, u8>) -> (Vec<Diagram<u8, u8>>, isize) {
    let fields = Arc::new(AtomicUsize::new(0));
    let seen = fields.clone();
    let diagrams = run(generator.filter(move |view| {
        if view.stage() == Stage::Fields {
            seen.fetch_add(1, Ordering::Relaxed);
        }
        true
    }));
    let edges: usize = diagrams
        .iter()
        .map(|d| {
            d.graph
                .iter_edges()
                .filter(|(pair, ..)| matches!(pair, HedgePair::Paired { .. }))
                .count()
        })
        .sum();
    let dead = fields.load(Ordering::Relaxed) as isize - edges as isize;
    (diagrams, dead)
}

#[test]
fn forward_checking_leaves_no_dead_field_assignments_on_gluon_trees() {
    // a quark or ghost line in a tree would have to end at a leg, so with gluon legs every
    // topology has only its all-gluon labelling; forward checking drops each quark and ghost
    // before the filter sees it, so the filter sees just the prefixes of the emitted labellings
    let (diagrams, dead) = dead_field_assignments(qcd(gluons(&[1, 2, 3, 4, 5, 6])).loops(0..=0));
    assert_eq!(diagrams.len(), 220);
    assert_eq!(dead, 0);
}

#[test]
fn topologies_with_an_edge_no_field_fits_get_no_field_assignment() {
    // phi^3 and an unrelated chi^4: no field fits an edge from a cubic to a quartic vertex, so the
    // vacuum topologies with both degrees are dropped before their first field, even where the
    // first bundle is far from such an edge; the others have one labelling each
    let signatures = [
        vec![Stub::undirected(PHI); 3],
        vec![Stub::undirected(CHI); 4],
    ];
    let both = GraphGenerator::<u8, u8>::new([], signatures)
        .allow_self_loops(true)
        .allow_zero_flow_edges(true)
        .loops(2..=4);
    let (diagrams, dead) = dead_field_assignments(both);
    let pure = |degree| run(phi(&[degree], &[]).loops(2..=4)).len();
    assert_eq!(diagrams.len(), pure(3) + pure(4));
    assert_eq!(dead, 0);
}

#[test]
fn filter_keeps_exactly_the_accepted_graphs() {
    let no_loops = |g: &HedgeGraph<u8, usize, Option<u8>>| {
        g.iter_edges().all(|(pair, ..)| match pair {
            HedgePair::Paired { source, sink } => g.node_id(source) != g.node_id(sink),
            _ => true,
        })
    };
    // on complete graphs only, and pruning partial topologies with self-loops
    let complete = run(phi(&[3], &[1, 2])
        .loops(4..=4)
        .filter(move |view| view.complete_graph().is_none_or(no_loops)));
    let pruned = run(phi(&[3], &[1, 2])
        .loops(4..=4)
        .filter(|view| (0..view.n_vertices()).all(|v| view.self_loops(v) == 0)));
    let disallowed = run(phi(&[3], &[1, 2]).allow_self_loops(false).loops(4..=4));
    let set = |d: &[Diagram<u8, u8>]| canonical_forms(d).into_iter().collect::<HashSet<_>>();
    assert_eq!(set(&complete), set(&disallowed));
    assert_eq!(set(&pruned), set(&disallowed));
}

fn canonical_forms(diagrams: &[Diagram<u8, u8>]) -> Vec<CanonicalForm<usize, u8, Option<u8>>> {
    diagrams.iter().map(|d| d.canonical_form.clone()).collect()
}

type View<'a> = GraphView<'a, u8, u8>;

/// No self-energy insertion on an external leg: no bridge cuts off a side with a single leg and a
/// loop. Decided in stage 1 as soon as such a side is final, and by the complete topology.
fn no_leg_self_energies(view: &View<'_>) -> bool {
    view.stage() != Stage::Topology
        || view.bridges().iter().all(|bridge| {
            let legs: usize = bridge.side.iter().map(|v| view.legs(v).count()).sum();
            legs != 1 || view.loops_in(bridge.side) == 0
        })
}

/// No internal edge carries `field` (a particle veto): decided edge by edge in stage 2.
fn without_field(field: u8) -> impl Fn(&View<'_>) -> bool + Clone + Send + Sync {
    move |view| {
        view.edges()
            .all(|e| e.stub.is_none_or(|stub| *stub.field != field))
    }
}

/// No vertex uses signature `rule` (a vertex veto): decided once every rule that fits a vertex is
/// vetoed, in stage 1 for the four-gluon vertex (the only one of degree four).
fn without_rule(rule: u32) -> impl Fn(&View<'_>) -> bool + Clone + Send + Sync {
    move |view| {
        (0..view.n_vertices()).all(|v| match view.rules(v) {
            Some(rules) => rules.iter().any(|&r| r != rule),
            None => rule != 2 || view.degree(v) != 4,
        })
    }
}

/// QCD coupling order at most `max` in [`qcd_qed`]: the sum over vertices of the smallest order
/// of the signatures that can still fit, a lower bound that only grows as the graph is built.
fn qcd_order_at_most(max: usize) -> impl Fn(&View<'_>) -> bool + Clone + Send + Sync {
    const ORDER: [usize; 6] = [1, 1, 2, 1, 0, 0];
    move |view| {
        let lower: usize = (0..view.n_vertices())
            .map(|v| match view.rules(v) {
                Some(rules) => rules.iter().map(|&r| ORDER[r as usize]).min().unwrap(),
                None if view.degree(v) == 4 => 2,
                None => 0,
            })
            .sum();
        lower <= max
    }
}

/// Graphs of a run with `filter`, and the number of states the search visited.
fn run_counting(
    generator: GraphGenerator<u8, u8>,
    filter: impl Fn(&View<'_>) -> bool + Send + Sync + 'static,
) -> (Vec<Diagram<u8, u8>>, usize) {
    let states = Arc::new(AtomicUsize::new(0));
    let visits = states.clone();
    let diagrams = run(generator.filter(move |view| {
        visits.fetch_add(1, Ordering::Relaxed);
        filter(view)
    }));
    (diagrams, states.load(Ordering::Relaxed))
}

/// Pruning with a hereditary filter emits the graphs that the filter accepts when only asked about
/// complete topologies and complete graphs, in the same order, and the heredity check accepts the
/// filter. Returns the states visited with and without pruning.
fn assert_prunes_like_post_filter(
    generator: impl Fn() -> GraphGenerator<u8, u8>,
    filter: impl Fn(&View<'_>) -> bool + Clone + Send + Sync + 'static,
) -> (usize, usize) {
    let (pruned, pruned_states) = run_counting(generator(), filter.clone());
    let post = filter.clone();
    let (post_filtered, states) = run_counting(generator(), move |view| match view.stage() {
        Stage::Topology => !view.is_topology_complete() || post(view),
        Stage::Fields => true,
        Stage::Complete => post(view),
    });
    let checked = run(generator().filter(filter.clone()).check_heredity(true));
    let parallel = run(generator().filter(filter).parallel(true));
    assert_eq!(canonical_forms(&pruned), canonical_forms(&post_filtered));
    assert_eq!(canonical_forms(&checked), canonical_forms(&post_filtered));
    assert_eq!(canonical_forms(&parallel), canonical_forms(&post_filtered));
    (pruned_states, states)
}

#[test]
fn hereditary_filters_prune_to_the_post_filtered_graphs() {
    let gggg = || qcd(gluons(&[1, 2, 3, 4])).loops(1..=1);
    let ee_qq = |loops| {
        move || {
            let (e, q) = (Stub::incoming(ELECTRON), Stub::incoming(QUARK));
            let legs = vec![
                (1, e.clone()),
                (2, e.flipped()),
                (3, q.clone()),
                (4, q.flipped()),
            ];
            qcd_qed(legs, true).loops(loops..=loops)
        }
    };
    // stage 1: self-energies on legs are cut before the topology is complete
    let (pruned, all) = assert_prunes_like_post_filter(gggg, no_leg_self_energies);
    assert!(pruned < all, "{pruned} < {all}");
    // stage 1 for the four-gluon vertex, stage 2 for the others
    let (pruned, all) = assert_prunes_like_post_filter(gggg, without_rule(2));
    assert!(pruned < all, "{pruned} < {all}");
    assert_prunes_like_post_filter(gggg, without_rule(0));
    // stage 2 only, edge by edge
    let (pruned, all) = assert_prunes_like_post_filter(gggg, without_field(GHOST));
    assert!(pruned < all, "{pruned} < {all}");
    let ghostless = run(gggg())
        .into_iter()
        .filter(|d| d.graph.iter_edges().all(|(_, _, f)| *f.data != GHOST))
        .collect::<Vec<_>>();
    assert_eq!(
        canonical_forms(&run(gggg().filter(without_field(GHOST)))),
        canonical_forms(&ghostless)
    );
    // stages 1 and 2
    let (pruned, all) = assert_prunes_like_post_filter(ee_qq(2), qcd_order_at_most(0));
    assert!(pruned < all, "{pruned} < {all}");
    for max in 1..=4 {
        assert_prunes_like_post_filter(ee_qq(2), qcd_order_at_most(max));
    }
    // a conjunction of hereditary filters is hereditary
    let (bridges, order) = (no_leg_self_energies, qcd_order_at_most(2));
    assert_prunes_like_post_filter(ee_qq(2), move |view| bridges(view) && order(view));
}

/// A filter that vetoes the eventual vertex rule while it is undecided in stage 2.
fn vetoes_undecided_rules(view: &View<'_>) -> bool {
    view.stage() != Stage::Fields
        || (0..view.n_vertices()).all(|v| view.rules(v).is_none_or(|rules| rules != [1]))
}

/// Rejects the partial topologies, which are completed by every topology with two vertices or
/// more.
fn rejects_partial_topologies(view: &View<'_>) -> bool {
    view.stage() != Stage::Topology || view.is_topology_complete()
}

fn gggg_tree() -> GraphGenerator<u8, u8> {
    qcd(gluons(&[1, 2, 3, 4])).loops(0..=0)
}

#[test]
fn filters_that_reject_too_early_lose_graphs() {
    // the s, t and u channels and the contact graph; the rule veto keeps only the contact graph
    assert_eq!(run(gggg_tree()).len(), 4);
    assert_eq!(run(gggg_tree().filter(vetoes_undecided_rules)).len(), 1);
    // the four distinct gluons are labelled after the contact vertex is placed, so the contact
    // graph has partial topologies too; with two pairs of equal labels, it is placed at once
    assert_eq!(run(gggg_tree().filter(rejects_partial_topologies)).len(), 0);
    let pairs = qcd(gluons(&[1, 1, 2, 2])).loops(0..=0);
    assert_eq!(run(pairs.filter(rejects_partial_topologies)).len(), 1);
}

#[test]
#[should_panic(expected = "not hereditary")]
fn heredity_check_reports_rejected_partial_labellings() {
    run(gggg_tree()
        .filter(vetoes_undecided_rules)
        .check_heredity(true));
}

#[test]
#[should_panic(expected = "not hereditary")]
fn heredity_check_reports_rejected_partial_topologies() {
    run(gggg_tree()
        .filter(rejects_partial_topologies)
        .check_heredity(true));
}

#[test]
#[should_panic(expected = "not hereditary")]
fn heredity_check_follows_parallel_subtrees() {
    run(gggg_tree()
        .filter(rejects_partial_topologies)
        .check_heredity(true)
        .parallel(true));
}

/// The invariants every view satisfies, at every stage, with `n_legs` legs in the process.
fn check_view(view: &View<'_>, n_legs: usize) {
    let n = view.n_vertices();
    let edges: Vec<_> = view.edges().collect();
    let legs: usize = (0..n).map(|v| view.legs(v).count()).sum();
    // stubs are legs, edge ends or open
    let ends = |v| {
        edges
            .iter()
            .map(|e| e.vertices.iter().filter(|&&w| w == v).count())
            .sum::<usize>()
    };
    assert!((0..n).all(|v| view.degree(v) == view.legs(v).count() + ends(v) + view.open_stubs(v)));
    assert_eq!(legs + view.unplaced_legs(), n_legs);
    assert_eq!(edges.len() + 1, n + view.loops());
    let graph = view.graph();
    assert_eq!((graph.n_live(), graph.n_edges()), (n, edges.len() + legs));
    // the placed vertices are connected, pairs exactly when adjacent, the empty set never
    assert!(n == 0 || view.is_connected((0..n).collect()));
    assert!(!view.is_connected(VertexSet::default()));
    for u in 0..n {
        for v in u + 1..n {
            let pair = [u, v].into_iter().collect();
            assert_eq!(view.is_connected(pair), view.multiplicity(u, v) > 0);
        }
    }
    for bridge in view.bridges() {
        let [outside, inside] = bridge.ends;
        assert!(bridge.side.contains(inside) && !bridge.side.contains(outside));
        assert_eq!(bridge.side.iter().collect::<VertexSet>(), bridge.side);
        assert!(view.is_connected(bridge.side));
        assert!(bridge.side.iter().all(|v| view.open_stubs(v) == 0));
        // the bridge is the only edge leaving the side
        let leaving = |v| {
            (0..n)
                .filter(|&w| !bridge.side.contains(w))
                .map(|w| view.multiplicity(v, w))
                .sum::<usize>()
        };
        assert!(bridge
            .side
            .iter()
            .all(|v| leaving(v) == usize::from(v == inside)));
    }
    match view.stage() {
        Stage::Topology => assert!(edges.iter().all(|e| e.stub.is_none())),
        Stage::Fields => {
            assert!(view.is_topology_complete());
            // fields are assigned in edge order; a vertex fits rules once its stubs are known
            let assigned = edges.iter().take_while(|e| e.stub.is_some()).count();
            assert!(assigned > 0 && edges[assigned..].iter().all(|e| e.stub.is_none()));
            for v in 0..n {
                let known = edges[assigned..].iter().all(|e| !e.vertices.contains(&v));
                assert_eq!(view.rules(v).is_some(), known);
            }
        }
        Stage::Complete => {
            assert!(edges.iter().all(|e| e.stub.is_some()));
            assert!((0..n).all(|v| view.rules(v).unwrap().len() == 1));
            let complete = view.complete_graph().unwrap();
            assert_eq!(complete.n_edges(), graph.n_edges());
        }
    }
}

#[test]
fn views_describe_the_search_states() {
    let quarks = || {
        vec![
            (1, Stub::undirected(GLUON)),
            (2, Stub::incoming(QUARK)),
            (3, Stub::outgoing(QUARK)),
        ]
    };
    // four distinct gluons are labelled after their topology is shaped: until then they are open
    // stubs and unplaced legs
    for externals in [quarks(), gluons(&[1, 2, 3, 4])] {
        let stages = Arc::new(std::sync::Mutex::new(HashSet::new()));
        let seen = stages.clone();
        let n_legs = externals.len();
        let generator = qcd(externals.clone()).loops(1..=1).filter(move |view| {
            // some legs placed and some not: for the gluons, only while labelling
            let placed = (0..view.n_vertices()).any(|v| view.legs(v).count() > 0);
            let partial = view.stage() == Stage::Topology && placed && view.unplaced_legs() > 0;
            seen.lock().unwrap().insert((view.stage(), partial));
            check_view(view, n_legs);
            true
        });
        let filtered = run(generator);
        let stages = stages.lock().unwrap();
        assert_eq!(
            stages.iter().map(|(s, _)| s).collect::<HashSet<_>>().len(),
            3
        );
        assert!(stages.contains(&(Stage::Topology, true)));
        let unfiltered = run(qcd(externals).loops(1..=1));
        assert_eq!(canonical_forms(&filtered), canonical_forms(&unfiltered));
        check_canonical(&filtered);
    }
}

#[test]
fn views_of_complete_graphs_match_the_emitted_graphs() {
    // charged fields are directed, and both conversions of a complete view agree
    let electron = Stub::incoming(ELECTRON);
    let legs = vec![
        (1, electron.clone()),
        (2, electron.flipped()),
        (3, Stub::undirected(PHOTON)),
    ];
    let generator = qed(legs.clone()).loops(1..=1).filter(|view| {
        if let Some(complete) = view.complete_graph() {
            for e in view.edges() {
                let stub = e.stub.unwrap();
                assert_eq!(stub.flow.is_some(), *stub.field == ELECTRON);
            }
            let graph = view.graph();
            assert!(graph.iter_nodes().all(|(_, _, &rule)| rule == Some(0)));
            assert!(graph
                .iter_edges()
                .zip(complete.iter_edges())
                .all(|((_, _, a), (_, _, b))| a.orientation == b.orientation
                    && *a.data == Some(*b.data)));
        }
        true
    });
    let diagrams = run(generator);
    assert_eq!(
        canonical_forms(&diagrams),
        canonical_forms(&run(qed(legs).loops(1..=1)))
    );
    check_canonical(&diagrams);
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(48))]

    /// Conjunctions of hereditary filters prune to the post-filtered graphs on processes of QCD
    /// with photons and electrons.
    #[test]
    fn pruning_matches_post_filtering(
        process in 0usize..5,
        vetoed_fields in proptest::collection::vec(prop_oneof![Just(QUARK), Just(GHOST), Just(PHOTON)], 0..2),
        vetoed_rules in proptest::collection::vec(0u32..6, 0..2),
        max_order in proptest::option::of(0usize..4),
        self_energies in any::<bool>(),
    ) {
        let (g, q, e, a) = (
            Stub::undirected(GLUON),
            Stub::incoming(QUARK),
            Stub::incoming(ELECTRON),
            Stub::undirected(PHOTON),
        );
        let (legs, loops) = [
            (vec![g.clone(), g.clone(), g.clone(), g.clone()], 1),
            (vec![e.clone(), e.clone().flipped(), q.clone(), q.clone().flipped()], 2),
            (vec![q.clone(), q.clone().flipped(), g.clone(), a.clone()], 1),
            (vec![e.clone(), e.clone().flipped(), a.clone(), a.clone()], 2),
            (vec![g.clone(), g], 2),
        ][process].clone();
        let generator = move || {
            let legs = legs.iter().cloned().enumerate().map(|(i, s)| (i as u8, s)).collect();
            qcd_qed(legs, true).loops(loops..=loops)
        };
        let filter = move |view: &View<'_>| {
            vetoed_fields.iter().all(|&f| without_field(f)(view))
                && vetoed_rules.iter().all(|&r| without_rule(r)(view))
                && max_order.is_none_or(|max| qcd_order_at_most(max)(view))
                && (!self_energies || no_leg_self_energies(view))
        };
        assert_prunes_like_post_filter(generator, filter);
    }
}

/// A budget as the tests state it: a weight per signature or a count, with an inclusive range.
#[derive(Clone, Debug)]
enum Budget {
    Weight(Vec<i32>, [i32; 2]),
    Count(Countable<u8>, [usize; 2]),
}

/// The paired edges of a graph as `(node, node, field)` and its legs as `(node, field)`.
type Edges = (Vec<(usize, usize, u8)>, Vec<(usize, u8)>);

fn edges_and_legs(g: &HedgeGraph<u8, usize, Option<u8>>) -> Edges {
    let (mut edges, mut legs) = (vec![], vec![]);
    for (pair, _, data) in g.iter_edges() {
        match pair {
            HedgePair::Paired { source, sink } => {
                edges.push((g.node_id(source).0, g.node_id(sink).0, *data.data))
            }
            HedgePair::Unpaired { hedge, .. } => legs.push((g.node_id(hedge).0, *data.data)),
            HedgePair::Split { .. } => unreachable!("generated graphs have no split edges"),
        }
    }
    (edges, legs)
}

/// Edge `i` is a bridge: its ends are disconnected without it.
fn is_bridge(n: usize, edges: &[(usize, usize, u8)], i: usize) -> bool {
    let (u, v, _) = edges[i];
    let mut seen = vec![false; n];
    let mut stack = vec![u];
    seen[u] = true;
    while let Some(x) = stack.pop() {
        for (j, &(a, b, _)) in edges.iter().enumerate() {
            for (p, q) in [(a, b), (b, a)] {
                if j != i && p == x && !seen[q] {
                    seen[q] = true;
                    stack.push(q);
                }
            }
        }
    }
    !seen[v]
}

/// The closed lines of `fields` as feynkit counts closed fermion loops: the components of their
/// edges in which every vertex has two edge stubs of them and no leg of them; `None` if a vertex
/// has three or more stubs of them, legs included.
fn closed_lines(n: usize, (edges, legs): &Edges, fields: &[u8]) -> Option<usize> {
    let (mut stubs, mut leg_stubs) = (vec![0; n], vec![0; n]);
    let line: Vec<_> = edges.iter().filter(|e| fields.contains(&e.2)).collect();
    for &&(u, v, _) in &line {
        stubs[u] += 1;
        stubs[v] += 1;
    }
    for &(v, f) in legs {
        leg_stubs[v] += fields.contains(&f) as usize;
    }
    if (0..n).any(|v| stubs[v] + leg_stubs[v] > 2) {
        return None;
    }
    let mut component: Vec<usize> = (0..n).collect();
    let find = |c: &[usize], mut v: usize| {
        while c[v] != v {
            v = c[v];
        }
        v
    };
    for &&(u, v, _) in &line {
        let (a, b) = (find(&component, u), find(&component, v));
        component[a] = b;
    }
    let roots: HashSet<usize> = line.iter().map(|e| find(&component, e.0)).collect();
    Some(
        roots
            .into_iter()
            .filter(|&r| {
                (0..n)
                    .filter(|&v| find(&component, v) == r)
                    .all(|v| stubs[v] == 2 && leg_stubs[v] == 0)
            })
            .count(),
    )
}

impl Budget {
    fn apply(&self, generator: GraphGenerator<u8, u8>) -> GraphGenerator<u8, u8> {
        match self.clone() {
            Budget::Weight(weights, [min, max]) => generator.weight(weights, min..=max),
            Budget::Count(count, [min, max]) => generator.count(count, min..=max),
        }
    }

    /// The emitted graph `d` meets the budget, computed on the graph alone.
    fn admits(&self, d: &Diagram<u8, u8>) -> bool {
        let g = &d.graph;
        let n = g.node_bound();
        let rules: Vec<usize> = g.iter_nodes().map(|(_, _, &rule)| rule).collect();
        let lines = edges_and_legs(g);
        let edges = &lines.0;
        match self {
            Budget::Weight(weights, [min, max]) => {
                let sum: i32 = rules.iter().map(|&r| weights[r]).sum();
                (*min..=*max).contains(&sum)
            }
            Budget::Count(count, [min, max]) => {
                let value = match count {
                    Countable::Vertices(signatures) => {
                        rules.iter().filter(|r| signatures.contains(r)).count()
                    }
                    Countable::Propagators(fields) => {
                        edges.iter().filter(|e| fields.contains(&e.2)).count()
                    }
                    Countable::Chords(fields) => (0..edges.len())
                        .filter(|&i| fields.contains(&edges[i].2) && !is_bridge(n, edges, i))
                        .count(),
                    Countable::ClosedLines(fields) => match closed_lines(n, &lines, fields) {
                        Some(closed) => closed,
                        None => return true,
                    },
                };
                (*min..=*max).contains(&value)
            }
        }
    }
}

/// Budgets prune to the graphs of the unbudgeted run that meet them, in its order, sequentially
/// and in parallel. Returns the states the filter sees with and without the budgets.
fn assert_budgets_prune_like_post_filter(
    generator: impl Fn() -> GraphGenerator<u8, u8>,
    budgets: &[Budget],
) -> (usize, usize) {
    assert_filtered_budgets_prune_like_post_filter(generator, |_| true, budgets)
}

/// [`assert_budgets_prune_like_post_filter`] with a hereditary `filter` pruning both runs.
fn assert_filtered_budgets_prune_like_post_filter(
    generator: impl Fn() -> GraphGenerator<u8, u8>,
    filter: impl Fn(&View<'_>) -> bool + Clone + Send + Sync + 'static,
    budgets: &[Budget],
) -> (usize, usize) {
    let budgeted = || budgets.iter().fold(generator(), |g, b| b.apply(g));
    let (pruned, pruned_states) = run_counting(budgeted(), filter.clone());
    let (all, states) = run_counting(generator(), filter.clone());
    let post: Vec<_> = all
        .into_iter()
        .filter(|d| budgets.iter().all(|b| b.admits(d)))
        .collect();
    assert_eq!(
        canonical_forms(&pruned),
        canonical_forms(&post),
        "{budgets:?}"
    );
    let parallel = run(budgeted().filter(filter).parallel(true));
    assert_eq!(
        canonical_forms(&parallel),
        canonical_forms(&post),
        "{budgets:?}"
    );
    (pruned_states, states)
}

/// e+ e- -> q q~ in [`qcd_qed`] at `loops` loops, and its QCD and QED orders per signature.
fn ee_qq(loops: usize) -> impl Fn() -> GraphGenerator<u8, u8> {
    move || {
        let (e, q) = (Stub::incoming(ELECTRON), Stub::incoming(QUARK));
        let legs = vec![
            (1, e.clone()),
            (2, e.flipped()),
            (3, q.clone()),
            (4, q.flipped()),
        ];
        qcd_qed(legs, true).loops(loops..=loops)
    }
}
const QCD_ORDER: [i32; 6] = [1, 1, 2, 1, 0, 0];
const QED_ORDER: [i32; 6] = [0, 0, 0, 0, 1, 1];

#[test]
fn coupling_orders_prune_to_the_post_filtered_graphs() {
    let order = |weights: [i32; 6], min, max| Budget::Weight(weights.to_vec(), [min, max]);
    for range in [[0, 0], [0, 1], [2, 2], [2, i32::MAX], [1, 3], [3, 2]] {
        let (pruned, all) = assert_budgets_prune_like_post_filter(
            ee_qq(2),
            &[order(QCD_ORDER, range[0], range[1])],
        );
        assert!(pruned < all, "{range:?}: {pruned} < {all}");
    }
    // both orders, which fix the loop order: QCD + QED = 2 L + 2 here
    assert_budgets_prune_like_post_filter(
        ee_qq(2),
        &[order(QCD_ORDER, 2, 2), order(QED_ORDER, 4, 4)],
    );
    // signed weights: QCD less QED
    let signed: Vec<i32> = (0..6).map(|s| QCD_ORDER[s] - QED_ORDER[s]).collect();
    for range in [[-2, -2], [-1, 0], [0, 6]] {
        assert_budgets_prune_like_post_filter(ee_qq(1), &[Budget::Weight(signed.clone(), range)]);
    }
}

#[test]
fn counts_prune_to_the_post_filtered_graphs() {
    let gggg = || qcd(gluons(&[1, 2, 3, 4])).loops(1..=1);
    let count = |count, min, max| Budget::Count(count, [min, max]);
    let cases = [
        count(Countable::Propagators(vec![GHOST]), 0, 0),
        count(Countable::Propagators(vec![GHOST, QUARK]), 2, 3),
        count(Countable::Chords(vec![GLUON]), 1, 2),
        count(Countable::Chords(vec![GLUON]), 0, 0),
        count(Countable::Vertices(vec![2]), 1, 1),
        count(Countable::ClosedLines(vec![QUARK]), 1, 1),
        count(Countable::ClosedLines(vec![QUARK, GHOST]), 0, 0),
    ];
    for budget in &cases {
        let (pruned, all) =
            assert_budgets_prune_like_post_filter(gggg, std::slice::from_ref(budget));
        assert!(pruned < all, "{budget:?}: {pruned} < {all}");
    }
    // the electron line is open, the quark lines may close
    for (min, max) in [(0, 0), (1, 1), (0, 1)] {
        assert_budgets_prune_like_post_filter(
            ee_qq(2),
            &[count(
                Countable::ClosedLines(vec![ELECTRON, QUARK]),
                min,
                max,
            )],
        );
    }
    // the contact graph of the gluon tree has no edge: its legs alone decide its counts
    let tree = || qcd(gluons(&[1, 2, 3, 4])).loops(0..=0);
    for budget in [
        count(Countable::Propagators(vec![GLUON]), 1, 1),
        count(Countable::ClosedLines(vec![QUARK]), 1, 1),
        count(Countable::Vertices(vec![2]), 0, 0),
    ] {
        assert_budgets_prune_like_post_filter(tree, &[budget]);
    }
    // a conjunction of counts and orders
    assert_budgets_prune_like_post_filter(
        ee_qq(2),
        &[
            count(Countable::Chords(vec![PHOTON]), 0, 1),
            Budget::Weight(QCD_ORDER.to_vec(), [0, 2]),
            count(Countable::Propagators(vec![QUARK]), 2, 4),
        ],
    );
}

#[test]
fn branching_lines_are_counted_on_complete_graphs() {
    // a four-quark vertex lets quark lines branch: graphs using it are not constrained by the
    // closed quark lines, the others are, and only complete graphs can tell them apart
    let (g, q) = (Stub::undirected(GLUON), Stub::incoming(QUARK));
    let generator = move || {
        let signatures = [
            vec![q.clone().flipped(), q.clone(), g.clone()],
            vec![g.clone(); 3],
            vec![
                q.clone().flipped(),
                q.clone(),
                q.clone().flipped(),
                q.clone(),
            ],
        ];
        let legs = [
            (1, q.clone()),
            (2, q.clone().flipped()),
            (3, g.clone()),
            (4, g.clone()),
        ];
        GraphGenerator::new(legs, signatures).loops(2..=2)
    };
    for (min, max) in [(0, 0), (1, 1)] {
        let budget = Budget::Count(Countable::ClosedLines(vec![QUARK]), [min, max]);
        assert_budgets_prune_like_post_filter(generator.clone(), &[budget]);
    }
}

#[test]
fn lines_follow_only_the_field_assignments_the_filter_keeps() {
    // the filters reject counted edges: the ghost edges of the quark and ghost lines, and the
    // third quark edge
    let few_quarks = |view: &View<'_>| {
        let quark = |e: &EdgeView<'_, u8>| e.stub.as_ref().is_some_and(|s| *s.field == QUARK);
        view.edges().filter(quark).count() <= 2
    };
    for (labels, loops) in [
        (&[1, 2, 3, 4][..], 2),
        (&[1, 1, 1, 1][..], 2),
        (&[1, 2][..], 3),
    ] {
        let generator = move || qcd(gluons(labels)).loops(loops..=loops);
        for range in [[0, 0], [1, 1], [2, 2], [0, 1]] {
            let lines = |fields| [Budget::Count(Countable::ClosedLines(fields), range)];
            let ghosts = lines(vec![QUARK, GHOST]);
            assert_filtered_budgets_prune_like_post_filter(
                generator,
                without_field(GHOST),
                &ghosts,
            );
            assert_filtered_budgets_prune_like_post_filter(
                generator,
                few_quarks,
                &lines(vec![QUARK]),
            );
        }
    }
}

#[test]
fn single_diagram_models_meet_budgets() {
    // phi^3 + phi^4 has one diagram per topology
    let two_point = || phi(&[3, 4], &[1, 2]).loops(3..=3);
    for budget in [
        Budget::Count(Countable::Vertices(vec![1]), [1, 1]),
        Budget::Weight(vec![1, 2], [5, 6]),
        Budget::Count(Countable::Chords(vec![PHI]), [4, 5]),
        Budget::Count(Countable::Propagators(vec![PHI]), [0, 4]),
    ] {
        assert_budgets_prune_like_post_filter(two_point, &[budget]);
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(48))]

    /// Random weights (of either sign) and counts prune to the post-filtered graphs on processes
    /// of QCD with photons and electrons.
    #[test]
    fn budgets_match_post_filtering(
        process in 0usize..5,
        weights in proptest::collection::vec(
            (proptest::collection::vec(-1i32..=2, 6), -2i32..5, 0i32..4), 0..3),
        counts in proptest::collection::vec(
            (0usize..4,
             proptest::collection::vec(
                prop_oneof![Just(GLUON), Just(QUARK), Just(GHOST), Just(PHOTON), Just(ELECTRON)],
                1..3),
             0usize..3, 0usize..3),
            0..3),
    ) {
        let (g, q, e, a) = (
            Stub::undirected(GLUON),
            Stub::incoming(QUARK),
            Stub::incoming(ELECTRON),
            Stub::undirected(PHOTON),
        );
        let (legs, loops) = [
            (vec![g.clone(), g.clone(), g.clone(), g.clone()], 1),
            (vec![e.clone(), e.clone().flipped(), q.clone(), q.clone().flipped()], 2),
            (vec![q.clone(), q.clone().flipped(), g.clone(), a.clone()], 1),
            (vec![e.clone(), e.clone().flipped(), a.clone(), a.clone()], 2),
            (vec![g.clone(), g], 2),
        ][process].clone();
        let generator = move || {
            let legs = legs.iter().cloned().enumerate().map(|(i, s)| (i as u8, s)).collect();
            qcd_qed(legs, true).loops(loops..=loops)
        };
        let mut budgets: Vec<Budget> = weights
            .into_iter()
            .map(|(w, min, width)| Budget::Weight(w, [min, min + width]))
            .collect();
        budgets.extend(counts.into_iter().map(|(kind, fields, min, width)| {
            let count = match kind {
                0 => Countable::Propagators(fields),
                1 => Countable::Chords(fields),
                2 => Countable::ClosedLines(fields),
                _ => Countable::Vertices(fields.iter().map(|&f| f as usize % 6).collect()),
            };
            Budget::Count(count, [min, min + width])
        }));
        assert_budgets_prune_like_post_filter(generator, &budgets);
    }
}

/// The states that the filter of `generator` sees, and the complete topologies it reports.
fn searched(generator: GraphGenerator<u8, u8>) -> (usize, usize) {
    let topologies = Arc::new(AtomicUsize::new(0));
    let seen = topologies.clone();
    let generator = generator.progress(move |progress| {
        seen.store(progress.topologies, Ordering::Relaxed);
        true
    });
    let (diagrams, states) = run_counting(generator, |_| true);
    assert!(diagrams.is_empty());
    (states, topologies.load(Ordering::Relaxed))
}

#[test]
fn requests_proven_empty_are_not_searched() {
    let e = Stub::incoming(ELECTRON);
    // charge: two incoming electrons
    let charge = qed(vec![(1, e.clone()), (2, e.clone())]).loops(0..=2);
    assert_eq!(searched(charge), (0, 0));
    // parity: one leg of a field that only quartic vertices carry
    let parity = GraphGenerator::<u8, u8>::new(
        [(1, Stub::undirected(CHI)), (2, Stub::undirected(PHI))],
        [
            vec![Stub::undirected(PHI); 3],
            vec![
                Stub::undirected(CHI),
                Stub::undirected(CHI),
                Stub::undirected(PHI),
                Stub::undirected(PHI),
            ],
        ],
    )
    .loops(0..=2);
    assert_eq!(searched(parity), (0, 0));
    // weights: e+ e- -> e+ e- needs QED order 2 + 2 L
    let legs = vec![
        (1, e.clone()),
        (2, e.clone().flipped()),
        (3, e.clone()),
        (4, e.flipped()),
    ];
    let order = qed(legs.clone()).loops(0..=1).weight([1], 3..=3);
    assert_eq!(searched(order), (0, 0));
    let vertices = qed(legs.clone())
        .loops(1..=1)
        .count(Countable::Vertices(vec![0]), ..4);
    assert_eq!(searched(vertices), (0, 0));
    // an empty range
    let empty = qed(legs.clone())
        .loops(0..=1)
        .count(Countable::Propagators(vec![PHOTON]), 2..2);
    assert_eq!(searched(empty), (0, 0));
    // the loop orders that the weights admit are searched
    assert_eq!(
        canonical_forms(&run(qed(legs.clone()).loops(0..=2).weight([1], 4..=4))),
        canonical_forms(&run(qed(legs).loops(1..=1)))
    );
}

/// The graphs of `generator` with `filter`, the field assignments the filter sees and those it
/// rejects.
fn field_assignments(
    generator: GraphGenerator<u8, u8>,
    filter: impl Fn(&View<'_>) -> bool + Send + Sync + 'static,
) -> (Vec<Diagram<u8, u8>>, usize, usize) {
    let (seen, rejected) = (Arc::new(AtomicUsize::new(0)), Arc::new(AtomicUsize::new(0)));
    let (s, r) = (seen.clone(), rejected.clone());
    let diagrams = run(generator.filter(move |view| {
        let keeps = filter(view);
        if view.stage() == Stage::Fields {
            s.fetch_add(1, Ordering::Relaxed);
            r.fetch_add(!keeps as usize, Ordering::Relaxed);
        }
        keeps
    }));
    (
        diagrams,
        seen.load(Ordering::Relaxed),
        rejected.load(Ordering::Relaxed),
    )
}

#[test]
fn budgets_prune_field_assignments_before_the_filter() {
    // ee -> qq with QCD order at most `max`: the budget bounds each partial vertex by the signatures
    // it can still become, so the filter that bounds the order by the vertices whose stubs are all
    // assigned (and four-gluon vertices) never rejects what the budget shows it, and sees fewer
    // assignments than when it prunes alone
    for (loops, max) in [(1, 0), (2, 0), (2, 1), (2, 2)] {
        let filter = qcd_order_at_most(max);
        let budgeted = ee_qq(loops)().weight(QCD_ORDER, ..=max as i32);
        let (diagrams, seen, rejected) = field_assignments(budgeted, filter.clone());
        let (alone, seen_alone, _) = field_assignments(ee_qq(loops)(), filter);
        assert_eq!(canonical_forms(&diagrams), canonical_forms(&alone));
        assert_eq!(rejected, 0, "L{loops} QCD <= {max}");
        assert!(
            seen < seen_alone,
            "L{loops} QCD <= {max}: {seen} < {seen_alone}"
        );
    }
}

#[test]
fn budget_errors_and_limits() {
    let wrong = qed(vec![]).loops(1..=1).weight([1, 2], ..);
    assert_eq!(
        wrong.iter().err(),
        Some(GenerationError::Weights {
            signatures: 1,
            weights: 2
        })
    );
    // at most 64 edge and line counts, not counting those that every graph meets
    let e = Stub::incoming(ELECTRON);
    let legs = vec![
        (1, e.clone()),
        (2, e.clone().flipped()),
        (3, e.clone()),
        (4, e.flipped()),
    ];
    let photons = |n, max| {
        (0..n).fold(qed(legs.clone()).loops(1..=1), |g, _| {
            g.count(Countable::Propagators(vec![PHOTON]), 0..=max)
        })
    };
    assert_eq!(
        photons(65, 9).iter().err(),
        Some(GenerationError::TooManyCounts)
    );
    let plain = canonical_forms(&run(photons(0, 0)));
    assert_eq!(canonical_forms(&run(photons(64, 9))), plain);
    assert_eq!(canonical_forms(&run(photons(65, usize::MAX))), plain);
    // budgets that every graph meets are left out, the others kept
    let budgets = |generator: GraphGenerator<u8, u8>| {
        let request = Compiled::request(&generator, &[PHOTON, ELECTRON]).unwrap();
        (
            request.weights.len(),
            request.edges.len(),
            request.lines.len(),
        )
    };
    let vacuous = qed(legs.clone())
        .weight([1], 0..)
        .weight([-1], ..=2)
        .count(Countable::Vertices(vec![0]), ..)
        .count(Countable::ClosedLines(vec![ELECTRON]), 0..);
    assert_eq!(budgets(vacuous), (0, 0, 0));
    let bounding = qed(legs)
        .weight([1], 0..=9)
        .weight([-1], 0..)
        .count(Countable::Vertices(vec![0]), 1..)
        .count(Countable::ClosedLines(vec![ELECTRON]), ..=1);
    assert_eq!(budgets(bounding), (3, 0, 1));
}

#[test]
fn budgets_do_not_outlive_their_generation() {
    // decorators return to a pool of the thread and serve later generations
    let gggg = || qcd(gluons(&[1, 2, 3, 4])).loops(1..=1);
    let plain = canonical_forms(&run(gggg()));
    let lines = gggg().count(Countable::ClosedLines(vec![QUARK]), 1..=1);
    assert!(run(lines).len() < plain.len());
    assert_eq!(canonical_forms(&run(gggg())), plain);
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

    // progress sees every complete topology, so it stops a run whose graphs are all filtered out
    let seen = topologies.clone();
    let generator = phi(&[3, 4], &[1, 2])
        .loops(3..=3)
        .filter(|view| !view.is_complete())
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
            .map(|d| (d.loops, d.graph.n_live(), d.automorphism_count))
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
            phi(&[3], &[1, 2, 3, 4, 5, 6, 7]).loops(0..=0),
            // several loop orders, of which only the first starts with the oldest subtree alone
            phi(&[3, 4], &[1, 2, 3, 4]).loops(0..=2),
        ]
    };
    for (sequential, parallel) in generators().into_iter().zip(generators()) {
        let order =
            |d: Vec<Diagram<u8, u8>>| d.into_iter().map(|d| d.canonical_form).collect::<Vec<_>>();
        let expected = order(run(sequential));
        let parallel = parallel.parallel(true);
        for threads in [1, 2, 4, 8, 32] {
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
        canon::CanonicalForm,
        involution::{Flow, Orientation},
    };

    type Plain = Graph<u8, u8>;
    type Key = CanonicalForm<usize, u8, Option<u8>>;

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
            Self::forms(self.generator())
        }

        fn forms(generator: GraphGenerator<u8, u8>) -> HashMap<Key, u128> {
            generator
                .iter()
                .unwrap()
                .map(|d| (d.canonical_form, d.automorphism_count.unwrap()))
                .collect()
        }

        fn generator(&self) -> GraphGenerator<u8, u8> {
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
        }

        /// Graphica's graphs that `keep` accepts.
        fn graphica(
            &self,
            keep: impl Fn(&HedgeGraph<u8, usize, Option<u8>>) -> bool,
        ) -> HashMap<Key, u128> {
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
                .filter_map(|(g, aut)| {
                    let graph = self.to_linnet(&g);
                    let form = graph.canonical_form_by(|&s| s, |&f| f, |&h| h);
                    keep(&graph).then(|| (form, aut.to_string().parse().unwrap()))
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

    /// QCD signatures in the order of [`qcd`](super::qcd).
    fn qcd_signatures() -> Vec<Vec<HalfEdge<u8>>> {
        let gluon = HalfEdge::undirected(GLUON);
        vec![
            vec![HalfEdge::outgoing(QUARK), HalfEdge::incoming(QUARK), gluon],
            vec![gluon; 3],
            vec![gluon; 4],
            vec![HalfEdge::outgoing(GHOST), HalfEdge::incoming(GHOST), gluon],
        ]
    }

    fn labelled(legs: &[HalfEdge<u8>]) -> Vec<(u8, HalfEdge<u8>)> {
        legs.iter()
            .enumerate()
            .map(|(i, &h)| (i as u8 + 1, h))
            .collect()
    }

    #[test]
    fn sets_and_symmetry_factors_match_graphica() {
        let phi = HalfEdge::undirected(PHI);
        let photon = HalfEdge::undirected(PHOTON);
        let gluon = HalfEdge::undirected(GLUON);
        let electron = HalfEdge::incoming(ELECTRON);
        let qed = vec![vec![electron, electron.flip(), photon]];
        let qcd = qcd_signatures();
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
                signatures: qed.clone(),
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
                signatures: qcd.clone(),
                loops: 2,
                max_bridges: Some(0),
                zero_flow_edges: false,
            },
            // distinct legs with a common stub: labelled on shaped topologies
            Case {
                externals: labelled(&[phi; 7]),
                signatures: vec![vec![phi; 3]],
                loops: 0,
                max_bridges: None,
                zero_flow_edges: false,
            },
            Case {
                externals: labelled(&[phi; 6]),
                signatures: vec![vec![phi; 3], vec![phi; 4]],
                loops: 0,
                max_bridges: None,
                zero_flow_edges: false,
            },
            Case {
                externals: labelled(&[phi; 5]),
                signatures: vec![vec![phi; 3]],
                loops: 2,
                max_bridges: None,
                zero_flow_edges: false,
            },
            Case {
                externals: [1, 1, 2, 3, 4].map(|l| (l, gluon)).to_vec(),
                signatures: qcd,
                loops: 1,
                max_bridges: None,
                zero_flow_edges: false,
            },
            Case {
                externals: labelled(&[electron.flip(), electron, photon, photon, photon]),
                signatures: qed.clone(),
                loops: 1,
                max_bridges: None,
                zero_flow_edges: false,
            },
        ];
        for (i, case) in cases.iter().enumerate() {
            let (native, graphica) = (case.native(), case.graphica(|_| true));
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
            assert_eq!(all.native(), all.graphica(|_| true), "L{loops}");
            for restricted in [case(loops, None, false), case(loops, Some(0), true)] {
                assert!(restricted.graphica(|_| true).is_empty(), "L{loops}");
                let native = restricted.native();
                assert!(!native.is_empty(), "L{loops}");
                assert!(native.keys().all(|k| all.native().contains_key(k)));
            }
        }
    }

    /// Pruning with a particle or vertex veto gives the Graphica graphs without that particle or
    /// vertex.
    #[test]
    fn pruned_sets_match_filtered_graphica() {
        let gluon = HalfEdge::undirected(GLUON);
        let quark = HalfEdge::incoming(QUARK);
        for externals in [
            labelled(&[gluon, gluon, gluon, gluon]),
            labelled(&[quark, quark.flip(), gluon, gluon]),
        ] {
            let case = Case {
                externals,
                signatures: qcd_signatures(),
                loops: 1,
                max_bridges: None,
                zero_flow_edges: false,
            };
            let native = Case::forms(case.generator().filter(without_field(GHOST)));
            let graphica = case.graphica(|g| g.iter_edges().all(|(_, _, f)| *f.data != GHOST));
            assert!(!native.is_empty() && native.len() < case.native().len());
            assert_eq!(native, graphica);
            let native = Case::forms(case.generator().filter(without_rule(2)));
            let graphica = case.graphica(|g| g.iter_nodes().all(|(_, _, &rule)| rule != 2));
            assert!(!native.is_empty() && native.len() < case.native().len());
            assert_eq!(native, graphica);
        }
    }
}
