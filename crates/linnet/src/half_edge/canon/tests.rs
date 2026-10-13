use rand::{rngs::SmallRng, seq::SliceRandom, Rng, SeedableRng};

use super::{ignore, Canonization};
use crate::half_edge::{
    builder::HedgeGraphBuilder,
    involution::{Flow, Hedge, Orientation},
    remap::{HedgeRemap, NodeRemap},
    HedgeGraph, NodeIndex,
};
use crate::permutation::Permutation;

/// Edge data, node data and half-edge data are all keys.
type G = HedgeGraph<u8, u8, u8>;

#[derive(Clone, Copy, Debug)]
enum Element {
    Edge {
        a: usize,
        ka: u8,
        b: usize,
        kb: u8,
        data: u8,
        orientation: Orientation,
    },
    Leg {
        node: usize,
        key: u8,
        data: u8,
        orientation: Orientation,
        flow: Flow,
    },
}

/// A graph description that can be rebuilt under random relabellings.
#[derive(Clone, Debug)]
struct Spec {
    nodes: Vec<u8>,
    elements: Vec<Element>,
}

impl Spec {
    fn new(nodes: usize) -> Self {
        Spec {
            nodes: vec![0; nodes],
            elements: vec![],
        }
    }

    fn edge(
        mut self,
        a: usize,
        ka: u8,
        b: usize,
        kb: u8,
        data: u8,
        orientation: Orientation,
    ) -> Self {
        self.elements.push(Element::Edge {
            a,
            ka,
            b,
            kb,
            data,
            orientation,
        });
        self
    }

    fn undirected(self, a: usize, ka: u8, b: usize, kb: u8) -> Self {
        self.edge(a, ka, b, kb, 0, Orientation::Undirected)
    }

    fn leg(mut self, node: usize, key: u8) -> Self {
        self.elements.push(Element::Leg {
            node,
            key,
            data: 0,
            orientation: Orientation::Undirected,
            flow: Flow::Source,
        });
        self
    }

    fn build(&self) -> G {
        let mut builder = HedgeGraphBuilder::new();
        let nodes: Vec<NodeIndex> = self.nodes.iter().map(|&c| builder.add_node(c)).collect();
        for &e in &self.elements {
            match e {
                Element::Edge {
                    a,
                    ka,
                    b,
                    kb,
                    data,
                    orientation,
                } => builder.add_edge(
                    nodes[a].add_data(ka),
                    nodes[b].add_data(kb),
                    data,
                    orientation,
                ),
                Element::Leg {
                    node,
                    key,
                    data,
                    orientation,
                    flow,
                } => builder.add_external_edge(nodes[node].add_data(key), data, orientation, flow),
            }
        }
        builder.build()
    }

    /// Random node order, element order and stored roles, with the same physical graph.
    fn relabel(&self, rng: &mut SmallRng) -> Spec {
        let mut perm: Vec<usize> = (0..self.nodes.len()).collect();
        perm.shuffle(rng);
        let mut nodes = vec![0; self.nodes.len()];
        for (v, &p) in perm.iter().enumerate() {
            nodes[p] = self.nodes[v];
        }
        let mut elements: Vec<Element> = self
            .elements
            .iter()
            .map(|&e| match e {
                Element::Edge {
                    a,
                    ka,
                    b,
                    kb,
                    data,
                    orientation,
                } if rng.gen_bool(0.5) => Element::Edge {
                    a: perm[b],
                    ka: kb,
                    b: perm[a],
                    kb: ka,
                    data,
                    orientation: orientation.reverse(),
                },
                Element::Edge {
                    a,
                    ka,
                    b,
                    kb,
                    data,
                    orientation,
                } => Element::Edge {
                    a: perm[a],
                    ka,
                    b: perm[b],
                    kb,
                    data,
                    orientation,
                },
                Element::Leg {
                    node,
                    key,
                    data,
                    orientation,
                    flow,
                } => {
                    let flip = rng.gen_bool(0.5);
                    Element::Leg {
                        node: perm[node],
                        key,
                        data,
                        orientation: if flip {
                            orientation.reverse()
                        } else {
                            orientation
                        },
                        flow: if flip { -flow } else { flow },
                    }
                }
            })
            .collect();
        elements.shuffle(rng);
        Spec { nodes, elements }
    }

    /// [`Self::relabel`], built, with the hedges then permuted at random as well (the builder
    /// numbers the hedges of each edge consecutively).
    fn relabelled_graph(&self, rng: &mut SmallRng) -> G {
        let mut g = self.relabel(rng).build();
        let mut map: Vec<usize> = (0..g.n_hedges()).collect();
        map.shuffle(rng);
        g.permute_hedges(&HedgeRemap::from_permutation(&Permutation::from_map(map)))
            .unwrap();
        g
    }

    fn random(rng: &mut SmallRng, max_nodes: usize, max_edges: usize, max_legs: usize) -> Spec {
        let n = rng.gen_range(1..=max_nodes);
        let node_colours = rng.gen_range(1..=2);
        let keys = rng.gen_range(1..=2);
        let data = rng.gen_range(1..=2);
        let orientation = |rng: &mut SmallRng| match rng.gen_range(0..4) {
            0 => Orientation::Default,
            1 => Orientation::Reversed,
            _ => Orientation::Undirected,
        };
        let mut spec = Spec {
            nodes: (0..n).map(|_| rng.gen_range(0..node_colours)).collect(),
            elements: vec![],
        };
        for _ in 0..rng.gen_range(0..=max_edges) {
            let o = orientation(rng);
            spec = spec.edge(
                rng.gen_range(0..n),
                rng.gen_range(0..keys),
                rng.gen_range(0..n),
                rng.gen_range(0..keys),
                rng.gen_range(0..data),
                o,
            );
        }
        for _ in 0..rng.gen_range(0..=max_legs) {
            spec.elements.push(Element::Leg {
                node: rng.gen_range(0..n),
                key: rng.gen_range(0..keys),
                data: rng.gen_range(0..data),
                orientation: orientation(rng),
                flow: if rng.gen_bool(0.5) {
                    Flow::Source
                } else {
                    Flow::Sink
                },
            });
        }
        spec
    }
}

/// Everything an isomorphism must preserve about a hedge, read through the public API.
type Colour = (u8, u8, Option<Flow>, bool);

fn colour(g: &G, h: Hedge) -> Colour {
    (
        *g.get_edge_data(h),
        g[h],
        g.superficial_hedge_orientation(h),
        g.is_dangling(h),
    )
}

/// The canonically relabelled graph: node colours by canonical node, and per canonical hedge its
/// canonical node, colour and canonical opposite.
type Relabelled = (Vec<u8>, Vec<(usize, Colour, usize)>);

fn relabelled(g: &G, c: &Canonization<&u8, &u8, &u8>) -> Relabelled {
    let mut nodes = vec![0; g.n_live()];
    for v in g.iter_node_ids() {
        nodes[c.canonical_node(v).0] = g[v];
    }
    let mut hedges: Vec<_> = (0..g.n_hedges())
        .map(Hedge)
        .map(|h| {
            (
                c.canonical_hedge(h).0,
                c.canonical_node(g.node_id(h)).0,
                colour(g, h),
                c.canonical_hedge(g.inv(h)).0,
            )
        })
        .collect();
    hedges.sort_unstable();
    let hedges = hedges.into_iter().map(|(_, v, c, o)| (v, c, o)).collect();
    (nodes, hedges)
}

/// Brute-force count of isomorphisms `a -> b` (with `pin` mapped as given), over hedge
/// bijections and colour-preserving bijections of isolated nodes.
fn count_isomorphisms(a: &G, b: &G, pin: Option<(Hedge, Hedge)>) -> u128 {
    struct State<'g> {
        a: &'g G,
        b: &'g G,
        image: Vec<usize>,
        used: Vec<bool>,
        node_map: Vec<usize>,
        node_used: Vec<bool>,
        pin: Option<(Hedge, Hedge)>,
    }
    fn rec(s: &mut State<'_>, x: usize) -> u128 {
        if x == s.a.n_hedges() {
            return 1;
        }
        let ix = s.a.inv(Hedge(x)).0;
        let mut total = 0;
        for y in 0..s.b.n_hedges() {
            if s.used[y]
                || colour(s.a, Hedge(x)) != colour(s.b, Hedge(y))
                || s.pin.is_some_and(|(p, q)| p.0 == x && q.0 != y)
            {
                continue;
            }
            let iy = s.b.inv(Hedge(y)).0;
            if (ix < x && s.image[ix] != iy) || (ix > x && s.used[iy]) {
                continue;
            }
            let (n, m) = (s.a.node_id(Hedge(x)).0, s.b.node_id(Hedge(y)).0);
            let fresh = s.node_map[n] == usize::MAX;
            if fresh && (s.node_used[m] || s.a[NodeIndex(n)] != s.b[NodeIndex(m)]) {
                continue;
            }
            if !fresh && s.node_map[n] != m {
                continue;
            }
            s.image[x] = y;
            s.used[y] = true;
            if fresh {
                s.node_map[n] = m;
                s.node_used[m] = true;
            }
            total += rec(s, x + 1);
            s.used[y] = false;
            if fresh {
                s.node_map[n] = usize::MAX;
                s.node_used[m] = false;
            }
        }
        total
    }
    if a.n_hedges() != b.n_hedges() || a.n_live() != b.n_live() {
        return 0;
    }
    let isolated = |g: &G| {
        let mut colours: Vec<u8> = g
            .iter_nodes()
            .filter(|(_, hedges, _)| hedges.len() == 0)
            .map(|(_, _, &c)| c)
            .collect();
        colours.sort_unstable();
        colours
    };
    let (ia, ib) = (isolated(a), isolated(b));
    if ia != ib {
        return 0;
    }
    let isolated_maps: u128 = ia
        .chunk_by(|x, y| x == y)
        .map(|run| (1..=run.len() as u128).product::<u128>())
        .product();
    let mut s = State {
        a,
        b,
        image: vec![0; a.n_hedges()],
        used: vec![false; b.n_hedges()],
        node_map: vec![usize::MAX; a.node_bound()],
        node_used: vec![false; b.node_bound()],
        pin,
    };
    rec(&mut s, 0) * isolated_maps
}

/// Group order from generators acting on `n` points, by deterministic Schreier-Sims.
fn schreier_sims_order(gens: &[Vec<usize>], n: usize) -> u128 {
    struct Level {
        base: usize,
        gens: Vec<Vec<usize>>,
        /// point -> u with u(base) = point
        transversal: Vec<Option<Vec<usize>>>,
    }
    fn compose(a: &[usize], b: &[usize]) -> Vec<usize> {
        b.iter().map(|&x| a[x]).collect()
    }
    fn inverse(a: &[usize]) -> Vec<usize> {
        let mut r = vec![0; a.len()];
        for (i, &x) in a.iter().enumerate() {
            r[x] = i;
        }
        r
    }
    fn is_identity(a: &[usize]) -> bool {
        a.iter().enumerate().all(|(i, &x)| i == x)
    }
    fn rebuild(l: &mut Level, n: usize) {
        l.transversal = vec![None; n];
        l.transversal[l.base] = Some((0..n).collect());
        let mut queue = vec![l.base];
        let mut i = 0;
        while i < queue.len() {
            let p = queue[i];
            i += 1;
            for g in &l.gens {
                let q = g[p];
                if l.transversal[q].is_none() {
                    l.transversal[q] = Some(compose(g, l.transversal[p].as_ref().unwrap()));
                    queue.push(q);
                }
            }
        }
    }
    fn sift(levels: &[Level], from: usize, mut g: Vec<usize>) -> Vec<usize> {
        for l in &levels[from..] {
            match &l.transversal[g[l.base]] {
                None => return g,
                Some(u) => g = compose(&inverse(u), &g),
            }
        }
        g
    }
    fn add(levels: &mut Vec<Level>, i: usize, g: Vec<usize>, n: usize) {
        if i == levels.len() {
            let base = (0..n).find(|&k| g[k] != k).unwrap();
            levels.push(Level {
                base,
                gens: vec![],
                transversal: vec![],
            });
        }
        levels[i].gens.push(g);
        rebuild(&mut levels[i], n);
        'restart: loop {
            for p in 0..n {
                let Some(up) = levels[i].transversal[p].clone() else {
                    continue;
                };
                for gi in 0..levels[i].gens.len() {
                    let s = levels[i].gens[gi].clone();
                    let usp = levels[i].transversal[s[p]].clone().unwrap();
                    let schreier = compose(&inverse(&usp), &compose(&s, &up));
                    if !is_identity(&sift(levels, i + 1, schreier.clone())) {
                        add(levels, i + 1, schreier, n);
                        continue 'restart;
                    }
                }
            }
            break;
        }
    }
    let mut levels: Vec<Level> = vec![];
    for g in gens {
        if !is_identity(&sift(&levels, 0, g.clone())) {
            add(&mut levels, 0, g.clone(), n);
        }
    }
    levels
        .iter()
        .map(|l| l.transversal.iter().filter(|t| t.is_some()).count() as u128)
        .product()
}

/// Generators are automorphisms and generate a group of the reported order.
fn check_generators(g: &G, c: &Canonization<&u8, &u8, &u8>) {
    let (nh, nn) = (g.n_hedges(), g.node_bound());
    let mut points = vec![];
    for gen in c.generators() {
        for h in (0..nh).map(Hedge) {
            let image = gen.hedge(h);
            assert_eq!(
                colour(g, image),
                colour(g, h),
                "generator preserves colours"
            );
            assert_eq!(
                g.inv(image),
                gen.hedge(g.inv(h)),
                "generator commutes with inv"
            );
            assert_eq!(
                g.node_id(image),
                gen.node(g.node_id(h)),
                "generator moves hedges with nodes"
            );
        }
        for v in (0..nn).map(NodeIndex) {
            assert_eq!(g[gen.node(v)], g[v], "generator preserves node colours");
        }
        let hp = gen.hedge_permutation();
        let np = gen.node_permutation();
        points.push(
            hp.map()
                .iter()
                .copied()
                .chain(np.map().iter().map(|&v| v + nh))
                .collect(),
        );
    }
    assert_eq!(
        Some(schreier_sims_order(&points, nh + nn)),
        c.automorphism_count(),
        "Schreier-Sims order of the generators"
    );
}

/// Full check of one graph: |Aut| against brute force (and `expected`), generators, orbit
/// sizes by orbit-stabilizer, and invariance of the canonical form and relabelled graph under
/// random relabellings, with valid isomorphisms between the copies.
fn check(spec: &Spec, expected: Option<u128>, rng: &mut SmallRng) {
    let g = spec.build();
    let c = g.canonize();
    let brute = count_isomorphisms(&g, &g, None);
    assert_eq!(c.automorphism_count(), Some(brute), "{spec:?}");
    if let Some(expected) = expected {
        assert_eq!(brute, expected, "{spec:?}");
    }
    check_generators(&g, &c);
    let orbits = c.hedge_orbits();
    for h in (0..g.n_hedges()).map(Hedge).take(3) {
        let orbit = orbits.iter().filter(|(_, &r)| r == orbits[h]).count() as u128;
        let stabilizer = count_isomorphisms(&g, &g, Some((h, h)));
        assert_eq!(
            orbit * stabilizer,
            brute,
            "orbit-stabilizer at {h:?} for {spec:?}"
        );
    }
    let form = relabelled(&g, &c);
    for _ in 0..6 {
        let g2 = spec.relabelled_graph(rng);
        let c2 = g2.canonize();
        assert_eq!(c2.canonical_form(), c.canonical_form(), "{spec:?}");
        assert_eq!(c2.automorphism_count(), c.automorphism_count());
        assert_eq!(relabelled(&g2, &c2), form, "{spec:?}");
        check_isomorphism(&g, &g2);
    }
}

/// [`HedgeGraph::isomorphism_by`] finds an isomorphism `a -> b` that preserves colours,
/// pairing, incidence, edges and node keys.
fn check_isomorphism(a: &G, b: &G) {
    let iso = a
        .isomorphism_by(b, |v| v, |e| e, |h| h)
        .expect("isomorphic");
    for h in (0..a.n_hedges()).map(Hedge) {
        let image = iso.hedges[h];
        assert_eq!(colour(b, image), colour(a, h));
        assert_eq!(b.inv(image), iso.hedges[a.inv(h)]);
        assert_eq!(b.node_id(image), iso.nodes[a.node_id(h)]);
        assert_eq!(b[&image], iso.edges[a[&h]]);
    }
    for v in a.iter_node_ids() {
        assert_eq!(a[v], b[iso.nodes[v]], "isomorphisms preserve node keys");
    }
}

fn hand_cases() -> Vec<(&'static str, Spec, u128)> {
    use Orientation::{Default, Reversed};
    let s = Spec::new;
    vec![
        (
            "tadpole loop, equal keys, 2 legs",
            s(1).undirected(0, 0, 0, 0).leg(0, 1).leg(0, 2),
            2,
        ),
        (
            "tadpole loop, keys 1|2",
            s(1).undirected(0, 1, 0, 2).leg(0, 1).leg(0, 2),
            1,
        ),
        (
            "figure-eight vacuum",
            s(1).undirected(0, 0, 0, 0).undirected(0, 0, 0, 0),
            8,
        ),
        (
            "figure-eight, one loop keyed 1|2",
            s(1).undirected(0, 1, 0, 2).undirected(0, 0, 0, 0),
            2,
        ),
        (
            "figure-eight, both loops keyed 1|2",
            s(1).undirected(0, 1, 0, 2).undirected(0, 1, 0, 2),
            2,
        ),
        (
            "figure-eight, loops keyed 1|2 and 2|1",
            s(1).undirected(0, 1, 0, 2).undirected(0, 2, 0, 1),
            2,
        ),
        (
            "directed self-loop + 2 legs",
            s(1).edge(0, 0, 0, 0, 2, Default).leg(0, 1).leg(0, 2),
            1,
        ),
        (
            "two directed self-loops",
            s(1).edge(0, 0, 0, 0, 2, Default)
                .edge(0, 0, 0, 0, 2, Reversed),
            2,
        ),
        (
            "theta",
            s(2).undirected(0, 0, 1, 0)
                .undirected(0, 0, 1, 0)
                .undirected(0, 0, 1, 0),
            12,
        ),
        (
            "theta keys 1-1,2-2,3-3",
            s(2).undirected(0, 1, 1, 1)
                .undirected(0, 2, 1, 2)
                .undirected(0, 3, 1, 3),
            2,
        ),
        (
            "theta keys 1-2,2-3,3-1",
            s(2).undirected(0, 1, 1, 2)
                .undirected(0, 2, 1, 3)
                .undirected(0, 3, 1, 1),
            1,
        ),
        (
            "theta keys 1-2,2-1,3-3",
            s(2).undirected(0, 1, 1, 2)
                .undirected(0, 2, 1, 1)
                .undirected(0, 3, 1, 3),
            2,
        ),
        (
            "bubble 2-pt",
            s(2).undirected(0, 0, 1, 0)
                .undirected(0, 0, 1, 0)
                .leg(0, 1)
                .leg(1, 2),
            2,
        ),
        (
            "bubble 2-pt keys (1,2)/(1,2)",
            s(2).undirected(0, 1, 1, 1)
                .undirected(0, 2, 1, 2)
                .leg(0, 1)
                .leg(1, 2),
            1,
        ),
        (
            "bubble 2-pt, symmetrized legs",
            s(2).undirected(0, 0, 1, 0)
                .undirected(0, 0, 1, 0)
                .leg(0, 1)
                .leg(1, 1),
            4,
        ),
        (
            "fermion bubble",
            s(2).edge(0, 0, 1, 0, 2, Default)
                .edge(1, 0, 0, 0, 2, Default)
                .leg(0, 1)
                .leg(1, 2),
            1,
        ),
        (
            "fermion bubble, symmetrized legs",
            s(2).edge(0, 0, 1, 0, 2, Default)
                .edge(1, 0, 0, 0, 2, Default)
                .leg(0, 1)
                .leg(1, 1),
            2,
        ),
        (
            "melon-4 vacuum",
            s(2).undirected(0, 0, 1, 0)
                .undirected(0, 0, 1, 0)
                .undirected(0, 0, 1, 0)
                .undirected(0, 0, 1, 0),
            48,
        ),
        (
            "dumbbell",
            s(2).undirected(0, 0, 0, 0)
                .undirected(1, 0, 1, 0)
                .undirected(0, 0, 1, 0),
            8,
        ),
        (
            "dumbbell, one loop keyed",
            s(2).undirected(0, 1, 0, 2)
                .undirected(1, 0, 1, 0)
                .undirected(0, 0, 1, 0),
            2,
        ),
        (
            "K4",
            s(4).undirected(0, 0, 1, 0)
                .undirected(0, 0, 2, 0)
                .undirected(0, 0, 3, 0)
                .undirected(1, 0, 2, 0)
                .undirected(1, 0, 3, 0)
                .undirected(2, 0, 3, 0),
            24,
        ),
        (
            "4 identical legs on one node",
            s(1).leg(0, 5).leg(0, 5).leg(0, 5).leg(0, 5),
            24,
        ),
        ("empty graph", s(0), 1),
        ("three isolated nodes", s(3), 6),
        (
            "isolated nodes of two colours",
            Spec {
                nodes: vec![0, 1, 0],
                elements: vec![],
            },
            2,
        ),
        (
            "isolated node next to a tadpole",
            s(2).undirected(0, 0, 0, 0),
            2,
        ),
        (
            "directed 2-cycle",
            s(2).edge(0, 0, 1, 0, 0, Default)
                .edge(1, 0, 0, 0, 0, Default),
            2,
        ),
        (
            "antiparallel vs parallel",
            s(2).edge(0, 0, 1, 0, 0, Default)
                .edge(0, 0, 1, 0, 0, Reversed)
                .edge(0, 0, 1, 0, 0, Default),
            2,
        ),
        (
            "leg directions",
            s(1).edge(0, 0, 0, 0, 0, Orientation::Undirected)
                .leg(0, 0)
                .leg(0, 0),
            4,
        ),
    ]
}

#[test]
fn hand_cases_match_expected_brute_force_and_are_invariant() {
    let mut rng = SmallRng::seed_from_u64(42);
    for (label, spec, expected) in hand_cases() {
        eprintln!("{label}");
        check(&spec, Some(expected), &mut rng);
    }
}

#[test]
fn leg_directions_distinguish() {
    let out = Spec::new(1).leg(0, 0);
    let mut incoming = out.clone();
    let Element::Leg { flow, .. } = &mut incoming.elements[0] else {
        unreachable!()
    };
    *flow = Flow::Sink;
    let mut directed = out.clone();
    let Element::Leg { orientation, .. } = &mut directed.elements[0] else {
        unreachable!()
    };
    *orientation = Orientation::Default;
    let mut reversed = directed.clone();
    let Element::Leg {
        orientation, flow, ..
    } = &mut reversed.elements[0]
    else {
        unreachable!()
    };
    (*orientation, *flow) = (Orientation::Reversed, Flow::Sink);
    let form = |s: &Spec| s.build().canonical_form_by(ignore, ignore, ignore);
    // an undirected leg ignores its stored flow; a directed leg keeps its physical direction
    assert_eq!(form(&out), form(&incoming));
    assert_ne!(form(&out), form(&directed));
    assert_eq!(form(&directed), form(&reversed));
}

#[test]
fn canonical_form_records_isolated_nodes() {
    let form = |nodes: Vec<u8>| {
        Spec {
            nodes,
            elements: vec![],
        }
        .undirected(0, 0, 0, 0)
        .build()
        .canonize_by(|&v| v, ignore, ignore)
        .into_canonical_form()
    };
    assert_ne!(form(vec![0]), form(vec![0, 0]));
    assert_ne!(form(vec![0, 0, 1]), form(vec![0, 1, 1]));
    assert_eq!(form(vec![0, 1, 0]), form(vec![0, 0, 1]));
}

#[test]
fn canonical_form_parts_are_its_words_and_key_tables() {
    let g = Spec::new(2).undirected(0, 1, 1, 2).leg(0, 3).build();
    let form = g.canonical_form();
    let (words, nodes, edges, hedges) = form.clone().into_parts();
    assert_eq!(words, form.words());
    assert_eq!(nodes, form.node_keys());
    assert_eq!(edges, form.edge_keys());
    assert_eq!(hedges, form.hedge_keys());
}

#[test]
fn ignored_payloads_do_not_distinguish() {
    let a = Spec::new(2)
        .edge(0, 1, 1, 2, 3, Orientation::Undirected)
        .build();
    let b = Spec::new(2)
        .edge(0, 7, 1, 7, 5, Orientation::Undirected)
        .build();
    assert!(a.is_isomorphic_by(&b, ignore, ignore, ignore));
    assert!(!a.is_isomorphic_by(&b, ignore, |&e| e, ignore));
    assert!(!a.is_isomorphic_by(&b, ignore, ignore, |&h| h));
    let keyed = a.canonize_by(ignore, ignore, |&h| h);
    let unkeyed = a.canonize_by(ignore, ignore, ignore);
    assert_eq!(keyed.automorphism_count(), Some(1));
    assert_eq!(unkeyed.automorphism_count(), Some(2));
}

#[test]
fn symmetric_graphs() {
    let mut rng = SmallRng::seed_from_u64(7);
    let from_edges = |n: usize, edges: &[(usize, usize)]| {
        edges
            .iter()
            .fold(Spec::new(n), |s, &(a, b)| s.undirected(a, 0, b, 0))
    };
    let cycle: Vec<_> = (0..8).map(|i| (i, (i + 1) % 8)).collect();
    let petersen: Vec<_> = (0..5)
        .flat_map(|i| [(i, (i + 1) % 5), (i, i + 5), (i + 5, (i + 2) % 5 + 5)])
        .collect();
    let cube: Vec<_> = (0..8usize)
        .flat_map(|i| (0..3).map(move |b| (i, i ^ (1 << b))))
        .filter(|(a, b)| a < b)
        .collect();
    let k5: Vec<_> = (0..5)
        .flat_map(|a| (a + 1..5).map(move |b| (a, b)))
        .collect();
    let k33: Vec<_> = (0..3).flat_map(|a| (3..6).map(move |b| (a, b))).collect();
    let disjoint_triangles: Vec<_> = (0..3)
        .flat_map(|t| (0..3).map(move |i| (3 * t + i, 3 * t + (i + 1) % 3)))
        .collect();
    for (spec, order) in [
        (from_edges(8, &cycle), 16),
        (from_edges(10, &petersen), 120),
        (from_edges(8, &cube), 48),
        (from_edges(5, &k5), 120),
        (from_edges(6, &k33), 72),
        (from_edges(9, &disjoint_triangles), 6 * 6 * 6 * 6),
    ] {
        let g = spec.build();
        let c = g.canonize();
        assert_eq!(c.automorphism_count(), Some(order), "{spec:?}");
        check_generators(&g, &c);
        let form = relabelled(&g, &c);
        for _ in 0..4 {
            let g2 = spec.relabelled_graph(&mut rng);
            let c2 = g2.canonize();
            assert_eq!(relabelled(&g2, &c2), form);
            assert_eq!(c2.canonical_form(), c.canonical_form());
        }
    }
}

#[test]
fn random_coloured_multigraphs_match_brute_force() {
    let mut rng = SmallRng::seed_from_u64(0x5EED);
    for _ in 0..600 {
        let spec = Spec::random(&mut rng, 5, 6, 3);
        check(&spec, None, &mut rng);
    }
}

#[test]
fn canonical_forms_separate_exactly_the_isomorphism_classes() {
    let mut rng = SmallRng::seed_from_u64(0xC1A55);
    // few nodes and many edges make non-isomorphic graphs with equal local data common
    let specs: Vec<Spec> = (0..400)
        .map(|i| Spec::random(&mut rng, 2 + i % 2, 3 + i % 2, 1))
        .collect();
    let graphs: Vec<G> = specs.iter().map(Spec::build).collect();
    let forms: Vec<_> = graphs
        .iter()
        .map(|g| {
            // the form alone is the form of the full canonization
            let form = g.canonical_form();
            assert_eq!(form, *g.canonize().canonical_form());
            form
        })
        .collect();
    for i in 0..graphs.len() {
        for j in 0..i {
            let isomorphic = count_isomorphisms(&graphs[i], &graphs[j], None) > 0;
            assert_eq!(
                forms[i] == forms[j],
                isomorphic,
                "{:?} vs {:?}",
                specs[i],
                specs[j]
            );
        }
    }
}

#[test]
fn many_identical_legs_overflow_to_none() {
    let spec = (0..40).fold(Spec::new(1), |s, _| s.leg(0, 0));
    let g = spec.build();
    let c = g.canonize();
    assert_eq!(c.automorphism_count(), None);
    assert_eq!(c.generators().len(), 39);
}

#[test]
fn parallel_edge_end_pairing_matters() {
    let straight = Spec::new(2).undirected(0, 1, 1, 1).undirected(0, 2, 1, 2);
    let crossed = Spec::new(2).undirected(0, 1, 1, 2).undirected(0, 2, 1, 1);
    let (a, b) = (straight.build(), crossed.build());
    assert!(!a.is_isomorphic_by(&b, ignore, ignore, |&h| h));
    assert!(a.is_isomorphic_by(&b, ignore, ignore, ignore));
}

/// Identifies random groups of nodes without compacting. Requests may repeat nodes or name
/// absorbed ones, and survivors get random keys, so absorbed slots may keep keys that no live
/// node carries.
fn identify_some(g: &mut G, rng: &mut SmallRng) {
    for _ in 0..rng.gen_range(1..=3) {
        let nodes: Vec<NodeIndex> = (0..rng.gen_range(1..=3))
            .map(|_| NodeIndex(rng.gen_range(0..g.node_bound())))
            .collect();
        g.identify_nodes(&nodes, rng.gen_range(0..3));
    }
}

/// A graph with pending identifications canonizes exactly like its compacted copy, through the
/// remap of [`HedgeGraph::compact`]: node lookups resolve absorbed nodes to their
/// representative, the relabelled graphs agree, the shared generators pass
/// [`check_generators`] on the compacted copy, and isomorphisms to relabelled copies are valid
/// both ways.
fn check_pending(g: &G, rng: &mut SmallRng) {
    let mut compacted = g.clone();
    let remap = compacted.compact();
    let to = |v: NodeIndex| remap.get(v).expect("compaction keeps or merges every node");
    let slots = || (0..g.node_bound()).map(NodeIndex);
    let (c, d) = (g.canonize(), compacted.canonize());
    assert_eq!(c.canonical_form(), d.canonical_form());
    assert_eq!(c.automorphism_count(), d.automorphism_count());
    assert_eq!(c.automorphism_count(), Some(count_isomorphisms(g, g, None)));
    check_generators(&compacted, &d);
    assert_eq!(relabelled(g, &c), relabelled(&compacted, &d));
    assert_eq!(c.hedge_permutation(), d.hedge_permutation());
    assert_eq!(c.edge_permutation(), d.edge_permutation());
    assert_eq!(c.hedge_orbits(), d.hedge_orbits());
    assert_eq!(c.edge_orbits(), d.edge_orbits());
    for v in slots() {
        assert_eq!(c.canonical_node(v), d.canonical_node(to(v)));
    }
    assert_eq!(c.generators().len(), d.generators().len());
    for (a, b) in c.generators().zip(d.generators()) {
        assert_eq!(a.hedge_permutation(), b.hedge_permutation());
        for v in slots() {
            assert!(g.node_store.is_live(a.node(v)));
            assert_eq!(to(a.node(v)), b.node(to(v)));
        }
    }
    let (orbits, compacted_orbits) = (c.node_orbits(), d.node_orbits());
    for u in slots() {
        assert!(
            g.node_store.is_live(orbits[u]),
            "live nodes name the orbits"
        );
        for v in slots() {
            assert_eq!(
                orbits[u] == orbits[v],
                compacted_orbits[to(u)] == compacted_orbits[to(v)]
            );
        }
    }
    assert!(g.is_isomorphic_by(&compacted, |v| v, |e| e, |h| h));
    let iso = g
        .isomorphism_by(&compacted, |v| v, |e| e, |h| h)
        .expect("isomorphic");
    assert!(slots().all(|v| iso.nodes[v] == to(v)));
    assert!(iso.hedges.iter().all(|(h, &image)| image == h));
    let back = compacted
        .isomorphism_by(g, |v| v, |e| e, |h| h)
        .expect("isomorphic");
    assert!(g.iter_node_ids().all(|v| back.nodes[to(v)] == v));

    let shuffled = |n: usize, rng: &mut SmallRng| {
        let mut map: Vec<usize> = (0..n).collect();
        map.shuffle(rng);
        Permutation::from_map(map)
    };
    let mut relabelled = compacted;
    let (hedges, nodes) = (
        shuffled(relabelled.n_hedges(), rng),
        shuffled(relabelled.node_bound(), rng),
    );
    relabelled
        .permute_hedges(&HedgeRemap::from_permutation(&hedges))
        .unwrap();
    relabelled
        .permute_nodes(&NodeRemap::from_permutation(&nodes))
        .unwrap();
    check_isomorphism(g, &relabelled);
    check_isomorphism(&relabelled, g);
}

#[test]
fn absorbed_nodes_are_not_isolated_nodes() {
    // a square whose opposite corners merge: two nodes joined by a double edge, each with a
    // self-loop, and two absorbed slots
    let mut g = Spec::new(4)
        .undirected(0, 0, 1, 0)
        .undirected(1, 0, 2, 0)
        .undirected(2, 0, 3, 0)
        .undirected(3, 0, 0, 0)
        .build();
    let survivors = [
        g.identify_nodes(&[NodeIndex(1), NodeIndex(2)], 0),
        g.identify_nodes(&[NodeIndex(3), NodeIndex(0)], 0),
    ];
    assert_eq!((g.n_live(), g.node_bound()), (2, 4));
    let c = g.canonize();
    // swap the nodes, the parallel edges and the ends of either self-loop
    assert_eq!(c.automorphism_count(), Some(16));
    for (absorbed, survivor) in [NodeIndex(2), NodeIndex(0)].into_iter().zip(survivors) {
        assert_eq!(c.canonical_node(absorbed), c.canonical_node(survivor));
    }
    check_pending(&g, &mut SmallRng::seed_from_u64(1));
}

#[test]
fn pending_identifications_canonize_like_the_compacted_graph() {
    let mut rng = SmallRng::seed_from_u64(0xAB50);
    let mut specs: Vec<Spec> = hand_cases().into_iter().map(|(_, spec, _)| spec).collect();
    specs.extend((0..300).map(|_| Spec::random(&mut rng, 5, 6, 3)));
    let mut pending = 0;
    for spec in specs.iter().filter(|spec| !spec.nodes.is_empty()) {
        let mut g = spec.build();
        identify_some(&mut g, &mut rng);
        pending += usize::from(g.n_live() < g.node_bound());
        check_pending(&g, &mut rng);
    }
    assert!(2 * pending > specs.len(), "most graphs keep absorbed nodes");
}

/// Two nodes joined by an edge, identified into one with a self-loop.
fn merged_edge() -> G {
    let mut g = Spec::new(2).undirected(0, 0, 1, 0).build();
    g.identify_nodes(&[NodeIndex(0), NodeIndex(1)], 0);
    g
}

#[test]
#[should_panic(expected = "compact")]
fn node_permutations_need_compacted_graphs() {
    merged_edge().canonize().node_permutation();
}

#[test]
#[should_panic(expected = "compact")]
fn automorphism_node_permutations_need_compacted_graphs() {
    let g = merged_edge();
    let c = g.canonize();
    let flip = c.generators().next().expect("the self-loop flips");
    flip.node_permutation();
}

/// Oracle: Graphica on generated Feynman-graph sets, plain graphs and the exact incidence
/// encoding for half-edge keys.
#[cfg(feature = "symbolica")]
mod graphica {
    use rand::{rngs::SmallRng, seq::SliceRandom, Rng, SeedableRng};
    use symbolica::graph::{GenerationSettings, Graph, HalfEdge};

    use crate::half_edge::{
        builder::HedgeGraphBuilder,
        canon::ignore,
        involution::{Flow, HedgePair, Orientation},
        HedgeGraph,
    };

    /// Node data: 0 for vertices, the leg label for external nodes.
    type Plain = Graph<usize, u8>;
    type L = HedgeGraph<u8, usize, usize>;

    fn aut(g: &Plain) -> u128 {
        g.canonize()
            .automorphism_group_size
            .to_string()
            .parse()
            .unwrap()
    }

    fn generate(
        externals: &[(usize, HalfEdge<u8>)],
        vertices: &[Vec<HalfEdge<u8>>],
        loops: usize,
        one_pi: bool,
    ) -> Vec<Plain> {
        let min_degree = vertices.iter().map(Vec::len).min().unwrap();
        let max_vertices = (2 * loops + externals.len() - 2) / (min_degree - 2);
        let mut settings = GenerationSettings::new()
            .max_loops(loops)
            .max_vertices(externals.len() + max_vertices)
            .allow_self_loops(true)
            .allow_zero_flow_edges(!one_pi);
        if one_pi {
            settings = settings.max_bridges(0);
        }
        let mut graphs: Vec<Plain> = Plain::generate(externals, vertices, settings)
            .unwrap()
            .into_keys()
            .filter(|g| g.num_loops() == loops)
            .collect();
        graphs.sort();
        graphs
    }

    fn sets() -> Vec<(&'static str, Vec<Plain>)> {
        let phi = HalfEdge::undirected(0);
        let (photon, electron) = (HalfEdge::undirected(1), HalfEdge::incoming(2));
        let (gluon, quark, ghost) = (
            HalfEdge::undirected(3),
            HalfEdge::incoming(4),
            HalfEdge::incoming(7),
        );
        let qcd = vec![
            vec![gluon; 3],
            vec![gluon; 4],
            vec![quark, quark.flip(), gluon],
            vec![ghost, ghost.flip(), gluon],
        ];
        vec![
            ("phi3 vacuum L4", generate(&[], &[vec![phi; 3]], 4, false)),
            ("phi3 vacuum L5", generate(&[], &[vec![phi; 3]], 5, false)),
            ("phi4 vacuum L4", generate(&[], &[vec![phi; 4]], 4, false)),
            (
                "phi34 2-point L3",
                generate(
                    &[(1, phi), (2, phi)],
                    &[vec![phi; 3], vec![phi; 4]],
                    3,
                    false,
                ),
            ),
            (
                "QED photon self-energy L3",
                generate(
                    &[(1, photon), (2, photon)],
                    &[vec![electron, electron.flip(), photon]],
                    3,
                    true,
                ),
            ),
            (
                "QCD gluon self-energy L2",
                generate(&[(1, gluon), (2, gluon)], &qcd, 2, true),
            ),
            (
                "gg -> gg L1",
                generate(
                    &[(1, gluon), (2, gluon), (3, gluon), (4, gluon)],
                    &qcd,
                    1,
                    false,
                ),
            ),
        ]
    }

    /// Random node order, edge order and endpoint order of undirected edges.
    fn permute(g: &Plain, rng: &mut SmallRng) -> Plain {
        let mut perm: Vec<usize> = (0..g.nodes().len()).collect();
        perm.shuffle(rng);
        let mut nodes = vec![0; perm.len()];
        for (v, &p) in perm.iter().enumerate() {
            nodes[p] = g.nodes()[v].data;
        }
        let mut h = Plain::new();
        for d in nodes {
            h.add_node(d);
        }
        let mut edges: Vec<_> = g.edges().iter().collect();
        edges.shuffle(rng);
        for e in edges {
            let (mut a, mut b) = (perm[e.vertices.0], perm[e.vertices.1]);
            if !e.directed && rng.gen_bool(0.5) {
                std::mem::swap(&mut a, &mut b);
            }
            h.add_edge(a, b, e.directed, e.data).unwrap();
        }
        h
    }

    /// External nodes become dangling hedges keyed by their label.
    fn to_linnet(g: &Plain) -> L {
        let mut builder = HedgeGraphBuilder::new();
        let nodes: Vec<_> = g
            .nodes()
            .iter()
            .map(|n| (n.data == 0).then(|| builder.add_node(0)))
            .collect();
        let orientation = |directed| {
            if directed {
                Orientation::Default
            } else {
                Orientation::Undirected
            }
        };
        for e in g.edges() {
            let (x, y) = e.vertices;
            let o = orientation(e.directed);
            match (nodes[x], nodes[y]) {
                (Some(a), Some(b)) => builder.add_edge(a.add_data(0), b.add_data(0), e.data, o),
                (Some(a), None) => builder.add_external_edge(
                    a.add_data(g.nodes()[y].data),
                    e.data,
                    o,
                    Flow::Source,
                ),
                (None, Some(b)) => {
                    builder.add_external_edge(b.add_data(g.nodes()[x].data), e.data, o, Flow::Sink)
                }
                (None, None) => unreachable!("legs join vertices"),
            }
        }
        builder.build()
    }

    /// Exact incidence encoding: a node per vertex, a centre per edge and a port per hedge
    /// coloured by its role and key.
    fn incidence(g: &L) -> Graph<(u8, usize), ()> {
        let mut out = Graph::new();
        let nodes: Vec<usize> = g
            .iter_nodes()
            .map(|(_, _, &c)| out.add_node((0, c)))
            .collect();
        for (pair, _, data) in g.iter_edges() {
            let centre = out.add_node((1, *data.data as usize));
            let hedges = match pair {
                HedgePair::Paired { source, sink } => vec![source, sink],
                HedgePair::Unpaired { hedge, .. } => vec![hedge],
                HedgePair::Split { source, sink, .. } => vec![source, sink],
            };
            for h in hedges {
                let role = match g.superficial_hedge_orientation(h) {
                    None => 2,
                    Some(Flow::Source) => 3,
                    Some(Flow::Sink) => 4,
                };
                let port = out.add_node((role, g[h]));
                out.add_edge(centre, port, false, ()).unwrap();
                out.add_edge(port, nodes[g.node_id(h).0], false, ())
                    .unwrap();
            }
        }
        out
    }

    #[test]
    fn automorphism_counts_and_classes_match_graphica() {
        let mut rng = SmallRng::seed_from_u64(0x6A9);
        for (name, graphs) in sets() {
            assert!(!graphs.is_empty(), "{name}");
            let mut forms = std::collections::HashSet::new();
            for g in &graphs {
                let l = to_linnet(g);
                let c = l.canonize_by(|&v| v, |&e| e, |&h| h);
                assert_eq!(c.automorphism_count(), Some(aut(g)), "{name}");
                // symmetrized legs: equal labels
                let mut sym = Plain::new();
                for n in g.nodes() {
                    sym.add_node((n.data != 0) as usize);
                }
                for e in g.edges() {
                    sym.add_edge(e.vertices.0, e.vertices.1, e.directed, e.data)
                        .unwrap();
                }
                let sym_count = to_linnet(&sym)
                    .canonize_by(ignore, |&e| e, |&h| h)
                    .automorphism_count();
                assert_eq!(sym_count, Some(aut(&sym)), "{name} symmetrized");
                let form = c.into_canonical_form();
                for _ in 0..2 {
                    let p = to_linnet(&permute(g, &mut rng));
                    assert_eq!(p.canonical_form_by(|&v| v, |&e| e, |&h| h), form, "{name}");
                }
                forms.insert(form);
            }
            assert_eq!(
                forms.len(),
                graphs.len(),
                "{name}: distinct graphs, distinct canonical forms"
            );
        }
    }

    #[test]
    fn hedge_keyed_variants_match_the_incidence_encoding() {
        let mut rng = SmallRng::seed_from_u64(0xBADC0DE);
        for (name, graphs) in sets() {
            for g in graphs.iter().step_by(graphs.len().div_ceil(40)) {
                let base = to_linnet(g);
                for _ in 0..3 {
                    let keyed = base.map_data_ref(
                        |_, _, &v| v,
                        |_, _, _, d| d.map(|&e| e),
                        |h, &k| {
                            if base.is_dangling(h) {
                                k
                            } else {
                                rng.gen_range(0..2)
                            }
                        },
                    );
                    let count = keyed
                        .canonize_by(|&v| v, |&e| e, |&h| h)
                        .automorphism_count();
                    let oracle: u128 = incidence(&keyed)
                        .canonize()
                        .automorphism_group_size
                        .to_string()
                        .parse()
                        .unwrap();
                    assert_eq!(count, Some(oracle), "{name}");
                }
            }
        }
    }
}

#[test]
fn many_distinct_keys_use_the_dense_colour_table() {
    // 20k parallel edges with distinct edge and half-edge keys overflow one-word colours
    let n = 20_000u32;
    let build = |swap: bool| -> HedgeGraph<u32, u8, u32> {
        let mut builder = HedgeGraphBuilder::new();
        let (a, b) = (builder.add_node(0u8), builder.add_node(0u8));
        for e in 0..n {
            let (ka, kb) = if swap && e == 7 {
                (2 * e + 1, 2 * e)
            } else {
                (2 * e, 2 * e + 1)
            };
            builder.add_edge(a.add_data(ka), b.add_data(kb), e, Orientation::Undirected);
        }
        builder.build()
    };
    let (g, swapped) = (build(false), build(true));
    let mut flat = super::flat::Flat::default();
    flat.load(&g, |&v| v, |&e| e, |&h| h);
    assert!(!flat.colour_table.is_empty());
    let c = g.canonize();
    assert_eq!(c.automorphism_count(), Some(1));
    let mut map: Vec<usize> = (0..g.n_hedges()).collect();
    map.shuffle(&mut SmallRng::seed_from_u64(3));
    let mut relabelled = build(false);
    relabelled
        .permute_hedges(&HedgeRemap::from_permutation(&Permutation::from_map(map)))
        .unwrap();
    assert_eq!(relabelled.canonical_form(), *c.canonical_form());
    assert_ne!(swapped.canonical_form(), *c.canonical_form());
}
