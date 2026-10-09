//! Graphs as nauty sees them, their half-edge form, and the checks every translated run shares:
//! valid generators of a group of the reported order, and certificates that do not depend on
//! the labelling.
//!
//! nauty's simple graphs become undirected paired edges, digraphs directed edges `from -> to`,
//! multigraphs parallel edges and vertex colourings (ordered partitions) node keys. None of the
//! suite's inputs has loops or digons, so `|Aut|` and the orbits agree exactly with nauty's except
//! for multigraphs, where Linnet also counts the swaps of parallel edges.

use std::collections::BTreeMap;

use linnet::half_edge::{
    builder::HedgeGraphBuilder,
    canon::{ignore, Canonization, Certificate},
    involution::{Hedge, Orientation},
    HedgeGraph, NoData, NodeIndex,
};
use rand::{rngs::SmallRng, seq::SliceRandom, Rng};

pub type Linnet = HedgeGraph<(), u32>;
pub type Canon = Canonization<u32, NoData, NoData>;
pub type Cert = Certificate<u32, NoData, NoData>;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Graph {
    pub n: usize,
    pub directed: bool,
    /// Edges, or arcs `from -> to`; parallel edges are repeated.
    pub edges: Vec<(usize, usize)>,
    /// Node keys: a vertex colouring, or the cell of each vertex in an ordered partition.
    pub colours: Vec<u32>,
}

impl Graph {
    pub fn undirected(n: usize, edges: Vec<(usize, usize)>) -> Self {
        Graph {
            n,
            directed: false,
            edges,
            colours: vec![0; n],
        }
    }

    pub fn directed(n: usize, arcs: Vec<(usize, usize)>) -> Self {
        Graph {
            directed: true,
            ..Graph::undirected(n, arcs)
        }
    }

    /// The undirected graph of a generated half-edge graph without legs, payloads dropped.
    pub fn of_hedge_graph<E, V, H>(g: &HedgeGraph<E, V, H>) -> Self {
        let edges = (0..g.n_hedges())
            .map(Hedge)
            .filter(|&h| g.inv(h) > h)
            .map(|h| (g.node_id(h).0, g.node_id(g.inv(h)).0))
            .collect();
        Graph::undirected(g.node_bound(), edges)
    }

    pub fn with_colours(mut self, colours: Vec<u32>) -> Self {
        assert_eq!(colours.len(), self.n);
        self.colours = colours;
        self
    }

    pub fn hedge_graph(&self) -> Linnet {
        let mut builder = HedgeGraphBuilder::new();
        let nodes: Vec<NodeIndex> = self.colours.iter().map(|&c| builder.add_node(c)).collect();
        let orientation = if self.directed {
            Orientation::Default
        } else {
            Orientation::Undirected
        };
        for &(a, b) in &self.edges {
            builder.add_edge(nodes[a], nodes[b], (), orientation);
        }
        builder.build()
    }

    pub fn canonize(&self) -> Canon {
        self.hedge_graph().canonize_by(|&c| c, ignore, ignore)
    }

    pub fn certificate(&self) -> Cert {
        self.canonize().into_certificate()
    }

    /// The edges as a sorted multiset, ends sorted for undirected graphs.
    fn edge_multiset(&self, map: impl Fn(usize) -> usize) -> Vec<(usize, usize)> {
        let mut edges: Vec<_> = self
            .edges
            .iter()
            .map(|&(a, b)| {
                let (a, b) = (map(a), map(b));
                if self.directed || a <= b {
                    (a, b)
                } else {
                    (b, a)
                }
            })
            .collect();
        edges.sort_unstable();
        edges
    }

    /// Whether the vertex map `image` preserves the colours and the edges.
    pub fn is_automorphism(&self, image: &[usize]) -> bool {
        image.len() == self.n
            && (0..self.n).all(|v| self.colours[image[v]] == self.colours[v])
            && self.edge_multiset(|v| image[v]) == self.edge_multiset(|v| v)
    }

    /// nauty's `relabel(g, perm)`: vertex `i` of the result is vertex `perm[i]` of `self`.
    pub fn relabelled_by(&self, perm: &[usize]) -> Graph {
        let mut position = vec![0; self.n];
        for (i, &p) in perm.iter().enumerate() {
            position[p] = i;
        }
        let edges = self
            .edges
            .iter()
            .map(|&(a, b)| (position[a], position[b]))
            .collect();
        let colours = perm.iter().map(|&p| self.colours[p]).collect();
        Graph {
            edges,
            colours,
            ..self.clone()
        }
    }

    /// A random relabelling, with the edges listed in random order and undirected edges stored
    /// with random ends first.
    pub fn shuffled(&self, rng: &mut SmallRng) -> Graph {
        let mut perm: Vec<usize> = (0..self.n).collect();
        perm.shuffle(rng);
        let mut g = self.relabelled_by(&perm);
        g.edges.shuffle(rng);
        if !self.directed {
            for e in &mut g.edges {
                if rng.gen_bool(0.5) {
                    *e = (e.1, e.0);
                }
            }
        }
        g
    }

    /// The complement of a loopless simple graph.
    pub fn complement(&self) -> Graph {
        let edges = self.edge_multiset(|v| v);
        let pairs: Vec<_> = if self.directed {
            (0..self.n)
                .flat_map(|a| (0..self.n).map(move |b| (a, b)))
                .filter(|(a, b)| a != b)
                .collect()
        } else {
            (0..self.n)
                .flat_map(|a| (a + 1..self.n).map(move |b| (a, b)))
                .collect()
        };
        let edges = pairs
            .into_iter()
            .filter(|e| edges.binary_search(e).is_err())
            .collect();
        Graph {
            edges,
            ..self.clone()
        }
    }

    pub fn converse(&self) -> Graph {
        Graph {
            edges: self.edges.iter().map(|&(a, b)| (b, a)).collect(),
            ..self.clone()
        }
    }

    /// Neighbours in the underlying undirected graph.
    pub fn neighbours(&self) -> Vec<Vec<usize>> {
        let mut adjacency = vec![vec![]; self.n];
        for &(a, b) in &self.edges {
            adjacency[a].push(b);
            adjacency[b].push(a);
        }
        adjacency
    }

    /// Distances in the underlying undirected graph (`usize::MAX` when disconnected).
    pub fn distances(&self) -> Vec<Vec<usize>> {
        let adjacency = self.neighbours();
        (0..self.n)
            .map(|s| {
                let mut dist = vec![usize::MAX; self.n];
                dist[s] = 0;
                let mut queue = std::collections::VecDeque::from([s]);
                while let Some(v) = queue.pop_front() {
                    for &w in &adjacency[v] {
                        if dist[w] == usize::MAX {
                            dist[w] = dist[v] + 1;
                            queue.push_back(w);
                        }
                    }
                }
                dist
            })
            .collect()
    }

    /// The vertex sets of the connected components, by smallest vertex.
    pub fn components(&self) -> Vec<Vec<usize>> {
        let distances = self.distances();
        partition(
            distances
                .iter()
                .map(|d| d.iter().position(|&x| x != usize::MAX).unwrap()),
        )
    }

    /// Radius and diameter of a connected graph.
    pub fn radius_diameter(&self) -> (usize, usize) {
        let eccentricities: Vec<usize> = self
            .distances()
            .iter()
            .map(|d| *d.iter().max().unwrap())
            .collect();
        (
            *eccentricities.iter().min().unwrap(),
            *eccentricities.iter().max().unwrap(),
        )
    }

    /// Linnet's automorphism group, after checking that every generator is an automorphism of
    /// the half-edge graph and of `self` and that the generators generate a group of the
    /// reported order.
    pub fn symmetry(&self) -> Symmetry {
        let g = self.hedge_graph();
        let c = g.canonize_by(|&c| c, ignore, ignore);
        let order = c.automorphism_count().expect("order fits in u128");
        let (n_hedges, n) = (g.n_hedges(), self.n);
        // without loops and parallel edges the group acts faithfully on the nodes alone
        let edges = self.edge_multiset(|v| v);
        let simple = edges.windows(2).all(|w| w[0] != w[1]) && edges.iter().all(|(a, b)| a != b);
        let mut points = vec![];
        for gen in c.generators() {
            for h in (0..n_hedges).map(Hedge) {
                let image = gen.hedge(h);
                assert_eq!(g.inv(image), gen.hedge(g.inv(h)), "commutes with pairing");
                assert_eq!(
                    g.node_id(image),
                    gen.node(g.node_id(h)),
                    "moves hedges with nodes"
                );
                assert_eq!(
                    g.superficial_hedge_orientation(image),
                    g.superficial_hedge_orientation(h),
                    "keeps directions"
                );
            }
            let nodes: Vec<usize> = (0..n).map(|v| gen.node(NodeIndex(v)).0).collect();
            assert!(
                self.is_automorphism(&nodes),
                "generator {nodes:?} of {self:?}"
            );
            points.push(if simple {
                nodes
            } else {
                (0..n_hedges)
                    .map(|h| gen.hedge(Hedge(h)).0)
                    .chain(nodes.iter().map(|&v| v + n_hedges))
                    .collect()
            });
        }
        assert_eq!(
            group_order(&points, if simple { n } else { n_hedges + n }),
            order,
            "Schreier-Sims order of the generators of {self:?}"
        );
        let node_orbits = c.node_orbits();
        let edge_orbits = c.edge_orbits();
        let mut edge_roots: Vec<_> = edge_orbits.iter().map(|(_, r)| *r).collect();
        edge_roots.sort_unstable();
        edge_roots.dedup();
        Symmetry {
            order,
            node_orbits: partition((0..n).map(|v| node_orbits[NodeIndex(v)].0)),
            edge_orbits: edge_roots.len(),
            certificate: c.into_certificate(),
        }
    }

    /// The certificate is unchanged by `copies` random relabellings, and the isomorphisms Linnet
    /// returns map the edges of `self` onto those of the copy.
    pub fn check_invariance(&self, certificate: &Cert, copies: usize, rng: &mut SmallRng) {
        let g = self.hedge_graph();
        for _ in 0..copies {
            let copy = self.shuffled(rng);
            let h = copy.hedge_graph();
            assert_eq!(&copy.certificate(), certificate, "relabelled {self:?}");
            let iso = g
                .isomorphism_by(&h, |&c| c, ignore, ignore)
                .expect("isomorphic");
            let image: Vec<usize> = (0..self.n).map(|v| iso.nodes[NodeIndex(v)].0).collect();
            assert_eq!(
                self.edge_multiset(|v| image[v]),
                copy.edge_multiset(|v| v),
                "isomorphism of {self:?}"
            );
            assert!((0..self.n).all(|v| self.colours[v] == copy.colours[image[v]]));
        }
    }
}

/// Linnet's automorphism group of a [`Graph`].
#[derive(Clone, Debug)]
pub struct Symmetry {
    pub order: u128,
    /// Node orbits, each sorted, ordered by smallest element.
    pub node_orbits: Vec<Vec<usize>>,
    pub edge_orbits: usize,
    pub certificate: Cert,
}

/// The classes of `0..n` under a map to class representatives, each sorted, ordered by
/// smallest element.
pub fn partition(roots: impl Iterator<Item = usize>) -> Vec<Vec<usize>> {
    let mut classes: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (v, root) in roots.enumerate() {
        classes.entry(root).or_default().push(v);
    }
    let mut classes: Vec<_> = classes.into_values().collect();
    classes.sort_unstable();
    classes
}

/// Order of the permutation group generated by `generators` on `0..n`, by deterministic
/// Schreier-Sims (Holt, Handbook of Computational Group Theory, SCHREIERSIMS): every Schreier
/// generator is sifted once, and a non-trivial residue becomes a strong generator. Linnet's own
/// helper is private to its unit tests.
pub fn group_order(generators: &[Vec<usize>], n: usize) -> u128 {
    /// `u` and `u^-1` with `u(base) = point`, for the points of the basic orbit.
    type Transversal = Vec<Option<(Vec<usize>, Vec<usize>)>>;
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
    let is_identity = |a: &[usize]| a.iter().enumerate().all(|(i, &x)| i == x);
    let mut base: Vec<usize> = vec![];
    // strong generators with their depth: each fixes the base points before its depth and lies
    // in the stabilisers of all levels up to it
    let mut strong: Vec<(Vec<usize>, usize)> = vec![];
    let mut transversals: Vec<Transversal> = vec![];
    let mut tested: Vec<std::collections::HashSet<(usize, usize)>> = vec![];
    let new_level = |point: usize| {
        let mut t: Transversal = vec![None; n];
        let id: Vec<usize> = (0..n).collect();
        t[point] = Some((id.clone(), id));
        t
    };
    // extend a basic orbit, keeping the transversal elements already chosen
    let extend = |t: &mut Transversal, strong: &[(Vec<usize>, usize)], level: usize| {
        let mut queue: Vec<usize> = (0..n).filter(|&p| t[p].is_some()).collect();
        while let Some(p) = queue.pop() {
            for (g, _) in strong.iter().filter(|(_, d)| *d >= level) {
                let q = g[p];
                if t[q].is_none() {
                    let u = compose(g, &t[p].as_ref().unwrap().0);
                    let inv = inverse(&u);
                    t[q] = Some((u, inv));
                    queue.push(q);
                }
            }
        }
    };
    for g in generators.iter().filter(|g| !is_identity(g)) {
        if base.iter().all(|&b| g[b] == b) {
            let point = (0..n).find(|&k| g[k] != k).unwrap();
            base.push(point);
            transversals.push(new_level(point));
            tested.push(Default::default());
        }
        strong.push((g.clone(), 0));
    }
    for (level, t) in transversals.iter_mut().enumerate() {
        extend(t, &strong, level);
    }
    let mut level = base.len();
    while level > 0 {
        let i = level - 1;
        let mut descend = None;
        'pairs: for p in 0..n {
            let Some((up, _)) = transversals[i][p].clone() else {
                continue;
            };
            for gi in 0..strong.len() {
                if strong[gi].1 < i || !tested[i].insert((p, gi)) {
                    continue;
                }
                let s = &strong[gi].0;
                let inv_usp = &transversals[i][s[p]].as_ref().unwrap().1;
                let mut residue = compose(inv_usp, &compose(s, &up));
                // sift through the deeper levels
                let mut j = i + 1;
                while j < base.len() {
                    match &transversals[j][residue[base[j]]] {
                        Some((_, inv)) => residue = compose(inv, &residue),
                        None => break,
                    }
                    j += 1;
                }
                if !is_identity(&residue) {
                    if j == base.len() {
                        let point = (0..n).find(|&k| residue[k] != k).unwrap();
                        base.push(point);
                        transversals.push(new_level(point));
                        tested.push(Default::default());
                    }
                    strong.push((residue, j));
                    for (l, t) in transversals.iter_mut().enumerate().take(j + 1).skip(i + 1) {
                        extend(t, &strong, l);
                    }
                    descend = Some(j + 1);
                    break 'pairs;
                }
            }
        }
        match descend {
            Some(deeper) => level = deeper,
            None => {
                level -= 1;
                if level > 0 {
                    // the level above sees the strong generators added below it
                    extend(&mut transversals[level - 1], &strong, level - 1);
                }
            }
        }
    }
    transversals
        .iter()
        .map(|t| t.iter().filter(|u| u.is_some()).count() as u128)
        .product()
}

#[test]
fn schreier_sims_orders() {
    // S_5 from a transposition and a 5-cycle, the dihedral group of the 6-gon, the trivial group
    assert_eq!(
        group_order(&[vec![1, 0, 2, 3, 4], vec![1, 2, 3, 4, 0]], 5),
        120
    );
    assert_eq!(
        group_order(&[vec![1, 2, 3, 4, 5, 0], vec![0, 5, 4, 3, 2, 1]], 6),
        12
    );
    assert_eq!(group_order(&[], 4), 1);
}

#[test]
fn schreier_sims_on_large_groups() {
    // S_3 wr S_11 on 11 disjoint triangles, generated by one rotation, one reflection and
    // swaps of neighbouring triangles: 6^11 11!
    let n = 33;
    let swap = |k: usize| {
        let mut p: Vec<usize> = (0..n).collect();
        for t in 0..3 {
            p.swap(3 * k + t, 3 * k + 3 + t);
        }
        p
    };
    let mut generators = vec![
        (0..n)
            .map(|v| if v < 3 { (v + 1) % 3 } else { v })
            .collect(),
        (0..n)
            .map(|v| [0, 2, 1].get(v).copied().unwrap_or(v))
            .collect(),
    ];
    generators.extend((0..10).map(swap));
    assert_eq!(group_order(&generators, n), 6u128.pow(11) * 39916800);
}
