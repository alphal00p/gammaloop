//! # Canonical forms and automorphisms of half-edge graphs
//!
//! Symmetries are defined on half-edges, as in Borinsky's `(H, nu, iota)` model: an automorphism
//! is a permutation of the hedges that commutes with the involution (dangling hedges stay
//! dangling), preserves the node partition, and preserves the caller-supplied node, edge and
//! half-edge keys, together with the induced node permutation (which also permutes isolated
//! nodes of equal key). Edge orientation is physical: the stored Source/Sink role of an
//! undirected edge never matters and [`Orientation::Reversed`](super::involution::Orientation)
//! resolves to the actual direction; dangling legs carry their direction too. Payloads not
//! covered by a key do not distinguish anything, so dangling legs with equal keys at one node are
//! interchangeable.
//!
//! The automorphism count therefore includes `m!` for `m` interchangeable parallel edges or legs
//! and `2` for every undirected self-loop whose two ends are indistinguishable, matching
//! Graphica's convention for the equivalent plain graph.
//!
//! The algorithm is individualization-refinement on the node partition (equitable refinement by
//! edge-end colours, target cells chosen like nauty's `bestcell`, first-path and best-leaf
//! search, trace and orbit pruning). Hedge orders are derived at leaves and symmetries fixing
//! every node (twin hedges) are written down directly.
//!
//! ```
//! use linnet::half_edge::{
//!     builder::HedgeGraphBuilder, canon::ignore, involution::Orientation, HedgeGraph,
//! };
//!
//! // a theta graph: two nodes joined by three parallel edges
//! let mut builder = HedgeGraphBuilder::new();
//! let (a, b) = (builder.add_node(()), builder.add_node(()));
//! for _ in 0..3 {
//!     builder.add_edge(a, b, (), Orientation::Undirected);
//! }
//! let theta: HedgeGraph<(), ()> = builder.build();
//!
//! let canon = theta.canonize();
//! // swap the two nodes, permute the three edges
//! assert_eq!(canon.automorphism_count(), Some(12));
//! assert!(theta.is_isomorphic_by(&theta, ignore, ignore, ignore));
//! ```

mod flat;
mod search;
#[cfg(test)]
mod tests;

use std::cell::RefCell;

use flat::{Flat, Twins};
pub(crate) use flat::{KeyTables, RawHedge, Role};
pub(crate) use search::mix;
use search::Search;

use super::{
    involution::{EdgeIndex, EdgeVec, Flow, Hedge, HedgeVec, Orientation},
    nodestore::Topology,
    HedgeGraph, NoData, NodeIndex, NodeVec,
};
use crate::permutation::Permutation;

/// Key that ignores its payload, for [`HedgeGraph::canonize_by`] and friends.
pub fn ignore<T: ?Sized>(_: &T) -> NoData {
    NoData {}
}

/// Canonical labelling, canonical form and automorphism group of a graph, from
/// [`HedgeGraph::canonize_by`].
///
/// Applying [`Self::hedge_permutation`] and [`Self::node_permutation`] to isomorphic graphs gives
/// identical relabelled graphs (same incidence, pairing, keys and physical orientations).
#[derive(Clone, Debug)]
pub struct Canonization<KV, KE, KH> {
    canonical_form: CanonicalForm<KV, KE, KH>,
    hedge_map: Vec<u32>,
    node_map: Vec<u32>,
    edge_map: Vec<u32>,
    edge_of: Vec<u32>,
    generator_hedges: Vec<u32>,
    generator_nodes: Vec<u32>,
    automorphism_count: Option<u128>,
}

/// Complete isomorphism invariant: equal iff the graphs are isomorphic under the same keys.
///
/// `words` encodes the canonically relabelled structure, with colours that are ranks into the
/// sorted distinct key tables, which are part of the canonical form. Forms are comparable only
/// within one version of the canonizer.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CanonicalForm<KV, KE, KH> {
    words: Vec<u32>,
    node_keys: Vec<KV>,
    edge_keys: Vec<KE>,
    hedge_keys: Vec<KH>,
}

/// One edge in a compact canonicalization input, with caller keys rather than ranks.
///
/// Elements follow insertion order. A paired edge contributes its source hedge then its sink;
/// a dangling edge contributes one hedge, exactly as [`super::builder::HedgeGraphBuilder`].
#[derive(Clone, Debug)]
pub enum KeyedEdge<KE, KH> {
    Paired {
        source: (NodeIndex, KH),
        sink: (NodeIndex, KH),
        key: KE,
        orientation: Orientation,
    },
    Dangling {
        node: NodeIndex,
        hedge_key: KH,
        key: KE,
        orientation: Orientation,
        flow: Flow,
    },
}

/// Invalid compact input. Pairing and dense edge/hedge identifiers are supplied by the loader.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
pub enum CanonicalInputError {
    #[error("canonicalization node {node:?} is outside {nodes} nodes")]
    UnknownNode { node: NodeIndex, nodes: usize },
    #[error("canonicalization input exceeds the u32 index space")]
    IndexOverflow,
}

impl<KV, KE, KH> CanonicalForm<KV, KE, KH> {
    /// The canonical form and automorphism count of a graph given by its node colours (ranks of
    /// the node keys) and hedges; equal to those of [`HedgeGraph::canonize_by`] on the same graph
    /// with keys whose sorted distinct values are `keys`.
    pub(crate) fn of_hedges(
        node_colours: &[u32],
        hedges: &[RawHedge],
        n_edges: usize,
        keys: KeyTables<KV, KE, KH>,
    ) -> (Self, Option<u128>) {
        Workspace::with(|ws| {
            ws.flat.load_hedges(
                node_colours,
                hedges,
                n_edges,
                keys.edges.len(),
                keys.hedges.len(),
            );
            ws.certify(keys)
        })
    }

    /// Canonical form and input-node -> canonical-node map of a compact keyed graph.
    ///
    /// This uses the same search, key ordering and hedge insertion order as
    /// [`HedgeGraph::canonize_by`], without constructing a graph, hedge/edge maps or the lifted
    /// half-edge automorphism group. Node keys are in input node order, including isolated nodes.
    pub fn label_edges(
        node_keys: impl IntoIterator<Item = KV, IntoIter: ExactSizeIterator>,
        edges: impl IntoIterator<Item = KeyedEdge<KE, KH>>,
    ) -> Result<(Self, NodeVec<NodeIndex>), CanonicalInputError>
    where
        KV: Ord,
        KE: Ord,
        KH: Ord,
    {
        Workspace::with(|ws| {
            let keys = ws.flat.load_edges(node_keys, edges)?;
            let form = ws.canonical_form(keys);
            let mut nodes = vec![NodeIndex(0); ws.flat.n_nodes()];
            for (canonical, &input) in ws.search.best_lab().iter().enumerate() {
                nodes[input as usize] = NodeIndex(canonical);
            }
            Ok((form, nodes.into_iter().collect()))
        })
    }

    /// The structural part, for hashing or storage next to the key tables.
    pub fn words(&self) -> &[u32] {
        &self.words
    }

    pub fn node_keys(&self) -> &[KV] {
        &self.node_keys
    }

    pub fn edge_keys(&self) -> &[KE] {
        &self.edge_keys
    }

    pub fn hedge_keys(&self) -> &[KH] {
        &self.hedge_keys
    }

    /// The words and the node, edge and half-edge key tables, by value.
    pub fn into_parts(self) -> (Vec<u32>, Vec<KV>, Vec<KE>, Vec<KH>) {
        (self.words, self.node_keys, self.edge_keys, self.hedge_keys)
    }
}

/// One automorphism: a hedge permutation and the node permutation it induces.
#[derive(Clone, Copy, Debug)]
pub struct Automorphism<'a> {
    hedges: &'a [u32],
    nodes: &'a [u32],
}

impl Automorphism<'_> {
    pub fn hedge(&self, hedge: Hedge) -> Hedge {
        Hedge(self.hedges[hedge.0] as usize)
    }

    pub fn node(&self, node: NodeIndex) -> NodeIndex {
        NodeIndex(self.nodes[node.0] as usize)
    }

    pub fn hedge_permutation(&self) -> Permutation {
        permutation(self.hedges)
    }

    pub fn node_permutation(&self) -> Permutation {
        permutation(self.nodes)
    }
}

/// An isomorphism between two graphs, as maps from the first graph's elements to the second's.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Isomorphism {
    pub hedges: HedgeVec<Hedge>,
    pub nodes: NodeVec<NodeIndex>,
    pub edges: EdgeVec<EdgeIndex>,
}

fn permutation(map: &[u32]) -> Permutation {
    Permutation::from_map(map.iter().map(|&x| x as usize).collect())
}

impl<KV, KE, KH> Canonization<KV, KE, KH> {
    pub fn canonical_form(&self) -> &CanonicalForm<KV, KE, KH> {
        &self.canonical_form
    }

    pub fn into_canonical_form(self) -> CanonicalForm<KV, KE, KH> {
        self.canonical_form
    }

    /// Position of `hedge` in the canonical labelling.
    pub fn canonical_hedge(&self, hedge: Hedge) -> Hedge {
        Hedge(self.hedge_map[hedge.0] as usize)
    }

    pub fn canonical_node(&self, node: NodeIndex) -> NodeIndex {
        NodeIndex(self.node_map[node.0] as usize)
    }

    /// Edges are ordered by their first canonical hedge.
    pub fn canonical_edge(&self, edge: EdgeIndex) -> EdgeIndex {
        EdgeIndex(self.edge_map[edge.0] as usize)
    }

    /// Input hedge -> canonical hedge.
    pub fn hedge_permutation(&self) -> Permutation {
        permutation(&self.hedge_map)
    }

    pub fn node_permutation(&self) -> Permutation {
        permutation(&self.node_map)
    }

    pub fn edge_permutation(&self) -> Permutation {
        permutation(&self.edge_map)
    }

    /// Order of the automorphism group, or `None` if it exceeds `u128`.
    pub fn automorphism_count(&self) -> Option<u128> {
        self.automorphism_count
    }

    /// Generators of the automorphism group: twin swaps and self-loop flips first, then the
    /// automorphisms found by the search, then swaps of isolated nodes.
    pub fn generators(&self) -> impl ExactSizeIterator<Item = Automorphism<'_>> {
        let (n_hedges, n_nodes) = (self.hedge_map.len(), self.node_map.len());
        // hedges need nodes, so a graph without nodes has no generators
        let count = self.generator_nodes.len().checked_div(n_nodes).unwrap_or(0);
        (0..count).map(move |i| Automorphism {
            hedges: &self.generator_hedges[i * n_hedges..(i + 1) * n_hedges],
            nodes: &self.generator_nodes[i * n_nodes..(i + 1) * n_nodes],
        })
    }

    /// Hedge orbits, each hedge mapped to the smallest hedge of its orbit.
    pub fn hedge_orbits(&self) -> HedgeVec<Hedge> {
        let mut orbits = Orbits::default();
        orbits.reset(self.hedge_map.len());
        for g in self.generators() {
            orbits.join(g.hedges);
        }
        (0..self.hedge_map.len() as u32)
            .map(|h| Hedge(orbits.find(h) as usize))
            .collect()
    }

    /// Node orbits, each node mapped to the smallest node of its orbit.
    pub fn node_orbits(&self) -> NodeVec<NodeIndex> {
        let mut orbits = Orbits::default();
        orbits.reset(self.node_map.len());
        for g in self.generators() {
            orbits.join(g.nodes);
        }
        (0..self.node_map.len() as u32)
            .map(|v| NodeIndex(orbits.find(v) as usize))
            .collect()
    }

    /// Edge orbits, each edge mapped to the smallest edge of its orbit.
    pub fn edge_orbits(&self) -> EdgeVec<EdgeIndex> {
        let mut orbits = Orbits::default();
        orbits.reset(self.edge_map.len());
        for g in self.generators() {
            for (h, &image) in g.hedges.iter().enumerate() {
                orbits.union(self.edge_of[h], self.edge_of[image as usize]);
            }
        }
        (0..self.edge_map.len() as u32)
            .map(|e| EdgeIndex(orbits.find(e) as usize))
            .collect()
    }

    /// The isomorphism from this graph to `other`, if their canonical forms agree.
    pub fn isomorphism_to(&self, other: &Self) -> Option<Isomorphism>
    where
        KV: PartialEq,
        KE: PartialEq,
        KH: PartialEq,
    {
        if self.canonical_form != other.canonical_form {
            return None;
        }
        let compose = |from: &[u32], to: &[u32]| {
            let mut inverse = vec![0usize; to.len()];
            for (x, &c) in to.iter().enumerate() {
                inverse[c as usize] = x;
            }
            from.iter()
                .map(move |&c| inverse[c as usize])
                .collect::<Vec<_>>()
        };
        Some(Isomorphism {
            hedges: compose(&self.hedge_map, &other.hedge_map)
                .into_iter()
                .map(Hedge)
                .collect(),
            nodes: compose(&self.node_map, &other.node_map)
                .into_iter()
                .map(NodeIndex)
                .collect(),
            edges: compose(&self.edge_map, &other.edge_map)
                .into_iter()
                .map(EdgeIndex)
                .collect(),
        })
    }
}

/// Union-find over `0..n` whose roots are the smallest element of their set, on a reusable
/// buffer: the search recomputes node orbits in its inner loop, which the data-carrying
/// [`crate::union_find::UnionFind`] is not built for.
#[derive(Debug, Default)]
pub(crate) struct Orbits {
    parent: Vec<u32>,
}

impl Orbits {
    pub(crate) fn reset(&mut self, n: usize) {
        self.parent.clear();
        self.parent.extend(0..n as u32);
    }

    pub(crate) fn find(&mut self, mut x: u32) -> u32 {
        let p = &mut self.parent;
        while p[x as usize] != x {
            p[x as usize] = p[p[x as usize] as usize];
            x = p[x as usize];
        }
        x
    }

    pub(crate) fn union(&mut self, a: u32, b: u32) {
        let (a, b) = (self.find(a), self.find(b));
        self.parent[a.max(b) as usize] = a.min(b);
    }

    /// Merge every point with its image under `permutation`.
    pub(crate) fn join(&mut self, permutation: &[u32]) {
        for (x, &y) in permutation.iter().enumerate() {
            self.union(x as u32, y);
        }
    }
}

/// Node canonization of a graph given by node colours and coloured undirected edges, for searches
/// that keep their own light graph state (the graph generator): the canonical node order and the
/// node automorphisms found, without hedge maps, canonical forms or twin handling. Nodes without
/// edges are not permuted among themselves.
#[derive(Debug, Default)]
pub(crate) struct NodeCanon {
    flat: Flat,
    search: Search,
}

impl NodeCanon {
    pub(crate) fn run(&mut self, node_colours: &[u32], edges: &[[u32; 3]]) {
        self.flat.load_coloured(node_colours, edges);
        self.search.run(&self.flat);
    }

    /// Nodes in canonical order.
    pub(crate) fn lab(&self) -> &[u32] {
        self.search.best_lab()
    }

    /// Order of the node group for the connected graphs used by generation; `None` on overflow.
    pub(crate) fn group_order(&self) -> Option<u128> {
        self.search.group_order
    }

    /// Generators of the node automorphism group, each as the image of every node.
    pub(crate) fn generators(&self) -> impl Iterator<Item = &[u32]> {
        (0..self.search.n_generators()).map(|i| self.search.generator(i))
    }
}

/// `n!`, or `None` beyond `u128`.
pub(crate) fn factorial(n: usize) -> Option<u128> {
    (2..=n as u128).try_fold(1u128, |acc, k| acc.checked_mul(k))
}

/// A group order times a factor, `None` once either overflows.
pub(crate) fn scaled(order: Option<u128>, factor: Option<u128>) -> Option<u128> {
    order.zip(factor).and_then(|(a, b)| a.checked_mul(b))
}

/// Automorphism generators stored back to back, and the group order (`None` beyond `u128`).
struct Group<'a> {
    n_hedges: usize,
    n_nodes: usize,
    hedges: &'a mut Vec<u32>,
    nodes: &'a mut Vec<u32>,
    order: Option<u128>,
}

impl Group<'_> {
    /// A generator that swaps the given hedge pairs and the given node pairs.
    fn push_swaps(&mut self, hedges: &[(u32, u32)], nodes: &[(u32, u32)]) {
        for (out, n, swaps) in [
            (&mut *self.hedges, self.n_hedges, hedges),
            (&mut *self.nodes, self.n_nodes, nodes),
        ] {
            let start = out.len();
            out.extend(0..n as u32);
            for &(a, b) in swaps {
                out.swap(start + a as usize, start + b as usize);
            }
        }
    }

    fn scale(&mut self, factor: Option<u128>) {
        self.order = scaled(self.order, factor);
    }

    /// Twin hedges: the automorphisms fixing every node.
    fn add_twins(&mut self, class: Twins<'_>, inv: &[u32]) {
        let paired = |a: u32, b: u32| [(a, b), (inv[a as usize], inv[b as usize])];
        self.scale(class.order());
        match class {
            Twins::Legs(hs) => {
                for w in hs.windows(2) {
                    self.push_swaps(&[(w[0], w[1])], &[]);
                }
            }
            Twins::Edges(hs) => {
                for w in hs.windows(2) {
                    self.push_swaps(&paired(w[0], w[1]), &[]);
                }
            }
            Twins::SymmetricLoops(hs) => {
                // one flip and the loop swaps generate the wreath product Z2 wr S_m
                self.push_swaps(&[(hs[0], inv[hs[0] as usize])], &[]);
                for w in hs.windows(2) {
                    self.push_swaps(&paired(w[0], w[1]), &[]);
                }
            }
        }
    }
}

/// Reusable buffers, one per thread.
#[derive(Debug, Default)]
struct Workspace {
    flat: Flat,
    search: Search,
    canonical_index: Vec<u32>,
    scratch: Vec<([u32; 3], u32)>,
    twins: Vec<u32>,
    pos: Vec<u32>,
    image_lab: Vec<u32>,
    image_pos: Vec<u32>,
    at: Vec<u32>,
    generator_hedges: Vec<u32>,
    generator_nodes: Vec<u32>,
}

thread_local! {
    static WORKSPACE: RefCell<Workspace> = RefCell::default();
}

impl Workspace {
    fn with<T>(f: impl FnOnce(&mut Workspace) -> T) -> T {
        WORKSPACE.with(|ws| match ws.try_borrow_mut() {
            Ok(mut ws) => f(&mut ws),
            // a key callback that canonizes again gets its own buffers
            Err(_) => f(&mut Workspace::default()),
        })
    }

    fn canonize<KV, KE, KH>(&mut self, keys: KeyTables<KV, KE, KH>) -> Canonization<KV, KE, KH> {
        let Workspace {
            flat: g,
            search,
            canonical_index,
            scratch,
            twins,
            pos,
            image_lab,
            image_pos,
            at,
            generator_hedges,
            generator_nodes,
        } = self;
        let g: &Flat = g;
        search.run(g);
        let (n, n_hedges) = (g.n_nodes(), g.n_hedges());
        let lab = search.best_lab();

        let words = Self::words(g, search);

        let mut node_map = vec![0; n];
        for (q, &v) in lab.iter().enumerate() {
            node_map[v as usize] = q as u32;
        }
        g.hedge_order(lab, canonical_index, scratch, pos);
        at.clear();
        at.resize(n_hedges, 0);
        for (h, &p) in pos.iter().enumerate() {
            at[p as usize] = h as u32;
        }
        let mut edge_map = vec![u32::MAX; g.n_edges()];
        let mut next = 0;
        for &h in at.iter() {
            let e = &mut edge_map[g.edge_of[h as usize] as usize];
            if *e == u32::MAX {
                *e = next;
                next += 1;
            }
        }

        generator_hedges.clear();
        generator_nodes.clear();
        let mut group = Group {
            n_hedges,
            n_nodes: n,
            hedges: generator_hedges,
            nodes: generator_nodes,
            order: search.group_order,
        };
        g.twin_classes(scratch, twins, |class| group.add_twins(class, &g.inv));
        // node automorphisms from the search, lifted to hedges through the canonical hedge
        // order of the image labelling
        for i in 0..search.n_generators() {
            let gamma = search.generator(i);
            image_lab.clear();
            image_lab.extend(lab.iter().map(|&v| gamma[v as usize]));
            g.hedge_order(image_lab, canonical_index, scratch, image_pos);
            for (h, &p) in image_pos.iter().enumerate() {
                at[p as usize] = h as u32;
            }
            group.hedges.extend(pos.iter().map(|&p| at[p as usize]));
            group.nodes.extend_from_slice(gamma);
        }
        for run in g.isolated_runs(lab) {
            group.scale(factorial(run.len()));
            for w in run.windows(2) {
                group.push_swaps(&[], &[(w[0], w[1])]);
            }
        }

        Canonization {
            canonical_form: CanonicalForm {
                words,
                node_keys: keys.nodes,
                edge_keys: keys.edges,
                hedge_keys: keys.hedges,
            },
            hedge_map: pos.clone(),
            node_map,
            edge_map,
            edge_of: g.edge_of.clone(),
            generator_hedges: group.hedges.clone(),
            generator_nodes: group.nodes.clone(),
            automorphism_count: group.order,
        }
    }

    /// The canonical-form words of the searched graph.
    fn words(g: &Flat, search: &Search) -> Vec<u32> {
        let lab = search.best_lab();
        let mut words =
            Vec::with_capacity(2 + g.colour_table.len() + lab.len() + search.best_cert().len());
        words.extend([g.n_nodes() as u32, g.n_hedges() as u32]);
        words.extend_from_slice(&g.colour_table);
        words.extend(lab.iter().map(|&v| g.ncol[v as usize]));
        words.extend_from_slice(search.best_cert());
        words
    }

    /// The canonical form only: no labelling maps, generators or automorphism count.
    fn canonical_form<KV, KE, KH>(
        &mut self,
        keys: KeyTables<KV, KE, KH>,
    ) -> CanonicalForm<KV, KE, KH> {
        self.search.run(&self.flat);
        CanonicalForm {
            words: Self::words(&self.flat, &self.search),
            node_keys: keys.nodes,
            edge_keys: keys.edges,
            hedge_keys: keys.hedges,
        }
    }

    /// The canonical form and the automorphism count, without labelling maps and generators.
    fn certify<KV, KE, KH>(
        &mut self,
        keys: KeyTables<KV, KE, KH>,
    ) -> (CanonicalForm<KV, KE, KH>, Option<u128>) {
        let form = self.canonical_form(keys);
        let Workspace {
            flat: g,
            search,
            scratch,
            twins,
            ..
        } = self;
        let mut order = search.group_order;
        g.twin_classes(scratch, twins, |class| order = scaled(order, class.order()));
        for run in g.isolated_runs(search.best_lab()) {
            order = scaled(order, factorial(run.len()));
        }
        (form, order)
    }
}

impl<E, V, H, N: Topology> HedgeGraph<E, V, H, N> {
    /// Canonize the graph under the given node, edge and half-edge keys.
    ///
    /// Elements with equal keys are interchangeable; pass [`ignore`] to disregard a payload.
    /// Dangling legs are distinguished only by their (edge and half-edge) keys and direction.
    pub fn canonize_by<'a, KV: Ord, KE: Ord, KH: Ord>(
        &'a self,
        node_key: impl FnMut(&'a V) -> KV,
        edge_key: impl FnMut(&'a E) -> KE,
        hedge_key: impl FnMut(&'a H) -> KH,
    ) -> Canonization<KV, KE, KH> {
        Workspace::with(|ws| {
            let keys = ws.flat.load(self, node_key, edge_key, hedge_key);
            ws.canonize(keys)
        })
    }

    /// Canonize the graph with every node, edge and half-edge payload as its key.
    pub fn canonize(&self) -> Canonization<&V, &E, &H>
    where
        V: Ord,
        E: Ord,
        H: Ord,
    {
        self.canonize_by(|v| v, |e| e, |h| h)
    }

    /// The canonical form of [`Self::canonize_by`], without the labelling maps, automorphism
    /// generators and automorphism count.
    pub fn canonical_form_by<'a, KV: Ord, KE: Ord, KH: Ord>(
        &'a self,
        node_key: impl FnMut(&'a V) -> KV,
        edge_key: impl FnMut(&'a E) -> KE,
        hedge_key: impl FnMut(&'a H) -> KH,
    ) -> CanonicalForm<KV, KE, KH> {
        Workspace::with(|ws| {
            let keys = ws.flat.load(self, node_key, edge_key, hedge_key);
            ws.canonical_form(keys)
        })
    }

    /// The canonical form with every node, edge and half-edge payload as its key.
    pub fn canonical_form(&self) -> CanonicalForm<&V, &E, &H>
    where
        V: Ord,
        E: Ord,
        H: Ord,
    {
        self.canonical_form_by(|v| v, |e| e, |h| h)
    }

    /// Whether `other` is isomorphic to this graph under the given keys.
    pub fn is_isomorphic_by<'a, KV: Ord, KE: Ord, KH: Ord, N2: Topology>(
        &'a self,
        other: &'a HedgeGraph<E, V, H, N2>,
        mut node_key: impl FnMut(&'a V) -> KV,
        mut edge_key: impl FnMut(&'a E) -> KE,
        mut hedge_key: impl FnMut(&'a H) -> KH,
    ) -> bool {
        self.n_live() == other.n_live()
            && self.n_hedges() == other.n_hedges()
            && self.canonical_form_by(&mut node_key, &mut edge_key, &mut hedge_key)
                == other.canonical_form_by(&mut node_key, &mut edge_key, &mut hedge_key)
    }

    /// An isomorphism from this graph to `other` under the given keys, if one exists.
    pub fn isomorphism_by<'a, KV: Ord, KE: Ord, KH: Ord, N2: Topology>(
        &'a self,
        other: &'a HedgeGraph<E, V, H, N2>,
        mut node_key: impl FnMut(&'a V) -> KV,
        mut edge_key: impl FnMut(&'a E) -> KE,
        mut hedge_key: impl FnMut(&'a H) -> KH,
    ) -> Option<Isomorphism> {
        if self.n_live() != other.n_live()
            || self.n_hedges() != other.n_hedges()
            || self.n_edges() != other.n_edges()
        {
            return None;
        }
        let a = self.canonize_by(&mut node_key, &mut edge_key, &mut hedge_key);
        let b = other.canonize_by(&mut node_key, &mut edge_key, &mut hedge_key);
        a.isomorphism_to(&b)
    }
}
