//! Role-normalized flat view of a [`HedgeGraph`] used by the canonizer.
//!
//! Caller keys are interned to dense ranks (positions in the sorted list of distinct keys),
//! which are isomorphism invariant. Every hedge gets one colour combining the rank of its edge
//! key, the rank of its own key and its physical role, so the arbitrary stored Source/Sink role
//! of an undirected edge never matters and `Orientation::Reversed` resolves to the physical
//! direction.

use crate::half_edge::{
    involution::{Flow, InvolutiveMapping, Orientation},
    nodestore::NodeStorageOps,
    HedgeGraph,
};

/// Marks the missing opposite of a dangling hedge.
pub(crate) const NONE: u32 = u32::MAX;

const UNDIRECTED: u32 = 0;
const TAIL: u32 = 1;
const HEAD: u32 = 2;
const LEG: u32 = 3;
const LEG_OUT: u32 = 4;
const LEG_IN: u32 = 5;
const ROLES: u32 = 6;

/// Hedges, their pairing, incidence, colours and a node -> hedges index, in `u32` words.
#[derive(Debug, Default)]
pub(crate) struct Flat {
    /// Opposite hedge; dangling hedges are fixed points.
    pub(crate) inv: Vec<u32>,
    pub(crate) node_of: Vec<u32>,
    pub(crate) edge_of: Vec<u32>,
    pub(crate) hcol: Vec<u32>,
    pub(crate) ncol: Vec<u32>,
    /// Hedges of node `v` are `adj[off[v]..off[v + 1]]`.
    pub(crate) off: Vec<u32>,
    pub(crate) adj: Vec<u32>,
    /// Certificate prefix fixing the meaning of `hcol`: empty for the mixed-radix colours,
    /// otherwise the sorted `(edge key, hedge key, role)` table that dense colours index.
    pub(crate) colour_table: Vec<u32>,
    edge_key: Vec<u32>,
    hedge_key: Vec<u32>,
}

/// Sorted distinct caller keys; the canonical colours are ranks into these tables.
pub(crate) struct KeyTables<KV, KE, KH> {
    pub(crate) nodes: Vec<KV>,
    pub(crate) edges: Vec<KE>,
    pub(crate) hedges: Vec<KH>,
}

impl Flat {
    pub(crate) fn n_nodes(&self) -> usize {
        self.ncol.len()
    }

    pub(crate) fn n_hedges(&self) -> usize {
        self.inv.len()
    }

    pub(crate) fn n_edges(&self) -> usize {
        self.edge_key.len()
    }

    pub(crate) fn hedges_of(&self, v: usize) -> &[u32] {
        &self.adj[self.off[v] as usize..self.off[v + 1] as usize]
    }

    pub(crate) fn load<'a, E, V, H, N, KV: Ord, KE: Ord, KH: Ord>(
        &mut self,
        graph: &'a HedgeGraph<E, V, H, N>,
        mut node_key: impl FnMut(&'a V) -> KV,
        mut edge_key: impl FnMut(&'a E) -> KE,
        mut hedge_key: impl FnMut(&'a H) -> KH,
    ) -> KeyTables<KV, KE, KH>
    where
        N: NodeStorageOps<NodeData = V>,
    {
        let keys = KeyTables {
            nodes: rank(
                graph.node_store.iter().map(|(v, d)| (v.0, node_key(d))),
                &mut self.ncol,
            ),
            edges: rank(
                graph
                    .edge_store
                    .data
                    .iter()
                    .map(|(e, (d, _))| (e.0, edge_key(d))),
                &mut self.edge_key,
            ),
            hedges: rank(
                graph.hedge_data.iter().map(|(h, d)| (h.0, hedge_key(d))),
                &mut self.hedge_key,
            ),
        };

        let n_hedges = graph.n_hedges();
        self.inv.clear();
        self.node_of.clear();
        self.edge_of.clear();
        self.hcol.clear();
        let mapping = &graph.edge_store.involution.inv;
        for (h, m) in mapping.iter() {
            let (opposite, edge, role) = match m {
                InvolutiveMapping::Identity { data, underlying } => {
                    let role = match (data.orientation, underlying) {
                        (Orientation::Undirected, _) => LEG,
                        (Orientation::Default, Flow::Source)
                        | (Orientation::Reversed, Flow::Sink) => LEG_OUT,
                        _ => LEG_IN,
                    };
                    (h, data.data, role)
                }
                InvolutiveMapping::Source { data, sink_idx } => {
                    (*sink_idx, data.data, edge_role(data.orientation, true))
                }
                InvolutiveMapping::Sink { source_idx } => {
                    let InvolutiveMapping::Source { data, .. } = &mapping[*source_idx] else {
                        unreachable!("a sink hedge pairs with a source hedge")
                    };
                    (*source_idx, data.data, edge_role(data.orientation, false))
                }
            };
            self.inv.push(opposite.0 as u32);
            self.node_of.push(graph.node_id(h).0 as u32);
            self.edge_of.push(edge.0 as u32);
            // provisional colour: the role; `colour` folds in the keys below
            self.hcol.push(role);
        }
        self.colour(keys.edges.len(), keys.hedges.len());

        let n_nodes = self.ncol.len();
        self.off.clear();
        self.off.resize(n_nodes + 1, 0);
        for &v in &self.node_of {
            self.off[v as usize + 1] += 1;
        }
        for v in 0..n_nodes {
            self.off[v + 1] += self.off[v];
        }
        self.adj.clear();
        self.adj.resize(n_hedges, 0);
        for (h, &v) in self.node_of.iter().enumerate() {
            let slot = &mut self.off[v as usize];
            self.adj[*slot as usize] = h as u32;
            *slot += 1;
        }
        // the fill advanced every offset by its node's degree: shift back
        for v in (1..=n_nodes).rev() {
            self.off[v] = self.off[v - 1];
        }
        if n_nodes > 0 {
            self.off[0] = 0;
        }
        keys
    }

    /// Replace the provisional roles in `hcol` by full hedge colours.
    fn colour(&mut self, n_edge_keys: usize, n_hedge_keys: usize) {
        self.colour_table.clear();
        let radix = (n_hedge_keys as u64) * ROLES as u64;
        if (n_edge_keys as u64).saturating_mul(radix) <= u32::MAX as u64 {
            for ((c, &e), &k) in self.hcol.iter_mut().zip(&self.edge_of).zip(&self.hedge_key) {
                *c += self.edge_key[e as usize] * radix as u32 + k * ROLES;
            }
            return;
        }
        // Too many distinct keys for one word: index the sorted table of present colours, which
        // then becomes part of the certificate.
        let packed = |c: u32, e: u32, k: u32| ((e as u128) << 64) | ((k as u128) << 32) | c as u128;
        let mut present: Vec<u128> = (0..self.hcol.len())
            .map(|h| {
                packed(
                    self.hcol[h],
                    self.edge_key[self.edge_of[h] as usize],
                    self.hedge_key[h],
                )
            })
            .collect();
        present.sort_unstable();
        present.dedup();
        for h in 0..self.hcol.len() {
            let p = packed(
                self.hcol[h],
                self.edge_key[self.edge_of[h] as usize],
                self.hedge_key[h],
            );
            self.hcol[h] = present.binary_search(&p).expect("colour is present") as u32;
        }
        self.colour_table.push(present.len() as u32);
        for p in present {
            self.colour_table
                .extend([(p >> 64) as u32, (p >> 32) as u32, p as u32]);
        }
    }
}

/// One twin class of hedges at a node, for the kernel (automorphisms fixing every node).
pub(crate) enum Twins<'a> {
    /// Identical dangling legs.
    Legs(&'a [u32]),
    /// Parallel edges or asymmetric self-loops, by their hedge at this node; counted once.
    Edges(&'a [u32]),
    /// Self-loops with equal end colours, each by its smaller hedge; `2^m m!` symmetries.
    SymmetricLoops(&'a [u32]),
}

impl Flat {
    /// Hedges at `v` sorted by `(colour, opposite node, opposite colour, hedge)`, with the
    /// opposite node given by `rank` (identity for input labels).
    fn sorted_hedges(&self, v: usize, rank: impl Fn(u32) -> u32, out: &mut Vec<([u32; 3], u32)>) {
        out.clear();
        for &h in self.hedges_of(v) {
            let o = self.inv[h as usize];
            let key = if o == h {
                [self.hcol[h as usize], NONE, NONE]
            } else {
                [
                    self.hcol[h as usize],
                    rank(self.node_of[o as usize]),
                    self.hcol[o as usize],
                ]
            };
            out.push((key, h));
        }
        out.sort_unstable();
    }

    /// Visit every twin class with at least two members.
    pub(crate) fn twin_classes(
        &self,
        scratch: &mut Vec<([u32; 3], u32)>,
        hedges: &mut Vec<u32>,
        mut visit: impl FnMut(Twins<'_>),
    ) {
        for v in 0..self.n_nodes() {
            self.sorted_hedges(v, |w| w, scratch);
            for run in scratch.chunk_by(|a, b| a.0 == b.0) {
                if run.len() < 2 {
                    continue;
                }
                let [c, w, c2] = run[0].0;
                hedges.clear();
                hedges.extend(run.iter().map(|&(_, h)| h));
                if w == NONE {
                    visit(Twins::Legs(hedges));
                } else if w as usize == v && c == c2 {
                    hedges.retain(|&h| h < self.inv[h as usize]);
                    visit(Twins::SymmetricLoops(hedges));
                } else if (v as u32, c) < (w, c2) {
                    visit(Twins::Edges(hedges));
                }
            }
        }
    }

    /// Canonical hedge positions for the node order `lab`: hedges sorted by (canonical node,
    /// colour, opposite canonical node, opposite colour). Twins are placed so that the relabelled
    /// involution is determined by the sorted keys: a twin run whose partners are already placed
    /// follows their order, and symmetric self-loops are placed as adjacent pairs.
    pub(crate) fn hedge_order(
        &self,
        lab: &[u32],
        canonical_index: &mut Vec<u32>,
        scratch: &mut Vec<([u32; 3], u32)>,
        pos: &mut Vec<u32>,
    ) {
        canonical_index.resize(self.n_nodes(), 0);
        for (q, &v) in lab.iter().enumerate() {
            canonical_index[v as usize] = q as u32;
        }
        pos.clear();
        pos.resize(self.n_hedges(), NONE);
        let mut next = 0;
        for (q, &v) in lab.iter().enumerate() {
            self.sorted_hedges(v as usize, |w| canonical_index[w as usize], scratch);
            for run in scratch.chunk_by_mut(|a, b| a.0 == b.0) {
                let [c, w, c2] = run[0].0;
                let q = q as u32;
                if w != NONE && w == q && c == c2 {
                    for &(_, h) in run.iter() {
                        if pos[h as usize] == NONE {
                            pos[h as usize] = next;
                            pos[self.inv[h as usize] as usize] = next + 1;
                            next += 2;
                        }
                    }
                    continue;
                }
                if w != NONE && (w < q || (w == q && c > c2)) {
                    run.sort_unstable_by_key(|&(_, h)| pos[self.inv[h as usize] as usize]);
                }
                for &(_, h) in run.iter() {
                    pos[h as usize] = next;
                    next += 1;
                }
            }
        }
    }
}

fn edge_role(orientation: Orientation, source: bool) -> u32 {
    match (orientation, source) {
        (Orientation::Undirected, _) => UNDIRECTED,
        (Orientation::Default, true) | (Orientation::Reversed, false) => TAIL,
        _ => HEAD,
    }
}

/// Rank `keys` (given with their dense indices) into `ranks`; returns the sorted distinct keys.
fn rank<K: Ord>(keys: impl Iterator<Item = (usize, K)>, ranks: &mut Vec<u32>) -> Vec<K> {
    let mut sorted: Vec<(K, u32)> = keys.map(|(i, k)| (k, i as u32)).collect();
    sorted.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    ranks.clear();
    ranks.resize(sorted.len(), 0);
    let mut table: Vec<K> = Vec::new();
    for (k, i) in sorted {
        if table.last() != Some(&k) {
            table.push(k);
        }
        ranks[i as usize] = table.len() as u32 - 1;
    }
    table
}
