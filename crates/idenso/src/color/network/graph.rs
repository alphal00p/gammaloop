//! Closed colour networks: generator traces, structure constants and
//! symmetric invariant tensors, with canonical labelling and edge cuts.

use super::poly::{CA, I2R, NA, NR, Poly};
#[allow(unused_imports)]
use super::poly::CR;

pub(crate) type Ix = u16;

/// The (single) non-adjoint representation R, or the adjoint A.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub(crate) enum Rep {
    R,
    A,
}

impl Rep {
    /// Tr(1).
    pub(crate) fn dimension(self) -> u16 {
        match self {
            Rep::R => NR,
            Rep::A => NA,
        }
    }
    /// Tr(T^a T^b) = index * delta^{ab}.
    pub(crate) fn index(self) -> u16 {
        match self {
            Rep::R => I2R,
            Rep::A => CA,
        }
    }
    /// T^a T^a = casimir * 1.
    pub(crate) fn casimir(self) -> u16 {
        match self {
            Rep::R => super::poly::CR,
            Rep::A => CA,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) enum Obj {
    /// f^{abc}, [T^a, T^b] = i f^{abc} T^c.
    F([Ix; 3]),
    /// Tr(T^{a1} ... T^{an}); for `Rep::A` the generators are (F^a)_{bc} = -i f^{abc}.
    Tr(Rep, Vec<Ix>),
    /// The normalized symmetrized trace Str(T^{a1} ... T^{an}).
    D(Rep, Vec<Ix>),
    /// Tr(S(T^{a1} ... T^{as}) T^{b1} ... T^{bk}): a trace whose first `s`
    /// generators are symmetrized, during the expansion of a trace.
    St(Rep, usize, Vec<Ix>),
}

impl Obj {
    pub(crate) fn indices(&self) -> &[Ix] {
        match self {
            Obj::F(i) => i,
            Obj::Tr(_, i) | Obj::D(_, i) | Obj::St(_, _, i) => i,
        }
    }
    fn indices_mut(&mut self) -> &mut [Ix] {
        match self {
            Obj::F(i) => i,
            Obj::Tr(_, i) | Obj::D(_, i) | Obj::St(_, _, i) => i,
        }
    }
}

/// Every index occurs exactly twice.
#[derive(Clone, Debug, Default)]
pub(crate) struct Net {
    pub(crate) objs: Vec<Obj>,
    pub(crate) next: Ix,
}

impl Net {
    pub(crate) fn new(objs: Vec<Obj>) -> Net {
        let next = objs
            .iter()
            .flat_map(|o| o.indices().iter().copied())
            .max()
            .map_or(0, |m| m + 1);
        Net { objs, next }
    }

    pub(crate) fn fresh(&mut self) -> Ix {
        let i = self.next;
        self.next += 1;
        i
    }

    /// (object, position) of both occurrences of every index.
    pub(crate) fn occurrences(&self) -> Vec<[(u16, u16); 2]> {
        let mut occ = vec![[(u16::MAX, u16::MAX); 2]; self.next as usize];
        for (o, obj) in self.objs.iter().enumerate() {
            for (p, &i) in obj.indices().iter().enumerate() {
                let slot = &mut occ[i as usize];
                if slot[0].0 == u16::MAX {
                    slot[0] = (o as u16, p as u16);
                } else {
                    debug_assert_eq!(slot[1].0, u16::MAX, "index {i} occurs more than twice");
                    slot[1] = (o as u16, p as u16);
                }
            }
        }
        occ
    }

    pub(crate) fn rename(&mut self, from: Ix, to: Ix) {
        for obj in &mut self.objs {
            for i in obj.indices_mut() {
                if *i == from {
                    *i = to;
                }
            }
        }
    }

    pub(crate) fn count_f(&self) -> usize {
        self.objs.iter().filter(|o| matches!(o, Obj::F(_))).count()
    }

    /// Remove trivial objects: empty, one- and two-generator traces and
    /// symmetric tensors, odd adjoint symmetric tensors and structure
    /// constants with a repeated index. `None` for a vanishing network.
    pub(crate) fn normalize(mut self) -> Option<(Poly, Net)> {
        let mut factor = Poly::one();
        let mut k = 0;
        while k < self.objs.len() {
            let (rep, idx, symmetric) = match &self.objs[k] {
                Obj::F(i) => {
                    if i[0] == i[1] || i[1] == i[2] || i[0] == i[2] {
                        return None;
                    }
                    k += 1;
                    continue;
                }
                Obj::Tr(rep, idx) => (*rep, idx.clone(), false),
                Obj::D(rep, idx) => (*rep, idx.clone(), true),
                Obj::St(rep, split, idx) => {
                    if *split == idx.len() {
                        self.objs[k] = Obj::D(*rep, idx.clone());
                        continue;
                    }
                    if idx.len() > 2 {
                        k += 1;
                        continue;
                    }
                    (*rep, idx.clone(), false)
                }
            };
            if symmetric && rep == Rep::A && idx.len() % 2 == 1 {
                return None;
            }
            match idx.len() {
                0 => {
                    factor = factor.mul(&Poly::var(rep.dimension()));
                    self.objs.swap_remove(k);
                }
                1 => return None,
                2 => {
                    self.objs.swap_remove(k);
                    if idx[0] == idx[1] {
                        // Tr(T^a T^a) = C_X N_X, as the local rules write it.
                        factor = factor.mul(&Poly::var(rep.casimir())).mul(&Poly::var(rep.dimension()));
                    } else {
                        factor = factor.mul(&Poly::var(rep.index()));
                        self.rename(idx[1], idx[0]);
                    }
                    k = 0;
                }
                _ => k += 1,
            }
        }
        Some((factor, self))
    }
}

/// Node kinds of the labelled graph.
const K_TR: u32 = 1;
const K_TA: u32 = 2;
const K_F: u32 = 3;
const K_DR: u32 = 16;
const K_DA: u32 = 64;
/// The symmetrized block of a partially symmetrized trace, plus its rank.
const K_SR: u32 = 128;
const K_SA: u32 = 192;
/// The role of a slot of a structure constant or symmetric tensor.
const ROLE_SYM: u8 = 3;

/// The graph of a network: one node per trace position, structure constant
/// and symmetric tensor. Trace nodes have slots [adjoint, in, out].
pub(crate) struct NodeGraph {
    pub(crate) kind: Vec<u32>,
    /// Per node and slot: the connected (node, slot).
    pub(crate) adj: Vec<Vec<(u32, u8)>>,
    /// Per node: (object, position).
    pub(crate) origin: Vec<(u16, u16)>,
    /// For F nodes, the original index labels of the slots (tie breaks).
    pub(crate) labels: Vec<Vec<Ix>>,
}

impl NodeGraph {
    pub(crate) fn of(net: &Net) -> NodeGraph {
        let mut kind = Vec::new();
        let mut origin = Vec::new();
        let mut first = Vec::with_capacity(net.objs.len());
        for (o, obj) in net.objs.iter().enumerate() {
            first.push(kind.len() as u32);
            match obj {
                Obj::F(_) => {
                    kind.push(K_F);
                    origin.push((o as u16, 0));
                }
                Obj::Tr(rep, idx) => {
                    for p in 0..idx.len() {
                        kind.push(if *rep == Rep::R { K_TR } else { K_TA });
                        origin.push((o as u16, p as u16));
                    }
                }
                Obj::D(rep, idx) => {
                    kind.push(if *rep == Rep::R { K_DR } else { K_DA } + idx.len() as u32);
                    origin.push((o as u16, 0));
                }
                Obj::St(rep, split, idx) => {
                    kind.push(if *rep == Rep::R { K_SR } else { K_SA } + *split as u32);
                    origin.push((o as u16, 0));
                    for p in *split..idx.len() {
                        kind.push(if *rep == Rep::R { K_TR } else { K_TA });
                        origin.push((o as u16, p as u16));
                    }
                }
            }
        }
        let n = kind.len();
        let mut adj: Vec<Vec<(u32, u8)>> = vec![Vec::new(); n];
        let mut labels: Vec<Vec<Ix>> = vec![Vec::new(); n];
        for (o, obj) in net.objs.iter().enumerate() {
            match obj {
                Obj::F(i) => {
                    adj[first[o] as usize] = vec![(u32::MAX, 0); 3];
                    labels[first[o] as usize] = i.to_vec();
                }
                Obj::D(_, i) => adj[first[o] as usize] = vec![(u32::MAX, 0); i.len()],
                Obj::Tr(_, i) => {
                    let len = i.len() as u32;
                    for p in 0..len {
                        let node = first[o] + p;
                        let next = first[o] + (p + 1) % len;
                        let prev = first[o] + (p + len - 1) % len;
                        adj[node as usize] = vec![(u32::MAX, 0), (prev, 2), (next, 1)];
                    }
                }
                Obj::St(_, split, i) => {
                    // Block slots: [in, out, legs...]; trace nodes follow.
                    let block = first[o];
                    let ord = (i.len() - split) as u32;
                    let mut slots = vec![(block + ord, 2), (block + 1, 1)];
                    slots.extend(std::iter::repeat_n((u32::MAX, 0), *split));
                    adj[block as usize] = slots;
                    for p in 0..ord {
                        let node = block + 1 + p;
                        let prev = if p == 0 { (block, 1) } else { (node - 1, 2) };
                        let next = if p + 1 == ord { (block, 0) } else { (node + 1, 1) };
                        adj[node as usize] = vec![(u32::MAX, 0), prev, next];
                    }
                }
            }
        }
        // Adjoint connections.
        let occ = net.occurrences();
        let node_slot = |o: u16, p: u16| -> (u32, u8) {
            match &net.objs[o as usize] {
                Obj::Tr(..) => (first[o as usize] + p as u32, 0),
                Obj::St(_, split, _) if (p as usize) < *split => (first[o as usize], 2 + p as u8),
                Obj::St(_, split, _) => (first[o as usize] + 1 + (p as u32 - *split as u32), 0),
                _ => (first[o as usize], p as u8),
            }
        };
        for pair in occ.iter() {
            if pair[0].0 == u16::MAX {
                continue;
            }
            let a = node_slot(pair[0].0, pair[0].1);
            let b = node_slot(pair[1].0, pair[1].1);
            adj[a.0 as usize][a.1 as usize] = b;
            adj[b.0 as usize][b.1 as usize] = a;
        }
        NodeGraph { kind, adj, origin, labels }
    }

    fn role(&self, node: u32, slot: u8) -> u8 {
        match self.kind[node as usize] {
            K_TR | K_TA => slot,
            k if k >= K_SR && slot < 2 => slot + 1,
            _ => ROLE_SYM,
        }
    }

    /// Colour refinement. Neighbour multisets are summarized by a
    /// commutative hash: a collision only refines less, the certificate of a
    /// leaf is exact, so the canonical form stays exact and invariant.
    fn refine(&self, colors: &mut [u32]) {
        let n = self.kind.len();
        let mut classes = count_classes(colors);
        let mut keyed: Vec<(u32, u64, u32)> = Vec::with_capacity(n);
        loop {
            keyed.clear();
            for v in 0..n {
                let mut aggregate = 0u64;
                for (s, &(w, t)) in self.adj[v].iter().enumerate() {
                    let term = ((self.role(v as u32, s as u8) as u64) << 40)
                        | ((self.role(w, t) as u64) << 32)
                        | colors[w as usize] as u64;
                    aggregate = aggregate.wrapping_add(mix(term));
                }
                keyed.push((colors[v], mix(aggregate ^ colors[v] as u64), v as u32));
            }
            keyed.sort_unstable();
            let mut rank = 0u32;
            for k in 0..n {
                if k > 0 && (keyed[k].0, keyed[k].1) != (keyed[k - 1].0, keyed[k - 1].1) {
                    rank += 1;
                }
                colors[keyed[k].2 as usize] = rank;
            }
            let new = rank as usize + 1;
            if new == classes {
                return;
            }
            classes = new;
        }
    }

    /// The certificate of a discrete colouring and the sign of the network
    /// relative to it.
    fn certificate(&self, colors: &[u32]) -> (Vec<u32>, i8) {
        let n = self.kind.len();
        let label = colors;
        let mut order = vec![0u32; n];
        for v in 0..n {
            order[label[v] as usize] = v as u32;
        }
        let mut cert = Vec::with_capacity(n * 5);
        let mut sign = 1i8;
        for &v in &order {
            let v = v as usize;
            cert.push(self.kind[v]);
            match self.kind[v] {
                K_TR | K_TA => {
                    for &(w, t) in &self.adj[v] {
                        cert.push((label[w as usize] << 2) | self.role(w, t) as u32);
                    }
                }
                K_F => {
                    let mut legs: Vec<(u32, Ix, usize)> = self.adj[v]
                        .iter()
                        .enumerate()
                        .map(|(s, &(w, t))| {
                            ((label[w as usize] << 2) | self.role(w, t) as u32, self.labels[v][s], s)
                        })
                        .collect();
                    legs.sort_unstable();
                    if permutation_parity([legs[0].2, legs[1].2, legs[2].2]) {
                        sign = -sign;
                    }
                    cert.extend(legs.iter().map(|l| l.0));
                }
                k => {
                    let skip = if k >= K_SR { 2 } else { 0 };
                    for &(w, t) in &self.adj[v][..skip] {
                        cert.push((label[w as usize] << 2) | self.role(w, t) as u32);
                    }
                    let mut legs: Vec<u32> = self.adj[v][skip..]
                        .iter()
                        .map(|&(w, t)| (label[w as usize] << 2) | self.role(w, t) as u32)
                        .collect();
                    legs.sort_unstable();
                    cert.extend(legs);
                }
            }
        }
        (cert, sign)
    }
}

fn mix(mut x: u64) -> u64 {
    // splitmix64
    x = x.wrapping_add(0x9e37_79b9_7f4a_7c15);
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^ (x >> 31)
}

fn count_classes(colors: &[u32]) -> usize {
    let mut c = colors.to_vec();
    c.sort_unstable();
    c.dedup();
    c.len()
}

/// Whether the permutation of three is odd.
fn permutation_parity(p: [usize; 3]) -> bool {
    let mut inversions = 0;
    for i in 0..3 {
        for j in i + 1..3 {
            if p[i] > p[j] {
                inversions += 1;
            }
        }
    }
    inversions % 2 == 1
}

/// Result of canonical labelling.
pub(crate) struct Canonical {
    pub(crate) certificate: Vec<u32>,
    /// value(net) = sign * value(canonical); 0 for an odd automorphism that
    /// moves nodes. Not 0 for an odd symmetry that only exchanges parallel
    /// edges (a structure constant with two legs on one symmetric tensor,
    /// which the rules set to zero), nor when an odd automorphism lies beyond
    /// LEAF_LIMIT; the sign is then exact for the chosen labelling.
    pub(crate) sign: i8,
    /// The node colouring of the minimal leaf (node -> canonical label).
    pub(crate) labels: Vec<u32>,
}

/// Upper bound on explored leaves; beyond it the best leaf so far is used,
/// which is still an exact relabelling, only possibly not canonical.
const LEAF_LIMIT: usize = 20000;

pub(crate) fn canonical(graph: &NodeGraph) -> Canonical {
    let n = graph.kind.len();
    let mut colors: Vec<u32> = graph.kind.clone();
    // Rank the initial colours.
    let mut kinds = graph.kind.clone();
    kinds.sort_unstable();
    kinds.dedup();
    for c in colors.iter_mut() {
        *c = kinds.binary_search(c).unwrap() as u32;
    }
    graph.refine(&mut colors);
    let mut state = Search { best: None, zero: false, leaves: 0, orbits: (0..n as u32).collect() };
    search(graph, colors, n, &mut state, true);
    let (certificate, sign, labels) = state.best.unwrap();
    Canonical { certificate, sign: if state.zero { 0 } else { sign }, labels }
}

struct Search {
    best: Option<(Vec<u32>, i8, Vec<u32>)>,
    zero: bool,
    leaves: usize,
    /// Union-find of the orbits of the automorphisms found so far.
    orbits: Vec<u32>,
}

fn orbit(orbits: &mut [u32], x: u32) -> u32 {
    let mut r = x;
    while orbits[r as usize] != r {
        r = orbits[r as usize];
    }
    let mut y = x;
    while orbits[y as usize] != r {
        let next = orbits[y as usize];
        orbits[y as usize] = r;
        y = next;
    }
    r
}

fn search(graph: &NodeGraph, colors: Vec<u32>, n: usize, state: &mut Search, root: bool) {
    if state.leaves >= LEAF_LIMIT {
        return;
    }
    // First non-singleton cell.
    let mut counts = vec![0u32; n];
    for &c in &colors {
        counts[c as usize] += 1;
    }
    let Some(cell) = (0..n).find(|&c| counts[c] > 1) else {
        state.leaves += 1;
        let (cert, sign) = graph.certificate(&colors);
        match &state.best {
            None => state.best = Some((cert, sign, colors)),
            Some((b, s, best_labels)) if cert == *b => {
                if *s != sign {
                    state.zero = true;
                }
                // An automorphism: node u and the node with u's label in the best leaf.
                let mut order = vec![0u32; n];
                for (v, &l) in best_labels.iter().enumerate() {
                    order[l as usize] = v as u32;
                }
                for u in 0..n {
                    let (a, b) = (orbit(&mut state.orbits, u as u32), orbit(&mut state.orbits, order[colors[u] as usize]));
                    if a != b {
                        state.orbits[a as usize] = b;
                    }
                }
            }
            Some((b, _, _)) if cert < *b => state.best = Some((cert, sign, colors)),
            _ => {}
        }
        return;
    };
    let mut explored: Vec<u32> = Vec::new();
    for v in 0..n {
        if colors[v] as usize != cell {
            continue;
        }
        // Children in the orbit of an explored child give the same leaves.
        if root && explored.iter().any(|&w| orbit(&mut state.orbits, w) == orbit(&mut state.orbits, v as u32)) {
            continue;
        }
        explored.push(v as u32);
        let mut next: Vec<u32> = colors
            .iter()
            .enumerate()
            .map(|(u, &c)| 2 * c + u32::from(c as usize == cell && u != v))
            .collect();
        // Re-rank.
        let mut sorted = next.clone();
        sorted.sort_unstable();
        sorted.dedup();
        for c in next.iter_mut() {
            *c = sorted.binary_search(c).unwrap() as u32;
        }
        graph.refine(&mut next);
        search(graph, next, n, state, false);
        if state.zero {
            return;
        }
    }
}

impl Net {
    /// The canonical form: (certificate, sign, canonical net).
    pub(crate) fn canonicalize(&self) -> (Vec<u32>, i8, Net) {
        let graph = NodeGraph::of(self);
        let canon = canonical(&graph);
        if canon.sign == 0 {
            return (canon.certificate, 0, Net::default());
        }
        let net = rebuild(self, &graph, &canon.labels);
        (canon.certificate, canon.sign, net)
    }
}

/// Rebuild the network in canonical node order with canonical index labels
/// and structure constants written in canonical slot order.
pub(crate) fn rebuild(net: &Net, graph: &NodeGraph, labels: &[u32]) -> Net {
    let n = graph.kind.len();
    let mut order = vec![0u32; n];
    for v in 0..n {
        order[labels[v] as usize] = v as u32;
    }
    // Slot order per node.
    let slot_order = |v: usize| -> Vec<usize> {
        match graph.kind[v] {
            K_TR | K_TA => vec![0],
            K_F => {
                let mut legs: Vec<(u32, Ix, usize)> = graph.adj[v]
                    .iter()
                    .enumerate()
                    .map(|(s, &(w, t))| ((labels[w as usize] << 2) | graph.role(w, t) as u32, graph.labels[v][s], s))
                    .collect();
                legs.sort_unstable();
                legs.into_iter().map(|l| l.2).collect()
            }
            k => {
                let skip = if k >= K_SR { 2 } else { 0 };
                let mut legs: Vec<(u32, usize)> = graph.adj[v]
                    .iter()
                    .enumerate()
                    .skip(skip)
                    .map(|(s, &(w, t))| ((labels[w as usize] << 2) | graph.role(w, t) as u32, s))
                    .collect();
                legs.sort_unstable();
                legs.into_iter().map(|l| l.1).collect()
            }
        }
    };
    // Assign index labels by first appearance.
    let mut edge_label: std::collections::HashMap<(u32, u8), Ix> = Default::default();
    let mut next_label: Ix = 0;
    let mut node_slots: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &v in &order {
        let v = v as usize;
        let slots = slot_order(v);
        for &s in &slots {
            let key = (v as u32, s as u8);
            if !edge_label.contains_key(&key) {
                let other = graph.adj[v][s];
                edge_label.insert(key, next_label);
                edge_label.insert(other, next_label);
                next_label += 1;
            }
        }
        node_slots[v] = slots;
    }
    let mut objs = Vec::new();
    let mut done_trace = vec![false; net.objs.len()];
    for &v in &order {
        let v = v as usize;
        let (o, _) = graph.origin[v];
        match &net.objs[o as usize] {
            Obj::F(_) => {
                let s = &node_slots[v];
                objs.push(Obj::F([
                    edge_label[&(v as u32, s[0] as u8)],
                    edge_label[&(v as u32, s[1] as u8)],
                    edge_label[&(v as u32, s[2] as u8)],
                ]));
            }
            Obj::D(rep, _) => {
                objs.push(Obj::D(*rep, node_slots[v].iter().map(|&s| edge_label[&(v as u32, s as u8)]).collect()));
            }
            Obj::St(rep, split, idx) => {
                if graph.kind[v] < K_SR {
                    continue;
                }
                let mut seq: Vec<Ix> = node_slots[v].iter().map(|&s| edge_label[&(v as u32, s as u8)]).collect();
                let mut w = graph.adj[v][1].0 as usize;
                for _ in *split..idx.len() {
                    seq.push(edge_label[&(w as u32, 0)]);
                    w = graph.adj[w][2].0 as usize;
                }
                objs.push(Obj::St(*rep, *split, seq));
            }
            Obj::Tr(rep, idx) => {
                if done_trace[o as usize] {
                    continue;
                }
                done_trace[o as usize] = true;
                // Start at this (lowest labelled) node and follow "out".
                let mut seq = Vec::with_capacity(idx.len());
                let mut w = v;
                for _ in 0..idx.len() {
                    seq.push(edge_label[&(w as u32, 0)]);
                    w = graph.adj[w][2].0 as usize;
                }
                objs.push(Obj::Tr(*rep, seq));
            }
        }
    }
    Net { objs, next: next_label }
}

pub(crate) enum SmallCut {
    Bridge,
    Two(u32, u32),
}

/// Edges of the node graph: (node, slot) pairs, adjoint or fundamental.
pub(crate) struct Edges {
    pub(crate) ends: Vec<[(u32, u8); 2]>,
    pub(crate) fundamental: Vec<bool>,
    /// node -> list of (edge, other node)
    pub(crate) incident: Vec<Vec<(u32, u32)>>,
}

impl Edges {
    pub(crate) fn of(graph: &NodeGraph) -> Edges {
        let n = graph.kind.len();
        let mut ends = Vec::new();
        let mut fundamental = Vec::new();
        let mut incident = vec![Vec::new(); n];
        for v in 0..n {
            for (s, &(w, t)) in graph.adj[v].iter().enumerate() {
                if (v as u32, s as u8) < (w, t) {
                    let e = ends.len() as u32;
                    ends.push([(v as u32, s as u8), (w, t)]);
                    let fund = (matches!(graph.kind[v], K_TR | K_TA) && s != 0)
                        || (graph.kind[v] >= K_SR && s < 2);
                    fundamental.push(fund);
                    incident[v].push((e, w));
                    incident[w as usize].push((e, v as u32));
                }
            }
        }
        Edges { ends, fundamental, incident }
    }

    /// Connected components of the nodes, ignoring the edge `skip`.
    pub(crate) fn components(&self, skip: &[u32]) -> (usize, Vec<u32>) {
        let n = self.incident.len();
        let mut comp = vec![u32::MAX; n];
        let mut count = 0u32;
        let mut stack = Vec::new();
        for s in 0..n {
            if comp[s] != u32::MAX {
                continue;
            }
            comp[s] = count;
            stack.push(s as u32);
            while let Some(v) = stack.pop() {
                for &(e, w) in &self.incident[v as usize] {
                    if skip.contains(&e) || comp[w as usize] != u32::MAX {
                        continue;
                    }
                    comp[w as usize] = count;
                    stack.push(w);
                }
            }
            count += 1;
        }
        (count as usize, comp)
    }

    /// A bridge, else the first pair of edges forming a 2-edge cut, of a
    /// connected graph: two edges form a cut iff their cycle-space cover
    /// sets coincide (a tree edge is covered by the non-tree edges whose
    /// fundamental cycle passes through it; a non-tree edge by itself).
    pub(crate) fn small_cut(&self) -> Option<SmallCut> {
        let n = self.incident.len();
        let m = self.ends.len();
        // DFS tree from node 0.
        let mut parent_edge = vec![u32::MAX; n];
        let mut order = Vec::with_capacity(n);
        let mut seen = vec![false; n];
        let mut stack = vec![0u32];
        seen[0] = true;
        let mut tree = vec![false; m];
        while let Some(v) = stack.pop() {
            order.push(v);
            for &(e, w) in &self.incident[v as usize] {
                if !seen[w as usize] {
                    seen[w as usize] = true;
                    parent_edge[w as usize] = e;
                    tree[e as usize] = true;
                    stack.push(w);
                }
            }
        }
        let non_tree: Vec<u32> = (0..m as u32).filter(|&e| !tree[e as usize]).collect();
        if non_tree.len() > 128 {
            return self.small_cut_slow();
        }
        let mut cover = vec![0u128; m];
        let mut val = vec![0u128; n];
        for (bit, &e) in non_tree.iter().enumerate() {
            let b = 1u128 << bit;
            cover[e as usize] = b;
            let [(u, _), (v, _)] = self.ends[e as usize];
            val[u as usize] ^= b;
            val[v as usize] ^= b;
        }
        // Subtree sums in reverse DFS order.
        for &v in order.iter().rev() {
            let e = parent_edge[v as usize];
            if e == u32::MAX {
                continue;
            }
            cover[e as usize] = val[v as usize];
            let [(a, _), (b, _)] = self.ends[e as usize];
            let p = if a == v { b } else { a };
            val[p as usize] ^= val[v as usize];
        }
        if (0..m).any(|e| tree[e] && cover[e] == 0) {
            return Some(SmallCut::Bridge);
        }
        for e1 in 0..m {
            for e2 in e1 + 1..m {
                if cover[e1] == cover[e2] {
                    return Some(SmallCut::Two(e1 as u32, e2 as u32));
                }
            }
        }
        None
    }

    fn small_cut_slow(&self) -> Option<SmallCut> {
        if !self.bridges(None).is_empty() {
            return Some(SmallCut::Bridge);
        }
        for e in 0..self.ends.len() as u32 {
            if let Some(&b) = self.bridges(Some(e)).first() {
                return Some(SmallCut::Two(e.min(b), e.max(b)));
            }
        }
        None
    }

    /// Bridges, ignoring edge `skip` (Tarjan, multigraph aware).
    pub(crate) fn bridges(&self, skip: Option<u32>) -> Vec<u32> {
        let n = self.incident.len();
        let mut disc = vec![u32::MAX; n];
        let mut low = vec![0u32; n];
        let mut time = 0u32;
        let mut out = Vec::new();
        // Iterative DFS: (node, parent edge, next incident position)
        for root in 0..n {
            if disc[root] != u32::MAX {
                continue;
            }
            let mut stack: Vec<(u32, u32, usize)> = vec![(root as u32, u32::MAX, 0)];
            disc[root] = time;
            low[root] = time;
            time += 1;
            while let Some(&mut (v, pe, ref mut pos)) = stack.last_mut() {
                if *pos < self.incident[v as usize].len() {
                    let (e, w) = self.incident[v as usize][*pos];
                    *pos += 1;
                    if Some(e) == skip || e == pe {
                        continue;
                    }
                    if disc[w as usize] == u32::MAX {
                        disc[w as usize] = time;
                        low[w as usize] = time;
                        time += 1;
                        stack.push((w, e, 0));
                    } else {
                        low[v as usize] = low[v as usize].min(disc[w as usize]);
                    }
                } else {
                    stack.pop();
                    if let Some(&(parent, _, _)) = stack.last() {
                        low[parent as usize] = low[parent as usize].min(low[v as usize]);
                        if low[v as usize] > disc[parent as usize] {
                            out.push(pe);
                        }
                    }
                }
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_form_respects_antisymmetry_and_relabelling() {
        // f(a,b,c) f(a,b,c) and f(b,a,c) f(a,b,c) differ by a sign.
        let x = Net::new(vec![Obj::F([0, 1, 2]), Obj::F([0, 1, 2])]);
        let y = Net::new(vec![Obj::F([1, 0, 2]), Obj::F([0, 1, 2])]);
        let (cx, sx, _) = x.canonicalize();
        let (cy, sy, _) = y.canonicalize();
        assert_eq!(cx, cy);
        assert_eq!(sx, -sy);
        // Relabelled traces are equal.
        let t1 = Net::new(vec![Obj::Tr(Rep::R, vec![0, 1, 2]), Obj::Tr(Rep::R, vec![0, 1, 2])]);
        let t2 = Net::new(vec![Obj::Tr(Rep::R, vec![5, 3, 4]), Obj::Tr(Rep::R, vec![4, 5, 3])]);
        let (c1, s1, _) = t1.canonicalize();
        let (c2, s2, _) = t2.canonicalize();
        assert_eq!((c1, s1), (c2, s2));
        // Tr(abc)Tr(acb) differs.
        let t3 = Net::new(vec![Obj::Tr(Rep::R, vec![0, 1, 2]), Obj::Tr(Rep::R, vec![0, 2, 1])]);
        assert_ne!(t3.canonicalize().0, t1.canonicalize().0);
    }

    #[test]
    fn odd_automorphism_is_detected() {
        // K_{3,3} of structure constants: exchanging two vertices of one side
        // reverses three of them (form_reference's six_f_k33 case).
        let (a, b, c, d, e, f, h, i, j) = (0, 1, 2, 3, 4, 5, 6, 7, 8);
        let x = Net::new(vec![
            Obj::F([a, b, c]),
            Obj::F([d, e, f]),
            Obj::F([h, i, j]),
            Obj::F([a, d, h]),
            Obj::F([b, e, i]),
            Obj::F([c, f, j]),
        ]);
        assert_eq!(x.canonicalize().1, 0);
        // A symmetric tensor network has no sign.
        let y = Net::new(vec![Obj::D(Rep::R, vec![0, 1, 2]), Obj::D(Rep::R, vec![0, 1, 2])]);
        assert_eq!(y.canonicalize().1, 1);
    }
}
