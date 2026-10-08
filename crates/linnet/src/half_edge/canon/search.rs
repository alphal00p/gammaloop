//! Individualization-refinement over the node partition (McKay's search tree).
//!
//! Cells are runs of `lab`; a cell is named by its start position. Refinement is a splitter
//! queue: for a splitter cell `W`, every node `v` adjacent to `W` gets the sum of weights of the
//! edge ends `(colour at v, colour at the W side)` joining it to `W`, and touched cells split by
//! that sum. Sums are commutative, so the input order of hedges never matters, and touched cells
//! and new fragments are processed in position order, so refinement is equivariant. Hash sums
//! may collide; that only weakens refinement, because leaves are compared by exact certificates.
//!
//! Only nodes are individualized, those of a target cell chosen like nauty's `bestcell` (see
//! [`Search::target_cell`]). A leaf is a partition whose non-singleton cells consist of
//! isolated nodes; its certificate lists, per node in leaf order, the degree and the sorted
//! `(colour, opposite node, opposite colour)` triples of its hedges. Twin hedges (automorphisms
//! fixing every node) and permutations of isolated nodes are handled analytically.

use std::cmp::Reverse;

use super::{
    flat::{Flat, NONE},
    Orbits,
};

#[inline]
pub(crate) fn mix(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^ (x >> 31)
}

fn pair(a: u32, b: u32) -> u64 {
    ((a as u64) << 32) | b as u64
}

/// Ordered partition of the nodes.
#[derive(Clone, Debug, Default)]
struct Part {
    lab: Vec<u32>,
    cell_of: Vec<u32>,
    /// End of the cell starting at a position; only meaningful at cell starts.
    cell_end: Vec<u32>,
}

impl Part {
    fn copy_from(&mut self, other: &Part) {
        self.lab.clone_from(&other.lab);
        self.cell_of.clone_from(&other.cell_of);
        self.cell_end.clone_from(&other.cell_end);
    }
}

#[derive(Debug, Default)]
struct Leaf {
    invariants: Vec<u64>,
    lab: Vec<u32>,
    cert: Vec<u32>,
    path: Vec<u32>,
}

impl Leaf {
    fn set(&mut self, invariants: &[u64], lab: &[u32], cert: &[u32], path: &[u32]) {
        assign(&mut self.invariants, invariants);
        assign(&mut self.lab, lab);
        assign(&mut self.cert, cert);
        assign(&mut self.path, path);
    }
}

fn assign<T: Copy>(dst: &mut Vec<T>, src: &[T]) {
    dst.clear();
    dst.extend_from_slice(src);
}

/// Per-level scratch: children of the node, explored children and an orbit table.
type LevelScratch = (Vec<u32>, Vec<u32>, Vec<u32>);

/// Reusable search state. After [`Search::run`], `best_lab` is the canonical node order,
/// `best_cert` the leaf certificate, `generators` the node automorphisms found (stride `n`) and
/// `group_order` the order of the induced node group without isolated nodes (`None` on overflow).
#[derive(Debug, Default)]
pub(crate) struct Search {
    n: usize,
    /// Neighbours across paired hedges, grouped like `Flat::adj`, with the weight of the edge end
    /// as seen from the neighbour.
    nbr_off: Vec<u32>,
    nbr: Vec<u32>,
    nbr_weight: Vec<u64>,
    parts: Vec<Part>,
    queue: Vec<u32>,
    in_queue: Vec<bool>,
    key: Vec<u64>,
    /// `node_stamp[v] == stamp` marks the nodes touched by the current splitter.
    node_stamp: Vec<u32>,
    stamp: u32,
    /// Per cell start, a count during one pass; zero between passes.
    cell_count: Vec<u32>,
    touched_cells: Vec<u32>,
    /// Keys and touched nodes of the cell being split.
    sorted: Vec<(u64, u32)>,
    /// Starts of the cells the target cell is chosen from.
    candidates: Vec<u32>,
    path: Vec<u32>,
    path_invariants: Vec<u64>,
    first: Leaf,
    best: Leaf,
    has_first: bool,
    pub(crate) generators: Vec<u32>,
    orbits: Orbits,
    canonical_index: Vec<u32>,
    cert: Vec<u32>,
    triples: Vec<[u32; 3]>,
    scratch: Vec<LevelScratch>,
    pub(crate) group_order: Option<u128>,
}

impl Search {
    pub(crate) fn best_lab(&self) -> &[u32] {
        &self.best.lab
    }

    pub(crate) fn best_cert(&self) -> &[u32] {
        &self.best.cert
    }

    pub(crate) fn n_generators(&self) -> usize {
        self.generators.len().checked_div(self.n).unwrap_or(0)
    }

    pub(crate) fn generator(&self, i: usize) -> &[u32] {
        &self.generators[i * self.n..(i + 1) * self.n]
    }

    /// Canonical node order, leaf certificate and node automorphisms of `g`.
    pub(crate) fn run(&mut self, g: &Flat) {
        let n = g.n_nodes();
        self.n = n;
        self.generators.clear();
        self.group_order = Some(1);
        self.has_first = false;
        self.first.set(&[], &[], &[], &[]);
        self.best.set(&[], &[], &[], &[]);
        self.path.clear();
        self.path_invariants.clear();
        self.queue.clear();
        self.in_queue.clear();
        self.in_queue.resize(n, false);
        self.key.clear();
        self.key.resize(n, 0);
        self.node_stamp.clear();
        self.node_stamp.resize(n, 0);
        self.stamp = 0;
        self.cell_count.clear();
        self.cell_count.resize(n, 0);
        self.canonical_index.resize(n, 0);
        if self.parts.is_empty() {
            self.parts.push(Part::default());
        }

        self.nbr_off.clear();
        self.nbr.clear();
        self.nbr_weight.clear();
        self.nbr_off.push(0);
        for v in 0..n {
            for &h in g.hedges_of(v) {
                let o = g.inv[h as usize];
                if o != h {
                    self.nbr.push(g.node_of[o as usize]);
                    self.nbr_weight
                        .push(mix(pair(g.hcol[o as usize], g.hcol[h as usize])) | 1);
                }
            }
            self.nbr_off.push(self.nbr.len() as u32);
        }

        self.initial_partition(g);
        let root = &self.parts[0];
        let mut s = 0;
        while s < n {
            self.queue.push(s as u32);
            self.in_queue[s] = true;
            s = root.cell_end[s] as usize;
        }
        self.refine(0, 0);
        if n > 0 {
            self.search(g, 0, true);
        }
    }

    /// Cells by node colour, then by the multiset of incident hedge colours (dangling legs only
    /// enter here); isolated nodes form their own cells.
    fn initial_partition(&mut self, g: &Flat) {
        let n = g.n_nodes();
        for v in 0..n {
            let hedges = g.hedges_of(v);
            let mut profile = 0u64;
            for &h in hedges {
                let o = g.inv[h as usize];
                let (opposite, self_loop) = if o == h {
                    (NONE, 0)
                } else {
                    (
                        g.hcol[o as usize],
                        (g.node_of[o as usize] as usize == v) as u64,
                    )
                };
                profile = profile.wrapping_add(mix(pair(g.hcol[h as usize], opposite) ^ self_loop));
            }
            self.key[v] = if hedges.is_empty() { 0 } else { profile | 1 };
        }
        let p = &mut self.parts[0];
        p.lab.clear();
        p.lab.extend(0..n as u32);
        p.cell_of.resize(n, 0);
        p.cell_end.resize(n, 0);
        let key = &self.key;
        p.lab
            .sort_unstable_by_key(|&v| (g.ncol[v as usize], key[v as usize]));
        let mut start = 0;
        for q in 1..=n {
            let split = q == n || {
                let (a, b) = (p.lab[q - 1] as usize, p.lab[q] as usize);
                (g.ncol[a], key[a]) != (g.ncol[b], key[b])
            };
            if split {
                p.cell_end[start] = q as u32;
                for &v in &p.lab[start..q] {
                    p.cell_of[v as usize] = start as u32;
                }
                start = q;
            }
        }
    }

    /// Refine `parts[level]` with the queued splitters; returns an invariant trace.
    ///
    /// Only the cells a splitter touches are visited, and only their touched nodes are sorted:
    /// the untouched nodes of a cell keep their order and become its first fragment.
    fn refine(&mut self, level: usize, seed: u64) -> u64 {
        let Search {
            nbr_off,
            nbr,
            nbr_weight,
            parts,
            queue,
            in_queue,
            key,
            node_stamp,
            stamp,
            cell_count,
            touched_cells,
            sorted,
            ..
        } = self;
        let Part {
            lab,
            cell_of,
            cell_end,
        } = &mut parts[level];
        // slices keep their bounds in registers across the stores below
        let (lab, cell_of, cell_end) = (&mut lab[..], &mut cell_of[..], &mut cell_end[..]);
        let (key, node_stamp, cell_count) =
            (&mut key[..], &mut node_stamp[..], &mut cell_count[..]);
        let mut trace = mix(seed);
        let mut qi = 0;
        while let Some(&w) = queue.get(qi) {
            let w = w as usize;
            qi += 1;
            in_queue[w] = false;
            if *stamp == u32::MAX {
                node_stamp.fill(0);
                *stamp = 0;
            }
            *stamp += 1;
            let st = *stamp;
            touched_cells.clear();
            for &u in &lab[w..cell_end[w] as usize] {
                let range = nbr_off[u as usize] as usize..nbr_off[u as usize + 1] as usize;
                for (&v, &weight) in nbr[range.clone()].iter().zip(&nbr_weight[range]) {
                    let v = v as usize;
                    if node_stamp[v] == st {
                        key[v] = key[v].wrapping_add(weight);
                        continue;
                    }
                    node_stamp[v] = st;
                    key[v] = weight;
                    let c = cell_of[v] as usize;
                    if cell_end[c] as usize - c > 1 {
                        if cell_count[c] == 0 {
                            touched_cells.push(c as u32);
                        }
                        cell_count[c] += 1;
                    }
                }
            }
            touched_cells.sort_unstable();
            for &c in touched_cells.iter() {
                let c = c as usize;
                let e = cell_end[c] as usize;
                let cell = &mut lab[c..e];
                let untouched = cell.len() - std::mem::take(&mut cell_count[c]) as usize;
                sorted.clear();
                if untouched == 0 {
                    let k0 = key[cell[0] as usize];
                    if cell.iter().all(|&v| key[v as usize] == k0) {
                        continue;
                    }
                    sorted.extend(cell.iter().map(|&v| (key[v as usize], v)));
                } else {
                    let mut q = 0;
                    for i in 0..cell.len() {
                        let v = cell[i];
                        if node_stamp[v as usize] == st {
                            sorted.push((key[v as usize], v));
                        } else {
                            cell[q] = v;
                            q += 1;
                        }
                    }
                }
                sorted.sort_unstable_by_key(|&(k, _)| k);
                for (slot, &(_, v)) in cell[untouched..].iter_mut().zip(sorted.iter()) {
                    *slot = v;
                }
                // fragments: the untouched nodes, then the touched ones by key
                let mut largest = (0, c);
                let mut fs = c;
                let untouched_run = (untouched > 0).then_some((untouched, 0));
                let runs = sorted
                    .chunk_by(|a, b| a.0 == b.0)
                    .map(|r| (r.len(), r[0].0));
                for (len, k) in untouched_run.into_iter().chain(runs) {
                    let q = fs + len;
                    cell_end[fs] = q as u32;
                    for &v in &lab[fs..q] {
                        cell_of[v as usize] = fs as u32;
                    }
                    trace = mix(trace ^ pair(c as u32, fs as u32) ^ k);
                    if len > largest.0 {
                        largest = (len, fs);
                    }
                    fs = q;
                }
                // Hopcroft: a split parent that is still queued covers its first fragment;
                // otherwise the counts into the largest fragment are implied by the others.
                let parent_queued = in_queue[c];
                let mut fs = c;
                while fs < e {
                    let skip = if parent_queued {
                        fs == c
                    } else {
                        fs == largest.1
                    };
                    if !skip {
                        in_queue[fs] = true;
                        queue.push(fs as u32);
                    }
                    fs = cell_end[fs] as usize;
                }
            }
        }
        queue.clear();
        trace
    }

    /// How many cells [`Self::target_cell`] scores. Scoring walks the edges of every candidate;
    /// scoring all non-singleton cells made each search node cost `O(n + m)` and deep searches
    /// through many small cells (twin pairs, unions of paths) up to 2.7 times slower at the
    /// same node count. Four keep the node counts of scoring all cells on projective planes and
    /// Latin squares.
    const CANDIDATES: usize = 4;

    /// The cell whose nodes the search individualizes. Among the first [`Self::CANDIDATES`]
    /// non-singleton cells that do not consist of isolated nodes, it is joined non-uniformly to
    /// the most other non-singleton cells, like nauty's `bestcell`, so leaves are reached sooner
    /// and differ earlier than under the first such cell. Ties go to the smallest cell, which has
    /// the fewest children, then to the first.
    ///
    /// A cell `d` counts as non-uniformly joined to `c` when the edge ends between them number
    /// neither zero nor `|c| |d|`; in an equitable partition, individualizing a node of `c` then
    /// splits `d`. Edge colours and multiplicities are not looked at, so a uniform count can
    /// hide such a split on coloured or multigraphs; telling those apart doubles the scan and
    /// cost up to 2% on Feynman graphs without saving search nodes.
    ///
    /// The choice is an isomorphism invariant of the partition, as the canonical form requires:
    /// an isomorphism maps each cell to the cell at the same position of the image partition,
    /// and the candidates and their scores depend only on cell positions, cell sizes and the
    /// number of edge ends between whole cells, so isomorphic graphs pick cells at the same
    /// position and their search trees stay isomorphic. Counting over whole cells rather than at
    /// one node, as nauty does, keeps this true when a hash collision in [`Self::refine`] leaves
    /// the partition short of equitable.
    fn target_cell(&mut self, g: &Flat, level: usize) -> Option<(usize, usize)> {
        let Search {
            n,
            nbr_off,
            nbr,
            parts,
            cell_count,
            candidates,
            touched_cells,
            ..
        } = self;
        let p = &parts[level];
        let size = |c: usize| p.cell_end[c] as usize - c;
        candidates.clear();
        let mut s = 0;
        while s < *n && candidates.len() < Self::CANDIDATES {
            if size(s) > 1 && !g.hedges_of(p.lab[s] as usize).is_empty() {
                candidates.push(s as u32);
            }
            s = p.cell_end[s] as usize;
        }
        let target = match candidates[..] {
            [] => return None,
            [c] => c as usize,
            _ => candidates
                .iter()
                .map(|&c| c as usize)
                .max_by_key(|&c| {
                    // edge ends from `c` into every other non-singleton cell
                    touched_cells.clear();
                    for &u in &p.lab[c..p.cell_end[c] as usize] {
                        let range = nbr_off[u as usize] as usize..nbr_off[u as usize + 1] as usize;
                        for &w in &nbr[range] {
                            let d = p.cell_of[w as usize] as usize;
                            if d != c && size(d) > 1 {
                                if cell_count[d] == 0 {
                                    touched_cells.push(d as u32);
                                }
                                cell_count[d] += 1;
                            }
                        }
                    }
                    let mut joined = 0;
                    for &d in touched_cells.iter() {
                        let d = d as usize;
                        let ends = std::mem::take(&mut cell_count[d]) as usize;
                        joined += usize::from(ends < size(c) * size(d));
                    }
                    (joined, Reverse(size(c)), Reverse(c))
                })
                .expect("several candidates"),
        };
        Some((target, p.cell_end[target] as usize))
    }

    /// `parts[level + 1]` = `parts[level]` with node `v` individualized; returns the seed of the
    /// refinement trace.
    fn individualize(&mut self, level: usize, v: u32) -> u64 {
        if self.parts.len() <= level + 1 {
            self.parts.push(Part::default());
        }
        let (done, next) = self.parts.split_at_mut(level + 1);
        let p = &mut next[0];
        p.copy_from(&done[level]);
        let s = p.cell_of[v as usize] as usize;
        let e = p.cell_end[s] as usize;
        let i = s + p.lab[s..e]
            .iter()
            .position(|&x| x == v)
            .expect("node in its cell");
        p.lab.swap(s, i);
        p.cell_end[s] = s as u32 + 1;
        p.cell_end[s + 1] = e as u32;
        for &x in &p.lab[s + 1..e] {
            p.cell_of[x as usize] = s as u32 + 1;
        }
        self.queue.push(s as u32);
        self.in_queue[s] = true;
        pair(level as u32, s as u32)
    }

    /// Node orbits under the generators fixing the first `level` nodes of the current path (or
    /// of the first path).
    fn orbits(&mut self, level: usize, on_first_path: bool, orbit: &mut Vec<u32>) {
        let n = self.n;
        let orbits = &mut self.orbits;
        orbits.reset(n);
        let prefix = if on_first_path {
            &self.first.path[..level]
        } else {
            &self.path[..level]
        };
        for gen in self.generators.chunks_exact(n) {
            if prefix.iter().any(|&v| gen[v as usize] != v) {
                continue;
            }
            orbits.join(gen);
        }
        orbit.clear();
        orbit.extend((0..n as u32).map(|v| orbits.find(v)));
    }

    /// Explore the subtree at `level`; returns the level to backjump to after an automorphism,
    /// or `usize::MAX`.
    fn search(&mut self, g: &Flat, level: usize, on_first_path: bool) -> usize {
        let Some((s, e)) = self.target_cell(g, level) else {
            return self.leaf(g, level);
        };
        if self.scratch.len() <= level {
            self.scratch.resize_with(level + 1, Default::default);
        }
        let (mut children, mut explored, mut orbit) = std::mem::take(&mut self.scratch[level]);
        children.clear();
        children.extend_from_slice(&self.parts[level].lab[s..e]);
        explored.clear();
        let jump = self.search_children(
            g,
            level,
            on_first_path,
            &children,
            &mut explored,
            &mut orbit,
        );
        self.scratch[level] = (children, explored, orbit);
        jump
    }

    fn search_children(
        &mut self,
        g: &Flat,
        level: usize,
        on_first_path: bool,
        children: &[u32],
        explored: &mut Vec<u32>,
        orbit: &mut Vec<u32>,
    ) -> usize {
        let mut orbit_generators = usize::MAX;
        for &v in children {
            if !explored.is_empty() {
                if orbit_generators != self.generators.len() {
                    self.orbits(level, false, orbit);
                    orbit_generators = self.generators.len();
                }
                let r = orbit[v as usize];
                if explored.iter().any(|&x| orbit[x as usize] == r) {
                    continue;
                }
            }
            explored.push(v);
            let seed = self.individualize(level, v);
            let invariant = self.refine(level + 1, seed);
            self.path.truncate(level);
            self.path.push(v);
            self.path_invariants.truncate(level);
            self.path_invariants.push(invariant);
            if self.has_first {
                let current = &self.path_invariants[..=level];
                let first = &self.first.invariants;
                let equals_first = first.len() > level && first[..=level] == *current;
                let best = &self.best.invariants;
                let m = best.len().min(level + 1);
                if !equals_first && current[..m] < best[..m] {
                    continue;
                }
            }
            let jump = self.search(g, level + 1, on_first_path && !self.has_first);
            if jump < level {
                return jump;
            }
        }
        if on_first_path {
            // orbit-stabilizer along the first path
            self.orbits(level, true, orbit);
            let r = orbit[self.first.path[level] as usize];
            let size = orbit.iter().filter(|&&x| x == r).count() as u128;
            self.group_order = self.group_order.and_then(|o| o.checked_mul(size));
        }
        usize::MAX
    }

    fn leaf(&mut self, g: &Flat, level: usize) -> usize {
        let lab = &self.parts[level].lab;
        for (q, &v) in lab.iter().enumerate() {
            self.canonical_index[v as usize] = q as u32;
        }
        self.cert.clear();
        for &v in lab {
            let hedges = g.hedges_of(v as usize);
            self.cert.push(hedges.len() as u32);
            self.triples.clear();
            for &h in hedges {
                let o = g.inv[h as usize];
                self.triples.push(if o == h {
                    [g.hcol[h as usize], NONE, NONE]
                } else {
                    [
                        g.hcol[h as usize],
                        self.canonical_index[g.node_of[o as usize] as usize],
                        g.hcol[o as usize],
                    ]
                });
            }
            self.triples.sort_unstable();
            self.cert.extend(self.triples.iter().flatten());
        }
        let invariants = &self.path_invariants[..level];
        let path = &self.path[..level];
        if !self.has_first {
            self.has_first = true;
            self.first.set(invariants, lab, &self.cert, path);
            self.best.set(invariants, lab, &self.cert, path);
            return usize::MAX;
        }
        let reference = if self.first.invariants == invariants && self.first.cert == self.cert {
            &self.first
        } else {
            match (invariants, &self.cert[..])
                .cmp(&(&self.best.invariants[..], &self.best.cert[..]))
            {
                std::cmp::Ordering::Equal => &self.best,
                std::cmp::Ordering::Greater => {
                    self.best.set(invariants, lab, &self.cert, path);
                    return usize::MAX;
                }
                std::cmp::Ordering::Less => return usize::MAX,
            }
        };
        // automorphism: reference leaf -> this leaf
        let start = self.generators.len();
        self.generators.resize(start + self.n, 0);
        for (&x, &y) in reference.lab.iter().zip(lab) {
            self.generators[start + x as usize] = y;
        }
        path.iter()
            .zip(&reference.path)
            .take_while(|(a, b)| a == b)
            .count()
    }
}
