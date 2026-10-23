//! SPQR decomposition primitive of the constrained embedding, and a
//! decomposition of the small blocks of diagram expansions by recursive
//! separation-pair splits.
use super::embedding::{require, Graph, Id, Ids, Result, Set};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet, VecDeque};

/// One block in incidence order, in the JSON shape of the native `decompose`
/// request.
#[derive(Serialize)]
pub(crate) struct SpqrRequest<'a> {
    pub(crate) nodes: &'a Set,
    pub(crate) edges: Vec<SpqrEdge<'a>>,
}

#[derive(Serialize)]
pub(crate) struct SpqrEdge<'a> {
    pub(crate) id: &'a str,
    pub(crate) source: &'a str,
    pub(crate) target: &'a str,
}

impl<'a> SpqrRequest<'a> {
    pub(crate) fn new(g: &'a Graph) -> Self {
        Self {
            nodes: &g.nodes,
            edges: g
                .edges
                .iter()
                .map(|(id, [source, target])| SpqrEdge { id, source, target })
                .collect(),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Deserialize, Serialize)]
pub(crate) enum SpqrKind {
    S,
    P,
    R,
}

/// SPQR tree of a planar biconnected block. Skeleton edge identities are
/// unique within the decomposition; real edges name their original edge and
/// virtual edges their twin in the adjacent skeleton. Rotations list skeleton
/// edges clockwise and form a planar embedding of every skeleton.
#[derive(Clone, Debug, Deserialize)]
pub(crate) struct Decomposition {
    pub(crate) root: Id,
    pub(crate) components: Vec<SpqrComponent>,
}

#[derive(Clone, Debug, Deserialize)]
pub(crate) struct SpqrComponent {
    pub(crate) id: Id,
    #[serde(rename = "type")]
    pub(crate) kind: SpqrKind,
    pub(crate) nodes: Ids,
    pub(crate) edges: Vec<SkeletonEdge>,
    pub(crate) rotation_cw: BTreeMap<Id, Ids>,
}

#[derive(Clone, Debug, Deserialize)]
pub(crate) struct SkeletonEdge {
    pub(crate) id: Id,
    pub(crate) source: Id,
    pub(crate) target: Id,
    pub(crate) real_edge: Option<Id>,
    pub(crate) twin: Option<SkeletonTwin>,
}

#[derive(Clone, Debug, Deserialize)]
pub(crate) struct SkeletonTwin {
    pub(crate) component: Id,
    pub(crate) edge: Id,
}

pub(crate) trait SpqrDecomposer {
    /// SPQR tree of a biconnected block with at least three edges and no
    /// self-loops, or `None` when the block is not planar.
    fn decompose(&self, block: &SpqrRequest<'_>) -> Result<Option<Decomposition>>;
}

/// Triconnected components by repeated splitting at separation pairs, merged
/// maximally into bonds (P) and polygons (S); the remaining triconnected
/// skeletons (R) are embedded by path addition. Blocks of diagram expansions
/// are small, so each split simply tests every vertex `a` for an articulation
/// vertex `b` of the graph without `a`.
///
/// Ties are broken by original identities. A skeleton edge's rank is the
/// smallest request position among the real edges it represents; skeletons
/// list edges by rank, and each rotation starts at its lowest-rank edge. The
/// root holds the first request edge and tree nodes are numbered breadth-first
/// from it in rank order. Bonds list their edges clockwise by rank around the
/// smaller pole. Of the two mirror images of a triconnected skeleton, the one
/// is chosen in which, around the smallest vertex, the clockwise successor of
/// its smallest neighbor is smaller than that neighbor's predecessor.
pub(crate) struct SplitDecomposer;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Arc {
    Real(usize),
    Virtual(usize),
}

/// Edge of a split component; virtual twins share their pair index.
#[derive(Clone, Copy)]
struct SplitEdge {
    ends: [usize; 2],
    arc: Arc,
}

impl SplitEdge {
    fn pair(&self) -> (usize, usize) {
        let [a, b] = self.ends;
        (a.min(b), a.max(b))
    }
}

struct Split {
    kind: SpqrKind,
    edges: Vec<SplitEdge>,
}

/// Simple graph on the block's vertex indices with sorted neighbor lists;
/// vertices outside the component have none.
struct Adjacency(Vec<Vec<usize>>);

impl Adjacency {
    fn new(edges: &[SplitEdge], vertices: usize) -> Self {
        let mut out = vec![Vec::new(); vertices];
        for e in edges {
            let [a, b] = e.ends;
            out[a].push(b);
            out[b].push(a);
        }
        for neighbors in &mut out {
            neighbors.sort_unstable();
            neighbors.dedup();
        }
        Self(out)
    }

    fn vertices(&self) -> impl Iterator<Item = usize> + '_ {
        (0..self.0.len()).filter(|&v| !self.0[v].is_empty())
    }

    /// Vertices reachable from `start` without entering `blocked` ones.
    fn reachable(&self, start: usize, blocked: impl Fn(usize) -> bool) -> Vec<bool> {
        let mut seen = vec![false; self.0.len()];
        seen[start] = true;
        let mut pending = vec![start];
        while let Some(u) = pending.pop() {
            for &v in &self.0[u] {
                if !blocked(v) && !seen[v] {
                    seen[v] = true;
                    pending.push(v);
                }
            }
        }
        seen
    }

    /// Articulation vertices of the graph without `removed`.
    fn articulations(&self, removed: Option<usize>) -> Vec<bool> {
        let mut search = Lowpoints {
            adjacency: &self.0,
            removed,
            discovery: vec![usize::MAX; self.0.len()],
            low: vec![0; self.0.len()],
            cuts: vec![false; self.0.len()],
            visited: 0,
        };
        if let Some(root) = self.vertices().find(|&v| Some(v) != removed) {
            search.visit(root, None);
        }
        search.cuts
    }

    /// First separation pair `{a, b}` of a simple biconnected graph and the
    /// vertices of one component of the graph without `a` and `b`.
    fn separation_pair(&self) -> Option<(usize, usize, Vec<bool>)> {
        self.vertices().find_map(|a| {
            let b = self.articulations(Some(a)).iter().position(|&cut| cut)?;
            let start = self.vertices().find(|&v| v != a && v != b)?;
            Some((a, b, self.reachable(start, |v| v == a || v == b)))
        })
    }
}

/// Tarjan low points of a depth-first search.
struct Lowpoints<'a> {
    adjacency: &'a [Vec<usize>],
    removed: Option<usize>,
    discovery: Vec<usize>,
    low: Vec<usize>,
    cuts: Vec<bool>,
    visited: usize,
}

impl Lowpoints<'_> {
    fn visit(&mut self, u: usize, parent: Option<usize>) {
        self.discovery[u] = self.visited;
        self.low[u] = self.visited;
        self.visited += 1;
        let mut children = 0;
        for &v in &self.adjacency[u] {
            if Some(v) == self.removed || Some(v) == parent {
                continue;
            }
            if self.discovery[v] != usize::MAX {
                self.low[u] = self.low[u].min(self.discovery[v]);
                continue;
            }
            children += 1;
            self.visit(v, Some(u));
            self.low[u] = self.low[u].min(self.low[v]);
            if parent.is_some() && self.low[v] >= self.discovery[u] {
                self.cuts[u] = true;
            }
        }
        if parent.is_none() && children > 1 {
            self.cuts[u] = true;
        }
    }
}

/// Split components of a biconnected multigraph: bonds, polygons and
/// triconnected simple graphs, connected by numbered virtual edge pairs.
#[derive(Default)]
struct Splitter {
    vertices: usize,
    pairs: usize,
    done: Vec<Split>,
}

impl Splitter {
    fn virtual_pair(&mut self, a: usize, b: usize) -> SplitEdge {
        self.pairs += 1;
        SplitEdge {
            ends: [a.min(b), a.max(b)],
            arc: Arc::Virtual(self.pairs - 1),
        }
    }

    fn finish(&mut self, kind: SpqrKind, edges: Vec<SplitEdge>) {
        self.done.push(Split { kind, edges });
    }

    fn split(&mut self, edges: Vec<SplitEdge>) {
        let mut multiplicity: BTreeMap<(usize, usize), usize> = BTreeMap::new();
        for e in &edges {
            *multiplicity.entry(e.pair()).or_default() += 1;
        }
        if multiplicity.len() == 1 {
            return self.finish(SpqrKind::P, edges);
        }
        // Split off a bundle of parallel edges as a bond.
        if let Some((&(u, v), _)) = multiplicity.iter().find(|(_, m)| **m > 1) {
            let (mut bundle, mut rest): (Vec<_>, Vec<_>) =
                edges.into_iter().partition(|e| e.pair() == (u, v));
            let twin = self.virtual_pair(u, v);
            bundle.push(twin);
            rest.push(twin);
            self.finish(SpqrKind::P, bundle);
            return self.split(rest);
        }
        let adjacency = Adjacency::new(&edges, self.vertices);
        if adjacency.vertices().all(|v| adjacency.0[v].len() == 2) {
            return self.finish(SpqrKind::S, edges);
        }
        let Some((a, b, side)) = adjacency.separation_pair() else {
            return self.finish(SpqrKind::R, edges);
        };
        let (mut first, mut second): (Vec<_>, Vec<_>) = edges
            .into_iter()
            .partition(|e| e.ends.iter().any(|&v| side[v]));
        let twin = self.virtual_pair(a, b);
        first.push(twin);
        second.push(twin);
        self.split(first);
        self.split(second);
    }
}

/// Merge adjacent bonds and adjacent polygons along their shared virtual
/// pairs; returns the merged components and the virtual pairs that remain
/// tree edges.
fn merge(splits: Vec<Split>, pairs: usize) -> (Vec<Split>, Vec<bool>) {
    let mut halves = vec![Vec::new(); pairs];
    for (i, s) in splits.iter().enumerate() {
        for e in &s.edges {
            if let Arc::Virtual(k) = e.arc {
                halves[k].push(i);
            }
        }
    }
    let mut parent: Vec<usize> = (0..splits.len()).collect();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    let mut tree_edge = vec![true; pairs];
    for (k, halves) in halves.iter().enumerate() {
        let [x, y] = halves[..] else {
            continue;
        };
        if splits[x].kind == splits[y].kind && splits[x].kind != SpqrKind::R {
            let (rx, ry) = (find(&mut parent, x), find(&mut parent, y));
            parent[rx.max(ry)] = rx.min(ry);
            tree_edge[k] = false;
        }
    }
    let mut groups: BTreeMap<usize, Split> = BTreeMap::new();
    for (i, s) in splits.into_iter().enumerate() {
        let group = groups.entry(find(&mut parent, i)).or_insert(Split {
            kind: s.kind,
            edges: Vec::new(),
        });
        group.edges.extend(
            s.edges
                .into_iter()
                .filter(|e| !matches!(e.arc, Arc::Virtual(k) if !tree_edge[k])),
        );
    }
    (groups.into_values().collect(), tree_edge)
}

/// Clockwise neighbor orders of a planar embedding of a simple biconnected
/// graph by Demoucron–Malgrange–Pertuiset path addition, or `None` if the
/// graph is not planar. Faces are traced as in `Graph::faces`: arriving at `b`
/// from `a`, leave towards the neighbor preceding `a` around `b`.
fn planar_rotation(adjacency: &Adjacency) -> Option<BTreeMap<usize, Vec<usize>>> {
    let v0 = adjacency.vertices().next()?;
    let v1 = adjacency.0[v0][0];
    // Initial cycle: the edge v0-v1 closed by a shortest path from v1 to v0.
    let mut pred = vec![usize::MAX; adjacency.0.len()];
    pred[v1] = v1;
    let mut queue = VecDeque::from([v1]);
    while let Some(u) = queue.pop_front() {
        for &w in &adjacency.0[u] {
            if (u, w) != (v1, v0) && pred[w] == usize::MAX {
                pred[w] = u;
                queue.push_back(w);
            }
        }
    }
    let mut cycle = vec![v0];
    while cycle[cycle.len() - 1] != v1 {
        let previous = pred[cycle[cycle.len() - 1]];
        if previous == usize::MAX {
            return None;
        }
        cycle.push(previous);
    }
    let mut rotation: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    let mut embedded: BTreeSet<(usize, usize)> = BTreeSet::new();
    for (i, &c) in cycle.iter().enumerate() {
        let next = cycle[(i + 1) % cycle.len()];
        rotation.insert(c, vec![cycle[(i + cycle.len() - 1) % cycle.len()], next]);
        embedded.insert((c.min(next), c.max(next)));
    }
    loop {
        let faces = trace_faces(&rotation);
        let mut choice: Option<(Vec<usize>, usize, usize)> = None;
        for bridge in Bridge::all(adjacency, &rotation, &embedded) {
            let admissible: Vec<usize> = faces
                .iter()
                .enumerate()
                .filter(|(_, face)| {
                    bridge
                        .attachments
                        .iter()
                        .all(|a| face.iter().any(|(_, v)| v == a))
                })
                .map(|(i, _)| i)
                .collect();
            let face = *admissible.first()?;
            if choice
                .as_ref()
                .is_none_or(|(_, _, count)| *count > 1 && admissible.len() == 1)
            {
                choice = Some((bridge.path, face, admissible.len()));
            }
        }
        let Some((path, face, _)) = choice else {
            return Some(rotation);
        };
        // The face corner at `v` lies between the neighbor it leaves towards
        // and the one it arrives from; the new edge follows the former.
        let corner = |v: usize| {
            let darts = &faces[face];
            let i = darts.iter().position(|&(_, b)| b == v)?;
            Some(darts[(i + 1) % darts.len()].1)
        };
        let (x, y) = (path[0], path[path.len() - 1]);
        let (wx, wy) = (corner(x)?, corner(y)?);
        for (v, w, new) in [(x, wx, path[1]), (y, wy, path[path.len() - 2])] {
            let r = rotation.get_mut(&v)?;
            let i = r.iter().position(|&n| n == w)?;
            r.insert(i + 1, new);
        }
        for window in path.windows(3) {
            rotation.insert(window[1], vec![window[0], window[2]]);
        }
        for pair in path.windows(2) {
            embedded.insert((pair[0].min(pair[1]), pair[0].max(pair[1])));
        }
    }
}

/// Faces of a rotation system as cyclic dart sequences.
fn trace_faces(rotation: &BTreeMap<usize, Vec<usize>>) -> Vec<Vec<(usize, usize)>> {
    let mut visited = BTreeSet::new();
    let mut out = Vec::new();
    for (&u, r) in rotation {
        for &v in r {
            let mut face = Vec::new();
            let (mut a, mut b) = (u, v);
            while visited.insert((a, b)) {
                face.push((a, b));
                let around = &rotation[&b];
                let i = around.iter().position(|&n| n == a).unwrap_or_default();
                (a, b) = (b, around[(i + around.len() - 1) % around.len()]);
            }
            if !face.is_empty() {
                out.push(face);
            }
        }
    }
    out
}

/// Bridge of the embedded subgraph: its embedded attachment vertices and a
/// path through it between two of them.
struct Bridge {
    attachments: BTreeSet<usize>,
    path: Vec<usize>,
}

impl Bridge {
    /// Chords first, then the components of the remaining vertices, each with
    /// a path starting at its smallest attachment.
    fn all(
        adjacency: &Adjacency,
        rotation: &BTreeMap<usize, Vec<usize>>,
        embedded: &BTreeSet<(usize, usize)>,
    ) -> Vec<Self> {
        let mut out = Vec::new();
        for u in adjacency.vertices() {
            for &v in &adjacency.0[u] {
                if u < v
                    && rotation.contains_key(&u)
                    && rotation.contains_key(&v)
                    && !embedded.contains(&(u, v))
                {
                    out.push(Self {
                        attachments: BTreeSet::from([u, v]),
                        path: vec![u, v],
                    });
                }
            }
        }
        let mut seen = vec![false; adjacency.0.len()];
        for root in adjacency.vertices() {
            if rotation.contains_key(&root) || seen[root] {
                continue;
            }
            let component = adjacency.reachable(root, |v| rotation.contains_key(&v));
            let members: Vec<usize> = (0..component.len()).filter(|&v| component[v]).collect();
            for &v in &members {
                seen[v] = true;
            }
            let attachments: BTreeSet<usize> = members
                .iter()
                .flat_map(|&v| &adjacency.0[v])
                .filter(|v| rotation.contains_key(v))
                .copied()
                .collect();
            let Some(&x) = attachments.first() else {
                continue;
            };
            // Breadth-first from x through the component to another attachment.
            let mut pred = vec![usize::MAX; adjacency.0.len()];
            let mut queue = VecDeque::new();
            for &v in &adjacency.0[x] {
                if component[v] && pred[v] == usize::MAX {
                    pred[v] = x;
                    queue.push_back(v);
                }
            }
            'search: while let Some(z) = queue.pop_front() {
                for &w in &adjacency.0[z] {
                    if w != x && attachments.contains(&w) {
                        let mut path = vec![w, z];
                        while path[path.len() - 1] != x {
                            path.push(pred[path[path.len() - 1]]);
                        }
                        path.reverse();
                        out.push(Self { attachments, path });
                        break 'search;
                    }
                    if component[w] && pred[w] == usize::MAX {
                        pred[w] = z;
                        queue.push_back(w);
                    }
                }
            }
        }
        out
    }
}

impl SpqrDecomposer for SplitDecomposer {
    fn decompose(&self, block: &SpqrRequest<'_>) -> Result<Option<Decomposition>> {
        let names: Vec<&Id> = block.nodes.iter().collect();
        let index: BTreeMap<&str, usize> = names
            .iter()
            .enumerate()
            .map(|(i, n)| (n.as_str(), i))
            .collect();
        let mut seen: BTreeSet<&str> = BTreeSet::new();
        let mut edges = Vec::new();
        for (i, e) in block.edges.iter().enumerate() {
            require(seen.insert(e.id), &format!("duplicate edge ID: {}", e.id))?;
            let (Some(&a), Some(&b)) = (index.get(e.source), index.get(e.target)) else {
                return Err(format!("unknown endpoint for edge: {}", e.id));
            };
            require(
                a != b,
                &format!("self-loop must be separated before SPQR: {}", e.id),
            )?;
            edges.push(SplitEdge {
                ends: [a, b],
                arc: Arc::Real(i),
            });
        }
        require(
            names.len() >= 2 && edges.len() >= 3,
            "SPQR requires a block with at least three edges",
        )?;
        let simple = Adjacency::new(&edges, names.len());
        require(
            simple.vertices().count() == names.len()
                && simple.reachable(0, |_| false).iter().all(|&seen| seen)
                && !simple.articulations(None).contains(&true),
            "SPQR requires a biconnected block",
        )?;
        let mut splitter = Splitter {
            vertices: names.len(),
            ..Splitter::default()
        };
        splitter.split(edges);
        let (components, tree_edge) = merge(splitter.done, splitter.pairs);
        SpqrSkeletons::new(block, &names, components, &tree_edge).decomposition()
    }
}

/// Merged components with their tree structure and edge ranks.
struct SpqrSkeletons<'a> {
    block: &'a SpqrRequest<'a>,
    names: &'a [&'a Id],
    components: Vec<Split>,
    /// Component and edge position of both halves of each virtual tree edge.
    halves: BTreeMap<usize, [(usize, usize); 2]>,
}

impl<'a> SpqrSkeletons<'a> {
    fn new(
        block: &'a SpqrRequest<'a>,
        names: &'a [&'a Id],
        components: Vec<Split>,
        tree_edge: &[bool],
    ) -> Self {
        let mut halves: BTreeMap<usize, Vec<(usize, usize)>> = BTreeMap::new();
        for (c, s) in components.iter().enumerate() {
            for (i, e) in s.edges.iter().enumerate() {
                if let Arc::Virtual(k) = e.arc {
                    debug_assert!(tree_edge[k]);
                    halves.entry(k).or_default().push((c, i));
                }
            }
        }
        let halves = halves
            .into_iter()
            .filter_map(|(k, h)| Some((k, <[_; 2]>::try_from(h).ok()?)))
            .collect();
        Self {
            block,
            names,
            components,
            halves,
        }
    }

    /// Twin half of a virtual edge.
    fn twin(&self, k: usize, own: usize) -> (usize, usize) {
        let [a, b] = self.halves[&k];
        if a.0 == own {
            b
        } else {
            a
        }
    }

    /// Smallest request position of the real edges behind every edge,
    /// looking away from its component.
    fn ranks(&self) -> Vec<Vec<usize>> {
        let mut ranks: Vec<Vec<Option<usize>>> = self
            .components
            .iter()
            .map(|s| {
                s.edges
                    .iter()
                    .map(|e| match e.arc {
                        Arc::Real(i) => Some(i),
                        Arc::Virtual(_) => None,
                    })
                    .collect()
            })
            .collect();
        // Rank of a virtual half: the smallest real edge across its tree edge.
        fn behind(
            skeletons: &SpqrSkeletons<'_>,
            c: usize,
            i: usize,
            ranks: &mut Vec<Vec<Option<usize>>>,
        ) -> usize {
            if let Some(rank) = ranks[c][i] {
                return rank;
            }
            let Arc::Virtual(k) = skeletons.components[c].edges[i].arc else {
                unreachable!("real edges are ranked");
            };
            let (d, j) = skeletons.twin(k, c);
            let rank = (0..skeletons.components[d].edges.len())
                .filter(|&m| m != j)
                .map(|m| behind(skeletons, d, m, ranks))
                .min()
                .unwrap_or(usize::MAX);
            ranks[c][i] = Some(rank);
            rank
        }
        for c in 0..self.components.len() {
            for i in 0..self.components[c].edges.len() {
                behind(self, c, i, &mut ranks);
            }
        }
        ranks
            .into_iter()
            .map(|r| r.into_iter().map(Option::unwrap_or_default).collect())
            .collect()
    }

    fn decomposition(&self) -> Result<Option<Decomposition>> {
        let ranks = self.ranks();
        let orders: Vec<Vec<usize>> = ranks
            .iter()
            .map(|r| {
                let mut order: Vec<usize> = (0..r.len()).collect();
                order.sort_unstable_by_key(|&i| r[i]);
                order
            })
            .collect();
        // Breadth-first numbering from the component holding the first edge.
        let root = self
            .components
            .iter()
            .position(|s| s.edges.iter().any(|e| e.arc == Arc::Real(0)))
            .ok_or("SPQR lost the first edge")?;
        let mut numbering = vec![None; self.components.len()];
        let mut sequence = Vec::new();
        let mut queue = VecDeque::from([root]);
        numbering[root] = Some(0);
        while let Some(c) = queue.pop_front() {
            sequence.push(c);
            for &i in &orders[c] {
                if let Arc::Virtual(k) = self.components[c].edges[i].arc {
                    let (d, _) = self.twin(k, c);
                    if numbering[d].is_none() {
                        numbering[d] = Some(sequence.len() + queue.len());
                        queue.push_back(d);
                    }
                }
            }
        }
        let cid = |c: usize| numbering[c].map(|n| n.to_string()).unwrap_or_default();
        let local: Vec<Vec<usize>> = orders
            .iter()
            .map(|order| {
                let mut local = vec![0; order.len()];
                for (position, &i) in order.iter().enumerate() {
                    local[i] = position;
                }
                local
            })
            .collect();
        let edge_id = |c: usize, i: usize| format!("{}:{}", cid(c), local[c][i]);
        let mut components = Vec::new();
        for &c in &sequence {
            let split = &self.components[c];
            let Some(rotation) = self.rotation(split, &ranks[c]) else {
                return Ok(None);
            };
            let edges = orders[c]
                .iter()
                .map(|&i| {
                    let e = &split.edges[i];
                    let (source, target, real_edge, twin) = match e.arc {
                        Arc::Real(r) => {
                            let edge = &self.block.edges[r];
                            let id = Some(edge.id.to_owned());
                            (edge.source.to_owned(), edge.target.to_owned(), id, None)
                        }
                        Arc::Virtual(k) => {
                            let (d, j) = self.twin(k, c);
                            let twin = SkeletonTwin {
                                component: cid(d),
                                edge: edge_id(d, j),
                            };
                            let [a, b] = e.ends.map(|v| self.names[v].clone());
                            (a, b, None, Some(twin))
                        }
                    };
                    SkeletonEdge {
                        id: edge_id(c, i),
                        source,
                        target,
                        real_edge,
                        twin,
                    }
                })
                .collect();
            let nodes: BTreeSet<usize> = split.edges.iter().flat_map(|e| e.ends).collect();
            components.push(SpqrComponent {
                id: cid(c),
                kind: split.kind,
                nodes: nodes.iter().map(|&v| self.names[v].clone()).collect(),
                edges,
                rotation_cw: rotation
                    .into_iter()
                    .map(|(v, r)| {
                        (
                            self.names[v].clone(),
                            r.into_iter().map(|i| edge_id(c, i)).collect(),
                        )
                    })
                    .collect(),
            });
        }
        Ok(Some(Decomposition {
            root: cid(root),
            components,
        }))
    }

    /// Clockwise edge positions around every vertex of a skeleton, each
    /// starting at its lowest-rank edge; `None` for a nonplanar skeleton.
    fn rotation(&self, split: &Split, ranks: &[usize]) -> Option<BTreeMap<usize, Vec<usize>>> {
        let mut incident: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for (i, e) in split.edges.iter().enumerate() {
            for v in e.ends {
                incident.entry(v).or_default().push(i);
            }
        }
        let mut rotation = match split.kind {
            SpqrKind::S => incident,
            SpqrKind::P => {
                let mut poles = incident.into_iter();
                let (low, mut order) = poles.next()?;
                let (high, _) = poles.next()?;
                order.sort_unstable_by_key(|&i| ranks[i]);
                let reversed = order.iter().rev().copied().collect();
                BTreeMap::from([(low, order), (high, reversed)])
            }
            SpqrKind::R => {
                let mut neighbors =
                    planar_rotation(&Adjacency::new(&split.edges, self.names.len()))?;
                let (_, around) = neighbors.iter().next()?;
                let first = around.iter().position(|n| Some(n) == around.iter().min())?;
                let d = around.len();
                if around[(first + 1) % d] > around[(first + d - 1) % d] {
                    for r in neighbors.values_mut() {
                        r.reverse();
                    }
                }
                let edge = |u: usize, v: usize| {
                    split
                        .edges
                        .iter()
                        .position(|e| e.pair() == (u.min(v), u.max(v)))
                };
                neighbors
                    .into_iter()
                    .map(|(u, r)| {
                        Some((
                            u,
                            r.into_iter()
                                .map(|v| edge(u, v))
                                .collect::<Option<Vec<_>>>()?,
                        ))
                    })
                    .collect::<Option<_>>()?
            }
        };
        for r in rotation.values_mut() {
            let start = (0..r.len()).min_by_key(|&i| ranks[r[i]])?;
            r.rotate_left(start);
        }
        Some(rotation)
    }
}
