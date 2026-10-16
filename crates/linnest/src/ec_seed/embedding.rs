//! Gutwenger–Klein–Mutzel EC expansion and two-sided SPQR insertion.
//! Edge-identity rotations are clockwise. Container order matches the reference
//! implementation: ties in shortest paths use the original incidence order.
use super::spqr::{Decomposition, SpqrDecomposer, SpqrKind, SpqrRequest};
use indexmap::IndexMap;
use serde::Serialize;
use std::collections::{hash_map::Entry, BTreeMap, BTreeSet, HashMap, HashSet, VecDeque};

pub(crate) type Id = String;
pub(crate) type Ids = Vec<Id>;
pub(crate) type Set = BTreeSet<Id>;
pub(crate) type Ends = [Id; 2];
type Dart = (Id, Id);
pub(crate) type Result<T> = std::result::Result<T, String>;

pub(crate) fn require(valid: bool, message: &str) -> Result<()> {
    if valid {
        Ok(())
    } else {
        Err(message.to_owned())
    }
}

/// Python dictionaries preserve insertion order, including for repeated edge
/// insertions. A sorted map here would silently change the certified seed.
/// Overwriting keeps a key's position; erasing shifts later keys forward.
#[derive(Clone, Debug, PartialEq, Serialize)]
#[serde(transparent)]
pub(crate) struct OrderedMap<T>(IndexMap<Id, T>);

impl<T> Default for OrderedMap<T> {
    fn default() -> Self {
        Self(IndexMap::new())
    }
}

impl<T> FromIterator<(Id, T)> for OrderedMap<T> {
    fn from_iter<I: IntoIterator<Item = (Id, T)>>(iter: I) -> Self {
        Self(iter.into_iter().collect())
    }
}

impl<T> OrderedMap<T> {
    pub(crate) fn get(&self, id: &str) -> Option<&T> {
        self.0.get(id)
    }

    pub(crate) fn at(&self, id: &str) -> Result<&T> {
        self.0
            .get(id)
            .ok_or_else(|| format!("Unknown identity: {id}"))
    }

    pub(crate) fn at_mut(&mut self, id: &str) -> Result<&mut T> {
        self.0
            .get_mut(id)
            .ok_or_else(|| format!("Unknown identity: {id}"))
    }

    pub(crate) fn contains(&self, id: &str) -> bool {
        self.0.contains_key(id)
    }

    pub(crate) fn insert(&mut self, id: Id, value: T) {
        self.0.insert(id, value);
    }

    pub(crate) fn erase(&mut self, id: &str) {
        self.0.shift_remove(id);
    }

    pub(crate) fn len(&self) -> usize {
        self.0.len()
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub(crate) fn iter(&self) -> indexmap::map::Iter<'_, Id, T> {
        self.0.iter()
    }

    pub(crate) fn iter_mut(&mut self) -> indexmap::map::IterMut<'_, Id, T> {
        self.0.iter_mut()
    }

    fn keys(&self) -> indexmap::map::Keys<'_, Id, T> {
        self.0.keys()
    }

    pub(crate) fn values(&self) -> indexmap::map::Values<'_, Id, T> {
        self.0.values()
    }

    pub(crate) fn key_set(&self) -> Set {
        self.0.keys().cloned().collect()
    }
}

impl<T: Clone> OrderedMap<T> {
    fn update(&mut self, other: &Self) {
        for (k, v) in other.iter() {
            self.insert(k.clone(), v.clone());
        }
    }
}

fn intersection(a: &Set, b: &Set) -> Set {
    a.intersection(b).cloned().collect()
}

pub(crate) fn as_set(v: &[Id]) -> Set {
    v.iter().cloned().collect()
}

fn index_of(v: &[Id], k: &str) -> Result<usize> {
    v.iter()
        .position(|x| x == k)
        .ok_or_else(|| format!("Missing incidence {k}"))
}

fn cyclic_equal(a: &[Id], b: &[Id]) -> bool {
    if a.len() != b.len() {
        return false;
    }
    if a.is_empty() {
        return true;
    }
    let Some(j) = b.iter().position(|x| *x == a[0]) else {
        return false;
    };
    (0..a.len()).all(|i| a[i] == b[(i + j) % b.len()])
}

fn after(v: &[Id], k: &str) -> Result<Ids> {
    let i = index_of(v, k)?;
    Ok([&v[i + 1..], &v[..i]].concat())
}

pub(crate) fn replace(v: &mut Ids, old: &str, replacement: &[Id]) -> Result<()> {
    let i = index_of(v, old)?;
    v.splice(i..=i, replacement.iter().cloned());
    Ok(())
}

/// Unordered pair of SPQR-tree nodes.
fn unordered(a: &str, b: &str) -> (Id, Id) {
    if a <= b {
        (a.to_owned(), b.to_owned())
    } else {
        (b.to_owned(), a.to_owned())
    }
}

/// Multigraph with clockwise edge rotations. Protected edges cannot be
/// crossed by insertions; oriented `hubs` fix the cyclic order of their
/// spokes, and faces at `wheel_hubs` must be protected triangles.
#[derive(Clone, Debug, Default)]
pub(crate) struct Graph {
    pub(crate) nodes: Set,
    pub(crate) protected_edges: Set,
    pub(crate) wheel_hubs: Set,
    pub(crate) edges: OrderedMap<Ends>,
    pub(crate) rotation: BTreeMap<Id, Ids>,
    pub(crate) hubs: BTreeMap<Id, Ids>,
}

/// Faces of a rotation system. A dart is an edge position with the side of
/// the endpoint it leaves from.
struct Faces<'a> {
    graph: &'a Graph,
    boundaries: Vec<Vec<(usize, usize)>>,
    face_of: Vec<[usize; 2]>,
}

impl<'a> Faces<'a> {
    fn len(&self) -> usize {
        self.boundaries.len()
    }

    /// Face of the dart of edge `e` at its endpoint `n`.
    fn of(&self, e: &str, n: &str) -> Result<usize> {
        self.graph
            .edges
            .position(e)
            .and_then(|i| Some(self.face_of[i][self.graph.edges.side(i, n)?]))
            .ok_or_else(|| format!("Unknown dart ({e}, {n})"))
    }

    /// `(edge, vertex)` darts around a face.
    fn darts(&self, face: usize) -> impl Iterator<Item = (&'a Id, &'a Id)> + '_ {
        self.boundaries[face].iter().map(|&(i, side)| {
            let (e, ends) = self.graph.edges.index(i);
            (e, &ends[side])
        })
    }
}

impl OrderedMap<Ends> {
    fn position(&self, e: &str) -> Option<usize> {
        self.0.get_index_of(e)
    }

    fn index(&self, i: usize) -> (&Id, &Ends) {
        self.0.get_index(i).expect("edge position")
    }

    /// Side of the edge at position `i` whose endpoint is `n`.
    fn side(&self, i: usize, n: &str) -> Option<usize> {
        self.index(i).1.iter().position(|m| m == n)
    }
}

impl Graph {
    pub(crate) fn new(nodes: Set, edges: OrderedMap<Ends>) -> Result<Self> {
        let mut rotation: BTreeMap<Id, Ids> =
            nodes.iter().map(|n| (n.clone(), Ids::new())).collect();
        for (e, [a, b]) in edges.iter() {
            require(
                nodes.contains(a) && nodes.contains(b) && a != b,
                "Invalid graph endpoints",
            )?;
            rotation.entry(a.clone()).or_default().push(e.clone());
            rotation.entry(b.clone()).or_default().push(e.clone());
        }
        Ok(Self {
            nodes,
            edges,
            rotation,
            ..Self::default()
        })
    }

    pub(crate) fn rotation_at(&self, n: &str) -> Result<&Ids> {
        self.rotation
            .get(n)
            .ok_or_else(|| format!("Unknown rotation vertex {n}"))
    }

    pub(crate) fn other(&self, e: &str, n: &str) -> Result<&Id> {
        let [a, b] = self.edges.at(e)?;
        if a == n {
            return Ok(b);
        }
        require(b == n, "Vertex not incident to edge")?;
        Ok(a)
    }

    /// The edges leaving the connected region `members`, in the cyclic order
    /// of its embedded boundary, starting with `first`.
    pub(crate) fn boundary(&self, members: &Set, first: &str) -> Result<Ids> {
        let internal: BTreeSet<&str> = self
            .edges
            .iter()
            .filter(|(_, p)| members.contains(&p[0]) && members.contains(&p[1]))
            .map(|(e, _)| e.as_str())
            .collect();
        let mut edge = first.to_owned();
        let p = self.edges.at(&edge)?;
        let mut point = if members.contains(&p[0]) {
            p[0].clone()
        } else {
            p[1].clone()
        };
        let mut visited: HashSet<Dart> = HashSet::new();
        let mut boundary = Ids::new();
        while visited.insert((edge.clone(), point.clone())) {
            let r = self.rotation_at(&point)?;
            edge = r[(index_of(r, &edge)? + 1) % r.len()].clone();
            if internal.contains(edge.as_str()) {
                point = self.other(&edge, &point)?.clone();
            } else {
                boundary.push(edge.clone());
            }
        }
        Ok(boundary)
    }

    /// Contract the embedded region `members` into its vertex `into`, whose
    /// rotation becomes the region's boundary order.
    pub(crate) fn contract(&mut self, members: &Set, into: &str) -> Result<()> {
        let leaving = |p: &Ends| members.contains(&p[0]) != members.contains(&p[1]);
        let first = self
            .edges
            .iter()
            .find(|(_, p)| leaving(p))
            .map(|(e, _)| e.clone())
            .ok_or("Contracted region has no boundary")?;
        let boundary = self.boundary(members, &first)?;
        let internal: Ids = self
            .edges
            .iter()
            .filter(|(_, p)| members.contains(&p[0]) && members.contains(&p[1]))
            .map(|(e, _)| e.clone())
            .collect();
        for e in &internal {
            self.edges.erase(e);
            self.protected_edges.remove(e);
        }
        for e in &boundary {
            for end in self.edges.at_mut(e)?.iter_mut() {
                if members.contains(end) {
                    *end = into.to_owned();
                }
            }
        }
        for n in members.iter().filter(|n| *n != into) {
            self.nodes.remove(n);
            self.rotation.remove(n);
        }
        self.rotation.insert(into.to_owned(), boundary);
        self.wheel_hubs.remove(into);
        self.hubs.remove(into);
        Ok(())
    }

    fn subgraph(&self, chosen: &Set, keep_nodes: bool) -> Result<Self> {
        let mut out = Self::default();
        if keep_nodes {
            out.nodes = self.nodes.clone();
        }
        for (e, p) in self.edges.iter() {
            if chosen.contains(e) {
                out.edges.insert(e.clone(), p.clone());
                out.nodes.extend(p.iter().cloned());
            }
        }
        for n in &out.nodes {
            let r = self.rotation_at(n)?;
            out.rotation.insert(
                n.clone(),
                r.iter().filter(|e| chosen.contains(*e)).cloned().collect(),
            );
        }
        out.protected_edges = intersection(&self.protected_edges, chosen);
        for (n, r) in &self.hubs {
            if out.nodes.contains(n) {
                out.hubs.insert(n.clone(), r.clone());
            }
        }
        out.wheel_hubs = intersection(&self.wheel_hubs, &out.nodes);
        Ok(out)
    }

    fn faces(&self) -> Result<Faces<'_>> {
        let mut incident: BTreeMap<&str, usize> = BTreeMap::new();
        for [a, b] in self.edges.values() {
            *incident.entry(a).or_default() += 1;
            if a != b {
                *incident.entry(b).or_default() += 1;
            }
        }
        // Edge preceding each dart's edge around the dart's vertex.
        let mut previous = vec![[usize::MAX; 2]; self.edges.len()];
        for n in &self.nodes {
            let invalid = || format!("Invalid rotation at {n}");
            let darts = self
                .rotation_at(n)?
                .iter()
                .map(|e| {
                    let i = self.edges.position(e)?;
                    Some((i, self.edges.side(i, n)?))
                })
                .collect::<Option<Vec<_>>>()
                .ok_or_else(invalid)?;
            if darts.len() != incident.get(n.as_str()).copied().unwrap_or_default() {
                return Err(invalid());
            }
            for (k, &(i, side)) in darts.iter().enumerate() {
                if previous[i][side] != usize::MAX {
                    return Err(invalid());
                }
                previous[i][side] = darts[(k + darts.len() - 1) % darts.len()].0;
            }
        }
        let mut out = Faces {
            graph: self,
            boundaries: Vec::new(),
            face_of: vec![[usize::MAX; 2]; self.edges.len()],
        };
        for start in (0..self.edges.len()).flat_map(|i| [(i, 0), (i, 1)]) {
            if out.face_of[start.0][start.1] != usize::MAX {
                continue;
            }
            let mut d = start;
            let mut boundary = Vec::new();
            while out.face_of[d.0][d.1] == usize::MAX {
                out.face_of[d.0][d.1] = out.boundaries.len();
                boundary.push(d);
                let (e, ends) = self.edges.index(d.0);
                let target = &ends[1 - d.1];
                let before = previous[d.0][1 - d.1];
                let side = (before != usize::MAX)
                    .then(|| self.edges.side(before, target))
                    .flatten()
                    .ok_or_else(|| format!("Unknown dart ({e}, {target})"))?;
                d = (before, side);
            }
            require(d == start, "Invalid face permutation")?;
            out.boundaries.push(boundary);
        }
        Ok(out)
    }

    fn valid_embedding(&self) -> Result<bool> {
        let fs = self.faces()?;
        let mut component_of: HashMap<&str, usize> = HashMap::new();
        let mut counts: Vec<[i64; 3]> = Vec::new();
        let mut adj: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
        for [a, b] in self.edges.values() {
            adj.entry(a).or_default().push(b);
            adj.entry(b).or_default().push(a);
        }
        for root in &self.nodes {
            if component_of.contains_key(root.as_str()) {
                continue;
            }
            let component = counts.len();
            let mut vertices = 0;
            let mut pending = vec![root.as_str()];
            while let Some(n) = pending.pop() {
                if let Entry::Vacant(entry) = component_of.entry(n) {
                    entry.insert(component);
                    vertices += 1;
                    pending.extend(adj.get(n).into_iter().flatten());
                }
            }
            counts.push([vertices, 0, 0]);
        }
        // Count each edge and face once, rather than rescanning the whole
        // graph for every component of a partly embedded edge prefix.
        for [a, _] in self.edges.values() {
            if let Some(&component) = component_of.get(a.as_str()) {
                counts[component][1] += 1;
            }
        }
        for face in 0..fs.len() {
            if let Some((_, n)) = fs.darts(face).next() {
                if let Some(&component) = component_of.get(n.as_str()) {
                    counts[component][2] += 1;
                }
            }
        }
        if counts
            .iter()
            .any(|&[nv, ne, nf]| ne != 0 && nv - ne + nf != 2)
        {
            return Ok(false);
        }
        for (hub, expected) in &self.hubs {
            if !cyclic_equal(self.rotation_at(hub)?, expected) {
                return Ok(false);
            }
        }
        for f in 0..fs.len() {
            let wheel = fs.darts(f).any(|(_, n)| self.wheel_hubs.contains(n));
            let protected_face = fs.darts(f).all(|(e, _)| self.protected_edges.contains(e));
            if wheel && (fs.boundaries[f].len() != 3 || !protected_face) {
                return Ok(false);
            }
        }
        Ok(true)
    }

    pub(crate) fn validate(&self) -> Result<()> {
        require(self.valid_embedding()?, "Invalid constrained embedding")
    }

    /// Request the SPQR decomposition of this block from the primitive.
    fn spqr(&self, spqr: &dyn SpqrDecomposer) -> Result<Option<Decomposition>> {
        spqr.decompose(&SpqrRequest::new(self))
    }

    /// Map SPQR skeletons onto graph identities and orient R skeletons by their
    /// oriented wheel hubs; `None` if one skeleton needs both orientations.
    fn feasible_skeletons(&self, tree: &Decomposition) -> Result<Option<OrderedMap<Self>>> {
        let reserved = self.edges.key_set();
        let mut out = OrderedMap::default();
        for c in &tree.components {
            let cid = &c.id;
            let mut names: BTreeMap<&str, Id> = BTreeMap::new();
            let mut real: BTreeMap<&str, Id> = BTreeMap::new();
            for e in &c.edges {
                let name = match &e.real_edge {
                    None => skeleton_id(cid, &e.id, &reserved),
                    Some(r) => r.clone(),
                };
                if let Some(r) = &e.real_edge {
                    real.insert(r, name.clone());
                }
                names.insert(&e.id, name);
            }
            let name = |e: &str| {
                names
                    .get(e)
                    .ok_or_else(|| format!("Unknown skeleton edge {e}"))
            };
            let mut s = Self {
                nodes: c.nodes.iter().cloned().collect(),
                ..Self::default()
            };
            for e in &c.edges {
                let p = match &e.real_edge {
                    None => {
                        let mut p = [e.source.clone(), e.target.clone()];
                        p.sort_unstable();
                        p
                    }
                    Some(r) => self.edges.at(r)?.clone(),
                };
                s.edges.insert(name(&e.id)?.clone(), p);
            }
            for (n, r) in &c.rotation_cw {
                for e in r {
                    s.rotation
                        .entry(n.clone())
                        .or_default()
                        .push(name(e)?.clone());
                }
            }
            for e in &self.protected_edges {
                if let Some(name) = real.get(e.as_str()) {
                    s.protected_edges.insert(name.clone());
                }
            }
            s.wheel_hubs = intersection(&self.wheel_hubs, &s.nodes);
            let mut orientation = BTreeSet::new();
            for (hub, expected) in &self.hubs {
                if !s.nodes.contains(hub) {
                    continue;
                }
                let mapped = expected
                    .iter()
                    .map(|e| {
                        real.get(e.as_str())
                            .cloned()
                            .ok_or_else(|| format!("Wheel spoke {e} is not real"))
                    })
                    .collect::<Result<Ids>>()?;
                s.hubs.insert(hub.clone(), mapped.clone());
                require(c.kind == SpqrKind::R, "Wheel hub is outside R skeleton")?;
                let mut actual = s.rotation_at(hub)?.clone();
                if cyclic_equal(&actual, &mapped) {
                    orientation.insert(false);
                } else {
                    actual.reverse();
                    require(cyclic_equal(&actual, &mapped), "SPQR changed wheel order")?;
                    orientation.insert(true);
                }
            }
            if orientation.len() > 1 {
                return Ok(None);
            }
            if orientation.contains(&true) {
                for r in s.rotation.values_mut() {
                    r.reverse();
                }
            }
            out.insert(cid.clone(), s);
        }
        Ok(Some(out))
    }

    /// Biconnected blocks as edge sets, in Hopcroft–Tarjan discovery order.
    fn blocks(&self) -> Vec<Set> {
        let mut search = BlockSearch::default();
        for (e, [a, b]) in self.edges.iter() {
            search.adj.entry(a).or_default().push((e, b));
            search.adj.entry(b).or_default().push((e, a));
        }
        for n in &self.nodes {
            if !search.discovery.contains_key(n.as_str()) {
                search.visit(n, "");
            }
        }
        search.result
    }

    /// Replace the virtual edge `pe` of this skeleton by `child` without its
    /// twin `ce`, splicing the child's pole rotations into the virtual slot.
    fn splice(self, pe: &str, child: &Self, ce: &str) -> Result<Self> {
        let poles: Set = self.edges.at(pe)?.iter().cloned().collect();
        require(
            poles == child.edges.at(ce)?.iter().cloned().collect::<Set>(),
            "Virtual edges have different poles",
        )?;
        require(
            child
                .edges
                .keys()
                .all(|e| e == pe || e == ce || !self.edges.contains(e)),
            "SPQR edge identities overlap",
        )?;
        require(
            self.nodes.intersection(&child.nodes).eq(&poles),
            "SPQR expansions share non-pole vertices",
        )?;
        let mut out = self;
        out.nodes.extend(child.nodes.iter().cloned());
        out.edges.erase(pe);
        for (e, p) in child.edges.iter() {
            if e != ce {
                out.edges.insert(e.clone(), p.clone());
            }
        }
        for n in &child.nodes {
            let r = child.rotation_at(n)?;
            if poles.contains(n) {
                replace(
                    out.rotation.entry(n.clone()).or_default(),
                    pe,
                    &after(r, ce)?,
                )?;
            } else {
                out.rotation.insert(n.clone(), r.clone());
            }
        }
        out.protected_edges.remove(pe);
        out.protected_edges
            .extend(child.protected_edges.iter().filter(|e| *e != ce).cloned());
        out.insert_hubs(child);
        out.wheel_hubs.extend(child.wheel_hubs.iter().cloned());
        Ok(out)
    }

    /// Insert the child's hubs without overwriting existing ones.
    fn insert_hubs(&mut self, child: &Self) {
        for (hub, spokes) in &child.hubs {
            self.hubs
                .entry(hub.clone())
                .or_insert_with(|| spokes.clone());
        }
    }

    /// Face path from a face at `start` across `crossings` to a face at `end`.
    fn witness(&self, start: &str, end: &str, crossings: &[Id]) -> Result<Witness> {
        let dual = Dual::new(self)?;
        let targets = dual.at(end)?;
        for source in dual.at(start)? {
            let mut face = source;
            let mut darts = Vec::new();
            let mut valid = true;
            for e in crossings {
                let i = self
                    .edges
                    .position(e)
                    .ok_or_else(|| format!("Unknown identity: {e}"))?;
                let p = self.edges.index(i).1;
                let [a, b] = dual.faces.face_of[i];
                if face == a {
                    darts.push((e.clone(), p[0].clone()));
                    face = b;
                } else if face == b {
                    darts.push((e.clone(), p[1].clone()));
                    face = a;
                } else {
                    valid = false;
                    break;
                }
            }
            if valid && targets.contains(&face) {
                let corner = |face: usize, vertex: &str| {
                    dual.faces
                        .darts(face)
                        .find(|(_, n)| *n == vertex)
                        .map(|(e, _)| e.clone())
                        .unwrap_or_default()
                };
                return Ok(Witness {
                    darts,
                    start_corner: corner(source, start),
                    end_corner: corner(face, end),
                });
            }
        }
        Err("Embedding does not realize insertion path".into())
    }

    /// First wheel-free corner at every vertex. Gluing blocks through these
    /// corners merges two wheel-free faces, so all other corner choices survive.
    fn free_corners(&self) -> Result<BTreeMap<Id, Id>> {
        let fs = self.faces()?;
        let wheel: Vec<_> = (0..fs.len())
            .map(|face| fs.darts(face).any(|(_, v)| self.wheel_hubs.contains(v)))
            .collect();
        let mut corners: BTreeMap<Id, Id> = BTreeMap::new();
        for (n, rotation) in &self.rotation {
            for e in rotation {
                if !wheel[fs.of(e, n)?] {
                    corners.insert(n.clone(), e.clone());
                    break;
                }
            }
        }
        Ok(corners)
    }

    /// Join a block at cut vertex `n`, inserting its rotation (starting after
    /// corner `cc`) after corner `pc` of this graph.
    fn join_vertex(self, child: &Self, n: &str, pc: &str, cc: &str) -> Result<Self> {
        require(
            self.nodes.intersection(&child.nodes).eq([n]),
            "Blocks must meet at one cut vertex",
        )?;
        let mut out = self;
        out.nodes.extend(child.nodes.iter().cloned());
        out.edges.update(&child.edges);
        for (v, r) in &child.rotation {
            if v != n {
                out.rotation.insert(v.clone(), r.clone());
            }
        }
        let mut order = after(child.rotation_at(n)?, cc)?;
        order.push(cc.to_owned());
        let r = out.rotation.entry(n.to_owned()).or_default();
        let i = index_of(r, pc)? + 1;
        r.splice(i..i, order);
        out.protected_edges
            .extend(child.protected_edges.iter().cloned());
        out.insert_hubs(child);
        out.wheel_hubs.extend(child.wheel_hubs.iter().cloned());
        Ok(out)
    }

    fn unite(self, b: &Self) -> Result<Self> {
        require(
            self.nodes.is_disjoint(&b.nodes),
            "Disjoint union shares vertices",
        )?;
        let mut out = self;
        out.nodes.extend(b.nodes.iter().cloned());
        out.edges.update(&b.edges);
        for (n, r) in &b.rotation {
            out.rotation.entry(n.clone()).or_insert_with(|| r.clone());
        }
        out.insert_hubs(b);
        out.protected_edges
            .extend(b.protected_edges.iter().cloned());
        out.wheel_hubs.extend(b.wheel_hubs.iter().cloned());
        Ok(out)
    }

    /// Attach embedded blocks at cut vertices through faces outside wheels.
    fn glue(&self, pending: Vec<Self>) -> Result<Self> {
        let mut pending = pending
            .into_iter()
            .map(|block| Ok((block.free_corners()?, block)))
            .collect::<Result<Vec<_>>>()?;
        let mut out = Self::default();
        let mut corners: BTreeMap<Id, Id> = BTreeMap::new();
        while !pending.is_empty() {
            let mut attached = false;
            for i in 0..pending.len() {
                let (child_corners, child) = &pending[i];
                let shared = intersection(&out.nodes, &child.nodes);
                let Some(n) = shared.first() else {
                    continue;
                };
                require(shared.len() == 1, "Blocks share multiple cut vertices")?;
                out = out.join_vertex(
                    child,
                    n,
                    corners.get(n).ok_or("Cannot attach block inside wheel")?,
                    child_corners
                        .get(n)
                        .ok_or("Cannot attach block inside wheel")?,
                )?;
                let (child_corners, _) = pending.remove(i);
                for (n, corner) in child_corners {
                    corners.entry(n).or_insert(corner);
                }
                attached = true;
                break;
            }
            if !attached {
                let (child_corners, child) = pending.remove(0);
                out = out.unite(&child)?;
                corners.extend(child_corners);
            }
        }
        out.nodes.extend(self.nodes.iter().cloned());
        for n in &self.nodes {
            out.rotation.entry(n.clone()).or_default();
        }
        Ok(out)
    }

    pub(crate) fn try_embed(&self, spqr: &dyn SpqrDecomposer) -> Result<Option<Self>> {
        let mut embedded = Vec::new();
        for es in self.blocks() {
            let mut block = self.subgraph(&es, false)?;
            if block.edges.len() > 2 {
                let Some(s) = block.spqr(spqr)? else {
                    return Ok(None);
                };
                let Some(skeletons) = block.feasible_skeletons(&s)? else {
                    return Ok(None);
                };
                let tree = Tree::new(&s, skeletons);
                block = tree.expand(&tree.root, &BTreeSet::new(), None)?;
            }
            embedded.push(block);
        }
        let out = self.glue(embedded)?;
        Ok(out.valid_embedding()?.then_some(out))
    }

    fn embed(&self, spqr: &dyn SpqrDecomposer) -> Result<Self> {
        self.try_embed(spqr)?
            .ok_or_else(|| "Expected an EC-planar graph".to_owned())
    }

    /// Optimal insertion path through one block: dynamic programming over the
    /// SPQR allocation path, choosing P permutations and R reflections.
    fn block_insertion(
        &self,
        start: &str,
        end: &str,
        spqr: &dyn SpqrDecomposer,
    ) -> Result<Insertion> {
        if self.edges.len() <= 2 {
            return Ok(Insertion::new(
                self.clone(),
                Ids::new(),
                start,
                end,
                Vec::new(),
            ));
        }
        let raw = self.spqr(spqr)?.ok_or("Insertion requires planar graph")?;
        let skeletons = self
            .feasible_skeletons(&raw)?
            .ok_or("Insertion constraints infeasible")?;
        let tree = Tree::new(&raw, skeletons);
        let path = tree.allocation_path(start, end)?;
        let boundaries: BTreeSet<(Id, Id)> =
            path.windows(2).map(|w| unordered(&w[0], &w[1])).collect();
        const INF: usize = usize::MAX / 4;
        let mut costs = [0, 0];
        let mut choices: Vec<[Choice; 2]> = Vec::new();
        let mut history = Vec::new();
        for (i, cid) in path.iter().enumerate() {
            let previous = i.checked_sub(1).map(|j| path[j].as_str());
            let following = path.get(i + 1).map(String::as_str);
            let kind = *tree.kinds.at(cid)?;
            let mut local: [Choice; 2] = Default::default();
            if kind == SpqrKind::P && path.len() > 1 {
                let (Some(previous), Some(following)) = (previous, following) else {
                    return Err("Shortest path ends at P node".into());
                };
                for (side, choice) in local.iter_mut().enumerate() {
                    let s = tree.p_embedding(cid, previous, following, side)?;
                    *choice = Choice {
                        entering: side,
                        crossed: Ids::new(),
                        graph: tree.expand(cid, &boundaries, Some(&s))?,
                        mirrored: false,
                    };
                }
            } else {
                let mut best = [INF, INF];
                for mirrored in [false, true] {
                    if mirrored && (kind != SpqrKind::R || !tree.skeletons.at(cid)?.hubs.is_empty())
                    {
                        continue;
                    }
                    let mut skeleton = tree.skeletons.at(cid)?.clone();
                    if mirrored {
                        for r in skeleton.rotation.values_mut() {
                            r.reverse();
                        }
                    }
                    let mut part = tree.expand(cid, &boundaries, Some(&skeleton))?;
                    let virtual_edge = |other: Option<&str>| {
                        other
                            .map(|other| Ok::<_, String>(tree.adjacent(cid)?.at(other)?[0].clone()))
                            .transpose()
                    };
                    let incoming: Option<Id> = virtual_edge(previous)?;
                    let outgoing: Option<Id> = virtual_edge(following)?;
                    part.protected_edges
                        .extend(incoming.iter().chain(&outgoing).cloned());
                    let dual = Dual::new(&part)?;
                    for (entering, &entry_cost) in costs.iter().enumerate() {
                        let sources = match &incoming {
                            Some(e) => dual.side(e, 1 - entering)?,
                            None => dual.at(start)?,
                        };
                        for leaving in 0..2 {
                            let targets = match &outgoing {
                                Some(e) => dual.side(e, leaving)?,
                                None => dual.at(end)?,
                            };
                            let Some(crossed) = dual.shortest(&sources, &targets) else {
                                continue;
                            };
                            if entry_cost == INF {
                                continue;
                            }
                            let cost = entry_cost + crossed.len();
                            if cost < best[leaving] {
                                best[leaving] = cost;
                                local[leaving] = Choice {
                                    entering,
                                    crossed,
                                    graph: part.clone(),
                                    mirrored,
                                };
                            }
                        }
                    }
                }
                costs = best;
            }
            choices.push(local);
            history.push(InsertionState {
                component: cid.clone(),
                kind,
                costs: costs.map(|c| (c != INF).then_some(c)),
                mirrored: false,
                entering: 0,
                leaving: 0,
            });
        }
        if costs[0].min(costs[1]) == INF {
            return Err("No EC insertion avoids protected edges".into());
        }
        require(costs[0] == costs[1], "Terminal insertion costs differ")?;
        let mut selected: Vec<Choice> = Vec::with_capacity(path.len());
        let mut side = 0;
        for (choice, state) in choices.iter_mut().zip(&mut history).rev() {
            let choice = std::mem::take(&mut choice[side]);
            state.mirrored = choice.mirrored;
            state.entering = choice.entering;
            state.leaving = side;
            side = choice.entering;
            selected.push(choice);
        }
        selected.reverse();
        let mut selected = selected.into_iter();
        let first = selected.next().ok_or("Empty SPQR allocation path")?;
        let mut embedding = first.graph;
        let mut crossings = first.crossed;
        for (i, choice) in selected.enumerate() {
            let pair = tree.adjacent(&path[i + 1])?.at(&path[i])?;
            embedding = choice.graph.splice(&pair[0], &embedding, &pair[1])?;
            crossings.extend(choice.crossed);
        }
        embedding.protected_edges = intersection(&embedding.protected_edges, &self.protected_edges);
        embedding.witness(start, end, &crossings)?;
        Ok(Insertion::new(embedding, crossings, start, end, history))
    }

    /// Optimal insertion path through the blocks between `start` and `end`.
    fn optimal_insertion(
        &self,
        start: &str,
        end: &str,
        spqr: &dyn SpqrDecomposer,
    ) -> Result<Insertion> {
        require(
            start != end && self.nodes.contains(start) && self.nodes.contains(end),
            "Insertion needs two distinct vertices",
        )?;
        let mut subgraphs = Vec::new();
        let mut incidence: BTreeMap<Id, Vec<usize>> = BTreeMap::new();
        for es in self.blocks() {
            let b = self.subgraph(&es, false)?;
            for n in &b.nodes {
                incidence
                    .entry(n.clone())
                    .or_default()
                    .push(subgraphs.len());
            }
            subgraphs.push(b);
        }
        // Breadth-first search in the block-cut tree.
        #[derive(Clone, PartialEq, Eq, Hash)]
        enum Key {
            Block(usize),
            Vertex(Id),
        }
        let source = Key::Vertex(start.to_owned());
        let target = Key::Vertex(end.to_owned());
        let mut queue = VecDeque::from([source.clone()]);
        let mut pred: HashMap<Key, Option<Key>> = HashMap::from([(source, None)]);
        while !pred.contains_key(&target) {
            let Some(current) = queue.pop_front() else {
                break;
            };
            let adjacent: Vec<Key> = match &current {
                Key::Vertex(n) => incidence
                    .get(n)
                    .into_iter()
                    .flatten()
                    .map(|&i| Key::Block(i))
                    .collect(),
                Key::Block(i) => subgraphs[*i]
                    .nodes
                    .iter()
                    .map(|n| Key::Vertex(n.clone()))
                    .collect(),
            };
            for n in adjacent {
                if !pred.contains_key(&n) {
                    pred.insert(n.clone(), Some(current.clone()));
                    queue.push_back(n);
                }
            }
        }
        require(
            pred.contains_key(&target),
            "Insertion endpoints are disconnected",
        )?;
        let mut sequence = Vec::new();
        let mut cursor = Some(target);
        while let Some(key) = cursor {
            cursor = pred[&key].clone();
            sequence.push(key);
        }
        sequence.reverse();
        let mut embedding: Option<Self> = None;
        let mut crossings = Ids::new();
        let mut used = BTreeSet::new();
        let mut history = Vec::new();
        for step in sequence.windows(3).step_by(2) {
            let [Key::Vertex(x), Key::Block(i), Key::Vertex(y)] = step else {
                return Err("Invalid block-cut path".into());
            };
            let step = subgraphs[*i].block_insertion(x, y, spqr)?;
            let joined = match embedding {
                None => step.embedding.clone(),
                Some(current) => {
                    let left = current.witness(start, x, &crossings)?;
                    let right = step.embedding.witness(x, y, &step.crossings)?;
                    current.join_vertex(
                        &step.embedding,
                        x,
                        &left.end_corner,
                        &right.start_corner,
                    )?
                }
            };
            crossings.extend(step.crossings);
            used.insert(*i);
            history.extend(step.states);
            joined.witness(start, y, &crossings)?;
            embedding = Some(joined);
        }
        let mut embedding = embedding.ok_or("Insertion endpoints are disconnected")?;
        let mut corners = embedding.free_corners()?;
        let mut pending: BTreeSet<usize> =
            (0..subgraphs.len()).filter(|i| !used.contains(i)).collect();
        while let Some(&first) = pending.first() {
            let mut progress = false;
            for i in pending.clone() {
                let current = &subgraphs[i];
                let shared = intersection(&embedding.nodes, &current.nodes);
                let Some(n) = shared.first() else {
                    continue;
                };
                require(shared.len() == 1, "Invalid block decomposition")?;
                let child_corners = current.free_corners()?;
                embedding = embedding.join_vertex(
                    current,
                    n,
                    corners.get(n).ok_or("Cannot attach block inside wheel")?,
                    child_corners
                        .get(n)
                        .ok_or("Cannot attach block inside wheel")?,
                )?;
                for (n, corner) in child_corners {
                    corners.entry(n).or_insert(corner);
                }
                pending.remove(&i);
                progress = true;
            }
            if !progress {
                embedding = embedding.unite(&subgraphs[first])?;
                corners.extend(subgraphs[first].free_corners()?);
                pending.remove(&first);
            }
        }
        embedding.nodes.extend(self.nodes.iter().cloned());
        for n in &self.nodes {
            embedding.rotation.entry(n.clone()).or_default();
        }
        embedding.witness(start, end, &crossings)?;
        Ok(Insertion::new(embedding, crossings, start, end, history))
    }

    /// Embed a maximal EC-planar subgraph greedily, then insert each remaining
    /// edge optimally while preserving earlier crossings as mirror wheels.
    pub(crate) fn planarize(&self, spqr: &dyn SpqrDecomposer) -> Result<Planarization> {
        if let Some(initial) = self.try_embed(spqr)? {
            return Ok(Planarization {
                embedding: initial,
                chains: self
                    .edges
                    .keys()
                    .map(|e| (e.clone(), vec![e.clone()]))
                    .collect(),
                ..Planarization::default()
            });
        }
        let mut chosen = self.protected_edges.clone();
        let mut embedded = self.subgraph(&chosen, true)?.embed(spqr)?;
        let mut removed = Ids::new();
        for e in self.edges.keys() {
            if chosen.contains(e) {
                continue;
            }
            let mut candidate = chosen.clone();
            candidate.insert(e.clone());
            if let Some(next) = self.subgraph(&candidate, true)?.try_embed(spqr)? {
                chosen = candidate;
                embedded = next;
            } else {
                removed.push(e.clone());
            }
        }
        let mut out = Planarization {
            embedding: embedded,
            chains: chosen
                .iter()
                .map(|e| (e.clone(), vec![e.clone()]))
                .collect(),
            ..Planarization::default()
        };
        for e in &removed {
            let owners = out.chains.owners();
            let constraints = out
                .crossing_nodes
                .keys()
                .map(|n| {
                    let children = out
                        .embedding
                        .rotation_at(n)?
                        .iter()
                        .cloned()
                        .map(ConstraintTree::Leaf)
                        .collect();
                    Ok((
                        n.clone(),
                        ConstraintTree::Node {
                            kind: ConstraintKind::Mirror,
                            children,
                        },
                    ))
                })
                .collect::<Result<OrderedMap<_>>>()?;
            let preserved = Expansion::new(&out.embedding, constraints)?;
            let constrained = preserved.graph.embed(spqr)?;
            let ends = self.edges.at(e)?;
            let path = constrained.optimal_insertion(&ends[0], &ends[1], spqr)?;
            let step = path.planarize(e)?;
            out.embedding = preserved.collapse(&step.embedding)?;
            for (n, record) in step.crossing_nodes.iter() {
                let owner = owners
                    .get(record[0].as_str())
                    .ok_or_else(|| format!("Unknown chain part {}", record[0]))?;
                out.crossing_nodes
                    .insert(n.clone(), [(*owner).to_owned(), e.clone()]);
            }
            for (_, chain) in out.chains.iter_mut() {
                let mut next = Ids::new();
                for part in chain.iter() {
                    next.extend(step.chains.at(part)?.iter().cloned());
                }
                *chain = next;
            }
            out.chains.insert(e.clone(), step.chains.at(e)?.clone());
            let owners = out.chains.owners();
            for (n, r) in out.crossing_nodes.iter() {
                let incidence = out
                    .embedding
                    .rotation_at(n)?
                    .iter()
                    .map(|part| {
                        owners
                            .get(part.as_str())
                            .map(|owner| (*owner).to_owned())
                            .ok_or_else(|| format!("Unknown chain part {part}"))
                    })
                    .collect::<Result<Ids>>()?;
                require(
                    incidence.len() == 4
                        && incidence[0] == incidence[2]
                        && incidence[1] == incidence[3]
                        && incidence[0] != incidence[1],
                    "Crossing alternation changed",
                )?;
                require(
                    as_set(&incidence) == as_set(r),
                    "Crossing provenance changed",
                )?;
            }
            let crossed = out
                .crossing_nodes
                .iter()
                .filter(|(n, _)| step.crossing_nodes.contains(n))
                .map(|(_, r)| r[0].clone())
                .collect();
            out.insertions.push(InsertionRecord {
                edge: e.clone(),
                crossed_edges: crossed,
                crossings: path.crossings.len(),
                states: path.states,
            });
        }
        require(
            out.chains.key_set() == self.edges.key_set(),
            "Planarization lost original edges",
        )?;
        out.embedding.validate()?;
        Ok(out)
    }
}

#[derive(Default)]
struct BlockSearch<'a> {
    adj: BTreeMap<&'a str, Vec<(&'a str, &'a str)>>,
    discovery: HashMap<&'a str, usize>,
    low: HashMap<&'a str, usize>,
    stack: Vec<&'a str>,
    result: Vec<Set>,
}

impl<'a> BlockSearch<'a> {
    fn visit(&mut self, u: &'a str, parent: &str) {
        let index = self.discovery.len();
        self.discovery.insert(u, index);
        self.low.insert(u, index);
        let neighbors = self.adj.get(u).cloned().unwrap_or_default();
        for (e, v) in neighbors {
            if e == parent {
                continue;
            }
            if let Some(&dv) = self.discovery.get(v) {
                if dv < self.discovery[u] {
                    self.stack.push(e);
                    self.low.insert(u, self.low[u].min(dv));
                }
                continue;
            }
            self.stack.push(e);
            self.visit(v, e);
            self.low.insert(u, self.low[u].min(self.low[v]));
            if self.low[v] >= self.discovery[u] {
                let mut component = Set::new();
                while let Some(current) = self.stack.pop() {
                    component.insert(current.to_owned());
                    if current == e {
                        break;
                    }
                }
                self.result.push(component);
            }
        }
    }
}

/// Graph identity of a virtual skeleton edge, avoiding real edge identities.
fn skeleton_id(component: &str, edge: &str, reserved: &Set) -> Id {
    let mut id = format!("@spqr:{component}:{edge}");
    while reserved.contains(&id) {
        id.insert(0, '@');
    }
    id
}

/// SPQR tree with skeletons expressed in graph identities.
struct Tree {
    kinds: OrderedMap<SpqrKind>,
    skeletons: OrderedMap<Graph>,
    /// Virtual edge pairs `[own, twin]` towards each neighboring tree node.
    adjacency: BTreeMap<Id, OrderedMap<Ends>>,
    root: Id,
}

impl Tree {
    fn new(s: &Decomposition, skeletons: OrderedMap<Graph>) -> Self {
        let mut kinds = OrderedMap::default();
        let mut reserved = Set::new();
        for c in &s.components {
            kinds.insert(c.id.clone(), c.kind);
            reserved.extend(c.edges.iter().filter_map(|e| e.real_edge.clone()));
        }
        let mut adjacency: BTreeMap<Id, OrderedMap<Ends>> = BTreeMap::new();
        for c in &s.components {
            for e in &c.edges {
                if let Some(t) = &e.twin {
                    adjacency.entry(c.id.clone()).or_default().insert(
                        t.component.clone(),
                        [
                            skeleton_id(&c.id, &e.id, &reserved),
                            skeleton_id(&t.component, &t.edge, &reserved),
                        ],
                    );
                }
            }
        }
        Self {
            kinds,
            skeletons,
            adjacency,
            root: s.root.clone(),
        }
    }

    fn adjacent(&self, n: &str) -> Result<&OrderedMap<Ends>> {
        self.adjacency
            .get(n)
            .ok_or_else(|| format!("SPQR node {n} has no neighbors"))
    }

    /// Expand the tree from `start`, optionally replacing its skeleton and not
    /// crossing the `blocked` tree edges.
    fn expand(
        &self,
        start: &str,
        blocked: &BTreeSet<(Id, Id)>,
        replacement: Option<&Graph>,
    ) -> Result<Graph> {
        self.expand_from(start, "", blocked, replacement)
    }

    fn expand_from(
        &self,
        n: &str,
        parent: &str,
        blocked: &BTreeSet<(Id, Id)>,
        replacement: Option<&Graph>,
    ) -> Result<Graph> {
        // Only the start node has no parent: tree recursion never revisits it.
        let mut out = match replacement {
            Some(g) if parent.is_empty() => g.clone(),
            _ => self.skeletons.at(n)?.clone(),
        };
        if let Some(adjacent) = self.adjacency.get(n) {
            for (next, p) in adjacent.iter() {
                if next != parent && !blocked.contains(&unordered(n, next)) {
                    let child = self.expand_from(next, n, blocked, None)?;
                    out = out.splice(&p[0], &child, &p[1])?;
                }
            }
        }
        Ok(out)
    }

    /// Shortest tree path from a skeleton containing `start` to one
    /// containing `end`.
    fn allocation_path(&self, start: &str, end: &str) -> Result<Ids> {
        let mut starts = Vec::new();
        let mut ends = BTreeSet::new();
        for (k, g) in self.skeletons.iter() {
            if g.nodes.contains(start) {
                starts.push(k.as_str());
            }
            if g.nodes.contains(end) {
                ends.insert(k.as_str());
            }
        }
        starts.sort_unstable();
        let mut pred: HashMap<&str, Option<&str>> = HashMap::new();
        let mut queue = VecDeque::new();
        for &k in &starts {
            pred.insert(k, None);
            queue.push_back(k);
        }
        while let Some(n) = queue.pop_front() {
            if ends.contains(n) {
                let mut path = Ids::new();
                let mut cursor = Some(n);
                while let Some(current) = cursor {
                    path.push(current.to_owned());
                    cursor = pred[current];
                }
                path.reverse();
                return Ok(path);
            }
            if let Some(adjacent) = self.adjacency.get(n) {
                for next in adjacent.keys() {
                    if !pred.contains_key(next.as_str()) {
                        pred.insert(next, Some(n));
                        queue.push_back(next);
                    }
                }
            }
        }
        Err("Missing SPQR endpoint".into())
    }

    /// P skeleton with the outgoing virtual edge first and the incoming one on
    /// the requested side of the remaining parallel edges.
    fn p_embedding(
        &self,
        cid: &str,
        previous: &str,
        following: &str,
        side: usize,
    ) -> Result<Graph> {
        let mut g = self.skeletons.at(cid)?.clone();
        let incoming = self.adjacent(cid)?.at(previous)?[0].clone();
        let outgoing = self.adjacent(cid)?.at(following)?[0].clone();
        let poles = g.edges.at(&outgoing)?.clone();
        let middle: Ids = g
            .rotation_at(&poles[0])?
            .iter()
            .filter(|e| **e != incoming && **e != outgoing)
            .cloned()
            .collect();
        let mut order = vec![outgoing];
        if side == 0 {
            order.push(incoming.clone());
        }
        order.extend(middle);
        if side == 1 {
            order.push(incoming);
        }
        g.rotation.insert(poles[0].clone(), order.clone());
        order.reverse();
        g.rotation.insert(poles[1].clone(), order);
        Ok(g)
    }
}

/// Faces of an embedding with arcs across its unprotected edges.
struct Dual<'a> {
    faces: Faces<'a>,
    arcs: Vec<Vec<(usize, &'a str)>>,
}

impl<'a> Dual<'a> {
    fn new(g: &'a Graph) -> Result<Self> {
        let faces = g.faces()?;
        let mut arcs = vec![Vec::new(); faces.len()];
        for (e, &[a, b]) in g.edges.keys().zip(&faces.face_of) {
            if !g.protected_edges.contains(e) && a != b {
                arcs[a].push((b, e.as_str()));
                arcs[b].push((a, e.as_str()));
            }
        }
        Ok(Self { faces, arcs })
    }

    fn at(&self, n: &str) -> Result<BTreeSet<usize>> {
        self.faces
            .graph
            .rotation_at(n)?
            .iter()
            .map(|e| self.faces.of(e, n))
            .collect()
    }

    fn side(&self, e: &str, s: usize) -> Result<BTreeSet<usize>> {
        Ok(BTreeSet::from([self
            .faces
            .of(e, &self.faces.graph.edges.at(e)?[s])?]))
    }

    /// Edges crossed by a breadth-first shortest face path between the sets.
    fn shortest(&self, sources: &BTreeSet<usize>, targets: &BTreeSet<usize>) -> Option<Ids> {
        let mut pred: HashMap<usize, Option<(usize, &str)>> = HashMap::new();
        let mut queue = VecDeque::new();
        for &n in sources {
            pred.insert(n, None);
            queue.push_back(n);
        }
        while let Some(mut n) = queue.pop_front() {
            if targets.contains(&n) {
                let mut path = Ids::new();
                while let Some((from, e)) = pred[&n] {
                    path.push(e.to_owned());
                    n = from;
                }
                path.reverse();
                return Some(path);
            }
            for &(next, e) in &self.arcs[n] {
                if let Entry::Vacant(slot) = pred.entry(next) {
                    slot.insert(Some((n, e)));
                    queue.push_back(next);
                }
            }
        }
        None
    }
}

/// Face path realizing an insertion: the crossed darts and the edges whose
/// corners at the endpoints open into the first and last faces.
struct Witness {
    darts: Vec<Dart>,
    start_corner: Id,
    end_corner: Id,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "lowercase")]
pub(crate) enum ConstraintKind {
    Group,
    Mirror,
    Oriented,
}

/// Cyclic embedding constraint at a vertex: leaves are incident edge identities.
#[derive(Clone, Debug, PartialEq, Serialize)]
#[serde(untagged)]
pub(crate) enum ConstraintTree {
    Leaf(Id),
    Node {
        kind: ConstraintKind,
        children: Vec<ConstraintTree>,
    },
}

impl ConstraintTree {
    /// A group of the children, or the only child itself.
    pub(crate) fn group(mut children: Vec<Self>) -> Self {
        if children.len() == 1 {
            children.remove(0)
        } else {
            Self::Node {
                kind: ConstraintKind::Group,
                children,
            }
        }
    }

    /// Children in a fixed cyclic orientation, or the only child itself.
    pub(crate) fn oriented(mut children: Vec<Self>) -> Self {
        if children.len() == 1 {
            children.remove(0)
        } else {
            Self::Node {
                kind: ConstraintKind::Oriented,
                children,
            }
        }
    }

    fn leaves(&self) -> Result<Ids> {
        match self {
            Self::Leaf(id) => Ok(vec![id.clone()]),
            Self::Node { children, .. } => {
                require(children.len() >= 2, "Constraint node needs two children")?;
                let mut out = Ids::new();
                for child in children {
                    out.extend(child.leaves()?);
                }
                Ok(out)
            }
        }
    }

    /// Whether some rotation of the cyclic order is a linear order allowed by
    /// the tree: children consecutive, oriented nodes forward, mirror nodes in
    /// either direction and groups permuted freely.
    fn admissible(&self, rotation: &[Id]) -> Result<bool> {
        let expected = self.leaves()?;
        if rotation.len() != expected.len() || as_set(rotation) != as_set(&expected) {
            return Ok(false);
        }
        for i in 0..rotation.len() {
            if self.linear(&[&rotation[i..], &rotation[..i]].concat())? {
                return Ok(true);
            }
        }
        Ok(false)
    }

    fn linear(&self, order: &[Id]) -> Result<bool> {
        let (kind, children) = match self {
            Self::Leaf(id) => return Ok(order.len() == 1 && order[0] == *id),
            Self::Node { kind, children } => (*kind, children),
        };
        let mut owner: BTreeMap<Id, usize> = BTreeMap::new();
        for (i, child) in children.iter().enumerate() {
            for e in child.leaves()? {
                owner.insert(e, i);
            }
        }
        let mut runs: Vec<(usize, Ids)> = Vec::new();
        for e in order {
            let i = *owner
                .get(e)
                .ok_or_else(|| format!("Unknown constraint leaf {e}"))?;
            match runs.last_mut() {
                Some((last, run)) if *last == i => run.push(e.clone()),
                _ => runs.push((i, vec![e.clone()])),
            }
        }
        let unique: BTreeSet<usize> = runs.iter().map(|r| r.0).collect();
        if runs.len() != children.len() || unique.len() != children.len() {
            return Ok(false);
        }
        let forward = runs.iter().enumerate().all(|(i, r)| r.0 == i);
        let reverse = runs
            .iter()
            .enumerate()
            .all(|(i, r)| r.0 == runs.len() - 1 - i);
        if kind == ConstraintKind::Oriented && !forward {
            return Ok(false);
        }
        if kind == ConstraintKind::Mirror && !forward && !reverse {
            return Ok(false);
        }
        for (i, r) in &runs {
            if !children[*i].linear(r)? {
                return Ok(false);
            }
        }
        Ok(true)
    }
}

/// Constrained vertices replaced by their expanded constraint trees: groups
/// become single vertices, mirror/oriented nodes wheels of protected edges.
pub(crate) struct Expansion {
    original: Graph,
    pub(crate) graph: Graph,
    constraints: OrderedMap<ConstraintTree>,
    pub(crate) members: BTreeMap<Id, Set>,
    pub(crate) leaf_ports: BTreeMap<Id, BTreeMap<Id, Id>>,
    /// Protected edges joining a constraint node to its parent node.
    pub(crate) tree_edges: Set,
    reserved: Set,
    counter: usize,
}

impl Expansion {
    pub(crate) fn new(g: &Graph, constraints: OrderedMap<ConstraintTree>) -> Result<Self> {
        let mut out = Self {
            original: g.clone(),
            graph: Graph::default(),
            constraints: OrderedMap::default(),
            members: BTreeMap::new(),
            leaf_ports: BTreeMap::new(),
            tree_edges: Set::new(),
            reserved: g.nodes.iter().chain(g.edges.keys()).cloned().collect(),
            counter: 0,
        };
        for n in &g.nodes {
            out.members.insert(n.clone(), Set::from([n.clone()]));
            out.leaf_ports.insert(n.clone(), BTreeMap::new());
            if !constraints.contains(n) {
                out.graph.nodes.insert(n.clone());
                out.graph.rotation.insert(n.clone(), Ids::new());
            }
        }
        out.graph.hubs = g.hubs.clone();
        out.graph.wheel_hubs = g.wheel_hubs.clone();
        for (n, t) in constraints.iter() {
            require(
                g.nodes.contains(n) && !g.wheel_hubs.contains(n),
                "Invalid constrained vertex",
            )?;
            let incidence: Set = g
                .edges
                .iter()
                .filter(|(_, p)| p[0] == *n || p[1] == *n)
                .map(|(e, _)| e.clone())
                .collect();
            let ls = t.leaves()?;
            require(
                ls.len() == as_set(&ls).len() && as_set(&ls) == incidence,
                "Constraint must contain all incidences once",
            )?;
            out.members.entry(n.clone()).or_default().clear();
            out.expand(n, t, None);
        }
        for (e, p) in g.edges.iter() {
            let port = |n: &Id| {
                out.leaf_ports
                    .get(n)
                    .and_then(|ports| ports.get(e))
                    .unwrap_or(n)
                    .clone()
            };
            let (a, b) = (port(&p[0]), port(&p[1]));
            out.edge(e.clone(), &a, &b, false);
        }
        out.graph
            .protected_edges
            .extend(g.protected_edges.iter().cloned());
        out.constraints = constraints;
        Ok(out)
    }

    fn fresh(&mut self, prefix: &str) -> Id {
        loop {
            self.counter += 1;
            let id = format!("__ec_{prefix}_{}", self.counter);
            if self.reserved.insert(id.clone()) {
                return id;
            }
        }
    }

    fn node(&mut self, owner: &str, prefix: &str) -> Id {
        let n = self.fresh(prefix);
        self.graph.nodes.insert(n.clone());
        self.graph.rotation.insert(n.clone(), Ids::new());
        self.members
            .entry(owner.to_owned())
            .or_default()
            .insert(n.clone());
        n
    }

    fn edge(&mut self, e: Id, a: &str, b: &str, protect: bool) -> Id {
        self.graph
            .edges
            .insert(e.clone(), [a.to_owned(), b.to_owned()]);
        self.graph
            .rotation
            .entry(a.to_owned())
            .or_default()
            .push(e.clone());
        self.graph
            .rotation
            .entry(b.to_owned())
            .or_default()
            .push(e.clone());
        if protect {
            self.graph.protected_edges.insert(e.clone());
        }
        e
    }

    fn set_port(&mut self, owner: &str, leaf: &str, port: Id) {
        self.leaf_ports
            .entry(owner.to_owned())
            .or_default()
            .insert(leaf.to_owned(), port);
    }

    fn expand(&mut self, owner: &str, tree: &ConstraintTree, parent: Option<&str>) {
        let (kind, children) = match tree {
            ConstraintTree::Leaf(leaf) => {
                let port = self.node(owner, "leaf");
                self.set_port(owner, leaf, port);
                return;
            }
            ConstraintTree::Node { kind, children } => (*kind, children),
        };
        let degree = children.len() + usize::from(parent.is_some());
        let ports: Ids = if kind == ConstraintKind::Group {
            vec![self.node(owner, "group"); degree]
        } else {
            let hub = self.node(owner, "hub");
            self.graph.wheel_hubs.insert(hub.clone());
            let rim: Ids = (0..2 * degree).map(|_| self.node(owner, "rim")).collect();
            let mut spokes = Ids::new();
            for (i, r) in rim.iter().enumerate() {
                let spoke = self.fresh("spoke");
                spokes.push(self.edge(spoke, &hub, r, true));
                let rim_edge = self.fresh("rimedge");
                self.edge(rim_edge, r, &rim[(i + 1) % rim.len()], true);
            }
            if kind == ConstraintKind::Oriented {
                self.graph.hubs.insert(hub, spokes);
            }
            rim.into_iter().step_by(2).collect()
        };
        if let (Some(parent), Some(port)) = (parent, ports.last()) {
            let tree_edge = self.fresh("treeedge");
            self.edge(tree_edge.clone(), parent, port, true);
            self.tree_edges.insert(tree_edge);
        }
        for (child, port) in children.iter().zip(&ports) {
            match child {
                ConstraintTree::Leaf(leaf) => self.set_port(owner, leaf, port.clone()),
                ConstraintTree::Node { .. } => self.expand(owner, child, Some(port)),
            }
        }
    }

    /// Contract expanded constraint trees back to their vertices, reading each
    /// vertex rotation from the boundary of its embedded expansion.
    pub(crate) fn collapse(&self, embedded: &Graph) -> Result<Graph> {
        let mut owner: BTreeMap<&str, &str> = BTreeMap::new();
        for (n, ns) in &self.members {
            for m in ns {
                owner.insert(m, n);
            }
        }
        let of = |n: &str| owner.get(n).copied().unwrap_or(n).to_owned();
        let ns: Set = embedded.nodes.iter().map(|n| of(n)).collect();
        let mut edges = OrderedMap::default();
        for (e, p) in embedded.edges.iter() {
            let (a, b) = (of(&p[0]), of(&p[1]));
            if a != b {
                edges.insert(e.clone(), [a, b]);
            }
        }
        let mut out = Graph::new(ns, edges.clone())?;
        out.protected_edges = intersection(&embedded.protected_edges, &edges.key_set());
        out.hubs = self.original.hubs.clone();
        out.wheel_hubs = self.original.wheel_hubs.clone();
        for n in &embedded.nodes {
            if !owner.contains_key(n.as_str()) {
                out.rotation
                    .insert(n.clone(), embedded.rotation_at(n)?.clone());
            }
        }
        for n in &self.original.nodes {
            let ms = self
                .members
                .get(n)
                .ok_or_else(|| format!("Unknown expansion owner {n}"))?;
            let incident: Ids = edges
                .iter()
                .filter(|(_, p)| p[0] == *n || p[1] == *n)
                .map(|(e, _)| e.clone())
                .collect();
            let Some(first) = incident.first() else {
                out.rotation.insert(n.clone(), Ids::new());
                continue;
            };
            let boundary = embedded.boundary(ms, first)?;
            require(
                boundary.len() == incident.len() && as_set(&boundary) == as_set(&incident),
                "Expansion boundary lost incidences",
            )?;
            out.rotation.insert(n.clone(), boundary);
        }
        out.validate()?;
        for (n, t) in self.constraints.iter() {
            let r = out.rotation_at(n)?;
            if as_set(r) == as_set(&t.leaves()?) {
                require(t.admissible(r)?, "Collapsed constraint violated")?;
            }
        }
        Ok(out)
    }
}

#[derive(Clone, Default)]
struct Choice {
    entering: usize,
    crossed: Ids,
    graph: Graph,
    mirrored: bool,
}

/// Dynamic-programming state of one SPQR node on the allocation path.
#[derive(Clone, Debug, Serialize)]
struct InsertionState {
    component: Id,
    #[serde(rename = "type")]
    kind: SpqrKind,
    costs: [Option<usize>; 2],
    mirrored: bool,
    entering: usize,
    leaving: usize,
}

struct Insertion {
    embedding: Graph,
    crossings: Ids,
    start: Id,
    end: Id,
    states: Vec<InsertionState>,
}

#[derive(Clone, Debug, Serialize)]
pub(crate) struct InsertionRecord {
    edge: Id,
    crossed_edges: Ids,
    crossings: usize,
    states: Vec<InsertionState>,
}

/// Planar drawing graph with each original edge's chain of drawing edges and
/// the two original edges meeting at every artificial crossing vertex.
#[derive(Default)]
pub(crate) struct Planarization {
    pub(crate) embedding: Graph,
    pub(crate) chains: OrderedMap<Ids>,
    pub(crate) crossing_nodes: OrderedMap<Ends>,
    pub(crate) insertions: Vec<InsertionRecord>,
}

impl OrderedMap<Ids> {
    /// Original edge owning each drawing edge of these chains.
    fn owners(&self) -> BTreeMap<&str, &str> {
        let mut owners = BTreeMap::new();
        for (owner, chain) in self.iter() {
            for part in chain {
                owners.insert(part.as_str(), owner.as_str());
            }
        }
        owners
    }
}

fn fresh_id(base: &str, existing: &Set) -> Id {
    let mut value = base.to_owned();
    let mut serial = 0;
    while existing.contains(&value) {
        serial += 1;
        value = format!("{base}:{serial}");
    }
    value
}

impl Insertion {
    fn new(
        embedding: Graph,
        crossings: Ids,
        start: &str,
        end: &str,
        states: Vec<InsertionState>,
    ) -> Self {
        Self {
            embedding,
            crossings,
            start: start.to_owned(),
            end: end.to_owned(),
            states,
        }
    }

    /// Insert edge `id` along the path, subdividing every crossed edge at a new
    /// degree-four crossing vertex.
    fn planarize(&self, id: &str) -> Result<Planarization> {
        let g = &self.embedding;
        require(!g.edges.contains(id), "Inserted edge exists")?;
        require(
            self.crossings.len() == as_set(&self.crossings).len(),
            "Path crosses edge twice",
        )?;
        let proof = g.witness(&self.start, &self.end, &self.crossings)?;
        let mut out = Planarization {
            embedding: g.clone(),
            chains: g
                .edges
                .keys()
                .map(|e| (e.clone(), vec![e.clone()]))
                .collect(),
            ..Planarization::default()
        };
        let mut nodes = vec![self.start.clone()];
        let mut split: Vec<Ends> = Vec::new();
        for (i, (e, origin)) in proof.darts.iter().enumerate() {
            require(
                !g.protected_edges.contains(e),
                "Insertion crosses protected edge",
            )?;
            let n = fresh_id(&format!("@cross:{id}:{i}"), &out.embedding.nodes);
            out.embedding.nodes.insert(n.clone());
            nodes.push(n.clone());
            let p = out.embedding.edges.at(e)?.clone();
            out.embedding.edges.erase(e);
            let a = fresh_id(&format!("{e}:part:0"), &out.embedding.edges.key_set());
            out.embedding
                .edges
                .insert(a.clone(), [p[0].clone(), n.clone()]);
            let b = fresh_id(&format!("{e}:part:1"), &out.embedding.edges.key_set());
            out.embedding
                .edges
                .insert(b.clone(), [n.clone(), p[1].clone()]);
            replace(
                out.embedding.rotation.entry(p[0].clone()).or_default(),
                e,
                std::slice::from_ref(&a),
            )?;
            replace(
                out.embedding.rotation.entry(p[1].clone()).or_default(),
                e,
                std::slice::from_ref(&b),
            )?;
            out.chains.insert(e.clone(), vec![a.clone(), b.clone()]);
            split.push(if *origin == p[0] { [a, b] } else { [b, a] });
            out.crossing_nodes.insert(n, [e.clone(), id.to_owned()]);
        }
        nodes.push(self.end.clone());
        let mut inserted = Ids::new();
        for (i, pair) in nodes.windows(2).enumerate() {
            let e = if self.crossings.is_empty() {
                id.to_owned()
            } else {
                fresh_id(&format!("{id}:part:{i}"), &out.embedding.edges.key_set())
            };
            out.embedding
                .edges
                .insert(e.clone(), [pair[0].clone(), pair[1].clone()]);
            inserted.push(e);
        }
        for (i, [a, b]) in split.iter().enumerate() {
            out.embedding.rotation.insert(
                nodes[i + 1].clone(),
                vec![
                    a.clone(),
                    inserted[i + 1].clone(),
                    b.clone(),
                    inserted[i].clone(),
                ],
            );
        }
        let corners = [
            (&self.start, &proof.start_corner, &inserted[0]),
            (&self.end, &proof.end_corner, &inserted[inserted.len() - 1]),
        ];
        for (v, c, e) in corners {
            let mut c = c.clone();
            if !out.embedding.edges.contains(&c) {
                for part in out.chains.at(&c)? {
                    let p = out.embedding.edges.at(part)?;
                    if p[0] == *v || p[1] == *v {
                        c = part.clone();
                        break;
                    }
                }
            }
            let r = out.embedding.rotation.entry(v.clone()).or_default();
            let i = index_of(r, &c)? + 1;
            r.insert(i, e.clone());
        }
        out.chains.insert(id.to_owned(), inserted);
        Ok(out)
    }
}

#[cfg(test)]
mod corner_tests {
    use super::*;

    fn wheel(prefix: &str, a: &str) -> Graph {
        let (h, b, c) = (
            format!("{prefix}:h"),
            format!("{prefix}:b"),
            format!("{prefix}:c"),
        );
        let edges = [
            ("ha", h.as_str(), a),
            ("hb", h.as_str(), b.as_str()),
            ("hc", h.as_str(), c.as_str()),
            ("ab", a, b.as_str()),
            ("bc", b.as_str(), c.as_str()),
            ("ca", c.as_str(), a),
        ]
        .map(|(e, u, v)| (format!("{prefix}:{e}"), [u.to_owned(), v.to_owned()]));
        let nodes = edges
            .iter()
            .flat_map(|(_, ends)| ends.iter().cloned())
            .collect();
        let mut graph = Graph::new(nodes, edges.into_iter().collect()).unwrap();
        for (vertex, order) in [
            (h.as_str(), ["ha", "hb", "hc"]),
            (a, ["ha", "ca", "ab"]),
            (b.as_str(), ["hb", "ab", "bc"]),
            (c.as_str(), ["hc", "bc", "ca"]),
        ] {
            graph.rotation.insert(
                vertex.to_owned(),
                order.map(|e| format!("{prefix}:{e}")).to_vec(),
            );
        }
        graph.wheel_hubs.insert(h);
        graph.protected_edges = graph.edges.key_set();
        graph.validate().unwrap();
        graph
    }

    #[test]
    fn first_free_corners_survive_block_joins() {
        let mut graph = wheel("root", "a");
        let mut corners = graph.free_corners().unwrap();
        assert!(!corners.contains_key("root:h"));
        for step in 0..24 {
            let prefix = format!("child:{step}");
            let at = if step % 3 == 0 { "a" } else { "root:b" };
            let child = if step % 2 == 0 {
                wheel(&prefix, at)
            } else {
                Graph::new(
                    [at.to_owned(), prefix.clone()].into_iter().collect(),
                    [(prefix.clone(), [at.to_owned(), prefix.clone()])]
                        .into_iter()
                        .collect(),
                )
                .unwrap()
            };
            let child_corners = child.free_corners().unwrap();
            graph = graph
                .join_vertex(&child, at, &corners[at], &child_corners[at])
                .unwrap();
            for (n, corner) in child_corners {
                corners.entry(n).or_insert(corner);
            }
            // Rewalking every face is the reference corner selection. Joining
            // at a cut vertex must retain the parent's first available edge.
            assert_eq!(corners, graph.free_corners().unwrap());
            graph.validate().unwrap();
        }
        let disjoint = wheel("disjoint", "other");
        corners.extend(disjoint.free_corners().unwrap());
        graph = graph.unite(&disjoint).unwrap();
        assert_eq!(corners, graph.free_corners().unwrap());
        graph.validate().unwrap();
    }
}
