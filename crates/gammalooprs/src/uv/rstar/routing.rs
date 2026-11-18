//! Opt-in routing of the root (boundary) momenta through an R* remainder
//! (root-routing lane, `GAMMALOOP_ROOT_ROUTING`).
//!
//! The remainder basis is [`LMBext::shrunken_sub_lmb`]: a depth-first
//! spanning tree of the remaining graph (insertions shrunk to points) rooted
//! at the node of the highest boundary half-edge. Every boundary momentum flows
//! along that tree to the root, so every tree edge on such a path carries a
//! root momentum in its `den` and is Taylor expanded by the outer Taylor (or
//! the IR hard Taylor in the soft edges). A depth-first tree makes those
//! paths long.
//!
//! With `GAMMALOOP_ROOT_ROUTING=shortest` the same `lmb_impl` is guided by a
//! spanning tree that contains a minimum Steiner tree of the boundary vertices
//! (for two boundary vertices: a shortest path). The boundary momenta then
//! flow only through the Steiner tree, i.e. through the minimal number of
//! propagators. `shortest-avoid-marked` / `shortest-prefer-marked` break ties
//! between minimal trees against / in favour of the marked (massive) edge.
//! Unset or `default`: the production basis, byte for byte.
//!
//! K R* is independent of the routing (a shift of loop momenta), so the
//! result is unchanged; only the size of the intermediate expressions is.
//!
//! `GAMMALOOP_ROOT_ROUTING_LOG=1|2` (diagnostics, no value depends on it):
//! `RR_BASIS` per remainder (propagators carrying a root momentum, against the
//! minimum), the stage timings `RR_CONTRACT`/`RR_OUTER`/`RR_PHASE`; with `2`
//! also the den-product Taylor shape `|D_j|` per remainder and, at the outer
//! Taylor, `RR_TAYLOR` with `|P_i|` and `|D_j|`.

use std::{
    collections::{BTreeMap, BinaryHeap},
    sync::LazyLock,
};

use eyre::{Result, eyre};
use linnet::half_edge::{
    involution::{EdgeIndex, HedgePair},
    subgraph::{Inclusion, InternalSubGraph, ModifySubSet, SuBitGraph, SubSetLike, SubSetOps},
};

use crate::graph::{Graph, LMBext, LoopMomentumBasis};

/// The routing rule for the boundary momenta of a remainder.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum RootRouting {
    Default,
    Shortest,
    AvoidMarked,
    PreferMarked,
    /// Keep the marked edge off the boundary path whenever the graph allows
    /// it, even through a longer path; shortest among those.
    StrictAvoidMarked,
    /// Spanning-tree enumeration: marked edge off the root path, then the
    /// smallest total routing (loop + root momentum terms over the dens).
    Compact,
    /// As `Compact`, keyed on the loop terms first.
    CompactLoops,
    /// As `Compact`, keyed on the root-carrying dens first.
    CompactRoots,
    /// Smallest total routing, the marked edge only as a tie-break.
    CompactFree,
}

impl RootRouting {
    pub(crate) fn as_str(self) -> &'static str {
        match self {
            RootRouting::Default => "default",
            RootRouting::Shortest => "shortest",
            RootRouting::AvoidMarked => "shortest-avoid-marked",
            RootRouting::PreferMarked => "shortest-prefer-marked",
            RootRouting::StrictAvoidMarked => "strict-avoid-marked",
            RootRouting::Compact => "compact",
            RootRouting::CompactLoops => "compact-loops",
            RootRouting::CompactRoots => "compact-roots",
            RootRouting::CompactFree => "compact-free",
        }
    }

    /// Edge weight: 1000 per propagator; the marked edge ±1 for tie-breaks.
    fn weight(self, edge: EdgeIndex, marked: Option<EdgeIndex>) -> u64 {
        match (self, Some(edge) == marked) {
            (RootRouting::AvoidMarked, true) => 1001,
            (RootRouting::PreferMarked, true) => 999,
            (RootRouting::StrictAvoidMarked, true) => 1_000_000,
            _ => 1000,
        }
    }
}

fn parse_mode(value: Option<&str>) -> std::result::Result<RootRouting, String> {
    match value {
        None | Some("" | "default" | "0" | "off") => Ok(RootRouting::Default),
        Some("shortest" | "1") => Ok(RootRouting::Shortest),
        Some("shortest-avoid-marked" | "avoid-marked") => Ok(RootRouting::AvoidMarked),
        Some("shortest-prefer-marked" | "prefer-marked") => Ok(RootRouting::PreferMarked),
        Some("strict-avoid-marked") => Ok(RootRouting::StrictAvoidMarked),
        Some("compact") => Ok(RootRouting::Compact),
        Some("compact-loops") => Ok(RootRouting::CompactLoops),
        Some("compact-roots") => Ok(RootRouting::CompactRoots),
        Some("compact-free") => Ok(RootRouting::CompactFree),
        Some(other) => Err(format!(
            "GAMMALOOP_ROOT_ROUTING={other:?}: expected default, shortest, shortest-avoid-marked, shortest-prefer-marked, strict-avoid-marked, compact, compact-loops, compact-roots or compact-free"
        )),
    }
}

static MODE: LazyLock<std::result::Result<RootRouting, String>> = LazyLock::new(|| {
    parse_mode(std::env::var("GAMMALOOP_ROOT_ROUTING").ok().as_deref())
});

/// The routing mode (read once); an invalid value is an error.
pub(crate) fn mode() -> Result<RootRouting> {
    MODE.clone().map_err(|e| eyre!(e))
}

/// `GAMMALOOP_ROOT_ROUTING_LOG` level (0 when unset).
pub(crate) fn log_level() -> u8 {
    static LEVEL: LazyLock<u8> = LazyLock::new(|| {
        std::env::var("GAMMALOOP_ROOT_ROUTING_LOG")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(0)
    });
    *LEVEL
}

/// A spanning-tree guide of the shrunken remainder.
struct Guide {
    /// The guide subgraph handed to `lmb_impl`: the remainder minus its
    /// non-tree propagators (self-loops and boundary hairs kept).
    guide: SuBitGraph,
    /// Propagators of the Steiner tree (they carry the boundary momenta).
    steiner_edges: Vec<EdgeIndex>,
    terminals: usize,
    exact: bool,
}

/// One propagator of the shrunken remainder between two (identified) nodes.
struct Arc {
    u: usize,
    v: usize,
    edge: EdgeIndex,
    pair: HedgePair,
    weight: u64,
}

struct UnionFind(Vec<usize>);

impl UnionFind {
    fn new(n: usize) -> Self {
        Self((0..n).collect())
    }
    fn find(&mut self, mut x: usize) -> usize {
        while self.0[x] != x {
            self.0[x] = self.0[self.0[x]];
            x = self.0[x];
        }
        x
    }
    fn union(&mut self, a: usize, b: usize) -> bool {
        let (a, b) = (self.find(a), self.find(b));
        if a == b {
            return false;
        }
        self.0[a] = b;
        true
    }
}

/// Minimum spanning tree of the nodes in `members` (a bitmask over node
/// positions, at most 64 nodes) using `arcs` sorted by (weight, edge);
/// `None` if they are not connected.
fn induced_mst(arcs: &[Arc], n: usize, members: &[bool]) -> Option<(u64, Vec<usize>)> {
    let mut uf = UnionFind::new(n);
    let mut chosen = Vec::new();
    let mut weight = 0;
    for (i, arc) in arcs.iter().enumerate() {
        if members[arc.u] && members[arc.v] && uf.union(arc.u, arc.v) {
            weight += arc.weight;
            chosen.push(i);
        }
    }
    let count = members.iter().filter(|m| **m).count();
    (chosen.len() + 1 == count).then_some((weight, chosen))
}

/// Shortest path (Dijkstra, ties by arc order) between two nodes: arc ids.
fn shortest_path(arcs: &[Arc], n: usize, from: usize, to: usize) -> Option<Vec<usize>> {
    let mut adjacency = vec![Vec::new(); n];
    for (i, arc) in arcs.iter().enumerate() {
        adjacency[arc.u].push((arc.v, i));
        adjacency[arc.v].push((arc.u, i));
    }
    let mut distance = vec![u64::MAX; n];
    let mut previous: Vec<Option<(usize, usize)>> = vec![None; n];
    let mut heap = BinaryHeap::new();
    distance[from] = 0;
    heap.push(std::cmp::Reverse((0_u64, from)));
    while let Some(std::cmp::Reverse((d, node))) = heap.pop() {
        if d > distance[node] {
            continue;
        }
        for &(next, arc) in &adjacency[node] {
            let candidate = d + arcs[arc].weight;
            if candidate < distance[next] {
                distance[next] = candidate;
                previous[next] = Some((node, arc));
                heap.push(std::cmp::Reverse((candidate, next)));
            }
        }
    }
    if distance[to] == u64::MAX {
        return None;
    }
    let mut path = Vec::new();
    let mut node = to;
    while let Some((prev, arc)) = previous[node] {
        path.push(arc);
        node = prev;
    }
    Some(path)
}

/// The non-terminal node count up to which the Steiner tree is exact
/// (enumeration of the Steiner node subsets).
const EXACT_STEINER_NODES: usize = 16;

/// A minimum Steiner tree of `terminals` (arc ids), exact for at most
/// [`EXACT_STEINER_NODES`] other nodes, else the shortest-path heuristic.
fn steiner_tree(arcs: &[Arc], n: usize, terminals: &[usize]) -> Option<(Vec<usize>, bool)> {
    if terminals.len() <= 1 {
        return Some((Vec::new(), true));
    }
    if terminals.len() == 2 {
        return shortest_path(arcs, n, terminals[0], terminals[1]).map(|p| (p, true));
    }
    let mut is_terminal = vec![false; n];
    for t in terminals {
        is_terminal[*t] = true;
    }
    let others: Vec<usize> = (0..n).filter(|v| !is_terminal[*v]).collect();
    if others.len() <= EXACT_STEINER_NODES {
        let mut best: Option<(u64, Vec<usize>)> = None;
        for mask in 0_u64..(1_u64 << others.len()) {
            let mut members = is_terminal.clone();
            for (bit, v) in others.iter().enumerate() {
                if mask >> bit & 1 == 1 {
                    members[*v] = true;
                }
            }
            if let Some((weight, chosen)) = induced_mst(arcs, n, &members)
                && best.as_ref().is_none_or(|(w, _)| weight < *w)
            {
                best = Some((weight, chosen));
            }
        }
        return best.map(|(_, chosen)| (chosen, true));
    }
    // Shortest-path heuristic: grow from the first terminal, adding the
    // nearest remaining terminal by a shortest path to the current tree.
    let mut in_tree = vec![false; n];
    in_tree[terminals[0]] = true;
    let mut chosen: Vec<usize> = Vec::new();
    let mut remaining: Vec<usize> = terminals[1..].to_vec();
    while !remaining.is_empty() {
        let mut best: Option<(u64, usize, Vec<usize>)> = None;
        for (index, t) in remaining.iter().enumerate() {
            for source in (0..n).filter(|v| in_tree[*v]) {
                if let Some(path) = shortest_path(arcs, n, source, *t) {
                    let w: u64 = path.iter().map(|a| arcs[*a].weight).sum();
                    if best.as_ref().is_none_or(|(bw, _, _)| w < *bw) {
                        best = Some((w, index, path));
                    }
                }
            }
        }
        let (_, index, path) = best?;
        for arc in path {
            in_tree[arcs[arc].u] = true;
            in_tree[arcs[arc].v] = true;
            if !chosen.contains(&arc) {
                chosen.push(arc);
            }
        }
        remaining.remove(index);
    }
    Some((chosen, false))
}

/// The guide of [`guided_shrunken_sub_lmb`] on the shrunken remainder, or
/// `None` when it does not apply (no boundary flow, disconnected remainder):
/// then the production basis is used.
fn guide<G>(
    contracted: &G,
    remainder: &SuBitGraph,
    externals: &SuBitGraph,
    weight: &dyn Fn(EdgeIndex) -> u64,
    preferred: &dyn Fn(EdgeIndex) -> bool,
) -> Option<Guide>
where
    G: GraphAccess,
{
    let mut nodes = BTreeMap::new();
    let node_of = |node: usize, nodes: &mut BTreeMap<usize, usize>| {
        let next = nodes.len();
        *nodes.entry(node).or_insert(next)
    };
    let mut arcs = Vec::new();
    for (pair, edge) in contracted.paired_edges(remainder) {
        let (a, b) = contracted.ends(&pair)?;
        let (u, v) = (node_of(a, &mut nodes), node_of(b, &mut nodes));
        if u == v {
            // a self-loop carries no boundary flow
            continue;
        }
        arcs.push(Arc { u, v, edge, pair, weight: weight(edge) });
    }
    let mut terminals = Vec::new();
    for hedge in externals.included_iter() {
        let t = node_of(contracted.hedge_node(hedge), &mut nodes);
        if !terminals.contains(&t) {
            terminals.push(t);
        }
    }
    let n = nodes.len();
    if terminals.len() <= 1 || n > 64 {
        return None;
    }
    arcs.sort_by_key(|arc| (arc.weight, usize::from(arc.edge)));
    let (steiner, exact) = steiner_tree(&arcs, n, &terminals)?;
    // Extend to a spanning tree: Steiner arcs, then the completion order of
    // GAMMALOOP_ROOT_ROUTING_TREE (keep: the production tree's arcs first;
    // kruskal: edge order; bfs: breadth first from the Steiner tree).
    let mut uf = UnionFind::new(n);
    let mut in_tree = vec![false; arcs.len()];
    for arc in &steiner {
        uf.union(arcs[*arc].u, arcs[*arc].v);
        in_tree[*arc] = true;
    }
    let order = completion_order(&arcs, n, &in_tree, preferred);
    for i in order {
        if !in_tree[i] && uf.union(arcs[i].u, arcs[i].v) {
            in_tree[i] = true;
        }
    }
    let tree_size = in_tree.iter().filter(|t| **t).count();
    if tree_size + 1 != n {
        return None;
    }
    let mut guide = remainder.clone();
    for (i, arc) in arcs.iter().enumerate() {
        if !in_tree[i] {
            guide.sub(arc.pair);
        }
    }
    let mut steiner_edges: Vec<EdgeIndex> = steiner.iter().map(|a| arcs[*a].edge).collect();
    steiner_edges.sort();
    Some(Guide {
        guide,
        steiner_edges,
        terminals: terminals.len(),
        exact,
    })
}

/// What [`guide`] needs of the shrunken graph.
trait GraphAccess {
    fn paired_edges(&self, subgraph: &SuBitGraph) -> Vec<(HedgePair, EdgeIndex)>;
    fn ends(&self, pair: &HedgePair) -> Option<(usize, usize)>;
    fn hedge_node(&self, hedge: linnet::half_edge::involution::Hedge) -> usize;
}

impl<E, V, H, N> GraphAccess for linnet::half_edge::HedgeGraph<E, V, H, N>
where
    N: linnet::half_edge::nodestore::NodeStorageOps<NodeData = V>,
{
    fn paired_edges(&self, subgraph: &SuBitGraph) -> Vec<(HedgePair, EdgeIndex)> {
        self.iter_edges_of(subgraph)
            .filter(|(pair, _, _)| pair.is_paired())
            .map(|(pair, edge, _)| (pair, edge))
            .collect()
    }
    fn ends(&self, pair: &HedgePair) -> Option<(usize, usize)> {
        match pair {
            HedgePair::Paired { source, sink } => {
                Some((self.node_id(*source).0, self.node_id(*sink).0))
            }
            _ => None,
        }
    }
    fn hedge_node(&self, hedge: linnet::half_edge::involution::Hedge) -> usize {
        self.node_id(hedge).0
    }
}

/// The boundary-flow report of one guided or measured basis.
#[derive(Clone, Debug, Default)]
pub(crate) struct SteinerReport {
    pub terminals: usize,
    pub min_edges: usize,
    pub exact: bool,
    pub contains_marked: bool,
}

/// [`LMBext::shrunken_sub_lmb`] with the forest guided by a spanning tree
/// that contains a minimum (weighted) Steiner tree of the boundary nodes.
/// `Ok(None)`: no guide applies; the caller uses the production basis.
pub(crate) fn guided_shrunken_sub_lmb(
    graph: &Graph,
    outer: &SuBitGraph,
    shrunken: &InternalSubGraph,
    externals: SuBitGraph,
    rule: RootRouting,
    marked: Option<EdgeIndex>,
) -> Result<Option<(LoopMomentumBasis, SteinerReport)>> {
    let g = &graph.underlying;
    let remainder = outer.subtract(&shrunken.filter);
    let contracted_externals = externals.subtract(&shrunken.filter);
    let mut contracted = g.to_ref();
    if !shrunken.is_empty() {
        for component in g.connected_components(shrunken) {
            let Some(root) = component.included_iter().next() else {
                continue;
            };
            let node_data = &g[g.node_id(root)];
            contracted.identify_nodes_of_subgraph_without_self_edges::<_, SuBitGraph>(
                &component, node_data,
            );
        }
    }
    let weight = |edge: EdgeIndex| rule.weight(edge, marked);
    let production = graph.shrunken_sub_lmb(outer, shrunken, externals.clone()).ok();
    let preferred = |edge: EdgeIndex| {
        production.as_ref().is_some_and(|lmb| {
            let (_, pair) = &g[&edge];
            matches!(pair, HedgePair::Paired { source, .. } if lmb.tree.includes(source))
        })
    };
    let found = if rule.enumerates() {
        enumerate_guide(&contracted, &remainder, &contracted_externals, marked, rule)
    } else {
        guide(&contracted, &remainder, &contracted_externals, &weight, &preferred)
    };
    let Some(guide) = found else {
        return Ok(None);
    };
    let mut lmb = contracted
        .lmb_impl(&remainder, &guide.guide, contracted_externals)
        .map_err(|e| eyre!("guided root routing: {e}"))?;
    graph.canonicalize_lmb_external_order(&mut lmb);
    let report = SteinerReport {
        terminals: guide.terminals,
        min_edges: guide.steiner_edges.len(),
        exact: guide.exact,
        contains_marked: marked.is_some_and(|m| guide.steiner_edges.contains(&m)),
    };
    Ok(Some((lmb, report)))
}

/// The minimum Steiner tree of the boundary nodes of the shrunken remainder
/// under `rule`'s weights, without building a basis (census).
pub(crate) fn minimal_report(
    graph: &Graph,
    outer: &SuBitGraph,
    shrunken: &InternalSubGraph,
    externals: &SuBitGraph,
    rule: RootRouting,
    marked: Option<EdgeIndex>,
) -> Option<SteinerReport> {
    let g = &graph.underlying;
    let remainder = outer.subtract(&shrunken.filter);
    let contracted_externals = externals.subtract(&shrunken.filter);
    let mut contracted = g.to_ref();
    if !shrunken.is_empty() {
        for component in g.connected_components(shrunken) {
            let Some(root) = component.included_iter().next() else {
                continue;
            };
            let node_data = &g[g.node_id(root)];
            contracted.identify_nodes_of_subgraph_without_self_edges::<_, SuBitGraph>(
                &component, node_data,
            );
        }
    }
    let rule = if rule == RootRouting::Default { RootRouting::Shortest } else { rule };
    let weight = |edge: EdgeIndex| rule.weight(edge, marked);
    let guide = guide(&contracted, &remainder, &contracted_externals, &weight, &|_| false)?;
    Some(SteinerReport {
        terminals: guide.terminals,
        min_edges: guide.steiner_edges.len(),
        exact: guide.exact,
        contains_marked: marked.is_some_and(|m| guide.steiner_edges.contains(&m)),
    })
}

/// The remainder basis under the current [`mode`]: the production
/// `shrunken_sub_lmb`, or the guided one.
pub(crate) fn remainder_basis(
    graph: &Graph,
    outer: &SuBitGraph,
    shrunken: &InternalSubGraph,
    externals: SuBitGraph,
    marked: Option<EdgeIndex>,
) -> Result<LoopMomentumBasis> {
    let rule = mode()?;
    if rule != RootRouting::Default
        && let Some((lmb, _)) =
            guided_shrunken_sub_lmb(graph, outer, shrunken, externals.clone(), rule, marked)?
    {
        return Ok(lmb);
    }
    graph
        .shrunken_sub_lmb(outer, shrunken, externals)
        .map_err(|e| eyre!("{e}"))
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn modes_parse() {
        assert_eq!(parse_mode(None), Ok(RootRouting::Default));
        assert_eq!(parse_mode(Some("default")), Ok(RootRouting::Default));
        assert_eq!(parse_mode(Some("shortest")), Ok(RootRouting::Shortest));
        assert_eq!(
            parse_mode(Some("shortest-avoid-marked")),
            Ok(RootRouting::AvoidMarked)
        );
        assert!(parse_mode(Some("bogus")).is_err());
    }

    fn arc(u: usize, v: usize, e: usize, w: u64) -> Arc {
        Arc {
            u,
            v,
            edge: EdgeIndex::from(e),
            pair: HedgePair::Unpaired {
                hedge: linnet::half_edge::involution::Hedge(0),
                flow: linnet::half_edge::involution::Flow::Source,
            },
            weight: w,
        }
    }

    #[test]
    fn steiner_star_and_path() {
        // A 4-cycle 0-1-2-3-0 plus centre 4 joined to 0,1,2,3 (spokes weight 1000,
        // rim weight 1000): terminals 0,1,2 -> 2 edges on the rim (0-1, 1-2) or
        // 3 spokes; min = 2 edges.
        let mut arcs = vec![
            arc(0, 1, 0, 1000),
            arc(1, 2, 1, 1000),
            arc(2, 3, 2, 1000),
            arc(3, 0, 3, 1000),
            arc(0, 4, 4, 1000),
            arc(1, 4, 5, 1000),
            arc(2, 4, 6, 1000),
            arc(3, 4, 7, 1000),
        ];
        arcs.sort_by_key(|a| (a.weight, usize::from(a.edge)));
        let (tree, exact) = steiner_tree(&arcs, 5, &[0, 1, 2]).unwrap();
        assert!(exact);
        assert_eq!(tree.len(), 2);
        // terminals 0 and 2: shortest path has 2 edges
        let path = shortest_path(&arcs, 5, 0, 2).unwrap();
        assert_eq!(path.len(), 2);
        // terminals 0,2 with 1 and 3 opposite: the opposite corners 1 and 3 need
        // 0-1? no: terminals 1,3 -> path 1-0-3 / 1-2-3 / 1-4-3 -> 2 edges
        let (tree, _) = steiner_tree(&arcs, 5, &[0, 1, 2, 3]).unwrap();
        assert_eq!(tree.len(), 3);
    }
}

fn root_carrying(basis: &LoopMomentumBasis, edge: EdgeIndex) -> bool {
    !basis
        .ext_atom::<symbolica::atom::Atom>(edge, crate::utils::GS.emr_mom, &[], true)
        .is_zero()
}

/// The paired (propagator) edges of `region.remaining` and how many of them
/// carry a boundary momentum in `basis`; whether `marked` is among them and
/// carries one.
fn carriers(
    graph: &Graph,
    region: &super::GraphRegion,
    basis: &LoopMomentumBasis,
    marked: Option<EdgeIndex>,
) -> (usize, usize, bool, bool) {
    let (mut dens, mut root_dens, mut marked_present, mut marked_root) = (0, 0, false, false);
    for (pair, edge, _) in graph.iter_edges_of(&region.remaining) {
        if !pair.is_paired() {
            continue;
        }
        dens += 1;
        let carries = root_carrying(basis, edge);
        root_dens += usize::from(carries);
        if Some(edge) == marked {
            marked_present = true;
            marked_root = carries;
        }
    }
    (dens, root_dens, marked_present, marked_root)
}

/// `RR_BASIS` (log level >= 1): per remainder construction, the propagators
/// carrying a boundary momentum in `basis` against the minimum; level >= 2
/// adds the den-product Taylor shape `|D_0..2|` (boundary edges scaled, the
/// reference mass compensated as in the outer Taylor).
pub(crate) fn log_basis(
    graph: &Graph,
    region: &super::GraphRegion,
    contracted: &InternalSubGraph,
    basis: &LoopMomentumBasis,
    marked: EdgeIndex,
    reference_mass_squared: &symbolica::atom::Atom,
    routed: &symbolica::atom::Atom,
) {
    use symbolica::atom::{Atom, AtomCore, AtomView};
    let rule = mode().unwrap_or(RootRouting::Default);
    let (dens, root_dens, marked_present, marked_root) =
        carriers(graph, region, basis, Some(marked));
    let mut boundary_edges: Vec<EdgeIndex> = graph
        .iter_edges_of(&region.boundary)
        .map(|(_, edge, _)| edge)
        .collect();
    boundary_edges.sort();
    boundary_edges.dedup();
    let report = |rule| {
        minimal_report(graph, &region.support, contracted, &region.boundary, rule, Some(marked))
    };
    let shortest = report(RootRouting::Shortest);
    let avoid = report(RootRouting::AvoidMarked);
    let prefer = report(RootRouting::PreferMarked);
    let paired = |edges: &mut dyn Iterator<Item = bool>| edges.filter(|p| *p).count();
    let (loop_terms, root_terms) = routing_terms(graph, region, basis);
    let root_scope = paired(&mut graph.iter_edges_of(&region.support).map(|(p, _, _)| p.is_paired()))
        == paired(&mut graph.underlying.iter_edges().map(|(p, _, _)| p.is_paired()));
    let mut line = format!(
        "RR_BASIS graph={} mode={} root_scope={root_scope} loops={} dens={dens} boundary_edges={} terminals={} root_dens={root_dens} min_root_dens={} exact={} marked_present={marked_present} marked_root={marked_root} min_marked_avoidable={} min_marked_possible={} loop_terms={loop_terms} root_terms={root_terms}",
        graph.name,
        rule.as_str(),
        region.loops,
        boundary_edges.len(),
        shortest.as_ref().map_or(0, |r| r.terminals),
        shortest
            .as_ref()
            .map_or("-".to_string(), |r| r.min_edges.to_string()),
        shortest.as_ref().is_none_or(|r| r.exact),
        avoid.as_ref().is_none_or(|r| !r.contains_marked),
        prefer.as_ref().is_some_and(|r| r.contains_marked),
    );
    if log_level() >= 2 {
        let parameter = symbolica::symbol!("rr_census_t"; Scalar);
        let auxiliary = Atom::var(symbolica::symbol!("rr_census_m2"));
        let den_factors: Vec<AtomView<'_>> = match routed.as_view() {
            AtomView::Mul(product) => product
                .iter()
                .filter(|factor| {
                    matches!(factor, AtomView::Pow(power)
                        if matches!(power.get_base_exp().0,
                            AtomView::Fun(f) if f.get_symbol() == crate::utils::GS.den))
                })
                .collect(),
            _ => Vec::new(),
        };
        if !den_factors.is_empty() {
            let product = Atom::mul_many(den_factors.iter().copied());
            let scaled = super::expression::RStarIntegrand::scale_edges_view(
                product.as_view(),
                &boundary_edges,
                parameter,
            );
            let compensated = if reference_mass_squared.is_zero() {
                scaled
            } else {
                scaled
                    .replace(reference_mass_squared.to_pattern())
                    .with((Atom::num(1) - Atom::var(parameter).pow(2)) * &auxiliary)
            };
            let shape =
                crate::uv::approx::counts::taylor_shape(compensated.as_view(), parameter, 2);
            line.push_str(&format!(
                " d_orders={} den_patterns={}",
                shape["d_orders"], shape["den_patterns"]
            ));
        }
    }
    eprintln!("{line}");
}

/// `RR_TAYLOR` (log level >= 2): at the outer Taylor, the propagators that
/// carry a root momentum and the `P_i x D_j` shape of the product.
pub(crate) fn log_taylor(
    graph: &Graph,
    region: &super::GraphRegion,
    basis: &LoopMomentumBasis,
    root_edges: &[EdgeIndex],
    degree: usize,
    shape: &serde_json::Value,
) {
    let (dens, root_dens, _, _) = carriers(graph, region, basis, None);
    eprintln!(
        "RR_TAYLOR graph={} mode={} loops={} degree={degree} roots={} dens={dens} root_dens={root_dens} p_orders={} d_orders={} den_patterns={} pair_estimate={} p_summands={} p_bytes={}",
        graph.name,
        mode().unwrap_or(RootRouting::Default).as_str(),
        region.loops,
        root_edges.len(),
        shape["p_orders"],
        shape["d_orders"],
        shape["den_patterns"],
        shape["pair_estimate"],
        shape["p_summands"],
        shape["p_bytes"],
    );
}

/// Σ over the remaining propagators of the number of loop momenta and of
/// boundary momenta in their routed momentum (the size of the routing).
fn routing_terms(graph: &Graph, region: &super::GraphRegion, basis: &LoopMomentumBasis) -> (usize, usize) {
    use symbolica::atom::Atom;
    let (mut loops, mut roots) = (0, 0);
    for (pair, edge, _) in graph.iter_edges_of(&region.remaining) {
        if !pair.is_paired() {
            continue;
        }
        let l = basis.loop_atom::<Atom>(edge, crate::utils::GS.emr_mom, &[], true);
        let e = basis.ext_atom::<Atom>(edge, crate::utils::GS.emr_mom, &[], true);
        loops += vakint::census::summands(l.as_view());
        roots += vakint::census::summands(e.as_view());
    }
    (loops, roots)
}

/// How a Steiner tree is completed to a spanning tree.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Completion {
    Keep,
    Kruskal,
    Bfs,
}

fn completion() -> Completion {
    static MODE: LazyLock<Completion> = LazyLock::new(|| {
        match std::env::var("GAMMALOOP_ROOT_ROUTING_TREE").as_deref() {
            Ok("kruskal") => Completion::Kruskal,
            Ok("bfs") => Completion::Bfs,
            _ => Completion::Keep,
        }
    });
    *MODE
}

/// The order in which the non-Steiner arcs are offered to the spanning tree.
fn completion_order(
    arcs: &[Arc],
    n: usize,
    in_tree: &[bool],
    preferred: &dyn Fn(EdgeIndex) -> bool,
) -> Vec<usize> {
    match completion() {
        Completion::Kruskal => (0..arcs.len()).collect(),
        Completion::Keep => {
            let (mut first, mut rest): (Vec<usize>, Vec<usize>) =
                (0..arcs.len()).partition(|i| preferred(arcs[*i].edge));
            first.sort_by_key(|i| usize::from(arcs[*i].edge));
            rest.sort_by_key(|i| usize::from(arcs[*i].edge));
            first.extend(rest);
            first
        }
        Completion::Bfs => {
            let mut reached = vec![false; n];
            let mut any = false;
            for (i, arc) in arcs.iter().enumerate() {
                if in_tree[i] {
                    reached[arc.u] = true;
                    reached[arc.v] = true;
                    any = true;
                }
            }
            if !any && n > 0 {
                reached[0] = true;
            }
            let mut order = Vec::new();
            let mut used = in_tree.to_vec();
            loop {
                let frontier: Vec<usize> = (0..arcs.len())
                    .filter(|i| !used[*i] && (reached[arcs[*i].u] != reached[arcs[*i].v]))
                    .collect();
                if frontier.is_empty() {
                    break;
                }
                for i in frontier {
                    if !(reached[arcs[i].u] && reached[arcs[i].v]) {
                        reached[arcs[i].u] = true;
                        reached[arcs[i].v] = true;
                        order.push(i);
                    }
                    used[i] = true;
                }
            }
            order.extend((0..arcs.len()).filter(|i| !used[*i]));
            order
        }
    }
}

impl RootRouting {
    fn enumerates(self) -> bool {
        matches!(
            self,
            RootRouting::Compact
                | RootRouting::CompactLoops
                | RootRouting::CompactRoots
                | RootRouting::CompactFree
        )
    }
}

/// The maximum number of spanning trees [`enumerate_guide`] evaluates; the
/// best one seen so far is used beyond it.
const MAX_TREES: usize = 500_000;

/// The routing size of one spanning tree (arc ids `tree`) rooted at `root`.
struct TreeCost {
    marked_root: bool,
    root_dens: usize,
    root_terms: usize,
    loop_terms: usize,
}

fn tree_cost(
    arcs: &[Arc],
    n: usize,
    tree: &[usize],
    root: usize,
    sources: &[usize],
    marked: Option<EdgeIndex>,
) -> TreeCost {
    let mut adjacency = vec![Vec::new(); n];
    for &a in tree {
        adjacency[arcs[a].u].push((arcs[a].v, a));
        adjacency[arcs[a].v].push((arcs[a].u, a));
    }
    let mut parent: Vec<Option<(usize, usize)>> = vec![None; n];
    let mut depth = vec![0usize; n];
    let mut seen = vec![false; n];
    let mut queue = std::collections::VecDeque::from([root]);
    seen[root] = true;
    while let Some(x) = queue.pop_front() {
        for &(y, a) in &adjacency[x] {
            if !seen[y] {
                seen[y] = true;
                parent[y] = Some((x, a));
                depth[y] = depth[x] + 1;
                queue.push_back(y);
            }
        }
    }
    let mut on_path = vec![false; arcs.len()];
    let mut root_terms = 0;
    for &s in sources {
        let mut x = s;
        while let Some((p, a)) = parent[x] {
            on_path[a] = true;
            root_terms += 1;
            x = p;
        }
    }
    let mut in_tree = vec![false; arcs.len()];
    for &a in tree {
        in_tree[a] = true;
    }
    let mut loop_terms = 0;
    for (i, arc) in arcs.iter().enumerate() {
        if in_tree[i] {
            continue;
        }
        let (mut x, mut y) = (arc.u, arc.v);
        let mut length = 1;
        while x != y {
            if depth[x] >= depth[y] {
                x = parent[x].map_or(x, |(p, _)| p);
            } else {
                y = parent[y].map_or(y, |(p, _)| p);
            }
            length += 1;
        }
        loop_terms += length;
    }
    TreeCost {
        marked_root: (0..arcs.len()).any(|i| on_path[i] && Some(arcs[i].edge) == marked),
        root_dens: on_path.iter().filter(|p| **p).count(),
        root_terms,
        loop_terms,
    }
}

/// The guide of the best spanning tree under `rule`'s key, by enumeration.
fn enumerate_guide<G: GraphAccess>(
    contracted: &G,
    remainder: &SuBitGraph,
    externals: &SuBitGraph,
    marked: Option<EdgeIndex>,
    rule: RootRouting,
) -> Option<Guide> {
    let mut nodes = BTreeMap::new();
    let node_of = |node: usize, nodes: &mut BTreeMap<usize, usize>| {
        let next = nodes.len();
        *nodes.entry(node).or_insert(next)
    };
    let mut arcs = Vec::new();
    for (pair, edge) in contracted.paired_edges(remainder) {
        let (a, b) = contracted.ends(&pair)?;
        let (u, v) = (node_of(a, &mut nodes), node_of(b, &mut nodes));
        if u != v {
            arcs.push(Arc { u, v, edge, pair, weight: 1000 });
        }
    }
    let hedges: Vec<_> = externals.included_iter().collect();
    let root_hedge = *hedges.last()?;
    let root = node_of(contracted.hedge_node(root_hedge), &mut nodes);
    let sources: Vec<usize> = hedges[..hedges.len() - 1]
        .iter()
        .map(|h| node_of(contracted.hedge_node(*h), &mut nodes))
        .collect();
    let n = nodes.len();
    let mut terminals: Vec<usize> = sources.clone();
    terminals.push(root);
    terminals.sort();
    terminals.dedup();
    if terminals.len() <= 1 || n < 2 {
        return None;
    }
    arcs.sort_by_key(|arc| usize::from(arc.edge));
    let key = |c: &TreeCost| -> (usize, usize, usize, usize) {
        let m = usize::from(c.marked_root);
        let total = c.loop_terms + c.root_terms;
        match rule {
            RootRouting::CompactLoops => (m, c.loop_terms, c.root_terms, c.root_dens),
            RootRouting::CompactRoots => (m, c.root_dens, c.loop_terms, c.root_terms),
            RootRouting::CompactFree => (total, m, c.root_dens, c.loop_terms),
            _ => (m, total, c.root_dens, c.loop_terms),
        }
    };
    let mut best: Option<((usize, usize, usize, usize), Vec<usize>, usize)> = None;
    let mut count = 0usize;
    // Backtracking over arcs: include (if acyclic) or exclude.
    fn walk(
        i: usize,
        chosen: &mut Vec<usize>,
        uf: &UnionFind,
        arcs: &[Arc],
        n: usize,
        count: &mut usize,
        visit: &mut dyn FnMut(&[usize]),
    ) {
        if *count >= MAX_TREES {
            return;
        }
        if chosen.len() + 1 == n {
            *count += 1;
            visit(chosen);
            return;
        }
        if i == arcs.len() || chosen.len() + (arcs.len() - i) + 1 < n {
            return;
        }
        let mut with = UnionFind(uf.0.clone());
        if with.union(arcs[i].u, arcs[i].v) {
            chosen.push(i);
            walk(i + 1, chosen, &with, arcs, n, count, visit);
            chosen.pop();
        }
        walk(i + 1, chosen, uf, arcs, n, count, visit);
    }
    let mut visit = |tree: &[usize]| {
        let cost = tree_cost(&arcs, n, tree, root, &sources, marked);
        let k = key(&cost);
        if best.as_ref().is_none_or(|(b, _, _)| k < *b) {
            best = Some((k, tree.to_vec(), cost.root_dens));
        }
    };
    walk(0, &mut Vec::new(), &UnionFind::new(n), &arcs, n, &mut count, &mut visit);
    let (_, tree, root_dens) = best?;
    let mut in_tree = vec![false; arcs.len()];
    for a in &tree {
        in_tree[*a] = true;
    }
    let mut guide = remainder.clone();
    for (i, arc) in arcs.iter().enumerate() {
        if !in_tree[i] {
            guide.sub(arc.pair);
        }
    }
    let _ = root_dens;
    Some(Guide {
        guide,
        steiner_edges: Vec::new(),
        terminals: terminals.len(),
        exact: count < MAX_TREES,
    })
}
