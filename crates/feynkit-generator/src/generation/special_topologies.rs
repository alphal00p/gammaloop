//! Self-energy, tadpole and snail vetoes and the factorized-loop count.
//!
//! The conditions are defined on the blocks (biconnected components) of the graph in which every
//! external leg ends on a degree-one node, so they do not depend on how the graph is labelled.
//! Removing the edges of a block splits the graph into one part per block vertex. The
//! momentum-carrying *core* is formed by the blocks whose parts hold external legs at two or more
//! of their vertices; with a single leg, whose momentum is conserved everywhere, by every block
//! but the self-loops; and by nothing in vacuum graphs. The other cyclic blocks are *vacuum
//! cycles*:
//!
//! - a vacuum cycle touching a core vertex or another cycle is a *snail*, and one reached from
//!   such a vertex through bridges a *tadpole*. Both are classified by the line they sit on: the
//!   core lines at that vertex, or else the lines of the block leading from it to the core,
//!   massive if one of them is massive. A vacuum graph with a loop is a tadpole on a massive
//!   line;
//! - a *self-energy insertion* is a core bridge with a core cycle in a part holding exactly one
//!   external leg, the others being in the other part. With two legs, both parts need a core
//!   cycle and both legs the vetoed mass, since each leg sees the other part as the insertion;
//! - the *factorized loop topologies* are the cyclic blocks through which at most one external
//!   momentum flows: their legs are in at most one part, or in two parts one of which holds a
//!   single leg.
//!
//! Edge masses may be unknown (untyped topologies); a mass-dependent veto then applies only when
//! it vetoes both masses, so an unknown mass never vetoes a graph that some typing would keep.

use crate::{Range, Veto};

/// The special-topology filters of a generation request.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct SpecialTopologyVetoes {
    pub(crate) self_energy: Option<Veto>,
    pub(crate) tadpoles: Option<Veto>,
    pub(crate) snails: Option<Veto>,
    pub(crate) factorized_loop_range: Range,
}

/// Whether vetoes for massive and for massless lines apply to a line of the given
/// masslessness.
fn applies(massless: Option<bool>, veto_massive: bool, veto_massless: bool) -> bool {
    match massless {
        Some(true) => veto_massless,
        Some(false) => veto_massive,
        None => veto_massive && veto_massless,
    }
}

impl SpecialTopologyVetoes {
    pub(crate) fn is_empty(&self) -> bool {
        self.self_energy.is_none()
            && self.tadpoles.is_none()
            && self.snails.is_none()
            && self.factorized_loop_range == Range::Unrestricted
    }

    pub(crate) fn vetoes(&self, graph: &PlainGraph) -> bool {
        if self.is_empty() {
            return false;
        }
        let blocks = Blocks::new(graph);
        if !self
            .factorized_loop_range
            .contains(blocks.factorized_loops())
        {
            return true;
        }

        if let Some(options) = self.self_energy
            && blocks.self_energy_legs(graph).iter().any(|legs| {
                legs.iter()
                    .all(|&massless| applies(Some(massless), options.massive, options.massless))
            })
        {
            return true;
        }
        blocks
            .vacuum_attachments(graph)
            .into_iter()
            .any(|attachment| {
                let options = match attachment.kind {
                    Attachment::Tadpole => self
                        .tadpoles
                        .map(|options| (options.massive, options.massless)),
                    Attachment::Snail => self
                        .snails
                        .map(|options| (options.massive, options.massless)),
                };
                options.is_some_and(|(veto_massive, veto_massless)| {
                    applies(attachment.massless, veto_massive, veto_massless)
                })
            })
    }
}

/// A multigraph whose external legs are degree-one nodes.
#[derive(Debug, Default)]
pub(crate) struct PlainGraph {
    edges: Vec<PlainEdge>,
    /// Incident edges of each node; a self-loop is listed once.
    edges_of: Vec<Vec<usize>>,
    /// The masslessness of the particle of each leg node.
    legs: Vec<Option<bool>>,
}

#[derive(Debug, Clone, Copy)]
struct PlainEdge {
    ends: (usize, usize),
    /// `None` when the edge is not typed yet.
    massless: Option<bool>,
}

impl PlainEdge {
    fn other_end(&self, node: usize) -> usize {
        if self.ends.0 == node {
            self.ends.1
        } else {
            self.ends.0
        }
    }

    fn touches(&self, node: usize) -> bool {
        self.ends.0 == node || self.ends.1 == node
    }
}

impl PlainGraph {
    pub(crate) fn add_vertex(&mut self) -> usize {
        self.edges_of.push(vec![]);
        self.legs.push(None);
        self.edges_of.len() - 1
    }

    pub(crate) fn add_edge(&mut self, ends: (usize, usize), massless: Option<bool>) {
        let edge = self.edges.len();
        self.edges.push(PlainEdge { ends, massless });
        self.edges_of[ends.0].push(edge);
        if ends.1 != ends.0 {
            self.edges_of[ends.1].push(edge);
        }
    }

    /// Attach an external leg of a particle with the given masslessness to `vertex`.
    pub(crate) fn add_leg(&mut self, vertex: usize, massless: bool) {
        let node = self.add_vertex();
        self.legs[node] = Some(massless);
        self.add_edge((node, vertex), Some(massless));
    }

    fn n_nodes(&self) -> usize {
        self.edges_of.len()
    }

    /// The edge sets of the blocks: Hopcroft-Tarjan with an edge stack, distinguishing
    /// parallel edges by the tree edge that entered a node. Every self-loop is its own block.
    fn blocks(&self) -> Vec<Vec<usize>> {
        let n = self.n_nodes();
        let mut blocks = vec![];
        let mut discovery = vec![usize::MAX; n];
        let mut low = vec![0; n];
        let mut time = 0;
        let mut edge_stack = vec![];
        for root in 0..n {
            if discovery[root] != usize::MAX {
                continue;
            }
            discovery[root] = time;
            low[root] = time;
            time += 1;
            // (node, entering edge, next incident position)
            let mut stack = vec![(root, usize::MAX, 0)];
            while let Some((node, entering, next)) =
                stack.last_mut().map(|top| (top.0, top.1, &mut top.2))
            {
                if let Some(&edge) = self.edges_of[node].get(*next) {
                    *next += 1;
                    let other = self.edges[edge].other_end(node);
                    if edge == entering || other == node {
                        continue;
                    }
                    if discovery[other] == usize::MAX {
                        edge_stack.push(edge);
                        discovery[other] = time;
                        low[other] = time;
                        time += 1;
                        stack.push((other, edge, 0));
                    } else if discovery[other] < discovery[node] {
                        edge_stack.push(edge);
                        low[node] = low[node].min(discovery[other]);
                    }
                    continue;
                }
                stack.pop();
                let Some(&(parent, _, _)) = stack.last() else {
                    continue;
                };
                low[parent] = low[parent].min(low[node]);
                if low[node] >= discovery[parent] {
                    let mut block = vec![];
                    while let Some(edge) = edge_stack.pop() {
                        block.push(edge);
                        if edge == entering {
                            break;
                        }
                    }
                    blocks.push(block);
                }
            }
        }
        blocks.extend(
            (0..self.edges.len())
                .filter(|edge| self.edges[*edge].ends.0 == self.edges[*edge].ends.1)
                .map(|edge| vec![edge]),
        );
        blocks
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Attachment {
    Tadpole,
    Snail,
}

/// A vacuum cycle hanging off a line.
#[derive(Debug, Clone, Copy)]
struct VacuumAttachment {
    kind: Attachment,
    /// The masslessness of the line it sits on.
    massless: Option<bool>,
}

/// Massive if one of the lines is massive, unknown if one is untyped.
fn masslessness(lines: impl IntoIterator<Item = Option<bool>>) -> Option<bool> {
    let mut massless = Some(true);
    for line in lines {
        match line {
            Some(false) => return Some(false),
            None => massless = None,
            Some(true) => {}
        }
    }
    massless
}

#[derive(Debug, Clone)]
struct Block {
    edges: Vec<usize>,
    vertices: Vec<usize>,
    cyclic: bool,
    self_loop: bool,
    /// The component of every node in the graph without the block's edges.
    component: Vec<usize>,
    /// The leg nodes of each component.
    legs: Vec<Vec<usize>>,
    core: bool,
}

impl Block {
    fn new(graph: &PlainGraph, edges: Vec<usize>) -> Self {
        let n = graph.n_nodes();
        let mut root: Vec<usize> = (0..n).collect();
        fn find(root: &mut [usize], node: usize) -> usize {
            let mut top = node;
            while root[top] != top {
                top = root[top];
            }
            let mut node = node;
            while root[node] != top {
                let next = root[node];
                root[node] = top;
                node = next;
            }
            top
        }
        let mut in_block = vec![false; graph.edges.len()];
        for &edge in &edges {
            in_block[edge] = true;
        }
        for (edge, data) in graph.edges.iter().enumerate() {
            if !in_block[edge] {
                let (a, b) = (find(&mut root, data.ends.0), find(&mut root, data.ends.1));
                root[a] = b;
            }
        }
        let component: Vec<usize> = (0..n).map(|node| find(&mut root, node)).collect();
        let mut legs = vec![vec![]; n];
        for node in 0..n {
            if graph.legs[node].is_some() {
                legs[component[node]].push(node);
            }
        }
        let mut vertices: Vec<usize> = edges
            .iter()
            .flat_map(|edge| [graph.edges[*edge].ends.0, graph.edges[*edge].ends.1])
            .collect();
        vertices.sort_unstable();
        vertices.dedup();
        let self_loop = vertices.len() == 1;
        Self {
            cyclic: edges.len() > 1 || self_loop,
            self_loop,
            edges,
            vertices,
            component,
            legs,
            core: false,
        }
    }

    /// The leg nodes in the part of the graph hanging at `vertex`, `vertex` included.
    fn legs_at(&self, vertex: usize) -> &[usize] {
        &self.legs[self.component[vertex]]
    }

    /// The numbers of legs of the parts at the vertices holding some.
    fn leg_sides(&self) -> impl Iterator<Item = usize> + '_ {
        self.vertices
            .iter()
            .map(|vertex| self.legs_at(*vertex).len())
            .filter(|legs| *legs > 0)
    }

    fn lines_at<'a>(
        &'a self,
        graph: &'a PlainGraph,
        vertex: usize,
    ) -> impl Iterator<Item = Option<bool>> + 'a {
        self.edges
            .iter()
            .map(|edge| graph.edges[*edge])
            .filter(move |edge| edge.touches(vertex))
            .map(|edge| edge.massless)
    }
}

/// The blocks of a [`PlainGraph`], the core marked.
struct Blocks {
    blocks: Vec<Block>,
    total_legs: usize,
}

impl Blocks {
    fn new(graph: &PlainGraph) -> Self {
        let total_legs = graph.legs.iter().flatten().count();
        let mut blocks: Vec<Block> = graph
            .blocks()
            .into_iter()
            .map(|edges| Block::new(graph, edges))
            .collect();
        for block in &mut blocks {
            // a single leg's momentum is conserved through every block, self-loops aside
            block.core = if total_legs == 1 {
                !block.self_loop
            } else {
                block.leg_sides().count() >= 2
            };
        }
        Self { blocks, total_legs }
    }

    /// The other blocks in the part of `block` hanging at `vertex`.
    fn blocks_at<'a>(&'a self, block: &'a Block, vertex: usize) -> impl Iterator<Item = &'a Block> {
        let part = block.component[vertex];
        self.blocks.iter().filter(move |other| {
            !std::ptr::eq(*other, block) && block.component[other.vertices[0]] == part
        })
    }

    /// The blocks containing `vertex`.
    fn containing(&self, vertex: usize) -> impl Iterator<Item = &Block> {
        self.blocks
            .iter()
            .filter(move |block| block.vertices.contains(&vertex))
    }

    /// Cyclic blocks through which at most one external momentum flows.
    fn factorized_loops(&self) -> usize {
        self.blocks
            .iter()
            .filter(|block| block.cyclic)
            .filter(|block| {
                let mut sides = block.leg_sides();
                match (sides.next(), sides.next(), sides.next()) {
                    (_, None, _) => true,
                    (Some(a), Some(b), None) => a == 1 || b == 1,
                    _ => false,
                }
            })
            .count()
    }

    /// For every self-energy insertion, the masslessness of the legs whose mass decides it.
    fn self_energy_legs(&self, graph: &PlainGraph) -> Vec<Vec<bool>> {
        if self.total_legs < 2 {
            return vec![];
        }
        let massless = |legs: &[usize]| -> Vec<bool> {
            legs.iter()
                .map(|leg| graph.legs[*leg].expect("a leg node"))
                .collect()
        };
        let mut insertions = vec![];
        for block in &self.blocks {
            let [left, right] = block.vertices[..] else {
                continue;
            };
            let sides = [left, right].map(|vertex| block.legs_at(vertex));
            // legs are not insertions, and the insertion must carry the leg momentum
            if block.cyclic
                || !block.core
                || [left, right]
                    .iter()
                    .any(|vertex| graph.legs[*vertex].is_some())
            {
                continue;
            }
            let has_loop = |vertex: usize, core: bool| {
                self.blocks_at(block, vertex)
                    .any(|other| other.cyclic && (other.core || !core))
            };
            if self.total_legs == 2 {
                if (has_loop(left, true) && has_loop(right, false))
                    || (has_loop(right, true) && has_loop(left, false))
                {
                    insertions.push([massless(sides[0]), massless(sides[1])].concat());
                }
                continue;
            }
            for (vertex, legs) in [left, right].into_iter().zip(sides) {
                if legs.len() == 1 && has_loop(vertex, true) {
                    insertions.push(massless(legs));
                }
            }
        }
        insertions
    }

    /// The block from `vertex` towards the core, other than `from`, and the masslessness of
    /// its lines at `vertex`; at a core vertex, the core lines.
    fn line_towards_core(&self, graph: &PlainGraph, vertex: usize, from: &Block) -> Option<bool> {
        let core_lines: Vec<_> = self
            .containing(vertex)
            .filter(|block| block.core)
            .flat_map(|block| block.lines_at(graph, vertex))
            .collect();
        if !core_lines.is_empty() {
            return masslessness(core_lines);
        }
        let towards = self
            .containing(vertex)
            .find(|block| !std::ptr::eq(*block, from) && block.legs_at(vertex).is_empty())
            .expect("a path from a vacuum part to the core");
        masslessness(towards.lines_at(graph, vertex))
    }

    /// Every vacuum cycle with the line it sits on: a snail when it touches a core vertex or
    /// another cycle, a tadpole when bridges lead to that line. A vacuum graph with a loop is
    /// a tadpole attached to a massive line.
    fn vacuum_attachments(&self, graph: &PlainGraph) -> Vec<VacuumAttachment> {
        if self.total_legs == 0 {
            return self
                .blocks
                .iter()
                .any(|block| block.cyclic)
                .then_some(VacuumAttachment {
                    kind: Attachment::Tadpole,
                    massless: Some(false),
                })
                .into_iter()
                .collect();
        }
        let on_line = |vertex: usize, from: &Block| {
            self.containing(vertex)
                .any(|block| !std::ptr::eq(block, from) && (block.cyclic || block.core))
        };
        let mut attachments = vec![];
        for cycle in self
            .blocks
            .iter()
            .filter(|block| block.cyclic && !block.core)
        {
            let attachment = *cycle
                .vertices
                .iter()
                .find(|vertex| !cycle.legs_at(**vertex).is_empty())
                .expect("the legs hang at one vertex of a vacuum cycle");
            if on_line(attachment, cycle) {
                attachments.push(VacuumAttachment {
                    kind: Attachment::Snail,
                    massless: self.line_towards_core(graph, attachment, cycle),
                });
                continue;
            }
            // follow the bridges towards the core
            let (mut vertex, mut from) = (attachment, cycle);
            loop {
                let bridge = self
                    .containing(vertex)
                    .find(|block| !std::ptr::eq(*block, from) && block.legs_at(vertex).is_empty())
                    .expect("a path from a vacuum part to the core");
                let next = bridge.vertices.iter().copied().find(|&end| end != vertex);
                let next = next.expect("a bridge has two ends");
                if bridge.core || on_line(next, bridge) {
                    attachments.push(VacuumAttachment {
                        kind: Attachment::Tadpole,
                        massless: self.line_towards_core(graph, next, bridge),
                    });
                    break;
                }
                (vertex, from) = (next, bridge);
            }
        }
        attachments
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn all_vetoes(snails: bool) -> SpecialTopologyVetoes {
        SpecialTopologyVetoes {
            self_energy: Some(Veto {
                massive: true,
                massless: true,
                ..Default::default()
            }),
            tadpoles: Some(Veto {
                massive: true,
                massless: true,
                ..Default::default()
            }),
            snails: snails.then_some(Veto {
                massive: true,
                massless: true,
                ..Default::default()
            }),
            factorized_loop_range: Range::Unrestricted,
        }
    }

    /// A propagator between two legs with a bubble inserted on it.
    fn bubble_on_line() -> PlainGraph {
        let mut graph = PlainGraph::default();
        let (a, b) = (graph.add_vertex(), graph.add_vertex());
        graph.add_leg(a, true);
        graph.add_leg(b, true);
        graph.add_edge((a, b), Some(true));
        graph.add_edge((a, b), Some(true));
        graph
    }

    #[test]
    fn one_particle_irreducible_two_point_graphs_are_not_self_energy_insertions() {
        assert!(!all_vetoes(true).vetoes(&bubble_on_line()));
    }

    #[test]
    fn reducible_two_point_graphs_are_self_energy_insertions() {
        // leg - bubble - bridge - bubble - leg
        let mut graph = PlainGraph::default();
        let [a, b, c, d] = [(); 4].map(|_| graph.add_vertex());
        graph.add_leg(a, true);
        graph.add_edge((a, b), Some(true));
        graph.add_edge((a, b), Some(true));
        graph.add_edge((b, c), Some(true));
        graph.add_edge((c, d), Some(true));
        graph.add_edge((c, d), Some(true));
        graph.add_leg(d, false);
        assert!(all_vetoes(false).vetoes(&graph));
        let massless_only = SpecialTopologyVetoes {
            self_energy: Some(Veto {
                massive: false,
                massless: true,
                ..Default::default()
            }),
            ..Default::default()
        };
        // each leg sees the other part as the insertion: both legs must be vetoed
        assert!(!massless_only.vetoes(&graph));
    }

    #[test]
    fn self_energies_on_external_legs_are_vetoed() {
        // leg 0 - a = (bubble) = b - c with legs 1 and 2 at c
        let mut graph = PlainGraph::default();
        let (a, b, c) = (graph.add_vertex(), graph.add_vertex(), graph.add_vertex());
        graph.add_leg(a, true);
        graph.add_edge((a, b), Some(true));
        graph.add_edge((a, b), Some(true));
        graph.add_edge((b, c), Some(true));
        graph.add_leg(c, true);
        graph.add_leg(c, true);
        assert!(all_vetoes(false).vetoes(&graph));
        let massive_only = SpecialTopologyVetoes {
            self_energy: Some(Veto {
                massive: true,
                massless: false,
                ..Default::default()
            }),
            ..Default::default()
        };
        assert!(!massive_only.vetoes(&graph));
    }

    #[test]
    fn unknown_attachment_masses_veto_only_when_both_masses_are_vetoed() {
        // legs at a and b, a tadpole hanging from a through c
        let mut graph = PlainGraph::default();
        let (a, b, c) = (graph.add_vertex(), graph.add_vertex(), graph.add_vertex());
        graph.add_leg(a, true);
        graph.add_edge((a, b), None);
        graph.add_leg(b, true);
        graph.add_edge((a, c), None);
        graph.add_edge((c, c), None);
        assert!(all_vetoes(false).vetoes(&graph));
        let massive_only = SpecialTopologyVetoes {
            tadpoles: Some(Veto {
                massive: true,
                massless: false,
                ..Default::default()
            }),
            ..Default::default()
        };
        assert!(!massive_only.vetoes(&graph));
    }

    #[test]
    fn snails_are_classified_by_the_line_they_sit_on() {
        // a massless self-loop on a vertex of a massive line
        let mut graph = PlainGraph::default();
        let (a, b) = (graph.add_vertex(), graph.add_vertex());
        graph.add_leg(a, false);
        graph.add_leg(a, false);
        graph.add_edge((a, b), Some(false));
        graph.add_edge((b, b), Some(true));
        graph.add_leg(b, false);
        let snails = |veto_attached_to_massive, veto_attached_to_massless| SpecialTopologyVetoes {
            snails: Some(Veto {
                massive: veto_attached_to_massive,
                massless: veto_attached_to_massless,
                ..Default::default()
            }),
            ..Default::default()
        };
        assert!(!snails(false, true).vetoes(&graph));
        assert!(snails(true, false).vetoes(&graph));
    }

    #[test]
    fn one_particle_irreducible_cubic_graphs_have_no_special_topologies() {
        // the two-loop vertex ladder: legs at a, d and e
        let mut graph = PlainGraph::default();
        let [a, b, c, d, e] = [(); 5].map(|_| graph.add_vertex());
        graph.add_leg(a, true);
        graph.add_leg(d, true);
        graph.add_leg(e, true);
        for ends in [(a, b), (a, c), (b, c), (b, d), (c, e), (d, e)] {
            graph.add_edge(ends, Some(true));
        }
        let mut vetoes = all_vetoes(true);
        vetoes.factorized_loop_range = Range::exact(0);
        assert!(!vetoes.vetoes(&graph));
    }

    #[test]
    fn factorized_loops_count_the_blocks_carrying_at_most_one_momentum() {
        let count = |graph: &PlainGraph| Blocks::new(graph).factorized_loops();
        // theta: one block
        let mut theta = PlainGraph::default();
        let (a, b) = (theta.add_vertex(), theta.add_vertex());
        for _ in 0..3 {
            theta.add_edge((a, b), Some(true));
        }
        assert_eq!(count(&theta), 1);
        // glasses: two bubbles joined by a bridge
        let mut glasses = PlainGraph::default();
        let [a, b, c, d] = [(); 4].map(|_| glasses.add_vertex());
        glasses.add_edge((a, b), Some(true));
        glasses.add_edge((a, b), Some(true));
        glasses.add_edge((b, c), Some(true));
        glasses.add_edge((c, d), Some(true));
        glasses.add_edge((c, d), Some(true));
        assert_eq!(count(&glasses), 2);
        // a two-point bubble carries one momentum
        assert_eq!(count(&bubble_on_line()), 1);
    }

    #[test]
    fn a_tadpole_on_a_one_particle_irreducible_two_point_core_is_not_a_self_energy() {
        // the triangle a, b, c with legs at a and b, and a bridge from c to a self-loop at d
        let mut graph = PlainGraph::default();
        let [a, d, b, c] = [(); 4].map(|_| graph.add_vertex());
        graph.add_leg(a, true);
        graph.add_leg(b, true);
        for ends in [(a, b), (a, c), (b, c), (c, d), (d, d)] {
            graph.add_edge(ends, Some(true));
        }
        let self_energy = SpecialTopologyVetoes {
            self_energy: all_vetoes(false).self_energy,
            ..Default::default()
        };
        assert!(!self_energy.vetoes(&graph));
        let tadpoles = SpecialTopologyVetoes {
            tadpoles: all_vetoes(false).tadpoles,
            ..Default::default()
        };
        assert!(tadpoles.vetoes(&graph));
    }

    #[test]
    fn a_bubble_inside_a_one_particle_irreducible_vertex_is_not_a_factorized_loop() {
        // legs at b, c and a; the bubble d = e lies on the internal line c - d - e - b
        let mut graph = PlainGraph::default();
        let [a, b, c, d, e] = [(); 5].map(|_| graph.add_vertex());
        graph.add_leg(b, true);
        graph.add_leg(c, true);
        graph.add_leg(a, true);
        for ends in [(a, b), (a, c), (b, e), (c, d), (d, e), (d, e)] {
            graph.add_edge(ends, Some(true));
        }
        assert_eq!(Blocks::new(&graph).factorized_loops(), 0);
        let mut vetoes = all_vetoes(true);
        vetoes.factorized_loop_range = Range::exact(0);
        assert!(!vetoes.vetoes(&graph));
    }

    #[test]
    fn a_one_particle_irreducible_one_point_graph_has_no_special_topologies() {
        // the leg vertex a joined to both ends of the bubble b = c
        let mut graph = PlainGraph::default();
        let [a, b, c] = [(); 3].map(|_| graph.add_vertex());
        graph.add_leg(a, true);
        for ends in [(a, b), (a, c), (b, c), (b, c)] {
            graph.add_edge(ends, Some(true));
        }
        assert!(!all_vetoes(true).vetoes(&graph));
        // its single block carries the one external momentum
        assert_eq!(Blocks::new(&graph).factorized_loops(), 1);
    }
}
