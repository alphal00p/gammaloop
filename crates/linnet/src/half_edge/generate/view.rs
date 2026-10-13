//! What the filter of [`super::GraphGenerator`] sees: read-only views of the search states of
//! both stages and of the complete graphs.

use std::{cell::OnceCell, iter};

use itertools::Either;

use super::{
    decorate::{Decoration, Edge},
    model::Direction,
    topology::Topology,
    Compiled, Filter, Stub,
};
use crate::half_edge::{
    builder::HedgeGraphBuilder,
    involution::{Flow, Orientation},
    HedgeGraph, NodeIndex,
};

/// A search state as the stages pass it to the filter.
#[derive(Clone, Copy)]
pub(crate) enum State<'a> {
    Topology(&'a Topology),
    Decoration(Decoration<'a>),
}

/// The filter as the stages call it on their partial states.
pub(crate) trait Prune: Sync {
    fn keeps(&self, state: State<'_>) -> bool;

    /// Search without pruning partial states, marking the rejected ones for the heredity check.
    fn checks_heredity(&self) -> bool;
}

/// The user filter with what its views need.
pub(crate) struct Pruner<'a, F, L> {
    pub(super) compiled: &'a Compiled<F, L>,
    pub(super) filter: &'a Filter<F, L>,
    pub(super) loops: usize,
    pub(super) check_heredity: bool,
}

impl<F: Send + Sync, L: Send + Sync> Prune for Pruner<'_, F, L> {
    fn keeps(&self, state: State<'_>) -> bool {
        (self.filter)(&GraphView::new(self.compiled, self.loops, state))
    }

    fn checks_heredity(&self) -> bool {
        self.check_heredity
    }
}

/// How far a [`GraphView`] is built.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Stage {
    /// An untyped, connected topology, partial or complete: vertices with their degrees, attached
    /// legs, self-loops and edge multiplicities, and open stubs that later vertices or legs will
    /// take. Legs not attached yet go to later vertices or to open stubs.
    ///
    /// Three or more legs with distinct labels and a common stub are placed as anonymous stubs
    /// while vertices and edges are added, and labelled one at a time once every vertex and edge is
    /// placed: until its label is assigned, such a leg is an open stub of its vertex and an
    /// unplaced leg.
    Topology,
    /// A complete topology with fields and orientations on some of its edges.
    Fields,
    /// A complete graph: every edge carries a field and every vertex a rule.
    Complete,
}

/// A state of the search or a complete graph, as the filter of
/// [`super::GraphGenerator::filter`] sees it.
///
/// Vertices are numbered `0..n_vertices()` and keep their numbers in every extension of a state:
/// stage 1 adds vertices with higher numbers, stage 2 adds fields and rules.
pub struct GraphView<'a, F, L> {
    compiled: &'a Compiled<F, L>,
    loops: usize,
    state: State<'a>,
    complete: OnceCell<HedgeGraph<F, usize, Option<L>>>,
}

/// An edge of a [`GraphView`].
#[derive(Debug, PartialEq, Eq)]
pub struct EdgeView<'a, F> {
    /// Its vertices, equal for a self-loop.
    pub vertices: [usize; 2],
    /// Its field and direction at `vertices[0]` once assigned: an outgoing stub makes the edge
    /// point from `vertices[0]` to `vertices[1]` (directed self-loops have an outgoing stub).
    pub stub: Option<Stub<&'a F>>,
}

/// A bridge of a [`GraphView`] together with the side of it that is final: that side has no open
/// stubs, so no later vertex or edge attaches to it and the bridge stays a bridge.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Bridge {
    /// The ends of the bridge; `ends[1]` is in `side`.
    pub ends: [usize; 2],
    /// The vertices on the final side.
    pub side: VertexSet,
}

/// A set of vertices of a [`GraphView`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct VertexSet(pub(super) u64);

impl VertexSet {
    pub fn contains(self, v: usize) -> bool {
        v < 64 && self.0 >> v & 1 == 1
    }

    pub fn len(self) -> usize {
        self.0.count_ones() as usize
    }

    pub fn is_empty(self) -> bool {
        self.0 == 0
    }

    /// The vertices in increasing order.
    pub fn iter(self) -> impl Iterator<Item = usize> {
        let mut bits = self.0;
        iter::from_fn(move || {
            (bits != 0).then(|| {
                let v = bits.trailing_zeros() as usize;
                bits &= bits - 1;
                v
            })
        })
    }
}

impl FromIterator<usize> for VertexSet {
    /// The set of the given vertices; panics on a vertex of 64 or more, which no view has.
    fn from_iter<I: IntoIterator<Item = usize>>(vertices: I) -> Self {
        VertexSet(vertices.into_iter().fold(0, |bits, v| {
            assert!(v < 64, "vertex {v} is out of range");
            bits | 1 << v
        }))
    }
}

impl<'a, F, L> GraphView<'a, F, L> {
    pub(super) fn new(compiled: &'a Compiled<F, L>, loops: usize, state: State<'a>) -> Self {
        GraphView {
            compiled,
            loops,
            state,
            complete: OnceCell::new(),
        }
    }

    fn topology(&self) -> &'a Topology {
        match self.state {
            State::Topology(t) => t,
            State::Decoration(d) => d.topology,
        }
    }

    pub fn stage(&self) -> Stage {
        match self.state {
            State::Topology(_) => Stage::Topology,
            State::Decoration(d) if d.rules.is_empty() => Stage::Fields,
            State::Decoration(_) => Stage::Complete,
        }
    }

    /// The view is a complete graph.
    pub fn is_complete(&self) -> bool {
        self.stage() == Stage::Complete
    }

    /// Every vertex, edge and leg is placed: no open stubs and no unattached legs. Always true
    /// from [`Stage::Fields`] on.
    pub fn is_topology_complete(&self) -> bool {
        self.topology().is_complete()
    }

    /// The loop order being generated.
    pub fn loop_order(&self) -> usize {
        self.loops
    }

    pub fn n_vertices(&self) -> usize {
        self.topology().n
    }

    pub fn degree(&self, v: usize) -> usize {
        self.topology().vertices[v].degree as usize
    }

    /// Stubs of `v` not taken by an edge or a leg yet; zero from [`Stage::Fields`] on.
    pub fn open_stubs(&self, v: usize) -> usize {
        self.topology().free_stubs(&self.compiled.model, v) as usize
    }

    pub fn self_loops(&self, v: usize) -> usize {
        self.topology().vertices[v].loops as usize
    }

    /// Number of edges between the distinct vertices `u` and `v`.
    pub fn multiplicity(&self, u: usize, v: usize) -> usize {
        self.topology().multiplicity(u, v) as usize
    }

    /// The legs attached to `v`, with their stubs at `v`.
    pub fn legs(&self, v: usize) -> impl Iterator<Item = (&'a L, Stub<&'a F>)> + 'a {
        let (compiled, packed) = (self.compiled, self.topology().vertices[v].labelled);
        compiled
            .tables
            .legs
            .iter()
            .enumerate()
            .flat_map(move |(c, leg)| {
                let stub = Stub {
                    field: &compiled.tables.fields[leg.field as usize],
                    flow: leg.flow,
                };
                let count = compiled.model.legs.count(packed, c) as usize;
                iter::repeat_n((&compiled.tables.labels[leg.label as usize], stub), count)
            })
    }

    /// Legs not attached to a vertex yet; zero from [`Stage::Fields`] on.
    pub fn unplaced_legs(&self) -> usize {
        let t = self.topology();
        let placed: usize = t.vertices[..t.n]
            .iter()
            .map(|v| v.n_legs as usize - v.unlabelled(&self.compiled.model) as usize)
            .sum();
        self.compiled.model.legs.total as usize - placed
    }

    /// Independent cycles of the edges placed so far.
    pub fn loops(&self) -> usize {
        let t = self.topology();
        let paired: usize = t.vertices[..t.n]
            .iter()
            .map(|v| (v.degree - v.n_legs - v.open) as usize)
            .sum();
        paired / 2 + 1 - t.n
    }

    /// Whether `side` is non-empty and connected by the edges placed so far.
    pub fn is_connected(&self, side: VertexSet) -> bool {
        let Some(first) = side.iter().next() else {
            return false;
        };
        let (mut reached, mut stack) = (VertexSet(1 << first), vec![first]);
        while let Some(u) = stack.pop() {
            for v in side.iter() {
                if !reached.contains(v) && self.multiplicity(u, v) > 0 {
                    reached.0 |= 1 << v;
                    stack.push(v);
                }
            }
        }
        reached == side
    }

    /// Independent cycles of the subgraph on the connected set `side` (see
    /// [`is_connected`](Self::is_connected)), such as the side of a [`Bridge`].
    pub fn loops_in(&self, side: VertexSet) -> usize {
        let edges: usize = side
            .iter()
            .map(|u| {
                self.self_loops(u)
                    + side
                        .iter()
                        .filter(|&v| v > u)
                        .map(|v| self.multiplicity(u, v))
                        .sum::<usize>()
            })
            .sum();
        (edges + 1).saturating_sub(side.len())
    }

    /// Every edge between vertices. From [`Stage::Fields`] on, edges are listed in the order
    /// stage 2 assigns their fields: the edges with a field come first.
    pub fn edges(&self) -> impl Iterator<Item = EdgeView<'a, F>> + 'a {
        match self.state {
            State::Topology(t) => Either::Left((0..t.n).flat_map(move |u| {
                (u..t.n).flat_map(move |v| {
                    let m = if u == v {
                        t.vertices[u].loops
                    } else {
                        t.multiplicity(u, v)
                    };
                    (0..m).map(move |_| EdgeView {
                        vertices: [u, v],
                        stub: None,
                    })
                })
            })),
            State::Decoration(d) => Either::Right(self.compiled.edges(d)),
        }
    }

    /// The signatures that fit `v`: once all its stubs have fields (from [`Stage::Fields`] on),
    /// those with exactly its stubs, and in a complete graph the rule of `v`.
    pub fn rules(&self, v: usize) -> Option<&'a [u32]> {
        let State::Decoration(d) = self.state else {
            return None;
        };
        if !d.rules.is_empty() {
            return Some(&d.rules[v..=v]);
        }
        // only the full states of the vertex automaton have rules
        let rules = self.compiled.model.rules(d.states[v]);
        (!rules.is_empty()).then_some(rules)
    }

    /// The bridges with a final side; a bridge of a complete topology is listed once per side.
    pub fn bridges(&self) -> Vec<Bridge> {
        let mut bridges = vec![];
        self.topology()
            .closed_sides(&self.compiled.model, |ends, side| {
                bridges.push(Bridge {
                    ends,
                    side: VertexSet(side),
                })
            });
        bridges
    }

    /// The state as a graph: the edges placed so far, with their field once assigned (undirected
    /// before), and the legs; nodes carry their rule once it is the only one that fits. Open
    /// stubs are left out.
    pub fn graph(&self) -> HedgeGraph<Option<F>, Option<usize>, Option<L>>
    where
        F: Clone,
        L: Clone,
    {
        Self::build(
            self.n_vertices(),
            self.compiled.model.legs.total as usize,
            (0..self.n_vertices())
                .flat_map(|v| self.legs(v).map(move |(label, stub)| (v, label, stub))),
            self.edges(),
            |v| match self.rules(v) {
                Some(&[rule]) => Some(rule as usize),
                _ => None,
            },
            |stub| stub.map(|s| s.field.clone()),
        )
    }

    /// The complete graph, as the generator emits it; `None` before [`Stage::Complete`].
    pub fn complete_graph(&self) -> Option<&HedgeGraph<F, usize, Option<L>>>
    where
        F: Clone,
        L: Clone,
    {
        let State::Decoration(d) = self.state else {
            return None;
        };
        (!d.rules.is_empty()).then(|| {
            self.complete.get_or_init(|| {
                Self::build(
                    self.n_vertices(),
                    self.compiled.model.legs.total as usize,
                    (0..self.n_vertices())
                        .flat_map(|v| self.legs(v).map(move |(label, stub)| (v, label, stub))),
                    self.compiled.edges(d),
                    |v| d.rules[v] as usize,
                    |stub| {
                        stub.expect("a complete graph has every field")
                            .field
                            .clone()
                    },
                )
            })
        })
    }

    /// Transfer only a complete graph already requested by the final filter.
    pub(crate) fn take_complete_graph(&mut self) -> Option<HedgeGraph<F, usize, Option<L>>> {
        self.complete.take()
    }

    /// The graph of the placed `edges` and the legs.
    pub(super) fn build<E, N>(
        n_vertices: usize,
        n_legs: usize,
        legs: impl Iterator<Item = (usize, &'a L, Stub<&'a F>)>,
        edges: impl Iterator<Item = EdgeView<'a, F>>,
        node: impl Fn(usize) -> N,
        edge: impl Fn(Option<Stub<&'a F>>) -> E,
    ) -> HedgeGraph<E, N, Option<L>>
    where
        F: 'a,
        L: Clone + 'a,
    {
        let n_hedges = n_legs + 2 * edges.size_hint().0;
        let mut builder = HedgeGraphBuilder::with_capacity(n_vertices, n_hedges);
        for v in 0..n_vertices {
            builder.add_node(node(v));
        }
        for (vertex, label, stub) in legs {
            let (orientation, underlying) = match stub.flow {
                None => (Orientation::Undirected, Flow::Source),
                Some(flow) => (Orientation::Default, flow),
            };
            builder.add_external_edge(
                NodeIndex(vertex).add_data(Some(label.clone())),
                edge(Some(stub)),
                orientation,
                underlying,
            );
        }
        for e in edges {
            let [u, v] = e.vertices.map(NodeIndex);
            match e.stub.as_ref().and_then(|s| s.flow) {
                None => builder.add_edge(u, v, edge(e.stub), Orientation::Undirected),
                Some(Flow::Source) => builder.add_edge(u, v, edge(e.stub), Orientation::Default),
                Some(Flow::Sink) => builder.add_edge(v, u, edge(e.stub), Orientation::Default),
            }
        }
        builder.build()
    }
}

impl<F, L> Compiled<F, L> {
    /// The edges of a decoration, with the stubs of the labelled ones.
    fn edges<'a>(&'a self, d: Decoration<'a>) -> impl Iterator<Item = EdgeView<'a, F>> + 'a {
        d.edges.iter().enumerate().map(move |(i, e)| EdgeView {
            vertices: [e.u as usize, e.v as usize],
            stub: d.labels.get(i).map(|&l| self.stub(e, l)),
        })
    }

    /// The field rank and direction of a decorated edge at its first vertex.
    pub(super) fn edge_stub(&self, e: &Edge, l: u8) -> Stub<u8> {
        let model = &self.model;
        let (field, flow) = if e.u == e.v {
            let label = model.loop_labels[l as usize];
            (label.field, label.directed.then_some(Flow::Source))
        } else {
            let label = model.labels[l as usize];
            let flow = match label.direction {
                Direction::Undirected => None,
                Direction::Forward => Some(Flow::Source),
                Direction::Backward => Some(Flow::Sink),
            };
            (label.field, flow)
        };
        Stub { field, flow }
    }

    fn stub(&self, e: &Edge, label: u8) -> Stub<&F> {
        let stub = self.edge_stub(e, label);
        Stub {
            field: &self.tables.fields[stub.field as usize],
            flow: stub.flow,
        }
    }
}
