//! Read-only typed views of skeletons, dressing states and complete assignments.

use std::{cell::OnceCell, iter, sync::OnceLock};

use itertools::Either;

use super::{
    decorate::{Decoration, Edge},
    model::{Direction, Model},
    output::SkeletonData,
    topology::Topology,
    CanonicalForm, Compiled, DiagramEdge, Hooks, OutputTables, Stub,
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
    pub(super) hooks: &'a Hooks<F, L>,
    pub(super) loops: usize,
    pub(super) check_heredity: bool,
}

impl<F: Send + Sync, L: Send + Sync> Prune for Pruner<'_, F, L> {
    fn keeps(&self, state: State<'_>) -> bool {
        let view = View::new(self.compiled, self.loops, state);
        match state {
            State::Topology(_) => self
                .hooks
                .skeletons
                .as_ref()
                .is_none_or(|hook| hook(&SkeletonView(view))),
            State::Decoration(_) => self
                .hooks
                .dressings
                .as_ref()
                .is_none_or(|hook| hook(&DressingView(view))),
        }
    }

    fn checks_heredity(&self) -> bool {
        self.check_heredity
    }
}

/// How far a [`View`] is built.
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

/// Shared private representation of the typed borrowed views.
///
/// Vertices are numbered `0..n_vertices()` and keep their numbers in every extension of a state:
/// stage 1 adds vertices with higher numbers, stage 2 adds fields and rules.
pub(super) struct View<'a, F, L> {
    compiled: Option<&'a Compiled<F, L>>,
    tables: &'a OutputTables<F, L>,
    skeleton: Option<&'a SkeletonData<F, L>>,
    loops: usize,
    state: State<'a>,
    complete: OnceCell<HedgeGraph<F, usize, Option<L>>>,
    /// Complete batch assignments share this cache across scoped views and hook stages.
    complete_cache: Option<&'a OnceLock<HedgeGraph<F, usize, Option<L>>>>,
    skeleton_id: OnceCell<CanonicalForm<Option<usize>, Option<F>, Option<L>>>,
}

// Typed boundaries share the compact representation and accessor semantics. The private View
// cannot be supplied as a public callback, so a hook cannot accidentally run at another stage.
macro_rules! typed_views {
    ($($(#[$doc:meta])* $name:ident),+ $(,)?) => {$ (
        $(#[$doc])*
        pub struct $name<'a, F, L>(pub(super) View<'a, F, L>);

        impl<'a, F, L> $name<'a, F, L> {
            pub fn stage(&self) -> Stage { self.0.stage() }
            pub fn is_complete(&self) -> bool { self.0.is_complete() }
            pub fn is_topology_complete(&self) -> bool { self.0.is_topology_complete() }
            pub fn loop_order(&self) -> usize { self.0.loop_order() }
            pub fn n_vertices(&self) -> usize { self.0.n_vertices() }
            pub fn degree(&self, vertex: usize) -> usize { self.0.degree(vertex) }
            pub fn open_stubs(&self, vertex: usize) -> usize { self.0.open_stubs(vertex) }
            pub fn self_loops(&self, vertex: usize) -> usize { self.0.self_loops(vertex) }
            pub fn multiplicity(&self, source: usize, target: usize) -> usize {
                self.0.multiplicity(source, target)
            }
            pub fn legs(&self, vertex: usize) -> impl Iterator<Item = (&'a L, Stub<&'a F>)> + 'a {
                self.0.legs(vertex)
            }
            pub fn unplaced_legs(&self) -> usize { self.0.unplaced_legs() }
            pub fn loops(&self) -> usize { self.0.loops() }
            pub fn is_connected(&self, side: VertexSet) -> bool { self.0.is_connected(side) }
            pub fn loops_in(&self, side: VertexSet) -> usize { self.0.loops_in(side) }
            pub fn edges(&self) -> impl Iterator<Item = EdgeView<'a, F>> + 'a { self.0.edges() }
            /// Assigned class or particle of an edge. Class members retain the original field order.
            pub fn field(&self, edge: usize) -> Option<Dressed<'a, F>> { self.0.field(edge) }
            /// Known direction at the first endpoint, including unresolved directed classes.
            /// None means undirected or not yet assigned; field() distinguishes those cases.
            pub fn flow(&self, edge: usize) -> Option<Flow> { self.0.flow(edge) }
            pub fn rules(&self, vertex: usize) -> Option<&'a [u32]> { self.0.rules(vertex) }
            /// Rules compatible with the assigned stubs at this vertex, in declaration order.
            pub fn domains(&self, vertex: usize) -> Vec<u32> { self.0.domains(vertex) }
            /// Bounds on declared theory weights for a complete skeleton. `None` before its
            /// vertex set is final, or when a vertex has no possible rule.
            pub fn order_bounds(&self) -> Option<Vec<OrderBound<'a>>> { self.0.order_bounds() }
            /// Stable topology identity once every edge and external label is placed.
            pub fn skeleton_id(&self) -> Option<&CanonicalForm<Option<usize>, Option<F>, Option<L>>>
            where F: Clone + Ord, L: Clone + Ord { self.0.skeleton_id() }
            pub fn bridges(&self) -> Vec<Bridge> { self.0.bridges() }
            pub fn graph(&self) -> HedgeGraph<Option<F>, Option<usize>, Option<L>>
            where F: Clone, L: Clone {
                self.0.graph()
            }
            pub fn complete_graph(&self) -> Option<&HedgeGraph<F, usize, Option<L>>>
            where F: Clone, L: Clone {
                self.0.complete_graph()
            }
        }
    )+};
}

typed_views! {
    /// A partial or complete untyped topology. Vertex numbers are retained by its refinements.
    /// Open stubs and unplaced legs describe the unfinished part; edge fields are unassigned.
    SkeletonView,
    /// A complete topology with a partial or complete assignment of fields and vertex rules.
    /// Pruning must remain hereditary as fields and rules become known.
    DressingView,
    /// A fully assigned diagram, after built-in restrictions and dressing pruning.
    /// Final selection may use any predicate; graph construction remains an explicit demand.
    DiagramView,
}

/// A field at the current dressing stage, borrowing the original immutable field table.
/// A class is not a representative particle: predicates must hold for its refinements.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Dressed<'a, F> {
    Class(&'a [F]),
    Particle(&'a F),
}

/// An inclusive range of the possible sum of one declared theory weight.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct OrderBound<'a> {
    pub name: &'a str,
    pub range: [i64; 2],
}

/// An edge of a [`View`].
#[derive(Debug, PartialEq, Eq)]
pub struct EdgeView<'a, F> {
    /// Its vertices, equal for a self-loop.
    pub vertices: [usize; 2],
    /// Its field and direction at `vertices[0]` once assigned: an outgoing stub makes the edge
    /// point from `vertices[0]` to `vertices[1]` (directed self-loops have an outgoing stub).
    pub stub: Option<Stub<&'a F>>,
    /// Direction is known even while a class has no concrete particle stub. None means
    /// undirected or still unassigned, as distinguished by the view's field accessor.
    pub flow: Option<Flow>,
}

/// A bridge of a [`View`] together with the side of it that is final: that side has no open
/// stubs, so no later vertex or edge attaches to it and the bridge stays a bridge.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Bridge {
    /// The ends of the bridge; `ends[1]` is in `side`.
    pub ends: [usize; 2],
    /// The vertices on the final side.
    pub side: VertexSet,
}

/// A set of vertices of a [`View`].
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

impl<'a, F, L> View<'a, F, L> {
    pub(super) fn new(compiled: &'a Compiled<F, L>, loops: usize, state: State<'a>) -> Self {
        View {
            compiled: Some(compiled),
            tables: &compiled.tables,
            skeleton: None,
            loops,
            state,
            complete: OnceCell::new(),
            complete_cache: None,
            skeleton_id: OnceCell::new(),
        }
    }

    pub(super) fn stored(skeleton: &'a SkeletonData<F, L>) -> Self {
        Self {
            compiled: None,
            tables: &skeleton.tables,
            skeleton: Some(skeleton),
            loops: skeleton.loops,
            state: State::Topology(&skeleton.topology),
            complete: OnceCell::new(),
            complete_cache: None,
            skeleton_id: OnceCell::new(),
        }
    }

    pub(super) fn with_complete_cache(
        mut self,
        cache: Option<&'a OnceLock<HedgeGraph<F, usize, Option<L>>>>,
    ) -> Self {
        self.complete_cache = cache;
        self
    }

    /// Restore a demanded batch graph to the same owner used by scalar output freezing.
    pub(super) fn with_complete_graph(
        mut self,
        graph: Option<HedgeGraph<F, usize, Option<L>>>,
    ) -> Self {
        if let Some(graph) = graph {
            self.complete = OnceCell::from(graph);
        }
        self
    }

    fn model(&self) -> &'a Model {
        let model = &self
            .compiled
            .expect("a search view borrows its compiled plan")
            .model;
        match self.state {
            State::Decoration(d) if d.particles => model.particles(),
            _ => model,
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
        {
            let vertex = &self.topology().vertices[v];
            vertex.open as usize + vertex.n_legs as usize
                - self.tables.packing.size(vertex.labelled) as usize
        }
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
        let (tables, packed) = (self.tables, self.topology().vertices[v].labelled);
        tables
            .legs
            .iter()
            .enumerate()
            .flat_map(move |(class, leg)| {
                let stub = Stub {
                    field: &tables.fields[leg.field as usize],
                    flow: leg.flow,
                };
                let count = tables.packing.count(packed, class) as usize;
                iter::repeat_n((&tables.labels[leg.label as usize], stub), count)
            })
    }

    /// Legs not attached to a vertex yet; zero from [`Stage::Fields`] on.
    pub fn unplaced_legs(&self) -> usize {
        let t = self.topology();
        let placed: usize = t.vertices[..t.n]
            .iter()
            .map(|v| self.tables.packing.size(v.labelled) as usize)
            .sum();
        self.tables.packing.total as usize - placed
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
    /// stage 2 assigns their classes. Particle refinement assigns a whole directed line at once,
    /// so resolved particle fields need not form a prefix.
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
                        flow: None,
                    })
                })
            })),
            State::Decoration(d) => Either::Right(self.compiled.expect("decoration view").edges(d)),
        }
    }

    fn field(&self, edge: usize) -> Option<Dressed<'a, F>> {
        match self.state {
            State::Topology(_) => None,
            State::Decoration(d) => self.compiled.expect("decoration view").dressed(d, edge),
        }
    }

    fn flow(&self, edge: usize) -> Option<Flow> {
        match self.state {
            State::Topology(_) => None,
            State::Decoration(d) => self
                .compiled
                .expect("decoration view")
                .decoration_edge(d, edge)
                .and_then(|e| e.stub.flow),
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
        let rules = self.model().rules(d.states[v]);
        (!rules.is_empty()).then_some(rules)
    }

    fn domains(&self, vertex: usize) -> Vec<u32> {
        let t = self.topology();
        if let Some(compiled) = self.compiled {
            let state = match self.state {
                State::Decoration(d) if !d.rules.is_empty() => return vec![d.rules[vertex]],
                State::Decoration(d) => Some(d.states[vertex]),
                State::Topology(_) => compiled
                    .model
                    .leg_state(t.vertices[vertex].degree, t.vertices[vertex].legs),
            };
            return state.map(|s| self.model().domain(s)).unwrap_or_default();
        }
        // Stored skeletons keep rule declarations but no automaton. Match the attached stubs
        // directly, preserving duplicate rules and their original insertion indices.
        let mut stubs = Vec::new();
        let v = &t.vertices[vertex];
        for (class, leg) in self.tables.legs.iter().enumerate() {
            let direction = Stub {
                field: leg.field,
                flow: leg.flow,
            }
            .direction();
            stubs.extend(iter::repeat_n(
                3 * leg.field + direction,
                self.tables.packing.count(v.labelled, class) as usize,
            ));
        }
        self.tables
            .signatures
            .iter()
            .enumerate()
            .filter_map(|(rule, signature)| {
                if signature.len() != v.degree as usize {
                    return None;
                }
                let mut remaining = signature.clone();
                for &stub in &stubs {
                    let index = remaining.iter().position(|&s| s == stub)?;
                    remaining.swap_remove(index);
                }
                Some(rule as u32)
            })
            .collect()
    }

    fn order_bounds(&self) -> Option<Vec<OrderBound<'a>>> {
        if !self.is_topology_complete() {
            return None;
        }
        let domains: Vec<_> = (0..self.n_vertices()).map(|v| self.domains(v)).collect();
        if domains.iter().any(Vec::is_empty) {
            return None;
        }
        Some(
            self.tables
                .weights
                .iter()
                .map(|(name, weights)| {
                    let mut range = [0i64; 2];
                    for domain in &domains {
                        let values = domain.iter().map(|&rule| i64::from(weights[rule as usize]));
                        range[0] += values.clone().min().expect("nonempty domain");
                        range[1] += values.max().expect("nonempty domain");
                    }
                    OrderBound { name, range }
                })
                .collect(),
        )
    }

    fn skeleton_id(&self) -> Option<&CanonicalForm<Option<usize>, Option<F>, Option<L>>>
    where
        F: Clone + Ord,
        L: Clone + Ord,
    {
        if !self.is_topology_complete() {
            return None;
        }
        Some(match self.skeleton {
            Some(skeleton) => skeleton.id(),
            None => self.skeleton_id.get_or_init(|| {
                let view = Self::new(
                    self.compiled.expect("search view"),
                    self.loops,
                    State::Topology(self.topology()),
                );
                view.graph()
                    .canonize_by(Clone::clone, Clone::clone, Clone::clone)
                    .canonical_form()
                    .clone()
            }),
        })
    }

    /// The bridges with a final side; a bridge of a complete topology is listed once per side.
    pub fn bridges(&self) -> Vec<Bridge> {
        let mut bridges = vec![];
        self.topology()
            .closed_sides(&self.tables.packing, |ends, side| {
                bridges.push(Bridge {
                    ends,
                    side: VertexSet(side),
                })
            });
        bridges
    }

    /// The state as a graph: assigned orientations and resolved particle fields, and the legs.
    /// A class has no concrete particle payload yet. Nodes carry their rule once unique. Open
    /// stubs are left out.
    pub fn graph(&self) -> HedgeGraph<Option<F>, Option<usize>, Option<L>>
    where
        F: Clone,
        L: Clone,
    {
        Self::build(
            self.n_vertices(),
            self.tables.packing.total as usize,
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
            let build = || {
                Self::build(
                    self.n_vertices(),
                    self.tables.packing.total as usize,
                    (0..self.n_vertices())
                        .flat_map(|v| self.legs(v).map(move |(label, stub)| (v, label, stub))),
                    self.compiled.expect("decoration view").edges(d),
                    |v| d.rules[v] as usize,
                    |stub| {
                        stub.expect("a complete graph has every field")
                            .field
                            .clone()
                    },
                )
            };
            match self.complete_cache {
                Some(cache) => cache.get_or_init(build),
                None => self.complete.get_or_init(build),
            }
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
            match e.flow {
                None => builder.add_edge(u, v, edge(e.stub), Orientation::Undirected),
                Some(Flow::Source) => builder.add_edge(u, v, edge(e.stub), Orientation::Default),
                Some(Flow::Sink) => builder.add_edge(v, u, edge(e.stub), Orientation::Default),
            }
        }
        builder.build()
    }
}

impl<F, L> Compiled<F, L> {
    /// The edges retain orientation before an actual class member is selected.
    fn edges<'a>(&'a self, d: Decoration<'a>) -> impl Iterator<Item = EdgeView<'a, F>> + 'a {
        d.edges.iter().enumerate().map(move |(i, e)| {
            let flow = self.decoration_edge(d, i).and_then(|edge| edge.stub.flow);
            let stub = match self.dressed(d, i) {
                Some(Dressed::Particle(field)) => Some(Stub { field, flow }),
                _ => None,
            };
            EdgeView {
                vertices: [e.u as usize, e.v as usize],
                stub,
                flow,
            }
        })
    }

    fn decoration_edge(&self, d: Decoration<'_>, index: usize) -> Option<DiagramEdge> {
        let edge = d.edges.get(index)?;
        let label = d
            .labels
            .get(index)
            .copied()
            .filter(|&label| label != u8::MAX);
        if d.particles {
            if let Some(label) = label {
                return Some(Self::edge_with(self.model.particles(), edge, label));
            }
            Some(Self::edge_with(
                &self.model,
                edge,
                *d.coarse_labels?.get(index)?,
            ))
        } else {
            Some(Self::edge_with(&self.model, edge, label?))
        }
    }

    fn dressed<'a>(&'a self, d: Decoration<'_>, index: usize) -> Option<Dressed<'a, F>> {
        let edge = self.decoration_edge(d, index)?;
        let fine = d.particles && d.labels.get(index).is_some_and(|&label| label != u8::MAX);
        if fine || self.model.refinement.is_none() {
            Some(Dressed::Particle(
                &self.tables.fields[edge.stub.field as usize],
            ))
        } else {
            let members = &self.tables.class_members[edge.stub.field as usize];
            Some(match members.as_ref() {
                [field] => Dressed::Particle(field),
                fields => Dressed::Class(fields),
            })
        }
    }

    /// The decorated edge, with its field rank and direction at the first vertex.
    pub(super) fn edge(&self, e: &Edge, l: u8) -> DiagramEdge {
        Self::edge_with(self.model.particles(), e, l)
    }

    fn edge_with(model: &Model, e: &Edge, l: u8) -> DiagramEdge {
        let (field, flow) = if e.u == e.v {
            let label = model.theory.loop_labels[l as usize];
            (label.field, label.directed.then_some(Flow::Source))
        } else {
            let label = model.theory.labels[l as usize];
            let flow = match label.direction {
                Direction::Undirected => None,
                Direction::Forward => Some(Flow::Source),
                Direction::Backward => Some(Flow::Sink),
            };
            (label.field, flow)
        };
        DiagramEdge {
            vertices: [e.u, e.v],
            stub: Stub { field, flow },
        }
    }
}
