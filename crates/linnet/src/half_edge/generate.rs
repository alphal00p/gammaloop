//! # Feynman-graph generation
//!
//! [`GraphGenerator`] enumerates the connected graphs built from vertex signatures (multisets of
//! typed stubs) and a set of external legs, one per isomorphism class, each with its exact
//! automorphism count and the certificate of the native canonizer.
//!
//! A stub is a field with a direction relative to its vertex: undirected (self-conjugate
//! fields), incoming or outgoing (fermions, ghosts, any charged field). Two stubs pair into an
//! edge when they carry the same field and are both undirected, or one is outgoing and the other
//! incoming. External legs are stubs at the vertex they attach to, with a label; they become
//! dangling half-edges. Legs with equal label and stub are interchangeable, so symmetrized
//! external states are expressed by repeating a label.
//!
//! Output graphs are `HedgeGraph<F, usize, Option<L>>`: edges and legs carry their field, nodes
//! the index of their signature, leg half-edges `Some(label)` and paired half-edges `None`.
//! Directed edges point from the outgoing to the incoming stub; legs keep the direction of their
//! stub. The certificate and automorphism count are those of
//! `graph.canonize_by(|&s| s, |f| f.clone(), |h| h.clone())`, so symmetry factors follow the
//! canonizer's half-edge convention: `m!` for `m` interchangeable parallel edges or legs and `2`
//! for an undirected self-loop. Signatures with identical stubs are distinct vertex rules and give
//! distinct graphs.
//!
//! ```
//! use linnet::half_edge::generate::{GraphGenerator, Stub};
//!
//! // QED photon self-energy at two loops, one-particle irreducible
//! let (photon, electron) = (Stub::undirected('a'), Stub::incoming('e'));
//! let generator = GraphGenerator::new(
//!     [(1, photon.clone()), (2, photon.clone())],
//!     [vec![electron.clone(), electron.flipped(), photon]],
//! )
//! .loops(2..=2)
//! .max_bridges(0);
//! let diagrams: Vec<_> = generator.iter().unwrap().collect();
//! assert_eq!(diagrams.len(), 3);
//! assert!(diagrams.iter().all(|d| d.automorphism_count == Some(1)));
//! ```
//!
//! ## Algorithm
//!
//! Generation runs in two stages per loop order. Untyped topologies (vertex degrees from the
//! signatures, legs attached, loop order, vertex and bridge limits) are built by canonical
//! augmentation one vertex at a time ([`topology`]), except that the pendant vertices completing a
//! topology are added at once. Three or more legs with distinct labels and a common stub are
//! placed as anonymous legs and labelled once the topology is shaped, each label up to the
//! stabilizer of the previous ones, so many distinct legs cost few topologies. Every complete
//! topology is then decorated with fields, orientations and vertex rules, one decoration per orbit
//! of its automorphism group ([`decorate`]). Only accepted diagrams are materialized as
//! [`HedgeGraph`]s and canonized.
//!
//! Both stages are resumable depth-first searches, so [`Diagrams`] streams: each `next` call
//! searches only up to the next diagram and nothing is buffered. Output order is the depth-first
//! order of the search and does not depend on parallelism: with the `rayon` feature,
//! [`GraphGenerator::parallel`] runs independent subtrees on the current thread pool, each at most
//! a fixed number of diagrams ahead of the output, and emits their results in that order. The
//! search starts from a few subtrees per thread; while fewer subtrees than threads can advance,
//! the unsearched rest of the oldest ones is split off into new subtrees. The first step runs only
//! the oldest subtree, up to its first diagram, and further steps buffer one diagram per subtree
//! until the first is emitted, so the first graph does not wait for a step of the whole window.
//!
//! ## Filtering
//!
//! [`GraphGenerator::filter`] prunes both searches. It sees each search state as a [`GraphView`]:
//! partial untyped topologies in stage 1, partial field assignments in stage 2 and the complete
//! graphs, and a rejected state is dropped with everything below it. A hereditary filter (one that
//! rejects every completion of a state it rejects) thus yields exactly the complete graphs it
//! accepts, and cuts the search as early as its conditions are decided:
//!
//! ```
//! use linnet::half_edge::generate::{GraphGenerator, Stub};
//!
//! // gg -> gg at one loop without ghosts: ghost fields are decided edge by edge in stage 2
//! let (g, q, c) = (Stub::undirected('g'), Stub::incoming('q'), Stub::incoming('c'));
//! let signatures = [
//!     vec![q.clone().flipped(), q, g.clone()],
//!     vec![g.clone(); 3],
//!     vec![g.clone(); 4],
//!     vec![c.clone().flipped(), c, g.clone()],
//! ];
//! let generator = GraphGenerator::new((1..=4).map(|l| (l, g.clone())), signatures)
//!     .allow_self_loops(true)
//!     .loops(1..=1);
//! let no_ghosts = generator.filter(|view| {
//!     view.edges().all(|e| e.stub.is_none_or(|s| *s.field != 'c'))
//! });
//! let diagrams: Vec<_> = no_ghosts.iter().unwrap().collect();
//! assert_eq!(diagrams.len(), 111);
//! ```

mod decorate;
mod model;
#[cfg(test)]
mod tests;
mod topology;
mod view;

use std::{
    cell::RefCell,
    collections::VecDeque,
    ops::{Bound, RangeBounds},
    sync::atomic::{AtomicBool, Ordering},
};

use decorate::{Decoration, Decorator, Edge};
use model::{Direction, Model, INCOMING, OUTGOING, UNDIRECTED};
use topology::{Context, Limits, Search, Topology, MAX_VERTICES};
pub use view::{Bridge, EdgeView, GraphView, Stage, VertexSet};
use view::{Prune, Pruner, State};

use super::{
    canon::{Certificate, KeyTables, RawHedge, Role},
    involution::Flow,
    HedgeGraph,
};

/// A typed half-edge slot of a vertex or an external leg.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Stub<F> {
    pub field: F,
    /// `Some(Flow::Source)` for outgoing, `Some(Flow::Sink)` for incoming, `None` for undirected.
    pub flow: Option<Flow>,
}

impl<F> Stub<F> {
    /// A self-conjugate field: pairs with another undirected stub of the same field.
    pub fn undirected(field: F) -> Self {
        Stub { field, flow: None }
    }

    /// A charged field flowing into the vertex: pairs with an outgoing stub of the same field.
    pub fn incoming(field: F) -> Self {
        Stub {
            field,
            flow: Some(Flow::Sink),
        }
    }

    /// A charged field flowing out of the vertex.
    pub fn outgoing(field: F) -> Self {
        Stub {
            field,
            flow: Some(Flow::Source),
        }
    }

    /// The same field with the opposite direction.
    pub fn flipped(self) -> Self {
        Stub {
            field: self.field,
            flow: self.flow.map(|f| -f),
        }
    }

    fn direction(&self) -> u8 {
        match self.flow {
            None => UNDIRECTED,
            Some(Flow::Sink) => INCOMING,
            Some(Flow::Source) => OUTGOING,
        }
    }
}

/// A generated graph.
#[derive(Clone, Debug)]
pub struct Diagram<F, L> {
    /// Fields on edges and legs, signature indices on nodes, `Some(label)` on leg half-edges.
    pub graph: HedgeGraph<F, usize, Option<L>>,
    /// Loop order.
    pub loops: usize,
    /// Order of the half-edge automorphism group (the symmetry factor); `None` if it exceeds
    /// `u128`.
    pub automorphism_count: Option<u128>,
    /// Equal for isomorphic graphs: the certificate of
    /// `graph.canonize_by(|&s| s, |f| f.clone(), |h| h.clone())`.
    pub certificate: Certificate<usize, F, Option<L>>,
}

/// Counts passed to the progress callback.
#[derive(Clone, Debug, Default)]
pub struct Progress {
    /// Loop order being generated.
    pub loops: usize,
    /// Complete topologies decorated so far.
    pub topologies: usize,
    /// Graphs emitted so far.
    pub diagrams: usize,
}

#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
pub enum GenerationError {
    #[error("bound the loop order or the number of vertices (vertices of degree two or less always need a vertex bound)")]
    Unbounded,
    #[error("signatures need between 1 and 15 stubs, found {0}")]
    UnsupportedDegree(usize),
    #[error("at most 84 distinct fields are supported")]
    TooManyFields,
    #[error("too many external legs to pack")]
    TooManyLegs,
    #[error("graphs of this size would need more than 64 vertices")]
    TooManyVertices,
}

type Filter<F, L> = dyn Fn(&GraphView<'_, F, L>) -> bool + Send + Sync;
type ProgressFn = dyn Fn(&Progress) -> bool + Send + Sync;
type AbortCheck = dyn Fn() -> bool + Send + Sync;

/// Generation settings; [`Self::iter`] runs the generation.
///
/// At least one of a loop bound and a vertex bound is required. Defaults: no self-loops, no
/// zero-flow edges (bridges with no external leg on one side), any number of bridges.
pub struct GraphGenerator<F, L> {
    externals: Vec<(L, Stub<F>)>,
    signatures: Vec<Vec<Stub<F>>>,
    loops: (usize, Option<usize>),
    max_vertices: Option<usize>,
    max_bridges: Option<usize>,
    allow_self_loops: bool,
    allow_zero_flow_edges: bool,
    #[cfg(feature = "rayon")]
    parallel: bool,
    filter: Option<Box<Filter<F, L>>>,
    check_heredity: bool,
    progress: Option<Box<ProgressFn>>,
    abort_check: Option<Box<AbortCheck>>,
}

impl<F, L> GraphGenerator<F, L>
where
    F: Clone + Ord + Send + Sync,
    L: Clone + Ord + Send + Sync,
{
    /// External legs as `(label, stub at the vertex)` and the vertex signatures; a node of a
    /// generated graph carries the index of its signature.
    pub fn new(
        externals: impl IntoIterator<Item = (L, Stub<F>)>,
        signatures: impl IntoIterator<Item = Vec<Stub<F>>>,
    ) -> Self {
        GraphGenerator {
            externals: externals.into_iter().collect(),
            signatures: signatures.into_iter().collect(),
            loops: (0, None),
            max_vertices: None,
            max_bridges: None,
            allow_self_loops: false,
            allow_zero_flow_edges: false,
            #[cfg(feature = "rayon")]
            parallel: false,
            filter: None,
            check_heredity: false,
            progress: None,
            abort_check: None,
        }
    }

    /// Loop orders to generate, e.g. `2..=2` for exactly two loops or `..=3` up to three.
    pub fn loops(mut self, loops: impl RangeBounds<usize>) -> Self {
        let low = match loops.start_bound() {
            Bound::Included(&l) => l,
            Bound::Excluded(&l) => l + 1,
            Bound::Unbounded => 0,
        };
        let high = match loops.end_bound() {
            Bound::Included(&h) => Some(h),
            Bound::Excluded(&h) => Some(h.saturating_sub(1)),
            Bound::Unbounded => None,
        };
        self.loops = (low, high);
        self
    }

    /// Maximum number of vertices (legs are not vertices).
    pub fn max_vertices(mut self, max_vertices: usize) -> Self {
        self.max_vertices = Some(max_vertices);
        self
    }

    /// Maximum number of bridges; `0` keeps one-particle-irreducible graphs. Bridges are paired
    /// edges whose removal disconnects the graph; legs are never bridges.
    pub fn max_bridges(mut self, max_bridges: usize) -> Self {
        self.max_bridges = Some(max_bridges);
        self
    }

    /// Allow edges from a vertex to itself.
    pub fn allow_self_loops(mut self, allow: bool) -> Self {
        self.allow_self_loops = allow;
        self
    }

    /// Allow bridges that carry no momentum: one side has no external leg (tadpoles, and every
    /// bridge of a vacuum graph).
    ///
    /// Only paired edges count, so a one-leg process keeps its graphs without bridges (at one
    /// loop in phi^3, the vertex with the leg and a self-loop) although its leg carries no
    /// momentum either. Graphica counts that leg as a zero-flow bridge and returns no one-leg
    /// graph without zero-flow edges or with `max_bridges(0)`.
    pub fn allow_zero_flow_edges(mut self, allow: bool) -> Self {
        self.allow_zero_flow_edges = allow;
        self
    }

    /// Keep only the graphs accepted by `filter`, and prune the search with it.
    ///
    /// `filter` sees every state of the search as a [`GraphView`]. Stage 1 builds untyped
    /// topologies ([`Stage::Topology`]) up to each complete topology: a vertex at a time, the
    /// pendant vertices that complete a topology at once, then the labels of anonymous legs one at
    /// a time (see [`Stage::Topology`]). Stage 2 assigns fields to the edges of an accepted
    /// complete topology ([`Stage::Fields`]) up to each complete graph ([`Stage::Complete`]), seen
    /// before it is canonized. The output is the graphs whose complete topology and complete graph
    /// `filter` accepts.
    ///
    /// A rejected state is dropped with everything that completes it, so `filter` must be
    /// hereditary within each stage: when it rejects a partial topology, it rejects every complete
    /// topology completing it (with more vertices, edges and legs), and when it rejects a partial
    /// field assignment, it rejects every complete graph completing it (with more fields and the
    /// vertex rules). Each state costs a call, so a condition is best tested where it is decided:
    /// one on the topology alone in stage 1, accepting every state from [`Stage::Fields`] on, and
    /// one on complete graphs only at [`Stage::Complete`], accepting every earlier state.
    /// [`Self::check_heredity`] validates a filter.
    ///
    /// The view shows each graph in one vertex numbering, so the verdict on complete graphs
    /// should not depend on it.
    pub fn filter(
        mut self,
        filter: impl Fn(&GraphView<'_, F, L>) -> bool + Send + Sync + 'static,
    ) -> Self {
        self.filter = Some(Box::new(filter));
        self
    }

    /// Validate the heredity of the filter instead of pruning with it: search every partial state,
    /// and panic when the filter accepts a complete topology or a complete graph although it
    /// rejected a partial state of the same stage that it completes. The output is that of a
    /// hereditary filter, at the cost of the search without the partial prunes.
    pub fn check_heredity(mut self, check: bool) -> Self {
        self.check_heredity = check;
        self
    }

    /// Called after every complete topology (before its graphs are decorated) and every emitted
    /// graph; returning `false` stops the generation. A parallel search counts the topologies of
    /// a subtree when its graphs are handed out.
    pub fn progress(
        mut self,
        progress: impl Fn(&Progress) -> bool + Send + Sync + 'static,
    ) -> Self {
        self.progress = Some(Box::new(progress));
        self
    }

    /// Polled during the search; returning `true` stops the generation.
    pub fn abort_check(mut self, abort_check: impl Fn() -> bool + Send + Sync + 'static) -> Self {
        self.abort_check = Some(Box::new(abort_check));
        self
    }

    /// Search independent subtrees on the current rayon thread pool. The output order is the
    /// same as without parallelism.
    #[cfg(feature = "rayon")]
    pub fn parallel(mut self, parallel: bool) -> Self {
        self.parallel = parallel;
        self
    }

    /// Generate lazily: each call to `next` searches up to the next graph.
    pub fn iter(&self) -> Result<Diagrams<'_, F, L>, GenerationError> {
        let compiled = Compiled::new(self)?;
        let n_legs = compiled.model.legs.total as usize;
        let degrees = &compiled.model.degrees;
        let mut orders = VecDeque::new();
        // without signatures there is no vertex, hence no graph
        if let (Some(&min_degree), Some(&max_degree)) = (degrees.first(), degrees.last()) {
            let (min_degree, max_degree) = (min_degree as usize, max_degree as usize);
            if min_degree <= 2 && self.max_vertices.is_none() {
                return Err(GenerationError::Unbounded);
            }
            let high = match (self.loops.1, self.max_vertices) {
                (Some(high), _) => high,
                // 2 L - 2 + legs = sum (degree - 2) <= max_vertices (max_degree - 2)
                (None, Some(v)) => {
                    (v * max_degree.saturating_sub(2) + 2).saturating_sub(n_legs) / 2
                }
                (None, None) => return Err(GenerationError::Unbounded),
            };
            for loops in self.loops.0..=high {
                let budget = 2 * loops as i32 - 2 + n_legs as i32;
                let by_budget = (min_degree > 2).then(|| budget.max(0) as usize / (min_degree - 2));
                // one of them is set: low degrees need a vertex bound
                let max_vertices = [self.max_vertices, by_budget]
                    .into_iter()
                    .flatten()
                    .min()
                    .ok_or(GenerationError::Unbounded)?;
                if max_vertices > MAX_VERTICES {
                    return Err(GenerationError::TooManyVertices);
                }
                let limits = Limits::new(
                    &compiled.model,
                    budget,
                    max_vertices,
                    self.max_bridges,
                    self.allow_zero_flow_edges,
                    self.allow_self_loops,
                );
                orders.push_back((loops, limits));
            }
        }
        Ok(Diagrams {
            generator: self,
            compiled,
            orders,
            current: None,
            run: None,
            scratch: None,
            tasks: VecDeque::new(),
            window: VecDeque::new(),
            alone: true,
            ready: VecDeque::new(),
            stop: AtomicBool::new(false),
            progress: Progress::default(),
        })
    }
}

/// Generation input in dense form, with the payloads needed to materialize graphs.
struct Compiled<F, L> {
    model: Model,
    /// Distinct fields, sorted: field ids are ranks.
    fields: Vec<F>,
    /// Distinct leg labels, sorted.
    labels: Vec<L>,
    /// Per leg class, in the model's class order.
    legs: Vec<LegClass>,
}

/// Interchangeable external legs: equal label and stub.
#[derive(Clone, Copy, Debug)]
struct LegClass {
    /// Index into the sorted distinct labels.
    label: u32,
    field: u8,
    flow: Option<Flow>,
}

/// Buffers for the flat form of one diagram.
#[derive(Debug, Default)]
struct RawDiagram {
    colours: Vec<u32>,
    hedges: Vec<RawHedge>,
    rules: Vec<u32>,
    fields: Vec<u8>,
}

impl<F: Clone + Ord, L: Clone + Ord> Compiled<F, L> {
    fn new(generator: &GraphGenerator<F, L>) -> Result<Self, GenerationError> {
        let mut fields: Vec<F> = generator
            .signatures
            .iter()
            .flatten()
            .chain(generator.externals.iter().map(|(_, s)| s))
            .map(|s| s.field.clone())
            .collect();
        fields.sort();
        fields.dedup();
        if fields.len() > model::MAX_FIELDS {
            return Err(GenerationError::TooManyFields);
        }
        let field = |s: &Stub<F>| fields.binary_search(&s.field).expect("interned field") as u8;
        let stub = |s: &Stub<F>| 3 * field(s) + s.direction();
        let signatures: Vec<Vec<u8>> = generator
            .signatures
            .iter()
            .map(|sig| sig.iter().map(stub).collect())
            .collect();
        let mut externals: Vec<(&L, u8, &Stub<F>)> = generator
            .externals
            .iter()
            .map(|(l, s)| (l, stub(s), s))
            .collect();
        externals.sort_by(|a, b| (a.0, a.1).cmp(&(b.0, b.1)));
        let mut labels: Vec<L> = externals.iter().map(|(l, ..)| (*l).clone()).collect();
        labels.dedup();
        let mut classes = vec![];
        let mut legs = vec![];
        for run in externals.chunk_by(|a, b| (a.0, a.1) == (b.0, b.1)) {
            let (label, stub_type, s) = run[0];
            let multiplicity = u8::try_from(run.len()).map_err(|_| GenerationError::TooManyLegs)?;
            classes.push((stub_type, multiplicity));
            legs.push(LegClass {
                label: labels.binary_search(label).expect("interned label") as u32,
                field: field(s),
                flow: s.flow,
            });
        }
        Ok(Compiled {
            model: Model::new(fields.len(), &signatures, &classes)?,
            fields,
            labels,
            legs,
        })
    }

    /// Filter, materialize and certify a decorated topology.
    fn diagram(
        &self,
        d: Decoration<'_>,
        loops: usize,
        filter: Option<&Filter<F, L>>,
        raw: &mut RawDiagram,
    ) -> Option<Diagram<F, L>> {
        let mut view = GraphView::new(self, loops, State::Decoration(d));
        if let Some(filter) = filter {
            if !filter(&view) {
                return None;
            }
            assert!(
                !d.vetoed,
                "the filter is not hereditary: it accepts a graph but rejected a partial field \
                 assignment it completes"
            );
        }
        let graph = view.take_complete_graph();
        let (certificate, automorphism_count) = self.certify(&d, raw);
        Some(Diagram {
            graph,
            loops,
            automorphism_count,
            certificate,
        })
    }

    /// The certificate of the materialized graph under `canonize_by(|&s| s, F::clone,
    /// Option::clone)`, from the decoration directly: node keys are signatures, edge keys fields
    /// and half-edge keys `None` for paired hedges and the label for legs.
    fn certify(
        &self,
        d: &Decoration<'_>,
        raw: &mut RawDiagram,
    ) -> (Certificate<usize, F, Option<L>>, Option<u128>) {
        let t = d.topology;
        let model = &self.model;
        raw.rules.clear();
        raw.rules.extend_from_slice(d.rules);
        raw.rules.sort_unstable();
        raw.rules.dedup();
        raw.colours.clear();
        raw.colours.extend(
            d.rules
                .iter()
                .map(|r| raw.rules.binary_search(r).expect("present rule") as u32),
        );
        let edge_field = |(e, &l): (&Edge, &u8)| {
            if e.u == e.v {
                model.loop_labels[l as usize].field
            } else {
                model.labels[l as usize].field
            }
        };
        raw.fields.clear();
        raw.fields.extend(self.legs.iter().map(|leg| leg.field));
        raw.fields
            .extend(d.edges.iter().zip(d.labels).map(edge_field));
        raw.fields.sort_unstable();
        raw.fields.dedup();
        let fields = &raw.fields;
        let field_rank = |f: u8| fields.binary_search(&f).expect("present field") as u32;
        // `None` (paired hedges) sorts before every label
        let paired = !d.edges.is_empty();
        let hedges = &mut raw.hedges;
        hedges.clear();
        let mut n_edges = 0;
        for v in 0..t.n {
            let labelled = t.vertices[v].labelled;
            for (c, leg) in self.legs.iter().enumerate() {
                let role = match leg.flow {
                    None => Role::Leg,
                    Some(Flow::Sink) => Role::LegIn,
                    Some(Flow::Source) => Role::LegOut,
                };
                for _ in 0..model.legs.count(labelled, c) {
                    let h = hedges.len() as u32;
                    hedges.push(RawHedge {
                        node: v as u32,
                        opposite: h,
                        edge: n_edges,
                        role,
                        edge_key: field_rank(leg.field),
                        hedge_key: leg.label + paired as u32,
                    });
                    n_edges += 1;
                }
            }
        }
        for (edge, &l) in d.edges.iter().zip(d.labels) {
            let (field, roles) = if edge.u == edge.v {
                let label = model.loop_labels[l as usize];
                let roles = if label.directed {
                    (Role::Tail, Role::Head)
                } else {
                    (Role::Undirected, Role::Undirected)
                };
                (label.field, roles)
            } else {
                let label = model.labels[l as usize];
                let roles = match label.direction {
                    Direction::Undirected => (Role::Undirected, Role::Undirected),
                    Direction::Forward => (Role::Tail, Role::Head),
                    Direction::Backward => (Role::Head, Role::Tail),
                };
                (label.field, roles)
            };
            let h = hedges.len() as u32;
            let edge_key = field_rank(field);
            for (node, opposite, role) in [(edge.u, h + 1, roles.0), (edge.v, h, roles.1)] {
                hedges.push(RawHedge {
                    node: node as u32,
                    opposite,
                    edge: n_edges,
                    role,
                    edge_key,
                    hedge_key: 0,
                });
            }
            n_edges += 1;
        }
        let keys = KeyTables {
            nodes: raw.rules.iter().map(|&r| r as usize).collect(),
            edges: raw
                .fields
                .iter()
                .map(|&f| self.fields[f as usize].clone())
                .collect(),
            hedges: paired
                .then_some(None)
                .into_iter()
                .chain(self.labels.iter().cloned().map(Some))
                .collect(),
        };
        Certificate::of_hedges(&raw.colours, &raw.hedges, n_edges as usize, keys)
    }
}

/// Search and decoration buffers.
#[derive(Debug, Default)]
struct Scratch {
    search: topology::Buffers,
    decorator: Decorator,
    raw: RawDiagram,
}

/// Spare buffers per thread, for later generations and the subtrees of a parallel search.
const POOL: usize = 16;

thread_local! {
    static SPARE: RefCell<Vec<Scratch>> = RefCell::default();
}

impl Scratch {
    fn take() -> Scratch {
        SPARE
            .with(|spare| spare.borrow_mut().pop())
            .unwrap_or_default()
    }

    fn give(self) {
        SPARE.with(|spare| {
            let mut spare = spare.borrow_mut();
            if spare.len() < POOL {
                spare.push(self);
            }
        });
    }
}

/// What a generation step needs besides its own state.
struct RunContext<'a, F, L> {
    compiled: &'a Compiled<F, L>,
    limits: &'a Limits,
    stop: &'a AtomicBool,
    abort_check: Option<&'a AbortCheck>,
    loops: usize,
    filter: Option<Pruner<'a, F, L>>,
}

impl<'a, F: Sync, L: Sync> RunContext<'a, F, L> {
    fn new(
        generator: &'a GraphGenerator<F, L>,
        compiled: &'a Compiled<F, L>,
        current: &'a Option<(usize, Limits)>,
        stop: &'a AtomicBool,
    ) -> Self {
        let (loops, limits) = current.as_ref().expect("a loop order is running");
        RunContext {
            compiled,
            limits,
            stop,
            abort_check: generator.abort_check.as_deref(),
            loops: *loops,
            filter: generator.filter.as_deref().map(|filter| Pruner {
                compiled,
                filter,
                loops: *loops,
                check_heredity: generator.check_heredity,
            }),
        }
    }

    /// The filter as both stages prune with it.
    fn prune(&self) -> Option<&dyn Prune> {
        self.filter.as_ref().map(|filter| filter as &dyn Prune)
    }

    fn search(&self) -> Context<'_> {
        Context {
            model: &self.compiled.model,
            limits: self.limits,
            stop: self.stop,
            abort_check: self.abort_check,
            filter: self.prune(),
        }
    }
}

/// Call the progress callback; `false` stops the generation.
fn report(callback: Option<&ProgressFn>, progress: &Progress, stop: &AtomicBool) {
    if callback.is_some_and(|callback| !callback(progress)) {
        stop.store(true, Ordering::Relaxed);
    }
}

/// A resumable generation under one search object: the topology search and the decoration of
/// its current complete topology.
struct Run {
    search: Search,
    decorator: Decorator,
    raw: RawDiagram,
    decorating: bool,
}

impl Run {
    fn new(model: &Model, task: Topology, scratch: Scratch) -> Self {
        Run {
            search: Search::new(model, task, scratch.search),
            decorator: scratch.decorator,
            raw: scratch.raw,
            decorating: false,
        }
    }

    fn into_scratch(self) -> Scratch {
        Scratch {
            search: self.search.into_buffers(),
            decorator: self.decorator,
            raw: self.raw,
        }
    }

    /// The next diagram; `None` when the subtree is exhausted or the generation stops.
    /// `topology` is called on every complete topology before it is decorated.
    fn next<F: Clone + Ord + Sync, L: Clone + Ord + Sync>(
        &mut self,
        cx: &RunContext<'_, F, L>,
        topology: &mut dyn FnMut(),
    ) -> Option<Diagram<F, L>> {
        let model = &cx.compiled.model;
        let filter = cx.filter.as_ref().map(|pruner| pruner.filter);
        loop {
            while self.decorating {
                let t = self.search.topology();
                let Some(d) = self.decorator.next(model, t, cx.prune()) else {
                    self.decorating = false;
                    break;
                };
                if let Some(diagram) = cx.compiled.diagram(d, cx.loops, filter, &mut self.raw) {
                    return Some(diagram);
                }
            }
            if !self.search.next(&cx.search()) {
                return None;
            }
            topology();
            if cx.stop.load(Ordering::Relaxed) {
                return None;
            }
            self.decorator
                .start(model, self.search.topology(), self.search.generators());
            self.decorating = true;
        }
    }
}

/// Subtrees per thread a parallel search starts from (splitting balances the rest); a window of
/// [`WINDOW_PER_THREAD`] per thread runs at once, each buffering at most [`AHEAD`] diagrams before
/// it is the oldest.
const TASKS_PER_THREAD: usize = 8;
const WINDOW_PER_THREAD: usize = 16;
const AHEAD: usize = 1024;

/// A subtree of a parallel search and the diagrams it found that are not emitted yet.
struct Slot<F, L> {
    task: Option<Topology>,
    run: Option<Run>,
    diagrams: Vec<Diagram<F, L>>,
    topologies: usize,
    done: bool,
}

impl<F: Clone + Ord + Sync, L: Clone + Ord + Sync> Slot<F, L> {
    fn new(task: Topology) -> Self {
        Slot {
            task: Some(task),
            run: None,
            diagrams: vec![],
            topologies: 0,
            done: false,
        }
    }

    /// Whether a step of the window that buffers up to `ahead` diagrams advances it.
    fn active(&self, ahead: usize) -> bool {
        !self.done && self.diagrams.len() < ahead
    }

    /// Hand the rest of the running search's shallowest unfinished level to new subtrees, which
    /// come right after this one in the output.
    fn split(&mut self, cx: &RunContext<'_, F, L>) -> Vec<Topology> {
        let Some(run) = &mut self.run else {
            return vec![];
        };
        let mut scratch = Scratch::take();
        let split = run.search.split(&cx.search(), &mut scratch.search);
        scratch.give();
        split
    }

    /// Search until `ahead` diagrams are buffered or the subtree is exhausted.
    fn advance(&mut self, cx: &RunContext<'_, F, L>, ahead: usize) {
        if !self.active(ahead) {
            return;
        }
        let run = self.run.get_or_insert_with(|| {
            let task = self.task.take().expect("an unstarted slot has its task");
            Run::new(&cx.compiled.model, task, Scratch::take())
        });
        while self.diagrams.len() < ahead {
            let Some(diagram) = run.next(cx, &mut || self.topologies += 1) else {
                self.done = true;
                self.run.take().expect("running").into_scratch().give();
                return;
            };
            self.diagrams.push(diagram);
        }
    }
}

/// Graphs in generation order; see [`GraphGenerator::iter`].
pub struct Diagrams<'g, F, L> {
    generator: &'g GraphGenerator<F, L>,
    compiled: Compiled<F, L>,
    /// Loop orders not started yet.
    orders: VecDeque<(usize, Limits)>,
    current: Option<(usize, Limits)>,
    /// The sequential search of the current loop order, and its buffers between orders.
    run: Option<Run>,
    scratch: Option<Scratch>,
    /// The parallel search of the current loop order: subtrees not started, running subtrees in
    /// order, and diagrams of finished ones not emitted yet.
    tasks: VecDeque<Topology>,
    window: VecDeque<Slot<F, L>>,
    /// Whether the next step runs only the oldest subtree, up to its first diagram: the first step,
    /// so that the first graph does not wait for a step of the whole window.
    alone: bool,
    ready: VecDeque<Diagram<F, L>>,
    stop: AtomicBool,
    progress: Progress,
}

impl<F, L> Diagrams<'_, F, L>
where
    F: Clone + Ord + Send + Sync,
    L: Clone + Ord + Send + Sync,
{
    /// Whether the abort check or the progress callback stopped the generation.
    pub fn aborted(&self) -> bool {
        self.stop.load(Ordering::Relaxed)
    }

    fn parallel(&self) -> bool {
        #[cfg(feature = "rayon")]
        return self.generator.parallel;
        #[cfg(not(feature = "rayon"))]
        false
    }

    /// Map `items` in order, on the current rayon thread pool in a parallel search.
    fn map_ordered<T: Send, R: Send>(
        &self,
        items: Vec<T>,
        f: impl Fn(T) -> R + Send + Sync,
    ) -> Vec<R> {
        #[cfg(feature = "rayon")]
        if self.generator.parallel {
            use rayon::prelude::*;
            return items.into_par_iter().map(f).collect();
        }
        items.into_iter().map(f).collect()
    }

    fn threads(&self) -> usize {
        #[cfg(feature = "rayon")]
        if self.generator.parallel {
            return rayon::current_num_threads();
        }
        1
    }

    fn report(&self) {
        report(
            self.generator.progress.as_deref(),
            &self.progress,
            &self.stop,
        );
    }

    /// Start the next loop order; `false` when all are done.
    fn next_loop_order(&mut self) -> bool {
        let Some((loops, limits)) = self.orders.pop_front() else {
            return false;
        };
        self.progress.loops = loops;
        let root = Topology::empty(&self.compiled.model, &limits);
        self.current = Some((loops, limits));
        if self.parallel() {
            self.tasks = self.frontier(root).into();
        } else {
            let scratch = self.scratch.take().unwrap_or_else(Scratch::take);
            self.run = Some(Run::new(&self.compiled.model, root, scratch));
        }
        true
    }

    /// Split the search into subtrees, breadth-first and in search order, until there are enough
    /// to balance the threads.
    fn frontier(&self, root: Topology) -> Vec<Topology> {
        let target = TASKS_PER_THREAD * self.threads();
        let run = RunContext::new(self.generator, &self.compiled, &self.current, &self.stop);
        let cx = run.search();
        let mut items = vec![root];
        while items.len() < target
            && !self.stop.load(Ordering::Relaxed)
            && items.iter().any(|item| !item.is_complete())
        {
            let children = self.map_ordered(items, |item| {
                if item.is_complete() {
                    return vec![item];
                }
                let mut scratch = Scratch::take();
                let children = Search::children(&cx, item, 0, &mut scratch.search);
                scratch.give();
                children
            });
            items = children.into_iter().flatten().collect();
        }
        items
    }

    /// Run the window of subtrees one step and move the diagrams of its oldest subtrees to
    /// `ready`; `false` when the loop order is done.
    fn step_window(&mut self) -> bool {
        let window = WINDOW_PER_THREAD * self.threads();
        while self.window.len() < window {
            let Some(task) = self.tasks.pop_front() else {
                break;
            };
            self.window.push_back(Slot::new(task));
        }
        if self.window.is_empty() {
            return false;
        }
        let cx = RunContext::new(self.generator, &self.compiled, &self.current, &self.stop);
        // until the first diagram is out, a step buffers one diagram per subtree
        let ahead = if self.progress.diagrams == 0 {
            1
        } else {
            AHEAD
        };
        if std::mem::take(&mut self.alone) {
            self.window[0].advance(&cx, ahead);
        } else {
            // a subtree that holds the output back runs alone: split the rest of the running
            // searches, oldest first, until every thread has a subtree that is neither finished nor
            // full, or the window is full (each slot buffers diagrams; the last split may add the
            // children of one level beyond it)
            let threads = self.threads();
            let mut i = 0;
            while i < self.window.len()
                && self.window.len() < window
                && self.window.iter().filter(|slot| slot.active(AHEAD)).count() < threads
            {
                let split = self.window[i].split(&cx);
                i += 1;
                for (k, task) in split.into_iter().enumerate() {
                    self.window.insert(i + k, Slot::new(task));
                }
            }
            #[cfg(feature = "rayon")]
            {
                use rayon::prelude::*;
                self.window
                    .par_iter_mut()
                    .for_each(|slot| slot.advance(&cx, ahead));
            }
            #[cfg(not(feature = "rayon"))]
            self.window
                .iter_mut()
                .for_each(|slot| slot.advance(&cx, ahead));
        }
        while let Some(oldest) = self.window.front_mut() {
            self.ready.extend(oldest.diagrams.drain(..));
            let topologies = std::mem::take(&mut oldest.topologies);
            let done = oldest.done;
            if done {
                self.window.pop_front();
            }
            if topologies > 0 {
                // a subtree's topologies are counted when its diagrams are handed out
                self.progress.topologies += topologies;
                self.report();
            }
            if !done {
                break;
            }
        }
        true
    }
}

impl<F, L> Drop for Diagrams<'_, F, L> {
    /// Hand the buffers of the sequential search to the next generation on this thread.
    fn drop(&mut self) {
        if let Some(scratch) = self
            .run
            .take()
            .map(Run::into_scratch)
            .or(self.scratch.take())
        {
            scratch.give();
        }
    }
}

impl<F, L> Iterator for Diagrams<'_, F, L>
where
    F: Clone + Ord + Send + Sync,
    L: Clone + Ord + Send + Sync,
{
    type Item = Diagram<F, L>;

    fn next(&mut self) -> Option<Diagram<F, L>> {
        loop {
            if self.stop.load(Ordering::Relaxed) {
                return None;
            }
            if let Some(diagram) = self.ready.pop_front() {
                self.progress.diagrams += 1;
                self.report();
                return Some(diagram);
            }
            if let Some(run) = &mut self.run {
                let cx = RunContext::new(self.generator, &self.compiled, &self.current, &self.stop);
                let (callback, progress, stop) = (
                    self.generator.progress.as_deref(),
                    &mut self.progress,
                    &self.stop,
                );
                let diagram = run.next(&cx, &mut || {
                    progress.topologies += 1;
                    report(callback, progress, stop);
                });
                match diagram {
                    Some(diagram) => {
                        self.progress.diagrams += 1;
                        self.report();
                        return Some(diagram);
                    }
                    None => {
                        let run = self.run.take().expect("running");
                        self.scratch = Some(run.into_scratch());
                        self.current = None;
                    }
                }
            } else if self.current.is_none() || !self.step_window() {
                self.current = None;
                if !self.next_loop_order() {
                    return None;
                }
            }
        }
    }
}
