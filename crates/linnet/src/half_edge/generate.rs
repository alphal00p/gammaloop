//! # Feynman-graph generation
//!
//! [`Plan`] enumerates the connected graphs built from vertex signatures (multisets of
//! typed stubs) and a set of external legs, one per isomorphism class, each with its exact
//! automorphism count and the canonical form of the native canonizer.
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
//! stub. The canonical form and automorphism count are those of
//! `graph.canonize_by(|&s| s, |f| f.clone(), |h| h.clone())`, so symmetry factors follow the
//! canonizer's half-edge convention: `m!` for `m` interchangeable parallel edges or legs and `2`
//! for an undirected self-loop. Signatures with identical stubs are distinct vertex rules and give
//! distinct graphs.
//!
//! ```
//! use linnet::half_edge::generate::{Theory, Legs, Query, Plan, Run, Structure, Stub};
//!
//! // QED photon self-energy at two loops, one-particle irreducible
//! let (photon, electron) = (Stub::undirected('a'), Stub::incoming('e'));
//! let theory = Theory::builder()
//!     .vertex([electron.clone(), electron.flipped(), photon.clone()], []).build().unwrap();
//! let legs = Legs::new().leg(1, photon.clone()).leg(2, photon);
//! let query = Query::new(&legs).loops(2..=2).structure(Structure::one_pi());
//! let plan = Plan::new(&theory, &query).unwrap();
//! let run = Run::default();
//! let diagrams: Vec<_> = plan.diagrams(&run).unwrap().collect();
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
//! of its automorphism group ([`decorate`]). Budgets ([`Query::weight`],
//! [`Query::count`]) bound every partial decoration and drop those whose bounds miss their
//! ranges. Loop orders without a vertex multiset that meets the weights, and requests whose legs
//! carry a charge or parity that no vertex balances, are not searched. Accepted diagrams own compact data; their [`HedgeGraph`] and canonical form are
//! materialized only when requested (also checked independently in debug builds).
//!
//! Both stages are resumable depth-first searches, so [`Diagrams`] streams: each `next` call
//! searches only up to the next diagram and nothing is buffered. Output order is the depth-first
//! order of the search and does not depend on parallelism: with the `rayon` feature,
//! [`Run::parallel`] runs independent subtrees on the current thread pool, each at most
//! a fixed number of diagrams ahead of the output, and emits their results in that order. The
//! search starts from a few subtrees per thread; while fewer subtrees than threads can advance,
//! the unsearched rest of the oldest ones is split off into new subtrees. The first step runs only
//! the oldest subtree, up to its first diagram, and further steps buffer one diagram per subtree
//! until the first is emitted, so the first graph does not wait for a step of the whole window.
//!
//! ## Filtering
//!
//! [`Query::prune_skeletons`] and [`Query::prune_dressings`] see typed borrowed views
//! and prune each search independently. A pruning predicate must reject every refinement of a
//! rejected state. [`Query::select`] sees only completed assignments and may be arbitrary;
//! [`Run::check_heredity`] diagnoses pruning violations independently of final selection.
//!
//! ```
//! use linnet::half_edge::generate::{Theory, Legs, Query, Plan, Run, Structure, Stub};
//!
//! // gg -> gg at one loop without ghosts: ghost fields are decided edge by edge in stage 2
//! let (g, q, c) = (Stub::undirected('g'), Stub::incoming('q'), Stub::incoming('c'));
//! let signatures = [
//!     vec![q.clone().flipped(), q, g.clone()],
//!     vec![g.clone(); 3],
//!     vec![g.clone(); 4],
//!     vec![c.clone().flipped(), c, g.clone()],
//! ];
//! let theory = signatures.into_iter().fold(Theory::builder(), |t, rule| t.vertex(rule, []))
//!     .build().unwrap();
//! let legs = (1..=4).fold(Legs::new(), |legs, label| legs.leg(label, g.clone()));
//! let query = Query::new(&legs).loops(1..=1)
//!     .structure(Structure { self_loops: true, ..Structure::default() })
//!     .prune_dressings(|view| view.edges().all(|e| e.stub.is_none_or(|s| *s.field != 'c')));
//! let plan = Plan::new(&theory, &query).unwrap();
//! let run = Run::default();
//! let diagrams: Vec<_> = plan.diagrams(&run).unwrap().collect();
//! assert_eq!(diagrams.len(), 111);
//! ```

#[doc(hidden)]
pub mod batch;
mod budget;
mod decorate;
mod model;
mod output;
mod query;
mod suspend;
#[cfg(test)]
mod tests;
mod topology;
mod view;

use std::{
    cell::RefCell,
    collections::{BTreeMap, VecDeque},
    ops::ControlFlow,
    sync::{
        atomic::{AtomicBool, Ordering},
        Arc, OnceLock,
    },
};

use num_bigint::BigUint;
use num_traits::Zero;

pub use budget::Countable;
use budget::Request;
use decorate::{Decoration, Decorator};
use model::{Model, INCOMING, OUTGOING, UNDIRECTED};
pub use output::{Counts, Dressings, LoopCounts, Skeleton, SkeletonSnapshot, Skeletons};
use output::{DiagramOutput, Emit, SkeletonData};
use query::Hooks;
pub use query::{
    Bridges, Legs, PlanCriterion, PlanReport, PlanRole, PlanStage, Query, Run, Structure, Theory,
    TheoryBuilder,
};
use topology::{Context, Limits, Search, Topology, MAX_VERTICES};
pub use view::{
    Bridge, DiagramView, DressingView, EdgeView, OrderBound, SkeletonView, Stage, VertexSet,
};
use view::{Prune, Pruner, State, View};

use super::{
    canon::{CanonicalForm, KeyTables, RawHedge, Role},
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

/// A generated graph, with compact owned data and independently lazy graph/form caches.
#[derive(Clone)]
pub struct Diagram<F, L> {
    rules: Box<[u32]>,
    edges: Box<[DiagramEdge]>,
    /// `(vertex, leg class)` in native half-edge insertion order.
    legs: Box<[(u8, u32)]>,
    tables: Arc<OutputTables<F, L>>,
    skeleton: Arc<SkeletonData<F, L>>,
    graph: OnceLock<HedgeGraph<F, usize, Option<L>>>,
    canonical_form: OnceLock<CanonicalForm<usize, F, Option<L>>>,
    /// Loop order.
    pub loops: usize,
    /// Order of the half-edge automorphism group (the symmetry factor); `None` beyond `u128`.
    pub automorphism_count: Option<u128>,
}

/// Payload ranks for an internal edge; no search buffers or automaton are retained.
#[derive(Clone)]
struct DiagramEdge {
    vertices: [u8; 2],
    stub: Stub<u8>,
}

impl<F, L> std::fmt::Debug for Diagram<F, L> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Diagram")
            .field("rules", &self.rules)
            .field("edges", &self.edges.len())
            .field("legs", &self.legs.len())
            .field("loops", &self.loops)
            .field("automorphism_count", &self.automorphism_count)
            .field("graph_cached", &self.graph.get().is_some())
            .field(
                "canonical_form_cached",
                &self.canonical_form.get().is_some(),
            )
            .finish()
    }
}

impl<F, L> Diagram<F, L> {
    /// Order of the half-edge automorphism group. Its reciprocal contributes to
    /// [`Counts::sum_inv_aut`]; computing either does not materialize a graph.
    pub fn symmetry_factor(&self) -> Result<u128, GenerationError> {
        self.automorphism_count
            .ok_or(GenerationError::AutomorphismOverflow)
    }

    /// The complete untyped topology, independent of the plan and its search buffers.
    pub fn skeleton(&self) -> SkeletonView<'_, F, L> {
        self.skeleton.view()
    }

    /// Retain the same immutable skeleton independently of this diagram and its originating run.
    pub fn skeleton_snapshot(&self) -> SkeletonSnapshot<F, L> {
        SkeletonSnapshot {
            data: Arc::clone(&self.skeleton),
        }
    }

    /// Signature indices in native vertex order, without materializing the graph.
    pub fn rules(&self) -> &[u32] {
        &self.rules
    }

    /// Internal edges in decorator order; every stub is assigned, at the first endpoint.
    pub fn edges(&self) -> impl ExactSizeIterator<Item = EdgeView<'_, F>> {
        self.edges.iter().map(|edge| EdgeView {
            vertices: edge.vertices.map(usize::from),
            stub: Some(Stub {
                field: &self.tables.fields[edge.stub.field as usize],
                flow: edge.stub.flow,
            }),
        })
    }

    /// External slots in insertion order: vertex, leg class, then repetition within the class.
    pub fn legs(&self) -> impl ExactSizeIterator<Item = (usize, &L, Stub<&F>)> {
        self.legs.iter().map(|&(vertex, class)| {
            let leg = self.tables.legs[class as usize];
            (
                vertex as usize,
                &self.tables.labels[leg.label as usize],
                Stub {
                    field: &self.tables.fields[leg.field as usize],
                    flow: leg.flow,
                },
            )
        })
    }
}

impl<F: Clone, L: Clone> Diagram<F, L> {
    /// Fields on edges/legs, signature indices on nodes and `Some(label)` on dangling hedges.
    pub fn graph(&self) -> &HedgeGraph<F, usize, Option<L>> {
        self.graph.get_or_init(|| {
            View::build(
                self.rules.len(),
                self.legs.len(),
                self.legs(),
                self.edges(),
                |v| self.rules[v] as usize,
                |stub| stub.expect("complete edge").field.clone(),
            )
        })
    }
}

impl<F: Clone + Ord, L: Clone + Ord> Diagram<F, L> {
    /// The form of `graph().canonize_by(rule, field, optional_leg_label)`, without building a graph.
    pub fn canonical_form(&self) -> &CanonicalForm<usize, F, Option<L>> {
        self.canonical_form.get_or_init(|| self.certify().0)
    }

    fn certify(&self) -> (CanonicalForm<usize, F, Option<L>>, Option<u128>) {
        let mut raw = RawDiagram::default();
        raw.rules.extend_from_slice(&self.rules);
        raw.rules.sort_unstable();
        raw.rules.dedup();
        raw.colours.extend(
            self.rules
                .iter()
                .map(|rule| raw.rules.binary_search(rule).expect("present rule") as u32),
        );
        raw.fields
            .extend(self.tables.legs.iter().map(|leg| leg.field));
        raw.fields
            .extend(self.edges.iter().map(|edge| edge.stub.field));
        raw.fields.sort_unstable();
        raw.fields.dedup();
        let field_rank = |field| raw.fields.binary_search(&field).expect("present field") as u32;
        // `None` (paired hedges) sorts before every label.
        let paired = !self.edges.is_empty();
        let mut n_edges = 0;
        for &(vertex, class) in &self.legs {
            let leg = self.tables.legs[class as usize];
            let role = match leg.flow {
                None => Role::Leg,
                Some(Flow::Sink) => Role::LegIn,
                Some(Flow::Source) => Role::LegOut,
            };
            let hedge = raw.hedges.len() as u32;
            raw.hedges.push(RawHedge {
                node: vertex as u32,
                opposite: hedge,
                edge: n_edges,
                role,
                edge_key: field_rank(leg.field),
                hedge_key: leg.label + paired as u32,
            });
            n_edges += 1;
        }
        for edge in &self.edges {
            let roles = match edge.stub.flow {
                None => (Role::Undirected, Role::Undirected),
                Some(Flow::Source) => (Role::Tail, Role::Head),
                Some(Flow::Sink) => (Role::Head, Role::Tail),
            };
            let hedge = raw.hedges.len() as u32;
            let edge_key = field_rank(edge.stub.field);
            for (node, opposite, role) in [
                (edge.vertices[0], hedge + 1, roles.0),
                (edge.vertices[1], hedge, roles.1),
            ] {
                raw.hedges.push(RawHedge {
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
            nodes: raw.rules.iter().map(|&rule| rule as usize).collect(),
            edges: raw
                .fields
                .iter()
                .map(|&field| self.tables.fields[field as usize].clone())
                .collect(),
            hedges: paired
                .then_some(None)
                .into_iter()
                .chain(self.tables.labels.iter().cloned().map(Some))
                .collect(),
        };
        CanonicalForm::of_hedges(&raw.colours, &raw.hedges, n_edges as usize, keys)
    }
}

/// Counts passed to the progress callback.
#[derive(Clone, Debug, Default)]
pub struct Progress {
    /// Loop order being generated.
    pub loops: usize,
    /// Complete topologies decorated so far. Loop orders and requests proven empty before the
    /// search (see [`Query::weight`] and [`Plan::diagrams`]) add none.
    pub topologies: BigUint,
    /// Graphs emitted so far.
    pub diagrams: BigUint,
}

#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
pub enum GenerationError {
    #[error(
        "bound the loop order or the number of vertices (vertices of degree two or less always need a vertex bound)"
    )]
    Unbounded,
    #[error("the automorphism order exceeds u128; an exact reciprocal cannot be counted")]
    AutomorphismOverflow,
    #[error("the skeleton belongs to a different plan")]
    ForeignSkeleton,
    #[error("this plan has deferred hooks; execute it with the batch driver")]
    DeferredHooks,
    #[error("the callback reply does not match the pending batch")]
    InvalidReply,
    #[error("signatures need between 1 and 15 stubs, found {0}")]
    UnsupportedDegree(usize),
    #[error("at most 84 distinct fields are supported")]
    TooManyFields,
    #[error("too many external legs to pack")]
    TooManyLegs,
    #[error("graphs of this size would need more than 64 vertices")]
    TooManyVertices,
    #[error("at most 64 propagator, chord and closed-line counts are supported")]
    TooManyCounts,
    #[error("unknown weight column {0:?}")]
    UnknownWeight(String),
    #[error("generation request is infeasible: {0}")]
    Infeasible(PlanReport),
    #[error("weight name {0:?} is declared more than once")]
    DuplicateWeight(String),
    #[error("vertex rule {rule} has {weights} weights, but the theory declares {columns} columns")]
    VertexWeights {
        rule: usize,
        columns: usize,
        weights: usize,
    },
}

type ProgressFn<'a> = dyn Fn(&Progress) -> ControlFlow<()> + Send + Sync + 'a;
type AbortCheck<'a> = dyn Fn() -> bool + Send + Sync + 'a;

/// A validated request, reusable for independent executions against one compiled theory.
/// Output diagrams retain only their compact data and immutable field/label tables.
pub struct Plan<F, L> {
    compiled: Arc<Compiled<F, L>>,
    hooks: Hooks<F, L>,
    orders: Vec<(usize, Arc<Limits>)>,
    report: PlanReport,
}

impl<F: Clone + Ord, L: Clone + Ord> Plan<F, L> {
    pub fn new(theory: &Theory<F>, query: &Query<F, L>) -> Result<Self, GenerationError> {
        let mut report = PlanReport::default();
        report.criterion(
            "loops",
            format!("{:?}", query.loops),
            PlanStage::Planning,
            true,
            PlanRole::Restriction,
        );
        report.criterion(
            "bridges",
            format!("{:?}", query.structure.bridges),
            PlanStage::Skeletons,
            true,
            PlanRole::Restriction,
        );
        report.criterion(
            "self loops",
            query.structure.self_loops.to_string(),
            PlanStage::Skeletons,
            true,
            PlanRole::Restriction,
        );
        report.criterion(
            "zero-flow edges",
            query.structure.zero_flow_edges.to_string(),
            PlanStage::Skeletons,
            true,
            PlanRole::Restriction,
        );
        report.criterion(
            "vertex bound",
            format!("{:?}", query.structure.max_vertices),
            PlanStage::Skeletons,
            true,
            PlanRole::Restriction,
        );
        let compiled = Compiled::new(theory, query, &mut report)?;
        let degrees = &compiled.model.theory.degrees;
        let (Some(&min_degree), Some(&max_degree)) = (degrees.first(), degrees.last()) else {
            return Err(report.infeasible("the theory has no vertex rules"));
        };
        let (min_degree, max_degree) = (min_degree as usize, max_degree as usize);
        let vertex_bound = query.structure.max_vertices;
        if min_degree <= 2 && vertex_bound.is_none() {
            return Err(GenerationError::Unbounded);
        }
        let n_legs = compiled.model.legs.total as usize;
        // 2 L - 2 + legs = sum (degree - 2) <= max_vertices (max_degree - 2).
        let high = match (query.loops.1, vertex_bound) {
            (Some(high), _) => high,
            (None, Some(v)) => {
                v.checked_mul(max_degree.saturating_sub(2))
                    .and_then(|n| n.checked_add(2))
                    .ok_or(GenerationError::TooManyVertices)?
                    .saturating_sub(n_legs)
                    / 2
            }
            (None, None) => return Err(GenerationError::Unbounded),
        };
        if !compiled.model.balanced {
            return Err(report.infeasible(
                "the external legs violate a conserved charge, parity or stub direction",
            ));
        }
        let budgets = &compiled.model.budgets;
        if budgets.unsatisfiable() {
            return Err(report.infeasible("a budget range has no possible value"));
        }
        let mut orders = Vec::new();
        for loops in query.loops.0..=high {
            let budget = i32::try_from(loops)
                .ok()
                .and_then(|l| l.checked_mul(2))
                .and_then(|n| n.checked_add(n_legs as i32))
                .and_then(|n| n.checked_sub(2))
                .ok_or(GenerationError::TooManyVertices)?;
            let by_budget = (min_degree > 2).then(|| budget.max(0) as usize / (min_degree - 2));
            let max_vertices = [vertex_bound, by_budget]
                .into_iter()
                .flatten()
                .min()
                .ok_or(GenerationError::Unbounded)?;
            if max_vertices > MAX_VERTICES {
                return Err(GenerationError::TooManyVertices);
            }
            if max_vertices == 0 || !budgets.partition_exists(budget, max_vertices) {
                report.empty_partitions.push((
                    loops,
                    "no degree/weight partition satisfies the request".into(),
                ));
                continue;
            }
            orders.push((
                loops,
                Arc::new(Limits::new(
                    &compiled.model,
                    budget,
                    max_vertices,
                    query.structure.bridges.limit(loops),
                    query.structure.zero_flow_edges,
                    query.structure.self_loops,
                )),
            ));
        }
        if orders.is_empty() {
            return Err(report.infeasible("no loop partition satisfies the request"));
        }
        for (present, name, stage, hereditary, role) in [
            (
                query.hooks.skeletons.is_some() || query.hooks.deferred.skeletons,
                "prune_skeletons",
                PlanStage::Skeletons,
                true,
                PlanRole::Restriction,
            ),
            (
                query.hooks.dressings.is_some() || query.hooks.deferred.dressings,
                "prune_dressings",
                PlanStage::Dressing,
                true,
                PlanRole::Restriction,
            ),
            (
                query.hooks.select.is_some() || query.hooks.deferred.select,
                "select",
                PlanStage::Selection,
                false,
                PlanRole::FinalCheck,
            ),
        ] {
            if present {
                report.criterion(name, "user predicate", stage, hereditary, role);
            }
        }
        Ok(Self {
            compiled: Arc::new(compiled),
            hooks: query.hooks.clone(),
            orders,
            report,
        })
    }

    /// The report produced while lowering this request, without running another planner.
    pub fn explain(&self) -> &PlanReport {
        &self.report
    }
}

impl<F: Clone + Ord + Send + Sync, L: Clone + Ord + Send + Sync> Plan<F, L> {
    /// Enumerate lazily, borrowing transient controls only for this execution.
    pub fn diagrams<'a>(&'a self, run: &'a Run<'a>) -> Result<Diagrams<'a, F, L>, GenerationError> {
        self.direct_hooks()?;
        Ok(Diagrams {
            traversal: Traversal::new(self, run),
        })
    }

    fn direct_hooks(&self) -> Result<(), GenerationError> {
        if self.hooks.deferred.any() {
            return Err(GenerationError::DeferredHooks);
        }
        Ok(())
    }
}

/// Generation input in dense form; emitted diagrams retain only the immutable output tables.
struct Compiled<F, L> {
    model: Arc<Model>,
    tables: Arc<OutputTables<F, L>>,
    /// Bindings translate key ranks through Python-owned tables without cloning automata.
    ranked: OnceLock<Arc<Compiled<u32, u32>>>,
}

struct OutputTables<F, L> {
    /// Distinct fields and labels, sorted: their ids are ranks.
    fields: Arc<[F]>,
    labels: Box<[L]>,
    packing: model::Legs,
    signatures: Arc<[Vec<u8>]>,
    weights: Arc<[(String, Vec<i32>)]>,
    /// Per leg class, in the model's class order.
    legs: Box<[LegClass]>,
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
    fn new(
        theory: &Theory<F>,
        query: &Query<F, L>,
        report: &mut PlanReport,
    ) -> Result<Self, GenerationError> {
        let fields = &theory.fields;
        let mut externals = Vec::with_capacity(query.legs.len());
        for (label, stub) in &query.legs.slots {
            let Ok(field) = fields.binary_search(&stub.field) else {
                return Err(report.infeasible("an external field is absent from the theory"));
            };
            externals.push((label, 3 * field as u8 + stub.direction(), stub));
        }
        let field = |stub: &Stub<F>| {
            fields
                .binary_search(&stub.field)
                .expect("validated external field") as u8
        };
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
        let request = Self::request(theory, query, report)?;
        let model = Model::new(Arc::clone(&theory.dense), &classes, &request)?;
        let packing = model.legs.clone();
        Ok(Self {
            model: Arc::new(model),
            ranked: OnceLock::new(),
            tables: Arc::new(OutputTables {
                fields: Arc::clone(fields),
                labels: labels.into_boxed_slice(),
                packing,
                signatures: Arc::clone(&theory.dense.signatures),
                weights: Arc::clone(&theory.weights),
                legs: legs.into_boxed_slice(),
            }),
        })
    }

    /// Budgets with theory fields as ranks: vertex counts become weights,
    /// and fields that no stub carries are left out of the counts. Budgets that every graph meets
    /// are left out, so they cost nothing.
    fn request(
        theory: &Theory<F>,
        query: &Query<F, L>,
        report: &mut PlanReport,
    ) -> Result<Request, GenerationError> {
        let n_signatures = theory.n_rules;
        let fields = &theory.fields;
        // a sum of values of one sign, with the range unbounded on that side of zero
        let bounds_nothing = |values: &[i32], [min, max]: [i64; 2]| {
            (min <= 0 && max == i64::MAX && values.iter().all(|&v| v >= 0))
                || (max >= 0 && min == i64::MIN && values.iter().all(|&v| v <= 0))
        };
        let mut request = Request::default();
        for (name, range) in &query.weights {
            let (_, weights) = theory
                .weights
                .iter()
                .find(|(column, _)| column == name)
                .ok_or_else(|| GenerationError::UnknownWeight(name.clone()))?;
            report.criterion(
                format!("weight {name}"),
                format!("{}..={}", range[0], range[1]),
                PlanStage::Dressing,
                true,
                PlanRole::Propagator,
            );
            report.criterion(
                format!("weight {name}"),
                format!("{}..={}", range[0], range[1]),
                PlanStage::Rules,
                false,
                PlanRole::FinalCheck,
            );
            request.weights.push((weights.clone(), *range));
        }
        let set = |counted: &[F]| {
            (counted.iter())
                .filter_map(|f| fields.binary_search(f).ok())
                .fold(0u128, |set, field| set | 1 << field)
        };
        for (count, range) in &query.counts {
            let criterion = match count {
                Countable::Vertices(rules) => format!("vertices in {rules:?}"),
                Countable::Propagators(_) => "propagators".into(),
                Countable::Chords(_) => "chords".into(),
                Countable::ClosedLines(_) => "closed lines".into(),
            };
            report.criterion(
                criterion,
                format!("{}..={}", range[0], range[1]),
                PlanStage::Dressing,
                true,
                PlanRole::Propagator,
            );
            match count {
                Countable::Vertices(signatures) => request.weights.push((
                    (0..n_signatures)
                        .map(|s| signatures.contains(&s) as i32)
                        .collect(),
                    *range,
                )),
                // counts are sums of ones
                _ if bounds_nothing(&[1], *range) => {}
                Countable::Propagators(counted) => {
                    request.edges.push((false, set(counted), *range))
                }
                Countable::Chords(counted) => request.edges.push((true, set(counted), *range)),
                Countable::ClosedLines(counted) => request.lines.push((set(counted), *range)),
            }
        }
        request
            .weights
            .retain(|(weights, range)| !bounds_nothing(weights, *range));
        Ok(request)
    }

    /// Freeze only an accepted output, transferring a graph only if a hook requested it.
    fn diagram(
        &self,
        d: Decoration<'_>,
        loops: usize,
        mut view: View<'_, F, L>,
        skeleton: Arc<SkeletonData<F, L>>,
    ) -> Diagram<F, L> {
        let diagram = Diagram {
            rules: d.rules.into(),
            edges: d
                .edges
                .iter()
                .zip(d.labels)
                .map(|(edge, &label)| self.edge(edge, label))
                .collect(),
            legs: (0..d.topology.n)
                .flat_map(|vertex| {
                    self.tables
                        .legs
                        .iter()
                        .enumerate()
                        .flat_map(move |(class, _)| {
                            std::iter::repeat_n(
                                (vertex as u8, class as u32),
                                self.model
                                    .legs
                                    .count(d.topology.vertices[vertex].labelled, class)
                                    as usize,
                            )
                        })
                })
                .collect(),
            tables: Arc::clone(&self.tables),
            skeleton,
            graph: view
                .take_complete_graph()
                .map_or_else(OnceLock::new, OnceLock::from),
            canonical_form: OnceLock::new(),
            loops,
            automorphism_count: d.automorphism_count,
        };
        #[cfg(debug_assertions)]
        {
            // Independent native-graph oracle: neither the certificate nor group count is trusted.
            let canon = diagram
                .graph()
                .canonize_by(|&rule| rule, F::clone, Option::clone);
            assert!(canon.canonical_form() == diagram.canonical_form());
            assert_eq!(canon.automorphism_count(), diagram.automorphism_count);
        }
        diagram
    }
}

/// Search and decoration buffers.
#[derive(Debug, Default)]
struct Scratch {
    search: topology::Buffers,
    decorator: Decorator,
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
    abort_check: Option<&'a AbortCheck<'a>>,
    loops: usize,
    filter: Option<Pruner<'a, F, L>>,
    hooks: &'a Hooks<F, L>,
}

impl<'a, F: Send + Sync, L: Send + Sync> RunContext<'a, F, L> {
    fn new(
        plan: &'a Plan<F, L>,
        controls: &'a Run<'a>,
        current: &'a Option<(usize, Arc<Limits>)>,
        stop: &'a AtomicBool,
    ) -> Self {
        let (loops, limits) = current.as_ref().expect("a loop order is running");
        let compiled = &plan.compiled;
        RunContext {
            compiled,
            limits,
            stop,
            abort_check: controls.abort_check.as_deref(),
            loops: *loops,
            filter: (plan.hooks.skeletons.is_some()
                || plan.hooks.dressings.is_some()
                || plan.hooks.deferred.skeletons
                || plan.hooks.deferred.dressings)
                .then_some(Pruner {
                    compiled,
                    hooks: &plan.hooks,
                    loops: *loops,
                    check_heredity: controls.check_heredity,
                }),
            hooks: &plan.hooks,
        }
    }

    fn stopped(&self) -> bool {
        if self.stop.load(Ordering::Relaxed) {
            return true;
        }
        if self.abort_check.is_some_and(|check| check()) {
            self.stop.store(true, Ordering::Relaxed);
            return true;
        }
        false
    }

    /// Complete-assignment acceptance shared by diagram, count and retained-skeleton runs.
    /// Cancellation inside a predicate stops before another predicate or output is delivered.
    fn accept<'b>(&'b self, d: Decoration<'b>) -> Option<View<'b, F, L>> {
        if self.stopped() {
            return None;
        }
        let mut view = View::new(self.compiled, self.loops, State::Decoration(d));
        if let Some(prune) = &self.hooks.dressings {
            let dressing = DressingView(view);
            let keep = prune(&dressing);
            if self.stopped() {
                return None;
            }
            assert!(
                !(keep && d.vetoed),
                "prune_dressings is not hereditary: it accepts a diagram but rejected an ancestor assignment"
            );
            if !keep {
                return None;
            }
            view = dressing.0;
        }
        if let Some(select) = &self.hooks.select {
            let diagram = DiagramView(view);
            let keep = select(&diagram);
            if self.stopped() || !keep {
                return None;
            }
            view = diagram.0;
        }
        Some(view)
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

/// Call the progress callback; `Break` stops the generation.
fn report(callback: Option<&ProgressFn<'_>>, progress: &Progress, stop: &AtomicBool) {
    if callback.is_some_and(|callback| callback(progress).is_break()) {
        stop.store(true, Ordering::Relaxed);
    }
}

/// A resumable generation under one search object: the topology search and the decoration of
/// its current complete topology.
struct Cursor<F, L> {
    search: Search,
    decorator: Decorator,
    decorating: bool,
    skeleton: Option<Arc<SkeletonData<F, L>>>,
}

impl<F, L> Cursor<F, L> {
    fn new(model: &Model, task: Topology, scratch: Scratch) -> Self {
        Cursor {
            search: Search::new(model, task, scratch.search),
            decorator: scratch.decorator,
            decorating: false,
            skeleton: None,
        }
    }

    fn into_scratch(self) -> Scratch {
        Scratch {
            search: self.search.into_buffers(),
            decorator: self.decorator,
        }
    }

    /// The next diagram; `None` when the subtree is exhausted or the generation stops.
    /// `topology` is called on every complete topology before it is decorated.
    fn next<M: Emit<F, L>>(
        &mut self,
        cx: &RunContext<'_, F, L>,
        topology: &mut dyn FnMut(),
    ) -> Option<M::Item>
    where
        F: Clone + Ord + Send + Sync,
        L: Clone + Ord + Send + Sync,
    {
        let model = &cx.compiled.model;
        loop {
            if cx.stopped() {
                return None;
            }
            while self.decorating {
                if cx.stopped() {
                    return None;
                }
                let t = self.search.topology();
                let Some(d) = self.decorator.next(model, t, cx.prune(), || cx.stopped()) else {
                    self.decorating = false;
                    break;
                };
                if let Some(view) = cx.accept(d) {
                    return Some(M::diagram(
                        cx.compiled,
                        d,
                        cx.loops,
                        view,
                        self.skeleton.as_ref(),
                    ));
                }
            }
            if !self.search.next(&cx.search()) {
                return None;
            }
            topology();
            if cx.stop.load(Ordering::Relaxed) {
                return None;
            }
            if M::SNAPSHOT {
                self.skeleton = Some(Arc::new(SkeletonData::new(
                    cx.compiled,
                    &self.search,
                    cx.loops,
                )));
            }
            if !M::DECORATE {
                return M::skeleton(Arc::clone(
                    self.skeleton.as_ref().expect("skeleton output snapshots"),
                ));
            }
            self.decorator.start(
                model,
                self.search.topology(),
                self.search.generators(),
                self.search.node_group_order(),
            );
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
struct Slot<F, L, M: Emit<F, L>> {
    task: Option<Topology>,
    run: Option<Cursor<F, L>>,
    diagrams: Vec<M::Item>,
    topologies: BigUint,
    done: bool,
}

impl<F: Clone + Ord + Send + Sync, L: Clone + Ord + Send + Sync, M: Emit<F, L>> Slot<F, L, M> {
    fn new(task: Topology) -> Self {
        Slot {
            task: Some(task),
            run: None,
            diagrams: vec![],
            topologies: BigUint::zero(),
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
            Cursor::new(&cx.compiled.model, task, Scratch::take())
        });
        while self.diagrams.len() < ahead {
            let Some(diagram) = run.next::<M>(cx, &mut || self.topologies += 1u8) else {
                self.done = true;
                self.run.take().expect("running").into_scratch().give();
                return;
            };
            self.diagrams.push(diagram);
        }
    }
}

/// Common bounded traversal, parameterized only at output acceptance.
struct Traversal<'g, F, L, M: Emit<F, L>> {
    plan: &'g Plan<F, L>,
    controls: &'g Run<'g>,
    /// Loop orders not started yet.
    orders: VecDeque<(usize, Arc<Limits>)>,
    current: Option<(usize, Arc<Limits>)>,
    /// The sequential search of the current loop order, and its buffers between orders.
    run: Option<Cursor<F, L>>,
    scratch: Option<Scratch>,
    /// The parallel search of the current loop order: subtrees not started, running subtrees in
    /// order, and diagrams of finished ones not emitted yet.
    tasks: VecDeque<Topology>,
    window: VecDeque<Slot<F, L, M>>,
    /// Whether the next step runs only the oldest subtree, up to its first diagram: the first step,
    /// so that the first graph does not wait for a step of the whole window.
    alone: bool,
    ready: VecDeque<M::Item>,
    stop: AtomicBool,
    progress: Progress,
    topologies: BTreeMap<usize, BigUint>,
}

impl<'a, F, L, M: Emit<F, L>> Traversal<'a, F, L, M>
where
    F: Clone + Ord + Send + Sync,
    L: Clone + Ord + Send + Sync,
{
    fn new(plan: &'a Plan<F, L>, controls: &'a Run<'a>) -> Self {
        Self {
            plan,
            controls,
            orders: plan.orders.iter().cloned().collect(),
            current: None,
            run: None,
            scratch: None,
            tasks: VecDeque::new(),
            window: VecDeque::new(),
            alone: true,
            ready: VecDeque::new(),
            stop: AtomicBool::new(false),
            progress: Progress::default(),
            topologies: BTreeMap::new(),
        }
    }

    /// Whether the abort check or the progress callback stopped the generation.
    pub fn aborted(&self) -> bool {
        self.stop.load(Ordering::Relaxed)
    }

    fn parallel(&self) -> bool {
        #[cfg(all(feature = "rayon", not(target_arch = "wasm32")))]
        return self.controls.parallel;
        #[cfg(not(all(feature = "rayon", not(target_arch = "wasm32"))))]
        false
    }

    /// Map `items` in order, on the current rayon thread pool in a parallel search.
    fn map_ordered<T: Send, R: Send>(
        &self,
        items: Vec<T>,
        f: impl Fn(T) -> R + Send + Sync,
    ) -> Vec<R> {
        #[cfg(all(feature = "rayon", not(target_arch = "wasm32")))]
        if self.controls.parallel {
            use rayon::prelude::*;
            return items.into_par_iter().map(f).collect();
        }
        items.into_iter().map(f).collect()
    }

    fn threads(&self) -> usize {
        #[cfg(all(feature = "rayon", not(target_arch = "wasm32")))]
        if self.controls.parallel {
            return rayon::current_num_threads();
        }
        1
    }

    fn report(&self) {
        report(
            self.controls.progress.as_deref(),
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
        let root = Topology::empty(&self.plan.compiled.model, &limits);
        self.current = Some((loops, limits));
        if self.parallel() {
            self.tasks = self.frontier(root).into();
        } else {
            let scratch = self.scratch.take().unwrap_or_else(Scratch::take);
            self.run = Some(Cursor::new(&self.plan.compiled.model, root, scratch));
        }
        true
    }

    /// Split the search into subtrees, breadth-first and in search order, until there are enough
    /// to balance the threads.
    fn frontier(&self, root: Topology) -> Vec<Topology> {
        let target = TASKS_PER_THREAD * self.threads();
        let run = RunContext::new(self.plan, self.controls, &self.current, &self.stop);
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
        let cx = RunContext::new(self.plan, self.controls, &self.current, &self.stop);
        // until the first diagram is out, a step buffers one diagram per subtree
        let ahead = if self.progress.diagrams.is_zero() {
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
            #[cfg(all(feature = "rayon", not(target_arch = "wasm32")))]
            {
                use rayon::prelude::*;
                self.window
                    .par_iter_mut()
                    .for_each(|slot| slot.advance(&cx, ahead));
            }
            #[cfg(not(all(feature = "rayon", not(target_arch = "wasm32"))))]
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
            if !topologies.is_zero() {
                // a subtree's topologies are counted when its diagrams are handed out
                *self.topologies.entry(self.progress.loops).or_default() += &topologies;
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

impl<F, L, M: Emit<F, L>> Drop for Traversal<'_, F, L, M> {
    /// Hand the buffers of the sequential search to the next generation on this thread.
    fn drop(&mut self) {
        if let Some(scratch) = self
            .run
            .take()
            .map(Cursor::into_scratch)
            .or(self.scratch.take())
        {
            scratch.give();
        }
    }
}

impl<F, L, M: Emit<F, L>> Iterator for Traversal<'_, F, L, M>
where
    F: Clone + Ord + Send + Sync,
    L: Clone + Ord + Send + Sync,
{
    type Item = M::Item;

    fn next(&mut self) -> Option<M::Item> {
        loop {
            if self
                .controls
                .abort_check
                .as_ref()
                .is_some_and(|check| check())
            {
                self.stop.store(true, Ordering::Relaxed);
            }
            if self.stop.load(Ordering::Relaxed) {
                return None;
            }
            if let Some(diagram) = self.ready.pop_front() {
                if M::DECORATE {
                    self.progress.diagrams += 1u8;
                }
                self.report();
                return Some(diagram);
            }
            if let Some(run) = &mut self.run {
                let cx = RunContext::new(self.plan, self.controls, &self.current, &self.stop);
                let (callback, progress, stop) = (
                    self.controls.progress.as_deref(),
                    &mut self.progress,
                    &self.stop,
                );
                let topologies = &mut self.topologies;
                let diagram = run.next::<M>(&cx, &mut || {
                    *topologies.entry(progress.loops).or_default() += 1u8;
                    progress.topologies += 1u8;
                    report(callback, progress, stop);
                });
                match diagram {
                    Some(diagram) => {
                        if M::DECORATE {
                            self.progress.diagrams += 1u8;
                        }
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

/// Accepted diagrams in deterministic search order.
pub struct Diagrams<'a, F: Clone + Ord + Send + Sync, L: Clone + Ord + Send + Sync> {
    traversal: Traversal<'a, F, L, DiagramOutput>,
}

impl<F: Clone + Ord + Send + Sync, L: Clone + Ord + Send + Sync> Diagrams<'_, F, L> {
    pub fn aborted(&self) -> bool {
        self.traversal.aborted()
    }
}

impl<F: Clone + Ord + Send + Sync, L: Clone + Ord + Send + Sync> Iterator for Diagrams<'_, F, L> {
    type Item = Diagram<F, L>;
    fn next(&mut self) -> Option<Self::Item> {
        self.traversal.next()
    }
}
