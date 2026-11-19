//! Scoped graph preparation for the UV-counterterm form of the R* operation.
//!
//! All supports and boundary half-edges refer to the original typed graph. In
//! particular, contracting a region never renumbers an insertion's tensor slots.

use std::{
    collections::{BTreeMap, BTreeSet},
    sync::{Arc, LazyLock},
    time::Instant,
};

use color_eyre::Result;
use eyre::{WrapErr, ensure, eyre};
use linnet::half_edge::{
    HedgeGraph, PowersetIterator,
    involution::{EdgeIndex, HedgePair},
    subgraph::{Inclusion, ModifySubSet, SuBitGraph, SubSetLike, SubSetOps},
};
use symbolica::atom::{Atom, AtomCore};
use vakint::Vakint;

use crate::{
    debug_tags,
    graph::{FeynmanGraph, Graph},
    model::Model,
    utils::GS,
    uv::{
        ApproximationType, RenormalizationPrescriptionSettings, Spinney, UVgenerationSettings,
        UltravioletGraph,
        approx::CutStructure,
        counterterm_cache::{CountertermCache, CountertermRequest},
        hedge_poset::{Forests, Wood},
        literal_zero_if_vanishing,
    },
};

/// How the proper UV children `Delta(S)` of an R* plan are computed.
///
/// `Massive` is the established route: the massive R operation (Bogoliubov
/// preparation, Taylor expansion to massive vacuum integrals). `Rstar` applies
/// the complete-graph R* procedure of this module to the child itself: the
/// child's support is the root of its own plan (automatic marked edge, IR
/// corrections, outer Taylor expansion in all of its boundary momenta up to its
/// local degree, single-mass closures onto massless p-integrals) and its own
/// children follow the same route recursively. Children with fewer than
/// `min_loops` loops keep the massive route. The returned polynomial has the
/// massive route's contract: open boundary slots, every boundary momentum, no
/// projector and no graph prefactor.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum UvChildrenRoute {
    Massive,
    Rstar { min_loops: usize },
}

impl std::fmt::Display for UvChildrenRoute {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Massive => write!(f, "massive"),
            Self::Rstar { min_loops } => write!(f, "rstar:{min_loops}"),
        }
    }
}

impl UvChildrenRoute {
    /// The process-wide route selected by `GAMMALOOP_UV_CHILDREN`.
    pub fn from_environment() -> Self {
        *UV_CHILDREN_ROUTE
    }

    /// `massive` (also unset or empty), `rstar` (every child) or `rstar:N`
    /// (children with at least N loops; smaller ones stay massive).
    pub fn parse(value: &str) -> Result<Self> {
        let value = value.trim();
        if value.is_empty() || value == "massive" {
            return Ok(Self::Massive);
        }
        let (name, threshold) = match value.split_once(':') {
            Some((name, threshold)) => (name.trim(), Some(threshold.trim())),
            None => (value, None),
        };
        ensure!(
            name == "rstar",
            "GAMMALOOP_UV_CHILDREN must be `massive`, `rstar` or `rstar:N`, not `{value}`"
        );
        let min_loops = match threshold {
            None => 1,
            Some(threshold) => threshold.parse::<usize>().wrap_err_with(|| {
                format!("GAMMALOOP_UV_CHILDREN loop threshold `{threshold}` is not an integer")
            })?,
        };
        Ok(Self::Rstar {
            min_loops: min_loops.max(1),
        })
    }

    fn uses_rstar(&self, loops: usize) -> bool {
        matches!(self, Self::Rstar { min_loops } if loops >= *min_loops)
    }
}

/// `GAMMALOOP_UV_CHILDREN`, read once. Unset keeps the massive route, so the
/// default evaluation is unchanged. An invalid value fails closed.
static UV_CHILDREN_ROUTE: LazyLock<UvChildrenRoute> = LazyLock::new(|| {
    let value = std::env::var("GAMMALOOP_UV_CHILDREN").unwrap_or_default();
    UvChildrenRoute::parse(&value).unwrap_or_else(|error| panic!("{error:#}"))
});

/// How an R*-route proper UV child with a closed tensor structure is
/// contracted. `Off` (the default) contracts the open child tensor.
/// `Isotropic { min_loops }` computes a logarithmically divergent child with
/// four open Lorentz slots and at least `min_loops` loops through its three
/// scalar metric channels ([`RStarPlan::isotropic_counterterm`]).
/// `Dirac { isotropic, dirac }` adds the Dirac closure
/// ([`RStarPlan::dirac_counterterm`]): a child with one open Dirac line whose
/// counterterm has a unique Dirac structure (a degree-0 vertex `c gamma^mu`,
/// a degree-1 self-energy `c slash(p)`) is computed as a closed trace, and the
/// isotropic closure as well when `isotropic` is set. Exact either way; which
/// is faster depends on the child class, so it is opt-in.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum UvChildClosure {
    #[default]
    Off,
    Isotropic {
        min_loops: usize,
    },
    Dirac {
        isotropic: Option<usize>,
        dirac: DiracClosure,
    },
}

/// The Dirac child closure of [`UvChildClosure::Dirac`]: children of at least
/// `min_loops` loops; `check` also computes the open tensor of every closed
/// child and requires it to equal the closed one (`DIRAC_CLOSURE_CHECK`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct DiracClosure {
    pub min_loops: usize,
    pub check: bool,
}

impl std::fmt::Display for DiracClosure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let name = if self.check { "dirac-check" } else { "dirac" };
        write!(f, "{name}:{}", self.min_loops)
    }
}

impl std::fmt::Display for UvChildClosure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Off => write!(f, "off"),
            Self::Isotropic { min_loops } => write!(f, "isotropic:{min_loops}"),
            Self::Dirac {
                isotropic: Some(min_loops),
                dirac,
            } => write!(f, "isotropic:{min_loops},{dirac}"),
            Self::Dirac {
                isotropic: None,
                dirac,
            } => write!(f, "{dirac}"),
        }
    }
}

impl UvChildClosure {
    /// `GAMMALOOP_UV_CHILD_CLOSURE`: unset, empty or `off`; otherwise a comma
    /// list of `isotropic`, `dirac` or `dirac-check` (each at most once, and
    /// not both Dirac modes), each optionally `:N` for children of at least
    /// N loops. Anything else is an error.
    pub fn parse(value: &str) -> Result<Self> {
        let value = value.trim();
        if value.is_empty() || value == "off" {
            return Ok(Self::Off);
        }
        let syntax = "GAMMALOOP_UV_CHILD_CLOSURE must be `off` or a comma list of \
                      `isotropic[:N]`, `dirac[:N]` or `dirac-check[:N]`";
        let mut isotropic = None;
        let mut dirac = None;
        for item in value.split(',') {
            let item = item.trim();
            let (name, threshold) = match item.split_once(':') {
                Some((name, threshold)) => (name.trim(), Some(threshold.trim())),
                None => (item, None),
            };
            let min_loops = match threshold {
                None => 1,
                Some(threshold) => threshold.parse::<usize>().wrap_err_with(|| {
                    format!(
                        "GAMMALOOP_UV_CHILD_CLOSURE loop threshold `{threshold}` is not an integer"
                    )
                })?,
            }
            .max(1);
            match name {
                "isotropic" => {
                    ensure!(
                        isotropic.is_none(),
                        "{syntax}, not `{value}` (isotropic twice)"
                    );
                    isotropic = Some(min_loops);
                }
                "dirac" | "dirac-check" => {
                    ensure!(
                        dirac.is_none(),
                        "{syntax}, not `{value}` (two Dirac closures)"
                    );
                    dirac = Some(DiracClosure {
                        min_loops,
                        check: name == "dirac-check",
                    });
                }
                _ => return Err(eyre!("{syntax}, not `{value}`")),
            }
        }
        Ok(match (isotropic, dirac) {
            (isotropic, Some(dirac)) => Self::Dirac { isotropic, dirac },
            (Some(min_loops), None) => Self::Isotropic { min_loops },
            (None, None) => unreachable!("a nonempty list names a closure"),
        })
    }

    pub fn from_environment() -> Self {
        *UV_CHILD_CLOSURE
    }

    /// Whether the isotropic closure applies to a child of `loops` loops.
    fn applies(&self, loops: usize) -> bool {
        match self {
            Self::Isotropic { min_loops }
            | Self::Dirac {
                isotropic: Some(min_loops),
                ..
            } => loops >= *min_loops,
            _ => false,
        }
    }

    /// The Dirac closure for a child of `loops` loops, when it applies.
    fn dirac(&self, loops: usize) -> Option<DiracClosure> {
        match self {
            Self::Dirac { dirac, .. } if loops >= dirac.min_loops => Some(*dirac),
            _ => None,
        }
    }
}

static UV_CHILD_CLOSURE: LazyLock<UvChildClosure> = LazyLock::new(|| {
    let value = std::env::var("GAMMALOOP_UV_CHILD_CLOSURE").unwrap_or_default();
    UvChildClosure::parse(&value).unwrap_or_else(|error| panic!("{error:#}"))
});

/// Children computed through [`RStarPlan::isotropic_counterterm`], for tests.
#[cfg(test)]
pub(crate) static ISOTROPIC_CHILDREN: std::sync::atomic::AtomicUsize =
    std::sync::atomic::AtomicUsize::new(0);

/// Children computed through [`RStarPlan::dirac_counterterm`], for tests.
#[cfg(test)]
pub(crate) static DIRAC_CHILDREN: std::sync::atomic::AtomicUsize =
    std::sync::atomic::AtomicUsize::new(0);

/// `dirac-check` outcomes, for tests: checked, identical, equal.
#[cfg(test)]
pub(crate) static DIRAC_CHECKS: [std::sync::atomic::AtomicUsize; 3] = [
    std::sync::atomic::AtomicUsize::new(0),
    std::sync::atomic::AtomicUsize::new(0),
    std::sync::atomic::AtomicUsize::new(0),
];

/// `GAMMALOOP_UV_CHILDREN_LOG=1`: one stderr line per proper UV child request
/// of an R*-route plan (route, loops, seconds, the child's edge ids); for smoke
/// tests, timing and per-child views.
fn uv_children_log() -> bool {
    static LOG: LazyLock<bool> =
        LazyLock::new(|| std::env::var_os("GAMMALOOP_UV_CHILDREN_LOG").is_some());
    *LOG
}

/// The internal (paired) edge ids of `support`, sorted. Edge ids are the
/// graph's `EdgeIndex`, which a DOT-parsed graph orders by the DOT edges'
/// explicit `id`; a child's own R* run and the decorated UV quotient keep them.
fn uv_child_edge_ids(graph: &Graph, support: &SuBitGraph) -> Vec<usize> {
    let mut ids = graph
        .underlying
        .iter_edges_of(support)
        .filter(|(pair, _, _)| matches!(pair, HedgePair::Paired { .. }))
        .map(|(_, id, _)| usize::from(id))
        .collect::<Vec<_>>();
    ids.sort_unstable();
    ids.dedup();
    ids
}

fn join_edge_ids(ids: &[usize]) -> String {
    ids.iter()
        .map(usize::to_string)
        .collect::<Vec<_>>()
        .join(",")
}

/// The log-only provenance fields of a `UV_CHILD` line besides its edges:
/// - `nested=1 parent_edges=… top_graph=…` when the plan's root is not every
///   non-dummy internal edge of `graph` (a child's own R* run, or the
///   remainder of a decorated UV quotient, whose consumed edges are dummies);
///   the ids stay in `graph`'s, i.e. the top-level diagram's, numbering;
/// - `quotient=1` when `graph` has dummy internal edges (a decorated quotient:
///   its contracted UV subgraph is a vertex, its other edges keep their ids);
/// - `surrogate=1` for a `_scalar_insertion` surrogate, whose edge ids are the
///   surrogate diagram's own (its bubble chain has no image in the original).
fn uv_child_context(graph: &Graph, root: &SuBitGraph) -> String {
    let mut all = Vec::new();
    let mut dummy = false;
    for (pair, id, data) in graph.underlying.iter_edges() {
        if matches!(pair, HedgePair::Paired { .. }) {
            if data.data.is_dummy {
                dummy = true;
            } else {
                all.push(usize::from(id));
            }
        }
    }
    all.sort_unstable();
    let root_edges = uv_child_edge_ids(graph, root);
    let mut fields = String::new();
    if root_edges != all {
        fields.push_str(&format!(
            " nested=1 parent_edges={} top_graph={}",
            join_edge_ids(&root_edges),
            graph.name
        ));
    }
    if dummy {
        fields.push_str(" quotient=1");
    }
    if graph.name.ends_with("_scalar_insertion") {
        fields.push_str(" surrogate=1");
    }
    fields
}

/// Which values of the reference mass a consumer of the proper UV children
/// needs. `Zero`: only the massless limit (the physical remainders of
/// `counterterm_from_pole_part` and `prepare_physical_uv_terms`). `Reference`:
/// the complete polynomial in the marked line's reference mass, as the
/// compensated root Taylor operation requires.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ChildMass {
    Reference,
    Zero,
}

/// A local polynomial replaces this support, with its original open indices.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct InsertionSlot {
    pub support: SuBitGraph,
    pub boundary: SuBitGraph,
}

/// An explicit, possibly decorated, counterterm argument or remaining graph.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct GraphRegion {
    /// Original support before the listed contractions.
    pub support: SuBitGraph,
    /// Propagators whose loop integrations have not been consumed.
    pub remaining: SuBitGraph,
    pub insertions: Vec<InsertionSlot>,
    pub boundary: SuBitGraph,
    pub loops: usize,
}

impl GraphRegion {
    fn new(graph: &Graph, support: &SuBitGraph, contractions: &[SuBitGraph]) -> Self {
        let mut remaining = support.clone();
        let insertions = contractions
            .iter()
            .map(|contracted| {
                remaining.subtract_with(contracted);
                InsertionSlot {
                    support: contracted.clone(),
                    boundary: graph.dummy_less_full_crown(contracted),
                }
            })
            .collect();
        let quotient = Self::quotient_topology(graph, contractions);
        Self {
            support: support.clone(),
            loops: quotient.cyclotomatic_number(&remaining),
            remaining,
            insertions,
            boundary: graph.dummy_less_full_crown(support),
        }
    }

    fn quotient_topology(graph: &Graph, contractions: &[SuBitGraph]) -> HedgeGraph<(), (), ()> {
        let mut quotient = graph.underlying.just_structure();
        for contracted in contractions {
            let nodes = quotient
                .iter_nodes_of(contracted)
                .map(|(node, _, _)| node)
                .collect::<Vec<_>>();
            if !nodes.is_empty() {
                // Identify, rather than delete, so original half-edge labels survive.
                quotient.identify_nodes(&nodes, ());
            }
        }
        quotient
    }

    fn is_one_particle_irreducible(&self, graph: &Graph) -> bool {
        // A disconnected ancestor must not appear connected merely because its
        // contracted hard component has become an isolated, invisible vertex.
        if self.loops == 0 || !graph.is_connected(&self.support) {
            return false;
        }
        let quotient = self.contracted_topology(graph);
        quotient.is_connected(&self.remaining) && quotient.bridges_of(&self.remaining).is_empty()
    }

    pub(crate) fn contracted_topology(&self, graph: &Graph) -> HedgeGraph<(), (), ()> {
        let contractions = self
            .insertions
            .iter()
            .map(|slot| slot.support.clone())
            .collect::<Vec<_>>();
        Self::quotient_topology(graph, &contractions)
    }

    pub fn contains_edge(&self, graph: &Graph, edge: EdgeIndex) -> bool {
        self.remaining.includes(&graph[&edge].1)
    }

    /// A swallowed mass kills an unintegrated vacuum loop, not a local term.
    pub(crate) fn is_scaleless(&self, graph: &Graph, marked_edge: EdgeIndex) -> bool {
        self.loops > 0 && !self.contains_edge(graph, marked_edge)
    }

    fn largest_cut_component(&self, graph: &Graph, marked_edge: EdgeIndex) -> usize {
        if self.loops == 0 {
            return 0;
        }
        let mut cut = self.remaining.clone();
        cut.sub(graph[&marked_edge].1);
        let quotient = self.contracted_topology(graph);
        let cyclic = quotient.non_bridges_of(&cut);
        quotient
            .connected_components(&cyclic)
            .iter()
            .map(|component| quotient.cyclotomatic_number(component))
            .max()
            .unwrap_or(0)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum CountertermScope {
    /// Child coefficients belong inside this prepared amplitude's own -K.
    PreparedAmplitude,
    /// Evaluate each inserted monomial; its coefficient stays outside this -K.
    InsertionMonomial,
}

/// One UV counterterm call, including its own integration and pole projection.
///
/// The subtraction sign -K belongs to the call, independently of an IR forest's
/// combinatorial sign. No projector from an enclosing graph is included here.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct CountertermCall {
    pub argument: GraphRegion,
    /// Upper bound from the original typed graph. Insertions can lower it.
    pub taylor_order: i32,
    pub scope: CountertermScope,
}

impl CountertermCall {
    fn original(graph: &Graph, support: &SuBitGraph, taylor_order: i32) -> Self {
        Self {
            argument: GraphRegion::new(graph, support, &[]),
            taylor_order,
            scope: CountertermScope::PreparedAmplitude,
        }
    }

    fn quotient(graph: &Graph, support: &SuBitGraph, contracted: &SuBitGraph) -> Self {
        Self {
            argument: GraphRegion::new(graph, support, std::slice::from_ref(contracted)),
            taylor_order: graph.local_dod(support),
            scope: CountertermScope::InsertionMonomial,
        }
    }
}

/// One spinney contribution: independent child calls inserted into a remainder.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UvPreparationTerm {
    pub remainder: GraphRegion,
    pub counterterms: Vec<CountertermCall>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IrForestTerm {
    /// -(-1)^|F|; each step also carries its own UV counterterm sign.
    pub sign: i8,
    /// Smaller ancestors act first. The last call contracts the whole quotient.
    /// Coefficients of the current insertion remain outside each next call.
    pub steps: Vec<CountertermCall>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IrCorrection {
    pub hard: SuBitGraph,
    /// Total momentum degree, including forced external derivative factors.
    pub taylor_order: i32,
    /// Full R_U preparation, including the hard component's overall counterterm.
    /// Remaining hard loops stay unintegrated until the root Taylor operation.
    pub preparation: Vec<UvPreparationTerm>,
    pub forest: Vec<IrForestTerm>,
}

/// Structural work before algebraic zeros, tensor moments and backend reduction.
/// These counts predict preparation size, not elapsed integration time.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct RStarCost {
    /// Largest edge-dependent reduction: a cut loop component or an IR UV map.
    /// Proper UV child loop counts do not distinguish the marked lines.
    pub max_reduction_loops: usize,
    pub ordinary_terms: usize,
    pub ir_terms: usize,
    /// UV operations across surviving IR summands, before cache reuse.
    pub quotient_calls: usize,
    /// Degree-weighted hard preparations; tensor ranks can increase actual work.
    pub hard_taylor_work: usize,
    /// Cutting outside triangles preserves potential propagator reductions.
    /// This does not certify that a backend triangle rule can be applied.
    pub cuts_triangle: bool,
}

impl RStarCost {
    pub fn terms(&self) -> usize {
        self.ordinary_terms + self.ir_terms
    }

    fn selection_key(&self) -> (usize, usize, usize, usize, bool) {
        (
            self.max_reduction_loops,
            self.terms(),
            self.quotient_calls,
            self.hard_taylor_work,
            self.cuts_triangle,
        )
    }

    fn estimate(
        graph: &Graph,
        support: &SuBitGraph,
        marked_edge: EdgeIndex,
        ordinary: &[UvPreparationTerm],
        corrections: &[IrCorrection],
    ) -> Self {
        let mut cost = Self {
            max_reduction_loops: ordinary
                .iter()
                .filter(|term| !term.remainder.is_scaleless(graph, marked_edge))
                .map(|term| term.remainder.largest_cut_component(graph, marked_edge))
                .max()
                .unwrap_or(0),
            ordinary_terms: ordinary
                .iter()
                .filter(|term| !term.remainder.is_scaleless(graph, marked_edge))
                .count(),
            cuts_triangle: Self::cuts_triangle(graph, support, marked_edge),
            ..Self::default()
        };
        for correction in corrections
            .iter()
            .filter(|correction| correction.hard.includes(&graph[&marked_edge].1))
        {
            let preparations = correction
                .preparation
                .iter()
                .filter(|term| !term.remainder.is_scaleless(graph, marked_edge))
                .collect::<Vec<_>>();
            let surviving = preparations.len();
            if surviving > 0 {
                for term in preparations {
                    cost.max_reduction_loops = cost
                        .max_reduction_loops
                        .max(term.remainder.largest_cut_component(graph, marked_edge));
                }
                for call in correction.forest.iter().flat_map(|term| &term.steps) {
                    cost.max_reduction_loops = cost.max_reduction_loops.max(call.argument.loops);
                }
            }
            cost.ir_terms += surviving * correction.forest.len();
            cost.quotient_calls += surviving
                * correction
                    .forest
                    .iter()
                    .map(|term| term.steps.len())
                    .sum::<usize>();
            cost.hard_taylor_work += surviving * (correction.taylor_order.max(0) as usize + 1);
        }
        cost
    }

    fn cuts_triangle(graph: &Graph, support: &SuBitGraph, edge: EdgeIndex) -> bool {
        let HedgePair::Paired { source, sink } = graph[&edge].1 else {
            return false;
        };
        let source = graph.node_id(source);
        let sink = graph.node_id(sink);
        if source == sink {
            return false;
        }
        let neighbors = |node| {
            graph
                .iter_edges_of(support)
                .filter_map(|(pair, _, _)| {
                    let HedgePair::Paired { source: a, sink: b } = pair else {
                        return None;
                    };
                    let a = graph.node_id(a);
                    let b = graph.node_id(b);
                    let other = if a == node {
                        b
                    } else if b == node {
                        a
                    } else {
                        return None;
                    };
                    (other != source && other != sink).then_some(other)
                })
                .collect::<BTreeSet<_>>()
        };
        !neighbors(source).is_disjoint(&neighbors(sink))
    }
}

/// Topology of -K I_rem T_compensated [R'_U + sum_gamma C^U B_gamma] G.
///
/// This catalogue is independent of Feynman-rule evaluation. Hard coefficients
/// remain outside every quotient call, while the root -K encloses the complete
/// assembled expression. The ordinary and hard preparations deliberately retain
/// contractions that swallow the marked line: only an actual remaining scaleless
/// loop integral may be discarded, never an already local mixed contribution.
pub struct RStarPlan {
    pub support: SuBitGraph,
    pub marked_edge: EdgeIndex,
    /// Total degree of the local polynomial, including external vertex momenta.
    /// This excludes external projectors, as does the graph numerator.
    pub taylor_order: i32,
    pub ordinary_uv: Vec<UvPreparationTerm>,
    pub ir: Vec<IrCorrection>,
    /// Undecorated UV children, ordered from fewer loops to more loops.
    /// Decorated quotient dependencies have their explicit order in `forest`.
    pub counterterms: Vec<CountertermCall>,
    catalogue: Vec<Spinney>,
    pub(super) counterterm_cache: Option<Arc<CountertermCache>>,
    /// Route of the proper UV children, `GAMMALOOP_UV_CHILDREN` by default.
    pub(super) uv_children: UvChildrenRoute,
    pub(super) uv_child_closure: UvChildClosure,
    /// `Some(crown)` when the root is a proper support of the graph (a UV
    /// child evaluated by R*): its boundary carriers are the full crown of
    /// `support` rather than the graph's external edges.
    pub(super) root_boundary: Option<SuBitGraph>,
}

impl RStarPlan {
    /// Reuse completed proper UV and decorated quotient counterterms across
    /// evaluations. Root projectors remain outside this shared cache.
    pub fn with_counterterm_cache(mut self, cache: Arc<CountertermCache>) -> Self {
        self.counterterm_cache = Some(cache);
        self
    }

    /// Select how this plan's proper UV children are computed (the default is
    /// the process-wide `GAMMALOOP_UV_CHILDREN`).
    /// Contract eligible R*-route children through their scalar channels.
    pub fn with_uv_child_closure(mut self, closure: UvChildClosure) -> Self {
        self.uv_child_closure = closure;
        self
    }

    pub fn with_uv_children(mut self, route: UvChildrenRoute) -> Self {
        self.uv_children = route;
        self
    }

    pub fn uv_children(&self) -> UvChildrenRoute {
        self.uv_children
    }

    /// An R* plan whose root is the proper support `support` of `graph`, with
    /// the automatic marked edge. The boundary of the root is the full crown
    /// of `support`; all edge and half-edge labels are those of `graph`.
    /// [`Self::counterterm`] of this plan is the local MS counterterm of the
    /// subgraph, with the contract of the massive R operation's `Delta(S)`.
    pub fn for_support(
        graph: &Graph,
        model: &Model,
        support: &SuBitGraph,
        settings: &UVgenerationSettings,
    ) -> Result<Self> {
        Self::for_support_with_edge(graph, model, support, None, settings)
    }

    /// [`Self::for_support`] with an explicit marked edge (automatic if `None`).
    pub fn for_support_with_edge(
        graph: &Graph,
        model: &Model,
        support: &SuBitGraph,
        marked_edge: Option<EdgeIndex>,
        settings: &UVgenerationSettings,
    ) -> Result<Self> {
        let mut plan =
            Self::build_for_support(graph, model, support.clone(), marked_edge, settings)?;
        plan.root_boundary = Some(graph.dummy_less_full_crown(support));
        Ok(plan)
    }

    pub fn new(
        graph: &Graph,
        model: &Model,
        marked_edge: EdgeIndex,
        settings: &UVgenerationSettings,
    ) -> Result<Self> {
        Self::build(graph, model, Some(marked_edge), settings)
    }

    /// Prefer smaller loop reductions, then fewer surviving graph terms,
    /// quotient calls and hard Taylor work, then a cut outside every triangle. Remaining ties
    /// use the original edge index. Classification is shared
    /// across candidates, and each hard component's forest is generated once.
    /// This graph-only heuristic does not predict numerator zeros or Forcer's
    /// reduction difficulty; an explicit mark remains available through `new`.
    pub fn with_automatic_edge(
        graph: &Graph,
        model: &Model,
        settings: &UVgenerationSettings,
    ) -> Result<Self> {
        Self::build(graph, model, None, settings)
    }

    pub fn estimated_cost(&self, graph: &Graph) -> RStarCost {
        if graph.compute_dod(&self.support) < 0 {
            return RStarCost::default();
        }
        RStarCost::estimate(
            graph,
            &self.support,
            self.marked_edge,
            &self.ordinary_uv,
            &self.ir,
        )
    }

    fn build(
        graph: &Graph,
        model: &Model,
        marked_edge: Option<EdgeIndex>,
        settings: &UVgenerationSettings,
    ) -> Result<Self> {
        let support = graph
            .no_dummy()
            .subtract(&graph.external_filter::<SuBitGraph>());
        Self::build_for_support(graph, model, support, marked_edge, settings)
    }

    fn build_for_support(
        graph: &Graph,
        model: &Model,
        support: SuBitGraph,
        marked_edge: Option<EdgeIndex>,
        settings: &UVgenerationSettings,
    ) -> Result<Self> {
        ensure!(
            marked_edge.is_none_or(|marked| graph
                .iter_edges_of(&support)
                .any(|(pair, edge, _)| { pair.is_paired() && edge == marked })),
            "the R* marked edge must be an internal, non-dummy propagator"
        );
        ensure!(
            graph.is_connected(&support)
                && graph.n_loops(&support) > 0
                && graph.bridges_of(&support).is_empty(),
            "an R* counterterm task requires a connected 1PI graph with loops"
        );
        let overall_divergent = graph.compute_dod(&support) >= 0;
        // R* scales the complete boundary numerator. Forced external momenta
        // count toward this Taylor degree even though they do not contribute
        // to loop UV power counting (for example the antighost derivative).
        let taylor_order = graph.local_dod(&support);

        // R* computes the MS UV map. Reuse the established spinney classification
        // and routing, without importing MUV finite terms or user exclusions.
        let ms_settings = UVgenerationSettings {
            subtract_uv: true,
            renormalization_prescription: RenormalizationPrescriptionSettings {
                log_divergent: ApproximationType::PolePart,
                massive_power_divergent: ApproximationType::PolePart,
                massless_power_divergent: ApproximationType::PolePart,
                overrides: Vec::new(),
            },
            ..settings.clone()
        };
        let mut catalogue =
            graph.classified_spinneys(&support, model, &ms_settings, &graph.loop_momentum_basis);
        catalogue.sort_by(|left, right| {
            left.filter()
                .n_included()
                .cmp(&right.filter().n_included())
                .then_with(|| left.filter().cmp(right.filter()))
        });
        ensure!(
            !overall_divergent || catalogue.iter().any(|spinney| spinney.filter() == &support),
            "the R* graph has no UV spinney with a compatible loop momentum basis"
        );

        let ordinary_uv = Self::preparation(graph, &support, &catalogue, false);
        // An overall convergent graph can still have subdivergence poles.
        // Its local UV map vanishes, so only the physical UV preparation is
        // needed for reconstructing its bare pole part.
        let mut ir = if overall_divergent {
            Self::ir_corrections(graph, &support, marked_edge, &catalogue)?
        } else {
            Vec::new()
        };
        let marked_edge = if let Some(marked_edge) = marked_edge {
            marked_edge
        } else {
            graph
                .iter_edges_of(&support)
                .filter(|(pair, _, _)| pair.is_paired())
                .map(|(_, edge, _)| {
                    let cost = if overall_divergent {
                        RStarCost::estimate(graph, &support, edge, &ordinary_uv, &ir)
                    } else {
                        RStarCost::default()
                    };
                    debug_tags!(#uv, #integrated, #rstar;
                        stage = "massive_edge_candidate",
                        marked_edge = usize::from(edge),
                        max_reduction_loops = cost.max_reduction_loops,
                        ordinary_terms = cost.ordinary_terms,
                        ir_terms = cost.ir_terms,
                        quotient_calls = cost.quotient_calls,
                        hard_taylor_work = cost.hard_taylor_work,
                        cuts_triangle = cost.cuts_triangle,
                        "Estimated R* work for an internal massive line"
                    );
                    (cost.selection_key(), edge)
                })
                .min()
                .ok_or_else(|| eyre!("the R* graph has no internal propagator"))?
                .1
        };
        ir.retain(|correction| correction.hard.includes(&graph[&marked_edge].1));
        let mut counterterms = ordinary_uv
            .iter()
            .chain(ir.iter().flat_map(|correction| &correction.preparation))
            .flat_map(|term| term.counterterms.iter().cloned())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        counterterms.sort_by(|left, right| {
            left.argument
                .loops
                .cmp(&right.argument.loops)
                .then_with(|| left.cmp(right))
        });
        Ok(Self {
            support,
            marked_edge,
            taylor_order,
            ordinary_uv,
            ir,
            counterterms,
            catalogue,
            counterterm_cache: None,
            uv_children: UvChildrenRoute::from_environment(),
            uv_child_closure: UvChildClosure::from_environment(),
            root_boundary: None,
        })
    }

    /// Reuse the existing UV forest for evaluations of the original UV map.
    /// Its all-mass Taylor kernel is not the R* root's single-line kernel.
    pub fn uv_forest(&self, graph: &Graph, settings: &UVgenerationSettings) -> Forests {
        Wood::from_spinneys(
            self.catalogue.clone(),
            graph,
            CutStructure::empty(graph),
            &settings.vakint,
        )
        .unfold()
    }

    /// Exact undecorated UV requests used by this plan's ordinary R children.
    /// The family mass and settings are identical to normal R* evaluation.
    /// Dynamic IR quotient insertions are not known at this stage.
    pub fn uv_child_requests<'a>(
        &self,
        graph: &Graph,
        model: &'a Model,
        settings: &UVgenerationSettings,
    ) -> Vec<CountertermRequest<'a>> {
        if self.counterterms.is_empty() {
            return Vec::new();
        }
        let reference_mass = Atom::var(*super::REFERENCE_MASS_SQUARED).sqrt();
        let (family, settings) = self.uv_family(graph, &reference_mass, settings);
        let requests = self.uv_family_requests(family, model, settings);
        if self.uv_children == UvChildrenRoute::Massive {
            return requests;
        }
        // R*-route children are computed on demand; only the children below
        // the loop threshold certainly use the massive route.
        self.counterterms
            .iter()
            .zip(requests)
            .filter(|(call, _)| !self.uv_children.uses_rstar(call.argument.loops))
            .map(|(_, request)| request)
            .collect()
    }

    fn uv_family(
        &self,
        graph: &Graph,
        reference_mass: &Atom,
        settings: &UVgenerationSettings,
    ) -> (Graph, UVgenerationSettings) {
        let mut family = graph.clone();
        family.underlying[self.marked_edge].particle = family[self.marked_edge]
            .particle
            .clone()
            .override_mass(Some(reference_mass.clone()));
        let settings = UVgenerationSettings {
            generate_integrated: true,
            // Diagnostic wrappers are not part of an inserted polynomial.
            add_marker: false,
            ..settings.clone()
        };
        (family, settings)
    }

    fn uv_family_requests<'a>(
        &self,
        family: Graph,
        model: &'a Model,
        settings: UVgenerationSettings,
    ) -> Vec<CountertermRequest<'a>> {
        // Match the MS catalogue fixed by `build`, independently of the
        // caller's ordinary-R prescription. Shared by lazy and eager requests.
        let settings = Arc::new(UVgenerationSettings {
            subtract_uv: true,
            renormalization_prescription: RenormalizationPrescriptionSettings {
                log_divergent: ApproximationType::PolePart,
                massive_power_divergent: ApproximationType::PolePart,
                massless_power_divergent: ApproximationType::PolePart,
                overrides: Vec::new(),
            },
            ..settings
        });
        let family = Arc::new(family);
        self.counterterms
            .iter()
            .map(|call| {
                CountertermRequest::new(
                    Arc::clone(&family),
                    call.argument.support.clone(),
                    model,
                    Arc::clone(&settings),
                )
            })
            .collect()
    }

    /// Evaluate only proper UV children with the existing massive R operation.
    /// `reference_mass` is a mass, rather than a mass squared. The returned
    /// polynomials retain its powers and their original open tensor slots.
    /// The N-loop root is excluded before forest construction, so this never
    /// computes the massive vacuum integral that the R* rearrangement replaces.
    ///
    /// Under [`UvChildrenRoute::Rstar`] eligible children are instead computed
    /// by R* on their own support (see [`Self::evaluate_uv_children_for`]).
    pub fn evaluate_uv_children(
        &self,
        graph: &Graph,
        model: &Model,
        reference_mass: &Atom,
        vakint: &Vakint,
        settings: &UVgenerationSettings,
    ) -> Result<BTreeMap<SuBitGraph, Atom>> {
        self.evaluate_uv_children_for(
            graph,
            model,
            reference_mass,
            vakint,
            settings,
            ChildMass::Reference,
        )
    }

    /// [`Self::evaluate_uv_children`] for a consumer that needs the children
    /// only at the given reference mass. The massive route ignores
    /// `child_mass` and is unchanged.
    ///
    /// The R* route computes a child `S` by R* on the massless `graph`
    /// restricted to `S` when `S` has at least the route's loop threshold and
    /// is a connected 1PI support. The reference mass of the family sits in a
    /// denominator only, so `Delta(S)` depends on it through even powers,
    /// which require `S` to contain the marked line and to have a degree of at
    /// least two. Only such children, when the consumer needs the reference
    /// mass, keep the massive route (as do children below the threshold).
    pub(crate) fn evaluate_uv_children_for(
        &self,
        graph: &Graph,
        model: &Model,
        reference_mass: &Atom,
        vakint: &Vakint,
        settings: &UVgenerationSettings,
        child_mass: ChildMass,
    ) -> Result<BTreeMap<SuBitGraph, Atom>> {
        if self.uv_children == UvChildrenRoute::Massive {
            return self.evaluate_massive_uv_children(
                graph,
                model,
                reference_mass,
                vakint,
                settings,
            );
        }
        let cache = self
            .counterterm_cache
            .clone()
            .unwrap_or_else(|| Arc::new(CountertermCache::default()));
        let marked = graph[&self.marked_edge].1;
        let (family, family_settings) = self.uv_family(graph, reference_mass, settings);
        let requests = self.uv_family_requests(family, model, family_settings);
        let mut children = BTreeMap::new();
        for (call, request) in self.counterterms.iter().zip(requests) {
            let support = &call.argument.support;
            let started = Instant::now();
            let route = if !self.uv_children.uses_rstar(call.argument.loops) {
                "massive"
            } else if child_mass == ChildMass::Reference
                && support.includes(&marked)
                && call.taylor_order >= 2
            {
                "massive_reference_mass"
            } else if !Self::admits_support_root(graph, support) {
                "massive_structure"
            } else {
                "rstar"
            };
            let polynomial = if route == "rstar" {
                Self::rstar_child_counterterm(
                    &cache,
                    graph,
                    support,
                    model,
                    vakint,
                    settings,
                    self.uv_children,
                    self.uv_child_closure,
                )?
            } else {
                request
                    .evaluate(&cache, vakint)?
                    .replace(GS.m_uv_expansion)
                    .with(GS.m_uv_vacuum)
            };
            tracing::info!(
                uv = true,
                rstar = true,
                timing = true,
                stage = "rstar_uv_child",
                graph = %graph.name,
                route,
                loops = call.argument.loops,
                degree = call.taylor_order,
                seconds = started.elapsed().as_secs_f64(),
                "Proper UV child of an R*-route plan"
            );
            if uv_children_log() {
                let seconds = started.elapsed().as_secs_f64();
                eprintln!(
                    "UV_CHILD graph={} route={route} loops={} degree={} parent_loops={} seconds={:.3} edges={}{}",
                    graph.name,
                    call.argument.loops,
                    call.taylor_order,
                    graph.n_loops(&self.support),
                    seconds,
                    join_edge_ids(&uv_child_edge_ids(graph, support)),
                    uv_child_context(graph, &self.support)
                );
            }
            children.insert(support.clone(), polynomial);
        }
        Ok(children)
    }

    /// The structural requirements of an R* root on a proper support.
    pub(super) fn admits_support_root(graph: &Graph, support: &SuBitGraph) -> bool {
        graph.is_connected(support)
            && graph.n_loops(support) > 0
            && graph.bridges_of(support).is_empty()
    }

    /// The Δ-cache settings key of an R* child: the serialised settings, the
    /// route, the child closure when on, and the root-truncation suffix when
    /// `GAMMALOOP_UV_CHILD_TRUNCATE` is on (absent otherwise, so the key is
    /// byte-identical to the key without the switch).
    pub(super) fn rstar_child_key(
        settings: String,
        route: UvChildrenRoute,
        closure: UvChildClosure,
        truncation: Option<String>,
    ) -> String {
        let mut key = match closure {
            UvChildClosure::Off => format!("{settings}#uv-children={route}"),
            closure => format!("{settings}#uv-children={route}#child-closure={closure}"),
        };
        if let Some(suffix) = truncation {
            key.push_str(&suffix);
        }
        key
    }

    /// `Delta(S)` of the support `S` of the massless `graph`, computed by the
    /// complete-graph R* procedure applied to `S` and memoised in `cache`
    /// under a key distinct from the massive route's raw entries. The result
    /// is identified (it carries no auxiliary UV mass) and is free of every
    /// reference mass.
    pub(crate) fn rstar_child_counterterm(
        cache: &Arc<CountertermCache>,
        graph: &Graph,
        support: &SuBitGraph,
        model: &Model,
        vakint: &Vakint,
        settings: &UVgenerationSettings,
        route: UvChildrenRoute,
        closure: UvChildClosure,
    ) -> Result<Atom> {
        // GAMMALOOP_UV_CHILD_TRUNCATE (when on): the child's Taylor order, in
        // the key as `#uv-child-truncate={mode}@{d}`.
        let mut truncation_degree = None;
        let truncation = super::expression::uv_child_truncation_key(|| {
            *truncation_degree.insert(graph.local_dod(support))
        });
        let key =
            Self::rstar_child_key(serde_json::to_string(settings)?, route, closure, truncation);
        // GAMMALOOP_IR_DIRAC_CLOSURE (when on): the child's IR hard parts are
        // closed, in the key as `#ir-dirac-closure={mode}`.
        let key = super::ir_dirac::child_key(key, super::ir_dirac::ir_dirac_mode()?);
        cache.counterterm_with(graph, support, model, key, || {
            let started = Instant::now();
            let plan = Self::for_support(graph, model, support, settings)?
                .with_uv_children(route)
                .with_uv_child_closure(closure)
                .with_counterterm_cache(Arc::clone(cache));
            debug_assert!(
                truncation_degree.is_none_or(|degree| degree == plan.taylor_order),
                "the truncation degree of the child key differs from the plan's Taylor order"
            );
            // A logarithmically divergent child with four open Lorentz slots
            // is a constant isotropic tensor: compute its three scalar
            // channels in its own root scope instead of the open tensor.
            // Decided on a cache miss only.
            let isotropic = if plan.taylor_order == 0 && closure.applies(graph.n_loops(support)) {
                super::expression::RStarIntegrand::isotropic_boundary(graph, support)?
            } else {
                None
            };
            // A child with one open Dirac line and a unique Dirac structure
            // (degree-0 vertex, degree-1 self-energy): a closed trace.
            let dirac = match closure.dirac(graph.n_loops(support)) {
                Some(mode) if isotropic.is_none() => {
                    super::expression::RStarIntegrand::dirac_boundary(
                        graph,
                        support,
                        plan.taylor_order,
                    )?
                    .map(|boundary| (mode, boundary))
                }
                _ => None,
            };
            let delta = super::expression::in_uv_child_scope(plan.taylor_order, || match isotropic
                .as_ref()
            {
                Some(pairings) => {
                    #[cfg(test)]
                    ISOTROPIC_CHILDREN.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    if uv_children_log() {
                        eprintln!(
                            "UV_CHILD_CLOSURE graph={} loops={} closure=isotropic edges={}",
                            graph.name,
                            graph.n_loops(support),
                            join_edge_ids(&uv_child_edge_ids(graph, support))
                        );
                    }
                    plan.isotropic_counterterm(graph, model, vakint, settings, pairings)
                }
                None => match dirac.as_ref() {
                    Some((mode, boundary)) => Self::dirac_child_counterterm(
                        &plan, graph, support, model, vakint, settings, *mode, boundary,
                    ),
                    None => plan.counterterm(graph, model, vakint, settings),
                },
            })
            .wrap_err_with(|| {
                format!(
                    "while computing the {}-loop UV child {} with R*",
                    graph.n_loops(support),
                    support.string_label()
                )
            })?;
            tracing::info!(
                uv = true,
                rstar = true,
                timing = true,
                stage = "rstar_uv_child_computed",
                graph = %graph.name,
                loops = graph.n_loops(support),
                marked_edge = usize::from(plan.marked_edge),
                taylor_order = plan.taylor_order,
                children = plan.counterterms.len(),
                seconds = started.elapsed().as_secs_f64(),
                "Computed a proper UV child with R*"
            );
            Ok(literal_zero_if_vanishing(delta))
        })
    }

    /// `Delta` of a child through [`Self::dirac_counterterm`]; under
    /// `dirac-check` also the open tensor, which must equal it
    /// (`DIRAC_CLOSURE_CHECK`, an error otherwise). The closed `Delta` is
    /// returned either way. Called in the child's root scope.
    #[allow(clippy::too_many_arguments)]
    fn dirac_child_counterterm(
        plan: &Self,
        graph: &Graph,
        support: &SuBitGraph,
        model: &Model,
        vakint: &Vakint,
        settings: &UVgenerationSettings,
        mode: DiracClosure,
        boundary: &super::expression::DiracBoundary,
    ) -> Result<Atom> {
        #[cfg(test)]
        DIRAC_CHILDREN.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let loops = graph.n_loops(support);
        let kind = match boundary.kind {
            super::expression::DiracKind::Vertex => "vertex",
            super::expression::DiracKind::SelfEnergy => "self_energy",
        };
        if uv_children_log() {
            eprintln!(
                "UV_CHILD_CLOSURE graph={} loops={loops} closure=dirac kind={kind}",
                graph.name
            );
        }
        let started = Instant::now();
        let closed = plan.dirac_counterterm(graph, model, vakint, settings, boundary)?;
        if !mode.check {
            return Ok(closed);
        }
        let closed_seconds = started.elapsed().as_secs_f64();
        let started = Instant::now();
        let open = plan.counterterm(graph, model, vakint, settings)?;
        let open_seconds = started.elapsed().as_secs_f64();
        let (identical, equal) = Self::dirac_closure_agrees(graph, support, model, &closed, &open)?;
        #[cfg(test)]
        {
            use std::sync::atomic::Ordering;
            DIRAC_CHECKS[0].fetch_add(1, Ordering::Relaxed);
            DIRAC_CHECKS[1].fetch_add(usize::from(identical), Ordering::Relaxed);
            DIRAC_CHECKS[2].fetch_add(usize::from(equal), Ordering::Relaxed);
        }
        eprintln!(
            "DIRAC_CLOSURE_CHECK graph={} support={} loops={loops} degree={} kind={kind} identical={identical} equal={equal} zero={} closed_seconds={closed_seconds:.3} open_seconds={open_seconds:.3}",
            graph.name,
            support.string_label(),
            plan.taylor_order,
            open.expand().is_zero(),
        );
        ensure!(
            equal,
            "DIRAC_CLOSURE_CHECK: the Dirac closure of the {loops}-loop {kind} child {} differs from the open tensor\nclosed: {closed}\nopen: {open}",
            support.string_label()
        );
        Ok(closed)
    }

    /// Whether the closed and open `Delta` of one child agree: `identical`
    /// as cached `Delta` (the cache's canonical form), `equal` after gamma and
    /// colour simplification of their difference in canonical labels.
    fn dirac_closure_agrees(
        graph: &Graph,
        support: &SuBitGraph,
        model: &Model,
        closed: &Atom,
        open: &Atom,
    ) -> Result<(bool, bool)> {
        use idenso::{
            CookMode, CookSettings,
            color::ColorSimplifySettings,
            dirac::GammaSimplifySettings,
            tensor::{AlgebraSettings, SymbolicTensor},
        };
        let identity =
            crate::uv::counterterm_identity::CountertermIdentity::new(graph, support, model)?;
        let identical = (identity.canonicalize_delta(closed)?
            - identity.canonicalize_delta(open)?)
        .expand()
        .is_zero();
        if identical {
            return Ok((true, true));
        }
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let difference = (identity.canonicalize_integrand(closed)?
            - identity.canonicalize_integrand(open)?)
        .expand();
        if difference.is_zero() {
            return Ok((false, true));
        }
        let tensor = SymbolicTensor::infer_indexed(cooking.try_cook(difference.as_view())?)?
            .undo_trace()?
            .undo_chain()?
            .simplify_algebra(&AlgebraSettings {
                gamma: Some(GammaSimplifySettings::default()),
                color: Some(ColorSimplifySettings::default()),
                ..Default::default()
            })?
            .to_dots()?;
        let reduced = cooking.uncook(tensor.expression().as_view()).expand();
        if reduced.is_zero() {
            return Ok((false, true));
        }
        let canonical = SymbolicTensor::infer(cooking.try_cook(reduced.as_view())?)?.canonize()?;
        let canonical = cooking.uncook(canonical.expression().as_view()).expand();
        Ok((false, canonical.together().is_zero()))
    }

    fn evaluate_massive_uv_children(
        &self,
        graph: &Graph,
        model: &Model,
        reference_mass: &Atom,
        vakint: &Vakint,
        settings: &UVgenerationSettings,
    ) -> Result<BTreeMap<SuBitGraph, Atom>> {
        let (family, settings) = self.uv_family(graph, reference_mass, settings);
        if let Some(cache) = &self.counterterm_cache {
            return self
                .counterterms
                .iter()
                .zip(self.uv_family_requests(family, model, settings))
                .map(|(call, request)| {
                    let polynomial = request
                        .evaluate(cache, vakint)?
                        .replace(GS.m_uv_expansion)
                        .with(GS.m_uv_vacuum);
                    Ok((call.argument.support.clone(), polynomial))
                })
                .collect();
        }
        let mut forest = Wood::from_spinneys(
            self.catalogue
                .iter()
                .filter(|spinney| spinney.filter() != &self.support)
                .cloned(),
            &family,
            CutStructure::empty(&family),
            &settings.vakint,
        )
        .unfold();
        forest.integrate(&family, model, vakint, &settings)?;
        self.counterterms
            .iter()
            .map(|call| {
                let support = &call.argument.support;
                let polynomial = forest
                    .counterterm_for_support(support)?
                    .replace(GS.m_uv_expansion)
                    .with(GS.m_uv_vacuum);
                Ok((support.clone(), polynomial))
            })
            .collect()
    }

    fn preparation(
        graph: &Graph,
        support: &SuBitGraph,
        catalogue: &[Spinney],
        include_overall: bool,
    ) -> Vec<UvPreparationTerm> {
        catalogue
            .iter()
            .filter(|spinney| {
                support.includes(spinney.filter())
                    && (include_overall || spinney.filter() != support)
            })
            .map(|spinney| {
                let mut components = graph.connected_components(spinney.filter());
                components.sort();
                UvPreparationTerm {
                    remainder: GraphRegion::new(graph, support, &components),
                    counterterms: components
                        .iter()
                        .map(|component| {
                            CountertermCall::original(
                                graph,
                                component,
                                graph.compute_dod(component),
                            )
                        })
                        .collect(),
                }
            })
            .collect()
    }

    /// The IR kernel is a sum over strict connected ancestor chains, separate
    /// from ordinary UV spinneys. Candidates are tested on the quotient, not
    /// restricted to the UV catalogue of the original graph.
    fn ir_corrections(
        graph: &Graph,
        support: &SuBitGraph,
        marked_edge: Option<EdgeIndex>,
        catalogue: &[Spinney],
    ) -> Result<Vec<IrCorrection>> {
        let mut corrections = Vec::new();
        for hard in catalogue.iter().filter(|spinney| {
            spinney.n_components() == 1
                && spinney.filter() != support
                && marked_edge.is_none_or(|edge| spinney.filter().includes(&graph[&edge].1))
                && graph.bridges_of(spinney.filter()).is_empty()
        }) {
            let mut forest = Vec::new();
            Self::ancestor_chains(graph, support, hard.filter(), &mut Vec::new(), &mut forest)?;
            corrections.push(IrCorrection {
                hard: hard.filter().clone(),
                taylor_order: graph.local_dod(hard.filter()),
                preparation: Self::preparation(graph, hard.filter(), catalogue, true),
                forest,
            });
        }
        Ok(corrections)
    }

    fn ancestor_chains(
        graph: &Graph,
        full: &SuBitGraph,
        contracted: &SuBitGraph,
        prefix: &mut Vec<CountertermCall>,
        terms: &mut Vec<IrForestTerm>,
    ) -> Result<()> {
        let overall = CountertermCall::quotient(graph, full, contracted);
        if overall.argument.is_one_particle_irreducible(graph) {
            let mut steps = prefix.clone();
            steps.push(overall);
            terms.push(IrForestTerm {
                sign: if prefix.len().is_multiple_of(2) {
                    -1
                } else {
                    1
                },
                steps,
            });
        }

        let available = graph
            .iter_edges_of(&full.subtract(contracted))
            .filter_map(|(pair, _, _)| pair.is_paired().then_some(pair))
            .collect::<Vec<_>>();
        ensure!(
            available.len() < usize::BITS as usize,
            "the R* ancestor powerset has too many edges to enumerate"
        );
        let edge_count = u8::try_from(available.len())
            .map_err(|_| eyre!("the R* ancestor edge count exceeds the powerset limit"))?;
        for subset in PowersetIterator::<usize>::new(edge_count).skip(1) {
            let mut ancestor = contracted.clone();
            for edge in subset.included_iter() {
                ancestor.add(available[edge]);
            }
            // Contracting G here would leave a repeated overall c-operation.
            // An insertion has total degree at most local_dod(contracted),
            // including inherited masses and forced external momenta. Thus
            // the decorated quotient's degree is at most local_dod(ancestor).
            // A negative bound makes this UV map, and its entire branch, zero.
            // Use the full local degree, not original loop power counting:
            // forced external momenta also belong to the inserted polynomial.
            if &ancestor == full || graph.local_dod(&ancestor) < 0 {
                continue;
            }
            let call = CountertermCall::quotient(graph, &ancestor, contracted);
            if !call.argument.is_one_particle_irreducible(graph) {
                continue;
            }
            prefix.push(call);
            Self::ancestor_chains(graph, full, &ancestor, prefix, terms)?;
            prefix.pop();
        }
        Ok(())
    }
}

#[cfg(test)]
pub(super) mod tests {
    use super::*;
    use crate::{
        finalized_runtime_dot, graph::parse::IntoFinalizedRuntimeGraph,
        initialisation::test_initialise, utils::load_generic_model,
    };

    fn edges(graph: &Graph, selected: &[usize]) -> SuBitGraph {
        let mut support = graph.empty_subgraph::<SuBitGraph>();
        for edge in selected {
            support.add(graph[&EdgeIndex::from(*edge)].1);
        }
        support
    }

    pub(in crate::uv::rstar) fn logarithmic_graph() -> Result<Graph> {
        let mut graph: Graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=1];
                node [num=1];
                edge [particle="scalar_1" num=1];
                A -> B [id=0 source="{ufo_order:0}" sink="{ufo_order:0}"];
                A -> B [id=1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
                A -> C [id=2 source="{ufo_order:2}" sink="{ufo_order:0}"];
                C -> B [id=3 lmb_id=1 source="{ufo_order:1}" sink="{ufo_order:2}"];
            }, "scalars"
        )?;
        let edges = graph
            .iter_edges()
            .map(|(_, edge, _)| edge)
            .collect::<Vec<_>>();
        for edge in edges {
            graph.underlying[edge].particle =
                graph[edge].particle.clone().override_mass(Some(Atom::Zero));
        }
        Ok(graph)
    }

    #[test]
    fn hard_preparation_keeps_local_mixed_term_and_original_boundary() -> Result<()> {
        test_initialise()?;
        let graph = logarithmic_graph()?;
        let model = load_generic_model("scalars");
        let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(0), &Default::default())?;
        assert_eq!(plan.taylor_order, 0);
        assert_eq!(plan.ordinary_uv.len(), 2);
        assert_eq!(plan.ir.len(), 1);
        let correction = &plan.ir[0];
        assert_eq!(correction.hard, edges(&graph, &[0, 1]));
        assert_eq!(correction.taylor_order, 0);
        assert_eq!(correction.preparation.len(), 2);
        let local = correction
            .preparation
            .iter()
            .find(|term| term.remainder.loops == 0)
            .expect("full R_U includes the local hard counterterm");
        assert!(local.remainder.remaining.is_empty());
        assert_eq!(local.counterterms.len(), 1);
        assert_eq!(local.remainder.insertions[0].support, correction.hard);
        assert_eq!(
            local.remainder.insertions[0].boundary,
            graph.dummy_less_full_crown(&correction.hard)
        );
        assert_eq!(correction.forest.len(), 1);
        assert_eq!(correction.forest[0].sign, -1);
        let steps = &correction.forest[0].steps;
        assert_eq!(steps.len(), 1);
        assert_eq!(steps[0].scope, CountertermScope::InsertionMonomial);
        assert_eq!(steps[0].argument.remaining, edges(&graph, &[2, 3]));
        assert_eq!(steps[0].argument.loops, 1);
        // Before contraction these propagators form a path, not a 1PI graph.
        assert!(!graph.bridges_of(&steps[0].argument.remaining).is_empty());
        assert!(steps[0].argument.is_one_particle_irreducible(&graph));
        assert!(
            !plan.ordinary_uv[1]
                .remainder
                .contains_edge(&graph, plan.marked_edge)
        );

        let alternate = RStarPlan::new(&graph, &model, EdgeIndex::from(2), &Default::default())?;
        assert!(alternate.ir.is_empty());
        assert!(
            alternate.ordinary_uv[1]
                .remainder
                .contains_edge(&graph, alternate.marked_edge)
        );
        assert_eq!(plan.counterterms.len(), 1);
        assert_eq!(plan.counterterms[0].argument.loops, 1);
        assert!(plan.uv_forest(&graph, &Default::default()).graph.n_nodes() > 1);
        Ok(())
    }

    #[test]
    fn automatic_edge_avoids_ir_forests_and_counts_local_mixed_terms() -> Result<()> {
        test_initialise()?;
        let graph = logarithmic_graph()?;
        let model = load_generic_model("scalars");
        let settings = UVgenerationSettings::default();
        let automatic = RStarPlan::with_automatic_edge(&graph, &model, &settings)?;
        assert_eq!(automatic.marked_edge, EdgeIndex::from(2));
        assert!(automatic.ir.is_empty());
        assert_eq!(automatic.estimated_cost(&graph).terms(), 2);
        for (_, edge, _) in graph.iter_edges_of(&automatic.support) {
            let explicit = RStarPlan::new(&graph, &model, edge, &settings)?;
            let cost = explicit.estimated_cost(&graph);
            assert!(automatic.estimated_cost(&graph).selection_key() <= cost.selection_key());
            if edge == EdgeIndex::from(0) {
                // One ordinary term is scaleless, but BOTH hard-preparation
                // terms survive: the bubble and its local counterterm.
                assert_eq!(cost.ordinary_terms, 1);
                assert_eq!(cost.ir_terms, 2);
                assert_eq!(cost.quotient_calls, 2);
            }
            if edge == automatic.marked_edge {
                assert_eq!(explicit.ordinary_uv, automatic.ordinary_uv);
                assert_eq!(explicit.ir, automatic.ir);
                assert_eq!(explicit.counterterms, automatic.counterterms);
            }
        }
        Ok(())
    }

    #[test]
    fn automatic_edge_ranks_unavoidable_ir_corrections_deterministically() -> Result<()> {
        test_initialise()?;
        let graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=1];
                node [num=1];
                edge [particle="scalar_1" num=1];
                A -> B [id=0 source="{ufo_order:0}" sink="{ufo_order:0}"];
                A -> B [id=1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
                A -> B [id=2 lmb_id=1 source="{ufo_order:2}" sink="{ufo_order:2}"];
            }, "scalars"
        )?;
        let model = load_generic_model("scalars");
        let settings = UVgenerationSettings::default();
        let automatic = RStarPlan::with_automatic_edge(&graph, &model, &settings)?;
        assert_eq!(automatic.marked_edge, EdgeIndex::from(0));
        assert!(!automatic.ir.is_empty());
        for (_, edge, _) in graph.iter_edges_of(&automatic.support) {
            let explicit = RStarPlan::new(&graph, &model, edge, &settings)?;
            assert_eq!(
                explicit.estimated_cost(&graph),
                automatic.estimated_cost(&graph)
            );
        }
        Ok(())
    }

    #[test]
    fn triangle_tie_break_ignores_orientation_parallel_edges_and_self_loops() -> Result<()> {
        test_initialise()?;
        let graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=1]; node [num=1];
                edge [particle="scalar_1" num=1];
                A -> B [id=0 source="{ufo_order:0}" sink="{ufo_order:0}"];
                B -> A [id=1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
                A -> C [id=2 source="{ufo_order:2}" sink="{ufo_order:0}"];
                C -> B [id=3 lmb_id=1 source="{ufo_order:1}" sink="{ufo_order:2}"];
                A -> A [id=4 lmb_id=2 source="{ufo_order:3}" sink="{ufo_order:4}"];
            }, "scalars"
        )?;
        let full = edges(&graph, &[0, 1, 2, 3, 4]);
        for edge in [0, 1, 2, 3].map(EdgeIndex::from) {
            assert!(RStarCost::cuts_triangle(&graph, &full, edge));
        }
        assert!(!RStarCost::cuts_triangle(&graph, &full, EdgeIndex::from(4)));
        let bubble = edges(&graph, &[0, 1, 4]);
        assert!(!RStarCost::cuts_triangle(
            &graph,
            &bubble,
            EdgeIndex::from(0)
        ));
        Ok(())
    }

    #[test]
    fn restricted_ancestors_are_coefficientwise_and_smallest_first() -> Result<()> {
        test_initialise()?;
        // Four parallel edges produce two strict, incomparable one-edge
        // ancestors of the hard bubble. Each remainder becomes a self-loop.
        let graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=1];
                node [num=1];
                edge [particle="scalar_1" num=1];
                A -> B [id=0 source="{ufo_order:0}" sink="{ufo_order:0}"];
                A -> B [id=1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
                A -> B [id=2 lmb_id=1 source="{ufo_order:2}" sink="{ufo_order:2}"];
                A -> B [id=3 lmb_id=2 source="{ufo_order:3}" sink="{ufo_order:3}"];
            }, "scalars"
        )?;
        let full = edges(&graph, &[0, 1, 2, 3]);
        let hard = edges(&graph, &[0, 1]);
        let mut terms = Vec::new();
        RStarPlan::ancestor_chains(&graph, &full, &hard, &mut Vec::new(), &mut terms)?;
        assert_eq!(terms.len(), 3);
        assert_eq!(terms[0].sign, -1);
        assert_eq!(terms[0].steps.len(), 1);
        for term in &terms[1..] {
            assert_eq!(term.sign, 1);
            assert_eq!(term.steps.len(), 2);
            let first = &term.steps[0];
            let second = &term.steps[1];
            assert_eq!(first.argument.insertions[0].support, hard);
            assert_eq!(first.argument.remaining.n_included(), 2);
            assert_eq!(
                second.argument.insertions[0].support,
                first.argument.support
            );
            assert_eq!(second.argument.support, full);
            assert_eq!(second.argument.remaining.n_included(), 2);
            assert!(term.steps.iter().all(|call| {
                call.scope == CountertermScope::InsertionMonomial && call.argument.loops > 0
            }));
        }
        Ok(())
    }

    #[test]
    fn ir_forests_skip_negative_degree_ancestors_before_expansion() -> Result<()> {
        test_initialise()?;
        let graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=1]; node [num=1];
                edge [particle="scalar_1" num=1];
                A -> B [id=0 source="{ufo_order:0}" sink="{ufo_order:0}"];
                A -> B [id=1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
                A -> C [id=2 source="{ufo_order:2}" sink="{ufo_order:0}"];
                C -> D [id=3 source="{ufo_order:1}" sink="{ufo_order:0}"];
                D -> B [id=4 lmb_id=1 source="{ufo_order:1}" sink="{ufo_order:2}"];
                C -> B [id=5 lmb_id=2 source="{ufo_order:2}" sink="{ufo_order:3}"];
            }, "scalars"
        )?;
        let hard = edges(&graph, &[0, 1]);
        let convergent = edges(&graph, &[0, 1, 2, 3, 4]);
        let divergent = edges(&graph, &[0, 1, 2, 5]);
        let full = edges(&graph, &[0, 1, 2, 3, 4, 5]);
        assert_eq!(graph.local_dod(&convergent), -2);
        assert_eq!(graph.local_dod(&divergent), 0);
        assert!(
            GraphRegion::new(&graph, &convergent, std::slice::from_ref(&hard))
                .is_one_particle_irreducible(&graph)
        );
        let mut terms = Vec::new();
        RStarPlan::ancestor_chains(&graph, &full, &hard, &mut Vec::new(), &mut terms)?;
        let calls = terms
            .iter()
            .flat_map(|term| &term.steps)
            .collect::<Vec<_>>();
        assert!(calls.iter().all(|call| call.argument.support != convergent));
        assert!(calls.iter().any(|call| call.argument.support == divergent));
        assert!(calls.iter().all(|call| call.taylor_order >= 0));
        Ok(())
    }

    #[test]
    fn uv_spinneys_reject_a_convergent_component_despite_nonnegative_total_degree() -> Result<()> {
        test_initialise()?;
        // A quark bubble with nested gluon chords. The vertex-disjoint
        // one-loop qqBB box and quark self-energy have degrees -1 and +1.
        // Linear quark numerators retain precisely this power counting.
        let graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=1];
                node [num=1];
                edge [particle="g" num=1 dod=-2];
                ext [style=invis];
                ext -> v1 [id=0 sink="{ufo_order:2}"];
                v0 -> ext [id=1 source="{ufo_order:2}"];
                v1 -> v0 [id=2 particle="d" num="Q(2,spenso::cind(1))" dod=-1 lmb_id=0 source="{ufo_order:0}" sink="{ufo_order:1}"];
                v0 -> v5 [id=3 particle="d" num="Q(3,spenso::cind(1))" dod=-1 source="{ufo_order:0}" sink="{ufo_order:1}"];
                v4 -> v1 [id=4 particle="d" num="Q(4,spenso::cind(1))" dod=-1 source="{ufo_order:0}" sink="{ufo_order:1}"];
                v2 -> v3 [id=5 particle="d" num="Q(5,spenso::cind(1))" dod=-1 lmb_id=1 source="{ufo_order:0}" sink="{ufo_order:1}"];
                v2 -> v3 [id=6 source="{ufo_order:2}" sink="{ufo_order:2}"];
                v5 -> v2 [id=7 particle="d" num="Q(7,spenso::cind(1))" dod=-1 source="{ufo_order:0}" sink="{ufo_order:1}"];
                v3 -> v4 [id=8 particle="d" num="Q(8,spenso::cind(1))" dod=-1 lmb_id=2 source="{ufo_order:0}" sink="{ufo_order:1}"];
                v4 -> v5 [id=9 source="{ufo_order:2}" sink="{ufo_order:2}"];
            }, "sm"
        )?;
        let convergent = edges(&graph, &[2, 3, 4, 9]);
        let divergent = edges(&graph, &[5, 6]);
        let mixed = convergent.union(&divergent);
        let support = edges(&graph, &[2, 3, 4, 5, 6, 7, 8, 9]);
        assert_eq!(graph.compute_dod(&convergent), -1);
        assert_eq!(graph.compute_dod(&divergent), 1);
        assert_eq!(graph.compute_dod(&mixed), 0);
        assert_eq!(graph.connected_components(&mixed).len(), 2);
        let candidate = graph
            .spinneys(&support)
            .into_iter()
            .find(|spinney| spinney.filter == mixed)
            .expect("the summed degree admits the mixed cycle union before classification");
        assert!(Spinney::new(candidate, &graph, &graph.loop_momentum_basis).is_none());

        let model = load_generic_model("sm");
        let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(7), &Default::default())?;
        assert!(
            plan.catalogue
                .iter()
                .all(|spinney| spinney.filter() != &mixed)
        );
        assert!(
            plan.counterterms
                .iter()
                .all(|call| { call.taylor_order >= 0 && call.argument.support != convergent })
        );
        assert!(
            plan.counterterms
                .iter()
                .any(|call| call.argument.support == divergent)
        );
        assert!(
            plan.counterterms
                .iter()
                .any(|call| { call.argument.support == edges(&graph, &[5, 6, 7, 8, 9]) })
        );
        Ok(())
    }

    #[test]
    fn disjoint_uv_children_keep_independent_projections() -> Result<()> {
        test_initialise()?;
        let graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=1];
                node [num=1];
                edge [particle="scalar_1" num=1];
                A -> B [id=0 source="{ufo_order:0}" sink="{ufo_order:0}"];
                A -> B [id=1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
                C -> D [id=2 source="{ufo_order:0}" sink="{ufo_order:0}"];
                C -> D [id=3 lmb_id=1 source="{ufo_order:1}" sink="{ufo_order:1}"];
                B -> C [id=4 source="{ufo_order:2}" sink="{ufo_order:2}"];
                D -> A [id=5 lmb_id=2 source="{ufo_order:2}" sink="{ufo_order:2}"];
            }, "scalars"
        )?;
        let model = load_generic_model("scalars");
        let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(0), &Default::default())?;
        let disjoint = plan
            .ordinary_uv
            .iter()
            .find(|term| term.counterterms.len() == 2)
            .expect("the two vertex-disjoint bubbles form a spinney");
        assert_eq!(disjoint.remainder.insertions.len(), 2);
        assert_eq!(disjoint.remainder.loops, 1);
        assert!(disjoint.counterterms.iter().all(|call| {
            call.argument.insertions.is_empty()
                && call.scope == CountertermScope::PreparedAmplitude
                && call.argument.loops == 1
        }));
        let disconnected = edges(&graph, &[0, 1, 2, 3]);
        assert!(plan.ir.iter().flat_map(|term| &term.forest).all(|term| {
            term.steps
                .iter()
                .all(|call| call.argument.support != disconnected)
        }));
        assert!(
            plan.counterterms
                .windows(2)
                .all(|calls| { calls[0].argument.loops <= calls[1].argument.loops })
        );
        let automatic = RStarPlan::with_automatic_edge(&graph, &model, &Default::default())?;
        assert_eq!(automatic.marked_edge, EdgeIndex::from(4));
        // Cutting a joining line separates two one-loop p-integrals. Cutting
        // inside a bubble leaves a two-loop block and a two-loop IR UV map.
        assert_eq!(automatic.estimated_cost(&graph).max_reduction_loops, 1);
        assert_eq!(plan.estimated_cost(&graph).max_reduction_loops, 2);
        Ok(())
    }

    #[test]
    fn unfolded_example_keeps_original_uv_catalogue_for_both_marked_lines() -> Result<()> {
        test_initialise()?;
        // The five-edge topology in the unfolded formula. Linear vertex
        // numerators give its original cubic-theory bounds without imposing a
        // projector or performing any tensor contraction in the graph plan.
        let graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=1];
                edge [particle="scalar_1" num=1];
                e [style=invis];
                A [num="Q(0,spenso::cind(1))" dod=1];
                B [num="Q(0,spenso::cind(2))" dod=1];
                C [num="Q(2,spenso::cind(3))+2*Q(3,spenso::cind(3))" dod=1];
                D [num="Q(2,spenso::cind(4))+3*Q(3,spenso::cind(4))" dod=1];
                A -> B [id=0 source="{ufo_order:0}" sink="{ufo_order:0}"];
                A -> C [id=1 source="{ufo_order:1}" sink="{ufo_order:0}"];
                C -> D [id=2 source="{ufo_order:1}" sink="{ufo_order:0}"];
                C -> D [id=3 lmb_id=0 source="{ufo_order:2}" sink="{ufo_order:1}"];
                D -> B [id=4 lmb_id=1 source="{ufo_order:2}" sink="{ufo_order:1}"];
                e -> A [id=5 sink="{ufo_order:2}"];
                B -> e [id=6 source="{ufo_order:2}"];
            }, "scalars"
        )?;
        let model = load_generic_model("scalars");
        let outside = RStarPlan::new(&graph, &model, EdgeIndex::from(1), &Default::default())?;
        let inside = RStarPlan::new(&graph, &model, EdgeIndex::from(2), &Default::default())?;
        assert_eq!(outside.taylor_order, 2);
        assert_eq!(outside.ordinary_uv, inside.ordinary_uv);
        assert_eq!(outside.ordinary_uv.len(), 4);
        let hard_supports = |plan: &RStarPlan| {
            plan.ir
                .iter()
                .map(|term| term.hard.clone())
                .collect::<BTreeSet<_>>()
        };
        assert_eq!(
            hard_supports(&outside),
            BTreeSet::from([edges(&graph, &[0, 1, 3, 4]), edges(&graph, &[0, 1, 2, 4])])
        );
        assert_eq!(
            hard_supports(&inside),
            BTreeSet::from([edges(&graph, &[2, 3]), edges(&graph, &[0, 1, 2, 4])])
        );
        let bubble = inside
            .ir
            .iter()
            .find(|term| term.hard == edges(&graph, &[2, 3]))
            .unwrap();
        assert_eq!(bubble.taylor_order, 2);
        assert_eq!(bubble.forest.len(), 1);
        assert_eq!(bubble.forest[0].steps.len(), 1);
        assert_eq!(bubble.preparation.len(), 2);
        Ok(())
    }

    #[test]
    fn uv_children_preserve_reference_mass_before_parent_taylor() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=0];
                node [num=1];
                edge [particle="scalar_1" num=1];
                A -> A [id=0 lmb_id=0 source="{ufo_order:0}" sink="{ufo_order:1}"];
                A -> B [id=1 source="{ufo_order:2}" sink="{ufo_order:0}"];
                A -> B [id=2 lmb_id=1 source="{ufo_order:3}" sink="{ufo_order:1}"];
            }, "scalars"
        )?;
        for edge in [0, 1, 2].map(EdgeIndex::from) {
            graph.underlying[edge].particle =
                graph[edge].particle.clone().override_mass(Some(Atom::Zero));
        }
        let model = load_generic_model("scalars");
        let settings = UVgenerationSettings::default();
        let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(0), &settings)?;
        let mass = symbolica::symbol!("rstar_plan_reference_mass"; Scalar);
        let children = plan.evaluate_uv_children(
            &graph,
            &model,
            &Atom::var(mass),
            crate::utils::vakint()?,
            &settings,
        )?;
        assert!(!children.contains_key(&plan.support));
        let tadpole = &children[&edges(&graph, &[0])];
        assert!(
            !tadpole.is_zero(),
            "the enclosing zero projector must not act"
        );
        assert!(tadpole.replace(mass).with(0).expand().is_zero());
        assert!(
            (tadpole - tadpole.replace(mass).with(1) * Atom::var(mass).pow(2))
                .expand()
                .is_zero(),
            "the complete UV polynomial retains its inherited mass squared"
        );
        Ok(())
    }
}
