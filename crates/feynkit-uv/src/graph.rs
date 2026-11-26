//! The graph interface of the UV combinatorics, and spinney enumeration.

use std::cell::RefCell;
use std::collections::BTreeSet;
use std::error::Error;
use std::fmt::{self, Display, Formatter};

use ahash::{AHashMap, AHashSet};
use feynkit_model::Model;
use linnet::half_edge::HedgeGraph;
use linnet::half_edge::involution::{EdgeIndex, Flow, Hedge};
use linnet::half_edge::subgraph::{
    Cycle, InternalSubGraph, ModifySubSet, PairwiseSubSetOps, SuBitGraph, SubGraphLike,
    SubGraphOps, SubSetLike, SubSetOps,
};
use thiserror::Error;
use tracing::debug;

use crate::{ApproximationType, CTIdentifier, RenormalizationPrescriptionSettings, Spinney};

/// Which check rejected a spinney or selection as belonging to another graph.
///
/// ```
/// use feynkit_uv::ForeignReason;
///
/// assert_eq!(ForeignReason::SubBasis.to_string(), "its routed sub-basis differs");
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ForeignReason {
    /// The selection has another number of half-edges than the graph.
    Width,
    /// The stored connected components differ from the graph's.
    Components,
    /// The stored loop count differs from the graph's.
    LoopCount,
    /// The stored basis differs from the sub-basis the graph routes for the
    /// selection from the given parent basis.
    SubBasis,
}

impl Display for ForeignReason {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            ForeignReason::Width => "its width differs",
            ForeignReason::Components => "its connected components differ",
            ForeignReason::LoopCount => "its loop count differs",
            ForeignReason::SubBasis => "its routed sub-basis differs",
        })
    }
}

/// Failure of a graph query made by the UV combinatorics.
///
/// The source is the implementing graph's own error.
#[derive(Debug, Error)]
pub enum UvGraphError<E> {
    #[error(
        "failed to select a compatible loop-momentum route for spinney {subgraph} (scheme {scheme}, DOD {dod})"
    )]
    CompatibleBasis {
        subgraph: String,
        scheme: ApproximationType,
        dod: i32,
        #[source]
        source: E,
    },
    #[error("failed to compute the degree of divergence of subgraph {subgraph}")]
    DegreeOfDivergence {
        subgraph: String,
        #[source]
        source: E,
    },
    #[error("subgraph {subgraph} is not within the graph's UV domain")]
    OutsideDomain { subgraph: String },
    #[error("spinney or selection {subgraph} does not belong to this graph: {reason}")]
    ForeignSpinney {
        subgraph: String,
        reason: ForeignReason,
    },
    #[error("failed to resolve the particles of subgraph {subgraph}")]
    Particle {
        subgraph: String,
        #[source]
        source: E,
    },
}

/// What the UV spinney, wood and forest combinatorics needs from a Feynman graph.
///
/// Subgraphs are half-edge selections of [`Self::hedge_graph`]. Implementors
/// supply the UV domain, power counting, particle identities and loop-momentum
/// bases; the provided methods enumerate and classify the UV-divergent subgraphs
/// (spinneys) from which [`crate::Wood`] is built.
///
/// # Domain
///
/// Enumeration and cycle unions act on [`Self::uv_domain`] of the requested
/// subgraph only: the propagators that may become UV hard, without dummy
/// attachments or anything else that carries no loop momentum. Single-spinney
/// constructors such as [`Spinney::new`] reject selections outside it.
///
/// # Power counting
///
/// Both degrees are taken in four spacetime dimensions, for a subgraph `γ`
/// that is a union of cycles of the domain (or empty, which must give zero).
/// Positive is power divergent, zero logarithmic and negative convergent.
///
/// - [`Self::degree_of_divergence`] classifies: `γ` is a spinney when it is
///   nonnegative, and it selects the scheme. It is the leading power of `1/t`
///   when the independent loop momenta of `γ`, routed in a basis of `γ` whose
///   external coordinates are the momenta crossing its boundary, are rescaled
///   as `k → k/t` with those external momenta fixed. It includes the measure
///   `t^(-4L)`, the local numerator of `γ` (its edge and vertex factors, without
///   graph-wide prefactors, projectors or dummy edges) and one massless inverse
///   quadratic per internal line; masses do not change it. Implementations may
///   track hard momenta differently, and no ordering between them is
///   guaranteed. FeynKit routes every line to common loop coordinates, so hard
///   parts of distinct lines can cancel in a local factor; GammaLoop tags each
///   line's hard momentum and does not cancel them. On FeynKit a local
///   numerator has degree zero when it is zero after routing or one routed edge
///   or vertex factor expands to zero. Indices are not contracted, so a
///   numerator that vanishes only by contraction across factors is not detected
///   and can come out at any degree, including a negative one: such a spinney
///   may be dropped although its integrand is zero. GammaLoop also gives an
///   identically zero numerator degree zero,
///   while it counts a numerator that only vanishes after routing, such as
///   `q₀⁰ + q₁⁰ − p⁰` on a triangle, by its tagged hard momenta.
/// - [`Self::local_degree_of_divergence`] prunes enumeration: a union is
///   dropped when it is negative. It must not be below the classifying degree
///   on any cycle union, or a divergent union is discarded. FeynKit guarantees
///   this by construction: its pruning degree is the larger of the superficial
///   local degree and the classifying degree, and its [`Self::passes_pruning`]
///   skips the classifying degree only when the superficial degree already
///   passes. Every union FeynKit classifies as divergent therefore passes
///   pruning (a numerator detected as vanishing is kept with degree zero on a bubble and a
///   triangle alike); whether it is reached is the separate limit described
///   under [enumeration](UvGraph#enumeration). GammaLoop prunes with its stored edge and
///   vertex degrees and does not guarantee the bound: a factor in another line's
///   momentum, an identically zero numerator on a triangle, or inverse powers of
///   external momenta (`1/(1 + E)` on a bubble) can make its pruning degree lower
///   than its classifying degree, and such a union is then not enumerated.
///
/// # Enumeration
///
/// [`Self::spinneys`] seeds with every nonempty element of the cycle space of the
/// domain (every edge-disjoint union of cycles). It first forms the union of each
/// seed with each seed and keeps the unions whose pruning degree is nonnegative;
/// it then repeatedly unites retained unions with seeds, keeping again only
/// nonnegative ones, until nothing new appears. A union that passes pruning is
/// therefore found whenever it is the union of two seeds, whatever the degrees of
/// the seeds. A union that needs more seeds is found only along some growth path
/// of retained unions, and is lost exactly when every such path passes through a
/// pruned union. The seeds are generated as subsets of a cycle basis of the
/// domain, so the domain's total cycle rank must stay below 64 (the limit of the
/// powerset iterator; for example 64 disconnected tadpoles exceed it).
///
/// A bridgeless graph is the union of two even subgraphs exactly when it has a
/// nowhere-zero 4-flow (Jaeger), and two-seed covers combine componentwise, so a
/// union that passes pruning is guaranteed when every connected component has
/// such a flow. This holds for every planar component (four-colour theorem) and
/// for every component with at most five independent loops: self-loops are even
/// subgraphs on their own, degree-two paths can be suppressed, and a
/// bridgeless-preserving vertex split (Fleischner's splitting lemma) reduces the
/// rest to bridgeless loopless cubic multigraphs with at most eight vertices,
/// all of which are 3-edge-colourable (the Petersen graph, with ten vertices and
/// six loops, is the smallest exception). The first possible loss is therefore a
/// nonplanar component with at least six loops and no such flow, every growth
/// path of which is pruned: for example the scalar Petersen graph with six
/// weighted vertices. Of physical numerators only the pure-gluon and the QED
/// Petersen vacua have been checked; neither loses a union.
///
/// ```
/// use feynkit_graph::FeynmanDiagram;
/// use feynkit_model::Model;
/// use feynkit_uv::{RenormalizationPrescriptionSettings, UvGraph};
///
/// let diagram = FeynmanDiagram::from_dot(
///     Model::phi4(),
///     "digraph basketball { edge [particle=phi]; a -> b; a -> b; a -> b; a -> b; }",
/// )?;
/// let spinneys = diagram.classify_spinneys(
///     &diagram.hedge_graph().full_filter(),
///     diagram.model(),
///     &RenormalizationPrescriptionSettings::default(),
///     diagram.loop_momentum_basis(),
/// )?;
/// // Six one-loop pairs, four two-loop triples and the full graph.
/// assert_eq!(spinneys.iter().filter(|spinney| !spinney.is_empty()).count(), 11);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub trait UvGraph {
    type Edge;
    type Vertex;
    type HedgeData;
    /// Loop-momentum basis stored on each spinney.
    type Basis;
    /// Failure of a power-counting, particle or basis query.
    type Error: Error + Send + Sync + 'static;

    /// The underlying half-edge graph.
    fn hedge_graph(&self) -> &HedgeGraph<Self::Edge, Self::Vertex, Self::HedgeData>;

    /// The part of `subgraph` on which UV subgraphs are enumerated: a subset of
    /// `subgraph` in the same half-edge index space.
    ///
    /// It keeps the propagators that may become UV hard. Half-edges outside it
    /// (dummy attachments, external legs, and any physical cut or sewn initial
    /// state the implementor excludes) never enter a spinney; boundary carriers
    /// that cross a spinney's crown stay soft external momenta of its basis.
    /// Callers of the single-spinney APIs must pass selections within it.
    ///
    /// `subgraph` must be a selection of this graph (same width); not checked.
    fn uv_domain(&self, subgraph: &SuBitGraph) -> SuBitGraph;

    /// Whether `edge` is a dummy attachment that carries no momentum.
    fn is_dummy_edge(&self, edge: EdgeIndex) -> bool;

    /// PDG code of the particle on `edge` in `model`, if the edge carries one.
    ///
    /// # Errors
    ///
    /// Fails when the edge's particle cannot be resolved in `model`.
    fn edge_pdg_code(&self, edge: EdgeIndex, model: &Model) -> Result<Option<isize>, Self::Error>;

    /// Whether the particle on `edge` is massive in `model`.
    ///
    /// # Errors
    ///
    /// Fails when the edge's particle cannot be resolved in `model`.
    fn is_massive_edge(&self, edge: EdgeIndex, model: &Model) -> Result<bool, Self::Error>;

    /// Pruning degree of divergence; see [the trait documentation](UvGraph#power-counting).
    ///
    /// # Errors
    ///
    /// Fails when the graph cannot power count `subgraph`.
    fn local_degree_of_divergence(&self, subgraph: &SuBitGraph) -> Result<i32, Self::Error>;

    /// Whether `subgraph` survives pruning, that is whether its pruning degree is
    /// nonnegative. Enumeration only needs this sign; implementations may decide
    /// it without computing the full degree.
    ///
    /// # Errors
    ///
    /// Fails when the graph cannot power count `subgraph`.
    fn passes_pruning(&self, subgraph: &SuBitGraph) -> Result<bool, Self::Error> {
        Ok(self.local_degree_of_divergence(subgraph)? >= 0)
    }

    /// Classifying degree of divergence; see
    /// [the trait documentation](UvGraph#power-counting).
    ///
    /// # Errors
    ///
    /// Fails when the graph cannot power count `subgraph`.
    fn degree_of_divergence(&self, subgraph: &SuBitGraph) -> Result<i32, Self::Error>;

    /// The basis with no loop or external generators.
    fn empty_basis(&self) -> Self::Basis;

    /// A basis of `subgraph` that reuses the loop edges of `parent` where possible,
    /// with the momentum-carrying crown of `subgraph` as its external carriers.
    ///
    /// # Errors
    ///
    /// Fails when no compatible basis exists.
    fn compatible_sub_basis(
        &self,
        subgraph: &InternalSubGraph,
        parent: &Self::Basis,
    ) -> Result<Self::Basis, Self::Error>;

    /// Number of independent loops of `subgraph`.
    ///
    /// `subgraph` must be a selection of this graph (same width); not checked.
    fn n_loops<S: SubGraphLike>(&self, subgraph: &S) -> usize {
        self.hedge_graph().cyclotomatic_number(subgraph)
    }

    /// The full crown of `subgraph` without dummy edges.
    ///
    /// `subgraph` must be a selection of this graph (same width); not checked.
    fn dummy_less_full_crown<S: SubGraphLike>(&self, subgraph: &S) -> S::Base
    where
        S::Base: ModifySubSet<Hedge> + SubGraphOps,
    {
        let graph = self.hedge_graph();
        let crown = graph.full_crown(subgraph);
        let mut dummy_less = crown.clone();
        for hedge in crown.included_iter() {
            if self.is_dummy_edge(graph[&hedge]) {
                dummy_less.sub(hedge);
            }
        }
        dummy_less
    }

    /// Signed PDG codes on the full crown of `subgraph`; outgoing flow is negated.
    ///
    /// `subgraph` must be a selection of this graph (same width); not checked.
    ///
    /// # Errors
    ///
    /// Fails when a boundary particle cannot be resolved in `model`.
    fn boundary_pdg_set<S: SubGraphLike<Base = SuBitGraph>>(
        &self,
        subgraph: &S,
        model: &Model,
    ) -> Result<BTreeSet<isize>, Self::Error> {
        let graph = self.hedge_graph();
        let mut codes = BTreeSet::new();
        for hedge in graph.full_crown(subgraph).included_iter() {
            if let Some(pdg) = self.edge_pdg_code(graph[&hedge], model)? {
                codes.insert(if graph.flow(hedge) == Flow::Source {
                    -pdg
                } else {
                    pdg
                });
            }
        }
        Ok(codes)
    }

    /// PDG codes of the internal edges of `subgraph`.
    ///
    /// `subgraph` must be a selection of this graph (same width); not checked.
    ///
    /// # Errors
    ///
    /// Fails when an internal particle cannot be resolved in `model`.
    fn internal_pdg_set<S: SubGraphLike>(
        &self,
        subgraph: &S,
        model: &Model,
    ) -> Result<BTreeSet<isize>, Self::Error> {
        let mut codes = BTreeSet::new();
        for (pair, edge_id, _) in self.hedge_graph().iter_edges_of(subgraph) {
            if pair.is_paired()
                && let Some(pdg) = self.edge_pdg_code(edge_id, model)?
            {
                codes.insert(pdg);
            }
        }
        Ok(codes)
    }

    /// Whether a massive particle crosses the boundary of `subgraph`.
    ///
    /// `subgraph` must be a selection of this graph (same width); not checked.
    ///
    /// # Errors
    ///
    /// Fails when a boundary particle cannot be resolved in `model`.
    fn has_massive_boundary_external<S: SubGraphLike<Base = SuBitGraph>>(
        &self,
        subgraph: &S,
        model: &Model,
    ) -> Result<bool, Self::Error> {
        let graph = self.hedge_graph();
        for hedge in graph.full_crown(subgraph).included_iter() {
            if self.is_massive_edge(graph[&hedge], model)? {
                return Ok(true);
            }
        }
        Ok(false)
    }

    /// The particle content that selects scheme overrides for `subgraph`.
    ///
    /// `subgraph` must be a selection of this graph (same width); not checked.
    ///
    /// # Errors
    ///
    /// Fails when a particle of `subgraph` cannot be resolved in `model`.
    fn ct_identifier<S: SubGraphLike<Base = SuBitGraph>>(
        &self,
        subgraph: &S,
        model: &Model,
    ) -> Result<CTIdentifier, Self::Error> {
        Ok(CTIdentifier::new(
            self.boundary_pdg_set(subgraph, model)?,
            Some(self.internal_pdg_set(subgraph, model)?),
        ))
    }

    /// The scheme `prescription` assigns to `subgraph` with degree `dod`.
    ///
    /// `subgraph` must be a selection of this graph (same width); not checked.
    ///
    /// # Errors
    ///
    /// Fails when a particle of `subgraph` cannot be resolved in `model`.
    fn approximation_scheme<S: SubGraphLike<Base = SuBitGraph>>(
        &self,
        subgraph: &S,
        model: &Model,
        prescription: &RenormalizationPrescriptionSettings,
        dod: i32,
    ) -> Result<ApproximationType, Self::Error> {
        Ok(prescription.approximation_scheme_for(
            &self.ct_identifier(subgraph, model)?,
            dod,
            self.has_massive_boundary_external(subgraph, model)?,
        ))
    }

    /// Classify one connected cycle union.
    ///
    /// Returns `Ok(None)` when it is convergent or its scheme is
    /// [`ApproximationType::Unsubtracted`].
    ///
    /// # Errors
    ///
    /// Fails when power counting, particle lookup in `model` or basis selection
    /// fails, or when `spinney` is not within the UV domain.
    fn classify_spinney(
        &self,
        spinney: InternalSubGraph,
        model: &Model,
        prescription: &RenormalizationPrescriptionSettings,
        parent_basis: &Self::Basis,
    ) -> Result<Option<Spinney<Self::Basis>>, UvGraphError<Self::Error>> {
        Spinney::<Self::Basis>::ensure_in_domain(&spinney, self)?;
        let dod = self
            .degree_of_divergence(&spinney.filter)
            .map_err(|source| UvGraphError::DegreeOfDivergence {
                subgraph: spinney.string_label(),
                source,
            })?;
        if dod < 0 {
            return Ok(None);
        }

        let renormalization_scheme = self
            .approximation_scheme(&spinney.filter, model, prescription, dod)
            .map_err(|source| UvGraphError::Particle {
                subgraph: spinney.string_label(),
                source,
            })?;

        if renormalization_scheme != ApproximationType::Unsubtracted {
            Spinney::with_scheme(spinney, self, parent_basis, renormalization_scheme, dod)
        } else {
            Ok(None)
        }
    }

    /// All classified spinneys of `subgraph`. The empty subgraph is included when
    /// its degree of divergence is nonnegative; [`crate::Wood`] merges it with the root.
    ///
    /// Connected cycle unions are classified individually. A disconnected union
    /// is retained, with the neutral `MUV` tag, exactly when each of its
    /// connected components was retained.
    ///
    /// # Errors
    ///
    /// Fails when power counting, particle lookup in `model` or basis selection
    /// fails.
    fn classify_spinneys<S: SubGraphLike<Base = SuBitGraph>>(
        &self,
        subgraph: &S,
        model: &Model,
        prescription: &RenormalizationPrescriptionSettings,
        parent_basis: &Self::Basis,
    ) -> Result<Vec<Spinney<Self::Basis>>, UvGraphError<Self::Error>> {
        let graph = self.hedge_graph();
        let (connected, disconnected): (Vec<_>, Vec<_>) = self
            .spinneys(subgraph)?
            .into_iter()
            .partition(|spinney| graph.connected_components(spinney).len() <= 1);

        let mut classified: Vec<_> = connected
            .into_iter()
            .filter_map(|spinney| {
                self.classify_spinney(spinney, model, prescription, parent_basis)
                    .transpose()
            })
            .collect::<Result<_, _>>()?;
        let retained_factors: AHashSet<_> = classified
            .iter()
            .map(|spinney| spinney.filter().clone())
            .collect();

        for spinney in disconnected {
            let components = graph.connected_components(&spinney);
            if !components.into_iter().all(|component| {
                let factor = InternalSubGraph::cleaned_filter_optimist(component, graph);
                retained_factors.contains(&factor.filter)
            }) {
                continue;
            }

            let dod = self
                .degree_of_divergence(&spinney.filter)
                .map_err(|source| UvGraphError::DegreeOfDivergence {
                    subgraph: spinney.string_label(),
                    source,
                })?;
            if let Some(spinney) =
                Spinney::with_scheme(spinney, self, parent_basis, ApproximationType::MUV, dod)?
            {
                classified.push(spinney);
            }
        }
        Ok(classified)
    }

    /// Every union of cycles of the UV domain of `subgraph`, without power counting.
    ///
    /// `subgraph` must be a selection of this graph (same width); not checked.
    fn all_cycle_unions<S: SubGraphLike<Base = SuBitGraph>>(
        &self,
        subgraph: &S,
    ) -> AHashSet<InternalSubGraph> {
        let graph = self.hedge_graph();
        let domain = self.uv_domain(subgraph.included());
        let all_subcycles: Vec<_> =
            Cycle::all_sum_powerset_filter_map(&graph.cycle_basis_of(&domain).0, &Some)
                .map(|cycles| {
                    cycles
                        .into_iter()
                        .map(|cycle| cycle.internal_graph(graph))
                        .collect()
                })
                .expect("the cycle powerset of a finite graph must exist");

        InternalSubGraph::all_unions_iterative(&all_subcycles)
    }

    /// The cycle unions of the UV domain of `subgraph` that survive pruning by the
    /// local degree of divergence, together with the empty subgraph.
    ///
    /// # Errors
    ///
    /// Fails when `subgraph` is not a selection of this graph, or when the local
    /// degree of divergence cannot be computed.
    fn spinneys<S: SubGraphLike<Base = SuBitGraph>>(
        &self,
        subgraph: &S,
    ) -> Result<AHashSet<InternalSubGraph>, UvGraphError<Self::Error>> {
        let graph = self.hedge_graph();
        if subgraph.included().size() != graph.n_hedges() {
            return Err(UvGraphError::ForeignSpinney {
                subgraph: subgraph.string_label(),
                reason: ForeignReason::Width,
            });
        }
        let domain = self.uv_domain(subgraph.included());
        debug!(subgraph = %graph.dot(&domain), "Spinneys of subgraph");

        if domain.is_empty() {
            let mut spinneys = AHashSet::new();
            spinneys.insert(graph.empty_subgraph());
            return Ok(spinneys);
        }

        let cycles = graph.cycle_basis_of(&domain).0;

        let all_subcycles: Vec<_> = Cycle::all_sum_powerset_filter_map(&cycles, &Some)
            .map(|cycles| {
                cycles
                    .into_iter()
                    .map(|cycle| cycle.internal_graph(graph))
                    .collect()
            })
            .expect("the cycle powerset of a finite graph must exist");

        let dod_cache = RefCell::new(AHashMap::new());
        let failure = RefCell::new(None);

        let mut spinneys: AHashSet<_> = InternalSubGraph::all_ops_iterative_filter_map(
            &all_subcycles,
            &|a, b| a.union(b),
            &|union| {
                let cached = dod_cache.borrow().get(&union).copied();
                let keep = match cached {
                    Some(keep) => keep,
                    None => match self.passes_pruning(&union.filter) {
                        Ok(keep) => {
                            dod_cache.borrow_mut().insert(union.clone(), keep);
                            keep
                        }
                        Err(source) => {
                            failure.borrow_mut().get_or_insert_with(|| {
                                UvGraphError::DegreeOfDivergence {
                                    subgraph: union.string_label(),
                                    source,
                                }
                            });
                            false
                        }
                    },
                };
                keep.then_some(union)
            },
        );
        if let Some(error) = failure.into_inner() {
            return Err(error);
        }

        spinneys.insert(graph.empty_subgraph());

        Ok(spinneys)
    }
}
