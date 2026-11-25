//! UV-divergent subgraphs (spinneys) and their renormalization data.

use std::cmp::Ordering;
use std::fmt::{self, Display, Formatter};
use std::hash::{Hash, Hasher};

use linnet::half_edge::involution::EdgeIndex;
use linnet::half_edge::subgraph::{Inclusion, InternalSubGraph, SuBitGraph, SubSetLike};
use tracing::debug;

use crate::{ApproximationType, ForeignReason, UvGraph, UvGraphError};

/// A UV-divergent subgraph: a union of cycles with its connected components,
/// degree of divergence, renormalization scheme and a loop-momentum basis of the
/// subgraph that is compatible with the enclosing graph's basis.
///
/// Spinneys are identified by their subgraph alone: equality and hashing ignore
/// the derived data, and the partial order is subgraph inclusion.
///
/// The basis type `B` is the host graph's [`UvGraph::Basis`].
///
/// ```
/// use feynkit_graph::FeynmanDiagram;
/// use feynkit_model::Model;
/// use feynkit_uv::{Spinney, UvGraph};
/// use linnet::half_edge::subgraph::InternalSubGraph;
///
/// let sunrise = FeynmanDiagram::from_dot(
///     Model::phi4(),
///     "digraph sunrise { edge [particle=phi]; ext [style=invis]; \
///          ext -> a; a -> b; a -> b; a -> b; b -> ext; }",
/// )?;
/// let graph = sunrise.hedge_graph();
/// let loops = InternalSubGraph::cleaned_filter_optimist(sunrise.internal_subgraph(), graph);
/// let spinney = Spinney::new(loops, &sunrise, sunrise.loop_momentum_basis())?
///     .expect("the quartic sunrise is quadratically divergent");
/// assert_eq!((spinney.loop_count(), spinney.dod(), spinney.n_components()), (2, 2, 1));
///
/// let root = Spinney::empty(&sunrise);
/// assert!(root < spinney);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Clone, Debug)]
pub struct Spinney<B> {
    subgraph: InternalSubGraph,
    components: Vec<SuBitGraph>,
    loop_count: usize,
    max_comp_loop_count: usize,
    dod: i32,
    lmb: B,
    renormalization_scheme: ApproximationType,
}

impl<B> Spinney<B> {
    /// The root of every wood: the empty subgraph with an empty basis.
    ///
    /// See [`Spinney`] for an example.
    pub fn empty<G: UvGraph<Basis = B> + ?Sized>(graph: &G) -> Self {
        Self {
            subgraph: InternalSubGraph::empty(graph.hedge_graph().n_hedges()),
            components: vec![],
            loop_count: 0,
            max_comp_loop_count: 0,
            dod: 0,
            lmb: graph.empty_basis(),
            renormalization_scheme: ApproximationType::MUV,
        }
    }

    /// Classify `subgraph` with the graph's degree of divergence and the neutral
    /// `MUV` scheme.
    ///
    /// Returns `Ok(None)` for a superficially convergent connected subgraph.
    ///
    /// # Errors
    ///
    /// Fails when `subgraph` is not within the graph's [UV domain](UvGraph::uv_domain),
    /// or when the degree of divergence or a compatible basis cannot be computed.
    pub fn new<G: UvGraph<Basis = B> + ?Sized>(
        subgraph: InternalSubGraph,
        graph: &G,
        parent_basis: &B,
    ) -> Result<Option<Self>, UvGraphError<G::Error>> {
        Self::ensure_in_domain(&subgraph, graph)?;
        let dod = graph
            .degree_of_divergence(&subgraph.filter)
            .map_err(|source| UvGraphError::DegreeOfDivergence {
                subgraph: subgraph.string_label(),
                source,
            })?;
        Self::with_scheme(subgraph, graph, parent_basis, ApproximationType::MUV, dod)
    }

    /// Build a spinney with an already selected scheme and degree of divergence.
    ///
    /// Returns `Ok(None)` for a superficially convergent connected subgraph.
    /// A disconnected union is retained even when its aggregate degree is
    /// negative: divergence is a property of each connected factor.
    ///
    /// # Errors
    ///
    /// Fails when `subgraph` is not within the graph's [UV domain](UvGraph::uv_domain)
    /// or no basis of the subgraph is compatible with `parent_basis`.
    pub fn with_scheme<G: UvGraph<Basis = B> + ?Sized>(
        subgraph: InternalSubGraph,
        graph: &G,
        parent_basis: &B,
        renormalization_scheme: ApproximationType,
        dod: i32,
    ) -> Result<Option<Self>, UvGraphError<G::Error>> {
        Self::ensure_in_domain(&subgraph, graph)?;
        let hedge_graph = graph.hedge_graph();
        let components = hedge_graph.connected_components(&subgraph);
        let max_comp_loop_count = components
            .iter()
            .map(|component| hedge_graph.cyclotomatic_number(component))
            .max()
            .unwrap_or(0);

        // Divergence is a property of each connected factor of a spinney.
        // A disconnected union may have a negative aggregate DOD even though
        // every factor was selected as divergent and must still be retained.
        if dod < 0 && components.len() <= 1 {
            return Ok(None);
        }

        let lmb = graph
            .compatible_sub_basis(&subgraph, parent_basis)
            .map_err(|source| UvGraphError::CompatibleBasis {
                subgraph: subgraph.string_label(),
                scheme: renormalization_scheme,
                dod,
                source,
            })?;

        debug!(
            dod = %dod,
            renormalization_scheme = %renormalization_scheme,
            string_label = %subgraph.string_label(),
            "Constructed spinney with scheme: {:?}",
            renormalization_scheme
        );

        Ok(Some(Self {
            loop_count: hedge_graph.cyclotomatic_number(&subgraph),
            components,
            dod,
            lmb,
            subgraph,
            renormalization_scheme,
            max_comp_loop_count,
        }))
    }

    /// Check that this spinney was built for `graph` and `parent_basis`: within
    /// the UV domain, with the components and loop count the graph gives it, and
    /// with the sub-basis the graph routes for it from `parent_basis`.
    pub(crate) fn ensure_belongs_to<G: UvGraph<Basis = B> + ?Sized>(
        &self,
        graph: &G,
        parent_basis: &B,
    ) -> Result<(), UvGraphError<G::Error>>
    where
        B: PartialEq,
    {
        Self::ensure_in_domain(&self.subgraph, graph)?;
        let hedge_graph = graph.hedge_graph();
        let foreign = |reason| UvGraphError::ForeignSpinney {
            subgraph: self.subgraph.string_label(),
            reason,
        };
        if hedge_graph.connected_components(&self.subgraph) != self.components {
            return Err(foreign(ForeignReason::Components));
        }
        if hedge_graph.cyclotomatic_number(&self.subgraph) != self.loop_count {
            return Err(foreign(ForeignReason::LoopCount));
        }
        // The selection alone cannot certify the routing: a line reversed in the
        // destination flips its loop signature. Recompute the sub-basis there.
        if !self.is_empty() {
            let routed = graph
                .compatible_sub_basis(&self.subgraph, parent_basis)
                .map_err(|source| UvGraphError::CompatibleBasis {
                    subgraph: self.subgraph.string_label(),
                    scheme: self.renormalization_scheme,
                    dod: self.dod,
                    source,
                })?;
            if routed != self.lmb {
                return Err(foreign(ForeignReason::SubBasis));
            }
        }
        Ok(())
    }

    /// Reject selections of another width or with half-edges outside the graph's
    /// UV domain, such as dummy attachments: their loops would be counted without
    /// being routed.
    pub(crate) fn ensure_in_domain<G: UvGraph<Basis = B> + ?Sized>(
        subgraph: &InternalSubGraph,
        graph: &G,
    ) -> Result<(), UvGraphError<G::Error>> {
        if subgraph.filter.size() != graph.hedge_graph().n_hedges() {
            return Err(UvGraphError::ForeignSpinney {
                subgraph: subgraph.string_label(),
                reason: ForeignReason::Width,
            });
        }
        if graph.uv_domain(&subgraph.filter) == subgraph.filter {
            Ok(())
        } else {
            Err(UvGraphError::OutsideDomain {
                subgraph: subgraph.string_label(),
            })
        }
    }

    /// The subgraph.
    pub fn subgraph(&self) -> &InternalSubGraph {
        &self.subgraph
    }

    /// Classifying degree of divergence of the whole subgraph.
    pub fn dod(&self) -> i32 {
        self.dod
    }

    /// Loop-momentum basis of the subgraph, compatible with the parent basis.
    pub fn lmb(&self) -> &B {
        &self.lmb
    }

    /// The renormalization scheme. Disconnected unions use `MUV` as a neutral
    /// tag; their component factors carry the schemes that determine the
    /// approximation.
    pub fn renormalization_scheme(&self) -> ApproximationType {
        self.renormalization_scheme
    }

    /// The half-edges of the subgraph.
    pub fn filter(&self) -> &SuBitGraph {
        &self.subgraph.filter
    }

    /// Whether this is the empty root spinney.
    pub fn is_empty(&self) -> bool {
        self.subgraph.is_empty()
    }

    /// Number of independent loops of the whole subgraph.
    pub fn loop_count(&self) -> usize {
        self.loop_count
    }

    /// Largest loop count of a connected component.
    pub fn max_comp_loop_count(&self) -> usize {
        self.max_comp_loop_count
    }

    /// The connected components, as half-edge selections.
    pub fn components(&self) -> &[SuBitGraph] {
        &self.components
    }

    /// Number of connected components; more than one marks a disconnected union.
    pub fn n_components(&self) -> usize {
        self.components.len()
    }

    /// The internal edges of the subgraph, in increasing order. `graph` must be
    /// the graph the spinney was built for; it is not checked here.
    pub fn edges<G: UvGraph<Basis = B> + ?Sized>(&self, graph: &G) -> Vec<EdgeIndex> {
        let mut edges: Vec<_> = graph
            .hedge_graph()
            .iter_edges_of(&self.subgraph)
            .filter_map(|(pair, edge, _)| pair.is_paired().then_some(edge))
            .collect();
        edges.sort();
        edges
    }
}

impl<B> Display for Spinney<B> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Spinney {{ dod: {}, components: {}, max_comp_loop_count: {},subgraph: {} }}",
            self.dod,
            self.components.len(),
            self.max_comp_loop_count,
            self.subgraph.string_label(),
        )
    }
}

impl<B> Hash for Spinney<B> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.subgraph.hash(state);
    }
}

impl<B> PartialEq for Spinney<B> {
    fn eq(&self, other: &Self) -> bool {
        self.subgraph == other.subgraph
    }
}

impl<B> Eq for Spinney<B> {}

impl<B> PartialOrd for Spinney<B> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.subgraph == other.subgraph {
            Some(Ordering::Equal)
        } else if self.subgraph.includes(&other.subgraph) {
            Some(Ordering::Greater)
        } else if other.subgraph.includes(&self.subgraph) {
            Some(Ordering::Less)
        } else {
            None
        }
    }
}
