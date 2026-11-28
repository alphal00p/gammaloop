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

    /// Whether every connected component has a nonnegative classifying degree
    /// in `graph`: the condition under which classification retains a
    /// spinney. [`Self::new`] also builds collective regions with a convergent
    /// factor, which hosts use as atomic Taylor sources; this tells them apart.
    ///
    /// `graph` must be the graph the spinney was built for; it is not checked.
    ///
    /// # Errors
    ///
    /// Fails when the graph cannot power count a component.
    ///
    /// ```
    /// use feynkit_graph::FeynmanDiagram;
    /// use feynkit_model::Model;
    /// use feynkit_uv::{Spinney, UvGraph};
    /// use linnet::half_edge::subgraph::{InternalSubGraph, ModifySubSet, SubSetLike};
    ///
    /// // A convergent cubic triangle and a divergent tadpole hanging off it.
    /// let kite = FeynmanDiagram::from_dot(
    ///     Model::phi3(),
    ///     "digraph kite { edge [particle=phi]; ext [style=invis]; \
    ///          a -> b; b -> c; c -> a; a -> t; t -> t; ext -> b; c -> ext; }",
    /// )?;
    /// let graph = kite.hedge_graph();
    /// let mut loops = kite.internal_subgraph();
    /// for bridge in graph.bridges_of(&loops).included_iter() {
    ///     loops.sub(bridge);
    /// }
    /// let region = InternalSubGraph::cleaned_filter_optimist(loops, graph);
    /// let collective = Spinney::new(region, &kite, kite.loop_momentum_basis())?
    ///     .expect("the union counts as logarithmic");
    /// assert_eq!(collective.n_components(), 2);
    /// assert!(!collective.has_divergent_factors(&kite)?);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn has_divergent_factors<G: UvGraph<Basis = B> + ?Sized>(
        &self,
        graph: &G,
    ) -> Result<bool, UvGraphError<G::Error>> {
        for component in &self.components {
            let degree = graph.degree_of_divergence(component).map_err(|source| {
                UvGraphError::DegreeOfDivergence {
                    subgraph: component.string_label(),
                    source,
                }
            })?;
            if degree < 0 {
                return Ok(false);
            }
        }
        Ok(true)
    }

    /// How this spinney relates to `other` in `graph`.
    ///
    /// Both spinneys must belong to `graph`; selections of another width are not
    /// compared and the result is then meaningless.
    ///
    /// ```
    /// use feynkit_graph::FeynmanDiagram;
    /// use feynkit_model::Model;
    /// use feynkit_uv::{RenormalizationPrescriptionSettings, SpinneyRelation, UvGraph};
    ///
    /// let sunrise = FeynmanDiagram::from_dot(
    ///     Model::phi4(),
    ///     "digraph sunrise { edge [particle=phi]; ext [style=invis]; \
    ///          ext -> a; a -> b; a -> b; a -> b; b -> ext; }",
    /// )?;
    /// let mut spinneys = sunrise.classify_spinneys(
    ///     &sunrise.hedge_graph().full_filter(),
    ///     sunrise.model(),
    ///     &RenormalizationPrescriptionSettings::default(),
    ///     sunrise.loop_momentum_basis(),
    /// )?;
    /// spinneys.retain(|spinney| !spinney.is_empty());
    /// spinneys.sort_by_key(|spinney| spinney.loop_count());
    /// let (bubble, other, full) = (&spinneys[0], &spinneys[1], &spinneys[3]);
    /// assert_eq!(bubble.relation(full, &sunrise), SpinneyRelation::ContainedIn);
    /// assert_eq!(full.relation(bubble, &sunrise), SpinneyRelation::Contains);
    /// // Two bubbles of the sunrise share a line.
    /// assert_eq!(bubble.relation(other, &sunrise), SpinneyRelation::Overlapping);
    /// assert!(!bubble.relation(other, &sunrise).is_compatible());
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn relation<G: UvGraph + ?Sized, C>(
        &self,
        other: &Spinney<C>,
        graph: &G,
    ) -> SpinneyRelation {
        let (this, that) = (self.filter(), other.filter());
        if this == that {
            SpinneyRelation::Equal
        } else if self.subgraph.includes(&other.subgraph) {
            SpinneyRelation::Contains
        } else if other.subgraph.includes(&self.subgraph) {
            SpinneyRelation::ContainedIn
        } else if this.intersects(that) {
            SpinneyRelation::Overlapping
        } else if graph
            .hedge_graph()
            .full_crown(&self.subgraph)
            .intersects(that)
        {
            SpinneyRelation::Touching
        } else {
            SpinneyRelation::Disjoint
        }
    }

    /// A short name listing the internal edge IDs in `graph`, such as `γ{0,1}`.
    /// `graph` must be the graph the spinney was built for; it is not checked.
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
    /// let loops = InternalSubGraph::cleaned_filter_optimist(sunrise.internal_subgraph(), sunrise.hedge_graph());
    /// let spinney = Spinney::new(loops, &sunrise, sunrise.loop_momentum_basis())?.unwrap();
    /// assert_eq!(spinney.label(&sunrise), "γ{0,1,2}");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn label<G: UvGraph<Basis = B> + ?Sized>(&self, graph: &G) -> String {
        let edges = self
            .edges(graph)
            .into_iter()
            .map(|edge| edge.0.to_string())
            .collect::<Vec<_>>();
        format!("γ{{{}}}", edges.join(","))
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

/// How two spinneys of one graph relate.
///
/// A forest of the forest formula holds only pairwise [compatible](Self::is_compatible)
/// spinneys: nested or disjoint, where disjoint spinneys share no vertex. Two
/// spinneys that share a vertex but no line touch; their union is connected and
/// they are not composed as a product. See [`Spinney::relation`] for an example.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SpinneyRelation {
    /// Both select the same half-edges.
    Equal,
    /// The first strictly includes the second.
    Contains,
    /// The first is strictly included in the second.
    ContainedIn,
    /// They share no vertex.
    Disjoint,
    /// They share a vertex but no line.
    Touching,
    /// They share a line and neither includes the other.
    Overlapping,
}

impl SpinneyRelation {
    /// Whether the two spinneys can belong to one forest: nested or disjoint.
    ///
    /// See [`Spinney::relation`] for an example.
    pub fn is_compatible(self) -> bool {
        !matches!(self, Self::Touching | Self::Overlapping)
    }
}

impl Display for SpinneyRelation {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Equal => "equal",
            Self::Contains => "contains",
            Self::ContainedIn => "contained in",
            Self::Disjoint => "disjoint",
            Self::Touching => "touching",
            Self::Overlapping => "overlapping",
        })
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
