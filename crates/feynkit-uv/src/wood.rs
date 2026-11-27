//! The hedge-poset wood: spinneys ordered by inclusion.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::{self, Display, Formatter};

use linnet::half_edge::algorithms::trace_unfold::{HiddenData, Independence, TraceUnfold};
use linnet::half_edge::involution::{EdgeIndex, Flow, HedgePair};
use linnet::half_edge::nodestore::{NodeStorageOps, NodeStorageVec};
use linnet::half_edge::subgraph::{
    Inclusion, InternalSubGraph, ModifySubSet, SuBitGraph, SubSetLike, SubSetOps,
};
use linnet::half_edge::{HedgeGraph, NoData, NodeIndex};

use crate::{ApproximationType, Forests, Spinney, UvGraph, UvGraphError};

/// Spinneys ordered by inclusion, rooted at the empty spinney.
///
/// Nodes are spinneys. An edge from a smaller to a larger spinney is one
/// transition of the forest formula; its data is the Taylor scope of that
/// transition. Where the larger spinney is cycle-independent of the smaller one,
/// the scope is only the reduced complement, so the transition commutes with
/// operations on disjoint subgraphs. A disconnected union whose components were
/// all classified keeps only the edges from its components, so trace unfolding
/// composes it as a product.
///
/// [`Self::unfold`] turns the wood into its [`Forests`].
///
/// ```rust
/// use feynkit_graph::FeynmanDiagram;
/// use feynkit_model::Model;
/// use feynkit_uv::{RenormalizationPrescriptionSettings, UvGraph, Wood};
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
/// let wood = Wood::from_spinneys(spinneys, &diagram, diagram.loop_momentum_basis())?;
/// assert_eq!(wood.n_spinneys(), 12);
/// assert!(wood.root_spinney().is_empty());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Clone, Debug)]
pub struct Wood<B> {
    graph: HedgeGraph<SuBitGraph, Spinney<B>>,
    root: NodeIndex,
}

/// One edge of a [`Wood`]: an inclusion `smaller ⊂ larger` and the Taylor scope
/// of the corresponding transition.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct WoodInclusion<'a> {
    /// The wood edge.
    pub edge: EdgeIndex,
    /// Node of the included spinney.
    pub smaller: NodeIndex,
    /// Node of the including spinney.
    pub larger: NodeIndex,
    /// The half-edges the transition's Taylor operation acts on.
    pub scope: &'a SuBitGraph,
}

impl<B> Independence<HiddenData<SuBitGraph, EdgeIndex>> for Wood<B> {
    fn independent(
        &self,
        a: &HiddenData<SuBitGraph, EdgeIndex>,
        b: &HiddenData<SuBitGraph, EdgeIndex>,
    ) -> bool {
        !a.order.intersects(&b.order)
    }
}

impl<B> TraceUnfold<SuBitGraph> for Wood<B> {
    type EdgeData = SuBitGraph;
    type HedgeData = NoData;
    type NodeData = Spinney<B>;
    type NodeStorage = NodeStorageVec<Spinney<B>>;

    fn graph(
        &self,
    ) -> &HedgeGraph<Self::EdgeData, Self::NodeData, Self::HedgeData, Self::NodeStorage> {
        &self.graph
    }

    fn key(&self, e: EdgeIndex) -> SuBitGraph {
        self.graph[e].clone()
    }

    /// Treats a wood node as a factorized join when its incoming sink edges are exactly the
    /// disjoint connected components of the target spinney.
    ///
    /// The returned `SuBitGraph`s are the required branch factors for the generic trace unfold:
    /// their union must be the target filter, they must be pairwise disjoint, and there must be
    /// one factor per connected component.
    fn join_factors(&self, target: NodeIndex) -> Option<BTreeSet<SuBitGraph>> {
        if self.graph[target].n_components() < 2 {
            return None;
        }

        let factors = self
            .graph
            .iter_crown(target)
            .filter(|hedge| self.graph.flow(*hedge) == Flow::Sink)
            .map(|hedge| self.graph[self.graph[&hedge]].clone())
            .collect::<BTreeSet<_>>();
        if factors.len() != self.graph[target].n_components() {
            return None;
        }
        if factors.iter().any(|factor| {
            !self
                .graph
                .iter_nodes()
                .any(|(_, _, spinney)| spinney.n_components() == 1 && spinney.filter() == factor)
        }) {
            return None;
        }

        let mut cover: Option<SuBitGraph> = None;
        for factor in &factors {
            if let Some(acc) = &mut cover {
                if acc.intersects(factor) {
                    return None;
                }
                acc.union_with(factor);
            } else {
                cover = Some(factor.clone());
            }
        }

        (cover.as_ref() == Some(self.graph[target].filter())).then_some(factors)
    }
}

impl<B> Wood<B> {
    /// Build the wood of `spinneys`, adding the empty root if it is missing.
    /// Spinneys with equal subgraphs are merged, keeping the first.
    ///
    /// Each spinney's routing is checked against `graph` and `parent_basis`: its
    /// selection must have the graph's width and lie in the UV domain, its stored
    /// components and loop count must match the graph, and its basis must equal
    /// the compatible sub-basis the graph routes for it from `parent_basis`. This
    /// recomputes one sub-basis per spinney and rejects a spinney routed for
    /// another diagram even when its selection is structurally identical.
    ///
    /// `parent_basis` should be the basis the spinneys were classified with. A
    /// spinney classified from another basis of the same diagram is rejected
    /// when its routed sub-basis differs, as it does for the full graph; one
    /// whose sub-basis coincides is accepted. The classification itself, the
    /// stored degree of divergence and scheme, is taken as given and not
    /// recomputed, since spinneys may carry an explicitly chosen scheme and
    /// degree ([`Spinney::with_scheme`]); a spinney classified with another
    /// numerator on the same topology is accepted.
    ///
    /// # Errors
    ///
    /// Fails with [`UvGraphError::ForeignSpinney`] when a spinney does not
    /// belong to `graph`, and with [`UvGraphError::CompatibleBasis`] when the
    /// graph cannot route a spinney's sub-basis from `parent_basis`.
    ///
    /// See [`Wood`] for an example.
    pub fn from_spinneys<G, I>(
        spinneys: I,
        graph: &G,
        parent_basis: &B,
    ) -> Result<Self, UvGraphError<G::Error>>
    where
        G: UvGraph<Basis = B> + ?Sized,
        I: IntoIterator<Item = Spinney<B>>,
        B: PartialEq,
    {
        let hedge_graph = graph.hedge_graph();
        let mut by_filter = BTreeMap::new();
        for spinney in spinneys {
            spinney.ensure_belongs_to(graph, parent_basis)?;
            by_filter.entry(spinney.filter().clone()).or_insert(spinney);
        }
        let empty = Spinney::empty(graph);
        by_filter.entry(empty.filter().clone()).or_insert(empty);

        let mut unions = BTreeSet::new();
        let g: HedgeGraph<_, _> = HedgeGraph::poset(by_filter.into_values());
        let mut poset = g.map(
            |_, _, v| v,
            |_, n, pair, _, e| {
                let HedgePair::Paired { source, sink } = pair else {
                    return e.map(|_| hedge_graph.empty_subgraph());
                };
                let nsource = n.node_id_ref(source);
                let nsink = n.node_id_ref(sink);

                let source_subgraph = n.get_node_data(nsource).subgraph();
                let sink_subgraph = n.get_node_data(nsink).subgraph();
                let reduced_subgraph = sink_subgraph.subtract(source_subgraph).filter;
                let hairy_source = hedge_graph.full_crown(source_subgraph);

                if hedge_graph.bridges_of(&reduced_subgraph).is_empty()
                    && !hairy_source.intersects(&reduced_subgraph)
                {
                    // if the reduced graph is bridgless,and has no node overlap with the source, then the source subgraph is cycle independent of reduced subgraph
                    // if the source is not empty then this is a disjoint union
                    if !source_subgraph.is_empty() {
                        unions.insert(nsink);
                    }
                    e.map(|_| reduced_subgraph)
                } else {
                    e.map(|_| sink_subgraph.filter.clone())
                }
            },
            |_, d| d,
        );

        let mut to_remove: SuBitGraph = poset.empty_subgraph();
        let classified_filters = poset
            .iter_nodes()
            .map(|(_, _, spinney)| spinney.filter().clone())
            .collect::<BTreeSet<_>>();

        // Not quite transitive closure. A disconnected target factorizes only when every
        // connected component is both an independently classified UV spinney and supplied by
        // an incoming edge. In that case, keep only those component edges so trace unfolding
        // builds the canonical product. If either condition fails, the target is a genuinely
        // collective UV region whose atomic Taylor term and nested forest terms are both
        // required. Label every retained transition by that full Taylor scope: the reduced
        // complement remains available from the edge's source/sink nodes, while using it as
        // the trace order would incorrectly make a collective operation commute with its
        // nested prefix and could coalesce distinct forests into a false factorized join.
        // A user-given disconnected spinney is itself a product: a target also
        // factorizes when the maximal supplied factors, each a union of the
        // target's components with no convergent component, partition it.
        for u in unions {
            let collective_scope = poset[u].subgraph().filter.clone();
            let comps: BTreeSet<_> = hedge_graph
                .connected_components(&collective_scope)
                .into_iter()
                .collect();
            let incoming = poset
                .iter_crown(u)
                .filter(|c| poset.flow(*c) == Flow::Sink)
                .collect::<Vec<_>>();
            let supplied = incoming
                .iter()
                .map(|c| &poset[poset[c]])
                .filter(|factor| comps.contains(*factor))
                .cloned()
                .collect::<BTreeSet<_>>();
            let independently_classified = comps
                .iter()
                .all(|component| classified_filters.contains(component));
            let partition = if independently_classified && supplied == comps {
                None
            } else {
                let mut candidates = Vec::new();
                // Each incoming edge's scope is the factor its source leaves.
                for c in &incoming {
                    let factor = &poset[poset[c]];
                    let own: BTreeSet<_> = hedge_graph
                        .connected_components(factor)
                        .into_iter()
                        .collect();
                    if factor.is_empty() || factor == &collective_scope || !own.is_subset(&comps) {
                        continue;
                    }
                    let mut divergent = true;
                    for component in &own {
                        let degree = graph.degree_of_divergence(component).map_err(|source| {
                            UvGraphError::DegreeOfDivergence {
                                subgraph: component.string_label(),
                                source,
                            }
                        })?;
                        divergent &= degree >= 0;
                    }
                    if divergent {
                        candidates.push((factor.clone(), own));
                    }
                }
                let maximal: Vec<_> = candidates
                    .iter()
                    .filter(|(factor, _)| {
                        !candidates
                            .iter()
                            .any(|(other, _)| other != factor && other.includes(factor))
                    })
                    .collect();
                let covered: Vec<_> = maximal.iter().flat_map(|(_, own)| own).collect();
                let distinct: BTreeSet<_> = covered.iter().collect();
                (maximal.len() > 1
                    && covered.len() == distinct.len()
                    && distinct.len() == comps.len())
                .then(|| {
                    maximal
                        .iter()
                        .map(|(factor, _)| factor.clone())
                        .collect::<BTreeSet<_>>()
                })
            };
            if independently_classified && supplied == comps {
                for c in incoming {
                    if comps.contains(&poset[poset[&c]]) {
                        continue;
                    }
                    to_remove.add(c);
                    to_remove.add(poset.inv(c));
                }
            } else if let Some(factors) = partition {
                for c in incoming {
                    if factors.contains(&poset[poset[&c]]) {
                        continue;
                    }
                    to_remove.add(c);
                    to_remove.add(poset.inv(c));
                }
            } else {
                for c in incoming {
                    let edge = poset[&c];
                    poset[edge] = collective_scope.clone();
                }
            }
        }

        poset.delete_hedges(&to_remove);
        let root = poset
            .iter_nodes()
            .find(|(_, _, s)| s.subgraph().is_empty())
            .map(|(n, _, _)| n)
            .expect("the empty spinney was inserted");

        Ok(Wood { graph: poset, root })
    }

    /// The unions of two or more pairwise vertex-disjoint spinneys of
    /// `spinneys` that are missing from it, built as classification builds
    /// them: the neutral `MUV` scheme and the union's own degree. A spinney may
    /// itself be disconnected; equal selections count once and each missing
    /// union is built once.
    ///
    /// [`Self::from_spinneys`] builds the given spinneys only, as hosts that
    /// deliberately prune unions need. Add these first to give user-given
    /// vertex-disjoint spinneys the product forests that a classified wood has.
    ///
    /// # Errors
    ///
    /// Fails when a union cannot be power counted or routed.
    ///
    /// ```
    /// use feynkit_graph::FeynmanDiagram;
    /// use feynkit_model::Model;
    /// use feynkit_uv::{RenormalizationPrescriptionSettings, UvGraph, Wood};
    ///
    /// let spectacles = FeynmanDiagram::from_dot(
    ///     Model::phi3(),
    ///     "digraph spectacles { edge [particle=phi]; a -> a; a -> b; b -> b; }",
    /// )?;
    /// let basis = spectacles.loop_momentum_basis();
    /// let mut tadpoles = spectacles.classify_spinneys(
    ///     &spectacles.hedge_graph().full_filter(),
    ///     spectacles.model(),
    ///     &RenormalizationPrescriptionSettings::default(),
    ///     basis,
    /// )?;
    /// tadpoles.retain(|spinney| spinney.n_components() == 1 && !spinney.is_empty());
    /// let literal = Wood::from_spinneys(tadpoles.clone(), &spectacles, basis)?;
    /// assert_eq!(literal.unfold().n_forests(), 3);
    /// let products = Wood::disjoint_products(&tadpoles, &spectacles, basis)?;
    /// tadpoles.extend(products);
    /// let closed = Wood::from_spinneys(tadpoles, &spectacles, basis)?;
    /// assert_eq!(closed.unfold().n_forests(), 4);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn disjoint_products<G>(
        spinneys: &[Spinney<B>],
        graph: &G,
        parent_basis: &B,
    ) -> Result<Vec<Spinney<B>>, UvGraphError<G::Error>>
    where
        G: UvGraph<Basis = B> + ?Sized,
    {
        let hedge_graph = graph.hedge_graph();
        let present: BTreeSet<&SuBitGraph> = spinneys.iter().map(Spinney::filter).collect();
        // Equal selections would only repeat the same families.
        let mut seen = BTreeSet::new();
        let distinct: Vec<&Spinney<B>> = spinneys
            .iter()
            .filter(|spinney| !spinney.is_empty() && seen.insert(spinney.filter()))
            .collect();
        let mut built: BTreeSet<SuBitGraph> = BTreeSet::new();
        // Grow families by later spinneys that share no vertex with the family.
        let mut families: Vec<(InternalSubGraph, usize, usize)> = distinct
            .iter()
            .enumerate()
            .map(|(index, spinney)| (spinney.subgraph().clone(), index, 1))
            .collect();
        let mut products = Vec::new();
        while let Some((union, last, size)) = families.pop() {
            if size > 1 && !present.contains(&union.filter) && built.insert(union.filter.clone()) {
                let dod = graph
                    .degree_of_divergence(&union.filter)
                    .map_err(|source| UvGraphError::DegreeOfDivergence {
                        subgraph: union.string_label(),
                        source,
                    })?;
                if let Some(product) = Spinney::with_scheme(
                    union.clone(),
                    graph,
                    parent_basis,
                    ApproximationType::MUV,
                    dod,
                )? {
                    products.push(product);
                }
            }
            let crown = hedge_graph.full_crown(&union);
            for (index, candidate) in distinct.iter().enumerate().skip(last + 1) {
                if !crown.intersects(candidate.filter())
                    && !union.filter.intersects(candidate.filter())
                {
                    let mut grown = union.clone();
                    grown.union_with(candidate.subgraph());
                    families.push((grown, index, size + 1));
                }
            }
        }
        Ok(products)
    }

    /// The underlying poset graph.
    pub fn graph(&self) -> &HedgeGraph<SuBitGraph, Spinney<B>> {
        &self.graph
    }

    /// The node of the empty root spinney.
    pub fn root(&self) -> NodeIndex {
        self.root
    }

    /// The empty root spinney.
    pub fn root_spinney(&self) -> &Spinney<B> {
        &self.graph[self.root]
    }

    /// Number of spinneys, including the empty root.
    pub fn n_spinneys(&self) -> usize {
        self.graph.n_nodes()
    }

    /// Every spinney with its node, including the empty root.
    pub fn spinneys(&self) -> impl Iterator<Item = (NodeIndex, &Spinney<B>)> {
        self.graph
            .iter_nodes()
            .map(|(node, _, spinney)| (node, spinney))
    }

    /// The spinney at `node`.
    pub fn spinney(&self, node: NodeIndex) -> &Spinney<B> {
        &self.graph[node]
    }

    /// Largest loop count of any spinney.
    pub fn max_loop_count(&self) -> usize {
        self.graph
            .iter_nodes()
            .map(|(_, _, spinney)| spinney.loop_count())
            .max()
            .unwrap_or(0)
    }

    /// The edges of the wood: every strict inclusion that remains a transition
    /// after disconnected unions are factorized, with its Taylor scope. This is
    /// the order relation, not only its cover relation.
    pub fn inclusions(&self) -> impl Iterator<Item = WoodInclusion<'_>> {
        self.graph
            .iter_edges()
            .filter_map(|(pair, edge, scope)| match pair {
                HedgePair::Paired { source, sink } => Some(WoodInclusion {
                    edge,
                    smaller: self.graph.node_id(source),
                    larger: self.graph.node_id(sink),
                    scope: scope.data,
                }),
                _ => None,
            })
    }

    /// The transition of wood edge `edge`, or `None` when it is not a paired
    /// edge of this wood.
    ///
    /// ```
    /// use feynkit_graph::FeynmanDiagram;
    /// use feynkit_model::Model;
    /// use feynkit_uv::{RenormalizationPrescriptionSettings, UvGraph, Wood};
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
    /// let wood = Wood::from_spinneys(spinneys, &diagram, diagram.loop_momentum_basis())?;
    /// for inclusion in wood.inclusions() {
    ///     assert_eq!(wood.inclusion(inclusion.edge), Some(inclusion));
    /// }
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn inclusion(&self, edge: EdgeIndex) -> Option<WoodInclusion<'_>> {
        if edge.0 >= self.graph.n_edges() {
            return None;
        }
        let (scope, pair) = &self.graph[&edge];
        match *pair {
            HedgePair::Paired { source, sink } => Some(WoodInclusion {
                edge,
                smaller: self.graph.node_id(source),
                larger: self.graph.node_id(sink),
                scope,
            }),
            _ => None,
        }
    }

    /// The forest elements that `transition` adds: its larger spinney, or, for
    /// the step of a factorized product, the factors its scope selects.
    ///
    /// A disconnected union `A∪B` of classified spinneys is reached as `∅ → A`
    /// then `A → A∪B`, and the second step acts on `B` only: its Zimmermann
    /// element is `B`, not the union. A factor given as a disconnected spinney
    /// itself, such as a user-given `B = B₁∪B₂`, stays one element.
    ///
    /// # Panics
    ///
    /// Panics if `transition` is not a transition of this wood.
    ///
    /// ```
    /// use feynkit_graph::FeynmanDiagram;
    /// use feynkit_model::Model;
    /// use feynkit_uv::{RenormalizationPrescriptionSettings, UvGraph, Wood};
    ///
    /// let spectacles = FeynmanDiagram::from_dot(
    ///     Model::phi3(),
    ///     "digraph spectacles { edge [particle=phi]; a -> a; a -> b; b -> b; }",
    /// )?;
    /// let spinneys = spectacles.classify_spinneys(
    ///     &spectacles.hedge_graph().full_filter(),
    ///     spectacles.model(),
    ///     &RenormalizationPrescriptionSettings::default(),
    ///     spectacles.loop_momentum_basis(),
    /// )?;
    /// let wood = Wood::from_spinneys(spinneys, &spectacles, spectacles.loop_momentum_basis())?;
    /// for transition in wood.inclusions() {
    ///     // Every element is one connected tadpole, never their union.
    ///     for element in wood.transition_elements(&transition) {
    ///         assert_eq!(element.n_components(), 1);
    ///         assert_eq!(element.filter(), transition.scope);
    ///     }
    /// }
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn transition_elements(&self, transition: &WoodInclusion<'_>) -> Vec<&Spinney<B>> {
        let larger = self.spinney(transition.larger);
        if transition.scope == larger.filter() {
            return vec![larger];
        }
        // A given disconnected factor is one element, as given.
        if let Some((_, factor)) = self
            .spinneys()
            .find(|(_, spinney)| spinney.filter() == transition.scope)
        {
            return vec![factor];
        }
        larger
            .components()
            .iter()
            .filter(|component| transition.scope.includes(*component))
            .map(|component| {
                self.spinneys()
                    .map(|(_, spinney)| spinney)
                    .find(|spinney| spinney.filter() == component)
                    .expect("the factors of a factorized product are classified spinneys")
            })
            .collect()
    }

    /// The cover relation of spinney inclusion (its Hasse diagram): the pairs
    /// `(smaller, larger)` with no spinney strictly between them.
    pub fn cover_relation(&self) -> Vec<(NodeIndex, NodeIndex)> {
        let spinneys = self.spinneys().collect::<Vec<_>>();
        let mut covers = Vec::new();
        for &(smaller, lower) in &spinneys {
            for &(larger, upper) in &spinneys {
                let between = || {
                    spinneys
                        .iter()
                        .any(|(_, middle)| lower < *middle && *middle < upper)
                };
                if lower < upper && !between() {
                    covers.push((smaller, larger));
                }
            }
        }
        covers
    }

    /// The forest nodes for the transition along wood edge `edge_id`: the larger
    /// spinney being approximated (current) and the smaller, already treated one
    /// (given), both at `order`.
    ///
    /// # Panics
    ///
    /// Panics if `edge_id` is not a paired wood edge.
    pub fn current_given_pair(
        &self,
        edge_id: EdgeIndex,
        order: usize,
    ) -> (ForestNode<'_, B>, ForestNode<'_, B>) {
        let HedgePair::Paired { source, sink } = self.graph[&edge_id].1 else {
            panic!("edge in self is not paired");
        };

        // The given node has already been computed: it is the parent of this edge.
        let given = ForestNode {
            spinney: &self.graph[self.graph.node_id(source)],
            topo_order: order,
        };

        // The current node is shared by all union edges, since they have the same sink.
        let current = ForestNode {
            spinney: &self.graph[self.graph.node_id(sink)],
            topo_order: order,
        };

        (current, given)
    }

    /// Unfold the commuting transitions into the forests of this wood.
    pub fn unfold(self) -> Forests<B> {
        Forests::new(self)
    }
}

impl<B> Display for Wood<B> {
    /// DOT rendering with the spinney subgraph labels on nodes and transition
    /// scopes on edges.
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        self.graph.dot_impl_fmt(
            f,
            &self.graph.full_filter(),
            "start=2;\n",
            &|_| None,
            &|a| Some(format!("label=\"{}\"", a.string_label())),
            &|v| Some(format!("label=\"{}\"", v.subgraph().string_label())),
        )
    }
}

/// A spinney at a position of a topological order of the forests.
#[derive(Clone, Copy, Debug)]
pub struct ForestNode<'a, B> {
    /// The spinney.
    pub spinney: &'a Spinney<B>,
    /// Position of the operation in the forest's topological order.
    pub topo_order: usize,
}

/// An owned [`ForestNode`], for spinneys built from an operation cover.
#[derive(Clone, Debug)]
pub struct OwnedForestNode<B> {
    /// The spinney.
    pub spinney: Spinney<B>,
    /// Position of the operation in the forest's topological order.
    pub topo_order: usize,
}

impl<B> ForestNode<'_, B> {
    /// The half-edges of `self` that are not in `given`.
    pub fn reduced_subgraph(&self, given: &Self) -> SuBitGraph {
        self.spinney
            .subgraph()
            .subtract(given.spinney.subgraph())
            .filter
    }
}

impl<B> OwnedForestNode<B> {
    /// The half-edges of `self` that are not in `given`.
    pub fn reduced_subgraph(&self, given: &Self) -> SuBitGraph {
        self.spinney
            .subgraph()
            .subtract(given.spinney.subgraph())
            .filter
    }
}

impl<B> Display for ForestNode<'_, B> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "subgraph={}, topo_order={}, dod={}",
            self.spinney.filter().string_label(),
            self.topo_order,
            self.spinney.dod()
        )
    }
}

impl<B> Display for OwnedForestNode<B> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "subgraph={}, topo_order={}, dod={}",
            self.spinney.filter().string_label(),
            self.topo_order,
            self.spinney.dod()
        )
    }
}
