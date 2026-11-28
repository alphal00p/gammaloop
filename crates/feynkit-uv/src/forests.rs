//! Forests: the unfolded traces of commuting UV operations of a wood.

use std::fmt::{self, Display, Formatter, Write};

use ahash::AHashMap;
use itertools::Itertools;
use linnet::half_edge::NodeIndex;
use linnet::half_edge::algorithms::DirectionBasis;
use linnet::half_edge::algorithms::topological_order::TopoError;
use linnet::half_edge::algorithms::trace_unfold::{
    HiddenData, TraceKey, TraceUnfold, UnfoldedTraceGraph,
};
use linnet::half_edge::involution::{EdgeIndex, HedgePair};
use linnet::half_edge::nodestore::NodeStorageVec;
use linnet::half_edge::subgraph::{Inclusion, InternalSubGraph, SuBitGraph, SubSetLike, SubSetOps};
use thiserror::Error;

use crate::{
    ForeignReason, ForestNode, OwnedForestNode, Spinney, UvGraph, UvGraphError, Wood, WoodInclusion,
};

/// One UV operation: a wood transition, ordered by its Taylor scope.
pub type Operation = HiddenData<SuBitGraph, EdgeIndex>;

/// A forest of the forest formula: a trace of commuting UV operations in Foata
/// normal form. Each level holds pairwise commuting operations; later levels act
/// on the result of earlier ones.
///
/// Each operation is a wood transition: its data is an edge index of the
/// [`Wood`] the forest was unfolded from, and its order is that transition's
/// Taylor scope, a half-edge selection of the source graph. Both are local to
/// that wood and graph, so operation nodes from different woods are not
/// comparable. Equality and hashing follow the trace key, which is how
/// applications key per-forest values. The empty trace is the unsubtracted root.
///
/// ```
/// use feynkit_graph::FeynmanDiagram;
/// use feynkit_model::Model;
/// use feynkit_uv::{RenormalizationPrescriptionSettings, UvGraph, Wood};
///
/// let sunrise = FeynmanDiagram::from_dot(
///     Model::phi4(),
///     "digraph sunrise { edge [particle=phi]; ext [style=invis]; \
///          ext -> a; a -> b; a -> b; a -> b; b -> ext; }",
/// )?;
/// let spinneys = sunrise.classify_spinneys(
///     &sunrise.hedge_graph().full_filter(),
///     sunrise.model(),
///     &RenormalizationPrescriptionSettings::default(),
///     sunrise.loop_momentum_basis(),
/// )?;
/// let forests = Wood::from_spinneys(spinneys, &sunrise, sunrise.loop_momentum_basis())?.unfold();
/// // A one-loop bubble followed by the full sunrise.
/// let (_, nested) = forests
///     .operations()
///     .find(|(_, forest)| forest.depth() == 2)
///     .expect("the sunrise nests each bubble in the full graph");
/// assert_eq!(nested.levels().len(), 2);
/// assert_eq!(nested.sign(), 1);
/// assert_eq!(nested.covers().as_ref(), Some(&sunrise.internal_subgraph()));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Clone, Debug, Hash, PartialEq, Eq)]
pub struct OperationNode {
    /// The trace of operations, in Foata normal form.
    pub key: TraceKey<SuBitGraph, EdgeIndex>,
}

/// Failure to locate a related node of the unfolded forests.
#[derive(Debug, Error)]
pub enum ForestError {
    #[error("{operation} has no disconnected component factors")]
    NoComponentFactors { operation: String },
    #[error("operation {op} in {operation} belongs to {owners} disconnected factors")]
    AmbiguousComponent {
        op: String,
        operation: String,
        owners: usize,
    },
    #[error("component {component} of {operation} is not in the unfolded forest")]
    MissingComponent {
        component: String,
        operation: String,
    },
    #[error("{operation} is not in the unfolded forest")]
    MissingOperation { operation: String },
    #[error("{operation} has no unique replay parent")]
    NoReplayParent { operation: String },
    #[error("cannot join integrated prefixes for {operation}")]
    IncompatiblePrefixes { operation: String },
    #[error("the forest selection has {found} half-edges, the forests have {expected}")]
    SelectionWidth { found: usize, expected: usize },
    #[error("the unfolded forests are not acyclic")]
    TopologicalOrder(#[from] TopoError),
}

/// Failure to build the spinney covered by an operation.
#[derive(Debug, Error)]
pub enum OperationSpinneyError<E> {
    #[error("failed to construct forest node for operation {operation}")]
    Spinney {
        operation: String,
        #[source]
        source: UvGraphError<E>,
    },
    #[error("operation {operation} cover has negative degree of divergence")]
    Convergent { operation: String },
}

impl OperationNode {
    /// Whether no operation of this forest touches the half-edges in `cut`.
    pub fn is_compatible_with(&self, cut: &SuBitGraph) -> bool {
        self.covers().is_none_or(|cover| !cover.intersects(cut))
    }

    /// The forest nodes of the last applied (leaf) operations, at `topo_order`.
    /// Returns `None` for the empty trace.
    ///
    /// # Panics
    ///
    /// Panics if an operation does not name a paired edge of `wood`.
    pub fn current<'a, B>(
        &'a self,
        wood: &'a Wood<B>,
        topo_order: usize,
    ) -> Option<Vec<ForestNode<'a, B>>> {
        if self.key.is_empty() {
            return None;
        }

        Some(
            self.key
                .iter_leaf_ops()
                .map(|op| {
                    let HedgePair::Paired { sink, .. } = wood.graph()[&op.data].1 else {
                        panic!("edge in trace key is not paired");
                    };
                    ForestNode {
                        spinney: wood.spinney(wood.graph().node_id(sink)),
                        topo_order,
                    }
                })
                .collect::<Vec<_>>(),
        )
    }

    /// The Foata levels, first applied first, each listing its operation scopes.
    pub fn levels(&self) -> Vec<Vec<&SuBitGraph>> {
        self.key
            .iter_levels_top_down()
            .map(|level| level.iter_leaf_ops().map(|op| &op.order).collect())
            .collect()
    }

    /// The Foata levels, first applied first, each listing the wood transitions
    /// of its operations. Each transition approximates its larger spinney after
    /// its smaller one was treated. [`Wood::transition_elements`] gives the
    /// forest elements each adds: the larger spinney, or the factor a product
    /// step acts on.
    ///
    /// # Panics
    ///
    /// Panics if an operation does not name a paired edge of `wood`, that is if
    /// this forest was unfolded from another wood.
    ///
    /// ```
    /// use feynkit_graph::FeynmanDiagram;
    /// use feynkit_model::Model;
    /// use feynkit_uv::{RenormalizationPrescriptionSettings, UvGraph, Wood};
    ///
    /// let eye = FeynmanDiagram::from_dot(
    ///     Model::phi4(),
    ///     "digraph eye { edge [particle=phi]; ext [style=invis]; ext -> a; ext -> a; \
    ///          a -> b; a -> c; b -> c; b -> c; b -> ext; c -> ext; }",
    /// )?;
    /// let spinneys = eye.classify_spinneys(
    ///     &eye.hedge_graph().full_filter(),
    ///     eye.model(),
    ///     &RenormalizationPrescriptionSettings::default(),
    ///     eye.loop_momentum_basis(),
    /// )?;
    /// let forests = Wood::from_spinneys(spinneys, &eye, eye.loop_momentum_basis())?.unfold();
    /// let (_, nested) = forests
    ///     .operations()
    ///     .find(|(_, forest)| forest.depth() == 2)
    ///     .expect("the bubble nests in the full graph");
    /// let levels = nested.transitions(forests.wood());
    /// let (inner, outer) = (levels[0][0], levels[1][0]);
    /// // The full graph is approximated after the bubble it contains.
    /// assert_eq!(outer.smaller, inner.larger);
    /// assert!(forests.wood().spinney(inner.smaller).is_empty());
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn transitions<'a, B>(&self, wood: &'a Wood<B>) -> Vec<Vec<WoodInclusion<'a>>> {
        self.key
            .iter_levels_top_down()
            .map(|level| {
                level
                    .iter_leaf_ops()
                    .map(|op| {
                        wood.inclusion(op.data)
                            .expect("forest operations name paired edges of their wood")
                    })
                    .collect()
            })
            .collect()
    }

    /// Number of operations in the forest.
    pub fn depth(&self) -> usize {
        self.key.op_count()
    }

    /// Sign of the forest term: each operation contributes one subtraction sign.
    pub fn sign(&self) -> i8 {
        if self.depth().is_multiple_of(2) {
            1
        } else {
            -1
        }
    }

    /// The Foata levels as `;`-separated levels of `,`-separated subgraph labels.
    pub fn foata_level_labels(&self) -> String {
        self.key
            .iter_levels_top_down()
            .map(|level| {
                level
                    .iter_leaf_ops()
                    .map(|op| op.order.string_label())
                    .join(",")
            })
            .join(";")
    }

    /// Union of the operation scopes; `None` for the empty trace.
    pub fn covers(&self) -> Option<SuBitGraph> {
        let mut acc: Option<SuBitGraph> = None;

        for level in self.key.iter_levels_top_down() {
            for op in level.iter_leaf_ops() {
                if let Some(a) = &mut acc {
                    a.union_with(&op.order);
                } else {
                    acc = Some(op.order.clone());
                }
            }
        }

        acc
    }

    /// The spinney covered by this forest, at `topo_order`. The empty trace
    /// covers the empty spinney.
    ///
    /// # Errors
    ///
    /// Fails when the cover cannot be classified or is convergent.
    pub fn forest_node<G: UvGraph + ?Sized>(
        &self,
        graph: &G,
        parent_basis: &G::Basis,
        topo_order: usize,
    ) -> Result<OwnedForestNode<G::Basis>, OperationSpinneyError<G::Error>> {
        let spinney = match self.covers() {
            Some(cover) if !cover.is_empty() => {
                // Cleaning indexes the cover against the graph, so reject an
                // operation of a differently sized diagram first.
                if cover.size() != graph.hedge_graph().n_hedges() {
                    return Err(OperationSpinneyError::Spinney {
                        operation: self.to_string(),
                        source: UvGraphError::ForeignSpinney {
                            subgraph: cover.string_label(),
                            reason: ForeignReason::Width,
                        },
                    });
                }
                let subgraph =
                    InternalSubGraph::cleaned_filter_optimist(cover, graph.hedge_graph());
                Spinney::new(subgraph, graph, parent_basis)
                    .map_err(|source| OperationSpinneyError::Spinney {
                        operation: self.to_string(),
                        source,
                    })?
                    .ok_or_else(|| OperationSpinneyError::Convergent {
                        operation: self.to_string(),
                    })?
            }
            _ => Spinney::empty(graph),
        };

        Ok(OwnedForestNode {
            spinney,
            topo_order,
        })
    }
}

impl AsRef<TraceKey<SuBitGraph, EdgeIndex>> for OperationNode {
    fn as_ref(&self) -> &TraceKey<SuBitGraph, EdgeIndex> {
        &self.key
    }
}

impl Display for OperationNode {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        if self.key.is_empty() {
            write!(f, "∅")
        } else {
            self.key.write_foata_like(f, |op| op.string_label())
        }
    }
}

/// One way to reach a forest from its integrated prefix: the node whose
/// integrated value is reused, followed by the wood edges applied locally.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UnionReplayState {
    /// Node of the forests whose integrated value is reused.
    pub integrated: NodeIndex,
    /// Wood edges applied locally after it, in order.
    pub local_edges: Vec<EdgeIndex>,
}

/// An edge of the unfolded forests: applying wood transition `operation` to
/// the forest at `parent` yields the forest at `child`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ForestEdge {
    /// Node of the forest before the operation.
    pub parent: NodeIndex,
    /// Node of the forest after the operation.
    pub child: NodeIndex,
    /// The wood edge of the applied transition.
    pub operation: EdgeIndex,
}

/// The unfolded forests of a [`Wood`].
///
/// Each node is one forest, an [`OperationNode`], and remembers the wood
/// spinney whose transitions produced it. An edge appends one operation.
/// Per-node computed values belong to the consumer, keyed by node or by
/// [`OperationNode`].
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
/// let forests = Wood::from_spinneys(spinneys, &diagram, diagram.loop_momentum_basis())?.unfold();
/// assert_eq!(forests.n_forests(), 46);
/// let root = forests.operation(forests.root());
/// assert_eq!((root.depth(), root.sign()), (0, 1));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub struct Forests<B> {
    graph: UnfoldedTraceGraph<Wood<B>, OperationNode, NodeStorageVec<OperationNode>>,
    root: NodeIndex,
    wood: Wood<B>,
}

impl<B> Forests<B> {
    pub(crate) fn new(wood: Wood<B>) -> Self {
        let unfolded = wood.trace_unfold::<NodeStorageVec<_>>(wood.root());
        let graph = unfolded.map(|_, _, key| OperationNode { key });
        let root = graph
            .iter_nodes()
            .find(|(_, _, operation)| operation.key.is_empty())
            .map(|(node, _, _)| node)
            .expect("no empty trace key found in unfolded hedge-poset forest");
        Self { graph, root, wood }
    }

    /// The unfolded trace graph.
    pub fn graph(
        &self,
    ) -> &UnfoldedTraceGraph<Wood<B>, OperationNode, NodeStorageVec<OperationNode>> {
        &self.graph
    }

    /// The wood these forests were unfolded from.
    pub fn wood(&self) -> &Wood<B> {
        &self.wood
    }

    /// The node of the empty (unsubtracted) forest.
    pub fn root(&self) -> NodeIndex {
        self.root
    }

    /// Number of forests, including the empty one.
    pub fn n_forests(&self) -> usize {
        self.graph.n_nodes()
    }

    /// The forest at `node`.
    pub fn operation(&self, node: NodeIndex) -> &OperationNode {
        &self.graph[node]
    }

    /// Every forest with its node.
    pub fn operations(&self) -> impl Iterator<Item = (NodeIndex, &OperationNode)> {
        self.graph
            .iter_nodes()
            .map(|(node, _, operation)| (node, operation))
    }

    /// The elements of the forest at `node` by Foata level, first applied
    /// first: pairwise nested or vertex-disjoint connected spinneys, plus any
    /// collective disconnected union the wood does not factorize. A product of
    /// disjoint spinneys lists each factor. See [`Wood::transition_elements`].
    ///
    /// # Panics
    ///
    /// Panics if `node` is not a forest of these forests.
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
    /// let forests = Wood::from_spinneys(spinneys, &spectacles, spectacles.loop_momentum_basis())?.unfold();
    /// let (product, _) = forests
    ///     .operations()
    ///     .find(|(node, _)| forests.is_disjoint_union(*node))
    ///     .expect("the two tadpoles factorize");
    /// // One commuting level holding both tadpoles.
    /// let levels = forests.forest_levels(product);
    /// assert_eq!(levels.len(), 1);
    /// let mut loops: Vec<_> = levels[0].iter().map(|s| s.edges(&spectacles)).collect();
    /// loops.sort();
    /// assert_eq!(loops.concat().iter().map(|e| e.0).collect::<Vec<_>>(), [0, 2]);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn forest_levels(&self, node: NodeIndex) -> Vec<Vec<&Spinney<B>>> {
        self.graph[node]
            .transitions(&self.wood)
            .iter()
            .map(|level| {
                level
                    .iter()
                    .flat_map(|transition| self.wood.transition_elements(transition))
                    .collect()
            })
            .collect()
    }

    /// The elements of the forest at `node`, in application order: the
    /// flattened [`Self::forest_levels`].
    ///
    /// # Panics
    ///
    /// Panics if `node` is not a forest of these forests.
    pub fn forest_spinneys(&self, node: NodeIndex) -> Vec<&Spinney<B>> {
        self.forest_levels(node).into_iter().flatten().collect()
    }

    /// The edges of the unfolded forests.
    pub fn edges(&self) -> impl Iterator<Item = ForestEdge> + '_ {
        self.graph
            .iter_edges()
            .filter_map(|(pair, _, operation)| match pair {
                HedgePair::Paired { source, sink } => Some(ForestEdge {
                    parent: self.graph.node_id(source),
                    child: self.graph.node_id(sink),
                    operation: *operation.data,
                }),
                _ => None,
            })
    }

    /// Forests that no further operation extends.
    pub fn terminal_nodes(&self) -> impl Iterator<Item = NodeIndex> + '_ {
        self.graph.iter_nodes().filter_map(|(node, mut crown, _)| {
            (!crown.any(|hedge| self.graph.flow(hedge).is_source())).then_some(node)
        })
    }

    /// The wood spinney whose transitions produced the forest at `node`.
    pub fn source_spinney(&self, node: NodeIndex) -> &Spinney<B> {
        self.wood.spinney(self.graph.source_node(node))
    }

    /// Whether the forest at `node` is the product of disjoint component forests.
    pub fn is_disjoint_union(&self, node: NodeIndex) -> bool {
        self.graph.is_disjoint_union(node)
    }

    /// The forest node of the root-to-`node` operation along wood edge `edge`
    /// applied at `order`: see [`Wood::current_given_pair`].
    pub fn current_given_pair(
        &self,
        edge: EdgeIndex,
        order: usize,
    ) -> (ForestNode<'_, B>, ForestNode<'_, B>) {
        self.wood.current_given_pair(edge, order)
    }

    /// The leaf operations of `node` with the forest each depends on.
    pub fn leaf_op_dependency_frontiers(
        &self,
        node: NodeIndex,
    ) -> impl Iterator<Item = (&Operation, NodeIndex)> + '_ {
        self.graph.leaf_op_dependency_frontiers(node, &self.wood)
    }

    /// A topological order of the forests in `subset`, starting at the root.
    ///
    /// `subset` must be a selection of the forest graph's half-edges (same
    /// width); this is checked.
    ///
    /// # Errors
    ///
    /// Fails if `subset` has another width or the selection is not acyclic.
    pub fn compatible_topological_order(
        &self,
        subset: &SuBitGraph,
    ) -> Result<Vec<NodeIndex>, ForestError> {
        let expected = self.graph.n_hedges();
        if subset.size() != expected {
            return Err(ForestError::SelectionWidth {
                found: subset.size(),
                expected,
            });
        }
        let mut order = self
            .graph
            .topo_sort_kahn_of(subset, DirectionBasis::Underlying)?;

        if let Some(root_position) = order.iter().position(|node| *node == self.root) {
            if root_position != 0 {
                let root = order.remove(root_position);
                order.insert(0, root);
            }
        } else {
            order.insert(0, self.root);
        }

        Ok(order)
    }

    /// A topological order of all forests.
    ///
    /// # Errors
    ///
    /// Fails if the forests are not acyclic.
    pub fn topological_order(&self) -> Result<Vec<NodeIndex>, ForestError> {
        Ok(self.graph.topo_sort_kahn(DirectionBasis::Underlying)?)
    }

    /// The component forests whose product is the disjoint union at `node`.
    ///
    /// # Errors
    ///
    /// Fails when `node` is not a factorized union of forests in this graph.
    pub fn disconnected_component_nodes(
        &self,
        node: NodeIndex,
    ) -> Result<Vec<NodeIndex>, ForestError> {
        let operation = &self.graph[node];
        let factors = self
            .wood
            .join_factors(self.graph.source_node(node))
            .ok_or_else(|| ForestError::NoComponentFactors {
                operation: operation.to_string(),
            })?;
        let target_ops = operation
            .key
            .iter_levels_top_down()
            .flat_map(|level| level.iter_leaf_ops())
            .cloned()
            .collect::<Vec<_>>();

        for op in &target_ops {
            let owners = factors
                .iter()
                .filter(|factor| factor.includes(&op.order))
                .count();
            if owners != 1 {
                return Err(ForestError::AmbiguousComponent {
                    op: op.order.string_label(),
                    operation: operation.to_string(),
                    owners,
                });
            }
        }

        factors
            .into_iter()
            .map(|factor| {
                let key = target_ops
                    .iter()
                    .filter(|op| factor.includes(&op.order))
                    .fold(TraceKey::empty(), |key, op| {
                        key.push(&self.wood, op.clone())
                    });
                self.graph
                    .iter_nodes()
                    .find_map(|(component, _, candidate)| {
                        (candidate.key == key && self.source_spinney(component).filter() == &factor)
                            .then_some(component)
                    })
                    .ok_or_else(|| ForestError::MissingComponent {
                        component: factor.string_label(),
                        operation: operation.to_string(),
                    })
            })
            .collect()
    }

    /// The node holding `operation` produced from the spinney it covers.
    ///
    /// # Errors
    ///
    /// Fails when no such node exists.
    pub fn operation_node_index(
        &self,
        operation: &OperationNode,
    ) -> Result<NodeIndex, ForestError> {
        let cover = operation.covers();
        self.graph
            .iter_nodes()
            .find_map(|(node, _, candidate)| {
                let source_matches = cover.as_ref().map_or_else(
                    || self.source_spinney(node).filter().is_empty(),
                    |cover| self.source_spinney(node).filter() == cover,
                );
                (candidate == operation && source_matches).then_some(node)
            })
            .ok_or_else(|| ForestError::MissingOperation {
                operation: operation.to_string(),
            })
    }

    /// The ways to replay `node` from integrated prefixes.
    ///
    /// `node` and every returned `integrated` index are nodes of
    /// [`Self::graph`]; `local_edges` are edges of [`Self::wood`]. Each state
    /// says: start from the integrated value already computed at `integrated`,
    /// then apply the listed wood transitions locally to reach `node`.
    ///
    /// For a connected forest these are the prefixes of its unique root path,
    /// ending with `node` itself and no local edges. For a disjoint union, one
    /// state of each component is chosen in every combination and their
    /// integrated prefixes are joined in Foata normal form; the joined prefix
    /// must itself be a forest of this graph. The empty forest replays from
    /// itself. A consumer must therefore have integrated values for every
    /// returned `integrated` node before using the states.
    ///
    /// # Errors
    ///
    /// Fails when a connected forest has no unique parent, when a union's
    /// components or a joined prefix are not forests of this graph.
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
    /// let forests = Wood::from_spinneys(spinneys, &spectacles, spectacles.loop_momentum_basis())?.unfold();
    /// let union = forests
    ///     .operations()
    ///     .find(|(node, _)| forests.is_disjoint_union(*node))
    ///     .expect("the two tadpoles factorize")
    ///     .0;
    /// // Each tadpole is integrated or applied locally: four combinations.
    /// let states = forests.union_replay_states(union)?;
    /// assert_eq!(states.len(), 4);
    /// assert!(states.iter().any(|state| state.integrated == forests.root()));
    /// assert!(states.iter().any(|state| state.integrated == union));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn union_replay_states(
        &self,
        node: NodeIndex,
    ) -> Result<Vec<UnionReplayState>, ForestError> {
        let operation = &self.graph[node];
        if operation.key.is_empty() {
            return Ok(vec![UnionReplayState {
                integrated: node,
                local_edges: Vec::new(),
            }]);
        }

        if !self.graph.is_disjoint_union(node) {
            let (parent, edge) =
                self.graph
                    .unique_parent(node)
                    .ok_or_else(|| ForestError::NoReplayParent {
                        operation: operation.to_string(),
                    })?;
            let mut states = self.union_replay_states(parent)?;
            for state in &mut states {
                state.local_edges.push(edge);
            }
            states.push(UnionReplayState {
                integrated: node,
                local_edges: Vec::new(),
            });
            return Ok(states);
        }

        let component_states = self
            .disconnected_component_nodes(node)?
            .into_iter()
            .map(|component| self.union_replay_states(component))
            .collect::<Result<Vec<_>, _>>()?;
        let mut states = Vec::new();
        for components in component_states
            .iter()
            .map(|states| states.iter())
            .multi_cartesian_product()
        {
            let key = TraceKey::try_foata_join(
                components
                    .iter()
                    .map(|component| &self.graph[component.integrated].key),
                &self.wood,
            )
            .ok_or_else(|| ForestError::IncompatiblePrefixes {
                operation: operation.to_string(),
            })?;
            let mut local_edges = Vec::new();
            for component in components {
                local_edges.extend_from_slice(&component.local_edges);
            }
            let integrated = self.operation_node_index(&OperationNode { key })?;
            states.push(UnionReplayState {
                integrated,
                local_edges,
            });
        }

        Ok(states)
    }

    /// Write the forests as DOT. Each node carries `label`, the Foata levels in
    /// `foata` and the covered subgraph in `cover`; `label(node)` supplies the
    /// label text (Typst code for the `tree-auto-root` template).
    ///
    /// # Errors
    ///
    /// Returns the writer's error.
    pub fn dot_serialize_with(
        &self,
        writer: &mut impl Write,
        mut label: impl FnMut(NodeIndex) -> String,
    ) -> fmt::Result {
        let empty = self.graph.empty_subgraph::<SuBitGraph>();
        let attrs = self
            .graph
            .iter_nodes()
            .map(|(node, _, operation)| {
                let cover = operation
                    .covers()
                    .unwrap_or_else(|| empty.clone())
                    .string_label();
                (
                    operation.clone(),
                    format!(
                        "label={} foata={} cover={}",
                        dot_attr_value(&label(node)),
                        dot_attr_value(&operation.foata_level_labels()),
                        dot_attr_value(&cover),
                    ),
                )
            })
            .collect::<AHashMap<_, _>>();

        self.graph.dot_impl_fmt(
            writer,
            &self.graph.full_filter(),
            "start=2;\n",
            &|_| None,
            &|_| None,
            &|operation| Some(attrs[operation].clone()),
        )
    }

    /// [`Self::dot_serialize_with`], labelling each forest by its Foata form as a
    /// Typst string.
    pub fn dot_serialize(&self) -> String {
        let mut output = String::new();
        self.dot_serialize_with(&mut output, |node| {
            format!("\"{}\"", typst_string_escape(&self.graph[node].to_string()))
        })
        .expect("writing hedge-poset forest DOT into a string must succeed");
        output
    }
}

/// Quote raw text as a complete DOT attribute value.
fn dot_attr_value(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len() + 2);
    escaped.push('"');
    for c in value.chars() {
        match c {
            '\\' => escaped.push_str("\\\\"),
            '"' => escaped.push_str("\\\""),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            _ => escaped.push(c),
        }
    }
    escaped.push('"');
    escaped
}

fn typst_string_escape(value: &str) -> String {
    value.replace('\\', "\\\\").replace('"', "\\\"")
}
