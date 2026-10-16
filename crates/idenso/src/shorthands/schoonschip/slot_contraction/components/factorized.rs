//! State merging on the collector's existing port graph.
//!
//! Remaining factors are opaque terminals of the same component reducer used
//! for a completed monomial. Only the selected factor contributes alternatives;
//! substitutions remain port bindings until a distinct output variable is emitted.

use super::*;
use crate::tensor::{
    ContractSettings, SymbolicNet, simplification::observation::DomainObservations,
};
use ahash::AHashSet;
use linnet::{
    half_edge::{
        NodeIndex,
        subgraph::{Inclusion, SuBitGraph},
        tree::SimpleTraversalTree,
    },
    tree::child_vec::ChildVecStore,
};
use spenso::{
    network::{
        graph::{NetworkEdge, NetworkLeaf, NetworkNode, NetworkOp},
        library::DummyLibrary,
        parsing::{LeafInterfaceCache, ParseSettings, ShorthandParsing},
        store::{TensorScalarStore, TensorScalarStoreMapping},
    },
    structure::{
        OrderedStructure, TensorStructure,
        abstract_index::AbstractIndex,
        partial::{PartialIndex, PartialStructureExt},
        slot::IsAbstractSlot,
    },
};
use std::sync::Arc;
use symbolica::{atom::AtomOrView, coefficient::CoefficientView};

// A state contains occurrence-local port bindings, not rebuilt expressions.
// Variables and the incidence scratch still belong to ComponentSum.
#[derive(Clone, PartialEq, Eq, Hash)]
struct Remaining<'a> {
    factors: Vec<(usize, Vec<Argument<'a>>)>,
    residual: Vec<(usize, u32)>,
}

impl Remaining<'_> {
    #[cfg(feature = "reference-cases")]
    fn storage_bytes(&self) -> usize {
        std::mem::size_of::<Self>()
            + self.factors.capacity() * std::mem::size_of::<(usize, Vec<Argument<'_>>)>()
            + self
                .factors
                .iter()
                .map(|(_, arguments)| arguments.capacity() * std::mem::size_of::<Argument<'_>>())
                .sum::<usize>()
            + self.residual.capacity() * std::mem::size_of::<(usize, u32)>()
    }
}

type FactorTerm = (Vec<(usize, u32)>, Rational);
// Selected expression/interface and retained factors.
type PrerequisiteFactors = (Atom, PartialStructure, Vec<Atom>);

pub(crate) use crate::tensor::simplification::ReductionStatus;

/// Exact output may retain unsupported contractions for subsequent operations.
/// Algebraic validity and completion are separate facts.
#[derive(Clone)]
pub(crate) struct FactorizedContraction<Root = Atom> {
    pub(crate) root: Root,
    pub(crate) status: ReductionStatus,
    pub(crate) observations: Option<DomainObservations>,
}

// A parse belongs to one arithmetic scope. Opening an occurrence is cached on
// that scope, so borrowed payloads and its port atoms outlive every DP state.
// This is parser scratch, not another symbolic tensor representation.
struct ContractionScope<'src> {
    source: AtomView<'src>,
    network: SymbolicNet<AbstractIndex, AtomOrView<'src>>,
    tree: SimpleTraversalTree<ChildVecStore<()>>,
    root: NodeIndex,
    ports: AHashMap<NodeIndex, (Vec<Atom>, PartialStructure)>,
    opened: AHashMap<NodeIndex, std::cell::OnceCell<Option<Box<Self>>>>,
    reductions: AHashMap<NodeIndex, std::cell::OnceCell<Option<FactorizedContraction>>>,
    // Each scope serves one fixed contraction or admission policy.
    work: AHashMap<NodeIndex, std::cell::OnceCell<bool>>,
    interfaces: Option<LeafInterfaceCache>,
}

impl<'src> ContractionScope<'src> {
    fn parse(source: AtomView<'src>) -> Option<Self> {
        let network = SlotContraction::shallow_network(source, None).ok()?;
        Self::from_network(source, network, None)
    }

    fn from_network(
        source: AtomView<'src>,
        network: SymbolicNet<AbstractIndex, AtomOrView<'src>>,
        interfaces: Option<LeafInterfaceCache>,
    ) -> Option<Self> {
        let tree = network.graph.expr_tree().cast::<ChildVecStore<()>>();
        let root = network.graph.graph.node_id(network.graph.head());
        let nodes = tree
            .iter_preorder_tree_nodes(&network.graph.graph, root)
            .collect::<Vec<_>>();
        let external: SuBitGraph = network.graph.graph.external_filter();
        // Preorder makes each arithmetic subtree a contiguous interval. Reuse
        // those intervals for every boundary instead of allocating a descendant
        // set and traversing the same tree for each node.
        let positions = nodes
            .iter()
            .copied()
            .enumerate()
            .map(|(position, node)| (node, position))
            .collect::<AHashMap<_, _>>();
        let mut ends = (1..=nodes.len()).collect::<Vec<_>>();
        for (position, &node) in nodes.iter().enumerate().rev() {
            for child in tree.iter_children(node, &network.graph.graph) {
                ends[position] = ends[position].max(ends[positions[&child]]);
            }
        }
        // The parser already observed each leaf's written ports. Retain their
        // literal index payloads when describing connected boundaries instead
        // of rebuilding names through user-defined normalizers.
        let mut written_slots = AHashMap::new();
        if let Some(interfaces) = &interfaces {
            let interfaces = interfaces.lock().unwrap();
            let mut matcher = SlotMatcher::default();
            for &node in &nodes {
                if let NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) =
                    &network.graph.graph[node]
                    && let Some((_, ports)) = interfaces.get(
                        network.store.tensors[*index]
                            .expression
                            .as_atom_view()
                            .get_data(),
                    )
                {
                    for atom in ports {
                        if let Ok(slot) = matcher.parse::<LibraryRep, AbstractIndex>(atom.as_view())
                        {
                            written_slots.entry(slot).or_insert_with(|| atom.clone());
                        }
                    }
                }
            }
        }
        let mut ports = AHashMap::with_capacity(nodes.len());
        for (begin, &node) in nodes.iter().enumerate() {
            let end = ends[begin];
            let mut boundary = Vec::new();
            for &member in &nodes[begin..end] {
                for hedge in network.graph.graph.iter_crown(member) {
                    let NetworkEdge::Slot(slot) = network.graph.graph[[&hedge]] else {
                        continue;
                    };
                    if !external.includes(&hedge)
                        && positions
                            .get(&network.graph.graph.node_id(network.graph.graph.inv(hedge)))
                            .is_some_and(|&other| begin <= other && other < end)
                    {
                        continue;
                    }
                    // Sewing shares one edge label; the occurrence's storage
                    // axis still owns this endpoint's variance.
                    let slot = match &network.graph.graph[member] {
                        NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => {
                            *network.store.tensors[*index]
                                .structure
                                .external_structure()
                                .get(usize::from(network.graph.slot_order[hedge.0]))?
                        }
                        _ => slot,
                    };
                    boundary.push(slot);
                }
            }
            boundary.sort();
            ports.insert(
                node,
                (
                    boundary
                        .iter()
                        .map(|slot| {
                            written_slots
                                .get(slot)
                                .cloned()
                                .unwrap_or_else(|| slot.to_atom())
                        })
                        .collect(),
                    PartialStructure::from_logical_slots(
                        boundary
                            .into_iter()
                            .map(|slot| slot.rep().slot(PartialIndex::Explicit(slot.aind()))),
                    ),
                ),
            );
        }
        Some(Self {
            source,
            network,
            tree,
            root,
            ports,
            opened: nodes
                .iter()
                .map(|&node| (node, Default::default()))
                .collect(),
            reductions: nodes
                .iter()
                .map(|&node| (node, Default::default()))
                .collect(),
            work: nodes
                .iter()
                .map(|&node| (node, Default::default()))
                .collect(),
            interfaces,
        })
    }

    fn payload(&self, node: NodeIndex) -> Option<&AtomOrView<'src>> {
        match self.network.graph.graph[node] {
            NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => {
                Some(&self.network.store.tensors[index].expression)
            }
            NetworkNode::Leaf(NetworkLeaf::Scalar(index)) => {
                Some(self.network.store.get_scalar_ref(index))
            }
            _ => None,
        }
    }

    fn literal(&self, node: NodeIndex) -> Option<AtomView<'_>> {
        self.payload(node)
            .map(AtomOrView::as_view)
            .or_else(|| (node == self.root).then_some(self.source))
    }

    fn open(&self, node: NodeIndex) -> Option<&Self> {
        self.opened[&node]
            .get_or_init(|| {
                // Composite arithmetic leaves retain their original view. The
                // parser's owned scalar coefficients contain no tensor interior.
                let AtomOrView::View(value) = self.payload(node)? else {
                    return None;
                };
                #[cfg(test)]
                crate::tensor::SELECTED_OPENS.with(|count| count.set(count.get() + 1));
                let network =
                    SlotContraction::shallow_network(*value, self.interfaces.as_ref()).ok()?;
                let opened = Self::from_network(*value, network, self.interfaces.clone())?;
                if matches!(
                    opened.network.graph.graph[opened.root],
                    NetworkNode::Leaf(_)
                ) && opened.literal(opened.root)? == *value
                {
                    // Unsupported powers and opaque scalar syntax do not expose
                    // arithmetic. Decline instead of reopening the same leaf.
                    return None;
                }
                Some(Box::new(opened))
            })
            .as_deref()
    }

    fn contains_work(
        &self,
        node: NodeIndex,
        contractor: &SlotContraction,
        settings: ContractSettings<'_>,
    ) -> bool {
        *self.work[&node].get_or_init(|| {
            let Some(value) = self.literal(node) else {
                return false;
            };
            if let Some(observed) = &contractor.observations {
                observed.region(value).is_none_or(|region| {
                    // Internal colour/gamma connections are not metric/vector
                    // work. Keep these regions opaque, including their open
                    // interfaces: incident sources can bind those ports without
                    // opening the payload. Minimal notation admission still
                    // needs its independent chain/trace candidates. Expanding
                    // collection also merges alpha-equivalent tensor terms.
                    (settings.expand
                        || settings.collect_chains
                        || settings.collect_traces
                        || !region.excludes_contraction_sources(settings))
                        && (region.has_indices() || region.counts[7] > 0)
                })
            } else {
                let observed = SimplificationCandidates::scan(value, [contractor.metric], || true);
                observed.symbols[0] || observed.dots || observed.repeated_indices
            }
        })
    }

    fn factors(&self, root: NodeIndex) -> Vec<NodeIndex> {
        if matches!(
            self.network.graph.graph[root],
            NetworkNode::Op(NetworkOp::Product)
        ) {
            let mut nodes = self
                .tree
                .iter_children(root, &self.network.graph.graph)
                .collect::<Vec<_>>();
            nodes.reverse();
            nodes
        } else {
            vec![root]
        }
    }
}

impl SlotContraction {
    fn contraction_settings(&self, rank_one: bool) -> ContractSettings<'_> {
        ContractSettings {
            metrics: self.metrics,
            representations: self.representations.as_deref(),
            rank_one,
            collect_chains: false,
            collect_traces: false,
            expand: self.expand,
            ..Default::default()
        }
    }

    /// Budget admission uses the same shallow scopes and boundary edges as the
    /// contraction. It does not normalize a tensor or invoke user callbacks.
    pub(crate) fn minimal_work_permitted(
        &self,
        source: AtomView<'_>,
        settings: crate::tensor::ContractSettings<'_>,
    ) -> Option<bool> {
        let scope = self.scope(source)?;
        self.minimal_scope_work(&scope, scope.root, settings)
    }

    fn minimal_scope_work(
        &self,
        scope: &ContractionScope<'_>,
        root: NodeIndex,
        settings: crate::tensor::ContractSettings<'_>,
    ) -> Option<bool> {
        match scope.network.graph.graph[root] {
            NetworkNode::Op(NetworkOp::Sum) => {
                for child in scope.tree.iter_children(root, &scope.network.graph.graph) {
                    if self.minimal_scope_work(scope, child, settings)? {
                        return Some(true);
                    }
                }
                return Some(false);
            }
            NetworkNode::Op(NetworkOp::Power(_)) => {
                let child = scope
                    .tree
                    .iter_children(root, &scope.network.graph.graph)
                    .next()?;
                if !scope.ports[&child].0.is_empty()
                    && !self.tensor_alternatives(scope.literal(child)?)
                {
                    // An atomic tensor power can connect copies directly.
                    let sources = self
                        .observations
                        .as_ref()
                        .and_then(|observed| observed.region(scope.literal(child)?))
                        .map(|region| region.contraction_sources(settings))?;
                    return Some(
                        (sources || (settings.collect_traces && scope.ports[&child].0.len() == 2))
                            && scope.ports[&child]
                                .1
                                .logical_slots()
                                .iter()
                                .any(|slot| self.permits(slot.rep().rep)),
                    );
                }
                return self.minimal_scope_work(scope, child, settings);
            }
            _ => {}
        }
        let nodes = scope.factors(root);
        let external: SuBitGraph = scope.network.graph.graph.external_filter();
        let components = scope.network.graph.slot_components(&scope.tree, &nodes);
        for &node in &nodes {
            let literal = scope.literal(node)?;
            let component = components
                .iter()
                .find(|component| component.contains(&node))?;
            let alternatives = component
                .iter()
                .copied()
                .filter(|&member| {
                    !scope.ports[&member].0.is_empty()
                        && scope
                            .literal(member)
                            .is_some_and(|source| self.tensor_alternatives(source))
                })
                .collect::<Vec<_>>();
            let open_component = component.iter().any(|&member| {
                scope.network.graph.graph.iter_crown(member).any(|edge| {
                    matches!(scope.network.graph.graph[[&edge]], NetworkEdge::Slot(_))
                        && external.includes(&edge)
                })
            });
            let matrix = scope.ports[&node].0.len() == 2
                && (settings.collect_traces || (settings.collect_chains && open_component));
            let connected = scope.network.graph.graph.iter_crown(node).any(|edge| {
                matches!(scope.network.graph.graph[[&edge]], NetworkEdge::Slot(slot)
                    if self.permits(slot.rep().rep))
                    && !external.includes(&edge)
                    && component.contains(
                        &scope
                            .network
                            .graph
                            .graph
                            .node_id(scope.network.graph.graph.inv(edge)),
                    )
            });
            if alternatives.len() == 1
                && alternatives[0] == node
                && connected
                && (matrix
                    || self
                        .observations
                        .as_ref()
                        .and_then(|observed| observed.region(literal))?
                        .contraction_sources(settings))
            {
                // One selected existing sum may supply a metric/vector to an
                // atomic partner. It never combines independent alternatives.
                return Some(true);
            }
            if let AtomView::Fun(function) = literal {
                let head = function.get_symbol();
                let source = (settings.metrics && head == self.metric)
                    || (settings.rank_one && head.has_tag(&self.tags.rank1));
                if (source || matrix) && connected {
                    return Some(true);
                }
                if !source
                    && scope.ports[&node].0.is_empty()
                    && settings.collect_traces
                    && self
                        .observations
                        .as_ref()
                        .and_then(|observed| observed.region(literal))?
                        .has_indices()
                {
                    // An internally closed matrix leaf has no graph boundary;
                    // its notation owner establishes trace eligibility.
                    return None;
                }
                if head == self.tags.chain || head == self.tags.trace || head == self.tags.bracket {
                    // These notation owners carry ports inside their words.
                    // Keep their existing admission and checked rewrite rules.
                    return None;
                }
            } else if scope.contains_work(node, self, settings) {
                if !matches!(scope.network.graph.graph[node], NetworkNode::Leaf(_)) {
                    if self.minimal_scope_work(scope, node, settings)? {
                        return Some(true);
                    }
                } else {
                    let opened = scope.open(node)?;
                    if self.minimal_scope_work(opened, opened.root, settings)? {
                        return Some(true);
                    }
                }
            }
        }
        Some(false)
    }

    fn tensor_alternatives(&self, source: AtomView<'_>) -> bool {
        if self
            .observations
            .as_ref()
            .and_then(|observed| observed.region(source))
            .is_some_and(|region| !region.has_indices())
        {
            return false;
        }
        match source {
            AtomView::Add(_) => true,
            AtomView::Mul(product) => product
                .iter()
                .any(|factor| self.tensor_alternatives(factor)),
            AtomView::Pow(power) => self.tensor_alternatives(power.get_base_exp().0),
            _ => false,
        }
    }

    fn scope<'a>(&'a self, source: AtomView<'a>) -> Option<ContractionScope<'a>> {
        if let Some((planned, graph)) = &self.graph
            && planned.as_view() == source
        {
            let network = graph.map_ref(
                |scalar| AtomOrView::View(scalar.as_view()),
                |tensor| SymbolicTensor {
                    expression: AtomOrView::View(tensor.expression.as_view()),
                    structure: tensor.structure.clone(),
                    is_metric: tensor.is_metric,
                    is_composite: tensor.is_composite,
                    proofs: tensor.proofs.clone(),
                },
            );
            ContractionScope::from_network(source, network, self.interfaces.clone())
        } else {
            ContractionScope::parse(source)
        }
    }

    /// The shared depth-one structural view keeps arithmetic leaves borrowed.
    /// Family collection opens selected leaves through its existing tape reader.
    pub(crate) fn shallow_network<'src>(
        source: AtomView<'src>,
        interfaces: Option<&LeafInterfaceCache>,
    ) -> Result<
        SymbolicNet<AbstractIndex, AtomOrView<'src>>,
        crate::tensor::inference::TensorInferenceError,
    > {
        #[cfg(feature = "reference-cases")]
        let _phase = crate::reference_cases::timing::scope(
            crate::reference_cases::timing::Phase::ShallowGraph,
        );
        let settings = ParseSettings {
            precontract_scalars: false,
            depth_limit: Some(1),
            shorthand_parsing: ShorthandParsing::Opaque,
            parse_composite_scalars_as_tensors: true,
            ..ParseSettings::default()
        };
        type Tensor<'src> = SymbolicTensor<OrderedStructure, AtomOrView<'src>>;
        let library = DummyLibrary::<Tensor<'src>>::new();
        if let Some(interfaces) = interfaces {
            SymbolicNet::<AbstractIndex, AtomOrView<'src>>::try_from_admitted_view::<
                OrderedStructure,
                _,
            >(source, &library, &settings, Arc::clone(interfaces))
        } else {
            SymbolicNet::<AbstractIndex, AtomOrView<'src>>::try_from_view::<OrderedStructure, _>(
                source, &library, &settings,
            )
        }
        .map_err(|error| crate::tensor::inference::TensorInferenceError::Invalid(error.to_string()))
    }

    /// Collect notation on connected scopes of the same shallow graph. Scalar
    /// spectators and components without a selected representation are reused.
    pub(crate) fn collect_matrix_notation(
        &self,
        source: AtomView<'_>,
        representations: &[LibraryRep],
        collect: bool,
        traces: bool,
    ) -> Option<Atom> {
        use crate::shorthands::chain::Chain;
        let scope = self.scope(source)?;
        let nodes = scope.factors(scope.root);
        let components = scope.network.graph.slot_components(&scope.tree, &nodes);
        let external: SuBitGraph = scope.network.graph.graph.external_filter();
        let mut output = Vec::new();
        let mut changed = false;
        for component in components {
            let factors = component
                .iter()
                .map(|&node| scope.literal(node))
                .collect::<Option<Vec<_>>>()?;
            let mut expression = if factors.len() == 1 {
                factors[0].to_owned()
            } else {
                Atom::mul_many(factors.iter())
            };
            for &representation in representations {
                let connection = component.len() > 1
                    // An opaque arithmetic leaf can own connections that the
                    // outer graph intentionally does not expose. Regional
                    // representation facts below select its local notation.
                    || factors.iter().any(|factor| {
                        matches!(factor, AtomView::Add(_) | AtomView::Mul(_) | AtomView::Pow(_))
                            || matches!(factor, AtomView::Fun(function)
                                if function.get_symbol() == self.tags.chain)
                    })
                    || (traces
                        && component.iter().any(|&node| {
                            scope
                                .tree
                                .iter_preorder_tree_nodes(&scope.network.graph.graph, node)
                                .any(|member| {
                                    scope.network.graph.graph.iter_crown(member).any(|edge| {
                                        let NetworkEdge::Slot(slot) =
                                            scope.network.graph.graph[[&edge]]
                                        else {
                                            return false;
                                        };
                                        slot.rep().rep.base() == representation.base()
                                            && !external.includes(&edge)
                                            && scope
                                                .network
                                                .graph
                                                .graph
                                                .node_id(scope.network.graph.graph.inv(edge))
                                                == member
                                    })
                                })
                        }));
                let relevant = connection
                    && factors.iter().any(|&factor| {
                        self.observations
                            .as_ref()
                            .and_then(|observed| observed.region(factor))
                            .is_none_or(|region| region.has_chain_channel(representation))
                    });
                if relevant {
                    let next = expression.collect_chains(representation, collect, traces, false);
                    changed |= next != expression;
                    expression = next;
                }
            }
            output.push(expression);
        }
        Some(if changed {
            Atom::mul_many(output)
        } else {
            source.to_owned()
        })
    }

    /// Select only a family's factors and metric/vector paths incident to them.
    /// Connectivity authorizes these prerequisites, never a foreign tensor body.
    pub(crate) fn prerequisite_factors(
        &self,
        source: AtomView<'_>,
        mut selected: impl FnMut(AtomView<'_>) -> bool,
        allowed_vector: &dyn Fn(LibraryRep) -> bool,
    ) -> Option<PrerequisiteFactors> {
        let scope = self.scope(source)?;
        let nodes = scope.factors(scope.root);
        let graph = &scope.network.graph.graph;
        let mut owners = AHashMap::new();
        for (position, &node) in nodes.iter().enumerate() {
            for member in scope.tree.iter_preorder_tree_nodes(graph, node) {
                owners.insert(member, position);
            }
        }
        let mut active = nodes
            .iter()
            .map(|&node| scope.literal(node).is_some_and(&mut selected))
            .collect::<Vec<_>>();
        if !active.iter().any(|active| *active) {
            return None;
        }
        let eligible = nodes
            .iter()
            .map(|&node| {
                matches!(scope.literal(node), Some(AtomView::Fun(function))
                if (self.metrics && function.get_symbol() == self.metric)
                    || (function.get_symbol().has_tag(&self.tags.rank1)
                        && scope.ports[&node].1.logical_slots().iter()
                            .all(|slot| allowed_vector(slot.rep().rep))))
            })
            .collect::<Vec<_>>();
        let mut joined = false;
        let mut touched = vec![false; nodes.len()];
        loop {
            let mut changed = false;
            for (position, &node) in nodes.iter().enumerate() {
                if active[position] || !eligible[position] {
                    continue;
                }
                let mut attached = false;
                for member in scope.tree.iter_preorder_tree_nodes(graph, node) {
                    for edge in graph.iter_crown(member) {
                        let NetworkEdge::Slot(slot) = graph[[&edge]] else {
                            continue;
                        };
                        if self.permits(slot.rep().rep)
                            && let Some(&other) = owners.get(&graph.node_id(graph.inv(edge)))
                            && other != position
                            && active[other]
                        {
                            touched[other] = true;
                            attached = true;
                        }
                    }
                }
                if attached {
                    active[position] = true;
                    touched[position] = true;
                    changed = true;
                    joined = true;
                }
            }
            if !changed {
                break;
            }
        }
        if !joined {
            return None;
        }
        let mut factors = Vec::new();
        let mut interfaces = Vec::new();
        let mut spectators = Vec::new();
        for (position, &node) in nodes.iter().enumerate() {
            let value = scope.literal(node)?.to_owned();
            if active[position] && touched[position] {
                // A selected but disconnected scalar child is still opaque:
                // another family's occurrence is no reason to open it here.
                factors.push(value);
                interfaces.push(scope.ports[&node].1.clone());
            } else {
                spectators.push(value);
            }
        }
        let interface = InterfaceInference::merge_explicit_interface_sequence(&interfaces).ok()?;
        Some((Atom::mul_many(factors), interface, spectators))
    }

    /// Plan an already admitted interface; the typed caller owns source validation.
    pub(crate) fn contract_factorized(
        &self,
        source: AtomView<'_>,
        order: Option<&[usize]>,
        rank_one: bool,
    ) -> Option<FactorizedContraction> {
        if rank_one && self.observations.is_none() && source.contains_symbol(self.tags.chain) {
            self.dummy_state().reserve_indices(source);
        }
        let intrinsic = self
            .graph
            .as_ref()
            .filter(|(planned, _)| planned.as_view() == source)
            .and(self.observations.as_ref())
            .map_or_else(
                || InterfaceInference::normalization_is_intrinsic(source),
                |observed| {
                    observed.candidates.intrinsic
                        || (!observed.candidates.complete
                            && InterfaceInference::normalization_is_intrinsic(source))
                },
            );
        if source.needs_normalization() || self.metric.get_evaluation_info().is_some() || !intrinsic
        {
            return None;
        }
        let scope = self.scope(source)?;
        self.contract_graph(&scope, scope.root, order, rank_one)
    }

    fn power_scope(
        &self,
        scope: &ContractionScope<'_>,
        node: NodeIndex,
    ) -> Option<(NodeIndex, i8, bool)> {
        let NetworkNode::Op(NetworkOp::Power(power)) = scope.network.graph.graph[node] else {
            return None;
        };
        let mut children = scope.tree.iter_children(node, &scope.network.graph.graph);
        let child = children.next()?;
        if children.next().is_some() {
            return None;
        }
        let source = scope.literal(child)?;
        if !scope.ports[&child].0.is_empty() && (self.expand || !self.tensor_alternatives(source)) {
            return None;
        }
        use spenso::network::parsing::AtomStructureExt;
        let internal = self
            .observations
            .as_ref()
            .and_then(|observed| observed.region(source))
            .map_or_else(
                || source.has_repeated_explicit_indices(),
                |region| region.has_indices(),
            );
        Some((child, power, internal))
    }

    fn scoped_power<'a>(
        &self,
        scope: &'a ContractionScope<'_>,
        node: NodeIndex,
        rank_one: bool,
    ) -> Option<&'a FactorizedContraction> {
        let (child, power, internal) = self.power_scope(scope, node)?;
        scope.reductions[&node]
            .get_or_init(|| {
                if !internal {
                    // A transparent bracket can remain after local dot pairing
                    // has removed its last explicit index. Lower just that
                    // selected wrapper, retaining ordinary scalar powers verbatim.
                    if matches!(scope.literal(child)?, AtomView::Fun(function)
                        if function.get_symbol() == self.tags.bracket)
                        && scope.contains_work(child, self, self.contraction_settings(rank_one))
                    {
                        let opened = scope.open(child)?;
                        let mut result =
                            self.contract_graph(opened, opened.root, None, rank_one)?;
                        result.root = result.root.pow(power);
                        return Some(result);
                    }
                    return Some(FactorizedContraction {
                        root: scope.literal(node)?.to_owned(),
                        status: ReductionStatus::Complete,
                        observations: self.observations.as_ref().and_then(|observed| {
                            observed.certify_scalar_region(scope.literal(node)?)
                        }),
                    });
                }
                // Contract the base in its own scope before applying the
                // exponent. Closed bases keep independent dummy copies; minimal
                // open sums retain their alternatives and external connections.
                let mut result = self.contract_graph(scope, child, None, rank_one)?;
                result.root = result.root.pow(power);
                Some(result)
            })
            .as_ref()
    }

    fn factor_order(
        &self,
        scope: &ContractionScope<'_>,
        nodes: &[NodeIndex],
        mut selected: Vec<usize>,
        order: Option<&[usize]>,
    ) -> Option<Vec<usize>> {
        if let Some(order) = order {
            let factors = match scope.source {
                AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
                value => vec![value],
            };
            let mut positions = order.to_vec();
            positions.sort_unstable();
            if positions != (0..factors.len()).collect::<Vec<_>>() {
                return None;
            }
            let ranks = selected
                .iter()
                .map(|&position| {
                    let literal = scope.literal(nodes[position])?;
                    let ordinal = factors.iter().position(|&factor| factor == literal)?;
                    Some((position, order.iter().position(|&entry| entry == ordinal)?))
                })
                .collect::<Option<AHashMap<_, _>>>()?;
            selected.sort_by_key(|position| ranks[position]);
            return Some(selected);
        }
        // Minimize the live boundary of the accumulated factor set. Incidence
        // comes entirely from the depth-one graph; scoring never opens a leaf.
        let graph = &scope.network.graph.graph;
        let edges = nodes
            .iter()
            .map(|&node| {
                scope
                    .tree
                    .iter_preorder_tree_nodes(graph, node)
                    .flat_map(|member| graph.iter_crown(member))
                    .filter(|hedge| scope.network.graph.graph[[hedge]].is_slot())
                    .filter(|&hedge| {
                        graph.inv(hedge) == hedge || graph.node_id(graph.inv(hedge)) != node
                    })
                    .map(|hedge| hedge.0.min(graph.inv(hedge).0))
                    .collect::<std::collections::BTreeSet<_>>()
            })
            .collect::<Vec<_>>();
        let mut boundary = std::collections::BTreeSet::new();
        let mut planned = Vec::with_capacity(selected.len());
        while !selected.is_empty() {
            let best = (0..selected.len()).min_by_key(|&i| {
                let position = selected[i];
                (
                    boundary.symmetric_difference(&edges[position]).count(),
                    std::cmp::Reverse(boundary.intersection(&edges[position]).count()),
                    position,
                )
            })?;
            let position = selected.remove(best);
            boundary = boundary
                .symmetric_difference(&edges[position])
                .copied()
                .collect();
            planned.push(position);
        }
        Some(planned)
    }

    /// Compile selected, admitted arithmetic on the existing contraction tape.
    /// The outer shallow graph already supplied topology and factor order.
    /// Selected sums and products use the same tape as their parsed graph;
    /// indexed powers and transparent brackets retain the scoped graph path.
    fn compile_arithmetic<'a>(
        &self,
        value: AtomView<'a>,
        tape: &mut TermTape<InputLeaf<'a>>,
        sum: &mut ComponentSum<'a, '_>,
    ) -> Option<(usize, (usize, usize))> {
        if let AtomView::Num(number) = value {
            return if let Ok(coefficient) = Rational::try_from(value) {
                tape.number(InputLeaf::Atom(value), coefficient)
            } else if matches!(
                number.get_coeff_view(),
                CoefficientView::Natural(..) | CoefficientView::Large(..)
            ) {
                tape.leaf(InputLeaf::Scalar(value))
            } else {
                None
            };
        }
        if self
            .observations
            .as_ref()
            .and_then(|observed| observed.region(value))
            .is_some_and(|region| region.certifies_scalar_interface_syntax())
        {
            sum.scalar_interfaces.insert(value, true);
            return tape.leaf(InputLeaf::Scalar(value));
        }
        match value {
            AtomView::Add(add) => {
                let children = add
                    .iter()
                    .map(|child| self.compile_arithmetic(child, tape, sum))
                    .collect::<Option<Vec<_>>>()?;
                tape.group(children.into_iter(), true)
            }
            AtomView::Mul(product) => {
                let children = product
                    .iter()
                    .map(|child| self.compile_arithmetic(child, tape, sum))
                    .collect::<Option<Vec<_>>>()?;
                tape.group(children.into_iter(), false)
            }
            AtomView::Pow(_) => None,
            _ => tape.leaf(InputLeaf::Atom(value)),
        }
    }

    fn compile_scope<'a>(
        &self,
        scope: &'a ContractionScope<'_>,
        node: NodeIndex,
        tape: &mut TermTape<InputLeaf<'a>>,
        sum: &mut ComponentSum<'a, '_>,
        positions: &mut AHashMap<(usize, NodeIndex), usize>,
        scoped: &mut AHashMap<(usize, NodeIndex), &'a FactorizedContraction>,
    ) -> Option<(usize, (usize, usize))> {
        tape.compile_graph(
            &scope.network.graph,
            &scope.tree,
            node,
            &mut |tape, node, value| {
                if let Some(literal) = scope.literal(node)
                    && matches!(literal, AtomView::Add(_) | AtomView::Mul(_) | AtomView::Pow(_))
                    && !scope.contains_work(node, self, self.contraction_settings(sum.rank_one))
                {
                    // The graph already established this boundary. A completed
                    // structural region can still contain identity candidates;
                    // retain its payload without inventing terminal algebra facts.
                    sum.scalar_interfaces.insert(literal, scope.ports[&node].0.is_empty());
                    return Some(if scope.ports[&node].0.is_empty() {
                        TermLeaf::Value(InputLeaf::Scalar(literal))
                    } else {
                        let position = sum.register_factor(scope, node, positions)?;
                        TermLeaf::Value(InputLeaf::Subtree(position))
                    });
                }
                if let Some((child, power, _)) = self.power_scope(scope, node) {
                    let result = self.scoped_power(scope, node, sum.rank_one)?;
                    let unchanged = scope.literal(node).map_or_else(
                        || {
                            scope
                                .literal(child)
                                .is_some_and(|base| result.root == base.pow(power))
                        },
                        |source| result.root.as_view() == source,
                    );
                    // Visiting an unchanged scalar power is bookkeeping, not
                    // a contraction. It must not disable exact source reuse.
                    if !unchanged
                        || result.status != ReductionStatus::Complete
                    {
                        scoped.insert((scope as *const _ as usize, node), result);
                    }
                    if let Some(observed) = &result.observations {
                        sum.scalar_forms.insert(
                            result.root.as_view(),
                            observed.needs_dot_normalization(),
                        );
                    }
                    return Some(if scope.ports[&node].0.is_empty() {
                        TermLeaf::Value(InputLeaf::Scalar(result.root.as_view()))
                    } else {
                        let position = sum.register_factor(scope, node, positions)?;
                        sum.opaque_factors[position].0 = result.root.as_view();
                        TermLeaf::Value(InputLeaf::Subtree(position))
                    });
                }
                if !matches!(value, NetworkNode::Leaf(_)) {
                    return None;
                }
                let literal = scope.literal(node)?;
                sum.scalar_interfaces
                    .insert(literal, scope.ports[&node].0.is_empty());
                if matches!(
                    literal,
                    AtomView::Add(_) | AtomView::Mul(_) | AtomView::Pow(_)
                ) || matches!(literal, AtomView::Fun(function) if function.get_symbol() == self.tags.bracket) {
                    if scope.contains_work(node, self, self.contraction_settings(sum.rank_one)) {
                        let region = self.observations.as_ref()
                            .and_then(|observed| observed.region(literal));
                        if self.expand
                            && matches!(literal, AtomView::Add(_) | AtomView::Mul(_))
                            && region.as_ref()
                                .is_some_and(|region| !region.has_indexed_powers() && region.counts[7] == 0)
                        {
                            let (node, size) = self.compile_arithmetic(literal, tape, sum)?;
                            return Some(TermLeaf::Reference(node, size));
                        }
                        #[cfg(feature = "reference-cases")]
                        {
                            use crate::reference_cases::timing::{ScopeOpenReason, scope_open};
                            let reason = match region.as_ref() {
                                None => ScopeOpenReason::Unobserved,
                                Some(region) if region.contraction_sources(self.contraction_settings(sum.rank_one)) => ScopeOpenReason::MetricOrVector,
                                Some(_) => ScopeOpenReason::OpaqueOrBracket,
                            };
                            scope_open(reason, region.as_ref().is_some_and(|region| region.has_indexed_powers()));
                        }
                        let opened = scope.open(node)?;
                        if !self.expand && sum.exposed_sum != Some((scope as *const _ as usize, node)) && self.tensor_alternatives(literal) {
                            // Reduce this existing scope independently; its
                            // alternatives remain inside one opaque factor of
                            // the surrounding product. Boundary bindings are
                            // applied by the same positional emitter below.
                            let result = scope.reductions[&node]
                                .get_or_init(|| self.contract_graph(opened, opened.root, None, sum.rank_one))
                                .as_ref()?;
                            if result.root.as_view() != literal || result.status != ReductionStatus::Complete {
                                scoped.insert((scope as *const _ as usize, node), result);
                            }
                            if let Some(observed) = &result.observations {
                                sum.scalar_forms.insert(result.root.as_view(), observed.needs_dot_normalization());
                            }
                            return Some(if scope.ports[&node].0.is_empty() {
                                TermLeaf::Value(InputLeaf::Scalar(result.root.as_view()))
                            } else {
                                let position = sum.register_factor(scope, node, positions)?;
                                sum.opaque_factors[position].0 = result.root.as_view();
                                TermLeaf::Value(InputLeaf::Subtree(position))
                            });
                        }
                        let (node, size) =
                            self.compile_scope(opened, opened.root, tape, sum, positions, scoped)?;
                        Some(TermLeaf::Reference(node, size))
                    } else if scope.ports[&node].0.is_empty() {
                        // Equal finished coefficients share expression identity
                        // across branches, independently of graph occurrence.
                        Some(TermLeaf::Value(InputLeaf::Scalar(literal)))
                    } else {
                        let position = sum.register_factor(scope, node, positions)?;
                        Some(TermLeaf::Value(InputLeaf::Subtree(position)))
                    }
                } else if let AtomView::Num(number) = literal {
                    if let Ok(coefficient) = Rational::try_from(literal) {
                        Some(TermLeaf::Number(InputLeaf::Atom(literal), coefficient))
                    } else if matches!(
                        number.get_coeff_view(),
                        CoefficientView::Natural(..) | CoefficientView::Large(..)
                    ) {
                        Some(TermLeaf::Value(InputLeaf::Scalar(literal)))
                    } else {
                        None
                    }
                } else {
                    Some(TermLeaf::Value(InputLeaf::Atom(literal)))
                }
            },
        )
    }

    fn contract_graph(
        &self,
        scope: &ContractionScope<'_>,
        root: NodeIndex,
        order: Option<&[usize]>,
        rank_one: bool,
    ) -> Option<FactorizedContraction> {
        if !self.expand
            && let Some(source) = scope.literal(root)
            && self
                .observations
                .as_ref()
                .and_then(|observed| observed.region(source))
                .is_some_and(|region| {
                    region.excludes_contraction_sources(self.contraction_settings(rank_one))
                })
        {
            return Some(FactorizedContraction {
                root: source.to_owned(),
                status: ReductionStatus::Complete,
                observations: self
                    .observations
                    .as_ref()
                    .and_then(|observed| observed.certify_scalar_region(source)),
            });
        }
        if self.power_scope(scope, root).is_some() {
            return self.scoped_power(scope, root, rank_one).cloned();
        }
        let scalar_boundary = scope.ports[&root].0.is_empty();
        if !self.expand
            && matches!(
                scope.network.graph.graph[root],
                NetworkNode::Op(NetworkOp::Sum)
            )
        {
            // Existing additive branches already own their dummy scopes. Visit
            // them independently instead of exposing them to a product tape.
            let mut terms = Vec::new();
            let mut status = ReductionStatus::Complete;
            let mut scalar_form = Some(false);
            for node in scope.tree.iter_children(root, &scope.network.graph.graph) {
                let result = self.contract_graph(scope, node, None, rank_one)?;
                status = status.max(result.status);
                scalar_form = scalar_form
                    .zip(result.observations.as_ref())
                    .map(|(needs_dots, observed)| needs_dots || observed.needs_dot_normalization());
                terms.push(result.root);
            }
            return Some(FactorizedContraction {
                root: Atom::add_many(terms),
                status,
                observations: scalar_form
                    .filter(|_| scalar_boundary)
                    .map(DomainObservations::closed_scalar),
            });
        }
        let nodes = scope.factors(root);
        let selected = nodes
            .iter()
            .enumerate()
            .filter_map(|(i, &node)| {
                scope
                    .contains_work(node, self, self.contraction_settings(rank_one))
                    .then_some(i)
            })
            .collect();
        let mut selected = self.factor_order(scope, &nodes, selected, order)?;
        if !self.expand && order.is_none() {
            // Apply outer atomic sources before visiting opaque arithmetic
            // scopes. A scalar or tensor sum is never a list of alternatives
            // in the surrounding contraction state.
            selected.sort_by_key(|&position| {
                matches!(
                    scope.literal(nodes[position]),
                    Some(AtomView::Add(_) | AtomView::Mul(_) | AtomView::Pow(_))
                )
            });
        }
        let components = scope.network.graph.slot_components(&scope.tree, &nodes);
        let mut exposed_sums = AHashSet::new();
        if !self.expand {
            for component in &components {
                let alternatives = component
                    .iter()
                    .copied()
                    .filter(|&node| {
                        !scope.ports[&node].0.is_empty()
                            && scope
                                .literal(node)
                                .is_some_and(|source| self.tensor_alternatives(source))
                    })
                    .collect::<Vec<_>>();
                if let [node] = alternatives.as_slice() {
                    exposed_sums.insert(*node);
                }
            }
        }
        let ranks = selected
            .into_iter()
            .enumerate()
            .map(|(rank, position)| (nodes[position], rank))
            .collect::<AHashMap<_, _>>();
        let mut roots = Vec::with_capacity(components.len());
        let mut status = ReductionStatus::Complete;
        let mut scalar_form = Some(false);
        let mut changed = false;
        let mut refused = Vec::new();
        for component in components {
            let mut selected = component
                .iter()
                .enumerate()
                .filter_map(|(position, node)| ranks.get(node).map(|&rank| (position, rank)))
                .collect::<Vec<_>>();
            selected.sort_unstable_by_key(|&(_, rank)| rank);
            let selected = selected
                .into_iter()
                .map(|(position, _)| position)
                .collect::<Vec<_>>();
            let result =
                self.contract_component(scope, &component, &selected, &exposed_sums, rank_one);
            let (result, component_changed) = match result {
                Some(result) => result,
                None => {
                    // A refusal can follow tentative port substitutions or
                    // coefficient emission. Retain the exact source, discarding
                    // this component's scratch and tentative observations.
                    let root = Atom::mul_many(
                        component
                            .iter()
                            .map(|&node| scope.literal(node).map(|value| value.to_owned()))
                            .collect::<Option<Vec<_>>>()?,
                    );
                    refused.push(roots.len());
                    (
                        FactorizedContraction {
                            root,
                            status: ReductionStatus::Deferred,
                            observations: None,
                        },
                        false,
                    )
                }
            };
            changed |= component_changed;
            status = status.max(result.status);
            scalar_form = scalar_form
                .zip(result.observations.as_ref())
                .map(|(needs_dots, observed)| needs_dots || observed.needs_dot_normalization());
            roots.push(result.root);
        }
        if !refused.is_empty() {
            if !changed {
                // Preserve the checked whole-source fallback when factorized
                // planning made no independent progress. In particular, an
                // inexact coefficient must not become an admitted polynomial.
                return None;
            }
            for position in refused {
                // Independent work is already exact. Rewrite only the refused
                // components, with the same dummy owner as the whole domain.
                roots[position] = self.contract_ordered(
                    roots[position].as_view(),
                    rank_one,
                    self.observations
                        .as_ref()
                        .map(|observed| &observed.candidates),
                );
            }
        }
        let root = if status == ReductionStatus::Complete && !changed {
            let source = scope.literal(root)?;
            // The bracket normalizer owns whether a product can be exposed:
            // closed scalar wrappers disappear, while AUTO operand order stays.
            if matches!(source, AtomView::Fun(function)
                if function.get_symbol() == self.tags.bracket)
            {
                crate::shorthands::bracket::BracketNormalizer::normalize(source)
            } else {
                source.to_owned()
            }
        } else {
            Atom::mul_many(roots)
        };
        Some(FactorizedContraction {
            root,
            status,
            observations: scalar_form
                .filter(|_| status == ReductionStatus::Complete && scalar_boundary)
                .map(DomainObservations::closed_scalar),
        })
    }

    /// Execute one connected component on the admitted graph. Scratch bindings,
    /// compiled factors and tentative proofs cannot escape a refused component.
    fn contract_component<'a>(
        &self,
        scope: &'a ContractionScope<'_>,
        nodes: &[NodeIndex],
        selected: &[usize],
        exposed_sums: &AHashSet<NodeIndex>,
        rank_one: bool,
    ) -> Option<(FactorizedContraction, bool)> {
        let mut slots = SlotMatcher::default();
        let mut sum = ComponentSum::new(self, Intake::Contraction, &mut slots);
        sum.rank_one = rank_one;
        let mut positions = AHashMap::new();
        let mut remaining = Vec::with_capacity(nodes.len());
        for &node in nodes {
            let position = sum.register_factor(scope, node, &mut positions)?;
            remaining.push((
                position,
                sum.opaque_factors[position]
                    .1
                    .iter()
                    .copied()
                    .map(Argument::Original)
                    .collect::<Vec<_>>(),
            ));
        }
        if selected.is_empty() {
            // Untouched components need neither a polynomial frontier nor a
            // unit coefficient. Keep the admitted scalar facts and exact leaves.
            let mut scalar_form = Some(false);
            let mut factors = Vec::with_capacity(remaining.len());
            for (position, arguments) in remaining {
                scalar_form = scalar_form
                    .zip(sum.scalar_form(sum.opaque_factors[position].0))
                    .filter(|_| sum.opaque_factors[position].1.is_empty())
                    .map(|(first, second)| first || second);
                factors.push(sum.emit_factor(position, &arguments)?);
            }
            return Some((
                FactorizedContraction {
                    root: Atom::mul_many(factors),
                    status: ReductionStatus::Complete,
                    observations: scalar_form.map(DomainObservations::closed_scalar),
                },
                false,
            ));
        }
        let mut scoped = AHashMap::new();
        let mut compile = |sum: &mut ComponentSum<'a, '_>, position: usize| {
            let mut tape = std::mem::take(&mut sum.input);
            sum.exposed_sum = exposed_sums
                .contains(&nodes[position])
                .then_some((scope as *const _ as usize, nodes[position]));
            let result = self.compile_scope(
                scope,
                nodes[position],
                &mut tape,
                sum,
                &mut positions,
                &mut scoped,
            );
            sum.exposed_sum = None;
            sum.input = tape;
            result.map(|(node, _)| node)
        };
        let mut result = sum.contract_states(
            Remaining {
                factors: remaining,
                residual: Vec::new(),
            },
            selected,
            &mut compile,
        )?;
        for replacement in scoped.values() {
            result.status = result.status.max(replacement.status);
        }
        Some((result, sum.contracted || !scoped.is_empty()))
    }

    /// The ordered rewrite is shared by refused intrinsic components and the
    /// checked callback path. Keep existing additive branches independently scoped.
    pub(crate) fn contract_ordered<const N: usize>(
        &self,
        source: AtomView<'_>,
        rank_one: bool,
        observed: Option<&SimplificationCandidates<N>>,
    ) -> Atom {
        use crate::shorthands::schoonschip::{SchoonschipSettings, SchoonschipWithSettings};
        let mut settings = SchoonschipSettings::default().with_chain_like_functions();
        settings.schoonschip_rank1_tensors = rank_one;
        settings.metrics = self.metrics;
        settings.representations = self.representations.clone();
        let rewrite = SchoonschipWithSettings {
            settings: &settings,
        };
        match source {
            AtomView::Add(sum) => Atom::add_many(
                sum.iter()
                    .map(|term| rewrite.run_observed(term, observed, self)),
            ),
            _ => rewrite.run_observed(source, observed, self),
        }
    }
}

impl<'a> ComponentSum<'a, '_> {
    /// Reuse the admitted inventory; never inspect a completed coefficient to
    /// rediscover whether its source contains an identity or an index.
    fn scalar_form(&self, value: AtomView<'a>) -> Option<bool> {
        if let Some(&needs_dots) = self.scalar_forms.get(&value) {
            return Some(needs_dots);
        }
        if matches!(value, AtomView::Num(_)) {
            return Some(false);
        }
        if let AtomView::Var(variable) = value {
            let head = variable.get_symbol();
            let tags = self.contractor.tags;
            if head.get_wildcard_level() == 0
                && !head.has_tag(&tags.tensor)
                && !head.has_tag(&tags.rank1)
                && !head.has_tag(&tags.representation)
                && head != tags.chain_in
                && head != tags.chain_out
            {
                return Some(false);
            }
        }
        let observed = self.contractor.observations.as_ref();
        let known = observed.is_some_and(|observed| observed.region(value).is_some());
        if known {
            observed
                .and_then(|observed| observed.certify_scalar_region(value))
                .map(|observed| observed.needs_dot_normalization())
        } else {
            // The shallow parser can compose scalar coefficients that were
            // separate factors at admission. Combine only these new arithmetic
            // connectors, stopping at every already-observed subtree.
            match value {
                AtomView::Add(sum) => sum.iter().try_fold(false, |needs_dots, term| {
                    self.scalar_form(term).map(|next| needs_dots || next)
                }),
                AtomView::Mul(product) => product.iter().try_fold(false, |needs_dots, factor| {
                    self.scalar_form(factor).map(|next| needs_dots || next)
                }),
                AtomView::Pow(power) => {
                    let (base, exponent) = power.get_base_exp();
                    self.scalar_form(base)
                        .zip(self.scalar_form(exponent))
                        .map(|(base, exponent)| base || exponent)
                }
                _ => None,
            }
        }
    }

    fn variable_scalar_form(&mut self, variable: &Variable<'a>) -> Option<bool> {
        match variable {
            Variable::Scalar(value @ AtomView::Fun(function))
                if self.scalar_form(*value).is_none() =>
            {
                // A compact dot can come from the shallow parser's scalar
                // store, after notation normalization changed its source key.
                // Reuse the same admitted vector/space parser as dot emission.
                let (space, first, second) = self.compact_dot(*function)?;
                let needs_dots =
                    self.variable_scalar_form(&Variable::Dot(space, [first, second]))?;
                self.scalar_forms.insert(*value, needs_dots);
                Some(needs_dots)
            }
            Variable::Scalar(value) => self.scalar_form(*value),
            Variable::Subtree(position) if self.opaque_factors[*position].1.is_empty() => {
                self.scalar_form(self.opaque_factors[*position].0)
            }
            Variable::Dot(_, vectors) => {
                // Vector admission already proved intrinsic heads, compatible
                // atomic dimensions and scalar metadata. Require the stronger
                // identity-free certificate before skipping future family work.
                for &vector in vectors {
                    let source = self.vectors[vector];
                    for argument in source.iter().take(source.get_nargs() - 1) {
                        self.scalar_form(argument)?;
                    }
                }
                Some(true)
            }
            _ => None,
        }
    }

    fn register_factor(
        &mut self,
        scope: &'a ContractionScope<'_>,
        node: NodeIndex,
        positions: &mut AHashMap<(usize, NodeIndex), usize>,
    ) -> Option<usize> {
        let key = (scope as *const _ as usize, node);
        if let Some(&position) = positions.get(&key) {
            return Some(position);
        }
        let position = self.opaque_factors.len();
        let (ports, interface) = &scope.ports[&node];
        self.scalar_interfaces
            .insert(scope.literal(node)?, ports.is_empty());
        self.opaque_factors.push((
            scope.literal(node)?,
            ports.iter().map(Atom::as_view).collect(),
            interface.clone(),
        ));
        self.factor_roots.push(None);
        positions.insert(key, position);
        Some(position)
    }

    fn opaque_factor(&mut self, position: usize, arguments: &[Argument<'a>]) -> Option<()> {
        let tensor = self.tensors.len();
        self.tensors
            .push((TensorSource::Factor(position), arguments.to_vec()));
        for (port, &argument) in arguments.iter().enumerate() {
            if let Argument::Original(slot) = argument {
                if !self.permits_slot(slot) {
                    continue;
                }
                let (space, index) = self.resolve_endpoint(slot)?;
                let node = self.endpoint(slot, space, index)?;
                self.nodes[node]
                    .terminals
                    .push(Terminal::Tensor(tensor, port));
            }
        }
        Some(())
    }

    fn residual_factor(&mut self, variable: usize, exponent: u32) -> Option<()> {
        match self.variables[variable].clone() {
            Variable::Vector(head, slot) => {
                let (space, index) = self.resolve_endpoint(slot)?;
                if exponent > 1 {
                    self.dot(space, head, head, exponent / 2);
                }
                if exponent % 2 == 1 {
                    let node = self.endpoint(slot, space, index)?;
                    self.nodes[node].terminals.push(Terminal::Vector(head));
                }
            }
            Variable::Metric([first, second]) => {
                self.metric(
                    Argument::Original(first),
                    Argument::Original(second),
                    exponent,
                )?;
            }
            Variable::Tensor(source, arguments) => {
                if exponent != 1 {
                    return None;
                }
                self.tensor_with_arguments(source, arguments)?;
            }
            value => self.variable(value, exponent),
        }
        Some(())
    }

    fn reset_components(&mut self, coefficient: &Rational) {
        self.nodes.clear();
        self.occurrences.clear();
        self.factors.clear();
        self.tensors.clear();
        self.metrics = 0;
        self.coefficient = coefficient.clone();
    }

    fn factor_terms(&mut self, position: usize) -> Option<Vec<FactorTerm>> {
        let node = self.factor_roots.get(position).copied().flatten()?;
        self.input.clear_terms();
        self.input
            .distribute(&mut vec![(node, 1)], &Rational::one())?;
        Some(self.input.take_terms())
    }

    fn contract_states(
        &mut self,
        remaining: Remaining<'a>,
        order: &[usize],
        compile: &mut impl FnMut(&mut Self, usize) -> Option<usize>,
    ) -> Option<FactorizedContraction> {
        #[cfg(feature = "reference-cases")]
        let _states_timing = crate::reference_cases::timing::scope(
            crate::reference_cases::timing::Phase::ContractionStates,
        );
        let mut states = vec![(remaining, ScalarPolynomial::new_one(&Q))];
        let mut emitted: AHashMap<usize, Atom> = AHashMap::new();
        #[cfg(feature = "reference-cases")]
        let mut emitted_bytes = 0usize;
        let mut scalar_variables: AHashMap<Atom, usize> = AHashMap::new();
        let mut scalar_form = Some(false);
        #[cfg(feature = "reference-cases")]
        let mut generated = 0usize;
        for &selected in order {
            if self.factor_roots[selected].is_none() {
                self.factor_roots[selected] = Some(compile(self, selected)?);
            }
            let terms = self.factor_terms(selected)?;
            // Measure generated coefficient storage without placing an
            // automatic work limit on the selected contraction.
            #[cfg(feature = "reference-cases")]
            let state_bytes = states
                .iter()
                .map(|(remaining, coefficient)| {
                    Self::coefficient_bytes(coefficient) + remaining.storage_bytes()
                })
                .fold(0usize, usize::saturating_add);
            #[cfg(feature = "reference-cases")]
            let count = states
                .iter()
                .map(|(_, weight)| weight.nterms())
                .fold(0usize, usize::saturating_add)
                .saturating_mul(terms.len());
            #[cfg(feature = "reference-cases")]
            let predicted = state_bytes.saturating_add(emitted_bytes.saturating_mul(2));
            #[cfg(feature = "reference-cases")]
            crate::reference_cases::timing::frontier(crate::reference_cases::timing::Frontier {
                selected,
                states: states.len(),
                local_terms: terms.len(),
                generated,
                predicted_bytes: predicted,
                coefficient_bytes: states
                    .iter()
                    .map(|(_, coefficient)| Self::coefficient_bytes(coefficient))
                    .sum(),
                max_coefficient_bytes: states
                    .iter()
                    .map(|(_, coefficient)| Self::coefficient_bytes(coefficient))
                    .max()
                    .unwrap_or(0),
            });
            #[cfg(feature = "reference-cases")]
            {
                generated = generated.saturating_add(count);
            }
            let mut positions: AHashMap<Remaining<'a>, usize> = AHashMap::new();
            let mut next: Vec<(Remaining<'a>, ScalarPolynomial)> = Vec::new();
            #[cfg(feature = "reference-cases")]
            let mut next_bytes = 0usize;
            for (remaining, weight) in &states {
                let source = remaining
                    .factors
                    .iter()
                    .find(|(position, _)| *position == selected)?;
                self.overrides = self.opaque_factors[selected]
                    .1
                    .iter()
                    .copied()
                    .zip(source.1.iter().copied())
                    .collect();
                for (factors, coefficient) in &terms {
                    self.reset_components(coefficient);
                    for &(factor, exponent) in factors {
                        let leaf = self.input.leaves()[factor];
                        match leaf {
                            InputLeaf::Scalar(value) => {
                                self.variable(Variable::Scalar(value), exponent);
                            }
                            InputLeaf::Atom(value) if self.scalar_spectator(value) => {
                                self.variable(Variable::Scalar(value), exponent);
                            }
                            InputLeaf::Atom(value) => self.factor(value, exponent)?,
                            InputLeaf::Subtree(position) => {
                                let ports = &self.opaque_factors[position].1;
                                if ports.is_empty() {
                                    self.variable(Variable::Subtree(position), exponent);
                                } else {
                                    if exponent != 1 {
                                        return None;
                                    }
                                    let arguments = ports
                                        .iter()
                                        .map(|&slot| self.argument(slot))
                                        .collect::<Vec<_>>();
                                    self.opaque_factor(position, &arguments)?;
                                }
                            }
                        }
                    }
                    for &(variable, exponent) in &remaining.residual {
                        self.residual_factor(variable, exponent)?;
                    }
                    for (position, arguments) in &remaining.factors {
                        if *position != selected {
                            self.opaque_factor(*position, arguments)?;
                        }
                    }
                    self.reduce_components()?;
                    let mut key = Remaining {
                        factors: Vec::new(),
                        residual: Vec::new(),
                    };
                    let mut variables = Vec::new();
                    let mut exponents: Vec<u32> = Vec::new();
                    for factor in 0..self.factors.len() {
                        let (variable, exponent) = self.factors[factor];
                        match self.variables[variable].clone() {
                            Variable::Tensor(TensorSource::Factor(position), arguments) => {
                                key.factors.push((position, arguments));
                            }
                            value @ (Variable::Scalar(_)
                            | Variable::Subtree(_)
                            | Variable::Dot(_, _)) => {
                                let atom = if let Some(value) = emitted.get(&variable) {
                                    value.clone()
                                } else {
                                    scalar_form = scalar_form
                                        .zip(self.variable_scalar_form(&value))
                                        .map(|(first, second)| first || second);
                                    let atom = self.emit_variable(&value)?;
                                    #[cfg(feature = "reference-cases")]
                                    {
                                        emitted_bytes = emitted_bytes
                                            .saturating_add(atom.as_view().get_byte_size());
                                    }
                                    emitted.insert(variable, atom.clone());
                                    atom
                                };
                                let index = *scalar_variables.entry(atom).or_insert(variable);
                                let variable = PolyVariable::Temporary(index);
                                if let Some(position) =
                                    variables.iter().position(|entry| entry == &variable)
                                {
                                    exponents[position] =
                                        exponents[position].checked_add(exponent)?;
                                } else {
                                    variables.push(variable);
                                    exponents.push(exponent);
                                }
                            }
                            _ => key.residual.push((variable, exponent)),
                        }
                    }
                    key.factors.sort_unstable_by_key(|(position, _)| *position);
                    key.residual.sort_unstable();
                    for (variable, exponent) in variables.iter().zip(&exponents) {
                        if let Some(position) = weight
                            .variables()
                            .iter()
                            .position(|entry| entry == variable)
                        {
                            // Machine exponent capacity is not a work budget.
                            // A decline preserves the caller's exact source.
                            weight.degree(position).checked_add(*exponent)?;
                        }
                    }
                    #[cfg(feature = "reference-cases")]
                    let incoming_bytes = weight.nterms().saturating_mul(
                        std::mem::size_of::<Rational>()
                            + (weight.nvars() + variables.len()) * std::mem::size_of::<u32>(),
                    );
                    let position = positions.get(&key).copied();
                    // The position map and frontier both own a key.
                    #[cfg(feature = "reference-cases")]
                    let key_bytes = if position.is_none() {
                        key.storage_bytes().saturating_mul(2)
                    } else {
                        0
                    };
                    #[cfg(feature = "reference-cases")]
                    let target_bytes = position.map_or(0, |position| {
                        let target = &next[position].1;
                        target.nterms().saturating_mul(
                            std::mem::size_of::<Rational>()
                                + (target.nvars() + weight.nvars() + variables.len())
                                    * std::mem::size_of::<u32>(),
                        )
                    });
                    // Existing polynomial addition owns the only workspace:
                    // two unified inputs and the merged output.
                    #[cfg(feature = "reference-cases")]
                    let live_bytes = state_bytes
                        .saturating_add(next_bytes)
                        .saturating_add(key_bytes)
                        .saturating_add(emitted_bytes.saturating_mul(2))
                        .saturating_add(
                            incoming_bytes
                                .saturating_add(target_bytes)
                                .saturating_mul(2),
                        );
                    #[cfg(feature = "reference-cases")]
                    crate::reference_cases::timing::frontier_workspace(live_bytes);
                    let monomial = ScalarPolynomial::new(&Q, Some(1), Arc::new(variables))
                        .monomial(self.coefficient.clone(), exponents);
                    let coefficient = weight * &monomial;
                    if coefficient.is_zero() {
                        continue;
                    }
                    #[cfg(feature = "reference-cases")]
                    let _merge_timing = crate::reference_cases::timing::scope(
                        crate::reference_cases::timing::Phase::CoefficientMerge,
                    );
                    if let Some(position) = position {
                        let old = std::mem::replace(
                            &mut next[position].1,
                            ScalarPolynomial::new_zero(&Q),
                        );
                        #[cfg(feature = "reference-cases")]
                        {
                            next_bytes -= Self::coefficient_bytes(&old);
                        }
                        let merged = old + coefficient;
                        #[cfg(feature = "reference-cases")]
                        {
                            next_bytes += Self::coefficient_bytes(&merged);
                        }
                        next[position].1 = merged;
                    } else {
                        #[cfg(feature = "reference-cases")]
                        {
                            next_bytes += Self::coefficient_bytes(&coefficient) + key_bytes;
                        }
                        positions.insert(key.clone(), next.len());
                        next.push((key, coefficient));
                    }
                }
            }
            self.overrides.clear();
            states = next
                .into_iter()
                .filter(|(_, weight)| !weight.is_zero())
                .collect();
        }
        self.overrides.clear();
        #[cfg(feature = "reference-cases")]
        let _emission_timing = crate::reference_cases::timing::scope(
            crate::reference_cases::timing::Phase::ContractionEmission,
        );
        let mut terms = Vec::with_capacity(states.len());
        for (remaining, weight) in states {
            let mut factors = vec![self.emit_coefficient(&weight)?];
            for (position, arguments) in remaining.factors {
                scalar_form = scalar_form
                    .zip(self.scalar_form(self.opaque_factors[position].0))
                    .filter(|_| self.opaque_factors[position].1.is_empty())
                    .map(|(first, second)| first || second);
                factors.push(self.emit_factor(position, &arguments)?);
            }
            for (variable, exponent) in remaining.residual {
                let atom = if let Some(value) = emitted.get(&variable) {
                    value.clone()
                } else {
                    scalar_form = scalar_form
                        .zip(self.variable_scalar_form(&self.variables[variable].clone()))
                        .map(|(first, second)| first || second);
                    let atom = self.emit_variable(&self.variables[variable])?;
                    emitted.insert(variable, atom.clone());
                    atom
                };
                factors.push(atom.pow(exponent));
            }
            terms.push(Atom::mul_many(factors));
        }
        Some(FactorizedContraction {
            root: Atom::add_many(terms),
            status: ReductionStatus::Complete,
            observations: scalar_form.map(DomainObservations::closed_scalar),
        })
    }

    pub(super) fn argument(&self, value: AtomView<'a>) -> Argument<'a> {
        self.overrides
            .iter()
            .find_map(|&(source, target)| (source == value).then_some(target))
            .unwrap_or(Argument::Original(value))
    }

    fn emit_argument(&self, argument: Argument<'a>) -> Atom {
        match argument {
            Argument::Original(value) => value.to_owned(),
            Argument::Vector(space, head) => {
                let (representation, dimension) = self.spaces[space];
                self.emit_vector(head, representation.to_symbolic([dimension]).as_view())
            }
        }
    }

    pub(super) fn emit_factor(&self, position: usize, arguments: &[Argument<'a>]) -> Option<Atom> {
        let (source, ports, interface) = &self.opaque_factors[position];
        let source = source.to_owned();
        let replacements = ports
            .iter()
            .zip(arguments)
            .enumerate()
            .filter(|(_, (port, argument))| **argument != Argument::Original(**port))
            .map(|(position, (_, &argument))| (position, self.emit_argument(argument)))
            .collect::<std::collections::HashMap<_, _>>();
        if replacements.is_empty() {
            return Some(source);
        }
        let source = SymbolicTensor::from_normalized_parts(source, interface.clone());
        let expression = crate::tensor::composition::rewrite_interface_ports(
            &source,
            &replacements,
            Some(self.contractor.dummy_state()),
        )
        .ok()?;
        Some(
            crate::shorthands::schoonschip::DotNormalizer::with_settings(
                expression.as_view(),
                crate::tensor::ContractSettings {
                    metrics: self.contractor.metrics,
                    rank_one: self.rank_one,
                    representations: self.contractor.representations.as_deref(),
                    ..Default::default()
                },
            ),
        )
    }

    pub(super) fn metric(
        &mut self,
        first: Argument<'a>,
        second: Argument<'a>,
        exponent: u32,
    ) -> Option<()> {
        match (first, second) {
            (Argument::Vector(space, a), Argument::Vector(other, b)) => {
                if space != other {
                    return None;
                }
                self.dot(space, a, b, exponent);
            }
            (Argument::Vector(space, head), Argument::Original(slot))
            | (Argument::Original(slot), Argument::Vector(space, head)) => {
                let (other, index) = self.resolve_endpoint(slot)?;
                if space != other {
                    return None;
                }
                if exponent > 1 {
                    self.dot(space, head, head, exponent / 2);
                    self.contracted = true;
                }
                if exponent % 2 == 1 {
                    let node = self.endpoint(slot, space, index)?;
                    self.nodes[node].terminals.push(Terminal::Vector(head));
                }
            }
            (Argument::Original(first), Argument::Original(second)) => {
                let (space, first_index) = self.resolve_endpoint(first)?;
                let (other, second_index) = self.resolve_endpoint(second)?;
                if space != other {
                    return None;
                }
                if exponent > 1 {
                    if first_index == second_index {
                        return None;
                    }
                    self.factor(self.spaces[space].1, exponent / 2)?;
                    self.contracted = true;
                }
                if exponent % 2 == 1 {
                    let a = self.endpoint(first, space, first_index)?;
                    let b = self.endpoint(second, space, second_index)?;
                    let a = self.root(a);
                    let b = self.root(b);
                    self.metrics += 1;
                    self.nodes[a].metrics += 1;
                    self.nodes[b].parent = a;
                }
            }
        }
        Some(())
    }
}

#[cfg(test)]
mod tests {
    use super::super::tests::{input, setup};
    use super::*;
    use crate::shorthands::schoonschip::Schoonschip;

    #[test]
    fn source_free_color_scopes_stay_opaque_beside_real_contractions() {
        let reps = crate::test_support::test_initialize();
        let [a, b, c, d, e, f] = [97801, 97802, 97803, 97804, 97805, 97806]
            .map(|index| reps.coad_da.to_symbolic([Atom::num(index)]));
        let rep = reps.cof_nc.to_symbolic([] as [Atom; 0]);
        let symmetric = |x: &Atom| {
            spenso::trace_sym!(
                &rep,
                crate::color_t!(x),
                crate::color_t!(&c),
                crate::color_t!(&d)
            )
        };
        let invariant = crate::color_f!(&a, &b, &e)
            * symmetric(&e)
            * crate::color_f!(&a, &b, &f)
            * symmetric(&f);
        let spectator = Atom::mul_many((0..32).map(|index| {
            let coefficient =
                Atom::var(symbolica::symbol!(&format!("opaque_color_scope::x{index}")));
            (coefficient + &invariant).pow(2)
        }));
        let mu = spenso::mink!(4, 97807);
        let nu = spenso::mink!(4, 97808);
        let source =
            SymbolicTensor::infer(&spectator * spenso::g!(&mu, &nu) * spenso::p!(&nu)).unwrap();
        let contractor = SlotContraction::new()
            .expanding(false)
            .observed(source.reduction_observations());
        crate::tensor::SELECTED_OPENS.with(|count| count.set(0));
        let result = contractor
            .contract_factorized(source.expression.as_view(), None, true)
            .unwrap();
        assert_eq!(result.root, spectator * spenso::p!(&mu));
        assert_eq!(result.status, ReductionStatus::Complete);
        assert!(
            result.observations.is_none(),
            "remaining colour identities are not terminal scalar algebra"
        );
        assert_eq!(crate::tensor::SELECTED_OPENS.with(|count| count.get()), 0);
    }

    #[test]
    fn source_free_open_color_sum_keeps_positional_metric_bindings() {
        let reps = crate::test_support::test_initialize();
        let [a, b, c, d, e, z] = [97811, 97812, 97813, 97814, 97815, 97816]
            .map(|index| reps.coad_da.to_symbolic([Atom::num(index)]));
        let rep = reps.cof_nc.to_symbolic([] as [Atom; 0]);
        let sum = |port: &Atom| {
            crate::color_f!(port, &b, &e)
                * spenso::trace_sym!(
                    &rep,
                    crate::color_t!(&e),
                    crate::color_t!(&c),
                    crate::color_t!(&d)
                )
                + crate::color_f!(port, &c, &e)
                    * spenso::trace_sym!(
                        &rep,
                        crate::color_t!(&e),
                        crate::color_t!(&b),
                        crate::color_t!(&d)
                    )
        };
        let source = SymbolicTensor::infer(spenso::g!(&a, &z) * sum(&a)).unwrap();
        let contractor = SlotContraction::new()
            .expanding(false)
            .observed(source.reduction_observations());
        crate::tensor::SELECTED_OPENS.with(|count| count.set(0));
        let result = contractor
            .contract_factorized(source.expression.as_view(), None, true)
            .unwrap();
        assert_eq!(result.root, sum(&z));
        assert_eq!(result.status, ReductionStatus::Complete);
        assert!(result.observations.is_none());
        assert_eq!(crate::tensor::SELECTED_OPENS.with(|count| count.get()), 0);
        assert_eq!(
            SymbolicTensor::infer(result.root).unwrap().structure,
            source.structure
        );
    }

    #[test]
    fn source_free_callback_region_still_requires_checked_contraction() {
        crate::test_support::test_initialize();
        let a = spenso::mink!(4, 97821);
        let b = spenso::mink!(4, 97822);
        let target = b.clone();
        let callback = spenso::tensor_symbol!(
            "source_proof_callback",
            norm = move |value, output| {
                if let AtomView::Fun(function) = value
                    && function.iter().any(|argument| argument == target.as_view())
                {
                    **output = Atom::one();
                }
            }
        );
        let sibling = spenso::tensor_symbol!("source_proof_callback_sibling");
        let source = SymbolicTensor::infer(
            spenso::g!(&a, &b)
                * (symbolica::function!(callback, &a) + symbolica::function!(sibling, &a)),
        )
        .unwrap();
        let contractor = SlotContraction::new().observed(source.reduction_observations());
        assert!(
            contractor
                .contract_factorized(source.expression.as_view(), None, true)
                .is_none()
        );
        assert!(source.contract(Default::default()).is_err());
    }

    #[test]
    fn negative_power_contracts_its_existing_graph_scope() {
        let contractor = setup();
        for (source, expected) in [
            (
                "spenso::bracket(p(spenso::mink(4,a))*q(spenso::mink(4,a)))^-1",
                "spenso::g(p(spenso::mink(4)),q(spenso::mink(4)))^-1",
            ),
            (
                "x+spenso::bracket(p(spenso::mink(4,a))*q(spenso::mink(4,a)))^-2",
                "x+spenso::g(p(spenso::mink(4)),q(spenso::mink(4)))^-2",
            ),
        ] {
            let source = input(source);
            let result = contractor
                .contract_factorized(source.as_view(), None, true)
                .unwrap();
            assert_eq!(result.root, input(expected), "{source}");
        }
    }

    #[test]
    fn factorized_states_retain_exact_complex_coefficients() {
        let contractor = setup();
        let source = input("1i*p(mink(4,a))*q(mink(4,a))+2i*r(mink(4,b))*s(mink(4,b))");
        let result = contractor
            .contract_factorized(source.as_view(), None, true)
            .unwrap();
        assert!(result.status == ReductionStatus::Complete);
        assert_eq!(
            result.root,
            input("1i*g(p(mink(4)),q(mink(4)))+2i*g(r(mink(4)),s(mink(4)))")
        );
        let inexact = input("1.25*p(mink(4,a))*q(mink(4,a))+2.5*r(mink(4,b))*s(mink(4,b))");
        assert!(
            contractor
                .contract_factorized(inexact.as_view(), None, true)
                .is_none()
        );
    }

    #[test]
    fn factorized_states_use_the_component_reducer_in_both_orders() {
        let contractor = setup();
        for source in [
            "g(mink(4,a),mink(4,b))*(p(mink(4,a))+q(mink(4,a)))*r(mink(4,b))",
            "(g(mink(4,a),mink(4,b))*p(mink(4,c))+g(mink(4,a),mink(4,c))*q(mink(4,b)))*(p(mink(4,a))*q(mink(4,b))*r(mink(4,c))+q(mink(4,a))*r(mink(4,b))*s(mink(4,c)))",
            "(g(mink(4,a),mink(4,b))+p(mink(4,a))*q(mink(4,b)))*(g(mink(4,a),mink(4,b))+r(mink(4,a))*s(mink(4,b)))",
            "(g(mink(4,a),mink(4,b))+t(mink(4,a),mink(4,b)))*p(mink(4,a))*q(mink(4,b))",
        ] {
            let source = input(source);
            let expected = source.expand().schoonschip().expand();
            let count = match source.as_view() {
                AtomView::Mul(product) => product.iter().len(),
                _ => 1,
            };
            let reversed = (0..count).rev().collect::<Vec<_>>();
            for order in [None, Some(reversed.as_slice())] {
                let contracted = contractor
                    .contract_factorized(source.as_view(), order, true)
                    .unwrap();
                SymbolicTensor::infer(contracted.root.clone()).unwrap();
                assert_eq!(contracted.root.expand(), expected, "{source}");
            }
        }
    }

    #[test]
    fn factorized_states_keep_foreign_and_disconnected_sums() {
        let contractor = setup();
        let spectator = input("(1+x)*(t(mink(4,c))+routing(z)*u(mink(4,c)))");
        let first = input("(p(mink(4,a))+q(mink(4,a)))*r(mink(4,a))");
        let second = input("(p(mink(4,b))+s(mink(4,b)))*q(mink(4,b))");
        let source = Atom::mul_many([&spectator, &first, &second]);
        let contracted = contractor
            .contract_factorized(source.as_view(), None, true)
            .unwrap();
        let result = contracted.root;
        let expected = Atom::mul_many([
            spectator,
            first.expand().schoonschip(),
            second.expand().schoonschip(),
        ]);
        assert_eq!(result, expected);
    }

    #[test]
    fn factorized_port_bindings_remain_occurrence_local() {
        let contractor = setup();
        let source =
            input("(p(mink(4,a))+q(mink(4,a)))*(p(mink(4,b))+q(mink(4,b)))*t(mink(4,a),mink(4,b))");
        let contracted = contractor
            .contract_factorized(source.as_view(), None, true)
            .unwrap();
        assert_eq!(
            contracted.root.expand(),
            source.expand().schoonschip().expand()
        );
        assert!(
            contractor
                .contract_factorized(source.as_view(), Some(&[0, 0, 1]), true)
                .is_none()
        );
    }

    #[test]
    fn polynomial_exponents_beyond_u16_remain_exact_and_complete() {
        let contractor = setup();
        let x = input("x");
        let spectator = input("(y+z)^17");
        let first = &spectator * x.pow(u16::MAX);
        let second = x.clone();
        let mut slots = SlotMatcher::default();
        let mut sum = ComponentSum::new(&contractor, Intake::Contraction, &mut slots);
        let leaf = sum.input.leaf(InputLeaf::Scalar(x.as_view())).unwrap();
        let power = sum.input.power(leaf, u32::from(u16::MAX)).unwrap();
        let scalar = sum
            .input
            .leaf(InputLeaf::Scalar(spectator.as_view()))
            .unwrap();
        let weighted = sum.input.group([power, scalar].into_iter(), false).unwrap();
        sum.opaque_factors = vec![
            (
                first.as_view(),
                Vec::new(),
                PartialStructure::from_logical_slots([]),
            ),
            (
                second.as_view(),
                Vec::new(),
                PartialStructure::from_logical_slots([]),
            ),
        ];
        sum.factor_roots = vec![Some(weighted.0), Some(leaf.0)];
        let result = sum
            .contract_states(
                Remaining {
                    factors: vec![(0, Vec::new()), (1, Vec::new())],
                    residual: Vec::new(),
                },
                &[0, 1],
                &mut |_, _| unreachable!("the two selected factors are already compiled"),
            )
            .unwrap();
        assert_eq!(result.status, ReductionStatus::Complete);
        assert_eq!(result.root, &spectator * x.pow(u32::from(u16::MAX) + 1));
        assert!(result.observations.is_some());
        assert!(sum.overrides.is_empty());
    }

    #[test]
    fn skewed_frontier_counts_each_opaque_coefficient_once() {
        let contractor = setup();
        let scalar = Atom::add_many(
            (0..512).map(|i| Atom::var(symbolica::symbol!(&format!("skewed_frontier::s{i}")))),
        );
        let source = (&scalar * input("p(mink(4,a))") + input("q(mink(4,a))"))
            * input(
                "(p(mink(4,b))+q(mink(4,b)))*(p(mink(4,c))+q(mink(4,c)))*t(mink(4,a),mink(4,b),mink(4,c))",
            );
        let result = contractor
            .contract_factorized(source.as_view(), None, true)
            .unwrap();
        assert_eq!(result.status, ReductionStatus::Complete);
        let mut expected = Vec::new();
        for a in ["p", "q"] {
            for b in ["p", "q"] {
                for c in ["p", "q"] {
                    let tensor = input(&format!("t({a}(mink(4)),{b}(mink(4)),{c}(mink(4)))"));
                    expected.push(if a == "p" { &scalar * tensor } else { tensor });
                }
            }
        }
        assert_eq!(result.root, Atom::add_many(expected));
    }

    #[test]
    fn factorized_frontier_completes_all_selected_contractions() {
        let contractor = setup();
        let slots = (0..9).map(|i| format!("mink(4,a{i})")).collect::<Vec<_>>();
        let source = input(&format!(
            "t({})*{}",
            slots.join(","),
            slots
                .iter()
                .map(|slot| format!("(x*p({slot})+y*q({slot}))"))
                .collect::<Vec<_>>()
                .join("*")
        ));
        let contracted = contractor
            .contract_factorized(source.as_view(), None, true)
            .unwrap();
        assert!(!contracted.root.is_zero());
        assert_eq!(contracted.status, ReductionStatus::Complete);
        use spenso::network::parsing::AtomStructureExt;
        assert!(
            !contracted.root.has_repeated_explicit_indices(),
            "all selected connections must be contracted"
        );
        let result = contracted.root;
        assert_eq!(
            result.expand().schoonschip().expand(),
            source.expand().schoonschip().expand()
        );
    }
    #[test]
    fn scalar_power_without_graph_incidence_keeps_its_original_source() {
        crate::test_support::test_initialize();
        let x = Atom::var(symbolica::symbol!("no_work_power_x"));
        let y = Atom::var(symbolica::symbol!("no_work_power_y"));
        let z = Atom::var(symbolica::symbol!("no_work_power_z"));
        let source = spenso::bracket!(x + y, z).pow(3);
        let result = SlotContraction::new()
            .contract_factorized(source.as_view(), None, true)
            .unwrap();
        assert!(result.status == ReductionStatus::Complete);
        assert_eq!(result.root, source);
    }

    #[test]
    fn scalar_rational_sums_keep_their_factorization() {
        let contractor = setup();
        for scalar in ["z", "q(0,spenso::cind(0))", "g(p(mink(4)),q(mink(4)))"] {
            let source = input(&format!(
                "n*(1/({scalar}+x)+1/({scalar}+y))+m*(1/({scalar}-x)+1/({scalar}-y))"
            ));
            for rank_one in [false, true] {
                let result = contractor
                    .contract_factorized(source.as_view(), None, rank_one)
                    .unwrap();
                assert!(result.status == ReductionStatus::Complete);
                // An exact Atom comparison protects the factored scalar input;
                // expanding both sides would hide this regression.
                assert_eq!(result.root, source, "{scalar}, rank_one={rank_one}");

                let tensor = SymbolicTensor::infer(source.clone()).unwrap();
                let settings = crate::tensor::ContractSettings::default();
                let settings = if rank_one {
                    settings
                } else {
                    settings.without_rank_one_tensors()
                };
                for settings in [settings, settings.with_order(&[0])] {
                    let result = tensor.contract(settings).unwrap();
                    assert_eq!(result, tensor);
                    assert_eq!(result.expression, source);
                }
            }
        }
    }
}

#[cfg(test)]
mod lazy_tests {
    use super::*;
    use spenso::network::tags::SPENSO_TAG;

    #[test]
    fn partial_arithmetic_scopes_have_only_one_level() {
        let _contractor = super::super::tests::setup();
        for (source, nodes) in [
            ("p(spenso::mink(4,a))*q(spenso::mink(4,a))+x*y", 3),
            ("(p(spenso::mink(4,a))*q(spenso::mink(4,a))+x)^2", 2),
            (
                "spenso::bracket(p(spenso::mink(4,a))*q(spenso::mink(4,a)))^-1",
                2,
            ),
        ] {
            let source = super::super::tests::input(source);
            let scope = ContractionScope::parse(source.as_view()).unwrap();
            assert_eq!(scope.ports.len(), nodes, "{source}");
            assert_eq!(scope.network.graph.graph.n_nodes(), nodes, "{source}");
            assert!(scope.opened.values().all(|entry| entry.get().is_none()));
        }
    }

    #[test]
    fn selected_arithmetic_reuses_admission_without_nested_graphs() {
        let setup = super::super::tests::setup();
        let input = super::super::tests::input;
        let scalar = Atom::add_many(
            (0..512).map(|i| Atom::var(symbolica::symbol!(&format!("selected_tape::s{i}")))),
        );
        let source = (&scalar * input("p(mink(4,a))") + input("q(mink(4,a))"))
            * input("t(mink(4,a),mink(4,b))");
        let tensor = SymbolicTensor::infer(source.clone()).unwrap();
        let contractor = setup.observed(tensor.reduction_observations());
        let scope = ContractionScope::parse(source.as_view()).unwrap();
        let result = contractor
            .contract_graph(&scope, scope.root, None, true)
            .unwrap();
        assert_eq!(result.status, ReductionStatus::Complete);
        assert_eq!(
            result.root,
            &scalar * input("t(p(mink(4)),mink(4,b))") + input("t(q(mink(4)),mink(4,b))")
        );
        assert!(
            scope.opened.values().all(|entry| entry.get().is_none()),
            "admitted arithmetic needs a tape traversal, not new port graphs"
        );
    }

    #[test]
    fn direct_selected_arithmetic_keeps_independent_dummy_scopes() {
        let setup = super::super::tests::setup();
        let input = super::super::tests::input;
        let source = input(
            "(p(mink(4,a))*q(mink(4,a))*r(mink(4,b))+s(mink(4,b)))*(p(mink(4,a))*q(mink(4,a))*t(mink(4,b),mink(4,c))+u(mink(4,b),mink(4,c)))",
        );
        let tensor = SymbolicTensor::infer(source.clone()).unwrap();
        let contractor = setup.observed(tensor.reduction_observations());
        let scope = ContractionScope::parse(source.as_view()).unwrap();
        let result = contractor
            .contract_graph(&scope, scope.root, None, true)
            .unwrap();
        let dot = input("g(p(mink(4)),q(mink(4)))");
        let expected = dot.pow(2) * input("t(r(mink(4)),mink(4,c))")
            + &dot * input("t(s(mink(4)),mink(4,c))")
            + &dot * input("u(r(mink(4)),mink(4,c))")
            + input("u(s(mink(4)),mink(4,c))");
        assert_eq!(result.status, ReductionStatus::Complete);
        assert_eq!(result.root.expand(), expected.expand());
        assert!(scope.opened.values().all(|entry| entry.get().is_none()));
    }

    #[test]
    fn scalar_powers_in_selected_sums_do_not_rebuild_graphs() {
        let setup = super::super::tests::setup();
        let input = super::super::tests::input;
        let source = input(
            "(m^4*p(mink(4,a))+m^2*g(q(mink(4)),s(mink(4)))^2*r(mink(4,a)))*t(mink(4,a),mink(4,b))",
        );
        let tensor = SymbolicTensor::infer(source.clone()).unwrap();
        let contractor = setup.observed(tensor.reduction_observations());
        let scope = ContractionScope::parse(source.as_view()).unwrap();
        let result = contractor
            .contract_graph(&scope, scope.root, None, true)
            .unwrap();
        assert_eq!(result.status, ReductionStatus::Complete);
        assert_eq!(
            result.root,
            input(
                "m^4*t(p(mink(4)),mink(4,b))+m^2*g(q(mink(4)),s(mink(4)))^2*t(r(mink(4)),mink(4,b))",
            )
        );
        assert!(scope.opened.values().all(|entry| entry.get().is_none()));
    }

    #[test]
    fn compact_tensor_sums_do_not_enter_the_scalar_coefficient_tape() {
        let setup = super::super::tests::setup();
        let input = super::super::tests::input;
        let source = input("p(mink(4))+q(mink(4))");
        let tensor = SymbolicTensor::infer(source.clone()).unwrap();
        assert_eq!(tensor.structure.logical_slots().len(), 1);
        let contractor = setup.observed(tensor.reduction_observations());
        let mut slots = SlotMatcher::default();
        let mut sum = ComponentSum::new(&contractor, Intake::Contraction, &mut slots);
        let mut tape = TermTape::default();
        contractor
            .compile_arithmetic(source.as_view(), &mut tape, &mut sum)
            .unwrap();
        // The rank-one alternatives keep their tensor ports. Zero written
        // indices does not establish the scalar interface required by a
        // coefficient, even before unresolved ports have been materialized.
        assert!(sum.scalar_interfaces.is_empty());
        assert_eq!(tape.leaves().len(), 2);
        assert!(
            tape.leaves()
                .iter()
                .all(|leaf| matches!(leaf, InputLeaf::Atom(_)))
        );
    }

    #[test]
    fn indexed_powers_in_selected_sums_keep_the_scoped_graph() {
        let setup = super::super::tests::setup();
        let input = super::super::tests::input;
        let source = input(
            "((p(mink(4,a))*q(mink(4,a))+x)^2*r(mink(4,b))+s(mink(4,b)))*t(mink(4,b),mink(4,c))",
        );
        let tensor = SymbolicTensor::infer(source.clone()).unwrap();
        let contractor = setup.observed(tensor.reduction_observations());
        let scope = ContractionScope::parse(source.as_view()).unwrap();
        let result = contractor
            .contract_graph(&scope, scope.root, None, true)
            .unwrap();
        assert_eq!(result.status, ReductionStatus::Complete);
        assert_eq!(
            result.root.expand(),
            input(
                "(g(p(mink(4)),q(mink(4)))+x)^2*t(r(mink(4)),mink(4,c))+t(s(mink(4)),mink(4,c))",
            )
            .expand()
        );
        assert!(scope.opened.values().any(|entry| entry.get().is_some()));
    }

    #[test]
    fn fixed_component_scalar_sums_stay_closed() {
        let contractor = super::super::tests::setup();
        let source =
            super::super::tests::input("n*(1/(q(0,spenso::cind(0))+x)+1/(q(0,spenso::cind(0))+y))");
        let scope = ContractionScope::parse(source.as_view()).unwrap();
        let result = contractor
            .contract_graph(&scope, scope.root, None, true)
            .unwrap();
        assert!(result.status == ReductionStatus::Complete);
        assert_eq!(result.root, source);
        assert!(scope.opened.values().all(|entry| entry.get().is_none()));
    }

    #[test]
    fn graph_order_starts_with_a_small_boundary_and_accepts_an_override() {
        let contractor = super::super::tests::setup();
        let source = super::super::tests::input(
            "spenso::g(spenso::mink(4,a),spenso::mink(4,b))*p(spenso::mink(4,a))*q(spenso::mink(4,b))",
        );
        let scope = ContractionScope::parse(source.as_view()).unwrap();
        let nodes = scope.factors(scope.root);
        let selected = (0..nodes.len()).collect::<Vec<_>>();
        let order = contractor
            .factor_order(&scope, &nodes, selected.clone(), None)
            .unwrap();
        assert_eq!(scope.ports[&nodes[order[0]]].0.len(), 1);
        let explicit = selected.iter().rev().copied().collect::<Vec<_>>();
        assert_eq!(
            contractor
                .factor_order(&scope, &nodes, selected, Some(&explicit))
                .unwrap(),
            explicit
        );
        let expected =
            super::super::tests::input("spenso::g(p(spenso::mink(4)),q(spenso::mink(4)))");
        for order in [None, Some(explicit.as_slice())] {
            let result = contractor
                .contract_graph(&scope, scope.root, order, true)
                .unwrap();
            assert!(result.status == ReductionStatus::Complete);
            assert_eq!(result.root, expected);
        }
    }

    #[test]
    fn unselected_thousand_term_factors_stay_closed() {
        crate::test_support::test_initialize();
        SPENSO_TAG.rank_one_tensor_symbol("lazy_contraction::p");
        SPENSO_TAG.tensor_symbol("lazy_contraction::t");
        SPENSO_TAG.tensor_symbol("lazy_contraction::u");
        let parse = |source: &str| {
            Atom::parse(
                source,
                "lazy_contraction",
                symbolica::parser::ParseSettings::symbolica(),
            )
            .unwrap()
        };
        for foreign in [false, true] {
            let first = Atom::add_many((0..1000).map(|i| {
                parse(&if foreign {
                    format!("t({i},spenso::mink(4,c))")
                } else {
                    format!("x({i})")
                })
            }));
            let second = Atom::add_many((0..1000).map(|i| {
                parse(&if foreign {
                    format!("u({i},spenso::mink(4,d))")
                } else {
                    format!("y({i})")
                })
            }));
            let core = parse("spenso::g(spenso::mink(4,a),spenso::mink(4,b))*p(spenso::mink(4,a))");
            let source = Atom::mul_many([core, first.clone(), second.clone()]);
            let scope = ContractionScope::parse(source.as_view()).unwrap();
            assert_eq!(scope.ports.len(), 5, "one product and four factor leaves");
            assert_eq!(scope.network.graph.graph.n_nodes(), 5);
            let protected = scope
                .factors(scope.root)
                .into_iter()
                .filter(|&node| {
                    let value = scope.literal(node).unwrap();
                    value == first.as_view() || value == second.as_view()
                })
                .collect::<Vec<_>>();
            assert_eq!(protected.len(), 2);
            let result = SlotContraction::new()
                .contract_graph(&scope, scope.root, None, true)
                .unwrap();
            assert!(result.status == ReductionStatus::Complete);
            assert_eq!(
                result.root,
                Atom::mul_many([parse("p(spenso::mink(4,b))"), first, second])
            );
            for node in protected {
                assert!(
                    scope.opened[&node].get().is_none(),
                    "unselected factor was opened"
                );
            }
            assert!(scope.opened.values().all(|opened| opened.get().is_none()));
        }
    }
}
