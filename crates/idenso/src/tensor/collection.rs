//! Selected algebra on the existing parsed graph and arithmetic tape.

use std::collections::hash_map::Entry;
use std::collections::{HashMap, HashSet};
use std::sync::{Arc, OnceLock};

use linnet::half_edge::NodeIndex;
use linnet::half_edge::involution::Hedge;
use linnet::half_edge::subgraph::{ModifySubSet, SuBitGraph};
use linnet::half_edge::tree::SimpleTraversalTree;
use linnet::tree::child_vec::ChildVecStore;
use spenso::network::TensorNetworkError;
use spenso::network::graph::{NetworkEdge, NetworkLeaf, NetworkNode, NetworkOp};
use spenso::network::parsing::{ParseSettings, ParseState};
use spenso::network::store::{TensorScalarStore, TensorScalarStoreMapping};
use spenso::network::tags::SPENSO_TAG;
use spenso::shadowing::{TensorCollectFilter, TermLeaf, TermTape};
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::partial::{
    PartialIndex, PartialSlot, PartialStructure, PartialStructureExt,
};
use spenso::structure::representation::{LibraryRep, LibrarySlot, Representation};
use spenso::structure::slot::{IsAbstractSlot, SlotMatcher};
use spenso::structure::{OrderedStructure, TensorStructure};
use symbolica::atom::{AddView, Atom, AtomCore, AtomOrView, AtomView};
use symbolica::domains::rational::Rational;
use symbolica::id::ConditionResult;

use super::composition::{self, ExplicitIndexOccurrences};
use super::inference::{InterfaceInference, LeafInference, TensorInferenceError};
use super::planning::{self, EnclosingLeaf, PlanningLeafSources, PlanningNetwork, PlanningSource};
use super::simplification::observation::DomainObservations;
use super::{SymbolicTensor, TensorProofs};

mod collected;
pub(crate) use collected::{CollectedOutput, CollectedPolynomial};

type Result<T> = std::result::Result<T, TensorInferenceError>;

#[cfg(test)]
thread_local! {
    pub(crate) static SHALLOW_BUILDS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    pub(crate) static SHALLOW_VIEWS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    pub(crate) static SELECTED_OPENS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    static PUBLISHED_TAPE_REGIONS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    static STREAM_THRESHOLD: std::cell::Cell<usize> = const { std::cell::Cell::new(STREAM_SUM_BYTES) };
    static PEAK_SUMMAND_NODES: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    static PEAK_OPEN_NODES: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    static ANALYZED_NODES: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

/// Whether the selected kernel needs monomials or can choose its own local
/// distribution while retaining arithmetic leaves. Both modes share the same
/// graph, copy scopes, tape and result publication.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum CollectionMode {
    Monomials,
    Factored,
}

type CollectionStates = Vec<(
    SymbolicTensor<PartialStructure>,
    SymbolicTensor<PartialStructure>,
)>;
type CollectionFinalizer<'a> =
    dyn FnMut(&mut CollectionStates, &[SymbolicTensor<PartialStructure>]) -> Result<()> + 'a;

enum Finalizer<'a, 'f> {
    Rows(&'a mut dyn FnMut(&mut CollectionStates) -> Result<()>),
    WithSpectators(&'a mut CollectionFinalizer<'f>),
}

impl Finalizer<'_, '_> {
    fn finish(
        &mut self,
        states: &mut CollectionStates,
        spectators: &[SymbolicTensor<PartialStructure>],
    ) -> Result<()> {
        match self {
            Self::Rows(finish) => finish(states),
            Self::WithSpectators(finish) => finish(states, spectators),
        }
    }
}

// Per-call intake state only. Graphs and arithmetic remain owned by Network and
// TermTape; the call retains its selected leaves and interfaces.
struct CollectionInput<Select> {
    observations: Arc<super::simplification::observation::DomainObservations>,
    state: ParseState<AbstractIndex>,
    select: Select,
    mode: CollectionMode,
    values: Vec<(SymbolicTensor<PartialStructure>, bool)>,
    positions: HashMap<(Atom, Vec<PartialSlot>, bool), usize>,
}

// Above this size the full sum graph can dominate the arithmetic tape.
const STREAM_SUM_BYTES: usize = 1 << 20;

fn stream_sum_bytes() -> usize {
    #[cfg(test)]
    {
        STREAM_THRESHOLD.get()
    }
    #[cfg(not(test))]
    {
        STREAM_SUM_BYTES
    }
}

type CollectionTree = SimpleTraversalTree<ChildVecStore<()>>;
type CollectionAnalysis = (
    HashMap<NodeIndex, bool>,
    HashMap<NodeIndex, PartialStructure>,
    HashMap<NodeIndex, Arc<ExplicitIndexOccurrences>>,
);

impl<Select: FnMut(AtomView<'_>) -> bool> CollectionInput<Select> {
    fn new(
        source: &SymbolicTensor<PartialStructure>,
        select: Select,
        mode: CollectionMode,
    ) -> Self {
        Self {
            observations: Arc::clone(source.reduction_observations()),
            state: SymbolicTensor::<PartialStructure>::reserved_dummies(std::iter::once(source)),
            select,
            mode,
            values: Vec::new(),
            positions: HashMap::new(),
        }
    }

    fn parse(
        &self,
        source: &SymbolicTensor<PartialStructure>,
    ) -> Result<Arc<PlanningNetwork<Atom>>> {
        #[cfg(test)]
        SHALLOW_VIEWS.with(|count| count.set(count.get() + 1));
        source.shallow_graph()
    }

    fn open<'src>(
        &self,
        source: AtomView<'src>,
        admission: PlanningSource<'src>,
    ) -> Result<PlanningNetwork<AtomOrView<'src>>> {
        #[cfg(test)]
        SELECTED_OPENS.with(|count| count.set(count.get() + 1));
        let opened = planning::parse(
            source,
            None,
            Some(admission),
            &ParseSettings::opaque_arithmetic(),
        )?;
        #[cfg(test)]
        PEAK_OPEN_NODES.set(PEAK_OPEN_NODES.get().max(opened.graph.graph.n_nodes()));
        Ok(opened)
    }

    fn selected(&mut self, value: AtomView<'_>) -> bool {
        if (self.select)(value) {
            return true;
        }
        if matches!(
            value,
            AtomView::Add(_) | AtomView::Mul(_) | AtomView::Pow(_)
        ) {
            if let Some(region) = self.observations.region(value)
                && let Some(selected) = region.any_leaf(&mut self.select)
            {
                return selected;
            }
            // Arbitrary user selection and newly exposed replacement regions
            // have no reusable family observation. Inspect only this leaf.
            let mut selected = false;
            value.visitor(&mut |node| {
                if selected {
                    return false;
                }
                if matches!(node, AtomView::Fun(_)) {
                    selected = (self.select)(node);
                    return false;
                }
                true
            });
            return selected;
        }
        false
    }

    fn emit<E: AtomCore>(
        &self,
        network: &PlanningNetwork<E>,
        tree: &CollectionTree,
        node: NodeIndex,
    ) -> Result<Atom> {
        network
            .graph
            .to_expression_at_with_operation(
                tree,
                node,
                &mut |_, value| match value {
                    NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => Ok(Some(
                        network.store.tensors[*index]
                            .expression
                            .as_atom_view()
                            .to_owned(),
                    )),
                    NetworkNode::Leaf(NetworkLeaf::Scalar(index)) => Ok(Some(
                        network
                            .store
                            .get_scalar_ref(*index)
                            .as_atom_view()
                            .to_owned(),
                    )),
                    NetworkNode::Op(_) => Ok(None),
                    _ => unreachable!("symbolic graph admission requires stored leaves"),
                },
                &mut |node, operation, children| {
                    // A power without slot edges raises a closed scalar. Keep
                    // its contraction atomic, as `pow` does: Symbolica would
                    // otherwise distribute the exponent over the factors and
                    // turn `(p·q)^2` into `p^2 q^2`.
                    let NetworkOp::Power(power) = operation else {
                        return Ok(None);
                    };
                    if network
                        .graph
                        .graph
                        .iter_crown(node)
                        .any(|hedge| network.graph.graph[[&hedge]].is_slot())
                    {
                        return Ok(None);
                    }
                    Ok(Some(
                        SymbolicTensor::<PartialStructure>::scalar_power_expression(
                            children[0].clone(),
                            None,
                        )
                        .pow(Atom::num(*power)),
                    ))
                },
            )
            .map_err(|error| TensorInferenceError::invalid(error.to_string()))
    }

    fn emitted_tensor<E: AtomCore>(
        &self,
        network: &PlanningNetwork<E>,
        tree: &CollectionTree,
        node: NodeIndex,
        interface: PartialStructure,
        history: Arc<ExplicitIndexOccurrences>,
    ) -> Result<SymbolicTensor<PartialStructure>> {
        Ok(Self::tensor_for_emission(
            self.emit(network, tree, node)?,
            network,
            tree,
            node,
            interface,
            history,
        ))
    }

    fn tensor_for_emission<E: AtomCore>(
        expression: Atom,
        network: &PlanningNetwork<E>,
        tree: &CollectionTree,
        node: NodeIndex,
        interface: PartialStructure,
        history: Arc<ExplicitIndexOccurrences>,
    ) -> SymbolicTensor<PartialStructure> {
        #[cfg(test)]
        PUBLISHED_TAPE_REGIONS.with(|count| count.set(count.get() + 1));
        let mut tensor = SymbolicTensor::from_validated_parts(expression, interface, history);
        network.inherit_leaf_interfaces(
            &mut tensor.proofs,
            tree.iter_preorder_tree_nodes(&network.graph.graph, node)
                .filter_map(|node| match &network.graph.graph[node] {
                    NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => {
                        Some(&network.store.tensors[*index])
                    }
                    _ => None,
                }),
        );
        tensor
    }

    fn analyze<E: AtomCore>(
        &mut self,
        network: &PlanningNetwork<E>,
        tree: &CollectionTree,
        root: NodeIndex,
    ) -> Result<CollectionAnalysis> {
        let traversal = tree
            .iter_preorder_tree_nodes(&network.graph.graph, root)
            .collect::<Vec<_>>();
        let mut owned = HashMap::new();
        let mut interfaces = HashMap::<_, PartialStructure>::new();
        let mut histories = HashMap::<_, Arc<ExplicitIndexOccurrences>>::new();
        for &node in traversal.iter().rev() {
            #[cfg(test)]
            ANALYZED_NODES.with(|count| count.set(count.get() + 1));
            let children = tree
                .iter_children(node, &network.graph.graph)
                .collect::<Vec<_>>();
            let (selected, interface) = match &network.graph.graph[node] {
                NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => {
                    let tensor = &network.store.tensors[*index];
                    let canonical = tensor.structure.external_structure();
                    let positions = network
                        .graph
                        .logical_port_order(node)
                        .map_err(|error| TensorInferenceError::invalid(error.to_string()))?;
                    let interface = PartialStructure::from_logical_slots(
                        positions.into_iter().map(|position| {
                            let slot = canonical[position];
                            slot.rep().slot(PartialIndex::Explicit(slot.aind()))
                        }),
                    );
                    (
                        self.selected(tensor.expression.as_atom_view()),
                        InterfaceInference::merge_explicit_interface_sequence(&[interface])?,
                    )
                }
                NetworkNode::Leaf(NetworkLeaf::Scalar(index)) => (
                    self.selected(network.store.get_scalar_ref(*index).as_atom_view()),
                    PartialStructure::from_logical_slots([]),
                ),
                NetworkNode::Op(operation) => {
                    let selected = children.iter().any(|child| owned[child]);
                    let child_interfaces = children
                        .iter()
                        .map(|child| interfaces[child].clone())
                        .collect::<Vec<_>>();
                    let interface = match operation {
                        NetworkOp::Product => {
                            InterfaceInference::merge_explicit_interface_sequence(
                                &child_interfaces,
                            )?
                        }
                        NetworkOp::Power(power) if power % 2 == 0 => {
                            PartialStructure::from_logical_slots([])
                        }
                        _ => child_interfaces
                            .first()
                            .cloned()
                            .unwrap_or_else(|| PartialStructure::from_logical_slots([])),
                    };
                    (selected, interface)
                }
                _ => {
                    return Err(TensorInferenceError::invalid(
                        "selected collection requires stored symbolic leaves",
                    ));
                }
            };
            let mut history = match &network.graph.graph[node] {
                NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => {
                    network.store.tensors[*index].history.as_ref().clone()
                }
                NetworkNode::Leaf(NetworkLeaf::Scalar(_)) => ExplicitIndexOccurrences::default(),
                NetworkNode::Op(NetworkOp::Product) => {
                    let mut history = ExplicitIndexOccurrences::default();
                    for child in &children {
                        // Each graph child owns its internal contractions. In
                        // particular, parser-lowered power copies can spell
                        // the same internal label before scope exposure runs.
                        // The enclosing carrier retains its complete history.
                        history.append(&histories[child].scoped(&interfaces[child]));
                    }
                    history
                }
                NetworkNode::Op(NetworkOp::Sum) => {
                    let mut history = ExplicitIndexOccurrences::default();
                    for child in &children {
                        history.branch_max(&histories[child].scoped(&interfaces[child]));
                    }
                    history
                }
                NetworkNode::Op(NetworkOp::Power(power)) => {
                    let mut history = ExplicitIndexOccurrences::default();
                    if let Some(child) = children.first() {
                        if interfaces[child].canonical().is_scalar() {
                            history = histories[child].scoped(&interfaces[child]);
                        } else {
                            let boundary = histories[child].scoped(&interfaces[child]);
                            for _ in 0..(*power).max(0) {
                                history.append(&boundary);
                            }
                        }
                    }
                    history
                }
                NetworkNode::Op(_) => children
                    .first()
                    .map(|child| histories[child].scoped(&interface))
                    .unwrap_or_default(),
                _ => unreachable!("stored symbolic leaves checked at admission"),
            };
            history.branch_max(&ExplicitIndexOccurrences::from_interface(&interface));
            history.validate()?;
            owned.insert(node, selected);
            interfaces.insert(node, interface);
            histories.insert(node, Arc::new(history));
        }
        Ok((owned, interfaces, histories))
    }

    fn leaf(&mut self, value: SymbolicTensor<PartialStructure>, owned: bool) -> TermLeaf<usize> {
        let key = (
            value.expression.clone(),
            value.structure.logical_slots(),
            owned,
        );
        let position = if let Some(&position) = self.positions.get(&key) {
            position
        } else {
            let position = self.values.len();
            self.positions.insert(key, position);
            self.values.push((value, owned));
            position
        };
        self.tape_leaf(position)
    }

    fn tape_leaf(&self, position: usize) -> TermLeaf<usize> {
        if let Ok(number) = Rational::try_from(self.values[position].0.expression.as_atom_view()) {
            TermLeaf::Number(position, number)
        } else {
            TermLeaf::Value(position)
        }
    }

    fn scope_is_intrinsic<E: AtomCore>(
        &self,
        network: &PlanningNetwork<E>,
        tree: &CollectionTree,
        node: NodeIndex,
    ) -> bool {
        tree.iter_preorder_tree_nodes(&network.graph.graph, node)
            .all(|node| match &network.graph.graph[node] {
                NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => {
                    let value = network.store.tensors[*index].expression.as_atom_view();
                    InterfaceInference::normalization_is_intrinsic(value)
                }
                NetworkNode::Leaf(NetworkLeaf::Scalar(index)) => {
                    let value = network.store.get_scalar_ref(*index).as_atom_view();
                    InterfaceInference::normalization_is_intrinsic(value)
                }
                NetworkNode::Op(NetworkOp::Function(symbol)) => {
                    InterfaceInference::intrinsic_normalization_head(*symbol)
                }
                NetworkNode::Op(_) => true,
                _ => false,
            })
    }

    /// Whether distributing this selected product joins a pair which one of its
    /// selected sum factors closes with another factor using the same name.
    /// Each sum binds its own dummies, so such a name may legitimately recur.
    fn distribution_joins_closed_pair<E: AtomCore>(
        &mut self,
        network: &PlanningNetwork<E>,
        tree: &CollectionTree,
        product: NodeIndex,
        analysis: &CollectionAnalysis,
    ) -> bool {
        let (owned, interfaces, _) = analysis;
        let children = tree
            .iter_children(product, &network.graph.graph)
            .collect::<Vec<_>>();
        let names = children
            .iter()
            .map(|&child| self.distributed_names(network, tree, child, analysis))
            .collect::<Vec<_>>();
        children.iter().enumerate().any(|(position, &child)| {
            let sum = match &network.graph.graph[child] {
                NetworkNode::Op(NetworkOp::Sum) => true,
                NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => matches!(
                    network.store.tensors[*index].expression.as_atom_view(),
                    AtomView::Add(_)
                ),
                _ => false,
            };
            if !owned[&child] || !sum {
                return false;
            }
            let boundary = ExplicitIndexOccurrences::from_interface(&interfaces[&child]);
            names[position].iter().any(|slot| {
                boundary.count(*slot) == 0
                    && names
                        .iter()
                        .enumerate()
                        .any(|(other, written)| other != position && written.contains(slot))
            })
        })
    }

    /// Identities which Monomial distribution of `node` places into its rows,
    /// keyed by base orientation. Only selected sums and products are opened,
    /// as compilation opens them; every other node contributes its established
    /// history, which exports only the boundary ports of a closed sum.
    fn distributed_names<E: AtomCore>(
        &mut self,
        network: &PlanningNetwork<E>,
        tree: &CollectionTree,
        node: NodeIndex,
        analysis: &CollectionAnalysis,
    ) -> HashSet<LibrarySlot<AbstractIndex>> {
        let (owned, _, histories) = analysis;
        let mut names = HashSet::new();
        let mut slots = SlotMatcher::default();
        let mut pending = vec![node];
        while let Some(node) = pending.pop() {
            match &network.graph.graph[node] {
                NetworkNode::Op(NetworkOp::Sum | NetworkOp::Product) if owned[&node] => {
                    pending.extend(tree.iter_children(node, &network.graph.graph));
                }
                NetworkNode::Leaf(NetworkLeaf::LocalTensor(index))
                    if owned[&node]
                        && matches!(
                            network.store.tensors[*index].expression.as_atom_view(),
                            AtomView::Add(_) | AtomView::Mul(_)
                        ) =>
                {
                    self.selected_arithmetic_names(
                        network.store.tensors[*index].expression.as_atom_view(),
                        &mut slots,
                        &mut names,
                    );
                }
                _ => names.extend(
                    histories[&node]
                        .entries()
                        .map(|(slot, _)| slot.rep().base().slot(slot.aind())),
                ),
            }
        }
        names
    }

    /// The same frontier inside a selected composite leaf, which compilation
    /// opens: a selected sum or product is distributed into its parts, and any
    /// other part keeps the occurrences its own scope exports.
    fn selected_arithmetic_names(
        &mut self,
        value: AtomView<'_>,
        slots: &mut SlotMatcher,
        names: &mut HashSet<LibrarySlot<AbstractIndex>>,
    ) {
        match value {
            AtomView::Add(sum) if self.selected(value) => {
                for term in sum.iter() {
                    self.selected_arithmetic_names(term, slots, names);
                }
            }
            AtomView::Mul(product) if self.selected(value) => {
                for factor in product.iter() {
                    self.selected_arithmetic_names(factor, slots, names);
                }
            }
            _ => names.extend(
                ExplicitIndexOccurrences::from_atom(value, slots)
                    .entries()
                    .map(|(slot, _)| slot.rep().base().slot(slot.aind())),
            ),
        }
    }

    fn expose_scope<E: AtomCore + Clone + From<E::Output>>(
        &mut self,
        network: &PlanningNetwork<E>,
        tree: &CollectionTree,
        node: NodeIndex,
    ) -> Result<Option<PlanningNetwork<Atom>>> {
        // Prove this graph scope intrinsic before emitting a temporary body:
        // function reconstruction can itself invoke a registered normalizer.
        if !self.scope_is_intrinsic(network, tree, node) {
            return Ok(None);
        }
        if let NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) = &network.graph.graph[node] {
            let expression = network.store.tensors[*index].expression.as_atom_view();
            if matches!(expression, AtomView::Add(_) | AtomView::Mul(_)) {
                // Expose the selected arithmetic before assigning this copy's
                // internal labels; an opaque leaf has only its public ports.
                let opened = self.open(
                    expression,
                    network.store.tensors[*index].source(expression, network.admission.as_ref()),
                )?;
                let opened_tree = opened.graph.expr_tree().cast::<ChildVecStore<()>>();
                let opened_root = opened.graph.graph.node_id(opened.graph.head());
                if !matches!(opened.graph.graph[opened_root], NetworkNode::Leaf(_)) {
                    return self.expose_scope(&opened, &opened_tree, opened_root);
                }
            }
        }
        if network.graph.has_bound_ports() {
            return Ok(None);
        }
        let mut subset: SuBitGraph = network.graph.graph.empty_subgraph();
        for child in tree.iter_preorder_tree_nodes(&network.graph.graph, node) {
            if let NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) = &network.graph.graph[child]
            {
                let tensor = &network.store.tensors[*index];
                for (slot, _) in tensor.history.entries() {
                    self.state.reserve_index(slot.aind());
                }
                for index in tensor.structure.external_indices_iter() {
                    self.state.reserve_index(index);
                }
            }
            for hedge in network.graph.graph.iter_crown(child) {
                subset.add(hedge);
            }
        }
        let mut graph = network.map_ref(|scalar| scalar, |tensor| tensor);
        graph.graph = graph.graph.extract(&subset);
        let mut boundary = Vec::new();
        for index in 0..graph.graph.graph.n_hedges() {
            let hedge = Hedge(index);
            if let NetworkEdge::Slot(slot) = graph.graph.graph[[&hedge]] {
                self.state.reserve_index(slot.aind());
                if graph.graph.graph.inv(hedge) == hedge {
                    boundary.push((hedge, slot));
                }
            }
        }
        let mut visited = graph
            .graph
            .relabel_slot_components(&boundary)
            .map_err(|error| TensorInferenceError::invalid(error.to_string()))?;
        for index in 0..graph.graph.graph.n_hedges() {
            let hedge = Hedge(index);
            if !visited.contains(&hedge)
                && let NetworkEdge::Slot(slot) = graph.graph.graph[[&hedge]]
            {
                visited.extend(
                    graph
                        .graph
                        .relabel_slot_components(&[(
                            hedge,
                            slot.rep().slot(self.state.fresh_index()),
                        )])
                        .map_err(|error| TensorInferenceError::invalid(error.to_string()))?,
                );
            }
        }
        let mut unavailable = false;
        let mapped = graph.map_occurrences_with_structure(
            |tensor| &tensor.structure,
            |scalar| Ok(scalar.as_atom_view().to_owned()),
            |tensor, ports, logical_order| {
                let fail = |message: &str| TensorNetworkError::Other(eyre::eyre!("{message}"));
                let Some(logical_order) = logical_order else {
                    return Err(fail("parsed occurrence has no logical layout"));
                };
                if ports.iter().any(|(_, _, bound)| bound.is_some()) {
                    unavailable = true;
                    return Err(fail("scope exposure requires an unbound occurrence"));
                }
                let source = tensor.structure.external_structure();
                let source_interface =
                    PartialStructure::from_logical_slots(logical_order.iter().map(|&position| {
                        let slot = source[position];
                        let index = match slot.aind() {
                            AbstractIndex::Open { axis, .. } => PartialIndex::open(axis),
                            index => PartialIndex::Explicit(index),
                        };
                        slot.rep().slot(index)
                    }));
                let targets = ports.iter().map(|(_, slot, _)| *slot).collect::<Vec<_>>();
                let replacements = logical_order
                    .iter()
                    .enumerate()
                    .map(|(logical, &storage)| (logical, targets[storage].to_atom()))
                    .collect::<HashMap<_, _>>();
                let history = &tensor.history;
                let stored_keys = source
                    .iter()
                    .map(|slot| slot.rep().base().slot::<AbstractIndex, _>(slot.aind()))
                    .collect::<HashSet<_>>();
                let mut internal = HashMap::new();
                for (slot, _) in history.entries() {
                    let key = slot.rep().base().slot::<AbstractIndex, _>(slot.aind());
                    if !stored_keys.contains(&key) {
                        internal.insert(key, self.state.fresh_index());
                    }
                }
                let mut source_tensor = SymbolicTensor::from_validated_parts(
                    tensor.expression.as_atom_view().to_owned(),
                    source_interface,
                    Arc::clone(history),
                );
                network.inherit_leaf_interfaces(&mut source_tensor.proofs, [*tensor]);
                if !internal.is_empty() {
                    // Opaque leaves may retain closed pairs absent from storage
                    // ports. Map those identities and their regional witnesses
                    // together, including uses erased by normalization.
                    let map_slot = |slot: LibrarySlot<AbstractIndex>| {
                        let key = slot.rep().base().slot::<AbstractIndex, _>(slot.aind());
                        slot.rep()
                            .slot(internal.get(&key).copied().unwrap_or(slot.aind()))
                    };
                    let (expression, mapped) = source_tensor
                        .map_admitted_indices(
                            |value| {
                                if let Ok(slot) = LibrarySlot::<AbstractIndex>::try_from(value) {
                                    return Ok(Some(map_slot(slot).to_atom()));
                                }
                                Ok(Representation::<LibraryRep>::try_from(value)
                                    .is_ok()
                                    .then(|| value.to_owned()))
                            },
                            |interface| {
                                Ok(PartialStructure::from_logical_slots(
                                    interface.logical_slots().into_iter().map(|slot| {
                                        match slot.aind {
                                            PartialIndex::Explicit(index) => {
                                                let mapped = map_slot(slot.rep().slot(index));
                                                mapped
                                                    .rep()
                                                    .slot(PartialIndex::Explicit(mapped.aind()))
                                            }
                                            PartialIndex::Open(_) => slot,
                                        }
                                    }),
                                ))
                            },
                        )
                        .map_err(|error| fail(&error.to_string()))?;
                    let history = history
                        .mapped(|slot| Ok::<_, TensorInferenceError>(map_slot(slot)))
                        .map_err(|error| fail(&error.to_string()))?;
                    let mut renamed = SymbolicTensor::from_validated_parts(
                        expression,
                        source_tensor.structure,
                        Arc::new(history),
                    );
                    renamed.proofs.inherit_interfaces([&mapped]);
                    source_tensor = renamed;
                }
                let (expression, rewritten) = composition::rewrite_interface_parts(
                    &source_tensor,
                    &replacements,
                    Some(&self.state),
                )
                .map_err(|error| fail(&error.to_string()))?;
                let expected =
                    PartialStructure::from_logical_slots(logical_order.iter().map(|&position| {
                        let slot = targets[position];
                        slot.rep().slot(PartialIndex::Explicit(slot.aind()))
                    }));
                #[cfg(debug_assertions)]
                TensorProofs::assert_normalization_preserves_interface(
                    "collection scope exposure",
                    &[&source_tensor.expression],
                    &expression,
                    Ok(&expected),
                    None,
                );
                let history = source_tensor
                    .proofs
                    .occurrences
                    .get()
                    .expect("scope source retains mapped construction history")
                    .mapped(|slot| {
                        Ok::<_, TensorInferenceError>(
                            source
                                .iter()
                                .position(|source| *source == slot)
                                .map_or(slot, |position| targets[position]),
                        )
                    })
                    .map_err(|error| fail(&error.to_string()))?;
                let mut proofs = source_tensor.proofs.clone();
                proofs.occurrences = OnceLock::from(Arc::new(history));
                proofs.publish_interfaces(expression.as_view(), expected, [&rewritten]);
                Ok(SymbolicTensor {
                    expression,
                    structure: OrderedStructure::new(targets).into_canonical(),
                    is_metric: tensor.is_metric,
                    is_composite: tensor.is_composite,
                    work: Default::default(),
                    proofs,
                })
            },
        );
        match mapped {
            Ok(mapped) => Ok(Some(PlanningNetwork::from_materialized(mapped))),
            Err(_) if unavailable => Ok(None),
            Err(error) => Err(TensorInferenceError::invalid(error.to_string())),
        }
    }

    /// Compile parser-lowered power arms after each copy owns its internal labels.
    fn compile_lowered_power<E: AtomCore + Clone + From<E::Output>>(
        &mut self,
        tape: &mut TermTape<usize>,
        network: &PlanningNetwork<E>,
        tree: &CollectionTree,
        root: NodeIndex,
    ) -> Result<Option<(usize, (usize, usize))>> {
        let mut copies = Vec::new();
        for child in tree.iter_children(root, &network.graph.graph) {
            let Some(graph) = self.expose_scope(network, tree, child)? else {
                return Ok(None);
            };
            let tree = graph.graph.expr_tree().cast::<ChildVecStore<()>>();
            let root = graph.graph.graph.node_id(graph.graph.head());
            let analysis = self.analyze(&graph, &tree, root)?;
            let Some(copy) = self.compile(tape, &graph, &tree, root, false, &analysis)? else {
                return Ok(None);
            };
            copies.push(copy);
        }
        Ok(tape.group(copies.into_iter(), false))
    }

    /// A sum or product node, or a leaf whose payload is one: what Monomial
    /// distribution opens when it is selected.
    fn is_composite<E: AtomCore>(network: &PlanningNetwork<E>, node: NodeIndex) -> bool {
        match &network.graph.graph[node] {
            NetworkNode::Op(NetworkOp::Sum | NetworkOp::Product) => true,
            NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => matches!(
                network.store.tensors[*index].expression.as_atom_view(),
                AtomView::Add(_) | AtomView::Mul(_)
            ),
            _ => false,
        }
    }

    /// Whether `value` writes an explicit index outside scalar metadata,
    /// including indices of pairs it closes internally.
    fn writes_explicit_index(value: AtomView<'_>) -> bool {
        let mut slots = SlotMatcher::default();
        let mut written = false;
        value.visitor(&mut |node| {
            written |= slots.parse::<LibraryRep, AbstractIndex>(node).is_ok();
            !written && !node.get_symbol().is_some_and(|symbol| symbol.is_scalar())
        });
        written
    }

    /// Whether the selector picks a variable or a function head inside
    /// `value`; scalar variables are selectable leaves too.
    fn selects_within(&mut self, value: AtomView<'_>) -> bool {
        if self.selected(value) {
            return true;
        }
        let mut selected = false;
        value.visitor(&mut |node| {
            if selected {
                return false;
            }
            if matches!(node, AtomView::Var(_) | AtomView::Fun(_)) {
                selected = (self.select)(node);
                return false;
            }
            true
        });
        selected
    }

    /// Names which `value` places into the product of a Monomial row. Selected
    /// sums and products are distributed into their parts; an unselected sum
    /// is a closed scope that exports only its free names; a closed power
    /// exports nothing, while an open square keeps its name in the row.
    fn row_names(
        &mut self,
        value: AtomView<'_>,
        slots: &mut SlotMatcher,
        names: &mut HashSet<LibrarySlot<AbstractIndex>>,
    ) {
        let key = |slot: LibrarySlot<AbstractIndex>| slot.rep().base().slot(slot.aind());
        match value {
            AtomView::Add(sum) if self.selects_within(value) => {
                for term in sum.iter() {
                    self.row_names(term, slots, names);
                }
            }
            AtomView::Mul(product) => {
                for factor in product.iter() {
                    self.row_names(factor, slots, names);
                }
            }
            AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.bracket => {
                for argument in function.iter() {
                    self.row_names(argument, slots, names);
                }
            }
            AtomView::Fun(function) if function.get_symbol().is_scalar() => {}
            AtomView::Fun(function) => {
                for argument in function.iter() {
                    if let Ok(slot) = slots.parse::<LibraryRep, AbstractIndex>(argument) {
                        names.insert(key(slot));
                    }
                }
            }
            // An open base keeps its free names in the row (an open square
            // uses them twice); a closed scalar power binds its own pairs.
            AtomView::Pow(power) => names.extend(
                ExplicitIndexOccurrences::from_atom(power.get_base_exp().0, slots)
                    .entries()
                    .filter(|&(_, count)| count == 1)
                    .map(|(slot, _)| key(slot)),
            ),
            AtomView::Num(_) | AtomView::Var(_) => {}
            _ => names.extend(
                ExplicitIndexOccurrences::from_atom(value, slots)
                    .entries()
                    .map(|(slot, _)| key(slot)),
            ),
        }
    }

    /// The number of selected factors of this selected product, and whether
    /// distributing it would join names: a selected sum or composite leaf
    /// closes a name which another factor also writes. Parser-lowered square
    /// arms with internal pairs always do. Names are read from the emitted
    /// factors, so pairs closed by an open square inside a term count too.
    fn selected_factors_share_closed_names<E: AtomCore + Clone + From<E::Output>>(
        &mut self,
        network: &PlanningNetwork<E>,
        tree: &CollectionTree,
        product: NodeIndex,
        analysis: &CollectionAnalysis,
    ) -> Result<(usize, bool)> {
        let (owned, interfaces, _) = analysis;
        let children = tree
            .iter_children(product, &network.graph.graph)
            .collect::<Vec<_>>();
        let selected = children.iter().filter(|child| owned[child]).count();
        let distributed = |child: NodeIndex| owned[&child] && Self::is_composite(network, child);
        if !children.iter().any(|&child| distributed(child)) {
            return Ok((selected, false));
        }
        let mut slots = SlotMatcher::default();
        let mut written = Vec::new();
        for &child in &children {
            let expression = self.emit(network, tree, child)?;
            let mut names = HashSet::new();
            self.row_names(expression.as_view(), &mut slots, &mut names);
            written.push(names);
        }
        let joins = children.iter().enumerate().any(|(position, &child)| {
            if !distributed(child) {
                return false;
            }
            let boundary = ExplicitIndexOccurrences::from_interface(&interfaces[&child]);
            written[position].iter().any(|name| {
                boundary.count(*name) == 0
                    && written
                        .iter()
                        .enumerate()
                        .any(|(other, names)| other != position && names.contains(name))
            })
        });
        Ok((selected, joins))
    }

    /// Compile a selected product whose selected factors are distributed
    /// against each other. Each selected composite factor with explicit
    /// indices receives independent internal labels first; free boundary
    /// names still join across factors.
    fn compile_independent_factors<E: AtomCore + Clone + From<E::Output>>(
        &mut self,
        tape: &mut TermTape<usize>,
        network: &PlanningNetwork<E>,
        tree: &CollectionTree,
        product: NodeIndex,
        analysis: &CollectionAnalysis,
    ) -> Result<Option<(usize, (usize, usize))>> {
        let (owned, _, _) = analysis;
        let mut factors = Vec::new();
        for child in tree.iter_children(product, &network.graph.graph) {
            let composite = Self::is_composite(network, child);
            let indexed = tree
                .iter_preorder_tree_nodes(&network.graph.graph, child)
                .any(|member| {
                    network
                        .graph
                        .graph
                        .iter_crown(member)
                        .any(|hedge| network.graph.graph[[&hedge]].is_slot())
                        || matches!(&network.graph.graph[member],
                        NetworkNode::Leaf(NetworkLeaf::LocalTensor(index))
                            if Self::writes_explicit_index(
                                network.store.tensors[*index].expression.as_atom_view(),
                            ))
                });
            let compiled = if owned[&child] && composite && indexed {
                let Some(graph) = self.expose_scope(network, tree, child)? else {
                    return Ok(None);
                };
                let scoped_tree = graph.graph.expr_tree().cast::<ChildVecStore<()>>();
                let scoped_root = graph.graph.graph.node_id(graph.graph.head());
                let scoped_analysis = self.analyze(&graph, &scoped_tree, scoped_root)?;
                self.compile(
                    tape,
                    &graph,
                    &scoped_tree,
                    scoped_root,
                    false,
                    &scoped_analysis,
                )?
            } else {
                self.compile(tape, network, tree, child, false, analysis)?
            };
            let Some(compiled) = compiled else {
                return Ok(None);
            };
            factors.push(compiled);
        }
        Ok(tape.group(factors.into_iter(), false))
    }

    /// Only an opened sum with a selected descendant is distributed. Until
    /// the first such summand, keep opaque operands off the tape. A predicate
    /// may have selected just the original sum, in which case the whole
    /// opened sum must be emitted as one unselected leaf.
    fn compile_sum(
        &mut self,
        tape: &mut TermTape<usize>,
        sum: AddView<'_>,
        proof: &TensorProofs,
    ) -> Result<Option<TermLeaf<usize>>> {
        let mut children = Vec::with_capacity(sum.get_nargs());
        let mut pending = Vec::new();
        let mut refused = false;
        for network in planning::parse_summands(sum, PlanningSource::root(proof)) {
            let network = network?;
            #[cfg(test)]
            PEAK_SUMMAND_NODES.set(PEAK_SUMMAND_NODES.get().max(network.graph.graph.n_nodes()));
            let tree = network.graph.expr_tree().cast::<ChildVecStore<()>>();
            let root = network.graph.graph.node_id(network.graph.head());
            let analysis = self.analyze(&network, &tree, root)?;
            if refused {
                // The whole graph validates every branch before compiling
                // the tape; a late invalid branch must still report its error.
                continue;
            }
            if children.is_empty() && !analysis.0[&root] {
                let mut sources = PlanningLeafSources::default();
                sources.extend(
                    &network,
                    tree.iter_preorder_tree_nodes(&network.graph.graph, root)
                        .filter_map(|node| match network.graph.graph[node] {
                            NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => {
                                Some(&network.store.tensors[index])
                            }
                            _ => None,
                        }),
                );
                pending.push((
                    self.emit(&network, &tree, root)?,
                    analysis.1[&root].clone(),
                    Arc::clone(&analysis.2[&root]),
                    sources,
                ));
                continue;
            }
            for (expression, interface, history, sources) in pending.drain(..) {
                let mut value =
                    SymbolicTensor::from_validated_parts(expression, interface, history);
                sources.publish(&mut value.proofs);
                let Some(child) = tape.push_leaf(self.leaf(value, false)) else {
                    refused = true;
                    break;
                };
                children.push(child);
            }
            if refused {
                continue;
            }
            let Some(child) = self.compile(tape, &network, &tree, root, false, &analysis)? else {
                refused = true;
                continue;
            };
            children.push(child);
        }
        if refused {
            return Ok(None);
        }
        if let Some((_, interface, _, _)) = pending.first() {
            let interface = interface.clone();
            let mut history = ExplicitIndexOccurrences::from_interface(&interface);
            let mut sources = PlanningLeafSources::default();
            let mut operands = Vec::with_capacity(pending.len());
            for (expression, interface, term_history, term_sources) in pending {
                operands.push(expression);
                history.branch_max(&term_history.scoped(&interface));
                sources.append(term_sources);
            }
            history.validate()?;
            let mut value = SymbolicTensor::from_validated_parts(
                Atom::add_many(&operands),
                interface,
                Arc::new(history),
            );
            sources.publish(&mut value.proofs);
            return Ok(Some(self.leaf(value, false)));
        }
        Ok(tape
            .group(children.into_iter(), true)
            .map(|(node, size)| TermLeaf::Reference(node, size)))
    }

    fn compile<E: AtomCore + Clone + From<E::Output>>(
        &mut self,
        tape: &mut TermTape<usize>,
        network: &PlanningNetwork<E>,
        tree: &CollectionTree,
        root: NodeIndex,
        split_root_sum: bool,
        analysis: &CollectionAnalysis,
    ) -> Result<Option<(usize, (usize, usize))>> {
        // These facts belong to this immutable graph. Selecting a smaller root
        // does not change its leaf interfaces or require another selection walk.
        let (owned, interfaces, histories) = analysis;
        let mut failure = None;
        let mut refused = false;
        let compiled = tape.compile_graph(&network.graph, tree, root, &mut |tape, node, value| {
            if failure.is_some() || refused {
                return None;
            }
            let scoped_power = match value {
                NetworkNode::Op(NetworkOp::Power(power)) if *power > 1 && owned[&node] => {
                    Some(*power)
                }
                _ => None,
            };
            // Existing outer sum terms are independent kernel regions. Sums
            // inside a product or a powered copy stay opaque: enumerating those
            // rows would distribute arithmetic that was still factored.
            let factored = self.mode == CollectionMode::Factored
                && (matches!(value, NetworkNode::Op(NetworkOp::Product))
                    || (matches!(value, NetworkNode::Op(NetworkOp::Sum))
                        && !(split_root_sum && node == root)));
            let function = matches!(value, NetworkNode::Op(NetworkOp::Function(_)));
            // Closed inverse powers belong to a scalar scope, not to the
            // positive-power copy frontier. Inspect the base interface: an
            // even power alone does not prove that its base is closed.
            let scalar_power = matches!(value, NetworkNode::Op(NetworkOp::Power(power)) if *power <= 0)
                && tree
                    .iter_children(node, &network.graph.graph)
                    .next()
                    .is_some_and(|base| interfaces[&base].canonical().is_scalar());
            // Selected factors of a nested product are distributed against
            // each other, as at the collection root, so each composite one
            // needs its own internal labels when a closed pair would meet a
            // name of another factor. The parser lowers an open-base square
            // into two identical product arms: without this, their closed
            // pairs share names and coalesce on the tape.
            if self.mode == CollectionMode::Monomials
                && owned[&node]
                && matches!(value, NetworkNode::Op(NetworkOp::Product))
            {
                // With one selected factor, the nested sum's closed pair would
                // meet an unselected neighbour: refuse, as documented for
                // `collect`, before any row or domain callback exists.
                let compiled =
                    match self.selected_factors_share_closed_names(network, tree, node, analysis) {
                        Ok((selected, true)) if selected > 1 => self
                            .compile_independent_factors(tape, network, tree, node, analysis)
                            .map(Some),
                        Ok((_, true)) => Ok(Some(None)),
                        Ok((_, false)) => Ok(None),
                        Err(error) => Err(error),
                    };
                match compiled {
                    Ok(None) => {}
                    Ok(Some(Some((compiled, size)))) => {
                        return Some(TermLeaf::Reference(compiled, size));
                    }
                    Ok(Some(None)) => {
                        refused = true;
                        return None;
                    }
                    Err(error) => {
                        failure = Some(error);
                        return None;
                    }
                }
            }
            if owned[&node]
                && !matches!(value, NetworkNode::Leaf(_))
                && scoped_power.is_none()
                && !function
                && !factored
                && !scalar_power
            {
                // Distribution would join a pair which a nested sum closes
                // with a neighbouring factor of the same name. Refuse before
                // any row or domain callback exists.
                if matches!(value, NetworkNode::Op(NetworkOp::Product))
                    && self.distribution_joins_closed_pair(network, tree, node, analysis)
                {
                    refused = true;
                }
                return None;
            }
            // A function or closed inverse power owns its complete scope.
            // Emit one opaque leaf without distributing the payload or
            // replaying a user callback during reconstruction.
            if owned[&node]
                && (function || scalar_power)
                && !self.scope_is_intrinsic(network, tree, node)
            {
                refused = true;
                return None;
            }
            let result = (|| {
                if let Some(power) = scoped_power {
                    let mut children = tree.iter_children(node, &network.graph.graph);
                    let base = children.next().ok_or_else(|| {
                        TensorInferenceError::invalid("selected power has no base")
                    })?;
                    if children.next().is_some() {
                        return Err(TensorInferenceError::invalid(
                            "selected power has multiple bases",
                        ));
                    }
                    if matches!(network.graph.graph[base], NetworkNode::Leaf(_)) {
                        let base_value = self.emit(network, tree, base)?;
                        if matches!(base_value.as_view(), AtomView::Add(_) | AtomView::Mul(_)) {
                            // Open the selected power before copying: internal
                            // dummy pairs belong to each copy, not its opaque leaf.
                            let value = self.emitted_tensor(
                                network, tree, node, interfaces[&node].clone(), Arc::clone(&histories[&node]),
                            )?;
                            let opened = self.open(value.expression.as_view(), PlanningSource::root(&value.proofs))?;
                            let tree = opened.graph.expr_tree().cast::<ChildVecStore<()>>();
                            let root = opened.graph.graph.node_id(opened.graph.head());
                            if matches!(value.expression.as_view(), AtomView::Pow(_))
                                && matches!(opened.graph.graph[root], NetworkNode::Op(NetworkOp::Product))
                            {
                                let Some((compiled, size)) = self.compile_lowered_power(tape, &opened, &tree, root)? else {
                                    refused = true;
                                    return Ok(None);
                                };
                                return Ok(Some(TermLeaf::Reference(compiled, size)));
                            }
                            let analysis = self.analyze(&opened, &tree, root)?;
                            return self.compile(tape, &opened, &tree, root, false, &analysis)
                                .map(|compiled| compiled.map(|(node, size)| TermLeaf::Reference(node, size)));
                        }
                    }
                    let mut copies = Vec::new();
                    for _ in 0..power {
                        let Some(graph) =
                            self.expose_scope(network, tree, base)?
                        else {
                            refused = true;
                            return Ok(None);
                        };
                        let tree = graph.graph.expr_tree().cast::<ChildVecStore<()>>();
                        let root = graph.graph.graph.node_id(graph.graph.head());
                        let analysis = self.analyze(&graph, &tree, root)?;
                        let Some(copy) =
                            self.compile(tape, &graph, &tree, root, false, &analysis)?
                        else {
                            refused = true;
                            return Ok(None);
                        };
                        copies.push(copy);
                    }
                    if let Some((node, size)) = tape.group(copies.into_iter(), false) {
                        return Ok(Some(TermLeaf::Reference(node, size)));
                    }
                    refused = true;
                    return Ok(None);
                }
                let expression = self.emit(network, tree, node)?;
                let reopen = owned[&node] && matches!(value, NetworkNode::Leaf(_))
                    && (self.mode == CollectionMode::Monomials
                        || matches!(expression.as_view(), AtomView::Pow(_))
                        || (split_root_sum && node == root && matches!(expression.as_view(), AtomView::Add(_))))
                    && matches!(expression.as_view(), AtomView::Add(_) | AtomView::Mul(_) | AtomView::Pow(_));
                // Emission may normalize a graph operation. Once its exact output
                // is known, repeated opaque occurrences already refer to the same
                // tape slot; only the first one supplies that slot's history/facts.
                if !reopen {
                    let interface = SymbolicTensor::canonical_structure(interfaces[&node].clone());
                    let key = (expression, interface.logical_slots(), owned[&node]);
                    let position = match self.positions.entry(key) {
                        Entry::Occupied(entry) => *entry.get(),
                        Entry::Vacant(entry) => {
                            let value = Self::tensor_for_emission(
                                entry.key().0.clone(), network, tree, node, interface, Arc::clone(&histories[&node]),
                            );
                            let position = self.values.len();
                            entry.insert(position);
                            self.values.push((value, owned[&node]));
                            position
                        }
                    };
                    return Ok(Some(self.tape_leaf(position)));
                }
                let emitted = Self::tensor_for_emission(
                    expression, network, tree, node, interfaces[&node].clone(), Arc::clone(&histories[&node]),
                );
                let expression = emitted.expression.as_view();
                {
                    // This leaf was selected from its retained head inventory.
                    // Preserve the established copy-local dummy machinery when
                    // opening its arithmetic; unrelated leaves remain opaque.
                    if let AtomView::Add(sum) = expression
                        && expression.get_byte_size() >= stream_sum_bytes()
                        && emitted.normalization_is_intrinsic()
                    {
                        let compiled = self.compile_sum(tape, sum, &emitted.proofs)?;
                        refused = compiled.is_none();
                        return Ok(compiled);
                    }
                    let opened = self.open(expression, PlanningSource::root(&emitted.proofs))?;
                    let opened_tree = opened.graph.expr_tree().cast::<ChildVecStore<()>>();
                    let opened_root = opened.graph.graph.node_id(opened.graph.head());
                    if !matches!(opened.graph.graph[opened_root], NetworkNode::Leaf(_)) {
                        if matches!(expression, AtomView::Pow(_))
                            && matches!(opened.graph.graph[opened_root], NetworkNode::Op(NetworkOp::Product))
                        {
                            let Some((compiled, size)) = self.compile_lowered_power(tape, &opened, &opened_tree, opened_root)? else {
                                refused = true;
                                return Ok(None);
                            };
                            return Ok(Some(TermLeaf::Reference(compiled, size)));
                        }
                        let split_root_sum = split_root_sum && node == root
                            && matches!(expression, AtomView::Add(_));
                        let analysis = self.analyze(&opened, &opened_tree, opened_root)?;
                        return self.compile(tape, &opened, &opened_tree, opened_root, split_root_sum, &analysis)
                            .map(|compiled| compiled.map(|(node, size)| TermLeaf::Reference(node, size)));
                    }
                }
                Ok(Some(self.leaf(emitted, owned[&node])))
            })();
            match result {
                Ok(leaf) => leaf,
                Err(error) => {
                    failure = Some(error);
                    None
                }
            }
        });
        if let Some(error) = failure {
            return Err(error);
        }
        Ok(if refused { None } else { compiled })
    }
}

impl SymbolicTensor<PartialStructure> {
    /// Check whether the factor's normalized expression and logical boundary
    /// survive singleton multiplication unchanged. Open ports and arbitrary
    /// normalization use the ordinary checked product path.
    fn private_collection_factor<'plan, E: AtomCore>(
        &self,
        network: &'plan PlanningNetwork<E>,
        node: NodeIndex,
        interface: &PartialStructure,
    ) -> Option<EnclosingLeaf<'plan>> {
        let leaf = network.enclosing_leaf(node, interface, &self.proofs)?;
        if leaf.expression.needs_normalization()
            || !interface.canonical().external_structure_iter().all(|slot| matches!(
                slot.aind, PartialIndex::Explicit(index)
                    if InterfaceInference::partial_index(index, LeafInference::Observe) == slot.aind
            ))
            || !self.normalization_is_intrinsic()
        {
            return None;
        }
        Some(leaf)
    }

    /// Retain the existing shallow Network and its discovered leaf interfaces
    /// with this immutable domain. A changed payload rebuilds topology while
    /// intrinsic replacements reuse the exact unchanged boundary interfaces.
    pub(crate) fn shallow_graph(&self) -> Result<Arc<PlanningNetwork<Atom>>> {
        self.work
            .get_or_init()
            .shallow_graph
            .get_or_init(|| {
                #[cfg(test)]
                SHALLOW_BUILDS.with(|count| count.set(count.get() + 1));
                let network = crate::shorthands::schoonschip::SlotContraction::shallow_network(
                    self.expression.as_view(),
                    self.proofs
                        .is_admitted()
                        .then(|| self.proofs.leaf_interfaces.read_cache().unwrap_or_default())
                        .as_ref(),
                    Some(PlanningSource::root(&self.proofs)),
                )
                .map_err(|error| error.to_string())?;
                Ok(Arc::new(network.into_owned()))
            })
            .as_ref()
            .map(Arc::clone)
            .map_err(|error| TensorInferenceError::invalid(error.clone()))
    }

    /// Compare selected coefficients without distributing either tensor payload.
    ///
    /// Exact matching operands return after checking their interfaces. Remaining
    /// identities use the same complete coefficient visitor and conservative
    /// zero proof; `Inconclusive` never certifies equality.
    pub fn coefficients_equal<const N: usize>(
        &self,
        right: &Self,
        filter: TensorCollectFilter<N>,
    ) -> Result<ConditionResult> {
        let difference = self.try_sub(right)?;
        if self.expression == right.expression {
            return Ok(ConditionResult::True);
        }
        // Select tensor scopes before factoring their coefficients. Factoring
        // inverse tensor bodies here could move copy-local dummy indices into
        // the surrounding numerator.
        difference.coefficients_are_zero(filter)
    }

    /// Conservatively prove that every coefficient of the selected tensors is zero.
    ///
    /// `True` certifies a zero tensor. `False` means a nonzero coefficient was
    /// found, not that the selected tensor structures form an independent basis.
    /// Unproved scalar identities and incomplete collection return `Inconclusive`.
    /// Coefficients retain their factorization; no numerical samples are taken.
    pub fn coefficients_are_zero<const N: usize>(
        &self,
        filter: TensorCollectFilter<N>,
    ) -> Result<ConditionResult> {
        let mut proof = ConditionResult::Inconclusive;
        self.collect_with_map(
            CollectionMode::Monomials,
            Some(&mut |states, spectators| {
                let spectator = Self::collection_product(spectators)?;
                proof = ConditionResult::True;
                for (selected, coefficient) in states.iter() {
                    if selected.expression.is_zero() {
                        continue;
                    }
                    let coefficient =
                        Self::collection_product(&[coefficient.clone(), spectator.clone()])?;
                    // Symbolica's zero-iteration check keeps unresolved sums
                    // inconclusive. Only exact factor cancellation can prove zero.
                    match coefficient.expression.collect_factors().check_zero_with(
                        &symbolica::atom::ZeroTestSettings::default()
                            .finite_field_samples(0)
                            .numerical_samples(0)
                            .exact_fallback(false),
                    ) {
                        ConditionResult::True => {}
                        ConditionResult::False => {
                            proof = ConditionResult::False;
                            break;
                        }
                        ConditionResult::Inconclusive => proof = ConditionResult::Inconclusive,
                    }
                }
                Ok(())
            }),
            |value| filter.matches(value),
            |value, _, _| Ok(value),
        )?;
        Ok(proof)
    }

    /// Return complete selected sectors and their factored coefficients.
    ///
    /// This is an explicit output of the same graph/tape collection used by
    /// `collect`. Scalar
    /// spectators are retained in the coefficient without polynomial expansion.
    /// A refused or capped traversal returns an error, never partial rows.
    pub fn coefficient_list<const N: usize>(
        &self,
        filter: TensorCollectFilter<N>,
    ) -> Result<Vec<(Self, Self)>> {
        self.coefficient_list_with(|value| filter.matches(value))
    }

    /// Collect leaves selected by a predicate, retaining all unselected tensor
    /// blocks in the coefficients with their logical interfaces intact.
    /// The predicate must return the same decision for an unchanged value;
    /// selection facts can be reused without invoking it again. The order
    /// and number of predicate visits are unspecified.
    pub fn coefficient_list_with(
        &self,
        select: impl FnMut(AtomView<'_>) -> bool,
    ) -> Result<Vec<(Self, Self)>> {
        let mut rows = None;
        self.collect_with_map(
            CollectionMode::Monomials,
            Some(&mut |states, spectators| {
                let spectator = Self::collection_product(spectators)?;
                rows = Some(
                    states
                        .iter()
                        .map(|(selected, coefficient)| {
                            Ok((
                                selected.clone(),
                                Self::collection_product(&[
                                    coefficient.clone(),
                                    spectator.clone(),
                                ])?,
                            ))
                        })
                        .collect::<Result<Vec<_>>>()?,
                );
                Ok(())
            }),
            select,
            |value, _, _| Ok(value),
        )?;
        rows.ok_or_else(|| {
            TensorInferenceError::invalid("selected coefficient collection did not complete")
        })
    }

    /// Collect selected tensors with factored coefficients.
    /// Unselected scalar and tensor subtrees remain factored. A dummy pair
    /// which a sum closes keeps its scope: when a neighbouring factor reuses
    /// its name, the collected payload gives that pair a fresh dummy name.
    /// Monomial distribution which would join a nested sum's closed pair with
    /// a neighbouring factor is refused: the input is returned as deferred.
    pub fn collect<const N: usize>(&self, filter: TensorCollectFilter<N>) -> Result<Self> {
        self.collect_with_map(
            CollectionMode::Monomials,
            None,
            |value| filter.matches(value),
            |value, _, _| Ok(value),
        )
    }

    /// Visit selected arithmetic through the graph's established factor order.
    /// Monomial consumers distribute selected sums; factored consumers retain
    /// them as leaves and decide which local identities require distribution.
    /// The boolean is true only when every selected factor has entered this
    /// frontier. A kernel may consume incident coefficient ports; the combined
    /// logical boundary remains invariant. Refused and capped frontiers remain exact.
    /// The finalizer can combine complete selected states before materialization.
    /// Like the kernel, it must publish trusted typed identity results preserving
    /// the combined ordered boundary. Spectators remain immutable. With no selected
    /// work it receives the original coefficient for inspection only.
    pub(crate) fn collect_with_map(
        &self,
        mode: CollectionMode,
        finalize: Option<&mut CollectionFinalizer<'_>>,
        select: impl FnMut(AtomView<'_>) -> bool,
        step: impl FnMut(Self, bool, &mut Self) -> Result<Self>,
    ) -> Result<Self> {
        self.collect_rows(mode, finalize.map(Finalizer::WithSpectators), select, step)
    }

    /// Row kernels never receive the graph's unchanged spectator factors.
    pub(crate) fn collect_rows_with_map(
        &self,
        mode: CollectionMode,
        finalize: &mut dyn FnMut(&mut CollectionStates) -> Result<()>,
        select: impl FnMut(AtomView<'_>) -> bool,
        step: impl FnMut(Self, bool, &mut Self) -> Result<Self>,
    ) -> Result<Self> {
        self.collect_rows(mode, Some(Finalizer::Rows(finalize)), select, step)
    }

    fn collect_rows(
        &self,
        mode: CollectionMode,
        mut finalize: Option<Finalizer<'_, '_>>,
        mut select: impl FnMut(AtomView<'_>) -> bool,
        mut step: impl FnMut(Self, bool, &mut Self) -> Result<Self>,
    ) -> Result<Self> {
        self.explicit_index_occurrences()?;
        let mut input = CollectionInput::new(self, &mut select, mode);
        let network = input.parse(self)?;
        let tree = network.graph.expr_tree().cast::<ChildVecStore<()>>();
        let root = network.graph.graph.node_id(network.graph.head());
        let analysis = input.analyze(&network, &tree, root)?;
        let (owned, interfaces, histories) = &analysis;
        if !owned[&root] {
            if let Some(finalize) = finalize.as_mut() {
                let one = Self::from_validated_parts(
                    Atom::one(),
                    PartialStructure::from_logical_slots([]),
                    Arc::default(),
                );
                finalize.finish(&mut vec![(one, self.clone())], &[])?;
            }
            return Ok(self.clone());
        }
        let nodes = if matches!(
            network.graph.graph[root],
            NetworkNode::Op(NetworkOp::Product)
        ) {
            let mut nodes = tree
                .iter_children(root, &network.graph.graph)
                .collect::<Vec<_>>();
            nodes.reverse();
            nodes
        } else {
            vec![root]
        };
        let mut tape = TermTape::<usize>::default();
        // A finalizer that receives these tensors can clone and retain them.
        // Only row-only kernels may leave descendant facts on the plan.
        let private_spectators = !matches!(finalize, Some(Finalizer::WithSpectators(_)));
        let mut factors = Vec::new();
        let mut spectators = Vec::new();
        let mut spectator_sources = Vec::new();
        let selected_count = nodes.iter().filter(|node| owned[node]).count();
        let expose = mode == CollectionMode::Monomials && selected_count > 1;
        let lowered_power = matches!(self.expression.as_view(), AtomView::Pow(_))
            && matches!(
                network.graph.graph[root],
                NetworkNode::Op(NetworkOp::Product)
            );
        for node in nodes {
            if !owned[&node] {
                let interface = &interfaces[&node];
                let spectator = if private_spectators
                    && let Some(leaf) = self.private_collection_factor(&network, node, interface)
                {
                    if !leaf.boundary_published {
                        // Publish a new scalar root. The enclosing tensor
                        // already supplies all of this factor's descendants.
                        spectator_sources.push(spectators.len());
                    }
                    let value = Self::from_validated_parts(
                        leaf.expression.to_owned(),
                        interface.clone(),
                        Arc::clone(&histories[&node]),
                    );
                    let _ = value.work.intrinsic.set(true);
                    value.observe_replacement_work(self);
                    value
                } else {
                    let value = Self::collection_product(&[input.emitted_tensor(
                        &network,
                        &tree,
                        node,
                        interface.clone(),
                        Arc::clone(&histories[&node]),
                    )?])?;
                    value.observe_replacement(self);
                    spectator_sources.push(spectators.len());
                    value
                };
                if spectator.is_scalar()
                    && let Some(facts) = input
                        .observations
                        .certify_scalar_region(spectator.expression.as_view())
                {
                    let _ = spectator.work.intrinsic.set(true);
                    let _ = spectator
                        .work
                        .get_or_init()
                        .observations
                        .set(Arc::new(facts));
                }
                spectators.push(spectator);
                continue;
            }
            // A selected composite leaf is opened and distributed like a
            // selected sum or product. Its internal pairs have no graph edge,
            // so its payload decides whether it needs independent labels.
            let composite_leaf = match &network.graph.graph[node] {
                NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => {
                    let expression = network.store.tensors[*index].expression.as_atom_view();
                    let mut slots = SlotMatcher::default();
                    let mut written = false;
                    expression.visitor(&mut |value| {
                        written |= slots.parse::<LibraryRep, AbstractIndex>(value).is_ok();
                        !written && !value.get_symbol().is_some_and(|symbol| symbol.is_scalar())
                    });
                    matches!(expression, AtomView::Add(_) | AtomView::Mul(_)) && written
                }
                _ => false,
            };
            let compiled = if (expose && composite_leaf)
                || ((lowered_power
                    || (expose
                        && matches!(
                            network.graph.graph[node],
                            NetworkNode::Op(NetworkOp::Sum | NetworkOp::Product)
                        )))
                    && tree
                        .iter_preorder_tree_nodes(&network.graph.graph, node)
                        .any(|member| {
                            network
                                .graph
                                .graph
                                .iter_crown(member)
                                .any(|hedge| network.graph.graph[[&hedge]].is_slot())
                        })) {
                // A single selected subtree retains its existing dummy scope.
                // Separately expanded factors and parser-lowered power arms
                // need independent internal labels in either collection mode.
                // Keep its selected leaves available even when their local
                // domain step uses a checked normalization callback.
                // The parser can lower a tensor square into two product arms.
                // Scope each already parsed factor's internal incidences before
                // tape coalescence; free boundary names still join across arms.
                let Some(graph) = input.expose_scope(&network, &tree, node)? else {
                    let mut deferred = self.clone();
                    deferred.work.get_mut().frontier =
                        super::simplification::ReductionStatus::Deferred;
                    return Ok(deferred);
                };
                let scoped_tree = graph.graph.expr_tree().cast::<ChildVecStore<()>>();
                let scoped_root = graph.graph.graph.node_id(graph.graph.head());
                let scoped_analysis = input.analyze(&graph, &scoped_tree, scoped_root)?;
                input.compile(
                    &mut tape,
                    &graph,
                    &scoped_tree,
                    scoped_root,
                    false,
                    &scoped_analysis,
                )?
            } else {
                input.compile(
                    &mut tape,
                    &network,
                    &tree,
                    node,
                    // Spectators are multiplied back after the selected rows
                    // are collected. Splitting a sole selected sum therefore
                    // preserves their factorization as well as the sum's.
                    selected_count == 1,
                    &analysis,
                )?
            };
            let Some((compiled, _)) = compiled else {
                // No domain callback has run; discard all refused tape scratch.
                let mut deferred = self.clone();
                deferred.work.get_mut().frontier = super::simplification::ReductionStatus::Deferred;
                return Ok(deferred);
            };
            factors.push((node, compiled));
        }
        let values = &input.values;
        let scalar = PartialStructure::from_logical_slots([]);
        let one = Self::from_validated_parts(Atom::one(), scalar, Arc::default());
        let mut states = vec![(one.clone(), one)];
        let mut remaining = Vec::new();
        let mut frontier = super::simplification::ReductionStatus::Complete;
        for (level, &(_, factor)) in factors.iter().enumerate() {
            tape.clear_terms();
            if tape
                .distribute(&mut vec![(factor, 1)], &Rational::one())
                .is_none()
            {
                frontier = frontier.max(super::simplification::ReductionStatus::Deferred);
                remaining.extend(factors[level..].iter().map(|(node, _)| *node));
                break;
            }
            let rows = tape.take_terms();
            let mut next = Vec::<(Self, Vec<Self>)>::new();
            let mut keys = HashMap::new();
            for (selected, coefficient) in states {
                for (powers, number) in &rows {
                    let mut selected_factors = vec![selected.clone()];
                    let mut coefficient_factors = vec![coefficient.clone()];
                    for &(position, power) in powers {
                        let (value, owned) = &values[tape.leaves()[position]];
                        if power == 1 {
                            if *owned {
                                selected_factors.push(value.clone());
                            } else {
                                coefficient_factors.push(value.clone());
                            }
                            continue;
                        }
                        let value = value.collection_power(power)?;
                        if *owned {
                            selected_factors.push(value);
                        } else {
                            coefficient_factors.push(value);
                        }
                    }
                    let monomial = Self::collection_product(&selected_factors)?;
                    monomial.observe_replacement(self);
                    let mut coefficient = Self::collection_product(&coefficient_factors)?;
                    let original_coefficient = coefficient.clone();
                    let mapped = if frontier == super::simplification::ReductionStatus::Capped {
                        // Finish retaining the already distributed row frontier,
                        // but do not invoke another kernel after its budget ended.
                        monomial.clone()
                    } else {
                        step(
                            monomial.clone(),
                            level + 1 == factors.len(),
                            &mut coefficient,
                        )?
                    };
                    frontier = frontier.max(mapped.work.frontier());
                    // Internal domain owners publish certified typed results. Only
                    // unchecked callbacks need to have their payload observed again.
                    let expected = if coefficient.structure == original_coefficient.structure {
                        monomial.structure.clone()
                    } else {
                        // A kernel may consume incident coefficient ports while
                        // reducing its selected word. Their combined boundary is
                        // the invariant, including logical port order.
                        let before = InterfaceInference::merge_explicit_interface_sequence(&[
                            monomial.structure.clone(),
                            original_coefficient.structure,
                        ])?;
                        let after = InterfaceInference::merge_explicit_interface_sequence(&[
                            mapped.structure.clone(),
                            coefficient.structure.clone(),
                        ])?;
                        if after != before {
                            return Err(TensorInferenceError::invalid(
                                "selected domain changes its combined coefficient interface",
                            ));
                        }
                        mapped.structure.clone()
                    };
                    if mapped.structure != expected {
                        return Err(TensorInferenceError::invalid(
                            "selected domain changes its established logical interface",
                        ));
                    }
                    let mapped = if mapped.proofs.is_admitted() {
                        mapped
                    } else if monomial.structure.open_positions().is_empty() {
                        monomial.with_checked_expression(mapped.expression)?
                    } else {
                        monomial.with_rewritten_expression(mapped.expression)?
                    };
                    let facts = Self::collection_observations(
                        std::slice::from_ref(&coefficient),
                        &coefficient.structure,
                        true,
                    );
                    let coefficient = coefficient.with_identity_result(
                        Atom::num(number.clone()) * &coefficient.expression,
                        facts,
                    )?;
                    let key = (mapped.expression.clone(), mapped.structure.logical_slots());
                    if let Some(&position) = keys.get(&key) {
                        let (_, coefficients): &mut (Self, Vec<Self>) = &mut next[position];
                        coefficients.push(coefficient);
                    } else {
                        keys.insert(key, next.len());
                        next.push((mapped, vec![coefficient]));
                    }
                }
            }
            states = next
                .into_iter()
                .map(|(selected, coefficients)| {
                    let first = &coefficients[0];
                    if coefficients.iter().any(|value| {
                        !InterfaceInference::additive_interfaces_match(
                            &first.structure,
                            &value.structure,
                        )
                    }) {
                        return Err(TensorInferenceError::invalid(
                            "collected coefficients have different interfaces",
                        ));
                    }
                    let body = first.with_identity_result(
                        Atom::add_many(coefficients.iter().map(|value| &value.expression)),
                        Self::collection_observations(&coefficients, &first.structure, true),
                    )?;
                    Ok((selected, body))
                })
                .collect::<Result<Vec<_>>>()?;
            if frontier == super::simplification::ReductionStatus::Capped {
                remaining.extend(factors[level + 1..].iter().map(|(node, _)| *node));
                break;
            }
        }
        // Only complete coefficient rows can be finalized or certify queries.
        if remaining.is_empty()
            && frontier != super::simplification::ReductionStatus::Capped
            && let Some(finalize) = finalize.as_mut()
        {
            finalize.finish(&mut states, &spectators)?;
        }
        let mut terms = states
            .into_iter()
            .map(|(selected, coefficient)| Self::collection_product(&[selected, coefficient]))
            .collect::<Result<Vec<_>>>()?;
        let single_term = terms.len() == 1;
        if let [term] = terms.as_mut_slice() {
            // The rows of a sum have joined into one term, so the pairs which
            // that sum closed now meet the surrounding factors.
            let mut neighbours = spectators
                .iter()
                .map(Self::explicit_index_occurrences)
                .collect::<Result<Vec<_>>>()?;
            neighbours.extend(remaining.iter().map(|node| histories[node].as_ref()));
            let renaming = term.colliding_closed_pairs(&neighbours, &input.state)?;
            if !renaming.is_empty() {
                if !term.normalization_is_intrinsic() {
                    // Rebuilding the payload could replay a registered normalizer.
                    let mut deferred = self.clone();
                    deferred.work.get_mut().frontier =
                        super::simplification::ReductionStatus::Deferred;
                    return Ok(deferred);
                }
                *term = term.with_fresh_dummies(&renaming)?;
            }
        }
        let selected_interface = terms.first().map(|term| &term.structure);
        let mut facts = selected_interface
            .and_then(|interface| Self::collection_observations(&terms, interface, true));
        let selected_expression = Atom::add_many(terms.iter().map(|value| &value.expression));
        let selected_facts = facts.clone();
        let complete = remaining.is_empty();
        let mut output = vec![selected_expression.clone()];
        let mut remaining_proofs = Vec::with_capacity(remaining.len());
        if !remaining.is_empty() {
            facts = None;
        }
        for node in remaining {
            let interface = &interfaces[&node];
            let body = if let Some(leaf) = self.private_collection_factor(&network, node, interface)
            {
                if leaf.boundary_published {
                    output.push(leaf.expression.to_owned());
                    continue;
                }
                Self::from_validated_parts(
                    leaf.expression.to_owned(),
                    interface.clone(),
                    Arc::clone(&histories[&node]),
                )
            } else {
                Self::collection_product(&[input.emitted_tensor(
                    &network,
                    &tree,
                    node,
                    interface.clone(),
                    Arc::clone(&histories[&node]),
                )?])?
            };
            output.push(body.expression);
            remaining_proofs.push(body.proofs);
        }
        facts = facts.and_then(|selected| {
            let mut needs_dots = selected.needs_dot_normalization();
            let mut ports = selected_interface?.logical_slots().len();
            for spectator in &spectators {
                let observed = spectator.reduction_observations();
                if !(observed.is_terminal_polynomial()
                    || (spectator.is_scalar() && observed.can_certify_scalar()))
                {
                    return None;
                }
                ports += spectator.structure.logical_slots().len();
                needs_dots |= observed.needs_dot_normalization();
            }
            // A product can connect ports that were external in its factors.
            // Retain terminal facts only when the established boundary proves
            // that no such new connection was formed.
            (ports == self.structure.logical_slots().len())
                .then(|| DomainObservations::terminal_polynomial(&self.structure, needs_dots))
        });
        output.extend(spectators.iter().map(|value| value.expression.clone()));
        let expression = Atom::mul_many(output);
        if facts.is_none()
            && complete
            && let Some(selected_facts) = selected_facts
        {
            facts = DomainObservations::compose_terminal_product(
                expression.as_view(),
                &self.structure,
                selected_expression.as_view(),
                &selected_facts,
                &spectators,
            );
        }
        let structure = Self::canonical_structure(self.structure.clone());
        let mut proofs = self.proofs.clone();
        proofs.publish_interfaces(
            expression.as_view(),
            structure.clone(),
            terms
                .iter()
                .chain(spectator_sources.iter().map(|&index| &spectators[index]))
                .map(|value| &value.proofs)
                .chain(&remaining_proofs),
        );
        let mut history = InterfaceInference::produced_occurrences(
            expression.as_view(),
            &mut SlotMatcher::default(),
            &proofs,
        )?;
        if lowered_power {
            // The parser allocated independent dummies inside each copy.
            // Their local histories were available to the selected kernel;
            // publishing the power exposes only its established outer scope.
            history = history.scoped(&self.structure);
        }
        history.branch_max(self.explicit_index_occurrences()?);
        history.branch_max(&ExplicitIndexOccurrences::from_interface(&self.structure));
        if single_term
            && let Some(pairs) = ExplicitIndexOccurrences::root_product(expression.as_view())
        {
            // Record the pairs which the collected term contracts in the
            // outermost product, so later compositions cannot reuse them.
            history.branch_max(&pairs);
        }
        history.validate()?;
        let mut root = Self::from_published_parts(
            expression,
            structure,
            Arc::new(history),
            proofs.admitted_interfaces,
            proofs.leaf_interfaces,
        );
        if let Some(facts) = facts {
            root.work.set_observations(facts);
        }
        root.work.get_mut().frontier = frontier;
        // Collection gives the leaves of closed scopes and colliding pairs fresh
        // dummies; publish their raw ports from the source templates.
        self.retain_renamed_leaf_ports(&root)?;
        Ok(root)
    }

    /// Fresh names for the pairs which this collected term contracts while a
    /// neighbouring factor also uses their names. A sum closed these pairs
    /// before its rows joined into this term. Empty when nothing collides.
    fn colliding_closed_pairs(
        &self,
        neighbours: &[&ExplicitIndexOccurrences],
        state: &ParseState<AbstractIndex>,
    ) -> Result<HashMap<LibrarySlot<AbstractIndex>, AbstractIndex>> {
        let boundary = ExplicitIndexOccurrences::from_interface(&self.structure);
        Ok(self
            .explicit_index_occurrences()?
            .entries()
            .filter(|&(slot, _)| {
                boundary.count(slot) == 0
                    && neighbours
                        .iter()
                        .any(|neighbour| neighbour.count(slot) != 0)
            })
            .map(|(slot, _)| {
                (
                    slot.rep().base().slot::<AbstractIndex, _>(slot.aind()),
                    state.fresh_index(),
                )
            })
            .collect())
    }

    /// Rename closed identities in the payload and its history together.
    /// The renamed pairs are internal, so the logical interface is unchanged.
    fn with_fresh_dummies(
        &self,
        renaming: &HashMap<LibrarySlot<AbstractIndex>, AbstractIndex>,
    ) -> Result<Self> {
        let rename = |slot: &LibrarySlot<AbstractIndex>| {
            renaming
                .get(&slot.rep().base().slot::<AbstractIndex, _>(slot.aind()))
                .copied()
        };
        let expression =
            SlotMatcher::default().map_expression_indices(self.expression.as_view(), rename);
        let history = self.explicit_index_occurrences()?.mapped(|slot| {
            Ok::<_, TensorInferenceError>(
                rename(&slot).map_or(slot, |index| slot.rep().slot(index)),
            )
        })?;
        let value =
            Self::from_validated_parts(expression, self.structure.clone(), Arc::new(history));
        let _ = value.work.intrinsic.set(true);
        Ok(value)
    }

    // Sums preserve a common terminal boundary; products preserve terminal facts
    // only when their established interface excludes new inter-factor connections.
    fn collection_observations(
        values: &[Self],
        interface: &PartialStructure,
        additive: bool,
    ) -> Option<DomainObservations> {
        let mut needs_dots = false;
        let mut ports = 0;
        for value in values {
            if additive && value.structure != *interface {
                return None;
            }
            let observed = value.reduction_observations();
            if !(observed.is_terminal_polynomial()
                || (value.is_scalar() && observed.can_certify_scalar()))
            {
                return None;
            }
            ports += value.structure.logical_slots().len();
            needs_dots |= observed.needs_dot_normalization();
        }
        if !additive && ports != interface.logical_slots().len() {
            return None;
        }
        Some(DomainObservations::terminal_polynomial(
            interface, needs_dots,
        ))
    }

    /// Whether the payload's own product closes a pair: a tensor argument at
    /// its product or power level names an identity that the interface does
    /// not export. `pow` distributes over products and nested powers only;
    /// pairs inside a sum or a function's arguments stay in their scope.
    fn closes_explicit_pairs(&self) -> Result<bool> {
        let boundary = ExplicitIndexOccurrences::from_interface(&self.structure);
        let mut slots = SlotMatcher::default();
        let mut closed = false;
        self.expression.as_view().visitor(&mut |node| {
            if closed || matches!(node, AtomView::Add(_)) {
                return false;
            }
            if let AtomView::Fun(function) = node {
                if !function.get_symbol().is_scalar() {
                    for argument in function.iter() {
                        if let Ok(slot) = slots.parse::<LibraryRep, AbstractIndex>(argument) {
                            closed |= boundary.count(slot.rep().base().slot(slot.aind())) == 0;
                        }
                    }
                }
                return false;
            }
            true
        });
        Ok(closed)
    }

    /// Publish a coalesced tape leaf's power using its established scope.
    fn collection_power(&self, power: u32) -> Result<Self> {
        let interface = if power.is_multiple_of(2) {
            PartialStructure::from_logical_slots([])
        } else {
            self.structure.clone()
        };
        // A closed scalar keeps its contraction atomic, as `pow` does:
        // Symbolica would distribute the exponent over its factors and turn
        // `(p·q)^2` or `(p·p)^2` into `p^2 q^2` or `p^4`. An open base which
        // also closes pairs cannot be raised without renaming each copy.
        let expression = if self.is_scalar() {
            Self::scalar_power_expression(self.expression.clone(), None).pow(power)
        } else {
            if power > 1 && self.closes_explicit_pairs()? {
                return Err(TensorInferenceError::invalid(format!(
                    "cannot raise `{}` to the power {power}: its closed pairs need \
                     independent names in each copy",
                    self.expression
                )));
            }
            self.expression.pow(power)
        };
        let intrinsic = self.normalization_is_intrinsic();
        #[cfg(debug_assertions)]
        TensorProofs::assert_normalization_preserves_interface(
            "collection_power",
            &[&self.expression],
            &expression,
            Ok(&interface),
            None,
        );
        let mut history = ExplicitIndexOccurrences::default();
        // A scalar base owns its internal dummy scope. Tape coalescence must
        // not export repeated copies of that history into the outer product.
        if !self.is_scalar() {
            let boundary = self.explicit_index_occurrences()?.scoped(&self.structure);
            for _ in 0..power {
                history.append(&boundary);
            }
            history.validate()?;
        }
        let mut value = Self::from_validated_parts(expression, interface, Arc::new(history));
        value.proofs.inherit_interfaces([&self.proofs]);
        if intrinsic {
            let _ = value.work.intrinsic.set(true);
        }
        Ok(value)
    }

    pub(crate) fn collection_product(values: &[Self]) -> Result<Self> {
        let interface = InterfaceInference::merge_explicit_interface_sequence(
            &values
                .iter()
                .map(|value| value.structure.clone())
                .collect::<Vec<_>>(),
        )?;
        // Graph scratch retains encoded AUTO identities for incidence and
        // exposure. A domain callback receives the declared logical interface;
        // its checked finisher still compares encoded identities in the payload.
        let interface = PartialStructure::from_logical_slots(
            interface.logical_slots().into_iter().map(|slot| {
                let index = match slot.aind {
                    PartialIndex::Explicit(index) => {
                        InterfaceInference::partial_index(index, LeafInference::Observe)
                    }
                    open => open,
                };
                slot.rep().slot(index)
            }),
        );
        let occurrences = ExplicitIndexOccurrences::product(&values.iter().collect::<Vec<_>>())?;
        let expression = Atom::mul_many(values.iter().map(|value| &value.expression));
        if let Some(value) = values.iter().find(|value| {
            value.proofs.is_admitted()
                && value.expression == expression
                && value.structure == interface
        }) {
            // Multiplication by an identity did not change the exact payload or
            // its ordered interface. Keep its observations and lazy graph; an
            // encoded AUTO boundary must still pass the decoding above first.
            return Ok(value.clone());
        }
        let intrinsic = values.iter().all(Self::normalization_is_intrinsic);
        #[cfg(debug_assertions)]
        TensorProofs::assert_normalization_preserves_interface(
            "collection_product",
            &values
                .iter()
                .map(|value| &value.expression)
                .collect::<Vec<_>>(),
            &expression,
            Ok(&interface),
            None,
        );
        let mut value = Self::from_validated_parts(expression, interface, Arc::new(occurrences));
        value
            .proofs
            .inherit_interfaces(values.iter().map(|source| &source.proofs));
        if value.expression.is_zero() {
            // A vanished selected word can leave ports joined to its coefficient.
            // Those nominal connections determine the typed zero's boundary, but
            // there is no surviving syntax to contract or inspect for identities.
            let _ = value.work.intrinsic.set(true);
            let _ = value.work.get_or_init().observations.set(Arc::new(
                DomainObservations::terminal_polynomial(&value.structure, false),
            ));
        } else if intrinsic {
            let _ = value.work.intrinsic.set(true);
            if let Some(observed) = Self::collection_observations(values, &value.structure, false) {
                let _ = value
                    .work
                    .get_or_init()
                    .observations
                    .set(Arc::new(observed));
            }
        }
        Ok(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dirac::AGS;
    use crate::epsilon::EPSILON_SYMBOL;
    use crate::tensor::composition::TensorCompositionError;
    use crate::tensor::inference::tests::{INFERENCE_CALLS, INFERENCE_VISITS};
    use crate::tensor::simplification::ReductionStatus;
    use crate::tensor::{SymbolicNet, SymbolicNetParse};
    use crate::test_support;
    use spenso::mink;
    use spenso::network::library::DummyLibrary;
    use spenso::network::library::symbolic::ETS;
    use spenso::network::parsing::ShorthandParsing;
    use spenso::network::tags::SPENSO_TAG;
    use spenso::shadowing::SYM;
    use spenso::structure::representation::{Euclidean, Minkowski, RepName};
    use std::slice;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use symbolica::atom::{FunctionBuilder, Symbol};
    use symbolica::coefficient::CoefficientView;
    use symbolica::{parse_lit, symbol};

    fn tensor(head: Symbol, slots: &[Atom]) -> Atom {
        FunctionBuilder::new(head).add_args(slots).finish()
    }

    #[test]
    fn streamed_summands_match_whole_tape_and_admitted_leaves() {
        test_support::test_initialize();
        let left = SPENSO_TAG.tensor_symbol("stream_sum_left");
        let right = SPENSO_TAG.tensor_symbol("stream_sum_right");
        let free = mink!(4, 97132);
        let dummy = mink!(4, 97133);
        let loop_leaf = tensor(left, &[free.clone(), dummy.clone(), dummy]);
        let spectator = tensor(right, &[free]);
        let source = SymbolicTensor::infer(
            &loop_leaf * parse_lit!(stream_a + stream_b) + &spectator * parse_lit!(stream_c),
        )
        .unwrap();
        assert_streamed_sum_matches_whole(&source, &[Some(left), Some(right), None]);

        // AUTO slots are already bound to their admitted owner. Both branches
        // must preserve those identities as well as their written port order.
        let gamma = SymbolicTensor::<PartialStructure>::from_signature(
            &AGS.gamma_strct::<AbstractIndex>(4),
        )
        .unwrap();
        let positional = gamma
            .try_scalar_mul(&parse_lit!(stream_s))
            .unwrap()
            .try_add(&gamma.try_scalar_mul(&parse_lit!(stream_t)).unwrap())
            .unwrap();
        assert!(!positional.structure.open_positions().is_empty());
        assert_streamed_sum_matches_whole(&positional, &[gamma.expression.get_symbol(), None]);
    }

    #[test]
    fn streamed_summands_bound_open_graph_size_during_collection() {
        test_support::test_initialize();
        let left = SPENSO_TAG.tensor_symbol("stream_nested_left");
        let right = SPENSO_TAG.tensor_symbol("stream_nested_right");
        let dummy = mink!(4, 97146);
        let slot = mink!(4, 97147);
        // Every term contracts across a product using the same written dummy.
        // The free slot survives into each collected row.
        let tail = tensor(right, &[dummy.clone(), slot]);
        let sum = Atom::add_many(
            (0..96)
                .map(|i| tensor(left, &[Atom::num(i), dummy.clone()]) * &tail)
                .collect::<Vec<_>>(),
        );
        let source = SymbolicTensor::infer(parse_lit!(stream_spectator) * sum).unwrap();
        for mode in [CollectionMode::Monomials, CollectionMode::Factored] {
            let run = |limit| {
                let previous = STREAM_THRESHOLD.replace(limit);
                PEAK_SUMMAND_NODES.set(0);
                PEAK_OPEN_NODES.set(0);
                let value = match mode {
                    CollectionMode::Monomials => {
                        source.collect(TensorCollectFilter::<0>::TaggedTensors)
                    }
                    CollectionMode::Factored => source.collect_with_map(
                        mode,
                        None,
                        |value| TensorCollectFilter::<0>::TaggedTensors.matches(value),
                        |value, _, _| Ok(value),
                    ),
                };
                STREAM_THRESHOLD.set(previous);
                (
                    value.unwrap(),
                    PEAK_SUMMAND_NODES.get(),
                    PEAK_OPEN_NODES.get(),
                )
            };
            let (whole, _, whole_nodes) = run(usize::MAX);
            let (streamed, summand_nodes, other_nodes) = run(0);
            assert_eq!(streamed, whole);
            // No nested opens for these flat two-factor products: this is the
            // peak number of live nodes in the deep planning network itself.
            assert_eq!(summand_nodes, 3);
            assert_eq!(other_nodes, 0);
            assert_eq!(whole_nodes, 1 + 96 * summand_nodes);
        }
    }

    #[test]
    fn streamed_summands_preserve_tape_refusal() {
        test_support::test_initialize();
        let selected = SPENSO_TAG.tensor_symbol("stream_overflow_leaf");
        let term = tensor(selected, &[]);
        let source = SymbolicTensor::infer((&term + parse_lit!(stream_x)).pow(80) + &term).unwrap();
        let AtomView::Add(sum) = source.expression.as_view() else {
            unreachable!()
        };
        let mut input = CollectionInput::new(
            &source,
            |node| node.get_symbol() == Some(selected),
            CollectionMode::Monomials,
        );
        let whole = input
            .open(sum.as_view(), PlanningSource::root(&source.proofs))
            .unwrap();
        let tree = whole.graph.expr_tree().cast::<ChildVecStore<()>>();
        let root = whole.graph.graph.node_id(whole.graph.head());
        let analysis = input.analyze(&whole, &tree, root).unwrap();
        assert!(
            input
                .compile(
                    &mut TermTape::default(),
                    &whole,
                    &tree,
                    root,
                    true,
                    &analysis
                )
                .unwrap()
                .is_none()
        );
        let mut input = CollectionInput::new(
            &source,
            |node| node.get_symbol() == Some(selected),
            CollectionMode::Monomials,
        );
        assert!(
            input
                .compile_sum(&mut TermTape::default(), sum, &source.proofs)
                .unwrap()
                .is_none()
        );
    }

    #[test]
    fn streamed_summands_validate_later_ports_after_refusal() {
        test_support::test_initialize();
        let left = SPENSO_TAG.tensor_symbol("stream_positional_left");
        let right = SPENSO_TAG.tensor_symbol("stream_positional_right");
        let selected = SPENSO_TAG.tensor_symbol("stream_refused_scalar");
        let scalar = tensor(selected, &[]);
        let source = SymbolicTensor::infer(
            (&scalar + parse_lit!(stream_x)).pow(80) * tensor(left, &[mink!(4)])
                + (&scalar + parse_lit!(stream_y)).pow(80) * tensor(right, &[mink!(4)]),
        )
        .unwrap();
        let AtomView::Add(sum) = source.expression.as_view() else {
            panic!("expected a sum in the positional fixture")
        };
        // Both possible first summands overflow the tape. Parsing the whole
        // sum still fails: its two AUTO occurrences allocate distinct labels.
        let owner = PlanningSource::root(&source.proofs);
        let mut input = CollectionInput::new(
            &source,
            |node| node.get_symbol() == Some(selected),
            CollectionMode::Monomials,
        );
        let first = planning::parse_summands(sum, owner.clone())
            .next()
            .unwrap()
            .unwrap();
        let tree = first.graph.expr_tree().cast::<ChildVecStore<()>>();
        let root = first.graph.graph.node_id(first.graph.head());
        let analysis = input.analyze(&first, &tree, root).unwrap();
        assert!(
            input
                .compile(
                    &mut TermTape::default(),
                    &first,
                    &tree,
                    root,
                    false,
                    &analysis
                )
                .unwrap()
                .is_none()
        );
        let whole = input.open(sum.as_view(), owner).unwrap_err();
        let streamed = input
            .compile_sum(&mut TermTape::default(), sum, &source.proofs)
            .err()
            .expect("later summand must be validated");
        for error in [whole, streamed] {
            assert!(matches!(error, TensorInferenceError::Planning(source)
                if matches!(*source, TensorNetworkError::IncompatibleSummand(_))));
        }
    }

    fn assert_streamed_sum_matches_whole(
        source: &SymbolicTensor<PartialStructure>,
        selections: &[Option<Symbol>],
    ) {
        let AtomView::Add(sum) = source.expression.as_view() else {
            panic!("expected a sum in the streaming fixture")
        };
        for mode in [CollectionMode::Monomials, CollectionMode::Factored] {
            // Each individual head selects a different subset of terms.
            // An absent head leaves the entire opened sum opaque.
            for &selection in selections {
                let run = |stream| {
                    let mut input = CollectionInput::new(
                        source,
                        |node| selection.is_some() && node.get_symbol() == selection,
                        mode,
                    );
                    let mut tape = TermTape::default();
                    let compiled = if stream {
                        input
                            .compile_sum(&mut tape, sum, &source.proofs)
                            .unwrap()
                            .and_then(|leaf| tape.push_leaf(leaf))
                    } else {
                        let graph = input
                            .open(sum.as_view(), PlanningSource::root(&source.proofs))
                            .unwrap();
                        let tree = graph.graph.expr_tree().cast::<ChildVecStore<()>>();
                        let root = graph.graph.graph.node_id(graph.graph.head());
                        let analysis = input.analyze(&graph, &tree, root).unwrap();
                        input
                            .compile(&mut tape, &graph, &tree, root, true, &analysis)
                            .unwrap()
                    };
                    let (node, size) = compiled.unwrap();
                    tape.distribute(&mut vec![(node, 1)], &Rational::one())
                        .unwrap();
                    let rows = tape.take_terms();
                    let leaves = input
                        .values
                        .into_iter()
                        .map(|(value, selected)| {
                            let history = value
                                .explicit_index_occurrences()
                                .unwrap()
                                .entries()
                                .collect::<HashMap<_, _>>();
                            let mut regions = Vec::new();
                            value.expression.visitor(&mut |node| {
                                regions.push((
                                    value.proofs.published_interface(node).cloned(),
                                    value.proofs.leaf_interfaces.fact(node),
                                ));
                                true
                            });
                            (
                                value.expression,
                                value.structure,
                                selected,
                                history,
                                regions,
                            )
                        })
                        .collect::<Vec<_>>();
                    (size, rows, leaves)
                };
                INFERENCE_CALLS.set(0);
                assert_eq!(run(true), run(false));
                assert_eq!(INFERENCE_CALLS.get(), 0);
            }
        }
    }

    #[test]
    fn selected_generic_scalar_powers_keep_admitted_leaf_facts() {
        test_support::test_initialize();
        let slot = mink!(4, 71921);
        let left = SPENSO_TAG.tensor_symbol("collection_admitted_power_left");
        let right = SPENSO_TAG.tensor_symbol("collection_admitted_power_right");
        let scalar = Atom::var(symbol!("collection_admitted_power_scalar"));
        let source = SymbolicTensor::infer(
            (tensor(left, slice::from_ref(&slot)) * tensor(right, &[slot]) + scalar).pow(2),
        )
        .unwrap();
        INFERENCE_CALLS.with(|count| count.set(0));
        INFERENCE_VISITS.with(|count| count.set(0));
        let collected = source
            .collect(TensorCollectFilter::<0>::TaggedTensors)
            .unwrap();
        assert!(collected.is_scalar());
        collected.expression.visitor(&mut |value| {
            if value
                .get_symbol()
                .is_some_and(|head| head == left || head == right)
            {
                assert_eq!(collected.admitted_region(value).unwrap().rank(), 1);
                assert!(
                    collected
                        .proofs
                        .leaf_interfaces
                        .get()
                        .unwrap()
                        .lock()
                        .unwrap()
                        .contains_key(value.get_data())
                );
                false
            } else {
                true
            }
        });
        let network = collected
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert!(network.graph.dangling_indices().is_empty());
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
    }

    #[test]
    fn materializing_generic_positional_ports_keeps_raw_order() {
        test_support::test_initialize();
        let head = SPENSO_TAG.tensor_symbol("collection_admitted_positional_ports");
        let source = SymbolicTensor::infer(tensor(head, &[mink!(4), mink!(4)])).unwrap();
        let targets = [AbstractIndex::Normal(71924), AbstractIndex::Normal(71923)];
        INFERENCE_CALLS.with(|count| count.set(0));
        INFERENCE_VISITS.with(|count| count.set(0));
        let materialized = source
            .reindex_interface_ports(&HashMap::from([(0, targets[0]), (1, targets[1])]))
            .unwrap();
        let ports = materialized
            .proofs
            .leaf_interfaces
            .get()
            .unwrap()
            .lock()
            .unwrap();
        let written = &ports
            .get(materialized.expression.as_view().get_data())
            .unwrap()
            .1;
        assert_eq!(written, &vec![mink!(4, 71924), mink!(4, 71923)]);
        drop(ports);
        let network = materialized
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert_eq!(network.graph.dangling_indices().len(), 2);
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
    }

    fn contract_selected(
        value: SymbolicTensor<PartialStructure>,
        _complete: bool,
        _coefficient: &mut SymbolicTensor<PartialStructure>,
    ) -> Result<SymbolicTensor<PartialStructure>> {
        let mut contracted = value.contract_parts(Default::default())?;
        contracted.root.work.get_mut().frontier = contracted.status;
        Ok(contracted.root)
    }

    /// Rank-one Euclidean vectors of dimension two, a scalar, and a selector
    /// for the closed-scope collection regressions.
    struct Components {
        heads: [Symbol; 5],
        indices: [Atom; 3],
        x: Atom,
    }

    impl Components {
        fn new(prefix: &str) -> Self {
            test_support::test_initialize();
            let rep = LibraryRep::from(Euclidean {}).new_rep(2);
            Self {
                heads: [0, 1, 2, 3, 4]
                    .map(|head| SPENSO_TAG.rank_one_tensor_symbol(&format!("{prefix}_v{head}"))),
                indices: [97_411, 97_412, 97_413].map(|index| {
                    rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                        .to_atom()
                }),
                x: Atom::var(symbol!(&format!("{prefix}_x"))),
            }
        }

        fn v(&self, head: usize, index: usize) -> Atom {
            tensor(self.heads[head], slice::from_ref(&self.indices[index]))
        }

        fn dot(&self, left: usize, right: usize, index: usize) -> Atom {
            self.v(left, index) * self.v(right, index)
        }

        fn selects(&self, value: AtomView<'_>) -> bool {
            value
                .get_symbol()
                .is_some_and(|head| head == self.heads[0] || head == self.heads[1])
        }

        /// Explicit indices which `value` exports to a surrounding product: a
        /// sum, a bracket or a closed power binds its own pairs.
        fn exported(&self, value: AtomView<'_>) -> Vec<Atom> {
            match value {
                AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.bracket => {
                    self.exported(function.iter().next().unwrap())
                }
                AtomView::Fun(function) => function.iter().map(|slot| slot.to_owned()).collect(),
                AtomView::Add(sum) => self.exported(sum.iter().next().unwrap()),
                AtomView::Mul(product) => {
                    let written = product
                        .iter()
                        .flat_map(|factor| self.exported(factor))
                        .collect::<Vec<_>>();
                    written
                        .iter()
                        .filter(|slot| written.iter().filter(|other| other == slot).count() == 1)
                        .cloned()
                        .collect()
                }
                _ => Vec::new(),
            }
        }

        /// The component of `value` at the bound index values. A product sums
        /// every name two of its factors export, shadowing outer bindings; an
        /// open square contracts its two copies.
        fn value(&self, value: AtomView<'_>, bound: &HashMap<Atom, usize>) -> i128 {
            let summed = |names: Vec<Atom>, term: &dyn Fn(&HashMap<Atom, usize>) -> i128| {
                (0..1_usize << names.len())
                    .map(|code| {
                        let mut inner = bound.clone();
                        for (bit, name) in names.iter().enumerate() {
                            inner.insert(name.clone(), (code >> bit) & 1);
                        }
                        term(&inner)
                    })
                    .sum::<i128>()
            };
            match value {
                AtomView::Num(number) => match number.get_coeff_view() {
                    CoefficientView::Natural(numerator, 1, 0, 1) => i128::from(numerator),
                    _ => panic!("non-integer coefficient {value}"),
                },
                AtomView::Var(_) => 3,
                AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.bracket => {
                    self.value(function.iter().next().unwrap(), bound)
                }
                AtomView::Fun(function) => {
                    let head = self
                        .heads
                        .iter()
                        .position(|head| *head == function.get_symbol())
                        .unwrap() as i128;
                    let index = bound[&function.iter().next().unwrap().to_owned()] as i128;
                    [2, -3, 5, 7, -2][head as usize] + 11 * index * (head + 1)
                }
                AtomView::Add(sum) => sum.iter().map(|term| self.value(term, bound)).sum(),
                AtomView::Mul(product) => {
                    let written = product
                        .iter()
                        .flat_map(|factor| self.exported(factor))
                        .collect::<Vec<_>>();
                    let mut names = written
                        .iter()
                        .filter(|slot| written.iter().filter(|other| other == slot).count() == 2)
                        .cloned()
                        .collect::<Vec<_>>();
                    names.sort();
                    names.dedup();
                    summed(names, &|inner| {
                        product
                            .iter()
                            .map(|factor| self.value(factor, inner))
                            .product()
                    })
                }
                AtomView::Pow(power) => {
                    let (base, exponent) = power.get_base_exp();
                    let AtomView::Num(number) = exponent else {
                        panic!("symbolic exponent {value}");
                    };
                    let CoefficientView::Natural(exponent, 1, 0, 1) = number.get_coeff_view()
                    else {
                        panic!("non-integer exponent {value}");
                    };
                    let exponent = u32::try_from(exponent).unwrap();
                    let open = self.exported(base);
                    if open.is_empty() {
                        self.value(base, bound).pow(exponent)
                    } else {
                        assert_eq!(exponent, 2);
                        summed(open, &|inner| self.value(base, inner).pow(2))
                    }
                }
            }
        }

        /// Every component of `value` over its free index `free`.
        fn components(&self, value: &Atom, free: usize) -> Vec<i128> {
            (0..2)
                .map(|component| {
                    let bound = HashMap::from([(self.indices[free].clone(), component)]);
                    self.value(value.as_view(), &bound)
                })
                .collect()
        }
    }

    #[test]
    fn collection_keeps_closed_scalar_powers_bracketed() {
        let c = Components::new("bracketed_power");
        let bracket = |value: Atom| spenso::bracket!(value);
        let (i, j) = (0, 1);
        let x = &c.x;
        let sources = [
            // Rows whose selected sum carries bracketed closed powers.
            (x * bracket(c.dot(0, 1, i)).pow(2) + x * bracket(c.dot(0, 1, j)).pow(2)) * c.v(2, j),
            // An unselected bracketed power is a spectator.
            (x + c.dot(0, 1, i)) * bracket(c.dot(2, 3, i)).pow(2) * c.v(4, j),
            // A selected bracketed power beside a factor reusing its name.
            x * bracket(c.dot(0, 1, j)).pow(3) * c.v(2, j),
        ];
        for source in sources {
            let typed = SymbolicTensor::infer(source.clone()).unwrap();
            let expected = c.components(&source, j);
            for mode in [CollectionMode::Monomials, CollectionMode::Factored] {
                let collected = typed
                    .collect_with_map(mode, None, |value| c.selects(value), keep_selected)
                    .unwrap();
                assert_eq!(collected.structure, typed.structure);
                // The result is admissible: no name occurs more than twice
                // in one product scope, and every component is unchanged.
                SymbolicTensor::infer(collected.expression.clone()).unwrap();
                assert_eq!(
                    c.components(&collected.expression, j),
                    expected,
                    "{source} -> {}",
                    collected.expression
                );
            }
        }
        // The factored rows are kept exactly, brackets included.
        let source =
            (x * bracket(c.dot(0, 1, i)).pow(2) + x * bracket(c.dot(0, 1, j)).pow(2)) * c.v(2, j);
        let collected = SymbolicTensor::infer(source.clone())
            .unwrap()
            .collect_with_map(
                CollectionMode::Factored,
                None,
                |value| c.selects(value),
                keep_selected,
            )
            .unwrap();
        assert_eq!(collected.expression, source);
    }

    #[test]
    fn monomial_collection_freshens_square_arms_inside_a_selected_leaf() {
        let c = Components::new("lowered_arms");
        let (i, j, k) = (0, 1, 2);
        let x = &c.x;
        let y = Atom::var(symbol!("lowered_arms_y"));
        // The parser lowers this open-base square into two product arms. The
        // open square `v4(i)^2` inside each arm closes its own pair, so the
        // arms must not share it: (v4·v4)^2, never `v4(i)^4`.
        let square = (c.v(2, j) * c.v(4, i).pow(2) + c.v(0, j)).pow(2);
        let sources = [
            &y + x * &square,
            c.v(1, k) + &square * c.v(3, k),
            &y * c.v(1, k) + x * &square * c.v(3, k) * (x + c.dot(0, 1, i)),
            // Two selected factors: the open square closing `i` in the
            // selected sum meets the selected neighbour `v1(i)`.
            &y + x * (Atom::num(2) * c.v(2, i).pow(2) + c.dot(0, 1, j)) * c.v(1, i) * c.v(4, i),
        ];
        for source in sources {
            let typed = SymbolicTensor::infer(source.clone()).unwrap();
            let expected = c.components(&source, k);
            for mode in [CollectionMode::Monomials, CollectionMode::Factored] {
                let collected = typed
                    .collect_with_map(mode, None, |value| c.selects(value), keep_selected)
                    .unwrap();
                assert_eq!(collected.structure, typed.structure);
                SymbolicTensor::infer(collected.expression.clone()).unwrap_or_else(|error| {
                    panic!("{source} -> {}: {error}", collected.expression)
                });
                collected.expression.visitor(&mut |node| {
                    if let AtomView::Pow(power) = node {
                        assert!(
                            !matches!(power.get_base_exp(), (AtomView::Fun(_), AtomView::Num(n))
                                if n.get_coeff_view() != CoefficientView::Natural(2, 1, 0, 1)),
                            "{source} -> {}",
                            collected.expression
                        );
                    }
                    true
                });
                assert_eq!(
                    c.components(&collected.expression, k),
                    expected,
                    "{source} -> {}",
                    collected.expression
                );
            }
        }
    }

    #[test]
    fn monomial_collection_keeps_open_squares_of_a_selected_scalar_sum_closed() {
        let c = Components::new("scalar_selected_sum");
        let (i, j) = (0, 1);
        let [y, z] =
            ["y", "z"].map(|name| Atom::var(symbol!(&format!("scalar_selected_sum_{name}"))));
        let (x, v2) = (c.x.get_symbol().unwrap(), c.heads[2]);
        // The scalar `x` and `v2` are selected. The selected sum's open
        // square closes `i`, which its unselected neighbours `v1(i) v4(i)`
        // also use: distributing it must not write `v1(i)^3`.
        let source = &y
            + &z * (Atom::num(2) * c.v(1, i).pow(2) + Atom::num(3) * &c.x)
                * c.v(1, i)
                * c.v(4, i)
                * c.v(3, j)
                * c.v(2, j);
        let typed = SymbolicTensor::infer(source.clone()).unwrap();
        let collected = typed
            .collect_with_map(
                CollectionMode::Monomials,
                None,
                |value| {
                    value
                        .get_symbol()
                        .is_some_and(|head| head == x || head == v2)
                },
                keep_selected,
            )
            .unwrap();
        SymbolicTensor::infer(collected.expression.clone())
            .unwrap_or_else(|error| panic!("{source} -> {}: {error}", collected.expression));
        assert_eq!(
            c.components(&collected.expression, j),
            c.components(&source, j),
            "{source} -> {}",
            collected.expression
        );
    }

    #[test]
    fn coalesced_closed_products_are_raised_as_one_scalar() {
        let c = Components::new("coalesced_power");
        let (i, j) = (0, 1);
        // A closed product and a closed open square keep their contraction
        // atomic under the power, as `pow` does.
        for closed in [c.dot(0, 1, i), c.v(4, i).pow(2)] {
            let value = SymbolicTensor::infer(closed.clone()).unwrap();
            let squared = value.collection_power(2).unwrap();
            assert_eq!(squared.expression, spenso::bracket!(closed.clone()).pow(2));
            assert_eq!(
                c.components(&squared.expression, j),
                c.components(&closed, j)
                    .into_iter()
                    .map(|component| component * component)
                    .collect::<Vec<_>>()
            );
        }
        // An open tensor is still raised to its open square.
        let open = SymbolicTensor::infer(c.v(2, j)).unwrap();
        assert_eq!(
            open.collection_power(2).unwrap().expression,
            c.v(2, j).pow(2)
        );
        // An open leaf which also closes a pair cannot be raised without
        // renaming each copy: refuse instead of writing `v4(i)^4`.
        let mixed = SymbolicTensor::infer(c.v(4, i).pow(2) * c.v(2, j)).unwrap();
        assert!(mixed.collection_power(2).is_err());
        assert!(mixed.collection_power(1).is_ok());
        // Pairs closed inside a sum stay in their scope under the power.
        let scoped = SymbolicTensor::infer((&c.x + c.dot(0, 1, i)) * c.v(2, j)).unwrap();
        assert_eq!(
            scoped.collection_power(2).unwrap().expression,
            ((&c.x + c.dot(0, 1, i)) * c.v(2, j)).pow(2)
        );
    }

    #[test]
    fn coalesced_tape_powers_close_scalar_history_and_reject_non_scalar_excess() {
        test_support::test_initialize();
        let slot = LibraryRep::from(Minkowski {})
            .new_rep(4)
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(99601));
        let left = SymbolicTensor::infer(tensor(
            SPENSO_TAG.rank_one_tensor_symbol("collection_power_scope_left"),
            &[slot.to_atom()],
        ))
        .unwrap();
        let right = SymbolicTensor::infer(tensor(
            SPENSO_TAG.rank_one_tensor_symbol("collection_power_scope_right"),
            &[slot.to_atom()],
        ))
        .unwrap();
        let scalar = left.multiply(&right).unwrap();
        assert!(scalar.is_scalar());
        assert!(
            matches!(scalar.expression.as_view(), AtomView::Fun(function)
            if function.get_symbol() == SPENSO_TAG.dot)
        );
        assert_eq!(scalar.explicit_index_occurrences().unwrap().count(slot), 2);

        for (factor, repetitions) in [(&scalar, 2), (&left, 3)] {
            let mut input = CollectionInput::new(factor, |_| true, CollectionMode::Monomials);
            let mut tape = TermTape::<usize>::default();
            let mut copies = Vec::new();
            for _ in 0..repetitions {
                let TermLeaf::Value(position) = input.leaf(factor.clone(), true) else {
                    panic!("a tensor factor must enter the tape as a value");
                };
                copies.push(tape.leaf(position).unwrap());
            }
            assert_eq!(input.values.len(), 1);
            assert_eq!(tape.leaves().len(), 1);
            let (product, _) = tape.group(copies.into_iter(), false).unwrap();
            tape.distribute(&mut vec![(product, 1)], &Rational::one())
                .unwrap();
            let rows = tape.take_terms();
            assert_eq!(rows.len(), 1);
            assert_eq!(rows[0].0.len(), 1);
            let (position, power) = rows[0].0[0];
            assert_eq!(power, repetitions);
            let value = &input.values[tape.leaves()[position]].0;
            let powered = value.collection_power(power);
            if factor.is_scalar() {
                let powered = powered.unwrap();
                assert_eq!(powered.expression, scalar.expression.pow(2));
                assert!(powered.explicit_index_occurrences().unwrap().is_empty());
                let reused = powered.multiply(&left).unwrap();
                assert_eq!(reused.structure, left.structure);
                assert_eq!(reused.explicit_index_occurrences().unwrap().count(slot), 1);
            } else {
                assert!(matches!(
                    powered,
                    Err(TensorInferenceError::Composition(
                        TensorCompositionError::InvalidExplicitMultiplicity { occurrences: 3, .. }
                    ))
                ));
            }
        }

        let free_slot = slot.rep().slot(AbstractIndex::Normal(99603));
        let free = SymbolicTensor::infer(tensor(
            SPENSO_TAG.rank_one_tensor_symbol("collection_power_scope_free"),
            &[free_slot.to_atom()],
        ))
        .unwrap();
        let open = scalar.multiply(&free).unwrap();
        assert_eq!(open.explicit_index_occurrences().unwrap().count(slot), 2);
        assert_eq!(
            open.explicit_index_occurrences().unwrap().count(free_slot),
            1
        );
        let squared = open.collection_power(2).unwrap();
        assert_eq!(squared.explicit_index_occurrences().unwrap().count(slot), 0);
        assert_eq!(
            squared
                .explicit_index_occurrences()
                .unwrap()
                .count(free_slot),
            2
        );
        assert!(squared.multiply(&left).is_ok());
        assert!(matches!(squared.multiply(&free),
            Err(TensorCompositionError::InvalidExplicitMultiplicity { index, occurrences: 3 })
                if index == free_slot.aind()));
    }

    #[test]
    fn vanished_selected_word_retains_terminal_facts_after_coefficient_connections() {
        crate::test_support::test_initialize();
        let a = spenso::mink!(4, 98972);
        let b = spenso::mink!(4, 98973);
        let selected = SymbolicTensor::infer(spenso::p!(&a) * spenso::q!(&b)).unwrap();
        let history = Arc::clone(selected.proofs.occurrences.get().unwrap());
        let zero = SymbolicTensor::from_validated_parts(Atom::Zero, selected.structure, history);
        let coefficient = SymbolicTensor::infer(spenso::p!(&a)).unwrap();
        let result = SymbolicTensor::collection_product(&[zero, coefficient]).unwrap();
        let survivor = SymbolicTensor::infer(spenso::q!(&b)).unwrap();
        let survivor_facts = DomainObservations::terminal_polynomial(&survivor.structure, false);
        let _ = survivor
            .work
            .get_or_init()
            .observations
            .set(Arc::new(survivor_facts));

        assert!(result.expression.is_zero());
        assert_eq!(
            result.structure.logical_slots(),
            survivor.structure.logical_slots()
        );
        assert!(result.reduction_observations().is_terminal_polynomial());
        assert!(
            result
                .reduction_observations()
                .excludes_internal_connections()
        );
        assert!(
            SymbolicTensor::collection_observations(
                &[result, survivor.clone()],
                &survivor.structure,
                true,
            )
            .is_some()
        );
    }

    #[test]
    fn terminal_selected_polynomial_keeps_open_colour_spectator_connection_proof() {
        crate::test_support::test_initialize();
        let ports = [98974, 98975, 98976, 98977].map(|index| spenso::mink!(4, Atom::num(index)));
        assert_eq!(
            ports.iter().collect::<std::collections::HashSet<_>>().len(),
            4
        );
        let [a, b, c, d] = ports;
        let head = spenso::tensor_symbol!("terminal_collection_selected");
        let selected = tensor(head, &[a.clone(), b.clone(), c.clone(), d.clone()]);
        let spectator = crate::color_t!([8, 98978], [3, 98979], [3, 98980]);
        let source = SymbolicTensor::infer(selected * &spectator).unwrap();
        let polynomial =
            spenso::g!(&a, &b) * spenso::g!(&c, &d) + spenso::g!(&a, &c) * spenso::g!(&b, &d);
        let result = source
            .collect_with_map(
                CollectionMode::Monomials,
                None,
                |value| matches!(value, AtomView::Fun(function) if function.get_symbol() == head),
                |selected, _, _| {
                    let facts = DomainObservations::terminal_polynomial(&selected.structure, false);
                    selected.with_identity_result(polynomial.clone(), Some(facts))
                },
            )
            .unwrap();
        assert_eq!(result.expression, polynomial * spectator);
        assert_eq!(
            result.structure.logical_slots(),
            source.structure.logical_slots()
        );
        let observed = result.reduction_observations();
        assert_eq!(observed.candidates.counts[18], 1);
        assert!(observed.excludes_internal_connections());
        assert!(!observed.is_terminal_polynomial());
        let planning_costs = || {
            use super::super::simplification::observation::{INITIAL_SCANS, REPLACEMENT_SCANS};
            (
                INITIAL_SCANS.with(|count| count.get()),
                REPLACEMENT_SCANS.with(|count| count.get()),
                SHALLOW_BUILDS.with(|count| count.get()),
                SELECTED_OPENS.with(|count| count.get()),
            )
        };
        let before = planning_costs();
        let contracted = result.contract_parts(Default::default()).unwrap();
        assert_eq!(
            contracted.status,
            super::super::simplification::ReductionStatus::Complete
        );
        assert_eq!(contracted.root.expression, result.expression);
        assert_eq!(planning_costs(), before);
    }

    #[test]
    fn selected_kernel_can_consume_incident_coefficient_ports() {
        crate::test_support::test_initialize();
        let mu = spenso::mink!(4, 98971);
        let compact = spenso::mink!(4);
        let head = spenso::tensor_symbol!("coefficient_connection_t");
        let t = tensor(head, std::slice::from_ref(&mu));
        let spectator =
            symbolica::parse_lit!((coefficient_connection_x + coefficient_connection_y) ^ 12);
        let source =
            SymbolicTensor::infer(&spectator * (spenso::p!(&mu) * &t + spenso::q!(&mu) * &t))
                .unwrap();
        let result = source
            .collect_with_map(
                CollectionMode::Monomials,
                None,
                |value| matches!(value, AtomView::Fun(function) if function.get_symbol() == head),
                |selected, complete, coefficient| {
                    if !complete {
                        return Ok(selected);
                    }
                    let joined =
                        SymbolicTensor::collection_product(&[selected, coefficient.clone()])?;
                    let result = joined.contract(Default::default())?;
                    *coefficient = SymbolicTensor::from_validated_parts(
                        Atom::one(),
                        PartialStructure::from_logical_slots([]),
                        Arc::default(),
                    );
                    Ok(result)
                },
            )
            .unwrap();
        assert_eq!(result.structure, source.structure);
        assert_eq!(
            result.expression,
            spectator
                * (tensor(head, &[spenso::p!(&compact)]) + tensor(head, &[spenso::q!(&compact)]))
        );
    }

    #[test]
    fn finalizer_combines_selected_states_before_materialization() {
        crate::test_support::test_initialize();
        let port = spenso::mink!(4, 98981);
        let heads = [
            spenso::tensor_symbol!("finalized_collection_a"),
            spenso::tensor_symbol!("finalized_collection_b"),
        ];
        let foreign = tensor(
            spenso::tensor_symbol!("finalized_collection_foreign"),
            &[spenso::euc!(4, 98982)],
        );
        let spectator =
            symbolica::parse_lit!((finalized_collection_x + finalized_collection_y) ^ 12);
        let selected = tensor(heads[0], std::slice::from_ref(&port))
            + Atom::num(2) * tensor(heads[1], &[port]);
        let source = SymbolicTensor::infer(&selected * &foreign * &spectator).unwrap();
        let mut calls = 0;
        let result = source.collect_with_map(CollectionMode::Monomials,
            Some(&mut |states, spectators| {
                calls += 1;
                assert_eq!(states.len(), 2);
                assert!(spectators.iter().any(|value| value.expression == spectator));
                assert!(spectators.iter().any(|value| value.expression == foreign));
                let terms = states.iter().map(|(selected, coefficient)| {
                    SymbolicTensor::collection_product(&[selected.clone(), coefficient.clone()])
                }).collect::<Result<Vec<_>>>()?;
                let combined = terms[0].with_identity_result(
                    Atom::add_many(terms.iter().map(|term| &term.expression)), None,
                )?;
                *states = vec![(combined, SymbolicTensor::from_validated_parts(
                    Atom::one(), PartialStructure::from_logical_slots([]), Arc::default()))];
                Ok(())
            }),
            |value| matches!(value, AtomView::Fun(function) if heads.contains(&function.get_symbol())),
            |selected, _, _| Ok(selected),
        ).unwrap();
        assert_eq!(calls, 1);
        assert_eq!(result.expression, source.expression);
        assert_eq!(
            result.structure.logical_slots(),
            source.structure.logical_slots()
        );
    }

    #[test]
    fn collection_retains_scalar_spectators_and_is_stable() {
        crate::test_support::test_initialize();
        let a = spenso::mink!(4, 97001);
        let t = spenso::tensor_symbol!("selected_collection_t");
        let u = spenso::tensor_symbol!("selected_collection_u");
        let x = Atom::var(symbolica::symbol!("selected_collection_x"));
        let y = Atom::var(symbolica::symbol!("selected_collection_y"));
        let spectator = (&x + &y).pow(3);
        let source = SymbolicTensor::infer(
            &spectator * (tensor(t, std::slice::from_ref(&a)) + tensor(u, &[a])),
        )
        .unwrap();
        let result = source
            .collect(TensorCollectFilter::<0>::TaggedTensors)
            .unwrap();
        assert_eq!(result.expression.expand(), source.expression.expand());
        assert!(
            matches!(result.expression.as_atom_view(), AtomView::Mul(product)
            if product.iter().any(|factor| factor == spectator.as_atom_view()))
        );
        let repeated = result
            .collect(TensorCollectFilter::<0>::TaggedTensors)
            .unwrap();
        assert_eq!(result, repeated);
    }

    #[test]
    fn selected_callback_receives_monomials_and_preserves_foreign_tensors() {
        crate::test_support::test_initialize();
        let a = spenso::mink!(4, 97002);
        let b = spenso::mink!(4, 97003);
        let selected = spenso::tensor_symbol!("selected_collection_owned");
        let foreign = spenso::tensor_symbol!("selected_collection_foreign");
        let x = Atom::var(symbolica::symbol!("selected_collection_weight"));
        let source = SymbolicTensor::infer(
            (tensor(selected, std::slice::from_ref(&a))
                + &x * tensor(selected, std::slice::from_ref(&a)))
                * tensor(foreign, &[b]),
        )
        .unwrap();
        let mut seen = Vec::new();
        let root = source
            .collect_with_map(
                CollectionMode::Monomials,
                None,
                |value| matches!(value, AtomView::Fun(f) if f.get_symbol() == selected),
                |value, _complete, _| {
                    assert!(!matches!(value.expression.as_atom_view(), AtomView::Add(_)));
                    seen.push(value.expression.clone());
                    Ok(value)
                },
            )
            .unwrap();
        assert!(!seen.is_empty());
        assert!(
            seen.iter()
                .all(|value| value == &tensor(selected, std::slice::from_ref(&a)))
        );
        assert_eq!(root.expression.expand(), source.expression.expand());
    }

    #[test]
    fn factored_selected_sum_reuses_callback_leaves_without_replaying_them() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };
        crate::test_support::test_initialize();
        let calls = Arc::new(AtomicUsize::new(0));
        let seen = Arc::clone(&calls);
        let head = spenso::tensor_symbol!(
            "factored_selected_callback",
            norm = move |_, _| {
                seen.fetch_add(1, Ordering::Relaxed);
            }
        );
        let port = spenso::mink!(4, 97003);
        let body = tensor(head, std::slice::from_ref(&port)) + spenso::p!(port);
        let coefficient = Atom::var(symbolica::symbol!("factored_callback_coefficient"));
        let source = SymbolicTensor::infer(&coefficient * &body).unwrap();
        calls.store(0, Ordering::Relaxed);
        let mut steps = 0;
        let result = source
            .collect_with_map(
                CollectionMode::Factored,
                None,
                |value| TensorCollectFilter::<0>::TaggedTensors.matches(value),
                |value, complete, _| {
                    steps += 1;
                    assert!(complete);
                    let AtomView::Add(terms) = body.as_view() else {
                        unreachable!()
                    };
                    assert!(terms.iter().any(|term| term == value.expression.as_view()));
                    Ok(value)
                },
            )
            .unwrap();
        assert_eq!(result, source);
        assert_eq!(steps, 2);
        assert_eq!(calls.load(Ordering::Relaxed), 0);
        assert!(
            source
                .collect_with_map(
                    CollectionMode::Factored,
                    None,
                    |value| TensorCollectFilter::<0>::TaggedTensors.matches(value),
                    |value, _, _| SymbolicTensor::from_declared_parts(Atom::one(), value.structure),
                )
                .is_err()
        );
    }

    #[test]
    fn factored_outer_sum_reuses_settled_terms_without_distributing_products() {
        use crate::tensor::simplification::observation::SettledRegions;
        crate::test_support::test_initialize();
        let port = spenso::mink!(4, 97033);
        let coefficient = symbolica::parse_lit!(factored_region_x + factored_region_y);
        let other_coefficient = symbolica::parse_lit!(factored_region_u + factored_region_v);
        let left = spenso::p!(&port) * &coefficient;
        let replacement = tensor(
            SPENSO_TAG.rank_one_tensor_symbol("factored_region_r"),
            std::slice::from_ref(&port),
        ) * &coefficient;
        let right = spenso::q!(&port) * &other_coefficient;
        let spectator = symbolica::parse_lit!((factored_spectator_x + factored_spectator_y) ^ 17);
        for spectator in [Atom::one(), spectator] {
            let source = SymbolicTensor::infer(&spectator * (&left + &right)).unwrap();
            let mut settled = SettledRegions::default();
            let mut calls = 0;
            let mut current = source.clone();
            for expected_calls in [2, 3] {
                current = current
                    .collect_with_map(
                        CollectionMode::Factored,
                        None,
                        |value| TensorCollectFilter::<0>::TaggedTensors.matches(value),
                        |value, complete, _| {
                            settled.run(value, complete, |value| {
                                calls += 1;
                                assert!(complete);
                                assert!([&left, &replacement, &right].contains(&&value.expression));
                                if value.expression == left {
                                    value.with_identity_result(replacement.clone(), None)
                                } else {
                                    Ok(value)
                                }
                            })
                        },
                    )
                    .unwrap();
                assert_eq!(calls, expected_calls);
                assert_eq!(current.expression, &spectator * (&replacement + &right));
                assert_eq!(current.structure, source.structure);
            }
        }
    }

    #[test]
    fn collection_product_reuses_a_validated_identity_operand() {
        crate::test_support::test_initialize();
        let source = SymbolicTensor::infer(spenso::p!(spenso::mink!(4, 97034))).unwrap();
        let observations = Arc::clone(source.reduction_observations());
        let graph = source.shallow_graph().unwrap();
        let one = SymbolicTensor::from_validated_parts(
            Atom::one(),
            PartialStructure::from_logical_slots([]),
            Arc::default(),
        );
        let result = SymbolicTensor::collection_product(&[one, source.clone()]).unwrap();
        assert_eq!(result, source);
        assert!(Arc::ptr_eq(&observations, result.reduction_observations()));
        assert!(Arc::ptr_eq(&graph, &result.shallow_graph().unwrap()));
    }

    #[test]
    fn selected_domain_cannot_discard_an_established_tensor_port() {
        crate::test_support::test_initialize();
        let a = spenso::mink!(4, 97004);
        let source = SymbolicTensor::infer(tensor(
            spenso::tensor_symbol!("selected_collection_rank"),
            &[a],
        ))
        .unwrap();
        assert!(
            source
                .collect_with_map(
                    CollectionMode::Monomials,
                    None,
                    |_| true,
                    |value, _complete, _| {
                        SymbolicTensor::from_declared_parts(Atom::one(), value.structure)
                    }
                )
                .is_err()
        );
    }

    #[test]
    fn selected_metric_step_checks_callback_induced_rank_loss() {
        crate::test_support::test_initialize();
        let a = spenso::mink!(4, 97005);
        let b = spenso::mink!(4, 97006);
        let target = b.clone();
        let head = spenso::tensor_symbol!(
            "selected_collection_hook",
            norm = move |node, output| {
                if let AtomView::Fun(function) = node
                    && function
                        .iter()
                        .any(|argument| argument == target.as_atom_view())
                {
                    **output = Atom::one();
                }
            }
        );
        let source =
            SymbolicTensor::infer(tensor(ETS.metric, &[a.clone(), b]) * tensor(head, &[a]))
                .unwrap();
        assert!(
            source
                .collect_with_map(
                    CollectionMode::Monomials,
                    None,
                    |value| TensorCollectFilter::<0>::TaggedTensors.matches(value),
                    contract_selected
                )
                .is_err()
        );
    }

    #[test]
    fn selected_dual_metadata_keeps_nested_callback_payload_opaque() {
        use crate::representations::ColorAntiFundamental;
        use spenso::structure::{abstract_index::AIND_SYMBOLS, representation::RepName};
        use std::sync::atomic::{AtomicUsize, Ordering};
        crate::test_support::test_initialize();
        let port = ColorAntiFundamental {}
            .new_rep(3)
            .slot::<AbstractIndex, _>(97113)
            .to_atom();
        let calls = Arc::new(AtomicUsize::new(0));
        let seen = Arc::clone(&calls);
        let hook = symbolica::symbol!(
            "selected_dual_payload_hook",
            norm = move |_, _| {
                seen.fetch_add(1, Ordering::Relaxed);
            }
        );
        let nested = FunctionBuilder::new(AIND_SYMBOLS.dind)
            .add_arg(FunctionBuilder::new(hook).add_arg(&port).finish())
            .finish();
        let metadata = FunctionBuilder::new(SPENSO_TAG.scalar)
            .add_arg(nested)
            .finish();
        let body = SymbolicTensor::infer(tensor(
            spenso::tensor_symbol!("selected_dual_callback_body"),
            &[metadata, port],
        ))
        .unwrap();
        calls.store(0, Ordering::Relaxed);
        let mut steps = 0;
        let root = body
            .collect_with_map(
                CollectionMode::Monomials,
                None,
                |_| true,
                |value, _complete, _| {
                    steps += 1;
                    Ok(value)
                },
            )
            .unwrap();
        assert_eq!(root, body);
        assert_eq!(steps, 1);
        assert_eq!(calls.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn selected_body_is_visible_beside_another_selected_factor() {
        crate::test_support::test_initialize();
        let a = spenso::mink!(4, 97007);
        let b = spenso::mink!(4, 97008);
        let head = SPENSO_TAG.tensor_symbol("selected_collection_adjacent");
        let body = SymbolicTensor::infer(tensor(head, &[a])).unwrap();
        let right = tensor(head, &[b]);
        let source = SymbolicTensor::infer(&body.expression * &right).unwrap();
        let mut seen = Vec::new();
        let root = source
            .collect_with_map(
                CollectionMode::Monomials,
                None,
                |value| matches!(value, AtomView::Fun(f) if f.get_symbol() == head),
                |value, complete, _| {
                    seen.push((value.expression.clone(), complete));
                    Ok(value)
                },
            )
            .unwrap();
        assert!(
            seen.iter()
                .any(|(value, complete)| *complete && value == &(&body.expression * &right))
        );
        assert!(seen.iter().any(|(_, complete)| !complete));
        assert!(
            seen.iter()
                .filter(|(_, complete)| *complete)
                .all(|(value, _)| value == &(&body.expression * &right))
        );
        assert_eq!(root.expression, &body.expression * right);
    }
    #[test]
    fn selected_power_copies_preserve_independent_internal_dummy_scopes() {
        crate::test_support::test_initialize();
        let p = SPENSO_TAG.rank_one_tensor_symbol("selected_power_p");
        let q = SPENSO_TAG.rank_one_tensor_symbol("selected_power_q");
        let a = spenso::mink!(4, 97011);
        let b = spenso::mink!(4, 97012);
        let x = Atom::var(symbol!("selected_power_x"));
        let dot = tensor(
            ETS.metric,
            &[
                tensor(p, &[spenso::mink!(4)]),
                tensor(q, &[spenso::mink!(4)]),
            ],
        );
        let base = tensor(p, std::slice::from_ref(&a)) * tensor(q, std::slice::from_ref(&a)) + &x;
        for (mode, power) in [CollectionMode::Monomials, CollectionMode::Factored]
            .into_iter()
            .flat_map(|mode| [2, 3].map(|power| (mode, power)))
        {
            let source = SymbolicTensor::infer(base.pow(power)).unwrap();
            let root = source
                .collect_with_map(
                    mode,
                    None,
                    |value| TensorCollectFilter::<0>::TaggedTensors.matches(value),
                    contract_selected,
                )
                .unwrap();
            let actual = root;
            assert_eq!(actual.expression.expand(), (&dot + &x).pow(power).expand());
        }
        let internal = LibrarySlot::<AbstractIndex>::try_from(a.as_view()).unwrap();
        let boundary = LibrarySlot::<AbstractIndex>::try_from(b.as_view()).unwrap();
        let open = (tensor(p, std::slice::from_ref(&a))
            * tensor(q, &[a])
            * tensor(p, std::slice::from_ref(&b))
            + &x * tensor(q, &[b]))
        .pow(2);
        let expected = dot.pow(2)
            * tensor(
                ETS.metric,
                &[
                    tensor(p, &[spenso::mink!(4)]),
                    tensor(p, &[spenso::mink!(4)]),
                ],
            )
            + Atom::num(2) * &x * dot.pow(2)
            + x.pow(2)
                * tensor(
                    ETS.metric,
                    &[
                        tensor(q, &[spenso::mink!(4)]),
                        tensor(q, &[spenso::mink!(4)]),
                    ],
                );
        let spectator = Atom::var(symbol!("selected_power_y"));
        for spectator in [Atom::one(), spectator] {
            let source = SymbolicTensor::infer(&spectator * &open).unwrap();
            let root = source
                .collect_with_map(
                    CollectionMode::Monomials,
                    None,
                    |value| TensorCollectFilter::<0>::TaggedTensors.matches(value),
                    contract_selected,
                )
                .unwrap();
            assert_eq!(root.expression.expand(), (&spectator * &expected).expand());
            assert_eq!(root.structure, source.structure);
            let history = root.explicit_index_occurrences().unwrap();
            let source_history = source.explicit_index_occurrences().unwrap();
            assert_eq!(history.count(internal), 0);
            assert_eq!(history.count(boundary), 2);
            assert_eq!(
                history.oriented_entries().collect::<Vec<_>>(),
                source_history.oriented_entries().collect::<Vec<_>>()
            );
            assert_eq!(history.selected_pair, source_history.selected_pair);
        }
    }
    #[test]
    fn selected_trace_power_freshens_slots_and_preserves_opaque_metadata() {
        test_support::test_initialize();
        let a = mink!(4, 97011);
        let b = mink!(4, 97012);
        let internal = LibrarySlot::<AbstractIndex>::try_from(a.as_view()).unwrap();
        let boundary = LibrarySlot::<AbstractIndex>::try_from(b.as_view()).unwrap();
        let left = SPENSO_TAG.tensor_symbol("selected_trace_metadata_left");
        let right = SPENSO_TAG.tensor_symbol("selected_trace_metadata_right");
        let input = Atom::var(SPENSO_TAG.chain_in);
        let output = Atom::var(SPENSO_TAG.chain_out);
        let vector = tensor(
            SPENSO_TAG.rank_one_tensor_symbol("selected_trace_metadata_vector"),
            slice::from_ref(&a),
        );
        let exponent = tensor(
            SPENSO_TAG.rank_one_tensor_symbol("selected_trace_exponent_left"),
            slice::from_ref(&a),
        ) * tensor(
            SPENSO_TAG.rank_one_tensor_symbol("selected_trace_exponent_right"),
            slice::from_ref(&a),
        );
        let metadata = [
            tensor(SPENSO_TAG.scalar, slice::from_ref(&a)),
            tensor(SPENSO_TAG.pure_scalar, slice::from_ref(&a)),
            tensor(symbol!("selected_trace_opaque_metadata"), &[vector]),
            Atom::var(symbol!("selected_trace_metadata_base")).pow(exponent),
        ];
        for metadata in metadata {
            let body = |index: &Atom| {
                tensor(
                    SPENSO_TAG.trace,
                    &[
                        mink!(4),
                        tensor(
                            *SYM,
                            &[
                                tensor(
                                    left,
                                    &[
                                        metadata.clone(),
                                        index.clone(),
                                        input.clone(),
                                        output.clone(),
                                    ],
                                ),
                                tensor(
                                    right,
                                    &[index.clone(), b.clone(), input.clone(), output.clone()],
                                ),
                            ],
                        ),
                    ],
                )
            };
            let one_copy = SymbolicTensor::infer(body(&a)).unwrap();
            let one_history = one_copy.explicit_index_occurrences().unwrap();
            assert_eq!(one_history.count(internal), 2);
            assert_eq!(one_history.count(boundary), 1);
            assert_eq!(one_copy.rank(), 1);
            let source = SymbolicTensor::infer(one_copy.expression.pow(2)).unwrap();
            assert!(source.is_scalar());
            let mut seen = Vec::new();
            let result = source
                .collect_with_map(
                    CollectionMode::Monomials,
                    None,
                    |_| true,
                    |value, _, _| {
                        seen.push(value.clone());
                        Ok(value)
                    },
                )
                .unwrap();
            let mut maximum_copies = 0;
            for value in &seen {
                let mut internals = Vec::new();
                let mut partners = Vec::new();
                value.expression.as_view().visitor(&mut |node| {
                    let AtomView::Fun(function) = node else {
                        return true;
                    };
                    let arguments = function.iter().collect::<Vec<_>>();
                    if function.get_symbol() == left {
                        assert_eq!(arguments.len(), 4);
                        assert_eq!(arguments[0], metadata.as_view());
                        internals
                            .push(LibrarySlot::<AbstractIndex>::try_from(arguments[1]).unwrap());
                        assert_eq!(arguments[2], input.as_view());
                        assert_eq!(arguments[3], output.as_view());
                        false
                    } else if function.get_symbol() == right {
                        assert_eq!(arguments.len(), 4);
                        partners
                            .push(LibrarySlot::<AbstractIndex>::try_from(arguments[0]).unwrap());
                        assert_eq!(arguments[1], b.as_view());
                        assert_eq!(arguments[2], input.as_view());
                        assert_eq!(arguments[3], output.as_view());
                        false
                    } else {
                        true
                    }
                });
                assert!((1..=2).contains(&internals.len()));
                maximum_copies = maximum_copies.max(internals.len());
                let unique = internals.iter().copied().collect::<HashSet<_>>();
                assert_eq!(unique.len(), internals.len());
                assert_eq!(partners.len(), internals.len());
                assert_eq!(unique, partners.into_iter().collect::<HashSet<_>>());
                let mut expected_history = ExplicitIndexOccurrences::default();
                for &fresh in &internals {
                    assert_eq!(fresh.rep(), internal.rep());
                    assert_ne!(fresh.aind(), internal.aind());
                    assert_ne!(fresh.aind(), boundary.aind());
                    let copy_history = one_history
                        .mapped(|slot| {
                            Ok::<_, TensorInferenceError>(if slot == internal {
                                slot.rep().slot(fresh.aind())
                            } else {
                                slot
                            })
                        })
                        .unwrap();
                    expected_history.append(&copy_history);
                }
                assert_eq!(
                    value.expression,
                    Atom::mul_many(internals.iter().map(|slot| body(&slot.to_atom())))
                );
                let expected_interface = if internals.len() == 1 {
                    one_copy.structure.clone()
                } else {
                    source.structure.clone()
                };
                assert_eq!(value.structure, expected_interface);
                let history = value.explicit_index_occurrences().unwrap();
                assert_eq!(history.count(internal), 0);
                assert_eq!(history.count(boundary), internals.len());
                assert_eq!(
                    history.oriented_entries().collect::<Vec<_>>(),
                    expected_history.oriented_entries().collect::<Vec<_>>()
                );
                assert_eq!(history.selected_pair, expected_history.selected_pair);
            }
            assert_eq!(maximum_copies, 2);
            assert_eq!(result.expression, seen.last().unwrap().expression);
            assert_eq!(result.structure, source.structure);
            assert_eq!(
                result
                    .explicit_index_occurrences()
                    .unwrap()
                    .oriented_entries()
                    .collect::<Vec<_>>(),
                seen.last()
                    .unwrap()
                    .explicit_index_occurrences()
                    .unwrap()
                    .oriented_entries()
                    .collect::<Vec<_>>()
            );
            assert_eq!(
                result.explicit_index_occurrences().unwrap().selected_pair,
                seen.last()
                    .unwrap()
                    .explicit_index_occurrences()
                    .unwrap()
                    .selected_pair
            );
            let original_history = source.explicit_index_occurrences().unwrap();
            let published_history = result.explicit_index_occurrences().unwrap();
            let (fresh, _) = published_history
                .entries()
                .find(|(slot, count)| *count == 2 && original_history.count(*slot) == 0)
                .expect("the materialized product retains a fresh contracted pair");
            let third = SymbolicTensor::infer(tensor(
                SPENSO_TAG.rank_one_tensor_symbol("published_power_pair_third"),
                &[fresh.to_atom()],
            ))
            .unwrap();
            assert!(matches!(
                result.multiply(&third),
                Err(TensorCompositionError::InvalidExplicitMultiplicity { occurrences: 3, .. })
            ));
        }
    }
    #[test]
    fn selected_builtin_operand_powers_freshen_all_retained_pairs() {
        test_support::test_initialize();
        let r = mink!(4);
        let a = mink!(4, 97011);
        let b = mink!(4, 97012);
        let c = mink!(4, 97013);
        let internal_a = LibrarySlot::<AbstractIndex>::try_from(a.as_view()).unwrap();
        let internal_c = LibrarySlot::<AbstractIndex>::try_from(c.as_view()).unwrap();
        let boundary = LibrarySlot::<AbstractIndex>::try_from(b.as_view()).unwrap();
        let p = SPENSO_TAG.tensor_symbol("selected_operand_p");
        let q = SPENSO_TAG.tensor_symbol("selected_operand_q");
        let left = SPENSO_TAG.tensor_symbol("selected_operand_left");
        let right = SPENSO_TAG.tensor_symbol("selected_operand_right");
        let vectors = [
            "selected_operand_v",
            "selected_operand_w",
            "selected_operand_x",
            "selected_operand_y",
        ]
        .map(|name| tensor(SPENSO_TAG.rank_one_tensor_symbol(name), slice::from_ref(&r)));
        for builtin in [ETS.metric, *EPSILON_SYMBOL] {
            let body = |a: &Atom, c: &Atom| {
                let consumed = if builtin == ETS.metric {
                    tensor(
                        builtin,
                        &[
                            tensor(p, &[a.clone(), r.clone()]),
                            tensor(q, &[a.clone(), r.clone()]),
                        ],
                    )
                } else {
                    let compound =
                        tensor(p, slice::from_ref(a)) * tensor(q, slice::from_ref(a)) * &vectors[0];
                    tensor(
                        builtin,
                        &[
                            compound,
                            vectors[1].clone(),
                            vectors[2].clone(),
                            vectors[3].clone(),
                        ],
                    )
                };
                tensor(
                    *SYM,
                    &[
                        consumed,
                        tensor(left, slice::from_ref(c)),
                        tensor(right, &[c.clone(), b.clone()]),
                    ],
                )
            };
            let one_copy = SymbolicTensor::infer(body(&a, &c)).unwrap();
            let one_history = one_copy.explicit_index_occurrences().unwrap();
            assert_eq!(one_copy.rank(), 1);
            assert_eq!(one_history.count(internal_a), 2);
            assert_eq!(one_history.count(internal_c), 2);
            assert_eq!(one_history.count(boundary), 1);
            let source = SymbolicTensor::infer(one_copy.expression.pow(2)).unwrap();
            assert!(source.is_scalar());
            let mut seen = Vec::new();
            let result = source
                .collect_with_map(
                    CollectionMode::Monomials,
                    None,
                    |_| true,
                    |value, _, _| {
                        seen.push(value.clone());
                        Ok(value)
                    },
                )
                .unwrap();
            let mut maximum_copies = 0;
            for value in &seen {
                let mut copies = Vec::new();
                value.expression.as_view().visitor(&mut |node| {
                    let AtomView::Fun(projector) = node else {
                        return true;
                    };
                    if projector.get_symbol() != *SYM {
                        return true;
                    }
                    let mut pairs = [Vec::new(), Vec::new()];
                    node.visitor(&mut |factor| {
                        let AtomView::Fun(function) = factor else {
                            return true;
                        };
                        let symbol = function.get_symbol();
                        let pair = if symbol == p || symbol == q {
                            0
                        } else if symbol == left || symbol == right {
                            1
                        } else {
                            return true;
                        };
                        let arguments = function.iter().collect::<Vec<_>>();
                        pairs[pair]
                            .push(LibrarySlot::<AbstractIndex>::try_from(arguments[0]).unwrap());
                        if symbol == right {
                            assert_eq!(arguments[1], b.as_view());
                        }
                        false
                    });
                    let [a, c] = pairs.map(|pair| {
                        assert_eq!(pair.len(), 2);
                        assert_eq!(pair[0], pair[1]);
                        pair[0]
                    });
                    assert_eq!(node, body(&a.to_atom(), &c.to_atom()).as_view());
                    copies.push((a, c));
                    false
                });
                assert!((1..=2).contains(&copies.len()));
                maximum_copies = maximum_copies.max(copies.len());
                let mut fresh = HashSet::new();
                let mut expected_history = ExplicitIndexOccurrences::default();
                for &(a, c) in &copies {
                    for slot in [a, c] {
                        assert_eq!(slot.rep(), internal_a.rep());
                        assert!(![internal_a, internal_c, boundary].contains(&slot));
                        assert!(fresh.insert(slot));
                    }
                    let copy_history = one_history
                        .mapped(|slot| {
                            Ok::<_, TensorInferenceError>(if slot == internal_a {
                                a
                            } else if slot == internal_c {
                                c
                            } else {
                                slot
                            })
                        })
                        .unwrap();
                    expected_history.append(&copy_history);
                }
                assert_eq!(
                    value.expression,
                    Atom::mul_many(copies.iter().map(|(a, c)| body(&a.to_atom(), &c.to_atom())))
                );
                let expected_interface = if copies.len() == 1 {
                    one_copy.structure.clone()
                } else {
                    source.structure.clone()
                };
                assert_eq!(value.structure, expected_interface);
                let history = value.explicit_index_occurrences().unwrap();
                assert_eq!(history.count(internal_a), 0);
                assert_eq!(history.count(internal_c), 0);
                assert_eq!(history.count(boundary), copies.len());
                assert_eq!(
                    history.oriented_entries().collect::<Vec<_>>(),
                    expected_history.oriented_entries().collect::<Vec<_>>()
                );
                assert_eq!(history.selected_pair, expected_history.selected_pair);
            }
            assert_eq!(maximum_copies, 2, "builtin={builtin:?}");
            assert_eq!(result.expression, seen.last().unwrap().expression);
            assert_eq!(result.structure, source.structure);
            assert_eq!(
                result
                    .explicit_index_occurrences()
                    .unwrap()
                    .oriented_entries()
                    .collect::<Vec<_>>(),
                seen.last()
                    .unwrap()
                    .explicit_index_occurrences()
                    .unwrap()
                    .oriented_entries()
                    .collect::<Vec<_>>()
            );
            assert_eq!(
                result.explicit_index_occurrences().unwrap().selected_pair,
                seen.last()
                    .unwrap()
                    .explicit_index_occurrences()
                    .unwrap()
                    .selected_pair
            );
            let original_history = source.explicit_index_occurrences().unwrap();
            let published_history = result.explicit_index_occurrences().unwrap();
            let (fresh, _) = published_history
                .entries()
                .find(|(slot, count)| *count == 2 && original_history.count(*slot) == 0)
                .expect("the materialized product retains a fresh contracted pair");
            let third = SymbolicTensor::infer(tensor(
                SPENSO_TAG.rank_one_tensor_symbol("published_power_pair_third"),
                &[fresh.to_atom()],
            ))
            .unwrap();
            assert!(matches!(
                result.multiply(&third),
                Err(TensorCompositionError::InvalidExplicitMultiplicity { occurrences: 3, .. })
            ));
        }
    }
    #[test]
    fn scalar_selected_sum_does_not_require_dummy_scope_exposure() {
        crate::test_support::test_initialize();
        let hook = symbolica::symbol!("selected_scalar_scope_hook"; Scalar; norm=|_,_| {});
        let selected = tensor(hook, &[Atom::var(symbolica::symbol!("scope_payload"))]);
        assert!(!InterfaceInference::normalization_is_intrinsic(
            selected.as_view()
        ));
        let x = Atom::var(symbolica::symbol!("scope_scalar_spectator"));
        let source = SymbolicTensor::infer(&x * (&selected + Atom::one())).unwrap();
        let mut seen = Vec::new();
        let root = source
            .collect_with_map(
                CollectionMode::Monomials,
                None,
                |value| matches!(value, AtomView::Fun(fun) if fun.get_symbol() == hook),
                |value, _complete, _| {
                    seen.push(value.expression.clone());
                    if value.expression == selected {
                        value.with_checked_expression(Atom::num(2))
                    } else {
                        assert!(value.expression.is_one());
                        Ok(value)
                    }
                },
            )
            .unwrap();
        assert_eq!(seen.iter().filter(|value| **value == selected).count(), 1);
        assert_eq!(root.expression, Atom::num(3) * x);
    }
    #[test]
    fn one_selected_sum_preserves_its_existing_open_scope() {
        crate::test_support::test_initialize();
        let a = spenso::mink!(4, 97021);
        let hook = spenso::tensor_symbol!("selected_open_scope_hook", norm = |_, _| {});
        let target = tensor(
            spenso::tensor_symbol!("selected_open_scope_target"),
            std::slice::from_ref(&a),
        );
        let selected = tensor(hook, &[a]);
        assert!(!InterfaceInference::normalization_is_intrinsic(
            selected.as_view()
        ));
        let x = Atom::var(symbolica::symbol!("open_scope_scalar_spectator"));
        let source = SymbolicTensor::infer(&x * (&selected + &target)).unwrap();
        let mut selected_calls = 0;
        let root = source
            .collect_with_map(
                CollectionMode::Monomials,
                None,
                |value| TensorCollectFilter::<0>::TaggedTensors.matches(value),
                |value, _complete, _| {
                    if value.expression == selected {
                        selected_calls += 1;
                        value.with_checked_expression(target.clone())
                    } else {
                        assert_eq!(value.expression, target);
                        Ok(value)
                    }
                },
            )
            .unwrap();
        assert_eq!(selected_calls, 1);
        assert_eq!(root.expression, Atom::num(2) * x * target);
    }

    #[test]
    fn coefficient_equality_checks_operands_before_forming_their_difference() {
        crate::test_support::test_initialize();
        let x = Atom::var(symbolica::symbol!("coefficient_equality_x"));
        let y = Atom::var(symbolica::symbol!("coefficient_equality_y"));
        let z = Atom::var(symbolica::symbol!("coefficient_equality_z"));
        let filter = TensorCollectFilter::<0>::TaggedTensors;
        let sum = SymbolicTensor::infer(&x + &y).unwrap();
        assert_eq!(
            sum.coefficients_equal(&sum, filter).unwrap(),
            ConditionResult::True
        );
        let left = SymbolicTensor::infer(&x * (&y + &z)).unwrap();
        let right = SymbolicTensor::infer(&x * &y + &x * &z).unwrap();
        assert_eq!(
            left.coefficients_equal(&right, filter).unwrap(),
            ConditionResult::True
        );
        let squared = SymbolicTensor::infer((&x + &y).pow(2)).unwrap();
        let distributed =
            SymbolicTensor::infer(x.pow(2) + Atom::num(2) * &x * &y + y.pow(2)).unwrap();
        assert_eq!(
            squared.coefficients_equal(&distributed, filter).unwrap(),
            ConditionResult::Inconclusive
        );
        let vector = SymbolicTensor::infer(tensor(
            spenso::tensor_symbol!("coefficient_equality_vector"),
            &[spenso::mink!(4, 97203)],
        ))
        .unwrap();
        assert!(vector.coefficients_equal(&sum, filter).is_err());
        assert_eq!(squared.expression, (&x + &y).pow(2));
    }

    #[test]
    fn coefficient_equality_keeps_closed_denominator_dummy_scopes() {
        crate::test_support::test_initialize();
        let [p, q, r, s, u, v] = [
            "coefficient_scope_p",
            "coefficient_scope_q",
            "coefficient_scope_r",
            "coefficient_scope_s",
            "coefficient_scope_u",
            "coefficient_scope_v",
        ]
        .map(|name| {
            tensor(
                SPENSO_TAG.rank_one_tensor_symbol(name),
                &[spenso::mink!(4, 97240)],
            )
        });
        let open = tensor(
            SPENSO_TAG.rank_one_tensor_symbol("coefficient_scope_open"),
            &[spenso::mink!(4, 97241)],
        );
        let first = &p * &q / (Atom::one() + &r * &s);
        let second = &r * &s / (Atom::one() + &u * &v);
        let sum = SymbolicTensor::infer(&open * &first + &open * &second).unwrap();
        let factored = SymbolicTensor::infer(&open * (&first + &second)).unwrap();
        let before = (sum.clone(), factored.clone());
        assert_eq!(
            sum.coefficients_equal(&factored, TensorCollectFilter::<0>::Tensors)
                .unwrap(),
            ConditionResult::True
        );
        assert_eq!((sum, factored), before);
    }

    #[test]
    fn coefficient_zero_proof_resolves_collected_rows_and_preserves_interfaces() {
        crate::test_support::test_initialize();
        let a = spenso::mink!(4, 97201);
        let selected = tensor(spenso::tensor_symbol!("zero_proof_tensor"), &[a]);
        let x = Atom::var(symbolica::symbol!("zero_proof_x"));
        let y = Atom::var(symbolica::symbol!("zero_proof_y"));
        let z = Atom::var(symbolica::symbol!("zero_proof_z"));
        let source = SymbolicTensor::infer(
            &x * (&y + &z) * &selected - &x * &y * &selected - &x * &z * &selected,
        )
        .unwrap();
        assert!(!source.expression.is_zero());
        let before = source.clone();
        assert_eq!(
            source
                .coefficients_are_zero(TensorCollectFilter::<0>::TaggedTensors)
                .unwrap(),
            ConditionResult::True,
        );
        assert_eq!(source, before);
        let nonzero = SymbolicTensor::infer(&x * selected).unwrap();
        assert_eq!(
            nonzero
                .coefficients_are_zero(TensorCollectFilter::<0>::TaggedTensors)
                .unwrap(),
            ConditionResult::False,
        );
    }

    #[test]
    fn coefficient_zero_proof_keeps_distributive_identities_inconclusive() {
        crate::test_support::test_initialize();
        let x = Atom::var(symbolica::symbol!("unproved_zero_x"));
        let y = Atom::var(symbolica::symbol!("unproved_zero_y"));
        let identity = (&x + &y).pow(2) - x.pow(2) - Atom::num(2) * &x * &y - y.pow(2);
        for expression in [
            identity.clone(),
            identity
                * tensor(
                    spenso::tensor_symbol!("unproved_zero_tensor"),
                    &[spenso::mink!(4, 97202)],
                ),
        ] {
            let source = SymbolicTensor::infer(expression).unwrap();
            let before = source.clone();
            assert_eq!(
                source
                    .coefficients_are_zero(TensorCollectFilter::<0>::TaggedTensors)
                    .unwrap(),
                ConditionResult::Inconclusive,
            );
            assert_eq!(source, before);
        }
    }

    #[test]
    fn coefficient_zero_proof_preserves_large_foreign_factors() {
        crate::test_support::test_initialize();
        let t = tensor(spenso::tensor_symbol!("opaque_proof_t"), &[]);
        let spectator = Atom::mul_many((0..30).map(|index| {
            Atom::var(symbolica::symbol!(&format!("zero_spectator_{index}"))) + Atom::one()
        }));
        let source = SymbolicTensor::infer(&spectator * &t).unwrap();
        let mut observed = false;
        source
            .collect_with_map(
                CollectionMode::Monomials,
                Some(&mut |_, spectators| {
                    assert_eq!(spectators.len(), 30);
                    assert_eq!(
                        Atom::mul_many(spectators.iter().map(|value| &value.expression)),
                        spectator,
                    );
                    observed = true;
                    Ok(())
                }),
                |value| TensorCollectFilter::<0>::TaggedTensors.matches(value),
                |value, _, _| Ok(value),
            )
            .unwrap();
        assert!(observed);
        assert_eq!(
            source
                .coefficients_are_zero(TensorCollectFilter::<0>::TaggedTensors)
                .unwrap(),
            // Symbolica proves this product is structurally nonzero without
            // opening its factors, before its disabled sampling/fallback stages.
            ConditionResult::False,
        );
    }

    #[test]
    fn row_only_collection_decodes_auto_spectators_before_publication() {
        test_support::test_initialize();
        let spectator = SymbolicTensor::<PartialStructure>::from_signature(
            &crate::dirac::AGS.gamma_strct::<AbstractIndex>(4),
        )
        .unwrap()
        .reindex_interface_ports(&HashMap::from([(2, AbstractIndex::Normal(99514))]))
        .unwrap();
        let head = SPENSO_TAG.tensor_symbol("collection_auto_factor");
        let selected = SymbolicTensor::infer(tensor(head, &[])).unwrap();
        let source = spectator.multiply(&selected).unwrap();
        let select = |value: AtomView<'_>| value.get_symbol() == Some(head);

        PUBLISHED_TAPE_REGIONS.set(0);
        let eager = source
            .collect_with_map(
                CollectionMode::Monomials,
                Some(&mut |_, spectators| {
                    assert_eq!(spectators.len(), 1);
                    assert_eq!(spectators[0].structure.open_positions().len(), 2);
                    Ok(())
                }),
                select,
                |value, _, _| Ok(value),
            )
            .unwrap();
        let published = PUBLISHED_TAPE_REGIONS.replace(0);
        let private = source
            .collect_rows_with_map(
                CollectionMode::Monomials,
                &mut |_| Ok(()),
                select,
                |value, _, _| Ok(value),
            )
            .unwrap();
        assert!(published >= 2);
        assert_eq!(PUBLISHED_TAPE_REGIONS.get(), published);
        assert_eq!(private, eager);
        assert_eq!(private.structure.open_positions().len(), 2);
    }

    #[test]
    fn encoded_auto_scratch_is_converted_before_certified_publication() {
        crate::test_support::test_initialize();
        let source = SymbolicTensor::<PartialStructure>::from_signature(
            &crate::dirac::AGS.gamma_strct::<AbstractIndex>(4),
        )
        .unwrap()
        .reindex_interface_ports(&HashMap::from([(2, AbstractIndex::Normal(99513))]))
        .unwrap();
        let mut input = CollectionInput::new(&source, |_| true, CollectionMode::Monomials);
        let network = input.parse(&source).unwrap();
        let tree = network.graph.expr_tree().cast::<ChildVecStore<()>>();
        let root = network.graph.graph.node_id(network.graph.head());
        let mut tape = TermTape::<usize>::default();
        let analysis = input.analyze(&network, &tree, root).unwrap();
        assert!(
            input
                .compile(&mut tape, &network, &tree, root, false, &analysis)
                .unwrap()
                .is_some()
        );
        assert_eq!(input.values.len(), 1);
        let scratch = &input.values[0].0;
        assert_eq!(scratch.expression, source.expression);
        assert!(scratch.structure.logical_slots().iter().any(|slot| {
            matches!(
                slot.aind,
                PartialIndex::Explicit(AbstractIndex::Open { .. })
            )
        }));
        assert!(scratch.proofs.is_admitted());
        assert!(scratch.proofs.occurrences.get().is_some());
        // Encoded graph incidences are not the public logical AUTO interface.
        assert!(
            SymbolicTensor::validate_interface(&scratch.expression, &scratch.structure).is_err()
        );

        // This is also the publication boundary for an unprocessed tape factor.
        let body = SymbolicTensor::collection_product(std::slice::from_ref(scratch)).unwrap();
        assert_eq!(body.expression, source.expression);
        assert_eq!(body.structure, source.structure);
        assert_eq!(body.structure.open_positions(), vec![0, 1]);
        assert!(body.proofs.is_admitted());
        SymbolicTensor::validate_interface(&body.expression, &body.structure).unwrap();
        assert_eq!(body, source);
    }

    #[test]
    fn selected_callback_receives_declared_auto_ports_in_logical_order() {
        crate::test_support::test_initialize();
        let source = SymbolicTensor::<PartialStructure>::from_signature(
            &crate::dirac::AGS.gamma_strct::<AbstractIndex>(4),
        )
        .unwrap()
        .reindex_interface_ports(&HashMap::from([(2, AbstractIndex::Normal(99511))]))
        .unwrap();
        let mut calls = 0;
        let result = source
            .collect_with_map(
                CollectionMode::Monomials,
                None,
                |_| true,
                |selected, _complete, _| {
                    calls += 1;
                    assert_eq!(
                        selected.structure.logical_slots(),
                        source.structure.logical_slots()
                    );
                    assert_eq!(selected.structure.open_positions(), vec![0, 1]);
                    Ok(selected)
                },
            )
            .unwrap();
        assert_eq!(calls, 1);
        assert_eq!(result, source);
    }
    #[test]
    fn coefficient_list_keeps_foreign_tensor_and_scalar_spectators_factored() {
        crate::test_support::test_initialize();
        let selected_rep = spenso::mink!(4, 97501);
        let foreign_rep = spenso::euc!(4, 97502);
        let a = tensor(
            spenso::tensor_symbol!("coefficient_output_a"),
            std::slice::from_ref(&selected_rep),
        );
        let b = tensor(
            spenso::tensor_symbol!("coefficient_output_b"),
            &[selected_rep],
        );
        let foreign = tensor(
            spenso::tensor_symbol!("coefficient_output_foreign"),
            &[foreign_rep],
        );
        let x = Atom::var(symbolica::symbol!("coefficient_output_x"));
        let y = Atom::var(symbolica::symbol!("coefficient_output_y"));
        let spectator = (&x + &y).pow(8) * foreign;
        let expression = (&a + &b) * &spectator;
        let source = SymbolicTensor::infer(expression.clone()).unwrap();
        let original = source.clone();
        let rows = source
            .coefficient_list(TensorCollectFilter::Reps([
                spenso::structure::representation::LibraryRep::from(
                    spenso::structure::representation::Minkowski {},
                ),
            ]))
            .unwrap();
        assert_eq!(rows.len(), 2);
        assert!(rows.iter().any(|(selected, _)| selected.expression == a));
        assert!(rows.iter().any(|(selected, _)| selected.expression == b));
        for (_, coefficient) in &rows {
            assert_eq!(coefficient.expression, spectator);
            assert_eq!(coefficient.structure.canonical().order(), 1);
        }
        let rebuilt = Atom::add_many(
            rows.iter()
                .map(|(selected, coefficient)| &selected.expression * &coefficient.expression),
        );
        assert_eq!(rebuilt.collect_factors(), expression.collect_factors());
        assert_eq!(source, original);
    }

    #[test]
    fn capped_kernel_stops_later_rows_and_retains_exact_remaining_factors() {
        crate::test_support::test_initialize();
        let heads = [
            "capped_kernel_frontier_a",
            "capped_kernel_frontier_b",
            "capped_kernel_frontier_c",
            "capped_kernel_frontier_d",
        ]
        .map(|name| spenso::network::tags::SPENSO_TAG.tensor_symbol(name));
        let source = SymbolicTensor::infer(
            (tensor(heads[0], &[]) + tensor(heads[1], &[]))
                * (tensor(heads[2], &[]) + tensor(heads[3], &[])),
        )
        .unwrap();
        let mut calls = 0;
        let mut complete_rows = false;
        let root = source.collect_with_map(CollectionMode::Monomials,
            Some(&mut |_, _| { complete_rows = true; Ok(()) }),
            |value| matches!(value, AtomView::Fun(function) if heads.contains(&function.get_symbol())),
            |mut selected, _, _| {
                calls += 1;
                selected.work.get_mut().frontier = super::super::simplification::ReductionStatus::Capped;
                Ok(selected)
            },
        ).unwrap();
        assert_eq!(calls, 1);
        assert!(!complete_rows);
        assert_eq!(
            root.work.frontier(),
            super::super::simplification::ReductionStatus::Capped
        );
        let result = root;
        assert_eq!(result.expression.expand(), source.expression.expand());
        assert_eq!(result.structure(), source.structure());
    }

    #[test]
    fn coefficient_list_preserves_empty_selection() {
        crate::test_support::test_initialize();
        let x = Atom::var(symbolica::symbol!("coefficient_output_scalar"));
        let source = SymbolicTensor::from_produced_parts(
            (&x + Atom::one()).pow(8),
            PartialStructure::from_logical_slots([]),
        )
        .unwrap();
        let rows = source
            .coefficient_list(TensorCollectFilter::<0>::TaggedTensors)
            .unwrap();
        assert_eq!(rows.len(), 1);
        assert!(rows[0].0.expression.is_one());
        assert_eq!(rows[0].1, source);
    }

    #[test]
    fn coefficient_list_keeps_selected_broadcast_as_one_open_leaf() {
        crate::test_support::test_initialize();
        let port = spenso::mink!(4, 97310);
        let value = tensor(spenso::vector_symbol!("coefficient_broadcast_p"), &[port]);
        let wrapper = tensor(spenso::broadcast_symbol!("coefficient_broadcast"), &[value]);
        let weight =
            (Atom::one() + Atom::var(symbolica::symbol!("coefficient_broadcast_x"))).pow(7);
        let source = SymbolicTensor::infer(&wrapper * &weight).unwrap();
        let rows = source
            .coefficient_list(TensorCollectFilter::Reps([
                spenso::structure::representation::LibraryRep::from(
                    spenso::structure::representation::Minkowski {},
                ),
            ]))
            .unwrap();
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].0.expression, wrapper);
        assert_eq!(rows[0].0.structure, source.structure);
        assert_eq!(rows[0].1.expression, weight);
        assert!(rows[0].1.structure.canonical().is_scalar());
    }

    #[test]
    fn coefficient_list_refuses_callback_inside_selected_broadcast_without_replay() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        crate::test_support::test_initialize();
        let calls = Arc::new(AtomicUsize::new(0));
        let seen = Arc::clone(&calls);
        let hook = spenso::tensor_symbol!(
            "coefficient_broadcast_hook",
            norm = move |_, _| {
                seen.fetch_add(1, Ordering::Relaxed);
            }
        );
        let value = tensor(hook, &[spenso::mink!(4, 97311)]);
        let wrapper = tensor(
            spenso::broadcast_symbol!("coefficient_sensitive_broadcast"),
            &[value],
        );
        let source = SymbolicTensor::infer(wrapper).unwrap();
        calls.store(0, Ordering::Relaxed);
        assert!(
            source
                .coefficient_list(TensorCollectFilter::<0>::TaggedTensors)
                .is_err()
        );
        assert_eq!(calls.load(Ordering::Relaxed), 0);
    }
    #[test]
    fn coefficient_list_keeps_closed_inverse_and_fractional_powers_selected() {
        crate::test_support::test_initialize();
        let p = tensor(
            spenso::vector_symbol!("coefficient_power_p"),
            &[spenso::mink!(4)],
        );
        let q = tensor(
            spenso::vector_symbol!("coefficient_power_q"),
            &[spenso::mink!(4)],
        );
        let dot = tensor(spenso::network::tags::SPENSO_TAG.dot, &[p, q]);
        let open = tensor(
            spenso::vector_symbol!("coefficient_power_open"),
            &[spenso::mink!(4, 97312)],
        );
        let coefficient =
            (Atom::one() + Atom::var(symbolica::symbol!("coefficient_power_x"))).pow(7);
        let filter =
            TensorCollectFilter::Reps([spenso::structure::representation::LibraryRep::from(
                spenso::structure::representation::Minkowski {},
            )]);
        for exponent in [Atom::num(-2), Atom::num((1, 2)), Atom::num((-1, 2))] {
            let power = (Atom::one() + &dot).pow(exponent);
            let selected = &power * &open;
            let source = SymbolicTensor::infer(&selected * &coefficient).unwrap();
            let rows = source.coefficient_list(filter).unwrap();
            assert_eq!(rows.len(), 1);
            assert_eq!(rows[0].0.expression, selected);
            assert_eq!(rows[0].0.structure, source.structure);
            assert_eq!(rows[0].1.expression, coefficient);
        }
    }

    #[test]
    fn coefficient_list_preserves_fixed_selector_powers_and_scalar_spectators() {
        crate::test_support::test_initialize();
        let x = Atom::var(symbolica::symbol!("fixed_coefficient_x"; Scalar));
        let y = Atom::var(symbolica::symbol!("fixed_coefficient_y"; Scalar));
        let spectator = (&x + &y).pow(9)
            * tensor(
                symbolica::symbol!("fixed_coefficient_opaque"; Scalar),
                &[(&x - &y).pow(11)],
            )
            / (&x + Atom::num(3));
        for dimension in [
            Atom::num(4),
            Atom::var(symbolica::symbol!("fixed_coefficient_d")),
        ] {
            let rep = spenso::mink!(
                spenso::structure::dimension::Dimension::try_from(dimension.as_view()).unwrap()
            );
            let p = tensor(
                spenso::vector_symbol!("fixed_coefficient_p"),
                std::slice::from_ref(&rep),
            );
            let q = tensor(
                spenso::vector_symbol!("fixed_coefficient_q"),
                std::slice::from_ref(&rep),
            );
            let ordinary = spenso::dot!(&p, &q);
            for component in 0..4 {
                let basis = tensor(
                    ETS.delta,
                    &[
                        tensor(
                            spenso::structure::abstract_index::AIND_SYMBOLS.cind,
                            &[Atom::num(component)],
                        ),
                        rep.clone(),
                    ],
                );
                let fixed = spenso::dot!(&q, &basis);
                for coefficient in [Atom::one(), Atom::num(2), x.clone(), spectator.clone()] {
                    for (body, expected) in [
                        (fixed.pow(2), vec![(fixed.pow(2), coefficient.clone())]),
                        (fixed.pow(3), vec![(fixed.pow(3), coefficient.clone())]),
                        (
                            ordinary.pow(2),
                            vec![(ordinary.pow(2), coefficient.clone())],
                        ),
                        (
                            (&ordinary - &fixed).pow(2),
                            vec![
                                (ordinary.pow(2), coefficient.clone()),
                                (&ordinary * &fixed, -Atom::num(2) * &coefficient),
                                (fixed.pow(2), coefficient.clone()),
                            ],
                        ),
                    ] {
                        let source = SymbolicTensor::infer(&coefficient * body).unwrap();
                        let before = source.clone();
                        let rows = source
                            .coefficient_list(TensorCollectFilter::<0>::Tensors)
                            .unwrap();
                        assert_eq!(rows.len(), expected.len());
                        assert!(
                            rows.iter()
                                .all(|(key, weight)| key.is_scalar() && weight.is_scalar())
                        );
                        let actual = rows
                            .into_iter()
                            .map(|(key, weight)| (key.expression, weight.expression))
                            .collect::<std::collections::BTreeMap<_, _>>();
                        assert_eq!(
                            actual,
                            expected
                                .into_iter()
                                .collect::<std::collections::BTreeMap<_, _>>()
                        );
                        assert_eq!(source, before);
                    }
                }
                let open = tensor(
                    spenso::tensor_symbol!("fixed_coefficient_open"),
                    &[spenso::mink!(
                        spenso::structure::dimension::Dimension::try_from(dimension.as_view())
                            .unwrap(),
                        97601
                    )],
                );
                let selected = &open * fixed.pow(2);
                let source = SymbolicTensor::infer(&selected * &spectator).unwrap();
                let rows = source
                    .coefficient_list(TensorCollectFilter::<0>::Tensors)
                    .unwrap();
                assert_eq!(rows.len(), 1);
                assert_eq!(rows[0].0.expression, selected);
                assert_eq!(rows[0].0.structure, source.structure);
                assert_eq!(rows[0].1.expression, spectator);
            }
        }
    }

    #[test]
    fn coefficient_list_refuses_inverse_callback_scope_without_replay() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        crate::test_support::test_initialize();
        let calls = Arc::new(AtomicUsize::new(0));
        let seen = Arc::clone(&calls);
        let hook = spenso::tensor_symbol!(
            "coefficient_inverse_hook",
            norm = move |_, _| {
                seen.fetch_add(1, Ordering::Relaxed);
            }
        );
        let value = tensor(hook, &[]);
        let source = SymbolicTensor::infer((value + Atom::one()).pow(-2)).unwrap();
        calls.store(0, Ordering::Relaxed);
        assert!(
            source
                .coefficient_list(TensorCollectFilter::<0>::TaggedTensors)
                .is_err()
        );
        assert_eq!(calls.load(Ordering::Relaxed), 0);
    }

    /// Heads and payloads shared by the closed-scope collection tests.
    struct ClosedScope {
        mu: Atom,
        nu: Atom,
        a: Symbol,
        b: Symbol,
        s: Symbol,
        t: Symbol,
        c: Atom,
        d: Atom,
    }

    impl ClosedScope {
        fn new() -> Self {
            test_support::test_initialize();
            Self {
                mu: mink!(4, 97401),
                nu: mink!(4, 97402),
                a: SPENSO_TAG.rank_one_tensor_symbol("collect_scope_a"),
                b: SPENSO_TAG.rank_one_tensor_symbol("collect_scope_b"),
                s: SPENSO_TAG.rank_one_tensor_symbol("collect_scope_s"),
                t: SPENSO_TAG.rank_one_tensor_symbol("collect_scope_t"),
                c: Atom::var(symbol!("collect_scope_c")),
                d: Atom::var(symbol!("collect_scope_d")),
            }
        }

        fn one(head: Symbol, index: &Atom) -> Atom {
            tensor(head, slice::from_ref(index))
        }

        fn pair(&self, index: &Atom) -> Atom {
            Self::one(self.a, index) * Self::one(self.b, index)
        }

        /// A sum which closes `mu`, beside a spectator which reuses `mu`.
        fn joined(&self) -> Atom {
            (&self.c * self.pair(&self.mu) + &self.d * self.pair(&self.mu))
                * Self::one(self.s, &self.mu)
        }

        fn joined_expected(&self) -> Atom {
            (&self.c + &self.d) * self.pair(&self.nu) * Self::one(self.s, &self.mu)
        }

        fn selects(&self, value: AtomView<'_>) -> bool {
            value
                .get_symbol()
                .is_some_and(|head| head == self.a || head == self.b)
        }

        /// The dummy name which the collected `a` factor carries.
        fn emitted_pair(&self, value: &Atom) -> Atom {
            let mut emitted = None;
            value.as_view().visitor(&mut |node| {
                if let AtomView::Fun(function) = node
                    && function.get_symbol() == self.a
                {
                    emitted = function.iter().next().map(|index| index.to_owned());
                }
                true
            });
            emitted.unwrap()
        }

        fn rejects_third_use(&self, value: &SymbolicTensor<PartialStructure>) -> bool {
            let third =
                SymbolicTensor::infer(Self::one(self.t, &self.emitted_pair(&value.expression)))
                    .unwrap();
            matches!(
                value.multiply(&third),
                Err(TensorCompositionError::InvalidExplicitMultiplicity { occurrences: 3, .. })
            )
        }
    }

    fn keep_selected(
        value: SymbolicTensor<PartialStructure>,
        _complete: bool,
        _coefficient: &mut SymbolicTensor<PartialStructure>,
    ) -> Result<SymbolicTensor<PartialStructure>> {
        Ok(value)
    }

    /// Expanded, canonically relabelled payload; panics unless it re-infers.
    /// Tensor expansion preserves every closed pair before canonical comparison.
    fn canonical_payload(value: &Atom) -> Atom {
        SymbolicTensor::infer(value.clone())
            .unwrap()
            .expanded(None, false)
            .unwrap()
            .canonize()
            .unwrap()
            .expression
    }

    #[test]
    fn monomial_collection_keeps_closed_pairs_beside_selected_neighbours() {
        test_support::test_initialize();
        let rep = LibraryRep::from(Euclidean {}).new_rep(2);
        let [i, j, k] = [97_421, 97_422, 97_423].map(|index| {
            rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                .to_atom()
        });
        let heads = [0, 1, 2, 3, 4]
            .map(|head| SPENSO_TAG.rank_one_tensor_symbol(&format!("selected_neighbour_v{head}")));
        let v = |head: usize, index: &Atom| tensor(heads[head], slice::from_ref(index));
        let x = Atom::var(symbol!("selected_neighbour_x"));
        let selects = |value: AtomView<'_>| {
            value
                .get_symbol()
                .is_some_and(|head| head == heads[0] || head == heads[1])
        };
        let collect = |source: &Atom, mode| {
            let calls = AtomicUsize::new(0);
            let result = SymbolicTensor::infer(source.clone())
                .unwrap()
                .collect_with_map(mode, None, selects, |value, complete, coefficient| {
                    calls.fetch_add(1, Ordering::Relaxed);
                    keep_selected(value, complete, coefficient)
                })
                .unwrap();
            (result, calls.load(Ordering::Relaxed))
        };
        let canonical = |value: &Atom| {
            SymbolicTensor::infer(value.clone())
                .unwrap()
                .expanded(None, false)
                .unwrap()
                .canonize()
                .unwrap()
                .expression
        };
        let renamed = [
            &x * (Atom::num(2) + Atom::num(2) * v(0, &k) * v(4, &k)) * v(0, &k) * v(3, &k),
            (Atom::num(3) * v(0, &k) * v(1, &j) * v(3, &k) + Atom::num(2) * v(4, &j))
                * v(0, &j)
                * v(1, &k),
        ];
        for source in &renamed {
            // The selected sum is opened beside a selected factor using the
            // name of its closed pair. That pair is renamed first, so every
            // row stays admissible and the collected value is unchanged.
            let (result, calls) = collect(source, CollectionMode::Monomials);
            assert!(calls > 0);
            assert_ne!(result.work.frontier(), ReductionStatus::Deferred);
            assert_ne!(&result.expression, source);
            assert_eq!(canonical(&result.expression), canonical(source), "{source}");
            // Factored collection does not distribute and is unaffected.
            let (result, _) = collect(source, CollectionMode::Factored);
            assert_eq!(&result.expression, source);
        }
        // Neighbours which are only spectators, or which use other names, are
        // still collected.
        for source in [
            (Atom::num(2) + v(0, &k) * v(4, &k)) * v(2, &k) * v(3, &k),
            (Atom::num(2) + v(0, &k) * v(4, &k)) * v(1, &i) * v(3, &i),
        ] {
            let (result, calls) = collect(&source, CollectionMode::Monomials);
            assert!(calls > 0, "{source}");
            assert_ne!(result.work.frontier(), ReductionStatus::Deferred);
            SymbolicTensor::infer(result.expression).unwrap();
        }
    }

    #[test]
    fn collected_rows_keep_closed_dummy_pairs_inside_their_sum() {
        let scope = ClosedScope::new();
        let (mu, nu) = (&scope.mu, &scope.nu);
        let rho = mink!(4, 97403);
        let x = SPENSO_TAG.rank_one_tensor_symbol("collect_scope_x");
        let y = SPENSO_TAG.rank_one_tensor_symbol("collect_scope_y");
        let e = Atom::var(symbol!("collect_scope_e"));
        let one = ClosedScope::one;
        let spectator = one(scope.s, mu);
        let selector = |value: AtomView<'_>| scope.selects(value);
        // Only some rows carry the pair; their collected sum still closes it,
        // so nothing is renamed.
        let partial = (&scope.c * scope.pair(mu) + &scope.d * scope.pair(mu) + &e * scope.pair(nu))
            * &spectator;
        let partial_expected =
            ((&scope.c + &scope.d) * scope.pair(mu) + &e * scope.pair(nu)) * &spectator;
        // The collected coefficient itself carries the reused pair `mu`.
        let carried = (one(x, mu) * one(y, mu) * scope.pair(nu)
            + one(x, mu) * one(y, mu) * scope.pair(&rho))
            * &spectator;
        let dot = tensor(
            ETS.metric,
            &[tensor(scope.a, &[mink!(4)]), tensor(scope.b, &[mink!(4)])],
        );
        let carried_expected = Atom::num(2) * one(x, nu) * one(y, nu) * dot * &spectator;
        for mode in [CollectionMode::Monomials, CollectionMode::Factored] {
            let source = SymbolicTensor::infer(scope.joined()).unwrap();
            let root = source
                .collect_with_map(mode, None, selector, keep_selected)
                .unwrap();
            assert_eq!(root.structure, source.structure);
            assert_eq!(
                canonical_payload(&root.expression),
                canonical_payload(&scope.joined_expected())
            );
            let source = SymbolicTensor::infer(partial.clone()).unwrap();
            let root = source
                .collect_with_map(mode, None, selector, keep_selected)
                .unwrap();
            assert_eq!(
                &root.expression,
                if mode == CollectionMode::Monomials {
                    &partial_expected
                } else {
                    &partial
                }
            );
            let source = SymbolicTensor::infer(carried.clone()).unwrap();
            let root = source
                .collect_with_map(mode, None, selector, contract_selected)
                .unwrap();
            assert_eq!(root.structure, source.structure);
            assert_eq!(
                canonical_payload(&root.expression),
                canonical_payload(&carried_expected)
            );
        }
        // Public collection by representation leaves the spectator unselected.
        let colour = LibraryRep::from(Euclidean {}).new_rep(3);
        let i = colour
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(97404))
            .to_atom();
        let u = SPENSO_TAG.tensor_symbol("collect_scope_u");
        let v = SPENSO_TAG.tensor_symbol("collect_scope_v");
        let colour_pair = |index: &Atom| {
            tensor(u, &[index.clone(), i.clone()]) * tensor(v, &[index.clone(), i.clone()])
        };
        let source = SymbolicTensor::infer(
            (&scope.c * colour_pair(mu) + &scope.d * colour_pair(mu)) * &spectator,
        )
        .unwrap();
        let root = source
            .collect(TensorCollectFilter::Reps([LibraryRep::from(Euclidean {})]))
            .unwrap();
        assert_eq!(
            canonical_payload(&root.expression),
            canonical_payload(&((&scope.c + &scope.d) * colour_pair(nu) * &spectator))
        );
    }

    #[test]
    fn nested_closed_sum_beside_its_neighbour_keeps_its_scope() {
        let scope = ClosedScope::new();
        let e = Atom::var(symbol!("collect_scope_e"));
        // The selected sum distributes a product whose own sum closes `mu`
        // beside a neighbour which reuses `mu`.
        let source =
            SymbolicTensor::infer(&e * (scope.joined() + ClosedScope::one(scope.t, &scope.mu)))
                .unwrap();
        for mode in [CollectionMode::Monomials, CollectionMode::Factored] {
            let mut calls = 0;
            let root = source
                .collect_with_map(
                    mode,
                    None,
                    |value| scope.selects(value),
                    |value, _, _| {
                        calls += 1;
                        Ok(value)
                    },
                )
                .unwrap();
            assert_eq!(root.expression, source.expression);
            assert_eq!(root.structure, source.structure);
            SymbolicTensor::infer(root.expression.clone()).unwrap();
            if mode == CollectionMode::Monomials {
                // Monomial rows would join the pair with its neighbour; the
                // collector refuses before any domain callback runs.
                assert_eq!(root.work.frontier(), ReductionStatus::Deferred);
                assert_eq!(calls, 0);
            }
        }
    }

    #[test]
    fn nested_unselected_coefficient_scope_keeps_the_baseline_frontier() {
        let scope = ClosedScope::new();
        let x = SPENSO_TAG.rank_one_tensor_symbol("collect_frontier_x");
        let y = SPENSO_TAG.rank_one_tensor_symbol("collect_frontier_y");
        let e = Atom::var(symbol!("collect_frontier_e"));
        let x_nu = ClosedScope::one(x, &scope.nu);
        let y_nu = ClosedScope::one(y, &scope.nu);
        // The mu pair belongs to an unselected coefficient sum. Selecting x/y
        // must not expose that pair beside the separate free use in s(mu).
        let expression = &e
            * ((&x_nu * &y_nu * (&scope.c + scope.pair(&scope.mu)) + &scope.d)
                * ClosedScope::one(scope.s, &scope.mu)
                + ClosedScope::one(scope.t, &scope.mu));
        let mut failures = Vec::new();
        for depth_limit in [None, Some(1), Some(2), Some(3)] {
            for mode in [CollectionMode::Monomials, CollectionMode::Factored] {
                for exact_selector in [false, true] {
                    let monomials = mode == CollectionMode::Monomials;
                    let row = format!(
                        "depth={depth_limit:?} monomials={monomials} exact={exact_selector}"
                    );
                    let mut source = SymbolicTensor::infer(expression.clone()).unwrap();
                    // Cache the shallow graph at this parse depth.
                    let settings = ParseSettings {
                        precontract_scalars: false,
                        depth_limit,
                        shorthand_parsing: ShorthandParsing::Opaque,
                        parse_composite_scalars_as_tensors: true,
                        ..ParseSettings::default()
                    };
                    type Tensor<'a> = SymbolicTensor<OrderedStructure, AtomOrView<'a>>;
                    let network = SymbolicNet::try_from_view::<OrderedStructure, _>(
                        source.expression.as_view(),
                        &DummyLibrary::<Tensor<'_>>::new(),
                        &settings,
                    )
                    .unwrap();
                    let network: SymbolicNet<AbstractIndex> = network.map_ref(
                        |scalar| scalar.as_view().to_owned(),
                        |tensor| SymbolicTensor {
                            expression: tensor.expression.as_view().to_owned(),
                            structure: tensor.structure.clone(),
                            is_metric: tensor.is_metric,
                            is_composite: tensor.is_composite,
                            work: Default::default(),
                            proofs: tensor.proofs.clone(),
                        },
                    );
                    source.work.get_mut().shallow_graph =
                        OnceLock::from(Ok(Arc::new(PlanningNetwork::from_materialized(network))));
                    let mut calls = 0;
                    let result = source.collect_with_map(
                        mode,
                        None,
                        |value| {
                            if exact_selector {
                                value == x_nu.as_view() || value == y_nu.as_view()
                            } else {
                                value
                                    .get_symbol()
                                    .is_some_and(|head| head == x || head == y)
                            }
                        },
                        |value, _, _| {
                            calls += 1;
                            Ok(value)
                        },
                    );
                    // PR-head baseline: Monomials distributes three rows,
                    // Factored keeps the input with two callbacks.
                    match result {
                        Ok(output) => {
                            if SymbolicTensor::infer(output.expression.clone()).is_err()
                                || output.structure != source.structure
                                || output.work.frontier() != ReductionStatus::Complete
                                || calls != if monomials { 3 } else { 2 }
                                || (output.expression == source.expression) == monomials
                            {
                                failures.push(format!(
                                    "{row}: frontier={:?} calls={calls} payload={}",
                                    output.work.frontier(),
                                    output.expression
                                ));
                            }
                        }
                        Err(error) => failures.push(format!("{row}: {error:?}")),
                    }
                }
            }
        }
        assert!(failures.is_empty(), "{failures:#?}");
    }

    #[test]
    fn unselected_collision_is_returned_unchanged_without_mapper_calls() {
        let scope = ClosedScope::new();
        let source = SymbolicTensor::infer(scope.joined()).unwrap();
        for mode in [CollectionMode::Monomials, CollectionMode::Factored] {
            let mut calls = 0;
            let root = source
                .collect_with_map(
                    mode,
                    None,
                    |_| false,
                    |value, _, _| {
                        calls += 1;
                        Ok(value)
                    },
                )
                .unwrap();
            assert_eq!(calls, 0);
            assert_eq!(root.expression, source.expression);
            assert_eq!(root.structure, source.structure);
            assert_eq!(
                root.explicit_index_occurrences()
                    .unwrap()
                    .entries()
                    .collect::<Vec<_>>(),
                source
                    .explicit_index_occurrences()
                    .unwrap()
                    .entries()
                    .collect::<Vec<_>>()
            );
        }
    }

    #[test]
    fn exact_selectors_match_original_names_beside_a_collision() {
        let scope = ClosedScope::new();
        let source = SymbolicTensor::infer(scope.joined()).unwrap();
        let a = ClosedScope::one(scope.a, &scope.mu);
        let b = ClosedScope::one(scope.b, &scope.mu);
        let mu = scope.mu.as_view();
        let by_atom = |value: AtomView<'_>| value == a.as_view() || value == b.as_view();
        // Select every tensor which uses the reused index, except the spectator.
        let by_index = |value: AtomView<'_>| {
            matches!(value, AtomView::Fun(function)
                if function.get_symbol() != scope.s && function.iter().any(|index| index == mu))
        };
        for root in [
            source.collect_with_map(CollectionMode::Monomials, None, by_atom, keep_selected),
            source.collect_with_map(CollectionMode::Monomials, None, by_index, keep_selected),
        ] {
            assert_eq!(
                canonical_payload(&root.unwrap().expression),
                canonical_payload(&scope.joined_expected())
            );
        }
    }

    #[test]
    fn collected_root_pairs_reject_a_third_occurrence_after_recollection() {
        let scope = ClosedScope::new();
        let selector = |value: AtomView<'_>| scope.selects(value);
        // With and without a spectator which forces a fresh dummy name.
        let alone = &scope.c * scope.pair(&scope.mu) + &scope.d * scope.pair(&scope.mu);
        for source in [scope.joined(), alone] {
            let source = SymbolicTensor::infer(source).unwrap();
            let once = source
                .collect_with_map(CollectionMode::Monomials, None, selector, keep_selected)
                .unwrap();
            assert!(!matches!(once.expression.as_view(), AtomView::Add(_)));
            assert!(scope.rejects_third_use(&once));
            let twice = once
                .collect_with_map(CollectionMode::Monomials, None, selector, keep_selected)
                .unwrap();
            assert!(scope.rejects_third_use(&twice));
        }
    }

    #[test]
    fn refused_inverse_scope_beside_a_collision_does_not_replay_callbacks() {
        let scope = ClosedScope::new();
        let calls = Arc::new(AtomicUsize::new(0));
        let seen = Arc::clone(&calls);
        let hook = spenso::tensor_symbol!(
            "collect_scope_inverse_hook",
            norm = move |_, _| {
                seen.fetch_add(1, Ordering::Relaxed);
            }
        );
        let inverse = (tensor(hook, &[]) + Atom::one()).pow(-2);
        let source = SymbolicTensor::infer(
            (&scope.c * scope.pair(&scope.mu) * inverse + &scope.d * scope.pair(&scope.mu))
                * ClosedScope::one(scope.s, &scope.mu),
        )
        .unwrap();
        calls.store(0, Ordering::Relaxed);
        // Monomial collection opens the selected sum before it reaches the
        // closed inverse scope, which it refuses without rebuilding it.
        let root = source
            .collect_with_map(
                CollectionMode::Monomials,
                None,
                |value| {
                    scope.selects(value)
                        || matches!(value, AtomView::Fun(function) if function.get_symbol() == hook)
                },
                keep_selected,
            )
            .unwrap();
        assert_eq!(root.work.frontier(), ReductionStatus::Deferred);
        assert_eq!(root.expression, source.expression);
        assert_eq!(calls.load(Ordering::Relaxed), 0);
    }
}

#[cfg(test)]
mod planning_tests {
    use super::*;
    use spenso::{g, mink};
    use symbolica::atom::FunctionBuilder;

    #[test]
    fn row_only_collection_keeps_exact_spectators_in_the_plan() {
        crate::test_support::test_initialize();
        let a = mink!(4, 98903);
        let b = mink!(4, 98904);
        let spectator = FunctionBuilder::new(spenso::tensor_symbol!("planning_spectator"))
            .add_arg(&a)
            .finish();
        let scalar = symbolica::parse_lit!(planning_scalar);
        let source = SymbolicTensor::infer(&scalar * &spectator * g!(&b, &a)).unwrap();
        let select = |value: AtomView<'_>| TensorCollectFilter::<0>::Metrics.matches(value);

        PUBLISHED_TAPE_REGIONS.set(0);
        let independent = source
            .collect_with_map(
                CollectionMode::Monomials,
                Some(&mut |_, spectators| {
                    assert!(spectators.iter().any(|value| value.expression == spectator));
                    assert!(spectators.iter().any(|value| value.expression == scalar));
                    Ok(())
                }),
                select,
                |value, _, _| Ok(value),
            )
            .unwrap();
        let including_spectators = PUBLISHED_TAPE_REGIONS.replace(0);
        let private = source
            .collect_rows_with_map(
                CollectionMode::Monomials,
                &mut |_| Ok(()),
                select,
                |value, _, _| Ok(value),
            )
            .unwrap();
        assert_eq!(including_spectators - PUBLISHED_TAPE_REGIONS.get(), 2);
        assert_eq!(private, independent);
        assert_eq!(
            private.proofs.published_interface(scalar.as_view()),
            Some(&PartialStructure::from_logical_slots([])),
        );
        assert_eq!(
            private.proofs.leaf_interfaces.fact(spectator.as_view()),
            Some((true, vec![a])),
        );
    }

    #[test]
    fn repeated_opaque_tape_occurrences_publish_once() {
        crate::test_support::test_initialize();
        let a = mink!(4, 98905);
        let atom = FunctionBuilder::new(spenso::tensor_symbol!("planning_repeated_leaf"))
            .add_arg(&a)
            .finish();
        let source = SymbolicTensor::infer(atom).unwrap();
        let mut input =
            CollectionInput::new(&source, |_: AtomView<'_>| true, CollectionMode::Monomials);
        let network = input.parse(&source).unwrap();
        let tree = network.graph.expr_tree().cast::<ChildVecStore<()>>();
        let root = network.graph.graph.node_id(network.graph.head());
        let analysis = input.analyze(&network, &tree, root).unwrap();

        PUBLISHED_TAPE_REGIONS.set(0);
        for _ in 0..2 {
            let mut tape = TermTape::<usize>::default();
            assert!(
                input
                    .compile(&mut tape, &network, &tree, root, false, &analysis)
                    .unwrap()
                    .is_some()
            );
        }
        assert_eq!(PUBLISHED_TAPE_REGIONS.get(), 1);
        assert_eq!(input.values.len(), 1);
    }

    #[test]
    fn unselected_arithmetic_leaves_remain_unopened() {
        crate::test_support::test_initialize();
        let vector = spenso::vector_symbol!("shallow_collection_vector");
        let a = mink!(4, 98901);
        let b = mink!(4, 98902);
        let spectator = symbolica::parse_lit!((shallow_spectator_x + shallow_spectator_y) ^ 1000);
        let expression =
            &spectator * g!(&a, &b) * FunctionBuilder::new(vector).add_arg(&a).finish();
        let source = SymbolicTensor::infer(expression.clone()).unwrap();
        SHALLOW_BUILDS.with(|count| count.set(0));
        SHALLOW_VIEWS.with(|count| count.set(0));
        SELECTED_OPENS.with(|count| count.set(0));
        ANALYZED_NODES.with(|count| count.set(0));
        let collected = source.collect(TensorCollectFilter::<0>::Metrics).unwrap();
        assert_eq!(collected.expression, expression);
        assert_eq!(SHALLOW_VIEWS.with(|count| count.get()), 1);
        assert_eq!(SELECTED_OPENS.with(|count| count.get()), 0);
        let graph = source.shallow_graph().unwrap();
        let tree = graph.graph.expr_tree().cast::<ChildVecStore<()>>();
        let root = graph.graph.graph.node_id(graph.graph.head());
        let nodes = tree
            .iter_preorder_tree_nodes(&graph.graph.graph, root)
            .count();
        assert_eq!(ANALYZED_NODES.with(|count| count.get()), nodes);
        source.collect(TensorCollectFilter::<0>::Metrics).unwrap();
        assert_eq!(ANALYZED_NODES.with(|count| count.get()), 2 * nodes);
        assert!(Arc::ptr_eq(&graph, &source.shallow_graph().unwrap()));
        assert_eq!(SHALLOW_BUILDS.with(|count| count.get()), 1);
    }
}
