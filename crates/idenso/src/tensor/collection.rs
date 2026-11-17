//! Selected algebra on the existing parsed graph and arithmetic tape.

use std::collections::{HashMap, HashSet};
use std::sync::{Arc, OnceLock};

use linnet::half_edge::NodeIndex;
use linnet::half_edge::involution::Hedge;
use linnet::half_edge::subgraph::{ModifySubSet, SuBitGraph};
use linnet::half_edge::tree::SimpleTraversalTree;
use linnet::tree::child_vec::ChildVecStore;
use spenso::network::TensorNetworkError;
use spenso::network::graph::{NetworkEdge, NetworkLeaf, NetworkNode, NetworkOp};
use spenso::network::library::DummyLibrary;
use spenso::network::library::function_lib::INBUILTS;
use spenso::network::parsing::{ParseSettings, ParseState, ShorthandParsing};
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
use symbolica::atom::{Atom, AtomCore, AtomOrView, AtomType, AtomView};
use symbolica::domains::rational::Rational;
use symbolica::id::ConditionResult;

use super::composition::{self, ExplicitIndexOccurrences};
use super::inference::{InterfaceInference, LeafInference, TensorInferenceError};
use super::simplification::observation::DomainObservations;
use super::{SymbolicNet, SymbolicTensor};

mod collected;
pub(crate) use collected::{CollectedOutput, CollectedPolynomial};

type Result<T> = std::result::Result<T, TensorInferenceError>;

#[cfg(test)]
thread_local! {
    pub(crate) static SHALLOW_BUILDS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    pub(crate) static SHALLOW_VIEWS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    pub(crate) static SELECTED_OPENS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
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

type CollectionFinalizer<'a> = dyn FnMut(
        &mut Vec<(
            SymbolicTensor<PartialStructure>,
            SymbolicTensor<PartialStructure>,
        )>,
        &[SymbolicTensor<PartialStructure>],
    ) -> Result<()>
    + 'a;

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
    ) -> Result<Arc<SymbolicNet<AbstractIndex>>> {
        #[cfg(test)]
        SHALLOW_VIEWS.with(|count| count.set(count.get() + 1));
        source.shallow_graph()
    }

    fn open<'src>(
        &self,
        source: AtomView<'src>,
    ) -> Result<SymbolicNet<AbstractIndex, AtomOrView<'src>>> {
        #[cfg(test)]
        SELECTED_OPENS.with(|count| count.set(count.get() + 1));
        let settings = ParseSettings {
            precontract_scalars: false,
            depth_limit: None,
            shorthand_parsing: ShorthandParsing::Opaque,
            parse_composite_scalars_as_tensors: true,
            ..ParseSettings::default()
        };
        type Tensor<'src> = SymbolicTensor<OrderedStructure, AtomOrView<'src>>;
        SymbolicNet::try_from_view::<OrderedStructure, _>(
            source,
            &DummyLibrary::<Tensor<'_>>::new(),
            &settings,
        )
        .map_err(|error| TensorInferenceError::invalid(error.to_string()))
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
        network: &SymbolicNet<AbstractIndex, E>,
        tree: &CollectionTree,
        node: NodeIndex,
    ) -> Result<Atom> {
        network
            .graph
            .to_expression_at(tree, node, &mut |_, value| match value {
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
            })
            .map_err(|error| TensorInferenceError::invalid(error.to_string()))
    }

    fn analyze<E: AtomCore>(
        &mut self,
        network: &SymbolicNet<AbstractIndex, E>,
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
                NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => network.store.tensors[*index]
                    .proofs
                    .occurrences
                    .get()
                    .map(|history| history.as_ref().clone())
                    .ok_or_else(|| {
                        TensorInferenceError::invalid(
                            "parsed tensor is missing construction history",
                        )
                    })?,
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
        if let Ok(number) = Rational::try_from(self.values[position].0.expression.as_atom_view()) {
            TermLeaf::Number(position, number)
        } else {
            TermLeaf::Value(position)
        }
    }

    fn scope_is_intrinsic<E: AtomCore>(
        &self,
        network: &SymbolicNet<AbstractIndex, E>,
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

    fn expose_scope<E: AtomCore + Clone + From<E::Output>>(
        &mut self,
        network: &SymbolicNet<AbstractIndex, E>,
        tree: &CollectionTree,
        node: NodeIndex,
    ) -> Result<Option<SymbolicNet<AbstractIndex, Atom>>> {
        // Prove this graph scope intrinsic before emitting a temporary body:
        // function reconstruction can itself invoke a registered normalizer.
        if !self.scope_is_intrinsic(network, tree, node) {
            return Ok(None);
        }
        if let NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) = &network.graph.graph[node] {
            let tensor = &network.store.tensors[*index];
            let expression = tensor.expression.as_atom_view();
            // Shallow storage does not expose a closed broadcast's arguments.
            // Retained internal pairs are not permission to rename that payload.
            if tensor.structure.order() == 0
                && matches!(expression, AtomView::Fun(function)
                    if function.get_symbol().has_tag(&SPENSO_TAG.broadcast)
                        && function.get_symbol() != INBUILTS.conj)
            {
                return Ok(None);
            }
            if matches!(expression, AtomView::Add(_) | AtomView::Mul(_)) {
                // Expose the selected arithmetic before assigning this copy's
                // internal labels; an opaque leaf has only its public ports.
                let opened = self.open(expression)?;
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
                let history = tensor.proofs.occurrences.get().ok_or_else(|| {
                    TensorInferenceError::invalid("parsed tensor is missing construction history")
                })?;
                for (slot, _) in history.entries() {
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
        let mut graph = network.map_ref(Clone::clone, Clone::clone);
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
        let mapped = graph.map_occurrences(
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
                let history = tensor
                    .proofs
                    .occurrences
                    .get()
                    .ok_or_else(|| fail("parsed tensor is missing construction history"))?;
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
                let (expression, history) =
                    if internal.is_empty() {
                        (
                            tensor.expression.as_atom_view().to_owned(),
                            Arc::clone(history),
                        )
                    } else {
                        // Projectors and other opaque leaves can close pairs before
                        // their public storage ports are built. Rename these known
                        // retained identities in the payload and history together,
                        // including historical uses whose syntax has vanished.
                        let mut parents = Vec::new();
                        let mut slots = SlotMatcher::default();
                        let expression = tensor.expression.as_atom_view().replace_map(
                            |value, context, output| {
                                // Only function descent increments function_level;
                                // arithmetic keeps the same enclosing function stack.
                                parents.truncate(context.function_level);
                                if context.parent_type == Some(AtomType::Fun)
                                    && !InterfaceInference::explicit_index_argument(
                                        parents[context.function_level - 1],
                                        value,
                                        &mut slots,
                                    )
                                {
                                    output.set_from_view(&value);
                                    return;
                                }
                                if let Ok(slot) = LibrarySlot::<AbstractIndex>::try_from(value) {
                                    let key =
                                        slot.rep().base().slot::<AbstractIndex, _>(slot.aind());
                                    if let Some(&index) = internal.get(&key) {
                                        **output =
                                            slot.rep().slot::<AbstractIndex, _>(index).to_atom();
                                    } else {
                                        output.set_from_view(&value);
                                    }
                                } else if Representation::<LibraryRep>::try_from(value).is_ok() {
                                    // Representations own their dimension metadata.
                                    output.set_from_view(&value);
                                } else if let AtomView::Fun(function) = value {
                                    let symbol = function.get_symbol();
                                    if SPENSO_TAG.is_scalar_metadata(symbol) {
                                        output.set_from_view(&value);
                                    } else {
                                        parents.push(symbol);
                                    }
                                }
                            },
                        );
                        let history = history
                            .mapped(|slot| {
                                let key = slot.rep().base().slot::<AbstractIndex, _>(slot.aind());
                                Ok::<_, TensorInferenceError>(
                                    slot.rep()
                                        .slot(internal.get(&key).copied().unwrap_or(slot.aind())),
                                )
                            })
                            .map_err(|error| fail(&error.to_string()))?;
                        (expression, Arc::new(history))
                    };
                let source_tensor = SymbolicTensor::from_validated_parts(
                    expression,
                    source_interface,
                    Arc::clone(&history),
                );
                let expression = composition::rewrite_interface_ports(
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
                if SymbolicTensor::validate_observed_interface(
                    &expression,
                    &expected,
                    LeafInference::ObserveStorage,
                )
                .is_err()
                {
                    unavailable = true;
                    return Err(fail("normalized occurrence changed its interface"));
                }
                let history = history
                    .mapped(|slot| {
                        Ok::<_, TensorInferenceError>(
                            source
                                .iter()
                                .position(|source| *source == slot)
                                .map_or(slot, |position| targets[position]),
                        )
                    })
                    .map_err(|error| fail(&error.to_string()))?;
                let mut proofs = tensor.proofs.algebra_result();
                proofs.occurrences = OnceLock::from(Arc::new(history));
                Ok(SymbolicTensor {
                    expression,
                    structure: OrderedStructure::new(targets).into_canonical(),
                    is_metric: tensor.is_metric,
                    is_composite: tensor.is_composite,
                    proofs,
                })
            },
        );
        match mapped {
            Ok(mapped) => Ok(Some(mapped)),
            Err(_) if unavailable => Ok(None),
            Err(error) => Err(TensorInferenceError::invalid(error.to_string())),
        }
    }

    /// Compile parser-lowered power arms after each copy owns its internal labels.
    fn compile_lowered_power<E: AtomCore + Clone + From<E::Output>>(
        &mut self,
        tape: &mut TermTape<usize>,
        network: &SymbolicNet<AbstractIndex, E>,
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

    fn compile<E: AtomCore + Clone + From<E::Output>>(
        &mut self,
        tape: &mut TermTape<usize>,
        network: &SymbolicNet<AbstractIndex, E>,
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
            if owned[&node]
                && !matches!(value, NetworkNode::Leaf(_))
                && scoped_power.is_none()
                && !function
                && !factored
                && !scalar_power
            {
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
                    // An unknown closed broadcast owns its argument scopes.
                    // Native conjugation is the explicitly supported exception;
                    // copying another function must not rename its closed payload.
                    if interfaces[&base].canonical().is_scalar()
                        && matches!(&network.graph.graph[base], NetworkNode::Op(NetworkOp::Function(symbol))
                            if *symbol != INBUILTS.conj)
                    {
                        refused = true;
                        return Ok(None);
                    }
                    if matches!(network.graph.graph[base], NetworkNode::Leaf(_)) {
                        let base_value = self.emit(network, tree, base)?;
                        let closed_conjugation = matches!(base_value.as_view(), AtomView::Fun(function)
                            if function.get_symbol() == INBUILTS.conj && function.get_nargs() == 1)
                            && interfaces[&base].canonical().is_scalar()
                            && !base_value.get_all_symbols(true).iter().any(|symbol| {
                                symbol.has_tag(&SPENSO_TAG.tensor_word)
                            });
                        if closed_conjugation && !self.scope_is_intrinsic(network, tree, base) {
                            refused = true;
                            return Ok(None);
                        }
                        if matches!(base_value.as_view(), AtomView::Add(_) | AtomView::Mul(_))
                            || matches!(base_value.as_view(), AtomView::Fun(function)
                                if function.get_symbol() == SPENSO_TAG.bracket)
                            || closed_conjugation
                        {
                            // Open the selected power before copying: internal
                            // dummy pairs belong to each copy, not its opaque leaf.
                            // A closed native conjugation can hide that pair at
                            // shallow depth; its intrinsic graph supplies the same
                            // copy witness without opening arbitrary functions.
                            let value = self.emit(network, tree, node)?;
                            let opened = self.open(value.as_view())?;
                            let tree = opened.graph.expr_tree().cast::<ChildVecStore<()>>();
                            let root = opened.graph.graph.node_id(opened.graph.head());
                            if matches!(value.as_view(), AtomView::Pow(_))
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
                if owned[&node] && matches!(value, NetworkNode::Leaf(_))
                    && (self.mode == CollectionMode::Monomials
                        || matches!(expression.as_view(), AtomView::Pow(_))
                        || (split_root_sum && node == root && matches!(expression.as_view(), AtomView::Add(_))))
                    && matches!(expression.as_view(), AtomView::Add(_) | AtomView::Mul(_) | AtomView::Pow(_))
                {
                    // This leaf was selected from its retained head inventory.
                    // Preserve the established copy-local dummy machinery when
                    // opening its arithmetic; unrelated leaves remain opaque.
                    let opened = self.open(expression.as_view())?;
                    let opened_tree = opened.graph.expr_tree().cast::<ChildVecStore<()>>();
                    let opened_root = opened.graph.graph.node_id(opened.graph.head());
                    if !matches!(opened.graph.graph[opened_root], NetworkNode::Leaf(_)) {
                        if matches!(expression.as_view(), AtomView::Pow(_))
                            && matches!(opened.graph.graph[opened_root], NetworkNode::Op(NetworkOp::Product))
                        {
                            let Some((compiled, size)) = self.compile_lowered_power(tape, &opened, &opened_tree, opened_root)? else {
                                refused = true;
                                return Ok(None);
                            };
                            return Ok(Some(TermLeaf::Reference(compiled, size)));
                        }
                        let split_root_sum = split_root_sum && node == root
                            && matches!(expression.as_view(), AtomView::Add(_));
                        let analysis = self.analyze(&opened, &opened_tree, opened_root)?;
                        return self.compile(tape, &opened, &opened_tree, opened_root, split_root_sum, &analysis)
                            .map(|compiled| compiled.map(|(node, size)| TermLeaf::Reference(node, size)));
                    }
                }
                let value =
                    SymbolicTensor::from_validated_parts(expression, interfaces[&node].clone(), Arc::clone(&histories[&node]));
                Ok(Some(self.leaf(value, owned[&node])))
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
    /// Retain the existing shallow Network and its discovered leaf interfaces
    /// with this immutable domain. A changed payload rebuilds topology while
    /// intrinsic replacements reuse the exact unchanged boundary interfaces.
    pub(crate) fn shallow_graph(&self) -> Result<Arc<super::SymbolicNet<AbstractIndex>>> {
        self.proofs
            .shallow_graph
            .get_or_init(|| {
                #[cfg(test)]
                SHALLOW_BUILDS.with(|count| count.set(count.get() + 1));
                let network = crate::shorthands::schoonschip::SlotContraction::shallow_network(
                    self.expression.as_view(),
                    (self.proofs.validated.get() == Some(&true))
                        .then(|| self.proofs.leaf_interfaces.get_or_init(Default::default)),
                )
                .map_err(|error| error.to_string())?;
                Ok(Arc::new(network.map_ref(
                    |scalar| scalar.as_view().to_owned(),
                    |tensor| SymbolicTensor {
                        expression: tensor.expression.as_view().to_owned(),
                        structure: tensor.structure.clone(),
                        is_metric: tensor.is_metric,
                        is_composite: tensor.is_composite,
                        proofs: tensor.proofs.clone(),
                    },
                )))
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
    /// selection facts can be reused without invoking it again.
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
    /// Unselected scalar and tensor subtrees remain factored.
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
        mut finalize: Option<&mut CollectionFinalizer<'_>>,
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
                finalize(&mut vec![(one, self.clone())], &[])?;
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
        let mut factors = Vec::new();
        let mut spectators = Vec::new();
        let selected_count = nodes.iter().filter(|node| owned[node]).count();
        let expose = mode == CollectionMode::Monomials && selected_count > 1;
        let lowered_power = matches!(self.expression.as_view(), AtomView::Pow(_))
            && matches!(
                network.graph.graph[root],
                NetworkNode::Op(NetworkOp::Product)
            );
        for node in nodes {
            if !owned[&node] {
                let spectator = Self::collection_product(&[Self::from_validated_parts(
                    input.emit(&network, &tree, node)?,
                    interfaces[&node].clone(),
                    Arc::clone(&histories[&node]),
                )])?;
                // These unchanged graph leaves retain their admitted boundary.
                // Use the same publication path as selected factors to decode
                // graph-local AUTO ports before reusing their regional facts.
                spectator.observe_replacement(self);
                if spectator.is_scalar()
                    && let Some(facts) = input
                        .observations
                        .certify_scalar_region(spectator.expression.as_view())
                {
                    let _ = spectator.proofs.intrinsic.set(true);
                    let _ = spectator.proofs.observations.set(Arc::new(facts));
                }
                spectators.push(spectator);
                continue;
            }
            let compiled = if (lowered_power
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
                    }) {
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
                    deferred.proofs.frontier = super::simplification::ReductionStatus::Deferred;
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
                deferred.proofs.frontier = super::simplification::ReductionStatus::Deferred;
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
                    frontier = frontier.max(mapped.proofs.frontier);
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
                    let mapped = if mapped.proofs.validated.get() == Some(&true) {
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
            finalize(&mut states, &spectators)?;
        }
        let terms = states
            .into_iter()
            .map(|(selected, coefficient)| Self::collection_product(&[selected, coefficient]))
            .collect::<Result<Vec<_>>>()?;
        let selected_interface = terms.first().map(|term| &term.structure);
        let mut facts = selected_interface
            .and_then(|interface| Self::collection_observations(&terms, interface, true));
        let selected_expression = Atom::add_many(terms.iter().map(|value| &value.expression));
        let selected_facts = facts.clone();
        let complete = remaining.is_empty();
        let mut output = vec![selected_expression.clone()];
        if !remaining.is_empty() {
            facts = None;
        }
        for node in remaining {
            let expression = network
                .graph
                .to_expression_at(&tree, node, &mut |_, value| match value {
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
                    _ => unreachable!("stored symbolic leaves checked at admission"),
                })
                .map_err(|error| TensorInferenceError::invalid(error.to_string()))?;
            // Convert the graph's occurrence-local ports before publishing a
            // remaining factor with its declared logical layout.
            let body = Self::collection_product(&[Self::from_validated_parts(
                expression,
                interfaces[&node].clone(),
                Arc::clone(&histories[&node]),
            )])?;
            output.push(body.expression);
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
        let mut root = self.with_identity_result(expression, facts)?;
        root.proofs.frontier = frontier;
        Ok(root)
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

    /// Publish a coalesced tape leaf's power using its established scope.
    fn collection_power(&self, power: u32) -> Result<Self> {
        let interface = if power.is_multiple_of(2) {
            PartialStructure::from_logical_slots([])
        } else {
            self.structure.clone()
        };
        let expression = self.expression.pow(power);
        let intrinsic = self.normalization_is_intrinsic();
        if !intrinsic {
            Self::validate_interface(&expression, &interface)?;
            Self::validate_atom(&expression)?;
        }
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
        let value = Self::from_validated_parts(expression, interface, Arc::new(history));
        if intrinsic {
            let _ = value.proofs.intrinsic.set(true);
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
            value.proofs.validated.get() == Some(&true)
                && value.expression == expression
                && value.structure == interface
        }) {
            // Multiplication by an identity did not change the exact payload or
            // its ordered interface. Keep its observations and lazy graph; an
            // encoded AUTO boundary must still pass the decoding above first.
            return Ok(value.clone());
        }
        let intrinsic = values.iter().all(Self::normalization_is_intrinsic);
        if !intrinsic {
            Self::validate_interface(&expression, &interface)?;
            Self::validate_atom(&expression)?;
        }
        let value = Self::from_validated_parts(expression, interface, Arc::new(occurrences));
        if value.expression.is_zero() {
            // A vanished selected word can leave ports joined to its coefficient.
            // Those nominal connections determine the typed zero's boundary, but
            // there is no surviving syntax to contract or inspect for identities.
            let _ = value.proofs.intrinsic.set(true);
            let _ =
                value
                    .proofs
                    .observations
                    .set(Arc::new(DomainObservations::terminal_polynomial(
                        &value.structure,
                        false,
                    )));
        } else if intrinsic {
            let _ = value.proofs.intrinsic.set(true);
            if let Some(observed) = Self::collection_observations(values, &value.structure, false) {
                let _ = value.proofs.observations.set(Arc::new(observed));
            }
        }
        Ok(value)
    }
}

#[cfg(test)]
mod conjugated_power_tests;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::epsilon::EPSILON_SYMBOL;
    use crate::tensor::composition::TensorCompositionError;
    use crate::test_support;
    use spenso::mink;
    use spenso::network::library::symbolic::ETS;
    use spenso::network::tags::SPENSO_TAG;
    use spenso::shadowing::SYM;
    use spenso::structure::representation::{Minkowski, RepName};
    use std::slice;
    use symbolica::atom::{FunctionBuilder, Symbol};
    use symbolica::symbol;

    fn tensor(head: Symbol, slots: &[Atom]) -> Atom {
        FunctionBuilder::new(head).add_args(slots).finish()
    }

    fn contract_selected(
        value: SymbolicTensor<PartialStructure>,
        _complete: bool,
        _coefficient: &mut SymbolicTensor<PartialStructure>,
    ) -> Result<SymbolicTensor<PartialStructure>> {
        let mut contracted = value.contract_parts(Default::default())?;
        contracted.root.proofs.frontier = contracted.status;
        Ok(contracted.root)
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
        let _ = survivor.proofs.observations.set(Arc::new(survivor_facts));

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
                source
                    .explicit_index_occurrences()
                    .unwrap()
                    .oriented_entries()
                    .collect::<Vec<_>>()
            );
            assert_eq!(
                result.explicit_index_occurrences().unwrap().selected_pair,
                source.explicit_index_occurrences().unwrap().selected_pair
            );
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
                source
                    .explicit_index_occurrences()
                    .unwrap()
                    .oriented_entries()
                    .collect::<Vec<_>>()
            );
            assert_eq!(
                result.explicit_index_occurrences().unwrap().selected_pair,
                source.explicit_index_occurrences().unwrap().selected_pair
            );
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
        assert_eq!(scratch.proofs.validated.get(), Some(&true));
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
        assert_eq!(body.proofs.validated.get(), Some(&true));
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
                selected.proofs.frontier = super::super::simplification::ReductionStatus::Capped;
                Ok(selected)
            },
        ).unwrap();
        assert_eq!(calls, 1);
        assert!(!complete_rows);
        assert_eq!(
            root.proofs.frontier,
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
}

#[cfg(test)]
mod planning_tests {
    use super::*;
    use spenso::{g, mink};
    use symbolica::atom::FunctionBuilder;

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
