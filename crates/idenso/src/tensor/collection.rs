//! Selected algebra on the existing parsed graph and arithmetic tape.

use std::{collections::HashMap, sync::Arc};

use linnet::{
    half_edge::{
        NodeIndex,
        subgraph::{ModifySubSet, SuBitGraph},
        tree::SimpleTraversalTree,
    },
    tree::child_vec::ChildVecStore,
};
use spenso::{
    network::{
        graph::{NetworkLeaf, NetworkNode, NetworkOp},
        library::DummyLibrary,
        parsing::{ParseSettings, ParseState, ShorthandParsing},
        store::{TensorScalarStore, TensorScalarStoreMapping},
    },
    shadowing::{TensorCollectFilter, TermLeaf, TermTape},
    structure::{
        OrderedStructure, TensorStructure,
        abstract_index::AbstractIndex,
        partial::{PartialIndex, PartialSlot, PartialStructure, PartialStructureExt},
        slot::IsAbstractSlot,
    },
};
use symbolica::{
    atom::{AliasedAtom, Atom, AtomCore, AtomOrView, AtomView},
    domains::rational::Rational,
    id::ConditionResult,
};

use super::{
    SymbolicNet, SymbolicTensor,
    aliases::{AliasInterfaces, Definition},
    inference::{InterfaceInference, LeafInference, TensorInferenceError},
};

type Result<T> = std::result::Result<T, TensorInferenceError>;
type CoefficientVisitor<'a> = dyn FnMut(
        &[(
            SymbolicTensor<PartialStructure>,
            SymbolicTensor<PartialStructure>,
        )],
        &[SymbolicTensor<PartialStructure>],
        &[Definition],
    ) -> Result<()>
    + 'a;

// Per-call intake state only. Graphs and arithmetic remain owned by Network and
// TermTape; cached templates retain borrowed source storage until this call ends.
struct CollectionInput<'a, Select> {
    definitions: &'a [Definition],
    templates: HashMap<usize, Arc<SymbolicNet<AbstractIndex, AtomOrView<'a>>>>,
    template_interfaces: HashMap<usize, PartialStructure>,
    template_owned: HashMap<usize, bool>,
    origins: HashMap<Atom, usize>,
    registry: HashMap<Atom, Definition>,
    state: ParseState<AbstractIndex>,
    definitions_intrinsic: bool,
    select: Select,
    values: Vec<(SymbolicTensor<PartialStructure>, bool)>,
    positions: HashMap<(Atom, Vec<PartialSlot>, bool), usize>,
}

type CollectionTree = SimpleTraversalTree<ChildVecStore<()>>;
type CollectionAnalysis = (
    HashMap<NodeIndex, bool>,
    HashMap<NodeIndex, PartialStructure>,
);

impl<'a, Select: FnMut(AtomView<'_>) -> bool> CollectionInput<'a, Select> {
    fn new(
        source: &SymbolicTensor<PartialStructure>,
        definitions: &'a [Definition],
        select: Select,
    ) -> Self {
        Self {
            definitions,
            templates: HashMap::new(),
            template_interfaces: HashMap::new(),
            template_owned: HashMap::new(),
            origins: definitions
                .iter()
                .enumerate()
                .map(|(i, pair)| (pair.0.expression.clone(), i))
                .collect(),
            registry: definitions
                .iter()
                .cloned()
                .map(|pair| (pair.0.expression.clone(), pair))
                .collect(),
            state: SymbolicTensor::reserved_dummies(
                std::iter::once(source).chain(definitions.iter().flat_map(|(a, b)| [a, b])),
            ),
            definitions_intrinsic: definitions.iter().all(|(_, body)| {
                InterfaceInference::normalization_is_intrinsic(body.expression.as_atom_view())
            }),
            select,
            values: Vec::new(),
            positions: HashMap::new(),
        }
    }

    fn parse(&self, source: AtomView<'a>) -> Result<SymbolicNet<AbstractIndex, AtomOrView<'a>>> {
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

    fn template(
        &mut self,
        index: usize,
    ) -> Result<Arc<SymbolicNet<AbstractIndex, AtomOrView<'a>>>> {
        if let Some(graph) = self.templates.get(&index) {
            return Ok(Arc::clone(graph));
        }
        let graph = Arc::new(self.parse(self.definitions[index].1.expression.as_atom_view())?);
        self.templates.insert(index, Arc::clone(&graph));
        Ok(graph)
    }

    fn selected(&mut self, value: AtomView<'_>, depth: usize) -> Result<bool> {
        let Some(&index) = self.origins.get(value.get_data()) else {
            return Ok((self.select)(value));
        };
        // A callback-sensitive registry stays opaque before any template
        // parsing or speculative construction. Domain mapping still visits its
        // original reachable bodies through the established checked schedule.
        if !self.definitions_intrinsic {
            return Ok(false);
        }
        if let Some(&owned) = self.template_owned.get(&index) {
            return Ok(owned);
        }
        if depth >= TermTape::<usize>::MAX_EXPANSION_DEPTH {
            return Ok(false);
        }
        let network = self.template(index)?;
        let tree = network.graph.expr_tree().cast::<ChildVecStore<()>>();
        let root = network.graph.graph.node_id(network.graph.head());
        let (owned, interfaces) = self.analyze(&network, &tree, root, depth + 1)?;
        self.template_owned.insert(index, owned[&root]);
        self.template_interfaces
            .insert(index, interfaces[&root].clone());
        Ok(owned[&root])
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
        depth: usize,
    ) -> Result<CollectionAnalysis> {
        let traversal = tree
            .iter_preorder_tree_nodes(&network.graph.graph, root)
            .collect::<Vec<_>>();
        let mut owned = HashMap::new();
        let mut interfaces = HashMap::<_, PartialStructure>::new();
        for &node in traversal.iter().rev() {
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
                        self.selected(tensor.expression.as_atom_view(), depth)?,
                        InterfaceInference::merge_explicit_interface_sequence(&[interface])?,
                    )
                }
                NetworkNode::Leaf(NetworkLeaf::Scalar(index)) => (
                    self.selected(network.store.get_scalar_ref(*index).as_atom_view(), depth)?,
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
            owned.insert(node, selected);
            interfaces.insert(node, interface);
        }
        Ok((owned, interfaces))
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
        self.definitions_intrinsic
            && tree
                .iter_preorder_tree_nodes(&network.graph.graph, node)
                .all(|node| match &network.graph.graph[node] {
                    NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => {
                        InterfaceInference::normalization_is_intrinsic(
                            network.store.tensors[*index].expression.as_atom_view(),
                        )
                    }
                    NetworkNode::Leaf(NetworkLeaf::Scalar(index)) => {
                        InterfaceInference::normalization_is_intrinsic(
                            network.store.get_scalar_ref(*index).as_atom_view(),
                        )
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
        interface: &PartialStructure,
    ) -> Result<Option<SymbolicNet<AbstractIndex, Atom>>> {
        // Prove this graph scope intrinsic before emitting a temporary body:
        // function reconstruction can itself invoke a registered normalizer.
        if !self.scope_is_intrinsic(network, tree, node) {
            return Ok(None);
        }
        let logical = PartialStructure::from_logical_slots(
            interface.logical_slots().into_iter().map(|slot| {
                let index = match slot.aind {
                    PartialIndex::Explicit(index) => {
                        InterfaceInference::partial_index(index, LeafInference::Observe)
                    }
                    index => index,
                };
                slot.rep().slot(index)
            }),
        );
        let body = SymbolicTensor::from_normalized_parts(self.emit(network, tree, node)?, logical);
        if !self.definitions_intrinsic
            || !InterfaceInference::normalization_is_intrinsic(body.expression.as_atom_view())
        {
            return Ok(None);
        }
        let handle = body.alias_handle()?;
        let replacements = interface
            .logical_slots()
            .into_iter()
            .enumerate()
            .map(|(position, slot)| {
                let PartialIndex::Explicit(index) = slot.aind else {
                    unreachable!("parsed graph boundary is encoded")
                };
                (
                    position,
                    slot.rep().slot::<AbstractIndex, _>(index).to_atom(),
                )
            })
            .collect::<HashMap<_, _>>();
        let literal =
            super::composition::rewrite_interface_ports(&handle, &replacements, &mut |_, _| {})
                .map_err(|error| TensorInferenceError::invalid(error.to_string()))?;
        let current = SymbolicTensor::from_normalized_parts(literal, interface.clone());
        let definition = (handle, body);
        let mut subset: SuBitGraph = network.graph.graph.empty_subgraph();
        for child in tree.iter_preorder_tree_nodes(&network.graph.graph, node) {
            for hedge in network.graph.graph.iter_crown(child) {
                subset.add(hedge);
            }
        }
        let mut template = network.map_ref(Clone::clone, Clone::clone);
        template.graph = template.graph.extract(&subset);
        let mut observations = Vec::new();
        let graph = SymbolicTensor::<AliasInterfaces, AliasedAtom>::expose_literal_graph(
            &template,
            &definition,
            interface,
            &current,
            &self.state,
            &mut self.registry,
            &mut |old, new| observations.push((old.clone(), new.clone())),
        )?;
        for (old, new) in observations {
            if let Some(&origin) = self.origins.get(&old) {
                self.origins.insert(new, origin);
            }
        }
        Ok(graph)
    }

    fn compile<E: AtomCore + Clone + From<E::Output>>(
        &mut self,
        tape: &mut TermTape<usize>,
        network: &SymbolicNet<AbstractIndex, E>,
        tree: &CollectionTree,
        root: NodeIndex,
        expose: bool,
        depth: usize,
    ) -> Result<Option<(usize, (usize, usize))>> {
        if depth >= TermTape::<usize>::MAX_EXPANSION_DEPTH {
            return Ok(None);
        }
        let (owned, interfaces) = self.analyze(network, tree, root, depth)?;
        // Exposure is contextual: a lone selected literal remains opaque and
        // stable. Adjacent selected factors can inspect its owned body; foreign
        // definitions never enter the selected tape.
        let mut exposure = HashMap::new();
        let mut pending = vec![(root, expose)];
        while let Some((node, inherited)) = pending.pop() {
            let children = tree
                .iter_children(node, &network.graph.graph)
                .collect::<Vec<_>>();
            let adjacent = inherited
                || match network.graph.graph[node] {
                    NetworkNode::Op(NetworkOp::Product) => {
                        children.iter().filter(|child| owned[child]).count() > 1
                    }
                    NetworkNode::Op(NetworkOp::Power(power)) => power > 1 && owned[&node],
                    _ => false,
                };
            exposure.insert(node, adjacent);
            pending.extend(children.into_iter().map(|child| (child, adjacent)));
        }
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
                    let mut copies = Vec::new();
                    for _ in 0..power {
                        let Some(graph) =
                            self.expose_scope(network, tree, base, &interfaces[&base])?
                        else {
                            refused = true;
                            return Ok(None);
                        };
                        let tree = graph.graph.expr_tree().cast::<ChildVecStore<()>>();
                        let root = graph.graph.graph.node_id(graph.graph.head());
                        let Some(copy) =
                            self.compile(tape, &graph, &tree, root, true, depth + 2)?
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
                let value =
                    SymbolicTensor::from_normalized_parts(expression, interfaces[&node].clone());
                if self.definitions_intrinsic
                    && owned[&node]
                    && exposure[&node]
                    && let Some(&index) = self.origins.get(&value.expression)
                {
                    let template = self.template(index)?;
                    let definition = &self.definitions[index];
                    let encoded = self.template_interfaces[&index].clone();
                    let mut observations = Vec::new();
                    let graph =
                        SymbolicTensor::<AliasInterfaces, AliasedAtom>::expose_literal_graph(
                            &template,
                            definition,
                            &encoded,
                            &value,
                            &self.state,
                            &mut self.registry,
                            &mut |old, new| observations.push((old.clone(), new.clone())),
                        )?;
                    for (old, new) in observations {
                        if let Some(&origin) = self.origins.get(&old) {
                            self.origins.insert(new, origin);
                        }
                    }
                    if let Some(graph) = graph {
                        let tree = graph.graph.expr_tree().cast::<ChildVecStore<()>>();
                        let root = graph.graph.graph.node_id(graph.graph.head());
                        if let Some((node, size)) =
                            self.compile(tape, &graph, &tree, root, true, depth + 1)?
                        {
                            return Ok(Some(TermLeaf::Reference(node, size)));
                        }
                        // A refused tape must be discarded by the caller.
                        refused = true;
                        return Ok(None);
                    }
                }
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
            Some(&mut |states, spectators, registry| {
                let spectator = Self::collection_product(spectators)?;
                proof = ConditionResult::True;
                for (selected, coefficient) in states {
                    if selected.expression.is_zero() {
                        continue;
                    }
                    let coefficient =
                        Self::collection_product(&[coefficient.clone(), spectator.clone()])?
                            .with_aliases(registry.iter().cloned())?
                            .resolved()?;
                    // Symbolica's zero-iteration check keeps unresolved sums
                    // inconclusive. Only exact factor cancellation can prove zero.
                    match coefficient.expression.collect_factors().zero_test(0, 0.0) {
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
            &[],
            |value| filter.matches(value),
            |value, _, _| Ok((value, vec![])),
        )?;
        Ok(proof)
    }

    /// Return complete selected sectors and their factored coefficients.
    ///
    /// This is an explicit output of the same graph/tape collection used by
    /// `collect`. Each pair resolves only its own alias dependencies; scalar
    /// spectators are retained in the coefficient without polynomial expansion.
    /// A refused or capped traversal returns an error, never partial rows.
    pub fn coefficient_list<const N: usize>(
        &self,
        filter: TensorCollectFilter<N>,
    ) -> Result<Vec<(Self, Self)>> {
        let mut rows = None;
        self.collect_with_map(
            Some(&mut |states, spectators, registry| {
                let spectator = Self::collection_product(spectators)?;
                rows = Some(
                    states
                        .iter()
                        .map(|(selected, coefficient)| {
                            Ok((
                                selected
                                    .clone()
                                    .with_aliases(registry.iter().cloned())?
                                    .resolved()?,
                                Self::collection_product(&[
                                    coefficient.clone(),
                                    spectator.clone(),
                                ])?
                                .with_aliases(registry.iter().cloned())?
                                .resolved()?,
                            ))
                        })
                        .collect::<Result<Vec<_>>>()?,
                );
                Ok(())
            }),
            &[],
            |value| filter.matches(value),
            |value, _, _| Ok((value, vec![])),
        )?;
        rows.ok_or_else(|| {
            TensorInferenceError::invalid("selected coefficient collection did not complete")
        })
    }

    /// Keep selected tensor coefficients in the existing literal alias registry.
    /// Unselected scalar and tensor subtrees remain factored.
    pub fn collect<const N: usize>(
        &self,
        filter: TensorCollectFilter<N>,
    ) -> Result<SymbolicTensor<AliasInterfaces, AliasedAtom>> {
        let mut handles = HashMap::<_, Self>::new();
        let (root, definitions) = self.collect_with_map(
            None,
            &[],
            |value| filter.matches(value),
            |value, _registry, _complete| {
                if value.expression.is_one() {
                    return Ok((value, vec![]));
                }
                let key = (value.expression.clone(), value.structure.logical_slots());
                if let Some(handle) = handles.get(&key) {
                    return Ok((handle.clone(), vec![]));
                }
                let handle = value.alias_handle()?;
                handles.insert(key, handle.clone());
                Ok((handle.clone(), vec![(handle, value)]))
            },
        )?;
        root.with_aliases(definitions)
    }

    /// Visit selected monomials through the graph's established factor order.
    /// Domain callbacks receive no outer product of sums. Generated definitions
    /// are returned to the alias owner, not resolved into the input expression.
    /// Each callback receives every literal definition retained so far, including
    /// earlier frontier results. It may use these handles without traversing the
    /// new bodies; original definition traversal belongs to the outer scheduler.
    /// The final argument is true only when every selected factor has entered
    /// this frontier. It does not certify normalization callbacks or mark a
    /// capped/refused prefix as complete.
    pub(crate) fn collect_with_map(
        &self,
        mut coefficients: Option<&mut CoefficientVisitor<'_>>,
        definitions: &[Definition],
        mut select: impl FnMut(AtomView<'_>) -> bool,
        mut step: impl FnMut(Self, &[Definition], bool) -> Result<(Self, Vec<Definition>)>,
    ) -> Result<(Self, Vec<Definition>)> {
        let mut input = CollectionInput::new(self, definitions, &mut select);
        let network = input.parse(self.expression.as_atom_view())?;
        let tree = network.graph.expr_tree().cast::<ChildVecStore<()>>();
        let root = network.graph.graph.node_id(network.graph.head());
        let (owned, interfaces) = input.analyze(&network, &tree, root, 0)?;
        if !owned[&root] {
            if let Some(visit) = coefficients.as_mut() {
                let one = Self::from_normalized_parts(
                    Atom::one(),
                    PartialStructure::from_logical_slots([]),
                );
                visit(&[(one, self.clone())], &[], definitions)?;
            }
            return Ok((self.clone(), vec![]));
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
        let expose = nodes.iter().filter(|node| owned[node]).count() > 1;
        for node in nodes {
            if !owned[&node] {
                spectators.push(Self::from_normalized_parts(
                    input.emit(&network, &tree, node)?,
                    interfaces[&node].clone(),
                ));
                continue;
            }
            let compiled = if expose
                && matches!(
                    network.graph.graph[node],
                    NetworkNode::Op(NetworkOp::Sum | NetworkOp::Product)
                )
                && tree
                    .iter_preorder_tree_nodes(&network.graph.graph, node)
                    .any(|member| {
                        network
                            .graph
                            .graph
                            .iter_crown(member)
                            .any(|hedge| network.graph.graph[[&hedge]].is_slot())
                    }) {
                // A single selected subtree retains its existing dummy scope;
                // only separately expanded factors can collide. Explicit powers
                // scope their copies in compile_power. A scalar subtree has no
                // graph dummy identities to freshen.
                // Keep its selected leaves available even when their local
                // domain step uses a checked normalization callback.
                // The parser can lower a tensor square into two product arms.
                // Scope each already parsed factor's internal incidences before
                // tape coalescence; free boundary names still join across arms.
                let Some(graph) = input.expose_scope(&network, &tree, node, &interfaces[&node])?
                else {
                    return Ok((self.clone(), vec![]));
                };
                let scoped_tree = graph.graph.expr_tree().cast::<ChildVecStore<()>>();
                let scoped_root = graph.graph.graph.node_id(graph.graph.head());
                input.compile(&mut tape, &graph, &scoped_tree, scoped_root, expose, 0)?
            } else {
                input.compile(&mut tape, &network, &tree, node, expose, 0)?
            };
            let Some((compiled, _)) = compiled else {
                // No domain callback has run; discard all refused tape scratch.
                return Ok((self.clone(), vec![]));
            };
            factors.push((node, compiled));
        }
        let values = &input.values;
        let scalar = PartialStructure::from_logical_slots([]);
        let one = Self::from_normalized_parts(Atom::one(), scalar);
        let mut states = vec![(one.clone(), one)];
        let mut additional = input
            .registry
            .values()
            .filter(|pair| {
                !definitions
                    .iter()
                    .any(|old| old.0.expression == pair.0.expression)
            })
            .cloned()
            .collect::<Vec<_>>();
        additional.sort_unstable_by(|a, b| {
            AtomView::cmp(&a.0.expression.as_view(), &b.0.expression.as_view())
        });
        let mut available = definitions.to_vec();
        available.extend(additional.iter().cloned());
        let mut register = |pair: Definition,
                            available: &mut Vec<Definition>,
                            additional: &mut Vec<Definition>| {
            if let Some(previous) = input.registry.get(&pair.0.expression) {
                if previous != &pair {
                    return Err(TensorInferenceError::invalid(
                        "conflicting selected-domain alias definitions",
                    ));
                }
            } else {
                input
                    .registry
                    .insert(pair.0.expression.clone(), pair.clone());
                available.push(pair.clone());
                additional.push(pair);
            }
            Ok(())
        };
        let mut generated = 0usize;
        let mut remaining = Vec::new();
        for (level, &(_, factor)) in factors.iter().enumerate() {
            tape.clear_terms();
            if tape
                .distribute(&mut vec![(factor, 1)], &Rational::one(), 0)
                .is_none()
            {
                remaining.extend(factors[level..].iter().map(|(node, _)| *node));
                break;
            }
            let rows = tape.take_terms();
            let count = states.len().saturating_mul(rows.len());
            let state_bytes = states
                .iter()
                .map(|(selected, coefficient)| {
                    selected
                        .expression
                        .as_atom_view()
                        .get_byte_size()
                        .saturating_add(coefficient.expression.as_atom_view().get_byte_size())
                        .saturating_add(2 * std::mem::size_of::<Self>())
                })
                .max()
                .unwrap_or(0);
            let row_bytes = rows
                .iter()
                .map(|(powers, _)| {
                    powers
                        .iter()
                        .map(|(position, _)| {
                            values[tape.leaves()[*position]]
                                .0
                                .expression
                                .as_atom_view()
                                .get_byte_size()
                        })
                        .fold(0usize, usize::saturating_add)
                })
                .max()
                .unwrap_or(0);
            let definition_bytes = additional
                .iter()
                .map(|(handle, body): &Definition| {
                    handle
                        .expression
                        .as_atom_view()
                        .get_byte_size()
                        .saturating_add(body.expression.as_atom_view().get_byte_size())
                })
                .fold(0usize, usize::saturating_add);
            // This predicts the next frontier and retained definitions. Atom
            // allocator overhead is not a hard heap bound. Refusal retains the
            // exact current frontier and unprocessed graph factors as aliases.
            if generated.saturating_add(count) > TermTape::<usize>::MAX_GENERATED_TERMS
                || count
                    .saturating_mul(state_bytes.saturating_add(row_bytes).saturating_mul(4))
                    .saturating_add(definition_bytes)
                    > TermTape::<usize>::MAX_GENERATED_FACTOR_BYTES
            {
                remaining.extend(factors[level..].iter().map(|(node, _)| *node));
                break;
            }
            generated += count;
            let mut next = Vec::<(Self, Vec<Self>)>::new();
            let mut keys = HashMap::new();
            for (selected, coefficient) in states {
                for (powers, number) in &rows {
                    let mut selected_factors = vec![selected.clone()];
                    let mut coefficient_factors = vec![coefficient.clone()];
                    for &(position, power) in powers {
                        let (value, owned) = &values[tape.leaves()[position]];
                        let interface = if power % 2 == 0 {
                            PartialStructure::from_logical_slots([])
                        } else {
                            value.structure.clone()
                        };
                        let expression = value.expression.pow(power);
                        if !InterfaceInference::normalization_is_intrinsic(
                            value.expression.as_atom_view(),
                        ) {
                            Self::validate_interface(&expression, &interface)?;
                            Self::validate_atom(&expression)?;
                        }
                        let value = Self::from_normalized_parts(expression, interface);
                        if *owned {
                            selected_factors.push(value);
                        } else {
                            coefficient_factors.push(value);
                        }
                    }
                    let monomial = Self::collection_product(&selected_factors)?;
                    let (mapped, emitted) =
                        step(monomial.clone(), &available, level + 1 == factors.len())?;
                    // The callback owns a domain identity, but its returned carrier
                    // may carry stale metadata. Check its payload against the trusted
                    // input, including encoded AUTO order, before retaining any result.
                    let mut mapped = if monomial.structure.open_positions().is_empty() {
                        monomial.with_checked_expression(mapped.expression)?
                    } else {
                        monomial.with_rewritten_expression(mapped.expression)?
                    };
                    for pair in emitted {
                        register(pair, &mut available, &mut additional)?;
                    }
                    // A domain can return a sum. Keep it in the same registry
                    // so the next local step still receives a monomial; newly
                    // generated definitions are visited on the next domain pass.
                    if matches!(mapped.expression.as_atom_view(), AtomView::Add(_)) {
                        let handle = mapped.alias_handle()?;
                        register((handle.clone(), mapped), &mut available, &mut additional)?;
                        mapped = handle;
                    }
                    let coefficient = Self::collection_product(&coefficient_factors)?;
                    let coefficient = coefficient
                        .with_algebra_result(Atom::num(number.clone()) * &coefficient.expression)?;
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
                    let body = first.with_algebra_result(Atom::add_many(
                        coefficients.iter().map(|value| &value.expression),
                    ))?;
                    if coefficients.len() == 1
                        || !matches!(body.expression.as_atom_view(), AtomView::Add(_))
                    {
                        return Ok((selected, body));
                    }
                    let handle = body.alias_handle()?;
                    register((handle.clone(), body), &mut available, &mut additional)?;
                    Ok((selected, handle))
                })
                .collect::<Result<Vec<_>>>()?;
        }
        // These are the complete coefficient rows. Definitions remain a
        // dependency registry; a capped frontier cannot certify its coefficients.
        if remaining.is_empty()
            && let Some(visit) = coefficients.as_mut()
        {
            visit(&states, &spectators, &available)?;
        }
        let terms = states
            .into_iter()
            .map(|(selected, coefficient)| {
                Self::collection_product(&[selected, coefficient]).map(|value| value.expression)
            })
            .collect::<Result<Vec<_>>>()?;
        let mut output = vec![Atom::add_many(terms)];
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
            let body = Self::from_normalized_parts(expression, interfaces[&node].clone());
            let handle = body.alias_handle()?;
            output.push(handle.expression.clone());
            register((handle, body), &mut available, &mut additional)?;
        }
        output.extend(spectators.into_iter().map(|value| value.expression));
        let expression = Atom::mul_many(output);
        let root = self.with_rewritten_expression(expression)?;
        Ok((root, additional))
    }

    fn collection_product(values: &[Self]) -> Result<Self> {
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
        let intrinsic = values.iter().all(|value| {
            InterfaceInference::normalization_is_intrinsic(value.expression.as_atom_view())
        });
        let expression = Atom::mul_many(values.iter().map(|value| &value.expression));
        if !intrinsic {
            Self::validate_interface(&expression, &interface)?;
            Self::validate_atom(&expression)?;
        }
        Ok(Self::from_normalized_parts(expression, interface))
    }
}

impl SymbolicTensor<AliasInterfaces, AliasedAtom> {
    pub fn collect<const N: usize>(
        self: &Arc<Self>,
        filter: TensorCollectFilter<N>,
    ) -> Result<Arc<Self>> {
        let definitions = self.aliases()?;
        let registered = definitions
            .iter()
            .map(|pair| pair.0.expression.clone())
            .collect::<std::collections::HashSet<_>>();
        let mut handles = definitions
            .iter()
            .map(|(handle, body)| {
                (
                    (body.expression.clone(), body.structure.logical_slots()),
                    handle.clone(),
                )
            })
            .collect::<HashMap<_, _>>();
        let (root, additional) = self.root().collect_with_map(
            None,
            &definitions,
            |value| filter.matches(value),
            |value, _registry, _complete| {
                if value.expression.is_one() || registered.contains(&value.expression) {
                    return Ok((value, vec![]));
                }
                let key = (value.expression.clone(), value.structure.logical_slots());
                if let Some(handle) = handles.get(&key) {
                    return Ok((handle.clone(), vec![]));
                }
                let handle = value.alias_handle()?;
                handles.insert(key, handle.clone());
                Ok((handle.clone(), vec![(handle, value)]))
            },
        )?;
        if root == self.root() && additional.is_empty() {
            return Ok(Arc::clone(self));
        }
        Ok(Arc::new(
            root.with_aliases(definitions.into_iter().chain(additional))?,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use spenso::network::{library::symbolic::ETS, tags::SPENSO_TAG};
    use symbolica::atom::{FunctionBuilder, Symbol};

    fn tensor(head: Symbol, slots: &[Atom]) -> Atom {
        FunctionBuilder::new(head).add_args(slots).finish()
    }

    fn contract_selected(
        value: SymbolicTensor<PartialStructure>,
        registry: &[Definition],
        _complete: bool,
    ) -> Result<(SymbolicTensor<PartialStructure>, Vec<Definition>)> {
        let contracted = value.contract_parts(Default::default())?;
        let mut registry = registry
            .iter()
            .cloned()
            .map(|pair| (pair.0.expression.clone(), pair))
            .collect();
        SymbolicTensor::<AliasInterfaces, AliasedAtom>::retain_contracted_definitions(
            &contracted,
            &mut registry,
        )?;
        Ok((contracted.root, registry.into_values().collect()))
    }

    #[test]
    fn collection_retains_scalar_spectators_and_reuses_literal_result() {
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
        let result = Arc::new(
            source
                .collect(TensorCollectFilter::<0>::TaggedTensors)
                .unwrap(),
        );
        assert_eq!(
            result.resolved().unwrap().expression.expand(),
            source.expression.expand()
        );
        assert!(!result.expression.get_aliases().is_empty());
        assert!(
            matches!(result.expression.get_root().as_atom_view(), AtomView::Mul(product)
            if product.iter().any(|factor| factor == spectator.as_atom_view()))
        );
        let repeated = result
            .collect(TensorCollectFilter::<0>::TaggedTensors)
            .unwrap();
        assert!(Arc::ptr_eq(&result, &repeated));
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
        let (root, definitions) = source
            .collect_with_map(
                None,
                &[],
                |value| matches!(value, AtomView::Fun(f) if f.get_symbol() == selected),
                |value, _registry, _complete| {
                    assert!(!matches!(value.expression.as_atom_view(), AtomView::Add(_)));
                    seen.push(value.expression.clone());
                    Ok((value, vec![]))
                },
            )
            .unwrap();
        assert!(!seen.is_empty());
        assert!(
            seen.iter()
                .all(|value| value == &tensor(selected, std::slice::from_ref(&a)))
        );
        assert_eq!(
            root.with_aliases(definitions)
                .unwrap()
                .resolved()
                .unwrap()
                .expression
                .expand(),
            source.expression.expand()
        );
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
                    None,
                    &[],
                    |_| true,
                    |value, _registry, _complete| {
                        Ok((
                            SymbolicTensor::from_normalized_parts(Atom::one(), value.structure),
                            vec![],
                        ))
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
                    None,
                    &[],
                    |value| TensorCollectFilter::<0>::TaggedTensors.matches(value),
                    contract_selected
                )
                .is_err()
        );
    }

    #[test]
    fn selected_dual_alias_keeps_nested_callback_templates_opaque() {
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
        let handle = body.alias_handle().unwrap();
        calls.store(0, Ordering::Relaxed);
        let mut steps = 0;
        let (root, additional) = handle
            .collect_with_map(
                None,
                &[(handle.clone(), body)],
                |_| true,
                |value, _, _complete| {
                    steps += 1;
                    Ok((value, vec![]))
                },
            )
            .unwrap();
        assert_eq!(root, handle);
        assert!(additional.is_empty());
        assert_eq!(steps, 0);
        assert_eq!(calls.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn selected_alias_body_is_visible_beside_another_selected_factor() {
        crate::test_support::test_initialize();
        let a = spenso::mink!(4, 97007);
        let b = spenso::mink!(4, 97008);
        let head = SPENSO_TAG.tensor_symbol("selected_collection_adjacent");
        let body = SymbolicTensor::infer(tensor(head, &[a])).unwrap();
        let handle = body.alias_handle().unwrap();
        let right = tensor(head, &[b]);
        let source = SymbolicTensor::infer(&handle.expression * &right).unwrap();
        let mut seen = Vec::new();
        let (root, definitions) = source
            .collect_with_map(
                None,
                &[(handle, body.clone())],
                |value| matches!(value, AtomView::Fun(f) if f.get_symbol() == head),
                |value, _registry, complete| {
                    seen.push((value.expression.clone(), complete));
                    Ok((value, vec![]))
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
        assert_eq!(
            root.with_aliases(definitions)
                .unwrap()
                .resolved()
                .unwrap()
                .expression,
            &body.expression * right
        );
    }
    #[test]
    fn selected_power_copies_preserve_independent_internal_dummy_scopes() {
        crate::test_support::test_initialize();
        let p = SPENSO_TAG.rank_one_tensor_symbol("selected_power_p");
        let q = SPENSO_TAG.rank_one_tensor_symbol("selected_power_q");
        let a = spenso::mink!(4, 97011);
        let b = spenso::mink!(4, 97012);
        let x = Atom::var(symbolica::symbol!("selected_power_x"));
        let dot = tensor(
            ETS.metric,
            &[
                tensor(p, &[spenso::mink!(4)]),
                tensor(q, &[spenso::mink!(4)]),
            ],
        );
        let base = tensor(p, std::slice::from_ref(&a)) * tensor(q, std::slice::from_ref(&a)) + &x;
        for power in [2, 3] {
            let source = SymbolicTensor::infer(base.pow(power)).unwrap();
            let (root, definitions) = source
                .collect_with_map(
                    None,
                    &[],
                    |value| TensorCollectFilter::<0>::TaggedTensors.matches(value),
                    contract_selected,
                )
                .unwrap();
            let actual = root.with_aliases(definitions).unwrap().resolved().unwrap();
            assert_eq!(actual.expression.expand(), (&dot + &x).pow(power).expand());
        }
        let open = (tensor(p, std::slice::from_ref(&a))
            * tensor(q, &[a])
            * tensor(p, std::slice::from_ref(&b))
            + &x * tensor(q, &[b]))
        .pow(2);
        let source = SymbolicTensor::infer(open).unwrap();
        let (root, definitions) = source
            .collect_with_map(
                None,
                &[],
                |value| TensorCollectFilter::<0>::TaggedTensors.matches(value),
                contract_selected,
            )
            .unwrap();
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
        assert_eq!(
            root.with_aliases(definitions)
                .unwrap()
                .resolved()
                .unwrap()
                .expression
                .expand(),
            expected.expand()
        );
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
        let (root, definitions) = source
            .collect_with_map(
                None,
                &[],
                |value| matches!(value, AtomView::Fun(fun) if fun.get_symbol() == hook),
                |value, _, _complete| {
                    seen.push(value.expression.clone());
                    if value.expression == selected {
                        Ok((value.with_checked_expression(Atom::num(2))?, vec![]))
                    } else {
                        assert!(value.expression.is_one());
                        Ok((value, vec![]))
                    }
                },
            )
            .unwrap();
        assert_eq!(seen.iter().filter(|value| **value == selected).count(), 1);
        assert_eq!(
            root.with_aliases(definitions)
                .unwrap()
                .resolved()
                .unwrap()
                .expression,
            Atom::num(3) * x
        );
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
        let (root, definitions) = source
            .collect_with_map(
                None,
                &[],
                |value| TensorCollectFilter::<0>::TaggedTensors.matches(value),
                |value, _, _complete| {
                    if value.expression == selected {
                        selected_calls += 1;
                        Ok((value.with_checked_expression(target.clone())?, vec![]))
                    } else {
                        assert_eq!(value.expression, target);
                        Ok((value, vec![]))
                    }
                },
            )
            .unwrap();
        assert_eq!(selected_calls, 1);
        assert_eq!(
            root.with_aliases(definitions)
                .unwrap()
                .resolved()
                .unwrap()
                .expression,
            Atom::num(2) * x * target
        );
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
    fn coefficient_zero_proof_preserves_large_foreign_factors_and_refuses_capped_selection() {
        crate::test_support::test_initialize();
        let t = tensor(spenso::tensor_symbol!("capped_proof_t"), &[]);
        let u = tensor(spenso::tensor_symbol!("capped_proof_u"), &[]);
        let spectator = Atom::mul_many((0..30).map(|index| {
            Atom::var(symbolica::symbol!(&format!("zero_spectator_{index}"))) + Atom::one()
        }));
        let source = SymbolicTensor::infer(&spectator * &t).unwrap();
        let mut observed = false;
        source
            .collect_with_map(
                Some(&mut |_, spectators, _| {
                    assert_eq!(spectators.len(), 30);
                    assert_eq!(
                        Atom::mul_many(spectators.iter().map(|value| &value.expression)),
                        spectator,
                    );
                    observed = true;
                    Ok(())
                }),
                &[],
                |value| TensorCollectFilter::<0>::TaggedTensors.matches(value),
                |value, _, _| Ok((value, vec![])),
            )
            .unwrap();
        assert!(observed);
        assert_eq!(
            source
                .coefficients_are_zero(TensorCollectFilter::<0>::TaggedTensors)
                .unwrap(),
            ConditionResult::Inconclusive,
        );
        let capped = SymbolicTensor::infer((t + u).pow(32)).unwrap();
        assert_eq!(
            capped
                .coefficients_are_zero(TensorCollectFilter::<0>::TaggedTensors)
                .unwrap(),
            ConditionResult::Inconclusive,
        );
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
        let (result, definitions) = source
            .collect_with_map(
                None,
                &[],
                |_| true,
                |selected, _, _complete| {
                    calls += 1;
                    assert_eq!(
                        selected.structure.logical_slots(),
                        source.structure.logical_slots()
                    );
                    assert_eq!(selected.structure.open_positions(), vec![0, 1]);
                    Ok((selected, vec![]))
                },
            )
            .unwrap();
        assert_eq!(calls, 1);
        assert!(definitions.is_empty());
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
    fn coefficient_list_preserves_empty_selection_and_refuses_incomplete_rows() {
        crate::test_support::test_initialize();
        let x = Atom::var(symbolica::symbol!("coefficient_output_scalar"));
        let source = SymbolicTensor::checked_parts(
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
        let t = tensor(spenso::tensor_symbol!("coefficient_output_capped_t"), &[]);
        let u = tensor(spenso::tensor_symbol!("coefficient_output_capped_u"), &[]);
        let capped = SymbolicTensor::infer((t + u).pow(32)).unwrap();
        let error = capped
            .coefficient_list(TensorCollectFilter::<0>::TaggedTensors)
            .unwrap_err();
        assert!(error.to_string().contains("did not complete"));
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
