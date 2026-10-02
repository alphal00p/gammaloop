//! Operation-local assembly for the existing parser. Recognition and payload
//! construction still happen in source order. Only graph/store emission is
//! deferred, so discarded scalar temporaries cannot change the final layout.

use crate::{
    network::{
        Network, NetworkState, TensorNetworkError,
        graph::{
            NetworkEdge, NetworkGraph, NetworkGraphBuilder, NetworkLeaf, NetworkNode, NetworkOp,
        },
        store::TensorScalarStore,
    },
    structure::{
        CanonicalLayout, Canonicalized, StructureError, TensorStructure,
        representation::{LibrarySlot, RepName},
        slot::{AbsInd, DualSlotTo, IsAbstractSlot},
    },
};
use linnet::half_edge::{
    NodeIndex,
    involution::{Flow, Hedge},
};
use std::{collections::BTreeMap, fmt::Debug, ops::Range};
use symbolica::atom::Symbol;

/// Bookkeeping for a subtree already recognized by the parser. This is not a
/// second expression graph: all topology is in the one existing Linnet builder.
struct Subtree<Aind> {
    state: NetworkState,
    children: Vec<NodeIndex>,
    // The existing product constructor assigns child inputs from last to first.
    // Other operations retain their child order; this records that topology.
    reversed_head_inputs: bool,
    hedges: Range<usize>,
    slots: Vec<(LibrarySlot<Aind>, Hedge)>,
    tensor_slots: Option<Vec<LibrarySlot<Aind>>>,
    layout: Option<CanonicalLayout>,
}

pub(super) struct Construction<S, K, Aind> {
    builder: Option<NetworkGraphBuilder<K, Symbol, Aind>>,
    store: S,
    subtrees: Vec<Subtree<Aind>>,
    slot_connections: Vec<(Hedge, Hedge)>,
    n_hedges: usize,
}

impl<S: TensorScalarStore, K: Debug, Aind: AbsInd> Construction<S, K, Aind> {
    pub(super) fn new() -> Self {
        Self {
            builder: Some(NetworkGraphBuilder::new()),
            ..Self::structure_only()
        }
    }

    /// Observe the existing parser's interface without emitting graph topology.
    pub(super) fn structure_only() -> Self {
        Self {
            builder: None,
            store: S::default(),
            subtrees: Vec::new(),
            slot_connections: Vec::new(),
            n_hedges: 0,
        }
    }

    pub(super) fn emits_graph(&self) -> bool {
        self.builder.is_some()
    }

    pub(super) fn state(&self, node: NodeIndex) -> NetworkState {
        self.subtrees[node.0].state
    }
    pub(super) fn slots(&self, node: NodeIndex) -> Vec<LibrarySlot<Aind>> {
        self.subtrees[node.0]
            .slots
            .iter()
            .map(|(slot, _)| *slot)
            .collect()
    }
    fn edge(&mut self, node: NodeIndex, edge: NetworkEdge<Aind>, flow: Flow) -> Hedge {
        let orientation = match edge {
            NetworkEdge::Slot(slot) => slot.rep_name().orientation(),
            _ => true.into(),
        };
        let hedge = Hedge(self.n_hedges);
        self.n_hedges += 1;
        if let Some(builder) = &mut self.builder {
            builder.add_external_edge(node, edge, orientation, flow);
        }
        hedge
    }
    fn node(
        &mut self,
        node: NetworkNode<K, Symbol, Aind>,
        state: NetworkState,
        children: Vec<NodeIndex>,
    ) -> NodeIndex {
        let id = NodeIndex(self.subtrees.len());
        if let Some(builder) = &mut self.builder {
            assert_eq!(builder.add_node(node), id);
        }
        self.subtrees.push(Subtree {
            state,
            children,
            reversed_head_inputs: false,
            hedges: self.n_hedges..self.n_hedges,
            slots: Vec::new(),
            tensor_slots: None,
            layout: None,
        });
        self.edge(id, NetworkEdge::Head, Flow::Source);
        id
    }
    fn close(&mut self, node: NodeIndex) -> NodeIndex {
        self.subtrees[node.0].hedges.end = self.n_hedges;
        node
    }
    fn slot_state(slots: &[(LibrarySlot<Aind>, Hedge)]) -> NetworkState {
        if slots.is_empty() {
            NetworkState::Scalar
        } else if slots.iter().all(|(s, _)| s.rep.rep.is_self_dual()) {
            NetworkState::SelfDualTensor
        } else {
            NetworkState::Tensor
        }
    }
    pub(super) fn scalar(&mut self, scalar: S::Scalar) -> NodeIndex {
        let index = self.store.add_scalar(scalar);
        let node = self.node(
            NetworkNode::Leaf(NetworkLeaf::Scalar(index.into())),
            NetworkState::PureScalar,
            Vec::new(),
        );
        self.close(node)
    }
    pub(super) fn tensor(
        &mut self,
        tensor: S::Tensor,
        layout: CanonicalLayout,
    ) -> Result<NodeIndex, StructureError>
    where
        S::Tensor: TensorStructure,
        <S::Tensor as TensorStructure>::Slot: IsAbstractSlot<Aind = Aind>,
    {
        let slots = tensor
            .external_structure_iter()
            .map(|s| s.to_lib())
            .collect::<Vec<_>>();
        if layout.order() != slots.len() {
            return Err(StructureError::WrongNumberOfArguments(
                layout.order(),
                slots.len(),
            ));
        }
        let index = self.store.add_tensor(tensor);
        Ok(self.tensor_node(NetworkLeaf::LocalTensor(index), slots, Some(layout)))
    }
    pub(super) fn library_tensor<T>(&mut self, tensor: &T, key: Canonicalized<K>) -> NodeIndex
    where
        T: TensorStructure,
        T::Slot: IsAbstractSlot<Aind = Aind>,
    {
        let slots = tensor
            .external_structure_iter()
            .map(|s| s.to_lib())
            .collect::<Vec<_>>();
        let layout = (key.layout().order() == slots.len()).then(|| key.layout().clone());
        self.tensor_node(NetworkLeaf::library_key(key), slots, layout)
    }
    fn tensor_node(
        &mut self,
        leaf: NetworkLeaf<K, Aind>,
        slots: Vec<LibrarySlot<Aind>>,
        layout: Option<CanonicalLayout>,
    ) -> NodeIndex {
        assert!(
            slots.len() <= usize::from(u8::MAX) + 1,
            "tensor slots must match their new graph crown"
        );
        let node = self.node(NetworkNode::Leaf(leaf), NetworkState::Scalar, Vec::new());
        let mut exposed = slots
            .iter()
            .map(|slot| {
                let hedge = self.edge(node, NetworkEdge::Slot(*slot), Flow::Source);
                Some((*slot, hedge))
            })
            .collect::<Vec<_>>();
        // The existing sewing owner chooses the last matching old endpoint and
        // the first greater endpoint, then repeats with those endpoints removed.
        loop {
            let mut pair = None;
            for i in 0..exposed.len() {
                if let Some((left, _)) = exposed[i]
                    && let Some(j) = (i + 1..exposed.len())
                        .find(|&j| exposed[j].is_some_and(|(right, _)| left.matches(&right)))
                {
                    pair = Some((i, j));
                }
            }
            let Some((i, j)) = pair else { break };
            self.slot_connections
                .push((exposed[i].take().unwrap().1, exposed[j].take().unwrap().1));
        }
        let exposed = exposed.into_iter().flatten().collect::<Vec<_>>();
        self.subtrees[node.0].state = Self::slot_state(&exposed);
        self.subtrees[node.0].slots = exposed;
        self.subtrees[node.0].tensor_slots = Some(slots);
        self.subtrees[node.0].layout = layout;
        self.close(node)
    }
    pub(super) fn product(&mut self, children: Vec<NodeIndex>) -> NodeIndex {
        assert!(!children.is_empty());
        let mut state = self.state(children[0]);
        for &child in &children[1..] {
            state *= self.state(child);
        }
        let node = self.node(NetworkNode::Op(NetworkOp::Product), state, children.clone());
        self.subtrees[node.0].reversed_head_inputs = true;
        for _ in &children {
            self.edge(node, NetworkEdge::Head, Flow::Sink);
        }
        let mut pending: Vec<(LibrarySlot<Aind>, Hedge)> = Vec::new();
        for child in children {
            let mut slots = self.subtrees[child.0]
                .slots
                .iter()
                .copied()
                .map(Some)
                .collect::<Vec<_>>();
            for index in (0..pending.len()).rev() {
                let (left, source) = pending[index];
                if let Some(right) = slots
                    .iter_mut()
                    .find(|right| right.as_ref().is_some_and(|(right, _)| left.matches(right)))
                {
                    self.slot_connections
                        .push((source, right.take().unwrap().1));
                    pending.remove(index);
                }
            }
            pending.extend(slots.into_iter().flatten());
        }
        if state.is_tensor() {
            self.subtrees[node.0].state = Self::slot_state(&pending);
        }
        self.subtrees[node.0].slots = pending;
        self.close(node)
    }
    pub(super) fn sum(
        &mut self,
        children: Vec<NodeIndex>,
    ) -> Result<NodeIndex, TensorNetworkError<K, Symbol>>
    where
        K: std::fmt::Display,
    {
        assert!(!children.is_empty());
        let slots = self.slots(children[0]);
        let mut state = self.state(children[0]);
        for &child in &children[1..] {
            if !state.is_compatible(&self.state(child)) {
                return Err(TensorNetworkError::IncompatibleSummand(format!(
                    "sum operand states {state:?} and {:?} are incompatible",
                    self.state(child),
                )));
            }
            state += self.state(child);
        }
        let node = self.node(NetworkNode::Op(NetworkOp::Sum), state, children.clone());
        for _ in &children {
            self.edge(node, NetworkEdge::Head, Flow::Sink);
        }
        let slot_count = slots.len();
        let mut slot_inputs: BTreeMap<LibrarySlot<Aind>, (Vec<Hedge>, usize)> = BTreeMap::new();
        for slot in slots {
            let hedge = self.edge(node, NetworkEdge::Slot(slot), Flow::Source);
            self.subtrees[node.0].slots.push((slot, hedge));
            let (inputs, multiplicity) = slot_inputs.entry(slot).or_default();
            *multiplicity += 1;
            for _ in &children {
                let input = self.edge(node, NetworkEdge::Slot(slot), Flow::Sink);
                inputs.push(input);
            }
        }
        let mut remaining_children = children.len();
        for child in children {
            remaining_children -= 1;
            if self.subtrees[child.0].slots.len() != slot_count {
                return Err(TensorNetworkError::IncompatibleSummand(
                    "sum operands have different numbers of exposed ports".to_owned(),
                ));
            }
            for (slot, hedge) in &self.subtrees[child.0].slots {
                self.slot_connections.push((
                    slot_inputs
                        .get_mut(slot)
                        // Reserve each later child's share. Equal total rank alone
                        // does not establish equal port identities or multiplicities.
                        .filter(|(inputs, multiplicity)| inputs.len() > remaining_children * *multiplicity)
                        .and_then(|(inputs, _)| inputs.pop())
                        .ok_or_else(|| TensorNetworkError::IncompatibleSummand(
                            "sum operands have different exposed port identities or multiplicities".to_owned(),
                        ))?,
                    *hedge,
                ));
            }
        }
        Ok(self.close(node))
    }
    pub(super) fn power(&mut self, child: NodeIndex, power: i8) -> NodeIndex {
        let state = self.state(child).pow(power);
        let node = self.node(NetworkNode::Op(NetworkOp::Power(power)), state, vec![child]);
        self.edge(node, NetworkEdge::Head, Flow::Sink);
        if power % 2 == 0 {
            for (slot, hedge) in self.subtrees[child.0].slots.clone() {
                let input = self.edge(node, NetworkEdge::Slot(slot), Flow::Sink);
                self.slot_connections.push((input, hedge));
            }
        } else {
            self.subtrees[node.0].slots = self.subtrees[child.0].slots.clone();
        }
        self.close(node)
    }
    pub(super) fn function(&mut self, child: NodeIndex, function: Symbol) -> NodeIndex {
        let slots = self.subtrees[child.0].slots.clone();
        let node = self.node(
            NetworkNode::Op(NetworkOp::Function(function)),
            Self::slot_state(&slots),
            vec![child],
        );
        self.edge(node, NetworkEdge::Head, Flow::Sink);
        self.subtrees[node.0].slots = slots;
        self.close(node)
    }
    pub(super) fn finish(self, root: NodeIndex) -> Network<S, K, Symbol, Aind> {
        let mut nodes = Vec::new();
        let mut stack = vec![root];
        while let Some(node) = stack.pop() {
            nodes.push(node);
            stack.extend(self.subtrees[node.0].children.iter().rev());
        }
        let retained_nodes = nodes.len();
        let (mut builder, hedges) = self
            .builder
            .expect("graph emission requested")
            .retain_nodes_in_order(&nodes);
        // Retention only has to permute dangling Head endpoints. Their known
        // seams connect here, before edge storage, in the same operand order.
        for &node in &nodes {
            let subtree = &self.subtrees[node.0];
            for (position, &child) in subtree.children.iter().enumerate() {
                let input = if subtree.reversed_head_inputs {
                    subtree.children.len() - 1 - position
                } else {
                    position
                };
                let input = Hedge(subtree.hedges.start + 1 + input);
                let output = Hedge(self.subtrees[child.0].hedges.start);
                builder
                    .connect_identities(
                        hedges[input].expect("retained operation input"),
                        hedges[output].expect("retained child head"),
                        NetworkGraph::<K, Symbol, Aind>::join_heads,
                    )
                    .expect("new expression heads are dangling");
            }
        }
        let mut graph = NetworkGraph::from(builder);
        for (index, &old) in nodes.iter().enumerate() {
            let subtree = &self.subtrees[old.0];
            if let Some(slots) = &subtree.tensor_slots {
                graph
                    .set_tensor_slot_order(NodeIndex(index), slots, slots)
                    .expect("new tensor crown");
                if let Some(layout) = &subtree.layout {
                    graph
                        .set_logical_layout(NodeIndex(index), layout)
                        .expect("checked tensor layout");
                }
            }
        }
        // Slot seams stay deferred until both endpoint layouts are installed.
        // Unlike structural heads, dual tensor ports retain distinct witnesses.
        for (source, sink) in self.slot_connections {
            match (hedges[source], hedges[sink]) {
                (Some(source), Some(sink)) => graph.connect_identities(source, sink),
                (None, None) => {}
                _ => unreachable!("a parsed subtree owns both seam endpoints"),
            }
        }
        let mut tensors = Vec::new();
        let mut scalars = Vec::new();
        let mut tensor_ids = vec![None; self.store.n_tensors()];
        let mut scalar_ids = vec![None; self.store.n_scalars()];
        for i in 0..retained_nodes {
            if let NetworkNode::Leaf(leaf) = &mut graph.graph[NodeIndex(i)] {
                leaf.map_tensor_refs(|id| {
                    *tensor_ids[id].get_or_insert_with(|| {
                        let next = tensors.len();
                        tensors.push(id);
                        next
                    })
                });
                leaf.map_scalar_refs(|scalar| {
                    scalar.map_index(|id| {
                        *scalar_ids[id].get_or_insert_with(|| {
                            let next = scalars.len();
                            scalars.push(id);
                            next
                        })
                    })
                });
            }
        }
        Network {
            graph,
            store: self.store.retain_in_order(&tensors, &scalars),
            state: self.subtrees[root.0].state,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        network::{
            graph::{NAdd, NMul},
            store::NetworkStore,
        },
        structure::{
            OrderedStructure,
            abstract_index::AbstractIndex,
            representation::{Lorentz, Minkowski},
        },
    };
    use linnet::half_edge::EdgeAccessors;

    type Store = NetworkStore<OrderedStructure, i32>;
    type Net = Network<Store, i8, Symbol>;
    type Builder = Construction<Store, i8, AbstractIndex>;
    type Pair = (NodeIndex, Net);

    fn tensor(builder: &mut Builder, slots: &[LibrarySlot<AbstractIndex>]) -> Pair {
        let structure = Canonicalized::<OrderedStructure>::from_iter(slots.iter().copied());
        let layout = structure.layout().clone();
        let payload = structure.into_canonical();
        let node = builder.tensor(payload.clone(), layout.clone()).unwrap();
        let mut old = Net::from_tensor(payload);
        old.graph
            .set_logical_layout(old.graph.graph.node_id(old.graph.head()), &layout)
            .unwrap();
        (node, old)
    }

    fn scalar(builder: &mut Builder, value: i32) -> Pair {
        (builder.scalar(value), Net::from_scalar(value))
    }

    fn product(builder: &mut Builder, children: Vec<Pair>) -> Pair {
        let (nodes, nets): (Vec<_>, Vec<_>) = children.into_iter().unzip();
        let mut nets = nets.into_iter();
        (builder.product(nodes), nets.next().unwrap().n_mul(nets))
    }

    fn sum(builder: &mut Builder, children: Vec<Pair>) -> Pair {
        let (nodes, nets): (Vec<_>, Vec<_>) = children.into_iter().unzip();
        let mut nets = nets.into_iter();
        (
            builder.sum(nodes).unwrap(),
            nets.next().unwrap().n_add(nets),
        )
    }

    fn identical(actual: Net, expected: Net) {
        actual.graph.graph.check().unwrap();
        assert_eq!(actual.state, expected.state);
        assert_eq!(actual.store, expected.store);
        assert_eq!(actual.graph.slot_order, expected.graph.slot_order);
        assert_eq!(
            actual.graph.logical_slot_order,
            expected.graph.logical_slot_order
        );
        let actual = actual.graph.graph;
        let expected = expected.graph.graph;
        assert_eq!(actual.n_nodes(), expected.n_nodes());
        assert_eq!(actual.n_hedges(), expected.n_hedges());
        for i in 0..actual.n_nodes() {
            assert_eq!(actual[NodeIndex(i)], expected[NodeIndex(i)]);
            assert_eq!(
                actual.iter_crown(NodeIndex(i)).collect::<Vec<_>>(),
                expected.iter_crown(NodeIndex(i)).collect::<Vec<_>>()
            );
        }
        let mut edge_ids = BTreeMap::new();
        let mut inverse = BTreeMap::new();
        for i in 0..actual.n_hedges() {
            let h = Hedge(i);
            assert_eq!(actual.inv(h), expected.inv(h), "involution at {i}");
            assert_eq!(actual.flow(h), expected.flow(h), "flow at {i}");
            assert_eq!(
                actual.orientation(h),
                expected.orientation(h),
                "orientation at {i}"
            );
            let edge: &NetworkEdge<_> = &actual[[&h]];
            assert_eq!(edge, &expected[[&h]], "edge at {i}");
            assert_eq!(actual.node_id(h), expected.node_id(h), "node at {i}");
            if let Some(previous) = edge_ids.insert(expected[&h], actual[&h]) {
                assert_eq!(previous, actual[&h]);
            }
            if let Some(previous) = inverse.insert(actual[&h], expected[&h]) {
                assert_eq!(previous, expected[&h]);
            }
        }
        assert_eq!(edge_ids.len(), expected.iter_edges().count());
        assert_eq!(inverse.len(), actual.iter_edges().count());
        for (pair, id, data) in expected.iter_edges() {
            let (actual_pair, _, actual_data) = actual
                .iter_edges()
                .find(|(_, actual_id, _)| *actual_id == edge_ids[&id])
                .unwrap();
            assert_eq!(actual_pair, pair);
            assert_eq!(actual_data, data);
        }
    }

    #[test]
    fn structure_observation_matches_retained_graph_port_order() {
        let slot = |i| Minkowski {}.new_slot(4, i).to_lib();
        let function = symbolica::symbol!("structure_only_function");
        for observe in [false, true] {
            for exponent in [-2, 2, 3] {
                let mut builder = if observe {
                    Builder::structure_only()
                } else {
                    Builder::new()
                };
                let unused = tensor(&mut builder, &[slot(11), slot(11)]);
                let _ = builder.power(unused.0, 2);
                let left = tensor(&mut builder, &[slot(7), slot(3)]);
                let right = tensor(&mut builder, &[slot(7), slot(3)]);
                let (sum, old) = sum(&mut builder, vec![left, right]);
                let sum = builder.function(sum, function);
                let sum = builder.power(sum, exponent);
                let other = tensor(&mut builder, &[slot(3), slot(19)]);
                let (root, expected) = product(
                    &mut builder,
                    vec![(sum, old.fun(function).pow(exponent)), other],
                );
                assert_eq!(builder.slots(root), expected.graph.dangling_indices());
                assert_eq!(builder.emits_graph(), !observe);
                if !observe {
                    identical(builder.finish(root), expected);
                }
            }
        }
    }

    #[test]
    fn arena_preserves_existing_assembly_and_discards_only_unused_payloads() {
        let slot = |i| Minkowski {}.new_slot(4, i).to_lib();
        let upper = Lorentz {}.new_slot(4, 47).to_lib();
        let lower = upper.dual();
        let fixtures = [
            vec![vec![]],
            vec![
                vec![slot(1), slot(2)],
                vec![slot(2), slot(3)],
                vec![slot(3), slot(1)],
            ],
            vec![vec![slot(1)], vec![slot(1)], vec![slot(1)]],
            vec![vec![slot(5), slot(5), slot(5), slot(5)], vec![slot(7)]],
            vec![vec![upper, lower], vec![upper], vec![lower]],
            vec![vec![slot(19), slot(11)], vec![slot(11), slot(23)]],
        ];
        for slots in fixtures {
            let mut builder = Builder::new();
            // Parsing may evaluate and then precontract scalar temporaries.
            // Their nodes and payloads must not survive the finished root.
            let ignored_a = scalar(&mut builder, 101);
            let ignored_b = scalar(&mut builder, 103);
            let _ignored = product(&mut builder, vec![ignored_a, ignored_b]);
            let operands = slots
                .iter()
                .map(|s| tensor(&mut builder, s))
                .collect::<Vec<_>>();
            let coefficient = scalar(&mut builder, 7);
            let mut reordered = vec![coefficient];
            reordered.extend(operands);
            let (root, expected) = product(&mut builder, reordered);
            identical(builder.finish(root), expected);
        }
        for power in [-2, 2, 3] {
            let mut builder = Builder::new();
            let a = tensor(&mut builder, &[slot(3), slot(1)]);
            let b = tensor(&mut builder, &[slot(3), slot(1)]);
            let (node, old) = sum(&mut builder, vec![a, b]);
            let node = builder.power(node, power);
            let function = symbolica::symbol!("arena_assembly_function");
            let node = builder.function(node, function);
            identical(builder.finish(node), old.pow(power).fun(function));
        }
        // Retain a sole leaf even when earlier discarded subtrees had edges.
        let mut builder = Builder::new();
        let unused = tensor(&mut builder, &[slot(1), slot(1)]);
        let _ = builder.power(unused.0, 2);
        let (root, old) = scalar(&mut builder, 31);
        identical(builder.finish(root), old);
    }
    #[test]
    fn arena_connects_nested_head_seams_without_changing_branch_order() {
        let function = symbolica::symbol!("arena_nested_head_function");
        for exponent in [-2, 2, 3] {
            let mut builder = Builder::new();
            let mut branches = Vec::new();
            for coefficient in [2, 3, 5] {
                let left = scalar(&mut builder, coefficient);
                let right = scalar(&mut builder, coefficient + 1);
                let (node, old) = sum(&mut builder, vec![left, right]);
                let node = builder.function(node, function);
                let node = builder.power(node, exponent);
                let value = (node, old.fun(function).pow(exponent));
                let tensor = tensor(&mut builder, &[]);
                branches.push(product(&mut builder, vec![value, tensor]));
            }
            let (root, old) = sum(&mut builder, branches);
            identical(builder.finish(root), old);
        }
    }

    #[test]
    fn arena_retains_a_scalar_child_and_discards_failed_subtrees() {
        let mut builder = Builder::new();
        let child = scalar(&mut builder, 29);
        let root = child.0;
        let expected = child.1.clone();
        let other = scalar(&mut builder, 31);
        let _discarded = product(&mut builder, vec![child, other]);
        identical(builder.finish(root), expected);

        let mut builder = Builder::new();
        let slot = |i| Lorentz {}.new_slot(4, i).to_lib();
        let child = tensor(&mut builder, &[slot(1), slot(2)]);
        let invalid = tensor(&mut builder, &[slot(1), slot(3)]);
        assert!(matches!(
            builder.sum(vec![child.0, invalid.0]),
            Err(TensorNetworkError::IncompatibleSummand(_))
        ));
        // The failed sum and all of its children are unreachable. Retaining a
        // tensor child alone would leave a Slot seam crossing that boundary,
        // which the existing finish contract deliberately does not admit.
        let (root, expected) = scalar(&mut builder, 37);
        identical(builder.finish(root), expected);
    }

    #[test]
    fn sum_rejects_missing_extra_and_duplicate_ports() {
        let slot = |index| Lorentz {}.new_slot(4, index).to_lib();
        for right in [vec![], vec![1], vec![1, 2, 3], vec![1, 3], vec![1, 1]] {
            let mut builder = Builder::new();
            let left = tensor(&mut builder, &[slot(1), slot(2)]).0;
            let right = tensor(
                &mut builder,
                &right.into_iter().map(slot).collect::<Vec<_>>(),
            )
            .0;
            assert!(matches!(
                builder.sum(vec![left, right]),
                Err(TensorNetworkError::IncompatibleSummand(_))
            ));
        }
    }

    #[test]
    fn arena_keeps_function_execution_order_and_results() {
        crate::structure::representation::initialize();
        use crate::{
            network::{
                ExecutionResult, Sequential, SequentialExtract, SequentialRef, SmallestDegree,
                library::{
                    DummyKey, DummyLibrary, DummyLibraryTensor, FunctionLibrary,
                    FunctionLibraryError,
                },
            },
            tensors::data::DenseTensor,
        };
        use std::cell::RefCell;
        type Tensor = DenseTensor<f64, OrderedStructure>;
        type Values = NetworkStore<Tensor, f64>;
        type NumericNet = Network<Values, DummyKey, Symbol>;
        type LibTensor = DummyLibraryTensor<Tensor>;
        type Lib = DummyLibrary<Tensor, DummyKey>;
        #[derive(Default)]
        struct Calls(RefCell<Vec<(Symbol, Vec<f64>)>>);
        impl FunctionLibrary<Tensor, f64> for Calls {
            type Key = Symbol;
            fn apply(
                &self,
                key: &Symbol,
                mut value: Tensor,
            ) -> Result<Tensor, FunctionLibraryError<Symbol>> {
                let mut calls = self.0.borrow_mut();
                calls.push((*key, value.data.clone()));
                for component in &mut value.data {
                    *component += calls.len() as f64;
                }
                Ok(value)
            }
            fn apply_scalar(
                &self,
                key: &Symbol,
                value: f64,
            ) -> Result<f64, FunctionLibraryError<Symbol>> {
                let mut calls = self.0.borrow_mut();
                calls.push((*key, vec![value]));
                Ok(value + calls.len() as f64)
            }
        }
        let mut builder = Construction::<Values, DummyKey, AbstractIndex>::new();
        let mut branches = Vec::new();
        let mut old = Vec::new();
        for (index, hook) in [
            symbolica::symbol!("arena_callback_a"),
            symbolica::symbol!("arena_callback_b"),
        ]
        .into_iter()
        .enumerate()
        {
            let slot = Minkowski {}.new_slot(2, index).to_lib();
            let structure = OrderedStructure::new(vec![slot]);
            let tensor = |values| {
                DenseTensor::from_storage_data(values, structure.canonical().clone()).unwrap()
            };
            let left = tensor(vec![1.0, 2.0]);
            let right = tensor(vec![3.0, 4.0]);
            let left_node = builder
                .tensor(left.clone(), structure.layout().clone())
                .unwrap();
            let left_node = builder.function(left_node, hook);
            let right_node = builder
                .tensor(right.clone(), structure.layout().clone())
                .unwrap();
            branches.push(builder.product(vec![left_node, right_node]));
            old.push(
                NumericNet::from_tensor(left)
                    .fun(hook)
                    .n_mul([NumericNet::from_tensor(right)]),
            );
        }
        let root = builder.sum(branches).unwrap();
        let hook = symbolica::symbol!("arena_callback_result");
        let root = builder.function(root, hook);
        let actual = builder.finish(root);
        let mut old = old.into_iter();
        let expected = old.next().unwrap().n_add(old).fun(hook);
        let lib = Lib::new();
        macro_rules! compare {
            ($strategy:ty) => {{
                let (mut actual, mut expected) = (actual.clone(), expected.clone());
                let (actual_calls, expected_calls) = (Calls::default(), Calls::default());
                actual
                    .execute::<$strategy, SmallestDegree, LibTensor, Lib, Calls>(
                        &lib,
                        &actual_calls,
                    )
                    .unwrap();
                expected
                    .execute::<$strategy, SmallestDegree, LibTensor, Lib, Calls>(
                        &lib,
                        &expected_calls,
                    )
                    .unwrap();
                let ExecutionResult::Val(actual) = actual.result_scalar().unwrap() else {
                    panic!("closed scalar result")
                };
                let ExecutionResult::Val(expected) = expected.result_scalar().unwrap() else {
                    panic!("closed scalar result")
                };
                assert_eq!(actual, expected);
                assert_eq!(actual_calls.0.borrow().len(), 3);
                assert_eq!(*actual_calls.0.borrow(), *expected_calls.0.borrow());
            }};
        }
        compare!(Sequential);
        compare!(SequentialRef);
        compare!(SequentialExtract);
    }
}
