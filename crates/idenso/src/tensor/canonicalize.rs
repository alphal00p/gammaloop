//! Signed tensor canonicalization on Spenso's symbolic network.
//!
//! Tensor syntax, operation grouping, and contraction topology belong to the
//! network parser. This module only projects a parsed monomial into Graphica's
//! incidence graph so that a negative automorphism of antisymmetric tensor
//! slots can be detected. Symbolica still owns the final dummy-index names and
//! expression rebuilding.

use std::collections::BTreeMap;

use linnet::half_edge::{
    NodeIndex,
    involution::{Hedge, HedgePair},
    subgraph::{ModifySubSet, SuBitGraph},
    tree::SimpleTraversalTree,
};
use linnet::permutation::Permutation;
use linnet::tree::child_pointer::ParentChildStore;
use spenso::{
    network::{
        graph::{NetworkEdge, NetworkLeaf, NetworkNode, NetworkOp},
        store::NetworkStoreAccess,
    },
    structure::{
        HasName, OrderedStructure, TensorStructure,
        representation::{LibraryRep, LibrarySlot, RepName, Representation},
        slot::{AbsInd, DualSlotTo, IsAbstractSlot, ParseableAind},
    },
};
use symbolica::{
    atom::{Atom, AtomView, Symbol},
    graph::{Graph, HiddenData},
};

use super::{SymbolicNet, SymbolicTensor};

const CYCLE_SYMMETRY_MARKER: usize = usize::MAX - 2;

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct TensorColor {
    head: Symbol,
    arguments: Vec<Atom>,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum TensorGraphNode<Aind> {
    Product,
    Tensor(TensorColor),
    /// A nonlinear function, sum, or unresolved leaf, identified per occurrence.
    Opaque(NodeIndex),
    Slot(SlotColor<Aind>),
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum SlotColor<Aind> {
    Internal(Representation<LibraryRep>),
    External(LibrarySlot<Aind>),
}

type TensorGraph<Aind> =
    Graph<TensorGraphNode<Aind>, HiddenData<(usize, Option<Representation<LibraryRep>>), usize>>;

struct SlotCopies<Aind> {
    slot: LibrarySlot<Aind>,
    power_path: Vec<NodeIndex>,
    vertices: Vec<usize>,
}

type SlotsByHedge<Aind> = BTreeMap<Hedge, SlotCopies<Aind>>;

/// Prune signed zeros on the already parsed network, without distributing sums.
///
/// Each mutation removes an operation/subtree, or replaces a nonzero tensor by
/// zero. Rebuild only the graph traversal after compacting node/hedge storage;
/// tensor syntax and interfaces are never parsed again. Nonlinear functions
/// and nonpositive powers retain their opaque boundary.
pub(super) fn remove_antisymmetric_zero_terms<Aind: AbsInd + ParseableAind>(
    network: &mut SymbolicNet<Aind>,
) -> bool {
    if !contains_antisymmetric_tensor(network) {
        return false;
    }
    let mut changed = false;
    loop {
        let tree: SimpleTraversalTree<ParentChildStore<()>> = network.graph.expr_tree().cast();
        let root = network.graph.graph.node_id(network.graph.head());
        if !prune_one(network, &tree, root) {
            return changed;
        }
        changed = true;
    }
}

fn is_zero<Aind: AbsInd>(network: &SymbolicNet<Aind>, node: NodeIndex) -> bool {
    match &network.graph.graph[node] {
        NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => {
            network.store.tensors[*index].expression.is_zero()
        }
        NetworkNode::Leaf(NetworkLeaf::Scalar(index)) => network.store.scalar_ref(*index).is_zero(),
        _ => false,
    }
}

fn prune_one<Aind: AbsInd + ParseableAind>(
    network: &mut SymbolicNet<Aind>,
    tree: &SimpleTraversalTree<ParentChildStore<()>>,
    node: NodeIndex,
) -> bool {
    if matches!(
        network.graph.graph[node],
        NetworkNode::Op(NetworkOp::Function(_) | NetworkOp::Power(..=0))
    ) || is_zero(network, node)
    {
        return false;
    }
    let children = tree
        .iter_children(node, network.graph.graph.as_ref())
        .collect::<Vec<_>>();
    for child in &children {
        if prune_one(network, tree, *child) {
            return true;
        }
    }
    let zeros = children
        .iter()
        .filter(|child| is_zero(network, **child))
        .count();
    match network.graph.graph[node] {
        NetworkNode::Op(NetworkOp::Sum) if zeros != 0 => {
            if zeros == children.len() {
                replace_with_zero(network, tree, node);
            } else {
                // A zero child has already become a leaf. Remove both halves
                // of its head and Sum-input slots, leaving the Sum output and
                // every surviving arm's incidence unchanged.
                let mut removed: SuBitGraph = network.graph.graph.empty_subgraph();
                for child in children.iter().filter(|child| is_zero(network, **child)) {
                    for hedge in network.graph.graph.iter_crown(*child) {
                        removed.add(hedge);
                        removed.add(network.graph.graph.inv(hedge));
                    }
                }
                network.graph.delete(&removed);
            }
            return true;
        }
        NetworkNode::Op(NetworkOp::Sum) if children.len() == 1 => {
            return remove_singleton_sum(network, tree, node, children[0]);
        }
        NetworkNode::Op(NetworkOp::Product | NetworkOp::Neg | NetworkOp::Power(1..))
            if zeros != 0 =>
        {
            replace_with_zero(network, tree, node);
            return true;
        }
        _ => {}
    }
    if has_odd_automorphism(network, tree, node) {
        replace_with_zero(network, tree, node);
        return true;
    }
    false
}

fn replace_with_zero<Aind: AbsInd>(
    network: &mut SymbolicNet<Aind>,
    tree: &SimpleTraversalTree<ParentChildStore<()>>,
    root: NodeIndex,
) {
    let nodes = tree
        .iter_preorder_tree_nodes(network.graph.graph.as_ref(), root)
        .collect::<Vec<_>>();
    let inside = nodes
        .iter()
        .copied()
        .collect::<std::collections::BTreeSet<_>>();
    let mut boundary = Vec::new();
    for node in &nodes {
        for hedge in network.graph.graph.iter_crown(*node) {
            let NetworkEdge::Slot(edge_slot) = network.graph.graph[[&hedge]] else {
                continue;
            };
            let other = network.graph.graph.inv(hedge);
            if other != hedge && inside.contains(&network.graph.graph.node_id(other)) {
                continue;
            }
            // Sewing a dual pair can retain the opposite endpoint's descriptor.
            // Recover this endpoint from its already established storage axis.
            let slot = match &network.graph.graph[*node] {
                NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => network.store.tensors[*index]
                    .structure
                    .external_structure_iter()
                    .nth(usize::from(network.graph.slot_order[hedge.0]))
                    .expect("a parsed tensor endpoint has a storage axis"),
                NetworkNode::Op(NetworkOp::Sum) => {
                    // A Sum output forwards each arm's identical interface.
                    // Its input seam still carries that descriptor even when
                    // an outer product retained the dual output descriptor.
                    network
                        .graph
                        .graph
                        .iter_crown(*node)
                        .find_map(|input| {
                            let NetworkEdge::Slot(slot) = network.graph.graph[[&input]] else {
                                return None;
                            };
                            let other = network.graph.graph.inv(input);
                            (other != input
                                && inside.contains(&network.graph.graph.node_id(other))
                                && (slot == edge_slot || slot.matches(&edge_slot)))
                            .then_some(slot)
                        })
                        .expect("a parsed Sum output has a matching input port")
                }
                _ => edge_slot,
            };
            boundary.push((hedge, slot));
        }
    }
    let layout = OrderedStructure::new(boundary.iter().map(|(_, slot)| *slot).collect());
    let slots = layout.canonical().external_structure();
    for (position, slot) in slots.iter().enumerate() {
        let found = boundary
            .iter()
            .position(|(_, current)| current == slot)
            .expect("the canonical zero interface permutes its boundary ports");
        let (hedge, _) = boundary.remove(found);
        network.graph.slot_order[hedge.0] =
            u8::try_from(position).expect("parsed slot position fits u8");
    }
    let index = network.store.tensors.len();
    network.graph.identify_nodes_without_self_edges(
        &nodes,
        NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)),
    );
    // Identification intentionally preserves pre-existing self loops. A zero
    // replacement has already consumed those internal contractions, so remove
    // their paired ports too. Recover its occurrence after node compaction.
    let node = network
        .graph
        .graph
        .iter_nodes()
        .find_map(|(node, _, data)| {
            matches!(data, NetworkNode::Leaf(NetworkLeaf::LocalTensor(found)) if *found == index)
                .then_some(node)
        })
        .expect("the replacement occurrence survives graph compaction");
    let mut internal: SuBitGraph = network.graph.graph.empty_subgraph();
    for hedge in network.graph.graph.iter_crown(node) {
        let other = network.graph.graph.inv(hedge);
        if other != hedge && network.graph.graph.node_id(other) == node {
            internal.add(hedge);
        }
    }
    network.graph.delete(&internal);
    // Keep the proven boundary until the existing execution owner consumes it;
    // deleting an open zero's ports here would corrupt the enclosing Sum.
    network.store.tensors.push(SymbolicTensor {
        proofs: Default::default(),
        expression: Atom::zero(),
        structure: layout.into_canonical(),
        is_metric: false,
        is_composite: true,
    });
}

fn remove_singleton_sum<Aind: AbsInd>(
    network: &mut SymbolicNet<Aind>,
    tree: &SimpleTraversalTree<ParentChildStore<()>>,
    node: NodeIndex,
    child: NodeIndex,
) -> bool {
    let descendants = tree
        .iter_preorder_tree_nodes(network.graph.graph.as_ref(), child)
        .collect::<std::collections::BTreeSet<_>>();
    let crown = network.graph.graph.iter_crown(node).collect::<Vec<_>>();
    let mut seams = Vec::new();
    for input in &crown {
        let inside = network.graph.graph.inv(*input);
        if inside == *input || !descendants.contains(&network.graph.graph.node_id(inside)) {
            continue;
        }
        let input_data = network.graph.graph[[input]];
        let Some(output) = crown.iter().copied().find(|output| {
            let outside = network.graph.graph.inv(*output);
            let same_port = match (network.graph.graph[[output]], input_data) {
                (NetworkEdge::Head, NetworkEdge::Head) => true,
                (NetworkEdge::Slot(output), NetworkEdge::Slot(input)) => {
                    output == input || output.matches(&input)
                }
                _ => false,
            };
            same_port
                && (outside == *output
                    || !descendants.contains(&network.graph.graph.node_id(outside)))
        }) else {
            return false;
        };
        seams.push((output, *input));
    }
    if seams.len() * 2 != crown.len() {
        return false;
    }
    // Bypass the existing Sum seams, retaining the outside edge's descriptor
    // and orientation. No tensor storage axis or source logical witness moves.
    for (output, input) in seams {
        let outside = network.graph.graph.inv(output);
        let inside = network.graph.graph.inv(input);
        let input_data = network
            .graph
            .graph
            .get_edge_data_full(input)
            .map(|data| *data);
        network
            .graph
            .graph
            .split_edge(input, input_data)
            .expect("Sum input is paired");
        if outside != output {
            let output_data = network
                .graph
                .graph
                .get_edge_data_full(output)
                .map(|data| *data);
            network
                .graph
                .graph
                .split_edge(output, output_data)
                .expect("Sum output is paired");
            network
                .graph
                .graph
                .connect_identities(outside, inside, |flow, data, _, _| (flow, data));
        }
    }
    let mut removed: SuBitGraph = network.graph.graph.empty_subgraph();
    for hedge in crown {
        removed.add(hedge);
    }
    network.graph.delete(&removed);
    true
}

fn contains_antisymmetric_tensor<Aind: AbsInd>(network: &SymbolicNet<Aind>) -> bool {
    network.graph.graph.iter_nodes().any(|(_, _, node)| {
        let NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) = node else {
            return false;
        };
        network.store.tensors[*index]
            .name()
            .is_some_and(|head| head.is_antisymmetric())
    })
}

fn has_odd_automorphism<Aind: AbsInd + ParseableAind>(
    network: &SymbolicNet<Aind>,
    tree: &SimpleTraversalTree<ParentChildStore<()>>,
    root: NodeIndex,
) -> bool {
    let Some(graph) = project_network(network, tree, root) else {
        return false;
    };
    if !graph.nodes().iter().any(|node| {
        matches!(&node.data,
        TensorGraphNode::Tensor(color) if color.head.is_antisymmetric())
    }) {
        return false;
    }
    // An odd stabilizer generator proves that the monomial equals its negative.
    let canonical = graph.canonize();
    canonical
        .orbit_generators
        .iter()
        .any(|generator| generator_is_odd(&canonical.graph, generator))
}

fn project_network<Aind: AbsInd + ParseableAind>(
    network: &SymbolicNet<Aind>,
    tree: &SimpleTraversalTree<ParentChildStore<()>>,
    root: NodeIndex,
) -> Option<TensorGraph<Aind>> {
    let mut graph = Graph::new();
    let mut slot_copies = BTreeMap::new();
    project_expression(network, tree, root, &mut graph, &mut slot_copies)?;
    connect_slots(network, &mut graph, &slot_copies)?;
    Some(graph)
}

fn project_expression<Aind: AbsInd + ParseableAind>(
    network: &SymbolicNet<Aind>,
    tree: &SimpleTraversalTree<ParentChildStore<()>>,
    node: NodeIndex,
    graph: &mut TensorGraph<Aind>,
    slot_copies: &mut SlotsByHedge<Aind>,
) -> Option<usize> {
    match &network.graph.graph[node] {
        NetworkNode::Op(NetworkOp::Product) => {
            let header = graph.add_node(TensorGraphNode::Product);
            for child in tree.iter_children(node, network.graph.graph.as_ref()) {
                let child = project_expression(network, tree, child, graph, slot_copies)?;
                add_incidence(graph, header, child);
            }
            Some(header)
        }
        NetworkNode::Op(NetworkOp::Power(power)) if *power > 0 => {
            let child = tree
                .iter_children(node, network.graph.graph.as_ref())
                .next()?;
            let header = graph.add_node(TensorGraphNode::Product);
            for _ in 0..*power {
                let copy = project_expression(network, tree, child, graph, slot_copies)?;
                add_incidence(graph, header, copy);
            }
            Some(header)
        }
        NetworkNode::Op(NetworkOp::Function(_)) => {
            Some(graph.add_node(TensorGraphNode::Opaque(node)))
        }
        NetworkNode::Op(NetworkOp::Neg) => {
            let child = tree
                .iter_children(node, network.graph.graph.as_ref())
                .next()?;
            project_expression(network, tree, child, graph, slot_copies)
        }
        // Sums stay factorized, and function arguments stay opaque. An
        // automorphism confined to one summand does not multiply its parent,
        // but independent factors remain analyzable.
        NetworkNode::Op(NetworkOp::Sum) => Some(graph.add_node(TensorGraphNode::Opaque(node))),
        NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => add_tensor(
            network,
            node,
            &network.store.tensors[*index],
            tree,
            graph,
            slot_copies,
        ),
        NetworkNode::Leaf(_) => Some(graph.add_node(TensorGraphNode::Opaque(node))),
        NetworkNode::Op(NetworkOp::Power(_)) => None,
    }
}

fn add_tensor<Aind: AbsInd + ParseableAind>(
    network: &SymbolicNet<Aind>,
    network_node: NodeIndex,
    tensor: &SymbolicTensor<OrderedStructure<LibraryRep, Aind>>,
    tree: &SimpleTraversalTree<ParentChildStore<()>>,
    graph: &mut TensorGraph<Aind>,
    slot_copies: &mut SlotsByHedge<Aind>,
) -> Option<usize> {
    let slots = tensor
        .structure
        .external_structure_iter()
        .collect::<Vec<_>>();
    let color = tensor_color(tensor, &slots)?;
    let head = color.head;
    let unordered = head.is_symmetric() || head.is_antisymmetric() || head.is_cyclesymmetric();
    let header = graph.add_node(TensorGraphNode::Tensor(color));
    let mut network_slots = network
        .graph
        .graph
        .iter_crown(network_node)
        .filter(|hedge| matches!(network.graph.graph[[hedge]], NetworkEdge::Slot(_)))
        .collect::<Vec<_>>();
    network_slots.sort_unstable_by_key(|hedge| network.graph.slot_order[hedge.0]);
    if network_slots.len() != slots.len() {
        return None;
    }
    let power_path = tree
        .ancestor_iter_node(network_node, network.graph.graph.as_ref())
        .skip(1)
        .filter(|ancestor| {
            matches!(
                network.graph.graph[*ancestor],
                NetworkNode::Op(NetworkOp::Power(1..))
            )
        })
        .collect::<Vec<_>>();

    let mut cycle_slots = Vec::with_capacity(slots.len());
    for (position, (slot, network_slot)) in slots.iter().copied().zip(network_slots).enumerate() {
        let vertex = graph.add_node(TensorGraphNode::Slot(SlotColor::External(slot)));
        let visible_position = if unordered { 0 } else { position };
        graph
            .add_edge(
                header,
                vertex,
                true,
                HiddenData::new((visible_position, None), position),
            )
            .unwrap();
        slot_copies
            .entry(network_slot)
            .or_insert_with(|| SlotCopies {
                slot,
                power_path: power_path.clone(),
                vertices: vec![],
            })
            .vertices
            .push(vertex);
        cycle_slots.push(vertex);
    }

    if head.is_cyclesymmetric() && cycle_slots.len() > 1 {
        for (&left, &right) in cycle_slots
            .iter()
            .zip(cycle_slots.iter().cycle().skip(1))
            .take(cycle_slots.len())
        {
            graph
                .add_edge(
                    left,
                    right,
                    true,
                    HiddenData::new((CYCLE_SYMMETRY_MARKER, None), 0),
                )
                .unwrap();
        }
    }
    Some(header)
}

fn tensor_color<Aind: AbsInd + ParseableAind>(
    tensor: &SymbolicTensor<OrderedStructure<LibraryRep, Aind>>,
    slots: &[LibrarySlot<Aind>],
) -> Option<TensorColor> {
    let AtomView::Fun(function) = tensor.expression.as_view() else {
        return None;
    };
    let mut slot_atoms = slots
        .iter()
        .map(IsAbstractSlot::to_atom)
        .collect::<Vec<_>>();
    let mut arguments = Vec::new();
    for argument in function.iter() {
        if let Some(position) = slot_atoms
            .iter()
            .position(|slot| slot.as_view() == argument)
        {
            let _ = slot_atoms.swap_remove(position);
        } else {
            arguments.push(argument.to_owned());
        }
    }
    if !slot_atoms.is_empty() {
        return None;
    }

    Some(TensorColor {
        head: function.get_symbol(),
        arguments,
    })
}

fn add_incidence<Aind: AbsInd>(graph: &mut TensorGraph<Aind>, parent: usize, child: usize) {
    graph
        .add_edge(parent, child, true, HiddenData::new((0, None), 0))
        .unwrap();
}

fn connect_slots<Aind: AbsInd>(
    network: &SymbolicNet<Aind>,
    graph: &mut TensorGraph<Aind>,
    slot_copies: &SlotsByHedge<Aind>,
) -> Option<()> {
    for (pair, _, data) in network.graph.graph.iter_edges() {
        let NetworkEdge::Slot(slot) = data.data else {
            continue;
        };
        let HedgePair::Paired { source, sink } = pair else {
            continue;
        };
        let source_copies = slot_copies.get(&source);
        let sink_copies = slot_copies.get(&sink);
        let group = if slot.rep_name().is_dual() {
            slot.rep().dual()
        } else {
            slot.rep()
        };

        match (source_copies, sink_copies) {
            (Some(source), Some(sink))
                if source.vertices.len() == sink.vertices.len()
                    && (source.vertices.len() == 1 || source.power_path == sink.power_path) =>
            {
                for (&source_vertex, &sink_vertex) in source.vertices.iter().zip(&sink.vertices) {
                    connect_slot_pair(
                        graph,
                        source_vertex,
                        sink_vertex,
                        source.slot.rep(),
                        sink.slot.rep(),
                        &group,
                    );
                }
            }
            (Some(copies), None) if is_power_node(network, sink) => {
                let [left, right] = copies.vertices.as_slice() else {
                    return None;
                };
                connect_slot_pair(
                    graph,
                    *left,
                    *right,
                    copies.slot.rep(),
                    copies.slot.rep(),
                    &group,
                );
            }
            (None, Some(copies)) if is_power_node(network, source) => {
                let [left, right] = copies.vertices.as_slice() else {
                    return None;
                };
                connect_slot_pair(
                    graph,
                    *left,
                    *right,
                    copies.slot.rep(),
                    copies.slot.rep(),
                    &group,
                );
            }
            // An opaque subtree can own the other endpoint. Keeping this
            // endpoint fixed is conservative and still permits independent factors.
            (Some(_), None) | (None, Some(_)) => {}
            (None, None) => {}
            _ => return None,
        }
    }
    if slot_copies.iter().any(|(hedge, copies)| {
        copies.vertices.len() != 1 && network.graph.graph.as_ref().is_identity(*hedge)
    }) {
        return None;
    }
    Some(())
}

fn is_power_node<Aind: AbsInd>(network: &SymbolicNet<Aind>, hedge: Hedge) -> bool {
    matches!(
        network.graph.graph[network.graph.graph.node_id(hedge)],
        NetworkNode::Op(NetworkOp::Power(_))
    )
}

fn connect_slot_pair<Aind: AbsInd>(
    graph: &mut TensorGraph<Aind>,
    source: usize,
    sink: usize,
    source_rep: Representation<LibraryRep>,
    sink_rep: Representation<LibraryRep>,
    group: &Representation<LibraryRep>,
) {
    graph.set_node_data(
        source,
        TensorGraphNode::Slot(SlotColor::Internal(source_rep)),
    );
    graph.set_node_data(sink, TensorGraphNode::Slot(SlotColor::Internal(sink_rep)));
    graph
        .add_edge(source, sink, false, HiddenData::new((0, Some(*group)), 0))
        .unwrap();
}

fn generator_is_odd<Aind: AbsInd>(graph: &TensorGraph<Aind>, cycles: &[Vec<usize>]) -> bool {
    let mut automorphism = (0..graph.nodes().len()).collect::<Vec<_>>();
    for cycle in cycles {
        for (source, target) in cycle
            .iter()
            .zip(cycle.iter().cycle().skip(1))
            .take(cycle.len())
        {
            automorphism[*source] = *target;
        }
    }

    let mut odd = false;
    for node_index in 0..graph.nodes().len() {
        let TensorGraphNode::Tensor(color) = &graph.node(node_index).data else {
            continue;
        };
        if !color.head.is_antisymmetric() {
            continue;
        }

        let target_header = automorphism[node_index];
        let mut induced_permutation = vec![];
        for &edge_index in &graph.node(node_index).edges {
            let edge = graph.edge(edge_index);
            if !edge.directed || edge.vertices.0 != node_index {
                continue;
            }

            let target_slot = automorphism[edge.vertices.1];
            let target_position = graph
                .node(target_header)
                .edges
                .iter()
                .map(|edge_index| graph.edge(*edge_index))
                .find(|target_edge| {
                    target_edge.directed && target_edge.vertices == (target_header, target_slot)
                })
                .expect("an automorphism preserves tensor-slot incidence")
                .data
                .hidden;
            induced_permutation.push((edge.data.hidden, target_position));
        }

        induced_permutation.sort_unstable_by_key(|(source, _)| *source);
        odd ^= Permutation::from_map(
            induced_permutation
                .into_iter()
                .map(|(_, target)| target)
                .collect(),
        )
        .sign()
            < 0;
    }
    odd
}
