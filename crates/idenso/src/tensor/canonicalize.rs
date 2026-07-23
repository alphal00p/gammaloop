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
    tree::SimpleTraversalTree,
};
use linnet::permutation::Permutation;
use linnet::tree::child_pointer::ParentChildStore;
use spenso::{
    network::{
        graph::{NetworkEdge, NetworkLeaf, NetworkNode, NetworkOp},
        parsing::ParseSettings,
    },
    structure::{
        HasName, TensorStructure,
        representation::{LibraryRep, LibrarySlot, RepName, Representation},
        slot::{AbsInd, DummyAind, IsAbstractSlot, ParseableAind},
    },
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    graph::{Graph, HiddenData},
};

use crate::Cookable;

use super::{SymbolicNet, SymbolicNetParse, SymbolicTensor};

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

/// Remove Spenso tensor monomials that are fixed by a negative automorphism.
///
/// Products and positive powers are expanded only when they contain a sum.
/// Sums that remain inside functions stay opaque. If parsing fails or no term
/// vanishes, the original expression is returned so its factorization is
/// preserved.
pub(crate) fn remove_antisymmetric_zero_terms<Aind: AbsInd + DummyAind + ParseableAind>(
    expression: AtomView<'_>,
) -> Atom {
    // Structured index payloads such as `hedge(1)` are not abstract indices
    // themselves. Cook a probe for parsing, but retain the original expression.
    let cooked = expression.cook_indices();
    let Ok(network) = cooked
        .as_view()
        .parse_to_symbolic_net::<Aind>(&ParseSettings::default())
    else {
        return expression.to_owned();
    };
    if !contains_antisymmetric_tensor(&network) {
        return expression.to_owned();
    }

    let expanded = requires_term_expansion(&network).then(|| expression.expand());
    let candidate = expanded
        .as_ref()
        .map_or(expression, |expanded| expanded.as_view());
    let terms = match candidate {
        AtomView::Add(add) => add.iter().collect::<Vec<_>>(),
        _ => vec![candidate],
    };

    let mut removed = false;
    let retained = terms
        .into_iter()
        .filter(|term| {
            let cooked = term.cook_indices();
            let vanishes = cooked
                .as_view()
                .parse_to_symbolic_net::<Aind>(&ParseSettings::default())
                .is_ok_and(|network| has_odd_automorphism(&network));
            removed |= vanishes;
            !vanishes
        })
        .map(|term| term.to_owned())
        .collect::<Vec<_>>();

    if removed {
        retained.into_iter().sum()
    } else {
        expression.to_owned()
    }
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

fn requires_term_expansion<Aind: AbsInd>(network: &SymbolicNet<Aind>) -> bool {
    let tree = network.graph.expr_tree();
    network.graph.graph.iter_nodes().any(|(node, _, data)| {
        if !matches!(data, NetworkNode::Op(NetworkOp::Sum)) {
            return false;
        }

        let mut below_product_or_power = false;
        for ancestor in tree
            .ancestor_iter_node(node, network.graph.graph.as_ref())
            .skip(1)
        {
            match network.graph.graph[ancestor] {
                NetworkNode::Op(NetworkOp::Function(_)) => return false,
                NetworkNode::Op(NetworkOp::Product | NetworkOp::Power(1..)) => {
                    below_product_or_power = true;
                }
                _ => {}
            }
        }
        below_product_or_power
    })
}

fn has_odd_automorphism<Aind: AbsInd + ParseableAind>(network: &SymbolicNet<Aind>) -> bool {
    if !contains_antisymmetric_tensor(network) {
        return false;
    }

    let Some(graph) = project_network(network) else {
        return false;
    };
    // An odd stabilizer generator proves that the monomial equals its negative.
    let canonical = graph.canonize();
    canonical
        .orbit_generators
        .iter()
        .any(|generator| generator_is_odd(&canonical.graph, generator))
}

fn project_network<Aind: AbsInd + ParseableAind>(
    network: &SymbolicNet<Aind>,
) -> Option<TensorGraph<Aind>> {
    let tree: SimpleTraversalTree<ParentChildStore<()>> = network.graph.expr_tree().cast();
    let root = network.graph.graph.node_id(network.graph.head());
    let mut graph = Graph::new();
    let mut slot_copies = BTreeMap::new();
    project_expression(network, &tree, root, &mut graph, &mut slot_copies)?;
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
        // `expand` deliberately leaves function arguments opaque. An
        // automorphism confined to one summand there does not multiply the
        // enclosing function, but independent factors remain analyzable.
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
    tensor: &SymbolicTensor<Aind>,
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
    tensor: &SymbolicTensor<Aind>,
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
