use std::collections::{BTreeMap, BTreeSet};

use linnet::{
    half_edge::{
        NodeIndex,
        involution::{Hedge, HedgePair},
        tree::SimpleTraversalTree,
    },
    tree::child_pointer::ParentChildStore,
    union_find::UnionFind,
};
use spenso::{
    network::{
        graph::{NetworkEdge, NetworkLeaf, NetworkNode, NetworkOp},
        store::TensorScalarStore,
    },
    structure::{
        TensorStructure,
        representation::{LibraryRep, LibrarySlot, Representation},
        slot::{AbsInd, DummyAind, IsAbstractSlot, ParseableAind},
    },
};
use symbolica::{
    atom::{Atom, Symbol},
    graph::{Graph, HiddenData},
};

use super::{
    CanonicalPolicyNet, CanonicalizationError, GroupKey, IncidenceRole, SymmetryKind, TensorColor,
    TensorLayout, group::SignedGroup, semantic,
};
use crate::tensor::{SymbolicNet, SymbolicTensor};

/// The literal-copy projection is deliberately bounded: Graphica cost is a
/// property of the complete expanded problem, not just of the input network.
/// The selected limits are recorded in the architecture plan's Phase 6 table.
pub(super) const DEFAULT_GRAPH_BUDGET: GraphBudget = GraphBudget {
    vertices: 128,
    edges: 160,
};

#[cfg(test)]
thread_local! {
    static GRAPHICA_CALLS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

#[cfg(test)]
pub(super) fn reset_graphica_calls() {
    GRAPHICA_CALLS.with(|calls| calls.set(0));
}

#[cfg(test)]
pub(super) fn graphica_calls() -> usize {
    GRAPHICA_CALLS.with(std::cell::Cell::get)
}

#[derive(Clone, Copy, Debug)]
pub(super) struct GraphBudget {
    pub(super) vertices: usize,
    pub(super) edges: usize,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(super) enum UnifiedNodeColor {
    Root,
    Product,
    Sum,
    Neg,
    Function(semantic::SemanticSymbolKey),
    PowerResult(i8),
    PowerMagnitude,
    PowerCopy,
    Scalar(semantic::SemanticAtomKey),
    Tensor(TensorColor),
    Port(semantic::SemanticAtomKey),
    InternalLine(semantic::SemanticAtomKey),
    ExternalLine(semantic::SemanticAtomKey),
}

/// Hidden provenance retained for reconstruction and invariant diagnostics.
/// It deliberately does not participate in graph identity.
#[allow(dead_code)]
#[derive(Clone, Debug)]
pub(super) enum UnifiedNodeOrigin {
    Root,
    Expression(usize),
    Magnitude(usize),
    Copy { power: usize, copy: usize },
    Port { occurrence: usize, flat_slot: usize },
    Line(usize),
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(super) enum UnifiedEdgeColor {
    Root,
    Child,
    Argument,
    Operand,
    Magnitude,
    Copy,
    Base,
    Incidence(IncidenceRole),
    Cyclic(GroupKey),
    PortLine,
}

impl UnifiedEdgeColor {
    fn directed(&self) -> bool {
        !matches!(self, Self::PortLine)
    }
}

/// Hidden incidence provenance retained independently of visible edge color.
#[allow(dead_code)]
#[derive(Clone, Debug)]
pub(super) enum UnifiedEdgeOrigin {
    None,
    Incidence { flat_slot: usize, member: usize },
}

pub(super) type UnifiedGraph = Graph<
    HiddenData<UnifiedNodeColor, UnifiedNodeOrigin>,
    HiddenData<UnifiedEdgeColor, UnifiedEdgeOrigin>,
>;

#[derive(Clone, Debug)]
pub(super) enum ExpressionKind {
    Product,
    Sum,
    Neg,
    Function(Symbol),
    Power {
        exponent: i8,
        magnitude: usize,
        copies: Vec<usize>,
    },
    Scalar(Atom),
    Tensor(usize),
}

#[derive(Clone, Debug)]
pub(super) struct ExpressionOccurrence {
    pub(super) root: usize,
    pub(super) kind: ExpressionKind,
    pub(super) children: Vec<usize>,
}

#[derive(Clone, Debug)]
pub(super) struct TensorPort<Aind: AbsInd> {
    pub(super) flat_slot: usize,
    pub(super) member: usize,
    pub(super) vertex: usize,
    pub(super) line: usize,
    pub(super) slot: LibrarySlot<Aind>,
    pub(super) role: IncidenceRole,
}

#[derive(Clone, Debug)]
pub(super) struct TensorOccurrence<Aind: AbsInd> {
    pub(super) expression: usize,
    pub(super) header: usize,
    pub(super) layout: TensorLayout,
    pub(super) tensor: SymbolicTensor<Aind>,
    pub(super) ports: Vec<TensorPort<Aind>>,
}

#[derive(Clone, Debug)]
pub(super) struct LineOccurrence<Aind: AbsInd> {
    pub(super) vertex: usize,
    pub(super) group: Representation<LibraryRep>,
    pub(super) external: Option<LibrarySlot<Aind>>,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(super) struct ProblemIdentity {
    nodes: Vec<UnifiedNodeColor>,
    edges: Vec<(usize, usize, UnifiedEdgeColor)>,
}

pub(super) struct CanonicalProjection<Aind: AbsInd> {
    pub(super) graph: UnifiedGraph,
    pub(super) vertex_map: Vec<usize>,
    pub(super) orbit_generators: Vec<Vec<Vec<usize>>>,
    pub(super) expressions: Vec<ExpressionOccurrence>,
    pub(super) tensors: Vec<TensorOccurrence<Aind>>,
    pub(super) lines: Vec<LineOccurrence<Aind>>,
    pub(super) powers: Vec<PowerDescriptor>,
    pub(super) root_expression: usize,
    pub(super) identity: ProblemIdentity,
}

impl<Aind: AbsInd> CanonicalProjection<Aind> {
    /// Return the checked complete-copy seam for a materialized Power.
    pub(super) fn power_descriptor(&self, expression: usize) -> Option<&PowerDescriptor> {
        self.powers
            .iter()
            .find(|power| power.expression == expression)
    }
}

/// One canonical line at the exposed boundary of a complete Power copy.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(super) struct PowerBoundaryDescriptor {
    pub(super) line: usize,
    pub(super) ports: Vec<usize>,
    pub(super) target: PowerBoundaryTarget,
}

/// The other side of a Power-copy boundary line.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(super) enum PowerBoundaryTarget {
    Pair {
        partner_copy_root: usize,
        partner_ports: Vec<usize>,
    },
    Interface {
        outside_ports: Vec<usize>,
        externally_open: bool,
    },
}

/// Canonical membership and interface data for one complete cloned base.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct PowerCopyDescriptor {
    pub(super) expression: usize,
    pub(super) root: usize,
    pub(super) base_root: usize,
    pub(super) vertices: BTreeSet<usize>,
    pub(super) boundaries: Vec<PowerBoundaryDescriptor>,
}

/// Checked canonical metadata for one signed Power occurrence.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct PowerDescriptor {
    pub(super) expression: usize,
    pub(super) exponent: i8,
    pub(super) result: usize,
    pub(super) magnitude: usize,
    pub(super) copies: Vec<PowerCopyDescriptor>,
}

#[derive(Clone, Debug)]
struct Component<Aind: AbsInd> {
    group: Representation<LibraryRep>,
    external: Option<LibrarySlot<Aind>>,
    nodes: BTreeSet<NodeIndex>,
}

type LineComponents<Aind> = (Vec<Option<usize>>, Vec<Component<Aind>>);

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum LineStep {
    Copy { power: usize, copy: usize },
    Pair { power: usize, pair: usize },
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct LineKey {
    component: usize,
    steps: Vec<LineStep>,
}

#[derive(Clone)]
struct PowerFrame {
    power: usize,
    copy: usize,
    magnitude: usize,
    base_nodes: BTreeSet<NodeIndex>,
}

#[derive(Clone, Debug)]
struct PowerDraft {
    power: usize,
    expression: usize,
    exponent: i8,
    result: usize,
    magnitude: usize,
    copies: Vec<PowerCopyDraft>,
}

#[derive(Clone, Debug)]
struct PowerCopyDraft {
    copy: usize,
    expression: usize,
    root: usize,
    base_root: usize,
}

impl PowerBoundaryDescriptor {
    fn mapped<E>(&self, map_vertex: &impl Fn(usize) -> Result<usize, E>) -> Result<Self, E> {
        let mut ports = self
            .ports
            .iter()
            .map(|&port| map_vertex(port))
            .collect::<Result<Vec<_>, _>>()?;
        ports.sort_unstable();
        let target = match &self.target {
            PowerBoundaryTarget::Pair {
                partner_copy_root,
                partner_ports,
            } => {
                let mut partner_ports = partner_ports
                    .iter()
                    .map(|&port| map_vertex(port))
                    .collect::<Result<Vec<_>, _>>()?;
                partner_ports.sort_unstable();
                PowerBoundaryTarget::Pair {
                    partner_copy_root: map_vertex(*partner_copy_root)?,
                    partner_ports,
                }
            }
            PowerBoundaryTarget::Interface {
                outside_ports,
                externally_open,
            } => {
                let mut outside_ports = outside_ports
                    .iter()
                    .map(|&port| map_vertex(port))
                    .collect::<Result<Vec<_>, _>>()?;
                outside_ports.sort_unstable();
                PowerBoundaryTarget::Interface {
                    outside_ports,
                    externally_open: *externally_open,
                }
            }
        };
        Ok(Self {
            line: map_vertex(self.line)?,
            ports,
            target,
        })
    }

    fn mapped_infallible(&self, vertex_map: &[usize]) -> Self {
        self.mapped(&|vertex| Ok::<_, std::convert::Infallible>(vertex_map[vertex]))
            .unwrap()
    }
}

impl PowerDescriptor {
    fn validate_local(&self) -> Result<(), CanonicalizationError> {
        if matches!(self.exponent, 0 | 1)
            || self.copies.len() != usize::from(self.exponent.unsigned_abs())
        {
            return Err(CanonicalizationError::Projection(format!(
                "Power {} has copy count {} inconsistent with exponent {}",
                self.expression,
                self.copies.len(),
                self.exponent
            )));
        }
        let roots = self
            .copies
            .iter()
            .enumerate()
            .map(|(copy, descriptor)| (descriptor.root, copy))
            .collect::<BTreeMap<_, _>>();
        if roots.len() != self.copies.len() {
            return Err(CanonicalizationError::Projection(format!(
                "Power {} repeats a copy root",
                self.expression
            )));
        }
        let mut complete_vertices = BTreeSet::new();
        for copy in &self.copies {
            if !copy.vertices.contains(&copy.root) || !copy.vertices.contains(&copy.base_root) {
                return Err(CanonicalizationError::Projection(format!(
                    "Power {} copy root {} has incomplete rooted membership",
                    self.expression, copy.root
                )));
            }
            if copy
                .vertices
                .iter()
                .any(|&vertex| !complete_vertices.insert(vertex))
            {
                return Err(CanonicalizationError::Projection(format!(
                    "Power {} complete-copy vertex sets overlap",
                    self.expression
                )));
            }
            let mut lines = BTreeSet::new();
            let mut boundary_partner = None;
            let mut reaches_interface = false;
            for boundary in &copy.boundaries {
                if boundary.ports.is_empty()
                    || !lines.insert(boundary.line)
                    || copy.vertices.contains(&boundary.line)
                    || boundary
                        .ports
                        .iter()
                        .any(|port| !copy.vertices.contains(port))
                {
                    return Err(CanonicalizationError::Projection(format!(
                        "Power {} copy root {} has an invalid boundary descriptor",
                        self.expression, copy.root
                    )));
                }
                match &boundary.target {
                    PowerBoundaryTarget::Pair {
                        partner_copy_root,
                        partner_ports,
                    } => {
                        if reaches_interface
                            || boundary_partner
                                .replace(*partner_copy_root)
                                .is_some_and(|expected| expected != *partner_copy_root)
                        {
                            return Err(CanonicalizationError::Projection(format!(
                                "Power {} copy root {} has inconsistent partner copies across boundary slots",
                                self.expression, copy.root
                            )));
                        }
                        let Some(&partner) = roots.get(partner_copy_root) else {
                            return Err(CanonicalizationError::Projection(format!(
                                "Power {} copy root {} pairs outside its magnitude scope",
                                self.expression, copy.root
                            )));
                        };
                        if partner_copy_root == &copy.root
                            || partner_ports.is_empty()
                            || partner_ports
                                .iter()
                                .any(|port| !self.copies[partner].vertices.contains(port))
                        {
                            return Err(CanonicalizationError::Projection(format!(
                                "Power {} copy root {} has an incomplete partner boundary",
                                self.expression, copy.root
                            )));
                        }
                        let reciprocal = self.copies[partner].boundaries.iter().any(|candidate| {
                            candidate.line == boundary.line
                                && candidate.ports == *partner_ports
                                && matches!(
                                    &candidate.target,
                                    PowerBoundaryTarget::Pair {
                                        partner_copy_root,
                                        partner_ports,
                                    } if *partner_copy_root == copy.root
                                        && *partner_ports == boundary.ports
                                )
                        });
                        if !reciprocal {
                            return Err(CanonicalizationError::Projection(format!(
                                "Power {} boundary pairing at line {} is not reciprocal",
                                self.expression, boundary.line
                            )));
                        }
                    }
                    PowerBoundaryTarget::Interface { outside_ports, .. } => {
                        if boundary_partner.is_some()
                            || self.copies.len().is_multiple_of(2)
                            || copy.root != self.copies.last().unwrap().root
                            || outside_ports
                                .iter()
                                .any(|port| copy.vertices.contains(port))
                        {
                            return Err(CanonicalizationError::Projection(format!(
                                "Power {} has an invalid odd-copy result interface",
                                self.expression
                            )));
                        }
                        reaches_interface = true;
                    }
                }
            }
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct Estimate {
    vertices: usize,
    edges: usize,
    largest_power: Option<(i8, u8)>,
}

#[derive(Default)]
struct EstimateState {
    size: Estimate,
}

#[derive(Clone, Copy, Debug, Default)]
struct LineSummary {
    plain: bool,
    stepped: usize,
}

pub(super) fn project<Aind>(
    policy: &CanonicalPolicyNet<Aind>,
    budget: GraphBudget,
) -> Result<CanonicalProjection<Aind>, CanonicalizationError>
where
    Aind: AbsInd + DummyAind + ParseableAind,
{
    ProjectionBuilder::new(policy.network(), policy.normalized_atom(), budget)?.build()
}

struct ProjectionBuilder<'a, Aind: AbsInd> {
    network: &'a SymbolicNet<Aind>,
    source: &'a Atom,
    budget: GraphBudget,
    tree: SimpleTraversalTree<ParentChildStore<()>>,
    hedge_component: Vec<Option<usize>>,
    components: Vec<Component<Aind>>,
    graph: UnifiedGraph,
    expressions: Vec<ExpressionOccurrence>,
    tensors: Vec<TensorOccurrence<Aind>>,
    lines: Vec<LineOccurrence<Aind>>,
    line_by_key: BTreeMap<LineKey, usize>,
    powers: Vec<PowerDraft>,
    next_power: usize,
}

impl<'a, Aind> ProjectionBuilder<'a, Aind>
where
    Aind: AbsInd + ParseableAind,
{
    fn new(
        network: &'a SymbolicNet<Aind>,
        source: &'a Atom,
        budget: GraphBudget,
    ) -> Result<Self, CanonicalizationError> {
        let tree: SimpleTraversalTree<ParentChildStore<()>> = network.graph.expr_tree().cast();
        let (hedge_component, components) = Self::line_components(network, &tree)?;
        Ok(Self {
            network,
            source,
            budget,
            tree,
            hedge_component,
            components,
            graph: Graph::new(),
            expressions: Vec::new(),
            tensors: Vec::new(),
            lines: Vec::new(),
            line_by_key: BTreeMap::new(),
            powers: Vec::new(),
            next_power: 0,
        })
    }

    fn build(mut self) -> Result<CanonicalProjection<Aind>, CanonicalizationError> {
        let root_expression = self.project_visible_graph()?;

        #[cfg(test)]
        GRAPHICA_CALLS.with(|calls| calls.set(calls.get() + 1));
        let canonical = self.graph.canonize();
        let powers = self.canonical_power_descriptors(&canonical.vertex_map)?;
        Self::validate_power_automorphisms(
            canonical.graph.nodes().len(),
            &powers,
            &canonical.orbit_generators,
        )?;
        let identity = Self::problem_identity(&canonical.graph);
        Ok(CanonicalProjection {
            graph: canonical.graph,
            vertex_map: canonical.vertex_map,
            orbit_generators: canonical.orbit_generators,
            expressions: self.expressions,
            tensors: self.tensors,
            lines: self.lines,
            powers,
            root_expression,
            identity,
        })
    }

    fn project_visible_graph(&mut self) -> Result<usize, CanonicalizationError> {
        let root = self.network.graph.graph.node_id(self.network.graph.head());
        let estimate = self.estimate_visible_graph(root)?;

        if estimate.vertices > self.budget.vertices || estimate.edges > self.budget.edges {
            return Err(self.graph_size_error(
                estimate.largest_power,
                estimate.vertices,
                estimate.edges,
            ));
        }

        let root_marker = self.graph.add_node(HiddenData::new(
            UnifiedNodeColor::Root,
            UnifiedNodeOrigin::Root,
        ));
        let root_expression = self.project_node(root, &[])?;
        self.add_edge(
            root_marker,
            self.expressions[root_expression].root,
            UnifiedEdgeColor::Root,
            UnifiedEdgeOrigin::None,
        )?;

        if self.graph.nodes().len() > self.budget.vertices
            || self.graph.edges().len() > self.budget.edges
        {
            return Err(self.graph_size_error(
                estimate.largest_power,
                self.graph.nodes().len(),
                self.graph.edges().len(),
            ));
        }
        if self.graph.nodes().len() != estimate.vertices
            || self.graph.edges().len() != estimate.edges
        {
            return Err(CanonicalizationError::Projection(format!(
                "canonical graph size estimator predicted {}/{} vertices/edges but constructed {}/{}",
                estimate.vertices,
                estimate.edges,
                self.graph.nodes().len(),
                self.graph.edges().len()
            )));
        }

        Ok(root_expression)
    }

    fn line_components(
        network: &SymbolicNet<Aind>,
        tree: &SimpleTraversalTree<ParentChildStore<()>>,
    ) -> Result<LineComponents<Aind>, CanonicalizationError> {
        let graph = &network.graph.graph;
        let mut union = UnionFind::new(vec![(); graph.n_hedges()]);
        for (pair, _, data) in graph.iter_edges() {
            if matches!(data.data, NetworkEdge::Slot(_))
                && let HedgePair::Paired { source, sink } = pair
            {
                union.union(source, sink, |(), ()| ());
            }
        }

        let root = graph.node_id(network.graph.head());
        for node in tree.iter_preorder_tree_nodes(graph.as_ref(), root) {
            if !matches!(graph[node], NetworkNode::Op(_)) {
                continue;
            }
            let mut by_slot = BTreeMap::new();
            for hedge in graph.iter_crown(node) {
                let NetworkEdge::Slot(slot) = graph[[&hedge]] else {
                    continue;
                };
                if let Some(previous) = by_slot.insert(slot, hedge) {
                    union.union(previous, hedge, |(), ()| ());
                }
            }
        }

        let mut roots = BTreeMap::<Hedge, usize>::new();
        let mut hedge_component = vec![None; graph.n_hedges()];
        let mut components = Vec::<Component<Aind>>::new();
        for hedge in (0..graph.n_hedges()).map(Hedge) {
            let NetworkEdge::Slot(slot) = graph[[&hedge]] else {
                continue;
            };
            let representative = union.find(hedge);
            let component = *roots.entry(representative).or_insert_with(|| {
                let component = components.len();
                components.push(Component {
                    group: slot.rep().base(),
                    external: None,
                    nodes: BTreeSet::new(),
                });
                component
            });
            hedge_component[hedge.0] = Some(component);
            let entry = &mut components[component];
            if entry.group != slot.rep().base() {
                return Err(CanonicalizationError::Projection(format!(
                    "line component mixes representation groups {} and {}",
                    entry.group,
                    slot.rep().base()
                )));
            }
            entry.nodes.insert(graph.node_id(hedge));
            if graph.inv(hedge) == hedge {
                match entry.external {
                    Some(external) if external != slot => {
                        return Err(CanonicalizationError::Projection(format!(
                            "line component mixes external slots {external} and {slot}"
                        )));
                    }
                    None => entry.external = Some(slot),
                    Some(_) => {}
                }
            }
        }

        for (node, _, data) in graph.iter_nodes() {
            let NetworkNode::Leaf(NetworkLeaf::LocalTensor(tensor_index)) = data else {
                continue;
            };
            let tensor = &network.store.tensors[*tensor_index];
            let structural_slots = tensor
                .structure
                .external_structure_iter()
                .collect::<Vec<_>>();
            let mut hedges = graph
                .iter_crown(node)
                .filter(|hedge| matches!(graph[[hedge]], NetworkEdge::Slot(_)))
                .collect::<Vec<_>>();
            hedges.sort_unstable_by_key(|hedge| network.graph.slot_order[hedge.0]);
            if structural_slots.len() != hedges.len() {
                return Err(CanonicalizationError::StructureMismatch {
                    expression: tensor.expression.clone(),
                });
            }
            for (slot, hedge) in structural_slots.into_iter().zip(hedges) {
                let component = hedge_component[hedge.0].ok_or_else(|| {
                    CanonicalizationError::Projection(format!(
                        "tensor slot hedge {} is absent from line closure",
                        hedge.0
                    ))
                })?;
                let edge_slot = match graph[[&hedge]] {
                    NetworkEdge::Slot(edge_slot) => edge_slot,
                    NetworkEdge::Head => unreachable!("slot hedges were filtered above"),
                };
                if slot.aind() != edge_slot.aind()
                    || slot.rep().base() != components[component].group
                {
                    return Err(CanonicalizationError::Projection(format!(
                        "tensor endpoint {slot} conflicts with line-component declaration {edge_slot}"
                    )));
                }
            }
        }
        Ok((hedge_component, components))
    }

    fn estimate_visible_graph(&self, root: NodeIndex) -> Result<Estimate, CanonicalizationError> {
        let mut state = EstimateState::default();
        // The unique root marker and its directed edge are present around every
        // projected expression, including normalized nullary operations.
        self.add_estimated_size(&mut state, 1, 1)?;
        self.estimate_node(root, &mut state)?;
        for component in 0..self.components.len() {
            let lines = self.estimate_component_lines(root, component, state.size.largest_power)?;
            let vertices = lines
                .stepped
                .checked_add(usize::from(lines.plain))
                .ok_or_else(|| {
                    self.graph_size_error(state.size.largest_power, usize::MAX, usize::MAX)
                })?;
            self.add_estimated_size(&mut state, vertices, 0)?;
        }
        Ok(state.size)
    }

    fn estimate_node(
        &self,
        node: NodeIndex,
        state: &mut EstimateState,
    ) -> Result<(), CanonicalizationError> {
        match &self.network.graph.graph[node] {
            NetworkNode::Leaf(NetworkLeaf::Scalar(_)) => self.add_estimated_size(state, 1, 0),
            NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => {
                self.estimate_tensor(node, *index, state)
            }
            NetworkNode::Leaf(_) => Err(CanonicalizationError::UnsupportedLeaf),
            NetworkNode::Op(NetworkOp::Product) => self.estimate_associative(node, true, state),
            NetworkNode::Op(NetworkOp::Sum) => self.estimate_associative(node, false, state),
            NetworkNode::Op(NetworkOp::Neg | NetworkOp::Function(_)) => {
                let child = self.only_child(node)?;
                self.add_estimated_size(state, 1, 1)?;
                self.estimate_node(child, state)
            }
            NetworkNode::Op(NetworkOp::Power(0)) => self.add_estimated_size(state, 1, 0),
            NetworkNode::Op(NetworkOp::Power(1)) => {
                self.estimate_node(self.only_child(node)?, state)
            }
            NetworkNode::Op(NetworkOp::Power(exponent)) => {
                let magnitude = exponent.unsigned_abs();
                if state
                    .size
                    .largest_power
                    .is_none_or(|(_, existing)| magnitude > existing)
                {
                    state.size.largest_power = Some((*exponent, magnitude));
                }
                let copies = usize::from(magnitude);
                let copy_edges = copies.checked_mul(2).ok_or_else(|| {
                    self.graph_size_error(state.size.largest_power, usize::MAX, usize::MAX)
                })?;
                // Result and magnitude roots, every copy root, and their
                // magnitude/copy, copy/base, and result/magnitude edges.
                self.add_estimated_size(state, 2 + copies, 1 + copy_edges)?;
                let child = self.only_child(node)?;
                let before = state.size;
                self.estimate_node(child, state)?;
                let child_vertices = state.size.vertices - before.vertices;
                let child_edges = state.size.edges - before.edges;
                if copies > 1 {
                    let repeated = copies - 1;
                    let vertices = child_vertices.checked_mul(repeated).ok_or_else(|| {
                        self.graph_size_error(state.size.largest_power, usize::MAX, usize::MAX)
                    })?;
                    let edges = child_edges.checked_mul(repeated).ok_or_else(|| {
                        self.graph_size_error(state.size.largest_power, usize::MAX, usize::MAX)
                    })?;
                    self.add_estimated_size(state, vertices, edges)?;
                }
                Ok(())
            }
        }
    }

    fn estimate_associative(
        &self,
        source: NodeIndex,
        product: bool,
        state: &mut EstimateState,
    ) -> Result<(), CanonicalizationError> {
        let mut sources = Vec::new();
        self.flatten_associative(source, product, &mut sources);
        match sources.as_slice() {
            [] => self.add_estimated_size(state, 1, 0),
            [source] => self.estimate_node(*source, state),
            _ => {
                self.add_estimated_size(state, 1, sources.len())?;
                for source in sources {
                    self.estimate_node(source, state)?;
                }
                Ok(())
            }
        }
    }

    fn estimate_tensor(
        &self,
        source: NodeIndex,
        tensor_index: usize,
        state: &mut EstimateState,
    ) -> Result<(), CanonicalizationError> {
        let tensor = &self.network.store.tensors[tensor_index];
        let layout = TensorLayout::scan(tensor)?;
        let mut hedges = self
            .network
            .graph
            .graph
            .iter_crown(source)
            .filter(|hedge| matches!(self.network.graph.graph[[hedge]], NetworkEdge::Slot(_)))
            .collect::<Vec<_>>();
        hedges.sort_unstable_by_key(|hedge| self.network.graph.slot_order[hedge.0]);
        if tensor.structure.external_structure_iter().count() != layout.slot_count
            || hedges.len() != layout.slot_count
        {
            return Err(CanonicalizationError::StructureMismatch {
                expression: tensor.expression.clone(),
            });
        }

        // Tensor header, followed by one port and two incidence edges per slot.
        self.add_estimated_size(state, 1, 0)?;
        let mut cyclic = BTreeMap::<GroupKey, usize>::new();
        for (structural_position, hedge) in hedges.into_iter().enumerate() {
            let flat_slot = layout.structural_holes[structural_position];
            self.hedge_component[hedge.0].ok_or_else(|| {
                CanonicalizationError::Projection(format!(
                    "tensor slot hedge {} is absent from line closure",
                    hedge.0
                ))
            })?;
            self.add_estimated_size(state, 1, 2)?;
            if let IncidenceRole::Group {
                key,
                kind: SymmetryKind::Cyclic,
            } = layout.incidence_role(flat_slot)
            {
                let members = cyclic.entry(key).or_default();
                *members = members.checked_add(1).ok_or_else(|| {
                    self.graph_size_error(state.size.largest_power, usize::MAX, usize::MAX)
                })?;
            }
        }
        for members in cyclic.into_values().filter(|members| *members > 1) {
            self.add_estimated_size(state, 0, members)?;
        }
        Ok(())
    }

    /// Count the distinct line keys for one original component without
    /// enumerating literal Power copies. `plain` records the key with no Power
    /// step; all other keys are counted by `stepped`.
    fn estimate_component_lines(
        &self,
        node: NodeIndex,
        component: usize,
        largest_power: Option<(i8, u8)>,
    ) -> Result<LineSummary, CanonicalizationError> {
        match &self.network.graph.graph[node] {
            NetworkNode::Leaf(NetworkLeaf::Scalar(_)) => Ok(LineSummary::default()),
            NetworkNode::Leaf(NetworkLeaf::LocalTensor(_)) => Ok(LineSummary {
                plain: self.components[component].nodes.contains(&node),
                stepped: 0,
            }),
            NetworkNode::Leaf(_) => Err(CanonicalizationError::UnsupportedLeaf),
            NetworkNode::Op(NetworkOp::Power(0)) => Ok(LineSummary::default()),
            NetworkNode::Op(NetworkOp::Power(1)) => {
                self.estimate_component_lines(self.only_child(node)?, component, largest_power)
            }
            NetworkNode::Op(NetworkOp::Power(exponent)) => {
                let child = self.only_child(node)?;
                let child_lines = self.estimate_component_lines(child, component, largest_power)?;
                let magnitude = usize::from(exponent.unsigned_abs());
                let mut stepped = child_lines
                    .stepped
                    .checked_mul(magnitude)
                    .ok_or_else(|| self.graph_size_error(largest_power, usize::MAX, usize::MAX))?;
                let base_nodes = self
                    .tree
                    .iter_preorder_tree_nodes(self.network.graph.graph.as_ref(), child)
                    .collect::<BTreeSet<_>>();
                let metadata = &self.components[component];
                let exposed = metadata.external.is_some()
                    || metadata.nodes.iter().any(|node| !base_nodes.contains(node));
                if child_lines.plain {
                    let added = if exposed { magnitude / 2 } else { magnitude };
                    stepped = stepped.checked_add(added).ok_or_else(|| {
                        self.graph_size_error(largest_power, usize::MAX, usize::MAX)
                    })?;
                }
                Ok(LineSummary {
                    plain: child_lines.plain && exposed && !magnitude.is_multiple_of(2),
                    stepped,
                })
            }
            NetworkNode::Op(_) => {
                let mut combined = LineSummary::default();
                for child in self.children(node) {
                    let child_lines =
                        self.estimate_component_lines(child, component, largest_power)?;
                    combined.plain |= child_lines.plain;
                    combined.stepped = combined
                        .stepped
                        .checked_add(child_lines.stepped)
                        .ok_or_else(|| {
                            self.graph_size_error(largest_power, usize::MAX, usize::MAX)
                        })?;
                }
                Ok(combined)
            }
        }
    }

    fn add_estimated_size(
        &self,
        state: &mut EstimateState,
        vertices: usize,
        edges: usize,
    ) -> Result<(), CanonicalizationError> {
        let requested_vertices = state.size.vertices.checked_add(vertices);
        let requested_edges = state.size.edges.checked_add(edges);
        let (Some(requested_vertices), Some(requested_edges)) =
            (requested_vertices, requested_edges)
        else {
            return Err(self.graph_size_error(
                state.size.largest_power,
                requested_vertices.unwrap_or(usize::MAX),
                requested_edges.unwrap_or(usize::MAX),
            ));
        };
        state.size.vertices = requested_vertices;
        state.size.edges = requested_edges;
        Ok(())
    }

    fn project_node(
        &mut self,
        source: NodeIndex,
        power_path: &[PowerFrame],
    ) -> Result<usize, CanonicalizationError> {
        match &self.network.graph.graph[source] {
            NetworkNode::Leaf(NetworkLeaf::Scalar(reference)) => {
                let atom = self.network.store.get_scalar_ref(*reference).clone();
                let expression = self.expressions.len();
                let root = self.graph.add_node(HiddenData::new(
                    UnifiedNodeColor::Scalar(semantic::SemanticAtomKey::new(atom.as_view())),
                    UnifiedNodeOrigin::Expression(expression),
                ));
                self.expressions.push(ExpressionOccurrence {
                    root,
                    kind: ExpressionKind::Scalar(atom),
                    children: Vec::new(),
                });
                Ok(expression)
            }
            NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => {
                self.project_tensor(source, *index, power_path)
            }
            NetworkNode::Leaf(_) => Err(CanonicalizationError::UnsupportedLeaf),
            NetworkNode::Op(NetworkOp::Product) => {
                self.project_associative(source, true, power_path)
            }
            NetworkNode::Op(NetworkOp::Sum) => self.project_associative(source, false, power_path),
            NetworkNode::Op(NetworkOp::Neg) => {
                let expression = self.add_operation(UnifiedNodeColor::Neg, ExpressionKind::Neg);
                let child = self.project_node(self.only_child(source)?, power_path)?;
                self.expressions[expression].children.push(child);
                self.add_edge(
                    self.expressions[expression].root,
                    self.expressions[child].root,
                    UnifiedEdgeColor::Operand,
                    UnifiedEdgeOrigin::None,
                )?;
                Ok(expression)
            }
            NetworkNode::Op(NetworkOp::Function(function)) => {
                let expression = self.add_operation(
                    UnifiedNodeColor::Function((*function).into()),
                    ExpressionKind::Function(*function),
                );
                let child = self.project_node(self.only_child(source)?, power_path)?;
                self.expressions[expression].children.push(child);
                self.add_edge(
                    self.expressions[expression].root,
                    self.expressions[child].root,
                    UnifiedEdgeColor::Argument,
                    UnifiedEdgeOrigin::None,
                )?;
                Ok(expression)
            }
            NetworkNode::Op(NetworkOp::Power(0)) => {
                let expression = self.expressions.len();
                let root = self.graph.add_node(HiddenData::new(
                    UnifiedNodeColor::Scalar(semantic::SemanticAtomKey::new(
                        Atom::num(1).as_view(),
                    )),
                    UnifiedNodeOrigin::Expression(expression),
                ));
                self.expressions.push(ExpressionOccurrence {
                    root,
                    kind: ExpressionKind::Scalar(Atom::num(1)),
                    children: Vec::new(),
                });
                Ok(expression)
            }
            NetworkNode::Op(NetworkOp::Power(1)) => {
                self.project_node(self.only_child(source)?, power_path)
            }
            NetworkNode::Op(NetworkOp::Power(exponent)) => {
                self.project_power(source, *exponent, power_path)
            }
        }
    }

    fn project_associative(
        &mut self,
        source: NodeIndex,
        product: bool,
        power_path: &[PowerFrame],
    ) -> Result<usize, CanonicalizationError> {
        let mut sources = Vec::new();
        self.flatten_associative(source, product, &mut sources);
        if sources.is_empty() {
            let atom = if product { Atom::num(1) } else { Atom::Zero };
            let expression = self.expressions.len();
            let root = self.graph.add_node(HiddenData::new(
                UnifiedNodeColor::Scalar(semantic::SemanticAtomKey::new(atom.as_view())),
                UnifiedNodeOrigin::Expression(expression),
            ));
            self.expressions.push(ExpressionOccurrence {
                root,
                kind: ExpressionKind::Scalar(atom),
                children: Vec::new(),
            });
            return Ok(expression);
        }
        if sources.len() == 1 {
            return self.project_node(sources[0], power_path);
        }

        let expression = self.add_operation(
            if product {
                UnifiedNodeColor::Product
            } else {
                UnifiedNodeColor::Sum
            },
            if product {
                ExpressionKind::Product
            } else {
                ExpressionKind::Sum
            },
        );
        for source in sources {
            let child = self.project_node(source, power_path)?;
            self.expressions[expression].children.push(child);
            self.add_edge(
                self.expressions[expression].root,
                self.expressions[child].root,
                UnifiedEdgeColor::Child,
                UnifiedEdgeOrigin::None,
            )?;
        }
        Ok(expression)
    }

    fn project_power(
        &mut self,
        source: NodeIndex,
        exponent: i8,
        power_path: &[PowerFrame],
    ) -> Result<usize, CanonicalizationError> {
        let child_source = self.only_child(source)?;
        let magnitude = usize::from(exponent.unsigned_abs());
        let power = self.next_power;
        self.next_power += 1;

        let expression = self.add_operation(
            UnifiedNodeColor::PowerResult(exponent),
            ExpressionKind::Power {
                exponent,
                magnitude: 0,
                copies: Vec::new(),
            },
        );
        let magnitude_root = self.graph.add_node(HiddenData::new(
            UnifiedNodeColor::PowerMagnitude,
            UnifiedNodeOrigin::Magnitude(expression),
        ));
        self.add_edge(
            self.expressions[expression].root,
            magnitude_root,
            UnifiedEdgeColor::Magnitude,
            UnifiedEdgeOrigin::None,
        )?;

        let base_nodes = self
            .tree
            .iter_preorder_tree_nodes(self.network.graph.graph.as_ref(), child_source)
            .collect::<BTreeSet<_>>();
        let mut copies = Vec::with_capacity(magnitude);
        let mut copy_drafts = Vec::with_capacity(magnitude);
        for copy in 0..magnitude {
            let copy_root = self.graph.add_node(HiddenData::new(
                UnifiedNodeColor::PowerCopy,
                UnifiedNodeOrigin::Copy { power, copy },
            ));
            self.add_edge(
                magnitude_root,
                copy_root,
                UnifiedEdgeColor::Copy,
                UnifiedEdgeOrigin::None,
            )?;
            let mut nested_path = power_path.to_vec();
            nested_path.push(PowerFrame {
                power,
                copy,
                magnitude,
                base_nodes: base_nodes.clone(),
            });
            let base = self.project_node(child_source, &nested_path)?;
            self.add_edge(
                copy_root,
                self.expressions[base].root,
                UnifiedEdgeColor::Base,
                UnifiedEdgeOrigin::None,
            )?;
            copy_drafts.push(PowerCopyDraft {
                copy,
                expression: base,
                root: copy_root,
                base_root: self.expressions[base].root,
            });
            copies.push(base);
        }
        self.expressions[expression].kind = ExpressionKind::Power {
            exponent,
            magnitude: magnitude_root,
            copies,
        };
        self.powers.push(PowerDraft {
            power,
            expression,
            exponent,
            result: self.expressions[expression].root,
            magnitude: magnitude_root,
            copies: copy_drafts,
        });
        Ok(expression)
    }

    fn project_tensor(
        &mut self,
        source: NodeIndex,
        tensor_index: usize,
        power_path: &[PowerFrame],
    ) -> Result<usize, CanonicalizationError> {
        let tensor = self.network.store.tensors[tensor_index].clone();
        let layout = TensorLayout::scan(&tensor)?;
        let structural_slots = tensor
            .structure
            .external_structure_iter()
            .collect::<Vec<_>>();
        let mut hedges = self
            .network
            .graph
            .graph
            .iter_crown(source)
            .filter(|hedge| matches!(self.network.graph.graph[[hedge]], NetworkEdge::Slot(_)))
            .collect::<Vec<_>>();
        hedges.sort_unstable_by_key(|hedge| self.network.graph.slot_order[hedge.0]);
        if structural_slots.len() != layout.slot_count || hedges.len() != layout.slot_count {
            return Err(CanonicalizationError::StructureMismatch {
                expression: tensor.expression,
            });
        }

        let expression = self.expressions.len();
        let occurrence = self.tensors.len();
        let header = self.graph.add_node(HiddenData::new(
            UnifiedNodeColor::Tensor(layout.color.clone()),
            UnifiedNodeOrigin::Expression(expression),
        ));
        let mut ports = Vec::with_capacity(layout.slot_count);
        let mut cyclic = BTreeMap::<GroupKey, Vec<(usize, usize)>>::new();
        for (structural_position, hedge) in hedges.into_iter().enumerate() {
            let flat_slot = layout.structural_holes[structural_position];
            // A paired network edge stores one shared slot payload; the tensor
            // structure retains the orientation of this particular endpoint.
            let slot = structural_slots[structural_position];
            let component = self.hedge_component[hedge.0].ok_or_else(|| {
                CanonicalizationError::Projection(format!(
                    "tensor slot hedge {} is absent from line closure",
                    hedge.0
                ))
            })?;
            let line = self.line(component, power_path);
            let role = layout.incidence_role(flat_slot);
            let member = layout.member_position(flat_slot);
            let port = self.graph.add_node(HiddenData::new(
                UnifiedNodeColor::Port(semantic::representation_key(slot.rep())),
                UnifiedNodeOrigin::Port {
                    occurrence,
                    flat_slot,
                },
            ));
            self.add_edge(
                header,
                port,
                UnifiedEdgeColor::Incidence(role.clone()),
                UnifiedEdgeOrigin::Incidence { flat_slot, member },
            )?;
            self.add_edge(
                port,
                self.lines[line].vertex,
                UnifiedEdgeColor::PortLine,
                UnifiedEdgeOrigin::None,
            )?;
            if let IncidenceRole::Group {
                key,
                kind: SymmetryKind::Cyclic,
            } = role
            {
                cyclic.entry(key).or_default().push((member, port));
            }
            ports.push(TensorPort {
                flat_slot,
                member,
                vertex: port,
                line,
                slot,
                role,
            });
        }

        for (key, mut members) in cyclic {
            members.sort_unstable_by_key(|(member, _)| *member);
            if members.len() > 1 {
                for (&(_, source), &(_, target)) in members
                    .iter()
                    .zip(members.iter().cycle().skip(1))
                    .take(members.len())
                {
                    self.add_edge(
                        source,
                        target,
                        UnifiedEdgeColor::Cyclic(key),
                        UnifiedEdgeOrigin::None,
                    )?;
                }
            }
        }

        self.tensors.push(TensorOccurrence {
            expression,
            header,
            layout,
            tensor,
            ports,
        });
        self.expressions.push(ExpressionOccurrence {
            root: header,
            kind: ExpressionKind::Tensor(occurrence),
            children: Vec::new(),
        });
        Ok(expression)
    }

    fn line(&mut self, component: usize, power_path: &[PowerFrame]) -> usize {
        let key = self.line_key(component, power_path);
        if let Some(line) = self.line_by_key.get(&key) {
            return *line;
        }
        let metadata = &self.components[component];
        let external = key.steps.is_empty().then_some(metadata.external).flatten();
        let color = external.map_or_else(
            || UnifiedNodeColor::InternalLine(semantic::representation_key(metadata.group)),
            |slot| UnifiedNodeColor::ExternalLine(semantic::slot_key(slot)),
        );
        let line = self.lines.len();
        let vertex = self
            .graph
            .add_node(HiddenData::new(color, UnifiedNodeOrigin::Line(line)));
        self.lines.push(LineOccurrence {
            vertex,
            group: metadata.group,
            external,
        });
        self.line_by_key.insert(key, line);
        line
    }

    fn line_key(&self, component: usize, power_path: &[PowerFrame]) -> LineKey {
        let metadata = &self.components[component];
        let mut steps = Vec::new();
        for frame in power_path {
            let exposed = metadata.external.is_some()
                || metadata
                    .nodes
                    .iter()
                    .any(|node| !frame.base_nodes.contains(node));
            if exposed {
                if frame.magnitude % 2 == 0 || frame.copy + 1 < frame.magnitude {
                    steps.push(LineStep::Pair {
                        power: frame.power,
                        pair: frame.copy / 2,
                    });
                }
            } else {
                steps.push(LineStep::Copy {
                    power: frame.power,
                    copy: frame.copy,
                });
            }
        }
        LineKey { component, steps }
    }

    fn add_operation(&mut self, color: UnifiedNodeColor, kind: ExpressionKind) -> usize {
        let expression = self.expressions.len();
        let root = self.graph.add_node(HiddenData::new(
            color,
            UnifiedNodeOrigin::Expression(expression),
        ));
        self.expressions.push(ExpressionOccurrence {
            root,
            kind,
            children: Vec::new(),
        });
        expression
    }

    fn add_edge(
        &mut self,
        source: usize,
        target: usize,
        color: UnifiedEdgeColor,
        origin: UnifiedEdgeOrigin,
    ) -> Result<(), CanonicalizationError> {
        self.graph
            .add_edge(
                source,
                target,
                color.directed(),
                HiddenData::new(color, origin),
            )
            .map_err(|error| CanonicalizationError::Projection(error.to_string()))?;
        Ok(())
    }

    fn flatten_associative(&self, node: NodeIndex, product: bool, children: &mut Vec<NodeIndex>) {
        let same = matches!(
            (&self.network.graph.graph[node], product),
            (NetworkNode::Op(NetworkOp::Product), true) | (NetworkNode::Op(NetworkOp::Sum), false)
        );
        if same {
            for child in self.children(node) {
                self.flatten_associative(child, product, children);
            }
        } else {
            children.push(node);
        }
    }

    fn children(&self, node: NodeIndex) -> Vec<NodeIndex> {
        self.tree
            .iter_children(node, self.network.graph.graph.as_ref())
            .collect()
    }

    fn only_child(&self, node: NodeIndex) -> Result<NodeIndex, CanonicalizationError> {
        match self.children(node).as_slice() {
            [child] => Ok(*child),
            _ => Err(CanonicalizationError::Projection(format!(
                "unary operation at node {node:?} does not have exactly one child"
            ))),
        }
    }

    fn graph_size_error(
        &self,
        largest_power: Option<(i8, u8)>,
        requested_vertices: usize,
        requested_edges: usize,
    ) -> CanonicalizationError {
        match largest_power {
            Some((power, magnitude)) => CanonicalizationError::GraphSizeLimit {
                expression: self.source.clone(),
                power,
                magnitude,
                requested_vertices,
                requested_edges,
                vertex_limit: self.budget.vertices,
                edge_limit: self.budget.edges,
            },
            None => CanonicalizationError::WholeGraphSizeLimit {
                expression: self.source.clone(),
                requested_vertices,
                requested_edges,
                vertex_limit: self.budget.vertices,
                edge_limit: self.budget.edges,
            },
        }
    }

    fn canonical_power_descriptors(
        &self,
        vertex_map: &[usize],
    ) -> Result<Vec<PowerDescriptor>, CanonicalizationError> {
        if vertex_map.len() != self.graph.nodes().len() {
            return Err(CanonicalizationError::Projection(
                "canonical vertex map has the wrong domain for Power validation".into(),
            ));
        }

        let powers_by_expression = self
            .powers
            .iter()
            .map(|power| (power.expression, power))
            .collect::<BTreeMap<_, _>>();
        let mut line_keys = vec![None; self.lines.len()];
        for (key, &line) in &self.line_by_key {
            line_keys[line] = Some(key);
        }
        let mut ports_by_line = vec![Vec::<(usize, usize)>::new(); self.lines.len()];
        for tensor in &self.tensors {
            for port in &tensor.ports {
                ports_by_line[port.line].push((tensor.expression, port.vertex));
            }
        }
        let tensors_by_expression = self
            .tensors
            .iter()
            .enumerate()
            .map(|(tensor, occurrence)| (occurrence.expression, tensor))
            .collect::<BTreeMap<_, _>>();
        let mut internal_lines = BTreeMap::<(usize, usize), Vec<usize>>::new();
        for (key, &line) in &self.line_by_key {
            for step in &key.steps {
                if let LineStep::Copy { power, copy } = step {
                    internal_lines
                        .entry((*power, *copy))
                        .or_default()
                        .push(line);
                }
            }
        }

        self.powers
            .iter()
            .map(|power| {
                let mut expressions_by_copy = Vec::with_capacity(power.copies.len());
                let mut copy_by_expression = BTreeMap::new();
                for copy in &power.copies {
                    let expressions = self.expression_descendants(copy.expression)?;
                    for &expression in &expressions {
                        if copy_by_expression.insert(expression, copy.copy).is_some() {
                            return Err(CanonicalizationError::Projection(format!(
                                "Power {} has overlapping complete-copy expression membership",
                                power.expression
                            )));
                        }
                    }
                    expressions_by_copy.push(expressions);
                }

                let mut copies = Vec::with_capacity(power.copies.len());
                for copy in &power.copies {
                    let expressions = &expressions_by_copy[copy.copy];
                    let mut vertices = BTreeSet::from([copy.root]);
                    for &expression in expressions {
                        let occurrence = &self.expressions[expression];
                        vertices.insert(occurrence.root);
                        if let Some(&tensor) = tensors_by_expression.get(&expression) {
                            vertices.extend(
                                self.tensors[tensor].ports.iter().map(|port| port.vertex),
                            );
                        }
                        if let Some(nested) = powers_by_expression.get(&expression) {
                            vertices.insert(nested.magnitude);
                            vertices.extend(nested.copies.iter().map(|copy| copy.root));
                        }
                    }
                    if let Some(lines) = internal_lines.get(&(power.power, copy.copy)) {
                        for &line in lines {
                            vertices.insert(self.lines[line].vertex);
                        }
                    }

                    let mut ports_by_boundary = BTreeMap::<usize, Vec<usize>>::new();
                    for &expression in expressions {
                        if let Some(&tensor) = tensors_by_expression.get(&expression) {
                            for port in &self.tensors[tensor].ports {
                                if !vertices.contains(&self.lines[port.line].vertex) {
                                    ports_by_boundary
                                        .entry(port.line)
                                        .or_default()
                                        .push(port.vertex);
                                }
                            }
                        }
                    }
                    let mut boundaries = Vec::with_capacity(ports_by_boundary.len());
                    for (line, mut ports) in ports_by_boundary {
                        ports.sort_unstable();
                        ports.dedup();
                        let key = line_keys[line].ok_or_else(|| {
                            CanonicalizationError::Projection(format!(
                                "Power {} boundary line {line} has no line-closure key",
                                power.expression
                            ))
                        })?;
                        let step = key.steps.iter().find(|step| match step {
                            LineStep::Copy { power: candidate, .. }
                            | LineStep::Pair { power: candidate, .. } => {
                                *candidate == power.power
                            }
                        });
                        let target = match step {
                            Some(LineStep::Pair { pair, .. }) => {
                                if *pair != copy.copy / 2 {
                                    return Err(CanonicalizationError::Projection(format!(
                                        "Power {} copy {} has inconsistent boundary pair {pair}",
                                        power.expression, copy.copy
                                    )));
                                }
                                let expected_partner = copy.copy ^ 1;
                                let mut partner_ports = Vec::new();
                                let mut partner = None;
                                for &(expression, port) in &ports_by_line[line] {
                                    match copy_by_expression.get(&expression).copied() {
                                        Some(owner) if owner == copy.copy => {}
                                        Some(owner) => {
                                            if owner != expected_partner
                                                || partner.replace(owner).is_some_and(|old| old != owner)
                                            {
                                                return Err(CanonicalizationError::Projection(
                                                    format!(
                                                        "Power {} copy {} boundary line {line} pairs with multiple or non-partner copies",
                                                        power.expression, copy.copy
                                                    ),
                                                ));
                                            }
                                            partner_ports.push(port);
                                        }
                                        None => {
                                            return Err(CanonicalizationError::Projection(format!(
                                                "Power {} paired boundary line {line} escapes its magnitude scope",
                                                power.expression
                                            )));
                                        }
                                    }
                                }
                                if partner != Some(expected_partner) || partner_ports.is_empty() {
                                    return Err(CanonicalizationError::Projection(format!(
                                        "Power {} copy {} has an incomplete boundary pairing on line {line}",
                                        power.expression, copy.copy
                                    )));
                                }
                                partner_ports.sort_unstable();
                                partner_ports.dedup();
                                PowerBoundaryTarget::Pair {
                                    partner_copy_root: power.copies[expected_partner].root,
                                    partner_ports,
                                }
                            }
                            Some(LineStep::Copy { .. }) => {
                                return Err(CanonicalizationError::Projection(format!(
                                    "Power {} internal copy line {line} was classified as a boundary",
                                    power.expression
                                )));
                            }
                            None => {
                                if power.copies.len() % 2 == 0
                                    || copy.copy + 1 != power.copies.len()
                                {
                                    return Err(CanonicalizationError::Projection(format!(
                                        "Power {} copy {} reaches the result interface outside the unique odd remainder",
                                        power.expression, copy.copy
                                    )));
                                }
                                if ports_by_line[line].iter().any(|(expression, _)| {
                                    copy_by_expression
                                        .get(expression)
                                        .is_some_and(|owner| *owner != copy.copy)
                                }) {
                                    return Err(CanonicalizationError::Projection(format!(
                                        "Power {} odd result interface line {line} reaches another materialized copy",
                                        power.expression
                                    )));
                                }
                                let mut outside_ports = ports_by_line[line]
                                    .iter()
                                    .filter_map(|&(expression, port)| {
                                        (!expressions.contains(&expression)).then_some(port)
                                    })
                                    .collect::<Vec<_>>();
                                outside_ports.sort_unstable();
                                outside_ports.dedup();
                                PowerBoundaryTarget::Interface {
                                    outside_ports,
                                    externally_open: self.lines[line].external.is_some(),
                                }
                            }
                        };
                        boundaries.push(PowerBoundaryDescriptor {
                            line: self.lines[line].vertex,
                            ports,
                            target,
                        });
                    }

                    let map_vertex = |vertex: usize| {
                        vertex_map.get(vertex).copied().ok_or_else(|| {
                            CanonicalizationError::Projection(format!(
                                "Power {} references vertex {vertex} outside the canonical map",
                                power.expression
                            ))
                        })
                    };
                    let mut canonical_boundaries = boundaries
                        .iter()
                        .map(|boundary| boundary.mapped(&map_vertex))
                        .collect::<Result<Vec<_>, _>>()?;
                    canonical_boundaries.sort();
                    copies.push(PowerCopyDescriptor {
                        expression: copy.expression,
                        root: map_vertex(copy.root)?,
                        base_root: map_vertex(copy.base_root)?,
                        vertices: vertices
                            .into_iter()
                            .map(map_vertex)
                            .collect::<Result<_, _>>()?,
                        boundaries: canonical_boundaries,
                    });
                }

                let descriptor = PowerDescriptor {
                    expression: power.expression,
                    exponent: power.exponent,
                    result: vertex_map[power.result],
                    magnitude: vertex_map[power.magnitude],
                    copies,
                };
                descriptor.validate_local()?;
                Ok(descriptor)
            })
            .collect()
    }

    fn expression_descendants(
        &self,
        root: usize,
    ) -> Result<BTreeSet<usize>, CanonicalizationError> {
        let mut descendants = BTreeSet::new();
        let mut pending = vec![root];
        while let Some(expression) = pending.pop() {
            if !descendants.insert(expression) {
                continue;
            }
            let occurrence = self.expressions.get(expression).ok_or_else(|| {
                CanonicalizationError::Projection(format!(
                    "Power copy references missing expression {expression}"
                ))
            })?;
            pending.extend(occurrence.children.iter().copied());
            if let ExpressionKind::Power { copies, .. } = &occurrence.kind {
                pending.extend(copies.iter().copied());
            }
        }
        Ok(descendants)
    }

    fn validate_power_automorphisms(
        vertex_count: usize,
        powers: &[PowerDescriptor],
        orbit_generators: &[Vec<Vec<usize>>],
    ) -> Result<(), CanonicalizationError> {
        let roots = powers
            .iter()
            .enumerate()
            .flat_map(|(power, descriptor)| {
                descriptor
                    .copies
                    .iter()
                    .enumerate()
                    .map(move |(copy, descriptor)| (descriptor.root, (power, copy)))
            })
            .collect::<BTreeMap<_, _>>();
        if roots.len() != powers.iter().map(|power| power.copies.len()).sum::<usize>() {
            return Err(CanonicalizationError::Projection(
                "Power descriptors repeat a canonical copy root".into(),
            ));
        }

        let group = SignedGroup::from_graphica(vertex_count, &[], orbit_generators)
            .map_err(|error| CanonicalizationError::Projection(error.to_string()))?;
        // Whole-copy membership and boundary descriptors are closed under
        // composition, so checking Graphica's complete generator set proves
        // the property for every generated canonical automorphism.
        for generator in group.generators() {
            let vertex_map = generator.vertex_map();
            for source_power in powers {
                for source_copy in &source_power.copies {
                    let mapped_root = vertex_map[source_copy.root];
                    let Some(&(target_power_index, target_copy_index)) = roots.get(&mapped_root)
                    else {
                        return Err(CanonicalizationError::Projection(format!(
                            "canonical automorphism maps Power copy root {} outside the complete-copy set",
                            source_copy.root
                        )));
                    };
                    let target_power = &powers[target_power_index];
                    let target_copy = &target_power.copies[target_copy_index];
                    let mapped_vertices = source_copy
                        .vertices
                        .iter()
                        .map(|&vertex| vertex_map[vertex])
                        .collect::<BTreeSet<_>>();
                    if mapped_vertices != target_copy.vertices {
                        return Err(CanonicalizationError::Projection(format!(
                            "canonical automorphism maps only part of the complete Power copy rooted at {}",
                            source_copy.root
                        )));
                    }
                    if source_power.exponent != target_power.exponent
                        || vertex_map[source_power.result] != target_power.result
                        || vertex_map[source_power.magnitude] != target_power.magnitude
                        || vertex_map[source_copy.base_root] != target_copy.base_root
                    {
                        return Err(CanonicalizationError::Projection(format!(
                            "canonical automorphism does not preserve signed Power polarity and rooted copy scope at root {}",
                            source_copy.root
                        )));
                    }
                    let mut mapped_boundaries = source_copy
                        .boundaries
                        .iter()
                        .map(|boundary| boundary.mapped_infallible(vertex_map))
                        .collect::<Vec<_>>();
                    mapped_boundaries.sort();
                    if mapped_boundaries != target_copy.boundaries {
                        return Err(CanonicalizationError::Projection(format!(
                            "canonical automorphism changes the Power boundary pairing at copy root {}",
                            source_copy.root
                        )));
                    }
                }
            }
        }
        Ok(())
    }

    fn problem_identity(graph: &UnifiedGraph) -> ProblemIdentity {
        let nodes = graph
            .nodes()
            .iter()
            .map(|node| node.data.data.clone())
            .collect();
        let mut edges = graph
            .edges()
            .iter()
            .map(|edge| {
                let mut vertices = edge.vertices;
                if !edge.data.data.directed() && vertices.0 > vertices.1 {
                    vertices = (vertices.1, vertices.0);
                }
                (vertices.0, vertices.1, edge.data.data.clone())
            })
            .collect::<Vec<_>>();
        edges.sort();
        ProblemIdentity { nodes, edges }
    }
}

#[cfg(test)]
mod tests {
    use std::{
        collections::{BTreeMap, BTreeSet},
        time::{Duration, Instant},
    };

    use spenso::{
        bracket, cyclic,
        network::graph::{NetworkEdge, NetworkLeaf, NetworkNode, NetworkOp},
        slot,
        structure::{abstract_index::AbstractIndex, slot::IsAbstractSlot},
        tensor_symbol,
    };
    use symbolica::{
        atom::{Atom, AtomCore},
        function, symbol,
    };

    use super::{
        CanonicalPolicyNet, CanonicalizationError, DEFAULT_GRAPH_BUDGET, Estimate, EstimateState,
        ExpressionKind, GraphBudget, PowerBoundaryDescriptor, PowerBoundaryTarget,
        PowerCopyDescriptor, PowerDescriptor, ProblemIdentity, ProjectionBuilder, UnifiedEdgeColor,
        UnifiedGraph, UnifiedNodeColor, graphica_calls, project, reset_graphica_calls,
    };
    use crate::tensor::canonicalize::driver::DEFAULT_ITERATION_LIMIT;
    use crate::test_support::test_initialize;

    struct Phase6Measurement {
        estimate: Estimate,
        graphica_search: Option<Duration>,
        iterations: usize,
        canonicalization: Duration,
    }

    #[test]
    fn ordered_partial_group_is_one_flat_tensor_occurrence() {
        let rep = test_initialize().mink4;
        let first = slot!(rep, projection_partial_group_first);
        let second = slot!(rep, projection_partial_group_second);
        let tensor = tensor_symbol!(projection_partial_group_tensor);
        let expression = function!(
            tensor,
            Atom::var(symbol!("projection_partial_group_parameter")),
            cyclic!(first, second)
        );
        let policy = CanonicalPolicyNet::<AbstractIndex>::parse(expression).unwrap();
        assert_eq!(
            policy
                .network()
                .graph
                .graph
                .iter_nodes()
                .filter(|(_, _, node)| matches!(
                    node,
                    NetworkNode::Leaf(NetworkLeaf::LocalTensor(_))
                ))
                .count(),
            1
        );

        let projection = project(&policy, DEFAULT_GRAPH_BUDGET).unwrap();
        assert_eq!(projection.tensors.len(), 1);
        assert_eq!(projection.tensors[0].ports.len(), 2);
        let colors = projection
            .graph
            .nodes()
            .iter()
            .map(|node| &node.data.data)
            .collect::<Vec<_>>();
        assert_eq!(
            colors
                .iter()
                .filter(|color| matches!(color, UnifiedNodeColor::Tensor(_)))
                .count(),
            1
        );
        assert_eq!(
            colors
                .iter()
                .filter(|color| matches!(color, UnifiedNodeColor::Port(_)))
                .count(),
            2
        );
        assert_eq!(
            colors
                .iter()
                .filter(|color| matches!(color, UnifiedNodeColor::ExternalLine(_)))
                .count(),
            2
        );
    }

    #[test]
    fn scalar_multiplier_is_a_visible_product_child() {
        let rep = test_initialize().mink4;
        let index = slot!(rep, projection_visible_scalar_index);
        let tensor = tensor_symbol!(projection_visible_scalar_tensor);
        let expression = Atom::num(2) * function!(tensor, index.to_atom());
        let policy = CanonicalPolicyNet::<AbstractIndex>::parse(expression).unwrap();
        let projection = project(&policy, DEFAULT_GRAPH_BUDGET).unwrap();
        let only_vertex = |predicate: fn(&UnifiedNodeColor) -> bool| {
            let mut vertices = projection
                .graph
                .nodes()
                .iter()
                .enumerate()
                .filter_map(|(vertex, node)| predicate(&node.data.data).then_some(vertex));
            let vertex = vertices.next().unwrap();
            assert!(vertices.next().is_none());
            vertex
        };
        let root = only_vertex(|color| matches!(color, UnifiedNodeColor::Root));
        let product = only_vertex(|color| matches!(color, UnifiedNodeColor::Product));
        let scalar = only_vertex(|color| matches!(color, UnifiedNodeColor::Scalar(_)));
        let tensor = only_vertex(|color| matches!(color, UnifiedNodeColor::Tensor(_)));
        let has_edge = |source, target, color| {
            projection.graph.edges().iter().any(|edge| {
                edge.vertices == (source, target) && edge.directed && edge.data.data == color
            })
        };

        assert!(has_edge(root, product, UnifiedEdgeColor::Root));
        assert!(has_edge(product, scalar, UnifiedEdgeColor::Child));
        assert!(has_edge(product, tensor, UnifiedEdgeColor::Child));
        assert!(projection.expressions.iter().any(|occurrence| {
            matches!(&occurrence.kind, ExpressionKind::Scalar(atom) if atom == &Atom::num(2))
        }));
    }

    #[test]
    fn root_and_equivalent_nested_operators_have_distinct_roles() {
        let rep = test_initialize().mink4;
        let index = slot!(rep, projection_root_nested_index);
        let tensor = tensor_symbol!(projection_root_nested_tensor);
        let wrapper = symbol!("projection_root_nested_wrapper");
        let expression = function!(
            wrapper,
            function!(wrapper, function!(tensor, index.to_atom()))
        );
        let policy = CanonicalPolicyNet::<AbstractIndex>::parse(expression).unwrap();
        let projection = project(&policy, DEFAULT_GRAPH_BUDGET).unwrap();
        let root = projection
            .graph
            .nodes()
            .iter()
            .position(|node| matches!(&node.data.data, UnifiedNodeColor::Root))
            .unwrap();
        let functions = projection
            .graph
            .nodes()
            .iter()
            .enumerate()
            .filter_map(|(vertex, node)| {
                matches!(&node.data.data, UnifiedNodeColor::Function(_)).then_some(vertex)
            })
            .collect::<Vec<_>>();
        assert_eq!(functions.len(), 2);
        assert_eq!(
            projection.graph.node(functions[0]).data.data,
            projection.graph.node(functions[1]).data.data
        );
        let outer = projection
            .graph
            .edges()
            .iter()
            .find_map(|edge| {
                (edge.vertices.0 == root && edge.data.data == UnifiedEdgeColor::Root)
                    .then_some(edge.vertices.1)
            })
            .unwrap();
        assert!(functions.contains(&outer));
        let inner = functions
            .into_iter()
            .find(|&function| function != outer)
            .unwrap();

        assert!(projection.graph.edges().iter().any(|edge| {
            edge.vertices == (outer, inner)
                && edge.directed
                && edge.data.data == UnifiedEdgeColor::Argument
        }));
    }

    #[test]
    fn tensor_endpoint_conflicting_with_line_component_is_rejected() {
        let reps = test_initialize();
        let index = slot!(reps.mink4, projection_conflicting_endpoint);
        let tensor = tensor_symbol!(projection_conflicting_endpoint_tensor);
        let expression = function!(tensor, index.to_atom());
        let policy = CanonicalPolicyNet::<AbstractIndex>::parse(expression.clone()).unwrap();
        let mut network = policy.into_network();
        let graph = &network.graph.graph;
        let tensor_node = graph
            .iter_nodes()
            .find_map(|(node, _, data)| {
                matches!(data, NetworkNode::Leaf(NetworkLeaf::LocalTensor(_))).then_some(node)
            })
            .unwrap();
        let hedge = graph
            .iter_crown(tensor_node)
            .find(|hedge| matches!(graph[[hedge]], NetworkEdge::Slot(_)))
            .unwrap();
        let conflicting = reps.bis4.slot::<AbstractIndex, _>(index.aind()).to_lib();
        network.graph.graph[[&hedge]] = NetworkEdge::Slot(conflicting);

        let Err(error) = ProjectionBuilder::new(&network, &expression, DEFAULT_GRAPH_BUDGET) else {
            panic!("conflicting tensor endpoint must fail projection setup");
        };
        assert!(matches!(error, CanonicalizationError::Projection(message)
            if message.contains("tensor endpoint")
                && message.contains("line-component declaration")));
    }

    #[test]
    fn same_abstract_index_in_different_representation_groups_stays_distinct() {
        let reps = test_initialize();
        let minkowski = slot!(reps.mink4, projection_shared_aind_different_groups);
        let bispinor = reps.bis4.slot::<AbstractIndex, _>(minkowski.aind());
        let left = tensor_symbol!(projection_shared_aind_minkowski);
        let right = tensor_symbol!(projection_shared_aind_bispinor);
        let expression =
            function!(left, minkowski.to_atom()) * function!(right, bispinor.to_atom());
        let policy = CanonicalPolicyNet::<AbstractIndex>::parse(expression).unwrap();
        let projection = project(&policy, DEFAULT_GRAPH_BUDGET).unwrap();
        let groups = projection
            .lines
            .iter()
            .map(|line| line.group)
            .collect::<BTreeSet<_>>();

        assert_eq!(projection.lines.len(), 2);
        assert_eq!(groups.len(), 2);
        assert!(groups.contains(&reps.mink4.to_lib().base()));
        assert!(groups.contains(&reps.bis4.to_lib().base()));
    }

    fn phase6_corpus() -> Vec<(&'static str, Atom)> {
        let rep = test_initialize().mink4;
        let slots = (0..27)
            .map(|index| slot!(rep, AbstractIndex::Dummy(index)))
            .collect::<Vec<_>>();

        let antisymmetric = tensor_symbol!(phase6_benchmark_antisymmetric; Antisymmetric);
        let f = |a: usize, b: usize, c: usize| {
            function!(
                antisymmetric,
                slots[a].to_atom(),
                slots[b].to_atom(),
                slots[c].to_atom()
            )
        };
        let k33_component = |offset: usize| {
            f(offset, offset + 1, offset + 2)
                * f(offset + 3, offset + 4, offset + 5)
                * f(offset + 6, offset + 7, offset + 8)
                * f(offset, offset + 3, offset + 6)
                * f(offset + 1, offset + 4, offset + 7)
                * f(offset + 2, offset + 5, offset + 8)
        };
        let k33 = k33_component(0);
        let three_k33 = k33.clone() * k33_component(9) * k33_component(18);

        let symmetric = tensor_symbol!(phase6_benchmark_symmetric; Symmetric);
        let left = tensor_symbol!(phase6_benchmark_left);
        let right = tensor_symbol!(phase6_benchmark_right);
        let branch = |a: usize, b: usize| {
            function!(symmetric, slots[a].to_atom(), slots[b].to_atom())
                * function!(left, slots[a].to_atom())
                * function!(right, slots[b].to_atom())
        };
        let repeated_sum = branch(0, 1) + branch(2, 3) + branch(4, 5) + branch(6, 7);

        let power_tensor = tensor_symbol!(phase6_benchmark_power; Symmetric);
        let power_base = function!(power_tensor, slots[0].to_atom(), slots[1].to_atom());

        vec![
            ("signed_k33", k33),
            ("three_signed_k33_components", three_k33),
            ("repeated_symmetric_sum", repeated_sum),
            ("positive_power_4", power_base.clone().pow(4)),
            ("negative_power_4", power_base.clone().pow(-4)),
            ("positive_power_24", power_base.pow(24)),
        ]
    }

    fn measure_phase6_case(expression: Atom, measure_search: bool) -> Phase6Measurement {
        let policy = CanonicalPolicyNet::<AbstractIndex>::parse(expression).unwrap();
        let mut builder = ProjectionBuilder::new(
            policy.network(),
            policy.normalized_atom(),
            GraphBudget {
                vertices: usize::MAX,
                edges: usize::MAX,
            },
        )
        .unwrap();
        let root = builder
            .network
            .graph
            .graph
            .node_id(builder.network.graph.head());
        let estimate = builder.estimate_visible_graph(root).unwrap();
        builder.project_visible_graph().unwrap();
        assert_eq!(builder.graph.nodes().len(), estimate.vertices);
        assert_eq!(builder.graph.edges().len(), estimate.edges);
        let graphica_search = measure_search.then(|| {
            let start = Instant::now();
            drop(builder.graph.canonize());
            start.elapsed()
        });
        drop(builder);

        reset_graphica_calls();
        let start = Instant::now();
        policy.canonize(AbstractIndex::Dummy).unwrap();
        let canonicalization = start.elapsed();
        let iterations = graphica_calls();
        assert!(iterations > 0);

        Phase6Measurement {
            estimate,
            graphica_search,
            iterations,
            canonicalization,
        }
    }

    fn assert_exact_estimate(expression: Atom) -> Estimate {
        let policy = CanonicalPolicyNet::<AbstractIndex>::parse(expression).unwrap();
        let mut builder = ProjectionBuilder::new(
            policy.network(),
            policy.normalized_atom(),
            GraphBudget {
                vertices: usize::MAX,
                edges: usize::MAX,
            },
        )
        .unwrap();
        let root = builder
            .network
            .graph
            .graph
            .node_id(builder.network.graph.head());
        let estimate = builder.estimate_visible_graph(root).unwrap();
        builder.project_visible_graph().unwrap();
        assert_eq!(builder.graph.nodes().len(), estimate.vertices);
        assert_eq!(builder.graph.edges().len(), estimate.edges);
        estimate
    }

    fn paired_power() -> PowerDescriptor {
        PowerDescriptor {
            expression: 0,
            exponent: 2,
            result: 0,
            magnitude: 1,
            copies: vec![
                PowerCopyDescriptor {
                    expression: 1,
                    root: 2,
                    base_root: 4,
                    vertices: [2, 4, 6].into_iter().collect(),
                    boundaries: vec![PowerBoundaryDescriptor {
                        line: 8,
                        ports: vec![6],
                        target: PowerBoundaryTarget::Pair {
                            partner_copy_root: 3,
                            partner_ports: vec![7],
                        },
                    }],
                },
                PowerCopyDescriptor {
                    expression: 2,
                    root: 3,
                    base_root: 5,
                    vertices: [3, 5, 7].into_iter().collect(),
                    boundaries: vec![PowerBoundaryDescriptor {
                        line: 8,
                        ports: vec![7],
                        target: PowerBoundaryTarget::Pair {
                            partner_copy_root: 2,
                            partner_ports: vec![6],
                        },
                    }],
                },
            ],
        }
    }

    fn shifted_power(mut power: PowerDescriptor, offset: usize, exponent: i8) -> PowerDescriptor {
        power.expression += 10;
        power.exponent = exponent;
        power.result += offset;
        power.magnitude += offset;
        for copy in &mut power.copies {
            copy.expression += 10;
            copy.root += offset;
            copy.base_root += offset;
            copy.vertices = copy
                .vertices
                .iter()
                .map(|vertex| *vertex + offset)
                .collect();
            for boundary in &mut copy.boundaries {
                boundary.line += offset;
                boundary.ports.iter_mut().for_each(|port| *port += offset);
                match &mut boundary.target {
                    PowerBoundaryTarget::Pair {
                        partner_copy_root,
                        partner_ports,
                    } => {
                        *partner_copy_root += offset;
                        partner_ports.iter_mut().for_each(|port| *port += offset);
                    }
                    PowerBoundaryTarget::Interface { outside_ports, .. } => {
                        outside_ports.iter_mut().for_each(|port| *port += offset)
                    }
                }
            }
        }
        power
    }

    #[test]
    fn phase6_production_budgets_cover_representative_corpus() {
        for (name, expression) in phase6_corpus() {
            let measurement = measure_phase6_case(expression, false);
            assert!(
                measurement.estimate.vertices <= DEFAULT_GRAPH_BUDGET.vertices,
                "{name} needs {} vertices, exceeding the production limit {}",
                measurement.estimate.vertices,
                DEFAULT_GRAPH_BUDGET.vertices,
            );
            assert!(
                measurement.estimate.edges <= DEFAULT_GRAPH_BUDGET.edges,
                "{name} needs {} edges, exceeding the production limit {}",
                measurement.estimate.edges,
                DEFAULT_GRAPH_BUDGET.edges,
            );
            assert!(
                measurement.iterations <= DEFAULT_ITERATION_LIMIT,
                "{name} needs {} iterations, exceeding the production limit {}",
                measurement.iterations,
                DEFAULT_ITERATION_LIMIT,
            );
        }
    }

    #[test]
    #[ignore = "manual Phase-6 canonicalization budget report"]
    fn phase6_canonicalization_budget_report() {
        println!("case\tvertices\tedges\tgraphica_ms\titerations\tcanonicalization_ms");
        for (name, expression) in phase6_corpus() {
            let measurement = measure_phase6_case(expression, true);
            let graphica_ms = measurement
                .graphica_search
                .expect("the benchmark measures Graphica directly")
                .as_secs_f64()
                * 1_000.0;
            println!(
                "{}\t{}\t{}\t{:.3}\t{}\t{:.3}",
                name,
                measurement.estimate.vertices,
                measurement.estimate.edges,
                graphica_ms,
                measurement.iterations,
                measurement.canonicalization.as_secs_f64() * 1_000.0,
            );
        }
        println!(
            "production_limits\t{}\t{}\t-\t{}\t-",
            DEFAULT_GRAPH_BUDGET.vertices, DEFAULT_GRAPH_BUDGET.edges, DEFAULT_ITERATION_LIMIT,
        );
    }

    #[test]
    fn exact_estimate_counts_shared_lines_and_intrinsic_cycle_edges() {
        let rep = test_initialize().mink4;
        let a = slot!(rep, exact_estimate_cycle_a);
        let b = slot!(rep, exact_estimate_cycle_b);
        let c = slot!(rep, exact_estimate_cycle_c);
        let cycle = tensor_symbol!(exact_estimate_cycle; Cyclesymmetric);
        let leaf = tensor_symbol!(exact_estimate_leaf);
        let expression = function!(cycle, a.to_atom(), b.to_atom(), c.to_atom())
            * function!(leaf, a.to_atom())
            * function!(leaf, b.to_atom())
            * function!(leaf, c.to_atom());

        assert_eq!(
            assert_exact_estimate(expression),
            Estimate {
                vertices: 15,
                edges: 20,
                largest_power: None,
            }
        );
    }

    #[test]
    fn exact_estimate_matches_power_line_topologies() {
        let rep = test_initialize().mink4;
        let internal_index = slot!(rep, exact_estimate_internal_power_line);
        let external_index = slot!(rep, exact_estimate_external_power_line);
        let nested_index = slot!(rep, exact_estimate_nested_power_line);
        let internal_left = tensor_symbol!(exact_estimate_internal_power_left);
        let internal_right = tensor_symbol!(exact_estimate_internal_power_right);
        let external = tensor_symbol!(exact_estimate_external_power_tensor);
        let nested = tensor_symbol!(exact_estimate_nested_power_tensor);
        let internal_base = function!(internal_left, internal_index.to_atom())
            * function!(internal_right, internal_index.to_atom());
        let external_base = function!(external, external_index.to_atom());
        let nested_base = function!(nested, nested_index.to_atom());
        let mut cases = Vec::new();

        for magnitude in 2..=5 {
            cases.push(bracket!(internal_base.clone()).pow(magnitude));
            cases.push(bracket!(external_base.clone()).pow(magnitude));
        }
        cases.push(bracket!(nested_base.pow(2)).pow(3));

        for expression in cases {
            assert_exact_estimate(expression);
        }
    }

    #[test]
    fn signed_power_projection_matrix_has_complete_copy_and_pairing_structure() {
        fn unsigned_identity(mut graph: UnifiedGraph) -> ProblemIdentity {
            for vertex in 0..graph.nodes().len() {
                let mut data = graph.node(vertex).data.clone();
                if let UnifiedNodeColor::PowerResult(exponent) = &mut data.data {
                    *exponent = exponent.unsigned_abs() as i8;
                }
                graph.set_node_data(vertex, data);
            }
            ProjectionBuilder::<'_, AbstractIndex>::problem_identity(&graph.canonize().graph)
        }

        let rep = test_initialize().mink4;
        let first = slot!(rep, projection_power_matrix_first);
        let second = slot!(rep, projection_power_matrix_second);
        let tensor = tensor_symbol!(projection_power_matrix_tensor; Symmetric);
        let base = function!(tensor, first.to_atom(), second.to_atom());
        let mut unsigned_identities = BTreeMap::new();

        for exponent in [-4_i8, -2, 2, 3, 4, 5] {
            let expression = base.clone().pow(exponent);
            let policy = CanonicalPolicyNet::<AbstractIndex>::parse(expression).unwrap();
            let projection = project(&policy, DEFAULT_GRAPH_BUDGET).unwrap();
            assert_eq!(projection.powers.len(), 1);
            let power = &projection.powers[0];
            let magnitude = usize::from(exponent.unsigned_abs());

            assert_eq!(power.exponent, exponent);
            assert_eq!(power.copies.len(), magnitude);
            assert_eq!(
                projection
                    .graph
                    .nodes()
                    .iter()
                    .filter(|node| matches!(&node.data.data, UnifiedNodeColor::PowerMagnitude))
                    .count(),
                1
            );
            assert_eq!(
                projection
                    .graph
                    .nodes()
                    .iter()
                    .filter(|node| matches!(&node.data.data, UnifiedNodeColor::PowerCopy))
                    .count(),
                magnitude
            );
            assert!(matches!(
                &projection.graph.node(power.result).data.data,
                UnifiedNodeColor::PowerResult(projected) if *projected == exponent
            ));
            assert!(matches!(
                &projection.graph.node(power.magnitude).data.data,
                UnifiedNodeColor::PowerMagnitude
            ));

            let mut partners = BTreeMap::new();
            let mut interface_copies = BTreeSet::new();
            for copy in &power.copies {
                assert_eq!(copy.boundaries.len(), 2);
                let pair_roots = copy
                    .boundaries
                    .iter()
                    .filter_map(|boundary| match &boundary.target {
                        PowerBoundaryTarget::Pair {
                            partner_copy_root, ..
                        } => Some(*partner_copy_root),
                        PowerBoundaryTarget::Interface { .. } => None,
                    })
                    .collect::<BTreeSet<_>>();
                let interface_count = copy
                    .boundaries
                    .iter()
                    .filter(|boundary| {
                        matches!(&boundary.target, PowerBoundaryTarget::Interface { .. })
                    })
                    .count();
                match (pair_roots.len(), interface_count) {
                    (1, 0) => {
                        partners.insert(copy.root, *pair_roots.first().unwrap());
                    }
                    (0, 2) => {
                        interface_copies.insert(copy.root);
                    }
                    state => panic!(
                        "copy {} has inconsistent boundary pairing state {state:?}",
                        copy.root
                    ),
                }
            }
            for (&copy, &partner) in &partners {
                assert_eq!(partners.get(&partner), Some(&copy));
            }
            if magnitude.is_multiple_of(2) {
                assert!(interface_copies.is_empty());
                assert_eq!(partners.len(), magnitude);
            } else {
                assert_eq!(interface_copies.len(), 1);
                assert_eq!(partners.len(), magnitude - 1);
            }

            unsigned_identities.insert(exponent, unsigned_identity(projection.graph));
        }

        for magnitude in [2, 4] {
            assert_eq!(
                unsigned_identities[&magnitude],
                unsigned_identities[&-magnitude]
            );
        }
    }

    #[test]
    fn exact_nested_power_estimate_enforces_injected_limits_before_projection() {
        let rep = test_initialize().mink4;
        let inner_index = slot!(rep, exact_estimate_inner);
        let outer_index = slot!(rep, exact_estimate_outer);
        let inner = tensor_symbol!(exact_estimate_inner_tensor);
        let outer = tensor_symbol!(exact_estimate_outer_tensor);
        let nested = function!(inner, inner_index.to_atom()).pow(2);
        let expression = bracket!(nested * function!(outer, outer_index.to_atom())).pow(3);
        let estimate = assert_exact_estimate(expression.clone());
        assert_eq!(estimate.largest_power, Some((3, 3)));

        for budget in [
            GraphBudget {
                vertices: estimate.vertices - 1,
                edges: estimate.edges,
            },
            GraphBudget {
                vertices: estimate.vertices,
                edges: estimate.edges - 1,
            },
        ] {
            let policy = CanonicalPolicyNet::<AbstractIndex>::parse(expression.clone()).unwrap();
            let mut builder =
                ProjectionBuilder::new(policy.network(), policy.normalized_atom(), budget).unwrap();
            let error = builder.project_visible_graph().unwrap_err();
            assert!(matches!(
                error,
                CanonicalizationError::GraphSizeLimit {
                    power: 3,
                    magnitude: 3,
                    requested_vertices,
                    requested_edges,
                    ..
                } if requested_vertices == estimate.vertices
                    && requested_edges == estimate.edges
            ));
            assert!(builder.graph.nodes().is_empty());
            assert!(builder.graph.edges().is_empty());
        }
    }

    #[test]
    fn minimum_signed_power_uses_unsigned_magnitude_and_preflights_before_projection() {
        let rep = test_initialize().mink4;
        let first = slot!(rep, minimum_signed_power_first);
        let second = slot!(rep, minimum_signed_power_second);
        let tensor = tensor_symbol!(minimum_signed_power_tensor; Symmetric);
        let expression = function!(tensor, first.to_atom(), second.to_atom()).pow(i8::MIN);
        let estimate = assert_exact_estimate(expression.clone());
        assert_eq!((estimate.vertices, estimate.edges), (643, 770));
        assert_eq!(estimate.largest_power, Some((i8::MIN, 128)));

        for budget in [
            DEFAULT_GRAPH_BUDGET,
            GraphBudget {
                vertices: estimate.vertices - 1,
                edges: estimate.edges,
            },
            GraphBudget {
                vertices: estimate.vertices,
                edges: estimate.edges - 1,
            },
        ] {
            let policy = CanonicalPolicyNet::<AbstractIndex>::parse(expression.clone()).unwrap();
            let mut builder =
                ProjectionBuilder::new(policy.network(), policy.normalized_atom(), budget).unwrap();
            let error = builder.project_visible_graph().unwrap_err();
            assert!(matches!(
                error,
                CanonicalizationError::GraphSizeLimit {
                    power: i8::MIN,
                    magnitude: 128,
                    requested_vertices,
                    requested_edges,
                    ..
                } if requested_vertices == estimate.vertices
                    && requested_edges == estimate.edges
            ));
            assert!(builder.graph.nodes().is_empty());
            assert!(builder.graph.edges().is_empty());
        }
    }

    #[test]
    #[ignore = "manual multi-minute Graphica check for the full signed i8::MIN domain"]
    fn minimum_signed_power_projects_with_a_sufficient_injected_budget() {
        let rep = test_initialize().mink4;
        let first = slot!(rep, minimum_signed_power_graphica_first);
        let second = slot!(rep, minimum_signed_power_graphica_second);
        let tensor = tensor_symbol!(minimum_signed_power_graphica_tensor; Symmetric);
        let expression = function!(tensor, first.to_atom(), second.to_atom()).pow(i8::MIN);
        let estimate = assert_exact_estimate(expression.clone());
        let policy = CanonicalPolicyNet::<AbstractIndex>::parse(expression).unwrap();
        let start = Instant::now();

        let projection = super::project(
            &policy,
            GraphBudget {
                vertices: estimate.vertices,
                edges: estimate.edges,
            },
        )
        .unwrap();

        assert_eq!(projection.graph.nodes().len(), estimate.vertices);
        assert_eq!(projection.graph.edges().len(), estimate.edges);
        println!(
            "Power(i8::MIN)\t{}\t{}\t{:.3} ms",
            estimate.vertices,
            estimate.edges,
            start.elapsed().as_secs_f64() * 1_000.0,
        );
    }

    #[test]
    fn complete_pair_swap_preserves_power_descriptor() {
        let power = paired_power();
        power.validate_local().unwrap();
        ProjectionBuilder::<'_, AbstractIndex>::validate_power_automorphisms(
            9,
            &[power],
            &[vec![vec![2, 3], vec![4, 5], vec![6, 7]]],
        )
        .unwrap();
    }

    #[test]
    fn complete_pair_of_pairs_swap_ignores_descriptor_copy_order() {
        let mut power = paired_power();
        power.exponent = 4;
        power
            .copies
            .extend(shifted_power(paired_power(), 7, 4).copies);
        let generator = vec![
            vec![2, 9],
            vec![3, 10],
            vec![4, 11],
            vec![5, 12],
            vec![6, 13],
            vec![7, 14],
            vec![8, 15],
        ];
        power.validate_local().unwrap();
        ProjectionBuilder::<'_, AbstractIndex>::validate_power_automorphisms(
            16,
            std::slice::from_ref(&power),
            std::slice::from_ref(&generator),
        )
        .unwrap();

        let mut reordered = power;
        reordered.copies.reverse();
        reordered.validate_local().unwrap();
        ProjectionBuilder::<'_, AbstractIndex>::validate_power_automorphisms(
            16,
            &[reordered],
            &[generator],
        )
        .unwrap();
    }

    #[test]
    fn different_copy_partners_across_boundary_slots_are_rejected() {
        let mut power = paired_power();
        power.exponent = 4;
        power.copies[0].vertices.insert(16);
        power.copies[0].boundaries.push(PowerBoundaryDescriptor {
            line: 15,
            ports: vec![16],
            target: PowerBoundaryTarget::Pair {
                partner_copy_root: 9,
                partner_ports: vec![13],
            },
        });
        power.copies.push(PowerCopyDescriptor {
            expression: 3,
            root: 9,
            base_root: 11,
            vertices: [9, 11, 13].into_iter().collect(),
            boundaries: vec![PowerBoundaryDescriptor {
                line: 15,
                ports: vec![13],
                target: PowerBoundaryTarget::Pair {
                    partner_copy_root: 2,
                    partner_ports: vec![16],
                },
            }],
        });
        power.copies.push(PowerCopyDescriptor {
            expression: 4,
            root: 10,
            base_root: 12,
            vertices: [10, 12, 14].into_iter().collect(),
            boundaries: Vec::new(),
        });

        let error = power.validate_local().unwrap_err();
        assert!(matches!(error, CanonicalizationError::Projection(message)
            if message.contains("inconsistent partner copies across boundary slots")));
    }

    #[test]
    fn partial_copy_automorphism_is_rejected() {
        let error = ProjectionBuilder::<'_, AbstractIndex>::validate_power_automorphisms(
            9,
            &[paired_power()],
            &[vec![vec![2, 3]]],
        )
        .unwrap_err();
        assert!(matches!(error, CanonicalizationError::Projection(message)
            if message.contains("only part of the complete Power copy")));
    }

    #[test]
    fn changed_boundary_pairing_is_rejected() {
        let error = ProjectionBuilder::<'_, AbstractIndex>::validate_power_automorphisms(
            10,
            &[paired_power()],
            &[vec![vec![2, 3], vec![4, 5], vec![6, 7], vec![8, 9]]],
        )
        .unwrap_err();
        assert!(matches!(error, CanonicalizationError::Projection(message)
            if message.contains("changes the Power boundary pairing")));
    }

    #[test]
    fn signed_power_polarity_is_not_exchangeable() {
        let positive = paired_power();
        let negative = shifted_power(positive.clone(), 9, -2);
        let error = ProjectionBuilder::<'_, AbstractIndex>::validate_power_automorphisms(
            18,
            &[positive, negative],
            &[vec![
                vec![0, 9],
                vec![1, 10],
                vec![2, 11],
                vec![3, 12],
                vec![4, 13],
                vec![5, 14],
                vec![6, 15],
                vec![7, 16],
                vec![8, 17],
            ]],
        )
        .unwrap_err();
        assert!(matches!(error, CanonicalizationError::Projection(message)
            if message.contains("signed Power polarity")));
    }

    #[test]
    fn power_free_graph_limit_has_no_synthetic_power_context() {
        let expression = Atom::num(2);
        let policy = CanonicalPolicyNet::<AbstractIndex>::parse(expression.clone()).unwrap();
        let mut builder = ProjectionBuilder::new(
            policy.network(),
            policy.normalized_atom(),
            GraphBudget {
                vertices: 1,
                edges: 1,
            },
        )
        .unwrap();

        assert!(matches!(
            builder.project_visible_graph(),
            Err(CanonicalizationError::WholeGraphSizeLimit {
                expression: limited,
                requested_vertices: 2,
                requested_edges: 1,
                ..
            }) if limited == expression
        ));
        assert!(builder.graph.nodes().is_empty());
        assert!(builder.graph.edges().is_empty());
    }

    #[test]
    fn estimator_overflow_is_a_limit_at_usize_max_budget() {
        let expression = Atom::num(2);
        let policy = CanonicalPolicyNet::<AbstractIndex>::parse(expression).unwrap();
        let builder = ProjectionBuilder::new(
            policy.network(),
            policy.normalized_atom(),
            GraphBudget {
                vertices: usize::MAX,
                edges: usize::MAX,
            },
        )
        .unwrap();

        for (vertices, edges, added_vertices, added_edges) in
            [(usize::MAX, 0, 1, 0), (0, usize::MAX, 0, 1)]
        {
            let mut state = EstimateState::default();
            state.size.vertices = vertices;
            state.size.edges = edges;
            assert!(matches!(
                builder.add_estimated_size(&mut state, added_vertices, added_edges),
                Err(CanonicalizationError::WholeGraphSizeLimit { .. })
            ));
        }
    }

    #[test]
    fn deeply_nested_minimum_powers_stop_at_the_injected_budget() {
        let rep = test_initialize().mink4;
        let index = slot!(rep, nested_minimum_power_index);
        let tensor = tensor_symbol!(nested_minimum_power_tensor; Symmetric);
        let wrapper = symbol!("nested_minimum_power_wrapper");
        let mut expression = function!(tensor, index.to_atom(), index.to_atom()).pow(i8::MIN);
        for _ in 0..4 {
            expression = function!(wrapper, expression).pow(i8::MIN);
        }
        let policy = CanonicalPolicyNet::<AbstractIndex>::parse(expression).unwrap();
        assert_eq!(
            policy
                .network()
                .graph
                .graph
                .iter_nodes()
                .filter(|(_, _, node)| {
                    matches!(node, NetworkNode::Op(NetworkOp::Power(i8::MIN)))
                })
                .count(),
            5
        );
        let budget = GraphBudget {
            vertices: 64,
            edges: usize::MAX,
        };
        let mut builder =
            ProjectionBuilder::new(policy.network(), policy.normalized_atom(), budget).unwrap();

        assert!(matches!(
            builder.project_visible_graph(),
            Err(CanonicalizationError::GraphSizeLimit {
                power: i8::MIN,
                magnitude: 128,
                requested_vertices: 172_880_888_323,
                requested_edges: 207_240_626_690,
                ..
            })
        ));
        assert!(builder.graph.nodes().is_empty());
        assert!(builder.graph.edges().is_empty());
    }
}
