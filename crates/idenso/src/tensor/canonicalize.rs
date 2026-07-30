//! Signed canonicalization for Spenso symbolic tensor networks.
//!
//! The network remains the semantic owner of products, sums, powers, and
//! enclosing functions.  Within one multiplicative scope this module projects
//! every local tensor to one header and one vertex per flat structural slot.
//! Graphica then supplies both the canonical labeling and stabilizer
//! generators.  Hidden graph data records only reconstruction origins; it does
//! not participate in graph identity.

use std::collections::{BTreeMap, BTreeSet};

use linnet::{
    half_edge::{
        NodeIndex,
        involution::{Hedge, HedgePair},
        tree::SimpleTraversalTree,
    },
    permutation::Permutation,
    tree::child_pointer::ParentChildStore,
    union_find::UnionFind,
};
use spenso::{
    network::{
        Network,
        graph::{NetworkEdge, NetworkLeaf, NetworkNode, NetworkOp, ScalarRef},
        store::TensorScalarStore,
        tags::SPENSO_TAG,
    },
    shadowing::{ANTISYM, CYCLIC, SYM},
    structure::{
        OrderedStructure, TensorStructure,
        abstract_index::AIND_SYMBOLS,
        representation::{LibraryRep, LibrarySlot, Representation},
        slot::{AbsInd, DummyAind, IsAbstractSlot, ParseableAind, Slot},
    },
};
use symbolica::{
    atom::{Atom, AtomView, FunctionBuilder, Symbol},
    domains::rational::Rational,
    graph::{Graph, HiddenData},
};
use thiserror::Error;

use crate::shorthands::{metric::MetricSimplifier, schoonschip::Schoonschip};

use super::{SymbolicNet, SymbolicTensor};

/// Failures detected while validating or canonicalizing a symbolic network.
#[derive(Clone, Debug, Error, PartialEq, Eq)]
pub enum CanonicalizationError {
    #[error(
        "intrinsically symmetric tensor {head} requires a direct slot at argument {argument}: {expression}"
    )]
    InvalidIntrinsicArgument {
        head: Symbol,
        argument: usize,
        expression: Atom,
    },
    #[error(
        "structural symmetry group {group} requires a direct slot at member {member}: {expression}"
    )]
    InvalidPartialGroup {
        group: Symbol,
        member: usize,
        expression: Atom,
    },
    #[error("tensor layout does not match its flat network structure: {expression}")]
    StructureMismatch { expression: Atom },
    #[error("unsupported tensor power {power} in canonicalization")]
    UnsupportedPower { power: i8 },
    #[error("unsupported symbolic-network leaf in canonicalization")]
    UnsupportedLeaf,
    #[error("canonical sign cannot be assigned to one nonlinear scope")]
    AmbiguousSignScope,
    #[error("failed to parse symbolic tensor network: {0}")]
    Network(String),
    #[error("failed to reversibly encode a structured index: {0}")]
    StructuredIndex(String),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum SymmetryKind {
    Symmetric,
    Antisymmetric,
    Cyclic,
}

impl SymmetryKind {
    fn of(symbol: Symbol) -> Option<Self> {
        if symbol.is_symmetric() {
            Some(Self::Symmetric)
        } else if symbol.is_antisymmetric() {
            Some(Self::Antisymmetric)
        } else if symbol.is_cyclesymmetric() {
            Some(Self::Cyclic)
        } else {
            None
        }
    }

    fn structural(symbol: Symbol) -> Option<Self> {
        if symbol == *SYM {
            Some(Self::Symmetric)
        } else if symbol == *ANTISYM {
            Some(Self::Antisymmetric)
        } else if symbol == *CYCLIC {
            Some(Self::Cyclic)
        } else {
            None
        }
    }
}

#[derive(Clone, Debug)]
enum LayoutArgument {
    Opaque(Atom),
    DirectSlot(usize),
    SlotBundle {
        holes: Vec<usize>,
    },
    Group {
        symbol: Symbol,
        kind: SymmetryKind,
        prefix: Atom,
        holes: Vec<usize>,
    },
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum LayoutColorArgument {
    Opaque(Atom),
    DirectSlot,
    SlotBundle(usize),
    Group {
        kind: SymmetryKind,
        members: usize,
        prefix: Atom,
    },
}

#[derive(Clone, Debug)]
struct TensorLayout {
    head: Symbol,
    arguments: Vec<LayoutArgument>,
    color: TensorColor,
    slot_count: usize,
    intrinsic: Option<SymmetryKind>,
    outer_linear: bool,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct TensorColor {
    head: Symbol,
    arguments: Vec<LayoutColorArgument>,
}

type PartialGroup<Aind> = (Symbol, SymmetryKind, Atom, Vec<LibrarySlot<Aind>>);

impl TensorLayout {
    fn scan<Aind: AbsInd + ParseableAind>(
        tensor: &SymbolicTensor<Aind>,
    ) -> Result<Self, CanonicalizationError> {
        let AtomView::Fun(function) = tensor.expression.as_view() else {
            return Err(CanonicalizationError::StructureMismatch {
                expression: tensor.expression.clone(),
            });
        };
        let head = function.get_symbol();
        let tagged = head.has_tag(&SPENSO_TAG.tensor);
        let intrinsic = tagged.then(|| SymmetryKind::of(head)).flatten();
        let root_group = SymmetryKind::structural(head);

        let mut arguments = Vec::with_capacity(function.get_nargs());
        let mut colors = Vec::with_capacity(function.get_nargs());
        let mut expression_slots = Vec::new();
        let mut slot_count = 0;

        if intrinsic.is_some() || root_group.is_some() {
            for (argument, value) in function.iter().enumerate() {
                let Ok(slot) = Slot::<LibraryRep, Aind>::try_from(value) else {
                    return Err(CanonicalizationError::InvalidIntrinsicArgument {
                        head,
                        argument,
                        expression: tensor.expression.clone(),
                    });
                };
                arguments.push(LayoutArgument::DirectSlot(slot_count));
                colors.push(LayoutColorArgument::DirectSlot);
                expression_slots.push(slot);
                slot_count += 1;
            }
        } else {
            for value in function.iter() {
                if let Ok(slot) = Slot::<LibraryRep, Aind>::try_from(value) {
                    arguments.push(LayoutArgument::DirectSlot(slot_count));
                    colors.push(LayoutColorArgument::DirectSlot);
                    expression_slots.push(slot);
                    slot_count += 1;
                    continue;
                }

                if let AtomView::Fun(bundle) = value
                    && bundle.get_symbol() == AIND_SYMBOLS.aind
                {
                    let start = slot_count;
                    for member in bundle.iter() {
                        let Ok(slot) = Slot::<LibraryRep, Aind>::try_from(member) else {
                            return Err(CanonicalizationError::StructureMismatch {
                                expression: tensor.expression.clone(),
                            });
                        };
                        expression_slots.push(slot);
                        slot_count += 1;
                    }
                    let holes = (start..slot_count).collect::<Vec<_>>();
                    colors.push(LayoutColorArgument::SlotBundle(holes.len()));
                    arguments.push(LayoutArgument::SlotBundle { holes });
                    continue;
                }

                if tagged
                    && let Some((symbol, kind, prefix, members)) =
                        Self::partial_group::<Aind>(value, &tensor.expression)?
                {
                    let start = slot_count;
                    slot_count += members.len();
                    expression_slots.extend(members);
                    let holes = (start..slot_count).collect::<Vec<_>>();
                    colors.push(LayoutColorArgument::Group {
                        kind,
                        members: holes.len(),
                        prefix: prefix.clone(),
                    });
                    arguments.push(LayoutArgument::Group {
                        symbol,
                        kind,
                        prefix,
                        holes,
                    });
                    continue;
                }

                arguments.push(LayoutArgument::Opaque(value.to_owned()));
                colors.push(LayoutColorArgument::Opaque(value.to_owned()));
            }
        }

        if slot_count != tensor.structure.order() {
            return Err(CanonicalizationError::StructureMismatch {
                expression: tensor.expression.clone(),
            });
        }
        let mut structural_slots = tensor
            .structure
            .external_structure_iter()
            .enumerate()
            .collect::<Vec<_>>();
        let mut expression_to_structure = Vec::with_capacity(expression_slots.len());
        for slot in expression_slots {
            let Some(position) = structural_slots
                .iter()
                .position(|(_, structural)| structural == &slot)
            else {
                return Err(CanonicalizationError::StructureMismatch {
                    expression: tensor.expression.clone(),
                });
            };
            expression_to_structure.push(structural_slots.remove(position).0);
        }
        for argument in &mut arguments {
            match argument {
                LayoutArgument::DirectSlot(hole) => *hole = expression_to_structure[*hole],
                LayoutArgument::SlotBundle { holes } | LayoutArgument::Group { holes, .. } => {
                    for hole in holes {
                        *hole = expression_to_structure[*hole];
                    }
                }
                LayoutArgument::Opaque(_) => {}
            }
        }

        Ok(Self {
            head,
            arguments,
            color: TensorColor {
                head,
                arguments: colors,
            },
            slot_count,
            intrinsic: intrinsic.or(root_group),
            outer_linear: head.is_linear(),
        })
    }

    fn partial_group<Aind: AbsInd + ParseableAind>(
        value: AtomView<'_>,
        expression: &Atom,
    ) -> Result<Option<PartialGroup<Aind>>, CanonicalizationError> {
        let (prefix, group) = match value {
            AtomView::Fun(group)
                if group.get_nargs() > 0
                    && SymmetryKind::structural(group.get_symbol()).is_some() =>
            {
                (Atom::num(1), group)
            }
            AtomView::Mul(product) => {
                let mut factors = product.iter();
                let Some(first) = factors.next() else {
                    return Ok(None);
                };
                let Some(second) = factors.next() else {
                    return Ok(None);
                };
                if factors.next().is_some() {
                    return Ok(None);
                }

                let (coefficient, group) = match (first, second) {
                    (coefficient @ AtomView::Num(_), AtomView::Fun(group))
                    | (AtomView::Fun(group), coefficient @ AtomView::Num(_)) => {
                        (coefficient, group)
                    }
                    _ => return Ok(None),
                };
                if Rational::try_from(coefficient).ok() != Some(Rational::from(-1))
                    || group.get_symbol() != *ANTISYM
                    || group.get_nargs() == 0
                {
                    return Ok(None);
                }
                (Atom::num(-1), group)
            }
            _ => return Ok(None),
        };

        let symbol = group.get_symbol();
        let kind = SymmetryKind::structural(symbol).expect("structural group was classified");
        let mut slots = Vec::with_capacity(group.get_nargs());
        for (member, value) in group.iter().enumerate() {
            let Ok(slot) = Slot::<LibraryRep, Aind>::try_from(value) else {
                return Err(CanonicalizationError::InvalidPartialGroup {
                    group: symbol,
                    member,
                    expression: expression.clone(),
                });
            };
            slots.push(slot);
        }
        Ok(Some((symbol, kind, prefix, slots)))
    }

    fn incidence_role(&self, flat_slot: usize) -> IncidenceRole {
        if let Some(kind) = self.intrinsic {
            return IncidenceRole::Group {
                key: GroupKey::Intrinsic,
                kind,
            };
        }
        for (argument, layout) in self.arguments.iter().enumerate() {
            match layout {
                LayoutArgument::DirectSlot(hole) if *hole == flat_slot => {
                    return IncidenceRole::Ordered { argument };
                }
                LayoutArgument::SlotBundle { holes } => {
                    if let Some(member) = holes.iter().position(|hole| *hole == flat_slot) {
                        return IncidenceRole::OrderedBundle { argument, member };
                    }
                }
                LayoutArgument::Group { kind, holes, .. } => {
                    if holes.contains(&flat_slot) {
                        return IncidenceRole::Group {
                            key: GroupKey::Argument(argument),
                            kind: *kind,
                        };
                    }
                }
                LayoutArgument::Opaque(_) | LayoutArgument::DirectSlot(_) => {}
            }
        }

        unreachable!("every flat slot belongs to one reversible layout hole")
    }

    fn member_position(&self, flat_slot: usize) -> usize {
        if self.intrinsic.is_some() {
            return self
                .arguments
                .iter()
                .position(
                    |layout| matches!(layout, LayoutArgument::DirectSlot(hole) if *hole == flat_slot),
                )
                .expect("every intrinsic slot is one direct argument");
        }
        for layout in &self.arguments {
            match layout {
                LayoutArgument::Group { holes, .. } | LayoutArgument::SlotBundle { holes } => {
                    if let Some(position) = holes.iter().position(|hole| *hole == flat_slot) {
                        return position;
                    }
                }
                LayoutArgument::DirectSlot(hole) if *hole == flat_slot => return 0,
                _ => {}
            }
        }
        flat_slot
    }

    fn antisymmetric_groups(&self) -> impl Iterator<Item = (GroupKey, bool)> + '_ {
        let intrinsic = (self.intrinsic == Some(SymmetryKind::Antisymmetric))
            .then_some((GroupKey::Intrinsic, true));
        intrinsic
            .into_iter()
            .chain(self.arguments.iter().enumerate().filter_map(
                |(argument, layout)| match layout {
                    LayoutArgument::Group {
                        kind: SymmetryKind::Antisymmetric,
                        ..
                    } => Some((GroupKey::Argument(argument), self.outer_linear)),
                    _ => None,
                },
            ))
    }

    fn sign_lifts(&self, key: GroupKey) -> bool {
        key == GroupKey::Intrinsic || self.outer_linear
    }

    fn rebuild<Aind: AbsInd + ParseableAind>(
        &self,
        slots: &[LibrarySlot<Aind>],
        zero_groups: &BTreeSet<GroupKey>,
        negative_groups: &BTreeSet<GroupKey>,
    ) -> (Atom, Vec<LibrarySlot<Aind>>) {
        let mut builder = FunctionBuilder::new(self.head);
        let mut active_slots = Vec::new();
        for (argument, layout) in self.arguments.iter().enumerate() {
            let value = match layout {
                LayoutArgument::Opaque(value) => value.clone(),
                LayoutArgument::DirectSlot(hole) => {
                    active_slots.push(slots[*hole]);
                    slots[*hole].to_atom()
                }
                LayoutArgument::SlotBundle { holes } => {
                    let mut bundle = FunctionBuilder::new(AIND_SYMBOLS.aind);
                    for hole in holes {
                        active_slots.push(slots[*hole]);
                        bundle = bundle.add_arg(slots[*hole].to_atom());
                    }
                    bundle.finish()
                }
                LayoutArgument::Group {
                    symbol,
                    prefix,
                    holes,
                    ..
                } => {
                    if zero_groups.contains(&GroupKey::Argument(argument)) {
                        Atom::Zero
                    } else {
                        active_slots.extend(holes.iter().map(|hole| slots[*hole]));
                        let mut group = FunctionBuilder::new(*symbol);
                        for hole in holes {
                            group = group.add_arg(slots[*hole].to_atom());
                        }
                        let group = group.finish();
                        match (
                            prefix == &Atom::num(1),
                            negative_groups.contains(&GroupKey::Argument(argument)),
                        ) {
                            (true, false) | (false, true) => group,
                            (true, true) => &Atom::num(-1) * group,
                            (false, false) => prefix * group,
                        }
                    }
                }
            };
            builder = builder.add_arg(value);
        }
        (builder.finish(), active_slots)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum GroupKey {
    Intrinsic,
    Argument(usize),
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum IncidenceRole {
    Ordered { argument: usize },
    OrderedBundle { argument: usize, member: usize },
    Group { key: GroupKey, kind: SymmetryKind },
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum SlotColor<Aind> {
    Boundary(Representation<LibraryRep>, usize),
    Internal(Representation<LibraryRep>),
    External(LibrarySlot<Aind>),
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum NodeColor<Aind> {
    Product,
    Tensor(TensorColor),
    Slot(SlotColor<Aind>),
}

#[derive(Clone, Debug)]
enum NodeOrigin<Aind: AbsInd> {
    None,
    Tensor(usize),
    Slot(LibrarySlot<Aind>),
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum EdgeColor {
    Factor,
    Incidence(IncidenceRole),
    Cyclic(GroupKey),
    Contraction(Representation<LibraryRep>),
}

#[derive(Clone, Debug)]
enum EdgeOrigin {
    Incidence { flat_slot: usize, member: usize },
    None,
}

type TensorGraph<Aind> =
    Graph<HiddenData<NodeColor<Aind>, NodeOrigin<Aind>>, HiddenData<EdgeColor, EdgeOrigin>>;

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum ContextNodeColor<Aind> {
    Sum,
    Neg,
    Product,
    Power(i8),
    Function(Symbol),
    Scalar(Atom),
    Tensor(TensorColor),
    Port(Representation<LibraryRep>),
    Line {
        group: Representation<LibraryRep>,
        boundary: bool,
    },
    External(LibrarySlot<Aind>),
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum ContextEdgeColor {
    Child,
    Incidence(IncidenceRole),
    Cyclic(GroupKey),
    Line,
}

#[derive(Clone)]
struct Occurrence<Aind: AbsInd> {
    header: usize,
    layout: TensorLayout,
    tensor: SymbolicTensor<Aind>,
}

#[derive(Clone)]
struct GroupMember<Aind: AbsInd> {
    flat_slot: usize,
    member: usize,
    vertex: usize,
    slot: LibrarySlot<Aind>,
    atom: Atom,
}

struct SlotCopies<Aind: AbsInd> {
    slot: LibrarySlot<Aind>,
    power_path: Vec<NodeIndex>,
    vertices: Vec<usize>,
}

struct Projection<Aind: AbsInd> {
    graph: TensorGraph<Aind>,
    occurrences: Vec<Occurrence<Aind>>,
    slot_copies: BTreeMap<Hedge, SlotCopies<Aind>>,
    scalar: Atom,
}

impl<Aind: AbsInd> Default for Projection<Aind> {
    fn default() -> Self {
        Self {
            graph: Graph::new(),
            occurrences: Vec::new(),
            slot_copies: BTreeMap::new(),
            scalar: Atom::num(1),
        }
    }
}

struct Canonicalizer<'a, Aind: AbsInd, F> {
    network: &'a SymbolicNet<Aind>,
    tree: SimpleTraversalTree<ParentChildStore<()>>,
    external_hedges: BTreeSet<Hedge>,
    boundary_ids: BTreeMap<Hedge, usize>,
    boundary_groups: BTreeMap<usize, Representation<LibraryRep>>,
    boundary_order: Vec<usize>,
    boundary_dummies: BTreeMap<usize, Aind>,
    dummy_pools: BTreeMap<Representation<LibraryRep>, Vec<Aind>>,
    next_dummy: BTreeMap<Representation<LibraryRep>, usize>,
    new_dummy: &'a mut F,
}

struct CanonicalSubnetwork<Aind: AbsInd> {
    network: SymbolicNet<Aind>,
    zero: bool,
}

impl<'a, Aind, F> Canonicalizer<'a, Aind, F>
where
    Aind: AbsInd + DummyAind + ParseableAind,
    F: FnMut(usize) -> Aind,
{
    fn new(
        network: &'a SymbolicNet<Aind>,
        new_dummy: &'a mut F,
    ) -> Result<Self, CanonicalizationError> {
        let tree: SimpleTraversalTree<ParentChildStore<()>> = network.graph.expr_tree().cast();
        let root = network.graph.graph.node_id(network.graph.head());
        let nodes = tree
            .iter_preorder_tree_nodes(network.graph.graph.as_ref(), root)
            .collect::<Vec<_>>();

        let mut tensor_scopes = BTreeMap::new();
        let mut pending = vec![root];
        while let Some(scope) = pending.pop() {
            let has_barrier = tree
                .iter_preorder_tree_nodes(network.graph.graph.as_ref(), scope)
                .any(|node| {
                    matches!(
                        network.graph.graph[node],
                        NetworkNode::Op(NetworkOp::Sum | NetworkOp::Function(_))
                    )
                });
            if has_barrier {
                pending.extend(tree.iter_children(scope, network.graph.graph.as_ref()));
            } else {
                for node in tree.iter_preorder_tree_nodes(network.graph.graph.as_ref(), scope) {
                    if matches!(
                        network.graph.graph[node],
                        NetworkNode::Leaf(NetworkLeaf::LocalTensor(_))
                    ) {
                        tensor_scopes.insert(node, scope);
                    }
                }
            }
        }

        let mut components = UnionFind::new(vec![(); network.graph.graph.n_hedges()]);
        for (pair, _, data) in network.graph.graph.iter_edges() {
            if !matches!(data.data, NetworkEdge::Slot(_)) {
                continue;
            }
            if let HedgePair::Paired { source, sink } = pair {
                components.union(source, sink, |(), ()| ());
            }
        }
        for &node in &nodes {
            if !matches!(network.graph.graph[node], NetworkNode::Op(_)) {
                continue;
            }
            let mut slots = BTreeMap::new();
            for hedge in network.graph.graph.iter_crown(node) {
                let NetworkEdge::Slot(slot) = network.graph.graph[[&hedge]] else {
                    continue;
                };
                if let Some(previous) = slots.insert(slot, hedge) {
                    components.union(previous, hedge, |(), ()| ());
                }
            }
        }

        let mut component_slots = BTreeMap::<
            Hedge,
            (
                Representation<LibraryRep>,
                Option<LibrarySlot<Aind>>,
                BTreeSet<NodeIndex>,
            ),
        >::new();
        for hedge in (0..network.graph.graph.n_hedges()).map(Hedge) {
            let NetworkEdge::Slot(slot) = network.graph.graph[[&hedge]] else {
                continue;
            };
            let component = components.find(hedge);
            let entry = component_slots.entry(component).or_insert((
                slot.rep().base(),
                None,
                BTreeSet::new(),
            ));
            if network.graph.graph.inv(hedge) == hedge {
                entry.1 = Some(slot);
            }
            let node = network.graph.graph.node_id(hedge);
            if let Some(scope) = tensor_scopes.get(&node) {
                entry.2.insert(*scope);
            }
        }

        let mut context = Graph::<ContextNodeColor<Aind>, ContextEdgeColor>::new();
        let mut context_nodes = BTreeMap::new();
        for &node in &nodes {
            let color = match &network.graph.graph[node] {
                NetworkNode::Op(NetworkOp::Sum) => ContextNodeColor::Sum,
                NetworkNode::Op(NetworkOp::Neg) => ContextNodeColor::Neg,
                NetworkNode::Op(NetworkOp::Product) => ContextNodeColor::Product,
                NetworkNode::Op(NetworkOp::Power(power)) => ContextNodeColor::Power(*power),
                NetworkNode::Op(NetworkOp::Function(function)) => {
                    if function.has_tag(&SPENSO_TAG.tensor) && SymmetryKind::of(*function).is_some()
                    {
                        return Err(CanonicalizationError::InvalidIntrinsicArgument {
                            head: *function,
                            argument: 0,
                            expression: FunctionBuilder::new(*function)
                                .add_arg(Atom::Zero)
                                .finish(),
                        });
                    }
                    ContextNodeColor::Function(*function)
                }
                NetworkNode::Leaf(NetworkLeaf::Scalar(reference)) => {
                    let scalar = network.store.get_scalar_ref(*reference);
                    validate_tensor_symmetry::<Aind>(scalar.as_view())?;
                    ContextNodeColor::Scalar(scalar.clone())
                }
                NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => ContextNodeColor::Tensor(
                    TensorLayout::scan(&network.store.tensors[*index])?.color,
                ),
                NetworkNode::Leaf(_) => return Err(CanonicalizationError::UnsupportedLeaf),
            };
            context_nodes.insert(node, context.add_node(color));
        }
        for &node in &nodes {
            for child in tree.iter_children(node, network.graph.graph.as_ref()) {
                context
                    .add_edge(
                        context_nodes[&node],
                        context_nodes[&child],
                        true,
                        ContextEdgeColor::Child,
                    )
                    .expect("an expression tree has one edge per child");
            }
        }

        let mut context_lines = BTreeMap::new();
        for (&component, (group, external, scopes)) in &component_slots {
            let color = external.map_or(
                ContextNodeColor::Line {
                    group: *group,
                    boundary: scopes.len() > 1,
                },
                ContextNodeColor::External,
            );
            context_lines.insert(component, context.add_node(color));
        }

        for &node in &nodes {
            let NetworkNode::Leaf(NetworkLeaf::LocalTensor(tensor_index)) =
                network.graph.graph[node]
            else {
                continue;
            };
            let tensor = &network.store.tensors[tensor_index];
            let layout = TensorLayout::scan(tensor)?;
            let mut hedges = network
                .graph
                .graph
                .iter_crown(node)
                .filter(|hedge| matches!(network.graph.graph[[hedge]], NetworkEdge::Slot(_)))
                .collect::<Vec<_>>();
            hedges.sort_unstable_by_key(|hedge| network.graph.slot_order[hedge.0]);
            if hedges.len() != layout.slot_count {
                return Err(CanonicalizationError::StructureMismatch {
                    expression: tensor.expression.clone(),
                });
            }
            let mut cycle_groups = BTreeMap::<GroupKey, Vec<(usize, usize)>>::new();
            for (flat_slot, hedge) in hedges.into_iter().enumerate() {
                let NetworkEdge::Slot(slot) = network.graph.graph[[&hedge]] else {
                    unreachable!("tensor slot hedges carry slot data")
                };
                let port = context.add_node(ContextNodeColor::Port(slot.rep()));
                let role = layout.incidence_role(flat_slot);
                context
                    .add_edge(
                        context_nodes[&node],
                        port,
                        true,
                        ContextEdgeColor::Incidence(role.clone()),
                    )
                    .expect("a tensor has one context edge per slot");
                context
                    .add_edge(
                        port,
                        context_lines[&components.find(hedge)],
                        false,
                        ContextEdgeColor::Line,
                    )
                    .expect("a tensor slot belongs to one context line");
                if let IncidenceRole::Group {
                    key,
                    kind: SymmetryKind::Cyclic,
                } = role
                {
                    cycle_groups
                        .entry(key)
                        .or_default()
                        .push((layout.member_position(flat_slot), port));
                }
            }
            for (key, mut ports) in cycle_groups {
                ports.sort_unstable_by_key(|(member, _)| *member);
                for (&(_, source), &(_, target)) in ports
                    .iter()
                    .zip(ports.iter().cycle().skip(1))
                    .take(ports.len())
                {
                    context
                        .add_edge(source, target, true, ContextEdgeColor::Cyclic(key))
                        .expect("cyclic context markers are unique");
                }
            }
        }

        let canonical_context = context.canonize();
        let mut boundary_components = component_slots
            .iter()
            .filter(|(_, (_, external, scopes))| external.is_none() && scopes.len() > 1)
            .map(|(&component, (group, _, _))| {
                (
                    canonical_context.vertex_map[context_lines[&component]],
                    component,
                    *group,
                )
            })
            .collect::<Vec<_>>();
        boundary_components.sort_by_key(|(canonical, _, _)| *canonical);

        let mut boundary_order = Vec::new();
        let mut seen_boundaries = BTreeSet::new();
        let mut stack = vec![canonical_context.vertex_map[context_nodes[&root]]];
        while let Some(node) = stack.pop() {
            match &canonical_context.graph.node(node).data {
                ContextNodeColor::Tensor(_) => {
                    for edge_index in &canonical_context.graph.node(node).edges {
                        let edge = canonical_context.graph.edge(*edge_index);
                        if edge.vertices.0 != node
                            || !matches!(edge.data, ContextEdgeColor::Incidence(_))
                        {
                            continue;
                        }
                        let port = edge.vertices.1;
                        let line = canonical_context.graph.node(port).edges.iter().find_map(
                            |line_edge_index| {
                                let line_edge = canonical_context.graph.edge(*line_edge_index);
                                matches!(line_edge.data, ContextEdgeColor::Line).then(|| {
                                    if line_edge.vertices.0 == port {
                                        line_edge.vertices.1
                                    } else {
                                        line_edge.vertices.0
                                    }
                                })
                            },
                        );
                        if let Some(line) = line
                            && matches!(
                                canonical_context.graph.node(line).data,
                                ContextNodeColor::Line { boundary: true, .. }
                            )
                            && seen_boundaries.insert(line)
                        {
                            boundary_order.push(line);
                        }
                    }
                }
                _ => {
                    let mut children = canonical_context
                        .graph
                        .node(node)
                        .edges
                        .iter()
                        .filter_map(|edge_index| {
                            let edge = canonical_context.graph.edge(*edge_index);
                            (edge.vertices.0 == node
                                && matches!(edge.data, ContextEdgeColor::Child))
                            .then_some(edge.vertices.1)
                        })
                        .collect::<Vec<_>>();
                    children.reverse();
                    stack.extend(children);
                }
            }
        }
        for (boundary, _, _) in &boundary_components {
            if seen_boundaries.insert(*boundary) {
                boundary_order.push(*boundary);
            }
        }

        let mut boundary_ids = BTreeMap::new();
        let mut boundary_groups = BTreeMap::new();
        for (canonical, component, group) in &boundary_components {
            boundary_groups.insert(*canonical, *group);
            for hedge in (0..network.graph.graph.n_hedges()).map(Hedge) {
                if matches!(network.graph.graph[[&hedge]], NetworkEdge::Slot(_))
                    && components.find(hedge) == *component
                {
                    boundary_ids.insert(hedge, *canonical);
                }
            }
        }
        let external_components = component_slots
            .iter()
            .filter_map(|(&component, (_, external, _))| external.map(|_| component))
            .collect::<BTreeSet<_>>();
        let external_hedges = (0..network.graph.graph.n_hedges())
            .map(Hedge)
            .filter(|hedge| {
                matches!(network.graph.graph[[hedge]], NetworkEdge::Slot(_))
                    && external_components.contains(&components.find(*hedge))
            })
            .collect();

        Ok(Self {
            tree,
            external_hedges,
            boundary_ids,
            boundary_groups,
            boundary_order,
            boundary_dummies: BTreeMap::new(),
            dummy_pools: BTreeMap::new(),
            next_dummy: BTreeMap::new(),
            network,
            new_dummy,
        })
    }

    fn run(&mut self) -> Result<SymbolicNet<Aind>, CanonicalizationError> {
        let root = self.network.graph.graph.node_id(self.network.graph.head());
        let mut capacities = self.dummy_usage(root)?;
        let mut boundary_counts = BTreeMap::new();
        for group in self.boundary_groups.values() {
            *boundary_counts.entry(*group).or_insert(0) += 1;
            *capacities.entry(*group).or_insert(0) += 1;
        }
        let mut groups = capacities.into_iter().collect::<Vec<_>>();
        groups.sort_by(|(left, _), (right, _)| {
            left.to_symbolic([Atom::num(0)])
                .as_view()
                .cmp(&right.to_symbolic([Atom::num(0)]).as_view())
        });
        let mut position = 0;
        for (group, capacity) in groups {
            let mut dummies = Vec::with_capacity(capacity);
            for _ in 0..capacity {
                dummies.push((self.new_dummy)(position));
                position += 1;
            }
            self.dummy_pools.insert(group, dummies);
        }
        let mut boundary_offsets = BTreeMap::new();
        for boundary in &self.boundary_order {
            let group = self.boundary_groups[boundary];
            let offset = boundary_offsets.entry(group).or_insert(0);
            self.boundary_dummies
                .insert(*boundary, self.dummy_pools[&group][*offset]);
            *offset += 1;
        }
        self.next_dummy = boundary_counts;
        Ok(self.canonize_node(root)?.network)
    }

    fn canonize_node(
        &mut self,
        node: NodeIndex,
    ) -> Result<CanonicalSubnetwork<Aind>, CanonicalizationError> {
        if !self.has_barrier(node) {
            return self.canonize_flat(node);
        }

        match &self.network.graph.graph[node] {
            NetworkNode::Op(NetworkOp::Sum) => {
                let children = self.children(node);
                let mut terms = Vec::new();
                let start_dummy = self.next_dummy.clone();
                let mut end_dummy = start_dummy.clone();
                for child in children {
                    self.next_dummy.clone_from(&start_dummy);
                    let term = self.canonize_node(child)?;
                    for (group, next) in &self.next_dummy {
                        let end = end_dummy.entry(*group).or_insert(0);
                        *end = (*end).max(*next);
                    }
                    if !term.zero {
                        terms.push(term.network);
                    }
                }
                self.next_dummy = end_dummy;
                let Some(first) = terms.pop() else {
                    return Ok(CanonicalSubnetwork {
                        network: Network::zero(),
                        zero: true,
                    });
                };
                Ok(CanonicalSubnetwork {
                    network: terms.into_iter().fold(first, |sum, term| sum + term),
                    zero: false,
                })
            }
            NetworkNode::Op(NetworkOp::Function(function)) => {
                let child = self.only_child(node)?;
                let inner = self.canonize_node(child)?;
                if inner.zero && function.is_linear() {
                    Ok(inner)
                } else {
                    Ok(CanonicalSubnetwork {
                        network: inner.network.fun(*function),
                        zero: false,
                    })
                }
            }
            NetworkNode::Op(NetworkOp::Product) => {
                let children = self.children(node);
                let mut factors = Vec::new();
                for child in children {
                    let factor = self.canonize_node(child)?;
                    if factor.zero {
                        return Ok(CanonicalSubnetwork {
                            network: Network::zero(),
                            zero: true,
                        });
                    }
                    factors.push(factor.network);
                }
                let Some(first) = factors.pop() else {
                    return Ok(CanonicalSubnetwork {
                        network: Network::one(),
                        zero: false,
                    });
                };
                Ok(CanonicalSubnetwork {
                    network: factors
                        .into_iter()
                        .fold(first, |product, factor| product * factor),
                    zero: false,
                })
            }
            NetworkNode::Op(NetworkOp::Neg) => {
                let child = self.canonize_node(self.only_child(node)?)?;
                if child.zero {
                    Ok(child)
                } else {
                    Ok(CanonicalSubnetwork {
                        network: -child.network,
                        zero: false,
                    })
                }
            }
            NetworkNode::Op(NetworkOp::Power(power)) => {
                if *power < 0 {
                    return Err(CanonicalizationError::UnsupportedPower { power: *power });
                }
                if *power == 0 {
                    return Ok(CanonicalSubnetwork {
                        network: Network::one(),
                        zero: false,
                    });
                }
                let base = self.canonize_node(self.only_child(node)?)?;
                if base.zero {
                    Ok(base)
                } else {
                    Ok(CanonicalSubnetwork {
                        network: base.network.pow(*power),
                        zero: false,
                    })
                }
            }
            _ => self.canonize_flat(node),
        }
    }

    fn children(&self, node: NodeIndex) -> Vec<NodeIndex> {
        self.tree
            .iter_children(node, self.network.graph.graph.as_ref())
            .collect()
    }

    fn only_child(&self, node: NodeIndex) -> Result<NodeIndex, CanonicalizationError> {
        let children = self.children(node);
        match children.as_slice() {
            [child] => Ok(*child),
            _ => Err(CanonicalizationError::UnsupportedLeaf),
        }
    }

    fn has_barrier(&self, root: NodeIndex) -> bool {
        self.tree
            .iter_preorder_tree_nodes(self.network.graph.graph.as_ref(), root)
            .any(|node| {
                matches!(
                    self.network.graph.graph[node],
                    NetworkNode::Op(NetworkOp::Sum | NetworkOp::Function(_))
                )
            })
    }

    fn dummy_usage(
        &self,
        root: NodeIndex,
    ) -> Result<BTreeMap<Representation<LibraryRep>, usize>, CanonicalizationError> {
        if !self.has_barrier(root) {
            let mut projection = Projection::default();
            self.project_expression(root, &mut Vec::new(), &mut projection)?;
            self.connect_slots(&mut projection)?;
            let mut usage = BTreeMap::new();
            for edge in projection.graph.edges() {
                if let EdgeColor::Contraction(group) = edge.data.data {
                    *usage.entry(group).or_insert(0) += 1;
                }
            }
            return Ok(usage);
        }

        let mut usage = BTreeMap::new();
        for child in self.children(root) {
            for (group, count) in self.dummy_usage(child)? {
                let current = usage.entry(group).or_insert(0);
                if matches!(
                    self.network.graph.graph[root],
                    NetworkNode::Op(NetworkOp::Sum)
                ) {
                    *current = (*current).max(count);
                } else {
                    *current += count;
                }
            }
        }
        Ok(usage)
    }

    fn canonize_flat(
        &mut self,
        root: NodeIndex,
    ) -> Result<CanonicalSubnetwork<Aind>, CanonicalizationError> {
        let mut projection = Projection::default();
        self.project_expression(root, &mut Vec::new(), &mut projection)?;
        self.connect_slots(&mut projection)?;

        if projection.occurrences.is_empty() {
            let zero = projection.scalar.as_view().is_zero();
            return Ok(CanonicalSubnetwork {
                network: Network::from_scalar(projection.scalar),
                zero,
            });
        }

        if projection.occurrences.len() > 1 {
            let product = projection
                .graph
                .add_node(HiddenData::new(NodeColor::Product, NodeOrigin::None));
            for occurrence in &projection.occurrences {
                projection
                    .graph
                    .add_edge(
                        product,
                        occurrence.header,
                        true,
                        HiddenData::new(EdgeColor::Factor, EdgeOrigin::None),
                    )
                    .expect("a product has one edge per tensor factor");
            }
        }

        let canonical = projection.graph.canonize();
        let zero_groups = self.stabilizer_zeros(
            &canonical.graph,
            &canonical.orbit_generators,
            &projection.occurrences,
        )?;
        if zero_groups.is_none() {
            return Ok(CanonicalSubnetwork {
                network: Network::zero(),
                zero: true,
            });
        }
        let zero_groups = zero_groups.expect("global zero was returned above");

        let assigned = self.assign_slots(&canonical.graph);
        let mut factors = Vec::new();
        let mut global_negative = false;
        if projection.scalar != Atom::num(1) {
            factors.push(Network::from_scalar(projection.scalar));
        }

        for node_index in 0..canonical.graph.nodes().len() {
            let node = canonical.graph.node(node_index);
            let (NodeColor::Tensor(_), NodeOrigin::Tensor(origin)) =
                (&node.data.data, &node.data.hidden)
            else {
                continue;
            };
            let occurrence = &projection.occurrences[*origin];
            let mut slots = vec![None; occurrence.layout.slot_count];
            let mut members = BTreeMap::<(GroupKey, SymmetryKind), Vec<GroupMember<Aind>>>::new();
            for edge_index in &node.edges {
                let edge = canonical.graph.edge(*edge_index);
                if edge.vertices.0 != node_index {
                    continue;
                }
                let EdgeOrigin::Incidence { flat_slot, member } = edge.data.hidden else {
                    continue;
                };
                let Some(slot) = assigned[edge.vertices.1] else {
                    continue;
                };
                if let EdgeColor::Incidence(IncidenceRole::Group { key, kind }) = edge.data.data {
                    members.entry((key, kind)).or_default().push(GroupMember {
                        flat_slot,
                        member,
                        vertex: edge.vertices.1,
                        slot,
                        atom: slot.to_atom(),
                    });
                } else {
                    slots[flat_slot] = Some(slot);
                }
            }

            let mut negative_groups = BTreeSet::new();
            for ((key, kind), group) in members {
                let holes = group
                    .iter()
                    .map(|member| (member.member, member.flat_slot))
                    .collect::<BTreeMap<_, _>>();
                let ordered = Self::canonical_group_members(&canonical.graph, key, kind, group)?;
                if kind == SymmetryKind::Antisymmetric
                    && !zero_groups
                        .get(origin)
                        .is_some_and(|groups| groups.contains(&key))
                    && Permutation::from_map(ordered.iter().map(|member| member.member).collect())
                        .sign()
                        < 0
                {
                    if occurrence.layout.sign_lifts(key) {
                        global_negative ^= true;
                    } else {
                        negative_groups.insert(key);
                    }
                }
                for (canonical_member, member) in ordered.into_iter().enumerate() {
                    let hole = holes[&canonical_member];
                    slots[hole] = Some(member.slot);
                }
            }
            let slots = slots
                .into_iter()
                .collect::<Option<Vec<_>>>()
                .ok_or_else(|| CanonicalizationError::StructureMismatch {
                    expression: occurrence.tensor.expression.clone(),
                })?;
            let groups = zero_groups.get(origin).cloned().unwrap_or_default();
            let (expression, active_slots) =
                occurrence.layout.rebuild(&slots, &groups, &negative_groups);
            if active_slots.is_empty() {
                factors.push(Network::from_scalar(expression));
            } else {
                factors.push(Network::from_tensor(SymbolicTensor {
                    structure: OrderedStructure::new(active_slots).structure,
                    is_metric: occurrence.tensor.is_metric,
                    is_composite: occurrence.tensor.is_composite,
                    expression,
                }));
            }
        }

        if global_negative {
            factors.push(Network::from_scalar(Atom::num(-1)));
        }

        let Some(first) = factors.pop() else {
            return Ok(CanonicalSubnetwork {
                network: Network::one(),
                zero: false,
            });
        };
        Ok(CanonicalSubnetwork {
            network: factors
                .into_iter()
                .fold(first, |product, factor| product * factor),
            zero: false,
        })
    }

    fn project_expression(
        &self,
        node: NodeIndex,
        power_path: &mut Vec<NodeIndex>,
        projection: &mut Projection<Aind>,
    ) -> Result<(), CanonicalizationError> {
        match &self.network.graph.graph[node] {
            NetworkNode::Op(NetworkOp::Product) => {
                for child in self.children(node) {
                    self.project_expression(child, power_path, projection)?;
                }
            }
            NetworkNode::Op(NetworkOp::Power(power)) if *power > 0 => {
                let child = self.only_child(node)?;
                power_path.push(node);
                for _ in 0..*power {
                    self.project_expression(child, power_path, projection)?;
                }
                power_path.pop();
            }
            NetworkNode::Op(NetworkOp::Power(0)) => {}
            NetworkNode::Op(NetworkOp::Power(power)) => {
                return Err(CanonicalizationError::UnsupportedPower { power: *power });
            }
            NetworkNode::Op(NetworkOp::Neg) => {
                projection.scalar = -projection.scalar.clone();
                self.project_expression(self.only_child(node)?, power_path, projection)?;
            }
            NetworkNode::Leaf(NetworkLeaf::Scalar(reference)) => {
                projection.scalar *= self.scalar(*reference).clone();
            }
            NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => {
                self.add_tensor(node, *index, power_path, projection)?;
            }
            NetworkNode::Op(NetworkOp::Sum | NetworkOp::Function(_)) => {
                return Err(CanonicalizationError::UnsupportedLeaf);
            }
            NetworkNode::Leaf(_) => return Err(CanonicalizationError::UnsupportedLeaf),
        }
        Ok(())
    }

    fn scalar(&self, reference: ScalarRef) -> &Atom {
        self.network.store.get_scalar_ref(reference)
    }

    fn add_tensor(
        &self,
        network_node: NodeIndex,
        tensor_index: usize,
        power_path: &[NodeIndex],
        projection: &mut Projection<Aind>,
    ) -> Result<(), CanonicalizationError> {
        let tensor = self.network.store.tensors[tensor_index].clone();
        let layout = TensorLayout::scan(&tensor)?;
        let structural_slots = tensor
            .structure
            .external_structure_iter()
            .collect::<Vec<_>>();
        let mut network_slots = self
            .network
            .graph
            .graph
            .iter_crown(network_node)
            .filter(|hedge| matches!(self.network.graph.graph[[hedge]], NetworkEdge::Slot(_)))
            .collect::<Vec<_>>();
        network_slots.sort_unstable_by_key(|hedge| self.network.graph.slot_order[hedge.0]);
        if structural_slots.len() != network_slots.len()
            || structural_slots.len() != layout.slot_count
        {
            return Err(CanonicalizationError::StructureMismatch {
                expression: tensor.expression,
            });
        }

        let occurrence_index = projection.occurrences.len();
        let header = projection.graph.add_node(HiddenData::new(
            NodeColor::Tensor(layout.color.clone()),
            NodeOrigin::Tensor(occurrence_index),
        ));
        let mut cycle_groups = BTreeMap::<GroupKey, Vec<(usize, usize)>>::new();

        for (flat_slot, (slot, hedge)) in structural_slots
            .iter()
            .copied()
            .zip(network_slots)
            .enumerate()
        {
            let color = if self.external_hedges.contains(&hedge) {
                SlotColor::External(slot)
            } else if let Some(boundary) = self.boundary_ids.get(&hedge) {
                SlotColor::Boundary(slot.rep(), *boundary)
            } else {
                SlotColor::Internal(slot.rep())
            };
            let vertex = projection.graph.add_node(HiddenData::new(
                NodeColor::Slot(color),
                NodeOrigin::Slot(slot),
            ));
            let role = layout.incidence_role(flat_slot);
            let member = layout.member_position(flat_slot);
            projection
                .graph
                .add_edge(
                    header,
                    vertex,
                    true,
                    HiddenData::new(
                        EdgeColor::Incidence(role.clone()),
                        EdgeOrigin::Incidence { flat_slot, member },
                    ),
                )
                .expect("tensor incidence cannot be a duplicate edge");
            if let IncidenceRole::Group {
                key,
                kind: SymmetryKind::Cyclic,
            } = role
            {
                cycle_groups.entry(key).or_default().push((member, vertex));
            }
            projection
                .slot_copies
                .entry(hedge)
                .or_insert_with(|| SlotCopies {
                    slot,
                    power_path: power_path.to_vec(),
                    vertices: Vec::new(),
                })
                .vertices
                .push(vertex);
        }

        for (key, mut vertices) in cycle_groups {
            vertices.sort_unstable_by_key(|(member, _)| *member);
            if vertices.len() > 1 {
                for (&(_, source), &(_, target)) in vertices
                    .iter()
                    .zip(vertices.iter().cycle().skip(1))
                    .take(vertices.len())
                {
                    projection
                        .graph
                        .add_edge(
                            source,
                            target,
                            true,
                            HiddenData::new(EdgeColor::Cyclic(key), EdgeOrigin::None),
                        )
                        .expect("cyclic marker cannot be a duplicate edge");
                }
            }
        }

        projection.occurrences.push(Occurrence {
            header,
            layout,
            tensor,
        });
        Ok(())
    }

    fn connect_slots(
        &self,
        projection: &mut Projection<Aind>,
    ) -> Result<(), CanonicalizationError> {
        for (pair, _, data) in self.network.graph.graph.iter_edges() {
            let NetworkEdge::Slot(slot) = data.data else {
                continue;
            };
            let HedgePair::Paired { source, sink } = pair else {
                continue;
            };
            let source_copies = projection.slot_copies.get(&source);
            let sink_copies = projection.slot_copies.get(&sink);
            let group = slot.rep().base();

            match (source_copies, sink_copies) {
                (Some(source), Some(sink))
                    if source.vertices.len() == sink.vertices.len()
                        && (source.vertices.len() == 1 || source.power_path == sink.power_path) =>
                {
                    let pairs = source
                        .vertices
                        .iter()
                        .copied()
                        .zip(sink.vertices.iter().copied())
                        .map(|(source_vertex, sink_vertex)| {
                            (
                                source_vertex,
                                sink_vertex,
                                source.slot.rep(),
                                sink.slot.rep(),
                            )
                        })
                        .collect::<Vec<_>>();
                    for (source_vertex, sink_vertex, source_rep, sink_rep) in pairs {
                        Self::connect_slot_pair(
                            &mut projection.graph,
                            source_vertex,
                            sink_vertex,
                            source_rep,
                            sink_rep,
                            group,
                        );
                    }
                }
                (Some(source), Some(sink)) => {
                    let source_remainder =
                        Self::connect_power_copies(&mut projection.graph, source, group);
                    let sink_remainder =
                        Self::connect_power_copies(&mut projection.graph, sink, group);
                    let (Some(source_vertex), Some(sink_vertex)) =
                        (source_remainder, sink_remainder)
                    else {
                        return Err(CanonicalizationError::UnsupportedLeaf);
                    };
                    Self::connect_slot_pair(
                        &mut projection.graph,
                        source_vertex,
                        sink_vertex,
                        source.slot.rep(),
                        sink.slot.rep(),
                        group,
                    );
                }
                (Some(copies), None) if self.is_power_node(sink) => {
                    let _ = Self::connect_power_copies(&mut projection.graph, copies, group);
                }
                (None, Some(copies)) if self.is_power_node(source) => {
                    let _ = Self::connect_power_copies(&mut projection.graph, copies, group);
                }
                (Some(_), None) | (None, Some(_)) | (None, None) => {}
            }
        }

        let dangling_powers = projection
            .slot_copies
            .values()
            .filter(|copies| copies.vertices.len() > 1)
            .filter(|copies| {
                copies.vertices.iter().all(|vertex| {
                    matches!(
                        projection.graph.node(*vertex).data.data,
                        NodeColor::Slot(SlotColor::External(_))
                    )
                })
            })
            .map(|copies| {
                (
                    copies.vertices.clone(),
                    copies.slot.rep(),
                    copies.slot.rep().base(),
                )
            })
            .collect::<Vec<_>>();
        for (vertices, representation, group) in dangling_powers {
            for pair in vertices.chunks_exact(2) {
                Self::connect_slot_pair(
                    &mut projection.graph,
                    pair[0],
                    pair[1],
                    representation,
                    representation,
                    group,
                );
            }
        }
        Ok(())
    }

    fn connect_power_copies(
        graph: &mut TensorGraph<Aind>,
        copies: &SlotCopies<Aind>,
        group: Representation<LibraryRep>,
    ) -> Option<usize> {
        for pair in copies.vertices.chunks_exact(2) {
            Self::connect_slot_pair(
                graph,
                pair[0],
                pair[1],
                copies.slot.rep(),
                copies.slot.rep(),
                group,
            );
        }
        (copies.vertices.len() % 2 == 1).then(|| *copies.vertices.last().unwrap())
    }

    fn is_power_node(&self, hedge: Hedge) -> bool {
        matches!(
            self.network.graph.graph[self.network.graph.graph.node_id(hedge)],
            NetworkNode::Op(NetworkOp::Power(_))
        )
    }

    fn connect_slot_pair(
        graph: &mut TensorGraph<Aind>,
        source: usize,
        sink: usize,
        source_rep: Representation<LibraryRep>,
        sink_rep: Representation<LibraryRep>,
        group: Representation<LibraryRep>,
    ) {
        let source_origin = graph.node(source).data.hidden.clone();
        let sink_origin = graph.node(sink).data.hidden.clone();
        graph.set_node_data(
            source,
            HiddenData::new(
                NodeColor::Slot(SlotColor::Internal(source_rep)),
                source_origin,
            ),
        );
        graph.set_node_data(
            sink,
            HiddenData::new(NodeColor::Slot(SlotColor::Internal(sink_rep)), sink_origin),
        );
        graph
            .add_edge(
                source,
                sink,
                false,
                HiddenData::new(EdgeColor::Contraction(group), EdgeOrigin::None),
            )
            .expect("network contraction cannot be a duplicate edge");
    }

    fn canonical_group_members(
        graph: &TensorGraph<Aind>,
        key: GroupKey,
        kind: SymmetryKind,
        mut members: Vec<GroupMember<Aind>>,
    ) -> Result<Vec<GroupMember<Aind>>, CanonicalizationError> {
        if members.len() < 2 {
            return Ok(members);
        }
        if kind != SymmetryKind::Cyclic {
            members.sort_by(|left, right| {
                left.atom
                    .as_view()
                    .cmp(&right.atom.as_view())
                    .then(left.vertex.cmp(&right.vertex))
            });
            return Ok(members);
        }

        let mut by_vertex = members
            .into_iter()
            .map(|member| (member.vertex, member))
            .collect::<BTreeMap<_, _>>();
        let Some(mut vertex) = by_vertex.keys().next().copied() else {
            return Ok(Vec::new());
        };
        let count = by_vertex.len();
        let mut cycle = Vec::with_capacity(count);
        for _ in 0..count {
            let member = by_vertex
                .remove(&vertex)
                .ok_or(CanonicalizationError::UnsupportedLeaf)?;
            cycle.push(member);
            vertex = graph
                .node(vertex)
                .edges
                .iter()
                .map(|edge| graph.edge(*edge))
                .find_map(|edge| {
                    (edge.vertices.0 == vertex
                        && matches!(edge.data.data, EdgeColor::Cyclic(edge_key) if edge_key == key))
                    .then_some(edge.vertices.1)
                })
                .ok_or(CanonicalizationError::UnsupportedLeaf)?;
        }

        let mut best = 0;
        for candidate in 1..cycle.len() {
            for offset in 0..cycle.len() {
                match cycle[(candidate + offset) % cycle.len()]
                    .atom
                    .as_view()
                    .cmp(&cycle[(best + offset) % cycle.len()].atom.as_view())
                {
                    std::cmp::Ordering::Less => {
                        best = candidate;
                        break;
                    }
                    std::cmp::Ordering::Greater => break,
                    std::cmp::Ordering::Equal => {}
                }
            }
        }
        cycle.rotate_left(best);
        Ok(cycle)
    }

    fn assign_slots(&mut self, graph: &TensorGraph<Aind>) -> Vec<Option<LibrarySlot<Aind>>> {
        let mut assigned = vec![None; graph.nodes().len()];
        for (node_index, node) in graph.nodes().iter().enumerate() {
            if let NodeColor::Slot(SlotColor::External(slot)) = node.data.data {
                assigned[node_index] = Some(slot);
            }
        }

        let product_headers = graph
            .nodes()
            .iter()
            .enumerate()
            .find_map(|(product, node)| {
                matches!(node.data.data, NodeColor::Product).then(|| {
                    node.edges
                        .iter()
                        .filter_map(|edge_index| {
                            let edge = graph.edge(*edge_index);
                            (edge.vertices.0 == product
                                && matches!(edge.data.data, EdgeColor::Factor))
                            .then_some(edge.vertices.1)
                        })
                        .collect::<Vec<_>>()
                })
            });
        let headers = product_headers.unwrap_or_else(|| {
            graph
                .nodes()
                .iter()
                .enumerate()
                .filter_map(|(node, data)| {
                    matches!(data.data.data, NodeColor::Tensor(_)).then_some(node)
                })
                .collect()
        });
        let canonical_slots = headers.into_iter().flat_map(|header| {
            graph
                .node(header)
                .edges
                .iter()
                .filter_map(move |edge_index| {
                    let edge = graph.edge(*edge_index);
                    (edge.vertices.0 == header && matches!(edge.data.data, EdgeColor::Incidence(_)))
                        .then_some(edge.vertices.1)
                })
        });

        for source in canonical_slots {
            if assigned[source].is_some() {
                continue;
            }
            if let NodeColor::Slot(SlotColor::Boundary(representation, boundary)) =
                graph.node(source).data.data
            {
                let dummy = self.boundary_dummies[&boundary];
                assigned[source] = Some(representation.slot(dummy));
                continue;
            }
            let Some(contraction) = graph.node(source).edges.iter().find_map(|edge_index| {
                let edge = graph.edge(*edge_index);
                matches!(edge.data.data, EdgeColor::Contraction(_)).then_some(edge)
            }) else {
                continue;
            };
            let sink = if contraction.vertices.0 == source {
                contraction.vertices.1
            } else {
                contraction.vertices.0
            };
            let EdgeColor::Contraction(group) = contraction.data.data else {
                unreachable!()
            };
            let next = self.next_dummy.entry(group).or_insert(0);
            let dummy = self.dummy_pools[&group][*next];
            *next += 1;
            let NodeColor::Slot(SlotColor::Internal(source_rep)) = graph.node(source).data.data
            else {
                unreachable!("a contraction endpoint is an internal slot")
            };
            let NodeColor::Slot(SlotColor::Internal(sink_rep)) = graph.node(sink).data.data else {
                unreachable!("a contraction endpoint is an internal slot")
            };
            assigned[source] = Some(source_rep.slot(dummy));
            assigned[sink] = Some(sink_rep.slot(dummy));
        }
        assigned
    }

    /// Return `None` for a global zero, otherwise the local group zeros keyed
    /// by the selected canonical occurrence origin.
    fn stabilizer_zeros(
        &self,
        graph: &TensorGraph<Aind>,
        generators: &[Vec<Vec<usize>>],
        occurrences: &[Occurrence<Aind>],
    ) -> Result<Option<BTreeMap<usize, BTreeSet<GroupKey>>>, CanonicalizationError> {
        let mut zeros = BTreeMap::<usize, BTreeSet<GroupKey>>::new();
        for cycles in generators {
            let mut automorphism = (0..graph.nodes().len()).collect::<Vec<_>>();
            for cycle in cycles {
                for (&source, &target) in cycle
                    .iter()
                    .zip(cycle.iter().cycle().skip(1))
                    .take(cycle.len())
                {
                    automorphism[source] = target;
                }
            }

            let mut global_odd = false;
            let mut local_odd = Vec::new();
            for header in 0..graph.nodes().len() {
                let (NodeColor::Tensor(_), NodeOrigin::Tensor(origin)) = (
                    &graph.node(header).data.data,
                    &graph.node(header).data.hidden,
                ) else {
                    continue;
                };
                let occurrence = &occurrences[*origin];
                for (key, lifts) in occurrence.layout.antisymmetric_groups() {
                    if !Self::group_is_odd(graph, &automorphism, header, key) {
                        continue;
                    }
                    if lifts {
                        global_odd ^= true;
                    } else if automorphism[header] == header
                        && Self::group_is_local(graph, header, key)
                    {
                        local_odd.push((*origin, key));
                    } else {
                        return Err(CanonicalizationError::AmbiguousSignScope);
                    }
                }
            }

            if local_odd.is_empty() {
                if global_odd {
                    return Ok(None);
                }
            } else {
                for (origin, key) in local_odd {
                    zeros.entry(origin).or_default().insert(key);
                }
            }
        }
        Ok(Some(zeros))
    }

    fn group_is_odd(
        graph: &TensorGraph<Aind>,
        automorphism: &[usize],
        header: usize,
        key: GroupKey,
    ) -> bool {
        let target_header = automorphism[header];
        let mut permutation = graph
            .node(header)
            .edges
            .iter()
            .filter_map(|edge_index| {
                let edge = graph.edge(*edge_index);
                let EdgeColor::Incidence(IncidenceRole::Group {
                    key: edge_key,
                    kind: SymmetryKind::Antisymmetric,
                }) = edge.data.data
                else {
                    return None;
                };
                if edge.vertices.0 != header || edge_key != key {
                    return None;
                }
                let EdgeOrigin::Incidence { member, .. } = edge.data.hidden else {
                    unreachable!("incidence edges retain member origins")
                };
                let target_slot = automorphism[edge.vertices.1];
                let target_member = graph
                    .node(target_header)
                    .edges
                    .iter()
                    .map(|edge_index| graph.edge(*edge_index))
                    .find_map(|target| {
                        let EdgeColor::Incidence(IncidenceRole::Group {
                            key: target_key,
                            kind: SymmetryKind::Antisymmetric,
                        }) = target.data.data
                        else {
                            return None;
                        };
                        if target.vertices == (target_header, target_slot) && target_key == key {
                            let EdgeOrigin::Incidence { member, .. } = target.data.hidden else {
                                unreachable!("incidence edges retain member origins")
                            };
                            Some(member)
                        } else {
                            None
                        }
                    })
                    .expect("an automorphism preserves typed tensor incidence");
                Some((member, target_member))
            })
            .collect::<Vec<_>>();
        permutation.sort_unstable_by_key(|(source, _)| *source);
        Permutation::from_map(permutation.into_iter().map(|(_, target)| target).collect()).sign()
            < 0
    }

    fn group_is_local(graph: &TensorGraph<Aind>, header: usize, key: GroupKey) -> bool {
        graph
            .node(header)
            .edges
            .iter()
            .filter_map(|edge_index| {
                let edge = graph.edge(*edge_index);
                matches!(
                edge.data.data,
                EdgeColor::Incidence(IncidenceRole::Group { key: edge_key, .. }) if edge_key == key
            )
            .then_some(edge.vertices.1)
            })
            .all(|slot| {
                graph
                    .node(slot)
                    .edges
                    .iter()
                    .find_map(|slot_edge_index| {
                        let slot_edge = graph.edge(*slot_edge_index);
                        matches!(slot_edge.data.data, EdgeColor::Contraction(_))
                            .then_some(slot_edge)
                    })
                    .is_some_and(|slot_edge| {
                        let other_slot = if slot_edge.vertices.0 == slot {
                            slot_edge.vertices.1
                        } else {
                            slot_edge.vertices.0
                        };
                        graph.node(other_slot).edges.iter().any(|other_edge_index| {
                    let other_edge = graph.edge(*other_edge_index);
                    other_edge.vertices.0 == header
                        && matches!(
                            other_edge.data.data,
                            EdgeColor::Incidence(IncidenceRole::Group { key: other_key, .. })
                                if other_key == key
                        )
                })
                    })
            })
    }
}

pub(super) fn canonize_network<Aind, F>(
    network: SymbolicNet<Aind>,
    mut new_dummy: F,
) -> Result<SymbolicNet<Aind>, CanonicalizationError>
where
    Aind: AbsInd + DummyAind + ParseableAind,
    F: FnMut(usize) -> Aind,
{
    Canonicalizer::new(&network, &mut new_dummy)?.run()
}

pub(crate) fn validate_tensor_symmetry<Aind: AbsInd + ParseableAind>(
    value: AtomView<'_>,
) -> Result<(), CanonicalizationError> {
    // Validate parser-recognized metric shorthand after its non-allocating
    // normalization; callers still parse or retain the original expression.
    let validation = value.normalize_dots().metric_shorthand_to_dot();
    let mut stack = vec![validation.as_view()];
    while let Some(value) = stack.pop() {
        match value {
            AtomView::Fun(function) => {
                let head = function.get_symbol();
                if head.has_tag(&SPENSO_TAG.tensor) && SymmetryKind::of(head).is_some() {
                    for (argument, value) in function.iter().enumerate() {
                        if Slot::<LibraryRep, Aind>::try_from(value).is_err() {
                            return Err(CanonicalizationError::InvalidIntrinsicArgument {
                                head,
                                argument,
                                expression: function.as_view().to_owned(),
                            });
                        }
                    }
                } else if head.has_tag(&SPENSO_TAG.tensor) {
                    for value in function.iter() {
                        let _ = TensorLayout::partial_group::<Aind>(
                            value,
                            &function.as_view().to_owned(),
                        )?;
                    }
                }
                stack.extend(function.iter());
            }
            AtomView::Add(add) => stack.extend(add.iter()),
            AtomView::Mul(product) => stack.extend(product.iter()),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                stack.push(base);
                stack.push(exponent);
            }
            AtomView::Num(_) | AtomView::Var(_) => {}
        }
    }
    Ok(())
}
