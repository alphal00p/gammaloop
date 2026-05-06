use std::{
    collections::{BTreeMap, BTreeSet},
    fmt::{Debug, Display},
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
    time::Instant,
};

use bincode::{Decode, Encode};
use linnet::{
    half_edge::{
        HedgeGraph, HedgeGraphError, NoData, NodeIndex,
        builder::HedgeGraphBuilder,
        involution::{EdgeData, Flow, Hedge},
        nodestore::NodeStorageOps,
        subgraph::{ModifySubSet, SuBitGraph, SubGraphLike, SubSetLike, SubSetOps},
        tree::SimpleTraversalTree,
    },
    permutation::Permutation,
    tree::{Forest, child_pointer::ParentChildStore, child_vec::ChildVecStore},
};
use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::{
    network::NetworkState,
    structure::{
        HasStructure, PermutedStructure, TensorStructure,
        abstract_index::AbstractIndex,
        permuted::{Perm, PermuteTensor},
        representation::{LibrarySlot, RepName},
        slot::{AbsInd, DualSlotTo, IsAbstractSlot},
    },
};

use super::{
    TensorNetworkError,
    library::{Library, LibraryTensor},
    profile::{self, Counter, Timer},
};

#[derive(
    Debug,
    Clone,
    Serialize,
    Deserialize,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
    PartialEq,
    Eq, // bincode_trait_derive::BorrowDecodeFromDecode,
)]
#[cfg_attr(
    feature = "shadowing",
    trait_decode(trait = symbolica::state::HasStateMap),
)]
pub struct NetworkGraph<K, FK = i8, Aind = AbstractIndex> {
    pub graph: NetworkHedgeGraph<K, FK, Aind>,
    pub slot_order: Vec<u8>,
    // #[bincode(with_serde)]
    // uncontracted: SuBitGraph,
}

type NetworkGraphBuilder<K, FK, Aind> =
    HedgeGraphBuilder<NetworkEdge<Aind>, NetworkNode<K, FK, Aind>>;

pub type NetworkNodeStore<K, FK, Aind> = Forest<NetworkNode<K, FK, Aind>, ParentChildStore<()>>;
pub type NetworkHedgeGraph<K, FK, Aind> =
    HedgeGraph<NetworkEdge<Aind>, NetworkNode<K, FK, Aind>, NoData, NetworkNodeStore<K, FK, Aind>>;

pub type ReadyNetworkOp<K, FK = i8, Aind = AbstractIndex> = (
    NetworkGraph<K, FK, Aind>,
    NetworkOp<FK>,
    Vec<NetworkLeaf<K, Aind>>,
);

#[derive(Debug, Clone)]
pub struct NetworkOperationRef<'a, K, FK = i8, Aind = AbstractIndex> {
    graph: &'a NetworkGraph<K, FK, Aind>,
    op_node: NodeIndex,
    op: &'a NetworkOp<FK>,
    children: Vec<NodeIndex>,
    subgraph: SuBitGraph,
}

impl<'a, K, FK, Aind> NetworkOperationRef<'a, K, FK, Aind> {
    pub fn graph(&self) -> &'a NetworkGraph<K, FK, Aind> {
        self.graph
    }

    pub fn op_node(&self) -> NodeIndex {
        self.op_node
    }

    pub fn op(&self) -> &'a NetworkOp<FK> {
        self.op
    }

    pub fn children(&self) -> &[NodeIndex] {
        &self.children
    }

    pub fn subgraph(&self) -> &SuBitGraph {
        &self.subgraph
    }

    pub fn leaf_count(&self) -> usize {
        self.children.len()
    }
}
#[derive(Debug, Clone)]
pub struct NetworkOperation<FK = i8> {
    op_node: NodeIndex,
    op: NetworkOp<FK>,
    children: Vec<NodeIndex>,
    subgraph: SuBitGraph,
}

impl<FK> NetworkOperation<FK> {
    pub fn op_node(&self) -> NodeIndex {
        self.op_node
    }

    pub fn op(&self) -> &NetworkOp<FK> {
        &self.op
    }

    pub fn children(&self) -> &[NodeIndex] {
        &self.children
    }

    pub fn subgraph(&self) -> &SuBitGraph {
        &self.subgraph
    }

    pub fn leaf_count(&self) -> usize {
        self.children.len()
    }
}

impl<K, FK: Clone, Aind> From<&NetworkOperationRef<'_, K, FK, Aind>> for NetworkOperation<FK> {
    fn from(op_ref: &NetworkOperationRef<'_, K, FK, Aind>) -> Self {
        Self {
            op_node: op_ref.op_node(),
            op: op_ref.op().clone(),
            children: op_ref.children().to_vec(),
            subgraph: op_ref.subgraph().clone(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct NetworkOperationBatchRef<'a, K, FK = i8, Aind = AbstractIndex> {
    operations: Vec<NetworkOperationRef<'a, K, FK, Aind>>,
    subgraph: SuBitGraph,
}

impl<'a, K, FK, Aind> NetworkOperationBatchRef<'a, K, FK, Aind> {
    pub fn operations(&self) -> &[NetworkOperationRef<'a, K, FK, Aind>] {
        &self.operations
    }

    pub fn iter(&self) -> impl Iterator<Item = &NetworkOperationRef<'a, K, FK, Aind>> {
        self.operations.iter()
    }

    pub fn subgraph(&self) -> &SuBitGraph {
        &self.subgraph
    }

    pub fn len(&self) -> usize {
        self.operations.len()
    }

    pub fn is_empty(&self) -> bool {
        self.operations.is_empty()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NetworkOperationReadiness {
    pub node_count: usize,
    pub hedge_count: usize,
    pub hidden_hedge_count: usize,
    pub cached_expression_node_count: usize,
    pub ready_operation_count: usize,
    pub batched_operation_count: usize,
    pub batched_subgraph_hedge_count: usize,
}

#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Serialize,
    Deserialize,
    Encode,
    bincode_trait_derive::Decode,
    // bincode_trait_derive::BorrowDecodeFromDecode,
)]
#[cfg_attr(
    feature = "shadowing",
    trait_decode(trait = symbolica::state::HasStateMap),
)]
pub enum NetworkEdge<Aind> {
    // Port,
    Head,
    Slot(LibrarySlot<Aind>),
}

impl<Aind: AbsInd> Display for NetworkEdge<Aind> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NetworkEdge::Head => write!(f, "Head"),
            NetworkEdge::Slot(slot) => write!(f, "Slot({})", slot),
        }
    }
}

impl<Aind> NetworkEdge<Aind> {
    pub fn is_head(&self) -> bool {
        matches!(self, NetworkEdge::Head)
    }

    pub fn is_slot(&self) -> bool {
        matches!(self, NetworkEdge::Slot(_))
    }
}

#[derive(
    Debug, Clone, PartialEq, Eq, Encode, bincode_trait_derive::Decode, Serialize, Deserialize,
)]
pub enum NetworkNode<LibKey, FunKey, Aind = AbstractIndex> {
    Leaf(NetworkLeaf<LibKey, Aind>),
    Op(NetworkOp<FunKey>),
    // Port,
}

impl<K, FK, Aind> NetworkNode<K, FK, Aind> {
    pub fn is_leaf(&self) -> bool {
        matches!(self, NetworkNode::Leaf(_))
    }

    pub fn is_op(&self) -> bool {
        matches!(self, NetworkNode::Op(_))
    }

    pub fn is_scalar(&self) -> bool {
        matches!(self, NetworkNode::Leaf(NetworkLeaf::Scalar(_)))
    }

    pub fn is_tensor(&self) -> bool {
        self.is_leaf() && !self.is_scalar()
    }
}

impl<K: Display, FK: Display, Aind> Display for NetworkNode<K, FK, Aind> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NetworkNode::Leaf(l) => write!(f, "{}", l),
            NetworkNode::Op(o) => write!(f, "{o}"),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Encode, Decode, Serialize, Deserialize)]
pub enum NetworkOp<FunKey> {
    Sum,
    Neg,
    Product,
    Power(i8),
    Function(FunKey),
}
impl<FK> NetworkOp<FK> {
    pub fn display_with(&self, function_map: impl Fn(&FK) -> String) -> String {
        match self {
            NetworkOp::Neg => "⊖".into(),
            NetworkOp::Product => "∏".into(),
            NetworkOp::Sum => "∑".into(),
            NetworkOp::Power(p) => format!("^( {p} )"),
            NetworkOp::Function(fun) => format!("Func({})", function_map(fun)),
        }
    }
}

impl<FunKey: Display> Display for NetworkOp<FunKey> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NetworkOp::Neg => write!(f, "⊖"),
            NetworkOp::Product => write!(f, "∏"),
            NetworkOp::Sum => write!(f, "∑"),
            NetworkOp::Power(p) => write!(f, "^( {p} )"),
            NetworkOp::Function(fun) => write!(f, "Func({})", fun),
        }
    }
}

#[derive(
    Debug, Clone, PartialEq, Eq, Encode, bincode_trait_derive::Decode, Serialize, Deserialize,
)]
pub enum NetworkLeaf<K, Aind = AbstractIndex> {
    LocalTensor(usize),
    LibraryKey {
        key: PermutedStructure<K>,
        indices: Vec<Aind>,
    },
    Scalar(usize),
}

impl<K, Aind> NetworkLeaf<K, Aind> {
    pub fn is_scalar(&self) -> bool {
        matches!(self, NetworkLeaf::Scalar(_))
    }
}

impl<K: Display, Aind> Display for NetworkLeaf<K, Aind> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NetworkLeaf::LibraryKey { key, .. } => write!(f, "Key:{key}"),
            NetworkLeaf::LocalTensor(l) => write!(f, "Tensor:{l}"),
            NetworkLeaf::Scalar(s) => write!(f, "Scalar:{s}"),
        }
    }
}

#[derive(Debug, Clone)]
pub enum NetworkLeafWithInds<K> {
    LocalTensor(usize),
    LibraryKey { key: K, inds: Vec<AbstractIndex> },
    Scalar(usize),
}

impl<K: Debug, FK: Debug, Aind: AbsInd>
    From<HedgeGraphBuilder<NetworkEdge<Aind>, NetworkNode<K, FK, Aind>>>
    for NetworkGraph<K, FK, Aind>
{
    fn from(builder: HedgeGraphBuilder<NetworkEdge<Aind>, NetworkNode<K, FK, Aind>>) -> Self {
        let _span = profile::span(Timer::BuildGraph);
        let graph: NetworkHedgeGraph<K, FK, Aind> = builder.build();
        let slot_order = vec![0; graph.n_hedges()];
        let mut g = Self {
            slot_order,
            graph,
            // uncontracted,
        };
        g.sync_order();
        g
    }
}
#[derive(Clone, Debug, Error)]
pub enum NetworkGraphError {
    #[error("Not head node")]
    NotHeadNode,
}

impl<K: Debug, FK: Debug, Aind: AbsInd> NetworkGraph<K, FK, Aind> {
    fn match_slots(
        _self_flow: Flow,
        self_data: EdgeData<&NetworkEdge<Aind>>,
        _other_flow: Flow,
        other_data: EdgeData<&NetworkEdge<Aind>>,
    ) -> bool {
        match (self_data.data, other_data.data) {
            (NetworkEdge::Slot(s), NetworkEdge::Slot(o)) => s.matches(o),
            _ => false,
        }
    }

    fn keep_first_edge(
        self_flow: Flow,
        self_data: EdgeData<NetworkEdge<Aind>>,
        _other_flow: Flow,
        _other_data: EdgeData<NetworkEdge<Aind>>,
    ) -> (Flow, EdgeData<NetworkEdge<Aind>>) {
        (self_flow, self_data)
    }

    fn sew_internal_tensor_slots(&mut self) {
        self.graph
            .sew(Self::match_slots, Self::keep_first_edge)
            .expect("sewing internal tensor slots should only connect dangling slot edges");
    }

    pub(crate) fn replace_node_deleting_self_loop_slots(
        &mut self,
        node: NodeIndex,
        node_data: NetworkNode<K, FK, Aind>,
    ) {
        let mut traced_slots: SuBitGraph = self.graph.empty_subgraph();
        for hedge in self.graph.iter_crown(node) {
            if self.graph[[&hedge]].is_slot() && self.graph.is_self_loop(hedge) {
                traced_slots.add(hedge);
            }
        }

        self.graph[node] = node_data;
        self.delete(&traced_slots);
        self.graph.node_store.check_and_set_nodes().unwrap();
    }

    fn set_tensor_slot_order(&mut self, node: NodeIndex, slots: &[LibrarySlot<Aind>]) {
        let mut slot_hedges = self
            .graph
            .iter_crown(node)
            .filter(|hedge| matches!(self.graph[[hedge]], NetworkEdge::Slot(_)))
            .collect::<Vec<_>>();

        for (order, slot) in slots.iter().enumerate() {
            let Some(position) = slot_hedges
                .iter()
                .position(|hedge| self.graph[[hedge]] == NetworkEdge::Slot(*slot))
            else {
                panic!(
                    "tensor graph is missing slot {slot} while assigning structural slot order; node crown: {:?}",
                    self.graph
                        .iter_crown(node)
                        .map(|hedge| self.graph[[&hedge]])
                        .collect::<Vec<_>>()
                );
            };
            let hedge = slot_hedges.remove(position);
            self.slot_order[hedge.0] = u8::try_from(order).unwrap_or_else(|_| {
                panic!(
                    "tensor graph slot order {order} for slot {slot} exceeds u8 storage; tensor has {} slots",
                    slots.len()
                )
            });
        }
    }

    pub fn slots(&self, nodeid: NodeIndex) -> Vec<LibrarySlot<Aind>> {
        let mut slots = Vec::new();
        let mut ord = Vec::new();
        if let NetworkNode::Leaf(_) = &self.graph[nodeid] {
            for n in self.graph.iter_crown(nodeid) {
                if let NetworkEdge::Slot(s) = self.graph[[&n]] {
                    slots.push(s);
                    ord.push(self.slot_order[n.0]);
                }
            }
        }

        let perm = Permutation::sort(&ord);
        perm.apply_slice_in_place(&mut slots);
        slots
    }

    pub fn state(&self) -> NetworkState {
        let dang = self.dangling_indices();
        if dang.is_empty() {
            NetworkState::Scalar
        } else if dang.iter().all(|a| a.rep.rep.is_self_dual()) {
            NetworkState::SelfDualTensor
        } else {
            NetworkState::Tensor
        }
    }

    pub fn sync_order(&mut self) {
        let _span = profile::span(Timer::SyncOrder);
        for (_, n, _) in self.graph.iter_nodes() {
            let slot_hedges = n
                .clone()
                .filter(|c| matches!(self.graph[[c]], NetworkEdge::Slot(_)))
                .collect::<Vec<_>>();
            if slot_hedges.iter().any(|h| self.slot_order[h.0] != 0) {
                continue;
            }

            let mut slots: BTreeMap<NetworkEdge<Aind>, Vec<Hedge>> = BTreeMap::new();
            for c in slot_hedges {
                slots
                    .entry(self.graph[self.graph[&c]])
                    .and_modify(|curr| curr.push(c))
                    .or_insert(vec![c]);
            }

            for (i, v) in slots.into_values().enumerate() {
                for h in v {
                    self.slot_order[h.0] = i as u8;
                }
            }
        }
    }

    pub fn inds(&self, nodeid: NodeIndex) -> Vec<Aind> {
        let mut slots = Vec::new();
        let mut ord = Vec::new();
        if let NetworkNode::Leaf(_) = &self.graph[nodeid] {
            for n in self.graph.iter_crown(nodeid) {
                if let NetworkEdge::Slot(s) = self.graph[[&n]] {
                    slots.push(s.aind);
                    ord.push(self.slot_order[n.0]);
                }
            }
        }

        let perm = Permutation::sort(&ord);
        perm.apply_slice_in_place(&mut slots);
        slots
    }

    pub fn get_lib_data<
        S,
        LT: LibraryTensor + Clone,
        L: Library<S, Key = K, Value = PermutedStructure<LT>>,
    >(
        &self,
        lib: &L,
        nodeid: NodeIndex,
    ) -> Option<LT::WithIndices>
    where
        K: Display + Debug,
        LT::WithIndices: PermuteTensor<Permuted = LT::WithIndices>,
        <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
            IsAbstractSlot<Aind = Aind>,
    {
        if let NetworkNode::Leaf(NetworkLeaf::LibraryKey { key, indices }) = &self.graph[nodeid] {
            let libt = lib.get(&key.structure).unwrap();
            let mappingperm = &key.index_permutation;
            let mut inds = indices.clone();

            // println!("Mapping perm: {mappingperm}");
            mappingperm.apply_slice_in_place_inv(&mut inds);
            // println!("Inds: {inds:?}");
            let libt_with_indices = libt.structure.with_indices(&inds).unwrap();
            // libt_with_indices.index_permutation = k.index_permutation.clone();

            Some(libt_with_indices.permute_inds())
        } else {
            None
        }
    }

    pub fn splice_descendents_of(&mut self, replacement: Self)
    where
        K: Clone + Debug,
    {
        let _span = profile::span(Timer::Splice);
        profile::bump(Counter::Splice, 1);
        // println!(
        //     "Joining \n{} with\n {}",
        //     self.dot_simple(),
        //     replacement.dot_simple()
        // );
        //
        self.slot_order.extend(replacement.slot_order);
        self.graph
            .join_mut(
                replacement.graph,
                |sf, sd, of, od| {
                    let flow_match = sf == -of;
                    let desc_match = sd == od;
                    // if desc_match{

                    // println!("looking at {sd} vs {od}, flow_match: {flow_match},sf  {sf:?},of {of:?}, desc_match: {desc_match}");
                    // }
                    // if flow_match && desc_match {
                    // sd.data
                    // println!("Splicing");
                    // }
                    flow_match && desc_match
                },
                |sf, sd, _, _| (sf, sd),
            )
            .unwrap();
    }

    pub fn extract<S: SubSetLike<Base = SuBitGraph>>(&mut self, subgraph: &S) -> Self
    where
        K: Clone + Display,
        FK: Clone + Display,
    {
        let _span = profile::span(Timer::GraphExtract);
        let mut left = Hedge(0);
        let mut extracted = Hedge(self.graph.n_hedges());
        while left < extracted {
            if !subgraph.includes(&left) {
                //left is in the right place
                left.0 += 1;
            } else {
                //left needs to be swapped
                extracted.0 -= 1;
                if !subgraph.includes(&extracted) {
                    //only with an extracted that is in the wrong spot
                    self.slot_order.swap(left.0, extracted.0);
                    left.0 += 1;
                }
            }
        }

        let extracted = self.graph.extract(
            subgraph,
            |a| a.map(Clone::clone),
            |a| a,
            |a| {
                // println!("a{a}");
                a.clone()
            },
            |a| a,
        );

        let slot_order = self.slot_order.split_off(left.0);

        Self {
            slot_order,
            graph: extracted,
        }
    }

    pub fn find_all_ready_ops(&mut self) -> Vec<ReadyNetworkOp<K, FK, Aind>>
    where
        K: Clone + Display,
        FK: Clone + Display,
    {
        let tt: SimpleTraversalTree<ParentChildStore<()>> = self.expr_tree().cast();
        let head = self.head();
        let mut out = vec![];
        let root_node = self.graph.node_id(head);
        let involution = self.graph.as_ref().clone();

        for nid in tt.iter_preorder_tree_nodes(&involution, root_node) {
            if let NetworkNode::Op(op) = &self.graph[nid] {
                let mut leaves = Vec::new();
                let mut subgraph: SuBitGraph = self.graph.empty_subgraph();
                let ok = tt.iter_children(nid, &self.graph).all(|child| {
                    for h in self.graph.iter_crown(child) {
                        subgraph.add(h);
                    }
                    if let NetworkNode::Leaf(l) = &self.graph[child] {
                        leaves.push(l.clone());
                        true
                    } else {
                        false
                    }
                });
                if ok {
                    let op = op.clone();
                    let extracted = self.extract(&subgraph);

                    out.push((extracted, op, leaves))
                }
            }
        }
        out
    }
    pub fn extract_next_ready_op(&mut self) -> Option<(Self, NetworkOp<FK>)>
    where
        K: Clone + Display,
        FK: Clone + Display,
    {
        // build a traversal over *all* internal edges
        let tt: SimpleTraversalTree<ChildVecStore<()>> = self.expr_tree().cast();
        let head = self.head();
        let root_node = self.graph.node_id(head);

        // println!("//tree:\n{}", self.graph.dot(&tt.tree_subgraph));
        // println!(
        //     "//tree:\n{}",
        //     self.graph.dot(&tt.tree_subgraph(self.graph.as_ref()))
        // );
        // println!("tree{:#?}", tt);
        // println!(
        //     "tree{}",
        //     tt.debug_draw(|a| match a {
        //         TTRoot::Child(a) => a.to_string(),
        //         TTRoot::Root => "Root".into(),
        //         _ => "".into(),
        //     })
        // );
        // look for the first op node whose children are all leaves
        for nid in tt.iter_preorder_tree_nodes(&self.graph, root_node) {
            if let NetworkNode::Op(op) = &self.graph[nid] {
                let mut subgraph: SuBitGraph = self.graph.empty_subgraph();

                let mut all_leaves = true;
                let mut has_children = false;

                for h in self.graph.iter_crown(nid) {
                    subgraph.add(h);
                }
                // println!("{nid}");

                for child in tt.iter_children(nid, &self.graph) {
                    // println!("{child}");
                    has_children = true;
                    match &self.graph[child] {
                        NetworkNode::Leaf(_) => {}
                        _ => {
                            all_leaves = false;
                            break;
                        }
                    }
                    for h in self.graph.iter_crown(child) {
                        subgraph.add(h);
                    }
                }
                if all_leaves && has_children {
                    let op = op.clone();

                    // println!(
                    //     "Extracting the: {}",
                    //     self.dot_impl(|a| a.to_string(), |_| "None".to_string(), |a| a.to_string())
                    // );

                    // println!(
                    //     "Extracting the: {}",
                    //     self.dot_impl_of(
                    //         &subgraph,
                    //         |a| a.to_string(),
                    //         |_| "None".to_string(),
                    //         |a| a.to_string()
                    //     )
                    // );

                    self.graph.check().unwrap();

                    let extracted = self.extract(&subgraph);

                    extracted.graph.check().unwrap();

                    return Some((extracted, op));
                }
            }
        }
        None
    }

    pub fn extract_next_ready_ref_op(&mut self) -> Option<(Self, NetworkOp<FK>)>
    where
        K: Clone + Display,
        FK: Clone + Display,
    {
        let (subgraph, op) = {
            let ready = self.ready_operation_ref_in_expr_tree_order()?;
            (ready.subgraph().clone(), ready.op().clone())
        };

        self.graph.check().unwrap();
        let extracted = self.extract(&subgraph);
        extracted.graph.check().unwrap();

        Some((extracted, op))
    }

    pub fn ready_operation_ref_in_expr_tree_order(
        &self,
    ) -> Option<NetworkOperationRef<'_, K, FK, Aind>> {
        let tt: SimpleTraversalTree<ChildVecStore<()>> = self.expr_tree().cast();
        let root_node = self.graph.node_id(self.head());

        for nid in tt.iter_preorder_tree_nodes(&self.graph, root_node) {
            let NetworkNode::Op(op) = &self.graph[nid] else {
                continue;
            };

            let mut children = Vec::new();
            let mut all_leaves = true;
            for child in tt.iter_children(nid, &self.graph) {
                if matches!(&self.graph[child], NetworkNode::Leaf(_)) {
                    children.push(child);
                } else {
                    all_leaves = false;
                    break;
                }
            }

            if all_leaves && !children.is_empty() {
                let subgraph = self.operation_subgraph(nid, &children);
                return Some(NetworkOperationRef {
                    graph: self,
                    op_node: nid,
                    op,
                    children,
                    subgraph,
                });
            }
        }

        None
    }

    pub fn cache_expr_tree_roots(&mut self) -> usize {
        let hidden: SuBitGraph = self.graph.empty_subgraph();
        self.cache_expr_tree_roots_ignoring(&hidden)
    }

    pub fn cache_expr_tree_roots_ignoring<S>(&mut self, hidden: &S) -> usize
    where
        S: SubSetLike<Base = SuBitGraph>,
    {
        let tt: SimpleTraversalTree<ChildVecStore<()>> = self.expr_tree_ignoring(hidden).cast();
        let head = self.head();
        let root_node = self.graph.node_id(head);
        let new_roots = tt
            .iter_preorder_tree_nodes(&self.graph, root_node)
            .map(|node| tt.root_hedge(node).into())
            .collect::<Vec<_>>();

        self.graph.node_store.reroot_many(new_roots)
    }

    pub fn cached_expr_children(&self, node: NodeIndex) -> Vec<NodeIndex> {
        let hidden: SuBitGraph = self.graph.empty_subgraph();
        self.cached_expr_children_ignoring(node, &hidden)
    }

    pub fn cached_expr_children_ignoring<S>(&self, node: NodeIndex, hidden: &S) -> Vec<NodeIndex>
    where
        S: SubSetLike<Base = SuBitGraph>,
    {
        let mut children = self
            .graph
            .iter_crown(node)
            .filter_map(|hedge| {
                if hidden.includes(&hedge) {
                    return None;
                }

                if self.graph[[&hedge]].is_slot() {
                    return None;
                }

                let other = self.graph.inv(hedge);
                if hidden.includes(&other) {
                    return None;
                }

                if other == hedge {
                    return None;
                }

                let other_node = self.graph.node_id(other);
                if other_node == node {
                    return None;
                }

                (self.graph.node_store.root_hedge_for_node(other_node) == other)
                    .then_some(other_node)
            })
            .collect::<Vec<_>>();

        children.sort();
        children.dedup();
        children
    }

    pub fn cached_expr_preorder_nodes(&self) -> Vec<NodeIndex> {
        let hidden: SuBitGraph = self.graph.empty_subgraph();
        self.cached_expr_preorder_nodes_ignoring(&hidden)
    }

    pub fn cached_expr_preorder_nodes_ignoring<S>(&self, hidden: &S) -> Vec<NodeIndex>
    where
        S: SubSetLike<Base = SuBitGraph>,
    {
        let root_node = self.graph.node_id(self.head());
        let mut stack = vec![root_node];
        let mut seen = BTreeSet::new();
        let mut nodes = Vec::new();

        while let Some(node) = stack.pop() {
            if !seen.insert(node) {
                continue;
            }

            nodes.push(node);
            let mut children = self.cached_expr_children_ignoring(node, hidden);
            children.reverse();
            stack.extend(children);
        }

        nodes
    }

    pub fn ready_operation_ref(&self) -> Option<NetworkOperationRef<'_, K, FK, Aind>> {
        let hidden: SuBitGraph = self.graph.empty_subgraph();
        self.ready_operation_ref_ignoring(&hidden)
    }

    pub fn ready_operation_ref_ignoring<S>(
        &self,
        hidden: &S,
    ) -> Option<NetworkOperationRef<'_, K, FK, Aind>>
    where
        S: SubSetLike<Base = SuBitGraph>,
    {
        for nid in self.cached_expr_preorder_nodes_ignoring(hidden) {
            let Some((op, children)) = self.ready_operation_parts_ignoring(nid, hidden) else {
                continue;
            };

            let subgraph = self.operation_subgraph_ignoring(nid, &children, hidden);
            return Some(NetworkOperationRef {
                graph: self,
                op_node: nid,
                op,
                children,
                subgraph,
            });
        }

        None
    }

    pub fn ready_operation_refs(&self) -> Vec<NetworkOperationRef<'_, K, FK, Aind>> {
        let hidden: SuBitGraph = self.graph.empty_subgraph();
        self.ready_operation_refs_ignoring(&hidden)
    }

    /// Return ready operations from the cached expression tree without building a
    /// second disjointness filter. Each selected operation is an expression-tree
    /// node plus its cached tree children, so one traversal cannot select
    /// overlapping operation subgraphs. Debug builds assert that invariant.
    pub fn ready_operation_tree_refs_ignoring<S>(
        &self,
        hidden: &S,
    ) -> Vec<NetworkOperationRef<'_, K, FK, Aind>>
    where
        S: SubSetLike<Base = SuBitGraph>,
    {
        let mut ready = Vec::new();

        #[cfg(debug_assertions)]
        let mut used: SuBitGraph = self.graph.empty_subgraph();

        let cached_nodes = {
            let _span = profile::span(Timer::ExecuteReadyPreorder);
            self.cached_expr_preorder_nodes_ignoring(hidden)
        };

        for nid in cached_nodes {
            let parts = {
                let _span = profile::span(Timer::ExecuteReadyParts);
                self.ready_operation_parts_ignoring(nid, hidden)
            };
            let Some((op, children)) = parts else {
                continue;
            };

            let subgraph = {
                let _span = profile::span(Timer::ExecuteReadySubgraph);
                self.operation_subgraph_ignoring(nid, &children, hidden)
            };

            #[cfg(debug_assertions)]
            {
                debug_assert!(subgraph.empty_intersection(&used));
                used.union_with(&subgraph);
            }

            ready.push(NetworkOperationRef {
                graph: self,
                op_node: nid,
                op,
                children,
                subgraph,
            });
        }

        ready
    }

    pub fn ready_operation_refs_ignoring<S>(
        &self,
        hidden: &S,
    ) -> Vec<NetworkOperationRef<'_, K, FK, Aind>>
    where
        S: SubSetLike<Base = SuBitGraph>,
    {
        let mut used: SuBitGraph = self.graph.empty_subgraph();
        let mut ready = Vec::new();

        let cached_nodes = {
            let _span = profile::span(Timer::ExecuteReadyPreorder);
            self.cached_expr_preorder_nodes_ignoring(hidden)
        };

        for nid in cached_nodes {
            let parts = {
                let _span = profile::span(Timer::ExecuteReadyParts);
                self.ready_operation_parts_ignoring(nid, hidden)
            };
            let Some((op, children)) = parts else {
                continue;
            };

            let subgraph = {
                let _span = profile::span(Timer::ExecuteReadySubgraph);
                self.operation_subgraph_ignoring(nid, &children, hidden)
            };
            let overlaps_used = {
                let _span = profile::span(Timer::ExecuteReadyIntersection);
                !subgraph.empty_intersection(&used)
            };
            if overlaps_used {
                continue;
            }

            {
                let _span = profile::span(Timer::ExecuteReadyUnion);
                used.union_with(&subgraph);
            }
            ready.push(NetworkOperationRef {
                graph: self,
                op_node: nid,
                op,
                children,
                subgraph,
            });
        }

        ready
    }

    pub fn ready_operation_batch(&self) -> NetworkOperationBatchRef<'_, K, FK, Aind> {
        let hidden: SuBitGraph = self.graph.empty_subgraph();
        self.ready_operation_batch_ignoring(&hidden)
    }

    pub fn ready_operation_batch_ignoring<S>(
        &self,
        hidden: &S,
    ) -> NetworkOperationBatchRef<'_, K, FK, Aind>
    where
        S: SubSetLike<Base = SuBitGraph>,
    {
        let operations = self.ready_operation_refs_ignoring(hidden);
        let mut subgraph: SuBitGraph = self.graph.empty_subgraph();
        for op_ref in &operations {
            let _span = profile::span(Timer::ExecuteReadyUnion);
            subgraph.union_with(op_ref.subgraph());
        }

        NetworkOperationBatchRef {
            operations,
            subgraph,
        }
    }

    fn ready_operation_parts_ignoring<'a, S>(
        &'a self,
        node: NodeIndex,
        hidden: &S,
    ) -> Option<(&'a NetworkOp<FK>, Vec<NodeIndex>)>
    where
        S: SubSetLike<Base = SuBitGraph>,
    {
        let NetworkNode::Op(op) = &self.graph[node] else {
            return None;
        };

        let children = self.cached_expr_children_ignoring(node, hidden);
        if children.is_empty()
            || children
                .iter()
                .any(|child| !matches!(&self.graph[*child], NetworkNode::Leaf(_)))
        {
            return None;
        }

        Some((op, children))
    }

    pub fn operation_readiness_diagnostics(&self) -> NetworkOperationReadiness {
        let hidden: SuBitGraph = self.graph.empty_subgraph();
        self.operation_readiness_diagnostics_ignoring(&hidden)
    }

    pub fn operation_readiness_diagnostics_ignoring<S>(
        &self,
        hidden: &S,
    ) -> NetworkOperationReadiness
    where
        S: SubSetLike<Base = SuBitGraph>,
    {
        let cached_expression_nodes = self.cached_expr_preorder_nodes_ignoring(hidden);
        let ready_operation_count = cached_expression_nodes
            .iter()
            .filter(|node| {
                self.ready_operation_parts_ignoring(**node, hidden)
                    .is_some()
            })
            .count();
        let ready_batch = self.ready_operation_batch_ignoring(hidden);

        NetworkOperationReadiness {
            node_count: self.graph.n_nodes(),
            hedge_count: self.graph.n_hedges(),
            hidden_hedge_count: hidden.n_included(),
            cached_expression_node_count: cached_expression_nodes.len(),
            ready_operation_count,
            batched_operation_count: ready_batch.len(),
            batched_subgraph_hedge_count: ready_batch.subgraph().n_included(),
        }
    }

    pub fn operation_subgraph(&self, op_node: NodeIndex, children: &[NodeIndex]) -> SuBitGraph {
        let hidden: SuBitGraph = self.graph.empty_subgraph();
        self.operation_subgraph_ignoring(op_node, children, &hidden)
    }

    pub fn operation_subgraph_ignoring<S>(
        &self,
        op_node: NodeIndex,
        children: &[NodeIndex],
        hidden: &S,
    ) -> SuBitGraph
    where
        S: SubSetLike<Base = SuBitGraph>,
    {
        let mut subgraph: SuBitGraph = self.graph.empty_subgraph();
        self.add_visible_crown_to_subgraph(op_node, hidden, &mut subgraph);
        for child in children {
            self.add_visible_crown_to_subgraph(*child, hidden, &mut subgraph);
        }
        subgraph
    }

    fn add_visible_crown_to_subgraph<S>(
        &self,
        node: NodeIndex,
        hidden: &S,
        subgraph: &mut SuBitGraph,
    ) where
        S: SubSetLike<Base = SuBitGraph>,
    {
        for hedge in self.graph.iter_crown(node) {
            let other = self.graph.inv(hedge);
            if !hidden.includes(&hedge) && !hidden.includes(&other) {
                subgraph.add(hedge);
            }
        }
    }

    pub fn sub_expression(&self, nid: NodeIndex) -> Result<SimpleTraversalTree, NetworkGraphError> {
        let include_hedge = self
            .graph
            .iter_crown(nid)
            .find(|i| !matches!(self.graph[[i]], NetworkEdge::Slot(_)));

        let headgraph: SuBitGraph = self.graph.from_filter(|a| matches!(a, NetworkEdge::Head));

        if include_hedge.is_some() {
            Ok(SimpleTraversalTree::depth_first_traverse(
                &self.graph,
                &headgraph,
                &nid,
                include_hedge,
            )
            .unwrap())
        } else {
            Err(NetworkGraphError::NotHeadNode)
        }
    }
    pub fn one() -> Self {
        Self::mul_graph(0)
    }

    pub fn shift_scalars(&mut self, shift: usize) {
        let _span = profile::span(Timer::ShiftScalars);
        profile::bump(Counter::ShiftScalars, 1);
        self.graph.iter_nodes_mut().for_each(|(_, _, d)| {
            if let NetworkNode::Leaf(NetworkLeaf::Scalar(s)) = d {
                *s += shift;
            }
        });
    }

    pub fn shift_tensors(&mut self, shift: usize) {
        let _span = profile::span(Timer::ShiftTensors);
        profile::bump(Counter::ShiftTensors, 1);
        self.graph.iter_nodes_mut().for_each(|(_, _, d)| {
            if let NetworkNode::Leaf(NetworkLeaf::LocalTensor(s)) = d {
                *s += shift;
            }
        });
    }

    pub fn delete<S: SubSetLike<Base = SuBitGraph>>(&mut self, subgraph: &S) {
        let mut left = Hedge(0);
        let mut extracted = Hedge(self.graph.n_hedges());
        while left < extracted {
            if !subgraph.includes(&left) {
                //left is in the right place
                left.0 += 1;
            } else {
                //left needs to be swapped
                extracted.0 -= 1;
                if !subgraph.includes(&extracted) {
                    // println!("{extracted}<=>{left}");
                    //only with an extracted that is in the wrong spot
                    self.slot_order.swap(left.0, extracted.0);
                    left.0 += 1;
                }
            }
        }
        self.graph.delete_hedges(subgraph);
    }
    pub fn identify_nodes_without_self_edges(
        &mut self,
        nodes: &[NodeIndex],
        node_data: NetworkNode<K, FK, Aind>,
    ) -> NodeIndex {
        let _span = profile::span(Timer::IdentifyNodes);
        let (n, sub) = self
            .graph
            .identify_nodes_without_self_edges::<SuBitGraph>(nodes, node_data);

        self.graph.forget_identification_history();
        self.graph.node_store.check_and_set_nodes().unwrap();
        self.delete(&sub);

        self.graph.node_store.check_and_set_nodes().unwrap();
        // println!(
        //     "identify_nodes_without_self_edges res:{}",
        //     self.dot_simple()
        // );
        n
    }

    pub fn identify_subgraph_nodes_without_deleting_self_edges<S>(
        &mut self,
        subgraph: &S,
        node_data: NetworkNode<K, FK, Aind>,
        ignored: &mut SuBitGraph,
    ) -> Option<NodeIndex>
    where
        S: SubSetLike<Base = SuBitGraph>,
    {
        let _span = profile::span(Timer::IdentifyNodes);
        self.graph
            .identify_nodes_of_subgraph_marking_self_edges(subgraph, node_data, ignored)
    }

    pub fn finish_deferred_node_identifications(&mut self) {
        self.graph.forget_identification_history();
        self.graph.node_store.check_and_set_nodes().unwrap();
    }

    pub fn identify_nodes_without_self_edges_merge_heads(
        &mut self,
        nodes: &[NodeIndex],
        node_data: NetworkNode<K, FK, Aind>,
    ) -> NodeIndex {
        // println!("Identifying:{:?}", nodes);
        let (n, mut sub) = self
            .graph
            .identify_nodes_without_self_edges::<SuBitGraph>(nodes, node_data);

        let mut first = true;
        for h in self.graph.iter_crown(n) {
            // println!("{h}");
            if self.graph[[&h]].is_head() && self.graph.inv(h) != h {
                if first {
                    first = false;
                } else {
                    sub.add(h);
                    sub.add(self.graph.inv(h));
                }
            }
        }

        // println!(
        //     "Deleting: res:{}",
        //     self.dot_impl_of(&sub, |i| i.to_string(), |_| "".into(), |i| i.to_string())
        // );
        self.graph.forget_identification_history();
        self.graph.node_store.check_and_set_nodes().unwrap();
        self.delete(&sub);
        self.graph.node_store.check_and_set_nodes().unwrap();
        // println!(
        //     "identify_nodes_without_self_edges_merge_heads res:{}",
        //     self.dot_simple()
        // );

        n
    }

    pub fn dot(&self) -> String
    where
        K: Display,
        FK: Display,
    {
        self.graph.dot_impl(
            &self.graph.full_filter(),
            "",
            &|_| None,
            &|e| {
                if let NetworkEdge::Slot(s) = e {
                    Some(format!("label=\"{s}\""))
                } else {
                    None
                }
            },
            &|n| match n {
                NetworkNode::Leaf(l) => match l {
                    NetworkLeaf::LibraryKey { key, .. } => Some(format!("label= \"L{key}\"")),
                    NetworkLeaf::LocalTensor(l) => Some(format!("label = \"T{l}\"")),
                    NetworkLeaf::Scalar(s) => Some(format!("label = \"S{s}\"")),
                },
                NetworkNode::Op(o) => Some(format!("label = \"{o}\"")),
            },
        )
    }

    pub fn dot_impl(
        &self,
        scalar_disp: impl Fn(usize) -> String,
        library_disp: impl Fn(&K) -> String,
        tensor_disp: impl Fn(usize) -> String,
        function_disp: impl Fn(&FK) -> String,
    ) -> String
// where
        // K: Display,
    {
        self.graph.dot_impl(
            &self.graph.full_filter(),
            "",
            &|_| None,
            &|e| {
                if let NetworkEdge::Slot(s) = e {
                    // #[cfg(feature = "shadowing")]
                    // {
                    //     use symbolica::atom::AtomCore;

                    //     use crate::shadowing::symbolica_utils::SpensoPrintSettings;

                    //     Some(format!(
                    //         "label=\"{s}={}\"",
                    //         s.to_atom()
                    //             .printer(SpensoPrintSettings::compact().nice_symbolica())
                    //     ))
                    // }
                    // #[cfg(not(feature = "shadowing"))]
                    Some(format!("label=\"{s}\""))
                } else {
                    None
                }
            },
            &|n| match n {
                NetworkNode::Leaf(l) => match l {
                    NetworkLeaf::LibraryKey { key, .. } => {
                        Some(format!("label= \"L:{}\"", library_disp(&key.structure)))
                    }
                    NetworkLeaf::LocalTensor(l) => {
                        Some(format!("label = \"T:{}\"", tensor_disp(*l)))
                    }
                    NetworkLeaf::Scalar(s) => Some(format!("label = \"S:{}\"", scalar_disp(*s))),
                },
                NetworkNode::Op(o) => {
                    Some(format!("label = \"{}\"", o.display_with(&function_disp)))
                }
            },
        )
    }
    pub fn dot_impl_of<S: SubGraphLike>(
        &self,
        subgraph: &S,
        scalar_disp: impl Fn(usize) -> String,
        library_disp: impl Fn(&K) -> String,
        tensor_disp: impl Fn(usize) -> String,
        function_disp: impl Fn(&FK) -> String,
    ) -> String
// where
        // K: Display,
    {
        self.graph.dot_impl(
            subgraph,
            "",
            &|_| None,
            &|e| {
                if let NetworkEdge::Slot(s) = e {
                    Some(format!("label=\"{s}\""))
                } else {
                    None
                }
            },
            &|n| match n {
                NetworkNode::Leaf(l) => match l {
                    NetworkLeaf::LibraryKey { key, .. } => {
                        Some(format!("label= \"L:{}\"", library_disp(&key.structure)))
                    }
                    NetworkLeaf::LocalTensor(l) => {
                        Some(format!("label = \"T:{}\"", tensor_disp(*l)))
                    }
                    NetworkLeaf::Scalar(s) => Some(format!("label = \"S:{}\"", scalar_disp(*s))),
                },
                NetworkNode::Op(o) => {
                    Some(format!("label = \"{}\"", o.display_with(&function_disp)))
                }
            },
        )
    }

    pub fn dot_simple(&self) -> String
// where
        // K: Display,
    {
        self.dot_impl(
            |i| i.to_string(),
            |_| "".into(),
            |i| i.to_string(),
            |_| "".into(),
        )
    }

    pub fn zero() -> Self {
        Self::add_graph(0, &[])
    }

    fn head_builder(
        node: NetworkNode<K, FK, Aind>,
    ) -> (NetworkGraphBuilder<K, FK, Aind>, NodeIndex) {
        let mut graph = HedgeGraphBuilder::new();
        let head = graph.add_node(node);
        graph.add_external_edge(head, NetworkEdge::Head, true, Flow::Source);
        (graph, head)
    }

    pub fn neg_graph() -> Self {
        let (mut graph, head) = Self::head_builder(NetworkNode::Op(NetworkOp::Neg));
        graph.add_external_edge(head, NetworkEdge::Head, true, Flow::Sink);
        graph.into()
    }

    pub fn pow_graph(pow: i8, slots: &[LibrarySlot<Aind>]) -> Self {
        let (mut graph, head) = Self::head_builder(NetworkNode::Op(NetworkOp::Power(pow)));
        graph.add_external_edge(head, NetworkEdge::Head, true, Flow::Sink);
        if pow % 2 == 0 {
            for s in slots {
                let orientation = s.rep_name().orientation();
                graph.add_external_edge(head, NetworkEdge::Slot(*s), orientation, Flow::Sink);
            }
        }
        graph.into()
    }

    pub fn fun_graph(key: FK) -> Self {
        let (mut graph, head) = Self::head_builder(NetworkNode::Op(NetworkOp::Function(key)));
        graph.add_external_edge(head, NetworkEdge::Head, true, Flow::Sink);
        graph.into()
    }

    pub fn mul_graph(n: usize) -> Self {
        let (mut graph, head) = Self::head_builder(NetworkNode::Op(NetworkOp::Product));
        for _ in 0..n {
            graph.add_external_edge(head, NetworkEdge::Head, true, Flow::Sink);
        }

        graph.into()
    }

    pub fn add_graph(n: usize, slots: &[LibrarySlot<Aind>]) -> Self {
        let (mut graph, head) = Self::head_builder(NetworkNode::Op(NetworkOp::Sum));
        for _ in 0..n {
            graph.add_external_edge(head, NetworkEdge::Head, true, Flow::Sink);
        }

        for s in slots {
            let orientation = s.rep_name().orientation();
            graph.add_external_edge(head, NetworkEdge::Slot(*s), orientation, Flow::Source);
            for _ in 0..n {
                graph.add_external_edge(head, NetworkEdge::Slot(*s), orientation, Flow::Sink);
            }
        }

        graph.into()
    }

    pub fn scalar(pos: usize) -> Self {
        let mut graph = HedgeGraphBuilder::new();

        let head = graph.add_node(NetworkNode::Leaf(NetworkLeaf::Scalar(pos)));
        graph.add_external_edge(head, NetworkEdge::Head, true, Flow::Source);

        graph.into()
    }

    pub fn key(key: PermutedStructure<K>) -> Self
    where
        K: TensorStructure,
        K::Slot: IsAbstractSlot<Aind = Aind>,
    {
        let slots = key
            .structure
            .external_structure_iter()
            .map(|a| a.to_lib())
            .collect::<Vec<_>>();

        let indices = slots.iter().map(|slot| slot.aind()).collect::<Vec<_>>();
        let (mut graph, head) =
            Self::head_builder(NetworkNode::Leaf(NetworkLeaf::LibraryKey { key, indices }));

        for lib in &slots {
            let orientation = lib.rep_name().orientation();
            graph.add_external_edge(head, NetworkEdge::Slot(*lib), orientation, Flow::Source);
        }
        let mut graph = Self::from(graph);
        graph.set_tensor_slot_order(head, &slots);
        graph.sew_internal_tensor_slots();
        graph
    }

    pub fn tensor<T: TensorStructure>(
        tensor: &T,
        node: NetworkLeaf<K, Aind>,
    ) -> NetworkGraph<K, FK, Aind>
    where
        T::Slot: IsAbstractSlot<Aind = Aind>,
    {
        let (mut graph, head) = Self::head_builder(NetworkNode::Leaf(node));
        let mut slots = Vec::new();

        for s in tensor.external_structure_iter() {
            let lib = s.to_lib();
            slots.push(lib);

            let orientation = lib.rep_name().orientation();
            graph.add_external_edge(head, NetworkEdge::Slot(lib), orientation, Flow::Source);
        }
        let mut graph = Self::from(graph);
        graph.set_tensor_slot_order(head, &slots);
        graph.sew_internal_tensor_slots();
        graph
    }

    fn match_heads(
        self_flow: Flow,
        self_data: EdgeData<&NetworkEdge<Aind>>,
        other_flow: Flow,
        other_data: EdgeData<&NetworkEdge<Aind>>,
    ) -> bool {
        if let (NetworkEdge::Head, NetworkEdge::Head) = (self_data.data, other_data.data) {
            self_flow == -other_flow
        } else {
            false
        }
    }

    // fn match_indices(
    //     _sf: Flow,
    //     self_data: EdgeData<&NetworkEdge>,
    //     _of: Flow,
    //     other_data: EdgeData<&NetworkEdge>,
    // ) -> bool {
    //     if let (NetworkEdge::Slot(s), NetworkEdge::Slot(o)) = (self_data.data, other_data.data) {
    //         s.matches(o)
    //     } else {
    //         false
    //     }
    // }

    fn prod_match(
        self_flow: Flow,
        self_data: EdgeData<&NetworkEdge<Aind>>,
        other_flow: Flow,
        other_data: EdgeData<&NetworkEdge<Aind>>,
    ) -> bool {
        match (self_data.data, other_data.data) {
            (NetworkEdge::Head, NetworkEdge::Head) => self_flow == -other_flow,
            (NetworkEdge::Slot(s), NetworkEdge::Slot(o)) => s.matches(o),
            _ => false,
        }
    }

    fn add_match(
        sf: Flow,
        sd: EdgeData<&NetworkEdge<Aind>>,
        of: Flow,
        od: EdgeData<&NetworkEdge<Aind>>,
    ) -> bool {
        match (sd.data, od.data) {
            (NetworkEdge::Head, NetworkEdge::Head) => sf == -of,
            (NetworkEdge::Slot(s), NetworkEdge::Slot(o)) => s == o && sf == -of,
            _ => false,
        }
    }

    pub fn head(&self) -> Hedge {
        let exts: SuBitGraph = self.graph.external_filter();
        let head = exts
            .included_iter()
            .find(|i| self.graph[[i]] == NetworkEdge::Head);
        head.unwrap()
    }

    pub fn dangling_indices(&self) -> Vec<LibrarySlot<Aind>> {
        let _span = profile::span(Timer::DanglingScan);
        profile::bump(Counter::DanglingScan, 1);
        let exts: SuBitGraph = self.graph.external_filter();
        exts.included_iter()
            .filter_map(|i| {
                if let NetworkEdge::Slot(s) = self.graph[[&i]] {
                    Some(s)
                } else {
                    None
                }
            })
            .collect()
    }

    pub fn n_dangling(&self) -> usize {
        let _span = profile::span(Timer::DanglingScan);
        profile::bump(Counter::DanglingScan, 1);
        self.graph
            .external_filter::<SuBitGraph>()
            .included_iter()
            .filter(|i| matches!(self.graph[[i]], NetworkEdge::Slot(_)))
            .count()
    }

    pub fn n_nodes(&self) -> usize {
        self.graph.n_nodes()
    }

    #[allow(clippy::result_large_err, clippy::type_complexity)]
    pub fn result(
        &self,
    ) -> Result<
        (&NetworkNode<K, FK, Aind>, NodeIndex, Vec<LibrarySlot<Aind>>),
        TensorNetworkError<K, FK>,
    >
    where
        K: Display,
        FK: Display,
    {
        let mut n_heads = 0;
        let _headgraph: SuBitGraph = self.graph.from_filter(|a| match a {
            NetworkEdge::Head => {
                n_heads += 1;
                true
            }
            // NetworkEdge::Port => true,
            _ => false,
        });

        match n_heads {
            0 => Err(TensorNetworkError::NoNodes),
            1 => {
                let head = self.head();
                let root_node = self.graph.node_id(head);

                let mut slots = vec![];
                for h in self.graph.iter_crown(root_node) {
                    if let NetworkEdge::Slot(s) = self.graph[[&h]] {
                        slots.push(s);
                    }
                }

                Ok((&self.graph[root_node], root_node, slots))
            }
            _ => Err(TensorNetworkError::MoreThanOneNode),
        }
    }

    fn join_heads(
        self_flow: Flow,
        self_data: EdgeData<NetworkEdge<Aind>>,
        _other_flow: Flow,
        _other_data: EdgeData<NetworkEdge<Aind>>,
    ) -> (Flow, EdgeData<NetworkEdge<Aind>>) {
        (self_flow, self_data)
    }

    pub fn expr_tree(&self) -> SimpleTraversalTree {
        let _span = profile::span(Timer::ExprTree);
        let headgraph = self.expression_subgraph();

        let head = self.head();
        let root_node = self.graph.node_id(head);
        SimpleTraversalTree::depth_first_traverse(&self.graph, &headgraph, &root_node, None)
            .unwrap()
    }

    pub fn expression_subgraph(&self) -> SuBitGraph {
        self.graph
            .from_filter(|a| !matches!(a, NetworkEdge::Slot(_)))
    }

    pub fn expression_subgraph_ignoring<S>(&self, hidden: &S) -> SuBitGraph
    where
        S: SubSetLike<Base = SuBitGraph>,
    {
        let mut expression = self.expression_subgraph();
        for hedge in hidden.included_iter() {
            expression.sub(hedge);
        }
        expression
    }

    pub fn expr_tree_ignoring<S>(&self, hidden: &S) -> SimpleTraversalTree
    where
        S: SubSetLike<Base = SuBitGraph>,
    {
        let _span = profile::span(Timer::ExprTree);
        let headgraph = self.expression_subgraph_ignoring(hidden);

        let head = self.head();
        let root_node = self.graph.node_id(head);
        SimpleTraversalTree::depth_first_traverse(&self.graph, &headgraph, &root_node, None)
            .unwrap()
    }

    pub fn merge_ops(&mut self)
    where
        K: Clone + Debug,
    {
        #[derive(Clone, Copy, Debug, PartialEq, Eq)]
        enum MergeOp {
            Sum,
            Product,
        }

        fn find(parents: &mut [usize], index: usize) -> usize {
            let parent = parents[index];
            if parent == index {
                index
            } else {
                let root = find(parents, parent);
                parents[index] = root;
                root
            }
        }

        fn union(parents: &mut [usize], left: usize, right: usize) {
            let left_root = find(parents, left);
            let right_root = find(parents, right);
            if left_root != right_root {
                parents[right_root] = left_root;
            }
        }

        if profile::enabled() {
            eprintln!(
                "spenso_profile merge_ops.start nodes={} hedges={}",
                self.graph.n_nodes(),
                self.graph.n_hedges()
            );
        }
        // build a traversal over *all* internal edges
        let tt: SimpleTraversalTree<ParentChildStore<()>> = self.expr_tree().cast();
        if profile::enabled() {
            eprintln!("spenso_profile merge_ops.after_expr_tree");
        }
        let head = self.head();
        let root_node = self.graph.node_id(head);

        let mut op_nodes = Vec::new();
        let mut op_index = BTreeMap::new();

        // look for repeated ops nodes along a chain
        for nid in tt.iter_preorder_tree_nodes(&self.graph, root_node) {
            if let NetworkNode::Op(op) = &self.graph[nid] {
                let merge_op = match op {
                    NetworkOp::Sum => Some(MergeOp::Sum),
                    NetworkOp::Product => Some(MergeOp::Product),
                    _ => None,
                };

                if let Some(merge_op) = merge_op {
                    op_index.insert(nid, op_nodes.len());
                    op_nodes.push((nid, merge_op));
                }
            }
        }
        let mut parents = (0..op_nodes.len()).collect::<Vec<_>>();

        for index in 0..op_nodes.len() {
            let (node, merge_op) = op_nodes[index];
            for child in tt.iter_children(node, &self.graph) {
                if let Some(child_index) = op_index.get(&child)
                    && op_nodes[*child_index].1 == merge_op
                {
                    union(&mut parents, index, *child_index);
                }
            }
        }

        if profile::enabled() {
            let sums = op_nodes
                .iter()
                .filter(|(_, op)| *op == MergeOp::Sum)
                .count();
            let prods = op_nodes
                .iter()
                .filter(|(_, op)| *op == MergeOp::Product)
                .count();
            eprintln!("spenso_profile merge_ops.after_scan sums={sums} prods={prods}");
        }

        let mut to_del: SuBitGraph = self.graph.empty_subgraph();
        let mut groups: BTreeMap<usize, Vec<NodeIndex>> = BTreeMap::new();
        for (index, (node, _)) in op_nodes.iter().enumerate() {
            let root = find(&mut parents, index);
            groups.entry(root).or_default().push(*node);
        }

        for (root, nodes) in groups {
            if nodes.len() > 1 {
                let op = match op_nodes[root].1 {
                    MergeOp::Sum => NetworkOp::Sum,
                    MergeOp::Product => NetworkOp::Product,
                };
                let (_, sub) = self
                    .graph
                    .identify_nodes_without_self_edges::<SuBitGraph>(&nodes, NetworkNode::Op(op));
                to_del.union_with(&sub);
            }
        }
        if profile::enabled() {
            eprintln!(
                "spenso_profile merge_ops.after_components to_del={}",
                to_del.n_included()
            );
        }

        // println!("{}", self.graph.dot(&to_del));

        self.graph.forget_identification_history();
        if profile::enabled() {
            eprintln!("spenso_profile merge_ops.after_forget");
        }
        self.graph.delete_hedges(&to_del);
        if profile::enabled() {
            eprintln!(
                "spenso_profile merge_ops.end nodes={} hedges={}",
                self.graph.n_nodes(),
                self.graph.n_hedges()
            );
        }
    }
    pub fn simplify_identity_ops(&mut self) {}

    fn join_mut(
        &mut self,
        other: Self,
        matching_fn: impl Fn(
            Flow,
            EdgeData<&NetworkEdge<Aind>>,
            Flow,
            EdgeData<&NetworkEdge<Aind>>,
        ) -> bool,
        merge_fn: impl Fn(
            Flow,
            EdgeData<NetworkEdge<Aind>>,
            Flow,
            EdgeData<NetworkEdge<Aind>>,
        ) -> (Flow, EdgeData<NetworkEdge<Aind>>),
    ) -> Result<(), HedgeGraphError> {
        let _span = profile::span(Timer::GraphJoin);
        profile::bump(Counter::GraphJoin, 1);
        let join_profile = profile::enabled().then(|| {
            (
                Instant::now(),
                self.graph.n_nodes(),
                self.graph.n_hedges(),
                other.graph.n_nodes(),
                other.graph.n_hedges(),
            )
        });
        self.graph.join_mut(other.graph, matching_fn, merge_fn)?;
        if let Some((start, self_nodes, self_hedges, other_nodes, other_hedges)) = join_profile {
            let elapsed = start.elapsed();
            if elapsed.as_millis() >= 50 {
                eprintln!(
                    "spenso_profile slow graph.join elapsed={elapsed:.3?} self_nodes={self_nodes} self_hedges={self_hedges} other_nodes={other_nodes} other_hedges={other_hedges} out_nodes={} out_hedges={}",
                    self.graph.n_nodes(),
                    self.graph.n_hedges(),
                );
            }
        }
        self.slot_order.extend(other.slot_order);
        // self.uncontracted.join_mut(other.uncontracted);
        Ok(())
    }

    pub fn pow(self, pow: i8) -> Self {
        let dangling = self.dangling_indices();
        let mut pow = NetworkGraph::pow_graph(pow, &dangling);

        pow.join_mut(
            self,
            NetworkGraph::<K, FK, Aind>::add_match,
            NetworkGraph::<K, FK, Aind>::join_heads,
        )
        .unwrap();

        pow
    }

    pub fn function(self, name: FK) -> Self {
        let mut fun = NetworkGraph::fun_graph(name);

        fun.join_mut(
            self,
            NetworkGraph::<K, FK, Aind>::match_heads,
            NetworkGraph::<K, FK, Aind>::join_heads,
        )
        .unwrap();

        fun
    }
}

pub trait NMul<Rhs = Self> {
    type Output;

    fn n_mul<I: IntoIterator<Item = Rhs>>(self, iter: I) -> Self::Output;
}

impl<K: Debug, FK: Debug, Aind: AbsInd> NMul for NetworkGraph<K, FK, Aind> {
    type Output = NetworkGraph<K, FK, Aind>;
    fn n_mul<I: IntoIterator<Item = Self>>(self, iter: I) -> Self::Output {
        let _span = profile::span(Timer::GraphNMul);
        profile::bump(Counter::GraphNMul, 1);
        let all = iter.into_iter().collect::<Vec<_>>();
        if profile::enabled() && (all.len() > 100 || self.graph.n_nodes() > 100) {
            eprintln!(
                "spenso_profile graph.n_mul.start inputs={} self_nodes={} self_hedges={}",
                all.len() + 1,
                self.graph.n_nodes(),
                self.graph.n_hedges(),
            );
        }
        let mut mul = Self::mul_graph(all.len() + 1);

        mul.join_mut(
            self,
            NetworkGraph::<K, FK, Aind>::match_heads,
            NetworkGraph::<K, FK, Aind>::join_heads,
        )
        .unwrap();

        for rhs in all {
            mul.join_mut(
                rhs,
                NetworkGraph::<K, FK, Aind>::prod_match,
                NetworkGraph::<K, FK, Aind>::join_heads,
            )
            .unwrap();
        }

        mul
    }
}

// impl<S: NMul<R>, R> Mul<R> for S {grrr
//     type Output = S::Output;
//     fn mul(self, rhs: R) -> Self::Output {
//         self.n_mul([rhs])
//     }
// }

pub trait NAdd<Rhs = Self> {
    type Output;

    fn n_add<I: IntoIterator<Item = Rhs>>(self, iter: I) -> Self::Output;
}

impl<K: Debug, FK: Debug, Aind: AbsInd> NAdd for NetworkGraph<K, FK, Aind> {
    type Output = NetworkGraph<K, FK, Aind>;

    fn n_add<I: IntoIterator<Item = Self>>(self, iter: I) -> Self::Output {
        let _span = profile::span(Timer::GraphNAdd);
        profile::bump(Counter::GraphNAdd, 1);
        let all = iter.into_iter().collect::<Vec<_>>();
        let slots = self.dangling_indices();
        if profile::enabled() && (all.len() > 100 || self.graph.n_nodes() > 100) {
            eprintln!(
                "spenso_profile graph.n_add.start inputs={} slots={} self_nodes={} self_hedges={}",
                all.len() + 1,
                slots.len(),
                self.graph.n_nodes(),
                self.graph.n_hedges(),
            );
        }

        let mut add = Self::add_graph(all.len() + 1, &slots);

        add.join_mut(
            self,
            NetworkGraph::<K, FK, Aind>::add_match,
            NetworkGraph::<K, FK, Aind>::join_heads,
        )
        .unwrap();

        for rhs in all {
            debug_assert!(
                slots.len() == rhs.n_dangling(),
                "Mismatched dangling edges in sum, Trying to add {} to {}",
                add.dot_simple(),
                rhs.dot_simple()
            );

            add.join_mut(
                rhs,
                NetworkGraph::<K, FK, Aind>::add_match,
                NetworkGraph::<K, FK, Aind>::join_heads,
            )
            .unwrap();
        }
        debug_assert!(
            slots.len() == add.n_dangling(),
            "addition with non matching dangling: \n{}",
            add.dot_impl(
                |_| "".to_string(),
                |_| "".to_string(),
                |_| "".to_string(),
                |_| "".to_string()
            )
        );

        add
    }
}

impl<K: Debug, FK: Debug, Aind: AbsInd> MulAssign for NetworkGraph<K, FK, Aind> {
    fn mul_assign(&mut self, rhs: Self) {
        let mul = Self::mul_graph(2);

        self.join_mut(
            mul,
            NetworkGraph::<K, FK, Aind>::match_heads,
            NetworkGraph::<K, FK, Aind>::join_heads,
        )
        .unwrap();

        self.join_mut(
            rhs,
            NetworkGraph::<K, FK, Aind>::prod_match,
            NetworkGraph::<K, FK, Aind>::join_heads,
        )
        .unwrap();
    }
}

impl<K: Clone + Debug, FK: Clone + Debug, Aind: AbsInd> MulAssign<&NetworkGraph<K, FK, Aind>>
    for NetworkGraph<K, FK, Aind>
{
    fn mul_assign(&mut self, rhs: &Self) {
        let mul = NetworkGraph::mul_graph(2);

        self.join_mut(
            mul,
            NetworkGraph::<K, FK, Aind>::match_heads,
            NetworkGraph::<K, FK, Aind>::join_heads,
        )
        .unwrap();

        self.join_mut(
            rhs.clone(),
            NetworkGraph::<K, FK, Aind>::prod_match,
            NetworkGraph::<K, FK, Aind>::join_heads,
        )
        .unwrap();
    }
}

impl<K: Debug, FK: Debug, Aind: AbsInd> Mul for NetworkGraph<K, FK, Aind> {
    type Output = NetworkGraph<K, FK, Aind>;
    fn mul(self, rhs: Self) -> Self::Output {
        let mut mul = Self::mul_graph(2);

        mul.join_mut(
            self,
            NetworkGraph::<K, FK, Aind>::match_heads,
            NetworkGraph::<K, FK, Aind>::join_heads,
        )
        .unwrap();

        mul.join_mut(
            rhs,
            NetworkGraph::<K, FK, Aind>::prod_match,
            NetworkGraph::<K, FK, Aind>::join_heads,
        )
        .unwrap();

        mul
    }
}
impl<K: Clone + Debug, FK: Clone + Debug, Aind: AbsInd> Mul<&NetworkGraph<K, FK, Aind>>
    for NetworkGraph<K, FK, Aind>
{
    type Output = NetworkGraph<K, FK, Aind>;
    fn mul(self, rhs: &Self) -> Self::Output {
        self * rhs.clone()
    }
}

impl<'b, K: Clone + Debug, FK: Clone + Debug, Aind: AbsInd> Mul<&'b NetworkGraph<K, FK, Aind>>
    for &NetworkGraph<K, FK, Aind>
{
    type Output = NetworkGraph<K, FK, Aind>;
    fn mul(self, rhs: &'b NetworkGraph<K, FK, Aind>) -> Self::Output {
        self.clone() * rhs
    }
}

impl<K: Clone + Debug, FK: Clone + Debug, Aind: AbsInd> Mul<NetworkGraph<K, FK, Aind>>
    for &NetworkGraph<K, FK, Aind>
{
    type Output = NetworkGraph<K, FK, Aind>;
    fn mul(self, rhs: NetworkGraph<K, FK, Aind>) -> Self::Output {
        rhs * self
    }
}

impl<K: Debug, FK: Debug, Aind: AbsInd> AddAssign for NetworkGraph<K, FK, Aind> {
    fn add_assign(&mut self, rhs: Self) {
        let slots = self.dangling_indices();
        debug_assert!(slots.len() == rhs.n_dangling());

        let add = Self::add_graph(2, &slots);

        self.join_mut(
            add,
            NetworkGraph::<K, FK, Aind>::add_match,
            NetworkGraph::<K, FK, Aind>::join_heads,
        )
        .unwrap();

        self.join_mut(
            rhs,
            NetworkGraph::<K, FK, Aind>::add_match,
            NetworkGraph::<K, FK, Aind>::join_heads,
        )
        .unwrap();
    }
}

impl<K: Clone + Debug, FK: Clone + Debug, Aind: AbsInd> AddAssign<&NetworkGraph<K, FK, Aind>>
    for NetworkGraph<K, FK, Aind>
{
    fn add_assign(&mut self, rhs: &Self) {
        let slots = self.dangling_indices();
        debug_assert!(slots.len() == rhs.n_dangling());

        let add = Self::add_graph(2, &slots);

        self.join_mut(
            add,
            NetworkGraph::<K, FK, Aind>::add_match,
            NetworkGraph::<K, FK, Aind>::join_heads,
        )
        .unwrap();

        self.join_mut(
            rhs.clone(),
            NetworkGraph::<K, FK, Aind>::add_match,
            NetworkGraph::<K, FK, Aind>::join_heads,
        )
        .unwrap();
    }
}

impl<K: Debug, FK: Debug, Aind: AbsInd> Add for NetworkGraph<K, FK, Aind> {
    type Output = NetworkGraph<K, FK, Aind>;
    fn add(self, rhs: Self) -> Self::Output {
        let slots = self.dangling_indices();
        debug_assert!(slots.len() == rhs.n_dangling());

        let mut add = Self::add_graph(2, &slots);

        add.join_mut(
            self,
            NetworkGraph::<K, FK, Aind>::add_match,
            NetworkGraph::<K, FK, Aind>::join_heads,
        )
        .unwrap();

        add.join_mut(
            rhs,
            NetworkGraph::<K, FK, Aind>::add_match,
            NetworkGraph::<K, FK, Aind>::join_heads,
        )
        .unwrap();

        debug_assert!(slots.len() == add.n_dangling());

        add
    }
}

impl<K: Clone + Debug, FK: Clone + Debug, Aind: AbsInd> Add<&NetworkGraph<K, FK, Aind>>
    for NetworkGraph<K, FK, Aind>
{
    type Output = NetworkGraph<K, FK, Aind>;
    fn add(self, rhs: &Self) -> Self::Output {
        self + rhs.clone()
    }
}

impl<'b, K: Clone + Debug, FK: Clone + Debug, Aind: AbsInd> Add<&'b NetworkGraph<K, FK, Aind>>
    for &NetworkGraph<K, FK, Aind>
{
    type Output = NetworkGraph<K, FK, Aind>;
    fn add(self, rhs: &'b NetworkGraph<K, FK, Aind>) -> Self::Output {
        self.clone() + rhs
    }
}

impl<K: Clone + Debug, FK: Clone + Debug, Aind: AbsInd> Add<NetworkGraph<K, FK, Aind>>
    for &NetworkGraph<K, FK, Aind>
{
    type Output = NetworkGraph<K, FK, Aind>;
    fn add(self, rhs: NetworkGraph<K, FK, Aind>) -> Self::Output {
        rhs + self
    }
}

impl<K: Clone + Debug, FK: Clone + Debug, Aind: AbsInd> Neg for NetworkGraph<K, FK, Aind> {
    type Output = NetworkGraph<K, FK, Aind>;
    fn neg(self) -> Self::Output {
        let mut neg = Self::neg_graph();

        neg.join_mut(
            self,
            NetworkGraph::<K, FK, Aind>::match_heads,
            NetworkGraph::<K, FK, Aind>::join_heads,
        )
        .unwrap();

        neg
    }
}

impl<K: Clone + Debug, FK: Clone + Debug, Aind: AbsInd> Neg for &NetworkGraph<K, FK, Aind> {
    type Output = NetworkGraph<K, FK, Aind>;
    fn neg(self) -> Self::Output {
        let mut neg = NetworkGraph::neg_graph();

        neg.join_mut(
            self.clone(),
            NetworkGraph::<K, FK, Aind>::match_heads,
            NetworkGraph::<K, FK, Aind>::join_heads,
        )
        .unwrap();

        neg
    }
}

impl<K: Clone + Debug, FK: Clone + Debug, Aind: AbsInd> SubAssign for NetworkGraph<K, FK, Aind> {
    fn sub_assign(&mut self, rhs: NetworkGraph<K, FK, Aind>) {
        *self += -rhs;
    }
}

impl<K: Clone + Debug, FK: Clone + Debug, Aind: AbsInd> SubAssign<&NetworkGraph<K, FK, Aind>>
    for NetworkGraph<K, FK, Aind>
{
    fn sub_assign(&mut self, rhs: &NetworkGraph<K, FK, Aind>) {
        *self += -rhs;
    }
}

impl<K: Clone + Debug, FK: Clone + Debug, Aind: AbsInd> Sub for NetworkGraph<K, FK, Aind> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + -rhs
    }
}
impl<K: Clone + Debug, FK: Clone + Debug, Aind: AbsInd> Sub<&NetworkGraph<K, FK, Aind>>
    for NetworkGraph<K, FK, Aind>
{
    type Output = Self;

    fn sub(self, rhs: &Self) -> Self::Output {
        self + -rhs
    }
}
impl<'b, K: Clone + Debug, FK: Clone + Debug, Aind: AbsInd> Sub<&'b NetworkGraph<K, FK, Aind>>
    for &NetworkGraph<K, FK, Aind>
{
    type Output = NetworkGraph<K, FK, Aind>;

    fn sub(self, rhs: &'b NetworkGraph<K, FK, Aind>) -> Self::Output {
        -rhs + self
    }
}
impl<K: Clone + Debug, FK: Clone + Debug, Aind: AbsInd> Sub<NetworkGraph<K, FK, Aind>>
    for &NetworkGraph<K, FK, Aind>
{
    type Output = NetworkGraph<K, FK, Aind>;

    fn sub(self, rhs: NetworkGraph<K, FK, Aind>) -> Self::Output {
        -rhs + self
    }
}

// pub trait TestContext {
//     fn get_statemap(&mut self) -> &mut StateMap;

//     fn get_model(&mut self) -> &mut Model;
// }

// pub struct Model {}

// pub struct Expr {}

// impl<T: TestContext> Decode<T> for Expr {
//     fn decode<D: bincode::de::Decoder<Context = T>>(
//         decoder: &mut D,
//     ) -> Result<Self, bincode::error::DecodeError> {
//         let context = decoder.context().get_statemap();
//         Ok(Expr {})
//     }
// }

// pub struct Test {
//     atom: Expr,
// }

// impl<T: TestContext> Decode<T> for Test {
//     fn decode<D: bincode::de::Decoder<Context = T>>(
//         decoder: &mut D,
//     ) -> Result<Self, bincode::error::DecodeError> {
//         let atom = Expr::decode(decoder)?;
//         Ok(Test { atom })
//     }
// }

#[cfg(test)]
pub mod test {

    use linnet::{
        half_edge::subgraph::{ModifySubSet, SuBitGraph, SubSetLike, SubSetOps},
        tree::child_vec::ChildVecStore,
    };

    use crate::{
        network::graph::NetworkLeaf,
        structure::{
            OrderedStructure, PermutedStructure,
            representation::{Euclidean, Lorentz, Minkowski, RepName},
            slot::IsAbstractSlot,
        },
    };

    use super::{NetworkGraph, NetworkNode, NetworkOp};

    #[test]
    fn addition() {
        let one = NetworkGraph::<i8>::one();
        let zero = NetworkGraph::<i8>::zero();

        let s = NetworkGraph::<i8>::scalar(2);

        let t = NetworkGraph::<i8>::tensor(
            &PermutedStructure::<OrderedStructure>::from_iter([
                Minkowski {}.new_slot(1, 2),
                Minkowski {}.new_slot(2, 2),
            ])
            .structure,
            NetworkLeaf::LocalTensor(1),
        );

        let t2 = NetworkGraph::<i8>::tensor(
            &PermutedStructure::<OrderedStructure>::from_iter([
                Minkowski {}.new_slot(1, 2),
                Minkowski {}.new_slot(2, 2),
            ])
            .structure,
            NetworkLeaf::LocalTensor(2),
        );

        let t3 = NetworkGraph::<i8>::tensor(
            &PermutedStructure::<OrderedStructure>::from_iter([
                Lorentz {}.new_slot(1, 2).to_lib(),
                Euclidean {}.new_slot(2, 2).to_lib(),
            ])
            .structure,
            NetworkLeaf::LocalTensor(2),
        );

        let t3b = NetworkGraph::<i8>::tensor(
            &PermutedStructure::<OrderedStructure>::from_iter([
                Lorentz {}.dual().new_slot(1, 2).to_lib(),
                Euclidean {}.new_slot(2, 1).to_lib(),
            ])
            .structure,
            NetworkLeaf::LocalTensor(2),
        );
        let s2 = NetworkGraph::<i8>::scalar(3);

        let mut expr = t3b * (&t + &t2) * s2 * (&s + &zero) * (t3 * (t + t2) * s * (one + zero));

        // println!("{}", expr.dot());
        expr.merge_ops();
        println!("{}", expr.dot());

        // expr.extract_next_ready_op();
        if let Some((a, _)) = expr.extract_next_ready_op() {
            println!("{}", a.dot());
        }
    }

    #[test]
    fn cached_expr_children_match_traversal_tree_after_root_alignment() {
        let scalar = NetworkGraph::<i8>::scalar(2);
        let scalar_b = NetworkGraph::<i8>::scalar(3);
        let tensor_a = NetworkGraph::<i8>::tensor(
            &PermutedStructure::<OrderedStructure>::from_iter([
                Minkowski {}.new_slot(1, 2),
                Minkowski {}.new_slot(2, 2),
            ])
            .structure,
            NetworkLeaf::LocalTensor(1),
        );
        let tensor_b = NetworkGraph::<i8>::tensor(
            &PermutedStructure::<OrderedStructure>::from_iter([
                Minkowski {}.new_slot(1, 2),
                Minkowski {}.new_slot(2, 2),
            ])
            .structure,
            NetworkLeaf::LocalTensor(2),
        );

        let sum = tensor_a + tensor_b;
        let mut expr = (sum.clone() * scalar) + (sum * scalar_b);
        expr.merge_ops();

        expr.cache_expr_tree_roots();

        let tt = expr.expr_tree().cast::<ChildVecStore<()>>();
        let root_node = expr.graph.node_id(expr.head());
        for node in tt.iter_preorder_tree_nodes(&expr.graph, root_node) {
            let mut expected = tt.iter_children(node, &expr.graph).collect::<Vec<_>>();
            expected.sort();
            expected.dedup();

            assert_eq!(expr.cached_expr_children(node), expected);
        }
    }

    #[test]
    fn hidden_expression_edges_are_ignored_by_cached_children() {
        let scalar = NetworkGraph::<i8>::scalar(2);
        let scalar_b = NetworkGraph::<i8>::scalar(3);
        let tensor_a = NetworkGraph::<i8>::tensor(
            &PermutedStructure::<OrderedStructure>::from_iter([
                Minkowski {}.new_slot(1, 2),
                Minkowski {}.new_slot(2, 2),
            ])
            .structure,
            NetworkLeaf::LocalTensor(1),
        );
        let tensor_b = NetworkGraph::<i8>::tensor(
            &PermutedStructure::<OrderedStructure>::from_iter([
                Minkowski {}.new_slot(1, 2),
                Minkowski {}.new_slot(2, 2),
            ])
            .structure,
            NetworkLeaf::LocalTensor(2),
        );

        let sum = tensor_a + tensor_b;
        let mut expr = (sum.clone() * scalar) + (sum * scalar_b);
        expr.merge_ops();
        expr.cache_expr_tree_roots();

        let tt = expr.expr_tree().cast::<ChildVecStore<()>>();
        let root_node = expr.graph.node_id(expr.head());
        let (parent, child) = tt
            .iter_preorder_tree_nodes(&expr.graph, root_node)
            .find_map(|node| {
                tt.iter_children(node, &expr.graph)
                    .next()
                    .map(|child| (node, child))
            })
            .unwrap();
        let child_root = tt.root_hedge(child);
        let parent_side = expr.graph.inv(child_root);

        let mut hidden: SuBitGraph = expr.graph.empty_subgraph();
        hidden.add(child_root);
        hidden.add(parent_side);

        expr.cache_expr_tree_roots_ignoring(&hidden);

        let children = expr.cached_expr_children_ignoring(parent, &hidden);
        assert!(!children.contains(&child));
    }

    #[test]
    fn ready_operation_ref_describes_ready_leaf_subgraph_without_extracting() {
        let scalar = NetworkGraph::<i8>::scalar(2);
        let scalar_b = NetworkGraph::<i8>::scalar(3);
        let tensor_a = NetworkGraph::<i8>::tensor(
            &PermutedStructure::<OrderedStructure>::from_iter([
                Minkowski {}.new_slot(1, 2),
                Minkowski {}.new_slot(2, 2),
            ])
            .structure,
            NetworkLeaf::LocalTensor(1),
        );
        let tensor_b = NetworkGraph::<i8>::tensor(
            &PermutedStructure::<OrderedStructure>::from_iter([
                Minkowski {}.new_slot(1, 2),
                Minkowski {}.new_slot(2, 2),
            ])
            .structure,
            NetworkLeaf::LocalTensor(2),
        );

        let sum = tensor_a + tensor_b;
        let mut expr = (sum.clone() * scalar) + (sum * scalar_b);
        expr.merge_ops();
        expr.cache_expr_tree_roots();

        let node_count = expr.n_nodes();
        let hedge_count = expr.graph.n_hedges();
        let op_ref = expr.ready_operation_ref().unwrap();

        assert!(op_ref.leaf_count() > 0);
        assert!(
            op_ref
                .children()
                .iter()
                .all(|child| matches!(&expr.graph[*child], NetworkNode::Leaf(_)))
        );

        let expected = expr.operation_subgraph(op_ref.op_node(), op_ref.children());
        assert_eq!(op_ref.subgraph(), &expected);
        assert!(op_ref.subgraph().n_included() >= op_ref.leaf_count());
        assert_eq!(op_ref.graph().n_nodes(), node_count);
        assert_eq!(expr.n_nodes(), node_count);
        assert_eq!(expr.graph.n_hedges(), hedge_count);
    }

    #[test]
    fn ready_operation_refs_batch_non_overlapping_leaf_operations() {
        let scalar = NetworkGraph::<i8>::scalar(2);
        let scalar_b = NetworkGraph::<i8>::scalar(3);
        let tensor_a = NetworkGraph::<i8>::tensor(
            &PermutedStructure::<OrderedStructure>::from_iter([
                Minkowski {}.new_slot(1, 2),
                Minkowski {}.new_slot(2, 2),
            ])
            .structure,
            NetworkLeaf::LocalTensor(1),
        );
        let tensor_b = NetworkGraph::<i8>::tensor(
            &PermutedStructure::<OrderedStructure>::from_iter([
                Minkowski {}.new_slot(1, 2),
                Minkowski {}.new_slot(2, 2),
            ])
            .structure,
            NetworkLeaf::LocalTensor(2),
        );

        let mut expr = (tensor_a + tensor_b) * (scalar + scalar_b);
        expr.merge_ops();
        expr.cache_expr_tree_roots();

        let ready = expr.ready_operation_refs();
        assert_eq!(ready.len(), 2);

        let mut used: SuBitGraph = expr.graph.empty_subgraph();
        for op_ref in &ready {
            assert!(matches!(op_ref.op(), NetworkOp::Sum));
            assert!(op_ref.subgraph().empty_intersection(&used));
            used.union_with(op_ref.subgraph());
            assert!(
                op_ref
                    .children()
                    .iter()
                    .all(|child| matches!(&expr.graph[*child], NetworkNode::Leaf(_)))
            );
        }
    }

    #[test]
    fn readiness_diagnostics_match_ready_extraction_boundary() {
        let scalar = NetworkGraph::<i8>::scalar(2);
        let scalar_b = NetworkGraph::<i8>::scalar(3);
        let tensor_a = NetworkGraph::<i8>::tensor(
            &PermutedStructure::<OrderedStructure>::from_iter([
                Minkowski {}.new_slot(1, 2),
                Minkowski {}.new_slot(2, 2),
            ])
            .structure,
            NetworkLeaf::LocalTensor(1),
        );
        let tensor_b = NetworkGraph::<i8>::tensor(
            &PermutedStructure::<OrderedStructure>::from_iter([
                Minkowski {}.new_slot(1, 2),
                Minkowski {}.new_slot(2, 2),
            ])
            .structure,
            NetworkLeaf::LocalTensor(2),
        );

        let mut expr = (tensor_a + tensor_b) * (scalar + scalar_b);
        expr.merge_ops();
        expr.cache_expr_tree_roots();

        let diagnostics = expr.operation_readiness_diagnostics();
        assert_eq!(
            diagnostics.cached_expression_node_count,
            expr.cached_expr_preorder_nodes().len()
        );
        assert_eq!(diagnostics.ready_operation_count, 2);
        assert_eq!(diagnostics.batched_operation_count, 2);
        assert!(diagnostics.batched_subgraph_hedge_count > 0);

        let op_ref = expr.ready_operation_ref().unwrap();
        let mut extracted_expr = expr.clone();
        let (extracted, op) = extracted_expr.extract_next_ready_op().unwrap();
        let extracted_leaf_count = extracted
            .graph
            .iter_nodes()
            .filter(|(_, _, node)| matches!(node, NetworkNode::Leaf(_)))
            .count();

        assert_eq!(op_ref.op(), &op);
        assert_eq!(op_ref.leaf_count(), extracted_leaf_count);
        assert_eq!(extracted.n_nodes(), op_ref.leaf_count() + 1);
    }

    #[test]
    fn ready_operation_batch_carries_parallel_scheduler_subgraph() {
        let scalar = NetworkGraph::<i8>::scalar(2);
        let scalar_b = NetworkGraph::<i8>::scalar(3);
        let tensor_a = NetworkGraph::<i8>::tensor(
            &PermutedStructure::<OrderedStructure>::from_iter([
                Minkowski {}.new_slot(1, 2),
                Minkowski {}.new_slot(2, 2),
            ])
            .structure,
            NetworkLeaf::LocalTensor(1),
        );
        let tensor_b = NetworkGraph::<i8>::tensor(
            &PermutedStructure::<OrderedStructure>::from_iter([
                Minkowski {}.new_slot(1, 2),
                Minkowski {}.new_slot(2, 2),
            ])
            .structure,
            NetworkLeaf::LocalTensor(2),
        );

        let mut expr = (tensor_a + tensor_b) * (scalar + scalar_b);
        expr.merge_ops();
        expr.cache_expr_tree_roots();

        let batch = expr.ready_operation_batch();
        assert_eq!(batch.len(), 2);
        assert!(!batch.is_empty());
        assert_eq!(batch.operations().len(), batch.len());

        let mut union: SuBitGraph = expr.graph.empty_subgraph();
        for op_ref in batch.iter() {
            assert!(op_ref.subgraph().empty_intersection(&union));
            union.union_with(op_ref.subgraph());
        }

        assert_eq!(batch.subgraph(), &union);
    }

    #[test]
    fn ready_ref_extraction_matches_legacy_ready_extraction_boundary() {
        let scalar = NetworkGraph::<i8>::scalar(2);
        let scalar_b = NetworkGraph::<i8>::scalar(3);
        let tensor_a = NetworkGraph::<i8>::tensor(
            &PermutedStructure::<OrderedStructure>::from_iter([
                Minkowski {}.new_slot(1, 2),
                Minkowski {}.new_slot(2, 2),
            ])
            .structure,
            NetworkLeaf::LocalTensor(1),
        );
        let tensor_b = NetworkGraph::<i8>::tensor(
            &PermutedStructure::<OrderedStructure>::from_iter([
                Minkowski {}.new_slot(1, 2),
                Minkowski {}.new_slot(2, 2),
            ])
            .structure,
            NetworkLeaf::LocalTensor(2),
        );

        let mut legacy_expr = (tensor_a + tensor_b) * (scalar + scalar_b);
        legacy_expr.merge_ops();
        let mut ref_expr = legacy_expr.clone();

        let (legacy_extracted, legacy_op) = legacy_expr.extract_next_ready_op().unwrap();
        let (ref_extracted, ref_op) = ref_expr.extract_next_ready_ref_op().unwrap();

        let legacy_leaf_count = legacy_extracted
            .graph
            .iter_nodes()
            .filter(|(_, _, node)| matches!(node, NetworkNode::Leaf(_)))
            .count();
        let ref_leaf_count = ref_extracted
            .graph
            .iter_nodes()
            .filter(|(_, _, node)| matches!(node, NetworkNode::Leaf(_)))
            .count();

        assert_eq!(legacy_op, ref_op);
        assert_eq!(legacy_leaf_count, ref_leaf_count);
        assert_eq!(legacy_extracted.n_nodes(), ref_extracted.n_nodes());
    }
}
