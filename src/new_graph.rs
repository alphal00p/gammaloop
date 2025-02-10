use std::path::{Path, PathBuf};

use bitvec::vec::BitVec;
use color_eyre::Result;
use hyperdual::Num;
use linnet::half_edge::{
    subgraph::{
        self, cycle::SignedCycle, Inclusion, InternalSubGraph, OrientedCut, SubGraph, SubGraphOps,
    },
    EdgeData, EdgeId, Hedge, HedgeGraph, HedgeVec, Involution, Orientation, Parent, TraversalTree,
};
use symbolica::atom::{representation::InlineNum, Atom};

use crate::{
    graph::{BareGraph, DerivedGraphData, Edge, LoopMomentumBasis, Vertex},
    momentum::{SignOrZero, Signature},
    numerator::{NumeratorState, PythonState, UnInit},
    ProcessSettings,
};

pub struct Graph<S: NumeratorState = PythonState> {
    multiplicity: Atom,
    underlying: HedgeGraph<Edge, Vertex>,
    loop_momentum_basis: HedgeLMB,
    derived_data: DerivedGraphData<S>,
}

impl<S: NumeratorState> Graph<S> {
    pub fn forget_type(self) -> Graph<PythonState> {
        Graph {
            multiplicity: self.multiplicity,
            underlying: self.underlying,
            loop_momentum_basis: self.loop_momentum_basis,
            derived_data: self.derived_data.forget_type(),
        }
    }
}

pub trait FeynmanGraph {
    fn new_lmb(&self) -> HedgeLMB;
    fn num_virtual_edges(&self, subgraph: BitVec) -> usize;
}

pub struct HedgeLMB {
    tree: TraversalTree,
    lmb_basis: Vec<Hedge>,
}

impl FeynmanGraph for HedgeGraph<Edge, Vertex> {
    fn new_lmb(&self) -> HedgeLMB {
        // root node should contain a dangling (external edge), that will be the dependent external
        let root = self
            .iter_nodes()
            .find(|(a, _)| {
                self.iter_egdes(*a)
                    .any(|(e, _)| matches!(e, EdgeId::Unpaired { .. }))
            })
            .unwrap_or(self.iter_nodes().next().unwrap())
            .0;

        let tree = TraversalTree::dfs(self, &self.full_filter(), root, None);

        let mut leaves = tree.leaf_edges();

        let mut external_leaves = self.empty_filter();
        let mut externals = vec![];

        for i in leaves.included_iter() {
            if self.involution.is_identity(i) {
                external_leaves.set(i.0, true);
                externals.push(i);
            }
        }

        let mut ext_signatures: HedgeVec<Signature> = self.new_derived_edge_data_empty();

        let empty_signature = Signature::from_iter(externals.iter().map(|e| SignOrZero::Zero));
        for (i, &h) in externals.iter().enumerate() {
            let mut signature = empty_signature.clone();
            signature.0[i] = SignOrZero::Plus;
            ext_signatures[h] = Some(signature);
        }

        let mut current_leaf_nodes = tree.leaf_nodes(&self);

        while let Some(leaf_node) = current_leaf_nodes.pop() {
            let hairs = &self.hairs_from_id(leaf_node).hairs;
            let mut root_pointer = None;
            let mut root_signature = empty_signature.clone();

            for h in hairs.included_iter() {
                match tree.parent(h) {
                    Parent::Root => {}
                    Parent::Hedge { hedge_to_root, .. } => {
                        if *hedge_to_root == h {
                            root_pointer = Some(h);
                            if self
                                .involved_node_hairs(h)
                                .unwrap()
                                .hairs
                                .included_iter()
                                .all(|a| a != h && ext_signatures.is_set(a))
                            {
                                current_leaf_nodes.push(self.involved_node_id(h).unwrap());
                            }
                        } else {
                            root_signature.sum(ext_signatures[h].as_ref().unwrap());
                        }
                    }
                    Parent::Unset => {}
                }
            }

            if let Some(root_pointer) = root_pointer {
                ext_signatures[root_pointer] = Some(root_signature);
            }
        }

        let tree_complement = tree.tree.complement(self);

        let mut cycle_basis = Vec::new();
        let mut lmb_basis = Vec::new();

        for (e, d) in self.iter_egdes(&tree_complement) {
            if let EdgeId::Paired { source, sink } = e {
                lmb_basis.push(source);
                cycle_basis.push(
                    SignedCycle::from_cycle(tree.cycle(source).unwrap(), source, self).unwrap(),
                );
            }
        }

        let signatures = self
            .involution
            .map_data_ref(|a| (), &|e| ())
            .map_edge_data(|e, d| {
                let e = match e {
                    EdgeId::Paired { source, sink } => {
                        let mut internal = vec![];
                        for (i, c) in cycle_basis.iter().enumerate() {
                            if c.filter.includes(&source) {
                                internal.push(SignOrZero::Plus);
                            } else if c.filter.includes(&sink) {
                                internal.push(SignOrZero::Minus);
                            } else {
                                internal.push(SignOrZero::Zero);
                            }
                        }

                        let internal_signature = Signature::from_iter(internal);

                        // return EdgeData::new(Signature::from_iter(iter), orientation)
                    }
                    EdgeId::Unpaired { hedge, flow } => {}
                    _ => {}
                };
                d
            });

        HedgeLMB { tree, lmb_basis }
    }

    fn num_virtual_edges(&self, subgraph: BitVec) -> usize {
        let internal_subgraph = InternalSubGraph::cleaned_filter_pessimist(subgraph, self);
        self.count_internal_edges(&internal_subgraph)
    }
}

impl Graph<UnInit> {
    pub fn new(multiplicity: Atom, underlying: HedgeGraph<Edge, Vertex>) -> Self {
        Self {
            multiplicity,
            loop_momentum_basis: underlying.new_lmb(),
            underlying,
            derived_data: DerivedGraphData::new_empty(),
        }
    }
}
