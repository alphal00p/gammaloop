use std::{
    collections::VecDeque,
    path::{Path, PathBuf},
};

use ahash::HashMap;
use bincode::{Decode, Encode};
use bitvec::vec::BitVec;
use color_eyre::{Report, Result};
use hyperdual::Num;
use itertools::Itertools;
use linnet::half_edge::{
    involution::EdgeIndex,
    subgraph::{
        self, cycle::SignedCycle, Inclusion, InternalSubGraph, OrientedCut, SubGraph, SubGraphOps,
    },
    EdgeData, EdgeId, Hedge, HedgeGraph, HedgeVec, Involution, Orientation, Parent, TraversalTree,
};
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{representation::InlineNum, Atom, AtomCore},
    id::Pattern,
};
use typed_index_collections::TiVec;

use crate::{
    gammaloop_integrand::BareSample,
    graph::{BareGraph, DerivedGraphData, Edge, EdgeType, Vertex},
    momentum::{FourMomentum, SignOrZero, Signature, ThreeMomentum},
    momentum_sample::LoopIndex,
    numerator::{NumeratorState, PythonState, UnInit},
    signature::{ExternalSignature, LoopExtSignature, LoopSignature},
    utils::{FloatLike, F},
    ProcessSettings,
};

pub struct Graph<S: NumeratorState = PythonState> {
    pub multiplicity: Atom,
    pub underlying: HedgeGraph<Edge, Vertex>,
    pub loop_momentum_basis: HedgeLMB,
    pub derived_data: DerivedGraphData<S>,
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
            let hairs = &self.hairs_from_id(leaf_no        if start == end {
                return Some(vec![]);
            }
    
            // Initialize BFS
            let mut queue = VecDeque::new();
            let mut visited: HashMap<usize, Option<(usize, usize, bool)>> = HashMap::new();
    
            queue.push_back(start);
            visited.insert(start, None);
    
            // Perform BFS
            while let Some(u) = queue.pop_front() {
                if u == end {
                    break;
                }
                if let Some(neighbors) = adjacency_list.get(&u) {
                    for &(v, edge_index, is_flipped) in neighbors {
                        #[allow(clippy::map_entry)]
                        if !visited.contains_key(&v) {
                            visited.insert(v, Some((u, edge_index, is_flipped)));
                            queue.push_back(v);
                        }
                    }
                }
            }
    
            // Reconstruct the path if end is reached
            if !visited.contains_key(&end) {
                return None;
            }
    
            let mut path = Vec::new();
            let mut current = end;
    
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

#[derive(Debug, Clone, Serialize, Deserialize /*, Encode, Decode*/)]

pub struct LoopMomentumBasis {
    pub basis: TiVec<LoopIndex, EdgeIndex>,
    pub edge_signatures: HedgeVec<LoopExtSignature>,
}

impl LoopMomentumBasis {
    pub fn spatial_emr<T: FloatLike>(&self, sample: &BareSample<T>) -> Vec<ThreeMomentum<F<T>>> {
        let three_externals = sample
            .external_moms
            .iter()
            .map(|m| m.spatial.clone())
            .collect_vec();
        self.edge_signatures
            .iter()
            .map(|sig| sig.compute_momentum(&sample.loop_moms, &three_externals))
            .collect()
    }

    pub fn to_massless_emr<T: FloatLike>(&self, sample: &BareSample<T>) -> Vec<FourMomentum<F<T>>> {
        self.edge_signatures
            .iter()
            .map(|sig| {
                sig.compute_four_momentum_from_three(&sample.loop_moms, &sample.external_moms)
            })
            .collect()
    }

    pub fn pattern(&self, edge_id: usize) -> Pattern {
        let signature = self.edge_signatures[edge_id].clone();

        let mut atom = Atom::new_num(0);

        for (i, sign) in signature.internal.into_iter().enumerate() {
            let k = sign * Atom::parse(&format!("K({},x_)", i)).unwrap();

            atom = &atom + &k;
        }

        for (i, sign) in signature.external.into_iter().enumerate() {
            let p = sign * Atom::parse(&format!("P({},x_)", i)).unwrap();
            atom = &atom + &p;
        }

        atom.to_pattern()
    }

    pub fn set_edge_signatures(&mut self, graph: &BareGraph) -> Result<(), Report> {
        // Initialize signature
        self.edge_signatures = vec![
            LoopExtSignature {
                internal: LoopSignature::from_iter(vec![SignOrZero::Zero; self.basis.len()]),
                external: ExternalSignature::from_iter(vec![
                    SignOrZero::Zero;
                    graph.external_edges.len()
                ])
            };
            graph.edges.len()
        ];

        // Build the adjacency list excluding vetoed edges
        let mut adj_list: HashMap<usize, Vec<(usize, usize, bool)>> = HashMap::new();
        for (edge_index, edge) in graph.edges.iter().enumerate() {
            if self.basis.contains(&edge_index) {
                continue;
            }
            let (u, v) = (edge.vertices[0], edge.vertices[1]);

            // Original orientation
            adj_list.entry(u).or_default().push((v, edge_index, false));
            // Flipped orientation
            adj_list.entry(v).or_default().push((u, edge_index, true));
        }

        // Route internal LMB momenta
        for (i_lmb, lmb_edge_id) in self.basis.iter().enumerate() {
            let edge = &graph.edges[*lmb_edge_id];
            let (u, v) = (edge.vertices[0], edge.vertices[1]);

            self.edge_signatures[*lmb_edge_id].internal.0[i_lmb] = SignOrZero::Plus;
            if let Some(path) = self.find_shortest_path(&adj_list, v, u) {
                for (edge_index, is_flipped) in path {
                    if self.edge_signatures[edge_index].internal.0[i_lmb] != SignOrZero::Zero {
                        return Err(eyre!(
                            "Inconsitency in edge momentum lmb signature assignment."
                        ));
                    }
                    self.edge_signatures[edge_index].internal.0[i_lmb] = if is_flipped {
                        SignOrZero::Minus
                    } else {
                        SignOrZero::Plus
                    };
                }
            } else {
                return Err(eyre!(
                    "No path found between vertices {} and {} for LMB: {:?}",
                    u,
                    v,
                    self.basis
                ));
            }
        }

        let sink_node = if let Some(last_external) = graph.external_edges.last() {
            match graph.edges[*last_external].edge_type {
                EdgeType::Outgoing => graph.edges[*last_external].vertices[1],
                EdgeType::Incoming => graph.edges[*last_external].vertices[0],
                _ => {
                    return Err(eyre!(
                        "External edge {} is not incoming or outgoing.",
                        graph.edges[*last_external].name
                    ))
                }
            }
        } else {
            0
        };

        // Route external momenta
        if graph.external_edges.len() >= 2 {
            for i_ext in 0..=(graph.external_edges.len() - 2) {
                let external_edge_index = graph.external_edges[i_ext];
                let external_edge = &graph.edges[external_edge_index];
                let (u, v) = match external_edge.edge_type {
                    EdgeType::Outgoing => (sink_node, external_edge.vertices[1]),
                    EdgeType::Incoming => (external_edge.vertices[0], sink_node),
                    _ => {
                        return Err(eyre!(
                            "External edge {} is not incoming or outgoing.",
                            external_edge.name
                        ))
                    }
                };

                if let Some(path) = self.find_shortest_path(&adj_list, u, v) {
                    //println!("External path from {}->{}: {} {:?}", u, v, i_ext, path);
                    for (edge_index, is_flipped) in path {
                        if self.edge_signatures[edge_index].external.0[i_ext] != SignOrZero::Zero {
                            return Err(eyre!(
                                "Inconsitency in edge momentum signature assignment."
                            ));
                        }
                        self.edge_signatures[edge_index].external.0[i_ext] = if is_flipped {
                            SignOrZero::Minus
                        } else {
                            SignOrZero::Plus
                        };
                    }
                } else {
                    return Err(eyre!(
                        "No path found between vertices {} and {} for LMB: {:?}",
                        u,
                        v,
                        self.basis
                    ));
                }
                if self.edge_signatures[external_edge_index].external.0[i_ext] != SignOrZero::Plus {
                    return Err(eyre!(
                        "Inconsitency in edge momentum external signature assignment."
                    ));
                }
            }
        }
        Ok(())
    }

    fn find_shortest_path(
        &self,
        adjacency_list: &HashMap<usize, Vec<(usize, usize, bool)>>,
        start: usize,
        end: usize,
    ) -> Option<Vec<(usize, bool)>> {
        if start == end {
            return Some(vec![]);
        }

        // Initialize BFS
        let mut queue = VecDeque::new();
        let mut visited: HashMap<usize, Option<(usize, usize, bool)>> = HashMap::new();

        queue.push_back(start);
        visited.insert(start, None);

        // Perform BFS
        while let Some(u) = queue.pop_front() {
            if u == end {
                break;
            }
            if let Some(neighbors) = adjacency_list.get(&u) {
                for &(v, edge_index, is_flipped) in neighbors {
                    #[allow(clippy::map_entry)]
                    if !visited.contains_key(&v) {
                        visited.insert(v, Some((u, edge_index, is_flipped)));
                        queue.push_back(v);
                    }
                }
            }
        }

        // Reconstruct the path if end is reached
        if !visited.contains_key(&end) {
            return None;
        }

        let mut path = Vec::new();
        let mut current = end;

        while let Some(Some((prev, edge_index, is_flipped))) = visited.get(&current) {
            path.push((*edge_index, *is_flipped));
            current = *prev;
        }

        path.reverse();
        Some(path)
    }
}
