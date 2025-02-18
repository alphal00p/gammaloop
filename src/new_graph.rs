use std::{
    borrow::Borrow,
    collections::VecDeque,
    path::{Path, PathBuf},
    sync::Arc,
};

use ahash::{HashMap, HashMapExt, HashSet};
use bincode::{Decode, Encode};
use bitvec::vec::BitVec;
use color_eyre::{Report, Result};
use eyre::eyre;
use hyperdual::Num;
use itertools::Itertools;
use linnet::half_edge::{
    hedgevec::HedgeVec,
    involution::{EdgeIndex, Flow, Hedge, HedgePair},
    subgraph::{
        self, cycle::SignedCycle, Inclusion, InternalSubGraph, OrientedCut, SubGraph, SubGraphOps,
    },
    tree::TraversalTree,
    HedgeGraph, NodeIndex,
};
use petgraph::graph;
use serde::{Deserialize, Serialize};
use smartstring::{LazyCompact, SmartString};
use spenso::structure::{
    abstract_index::AbstractIndex, representation::Minkowski, NamedStructure, ToSymbolic,
};
use symbolica::{
    atom::{representation::InlineNum, Atom, AtomCore, AtomView},
    coefficient::CoefficientView,
    graph::Node,
    id::{Pattern, Replacement},
    parse, with_default_namespace,
};
use typed_index_collections::TiVec;

use crate::{
    gammaloop_integrand::BareSample,
    graph::{DerivedGraphData, EdgeType, Vertex},
    model::{self, EdgeSlots},
    momentum::{FourMomentum, SignOrZero, Signature, ThreeMomentum},
    momentum_sample::{
        BareMomentumSample, ExternalFourMomenta, ExternalIndex, ExternalThreeMomenta, LoopIndex,
    },
    numerator::{ufo::preprocess_ufo_spin_wrapped, NumeratorState, PythonState, UnInit},
    signature::{ExternalSignature, LoopExtSignature, LoopSignature},
    utils::{FloatLike, F, GS},
    ProcessSettings, GAMMALOOP_NAMESPACE,
};

pub struct Graph<S: NumeratorState = PythonState> {
    pub multiplicity: Atom,
    pub underlying: HedgeGraph<Edge, Vertex>,
    pub loop_momentum_basis: LoopMomentumBasis,
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
    fn new_lmb(&self) -> Result<LoopMomentumBasis>;
    fn num_virtual_edges(&self, subgraph: BitVec) -> usize;
    fn is_incoming_to(&self, edge: EdgeIndex, vertex: NodeIndex) -> bool;
    fn denominator(&self, edge: EdgeIndex) -> (Atom, Atom);
    fn substitute_lmb(&self, edge: EdgeIndex, atom: Atom, lmb: &LoopMomentumBasis) -> Atom;
    fn in_slot(&self, edge: EdgeIndex) -> EdgeSlots<Minkowski>;
    fn out_slot(&self, edge: EdgeIndex) -> EdgeSlots<Minkowski>;
}

impl FeynmanGraph for HedgeGraph<Edge, Vertex> {
    fn new_lmb(&self) -> Result<LoopMomentumBasis> {
        // root node should contain a dangling (external edge), that will be the dependent external
        let root = self
            .iter_nodes()
            .find(|(a, _)| {
                self.iter_edges(*a)
                    .any(|(e, _, _)| matches!(e, HedgePair::Unpaired { .. }))
            })
            .unwrap_or(self.iter_nodes().next().unwrap())
            .0;

        let tree = TraversalTree::dfs(self, &self.full_filter(), root, None);
        let tree_complement = tree.tree.complement(self);

        let mut lmb_vec = Vec::new();
        let mut cycle_basis = Vec::new();
        let mut lmb_basis = Vec::new();

        for (pair, edge_index, _) in self.iter_edges(&tree_complement) {
            if let HedgePair::Paired { source, .. } = pair {
                lmb_vec.push(edge_index);
                lmb_basis.push(source);
                cycle_basis.push(
                    SignedCycle::from_cycle(tree.cycle(source).unwrap(), source, self).unwrap(),
                );
            }
        }

        let loop_number = lmb_vec.len();
        let mut lmb = LoopMomentumBasis {
            tree,
            basis: TiVec::from_iter(lmb_vec),
            edge_signatures: self.new_hedgevec(&|_, _| LoopExtSignature {
                internal: LoopSignature::from_iter(vec![SignOrZero::Zero; loop_number]),
                external: ExternalSignature::from_iter(vec![SignOrZero::Zero; self.n_externals()]),
            }),
        };

        lmb.set_edge_signatures(self)?;

        Ok(lmb)

        //        let mut leaves = tree.leaf_edges();
        //
        //        let mut external_leaves = self.empty_filter();
        //        let mut externals = vec![];
        //
        //        for i in leaves.included_iter() {
        //            if self.involution.is_identity(i) {
        //                external_leaves.set(i.0, true);
        //                externals.push(i);
        //            }
        //        }
        //
        //        let mut ext_signatures: HedgeVec<Signature> = self.new_derived_edge_data_empty();
        //
        //        let empty_signature = Signature::from_iter(externals.iter().map(|e| SignOrZero::Zero));
        //        for (i, &h) in externals.iter().enumerate() {
        //            let mut signature = empty_signature.clone();
        //            signature.0[i] = SignOrZero::Plus;
        //            ext_signatures[h] = Some(signature);
        //        }
        //
        //        let mut current_leaf_nodes = tree.leaf_nodes(&self);
        //
        //        while let Some(leaf_node) = current_leaf_nodes.pop() {
        //            let hairs = &self.hairs_from_id(leaf_node).hairs;
        //            let mut root_pointer = None;
        //            let mut root_signature = empty_signature.clone();
        //
        //            for h in hairs.included_iter() {
        //                match tree.parent(h) {
        //                    Parent::Root => {}
        //                    Parent::Hedge { hedge_to_root, .. } => {
        //                        if *hedge_to_root == h {
        //                            root_pointer = Some(h);
        //                            if self
        //                                .involved_node_hairs(h)
        //                                .unwrap()
        //                                .hairs
        //                                .included_iter()
        //                                .all(|a| a != h && ext_signatures.is_set(a))
        //                            {
        //                                current_leaf_nodes.push(self.involved_node_id(h).unwrap());
        //                            }
        //                        } else {
        //                            root_signature.sum(ext_signatures[h].as_ref().unwrap());
        //                        }
        //                    }
        //                    Parent::Unset => {}
        //                }
        //            }
        //
        //            if let Some(root_pointer) = root_pointer {
        //                ext_signatures[root_pointer] = Some(root_signature);
        //            }
        //        }
        //
        //        let signatures = self
        //            .involution
        //            .map_data_ref(|a| (), &|e| ())
        //            .map_edge_data(|e, d| {
        //                let e = match e {
        //                    EdgeId::Paired { source, sink } => {
        //                        let mut internal = vec![];
        //                        for (i, c) in cycle_basis.iter().enumerate() {
        //                            if c.filter.includes(&source) {
        //                                internal.push(SignOrZero::Plus);
        //                            } else if c.filter.includes(&sink) {
        //                                internal.push(SignOrZero::Minus);
        //                            } else {
        //                                internal.push(SignOrZero::Zero);
        //                            }
        //                        }
        //
        //                        let internal_signature = Signature::from_iter(internal);
        //
        //                        // return EdgeData::new(Signature::from_iter(iter), orientation)
        //                    }
        //                    EdgeId::Unpaired { hedge, flow } => {}
        //                    _ => {}
        //                };
        //                d
        //            });
        //
        //        HedgeLMB { tree, lmb_basis }
    }

    fn num_virtual_edges(&self, subgraph: BitVec) -> usize {
        let internal_subgraph = InternalSubGraph::cleaned_filter_pessimist(subgraph, self);
        self.count_internal_edges(&internal_subgraph)
    }

    fn is_incoming_to(&self, edge: EdgeIndex, vertex: NodeIndex) -> bool {
        let (_, pair) = self[&edge];
        match pair {
            HedgePair::Unpaired { hedge, flow } => {
                self.node_id(hedge) == vertex && matches!(flow, Flow::Sink)
            }
            HedgePair::Paired { source: _, sink } => self.node_id(sink) == vertex,
            HedgePair::Split {
                source: _,
                sink,
                split: _,
            } => self.node_id(sink) == vertex,
        }
    }

    fn denominator(&self, edge: EdgeIndex) -> (Atom, Atom) {
        let mom = parse!(&format!("Q{}", Into::<usize>::into(edge))).unwrap();
        let mass = self[edge]
            .particle
            .mass
            .expression
            .clone()
            .unwrap_or(Atom::new_num(0));

        (mom, mass)
    }

    fn substitute_lmb(&self, edge: EdgeIndex, atom: Atom, lmb: &LoopMomentumBasis) -> Atom {
        let num = Into::<usize>::into(edge);
        let mom = parse!(&format!("Q({num},x_)")).unwrap().to_pattern();
        let mom_rep = lmb.pattern(edge);
        atom.replace_all(&mom, &mom_rep, None, None)
    }

    fn in_slot(&self, edge: EdgeIndex) -> EdgeSlots<Minkowski> {
        let source = match self[&edge].1 {
            HedgePair::Unpaired { hedge, flow } => todo!(),
            HedgePair::Paired { source, sink } => self.node_id(source),
            HedgePair::Split {
                source,
                sink,
                split,
            } => self.node_id(source),
        };
        //let local_pos_in_sink_vertex =
        //    self[source].get_local_edge_position(&self[edge], self, false);
        todo!()
    }

    fn out_slot(&self, edge: EdgeIndex) -> EdgeSlots<Minkowski> {
        todo!()
    }
}

impl Graph<UnInit> {
    pub fn new(multiplicity: Atom, underlying: HedgeGraph<Edge, Vertex>) -> Result<Self> {
        Ok(Self {
            multiplicity,
            loop_momentum_basis: underlying.new_lmb()?,
            underlying,
            derived_data: DerivedGraphData::new_empty(),
        })
    }
}

#[derive(Debug, Clone)]
pub struct Edge {
    pub name: SmartString<LazyCompact>,
    pub edge_type: EdgeType,
    pub propagator: Arc<model::Propagator>,
    pub particle: Arc<model::Particle>,
    pub vertices: [Option<usize>; 2], // do we keep this redundant info?
    pub internal_index: Vec<AbstractIndex>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct SerializableEdge {
    name: SmartString<LazyCompact>,
    edge_type: EdgeType,
    particle: SmartString<LazyCompact>,
    propagator: SmartString<LazyCompact>,
    vertices: [SmartString<LazyCompact>; 2],
}

impl SerializableEdge {
    pub fn from_edge(graph: &HedgeGraph<Edge, Vertex>, edge: &Edge) -> SerializableEdge {
        //SerializableEdge {
        //    name: edge.name.clone(),
        //    edge_type: edge.edge_type,
        //    particle: edge.particle.name.clone(),
        //    propagator: edge.propagator.name.clone(),
        //    vertices: [
        //        graph.vertices[edge.vertices[0]].name.clone(),
        //        graph.vertices[edge.vertices[1]].name.clone(),
        //    ],
        //}

        todo!()
    }
}

impl Edge {
    pub fn n_dummies(&self) -> usize {
        5
    }
    pub fn from_serializable_edge(
        model: &model::Model,
        graph: &HedgeGraph<Edge, Vertex>,
        serializable_edge: &SerializableEdge,
    ) -> Edge {
        todo!()
        //Edge {
        //    name: serializable_edge.name.clone(),
        //    edge_type: serializable_edge.edge_type,
        //    particle: model.get_particle(&serializable_edge.particle),
        //    propagator: model.get_propagator(&serializable_edge.propagator),
        //    vertices: [
        //        graph
        //            .get_vertex_position(&serializable_edge.vertices[0])
        //            .unwrap(),
        //        graph
        //            .get_vertex_position(&serializable_edge.vertices[1])
        //            .unwrap(),
        //    ],
        //    internal_index: vec![],
        //}
    }

    // pub fn is_incoming_to(&self, vertex: usize) -> bool {
    //     self.vertices[1] == vertex
    // }

    //pub fn denominator(&self, graph: &BareGraph) -> (Atom, Atom) {
    //    let num = *graph.edge_name_to_position.get(&self.name).unwrap();
    //    let mom = parse!(&format!("Q{num}")).unwrap();
    //    let mass = self
    //        .particle
    //        .mass
    //        .expression
    //        .clone()
    //        .unwrap_or(Atom::new_num(0));

    //    (mom, mass)
    //}

    //pub fn substitute_lmb(&self, atom: Atom, graph: &BareGraph, lmb: &LoopMomentumBasis) -> Atom {
    //     let num = *graph.edge_name_to_position.get(&self.name).unwrap();
    //     let mom = parse!(&format!("Q({num},x_)")).unwrap().to_pattern();
    //     let mom_rep = lmb.pattern(num);
    //     atom.replace_all(&mom, &mom_rep, None, None)
    // }

    // pub fn edge_momentum_symbol(&self, graph: &Graph) -> Symbol {
    //     let num = *graph.edge_name_to_position.get(&self.name).unwrap();
    //     State::get_symbol(format!("Q{num}"))
    // }

    //pub fn in_slot(&self, graph: &BareGraph) -> EdgeSlots<Minkowski> {
    //    let local_pos_in_sink_vertex =
    //        graph.vertices[self.vertices[0]].get_local_edge_position(self, graph, false);

    //    graph.vertex_slots[self.vertices[0]][local_pos_in_sink_vertex].dual()
    //}

    //pub fn out_slot(&self, graph: &BareGraph) -> EdgeSlots<Minkowski> {
    //    let local_pos_in_sink_vertex = graph.vertices[self.vertices[1]].get_local_edge_position(
    //        self,
    //        graph,
    //        self.vertices[0] == self.vertices[1],
    //    );

    //    graph.vertex_slots[self.vertices[1]][local_pos_in_sink_vertex].dual()
    //}

    pub fn numerator(&self, graph: &HedgeGraph<Edge, Vertex>, edge_index: EdgeIndex) -> Atom {
        let [colorless, color] = self.color_separated_numerator(graph, edge_index);

        colorless * color
    }

    pub fn color_separated_numerator(
        &self,
        graph: &HedgeGraph<Edge, Vertex>,
        num: EdgeIndex,
    ) -> [Atom; 2] {
        // let num = *graph.edge_name_to_position.get(&self.name).unwrap();
        let in_slots = graph.in_slot(num);
        let out_slots = graph.out_slot(num);

        match self.edge_type {
            EdgeType::Incoming => {
                let [lorentz, spin, color] = in_slots.dual().kroneker(&out_slots);

                [lorentz * spin, color]
            }
            EdgeType::Outgoing => {
                let [lorentz, spin, color] = out_slots.dual().kroneker(&in_slots);

                [lorentz * spin, color]
            }
            EdgeType::Virtual => {
                let mut atom = self.propagator.numerator.clone();

                let pfun = parse!("P(x_)").unwrap().to_pattern();
                if self.particle.is_antiparticle() {
                    atom = atom.replace_all(
                        &pfun,
                        parse!(&format!(
                            "-Q({},mink(4,indexid(x_)))",
                            Into::<usize>::into(num)
                        ))
                        .unwrap()
                        .to_pattern(),
                        None,
                        None,
                    );
                } else {
                    atom = atom.replace_all(
                        &pfun,
                        parse!(&format!(
                            "Q({},mink(4,indexid(x_)))",
                            Into::<usize>::into(num)
                        ))
                        .unwrap()
                        .to_pattern(),
                        None,
                        None,
                    );
                }

                let pslashfun = parse!("PSlash(i_,j_)").unwrap().to_pattern();
                let pindex_num: usize = self.internal_index[0].into();
                if self.particle.is_antiparticle() {
                    atom = atom.replace_all(
                        &pslashfun,
                        parse!(&format!(
                            "-Q({},mink(4,{}))*Gamma({},i_,j_)",
                            Into::<usize>::into(num),
                            pindex_num,
                            pindex_num
                        ))
                        .unwrap()
                        .to_pattern(),
                        None,
                        None,
                    );
                } else {
                    atom = atom.replace_all(
                        &pslashfun,
                        parse!(&format!(
                            "Q({},mink(4,{}))*Gamma({},i_,j_)",
                            Into::<usize>::into(num),
                            pindex_num,
                            pindex_num
                        ))
                        .unwrap()
                        .to_pattern(),
                        None,
                        None,
                    );
                }

                atom = preprocess_ufo_spin_wrapped(atom);
                let indexidpat = parse!("indexid(x_)").unwrap().to_pattern();

                let dummies: HashSet<_> = atom
                    .pattern_match(&indexidpat, None, None)
                    .filter_map(|a| {
                        if let AtomView::Num(n) = a[&GS.x_].as_view() {
                            let e = if let CoefficientView::Natural(a, b) = n.get_coeff_view() {
                                if b == 1 {
                                    a
                                } else {
                                    0
                                }
                            } else {
                                0
                            };
                            if e < 0 {
                                Some(e)
                            } else {
                                None
                            }
                        } else {
                            None
                        }
                    })
                    .collect();

                let (replacements_in, replacements_out) = if self.particle.is_antiparticle() {
                    (in_slots.replacements(2), out_slots.replacements(1))
                } else {
                    (in_slots.replacements(1), out_slots.replacements(2))
                };

                // replacements_out.push(Replacement::new(
                //     parse!("indexid(x_)").unwrap().to_pattern(),
                //     parse!("x_").unwrap().to_pattern(),
                // ));

                let mut color_atom = Atom::new_num(1);
                for (&cin, &cout) in in_slots.color.iter().zip(out_slots.color.iter()) {
                    let id: NamedStructure<String, ()> =
                        NamedStructure::from_iter([cin, cout], "id".into(), None);
                    color_atom = color_atom * &id.to_symbolic().unwrap();
                }

                let reps: Vec<Replacement> = replacements_in
                    .into_iter()
                    .chain(replacements_out)
                    .collect();

                let atom = atom.replace_all_multiple(&reps);
                let color_atom = color_atom.replace_all_multiple(&reps);

                let indexid_reps: Vec<_> = dummies
                    .into_iter()
                    .enumerate()
                    .sorted()
                    .map(|(i, d)| {
                        Replacement::new(
                            parse!(&format!("indexid({})", d)).unwrap().to_pattern(),
                            parse!(&format!("{}", self.internal_index[i + 1]))
                                .unwrap()
                                .to_pattern(),
                        )
                    })
                    .collect();

                let atom = atom.replace_all_multiple(&indexid_reps);
                let color_atom = color_atom.replace_all_multiple(&indexid_reps);

                [
                    atom.replace_all(
                        &parse!("indexid(x_)").unwrap().to_pattern(),
                        Atom::new_var(GS.x_).to_pattern(),
                        None,
                        None,
                    ),
                    color_atom.replace_all(
                        &parse!("indexid(x_)").unwrap().to_pattern(),
                        Atom::new_var(GS.x_).to_pattern(),
                        None,
                        None,
                    ),
                ]
            }
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize /*, Encode, Decode*/)]

pub struct LoopMomentumBasis {
    pub tree: TraversalTree,
    pub basis: TiVec<LoopIndex, EdgeIndex>,
    pub edge_signatures: HedgeVec<LoopExtSignature>,
}

impl LoopMomentumBasis {
    pub fn spatial_emr<T: FloatLike>(
        &self,
        sample: &BareMomentumSample<T>,
    ) -> Vec<ThreeMomentum<F<T>>> {
        let three_externals: ExternalThreeMomenta<F<T>> = sample
            .external_moms
            .iter()
            .map(|m| m.spatial.clone())
            .collect();
        self.edge_signatures
            .borrow()
            .into_iter()
            .map(|(_, sig)| sig.compute_momentum(&sample.loop_moms, &three_externals))
            .collect()
    }

    pub fn to_massless_emr<T: FloatLike>(
        &self,
        sample: &BareMomentumSample<T>,
    ) -> Vec<FourMomentum<F<T>>> {
        self.edge_signatures
            .borrow()
            .into_iter()
            .map(|(_, sig)| {
                sig.compute_four_momentum_from_three(&sample.loop_moms, &sample.external_moms)
            })
            .collect()
    }

    pub fn pattern(&self, edge_id: EdgeIndex) -> Pattern {
        let signature = self.edge_signatures[edge_id].clone();

        let mut atom = Atom::new_num(0);

        for (i, sign) in signature.internal.into_iter().enumerate() {
            let k = sign
                * Atom::parse(with_default_namespace!(
                    &format!("K({},x_)", i),
                    GAMMALOOP_NAMESPACE
                ))
                .unwrap();

            atom = &atom + &k;
        }

        for (i, sign) in signature.external.into_iter().enumerate() {
            let p = sign
                * Atom::parse(with_default_namespace!(
                    &format!("P({},x_)", i),
                    GAMMALOOP_NAMESPACE
                ))
                .unwrap();
            atom = &atom + &p;
        }

        atom.to_pattern()
    }

    pub fn set_edge_signatures<E, V>(&mut self, graph: &HedgeGraph<E, V>) -> Result<(), Report> {
        self.edge_signatures = graph.new_hedgevec(&|_, _edge_index| LoopExtSignature {
            internal: LoopSignature::from_iter(vec![SignOrZero::Zero; self.basis.len()]),
            external: ExternalSignature::from_iter(vec![SignOrZero::Zero; graph.n_externals()]),
        });

        struct ExternalEdgeInfo {
            edge_index: EdgeIndex,
            fake_node: NodeIndex, // fake node is a hack to reuse the code from BareGraph
            real_node: NodeIndex,
            flow: Flow,
        }

        let mut current_extra_node = graph.n_nodes();
        let mut external_edge_info = TiVec::<ExternalIndex, ExternalEdgeInfo>::new();

        // Build the adjacency list excluding vetoed edges, we include "fake nodes" on the externals such that we do
        // not need to port too much of the code.
        let mut adj_list: HashMap<NodeIndex, Vec<(NodeIndex, EdgeIndex, bool)>> = HashMap::new();
        for (hedge_pair, edge_index, _edge_data) in graph.iter_all_edges() {
            if self.basis.contains(&edge_index) {
                continue;
            }
            // let (u, v) = (edge.vertices[0], edge.vertices[1]);

            match hedge_pair {
                HedgePair::Unpaired { hedge, flow } => {
                    let real_node = graph.node_id(hedge);
                    let extra_node = NodeIndex::from(current_extra_node);
                    external_edge_info.push(ExternalEdgeInfo {
                        edge_index,
                        fake_node: extra_node,
                        real_node,
                        flow,
                    });
                    match flow {
                        Flow::Sink => {
                            adj_list
                                .entry(extra_node)
                                .or_default()
                                .push((real_node, edge_index, false));
                            adj_list
                                .entry(real_node)
                                .or_default()
                                .push((extra_node, edge_index, true));
                        }
                        Flow::Source => {
                            adj_list
                                .entry(real_node)
                                .or_default()
                                .push((extra_node, edge_index, false));
                            adj_list
                                .entry(extra_node)
                                .or_default()
                                .push((real_node, edge_index, true));
                        }
                    }
                    current_extra_node += 1;
                }
                HedgePair::Paired { source, sink } => {
                    let u = graph.node_id(source);
                    let v = graph.node_id(sink);

                    // Original orientation
                    adj_list.entry(u).or_default().push((v, edge_index, false));
                    // Flipped orientation
                    adj_list.entry(v).or_default().push((u, edge_index, true));
                }
                HedgePair::Split { .. } => {
                    panic!("Can not set edge signatures for split edges yet")
                }
            }
        }

        // Route internal LMB momenta
        for (i_lmb, lmb_edge_id) in self.basis.iter_enumerated() {
            let (_, hedge_pair) = &graph[lmb_edge_id];

            let (u, v) = match hedge_pair {
                HedgePair::Paired { source, sink } => {
                    (graph.node_id(*source), graph.node_id(*sink))
                }
                _ => {
                    return Err(eyre!(
                        "Loop_momentum {} is an external edge. Edge: {:?}",
                        i_lmb,
                        lmb_edge_id
                    ))
                }
            };

            self.edge_signatures[*lmb_edge_id].internal[i_lmb] = SignOrZero::Plus;
            if let Some(path) = self.find_shortest_path(&adj_list, v, u) {
                for (edge_index, is_flipped) in path {
                    if self.edge_signatures[edge_index].internal[i_lmb] != SignOrZero::Zero {
                        return Err(eyre!(
                            "Inconsitency in edge momentum lmb signature assignment."
                        ));
                    }
                    self.edge_signatures[edge_index].internal[i_lmb] = if is_flipped {
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

        // sink node is the last external node
        let sink_node = external_edge_info
            .pop()
            .map(|edge_info| edge_info.fake_node)
            .unwrap_or(NodeIndex::from(0));

        // Route external momenta
        if graph.n_externals() >= 2 {
            for (external_index, external_edge_info) in external_edge_info.into_iter_enumerated() {
                let (u, v) = match external_edge_info.flow {
                    Flow::Source => (sink_node, external_edge_info.fake_node),
                    Flow::Sink => (external_edge_info.fake_node, sink_node),
                };

                if let Some(path) = self.find_shortest_path(&adj_list, u, v) {
                    //println!("External path from {}->{}: {} {:?}", u, v, i_ext, path);
                    for (edge_index, is_flipped) in path {
                        if self.edge_signatures[edge_index].external[external_index]
                            != SignOrZero::Zero
                        {
                            return Err(eyre!(
                                "Inconsitency in edge momentum signature assignment."
                            ));
                        }
                        self.edge_signatures[edge_index].external[external_index] = if is_flipped {
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
                if self.edge_signatures[external_edge_info.edge_index].external[external_index]
                    != SignOrZero::Plus
                {
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
        adjacency_list: &HashMap<NodeIndex, Vec<(NodeIndex, EdgeIndex, bool)>>,
        start: NodeIndex,
        end: NodeIndex,
    ) -> Option<Vec<(EdgeIndex, bool)>> {
        if start == end {
            return Some(vec![]);
        }

        // Initialize BFS
        let mut queue = VecDeque::new();
        let mut visited: HashMap<NodeIndex, Option<(NodeIndex, EdgeIndex, bool)>> = HashMap::new();

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
