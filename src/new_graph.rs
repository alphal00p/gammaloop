use std::{
    borrow::Borrow,
    collections::VecDeque,
    path::{Path, PathBuf},
    sync::Arc,
};

use ahash::{AHashMap, AHashSet, HashMap, HashMapExt, HashSet};
use bincode::{Decode, Encode};
use bitvec::vec::BitVec;
use color_eyre::{Report, Result};
use derive_more::{From, Into};
use eyre::eyre;
use hyperdual::Num;
use itertools::Itertools;
use linnet::half_edge::{
    builder::HedgeGraphBuilder,
    hedgevec::HedgeVec,
    involution::{EdgeData, EdgeIndex, Flow, Hedge, HedgePair, Orientation},
    subgraph::{
        self, cycle::SignedCycle, node, Inclusion, InternalSubGraph, OrientedCut, SubGraph,
        SubGraphOps,
    },
    HedgeGraph, NodeIndex,
};
use log::debug;
use nalgebra::DMatrix;
// use petgraph::Direction::Outgoing;
use serde::{de::value, Deserialize, Serialize};
use smartstring::{LazyCompact, SmartString};
use spenso::{
    contraction::IsZero,
    data::DataTensor,
    scalar::Scalar,
    structure::{
        abstract_index::AbstractIndex, representation::Minkowski, NamedStructure, ToSymbolic,
    },
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
    cff::generation::ShiftRewrite,
    disable,
    feyngen::diagram_generator::FeynGen,
    gammaloop_integrand::BareSample,
    graph::{
        BareEdge, BareGraph, BareVertex, DerivedGraphData, EdgeType, HasVertexInfo, Shifts,
        VertexInfo,
    },
    model::{self, ArcParticle, ArcPropagator, EdgeSlots, Model, Particle, VertexSlots},
    momentum::{FourMomentum, SignOrZero, Signature, ThreeMomentum},
    momentum_sample::{
        BareMomentumSample, ExternalFourMomenta, ExternalIndex, ExternalThreeMomenta, LoopIndex,
        LoopMomenta,
    },
    new_gammaloop_integrand::LmbMultiChannelingSetup,
    numerator::{ufo::preprocess_ufo_spin_wrapped, NumeratorState, PythonState, UnInit},
    signature::{ExternalSignature, LoopExtSignature, LoopSignature, SignatureLike},
    utils::{external_energy_atom_from_index, ose_atom_from_index, FloatLike, F, GS},
    GammaLoopContext, ProcessSettings, GAMMALOOP_NAMESPACE,
};

#[derive(Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = crate::GammaLoopContext)]
pub struct Graph {
    pub multiplicity: Atom,
    pub name: String,
    pub underlying: HedgeGraph<Edge, Vertex>,
    pub loop_momentum_basis: LoopMomentumBasis,
    pub vertex_slots: TiVec<NodeIndex, VertexSlots>,
    pub external_connections: Option<Vec<ExternalConnection>>,
}

impl From<BareGraph> for Graph {
    fn from(value: BareGraph) -> Self {
        let loop_momentum_basis = value.loop_momentum_basis.clone();
        let vertex_slots = value.vertex_slots.clone().into();
        let multiplicity = FeynGen::evaluate_overall_factor(value.overall_factor.as_view());
        let name = value.name.clone();

        // convert old external connections to new format
        let external_connections = value
            .external_connections
            .iter()
            .map(|(incoming_node, outgoing_node)| {
                let incoming_index = incoming_node.map(|node_id| {
                    let edge_id = value.vertices[node_id].edges[0];
                    let external_index = value.get_external_index(edge_id).unwrap();
                    external_index
                });

                let outgoing_index = outgoing_node.map(|node_id| {
                    let edge_id = value.vertices[node_id].edges[0];
                    let external_index = value.get_external_index(edge_id).unwrap();
                    external_index
                });

                if let (Some(incoming_index), Some(outgoing_index)) =
                    (incoming_index, outgoing_index)
                {
                    Some(ExternalConnection {
                        incoming_index,
                        outgoing_index,
                    })
                } else {
                    None
                }
            })
            .collect::<Option<Vec<ExternalConnection>>>();

        let underlying = value.into();
        Self {
            name: name.to_string(),
            external_connections,
            multiplicity,
            vertex_slots,
            loop_momentum_basis,
            underlying,
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
    fn get_local_edge_position(
        &self,
        node_id: NodeIndex,
        edge_id: EdgeIndex,
        skip_one: bool,
    ) -> usize;
    fn add_signs_to_edges(&self, node_id: NodeIndex) -> Vec<isize>;
    fn get_cff_inverse_energy_product(&self) -> Atom;
    fn get_loop_number(&self) -> usize;
    fn get_real_mass_vector<T: FloatLike>(&self) -> HedgeVec<F<T>>;
    fn get_energy_cache<T: FloatLike>(
        &self,
        loop_moms: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        lmb: &LoopMomentumBasis,
    ) -> HedgeVec<F<T>>;
    fn get_esurface_canonization(&self, lmb: &LoopMomentumBasis) -> Option<ShiftRewrite>;
    fn external_in_or_out_signature(&self) -> ExternalSignature;
    fn get_external_partcles(&self) -> Vec<ArcParticle>;
    fn get_external_signature(&self) -> SignatureLike<ExternalIndex>;
    fn get_energy_atoms(&self) -> Vec<Atom>;
}

impl FeynmanGraph for HedgeGraph<Edge, Vertex> {
    fn new_lmb(&self) -> Result<LoopMomentumBasis> {
        todo!();
        // root node should contain a dangling (external edge), that will be the dependent external
        disable! {
        let root = self
            .iter_nodes()
            .find(|(a, _, _)| {
                self.iter_edges(a)
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
                tree: Some(tree),
                basis: TiVec::from_iter(lmb_vec),
                edge_signatures: self.new_hedgevec(&|_, _| LoopExtSignature {
                    internal: LoopSignature::from_iter(vec![SignOrZero::Zero; loop_number]),
                    external: ExternalSignature::from_iter(vec![SignOrZero::Zero; self.n_externals()]),
                }),
            };

            lmb.set_edge_signatures(self)?;

            Ok(lmb)
        }
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
            .0
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
        atom.replace(mom).with(mom_rep)
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

    fn get_local_edge_position(
        &self,
        node_id: NodeIndex,
        edge_id: EdgeIndex,
        skip_one: bool,
    ) -> usize {
        unimplemented!()
    }

    fn add_signs_to_edges(&self, node_id: NodeIndex) -> Vec<isize> {
        let node_hairs: BitVec = self.iter_crown(node_id).into();

        self.iter_edges(&node_hairs)
            .map(|(_, edge_index, _)| {
                if !self.is_incoming_to(edge_index, node_id) {
                    -(Into::<usize>::into(edge_index) as isize)
                } else {
                    Into::<usize>::into(edge_index) as isize
                }
            })
            .collect()
    }

    /// This includes the factor 2 for each edge, inversion already performed
    fn get_cff_inverse_energy_product(&self) -> Atom {
        let full_subgraph = self.full_filter();
        get_cff_inverse_energy_product_impl(self, &full_subgraph)
    }

    fn get_loop_number(&self) -> usize {
        let internal_subgraph =
            InternalSubGraph::cleaned_filter_pessimist(self.full_filter(), self);
        self.cyclotomatic_number(&internal_subgraph)
    }

    fn get_real_mass_vector<T: FloatLike>(&self) -> HedgeVec<F<T>> {
        self.new_hedgevec(|edge, _edge_id, _| match edge.particle.0.mass.value {
            Some(mass) => F::from_ff64(mass.re),
            None => F::from_f64(0.0),
        })
    }

    fn get_energy_cache<T: FloatLike>(
        &self,
        loop_moms: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        lmb: &LoopMomentumBasis,
    ) -> HedgeVec<F<T>> {
        self.new_hedgevec_from_iter(
            lmb.edge_signatures
                .borrow()
                .into_iter()
                .map(|(_, sig)| sig.compute_four_momentum_from_three(loop_moms, external_moms))
                .zip(self.iter_all_edges())
                .map(|(emr_mom, (_, _, edge))| match edge.data.edge_type {
                    EdgeType::Virtual => {
                        emr_mom
                            .spatial
                            .on_shell_energy(edge.data.particle.0.mass.value.map(|m| {
                                if m.im.is_non_zero() {
                                    panic!("Complex masses not yet supported in gammaLoop")
                                }
                                F::<T>::from_ff64(m.re)
                            }))
                            .value
                    }
                    _ => emr_mom.temporal.value, // a wierd way of just obtaining the energy of the external particles
                }),
        )
        .unwrap()
    }

    fn get_esurface_canonization(&self, lmb: &LoopMomentumBasis) -> Option<ShiftRewrite> {
        let external_edges: TiVec<ExternalIndex, _> = self
            .iter_all_edges()
            .filter(|(pair, _, _)| matches!(pair, HedgePair::Unpaired { .. }))
            .collect();

        // find the external leg which does not appear in it's own signature
        external_edges
            .iter_enumerated()
            .find(|(external_index, (_, edge_id, _))| {
                lmb.edge_signatures[*edge_id].external[*external_index] == SignOrZero::Zero
            })
            .map(|(_, (_, dep_mom_edge_id, _))| {
                let dep_mom_signatrue = &lmb.edge_signatures[*dep_mom_edge_id].external;

                let external_shift = external_edges
                    .iter()
                    .zip(dep_mom_signatrue)
                    .filter(|(_, dep_mom_sign)| dep_mom_sign.is_sign())
                    .map(|((_, external_edge, _), dep_mom_sign)| {
                        (*external_edge, dep_mom_sign as i64)
                    })
                    .collect_vec();

                ShiftRewrite {
                    dependent_momentum: *dep_mom_edge_id,
                    dependent_momentum_expr: external_shift,
                }
            })
    }

    fn external_in_or_out_signature(&self) -> ExternalSignature {
        self.iter_all_edges()
            .filter_map(|(pair, _, _)| match pair {
                HedgePair::Unpaired { flow, .. } => match flow {
                    Flow::Sink => Some(1i8),
                    Flow::Source => Some(-1i8),
                },
                _ => None,
            })
            .collect()
    }

    fn get_external_partcles(&self) -> Vec<ArcParticle> {
        self.iter_all_edges()
            .filter_map(|(pair, _, data)| match pair {
                HedgePair::Unpaired { .. } => Some(data.data.particle.clone()),
                _ => None,
            })
            .collect()
    }

    fn get_external_signature(&self) -> SignatureLike<ExternalIndex> {
        SignatureLike::from_iter(self.iter_all_edges().filter_map(|(pair, _, _)| match pair {
            HedgePair::Unpaired { flow, .. } => match flow {
                Flow::Source => Some(SignOrZero::Minus),
                Flow::Sink => Some(SignOrZero::Plus),
            },
            _ => None,
        }))
    }

    fn get_energy_atoms(&self) -> Vec<Atom> {
        self.iter_all_edges()
            .map(|(pair, edge_id, _)| match pair {
                HedgePair::Paired { .. } => ose_atom_from_index(edge_id),
                HedgePair::Unpaired { .. } => external_energy_atom_from_index(edge_id),
                _ => unreachable!(),
            })
            .collect_vec()
    }
}

impl Graph {
    pub fn new(
        name: SmartString<LazyCompact>,
        multiplicity: Atom,
        underlying: HedgeGraph<Edge, Vertex>,
    ) -> Result<Self> {
        Ok(Self {
            name: name.to_string(),
            multiplicity,
            loop_momentum_basis: underlying.new_lmb()?,
            underlying,
            external_connections: None,
            vertex_slots: TiVec::new(),
        })
    }

    pub fn build_multi_channeling_channels(
        &self,
        lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
    ) -> LmbMultiChannelingSetup {
        let mut channels: Vec<LmbIndex> = Vec::new();

        // Filter out channels that are non-singular, or have the same singularities as another channel already included
        for (lmb_index, lmb) in lmbs.iter_enumerated() {
            let massless_edges = lmb
                .basis
                .iter()
                .filter(|&edge_id| self.underlying[*edge_id].particle.0.is_massless())
                .collect_vec();

            if massless_edges.is_empty() {
                continue;
            }

            if channels.iter().any(|channel| {
                let basis_of_included_channel = &lmbs[*channel].basis;
                massless_edges
                    .iter()
                    .all(|edge_id| basis_of_included_channel.contains(edge_id))
            }) {
                continue;
            }

            // only for 1 to n for now, assuming center of mass
            if self
                .external_connections
                .as_ref()
                .map(|external_connections| external_connections.len() == 1)
                .unwrap_or(false)
                && channels.iter().any(|channel| {
                    let massless_edges_of_included_channel = lmbs[*channel]
                        .basis
                        .iter()
                        .filter(|&edge_id| self.underlying[*edge_id].particle.0.is_massless())
                        .collect_vec();

                    let loop_signatures_of_massless_edges_of_included_channel =
                        massless_edges_of_included_channel
                            .iter()
                            .map(|edge_index| {
                                self.loop_momentum_basis.edge_signatures[**edge_index]
                                    .internal
                                    .first_abs()
                            })
                            .collect::<HashSet<_>>();

                    let loop_signatures_of_massless_edges_of_potential_channel = massless_edges
                        .iter()
                        .map(|edge_index| {
                            self.loop_momentum_basis.edge_signatures[**edge_index]
                                .internal
                                .first_abs()
                        })
                        .collect::<HashSet<_>>();

                    loop_signatures_of_massless_edges_of_included_channel
                        == loop_signatures_of_massless_edges_of_potential_channel
                })
            {
                continue;
            }

            channels.push(lmb_index);
        }

        let channels: TiVec<_, _> = channels.into_iter().sorted().collect();

        debug!(
            "number of lmbs: {}, number of channels: {}",
            lmbs.len(),
            channels.len()
        );

        LmbMultiChannelingSetup { channels }
    }

    pub fn iter_loop_edges(&self) -> impl Iterator<Item = (HedgePair, EdgeIndex, EdgeData<&Edge>)> {
        self.underlying
            .iter_all_edges()
            .filter(|(_, edge_index, _)| {
                self.loop_momentum_basis.edge_signatures[*edge_index]
                    .internal
                    .iter()
                    .any(|sign| sign.is_sign())
            })
    }

    pub fn iter_non_loop_edges(
        &self,
    ) -> impl Iterator<Item = (HedgePair, EdgeIndex, EdgeData<&Edge>)> {
        self.underlying
            .iter_all_edges()
            .filter(|(_, edge_index, _)| {
                self.loop_momentum_basis.edge_signatures[*edge_index]
                    .internal
                    .iter()
                    .all(|sign| sign.is_zero())
            })
    }
}

impl Graph {
    pub fn apply_vertex_rule(&self, node_id: NodeIndex) -> Option<[DataTensor<Atom>; 3]> {
        self.underlying[node_id].vertex_info.apply_vertex_rule(
            &self.underlying.add_signs_to_edges(node_id),
            Into::<usize>::into(node_id),
            &self.vertex_slots[node_id],
        )
    }
}

#[derive(Debug, Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = crate::GammaLoopContext)]
pub struct Edge {
    // #[bincode(with_serde)]
    pub name: String,
    pub edge_type: EdgeType,
    pub propagator: ArcPropagator,
    pub particle: ArcParticle,
    pub num: Atom,
    pub dod: i32,
    // #[bincode(with_serde)]
    pub internal_index: Vec<AbstractIndex>,
}

impl Edge {
    pub fn n_dummies(&self) -> usize {
        5
    }

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
                if self.particle.0.is_antiparticle() {
                    atom = atom.replace(&pfun).with(
                        parse!(&format!(
                            "-Q({},mink(4,indexid(x_)))",
                            Into::<usize>::into(num)
                        ))
                        .unwrap()
                        .to_pattern(),
                    );
                } else {
                    atom = atom.replace(&pfun).with(
                        parse!(&format!(
                            "Q({},mink(4,indexid(x_)))",
                            Into::<usize>::into(num)
                        ))
                        .unwrap()
                        .to_pattern(),
                    );
                }

                let pslashfun = parse!("PSlash(i_,j_)").unwrap().to_pattern();
                let pindex_num: usize = self.internal_index[0].into();
                if self.particle.0.is_antiparticle() {
                    atom = atom.replace(&pslashfun).with(
                        parse!(&format!(
                            "-Q({},mink(4,{}))*Gamma({},i_,j_)",
                            Into::<usize>::into(num),
                            pindex_num,
                            pindex_num
                        ))
                        .unwrap()
                        .to_pattern(),
                    );
                } else {
                    atom = atom.replace(&pslashfun).with(
                        parse!(&format!(
                            "Q({},mink(4,{}))*Gamma({},i_,j_)",
                            Into::<usize>::into(num),
                            pindex_num,
                            pindex_num
                        ))
                        .unwrap()
                        .to_pattern(),
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

                let (replacements_in, replacements_out) = if self.particle.0.is_antiparticle() {
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

                let atom = atom.replace_multiple(&reps);
                let color_atom = color_atom.replace_multiple(&reps);

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

                let atom = atom.replace_multiple(&indexid_reps);
                let color_atom = color_atom.replace_multiple(&indexid_reps);

                [
                    atom.replace(&parse!("indexid(x_)").unwrap().to_pattern())
                        .with(Atom::new_var(GS.x_).to_pattern()),
                    color_atom
                        .replace(&parse!("indexid(x_)").unwrap().to_pattern())
                        .with(Atom::new_var(GS.x_).to_pattern()),
                ]
            }
        }
    }
}

impl From<BareEdge> for Edge {
    fn from(value: BareEdge) -> Self {
        Self {
            edge_type: value.edge_type,
            internal_index: value.internal_index,
            name: value.name.into(),
            propagator: ArcPropagator(value.propagator),
            particle: value.particle,
            dod: -2,
            num: Atom::one(),
        }
    }
}

impl From<BareVertex> for Vertex {
    fn from(value: BareVertex) -> Self {
        Self {
            name: value.name.into(),
            vertex_info: value.vertex_info,
            dod: 0,
            num: Atom::one(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]

pub struct LoopMomentumBasis {
    #[bincode(with_serde)]
    pub tree: Option<()>,
    #[bincode(with_serde)]
    pub basis: TiVec<LoopIndex, EdgeIndex>,
    #[bincode(with_serde)]
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
        self.edge_signatures = graph.new_hedgevec(|_, _edge_index, _| LoopExtSignature {
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

    pub fn generate_loop_momentum_bases<E, V>(
        &self,
        graph: &HedgeGraph<E, V>,
    ) -> TiVec<LmbIndex, LoopMomentumBasis> {
        let loop_number = self.basis.len();
        let num_edges = graph.iter_all_edges().count();
        let external_signature_length = self.edge_signatures[EdgeIndex::from(0)].external.len();

        // the full virtual signature matrix in the form of a flattened vector
        let signature_matrix_flattened = self
            .edge_signatures
            .borrow()
            .into_iter()
            .flat_map(|(_, sig)| sig.internal.iter().map(|s| (*s as i8) as f64).collect_vec())
            .collect_vec();

        // convert to dmatrix
        let signature_matrix =
            DMatrix::from_row_slice(num_edges, loop_number, &signature_matrix_flattened);

        // the full external signature matrix in the form of a flattened vector
        let external_signature_matrix_flattened = self
            .edge_signatures
            .borrow()
            .into_iter()
            .flat_map(|(_, sig)| sig.external.iter().map(|&s| (s as i8) as f64).collect_vec())
            .collect_vec();

        // convert to dmatrix
        let external_signature_matrix = DMatrix::from_row_slice(
            num_edges,
            external_signature_length,
            &external_signature_matrix_flattened,
        );

        let possible_lmbs = graph
            .iter_all_edges()
            .filter(|(pair, _, _)| matches!(pair, HedgePair::Paired { .. }))
            .map(|(_, edge_index, _)| edge_index)
            .combinations(loop_number);

        let valid_lmbs = possible_lmbs
            .map(|basis| {
                let reduced_signature_matrix_flattened = basis
                    .iter()
                    .flat_map(|e| {
                        self.edge_signatures[*e]
                            .internal
                            .iter()
                            .map(|s| (*s as i8) as f64)
                    })
                    .collect_vec();

                (
                    basis,
                    DMatrix::from_row_slice(
                        loop_number,
                        loop_number,
                        &reduced_signature_matrix_flattened,
                    ),
                )
            })
            .filter(|(_basis, reduced_signature_matrix)| {
                reduced_signature_matrix.determinant() != 0. // nonzero determinant means valid lmb
            })
            .map(|(basis, reduced_signature_matrix)| {
                let mut sorted_basis = basis;
                sorted_basis.sort();
                (
                    Into::<TiVec<LoopIndex, EdgeIndex>>::into(sorted_basis),
                    reduced_signature_matrix,
                )
            });

        let lmbs = valid_lmbs
            .map(|(basis, reduced_signature_matrix)| {
                // construct the signatures

                // for this we need the reduced external signatures
                let reduced_external_signatures_vec = basis
                    .iter()
                    .flat_map(|&e| {
                        self.edge_signatures[e]
                            .external
                            .iter()
                            .map(|s| (*s as i8) as f64)
                    })
                    .collect_vec();

                let reduced_external_signatures = DMatrix::from_row_slice(
                    loop_number,
                    external_signature_length,
                    &reduced_external_signatures_vec,
                );

                let reduced_signature_matrix_inverse =
                    reduced_signature_matrix.clone().try_inverse().unwrap();

                let new_virtual_signatures =
                    signature_matrix.clone() * reduced_signature_matrix_inverse.clone();
                let new_external_signatures = external_signature_matrix.clone()
                    - signature_matrix.clone()
                        * reduced_signature_matrix_inverse.clone()
                        * reduced_external_signatures;

                let new_signatures = graph.new_hedgevec(|_, edge_id, _| {
                    let new_virtual_signature = new_virtual_signatures
                        .row(edge_id.into())
                        .iter()
                        .map(|s| s.round() as i8)
                        .collect();
                    let new_external_signature = new_external_signatures
                        .row(edge_id.into())
                        .iter()
                        .map(|s| s.round() as i8)
                        .collect();

                    LoopExtSignature {
                        internal: new_virtual_signature,
                        external: new_external_signature,
                    }
                });

                LoopMomentumBasis {
                    tree: None,
                    basis,
                    edge_signatures: new_signatures,
                }
            })
            .collect();

        lmbs
    }
}

#[derive(
    Debug,
    Clone,
    Serialize,
    Deserialize,
    Encode,
    Decode,
    Copy,
    Hash,
    From,
    Into,
    Eq,
    PartialEq,
    Ord,
    PartialOrd,
)]
pub struct LmbIndex(usize);

#[derive(Debug, Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct Vertex {
    // #[bincode(with_serde)]
    pub name: String,
    pub vertex_info: VertexInfo,
    pub num: Atom,
    pub dod: i32,
}

impl Vertex {
    pub fn generate_vertex_slots(&self, shifts: Shifts, model: &Model) -> (VertexSlots, Shifts) {
        self.vertex_info.generate_vertex_slots(shifts, model)
    }
}

impl From<BareGraph> for HedgeGraph<Edge, Vertex> {
    fn from(value: BareGraph) -> Self {
        let mut node_map = AHashMap::default();
        let mut external_nodes_to_be_dropped = AHashSet::default();
        let mut builder = HedgeGraphBuilder::new();

        for (bare_node_id, bare_node) in value.vertices.into_iter().enumerate() {
            let is_external = bare_node.edges.len() == 1;
            if is_external {
                external_nodes_to_be_dropped.insert(bare_node_id);
                continue;
            }

            let node_id = builder.add_node(bare_node.into());
            node_map.insert(bare_node_id, node_id);
        }

        for bare_edge in value.edges.into_iter() {
            match bare_edge.edge_type {
                EdgeType::Incoming => {
                    assert!(external_nodes_to_be_dropped.contains(&bare_edge.vertices[0]));
                    builder.add_external_edge(
                        node_map[&bare_edge.vertices[1]],
                        bare_edge.into(),
                        Orientation::Default,
                        Flow::Sink,
                    );
                }
                EdgeType::Outgoing => {
                    assert!(external_nodes_to_be_dropped.contains(&bare_edge.vertices[1]));
                    builder.add_external_edge(
                        node_map[&bare_edge.vertices[0]],
                        bare_edge.into(),
                        Orientation::Default,
                        Flow::Source,
                    );
                }
                EdgeType::Virtual => {
                    builder.add_edge(
                        node_map[&bare_edge.vertices[0]],
                        node_map[&bare_edge.vertices[1]],
                        bare_edge.into(),
                        Orientation::Default,
                    );
                }
            }
        }

        builder.build()
    }
}

#[derive(
    Debug, Copy, Clone, PartialEq, Eq, bincode_trait_derive::Encode, bincode_trait_derive::Decode,
)]
#[trait_decode(trait = symbolica::state::HasStateMap)]
pub struct ExternalConnection {
    pub incoming_index: ExternalIndex,
    pub outgoing_index: ExternalIndex,
}

pub fn get_cff_inverse_energy_product_impl<E, V, S: SubGraph>(
    graph: &HedgeGraph<E, V>,
    subgraph: &S,
) -> Atom {
    Atom::new_num(1)
        / graph
            .iter_edges(subgraph)
            .filter_map(|(pair, edge_index, _)| match pair {
                HedgePair::Paired { .. } => {
                    Some(Atom::new_num(2) * ose_atom_from_index(edge_index))
                }
                _ => None,
            })
            .reduce(|acc, x| acc * x)
            .unwrap_or_else(|| Atom::new_num(1))
}
