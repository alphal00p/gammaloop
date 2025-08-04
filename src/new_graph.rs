use std::borrow::Borrow;

use ahash::HashSet;
use bitvec::vec::BitVec;
use color_eyre::Result;
use itertools::Itertools;
use linnet::half_edge::{
    involution::{EdgeData, EdgeIndex, EdgeVec, Flow, HedgePair},
    subgraph::{InternalSubGraph, ModifySubgraph, SubGraph, SubGraphOps},
    HedgeGraph, NodeIndex,
};
use log::debug;
use momtrop::float::MomTropFloat;

use spenso::algebra::{algebraic_traits::IsZero, complex::Complex};
// use petgraph::Direction::Outgoing;
use idenso::metric::MS;
use smartstring::{LazyCompact, SmartString};
use spenso::structure::representation::Minkowski;
use symbolica::{
    atom::{Atom, AtomCore},
    function,
    id::Replacement,
    parse,
};
use typed_index_collections::TiVec;

use crate::{
    cff::generation::ShiftRewrite,
    graph::BareGraph,
    model::{ArcParticle, EdgeSlots, VertexSlots},
    momentum::SignOrZero,
    momentum_sample::{ExternalFourMomenta, ExternalIndex, LoopMomenta},
    new_gammaloop_integrand::LmbMultiChannelingSetup,
    numerator::GlobalPrefactor,
    signature::{ExternalSignature, SignatureLike},
    utils::{external_energy_atom_from_index, ose_atom_from_index, FloatLike, F, GS},
};

pub mod global;

#[derive(Clone, Copy, bincode_trait_derive::Encode, bincode_trait_derive::Decode, Default)]
pub struct VertexOrder(pub u8);

#[derive(Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = crate::GammaLoopContext)]
pub struct Graph {
    pub multiplicity: Atom,
    pub name: String,
    pub underlying: HedgeGraph<Edge, Vertex, NumHedgeData>,
    pub loop_momentum_basis: LoopMomentumBasis,
    pub global_prefactor: GlobalPrefactor,
}

// impl Deref for Graph {
//     type Target = HedgeGraph<Edge, Vertex>;
// }

impl From<BareGraph> for Graph {
    fn from(value: BareGraph) -> Self {
        todo!()
    }
}

pub trait FeynmanGraph {
    fn get_emr_vec_cache<T: FloatLike>(
        &self,
        loop_moms: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        lmb: &LoopMomentumBasis,
    ) -> Vec<F<T>>;

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
    fn get_real_mass_vector<T: FloatLike>(&self) -> EdgeVec<F<T>>;
    fn get_energy_cache<T: FloatLike>(
        &self,
        loop_moms: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        lmb: &LoopMomentumBasis,
    ) -> EdgeVec<F<T>>;
    fn get_esurface_canonization(&self, lmb: &LoopMomentumBasis) -> Option<ShiftRewrite>;
    fn external_in_or_out_signature(&self) -> ExternalSignature;
    fn get_external_partcles(&self) -> Vec<ArcParticle>;
    fn get_external_signature(&self) -> SignatureLike<ExternalIndex>;
    fn get_energy_atoms(&self) -> Vec<Atom>;
    fn get_external_energy_atoms(&self) -> Vec<Atom>;
    fn explicit_ose_atom(&self, edge: EdgeIndex) -> Atom;
    fn get_ose_replacements(&self) -> Vec<Replacement>;
    fn expected_scale(&self, e_cm: F<f64>) -> F<f64>;
    fn dummy_list(&self) -> Vec<EdgeIndex>;
    fn no_dummy(&self) -> BitVec;
}

impl FeynmanGraph for HedgeGraph<Edge, Vertex, NumHedgeData> {
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
        let mom = parse!(&format!("Q{}", Into::<usize>::into(edge)));
        let mass = self[edge]
            .particle
            .0
            .mass
            .expression
            .clone()
            .unwrap_or(Atom::num(0));

        (mom, mass)
    }

    fn substitute_lmb(&self, edge: EdgeIndex, atom: Atom, lmb: &LoopMomentumBasis) -> Atom {
        let num = Into::<usize>::into(edge);
        let mom = parse!(&format!("Q({num},x_)")).to_pattern();
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

        self.iter_edges_of(&node_hairs)
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
        get_cff_inverse_energy_product_impl(self, &full_subgraph, &[])
    }

    fn get_loop_number(&self) -> usize {
        let internal_subgraph =
            InternalSubGraph::cleaned_filter_pessimist(self.full_filter(), self);
        self.cyclotomatic_number(&internal_subgraph)
    }

    fn get_real_mass_vector<T: FloatLike>(&self) -> EdgeVec<F<T>> {
        self.new_edgevec(|edge, _edge_id, _| match edge.particle.0.mass.value {
            Some(mass) => F::from_ff64(mass.re),
            None => F::from_f64(0.0),
        })
    }

    fn get_energy_cache<T: FloatLike>(
        &self,
        loop_moms: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        lmb: &LoopMomentumBasis,
    ) -> EdgeVec<F<T>> {
        self.new_edgevec_from_iter(
            lmb.edge_signatures
                .borrow()
                .into_iter()
                .map(|(_, sig)| sig.compute_four_momentum_from_three(loop_moms, external_moms))
                .zip(self.iter_edges())
                .map(|(emr_mom, (p, _, edge))| {
                    if p.is_paired() {
                        emr_mom
                            .spatial
                            .on_shell_energy(edge.data.particle.0.mass.value.map(|m| {
                                if m.im.is_non_zero() {
                                    panic!("Complex masses not yet supported in gammaLoop")
                                }
                                F::<T>::from_ff64(m.re)
                            }))
                            .value
                    } else {
                        emr_mom.temporal.value // a wierd way of just obtaining the energy of the external particles
                    }
                }),
        )
        .unwrap()
    }

    fn get_emr_vec_cache<T: FloatLike>(
        &self,
        loop_moms: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        lmb: &LoopMomentumBasis,
    ) -> Vec<F<T>> {
        lmb.edge_signatures
            .iter()
            .zip(self.iter_edges())
            .filter(|(_, (pair, _, _))| matches!(pair, HedgePair::Paired { .. }))
            .map(|((_, sig), _)| sig.compute_three_momentum_from_four(loop_moms, external_moms))
            .flat_map(|emr_vec| [emr_vec.px, emr_vec.py, emr_vec.pz])
            .collect()
    }

    fn get_esurface_canonization(&self, lmb: &LoopMomentumBasis) -> Option<ShiftRewrite> {
        // let external_edges: TiVec<ExternalIndex, _> =
        // self
        // .iter_edges()
        // .filter(|(pair, _, _)| matches!(pair, HedgePair::Unpaired { .. }))
        // .collect();

        // find the external leg which does not appear in it's own signature
        lmb.ext_edges
            .iter_enumerated()
            .find(|(external_index, edge_id)| {
                lmb.edge_signatures[**edge_id].external[*external_index] == SignOrZero::Zero
            })
            .map(|(_, dep_mom_edge_id)| {
                let dep_mom_signatrue = &lmb.edge_signatures[*dep_mom_edge_id].external;

                let external_shift = lmb
                    .ext_edges
                    .iter()
                    .zip(dep_mom_signatrue)
                    .filter(|(_, dep_mom_sign)| dep_mom_sign.is_sign())
                    .map(|(external_edge, dep_mom_sign)| (*external_edge, dep_mom_sign as i64))
                    .collect_vec();

                ShiftRewrite {
                    dependent_momentum: *dep_mom_edge_id,
                    dependent_momentum_expr: external_shift,
                }
            })
    }

    fn external_in_or_out_signature(&self) -> ExternalSignature {
        self.iter_edges()
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
        self.iter_edges()
            .filter_map(|(pair, _, data)| match pair {
                HedgePair::Unpaired { .. } => Some(data.data.particle.clone()),
                _ => None,
            })
            .collect()
    }

    fn get_external_signature(&self) -> SignatureLike<ExternalIndex> {
        SignatureLike::from_iter(self.iter_edges().filter_map(|(pair, _, _)| match pair {
            HedgePair::Unpaired { flow, .. } => match flow {
                Flow::Source => Some(SignOrZero::Minus),
                Flow::Sink => Some(SignOrZero::Plus),
            },
            _ => None,
        }))
    }

    fn get_energy_atoms(&self) -> Vec<Atom> {
        self.iter_edges()
            .map(|(pair, edge_id, _)| match pair {
                HedgePair::Paired { .. } => ose_atom_from_index(edge_id),
                HedgePair::Unpaired { .. } => external_energy_atom_from_index(edge_id),
                _ => unreachable!(),
            })
            .collect_vec()
    }

    fn explicit_ose_atom(&self, edge: EdgeIndex) -> Atom {
        let mass = self[edge].particle.symbolic_mass();
        let mass2 = &mass * &mass;

        let emr_vec = function!(GS.emr_vec, Into::<usize>::into(edge) as i64);
        let dot = function!(MS.dot, emr_vec, emr_vec);

        (mass2 - dot).sqrt()
    }

    fn get_ose_replacements(&self) -> Vec<Replacement> {
        self.iter_edges()
            .filter(|(pair, _, _)| matches!(pair, HedgePair::Paired { .. }))
            .map(|(_, edge_id, _)| {
                let ose_atom = ose_atom_from_index(edge_id);
                let explicit = self.explicit_ose_atom(edge_id);
                Replacement::new(ose_atom.to_pattern(), explicit.to_pattern())
            })
            .collect()
    }

    fn get_external_energy_atoms(&self) -> Vec<Atom> {
        self.iter_edges()
            .filter_map(|(pair, edge_id, _)| match pair {
                HedgePair::Unpaired { .. } => Some(external_energy_atom_from_index(edge_id)),
                _ => None,
            })
            .collect_vec()
    }

    fn expected_scale(&self, e_cm: F<f64>) -> F<f64> {
        let mut scale = Complex::new_re(F(1.0));
        for (_, _, vertex) in self.iter_nodes() {
            // include the values of all couplings
            let coupling_value = vertex
                .vertex_rule
                .couplings
                .iter()
                .flat_map(|couplings| {
                    couplings.iter().map(|coupling| {
                        coupling
                            .as_ref()
                            .map(|coupling| {
                                coupling
                                    .value
                                    .map(|x| Complex::new(F(x.re), F(x.im)))
                                    .unwrap_or(Complex::new_re(F(1.0)))
                            })
                            .unwrap_or(Complex::new_re(F(1.0)))
                    })
                })
                .fold(Complex::new_re(F(1.0)), |product, term| product * term);

            scale *= Complex::new_re(e_cm.powi(vertex.dod)) * coupling_value;
        }

        for (_, _, edge_data) in
            self.iter_edges_of(&self.full_filter().subtract(&self.external_filter()))
        {
            scale *= Complex::new_re(e_cm.powi(edge_data.data.dod))
        }

        scale.norm_squared().sqrt()
    }

    fn no_dummy(&self) -> BitVec {
        let mut subgraph = self.full_filter();
        for (hedge_pair, _, edge) in self.iter_edges() {
            if edge.data.is_dummy {
                subgraph.sub(hedge_pair)
            }
        }
        subgraph
    }

    fn dummy_list(&self) -> Vec<EdgeIndex> {
        self.iter_edges()
            .filter_map(|(hedge_pair, edge_index, edge_data)| {
                if edge_data.data.is_dummy {
                    Some(edge_index)
                } else {
                    None
                }
            })
            .collect()
    }
}

impl Graph {
    pub fn new(
        name: SmartString<LazyCompact>,
        multiplicity: Atom,
        underlying: HedgeGraph<Edge, Vertex, NumHedgeData>,
    ) -> Result<Self> {
        Ok(Self {
            name: name.to_string(),
            multiplicity,
            loop_momentum_basis: underlying.lmb(&underlying.full_filter()),
            underlying,
            global_prefactor: GlobalPrefactor::default(),
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
                .loop_edges
                .iter()
                .filter(|&edge_id| self.underlying[*edge_id].particle.0.is_massless())
                .collect_vec();

            if massless_edges.is_empty() {
                continue;
            }

            if channels.iter().any(|channel| {
                let basis_of_included_channel = &lmbs[*channel].loop_edges;
                massless_edges
                    .iter()
                    .all(|edge_id| basis_of_included_channel.contains(edge_id))
            }) {
                continue;
            }

            // only for 1 to n for now, assuming center of mass
            if channels.iter().any(|channel| {
                let massless_edges_of_included_channel = lmbs[*channel]
                    .loop_edges
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
            }) {
                continue;
            }

            channels.push(lmb_index);
        }

        let channels: TiVec<_, _> = if !channels.is_empty() {
            channels.into_iter().sorted().collect()
        } else {
            let current_lmb_index = lmbs
                .iter_enumerated()
                .find(|(_lmb_index, lmb)| lmb.loop_edges == self.loop_momentum_basis.loop_edges)
                .map(|(lmb_index, _)| lmb_index)
                .unwrap();

            vec![current_lmb_index].into()
        };

        debug!(
            "number of lmbs: {}, number of channels: {}",
            lmbs.len(),
            channels.len()
        );

        LmbMultiChannelingSetup { channels }
    }

    pub fn iter_loop_edges(&self) -> impl Iterator<Item = (HedgePair, EdgeIndex, EdgeData<&Edge>)> {
        self.underlying.iter_edges().filter(|(_, edge_index, _)| {
            self.loop_momentum_basis.edge_signatures[*edge_index]
                .internal
                .iter()
                .any(|sign| sign.is_sign())
        })
    }

    pub fn iter_non_loop_edges(
        &self,
    ) -> impl Iterator<Item = (HedgePair, EdgeIndex, EdgeData<&Edge>)> {
        self.underlying.iter_edges().filter(|(_, edge_index, _)| {
            self.loop_momentum_basis.edge_signatures[*edge_index]
                .internal
                .iter()
                .all(|sign| sign.is_zero())
        })
    }
}

impl Graph {
    // pub fn apply_vertex_rule(&self, node_id: NodeIndex) -> Option<[DataTensor<Atom>; 3]> {
    //     // self.underlying[node_id].vertex_info.apply_vertex_rule(
    //     //     &self.underlying.add_signs_to_edges(node_id),
    //     //     Into::<usize>::into(node_id),
    //     //     &self.vertex_slots[node_id],
    //     // )
    // }
}

pub mod edge;
pub mod parse;
pub use edge::Edge;
pub mod hedge_data;
pub use hedge_data::NumHedgeData;
pub mod vertex;
pub use vertex::Vertex;

pub mod lmb;
pub use lmb::{LMBext, LmbIndex, LoopMomentumBasis};

// impl From<BareGraph> for HedgeGraph<Edge, Vertex, VertexOrder> {
//     fn from(value: BareGraph) -> Self {
//         let mut node_map = AHashMap::default();
//         let mut external_nodes_to_be_dropped = AHashSet::default();
//         let mut builder = HedgeGraphBuilder::new();

//         for (bare_node_id, bare_node) in value.vertices.iter().enumerate() {
//             let is_external = bare_node.edges.len() == 1;
//             if is_external {
//                 external_nodes_to_be_dropped.insert(bare_node_id);
//                 continue;
//             }

//             let node_id = builder.add_node(bare_node.clone().into());
//             node_map.insert(bare_node_id, node_id);
//         }

//         for (bare_edge_id, bare_edge) in value.edges.into_iter().enumerate() {
//             match bare_edge.edge_type {
//                 EdgeType::Incoming => {
//                     assert!(external_nodes_to_be_dropped.contains(&bare_edge.vertices[0]));
//                     builder.add_external_edge(
//                         node_map[&bare_edge.vertices[1]].add_data(VertexOrder(0)),
//                         bare_edge.into(),
//                         Orientation::Default,
//                         Flow::Sink,
//                     );
//                 }
//                 EdgeType::Outgoing => {
//                     assert!(external_nodes_to_be_dropped.contains(&bare_edge.vertices[1]));
//                     builder.add_external_edge(
//                         node_map[&bare_edge.vertices[0]].add_data(VertexOrder(0)),
//                         bare_edge.into(),
//                         Orientation::Default,
//                         Flow::Source,
//                     );
//                 }
//                 EdgeType::Virtual => {
//                     let source = &bare_edge.vertices[0];
//                     let (source_pos, _) = value.vertices[*source]
//                         .edges
//                         .iter()
//                         .find_position(|i| &bare_edge_id == *i)
//                         .unwrap();
//                     let sink = &bare_edge.vertices[1];
//                     let (sink_pos, _) = value.vertices[*sink]
//                         .edges
//                         .iter()
//                         .find_position(|i| &bare_edge_id == *i)
//                         .unwrap();

//                     builder.add_edge(
//                         node_map[&bare_edge.vertices[0]].add_data(VertexOrder(source_pos as u8)),
//                         node_map[&bare_edge.vertices[1]].add_data(VertexOrder(sink_pos as u8)),
//                         bare_edge.into(),
//                         Orientation::Default,
//                     );
//                 }
//             }
//         }

//         builder.build()
//     }
// }

#[derive(
    Debug, Copy, Clone, PartialEq, Eq, bincode_trait_derive::Encode, bincode_trait_derive::Decode,
)]
#[trait_decode(trait = symbolica::state::HasStateMap)]
pub struct ExternalConnection {
    pub incoming_index: ExternalIndex,
    pub outgoing_index: ExternalIndex,
}

pub fn get_cff_inverse_energy_product_impl<E, V, H, S: SubGraph>(
    graph: &HedgeGraph<E, V, H>,
    subgraph: &S,
    contract_edges: &[EdgeIndex],
) -> Atom {
    Atom::num(1)
        / graph
            .iter_edges_of(subgraph)
            .filter_map(|(pair, edge_index, _)| match pair {
                HedgePair::Paired { .. } => {
                    if contract_edges.contains(&edge_index) {
                        None
                    } else {
                        Some(-Atom::num(2) * ose_atom_from_index(edge_index))
                    }
                }
                _ => None,
            })
            .reduce(|acc, x| acc * x)
            .unwrap_or_else(|| Atom::num(1))
}
