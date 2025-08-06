use std::{borrow::Borrow, ops::Deref};

use bitvec::vec::BitVec;
use color_eyre::Result;
use itertools::Itertools;
use linnet::half_edge::{
    involution::{EdgeIndex, EdgeVec, Flow, HedgePair},
    subgraph::{InternalSubGraph, ModifySubgraph, SubGraphOps},
    HedgeGraph, NodeIndex,
};
use log::debug;
use momtrop::float::MomTropFloat;

use spenso::{
    algebra::{algebraic_traits::IsZero, complex::Complex},
    structure::concrete_index::ExpandedIndex,
};
// use petgraph::Direction::Outgoing;
use idenso::metric::MS;
use symbolica::{
    atom::{Atom, AtomCore},
    function,
    id::Replacement,
    parse,
};

use crate::{
    cff::generation::ShiftRewrite,
    model::ArcParticle,
    momentum::SignOrZero,
    momentum_sample::{ExternalFourMomenta, ExternalIndex, LoopMomenta},
    signature::{ExternalSignature, SignatureLike},
    utils::{external_energy_atom_from_index, ose_atom_from_index, FloatLike, F, GS},
    uv::uv_graph::UVE,
};

use super::{
    get_cff_inverse_energy_product_impl, Edge, Graph, LoopMomentumBasis, NumHedgeData, Vertex,
};

pub trait FeynmanGraph {
    fn get_emr_vec_cache<T: FloatLike>(
        &self,
        loop_moms: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        lmb: &LoopMomentumBasis,
    ) -> Vec<F<T>>;

    fn num_virtual_edges(&self, subgraph: BitVec) -> usize;
    fn is_incoming_to(&self, edge: EdgeIndex, vertex: NodeIndex) -> bool;
    // fn denominator(&self, edge: EdgeIndex) -> (Atom, Atom);
    fn substitute_lmb(&self, edge: EdgeIndex, atom: Atom, lmb: &LoopMomentumBasis) -> Atom;

    // fn get_local_edge_position(
    //     &self,
    //     node_id: NodeIndex,
    //     edge_id: EdgeIndex,
    //     skip_one: bool,
    // ) -> usize;
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

impl Deref for Graph {
    type Target = HedgeGraph<Edge, Vertex, NumHedgeData>;

    fn deref(&self) -> &Self::Target {
        &self.underlying
    }
}

impl FeynmanGraph for Graph {
    fn add_signs_to_edges(&self, node_id: NodeIndex) -> Vec<isize> {
        self.underlying.add_signs_to_edges(node_id)
    }

    fn dummy_list(&self) -> Vec<EdgeIndex> {
        self.underlying.dummy_list()
    }

    fn expected_scale(&self, e_cm: F<f64>) -> F<f64> {
        self.underlying.expected_scale(e_cm)
    }

    fn explicit_ose_atom(&self, edge: EdgeIndex) -> Atom {
        self.underlying.explicit_ose_atom(edge)
    }

    fn external_in_or_out_signature(&self) -> ExternalSignature {
        self.underlying.external_in_or_out_signature()
    }

    fn get_cff_inverse_energy_product(&self) -> Atom {
        self.underlying.get_cff_inverse_energy_product()
    }

    fn get_emr_vec_cache<T: FloatLike>(
        &self,
        loop_moms: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        lmb: &LoopMomentumBasis,
    ) -> Vec<F<T>> {
        self.underlying
            .get_emr_vec_cache(loop_moms, external_moms, lmb)
    }

    fn get_energy_atoms(&self) -> Vec<Atom> {
        self.underlying.get_energy_atoms()
    }

    fn get_energy_cache<T: FloatLike>(
        &self,
        loop_moms: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        lmb: &LoopMomentumBasis,
    ) -> EdgeVec<F<T>> {
        self.underlying
            .get_energy_cache(loop_moms, external_moms, lmb)
    }

    fn get_esurface_canonization(&self, lmb: &LoopMomentumBasis) -> Option<ShiftRewrite> {
        self.underlying.get_esurface_canonization(lmb)
    }

    fn get_external_energy_atoms(&self) -> Vec<Atom> {
        self.underlying.get_external_energy_atoms()
    }

    fn get_external_partcles(&self) -> Vec<ArcParticle> {
        self.underlying.get_external_partcles()
    }

    fn get_external_signature(&self) -> SignatureLike<ExternalIndex> {
        self.underlying.get_external_signature()
    }

    // fn get_local_edge_position(
    //     &self,
    //     node_id: NodeIndex,
    //     edge_id: EdgeIndex,
    //     skip_one: bool,
    // ) -> usize {
    //     self.underlying
    //         .get_local_edge_position(node_id, edge_id, skip_one)
    // }

    fn get_loop_number(&self) -> usize {
        self.underlying.get_loop_number()
    }

    fn get_ose_replacements(&self) -> Vec<Replacement> {
        self.underlying.get_ose_replacements()
    }

    fn get_real_mass_vector<T: FloatLike>(&self) -> EdgeVec<F<T>> {
        self.underlying.get_real_mass_vector()
    }
    fn is_incoming_to(&self, edge: EdgeIndex, vertex: NodeIndex) -> bool {
        self.underlying.is_incoming_to(edge, vertex)
    }

    fn no_dummy(&self) -> BitVec {
        self.underlying.no_dummy()
    }

    fn num_virtual_edges(&self, subgraph: BitVec) -> usize {
        self.underlying.num_virtual_edges(subgraph)
    }

    fn substitute_lmb(&self, edge: EdgeIndex, atom: Atom, lmb: &LoopMomentumBasis) -> Atom {
        self.underlying.substitute_lmb(edge, atom, lmb)
    }
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

    // fn denominator(&self, edge: EdgeIndex) -> (Atom, Atom) {
    //     let mom = parse!(&format!("Q{}", Into::<usize>::into(edge)));
    //     let mass = self[edge]
    //         .particle
    //         .0
    //         .mass
    //         .expression
    //         .clone()
    //         .unwrap_or(Atom::num(0));

    //     (mom, mass)
    // }

    fn substitute_lmb(&self, edge: EdgeIndex, atom: Atom, lmb: &LoopMomentumBasis) -> Atom {
        let num = Into::<usize>::into(edge);
        let mom = parse!(&format!("Q({num},x_)")).to_pattern();
        let mom_rep = lmb.pattern(edge);
        atom.replace(mom).with(mom_rep)
    }

    // fn get_local_edge_position(
    //     &self,
    //     node_id: NodeIndex,
    //     edge_id: EdgeIndex,
    //     skip_one: bool,
    // ) -> usize {
    //     unimplemented!()
    // }

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
        self.new_edgevec(|edge, _edge_id, _| match edge.mass_value() {
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
                            .on_shell_energy(edge.data.mass_value().map(|m| {
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
                HedgePair::Unpaired { .. } => data.data.particle(),
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
        let mass = self[edge].mass_atom();
        let mass2 = &mass * &mass;

        let dot = GS.emr_mom(edge, Atom::from(ExpandedIndex::from_iter([1])))
            * GS.emr_mom(edge, Atom::from(ExpandedIndex::from_iter([1])))
            + GS.emr_mom(edge, Atom::from(ExpandedIndex::from_iter([2])))
                * GS.emr_mom(edge, Atom::from(ExpandedIndex::from_iter([2])))
            + GS.emr_mom(edge, Atom::from(ExpandedIndex::from_iter([3])))
                * GS.emr_mom(edge, Atom::from(ExpandedIndex::from_iter([3])));

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
                .as_ref()
                .map(|r| {
                    r.couplings
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
                        .fold(Complex::new_re(F(1.0)), |product, term| product * term)
                })
                .unwrap_or(Complex::new_re(F(1.0)));

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
