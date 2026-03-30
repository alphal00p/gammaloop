use std::{borrow::Borrow, ops::Deref};

use itertools::Itertools;
use linnet::half_edge::{
    HedgeGraph, NodeIndex,
    involution::{EdgeData, EdgeIndex, EdgeVec, Flow, HedgePair},
    subgraph::{
        HedgeNode, InternalSubGraph, ModifySubSet, OrientedCut, SuBitGraph, SubGraphLike,
        SubSetLike, SubSetOps,
    },
};

use spenso::{
    algebra::{algebraic_traits::IsZero, complex::Complex},
    structure::concrete_index::ExpandedIndex,
};
// use petgraph::Direction::Outgoing;
use symbolica::{
    atom::{Atom, AtomCore},
    id::Replacement,
};
use typed_index_collections::TiVec;

use crate::{
    cff::generation::ShiftRewrite,
    integrands::process::param_builder::{ParamBuilderGraph, SplitPolarizations},
    model::{ArcParticle, Model},
    momentum::sample::{ExternalFourMomenta, ExternalIndex, LoopMomenta},
    momentum::signature::{ExternalSignature, SignatureLike},
    momentum::{PolDef, SignOrZero},
    numerator::graph::ReversibleEdge,
    utils::{F, FloatLike, GS, external_energy_atom_from_index, ose_atom_from_index},
    uv::uv_graph::UVE,
};

use super::{
    Edge, Graph, LoopMomentumBasis, NumHedgeData, Vertex, get_cff_inverse_energy_product_impl,
};

pub trait FeynmanGraph {
    fn get_emr_vec_cache<T: FloatLike>(
        &self,
        loop_moms: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        lmb: &LoopMomentumBasis,
    ) -> Vec<F<T>>;

    fn num_virtual_edges(&self, subgraph: SuBitGraph) -> usize;
    fn is_incoming_to(&self, edge: EdgeIndex, vertex: NodeIndex) -> bool;
    // fn denominator(&self, edge: EdgeIndex) -> (Atom, Atom);
    // fn get_local_edge_position(
    //     &self,
    //     node_id: NodeIndex,
    //     edge_id: EdgeIndex,
    //     skip_one: bool,
    // ) -> usize;
    fn add_signs_to_edges(&self, node_id: NodeIndex) -> Vec<isize>;
    fn get_cff_inverse_energy_product(&self) -> Atom;
    fn get_loop_number(&self) -> usize;
    fn get_real_mass_vector<T: FloatLike>(&self, model: &Model) -> EdgeVec<F<T>>;
    fn get_external_masses<T: FloatLike>(&self, model: &Model) -> TiVec<ExternalIndex, F<T>>;
    fn get_energy_cache<T: FloatLike>(
        &self,
        model: &Model,

        loop_moms: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        lmb: &LoopMomentumBasis,
    ) -> EdgeVec<F<T>>;
    fn get_esurface_canonization(&self, lmb: &LoopMomentumBasis) -> Option<ShiftRewrite>;
    fn external_in_or_out_signature(&self) -> ExternalSignature;

    fn get_external_partcles(&self) -> Vec<ArcParticle>;
    fn get_external_signature(&self) -> SignatureLike<ExternalIndex>;
    fn get_energy_atoms(&self) -> Vec<Atom>;
    // fn get_external_energy_atoms(&self) -> Vec<Atom>;
    // fn explicit_ose_atom(&self, edge: EdgeIndex) -> Atom;
    // fn get_ose_replacements(&self) -> Vec<Replacement>;
    fn expected_scale(&self, e_cm: F<f64>, model: &Model) -> F<f64>;
    fn dummy_list(&self) -> Vec<EdgeIndex>;
    fn no_dummy(&self) -> SuBitGraph;
    fn all_st_cuts_for_cs(
        &self,
        source_nodes: HedgeNode,
        target_nodes: HedgeNode,
        initial_state_tree: &SuBitGraph,
    ) -> Vec<(SuBitGraph, OrientedCut, SuBitGraph)>;
}

impl Deref for Graph {
    type Target = HedgeGraph<Edge, Vertex, NumHedgeData>;

    fn deref(&self) -> &Self::Target {
        &self.underlying
    }
}

impl SplitPolarizations for Graph {
    fn polarizations(&self) -> Vec<Atom> {
        self.polarizations.iter().map(|a| a.1.clone()).collect()
    }
}

impl<'a, V, E: UVE> SplitPolarizations
    for (&'a Vec<(PolDef, Atom)>, &'a HedgeGraph<E, V, NumHedgeData>)
{
    fn polarizations(&self) -> Vec<Atom> {
        self.0.iter().map(|a| a.1.clone()).collect()
    }
}

impl<'a, V, E: UVE> ParamBuilderGraph
    for (&'a Vec<(PolDef, Atom)>, &'a HedgeGraph<E, V, NumHedgeData>)
where
    for<'b> EdgeData<&'b E>: ReversibleEdge,
{
    fn iter_edge_ids(&self) -> impl Iterator<Item = EdgeIndex> + '_ {
        self.1.iter_edge_ids()
    }

    fn get_external_energy_atoms(&self) -> Vec<Atom> {
        self.1.get_external_energy_atoms()
    }

    fn get_ose_replacements(&self) -> Vec<Replacement> {
        self.1.get_ose_replacements()
    }

    fn explicit_ose_atom(&self, edge: EdgeIndex) -> Atom {
        self.1.explicit_ose_atom(edge)
    }

    fn loop_mom_params(&self, lmb: &LoopMomentumBasis) -> Vec<Atom> {
        lmb.loop_edges
            .iter()
            .flat_map(|edge_id| {
                vec![
                    GS.emr_mom(*edge_id, Atom::from(ExpandedIndex::from_iter([1]))),
                    GS.emr_mom(*edge_id, Atom::from(ExpandedIndex::from_iter([2]))),
                    GS.emr_mom(*edge_id, Atom::from(ExpandedIndex::from_iter([3]))),
                ]
            })
            .collect()
    }

    fn external_spatial_params(&self) -> Vec<Atom> {
        self.1.external_spatial_params()
    }
}

impl<V, E: UVE> ParamBuilderGraph for HedgeGraph<E, V, NumHedgeData>
where
    for<'a> EdgeData<&'a E>: ReversibleEdge,
{
    fn iter_edge_ids(&self) -> impl Iterator<Item = EdgeIndex> + '_ {
        self.iter_edges().map(|(_, e, _)| e)
    }

    fn get_external_energy_atoms(&self) -> Vec<Atom> {
        self.iter_edges()
            .filter_map(|(pair, edge_id, _)| match pair {
                HedgePair::Unpaired { .. } => Some(external_energy_atom_from_index(edge_id)),
                _ => None,
            })
            .collect_vec()
    }

    fn explicit_ose_atom(&self, edge: EdgeIndex) -> Atom {
        let mass = self[edge].mass_atom();
        let mass2 = &mass * &mass;

        // println!("{}", mass2);

        let dot = GS.emr_mom(edge, Atom::from(ExpandedIndex::from_iter([1])))
            * GS.emr_mom(edge, Atom::from(ExpandedIndex::from_iter([1])))
            + GS.emr_mom(edge, Atom::from(ExpandedIndex::from_iter([2])))
                * GS.emr_mom(edge, Atom::from(ExpandedIndex::from_iter([2])))
            + GS.emr_mom(edge, Atom::from(ExpandedIndex::from_iter([3])))
                * GS.emr_mom(edge, Atom::from(ExpandedIndex::from_iter([3])));

        (dot + mass2).sqrt()
    }

    fn loop_mom_params(&self, lmb: &LoopMomentumBasis) -> Vec<Atom> {
        lmb.loop_edges
            .iter()
            .flat_map(|edge_id| {
                vec![
                    GS.emr_mom(*edge_id, Atom::from(ExpandedIndex::from_iter([1]))),
                    GS.emr_mom(*edge_id, Atom::from(ExpandedIndex::from_iter([2]))),
                    GS.emr_mom(*edge_id, Atom::from(ExpandedIndex::from_iter([3]))),
                ]
            })
            .collect()
    }

    fn external_spatial_params(&self) -> Vec<Atom> {
        self.iter_edges()
            .flat_map(|(pair, edge_id, _)| {
                if let HedgePair::Unpaired { .. } = pair {
                    vec![
                        GS.emr_mom(edge_id, Atom::from(ExpandedIndex::from_iter([1]))),
                        GS.emr_mom(edge_id, Atom::from(ExpandedIndex::from_iter([2]))),
                        GS.emr_mom(edge_id, Atom::from(ExpandedIndex::from_iter([3]))),
                    ]
                } else {
                    vec![]
                }
            })
            .collect()
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
}

impl ParamBuilderGraph for Graph {
    // fn polarizations(&self) -> Vec<Atom> {
    //     self.polarizations.iter().map(|a| a.1.clone()).collect()
    // }

    fn iter_edge_ids(&self) -> impl Iterator<Item = EdgeIndex> + '_ {
        self.underlying.iter_edge_ids()
    }

    fn get_external_energy_atoms(&self) -> Vec<Atom> {
        if self.initial_state_cut.nedges(&self.underlying) == 0 {
            self.underlying.get_external_energy_atoms()
        } else {
            self.underlying
                .iter_edges_of(&self.initial_state_cut)
                .sorted_by(|a, b| a.1.cmp(&b.1))
                .map(|(_, edge_id, _)| external_energy_atom_from_index(edge_id))
                .collect()
        }
    }

    fn get_ose_replacements(&self) -> Vec<Replacement> {
        if self.initial_state_cut.nedges(&self.underlying) == 0 {
            self.underlying.get_ose_replacements()
        } else {
            let underlying_without_is_cut = self
                .underlying
                .full_filter()
                .subtract(&self.initial_state_cut.left)
                .subtract(&self.initial_state_cut.right);

            self.underlying
                .iter_edges_of(&underlying_without_is_cut)
                .map(|(_, edge_id, _)| {
                    let ose_atom = ose_atom_from_index(edge_id);
                    let explicit = self.underlying.explicit_ose_atom(edge_id);
                    Replacement::new(ose_atom.to_pattern(), explicit.to_pattern())
                })
                .collect()
        }
    }

    fn explicit_ose_atom(&self, edge: EdgeIndex) -> Atom {
        self.underlying.explicit_ose_atom(edge)
    }

    #[allow(unused_variables)]
    fn loop_mom_params(&self, lmb: &LoopMomentumBasis) -> Vec<Atom> {
        self.loop_momentum_basis
            .loop_edges
            .iter()
            .flat_map(|edge_id| {
                vec![
                    GS.emr_mom(*edge_id, Atom::from(ExpandedIndex::from_iter([1]))),
                    GS.emr_mom(*edge_id, Atom::from(ExpandedIndex::from_iter([2]))),
                    GS.emr_mom(*edge_id, Atom::from(ExpandedIndex::from_iter([3]))),
                ]
            })
            .collect()
    }

    fn external_spatial_params(&self) -> Vec<Atom> {
        if self.initial_state_cut.nedges(&self.underlying) == 0 {
            self.underlying.external_spatial_params()
        } else {
            self.underlying
                .iter_edges_of(&self.initial_state_cut)
                .sorted_by(|a, b| a.1.cmp(&b.1))
                .flat_map(|(_, edge_id, _)| {
                    vec![
                        GS.emr_mom(edge_id, Atom::from(ExpandedIndex::from_iter([1]))),
                        GS.emr_mom(edge_id, Atom::from(ExpandedIndex::from_iter([2]))),
                        GS.emr_mom(edge_id, Atom::from(ExpandedIndex::from_iter([3]))),
                    ]
                })
                .collect()
        }
    }
}

// impl FeynmanGraph for Graph {
//     fn add_signs_to_edges(&self, node_id: NodeIndex) -> Vec<isize> {
//         self.underlying.add_signs_to_edges(node_id)
//     }

//     fn dummy_list(&self) -> Vec<EdgeIndex> {
//         self.underlying.dummy_list()
//     }

//     fn expected_scale(&self, e_cm: F<f64>, model: &Model) -> F<f64> {
//         self.underlying.expected_scale(e_cm, model)
//     }

//     // fn explicit_ose_atom(&self, edge: EdgeIndex) -> Atom {
//     //     self.underlying.explicit_ose_atom(edge)
//     // }

//     fn external_in_or_out_signature(&self) -> ExternalSignature {
//         self.underlying.external_in_or_out_signature()
//     }

//     fn get_cff_inverse_energy_product(&self) -> Atom {
//         self.underlying.get_cff_inverse_energy_product()
//     }

//     fn get_emr_vec_cache<T: FloatLike>(
//         &self,
//         loop_moms: &LoopMomenta<F<T>>,
//         external_moms: &ExternalFourMomenta<F<T>>,
//         lmb: &LoopMomentumBasis,
//     ) -> Vec<F<T>> {
//         self.underlying
//             .get_emr_vec_cache(loop_moms, external_moms, lmb)
//     }

//     fn get_energy_atoms(&self) -> Vec<Atom> {
//         self.underlying.get_energy_atoms()
//     }

//     fn get_energy_cache<T: FloatLike>(
//         &self,
//         model: &Model,
//         paramb: &ParamBuilder,
//         loop_moms: &LoopMomenta<F<T>>,
//         external_moms: &ExternalFourMomenta<F<T>>,
//         lmb: &LoopMomentumBasis,
//     ) -> EdgeVec<F<T>> {
//         self.underlying
//             .get_energy_cache(model, paramb, loop_moms, external_moms, lmb)
//     }

//     fn get_esurface_canonization(&self, lmb: &LoopMomentumBasis) -> Option<ShiftRewrite> {
//         self.underlying.get_esurface_canonization(lmb)
//     }

//     // fn get_external_energy_atoms(&self) -> Vec<Atom> {
//     //     self.underlying.get_external_energy_atoms()
//     // }

//     fn get_external_partcles(&self) -> Vec<ArcParticle> {
//         self.underlying.get_external_partcles()
//     }

//     fn get_external_signature(&self) -> SignatureLike<ExternalIndex> {
//         self.underlying.get_external_signature()
//     }

//     // fn get_local_edge_position(
//     //     &self,
//     //     node_id: NodeIndex,
//     //     edge_id: EdgeIndex,
//     //     skip_one: bool,
//     // ) -> usize {
//     //     self.underlying
//     //         .get_local_edge_position(node_id, edge_id, skip_one)
//     // }

//     fn get_loop_number(&self) -> usize {
//         self.underlying.get_loop_number()
//     }

//     fn get_real_mass_vector<T: FloatLike>(
//         &self,
//         model: &Model,
//         paramb: &ParamBuilder,
//     ) -> EdgeVec<F<T>> {
//         self.underlying.get_real_mass_vector(model, paramb)
//     }
//     fn is_incoming_to(&self, edge: EdgeIndex, vertex: NodeIndex) -> bool {
//         self.underlying.is_incoming_to(edge, vertex)
//     }

//     fn no_dummy(&self) -> SuBitGraph {
//         self.underlying.no_dummy()
//     }

//     fn num_virtual_edges(&self, subgraph: SuBitGraph) -> usize {
//         self.underlying.num_virtual_edges(subgraph)
//     }

//     fn substitute_lmb(&self, edge: EdgeIndex, atom: Atom, lmb: &LoopMomentumBasis) -> Atom {
//         self.underlying.substitute_lmb(edge, atom, lmb)
//     }
// }

impl FeynmanGraph for Graph {
    fn num_virtual_edges(&self, subgraph: SuBitGraph) -> usize {
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

    // fn get_local_edge_position(
    //     &self,
    //     node_id: NodeIndex,
    //     edge_id: EdgeIndex,
    //     skip_one: bool,
    // ) -> usize {
    //     unimplemented!()
    // }

    fn add_signs_to_edges(&self, node_id: NodeIndex) -> Vec<isize> {
        let node_hairs: SuBitGraph = self.iter_crown(node_id).into();

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
        let n_initial = self.initial_state_cut.nedges(self);
        self.cyclotomatic_number(&internal_subgraph) - n_initial
    }

    fn get_real_mass_vector<T: FloatLike>(&self, model: &Model) -> EdgeVec<F<T>> {
        self.new_edgevec(|edge, _edge_id, _| {
            let c = edge
                .mass_value(model, &self.param_builder)
                .unwrap_or(Complex {
                    re: F::from_f64(0.0),
                    im: F::from_f64(0.0),
                });

            if c.im.is_zero() {
                c.re
            } else {
                panic!(
                    "Complex masses not yet supported in gammaLoop for {}:{}",
                    edge.mass_atom(),
                    c
                )
            }
        })
    }

    fn get_external_masses<T: FloatLike>(&self, model: &Model) -> TiVec<ExternalIndex, F<T>> {
        let external_filter: SuBitGraph = self.external_filter();

        self.iter_edges_of(&external_filter)
            .map(|(_, _, edge)| {
                let c = edge
                    .data
                    .mass_value(model, &self.param_builder)
                    .unwrap_or(Complex {
                        re: F::from_f64(0.0),
                        im: F::from_f64(0.0),
                    });

                if c.im.is_zero() {
                    c.re
                } else {
                    panic!(
                        "Complex masses not yet supported in gammaLoop for {}:{}",
                        edge.data.mass_atom(),
                        c
                    )
                }
            })
            .collect()
    }

    fn get_energy_cache<T: FloatLike>(
        &self,
        model: &Model,
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
                            .on_shell_energy(edge.data.mass_value(model, &self.param_builder).map(
                                |m| {
                                    if m.im.is_non_zero() {
                                        panic!("Complex masses not yet supported in gammaLoop")
                                    }
                                    F::<T>::from_ff64(m.re)
                                },
                            ))
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
        if self.initial_state_cut.nedges(&self.underlying) == 0 {
            self.iter_edges()
                .flat_map(|(pair, edge_id, _)| {
                    if let HedgePair::Paired { .. } = pair {
                        let emr_vec = lmb.edge_signatures[edge_id]
                            .compute_three_momentum_from_four(loop_moms, external_moms);
                        vec![emr_vec.px, emr_vec.py, emr_vec.pz]
                    } else {
                        vec![]
                    }
                })
                .collect()
        } else {
            let underlying_without_is_cut = self
                .underlying
                .full_filter()
                .subtract(&self.initial_state_cut.left)
                .subtract(&self.initial_state_cut.right);

            self.underlying
                .iter_edges_of(&underlying_without_is_cut)
                .sorted_by(|a, b| a.1.cmp(&b.1))
                .flat_map(|(_pair, edge_id, _)| {
                    let emr_vec = lmb.edge_signatures[edge_id]
                        .compute_three_momentum_from_four(loop_moms, external_moms);
                    vec![emr_vec.px, emr_vec.py, emr_vec.pz]
                })
                .collect()
        }
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
        let externals: SuBitGraph = self.external_filter();

        SignatureLike::from_iter(externals.included_iter().map(|h| match self.flow(h) {
            Flow::Source => SignOrZero::Minus,
            Flow::Sink => SignOrZero::Plus,
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

    fn expected_scale(&self, e_cm: F<f64>, model: &Model) -> F<f64> {
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
                                        model.couplings[coupling]
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

        for (_, _, edge_data) in self.iter_edges_of(
            &self
                .full_filter()
                .subtract(&self.external_filter::<SuBitGraph>()),
        ) {
            scale *= Complex::new_re(e_cm.powi(edge_data.data.dod))
        }

        scale.norm_squared().sqrt()
    }

    fn no_dummy(&self) -> SuBitGraph {
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
            .filter_map(|(_, edge_index, edge_data)| {
                if edge_data.data.is_dummy {
                    Some(edge_index)
                } else {
                    None
                }
            })
            .collect()
    }

    fn all_st_cuts_for_cs(
        &self,
        source_nodes: HedgeNode,
        target_nodes: HedgeNode,
        initial_state_tree: &SuBitGraph,
    ) -> Vec<(SuBitGraph, OrientedCut, SuBitGraph)> {
        self.underlying
            .all_cuts(source_nodes, target_nodes)
            .into_iter()
            .map(|(mut l, mut c, mut r)| {
                // remove initial state cut edges from cut
                c.left.subtract_with(&self.initial_state_cut.left);
                c.right.subtract_with(&self.initial_state_cut.left);
                c.left.subtract_with(&self.initial_state_cut.right);
                c.right.subtract_with(&self.initial_state_cut.right);
                l.subtract_with(initial_state_tree);
                r.subtract_with(initial_state_tree);

                (l, c, r)
            })
            .collect_vec()
    }
}
