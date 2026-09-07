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
    function,
    id::Replacement,
    transcendental::{coth, csch, sech, tanh},
};
use typed_index_collections::TiVec;

use crate::{
    cff::generation::ShiftRewrite,
    integrands::process::param_builder::{ParamBuilderGraph, SplitPolarizations},
    model::{ArcParticle, Model, UFOSymbol},
    momentum::sample::{ExternalFourMomenta, ExternalIndex, LoopMomenta},
    momentum::signature::{ExternalSignature, SignatureLike},
    momentum::{PolDef, SignOrZero},
    numerator::{graph::ReversibleEdge, ufo::UFO},
    utils::{
        F, FloatLike, GS, external_energy_atom_from_index, ose_atom_from_index,
        symbols::ThermalDistributionLimit,
    },
    uv::uv_graph::UVE,
};

use super::{
    Edge, Graph, HedgeData, LoopMomentumBasis, Vertex, get_cff_inverse_energy_product_impl,
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

impl Graph {
    /// LU uses real on-shell energies and stable-particle cuts; width metadata alone
    /// does not turn a real mass into a complex-mass-scheme propagator.
    /// Used UFO denominators must agree with the standard pole reconstructed from the mass.
    pub(crate) fn validate_real_masses(&self, model: &Model) -> eyre::Result<()> {
        // DOT mass expressions use evaluator slots, which may predate a model update.
        let mut parameters = self.param_builder.clone();
        parameters.update_model_values(model);
        for (pair, edge_id, edge) in self.iter_edges() {
            if let Some(mass) = edge.data.mass_value::<f64>(model, &parameters)
                && !mass.im.is_zero()
            {
                return Err(eyre::eyre!(
                    "Cross-section graph '{}' has unsupported complex mass '{}' = {} on edge {}. LU requires real masses; the complex-mass scheme and finite-width propagators are not supported.",
                    self.name,
                    edge.data.mass_atom(),
                    mass,
                    edge_id
                ));
            }
            if matches!(pair, HedgePair::Paired { .. })
                && let Some(particle) = edge.data.particle()
            {
                let propagator = model.get_propagator_for_particle(&particle.name);
                let momentum = function!(UFO.momentum, 1);
                let expected = momentum.pow(2) - Atom::var(particle.mass.0.0).pow(2);
                // Normalize only the model denominator, never the graph numerator.
                // DOT mass overrides replace the declared pole after this compatibility check.
                let difference = (propagator
                    .denominator
                    .replace(function!(UFO.momentum, function!(UFO.idx, 1, 1)))
                    .with(momentum)
                    - &expected)
                    .replace(UFOSymbol::zero().0)
                    .with(Atom::Zero);
                if !difference.expand().is_zero() {
                    return Err(eyre::eyre!(
                        "Cross-section graph '{}' has unsupported propagator denominator '{}' for particle '{}' on edge {}. LU reconstructs the standard real pole '{}'; custom denominators, including explicit finite-width terms, are not supported.",
                        self.name,
                        propagator.denominator,
                        particle.name,
                        edge_id,
                        expected
                    ));
                }
            }
        }
        Ok(())
    }
}

impl Deref for Graph {
    type Target = HedgeGraph<Edge, Vertex, HedgeData>;

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
    for (&'a Vec<(PolDef, Atom)>, &'a HedgeGraph<E, V, HedgeData>)
{
    fn polarizations(&self) -> Vec<Atom> {
        self.0.iter().map(|a| a.1.clone()).collect()
    }
}

impl<'a, V, E: UVE> ParamBuilderGraph for (&'a Vec<(PolDef, Atom)>, &'a HedgeGraph<E, V, HedgeData>)
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

    fn explicit_thermal_distribution_atom(
        &self,
        edge: EdgeIndex,
        derivative_order: usize,
        thermal_sign: Atom,
        limit: ThermalDistributionLimit,
    ) -> Option<Atom> {
        self.1
            .explicit_thermal_distribution_atom(edge, derivative_order, thermal_sign, limit)
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

impl<V, E: UVE> ParamBuilderGraph for HedgeGraph<E, V, HedgeData>
where
    for<'a> EdgeData<&'a E>: ReversibleEdge,
{
    fn iter_edge_ids(&self) -> impl Iterator<Item = EdgeIndex> + '_ {
        self.iter_edges().map(|(_, e, _)| e)
    }

    fn get_external_energy_atoms(&self) -> Vec<Atom> {
        let external_filter: SuBitGraph = self.external_filter();

        external_filter
            .included_iter()
            .map(|hedge| external_energy_atom_from_index(self[&hedge]))
            .collect()
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

    fn explicit_thermal_distribution_atom(
        &self,
        edge: EdgeIndex,
        derivative_order: usize,
        thermal_sign: Atom,
        limit: ThermalDistributionLimit,
    ) -> Option<Atom> {
        match limit {
            ThermalDistributionLimit::Default => {
                let chemical_potential = self[edge].chemical_potential_atom();
                let shifted_ose = match chemical_potential {
                    Some(mu) => ose_atom_from_index(edge) - GS.sign(edge) * mu,
                    None => ose_atom_from_index(edge),
                };

                let beta = Atom::var(GS.inverse_temperature);
                let arg = beta.clone() * shifted_ose / Atom::num(2);
                match (self[edge].is_fermion(), derivative_order) {
                    (true, 0) => Some((thermal_sign + tanh().call_args([arg])) / Atom::num(2)),
                    (false, 0) => Some((thermal_sign + coth().call_args([arg])) / Atom::num(2)),
                    (true, 1) => Some(-beta * sech().call_args([arg]).pow(2) / Atom::num(4)),
                    (false, 1) => Some(beta * csch().call_args([arg]).pow(2) / Atom::num(4)),
                    (true, 2) => Some(
                        -beta.pow(2) * sech().call_args([&arg]).pow(2) * tanh().call_args([arg])
                            / Atom::num(4),
                    ),
                    (false, 2) => Some(
                        beta.pow(2) * csch().call_args([&arg]).pow(2) * coth().call_args([arg])
                            / Atom::num(4),
                    ),
                    (_, _) => None,
                }
            }
            ThermalDistributionLimit::ZeroTemperature => {
                let chemical_potential = self[edge].chemical_potential_atom();
                let shifted_ose = match chemical_potential {
                    Some(mu) => ose_atom_from_index(edge) - GS.sign(edge) * mu,
                    None => ose_atom_from_index(edge),
                };
                let chemical_potential = self[edge].chemical_potential_atom();
                match (chemical_potential, derivative_order) {
                    (Some(_), 0) => {
                        Some(thermal_sign.clone() * GS.heaviside(thermal_sign * shifted_ose))
                    }
                    (None, 0) => Some((Atom::num(1) + thermal_sign) / Atom::num(2)),
                    (None, _) => Some(Atom::num(0)),
                    // TODO: distribution derivatives with chemical potential not yet
                    // supported at zero temperature
                    _ => None,
                }
            }
            ThermalDistributionLimit::Vacuum => match derivative_order {
                0 => Some((Atom::num(1) + thermal_sign) / Atom::num(2)),
                _ => Some(Atom::num(0)),
            },
        }
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
        let external_filter: SuBitGraph = self.external_filter();

        external_filter
            .included_iter()
            .flat_map(|hedge| {
                let edge_id = self[&hedge];
                vec![
                    GS.emr_mom(edge_id, Atom::from(ExpandedIndex::from_iter([1]))),
                    GS.emr_mom(edge_id, Atom::from(ExpandedIndex::from_iter([2]))),
                    GS.emr_mom(edge_id, Atom::from(ExpandedIndex::from_iter([3]))),
                ]
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
        self.external_momentum_edge_order()
            .into_iter()
            .map(external_energy_atom_from_index)
            .collect()
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

    fn explicit_thermal_distribution_atom(
        &self,
        edge: EdgeIndex,
        derivative_order: usize,
        thermal_sign: Atom,
        limit: ThermalDistributionLimit,
    ) -> Option<Atom> {
        self.underlying.explicit_thermal_distribution_atom(
            edge,
            derivative_order,
            thermal_sign,
            limit,
        )
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
        self.external_momentum_edge_order()
            .into_iter()
            .flat_map(|edge_id| {
                vec![
                    GS.emr_mom(edge_id, Atom::from(ExpandedIndex::from_iter([1]))),
                    GS.emr_mom(edge_id, Atom::from(ExpandedIndex::from_iter([2]))),
                    GS.emr_mom(edge_id, Atom::from(ExpandedIndex::from_iter([3]))),
                ]
            })
            .collect()
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
            .sorted_by_key(|(pair, _, _)| match pair {
                HedgePair::Unpaired { hedge, .. } => *hedge,
                _ => unreachable!(),
            })
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
        let external_filter: SuBitGraph = self.external_filter();

        self.iter_edges_of(&external_filter)
            .sorted_by_key(|(pair, _, _)| match pair {
                HedgePair::Unpaired { hedge, .. } => *hedge,
                _ => unreachable!(),
            })
            .map(|(pair, _, _)| match pair {
                HedgePair::Unpaired {
                    flow: Flow::Sink, ..
                } => 1i8,
                HedgePair::Unpaired {
                    flow: Flow::Source, ..
                } => -1i8,
                _ => unreachable!(),
            })
            .collect()
    }

    fn get_external_partcles(&self) -> Vec<ArcParticle> {
        let external_filter: SuBitGraph = self.external_filter();

        self.iter_edges_of(&external_filter)
            .sorted_by_key(|(pair, _, _)| match pair {
                HedgePair::Unpaired { hedge, .. } => *hedge,
                _ => unreachable!(),
            })
            .filter_map(|(_, _, data)| data.data.particle())
            .collect()
    }

    fn get_external_signature(&self) -> SignatureLike<ExternalIndex> {
        let externals: SuBitGraph = self.external_filter();

        SignatureLike::from_iter(
            self.iter_edges_of(&externals)
                .sorted_by_key(|(pair, _, _)| match pair {
                    HedgePair::Unpaired { hedge, .. } => *hedge,
                    _ => unreachable!(),
                })
                .map(|(pair, _, _)| match pair {
                    HedgePair::Unpaired {
                        flow: Flow::Source, ..
                    } => SignOrZero::Minus,
                    HedgePair::Unpaired {
                        flow: Flow::Sink, ..
                    } => SignOrZero::Plus,
                    _ => unreachable!(),
                }),
        )
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

            scale *= Complex::new_re(e_cm.powi(vertex.dod.value)) * coupling_value;
        }

        for (_, _, edge_data) in self.iter_edges_of(
            &self
                .full_filter()
                .subtract(&self.external_filter::<SuBitGraph>()),
        ) {
            scale *= Complex::new_re(e_cm.powi(edge_data.data.dod.value))
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        dot,
        graph::parse::from_dot::IntoGraph,
        initialisation::test_initialise,
        utils::{ArbPrec, load_generic_model},
    };
    use symbolica::{
        evaluate::{CompileOptions, ExportSettings, FunctionMap, OptimizationSettings},
        symbol,
    };

    #[test]
    fn cross_section_real_masses_accept_zero_real_and_width_metadata() -> eyre::Result<()> {
        test_initialise()?;
        let model = load_generic_model("sm");
        assert!(!model.get_parameter("WT").value.unwrap().re.is_zero());
        let graph: Graph = r#"
            digraph real_masses {
                node [num=1];
                edge [num=1];
                A -> B [id=0 mass=0];
                A -> B [id=1 mass=2];
                A -> B [id=2 particle="t"];
            }
        "#
        .into_graph(&model)?;
        graph.validate_real_masses(&model)
    }

    #[test]
    fn cross_section_standard_denominators_accept_wrapping_and_factorization() -> eyre::Result<()> {
        test_initialise()?;
        let mut model = load_generic_model("sm");
        for name in ["g", "t"] {
            let position = model
                .propagators
                .iter()
                .position(|propagator| propagator.particle.name == name)
                .unwrap();
            let mass = Atom::var(model.get_particle(name).mass.0.0);
            for momentum in [
                function!(UFO.momentum, 1),
                function!(UFO.momentum, function!(UFO.idx, 1, 1)),
            ] {
                // The factored form certifies denominator-only normalization.
                std::sync::Arc::make_mut(&mut model.propagators[position]).denominator =
                    (&momentum - &mass) * (&momentum + &mass);
                let graph: Graph = format!(
                    r#"digraph standard_denominator {{
                        node [num=1]; edge [num=1];
                        A -> B [particle="{name}"];
                    }}"#
                )
                .into_graph(&model)?;
                graph.validate_real_masses(&model)?;
            }
        }
        Ok(())
    }

    #[test]
    fn cross_section_custom_denominators_reject_used_lines_only() -> eyre::Result<()> {
        test_initialise()?;
        let mut model = load_generic_model("sm");
        let position = model
            .propagators
            .iter()
            .position(|propagator| propagator.particle.name == "H")
            .unwrap();
        let width_term = Atom::i()
            * Atom::var(model.get_parameter("MH").name.0)
            * Atom::var(model.get_parameter("WH").name.0);
        std::sync::Arc::make_mut(&mut model.propagators[position]).denominator += width_term;
        let unused: Graph = r#"digraph unused_width {
            node [num=1]; edge [num=1];
            A -> B [particle="g"];
        }"#
        .into_graph(&model)?;
        unused.validate_real_masses(&model)?;
        let amputated: Graph = r#"digraph amputated_width {
            ext [style=invis]; vertex [num=1]; edge [num=1];
            ext -> vertex [particle="H"];
        }"#
        .into_graph(&model)?;
        assert!(amputated.iter_edges().all(|(pair, _, _)| !pair.is_paired()));
        amputated.validate_real_masses(&model)?;
        let used: Graph = r#"digraph explicit_width {
            node [num=1]; edge [num=1];
            A -> B [particle="H"];
        }"#
        .into_graph(&model)?;
        let error = used.validate_real_masses(&model).unwrap_err().to_string();
        assert!(error.contains("explicit_width"));
        assert!(error.contains("unsupported propagator denominator"));
        assert!(error.contains("finite-width"));
        Ok(())
    }

    #[test]
    fn cross_section_denominator_validation_uses_the_current_model() -> eyre::Result<()> {
        test_initialise()?;
        let mut model = load_generic_model("sm");
        let mass = Atom::var(model.get_parameter("MT").name.0).to_canonical_string();
        let graph: Graph = format!(
            r#"digraph updated_denominator {{
                node [num=1]; edge [num=1];
                A -> B [particle="t" mass="{mass}+1"];
            }}"#
        )
        .into_graph(&model)?;
        graph.validate_real_masses(&model)?;
        let position = model
            .propagators
            .iter()
            .position(|propagator| propagator.particle.name == "t")
            .unwrap();
        // A loaded graph must not keep validating its old propagator declaration.
        // The same owner runs again at direct integrand warm-up.
        std::sync::Arc::make_mut(&mut model.propagators[position]).denominator += Atom::one();
        let error = graph.validate_real_masses(&model).unwrap_err().to_string();
        assert!(error.contains("updated_denominator"));
        assert!(error.contains("unsupported propagator denominator"));
        Ok(())
    }

    #[test]
    fn cross_section_complex_mass_is_rejected_at_generation() -> eyre::Result<()> {
        test_initialise()?;
        let model = load_generic_model("sm");
        let graph: Graph = r#"
            digraph complex_mass {
                node [num=1];
                edge [num=1];
                A -> B [id=0 mass="1+1𝑖"];
                A -> B [id=1 mass=0];
            }
        "#
        .into_graph(&model)?;
        let error = graph.validate_real_masses(&model).unwrap_err().to_string();
        assert!(error.contains("complex_mass"));
        assert!(error.contains("unsupported complex mass"));
        assert!(error.contains("edge e0"));
        assert!(error.contains("complex-mass scheme"));
        Ok(())
    }

    #[test]
    fn cross_section_mass_validation_refreshes_runtime_model_values() -> eyre::Result<()> {
        test_initialise()?;
        for offset in ["", "+1"] {
            let mut model = load_generic_model("sm");
            let mass = Atom::var(model.get_parameter("MT").name.0).to_canonical_string();
            let graph: Graph = format!(
                r#"digraph runtime_mass {{
                    node [num=1];
                    edge [num=1];
                    A -> B [id=0 mass="{mass}{offset}"];
                    A -> B [id=1 mass=0];
                }}"#
            )
            .into_graph(&model)?;
            graph.validate_real_masses(&model)?;
            let cached_values = graph.param_builder.values.clone();

            model.get_parameter_mut("MT")?.value = Some(Complex::new(F(2.0), F(0.5)));
            if !offset.is_empty() {
                // The DOT evaluator still sees the old real value until its slots refresh.
                assert!(
                    graph[EdgeIndex(0)]
                        .mass_value::<f64>(&model, &graph.param_builder)
                        .unwrap()
                        .im
                        .is_zero()
                );
            }
            let error = graph.validate_real_masses(&model).unwrap_err().to_string();
            assert!(error.contains("runtime_mass"));
            assert!(error.contains("unsupported complex mass"));
            assert_eq!(graph.param_builder.values, cached_values);

            model.get_parameter_mut("MT")?.value = Some(Complex::new_re(F(3.0)));
            graph.validate_real_masses(&model)?;
        }
        Ok(())
    }

    #[test]
    fn native_thermal_distributions() {
        test_initialise().unwrap();
        let graph: Graph = dot!(digraph thermal {
            node [num=1]
            edge [num=1]
            A -> B [particle="d"]
            B -> A [particle="d"]
            A -> B [particle="g"]
        })
        .unwrap();
        let directory =
            std::env::temp_dir().join(format!("gammaloop-native-thermal-{}", std::process::id()));
        std::fs::create_dir_all(&directory).unwrap();
        let energy = symbol!("thermal_test_energy").to_atom();
        let mu = symbol!("thermal_test_mu").to_atom();
        let orientation = symbol!("thermal_test_orientation").to_atom();
        let thermal_sign = symbol!("thermal_test_sign").to_atom();
        for edge in [EdgeIndex(0), EdgeIndex(2)] {
            let fermion = graph[edge].is_fermion();
            assert_eq!(fermion, edge == EdgeIndex(0));
            let chemical_potential = graph[edge].chemical_potential_atom();
            let has_mu = chemical_potential.as_ref().is_some_and(|mu| !mu.is_zero());
            assert_eq!(has_mu, fermion);
            for order in 0..=2 {
                let mut expression = graph
                    .explicit_thermal_distribution_atom(
                        edge,
                        order,
                        thermal_sign.clone(),
                        ThermalDistributionLimit::Default,
                    )
                    .unwrap()
                    .replace(ose_atom_from_index(edge))
                    .with(energy.clone())
                    .replace(GS.inverse_temperature)
                    .with(Atom::num(2))
                    .replace(GS.sign(edge))
                    .with(orientation.clone());
                if has_mu {
                    expression = expression
                        .replace(chemical_potential.clone().unwrap())
                        .with(mu.clone());
                }
                let evaluator = expression
                    .as_view()
                    .to_evaluation_tree(
                        &FunctionMap::new(),
                        &[
                            energy.clone(),
                            mu.clone(),
                            orientation.clone(),
                            thermal_sign.clone(),
                        ],
                    )
                    .unwrap()
                    .linearize(&OptimizationSettings::default());
                let mut double = evaluator.clone().map_coeff(&|c| F::<f64>::from(&c.re));
                let mut precise = evaluator.clone().map_coeff(&|c| F::<ArbPrec>::from(&c.re));
                let path = directory.join(format!("thermal_{}_{}", edge.0, order));
                let mut compiled = evaluator
                    .export_cpp::<f64>(path.with_extension("cpp"), "thermal", ExportSettings::new())
                    .unwrap()
                    .compile(path.with_extension("so"), CompileOptions::default())
                    .unwrap()
                    .load()
                    .unwrap();
                for x in [
                    -1000.0_f64,
                    -40.0,
                    -1.0,
                    -2.0_f64.powi(-26),
                    0.0,
                    2.0_f64.powi(-26),
                    1.0,
                    40.0,
                    1000.0,
                ] {
                    if !fermion && x == 0.0 {
                        continue; // The Bose distribution has a pole at zero.
                    }
                    let q = (-2.0 * x.abs()).exp();
                    let denominator = if fermion {
                        1.0 + q
                    } else {
                        -(-2.0 * x.abs()).exp_m1()
                    };
                    let hyperbolic = x.signum()
                        * if fermion {
                            -(-2.0 * x.abs()).exp_m1() / denominator
                        } else {
                            (1.0 + q) / denominator
                        };
                    let squared = 4.0 * q / denominator.powi(2);
                    for sign in [-1.0, 1.0] {
                        for chemical_potential in [0.0, 0.5] {
                            let energy_value = x + if has_mu {
                                sign * chemical_potential
                            } else {
                                0.0
                            };
                            let expected = match order {
                                0 => (sign + hyperbolic) / 2.0,
                                1 => {
                                    if fermion {
                                        -squared / 2.0
                                    } else {
                                        squared / 2.0
                                    }
                                }
                                2 => {
                                    if fermion {
                                        -squared * hyperbolic
                                    } else {
                                        squared * hyperbolic
                                    }
                                }
                                _ => unreachable!(),
                            };
                            let input = [energy_value, chemical_potential, sign, sign];
                            let mut output = [0.0];
                            compiled.evaluate(&input, &mut output);
                            let results = [
                                double.evaluate_single(&input.map(F)).0,
                                precise
                                    .evaluate_single(&input.map(|v| F::<ArbPrec>::from_ff64(F(v))))
                                    .into_ff64()
                                    .0,
                                output[0],
                            ];
                            for result in results {
                                let tolerance = if order == 0 {
                                    2e-15 * expected.abs().max(1.0)
                                } else {
                                    2e-13 * expected.abs() + 1e-300
                                };
                                assert!(
                                    (result - expected).abs() <= tolerance,
                                    "fermion={fermion}, order={order}, x={x}, mu={chemical_potential}, sign={sign}: {result} != {expected}"
                                );
                            }
                        }
                    }
                }
            }
        }
        std::fs::remove_dir_all(directory).unwrap();
    }
}
