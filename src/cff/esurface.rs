use std::fmt::Display;

use bincode_trait_derive::{Decode, Encode};
use bitvec::vec::BitVec;
use derive_more::{From, Into};
use eyre::eyre;
use itertools::Itertools;
use linnet::half_edge::involution::{EdgeIndex, EdgeVec, Flow, HedgePair};
use linnet::half_edge::subgraph::{InternalSubGraph, ModifySubgraph, OrientedCut, SubGraphOps};
use linnet::half_edge::HedgeGraph;
use lorentz_vector::LorentzVector;
use ref_ops::RefNeg;
use serde::{Deserialize, Serialize};
use spenso::network::graph;
use symbolica::atom::{Atom, AtomCore};
use symbolica::domains::float::{NumericalFloatLike, Real};
use symbolica::id::Replacement;
use symbolica::{function, parse};
use typed_index_collections::TiVec;

use crate::cff::cff_graph::VertexSet;

use crate::define_index;
use crate::gammaloop_integrand::GammaloopIntegrand;
use crate::graph::{Graph, GraphGroupPosition, LmbIndex, LoopMomentumBasis};

use crate::momentum::ThreeMomentum;
use crate::momentum_sample::{
    ExternalFourMomenta, ExternalIndex, ExternalThreeMomenta, LoopIndex, LoopMomenta, Subspace,
    SubspaceData,
};
use crate::processes::CrossSectionCut;
use crate::utils::{
    compute_loop_part, compute_loop_part_subspace, compute_shift_part, compute_shift_part_subspace,
    compute_t_part_of_shift_part, cut_energy, external_energy_atom_from_index, ose_atom_from_index,
    FloatLike, F, GS,
};
use crate::uv::uv_graph::UVE;
use color_eyre::Result;

use super::generation::ShiftRewrite;

/// Core esurface struct
#[derive(Serialize, Deserialize, Debug, Clone, bincode::Encode, bincode::Decode)]
pub struct Esurface {
    pub energies: Vec<EdgeIndex>,
    pub external_shift: ExternalShift,
    pub vertex_set: VertexSet,
    //#[bincode(with_serde)]
    //pub subspace_graph: InternalSubGraph,
}

impl PartialEq for Esurface {
    fn eq(&self, other: &Self) -> bool {
        self.energies == other.energies && self.external_shift == other.external_shift
    }
}

impl Eq for Esurface {}

impl Esurface {
    pub(crate) fn to_atom(&self, cut_edges: &[EdgeIndex]) -> Atom {
        let symbolic_energies = self
            .energies
            .iter()
            .map(|i| {
                if cut_edges.contains(i) {
                    cut_energy(*i)
                } else {
                    ose_atom_from_index(*i)
                }
            })
            .collect_vec();

        let symbolic_shift = self
            .external_shift
            .iter()
            .fold(Atom::new(), |sum, (i, sign)| {
                external_energy_atom_from_index(*i) * &Atom::num(*sign) + &sum
            });

        let builder_atom = Atom::new();
        let energy_sum = symbolic_energies
            .iter()
            .fold(builder_atom, |acc, energy| acc + energy);

        energy_sum + &symbolic_shift
    }

    /// Compute the value of the esurface from the momenta, needed to check if an arbitrary point
    /// is inside the esurface
    #[inline]
    pub(crate) fn compute_from_momenta<T: FloatLike>(
        &self,
        lmb: &LoopMomentumBasis,
        real_mass_vector: &EdgeVec<F<T>>,
        loop_moms: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
    ) -> F<T> {
        let spatial_part_of_externals = external_moms
            .iter()
            .map(|mom| mom.spatial.clone())
            .collect::<TiVec<ExternalIndex, _>>();

        let energy_sum = self
            .energies
            .iter()
            .map(|index| {
                let signature = &lmb.edge_signatures[*index];
                let momentum = signature.compute_momentum(loop_moms, &spatial_part_of_externals);
                let mass = &real_mass_vector[*index];

                (momentum.norm_squared() + mass * mass).sqrt()
            })
            .reduce(|acc, x| acc + x)
            .unwrap_or_else(|| loop_moms[LoopIndex(0)].px.zero());

        let shift_part = self
            .external_shift
            .iter()
            .map(|(index, sign)| {
                let external_signature = &lmb.edge_signatures[*index].external;
                F::from_f64(*sign as f64)
                    * compute_t_part_of_shift_part(external_signature, external_moms)
            })
            .reduce(|acc, x| acc + x)
            .unwrap_or_else(|| energy_sum.zero());

        energy_sum + shift_part
    }

    #[inline]
    /// TODO: upgrade to a status type
    pub(crate) fn exists_subspace<T: FloatLike>(
        &self,
        loop_moms: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        subspace: &SubspaceData,
        all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
        graph: &Graph,
        real_mass_vector: &EdgeVec<F<T>>,
        e_cm: &F<T>,
    ) -> bool {
        //todo!("refactor for subspaces");
        if self.external_shift.is_empty() {
            return false;
        }

        let shift_part = self.compute_shift_part_from_momenta_in_subspace(
            loop_moms,
            external_moms,
            subspace,
            all_lmbs,
            graph,
            real_mass_vector,
        );

        if shift_part < -F::from_ff64(SHIFT_THRESHOLD) * e_cm {
            let lmb = subspace.get_lmb(all_lmbs);

            let mass_sum: F<T> = self
                .energies
                .iter()
                .map(|index| &real_mass_vector[*index])
                .fold(F::from_f64(0.0), |acc, x| acc + x);

            let zero_vector = ThreeMomentum::new(e_cm.zero(), e_cm.zero(), e_cm.zero());

            let graph_vector = self
                .external_shift
                .iter()
                .map(|(index, sign)| {
                    let external_signature = &lmb.edge_signatures[*index].external;
                    compute_shift_part(external_signature, external_moms).spatial
                        * F::from_f64(*sign as f64)
                })
                .reduce(|acc, x| acc + x)
                .unwrap_or_else(|| zero_vector.clone());

            let other_part = subspace
                .does_not_contain(&self.energies, graph)
                .into_iter()
                .map(|index| {
                    let signature = &lmb.edge_signatures[index];
                    let momentum = signature.compute_momentum(
                        loop_moms,
                        &external_moms
                            .iter()
                            .map(|mom| mom.spatial.clone())
                            .collect::<TiVec<ExternalIndex, _>>(),
                    );
                    momentum
                })
                .reduce(|acc, x| acc + x)
                .unwrap_or_else(|| zero_vector.clone());

            let shift_vector_sq = (&graph_vector + &other_part).norm_squared();

            &shift_part * &shift_part - shift_vector_sq - &mass_sum * &mass_sum
                > F::from_ff64(EXISTENCE_THRESHOLD) * e_cm * e_cm
        } else {
            false
        }
    }

    #[inline]
    /// TODO: upgrade to a status type
    pub(crate) fn exists<T: FloatLike>(
        &self,
        external_moms: &ExternalFourMomenta<F<T>>,
        lmb: &LoopMomentumBasis,
        real_mass_vector: &EdgeVec<F<T>>,
        e_cm: &F<T>,
    ) -> bool {
        //todo!("refactor for subspaces");
        if self.external_shift.is_empty() {
            return false;
        }

        let shift_part = self.compute_shift_part_from_momenta(external_moms, lmb);

        if shift_part < -F::from_ff64(SHIFT_THRESHOLD) * e_cm {
            let mass_sum: F<T> = self
                .energies
                .iter()
                .map(|index| &real_mass_vector[*index])
                .fold(F::from_f64(0.0), |acc, x| acc + x);

            let zero_vector = ThreeMomentum::new(e_cm.zero(), e_cm.zero(), e_cm.zero());

            let shift_vector = self
                .external_shift
                .iter()
                .map(|(index, sign)| {
                    let external_signature = &lmb.edge_signatures[*index].external;
                    compute_shift_part(external_signature, external_moms).spatial
                        * F::from_f64(*sign as f64)
                })
                .reduce(|acc, x| acc + x)
                .unwrap_or_else(|| zero_vector.clone());

            let shift_vector_sq = shift_vector.norm_squared();

            &shift_part * &shift_part - shift_vector_sq - &mass_sum * &mass_sum
                > F::from_ff64(EXISTENCE_THRESHOLD) * e_cm * e_cm
        } else {
            false
        }
    }

    /// Only compute the shift part, useful for center finding.
    pub(crate) fn compute_shift_part_from_momenta_in_subspace<T: FloatLike>(
        &self,
        loop_moms: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        subspace: &SubspaceData,
        all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
        graph: &Graph,
        masses: &EdgeVec<F<T>>,
    ) -> F<T> {
        let lmb = subspace.get_lmb(all_lmbs);

        let full_external_shift = self
            .external_shift
            .iter()
            .map(|(index, sign)| {
                let external_signature = &lmb.edge_signatures[*index].external;
                F::from_f64(*sign as f64)
                    * compute_t_part_of_shift_part(external_signature, external_moms)
            })
            .reduce(|acc, x| acc + x)
            .unwrap_or_else(|| external_moms[ExternalIndex(0)].temporal.value.zero());

        let spatial_externals = external_moms
            .iter()
            .map(|mom| mom.spatial.clone())
            .collect::<TiVec<ExternalIndex, _>>();

        let remaining_shift = subspace
            .does_not_contain(&self.energies, graph)
            .into_iter()
            .map(|index| {
                let signature = &lmb.edge_signatures[index];
                let momentum = signature.compute_momentum(loop_moms, &spatial_externals);
                let mass = &masses[index];
                (momentum.norm_squared() + mass * mass).sqrt()
            })
            .reduce(|acc, x| acc + x)
            .unwrap_or_else(|| full_external_shift.zero());

        full_external_shift + remaining_shift
    }

    /// Only compute the shift part, useful for center finding.
    pub(crate) fn compute_shift_part_from_momenta<T: FloatLike>(
        &self,
        external_moms: &ExternalFourMomenta<F<T>>,
        lmb: &LoopMomentumBasis,
    ) -> F<T> {
        self.external_shift
            .iter()
            .map(|(index, sign)| {
                let external_signature = &lmb.edge_signatures[*index].external;
                F::from_f64(*sign as f64)
                    * compute_t_part_of_shift_part(external_signature, external_moms)
            })
            .reduce(|acc, x| acc + x)
            .unwrap_or_else(|| external_moms[ExternalIndex(0)].temporal.value.zero())
    }

    #[inline]
    pub(crate) fn compute_self_and_r_derivative_subspace<T: FloatLike>(
        &self,
        radius: &F<T>,
        shifted_unit_loops_in_subspace: &LoopMomenta<F<T>>,
        center_in_subspace: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        real_mass_vector: &EdgeVec<F<T>>,
        subspace: &SubspaceData,
        all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
        graph: &Graph,
    ) -> (F<T>, F<T>) {
        let spatial_part_of_externals: ExternalThreeMomenta<F<T>> = external_moms
            .iter()
            .map(|mom| mom.spatial.clone())
            .collect();

        let loops: LoopMomenta<F<T>> = subspace
            .iter_lmb_indices()
            .map(|loop_index| {
                &shifted_unit_loops_in_subspace[loop_index] * radius
                    + &center_in_subspace[loop_index]
            })
            .collect();

        let shift = self.compute_shift_part_from_momenta_in_subspace(
            shifted_unit_loops_in_subspace,
            external_moms,
            subspace,
            all_lmbs,
            graph,
            real_mass_vector,
        );

        let lmb = subspace.get_lmb(all_lmbs);
        let (derivative, energy_sum) = subspace
            .contains(&self.energies, graph)
            .map(|index| {
                let signature = &lmb.edge_signatures[index];

                let momentum = signature.compute_momentum(&loops, &spatial_part_of_externals);
                let unit_loop_part = compute_loop_part_subspace(
                    &signature.internal,
                    shifted_unit_loops_in_subspace,
                    subspace,
                );

                let energy = (momentum.norm_squared()
                    + &real_mass_vector[index] * &real_mass_vector[index])
                    .sqrt();

                let numerator = momentum * &unit_loop_part;

                (numerator / &energy, energy)
            })
            .fold(
                (radius.zero(), radius.zero()),
                |(der_sum, en_sum), (der, en)| (der_sum + der, en_sum + en),
            );

        (energy_sum + shift, derivative)
    }

    #[inline]
    pub(crate) fn compute_self_and_r_derivative<T: FloatLike>(
        &self,
        radius: &F<T>,
        shifted_unit_loops: &LoopMomenta<F<T>>,
        center: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        real_mass_vector: &EdgeVec<F<T>>,
        lmb: &LoopMomentumBasis,
    ) -> (F<T>, F<T>) {
        let spatial_part_of_externals: ExternalThreeMomenta<F<T>> = external_moms
            .iter()
            .map(|mom| mom.spatial.clone())
            .collect();

        let loops: LoopMomenta<F<T>> = shifted_unit_loops
            .iter()
            .zip(center.iter())
            .map(|(momentum, center)| momentum * radius + center)
            .collect();

        let shift = self.compute_shift_part_from_momenta(external_moms, lmb);

        let (derivative, energy_sum) = self
            .energies
            .iter()
            .map(|&index| {
                let signature = &lmb.edge_signatures[index];

                let momentum = signature.compute_momentum(&loops, &spatial_part_of_externals);
                let unit_loop_part = compute_loop_part(&signature.internal, shifted_unit_loops);

                let energy = (momentum.norm_squared()
                    + &real_mass_vector[index] * &real_mass_vector[index])
                    .sqrt();

                let numerator = momentum * &unit_loop_part;

                (numerator / &energy, energy)
            })
            .fold(
                (radius.zero(), radius.zero()),
                |(der_sum, en_sum), (der, en)| (der_sum + der, en_sum + en),
            );

        (energy_sum + shift, derivative)
    }

    pub(crate) fn get_subgraph_components<E, V, H>(
        &self,
        graph: &HedgeGraph<E, V, H>,
    ) -> (BitVec, BitVec) {
        let vertex_subgraph = self.vertex_set.subgraph(graph);
        let complement = vertex_subgraph.complement(graph);
        (vertex_subgraph, complement)
    }

    pub(crate) fn bitvec<E, V, H>(&self, graph: &HedgeGraph<E, V, H>) -> BitVec {
        let mut result: BitVec = graph.empty_subgraph();
        for edge_id in self.energies.iter() {
            let (_, pair) = graph[edge_id];
            result.add(pair);
        }

        result
    }

    // #[inline]
    /// the "loops_unit_in_subspace" means that the loop momenta that are part of the subspace are jointly normalized to unit length
    pub(crate) fn get_radius_guess_subspace<T: FloatLike>(
        &self,
        loops_unit_in_subspace: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        subspace: &SubspaceData,
        all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
        graph: &Graph,
        masses: &EdgeVec<F<T>>,
    ) -> (F<T>, F<T>) {
        let const_builder = &loops_unit_in_subspace[LoopIndex(0)].px;

        let esurface_shift = self.compute_shift_part_from_momenta_in_subspace(
            loops_unit_in_subspace,
            external_moms,
            subspace,
            all_lmbs,
            graph,
            masses,
        );

        let mut radius_guess = const_builder.zero();
        let mut denominator = const_builder.zero();

        let lmb = subspace.get_lmb(all_lmbs);

        let external_3_momenta = external_moms.iter().map(|x| x.spatial.clone()).collect();

        //println!("got to energy loop");
        for energy in subspace.contains(&self.energies, graph).into_iter() {
            //println!("computing contribution for energy {:?}", energy);
            let signature = &lmb.edge_signatures[energy];
            //println!("signature {:?}", signature);

            let unit_loop_part =
                compute_loop_part_subspace(&signature.internal, loops_unit_in_subspace, &subspace);
            //println!("computed_loop_part {:?}", unit_loop_part);

            let three_shift = compute_shift_part_subspace(
                &signature.internal,
                &signature.external,
                &loops_unit_in_subspace,
                &external_3_momenta,
                subspace,
            );
            //./bprintln!("computed_shift {:?}", shift);

            let norm_unit_loop_part_squared = unit_loop_part.norm_squared();
            let loop_dot_shift = &unit_loop_part * three_shift;

            radius_guess += loop_dot_shift.abs() / &norm_unit_loop_part_squared;
            denominator += norm_unit_loop_part_squared.sqrt();
        }

        radius_guess += esurface_shift.abs() / denominator;
        let negative_radius = radius_guess.ref_neg();
        (radius_guess, negative_radius)
    }

    pub(crate) fn get_radius_guess<T: FloatLike>(
        &self,
        unit_loops: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        lmb: &LoopMomentumBasis,
    ) -> (F<T>, F<T>) {
        let const_builder = &unit_loops[LoopIndex(0)].px;

        let esurface_shift = self.compute_shift_part_from_momenta(external_moms, lmb);

        let mut radius_guess = const_builder.zero();
        let mut denominator = const_builder.zero();

        //println!("got to energy loop");
        for &energy in &self.energies {
            //println!("computing contribution for energy {:?}", energy);
            let signature = &lmb.edge_signatures[energy];
            //println!("signature {:?}", signature);

            let unit_loop_part = compute_loop_part(&signature.internal, unit_loops);
            //println!("computed_loop_part {:?}", unit_loop_part);

            let three_shift = compute_shift_part(&signature.external, external_moms).spatial;
            //./bprintln!("computed_shift {:?}", shift);

            let norm_unit_loop_part_squared = unit_loop_part.norm_squared();
            let loop_dot_shift = &unit_loop_part * three_shift;

            radius_guess += loop_dot_shift.abs() / &norm_unit_loop_part_squared;
            denominator += norm_unit_loop_part_squared.sqrt();
        }

        radius_guess += esurface_shift.abs() / denominator;
        let negative_radius = radius_guess.ref_neg();
        (radius_guess, negative_radius)
    }

    pub(crate) fn canonicalize_shift(&mut self, shift_rewrite: &ShiftRewrite) {
        if let Some(dep_mom_pos) = self
            .external_shift
            .iter()
            .position(|(index, _)| *index == shift_rewrite.dependent_momentum)
        {
            let (_, dep_mom_sign) = self.external_shift.remove(dep_mom_pos);

            let external_shift = shift_rewrite
                .dependent_momentum_expr
                .iter()
                .map(|(index, sign)| (*index, dep_mom_sign * sign))
                .collect();

            self.external_shift = add_external_shifts(&self.external_shift, &external_shift);
        }
    }

    pub(crate) fn new_from_cut_left<E, V, H>(
        graph: &HedgeGraph<E, V, H>,
        cut: &CrossSectionCut,
        initial_state_cut: Option<&OrientedCut>,
    ) -> Self {
        let edges = graph
            .iter_edges_of(&cut.cut)
            .map(|(_, id, _)| id)
            .sorted()
            .collect();

        let external_shift = if let Some(is_cut) = initial_state_cut {
            graph
                .iter_edges_of(is_cut)
                .map(|(_, edge_index, __)| (edge_index, -1))
                .sorted_by(|a, b| a.0.cmp(&b.0))
                .collect()
        } else {
            graph
                .iter_edges_of(&cut.left)
                .filter_map(|(hedge_pair, edge_index, _)| match hedge_pair {
                    HedgePair::Unpaired { flow, .. } => match flow {
                        Flow::Sink => Some((edge_index, -1)),
                        Flow::Source => Some((edge_index, 1)),
                    },
                    _ => None,
                })
                .sorted_by(|a, b| a.0.cmp(&b.0))
                .collect()
        };

        let vertex_set = graph
            .iter_nodes_of(&cut.left)
            .map(|(node_id, _, _)| VertexSet::from_usize(node_id.into()))
            .reduce(|acc, v| acc.join(&v))
            .unwrap();

        Self {
            energies: edges,
            external_shift,
            vertex_set,
            //subspace_graph: graph.full_graph(),
        }
    }

    pub(crate) fn lmb_atom(&self, graph: &Graph, lmb_reps: &[Replacement]) -> Atom {
        self.energies
            .iter()
            .map(|index| {
                let mass_symbol = graph.underlying[*index].mass_atom();
                let emr_symbols = (0..3)
                    .map(|i| function!(GS.emr_mom, usize::from(*index), i + 1))
                    .collect_vec();

                let atom = (&emr_symbols[0] * &emr_symbols[0]
                    + &emr_symbols[1] * &emr_symbols[1]
                    + &emr_symbols[2] * &emr_symbols[2]
                    + &mass_symbol * &mass_symbol)
                    .sqrt();
                atom
            })
            .chain(self.external_shift.iter().map(|(index, sign)| {
                function!(GS.emr_mom, usize::from(*index), 0) * Atom::num(*sign)
            }))
            .reduce(|sum, atom| sum + atom)
            .unwrap_or_else(|| Atom::new())
            .replace_multiple(lmb_reps)
            .replace(parse!("ZERO"))
            .with(Atom::new())
            .expand() // ensure canonical form
    }

    // more readable version for debugging, because it doesn't write out components
    pub(crate) fn lmb_atom_simplified(&self, graph: &Graph, lmb_reps: &[Replacement]) -> Atom {
        self.energies
            .iter()
            .map(|index| {
                let mass_symbol = graph.underlying[*index].mass_atom();
                let emr_symbol = function!(GS.emr_mom, usize::from(*index));
                let atom = (&emr_symbol * &emr_symbol + &mass_symbol * &mass_symbol).sqrt();
                atom
            })
            .chain(self.external_shift.iter().map(|(index, sign)| {
                function!(GS.emr_mom, usize::from(*index), 0) * Atom::num(*sign)
            }))
            .reduce(|sum, atom| sum + atom)
            .unwrap_or_else(|| Atom::new())
            .replace_multiple(lmb_reps)
            .replace(parse!("ZERO"))
            .with(Atom::new())
            .expand() // ensure canonical form
    }
}

define_index! {pub struct GroupEsurfaceId;}

pub type EsurfaceCollection = TiVec<EsurfaceID, Esurface>;

pub type EsurfaceCache<T> = TiVec<EsurfaceID, T>;

/// Index type for esurface, location of an esurface in the list of all esurfaces of a graph
#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq, From, Into, Eq, Encode, Decode)]
pub struct EsurfaceID(pub usize);

/// Container for esurfaces that exist at a given point in the phase space
pub type ExistingEsurfaces = TiVec<ExistingEsurfaceId, GroupEsurfaceId>;

pub(crate) fn get_representative(
    esurface_map: &TiVec<GraphGroupPosition, Option<EsurfaceID>>,
) -> Result<(GraphGroupPosition, EsurfaceID)> {
    for (group_pos, esurface_option) in esurface_map.iter_enumerated() {
        if let Some(esurface_id) = esurface_option {
            return Ok((group_pos, *esurface_id));
        }
    }

    Err(eyre!(
        "No representative esurface found, esurface map corrupted"
    ))
}

/// Index in the list of all existing esurfaces, essentially a pointer to a pointer to an esurface
#[derive(
    Debug,
    From,
    Into,
    Copy,
    Clone,
    Serialize,
    Deserialize,
    PartialEq,
    Eq,
    Hash,
    PartialOrd,
    Ord,
    Encode,
    Decode,
)]
pub struct ExistingEsurfaceId(usize);

impl Display for ExistingEsurfaceId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "ExistingEsurfaceID({})", self.0)
    }
}

const SHIFT_THRESHOLD: F<f64> = F(1.0e-13);
const EXISTENCE_THRESHOLD: F<f64> = F(1.0e-7);

pub type ExternalShift = Vec<(EdgeIndex, i64)>;

/// add two external shifts, eliminates zero signs and sorts
pub(crate) fn add_external_shifts(lhs: &ExternalShift, rhs: &ExternalShift) -> ExternalShift {
    let mut res = lhs.clone();

    for rhs_element in rhs.iter() {
        if let Some(lhs_element) = res
            .iter_mut()
            .find(|lhs_element| rhs_element.0 == lhs_element.0)
        {
            lhs_element.1 += rhs_element.1;
        } else {
            res.push(*rhs_element)
        }
    }

    res.retain(|(_index, sign)| *sign != 0);
    res.sort_by(|(index_1, _), (index_2, _)| index_1.cmp(index_2));
    res
}

impl From<EsurfaceID> for Atom {
    fn from(id: EsurfaceID) -> Self {
        parse!(&format!("η({})", Into::<usize>::into(id.0)))
    }
}

#[cfg(test)]
mod tests {
    use bitvec::vec::BitVec;
    use itertools::Itertools;
    use linnet::half_edge::builder::HedgeGraphBuilder;
    use linnet::half_edge::involution::{EdgeIndex, Flow, Orientation};
    use linnet::half_edge::nodestore::NodeStorageVec;
    use linnet::half_edge::subgraph::InternalSubGraph;
    use linnet::half_edge::HedgeGraph;
    use symbolica::atom::{Atom, AtomCore};
    use symbolica::parse;

    use crate::cff::cff_graph::VertexSet;
    use crate::processes::CrossSectionCut;
    use crate::{
        cff::{esurface::Esurface, generation::ShiftRewrite},
        utils::{dummy_hedge_graph, F},
    };

    use super::add_external_shifts;

    #[test]
    fn test_esurface() {
        let dummy_graph = dummy_hedge_graph(5);

        let energies_cache = dummy_graph
            .new_edgevec_from_iter([F(1.), F(2.), F(3.), F(4.), F(5.)])
            .unwrap();

        let energies = vec![EdgeIndex::from(0), EdgeIndex::from(1), EdgeIndex::from(2)];

        let external_shift = vec![(EdgeIndex::from(3), 1), (EdgeIndex::from(4), 1)];

        let mut esurface = Esurface {
            energies,
            external_shift,
            vertex_set: VertexSet::dummy(),
            //subspace_graph: dummy_graph.full_graph(),
        };

        let shift_rewrite = ShiftRewrite {
            dependent_momentum: EdgeIndex::from(4),
            dependent_momentum_expr: vec![
                (EdgeIndex::from(1), -1),
                (EdgeIndex::from(2), -1),
                (EdgeIndex::from(3), -1),
            ],
        };

        esurface.canonicalize_shift(&shift_rewrite);

        assert_eq!(
            esurface.external_shift,
            vec![(EdgeIndex::from(1), -1), (EdgeIndex::from(2), -1)]
        );

        let energies = vec![EdgeIndex::from(0), EdgeIndex::from(2)];

        let external_shift = vec![(EdgeIndex::from(1), -1)];

        let esurface = Esurface {
            energies,
            external_shift,
            vertex_set: VertexSet::dummy(),
            //subspace_graph: dummy_graph.full_graph(),
        };
    }

    #[test]
    fn test_to_atom() {
        let external_shift = vec![(EdgeIndex::from(1), -1)];

        let esurface = Esurface {
            energies: vec![EdgeIndex::from(2), EdgeIndex::from(3)],
            external_shift,
            vertex_set: VertexSet::dummy(),
            // subspace_graph: unsafe { InternalSubGraph::new_unchecked(BitVec::new()) },
        };

        let esurface_atom = esurface.to_atom(&[]);
        let expected_atom = parse!("Q(2, cind(0)) + Q(3, cind(0)) - P(1, cind(0))");

        let diff = esurface_atom - &expected_atom;
        let diff = diff.expand();
        assert_eq!(diff, Atom::new());
    }

    #[test]
    fn test_add_external_shifts() {
        let shift_1 = vec![
            (EdgeIndex::from(0), 1),
            (EdgeIndex::from(1), 1),
            (EdgeIndex::from(2), -1),
        ];
        let shift_2 = vec![(EdgeIndex::from(1), -1), (EdgeIndex::from(2), 1)];

        let add = add_external_shifts(&shift_1, &shift_2);

        assert_eq!(add, vec![(EdgeIndex::from(0), 1)]);

        let shift_3 = vec![(EdgeIndex::from(3), 1), (EdgeIndex::from(4), -1)];
        let shift_4 = vec![
            (EdgeIndex::from(0), 1),
            (EdgeIndex::from(1), 1),
            (EdgeIndex::from(2), 1),
            (EdgeIndex::from(4), 1),
        ];

        let add = add_external_shifts(&shift_3, &shift_4);

        assert_eq!(
            add,
            vec![
                (EdgeIndex::from(0), 1),
                (EdgeIndex::from(1), 1),
                (EdgeIndex::from(2), 1),
                (EdgeIndex::from(3), 1)
            ]
        );
    }

    #[test]
    fn test_esurface_equality() {
        let esurface_1 = Esurface {
            energies: vec![EdgeIndex::from(3), EdgeIndex::from(5)],
            external_shift: vec![(EdgeIndex::from(0), 1), (EdgeIndex::from(1), 1)],
            vertex_set: VertexSet::dummy(),
            //subspace_graph: unsafe { InternalSubGraph::new_unchecked(BitVec::new()) },
        };

        let esurface_2 = Esurface {
            energies: vec![EdgeIndex::from(3), EdgeIndex::from(5)],
            external_shift: vec![(EdgeIndex::from(0), 1), (EdgeIndex::from(1), 1)],
            vertex_set: VertexSet::dummy(),
            //subspace_graph: unsafe { InternalSubGraph::new_unchecked(BitVec::new()) },
        };

        assert_eq!(esurface_1, esurface_2);
    }

    #[test]
    fn test_from_cut_left_dt() {
        let mut hedge_graph_builder = HedgeGraphBuilder::new();
        let nodes = (0..4)
            .map(|_| hedge_graph_builder.add_node(()))
            .collect_vec();

        hedge_graph_builder.add_edge(nodes[0], nodes[1], (), Orientation::Undirected);
        hedge_graph_builder.add_edge(nodes[0], nodes[2], (), Orientation::Undirected);
        hedge_graph_builder.add_edge(nodes[1], nodes[2], (), Orientation::Undirected);
        hedge_graph_builder.add_edge(nodes[1], nodes[3], (), Orientation::Undirected);
        hedge_graph_builder.add_edge(nodes[2], nodes[3], (), Orientation::Undirected);

        hedge_graph_builder.add_external_edge(nodes[0], (), Orientation::Undirected, Flow::Sink);
        hedge_graph_builder.add_external_edge(nodes[3], (), Orientation::Undirected, Flow::Source);

        let double_triangle: HedgeGraph<(), (), ()> =
            hedge_graph_builder.build::<NodeStorageVec<()>>();
        let node_0 = double_triangle.iter_crown(nodes[0]).into();
        let node_3 = double_triangle.iter_crown(nodes[3]).into();

        let cuts = double_triangle.all_cuts(node_0, node_3);

        let cross_section_cuts = cuts
            .into_iter()
            .map(|(node_l, cut, node_r)| CrossSectionCut {
                cut,
                left: node_l,
                right: node_r,
            })
            .map(|cut| Esurface::new_from_cut_left(&double_triangle, &cut, None))
            .collect_vec();

        let expected_esurfaces = vec![
            Esurface {
                energies: vec![EdgeIndex::from(0), EdgeIndex::from(1)],
                external_shift: vec![(EdgeIndex::from(5), -1)],
                vertex_set: VertexSet::dummy(),
                //subspace_graph: double_triangle.full_graph(),
            },
            Esurface {
                energies: vec![EdgeIndex::from(0), EdgeIndex::from(2), EdgeIndex::from(4)],
                external_shift: vec![(EdgeIndex::from(5), -1)],
                vertex_set: VertexSet::dummy(),
                //subspace_graph: double_triangle.full_graph(),
            },
            Esurface {
                energies: vec![EdgeIndex::from(3), EdgeIndex::from(4)],
                external_shift: vec![(EdgeIndex::from(5), -1)],
                vertex_set: VertexSet::dummy(),
                //subspace_graph: double_triangle.full_graph(),
            },
            Esurface {
                energies: vec![EdgeIndex::from(1), EdgeIndex::from(2), EdgeIndex::from(3)],
                external_shift: vec![(EdgeIndex::from(5), -1)],
                vertex_set: VertexSet::dummy(),
                //subspace_graph: double_triangle.full_graph(),
            },
        ];

        for expected_esurface in expected_esurfaces {
            assert!(cross_section_cuts.contains(&expected_esurface));
        }
    }

    #[test]
    fn test_from_cut_left_box() {
        let mut hedge_graph_builder = HedgeGraphBuilder::new();
        let nodes = (0..4)
            .map(|_| hedge_graph_builder.add_node(()))
            .collect_vec();

        hedge_graph_builder.add_edge(nodes[0], nodes[1], (), Orientation::Undirected);
        hedge_graph_builder.add_edge(nodes[1], nodes[2], (), Orientation::Undirected);
        hedge_graph_builder.add_edge(nodes[2], nodes[3], (), Orientation::Undirected);
        hedge_graph_builder.add_edge(nodes[3], nodes[0], (), Orientation::Undirected);

        hedge_graph_builder.add_external_edge(nodes[0], (), Orientation::Undirected, Flow::Sink);
        hedge_graph_builder.add_external_edge(nodes[1], (), Orientation::Undirected, Flow::Source);
        hedge_graph_builder.add_external_edge(nodes[2], (), Orientation::Undirected, Flow::Source);
        hedge_graph_builder.add_external_edge(nodes[3], (), Orientation::Undirected, Flow::Sink);

        let box_graph: HedgeGraph<(), (), ()> = hedge_graph_builder.build::<NodeStorageVec<()>>();

        let node_0 = box_graph.iter_crown(nodes[0]).into();
        let node_2 = box_graph.iter_crown(nodes[2]).into();

        let cuts = box_graph.all_cuts(node_0, node_2);
        assert_eq!(cuts.len(), 4);

        let cross_section_cuts = cuts
            .into_iter()
            .map(|(node_l, cut, node_r)| CrossSectionCut {
                cut,
                left: node_l,
                right: node_r,
            })
            .map(|cut| Esurface::new_from_cut_left(&box_graph, &cut, None))
            .collect_vec();

        let expected_esurfaces = vec![
            Esurface {
                energies: vec![EdgeIndex::from(0), EdgeIndex::from(3)],
                external_shift: vec![(EdgeIndex::from(4), -1)],
                vertex_set: VertexSet::dummy(),
                //subspace_graph: box_graph.full_graph(),
            },
            Esurface {
                energies: vec![EdgeIndex::from(0), EdgeIndex::from(2)],
                external_shift: vec![(EdgeIndex::from(4), -1), (EdgeIndex::from(7), -1)],
                vertex_set: VertexSet::dummy(),
                //subspace_graph: box_graph.full_graph(),
            },
            Esurface {
                energies: vec![EdgeIndex::from(1), EdgeIndex::from(3)],
                external_shift: vec![(EdgeIndex::from(4), -1), (EdgeIndex::from(5), 1)],
                vertex_set: VertexSet::dummy(),
                //subspace_graph: box_graph.full_graph(),
            },
            Esurface {
                energies: vec![EdgeIndex::from(1), EdgeIndex::from(2)],
                external_shift: vec![
                    (EdgeIndex::from(4), -1),
                    (EdgeIndex::from(5), 1),
                    (EdgeIndex::from(7), -1),
                ],
                vertex_set: VertexSet::dummy(),
                //subspace_graph: box_graph.full_graph(),
            },
        ];

        for expected_esurface in expected_esurfaces {
            assert!(cross_section_cuts.contains(&expected_esurface));
        }
    }
}
