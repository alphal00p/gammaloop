use std::ops::Index;

use bincode_trait_derive::{Decode, Encode};
use bitvec::vec::BitVec;
use color_eyre::Report;
use colored::Colorize;
use derive_more::{From, Into};
use eyre::eyre;
use itertools::Itertools;
use linnet::half_edge::hedgevec::EdgeVec;
use linnet::half_edge::involution::{EdgeIndex, Flow, HedgePair};
use linnet::half_edge::subgraph::SubGraphOps;
use linnet::half_edge::HedgeGraph;
use lorentz_vector::LorentzVector;
use ref_ops::RefNeg;
use serde::{Deserialize, Serialize};
use symbolica::atom::Atom;
use symbolica::domains::float::{NumericalFloatLike, Real};
use symbolica::parse;
use typed_index_collections::TiVec;

use crate::cff::cff_graph::VertexSet;
use crate::debug_info::DEBUG_LOGGER;
use crate::momentum::FourMomentum;
use crate::momentum_sample::{
    ExternalFourMomenta, ExternalIndex, ExternalThreeMomenta, LoopIndex, LoopMomenta,
};
use crate::new_cs::CrossSectionCut;
use crate::new_graph::{self, LmbIndex, LoopMomentumBasis};
use crate::signature::ExternalSignature;
use crate::utils::{
    compute_loop_part, compute_shift_part, compute_t_part_of_shift_part, cut_energy,
    external_energy_atom_from_index, ose_atom_from_index, FloatLike, F,
};

use super::generation::ShiftRewrite;
use super::surface::{self, Surface};

/// Core esurface struct
#[derive(Serialize, Deserialize, Debug, Clone, Encode, Decode)]
pub struct Esurface {
    pub energies: Vec<EdgeIndex>,
    pub external_shift: ExternalShift,
    pub vertex_set: VertexSet,
}

impl Surface for Esurface {
    fn get_positive_energies(&self) -> impl Iterator<Item = &EdgeIndex> {
        self.energies.iter()
    }

    fn get_external_shift(&self) -> impl Iterator<Item = &(EdgeIndex, i64)> {
        self.external_shift.iter()
    }
}

impl PartialEq for Esurface {
    fn eq(&self, other: &Self) -> bool {
        self.energies == other.energies && self.external_shift == other.external_shift
    }
}

impl Eq for Esurface {}

impl Esurface {
    pub fn to_atom(&self, cut_edges: &[EdgeIndex]) -> Atom {
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

    /// Compute the value of the esurface from an energy cache that can be computed from the underlying graph
    /// This is the fastest way to compute the value of all esurfaces in a full evaluation
    #[inline]
    pub fn compute_value<T: FloatLike>(&self, energy_cache: &EdgeVec<F<T>>) -> F<T> {
        surface::compute_value(self, energy_cache)
    }

    /// Only compute the shift part, useful for existence checks
    #[inline]
    pub fn compute_shift_part<T: FloatLike>(&self, energy_cache: &EdgeVec<F<T>>) -> F<T> {
        surface::compute_shift_part(self, energy_cache)
    }

    /// Compute the value of the esurface from the momenta, needed to check if an arbitrary point
    /// is inside the esurface
    #[inline]
    pub fn compute_from_momenta<T: FloatLike>(
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

        energy_sum + self.compute_shift_part_from_momenta(lmb, external_moms)
    }

    /// Only compute the shift part, useful for center finding.
    #[inline]
    pub fn compute_shift_part_from_momenta<T: FloatLike>(
        &self,
        lmb: &LoopMomentumBasis,
        external_moms: &ExternalFourMomenta<F<T>>,
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
    pub fn compute_self_and_r_derivative<T: FloatLike>(
        &self,
        radius: &F<T>,
        shifted_unit_loops: &LoopMomenta<F<T>>,
        center: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        lmb: &LoopMomentumBasis,
        real_mass_vector: &EdgeVec<F<T>>,
    ) -> (F<T>, F<T>) {
        let spatial_part_of_externals: ExternalThreeMomenta<F<T>> = external_moms
            .iter()
            .map(|mom| mom.spatial.clone())
            .collect();

        let loops: LoopMomenta<F<T>> = shifted_unit_loops
            .iter()
            .zip(center.iter())
            .map(|(loop_mom, center)| loop_mom * radius + center)
            .collect();

        let shift = self.compute_shift_part_from_momenta(lmb, external_moms);

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

    pub fn get_subgraph_components<E, V>(&self, graph: &HedgeGraph<E, V>) -> [BitVec; 2] {
        let vertex_subgraph = self.vertex_set.subgraph(graph);
        let complement = vertex_subgraph.complement(graph);
        [vertex_subgraph, complement]
    }

    #[inline]
    pub fn get_radius_guess<T: FloatLike>(
        &self,
        unit_loops: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        lmb: &LoopMomentumBasis,
    ) -> (F<T>, F<T>) {
        let const_builder = &unit_loops[LoopIndex(0)].px;

        let esurface_shift = self.compute_shift_part_from_momenta(lmb, external_moms);

        let mut radius_guess = const_builder.zero();
        let mut denominator = const_builder.zero();

        //println!("got to energy loop");
        for energy in self.energies.iter() {
            //println!("computing contribution for energy {:?}", energy);
            let signature = &lmb.edge_signatures[*energy];
            //println!("signature {:?}", signature);

            let unit_loop_part = compute_loop_part(&signature.internal, unit_loops);
            //println!("computed_loop_part {:?}", unit_loop_part);

            let shift = compute_shift_part(&signature.external, external_moms);
            //./bprintln!("computed_shift {:?}", shift);

            let three_shift = shift.spatial;
            let norm_unit_loop_part_squared = unit_loop_part.norm_squared();
            let loop_dot_shift = &unit_loop_part * three_shift;

            radius_guess += loop_dot_shift.abs() / &norm_unit_loop_part_squared;
            denominator += norm_unit_loop_part_squared.sqrt();
        }

        radius_guess += esurface_shift.abs() / denominator;
        let negative_radius = radius_guess.ref_neg();
        (radius_guess, negative_radius)
    }

    /// Write out the esurface expression in a given lmb
    pub fn string_format_in_lmb(&self, lmb: &LoopMomentumBasis) -> String {
        let mut energy_sum = self
            .energies
            .iter()
            .map(|index| {
                let signature = &lmb.edge_signatures[*index];
                format!("|{}|", signature.format_momentum())
            })
            .join(" + ");

        let shift_part = self
            .external_shift
            .iter()
            .map(|(index, sign)| {
                let signature = &lmb.edge_signatures[*index];

                let sign = if *sign == 1 {
                    "+".to_owned()
                } else if *sign == -1 {
                    "-".to_owned()
                } else {
                    format!("+{}", sign)
                };
                format!(" {} ({})^0", sign, signature.format_momentum())
            })
            .join("");

        energy_sum.push_str(&shift_part);
        energy_sum
    }

    /// Write out the esurface expression in a generic way
    pub fn string_format(&self) -> String {
        surface::string_format(self)
    }

    pub fn get_point_inside(&self) -> Vec<LorentzVector<f64>> {
        todo!()
    }

    pub fn canonicalize_shift(&mut self, shift_rewrite: &ShiftRewrite) {
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

    pub fn new_from_cut_left<E, V, H>(graph: &HedgeGraph<E, V, H>, cut: &CrossSectionCut) -> Self {
        let edges = graph
            .iter_edges_of(&cut.cut)
            .map(|(_, id, _)| id)
            .sorted()
            .collect();

        let external_shift = graph
            .iter_edges_of(&cut.left)
            .filter_map(|(hedge_pair, edge_index, _)| match hedge_pair {
                HedgePair::Unpaired { flow, .. } => match flow {
                    Flow::Sink => Some((edge_index, -1)),
                    Flow::Source => Some((edge_index, 1)),
                },
                _ => None,
            })
            .sorted_by(|a, b| a.0.cmp(&b.0))
            .collect();

        let vertex_set = graph
            .iter_nodes_of(&cut.left)
            .map(|(node_id, _, _)| VertexSet::from_usize(node_id.into()))
            .reduce(|acc, v| acc.join(&v))
            .unwrap();

        Self {
            energies: edges,
            external_shift,
            vertex_set,
        }
    }
}

pub type EsurfaceCollection = TiVec<EsurfaceID, Esurface>;

pub fn compute_esurface_cache<T: FloatLike>(
    esurfaces: &EsurfaceCollection,
    energy_cache: &EdgeVec<F<T>>,
) -> EsurfaceCache<F<T>> {
    esurfaces
        .iter()
        .map(|esurface| esurface.compute_value(energy_cache))
        .collect::<Vec<F<T>>>()
        .into()
}

pub type EsurfaceCache<T> = TiVec<EsurfaceID, T>;

/// Index type for esurface, location of an esurface in the list of all esurfaces of a graph
#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq, From, Into, Eq, Encode, Decode)]
pub struct EsurfaceID(usize);

/// Container for esurfaces that exist at a given point in the phase space
pub type ExistingEsurfaces = TiVec<ExistingEsurfaceId, EsurfaceID>;

/// Index in the list of all existing esurfaces, essentially a pointer to a pointer to an esurface
#[derive(
    Debug, From, Into, Copy, Clone, Serialize, Deserialize, PartialEq, Eq, Hash, PartialOrd, Ord,
)]
pub struct ExistingEsurfaceId(usize);

const MAX_EXPECTED_CAPACITY: usize = 32; // Used to prevent reallocations during existence check
const SHIFT_THRESHOLD: F<f64> = F(1.0e-13);
const EXISTENCE_THRESHOLD: F<f64> = F(1.0e-7);

/// Returns the list of esurfaces which may exist, must be called each time at evaluation if externals are not fixed.
#[inline]
pub fn get_existing_esurfaces<T: FloatLike>(
    esurfaces: &EsurfaceCollection,
    esurface_derived_data: &EsurfaceDerivedData,
    externals: &ExternalFourMomenta<F<T>>,
    lmb: &LoopMomentumBasis,
    debug: usize,
    e_cm: F<f64>,
) -> ExistingEsurfaces {
    if lmb.loop_edges.is_empty() {
        return ExistingEsurfaces::new();
    }
    if debug > 1 {
        println!(
            "{}",
            "Determining all esurfaces which can satisfy the existence condition".green()
        )
    }

    let mut existing_esurfaces = ExistingEsurfaces::with_capacity(MAX_EXPECTED_CAPACITY);

    for orientation_pair in &esurface_derived_data.orientation_pairs {
        if let Some((esurface_to_check_id, shift_zero_sq)) = {
            let (esurface_id, other_esurface_id) = orientation_pair;

            let esurface = &esurfaces[*esurface_id];

            if debug > 1 {
                DEBUG_LOGGER.write("esurface_pair", &(esurface_id, other_esurface_id));
            }

            let shift_part = esurface.compute_shift_part_from_momenta(lmb, externals);
            let shift_zero_sq = &shift_part * &shift_part;

            if shift_part < -F::from_ff64(SHIFT_THRESHOLD * e_cm) {
                if debug > 1 {
                    DEBUG_LOGGER.write("negative_shift_esurface", &(esurface_id, shift_part));
                }

                Some((*esurface_id, shift_zero_sq))
            } else if shift_part > F::from_ff64(SHIFT_THRESHOLD * e_cm) {
                if debug > 1 {
                    DEBUG_LOGGER.write("negative_shift_esurface", &(other_esurface_id, shift_part));
                }

                Some((*other_esurface_id, shift_zero_sq))
            } else {
                None
            }
        } {
            if let Some(esurface_derived_data) = &esurface_derived_data[esurface_to_check_id] {
                let shift_signature = &esurface_derived_data.shift_signature;

                let esurface_shift = compute_shift_part(shift_signature, externals);

                let shift_spatial_sq = esurface_shift.spatial.norm_squared();

                let mass_sum_squared = esurface_derived_data.mass_sum_squared;

                let existence_condition =
                    &shift_zero_sq - &shift_spatial_sq - F::from_ff64(mass_sum_squared);

                if debug > 1 {
                    let helper_struct = ExistenceCheckDebug {
                        esurface_id: esurface_to_check_id,
                        shift_zero_sq: shift_zero_sq.into_ff64(),
                        shift_spatial_sq: shift_spatial_sq.into_ff64(),
                        mass_sum_sq: mass_sum_squared.into_ff64(),
                        existence_condition: existence_condition.into_ff64(),
                        threshold: F::from_ff64(
                            EXISTENCE_THRESHOLD * EXISTENCE_THRESHOLD * e_cm * e_cm,
                        ),
                    };

                    DEBUG_LOGGER.write("existence_check", &helper_struct);
                }

                if existence_condition
                    > F::from_ff64(EXISTENCE_THRESHOLD * EXISTENCE_THRESHOLD * e_cm * e_cm)
                {
                    existing_esurfaces.push(esurface_to_check_id);
                }
            }
        }
    }
    existing_esurfaces
}

#[derive(Serialize)]
struct ExistenceCheckDebug {
    esurface_id: EsurfaceID,
    shift_zero_sq: F<f64>,
    shift_spatial_sq: F<f64>,
    mass_sum_sq: F<f64>,
    existence_condition: F<f64>,
    threshold: F<f64>,
}

#[derive(Clone, Serialize, Deserialize, Debug, Encode, bincode_trait_derive::Decode)]
pub struct EsurfaceDerivedData {
    esurface_data: Vec<Option<EsurfaceData>>,
    orientation_pairs: Vec<(EsurfaceID, EsurfaceID)>,
}

impl Index<EsurfaceID> for EsurfaceDerivedData {
    type Output = Option<EsurfaceData>;

    fn index(&self, index: EsurfaceID) -> &Self::Output {
        &self.esurface_data[index.0]
    }
}

#[derive(Clone, Serialize, Deserialize, Debug, Encode, bincode_trait_derive::Decode)]
pub struct EsurfaceData {
    pub cut_momentum_basis: LmbIndex,
    pub mass_sum_squared: F<f64>,
    pub shift_signature: ExternalSignature,
}

impl EsurfaceData {
    #[allow(dead_code)]
    fn existence_condition<T: FloatLike>(&self, externals: &[FourMomentum<F<T>>]) -> (F<T>, F<T>) {
        let shift = self.shift_signature.apply(externals);

        let shift_squared = shift.square();

        (
            shift.temporal.value,
            shift_squared - F::from_ff64(self.mass_sum_squared),
        )
    }

    pub fn compute_shift_part_from_externals<T: FloatLike>(
        &self,
        externals: &[FourMomentum<F<T>>],
    ) -> F<T> {
        self.shift_signature
            .apply_iter::<_, F<T>>(externals.iter().map(|mom| &mom.temporal.value))
            .unwrap()
    }
}

pub fn generate_esurface_data(
    graph: &new_graph::Graph,
    lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
    esurfaces: &EsurfaceCollection,
) -> Result<EsurfaceDerivedData, Report> {
    let edge_masses = graph
        .underlying
        .new_edgevec(|edge, _, _| edge.particle.0.mass.value);

    let data = esurfaces
        .iter()
        .map(|esurface| {
            // find the cut momentum basis
            let energies = &esurface.energies;

            if let Some((lmb_index, lmb)) = lmbs.iter_enumerated().find(|(_i, lmb)| {
                lmb.loop_edges
                    .iter()
                    .filter(|&i| energies.contains(i))
                    .count()
                    == energies.len() - 1
            }) {
                let energy_not_in_cmb = *energies
                    .iter()
                    .find(|&i| !lmb.loop_edges.contains(i))
                    .ok_or_else(|| eyre!("No remaining edge in esurface"))?;

                let shift_signature = lmb.edge_signatures[energy_not_in_cmb].external.clone();

                let mass_sum: F<f64> = esurface
                    .energies
                    .iter()
                    .map(|&i| edge_masses[i])
                    .filter(|mass| mass.is_some())
                    .map(|mass| mass.unwrap_or_else(|| unreachable!()).re)
                    .reduce(|acc, x| acc + x)
                    .unwrap_or_else(|| F::from_f64(0.0));

                Ok(Some(EsurfaceData {
                    cut_momentum_basis: lmb_index,
                    mass_sum_squared: mass_sum * mass_sum,
                    shift_signature,
                }))
            } else {
                if esurface.external_shift.is_empty() {
                    Ok(None)
                } else {
                    Err(eyre!(
                        "Could not find cut momentum basis for esurface: {:?}",
                        esurface
                    ))
                }
            }
        })
        .collect::<Result<Vec<_>, Report>>()?;

    let mut esurface_ids = (0..esurfaces.len())
        .map(Into::<EsurfaceID>::into)
        .collect_vec();
    let mut orientation_pairs = Vec::with_capacity(esurface_ids.len() / 2);

    while let Some(esurface_id) = esurface_ids.pop() {
        let esurface = &esurfaces[esurface_id];

        if esurface.external_shift.is_empty() {
            // these do not always have a pair
            continue;
        }

        if let Some((position, other_esurface_id)) =
            esurface_ids
                .iter()
                .enumerate()
                .find(|(_pos, other_esurface_id)| {
                    let other_esurface = &esurfaces[**other_esurface_id];

                    let energies_match = esurface.energies == other_esurface.energies;

                    let sign_flipped = esurface
                        .external_shift
                        .iter()
                        .zip(other_esurface.external_shift.iter())
                        .all(|((index, signature), (other_index, other_signature))| {
                            index == other_index && *signature == -other_signature
                        })
                        && esurface.external_shift.len() == other_esurface.external_shift.len();

                    energies_match && sign_flipped
                })
        {
            orientation_pairs.push((esurface_id, *other_esurface_id));
            esurface_ids.remove(position);
        }
    }

    Ok(EsurfaceDerivedData {
        esurface_data: data,
        orientation_pairs,
    })
}

pub type ExternalShift = Vec<(EdgeIndex, i64)>;

/// add two external shifts, eliminates zero signs and sorts
pub fn add_external_shifts(lhs: &ExternalShift, rhs: &ExternalShift) -> ExternalShift {
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
    use itertools::Itertools;
    use linnet::half_edge::builder::HedgeGraphBuilder;
    use linnet::half_edge::involution::{EdgeIndex, Flow, Orientation};
    use linnet::half_edge::nodestore::NodeStorageVec;
    use linnet::half_edge::HedgeGraph;
    use symbolica::atom::{Atom, AtomCore};
    use symbolica::parse;

    use crate::cff::cff_graph::VertexSet;
    use crate::new_cs::CrossSectionCut;
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
        };

        let res = esurface.compute_value(&energies_cache);
        assert_eq!(res.0, 15.);

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
        };

        let res = esurface.compute_value(&energies_cache);
        assert_eq!(res.0, 2.);
    }

    #[test]
    fn test_to_atom() {
        let external_shift = vec![(EdgeIndex::from(1), -1)];

        let esurface = Esurface {
            energies: vec![EdgeIndex::from(2), EdgeIndex::from(3)],
            external_shift,
            vertex_set: VertexSet::dummy(),
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
        };

        let esurface_2 = Esurface {
            energies: vec![EdgeIndex::from(3), EdgeIndex::from(5)],
            external_shift: vec![(EdgeIndex::from(0), 1), (EdgeIndex::from(1), 1)],
            vertex_set: VertexSet::dummy(),
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
            .map(|cut| Esurface::new_from_cut_left(&double_triangle, &cut))
            .collect_vec();

        let expected_esurfaces = vec![
            Esurface {
                energies: vec![EdgeIndex::from(0), EdgeIndex::from(1)],
                external_shift: vec![(EdgeIndex::from(5), -1)],
                vertex_set: VertexSet::dummy(),
            },
            Esurface {
                energies: vec![EdgeIndex::from(0), EdgeIndex::from(2), EdgeIndex::from(4)],
                external_shift: vec![(EdgeIndex::from(5), -1)],
                vertex_set: VertexSet::dummy(),
            },
            Esurface {
                energies: vec![EdgeIndex::from(3), EdgeIndex::from(4)],
                external_shift: vec![(EdgeIndex::from(5), -1)],
                vertex_set: VertexSet::dummy(),
            },
            Esurface {
                energies: vec![EdgeIndex::from(1), EdgeIndex::from(2), EdgeIndex::from(3)],
                external_shift: vec![(EdgeIndex::from(5), -1)],
                vertex_set: VertexSet::dummy(),
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
            .map(|cut| Esurface::new_from_cut_left(&box_graph, &cut))
            .collect_vec();

        let expected_esurfaces = vec![
            Esurface {
                energies: vec![EdgeIndex::from(0), EdgeIndex::from(3)],
                external_shift: vec![(EdgeIndex::from(4), -1)],
                vertex_set: VertexSet::dummy(),
            },
            Esurface {
                energies: vec![EdgeIndex::from(0), EdgeIndex::from(2)],
                external_shift: vec![(EdgeIndex::from(4), -1), (EdgeIndex::from(7), -1)],
                vertex_set: VertexSet::dummy(),
            },
            Esurface {
                energies: vec![EdgeIndex::from(1), EdgeIndex::from(3)],
                external_shift: vec![(EdgeIndex::from(4), -1), (EdgeIndex::from(5), 1)],
                vertex_set: VertexSet::dummy(),
            },
            Esurface {
                energies: vec![EdgeIndex::from(1), EdgeIndex::from(2)],
                external_shift: vec![
                    (EdgeIndex::from(4), -1),
                    (EdgeIndex::from(5), 1),
                    (EdgeIndex::from(7), -1),
                ],
                vertex_set: VertexSet::dummy(),
            },
        ];

        for expected_esurface in expected_esurfaces {
            assert!(cross_section_cuts.contains(&expected_esurface));
        }
    }
}
