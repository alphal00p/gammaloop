use std::ops::Index;

use bincode::{Decode, Encode};
use color_eyre::Report;
use colored::Colorize;
use derive_more::{From, Into};
use eyre::eyre;
use itertools::Itertools;
use lorentz_vector::LorentzVector;
use ref_ops::RefNeg;
use serde::{Deserialize, Serialize};
use symbolica::atom::Atom;
use symbolica::domains::float::{NumericalFloatLike, Real};
use typed_index_collections::TiVec;

use crate::debug_info::DEBUG_LOGGER;
use crate::graph::{Graph, LoopMomentumBasis};
use crate::momentum::{FourMomentum, Signature, ThreeMomentum};
use crate::numerator::NumeratorState;
use crate::utils::{
    compute_loop_part, compute_shift_part, compute_t_part_of_shift_part, FloatLike, F,
};

use super::cff_graph::VertexSet;
use super::surface::{self, Surface};

/// Core esurface struct
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Esurface {
    pub energies: Vec<usize>,
    pub sub_orientation: Vec<bool>,
    pub external_shift: ExternalShift,
    pub circled_vertices: VertexSet,
}

impl Surface for Esurface {
    fn get_positive_energies(&self) -> impl Iterator<Item = &usize> {
        self.energies.iter()
    }

    fn get_external_shift(&self) -> impl Iterator<Item = &(usize, i64)> {
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
    pub fn to_atom(&self) -> Atom {
        let symbolic_energies = self
            .energies
            .iter()
            .map(|i| Atom::parse(&format!("E{}", i)).unwrap())
            .collect_vec();

        let symbolic_shift = self
            .external_shift
            .iter()
            .fold(Atom::new(), |sum, (i, sign)| {
                Atom::parse(&format!("p{}", i)).unwrap() * &Atom::new_num(*sign) + &sum
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
    pub fn compute_value<T: FloatLike>(&self, energy_cache: &[F<T>]) -> F<T> {
        surface::compute_value(self, energy_cache)
    }

    /// Only compute the shift part, useful for existence checks
    #[inline]
    pub fn compute_shift_part<T: FloatLike>(&self, energy_cache: &[F<T>]) -> F<T> {
        surface::compute_shift_part(self, energy_cache)
    }

    /// Compute the value of the esurface from the momenta, needed to check if an arbitrary point
    /// is inside the esurface
    #[inline]
    pub fn compute_from_momenta<T: FloatLike>(
        &self,
        lmb: &LoopMomentumBasis,
        real_mass_vector: &[F<T>],
        loop_moms: &[ThreeMomentum<F<T>>],
        external_moms: &[FourMomentum<F<T>>],
    ) -> F<T> {
        let spatial_part_of_externals = external_moms
            .iter()
            .map(|mom| mom.spatial.clone())
            .collect_vec();

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
            .unwrap_or_else(|| loop_moms[0].px.zero());

        energy_sum + self.compute_shift_part_from_momenta(lmb, external_moms)
    }

    /// Only compute the shift part, useful for center finding.
    #[inline]
    pub fn compute_shift_part_from_momenta<T: FloatLike>(
        &self,
        lmb: &LoopMomentumBasis,
        external_moms: &[FourMomentum<F<T>>],
    ) -> F<T> {
        self.external_shift
            .iter()
            .map(|(index, sign)| {
                let external_signature = &lmb.edge_signatures[*index].external;
                F::from_f64(*sign as f64)
                    * compute_t_part_of_shift_part(external_signature, external_moms)
            })
            .reduce(|acc, x| acc + x)
            .unwrap_or_else(|| external_moms[0].temporal.value.zero())
    }

    #[inline]
    pub fn compute_self_and_r_derivative<T: FloatLike>(
        &self,
        radius: &F<T>,
        shifted_unit_loops: &[ThreeMomentum<F<T>>],
        center: &[ThreeMomentum<F<T>>],
        external_moms: &[FourMomentum<F<T>>],
        lmb: &LoopMomentumBasis,
        real_mass_vector: &[F<T>],
    ) -> (F<T>, F<T>) {
        let spatial_part_of_externals = external_moms
            .iter()
            .map(|mom| mom.spatial.clone())
            .collect_vec();

        let loops = shifted_unit_loops
            .iter()
            .zip(center)
            .map(|(loop_mom, center)| loop_mom * radius + center)
            .collect_vec();

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

    #[inline]
    pub fn get_radius_guess<T: FloatLike>(
        &self,
        unit_loops: &[ThreeMomentum<F<T>>],
        external_moms: &[FourMomentum<F<T>>],
        lmb: &LoopMomentumBasis,
    ) -> (F<T>, F<T>) {
        let const_builder = &unit_loops[0].px;

        let esurface_shift = self.compute_shift_part_from_momenta(lmb, external_moms);

        let mut radius_guess = const_builder.zero();
        let mut denominator = const_builder.zero();

        for energy in self.energies.iter() {
            let signature = &lmb.edge_signatures[*energy];

            let unit_loop_part = compute_loop_part(&signature.internal, unit_loops);
            let shift = compute_shift_part(&signature.external, external_moms);
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

    pub fn canonicalize_shift(&mut self, dep_mom: usize, dep_mom_expr: &ExternalShift) {
        if let Some(dep_mom_pos) = self
            .external_shift
            .iter()
            .position(|(index, _)| *index == dep_mom)
        {
            let (_, dep_mom_sign) = self.external_shift.remove(dep_mom_pos);

            let external_shift = dep_mom_expr
                .iter()
                .map(|(index, sign)| (*index, dep_mom_sign * sign))
                .collect();

            self.external_shift = add_external_shifts(&self.external_shift, &external_shift);
        }
    }
}

pub type EsurfaceCollection = TiVec<EsurfaceID, Esurface>;

pub fn compute_esurface_cache<T: FloatLike>(
    esurfaces: &EsurfaceCollection,
    energy_cache: &[F<T>],
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
    externals: &[FourMomentum<F<T>>],
    lmb: &LoopMomentumBasis,
    debug: usize,
    e_cm: F<f64>,
) -> ExistingEsurfaces {
    if lmb.basis.is_empty() {
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
            let shift_signature = &esurface_derived_data[esurface_to_check_id].shift_signature;

            let esurface_shift = compute_shift_part(shift_signature, externals);

            let shift_spatial_sq = esurface_shift.spatial.norm_squared();

            let mass_sum_squared = esurface_derived_data[esurface_to_check_id].mass_sum_squared;

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

#[derive(Clone, Serialize, Deserialize, Debug, Encode, Decode)]
pub struct EsurfaceDerivedData {
    esurface_data: Vec<EsurfaceData>,
    orientation_pairs: Vec<(EsurfaceID, EsurfaceID)>,
}

impl Index<EsurfaceID> for EsurfaceDerivedData {
    type Output = EsurfaceData;

    fn index(&self, index: EsurfaceID) -> &Self::Output {
        &self.esurface_data[index.0]
    }
}

#[derive(Clone, Serialize, Deserialize, Debug, Encode, Decode)]
pub struct EsurfaceData {
    cut_momentum_basis: usize,
    mass_sum_squared: F<f64>,
    shift_signature: Signature,
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

pub fn generate_esurface_data<S: NumeratorState>(
    graph: &Graph<S>,
    esurfaces: &EsurfaceCollection,
) -> Result<EsurfaceDerivedData, Report> {
    let lmbs = graph
        .derived_data
        .as_ref()
        .unwrap()
        .loop_momentum_bases
        .as_ref()
        .ok_or_else(|| {
            eyre!("Could not generate esurface data, loop momentum bases not generated.")
        })?;

    let data = esurfaces
        .iter()
        .map(|esurface| {
            // find the cut momentum basis
            let energies = &esurface.energies;

            let (lmb_index, lmb) = lmbs
                .iter()
                .enumerate()
                .find(|(_i, lmb)| {
                    lmb.basis.iter().filter(|&i| energies.contains(i)).count() == energies.len() - 1
                })
                .ok_or_else(|| eyre!("Could not find a cut momentum basis for esurface"))?;

            let energy_not_in_cmb = *energies
                .iter()
                .find(|&i| !lmb.basis.contains(i))
                .ok_or_else(|| eyre!("No remaining edge in esurface"))?;

            let shift_signature = lmb.edge_signatures[energy_not_in_cmb].external.clone();

            let mass_sum: F<f64> = esurface
                .energies
                .iter()
                .map(|&i| graph.bare_graph.edges[i].particle.mass.value)
                .filter(|mass| mass.is_some())
                .map(|mass| mass.unwrap_or_else(|| unreachable!()).re)
                .reduce(|acc, x| acc + x)
                .unwrap_or_else(|| F::from_f64(0.0));

            Ok(EsurfaceData {
                cut_momentum_basis: lmb_index,
                mass_sum_squared: mass_sum * mass_sum,
                shift_signature,
            })
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

        let (position, other_esurface_id) = esurface_ids
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
            .ok_or_else(|| {
                eyre!(
                    "Could not find mirror esurface for esurface: {:?}, list of esurfaces: {:#?}",
                    esurface_id,
                    esurfaces
                )
            })?;

        orientation_pairs.push((esurface_id, *other_esurface_id));
        esurface_ids.remove(position);
    }

    Ok(EsurfaceDerivedData {
        esurface_data: data,
        orientation_pairs,
    })
}

pub type ExternalShift = Vec<(usize, i64)>;

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

#[cfg(test)]
mod tests {
    use symbolica::atom::{Atom, AtomCore};

    use crate::{
        cff::{cff_graph::VertexSet, esurface::Esurface},
        utils::F,
    };

    use super::add_external_shifts;

    #[test]
    fn test_esurface() {
        let energies_cache = [F(1.), F(2.), F(3.), F(4.), F(5.)];
        let energies = vec![0, 1, 2];

        let external_shift = vec![(3, 1), (4, 1)];

        let sub_orientation = vec![true, false, true];

        let dummy_circled_vertices = VertexSet::dummy();

        let mut esurface = Esurface {
            sub_orientation: sub_orientation.clone(),
            energies,
            external_shift,
            circled_vertices: dummy_circled_vertices,
        };

        let res = esurface.compute_value(&energies_cache);
        assert_eq!(res.0, 15.);

        let canon_shift = vec![(1, -1), (2, -1), (3, -1)];
        esurface.canonicalize_shift(4, &canon_shift);

        assert_eq!(esurface.external_shift, vec![(1, -1), (2, -1)]);

        let energies = vec![0, 2];

        let external_shift = vec![(1, -1)];

        let esurface = Esurface {
            energies,
            sub_orientation: sub_orientation.clone(),
            external_shift,
            circled_vertices: dummy_circled_vertices,
        };

        let res = esurface.compute_value(&energies_cache);
        assert_eq!(res.0, 2.);
    }

    #[test]
    fn test_to_atom() {
        let external_shift = vec![(1, -1)];

        let esurface = Esurface {
            sub_orientation: vec![true, true],
            energies: vec![2, 3],
            external_shift,
            circled_vertices: VertexSet::dummy(),
        };

        let esurface_atom = esurface.to_atom();
        let expected_atom = Atom::parse("E2 + E3 - p1").unwrap();

        let diff = esurface_atom - &expected_atom;
        let diff = diff.expand();
        assert_eq!(diff, Atom::new());
    }

    #[test]
    fn test_add_external_shifts() {
        let shift_1 = vec![(0, 1), (1, 1), (2, -1)];
        let shift_2 = vec![(1, -1), (2, 1)];

        let add = add_external_shifts(&shift_1, &shift_2);

        assert_eq!(add, vec![(0, 1)]);

        let shift_3 = vec![(3, 1), (4, -1)];
        let shift_4 = vec![(0, 1), (1, 1), (2, 1), (4, 1)];

        let add = add_external_shifts(&shift_3, &shift_4);

        assert_eq!(add, vec![(0, 1), (1, 1), (2, 1), (3, 1)]);
    }

    #[test]
    fn test_esurface_equality() {
        let esurface_1 = Esurface {
            energies: vec![3, 5],
            sub_orientation: vec![true, false],
            external_shift: vec![(0, 1), (1, 1)],
            circled_vertices: VertexSet::dummy(),
        };

        let esurface_2 = Esurface {
            energies: vec![3, 5],
            sub_orientation: vec![true, false],
            external_shift: vec![(0, 1), (1, 1)],
            circled_vertices: VertexSet::from_usize(1),
        };

        assert_eq!(esurface_1, esurface_2);
    }
}
