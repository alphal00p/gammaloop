use std::ops::Index;

use color_eyre::Report;
use derive_more::{From, Into};
use eyre::eyre;
use itertools::Itertools;
use lorentz_vector::LorentzVector;
use serde::{Deserialize, Serialize};
use symbolica::representations::Atom;
use typed_index_collections::TiVec;

use crate::graph::{Graph, LoopMomentumBasis};
use crate::utils::{
    compute_loop_part, compute_momentum, compute_shift_part, compute_t_part_of_shift_part,
    format_momentum, FloatLike,
};

use super::cff_graph::VertexSet;
use super::surface::{self, Surface};

/// Core esurface struct
#[derive(Serialize, Deserialize, Debug, Clone, Eq)]
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

// This equality is naive in the presence of raised propagators
impl PartialEq for Esurface {
    fn eq(&self, other: &Self) -> bool {
        self.energies == other.energies && self.sub_orientation == other.sub_orientation
    }
}

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
    pub fn compute_value<T: FloatLike>(&self, energy_cache: &[T]) -> T {
        surface::compute_value(self, energy_cache)
    }

    /// Only compute the shift part, useful for existence checks
    #[inline]
    pub fn compute_shift_part<T: FloatLike>(&self, energy_cache: &[T]) -> T {
        surface::compute_shift_part(self, energy_cache)
    }

    /// Compute the value of the esurface from the momenta, needed to check if an arbitrary point
    /// is inside the esurface
    #[inline]
    pub fn compute_from_momenta<T: FloatLike>(
        &self,
        lmb: &LoopMomentumBasis,
        real_mass_vector: &[T],
        loop_moms: &[LorentzVector<T>],
        external_moms: &[LorentzVector<T>],
    ) -> T {
        let energy_sum = self
            .energies
            .iter()
            .map(|index| {
                let signature = &lmb.edge_signatures[*index];
                let momentum = compute_momentum(signature, loop_moms, external_moms);
                let mass = real_mass_vector[*index];

                (momentum.spatial_squared() + mass * mass).sqrt()
            })
            .sum::<T>();

        energy_sum + self.compute_shift_part_from_momenta(lmb, external_moms)
    }

    /// Only compute the shift part, useful for center finding.
    #[inline]
    pub fn compute_shift_part_from_momenta<T: FloatLike>(
        &self,
        lmb: &LoopMomentumBasis,
        external_moms: &[LorentzVector<T>],
    ) -> T {
        self.external_shift
            .iter()
            .map(|(index, sign)| {
                let external_signature = &lmb.edge_signatures[*index].1;
                Into::<T>::into(*sign as f64)
                    * compute_t_part_of_shift_part(external_signature, external_moms)
            })
            .sum::<T>()
    }

    #[inline]
    pub fn compute_self_and_r_derivative<T: FloatLike>(
        &self,
        radius: T,
        shifted_unit_loops: &[LorentzVector<T>],
        center: &[LorentzVector<T>],
        external_moms: &[LorentzVector<T>],
        lmb: &LoopMomentumBasis,
        real_mass_vector: &[T],
    ) -> (T, T) {
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

                let momentum = compute_momentum(signature, &loops, external_moms);
                let unit_loop_part = compute_loop_part(&signature.0, shifted_unit_loops);

                let energy = (momentum.spatial_squared()
                    + real_mass_vector[index] * real_mass_vector[index])
                    .sqrt();

                let numerator = momentum.spatial_dot(&unit_loop_part);

                (numerator / energy, energy)
            })
            .fold((T::zero(), T::zero()), |(der_sum, en_sum), (der, en)| {
                (der_sum + der, en_sum + en)
            });

        (energy_sum + shift, derivative)
    }

    #[inline]
    pub fn get_radius_guess<T: FloatLike>(
        &self,
        _unit_loops: &[LorentzVector<T>],
        external_moms: &[LorentzVector<T>],
        lmb: &LoopMomentumBasis,
    ) -> T {
        //let mut radius_guess = T::zero();
        //let mut denominator = T::zero();

        let esurface_shift = self.compute_shift_part_from_momenta(lmb, external_moms);
        Into::<T>::into(2.0) * esurface_shift.abs()

        //for energy in self.energies.iter() {
        //    let signature = &lmb.edge_signatures[*energy];

        //    let unit_loop_part = compute_loop_part(&signature.0, unit_loops);
        //    let shift = compute_shift_part(&signature.1, external_moms);
        //    let norm_unit_loop_part_squared = unit_loop_part.spatial_squared();
        //    let loop_dot_shift = unit_loop_part.spatial_dot(&shift);

        //    radius_guess += loop_dot_shift.abs() / norm_unit_loop_part_squared;
        //    denominator += norm_unit_loop_part_squared.sqrt();
        //}

        //radius_guess += esurface_shift / denominator;
        //radius_guess
    }

    /// Write out the esurface expression in a given lmb
    pub fn string_format_in_lmb(&self, lmb: &LoopMomentumBasis) -> String {
        let mut energy_sum = self
            .energies
            .iter()
            .map(|index| {
                let signature = &lmb.edge_signatures[*index];
                format!("|{}|", format_momentum(signature))
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
                format!(" {} ({})^0", sign, format_momentum(signature))
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
}

pub type EsurfaceCollection = TiVec<EsurfaceID, Esurface>;

pub fn compute_esurface_cache<T: FloatLike>(
    esurfaces: &EsurfaceCollection,
    energy_cache: &[T],
) -> EsurfaceCache<T> {
    esurfaces
        .iter()
        .map(|esurface| esurface.compute_value(energy_cache))
        .collect::<Vec<T>>()
        .into()
}

pub type EsurfaceCache<T> = TiVec<EsurfaceID, T>;

/// Index type for esurface, location of an esurface in the list of all esurfaces of a graph
#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq, From, Into, Eq)]
pub struct EsurfaceID(usize);

pub type ExistingEsurfaces = TiVec<ExistingEsurfaceId, EsurfaceID>;

/// Container for esurfaces that exist at a given point in the phase space

/// Index in the list of all existing esurfaces, essentially a pointer to a pointer to an esurface
#[derive(
    Debug, From, Into, Copy, Clone, Serialize, Deserialize, PartialEq, Eq, Hash, PartialOrd, Ord,
)]
pub struct ExistingEsurfaceId(usize);

const MAX_EXPECTED_CAPACITY: usize = 32; // Used to prevent reallocations during existence check
const SHIFT_THRESHOLD: f64 = 1.0e-13;
const EXISTENCE_THRESHOLD: f64 = 1.0e-7;

/// Returns the list of esurfaces which may exist, must be called each time at evaluation if externals are not fixed.
#[inline]
pub fn get_existing_esurfaces<T: FloatLike>(
    esurfaces: &EsurfaceCollection,
    esurface_derived_data: &EsurfaceDerivedData,
    externals: &[LorentzVector<T>],
    lmb: &LoopMomentumBasis,
    debug: usize,
    e_cm: f64,
) -> ExistingEsurfaces {
    let mut existing_esurfaces = ExistingEsurfaces::with_capacity(MAX_EXPECTED_CAPACITY);

    for orientation_pair in &esurface_derived_data.orientation_pairs {
        if let Some((esurface_to_check_id, shift_zero_sq)) = {
            let (esurface_id, other_esurface_id) = orientation_pair;

            let esurface = &esurfaces[*esurface_id];

            let shift_part = esurface.compute_shift_part_from_momenta(lmb, externals);
            let shift_zero_sq = shift_part * shift_part;

            if shift_part < -Into::<T>::into(SHIFT_THRESHOLD * e_cm) {
                Some((*esurface_id, shift_zero_sq))
            } else if shift_part > Into::<T>::into(SHIFT_THRESHOLD * e_cm) {
                Some((*other_esurface_id, shift_zero_sq))
            } else {
                None
            }
        } {
            let shift_signature = &esurface_derived_data[esurface_to_check_id].shift_signature;

            let esurface_shift = compute_shift_part(shift_signature, externals);

            let shift_spatial_sq = esurface_shift.spatial_squared();
            let mass_sum_squared = esurface_derived_data[esurface_to_check_id].mass_sum_squared;

            let existence_condition =
                shift_zero_sq - shift_spatial_sq - Into::<T>::into(mass_sum_squared);

            if existence_condition
                > Into::<T>::into(EXISTENCE_THRESHOLD * EXISTENCE_THRESHOLD * e_cm * e_cm)
            {
                if debug > 1 {
                    println!(
                        "existing esurface: {:?}, existence_condition: {}, shift_part_squared: {}",
                        esurface_to_check_id, existence_condition, shift_zero_sq
                    );
                }
                existing_esurfaces.push(esurface_to_check_id);
            }
        }
    }
    existing_esurfaces
}

#[derive(Clone, Serialize, Deserialize, Debug)]
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

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct EsurfaceData {
    cut_momentum_basis: usize,
    mass_sum_squared: f64,
    shift_signature: Vec<isize>,
}

impl EsurfaceData {
    #[allow(dead_code)]
    fn existence_condition<T: FloatLike>(&self, externals: &[LorentzVector<T>]) -> (T, T) {
        let mut shift = LorentzVector::new();

        for (i, external) in externals.iter().enumerate() {
            match self.shift_signature[i] {
                1 => shift += external,
                -1 => shift -= external,
                0 => {}
                _ => unreachable!("Shift signature must be -1, 0 or 1"),
            }
        }

        (
            shift.t,
            shift.square() - Into::<T>::into(self.mass_sum_squared),
        )
    }

    pub fn compute_shift_part_from_externals<T: FloatLike>(
        &self,
        externals: &[LorentzVector<T>],
    ) -> T {
        let mut shift = T::zero();

        for (i, external) in externals.iter().enumerate() {
            match self.shift_signature[i] {
                1 => shift += external.t,
                -1 => shift -= external.t,
                0 => {}
                _ => unreachable!("Shift signature must be -1, 0 or 1"),
            }
        }

        shift
    }
}

pub fn generate_esurface_data(
    graph: &Graph,
    esurfaces: &EsurfaceCollection,
) -> Result<EsurfaceDerivedData, Report> {
    let lmbs = graph
        .derived_data
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

            let shift_signature = lmb.edge_signatures[energy_not_in_cmb].1.clone();

            let mass_sum: f64 = esurface
                .energies
                .iter()
                .map(|&i| graph.edges[i].particle.mass.value)
                .filter(|mass| mass.is_some())
                .map(|mass| mass.unwrap_or_else(|| unreachable!()).re)
                .sum();

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
        let (posistion, other_esurface_id) = esurface_ids
            .iter()
            .enumerate()
            .find(|(_pos, other_esurface_id)| {
                let esurface = &esurfaces[esurface_id];
                let other_esurface = &esurfaces[**other_esurface_id];

                let energies_match = esurface.energies == other_esurface.energies;
                let orientation_flipped =
                    esurface.sub_orientation != other_esurface.sub_orientation;

                energies_match && orientation_flipped
            })
            .ok_or_else(|| eyre!("Could not find mirror esurface"))?;

        orientation_pairs.push((esurface_id, *other_esurface_id));
        esurface_ids.remove(posistion);
    }

    Ok(EsurfaceDerivedData {
        esurface_data: data,
        orientation_pairs,
    })
}

pub type ExternalShift = Vec<(usize, i64)>;

/// writes into lhs while emptying rhs
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
    use symbolica::representations::Atom;

    use crate::cff::{cff_graph::VertexSet, esurface::Esurface};

    use super::add_external_shifts;

    #[test]
    fn test_esurface() {
        let energies_cache = [1., 2., 3., 4., 5.];
        let energies = vec![0, 1, 2];

        let external_shift = vec![(3, 1), (4, 1)];

        let sub_orientation = vec![true, false, true];

        let dummy_circled_vertices = VertexSet::dummy();

        let esurface = Esurface {
            sub_orientation: sub_orientation.clone(),
            energies,
            external_shift,
            circled_vertices: dummy_circled_vertices,
        };

        let res = esurface.compute_value(&energies_cache);
        assert_eq!(res, 15.);

        let energies = vec![0, 2];

        let external_shift = vec![(1, -1)];

        let esurface = Esurface {
            energies,
            sub_orientation: sub_orientation.clone(),
            external_shift,
            circled_vertices: dummy_circled_vertices,
        };

        let res = esurface.compute_value(&energies_cache);
        assert_eq!(res, 2.);
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
}
