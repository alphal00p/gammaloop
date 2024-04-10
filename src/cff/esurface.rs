use std::ops::Index;

use color_eyre::Report;
use eyre::eyre;
use itertools::Itertools;
use lorentz_vector::LorentzVector;
use serde::{Deserialize, Serialize};
use symbolica::representations::Atom;

use crate::graph::{Graph, LoopMomentumBasis};
use crate::utils::{
    compute_momentum, compute_shift_part, compute_t_part_of_shift_part, format_momentum, FloatLike,
};

/// Core esurface struct
#[derive(Serialize, Deserialize, Debug, Clone, Eq)]
pub struct Esurface {
    pub energies: Vec<usize>,
    pub sub_orientation: Vec<bool>,
    pub shift: Vec<usize>,
    pub shift_signature: Vec<bool>,
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
            .shift
            .iter()
            .map(|i| Atom::parse(&format!("p{}", i)).unwrap())
            .collect_vec();

        let builder_atom = Atom::new();
        let energy_sum = symbolic_energies
            .iter()
            .fold(builder_atom, |acc, energy| acc + energy);

        let esurf = symbolic_shift.iter().zip(self.shift_signature.iter()).fold(
            energy_sum,
            |acc, (shift, &shift_signature)| {
                if shift_signature {
                    acc + shift
                } else {
                    acc - shift
                }
            },
        );

        esurf
    }

    /// Compute the value of the esurface from an energy cache that can be computed from the underlying graph
    /// This is the fastest way to compute the value of all esurfaces in a full evaluation
    #[inline]
    pub fn compute_value<T: FloatLike>(&self, energy_cache: &[T]) -> T {
        let energy_sum = self
            .energies
            .iter()
            .map(|index| energy_cache[*index])
            .sum::<T>();

        energy_sum + self.compute_shift_part(energy_cache)
    }

    /// Only compute the shift part, useful for existence checks
    #[inline]
    pub fn compute_shift_part<T: FloatLike>(&self, energy_cache: &[T]) -> T {
        let shift_sum = self.shift.iter().zip(self.shift_signature.iter()).fold(
            T::zero(),
            |acc, (index, sign)| match sign {
                true => acc + energy_cache[*index],
                false => acc - energy_cache[*index],
            },
        );

        shift_sum
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
        self.shift
            .iter()
            .zip(self.shift_signature.iter())
            .map(|(index, sign)| {
                let external_signature = &lmb.edge_signatures[*index].1;
                match sign {
                    true => compute_t_part_of_shift_part(external_signature, external_moms),
                    false => -compute_t_part_of_shift_part(external_signature, external_moms),
                }
            })
            .sum::<T>()
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
            .shift_signature
            .iter()
            .zip(self.shift.iter())
            .map(|(sign, index)| {
                let signature = &lmb.edge_signatures[*index];
                let sign = if *sign { "+" } else { "-" };
                format!(" {} ({})^0", sign, format_momentum(signature))
            })
            .join("");

        energy_sum.push_str(&shift_part);
        energy_sum
    }

    /// Write out the esurface expression in a generic way
    pub fn string_format(&self) -> String {
        let mut energy_sum = self
            .energies
            .iter()
            .map(|index| format!("E{}", index))
            .join(" + ");

        let shift_part = self
            .shift_signature
            .iter()
            .zip(self.shift.iter())
            .map(|(sign, index)| {
                let sign = if *sign { "+" } else { "-" };
                format!(" {} p{}", sign, index)
            })
            .join("");

        energy_sum.push_str(&shift_part);
        energy_sum
    }

    pub fn get_point_inside(&self) -> Vec<LorentzVector<f64>> {
        todo!()
    }
}

/// Container for esurfaces, supposed to be used for the list of all esurface of a graph
#[derive(Debug, Clone)]
pub struct EsurfaceCollection {
    esurfaces: Vec<Esurface>,
}

#[derive(Debug, Clone)]
pub struct EsurfaceCache<T> {
    cache: Vec<T>,
}

impl<T> Index<EsurfaceId> for EsurfaceCache<T> {
    type Output = T;

    fn index(&self, index: EsurfaceId) -> &Self::Output {
        &self.cache[index.0]
    }
}

impl<T> FromIterator<T> for EsurfaceCache<T> {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        Self {
            cache: iter.into_iter().collect(),
        }
    }
}

/// Index type for esurface, location of an esurface in the list of all esurfaces of a graph
#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct EsurfaceId(usize);

impl From<EsurfaceId> for usize {
    fn from(id: EsurfaceId) -> Self {
        id.0
    }
}

impl From<usize> for EsurfaceId {
    fn from(id: usize) -> Self {
        Self(id)
    }
}

impl EsurfaceCollection {
    pub fn to_vec(self) -> Vec<Esurface> {
        self.esurfaces
    }

    pub fn from_vec(esurfaces: Vec<Esurface>) -> Self {
        Self { esurfaces }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.esurfaces.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.esurfaces.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = &Esurface> {
        self.esurfaces.iter()
    }

    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut Esurface> {
        self.esurfaces.iter_mut()
    }

    pub fn push(&mut self, esurface: Esurface) {
        self.esurfaces.push(esurface);
    }

    pub fn iterate_all_ids(&self) -> impl Iterator<Item = EsurfaceId> {
        (0..self.len()).map(EsurfaceId)
    }
}

impl Index<EsurfaceId> for EsurfaceCollection {
    type Output = Esurface;

    fn index(&self, index: EsurfaceId) -> &Self::Output {
        &self.esurfaces[index.0]
    }
}

/// Container for esurfaces that exist at a given point in the phase space
#[derive(Debug, Clone)]
pub struct ExistingEsurfaces {
    existing_esurfaces: Vec<EsurfaceId>,
}

impl ExistingEsurfaces {
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            existing_esurfaces: Vec::with_capacity(capacity),
        }
    }

    pub fn push(&mut self, esurface_id: EsurfaceId) {
        self.existing_esurfaces.push(esurface_id);
    }

    pub fn len(&self) -> usize {
        self.existing_esurfaces.len()
    }

    pub fn is_empty(&self) -> bool {
        self.existing_esurfaces.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = &EsurfaceId> {
        self.existing_esurfaces.iter()
    }

    pub fn to_vec(self) -> Vec<EsurfaceId> {
        self.existing_esurfaces
    }

    pub fn from_vec(existing_esurfaces: Vec<EsurfaceId>) -> Self {
        Self { existing_esurfaces }
    }

    pub fn as_slice(&self) -> &[EsurfaceId] {
        &self.existing_esurfaces
    }

    pub fn all_indices(&self) -> impl Iterator<Item = ExistingEsurfaceId> {
        (0..self.len()).map(ExistingEsurfaceId)
    }

    pub fn iter_enumerate(&self) -> impl Iterator<Item = (ExistingEsurfaceId, &EsurfaceId)> {
        self.iter()
            .enumerate()
            .map(|(i, id)| (ExistingEsurfaceId(i), id))
    }
}

/// Index in the list of all existing esurfaces, essentially a pointer to a pointer to an esurface
#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct ExistingEsurfaceId(usize);

impl Index<ExistingEsurfaceId> for ExistingEsurfaces {
    type Output = EsurfaceId;

    fn index(&self, index: ExistingEsurfaceId) -> &Self::Output {
        &self.existing_esurfaces[index.0]
    }
}

impl From<ExistingEsurfaceId> for usize {
    fn from(id: ExistingEsurfaceId) -> Self {
        id.0
    }
}

impl From<usize> for ExistingEsurfaceId {
    fn from(id: usize) -> Self {
        Self(id)
    }
}

const MAX_EXPECTED_CAPACITY: usize = 32; // Used to prevent reallocations during existence check
const SHIFT_THRESHOLD: f64 = 1.0e-13;
const EXISTENCE_THRESHOLD: f64 = 1.0e-7;

/// Returns the list of esurfaces which may exist, must be called each time at evaluation if externals are not fixed.
#[inline]
pub fn get_existing_esurfaces<T: FloatLike>(
    esurfaces: &EsurfaceCollection,
    esurface_derived_data: &EsurfaceDerivedData,
    externals: &[LorentzVector<T>],
    energy_cache: &[T],
    debug: usize,
    e_cm: f64,
) -> ExistingEsurfaces {
    let mut existing_esurfaces = ExistingEsurfaces::with_capacity(MAX_EXPECTED_CAPACITY);

    for orientation_pair in &esurface_derived_data.orientation_pairs {
        if let Some((esurface_to_check_id, shift_zero_sq)) = {
            let (esurface_id, other_esurface_id) = orientation_pair;

            let esurface = &esurfaces[*esurface_id];

            let shift_part = esurface.compute_shift_part(energy_cache);
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
    orientation_pairs: Vec<(EsurfaceId, EsurfaceId)>,
}

impl Index<EsurfaceId> for EsurfaceDerivedData {
    type Output = EsurfaceData;

    fn index(&self, index: EsurfaceId) -> &Self::Output {
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

pub fn generate_esurface_data(graph: &Graph) -> Result<EsurfaceDerivedData, Report> {
    let lmbs = graph
        .derived_data
        .loop_momentum_bases
        .as_ref()
        .ok_or_else(|| {
            eyre!("Could not generate esurface data, loop momentum bases not generated.")
        })?;
    let cff =
        graph.derived_data.cff_expression.as_ref().ok_or_else(|| {
            eyre!("Could not generate esurface data, cff expression not generated.")
        })?;

    let data = cff
        .esurfaces
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

    let mut esurface_ids = cff.esurfaces.iterate_all_ids().collect_vec();
    let mut orientation_pairs = Vec::with_capacity(esurface_ids.len() / 2);

    while let Some(esurface_id) = esurface_ids.pop() {
        let (posistion, other_esurface_id) = esurface_ids
            .iter()
            .enumerate()
            .find(|(_pos, other_esurface_id)| {
                let esurface = &cff.esurfaces[esurface_id];
                let other_esurface = &cff.esurfaces[**other_esurface_id];

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
