use crate::cff::surface::Surface;
use crate::utils::FloatLike;
use derive_more::{From, Into};
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use typed_index_collections::TiVec;

use super::surface;

#[derive(From, Into, Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct HsurfaceID(usize);
pub type HsurfaceCollection = TiVec<HsurfaceID, Hsurface>;
pub type HsurfaceCache<T> = TiVec<HsurfaceID, T>;

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Hsurface {
    pub positive_energies: Vec<usize>,
    pub negative_energies: Vec<usize>,
    pub shift: Vec<usize>,
    pub shift_signature: Vec<bool>,
}

impl Surface for Hsurface {
    fn get_positive_energies(&self) -> impl Iterator<Item = &usize> {
        self.positive_energies.iter()
    }

    fn get_negative_energies(&self) -> impl Iterator<Item = &usize> {
        self.negative_energies.iter()
    }

    fn get_external_shift(&self) -> impl Iterator<Item = (&usize, &bool)> {
        self.shift.iter().zip(&self.shift_signature)
    }
}

impl PartialEq for Hsurface {
    fn eq(&self, other: &Self) -> bool {
        self.positive_energies == other.positive_energies
            && self.negative_energies == other.negative_energies
    }
}

impl Hsurface {
    #[inline]
    pub fn compute_value<T: FloatLike>(&self, energy_cache: &[T]) -> T {
        surface::compute_value(self, energy_cache)
    }

    #[inline]
    pub fn compute_shift_part<T: FloatLike>(&self, energy_cache: &[T]) -> T {
        surface::compute_shift_part(self, energy_cache)
    }
}

pub fn compute_hsurface_cache<T: FloatLike>(
    hsurfaces: &HsurfaceCollection,
    energy_cache: &[T],
) -> HsurfaceCache<T> {
    hsurfaces
        .iter()
        .map(|hsurface| hsurface.compute_value(energy_cache))
        .collect_vec()
        .into()
}

#[cfg(test)]
mod tests {
    use super::Hsurface;

    #[test]
    fn test_compute_shift_part() {
        let h_surface = Hsurface {
            positive_energies: vec![0],
            negative_energies: vec![1],
            shift: vec![2, 3],
            shift_signature: vec![true, false],
        };

        let energy_cache = [1.0, 2.0, 3.0, 4.0];
        let shift_part = h_surface.compute_shift_part(&energy_cache);
        assert_eq!(shift_part, -1.0);
    }

    #[test]
    fn test_compute_value() {
        let h_surface = Hsurface {
            positive_energies: vec![0, 1],
            negative_energies: vec![2, 3],
            shift: vec![4, 5],
            shift_signature: vec![false, true],
        };

        let energy_cache = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let value = h_surface.compute_value(&energy_cache);
        assert_eq!(value, 1.0 + 2.0 - 3.0 - 4.0 - 5.0 + 6.0);
    }
}
