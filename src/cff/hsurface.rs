use serde::{Deserialize, Serialize};

use crate::utils::FloatLike;

use crate::cff::surface::Surface;

use super::surface;

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Hsurface {
    pub positive_energies: Vec<usize>,
    pub negative_energies: Vec<usize>,
    pub sub_orientation: Vec<bool>,
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
            && self.sub_orientation == other.sub_orientation
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

#[cfg(test)]
mod tests {
    use super::Hsurface;

    #[test]
    fn test_compute_shift_part() {
        let h_surface = Hsurface {
            positive_energies: vec![0],
            negative_energies: vec![1],
            sub_orientation: vec![true, true],
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
            sub_orientation: vec![true, true, true, true],
            shift: vec![4, 5],
            shift_signature: vec![false, true],
        };

        let energy_cache = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let value = h_surface.compute_value(&energy_cache);
        assert_eq!(value, 1.0 + 2.0 - 3.0 - 4.0 - 5.0 + 6.0);
    }
}
