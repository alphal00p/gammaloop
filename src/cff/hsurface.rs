use serde::{Deserialize, Serialize};

use crate::utils::FloatLike;

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct HSurface {
    pub positive_energies: Vec<usize>,
    pub negative_energies: Vec<usize>,
    pub sub_orientation: Vec<bool>,
    pub shift: Vec<usize>,
    pub shift_signature: Vec<bool>,
}

impl PartialEq for HSurface {
    fn eq(&self, other: &Self) -> bool {
        self.positive_energies == other.positive_energies
            && self.sub_orientation == other.sub_orientation
    }
}

impl HSurface {
    #[inline]
    pub fn compute_value<T: FloatLike>(&self, energy_cache: &[T]) -> T {
        let positive_energy_sum = self
            .positive_energies
            .iter()
            .map(|index| energy_cache[*index])
            .sum::<T>();

        let negative_energy_sum = self
            .negative_energies
            .iter()
            .map(|index| energy_cache[*index])
            .sum::<T>();

        let shift_part = self.compute_shift_part(energy_cache);

        positive_energy_sum - negative_energy_sum + shift_part
    }

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
}

#[cfg(test)]
mod tests {
    use super::HSurface;

    #[test]
    fn test_compute_shift_part() {
        let h_surface = HSurface {
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
        let h_surface = HSurface {
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
