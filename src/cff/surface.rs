use super::{
    esurface::{Esurface, EsurfaceID},
    hsurface::{Hsurface, HsurfaceID},
};
use crate::utils::FloatLike;
use derive_more::From;
use enum_dispatch::enum_dispatch;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use typed_index_collections::TiVec;

pub trait Surface {
    fn get_positive_energies(&self) -> impl Iterator<Item = &usize>;

    fn get_negative_energies(&self) -> impl Iterator<Item = &usize> {
        std::iter::empty::<&usize>()
    }

    fn get_external_shift(&self) -> impl Iterator<Item = (&usize, &bool)>;
}

#[inline]
pub fn compute_value<T: FloatLike, S: Surface>(surface: &S, energy_cache: &[T]) -> T {
    let positive_energies = surface.get_positive_energies();
    let negative_energies = surface.get_negative_energies();

    let positive_energy_sum = positive_energies
        .map(|&index| energy_cache[index])
        .sum::<T>();

    let negative_energy_sum = negative_energies
        .map(|&index| energy_cache[index])
        .sum::<T>();

    let shift_part = compute_shift_part(surface, energy_cache);

    positive_energy_sum - negative_energy_sum + shift_part
}

#[inline]
pub fn compute_shift_part<T: FloatLike, S: Surface>(surface: &S, energy_cache: &[T]) -> T {
    surface
        .get_external_shift()
        .map(|(&index, &sign)| match sign {
            true => energy_cache[index],
            false => -energy_cache[index],
        })
        .sum::<T>()
}

pub fn string_format<S: Surface>(surface: &S) -> String {
    let positive_energy_sum = surface
        .get_positive_energies()
        .map(|index| format!("E_{}", index))
        .join(" + ");

    let negative_energy_sum = surface
        .get_negative_energies()
        .map(|index| format!("E_{}", index))
        .join(" - ");

    let shift_part = surface
        .get_external_shift()
        .map(|(index, &sign)| {
            let sign = if sign { "+" } else { "-" };
            format!(" {} p{}", sign, index)
        })
        .join("");

    let mut res = String::new();
    res.push_str(&positive_energy_sum);
    res.push_str(" - ");
    res.push_str(&negative_energy_sum);
    res.push_str(&shift_part);

    res
}

#[enum_dispatch]
pub enum HybridSurface {
    Esurface(Esurface),
    Hsurface(Hsurface),
}

#[derive(From, Clone, Copy, Debug, Serialize, Deserialize)]
pub enum HybridSurfaceID {
    Esurface(EsurfaceID),
    Hsurface(HsurfaceID),
}

pub type HybridSurfaceCollection = TiVec<HybridSurfaceID, HybridSurface>;
pub type HybridSurfaceCache<T> = TiVec<HybridSurfaceID, T>;
