use super::{
    esurface::{Esurface, EsurfaceID},
    hsurface::{Hsurface, HsurfaceID},
};
use crate::utils::{FloatLike, F};
use derive_more::From;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use symbolica::domains::float::NumericalFloatLike;
use typed_index_collections::TiVec;

pub trait Surface {
    fn get_positive_energies(&self) -> impl Iterator<Item = &usize>;

    fn get_negative_energies(&self) -> impl Iterator<Item = &usize> {
        std::iter::empty::<&usize>()
    }

    fn get_external_shift(&self) -> impl Iterator<Item = &(usize, i64)>;
}

pub type UnitSurface = ();

impl Surface for UnitSurface {
    fn get_positive_energies(&self) -> impl Iterator<Item = &usize> {
        std::iter::empty::<&usize>()
    }

    fn get_external_shift(&self) -> impl Iterator<Item = &(usize, i64)> {
        std::iter::empty::<&(usize, i64)>()
    }
}

#[inline]
pub fn compute_value<T: FloatLike, S: Surface>(surface: &S, energy_cache: &[F<T>]) -> F<T> {
    let positive_energies = surface.get_positive_energies();
    let negative_energies = surface.get_negative_energies();

    let positive_energy_sum = positive_energies.fold(energy_cache[0].zero(), |acc, &index| {
        &acc + &energy_cache[index]
    });

    let negative_energy_sum = negative_energies.fold(energy_cache[0].zero(), |acc, &index| {
        &acc + &energy_cache[index]
    });

    let shift_part = compute_shift_part(surface, energy_cache);

    positive_energy_sum - negative_energy_sum + shift_part
}

#[inline]
pub fn compute_shift_part<T: FloatLike, S: Surface>(surface: &S, energy_cache: &[F<T>]) -> F<T> {
    surface
        .get_external_shift()
        .map(|&(index, sign)| F::from_f64(sign as f64) * &energy_cache[index])
        .reduce(|acc, x| &acc + &x)
        .unwrap_or_else(|| energy_cache[0].zero())
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
        .map(|&(index, sign)| {
            let sign = if sign == 1 {
                "+".to_owned()
            } else if sign == -1 {
                "-".to_owned()
            } else {
                format!("+{}", sign)
            };

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

pub enum HybridSurface {
    Esurface(Esurface),
    Hsurface(Hsurface),
    Unit(UnitSurface),
}

#[derive(From, Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub enum HybridSurfaceID {
    Esurface(EsurfaceID),
    Hsurface(HsurfaceID),
    Unit,
}

pub type HybridSurfaceCollection = TiVec<HybridSurfaceID, HybridSurface>;
pub type HybridSurfaceCache<T> = TiVec<HybridSurfaceID, T>;
