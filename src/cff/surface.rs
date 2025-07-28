use super::{
    esurface::{Esurface, EsurfaceID},
    hsurface::{Hsurface, HsurfaceID},
};
use crate::utils::{FloatLike, F};
use bincode_trait_derive::{Decode, Encode};
use derive_more::From;
use itertools::Itertools;
use linnet::half_edge::involution::{EdgeIndex, EdgeVec};
use serde::{Deserialize, Serialize};
use symbolica::{atom::Atom, domains::float::NumericalFloatLike};
use typed_index_collections::TiVec;

pub trait Surface {
    fn get_positive_energies(&self) -> impl Iterator<Item = &EdgeIndex>;

    fn get_negative_energies(&self) -> impl Iterator<Item = &EdgeIndex> {
        std::iter::empty::<&EdgeIndex>()
    }

    fn get_external_shift(&self) -> impl Iterator<Item = &(EdgeIndex, i64)>;
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct UnitSurface {}

impl Surface for UnitSurface {
    fn get_positive_energies(&self) -> impl Iterator<Item = &EdgeIndex> {
        std::iter::empty::<&EdgeIndex>()
    }

    fn get_external_shift(&self) -> impl Iterator<Item = &(EdgeIndex, i64)> {
        std::iter::empty::<&(EdgeIndex, i64)>()
    }
}

#[inline]
pub fn compute_value<T: FloatLike, S: Surface>(surface: &S, energy_cache: &EdgeVec<F<T>>) -> F<T> {
    let positive_energies = surface.get_positive_energies();
    let negative_energies = surface.get_negative_energies();

    let zero_index = EdgeIndex::from(0);

    let positive_energy_sum = positive_energies
        .fold(energy_cache[zero_index].zero(), |acc, &index| {
            &acc + &energy_cache[index]
        });

    let negative_energy_sum = negative_energies
        .fold(energy_cache[zero_index].zero(), |acc, &index| {
            &acc + &energy_cache[index]
        });

    let shift_part = compute_shift_part(surface, energy_cache);

    positive_energy_sum - negative_energy_sum + shift_part
}

#[inline]
pub fn compute_shift_part<T: FloatLike, S: Surface>(
    surface: &S,
    energy_cache: &EdgeVec<F<T>>,
) -> F<T> {
    surface
        .get_external_shift()
        .map(|&(index, sign)| F::from_f64(sign as f64) * &energy_cache[index])
        .reduce(|acc, x| &acc + &x)
        .unwrap_or_else(|| energy_cache[EdgeIndex::from(0)].zero())
}

pub fn string_format<S: Surface>(surface: &S) -> String {
    let positive_energy_sum = surface
        .get_positive_energies()
        .map(|index| format!("E_{}", Into::<usize>::into(*index)))
        .join(" + ");

    let negative_energy_sum = surface
        .get_negative_energies()
        .map(|index| format!("E_{}", Into::<usize>::into(*index)))
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

            format!(" {} p{}", sign, Into::<usize>::into(index))
        })
        .join("");

    let mut res = String::new();
    res.push_str(&positive_energy_sum);
    res.push_str(" - ");
    res.push_str(&negative_energy_sum);
    res.push_str(&shift_part);

    res
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum HybridSurface {
    Esurface(Esurface),
    Hsurface(Hsurface),
    Unit(UnitSurface),
}

impl HybridSurface {
    pub fn to_atom(&self, cut_edges: &[EdgeIndex]) -> Atom {
        match self {
            HybridSurface::Esurface(surface) => surface.to_atom(cut_edges),
            HybridSurface::Hsurface(surface) => surface.to_atom(cut_edges),
            HybridSurface::Unit(_) => Atom::num(1),
        }
    }
}

#[derive(Debug, Clone)]
pub enum HybridSurfaceRef<'a> {
    Esurface(&'a Esurface),
    Hsurface(&'a Hsurface),
    Unit(UnitSurface),
}

impl HybridSurfaceRef<'_> {
    pub fn to_atom(&self, cut_edges: &[EdgeIndex]) -> Atom {
        match self {
            HybridSurfaceRef::Esurface(surface) => surface.to_atom(cut_edges),
            HybridSurfaceRef::Hsurface(surface) => surface.to_atom(cut_edges),
            HybridSurfaceRef::Unit(_) => Atom::num(1),
        }
    }
}

#[derive(From, Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq, Encode, Decode)]
pub enum HybridSurfaceID {
    Esurface(EsurfaceID),
    Hsurface(HsurfaceID),
    Unit,
}

pub type HybridSurfaceCollection = TiVec<HybridSurfaceID, HybridSurface>;
pub type HybridSurfaceCache<T> = TiVec<HybridSurfaceID, T>;

impl From<HybridSurfaceID> for Atom {
    fn from(id: HybridSurfaceID) -> Atom {
        match id {
            HybridSurfaceID::Esurface(id) => Atom::from(id),
            HybridSurfaceID::Hsurface(id) => Atom::from(id),
            HybridSurfaceID::Unit => Atom::num(1),
        }
    }
}
