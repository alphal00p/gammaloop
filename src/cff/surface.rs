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

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct UnitSurface {}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum HybridSurface {
    Esurface(Esurface),
    Hsurface(Hsurface),
    Unit(UnitSurface),
}

impl HybridSurface {
    pub(crate) fn to_atom(&self, cut_edges: &[EdgeIndex]) -> Atom {
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
    pub(crate) fn to_atom(&self, cut_edges: &[EdgeIndex]) -> Atom {
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
