//! GammaLoop symbolic extensions for FeynKit-owned CFF surfaces.

use feynkit_cff::Surface;
use linnet::half_edge::involution::EdgeIndex;
use symbolica::{atom::Atom, parse};

use crate::{
    cff::{esurface::EnergySurfaceExt, hsurface::HSurfaceExt},
    graph::LoopMomentumBasis,
};

pub trait SurfaceExt {
    fn to_atom(&self, cut_edges: &[EdgeIndex]) -> Atom;
    fn to_atom_in_lmb(&self, cut_edges: &[EdgeIndex], lmb: &LoopMomentumBasis) -> Atom;
}

impl SurfaceExt for Surface {
    fn to_atom(&self, cut_edges: &[EdgeIndex]) -> Atom {
        match self {
            Self::Energy(surface) => surface.to_atom(cut_edges),
            Self::H(surface) => surface.to_atom(cut_edges),
            Self::Unit => Atom::num(1),
            Self::Infinite => parse!("η_inf"),
        }
    }

    fn to_atom_in_lmb(&self, cut_edges: &[EdgeIndex], lmb: &LoopMomentumBasis) -> Atom {
        match self {
            Self::Energy(surface) => surface.to_atom_in_lmb(cut_edges, lmb),
            Self::H(surface) => surface.to_atom_in_lmb(cut_edges, lmb),
            Self::Unit => Atom::num(1),
            Self::Infinite => parse!("η_inf"),
        }
    }
}
