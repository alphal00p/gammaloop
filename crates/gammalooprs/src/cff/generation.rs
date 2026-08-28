//! GammaLoop runtime extensions for FeynKit-owned CFF generation.

use color_eyre::Result;
use feynkit_cff::{CffResult, ShiftRewrite, SurfaceCache};
use linnet::half_edge::involution::EdgeIndex;
use symbolica::{
    atom::{Atom, AtomCore},
    id::{Pattern, Replacement},
};

use crate::{
    cff::{esurface::EnergySurfaceExt, hsurface::HSurfaceExt},
    graph::LoopMomentumBasis,
    settings::global::OrientationPattern,
};

/// GammaLoop entry points layered on the canonical FeynKit CFF generator.
pub(crate) trait GammaLoopGraphCffExt {
    fn generate_cff(
        &mut self,
        contract_edges: &[EdgeIndex],
        canonize_esurface: &Option<ShiftRewrite>,
        orientation_pattern: &OrientationPattern,
    ) -> Result<CffResult>;
}

/// GammaLoop symbolic lowering over the canonical FeynKit surface arena.
pub trait GammaLoopSurfaceCacheExt {
    fn substitute_energies(&self, atom: &Atom, cut_edges: &[EdgeIndex]) -> Atom;
    fn get_all_replacements(&self, cut_edges: &[EdgeIndex]) -> Vec<Replacement>;
    fn get_all_replacements_in_lmb(
        &self,
        cut_edges: &[EdgeIndex],
        lmb: &LoopMomentumBasis,
    ) -> Vec<Replacement>;
}

impl GammaLoopSurfaceCacheExt for SurfaceCache {
    fn substitute_energies(&self, atom: &Atom, cut_edges: &[EdgeIndex]) -> Atom {
        atom.replace_multiple(self.get_all_replacements(cut_edges))
    }

    fn get_all_replacements(&self, cut_edges: &[EdgeIndex]) -> Vec<Replacement> {
        let energy = self
            .energy_surfaces()
            .iter()
            .enumerate()
            .map(|(index, surface)| {
                Replacement::new(
                    Pattern::from(Atom::from(feynkit_cff::EnergySurfaceId::from(index))),
                    Pattern::from(surface.to_atom(cut_edges)),
                )
            });
        let h = self
            .h_surfaces()
            .iter()
            .enumerate()
            .map(|(index, surface)| {
                Replacement::new(
                    Pattern::from(Atom::from(feynkit_cff::HSurfaceId::from(index))),
                    Pattern::from(surface.to_atom(cut_edges)),
                )
            });
        energy.chain(h).collect()
    }

    fn get_all_replacements_in_lmb(
        &self,
        cut_edges: &[EdgeIndex],
        lmb: &LoopMomentumBasis,
    ) -> Vec<Replacement> {
        let energy = self
            .energy_surfaces()
            .iter()
            .enumerate()
            .map(|(index, surface)| {
                Replacement::new(
                    Pattern::from(Atom::from(feynkit_cff::EnergySurfaceId::from(index))),
                    Pattern::from(surface.to_atom_in_lmb(cut_edges, lmb)),
                )
            });
        let h = self
            .h_surfaces()
            .iter()
            .enumerate()
            .map(|(index, surface)| {
                Replacement::new(
                    Pattern::from(Atom::from(feynkit_cff::HSurfaceId::from(index))),
                    Pattern::from(surface.to_atom_in_lmb(cut_edges, lmb)),
                )
            });
        energy.chain(h).collect()
    }
}
