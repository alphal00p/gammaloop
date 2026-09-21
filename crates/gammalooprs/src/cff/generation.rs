//! GammaLoop runtime extensions for FeynKit-owned CFF generation.

use color_eyre::Result;
use feynkit_cff::{CffResult, ShiftRewrite, SurfaceCache};
use linnet::half_edge::involution::EdgeIndex;
use symbolica::{
    atom::{Atom, AtomCore},
    id::Replacement,
};

use crate::{graph::LoopMomentumBasis, settings::global::OrientationPattern};

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
        self.replacements_with(
            |edge| {
                if cut_edges.contains(&edge) {
                    crate::utils::cut_energy(edge)
                } else {
                    crate::utils::ose_atom_from_index(edge)
                }
            },
            crate::utils::external_energy_atom_from_index,
        )
    }

    fn get_all_replacements_in_lmb(
        &self,
        cut_edges: &[EdgeIndex],
        lmb: &LoopMomentumBasis,
    ) -> Vec<Replacement> {
        self.replacements_with(
            |edge| {
                if cut_edges.contains(&edge) {
                    crate::utils::cut_energy(edge)
                } else {
                    crate::utils::ose_atom_from_index(edge)
                }
            },
            |edge| {
                lmb.edge_signatures[edge].external.iter_enumerated().fold(
                    Atom::Zero,
                    |sum, (external_index, sign)| {
                        let energy = crate::utils::external_energy_atom_from_index(
                            lmb.ext_edges[external_index],
                        );
                        match sign {
                            crate::momentum::SignOrZero::Zero => sum,
                            crate::momentum::SignOrZero::Plus => sum + energy,
                            crate::momentum::SignOrZero::Minus => sum - energy,
                        }
                    },
                )
            },
        )
    }
}
