/// A surface whose inverse is zero, useful for setting surface terms to zero.
pub use three_dimensional_reps::surface::InfiniteSurface;
/// A surface equal to one, useful for representing a single vertex.
pub use three_dimensional_reps::surface::UnitSurface;
pub use three_dimensional_reps::surface::{
    EsurfaceCollection, EsurfaceID, HsurfaceCollection, HsurfaceID, HybridSurfaceID,
    LinearEnergyExpr, LinearSurface, LinearSurfaceCollection, LinearSurfaceID, LinearSurfaceKind,
    SurfaceAtom, SurfaceOrigin,
};

use linnet::{half_edge::involution::EdgeIndex, num_traits::SignOrZero};
use symbolica::{
    atom::{Atom, AtomCore},
    id::{Pattern, Replacement},
};

use super::{esurface::Esurface, hsurface::Hsurface};
use crate::{
    graph::LoopMomentumBasis,
    utils::{GS, cut_energy, external_energy_atom_from_index, ose_atom_from_index},
};

pub type HybridSurface = three_dimensional_reps::surface::HybridSurface<Esurface, Hsurface>;
pub type HybridSurfaceRef<'a> =
    three_dimensional_reps::surface::HybridSurfaceRef<'a, Esurface, Hsurface>;
pub type HybridSurfaceCollection =
    three_dimensional_reps::surface::HybridSurfaceCollection<Esurface, Hsurface>;
pub type HybridSurfaceCache<T> = three_dimensional_reps::surface::HybridSurfaceCache<T>;
pub type SurfaceCache = three_dimensional_reps::surface::SurfaceCache<Esurface, Hsurface>;

impl SurfaceAtom for Esurface {
    fn to_atom(&self, cut_edges: &[EdgeIndex]) -> Atom {
        self.to_atom(cut_edges)
    }
}

impl SurfaceAtom for Hsurface {
    fn to_atom(&self, cut_edges: &[EdgeIndex]) -> Atom {
        self.to_atom(cut_edges)
    }
}

pub trait GammaLoopLinearEnergyExpr {
    fn to_atom_gs(&self, cut_edges: &[EdgeIndex]) -> Atom;
}

impl GammaLoopLinearEnergyExpr for LinearEnergyExpr {
    fn to_atom_gs(&self, cut_edges: &[EdgeIndex]) -> Atom {
        let internal = self
            .internal_terms
            .iter()
            .fold(Atom::new(), |acc, (edge_id, coeff)| {
                let energy = if cut_edges.contains(edge_id) {
                    cut_energy(*edge_id)
                } else {
                    ose_atom_from_index(*edge_id)
                };
                acc + coeff.clone() * energy
            });

        let external = self
            .external_terms
            .iter()
            .fold(Atom::new(), |acc, (edge_id, coeff)| {
                acc + coeff.clone() * external_energy_atom_from_index(*edge_id)
            });

        let scale = if self.uniform_scale_coeff.is_zero() {
            Atom::new()
        } else {
            self.uniform_scale_coeff.clone() * Atom::var(GS.numerator_sampling_scale)
        };

        internal + external + scale + self.constant.clone()
    }
}

pub trait GammaLoopHybridSurfaceRef {
    fn to_atom_gs(&self, cut_edges: &[EdgeIndex]) -> Atom;
}

impl<E: SurfaceAtom, H: SurfaceAtom> GammaLoopHybridSurfaceRef
    for three_dimensional_reps::surface::HybridSurfaceRef<'_, E, H>
{
    fn to_atom_gs(&self, cut_edges: &[EdgeIndex]) -> Atom {
        match self {
            three_dimensional_reps::surface::HybridSurfaceRef::Esurface(surface) => {
                surface.to_atom(cut_edges)
            }
            three_dimensional_reps::surface::HybridSurfaceRef::Hsurface(surface) => {
                surface.to_atom(cut_edges)
            }
            three_dimensional_reps::surface::HybridSurfaceRef::Linear(surface) => {
                surface.expression.to_atom_gs(cut_edges)
            }
            three_dimensional_reps::surface::HybridSurfaceRef::Unit(_) => Atom::num(1),
            three_dimensional_reps::surface::HybridSurfaceRef::Infinite(_) => {
                symbolica::parse!("η_inf")
            }
        }
    }
}

pub trait GammaLoopSurfaceCache {
    fn get_all_replacements_gs(&self, cut_edges: &[EdgeIndex]) -> Vec<Replacement>;
    fn get_all_replacements_gs_in_lmb(
        &self,
        cut_edges: &[EdgeIndex],
        lmb: &LoopMomentumBasis,
    ) -> Vec<Replacement>;
}

impl<E: SurfaceAtom, H: SurfaceAtom> GammaLoopSurfaceCache
    for three_dimensional_reps::surface::SurfaceCache<E, H>
{
    fn get_all_replacements_gs(&self, cut_edges: &[EdgeIndex]) -> Vec<Replacement> {
        self.iter_all_surfaces()
            .map(|(id, surface)| {
                let id_atom = Pattern::from(Atom::from(id));
                let surface_atom = Pattern::from(surface.to_atom_gs(cut_edges));
                Replacement::new(id_atom, surface_atom)
            })
            .collect()
    }

    fn get_all_replacements_gs_in_lmb(
        &self,
        cut_edges: &[EdgeIndex],
        lmb: &LoopMomentumBasis,
    ) -> Vec<Replacement> {
        let external_shift_replacements = lmb
            .edge_signatures
            .iter()
            .map(|(edge_id, signature)| {
                let canonical = signature.external.iter_enumerated().fold(
                    Atom::Zero,
                    |sum, (external_index, sign)| {
                        let energy = external_energy_atom_from_index(lmb.ext_edges[external_index]);
                        match sign {
                            SignOrZero::Zero => sum,
                            SignOrZero::Plus => sum + energy,
                            SignOrZero::Minus => sum - energy,
                        }
                    },
                );
                Replacement::new(
                    external_energy_atom_from_index(edge_id).to_pattern(),
                    canonical,
                )
            })
            .collect::<Vec<_>>();

        self.iter_all_surfaces()
            .map(|(id, surface)| {
                let id_atom = Pattern::from(Atom::from(id));
                let surface_atom = Pattern::from(
                    surface
                        .to_atom_gs(cut_edges)
                        .replace_multiple(&external_shift_replacements),
                );
                Replacement::new(id_atom, surface_atom)
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use linnet::half_edge::{
        involution::EdgeIndex,
        subgraph::{SuBitGraph, SubSetLike},
    };
    use symbolica::atom::{Atom, AtomCore};

    use super::{EsurfaceID, GammaLoopSurfaceCache, SurfaceCache};
    use crate::{
        cff::{VertexSet, esurface::Esurface},
        graph::LoopMomentumBasis,
        momentum::signature::LoopExtSignature,
        utils::{external_energy_atom_from_index, test_utils::dummy_hedge_graph},
    };

    #[test]
    fn replacements_in_lmb_use_canonical_external_edges_not_carrier_edges() {
        let dummy_graph = dummy_hedge_graph(9);
        let mut edge_signatures = dummy_graph
            .new_edgevec_from_iter(
                (0..9).map(|_| LoopExtSignature::from((Vec::<isize>::new(), vec![0, 0]))),
            )
            .unwrap();
        edge_signatures[EdgeIndex::from(8)] =
            LoopExtSignature::from((Vec::<isize>::new(), vec![0, 1]));
        let lmb = LoopMomentumBasis {
            tree: SuBitGraph::empty(0),
            loop_edges: vec![].into(),
            ext_edges: vec![EdgeIndex::from(2), EdgeIndex::from(6)].into(),
            edge_signatures,
        };
        let mut surfaces = SurfaceCache::new();
        surfaces.esurface_cache.push(Esurface {
            energies: vec![],
            external_shift: vec![(EdgeIndex::from(8), -1)],
            vertex_set: VertexSet::dummy(),
        });

        let atom = Atom::from(EsurfaceID::from(0))
            .replace_multiple(surfaces.get_all_replacements_gs_in_lmb(&[], &lmb))
            .expand();
        let expected =
            (Atom::num(-1) * external_energy_atom_from_index(EdgeIndex::from(6))).expand();

        assert_eq!(atom.to_canonical_string(), expected.to_canonical_string());
    }
}
