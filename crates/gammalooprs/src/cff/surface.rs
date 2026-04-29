pub use three_dimensional_reps::surface::{
    EsurfaceCollection, EsurfaceID, HsurfaceCollection, HsurfaceID, HybridSurfaceID,
    InfiniteSurface, LinearEnergyExpr, LinearSurface, LinearSurfaceCollection, LinearSurfaceID,
    LinearSurfaceKind, RationalCoefficient, SurfaceAtom, SurfaceOrigin, UnitSurface,
};

use linnet::half_edge::involution::EdgeIndex;
use symbolica::{
    atom::{Atom, AtomCore},
    id::{Pattern, Replacement},
};

use super::{esurface::Esurface, hsurface::Hsurface};
use crate::utils::{GS, cut_energy, external_energy_atom_from_index, ose_atom_from_index};

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
                acc + coeff.to_atom() * energy
            });

        let external = self
            .external_terms
            .iter()
            .fold(Atom::new(), |acc, (edge_id, coeff)| {
                acc + coeff.to_atom() * external_energy_atom_from_index(*edge_id)
            });

        let scale = if self.uniform_scale_coeff.is_zero() {
            Atom::new()
        } else {
            self.uniform_scale_coeff.to_atom() * Atom::var(GS.numerator_sampling_scale)
        };

        internal + external + scale + self.constant.to_atom()
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
    fn substitute_energies_gs(&self, atom: &Atom, cut_edges: &[EdgeIndex]) -> Atom;
    fn get_all_replacements_gs(&self, cut_edges: &[EdgeIndex]) -> Vec<Replacement>;
}

impl<E: SurfaceAtom, H: SurfaceAtom> GammaLoopSurfaceCache
    for three_dimensional_reps::surface::SurfaceCache<E, H>
{
    fn substitute_energies_gs(&self, atom: &Atom, cut_edges: &[EdgeIndex]) -> Atom {
        let replacement_rules = self.get_all_replacements_gs(cut_edges);
        atom.replace_multiple(&replacement_rules)
    }

    fn get_all_replacements_gs(&self, cut_edges: &[EdgeIndex]) -> Vec<Replacement> {
        self.iter_all_surfaces()
            .map(|(id, surface)| {
                let id_atom = Pattern::from(Atom::from(id));
                let surface_atom = Pattern::from(surface.to_atom_gs(cut_edges));
                Replacement::new(id_atom, surface_atom)
            })
            .collect()
    }
}
