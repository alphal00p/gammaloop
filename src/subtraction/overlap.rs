use lorentz_vector::LorentzVector;

use crate::graph::Graph;
use crate::utils::FloatLike;

pub fn get_existing_esurfaces(graph: &Graph, independent_external_momenta: &[usize]) -> Vec<usize> {
    let esurfaces = &graph
        .derived_data
        .cff_expression
        .as_ref()
        .unwrap()
        .esurfaces;

    todo!()
}

struct EsurfaceDerivedData {
    esurface_data: Vec<EsurfaceData>,
}

struct EsurfaceData {
    cut_momentum_basis: usize,
    mirror_esurface: usize, // points to the esurface with opposite orientation
    mass_sum_squared: f64,
    shift_signature: Vec<isize>,
}

impl EsurfaceData {
    fn existence_condition<T: FloatLike>(&self, externals: &[LorentzVector<T>]) -> (T, T) {
        let mut shift = LorentzVector::new();

        for (i, external) in externals.iter().enumerate() {
            match self.shift_signature[i] {
                1 => shift += external,
                -1 => shift -= external,
                0 => {}
                _ => unreachable!("Shift signature must be -1, 0 or 1"),
            }
        }

        (
            shift.t,
            shift.square() - Into::<T>::into(self.mass_sum_squared),
        )
    }
}
