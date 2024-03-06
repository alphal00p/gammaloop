use lorentz_vector::LorentzVector;
use num::Float;
use serde::{Deserialize, Serialize};

use crate::graph::Graph;
use crate::utils::FloatLike;

pub fn get_existing_esurfaces<T: FloatLike>(
    graph: &Graph,
    esurface_derived_data: &EsurfaceDerivedData,
    externals: &[LorentzVector<T>],
) -> Vec<usize> {
    let esurfaces = &graph
        .derived_data
        .cff_expression
        .as_ref()
        .unwrap()
        .esurfaces;

    todo!()
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct EsurfaceDerivedData {
    esurface_data: Vec<EsurfaceData>,
}

#[derive(Clone, Serialize, Deserialize, Debug)]
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

pub fn generate_esurface_data(graph: &Graph) -> EsurfaceDerivedData {
    let lmbs = graph.derived_data.loop_momentum_bases.as_ref().unwrap();
    let cff = graph.derived_data.cff_expression.as_ref().unwrap();

    for esurface in cff.esurfaces.iter() {
        // find the cut momentum basis
        let energies = &esurface.energies;

        // find the lmb that has at least n-1 energies of the esurface
        let cut_momentum_basis = lmbs
            .iter()
            .position(|lmb| {
                lmb.basis.iter().filter(|&i| energies.contains(&i)).count() >= energies.len() - 1
            })
            .unwrap();

        todo!()
    }

    todo!()
}
