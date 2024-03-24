use core::panic;

use ahash::HashMap;
use ahash::HashMapExt;
use clarabel::algebra::*;
use clarabel::solver::*;
use lorentz_vector::LorentzVector;
use num::Complex;

use crate::cff::Esurface;
use crate::graph::LoopMomentumBasis;
use crate::utils::compute_shift_part;
use crate::utils::FloatLike;

/// Helper struct to construct the socp problem
struct PropagatorConstraint<'a> {
    propagator_id: usize,            // index into the graph's edges
    mass_pointer: Option<usize>,     // pointer to value of unique mass
    loop_signature: &'a [isize],     // loop signature of the propagator
    external_signature: &'a [isize], // external signature of the propagator
}

impl<'a> PropagatorConstraint<'a> {
    fn get_dimension(&self) -> usize {
        let mass_value = if self.mass_pointer.is_some() { 1 } else { 0 };

        3 * self.loop_signature.iter().filter(|&x| *x != 0).count() + mass_value + 1
    }
}

fn extract_center<T: FloatLike>(num_loops: usize, solution: &[T]) -> Vec<LorentzVector<T>> {
    let len = solution.len();
    let num_loop_vars = 3 * num_loops;

    solution[len - num_loop_vars..]
        .windows(3)
        .map(|window| LorentzVector::from_args(T::zero(), window[0], window[1], window[2]))
        .collect()
}

fn construct_solver(
    lmb: &LoopMomentumBasis,
    existing_esurface_ids: &[usize],
    esurfaces: &[Esurface],
    edge_masses: &[Option<Complex<f64>>],
    external_momenta: &[LorentzVector<f64>],
) -> DefaultSolver {
    let num_loops = lmb.basis.len();
    let num_loop_vars = 3 * num_loops;
    let num_edges = lmb.edge_signatures.len();

    // first we study the structure of the problem
    let mut propagator_constraints: Vec<PropagatorConstraint> = Vec::with_capacity(num_edges);

    let mut inequivalent_masses: Vec<Complex<f64>> = vec![];

    let mut esurface_constraints: Vec<Vec<usize>> = Vec::with_capacity(existing_esurface_ids.len());

    for esurface_id in existing_esurface_ids {
        let esurface = &esurfaces[*esurface_id];
        let mut esurface_constraint_indices: Vec<usize> = Vec::with_capacity(6);

        for &edge_id in &esurface.energies {
            if let Some(edge_position) = propagator_constraints
                .iter()
                .position(|constraint| constraint.propagator_id == edge_id)
            {
                esurface_constraint_indices.push(edge_position);
            } else {
                let propagator_id = edge_id;
                let mass_pointer = edge_masses[propagator_id].map(|m| {
                    if let Some(mass_position) = inequivalent_masses.iter().position(|&x| x == m) {
                        mass_position
                    } else {
                        inequivalent_masses.push(m);
                        inequivalent_masses.len() - 1
                    }
                });

                let (loop_signature, external_signature) = &lmb.edge_signatures[edge_id];

                let propagator_constraint = PropagatorConstraint {
                    propagator_id,
                    mass_pointer,
                    loop_signature,
                    external_signature,
                };

                propagator_constraints.push(propagator_constraint);
                esurface_constraint_indices.push(propagator_constraints.len() - 1);
            };
        }

        esurface_constraints.push(esurface_constraint_indices);
    }

    // now we know the structure, so we can put it in matrix form.
    // variables are stored as [r , x_p_0, ... x_p_m, k1_x, k1_y, ... k_n_x, k_n_y, k_n_z]
    let propagator_index_offset = 1;
    let loop_momentum_offset = propagator_index_offset + propagator_constraints.len();

    let num_primal_variables = loop_momentum_offset + num_loop_vars;

    let cone_dimension_sum = propagator_constraints
        .iter()
        .map(|constraint| constraint.get_dimension())
        .sum::<usize>();

    let num_constaints = 1 + esurface_constraints.len() + cone_dimension_sum;

    // quadratic part of the objective function, we don't need this for now
    let p_matrix: CscMatrix<f64> =
        CscMatrix::spalloc((num_primal_variables, num_primal_variables), 0);

    // objective function
    let mut q_vector = vec![0.0; num_primal_variables];
    q_vector[0] = 1.0;

    // construct the cones
    let mut cones: Vec<SupportedConeT<f64>> = Vec::with_capacity(1 + propagator_constraints.len());
    cones.push(NonnegativeConeT(1));
    cones.push(NonnegativeConeT(esurface_constraints.len()));

    for prop_constraint in &propagator_constraints {
        cones.push(SecondOrderConeT(prop_constraint.get_dimension()));
    }

    // write the constaint equations
    let mut a_matrix = vec![vec![0.0; num_primal_variables]; num_constaints];
    let mut b_vector = vec![0.0; num_constaints];

    a_matrix[0][0] = 1.0;
    // esurface constraints
    for (constraint_index, (esurface_id, esurface_constraint)) in existing_esurface_ids
        .iter()
        .zip(esurface_constraints.iter())
        .enumerate()
    {
        for prop_index in esurface_constraint {
            a_matrix[constraint_index + 1][*prop_index + propagator_index_offset] = 1.0;
        }

        let shift_part =
            esurfaces[*esurface_id].compute_shift_part_from_momenta(lmb, external_momenta);
        b_vector[constraint_index + 1] = -shift_part;
        a_matrix[constraint_index + 1][0] = -1.0;
    }

    // propagator constraints
    let mut vertical_offset = esurface_constraints.len() + 1;
    for (cone_index, propagator_constraint) in propagator_constraints.iter().enumerate() {
        a_matrix[vertical_offset][propagator_index_offset + cone_index] = -1.0;
        vertical_offset += 1;

        let spatial_shift =
            compute_shift_part(propagator_constraint.external_signature, external_momenta);

        b_vector[vertical_offset] = spatial_shift.x;
        b_vector[vertical_offset + 1] = spatial_shift.y;
        b_vector[vertical_offset + 2] = spatial_shift.z;

        for (loop_index, individual_loop_signature) in
            propagator_constraint.loop_signature.iter().enumerate()
        {
            if *individual_loop_signature != 0 {
                a_matrix[vertical_offset][loop_momentum_offset + 3 * loop_index] =
                    -*individual_loop_signature as f64;
                vertical_offset += 1;
                a_matrix[vertical_offset][loop_momentum_offset + 3 * loop_index + 1] =
                    -*individual_loop_signature as f64;
                vertical_offset += 1;
                a_matrix[vertical_offset][loop_momentum_offset + 3 * loop_index + 2] =
                    -*individual_loop_signature as f64;
                vertical_offset += 1;
            }
        }

        if let Some(mass_index) = propagator_constraint.mass_pointer {
            b_vector[vertical_offset] = inequivalent_masses[mass_index].re;
            vertical_offset += 1;
        }
    }

    let a_matrix_sparse = CscMatrix::from(&a_matrix);

    let settings = DefaultSettingsBuilder::default()
        .verbose(false)
        .build()
        .unwrap();

    DefaultSolver::new(
        &p_matrix,
        &q_vector,
        &a_matrix_sparse,
        &b_vector,
        &cones,
        settings,
    )
}

pub fn find_center(
    lmb: &LoopMomentumBasis,
    existing_esurfaces_ids: &[usize],
    esurfaces: &[Esurface],
    edge_masses: &[Option<Complex<f64>>],
    external_momenta: &[LorentzVector<f64>],
) -> Option<Vec<LorentzVector<f64>>> {
    let mut solver = construct_solver(
        lmb,
        existing_esurfaces_ids,
        esurfaces,
        edge_masses,
        external_momenta,
    );
    solver.solve();

    let loop_number = lmb.basis.len();

    if solver.solution.status == SolverStatus::Solved
        || solver.solution.status == SolverStatus::AlmostSolved
    {
        let center = extract_center(loop_number, &solver.solution.x);
        Some(center)
    } else {
        None
    }
}

#[allow(unused_variables)]
pub fn find_maximal_overlap(
    lmb: &LoopMomentumBasis,
    existing_esurface_ids: &[usize],
    esurfaces: &[Esurface],
    edge_masses: &[Option<Complex<f64>>],
    external_momenta: &[LorentzVector<f64>],
) -> Vec<(Vec<usize>, Vec<LorentzVector<f64>>)> {
    let mut res = vec![];
    let num_loops = lmb.basis.len();

    // first try if all esurfaces have a single center
    let option_center = find_center(
        lmb,
        existing_esurface_ids,
        esurfaces,
        edge_masses,
        external_momenta,
    );

    if let Some(center) = option_center {
        res.push((existing_esurface_ids.to_vec(), center));
        return res;
    }

    // construct pairs of esurfaces
    let len = existing_esurface_ids.len();
    let mut pair_centers = EsurfacePairs::new(len);
    let mut has_overlap_with = vec![vec![]; len];

    for i in 0..len {
        for j in i + 1..len {
            let option_center = find_center(
                lmb,
                &[existing_esurface_ids[i], existing_esurface_ids[j]],
                esurfaces,
                edge_masses,
                external_momenta,
            );

            if let Some(center) = option_center {
                has_overlap_with[i].push(j);
                has_overlap_with[j].push(i);
                pair_centers.insert((i, j), center);
            }
        }
    }

    println!("{:?}", has_overlap_with);
    println!("{:?}", pair_centers);

    let existing_esurfaces_for_stage_2 = existing_esurface_ids.to_vec();

    for (esurface_id_index, overlaps_with) in has_overlap_with.iter().enumerate() {
        if overlaps_with.is_empty() {}
    }

    panic!();
}

#[derive(Debug)]
struct EsurfacePairs {
    data: HashMap<(usize, usize), Vec<LorentzVector<f64>>>,
}

impl EsurfacePairs {
    #[allow(dead_code)]
    fn get_remove(&mut self, mut pair: (usize, usize)) -> Option<Vec<LorentzVector<f64>>> {
        if pair.0 > pair.1 {
            std::mem::swap(&mut pair.0, &mut pair.1);
        }

        self.data.remove(&pair)
    }

    fn insert(&mut self, pair: (usize, usize), center: Vec<LorentzVector<f64>>) {
        if pair.0 > pair.1 {
            self.data.insert((pair.1, pair.0), center);
        } else {
            self.data.insert(pair, center);
        }
    }

    fn new(num_esurfaces: usize) -> Self {
        Self {
            data: HashMap::with_capacity(num_esurfaces * (num_esurfaces - 1) / 2),
        }
    }
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;
    use lorentz_vector::LorentzVector;
    use num::Complex;

    use crate::{cff::Esurface, graph::LoopMomentumBasis};

    struct HelperBoxStructure {
        external_momenta: [LorentzVector<f64>; 3],
        lmb: LoopMomentumBasis,
        esurfaces: [Esurface; 4],
        edge_masses: Vec<Option<Complex<f64>>>,
    }

    impl HelperBoxStructure {
        fn new(masses: Option<[f64; 4]>) -> Self {
            let external_momenta = [
                LorentzVector::from_args(14.0, -6.6, -40.0, 0.0),
                LorentzVector::from_args(-43.0, 15.2, 33.0, 0.0),
                LorentzVector::from_args(-17.9, -50.0, 11.8, 0.0),
            ];

            let box_basis = vec![4];
            let box_signatures = vec![
                (vec![0], vec![1, 0, 0]),
                (vec![0], vec![0, 1, 0]),
                (vec![0], vec![0, 0, 1]),
                (vec![0], vec![-1, -1, -1]),
                (vec![1], vec![0, 0, 0]),
                (vec![1], vec![1, 0, 0]),
                (vec![1], vec![1, 1, 0]),
                (vec![1], vec![1, 1, 1]),
            ];

            let box_lmb = LoopMomentumBasis {
                basis: box_basis,
                edge_signatures: box_signatures,
            };

            let esurfaces = [
                Esurface {
                    energies: vec![5, 6],
                    sub_orientation: vec![],
                    shift: vec![1],
                    shift_signature: vec![true],
                },
                Esurface {
                    energies: vec![4, 7],
                    sub_orientation: vec![],
                    shift: vec![1, 2],
                    shift_signature: vec![true, true],
                },
                Esurface {
                    energies: vec![4, 6],
                    sub_orientation: vec![],
                    shift: vec![0, 1],
                    shift_signature: vec![true, true],
                },
                Esurface {
                    energies: vec![4, 7],
                    sub_orientation: vec![],
                    shift: vec![0, 1, 2],
                    shift_signature: vec![true, true, true],
                },
            ];

            let edge_masses = match masses {
                Some(masses) => {
                    let mut edge_masses = vec![None; 4];
                    let mut real_masses = masses
                        .iter()
                        .map(|&x| Some(Complex::new(x, 0.0)))
                        .collect_vec();

                    edge_masses.append(&mut real_masses);
                    edge_masses
                }
                None => vec![None; 8],
            };

            Self {
                external_momenta,
                lmb: box_lmb,
                esurfaces,
                edge_masses,
            }
        }
    }

    #[test]
    fn test_box_4e() {
        // massless variant
        let box4e = HelperBoxStructure::new(None);
    }
}
