use clarabel::algebra::*;
use clarabel::solver::*;
use lorentz_vector::LorentzVector;
use num::Complex;
use rayon::vec;
use serde_yaml::Index;

use crate::graph::Graph;
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
    graph: &Graph,
    esurface_ids: &[usize],
    external_momenta: &[LorentzVector<f64>],
) -> DefaultSolver {
    let lmb = &graph.loop_momentum_basis;
    let num_loops = lmb.basis.len();
    let num_loop_vars = 3 * num_loops;

    let esurfaces = &graph
        .derived_data
        .cff_expression
        .as_ref()
        .unwrap()
        .esurfaces;

    let _esurface_derived_data = graph.derived_data.esurface_derived_data.as_ref().unwrap();

    let temp_loop_momenta = vec![LorentzVector::default(); num_loops];
    let energy_cache = graph.compute_onshell_energies(&temp_loop_momenta, external_momenta);

    // first we study the structure of the problem
    let mut propagator_constraints: Vec<PropagatorConstraint> =
        Vec::with_capacity(graph.edges.len());

    let mut inequivalent_masses: Vec<Complex<f64>> = vec![];

    let mut esurface_constraints: Vec<Vec<usize>> = Vec::with_capacity(esurface_ids.len());

    for esurface_id in esurface_ids {
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
                let mass_pointer = graph.edges[propagator_id].particle.mass.value.map(|m| {
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
    for (constraint_index, (esurface_id, esurface_constraint)) in esurface_ids
        .iter()
        .zip(esurface_constraints.iter())
        .enumerate()
    {
        for prop_index in esurface_constraint {
            a_matrix[constraint_index + 1][*prop_index + propagator_index_offset] = 1.0;
        }

        let shift_part = esurfaces[*esurface_id].compute_shift_part(&energy_cache);
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
    graph: &Graph,
    esurface_ids: &[usize],
    external_momenta: &[LorentzVector<f64>],
) -> Option<Vec<LorentzVector<f64>>> {
    let mut solver = construct_solver(graph, esurface_ids, external_momenta);
    solver.solve();

    if solver.solution.status == SolverStatus::Solved
        || solver.solution.status == SolverStatus::AlmostSolved
    {
        let center = extract_center(graph.loop_momentum_basis.basis.len(), &solver.solution.x);
        Some(center)
    } else {
        None
    }
}

pub fn find_maximal_overlap(
    graph: &Graph,
    mut existing_esurface_ids: Vec<usize>,
    external_momenta: &[LorentzVector<f64>],
) -> Vec<(Vec<usize>, Vec<LorentzVector<f64>>)> {
    let mut res = vec![];
    let num_loops = graph.loop_momentum_basis.basis.len();

    let esurface_list = &graph
        .derived_data
        .cff_expression
        .as_ref()
        .unwrap()
        .esurfaces;

    // first try if all esurfaces have a single center
    let option_center = find_center(graph, &existing_esurface_ids, external_momenta);

    if let Some(center) = option_center {
        return vec![(existing_esurface_ids, center)];
    }

    // construct pairs of esurfaces
    let len = existing_esurface_ids.len();
    let mut pair_centers = Vec::with_capacity(len * (len - 1) / 2);
    let mut has_overlap_with = vec![vec![]; len];

    for i in 0..len {
        for j in i + 1..len {
            let center = find_center(
                graph,
                &[existing_esurface_ids[i], existing_esurface_ids[j]],
                external_momenta,
            );

            if center.is_some() {
                has_overlap_with[i].push(j);
                has_overlap_with[j].push(i);
            }

            pair_centers.push(center);
        }
    }

    for (index, esurfaces_paired) in has_overlap_with.iter().enumerate() {
        if esurfaces_paired.is_empty() {
            res.push((
                vec![existing_esurface_ids[index]],
                esurface_list[existing_esurface_ids[index]].get_point_inside(),
            ));

            let index_to_remove = existing_esurface_ids
                .iter()
                .position(|&x| x == index)
                .unwrap();

            existing_esurface_ids.remove(index_to_remove);
        } else if esurfaces_paired.len() == 1 {
            res.push((
                vec![
                    existing_esurface_ids[index],
                    existing_esurface_ids[esurfaces_paired[0]],
                ],
                pair_centers[index + len * esurfaces_paired[0]]
                    .as_ref()
                    .unwrap()
                    .clone(),
            ));
        }
    }

    res
}

struct EsurfacePairs {
    data: Vec<Option<Vec<LorentzVector<f64>>>>,
}
