use ahash::HashMap;
use ahash::HashMapExt;
use ahash::HashSet;
use clarabel::algebra::*;
use clarabel::solver::*;
use core::panic;
use itertools::Itertools;
use spenso::complex::Complex;

use crate::cff::esurface::EsurfaceCollection;
use crate::cff::esurface::ExistingEsurfaceId;
use crate::cff::esurface::ExistingEsurfaces;
use crate::graph::LoopMomentumBasis;
use crate::momentum::FourMomentum;
use crate::momentum::ThreeMomentum;
use crate::utils::compute_shift_part;
use crate::utils::F;

#[derive(Debug, Clone)]
pub struct OverlapGroup {
    pub existing_esurfaces: Vec<ExistingEsurfaceId>,
    pub center: Vec<ThreeMomentum<F<f64>>>,
}

#[derive(Debug, Clone)]
pub struct OverlapStructure {
    pub overlap_groups: Vec<OverlapGroup>,
}

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

        3 + mass_value + 1
    }
}

fn extract_center(num_loops: usize, solution: &[f64]) -> Vec<ThreeMomentum<F<f64>>> {
    let len = solution.len();
    let num_loop_vars = 3 * num_loops;

    solution[len - num_loop_vars..]
        .chunks(3)
        .map(|window| {
            ThreeMomentum::new(
                F::from_f64(window[0]),
                F::from_f64(window[1]),
                F::from_f64(window[2]),
            )
        })
        .collect()
}

fn construct_solver(
    lmb: &LoopMomentumBasis,
    esurfaces_to_consider: &[ExistingEsurfaceId],
    existing_esurfaces: &ExistingEsurfaces,
    esurfaces: &EsurfaceCollection,
    edge_masses: &[Option<Complex<F<f64>>>],
    external_momenta: &[FourMomentum<F<f64>>],
    verbose: bool,
) -> DefaultSolver {
    let num_loops = lmb.basis.len();
    let num_loop_vars = 3 * num_loops;
    let num_edges = lmb.edge_signatures.len();

    // first we study the structure of the problem
    let mut propagator_constraints: Vec<PropagatorConstraint> = Vec::with_capacity(num_edges);

    let mut inequivalent_masses: Vec<Complex<F<f64>>> = vec![];

    let mut esurface_constraints: Vec<Vec<usize>> = Vec::with_capacity(esurfaces_to_consider.len());

    for esurface_id in esurfaces_to_consider.iter() {
        let esurface = &esurfaces[existing_esurfaces[*esurface_id]];
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

    // perhaps it is faster to build the sparse matrix directly, but it is also extremely unreadable
    let mut a_matrix = vec![vec![0.0; num_primal_variables]; num_constaints];
    let mut b_vector = vec![0.0; num_constaints];

    a_matrix[0][0] = 1.0;
    // esurface constraints
    for (constraint_index, (esurface_id, esurface_constraint)) in esurfaces_to_consider
        .iter()
        .zip(esurface_constraints.iter())
        .enumerate()
    {
        for prop_index in esurface_constraint {
            a_matrix[constraint_index + 1][*prop_index + propagator_index_offset] = 1.0;
        }

        let shift_part = esurfaces[existing_esurfaces[*esurface_id]]
            .compute_shift_part_from_momenta(lmb, external_momenta);
        b_vector[constraint_index + 1] = -shift_part.0;
        a_matrix[constraint_index + 1][0] = -1.0;
    }

    // propagator constraints
    let mut vertical_offset = esurface_constraints.len() + 1;
    for (cone_index, propagator_constraint) in propagator_constraints.iter().enumerate() {
        a_matrix[vertical_offset][propagator_index_offset + cone_index] = -1.0;
        vertical_offset += 1;

        let spatial_shift =
            compute_shift_part(propagator_constraint.external_signature, external_momenta);

        b_vector[vertical_offset] = spatial_shift.spatial.px.0;
        b_vector[vertical_offset + 1] = spatial_shift.spatial.py.0;
        b_vector[vertical_offset + 2] = spatial_shift.spatial.pz.0;

        for (loop_index, individual_loop_signature) in
            propagator_constraint.loop_signature.iter().enumerate()
        {
            if *individual_loop_signature != 0 {
                a_matrix[vertical_offset][loop_momentum_offset + 3 * loop_index] =
                    -*individual_loop_signature as f64;
                a_matrix[vertical_offset + 1][loop_momentum_offset + 3 * loop_index + 1] =
                    -*individual_loop_signature as f64;
                a_matrix[vertical_offset + 2][loop_momentum_offset + 3 * loop_index + 2] =
                    -*individual_loop_signature as f64;
            }
        }

        vertical_offset += 3;

        if let Some(mass_index) = propagator_constraint.mass_pointer {
            b_vector[vertical_offset] = inequivalent_masses[mass_index].re.0;
            vertical_offset += 1;
        }
    }

    if verbose {
        println!("problem");
        for (row_a, row_b) in a_matrix.iter().zip(b_vector.iter()) {
            println!("{:?} - {:?}", row_a, row_b);
        }
    }

    let a_matrix_sparse = CscMatrix::from(&a_matrix);

    let settings = DefaultSettingsBuilder::default()
        .verbose(verbose)
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
    esurfaces_to_consider: &[ExistingEsurfaceId],
    existing_esurfaces: &ExistingEsurfaces,
    esurfaces: &EsurfaceCollection,
    edge_masses: &[Option<Complex<F<f64>>>],
    external_momenta: &[FourMomentum<F<f64>>],
    verbose: bool,
) -> Option<Vec<ThreeMomentum<F<f64>>>> {
    let mut solver = construct_solver(
        lmb,
        esurfaces_to_consider,
        existing_esurfaces,
        esurfaces,
        edge_masses,
        external_momenta,
        verbose,
    );

    solver.solve();

    let loop_number = lmb.basis.len();

    if solver.solution.status == SolverStatus::Solved {
        let center = extract_center(loop_number, &solver.solution.x);
        Some(center)
    } else if solver.solution.status == SolverStatus::AlmostSolved
        || solver.solution.status == SolverStatus::InsufficientProgress
    {
        // if the solver did not converge, we check if the solution is still valid
        let center = extract_center(loop_number, &solver.solution.x);

        let is_valid = esurfaces_to_consider.iter().all(|&existing_esurface_id| {
            let esurface = &esurfaces[existing_esurfaces[existing_esurface_id]];
            esurface.compute_from_momenta(
                lmb,
                &to_real_mass_vector(edge_masses),
                &center,
                external_momenta,
            ) < F::from_f64(0.0)
        });

        if is_valid {
            Some(center)
        } else {
            None
        }
    } else {
        if verbose {
            println!("{:?}", solver.solution.x);
        }

        None
    }
}

pub fn find_maximal_overlap(
    lmb: &LoopMomentumBasis,
    existing_esurfaces: &ExistingEsurfaces,
    esurfaces: &EsurfaceCollection,
    edge_masses: &[Option<Complex<F<f64>>>],
    external_momenta: &[FourMomentum<F<f64>>],
    debug: usize,
) -> OverlapStructure {
    let mut res = OverlapStructure {
        overlap_groups: vec![],
    };

    let all_existing_esurfaces = existing_esurfaces
        .iter_enumerated()
        .map(|a| a.0)
        .collect_vec();

    // first try if all esurfaces have a single center, we explitely seach a center instead of trying the
    // origin. This is because the origin might not be optimal.
    let option_center = find_center(
        lmb,
        &all_existing_esurfaces,
        existing_esurfaces,
        esurfaces,
        edge_masses,
        external_momenta,
        false,
    );

    if let Some(center) = option_center {
        let single_group = OverlapGroup {
            existing_esurfaces: all_existing_esurfaces,
            center,
        };
        res.overlap_groups.push(single_group);
        return res;
    }

    // if the center is not valid, create a table of all pairs
    let pair_benchmark_start = std::time::Instant::now();
    let esurface_pairs = EsurfacePairs::new(
        lmb,
        existing_esurfaces,
        esurfaces,
        edge_masses,
        external_momenta,
    );
    let pair_benchmark = pair_benchmark_start.elapsed();

    if debug > 3 {
        println!(
            "Time to generate {} pairs: {}",
            esurface_pairs.data.len(),
            pair_benchmark.as_millis(),
        );
    }

    let mut num_disconnected_surfaces = 0;

    for (existing_esurface_id, &esurface_id) in existing_esurfaces.iter_enumerated() {
        // if an esurface overlaps with no other esurface, it is part of the maximal overlap structure
        if esurface_pairs.has_no_overlap(existing_esurface_id) {
            let center = find_center(
                lmb,
                &[existing_esurface_id],
                existing_esurfaces,
                esurfaces,
                edge_masses,
                external_momenta,
                false,
            );

            res.overlap_groups.push(OverlapGroup {
                existing_esurfaces: vec![existing_esurface_id],
                center: center.unwrap_or_else(|| {
                    let esurface = &esurfaces[esurface_id];

                    let mut error_message = String::new();

                    error_message.push_str(&format!(
                        "Could not find center of esuface {:?}\n",
                        esurface_id
                    ));
                    error_message.push_str(&format!("edges: {:?}\n", esurface.energies));

                    error_message
                        .push_str(&format!("external shift: {:?}\n", esurface.external_shift));

                    error_message.push_str(&format!(
                        "Equation in LMB: {}\n",
                        esurface.string_format_in_lmb(lmb)
                    ));

                    error_message.push_str(&format!("External momenta: {:#?}\n", external_momenta));

                    panic!("{}", error_message)
                }),
            });
            num_disconnected_surfaces += 1;
        }
    }

    if debug > 3 {
        println!("num disconnected surfaces: {}", num_disconnected_surfaces);
    }

    if num_disconnected_surfaces == existing_esurfaces.len() {
        return res;
    }

    let mut subset_size =
        if let Some(size) = esurface_pairs.has_pair_with.iter().map(Vec::len).max() {
            size + 1
        } else {
            1
        };

    while subset_size > 1 {
        if debug > 3 {
            println!("current subset_size: {}", subset_size);
        }

        let possible_subsets =
            esurface_pairs.construct_possible_subsets_of_len(existing_esurfaces, subset_size, &res);

        if debug > 3 {
            println!("possible subsets at this size: {}", possible_subsets.len());
        }

        for subset in possible_subsets.iter() {
            let option_center = find_center(
                lmb,
                subset,
                existing_esurfaces,
                esurfaces,
                edge_masses,
                external_momenta,
                false,
            );

            if let Some(center) = option_center {
                res.overlap_groups.push(OverlapGroup {
                    existing_esurfaces: subset.clone(),
                    center,
                });
            }
        }

        subset_size -= 1;

        if debug > 3 {
            println!("current result: {:?}", res);
        }
    }

    res
}

fn is_subset_of_result(subset: &[ExistingEsurfaceId], result: &OverlapStructure) -> bool {
    result.overlap_groups.iter().any(|group| {
        subset
            .iter()
            .all(|&x| group.existing_esurfaces.contains(&x))
    })
}

#[derive(Debug)]
struct EsurfacePairs {
    data: HashMap<(ExistingEsurfaceId, ExistingEsurfaceId), Vec<ThreeMomentum<F<f64>>>>,
    has_pair_with: Vec<Vec<ExistingEsurfaceId>>,
}

impl EsurfacePairs {
    fn insert(
        &mut self,
        pair: (ExistingEsurfaceId, ExistingEsurfaceId),
        center: Vec<ThreeMomentum<F<f64>>>,
    ) {
        if pair.0 > pair.1 {
            self.data.insert((pair.1, pair.0), center);
        } else {
            self.data.insert(pair, center);
        }
    }

    fn pair_exists(&self, pair: (ExistingEsurfaceId, ExistingEsurfaceId)) -> bool {
        if pair.0 > pair.1 {
            self.data.contains_key(&(pair.1, pair.0))
        } else {
            self.data.contains_key(&pair)
        }
    }

    fn new_empty(num_existing_esurfaces: usize) -> Self {
        Self {
            data: HashMap::with_capacity(num_existing_esurfaces * (num_existing_esurfaces - 1) / 2),
            has_pair_with: vec![Vec::with_capacity(num_existing_esurfaces); num_existing_esurfaces],
        }
    }

    fn new(
        lmb: &LoopMomentumBasis,
        existing_esurfaces: &ExistingEsurfaces,
        esurfaces: &EsurfaceCollection,
        edge_masses: &[Option<Complex<F<f64>>>],
        external_momenta: &[FourMomentum<F<f64>>],
    ) -> Self {
        let mut res = Self::new_empty(existing_esurfaces.len());

        let all_existing_esurfaces = existing_esurfaces
            .iter_enumerated()
            .map(|a| a.0)
            .collect_vec();

        for (i, &esurface_id_1) in all_existing_esurfaces.iter().enumerate() {
            for &esurface_id_2 in all_existing_esurfaces.iter().skip(i + 1) {
                let center = find_center(
                    lmb,
                    &[esurface_id_1, esurface_id_2],
                    existing_esurfaces,
                    esurfaces,
                    edge_masses,
                    external_momenta,
                    false,
                );

                if let Some(center) = center {
                    res.insert((esurface_id_1, esurface_id_2), center);
                    res.has_pair_with[Into::<usize>::into(esurface_id_1)].push(esurface_id_2);
                    res.has_pair_with[Into::<usize>::into(esurface_id_2)].push(esurface_id_1);
                }
            }
        }

        res
    }

    fn has_no_overlap(&self, esurface_id: ExistingEsurfaceId) -> bool {
        self.has_pair_with[Into::<usize>::into(esurface_id)].is_empty()
    }

    fn construct_possible_subsets_of_len(
        &self,
        existing_esurfaces: &ExistingEsurfaces,
        subset_len: usize,
        result: &OverlapStructure,
    ) -> HashSet<Vec<ExistingEsurfaceId>> {
        let mut res = HashSet::default();
        let existing_esurfaces_not_in_overlap = existing_esurfaces
            .iter_enumerated()
            .map(|a| a.0)
            .filter(|&existing_esurface_id| {
                !result
                    .overlap_groups
                    .iter()
                    .any(|group| group.existing_esurfaces.contains(&existing_esurface_id))
            });

        let mut possible_options_from_esurfaces_not_in_overlap = HashSet::default();

        for esurface in existing_esurfaces_not_in_overlap.filter(|existing_esurface_id| {
            self.has_pair_with[Into::<usize>::into(*existing_esurface_id)].len() >= subset_len - 1
        }) {
            for possible_combination in self.has_pair_with[Into::<usize>::into(esurface)]
                .iter()
                .combinations(subset_len - 1)
            {
                let mut option = vec![esurface];
                option.extend(possible_combination.iter().copied());
                let mut is_valid = true;

                'pair_loop: for i in 0..subset_len - 1 {
                    for j in i + 1..subset_len - 1 {
                        let pair = (
                            Into::<ExistingEsurfaceId>::into(*possible_combination[i]),
                            Into::<ExistingEsurfaceId>::into(*possible_combination[j]),
                        );

                        if !self.pair_exists(pair) {
                            is_valid = false;
                            break 'pair_loop;
                        }
                    }
                }

                if is_valid {
                    option.sort_unstable();
                    possible_options_from_esurfaces_not_in_overlap.insert(option);
                }
            }
        }

        let existing_pairs_not_in_overlap = self.data.keys().filter(|(left, right)| {
            !is_subset_of_result(&[*left, *right], result)
                && is_subset_of_result(&[*left], result)
                && is_subset_of_result(&[*right], result)
                && self.has_pair_with[Into::<usize>::into(*left)].len() >= subset_len - 1
                && self.has_pair_with[Into::<usize>::into(*right)].len() >= subset_len - 1
        });

        let mut possible_options_from_pairs_not_in_overlap = HashSet::default();

        for pair in existing_pairs_not_in_overlap {
            let possible_additions = existing_esurfaces
                .iter_enumerated()
                .map(|a| a.0)
                .filter(|&id| {
                    id != pair.0
                        && id != pair.1
                        && self.has_pair_with[Into::<usize>::into(pair.0)].contains(&id)
                        && self.has_pair_with[Into::<usize>::into(pair.1)].contains(&id)
                })
                .combinations(subset_len - 2);

            for possible_addition in possible_additions {
                let mut option = vec![pair.0, pair.1];
                option.extend(possible_addition.iter().copied());
                option.sort_unstable();
                possible_options_from_pairs_not_in_overlap.insert(option);
            }
        }

        res.extend(possible_options_from_esurfaces_not_in_overlap);
        res.extend(possible_options_from_pairs_not_in_overlap);

        res
    }
}

fn to_real_mass_vector(edge_masses: &[Option<Complex<F<f64>>>]) -> Vec<F<f64>> {
    edge_masses
        .iter()
        .map(|mass| match mass {
            Some(m) => m.re,
            None => F::from_f64(0.0),
        })
        .collect_vec()
}

#[cfg(test)]
#[allow(dead_code, unused_variables)]
mod tests {
    use super::*;
    use itertools::Itertools;
    use spenso::complex::Complex;

    use crate::{
        cff::{
            cff_graph::VertexSet,
            esurface::{Esurface, EsurfaceID},
        },
        graph::LoopMomentumBasis,
    };

    struct HelperBoxStructure {
        external_momenta: [FourMomentum<F<f64>>; 3],
        lmb: LoopMomentumBasis,
        esurfaces: EsurfaceCollection,
        existing_esurfaces: ExistingEsurfaces,
        edge_masses: Vec<Option<Complex<F<f64>>>>,
    }

    struct HelperBananaStructure {
        external_momenta: [FourMomentum<F<f64>>; 1],
        lmb: LoopMomentumBasis,
        esurfaces: EsurfaceCollection,
        existing_esurfaces: ExistingEsurfaces,
        edge_masses: Vec<Option<Complex<F<f64>>>>,
    }

    impl HelperBoxStructure {
        fn new(masses: Option<[F<f64>; 4]>) -> Self {
            let external_momenta = [
                FourMomentum::from_args(F(14.0), F(-6.6), F(-40.0), F(0.0)),
                FourMomentum::from_args(F(-43.0), F(15.2), F(33.0), F(0.0)),
                FourMomentum::from_args(F(-17.9), F(-50.0), F(11.8), F(0.0)),
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

            let circled_vertices = VertexSet::dummy();

            let esurfaces_array = [
                Esurface {
                    energies: vec![5, 6],
                    sub_orientation: vec![],
                    external_shift: vec![(1, 1)],
                    circled_vertices,
                },
                Esurface {
                    energies: vec![5, 7],
                    sub_orientation: vec![],
                    external_shift: vec![(1, 1), (2, 1)],
                    circled_vertices,
                },
                Esurface {
                    energies: vec![4, 6],
                    sub_orientation: vec![],
                    external_shift: vec![(0, 1), (1, 1)],
                    circled_vertices,
                },
                Esurface {
                    energies: vec![4, 7],
                    sub_orientation: vec![],
                    external_shift: vec![(0, 1), (1, 1), (2, 1)],
                    circled_vertices,
                },
            ];

            let esurfaces = esurfaces_array.to_vec().into();

            let edge_masses = match masses {
                Some(masses) => {
                    let mut edge_masses = vec![None; 4];
                    let mut real_masses = masses
                        .iter()
                        .map(|&x| Some(Complex::new(x, F(0.0))))
                        .collect_vec();

                    edge_masses.append(&mut real_masses);
                    edge_masses
                }
                None => vec![None; 8],
            };

            let existing_esurfaces = (0..4).map(Into::<EsurfaceID>::into).collect_vec().into();

            Self {
                external_momenta,
                lmb: box_lmb,
                existing_esurfaces,
                esurfaces,
                edge_masses,
            }
        }
    }

    impl HelperBananaStructure {
        fn new() -> Self {
            let external_momenta = [FourMomentum::from_args(
                F(10.0),
                F(-10.00000000),
                F(0.0),
                F(0.0),
            )];
            let banana_basis = vec![2, 3];
            let banana_edge_sigs = vec![
                (vec![0, 0], vec![1]),
                (vec![0, 0], vec![-1]),
                (vec![1, 0], vec![0]),
                (vec![0, 1], vec![0]),
                (vec![1, 1], vec![-1]),
            ];

            let banana_lmb = LoopMomentumBasis {
                basis: banana_basis,
                edge_signatures: banana_edge_sigs,
            };

            let only_esurface = Esurface {
                energies: vec![2, 3, 4],
                sub_orientation: vec![],
                external_shift: vec![(0, -1)],
                circled_vertices: VertexSet::dummy(),
            };

            let esurfaces = vec![only_esurface].into();

            let existing_esurfaces = vec![Into::<EsurfaceID>::into(0)].into();
            let edge_masses = vec![None; 5];

            Self {
                external_momenta,
                lmb: banana_lmb,
                esurfaces,
                existing_esurfaces,
                edge_masses,
            }
        }
    }

    #[test]
    fn test_is_subset_of_result() {
        let fake_res = vec![
            (
                vec![
                    Into::<ExistingEsurfaceId>::into(1),
                    Into::<ExistingEsurfaceId>::into(2),
                    Into::<ExistingEsurfaceId>::into(3),
                ],
                vec![],
            ),
            (
                vec![
                    Into::<ExistingEsurfaceId>::into(1),
                    Into::<ExistingEsurfaceId>::into(2),
                    Into::<ExistingEsurfaceId>::into(4),
                ],
                vec![],
            ),
            (
                vec![
                    Into::<ExistingEsurfaceId>::into(2),
                    Into::<ExistingEsurfaceId>::into(3),
                    Into::<ExistingEsurfaceId>::into(4),
                ],
                vec![],
            ),
        ];

        let fake_res = OverlapStructure {
            overlap_groups: fake_res
                .into_iter()
                .map(|(group, center)| OverlapGroup {
                    existing_esurfaces: group,
                    center,
                })
                .collect_vec(),
        };

        let fake_subset = vec![
            Into::<ExistingEsurfaceId>::into(1),
            Into::<ExistingEsurfaceId>::into(2),
        ];

        assert!(is_subset_of_result(&fake_subset, &fake_res));

        let fake_subset_2 = vec![
            Into::<ExistingEsurfaceId>::into(0),
            Into::<ExistingEsurfaceId>::into(4),
        ];

        assert!(!is_subset_of_result(&fake_subset_2, &fake_res));
    }

    #[test]
    fn test_pair_creator() {
        let box4e = HelperBoxStructure::new(None);

        let esurface_pairs = EsurfacePairs::new(
            &box4e.lmb,
            &box4e.existing_esurfaces,
            &box4e.esurfaces,
            &box4e.edge_masses,
            &box4e.external_momenta,
        );

        assert_eq!(esurface_pairs.data.len(), 4);

        let box4e_massive = HelperBoxStructure::new(Some([F(10.5); 4]));
        let esurface_pairs_massive = EsurfacePairs::new(
            &box4e_massive.lmb,
            &box4e_massive.existing_esurfaces,
            &box4e_massive.esurfaces,
            &box4e_massive.edge_masses,
            &box4e_massive.external_momenta,
        );

        assert_eq!(esurface_pairs_massive.data.len(), 0);
    }

    #[test]
    fn test_subset_generator() {
        let box4e = HelperBoxStructure::new(None);

        let esurface_pairs = EsurfacePairs::new(
            &box4e.lmb,
            &box4e.existing_esurfaces,
            &box4e.esurfaces,
            &box4e.edge_masses,
            &box4e.external_momenta,
        );

        let res = OverlapStructure {
            overlap_groups: vec![],
        };
        let subsets_3 =
            esurface_pairs.construct_possible_subsets_of_len(&box4e.existing_esurfaces, 3, &res);

        assert_eq!(subsets_3.len(), 0);

        let subsets_2 =
            esurface_pairs.construct_possible_subsets_of_len(&box4e.existing_esurfaces, 2, &res);
        assert_eq!(subsets_2.len(), 4);
    }

    #[test]
    fn test_box_4e() {
        // massless variant
        let box4e = HelperBoxStructure::new(None);

        let maximal_overlap = find_maximal_overlap(
            &box4e.lmb,
            &box4e.existing_esurfaces,
            &box4e.esurfaces,
            &box4e.edge_masses,
            &box4e.external_momenta,
            0,
        );

        assert_eq!(maximal_overlap.overlap_groups.len(), 4);

        for overlap_group in maximal_overlap.overlap_groups.iter() {
            let esurfaces = &overlap_group.existing_esurfaces;
            let center = &overlap_group.center;

            assert_eq!(esurfaces.len(), 2);

            for esurface in esurfaces.iter() {
                let esurfaec_val = box4e.esurfaces[box4e.existing_esurfaces[*esurface]]
                    .compute_from_momenta(
                        &box4e.lmb,
                        &to_real_mass_vector(box4e.edge_masses.as_slice()),
                        center,
                        &box4e.external_momenta,
                    );

                assert!(esurfaec_val.0 < 0.0);
            }
        }
    }

    /// This test deforms the threshold structure into 4 pieces with no overlap
    #[test]
    fn test_disconnected_box_4e() {
        let box4e = HelperBoxStructure::new(Some([F(10.5); 4]));

        let maximal_overlap = find_maximal_overlap(
            &box4e.lmb,
            &box4e.existing_esurfaces,
            &box4e.esurfaces,
            &box4e.edge_masses,
            &box4e.external_momenta,
            0,
        );

        assert_eq!(maximal_overlap.overlap_groups.len(), 4);

        for overlap_group in maximal_overlap.overlap_groups.iter() {
            let esurfaces = &overlap_group.existing_esurfaces;
            let center = &overlap_group.center;

            assert_eq!(esurfaces.len(), 1);

            for esurface in esurfaces.iter() {
                let esurfaec_val = box4e.esurfaces[box4e.existing_esurfaces[*esurface]]
                    .compute_from_momenta(
                        &box4e.lmb,
                        &to_real_mass_vector(box4e.edge_masses.as_slice()),
                        center,
                        &box4e.external_momenta,
                    );

                assert!(esurfaec_val < F(0.0));
            }
        }
    }

    #[test]
    fn test_banana() {
        let banana = HelperBananaStructure::new();

        println!(
            "{}",
            banana.esurfaces[Into::<EsurfaceID>::into(0)].string_format_in_lmb(&banana.lmb)
        );

        let maximal_overlap = find_maximal_overlap(
            &banana.lmb,
            &banana.existing_esurfaces,
            &banana.esurfaces,
            &banana.edge_masses,
            &banana.external_momenta,
            0,
        );

        println!("center: {:#?}", maximal_overlap);

        assert_eq!(maximal_overlap.overlap_groups.len(), 1);
    }
}
