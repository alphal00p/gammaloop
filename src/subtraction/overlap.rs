use ahash::HashMap;
use ahash::HashMapExt;
use clarabel::algebra::*;
use clarabel::solver::*;
use itertools::Itertools;
use lorentz_vector::LorentzVector;
use num::Complex;
use num::Zero;
use rayon::result;

use crate::cff::esurface::ExistingEsurfaceId;
use crate::cff::esurface::ExistingEsurfaces;
use crate::cff::esurface::{EsurfaceCollection, EsurfaceId};
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
    esurfaces_to_consider: &[ExistingEsurfaceId],
    existing_esurfaces: &ExistingEsurfaces,
    esurfaces: &EsurfaceCollection,
    edge_masses: &[Option<Complex<f64>>],
    external_momenta: &[LorentzVector<f64>],
) -> DefaultSolver {
    let num_loops = lmb.basis.len();
    let num_loop_vars = 3 * num_loops;
    let num_edges = lmb.edge_signatures.len();

    // first we study the structure of the problem
    let mut propagator_constraints: Vec<PropagatorConstraint> = Vec::with_capacity(num_edges);

    let mut inequivalent_masses: Vec<Complex<f64>> = vec![];

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
    esurfaces_to_consider: &[ExistingEsurfaceId],
    existing_esurfaces: &ExistingEsurfaces,
    esurfaces: &EsurfaceCollection,
    edge_masses: &[Option<Complex<f64>>],
    external_momenta: &[LorentzVector<f64>],
) -> Option<Vec<LorentzVector<f64>>> {
    let mut solver = construct_solver(
        lmb,
        esurfaces_to_consider,
        existing_esurfaces,
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
    existing_esurfaces: &ExistingEsurfaces,
    esurfaces: &EsurfaceCollection,
    edge_masses: &[Option<Complex<f64>>],
    external_momenta: &[LorentzVector<f64>],
) -> Vec<(Vec<EsurfaceId>, Vec<LorentzVector<f64>>)> {
    let mut res = vec![];
    let num_loops = lmb.basis.len();

    let all_existing_esurfaces = existing_esurfaces.all_indices().collect_vec();

    // first try if all esurfaces have a single center, we explitely seach a center instead of trying the
    // origin. This is because the origin might not be optimal.
    let option_center = find_center(
        lmb,
        &all_existing_esurfaces,
        existing_esurfaces,
        esurfaces,
        edge_masses,
        external_momenta,
    );

    if let Some(center) = option_center {
        res.push((existing_esurfaces.clone().to_vec(), center));
        return res;
    }

    // if the center is not valid, create a table of all pairs
    let esurface_pairs = EsurfacePairs::new(
        lmb,
        existing_esurfaces,
        esurfaces,
        edge_masses,
        external_momenta,
    );

    let mut num_disconnected_surfaces = 0;

    for (existing_esurface_id, &esurface_id) in existing_esurfaces.iter_enumerate() {
        // if an esurface overlaps with no other esurface, it is part of the maximal overlap structure
        if esurface_pairs.has_no_overlap(existing_esurface_id) {
            let center = find_center(
                lmb,
                &[existing_esurface_id],
                existing_esurfaces,
                esurfaces,
                edge_masses,
                external_momenta,
            );

            res.push((
                vec![esurface_id],
                center.unwrap_or_else(|| {
                    panic!("Could not find center of esurface {:?}", esurface_id)
                }),
            ));
            num_disconnected_surfaces += 1;
        }
    }

    if num_disconnected_surfaces == existing_esurfaces.len() {
        return res;
    }

    let mut subset_size = existing_esurfaces.len() - num_disconnected_surfaces - 1;

    while subset_size > 1 {
        let possible_subsets =
            esurface_pairs.construct_possible_subsets_of_len(existing_esurfaces, subset_size);

        for subset in possible_subsets.iter() {
            if is_subset_of_result(subset, &res, existing_esurfaces) {
                continue;
            }

            let option_center = find_center(
                lmb,
                subset,
                existing_esurfaces,
                esurfaces,
                edge_masses,
                external_momenta,
            );

            if let Some(center) = option_center {
                res.push((
                    subset.iter().map(|&i| existing_esurfaces[i]).collect(),
                    center,
                ));
            }
        }

        subset_size -= 1;
    }

    res
}

fn is_subset_of_result(
    subset: &[ExistingEsurfaceId],
    result: &[(Vec<EsurfaceId>, Vec<LorentzVector<f64>>)],
    existing_esurfaces: &ExistingEsurfaces,
) -> bool {
    result.iter().any(|(result_subset, _)| {
        subset
            .iter()
            .all(|&x| result_subset.contains(&existing_esurfaces[x]))
    })
}

#[derive(Debug)]
struct EsurfacePairs {
    data: HashMap<(ExistingEsurfaceId, ExistingEsurfaceId), Vec<LorentzVector<f64>>>,
    has_pair_with: Vec<Vec<ExistingEsurfaceId>>,
}

impl EsurfacePairs {
    fn insert(
        &mut self,
        pair: (ExistingEsurfaceId, ExistingEsurfaceId),
        center: Vec<LorentzVector<f64>>,
    ) {
        if pair.0 > pair.1 {
            self.data.insert((pair.1, pair.0), center);
        } else {
            self.data.insert(pair, center);
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
        edge_masses: &[Option<Complex<f64>>],
        external_momenta: &[LorentzVector<f64>],
    ) -> Self {
        let mut res = Self::new_empty(existing_esurfaces.len());

        let all_existing_esurfaces = existing_esurfaces.all_indices().collect_vec();

        for (i, &esurface_id_1) in all_existing_esurfaces.iter().enumerate() {
            for &esurface_id_2 in all_existing_esurfaces.iter().skip(i + 1) {
                let center = find_center(
                    lmb,
                    &[esurface_id_1, esurface_id_2],
                    existing_esurfaces,
                    esurfaces,
                    edge_masses,
                    external_momenta,
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
    ) -> Vec<Vec<ExistingEsurfaceId>> {
        let mut res = vec![];
        let num_existing_esurfaces = self.has_pair_with.len();
        assert!(subset_len < num_existing_esurfaces);

        let subsets = existing_esurfaces.all_indices().combinations(subset_len);

        // horrendously inefficient, many subsets can be skipped once one is found to not be valid.
        for subset in subsets {
            let is_valid = subset.iter().combinations(2).all(|pair| {
                assert!(pair[0] < pair[1]);
                let pair = (
                    Into::<ExistingEsurfaceId>::into(*pair[0]),
                    Into::<ExistingEsurfaceId>::into(*pair[1]),
                );
                self.data.contains_key(&pair)
            });

            if is_valid {
                res.push(subset);
            }
        }

        res
    }
}

#[cfg(test)]
#[allow(dead_code, unused_variables)]
mod tests {
    use super::*;
    use itertools::Itertools;
    use lorentz_vector::LorentzVector;
    use num::Complex;

    fn to_real_mass_vector(edge_masses: &[Option<Complex<f64>>]) -> Vec<f64> {
        edge_masses
            .iter()
            .map(|mass| match mass {
                Some(m) => m.re,
                None => 0.0,
            })
            .collect_vec()
    }

    use crate::{
        cff::esurface::{self, Esurface},
        graph::LoopMomentumBasis,
    };

    struct HelperBoxStructure {
        external_momenta: [LorentzVector<f64>; 3],
        lmb: LoopMomentumBasis,
        esurfaces: EsurfaceCollection,
        existing_esurfaces: ExistingEsurfaces,
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

            let esurfaces_array = [
                Esurface {
                    energies: vec![5, 6],
                    sub_orientation: vec![],
                    shift: vec![1],
                    shift_signature: vec![true],
                },
                Esurface {
                    energies: vec![5, 7],
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

            let esurfaces = EsurfaceCollection::from_vec(esurfaces_array.to_vec());

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

            let existing_esurfaces =
                ExistingEsurfaces::from_vec((0..4).map(Into::<EsurfaceId>::into).collect());

            Self {
                external_momenta,
                lmb: box_lmb,
                existing_esurfaces,
                esurfaces,
                edge_masses,
            }
        }
    }

    #[test]
    fn test_is_subset_of_result() {
        let fake_res = vec![
            (
                vec![
                    Into::<EsurfaceId>::into(1),
                    Into::<EsurfaceId>::into(2),
                    Into::<EsurfaceId>::into(3),
                ],
                vec![],
            ),
            (
                vec![
                    Into::<EsurfaceId>::into(1),
                    Into::<EsurfaceId>::into(2),
                    Into::<EsurfaceId>::into(4),
                ],
                vec![],
            ),
            (
                vec![
                    Into::<EsurfaceId>::into(2),
                    Into::<EsurfaceId>::into(3),
                    Into::<EsurfaceId>::into(4),
                ],
                vec![],
            ),
        ];

        let existing_esurfaces =
            ExistingEsurfaces::from_vec((0..5).map(Into::<EsurfaceId>::into).collect());

        let fake_subset = vec![
            Into::<ExistingEsurfaceId>::into(1),
            Into::<ExistingEsurfaceId>::into(2),
        ];

        assert!(is_subset_of_result(
            &fake_subset,
            &fake_res,
            &existing_esurfaces
        ));

        let fake_subset_2 = vec![
            Into::<ExistingEsurfaceId>::into(0),
            Into::<ExistingEsurfaceId>::into(4),
        ];

        assert!(!is_subset_of_result(
            &fake_subset_2,
            &fake_res,
            &existing_esurfaces
        ));
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

        let box4e_massive = HelperBoxStructure::new(Some([10.5; 4]));
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

        let subsets_3 =
            esurface_pairs.construct_possible_subsets_of_len(&box4e.existing_esurfaces, 3);

        assert_eq!(subsets_3.len(), 0);

        let subsets_2 =
            esurface_pairs.construct_possible_subsets_of_len(&box4e.existing_esurfaces, 2);
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
        );

        assert_eq!(maximal_overlap.len(), 4);

        for (esurfaces, center) in maximal_overlap.iter() {
            assert_eq!(esurfaces.len(), 2);

            for esurface in esurfaces.iter() {
                let esurfaec_val = box4e.esurfaces[*esurface].compute_from_momenta(
                    &box4e.lmb,
                    &to_real_mass_vector(box4e.edge_masses.as_slice()),
                    center,
                    &box4e.external_momenta,
                );

                assert!(esurfaec_val < 0.0);
            }
        }
    }

    /// This test deforms the threshold structure into 4 pieces with no overlap
    #[test]
    fn test_disconnected_box_4e() {
        let box4e = HelperBoxStructure::new(Some([10.5; 4]));

        let maximal_overlap = find_maximal_overlap(
            &box4e.lmb,
            &box4e.existing_esurfaces,
            &box4e.esurfaces,
            &box4e.edge_masses,
            &box4e.external_momenta,
        );

        assert_eq!(maximal_overlap.len(), 4);

        for (esurfaces, center) in maximal_overlap.iter() {
            assert_eq!(esurfaces.len(), 1);

            for esurface in esurfaces.iter() {
                let esurfaec_val = box4e.esurfaces[*esurface].compute_from_momenta(
                    &box4e.lmb,
                    &to_real_mass_vector(box4e.edge_masses.as_slice()),
                    center,
                    &box4e.external_momenta,
                );

                assert!(esurfaec_val < 0.0);
            }
        }
    }
}
