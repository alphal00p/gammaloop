use crate::GammaLoopContext;
use crate::cff::esurface::EsurfaceCollection;
use crate::cff::esurface::ExistingEsurfaceId;
use crate::cff::esurface::ExistingThresholds;
use crate::cff::esurface::esurface_value_is_strictly_inside;
use crate::graph::FeynmanGraph;
use crate::graph::Graph;
use crate::graph::LmbIndex;
use crate::graph::LoopMomentumBasis;
use crate::momentum::sample::ExternalFourMomenta;
use crate::momentum::sample::LoopIndex;
use crate::momentum::sample::LoopMomenta;
use crate::momentum::sample::SubspaceData;
use crate::momentum::signature::LoopExtSignature;
use crate::momentum::{Rotation, ThreeMomentum};
use crate::settings::RuntimeSettings;
use crate::utils::F;
use crate::utils::compute_shift_part_subspace;
use ahash::HashMap;
use ahash::HashMapExt;
use ahash::HashSet;
use bincode_trait_derive::Decode;
use bincode_trait_derive::Encode;
use clarabel::algebra::*;
use clarabel::solver::*;
use eyre::{Result, eyre};
use itertools::Itertools;
use linnet::half_edge::involution::EdgeVec;
use spenso::algebra::algebraic_traits::IsZero;
use std::fmt::Display;
use typed_index_collections::TiVec;

#[derive(Debug, Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct OverlapGroup {
    pub existing_esurfaces: Vec<ExistingEsurfaceId>,
    pub complement: Vec<ExistingEsurfaceId>,
    /// LU overlap centers are stored in the current probe and cut-side LMB frame.
    /// Solver-derived centers therefore require no further rotation at consumption.
    pub center: LoopMomenta<F<f64>>,
}

#[derive(Debug, Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct OverlapStructure {
    pub overlap_groups: Vec<OverlapGroup>,
    pub existing_esurfaces: ExistingThresholds,
}

impl Display for OverlapStructure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let existing_esurfaces: Vec<_> = self.existing_esurfaces.iter().map(|id| id.0).collect();
        writeln!(f, "existing esurfaces: {:?}", existing_esurfaces)?;

        for (i, group) in self.overlap_groups.iter().enumerate() {
            writeln!(f, "Group {}:", i)?;
            let existing_esurfaces_in_group: Vec<_> = group
                .existing_esurfaces
                .iter()
                .map(|id| self.existing_esurfaces[*id].0)
                .collect();

            writeln!(f, "center:\n{}", group.center)?;
            writeln!(
                f,
                "existing esurfaces in group: {:?}",
                existing_esurfaces_in_group
            )?;
        }

        Ok(())
    }
}

impl OverlapStructure {
    pub fn fill_in_complements(&mut self) {
        for group in self.overlap_groups.iter_mut() {
            group.complement = self
                .existing_esurfaces
                .iter_enumerated()
                .map(|(existing_esurface_id, _)| existing_esurface_id)
                .filter(|&existing_esurface_id| {
                    !group.existing_esurfaces.contains(&existing_esurface_id)
                })
                .collect();
        }
    }

    pub fn new_empty() -> Self {
        Self {
            overlap_groups: vec![],
            existing_esurfaces: ExistingThresholds::new(),
        }
    }
}
/// Helper struct to construct the socp problem
struct PropagatorConstraint<'a> {
    mass_pointer: Option<usize>, // pointer to value of unique mass
    signature: &'a LoopExtSignature,
}

impl PropagatorConstraint<'_> {
    fn get_dimension(&self) -> usize {
        let mass_value = if self.mass_pointer.is_some() { 1 } else { 0 };

        3 + mass_value + 1
    }
}

fn extract_center(
    global_loop_num: usize,
    subspace: &SubspaceData,
    solution: &[f64],
) -> LoopMomenta<F<f64>> {
    let len = solution.len();
    let num_loop_vars = 3 * subspace.loopcount();

    let mut loop_chunks = solution[len - num_loop_vars..].chunks(3);

    (0..global_loop_num)
        .map(LoopIndex::from)
        .map(|loop_index| {
            if subspace.contains_loop_index(loop_index) {
                let window = loop_chunks.next().expect("not enough loop momenta");
                ThreeMomentum::new(
                    F::from_f64(window[0]),
                    F::from_f64(window[1]),
                    F::from_f64(window[2]),
                )
            } else {
                ThreeMomentum::new(F(0.0), F(0.0), F(0.0))
            }
        })
        .collect()
}

fn construct_solver(
    overlap_input: &OverlapInput,
    esurfaces_to_consider: &[ExistingEsurfaceId],
    existing_esurfaces: &ExistingThresholds,
    loop_moms: &LoopMomenta<F<f64>>,
    external_momenta: &ExternalFourMomenta<F<f64>>,
    verbose: bool,
) -> DefaultSolver {
    let num_loops = overlap_input.subspace.loopcount();

    let num_loop_vars = 3 * num_loops;

    // first we study the structure of the problem
    let mut propagator_constraints: Vec<PropagatorConstraint> = Vec::with_capacity(20);

    let mut inequivalent_masses: Vec<F<f64>> = vec![];

    let mut esurface_constraints: Vec<Vec<usize>> = Vec::with_capacity(esurfaces_to_consider.len());

    for existing_esurface_id in esurfaces_to_consider.iter() {
        let surface_id = existing_esurfaces[*existing_esurface_id];

        let esurface = &overlap_input.thresholds[surface_id];
        let lmb = overlap_input.subspace.get_lmb(overlap_input.lmbs);
        let edge_masses = &overlap_input.edge_masses;

        let mut esurface_constraint_indices: Vec<usize> = Vec::with_capacity(6);

        for edge_id in overlap_input
            .subspace
            .contains(&esurface.energies, overlap_input.graph)
        {
            if let Some(edge_position) = propagator_constraints
                .iter()
                .position(|constraint| *constraint.signature == lmb.edge_signatures[edge_id])
            {
                esurface_constraint_indices.push(edge_position);
            } else {
                let mass_pointer = if edge_masses[edge_id].is_zero() {
                    None
                } else {
                    Some(
                        if let Some(mass_position) = inequivalent_masses
                            .iter()
                            .position(|&x| x == edge_masses[edge_id])
                        {
                            mass_position
                        } else {
                            inequivalent_masses.push(edge_masses[edge_id]);
                            inequivalent_masses.len() - 1
                        },
                    )
                };

                let signature = &lmb.edge_signatures[edge_id];

                let propagator_constraint = PropagatorConstraint {
                    mass_pointer,
                    signature,
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
    for (constraint_index, (existing_esurface_id, esurface_constraint)) in esurfaces_to_consider
        .iter()
        .zip(esurface_constraints.iter())
        .enumerate()
    {
        for prop_index in esurface_constraint {
            a_matrix[constraint_index + 1][*prop_index + propagator_index_offset] = 1.0;
        }

        let esurface_id = existing_esurfaces[*existing_esurface_id];
        let esurface = &overlap_input.thresholds[esurface_id];

        let shift_part = esurface.compute_shift_part_from_momenta_in_subspace(
            loop_moms,
            external_momenta,
            overlap_input.subspace,
            overlap_input.lmbs,
            overlap_input.graph,
            &overlap_input.edge_masses,
        );
        b_vector[constraint_index + 1] = -shift_part.0;
        a_matrix[constraint_index + 1][0] = -1.0;
    }

    let spatial_part_of_externals = external_momenta.iter().map(|p| p.spatial).collect();

    // propagator constraints
    let mut vertical_offset = esurface_constraints.len() + 1;
    for (cone_index, propagator_constraint) in propagator_constraints.iter().enumerate() {
        a_matrix[vertical_offset][propagator_index_offset + cone_index] = -1.0;
        vertical_offset += 1;

        let spatial_shift = compute_shift_part_subspace(
            &propagator_constraint.signature.internal,
            &propagator_constraint.signature.external,
            loop_moms,
            &spatial_part_of_externals,
            overlap_input.subspace,
        );

        b_vector[vertical_offset] = spatial_shift.px.0;
        b_vector[vertical_offset + 1] = spatial_shift.py.0;
        b_vector[vertical_offset + 2] = spatial_shift.pz.0;

        for (subspace_loop_index, individual_loop_signature) in overlap_input
            .subspace
            .project_loop_signature_filtered(&propagator_constraint.signature.internal)
            .enumerate()
        {
            if individual_loop_signature.is_sign() {
                a_matrix[vertical_offset][loop_momentum_offset + 3 * subspace_loop_index] =
                    -(individual_loop_signature as i8) as f64;
                a_matrix[vertical_offset + 1][loop_momentum_offset + 3 * subspace_loop_index + 1] =
                    -(individual_loop_signature as i8) as f64;
                a_matrix[vertical_offset + 2][loop_momentum_offset + 3 * subspace_loop_index + 2] =
                    -(individual_loop_signature as i8) as f64;
            }
        }

        vertical_offset += 3;

        if let Some(mass_index) = propagator_constraint.mass_pointer {
            b_vector[vertical_offset] = inequivalent_masses[mass_index].0;
            vertical_offset += 1;
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
    .unwrap()
}

pub(crate) fn find_center(
    overlap_input: &OverlapInput,
    esurfaces_to_consider: &[ExistingEsurfaceId],
    existing_esurfaces: &ExistingThresholds,
    loop_moms: &LoopMomenta<F<f64>>,
    external_momenta: &ExternalFourMomenta<F<f64>>,
    verbose: bool,
) -> Option<LoopMomenta<F<f64>>> {
    let mut solver = construct_solver(
        overlap_input,
        esurfaces_to_consider,
        existing_esurfaces,
        loop_moms,
        external_momenta,
        verbose,
    );

    solver.solve();

    let global_loop_number = overlap_input.graph.get_loop_number();
    let esurfaces_to_check = esurfaces_to_consider
        .iter()
        .map(|existing_esurface_id| existing_esurfaces[*existing_esurface_id])
        .collect();

    if solver.solution.status == SolverStatus::Solved {
        let center = extract_center(
            global_loop_number,
            overlap_input.subspace,
            &solver.solution.x,
        );
        check_global_center(
            overlap_input,
            &esurfaces_to_check,
            &center,
            loop_moms,
            external_momenta,
        )
        .then_some(center)
    } else if solver.solution.status == SolverStatus::AlmostSolved
        || solver.solution.status == SolverStatus::InsufficientProgress
    {
        // if the solver did not converge, we check if the solution is still valid
        let center = extract_center(
            global_loop_number,
            overlap_input.subspace,
            &solver.solution.x,
        );

        let is_valid = check_global_center(
            overlap_input,
            &esurfaces_to_check,
            &center,
            loop_moms,
            external_momenta,
        );

        if is_valid { Some(center) } else { None }
    } else {
        if verbose {
            println!("{:?}", solver.solution.x);
        }

        None
    }
}

pub(crate) struct OverlapInput<'a> {
    pub graph: &'a Graph,
    pub settings: &'a RuntimeSettings,
    pub subspace: &'a SubspaceData,
    pub lmbs: &'a TiVec<LmbIndex, LoopMomentumBasis>,
    pub thresholds: &'a EsurfaceCollection,
    pub edge_masses: EdgeVec<F<f64>>,
}

pub(crate) fn check_global_center(
    overlap_input: &OverlapInput,
    existing_esurfaces: &ExistingThresholds,
    center: &LoopMomenta<F<f64>>,
    loop_moms: &LoopMomenta<F<f64>>,
    external_momenta: &ExternalFourMomenta<F<f64>>,
) -> bool {
    let mut center_with_fixed_complement = loop_moms.clone();
    for loop_index in overlap_input.subspace.iter_lmb_indices() {
        center_with_fixed_complement[loop_index] = center[loop_index];
    }

    existing_esurfaces.iter().all(|esurface_id| {
        let esurface = &overlap_input.thresholds[*esurface_id];

        let lmb = overlap_input.subspace.get_lmb(overlap_input.lmbs);
        let edge_masses = &overlap_input.edge_masses;

        let esurface_val = esurface.compute_from_momenta(
            lmb,
            edge_masses,
            &center_with_fixed_complement,
            external_momenta,
        );

        esurface_value_is_strictly_inside(&esurface_val, &F(overlap_input.settings.kinematics.e_cm))
    })
}

/// Runtime overlap failures are returned so the stability machinery can retry at higher precision.
/// Structural generation invariants are still asserted where malformed generated data is unrecoverable.
/// Solver-derived centers are already found in the current probe frame. `probe_rotation` is
/// needed only for a configured forced center, whose coordinates are defined in the identity
/// frame and must therefore be rotated exactly once before the cut-side LMB transform.
pub(crate) fn find_maximal_overlap(
    overlap_input: &OverlapInput,
    existing_esurfaces: &ExistingThresholds,
    loop_moms: &LoopMomenta<F<f64>>,
    external_momenta: &ExternalFourMomenta<F<f64>>,
    probe_rotation: &Rotation,
) -> Result<OverlapStructure> {
    let mut res = OverlapStructure {
        overlap_groups: vec![],
        existing_esurfaces: existing_esurfaces.clone(),
    };

    let settings = overlap_input.settings;

    let all_existing_esurfaces = existing_esurfaces
        .iter_enumerated()
        .map(|a| a.0)
        .collect_vec();

    if let Some(global_center) = &settings.subtraction.overlap_settings.force_global_center {
        let global_center_identity: LoopMomenta<F<f64>> = global_center
            .iter()
            .map(|coordinates| ThreeMomentum {
                px: F(coordinates[0]),
                py: F(coordinates[1]),
                pz: F(coordinates[2]),
            })
            .collect();
        let rotated_external_spatial = external_momenta
            .iter()
            .map(|momentum| momentum.spatial)
            .collect();
        // Forced centers are configured in the graph's identity-frame LMB. Apply the probe
        // rotation once, then express that rotated point in the cut-side LMB used by the solver.
        let global_center_probe = global_center_identity.rotate(probe_rotation).lmb_transform(
            &overlap_input.graph.loop_momentum_basis,
            overlap_input.subspace.get_lmb(overlap_input.lmbs),
            &rotated_external_spatial,
        );

        tracing::debug!(
            graph = %overlap_input.graph.name,
            rotation_id = %probe_rotation.method,
            center_provenance = "forced_identity_frame_rotated_and_lmb_transformed_once",
            identity_center = %global_center_identity,
            probe_cut_lmb_center = %global_center_probe,
            "prepared forced LU overlap center"
        );

        if !settings.subtraction.overlap_settings.check_global_center {
            tracing::warn!(
                graph = %overlap_input.graph.name,
                "overlap_settings.check_global_center=false is deprecated; forced centers are always validated"
            );
        }

        let is_valid = check_global_center(
            overlap_input,
            existing_esurfaces,
            &global_center_probe,
            loop_moms,
            external_momenta,
        );

        if !is_valid {
            return Err(eyre!(
                "Forced identity-frame center is not finite and strictly inside all existing esurfaces after applying probe rotation {} and the cut-side LMB transform",
                probe_rotation.method,
            ));
        }

        let single_group = OverlapGroup {
            existing_esurfaces: all_existing_esurfaces,
            center: global_center_probe,
            complement: vec![],
        };
        res.overlap_groups.push(single_group);

        res.fill_in_complements();
        return Ok(res);
    }

    if settings.subtraction.overlap_settings.try_origin {
        let global_loop_count = overlap_input.graph.get_loop_number();
        let origin = LoopMomenta::from_iter(
            (0..global_loop_count).map(|_| ThreeMomentum::new(F(0.0), F(0.0), F(0.0))),
        );

        let is_valid = check_global_center(
            overlap_input,
            existing_esurfaces,
            &origin,
            loop_moms,
            external_momenta,
        );

        if is_valid {
            let single_group = OverlapGroup {
                existing_esurfaces: all_existing_esurfaces,
                center: origin,
                complement: vec![],
            };
            res.overlap_groups.push(single_group);
            res.fill_in_complements();
            return Ok(res);
        }
    }

    if settings.subtraction.overlap_settings.try_origin_all_lmbs {
        todo!("Not all heuristics implemented")
    }

    // first try if all esurfaces have a single center, we explitely seach a center instead of trying the
    // origin. This is because the origin might not be optimal.
    let option_center = find_center(
        overlap_input,
        &all_existing_esurfaces,
        existing_esurfaces,
        loop_moms,
        external_momenta,
        false,
    );

    if let Some(center) = option_center {
        let single_group = OverlapGroup {
            existing_esurfaces: all_existing_esurfaces,
            center,
            complement: vec![],
        };
        res.overlap_groups.push(single_group);
        res.fill_in_complements();
        return Ok(res);
    }

    // if the center is not valid, create a table of all pairs
    let esurface_pairs = EsurfacePairs::new(
        overlap_input,
        existing_esurfaces,
        loop_moms,
        external_momenta,
    );

    // if settings.general.debug > 3 {
    //     DEBUG_LOGGER.write("overlap_pairs", &esurface_pairs);
    // }

    let mut num_disconnected_surfaces = 0;

    for (existing_esurface_id, &esurface_id) in existing_esurfaces.iter_enumerated() {
        // if an esurface overlaps with no other esurface, it is part of the maximal overlap structure
        if esurface_pairs.has_no_overlap(existing_esurface_id) {
            let center = find_center(
                overlap_input,
                &[existing_esurface_id],
                existing_esurfaces,
                loop_moms,
                external_momenta,
                false,
            )
            .ok_or_else(|| {
                let esurface = &overlap_input.thresholds[esurface_id];

                let mut error_message = String::new();

                error_message.push_str(&format!(
                    "Could not find center of esuface {:?}\n",
                    esurface_id
                ));
                error_message.push_str(&format!("edges: {:?}\n", esurface.energies));

                error_message.push_str(&format!("external shift: {:?}\n", esurface.external_shift));

                error_message.push_str(&format!("External momenta: {:#?}\n", external_momenta));

                eyre!("{}", error_message)
            })?;

            res.overlap_groups.push(OverlapGroup {
                existing_esurfaces: vec![existing_esurface_id],
                center,
                complement: vec![],
            });
            num_disconnected_surfaces += 1;
        }
    }

    // if settings.general.debug > 3 {
    //     DEBUG_LOGGER.write("num_disconnected_surfaces", &num_disconnected_surfaces);
    // }

    if num_disconnected_surfaces == existing_esurfaces.len() {
        res.fill_in_complements();
        return Ok(res);
    }

    let mut subset_size =
        if let Some(size) = esurface_pairs.has_pair_with.iter().map(Vec::len).max() {
            size + 1
        } else {
            1
        };

    while subset_size > 1 {
        let possible_subsets =
            esurface_pairs.construct_possible_subsets_of_len(existing_esurfaces, subset_size, &res);

        for subset in possible_subsets.iter() {
            let option_center = find_center(
                overlap_input,
                subset,
                existing_esurfaces,
                loop_moms,
                external_momenta,
                false,
            );

            if let Some(center) = option_center {
                res.overlap_groups.push(OverlapGroup {
                    existing_esurfaces: subset.clone(),
                    center,
                    complement: vec![],
                });
            }
        }

        // if settings.general.debug > 3 {
        //     DEBUG_LOGGER.write(
        //         "subset_size_and_num_possible_subsets_and_res",
        //         &(subset_size, possible_subsets.len(), &res),
        //     );
        // }

        subset_size -= 1;
    }

    res.fill_in_complements();
    Ok(res)
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
    data: HashMap<(ExistingEsurfaceId, ExistingEsurfaceId), LoopMomenta<F<f64>>>,
    has_pair_with: Vec<Vec<ExistingEsurfaceId>>,
}

impl EsurfacePairs {
    fn insert(
        &mut self,
        pair: (ExistingEsurfaceId, ExistingEsurfaceId),
        center: LoopMomenta<F<f64>>,
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
        let capacity = match num_existing_esurfaces {
            0 => 0,
            1 => 0,
            _ => num_existing_esurfaces * (num_existing_esurfaces - 1) / 2,
        };

        Self {
            data: HashMap::with_capacity(capacity),
            has_pair_with: vec![Vec::with_capacity(num_existing_esurfaces); num_existing_esurfaces],
        }
    }

    fn new(
        overlap_input: &OverlapInput,
        existing_esurfaces: &ExistingThresholds,
        loop_moms: &LoopMomenta<F<f64>>,
        external_momenta: &ExternalFourMomenta<F<f64>>,
    ) -> Self {
        let mut res = Self::new_empty(existing_esurfaces.len());

        let all_existing_esurfaces = existing_esurfaces
            .iter_enumerated()
            .map(|a| a.0)
            .collect_vec();

        for (i, &esurface_id_1) in all_existing_esurfaces.iter().enumerate() {
            for &esurface_id_2 in all_existing_esurfaces.iter().skip(i + 1) {
                let center = find_center(
                    overlap_input,
                    &[esurface_id_1, esurface_id_2],
                    existing_esurfaces,
                    loop_moms,
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
        existing_esurfaces: &ExistingThresholds,
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        cff::{
            cff_graph::VertexSet,
            esurface::{Esurface, EsurfaceExistence, EsurfaceID},
        },
        dot,
        graph::{LMBext, parse::from_dot::IntoGraph},
        initialisation::test_initialise,
        momentum::{FourMomentum, Rotatable, Rotation, RotationMethod},
    };
    use linnet::half_edge::involution::EdgeIndex;
    use linnet::half_edge::subgraph::{SuBitGraph, SubSetOps};
    use typed_index_collections::ti_vec;

    #[test]
    fn global_center_check_preserves_fixed_complement_for_multidimensional_subspace() {
        test_initialise().unwrap();
        let graph: Graph = dot!(digraph subspace_center {
            ext [style=invis]
            edge [num=1 mass=0]
            node [num=1]
            ext->v1:0 [id=0]
            v1->v2 [id=1]
            v2->v1 [id=2]
            v1->v2 [id=3]
            v2->v1 [id=4]
            ext->v2:1 [id=5]
        })
        .unwrap();

        assert_eq!(graph.loop_momentum_basis.loop_edges.len(), 3);
        let all_lmbs = ti_vec![graph.loop_momentum_basis.clone()];

        let mut parallel_edges_only: SuBitGraph = graph.empty_subgraph();
        for active_loop_index in [LoopIndex(0), LoopIndex(1)] {
            let active_loop_edge = graph.loop_momentum_basis.loop_edges[active_loop_index];
            parallel_edges_only.union_with(&graph.get_edge_subgraph(active_loop_edge));
        }
        let parallel_edge_subspace = SubspaceData::new_with_user_selected_lmb(
            parallel_edges_only,
            LmbIndex::from(0),
            &graph,
            &all_lmbs,
        )
        .unwrap();
        assert_eq!(
            parallel_edge_subspace.loopcount(),
            1,
            "two graph-parallel defining edges without their spanning support contain only one independent loop"
        );

        let raised_graph: Graph = dot!(digraph raised_signature_subspace {
            ext [style=invis]
            edge [num=1 mass=0]
            node [num=1]
            ext->a [id=0]
            a->b [id=1]
            b->c [id=2]
            c->a [id=3]
            ext->c [id=4]
        })
        .unwrap();
        let raised_group = raised_graph
            .get_raised_edge_groups()
            .into_iter()
            .find(|group| group.len() >= 2)
            .expect("test graph must contain a raised equal-signature edge group");
        assert!(raised_group.iter().tuple_windows().all(|(left, right)| {
            raised_graph.loop_momentum_basis.edge_signatures[*left]
                .equality_up_to_sign(&raised_graph.loop_momentum_basis.edge_signatures[*right])
        }));
        let mut raised_subgraph: SuBitGraph = raised_graph.empty_subgraph();
        for &edge in &raised_group {
            raised_subgraph.union_with(&raised_graph.get_edge_subgraph(edge));
        }
        let raised_lmbs = ti_vec![raised_graph.loop_momentum_basis.clone()];
        let raised_edges_only = SubspaceData::new_with_user_selected_lmb(
            raised_subgraph.clone(),
            LmbIndex::from(0),
            &raised_graph,
            &raised_lmbs,
        )
        .unwrap();
        assert_eq!(
            raised_edges_only.loopcount(),
            0,
            "a chain of raised equal-signature edges is not multiple independent loops"
        );

        let support_edge = raised_graph
            .iter_loop_edges()
            .map(|(_, edge, _)| edge)
            .find(|edge| !raised_group.contains(edge))
            .expect("test graph must contain the support edge closing the loop");
        raised_subgraph.union_with(&raised_graph.get_edge_subgraph(support_edge));
        let raised_cycle = SubspaceData::new_with_user_selected_lmb(
            raised_subgraph,
            LmbIndex::from(0),
            &raised_graph,
            &raised_lmbs,
        )
        .unwrap();
        assert_eq!(
            raised_cycle.loopcount(),
            1,
            "raised equal-signature edges must be counted by topology, not once per edge"
        );

        // Include the parent LMB spanning-tree support. Two graph-parallel
        // defining edges alone contain only one independent cycle.
        let mut subgraph = graph.loop_momentum_basis.tree.clone();
        for active_loop_index in [LoopIndex(0), LoopIndex(1)] {
            let active_loop_edge = graph.loop_momentum_basis.loop_edges[active_loop_index];
            subgraph.union_with(&graph.get_edge_subgraph(active_loop_edge));
        }
        let subspace = SubspaceData::new_with_user_selected_lmb(
            subgraph,
            LmbIndex::from(0),
            &graph,
            &all_lmbs,
        )
        .unwrap();
        assert_eq!(subspace.loopcount(), 2);

        let radial_surface = Esurface {
            energies: vec![graph.loop_momentum_basis.loop_edges[LoopIndex(0)]],
            external_shift: vec![(EdgeIndex::from(0), -1)],
            vertex_set: VertexSet::dummy(),
        };
        assert!(radial_surface.has_radial_dependence_in_subspace(&subspace, &all_lmbs, &graph,));
        let zero_loop_momenta = LoopMomenta::from_iter([
            ThreeMomentum::new(F(0.0), F(0.0), F(0.0)),
            ThreeMomentum::new(F(0.0), F(0.0), F(0.0)),
            ThreeMomentum::new(F(0.0), F(0.0), F(0.0)),
        ]);
        let classify_radial_surface = |energy| {
            let external_momenta = ExternalFourMomenta::from_iter([
                FourMomentum::from_args(F(energy), F(10.0), F(0.0), F(0.0)),
                FourMomentum::from_args(F(-energy), F(-10.0), F(0.0), F(0.0)),
            ]);
            radial_surface.classify_existence_subspace(
                &zero_loop_momenta,
                &external_momenta,
                &subspace,
                &all_lmbs,
                &graph,
                &graph.underlying.new_edgevec(|_, _, _| F(0.0)),
                &[],
                &F(10.0),
                &F(crate::utils::DEFAULT_ESURFACE_EXISTENCE_THRESHOLD),
            )
        };
        assert!(matches!(
            classify_radial_surface(11.0),
            EsurfaceExistence::Existing { .. }
        ));
        assert!(matches!(
            classify_radial_surface(10.0),
            EsurfaceExistence::Pinched { .. }
        ));
        assert!(matches!(
            classify_radial_surface(9.0),
            EsurfaceExistence::NonExisting { .. }
        ));

        let complement_edge = graph.loop_momentum_basis.loop_edges[LoopIndex(2)];
        let thresholds: crate::cff::esurface::EsurfaceCollection = vec![Esurface {
            energies: vec![complement_edge],
            external_shift: vec![(EdgeIndex::from(0), -1)],
            vertex_set: VertexSet::dummy(),
        }]
        .into();
        assert!(
            !thresholds[EsurfaceID::from(0)]
                .has_radial_dependence_in_subspace(&subspace, &all_lmbs, &graph,)
        );
        let masses = graph.underlying.new_edgevec(|_, _, _| F(0.0));
        let settings = RuntimeSettings::default();
        let overlap_input = OverlapInput {
            graph: &graph,
            settings: &settings,
            subspace: &subspace,
            lmbs: &all_lmbs,
            thresholds: &thresholds,
            edge_masses: masses,
        };

        let external_momenta = ExternalFourMomenta::from_iter([
            FourMomentum::from_args(F(1.0), F(0.0), F(0.0), F(0.0)),
            FourMomentum::from_args(F(-1.0), F(0.0), F(0.0), F(0.0)),
        ]);
        let center = LoopMomenta::from_iter([
            ThreeMomentum::new(F(0.0), F(0.0), F(0.0)),
            ThreeMomentum::new(F(0.0), F(0.0), F(0.0)),
            ThreeMomentum::new(F(0.0), F(0.0), F(0.0)),
        ]);
        let sampled_momenta = LoopMomenta::from_iter([
            ThreeMomentum::new(F(0.0), F(0.0), F(0.0)),
            ThreeMomentum::new(F(0.0), F(0.0), F(0.0)),
            ThreeMomentum::new(F(2.0), F(0.0), F(0.0)),
        ]);

        let origin_value = overlap_input.thresholds[EsurfaceID::from(0)].compute_from_momenta(
            subspace.get_lmb(&all_lmbs),
            &overlap_input.edge_masses,
            &center,
            &external_momenta,
        );
        assert!(origin_value < F(0.0));
        assert!(!check_global_center(
            &overlap_input,
            &ti_vec![EsurfaceID::from(0)],
            &center,
            &sampled_momenta,
            &external_momenta,
        ));
        assert!(check_global_center(
            &overlap_input,
            &ti_vec![EsurfaceID::from(0)],
            &center,
            &center,
            &external_momenta,
        ));

        let forced_center_coordinates = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let identity_frame_center =
            LoopMomenta::from_iter(forced_center_coordinates.iter().map(|coordinates| {
                ThreeMomentum::new(F(coordinates[0]), F(coordinates[1]), F(coordinates[2]))
            }));
        let probe_rotation = Rotation::new(RotationMethod::Pi2Z);
        let expected_probe_center = identity_frame_center.rotate(&probe_rotation);
        let mut forced_settings = RuntimeSettings::default();
        forced_settings
            .subtraction
            .overlap_settings
            .force_global_center = Some(forced_center_coordinates);
        let forced_overlap_input = OverlapInput {
            graph: &graph,
            settings: &forced_settings,
            subspace: &subspace,
            lmbs: &all_lmbs,
            thresholds: &thresholds,
            edge_masses: overlap_input.edge_masses.clone(),
        };

        let forced_overlap = find_maximal_overlap(
            &forced_overlap_input,
            &ti_vec![EsurfaceID::from(0)],
            &center,
            &external_momenta,
            &probe_rotation,
        )
        .unwrap();
        assert_eq!(forced_overlap.overlap_groups.len(), 1);
        assert_eq!(
            forced_overlap.overlap_groups[0].center, expected_probe_center,
            "an identity-frame forced center must be rotated exactly once into the probe frame"
        );

        let alternate_lmb = graph
            .generate_loop_momentum_bases()
            .into_iter()
            .find(|lmb| lmb.loop_edges != graph.loop_momentum_basis.loop_edges)
            .expect("test graph must admit a non-default parent LMB");
        let alternate_lmbs = ti_vec![graph.loop_momentum_basis.clone(), alternate_lmb];
        let full_subspace = SubspaceData::new_with_user_selected_lmb(
            graph.full_filter(),
            LmbIndex::from(1),
            &graph,
            &alternate_lmbs,
        )
        .unwrap();
        let empty_thresholds: EsurfaceCollection = Vec::new().into();
        let alternate_overlap_input = OverlapInput {
            graph: &graph,
            settings: &forced_settings,
            subspace: &full_subspace,
            lmbs: &alternate_lmbs,
            thresholds: &empty_thresholds,
            edge_masses: overlap_input.edge_masses.clone(),
        };
        let rotated_external_spatial = external_momenta
            .iter()
            .map(|momentum| momentum.spatial)
            .collect();
        let expected_alternate_lmb_center = expected_probe_center.lmb_transform(
            &graph.loop_momentum_basis,
            full_subspace.get_lmb(&alternate_lmbs),
            &rotated_external_spatial,
        );
        let alternate_overlap = find_maximal_overlap(
            &alternate_overlap_input,
            &ti_vec![],
            &center,
            &external_momenta,
            &probe_rotation,
        )
        .unwrap();
        assert_eq!(
            alternate_overlap.overlap_groups[0].center, expected_alternate_lmb_center,
            "a forced center must be transformed from the graph LMB into the selected parent LMB after its single probe rotation"
        );

        let support_edge = [1, 2, 3, 4]
            .into_iter()
            .map(EdgeIndex::from)
            .find(|edge| !graph.loop_momentum_basis.loop_edges.contains(edge))
            .unwrap();
        let covariant_thresholds = vec![Esurface {
            energies: vec![
                graph.loop_momentum_basis.loop_edges[LoopIndex(0)],
                support_edge,
            ],
            external_shift: vec![(EdgeIndex::from(0), -1)],
            vertex_set: VertexSet::dummy(),
        }]
        .into();
        let covariant_overlap_input = OverlapInput {
            graph: &graph,
            settings: &settings,
            subspace: &subspace,
            lmbs: &all_lmbs,
            thresholds: &covariant_thresholds,
            edge_masses: overlap_input.edge_masses.clone(),
        };
        let covariant_loop_momenta = LoopMomenta::from_iter([
            ThreeMomentum::new(F(0.3), F(-0.4), F(0.5)),
            ThreeMomentum::new(F(-0.2), F(0.7), F(0.1)),
            ThreeMomentum::new(F(0.6), F(-0.8), F(0.9)),
        ]);
        let covariant_externals = ExternalFourMomenta::from_iter([
            FourMomentum::from_args(F(20.0), F(6.0), F(8.0), F(0.0)),
            FourMomentum::from_args(F(-20.0), F(-6.0), F(-8.0), F(0.0)),
        ]);
        let existing = ti_vec![EsurfaceID::from(0)];
        let identity_center = find_center(
            &covariant_overlap_input,
            &[ExistingEsurfaceId::from(0)],
            &existing,
            &covariant_loop_momenta,
            &covariant_externals,
            false,
        )
        .unwrap();
        assert!(check_global_center(
            &covariant_overlap_input,
            &existing,
            &identity_center,
            &covariant_loop_momenta,
            &covariant_externals,
        ));

        for rotation_method in [
            RotationMethod::Identity,
            RotationMethod::Pi2X,
            RotationMethod::Pi2Y,
            RotationMethod::Pi2Z,
            RotationMethod::EulerAngles(0.1, 0.2, 0.3),
        ] {
            let stability_rotation = Rotation::new(rotation_method);
            let rotated_loop_momenta = covariant_loop_momenta.rotate(&stability_rotation);
            let rotated_externals = covariant_externals
                .iter()
                .map(|momentum| FourMomentum {
                    temporal: momentum.temporal,
                    spatial: momentum.spatial.rotate(&stability_rotation),
                })
                .collect();
            let rotated_center = find_center(
                &covariant_overlap_input,
                &[ExistingEsurfaceId::from(0)],
                &existing,
                &rotated_loop_momenta,
                &rotated_externals,
                false,
            )
            .unwrap();
            let expected_rotated_center = identity_center.rotate(&stability_rotation);

            for (actual, expected) in rotated_center.iter().zip(expected_rotated_center.iter()) {
                for (actual_component, expected_component) in actual.into_iter().zip(expected) {
                    assert!(
                        (actual_component.0 - expected_component.0).abs() < 1.0e-8,
                        "solver-derived center is not covariant under stability rotation {rotation_method}: actual={actual_component}, expected={expected_component}"
                    );
                }
            }
            assert!(
                check_global_center(
                    &covariant_overlap_input,
                    &existing,
                    &rotated_center,
                    &rotated_loop_momenta,
                    &rotated_externals,
                ),
                "solver-derived center is not strictly interior after stability rotation {rotation_method}",
            );
        }
    }
}
