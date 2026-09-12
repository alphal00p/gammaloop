use crate::GammaLoopContext;
use crate::cff::esurface::EsurfaceCollection;
use crate::cff::esurface::EsurfaceID;
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
use crate::utils::Length;
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
#[cfg(test)]
use std::collections::BTreeMap;
use std::fmt::Display;
use typed_index_collections::TiVec;

#[cfg(test)]
#[path = "overlap_subspace_tests.rs"]
mod regression_tests;

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

/// Kinematic data for one threshold instance. Different Cutkosky cuts may have different
/// solved `t*` values and external data even when they share a solve signature.
#[derive(Debug, Clone)]
pub(crate) struct OverlapKinematics {
    pub loop_moms: LoopMomenta<F<f64>>,
    pub external_momenta: ExternalFourMomenta<F<f64>>,
    /// Graph-group members can share routing while carrying different propagator masses.
    pub edge_masses: Option<EdgeVec<F<f64>>>,
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
    subspace: &'a SubspaceData,
    loop_moms: LoopMomenta<F<f64>>,
    external_momenta: ExternalFourMomenta<F<f64>>,
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
        let threshold_subspace = overlap_input.threshold_subspace(surface_id);
        let (surface_loop_moms, surface_external_momenta, edge_masses) = overlap_input
            .surface_kinematics
            .and_then(|kinematics| kinematics.get(surface_id.0))
            .map(|kinematics| {
                (
                    &kinematics.loop_moms,
                    &kinematics.external_momenta,
                    kinematics
                        .edge_masses
                        .as_ref()
                        .unwrap_or(&overlap_input.edge_masses),
                )
            })
            .unwrap_or((loop_moms, external_momenta, &overlap_input.edge_masses));
        let lmb = threshold_subspace.get_lmb(overlap_input.lmbs);

        let mut esurface_constraint_indices: Vec<usize> = Vec::with_capacity(6);

        for edge_id in threshold_subspace.contains(&esurface.energies, overlap_input.graph) {
            if let Some(edge_position) = propagator_constraints.iter().position(|constraint| {
                overlap_input.surface_kinematics.is_none()
                    && *constraint.signature == lmb.edge_signatures[edge_id]
                    && constraint.subspace.solve_signature(overlap_input.lmbs)
                        == threshold_subspace.solve_signature(overlap_input.lmbs)
            }) {
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
                    subspace: threshold_subspace,
                    loop_moms: surface_loop_moms.clone(),
                    external_momenta: surface_external_momenta.clone(),
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
        let threshold_subspace = overlap_input.threshold_subspace(esurface_id);
        let (surface_loop_moms, surface_external_momenta, edge_masses) = overlap_input
            .surface_kinematics
            .and_then(|kinematics| kinematics.get(esurface_id.0))
            .map(|kinematics| {
                (
                    &kinematics.loop_moms,
                    &kinematics.external_momenta,
                    kinematics
                        .edge_masses
                        .as_ref()
                        .unwrap_or(&overlap_input.edge_masses),
                )
            })
            .unwrap_or((loop_moms, external_momenta, &overlap_input.edge_masses));

        let shift_part = esurface.compute_shift_part_from_momenta_in_subspace(
            surface_loop_moms,
            surface_external_momenta,
            threshold_subspace,
            overlap_input.lmbs,
            overlap_input.graph,
            edge_masses,
        );
        b_vector[constraint_index + 1] = -shift_part.0;
        a_matrix[constraint_index + 1][0] = -1.0;
    }

    // propagator constraints
    let mut vertical_offset = esurface_constraints.len() + 1;
    for (cone_index, propagator_constraint) in propagator_constraints.iter().enumerate() {
        a_matrix[vertical_offset][propagator_index_offset + cone_index] = -1.0;
        vertical_offset += 1;

        let spatial_shift = compute_shift_part_subspace(
            &propagator_constraint.signature.internal,
            &propagator_constraint.signature.external,
            &propagator_constraint.loop_moms,
            &propagator_constraint
                .external_momenta
                .iter()
                .map(|p| p.spatial)
                .collect(),
            propagator_constraint.subspace,
        );

        b_vector[vertical_offset] = spatial_shift.px.0;
        b_vector[vertical_offset + 1] = spatial_shift.py.0;
        b_vector[vertical_offset + 2] = spatial_shift.pz.0;

        for (subspace_loop_index, loop_index) in
            overlap_input.subspace.iter_lmb_indices().enumerate()
        {
            if !propagator_constraint
                .subspace
                .contains_loop_index(loop_index)
            {
                continue;
            }
            let individual_loop_signature = propagator_constraint.signature.internal[loop_index];
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

    crate::debug_tags!(#subtraction, #threshold, #overlap, #socp;
        stage = "threshold_socp_input",
        graph = %overlap_input.graph.name,
        surfaces = ?esurfaces_to_consider,
        file.a = ?a_matrix,
        file.b = ?b_vector,
        file.q = ?q_vector,
        file.cones = ?cones,
        "constructed threshold overlap problem"
    );

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
) -> Result<Option<LoopMomenta<F<f64>>>> {
    let mut solver = construct_solver(
        overlap_input,
        esurfaces_to_consider,
        existing_esurfaces,
        loop_moms,
        external_momenta,
        verbose,
    );

    solver.solve();

    crate::debug_tags!(#subtraction, #threshold, #overlap, #socp;
        stage = "threshold_socp_result",
        graph = %overlap_input.graph.name,
        surfaces = ?esurfaces_to_consider,
        status = ?solver.solution.status,
        iterations = solver.solution.iterations,
        primal_objective = solver.solution.obj_val,
        dual_objective = solver.solution.obj_val_dual,
        primal_residual = solver.solution.r_prim,
        dual_residual = solver.solution.r_dual,
        file.solution = ?solver.solution,
        "finished threshold overlap solve"
    );

    let global_loop_number = overlap_input.graph.get_loop_number();
    let esurfaces_to_check = esurfaces_to_consider
        .iter()
        .map(|existing_esurface_id| existing_esurfaces[*existing_esurface_id])
        .collect();

    // Even if optimization did not converge, a physically valid center certifies overlap.
    // Conversely, a failed solve must not silently remove an overlap from the catalogue.
    let center = extract_center(
        global_loop_number,
        overlap_input.subspace,
        &solver.solution.x,
    );
    if check_global_center(
        overlap_input,
        &esurfaces_to_check,
        &center,
        loop_moms,
        external_momenta,
    ) {
        Ok(Some(center))
    } else if solver.solution.status == SolverStatus::PrimalInfeasible {
        Ok(None)
    } else {
        Err(eyre!(
            "Threshold SOCP for graph '{}' and E-surfaces {:?} returned {:?} without a strictly interior center: primal objective={}, dual objective={}, primal residual={}, dual residual={}, iterations={}",
            overlap_input.graph.name,
            esurfaces_to_check,
            solver.solution.status,
            solver.solution.obj_val,
            solver.solution.obj_val_dual,
            solver.solution.r_prim,
            solver.solution.r_dual,
            solver.solution.iterations,
        ))
    }
}

pub(crate) struct OverlapInput<'a> {
    pub graph: &'a Graph,
    pub settings: &'a RuntimeSettings,
    /// Union of the active coordinates used to represent a common center.
    pub subspace: &'a SubspaceData,
    /// Optional per-threshold active coordinates. `None` retains the homogeneous legacy path.
    pub threshold_subspaces: Option<&'a [SubspaceData]>,
    pub lmbs: &'a TiVec<LmbIndex, LoopMomentumBasis>,
    pub thresholds: &'a EsurfaceCollection,
    pub edge_masses: EdgeVec<F<f64>>,
    /// Optional per-threshold kinematics for a multi-cut solve. If absent, the single sample
    /// passed to the overlap routine is used for every threshold.
    pub surface_kinematics: Option<&'a [OverlapKinematics]>,
}

impl OverlapInput<'_> {
    fn threshold_subspace(&self, esurface_id: EsurfaceID) -> &SubspaceData {
        self.threshold_subspaces
            .map_or(self.subspace, |subspaces| &subspaces[esurface_id.0])
    }

    fn validate_subspaces(&self) -> Result<()> {
        let Some(threshold_subspaces) = self.threshold_subspaces else {
            return Ok(());
        };
        if threshold_subspaces.len() != self.thresholds.len() {
            return Err(eyre!(
                "Projected threshold overlap has {} subspaces for {} E-surfaces",
                threshold_subspaces.len(),
                self.thresholds.len(),
            ));
        }
        let common_parent = self.subspace.parent_lmb_index();
        for (esurface_index, subspace) in threshold_subspaces.iter().enumerate() {
            if subspace.parent_lmb_index() != common_parent {
                return Err(eyre!(
                    "Projected threshold E-surface instance {} uses parent LMB {}, while its common center uses parent LMB {}",
                    esurface_index,
                    usize::from(subspace.parent_lmb_index()),
                    usize::from(common_parent),
                ));
            }
            if subspace
                .iter_lmb_indices()
                .any(|loop_index| !self.subspace.contains_loop_index(loop_index))
            {
                return Err(eyre!(
                    "Projected threshold E-surface instance {} contains coordinates outside the common-center subspace",
                    esurface_index,
                ));
            }
        }
        Ok(())
    }
}

pub(crate) fn check_global_center(
    overlap_input: &OverlapInput,
    existing_esurfaces: &ExistingThresholds,
    center: &LoopMomenta<F<f64>>,
    loop_moms: &LoopMomenta<F<f64>>,
    external_momenta: &ExternalFourMomenta<F<f64>>,
) -> bool {
    existing_esurfaces.iter().all(|esurface_id| {
        let esurface = &overlap_input.thresholds[*esurface_id];
        let threshold_subspace = overlap_input.threshold_subspace(*esurface_id);
        let (surface_loop_moms, surface_external_momenta, edge_masses) = overlap_input
            .surface_kinematics
            .and_then(|kinematics| kinematics.get(esurface_id.0))
            .map(|kinematics| {
                (
                    &kinematics.loop_moms,
                    &kinematics.external_momenta,
                    kinematics
                        .edge_masses
                        .as_ref()
                        .unwrap_or(&overlap_input.edge_masses),
                )
            })
            .unwrap_or((loop_moms, external_momenta, &overlap_input.edge_masses));
        let mut center_with_fixed_complement = surface_loop_moms.clone();
        for loop_index in threshold_subspace.iter_lmb_indices() {
            center_with_fixed_complement[loop_index] = center[loop_index];
        }

        let lmb = threshold_subspace.get_lmb(overlap_input.lmbs);

        let esurface_val = esurface.compute_from_momenta(
            lmb,
            edge_masses,
            &center_with_fixed_complement,
            surface_external_momenta,
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
    if let Some(surface_kinematics) = overlap_input.surface_kinematics
        && surface_kinematics.len() != overlap_input.thresholds.len()
    {
        return Err(eyre!(
            "per-surface overlap kinematics has {} entries for {} threshold surfaces",
            surface_kinematics.len(),
            overlap_input.thresholds.len()
        ));
    }
    overlap_input.validate_subspaces()?;
    // An empty catalogue has no automatic partitions, but an explicit center still needs
    // its identity-frame rotation and parent-LMB projection.
    if existing_esurfaces.is_empty()
        && overlap_input
            .settings
            .subtraction
            .overlap_settings
            .force_global_center
            .is_some()
    {
        return find_maximal_overlap_single(
            overlap_input,
            existing_esurfaces,
            loop_moms,
            external_momenta,
            probe_rotation,
        );
    }
    let mut partitions: Vec<(SubspaceData, ExistingThresholds)> = Vec::new();
    for existing_id in existing_esurfaces.iter_enumerated().map(|(id, _)| id) {
        let surface_id = existing_esurfaces[existing_id];
        let threshold_subspace = overlap_input.threshold_subspace(surface_id);
        if let Some((_, members)) = partitions.iter_mut().find(|(representative, _)| {
            representative.solve_signature(overlap_input.lmbs)
                == threshold_subspace.solve_signature(overlap_input.lmbs)
        }) {
            members.push(surface_id);
        } else {
            let mut members = ExistingThresholds::new();
            members.push(surface_id);
            partitions.push((threshold_subspace.clone(), members));
        }
    }

    let mut combined = OverlapStructure {
        overlap_groups: Vec::new(),
        existing_esurfaces: existing_esurfaces.clone(),
    };
    let mut original_ids_by_partition = Vec::new();
    for (_, members) in &partitions {
        original_ids_by_partition.push(
            members
                .iter()
                .map(|surface_id| {
                    existing_esurfaces
                        .iter_enumerated()
                        .find_map(|(id, candidate)| (candidate == surface_id).then_some(id))
                        .expect("partition surface must belong to original threshold set")
                })
                .collect::<Vec<_>>(),
        );
    }
    for (partition_index, (representative, members)) in partitions.into_iter().enumerate() {
        let partition_input = OverlapInput {
            graph: overlap_input.graph,
            settings: overlap_input.settings,
            subspace: &representative,
            threshold_subspaces: None,
            lmbs: overlap_input.lmbs,
            thresholds: overlap_input.thresholds,
            edge_masses: overlap_input.edge_masses.clone(),
            surface_kinematics: overlap_input.surface_kinematics,
        };
        let local = find_maximal_overlap_single(
            &partition_input,
            &members,
            loop_moms,
            external_momenta,
            probe_rotation,
        )?;
        for mut group in local.overlap_groups {
            for id in group
                .existing_esurfaces
                .iter_mut()
                .chain(&mut group.complement)
            {
                *id = original_ids_by_partition[partition_index][usize::from(*id)];
            }
            combined.overlap_groups.push(group);
        }
    }
    Ok(combined)
}

/// Variant of [`find_maximal_overlap`] that keeps explicitly tagged thresholds in separate
/// solve problems. Unlabelled thresholds still use the automatic exact-subspace partition.
#[cfg(test)]
pub(crate) fn find_maximal_overlap_with_group_ids(
    overlap_input: &OverlapInput,
    existing_esurfaces: &ExistingThresholds,
    group_ids: &[Option<usize>],
    loop_moms: &LoopMomenta<F<f64>>,
    external_momenta: &ExternalFourMomenta<F<f64>>,
    probe_rotation: &Rotation,
) -> Result<OverlapStructure> {
    if let Some(surface_kinematics) = overlap_input.surface_kinematics
        && surface_kinematics.len() != overlap_input.thresholds.len()
    {
        return Err(eyre!(
            "per-surface overlap kinematics has {} entries for {} threshold surfaces",
            surface_kinematics.len(),
            overlap_input.thresholds.len()
        ));
    }
    if group_ids.len() < overlap_input.thresholds.len() {
        return Err(eyre!(
            "threshold solve-group metadata has {} entries for {} threshold surfaces",
            group_ids.len(),
            overlap_input.thresholds.len()
        ));
    }
    overlap_input.validate_subspaces()?;
    let mut labelled = BTreeMap::<Option<usize>, ExistingThresholds>::new();
    for existing_id in existing_esurfaces.iter_enumerated().map(|(id, _)| id) {
        let surface_id = existing_esurfaces[existing_id];
        let label = group_ids[surface_id.0];
        labelled.entry(label).or_default().push(surface_id);
    }
    for (label, members) in &labelled {
        let Some(group_id) = label else {
            continue;
        };
        let Some(first_surface_id) = members.first() else {
            continue;
        };
        let first_subspace = overlap_input.threshold_subspace(*first_surface_id);
        let first_signature = first_subspace.solve_signature(overlap_input.lmbs);
        let incompatible = members.iter().skip(1).any(|surface_id| {
            overlap_input
                .threshold_subspace(*surface_id)
                .solve_signature(overlap_input.lmbs)
                != first_signature
        });
        if incompatible {
            let members = members
                .iter()
                .map(|surface_id| {
                    let subspace = overlap_input.threshold_subspace(*surface_id);
                    format!(
                        "E-surface instance {}: parent {:?}, selected cycles {:?}",
                        surface_id.0,
                        subspace.get_lmb(overlap_input.lmbs).loop_edges,
                        subspace.solve_signature(overlap_input.lmbs)
                    )
                })
                .join("; ");
            return Err(eyre!(
                "explicit threshold group_id={} contains incompatible solve subspaces: {}",
                group_id,
                members,
            ));
        }
    }
    let mut combined = OverlapStructure {
        overlap_groups: Vec::new(),
        existing_esurfaces: existing_esurfaces.clone(),
    };
    for members in labelled.into_values() {
        let local = find_maximal_overlap(
            overlap_input,
            &members,
            loop_moms,
            external_momenta,
            probe_rotation,
        )?;
        for mut group in local.overlap_groups {
            for id in group
                .existing_esurfaces
                .iter_mut()
                .chain(&mut group.complement)
            {
                let surface_id = local.existing_esurfaces[*id];
                *id = existing_esurfaces
                    .iter_enumerated()
                    .find_map(|(original_id, candidate)| {
                        (candidate == &surface_id).then_some(original_id)
                    })
                    .expect("labelled overlap surface must belong to original threshold set");
            }
            combined.overlap_groups.push(group);
        }
    }
    Ok(combined)
}

fn find_maximal_overlap_single(
    overlap_input: &OverlapInput,
    existing_esurfaces: &ExistingThresholds,
    loop_moms: &LoopMomenta<F<f64>>,
    external_momenta: &ExternalFourMomenta<F<f64>>,
    probe_rotation: &Rotation,
) -> Result<OverlapStructure> {
    overlap_input.validate_subspaces()?;
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
        let center_externals = existing_esurfaces
            .first()
            .and_then(|id| {
                overlap_input
                    .surface_kinematics
                    .map(|kinematics| &kinematics[id.0].external_momenta)
            })
            .unwrap_or(external_momenta);
        let rotated_external_spatial = center_externals
            .iter()
            .map(|momentum| momentum.spatial)
            .collect();
        // Forced centers are configured in the graph's identity-frame LMB. Apply the probe
        // rotation once, then express that rotated point in the cut-side LMB used by the solver.
        let mut global_center_probe = global_center_identity.rotate(probe_rotation).lmb_transform(
            &overlap_input.graph.loop_momentum_basis,
            overlap_input.subspace.get_lmb(overlap_input.lmbs),
            &rotated_external_spatial,
        );
        if let Some(surface_kinematics) = overlap_input.surface_kinematics {
            for surface_id in existing_esurfaces {
                let surface_externals = &surface_kinematics[surface_id.0].external_momenta;
                if surface_externals == center_externals {
                    continue;
                }
                let surface_center = global_center_identity.rotate(probe_rotation).lmb_transform(
                    &overlap_input.graph.loop_momentum_basis,
                    overlap_input.subspace.get_lmb(overlap_input.lmbs),
                    &surface_externals
                        .iter()
                        .map(|momentum| momentum.spatial)
                        .collect(),
                );
                if overlap_input
                    .subspace
                    .iter_lmb_indices()
                    .any(|index| surface_center[index] != global_center_probe[index])
                {
                    return Err(eyre!(
                        "forced global threshold center gives incompatible active coordinates for E-surface instance {}: selected edges {:?}, center {}, expected {}",
                        surface_id.0,
                        overlap_input
                            .subspace
                            .iter_basis_edges(overlap_input.lmbs)
                            .collect_vec(),
                        surface_center,
                        global_center_probe,
                    ));
                }
            }
        }
        // A subspace center specifies only its active coordinates. The complementary loop
        // coordinates remain fixed at the sampled point throughout center validation, radial
        // solving, and final CT reconstruction. Affine LMB transforms can populate inactive
        // components even when the identity-frame center is zero, so project them out explicitly.
        for loop_index in (0..global_center_probe.len()).map(LoopIndex::from) {
            if !overlap_input.subspace.contains_loop_index(loop_index) {
                global_center_probe[loop_index] = ThreeMomentum::new(F(0.0), F(0.0), F(0.0));
            }
        }

        tracing::debug!(
            graph = %overlap_input.graph.name,
            rotation_id = %probe_rotation.method,
            center_provenance = "forced_identity_frame_rotated_lmb_transformed_and_projected_once",
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
    )?;

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
    )?;

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
            )?
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
            )?;

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
    ) -> Result<Self> {
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
                )?;

                if let Some(center) = center {
                    res.insert((esurface_id_1, esurface_id_2), center);
                    res.has_pair_with[Into::<usize>::into(esurface_id_1)].push(esurface_id_2);
                    res.has_pair_with[Into::<usize>::into(esurface_id_2)].push(esurface_id_1);
                }
            }
        }

        Ok(res)
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
        // A maximal overlap can have every pair covered by different larger overlaps.
        // Only containment of the complete candidate permits skipping its SOCP problem.
        for (first, _) in existing_esurfaces.iter_enumerated() {
            for remaining in self.has_pair_with[usize::from(first)]
                .iter()
                .copied()
                .filter(|&id| id > first)
                .combinations(subset_len - 1)
            {
                if !remaining
                    .iter()
                    .tuple_combinations()
                    .all(|(&left, &right)| self.pair_exists((left, right)))
                {
                    continue;
                }
                let mut candidate = vec![first];
                candidate.extend(remaining);
                candidate.sort_unstable();
                if !is_subset_of_result(&candidate, result) {
                    res.insert(candidate);
                }
            }
        }
        res
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        cff::{
            VertexSet,
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
    fn maximal_overlap_candidates_include_cliques_whose_pairs_are_already_covered() {
        let existing: ExistingThresholds = (0..9).map(EsurfaceID::from).collect();
        let center = LoopMomenta::from_iter([ThreeMomentum::new(F(0.0), F(0.0), F(0.0))]);
        // X,Y,Z have a common region. Two small balls in each pair's exclusive region
        // produce the larger overlaps XYuv, XZwx, YZyz, while XYZ remains maximal too.
        let larger = [[0, 1, 3, 4], [0, 2, 5, 6], [1, 2, 7, 8]];
        let result = OverlapStructure {
            existing_esurfaces: existing.clone(),
            overlap_groups: larger
                .map(|members| OverlapGroup {
                    existing_esurfaces: members.into_iter().map(ExistingEsurfaceId::from).collect(),
                    complement: vec![],
                    center: center.clone(),
                })
                .to_vec(),
        };
        let mut pairs = EsurfacePairs::new_empty(existing.len());
        for members in larger {
            for (left, right) in members
                .into_iter()
                .map(ExistingEsurfaceId::from)
                .tuple_combinations()
            {
                if !pairs.pair_exists((left, right)) {
                    pairs.insert((left, right), center.clone());
                    pairs.has_pair_with[usize::from(left)].push(right);
                    pairs.has_pair_with[usize::from(right)].push(left);
                }
            }
        }
        let central = (0..3).map(ExistingEsurfaceId::from).collect_vec();
        assert!(
            central
                .iter()
                .tuple_combinations()
                .all(|(&left, &right)| { is_subset_of_result(&[left, right], &result) })
        );
        assert_eq!(
            pairs.construct_possible_subsets_of_len(&existing, 3, &result),
            HashSet::from_iter([central]),
        );
    }

    #[test]
    fn overlap_partitions_remap_complements_and_keep_instance_kinematics() {
        test_initialise().unwrap();
        let graph: Graph = dot!(digraph partitioned_centers {
            ext [style=invis]
            edge [num=1 mass=0]
            node [num=1]
            ext->a:0 [id=0]
            a->b [id=1]
            b->a [id=2]
            a->b [id=3]
            b->a [id=4]
            ext->b:1 [id=5]
        })
        .unwrap();
        let lmbs = ti_vec![graph.loop_momentum_basis.clone()];
        let subspaces = [LoopIndex(0), LoopIndex(1)].map(|index| {
            SubspaceData::new_from_parent_basis_edges(
                &[graph.loop_momentum_basis.loop_edges[index]],
                &graph.full_filter(),
                LmbIndex::from(0),
                &graph,
                &lmbs,
            )
            .unwrap()
        });
        let common = SubspaceData::union_in_common_parent(subspaces.iter(), &graph, &lmbs).unwrap();
        let support = graph
            .iter_loop_edges()
            .map(|(_, edge, _)| edge)
            .find(|edge| !graph.loop_momentum_basis.loop_edges.contains(edge))
            .unwrap();
        let external =
            graph.loop_momentum_basis.ext_edges[crate::momentum::sample::ExternalIndex(0)];
        let thresholds: EsurfaceCollection = (0..4)
            .map(|_| Esurface {
                energies: vec![support],
                external_shift: vec![(external, -1)],
                vertex_set: VertexSet::dummy(),
            })
            .collect::<Vec<_>>()
            .into();
        let threshold_subspaces = [0, 1, 0, 1].map(|index| subspaces[index].clone());
        let masses = graph.underlying.new_edgevec(|_, _, _| F(0.0));
        let externals: ExternalFourMomenta<_> = [
            FourMomentum::from_args(F(1.0), F(0.0), F(0.0), F(0.0)),
            FourMomentum::from_args(F(-1.0), F(0.0), F(0.0), F(0.0)),
        ]
        .into_iter()
        .collect();
        let mut kinematics = (0..4)
            .map(|index| {
                let mut loop_moms = LoopMomenta::from_iter(
                    (0..3).map(|_| ThreeMomentum::new(F(0.0), F(0.0), F(0.0))),
                );
                loop_moms[LoopIndex(2)].px = F(10.0 * (index + 1) as f64);
                OverlapKinematics {
                    loop_moms,
                    external_momenta: externals.clone(),
                    edge_masses: None,
                }
            })
            .collect_vec();
        let mut member_masses = masses.clone();
        member_masses[support] = F(0.5);
        kinematics[0].edge_masses = Some(member_masses);
        let settings = RuntimeSettings::default();
        let input = OverlapInput {
            graph: &graph,
            settings: &settings,
            subspace: &common,
            threshold_subspaces: Some(&threshold_subspaces),
            lmbs: &lmbs,
            thresholds: &thresholds,
            edge_masses: masses.clone(),
            surface_kinematics: Some(&kinematics),
        };
        let existing: ExistingThresholds = thresholds.keys().collect();
        let rotation = Rotation::new(RotationMethod::Identity);
        let automatic = find_maximal_overlap(
            &input,
            &existing,
            &kinematics[0].loop_moms,
            &externals,
            &rotation,
        )
        .unwrap();
        let labelled = find_maximal_overlap_with_group_ids(
            &input,
            &existing,
            &[Some(0), Some(1), Some(0), Some(1)],
            &kinematics[0].loop_moms,
            &externals,
            &rotation,
        )
        .unwrap();
        let conflict = find_maximal_overlap_with_group_ids(
            &input,
            &existing,
            &[Some(0); 4],
            &kinematics[0].loop_moms,
            &externals,
            &rotation,
        )
        .unwrap_err()
        .to_string();
        assert!(conflict.contains("group_id=0"));
        for member in 0..4 {
            assert!(conflict.contains(&format!("E-surface instance {member}:")));
        }
        for overlap in [&automatic, &labelled] {
            assert_eq!(overlap.overlap_groups.len(), 4);
            for group in &overlap.overlap_groups {
                assert_eq!(group.existing_esurfaces.len(), 1);
                let member = usize::from(group.existing_esurfaces[0]);
                assert_eq!(
                    group.complement,
                    vec![ExistingEsurfaceId::from((member + 2) % 4)]
                );
                assert!(check_global_center(
                    &input,
                    &ti_vec![EsurfaceID(member)],
                    &group.center,
                    &kinematics[member].loop_moms,
                    &externals
                ));
                assert!(group.center.hyper_radius_squared(None) > F(1.0));
            }
        }
        let first = automatic
            .overlap_groups
            .iter()
            .find(|group| group.existing_esurfaces == vec![ExistingEsurfaceId::from(0)])
            .unwrap();
        let mut outside_massive_surface = first.center.clone();
        outside_massive_surface[LoopIndex(0)].px += F(0.95);
        assert!(!check_global_center(
            &input,
            &ti_vec![EsurfaceID(0)],
            &outside_massive_surface,
            &kinematics[0].loop_moms,
            &externals
        ));

        let mut forced_settings = settings.clone();
        forced_settings
            .subtraction
            .overlap_settings
            .force_global_center = Some(
            first
                .center
                .iter()
                .map(|momentum| [momentum.px.0, momentum.py.0, momentum.pz.0])
                .collect(),
        );
        // Other independent groups may have different external data; they must not participate
        // in validating this group's configured center.
        kinematics[1].external_momenta[crate::momentum::sample::ExternalIndex(0)]
            .temporal
            .value = F(2.0);
        let forced_input = OverlapInput {
            graph: &graph,
            settings: &forced_settings,
            subspace: &common,
            threshold_subspaces: Some(&threshold_subspaces),
            lmbs: &lmbs,
            thresholds: &thresholds,
            edge_masses: masses,
            surface_kinematics: Some(&kinematics),
        };
        assert!(
            find_maximal_overlap(
                &forced_input,
                &ti_vec![EsurfaceID(0)],
                &kinematics[0].loop_moms,
                &externals,
                &rotation
            )
            .is_ok()
        );
    }

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
            subgraph.clone(),
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
            threshold_subspaces: None,
            lmbs: &all_lmbs,
            thresholds: &thresholds,
            edge_masses: masses,
            surface_kinematics: None,
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

        let full_parent_subspace = SubspaceData::new_with_user_selected_lmb(
            graph.full_filter(),
            LmbIndex::from(0),
            &graph,
            &all_lmbs,
        )
        .unwrap();
        let common_subspace = SubspaceData::union_in_common_parent(
            [&subspace, &full_parent_subspace],
            &graph,
            &all_lmbs,
        )
        .unwrap();
        // The same supergraph geometry is deliberately presented twice as distinct projected
        // E-surface instances, each retaining its own fixed complement.
        let projected_thresholds: crate::cff::esurface::EsurfaceCollection = vec![
            thresholds[EsurfaceID::from(0)].clone(),
            thresholds[EsurfaceID::from(0)].clone(),
        ]
        .into();
        let projected_subspaces: TiVec<EsurfaceID, _> =
            ti_vec![subspace.clone(), full_parent_subspace];
        let projected_input = OverlapInput {
            graph: &graph,
            settings: &settings,
            subspace: &common_subspace,
            threshold_subspaces: Some(&projected_subspaces.raw),
            lmbs: &all_lmbs,
            thresholds: &projected_thresholds,
            edge_masses: overlap_input.edge_masses.clone(),
            surface_kinematics: overlap_input.surface_kinematics,
        };
        assert!(!check_global_center(
            &projected_input,
            &ti_vec![EsurfaceID::from(0)],
            &center,
            &sampled_momenta,
            &external_momenta,
        ));
        assert!(check_global_center(
            &projected_input,
            &ti_vec![EsurfaceID::from(1)],
            &center,
            &sampled_momenta,
            &external_momenta,
        ));

        let forced_center_coordinates = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let identity_frame_center =
            LoopMomenta::from_iter(forced_center_coordinates.iter().map(|coordinates| {
                ThreeMomentum::new(F(coordinates[0]), F(coordinates[1]), F(coordinates[2]))
            }));
        let probe_rotation = Rotation::new(RotationMethod::Pi2Z);
        let expected_probe_center = identity_frame_center.rotate(&probe_rotation);
        let mut expected_subspace_probe_center = expected_probe_center.clone();
        for loop_index in (0..expected_subspace_probe_center.len()).map(LoopIndex::from) {
            if !subspace.contains_loop_index(loop_index) {
                expected_subspace_probe_center[loop_index] =
                    ThreeMomentum::new(F(0.0), F(0.0), F(0.0));
            }
        }
        let mut forced_settings = RuntimeSettings::default();
        forced_settings
            .subtraction
            .overlap_settings
            .force_global_center = Some(forced_center_coordinates.clone());
        let forced_overlap_input = OverlapInput {
            graph: &graph,
            settings: &forced_settings,
            subspace: &subspace,
            threshold_subspaces: None,
            lmbs: &all_lmbs,
            thresholds: &thresholds,
            edge_masses: overlap_input.edge_masses.clone(),
            surface_kinematics: overlap_input.surface_kinematics,
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
            forced_overlap.overlap_groups[0].center, expected_subspace_probe_center,
            "an identity-frame forced center must be rotated exactly once and projected onto the active subspace"
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
            threshold_subspaces: None,
            lmbs: &alternate_lmbs,
            thresholds: &empty_thresholds,
            edge_masses: overlap_input.edge_masses.clone(),
            surface_kinematics: overlap_input.surface_kinematics,
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

        let proper_alternate_subspace = SubspaceData::new_with_user_selected_lmb(
            subgraph,
            LmbIndex::from(1),
            &graph,
            &alternate_lmbs,
        )
        .expect("the non-default parent LMB must support the same proper two-loop subspace");
        assert_eq!(proper_alternate_subspace.loopcount(), 2);
        let affine_external_momenta = ExternalFourMomenta::from_iter([
            FourMomentum::from_args(F(5.0), F(1.0), F(2.0), F(3.0)),
            FourMomentum::from_args(F(-5.0), F(-1.0), F(-2.0), F(-3.0)),
        ]);
        let affine_external_spatial = affine_external_momenta
            .iter()
            .map(|momentum| momentum.spatial)
            .collect();
        let unprojected_affine_center = identity_frame_center.lmb_transform(
            &graph.loop_momentum_basis,
            proper_alternate_subspace.get_lmb(&alternate_lmbs),
            &affine_external_spatial,
        );
        assert!(
            unprojected_affine_center
                .iter_enumerated()
                .any(|(loop_index, momentum)| {
                    !proper_alternate_subspace.contains_loop_index(loop_index)
                        && momentum.norm_squared() > F(0.0)
                }),
            "the fixture must exercise an affine LMB transform with a nonzero inactive component"
        );
        let mut expected_projected_affine_center = unprojected_affine_center.clone();
        for loop_index in (0..expected_projected_affine_center.len()).map(LoopIndex::from) {
            if !proper_alternate_subspace.contains_loop_index(loop_index) {
                expected_projected_affine_center[loop_index] =
                    ThreeMomentum::new(F(0.0), F(0.0), F(0.0));
            }
        }
        let mut affine_forced_settings = RuntimeSettings::default();
        affine_forced_settings
            .subtraction
            .overlap_settings
            .force_global_center = Some(forced_center_coordinates);
        let affine_overlap_input = OverlapInput {
            graph: &graph,
            settings: &affine_forced_settings,
            subspace: &proper_alternate_subspace,
            threshold_subspaces: None,
            lmbs: &alternate_lmbs,
            thresholds: &empty_thresholds,
            edge_masses: overlap_input.edge_masses.clone(),
            surface_kinematics: overlap_input.surface_kinematics,
        };
        let affine_overlap = find_maximal_overlap(
            &affine_overlap_input,
            &ti_vec![],
            &sampled_momenta,
            &affine_external_momenta,
            &Rotation::new(RotationMethod::Identity),
        )
        .unwrap();
        assert_eq!(
            affine_overlap.overlap_groups[0].center, expected_projected_affine_center,
            "an affine parent-LMB transform must not displace the fixed complement of a forced subspace center"
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
            threshold_subspaces: None,
            lmbs: &all_lmbs,
            thresholds: &covariant_thresholds,
            edge_masses: overlap_input.edge_masses.clone(),
            surface_kinematics: overlap_input.surface_kinematics,
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
        .unwrap()
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
            .unwrap()
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
