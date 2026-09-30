use crate::GammaLoopContext;
use crate::cff::esurface::EsurfaceCollection;
use crate::cff::esurface::EsurfaceID;
use crate::cff::esurface::ExistingEsurfaceId;
use crate::cff::esurface::ExistingEsurfaces;
use crate::cff::esurface::GroupEsurfaceId;
use crate::cff::esurface::RaisedEsurfaceData;
use crate::cff::esurface::RaisedEsurfaceId;
use crate::cff::esurface::{esurface_value_is_strictly_inside, get_representative};
use crate::graph::GraphGroupPosition;
use crate::graph::LoopMomentumBasis;
use crate::integrands::process::GenericEvaluator;
use crate::momentum::ThreeMomentum;
use crate::momentum::sample::ExternalFourMomenta;
use crate::momentum::sample::LoopMomenta;
use crate::momentum::signature::LoopExtSignature;
use crate::processes::EvaluatorSettings;
use crate::settings::RuntimeSettings;
use crate::settings::runtime::OverlapCenterObjective;
use crate::utils::F;
use crate::utils::GS;
use crate::utils::compute_shift_part;
use crate::utils::hyperdual_utils::simple_n_deriv_shape;
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
use std::cell::RefCell;
use std::fmt::Display;
use symbolica::atom::Atom;
use symbolica::atom::AtomCore;
use symbolica::evaluate::FunctionMap;
use symbolica::evaluate::OptimizationSettings;
use symbolica::function;
use typed_index_collections::TiVec;

#[derive(Debug, Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct OverlapGroup {
    pub existing_esurfaces: Vec<ExistingEsurfaceId>,
    pub complement: Vec<ExistingEsurfaceId>,
    /// Amplitude overlap centers are stored in the identity-probe frame.
    /// Counterterm evaluation rotates them exactly once into the current probe frame.
    pub center: LoopMomenta<F<f64>>,
    pub prefactor_evaluator: Option<Vec<RefCell<GenericEvaluator>>>,
}

#[derive(Debug, Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct OverlapStructure {
    pub overlap_groups: Vec<OverlapGroup>,
    pub existing_esurfaces: ExistingEsurfaces,
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

            writeln!(f, "center:\n {}", group.center)?;
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

    pub fn build_evaluators(
        &mut self,
        atoms: &TiVec<GroupEsurfaceId, Atom>,
        optimization_settings: &OptimizationSettings,
        num_loops: usize,
        num_externals: usize,
        model_params: Vec<Atom>,
        power: i32,
    ) -> Result<()> {
        let group_square_atoms = self
            .overlap_groups
            .iter()
            .map(|group| {
                group
                    .complement
                    .iter()
                    .map(|&existing_esurface_id| {
                        let esurface = self.existing_esurfaces[existing_esurface_id];
                        let atom = &atoms[esurface];
                        atom.pow(power)
                    })
                    .reduce(|prod, atom| prod * atom)
                    .unwrap_or_else(|| Atom::num(1))
            })
            .collect_vec();

        let denominator = group_square_atoms
            .iter()
            .fold(Atom::new(), |sum, atom| sum + atom);

        let params = (0..num_loops)
            .flat_map(|loop_index| {
                (1..=3).map(move |spatial_index| function!(GS.loop_mom, loop_index, spatial_index))
            })
            .chain((0..num_externals).flat_map(|external_index| {
                (0..=3).map(move |spatial_index| {
                    function!(GS.external_mom, external_index, spatial_index)
                })
            }))
            .chain(model_params)
            .collect_vec();

        for (group, square_atom) in self.overlap_groups.iter_mut().zip(group_square_atoms) {
            let atom = square_atom / &denominator;
            let num_orders = power.saturating_sub(1).max(1) as usize;

            let evaluators = (0..num_orders)
                .map(|order_index| {
                    GenericEvaluator::new_from_raw_params(
                        [atom.clone()],
                        &params,
                        &FunctionMap::new(),
                        vec![],
                        optimization_settings.clone(),
                        (order_index > 0)
                            .then(|| simple_n_deriv_shape(order_index))
                            .map(|shape| (shape, Vec::new())),
                        &EvaluatorSettings::default(),
                    )
                    .map(RefCell::new)
                })
                .collect::<Result<Vec<_>>>()?;

            group.prefactor_evaluator = Some(evaluators);
        }

        Ok(())
    }

    pub fn new_empty() -> Self {
        Self {
            overlap_groups: vec![],
            existing_esurfaces: ExistingEsurfaces::new(),
        }
    }

    pub(crate) fn localized_to_existing_surfaces(
        &self,
        local_esurface_exists: &TiVec<GroupEsurfaceId, bool>,
    ) -> Self {
        let mut remapped_existing_esurfaces = ExistingEsurfaces::new();
        let mut existing_esurface_map: TiVec<ExistingEsurfaceId, Option<ExistingEsurfaceId>> =
            TiVec::with_capacity(self.existing_esurfaces.len());

        for &group_esurface_id in self.existing_esurfaces.iter() {
            if local_esurface_exists[group_esurface_id] {
                let remapped_existing_esurface_id =
                    ExistingEsurfaceId::from(remapped_existing_esurfaces.len());
                remapped_existing_esurfaces.push(group_esurface_id);
                existing_esurface_map.push(Some(remapped_existing_esurface_id));
            } else {
                existing_esurface_map.push(None);
            }
        }

        let mut localized_overlap_groups = Vec::with_capacity(self.overlap_groups.len());
        for overlap_group in &self.overlap_groups {
            let mut localized_existing_esurfaces = overlap_group
                .existing_esurfaces
                .iter()
                .filter_map(|&existing_esurface_id| existing_esurface_map[existing_esurface_id])
                .collect_vec();
            localized_existing_esurfaces.sort_unstable();
            localized_existing_esurfaces.dedup();

            if localized_existing_esurfaces.is_empty()
                || localized_overlap_groups.iter().any(|group: &OverlapGroup| {
                    group.existing_esurfaces == localized_existing_esurfaces
                })
            {
                continue;
            }

            localized_overlap_groups.push(OverlapGroup {
                existing_esurfaces: localized_existing_esurfaces,
                complement: vec![],
                center: overlap_group.center.clone(),
                prefactor_evaluator: None,
            });
        }

        let mut localized = Self {
            overlap_groups: localized_overlap_groups,
            existing_esurfaces: remapped_existing_esurfaces,
        };
        localized.fill_in_complements();
        localized
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

fn extract_center(num_loops: usize, solution: &[f64]) -> LoopMomenta<F<f64>> {
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

impl OverlapCenterObjective {
    /// Refine a final overlap group's certified witness without changing its membership.
    /// The builders retain their own geometry; this policy bounds optional work and keeps
    /// a certified center whenever construction, optimization, or physical checks fail.
    pub(crate) fn refine_center(
        self,
        center: &mut LoopMomenta<F<f64>>,
        e_cm: f64,
        mut construct: impl FnMut(Self, f64) -> Result<DefaultSolver>,
        extract: impl Fn(&[f64]) -> LoopMomenta<F<f64>>,
        certify: impl Fn(&LoopMomenta<F<f64>>) -> Option<f64>,
    ) {
        if self == Self::MaxMinDepth {
            return;
        }
        let Some(mut clearance) = certify(center) else {
            tracing::error!("optional overlap refinement received an uncertified witness");
            return;
        };
        let mut minimum_radius = 0.0;
        let mut required_clearance = 0.0;
        for objective in [Self::RelaxedChebyshev, Self::MinSum] {
            let mut solver = match construct(objective, minimum_radius) {
                Ok(solver) => solver,
                Err(error) => {
                    crate::debug_tags!(#subtraction, #threshold, #overlap, #socp;
                        stage = "final_overlap_center_refinement",
                        objective = ?objective,
                        accepted = false,
                        error = %error,
                        "overlap refinement construction failed; retaining certified center"
                    );
                    return;
                }
            };
            solver.solve();
            let candidate = extract(&solver.solution.x);
            let candidate_clearance = certify(&candidate);
            let accepted = solver.solution.status == SolverStatus::Solved
                && candidate_clearance.is_some_and(|radius| radius >= required_clearance);
            crate::debug_tags!(#subtraction, #threshold, #overlap, #socp;
                stage = "final_overlap_center_refinement",
                objective = ?objective,
                status = ?solver.solution.status,
                iterations = solver.solution.iterations,
                previous_clearance = clearance,
                candidate_clearance = ?candidate_clearance,
                minimum_radius,
                required_clearance,
                accepted,
                "optional final-center optimization; rejected candidates retain the certified witness"
            );
            if !accepted {
                return;
            }
            let candidate_clearance = candidate_clearance.expect("accepted physical center");
            if objective == Self::MinSum || candidate_clearance > clearance {
                *center = candidate;
                clearance = candidate_clearance;
            }
            if self == Self::RelaxedChebyshev || !clearance.is_finite() {
                return;
            }
            // Interpret the solver's tolerance constants at the physical energy scale,
            // so the allowed loss of geometric clearance covaries with the input units.
            // The physical certificate remains authoritative even if the solver's own
            // unscaled absolute stopping criterion is weaker for a tiny overlap.
            let accuracy = (solver.settings().tol_gap_abs + solver.settings().tol_feas)
                * e_cm.abs()
                + solver.settings().tol_gap_rel * clearance.abs();
            let allowance = accuracy.min(0.5 * clearance);
            required_clearance = clearance - allowance;
            // Leave half the allowed accuracy for the second solver's constraint residual;
            // the physical certificate above enforces the full promised radius floor.
            minimum_radius = clearance - 0.5 * allowance;
        }
    }
}

fn construct_solver(
    overlap_input: &OverlapInput,
    esurfaces_to_consider: &[ExistingEsurfaceId],
    existing_esurfaces: &ExistingEsurfaces,
    external_momenta: &ExternalFourMomenta<F<f64>>,
    verbose: bool,
    objective: OverlapCenterObjective,
    minimum_radius: f64,
) -> Result<DefaultSolver> {
    let num_loops = overlap_input
        .graph_data
        .first()
        .expect("no graphs passed to overlap")
        .lmb
        .loop_edges
        .len();

    let num_loop_vars = 3 * num_loops;

    // first we study the structure of the problem
    let mut propagator_constraints: Vec<PropagatorConstraint> = Vec::with_capacity(20);

    let mut inequivalent_masses: Vec<F<f64>> = vec![];

    let mut esurface_constraints: Vec<Vec<usize>> = Vec::with_capacity(esurfaces_to_consider.len());
    let local_esurfaces_to_consider = esurfaces_to_consider
        .iter()
        .flat_map(|existing_esurface_id| {
            let group_esurface_id = existing_esurfaces[*existing_esurface_id];
            overlap_input.group_esurface_map[group_esurface_id]
                .iter_enumerated()
                .filter_map(move |(graph_group_pos, option_raised_esurface_id)| {
                    option_raised_esurface_id.and_then(|raised_esurface_id| {
                        overlap_input.local_esurface_exists[graph_group_pos][group_esurface_id]
                            .then_some((*existing_esurface_id, graph_group_pos, raised_esurface_id))
                    })
                })
        })
        .collect_vec();

    for (_, graph_group_pos, raised_esurface_id) in local_esurfaces_to_consider.iter().copied() {
        let esurface_id = representative_local_esurface_id(
            &overlap_input.graph_data[graph_group_pos],
            raised_esurface_id,
        );

        let esurface = &overlap_input.graph_data[graph_group_pos].esurfaces[esurface_id];
        let lmb = overlap_input.graph_data[graph_group_pos].lmb;
        let edge_masses = &overlap_input.graph_data[graph_group_pos].edge_masses;

        let mut esurface_constraint_indices: Vec<usize> = Vec::with_capacity(6);

        for &edge_id in &esurface.energies {
            if let Some(edge_position) = propagator_constraints.iter().position(|constraint| {
                *constraint.signature == lmb.edge_signatures[edge_id]
                    && constraint
                        .mass_pointer
                        .map_or(F(0.0), |index| inequivalent_masses[index])
                        == edge_masses[edge_id]
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
    if objective == OverlapCenterObjective::MinSum {
        // Group IDs identify inequivalent physical surfaces across graph copies. Keep every
        // local feasibility row but count each physical surface only once in the objective.
        // Repeated energy occurrences WITHIN one surface still contribute their multiplicity.
        let mut objective_surfaces = HashSet::default();
        for (row, (surface, _, _)) in esurface_constraints
            .iter()
            .zip(&local_esurfaces_to_consider)
        {
            if !objective_surfaces.insert(*surface) {
                continue;
            }
            for index in row {
                q_vector[propagator_index_offset + index] += 1.0;
            }
        }
    } else {
        q_vector[0] = 1.0;
    }

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
    b_vector[0] = -minimum_radius;
    // esurface constraints
    for (constraint_index, ((_, graph_group_pos, raised_esurface_id), esurface_constraint)) in
        local_esurfaces_to_consider
            .iter()
            .zip(esurface_constraints.iter())
            .enumerate()
    {
        for prop_index in esurface_constraint {
            a_matrix[constraint_index + 1][*prop_index + propagator_index_offset] += 1.0;
        }
        let esurface_id = representative_local_esurface_id(
            &overlap_input.graph_data[*graph_group_pos],
            *raised_esurface_id,
        );
        let lmb = overlap_input.graph_data[*graph_group_pos].lmb;
        let esurface = &overlap_input.graph_data[*graph_group_pos].esurfaces[esurface_id];

        let shift_part = esurface.compute_shift_part_from_momenta(external_momenta, lmb);
        b_vector[constraint_index + 1] = -shift_part.0;
        a_matrix[constraint_index + 1][0] = if objective != OverlapCenterObjective::MaxMinDepth {
            -overlap_input.surface_lipschitz(*graph_group_pos, *raised_esurface_id)
        } else {
            -1.0
        };
    }

    if objective == OverlapCenterObjective::RelaxedChebyshev
        && (1..=esurface_constraints.len()).all(|row| a_matrix[row][0] == 0.0)
    {
        q_vector[0] = 0.0;
    }

    // propagator constraints
    let mut vertical_offset = esurface_constraints.len() + 1;
    for (cone_index, propagator_constraint) in propagator_constraints.iter().enumerate() {
        a_matrix[vertical_offset][propagator_index_offset + cone_index] = -1.0;
        vertical_offset += 1;

        let spatial_shift =
            compute_shift_part(&propagator_constraint.signature.external, external_momenta);

        b_vector[vertical_offset] = spatial_shift.spatial.px.0;
        b_vector[vertical_offset + 1] = spatial_shift.spatial.py.0;
        b_vector[vertical_offset + 2] = spatial_shift.spatial.pz.0;

        for (loop_index, individual_loop_signature) in
            propagator_constraint.signature.internal.iter().enumerate()
        {
            if individual_loop_signature.is_sign() {
                a_matrix[vertical_offset][loop_momentum_offset + 3 * loop_index] =
                    -(*individual_loop_signature as i8) as f64;
                a_matrix[vertical_offset + 1][loop_momentum_offset + 3 * loop_index + 1] =
                    -(*individual_loop_signature as i8) as f64;
                a_matrix[vertical_offset + 2][loop_momentum_offset + 3 * loop_index + 2] =
                    -(*individual_loop_signature as i8) as f64;
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
        // Optional center refinement has a fixed iteration budget, not a scheduling-dependent
        // wall-clock limit. The baseline feasibility solve keeps its original settings.
        .max_iter(if objective == OverlapCenterObjective::MaxMinDepth {
            DefaultSettings::<f64>::default().max_iter
        } else {
            64
        })
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
    .map_err(Into::into)
}

pub(crate) fn find_center(
    overlap_input: &OverlapInput,
    esurfaces_to_consider: &[ExistingEsurfaceId],
    existing_esurfaces: &ExistingEsurfaces,
    external_momenta: &ExternalFourMomenta<F<f64>>,
    verbose: bool,
) -> Result<Option<LoopMomenta<F<f64>>>> {
    let mut solver = construct_solver(
        overlap_input,
        esurfaces_to_consider,
        existing_esurfaces,
        external_momenta,
        verbose,
        OverlapCenterObjective::MaxMinDepth,
        0.0,
    )?;

    solver.solve();

    crate::debug_tags!(#subtraction, #threshold, #overlap, #socp;
        stage = "threshold_socp_result",
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

    let loop_number = overlap_input
        .graph_data
        .first()
        .expect("no graphs passed to overlap")
        .lmb
        .loop_edges
        .len();

    let group_esurfaces_to_check = esurfaces_to_consider
        .iter()
        .map(|&existing_esurface_id| existing_esurfaces[existing_esurface_id])
        .collect_vec();

    // Even if the solver did not converge, check whether its candidate is still valid.
    let center = extract_center(loop_number, &solver.solution.x);
    if overlap_input
        .center_clearance(&group_esurfaces_to_check, &center, external_momenta)
        .is_some()
    {
        return Ok(Some(center));
    }

    if solver.solution.status == SolverStatus::PrimalInfeasible {
        return Ok(None);
    }

    if verbose {
        println!("{:?}", solver.solution.x);
    }

    Err(eyre!(
        "Threshold overlap solve has no valid interior center or primal-infeasibility certificate: status={:?}, existing_surface_ids={:?}, group_surface_ids={:?}, primal_objective={}, dual_objective={}, primal_residual={}, dual_residual={}, iterations={}",
        solver.solution.status,
        esurfaces_to_consider,
        group_esurfaces_to_check,
        solver.solution.obj_val,
        solver.solution.obj_val_dual,
        solver.solution.r_prim,
        solver.solution.r_dual,
        solver.solution.iterations,
    ))
}

pub struct SingleGraphOverlapData<'a> {
    pub lmb: &'a LoopMomentumBasis,
    pub esurfaces: &'a EsurfaceCollection,
    pub raised_data: &'a RaisedEsurfaceData,
    pub edge_masses: EdgeVec<F<f64>>,
}

pub struct OverlapInput<'a> {
    pub graph_data: TiVec<GraphGroupPosition, SingleGraphOverlapData<'a>>,
    pub settings: &'a RuntimeSettings,
    pub group_esurface_map:
        TiVec<GroupEsurfaceId, TiVec<GraphGroupPosition, Option<RaisedEsurfaceId>>>,
    pub local_esurface_exists: TiVec<GraphGroupPosition, TiVec<GroupEsurfaceId, bool>>,
}

fn representative_local_esurface_id(
    graph_data: &SingleGraphOverlapData,
    raised_esurface_id: RaisedEsurfaceId,
) -> EsurfaceID {
    graph_data.raised_data.raised_groups[raised_esurface_id].esurface_ids[0]
}

impl OverlapInput<'_> {
    fn refine_centers(
        &self,
        overlap: &mut OverlapStructure,
        external_momenta: &ExternalFourMomenta<F<f64>>,
    ) {
        let objective = self.settings.subtraction.overlap_settings.objective;
        if objective == OverlapCenterObjective::MaxMinDepth {
            return;
        }
        let loop_count = self
            .graph_data
            .first()
            .expect("overlap geometry")
            .lmb
            .loop_edges
            .len();
        for group in &mut overlap.overlap_groups {
            let selected = group.existing_esurfaces.clone();
            let surfaces = selected
                .iter()
                .map(|id| overlap.existing_esurfaces[*id])
                .collect_vec();
            objective.refine_center(
                &mut group.center,
                self.settings.kinematics.e_cm,
                |objective, radius| {
                    construct_solver(
                        self,
                        &selected,
                        &overlap.existing_esurfaces,
                        external_momenta,
                        false,
                        objective,
                        radius,
                    )
                },
                |coordinates| extract_center(loop_count, coordinates),
                |center| self.center_clearance(&surfaces, center, external_momenta),
            );
        }
    }

    /// Global Lipschitz bound in the Euclidean norm of the active LMB coordinates.
    fn surface_lipschitz(&self, graph_group: GraphGroupPosition, raised: RaisedEsurfaceId) -> f64 {
        let data = &self.graph_data[graph_group];
        data.esurfaces[representative_local_esurface_id(data, raised)]
            .energies
            .iter()
            .map(|edge| {
                data.lmb.edge_signatures[*edge]
                    .internal
                    .iter()
                    .map(|sign| f64::from(*sign as i8).powi(2))
                    .sum::<f64>()
                    .sqrt()
            })
            .sum()
    }

    /// Return the conservative common-ball radius only for a certified interior center.
    fn center_clearance(
        &self,
        group_esurfaces: &[GroupEsurfaceId],
        center: &LoopMomenta<F<f64>>,
        external_momenta: &ExternalFourMomenta<F<f64>>,
    ) -> Option<f64> {
        let mut depth = f64::INFINITY;
        for &group_esurface_id in group_esurfaces {
            let mut has_local_esurface = false;
            for (graph_group_pos, raised_esurface_id) in self.group_esurface_map[group_esurface_id]
                .iter_enumerated()
                .filter_map(|(graph_group_pos, raised)| {
                    raised.and_then(|id| {
                        self.local_esurface_exists[graph_group_pos][group_esurface_id]
                            .then_some((graph_group_pos, id))
                    })
                })
            {
                has_local_esurface = true;
                let data = &self.graph_data[graph_group_pos];
                let surface =
                    &data.esurfaces[representative_local_esurface_id(data, raised_esurface_id)];
                let value = surface.compute_from_momenta(
                    data.lmb,
                    &data.edge_masses,
                    center,
                    external_momenta,
                );
                if !esurface_value_is_strictly_inside(&value, &F(self.settings.kinematics.e_cm)) {
                    return None;
                }
                let lipschitz = self.surface_lipschitz(graph_group_pos, raised_esurface_id);
                if lipschitz > 0.0 {
                    depth = depth.min(-value.0 / lipschitz);
                }
            }
            if !has_local_esurface {
                return None;
            }
        }
        Some(depth)
    }
}

pub(crate) fn check_global_center(
    overlap_input: &OverlapInput,
    existing_esurfaces: &ExistingEsurfaces,
    center: &LoopMomenta<F<f64>>,
    external_momenta: &ExternalFourMomenta<F<f64>>,
) -> bool {
    let group_esurfaces = existing_esurfaces.iter().copied().collect_vec();
    overlap_input
        .center_clearance(&group_esurfaces, center, external_momenta)
        .is_some()
}

/// Runtime overlap failures are returned so the stability machinery can retry at higher precision.
/// Structural generation invariants are still asserted where malformed generated data is unrecoverable.
pub(crate) fn find_maximal_overlap(
    overlap_input: &OverlapInput,
    existing_esurfaces: &ExistingEsurfaces,
    external_momenta: &ExternalFourMomenta<F<f64>>,
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
        let global_center_f = global_center
            .iter()
            .map(|coordinates| ThreeMomentum {
                px: F(coordinates[0]),
                py: F(coordinates[1]),
                pz: F(coordinates[2]),
            })
            .collect();

        if !settings.subtraction.overlap_settings.check_global_center {
            tracing::warn!(
                "overlap_settings.check_global_center=false is deprecated; forced centers are always validated"
            );
        }

        let is_valid = check_global_center(
            overlap_input,
            existing_esurfaces,
            &global_center_f,
            external_momenta,
        );

        if !is_valid {
            return Err(eyre!(
                "Center provided is not finite and strictly inside all existing esurfaces"
            ));
        }

        let single_group = OverlapGroup {
            existing_esurfaces: all_existing_esurfaces,
            center: global_center_f,
            complement: vec![],
            prefactor_evaluator: None,
        };
        res.overlap_groups.push(single_group);

        res.fill_in_complements();
        return Ok(res);
    }

    if settings.subtraction.overlap_settings.enable_heuristics
        && settings.subtraction.overlap_settings.try_origin
    {
        let global_loop_count = overlap_input
            .graph_data
            .first()
            .unwrap()
            .lmb
            .loop_edges
            .len();
        let origin = LoopMomenta::from_iter(
            (0..global_loop_count).map(|_| ThreeMomentum::new(F(0.0), F(0.0), F(0.0))),
        );

        let is_valid =
            check_global_center(overlap_input, existing_esurfaces, &origin, external_momenta);

        if is_valid {
            let single_group = OverlapGroup {
                existing_esurfaces: all_existing_esurfaces,
                center: origin,
                complement: vec![],
                prefactor_evaluator: None,
            };
            res.overlap_groups.push(single_group);
            res.fill_in_complements();
            return Ok(res);
        }
    }

    if settings.subtraction.overlap_settings.enable_heuristics
        && settings.subtraction.overlap_settings.try_origin_all_lmbs
    {
        todo!("Not all heuristics implemented")
    }

    // first try if all esurfaces have a single center, we explitely seach a center instead of trying the
    // origin. This is because the origin might not be optimal.
    let option_center = find_center(
        overlap_input,
        &all_existing_esurfaces,
        existing_esurfaces,
        external_momenta,
        false,
    )?;

    if let Some(center) = option_center {
        let single_group = OverlapGroup {
            existing_esurfaces: all_existing_esurfaces,
            center,
            complement: vec![],
            prefactor_evaluator: None,
        };
        res.overlap_groups.push(single_group);
        res.fill_in_complements();
        overlap_input.refine_centers(&mut res, external_momenta);
        return Ok(res);
    }

    // If the full intersection is infeasible, create a table of all pairs.
    let esurface_pairs = EsurfacePairs::new(overlap_input, existing_esurfaces, external_momenta)?;

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
                external_momenta,
                false,
            )?
            .ok_or_else(|| {
                let (graph_group_pos, raised_esurface_id) =
                    get_representative(&overlap_input.group_esurface_map[esurface_id])
                        .expect("overlap corrupted");
                let esurface_id = representative_local_esurface_id(
                    &overlap_input.graph_data[graph_group_pos],
                    raised_esurface_id,
                );

                let esurface = &overlap_input.graph_data[graph_group_pos].esurfaces[esurface_id];

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
                prefactor_evaluator: None,
            });
            num_disconnected_surfaces += 1;
        }
    }

    // if settings.general.debug > 3 {
    //     DEBUG_LOGGER.write("num_disconnected_surfaces", &num_disconnected_surfaces);
    // }

    if num_disconnected_surfaces == existing_esurfaces.len() {
        res.fill_in_complements();
        overlap_input.refine_centers(&mut res, external_momenta);
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

        for subset in possible_subsets.iter().sorted() {
            let option_center = find_center(
                overlap_input,
                subset,
                existing_esurfaces,
                external_momenta,
                false,
            )?;

            if let Some(center) = option_center {
                res.overlap_groups.push(OverlapGroup {
                    existing_esurfaces: subset.clone(),
                    center,
                    complement: vec![],
                    prefactor_evaluator: None,
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
    overlap_input.refine_centers(&mut res, external_momenta);
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
        existing_esurfaces: &ExistingEsurfaces,
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
        existing_esurfaces: &ExistingEsurfaces,
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
#[allow(dead_code, unused_variables)]
mod tests {
    use super::*;
    use itertools::Itertools;
    use linnet::half_edge::{
        involution::EdgeIndex,
        subgraph::{SuBitGraph, SubSetLike},
    };
    use typed_index_collections::ti_vec;

    use crate::{
        cff::{
            VertexSet,
            esurface::{
                Esurface, EsurfaceExistence, EsurfaceID, RaisedEsurfaceData, RaisedEsurfaceGroup,
                RaisedEsurfaceId,
            },
        },
        graph::LoopMomentumBasis,
        momentum::FourMomentum,
        momentum::signature::LoopExtSignature,
        settings::RuntimeSettings,
        utils::test_utils::dummy_hedge_graph,
    };

    #[test]
    fn overlap_structure_localizes_to_graph_existing_surfaces() {
        let center = LoopMomenta::from_iter([ThreeMomentum::new(F(0.0), F(0.0), F(0.0))]);
        let overlap = OverlapStructure {
            existing_esurfaces: ti_vec![
                GroupEsurfaceId::from(0),
                GroupEsurfaceId::from(1),
                GroupEsurfaceId::from(2),
            ],
            overlap_groups: vec![
                OverlapGroup {
                    existing_esurfaces: vec![
                        ExistingEsurfaceId::from(0),
                        ExistingEsurfaceId::from(1),
                    ],
                    complement: vec![],
                    center: center.clone(),
                    prefactor_evaluator: None,
                },
                OverlapGroup {
                    existing_esurfaces: vec![
                        ExistingEsurfaceId::from(1),
                        ExistingEsurfaceId::from(2),
                    ],
                    complement: vec![],
                    center: center.clone(),
                    prefactor_evaluator: None,
                },
            ],
        };

        let localized = overlap.localized_to_existing_surfaces(&ti_vec![true, false, true]);

        assert_eq!(
            localized
                .existing_esurfaces
                .iter()
                .map(|group_esurface_id| group_esurface_id.0)
                .collect_vec(),
            vec![0, 2]
        );
        assert_eq!(localized.overlap_groups.len(), 2);
        assert_eq!(
            localized.overlap_groups[0].existing_esurfaces,
            vec![ExistingEsurfaceId::from(0)]
        );
        assert_eq!(
            localized.overlap_groups[0].complement,
            vec![ExistingEsurfaceId::from(1)]
        );
        assert_eq!(
            localized.overlap_groups[1].existing_esurfaces,
            vec![ExistingEsurfaceId::from(1)]
        );
        assert_eq!(
            localized.overlap_groups[1].complement,
            vec![ExistingEsurfaceId::from(0)]
        );
    }

    struct HelperBoxStructure {
        external_momenta: ExternalFourMomenta<F<f64>>,
        lmb: LoopMomentumBasis,
        esurfaces: EsurfaceCollection,
        raised_data: RaisedEsurfaceData,
        existing_esurfaces: ExistingEsurfaces,
        edge_masses: EdgeVec<F<f64>>,
    }

    struct HelperBananaStructure {
        external_momenta: ExternalFourMomenta<F<f64>>,
        lmb: LoopMomentumBasis,
        esurfaces: EsurfaceCollection,
        raised_data: RaisedEsurfaceData,
        existing_esurfaces: ExistingEsurfaces,
        edge_masses: EdgeVec<F<f64>>,
    }

    fn trivial_raised_data(num_esurfaces: usize) -> RaisedEsurfaceData {
        RaisedEsurfaceData {
            raised_groups: (0..num_esurfaces)
                .map(|index| RaisedEsurfaceGroup {
                    esurface_ids: vec![EsurfaceID::from(index)],
                    max_occurence: 1,
                })
                .collect(),
            pass_two_evaluator: None,
        }
    }

    impl HelperBoxStructure {
        fn new(masses: Option<[F<f64>; 4]>) -> Self {
            let external_momenta = ExternalFourMomenta::from_iter([
                FourMomentum::from_args(F(14.0), F(-6.6), F(-40.0), F(0.0)),
                FourMomentum::from_args(F(-43.0), F(15.2), F(33.0), F(0.0)),
                FourMomentum::from_args(F(-17.9), F(-50.0), F(11.8), F(0.0)),
            ]);

            let dummy_hedge_graph = dummy_hedge_graph(8);

            let box_basis = ti_vec![EdgeIndex::from(4)];
            let box_signatures = dummy_hedge_graph
                .new_edgevec_from_iter(vec![
                    (vec![0], vec![1, 0, 0]).into(),
                    (vec![0], vec![0, 1, 0]).into(),
                    (vec![0], vec![0, 0, 1]).into(),
                    (vec![0], vec![-1, -1, -1]).into(),
                    (vec![1], vec![0, 0, 0]).into(),
                    (vec![1], vec![1, 0, 0]).into(),
                    (vec![1], vec![1, 1, 0]).into(),
                    (vec![1], vec![1, 1, 1]).into(),
                ])
                .unwrap();

            let box_lmb = LoopMomentumBasis {
                tree: SuBitGraph::empty(0),
                ext_edges: vec![].into(),
                loop_edges: box_basis,
                edge_signatures: box_signatures,
            };

            let esurfaces_array = [
                Esurface {
                    energies: vec![EdgeIndex::from(5), EdgeIndex::from(6)],
                    external_shift: vec![(EdgeIndex::from(1), 1)],
                    vertex_set: VertexSet::dummy(),
                    // subspace_graph: dummy_hedge_graph.full_graph(),
                },
                Esurface {
                    energies: vec![EdgeIndex::from(5), EdgeIndex::from(7)],
                    external_shift: vec![(EdgeIndex::from(1), 1), (EdgeIndex::from(2), 1)],
                    vertex_set: VertexSet::dummy(),
                    //subspace_graph: dummy_hedge_graph.full_graph(),
                },
                Esurface {
                    energies: vec![EdgeIndex::from(4), EdgeIndex::from(6)],
                    external_shift: vec![(EdgeIndex::from(0), 1), (EdgeIndex::from(1), 1)],
                    vertex_set: VertexSet::dummy(),
                    //subspace_graph: dummy_hedge_graph.full_graph(),
                },
                Esurface {
                    energies: vec![EdgeIndex::from(4), EdgeIndex::from(7)],
                    external_shift: vec![
                        (EdgeIndex::from(0), 1),
                        (EdgeIndex::from(1), 1),
                        (EdgeIndex::from(2), 1),
                    ],
                    vertex_set: VertexSet::dummy(),
                    //subspace_graph: dummy_hedge_graph.full_graph(),
                },
            ];

            let esurfaces = esurfaces_array.to_vec().into();

            let edge_masses = match masses {
                Some(masses) => {
                    let mut edge_masses = vec![F(0.0); 4];
                    let mut real_masses = masses.iter().copied().collect_vec();

                    edge_masses.append(&mut real_masses);
                    edge_masses
                }
                None => vec![F(0.0); 8],
            };

            let existing_esurfaces = (0..4).map(Into::<GroupEsurfaceId>::into).collect();

            Self {
                external_momenta,
                lmb: box_lmb,
                existing_esurfaces,
                esurfaces,
                raised_data: trivial_raised_data(4),
                edge_masses: dummy_hedge_graph
                    .new_edgevec_from_iter(edge_masses)
                    .unwrap(),
            }
        }
    }

    impl HelperBananaStructure {
        fn new() -> Self {
            let external_momenta = ExternalFourMomenta::from_iter([FourMomentum::from_args(
                F(10.0),
                F(-10.00000000),
                F(0.0),
                F(0.0),
            )]);
            let banana_basis = ti_vec![EdgeIndex::from(2), EdgeIndex::from(3)];

            let dummy_hedge_graph = dummy_hedge_graph(5);

            let banana_edge_sigs = dummy_hedge_graph
                .new_edgevec_from_iter(vec![
                    LoopExtSignature {
                        internal: vec![0, 0].into(),
                        external: vec![1].into(),
                    },
                    LoopExtSignature {
                        internal: vec![0, 0].into(),
                        external: vec![-1].into(),
                    },
                    LoopExtSignature {
                        internal: vec![1, 0].into(),
                        external: vec![0].into(),
                    },
                    LoopExtSignature {
                        internal: vec![0, 1].into(),
                        external: vec![0].into(),
                    },
                    LoopExtSignature {
                        internal: vec![1, 1].into(),
                        external: vec![-1].into(),
                    },
                ])
                .unwrap();

            let banana_lmb = LoopMomentumBasis {
                tree: SuBitGraph::empty(0),
                loop_edges: banana_basis,
                ext_edges: vec![].into(),
                edge_signatures: banana_edge_sigs,
            };

            let only_esurface = Esurface {
                energies: vec![EdgeIndex::from(2), EdgeIndex::from(3), EdgeIndex::from(4)],
                external_shift: vec![(EdgeIndex::from(0), -1)],
                vertex_set: VertexSet::dummy(),
                //subspace_graph: dummy_hedge_graph.full_graph(),
            };

            let esurfaces = vec![only_esurface].into();

            let existing_esurfaces = vec![Into::<GroupEsurfaceId>::into(0)].into();
            let edge_masses = dummy_hedge_graph
                .new_edgevec_from_iter(vec![F(0.0); 5])
                .unwrap();

            Self {
                external_momenta,
                lmb: banana_lmb,
                esurfaces,
                raised_data: trivial_raised_data(1),
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
                LoopMomenta::from(vec![]),
            ),
            (
                vec![
                    Into::<ExistingEsurfaceId>::into(1),
                    Into::<ExistingEsurfaceId>::into(2),
                    Into::<ExistingEsurfaceId>::into(4),
                ],
                LoopMomenta::from(vec![]),
            ),
            (
                vec![
                    Into::<ExistingEsurfaceId>::into(2),
                    Into::<ExistingEsurfaceId>::into(3),
                    Into::<ExistingEsurfaceId>::into(4),
                ],
                LoopMomenta::from(vec![]),
            ),
        ];

        let fake_res = OverlapStructure {
            overlap_groups: fake_res
                .into_iter()
                .map(|(group, center)| OverlapGroup {
                    existing_esurfaces: group,
                    center,
                    complement: vec![],
                    prefactor_evaluator: None,
                })
                .collect_vec(),
            existing_esurfaces: ti_vec![],
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

        let massless_overlap_input = OverlapInput {
            graph_data: ti_vec![SingleGraphOverlapData {
                lmb: &box4e.lmb,
                esurfaces: &box4e.esurfaces,
                raised_data: &box4e.raised_data,
                edge_masses: box4e.edge_masses.clone(),
            }],
            settings: &RuntimeSettings::default(),
            group_esurface_map: (0..4)
                .map(|i| ti_vec![Some(Into::<RaisedEsurfaceId>::into(i))])
                .collect(),
            local_esurface_exists: ti_vec![ti_vec![true; 4]],
        };

        let esurface_pairs = EsurfacePairs::new(
            &massless_overlap_input,
            &box4e.existing_esurfaces,
            &box4e.external_momenta,
        )
        .unwrap();

        assert_eq!(esurface_pairs.data.len(), 4);

        let box4e_massive = HelperBoxStructure::new(Some([F(10.5); 4]));

        let massive_overlap_input = OverlapInput {
            graph_data: ti_vec![SingleGraphOverlapData {
                lmb: &box4e_massive.lmb,
                esurfaces: &box4e_massive.esurfaces,
                raised_data: &box4e_massive.raised_data,
                edge_masses: box4e_massive.edge_masses.clone(),
            }],
            settings: &RuntimeSettings::default(),
            group_esurface_map: (0..4)
                .map(|i| ti_vec![Some(Into::<RaisedEsurfaceId>::into(i))])
                .collect(),
            local_esurface_exists: ti_vec![ti_vec![true; 4]],
        };

        let esurface_pairs_massive = EsurfacePairs::new(
            &massive_overlap_input,
            &box4e_massive.existing_esurfaces,
            &box4e_massive.external_momenta,
        )
        .unwrap();

        assert_eq!(esurface_pairs_massive.data.len(), 0);
    }

    #[test]
    fn maximal_overlap_candidates_include_cliques_whose_pairs_are_already_covered() {
        let existing: ExistingEsurfaces = (0..9).map(GroupEsurfaceId::from).collect();
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
                    prefactor_evaluator: None,
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
    fn test_subset_generator() {
        let box4e = HelperBoxStructure::new(None);

        let overlap_input = OverlapInput {
            graph_data: ti_vec![SingleGraphOverlapData {
                lmb: &box4e.lmb,
                esurfaces: &box4e.esurfaces,
                raised_data: &box4e.raised_data,
                edge_masses: box4e.edge_masses.clone(),
            }],
            settings: &RuntimeSettings::default(),
            group_esurface_map: (0..4)
                .map(|i| ti_vec![Some(Into::<RaisedEsurfaceId>::into(i))])
                .collect(),
            local_esurface_exists: ti_vec![ti_vec![true; 4]],
        };

        let esurface_pairs = EsurfacePairs::new(
            &overlap_input,
            &box4e.existing_esurfaces,
            &box4e.external_momenta,
        )
        .unwrap();

        let res = OverlapStructure {
            overlap_groups: vec![],
            existing_esurfaces: box4e.existing_esurfaces.clone(),
        };
        let subsets_3 =
            esurface_pairs.construct_possible_subsets_of_len(&box4e.existing_esurfaces, 3, &res);

        assert_eq!(subsets_3.len(), 0);

        let subsets_2 =
            esurface_pairs.construct_possible_subsets_of_len(&box4e.existing_esurfaces, 2, &res);
        assert_eq!(subsets_2.len(), 4);
    }

    #[test]
    fn min_sum_center_counts_shared_physical_surfaces_once_across_graphs() {
        // The first loop fixes the best ball radius, while the signed sum fixes the
        // second-loop center. Duplicating one graph's E2 must not give E2 extra weight.
        let mut fixture = HelperBoxStructure::new(Some([F(0.0), F(1.0), F(1.0), F(0.0)]));
        fixture.lmb.loop_edges = ti_vec![EdgeIndex(4), EdgeIndex(7)];
        for (edge, internal, external) in [
            (0, vec![0, 0], vec![1, 0, 0]),
            (1, vec![0, 0], vec![0, 1, 0]),
            (2, vec![0, 0], vec![0, 0, 1]),
            (3, vec![0, 0], vec![-1, -1, -1]),
            (4, vec![1, 0], vec![0, 0, 0]),
            (5, vec![0, 1], vec![0, 1, 0]),
            (6, vec![0, 1], vec![0, 0, 1]),
            (7, vec![0, 1], vec![0, 0, 0]),
        ] {
            fixture.lmb.edge_signatures[EdgeIndex(edge)] = (internal, external).into();
        }
        fixture.external_momenta = [
            FourMomentum::from_args(F(1.0), F(0.0), F(0.0), F(0.0)),
            FourMomentum::from_args(F(10.0), F(-1.0), F(0.0), F(0.0)),
            FourMomentum::from_args(F(10.0), F(1.0), F(0.0), F(0.0)),
        ]
        .into_iter()
        .collect();
        fixture.esurfaces = (0..3)
            .map(|id| Esurface {
                energies: vec![EdgeIndex(4 + id)],
                external_shift: vec![(EdgeIndex(id), -1)],
                vertex_set: VertexSet::dummy(),
            })
            .collect();
        fixture.raised_data = trivial_raised_data(3);
        fixture.existing_esurfaces = (0..3).map(GroupEsurfaceId).collect();
        let mut settings = RuntimeSettings::default();
        settings.kinematics.e_cm = 10.0;
        settings.subtraction.overlap_settings.enable_heuristics = false;
        settings.subtraction.overlap_settings.objective = OverlapCenterObjective::MinSum;
        let mut previous: Option<f64> = None;
        for copies in [1, 2] {
            let input = OverlapInput {
                graph_data: (0..copies)
                    .map(|_| SingleGraphOverlapData {
                        lmb: &fixture.lmb,
                        esurfaces: &fixture.esurfaces,
                        raised_data: &fixture.raised_data,
                        edge_masses: fixture.edge_masses.clone(),
                    })
                    .collect(),
                settings: &settings,
                group_esurface_map: (0..3)
                    .map(|id| {
                        (0..copies)
                            .map(|copy| (copy == 0 || id == 1).then_some(RaisedEsurfaceId(id)))
                            .collect()
                    })
                    .collect(),
                local_esurface_exists: (0..copies)
                    .map(|copy| (0..3).map(|id| copy == 0 || id == 1).collect())
                    .collect(),
            };
            let result = find_maximal_overlap(
                &input,
                &fixture.existing_esurfaces,
                &fixture.external_momenta,
            )
            .unwrap();
            assert_eq!(result.overlap_groups.len(), 1);
            let center = &result.overlap_groups[0].center;
            let x = center[crate::momentum::sample::LoopIndex(1)].px.0;
            assert!(
                x.abs() < 5.0e-4,
                "duplicated graph must not weight one physical surface: {x}"
            );
            if let Some(previous) = previous {
                assert!((x - previous).abs() < 5.0e-4);
            }
            previous = Some(x);
            assert!(check_global_center(
                &input,
                &fixture.existing_esurfaces,
                center,
                &fixture.external_momenta
            ));
        }
    }

    #[test]
    fn overlap_refinement_keeps_certified_witness_on_every_optional_failure() {
        use std::cell::Cell;
        let original: LoopMomenta<_> = [ThreeMomentum::new(F(2.0), F(0.0), F(0.0))]
            .into_iter()
            .collect();
        let extract = |x: &[f64]| {
            [ThreeMomentum::new(F(x[0]), F(0.0), F(0.0))]
                .into_iter()
                .collect()
        };
        let certify = |center: &LoopMomenta<F<f64>>| {
            let x = center[crate::momentum::sample::LoopIndex(0)].px.0;
            (x.is_finite() && x < 3.0).then_some(3.0 - x)
        };
        for failure in 0..5 {
            let mut center = original.clone();
            let calls = Cell::new(0);
            OverlapCenterObjective::MinSum.refine_center(
                &mut center,
                10.0,
                |objective, _| {
                    calls.set(calls.get() + 1);
                    if failure == 0 || (failure == 2 && objective == OverlapCenterObjective::MinSum)
                    {
                        return Err(eyre!("deliberate optional-construction failure"));
                    }
                    DefaultSolver::new(
                        &CscMatrix::spalloc((1, 1), 0),
                        &[1.0],
                        &CscMatrix::from(&[[1.0], [-1.0]]),
                        &[2.0, -1.0],
                        &[NonnegativeConeT(2)],
                        DefaultSettingsBuilder::default()
                            .verbose(false)
                            .max_iter(if failure == 1 { 0 } else { 64 })
                            .build()
                            .unwrap(),
                    )
                    .map_err(Into::into)
                },
                |x| {
                    if failure == 4 && calls.get() == 2 {
                        // A physically interior phase-II candidate can still violate the
                        // certified phase-I radius floor and must be rejected.
                        extract(&[2.5])
                    } else {
                        extract(x)
                    }
                },
                |candidate| {
                    if failure == 3 && candidate != &original {
                        None
                    } else {
                        certify(candidate)
                    }
                },
            );
            if failure < 2 || failure == 3 {
                assert_eq!(
                    center, original,
                    "failed optimization must preserve the existing witness"
                );
                assert_eq!(calls.get(), 1);
            } else {
                assert!(
                    (center[crate::momentum::sample::LoopIndex(0)].px.0 - 1.0).abs() < 1.0e-7,
                    "phase-II failure keeps phase-I's certified improvement"
                );
                assert_eq!(calls.get(), 2);
            }
        }
        let mut center = original.clone();
        OverlapCenterObjective::MaxMinDepth.refine_center(
            &mut center,
            10.0,
            |_, _| panic!("default must perform no optional solve"),
            extract,
            certify,
        );
        assert_eq!(center, original);
    }

    #[test]
    fn overlap_refinement_preserves_the_maximal_catalogue() {
        let fixture = HelperBoxStructure::new(None);
        let mut previous = None;
        for objective in [
            OverlapCenterObjective::MaxMinDepth,
            OverlapCenterObjective::RelaxedChebyshev,
            OverlapCenterObjective::MinSum,
        ] {
            let mut settings = RuntimeSettings::default();
            settings.subtraction.overlap_settings.enable_heuristics = false;
            settings.subtraction.overlap_settings.objective = objective;
            let input = OverlapInput {
                graph_data: ti_vec![SingleGraphOverlapData {
                    lmb: &fixture.lmb,
                    esurfaces: &fixture.esurfaces,
                    raised_data: &fixture.raised_data,
                    edge_masses: fixture.edge_masses.clone()
                }],
                settings: &settings,
                group_esurface_map: (0..4).map(|i| ti_vec![Some(RaisedEsurfaceId(i))]).collect(),
                local_esurface_exists: ti_vec![ti_vec![true; 4]],
            };
            for _ in 0..3 {
                let result = find_maximal_overlap(
                    &input,
                    &fixture.existing_esurfaces,
                    &fixture.external_momenta,
                )
                .unwrap();
                let catalogue = result
                    .overlap_groups
                    .iter()
                    .map(|group| (group.existing_esurfaces.clone(), group.complement.clone()))
                    .collect_vec();
                if let Some(expected) = &previous {
                    assert_eq!(&catalogue, expected);
                }
                previous = Some(catalogue);
                for group in &result.overlap_groups {
                    let surfaces = group
                        .existing_esurfaces
                        .iter()
                        .map(|id| result.existing_esurfaces[*id])
                        .collect_vec();
                    assert!(
                        input
                            .center_clearance(&surfaces, &group.center, &fixture.external_momenta)
                            .is_some()
                    );
                }
            }
        }
    }

    #[test]
    fn overlap_objectives_solve_their_analytic_lens_and_preserve_energy_scaling() {
        // E1=|k|-2, E2=2|k-3 ex|-4. The unequal routing weights distinguish
        // energy max-min (x=4/3) from the relaxed Euclidean center (x=3/2).
        let mut fixture = HelperBoxStructure::new(None);
        fixture.lmb.edge_signatures[EdgeIndex(5)] = (vec![1], vec![0, 1, 0]).into();
        fixture.lmb.edge_signatures[EdgeIndex(6)] =
            fixture.lmb.edge_signatures[EdgeIndex(5)].clone();
        fixture.esurfaces = vec![
            Esurface {
                energies: vec![EdgeIndex(4)],
                external_shift: vec![(EdgeIndex(0), -1)],
                vertex_set: VertexSet::dummy(),
            },
            Esurface {
                energies: vec![EdgeIndex(5), EdgeIndex(6)],
                external_shift: vec![(EdgeIndex(1), -1)],
                vertex_set: VertexSet::dummy(),
            },
        ]
        .into();
        fixture.raised_data = trivial_raised_data(2);
        fixture.existing_esurfaces = ti_vec![GroupEsurfaceId(0), GroupEsurfaceId(1)];
        for scale in [1.0, 8.0, 0.125] {
            fixture.external_momenta = [
                FourMomentum::from_args(F(2.0 * scale), F(0.0), F(0.0), F(0.0)),
                FourMomentum::from_args(F(4.0 * scale), F(-3.0 * scale), F(0.0), F(0.0)),
                FourMomentum::from_args(F(0.0), F(0.0), F(0.0), F(0.0)),
            ]
            .into_iter()
            .collect();
            for (objective, expected) in [
                (OverlapCenterObjective::MaxMinDepth, 4.0 / 3.0),
                (OverlapCenterObjective::RelaxedChebyshev, 1.5),
                (OverlapCenterObjective::MinSum, 1.5),
            ] {
                let mut settings = RuntimeSettings::default();
                settings.kinematics.e_cm = 10.0 * scale;
                settings.subtraction.overlap_settings.objective = objective;
                settings.subtraction.overlap_settings.enable_heuristics = false;
                settings.subtraction.overlap_settings.try_origin_all_lmbs = true;
                let input = OverlapInput {
                    graph_data: ti_vec![SingleGraphOverlapData {
                        lmb: &fixture.lmb,
                        esurfaces: &fixture.esurfaces,
                        raised_data: &fixture.raised_data,
                        edge_masses: fixture.edge_masses.clone(),
                    }],
                    settings: &settings,
                    group_esurface_map: ti_vec![
                        ti_vec![Some(RaisedEsurfaceId(0))],
                        ti_vec![Some(RaisedEsurfaceId(1))]
                    ],
                    local_esurface_exists: ti_vec![ti_vec![true; 2]],
                };
                let result = find_maximal_overlap(
                    &input,
                    &fixture.existing_esurfaces,
                    &fixture.external_momenta,
                )
                .unwrap();
                assert_eq!(result.overlap_groups.len(), 1);
                let center = &result.overlap_groups[0].center;
                assert!(
                    (center[crate::momentum::sample::LoopIndex(0)].px.0 / scale - expected).abs()
                        < 2.0e-6,
                    "{objective:?}: {center}"
                );
                assert!(check_global_center(
                    &input,
                    &fixture.existing_esurfaces,
                    center,
                    &fixture.external_momenta
                ));
                let reversed: ExistingEsurfaces =
                    fixture.existing_esurfaces.iter().rev().copied().collect();
                let reordered =
                    find_maximal_overlap(&input, &reversed, &fixture.external_momenta).unwrap();
                assert!(
                    (reordered.overlap_groups[0].center[crate::momentum::sample::LoopIndex(0)]
                        .px
                        .0
                        / scale
                        - expected)
                        .abs()
                        < 2.0e-6
                );
            }
        }
        // A valid heuristic keeps its exact old result for every objective. Disabling
        // the master gate forces a solve without changing the individual heuristic options.
        fixture.external_momenta[crate::momentum::sample::ExternalIndex(0)]
            .temporal
            .value = F(5.0);
        fixture.external_momenta[crate::momentum::sample::ExternalIndex(1)] =
            FourMomentum::from_args(F(8.0), F(-3.0), F(0.0), F(0.0));
        for objective in [
            OverlapCenterObjective::MaxMinDepth,
            OverlapCenterObjective::RelaxedChebyshev,
            OverlapCenterObjective::MinSum,
        ] {
            let mut settings = RuntimeSettings::default();
            settings.subtraction.overlap_settings.objective = objective;
            for (enabled, forced) in [(true, false), (false, false), (false, true)] {
                settings.subtraction.overlap_settings.enable_heuristics = enabled;
                settings.subtraction.overlap_settings.force_global_center =
                    forced.then_some(vec![[1.0, 0.0, 0.0]]);
                let input = OverlapInput {
                    graph_data: ti_vec![SingleGraphOverlapData {
                        lmb: &fixture.lmb,
                        esurfaces: &fixture.esurfaces,
                        raised_data: &fixture.raised_data,
                        edge_masses: fixture.edge_masses.clone()
                    }],
                    settings: &settings,
                    group_esurface_map: ti_vec![
                        ti_vec![Some(RaisedEsurfaceId(0))],
                        ti_vec![Some(RaisedEsurfaceId(1))]
                    ],
                    local_esurface_exists: ti_vec![ti_vec![true; 2]],
                };
                let result = find_maximal_overlap(
                    &input,
                    &fixture.existing_esurfaces,
                    &fixture.external_momenta,
                )
                .unwrap();
                let x = result.overlap_groups[0].center[crate::momentum::sample::LoopIndex(0)]
                    .px
                    .0;
                if enabled {
                    assert_eq!(x, 0.0);
                } else if forced {
                    assert_eq!(x, 1.0);
                } else {
                    assert!(x > 0.1);
                }
            }
        }
    }

    #[test]
    fn overlap_epigraphs_preserve_repeated_energies_and_distinct_masses() {
        let mut fixture = HelperBoxStructure::new(Some([F(1.0), F(2.0), F(1.0), F(0.0)]));
        fixture.lmb.edge_signatures[EdgeIndex(5)] =
            fixture.lmb.edge_signatures[EdgeIndex(4)].clone();
        fixture.lmb.edge_signatures[EdgeIndex(6)] = (vec![1], vec![0, 0, 1]).into();
        fixture.external_momenta = [
            FourMomentum::from_args(F(3.0), F(0.0), F(0.0), F(0.0)),
            FourMomentum::from_args(F(2.1), F(0.0), F(0.0), F(0.0)),
            FourMomentum::from_args(F(2.0), F(-2.0), F(0.0), F(0.0)),
        ]
        .into_iter()
        .collect();
        fixture.esurfaces = (0..3)
            .map(|i| Esurface {
                energies: vec![EdgeIndex(4 + i)],
                external_shift: vec![(EdgeIndex(i), -1)],
                vertex_set: VertexSet::dummy(),
            })
            .collect();
        fixture.raised_data = trivial_raised_data(3);
        fixture.existing_esurfaces = (0..3).map(GroupEsurfaceId).collect();
        let mut settings = RuntimeSettings::default();
        settings.subtraction.overlap_settings.enable_heuristics = false;
        for objective in [
            OverlapCenterObjective::MaxMinDepth,
            OverlapCenterObjective::RelaxedChebyshev,
            OverlapCenterObjective::MinSum,
        ] {
            settings.subtraction.overlap_settings.objective = objective;
            let input = OverlapInput {
                graph_data: ti_vec![SingleGraphOverlapData {
                    lmb: &fixture.lmb,
                    esurfaces: &fixture.esurfaces,
                    raised_data: &fixture.raised_data,
                    edge_masses: fixture.edge_masses.clone()
                }],
                settings: &settings,
                group_esurface_map: (0..3).map(|i| ti_vec![Some(RaisedEsurfaceId(i))]).collect(),
                local_esurface_exists: ti_vec![ti_vec![true; 3]],
            };
            let result = find_maximal_overlap(
                &input,
                &fixture.existing_esurfaces,
                &fixture.external_momenta,
            )
            .unwrap();
            assert_eq!(result.overlap_groups.len(), 1);
            let center = result.overlap_groups[0].center.clone();
            assert!(check_global_center(
                &input,
                &fixture.existing_esurfaces,
                &center,
                &fixture.external_momenta
            ));
            assert!(
                center[crate::momentum::sample::LoopIndex(0)].px.0 < 0.65,
                "mass-1 epigraph must not replace mass-2 energy"
            );
        }
        fixture.edge_masses[EdgeIndex(5)] = F(1.0);
        fixture.esurfaces[EsurfaceID(0)].energies = vec![EdgeIndex(4), EdgeIndex(5)];
        fixture.external_momenta[crate::momentum::sample::ExternalIndex(0)]
            .temporal
            .value = F(1.5);
        let input = OverlapInput {
            graph_data: ti_vec![SingleGraphOverlapData {
                lmb: &fixture.lmb,
                esurfaces: &fixture.esurfaces,
                raised_data: &fixture.raised_data,
                edge_masses: fixture.edge_masses.clone()
            }],
            settings: &settings,
            group_esurface_map: (0..3).map(|i| ti_vec![Some(RaisedEsurfaceId(i))]).collect(),
            local_esurface_exists: ti_vec![ti_vec![true; 3]],
        };
        assert!(
            find_center(
                &input,
                &[ExistingEsurfaceId::from(0)],
                &fixture.existing_esurfaces,
                &fixture.external_momenta,
                false
            )
            .unwrap()
            .is_none(),
            "2sqrt(k²+1)-1.5 is everywhere positive"
        );
    }

    #[test]
    fn overlap_objectives_accept_stationary_massless_nonunique_centers() {
        let mut fixture = HelperBoxStructure::new(None);
        fixture.lmb.edge_signatures[EdgeIndex(5)] = (vec![1], vec![0, 1, 0]).into();
        fixture.lmb.edge_signatures[EdgeIndex(6)] = (vec![1], vec![0, 0, 1]).into();
        fixture.external_momenta = [
            FourMomentum::from_args(F(3.0), F(0.0), F(0.0), F(0.0)),
            FourMomentum::from_args(F(0.0), F(-1.0), F(0.0), F(0.0)),
            FourMomentum::from_args(F(0.0), F(1.0), F(0.0), F(0.0)),
        ]
        .into_iter()
        .collect();
        fixture.esurfaces[EsurfaceID(0)] = Esurface {
            energies: vec![EdgeIndex(5), EdgeIndex(6)],
            external_shift: vec![(EdgeIndex(0), -1)],
            vertex_set: VertexSet::dummy(),
        };
        for objective in [
            OverlapCenterObjective::MaxMinDepth,
            OverlapCenterObjective::RelaxedChebyshev,
            OverlapCenterObjective::MinSum,
        ] {
            let mut settings = RuntimeSettings::default();
            settings.subtraction.overlap_settings.objective = objective;
            settings.subtraction.overlap_settings.enable_heuristics = false;
            let input = OverlapInput {
                graph_data: ti_vec![SingleGraphOverlapData {
                    lmb: &fixture.lmb,
                    esurfaces: &fixture.esurfaces,
                    raised_data: &fixture.raised_data,
                    edge_masses: fixture.edge_masses.clone()
                }],
                settings: &settings,
                group_esurface_map: (0..4).map(|i| ti_vec![Some(RaisedEsurfaceId(i))]).collect(),
                local_esurface_exists: ti_vec![ti_vec![true; 4]],
            };
            let result = find_maximal_overlap(
                &input,
                &ti_vec![GroupEsurfaceId(0)],
                &fixture.external_momenta,
            )
            .unwrap();
            let center = result.overlap_groups[0].center.clone();
            // Every point on the focal segment is optimal: assert the physical contract,
            // not an invented uniqueness guarantee or one solver-specific bit pattern.
            assert!(
                input
                    .center_clearance(&[GroupEsurfaceId(0)], &center, &fixture.external_momenta)
                    .unwrap()
                    > 0.49999
            );
        }
    }

    #[test]
    fn box_4e_objectives_preserve_four_certified_groups_and_proportional_clearance() {
        // arXiv:1912.09291, Eq. (3.10). Each one-loop surface has L_s=2, so the
        // first two objectives are proportional; a nonunique optimum need not
        // produce bit-identical coordinates. MinSum selects along the optimal face.
        let fixture = HelperBoxStructure::new(None);
        let mut results = Vec::new();
        let scale = 100.0;
        for objective in [
            OverlapCenterObjective::MaxMinDepth,
            OverlapCenterObjective::RelaxedChebyshev,
            OverlapCenterObjective::MinSum,
        ] {
            let mut settings = RuntimeSettings::default();
            settings.kinematics.e_cm = scale;
            settings.subtraction.overlap_settings.enable_heuristics = false;
            settings.subtraction.overlap_settings.objective = objective;
            let input = OverlapInput {
                graph_data: ti_vec![SingleGraphOverlapData {
                    lmb: &fixture.lmb,
                    esurfaces: &fixture.esurfaces,
                    raised_data: &fixture.raised_data,
                    edge_masses: fixture.edge_masses.clone(),
                }],
                settings: &settings,
                group_esurface_map: (0..4)
                    .map(|id| ti_vec![Some(RaisedEsurfaceId(id))])
                    .collect(),
                local_esurface_exists: ti_vec![ti_vec![true; 4]],
            };
            for id in 0..4 {
                assert_eq!(
                    input.surface_lipschitz(GraphGroupPosition::from(0), RaisedEsurfaceId(id)),
                    2.0
                );
            }
            let overlap = find_maximal_overlap(
                &input,
                &fixture.existing_esurfaces,
                &fixture.external_momenta,
            )
            .unwrap();
            assert_eq!(overlap.overlap_groups.len(), 4);
            let clearances_and_sums = overlap
                .overlap_groups
                .iter()
                .map(|group| {
                    assert_eq!(group.existing_esurfaces.len(), 2);
                    assert_eq!(group.complement.len(), 2);
                    let surfaces = group
                        .existing_esurfaces
                        .iter()
                        .map(|id| overlap.existing_esurfaces[*id])
                        .collect_vec();
                    let clearance = input
                        .center_clearance(&surfaces, &group.center, &fixture.external_momenta)
                        .expect("every selected center must satisfy the physical interior guard");
                    let sum = surfaces
                        .iter()
                        .map(|id| {
                            fixture.esurfaces[EsurfaceID::from(id.0)]
                                .compute_from_momenta(
                                    &fixture.lmb,
                                    &fixture.edge_masses,
                                    &group.center,
                                    &fixture.external_momenta,
                                )
                                .0
                        })
                        .sum::<f64>();
                    (clearance, sum)
                })
                .collect_vec();
            results.push((overlap, clearances_and_sums));
        }
        let tolerances = DefaultSettings::<f64>::default();
        let mut changed_centers = 0;
        for index in 0..4 {
            let baseline = &results[0].0.overlap_groups[index];
            let chebyshev = &results[1].0.overlap_groups[index];
            let min_sum = &results[2].0.overlap_groups[index];
            assert_eq!(baseline.existing_esurfaces, chebyshev.existing_esurfaces);
            assert_eq!(baseline.existing_esurfaces, min_sum.existing_esurfaces);
            assert_eq!(baseline.complement, chebyshev.complement);
            assert_eq!(baseline.complement, min_sum.complement);
            let baseline_radius = results[0].1[index].0;
            let chebyshev_radius = results[1].1[index].0;
            let accuracy = (tolerances.tol_gap_abs + tolerances.tol_feas) * scale
                + tolerances.tol_gap_rel * chebyshev_radius;
            assert!((baseline_radius - chebyshev_radius).abs() <= accuracy);
            assert!(results[2].1[index].0 >= chebyshev_radius - accuracy);
            assert!(results[2].1[index].1 <= results[1].1[index].1 + accuracy);
            let displacement_squared = chebyshev
                .center
                .iter()
                .zip(min_sum.center.iter())
                .map(|(left, right)| (left - right).norm_squared().0)
                .sum::<f64>();
            changed_centers += usize::from(displacement_squared > 1.0e-6);
        }
        assert_eq!(
            changed_centers, 4,
            "MinSum must select visibly different centers"
        );
    }

    #[test]
    fn test_box_4e() {
        // massless variant
        let box4e = HelperBoxStructure::new(None);

        let massless_overlap_input = OverlapInput {
            graph_data: ti_vec![SingleGraphOverlapData {
                lmb: &box4e.lmb,
                esurfaces: &box4e.esurfaces,
                raised_data: &box4e.raised_data,
                edge_masses: box4e.edge_masses.clone(),
            }],
            settings: &RuntimeSettings::default(),
            group_esurface_map: (0..4)
                .map(|i| ti_vec![Some(Into::<RaisedEsurfaceId>::into(i))])
                .collect(),
            local_esurface_exists: ti_vec![ti_vec![true; 4]],
        };

        let maximal_overlap = find_maximal_overlap(
            &massless_overlap_input,
            &box4e.existing_esurfaces,
            &box4e.external_momenta,
        )
        .unwrap();

        assert_eq!(maximal_overlap.overlap_groups.len(), 4);

        for overlap_group in maximal_overlap.overlap_groups.iter() {
            let esurfaces = &overlap_group.existing_esurfaces;
            let center = &overlap_group.center;

            assert_eq!(esurfaces.len(), 2);
            assert_eq!(overlap_group.complement.len(), 2);

            for esurface in esurfaces.iter() {
                let raised_esurface_id = massless_overlap_input.group_esurface_map
                    [box4e.existing_esurfaces[*esurface]][GraphGroupPosition::from(0)]
                .unwrap();
                let esurface_id =
                    box4e.raised_data.raised_groups[raised_esurface_id].esurface_ids[0];
                let esurfaec_val = box4e.esurfaces[esurface_id].compute_from_momenta(
                    &box4e.lmb,
                    &box4e.edge_masses,
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

        let overlap_input = OverlapInput {
            graph_data: ti_vec![SingleGraphOverlapData {
                lmb: &box4e.lmb,
                esurfaces: &box4e.esurfaces,
                raised_data: &box4e.raised_data,
                edge_masses: box4e.edge_masses.clone(),
            }],
            settings: &RuntimeSettings::default(),
            group_esurface_map: (0..4)
                .map(|i| ti_vec![Some(Into::<RaisedEsurfaceId>::into(i))])
                .collect(),
            local_esurface_exists: ti_vec![ti_vec![true; 4]],
        };

        let maximal_overlap = find_maximal_overlap(
            &overlap_input,
            &box4e.existing_esurfaces,
            &box4e.external_momenta,
        )
        .unwrap();

        assert_eq!(maximal_overlap.overlap_groups.len(), 4);

        for overlap_group in maximal_overlap.overlap_groups.iter() {
            let esurfaces = &overlap_group.existing_esurfaces;
            let center = &overlap_group.center;

            assert_eq!(esurfaces.len(), 1);

            for esurface in esurfaces.iter() {
                let raised_esurface_id = overlap_input.group_esurface_map
                    [box4e.existing_esurfaces[*esurface]][GraphGroupPosition::from(0)]
                .unwrap();
                let esurface_id =
                    box4e.raised_data.raised_groups[raised_esurface_id].esurface_ids[0];
                let esurfaec_val = box4e.esurfaces[esurface_id].compute_from_momenta(
                    &box4e.lmb,
                    &box4e.edge_masses,
                    center,
                    &box4e.external_momenta,
                );

                assert!(esurfaec_val < F(0.0));
            }

            assert_eq!(overlap_group.complement.len(), 3);
        }
    }

    #[test]
    fn test_banana() {
        let banana = HelperBananaStructure::new();

        let classification = banana.esurfaces[EsurfaceID::from(0)].classify_existence(
            &banana.external_momenta,
            &banana.lmb,
            &banana.edge_masses,
            &F(10.0),
            &F(crate::utils::DEFAULT_ESURFACE_EXISTENCE_THRESHOLD),
        );
        assert!(matches!(classification, EsurfaceExistence::Pinched { .. }));

        let overlap_input = OverlapInput {
            graph_data: ti_vec![SingleGraphOverlapData {
                lmb: &banana.lmb,
                esurfaces: &banana.esurfaces,
                raised_data: &banana.raised_data,
                edge_masses: banana.edge_masses.clone(),
            }],
            settings: &RuntimeSettings::default(),
            group_esurface_map: ti_vec![ti_vec![Some(Into::<RaisedEsurfaceId>::into(0)),]],
            local_esurface_exists: ti_vec![ti_vec![true]],
        };

        let result = find_maximal_overlap(
            &overlap_input,
            &banana.existing_esurfaces,
            &banana.external_momenta,
        );

        assert!(
            result.is_err(),
            "a caller that manually marks a pinched surface as existing must be rejected"
        );

        let mut forced_settings = RuntimeSettings::default();
        forced_settings
            .subtraction
            .overlap_settings
            .force_global_center = Some(vec![[0.0, 0.0, 0.0]; 2]);
        forced_settings
            .subtraction
            .overlap_settings
            .check_global_center = false;
        let forced_overlap_input = OverlapInput {
            graph_data: ti_vec![SingleGraphOverlapData {
                lmb: &banana.lmb,
                esurfaces: &banana.esurfaces,
                raised_data: &banana.raised_data,
                edge_masses: banana.edge_masses.clone(),
            }],
            settings: &forced_settings,
            group_esurface_map: ti_vec![ti_vec![Some(Into::<RaisedEsurfaceId>::into(0)),]],
            local_esurface_exists: ti_vec![ti_vec![true]],
        };
        assert!(
            find_maximal_overlap(
                &forced_overlap_input,
                &banana.existing_esurfaces,
                &banana.external_momenta,
            )
            .is_err(),
            "check_global_center=false must not allow a forced center to bypass the invariant"
        );
    }
}
