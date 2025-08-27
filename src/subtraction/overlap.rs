use crate::cff::esurface::get_representative;
use crate::cff::esurface::EsurfaceCollection;
use crate::cff::esurface::EsurfaceID;
use crate::cff::esurface::ExistingEsurfaceId;
use crate::cff::esurface::ExistingEsurfaces;
use crate::cff::esurface::GroupEsurfaceId;
use crate::cff::esurface::OldExistingEsurfaces;
use crate::gammaloop_integrand::GenericEvaluator;
use crate::graph::Graph;
use crate::graph::GraphGroupPosition;
use crate::graph::LoopMomentumBasis;
use crate::momentum::ThreeMomentum;
use crate::momentum_sample::ExternalFourMomenta;
use crate::momentum_sample::LoopMomenta;
use crate::settings::RuntimeSettings;
use crate::signature::LoopExtSignature;
use crate::utils::compute_shift_part;
use crate::utils::F;
use crate::utils::GS;
use crate::GammaLoopContext;
use ahash::HashMap;
use ahash::HashMapExt;
use ahash::HashSet;
use bincode_trait_derive::Decode;
use bincode_trait_derive::Encode;
use clarabel::algebra::*;
use clarabel::solver::*;
use eyre::{eyre, Result};
use itertools::Itertools;
use linnet::half_edge::involution::EdgeVec;
use linnet::half_edge::swap::Swap;
use serde::Serialize;
use serde_with::serde_as;
use spenso::algebra::algebraic_traits::IsZero;
use spenso::algebra::complex::Complex;
use std::fmt::Display;
use symbolica::atom::Atom;
use symbolica::evaluate::FunctionMap;
use symbolica::evaluate::OptimizationSettings;
use symbolica::function;
use tracing::debug;
use typed_index_collections::TiVec;

#[derive(Debug, Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct OverlapGroup {
    pub existing_esurfaces: Vec<ExistingEsurfaceId>,
    pub complement: Vec<ExistingEsurfaceId>,
    pub center: LoopMomenta<F<f64>>,
    pub prefactor_evaluator: Option<GenericEvaluator>,
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
        mass_atoms: Vec<Atom>,
    ) -> Result<()> {
        let group_square_atoms = self
            .overlap_groups
            .iter()
            .map(|group| {
                group
                    .existing_esurfaces
                    .iter()
                    .map(|&existing_esurface_id| {
                        let esurface = self.existing_esurfaces[existing_esurface_id];
                        let atom = &atoms[esurface];
                        atom * atom
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
                (1..=3).map(move |spatial_index| {
                    function!(GS.loop_mom, loop_index as i32, spatial_index)
                })
            })
            .chain((0..num_externals).flat_map(|external_index| {
                (0..=3).map(move |spatial_index| {
                    function!(GS.external_mom, external_index as i32, spatial_index)
                })
            }))
            .chain(mass_atoms)
            .collect_vec();

        for (group, square_atom) in self.overlap_groups.iter_mut().zip(group_square_atoms) {
            let atom = square_atom / &denominator;

            let evalautor = GenericEvaluator::new_from_raw_params(
                [atom],
                &params,
                &FunctionMap::new(),
                optimization_settings.clone(),
            )
            .ok_or_else(|| eyre!("Could not build evaluator for overlap prefactor"))?;

            group.prefactor_evaluator = Some(evalautor);
        }

        Ok(())
    }

    pub fn new_empty() -> Self {
        Self {
            overlap_groups: vec![],
            existing_esurfaces: ExistingEsurfaces::new(),
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

fn construct_solver(
    overlap_input: &OverlapInput,
    esurfaces_to_consider: &[ExistingEsurfaceId],
    existing_esurfaces: &ExistingEsurfaces,
    external_momenta: &ExternalFourMomenta<F<f64>>,
    verbose: bool,
) -> DefaultSolver {
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

    for esurface_id in esurfaces_to_consider.iter() {
        let group_esurface_id = existing_esurfaces[*esurface_id];
        let (graph_group_pos, esurface_id) =
            get_representative(&overlap_input.group_esurface_map[group_esurface_id])
                .expect("esurface map corrupted");

        let esurface = &overlap_input.graph_data[graph_group_pos].esurfaces[esurface_id];
        let lmb = overlap_input.graph_data[graph_group_pos].lmb;
        let edge_masses = overlap_input.graph_data[graph_group_pos].edge_masses;

        let mut esurface_constraint_indices: Vec<usize> = Vec::with_capacity(6);

        for &edge_id in &esurface.energies {
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
    for (constraint_index, (esurface_id, esurface_constraint)) in esurfaces_to_consider
        .iter()
        .zip(esurface_constraints.iter())
        .enumerate()
    {
        for prop_index in esurface_constraint {
            a_matrix[constraint_index + 1][*prop_index + propagator_index_offset] = 1.0;
        }
        let (graph_group_pos, esurface_id) =
            get_representative(&overlap_input.group_esurface_map[existing_esurfaces[*esurface_id]])
                .expect("esurface map corrupted");
        let lmb = overlap_input.graph_data[graph_group_pos].lmb;
        let esurface = &overlap_input.graph_data[graph_group_pos].esurfaces[esurface_id];

        let shift_part = esurface.compute_shift_part_from_momenta(lmb, external_momenta);
        b_vector[constraint_index + 1] = -shift_part.0;
        a_matrix[constraint_index + 1][0] = -1.0;
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

    debug!(name:"overlap_a_matrix", value = ?a_matrix);
    debug!(name:"overlap_b_vector", value = ?b_vector);

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
    existing_esurfaces: &ExistingEsurfaces,
    external_momenta: &ExternalFourMomenta<F<f64>>,
    verbose: bool,
) -> Option<LoopMomenta<F<f64>>> {
    let mut solver = construct_solver(
        overlap_input,
        esurfaces_to_consider,
        existing_esurfaces,
        external_momenta,
        verbose,
    );

    solver.solve();

    let loop_number = overlap_input
        .graph_data
        .first()
        .expect("no graphs passed to overlap")
        .lmb
        .loop_edges
        .len();

    if solver.solution.status == SolverStatus::Solved {
        let center = extract_center(loop_number, &solver.solution.x);
        Some(center)
    } else if solver.solution.status == SolverStatus::AlmostSolved
        || solver.solution.status == SolverStatus::InsufficientProgress
    {
        // if the solver did not converge, we check if the solution is still valid
        let center = extract_center(loop_number, &solver.solution.x);

        let is_valid = esurfaces_to_consider.iter().all(|&existing_esurface_id| {
            let (graph_group_pos, esurface_id) = get_representative(
                &overlap_input.group_esurface_map[existing_esurfaces[existing_esurface_id]],
            )
            .expect("overlap corrupted");

            let lmb = overlap_input.graph_data[graph_group_pos].lmb;
            let edge_masses = overlap_input.graph_data[graph_group_pos].edge_masses;
            let esurface = &overlap_input.graph_data[graph_group_pos].esurfaces[esurface_id];
            esurface.compute_from_momenta(lmb, &edge_masses, &center, external_momenta)
                < F::from_f64(0.0)
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

pub struct SingleGraphOverlapData<'a> {
    pub lmb: &'a LoopMomentumBasis,
    pub esurfaces: &'a EsurfaceCollection,
    pub edge_masses: &'a EdgeVec<F<f64>>,
}

pub struct OverlapInput<'a> {
    pub graph_data: TiVec<GraphGroupPosition, SingleGraphOverlapData<'a>>,
    pub settings: &'a RuntimeSettings,
    pub group_esurface_map: TiVec<GroupEsurfaceId, TiVec<GraphGroupPosition, Option<EsurfaceID>>>,
}

/// TODO: When this function will be called at runtime, panics should be removed and this function should return result.
/// When the overlap finding fails, treat the point as unstable
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

        if settings.subtraction.overlap_settings.check_global_center {
            let is_valid = existing_esurfaces.iter().all(|existing_esurface_id| {
                let (graph_group_postition, esurface_id) =
                    get_representative(&overlap_input.group_esurface_map[*existing_esurface_id])
                        .expect("overlap corrupted");
                let esurface =
                    &overlap_input.graph_data[graph_group_postition].esurfaces[esurface_id];

                let lmb = overlap_input.graph_data[graph_group_postition].lmb;
                let edge_masses = overlap_input.graph_data[graph_group_postition].edge_masses;

                let esurface_val = esurface.compute_from_momenta(
                    lmb,
                    &edge_masses,
                    &global_center_f,
                    external_momenta,
                );

                esurface_val < F(0.0)
            });

            if !is_valid {
                return Err(eyre!(
                    "Center provided is not inside all existing esurfaces"
                ));
            }
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

    if settings.subtraction.overlap_settings.try_origin
        || settings.subtraction.overlap_settings.try_origin_all_lmbs
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
    );

    if let Some(center) = option_center {
        let single_group = OverlapGroup {
            existing_esurfaces: all_existing_esurfaces,
            center,
            complement: vec![],
            prefactor_evaluator: None,
        };
        res.overlap_groups.push(single_group);
        res.fill_in_complements();
        return Ok(res);
    }

    // if the center is not valid, create a table of all pairs
    let esurface_pairs = EsurfacePairs::new(overlap_input, existing_esurfaces, external_momenta);

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
            )
            .ok_or_else(|| {
                let (graph_group_pos, esurface_id) =
                    get_representative(&overlap_input.group_esurface_map[esurface_id])
                        .expect("overlap corrupted");

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
                center: center,
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
                external_momenta,
                false,
            );

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
    Ok(res)
}

fn is_subset_of_result(subset: &[ExistingEsurfaceId], result: &OverlapStructure) -> bool {
    result.overlap_groups.iter().any(|group| {
        subset
            .iter()
            .all(|&x| group.existing_esurfaces.contains(&x))
    })
}

#[serde_as]
#[derive(Debug, Serialize)]
struct EsurfacePairs {
    #[serde_as(as = "Vec<(_,_)>")]
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

#[cfg(test)]
#[allow(dead_code, unused_variables)]
mod tests {
    use super::*;
    use itertools::Itertools;
    use linnet::half_edge::involution::EdgeIndex;
    use typed_index_collections::ti_vec;

    use crate::{
        cff::{
            cff_graph::VertexSet,
            esurface::{Esurface, EsurfaceID},
        },
        graph::LoopMomentumBasis,
        momentum::FourMomentum,
        settings::RuntimeSettings,
        signature::LoopExtSignature,
        utils::dummy_hedge_graph,
    };

    struct HelperBoxStructure {
        external_momenta: ExternalFourMomenta<F<f64>>,
        lmb: LoopMomentumBasis,
        esurfaces: EsurfaceCollection,
        existing_esurfaces: ExistingEsurfaces,
        edge_masses: EdgeVec<F<f64>>,
    }

    struct HelperBananaStructure {
        external_momenta: ExternalFourMomenta<F<f64>>,
        lmb: LoopMomentumBasis,
        esurfaces: EsurfaceCollection,
        existing_esurfaces: ExistingEsurfaces,
        edge_masses: EdgeVec<F<f64>>,
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
                tree: None,
                ext_edges: vec![].into(),
                loop_edges: box_basis,
                edge_signatures: box_signatures,
            };

            let esurfaces_array = [
                Esurface {
                    energies: vec![EdgeIndex::from(5), EdgeIndex::from(6)],
                    external_shift: vec![(EdgeIndex::from(1), 1)],
                    vertex_set: VertexSet::dummy(),
                },
                Esurface {
                    energies: vec![EdgeIndex::from(5), EdgeIndex::from(7)],
                    external_shift: vec![(EdgeIndex::from(1), 1), (EdgeIndex::from(2), 1)],
                    vertex_set: VertexSet::dummy(),
                },
                Esurface {
                    energies: vec![EdgeIndex::from(4), EdgeIndex::from(6)],
                    external_shift: vec![(EdgeIndex::from(0), 1), (EdgeIndex::from(1), 1)],
                    vertex_set: VertexSet::dummy(),
                },
                Esurface {
                    energies: vec![EdgeIndex::from(4), EdgeIndex::from(7)],
                    external_shift: vec![
                        (EdgeIndex::from(0), 1),
                        (EdgeIndex::from(1), 1),
                        (EdgeIndex::from(2), 1),
                    ],
                    vertex_set: VertexSet::dummy(),
                },
            ];

            let esurfaces = esurfaces_array.to_vec().into();

            let edge_masses = match masses {
                Some(masses) => {
                    let mut edge_masses = vec![F(0.0); 4];
                    let mut real_masses = masses.iter().map(|&x| x).collect_vec();

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
                tree: None,
                loop_edges: banana_basis,
                ext_edges: vec![].into(),
                edge_signatures: banana_edge_sigs,
            };

            let only_esurface = Esurface {
                energies: vec![EdgeIndex::from(2), EdgeIndex::from(3), EdgeIndex::from(4)],
                external_shift: vec![(EdgeIndex::from(0), -1)],
                vertex_set: VertexSet::dummy(),
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
                edge_masses: &box4e.edge_masses,
            }],
            settings: &RuntimeSettings::default(),
            group_esurface_map: (0..4)
                .map(|i| ti_vec![Some(Into::<EsurfaceID>::into(i))])
                .collect(),
        };

        let esurface_pairs = EsurfacePairs::new(
            &massless_overlap_input,
            &box4e.existing_esurfaces,
            &box4e.external_momenta,
        );

        assert_eq!(esurface_pairs.data.len(), 4);

        let box4e_massive = HelperBoxStructure::new(Some([F(10.5); 4]));

        let massive_overlap_input = OverlapInput {
            graph_data: ti_vec![SingleGraphOverlapData {
                lmb: &box4e_massive.lmb,
                esurfaces: &box4e_massive.esurfaces,
                edge_masses: &box4e_massive.edge_masses,
            }],
            settings: &RuntimeSettings::default(),
            group_esurface_map: (0..4)
                .map(|i| ti_vec![Some(Into::<EsurfaceID>::into(i))])
                .collect(),
        };

        let esurface_pairs_massive = EsurfacePairs::new(
            &massive_overlap_input,
            &box4e_massive.existing_esurfaces,
            &box4e_massive.external_momenta,
        );

        assert_eq!(esurface_pairs_massive.data.len(), 0);
    }

    #[test]
    fn test_subset_generator() {
        let box4e = HelperBoxStructure::new(None);

        let overlap_input = OverlapInput {
            graph_data: ti_vec![SingleGraphOverlapData {
                lmb: &box4e.lmb,
                esurfaces: &box4e.esurfaces,
                edge_masses: &box4e.edge_masses,
            }],
            settings: &RuntimeSettings::default(),
            group_esurface_map: (0..4)
                .map(|i| ti_vec![Some(Into::<EsurfaceID>::into(i))])
                .collect(),
        };

        let esurface_pairs = EsurfacePairs::new(
            &overlap_input,
            &box4e.existing_esurfaces,
            &box4e.external_momenta,
        );

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
    fn test_box_4e() {
        // massless variant
        let box4e = HelperBoxStructure::new(None);

        let massless_overlap_input = OverlapInput {
            graph_data: ti_vec![SingleGraphOverlapData {
                lmb: &box4e.lmb,
                esurfaces: &box4e.esurfaces,
                edge_masses: &box4e.edge_masses,
            }],
            settings: &RuntimeSettings::default(),
            group_esurface_map: (0..4)
                .map(|i| ti_vec![Some(Into::<EsurfaceID>::into(i))])
                .collect(),
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
                let esurfaec_val = box4e.esurfaces[massless_overlap_input.group_esurface_map
                    [box4e.existing_esurfaces[*esurface]][GraphGroupPosition::from(0)]
                .unwrap()]
                .compute_from_momenta(
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
                edge_masses: &box4e.edge_masses,
            }],
            settings: &RuntimeSettings::default(),
            group_esurface_map: (0..4)
                .map(|i| ti_vec![Some(Into::<EsurfaceID>::into(i))])
                .collect(),
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
                let esurfaec_val = box4e.esurfaces[overlap_input.group_esurface_map
                    [box4e.existing_esurfaces[*esurface]][GraphGroupPosition::from(0)]
                .unwrap()]
                .compute_from_momenta(
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

        let overlap_input = OverlapInput {
            graph_data: ti_vec![SingleGraphOverlapData {
                lmb: &banana.lmb,
                esurfaces: &banana.esurfaces,
                edge_masses: &banana.edge_masses,
            }],
            settings: &RuntimeSettings::default(),
            group_esurface_map: ti_vec![ti_vec![Some(Into::<EsurfaceID>::into(0)),]],
        };

        let maximal_overlap = find_maximal_overlap(
            &overlap_input,
            &banana.existing_esurfaces,
            &banana.external_momenta,
        )
        .unwrap();

        assert_eq!(maximal_overlap.overlap_groups.len(), 1);
        assert_eq!(
            maximal_overlap.overlap_groups[0].existing_esurfaces.len(),
            1
        );
        assert!(maximal_overlap.overlap_groups[0].complement.is_empty());
    }
}
