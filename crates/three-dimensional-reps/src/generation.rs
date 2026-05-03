use std::{
    collections::{BTreeMap, BTreeSet, HashMap},
    ops::{Add, Neg, Sub},
};

use itertools::Itertools;
use linnet::half_edge::involution::{EdgeIndex, EdgeVec, Orientation};
use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::{
    ContourClosure, OrientationData, OrientationExpression, OrientationID, ParsedGraph,
    ThreeDExpression,
    cff_recursion::enumerate_cff_surface_chains,
    cut_structure::ltd_residues,
    energy_bounds::{energy_divergence_report, normalize_energy_degree_bounds},
    expression::{ResidualDenominator, assign_numerator_map_labels},
    graph_io::{
        ParsedGraphExternalEdge, ParsedGraphInitialStateCutEdge, ParsedGraphInternalEdge,
        ThreeDGraphSource, repeated_groups,
    },
    graph_signatures::{
        ExtractedSignatureExpression, MomentumSignature, ReconstructDotFormat,
        ReconstructDotOptions, reconstruct_parsed_graph,
    },
    surface::{
        HybridSurfaceID, LinearEnergyExpr, LinearSurface, LinearSurfaceID, LinearSurfaceKind,
        RationalAtomExt, SurfaceOrigin, rational_coeff_atom, rational_coeff_new,
        rational_coeff_one,
    },
    tree::{NodeId, Tree},
    utils::{
        Rational, RationalExt, binomial, determinant_i32_is_nonzero, factorial, multi_factorial,
        multiindices_leq, rank_i64, rank_rational, rational_pow_i64, rising, solve_rational_system,
    },
};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum RepresentationMode {
    Ltd,
    Cff,
    #[cfg(any(feature = "diagnostics", feature = "test-support"))]
    PureLtd,
    #[cfg(feature = "old_cff")]
    OldCff,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum NumeratorSamplingScaleMode {
    #[default]
    None,
    BeyondQuadratic,
    All,
}

impl NumeratorSamplingScaleMode {
    pub const fn is_active_for_degree(self, degree: usize) -> bool {
        match self {
            Self::None => false,
            Self::BeyondQuadratic => degree > 2,
            Self::All => degree > 0,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Generate3DExpressionOptions {
    pub representation: RepresentationMode,
    pub energy_degree_bounds: Vec<(usize, usize)>,
    pub numerator_sampling_scale: NumeratorSamplingScaleMode,
    #[serde(default)]
    pub preserve_internal_edges_as_four_d_denominators: Vec<usize>,
}

impl Default for Generate3DExpressionOptions {
    fn default() -> Self {
        Self {
            representation: RepresentationMode::Cff,
            energy_degree_bounds: Vec::new(),
            numerator_sampling_scale: NumeratorSamplingScaleMode::None,
            preserve_internal_edges_as_four_d_denominators: Vec::new(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GenerationDiagnostic {
    pub warnings: Vec<GenerationWarning>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum GenerationWarning {
    PureLtdRepeatedPropagatorsRequireNumeratorDerivatives,
}

#[derive(Debug, Error)]
pub enum GenerationError {
    #[error("generation mode {mode:?} is unavailable for this input API or feature set")]
    NotImplemented { mode: RepresentationMode },
    #[error("old_cff only supports the legacy affine/linear energy-numerator regime")]
    OldCffHigherEnergyPower,
    #[error(
        "this generalized CFF higher energy-numerator sector is not supported by the current Rust port"
    )]
    CffHigherEnergyPowerNotImplemented,
    #[error("cut-structure generation failed: {0}")]
    CutStructure(#[from] crate::cut_structure::CutStructureError),
    #[error("could not find a nonsingular loop-energy basis")]
    SingularBasis,
    #[error("loop-energy solve produced non-integral coefficients")]
    NonIntegralEnergyMap,
    #[error("rational coefficient is outside the supported i64 storage range")]
    CoefficientOutOfRange,
    #[error("physical numerator sample orbit is not affine-complete")]
    NumeratorSampleOrbitIncomplete,
    #[error("{0}")]
    EnergyBounds(#[from] crate::energy_bounds::EnergyBoundsError),
    #[error("{0}")]
    GraphIo(#[from] crate::graph_io::GraphIoError),
    #[error("energy-degree bounds leave a non-vanishing residue at infinity")]
    EnergyDegreeBoundsNonConvergent,
    #[error(
        "energy-degree bound was requested for edge {edge_id}, but this is not an internal edge of the selected graph"
    )]
    UnknownEnergyDegreeBoundEdge { edge_id: usize },
    #[error(
        "CFF preservation was requested for edge {edge_id}, but this is not an internal edge of the selected graph"
    )]
    UnknownCffPreservedInternalEdge { edge_id: usize },
    #[error(
        "energy-degree bound was requested for preserved 4D denominator edge {edge_id}; preserved edges stay in the four-dimensional denominator and cannot receive CFF energy-degree bounds"
    )]
    CffPreservedEdgeHasEnergyDegreeBound { edge_id: usize },
}

pub type Result<T> = std::result::Result<T, GenerationError>;

pub fn generate_3d_expression<G: ThreeDGraphSource + ?Sized>(
    graph: &G,
    options: &Generate3DExpressionOptions,
) -> Result<ThreeDExpression<OrientationID>> {
    let parsed = graph.to_three_d_parsed_graph()?;
    let Some(edge_map) = graph.energy_edge_index_map(&parsed) else {
        return generate_3d_expression_from_parsed(&parsed, options);
    };
    let mut local_options = options.clone();
    local_options.energy_degree_bounds = edge_map
        .remap_bounds_to_local(&options.energy_degree_bounds)
        .map_err(|edge_id| GenerationError::UnknownEnergyDegreeBoundEdge { edge_id })?;
    let reverse_internal = edge_map.internal_to_local();
    local_options.preserve_internal_edges_as_four_d_denominators = options
        .preserve_internal_edges_as_four_d_denominators
        .iter()
        .map(|edge_id| {
            reverse_internal
                .get(edge_id)
                .copied()
                .ok_or(GenerationError::UnknownCffPreservedInternalEdge { edge_id: *edge_id })
        })
        .collect::<Result<Vec<_>>>()?;
    let mut expression = generate_3d_expression_from_parsed(&parsed, &local_options)?
        .remap_energy_edge_indices(&edge_map)
        .fuse_compatible_variants();
    assign_numerator_map_labels(&mut expression.orientations);
    Ok(expression)
}

#[allow(unreachable_patterns)]
pub fn generate_3d_expression_from_parsed(
    parsed: &ParsedGraph,
    options: &Generate3DExpressionOptions,
) -> Result<ThreeDExpression<OrientationID>> {
    if options.representation == RepresentationMode::Cff
        && !options
            .preserve_internal_edges_as_four_d_denominators
            .is_empty()
    {
        return build_cff_expression_preserving_internal_edges(parsed, options);
    }

    if !options.energy_degree_bounds.is_empty() {
        let signatures = parsed
            .internal_edges
            .iter()
            .map(|edge| edge.signature.clone())
            .collect::<Vec<_>>();
        let report = energy_divergence_report(&signatures, &options.energy_degree_bounds)?;
        if !report.convergent {
            return Err(GenerationError::EnergyDegreeBoundsNonConvergent);
        }
    }

    let expression = match options.representation {
        RepresentationMode::Ltd => generate_ltd_expression_from_parsed(parsed, options),
        RepresentationMode::Cff if has_higher_energy_power(options) => {
            BoundedCffBuilder::new(parsed, options)?.build()
        }
        RepresentationMode::Cff => generate_pure_cff_expression_from_parsed(parsed),
        #[cfg(any(feature = "diagnostics", feature = "test-support"))]
        RepresentationMode::PureLtd => generate_pure_ltd_expression_from_parsed(parsed),
        #[cfg(feature = "old_cff")]
        RepresentationMode::OldCff if has_higher_energy_power(options) => {
            Err(GenerationError::OldCffHigherEnergyPower)
        }
        #[cfg(feature = "old_cff")]
        RepresentationMode::OldCff => generate_pure_cff_expression_from_parsed(parsed),
        mode => Err(GenerationError::NotImplemented { mode }),
    }?;
    let mut expression = expression.fuse_compatible_variants();
    assign_numerator_map_labels(&mut expression.orientations);
    Ok(expression)
}

fn generate_pure_cff_expression_from_parsed(
    parsed: &ParsedGraph,
) -> Result<ThreeDExpression<OrientationID>> {
    generate_pure_cff_expression_from_parsed_with_duplicate_sign(parsed, true)
}

fn build_cff_expression_preserving_internal_edges(
    parsed: &ParsedGraph,
    options: &Generate3DExpressionOptions,
) -> Result<ThreeDExpression<OrientationID>> {
    let preserved = options
        .preserve_internal_edges_as_four_d_denominators
        .iter()
        .copied()
        .collect::<BTreeSet<_>>();
    for edge_id in &preserved {
        if *edge_id >= parsed.internal_edges.len() {
            return Err(GenerationError::UnknownCffPreservedInternalEdge { edge_id: *edge_id });
        }
    }

    let (active_parsed, active_to_orig) = contract_preserved_parsed_edges(parsed, &preserved);
    if active_parsed.internal_edges.is_empty() {
        return cff_expression_with_only_preserved_edges(parsed);
    }

    let orig_to_active = active_to_orig
        .iter()
        .enumerate()
        .map(|(active_id, orig_id)| (*orig_id, active_id))
        .collect::<BTreeMap<_, _>>();
    let mut active_options = options.clone();
    active_options
        .preserve_internal_edges_as_four_d_denominators
        .clear();
    active_options.energy_degree_bounds = options
        .energy_degree_bounds
        .iter()
        .filter_map(|(edge_id, degree)| {
            if *edge_id >= parsed.internal_edges.len() {
                return Some(Err(GenerationError::UnknownEnergyDegreeBoundEdge {
                    edge_id: *edge_id,
                }));
            }
            if preserved.contains(edge_id) {
                if *degree > 0 {
                    return Some(Err(GenerationError::CffPreservedEdgeHasEnergyDegreeBound {
                        edge_id: *edge_id,
                    }));
                }
                return None;
            }
            Some(
                orig_to_active
                    .get(edge_id)
                    .copied()
                    .map(|active_id| (active_id, *degree))
                    .ok_or(GenerationError::UnknownEnergyDegreeBoundEdge { edge_id: *edge_id }),
            )
        })
        .collect::<Result<Vec<_>>>()?;

    let active_expression = generate_3d_expression_from_parsed(&active_parsed, &active_options)?;
    lift_cff_expression_to_preserved_graph(parsed, &active_expression, &active_to_orig, &preserved)
}

fn cff_expression_with_only_preserved_edges(
    parsed: &ParsedGraph,
) -> Result<ThreeDExpression<OrientationID>> {
    if !parsed.loop_names.is_empty() {
        return Err(GenerationError::SingularBasis);
    }

    let signatures = parsed
        .internal_edges
        .iter()
        .map(|edge| edge.signature.clone())
        .collect::<Vec<_>>();
    let edge_energy_map = edge_q0_from_loop_exprs(&signatures, &[]);
    let orientation =
        EdgeVec::from_iter((0..parsed.internal_edges.len()).map(|_| Orientation::Undirected));
    let mut data = OrientationData::new(orientation);
    data.label = Some("preserved_tree".to_string());

    let mut expression = ThreeDExpression::<OrientationID>::new_empty();
    expression.orientations.push(OrientationExpression {
        data,
        loop_energy_map: Vec::new(),
        edge_energy_map,
        variants: vec![crate::expression::CFFVariant {
            origin: Some("preserved_tree".to_string()),
            prefactor: rational_coeff_one(),
            half_edges: Vec::new(),
            uniform_scale_power: 0,
            numerator_surfaces: Vec::new(),
            denominator: Tree::from_root(HybridSurfaceID::Unit),
        }],
    });
    expression.residual_denominators = parsed
        .internal_edges
        .iter()
        .enumerate()
        .map(|(edge_id, edge)| {
            ResidualDenominator::new(EdgeIndex(edge_id), Some(edge.label.clone()))
        })
        .collect();
    assign_numerator_map_labels(&mut expression.orientations);
    Ok(expression)
}

fn contract_preserved_parsed_edges(
    parsed: &ParsedGraph,
    preserved: &BTreeSet<usize>,
) -> (ParsedGraph, Vec<usize>) {
    let mut parent = parsed
        .node_name_to_internal
        .values()
        .copied()
        .map(|node| (node, node))
        .collect::<BTreeMap<_, _>>();

    for (edge_index, edge) in parsed.internal_edges.iter().enumerate() {
        if preserved.contains(&edge_index) {
            union_nodes(&mut parent, edge.tail, edge.head);
        }
    }

    let parent_keys = parent.keys().copied().collect::<Vec<_>>();
    let roots = parent_keys
        .iter()
        .map(|node| find_node_root(&mut parent, *node))
        .collect::<BTreeSet<_>>();
    let root_to_new = roots
        .into_iter()
        .enumerate()
        .map(|(new_id, root)| (root, new_id))
        .collect::<BTreeMap<_, _>>();
    let old_to_new = parent_keys
        .iter()
        .map(|node| {
            let root = find_node_root(&mut parent, *node);
            (*node, root_to_new[&root])
        })
        .collect::<BTreeMap<_, _>>();

    let mut active_to_orig = Vec::new();
    let mut internal_edges = Vec::new();
    for (edge_index, edge) in parsed.internal_edges.iter().enumerate() {
        if preserved.contains(&edge_index) {
            continue;
        }
        let active_id = internal_edges.len();
        active_to_orig.push(edge_index);
        internal_edges.push(ParsedGraphInternalEdge {
            edge_id: active_id,
            tail: old_to_new[&edge.tail],
            head: old_to_new[&edge.head],
            label: edge.label.clone(),
            mass_key: edge.mass_key.clone(),
            signature: edge.signature.clone(),
            had_pow: edge.had_pow,
        });
    }

    let external_edges = parsed
        .external_edges
        .iter()
        .map(|edge| ParsedGraphExternalEdge {
            edge_id: edge.edge_id,
            source: edge.source.map(|source| old_to_new[&source]),
            destination: edge.destination.map(|destination| old_to_new[&destination]),
            label: edge.label.clone(),
            external_coefficients: edge.external_coefficients.clone(),
        })
        .collect::<Vec<_>>();
    let orig_to_active = active_to_orig
        .iter()
        .enumerate()
        .map(|(active_id, orig_id)| (*orig_id, active_id))
        .collect::<BTreeMap<_, _>>();
    let initial_state_cut_edges = parsed
        .initial_state_cut_edges
        .iter()
        .filter_map(|cut_edge| {
            orig_to_active
                .get(&cut_edge.edge_id)
                .copied()
                .map(|edge_id| ParsedGraphInitialStateCutEdge {
                    edge_id,
                    external_id: cut_edge.external_id,
                    external_sign: cut_edge.external_sign,
                })
        })
        .collect::<Vec<_>>();
    let node_name_to_internal = old_to_new
        .values()
        .copied()
        .collect::<BTreeSet<_>>()
        .into_iter()
        .map(|node| (format!("p{node}"), node))
        .collect();

    (
        ParsedGraph {
            internal_edges,
            external_edges,
            initial_state_cut_edges,
            loop_names: parsed.loop_names.clone(),
            external_names: parsed.external_names.clone(),
            node_name_to_internal,
        },
        active_to_orig,
    )
}

fn lift_cff_expression_to_preserved_graph(
    parsed: &ParsedGraph,
    source: &ThreeDExpression<OrientationID>,
    active_to_orig: &[usize],
    preserved: &BTreeSet<usize>,
) -> Result<ThreeDExpression<OrientationID>> {
    let active_edge_map = active_to_orig
        .iter()
        .enumerate()
        .map(|(active_id, orig_id)| (active_id, *orig_id))
        .collect::<BTreeMap<_, _>>();
    let signatures = parsed
        .internal_edges
        .iter()
        .map(|edge| edge.signature.clone())
        .collect::<Vec<_>>();

    let mut expression = ThreeDExpression::<OrientationID>::new_empty();
    expression.residual_denominators = preserved
        .iter()
        .map(|edge_id| {
            ResidualDenominator::new(
                EdgeIndex(*edge_id),
                Some(parsed.internal_edges[*edge_id].label.clone()),
            )
        })
        .collect();
    let surface_map = {
        let mut interner = LinearSurfaceInterner::new(&mut expression);
        source
            .surfaces
            .linear_surface_cache
            .iter_enumerated()
            .map(|(id, surface)| {
                (
                    HybridSurfaceID::Linear(id),
                    interner.intern_with_origin(
                        surface
                            .expression
                            .clone()
                            .remap_internal_edges(&active_edge_map),
                        surface.origin,
                        surface.numerator_only,
                    ),
                )
            })
            .collect::<HashMap<_, _>>()
    };

    for orientation in &source.orientations {
        let loop_energy_map = orientation
            .loop_energy_map
            .iter()
            .cloned()
            .map(|expr| expr.remap_internal_edges(&active_edge_map))
            .collect::<Vec<_>>();
        let mut edge_energy_map = edge_q0_from_loop_exprs(&signatures, &loop_energy_map);
        for (active_id, orig_id) in active_to_orig.iter().enumerate() {
            edge_energy_map[*orig_id] = orientation.edge_energy_map[active_id]
                .clone()
                .remap_internal_edges(&active_edge_map);
        }

        let (orientation_data, base_label) = orientation_from_edge_exprs(&edge_energy_map);
        let mut data = OrientationData::new(orientation_data);
        data.label = Some(base_label);
        let variants = orientation
            .variants
            .iter()
            .map(|variant| crate::expression::CFFVariant {
                origin: variant.origin.clone(),
                prefactor: variant.prefactor.clone(),
                half_edges: variant
                    .half_edges
                    .iter()
                    .map(|edge_id| EdgeIndex(active_edge_map[&edge_id.0]))
                    .collect(),
                uniform_scale_power: variant.uniform_scale_power,
                numerator_surfaces: variant
                    .numerator_surfaces
                    .iter()
                    .map(|surface_id| map_surface_id(*surface_id, &surface_map))
                    .collect(),
                denominator: variant
                    .denominator
                    .clone()
                    .map(|surface_id| map_surface_id(surface_id, &surface_map)),
            })
            .collect::<Vec<_>>();

        expression.orientations.push(OrientationExpression {
            data,
            loop_energy_map,
            edge_energy_map,
            variants,
        });
    }

    assign_numerator_map_labels(&mut expression.orientations);
    Ok(expression.fuse_compatible_variants())
}

fn generate_pure_cff_expression_from_parsed_with_duplicate_sign(
    parsed: &ParsedGraph,
    include_duplicate_excess_sign: bool,
) -> Result<ThreeDExpression<OrientationID>> {
    let signatures = parsed
        .internal_edges
        .iter()
        .map(|edge| edge.signature.clone())
        .collect::<Vec<_>>();
    let n_internal = signatures.len();
    let denominator_edge_ids = parsed.denominator_internal_edge_ids();
    let basis = choose_basis_indices_from_edges(&signatures, &denominator_edge_ids)?;
    let duplicate_excess = if include_duplicate_excess_sign {
        cff_duplicate_signature_excess(parsed)
    } else {
        0
    };
    let n_loops = signatures
        .first()
        .map(|signature| signature.loop_signature.len())
        .unwrap_or(0);
    let overall_sign = if (n_loops.saturating_sub(1) + duplicate_excess) % 2 == 0 {
        1
    } else {
        -1
    };

    let mut expression = ThreeDExpression::<OrientationID>::new_empty();
    let mut surface_index =
        HashMap::<(LinearSurfaceKind, LinearEnergyExpr), HybridSurfaceID>::new();

    for bitmask in 0..(1usize << denominator_edge_ids.len()) {
        let mut signs = vec![1; n_internal];
        for (bit, edge_index) in denominator_edge_ids.iter().copied().enumerate() {
            signs[edge_index] = if bitmask & (1usize << bit) == 0 {
                1
            } else {
                -1
            };
        }
        let surface_chains = enumerate_cff_surface_chains(parsed, &signs);
        if surface_chains.is_empty() {
            continue;
        }

        let mut edge_energy_map = vec![LinearEnergyExpr::zero(); n_internal];
        for edge_index in &denominator_edge_ids {
            edge_energy_map[*edge_index] =
                LinearEnergyExpr::ose(EdgeIndex(*edge_index), i64::from(signs[*edge_index]));
        }
        apply_initial_state_cut_edge_energy_exprs(parsed, &mut edge_energy_map);
        let loop_energy_map =
            solve_loop_energy_from_target_edge_exprs(&signatures, &basis, &edge_energy_map)?;
        let orientation = EdgeVec::from_iter((0..n_internal).map(|edge_index| {
            if parsed.is_initial_state_cut_edge(edge_index) || signs[edge_index] >= 0 {
                Orientation::Default
            } else {
                Orientation::Reversed
            }
        }));
        let data = OrientationData::new(orientation);
        let denominator_chains = surface_chains
            .into_iter()
            .map(|chain| {
                chain
                    .into_iter()
                    .map(|surface_expr| {
                        intern_linear_surface(
                            &mut expression,
                            &mut surface_index,
                            surface_expr,
                            false,
                        )
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let variants = vec![crate::expression::CFFVariant {
            origin: Some("pure_cff".to_string()),
            prefactor: rational_coeff_new(overall_sign, 1),
            half_edges: denominator_edge_ids
                .iter()
                .copied()
                .map(EdgeIndex)
                .collect(),
            uniform_scale_power: 0,
            numerator_surfaces: Vec::new(),
            denominator: denominator_tree_from_chains(&denominator_chains),
        }];

        expression.orientations.push(OrientationExpression {
            data,
            loop_energy_map,
            edge_energy_map,
            variants,
        });
    }

    Ok(expression)
}

fn generate_ltd_expression_from_parsed(
    parsed: &ParsedGraph,
    options: &Generate3DExpressionOptions,
) -> Result<ThreeDExpression<OrientationID>> {
    if !repeated_groups(parsed).is_empty() {
        return RepeatedLtdBuilder::new(parsed, options).build();
    }
    generate_pure_ltd_expression_from_parsed(parsed)
}

fn generate_pure_ltd_expression_from_parsed(
    parsed: &ParsedGraph,
) -> Result<ThreeDExpression<OrientationID>> {
    let signatures = parsed
        .internal_edges
        .iter()
        .map(|edge| edge.signature.clone())
        .collect::<Vec<_>>();
    let n_internal = signatures.len();
    let denominator_edge_ids = parsed.denominator_internal_edge_ids();
    let n_loops = signatures
        .first()
        .map(|signature| signature.loop_signature.len())
        .unwrap_or(0);
    let loop_line_signatures = denominator_edge_ids
        .iter()
        .map(|edge_id| signatures[*edge_id].loop_signature.clone())
        .collect::<Vec<_>>();
    let residues = ltd_residues(&loop_line_signatures, &vec![ContourClosure::Below; n_loops])?;

    let mut expression = ThreeDExpression::<OrientationID>::new_empty();
    let mut surface_index =
        HashMap::<(LinearSurfaceKind, LinearEnergyExpr), HybridSurfaceID>::new();

    for residue in residues {
        let basis = residue
            .basis
            .iter()
            .map(|edge| EdgeIndex(denominator_edge_ids[*edge]))
            .collect_vec();
        let cut_signs = residue.sigmas.clone();
        let basis_edges = residue
            .basis
            .iter()
            .map(|edge| denominator_edge_ids[*edge])
            .collect::<Vec<_>>();
        let (loop_energy_map, edge_energy_map) =
            solve_loop_energy_substitutions(parsed, &signatures, &basis_edges, &cut_signs)?;

        let basis_set = residue
            .basis
            .iter()
            .map(|edge| denominator_edge_ids[*edge])
            .collect::<std::collections::BTreeSet<_>>();
        let mut denominator_chain = Vec::new();
        for edge_index in &denominator_edge_ids {
            if basis_set.contains(edge_index) {
                continue;
            }
            let edge_expr = &edge_energy_map[*edge_index];
            let edge_id = EdgeIndex(*edge_index);
            let minus = edge_expr.clone() - LinearEnergyExpr::ose(edge_id, 1);
            let plus = edge_expr.clone() + LinearEnergyExpr::ose(edge_id, 1);
            denominator_chain.push(intern_linear_surface(
                &mut expression,
                &mut surface_index,
                minus,
                false,
            ));
            denominator_chain.push(intern_linear_surface(
                &mut expression,
                &mut surface_index,
                plus,
                false,
            ));
        }

        let mut orientation = EdgeVec::from_iter((0..n_internal).map(|_| Orientation::Undirected));
        for (edge_index, cut_sign) in basis_edges.iter().zip(&cut_signs) {
            orientation[EdgeIndex(*edge_index)] = if *cut_sign >= 0 {
                Orientation::Default
            } else {
                Orientation::Reversed
            };
        }
        let mut data = OrientationData::new(orientation);
        data.numerator_map_index = None;
        let denominator = denominator_tree_from_chain(&denominator_chain);
        let prefactor = if residue.sign >= 0 {
            rational_coeff_one()
        } else {
            rational_coeff_new(-1, 1)
        };
        let variant = crate::expression::CFFVariant {
            origin: Some("pure_ltd".to_string()),
            prefactor,
            half_edges: basis,
            uniform_scale_power: 0,
            numerator_surfaces: Vec::new(),
            denominator,
        };
        expression.orientations.push(OrientationExpression {
            data,
            loop_energy_map,
            edge_energy_map,
            variants: vec![variant],
        });
    }

    Ok(expression)
}

#[derive(Debug, Clone)]
struct LogicalChannel {
    rep_edge: usize,
    members: Vec<usize>,
    power: usize,
}

#[derive(Debug, Clone)]
struct PhysicalGroupOption {
    tau: Vec<i32>,
    chain_sign: i32,
    label: String,
}

#[derive(Debug, Clone)]
struct PhysicalSample {
    loop_exprs: Vec<LinearEnergyExpr>,
    edge_exprs: Vec<LinearEnergyExpr>,
    label: String,
}

#[derive(Debug, Clone)]
struct NumeratorSample {
    coeff: Rational,
    extra_half_edges: Vec<usize>,
    loop_exprs: Vec<LinearEnergyExpr>,
    edge_exprs: Vec<LinearEnergyExpr>,
    label: String,
    uniform_scale_power: usize,
}

#[derive(Debug, Clone)]
enum DenominatorFactor {
    Cut {
        edge: usize,
        power: usize,
        derivs: Vec<Rational>,
    },
    Surface {
        surface: HybridSurfaceID,
        power: usize,
        derivs: Vec<Rational>,
    },
}

#[derive(Debug, Clone)]
struct DenominatorTerm {
    coeff: Rational,
    half_edges: Vec<usize>,
    chain: Vec<HybridSurfaceID>,
}

#[derive(Debug, Clone)]
struct ContactComponent {
    sample: i32,
    prefactor: Rational,
    half_edges: Vec<usize>,
    numerator_surfaces: Vec<HybridSurfaceID>,
}

#[derive(Debug, Clone)]
struct ChannelNormalFormTerm {
    remaining_power: usize,
    parity: usize,
    cancelled_power: usize,
    inverse_power: usize,
    positive_ose_power: usize,
    coeff: Rational,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct KnownLinearExpr {
    var_terms: Vec<(usize, Rational)>,
    ose_terms: Vec<(usize, Rational)>,
    external_terms: Vec<(usize, Rational)>,
    constant: Rational,
    uniform_scale_coeff: Rational,
}

impl Default for KnownLinearExpr {
    fn default() -> Self {
        Self {
            var_terms: Vec::new(),
            ose_terms: Vec::new(),
            external_terms: Vec::new(),
            constant: Rational::zero(),
            uniform_scale_coeff: Rational::zero(),
        }
    }
}

impl KnownLinearExpr {
    fn zero() -> Self {
        Self::default()
    }

    fn var(edge_id: usize, coeff: i64) -> Self {
        Self::var_with_coeff(edge_id, Rational::from(coeff))
    }

    fn var_with_coeff(edge_id: usize, coeff: Rational) -> Self {
        if coeff.is_zero() {
            Self::zero()
        } else {
            Self {
                var_terms: vec![(edge_id, coeff)],
                ..Self::default()
            }
            .canonical()
        }
    }

    fn ose(edge_id: usize, coeff: i64) -> Self {
        Self::ose_with_coeff(edge_id, Rational::from(coeff))
    }

    fn ose_with_coeff(edge_id: usize, coeff: Rational) -> Self {
        if coeff.is_zero() {
            Self::zero()
        } else {
            Self {
                ose_terms: vec![(edge_id, coeff)],
                ..Self::default()
            }
            .canonical()
        }
    }

    fn external(edge_id: usize, coeff: i64) -> Self {
        Self::external_with_coeff(edge_id, Rational::from(coeff))
    }

    fn external_with_coeff(edge_id: usize, coeff: Rational) -> Self {
        if coeff.is_zero() {
            Self::zero()
        } else {
            Self {
                external_terms: vec![(edge_id, coeff)],
                ..Self::default()
            }
            .canonical()
        }
    }

    fn canonical(mut self) -> Self {
        self.var_terms = Self::collect_terms(self.var_terms);
        self.ose_terms = Self::collect_terms(self.ose_terms);
        self.external_terms = Self::collect_terms(self.external_terms);
        self
    }

    fn replace_var_with_ose(&self, edge_id: usize, sample: i32, ose_edge_id: usize) -> Self {
        let mut out = Self {
            ose_terms: self.ose_terms.clone(),
            external_terms: self.external_terms.clone(),
            constant: self.constant.clone(),
            uniform_scale_coeff: self.uniform_scale_coeff.clone(),
            ..Self::default()
        };
        for (var_edge_id, coeff) in &self.var_terms {
            if *var_edge_id == edge_id {
                out = out
                    + Self::ose_with_coeff(
                        ose_edge_id,
                        coeff.clone() * Rational::from(i64::from(sample)),
                    );
            } else {
                out = out + Self::var_with_coeff(*var_edge_id, coeff.clone());
            }
        }
        out.canonical()
    }

    fn mul_rational(&self, scale: Rational) -> Result<Self> {
        fn scale_terms(terms: &[(usize, Rational)], scale: &Rational) -> Vec<(usize, Rational)> {
            terms
                .iter()
                .filter_map(|(edge_id, coeff)| {
                    let value = coeff.clone() * scale.clone();
                    (!value.is_zero()).then_some((*edge_id, value))
                })
                .collect()
        }

        Ok(Self {
            var_terms: scale_terms(&self.var_terms, &scale),
            ose_terms: scale_terms(&self.ose_terms, &scale),
            external_terms: scale_terms(&self.external_terms, &scale),
            constant: self.constant.clone() * scale.clone(),
            uniform_scale_coeff: self.uniform_scale_coeff.clone() * scale,
        }
        .canonical())
    }

    fn is_zero(&self) -> bool {
        let item = self.clone().canonical();
        item.var_terms.is_empty()
            && item.ose_terms.is_empty()
            && item.external_terms.is_empty()
            && item.constant.is_zero()
            && item.uniform_scale_coeff.is_zero()
    }

    fn variable_edges(&self) -> BTreeSet<usize> {
        self.var_terms.iter().map(|(edge_id, _)| *edge_id).collect()
    }

    fn to_surface_expr(&self, edge_exprs: &[LinearEnergyExpr]) -> Result<LinearEnergyExpr> {
        let mut out = LinearEnergyExpr {
            external_terms: self
                .external_terms
                .iter()
                .map(|(edge_id, coeff)| (EdgeIndex(*edge_id), rational_coeff_atom(coeff.clone())))
                .collect(),
            uniform_scale_coeff: rational_coeff_atom(self.uniform_scale_coeff.clone()),
            constant: rational_to_coefficient(self.constant.clone())?,
            ..LinearEnergyExpr::zero()
        }
        .canonical();
        for (edge_id, coeff) in &self.ose_terms {
            out = out
                + LinearEnergyExpr::ose_with_coeff(
                    EdgeIndex(*edge_id),
                    rational_coeff_atom(coeff.clone()),
                );
        }
        for (edge_id, coeff) in &self.var_terms {
            out = out + edge_exprs[*edge_id].clone().scale_rational(coeff.clone());
        }
        Ok(out.canonical())
    }

    fn product_is_zero(factors: &[Self]) -> bool {
        factors.iter().any(Self::is_zero)
    }

    fn degrees(n_edges: usize, factors: &[Self]) -> Vec<usize> {
        let mut degrees = vec![0; n_edges];
        for factor in factors {
            for edge_id in factor.variable_edges() {
                degrees[edge_id] += 1;
            }
        }
        degrees
    }

    fn collect_terms(terms: Vec<(usize, Rational)>) -> Vec<(usize, Rational)> {
        let mut collected = BTreeMap::<usize, Rational>::new();
        for (edge_id, coeff) in terms {
            let entry = collected.entry(edge_id).or_insert_with(Rational::zero);
            *entry = entry.clone() + coeff;
        }
        collected
            .into_iter()
            .filter(|(_, coeff)| !coeff.is_zero())
            .collect()
    }
}

impl Add for KnownLinearExpr {
    type Output = Self;

    fn add(mut self, rhs: Self) -> Self::Output {
        self.var_terms.extend(rhs.var_terms);
        self.ose_terms.extend(rhs.ose_terms);
        self.external_terms.extend(rhs.external_terms);
        self.constant += rhs.constant;
        self.uniform_scale_coeff += rhs.uniform_scale_coeff;
        self.canonical()
    }
}

impl Neg for KnownLinearExpr {
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        for (_, coeff) in &mut self.var_terms {
            *coeff = -coeff.clone();
        }
        for (_, coeff) in &mut self.ose_terms {
            *coeff = -coeff.clone();
        }
        for (_, coeff) in &mut self.external_terms {
            *coeff = -coeff.clone();
        }
        self.constant = -self.constant;
        self.uniform_scale_coeff = -self.uniform_scale_coeff;
        self.canonical()
    }
}

impl Sub for KnownLinearExpr {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

#[derive(Debug, Clone)]
struct AccumulatedVariant {
    prefactor: symbolica::atom::Atom,
    half_edges: Vec<usize>,
    chains: Vec<Vec<HybridSurfaceID>>,
}

#[derive(Debug, Default)]
struct VariantAccumulator {
    variants: Vec<AccumulatedVariant>,
}

impl VariantAccumulator {
    fn add(
        &mut self,
        prefactor: symbolica::atom::Atom,
        half_edges: Vec<usize>,
        chain: Vec<HybridSurfaceID>,
    ) {
        if let Some(existing) = self
            .variants
            .iter_mut()
            .find(|variant| variant.prefactor == prefactor && variant.half_edges == half_edges)
        {
            existing.chains.push(chain);
            return;
        }
        self.variants.push(AccumulatedVariant {
            prefactor,
            half_edges,
            chains: vec![chain],
        });
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum CoeffSymbol {
    Constant,
    Internal(usize),
    External(usize),
    UniformScale,
}

type CoeffMap = BTreeMap<CoeffSymbol, Rational>;
type RationalMatrix = Vec<Vec<Rational>>;

struct LinearSurfaceInterner<'a> {
    expression: &'a mut ThreeDExpression<OrientationID>,
    index: HashMap<(LinearSurfaceKind, LinearEnergyExpr), HybridSurfaceID>,
}

impl<'a> LinearSurfaceInterner<'a> {
    fn new(expression: &'a mut ThreeDExpression<OrientationID>) -> Self {
        Self {
            expression,
            index: HashMap::new(),
        }
    }

    fn intern(&mut self, surface_expr: LinearEnergyExpr, numerator_only: bool) -> HybridSurfaceID {
        self.intern_with_origin(surface_expr, SurfaceOrigin::Physical, numerator_only)
    }

    fn intern_with_origin(
        &mut self,
        surface_expr: LinearEnergyExpr,
        origin: SurfaceOrigin,
        numerator_only: bool,
    ) -> HybridSurfaceID {
        let surface_expr = surface_expr.canonical();
        let kind = classify_surface_kind(&surface_expr);
        let key = (kind, surface_expr.clone());
        if let Some(surface_id) = self.index.get(&key) {
            if !numerator_only && let HybridSurfaceID::Linear(id) = *surface_id {
                self.expression.surfaces.linear_surface_cache[id].numerator_only = false;
                self.expression.surfaces.linear_surface_cache[id].origin = SurfaceOrigin::Physical;
            }
            return *surface_id;
        }

        let id = LinearSurfaceID(self.expression.surfaces.linear_surface_cache.len());
        let surface_id = HybridSurfaceID::Linear(id);
        self.expression
            .surfaces
            .linear_surface_cache
            .push(LinearSurface {
                kind,
                expression: surface_expr,
                origin,
                numerator_only,
            });
        self.index.insert(key, surface_id);
        surface_id
    }
}

struct EnergySolver<'a> {
    signatures: &'a [MomentumSignature],
}

impl<'a> EnergySolver<'a> {
    fn new(signatures: &'a [MomentumSignature]) -> Self {
        Self { signatures }
    }

    fn solve_from_target_edges(
        &self,
        basis: &[usize],
        target_edge_exprs: &[LinearEnergyExpr],
    ) -> Result<Vec<LinearEnergyExpr>> {
        solve_loop_energy_from_target_edge_exprs(self.signatures, basis, target_edge_exprs)
    }

    fn edge_q0_from_loop_exprs(&self, loop_exprs: &[LinearEnergyExpr]) -> Vec<LinearEnergyExpr> {
        edge_q0_from_loop_exprs(self.signatures, loop_exprs)
    }

    fn derivative_matrices(&self, basis: &[usize]) -> Result<(RationalMatrix, RationalMatrix)> {
        let n_loops = self
            .signatures
            .first()
            .map(|signature| signature.loop_signature.len())
            .unwrap_or(0);
        let n_basis = basis.len();
        let matrix = basis
            .iter()
            .map(|edge| {
                self.signatures[*edge]
                    .loop_signature
                    .iter()
                    .map(|value| Rational::from(*value))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let mut loop_by_basis = vec![vec![Rational::zero(); n_basis]; n_loops];

        for (basis_position, _) in basis.iter().enumerate() {
            let rhs = (0..n_basis)
                .map(|row| {
                    if row == basis_position {
                        Rational::one()
                    } else {
                        Rational::zero()
                    }
                })
                .collect::<Vec<_>>();
            let solution =
                solve_rational_system(matrix.clone(), rhs).ok_or(GenerationError::SingularBasis)?;
            for (loop_id, value) in solution.into_iter().enumerate() {
                loop_by_basis[loop_id][basis_position] = value;
            }
        }

        let edge_by_basis = self
            .signatures
            .iter()
            .map(|signature| {
                (0..n_basis)
                    .map(|basis_position| {
                        signature.loop_signature.iter().enumerate().fold(
                            Rational::zero(),
                            |acc, (loop_id, coeff)| {
                                acc + Rational::from(*coeff)
                                    * loop_by_basis[loop_id][basis_position].clone()
                            },
                        )
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        Ok((loop_by_basis, edge_by_basis))
    }
}

struct DenominatorDerivativeExpander<'a> {
    factors: &'a [DenominatorFactor],
}

impl<'a> DenominatorDerivativeExpander<'a> {
    fn new(factors: &'a [DenominatorFactor]) -> Self {
        Self { factors }
    }

    fn terms(&self, gamma: &[usize]) -> Vec<DenominatorTerm> {
        if self.factors.is_empty() {
            return if gamma.iter().all(|order| *order == 0) {
                vec![DenominatorTerm {
                    coeff: Rational::one(),
                    half_edges: Vec::new(),
                    chain: Vec::new(),
                }]
            } else {
                Vec::new()
            };
        }

        let total_factor = multi_factorial(gamma);
        let mut accum = BTreeMap::<(Vec<usize>, Vec<HybridSurfaceID>), Rational>::new();
        self.accumulate_terms(
            0,
            gamma.to_vec(),
            Rational::one(),
            Vec::new(),
            Vec::new(),
            Rational::one(),
            total_factor,
            &mut accum,
        );
        accum
            .into_iter()
            .filter_map(|((half_edges, chain), coeff)| {
                (!coeff.is_zero()).then_some(DenominatorTerm {
                    coeff,
                    half_edges,
                    chain,
                })
            })
            .collect()
    }

    #[allow(clippy::too_many_arguments)]
    fn accumulate_terms(
        &self,
        factor_index: usize,
        remaining: Vec<usize>,
        coeff: Rational,
        half_edges: Vec<usize>,
        chain: Vec<HybridSurfaceID>,
        denom_factorial: Rational,
        total_factor: Rational,
        accum: &mut BTreeMap<(Vec<usize>, Vec<HybridSurfaceID>), Rational>,
    ) {
        if factor_index == self.factors.len() {
            if remaining.iter().any(|order| *order != 0) {
                return;
            }
            let final_coeff = coeff * total_factor / denom_factorial;
            if final_coeff.is_zero() {
                return;
            }
            let mut half_edges = half_edges;
            half_edges.sort_unstable();
            let entry = accum.entry((half_edges, chain)).or_insert(Rational::zero());
            *entry = entry.clone() + final_coeff;
            return;
        }

        for delta in multiindices_leq(&remaining) {
            let Some(piece) = self.factors[factor_index].derivative_piece(&delta) else {
                continue;
            };
            let next_remaining = remaining
                .iter()
                .zip(&delta)
                .map(|(a, b)| a - b)
                .collect::<Vec<_>>();
            let mut next_half_edges = half_edges.clone();
            next_half_edges.extend(piece.half_edges);
            let mut next_chain = chain.clone();
            next_chain.extend(piece.chain);
            self.accumulate_terms(
                factor_index + 1,
                next_remaining,
                coeff.clone() * piece.coeff,
                next_half_edges,
                next_chain,
                denom_factorial.clone() * multi_factorial(&delta),
                total_factor.clone(),
                accum,
            );
        }
    }
}

struct DenominatorFactorDerivativePiece {
    coeff: Rational,
    half_edges: Vec<usize>,
    chain: Vec<HybridSurfaceID>,
}

impl DenominatorFactor {
    fn derivative_piece(&self, delta: &[usize]) -> Option<DenominatorFactorDerivativePiece> {
        let order = delta.iter().sum::<usize>();
        let (power, derivs) = match self {
            Self::Cut { power, derivs, .. } | Self::Surface { power, derivs, .. } => {
                (*power, derivs)
            }
        };
        let mut coeff = if order == 0 {
            Rational::one()
        } else {
            let sign = if order % 2 == 0 { 1 } else { -1 };
            let coeff = rising(power, order);
            if sign < 0 { -coeff } else { coeff }
        };
        for (derivative_order, derivative_coeff) in delta.iter().zip(derivs) {
            if *derivative_order != 0 && derivative_coeff.is_zero() {
                return None;
            }
            coeff *= derivative_coeff.pow_usize(*derivative_order);
        }
        let effective_power = power + order;
        match self {
            Self::Cut { edge, .. } => Some(DenominatorFactorDerivativePiece {
                coeff,
                half_edges: vec![*edge; effective_power],
                chain: Vec::new(),
            }),
            Self::Surface { surface, .. } => Some(DenominatorFactorDerivativePiece {
                coeff,
                half_edges: Vec::new(),
                chain: vec![*surface; effective_power],
            }),
        }
    }
}

struct BoundedCffBuilder<'a> {
    parsed: &'a ParsedGraph,
    bounds: Vec<usize>,
    sampling_scale_mode: NumeratorSamplingScaleMode,
    expression: ThreeDExpression<OrientationID>,
    surface_index: HashMap<(LinearSurfaceKind, LinearEnergyExpr), HybridSurfaceID>,
}

impl<'a> BoundedCffBuilder<'a> {
    fn new(parsed: &'a ParsedGraph, options: &'a Generate3DExpressionOptions) -> Result<Self> {
        let bounds = normalize_energy_degree_bounds(
            &options.energy_degree_bounds,
            parsed.internal_edges.len(),
        )?;
        Ok(Self {
            parsed,
            bounds,
            sampling_scale_mode: options.numerator_sampling_scale,
            expression: ThreeDExpression::new_empty(),
            surface_index: HashMap::new(),
        })
    }

    fn for_bounds(parsed: &'a ParsedGraph, bounds: Vec<usize>) -> Self {
        Self {
            parsed,
            bounds,
            sampling_scale_mode: NumeratorSamplingScaleMode::None,
            expression: ThreeDExpression::new_empty(),
            surface_index: HashMap::new(),
        }
    }

    fn build(mut self) -> Result<ThreeDExpression<OrientationID>> {
        if self.bounds.iter().all(|degree| *degree <= 1) {
            return generate_pure_cff_expression_from_parsed(self.parsed);
        }
        let uniform_sampling_for_nonlinear_degree = self
            .bounds
            .iter()
            .any(|degree| *degree > 1 && self.sampling_scale_mode.is_active_for_degree(*degree));
        if !uniform_sampling_for_nonlinear_degree && self.supports_quadratic_e_surface_only() {
            self.build_quadratic_e_surface_only()?;
            self.finalize_numerator_map_labels();
            return Ok(self.expression);
        }
        if !uniform_sampling_for_nonlinear_degree && self.supports_quadratic_recursive() {
            return self.build_quadratic_recursive(false);
        }
        if self.supports_known_factor_recursive() {
            return KnownFactorCffBuilder::new(self.parsed, self.bounds, self.sampling_scale_mode)
                .build();
        }
        Err(GenerationError::CffHigherEnergyPowerNotImplemented)
    }

    fn supports_quadratic_e_surface_only(&self) -> bool {
        self.parsed.loop_names.len() == 1
            && !self.has_duplicate_signature_ignoring_mass()
            && self.bounds.iter().all(|degree| *degree <= 2)
    }

    fn supports_quadratic_recursive(&self) -> bool {
        self.bounds.iter().all(|degree| *degree <= 2)
            && (self.parsed.loop_names.len() == 1
                || self.bounds.iter().filter(|degree| **degree > 1).count() == 1)
    }

    fn supports_known_factor_recursive(&self) -> bool {
        self.bounds.iter().any(|degree| *degree > 1)
            || RepeatedLtdBuilder::logical_channels(self.parsed)
                .iter()
                .any(|channel| {
                    channel
                        .members
                        .iter()
                        .map(|edge_id| self.bounds[*edge_id])
                        .sum::<usize>()
                        > 2
                })
    }

    fn has_duplicate_signature_ignoring_mass(&self) -> bool {
        let mut counts = BTreeMap::<MomentumSignature, usize>::new();
        for edge in &self.parsed.internal_edges {
            let (signature, _) = edge.signature.canonical_up_to_sign();
            *counts.entry(signature).or_default() += 1;
        }
        counts.values().any(|count| *count > 1)
    }

    fn build_quadratic_e_surface_only(&mut self) -> Result<()> {
        let quadratic_edges = self
            .bounds
            .iter()
            .enumerate()
            .filter_map(|(edge_id, degree)| (*degree > 1).then_some(edge_id))
            .collect::<Vec<_>>();

        for size in 0..=quadratic_edges.len() {
            for pinched in quadratic_edges.iter().copied().combinations(size) {
                self.lift_contracted_cff_terms(&pinched)?;
            }
        }
        Ok(())
    }

    fn build_quadratic_recursive(
        mut self,
        lower_sector_base: bool,
    ) -> Result<ThreeDExpression<OrientationID>> {
        let Some(active_edge) = self.bounds.iter().position(|degree| *degree > 1) else {
            return if lower_sector_base {
                self.lower_sector_base_expression()
            } else {
                generate_pure_cff_expression_from_parsed(self.parsed)
            };
        };

        let mut remainder_bounds = self.bounds.clone();
        remainder_bounds[active_edge] = 1;
        let remainder_expression = BoundedCffBuilder::for_bounds(self.parsed, remainder_bounds)
            .build_quadratic_recursive(lower_sector_base)?;
        self.append_recursive_remainder_terms(active_edge, &remainder_expression)?;

        let (subparsed, sub_to_orig) = self.delete_parsed_edges(&[active_edge]);
        if !subparsed.internal_edges.is_empty() {
            let sub_bounds = sub_to_orig
                .iter()
                .map(|orig_id| self.bounds[*orig_id])
                .collect::<Vec<_>>();
            let contact_expression = BoundedCffBuilder::for_bounds(&subparsed, sub_bounds)
                .build_quadratic_recursive(true)?;
            self.append_recursive_contact_terms(active_edge, &contact_expression, &sub_to_orig)?;
        }

        self.finalize_numerator_map_labels();
        Ok(self.expression)
    }

    fn lower_sector_base_expression(&self) -> Result<ThreeDExpression<OrientationID>> {
        LowerSectorCffBuilder::new(self.parsed).build()
    }

    fn append_recursive_remainder_terms(
        &mut self,
        edge_id: usize,
        source: &ThreeDExpression<OrientationID>,
    ) -> Result<()> {
        let edge_map = (0..self.parsed.internal_edges.len())
            .map(|id| (id, id))
            .collect::<BTreeMap<_, _>>();
        let surface_map = self.copy_expression_surfaces(source, &edge_map);

        for orientation in &source.orientations {
            for variant in &orientation.variants {
                let denominator = variant
                    .denominator
                    .clone()
                    .map(|surface_id| map_surface_id(surface_id, &surface_map));
                let base_half_edges = variant
                    .half_edges
                    .iter()
                    .map(|edge| edge.0)
                    .collect::<Vec<_>>();
                let base_num_surfaces = variant
                    .numerator_surfaces
                    .iter()
                    .map(|surface_id| map_surface_id(*surface_id, &surface_map))
                    .collect::<Vec<_>>();

                for component in self.finite_pole_remainder_components(
                    edge_id,
                    orientation.edge_energy_map[edge_id].clone(),
                ) {
                    let mut edge_exprs = orientation.edge_energy_map.clone();
                    edge_exprs[edge_id] =
                        LinearEnergyExpr::ose(EdgeIndex(edge_id), i64::from(component.sample));
                    let mut half_edges = base_half_edges.clone();
                    half_edges.extend(component.half_edges);
                    half_edges.sort_unstable();
                    let mut numerator_surfaces = base_num_surfaces.clone();
                    numerator_surfaces.extend(component.numerator_surfaces);
                    let prefactor =
                        rational_from_coefficient(&variant.prefactor) * component.prefactor.clone();
                    if prefactor.is_zero() {
                        continue;
                    }
                    self.push_variant_for_maps(
                        orientation.loop_energy_map.clone(),
                        edge_exprs,
                        crate::expression::CFFVariant {
                            origin: Some(format!(
                                "bounded_degree_quadratic_recursive_remainder:e{edge_id}={}",
                                if component.sample > 0 { "+" } else { "-" }
                            )),
                            prefactor: rational_to_coefficient(prefactor)?,
                            half_edges: half_edges.into_iter().map(EdgeIndex).collect(),
                            uniform_scale_power: variant.uniform_scale_power,
                            numerator_surfaces,
                            denominator: denominator.clone(),
                        },
                    );
                }
            }
        }
        Ok(())
    }

    fn append_recursive_contact_terms(
        &mut self,
        edge_id: usize,
        source: &ThreeDExpression<OrientationID>,
        sub_to_orig: &[usize],
    ) -> Result<()> {
        let edge_map = sub_to_orig
            .iter()
            .enumerate()
            .map(|(sub_id, orig_id)| (sub_id, *orig_id))
            .collect::<BTreeMap<_, _>>();
        let surface_map = self.copy_expression_surfaces(source, &edge_map);
        let signatures = self
            .parsed
            .internal_edges
            .iter()
            .map(|edge| edge.signature.clone())
            .collect::<Vec<_>>();

        for orientation in &source.orientations {
            let full_loop_exprs = orientation
                .loop_energy_map
                .iter()
                .cloned()
                .map(|expr| expr.remap_internal_edges(&edge_map))
                .collect::<Vec<_>>();
            let current_edge_exprs = edge_q0_from_loop_exprs(&signatures, &full_loop_exprs);
            let mut base_edge_exprs =
                vec![LinearEnergyExpr::zero(); self.parsed.internal_edges.len()];
            for (sub_id, orig_id) in sub_to_orig.iter().enumerate() {
                base_edge_exprs[*orig_id] = orientation.edge_energy_map[sub_id]
                    .clone()
                    .remap_internal_edges(&edge_map);
            }

            let components = self.finite_pole_contact_components(
                edge_id,
                self.bounds[edge_id],
                current_edge_exprs[edge_id].clone(),
            );
            for variant in &orientation.variants {
                let denominator = variant
                    .denominator
                    .clone()
                    .map(|surface_id| map_surface_id(surface_id, &surface_map));
                let base_half_edges = variant
                    .half_edges
                    .iter()
                    .map(|edge| edge_map.get(&edge.0).copied().unwrap_or(edge.0))
                    .collect::<Vec<_>>();
                let base_num_surfaces = variant
                    .numerator_surfaces
                    .iter()
                    .map(|surface_id| map_surface_id(*surface_id, &surface_map))
                    .collect::<Vec<_>>();

                for component in &components {
                    let mut edge_exprs = base_edge_exprs.clone();
                    edge_exprs[edge_id] = if component.sample == 0 {
                        LinearEnergyExpr::zero()
                    } else {
                        LinearEnergyExpr::ose(EdgeIndex(edge_id), i64::from(component.sample))
                    };
                    let mut half_edges = base_half_edges.clone();
                    half_edges.extend(component.half_edges.iter().copied());
                    half_edges.sort_unstable();
                    let mut numerator_surfaces = base_num_surfaces.clone();
                    numerator_surfaces.extend(component.numerator_surfaces.iter().copied());
                    let prefactor =
                        rational_from_coefficient(&variant.prefactor) * component.prefactor.clone();
                    if prefactor.is_zero() {
                        continue;
                    }
                    self.push_variant_for_maps(
                        full_loop_exprs.clone(),
                        edge_exprs,
                        crate::expression::CFFVariant {
                            origin: Some(format!(
                                "bounded_degree_quadratic_recursive_contact:e{edge_id}={}",
                                match component.sample {
                                    0 => "0",
                                    value if value > 0 => "+",
                                    _ => "-",
                                }
                            )),
                            prefactor: rational_to_coefficient(prefactor)?,
                            half_edges: half_edges.into_iter().map(EdgeIndex).collect(),
                            uniform_scale_power: variant.uniform_scale_power,
                            numerator_surfaces,
                            denominator: denominator.clone(),
                        },
                    );
                }
            }
        }
        Ok(())
    }

    fn lift_contracted_cff_terms(&mut self, pinched_edges: &[usize]) -> Result<()> {
        let (subparsed, sub_to_orig) = self.contract_parsed_edges(pinched_edges);
        if subparsed.internal_edges.is_empty() {
            return Ok(());
        }

        let use_terminal_ltd = subparsed.internal_edges.len() == self.parsed.loop_names.len();
        let sub_expression = if use_terminal_ltd {
            generate_pure_ltd_expression_from_parsed(&subparsed)?
        } else {
            generate_pure_cff_expression_from_parsed(&subparsed)?
        };
        let edge_map = sub_to_orig
            .iter()
            .enumerate()
            .map(|(sub_id, orig_id)| (sub_id, *orig_id))
            .collect::<BTreeMap<_, _>>();
        let surface_map = self.copy_expression_surfaces(&sub_expression, &edge_map);
        let signatures = self
            .parsed
            .internal_edges
            .iter()
            .map(|edge| edge.signature.clone())
            .collect::<Vec<_>>();

        for orientation in &sub_expression.orientations {
            let full_loop_exprs = orientation
                .loop_energy_map
                .iter()
                .cloned()
                .map(|expr| expr.remap_internal_edges(&edge_map))
                .collect::<Vec<_>>();
            let current_edge_exprs = edge_q0_from_loop_exprs(&signatures, &full_loop_exprs);
            let mut base_edge_exprs =
                vec![LinearEnergyExpr::zero(); self.parsed.internal_edges.len()];
            for (sub_id, orig_id) in sub_to_orig.iter().enumerate() {
                base_edge_exprs[*orig_id] = orientation.edge_energy_map[sub_id]
                    .clone()
                    .remap_internal_edges(&edge_map);
            }

            let per_edge_components = pinched_edges
                .iter()
                .map(|edge_id| {
                    self.finite_pole_contact_components(
                        *edge_id,
                        self.bounds[*edge_id],
                        current_edge_exprs[*edge_id].clone(),
                    )
                })
                .collect::<Vec<_>>();
            if per_edge_components.iter().any(Vec::is_empty) {
                continue;
            }
            let component_choices = if per_edge_components.is_empty() {
                vec![Vec::new()]
            } else {
                per_edge_components
                    .into_iter()
                    .multi_cartesian_product()
                    .collect::<Vec<_>>()
            };

            for variant in &orientation.variants {
                let remapped_denominator = variant
                    .denominator
                    .clone()
                    .map(|surface_id| map_surface_id(surface_id, &surface_map));
                let remapped_half_edges = variant
                    .half_edges
                    .iter()
                    .map(|edge_id| edge_map.get(&edge_id.0).copied().unwrap_or(edge_id.0))
                    .collect::<Vec<_>>();
                let remapped_numerator_surfaces = variant
                    .numerator_surfaces
                    .iter()
                    .map(|surface_id| map_surface_id(*surface_id, &surface_map))
                    .collect::<Vec<_>>();
                for components in &component_choices {
                    let mut prefactor = rational_from_coefficient(&variant.prefactor);
                    let mut half_edges = remapped_half_edges.clone();
                    let mut numerator_surfaces = remapped_numerator_surfaces.clone();
                    let mut edge_exprs = base_edge_exprs.clone();
                    let mut sample_labels = Vec::new();

                    for (edge_id, component) in pinched_edges.iter().zip(components) {
                        prefactor *= component.prefactor.clone();
                        half_edges.extend(component.half_edges.iter().copied());
                        numerator_surfaces.extend(component.numerator_surfaces.iter().copied());
                        if component.sample == 0 {
                            edge_exprs[*edge_id] = LinearEnergyExpr::zero();
                            sample_labels.push(format!("e{edge_id}=0"));
                        } else {
                            edge_exprs[*edge_id] = LinearEnergyExpr::ose(
                                EdgeIndex(*edge_id),
                                i64::from(component.sample),
                            );
                            sample_labels.push(format!(
                                "e{edge_id}={}",
                                if component.sample > 0 { "+" } else { "-" }
                            ));
                        }
                    }
                    if !pinched_edges.is_empty()
                        && !use_terminal_ltd
                        && pinched_edges.len() % 2 == 1
                    {
                        prefactor = -prefactor;
                    }
                    if prefactor.is_zero() {
                        continue;
                    }
                    half_edges.sort_unstable();

                    let origin = if sample_labels.is_empty() {
                        "bounded_degree_e_surface_pinch_cff".to_string()
                    } else {
                        format!(
                            "bounded_degree_e_surface_pinch_cff:{}",
                            sample_labels.join("|")
                        )
                    };
                    self.push_variant_for_maps(
                        full_loop_exprs.clone(),
                        edge_exprs,
                        crate::expression::CFFVariant {
                            origin: Some(match &variant.origin {
                                Some(source) => format!("{origin}:{source}"),
                                None => origin,
                            }),
                            prefactor: rational_to_coefficient(prefactor)?,
                            half_edges: half_edges.into_iter().map(EdgeIndex).collect(),
                            uniform_scale_power: variant.uniform_scale_power,
                            numerator_surfaces,
                            denominator: remapped_denominator.clone(),
                        },
                    );
                }
            }
        }

        Ok(())
    }

    fn finite_pole_contact_components(
        &mut self,
        edge_id: usize,
        bound: usize,
        current_expr: LinearEnergyExpr,
    ) -> Vec<ContactComponent> {
        let q_surface = (!current_expr.is_zero())
            .then(|| self.intern_surface(current_expr, SurfaceOrigin::Helper, true));
        let mut components = Vec::new();
        for (sample, poly) in contact_weight_polys(bound) {
            for (power, coeff) in poly.into_iter().enumerate() {
                if coeff.is_zero() || (power > 0 && q_surface.is_none()) {
                    continue;
                }
                let prefactor = coeff * rational_pow_i64(2, power + 2);
                if prefactor.is_zero() {
                    continue;
                }
                components.push(ContactComponent {
                    sample,
                    prefactor,
                    half_edges: std::iter::repeat_n(edge_id, power + 2).collect(),
                    numerator_surfaces: q_surface
                        .into_iter()
                        .cycle()
                        .take(power)
                        .collect::<Vec<_>>(),
                });
            }
        }
        components
    }

    fn finite_pole_remainder_components(
        &mut self,
        edge_id: usize,
        current_expr: LinearEnergyExpr,
    ) -> Vec<ContactComponent> {
        let mut components = Vec::new();
        let plus_expr = current_expr.clone() + LinearEnergyExpr::ose(EdgeIndex(edge_id), 1);
        if !plus_expr.is_zero() {
            let plus = self.intern_surface(plus_expr, SurfaceOrigin::Helper, true);
            components.push(ContactComponent {
                sample: 1,
                prefactor: Rational::one(),
                half_edges: vec![edge_id],
                numerator_surfaces: vec![plus],
            });
        }

        let minus_expr = LinearEnergyExpr::ose(EdgeIndex(edge_id), 1) - current_expr;
        if !minus_expr.is_zero() {
            let minus = self.intern_surface(minus_expr, SurfaceOrigin::Helper, true);
            components.push(ContactComponent {
                sample: -1,
                prefactor: Rational::one(),
                half_edges: vec![edge_id],
                numerator_surfaces: vec![minus],
            });
        }
        components
    }

    fn copy_expression_surfaces(
        &mut self,
        source: &ThreeDExpression<OrientationID>,
        edge_map: &BTreeMap<usize, usize>,
    ) -> HashMap<HybridSurfaceID, HybridSurfaceID> {
        source
            .surfaces
            .linear_surface_cache
            .iter_enumerated()
            .map(|(id, surface)| {
                (
                    HybridSurfaceID::Linear(id),
                    self.intern_surface(
                        surface.expression.clone().remap_internal_edges(edge_map),
                        surface.origin,
                        surface.numerator_only,
                    ),
                )
            })
            .collect()
    }

    fn intern_surface(
        &mut self,
        surface_expr: LinearEnergyExpr,
        origin: SurfaceOrigin,
        numerator_only: bool,
    ) -> HybridSurfaceID {
        let surface_expr = surface_expr.canonical();
        let kind = classify_surface_kind(&surface_expr);
        let key = (kind, surface_expr.clone());
        if let Some(surface_id) = self.surface_index.get(&key) {
            if !numerator_only && let HybridSurfaceID::Linear(id) = *surface_id {
                self.expression.surfaces.linear_surface_cache[id].origin = SurfaceOrigin::Physical;
                self.expression.surfaces.linear_surface_cache[id].numerator_only = false;
            }
            return *surface_id;
        }

        let id = LinearSurfaceID(self.expression.surfaces.linear_surface_cache.len());
        let surface_id = HybridSurfaceID::Linear(id);
        self.expression
            .surfaces
            .linear_surface_cache
            .push(LinearSurface {
                kind,
                expression: surface_expr,
                origin,
                numerator_only,
            });
        self.surface_index.insert(key, surface_id);
        surface_id
    }

    fn push_variant_for_maps(
        &mut self,
        loop_energy_map: Vec<LinearEnergyExpr>,
        edge_energy_map: Vec<LinearEnergyExpr>,
        variant: crate::expression::CFFVariant,
    ) {
        let (orientation, base_label) = orientation_from_edge_exprs(&edge_energy_map);
        if let Some(existing) = self
            .expression
            .orientations
            .iter_mut()
            .find(|orientation_expr| {
                orientation_expr.data.label.as_deref() == Some(base_label.as_str())
                    && orientation_expr.loop_energy_map == loop_energy_map
                    && orientation_expr.edge_energy_map == edge_energy_map
            })
        {
            existing.variants.push(variant);
            return;
        }
        let mut data = OrientationData::new(orientation);
        data.label = Some(base_label);
        self.expression.orientations.push(OrientationExpression {
            data,
            loop_energy_map,
            edge_energy_map,
            variants: vec![variant],
        });
    }

    fn finalize_numerator_map_labels(&mut self) {
        assign_numerator_map_labels(&mut self.expression.orientations);
    }

    fn contract_parsed_edges(&self, pinched_edges: &[usize]) -> (ParsedGraph, Vec<usize>) {
        let pinched = pinched_edges.iter().copied().collect::<BTreeSet<_>>();
        let mut parent = self
            .parsed
            .node_name_to_internal
            .values()
            .copied()
            .map(|node| (node, node))
            .collect::<BTreeMap<_, _>>();

        for edge in &self.parsed.internal_edges {
            if pinched.contains(&edge.edge_id) {
                union_nodes(&mut parent, edge.tail, edge.head);
            }
        }

        let parent_keys = parent.keys().copied().collect::<Vec<_>>();
        let roots = parent_keys
            .iter()
            .map(|node| find_node_root(&mut parent.clone(), *node))
            .collect::<BTreeSet<_>>();
        let root_to_new = roots
            .into_iter()
            .enumerate()
            .map(|(new_id, root)| (root, new_id))
            .collect::<BTreeMap<_, _>>();
        let old_to_new = parent_keys
            .iter()
            .map(|node| {
                let root = find_node_root(&mut parent, *node);
                (*node, root_to_new[&root])
            })
            .collect::<BTreeMap<_, _>>();

        let mut sub_to_orig = Vec::new();
        let mut internal_edges = Vec::new();
        for edge in &self.parsed.internal_edges {
            if pinched.contains(&edge.edge_id) {
                continue;
            }
            let sub_id = internal_edges.len();
            sub_to_orig.push(edge.edge_id);
            internal_edges.push(ParsedGraphInternalEdge {
                edge_id: sub_id,
                tail: old_to_new[&edge.tail],
                head: old_to_new[&edge.head],
                label: edge.label.clone(),
                mass_key: edge.mass_key.clone(),
                signature: edge.signature.clone(),
                had_pow: edge.had_pow,
            });
        }

        let external_edges = self
            .parsed
            .external_edges
            .iter()
            .map(|edge| ParsedGraphExternalEdge {
                edge_id: edge.edge_id,
                source: edge.source.map(|source| old_to_new[&source]),
                destination: edge.destination.map(|destination| old_to_new[&destination]),
                label: edge.label.clone(),
                external_coefficients: edge.external_coefficients.clone(),
            })
            .collect::<Vec<_>>();
        let node_name_to_internal = old_to_new
            .values()
            .copied()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .map(|node| (format!("c{node}"), node))
            .collect();

        (
            ParsedGraph {
                internal_edges,
                external_edges,
                initial_state_cut_edges: remap_initial_state_cut_edges(self.parsed, &sub_to_orig),
                loop_names: self.parsed.loop_names.clone(),
                external_names: self.parsed.external_names.clone(),
                node_name_to_internal,
            },
            sub_to_orig,
        )
    }

    fn delete_parsed_edges(&self, deleted_edges: &[usize]) -> (ParsedGraph, Vec<usize>) {
        let deleted = deleted_edges.iter().copied().collect::<BTreeSet<_>>();
        let mut sub_to_orig = Vec::new();
        let mut internal_edges = Vec::new();
        for edge in &self.parsed.internal_edges {
            if deleted.contains(&edge.edge_id) {
                continue;
            }
            let sub_id = internal_edges.len();
            sub_to_orig.push(edge.edge_id);
            internal_edges.push(ParsedGraphInternalEdge {
                edge_id: sub_id,
                tail: edge.tail,
                head: edge.head,
                label: edge.label.clone(),
                mass_key: edge.mass_key.clone(),
                signature: edge.signature.clone(),
                had_pow: edge.had_pow,
            });
        }

        (
            ParsedGraph {
                internal_edges,
                external_edges: self.parsed.external_edges.clone(),
                initial_state_cut_edges: remap_initial_state_cut_edges(self.parsed, &sub_to_orig),
                loop_names: self.parsed.loop_names.clone(),
                external_names: self.parsed.external_names.clone(),
                node_name_to_internal: self.parsed.node_name_to_internal.clone(),
            },
            sub_to_orig,
        )
    }
}

struct KnownFactorCffBuilder<'a> {
    original: &'a ParsedGraph,
    bounds: Vec<usize>,
    sampling_scale_mode: NumeratorSamplingScaleMode,
    expression: ThreeDExpression<OrientationID>,
    surface_index: HashMap<(LinearSurfaceKind, LinearEnergyExpr), HybridSurfaceID>,
}

impl<'a> KnownFactorCffBuilder<'a> {
    fn new(
        original: &'a ParsedGraph,
        bounds: Vec<usize>,
        sampling_scale_mode: NumeratorSamplingScaleMode,
    ) -> Self {
        Self {
            original,
            bounds,
            sampling_scale_mode,
            expression: ThreeDExpression::new_empty(),
            surface_index: HashMap::new(),
        }
    }

    fn build(mut self) -> Result<ThreeDExpression<OrientationID>> {
        let local_to_orig = (0..self.original.internal_edges.len()).collect::<Vec<_>>();
        self.channel_recursive_terms(
            self.original,
            &local_to_orig,
            &BTreeMap::new(),
            &[],
            &[],
            0,
            Rational::one(),
            false,
            0,
        )?;
        self.finalize_numerator_map_labels();
        Ok(self.expression)
    }

    #[allow(clippy::too_many_arguments)]
    fn channel_recursive_terms(
        &mut self,
        parsed: &ParsedGraph,
        local_to_orig: &[usize],
        replacements: &BTreeMap<usize, LinearEnergyExpr>,
        known_factors: &[KnownLinearExpr],
        extra_half_edges: &[usize],
        extra_uniform_scale_power: usize,
        prefactor: Rational,
        lower_sector_base: bool,
        depth: usize,
    ) -> Result<()> {
        if depth > 8 {
            return Err(GenerationError::CffHigherEnergyPowerNotImplemented);
        }
        let Some((channel, degree)) =
            self.active_repeated_channel(parsed, local_to_orig, replacements)
        else {
            return self.recursive_terms(
                parsed,
                local_to_orig,
                replacements,
                known_factors,
                extra_half_edges,
                extra_uniform_scale_power,
                prefactor,
                lower_sector_base,
                0,
            );
        };

        let rep_local = channel.rep_edge;
        let rep_orig = local_to_orig[rep_local];
        let rep_signature = parsed.internal_edges[rep_local].signature.clone();
        let nodes = Self::interpolation_nodes(degree);
        let use_uniform_scale = self.sampling_scale_mode.is_active_for_degree(degree);

        for (node_idx, sample) in nodes.iter().enumerate() {
            let basis_poly = lagrange_basis(&nodes, node_idx);
            let channel_terms = if use_uniform_scale {
                Self::channel_uniform_normal_form_terms(&basis_poly, channel.power)
            } else {
                Self::channel_normal_form_terms(&basis_poly, channel.power)
            };
            for term in channel_terms {
                if term.coeff.is_zero() {
                    continue;
                }
                let keep = channel
                    .members
                    .iter()
                    .take(term.remaining_power)
                    .copied()
                    .collect::<BTreeSet<_>>();
                let delete = channel
                    .members
                    .iter()
                    .copied()
                    .filter(|local_id| !keep.contains(local_id))
                    .collect::<Vec<_>>();
                let (subparsed, sub_to_local) =
                    BoundedCffBuilder::for_bounds(parsed, vec![0; parsed.internal_edges.len()])
                        .delete_parsed_edges(&delete);
                if subparsed.internal_edges.is_empty() {
                    continue;
                }
                let sub_local_to_orig = sub_to_local
                    .iter()
                    .map(|local_id| local_to_orig[*local_id])
                    .collect::<Vec<_>>();

                let mut sub_replacements = replacements.clone();
                for local_id in &channel.members {
                    let orig_id = local_to_orig[*local_id];
                    let rel_sign = Self::relative_signature_sign(
                        &rep_signature,
                        &parsed.internal_edges[*local_id].signature,
                    )?;
                    sub_replacements.insert(
                        orig_id,
                        Self::channel_replacement_expr(
                            orig_id,
                            *sample,
                            rel_sign,
                            use_uniform_scale,
                        ),
                    );
                }

                let mut sub_known = known_factors
                    .iter()
                    .map(|factor| self.remap_known_factor_to_sub(parsed, &subparsed, factor))
                    .collect::<Result<Vec<_>>>()?;
                if term.parity != 0 {
                    sub_known.push(self.known_signature_expr(&subparsed, &rep_signature)?);
                }
                if term.positive_ose_power != 0 {
                    sub_known.extend(
                        (0..term.positive_ose_power).map(|_| KnownLinearExpr::ose(rep_orig, 1)),
                    );
                }
                if term.cancelled_power != 0 {
                    let y_expr = self.known_signature_expr(&subparsed, &rep_signature)?;
                    let plus = (y_expr.clone() + KnownLinearExpr::ose(rep_orig, 1)).canonical();
                    let minus = (y_expr - KnownLinearExpr::ose(rep_orig, 1)).canonical();
                    for _ in 0..term.cancelled_power {
                        sub_known.push(minus.clone());
                        sub_known.push(plus.clone());
                    }
                }

                let mut next_half_edges = extra_half_edges.to_vec();
                let channel_uniform_power = if use_uniform_scale {
                    term.inverse_power
                } else {
                    next_half_edges.extend(std::iter::repeat_n(rep_orig, term.inverse_power));
                    0
                };
                let channel_prefactor = prefactor.clone()
                    * term.coeff
                    * if use_uniform_scale {
                        Rational::one()
                    } else {
                        rational_pow_i64(2, term.inverse_power)
                    };
                if channel_prefactor.is_zero() {
                    continue;
                }

                self.channel_recursive_terms(
                    &subparsed,
                    &sub_local_to_orig,
                    &sub_replacements,
                    &sub_known,
                    &next_half_edges,
                    extra_uniform_scale_power + channel_uniform_power,
                    channel_prefactor,
                    lower_sector_base || !delete.is_empty(),
                    depth + 1,
                )?;
            }
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn recursive_terms(
        &mut self,
        parsed: &ParsedGraph,
        local_to_orig: &[usize],
        replacements: &BTreeMap<usize, LinearEnergyExpr>,
        known_factors: &[KnownLinearExpr],
        extra_half_edges: &[usize],
        extra_uniform_scale_power: usize,
        prefactor: Rational,
        lower_sector_base: bool,
        depth: usize,
    ) -> Result<()> {
        if depth > 12 {
            return Err(GenerationError::CffHigherEnergyPowerNotImplemented);
        }
        let known_factors = known_factors
            .iter()
            .cloned()
            .map(KnownLinearExpr::canonical)
            .collect::<Vec<_>>();
        if KnownLinearExpr::product_is_zero(&known_factors) {
            return Ok(());
        }

        let total_bounds =
            self.known_total_bounds(parsed, local_to_orig, replacements, &known_factors);
        let active = total_bounds
            .iter()
            .enumerate()
            .find_map(|(idx, degree)| {
                let orig = local_to_orig[idx];
                (*degree > 1 && self.bounds[orig] > 2 && !replacements.contains_key(&orig))
                    .then_some(idx)
            })
            .or_else(|| {
                total_bounds
                    .iter()
                    .enumerate()
                    .find_map(|(idx, degree)| (*degree > 1).then_some(idx))
            });
        let Some(active) = active else {
            return self.append_base_terms(
                parsed,
                local_to_orig,
                replacements,
                &known_factors,
                extra_half_edges,
                extra_uniform_scale_power,
                prefactor.clone(),
                lower_sector_base,
            );
        };

        let orig_edge = local_to_orig[active];
        for (sample, weight) in [
            (
                1,
                KnownLinearExpr::var(active, 1) + KnownLinearExpr::ose(orig_edge, 1),
            ),
            (
                -1,
                KnownLinearExpr::ose(orig_edge, 1) - KnownLinearExpr::var(active, 1),
            ),
        ] {
            let mut sampled_factors = known_factors
                .iter()
                .map(|factor| factor.replace_var_with_ose(active, sample, orig_edge))
                .collect::<Vec<_>>();
            sampled_factors.push(weight);
            if KnownLinearExpr::product_is_zero(&sampled_factors) {
                continue;
            }
            let mut next_replacements = replacements.clone();
            next_replacements
                .entry(orig_edge)
                .or_insert_with(|| LinearEnergyExpr::ose(EdgeIndex(orig_edge), i64::from(sample)));
            let mut next_half_edges = extra_half_edges.to_vec();
            next_half_edges.push(orig_edge);
            self.recursive_terms(
                parsed,
                local_to_orig,
                &next_replacements,
                &sampled_factors,
                &next_half_edges,
                extra_uniform_scale_power,
                prefactor.clone(),
                lower_sector_base,
                depth + 1,
            )?;
        }

        let (subparsed, sub_to_local) =
            BoundedCffBuilder::for_bounds(parsed, vec![0; parsed.internal_edges.len()])
                .delete_parsed_edges(&[active]);
        if subparsed.internal_edges.is_empty() {
            return Ok(());
        }
        let sub_local_to_orig = sub_to_local
            .iter()
            .map(|local_id| local_to_orig[*local_id])
            .collect::<Vec<_>>();
        for (sample, poly) in contact_weight_polys(total_bounds[active]) {
            for (power, coeff) in poly.into_iter().enumerate() {
                if coeff.is_zero() {
                    continue;
                }
                let mut sampled_factors = known_factors
                    .iter()
                    .map(|factor| factor.replace_var_with_ose(active, sample, orig_edge))
                    .collect::<Vec<_>>();
                sampled_factors.extend((0..power).map(|_| KnownLinearExpr::var(active, 1)));
                if KnownLinearExpr::product_is_zero(&sampled_factors) {
                    continue;
                }
                let remapped_factors = sampled_factors
                    .iter()
                    .map(|factor| self.remap_known_factor_to_sub(parsed, &subparsed, factor))
                    .collect::<Result<Vec<_>>>()?;
                let mut next_replacements = replacements.clone();
                next_replacements.entry(orig_edge).or_insert_with(|| {
                    if sample == 0 {
                        LinearEnergyExpr::zero()
                    } else {
                        LinearEnergyExpr::ose(EdgeIndex(orig_edge), i64::from(sample))
                    }
                });
                let mut next_half_edges = extra_half_edges.to_vec();
                next_half_edges.extend(std::iter::repeat_n(orig_edge, power + 2));
                let branch_prefactor = prefactor.clone() * coeff * rational_pow_i64(2, power + 2);
                self.recursive_terms(
                    &subparsed,
                    &sub_local_to_orig,
                    &next_replacements,
                    &remapped_factors,
                    &next_half_edges,
                    extra_uniform_scale_power,
                    branch_prefactor,
                    true,
                    depth + 1,
                )?;
            }
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn append_base_terms(
        &mut self,
        parsed: &ParsedGraph,
        local_to_orig: &[usize],
        replacements: &BTreeMap<usize, LinearEnergyExpr>,
        known_factors: &[KnownLinearExpr],
        extra_half_edges: &[usize],
        extra_uniform_scale_power: usize,
        prefactor: Rational,
        lower_sector_base: bool,
    ) -> Result<()> {
        if KnownLinearExpr::product_is_zero(known_factors) {
            return Ok(());
        }
        let base_expression = if lower_sector_base {
            LowerSectorCffBuilder::new(parsed).build()?
        } else {
            generate_pure_cff_expression_from_parsed(parsed)?
        };
        let edge_map = local_to_orig
            .iter()
            .enumerate()
            .map(|(local_id, orig_id)| (local_id, *orig_id))
            .collect::<BTreeMap<_, _>>();
        let surface_map = self.copy_expression_surfaces(&base_expression, &edge_map);
        let n_original = self.original.internal_edges.len();

        for orientation in &base_expression.orientations {
            let local_edge_exprs = orientation
                .edge_energy_map
                .iter()
                .cloned()
                .map(|expr| expr.remap_internal_edges(&edge_map))
                .collect::<Vec<_>>();
            let loop_exprs = orientation
                .loop_energy_map
                .iter()
                .cloned()
                .map(|expr| expr.remap_internal_edges(&edge_map))
                .collect::<Vec<_>>();
            let mut full_edge_exprs = vec![LinearEnergyExpr::zero(); n_original];
            for (local_id, orig_id) in local_to_orig.iter().enumerate() {
                full_edge_exprs[*orig_id] = local_edge_exprs[local_id].clone();
            }
            for (orig_id, expr) in replacements {
                full_edge_exprs[*orig_id] = expr.clone();
            }

            for variant in &orientation.variants {
                let coeff = rational_from_coefficient(&variant.prefactor) * prefactor.clone();
                if coeff.is_zero() {
                    continue;
                }
                let mut numerator_surfaces = variant
                    .numerator_surfaces
                    .iter()
                    .map(|surface_id| map_surface_id(*surface_id, &surface_map))
                    .collect::<Vec<_>>();
                let mut skip = false;
                for factor in known_factors {
                    let surface_expr = factor.to_surface_expr(&local_edge_exprs)?;
                    if surface_expr.is_zero() {
                        skip = true;
                        break;
                    }
                    if surface_expr.is_one() {
                        continue;
                    }
                    numerator_surfaces.push(self.intern_surface(
                        surface_expr,
                        SurfaceOrigin::Helper,
                        true,
                    ));
                }
                if skip {
                    continue;
                }
                let mut half_edges = variant
                    .half_edges
                    .iter()
                    .map(|edge| local_to_orig[edge.0])
                    .chain(extra_half_edges.iter().copied())
                    .collect::<Vec<_>>();
                half_edges.sort_unstable();
                let denominator = variant
                    .denominator
                    .clone()
                    .map(|surface_id| map_surface_id(surface_id, &surface_map));
                self.push_variant_for_maps(
                    loop_exprs.clone(),
                    full_edge_exprs.clone(),
                    crate::expression::CFFVariant {
                        origin: Some("bounded_degree_known_factor_cff".to_string()),
                        prefactor: rational_to_coefficient(coeff)?,
                        half_edges: half_edges.into_iter().map(EdgeIndex).collect(),
                        uniform_scale_power: variant.uniform_scale_power
                            + extra_uniform_scale_power,
                        numerator_surfaces,
                        denominator,
                    },
                );
            }
        }
        Ok(())
    }

    fn known_total_bounds(
        &self,
        parsed: &ParsedGraph,
        local_to_orig: &[usize],
        replacements: &BTreeMap<usize, LinearEnergyExpr>,
        known_factors: &[KnownLinearExpr],
    ) -> Vec<usize> {
        let known = KnownLinearExpr::degrees(parsed.internal_edges.len(), known_factors);
        local_to_orig
            .iter()
            .enumerate()
            .map(|(local_id, orig_id)| {
                let blackbox = if replacements.contains_key(orig_id) {
                    0
                } else {
                    self.bounds[*orig_id]
                };
                blackbox + known[local_id]
            })
            .collect()
    }

    fn active_repeated_channel(
        &self,
        parsed: &ParsedGraph,
        local_to_orig: &[usize],
        replacements: &BTreeMap<usize, LinearEnergyExpr>,
    ) -> Option<(LogicalChannel, usize)> {
        RepeatedLtdBuilder::logical_channels(parsed)
            .into_iter()
            .filter(|channel| channel.power > 1)
            .filter_map(|channel| {
                let degree = channel
                    .members
                    .iter()
                    .map(|local_id| {
                        let orig_id = local_to_orig[*local_id];
                        if replacements.contains_key(&orig_id) {
                            0
                        } else {
                            self.bounds[orig_id]
                        }
                    })
                    .sum::<usize>();
                (degree > 2 || self.sampling_scale_mode.is_active_for_degree(degree))
                    .then_some((channel, degree))
            })
            .next()
    }

    fn channel_replacement_expr(
        edge_id: usize,
        sample: i32,
        rel_sign: i32,
        use_uniform_scale: bool,
    ) -> LinearEnergyExpr {
        let sample = sample * rel_sign;
        if sample == 0 {
            LinearEnergyExpr::zero()
        } else if use_uniform_scale {
            LinearEnergyExpr::uniform_scale(i64::from(sample))
        } else {
            LinearEnergyExpr::ose(EdgeIndex(edge_id), i64::from(sample))
        }
    }

    fn relative_signature_sign(
        reference: &MomentumSignature,
        candidate: &MomentumSignature,
    ) -> Result<i32> {
        if candidate == reference {
            return Ok(1);
        }
        let negated = MomentumSignature {
            loop_signature: reference
                .loop_signature
                .iter()
                .map(|value| -*value)
                .collect(),
            external_signature: reference
                .external_signature
                .iter()
                .map(|value| -*value)
                .collect(),
        };
        if *candidate == negated {
            Ok(-1)
        } else {
            Err(GenerationError::CffHigherEnergyPowerNotImplemented)
        }
    }

    fn interpolation_nodes(degree: usize) -> Vec<i32> {
        let mut nodes = vec![1, -1, 0];
        let mut next = 2;
        while nodes.len() < degree + 1 {
            nodes.push(next);
            if nodes.len() == degree + 1 {
                break;
            }
            nodes.push(-next);
            next += 1;
        }
        nodes.truncate(degree + 1);
        nodes
    }

    fn channel_normal_form_terms(
        poly: &[Rational],
        channel_power: usize,
    ) -> Vec<ChannelNormalFormTerm> {
        let mut combined = BTreeMap::<(usize, usize, usize, usize), Rational>::new();
        for (power, coeff) in poly.iter().cloned().enumerate() {
            if coeff.is_zero() {
                continue;
            }
            let quotient_power = power / 2;
            let parity = power % 2;
            for z_power in 0..=quotient_power {
                let term_coeff = coeff.clone() * binomial(quotient_power, z_power);
                let remaining = channel_power.saturating_sub(z_power);
                let cancelled = z_power.saturating_sub(channel_power);
                let inverse_power = 2 * z_power + parity;
                let key = (remaining, parity, cancelled, inverse_power);
                let entry = combined.entry(key).or_insert(Rational::zero());
                *entry = entry.clone() + term_coeff;
            }
        }
        combined
            .into_iter()
            .filter_map(
                |((remaining_power, parity, cancelled_power, inverse_power), coeff)| {
                    (!coeff.is_zero()).then_some(ChannelNormalFormTerm {
                        remaining_power,
                        parity,
                        cancelled_power,
                        inverse_power,
                        positive_ose_power: 0,
                        coeff,
                    })
                },
            )
            .collect()
    }

    fn channel_uniform_normal_form_terms(
        poly: &[Rational],
        channel_power: usize,
    ) -> Vec<ChannelNormalFormTerm> {
        let mut combined = BTreeMap::<(usize, usize, usize, usize, usize), Rational>::new();
        for (power, coeff) in poly.iter().cloned().enumerate() {
            if coeff.is_zero() {
                continue;
            }
            let quotient_power = power / 2;
            let parity = power % 2;
            for d_power in 0..=quotient_power {
                let term_coeff = coeff.clone() * binomial(quotient_power, d_power);
                let remaining = channel_power.saturating_sub(d_power);
                let cancelled = d_power.saturating_sub(channel_power);
                let inverse_scale_power = power;
                let positive_ose_power = 2 * (quotient_power - d_power);
                let key = (
                    remaining,
                    parity,
                    cancelled,
                    inverse_scale_power,
                    positive_ose_power,
                );
                let entry = combined.entry(key).or_insert(Rational::zero());
                *entry = entry.clone() + term_coeff;
            }
        }
        combined
            .into_iter()
            .filter_map(
                |(
                    (remaining_power, parity, cancelled_power, inverse_power, positive_ose_power),
                    coeff,
                )| {
                    (!coeff.is_zero()).then_some(ChannelNormalFormTerm {
                        remaining_power,
                        parity,
                        cancelled_power,
                        inverse_power,
                        positive_ose_power,
                        coeff,
                    })
                },
            )
            .collect()
    }

    fn remap_known_factor_to_sub(
        &self,
        parsed: &ParsedGraph,
        subparsed: &ParsedGraph,
        factor: &KnownLinearExpr,
    ) -> Result<KnownLinearExpr> {
        let mut out = KnownLinearExpr {
            ose_terms: factor.ose_terms.clone(),
            external_terms: factor.external_terms.clone(),
            constant: factor.constant.clone(),
            uniform_scale_coeff: factor.uniform_scale_coeff.clone(),
            ..KnownLinearExpr::zero()
        };
        for (edge_id, coeff) in &factor.var_terms {
            let expr =
                self.known_signature_expr(subparsed, &parsed.internal_edges[*edge_id].signature)?;
            out = out + expr.mul_rational(coeff.clone())?;
        }
        Ok(out.canonical())
    }

    fn known_signature_expr(
        &self,
        parsed: &ParsedGraph,
        signature: &MomentumSignature,
    ) -> Result<KnownLinearExpr> {
        let signatures = parsed
            .internal_edges
            .iter()
            .map(|edge| edge.signature.clone())
            .collect::<Vec<_>>();
        let edges = (0..signatures.len()).collect::<Vec<_>>();
        let lower_sector = LowerSectorCffBuilder::new(parsed);
        let basis = lower_sector.component_basis_edges(&signatures, &edges);
        let basis_rows = basis
            .iter()
            .map(|edge_id| signatures[*edge_id].loop_signature.clone())
            .collect::<Vec<_>>();
        let coords =
            lower_sector.row_coordinates_in_basis(&basis_rows, &signature.loop_signature)?;

        let mut reconstructed = vec![0i64; signature.loop_signature.len()];
        for (coord, row) in coords.iter().zip(&basis_rows) {
            for (idx, value) in row.iter().enumerate() {
                reconstructed[idx] += i64::from(*coord) * i64::from(*value);
            }
        }
        if reconstructed
            .iter()
            .zip(&signature.loop_signature)
            .any(|(lhs, rhs)| *lhs != i64::from(*rhs))
        {
            return Err(GenerationError::CffHigherEnergyPowerNotImplemented);
        }

        let mut out = KnownLinearExpr::zero();
        let mut external_terms = signature
            .external_signature
            .iter()
            .map(|value| i64::from(*value))
            .collect::<Vec<_>>();
        for (coord, edge_id) in coords.into_iter().zip(basis) {
            if coord == 0 {
                continue;
            }
            out = out + KnownLinearExpr::var(edge_id, i64::from(coord));
            for (external_id, basis_ext_coeff) in
                signatures[edge_id].external_signature.iter().enumerate()
            {
                external_terms[external_id] -= i64::from(coord) * i64::from(*basis_ext_coeff);
            }
        }
        for (external_id, coeff) in external_terms.into_iter().enumerate() {
            if coeff != 0 {
                out = out + KnownLinearExpr::external(external_id, coeff);
            }
        }
        Ok(out.canonical())
    }

    fn copy_expression_surfaces(
        &mut self,
        source: &ThreeDExpression<OrientationID>,
        edge_map: &BTreeMap<usize, usize>,
    ) -> HashMap<HybridSurfaceID, HybridSurfaceID> {
        source
            .surfaces
            .linear_surface_cache
            .iter_enumerated()
            .map(|(id, surface)| {
                (
                    HybridSurfaceID::Linear(id),
                    self.intern_surface(
                        surface.expression.clone().remap_internal_edges(edge_map),
                        surface.origin,
                        surface.numerator_only,
                    ),
                )
            })
            .collect()
    }

    fn intern_surface(
        &mut self,
        surface_expr: LinearEnergyExpr,
        origin: SurfaceOrigin,
        numerator_only: bool,
    ) -> HybridSurfaceID {
        let surface_expr = surface_expr.canonical();
        let kind = classify_surface_kind(&surface_expr);
        let key = (kind, surface_expr.clone());
        if let Some(surface_id) = self.surface_index.get(&key) {
            if !numerator_only && let HybridSurfaceID::Linear(id) = *surface_id {
                self.expression.surfaces.linear_surface_cache[id].origin = SurfaceOrigin::Physical;
                self.expression.surfaces.linear_surface_cache[id].numerator_only = false;
            }
            return *surface_id;
        }

        let id = LinearSurfaceID(self.expression.surfaces.linear_surface_cache.len());
        let surface_id = HybridSurfaceID::Linear(id);
        self.expression
            .surfaces
            .linear_surface_cache
            .push(LinearSurface {
                kind,
                expression: surface_expr,
                origin,
                numerator_only,
            });
        self.surface_index.insert(key, surface_id);
        surface_id
    }

    fn push_variant_for_maps(
        &mut self,
        loop_energy_map: Vec<LinearEnergyExpr>,
        edge_energy_map: Vec<LinearEnergyExpr>,
        variant: crate::expression::CFFVariant,
    ) {
        let (orientation, base_label) = orientation_from_edge_exprs(&edge_energy_map);
        if let Some(existing) = self
            .expression
            .orientations
            .iter_mut()
            .find(|orientation_expr| {
                orientation_expr.data.label.as_deref() == Some(base_label.as_str())
                    && orientation_expr.loop_energy_map == loop_energy_map
                    && orientation_expr.edge_energy_map == edge_energy_map
            })
        {
            existing.variants.push(variant);
            return;
        }
        let mut data = OrientationData::new(orientation);
        data.label = Some(base_label);
        self.expression.orientations.push(OrientationExpression {
            data,
            loop_energy_map,
            edge_energy_map,
            variants: vec![variant],
        });
    }

    fn finalize_numerator_map_labels(&mut self) {
        assign_numerator_map_labels(&mut self.expression.orientations);
    }
}

fn map_surface_id(
    surface_id: HybridSurfaceID,
    surface_map: &HashMap<HybridSurfaceID, HybridSurfaceID>,
) -> HybridSurfaceID {
    match surface_id {
        HybridSurfaceID::Linear(_) => surface_map[&surface_id],
        HybridSurfaceID::Unit | HybridSurfaceID::Infinite => surface_id,
        HybridSurfaceID::Esurface(_) | HybridSurfaceID::Hsurface(_) => surface_id,
    }
}

fn rational_from_coefficient(value: &symbolica::atom::Atom) -> Rational {
    value.rational_coeff()
}

fn union_nodes(parent: &mut BTreeMap<usize, usize>, lhs: usize, rhs: usize) {
    let lhs_root = find_node_root(parent, lhs);
    let rhs_root = find_node_root(parent, rhs);
    if lhs_root != rhs_root {
        parent.insert(rhs_root, lhs_root);
    }
}

fn find_node_root(parent: &mut BTreeMap<usize, usize>, node: usize) -> usize {
    let mut current = node;
    while parent[&current] != current {
        current = parent[&current];
    }
    let root = current;
    let mut current = node;
    while parent[&current] != current {
        let next = parent[&current];
        parent.insert(current, root);
        current = next;
    }
    root
}

fn contact_nodes(bound: usize) -> Vec<i32> {
    let needed = bound.saturating_sub(1);
    let mut nodes = vec![0];
    let mut next = 2;
    while nodes.len() < needed {
        nodes.push(next);
        if nodes.len() == needed {
            break;
        }
        nodes.push(-next);
        next += 1;
    }
    nodes.truncate(needed);
    nodes
}

fn lagrange_basis(nodes: &[i32], index: usize) -> Vec<Rational> {
    let mut poly = vec![Rational::one()];
    let mut denominator = Rational::one();
    let xj = Rational::from(nodes[index]);
    for (other_index, xk_int) in nodes.iter().enumerate() {
        if other_index == index {
            continue;
        }
        let xk = Rational::from(*xk_int);
        poly = poly_mul(&poly, &[-xk.clone(), Rational::one()]);
        denominator *= xj.clone() - xk;
    }
    poly_scale(&poly, Rational::one() / denominator)
}

fn contact_weight_polys(bound: usize) -> BTreeMap<i32, Vec<Rational>> {
    if bound <= 1 {
        return BTreeMap::new();
    }
    let nodes = contact_nodes(bound);
    let mut weights = BTreeMap::<i32, Vec<Rational>>::new();
    for (index, sample) in nodes.iter().enumerate() {
        let a = Rational::from(*sample);
        let denominator = a.clone() * a.clone() - Rational::one();
        let basis = lagrange_basis(&nodes, index);
        let contributions = [
            (*sample, Rational::one() / denominator.clone()),
            (
                1,
                -((a.clone() + Rational::one()) / Rational::from(2)) / denominator.clone(),
            ),
            (
                -1,
                -((Rational::one() - a) / Rational::from(2)) / denominator,
            ),
        ];
        for (node, scale) in contributions {
            let contribution = poly_scale(&basis, scale);
            let entry = weights.entry(node).or_default();
            *entry = poly_add(entry, &contribution);
        }
    }
    weights.retain(|_, poly| !poly.is_empty());
    weights
}

fn poly_add(lhs: &[Rational], rhs: &[Rational]) -> Vec<Rational> {
    let mut out = vec![Rational::zero(); lhs.len().max(rhs.len())];
    for (index, value) in lhs.iter().enumerate() {
        out[index] = out[index].clone() + value.clone();
    }
    for (index, value) in rhs.iter().enumerate() {
        out[index] = out[index].clone() + value.clone();
    }
    trim_poly(out)
}

fn poly_scale(poly: &[Rational], scale: Rational) -> Vec<Rational> {
    trim_poly(
        poly.iter()
            .map(|value| value.clone() * scale.clone())
            .collect(),
    )
}

fn poly_mul(lhs: &[Rational], rhs: &[Rational]) -> Vec<Rational> {
    if lhs.is_empty() || rhs.is_empty() {
        return Vec::new();
    }
    let mut out = vec![Rational::zero(); lhs.len() + rhs.len() - 1];
    for (lhs_index, lhs_value) in lhs.iter().enumerate() {
        for (rhs_index, rhs_value) in rhs.iter().enumerate() {
            out[lhs_index + rhs_index] =
                out[lhs_index + rhs_index].clone() + lhs_value.clone() * rhs_value.clone();
        }
    }
    trim_poly(out)
}

fn trim_poly(mut poly: Vec<Rational>) -> Vec<Rational> {
    while poly.last().is_some_and(|value| value.is_zero()) {
        poly.pop();
    }
    poly
}

#[derive(Debug, Clone)]
struct LowerSectorComponent {
    basis_edges: Vec<usize>,
    local_to_sub: Vec<usize>,
    expression: ThreeDExpression<OrientationID>,
    prefactor_correction: symbolica::atom::Atom,
}

#[derive(Debug, Clone)]
struct LowerSectorPartial {
    coeff: Rational,
    half_edges: Vec<usize>,
    chain: Vec<HybridSurfaceID>,
    numerator_surfaces: Vec<HybridSurfaceID>,
    targets: BTreeMap<usize, LinearEnergyExpr>,
    edge_exprs: BTreeMap<usize, LinearEnergyExpr>,
}

struct LowerSectorCffBuilder<'a> {
    parsed: &'a ParsedGraph,
    expression: ThreeDExpression<OrientationID>,
    surface_index: HashMap<(LinearSurfaceKind, LinearEnergyExpr), HybridSurfaceID>,
}

impl<'a> LowerSectorCffBuilder<'a> {
    fn new(parsed: &'a ParsedGraph) -> Self {
        Self {
            parsed,
            expression: ThreeDExpression::new_empty(),
            surface_index: HashMap::new(),
        }
    }

    fn build(mut self) -> Result<ThreeDExpression<OrientationID>> {
        if self.parsed.internal_edges.is_empty() {
            return Ok(self.expression);
        }
        let signatures = self.signatures();
        let rows = signatures
            .iter()
            .map(|signature| {
                signature
                    .loop_signature
                    .iter()
                    .map(|value| i64::from(*value))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        if rank_i64(&rows) == 0 {
            return Err(GenerationError::CffHigherEnergyPowerNotImplemented);
        }

        let components = self.component_bundles(&signatures)?;
        let mut partials = vec![LowerSectorPartial {
            coeff: Rational::one(),
            half_edges: Vec::new(),
            chain: Vec::new(),
            numerator_surfaces: Vec::new(),
            targets: BTreeMap::new(),
            edge_exprs: BTreeMap::new(),
        }];

        for component in &components {
            let edge_map = component
                .local_to_sub
                .iter()
                .enumerate()
                .map(|(local_id, sub_id)| (local_id, *sub_id))
                .collect::<BTreeMap<_, _>>();
            let surface_map = self.copy_expression_surfaces(&component.expression, &edge_map);
            let basis_sub = component
                .basis_edges
                .iter()
                .copied()
                .collect::<BTreeSet<_>>();
            let mut next_partials = Vec::new();
            for partial in &partials {
                for orientation in &component.expression.orientations {
                    for variant in &orientation.variants {
                        let mut item = partial.clone();
                        item.coeff = item.coeff
                            * rational_from_coefficient(&component.prefactor_correction)
                            * rational_from_coefficient(&variant.prefactor);
                        item.half_edges.extend(
                            variant
                                .half_edges
                                .iter()
                                .map(|edge| edge_map.get(&edge.0).copied().unwrap_or(edge.0)),
                        );
                        item.numerator_surfaces.extend(
                            variant
                                .numerator_surfaces
                                .iter()
                                .map(|surface| map_surface_id(*surface, &surface_map)),
                        );
                        for chain in denominator_tree_chains(&variant.denominator) {
                            let mut branched = item.clone();
                            branched.chain.extend(
                                chain
                                    .into_iter()
                                    .map(|sid| map_surface_id(sid, &surface_map)),
                            );
                            for (local_id, sub_id) in &edge_map {
                                let lifted = orientation.edge_energy_map[*local_id]
                                    .clone()
                                    .remap_internal_edges(&edge_map);
                                branched.edge_exprs.insert(*sub_id, lifted.clone());
                                if basis_sub.contains(sub_id) {
                                    branched.targets.insert(*sub_id, lifted);
                                }
                            }
                            next_partials.push(branched);
                        }
                    }
                }
            }
            partials = next_partials;
        }

        let global_basis = components
            .iter()
            .flat_map(|component| component.basis_edges.iter().copied())
            .collect::<Vec<_>>();
        let target_template = vec![LinearEnergyExpr::zero(); self.parsed.internal_edges.len()];
        for partial in partials {
            if partial.coeff.is_zero() {
                continue;
            }
            let mut targets = target_template.clone();
            for edge_id in &global_basis {
                if let Some(expr) = partial.targets.get(edge_id) {
                    targets[*edge_id] = expr.clone();
                }
            }
            let loop_exprs = solve_loop_energy_particular_from_target_edge_exprs(
                &signatures,
                &global_basis,
                &targets,
            )?;
            let mut edge_exprs = target_template.clone();
            for (edge_id, expr) in partial.edge_exprs {
                edge_exprs[edge_id] = expr;
            }
            let mut half_edges = partial.half_edges;
            half_edges.sort_unstable();
            self.push_variant_for_maps(
                loop_exprs,
                edge_exprs,
                crate::expression::CFFVariant {
                    origin: Some("lower_sector_cff_e_surface_component_product".to_string()),
                    prefactor: rational_to_coefficient(partial.coeff)?,
                    half_edges: half_edges.into_iter().map(EdgeIndex).collect(),
                    uniform_scale_power: 0,
                    numerator_surfaces: partial.numerator_surfaces,
                    denominator: denominator_tree_from_chain(&partial.chain),
                },
            );
        }
        self.finalize_numerator_map_labels();
        Ok(self.expression)
    }

    fn signatures(&self) -> Vec<MomentumSignature> {
        self.parsed
            .internal_edges
            .iter()
            .map(|edge| edge.signature.clone())
            .collect()
    }

    fn component_bundles(
        &self,
        signatures: &[MomentumSignature],
    ) -> Result<Vec<LowerSectorComponent>> {
        self.vector_matroid_components(signatures)
            .into_iter()
            .map(|edges| self.component_bundle(signatures, edges))
            .collect()
    }

    fn component_bundle(
        &self,
        signatures: &[MomentumSignature],
        edges: Vec<usize>,
    ) -> Result<LowerSectorComponent> {
        let rank = rank_i64(
            &edges
                .iter()
                .map(|edge_id| {
                    signatures[*edge_id]
                        .loop_signature
                        .iter()
                        .map(|value| i64::from(*value))
                        .collect()
                })
                .collect::<Vec<Vec<_>>>(),
        );
        if rank == 0 {
            return Err(GenerationError::CffHigherEnergyPowerNotImplemented);
        }
        let basis_edges = self.component_basis_edges(signatures, &edges);
        let (component_parsed, local_to_sub) =
            self.project_component_parsed(signatures, &edges, &basis_edges)?;
        let expression = if component_parsed.internal_edges.len() == rank {
            generate_pure_ltd_expression_from_parsed(&component_parsed)?
        } else {
            generate_pure_cff_expression_from_parsed_with_duplicate_sign(&component_parsed, false)?
        };
        let is_auxiliary = component_parsed.internal_edges.len() != rank;
        let prefactor_correction =
            if is_auxiliary && (component_parsed.internal_edges.len() - rank + 1) % 2 == 1 {
                rational_coeff_new(-1, 1)
            } else {
                rational_coeff_one()
            };

        Ok(LowerSectorComponent {
            basis_edges,
            local_to_sub,
            expression,
            prefactor_correction,
        })
    }

    fn project_component_parsed(
        &self,
        signatures: &[MomentumSignature],
        edges: &[usize],
        basis_edges: &[usize],
    ) -> Result<(ParsedGraph, Vec<usize>)> {
        let basis_rows = basis_edges
            .iter()
            .map(|edge_id| signatures[*edge_id].loop_signature.clone())
            .collect::<Vec<_>>();
        let projected = edges
            .iter()
            .map(|edge_id| {
                self.row_coordinates_in_basis(&basis_rows, &signatures[*edge_id].loop_signature)
            })
            .collect::<Result<Vec<_>>>()?;
        let rank = basis_edges.len();
        let loop_names = (0..rank).map(|idx| format!("ell{idx}")).collect::<Vec<_>>();

        if edges.len() == rank {
            let internal_edges = edges
                .iter()
                .zip(projected)
                .enumerate()
                .map(|(local_id, (edge_id, loop_coeffs))| {
                    let edge = &self.parsed.internal_edges[*edge_id];
                    ParsedGraphInternalEdge {
                        edge_id: local_id,
                        tail: local_id,
                        head: (local_id + 1) % edges.len(),
                        label: edge.label.clone(),
                        mass_key: edge.mass_key.clone(),
                        signature: MomentumSignature {
                            loop_signature: loop_coeffs,
                            external_signature: edge.signature.external_signature.clone(),
                        },
                        had_pow: edge.had_pow,
                    }
                })
                .collect::<Vec<_>>();
            let node_name_to_internal = (0..edges.len())
                .map(|node| (format!("aux{node}"), node))
                .collect::<BTreeMap<_, _>>();
            return Ok((
                ParsedGraph {
                    internal_edges,
                    external_edges: Vec::new(),
                    initial_state_cut_edges: Vec::new(),
                    loop_names,
                    external_names: self.parsed.external_names.clone(),
                    node_name_to_internal,
                },
                edges.to_vec(),
            ));
        }

        let projected_signatures = edges
            .iter()
            .zip(&projected)
            .map(|(edge_id, loop_coeffs)| MomentumSignature {
                loop_signature: loop_coeffs.clone(),
                external_signature: self.parsed.internal_edges[*edge_id]
                    .signature
                    .external_signature
                    .clone(),
            })
            .collect::<Vec<_>>();
        let marker_masses = (0..edges.len())
            .map(|local_id| format!("__lower_cff_edge_{local_id}"))
            .collect::<Vec<_>>();
        let extracted = ExtractedSignatureExpression {
            signatures: projected_signatures,
            loop_names: loop_names.clone(),
            external_names: self.parsed.external_names.clone(),
            masses: marker_masses.clone(),
        };
        let reconstructed = reconstruct_parsed_graph(
            &extracted,
            &ReconstructDotOptions {
                require_connected: true,
                minimize_external_legs: true,
                format: ReconstructDotFormat::Gammaloop,
                ..ReconstructDotOptions::default()
            },
        )
        .map_err(|_| GenerationError::CffHigherEnergyPowerNotImplemented)?;
        let marker_to_sub = marker_masses
            .iter()
            .zip(edges)
            .map(|(marker, edge_id)| (marker.clone(), *edge_id))
            .collect::<BTreeMap<_, _>>();
        let mut local_to_sub = Vec::new();
        let mut internal_edges = Vec::new();
        for (new_id, edge) in reconstructed.internal_edges.iter().enumerate() {
            let Some(marker) = &edge.mass_key else {
                return Err(GenerationError::CffHigherEnergyPowerNotImplemented);
            };
            let Some(sub_id) = marker_to_sub.get(marker).copied() else {
                return Err(GenerationError::CffHigherEnergyPowerNotImplemented);
            };
            let original = &self.parsed.internal_edges[sub_id];
            local_to_sub.push(sub_id);
            internal_edges.push(ParsedGraphInternalEdge {
                edge_id: new_id,
                tail: edge.tail,
                head: edge.head,
                label: original.label.clone(),
                mass_key: original.mass_key.clone(),
                signature: edge.signature.clone(),
                had_pow: original.had_pow,
            });
        }
        Ok((
            ParsedGraph {
                internal_edges,
                external_edges: reconstructed.external_edges,
                initial_state_cut_edges: Vec::new(),
                loop_names,
                external_names: self.parsed.external_names.clone(),
                node_name_to_internal: reconstructed.node_name_to_internal,
            },
            local_to_sub,
        ))
    }

    fn vector_matroid_components(&self, signatures: &[MomentumSignature]) -> Vec<Vec<usize>> {
        let rows = signatures
            .iter()
            .map(|signature| signature.loop_signature.clone())
            .collect::<Vec<_>>();
        let nonzero = rows
            .iter()
            .enumerate()
            .filter_map(|(edge_id, row)| row.iter().any(|value| *value != 0).then_some(edge_id))
            .collect::<Vec<_>>();
        let mut parent = nonzero
            .iter()
            .map(|edge_id| (*edge_id, *edge_id))
            .collect::<BTreeMap<_, _>>();
        let mut basis = Vec::<usize>::new();
        let mut basis_rows = Vec::<Vec<i32>>::new();
        for edge_id in &nonzero {
            let mut trial = basis_rows.clone();
            trial.push(rows[*edge_id].clone());
            if rank_i64(
                &trial
                    .iter()
                    .map(|row| row.iter().map(|value| i64::from(*value)).collect())
                    .collect::<Vec<Vec<_>>>(),
            ) > basis_rows.len()
            {
                basis.push(*edge_id);
                basis_rows.push(rows[*edge_id].clone());
                continue;
            }
            if let Ok(coords) = self.row_coordinates_in_basis(&basis_rows, &rows[*edge_id]) {
                for (basis_edge, coeff) in basis.iter().zip(coords) {
                    if coeff != 0 {
                        union_nodes(&mut parent, *edge_id, *basis_edge);
                    }
                }
            }
        }

        let mut grouped = BTreeMap::<usize, Vec<usize>>::new();
        for edge_id in nonzero {
            let root = find_node_root(&mut parent, edge_id);
            grouped.entry(root).or_default().push(edge_id);
        }
        let mut components = grouped
            .into_values()
            .map(|mut group| {
                group.sort_unstable();
                group
            })
            .collect::<Vec<_>>();
        components.extend(rows.iter().enumerate().filter_map(|(edge_id, row)| {
            (!row.iter().any(|value| *value != 0)).then_some(vec![edge_id])
        }));
        components.sort_by_key(|component| component[0]);
        components
    }

    fn component_basis_edges(
        &self,
        signatures: &[MomentumSignature],
        edges: &[usize],
    ) -> Vec<usize> {
        let mut selected = Vec::new();
        let mut selected_rows = Vec::<Vec<i32>>::new();
        let mut current_rank = 0usize;
        for edge_id in edges {
            let mut trial = selected_rows.clone();
            trial.push(signatures[*edge_id].loop_signature.clone());
            let trial_rank = rank_i64(
                &trial
                    .iter()
                    .map(|row| row.iter().map(|value| i64::from(*value)).collect())
                    .collect::<Vec<Vec<_>>>(),
            );
            if trial_rank > current_rank {
                selected.push(*edge_id);
                selected_rows.push(signatures[*edge_id].loop_signature.clone());
                current_rank = trial_rank;
            }
        }
        selected
    }

    fn row_coordinates_in_basis(&self, basis_rows: &[Vec<i32>], row: &[i32]) -> Result<Vec<i32>> {
        let rank = basis_rows.len();
        if rank == 0 {
            return if row.iter().all(|value| *value == 0) {
                Ok(Vec::new())
            } else {
                Err(GenerationError::SingularBasis)
            };
        }
        for columns in (0..row.len()).combinations(rank) {
            let square = columns
                .iter()
                .map(|column| {
                    basis_rows
                        .iter()
                        .map(|basis_row| Rational::from(basis_row[*column]))
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
            if rank_rational(&square) != rank {
                continue;
            }
            let rhs = columns
                .iter()
                .map(|column| Rational::from(row[*column]))
                .collect::<Vec<_>>();
            let coords =
                solve_rational_system(square, rhs).ok_or(GenerationError::SingularBasis)?;
            if coords.iter().any(|coord| {
                coord
                    .to_i64_pair()
                    .is_none_or(|(_, denominator)| denominator != 1)
            }) {
                return Err(GenerationError::NonIntegralEnergyMap);
            }
            return coords
                .into_iter()
                .map(|coord| {
                    let (numerator, _) = coord
                        .to_i64_pair()
                        .ok_or(GenerationError::CoefficientOutOfRange)?;
                    i32::try_from(numerator).map_err(|_| GenerationError::CoefficientOutOfRange)
                })
                .collect();
        }
        Err(GenerationError::SingularBasis)
    }

    fn copy_expression_surfaces(
        &mut self,
        source: &ThreeDExpression<OrientationID>,
        edge_map: &BTreeMap<usize, usize>,
    ) -> HashMap<HybridSurfaceID, HybridSurfaceID> {
        source
            .surfaces
            .linear_surface_cache
            .iter_enumerated()
            .map(|(id, surface)| {
                (
                    HybridSurfaceID::Linear(id),
                    self.intern_surface(
                        surface.expression.clone().remap_internal_edges(edge_map),
                        surface.origin,
                        surface.numerator_only,
                    ),
                )
            })
            .collect()
    }

    fn intern_surface(
        &mut self,
        surface_expr: LinearEnergyExpr,
        origin: SurfaceOrigin,
        numerator_only: bool,
    ) -> HybridSurfaceID {
        let surface_expr = surface_expr.canonical();
        let kind = classify_surface_kind(&surface_expr);
        let key = (kind, surface_expr.clone());
        if let Some(surface_id) = self.surface_index.get(&key) {
            if !numerator_only && let HybridSurfaceID::Linear(id) = *surface_id {
                self.expression.surfaces.linear_surface_cache[id].origin = SurfaceOrigin::Physical;
                self.expression.surfaces.linear_surface_cache[id].numerator_only = false;
            }
            return *surface_id;
        }

        let id = LinearSurfaceID(self.expression.surfaces.linear_surface_cache.len());
        let surface_id = HybridSurfaceID::Linear(id);
        self.expression
            .surfaces
            .linear_surface_cache
            .push(LinearSurface {
                kind,
                expression: surface_expr,
                origin,
                numerator_only,
            });
        self.surface_index.insert(key, surface_id);
        surface_id
    }

    fn push_variant_for_maps(
        &mut self,
        loop_energy_map: Vec<LinearEnergyExpr>,
        edge_energy_map: Vec<LinearEnergyExpr>,
        variant: crate::expression::CFFVariant,
    ) {
        let (orientation, base_label) = orientation_from_edge_exprs(&edge_energy_map);
        if let Some(existing) = self
            .expression
            .orientations
            .iter_mut()
            .find(|orientation_expr| {
                orientation_expr.data.label.as_deref() == Some(base_label.as_str())
                    && orientation_expr.loop_energy_map == loop_energy_map
                    && orientation_expr.edge_energy_map == edge_energy_map
            })
        {
            existing.variants.push(variant);
            return;
        }
        let mut data = OrientationData::new(orientation);
        data.label = Some(base_label);
        self.expression.orientations.push(OrientationExpression {
            data,
            loop_energy_map,
            edge_energy_map,
            variants: vec![variant],
        });
    }

    fn finalize_numerator_map_labels(&mut self) {
        assign_numerator_map_labels(&mut self.expression.orientations);
    }
}

struct RepeatedLtdBuilder<'a> {
    parsed: &'a ParsedGraph,
    options: &'a Generate3DExpressionOptions,
    signatures: Vec<MomentumSignature>,
    channels: Vec<LogicalChannel>,
    expression: ThreeDExpression<OrientationID>,
}

impl<'a> RepeatedLtdBuilder<'a> {
    fn new(parsed: &'a ParsedGraph, options: &'a Generate3DExpressionOptions) -> Self {
        let signatures = parsed
            .internal_edges
            .iter()
            .map(|edge| edge.signature.clone())
            .collect::<Vec<_>>();
        let channels = Self::logical_channels(parsed);
        Self {
            parsed,
            options,
            signatures,
            channels,
            expression: ThreeDExpression::new_empty(),
        }
    }

    fn build(mut self) -> Result<ThreeDExpression<OrientationID>> {
        let n_loops = self
            .signatures
            .first()
            .map(|signature| signature.loop_signature.len())
            .unwrap_or(0);
        let collapsed_signatures = self
            .channels
            .iter()
            .map(|channel| self.signatures[channel.rep_edge].loop_signature.clone())
            .collect::<Vec<_>>();
        let residues = ltd_residues(&collapsed_signatures, &vec![ContourClosure::Below; n_loops])?;
        let bounds = if self.options.energy_degree_bounds.is_empty() {
            None
        } else {
            Some(crate::energy_bounds::normalize_energy_degree_bounds(
                &self.options.energy_degree_bounds,
                self.signatures.len(),
            )?)
        };

        for residue in residues {
            self.add_residue(&residue, bounds.as_deref())?;
        }
        self.finalize_numerator_map_labels();
        Ok(self.expression)
    }

    fn add_residue(&mut self, residue: &crate::Residue, bounds: Option<&[usize]>) -> Result<()> {
        let basis_logical = residue.basis.clone();
        let cut_signs = residue.sigmas.clone();
        let basis_orig = basis_logical
            .iter()
            .map(|logical_idx| self.channels[*logical_idx].rep_edge)
            .collect::<Vec<_>>();
        let powers = basis_logical
            .iter()
            .map(|logical_idx| self.channels[*logical_idx].power)
            .collect::<Vec<_>>();
        let alpha = powers
            .iter()
            .map(|power| power.saturating_sub(1))
            .collect::<Vec<_>>();

        let solver = EnergySolver::new(&self.signatures);
        let mut targets = vec![LinearEnergyExpr::zero(); self.signatures.len()];
        for (edge, sigma) in basis_orig.iter().zip(&cut_signs) {
            targets[*edge] = LinearEnergyExpr::ose(EdgeIndex(*edge), i64::from(*sigma));
        }
        let loop_exprs = solver.solve_from_target_edges(&basis_orig, &targets)?;
        let mut edge_exprs = solver.edge_q0_from_loop_exprs(&loop_exprs);
        apply_initial_state_cut_edge_energy_exprs(self.parsed, &mut edge_exprs);
        let (loop_derivs, edge_derivs) = solver.derivative_matrices(&basis_orig)?;

        let mut surface_interner = LinearSurfaceInterner::new(&mut self.expression);
        let mut factors = Vec::new();
        for (basis_position, ((logical_idx, sigma), power)) in basis_logical
            .iter()
            .zip(&cut_signs)
            .zip(&powers)
            .enumerate()
        {
            let mut derivs = vec![Rational::zero(); basis_logical.len()];
            derivs[basis_position] = Rational::from(*sigma);
            factors.push(DenominatorFactor::Cut {
                edge: self.channels[*logical_idx].rep_edge,
                power: *power,
                derivs,
            });
        }

        let basis_set = basis_logical.iter().copied().collect::<BTreeSet<_>>();
        for (logical_idx, channel) in self.channels.iter().enumerate() {
            if basis_set.contains(&logical_idx) {
                continue;
            }
            let rep = channel.rep_edge;
            let q = edge_exprs[rep].clone();
            let minus = surface_interner
                .intern(q.clone() - LinearEnergyExpr::ose(EdgeIndex(rep), 1), false);
            let plus = surface_interner.intern(q + LinearEnergyExpr::ose(EdgeIndex(rep), 1), false);
            factors.push(DenominatorFactor::Surface {
                surface: minus,
                power: channel.power,
                derivs: edge_derivs[rep].clone(),
            });
            factors.push(DenominatorFactor::Surface {
                surface: plus,
                power: channel.power,
                derivs: edge_derivs[rep].clone(),
            });
        }
        drop(surface_interner);

        let mut residue_sign = Rational::from(if residue.sign >= 0 { 1 } else { -1 });
        for (sigma, order) in cut_signs.iter().zip(&alpha) {
            if order % 2 == 1 {
                residue_sign *= Rational::from(*sigma);
            }
        }
        let residue_norm = Rational::one() / multi_factorial(&alpha);
        let beta_candidates = if bounds.is_some() {
            multiindices_leq(&alpha)
        } else {
            multiindices_leq(&alpha)
                .into_iter()
                .filter(|beta| beta.iter().sum::<usize>() <= 1)
                .collect()
        };

        for beta in beta_candidates {
            let gamma = alpha
                .iter()
                .zip(&beta)
                .map(|(a, b)| a - b)
                .collect::<Vec<_>>();
            let leibniz = alpha
                .iter()
                .zip(&beta)
                .fold(Rational::one(), |acc, (a, b)| acc * binomial(*a, *b));
            let numerator_samples = if let Some(bounds) = bounds {
                self.bounded_numerator_samples(
                    &beta,
                    bounds,
                    &basis_logical,
                    &cut_signs,
                    &edge_derivs,
                )?
            } else {
                self.affine_numerator_samples(
                    &beta,
                    &basis_logical,
                    &cut_signs,
                    &alpha,
                    &loop_exprs,
                    &edge_exprs,
                    &loop_derivs,
                    &edge_derivs,
                )?
            };
            if numerator_samples.is_empty() {
                continue;
            }
            let denominator_terms = DenominatorDerivativeExpander::new(&factors).terms(&gamma);
            if denominator_terms.is_empty() {
                continue;
            }

            for numerator_sample in &numerator_samples {
                let mut accumulated = VariantAccumulator::default();
                for denominator_term in &denominator_terms {
                    let coeff = residue_sign.clone()
                        * residue_norm.clone()
                        * leibniz.clone()
                        * numerator_sample.coeff.clone()
                        * denominator_term.coeff.clone();
                    if coeff.is_zero() {
                        continue;
                    }
                    let mut half_edges = denominator_term.half_edges.clone();
                    half_edges.extend(numerator_sample.extra_half_edges.iter().copied());
                    half_edges.sort_unstable();
                    accumulated.add(
                        rational_to_coefficient(coeff)?,
                        half_edges,
                        denominator_term.chain.clone(),
                    );
                }

                for item in accumulated.variants {
                    let variant = crate::expression::CFFVariant {
                        origin: Some(format!(
                            "ltd_confluent:{}:beta={}:gamma={}",
                            numerator_sample.label,
                            beta.iter().map(usize::to_string).join(","),
                            gamma.iter().map(usize::to_string).join(",")
                        )),
                        prefactor: item.prefactor,
                        half_edges: item.half_edges.into_iter().map(EdgeIndex).collect(),
                        uniform_scale_power: numerator_sample.uniform_scale_power,
                        numerator_surfaces: Vec::new(),
                        denominator: denominator_tree_from_chains(&item.chains),
                    };
                    self.push_variant_for_numerator_sample(numerator_sample, variant);
                }
            }
        }

        Ok(())
    }

    fn logical_channels(parsed: &ParsedGraph) -> Vec<LogicalChannel> {
        let repeated = repeated_groups(parsed);
        let group_by_member = repeated
            .iter()
            .flat_map(|group| group.edge_ids.iter().map(move |edge_id| (*edge_id, group)))
            .collect::<HashMap<_, _>>();
        let mut seen_groups = BTreeSet::new();
        let mut channels = Vec::new();
        for edge in &parsed.internal_edges {
            if parsed.is_initial_state_cut_edge(edge.edge_id) {
                continue;
            }
            let Some(group) = group_by_member.get(&edge.edge_id) else {
                channels.push(LogicalChannel {
                    rep_edge: edge.edge_id,
                    members: vec![edge.edge_id],
                    power: 1,
                });
                continue;
            };
            if !seen_groups.insert(group.edge_ids.clone()) {
                continue;
            }
            channels.push(LogicalChannel {
                rep_edge: group.edge_ids[0],
                members: group.edge_ids.clone(),
                power: group.edge_ids.len(),
            });
        }
        channels
    }

    fn push_variant_for_numerator_sample(
        &mut self,
        sample: &NumeratorSample,
        variant: crate::expression::CFFVariant,
    ) {
        let (orientation, base_label) = orientation_from_edge_exprs(&sample.edge_exprs);
        if let Some(existing) = self
            .expression
            .orientations
            .iter_mut()
            .find(|orientation_expr| {
                orientation_expr.data.label.as_deref() == Some(base_label.as_str())
                    && orientation_expr.loop_energy_map == sample.loop_exprs
                    && orientation_expr.edge_energy_map == sample.edge_exprs
            })
        {
            existing.variants.push(variant);
            return;
        }
        let mut data = OrientationData::new(orientation);
        data.label = Some(base_label);
        self.expression.orientations.push(OrientationExpression {
            data,
            loop_energy_map: sample.loop_exprs.clone(),
            edge_energy_map: sample.edge_exprs.clone(),
            variants: vec![variant],
        });
    }

    fn finalize_numerator_map_labels(&mut self) {
        assign_numerator_map_labels(&mut self.expression.orientations);
    }

    #[allow(clippy::too_many_arguments)]
    fn affine_numerator_samples(
        &self,
        beta: &[usize],
        basis_logical: &[usize],
        cut_signs: &[i32],
        alpha: &[usize],
        base_loop_exprs: &[LinearEnergyExpr],
        base_edge_exprs: &[LinearEnergyExpr],
        loop_derivs: &[Vec<Rational>],
        edge_derivs: &[Vec<Rational>],
    ) -> Result<Vec<NumeratorSample>> {
        if beta.iter().any(|order| *order > 1) || beta.iter().sum::<usize>() > 1 {
            return Ok(Vec::new());
        }
        let samples = self.physical_sample_orbit(basis_logical, cut_signs, alpha)?;
        if samples.is_empty() {
            return Ok(Vec::new());
        }

        let internal_alias = self.internal_alias();
        let sample_component_maps = samples
            .iter()
            .map(|sample| {
                component_coeff_maps(&sample.loop_exprs, &sample.edge_exprs, &internal_alias)
            })
            .collect::<Vec<_>>();
        let n_components = sample_component_maps[0].len();
        let (constant_target, target_maps, extra_half_edges, sample_label) = if beta
            .iter()
            .all(|order| *order == 0)
        {
            (
                Rational::one(),
                component_coeff_maps(base_loop_exprs, base_edge_exprs, &internal_alias),
                Vec::new(),
                "n0".to_string(),
            )
        } else {
            let Some(basis_position) = beta.iter().position(|order| *order == 1) else {
                return Ok(Vec::new());
            };
            let logical_idx = basis_logical[basis_position];
            let rep_edge = self.channels[logical_idx].rep_edge;
            let aliased_rep = *internal_alias.get(&rep_edge).unwrap_or(&rep_edge);
            let mut target_maps = Vec::new();
            for row in loop_derivs {
                let val = row[basis_position].clone();
                target_maps.push(if val.is_zero() {
                    BTreeMap::new()
                } else {
                    BTreeMap::from([(CoeffSymbol::Internal(aliased_rep), val * Rational::from(2))])
                });
            }
            for row in edge_derivs {
                let val = row[basis_position].clone();
                target_maps.push(if val.is_zero() {
                    BTreeMap::new()
                } else {
                    BTreeMap::from([(CoeffSymbol::Internal(aliased_rep), val * Rational::from(2))])
                });
            }
            (
                Rational::zero(),
                target_maps,
                vec![rep_edge],
                format!("nD{basis_position}phys"),
            )
        };

        if target_maps.len() != n_components {
            return Err(GenerationError::NumeratorSampleOrbitIncomplete);
        }
        let weights =
            AffineWeightSystem::new(&sample_component_maps, &target_maps, constant_target)
                .solve()?;
        Ok(weights
            .into_iter()
            .zip(samples)
            .filter_map(|(coeff, sample)| {
                (!coeff.is_zero()).then_some(NumeratorSample {
                    coeff,
                    extra_half_edges: extra_half_edges.clone(),
                    loop_exprs: sample.loop_exprs,
                    edge_exprs: sample.edge_exprs,
                    label: format!("{sample_label}:{}", sample.label),
                    uniform_scale_power: 0,
                })
            })
            .collect())
    }

    fn bounded_numerator_samples(
        &self,
        beta: &[usize],
        bounds: &[usize],
        basis_logical: &[usize],
        cut_signs: &[i32],
        edge_derivs: &[Vec<Rational>],
    ) -> Result<Vec<NumeratorSample>> {
        let basis_orig = basis_logical
            .iter()
            .map(|logical_idx| self.channels[*logical_idx].rep_edge)
            .collect::<Vec<_>>();
        if beta.iter().all(|order| *order == 0) {
            let solver = EnergySolver::new(&self.signatures);
            let mut targets = vec![LinearEnergyExpr::zero(); self.signatures.len()];
            for (edge, sigma) in basis_orig.iter().zip(cut_signs) {
                targets[*edge] = LinearEnergyExpr::ose(EdgeIndex(*edge), i64::from(*sigma));
            }
            let loop_exprs = solver.solve_from_target_edges(&basis_orig, &targets)?;
            let edge_exprs = solver.edge_q0_from_loop_exprs(&loop_exprs);
            return Ok(vec![NumeratorSample {
                coeff: Rational::one(),
                extra_half_edges: Vec::new(),
                loop_exprs,
                edge_exprs,
                label: "n0:fd0".to_string(),
                uniform_scale_power: 0,
            }]);
        }

        let degree_by_basis = self.basis_variable_degree_bounds(bounds, edge_derivs);
        let channel_degree_by_basis = basis_logical
            .iter()
            .map(|logical_idx| {
                self.channels[*logical_idx]
                    .members
                    .iter()
                    .map(|member| bounds[*member])
                    .sum::<usize>()
            })
            .collect::<Vec<_>>();
        let mut per_axis = Vec::new();
        for (basis_position, order) in beta.iter().enumerate() {
            if *order == 0 {
                continue;
            }
            let degree = degree_by_basis[basis_position];
            if *order > degree {
                return Ok(Vec::new());
            }
            let nodes = derivative_nodes(degree);
            let weights = finite_difference_weights(&nodes, *order, degree)?;
            let active = self
                .options
                .numerator_sampling_scale
                .is_active_for_degree(channel_degree_by_basis[basis_position]);
            let choices = nodes
                .into_iter()
                .zip(weights)
                .filter(|(_, weight)| !weight.is_zero())
                .collect::<Vec<_>>();
            per_axis.push((basis_position, active, choices));
        }

        let mut out = Vec::new();
        let choices = per_axis
            .iter()
            .map(|(_, _, choices)| choices.clone())
            .multi_cartesian_product();
        for axis_choices in choices {
            let mut coeff = Rational::one();
            let mut offsets = BTreeMap::new();
            let mut axis_uniform = BTreeMap::new();
            for ((basis_position, use_uniform, _), (node, weight)) in
                per_axis.iter().zip(axis_choices)
            {
                coeff *= weight;
                offsets.insert(*basis_position, node);
                axis_uniform.insert(*basis_position, *use_uniform);
            }
            if coeff.is_zero() {
                continue;
            }

            let solver = EnergySolver::new(&self.signatures);
            let mut targets = vec![LinearEnergyExpr::zero(); self.signatures.len()];
            let mut labels = Vec::new();
            for (basis_position, (edge, sigma)) in basis_orig.iter().zip(cut_signs).enumerate() {
                let offset = *offsets.get(&basis_position).unwrap_or(&0);
                if *axis_uniform.get(&basis_position).unwrap_or(&false) {
                    targets[*edge] = LinearEnergyExpr::ose(EdgeIndex(*edge), i64::from(*sigma))
                        + LinearEnergyExpr::uniform_scale(i64::from(offset));
                    labels.push(format!("b{basis_position}{sigma:+}{offset:+}M"));
                } else {
                    let sample_coeff = *sigma + offset;
                    targets[*edge] =
                        LinearEnergyExpr::ose(EdgeIndex(*edge), i64::from(sample_coeff));
                    labels.push(format!("b{basis_position}{sample_coeff:+}"));
                }
            }
            let loop_exprs = solver.solve_from_target_edges(&basis_orig, &targets)?;
            let edge_exprs = solver.edge_q0_from_loop_exprs(&loop_exprs);
            let mut extra_half_edges = Vec::new();
            let mut uniform_scale_power = 0usize;
            let mut nonuniform_order_sum = 0usize;
            for (basis_position, order) in beta.iter().enumerate() {
                if *axis_uniform.get(&basis_position).unwrap_or(&false) {
                    uniform_scale_power += *order;
                    continue;
                }
                let rep_edge = self.channels[basis_logical[basis_position]].rep_edge;
                extra_half_edges.extend(std::iter::repeat_n(rep_edge, *order));
                nonuniform_order_sum += *order;
            }
            coeff *= rational_pow_i64(2, nonuniform_order_sum);
            out.push(NumeratorSample {
                coeff,
                extra_half_edges,
                loop_exprs,
                edge_exprs,
                label: format!(
                    "nD{}:fd{}",
                    beta.iter().map(usize::to_string).join(""),
                    labels.join(",")
                ),
                uniform_scale_power,
            });
        }
        Ok(out)
    }

    fn physical_sample_orbit(
        &self,
        basis_logical: &[usize],
        cut_signs: &[i32],
        alpha: &[usize],
    ) -> Result<Vec<PhysicalSample>> {
        let active_positions = alpha
            .iter()
            .enumerate()
            .filter_map(|(basis_position, order)| (*order > 0).then_some(basis_position))
            .collect::<Vec<_>>();
        let mut base_by_position = BTreeMap::new();
        let mut options_by_position = BTreeMap::new();
        for basis_position in &active_positions {
            let channel = &self.channels[basis_logical[*basis_position]];
            let sigma = cut_signs[*basis_position];
            let options = physical_group_options(channel);
            let base_tau = vec![sigma; channel.members.len()];
            let base = options
                .iter()
                .find(|option| option.tau == base_tau && option.chain_sign == sigma)
                .expect("base physical option must exist")
                .clone();
            base_by_position.insert(*basis_position, base);
            options_by_position.insert(*basis_position, options);
        }

        let mut specs = vec![base_by_position.clone()];
        for basis_position in &active_positions {
            for option in options_by_position
                .get(basis_position)
                .into_iter()
                .flatten()
            {
                let mut spec = base_by_position.clone();
                spec.insert(*basis_position, option.clone());
                specs.push(spec);
            }
        }

        let mut seen = BTreeSet::new();
        let mut samples = Vec::new();
        for spec in specs {
            let key = active_positions
                .iter()
                .filter_map(|basis_position| {
                    spec.get(basis_position)
                        .map(|option| (*basis_position, option.tau.clone(), option.chain_sign))
                })
                .collect::<Vec<_>>();
            if !seen.insert(key) {
                continue;
            }

            let basis_orig = basis_logical
                .iter()
                .map(|logical_idx| self.channels[*logical_idx].rep_edge)
                .collect::<Vec<_>>();
            let mut targets = vec![LinearEnergyExpr::zero(); self.signatures.len()];
            for (basis_position, (logical_idx, sigma)) in
                basis_logical.iter().zip(cut_signs).enumerate()
            {
                let channel = &self.channels[*logical_idx];
                let chain_sign = spec
                    .get(&basis_position)
                    .map(|option| option.chain_sign)
                    .unwrap_or(*sigma);
                targets[channel.rep_edge] =
                    LinearEnergyExpr::ose(EdgeIndex(channel.rep_edge), i64::from(chain_sign));
            }
            let solver = EnergySolver::new(&self.signatures);
            let loop_exprs = solver.solve_from_target_edges(&basis_orig, &targets)?;
            let mut edge_exprs = solver.edge_q0_from_loop_exprs(&loop_exprs);
            let mut labels = Vec::new();
            for (basis_position, logical_idx) in basis_logical.iter().enumerate() {
                let channel = &self.channels[*logical_idx];
                if let Some(option) = spec.get(&basis_position) {
                    labels.push(format!("G{}{}", channel.rep_edge, option.label));
                    for (member, sign) in channel.members.iter().zip(&option.tau) {
                        edge_exprs[*member] =
                            LinearEnergyExpr::ose(EdgeIndex(*member), i64::from(*sign));
                    }
                }
            }
            samples.push(PhysicalSample {
                loop_exprs,
                edge_exprs,
                label: if labels.is_empty() {
                    "phys0".to_string()
                } else {
                    format!("phys{}", labels.join("_"))
                },
            });
        }
        Ok(samples)
    }

    fn internal_alias(&self) -> BTreeMap<usize, usize> {
        self.channels
            .iter()
            .flat_map(|channel| {
                channel
                    .members
                    .iter()
                    .map(move |member| (*member, channel.rep_edge))
            })
            .collect()
    }

    fn basis_variable_degree_bounds(
        &self,
        bounds: &[usize],
        edge_derivs: &[Vec<Rational>],
    ) -> Vec<usize> {
        let n_basis = edge_derivs.first().map(Vec::len).unwrap_or(0);
        (0..n_basis)
            .map(|basis_position| {
                edge_derivs
                    .iter()
                    .enumerate()
                    .filter_map(|(edge_id, row)| {
                        (!row[basis_position].is_zero()).then_some(bounds[edge_id])
                    })
                    .sum()
            })
            .collect()
    }
}

struct AffineWeightSystem<'a> {
    sample_component_maps: &'a [Vec<CoeffMap>],
    target_maps: &'a [CoeffMap],
    constant_target: Rational,
}

impl<'a> AffineWeightSystem<'a> {
    fn new(
        sample_component_maps: &'a [Vec<CoeffMap>],
        target_maps: &'a [CoeffMap],
        constant_target: Rational,
    ) -> Self {
        Self {
            sample_component_maps,
            target_maps,
            constant_target,
        }
    }

    fn solve(&self) -> Result<Vec<Rational>> {
        let Some(first_sample) = self.sample_component_maps.first() else {
            return Ok(Vec::new());
        };
        let n_samples = self.sample_component_maps.len();
        let mut ordered_symbols = BTreeSet::new();
        for component_maps in self.sample_component_maps {
            for coeff_map in component_maps {
                ordered_symbols.extend(coeff_map.keys().copied());
            }
        }
        for coeff_map in self.target_maps {
            ordered_symbols.extend(coeff_map.keys().copied());
        }

        let mut rows = vec![vec![Rational::one(); n_samples]];
        let mut targets = vec![self.constant_target.clone()];
        for component_index in 0..first_sample.len() {
            for symbol in &ordered_symbols {
                let row = self
                    .sample_component_maps
                    .iter()
                    .map(|component_maps| {
                        component_maps[component_index]
                            .get(symbol)
                            .cloned()
                            .unwrap_or_else(Rational::zero)
                    })
                    .collect::<Vec<_>>();
                let target = self.target_maps[component_index]
                    .get(symbol)
                    .cloned()
                    .unwrap_or_else(Rational::zero);
                if row.iter().any(|value| !value.is_zero()) || !target.is_zero() {
                    rows.push(row);
                    targets.push(target);
                }
            }
        }
        Self::minimum_norm_affine_weights(&rows, &targets)
    }

    fn minimum_norm_affine_weights(
        rows: &[Vec<Rational>],
        targets: &[Rational],
    ) -> Result<Vec<Rational>> {
        if rows.is_empty() {
            return Ok(Vec::new());
        }
        let n_cols = rows[0].len();
        if n_cols == 0 {
            return Ok(Vec::new());
        }
        for (row, target) in rows.iter().zip(targets) {
            if row.iter().all(|value| value.is_zero()) && !target.is_zero() {
                return Err(GenerationError::NumeratorSampleOrbitIncomplete);
            }
        }

        let mut keep = Vec::new();
        let mut current = Vec::<Vec<Rational>>::new();
        let mut current_rank = 0;
        for (index, row) in rows.iter().enumerate() {
            if row.iter().all(|value| value.is_zero()) {
                continue;
            }
            let trial = current
                .iter()
                .cloned()
                .chain(std::iter::once(row.clone()))
                .collect::<Vec<_>>();
            let trial_rank = rank_rational(&trial);
            if trial_rank > current_rank {
                keep.push(index);
                current = trial;
                current_rank = trial_rank;
            }
        }
        if keep.is_empty() {
            return if targets.iter().all(|target| target.is_zero()) {
                Ok(vec![Rational::zero(); n_cols])
            } else {
                Err(GenerationError::NumeratorSampleOrbitIncomplete)
            };
        }
        let basis_rows = keep
            .iter()
            .map(|index| rows[*index].clone())
            .collect::<Vec<_>>();
        let basis_targets = keep
            .iter()
            .map(|index| targets[*index].clone())
            .collect::<Vec<_>>();
        let gram = basis_rows
            .iter()
            .map(|row_i| {
                basis_rows
                    .iter()
                    .map(|row_j| {
                        row_i
                            .iter()
                            .zip(row_j)
                            .fold(Rational::zero(), |acc, (a, b)| acc + a.clone() * b.clone())
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let dual =
            solve_rational_system(gram, basis_targets).ok_or(GenerationError::SingularBasis)?;
        let coeffs = (0..n_cols)
            .map(|column| {
                basis_rows
                    .iter()
                    .zip(&dual)
                    .fold(Rational::zero(), |acc, (row, dual_coeff)| {
                        acc + row[column].clone() * dual_coeff.clone()
                    })
            })
            .collect::<Vec<_>>();
        for (row, target) in rows.iter().zip(targets) {
            let got = row
                .iter()
                .zip(&coeffs)
                .fold(Rational::zero(), |acc, (a, b)| acc + a.clone() * b.clone());
            if got != target.clone() {
                return Err(GenerationError::NumeratorSampleOrbitIncomplete);
            }
        }
        Ok(coeffs)
    }
}

fn physical_group_options(channel: &LogicalChannel) -> Vec<PhysicalGroupOption> {
    (0..channel.power)
        .map(|_| [-1, 1])
        .multi_cartesian_product()
        .flat_map(|tau| {
            let chain_signs = if tau.iter().all(|sign| *sign > 0) {
                vec![1]
            } else if tau.iter().all(|sign| *sign < 0) {
                vec![-1]
            } else {
                vec![-1, 1]
            };
            let tau_label = tau
                .iter()
                .map(|sign| if *sign > 0 { "p" } else { "m" })
                .join("");
            chain_signs
                .into_iter()
                .map(move |chain_sign| PhysicalGroupOption {
                    tau: tau.clone(),
                    chain_sign,
                    label: format!("t{}c{}", tau_label, if chain_sign > 0 { "p" } else { "m" }),
                })
        })
        .collect()
}

fn component_coeff_maps(
    loop_exprs: &[LinearEnergyExpr],
    edge_exprs: &[LinearEnergyExpr],
    internal_alias: &BTreeMap<usize, usize>,
) -> Vec<CoeffMap> {
    loop_exprs
        .iter()
        .chain(edge_exprs)
        .map(|expr| coeff_map(expr, internal_alias))
        .collect()
}

fn coeff_map(expr: &LinearEnergyExpr, internal_alias: &BTreeMap<usize, usize>) -> CoeffMap {
    let mut out = BTreeMap::new();
    let constant = expr.constant.rational_coeff();
    if !constant.is_zero() {
        out.insert(CoeffSymbol::Constant, constant);
    }
    for (edge_id, coeff) in &expr.internal_terms {
        let aliased = internal_alias.get(&edge_id.0).copied().unwrap_or(edge_id.0);
        add_coeff(
            &mut out,
            CoeffSymbol::Internal(aliased),
            coeff.rational_coeff(),
        );
    }
    for (edge_id, coeff) in &expr.external_terms {
        add_coeff(
            &mut out,
            CoeffSymbol::External(edge_id.0),
            coeff.rational_coeff(),
        );
    }
    if !expr.uniform_scale_coeff.is_zero_coeff() {
        add_coeff(
            &mut out,
            CoeffSymbol::UniformScale,
            expr.uniform_scale_coeff.rational_coeff(),
        );
    }
    out.retain(|_, value| !value.is_zero());
    out
}

fn add_coeff(map: &mut CoeffMap, key: CoeffSymbol, value: Rational) {
    let entry = map.entry(key).or_insert(Rational::zero());
    *entry = entry.clone() + value;
}

fn orientation_from_edge_exprs(edge_exprs: &[LinearEnergyExpr]) -> (EdgeVec<Orientation>, String) {
    let mut label = String::new();
    let orientation = EdgeVec::from_iter(edge_exprs.iter().enumerate().map(|(edge, expr)| {
        if *expr == LinearEnergyExpr::ose(EdgeIndex(edge), 1) {
            label.push('+');
            Orientation::Default
        } else if *expr == LinearEnergyExpr::ose(EdgeIndex(edge), -1) {
            label.push('-');
            Orientation::Reversed
        } else if *expr == LinearEnergyExpr::zero() {
            label.push('0');
            Orientation::Undirected
        } else {
            label.push('x');
            Orientation::Undirected
        }
    }));
    (orientation, label)
}

fn rational_to_coefficient(value: Rational) -> Result<symbolica::atom::Atom> {
    Ok(rational_coeff_atom(value))
}

fn derivative_nodes(degree: usize) -> Vec<i32> {
    let mut nodes = vec![0];
    let mut step = 1;
    while nodes.len() < degree + 1 {
        nodes.push(step);
        if nodes.len() == degree + 1 {
            break;
        }
        nodes.push(-step);
        step += 1;
    }
    nodes
}

fn finite_difference_weights(
    nodes: &[i32],
    derivative_order: usize,
    degree: usize,
) -> Result<Vec<Rational>> {
    if derivative_order > degree {
        return Ok(vec![Rational::zero(); nodes.len()]);
    }
    let matrix = (0..=degree)
        .map(|power| {
            nodes
                .iter()
                .map(|node| Rational::from(*node).pow_usize(power))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let rhs = (0..=degree)
        .map(|power| {
            if power == derivative_order {
                factorial(derivative_order)
            } else {
                Rational::zero()
            }
        })
        .collect::<Vec<_>>();
    solve_rational_system(matrix, rhs).ok_or(GenerationError::SingularBasis)
}

fn intern_linear_surface(
    expression: &mut ThreeDExpression<OrientationID>,
    surface_index: &mut HashMap<(LinearSurfaceKind, LinearEnergyExpr), HybridSurfaceID>,
    surface_expr: LinearEnergyExpr,
    numerator_only: bool,
) -> HybridSurfaceID {
    let surface_expr = surface_expr.canonical();
    let kind = classify_surface_kind(&surface_expr);
    let key = (kind, surface_expr.clone());
    if let Some(surface_id) = surface_index.get(&key) {
        return *surface_id;
    }

    let id = crate::surface::LinearSurfaceID(expression.surfaces.linear_surface_cache.len());
    let surface_id = HybridSurfaceID::Linear(id);
    expression
        .surfaces
        .linear_surface_cache
        .push(LinearSurface {
            kind,
            expression: surface_expr,
            origin: SurfaceOrigin::Physical,
            numerator_only,
        });
    surface_index.insert(key, surface_id);
    surface_id
}

fn classify_surface_kind(expr: &LinearEnergyExpr) -> LinearSurfaceKind {
    let coeffs = expr
        .internal_terms
        .iter()
        .filter_map(|(_, coeff)| (!coeff.is_zero_coeff()).then_some(coeff))
        .collect::<Vec<_>>();
    if coeffs.len() <= 1
        || coeffs.iter().all(|coeff| !coeff.is_negative_coeff())
        || coeffs.iter().all(|coeff| coeff.is_negative_coeff())
    {
        LinearSurfaceKind::Esurface
    } else {
        LinearSurfaceKind::Hsurface
    }
}

fn denominator_tree_from_chain(chain: &[HybridSurfaceID]) -> Tree<HybridSurfaceID> {
    if chain.is_empty() {
        return Tree::from_root(HybridSurfaceID::Unit);
    }
    let mut tree = Tree::from_root(chain[0]);
    let mut parent = NodeId::root();
    for surface_id in chain.iter().copied().skip(1) {
        tree.insert_node(parent, surface_id);
        parent = NodeId(parent.0 + 1);
    }
    tree
}

fn denominator_tree_from_chains(chains: &[Vec<HybridSurfaceID>]) -> Tree<HybridSurfaceID> {
    if chains.is_empty() || chains.iter().all(Vec::is_empty) {
        return Tree::from_root(HybridSurfaceID::Unit);
    }

    let mut tree = Tree::from_root(HybridSurfaceID::Unit);
    for chain in chains {
        if chain.is_empty() {
            tree.insert_node(NodeId::root(), HybridSurfaceID::Unit);
            continue;
        }
        let mut parent = NodeId::root();
        for surface_id in chain {
            let existing_child = tree
                .get_node(parent)
                .children
                .iter()
                .copied()
                .find(|child| tree.get_node(*child).data == *surface_id);
            if let Some(child) = existing_child {
                parent = child;
                continue;
            }
            let child = NodeId(tree.get_num_nodes());
            tree.insert_node(parent, *surface_id);
            parent = child;
        }
    }
    tree
}

fn denominator_tree_chains(tree: &Tree<HybridSurfaceID>) -> Vec<Vec<HybridSurfaceID>> {
    fn walk(
        tree: &Tree<HybridSurfaceID>,
        node_id: NodeId,
        current: &mut Vec<HybridSurfaceID>,
        out: &mut Vec<Vec<HybridSurfaceID>>,
    ) {
        let node = tree.get_node(node_id);
        if node.data != HybridSurfaceID::Unit {
            current.push(node.data);
        }
        if node.children.is_empty() {
            out.push(current.clone());
        } else {
            for child in &node.children {
                walk(tree, *child, current, out);
            }
        }
        if node.data != HybridSurfaceID::Unit {
            current.pop();
        }
    }

    if tree.get_num_nodes() == 0 {
        return Vec::new();
    }
    let mut out = Vec::new();
    walk(tree, NodeId::root(), &mut Vec::new(), &mut out);
    if out.is_empty() {
        vec![Vec::new()]
    } else {
        out
    }
}

fn solve_loop_energy_substitutions(
    parsed: &ParsedGraph,
    signatures: &[MomentumSignature],
    basis: &[usize],
    cut_signs: &[i32],
) -> Result<(Vec<LinearEnergyExpr>, Vec<LinearEnergyExpr>)> {
    let mut target_edge_exprs = vec![LinearEnergyExpr::zero(); signatures.len()];
    for (edge_index, cut_sign) in basis.iter().zip(cut_signs) {
        target_edge_exprs[*edge_index] =
            LinearEnergyExpr::ose(EdgeIndex(*edge_index), i64::from(*cut_sign));
    }
    let loop_exprs =
        solve_loop_energy_from_target_edge_exprs(signatures, basis, &target_edge_exprs)?;
    let mut edge_exprs = edge_q0_from_loop_exprs(signatures, &loop_exprs);
    apply_initial_state_cut_edge_energy_exprs(parsed, &mut edge_exprs);
    Ok((loop_exprs, edge_exprs))
}

fn solve_loop_energy_particular_from_target_edge_exprs(
    signatures: &[MomentumSignature],
    basis: &[usize],
    target_edge_exprs: &[LinearEnergyExpr],
) -> Result<Vec<LinearEnergyExpr>> {
    let n_loops = signatures
        .first()
        .map(|signature| signature.loop_signature.len())
        .unwrap_or(0);
    if basis.len() == n_loops {
        return solve_loop_energy_from_target_edge_exprs(signatures, basis, target_edge_exprs);
    }
    if basis.is_empty() {
        return Ok(vec![LinearEnergyExpr::zero(); n_loops]);
    }

    let matrix = basis
        .iter()
        .map(|edge_index| signatures[*edge_index].loop_signature.clone())
        .collect::<Vec<_>>();
    if rank_i64(
        &matrix
            .iter()
            .map(|row| row.iter().map(|value| i64::from(*value)).collect())
            .collect::<Vec<Vec<_>>>(),
    ) != basis.len()
    {
        return Err(GenerationError::SingularBasis);
    }
    let rhs = basis
        .iter()
        .map(|edge_index| {
            let mut expr = target_edge_exprs[*edge_index].clone();
            for (external_id, coeff) in signatures[*edge_index]
                .external_signature
                .iter()
                .enumerate()
            {
                if *coeff != 0 {
                    expr = expr
                        - LinearEnergyExpr::external(EdgeIndex(external_id), i64::from(*coeff));
                }
            }
            expr
        })
        .collect::<Vec<_>>();

    for columns in (0..n_loops).combinations(basis.len()) {
        let square = matrix
            .iter()
            .map(|row| {
                columns
                    .iter()
                    .map(|column| row[*column])
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        if !determinant_i32_is_nonzero(&square) {
            continue;
        }
        let solved = solve_expr_system_unimodular(
            square
                .into_iter()
                .map(|row| row.into_iter().map(i64::from).collect())
                .collect(),
            rhs.clone(),
        )?;
        let mut out = vec![LinearEnergyExpr::zero(); n_loops];
        for (column, expr) in columns.into_iter().zip(solved) {
            out[column] = expr;
        }
        return Ok(out);
    }
    Err(GenerationError::SingularBasis)
}

fn choose_basis_indices_from_edges(
    signatures: &[MomentumSignature],
    edge_ids: &[usize],
) -> Result<Vec<usize>> {
    let n_loops = signatures
        .first()
        .map(|signature| signature.loop_signature.len())
        .unwrap_or(0);
    if n_loops == 0 {
        return Ok(Vec::new());
    }
    edge_ids
        .iter()
        .copied()
        .combinations(n_loops)
        .find(|basis| {
            let matrix = basis
                .iter()
                .map(|edge_index| signatures[*edge_index].loop_signature.clone())
                .collect::<Vec<_>>();
            determinant_i32_is_nonzero(&matrix)
        })
        .ok_or(GenerationError::SingularBasis)
}

fn solve_loop_energy_from_target_edge_exprs(
    signatures: &[MomentumSignature],
    basis: &[usize],
    target_edge_exprs: &[LinearEnergyExpr],
) -> Result<Vec<LinearEnergyExpr>> {
    let n_loops = signatures
        .first()
        .map(|signature| signature.loop_signature.len())
        .unwrap_or(0);
    if basis.len() != n_loops {
        return Err(GenerationError::SingularBasis);
    }

    let matrix = basis
        .iter()
        .map(|edge_index| {
            signatures[*edge_index]
                .loop_signature
                .iter()
                .map(|value| i64::from(*value))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let rhs = basis
        .iter()
        .map(|edge_index| {
            let mut expr = target_edge_exprs[*edge_index].clone();
            for (external_id, coeff) in signatures[*edge_index]
                .external_signature
                .iter()
                .enumerate()
            {
                if *coeff != 0 {
                    expr = expr
                        - LinearEnergyExpr::external(EdgeIndex(external_id), i64::from(*coeff));
                }
            }
            expr
        })
        .collect::<Vec<_>>();

    solve_expr_system_unimodular(matrix, rhs)
}

fn edge_q0_from_loop_exprs(
    signatures: &[MomentumSignature],
    loop_exprs: &[LinearEnergyExpr],
) -> Vec<LinearEnergyExpr> {
    signatures
        .iter()
        .map(|signature| {
            let mut expr = LinearEnergyExpr::zero();
            for (loop_id, coeff) in signature.loop_signature.iter().enumerate() {
                if *coeff != 0 {
                    expr = expr + loop_exprs[loop_id].clone() * i64::from(*coeff);
                }
            }
            for (external_id, coeff) in signature.external_signature.iter().enumerate() {
                if *coeff != 0 {
                    expr = expr
                        + LinearEnergyExpr::external(EdgeIndex(external_id), i64::from(*coeff));
                }
            }
            expr.canonical()
        })
        .collect()
}

fn apply_initial_state_cut_edge_energy_exprs(
    parsed: &ParsedGraph,
    edge_exprs: &mut [LinearEnergyExpr],
) {
    for cut_edge in &parsed.initial_state_cut_edges {
        if let Some(edge_expr) = edge_exprs.get_mut(cut_edge.edge_id) {
            *edge_expr = LinearEnergyExpr::external(
                EdgeIndex(cut_edge.external_id),
                i64::from(cut_edge.external_sign),
            );
        }
    }
}

fn cff_duplicate_signature_excess(parsed: &ParsedGraph) -> usize {
    let mut counts = BTreeMap::<MomentumSignature, usize>::new();
    for edge in &parsed.internal_edges {
        if parsed.is_initial_state_cut_edge(edge.edge_id) {
            continue;
        }
        let (signature, _) = edge.signature.canonical_up_to_sign();
        *counts.entry(signature).or_default() += 1;
    }
    counts.values().map(|count| count.saturating_sub(1)).sum()
}

fn remap_initial_state_cut_edges(
    parsed: &ParsedGraph,
    local_to_orig: &[usize],
) -> Vec<ParsedGraphInitialStateCutEdge> {
    let orig_to_local = local_to_orig
        .iter()
        .enumerate()
        .map(|(local_id, orig_id)| (*orig_id, local_id))
        .collect::<BTreeMap<_, _>>();
    parsed
        .initial_state_cut_edges
        .iter()
        .filter_map(|cut_edge| {
            orig_to_local
                .get(&cut_edge.edge_id)
                .copied()
                .map(|edge_id| ParsedGraphInitialStateCutEdge {
                    edge_id,
                    external_id: cut_edge.external_id,
                    external_sign: cut_edge.external_sign,
                })
        })
        .collect()
}

fn solve_expr_system_unimodular(
    matrix: Vec<Vec<i64>>,
    rhs: Vec<LinearEnergyExpr>,
) -> Result<Vec<LinearEnergyExpr>> {
    let n = matrix.len();
    if rhs.len() != n || matrix.iter().any(|row| row.len() != n) {
        return Err(GenerationError::SingularBasis);
    }
    let rational_matrix = matrix
        .iter()
        .map(|row| row.iter().map(|value| Rational::from(*value)).collect())
        .collect::<Vec<Vec<_>>>();
    let mut out = vec![LinearEnergyExpr::zero(); n];
    for (rhs_position, rhs_expr) in rhs.iter().enumerate().take(n) {
        let unit = (0..n)
            .map(|row| {
                if row == rhs_position {
                    Rational::one()
                } else {
                    Rational::zero()
                }
            })
            .collect::<Vec<_>>();
        let solution = solve_rational_system(rational_matrix.clone(), unit)
            .ok_or(GenerationError::SingularBasis)?;
        for (solution_position, coeff) in solution.into_iter().enumerate() {
            if coeff.is_zero() {
                continue;
            }
            out[solution_position] = out[solution_position].clone()
                + scale_linear_energy_expr_rational(rhs_expr, &coeff)?;
        }
    }
    Ok(out.into_iter().map(LinearEnergyExpr::canonical).collect())
}

fn scale_linear_energy_expr_rational(
    expr: &LinearEnergyExpr,
    scale: &Rational,
) -> Result<LinearEnergyExpr> {
    fn scale_terms(
        terms: &[(EdgeIndex, symbolica::atom::Atom)],
        scale: &Rational,
    ) -> Vec<(EdgeIndex, symbolica::atom::Atom)> {
        terms
            .iter()
            .filter_map(|(edge_id, coeff)| {
                let value = coeff.rational_coeff() * scale.clone();
                (!value.is_zero()).then_some((*edge_id, rational_coeff_atom(value)))
            })
            .collect()
    }

    let uniform = rational_coeff_atom(expr.uniform_scale_coeff.rational_coeff() * scale.clone());

    Ok(LinearEnergyExpr {
        internal_terms: scale_terms(&expr.internal_terms, scale),
        external_terms: scale_terms(&expr.external_terms, scale),
        uniform_scale_coeff: uniform,
        constant: rational_coeff_atom(expr.constant.rational_coeff() * scale.clone()),
    }
    .canonical())
}

fn has_higher_energy_power(options: &Generate3DExpressionOptions) -> bool {
    options
        .energy_degree_bounds
        .iter()
        .any(|(_, degree)| *degree > 1)
}

#[cfg(all(test, feature = "old_cff"))]
mod tests {
    use super::*;

    #[test]
    fn old_cff_rejects_higher_energy_powers() {
        let parsed = crate::graph_io::test_graphs::box_graph();
        let options = Generate3DExpressionOptions {
            representation: RepresentationMode::OldCff,
            energy_degree_bounds: vec![(0, 2)],
            ..Default::default()
        };

        assert!(matches!(
            generate_3d_expression(&parsed, &options),
            Err(GenerationError::OldCffHigherEnergyPower)
        ));
    }

    #[test]
    fn old_cff_parsed_mode_uses_migrated_pure_cff_generator() {
        let parsed = crate::graph_io::test_graphs::box_graph();
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::OldCff,
                energy_degree_bounds: Vec::new(),
                ..Default::default()
            },
        )
        .unwrap();

        assert!(!expression.orientations.is_empty());
    }
}

#[cfg(test)]
mod graph_source_tests {
    use std::collections::{BTreeMap, BTreeSet};

    use crate::graph_io::EnergyEdgeIndexMap;

    use super::*;

    #[derive(Clone)]
    struct RemappedSource {
        parsed: ParsedGraph,
        edge_map: EnergyEdgeIndexMap,
    }

    impl ThreeDGraphSource for RemappedSource {
        fn to_three_d_parsed_graph(&self) -> crate::graph_io::Result<ParsedGraph> {
            Ok(self.parsed.clone())
        }

        fn energy_edge_index_map(&self, _parsed: &ParsedGraph) -> Option<EnergyEdgeIndexMap> {
            Some(self.edge_map.clone())
        }
    }

    #[test]
    fn rich_graph_source_remaps_compact_energy_ids_back_to_source_edge_ids() {
        let parsed = crate::graph_io::test_graphs::box_graph();
        let source = RemappedSource {
            parsed,
            edge_map: EnergyEdgeIndexMap {
                internal: BTreeMap::from([(0, 10), (1, 11), (2, 12), (3, 13)]),
                external: BTreeMap::from([(0, 20), (1, 21), (2, 22)]),
                orientation_edge_count: 16,
            },
        };

        let expression = generate_3d_expression(
            &source,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Cff,
                energy_degree_bounds: vec![(10, 1), (11, 1)],
                ..Default::default()
            },
        )
        .unwrap();

        let half_edges = expression
            .orientations
            .iter()
            .flat_map(|orientation| &orientation.variants)
            .flat_map(|variant| variant.half_edges.iter().map(|edge_id| edge_id.0))
            .collect::<BTreeSet<_>>();
        assert_eq!(half_edges, BTreeSet::from([10, 11, 12, 13]));

        assert!(
            expression
                .orientations
                .iter()
                .all(|orientation| orientation.edge_energy_map.len() == 16)
        );
        assert!(
            expression
                .orientations
                .iter()
                .flat_map(|orientation| orientation.edge_energy_map.iter().enumerate())
                .all(|(edge_id, expr)| expr.is_zero() || edge_id >= 10)
        );

        assert!(
            expression
                .surfaces
                .linear_surface_cache
                .iter()
                .flat_map(|surface| &surface.expression.external_terms)
                .all(|(edge_id, _)| edge_id.0 >= 20)
        );
    }
}

#[cfg(test)]
mod initial_state_cut_tests {
    use super::*;

    fn generated_expression(
        representation: RepresentationMode,
        external_sign: i32,
    ) -> ThreeDExpression<OrientationID> {
        let parsed = crate::graph_io::test_graphs::initial_state_cut_line_graph(external_sign);
        generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                representation,
                energy_degree_bounds: Vec::new(),
                ..Default::default()
            },
        )
        .unwrap()
    }

    fn external_coeff(expr: &LinearEnergyExpr, external_id: usize) -> Rational {
        expr.external_terms
            .iter()
            .filter(|(edge_id, _)| edge_id.0 == external_id)
            .map(|(_, coeff)| coeff.rational_coeff())
            .fold(Rational::zero(), |acc, coeff| acc + coeff)
    }

    #[test]
    fn cff_initial_state_cut_edges_are_external_energy_shifts() {
        for external_sign in [1, -1] {
            let expression = generated_expression(RepresentationMode::Cff, external_sign);
            assert!(!expression.orientations.is_empty());
            assert!(expression.orientations.iter().all(|orientation| {
                orientation.edge_energy_map.first().is_some_and(|expr| {
                    expr.internal_terms.is_empty() && external_coeff(expr, 0) == external_sign
                })
            }));
            assert!(expression.orientations.iter().all(|orientation| {
                orientation.variants.iter().all(|variant| {
                    variant.half_edges == vec![EdgeIndex(1)]
                        && denominator_tree_chains(&variant.denominator)
                            .into_iter()
                            .flatten()
                            .all(|surface_id| match surface_id {
                                HybridSurfaceID::Linear(id) => {
                                    expression.surfaces.linear_surface_cache[id]
                                        .expression
                                        .internal_terms
                                        .iter()
                                        .all(|(edge_id, _)| edge_id.0 != 0)
                                }
                                _ => true,
                            })
                })
            }));

            let default_loop_orientation = expression
                .orientations
                .iter()
                .find(|orientation| {
                    orientation.data.orientation[EdgeIndex(1)] == Orientation::Default
                })
                .expect("default loop orientation should be present");
            let coeffs = default_loop_orientation
                .variants
                .iter()
                .flat_map(|variant| denominator_tree_chains(&variant.denominator))
                .flatten()
                .filter_map(|surface_id| match surface_id {
                    HybridSurfaceID::Linear(id) => Some(external_coeff(
                        &expression.surfaces.linear_surface_cache[id].expression,
                        0,
                    )),
                    _ => None,
                })
                .collect::<Vec<_>>();
            assert!(
                coeffs.contains(&Rational::from(external_sign)),
                "CFF should preserve the sign of the initial-state cut edge as an external shift"
            );
        }
    }

    #[test]
    fn ltd_initial_state_cut_edges_are_not_cut_denominators() {
        let expression = generated_expression(RepresentationMode::Ltd, -1);
        assert!(!expression.orientations.is_empty());
        assert!(expression.orientations.iter().all(|orientation| {
            orientation
                .edge_energy_map
                .first()
                .is_some_and(|expr| expr.internal_terms.is_empty() && external_coeff(expr, 0) == -1)
                && orientation.variants.iter().all(|variant| {
                    variant.half_edges == vec![EdgeIndex(1)]
                        && denominator_tree_chains(&variant.denominator)
                            .into_iter()
                            .flatten()
                            .all(|surface_id| match surface_id {
                                HybridSurfaceID::Linear(id) => {
                                    expression.surfaces.linear_surface_cache[id]
                                        .expression
                                        .internal_terms
                                        .iter()
                                        .all(|(edge_id, _)| edge_id.0 != 0)
                                }
                                _ => true,
                            })
                })
        }));
    }

    #[test]
    fn cff_duplicate_sign_correction_is_mass_insensitive() {
        let mut parsed = crate::graph_io::test_graphs::box_pow3_graph();
        parsed.internal_edges[4].mass_key = Some("different_mass".to_string());
        assert_eq!(cff_duplicate_signature_excess(&parsed), 2);
    }
}

#[cfg(test)]
mod ltd_tests {
    use super::*;

    fn parsed_fixture(name: &str) -> ParsedGraph {
        match name {
            "box.dot" => crate::graph_io::test_graphs::box_graph(),
            "box_pow3.dot" => crate::graph_io::test_graphs::box_pow3_graph(),
            "sunrise_pow4.dot" => crate::graph_io::test_graphs::sunrise_pow4_graph(),
            other => panic!("unknown parsed fixture {other}"),
        }
    }

    #[test]
    fn ltd_generation_builds_normal_box_structure() {
        let parsed = parsed_fixture("box.dot");
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Ltd,
                energy_degree_bounds: Vec::new(),
                ..Default::default()
            },
        )
        .unwrap();

        assert!(!expression.orientations.is_empty());
        assert!(!expression.surfaces.linear_surface_cache.is_empty());
        assert!(
            expression
                .orientations
                .iter()
                .all(|orientation| !orientation.variants.is_empty())
        );
    }

    #[test]
    fn ltd_generation_builds_repeated_propagator_confluent_structure() {
        let parsed = parsed_fixture("box_pow3.dot");
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Ltd,
                energy_degree_bounds: Vec::new(),
                ..Default::default()
            },
        )
        .unwrap();

        assert!(!expression.orientations.is_empty());
        assert!(expression.orientations.iter().all(|orientation| {
            orientation
                .data
                .label
                .as_deref()
                .is_some_and(|label| label.contains("|N"))
        }));
        assert!(
            expression
                .orientations
                .iter()
                .any(|orientation| orientation.variants.len() > 1)
        );
    }

    #[test]
    fn ltd_generation_builds_multiloop_repeated_confluent_structure() {
        let parsed = parsed_fixture("sunrise_pow4.dot");
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Ltd,
                energy_degree_bounds: Vec::new(),
                ..Default::default()
            },
        )
        .unwrap();

        assert!(!expression.orientations.is_empty());
        assert!(
            expression
                .orientations
                .iter()
                .any(|orientation| orientation.variants.len() > 1)
        );
        assert!(
            expression
                .orientations
                .iter()
                .all(|orientation| orientation.data.numerator_map_index.is_some())
        );
    }

    #[test]
    fn bounded_repeated_ltd_can_use_uniform_sampling_scale() {
        let parsed = parsed_fixture("box_pow3.dot");
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Ltd,
                energy_degree_bounds: vec![(3, 4)],
                numerator_sampling_scale: NumeratorSamplingScaleMode::BeyondQuadratic,
                preserve_internal_edges_as_four_d_denominators: Vec::new(),
            },
        )
        .unwrap();

        assert!(expression.orientations.iter().any(|orientation| {
            orientation
                .edge_energy_map
                .iter()
                .any(LinearEnergyExpr::uses_uniform_scale)
        }));
        assert!(
            expression
                .orientations
                .iter()
                .flat_map(|orientation| &orientation.variants)
                .any(|variant| variant.uniform_scale_power > 0)
        );
    }

    #[test]
    fn cff_generation_builds_repeated_box_with_branching_denominator_trees() {
        let parsed = parsed_fixture("box_pow3.dot");
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Cff,
                energy_degree_bounds: Vec::new(),
                ..Default::default()
            },
        )
        .unwrap();

        assert!(!expression.orientations.is_empty());
        assert!(
            expression
                .orientations
                .iter()
                .all(|orientation| orientation.variants.len() == 1)
        );
        assert!(
            expression
                .orientations
                .iter()
                .flat_map(|orientation| &orientation.variants)
                .any(|variant| variant.denominator.get_bottom_layer().len() > 1)
        );
        assert!(
            expression
                .surfaces
                .linear_surface_cache
                .iter()
                .any(|surface| !surface.expression.external_terms.is_empty())
        );
    }

    #[test]
    fn cff_generation_preserves_external_tree_edges() {
        let parsed = crate::graph_io::test_graphs::triangle_with_external_tree_graph();
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Cff,
                preserve_internal_edges_as_four_d_denominators: vec![3],
                ..Default::default()
            },
        )
        .unwrap();

        assert!(!expression.orientations.is_empty());
        assert!(expression.orientations.iter().all(|orientation| {
            orientation.data.orientation[EdgeIndex(3)] == Orientation::Undirected
                && orientation.edge_energy_map.len() == 4
                && orientation.edge_energy_map[3] == LinearEnergyExpr::external(EdgeIndex(1), 1)
        }));
        assert_eq!(expression.residual_denominators.len(), 1);
        assert_eq!(expression.residual_denominators[0].edge_id, EdgeIndex(3));
        assert_eq!(expression.residual_denominators[0].power, 1);
        assert!(
            expression
                .orientations
                .iter()
                .flat_map(|orientation| &orientation.edge_energy_map[3].internal_terms)
                .next()
                .is_none()
        );
        assert!(
            expression
                .orientations
                .iter()
                .flat_map(|orientation| &orientation.variants)
                .flat_map(|variant| &variant.half_edges)
                .all(|edge_id| edge_id.0 != 3)
        );
        assert!(
            expression
                .surfaces
                .linear_surface_cache
                .iter()
                .flat_map(|surface| &surface.expression.internal_terms)
                .all(|(edge_id, _)| edge_id.0 != 3)
        );
    }

    #[test]
    fn cff_generation_leaves_pure_trees_untouched() {
        let parsed = crate::graph_io::test_graphs::pure_tree_graph();
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Cff,
                preserve_internal_edges_as_four_d_denominators: vec![0, 1],
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(expression.orientations.len(), 1);
        assert!(expression.surfaces.linear_surface_cache.is_empty());
        assert_eq!(
            expression
                .residual_denominators
                .iter()
                .map(|denominator| denominator.edge_id)
                .collect_vec(),
            vec![EdgeIndex(0), EdgeIndex(1)]
        );
        let orientation = &expression.orientations[OrientationID(0)];
        assert!(
            orientation
                .data
                .orientation
                .iter()
                .all(|(_, value)| *value == Orientation::Undirected)
        );
        assert_eq!(
            orientation.edge_energy_map,
            vec![
                LinearEnergyExpr::external(EdgeIndex(0), 1),
                LinearEnergyExpr::external(EdgeIndex(0), 1)
                    + LinearEnergyExpr::external(EdgeIndex(1), 1),
            ]
        );
        assert_eq!(orientation.variants.len(), 1);
        assert!(orientation.variants[0].half_edges.is_empty());
        assert_eq!(
            orientation.variants[0]
                .denominator
                .get_node(NodeId::root())
                .data,
            HybridSurfaceID::Unit
        );
    }

    #[test]
    fn cff_generation_builds_multiloop_repeated_high_power_known_factor_completion() {
        let parsed = parsed_fixture("sunrise_pow4.dot");
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Cff,
                energy_degree_bounds: vec![(2, 5)],
                ..Default::default()
            },
        )
        .unwrap();

        assert!(!expression.orientations.is_empty());
        assert!(
            expression
                .surfaces
                .linear_surface_cache
                .iter()
                .all(|surface| surface.kind == LinearSurfaceKind::Esurface)
        );
        assert!(
            expression
                .orientations
                .iter()
                .flat_map(|orientation| &orientation.variants)
                .any(|variant| {
                    variant.origin.as_deref() == Some("bounded_degree_known_factor_cff")
                        && !variant.numerator_surfaces.is_empty()
                })
        );
    }

    #[test]
    fn cff_generation_builds_one_loop_quadratic_e_surface_completion() {
        let parsed = parsed_fixture("box.dot");
        let ordinary = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Cff,
                energy_degree_bounds: Vec::new(),
                ..Default::default()
            },
        )
        .unwrap();
        let bounded = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Cff,
                energy_degree_bounds: vec![(0, 2)],
                ..Default::default()
            },
        )
        .unwrap();

        assert!(bounded.orientations.len() > ordinary.orientations.len());
        assert!(
            bounded
                .surfaces
                .linear_surface_cache
                .iter()
                .all(|surface| surface.kind == LinearSurfaceKind::Esurface)
        );
        assert!(
            bounded
                .orientations
                .iter()
                .flat_map(|orientation| &orientation.variants)
                .any(|variant| variant.origin.as_deref().is_some_and(|origin| {
                    origin.starts_with("bounded_degree_e_surface_pinch_cff")
                }))
        );
    }

    #[test]
    fn cff_generation_builds_one_loop_single_high_power_completion() {
        let parsed = parsed_fixture("box.dot");
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Cff,
                energy_degree_bounds: vec![(0, 3)],
                ..Default::default()
            },
        )
        .unwrap();

        assert!(!expression.orientations.is_empty());
        assert!(
            expression
                .surfaces
                .linear_surface_cache
                .iter()
                .all(|surface| surface.kind == LinearSurfaceKind::Esurface)
        );
        assert!(
            expression
                .orientations
                .iter()
                .flat_map(|orientation| &orientation.variants)
                .any(|variant| !variant.numerator_surfaces.is_empty())
        );
    }

    #[test]
    fn cff_generation_builds_one_loop_mixed_high_power_known_factor_completion() {
        let parsed = parsed_fixture("box.dot");
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Cff,
                energy_degree_bounds: vec![(0, 1), (1, 1), (2, 1), (3, 3)],
                ..Default::default()
            },
        )
        .unwrap();

        assert!(!expression.orientations.is_empty());
        assert!(
            expression
                .surfaces
                .linear_surface_cache
                .iter()
                .all(|surface| surface.kind == LinearSurfaceKind::Esurface)
        );
        assert!(
            expression
                .orientations
                .iter()
                .flat_map(|orientation| &orientation.variants)
                .any(|variant| {
                    variant.origin.as_deref() == Some("bounded_degree_known_factor_cff")
                        && !variant.numerator_surfaces.is_empty()
                })
        );
    }

    #[test]
    fn cff_generation_builds_one_loop_quartic_known_factor_contact_completion() {
        let parsed = parsed_fixture("box.dot");
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Cff,
                energy_degree_bounds: vec![(0, 4), (1, 2)],
                ..Default::default()
            },
        )
        .unwrap();

        assert!(!expression.orientations.is_empty());
        assert!(
            expression
                .surfaces
                .linear_surface_cache
                .iter()
                .all(|surface| surface.kind == LinearSurfaceKind::Esurface)
        );
        assert!(
            expression
                .orientations
                .iter()
                .flat_map(|orientation| &orientation.variants)
                .any(|variant| {
                    variant.denominator.get_node(NodeId::root()).data == HybridSurfaceID::Unit
                })
        );
    }

    #[test]
    fn cff_generation_builds_one_loop_repeated_quadratic_recursive_completion() {
        let parsed = parsed_fixture("box_pow3.dot");
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Cff,
                energy_degree_bounds: vec![(0, 2), (1, 1), (2, 1), (3, 2)],
                ..Default::default()
            },
        )
        .unwrap();

        assert!(!expression.orientations.is_empty());
        assert!(
            expression
                .surfaces
                .linear_surface_cache
                .iter()
                .all(|surface| surface.kind == LinearSurfaceKind::Esurface)
        );
        assert!(
            expression
                .orientations
                .iter()
                .flat_map(|orientation| &orientation.variants)
                .any(|variant| variant.origin.as_deref().is_some_and(|origin| {
                    origin.starts_with("bounded_degree_quadratic_recursive_contact")
                }))
        );
        assert!(
            expression
                .orientations
                .iter()
                .flat_map(|orientation| &orientation.variants)
                .any(|variant| variant.origin.as_deref().is_some_and(|origin| {
                    origin.starts_with("bounded_degree_quadratic_recursive_remainder")
                }))
        );
    }

    #[test]
    fn cff_generation_builds_multiloop_single_quadratic_lower_sector_completion() {
        let parsed = parsed_fixture("sunrise_pow4.dot");
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Cff,
                energy_degree_bounds: vec![(0, 2)],
                ..Default::default()
            },
        )
        .unwrap();

        let denominator_surface_ids = expression
            .orientations
            .iter()
            .flat_map(|orientation| &orientation.variants)
            .flat_map(|variant| denominator_tree_chains(&variant.denominator))
            .flatten()
            .filter_map(|surface_id| match surface_id {
                HybridSurfaceID::Linear(id) => Some(id.0),
                _ => None,
            })
            .collect::<BTreeSet<_>>();
        assert!(!denominator_surface_ids.is_empty());
        assert!(denominator_surface_ids.into_iter().all(|surface_id| {
            expression.surfaces.linear_surface_cache[LinearSurfaceID(surface_id)].kind
                == LinearSurfaceKind::Esurface
        }));
    }

    #[test]
    fn cff_generation_builds_one_loop_repeated_high_power_channel_completion() {
        let parsed = parsed_fixture("box_pow3.dot");
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Cff,
                energy_degree_bounds: vec![(0, 1), (1, 1), (3, 4)],
                ..Default::default()
            },
        )
        .unwrap();

        assert!(!expression.orientations.is_empty());
        assert!(
            expression
                .surfaces
                .linear_surface_cache
                .iter()
                .all(|surface| surface.kind == LinearSurfaceKind::Esurface)
        );
        assert!(
            expression
                .orientations
                .iter()
                .flat_map(|orientation| &orientation.variants)
                .any(|variant| !variant.numerator_surfaces.is_empty())
        );
    }

    #[test]
    fn cff_generation_repeated_high_power_channel_can_use_uniform_scale() {
        let parsed = parsed_fixture("box_pow3.dot");
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Cff,
                energy_degree_bounds: vec![(0, 1), (1, 1), (3, 4)],
                numerator_sampling_scale: NumeratorSamplingScaleMode::All,
                preserve_internal_edges_as_four_d_denominators: Vec::new(),
            },
        )
        .unwrap();

        assert!(
            expression
                .orientations
                .iter()
                .flat_map(|orientation| &orientation.edge_energy_map)
                .any(LinearEnergyExpr::uses_uniform_scale)
        );
        assert!(
            expression
                .orientations
                .iter()
                .flat_map(|orientation| &orientation.variants)
                .any(|variant| variant.uniform_scale_power > 0)
        );
        assert!(
            expression
                .surfaces
                .linear_surface_cache
                .iter()
                .all(|surface| surface.kind == LinearSurfaceKind::Esurface)
        );
    }

    #[test]
    fn cff_generation_builds_one_loop_high_power_with_repeated_spectators() {
        let parsed = parsed_fixture("box_pow3.dot");
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Cff,
                energy_degree_bounds: vec![(0, 3), (1, 1), (2, 1), (3, 1)],
                ..Default::default()
            },
        )
        .unwrap();

        assert!(!expression.orientations.is_empty());
        assert!(
            expression
                .surfaces
                .linear_surface_cache
                .iter()
                .all(|surface| surface.kind == LinearSurfaceKind::Esurface)
        );
        assert!(
            expression
                .orientations
                .iter()
                .flat_map(|orientation| &orientation.variants)
                .any(|variant| {
                    variant.origin.as_deref() == Some("bounded_degree_known_factor_cff")
                })
        );
    }

    #[test]
    fn all_sampling_scale_mode_activates_for_linear_energy_degree() {
        assert!(!NumeratorSamplingScaleMode::None.is_active_for_degree(1));
        assert!(!NumeratorSamplingScaleMode::BeyondQuadratic.is_active_for_degree(2));
        assert!(NumeratorSamplingScaleMode::BeyondQuadratic.is_active_for_degree(3));
        assert!(NumeratorSamplingScaleMode::All.is_active_for_degree(1));
    }

    #[test]
    fn all_sampling_scale_mode_affects_quadratic_reconstruction() {
        let parsed = parsed_fixture("box_pow3.dot");
        let beyond_quadratic = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Cff,
                energy_degree_bounds: vec![(3, 2)],
                numerator_sampling_scale: NumeratorSamplingScaleMode::BeyondQuadratic,
                preserve_internal_edges_as_four_d_denominators: Vec::new(),
            },
        )
        .unwrap();
        let all = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Cff,
                energy_degree_bounds: vec![(3, 2)],
                numerator_sampling_scale: NumeratorSamplingScaleMode::All,
                preserve_internal_edges_as_four_d_denominators: Vec::new(),
            },
        )
        .unwrap();

        assert!(
            !beyond_quadratic
                .orientations
                .iter()
                .flat_map(|orientation| &orientation.edge_energy_map)
                .any(LinearEnergyExpr::uses_uniform_scale),
            "beyond-quadratic mode must not use M as a quadratic sampling node"
        );
        assert!(
            all.orientations
                .iter()
                .flat_map(|orientation| &orientation.edge_energy_map)
                .any(LinearEnergyExpr::uses_uniform_scale),
            "all mode must use M as a quadratic sampling node"
        );
    }
}
