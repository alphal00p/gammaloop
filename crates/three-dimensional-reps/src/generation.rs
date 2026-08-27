use std::{
    collections::{BTreeMap, BTreeSet, HashMap},
    ops::{Add, Neg, Sub},
};

use bincode_trait_derive::{Decode, Encode};
use itertools::Itertools;
use linnet::half_edge::involution::{EdgeIndex, EdgeVec, Orientation};
use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::{
    OrientationData, OrientationExpression, OrientationID, ParsedGraph, ThreeDExpression,
    cff_recursion::enumerate_cff_surface_chains,
    cut_structure::{ContourClosure, energy_residues},
    energy_bounds::normalize_energy_degree_bounds,
    expression::{ResidualDenominator, assign_numerator_map_labels},
    graph_io::{
        ParsedGraphExternalEdge, ParsedGraphInitialStateCutEdge, ParsedGraphInternalEdge,
        ThreeDGraphSource, repeated_groups,
    },
    graph_signatures::MomentumSignature,
    surface::{
        HybridSurfaceID, LinearEnergyExpr, LinearSurface, LinearSurfaceID, LinearSurfaceKind,
        RationalAtomExt, SurfaceOrigin, rational_coeff_atom, rational_coeff_new,
        rational_coeff_one,
    },
    tree::{NodeId, Tree},
    utils::{
        Rational, RationalExt, binomial, determinant_i32_is_nonzero, rank_i64, rank_rational,
        rational_pow_i64, solve_rational_system,
    },
};

#[derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize, Default,
)]
#[serde(rename_all = "snake_case")]
pub enum RepresentationMode {
    #[default]
    Cff,
    Ltd,
}

#[derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize, Default,
)]
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

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct Generate3DExpressionOptions {
    #[serde(default)]
    pub representation: RepresentationMode,
    /// `None` keeps the legacy numerator class, which is affine in every EMR
    /// edge energy. `Some` selects an explicit bounded class; omitted edges
    /// then have degree zero, including when the vector itself is empty.
    #[serde(default)]
    pub energy_degree_bounds: Option<Vec<(usize, usize)>>,
    pub numerator_sampling_scale: NumeratorSamplingScaleMode,
    /// Internal edge IDs whose unintegrated four-dimensional denominators stay
    /// outside the causal CFF graph. Their affine energy maps and typed
    /// residual factors remain in the generated result.
    #[serde(default)]
    pub preserve_internal_edges_as_four_d_denominators: Vec<usize>,
}

impl Default for Generate3DExpressionOptions {
    fn default() -> Self {
        Self {
            representation: RepresentationMode::Cff,
            energy_degree_bounds: None,
            numerator_sampling_scale: NumeratorSamplingScaleMode::None,
            preserve_internal_edges_as_four_d_denominators: Vec::new(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Encode, Decode)]
pub enum CffEnergyFactorOwnership {
    GlobalSourceProduct,
    VariantLocal,
}

/// Ownership of the on-shell energy factors for one denominator-connected
/// component. Edge IDs use the internal-energy namespace of the generated
/// expression. Component boundaries are retained when disconnected products
/// are formed.
#[derive(Debug, Clone, PartialEq, Eq, Encode, Decode)]
pub struct CffEnergyFactorComponent {
    pub internal_edge_ids: Vec<usize>,
    pub ownership: CffEnergyFactorOwnership,
}

impl CffEnergyFactorComponent {
    fn remap_internal_edges(&self, edge_map: &BTreeMap<usize, usize>) -> Self {
        let mut internal_edge_ids = self
            .internal_edge_ids
            .iter()
            .map(|edge_id| edge_map.get(edge_id).copied().unwrap_or(*edge_id))
            .collect::<Vec<_>>();
        internal_edge_ids.sort_unstable();
        internal_edge_ids.dedup();
        Self {
            internal_edge_ids,
            ownership: self.ownership,
        }
    }
}

/// The uniform prefactor sign inserted by CFF generation.
///
/// For each connected source this contains the shared core's
/// `(-1)^(L-1)` contour convention. Pure CFF also includes its uniform
/// duplicate-denominator sign. Duplicate signs introduced only inside
/// individual generalized variants are deliberately not represented here.
/// The metadata accompanies the generated expression so consumers can bridge
/// their own prefactor convention without reconstructing provenance from its
/// final algebra.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Encode, Decode)]
pub struct CffGlobalPrefactorSign {
    odd: bool,
}

impl CffGlobalPrefactorSign {
    fn from_exponent(exponent: usize) -> Self {
        Self {
            odd: !exponent.is_multiple_of(2),
        }
    }

    fn product(self, rhs: Self) -> Self {
        Self {
            odd: self.odd != rhs.odd,
        }
    }

    pub const fn factor(self) -> i64 {
        if self.odd { -1 } else { 1 }
    }
}

#[derive(Debug, Clone)]
pub struct GeneratedThreeDExpression<E = (), H = ()> {
    pub expression: ThreeDExpression<OrientationID, E, H>,
    pub energy_factor_ownership: CffEnergyFactorOwnership,
    /// Transient generation metadata. Persisted expressions retain the
    /// aggregate ownership above; exact-source consumers use these component
    /// boundaries before the converted expression is stored.
    pub energy_factor_components: Vec<CffEnergyFactorComponent>,
    /// Transient source-edge bounds supplied to this generation call. These always
    /// remain in the generator input's namespace; higher-level consumers that also
    /// need physical-parent bounds must retain them as separate metadata.
    pub source_energy_degree_bounds: Vec<(usize, usize)>,
    pub core_global_prefactor_sign: CffGlobalPrefactorSign,
}

// Keep the established persisted layout stable: component ownership and source
// bounds are generation-time metadata and are deliberately reconstructed only
// while a source is generated, rather than invalidating stored GammaLoop processes.
impl<E, H> bincode::Encode for GeneratedThreeDExpression<E, H>
where
    ThreeDExpression<OrientationID, E, H>: bincode::Encode,
{
    fn encode<Encoder: bincode::enc::Encoder>(
        &self,
        encoder: &mut Encoder,
    ) -> std::result::Result<(), bincode::error::EncodeError> {
        bincode::Encode::encode(&self.expression, encoder)?;
        bincode::Encode::encode(&self.energy_factor_ownership, encoder)?;
        bincode::Encode::encode(&self.core_global_prefactor_sign, encoder)
    }
}

impl<E, H, Context> bincode::Decode<Context> for GeneratedThreeDExpression<E, H>
where
    Context: symbolica::state::HasStateMap,
    ThreeDExpression<OrientationID, E, H>: bincode::Decode<Context>,
    CffEnergyFactorOwnership: bincode::Decode<Context>,
    CffGlobalPrefactorSign: bincode::Decode<Context>,
{
    fn decode<Decoder: bincode::de::Decoder<Context = Context>>(
        decoder: &mut Decoder,
    ) -> std::result::Result<Self, bincode::error::DecodeError> {
        Ok(Self {
            expression: bincode::Decode::decode(decoder)?,
            energy_factor_ownership: bincode::Decode::decode(decoder)?,
            energy_factor_components: Vec::new(),
            source_energy_degree_bounds: Vec::new(),
            core_global_prefactor_sign: bincode::Decode::decode(decoder)?,
        })
    }
}

#[cfg(test)]
mod generated_expression_persistence_tests {
    use symbolica::state::{HasStateMap, StateMap};

    use super::*;

    struct UnusedStateMap;

    impl HasStateMap for UnusedStateMap {
        fn get_state_map(&self) -> &StateMap {
            panic!("an empty expression must not request Symbolica state during decoding")
        }
    }

    #[derive(Encode)]
    struct LegacyGeneratedThreeDExpression<'a> {
        expression: &'a ThreeDExpression<OrientationID>,
        energy_factor_ownership: CffEnergyFactorOwnership,
        core_global_prefactor_sign: CffGlobalPrefactorSign,
    }

    #[test]
    fn transient_energy_factor_components_preserve_the_legacy_wire_layout() {
        let generated = GeneratedThreeDExpression {
            expression: ThreeDExpression::new_empty(),
            energy_factor_ownership: CffEnergyFactorOwnership::VariantLocal,
            energy_factor_components: vec![CffEnergyFactorComponent {
                internal_edge_ids: vec![3, 4],
                ownership: CffEnergyFactorOwnership::VariantLocal,
            }],
            source_energy_degree_bounds: vec![(3, 2)],
            core_global_prefactor_sign: CffGlobalPrefactorSign::from_exponent(0),
        };

        let legacy = LegacyGeneratedThreeDExpression {
            expression: &generated.expression,
            energy_factor_ownership: generated.energy_factor_ownership,
            core_global_prefactor_sign: generated.core_global_prefactor_sign,
        };
        let legacy_bytes = bincode::encode_to_vec(legacy, bincode::config::standard()).unwrap();
        assert_eq!(
            bincode::encode_to_vec(&generated, bincode::config::standard()).unwrap(),
            legacy_bytes,
        );

        let (decoded, bytes_read): (GeneratedThreeDExpression, _) =
            bincode::decode_from_slice_with_context(
                &legacy_bytes,
                bincode::config::standard(),
                UnusedStateMap,
            )
            .unwrap();
        assert_eq!(bytes_read, legacy_bytes.len());
        assert_eq!(
            decoded.energy_factor_ownership,
            generated.energy_factor_ownership
        );
        assert_eq!(
            decoded.core_global_prefactor_sign,
            generated.core_global_prefactor_sign
        );
        assert!(decoded.energy_factor_components.is_empty());
        assert!(decoded.source_energy_degree_bounds.is_empty());
        assert_eq!(
            bincode::encode_to_vec(decoded, bincode::config::standard()).unwrap(),
            legacy_bytes,
        );
    }
}

#[derive(Debug, Error)]
pub enum GenerationError {
    #[error(
        "three-dimensional representation mode {mode:?} is not implemented; only CFF is currently supported"
    )]
    NotImplemented { mode: RepresentationMode },
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
    #[error("{0}")]
    EnergyBounds(#[from] crate::energy_bounds::EnergyBoundsError),
    #[error("{0}")]
    GraphIo(#[from] crate::graph_io::GraphIoError),
    #[error(
        "energy-degree bound was requested for edge {edge_id}, but this is not an internal edge of the selected graph"
    )]
    UnknownEnergyDegreeBoundEdge { edge_id: usize },
    #[error(
        "CFF preservation was requested for edge {edge_id}, but this is not an internal edge of the selected graph"
    )]
    UnknownCffPreservedInternalEdge { edge_id: usize },
    #[error(
        "initial-state cut edge {edge_id} is an external alias and cannot be preserved as a residual four-dimensional denominator"
    )]
    CffPreservedInitialStateCutEdge { edge_id: usize },
    #[error(
        "preserved four-dimensional denominator edge {edge_id} retains loop-energy coefficients {loop_signature:?}; a CFF-external tree edge must have purely external momentum"
    )]
    CffPreservedEdgeCarriesLoopEnergy {
        edge_id: usize,
        loop_signature: Vec<i32>,
    },
    #[error(
        "contracting the last denominator carrying loop-energy variable {variable} leaves an unlocalized quotient of power {power}"
    )]
    UnlocalizedEnergyContact { variable: usize, power: usize },
    #[error(
        "disconnected denominator components {left:?} and {right:?} share loop-energy variables {variables:?}; component products require disjoint loop subspaces"
    )]
    DisconnectedComponentsShareLoopVariables {
        left: Vec<usize>,
        right: Vec<usize>,
        variables: Vec<usize>,
    },
}

pub type Result<T> = std::result::Result<T, GenerationError>;

pub fn generate_3d_expression<G: ThreeDGraphSource + ?Sized>(
    graph: &G,
    options: &Generate3DExpressionOptions,
) -> Result<GeneratedThreeDExpression> {
    if options.representation != RepresentationMode::Cff {
        return Err(GenerationError::NotImplemented {
            mode: options.representation,
        });
    }
    let parsed = graph.to_three_d_parsed_graph()?;
    let Some(edge_map) = graph.energy_edge_index_map(&parsed) else {
        let mut generated = generate_3d_expression_from_parsed_generated(&parsed, options)?;
        generated.source_energy_degree_bounds =
            options.energy_degree_bounds.clone().unwrap_or_default();
        return Ok(generated);
    };
    let mut local_options = options.clone();
    local_options.energy_degree_bounds = options
        .energy_degree_bounds
        .as_deref()
        .map(|bounds| -> Result<Vec<(usize, usize)>> {
            edge_map
                .remap_bounds_to_local(bounds)
                .map_err(|edge_id| GenerationError::UnknownEnergyDegreeBoundEdge { edge_id })
        })
        .transpose()?;
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
    let mut generated = generate_3d_expression_from_parsed_generated(&parsed, &local_options)?;
    generated.expression = generated
        .expression
        .remap_energy_edge_indices(&edge_map)
        .fuse_compatible_variants();
    generated.energy_factor_components = generated
        .energy_factor_components
        .iter()
        .map(|component| component.remap_internal_edges(&edge_map.internal))
        .collect();
    generated.source_energy_degree_bounds =
        options.energy_degree_bounds.clone().unwrap_or_default();
    assign_numerator_map_labels(&mut generated.expression.orientations);
    Ok(generated)
}

#[cfg(test)]
pub(crate) fn generate_3d_expression_from_parsed(
    parsed: &ParsedGraph,
    options: &Generate3DExpressionOptions,
) -> Result<ThreeDExpression<OrientationID>> {
    Ok(generate_3d_expression_from_parsed_generated(parsed, options)?.expression)
}

fn generate_3d_expression_from_parsed_generated(
    parsed: &ParsedGraph,
    options: &Generate3DExpressionOptions,
) -> Result<GeneratedThreeDExpression> {
    if options.representation != RepresentationMode::Cff {
        return Err(GenerationError::NotImplemented {
            mode: options.representation,
        });
    }
    if !options
        .preserve_internal_edges_as_four_d_denominators
        .is_empty()
    {
        return build_expression_preserving_internal_edges(parsed, options);
    }
    if let Some(mut generated) = generate_disconnected_component_product(parsed, options)? {
        generated.expression = generated.expression.fuse_compatible_variants();
        assign_numerator_map_labels(&mut generated.expression.orientations);
        return Ok(generated);
    }

    let (expression, energy_factor_ownership, core_global_prefactor_sign) =
        if cff_bounds_need_generalized_expression_from_options(parsed, options)? {
            (
                BoundedCffBuilder::new(parsed, options)?.build()?,
                CffEnergyFactorOwnership::VariantLocal,
                CffGlobalPrefactorSign::from_exponent(parsed.loop_names.len().saturating_sub(1)),
            )
        } else {
            let duplicate_excess = cff_duplicate_signature_excess(parsed);
            (
                generate_pure_cff_expression_from_parsed(parsed)?,
                CffEnergyFactorOwnership::GlobalSourceProduct,
                CffGlobalPrefactorSign::from_exponent(
                    parsed.loop_names.len().saturating_sub(1) + duplicate_excess,
                ),
            )
        };
    let mut expression = expression.fuse_compatible_variants();
    assign_numerator_map_labels(&mut expression.orientations);
    let internal_edge_ids = parsed.denominator_internal_edge_ids();
    let energy_factor_components = (!internal_edge_ids.is_empty())
        .then_some(CffEnergyFactorComponent {
            internal_edge_ids,
            ownership: energy_factor_ownership,
        })
        .into_iter()
        .collect();
    Ok(GeneratedThreeDExpression {
        expression,
        energy_factor_ownership,
        energy_factor_components,
        source_energy_degree_bounds: Vec::new(),
        core_global_prefactor_sign,
    })
}

fn generate_pure_cff_expression_from_parsed(
    parsed: &ParsedGraph,
) -> Result<ThreeDExpression<OrientationID>> {
    let signatures = parsed
        .internal_edges
        .iter()
        .map(|edge| edge.signature.clone())
        .collect::<Vec<_>>();
    let mut factorized = LowerSectorCffBuilder::new(parsed);
    if parsed.initial_state_cut_edges.is_empty()
        && factorized.vector_matroid_components(&signatures).len() > 1
    {
        // A graph may be vertex-connected while its rational denominator
        // factorizes into independent loop-energy components (for example a
        // tadpole attached at one vertex). Construct its causal product in
        // those independent coordinates instead of inventing surfaces which
        // mix components merely because they share an incidence vertex.
        factorized.force_component_factorization = true;
        return factorized.build();
    }
    generate_pure_cff_expression_from_parsed_with_duplicate_sign(parsed, true)
}

fn build_expression_preserving_internal_edges(
    parsed: &ParsedGraph,
    options: &Generate3DExpressionOptions,
) -> Result<GeneratedThreeDExpression> {
    let preserved = options
        .preserve_internal_edges_as_four_d_denominators
        .iter()
        .copied()
        .collect::<BTreeSet<_>>();
    for edge_id in &preserved {
        let Some(edge) = parsed.internal_edges.get(*edge_id) else {
            return Err(GenerationError::UnknownCffPreservedInternalEdge { edge_id: *edge_id });
        };
        if parsed.is_initial_state_cut_edge(*edge_id) {
            return Err(GenerationError::CffPreservedInitialStateCutEdge { edge_id: *edge_id });
        }
        if edge
            .signature
            .loop_signature
            .iter()
            .any(|coefficient| *coefficient != 0)
        {
            return Err(GenerationError::CffPreservedEdgeCarriesLoopEnergy {
                edge_id: *edge_id,
                loop_signature: edge.signature.loop_signature.clone(),
            });
        }
    }

    let (active_parsed, active_to_orig) = contract_preserved_parsed_edges(parsed, &preserved);
    if active_parsed.denominator_internal_edge_ids().is_empty() {
        return expression_with_only_preserved_edges(parsed, &preserved);
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
        .as_ref()
        .map(|bounds| {
            bounds
                .iter()
                .filter_map(|(edge_id, degree)| {
                    if *edge_id >= parsed.internal_edges.len() {
                        return Some(Err(GenerationError::UnknownEnergyDegreeBoundEdge {
                            edge_id: *edge_id,
                        }));
                    }
                    if preserved.contains(edge_id) {
                        return None;
                    }
                    Some(
                        orig_to_active
                            .get(edge_id)
                            .copied()
                            .map(|active_id| (active_id, *degree))
                            .ok_or(GenerationError::UnknownEnergyDegreeBoundEdge {
                                edge_id: *edge_id,
                            }),
                    )
                })
                .collect::<Result<Vec<_>>>()
        })
        .transpose()?;

    let generated = generate_3d_expression_from_parsed_generated(&active_parsed, &active_options)?;
    let active_edge_map = active_to_orig
        .iter()
        .enumerate()
        .map(|(active_id, orig_id)| (active_id, *orig_id))
        .collect::<BTreeMap<_, _>>();
    let energy_factor_components = generated
        .energy_factor_components
        .iter()
        .map(|component| component.remap_internal_edges(&active_edge_map))
        .collect();
    Ok(GeneratedThreeDExpression {
        expression: lift_expression_to_preserved_graph(
            parsed,
            &generated.expression,
            &active_to_orig,
            &preserved,
        )?,
        energy_factor_ownership: generated.energy_factor_ownership,
        energy_factor_components,
        source_energy_degree_bounds: Vec::new(),
        core_global_prefactor_sign: generated.core_global_prefactor_sign,
    })
}

fn expression_with_only_preserved_edges(
    parsed: &ParsedGraph,
    preserved: &BTreeSet<usize>,
) -> Result<GeneratedThreeDExpression> {
    if !parsed.loop_names.is_empty() {
        return Err(GenerationError::SingularBasis);
    }

    let signatures = parsed
        .internal_edges
        .iter()
        .map(|edge| edge.signature.clone())
        .collect::<Vec<_>>();
    let mut edge_energy_map = edge_q0_from_loop_exprs(&signatures, &[]);
    apply_initial_state_cut_edge_energy_exprs(parsed, &mut edge_energy_map);
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
            denominator_edges: Vec::new(),
            denominator_surface_signs: BTreeMap::new(),
            denominator_edge_support_signs: BTreeMap::new(),
            uniform_scale_power: 0,
            numerator_surfaces: Vec::new(),
            denominator: Tree::from_root(HybridSurfaceID::Unit),
        }],
    });
    expression.residual_denominators = preserved
        .iter()
        .map(|edge_id| {
            ResidualDenominator::new(
                EdgeIndex(*edge_id),
                Some(parsed.internal_edges[*edge_id].label.clone()),
            )
        })
        .collect();
    assign_numerator_map_labels(&mut expression.orientations);
    Ok(GeneratedThreeDExpression {
        expression,
        energy_factor_ownership: CffEnergyFactorOwnership::GlobalSourceProduct,
        energy_factor_components: Vec::new(),
        source_energy_degree_bounds: Vec::new(),
        core_global_prefactor_sign: CffGlobalPrefactorSign::default(),
    })
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

fn lift_expression_to_preserved_graph(
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
        apply_initial_state_cut_edge_energy_exprs(parsed, &mut edge_energy_map);
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
                denominator_edges: variant
                    .denominator_edges
                    .iter()
                    .map(|edge_id| EdgeIndex(active_edge_map[&edge_id.0]))
                    .collect(),
                denominator_surface_signs: variant
                    .denominator_surface_signs
                    .iter()
                    .map(|(surface_id, sign)| (map_surface_id(*surface_id, &surface_map), *sign))
                    .collect(),
                denominator_edge_support_signs: map_edge_support_signs(
                    &variant.denominator_edge_support_signs,
                    &active_edge_map,
                ),
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

#[derive(Debug, Clone)]
struct ParsedComponentEmbedding {
    local_to_orig_edge: Vec<usize>,
    local_to_orig_loop: Vec<usize>,
}

fn generate_disconnected_component_product(
    parsed: &ParsedGraph,
    options: &Generate3DExpressionOptions,
) -> Result<Option<GeneratedThreeDExpression>> {
    let components = denominator_connected_components(parsed);
    if components.len() <= 1 {
        return Ok(None);
    }

    let loop_supports = components
        .iter()
        .map(|component| {
            (0..parsed.loop_names.len())
                .filter(|loop_id| {
                    component.iter().any(|edge_id| {
                        parsed.internal_edges[*edge_id].signature.loop_signature[*loop_id] != 0
                    })
                })
                .collect::<BTreeSet<_>>()
        })
        .collect::<Vec<_>>();
    for (left, right) in (0..components.len()).tuple_combinations() {
        let variables = loop_supports[left]
            .intersection(&loop_supports[right])
            .copied()
            .collect::<Vec<_>>();
        if !variables.is_empty() {
            return Err(GenerationError::DisconnectedComponentsShareLoopVariables {
                left: components[left].clone(),
                right: components[right].clone(),
                variables,
            });
        }
    }

    let component_expressions = components
        .into_iter()
        .map(|component_edges| {
            let (component, embedding) = project_denominator_component(parsed, &component_edges)?;
            let component_options =
                project_component_options(options, parsed.internal_edges.len(), &embedding)?;
            let generated =
                generate_3d_expression_from_parsed_generated(&component, &component_options)?;
            Ok((generated, embedding))
        })
        .collect::<Result<Vec<_>>>()?;

    let energy_factor_ownership = if component_expressions.iter().any(|(generated, _)| {
        generated.energy_factor_ownership == CffEnergyFactorOwnership::VariantLocal
    }) {
        CffEnergyFactorOwnership::VariantLocal
    } else {
        CffEnergyFactorOwnership::GlobalSourceProduct
    };
    let expressions = component_expressions
        .iter()
        .map(|(generated, embedding)| (&generated.expression, embedding))
        .collect::<Vec<_>>();
    let core_global_prefactor_sign = component_expressions
        .iter()
        .fold(CffGlobalPrefactorSign::default(), |sign, (generated, _)| {
            sign.product(generated.core_global_prefactor_sign)
        });
    let mut energy_factor_components = Vec::new();
    for (generated, embedding) in &component_expressions {
        let edge_map = embedding
            .local_to_orig_edge
            .iter()
            .enumerate()
            .map(|(local_id, orig_id)| (local_id, *orig_id))
            .collect::<BTreeMap<_, _>>();
        energy_factor_components.extend(
            generated
                .energy_factor_components
                .iter()
                .map(|component| component.remap_internal_edges(&edge_map)),
        );
    }
    Ok(Some(GeneratedThreeDExpression {
        expression: lift_component_expression_product(parsed, &expressions)?,
        energy_factor_ownership,
        energy_factor_components,
        source_energy_degree_bounds: Vec::new(),
        core_global_prefactor_sign,
    }))
}

fn denominator_connected_components(parsed: &ParsedGraph) -> Vec<Vec<usize>> {
    let denominator_edges = parsed.denominator_internal_edge_ids();
    if denominator_edges.is_empty() {
        return Vec::new();
    }

    let mut parent = BTreeMap::<usize, usize>::new();
    for edge_id in &denominator_edges {
        let edge = &parsed.internal_edges[*edge_id];
        parent.entry(edge.tail).or_insert(edge.tail);
        parent.entry(edge.head).or_insert(edge.head);
        union_nodes(&mut parent, edge.tail, edge.head);
    }

    let mut grouped = BTreeMap::<usize, Vec<usize>>::new();
    for edge_id in denominator_edges {
        let edge = &parsed.internal_edges[edge_id];
        let root = find_node_root(&mut parent, edge.tail);
        grouped.entry(root).or_default().push(edge_id);
    }
    grouped
        .into_values()
        .map(|mut edges| {
            edges.sort_unstable();
            edges
        })
        .collect()
}

fn project_denominator_component(
    parsed: &ParsedGraph,
    component_edges: &[usize],
) -> Result<(ParsedGraph, ParsedComponentEmbedding)> {
    let component_edge_set = component_edges.iter().copied().collect::<BTreeSet<_>>();
    let active_loop_columns = parsed
        .loop_names
        .iter()
        .enumerate()
        .filter_map(|(loop_id, _)| {
            component_edges
                .iter()
                .any(|edge_id| {
                    parsed.internal_edges[*edge_id].signature.loop_signature[loop_id] != 0
                })
                .then_some(loop_id)
        })
        .collect::<Vec<_>>();
    let loop_names = active_loop_columns
        .iter()
        .map(|loop_id| parsed.loop_names[*loop_id].clone())
        .collect::<Vec<_>>();

    let component_nodes = component_edges
        .iter()
        .flat_map(|edge_id| {
            let edge = &parsed.internal_edges[*edge_id];
            [edge.tail, edge.head]
        })
        .collect::<BTreeSet<_>>();
    let old_to_new = component_nodes
        .iter()
        .copied()
        .enumerate()
        .map(|(new_id, old_id)| (old_id, new_id))
        .collect::<BTreeMap<_, _>>();

    let mut orig_to_local = BTreeMap::new();
    let mut local_to_orig_edge = Vec::new();
    let internal_edges = component_edges
        .iter()
        .copied()
        .enumerate()
        .map(|(local_id, orig_id)| {
            orig_to_local.insert(orig_id, local_id);
            local_to_orig_edge.push(orig_id);
            let edge = &parsed.internal_edges[orig_id];
            ParsedGraphInternalEdge {
                edge_id: local_id,
                tail: old_to_new[&edge.tail],
                head: old_to_new[&edge.head],
                label: edge.label.clone(),
                mass_key: edge.mass_key.clone(),
                signature: MomentumSignature {
                    loop_signature: active_loop_columns
                        .iter()
                        .map(|loop_id| edge.signature.loop_signature[*loop_id])
                        .collect(),
                    external_signature: edge.signature.external_signature.clone(),
                },
                had_pow: edge.had_pow,
            }
        })
        .collect::<Vec<_>>();

    let mut next_external_edge_id = 0usize;
    let mut external_edges = parsed
        .external_edges
        .iter()
        .filter_map(|edge| {
            let source = edge
                .source
                .and_then(|source| old_to_new.get(&source).copied());
            let destination = edge
                .destination
                .and_then(|destination| old_to_new.get(&destination).copied());
            (source.is_some() || destination.is_some()).then(|| {
                let edge_id = next_external_edge_id;
                next_external_edge_id += 1;
                ParsedGraphExternalEdge {
                    edge_id,
                    source,
                    destination,
                    label: edge.label.clone(),
                    external_coefficients: edge.external_coefficients.clone(),
                }
            })
        })
        .collect::<Vec<_>>();
    for edge in &parsed.internal_edges {
        if component_edge_set.contains(&edge.edge_id)
            || parsed.is_initial_state_cut_edge(edge.edge_id)
        {
            continue;
        }
        if let Some(source) = old_to_new.get(&edge.tail).copied() {
            external_edges.push(ParsedGraphExternalEdge {
                edge_id: next_external_edge_id,
                source: Some(source),
                destination: None,
                label: format!("{}_source_boundary", edge.label),
                external_coefficients: edge.signature.external_signature.clone(),
            });
            next_external_edge_id += 1;
        }
        if let Some(destination) = old_to_new.get(&edge.head).copied() {
            external_edges.push(ParsedGraphExternalEdge {
                edge_id: next_external_edge_id,
                source: None,
                destination: Some(destination),
                label: format!("{}_sink_boundary", edge.label),
                external_coefficients: edge.signature.external_signature.clone(),
            });
            next_external_edge_id += 1;
        }
    }

    let initial_state_cut_edges = parsed
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
        .collect::<Vec<_>>();

    Ok((
        ParsedGraph {
            internal_edges,
            external_edges,
            initial_state_cut_edges,
            loop_names,
            external_names: parsed.external_names.clone(),
            node_name_to_internal: old_to_new
                .into_iter()
                .map(|(old_id, new_id)| (format!("component_{old_id}"), new_id))
                .collect(),
        },
        ParsedComponentEmbedding {
            local_to_orig_edge,
            local_to_orig_loop: active_loop_columns,
        },
    ))
}

fn project_component_options(
    options: &Generate3DExpressionOptions,
    source_edge_count: usize,
    embedding: &ParsedComponentEmbedding,
) -> Result<Generate3DExpressionOptions> {
    let energy_degree_bounds = options
        .energy_degree_bounds
        .as_deref()
        .map(|bounds| -> Result<Vec<(usize, usize)>> {
            let bounds = normalize_energy_degree_bounds(bounds, source_edge_count)?;
            Ok(embedding
                .local_to_orig_edge
                .iter()
                .enumerate()
                .filter_map(|(local_edge, source_edge)| {
                    let degree = bounds[*source_edge];
                    (degree != 0).then_some((local_edge, degree))
                })
                .collect())
        })
        .transpose()?;
    Ok(Generate3DExpressionOptions {
        representation: options.representation,
        energy_degree_bounds,
        numerator_sampling_scale: options.numerator_sampling_scale,
        preserve_internal_edges_as_four_d_denominators: Vec::new(),
    })
}

fn lift_component_expression_product(
    parsed: &ParsedGraph,
    components: &[(&ThreeDExpression<OrientationID>, &ParsedComponentEmbedding)],
) -> Result<ThreeDExpression<OrientationID>> {
    let mut expression = ThreeDExpression::<OrientationID>::new_empty();
    let mut partials = vec![OrientationExpression {
        data: OrientationData {
            orientation: EdgeVec::from_iter(
                (0..parsed.internal_edges.len()).map(|_| Orientation::Undirected),
            ),
            label: Some("component_product_identity".to_string()),
            numerator_map_index: None,
        },
        loop_energy_map: vec![LinearEnergyExpr::zero(); parsed.loop_names.len()],
        edge_energy_map: vec![LinearEnergyExpr::zero(); parsed.internal_edges.len()],
        variants: vec![crate::expression::CFFVariant {
            origin: Some("component_product_identity".to_string()),
            prefactor: rational_coeff_one(),
            half_edges: Vec::new(),
            denominator_edges: Vec::new(),
            denominator_surface_signs: BTreeMap::new(),
            denominator_edge_support_signs: BTreeMap::new(),
            uniform_scale_power: 0,
            numerator_surfaces: Vec::new(),
            denominator: Tree::from_root(HybridSurfaceID::Unit),
        }],
    }];

    for (component_expression, embedding) in components {
        let edge_map = embedding
            .local_to_orig_edge
            .iter()
            .enumerate()
            .map(|(local_id, orig_id)| (local_id, *orig_id))
            .collect::<BTreeMap<_, _>>();
        let surface_map = {
            let mut interner = LinearSurfaceInterner::new(&mut expression);
            component_expression
                .surfaces
                .linear_surface_cache
                .iter_enumerated()
                .map(|(id, surface)| {
                    (
                        HybridSurfaceID::Linear(id),
                        interner.intern_with_origin(
                            surface.expression.clone().remap_internal_edges(&edge_map),
                            surface.origin,
                            surface.numerator_only,
                        ),
                    )
                })
                .collect::<HashMap<_, _>>()
        };

        let mut next_partials = Vec::new();
        for partial in &partials {
            for component_orientation in &component_expression.orientations {
                let mut merged = partial.clone();
                for (local_id, orig_id) in embedding.local_to_orig_edge.iter().enumerate() {
                    merged.data.orientation[EdgeIndex(*orig_id)] =
                        component_orientation.data.orientation[EdgeIndex(local_id)];
                    merged.edge_energy_map[*orig_id] = component_orientation.edge_energy_map
                        [local_id]
                        .clone()
                        .remap_internal_edges(&edge_map);
                }
                for (local_loop_id, orig_loop_id) in embedding.local_to_orig_loop.iter().enumerate()
                {
                    merged.loop_energy_map[*orig_loop_id] = component_orientation.loop_energy_map
                        [local_loop_id]
                        .clone()
                        .remap_internal_edges(&edge_map);
                }
                let (orientation_data, label) =
                    orientation_from_edge_exprs(&merged.edge_energy_map);
                merged.data.orientation = orientation_data;
                merged.data.label = Some(label);

                let mut variants = Vec::new();
                for partial_variant in &partial.variants {
                    for component_variant in &component_orientation.variants {
                        variants.push(product_variants(
                            partial_variant,
                            component_variant,
                            &edge_map,
                            &surface_map,
                        )?);
                    }
                }
                merged.variants = variants;
                next_partials.push(merged);
            }
        }
        partials = next_partials;
    }

    expression.orientations = partials.into_iter().collect();
    Ok(expression)
}

fn product_variants(
    lhs: &crate::expression::CFFVariant,
    rhs: &crate::expression::CFFVariant,
    rhs_edge_map: &BTreeMap<usize, usize>,
    rhs_surface_map: &HashMap<HybridSurfaceID, HybridSurfaceID>,
) -> Result<crate::expression::CFFVariant> {
    let mut half_edges = lhs.half_edges.clone();
    half_edges.extend(
        rhs.half_edges
            .iter()
            .map(|edge_id| EdgeIndex(rhs_edge_map[&edge_id.0])),
    );
    half_edges.sort_unstable();

    let mut denominator_edges = lhs.denominator_edges.clone();
    denominator_edges.extend(
        rhs.denominator_edges
            .iter()
            .map(|edge_id| EdgeIndex(rhs_edge_map[&edge_id.0])),
    );
    denominator_edges.sort_unstable();

    let mut numerator_surfaces = lhs.numerator_surfaces.clone();
    numerator_surfaces.extend(
        rhs.numerator_surfaces
            .iter()
            .map(|surface_id| map_surface_id(*surface_id, rhs_surface_map)),
    );

    let mut denominator_surface_signs = lhs.denominator_surface_signs.clone();
    for (surface_id, sign) in &rhs.denominator_surface_signs {
        let mapped = map_surface_id(*surface_id, rhs_surface_map);
        *denominator_surface_signs.entry(mapped).or_insert(1) *= *sign;
    }
    let mut denominator_edge_support_signs = lhs.denominator_edge_support_signs.clone();
    for (support_edges, sign) in
        map_edge_support_signs(&rhs.denominator_edge_support_signs, rhs_edge_map)
    {
        *denominator_edge_support_signs
            .entry(support_edges)
            .or_insert(1) *= sign;
    }

    Ok(crate::expression::CFFVariant {
        origin: Some(format!(
            "component_product:{}:{}",
            lhs.origin.as_deref().unwrap_or("lhs"),
            rhs.origin.as_deref().unwrap_or("rhs")
        )),
        prefactor: rational_to_coefficient(
            rational_from_coefficient(&lhs.prefactor) * rational_from_coefficient(&rhs.prefactor),
        )?,
        half_edges,
        denominator_edges,
        denominator_surface_signs,
        denominator_edge_support_signs,
        uniform_scale_power: lhs.uniform_scale_power + rhs.uniform_scale_power,
        numerator_surfaces,
        denominator: product_denominator_trees(&lhs.denominator, &rhs.denominator, rhs_surface_map),
    })
}

fn product_denominator_trees(
    lhs: &Tree<HybridSurfaceID>,
    rhs: &Tree<HybridSurfaceID>,
    rhs_surface_map: &HashMap<HybridSurfaceID, HybridSurfaceID>,
) -> Tree<HybridSurfaceID> {
    let lhs_chains = denominator_tree_chains(lhs);
    let rhs_chains = denominator_tree_chains(rhs);
    let chains = lhs_chains
        .iter()
        .cartesian_product(rhs_chains.iter())
        .map(|(lhs_chain, rhs_chain)| {
            lhs_chain
                .iter()
                .copied()
                .chain(
                    rhs_chain
                        .iter()
                        .map(|surface_id| map_surface_id(*surface_id, rhs_surface_map)),
                )
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    denominator_tree_from_chains(&chains)
}

fn generate_pure_cff_expression_from_parsed_with_duplicate_sign(
    parsed: &ParsedGraph,
    include_duplicate_excess_sign: bool,
) -> Result<ThreeDExpression<OrientationID>> {
    let duplicate_excess = if include_duplicate_excess_sign {
        cff_duplicate_signature_excess(parsed)
    } else {
        0
    };
    generate_pure_cff_expression_from_parsed_with_duplicate_excess(parsed, duplicate_excess)
}

fn generate_pure_cff_expression_from_parsed_with_duplicate_excess(
    parsed: &ParsedGraph,
    duplicate_excess: usize,
) -> Result<ThreeDExpression<OrientationID>> {
    let signatures = parsed
        .internal_edges
        .iter()
        .map(|edge| edge.signature.clone())
        .collect::<Vec<_>>();
    let n_internal = signatures.len();
    let denominator_edge_ids = parsed.denominator_internal_edge_ids();
    let basis = LowerSectorCffBuilder::component_basis_edges(&signatures, &denominator_edge_ids);
    let n_loops = parsed.loop_names.len();
    if basis.len() != n_loops {
        return Err(GenerationError::SingularBasis);
    }
    let overall_sign = if (n_loops.saturating_sub(1) + duplicate_excess).is_multiple_of(2) {
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

        let mut denominator_edge_energy_map = vec![LinearEnergyExpr::zero(); n_internal];
        for edge_index in &denominator_edge_ids {
            denominator_edge_energy_map[*edge_index] =
                LinearEnergyExpr::ose(EdgeIndex(*edge_index), i64::from(signs[*edge_index]));
        }
        apply_initial_state_cut_edge_energy_exprs(parsed, &mut denominator_edge_energy_map);
        let loop_energy_map = solve_loop_energy_from_target_edge_exprs(
            &signatures,
            &basis,
            &denominator_edge_energy_map,
        )?;
        let edge_energy_map = denominator_edge_energy_map;
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
            denominator_edges: denominator_edge_ids
                .iter()
                .copied()
                .map(EdgeIndex)
                .collect(),
            denominator_surface_signs: BTreeMap::new(),
            denominator_edge_support_signs: BTreeMap::new(),
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

fn generate_simple_residue_basis_expression_from_parsed(
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
    let residues = energy_residues(&loop_line_signatures, &vec![ContourClosure::Below; n_loops])?;

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
            origin: Some("residue_basis".to_string()),
            prefactor,
            half_edges: basis,
            denominator_edges: denominator_edge_ids
                .iter()
                .copied()
                .map(EdgeIndex)
                .collect(),
            denominator_surface_signs: BTreeMap::new(),
            denominator_edge_support_signs: BTreeMap::new(),
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
struct ContactComponent {
    sample: i32,
    prefactor: Rational,
    half_edges: Vec<usize>,
    numerator_surfaces: Vec<HybridSurfaceID>,
}

#[derive(Debug, Clone)]
struct LogicalChannel {
    rep_edge: usize,
    members: Vec<usize>,
    power: usize,
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum KnownPolynomialVar {
    UniformScale,
    External(usize),
    Ose(usize),
    Loop(usize),
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
struct KnownMonomial(Vec<(KnownPolynomialVar, usize)>);

impl KnownMonomial {
    fn one() -> Self {
        Self(Vec::new())
    }

    fn var(var: KnownPolynomialVar) -> Self {
        Self(vec![(var, 1)])
    }

    fn mul(&self, rhs: &Self) -> Self {
        let mut powers = BTreeMap::<KnownPolynomialVar, usize>::new();
        for (var, power) in self.0.iter().chain(&rhs.0) {
            *powers.entry(*var).or_default() += *power;
        }
        Self(
            powers
                .into_iter()
                .filter(|(_, power)| *power != 0)
                .collect(),
        )
    }

    fn divides(&self, rhs: &Self) -> bool {
        let rhs_powers = rhs.0.iter().copied().collect::<BTreeMap<_, _>>();
        self.0
            .iter()
            .all(|(var, power)| rhs_powers.get(var).copied().unwrap_or(0) >= *power)
    }

    fn quotient_after_dividing_by(&self, divisor: &Self) -> Option<Self> {
        if !divisor.divides(self) {
            return None;
        }
        let mut powers = self.0.iter().copied().collect::<BTreeMap<_, _>>();
        for (var, divisor_power) in &divisor.0 {
            let entry = powers.get_mut(var)?;
            *entry -= *divisor_power;
        }
        Some(Self(
            powers
                .into_iter()
                .filter(|(_, power)| *power != 0)
                .collect(),
        ))
    }

    fn loop_variables(&self) -> impl Iterator<Item = usize> + '_ {
        self.0.iter().filter_map(|(var, _)| match var {
            KnownPolynomialVar::Loop(loop_id) => Some(*loop_id),
            _ => None,
        })
    }

    fn numerator_surface_exprs(
        &self,
        loop_exprs: &[LinearEnergyExpr],
    ) -> Result<Option<Vec<LinearEnergyExpr>>> {
        let mut out = Vec::new();
        for (var, power) in &self.0 {
            let expr = match var {
                KnownPolynomialVar::UniformScale => LinearEnergyExpr::uniform_scale(1),
                KnownPolynomialVar::External(edge_id) => {
                    LinearEnergyExpr::external(EdgeIndex(*edge_id), 1)
                }
                KnownPolynomialVar::Ose(edge_id) => LinearEnergyExpr::ose(EdgeIndex(*edge_id), 1),
                KnownPolynomialVar::Loop(loop_id) => loop_exprs
                    .get(*loop_id)
                    .cloned()
                    .ok_or(GenerationError::CffHigherEnergyPowerNotImplemented)?,
            }
            .canonical();
            if expr.is_zero() {
                return Ok(None);
            }
            if expr.is_one() {
                continue;
            }
            out.extend(std::iter::repeat_n(expr, *power));
        }
        Ok(Some(out))
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
struct KnownPolynomial {
    terms: BTreeMap<KnownMonomial, Rational>,
}

impl KnownPolynomial {
    fn zero() -> Self {
        Self::default()
    }

    fn one() -> Self {
        Self::from_monomial(KnownMonomial::one(), Rational::one())
    }

    fn from_monomial(monomial: KnownMonomial, coeff: Rational) -> Self {
        let mut out = Self::zero();
        out.add_term(monomial, coeff);
        out
    }

    fn variable(var: KnownPolynomialVar) -> Self {
        Self::from_monomial(KnownMonomial::var(var), Rational::one())
    }

    fn constant(value: Rational) -> Self {
        Self::from_monomial(KnownMonomial::one(), value)
    }

    fn q0_from_signature(signature: &MomentumSignature) -> Self {
        let mut out = Self::zero();
        for (loop_id, coeff) in signature.loop_signature.iter().enumerate() {
            if *coeff != 0 {
                out = out
                    + Self::variable(KnownPolynomialVar::Loop(loop_id))
                        .scale(Rational::from(i64::from(*coeff)));
            }
        }
        for (external_id, coeff) in signature.external_signature.iter().enumerate() {
            if *coeff != 0 {
                out = out
                    + Self::variable(KnownPolynomialVar::External(external_id))
                        .scale(Rational::from(i64::from(*coeff)));
            }
        }
        out
    }

    fn from_known_factor(parsed: &ParsedGraph, factor: &KnownLinearExpr) -> Self {
        let mut out = Self::constant(factor.constant.clone())
            + Self::variable(KnownPolynomialVar::UniformScale)
                .scale(factor.uniform_scale_coeff.clone());
        for (edge_id, coeff) in &factor.external_terms {
            out = out + Self::variable(KnownPolynomialVar::External(*edge_id)).scale(coeff.clone());
        }
        for (edge_id, coeff) in &factor.ose_terms {
            out = out + Self::variable(KnownPolynomialVar::Ose(*edge_id)).scale(coeff.clone());
        }
        for (edge_id, coeff) in &factor.var_terms {
            out = out
                + Self::q0_from_signature(&parsed.internal_edges[*edge_id].signature)
                    .scale(coeff.clone());
        }
        out
    }

    fn product_from_known_factors(parsed: &ParsedGraph, factors: &[KnownLinearExpr]) -> Self {
        factors.iter().fold(Self::one(), |acc, factor| {
            acc * Self::from_known_factor(parsed, factor)
        })
    }

    fn denominator_for_edge(parsed: &ParsedGraph, local_to_orig: &[usize], edge_id: usize) -> Self {
        let q0 = Self::q0_from_signature(&parsed.internal_edges[edge_id].signature);
        let ose = Self::variable(KnownPolynomialVar::Ose(local_to_orig[edge_id]));
        q0.clone() * q0 - ose.clone() * ose
    }

    fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    fn add_term(&mut self, monomial: KnownMonomial, coeff: Rational) {
        if coeff.is_zero() {
            return;
        }
        let entry = self.terms.entry(monomial).or_insert_with(Rational::zero);
        *entry = entry.clone() + coeff;
        if entry.is_zero() {
            self.terms.retain(|_, value| !value.is_zero());
        }
    }

    fn scale(mut self, scale: Rational) -> Self {
        if scale.is_zero() {
            return Self::zero();
        }
        for coeff in self.terms.values_mut() {
            *coeff = coeff.clone() * scale.clone();
        }
        self
    }

    fn leading_term(&self) -> Option<(KnownMonomial, Rational)> {
        self.terms
            .iter()
            .next_back()
            .map(|(monomial, coeff)| (monomial.clone(), coeff.clone()))
    }

    fn div_rem(&self, divisor: &Self) -> Result<(Self, Self)> {
        let Some((divisor_monomial, divisor_coeff)) = divisor.leading_term() else {
            return Err(GenerationError::CffHigherEnergyPowerNotImplemented);
        };
        let mut quotient = Self::zero();
        let mut remainder = self.clone();
        while let Some((lead_monomial, lead_coeff)) = remainder.leading_term() {
            let Some(quotient_monomial) =
                lead_monomial.quotient_after_dividing_by(&divisor_monomial)
            else {
                break;
            };
            let quotient_coeff = lead_coeff / divisor_coeff.clone();
            let quotient_term =
                Self::from_monomial(quotient_monomial.clone(), quotient_coeff.clone());
            quotient.add_term(quotient_monomial, quotient_coeff);
            remainder = remainder - quotient_term * divisor.clone();
        }
        Ok((quotient, remainder))
    }

    fn loop_variables(&self) -> BTreeSet<usize> {
        self.terms
            .keys()
            .flat_map(KnownMonomial::loop_variables)
            .collect()
    }
}

impl Add for KnownPolynomial {
    type Output = Self;

    fn add(mut self, rhs: Self) -> Self::Output {
        for (monomial, coeff) in rhs.terms {
            self.add_term(monomial, coeff);
        }
        self
    }
}

impl Neg for KnownPolynomial {
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        for coeff in self.terms.values_mut() {
            *coeff = -coeff.clone();
        }
        self
    }
}

impl Sub for KnownPolynomial {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl std::ops::Mul for KnownPolynomial {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut out = Self::zero();
        for (lhs_monomial, lhs_coeff) in self.terms {
            for (rhs_monomial, rhs_coeff) in &rhs.terms {
                out.add_term(
                    lhs_monomial.mul(rhs_monomial),
                    lhs_coeff.clone() * rhs_coeff.clone(),
                );
            }
        }
        out
    }
}

#[derive(Debug, Clone)]
struct KnownPolynomialBranch {
    denominator_edges: Vec<usize>,
    numerator: KnownPolynomial,
}

struct KnownPolynomialNormalForm<'a> {
    parsed: &'a ParsedGraph,
    local_to_orig: &'a [usize],
}

impl<'a> KnownPolynomialNormalForm<'a> {
    fn new(parsed: &'a ParsedGraph, local_to_orig: &'a [usize]) -> Self {
        Self {
            parsed,
            local_to_orig,
        }
    }

    fn branches(&self, numerator: KnownPolynomial) -> Result<Vec<KnownPolynomialBranch>> {
        let edges = (0..self.parsed.internal_edges.len()).collect::<Vec<_>>();
        let mut branches = Vec::new();
        self.reduce(numerator, edges.clone(), edges, &mut branches)?;
        Ok(branches)
    }

    fn reduce(
        &self,
        numerator: KnownPolynomial,
        denominator_edges: Vec<usize>,
        candidates: Vec<usize>,
        branches: &mut Vec<KnownPolynomialBranch>,
    ) -> Result<()> {
        if numerator.is_zero() {
            return Ok(());
        }
        let Some((active, rest)) = candidates.split_first() else {
            if !denominator_edges.is_empty() {
                branches.push(KnownPolynomialBranch {
                    denominator_edges,
                    numerator,
                });
            }
            return Ok(());
        };

        let denominator =
            KnownPolynomial::denominator_for_edge(self.parsed, self.local_to_orig, *active);
        let (quotient, remainder) = numerator.div_rem(&denominator)?;
        if !quotient.is_zero() {
            let quotient_denominators = denominator_edges
                .iter()
                .copied()
                .filter(|edge_id| edge_id != active)
                .collect::<Vec<_>>();
            let quotient_candidates = rest
                .iter()
                .copied()
                .filter(|edge_id| edge_id != active)
                .collect::<Vec<_>>();
            self.reduce(
                quotient,
                quotient_denominators,
                quotient_candidates,
                branches,
            )?;
        }
        self.reduce(remainder, denominator_edges, rest.to_vec(), branches)
    }
}

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

struct BoundedCffBuilder<'a> {
    parsed: &'a ParsedGraph,
    bounds: Vec<usize>,
    sampling_scale_mode: NumeratorSamplingScaleMode,
    expression: ThreeDExpression<OrientationID>,
    surface_index: HashMap<(LinearSurfaceKind, LinearEnergyExpr), HybridSurfaceID>,
}

impl<'a> BoundedCffBuilder<'a> {
    fn new(parsed: &'a ParsedGraph, options: &Generate3DExpressionOptions) -> Result<Self> {
        let bounds = normalize_energy_degree_bounds(
            options.energy_degree_bounds.as_deref().unwrap_or(&[]),
            parsed.internal_edges.len(),
        )?;
        Ok(Self::for_bounds(parsed, bounds).with_options(options))
    }

    fn for_bounds(parsed: &'a ParsedGraph, mut bounds: Vec<usize>) -> Self {
        // A repeated denominator group is one algebraic on-shell-energy
        // channel. Canonicalize only this private structural bound; the
        // caller's occurrence caps and factor-local numerator assignment stay
        // minimax-distributed across equivalent occurrences.
        for channel in KnownFactorCffBuilder::logical_channels(parsed) {
            let degree = channel.members.iter().map(|edge_id| bounds[*edge_id]).sum();
            for edge_id in channel.members {
                bounds[edge_id] = 0;
            }
            bounds[channel.rep_edge] = degree;
        }
        Self {
            parsed,
            bounds,
            sampling_scale_mode: NumeratorSamplingScaleMode::None,
            expression: ThreeDExpression::new_empty(),
            surface_index: HashMap::new(),
        }
    }

    fn with_options(mut self, options: &Generate3DExpressionOptions) -> Self {
        self.sampling_scale_mode = options.numerator_sampling_scale;
        self
    }

    fn build(mut self) -> Result<ThreeDExpression<OrientationID>> {
        let needs_generalized_expression =
            cff_bounds_need_generalized_expression(self.parsed, &self.bounds);
        if !needs_generalized_expression {
            return generate_pure_cff_expression_from_parsed(self.parsed);
        }
        let uniform_sampling_for_nonlinear_degree = self
            .bounds
            .iter()
            .any(|degree| *degree > 1 && self.sampling_scale_mode.is_active_for_degree(*degree));
        // A repeated denominator channel shares one physical EMR energy. Reconstruct its
        // aggregate numerator degree before taking a per-edge quadratic shortcut.
        let repeated_channel_needs_known_factor =
            KnownFactorCffBuilder::logical_channels(self.parsed)
                .iter()
                .any(|channel| {
                    channel
                        .members
                        .iter()
                        .map(|edge_id| self.bounds[*edge_id])
                        .sum::<usize>()
                        > 2
                });
        if repeated_channel_needs_known_factor {
            return KnownFactorCffBuilder::new(self.parsed, self.bounds, self.sampling_scale_mode)
                .build();
        }
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
        self.bounds.iter().any(|degree| *degree > 1)
            && self.parsed.loop_names.len() == 1
            && !self.has_duplicate_signature_ignoring_mass()
            && self.bounds.iter().all(|degree| *degree <= 2)
    }

    fn supports_quadratic_recursive(&self) -> bool {
        self.bounds.iter().any(|degree| *degree > 1)
            && self.bounds.iter().all(|degree| *degree <= 2)
    }

    fn supports_known_factor_recursive(&self) -> bool {
        self.bounds.iter().any(|degree| *degree > 1)
            || KnownFactorCffBuilder::logical_channels(self.parsed)
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
        has_duplicate_signature_ignoring_mass(self.parsed)
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

        let (subparsed, sub_to_orig) = self.project_parsed_edges(&[active_edge]);
        let powered_edges = KnownFactorCffBuilder::logical_channels(self.parsed)
            .into_iter()
            .flat_map(|channel| channel.members)
            .collect::<BTreeSet<_>>();
        // Once this edge supplies the lower-sector contact, other ordinary
        // simple poles need only sample the remaining numerator. Carry a
        // higher bound onward solely on powered channels, where the raised
        // residue genuinely differentiates that numerator. Recursing another
        // simple-pole contact here would count the same pinch twice.
        let sub_bounds = sub_to_orig
            .iter()
            .map(|orig_id| {
                if powered_edges.contains(orig_id) {
                    self.bounds[*orig_id]
                } else {
                    self.bounds[*orig_id].min(1)
                }
            })
            .collect::<Vec<_>>();
        let contact_expression = BoundedCffBuilder::for_bounds(&subparsed, sub_bounds)
            .build_quadratic_recursive(true)?;
        self.append_recursive_contact_terms(
            active_edge,
            &contact_expression,
            &sub_to_orig,
            &contact_weight_polys(self.bounds[active_edge]),
        )?;

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
                    self.assign_repeated_channel_sample(
                        &mut edge_exprs,
                        edge_id,
                        component.sample,
                    )?;
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
                            denominator_edges: variant.denominator_edges.clone(),
                            denominator_surface_signs: variant.denominator_surface_signs.clone(),
                            denominator_edge_support_signs: variant
                                .denominator_edge_support_signs
                                .clone(),
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
        weight_polys: &BTreeMap<i32, Vec<Rational>>,
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
            let components = self.finite_pole_contact_components_from_weights(
                edge_id,
                weight_polys,
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
                let base_denominator_edges = variant
                    .denominator_edges
                    .iter()
                    .map(|edge| EdgeIndex(edge_map.get(&edge.0).copied().unwrap_or(edge.0)))
                    .chain(std::iter::once(EdgeIndex(edge_id)))
                    .sorted()
                    .dedup()
                    .collect::<Vec<_>>();

                for component in &components {
                    let mut edge_exprs = base_edge_exprs.clone();
                    self.assign_repeated_channel_sample(
                        &mut edge_exprs,
                        edge_id,
                        component.sample,
                    )?;
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
                            denominator_edges: base_denominator_edges.clone(),
                            denominator_surface_signs: variant.denominator_surface_signs.clone(),
                            denominator_edge_support_signs: variant
                                .denominator_edge_support_signs
                                .clone(),
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

    fn assign_repeated_channel_sample(
        &self,
        edge_exprs: &mut [LinearEnergyExpr],
        edge_id: usize,
        sample: i32,
    ) -> Result<()> {
        let channel = KnownFactorCffBuilder::logical_channels(self.parsed)
            .into_iter()
            .find(|channel| channel.members.contains(&edge_id));
        let members = channel
            .as_ref()
            .map(|channel| channel.members.as_slice())
            .unwrap_or(std::slice::from_ref(&edge_id));
        let reference = &self.parsed.internal_edges[edge_id].signature;

        // Repeated denominator occurrences are one algebraic q0 channel. A
        // quadratic interpolation sample therefore applies to every
        // occurrence, while the occurrence-local map retains only the routing
        // sign needed by the factorized numerator assignment.
        for member in members {
            let relative_sign = KnownFactorCffBuilder::relative_signature_sign(
                reference,
                &self.parsed.internal_edges[*member].signature,
            )?;
            edge_exprs[*member] = if sample == 0 {
                LinearEnergyExpr::zero()
            } else {
                LinearEnergyExpr::ose(EdgeIndex(*member), i64::from(sample * relative_sign))
            };
        }
        Ok(())
    }

    fn lift_contracted_cff_terms(&mut self, pinched_edges: &[usize]) -> Result<()> {
        let (subparsed, sub_to_orig) = self.project_parsed_edges(pinched_edges);
        if subparsed.internal_edges.is_empty() {
            let &[edge_id] = pinched_edges else {
                return Err(GenerationError::CffHigherEnergyPowerNotImplemented);
            };
            let scalar = LowerSectorCffBuilder::new(&subparsed).build()?;
            return self.append_recursive_contact_terms(
                edge_id,
                &scalar,
                &sub_to_orig,
                &contact_weight_polys(self.bounds[edge_id]),
            );
        }

        let use_terminal_residue_basis =
            subparsed.internal_edges.len() == self.parsed.loop_names.len();
        let sub_expression = if use_terminal_residue_basis {
            generate_simple_residue_basis_expression_from_parsed(&subparsed)?
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
                let remapped_denominator_edges = variant
                    .denominator_edges
                    .iter()
                    .map(|edge_id| {
                        EdgeIndex(edge_map.get(&edge_id.0).copied().unwrap_or(edge_id.0))
                    })
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
                    // A residue-basis lower sector already carries the contact
                    // quotient's algebraic sign; only a pure-CFF recursive
                    // pinch needs the extra residue-orientation sign.
                    if !pinched_edges.is_empty()
                        && !use_terminal_residue_basis
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
                            denominator_edges: remapped_denominator_edges.clone(),
                            denominator_surface_signs: variant.denominator_surface_signs.clone(),
                            denominator_edge_support_signs: variant
                                .denominator_edge_support_signs
                                .clone(),
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
        self.finite_pole_contact_components_from_weights(
            edge_id,
            &contact_weight_polys(bound),
            current_expr,
        )
    }

    fn finite_pole_contact_components_from_weights(
        &mut self,
        edge_id: usize,
        weight_polys: &BTreeMap<i32, Vec<Rational>>,
        current_expr: LinearEnergyExpr,
    ) -> Vec<ContactComponent> {
        let q_surface = (!current_expr.is_zero())
            .then(|| self.intern_surface(current_expr, SurfaceOrigin::Helper, true));
        let mut components = Vec::new();
        for (sample, poly) in weight_polys {
            for (power, coeff) in poly.iter().enumerate() {
                if coeff.is_zero() || (power > 0 && q_surface.is_none()) {
                    continue;
                }
                let prefactor = coeff.clone() * rational_pow_i64(2, power + 2);
                if prefactor.is_zero() {
                    continue;
                }
                components.push(ContactComponent {
                    sample: *sample,
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

    fn project_parsed_edges(&self, removed_edges: &[usize]) -> (ParsedGraph, Vec<usize>) {
        let removed = removed_edges.iter().copied().collect::<BTreeSet<_>>();
        let mut parent = self
            .parsed
            .node_name_to_internal
            .values()
            .copied()
            .map(|node| (node, node))
            .collect::<BTreeMap<_, _>>();

        for edge in &self.parsed.internal_edges {
            if removed.contains(&edge.edge_id) {
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
            if removed.contains(&edge.edge_id) {
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
}

struct KnownFactorCffBuilder<'a> {
    original: &'a ParsedGraph,
    bounds: Vec<usize>,
    sampling_scale_mode: NumeratorSamplingScaleMode,
    contact_only: bool,
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
            contact_only: false,
            expression: ThreeDExpression::new_empty(),
            surface_index: HashMap::new(),
        }
    }

    fn build(mut self) -> Result<ThreeDExpression<OrientationID>> {
        let local_to_orig = (0..self.original.internal_edges.len()).collect::<Vec<_>>();
        let recursion_budget = self.recursion_budget();
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
            recursion_budget,
        )?;
        self.finalize_numerator_map_labels();
        Ok(self.expression)
    }

    fn recursion_budget(&self) -> usize {
        2 * self.original.internal_edges.len()
            + self.bounds.iter().copied().sum::<usize>()
            + Self::logical_channels(self.original)
                .iter()
                .map(|channel| channel.members.len() * channel.power)
                .sum::<usize>()
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
        recursion_budget: usize,
    ) -> Result<()> {
        if depth > recursion_budget {
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
                recursion_budget,
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
                        .project_parsed_edges(&delete);
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
                    let relative_sign = Self::relative_signature_sign(
                        &rep_signature,
                        &parsed.internal_edges[*local_id].signature,
                    )?;
                    sub_replacements.insert(
                        orig_id,
                        Self::channel_replacement_expr(
                            orig_id,
                            *sample,
                            relative_sign,
                            use_uniform_scale,
                        ),
                    );
                }

                let mut next_known = known_factors.to_vec();
                if term.parity != 0 {
                    next_known.push(self.known_signature_expr(parsed, &rep_signature)?);
                }
                if term.positive_ose_power != 0 {
                    next_known.extend(
                        (0..term.positive_ose_power).map(|_| KnownLinearExpr::ose(rep_orig, 1)),
                    );
                }
                if term.cancelled_power != 0 {
                    let y_expr = self.known_signature_expr(parsed, &rep_signature)?;
                    let plus = (y_expr.clone() + KnownLinearExpr::ose(rep_orig, 1)).canonical();
                    let minus = (y_expr - KnownLinearExpr::ose(rep_orig, 1)).canonical();
                    for _ in 0..term.cancelled_power {
                        next_known.push(minus.clone());
                        next_known.push(plus.clone());
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

                match next_known
                    .iter()
                    .map(|factor| self.remap_known_factor_to_sub(parsed, &subparsed, factor))
                    .collect::<Result<Vec<_>>>()
                {
                    Ok(sub_known) => self.channel_recursive_terms(
                        &subparsed,
                        &sub_local_to_orig,
                        &sub_replacements,
                        &sub_known,
                        &next_half_edges,
                        extra_uniform_scale_power + channel_uniform_power,
                        channel_prefactor,
                        lower_sector_base || !delete.is_empty(),
                        depth + 1,
                        recursion_budget,
                    )?,
                    Err(GenerationError::SingularBasis)
                        if !sub_local_to_orig.iter().any(|orig_id| {
                            self.bounds[*orig_id] > 1 && !sub_replacements.contains_key(orig_id)
                        }) =>
                    {
                        self.append_generalized_contact_terms(
                            parsed,
                            &subparsed,
                            &sub_local_to_orig,
                            &sub_replacements,
                            &next_known,
                            &next_half_edges,
                            extra_uniform_scale_power + channel_uniform_power,
                            channel_prefactor,
                        )?;
                    }
                    Err(error) => return Err(error),
                }
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
        recursion_budget: usize,
    ) -> Result<()> {
        if depth > recursion_budget {
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
        // An unresolved black-box numerator bound requires another
        // interpolation step. Exact factors introduced by an earlier step do
        // so only on a genuinely powered denominator channel, where their
        // derivatives enter the raised residue. On ordinary edges they stay
        // factorized on the lower-sector CFF; interpreting their basis degrees
        // as fresh bounds would pinch the same cancellation a second time.
        let repeated_edges = Self::logical_channels(parsed)
            .into_iter()
            .flat_map(|channel| channel.members)
            .collect::<BTreeSet<_>>();
        let active = local_to_orig
            .iter()
            .enumerate()
            .find_map(|(local_id, orig_id)| {
                (self.bounds[*orig_id] > 1 && !replacements.contains_key(orig_id))
                    .then_some(local_id)
            })
            .or_else(|| {
                total_bounds
                    .iter()
                    .enumerate()
                    .find_map(|(local_id, degree)| {
                        (*degree > 1 && repeated_edges.contains(&local_id)).then_some(local_id)
                    })
            });
        let Some(active) = active else {
            if self.contact_only && !lower_sector_base {
                return Ok(());
            }
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
                recursion_budget,
            )?;
        }

        let (subparsed, sub_to_local) =
            BoundedCffBuilder::for_bounds(parsed, vec![0; parsed.internal_edges.len()])
                .project_parsed_edges(&[active]);
        let active_signature = &parsed.internal_edges[active].signature.loop_signature;
        let unlocalized_variable = (0..parsed.loop_names.len()).find(|variable| {
            active_signature
                .iter()
                .enumerate()
                .all(|(candidate, coefficient)| {
                    if candidate == *variable {
                        coefficient.abs() == 1
                    } else {
                        *coefficient == 0
                    }
                })
                && subparsed
                    .internal_edges
                    .iter()
                    .all(|edge| edge.signature.loop_signature[*variable] == 0)
        });
        if total_bounds[active] > 2
            && let Some(variable) = unlocalized_variable
        {
            return Err(GenerationError::UnlocalizedEnergyContact {
                variable,
                power: total_bounds[active] - 2,
            });
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
                match sampled_factors
                    .iter()
                    .map(|factor| self.remap_known_factor_to_sub(parsed, &subparsed, factor))
                    .collect::<Result<Vec<_>>>()
                {
                    Ok(remapped_factors) => self.recursive_terms(
                        &subparsed,
                        &sub_local_to_orig,
                        &next_replacements,
                        &remapped_factors,
                        &next_half_edges,
                        extra_uniform_scale_power,
                        branch_prefactor,
                        true,
                        depth + 1,
                        recursion_budget,
                    )?,
                    Err(
                        GenerationError::CffHigherEnergyPowerNotImplemented
                        | GenerationError::SingularBasis,
                    ) => {
                        self.append_generalized_contact_terms(
                            parsed,
                            &subparsed,
                            &sub_local_to_orig,
                            &next_replacements,
                            &sampled_factors,
                            &next_half_edges,
                            extra_uniform_scale_power,
                            branch_prefactor,
                        )?;
                    }
                    Err(err) => return Err(err),
                }
            }
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn append_generalized_contact_terms(
        &mut self,
        parent: &ParsedGraph,
        parsed: &ParsedGraph,
        local_to_orig: &[usize],
        replacements: &BTreeMap<usize, LinearEnergyExpr>,
        known_factors: &[KnownLinearExpr],
        extra_half_edges: &[usize],
        extra_uniform_scale_power: usize,
        prefactor: Rational,
    ) -> Result<()> {
        let numerator = KnownPolynomial::product_from_known_factors(parent, known_factors);
        let branches = KnownPolynomialNormalForm::new(parsed, local_to_orig).branches(numerator)?;
        for branch in branches {
            self.append_known_polynomial_branch(
                parsed,
                local_to_orig,
                replacements,
                branch,
                extra_half_edges,
                extra_uniform_scale_power,
                prefactor.clone(),
            )?;
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn append_known_polynomial_branch(
        &mut self,
        parsed: &ParsedGraph,
        local_to_orig: &[usize],
        replacements: &BTreeMap<usize, LinearEnergyExpr>,
        branch: KnownPolynomialBranch,
        extra_half_edges: &[usize],
        extra_uniform_scale_power: usize,
        prefactor: Rational,
    ) -> Result<()> {
        if branch.numerator.is_zero() || branch.denominator_edges.is_empty() {
            return Ok(());
        }
        if !self.branch_localizes_loop_variables(
            local_to_orig,
            &branch.denominator_edges,
            &branch.numerator,
            replacements,
        ) {
            return Err(GenerationError::CffHigherEnergyPowerNotImplemented);
        }

        let denominator_set = branch
            .denominator_edges
            .iter()
            .copied()
            .collect::<BTreeSet<_>>();
        let delete = (0..parsed.internal_edges.len())
            .filter(|edge_id| !denominator_set.contains(edge_id))
            .collect::<Vec<_>>();
        let (branch_parsed, branch_to_local) =
            BoundedCffBuilder::for_bounds(parsed, vec![0; parsed.internal_edges.len()])
                .project_parsed_edges(&delete);
        if branch_parsed.internal_edges.is_empty() {
            return Ok(());
        }
        let branch_local_to_orig = branch_to_local
            .iter()
            .map(|local_id| local_to_orig[*local_id])
            .collect::<Vec<_>>();
        let base_expression = LowerSectorCffBuilder::new(&branch_parsed).build()?;
        let edge_map = branch_local_to_orig
            .iter()
            .enumerate()
            .map(|(local_id, orig_id)| (local_id, *orig_id))
            .collect::<BTreeMap<_, _>>();
        let surface_map = self.copy_expression_surfaces(&base_expression, &edge_map);
        let original_signatures = self
            .original
            .internal_edges
            .iter()
            .map(|edge| edge.signature.clone())
            .collect::<Vec<_>>();
        for orientation in &base_expression.orientations {
            // A lower sector can lose a denominator direction which an earlier
            // generalized residue has already sampled. Reconstruct the full
            // loop-energy map from its surviving poles first, and use sampled
            // replacements only to complete directions those poles no longer
            // span. A contact sample must not displace a physical pole and
            // thereby change the loop map used by the polynomial remainder.
            let mut target_edge_exprs = vec![LinearEnergyExpr::zero(); original_signatures.len()];
            for (local_id, orig_id) in branch_local_to_orig.iter().copied().enumerate() {
                target_edge_exprs[orig_id] = orientation.edge_energy_map[local_id]
                    .clone()
                    .remap_internal_edges(&edge_map);
            }
            let mut candidate_edges = branch_local_to_orig.clone();
            candidate_edges.extend(
                replacements
                    .keys()
                    .filter(|edge_id| !branch_local_to_orig.contains(edge_id))
                    .copied(),
            );
            let basis_edges = LowerSectorCffBuilder::component_basis_edges(
                &original_signatures,
                &candidate_edges,
            );
            for orig_id in basis_edges
                .iter()
                .filter(|edge_id| !branch_local_to_orig.contains(edge_id))
            {
                target_edge_exprs[*orig_id] = replacements[orig_id].clone();
            }
            let loop_exprs = solve_loop_energy_particular_from_target_edge_exprs(
                &original_signatures,
                &basis_edges,
                &target_edge_exprs,
            )?;
            let mut full_edge_exprs = edge_q0_from_loop_exprs(&original_signatures, &loop_exprs);
            for (orig_id, expr) in replacements {
                full_edge_exprs[*orig_id] = expr.clone();
            }

            for variant in &orientation.variants {
                let base_coeff = rational_from_coefficient(&variant.prefactor) * prefactor.clone();
                if base_coeff.is_zero() {
                    continue;
                }
                let base_half_edges = variant
                    .half_edges
                    .iter()
                    .map(|edge| branch_local_to_orig[edge.0])
                    .chain(extra_half_edges.iter().copied())
                    .collect::<Vec<_>>();
                let base_denominator_edges = variant
                    .denominator_edges
                    .iter()
                    .map(|edge| EdgeIndex(branch_local_to_orig[edge.0]))
                    .collect::<Vec<_>>();
                let base_numerator_surfaces = variant
                    .numerator_surfaces
                    .iter()
                    .map(|surface_id| map_surface_id(*surface_id, &surface_map))
                    .collect::<Vec<_>>();
                let denominator = variant
                    .denominator
                    .clone()
                    .map(|surface_id| map_surface_id(surface_id, &surface_map));

                for (monomial, monomial_coeff) in &branch.numerator.terms {
                    let coeff = base_coeff.clone() * monomial_coeff.clone();
                    if coeff.is_zero() {
                        continue;
                    }
                    let Some(poly_surfaces) = monomial.numerator_surface_exprs(&loop_exprs)? else {
                        continue;
                    };
                    let mut numerator_surfaces = base_numerator_surfaces.clone();
                    for surface_expr in poly_surfaces {
                        numerator_surfaces.push(self.intern_surface(
                            surface_expr,
                            SurfaceOrigin::Helper,
                            true,
                        ));
                    }
                    let mut half_edges = base_half_edges.clone();
                    half_edges.sort_unstable();
                    self.push_variant_for_maps(
                        loop_exprs.clone(),
                        full_edge_exprs.clone(),
                        crate::expression::CFFVariant {
                            origin: Some(
                                if self.contact_only {
                                    "bounded_degree_known_factor_cff_contact_generalized"
                                } else {
                                    "bounded_degree_known_factor_cff_generalized_contact"
                                }
                                .to_string(),
                            ),
                            prefactor: rational_to_coefficient(coeff)?,
                            half_edges: half_edges.into_iter().map(EdgeIndex).collect(),
                            denominator_edges: base_denominator_edges.clone(),
                            denominator_surface_signs: variant.denominator_surface_signs.clone(),
                            denominator_edge_support_signs: variant
                                .denominator_edge_support_signs
                                .clone(),
                            uniform_scale_power: variant.uniform_scale_power
                                + extra_uniform_scale_power,
                            numerator_surfaces,
                            denominator: denominator.clone(),
                        },
                    );
                }
            }
        }
        Ok(())
    }

    fn branch_localizes_loop_variables(
        &self,
        local_to_orig: &[usize],
        denominator_edges: &[usize],
        numerator: &KnownPolynomial,
        replacements: &BTreeMap<usize, LinearEnergyExpr>,
    ) -> bool {
        let loop_variables = numerator.loop_variables();
        if loop_variables.is_empty() {
            return true;
        }
        let lower_sector = LowerSectorCffBuilder::new(self.original);
        let signatures = self
            .original
            .internal_edges
            .iter()
            .map(|edge| edge.signature.clone())
            .collect::<Vec<_>>();
        let candidate_edges = denominator_edges
            .iter()
            .map(|edge_id| local_to_orig[*edge_id])
            .chain(replacements.keys().copied())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let basis_edges =
            LowerSectorCffBuilder::component_basis_edges(&signatures, &candidate_edges);
        let basis_rows = basis_edges
            .iter()
            .map(|edge_id| signatures[*edge_id].loop_signature.clone())
            .collect::<Vec<_>>();
        loop_variables.into_iter().all(|loop_id| {
            let mut unit = vec![0; self.original.loop_names.len()];
            unit[loop_id] = 1;
            lower_sector
                .row_coordinates_in_basis(&basis_rows, &unit)
                .is_ok()
        })
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
        } else if !known_factors.is_empty() || !replacements.is_empty() {
            generate_pure_cff_expression_from_parsed_with_duplicate_excess(
                parsed,
                cff_duplicate_signature_excess(parsed),
            )?
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
                let denominator_edges = variant
                    .denominator_edges
                    .iter()
                    .map(|edge| EdgeIndex(local_to_orig[edge.0]))
                    .collect::<Vec<_>>();
                let denominator = variant
                    .denominator
                    .clone()
                    .map(|surface_id| map_surface_id(surface_id, &surface_map));
                self.push_variant_for_maps(
                    loop_exprs.clone(),
                    full_edge_exprs.clone(),
                    crate::expression::CFFVariant {
                        origin: Some(
                            if self.contact_only {
                                "bounded_degree_known_factor_cff_contact"
                            } else {
                                "bounded_degree_known_factor_cff"
                            }
                            .to_string(),
                        ),
                        prefactor: rational_to_coefficient(coeff)?,
                        half_edges: half_edges.into_iter().map(EdgeIndex).collect(),
                        denominator_edges,
                        denominator_surface_signs: variant.denominator_surface_signs.clone(),
                        denominator_edge_support_signs: variant
                            .denominator_edge_support_signs
                            .clone(),
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

    fn logical_channels(parsed: &ParsedGraph) -> Vec<LogicalChannel> {
        repeated_groups(parsed)
            .into_iter()
            .map(|group| LogicalChannel {
                rep_edge: group.edge_ids[0],
                power: group.edge_ids.len(),
                members: group.edge_ids,
            })
            .collect()
    }

    fn active_repeated_channel(
        &self,
        parsed: &ParsedGraph,
        local_to_orig: &[usize],
        replacements: &BTreeMap<usize, LinearEnergyExpr>,
    ) -> Option<(LogicalChannel, usize)> {
        Self::logical_channels(parsed)
            .into_iter()
            .find_map(|channel| {
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
    }

    fn channel_replacement_expr(
        edge_id: usize,
        sample: i32,
        relative_sign: i32,
        use_uniform_scale: bool,
    ) -> LinearEnergyExpr {
        let sample = sample * relative_sign;
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

    // Divide the sampled numerator polynomial by the powered physical
    // denominator before contracting any occurrence. This keeps cancellations
    // between repeated-pole orders visible to the causal completion.
    fn channel_normal_form_terms(
        polynomial: &[Rational],
        channel_power: usize,
    ) -> Vec<ChannelNormalFormTerm> {
        let mut combined = BTreeMap::<(usize, usize, usize, usize), Rational>::new();
        for (power, coefficient) in polynomial.iter().cloned().enumerate() {
            if coefficient.is_zero() {
                continue;
            }
            let quotient_power = power / 2;
            let parity = power % 2;
            for denominator_power in 0..=quotient_power {
                let coefficient = coefficient.clone() * binomial(quotient_power, denominator_power);
                let remaining = channel_power.saturating_sub(denominator_power);
                let cancelled = denominator_power.saturating_sub(channel_power);
                let inverse_power = 2 * denominator_power + parity;
                let entry = combined
                    .entry((remaining, parity, cancelled, inverse_power))
                    .or_insert_with(Rational::zero);
                *entry = entry.clone() + coefficient;
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
        polynomial: &[Rational],
        channel_power: usize,
    ) -> Vec<ChannelNormalFormTerm> {
        let mut combined = BTreeMap::<(usize, usize, usize, usize, usize), Rational>::new();
        for (power, coefficient) in polynomial.iter().cloned().enumerate() {
            if coefficient.is_zero() {
                continue;
            }
            let quotient_power = power / 2;
            let parity = power % 2;
            for denominator_power in 0..=quotient_power {
                let coefficient = coefficient.clone() * binomial(quotient_power, denominator_power);
                let remaining = channel_power.saturating_sub(denominator_power);
                let cancelled = denominator_power.saturating_sub(channel_power);
                let positive_ose_power = 2 * (quotient_power - denominator_power);
                let entry = combined
                    .entry((remaining, parity, cancelled, power, positive_ose_power))
                    .or_insert_with(Rational::zero);
                *entry = entry.clone() + coefficient;
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
        let basis = LowerSectorCffBuilder::component_basis_edges(&signatures, &edges);
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

fn map_edge_support_signs(
    signs: &BTreeMap<Vec<EdgeIndex>, i64>,
    edge_map: &BTreeMap<usize, usize>,
) -> BTreeMap<Vec<EdgeIndex>, i64> {
    signs
        .iter()
        .map(|(support_edges, sign)| {
            let mut mapped_support = support_edges
                .iter()
                .map(|edge_id| EdgeIndex(edge_map[&edge_id.0]))
                .collect::<Vec<_>>();
            mapped_support.sort_unstable();
            mapped_support.dedup();
            (mapped_support, *sign)
        })
        .collect()
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
    denominator_edges: Vec<usize>,
    chain: Vec<HybridSurfaceID>,
    numerator_surfaces: Vec<HybridSurfaceID>,
    targets: BTreeMap<usize, LinearEnergyExpr>,
    edge_exprs: BTreeMap<usize, LinearEnergyExpr>,
}

struct LowerSectorCffBuilder<'a> {
    parsed: &'a ParsedGraph,
    force_component_factorization: bool,
    expression: ThreeDExpression<OrientationID>,
    surface_index: HashMap<(LinearSurfaceKind, LinearEnergyExpr), HybridSurfaceID>,
}

impl<'a> LowerSectorCffBuilder<'a> {
    fn new(parsed: &'a ParsedGraph) -> Self {
        Self {
            parsed,
            force_component_factorization: false,
            expression: ThreeDExpression::new_empty(),
            surface_index: HashMap::new(),
        }
    }

    fn build(mut self) -> Result<ThreeDExpression<OrientationID>> {
        if self.parsed.internal_edges.is_empty() {
            // A fully contracted scalar contact is the zero-edge lower sector,
            // whose denominator is the multiplicative identity.
            let mut data = OrientationData::new(EdgeVec::new());
            data.label = Some("lower_sector_unit".to_string());
            let mut orientation = OrientationExpression::lower_sector_unit(data);
            orientation.loop_energy_map =
                vec![LinearEnergyExpr::zero(); self.parsed.loop_names.len()];
            self.expression.orientations.push(orientation);
            self.finalize_numerator_map_labels();
            return Ok(self.expression);
        }
        let signatures = self.signatures();
        let denominator_edge_ids = self.parsed.denominator_internal_edge_ids();
        let denominator_rows = denominator_edge_ids
            .iter()
            .map(|edge_id| {
                let signature = &signatures[*edge_id];
                signature
                    .loop_signature
                    .iter()
                    .map(|value| i64::from(*value))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let denominator_rank = rank_i64(&denominator_rows);
        if denominator_rank == 0 {
            return Err(GenerationError::CffHigherEnergyPowerNotImplemented);
        }
        if denominator_rank == self.parsed.loop_names.len()
            && denominator_connected_components(self.parsed).len() == 1
            && repeated_groups(self.parsed).is_empty()
            && !self.force_component_factorization
        {
            // A connected, unraised full-rank lower sector is already an
            // ordinary affine CFF problem. Keep it intact: vector-matroid
            // projection would turn loop energy carried by another signature
            // component into incomplete external boundary data. Disconnected
            // sectors and powered channels stay on the component path below
            // for their shared core sign and numerator derivatives.
            return generate_pure_cff_expression_from_parsed(self.parsed);
        }
        let components = self.component_bundles(&signatures)?;
        // A component product needs the relative core sign (-1)^(C-1): for a
        // particular lower-sector lift C counts genuinely disconnected
        // denominator factors, while a top-level factorization counts its
        // independent signature factors. Signature bundles used only as a
        // routing device inside one connected particular lift do not count.
        // An initial-state cut retains the already-fixed shared LU convention.
        let denominator_component_count = if self.force_component_factorization {
            components.len()
        } else {
            denominator_connected_components(self.parsed).len()
        };
        let component_product_sign = if self.parsed.initial_state_cut_edges.is_empty()
            && denominator_component_count % 2 == 0
        {
            Rational::from(-1)
        } else {
            Rational::one()
        };
        let mut partials = vec![LowerSectorPartial {
            coeff: component_product_sign,
            half_edges: Vec::new(),
            denominator_edges: Vec::new(),
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
                        item.denominator_edges.extend(
                            variant
                                .denominator_edges
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
            apply_initial_state_cut_edge_energy_exprs(self.parsed, &mut edge_exprs);
            let mut half_edges = partial.half_edges;
            half_edges.sort_unstable();
            self.push_variant_for_maps(
                loop_exprs,
                edge_exprs,
                crate::expression::CFFVariant {
                    origin: Some("lower_sector_cff_e_surface_component_product".to_string()),
                    prefactor: rational_to_coefficient(partial.coeff)?,
                    half_edges: half_edges.into_iter().map(EdgeIndex).collect(),
                    denominator_edges: partial
                        .denominator_edges
                        .into_iter()
                        .map(EdgeIndex)
                        .collect(),
                    denominator_surface_signs: BTreeMap::new(),
                    denominator_edge_support_signs: BTreeMap::new(),
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
        let basis_edges = Self::component_basis_edges(signatures, &edges);
        let (component_parsed, local_to_sub) =
            self.project_component_parsed(signatures, &edges, &basis_edges)?;
        // A basis-only bundle inside one connected affine graph supplies one
        // particular lower-sector contour closure used to reconstruct the
        // parent loop-energy map. A genuinely disconnected factor, like an
        // auxiliary-denominator bundle, instead requires its complete explicit
        // CFF orientation sum. That CFF omits its duplicate-factor sign here so
        // the product lift can restore the lower-sector parity exactly once.
        let genuinely_disconnected = denominator_connected_components(self.parsed).len() > 1;
        let expression = if component_parsed.internal_edges.len() == rank && !genuinely_disconnected
        {
            generate_simple_residue_basis_expression_from_parsed(&component_parsed)?
        } else {
            generate_pure_cff_expression_from_parsed_with_duplicate_sign(&component_parsed, false)?
        };
        // The component CFF deliberately omits its duplicate-denominator sign.
        // Restore precisely that parity for uncut factors: rank deficiency can
        // also come from distinct affine propagators and is not itself a sign.
        // A cut particular lift retains the established auxiliary-contour
        // parity, which already participates in the LU cut convention.
        let is_auxiliary = component_parsed.internal_edges.len() != rank;
        let correction_is_negative = if self.parsed.initial_state_cut_edges.is_empty() {
            cff_duplicate_signature_excess(&component_parsed) % 2 == 1
        } else {
            is_auxiliary && (component_parsed.internal_edges.len() - rank + 1) % 2 == 1
        };
        let prefactor_correction = if correction_is_negative {
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

        let component_nodes = edges
            .iter()
            .flat_map(|edge_id| {
                let edge = &self.parsed.internal_edges[*edge_id];
                [edge.tail, edge.head]
            })
            .collect::<BTreeSet<_>>();
        let old_to_new = component_nodes
            .iter()
            .copied()
            .enumerate()
            .map(|(new_id, old_id)| (old_id, new_id))
            .collect::<BTreeMap<_, _>>();
        let mut local_to_sub = Vec::new();
        let mut internal_edges = Vec::new();
        for (new_id, (sub_id, loop_coeffs)) in edges.iter().copied().zip(projected).enumerate() {
            let original = &self.parsed.internal_edges[sub_id];
            local_to_sub.push(sub_id);
            internal_edges.push(ParsedGraphInternalEdge {
                edge_id: new_id,
                tail: old_to_new[&original.tail],
                head: old_to_new[&original.head],
                label: original.label.clone(),
                mass_key: original.mass_key.clone(),
                signature: MomentumSignature {
                    loop_signature: loop_coeffs,
                    external_signature: original.signature.external_signature.clone(),
                },
                had_pow: original.had_pow,
            });
        }
        let orig_to_local = local_to_sub
            .iter()
            .enumerate()
            .map(|(local_id, sub_id)| (*sub_id, local_id))
            .collect::<BTreeMap<_, _>>();
        let mut next_external_edge_id = 0;
        let mut external_edges = self
            .parsed
            .external_edges
            .iter()
            .filter_map(|edge| {
                let source = edge
                    .source
                    .and_then(|source| old_to_new.get(&source).copied());
                let destination = edge
                    .destination
                    .and_then(|destination| old_to_new.get(&destination).copied());
                (source.is_some() || destination.is_some()).then(|| {
                    let edge_id = next_external_edge_id;
                    next_external_edge_id += 1;
                    ParsedGraphExternalEdge {
                        edge_id,
                        source,
                        destination,
                        label: edge.label.clone(),
                        external_coefficients: edge.external_coefficients.clone(),
                    }
                })
            })
            .collect::<Vec<_>>();
        let component_edge_set = edges.iter().copied().collect::<BTreeSet<_>>();
        for edge in &self.parsed.internal_edges {
            if component_edge_set.contains(&edge.edge_id) {
                continue;
            }
            if let Some(source) = old_to_new.get(&edge.tail).copied() {
                external_edges.push(ParsedGraphExternalEdge {
                    edge_id: next_external_edge_id,
                    source: Some(source),
                    destination: None,
                    label: format!("{}_source_boundary", edge.label),
                    external_coefficients: edge.signature.external_signature.clone(),
                });
                next_external_edge_id += 1;
            }
            if let Some(destination) = old_to_new.get(&edge.head).copied() {
                external_edges.push(ParsedGraphExternalEdge {
                    edge_id: next_external_edge_id,
                    source: None,
                    destination: Some(destination),
                    label: format!("{}_sink_boundary", edge.label),
                    external_coefficients: edge.signature.external_signature.clone(),
                });
                next_external_edge_id += 1;
            }
        }
        let initial_state_cut_edges = self
            .parsed
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
            .collect::<Vec<_>>();
        let node_name_to_internal = old_to_new
            .iter()
            .map(|(old_id, new_id)| (format!("lower_{old_id}"), *new_id))
            .collect::<BTreeMap<_, _>>();
        Ok((
            ParsedGraph {
                internal_edges,
                external_edges,
                initial_state_cut_edges,
                loop_names,
                external_names: self.parsed.external_names.clone(),
                node_name_to_internal,
            },
            local_to_sub,
        ))
    }

    fn vector_matroid_components(&self, signatures: &[MomentumSignature]) -> Vec<Vec<usize>> {
        let rows = signatures
            .iter()
            .map(|signature| signature.loop_signature.clone())
            .collect::<Vec<_>>();
        let denominator_edges = self
            .parsed
            .denominator_internal_edge_ids()
            .into_iter()
            .collect::<BTreeSet<_>>();
        let nonzero = rows
            .iter()
            .enumerate()
            .filter_map(|(edge_id, row)| {
                (denominator_edges.contains(&edge_id) && row.iter().any(|value| *value != 0))
                    .then_some(edge_id)
            })
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
            if let Ok(coords) = self.rational_row_coordinates_in_basis(&basis_rows, &rows[*edge_id])
            {
                for (basis_edge, coeff) in basis.iter().zip(coords) {
                    if !coeff.is_zero() {
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
            (denominator_edges.contains(&edge_id) && !row.iter().any(|value| *value != 0))
                .then_some(vec![edge_id])
        }));
        components.sort_by_key(|component| component[0]);
        components
    }

    fn component_basis_edges(signatures: &[MomentumSignature], edges: &[usize]) -> Vec<usize> {
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
        let coords = self.rational_row_coordinates_in_basis(basis_rows, row)?;
        if coords.iter().any(|coord| {
            coord
                .to_i64_pair()
                .is_none_or(|(_, denominator)| denominator != 1)
        }) {
            return Err(GenerationError::NonIntegralEnergyMap);
        }
        coords
            .into_iter()
            .map(|coord| {
                let (numerator, _) = coord
                    .to_i64_pair()
                    .ok_or(GenerationError::CoefficientOutOfRange)?;
                i32::try_from(numerator).map_err(|_| GenerationError::CoefficientOutOfRange)
            })
            .collect()
    }

    fn rational_row_coordinates_in_basis(
        &self,
        basis_rows: &[Vec<i32>],
        row: &[i32],
    ) -> Result<Vec<Rational>> {
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
            let reconstructs_row = (0..row.len()).all(|column| {
                basis_rows.iter().zip(&coords).fold(
                    Rational::zero(),
                    |value, (basis_row, coord)| {
                        value + Rational::from(basis_row[column]) * coord.clone()
                    },
                ) == row[column]
            });
            if reconstructs_row {
                return Ok(coords);
            }
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
            insert_terminal_unit_if_missing(&mut tree, NodeId::root());
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
        insert_terminal_unit_if_missing(&mut tree, parent);
    }
    tree
}

fn insert_terminal_unit_if_missing(tree: &mut Tree<HybridSurfaceID>, parent: NodeId) {
    let has_terminal_unit = tree
        .get_node(parent)
        .children
        .iter()
        .any(|child| tree.get_node(*child).data == HybridSurfaceID::Unit);
    if !has_terminal_unit {
        tree.insert_node(parent, HybridSurfaceID::Unit);
    }
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
    let mut counts = BTreeMap::<(MomentumSignature, Option<String>), usize>::new();
    for edge in &parsed.internal_edges {
        if parsed.is_initial_state_cut_edge(edge.edge_id) {
            continue;
        }
        let (signature, _) = edge.signature.canonical_up_to_sign();
        *counts
            .entry((signature, edge.mass_key.clone()))
            .or_default() += 1;
    }
    counts.values().map(|count| count.saturating_sub(1)).sum()
}

fn has_duplicate_signature_ignoring_mass(parsed: &ParsedGraph) -> bool {
    let mut counts = BTreeMap::<MomentumSignature, usize>::new();
    for edge in &parsed.internal_edges {
        if parsed.is_initial_state_cut_edge(edge.edge_id) {
            continue;
        }
        let (signature, _) = edge.signature.canonical_up_to_sign();
        *counts.entry(signature).or_default() += 1;
    }
    counts.values().any(|count| *count > 1)
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

fn cff_bounds_need_generalized_expression_from_options(
    parsed: &ParsedGraph,
    options: &Generate3DExpressionOptions,
) -> Result<bool> {
    let bounds = normalize_energy_degree_bounds(
        options.energy_degree_bounds.as_deref().unwrap_or(&[]),
        parsed.internal_edges.len(),
    )?;
    Ok(cff_bounds_need_generalized_expression(parsed, &bounds))
}

fn cff_bounds_need_generalized_expression(parsed: &ParsedGraph, bounds: &[usize]) -> bool {
    if bounds.iter().all(|degree| *degree == 0) {
        return false;
    }
    bounds.iter().any(|degree| *degree > 1)
        || KnownFactorCffBuilder::logical_channels(parsed)
            .iter()
            .any(|channel| {
                channel
                    .members
                    .iter()
                    .map(|edge_id| bounds[*edge_id])
                    .sum::<usize>()
                    > 1
            })
}

#[cfg(test)]
mod representation_tests {
    use super::*;

    #[test]
    fn ltd_setting_returns_explicit_not_implemented_error() {
        let error = generate_3d_expression_from_parsed(
            &crate::graph_io::test_graphs::box_graph(),
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Ltd,
                ..Default::default()
            },
        )
        .unwrap_err();

        assert!(matches!(
            error,
            GenerationError::NotImplemented {
                mode: RepresentationMode::Ltd
            }
        ));
    }
}

#[cfg(test)]
mod causal_generation_tests {
    use super::*;

    #[cfg(feature = "eval")]
    #[test]
    fn factorized_lower_sector_restores_global_core_sign_and_orientation_sum() {
        let parsed = ParsedGraph {
            internal_edges: vec![
                ParsedGraphInternalEdge {
                    edge_id: 0,
                    tail: 0,
                    head: 0,
                    label: "q0".to_string(),
                    mass_key: Some("m0".to_string()),
                    signature: MomentumSignature {
                        loop_signature: vec![1, 0],
                        external_signature: Vec::new(),
                    },
                    had_pow: false,
                },
                ParsedGraphInternalEdge {
                    edge_id: 1,
                    tail: 1,
                    head: 1,
                    label: "q1".to_string(),
                    mass_key: Some("m1".to_string()),
                    signature: MomentumSignature {
                        loop_signature: vec![0, 1],
                        external_signature: Vec::new(),
                    },
                    had_pow: false,
                },
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["q0".to_string(), "q1".to_string()],
            external_names: Vec::new(),
            node_name_to_internal: BTreeMap::from([("n0".to_string(), 0), ("n1".to_string(), 1)]),
        };
        let lower = LowerSectorCffBuilder::new(&parsed).build().unwrap();
        let input = crate::eval::EvaluationInput {
            external_momenta: Vec::new(),
            loop_spatial_momenta: vec![[0.0; 3]; 2],
            masses: vec![0.7, 1.1],
            uniform_scale: None,
        };
        let lower_value = crate::eval::evaluate_expression(&parsed, &lower, "1", &input)
            .unwrap()
            .value;
        // Each one-loop component contributes its two explicit orientations,
        // 1/E_i. The shared two-loop core supplies (-1)^(2-1).
        let expected = -1.0 / (0.7 * 1.1);

        assert!(
            (lower_value - expected).abs() <= 1.0e-13 * expected.abs(),
            "factorized lower-sector CFF changed the shared core sign or orientation sum: lower={lower_value:e}, expected={expected:e}"
        );
    }

    #[test]
    fn finite_pole_sampling_changes_only_the_selected_edge_map() {
        let parsed = crate::graph_io::test_graphs::box_pow3_graph();
        let source = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(0, 1), (1, 1), (2, 1), (3, 2)]),
                ..Default::default()
            },
        )
        .unwrap();
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(0, 2), (1, 1), (2, 1), (3, 2)]),
                ..Default::default()
            },
        )
        .unwrap();

        let remainder_orientations = expression
            .orientations
            .iter()
            .filter(|orientation| {
                orientation.variants.iter().any(|variant| {
                    variant.origin.as_deref().is_some_and(|origin| {
                        origin.starts_with("bounded_degree_quadratic_recursive_remainder:e0=")
                    })
                })
            })
            .collect::<Vec<_>>();
        assert!(!remainder_orientations.is_empty());

        for orientation in remainder_orientations {
            assert!(
                orientation.edge_energy_map[0] == LinearEnergyExpr::ose(EdgeIndex(0), 1)
                    || orientation.edge_energy_map[0] == LinearEnergyExpr::ose(EdgeIndex(0), -1)
            );
            assert!(source.orientations.iter().any(|source_orientation| {
                source_orientation.loop_energy_map == orientation.loop_energy_map
                    && source_orientation
                        .edge_energy_map
                        .iter()
                        .zip(&orientation.edge_energy_map)
                        .enumerate()
                        .all(|(edge_id, (source, sampled))| edge_id == 0 || source == sampled)
            }));
        }
    }

    #[test]
    fn lower_sector_joins_rationally_dependent_rows_but_rejects_nonintegral_projection() {
        let parsed = ParsedGraph {
            internal_edges: vec![
                ParsedGraphInternalEdge {
                    edge_id: 0,
                    tail: 0,
                    head: 1,
                    label: "q0".to_string(),
                    mass_key: Some("m0".to_string()),
                    signature: MomentumSignature {
                        loop_signature: vec![2, 0],
                        external_signature: Vec::new(),
                    },
                    had_pow: false,
                },
                ParsedGraphInternalEdge {
                    edge_id: 1,
                    tail: 1,
                    head: 0,
                    label: "q1".to_string(),
                    mass_key: Some("m1".to_string()),
                    signature: MomentumSignature {
                        loop_signature: vec![1, 0],
                        external_signature: Vec::new(),
                    },
                    had_pow: false,
                },
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["k0".to_string(), "k1".to_string()],
            external_names: Vec::new(),
            node_name_to_internal: BTreeMap::from([("n0".to_string(), 0), ("n1".to_string(), 1)]),
        };
        let lower_sector = LowerSectorCffBuilder::new(&parsed);

        assert_eq!(
            lower_sector.vector_matroid_components(&lower_sector.signatures()),
            vec![vec![0, 1]]
        );
        // Component connectivity is a rational-vector-matroid property, while
        // production signatures still require an integral projected lattice.
        assert!(matches!(
            lower_sector.build(),
            Err(GenerationError::NonIntegralEnergyMap)
        ));
    }

    #[test]
    fn bounded_cff_reconstructs_sampled_directions_in_theta_contacts() {
        let incidences = [
            (0, 2, [1, 0], "m0"),
            (2, 1, [1, 0], "m0"),
            (0, 3, [0, 1], "m1"),
            (3, 1, [0, 1], "m1"),
            (0, 1, [-1, -1], "m2"),
        ];
        let parsed = ParsedGraph {
            internal_edges: incidences
                .into_iter()
                .enumerate()
                .map(
                    |(edge_id, (tail, head, loop_signature, mass))| ParsedGraphInternalEdge {
                        edge_id,
                        tail,
                        head,
                        label: format!("q{edge_id}"),
                        mass_key: Some(mass.to_string()),
                        signature: MomentumSignature {
                            loop_signature: loop_signature.to_vec(),
                            external_signature: Vec::new(),
                        },
                        had_pow: false,
                    },
                )
                .collect(),
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["k0".to_string(), "k1".to_string()],
            external_names: Vec::new(),
            node_name_to_internal: (0..=3).map(|node| (format!("n{node}"), node)).collect(),
        };
        let expression = BoundedCffBuilder::new(
            &parsed,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(0, 3), (2, 3)]),
                ..Default::default()
            },
        )
        .unwrap()
        .build()
        .unwrap();

        assert!(!expression.orientations.is_empty());
    }

    #[test]
    fn top_level_cff_rejects_rank_deficiency_even_when_cut_rows_span_it() {
        let mut parsed = crate::graph_io::test_graphs::box_graph();
        parsed.loop_names.push("k_cut".to_string());
        for edge in &mut parsed.internal_edges {
            edge.signature.loop_signature.push(0);
        }
        let cut_edge_id = parsed.internal_edges.len();
        parsed.internal_edges.push(ParsedGraphInternalEdge {
            edge_id: cut_edge_id,
            tail: 0,
            head: 1,
            label: "p_cut_0".to_string(),
            mass_key: Some("m_cut_0".to_string()),
            signature: MomentumSignature {
                // These stored rows are deliberately outside the active
                // denominator span and must not add residue rank.
                loop_signature: vec![0, 1],
                external_signature: vec![1, 0, 0],
            },
            had_pow: false,
        });
        parsed.internal_edges.push(ParsedGraphInternalEdge {
            edge_id: cut_edge_id + 1,
            tail: 1,
            head: 0,
            label: "p_cut_1".to_string(),
            mass_key: Some("m_cut_1".to_string()),
            signature: MomentumSignature {
                loop_signature: vec![0, 1],
                external_signature: vec![0, 1, 0],
            },
            had_pow: false,
        });
        parsed.initial_state_cut_edges.extend([
            ParsedGraphInitialStateCutEdge {
                edge_id: cut_edge_id,
                external_id: 0,
                external_sign: 1,
            },
            ParsedGraphInitialStateCutEdge {
                edge_id: cut_edge_id + 1,
                external_id: 1,
                external_sign: 1,
            },
        ]);

        assert!(crate::validate_parsed_graph(&parsed).ok);
        let error = generate_3d_expression(
            &parsed,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(Vec::new()),
                ..Default::default()
            },
        )
        .unwrap_err();
        assert!(matches!(error, GenerationError::SingularBasis));
    }

    #[test]
    fn lower_sector_particular_lift_keeps_initial_cut_external() {
        let mut parsed = crate::graph_io::test_graphs::box_graph();
        parsed.loop_names.push("k_cut".to_string());
        for edge in &mut parsed.internal_edges {
            edge.signature.loop_signature.push(0);
        }
        let cut_edge_id = parsed.internal_edges.len();
        parsed.internal_edges.push(ParsedGraphInternalEdge {
            edge_id: cut_edge_id,
            tail: 0,
            head: 1,
            label: "p_cut_0".to_string(),
            mass_key: Some("m_cut_0".to_string()),
            signature: MomentumSignature {
                // These stored rows are deliberately outside the active
                // denominator span and must not add residue rank.
                loop_signature: vec![0, 1],
                external_signature: vec![1, 0, 0],
            },
            had_pow: false,
        });
        parsed.internal_edges.push(ParsedGraphInternalEdge {
            edge_id: cut_edge_id + 1,
            tail: 1,
            head: 0,
            label: "p_cut_1".to_string(),
            mass_key: Some("m_cut_1".to_string()),
            signature: MomentumSignature {
                loop_signature: vec![0, 1],
                external_signature: vec![0, 1, 0],
            },
            had_pow: false,
        });
        parsed.initial_state_cut_edges.extend([
            ParsedGraphInitialStateCutEdge {
                edge_id: cut_edge_id,
                external_id: 0,
                external_sign: 1,
            },
            ParsedGraphInitialStateCutEdge {
                edge_id: cut_edge_id + 1,
                external_id: 1,
                external_sign: 1,
            },
        ]);

        assert!(crate::validate_parsed_graph(&parsed).ok);
        let expression = LowerSectorCffBuilder::new(&parsed).build().unwrap();
        assert!(!expression.orientations.is_empty());
        assert!(expression.orientations.iter().all(|orientation| {
            orientation.loop_energy_map[1] == LinearEnergyExpr::zero()
                && orientation.edge_energy_map[cut_edge_id]
                    == LinearEnergyExpr::external(EdgeIndex(0), 1)
                && orientation.edge_energy_map[cut_edge_id + 1]
                    == LinearEnergyExpr::external(EdgeIndex(1), 1)
        }));
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

        let generated = generate_3d_expression(
            &source,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(10, 1), (11, 1)]),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            generated.energy_factor_components,
            vec![CffEnergyFactorComponent {
                internal_edge_ids: vec![10, 11, 12, 13],
                ownership: generated.energy_factor_ownership,
            }]
        );
        let expression = generated.expression;

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

    fn generated_expression(external_sign: i32) -> ThreeDExpression<OrientationID> {
        let parsed = crate::graph_io::test_graphs::initial_state_cut_line_graph(external_sign);
        generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(Vec::new()),
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
            let expression = generated_expression(external_sign);
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
    fn cff_duplicate_sign_correction_only_counts_identical_denominators() {
        let mut parsed = crate::graph_io::test_graphs::box_pow3_graph();
        parsed.internal_edges[4].mass_key = Some("different_mass".to_string());
        assert_eq!(cff_duplicate_signature_excess(&parsed), 1);
    }
}

#[cfg(test)]
mod cff_tests {
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
    fn cff_generation_preserves_external_tree_edges() {
        let parsed = crate::graph_io::test_graphs::triangle_with_external_tree_graph();
        let generated = generate_3d_expression_from_parsed_generated(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Cff,
                preserve_internal_edges_as_four_d_denominators: vec![3],
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            generated.energy_factor_components,
            vec![CffEnergyFactorComponent {
                internal_edge_ids: vec![0, 1, 2],
                ownership: generated.energy_factor_ownership,
            }]
        );
        let expression = generated.expression;

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
    fn cff_preserved_edge_embedding_remaps_component_metadata() {
        let mut parsed = crate::graph_io::test_graphs::triangle_with_external_tree_graph();
        parsed.internal_edges.swap(0, 3);
        for (edge_id, edge) in parsed.internal_edges.iter_mut().enumerate() {
            edge.edge_id = edge_id;
        }

        let generated = generate_3d_expression_from_parsed_generated(
            &parsed,
            &Generate3DExpressionOptions {
                preserve_internal_edges_as_four_d_denominators: vec![0],
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(
            generated.energy_factor_components,
            vec![CffEnergyFactorComponent {
                internal_edge_ids: vec![1, 2, 3],
                ownership: generated.energy_factor_ownership,
            }]
        );
        assert_eq!(
            generated.expression.residual_denominators[0].edge_id,
            EdgeIndex(0)
        );
    }

    #[test]
    fn cff_generation_leaves_pure_trees_untouched() {
        let parsed = crate::graph_io::test_graphs::pure_tree_graph();
        let generated = generate_3d_expression_from_parsed_generated(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Cff,
                preserve_internal_edges_as_four_d_denominators: vec![0, 1],
                ..Default::default()
            },
        )
        .unwrap();
        assert!(generated.energy_factor_components.is_empty());
        let expression = generated.expression;

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
    fn cff_zero_edge_expression_has_no_energy_factor_components() {
        let parsed = ParsedGraph {
            internal_edges: Vec::new(),
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: Vec::new(),
            external_names: Vec::new(),
            node_name_to_internal: BTreeMap::new(),
        };

        let generated = generate_3d_expression_from_parsed_generated(
            &parsed,
            &Generate3DExpressionOptions::default(),
        )
        .unwrap();

        assert!(generated.energy_factor_components.is_empty());
        assert_eq!(
            generated.energy_factor_ownership,
            CffEnergyFactorOwnership::GlobalSourceProduct
        );
    }

    fn penta_contact_graph() -> ParsedGraph {
        ParsedGraph {
            internal_edges: vec![
                ParsedGraphInternalEdge {
                    edge_id: 0,
                    tail: 0,
                    head: 1,
                    label: "q0".to_string(),
                    mass_key: Some("m0".to_string()),
                    signature: MomentumSignature {
                        loop_signature: vec![1, 1],
                        external_signature: vec![0, 0, 1],
                    },
                    had_pow: false,
                },
                ParsedGraphInternalEdge {
                    edge_id: 1,
                    tail: 1,
                    head: 2,
                    label: "q1".to_string(),
                    mass_key: Some("m1".to_string()),
                    signature: MomentumSignature {
                        loop_signature: vec![1, 0],
                        external_signature: vec![0, 0, 0],
                    },
                    had_pow: false,
                },
                ParsedGraphInternalEdge {
                    edge_id: 2,
                    tail: 2,
                    head: 0,
                    label: "q2".to_string(),
                    mass_key: Some("m2".to_string()),
                    signature: MomentumSignature {
                        loop_signature: vec![1, 0],
                        external_signature: vec![1, 0, 0],
                    },
                    had_pow: false,
                },
                ParsedGraphInternalEdge {
                    edge_id: 3,
                    tail: 0,
                    head: 3,
                    label: "q3".to_string(),
                    mass_key: Some("m3".to_string()),
                    signature: MomentumSignature {
                        loop_signature: vec![0, 1],
                        external_signature: vec![0, 0, 0],
                    },
                    had_pow: false,
                },
                ParsedGraphInternalEdge {
                    edge_id: 4,
                    tail: 3,
                    head: 1,
                    label: "q4".to_string(),
                    mass_key: Some("m4".to_string()),
                    signature: MomentumSignature {
                        loop_signature: vec![0, 1],
                        external_signature: vec![0, 1, 0],
                    },
                    had_pow: false,
                },
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["k1".to_string(), "k2".to_string()],
            external_names: vec!["p1".to_string(), "p2".to_string(), "p3".to_string()],
            node_name_to_internal: BTreeMap::from([
                ("n0".to_string(), 0),
                ("n1".to_string(), 1),
                ("n2".to_string(), 2),
                ("n3".to_string(), 3),
            ]),
        }
    }

    fn isolated_vertex_repeated_bubble_source() -> ParsedGraph {
        ParsedGraph {
            internal_edges: vec![
                ParsedGraphInternalEdge {
                    edge_id: 0,
                    tail: 1,
                    head: 2,
                    label: "q0".to_string(),
                    mass_key: Some("m".to_string()),
                    signature: MomentumSignature {
                        loop_signature: vec![1],
                        external_signature: Vec::new(),
                    },
                    had_pow: false,
                },
                ParsedGraphInternalEdge {
                    edge_id: 1,
                    tail: 2,
                    head: 1,
                    label: "q1".to_string(),
                    mass_key: Some("m".to_string()),
                    signature: MomentumSignature {
                        loop_signature: vec![-1],
                        external_signature: Vec::new(),
                    },
                    had_pow: false,
                },
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["k1".to_string()],
            external_names: Vec::new(),
            node_name_to_internal: BTreeMap::from([
                ("n0".to_string(), 0),
                ("n1".to_string(), 1),
                ("n2".to_string(), 2),
            ]),
        }
    }

    #[test]
    fn disconnected_mixed_energy_factor_ownership_keeps_both_component_factors() {
        let edge =
            |edge_id, tail, head, loop_signature: [i32; 2], external_signature, mass: &str| {
                ParsedGraphInternalEdge {
                    edge_id,
                    tail,
                    head,
                    label: format!("q{edge_id}"),
                    mass_key: Some(mass.to_string()),
                    signature: MomentumSignature {
                        loop_signature: loop_signature.to_vec(),
                        external_signature,
                    },
                    had_pow: false,
                }
            };
        let parsed = ParsedGraph {
            internal_edges: vec![
                edge(0, 0, 1, [1, 0], vec![0], "mr"),
                edge(1, 1, 0, [-1, 0], vec![0], "mr"),
                edge(2, 2, 3, [0, 1], vec![0], "m2"),
                edge(3, 3, 2, [0, 1], vec![1], "m3"),
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["ka".to_string(), "kb".to_string()],
            external_names: vec!["p".to_string()],
            node_name_to_internal: BTreeMap::from([
                ("a0".to_string(), 0),
                ("a1".to_string(), 1),
                ("b0".to_string(), 2),
                ("b1".to_string(), 3),
            ]),
        };
        let options = Generate3DExpressionOptions {
            energy_degree_bounds: Some(vec![(0, 2)]),
            ..Default::default()
        };
        let component_ownerships = denominator_connected_components(&parsed)
            .into_iter()
            .map(|component_edges| {
                let (component, embedding) =
                    project_denominator_component(&parsed, &component_edges)?;
                let component_options =
                    project_component_options(&options, parsed.internal_edges.len(), &embedding)?;
                Ok(
                    generate_3d_expression_from_parsed_generated(&component, &component_options)?
                        .energy_factor_ownership,
                )
            })
            .collect::<Result<Vec<_>>>()
            .unwrap();

        assert_eq!(
            component_ownerships,
            vec![
                CffEnergyFactorOwnership::VariantLocal,
                CffEnergyFactorOwnership::GlobalSourceProduct,
            ]
        );

        let generated = generate_3d_expression_from_parsed_generated(&parsed, &options).unwrap();
        assert_eq!(
            generated.energy_factor_ownership,
            CffEnergyFactorOwnership::VariantLocal
        );
        assert_eq!(
            generated.energy_factor_components,
            vec![
                CffEnergyFactorComponent {
                    internal_edge_ids: vec![0, 1],
                    ownership: CffEnergyFactorOwnership::VariantLocal,
                },
                CffEnergyFactorComponent {
                    internal_edge_ids: vec![2, 3],
                    ownership: CffEnergyFactorOwnership::GlobalSourceProduct,
                },
            ],
            "component-local ownership must survive the deterministic disconnected product order"
        );
        let variants = generated
            .expression
            .orientations
            .iter()
            .flat_map(|orientation| &orientation.variants)
            .collect::<Vec<_>>();
        assert!(!variants.is_empty());
        let repeated_component = BTreeSet::from([0, 1]);
        let ordinary_component = BTreeSet::from([2, 3]);
        assert!(variants.into_iter().all(|variant| {
            let half_edges = variant
                .half_edges
                .iter()
                .map(|edge| edge.0)
                .collect::<BTreeSet<_>>();
            !half_edges.is_disjoint(&repeated_component)
                && !half_edges.is_disjoint(&ordinary_component)
        }));
    }

    #[test]
    fn disconnected_metadata_composes_componentwise_global_prefactor_signs() {
        let edge =
            |edge_id, tail, head, loop_signature: [i32; 2], mass: &str| ParsedGraphInternalEdge {
                edge_id,
                tail,
                head,
                label: format!("q{edge_id}"),
                mass_key: Some(mass.to_string()),
                signature: MomentumSignature {
                    loop_signature: loop_signature.to_vec(),
                    external_signature: Vec::new(),
                },
                had_pow: false,
            };
        let parsed = ParsedGraph {
            internal_edges: vec![
                edge(0, 0, 1, [1, 0], "ma"),
                edge(1, 1, 0, [-1, 0], "ma"),
                edge(2, 2, 3, [0, 1], "mb"),
                edge(3, 3, 2, [0, -1], "mb"),
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["ka".to_string(), "kb".to_string()],
            external_names: Vec::new(),
            node_name_to_internal: BTreeMap::from([
                ("a0".to_string(), 0),
                ("a1".to_string(), 1),
                ("b0".to_string(), 2),
                ("b1".to_string(), 3),
            ]),
        };
        let options = Generate3DExpressionOptions {
            energy_degree_bounds: Some(vec![(0, 2)]),
            ..Default::default()
        };
        let components = denominator_connected_components(&parsed)
            .into_iter()
            .map(|component_edges| {
                let (component, embedding) =
                    project_denominator_component(&parsed, &component_edges)?;
                let component_options =
                    project_component_options(&options, parsed.internal_edges.len(), &embedding)?;
                generate_3d_expression_from_parsed_generated(&component, &component_options)
            })
            .collect::<Result<Vec<_>>>()
            .unwrap();

        assert_eq!(
            components
                .iter()
                .map(|generated| generated.energy_factor_ownership)
                .collect::<Vec<_>>(),
            vec![
                CffEnergyFactorOwnership::VariantLocal,
                CffEnergyFactorOwnership::GlobalSourceProduct,
            ]
        );
        assert_eq!(
            components
                .iter()
                .map(|generated| generated.core_global_prefactor_sign.factor())
                .collect::<Vec<_>>(),
            vec![1, -1],
            "branch-local generalized signs must not be inferred as a global bridge"
        );

        let generated = generate_3d_expression_from_parsed_generated(&parsed, &options).unwrap();
        assert_eq!(
            generated.energy_factor_ownership,
            CffEnergyFactorOwnership::VariantLocal
        );
        assert_eq!(generated.core_global_prefactor_sign.factor(), -1);
    }

    #[test]
    fn disconnected_component_product_rejects_shared_loop_variables() {
        let edge = |edge_id, tail, head, sign| ParsedGraphInternalEdge {
            edge_id,
            tail,
            head,
            label: format!("q{edge_id}"),
            mass_key: Some("m".to_string()),
            signature: MomentumSignature {
                loop_signature: vec![sign],
                external_signature: Vec::new(),
            },
            had_pow: false,
        };
        let parsed = ParsedGraph {
            internal_edges: vec![
                edge(0, 0, 1, 1),
                edge(1, 1, 0, -1),
                edge(2, 2, 3, 1),
                edge(3, 3, 2, -1),
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["k".to_string()],
            external_names: Vec::new(),
            node_name_to_internal: BTreeMap::from([
                ("a0".to_string(), 0),
                ("a1".to_string(), 1),
                ("b0".to_string(), 2),
                ("b1".to_string(), 3),
            ]),
        };

        assert!(matches!(
            generate_3d_expression_from_parsed_generated(
                &parsed,
                &Generate3DExpressionOptions::default(),
            ),
            Err(GenerationError::DisconnectedComponentsShareLoopVariables {
                left,
                right,
                variables,
            }) if left == vec![0, 1] && right == vec![2, 3] && variables == vec![0]
        ));
    }

    #[test]
    fn connected_multiloop_metadata_carries_pure_and_generalized_prefactor_signs() {
        let mut parsed = parsed_fixture("sunrise_pow4.dot");
        parsed.internal_edges.truncate(3);
        parsed.internal_edges[2].head = 1;

        let pure = generate_3d_expression_from_parsed_generated(
            &parsed,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(Vec::new()),
                ..Default::default()
            },
        )
        .unwrap();
        let generalized = generate_3d_expression_from_parsed_generated(
            &parsed,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(2, 2)]),
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(cff_duplicate_signature_excess(&parsed), 0);
        assert_eq!(pure.core_global_prefactor_sign.factor(), -1);
        assert_eq!(generalized.core_global_prefactor_sign.factor(), -1);
    }

    #[test]
    fn cff_generation_builds_repeated_box_with_branching_denominator_trees() {
        let parsed = parsed_fixture("box_pow3.dot");
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(Vec::new()),
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
    fn causal_repeated_generation_handles_energy_dependent_bubble_source() {
        let parsed = isolated_vertex_repeated_bubble_source();

        let generated = generate_3d_expression(
            &parsed,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(0, 1)]),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            generated.energy_factor_ownership,
            CffEnergyFactorOwnership::GlobalSourceProduct
        );
        let expression = generated.expression;

        assert!(!expression.orientations.is_empty());
        assert!(
            expression
                .orientations
                .iter()
                .flat_map(|orientation| &orientation.variants)
                .all(
                    |variant| variant.denominator_edges == vec![EdgeIndex(0), EdgeIndex(1)]
                        && !variant.half_edges.is_empty()
                )
        );
        assert!(
            expression
                .surfaces
                .linear_surface_cache
                .iter()
                .all(|surface| !surface.expression.is_zero())
        );
    }

    #[test]
    fn cff_generation_builds_multiloop_repeated_high_power_known_factor_completion() {
        let parsed = parsed_fixture("sunrise_pow4.dot");
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(2, 5)]),
                ..Default::default()
            },
        )
        .unwrap();

        assert!(!expression.orientations.is_empty());
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
        assert!(denominator_surface_ids.into_iter().all(|surface_id| {
            expression.surfaces.linear_surface_cache[LinearSurfaceID(surface_id)].kind
                == LinearSurfaceKind::Esurface
        }));
        let variants = expression
            .orientations
            .iter()
            .flat_map(|orientation| &orientation.variants)
            .collect::<Vec<_>>();
        assert!(variants.iter().any(|variant| {
            variant.origin.as_deref() == Some("bounded_degree_known_factor_cff")
                && !variant.numerator_surfaces.is_empty()
        }));
    }

    #[test]
    fn cff_generation_builds_one_loop_quadratic_e_surface_completion() {
        let parsed = parsed_fixture("box.dot");
        let ordinary = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(Vec::new()),
                ..Default::default()
            },
        )
        .unwrap();
        let bounded = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(0, 2)]),
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
                energy_degree_bounds: Some(vec![(0, 3)]),
                ..Default::default()
            },
        )
        .unwrap();

        assert!(!expression.orientations.is_empty());
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
        assert!(denominator_surface_ids.into_iter().all(|surface_id| {
            expression.surfaces.linear_surface_cache[LinearSurfaceID(surface_id)].kind
                == LinearSurfaceKind::Esurface
        }));
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
                energy_degree_bounds: Some(vec![(0, 1), (1, 1), (2, 1), (3, 3)]),
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
                energy_degree_bounds: Some(vec![(0, 4), (1, 2)]),
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
    fn cff_generation_builds_multiloop_high_contact_completion() {
        let parsed = penta_contact_graph();
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Cff,
                energy_degree_bounds: Some(vec![(0, 1), (2, 3), (4, 3)]),
                ..Default::default()
            },
        )
        .unwrap();

        assert!(!expression.orientations.is_empty());
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
        assert!(denominator_surface_ids.into_iter().all(|surface_id| {
            expression.surfaces.linear_surface_cache[LinearSurfaceID(surface_id)].kind
                == LinearSurfaceKind::Esurface
        }));
    }

    #[test]
    fn cff_generation_builds_one_loop_repeated_quadratic_recursive_completion() {
        let parsed = parsed_fixture("box_pow3.dot");
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Cff,
                energy_degree_bounds: Some(vec![(0, 2), (1, 1), (2, 1), (3, 2)]),
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
    fn cff_generation_samples_a_nonlinear_repeated_channel_as_one_emr_energy() {
        let parsed = parsed_fixture("box_pow3.dot");
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Cff,
                energy_degree_bounds: Some(vec![(3, 2), (5, 1)]),
                numerator_sampling_scale: NumeratorSamplingScaleMode::BeyondQuadratic,
                preserve_internal_edges_as_four_d_denominators: Vec::new(),
            },
        )
        .unwrap();

        assert!(!expression.orientations.is_empty());
        assert!(expression.orientations.iter().all(|orientation| {
            let channel_maps = &orientation.edge_energy_map[3..=5];
            channel_maps
                .iter()
                .all(|energy_map| energy_map == &channel_maps[0])
        }));
        assert!(expression.orientations.iter().any(|orientation| {
            orientation.edge_energy_map[3..=5]
                .iter()
                .all(LinearEnergyExpr::is_zero)
        }));
        assert!(expression.orientations.iter().any(|orientation| {
            orientation.edge_energy_map[3..=5]
                .iter()
                .all(LinearEnergyExpr::uses_uniform_scale)
        }));
    }

    #[test]
    fn cff_repeated_quadratic_channel_is_invariant_under_bound_ownership() {
        let mut parsed = parsed_fixture("box_pow3.dot");
        parsed.internal_edges[5]
            .signature
            .loop_signature
            .iter_mut()
            .for_each(|coefficient| *coefficient = -*coefficient);
        parsed.internal_edges[5]
            .signature
            .external_signature
            .iter_mut()
            .for_each(|coefficient| *coefficient = -*coefficient);

        let mut reference = None;
        for bounds in [
            vec![(3, 2)],
            vec![(4, 2)],
            vec![(5, 2)],
            vec![(3, 1), (4, 1)],
            vec![(3, 1), (5, 1)],
            vec![(4, 1), (5, 1)],
        ] {
            let generated = generate_3d_expression_from_parsed_generated(
                &parsed,
                &Generate3DExpressionOptions {
                    energy_degree_bounds: Some(bounds),
                    ..Default::default()
                },
            )
            .unwrap();
            assert_eq!(
                generated.energy_factor_ownership,
                CffEnergyFactorOwnership::VariantLocal
            );
            for orientation in &generated.expression.orientations {
                let representative = &orientation.edge_energy_map[3];
                let sample = if representative.is_zero() {
                    0
                } else if representative == &LinearEnergyExpr::ose(EdgeIndex(3), 1) {
                    1
                } else if representative == &LinearEnergyExpr::ose(EdgeIndex(3), -1) {
                    -1
                } else {
                    panic!(
                        "quadratic repeated-channel representative has a non-sample map: {representative:?}"
                    );
                };
                assert_eq!(
                    orientation.edge_energy_map[4],
                    LinearEnergyExpr::ose(EdgeIndex(4), sample),
                    "the same-routing occurrence must receive the representative sample"
                );
                assert_eq!(
                    orientation.edge_energy_map[5],
                    LinearEnergyExpr::ose(EdgeIndex(5), -sample),
                    "the reversed-routing occurrence must receive the signed representative sample"
                );
            }
            let encoded =
                bincode::encode_to_vec(&generated.expression, bincode::config::standard()).unwrap();
            if let Some(reference) = &reference {
                assert_eq!(
                    &encoded, reference,
                    "a common quadratic EMR channel must not depend on which Q/Q/-Q occurrence owns its bound"
                );
            } else {
                reference = Some(encoded);
            }
        }
    }

    #[test]
    fn cff_generation_builds_multiloop_single_quadratic_lower_sector_completion() {
        let parsed = parsed_fixture("sunrise_pow4.dot");
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(0, 2)]),
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
                energy_degree_bounds: Some(vec![(0, 1), (1, 1), (3, 4)]),
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
                energy_degree_bounds: Some(vec![(0, 1), (1, 1), (3, 4)]),
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
    fn cff_generation_keeps_distinct_maps_for_the_same_coarse_orientation() {
        let parsed = parsed_fixture("box_pow3.dot");
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Cff,
                energy_degree_bounds: Some(vec![(0, 1), (1, 1), (3, 4)]),
                numerator_sampling_scale: NumeratorSamplingScaleMode::All,
                preserve_internal_edges_as_four_d_denominators: Vec::new(),
            },
        )
        .unwrap();

        let (lhs, rhs) = expression
            .orientations
            .iter()
            .enumerate()
            .find_map(|(lhs_index, lhs)| {
                expression
                    .orientations
                    .iter()
                    .skip(lhs_index + 1)
                    .find(|rhs| {
                        lhs.data.orientation == rhs.data.orientation
                            && (lhs.loop_energy_map != rhs.loop_energy_map
                                || lhs.edge_energy_map != rhs.edge_energy_map)
                    })
                    .map(|rhs| (lhs, rhs))
            })
            .expect("raised CFF generation should retain exact maps sharing an orientation");

        assert_ne!(lhs.data.numerator_map_index, rhs.data.numerator_map_index);
        assert!(lhs.data.numerator_map_index.is_some());
        assert!(rhs.data.numerator_map_index.is_some());
        assert_eq!(
            lhs.data
                .label
                .as_deref()
                .and_then(|label| label.split_once('|'))
                .map(|(base, _)| base),
            rhs.data
                .label
                .as_deref()
                .and_then(|label| label.split_once('|'))
                .map(|(base, _)| base),
        );
    }

    #[test]
    fn cff_generation_builds_one_loop_high_power_with_repeated_spectators() {
        let parsed = parsed_fixture("box_pow3.dot");
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(0, 3), (1, 1), (2, 1), (3, 1)]),
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
        let variants = expression
            .orientations
            .iter()
            .flat_map(|orientation| &orientation.variants)
            .collect::<Vec<_>>();
        assert!(variants
            .iter()
            .any(|variant| variant.origin.as_deref() == Some("bounded_degree_known_factor_cff")));
        assert!(
            variants
                .iter()
                .any(|variant| !variant.numerator_surfaces.is_empty())
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
                energy_degree_bounds: Some(vec![(3, 2)]),
                numerator_sampling_scale: NumeratorSamplingScaleMode::BeyondQuadratic,
                preserve_internal_edges_as_four_d_denominators: Vec::new(),
            },
        )
        .unwrap();
        let all = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Cff,
                energy_degree_bounds: Some(vec![(3, 2)]),
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
