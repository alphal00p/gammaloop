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

/// Selects whether a CFF source is a standalone orientation catalogue or one
/// independently integrated CFF factor embedded in a larger product.
///
/// An embedded source keeps the ordinary causal sum whenever finite-pole
/// denominators remain. If its denominator set is itself a complete residue
/// basis, it retains the single deterministic Below residue instead of
/// reopening the two equivalent routings of every terminal contour.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize, Default,
)]
#[serde(rename_all = "snake_case")]
pub enum CffGenerationContext {
    #[default]
    Standalone,
    EmbeddedCffFactor,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct Generate3DExpressionOptions {
    #[serde(default)]
    pub representation: RepresentationMode,
    #[serde(default)]
    pub cff_generation_context: CffGenerationContext,
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
            cff_generation_context: CffGenerationContext::Standalone,
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

/// Ownership of the on-shell energy factors for one rational denominator
/// component. Edge IDs use the internal-energy namespace of the generated
/// expression. Component boundaries follow independently generated rational
/// source frames and survive their deterministic product. A vector-matroid
/// factorization internal to one Laurent functional remains in that single
/// frame because its numerator and variant-local factors are not independent.
#[derive(Debug, Clone, PartialEq, Eq, Encode, Decode)]
pub struct CffEnergyFactorComponent {
    pub internal_edge_ids: Vec<usize>,
    pub ownership: CffEnergyFactorOwnership,
    /// The scalar denominator source's uniform CFF convention in this
    /// independently generated rational component.
    pub denominator_only_global_prefactor_sign: CffGlobalPrefactorSign,
    /// The shared causal core's uniform contour convention in this
    /// independently generated rational component.
    pub core_global_prefactor_sign: CffGlobalPrefactorSign,
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
            denominator_only_global_prefactor_sign: self.denominator_only_global_prefactor_sign,
            core_global_prefactor_sign: self.core_global_prefactor_sign,
        }
    }
}

/// The uniform prefactor sign inserted by CFF generation.
///
/// For each source this contains the shared core's `(-1)^(L-1)` contour
/// convention. Pure CFF also includes its uniform duplicate-denominator sign;
/// when such a pure rational component is lifted into a generalized parent,
/// that component's uniform frame parity is retained here as well. Duplicate,
/// interpolation, closure, and contact signs introduced only inside individual
/// generalized variants remain variant-local and are deliberately not folded
/// into this record. The metadata accompanies the generated expression so
/// consumers can bridge their own prefactor convention without reconstructing
/// provenance from its final algebra.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Encode, Decode)]
pub struct CffGlobalPrefactorSign {
    odd: bool,
}

impl CffGlobalPrefactorSign {
    /// Construct the sign represented by `(-1)^exponent`.
    pub const fn from_exponent(exponent: usize) -> Self {
        Self {
            odd: !exponent.is_multiple_of(2),
        }
    }

    /// Compose two independent uniform prefactor conventions.
    pub const fn product(self, rhs: Self) -> Self {
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
    /// Rational-component metadata retained with persisted expressions.  The
    /// aggregate ownership above is insufficient for products mixing ordinary
    /// and generalized components, whose source-convention signs must be
    /// converted independently.
    pub energy_factor_components: Vec<CffEnergyFactorComponent>,
    /// Transient source-edge bounds supplied to this generation call. These always
    /// remain in the generator input's namespace; higher-level consumers that also
    /// need physical-parent bounds must retain them as separate metadata.
    pub source_energy_degree_bounds: Vec<(usize, usize)>,
    /// The scalar denominator source's uniform CFF convention, independent of
    /// numerator bounds and generalized interpolation. This is one ingredient
    /// of a consumer's componentwise convention bridge; it is not, by itself,
    /// the GammaLoop production bridge.
    pub denominator_only_global_prefactor_sign: CffGlobalPrefactorSign,
    pub core_global_prefactor_sign: CffGlobalPrefactorSign,
}

// Component boundaries and prefactor conventions are persisted because they
// define the source frame. Source bounds remain transient analysis metadata.
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
        bincode::Encode::encode(&self.energy_factor_components, encoder)?;
        bincode::Encode::encode(&self.denominator_only_global_prefactor_sign, encoder)?;
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
            energy_factor_components: bincode::Decode::decode(decoder)?,
            source_energy_degree_bounds: Vec::new(),
            denominator_only_global_prefactor_sign: bincode::Decode::decode(decoder)?,
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

    #[test]
    fn persisted_component_prefactor_metadata_roundtrips() {
        let generated: GeneratedThreeDExpression =
            GeneratedThreeDExpression {
                expression: ThreeDExpression::new_empty(),
                energy_factor_ownership: CffEnergyFactorOwnership::VariantLocal,
                energy_factor_components: vec![
                    CffEnergyFactorComponent {
                        internal_edge_ids: vec![3, 4],
                        ownership: CffEnergyFactorOwnership::VariantLocal,
                        denominator_only_global_prefactor_sign:
                            CffGlobalPrefactorSign::from_exponent(1),
                        core_global_prefactor_sign: CffGlobalPrefactorSign::from_exponent(0),
                    },
                    CffEnergyFactorComponent {
                        internal_edge_ids: vec![8],
                        ownership: CffEnergyFactorOwnership::GlobalSourceProduct,
                        denominator_only_global_prefactor_sign:
                            CffGlobalPrefactorSign::from_exponent(0),
                        core_global_prefactor_sign: CffGlobalPrefactorSign::from_exponent(1),
                    },
                ],
                source_energy_degree_bounds: vec![(3, 2)],
                denominator_only_global_prefactor_sign: CffGlobalPrefactorSign::from_exponent(1),
                core_global_prefactor_sign: CffGlobalPrefactorSign::from_exponent(0),
            };

        let bytes = bincode::encode_to_vec(&generated, bincode::config::standard()).unwrap();

        let (decoded, bytes_read): (GeneratedThreeDExpression, _) =
            bincode::decode_from_slice_with_context(
                &bytes,
                bincode::config::standard(),
                UnusedStateMap,
            )
            .unwrap();
        assert_eq!(bytes_read, bytes.len());
        assert_eq!(
            decoded.energy_factor_ownership,
            generated.energy_factor_ownership
        );
        assert_eq!(
            decoded.denominator_only_global_prefactor_sign,
            generated.denominator_only_global_prefactor_sign
        );
        assert_eq!(
            decoded.core_global_prefactor_sign,
            generated.core_global_prefactor_sign
        );
        assert_eq!(
            decoded.energy_factor_components,
            generated.energy_factor_components
        );
        assert!(decoded.source_energy_degree_bounds.is_empty());
        assert_eq!(
            bincode::encode_to_vec(decoded, bincode::config::standard()).unwrap(),
            bytes,
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
    if let Some(mut generated) = generate_rational_component_product(parsed, options)? {
        generated.expression = generated.expression.fuse_compatible_variants();
        assign_numerator_map_labels(&mut generated.expression.orientations);
        return Ok(generated);
    }

    let bounds = normalize_energy_degree_bounds(
        options.energy_degree_bounds.as_deref().unwrap_or(&[]),
        parsed.internal_edges.len(),
    )?;
    let uses_generalized_expression = cff_bounds_need_generalized_expression(parsed, &bounds);
    let denominator_edge_ids = parsed.denominator_internal_edge_ids();
    let duplicate_excess = cff_duplicate_signature_excess(parsed);
    // The scalar denominator convention is the product of the independent
    // rational contour frames. Incidence may join those factors at a vertex,
    // but it cannot replace their (-1)^(L-C) product by (-1)^(L-1).
    let denominator_only_global_prefactor_sign = denominator_contour_frame_sign(parsed);
    let embedded_terminal_basis = options.cff_generation_context
        == CffGenerationContext::EmbeddedCffFactor
        && denominator_set_is_complete_residue_basis(parsed);
    let (expression, energy_factor_ownership, core_global_prefactor_sign) =
        if embedded_terminal_basis {
            (
                KnownFactorCffBuilder::new(parsed, bounds, options.numerator_sampling_scale)
                    .build(true)?,
                if uses_generalized_expression {
                    CffEnergyFactorOwnership::VariantLocal
                } else {
                    CffEnergyFactorOwnership::GlobalSourceProduct
                },
                CffGlobalPrefactorSign::from_exponent(
                    parsed.loop_names.len().saturating_sub(1)
                        + if uses_generalized_expression {
                            0
                        } else {
                            duplicate_excess
                        },
                ),
            )
        } else if uses_generalized_expression {
            (
                // Keep the whole numerator in one Laurent functional. Rational
                // component decomposition remains available inside the lower-
                // sector builder, but is not a separate top-level generation
                // frame.
                BoundedCffBuilder::new(parsed, options)?.build()?,
                CffEnergyFactorOwnership::VariantLocal,
                CffGlobalPrefactorSign::from_exponent(parsed.loop_names.len().saturating_sub(1)),
            )
        } else {
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
    let energy_factor_components = (!denominator_edge_ids.is_empty())
        .then_some(CffEnergyFactorComponent {
            internal_edge_ids: denominator_edge_ids,
            ownership: energy_factor_ownership,
            denominator_only_global_prefactor_sign,
            core_global_prefactor_sign,
        })
        .into_iter()
        .collect();
    Ok(GeneratedThreeDExpression {
        expression,
        energy_factor_ownership,
        energy_factor_components,
        source_energy_degree_bounds: Vec::new(),
        denominator_only_global_prefactor_sign,
        core_global_prefactor_sign,
    })
}

/// Whether the denominator set closes every loop-energy contour without
/// leaving a finite-pole denominator behind.
///
/// Only this terminal case retains the deterministic Below residue when it is
/// embedded as one independently integrated factor of a larger product. A nonterminal or
/// repeated source still represents an ordinary generalized CFF sum.
fn denominator_set_is_complete_residue_basis(parsed: &ParsedGraph) -> bool {
    parsed.denominator_internal_edge_ids().len() == parsed.loop_names.len()
        && denominator_connected_components(parsed).len() == 1
        && repeated_groups(parsed).is_empty()
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
    if factorized.vector_matroid_components(&signatures).len() > 1 {
        // A graph may be vertex-connected while its rational denominator
        // factorizes into independent loop-energy components (for example a
        // tadpole attached at one vertex). Construct its causal product in
        // those independent coordinates instead of inventing surfaces which
        // mix components merely because they share an incidence vertex. The
        // component containing a structural cut retains that cut and its fixed
        // contour; every uncut component retains its complete public CFF sum.
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
        denominator_only_global_prefactor_sign: generated.denominator_only_global_prefactor_sign,
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
        denominator_only_global_prefactor_sign: CffGlobalPrefactorSign::default(),
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

fn generate_rational_component_product(
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
                        !parsed.is_initial_state_cut_edge(*edge_id)
                            && parsed.internal_edges[*edge_id].signature.loop_signature[*loop_id]
                                != 0
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
    let denominator_only_global_prefactor_sign = component_expressions
        .iter()
        .fold(CffGlobalPrefactorSign::default(), |sign, (generated, _)| {
            sign.product(generated.denominator_only_global_prefactor_sign)
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
        denominator_only_global_prefactor_sign,
        core_global_prefactor_sign,
    }))
}

fn denominator_connected_components(parsed: &ParsedGraph) -> Vec<Vec<usize>> {
    let denominator_edges = parsed.denominator_internal_edge_ids();
    if denominator_edges.is_empty() {
        return Vec::new();
    }
    let initial_state_cut_edges = parsed
        .initial_state_cut_edges
        .iter()
        .map(|cut_edge| cut_edge.edge_id)
        .collect::<BTreeSet<_>>();

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
    // A cut whose endpoints already lie in one denominator component is an
    // internal structural alias of that component and must survive projection.
    // A cut between independent denominator components remains an external
    // boundary flow for both projections; it must not merge their contour
    // subspaces. Cut-only components have no denominator contour and are not
    // emitted.
    for edge_id in initial_state_cut_edges {
        let edge = &parsed.internal_edges[edge_id];
        if parent.contains_key(&edge.tail) && parent.contains_key(&edge.head) {
            let tail_root = find_node_root(&mut parent, edge.tail);
            let head_root = find_node_root(&mut parent, edge.head);
            if tail_root == head_root {
                grouped.entry(tail_root).or_default().push(edge_id);
            }
        }
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
                    !parsed.is_initial_state_cut_edge(*edge_id)
                        && parsed.internal_edges[*edge_id].signature.loop_signature[loop_id] != 0
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
        cff_generation_context: options.cff_generation_context,
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

    for partial in &mut partials {
        // Cuts crossing components (and cut-only aliases) are deliberately not
        // owned by a denominator factor. Restore their fixed external-energy
        // maps after the product lift so they cannot remain at the identity's
        // zero placeholder.
        apply_initial_state_cut_edge_energy_exprs(parsed, &mut partial.edge_energy_map);
        let (orientation, label) = orientation_from_edge_exprs(&partial.edge_energy_map);
        partial.data.orientation = orientation;
        partial.data.label = Some(label);
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
    contour_closure: &[ContourClosure],
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
    let denominator_rows = loop_line_signatures
        .iter()
        .map(|row| {
            row.iter()
                .map(|value| i64::from(*value))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    // This constructor is the terminal D=L base of generalized CFF recursion,
    // where every denominator belongs to the residue basis and no finite-pole
    // factors remain. Keep the broader direct-residue expression out of the
    // production path.
    if denominator_edge_ids.len() != n_loops || rank_i64(&denominator_rows) != n_loops {
        return Err(GenerationError::CffHigherEnergyPowerNotImplemented);
    }
    let residues = energy_residues(&loop_line_signatures, contour_closure)?;

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
        let mut dividend = self.clone();
        let mut remainder = Self::zero();
        while let Some((lead_monomial, lead_coeff)) = dividend.leading_term() {
            if let Some(quotient_monomial) =
                lead_monomial.quotient_after_dividing_by(&divisor_monomial)
            {
                let quotient_coeff = lead_coeff / divisor_coeff.clone();
                let quotient_term =
                    Self::from_monomial(quotient_monomial.clone(), quotient_coeff.clone());
                quotient.add_term(quotient_monomial, quotient_coeff);
                dividend = dividend - quotient_term * divisor.clone();
            } else {
                let lead_term = Self::from_monomial(lead_monomial.clone(), lead_coeff.clone());
                remainder.add_term(lead_monomial, lead_coeff);
                dividend = dividend - lead_term;
            }
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
        let edges = self.parsed.denominator_internal_edge_ids();
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
    inherited_contour_rows: Vec<Vec<i32>>,
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
            inherited_contour_rows: Vec::new(),
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
                .build_with_inherited_contours(false, &self.inherited_contour_rows);
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
                .build_with_inherited_contours(false, &self.inherited_contour_rows);
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
                // The finite-pole remainder is the original denominator
                // sector. Its duplicate-line parity is therefore part of the
                // Hermite functional itself, rather than a child-to-parent
                // core convention to consume at the embedding boundary.
                generate_pure_cff_expression_from_parsed(self.parsed)
            };
        };

        let mut remainder_bounds = self.bounds.clone();
        remainder_bounds[active_edge] = 1;
        let mut remainder_builder = BoundedCffBuilder::for_bounds(self.parsed, remainder_bounds);
        remainder_builder.inherited_contour_rows = self.inherited_contour_rows.clone();
        remainder_builder.sampling_scale_mode = self.sampling_scale_mode;
        let remainder = remainder_builder.build_quadratic_recursive(lower_sector_base)?;
        self.append_recursive_remainder_terms(active_edge, &remainder)?;

        let (subparsed, sub_to_orig) = self.project_parsed_edges(&[active_edge]);
        // Projection consumes only the active edge's contact. Every surviving
        // pre-existing source bound remains immutable in the lower sector.
        let sub_bounds = sub_to_orig
            .iter()
            .map(|orig_id| self.bounds[*orig_id])
            .collect::<Vec<_>>();
        let mut contact_builder = BoundedCffBuilder::for_bounds(&subparsed, sub_bounds);
        contact_builder.inherited_contour_rows = self.inherited_contour_rows.clone();
        contact_builder.inherited_contour_rows.push(
            self.parsed.internal_edges[active_edge]
                .signature
                .loop_signature
                .clone(),
        );
        contact_builder.sampling_scale_mode = self.sampling_scale_mode;
        let mut contact = contact_builder.build_quadratic_recursive(true)?;
        let active_is_repeated_denominator = repeated_groups(self.parsed)
            .iter()
            .any(|group| group.edge_ids.contains(&active_edge));
        if active_is_repeated_denominator
            && subparsed.denominator_internal_edge_ids().len() == subparsed.loop_names.len()
            && denominator_connected_components(&subparsed).len() == 1
            && repeated_groups(&subparsed).is_empty()
        {
            // Removing one occurrence of an exact repeated denominator leaves
            // its terminal pole in the duplicate-channel frame. Consume that
            // frame transition; unequal-mass parallel lines are not repeated.
            for variant in contact
                .orientations
                .iter_mut()
                .flat_map(|orientation| &mut orientation.variants)
            {
                variant.prefactor =
                    rational_to_coefficient(-rational_from_coefficient(&variant.prefactor))?;
            }
        }
        self.append_recursive_contact_terms(
            active_edge,
            &contact,
            &sub_to_orig,
            &contact_weight_polys(self.bounds[active_edge]),
        )?;

        self.finalize_numerator_map_labels();
        Ok(self.expression)
    }

    fn lower_sector_base_expression(&self) -> Result<ThreeDExpression<OrientationID>> {
        // Once a recursive contact leaves exactly one denominator basis, its
        // pole inherits the parent's ordered Below contour. Retain that signed
        // residue instead of reopening the terminal with both free directions.
        // A powered or nonterminal lower sector keeps the standard component
        // construction; exact rational normalization remains encoded directly
        // in the generated expression coefficients.
        if self.parsed.denominator_internal_edge_ids().len() == self.parsed.loop_names.len()
            && denominator_connected_components(self.parsed).len() == 1
            && repeated_groups(self.parsed).is_empty()
        {
            generate_simple_residue_basis_expression_from_parsed(
                self.parsed,
                &vec![ContourClosure::Below; self.parsed.loop_names.len()],
            )
        } else {
            let mut lower = LowerSectorCffBuilder::new(self.parsed);
            lower.inherited_contour_rows = self.inherited_contour_rows.clone();
            lower.build()
        }
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
        // Lower-sector variants already carry their complete scalar
        // normalization in their coefficients.
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
            generate_simple_residue_basis_expression_from_parsed(
                &subparsed,
                &vec![ContourClosure::Below; subparsed.loop_names.len()],
            )?
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
                        && !pinched_edges.len().is_multiple_of(2)
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

    fn build(self, lower_sector_base: bool) -> Result<ThreeDExpression<OrientationID>> {
        self.build_with_inherited_contours(lower_sector_base, &[])
    }

    fn build_with_inherited_contours(
        mut self,
        lower_sector_base: bool,
        inherited_contour_rows: &[Vec<i32>],
    ) -> Result<ThreeDExpression<OrientationID>> {
        let local_to_orig = (0..self.original.internal_edges.len()).collect::<Vec<_>>();
        let recursion_budget = self.recursion_budget();
        self.channel_recursive_terms(
            self.original,
            &local_to_orig,
            &BTreeMap::new(),
            &[],
            &[],
            inherited_contour_rows,
            0,
            Rational::one(),
            lower_sector_base,
            false,
            0,
            recursion_budget,
        )?;
        // Pure CFF stores the shared connected-core contour sign directly in
        // every variant. Known-factor branches therefore bridge only a raw
        // inherited lower contour at their append boundary. A post-build flip
        // would apply the same sign twice to ordinary finite-pole branches.
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
        inherited_contour_rows: &[Vec<i32>],
        extra_uniform_scale_power: usize,
        prefactor: Rational,
        lower_sector_base: bool,
        repeated_channel_normal_form_consumed: bool,
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
                inherited_contour_rows,
                extra_uniform_scale_power,
                prefactor,
                lower_sector_base,
                repeated_channel_normal_form_consumed,
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
                // Removing serial copies of one repeated denominator channel
                // selects an algebraic quotient; it does not integrate out a
                // loop-energy direction. Preserve only genuine inherited
                // contour closures from the parent recursion.
                let next_contour_rows = inherited_contour_rows.to_vec();

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
                    next_known.push(KnownLinearExpr::var(rep_local, 1));
                }
                if term.positive_ose_power != 0 {
                    next_known.extend(
                        (0..term.positive_ose_power).map(|_| KnownLinearExpr::ose(rep_orig, 1)),
                    );
                }
                if term.cancelled_power != 0 {
                    let y_expr = KnownLinearExpr::var(rep_local, 1);
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
                    .map(|factor| {
                        self.remap_known_factor_to_sub(parsed, &subparsed, &sub_to_local, factor)
                    })
                    .collect::<Result<Vec<_>>>()
                {
                    Ok(sub_known) => self.channel_recursive_terms(
                        &subparsed,
                        &sub_local_to_orig,
                        &sub_replacements,
                        &sub_known,
                        &next_half_edges,
                        &next_contour_rows,
                        extra_uniform_scale_power + channel_uniform_power,
                        channel_prefactor,
                        lower_sector_base || !delete.is_empty(),
                        true,
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
                            &next_contour_rows,
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
        inherited_contour_rows: &[Vec<i32>],
        extra_uniform_scale_power: usize,
        prefactor: Rational,
        lower_sector_base: bool,
        repeated_channel_normal_form_consumed: bool,
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
        // factorized and are completed by the lower-sector polynomial normal
        // form below, rather than resampling the black-box numerator.
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
            // Once interpolation has resolved every black-box numerator
            // direction, keep derivative-generated factors in the exact base
            // CFF below. Their cut-aware polynomial normal form is only the
            // explicit singular/remap-error fallback; eagerly resampling them
            // here would make the answer depend on unused interpolation
            // capacity.
            return self.append_base_terms(
                parsed,
                local_to_orig,
                replacements,
                &known_factors,
                extra_half_edges,
                inherited_contour_rows,
                extra_uniform_scale_power,
                prefactor.clone(),
                lower_sector_base,
                repeated_channel_normal_form_consumed,
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
                inherited_contour_rows,
                extra_uniform_scale_power,
                prefactor.clone(),
                lower_sector_base,
                repeated_channel_normal_form_consumed,
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
        let mut next_contour_rows = inherited_contour_rows.to_vec();
        next_contour_rows.push(active_signature.clone());
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
                    .map(|factor| {
                        self.remap_known_factor_to_sub(parsed, &subparsed, &sub_to_local, factor)
                    })
                    .collect::<Result<Vec<_>>>()
                {
                    Ok(remapped_factors) => self.recursive_terms(
                        &subparsed,
                        &sub_local_to_orig,
                        &next_replacements,
                        &remapped_factors,
                        &next_half_edges,
                        &next_contour_rows,
                        extra_uniform_scale_power,
                        branch_prefactor,
                        true,
                        repeated_channel_normal_form_consumed,
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
                            &next_contour_rows,
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
        inherited_contour_rows: &[Vec<i32>],
        extra_uniform_scale_power: usize,
        prefactor: Rational,
    ) -> Result<()> {
        let numerator = KnownPolynomial::product_from_known_factors(parent, known_factors);
        // This path is reached only when a derivative-generated factor cannot
        // be remapped to the surviving graph. Keep its exact polynomial normal
        // form as the singular-projection fallback; ordinary lower contacts
        // stay factorized and obtain bounded CFF maps in `append_base_terms`.
        let branches = KnownPolynomialNormalForm::new(parsed, local_to_orig).branches(numerator)?;
        for branch in branches {
            self.append_known_polynomial_branch(
                parsed,
                local_to_orig,
                replacements,
                branch,
                extra_half_edges,
                inherited_contour_rows,
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
        inherited_contour_rows: &[Vec<i32>],
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
        )? {
            return Err(GenerationError::CffHigherEnergyPowerNotImplemented);
        }

        let denominator_set = branch
            .denominator_edges
            .iter()
            .copied()
            .collect::<BTreeSet<_>>();
        let delete = parsed
            .denominator_internal_edge_ids()
            .into_iter()
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
        let mut branch_contour_rows = inherited_contour_rows.to_vec();
        branch_contour_rows.extend(delete.iter().map(|local_id| {
            parsed.internal_edges[*local_id]
                .signature
                .loop_signature
                .clone()
        }));
        let original_signatures = self
            .original
            .internal_edges
            .iter()
            .map(|edge| edge.signature.clone())
            .collect::<Vec<_>>();
        let mut reconstruction_edges = branch_local_to_orig.clone();
        reconstruction_edges.extend(
            replacements
                .keys()
                .filter(|edge_id| !branch_local_to_orig.contains(edge_id))
                .copied(),
        );
        let reconstruction_basis = LowerSectorCffBuilder::component_basis_edges(
            &original_signatures,
            &reconstruction_edges,
        );
        // This is an explicitly selected rational quotient of the generalized
        // CFF, not a simplification of the physical 4D numerator. Construct it
        // through the same lower-sector CFF path as every other quotient; no
        // independent contour or residue-frame oracle is introduced.
        let mut lower_sector = LowerSectorCffBuilder::new(&branch_parsed);
        lower_sector.inherited_contour_rows = branch_contour_rows;
        let base_expression = lower_sector.build()?;
        let branch_prefactor = prefactor;
        let edge_map = branch_local_to_orig
            .iter()
            .enumerate()
            .map(|(local_id, orig_id)| (local_id, *orig_id))
            .collect::<BTreeMap<_, _>>();
        let surface_map = self.copy_expression_surfaces(&base_expression, &edge_map);
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
            for orig_id in reconstruction_basis
                .iter()
                .filter(|edge_id| !branch_local_to_orig.contains(edge_id))
            {
                target_edge_exprs[*orig_id] = replacements[orig_id].clone();
            }
            let loop_exprs = solve_loop_energy_particular_from_target_edge_exprs(
                &original_signatures,
                &reconstruction_basis,
                &target_edge_exprs,
            )?;
            let mut full_edge_exprs = edge_q0_from_loop_exprs(&original_signatures, &loop_exprs);
            for (local_id, orig_id) in branch_local_to_orig.iter().copied().enumerate() {
                // Loop reconstruction completes only directions deleted by the
                // quotient. Every surviving physical occurrence keeps the
                // public CFF map generated for that occurrence, including a
                // degenerate sibling with its own OSE identity.
                full_edge_exprs[orig_id] = orientation.edge_energy_map[local_id]
                    .clone()
                    .remap_internal_edges(&edge_map);
            }
            for (orig_id, expr) in replacements {
                full_edge_exprs[*orig_id] = expr.clone();
            }

            for variant in &orientation.variants {
                let base_coeff =
                    rational_from_coefficient(&variant.prefactor) * branch_prefactor.clone();
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
    ) -> Result<bool> {
        let loop_variables = numerator.loop_variables();
        if loop_variables.is_empty() {
            return Ok(true);
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
        Ok(loop_variables.into_iter().all(|loop_id| {
            let mut unit = vec![0; self.original.loop_names.len()];
            unit[loop_id] = 1;
            lower_sector
                .row_coordinates_in_basis(&basis_rows, &unit)
                .is_ok()
        }))
    }

    #[allow(clippy::too_many_arguments)]
    fn append_base_terms(
        &mut self,
        parsed: &ParsedGraph,
        local_to_orig: &[usize],
        replacements: &BTreeMap<usize, LinearEnergyExpr>,
        known_factors: &[KnownLinearExpr],
        extra_half_edges: &[usize],
        inherited_contour_rows: &[Vec<i32>],
        extra_uniform_scale_power: usize,
        prefactor: Rational,
        lower_sector_base: bool,
        repeated_channel_normal_form_consumed: bool,
    ) -> Result<()> {
        if KnownLinearExpr::product_is_zero(known_factors) {
            return Ok(());
        }
        let total_bounds =
            self.known_total_bounds(parsed, local_to_orig, replacements, known_factors);
        let bounded_lower = if repeated_channel_normal_form_consumed
            && cff_bounds_need_generalized_expression(parsed, &total_bounds)
        {
            let mut bounded = BoundedCffBuilder::for_bounds(parsed, total_bounds.clone());
            bounded.inherited_contour_rows = inherited_contour_rows.to_vec();
            bounded.sampling_scale_mode = self.sampling_scale_mode;
            Some(bounded.build()?)
        } else {
            None
        };
        let lower_sector_base_expression = || {
            let mut lower = LowerSectorCffBuilder::new(parsed);
            lower.inherited_contour_rows = inherited_contour_rows.to_vec();
            lower.build()
        };
        let direct_base = || -> Result<ThreeDExpression<OrientationID>> {
            // Interpolation has already separated every powered channel before
            // reaching this base.  What remains is the ordinary rational CFF
            // of this exact denominator product, not another lower-sector
            // component reconstruction. Keep its duplicate-line parity on the
            // direct variants in the inherited parent functional.
            let duplicate_excess = cff_duplicate_signature_excess(parsed);
            generate_pure_cff_expression_from_parsed_with_duplicate_excess(parsed, duplicate_excess)
        };
        let uses_terminal_residue_basis = lower_sector_base
            && parsed.denominator_internal_edge_ids().len() == parsed.loop_names.len()
            && denominator_connected_components(parsed).len() == 1
            && repeated_groups(parsed).is_empty();
        let base = if let Some(base) = bounded_lower {
            base
        } else if lower_sector_base
            && parsed.initial_state_cut_edges.is_empty()
            && !repeated_groups(parsed).is_empty()
        {
            // Reuse the component constructor. A repeated channel restores
            // its surviving duplicate parity directly on the generated
            // variants; the standard lower-sector constructor owns the
            // remaining component relation.
            lower_sector_base_expression()?
        } else if uses_terminal_residue_basis {
            // A fully reduced denominator basis inherits the ordered Below
            // contour which produced this contact. Retain that single residue
            // instead of reopening both standalone CFF directions. A root
            // factorized pole and the reduction of an exact repeated channel
            // consume the one-minus bridge to their public/duplicate core;
            // an ordinary unequal-mass contact already carries the signed
            // inherited residue in the parent frame.
            let root_factorized_pole = inherited_contour_rows.is_empty();
            let reduced_repeated_channel = repeated_groups(self.original).iter().any(|group| {
                let surviving = group
                    .edge_ids
                    .iter()
                    .filter(|edge_id| local_to_orig.contains(edge_id))
                    .count();
                surviving != 0 && surviving != group.edge_ids.len()
            });
            let needs_terminal_bridge = root_factorized_pole || reduced_repeated_channel;
            let mut terminal = if root_factorized_pole {
                generate_simple_residue_basis_expression_from_parsed(
                    parsed,
                    &vec![ContourClosure::Below; parsed.loop_names.len()],
                )?
            } else {
                lower_sector_base_expression()?
            };
            if needs_terminal_bridge {
                for variant in terminal
                    .orientations
                    .iter_mut()
                    .flat_map(|orientation| &mut orientation.variants)
                {
                    variant.prefactor =
                        rational_to_coefficient(-rational_from_coefficient(&variant.prefactor))?;
                }
            }
            terminal
        } else if lower_sector_base {
            lower_sector_base_expression()?
        } else if !known_factors.is_empty() || !replacements.is_empty() {
            // The surrounding generalized branch has selected this exact
            // denominator product. Once interpolation has consumed every
            // powered channel it is an ordinary direct CFF; otherwise retain
            // the established derivative-aware lower-sector construction.
            if repeated_channel_normal_form_consumed {
                direct_base()?
            } else {
                lower_sector_base_expression()?
            }
        } else {
            generate_pure_cff_expression_from_parsed(parsed)?
        };
        let base_expression = base;
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
        sub_to_parent: &[usize],
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
            let expr = if let Some(sub_edge_id) = sub_to_parent
                .iter()
                .position(|parent_edge_id| parent_edge_id == edge_id)
            {
                // Projection clones a surviving occurrence verbatim. Keep its
                // own generalized sample rather than rebuilding the same
                // physical source loop-coordinate momentum through unrelated
                // sampled edges.
                KnownLinearExpr::var(sub_edge_id, 1)
            } else {
                self.known_signature_expr(subparsed, &parsed.internal_edges[*edge_id].signature)?
            };
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
    requires_inherited_component_bridge: bool,
}

#[derive(Debug, Clone)]
struct LowerSectorPartial {
    coeff: Rational,
    half_edges: Vec<usize>,
    denominator_edges: Vec<usize>,
    chain: Vec<HybridSurfaceID>,
    numerator_surfaces: Vec<HybridSurfaceID>,
    denominator_surface_signs: BTreeMap<HybridSurfaceID, i64>,
    denominator_edge_support_signs: BTreeMap<Vec<EdgeIndex>, i64>,
    uniform_scale_power: usize,
    origins: Vec<String>,
    targets: BTreeMap<usize, LinearEnergyExpr>,
    edge_exprs: BTreeMap<usize, LinearEnergyExpr>,
}

struct LowerSectorCffBuilder<'a> {
    parsed: &'a ParsedGraph,
    force_component_factorization: bool,
    inherited_contour_rows: Vec<Vec<i32>>,
    expression: ThreeDExpression<OrientationID>,
    surface_index: HashMap<(LinearSurfaceKind, LinearEnergyExpr), HybridSurfaceID>,
}

impl<'a> LowerSectorCffBuilder<'a> {
    fn new(parsed: &'a ParsedGraph) -> Self {
        Self {
            parsed,
            force_component_factorization: false,
            inherited_contour_rows: Vec::new(),
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
            && self.inherited_contour_rows.is_empty()
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
        let needs_inherited_component_bridge = components
            .iter()
            .any(|component| component.requires_inherited_component_bridge);
        // `assemble_component_product` multiplies one independently closed
        // rational contour per vector-matroid component. A top-level product
        // needs the relative (-1)^(C-1) convention. A nested lower sector needs
        // that bridge exactly when its outermost deleted contour still closes
        // on a surviving denominator carrier. If that first contour has no
        // carrier, its parent contact already fixed the residual product frame;
        // later inherited closures act inside that frame and cannot reopen it.
        // Closure direction and its Jacobian remain component-local signs.
        let starts_component_product_frame = self.inherited_contour_rows.is_empty()
            || self.force_component_factorization
            || needs_inherited_component_bridge;
        let component_product_exponent =
            components.len().saturating_sub(1) * usize::from(starts_component_product_frame);
        let component_product_sign = if self.parsed.initial_state_cut_edges.is_empty()
            && !component_product_exponent.is_multiple_of(2)
        {
            Rational::from(-1)
        } else {
            Rational::one()
        };
        self.assemble_component_product(&components, component_product_sign)
    }

    fn assemble_component_product(
        mut self,
        components: &[LowerSectorComponent],
        initial_coeff: Rational,
    ) -> Result<ThreeDExpression<OrientationID>> {
        let signatures = self.signatures();
        let mut partials = vec![LowerSectorPartial {
            coeff: initial_coeff,
            half_edges: Vec::new(),
            denominator_edges: Vec::new(),
            chain: Vec::new(),
            numerator_surfaces: Vec::new(),
            denominator_surface_signs: BTreeMap::new(),
            denominator_edge_support_signs: BTreeMap::new(),
            uniform_scale_power: 0,
            origins: Vec::new(),
            targets: BTreeMap::new(),
            edge_exprs: BTreeMap::new(),
        }];

        for component in components {
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
                        item.coeff = item.coeff * rational_from_coefficient(&variant.prefactor);
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
                        for (surface, sign) in &variant.denominator_surface_signs {
                            let surface = map_surface_id(*surface, &surface_map);
                            *item.denominator_surface_signs.entry(surface).or_insert(1) *= sign;
                        }
                        for (support, sign) in map_edge_support_signs(
                            &variant.denominator_edge_support_signs,
                            &edge_map,
                        ) {
                            *item
                                .denominator_edge_support_signs
                                .entry(support)
                                .or_insert(1) *= sign;
                        }
                        item.uniform_scale_power += variant.uniform_scale_power;
                        item.origins.push(
                            variant
                                .origin
                                .clone()
                                .unwrap_or_else(|| "anonymous".to_string()),
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
            let origin = format!(
                "lower_sector_cff_e_surface_component_product:{}",
                partial.origins.join(":")
            );
            self.push_variant_for_maps(
                loop_exprs,
                edge_exprs,
                crate::expression::CFFVariant {
                    origin: Some(origin),
                    prefactor: rational_to_coefficient(partial.coeff)?,
                    half_edges: half_edges.into_iter().map(EdgeIndex).collect(),
                    denominator_edges: partial
                        .denominator_edges
                        .into_iter()
                        .map(EdgeIndex)
                        .collect(),
                    denominator_surface_signs: partial.denominator_surface_signs,
                    denominator_edge_support_signs: partial.denominator_edge_support_signs,
                    uniform_scale_power: partial.uniform_scale_power,
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
        let components = self.vector_matroid_components(signatures);
        let component_bases = components
            .iter()
            .map(|component| Self::component_basis_edges(signatures, component))
            .collect::<Vec<_>>();
        // Keep the canonical component order while cumulatively resolving an
        // inherited contour. Which component supplies the final independent
        // pivot is part of the residue construction: reordering a terminal
        // D=L component behind a nonterminal D>L component can replace the
        // latter's complete public CFF by one compact pole representative and
        // thereby lose half of a contact term.
        let component_elimination_order = 0..components.len();
        let mut inherited_closures = vec![Vec::new(); components.len()];
        let mut inherited_component_bridges = vec![false; components.len()];
        for (inherited_row_index, inherited_row) in self.inherited_contour_rows.iter().enumerate() {
            let mut ordered_basis_rows = Vec::new();
            let mut selected = None;
            for component_id in component_elimination_order.clone() {
                let basis_edges = &component_bases[component_id];
                for (local_basis_id, edge_id) in basis_edges.iter().copied().enumerate() {
                    let mut basis_row = signatures[edge_id].loop_signature.clone();
                    if basis_row
                        .iter()
                        .find(|coefficient| **coefficient != 0)
                        .is_some_and(|coefficient| *coefficient < 0)
                    {
                        basis_row
                            .iter_mut()
                            .for_each(|coefficient| *coefficient *= -1);
                    }
                    ordered_basis_rows.push(basis_row);
                    let Ok(coordinates) =
                        self.rational_row_coordinates_in_basis(&ordered_basis_rows, inherited_row)
                    else {
                        continue;
                    };
                    let Some((numerator, _)) = coordinates.last().and_then(Rational::to_i64_pair)
                    else {
                        return Err(GenerationError::CoefficientOutOfRange);
                    };
                    if numerator == 0 {
                        continue;
                    }
                    selected = Some((
                        component_id,
                        local_basis_id,
                        if numerator > 0 {
                            ContourClosure::Below
                        } else {
                            ContourClosure::Above
                        },
                        numerator.signum(),
                    ));
                    break;
                }
                if selected.is_some() {
                    break;
                }
            }
            // If the deleted contour coordinate has no surviving denominator
            // carrier, its contact is already fully localized and there is no
            // lower-sector component whose public orientation sum must be
            // restricted.
            let Some((component_id, local_basis_id, closure, jacobian_sign)) = selected else {
                continue;
            };
            // Contact rows are appended outermost to innermost. Only the
            // outermost row decides whether the parent's component-product
            // contour frame survives into this lower sector. Its routing sign
            // is deliberately not folded into this predicate: `closure` and
            // `jacobian_sign` already carry that oriented residue data.
            if inherited_row_index == 0 {
                inherited_component_bridges[component_id] = true;
            }
            inherited_closures[component_id].push((local_basis_id, closure, jacobian_sign));
        }
        components
            .into_iter()
            .zip(inherited_closures)
            .zip(inherited_component_bridges)
            .map(|((edges, closures), bridges)| {
                self.component_bundle(signatures, edges, closures, bridges)
            })
            .collect()
    }

    fn component_bundle(
        &self,
        signatures: &[MomentumSignature],
        edges: Vec<usize>,
        inherited_closures: Vec<(usize, ContourClosure, i64)>,
        inherited_component_bridge: bool,
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
        // Rationally independent components retain their complete explicit CFF
        // orientation sums even when incidence joins them at a vertex. Only a
        // structural cut internal to this component fixes its contour and may
        // select the compact residue basis. Otherwise compose the existing
        // public CFF, including its own duplicate-denominator parity.
        let denominator_count = component_parsed.denominator_internal_edge_ids().len();
        let mut contour_closure = vec![ContourClosure::Below; rank];
        let mut assigned_closures = vec![None; rank];
        for (loop_id, closure, jacobian_sign) in &inherited_closures {
            // Nested contacts are ordered from outermost to innermost. If
            // contractions make two deleted rows close on the same surviving
            // pivot, the innermost contact owns that terminal residue; the
            // outer contour was already consumed by its parent coefficient.
            assigned_closures[*loop_id] = Some((*closure, *jacobian_sign));
        }
        let mut inherited_jacobian_sign = 1i64;
        for (loop_id, assigned) in assigned_closures.into_iter().enumerate() {
            if let Some((closure, jacobian_sign)) = assigned {
                contour_closure[loop_id] = closure;
                inherited_jacobian_sign *= jacobian_sign;
            }
        }
        let uses_simple_residue_basis = denominator_count == rank
            && (!component_parsed.initial_state_cut_edges.is_empty()
                || !inherited_closures.is_empty());
        let mut expression = if uses_simple_residue_basis {
            generate_simple_residue_basis_expression_from_parsed(
                &component_parsed,
                &contour_closure,
            )?
        } else {
            generate_pure_cff_expression_from_parsed(&component_parsed)?
        };
        if uses_simple_residue_basis && inherited_jacobian_sign < 0 {
            for variant in expression
                .orientations
                .iter_mut()
                .flat_map(|orientation| &mut orientation.variants)
            {
                variant.prefactor =
                    rational_to_coefficient(-rational_from_coefficient(&variant.prefactor))?;
            }
        }
        Ok(LowerSectorComponent {
            basis_edges,
            local_to_sub,
            expression,
            requires_inherited_component_bridge: inherited_component_bridge,
        })
    }

    fn project_component_parsed(
        &self,
        signatures: &[MomentumSignature],
        edges: &[usize],
        basis_edges: &[usize],
    ) -> Result<(ParsedGraph, Vec<usize>)> {
        // Orient component coordinates in the source loop-coordinate basis,
        // independently of which physical edge happened to supply the basis
        // row. An oppositely routed edge then keeps a local -1 row and hence
        // its inherited pole.
        let basis_rows = basis_edges
            .iter()
            .map(|edge_id| {
                let mut row = signatures[*edge_id].loop_signature.clone();
                if row
                    .iter()
                    .find(|coefficient| **coefficient != 0)
                    .is_some_and(|coefficient| *coefficient < 0)
                {
                    row.iter_mut().for_each(|coefficient| *coefficient *= -1);
                }
                row
            })
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
        // A fixed cut belongs structurally to this rational component only if
        // both of its endpoints already lie in the component selected by the
        // denominator rows. A cut touching or crossing components is boundary
        // flow; it must neither merge their loop spaces nor change the parity
        // of the factor it merely touches.
        let retained_cut_edges = self
            .parsed
            .initial_state_cut_edges
            .iter()
            .filter_map(|cut_edge| {
                let edge = &self.parsed.internal_edges[cut_edge.edge_id];
                (component_nodes.contains(&edge.tail) && component_nodes.contains(&edge.head))
                    .then_some(cut_edge.edge_id)
            })
            .collect::<BTreeSet<_>>();
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
        for sub_id in retained_cut_edges.iter().copied() {
            let original = &self.parsed.internal_edges[sub_id];
            let new_id = internal_edges.len();
            local_to_sub.push(sub_id);
            internal_edges.push(ParsedGraphInternalEdge {
                edge_id: new_id,
                tail: old_to_new[&original.tail],
                head: old_to_new[&original.head],
                label: original.label.clone(),
                mass_key: original.mass_key.clone(),
                signature: MomentumSignature {
                    // The cut energy is the fixed external alias below. Its
                    // stored source row is provenance and cannot contribute a
                    // contour variable or component rank.
                    loop_signature: vec![0; rank],
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
        // Preserve genuine source-graph external insertions which touch this
        // component. A denominator in another vector-matroid component is a
        // multiplicative rational factor, not external boundary data here.
        // Likewise, a structural cut which was not retained above belongs to
        // another rational factor; merely touching this component at one
        // endpoint cannot anchor its otherwise independent contour.
        let external_edges = self
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

fn denominator_contour_frame_exponent(parsed: &ParsedGraph) -> usize {
    let signatures = parsed
        .internal_edges
        .iter()
        .map(|edge| edge.signature.clone())
        .collect::<Vec<_>>();
    let denominator_rows = parsed
        .denominator_internal_edge_ids()
        .into_iter()
        .map(|edge_id| {
            signatures[edge_id]
                .loop_signature
                .iter()
                .map(|coefficient| i64::from(*coefficient))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    // Projected lower sectors retain the parent's loop-coordinate namespace,
    // including directions already consumed by deleted contours. Their
    // denominator frame is set by the rank that remains in the rational
    // denominator, not by the number of names in that ambient namespace.
    let denominator_rank = rank_i64(&denominator_rows);
    let rational_component_count = LowerSectorCffBuilder::new(parsed)
        .vector_matroid_components(&signatures)
        .into_iter()
        .filter(|component| {
            component.iter().any(|edge_id| {
                signatures[*edge_id]
                    .loop_signature
                    .iter()
                    .any(|coefficient| *coefficient != 0)
            })
        })
        .count();
    let mut active_signature_counts = BTreeMap::<(MomentumSignature, Option<String>), usize>::new();
    for edge_id in parsed.denominator_internal_edge_ids() {
        let signature = &signatures[edge_id];
        if !signature
            .loop_signature
            .iter()
            .any(|coefficient| *coefficient != 0)
        {
            continue;
        }
        let (canonical, _) = signature.canonical_up_to_sign();
        *active_signature_counts
            .entry((canonical, parsed.internal_edges[edge_id].mass_key.clone()))
            .or_default() += 1;
    }
    let duplicate_excess = active_signature_counts
        .values()
        .map(|count| count.saturating_sub(1))
        .sum::<usize>();
    denominator_rank.saturating_sub(rational_component_count) + duplicate_excess
}

fn denominator_contour_frame_sign(parsed: &ParsedGraph) -> CffGlobalPrefactorSign {
    CffGlobalPrefactorSign::from_exponent(denominator_contour_frame_exponent(parsed))
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

    #[test]
    fn known_factor_projection_preserves_survivors_and_reconstructs_only_deleted_variables() {
        let edge = |edge_id, tail, head, loop_signature| ParsedGraphInternalEdge {
            edge_id,
            tail,
            head,
            label: format!("q{edge_id}"),
            mass_key: Some(format!("m{edge_id}")),
            signature: MomentumSignature {
                loop_signature,
                external_signature: Vec::new(),
            },
            had_pow: false,
        };
        let parsed = ParsedGraph {
            internal_edges: vec![
                edge(0, 0, 1, vec![1, 1]),
                edge(1, 0, 2, vec![1, 0]),
                edge(2, 2, 1, vec![0, 1]),
                edge(3, 0, 1, vec![1, 1]),
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["k0".to_string(), "k1".to_string()],
            external_names: Vec::new(),
            node_name_to_internal: (0..3).map(|node| (format!("v{node}"), node)).collect(),
        };
        let builder = KnownFactorCffBuilder::new(
            &parsed,
            vec![0; parsed.internal_edges.len()],
            NumeratorSamplingScaleMode::None,
        );
        let (subparsed, sub_to_parent) =
            BoundedCffBuilder::for_bounds(&parsed, vec![0; parsed.internal_edges.len()])
                .project_parsed_edges(&[0]);

        assert_eq!(sub_to_parent, vec![1, 2, 3]);
        assert_eq!(
            builder
                .remap_known_factor_to_sub(
                    &parsed,
                    &subparsed,
                    &sub_to_parent,
                    &KnownLinearExpr::var(3, 1),
                )
                .unwrap(),
            KnownLinearExpr::var(2, 1),
            "a surviving occurrence must retain its own generalized sample",
        );
        assert_eq!(
            builder
                .remap_known_factor_to_sub(
                    &parsed,
                    &subparsed,
                    &sub_to_parent,
                    &KnownLinearExpr::var(0, 1),
                )
                .unwrap(),
            KnownLinearExpr::var(0, 1) + KnownLinearExpr::var(1, 1),
            "only the deleted occurrence is reconstructed in the lower-sector loop-coordinate basis",
        );

        let mut opposite = parsed.clone();
        opposite.internal_edges[3].signature.loop_signature = vec![-1, -1];
        let opposite_builder = KnownFactorCffBuilder::new(
            &opposite,
            vec![0; opposite.internal_edges.len()],
            NumeratorSamplingScaleMode::None,
        );
        let (opposite_subparsed, opposite_sub_to_parent) =
            BoundedCffBuilder::for_bounds(&opposite, vec![0; opposite.internal_edges.len()])
                .project_parsed_edges(&[0]);
        assert_eq!(
            opposite_builder
                .remap_known_factor_to_sub(
                    &opposite,
                    &opposite_subparsed,
                    &opposite_sub_to_parent,
                    &KnownLinearExpr::var(3, 1),
                )
                .unwrap(),
            KnownLinearExpr::var(2, 1),
            "a surviving -Q occurrence carries its sign in its cloned energy map",
        );
    }

    #[test]
    fn polynomial_division_keeps_scanning_after_a_nondivisible_leading_term() {
        let loop_term = KnownPolynomial::variable(KnownPolynomialVar::Loop(0));
        let ose = KnownPolynomial::variable(KnownPolynomialVar::Ose(0));
        let dividend = loop_term.clone() + ose.clone() * ose.clone();

        let (quotient, remainder) = dividend.div_rem(&ose).unwrap();

        assert_eq!(quotient, ose.clone());
        assert_eq!(remainder, loop_term);
        assert_eq!(quotient * ose + remainder, dividend);
    }

    #[test]
    fn inherited_full_rank_terminal_keeps_single_signed_below_residue() {
        let parsed = ParsedGraph {
            internal_edges: vec![ParsedGraphInternalEdge {
                edge_id: 0,
                tail: 0,
                head: 0,
                label: "q0".to_string(),
                mass_key: Some("m0".to_string()),
                signature: MomentumSignature {
                    loop_signature: vec![1],
                    external_signature: Vec::new(),
                },
                had_pow: false,
            }],
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["q0".to_string()],
            external_names: Vec::new(),
            node_name_to_internal: BTreeMap::from([("n0".to_string(), 0)]),
        };
        let expression = BoundedCffBuilder::for_bounds(&parsed, vec![1])
            .build_quadratic_recursive(true)
            .unwrap();

        let orientation = expression.orientations.iter().exactly_one().unwrap();
        let variant = orientation.variants.iter().exactly_one().unwrap();
        assert_eq!(variant.origin.as_deref(), Some("residue_basis"));
        assert_eq!(
            rational_from_coefficient(&variant.prefactor),
            Rational::from(-1)
        );
    }

    #[test]
    fn embedded_full_rank_factor_keeps_single_below_residue() {
        let parsed = ParsedGraph {
            internal_edges: vec![ParsedGraphInternalEdge {
                edge_id: 0,
                tail: 0,
                head: 0,
                label: "q0".to_string(),
                mass_key: Some("m0".to_string()),
                signature: MomentumSignature {
                    loop_signature: vec![1],
                    external_signature: Vec::new(),
                },
                had_pow: false,
            }],
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["q0".to_string()],
            external_names: Vec::new(),
            node_name_to_internal: BTreeMap::from([("n0".to_string(), 0)]),
        };
        let generate = |cff_generation_context| {
            generate_3d_expression_from_parsed_generated(
                &parsed,
                &Generate3DExpressionOptions {
                    cff_generation_context,
                    energy_degree_bounds: Some(vec![(0, 1)]),
                    ..Default::default()
                },
            )
            .unwrap()
        };

        let standalone = generate(CffGenerationContext::Standalone);
        let factorized = generate(CffGenerationContext::EmbeddedCffFactor);

        assert_eq!(standalone.expression.orientations.len(), 2);
        let factorized_orientation = factorized
            .expression
            .orientations
            .iter()
            .exactly_one()
            .unwrap();
        assert!(
            standalone
                .expression
                .orientations
                .iter()
                .any(|orientation| orientation.to_atom() == factorized_orientation.to_atom()),
            "the factorized Below residue must be one existing standalone routing, not a new normalization"
        );
        assert_eq!(
            factorized.core_global_prefactor_sign,
            standalone.core_global_prefactor_sign
        );
    }

    #[test]
    fn factorized_generalized_core_does_not_duplicate_denominator_parity() {
        let parsed = ParsedGraph {
            internal_edges: (0..2)
                .map(|edge_id| ParsedGraphInternalEdge {
                    edge_id,
                    tail: edge_id,
                    head: (edge_id + 1) % 2,
                    label: format!("q{edge_id}"),
                    mass_key: Some("m".to_string()),
                    signature: MomentumSignature {
                        loop_signature: vec![1],
                        external_signature: Vec::new(),
                    },
                    had_pow: false,
                })
                .collect(),
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["q".to_string()],
            external_names: Vec::new(),
            node_name_to_internal: BTreeMap::from([("n0".to_string(), 0), ("n1".to_string(), 1)]),
        };
        let generate = |energy_degree_bounds| {
            generate_3d_expression_from_parsed_generated(
                &parsed,
                &Generate3DExpressionOptions {
                    cff_generation_context: CffGenerationContext::EmbeddedCffFactor,
                    energy_degree_bounds: Some(energy_degree_bounds),
                    ..Default::default()
                },
            )
            .unwrap()
        };

        let constant = generate(Vec::new());
        let generalized = generate(vec![(0, 2)]);

        assert_eq!(cff_duplicate_signature_excess(&parsed), 1);
        assert_eq!(constant.core_global_prefactor_sign.factor(), -1);
        assert_eq!(
            generalized.core_global_prefactor_sign.factor(),
            1,
            "variant-local generalized residues already carry the repeated-channel parity",
        );
        assert_eq!(
            generalized.denominator_only_global_prefactor_sign.factor(),
            -1,
            "denominator-only source metadata remains independent of numerator ownership",
        );
    }

    #[cfg(feature = "eval")]
    fn powered_identity_test_graph(
        signatures: Vec<MomentumSignature>,
        loop_count: usize,
        external_count: usize,
    ) -> ParsedGraph {
        let edge_count = signatures.len();
        ParsedGraph {
            internal_edges: signatures
                .into_iter()
                .enumerate()
                .map(|(edge_id, signature)| ParsedGraphInternalEdge {
                    edge_id,
                    tail: edge_id,
                    head: (edge_id + 1) % edge_count,
                    label: format!("q{edge_id}"),
                    mass_key: Some("m".to_string()),
                    signature,
                    had_pow: false,
                })
                .collect(),
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: (0..loop_count)
                .map(|loop_id| format!("k{loop_id}"))
                .collect(),
            external_names: (0..external_count)
                .map(|external_id| format!("p{external_id}"))
                .collect(),
            node_name_to_internal: (0..edge_count)
                .map(|node| (format!("n{node}"), node))
                .collect(),
        }
    }

    #[cfg(feature = "eval")]
    fn powered_identity_test_generate(
        parsed: &ParsedGraph,
        degree: usize,
    ) -> GeneratedThreeDExpression {
        generate_3d_expression_from_parsed_generated(
            parsed,
            &Generate3DExpressionOptions {
                cff_generation_context: CffGenerationContext::Standalone,
                energy_degree_bounds: (degree != 0).then_some(vec![(0, degree)]),
                ..Default::default()
            },
        )
        .unwrap()
    }

    #[cfg(feature = "eval")]
    fn powered_identity_test_value(
        parsed: &ParsedGraph,
        generated: &GeneratedThreeDExpression,
        numerator: &str,
        input: &crate::eval::EvaluationInput,
    ) -> f64 {
        crate::eval::evaluate_expression(parsed, &generated.expression, numerator, input)
            .unwrap()
            .value
    }

    #[cfg(feature = "eval")]
    #[test]
    fn same_routing_unequal_mass_quadratic_matches_partial_fractions() {
        let signature = || MomentumSignature {
            loop_signature: vec![1],
            external_signature: Vec::new(),
        };
        let with_mass_keys = |mut parsed: ParsedGraph, mass_keys: &[&str]| {
            for (edge, mass_key) in parsed.internal_edges.iter_mut().zip(mass_keys) {
                edge.mass_key = Some((*mass_key).to_string());
            }
            parsed
        };
        let parent = with_mass_keys(
            powered_identity_test_graph(vec![signature(), signature()], 1, 0),
            &["ma", "mb"],
        );
        let pole_a = with_mass_keys(
            powered_identity_test_graph(vec![signature()], 1, 0),
            &["ma"],
        );
        let pole_b = with_mass_keys(
            powered_identity_test_graph(vec![signature()], 1, 0),
            &["mb"],
        );
        let spatial = [0.31, -0.47, 0.83];
        let mass_a: f64 = 0.73;
        let mass_b: f64 = 1.19;
        let energy_a_squared = spatial
            .iter()
            .map(|component| component * component)
            .sum::<f64>()
            + mass_a.powi(2);
        let energy_b_squared = spatial
            .iter()
            .map(|component| component * component)
            .sum::<f64>()
            + mass_b.powi(2);
        let input = |masses| crate::eval::EvaluationInput {
            external_momenta: Vec::new(),
            loop_spatial_momenta: vec![spatial],
            masses,
            uniform_scale: None,
        };
        let generate = |parsed: &ParsedGraph, cff_generation_context, energy_degree_bounds| {
            generate_3d_expression_from_parsed_generated(
                parsed,
                &Generate3DExpressionOptions {
                    cff_generation_context,
                    energy_degree_bounds,
                    ..Default::default()
                },
            )
            .unwrap()
        };

        let simple_a = generate(
            &pole_a,
            CffGenerationContext::EmbeddedCffFactor,
            Some(Vec::new()),
        );
        let simple_b = generate(
            &pole_b,
            CffGenerationContext::EmbeddedCffFactor,
            Some(Vec::new()),
        );
        let value_a = powered_identity_test_value(&pole_a, &simple_a, "1", &input(vec![mass_a]));
        let value_b = powered_identity_test_value(&pole_b, &simple_b, "1", &input(vec![mass_b]));
        // Each factorized one-pole expression is the public positive pole
        // convention.  The reduced pole inside the two-denominator parent
        // instead retains its signed Below contour, hence the overall minus.
        let expected = -(energy_a_squared * value_a - energy_b_squared * value_b)
            / (energy_a_squared - energy_b_squared);
        let mut failures = Vec::new();
        for context in [
            CffGenerationContext::Standalone,
            CffGenerationContext::EmbeddedCffFactor,
        ] {
            let quadratic = generate(&parent, context, Some(vec![(0, 2)]));
            let actual = powered_identity_test_value(
                &parent,
                &quadratic,
                "edges[0][0]**2",
                &input(vec![mass_a, mass_b]),
            );
            let scale = actual.abs().max(expected.abs()).max(f64::MIN_POSITIVE);
            if (actual - expected).abs() > 1.0e-12 * scale {
                failures.push(format!(
                    "{context:?}: actual={actual:.17e}, expected={expected:.17e}"
                ));
            }
        }
        assert!(
            failures.is_empty(),
            "same-routing unequal-mass quadratic partial fractions failed: {}",
            failures.join("; "),
        );
    }

    #[cfg(feature = "eval")]
    #[test]
    fn powered_denominator_lowering_preserves_linear_factor_and_routing() {
        let graph = |signs: &[i32]| {
            powered_identity_test_graph(
                signs
                    .iter()
                    .map(|sign| MomentumSignature {
                        loop_signature: vec![*sign],
                        external_signature: Vec::new(),
                    })
                    .collect(),
                1,
                0,
            )
        };
        let input = |edge_count| crate::eval::EvaluationInput {
            external_momenta: Vec::new(),
            loop_spatial_momenta: vec![[0.31, -0.47, 0.83]],
            masses: vec![0.73; edge_count],
            uniform_scale: None,
        };
        let denominator = "dot(edges[0],edges[0])-0.5329";
        let mut reference: Option<f64> = None;
        for (parent_signs, child_signs) in [
            (&[1, 1, 1][..], &[1, 1][..]),
            (&[1, 1, -1][..], &[1, 1][..]),
            (&[-1, 1, 1][..], &[-1, 1][..]),
            (&[1, -1, 1][..], &[1, -1][..]),
        ] {
            let parent = graph(parent_signs);
            let child = graph(child_signs);
            let parent_generated = powered_identity_test_generate(&parent, 3);
            let child_generated = powered_identity_test_generate(&child, 1);
            let physical_q0 = |sign: i32| {
                if sign == 1 {
                    "edges[0][0]+0.37".to_string()
                } else {
                    "-edges[0][0]+0.37".to_string()
                }
            };
            let parent_numerator = format!("({denominator})*({})", physical_q0(parent_signs[0]));
            let child_numerator = physical_q0(child_signs[0]);
            let parent_value = powered_identity_test_value(
                &parent,
                &parent_generated,
                &parent_numerator,
                &input(parent_signs.len()),
            );
            let child_value = powered_identity_test_value(
                &child,
                &child_generated,
                &child_numerator,
                &input(child_signs.len()),
            );
            let scale = parent_value
                .abs()
                .max(child_value.abs())
                .max(f64::MIN_POSITIVE);
            assert!(
                (parent_value - child_value).abs() <= 1.0e-11 * scale,
                "D(Q)*(Q0+c)/D(Q)^3 changed under denominator lowering or routing {parent_signs:?}: parent={parent_value:e}, child={child_value:e}",
            );
            if let Some(reference) = reference {
                assert!(
                    (parent_value - reference).abs() <= 1.0e-11 * scale,
                    "the physical powered quotient changed under Q/-Q occurrence reordering: parent={parent_value:e}, reference={reference:e}",
                );
            } else {
                reference = Some(parent_value);
            }
        }
    }

    #[cfg(feature = "eval")]
    #[test]
    fn powered_denominator_lowering_preserves_quadratic_retained_factor() {
        let graph = |power| {
            powered_identity_test_graph(
                (0..power)
                    .map(|_| MomentumSignature {
                        loop_signature: vec![1],
                        external_signature: Vec::new(),
                    })
                    .collect(),
                1,
                0,
            )
        };
        let parent = graph(3);
        let child = graph(2);
        let parent_generated = powered_identity_test_generate(&parent, 4);
        let child_generated = powered_identity_test_generate(&child, 2);
        let parent_input = crate::eval::EvaluationInput {
            external_momenta: Vec::new(),
            loop_spatial_momenta: vec![[-0.59, 0.23, 1.07]],
            masses: vec![0.73; 3],
            uniform_scale: None,
        };
        let child_input = crate::eval::EvaluationInput {
            masses: vec![0.73; 2],
            ..parent_input.clone()
        };
        let retained = "(edges[0][0]+0.37)**2";
        let parent_value = powered_identity_test_value(
            &parent,
            &parent_generated,
            &format!("(dot(edges[0],edges[0])-0.5329)*{retained}"),
            &parent_input,
        );
        let child_value =
            powered_identity_test_value(&child, &child_generated, retained, &child_input);
        let scale = parent_value
            .abs()
            .max(child_value.abs())
            .max(f64::MIN_POSITIVE);
        assert!(
            (parent_value - child_value).abs() <= 1.0e-11 * scale,
            "the DOD-2 retained quadratic changed under D(Q)/D(Q)^3 -> 1/D(Q)^2: parent={parent_value:e}, child={child_value:e}",
        );
    }

    #[cfg(feature = "eval")]
    #[test]
    fn quadratic_repeated_channel_matches_half_simple_pole_under_physical_routing_reversal() {
        let signature = |sign| MomentumSignature {
            loop_signature: vec![sign],
            external_signature: Vec::new(),
        };
        let spatial = [0.11_f64, -0.17, 0.23];
        let mass = 0.73_f64;
        let input = |edge_count| crate::eval::EvaluationInput {
            external_momenta: Vec::new(),
            loop_spatial_momenta: vec![spatial],
            masses: vec![mass; edge_count],
            uniform_scale: None,
        };
        let simple = powered_identity_test_graph(vec![signature(1)], 1, 0);
        let simple_generated = generate_3d_expression_from_parsed_generated(
            &simple,
            &Generate3DExpressionOptions {
                cff_generation_context: CffGenerationContext::EmbeddedCffFactor,
                energy_degree_bounds: Some(Vec::new()),
                ..Default::default()
            },
        )
        .unwrap();
        let expected =
            powered_identity_test_value(&simple, &simple_generated, "1", &input(1)) / 2.0;

        for signs in [[1, -1], [-1, 1], [1, 1], [-1, -1]] {
            let mut parsed =
                powered_identity_test_graph(signs.into_iter().map(signature).collect(), 1, 0);
            for (edge, sign) in parsed.internal_edges.iter_mut().zip(signs) {
                if sign == -1 {
                    std::mem::swap(&mut edge.tail, &mut edge.head);
                }
            }
            assert!(crate::validate_parsed_graph(&parsed).ok);
            for context in [
                CffGenerationContext::Standalone,
                CffGenerationContext::EmbeddedCffFactor,
            ] {
                let generated = generate_3d_expression_from_parsed_generated(
                    &parsed,
                    &Generate3DExpressionOptions {
                        cff_generation_context: context,
                        energy_degree_bounds: Some(vec![(0, 2)]),
                        ..Default::default()
                    },
                )
                .unwrap();
                let actual =
                    powered_identity_test_value(&parsed, &generated, "edges[0][0]**2", &input(2));
                let scale = actual.abs().max(expected.abs()).max(f64::MIN_POSITIVE);
                assert!(
                    (actual - expected).abs() <= 1.0e-12 * scale,
                    "Q0^2/D(Q)^2 changed under the valid routing {signs:?} in {context:?}: actual={actual:e}, expected={expected:e}",
                );
            }
        }
    }

    #[cfg(feature = "eval")]
    #[test]
    fn quartic_repeated_channel_is_invariant_under_physical_routing_reversal() {
        let graph = |signs: [i32; 3]| {
            let mut parsed = powered_identity_test_graph(
                signs
                    .iter()
                    .map(|sign| MomentumSignature {
                        loop_signature: vec![*sign],
                        external_signature: Vec::new(),
                    })
                    .collect(),
                1,
                0,
            );
            // Reversing Q -> -Q for a physical propagator also reverses its
            // incidence. Keeping these two pieces synchronized is essential:
            // provenance/routing changes may not manufacture a graph which
            // violates momentum conservation at its vertices.
            for (edge, sign) in parsed.internal_edges.iter_mut().zip(signs) {
                if sign == -1 {
                    std::mem::swap(&mut edge.tail, &mut edge.head);
                }
            }
            assert!(crate::validate_parsed_graph(&parsed).ok);
            parsed
        };
        let spatial = [0.11_f64, -0.17, 0.23];
        let mass = 0.73_f64;
        let energy = (spatial
            .iter()
            .map(|component| component * component)
            .sum::<f64>()
            + mass * mass)
            .sqrt();
        let expected = 3.0 / (16.0 * energy);
        let mut reference: Option<f64> = None;

        for signs in [[1, 1, 1], [-1, 1, 1], [1, -1, -1], [-1, -1, -1]] {
            let parsed = graph(signs);
            let input = crate::eval::EvaluationInput {
                external_momenta: Vec::new(),
                loop_spatial_momenta: vec![spatial],
                masses: vec![mass; 3],
                uniform_scale: None,
            };
            for context in [
                CffGenerationContext::Standalone,
                CffGenerationContext::EmbeddedCffFactor,
            ] {
                let generated = generate_3d_expression_from_parsed_generated(
                    &parsed,
                    &Generate3DExpressionOptions {
                        cff_generation_context: context,
                        energy_degree_bounds: Some(vec![(0, 4)]),
                        ..Default::default()
                    },
                )
                .unwrap();
                let actual =
                    powered_identity_test_value(&parsed, &generated, "edges[0][0]**4", &input);
                let scale = actual.abs().max(expected.abs()).max(f64::MIN_POSITIVE);
                assert!(
                    (actual - expected).abs() <= 1.0e-12 * scale,
                    "Q0^4/D(Q)^3 changed under the valid routing {signs:?} in {context:?}: actual={actual:e}, expected={expected:e}",
                );
                if let Some(reference) = reference {
                    assert!(
                        (actual - reference).abs() <= 1.0e-12 * scale,
                        "Q0^4/D(Q)^3 is not routing/context invariant: actual={actual:e}, reference={reference:e}",
                    );
                } else {
                    reference = Some(actual);
                }
            }
        }
    }

    #[cfg(feature = "eval")]
    #[test]
    fn shifted_spectator_powered_quotient_matches_lower_source() {
        let signature = |loop_signature, external_signature| MomentumSignature {
            loop_signature,
            external_signature,
        };
        let parent = powered_identity_test_graph(
            vec![
                signature(vec![1], vec![0]),
                signature(vec![1], vec![0]),
                signature(vec![1], vec![1, 0]),
                signature(vec![1], vec![0, 1]),
            ],
            1,
            2,
        );
        let child = powered_identity_test_graph(
            vec![
                signature(vec![1], vec![1, 0]),
                signature(vec![1], vec![0, 1]),
            ],
            1,
            2,
        );
        let parent_generated = powered_identity_test_generate(&parent, 4);
        let child_generated = powered_identity_test_generate(&child, 0);
        let parent_input = crate::eval::EvaluationInput {
            external_momenta: vec![[1.4, 0.2, -0.1, 0.3], [0.9, -0.1, 0.4, -0.2]],
            loop_spatial_momenta: vec![[0.31, -0.47, 0.83]],
            masses: vec![0.73; 4],
            uniform_scale: None,
        };
        let child_input = crate::eval::EvaluationInput {
            masses: vec![0.73; 2],
            ..parent_input.clone()
        };
        let parent_value = powered_identity_test_value(
            &parent,
            &parent_generated,
            "(dot(edges[0],edges[0])-0.5329)**2",
            &parent_input,
        );
        let child_value = powered_identity_test_value(&child, &child_generated, "1", &child_input);
        let scale = parent_value
            .abs()
            .max(child_value.abs())
            .max(f64::MIN_POSITIVE);
        assert!(
            (parent_value - child_value).abs() <= 1.0e-11 * scale,
            "D(q)^2/[D(q)^2 D(q+p) D(q+r)] did not reduce to 1/[D(q+p) D(q+r)]: parent={parent_value:e}, child={child_value:e}",
        );
    }

    #[cfg(feature = "eval")]
    #[test]
    fn denominator_frame_metadata_uses_active_rank_after_loop_direction_deletion() {
        let signature = |loop_signature| MomentumSignature {
            loop_signature,
            external_signature: Vec::new(),
        };
        let parent = powered_identity_test_graph(
            vec![
                signature(vec![1, 0]),
                signature(vec![1, 0]),
                signature(vec![0, 1]),
            ],
            2,
            0,
        );
        let child =
            powered_identity_test_graph(vec![signature(vec![1, 0]), signature(vec![1, 0])], 2, 0);

        assert_eq!(denominator_contour_frame_sign(&parent).factor(), -1);
        assert_eq!(
            denominator_contour_frame_sign(&child).factor(),
            -1,
            "a projected child must not count a consumed ambient loop direction in its denominator metadata",
        );
    }

    #[test]
    fn embedded_nonterminal_factor_keeps_full_standalone_residue_sum() {
        let parsed = crate::graph_io::test_graphs::box_graph();
        let generate = |cff_generation_context| {
            generate_3d_expression_from_parsed_generated(
                &parsed,
                &Generate3DExpressionOptions {
                    cff_generation_context,
                    energy_degree_bounds: Some(vec![(0, 1)]),
                    ..Default::default()
                },
            )
            .unwrap()
        };

        let standalone = generate(CffGenerationContext::Standalone);
        let factorized = generate(CffGenerationContext::EmbeddedCffFactor);

        assert!(standalone.expression.orientations.len() > 1);
        assert_eq!(
            factorized.expression.orientations.len(),
            standalone.expression.orientations.len()
        );
        for (standalone, factorized) in standalone
            .expression
            .orientations
            .iter()
            .zip(&factorized.expression.orientations)
        {
            assert_eq!(factorized.data.orientation, standalone.data.orientation);
            assert_eq!(factorized.loop_energy_map, standalone.loop_energy_map);
            assert_eq!(factorized.edge_energy_map, standalone.edge_energy_map);
            assert_eq!(factorized.to_atom(), standalone.to_atom());
        }
    }

    #[test]
    fn two_loop_full_rank_terminal_keeps_all_below_residue_sign() {
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
                    tail: 0,
                    head: 0,
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
            node_name_to_internal: BTreeMap::from([("n0".to_string(), 0)]),
        };
        let expression = BoundedCffBuilder::for_bounds(&parsed, vec![1, 1])
            .build_quadratic_recursive(true)
            .unwrap();

        let orientation = expression.orientations.iter().exactly_one().unwrap();
        let variant = orientation.variants.iter().exactly_one().unwrap();
        assert_eq!(variant.origin.as_deref(), Some("residue_basis"));
        assert_eq!(
            rational_from_coefficient(&variant.prefactor),
            Rational::one(),
            "the two-loop all-Below residue carries (-1)^L, one minus relative to the pure-CFF core (-1)^(L-1)",
        );
    }

    #[test]
    fn projected_contact_preserves_unrelated_quadratic_source_bound() {
        let edge = |edge_id, tail, head, loop_signature, mass: &str| ParsedGraphInternalEdge {
            edge_id,
            tail,
            head,
            label: format!("q{edge_id}"),
            mass_key: Some(mass.to_string()),
            signature: MomentumSignature {
                loop_signature,
                external_signature: Vec::new(),
            },
            had_pow: false,
        };
        // Momentum conservation makes this a local two-loop graph: q0=z,
        // q1=z+q and q2,q3,q4=q. The factorized numerator q0^2*q2^2 has DOD
        // -1 in the z subcycle and -5 in the q subcycle. Projecting q0 may
        // consume only its own quadratic capacity; q2's pre-existing source
        // bound must still generate the deterministic zero contact sample.
        let parsed = ParsedGraph {
            internal_edges: vec![
                edge(0, 0, 1, vec![1, 0], "mz"),
                edge(1, 1, 0, vec![1, 1], "ms"),
                edge(2, 0, 2, vec![0, 1], "mq0"),
                edge(3, 2, 3, vec![0, 1], "mq1"),
                edge(4, 3, 1, vec![0, 1], "mq2"),
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["z".to_string(), "q".to_string()],
            external_names: Vec::new(),
            node_name_to_internal: (0..4).map(|node| (format!("v{node}"), node)).collect(),
        };
        let expression = BoundedCffBuilder::for_bounds(&parsed, vec![2, 0, 2, 0, 0])
            .build()
            .unwrap();

        assert!(expression.orientations.iter().any(|orientation| {
            orientation.edge_energy_map[0].is_zero()
                && orientation.edge_energy_map[2].is_zero()
                && orientation.variants.iter().any(|variant| {
                    variant.origin.as_deref()
                        == Some("bounded_degree_quadratic_recursive_contact:e0=0")
                })
        }));
    }

    #[test]
    fn nested_inherited_contours_use_latest_closure_on_a_shared_terminal_pivot() {
        let edge = |edge_id, loop_signature| ParsedGraphInternalEdge {
            edge_id,
            tail: 0,
            head: 0,
            label: format!("q{edge_id}"),
            mass_key: Some(format!("m{edge_id}")),
            signature: MomentumSignature {
                loop_signature,
                external_signature: Vec::new(),
            },
            had_pow: false,
        };
        let parsed = ParsedGraph {
            internal_edges: vec![edge(0, vec![1, 1]), edge(1, vec![0, 1])],
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["z".to_string(), "q".to_string()],
            external_names: Vec::new(),
            node_name_to_internal: BTreeMap::from([("v0".to_string(), 0)]),
        };
        let terminal_component = |inherited_contour_rows| {
            let mut builder = LowerSectorCffBuilder::new(&parsed);
            builder.inherited_contour_rows = inherited_contour_rows;
            let signatures = builder.signatures();
            builder
                .component_bundles(&signatures)
                .unwrap()
                .into_iter()
                .find(|component| component.local_to_sub == vec![1])
                .unwrap()
        };

        // z=(z+q)-q first closes the q pivot Above. The later, inner q row
        // closes that same pivot Below and must shadow it because the outer
        // residue/Jacobian has already been consumed by the parent contact.
        let inner_below = terminal_component(vec![vec![1, 0], vec![0, 1]]);
        let inner_above = terminal_component(vec![vec![0, 1], vec![1, 0]]);
        assert!(inner_below.requires_inherited_component_bridge);
        assert!(inner_above.requires_inherited_component_bridge);
        assert_eq!(
            inner_below
                .expression
                .orientations
                .iter()
                .exactly_one()
                .unwrap()
                .edge_energy_map[0],
            LinearEnergyExpr::ose(EdgeIndex(0), 1),
        );
        assert_eq!(
            inner_above
                .expression
                .orientations
                .iter()
                .exactly_one()
                .unwrap()
                .edge_energy_map[0],
            LinearEnergyExpr::ose(EdgeIndex(0), -1),
        );

        let component_bridges = |parsed: &ParsedGraph, inherited_contour_rows| {
            let mut builder = LowerSectorCffBuilder::new(parsed);
            builder.inherited_contour_rows = inherited_contour_rows;
            let signatures = builder.signatures();
            builder
                .component_bundles(&signatures)
                .unwrap()
                .into_iter()
                .map(|component| component.requires_inherited_component_bridge)
                .collect::<Vec<_>>()
        };
        let mut powered = parsed.clone();
        powered.internal_edges.push(edge(2, vec![0, 1]));
        assert!(
            component_bridges(&powered, vec![vec![0, 1]])
                .into_iter()
                .any(|bridges| bridges),
            "an outermost deleted contour with a surviving powered carrier retains the parent component-product frame",
        );
        assert!(
            component_bridges(&powered, vec![vec![0, -1]])
                .into_iter()
                .any(|bridges| bridges),
            "reversing the carrier changes its closure and Jacobian, not whether the component-product frame survives",
        );
        assert!(
            component_bridges(&powered, vec![vec![1, 0]])
                .into_iter()
                .any(|bridges| bridges),
            "a row representable only across components needs their product bridge",
        );
        let mut no_carrier = powered;
        no_carrier.loop_names.push("u".to_string());
        for edge in &mut no_carrier.internal_edges {
            edge.signature.loop_signature.push(0);
        }
        assert!(
            !component_bridges(&no_carrier, vec![vec![0, 0, 1]])
                .into_iter()
                .any(|bridges| bridges),
            "a deleted row with no surviving carrier creates no component anchor",
        );
        assert!(
            !component_bridges(&no_carrier, vec![vec![0, 0, 1], vec![0, 1, 0]])
                .into_iter()
                .any(|bridges| bridges),
            "an inner surviving carrier cannot reopen a product frame consumed by the outer contact",
        );
    }

    #[test]
    fn terminal_duplicate_contact_sign_is_routing_and_order_invariant() {
        let expected = BTreeMap::from([
            ("+".to_string(), Rational::from(2)),
            ("-".to_string(), Rational::from(2)),
            ("0".to_string(), Rational::from(-4)),
        ]);

        for (label, edge_data, active_edge) in [
            ("Q/Q", [(0, 1, 1), (1, 0, 1)], 0),
            ("-Q/-Q reordered", [(1, 0, -1), (0, 1, -1)], 1),
            ("Q/-Q", [(0, 1, 1), (0, 1, -1)], 0),
            ("-Q/Q reordered", [(0, 1, -1), (0, 1, 1)], 1),
        ] {
            let parsed = ParsedGraph {
                internal_edges: edge_data
                    .into_iter()
                    .enumerate()
                    .map(|(edge_id, (tail, head, sign))| ParsedGraphInternalEdge {
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
                    })
                    .collect(),
                external_edges: Vec::new(),
                initial_state_cut_edges: Vec::new(),
                loop_names: vec!["q".to_string()],
                external_names: Vec::new(),
                node_name_to_internal: BTreeMap::from([
                    ("n0".to_string(), 0),
                    ("n1".to_string(), 1),
                ]),
            };
            let mut bounds = vec![0; 2];
            bounds[active_edge] = 2;
            let expression = BoundedCffBuilder::for_bounds(&parsed, bounds)
                .build_quadratic_recursive(false)
                .unwrap();
            let contact_variants = expression
                .orientations
                .iter()
                .flat_map(|orientation| &orientation.variants)
                .filter(|variant| {
                    variant.origin.as_deref().is_some_and(|origin| {
                        origin.starts_with("bounded_degree_quadratic_recursive_contact")
                    })
                })
                .collect::<Vec<_>>();
            assert_eq!(
                contact_variants.len(),
                3,
                "unexpected contact count for {label}"
            );
            let coefficients = contact_variants
                .into_iter()
                .map(|variant| {
                    let origin = variant.origin.as_deref().unwrap();
                    let sample = origin.rsplit_once('=').unwrap().1.to_string();
                    (sample, rational_from_coefficient(&variant.prefactor))
                })
                .collect::<BTreeMap<_, _>>();
            assert_eq!(
                coefficients, expected,
                "terminal duplicate contact sign changed for {label}"
            );
        }
    }

    #[test]
    fn quadratic_hermite_remainder_retains_powered_denominator_parity() {
        let repeated_cycle = |power: usize| ParsedGraph {
            internal_edges: (0..power)
                .map(|edge_id| ParsedGraphInternalEdge {
                    edge_id,
                    tail: edge_id,
                    head: (edge_id + 1) % power,
                    label: format!("q{edge_id}"),
                    mass_key: Some("m".to_string()),
                    signature: MomentumSignature {
                        loop_signature: vec![1],
                        external_signature: Vec::new(),
                    },
                    had_pow: false,
                })
                .collect(),
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["q".to_string()],
            external_names: Vec::new(),
            node_name_to_internal: (0..power).map(|node| (format!("n{node}"), node)).collect(),
        };
        let powered = repeated_cycle(2);
        let expression = BoundedCffBuilder::for_bounds(&powered, vec![2, 0])
            .build_quadratic_recursive(false)
            .unwrap();
        let sector = |prefix: &str| {
            let mut selected = expression.clone();
            for orientation in &mut selected.orientations {
                orientation.variants.retain(|variant| {
                    variant
                        .origin
                        .as_deref()
                        .is_some_and(|origin| origin.starts_with(prefix))
                });
            }
            selected
                .orientations
                .retain(|orientation| !orientation.variants.is_empty());
            selected
        };
        let input = crate::eval::EvaluationInput {
            external_momenta: Vec::new(),
            loop_spatial_momenta: vec![[0.31, -0.47, 0.83]],
            masses: vec![0.73; 2],
            uniform_scale: None,
        };
        let numerator = "edges[0][0]**2";
        let contact = crate::eval::evaluate_expression(
            &powered,
            &sector("bounded_degree_quadratic_recursive_contact"),
            numerator,
            &input,
        )
        .unwrap()
        .value;
        let remainder = crate::eval::evaluate_expression(
            &powered,
            &sector("bounded_degree_quadratic_recursive_remainder"),
            numerator,
            &input,
        )
        .unwrap()
        .value;
        let total = crate::eval::evaluate_expression(&powered, &expression, numerator, &input)
            .unwrap()
            .value;

        // x^2 = E^2 + (x^2-E^2): the lower contact is 1/D, while the
        // repeated-denominator remainder E^2/D^2 contributes -1/2 of it.
        // Consequently x^2/D^2 is exactly one half of the lower contact.
        let scale = contact.abs().max(f64::MIN_POSITIVE);
        assert!((remainder + contact / 2.0).abs() <= 1.0e-13 * scale);
        assert!((total - contact / 2.0).abs() <= 1.0e-13 * scale);
    }

    #[test]
    fn duplicate_contact_records_energy_factor_ownership_without_a_second_frame() {
        let repeated_cycle = |power: usize| ParsedGraph {
            internal_edges: (0..power)
                .map(|edge_id| ParsedGraphInternalEdge {
                    edge_id,
                    tail: edge_id,
                    head: (edge_id + 1) % power,
                    label: format!("q{edge_id}"),
                    mass_key: Some("m".to_string()),
                    signature: MomentumSignature {
                        loop_signature: vec![1],
                        external_signature: Vec::new(),
                    },
                    had_pow: false,
                })
                .collect(),
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["q".to_string()],
            external_names: Vec::new(),
            node_name_to_internal: (0..power).map(|node| (format!("n{node}"), node)).collect(),
        };
        let parent = generate_3d_expression_from_parsed_generated(
            &repeated_cycle(2),
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(0, 2)]),
                ..Default::default()
            },
        )
        .unwrap();
        let lower = generate_3d_expression_from_parsed_generated(
            &repeated_cycle(1),
            &Generate3DExpressionOptions::default(),
        )
        .unwrap();

        assert_eq!(
            parent.energy_factor_ownership,
            CffEnergyFactorOwnership::VariantLocal
        );
        assert_eq!(
            lower.energy_factor_ownership,
            CffEnergyFactorOwnership::GlobalSourceProduct
        );
        assert_eq!(parent.energy_factor_components.len(), 1);
        assert_eq!(lower.energy_factor_components.len(), 1);
    }

    #[cfg(feature = "eval")]
    #[test]
    fn nonterminal_quadratic_contact_matches_independent_lower_cff_frame() {
        let edge = |edge_id, tail, head, loop_signature, external_shift| ParsedGraphInternalEdge {
            edge_id,
            tail,
            head,
            label: format!("q{edge_id}"),
            mass_key: Some("m".to_string()),
            signature: MomentumSignature {
                loop_signature,
                external_signature: vec![external_shift],
            },
            had_pow: false,
        };
        // This is the smallest raw rational skeleton of the GL0 cograph:
        // q2 and q3 are the repeated outer channel, while q4 is its affine
        // spectator. The q2^2 contact removes one occurrence but leaves a
        // connected, full-rank, nonterminal two-loop lower sector.
        let parsed = ParsedGraph {
            internal_edges: vec![
                edge(0, 0, 1, vec![1, 0], 0),
                edge(1, 0, 1, vec![0, 1], 0),
                edge(2, 3, 0, vec![1, 1], 0),
                edge(3, 1, 2, vec![1, 1], 0),
                edge(4, 2, 3, vec![1, 1], -1),
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["k0".to_string(), "k1".to_string()],
            external_names: vec!["p".to_string()],
            node_name_to_internal: (0..4).map(|node| (format!("v{node}"), node)).collect(),
        };
        let bounds = vec![(2, 2)];
        let quadratic = generate_3d_expression_from_parsed_generated(
            &parsed,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(bounds),
                ..Default::default()
            },
        )
        .unwrap();
        let scalar = generate_3d_expression_from_parsed_generated(
            &parsed,
            &Generate3DExpressionOptions::default(),
        )
        .unwrap();
        let (lower, lower_to_orig) =
            BoundedCffBuilder::for_bounds(&parsed, vec![0; parsed.internal_edges.len()])
                .project_parsed_edges(&[2]);
        let lower_scalar = generate_3d_expression_from_parsed_generated(
            &lower,
            &Generate3DExpressionOptions::default(),
        )
        .unwrap();
        let unit_parent = ParsedGraph {
            internal_edges: vec![
                edge(0, 0, 2, vec![1], 0),
                edge(1, 1, 0, vec![1], -1),
                edge(2, 2, 1, vec![1], 0),
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["q".to_string()],
            external_names: vec!["p".to_string()],
            node_name_to_internal: (0..3).map(|node| (format!("u{node}"), node)).collect(),
        };
        let unit_quadratic = generate_3d_expression_from_parsed_generated(
            &unit_parent,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(0, 2)]),
                ..Default::default()
            },
        )
        .unwrap();
        let unit_scalar = generate_3d_expression_from_parsed_generated(
            &unit_parent,
            &Generate3DExpressionOptions::default(),
        )
        .unwrap();
        let (unit_lower, unit_lower_to_orig) =
            BoundedCffBuilder::for_bounds(&unit_parent, vec![0; 3]).project_parsed_edges(&[0]);
        let unit_lower_scalar = generate_3d_expression_from_parsed_generated(
            &unit_lower,
            &Generate3DExpressionOptions::default(),
        )
        .unwrap();

        let rank = |graph: &ParsedGraph| {
            rank_i64(
                &graph
                    .denominator_internal_edge_ids()
                    .iter()
                    .map(|edge_id| {
                        graph.internal_edges[*edge_id]
                            .signature
                            .loop_signature
                            .iter()
                            .map(|coefficient| i64::from(*coefficient))
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>(),
            )
        };
        let rational_component_count = |graph: &ParsedGraph| {
            let signatures = graph
                .internal_edges
                .iter()
                .map(|edge| edge.signature.clone())
                .collect::<Vec<_>>();
            LowerSectorCffBuilder::new(graph)
                .vector_matroid_components(&signatures)
                .len()
        };
        assert_eq!(repeated_groups(&parsed)[0].edge_ids, [2, 3]);
        assert_eq!(repeated_groups(&unit_parent)[0].edge_ids, [0, 2]);
        assert_eq!(cff_duplicate_signature_excess(&parsed), 1);
        assert_eq!(cff_duplicate_signature_excess(&lower), 0);
        assert_eq!(cff_duplicate_signature_excess(&unit_parent), 1);
        assert_eq!(cff_duplicate_signature_excess(&unit_lower), 0);
        assert_eq!((rank(&parsed), rank(&lower)), (2, 2));
        assert_eq!(
            (
                denominator_connected_components(&parsed).len(),
                denominator_connected_components(&lower).len(),
                rational_component_count(&parsed),
                rational_component_count(&lower),
            ),
            (1, 1, 1, 1),
            "pinching preserves rank and component structure; only the lower-CFF constructor frame changes"
        );

        let expression_value =
            |generated: &GeneratedThreeDExpression,
             graph: &ParsedGraph,
             numerator: &str,
             input: &crate::eval::EvaluationInput| {
                crate::eval::evaluate_expression(graph, &generated.expression, numerator, input)
                    .unwrap()
                    .value
            };
        let energy_squared = |graph: &ParsedGraph,
                              edge_id: usize,
                              input: &crate::eval::EvaluationInput| {
            let active = &graph.internal_edges[edge_id].signature;
            input.masses[edge_id].powi(2)
                + (0..3)
                    .map(|axis| {
                        active
                            .loop_signature
                            .iter()
                            .zip(&input.loop_spatial_momenta)
                            .map(|(coefficient, momentum)| f64::from(*coefficient) * momentum[axis])
                            .sum::<f64>()
                            + active
                                .external_signature
                                .iter()
                                .zip(&input.external_momenta)
                                .map(|(coefficient, momentum)| {
                                    f64::from(*coefficient) * momentum[axis + 1]
                                })
                                .sum::<f64>()
                    })
                    .map(|component| component * component)
                    .sum::<f64>()
        };
        for seed in [11, 29, 47] {
            let input = crate::eval::EvaluationInput::deterministic(
                &parsed,
                seed,
                &BTreeMap::from([("m".to_string(), 0.61)]),
                None,
            )
            .unwrap();
            let lower_input = crate::eval::EvaluationInput {
                external_momenta: input.external_momenta.clone(),
                loop_spatial_momenta: input.loop_spatial_momenta.clone(),
                masses: lower_to_orig
                    .iter()
                    .map(|orig_id| input.masses[*orig_id])
                    .collect(),
                uniform_scale: None,
            };
            let active_energy_squared = energy_squared(&parsed, 2, &input);
            // Each generated expression already contains the global sign
            // recorded by `core_global_prefactor_sign`. Multiplying the
            // independently generated children by that metadata here would
            // apply different signs to the two pieces of one Hermite identity.
            // Compare their embedded expressions directly; a consumer may
            // strip or replace one common parent sign only after recombination.
            let actual = expression_value(&quadratic, &parsed, "edges[2][0]**2", &input);
            let expected = active_energy_squared * expression_value(&scalar, &parsed, "1", &input)
                + expression_value(&lower_scalar, &lower, "1", &lower_input);
            let scale = actual.abs().max(expected.abs()).max(f64::MIN_POSITIVE);
            assert!(
                (actual - expected).abs() <= 1.0e-11 * scale,
                "quadratic repeated-channel contact violates q0^2=D(q)+E^2 at seed {seed}: actual={actual:.17e}, expected={expected:.17e}"
            );
            let unit_input = crate::eval::EvaluationInput::deterministic(
                &unit_parent,
                seed,
                &BTreeMap::from([("m".to_string(), 0.61)]),
                None,
            )
            .unwrap();
            let unit_lower_input = crate::eval::EvaluationInput {
                external_momenta: unit_input.external_momenta.clone(),
                loop_spatial_momenta: unit_input.loop_spatial_momenta.clone(),
                masses: unit_lower_to_orig
                    .iter()
                    .map(|orig_id| unit_input.masses[*orig_id])
                    .collect(),
                uniform_scale: None,
            };
            let unit_actual =
                expression_value(&unit_quadratic, &unit_parent, "edges[0][0]**2", &unit_input);
            let unit_expected = energy_squared(&unit_parent, 0, &unit_input)
                * expression_value(&unit_scalar, &unit_parent, "1", &unit_input)
                + expression_value(&unit_lower_scalar, &unit_lower, "1", &unit_lower_input);
            let unit_scale = unit_actual
                .abs()
                .max(unit_expected.abs())
                .max(f64::MIN_POSITIVE);
            assert!(
                (unit_actual - unit_expected).abs() <= 1.0e-11 * unit_scale,
                "one-loop quadratic contact violates q0^2=D(q)+E^2 at seed {seed}: actual={unit_actual:.17e}, expected={unit_expected:.17e}"
            );
        }
    }

    #[test]
    fn lower_sector_constructors_return_direct_expressions_for_all_bases() {
        let edge = |edge_id, tail, head, loop_signature, external_signature, mass: &str| {
            ParsedGraphInternalEdge {
                edge_id,
                tail,
                head,
                label: format!("q{edge_id}"),
                mass_key: Some(mass.to_string()),
                signature: MomentumSignature {
                    loop_signature,
                    external_signature,
                },
                had_pow: false,
            }
        };
        let pure_repeated = ParsedGraph {
            internal_edges: vec![
                edge(0, 0, 1, vec![1], Vec::new(), "m0"),
                edge(1, 1, 0, vec![-1], Vec::new(), "m0"),
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["k".to_string()],
            external_names: Vec::new(),
            node_name_to_internal: BTreeMap::from([("n0".to_string(), 0), ("n1".to_string(), 1)]),
        };
        let pure_component = LowerSectorCffBuilder::new(&pure_repeated).build().unwrap();
        assert!(!pure_component.orientations.is_empty());

        let affine = ParsedGraph {
            internal_edges: vec![
                edge(0, 0, 1, vec![1], vec![0], "m0"),
                edge(1, 1, 0, vec![-1], vec![0], "m0"),
                edge(2, 0, 1, vec![1], vec![-1], "m1"),
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["k".to_string()],
            external_names: vec!["p".to_string()],
            node_name_to_internal: pure_repeated.node_name_to_internal.clone(),
        };
        let affine = LowerSectorCffBuilder::new(&affine).build().unwrap();
        assert!(!affine.orientations.is_empty());

        let disconnected = ParsedGraph {
            internal_edges: vec![
                edge(0, 0, 1, vec![1, 0], vec![0], "m0"),
                edge(1, 1, 0, vec![-1, 0], vec![0], "m0"),
                edge(2, 2, 3, vec![0, 1], vec![0], "m1"),
                edge(3, 3, 2, vec![0, 1], vec![-1], "m2"),
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["k0".to_string(), "k1".to_string()],
            external_names: vec!["p".to_string()],
            node_name_to_internal: (0..4).map(|node| (format!("n{node}"), node)).collect(),
        };
        let disconnected = LowerSectorCffBuilder::new(&disconnected).build().unwrap();
        assert!(!disconnected.orientations.is_empty());
    }

    #[test]
    fn inherited_connected_full_rank_base_preserves_causal_structure() {
        let parsed = ParsedGraph {
            internal_edges: vec![
                ParsedGraphInternalEdge {
                    edge_id: 0,
                    tail: 0,
                    head: 1,
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
                    tail: 0,
                    head: 1,
                    label: "q1".to_string(),
                    mass_key: Some("m1".to_string()),
                    signature: MomentumSignature {
                        loop_signature: vec![0, 1],
                        external_signature: Vec::new(),
                    },
                    had_pow: false,
                },
                ParsedGraphInternalEdge {
                    edge_id: 2,
                    tail: 0,
                    head: 1,
                    label: "q2".to_string(),
                    mass_key: Some("m2".to_string()),
                    signature: MomentumSignature {
                        loop_signature: vec![-1, -1],
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
        let standalone = generate_pure_cff_expression_from_parsed(&parsed).unwrap();
        let builder = LowerSectorCffBuilder::new(&parsed);
        let inherited = builder.build().unwrap();

        // The fast path is structurally byte-for-byte the standalone causal
        // expression, including its normalization. No parent/child frame
        // conversion is attached to it.
        assert_eq!(
            inherited.surfaces.linear_surface_cache,
            standalone.surfaces.linear_surface_cache
        );
        assert_eq!(inherited.orientations.len(), standalone.orientations.len());
        for (inherited, standalone) in inherited.orientations.iter().zip(&standalone.orientations) {
            assert_eq!(inherited.data.orientation, standalone.data.orientation);
            assert_eq!(inherited.loop_energy_map, standalone.loop_energy_map);
            assert_eq!(inherited.edge_energy_map, standalone.edge_energy_map);
            assert_eq!(inherited.variants.len(), standalone.variants.len());
            for (inherited, standalone) in inherited.variants.iter().zip(&standalone.variants) {
                assert_eq!(inherited.denominator, standalone.denominator);
                assert_eq!(inherited.numerator_surfaces, standalone.numerator_surfaces);
            }
        }
    }

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

    #[cfg(feature = "eval")]
    #[test]
    fn cut_alias_does_not_anchor_attached_tadpole_orientation() {
        let edge = |edge_id, tail, head, loop_signature, external_signature, mass: &str| {
            ParsedGraphInternalEdge {
                edge_id,
                tail,
                head,
                label: format!("q{edge_id}"),
                mass_key: Some(mass.to_string()),
                signature: MomentumSignature {
                    loop_signature,
                    external_signature,
                },
                had_pow: false,
            }
        };
        let parsed = ParsedGraph {
            internal_edges: vec![
                edge(0, 0, 1, vec![0, 0], vec![1], "m_cut"),
                edge(1, 0, 1, vec![1, 0], vec![0], "m1"),
                edge(2, 1, 0, vec![-1, 0], vec![1], "m2"),
                edge(3, 0, 0, vec![0, 1], vec![0], "m3"),
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: vec![ParsedGraphInitialStateCutEdge {
                edge_id: 0,
                external_id: 0,
                external_sign: 1,
            }],
            loop_names: vec!["k0".to_string(), "k1".to_string()],
            external_names: vec!["p0".to_string()],
            node_name_to_internal: BTreeMap::from([("n0".to_string(), 0), ("n1".to_string(), 1)]),
        };

        let mut builder = LowerSectorCffBuilder::new(&parsed);
        builder.force_component_factorization = true;
        let lower = builder.build().unwrap();
        assert!(!lower.orientations.is_empty());
        let public = generate_pure_cff_expression_from_parsed(&parsed).unwrap();
        for seed in [17, 41] {
            let input = crate::eval::EvaluationInput::deterministic(
                &parsed,
                seed,
                &BTreeMap::from([
                    ("m0".to_string(), 0.17),
                    ("m1".to_string(), 0.23),
                    ("m2".to_string(), 0.31),
                    ("m3".to_string(), 0.37),
                    ("m_cut".to_string(), 1.0),
                ]),
                None,
            )
            .unwrap();
            let public_value = crate::eval::evaluate_expression(&parsed, &public, "1", &input)
                .unwrap()
                .value;
            let factorized_value = crate::eval::evaluate_expression(&parsed, &lower, "1", &input)
                .unwrap()
                .value;
            let scale = public_value
                .abs()
                .max(factorized_value.abs())
                .max(f64::MIN_POSITIVE);
            assert!(
                (public_value - factorized_value).abs() <= 1.0e-13 * scale,
                "a cut touching only the attached component's incidence vertex changed its independent public CFF sum at seed {seed}: public={public_value:e}, factorized={factorized_value:e}"
            );
        }
        assert!(lower.orientations.iter().any(|orientation| {
            orientation.edge_energy_map[3] == LinearEnergyExpr::ose(EdgeIndex(3), 1)
        }));
        assert!(lower.orientations.iter().any(|orientation| {
            orientation.edge_energy_map[3] == LinearEnergyExpr::ose(EdgeIndex(3), -1)
        }));
    }

    #[cfg(feature = "eval")]
    #[test]
    fn cut_aware_vector_matroid_projection_anchors_only_internal_cut() {
        let edge = |edge_id, tail, head, loop_signature, external_signature, mass: &str| {
            ParsedGraphInternalEdge {
                edge_id,
                tail,
                head,
                label: format!("q{edge_id}"),
                mass_key: Some(mass.to_string()),
                signature: MomentumSignature {
                    loop_signature,
                    external_signature,
                },
                had_pow: false,
            }
        };
        // Incidence sees two bubbles joined at vertex 0, whereas the rational
        // energy rows split into independent k0 and k1 factors. The cut remains
        // internal to the first factor after vector-matroid projection.
        let parsed = ParsedGraph {
            internal_edges: vec![
                edge(0, 0, 1, vec![1, 0], vec![0], "m0"),
                edge(1, 1, 0, vec![1, 0], vec![1], "m1"),
                edge(2, 0, 2, vec![0, 1], vec![0], "m2"),
                edge(3, 2, 0, vec![0, 1], vec![0], "m3"),
                edge(4, 0, 1, vec![0, 0], vec![1], "m_cut"),
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: vec![ParsedGraphInitialStateCutEdge {
                edge_id: 4,
                external_id: 0,
                external_sign: 1,
            }],
            loop_names: vec!["k0".to_string(), "k1".to_string()],
            external_names: vec!["p0".to_string()],
            node_name_to_internal: (0..3).map(|node| (format!("n{node}"), node)).collect(),
        };
        let mut builder = LowerSectorCffBuilder::new(&parsed);
        builder.force_component_factorization = true;
        let signatures = builder.signatures();
        let components = builder.vector_matroid_components(&signatures);
        assert_eq!(components, vec![vec![0, 1], vec![2, 3]]);
        let (cut_component, cut_embedding) = builder
            .project_component_parsed(&signatures, &components[0], &[0])
            .unwrap();
        assert_eq!(cut_embedding, vec![0, 1, 4]);
        assert_eq!(
            cut_component.initial_state_cut_edges,
            vec![ParsedGraphInitialStateCutEdge {
                edge_id: 2,
                external_id: 0,
                external_sign: 1,
            }]
        );
        let (cut_disjoint_component, cut_disjoint_embedding) = builder
            .project_component_parsed(&signatures, &components[1], &[2])
            .unwrap();
        assert_eq!(cut_disjoint_embedding, vec![2, 3]);
        assert!(cut_disjoint_component.initial_state_cut_edges.is_empty());
        assert!(cut_disjoint_component.external_edges.is_empty());
        let lower = builder.build().unwrap();
        assert!(!lower.orientations.is_empty());
    }

    #[cfg(feature = "eval")]
    #[test]
    fn shared_cut_alias_is_retained_in_each_rational_factor() {
        let edge = |edge_id, tail, head, loop_signature, external_signature, mass: &str| {
            ParsedGraphInternalEdge {
                edge_id,
                tail,
                head,
                label: format!("q{edge_id}"),
                mass_key: Some(mass.to_string()),
                signature: MomentumSignature {
                    loop_signature,
                    external_signature,
                },
                had_pow: false,
            }
        };
        let parsed = ParsedGraph {
            internal_edges: vec![
                edge(0, 0, 1, vec![1, 0], vec![0], "m0"),
                edge(1, 1, 0, vec![1, 0], vec![1], "m1"),
                edge(2, 0, 1, vec![0, 1], vec![0], "m2"),
                edge(3, 1, 0, vec![0, 1], vec![1], "m3"),
                edge(4, 0, 1, vec![0, 0], vec![1], "m_cut"),
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: vec![ParsedGraphInitialStateCutEdge {
                edge_id: 4,
                external_id: 0,
                external_sign: 1,
            }],
            loop_names: vec!["k0".to_string(), "k1".to_string()],
            external_names: vec!["p0".to_string()],
            node_name_to_internal: BTreeMap::from([("n0".to_string(), 0), ("n1".to_string(), 1)]),
        };
        let mut builder = LowerSectorCffBuilder::new(&parsed);
        builder.force_component_factorization = true;
        let signatures = builder.signatures();
        let components = builder.vector_matroid_components(&signatures);
        assert_eq!(components, vec![vec![0, 1], vec![2, 3]]);
        for (component, basis) in components.iter().zip([0, 2]) {
            let (projected, embedding) = builder
                .project_component_parsed(&signatures, component, &[basis])
                .unwrap();
            assert_eq!(embedding.last(), Some(&4));
            assert_eq!(projected.initial_state_cut_edges.len(), 1);
            assert_eq!(
                projected.denominator_internal_edge_ids().len(),
                component.len()
            );
        }

        let lower = builder.build().unwrap();
        assert!(lower.orientations.iter().all(|orientation| {
            orientation.edge_energy_map[4] == LinearEnergyExpr::external(EdgeIndex(0), 1)
                && orientation.variants.iter().all(|variant| {
                    !variant.half_edges.contains(&EdgeIndex(4))
                        && !variant.denominator_edges.contains(&EdgeIndex(4))
                })
        }));
    }

    #[cfg(feature = "eval")]
    #[test]
    fn attached_repeated_factor_selector_branches_match_lower_sector_cff() {
        let edge = |edge_id, tail, head, loop_signature, external_signature, mass: &str| {
            ParsedGraphInternalEdge {
                edge_id,
                tail,
                head,
                label: format!("q{edge_id}"),
                mass_key: Some(mass.to_string()),
                signature: MomentumSignature {
                    loop_signature,
                    external_signature: vec![external_signature],
                },
                had_pow: false,
            }
        };
        let parsed = ParsedGraph {
            internal_edges: vec![
                edge(0, 0, 3, vec![1, 0, 0], 0, "m0"),
                edge(1, 4, 0, vec![1, 0, 0], 1, "m0"),
                edge(2, 1, 2, vec![0, 1, 0], 0, "m0"),
                edge(3, 3, 1, vec![0, 1, 0], -1, "m0"),
                edge(4, 3, 4, vec![1, -1, 0], 1, "m0"),
                edge(5, 2, 4, vec![0, 1, 0], 0, "m1"),
                edge(6, 0, 1, vec![0, 0, 0], 1, "m_cut"),
                edge(7, 0, 5, vec![0, 0, 1], 0, "m_uv"),
                edge(8, 5, 0, vec![0, 0, 1], 0, "m_uv"),
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: vec![ParsedGraphInitialStateCutEdge {
                edge_id: 6,
                external_id: 0,
                external_sign: 1,
            }],
            loop_names: vec!["q1".to_string(), "q3".to_string(), "q_uv".to_string()],
            external_names: vec!["p0".to_string()],
            node_name_to_internal: (0..6).map(|node| (format!("n{node}"), node)).collect(),
        };
        let lower_sector = LowerSectorCffBuilder::new(&parsed).build().unwrap();
        let input = crate::eval::EvaluationInput::deterministic(
            &parsed,
            17,
            &BTreeMap::from([
                ("m0".to_string(), 0.0),
                ("m1".to_string(), 0.1),
                ("m_cut".to_string(), 1.0),
                ("m_uv".to_string(), 20.0),
            ]),
            None,
        )
        .unwrap();
        let edge_zero_energy = input.loop_spatial_momenta[0]
            .iter()
            .map(|component| component * component)
            .sum::<f64>()
            .sqrt();
        let evaluate_lower = |build: &ThreeDExpression<OrientationID>, numerator: &str| {
            crate::eval::evaluate_expression(&parsed, build, numerator, &input)
                .unwrap()
                .value
        };
        let standalone_scalar = evaluate_lower(&lower_sector, "1");
        let standalone_linear = evaluate_lower(&lower_sector, "edges[0][0]");
        for sample in [1, -1] {
            let mut branch = KnownFactorCffBuilder::new(
                &parsed,
                vec![0; parsed.internal_edges.len()],
                NumeratorSamplingScaleMode::None,
            );
            let replacement = LinearEnergyExpr::ose(EdgeIndex(0), i64::from(sample));
            let selector = if sample == 1 {
                KnownLinearExpr::var(0, 1) + KnownLinearExpr::ose(0, 1)
            } else {
                KnownLinearExpr::ose(0, 1) - KnownLinearExpr::var(0, 1)
            };
            branch
                .append_base_terms(
                    &parsed,
                    &(0..parsed.internal_edges.len()).collect::<Vec<_>>(),
                    &BTreeMap::from([(0, replacement)]),
                    &[selector],
                    &[0],
                    &[],
                    0,
                    Rational::one(),
                    false,
                    false,
                )
                .unwrap();
            let actual = crate::eval::evaluate_expression(
                &parsed,
                &branch.expression,
                "edges[0][0]**4",
                &input,
            )
            .unwrap()
            .value;
            let selected = if sample == 1 {
                standalone_linear + edge_zero_energy * standalone_scalar
            } else {
                edge_zero_energy * standalone_scalar - standalone_linear
            } / (2.0 * edge_zero_energy);
            let expected = edge_zero_energy.powi(4) * selected;
            let scale = actual.abs().max(expected.abs()).max(f64::MIN_POSITIVE);
            assert!(
                (actual - expected).abs() <= 1.0e-12 * scale,
                "factorized append-base selector branch {sample:+} differs from its direct affine LowerSector action: actual={actual:e}, expected={expected:e}"
            );
        }
    }
    #[cfg(feature = "eval")]
    #[test]
    fn gl24_lower_sector_cff_is_invariant_under_edge_reversal() {
        let edge =
            |edge_id, tail, head, loop_signature, external_signature| ParsedGraphInternalEdge {
                edge_id,
                tail,
                head,
                label: format!("q{edge_id}"),
                mass_key: Some(format!("m{edge_id}")),
                signature: MomentumSignature {
                    loop_signature,
                    external_signature: vec![external_signature],
                },
                had_pow: false,
            };
        // This is the smallest denominator-deleting GL24 quotient: the
        // q3/q5 component carries the selected LU cut and q1-q3+p is its
        // independent one-denominator component.
        let parsed = ParsedGraph {
            internal_edges: vec![
                edge(0, 0, 1, vec![0, 0, 0], 1),
                edge(1, 1, 2, vec![0, 1, 0], 0),
                edge(2, 1, 0, vec![0, -1, 0], 1),
                edge(3, 2, 0, vec![0, 0, 1], 0),
                edge(4, 2, 0, vec![0, 1, -1], 0),
                edge(5, 0, 0, vec![1, -1, 0], 1),
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: vec![ParsedGraphInitialStateCutEdge {
                edge_id: 0,
                external_id: 0,
                external_sign: 1,
            }],
            loop_names: vec!["q1".to_string(), "q3".to_string(), "q5".to_string()],
            external_names: vec!["p0".to_string()],
            node_name_to_internal: (0..3).map(|node| (format!("v{node}"), node)).collect(),
        };
        let input = crate::eval::EvaluationInput {
            external_momenta: vec![[1.37, 0.0, 0.0, 0.0]],
            loop_spatial_momenta: vec![
                [0.13, -0.09, 0.21],
                [-0.14160937677509397, -0.43582884761035434, 0.2],
                [-0.16, 0.19, 0.08],
            ],
            masses: (0..parsed.internal_edges.len())
                .map(|edge_id| 0.04 * edge_id as f64)
                .collect(),
            uniform_scale: None,
        };

        let build_and_evaluate = |label: &str, graph: &ParsedGraph| {
            let signatures = graph
                .internal_edges
                .iter()
                .map(|edge| edge.signature.clone())
                .collect::<Vec<_>>();
            let builder = LowerSectorCffBuilder::new(graph);
            assert_eq!(
                builder.vector_matroid_components(&signatures),
                vec![vec![1, 2, 3, 4], vec![5]],
            );
            let expression = builder.build().unwrap();
            assert!(
                !expression.orientations.is_empty(),
                "{label} GL24 quotient produced no CFF orientations"
            );
            crate::eval::evaluate_expression(graph, &expression, "1", &input)
                .unwrap()
                .value
        };

        let positive = build_and_evaluate("Q", &parsed);
        let mut reversed = parsed.clone();
        let isolated = &mut reversed.internal_edges[5];
        isolated
            .signature
            .loop_signature
            .iter_mut()
            .for_each(|value| *value *= -1);
        isolated
            .signature
            .external_signature
            .iter_mut()
            .for_each(|value| *value *= -1);

        // D(Q)=D(-Q). Component projection must retain the source
        // loop-coordinate reversal while the complete CFF remains invariant.
        let signatures = reversed
            .internal_edges
            .iter()
            .map(|edge| edge.signature.clone())
            .collect::<Vec<_>>();
        let builder = LowerSectorCffBuilder::new(&reversed);
        let (isolated, local_to_parent) = builder
            .project_component_parsed(&signatures, &[5], &[5])
            .unwrap();
        assert_eq!(local_to_parent, vec![5]);
        assert_eq!(
            isolated.internal_edges[0].signature.loop_signature,
            vec![-1],
        );
        let negative = build_and_evaluate("-Q", &reversed);
        assert!((positive - negative).abs() <= 1.0e-12 * positive.abs());
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
    fn finite_pole_sampling_updates_entire_repeated_logical_channel() {
        let mut parsed = crate::graph_io::test_graphs::box_pow3_graph();
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
        let builder = BoundedCffBuilder::new(
            &parsed,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(3, 2)]),
                ..Default::default()
            },
        )
        .unwrap();
        let mut edge_maps = (0..parsed.internal_edges.len())
            .map(|edge_id| LinearEnergyExpr::ose(EdgeIndex(edge_id), 1))
            .collect::<Vec<_>>();
        edge_maps[5] = LinearEnergyExpr::ose(EdgeIndex(5), -1);

        builder
            .assign_repeated_channel_sample(&mut edge_maps, 3, 0)
            .unwrap();

        assert!(edge_maps[3].is_zero());
        assert!(edge_maps[4].is_zero());
        assert!(edge_maps[5].is_zero());

        builder
            .assign_repeated_channel_sample(&mut edge_maps, 3, -1)
            .unwrap();

        assert_eq!(edge_maps[3], LinearEnergyExpr::ose(EdgeIndex(3), -1));
        assert_eq!(edge_maps[4], LinearEnergyExpr::ose(EdgeIndex(4), -1));
        assert_eq!(edge_maps[5], LinearEnergyExpr::ose(EdgeIndex(5), 1));
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

    #[cfg(feature = "eval")]
    #[test]
    fn structural_cut_cff_is_invariant_under_conservative_inner_bound() {
        let edge =
            |edge_id, tail, head, loop_signature, external_signature| ParsedGraphInternalEdge {
                edge_id,
                tail,
                head,
                label: format!("e{edge_id}"),
                mass_key: Some("m".to_string()),
                signature: MomentumSignature {
                    loop_signature,
                    external_signature,
                },
                had_pow: false,
            };
        let parsed = ParsedGraph {
            internal_edges: vec![
                edge(0, 2, 3, vec![0, 0], vec![1]),
                edge(1, 0, 1, vec![1, 0], vec![0]),
                edge(2, 0, 1, vec![0, 1], vec![0]),
                edge(3, 3, 0, vec![1, 1], vec![0]),
                edge(4, 1, 2, vec![1, 1], vec![0]),
                edge(5, 2, 3, vec![1, 1], vec![-1]),
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: vec![ParsedGraphInitialStateCutEdge {
                edge_id: 0,
                external_id: 0,
                external_sign: 1,
            }],
            loop_names: vec!["q1".to_string(), "q2".to_string()],
            external_names: vec!["p0".to_string()],
            node_name_to_internal: (0..4).map(|node| (format!("n{node}"), node)).collect(),
        };
        let generate = |bounds| {
            generate_3d_expression(
                &parsed,
                &Generate3DExpressionOptions {
                    energy_degree_bounds: Some(bounds),
                    ..Default::default()
                },
            )
            .unwrap()
        };
        let exact = generate(vec![(1, 1), (3, 2), (4, 1)]);
        let conservative = generate(vec![(1, 2), (3, 2), (4, 1)]);
        assert_eq!(
            exact.core_global_prefactor_sign,
            conservative.core_global_prefactor_sign
        );
        for (label, generated) in [("exact", &exact), ("conservative", &conservative)] {
            let origins = generated
                .expression
                .orientations
                .iter()
                .flat_map(|orientation| &orientation.variants)
                .filter_map(|variant| variant.origin.as_deref())
                .collect::<Vec<_>>();
            assert!(
                origins.contains(&"bounded_degree_known_factor_cff"),
                "{label} structural-cut CFF did not retain its factorized lower-sector branch: {origins:?}",
            );
            // A conservative quadratic capacity may add internal interpolation
            // contacts. Those variants still sample the same runtime
            // factorized numerator; their presence is not numerator expansion.
        }

        let numerator = "(edges[1][0]*edges[3][0]-edges[1][1]*edges[3][1]-edges[1][2]*edges[3][2]-edges[1][3]*edges[3][3])*edges[3][0]*edges[4][0]";
        for seed in [17, 41, 73] {
            let input = crate::eval::EvaluationInput::deterministic(
                &parsed,
                seed,
                &BTreeMap::from([("m".to_string(), 0.71)]),
                None,
            )
            .unwrap();
            let evaluate = |generated: &GeneratedThreeDExpression| {
                crate::eval::evaluate_expression(&parsed, &generated.expression, numerator, &input)
                    .unwrap()
                    .value
                    * generated.core_global_prefactor_sign.factor() as f64
            };
            let exact_value = evaluate(&exact);
            let conservative_value = evaluate(&conservative);
            let scale = exact_value
                .abs()
                .max(conservative_value.abs())
                .max(f64::MIN_POSITIVE);
            assert!(
                (exact_value - conservative_value).abs() <= 1.0e-11 * scale,
                "conservative inner bound changes the full unselected CFF value at seed {seed}: exact={exact_value:.17e}, conservative={conservative_value:.17e}",
            );
        }
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
                denominator_only_global_prefactor_sign: generated
                    .denominator_only_global_prefactor_sign,
                core_global_prefactor_sign: generated.core_global_prefactor_sign,
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
    fn repeated_p6_and_p5_share_one_variant_local_energy_frame() {
        let repeated_cycle = |power: usize| ParsedGraph {
            internal_edges: (0..power)
                .map(|edge_id| ParsedGraphInternalEdge {
                    edge_id,
                    tail: edge_id,
                    head: (edge_id + 1) % power,
                    label: format!("q{edge_id}"),
                    mass_key: Some("m".to_string()),
                    signature: MomentumSignature {
                        loop_signature: vec![-1],
                        external_signature: Vec::new(),
                    },
                    had_pow: false,
                })
                .collect(),
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["ell".to_string()],
            external_names: Vec::new(),
            node_name_to_internal: (0..power).map(|node| (format!("v{node}"), node)).collect(),
        };
        let generated = [
            (
                repeated_cycle(6),
                vec![(1, 2), (2, 1), (3, 1), (4, 1), (5, 1)],
            ),
            (repeated_cycle(5), vec![(1, 1), (2, 1), (3, 1), (4, 1)]),
        ]
        .map(|(parsed, bounds)| {
            generate_3d_expression_from_parsed_generated(
                &parsed,
                &Generate3DExpressionOptions {
                    energy_degree_bounds: Some(bounds),
                    ..Default::default()
                },
            )
            .unwrap()
        });

        for expression in generated {
            assert_eq!(
                expression.energy_factor_ownership,
                CffEnergyFactorOwnership::VariantLocal
            );
            assert_eq!(expression.energy_factor_components.len(), 1);
            assert_eq!(expression.core_global_prefactor_sign.factor(), 1);
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
                denominator_only_global_prefactor_sign: generated
                    .denominator_only_global_prefactor_sign,
                core_global_prefactor_sign: generated.core_global_prefactor_sign,
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
                denominator_only_global_prefactor_sign: generated
                    .denominator_only_global_prefactor_sign,
                core_global_prefactor_sign: generated.core_global_prefactor_sign,
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
    fn disconnected_cograph_projection_retains_its_initial_state_cut() {
        let edge =
            |edge_id, tail, head, loop_signature, external_signature| ParsedGraphInternalEdge {
                edge_id,
                tail,
                head,
                label: format!("q{edge_id}"),
                mass_key: Some("m".to_string()),
                signature: MomentumSignature {
                    loop_signature,
                    external_signature,
                },
                had_pow: false,
            };
        let parsed = ParsedGraph {
            internal_edges: vec![
                // The cut alias belongs to the cograph. Its stored second-loop
                // row deliberately cannot add a contour integration variable.
                edge(0, 0, 1, vec![0, 1], vec![1, 0, 0]),
                edge(1, 0, 2, vec![1, 0], vec![0, 0, 0]),
                edge(2, 2, 1, vec![-1, 0], vec![1, 0, 0]),
                // A genuinely disconnected repeated vacuum bubble.
                edge(3, 3, 4, vec![0, 1], vec![0, 0, 0]),
                edge(4, 4, 3, vec![0, -1], vec![0, 0, 0]),
                // This cut crosses the two denominator components. It remains
                // external boundary flow and cannot merge their loop spaces.
                edge(5, 1, 3, vec![1, 1], vec![0, 1, 0]),
                // A cut-only structural component has no denominator contour
                // and must not be emitted as a standalone CFF factor.
                edge(6, 5, 6, vec![1, 0], vec![0, 0, 1]),
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: vec![
                ParsedGraphInitialStateCutEdge {
                    edge_id: 0,
                    external_id: 0,
                    external_sign: 1,
                },
                ParsedGraphInitialStateCutEdge {
                    edge_id: 5,
                    external_id: 1,
                    external_sign: 1,
                },
                ParsedGraphInitialStateCutEdge {
                    edge_id: 6,
                    external_id: 2,
                    external_sign: 1,
                },
            ],
            loop_names: vec!["kc".to_string(), "kuv".to_string()],
            external_names: vec!["p0".to_string(), "p1".to_string(), "p2".to_string()],
            node_name_to_internal: (0..=6).map(|node| (format!("n{node}"), node)).collect(),
        };

        let components = denominator_connected_components(&parsed);
        assert_eq!(components, vec![vec![0, 1, 2], vec![3, 4]]);

        let (cograph, cograph_embedding) =
            project_denominator_component(&parsed, &components[0]).unwrap();
        assert_eq!(cograph_embedding.local_to_orig_edge, vec![0, 1, 2]);
        assert_eq!(cograph_embedding.local_to_orig_loop, vec![0]);
        assert_eq!(
            cograph.initial_state_cut_edges,
            vec![ParsedGraphInitialStateCutEdge {
                edge_id: 0,
                external_id: 0,
                external_sign: 1,
            }]
        );
        assert_eq!(cograph.external_edges.len(), 1);
        assert_eq!(cograph.external_edges[0].label, "q5_source_boundary");

        let (vacuum_bubble, vacuum_embedding) =
            project_denominator_component(&parsed, &components[1]).unwrap();
        assert_eq!(vacuum_embedding.local_to_orig_edge, vec![3, 4]);
        assert_eq!(vacuum_embedding.local_to_orig_loop, vec![1]);
        assert!(vacuum_bubble.initial_state_cut_edges.is_empty());
        assert_eq!(vacuum_bubble.external_edges.len(), 1);
        assert_eq!(vacuum_bubble.external_edges[0].label, "q5_sink_boundary");

        let generated = generate_3d_expression_from_parsed_generated(
            &parsed,
            &Generate3DExpressionOptions::default(),
        )
        .unwrap();
        assert!(generated.expression.orientations.iter().all(|orientation| {
            orientation.edge_energy_map[0] == LinearEnergyExpr::external(EdgeIndex(0), 1)
                && orientation.edge_energy_map[5] == LinearEnergyExpr::external(EdgeIndex(1), 1)
                && orientation.edge_energy_map[6] == LinearEnergyExpr::external(EdgeIndex(2), 1)
        }));
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
                    denominator_only_global_prefactor_sign: CffGlobalPrefactorSign::from_exponent(
                        1
                    ),
                    core_global_prefactor_sign: CffGlobalPrefactorSign::from_exponent(0),
                },
                CffEnergyFactorComponent {
                    internal_edge_ids: vec![2, 3],
                    ownership: CffEnergyFactorOwnership::GlobalSourceProduct,
                    denominator_only_global_prefactor_sign: CffGlobalPrefactorSign::from_exponent(
                        0,
                    ),
                    core_global_prefactor_sign: CffGlobalPrefactorSign::from_exponent(0),
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

    #[cfg(feature = "eval")]
    #[test]
    fn connected_odd_triangle_bubble_keeps_rational_component_product() {
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
                        external_signature: vec![external_signature],
                    },
                    had_pow: false,
                }
            };
        // Incidence joins the odd triangle and bubble at node zero. The
        // bounded source therefore remains one factorized numerator evaluated
        // in one variant-local rational frame.
        let parsed = ParsedGraph {
            internal_edges: vec![
                edge(0, 0, 1, [1, 0], 0, "m0"),
                edge(1, 1, 2, [1, 0], 1, "m1"),
                edge(2, 2, 0, [1, 0], -1, "m2"),
                edge(3, 0, 3, [0, 1], 0, "m3"),
                edge(4, 3, 0, [0, 1], 1, "m4"),
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["x".to_string(), "y".to_string()],
            external_names: vec!["p".to_string()],
            node_name_to_internal: (0..4).map(|node| (format!("v{node}"), node)).collect(),
        };
        let quadratic = generate_3d_expression_from_parsed_generated(
            &parsed,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(0, 2)]),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            quadratic.energy_factor_components,
            vec![CffEnergyFactorComponent {
                internal_edge_ids: vec![0, 1, 2, 3, 4],
                ownership: CffEnergyFactorOwnership::VariantLocal,
                denominator_only_global_prefactor_sign: quadratic
                    .denominator_only_global_prefactor_sign,
                core_global_prefactor_sign: quadratic.core_global_prefactor_sign,
            }]
        );
    }

    #[test]
    fn connected_incidence_does_not_use_disconnected_component_product() {
        let edge = |edge_id, tail, head, loop_signature: [i32; 2], external_signature| {
            ParsedGraphInternalEdge {
                edge_id,
                tail,
                head,
                label: format!("q{edge_id}"),
                mass_key: Some(format!("m{edge_id}")),
                signature: MomentumSignature {
                    loop_signature: loop_signature.to_vec(),
                    external_signature: vec![external_signature],
                },
                had_pow: false,
            }
        };
        let unique = ParsedGraph {
            internal_edges: vec![
                edge(0, 0, 1, [1, 0], 0),
                edge(1, 1, 0, [1, 0], 1),
                edge(2, 0, 2, [0, 1], 0),
                edge(3, 2, 0, [0, 1], 1),
                edge(4, 0, 1, [0, 0], 1),
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: vec![ParsedGraphInitialStateCutEdge {
                edge_id: 4,
                external_id: 0,
                external_sign: 1,
            }],
            loop_names: vec!["x".to_string(), "y".to_string()],
            external_names: vec!["p".to_string()],
            node_name_to_internal: (0..3).map(|node| (format!("n{node}"), node)).collect(),
        };
        let options = Generate3DExpressionOptions {
            energy_degree_bounds: Some(vec![(0, 2)]),
            ..Default::default()
        };
        assert!(
            generate_rational_component_product(&unique, &options)
                .unwrap()
                .is_none(),
            "a connected denominator incidence must stay on the whole-graph construction"
        );

        let mut ambiguous = unique;
        ambiguous.internal_edges[2].head = 1;
        ambiguous.internal_edges[3].tail = 1;
        assert!(
            generate_rational_component_product(&ambiguous, &options)
                .unwrap()
                .is_none(),
            "changing the cut ownership cannot turn connected denominator incidence into a disconnected product"
        );
    }

    #[test]
    fn connected_generalized_product_tracks_lifted_global_duplicates() {
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
                        external_signature: vec![external_signature],
                    },
                    had_pow: false,
                }
            };
        // Incidence joins the two factors at node zero, but their loop rows
        // split into k and q components. The k bubble has two distinct affine
        // denominators, so it carries no standalone duplicate-line parity; the
        // q component alone needs generalized numerator sampling.
        let parsed = ParsedGraph {
            internal_edges: vec![
                edge(0, 0, 1, [1, 0], 0, "mk"),
                edge(1, 1, 0, [1, 0], 1, "mk"),
                edge(2, 0, 2, [0, 1], 0, "mq"),
                edge(3, 2, 3, [0, 1], 0, "mq"),
                edge(4, 3, 0, [0, 1], 0, "mq"),
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["k".to_string(), "q".to_string()],
            external_names: vec!["p".to_string()],
            node_name_to_internal: (0..4).map(|node| (format!("n{node}"), node)).collect(),
        };
        let signatures = parsed
            .internal_edges
            .iter()
            .map(|edge| edge.signature.clone())
            .collect::<Vec<_>>();
        let builder = LowerSectorCffBuilder::new(&parsed);
        let components = builder.vector_matroid_components(&signatures);
        assert_eq!(components, vec![vec![0, 1], vec![2, 3, 4]]);
        let (ordinary_component, _) = builder
            .project_component_parsed(&signatures, &components[0], &[0])
            .unwrap();
        assert_eq!(cff_duplicate_signature_excess(&ordinary_component), 0);

        let generated = generate_3d_expression_from_parsed_generated(
            &parsed,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(2, 2)]),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            generated
                .energy_factor_components
                .iter()
                .map(|component| (component.internal_edge_ids.len(), component.ownership))
                .collect::<Vec<_>>(),
            vec![(5, CffEnergyFactorOwnership::VariantLocal)],
        );
        assert_eq!(
            generated.core_global_prefactor_sign.factor(),
            -1,
            "without a component-local duplicate parity, the connected two-loop product keeps only its established (-1)^(L-1) bridge",
        );

        let mut parsed_with_pure_duplicate = parsed;
        parsed_with_pure_duplicate.internal_edges[1]
            .signature
            .external_signature = vec![0];
        let signatures = parsed_with_pure_duplicate
            .internal_edges
            .iter()
            .map(|edge| edge.signature.clone())
            .collect::<Vec<_>>();
        let builder = LowerSectorCffBuilder::new(&parsed_with_pure_duplicate);
        let components = builder.vector_matroid_components(&signatures);
        let (ordinary_component, _) = builder
            .project_component_parsed(&signatures, &components[0], &[0])
            .unwrap();
        assert_eq!(cff_duplicate_signature_excess(&ordinary_component), 1);

        let generated = generate_3d_expression_from_parsed_generated(
            &parsed_with_pure_duplicate,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(2, 2)]),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            generated
                .energy_factor_components
                .iter()
                .map(|component| (component.internal_edge_ids.len(), component.ownership))
                .collect::<Vec<_>>(),
            vec![(5, CffEnergyFactorOwnership::VariantLocal)],
        );
        assert_eq!(
            generated.core_global_prefactor_sign.factor(),
            -1,
            "a connected whole-graph construction keeps its established (-1)^(L-1) bridge even when one rational subspace contains duplicate signatures",
        );
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
    fn denominator_only_prefactor_is_bound_invariant_across_duplicate_parity() {
        let repeated_cycle = |power: usize| ParsedGraph {
            internal_edges: (0..power)
                .map(|edge_id| ParsedGraphInternalEdge {
                    edge_id,
                    tail: edge_id,
                    head: (edge_id + 1) % power,
                    label: format!("q{edge_id}"),
                    mass_key: Some("m".to_string()),
                    signature: MomentumSignature {
                        loop_signature: vec![1],
                        external_signature: Vec::new(),
                    },
                    had_pow: false,
                })
                .collect(),
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["q".to_string()],
            external_names: Vec::new(),
            node_name_to_internal: (0..power).map(|node| (format!("n{node}"), node)).collect(),
        };

        for (power, expected_denominator_sign) in [(2, -1), (3, 1)] {
            let parsed = repeated_cycle(power);
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
                    energy_degree_bounds: Some(vec![(0, 2)]),
                    ..Default::default()
                },
            )
            .unwrap();

            assert_eq!(
                pure.denominator_only_global_prefactor_sign.factor(),
                expected_denominator_sign
            );
            assert_eq!(
                generalized.denominator_only_global_prefactor_sign.factor(),
                expected_denominator_sign,
                "the physical denominator bridge changed when only the numerator bound changed for power {power}",
            );
            assert_eq!(
                generalized.core_global_prefactor_sign.factor(),
                1,
                "the generalized Laurent-functional frame remains independent metadata",
            );
        }
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
                cff_generation_context: CffGenerationContext::Standalone,
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
                cff_generation_context: CffGenerationContext::Standalone,
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
                cff_generation_context: CffGenerationContext::Standalone,
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
                cff_generation_context: CffGenerationContext::Standalone,
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
                cff_generation_context: CffGenerationContext::Standalone,
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
