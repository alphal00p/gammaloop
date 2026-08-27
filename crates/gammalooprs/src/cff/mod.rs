use std::{collections::BTreeMap, fmt::Display, sync::Arc};

use bincode_trait_derive::{Decode, Encode};
use linnet::half_edge::{
    involution::EdgeIndex,
    subgraph::{SubGraphLike, SubSetLike, SubSetOps},
};
use serde::{Deserialize, Serialize};
use symbolica::atom::{Atom, AtomCore};

use crate::{
    cff::{
        expression::{
            GammaLoopOrientationExpression, OrientationExpression, OrientationID,
            OrientationSelector, ThreeDExpression,
            normalize_three_d_expression_cut_support_with_raised_edge_groups,
        },
        orientations::GraphOrientation,
        surface::GammaLoopSurfaceCache,
    },
    graph::{
        FeynmanGraph, FourDDenominator, Graph, GraphThreeDSource, cuts::CutSet,
        get_cff_inverse_energy_product_impl, three_d_source::ExactSourceEnergyMapper,
    },
    numerator::energy_degree::EnergyPowerAssignmentPlan,
    settings::global::OrientationPattern,
    utils::GS,
    uv::Integrands,
};
use color_eyre::Result;
use three_dimensional_reps::{
    CffEnergyFactorOwnership, CffGlobalPrefactorSign, Generate3DExpressionOptions,
    GeneratedThreeDExpression,
};

pub mod orientations;
//pub mod cut_expression;
pub mod esurface;
pub mod expression;
pub mod generation;
pub mod hsurface;
pub mod surface;
pub mod tree;
mod vertex_set;
pub(crate) use vertex_set::VertexSet;

pub(crate) struct CFFOrientationTerm {
    pub(crate) expression: Atom,
    pub(crate) orientation: OrientationExpression,
    pub(crate) production_orientation_id: Option<OrientationID>,
}

pub struct CFFTerm {
    // Ordinary CFF maps retain production identity for direct-3D localization;
    // exact 4D maps instead remain source-local and share their parent mapper.
    pub(crate) orientations: Vec<CFFOrientationTerm>,
    exact_source_numerator: Option<Arc<PlannedExactSourceNumerator>>,
}

struct PlannedExactSourceNumerator {
    mapper: ExactSourceEnergyMapper,
    assignment: EnergyPowerAssignmentPlan,
}

impl CFFTerm {
    pub(crate) fn map_exact_source_numerator(
        &self,
        orientation: &OrientationExpression,
    ) -> Result<Atom> {
        let planned = self
            .exact_source_numerator
            .as_ref()
            .ok_or_else(|| eyre::eyre!("ordinary CFF term has no exact-source numerator plan"))?;
        planned.mapper.map_planned_numerator(
            &orientation.loop_energy_map,
            &orientation.edge_energy_map,
            &planned.assignment,
        )
    }

    pub fn expression_with_selectors(&self) -> Atom {
        self.orientations
            .iter()
            .map(|term| {
                term.expression.clone() * term.orientation.data.orientation.orientation_thetas()
            })
            .reduce(|left, right| left + right)
            .unwrap_or(Atom::Zero)
    }
}

#[derive(
    Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Hash, Encode, Decode, Serialize, Deserialize,
)]
// This describes the combinations of residues that are selected.
pub struct CutCFFIndex {
    pub left_threshold_order: Option<usize>,
    pub right_threshold_order: Option<usize>,
    pub lu_cut_order: Option<usize>,
}

impl CutCFFIndex {
    pub fn new_all_none() -> Self {
        Self {
            left_threshold_order: None,
            right_threshold_order: None,
            lu_cut_order: None,
        }
    }
}

impl Display for CutCFFIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut parts = vec![];
        if let Some(order) = self.lu_cut_order {
            parts.push(format!("lu_cut_{}", order));
        }

        if let Some(order) = self.left_threshold_order {
            parts.push(format!("left_th_{}", order));
        }

        if let Some(order) = self.right_threshold_order {
            parts.push(format!("right_th_{}", order));
        }

        if parts.is_empty() {
            write!(f, "")
        } else {
            write!(f, "{}", parts.join("_"))
        }
    }
}

/// Namespace used by the energy-degree bounds supplied to a CFF source.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum CffEnergyBoundSourceKind {
    /// Source IDs are physical graph-edge IDs.
    PhysicalGraph,
    /// Source IDs are exact four-dimensional denominator-occurrence IDs.
    ExactFourD,
}

/// Transient diagnostic record for one CFF generation source.
///
/// Exact four-dimensional Taylor terms can contain several denominator
/// occurrences owned by one physical edge. The CFF generator must receive the
/// minimax assignment in that occurrence namespace, while the physical-parent
/// bounds remain useful for checking the numerator analysis. Keeping both
/// explicitly prevents occurrence IDs from being mistaken for physical edges.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct CffEnergyDegreeBoundReport {
    pub source_kind: CffEnergyBoundSourceKind,
    pub physical_parent_bounds: Vec<(usize, usize)>,
    pub assigned_cff_source_bounds: Vec<(usize, usize)>,
}

pub struct CutCFF {
    pub terms: BTreeMap<CutCFFIndex, CFFTerm>,
    pub(crate) energy_degree_bound_report: CffEnergyDegreeBoundReport,
    // Terms retain the shared-core contour convention for exact-CFF users and
    // normalization oracles. GammaLoop production consumes this typed bridge
    // only when localizing either the direct-3D or exact-4D route.
    production_prefactor_bridge: CffGlobalPrefactorSign,
}

impl CutCFF {
    pub(crate) const fn production_prefactor_factor(&self) -> i64 {
        self.production_prefactor_bridge.factor()
    }

    pub fn expression_with_selectors(&self) -> Integrands {
        let production_prefactor = Atom::num(self.production_prefactor_factor());
        self.terms
            .iter()
            .map(|(index, term)| {
                (
                    *index,
                    term.expression_with_selectors() * &production_prefactor,
                )
            })
            .collect()
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum CutCffResidueAxis {
    RightThreshold,
    LeftThreshold,
    LuCut,
}

#[derive(Clone, Copy, Default)]
struct SelectedCffResidueAxisParity {
    all: bool,
    variant_local: bool,
}

impl SelectedCffResidueAxisParity {
    fn select(self, variant_local: bool) -> Self {
        Self {
            all: !self.all,
            variant_local: self.variant_local != variant_local,
        }
    }
}

impl CutCffResidueAxis {
    fn set_order(self, index: &mut CutCFFIndex, order: usize) {
        match self {
            Self::RightThreshold => index.right_threshold_order = Some(order),
            Self::LeftThreshold => index.left_threshold_order = Some(order),
            Self::LuCut => index.lu_cut_order = Some(order),
        }
    }
}

fn apply_indexed_residue_selection<F>(
    residues: Vec<(
        CutCFFIndex,
        ThreeDExpression<OrientationID>,
        SelectedCffResidueAxisParity,
    )>,
    axis: CutCffResidueAxis,
    variant_local: bool,
    mut select: F,
) -> Vec<(
    CutCFFIndex,
    ThreeDExpression<OrientationID>,
    SelectedCffResidueAxisParity,
)>
where
    F: FnMut(ThreeDExpression<OrientationID>) -> Vec<ThreeDExpression<OrientationID>>,
{
    residues
        .into_iter()
        .flat_map(|(index, expression, selected_axis_parity)| {
            select(expression)
                .into_iter()
                .enumerate()
                .map(move |(i, residue)| {
                    let mut new_index = index;
                    axis.set_order(&mut new_index, i + 1);
                    // One causal coordinate has been consumed regardless of
                    // which raised-order branch `i` labels.
                    (
                        new_index,
                        residue,
                        selected_axis_parity.select(variant_local),
                    )
                })
        })
        .collect()
}

fn select_indexed_cff_residues_with_variant_local_axes(
    cff: ThreeDExpression<OrientationID>,
    cutset: &CutSet,
    variant_local_axes: &BTreeMap<CutCffResidueAxis, bool>,
) -> Result<
    Vec<(
        CutCFFIndex,
        ThreeDExpression<OrientationID>,
        SelectedCffResidueAxisParity,
    )>,
> {
    if let Some(lu) = cutset.residue_selector.lu.as_ref()
        && (lu.cut_edge_alternatives.is_empty()
            || lu.cut_edge_alternatives.iter().any(Vec::is_empty))
    {
        return Err(eyre::eyre!(
            "LU residue selection requires at least one non-empty physical cut alternative"
        ));
    }
    // Carry selected causal coordinates through the operation pipeline rather
    // than reconstructing them later from raised-order labels in
    // `CutCFFIndex`. The component-local parity distinguishes coordinates
    // whose on-shell-energy factors remain attached to generalized variants.
    let mut residues = vec![(
        CutCFFIndex::new_all_none(),
        cff,
        SelectedCffResidueAxisParity::default(),
    )];

    if let Some(right_threshold) = cutset.residue_selector.right_th_cut.as_ref() {
        residues = apply_indexed_residue_selection(
            residues,
            CutCffResidueAxis::RightThreshold,
            variant_local_axes
                .get(&CutCffResidueAxis::RightThreshold)
                .copied()
                .unwrap_or(false),
            |expression| expression.select_esurface_residue(right_threshold),
        );
    }

    if let Some(left_threshold) = cutset.residue_selector.left_th_cut.as_ref() {
        residues = apply_indexed_residue_selection(
            residues,
            CutCffResidueAxis::LeftThreshold,
            variant_local_axes
                .get(&CutCffResidueAxis::LeftThreshold)
                .copied()
                .unwrap_or(false),
            // Powered poles are completed in the canonical causal basis, so
            // left- and right-threshold selections use the same positive-
            // energy Cutkosky convention as ordinary CFF terms.
            |expression| expression.select_esurface_residue(left_threshold),
        );
    }

    if let Some(lu) = cutset.residue_selector.lu.as_ref() {
        residues = apply_indexed_residue_selection(
            residues,
            CutCffResidueAxis::LuCut,
            variant_local_axes
                .get(&CutCffResidueAxis::LuCut)
                .copied()
                .unwrap_or(false),
            |expression| {
                // Physical support selection is independent of the 3D
                // representation. CFF then consumes the selected surface in
                // GammaLoop's positive-energy Cutkosky convention; future LTD
                // signs and local-series coordinates remain LTD-owned.
                expression
                    .restrict_to_cut_alternatives(&lu.cut_edge_alternatives)
                    .select_esurface_residue(&lu.raised_group)
            },
        );
    }

    Ok(residues)
}

fn select_indexed_cff_residues(
    cff: ThreeDExpression<OrientationID>,
    cutset: &CutSet,
) -> Result<Vec<(CutCFFIndex, ThreeDExpression<OrientationID>, bool)>> {
    Ok(
        select_indexed_cff_residues_with_variant_local_axes(cff, cutset, &BTreeMap::new())?
            .into_iter()
            .map(|(index, expression, parity)| (index, expression, parity.all))
            .collect(),
    )
}

fn variant_local_residue_axes(
    cutset: &CutSet,
    surface_ownership: &BTreeMap<esurface::EsurfaceID, CffEnergyFactorOwnership>,
) -> Result<BTreeMap<CutCffResidueAxis, bool>> {
    let group_is_variant_local = |axis: CutCffResidueAxis,
                                  group: &esurface::RaisedEsurfaceGroup|
     -> Result<bool> {
        let mut ownership = None;
        for surface_id in &group.esurface_ids {
            let Some(candidate) = surface_ownership.get(surface_id).copied() else {
                continue;
            };
            if let Some(previous) = ownership
                && previous != candidate
            {
                return Err(eyre::eyre!(
                    "exact CFF residue axis {axis:?} spans energy-factor components with incompatible ownership"
                ));
            }
            ownership = Some(candidate);
        }
        Ok(ownership == Some(CffEnergyFactorOwnership::VariantLocal))
    };

    let mut axes = BTreeMap::new();
    if let Some(group) = cutset.residue_selector.right_th_cut.as_ref() {
        axes.insert(
            CutCffResidueAxis::RightThreshold,
            group_is_variant_local(CutCffResidueAxis::RightThreshold, group)?,
        );
    }
    if let Some(group) = cutset.residue_selector.left_th_cut.as_ref() {
        axes.insert(
            CutCffResidueAxis::LeftThreshold,
            group_is_variant_local(CutCffResidueAxis::LeftThreshold, group)?,
        );
    }
    if let Some(lu) = cutset.residue_selector.lu.as_ref() {
        axes.insert(
            CutCffResidueAxis::LuCut,
            group_is_variant_local(CutCffResidueAxis::LuCut, &lu.raised_group)?,
        );
    }
    Ok(axes)
}

impl Graph {
    pub(crate) fn cff_from_production_expression(
        &self,
        production: &GeneratedThreeDExpression<esurface::Esurface, hsurface::Hsurface>,
        cutset: &CutSet,
        orientation_pattern: &OrientationPattern,
    ) -> Result<CutCFF> {
        let production_prefactor_bridge = production.core_global_prefactor_sign;
        let mut cff = production.expression.clone();
        normalize_three_d_expression_cut_support_with_raised_edge_groups(
            &mut cff,
            &self.get_raised_edge_groups(),
        );
        let residues = select_indexed_cff_residues(cff, cutset)?;
        let contract_subgraph = self.tree_edges.subtract(&self.initial_state_cut);
        let contract_edges = self
            .iter_edges_of(&contract_subgraph)
            .filter_map(|(pair, edge_id, _)| pair.is_paired().then_some(edge_id))
            .collect::<Vec<_>>();
        let graph_without_is_cut = self
            .underlying
            .full_filter()
            .subtract(&self.initial_state_cut.left)
            .subtract(&self.initial_state_cut.right);
        let cff_loop_number = self
            .get_loop_number()
            .saturating_sub(self.cyclotomatic_number(&contract_subgraph));
        let cff_phase = (-Atom::i()).pow(cff_loop_number as i64);
        let cff_normalization = cff_phase / (Atom::var(GS.pi) * 2).pow(3 * cff_loop_number as i64);
        let cff_energy_factor = match production.energy_factor_ownership {
            CffEnergyFactorOwnership::GlobalSourceProduct => {
                get_cff_inverse_energy_product_impl(self, &graph_without_is_cut, &contract_edges)
            }
            CffEnergyFactorOwnership::VariantLocal => Atom::one(),
        };

        let mut terms = BTreeMap::new();
        for (cut_cff_index, expr, _) in residues {
            let replacement_rules = if cutset.canonicalize_external_shifts {
                expr.surfaces
                    .get_all_replacements_gs_in_lmb(&[], &self.loop_momentum_basis)
            } else {
                expr.surfaces.get_all_replacements_gs(&[])
            };
            let mut cff_term = CFFTerm {
                orientations: Vec::new(),
                exact_source_numerator: None,
            };
            for (orientation_index, orientation) in expr.orientations.into_iter().enumerate() {
                if !orientation_pattern.filter_orientation(&orientation.data.orientation) {
                    continue;
                }
                let expression = orientation
                    .to_atom_gs()
                    .replace_multiple(&replacement_rules)
                    * &cff_energy_factor
                    * &cff_normalization;
                cff_term.orientations.push(CFFOrientationTerm {
                    expression,
                    orientation,
                    production_orientation_id: Some(OrientationID(orientation_index)),
                });
            }
            terms.insert(cut_cff_index, cff_term);
        }

        Ok(CutCFF {
            terms,
            energy_degree_bound_report: CffEnergyDegreeBoundReport {
                source_kind: CffEnergyBoundSourceKind::PhysicalGraph,
                physical_parent_bounds: production.source_energy_degree_bounds.clone(),
                assigned_cff_source_bounds: production.source_energy_degree_bounds.clone(),
            },
            production_prefactor_bridge,
        })
    }

    #[cfg(test)]
    pub(crate) fn cff_from_4d_denominators(
        &mut self,
        denominators: &[FourDDenominator],
        cutset: &CutSet,
        options: &Generate3DExpressionOptions,
        analysis_numerator: &Atom,
    ) -> Result<(CutCFF, linnet::half_edge::subgraph::SuBitGraph)> {
        self.cff_from_4d_denominators_in_uv_edges(
            denominators,
            [],
            cutset,
            options,
            analysis_numerator,
            None,
        )
    }

    pub(crate) fn cff_from_4d_denominators_in_uv_edges(
        &mut self,
        denominators: &[FourDDenominator],
        uv_edges: impl IntoIterator<Item = EdgeIndex>,
        cutset: &CutSet,
        options: &Generate3DExpressionOptions,
        analysis_numerator: &Atom,
        generation_cache: Option<&mut generation::ExactCffGenerationCache>,
    ) -> Result<(CutCFF, linnet::half_edge::subgraph::SuBitGraph)> {
        let (
            generated,
            physical_surfaces,
            physical_energy_edges,
            exact_source_numerator,
            inverse_energy_product,
            cff_loop_number,
            contract_subgraph,
            energy_degree_bound_report,
        ) = {
            let source = GraphThreeDSource::from_exact_denominators_in_uv_edges(
                self,
                denominators,
                uv_edges,
            )?;
            let (
                generated,
                exact_source_energy_mapper,
                energy_assignment,
                energy_degree_bound_report,
            ) = self.generate_3d_expression_for_4d_term(
                &source,
                options,
                analysis_numerator,
                generation_cache,
            )?;
            let physical_surfaces = generated
                .expression
                .surfaces
                .linear_surface_cache
                .iter()
                .map(|surface| source.physical_linear_surface(surface))
                .collect::<Vec<_>>();
            (
                generated,
                physical_surfaces,
                source
                    .physical_energy_edge_index_map()
                    .expect("exact 4D source has a physical energy-edge projection"),
                Arc::new(PlannedExactSourceNumerator {
                    mapper: exact_source_energy_mapper,
                    assignment: energy_assignment,
                }),
                source
                    .exact_inverse_energy_product()
                    .expect("exact 4D source has an occurrence-local energy product"),
                source.active_loop_count(),
                source.contract_subgraph(),
                energy_degree_bound_report,
            )
        };
        let (generated, surface_ownership) =
            self.convert_4d_expression_surfaces(generated, &physical_surfaces)?;
        let energy_factor_ownership = generated.energy_factor_ownership;
        let production_prefactor_bridge = generated.core_global_prefactor_sign;
        crate::debug_tags!(#generation, #cff, #inspect;
            denominator_count = denominators.len(),
            aggregate_ownership = ?energy_factor_ownership,
            components = ?generated.energy_factor_components,
            surface_ownership = ?surface_ownership,
            "Exact CFF energy-factor component metadata"
        );
        let mut cff = generated.expression;
        // Residue support belongs to physical Cutkosky alternatives even
        // though exact numerator maps and half-edge energies remain
        // occurrence-local. Project only that neutral provenance across the
        // dual-ID boundary.
        let physical_edge = |edge: linnet::half_edge::involution::EdgeIndex| {
            let edge_id = usize::from(edge);
            if edge_id < physical_energy_edges.orientation_edge_count {
                return Ok(edge);
            }
            physical_energy_edges
                .internal
                .get(&edge_id)
                .copied()
                .map(linnet::half_edge::involution::EdgeIndex)
                .ok_or_else(|| {
                    eyre::eyre!("exact CFF cut support contains unmapped occurrence edge {edge_id}")
                })
        };
        for orientation in cff.orientations.iter_mut() {
            for variant in &mut orientation.variants {
                variant.denominator_edges = variant
                    .denominator_edges
                    .iter()
                    .copied()
                    .map(&physical_edge)
                    .collect::<Result<Vec<_>>>()?;
                variant.denominator_edge_support_signs =
                    std::mem::take(&mut variant.denominator_edge_support_signs)
                        .into_iter()
                        .try_fold(BTreeMap::new(), |mut mapped, (support, sign)| {
                            let mut support = support
                                .into_iter()
                                .map(&physical_edge)
                                .collect::<Result<Vec<_>>>()?;
                            support.sort_unstable();
                            support.dedup();
                            *mapped.entry(support).or_insert(1) *= sign;
                            Ok::<_, color_eyre::Report>(mapped)
                        })?;
            }
        }
        normalize_three_d_expression_cut_support_with_raised_edge_groups(
            &mut cff,
            &self.get_raised_edge_groups(),
        );
        let variant_local_axes = variant_local_residue_axes(cutset, &surface_ownership)?;
        crate::debug_tags!(#generation, #cff, #inspect;
            variant_local_axes = ?variant_local_axes,
            residue_selector = ?cutset.residue_selector,
            "Exact CFF component-local residue ownership"
        );
        let residues =
            select_indexed_cff_residues_with_variant_local_axes(cff, cutset, &variant_local_axes)?;
        let cff_phase = (-Atom::i()).pow(cff_loop_number as i64);
        let cff_normalization = cff_phase / (Atom::var(GS.pi) * 2).pow(3 * cff_loop_number as i64);
        let mut terms = BTreeMap::new();
        for (cut_cff_index, expr, selected_residue_axis_parity) in residues {
            let cff_energy_factor = match energy_factor_ownership {
                CffEnergyFactorOwnership::GlobalSourceProduct => inverse_energy_product.clone(),
                CffEnergyFactorOwnership::VariantLocal => {
                    // Exact 4D sources define every active occurrence through
                    // 1/(-2E). Generalized CFF variants own the corresponding
                    // magnitudes locally as 1/(+2E), so retain the parity of
                    // the otherwise-discarded parent source product here.
                    // Each selected causal residue consumes one linear
                    // denominator sign, independently of its raised order,
                    // only in the component whose energy factors remain
                    // variant-local. A residue in a separate ordinary CFF
                    // component has already been handled by that component's
                    // pure-CFF selection and must not cancel this source sign.
                    // This also makes a numerator cancellation agree with the
                    // separately generated lower source, whose occurrence
                    // parity changes.
                    let source_parity = !denominators.len().is_multiple_of(2);
                    Atom::num(
                        if source_parity == selected_residue_axis_parity.variant_local {
                            1
                        } else {
                            -1
                        },
                    )
                }
            };
            crate::debug_tags!(#generation, #cff, #inspect;
                cut_index = ?cut_cff_index,
                selected_all_axis_parity = selected_residue_axis_parity.all,
                selected_variant_local_axis_parity = selected_residue_axis_parity.variant_local,
                log.cff_energy_factor = cff_energy_factor,
                "Exact CFF component-local energy-factor bridge"
            );
            let replacement_rules = if cutset.canonicalize_external_shifts {
                expr.surfaces
                    .get_all_replacements_gs_in_lmb(&[], &self.loop_momentum_basis)
            } else {
                expr.surfaces.get_all_replacements_gs(&[])
            };
            let mut cff_term = CFFTerm {
                orientations: Vec::new(),
                exact_source_numerator: Some(exact_source_numerator.clone()),
            };
            for orientation in expr.orientations {
                let expression = orientation
                    .to_atom_gs()
                    .replace_multiple(&replacement_rules)
                    .replace_multiple(exact_source_numerator.mapper.exact_ose_replacements())
                    * &cff_energy_factor
                    * &cff_normalization;
                // Distinct exact factors remain in `expression`. The
                // orientation is deliberately source-local: its affine map
                // evaluates the complete parent numerator directly, rather
                // than being remapped to a production OrientationID.
                cff_term.orientations.push(CFFOrientationTerm {
                    expression,
                    orientation,
                    production_orientation_id: None,
                });
            }
            terms.insert(cut_cff_index, cff_term);
        }
        Ok((
            CutCFF {
                terms,
                energy_degree_bound_report,
                production_prefactor_bridge,
            },
            contract_subgraph,
        ))
    }

    pub fn cff<S: SubGraphLike + SubSetLike>(
        &mut self,
        contract_subgraph: &S,
        cutset: &CutSet,
        orientation_pattern: &OrientationPattern,
        options: &Generate3DExpressionOptions,
        analysis_numerator: Option<&Atom>,
    ) -> Result<CutCFF> {
        let mut contract_edges = vec![];

        for (p, eid, _) in self.iter_edges_of(contract_subgraph) {
            if p.is_paired() {
                contract_edges.push(eid);
            }
        }
        contract_edges.sort_unstable();
        contract_edges.dedup();

        let canonize_esurface = self.get_esurface_canonization(&self.loop_momentum_basis);

        // Reduced UV graphs use the same CFF dispatcher and representation family as the
        // production graph. Switching only this side to a residue coordinate would make its
        // surviving edge-energy maps different numerator samples rather than exact restrictions
        // of production maps. Proper LTD remains a separate deferred representation.
        let generated = self.generate_3d_expression_for_integrand(
            &contract_edges,
            &canonize_esurface,
            options,
            analysis_numerator,
        )?;
        let energy_factor_ownership = generated.energy_factor_ownership;
        let source_energy_degree_bounds = generated.source_energy_degree_bounds.clone();
        let production_prefactor_bridge = generated.core_global_prefactor_sign;
        let mut cff = generated.expression;
        normalize_three_d_expression_cut_support_with_raised_edge_groups(
            &mut cff,
            &self.get_raised_edge_groups(),
        );

        let residues = select_indexed_cff_residues(cff, cutset)?;

        // println!("residue orders: {}", residue.len());

        let graph_without_is_cut = self
            .underlying
            .full_filter()
            .subtract(&self.initial_state_cut.left)
            .subtract(&self.initial_state_cut.right);

        // The CFF carries the measure normalization for the loop variables that remain after
        // contracting a UV subgraph. Fully contracted integrated CTs therefore get no extra CFF
        // measure factor, while ordinary root terms get the full graph-loop factor.
        let cff_loop_number = self
            .get_loop_number()
            .saturating_sub(self.cyclotomatic_number(contract_subgraph));
        let cff_phase = (-Atom::i()).pow(cff_loop_number as i64);
        let cff_normalization = cff_phase / (Atom::var(GS.pi) * 2).pow(3 * cff_loop_number as i64);
        let cff_energy_factor = match energy_factor_ownership {
            CffEnergyFactorOwnership::GlobalSourceProduct => {
                get_cff_inverse_energy_product_impl(self, &graph_without_is_cut, &contract_edges)
            }
            CffEnergyFactorOwnership::VariantLocal => Atom::one(),
        };
        crate::debug_tags!(#cff, #trace;
            stage = "graph_cff_normalization",
            graph = %self.name,
            cff_loop_number = cff_loop_number,
            log.cff_normalization = cff_normalization,
            "Graph CFF normalization"
        );

        let mut terms = BTreeMap::new();

        for (cut_cff_index, expr, _) in residues {
            let replacement_rules = if cutset.canonicalize_external_shifts {
                expr.surfaces
                    .get_all_replacements_gs_in_lmb(&[], &self.loop_momentum_basis)
            } else {
                expr.surfaces.get_all_replacements_gs(&[])
            };
            let mut cff_term = CFFTerm {
                orientations: vec![],
                exact_source_numerator: None,
            };
            for orientation in expr.orientations.iter().filter(|orientation| {
                orientation_pattern.filter_orientation(&orientation.data.orientation)
            }) {
                let eta_expr = orientation.to_atom_gs();
                let mut ose_expr = eta_expr.replace_multiple(&replacement_rules);
                ose_expr *= &cff_energy_factor;
                ose_expr *= cff_normalization.clone();

                crate::debug_tags!(#cff, #trace;
                    stage = "graph_cff_term_expr",
                    graph = %self.name,
                    cut_index = ?cut_cff_index,
                    log.expr = ose_expr,
                    "Graph CFF term expression"
                );
                // println!("ose expr :{}", ose_expr);
                cff_term.orientations.push(CFFOrientationTerm {
                    expression: ose_expr,
                    orientation: orientation.clone(),
                    production_orientation_id: None,
                });
            }
            terms.insert(cut_cff_index, cff_term);
        }

        let cut_cff = CutCFF {
            terms,
            energy_degree_bound_report: CffEnergyDegreeBoundReport {
                source_kind: CffEnergyBoundSourceKind::PhysicalGraph,
                physical_parent_bounds: source_energy_degree_bounds.clone(),
                assigned_cff_source_bounds: source_energy_degree_bounds,
            },
            production_prefactor_bridge,
        };
        Ok(cut_cff)
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use super::*;
    use crate::{dot, graph::parse::IntoGraph, initialisation::test_initialise, utils::W_};
    use itertools::Itertools;
    use linnet::half_edge::involution::{EdgeIndex, Orientation};
    use linnet::half_edge::subgraph::{ModifySubSet, SuBitGraph};
    use symbolica::atom::FunctionBuilder;
    use three_dimensional_reps::ThreeDGraphSource;

    #[test]
    fn lu_residue_selection_rejects_missing_physical_cut_support() {
        let mut cutset = CutSet::empty(0);
        cutset.residue_selector.lu = Some(crate::graph::cuts::LuCutSelection {
            raised_group: crate::cff::esurface::RaisedEsurfaceGroup {
                esurface_ids: Vec::new(),
                max_occurence: 1,
            },
            cut_edge_alternatives: Vec::new(),
        });

        assert!(select_indexed_cff_residues(ThreeDExpression::new_empty(), &cutset).is_err());
    }

    #[test]
    fn indexed_residue_selection_tracks_axis_parity_independent_of_order() -> Result<()> {
        let raised_group = crate::cff::esurface::RaisedEsurfaceGroup {
            esurface_ids: vec![crate::cff::esurface::EsurfaceID::from(0)],
            max_occurence: 2,
        };
        for selected_axis_count in 0usize..=3 {
            let mut cutset = CutSet::empty(1);
            if selected_axis_count >= 1 {
                cutset.residue_selector.right_th_cut = Some(raised_group.clone());
            }
            if selected_axis_count >= 2 {
                cutset.residue_selector.left_th_cut = Some(raised_group.clone());
            }
            if selected_axis_count >= 3 {
                cutset.residue_selector.lu = Some(crate::graph::cuts::LuCutSelection {
                    raised_group: raised_group.clone(),
                    cut_edge_alternatives: vec![vec![EdgeIndex(0)]],
                });
            }

            let residues = select_indexed_cff_residues(ThreeDExpression::new_empty(), &cutset)?;
            assert_eq!(residues.len(), 1usize << selected_axis_count);
            assert!(residues.iter().all(|(index, _, parity)| {
                let present_axes = [
                    index.right_threshold_order,
                    index.left_threshold_order,
                    index.lu_cut_order,
                ]
                .into_iter()
                .flatten()
                .count();
                present_axes == selected_axis_count
                    && *parity == !selected_axis_count.is_multiple_of(2)
            }));
            if selected_axis_count != 0 {
                assert!(
                    residues.iter().any(|(index, _, _)| {
                        [
                            index.right_threshold_order,
                            index.left_threshold_order,
                            index.lu_cut_order,
                        ]
                        .into_iter()
                        .flatten()
                        .any(|order| order == 2)
                    }),
                    "raised order two must retain the same selected-axis parity"
                );
            }
        }
        Ok(())
    }

    #[test]
    fn indexed_residue_selection_tracks_only_variant_local_component_axes() -> Result<()> {
        let raised_group = crate::cff::esurface::RaisedEsurfaceGroup {
            esurface_ids: vec![crate::cff::esurface::EsurfaceID::from(0)],
            max_occurence: 2,
        };
        let mut cutset = CutSet::empty(1);
        cutset.residue_selector.right_th_cut = Some(raised_group.clone());
        cutset.residue_selector.left_th_cut = Some(raised_group.clone());
        cutset.residue_selector.lu = Some(crate::graph::cuts::LuCutSelection {
            raised_group,
            cut_edge_alternatives: vec![vec![EdgeIndex(0)]],
        });
        let variant_local_axes = BTreeMap::from([
            (CutCffResidueAxis::RightThreshold, true),
            (CutCffResidueAxis::LeftThreshold, false),
            (CutCffResidueAxis::LuCut, true),
        ]);

        let residues = select_indexed_cff_residues_with_variant_local_axes(
            ThreeDExpression::new_empty(),
            &cutset,
            &variant_local_axes,
        )?;
        assert_eq!(residues.len(), 8);
        assert!(
            residues
                .iter()
                .all(|(_, _, parity)| { parity.all && !parity.variant_local })
        );
        Ok(())
    }

    #[test]
    fn exact_residue_axis_ownership_rejects_cross_component_groups() -> Result<()> {
        let global = crate::cff::esurface::EsurfaceID::from(0);
        let variant_local = crate::cff::esurface::EsurfaceID::from(1);
        let group = |esurface_ids| crate::cff::esurface::RaisedEsurfaceGroup {
            esurface_ids,
            max_occurence: 1,
        };
        let mut cutset = CutSet::empty(1);
        cutset.residue_selector.right_th_cut = Some(group(vec![global]));
        cutset.residue_selector.lu = Some(crate::graph::cuts::LuCutSelection {
            raised_group: group(vec![variant_local]),
            cut_edge_alternatives: vec![vec![EdgeIndex(0)]],
        });
        let ownership = BTreeMap::from([
            (global, CffEnergyFactorOwnership::GlobalSourceProduct),
            (variant_local, CffEnergyFactorOwnership::VariantLocal),
        ]);

        let axes = variant_local_residue_axes(&cutset, &ownership)?;
        assert_eq!(axes[&CutCffResidueAxis::RightThreshold], false);
        assert_eq!(axes[&CutCffResidueAxis::LuCut], true);

        cutset.residue_selector.left_th_cut = Some(group(vec![global, variant_local]));
        assert!(variant_local_residue_axes(&cutset, &ownership).is_err());
        Ok(())
    }

    #[test]
    fn exact_cff_keeps_dotted_same_edge_occurrences() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph exact_dotted {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]

            incoming -> a [id=2]
            a -> b [id=0 lmb_id=0]
            a -> b [id=1]
            b -> outgoing [id=3]
        })?;
        let edge = EdgeIndex::from(0);
        let momentum = FunctionBuilder::new(GS.emr_mom)
            .add_arg(usize::from(edge))
            .finish();
        let mass_squared = graph.underlying[edge].particle.mass_atom().pow(2);
        let denominators = [
            FourDDenominator {
                source_edge: edge,
                momentum: momentum.clone(),
                mass_squared: mass_squared.clone(),
                full_expr: Atom::var(symbolica::symbol!("exact_cff_test::first")),
            },
            FourDDenominator {
                source_edge: edge,
                momentum,
                mass_squared: mass_squared.clone(),
                full_expr: Atom::var(symbolica::symbol!("exact_cff_test::second")),
            },
        ];
        let options = graph.denominator_only_cff_3d_expression_options();
        let cutset = CutSet::empty(graph.n_hedges());

        let (cff, _) =
            graph.cff_from_4d_denominators(&denominators, &cutset, &options, &Atom::one())?;
        assert_eq!(cff.terms.len(), 1);
        let on_shell_energy = (1..=3)
            .fold(mass_squared, |norm_squared, spatial_index| {
                norm_squared + GS.emr_mom(edge, GS.cind(spatial_index)).pow(2)
            })
            .sqrt();
        // Ordinary CFF distributes the double-pole residue across orientations,
        // so only the complete explicit sum has the contour normalization.
        let expected_orientation_sum =
            Atom::i() / (Atom::num(32) * Atom::var(GS.pi).pow(3) * on_shell_energy.pow(3));
        let occurrence_offset = graph.underlying.n_edges();
        let first_occurrence = EdgeIndex(occurrence_offset);
        let second_occurrence = EdgeIndex(occurrence_offset + 1);
        let terms = &cff
            .terms
            .values()
            .next()
            .expect("the empty cutset has one CFF term")
            .orientations;
        assert_eq!(terms.len(), 2);
        for term in terms {
            let orientation = &term.orientation.data.orientation;
            assert_ne!(orientation[first_occurrence], Orientation::Undirected);
            assert_ne!(orientation[second_occurrence], Orientation::Undirected);
            assert_ne!(
                orientation[second_occurrence],
                orientation[first_occurrence]
            );
            assert_eq!(orientation[edge], Orientation::Undirected);
        }
        assert_ne!(
            terms[0].orientation.data.orientation[first_occurrence],
            terms[1].orientation.data.orientation[first_occurrence]
        );
        let orientation_sum = terms
            .iter()
            .fold(Atom::zero(), |sum, term| sum + &term.expression)
            .replace(GS.ose(edge))
            .with(on_shell_energy);
        assert!(
            (orientation_sum - expected_orientation_sum)
                .together()
                .is_zero()
        );
        Ok(())
    }

    #[test]
    fn exact_cff_uncancelled_powered_denominator_matches_lower_source() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph exact_powered_identity {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0 lmb_id=0]
            a -> b [id=1]
        })?;
        let edge = EdgeIndex(0);
        let momentum = FunctionBuilder::new(GS.emr_mom)
            .add_arg(usize::from(edge))
            .finish();
        let mass_squared = graph.underlying[edge].particle.mass_atom().pow(2);
        let full_expr = GS.emr_mom(edge, GS.cind(0)).pow(2)
            - (1..=3).fold(mass_squared.clone(), |norm_squared, spatial_index| {
                norm_squared + GS.emr_mom(edge, GS.cind(spatial_index)).pow(2)
            });
        let denominator = FourDDenominator {
            source_edge: edge,
            momentum: momentum.clone(),
            mass_squared: mass_squared.clone(),
            full_expr: full_expr.clone(),
        };
        let spectator_edge = EdgeIndex(1);
        let spectator_momentum = FunctionBuilder::new(GS.emr_mom)
            .add_arg(usize::from(spectator_edge))
            .finish();
        let spectator_mass_squared = graph.underlying[spectator_edge].particle.mass_atom().pow(2);
        let spectator = FourDDenominator {
            source_edge: spectator_edge,
            momentum: spectator_momentum,
            mass_squared: spectator_mass_squared.clone(),
            full_expr: GS.emr_mom(spectator_edge, GS.cind(0)).pow(2)
                - (1..=3).fold(spectator_mass_squared, |norm_squared, spatial_index| {
                    norm_squared + GS.emr_mom(spectator_edge, GS.cind(spatial_index)).pow(2)
                }),
        };
        let denominators = [spectator, denominator.clone(), denominator];
        let retained_constant = Atom::var(symbolica::symbol!("exact_cff_test::retained_factor"));
        let retained_factor = GS.emr_mom(edge, GS.cind(0)) + &retained_constant;
        let numerator =
            GS.den(usize::from(edge), momentum, mass_squared, full_expr) * &retained_factor;
        let options = graph.denominator_only_cff_3d_expression_options();
        let cutset = CutSet::empty(graph.n_hedges());

        let (powered, _) =
            graph.cff_from_4d_denominators(&denominators, &cutset, &options, &numerator)?;
        let powered_prefactor = Atom::num(powered.production_prefactor_factor());
        let powered_sum = powered
            .terms
            .values()
            .flat_map(|term| {
                term.orientations.iter().map(|orientation| {
                    Ok(orientation.expression.clone()
                        * term.map_exact_source_numerator(&orientation.orientation)?)
                })
            })
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .fold(Atom::Zero, |sum, term| sum + term)
            * powered_prefactor;
        let (lower, _) = graph.cff_from_4d_denominators(
            &denominators[..2],
            &cutset,
            &options,
            &retained_factor,
        )?;
        let lower_prefactor = Atom::num(lower.production_prefactor_factor());
        let lower_sum = lower
            .terms
            .values()
            .flat_map(|term| {
                term.orientations.iter().map(|orientation| {
                    Ok(orientation.expression.clone()
                        * term.map_exact_source_numerator(&orientation.orientation)?)
                })
            })
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .fold(Atom::Zero, |sum, term| sum + term)
            * lower_prefactor;
        let contract: SuBitGraph = graph.empty_subgraph();
        let ordinary = graph.cff(
            &contract,
            &cutset,
            &OrientationPattern::default(),
            &options,
            Some(&retained_factor),
        )?;
        assert_eq!(
            lower.production_prefactor_factor(),
            ordinary.production_prefactor_factor(),
            "exact and ordinary lower sources must use the same production prefactor bridge"
        );
        let ordinary_prefactor = Atom::num(ordinary.production_prefactor_factor());
        let ordinary_sum = ordinary
            .terms
            .values()
            .flat_map(|term| {
                term.orientations.iter().map(|orientation| {
                    orientation.expression.clone()
                        * retained_factor.replace_multiple(
                            orientation.orientation.energy_replacements_gs(&graph),
                        )
                })
            })
            .fold(Atom::Zero, |sum, term| sum + term)
            * ordinary_prefactor;
        let powered_sum = powered_sum
            .replace(GS.den(W_.a_, W_.b_, W_.c_, W_.d_))
            .with(W_.d_);

        let fixed_point = |mut expression: Atom, spatial: &Atom, energy: &Atom, constant: &Atom| {
            for source_edge in [edge, spectator_edge] {
                expression = expression
                    .replace(GS.emr_mom(source_edge, GS.cind(1)))
                    .with(spatial.clone())
                    .replace(GS.ose(source_edge))
                    .with(energy.clone());
                for spatial_index in 2..=3 {
                    expression = expression
                        .replace(GS.emr_mom(source_edge, GS.cind(spatial_index)))
                        .with(Atom::Zero);
                }
            }
            expression
                .replace(retained_constant.clone())
                .with(constant.clone())
        };
        let points = [
            (Atom::Zero, Atom::one(), Atom::num(2)),
            (
                Atom::num(symbolica::domains::rational::Rational::from((3, 4))),
                Atom::num(symbolica::domains::rational::Rational::from((5, 4))),
                Atom::num(symbolica::domains::rational::Rational::from((7, 3))),
            ),
        ];
        for (spatial, energy, constant) in points {
            let difference = (fixed_point(powered_sum.clone(), &spatial, &energy, &constant)
                - fixed_point(lower_sum.clone(), &spatial, &energy, &constant))
            .together();
            assert!(
                difference.is_zero(),
                "uncancelled D8*(Q0+c)/(D7*D8^2) exact CFF differs from (Q0+c)/(D7*D8) at spatial={spatial}, c={constant}: {difference}"
            );
            let normalization_difference =
                (fixed_point(lower_sum.clone(), &spatial, &energy, &constant)
                    - fixed_point(ordinary_sum.clone(), &spatial, &energy, &constant))
                .together();
            assert!(
                normalization_difference.is_zero(),
                "exact and ordinary (Q0+c)/(D7*D8) CFF sums differ at spatial={spatial}, c={constant}: {normalization_difference}"
            );
        }
        Ok(())
    }

    #[test]
    fn exact_cff_uncancelled_powered_denominator_matches_lower_lu_residues() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph exact_powered_lu_identity {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]

            incoming -> v1 [id=0]
            v1 -> v2 [id=1 lmb_id=0]
            v2 -> v3 [id=2]
            v1 -> v3 [id=3]
            v3 -> outgoing [id=4]
        })?;
        let edge = EdgeIndex(1);
        let momentum = FunctionBuilder::new(GS.emr_mom)
            .add_arg(usize::from(edge))
            .finish();
        let mass_squared = graph.underlying[edge].particle.mass_atom().pow(2);
        let full_expr = GS.emr_mom(edge, GS.cind(0)).pow(2)
            - (1..=3).fold(mass_squared.clone(), |norm_squared, spatial_index| {
                norm_squared + GS.emr_mom(edge, GS.cind(spatial_index)).pow(2)
            });
        let denominator = FourDDenominator {
            source_edge: edge,
            momentum: momentum.clone(),
            mass_squared: mass_squared.clone(),
            full_expr: full_expr.clone(),
        };
        let spectator_edge = EdgeIndex(3);
        let spectator_mass_squared = graph.underlying[spectator_edge].particle.mass_atom().pow(2);
        let spectator = FourDDenominator {
            source_edge: spectator_edge,
            momentum: FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(spectator_edge))
                .finish(),
            mass_squared: spectator_mass_squared.clone(),
            full_expr: GS.emr_mom(spectator_edge, GS.cind(0)).pow(2)
                - (1..=3).fold(spectator_mass_squared, |norm_squared, spatial_index| {
                    norm_squared + GS.emr_mom(spectator_edge, GS.cind(spatial_index)).pow(2)
                }),
        };
        let denominators = [spectator, denominator.clone(), denominator];
        let retained_constant = Atom::var(symbolica::symbol!("exact_cff_test::lu_retained_factor"));
        let retained_factor = GS.emr_mom(edge, GS.cind(0)) + &retained_constant;
        let numerator =
            GS.den(usize::from(edge), momentum, mass_squared, full_expr) * &retained_factor;
        let options = graph.denominator_only_cff_3d_expression_options();
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let production = graph.generate_3d_expression_for_integrand(
            &[EdgeIndex(2)],
            &canonization,
            &options,
            Some(&retained_factor),
        )?;
        let mut lu_cut = graph
            .determine_raised_esurfaces_from_expression(&production.expression)
            .raised_groups
            .into_iter()
            .find(|group| {
                group.esurface_ids.iter().any(|esurface_id| {
                    !production.expression.surfaces.esurface_cache[*esurface_id]
                        .external_shift
                        .is_empty()
                })
            })
            .expect("the lower production CFF contains an ordinary physical LU surface");
        lu_cut.max_occurence = 2;
        let mut cutset = CutSet::empty(graph.n_hedges());
        cutset.residue_selector.lu = Some(crate::graph::cuts::LuCutSelection {
            raised_group: lu_cut.clone(),
            cut_edge_alternatives: lu_cut
                .esurface_ids
                .iter()
                .map(|esurface_id| {
                    production.expression.surfaces.esurface_cache[*esurface_id]
                        .energies
                        .clone()
                })
                .collect(),
        });

        let (powered, _) =
            graph.cff_from_4d_denominators(&denominators, &cutset, &options, &numerator)?;
        let (lower, _) = graph.cff_from_4d_denominators(
            &denominators[..2],
            &cutset,
            &options,
            &retained_factor,
        )?;
        let mut contract: SuBitGraph = graph.empty_subgraph();
        contract.add(graph[&EdgeIndex(2)].1);
        let ordinary = graph.cff(
            &contract,
            &cutset,
            &OrientationPattern::default(),
            &options,
            Some(&retained_factor),
        )?;
        assert_eq!(
            powered.terms.keys().collect::<Vec<_>>(),
            lower.terms.keys().collect::<Vec<_>>(),
            "powered and lower exact sources must expose the same LU residues"
        );
        assert_eq!(
            lower.terms.keys().collect::<Vec<_>>(),
            ordinary.terms.keys().collect::<Vec<_>>(),
            "exact and ordinary sources must expose the same LU residues"
        );

        let points = [
            (
                Atom::Zero,
                Atom::one(),
                Atom::num(3) / 4,
                Atom::num(5) / 4,
                Atom::num(7),
                Atom::num(2),
            ),
            (
                Atom::num(3) / 4,
                Atom::num(5) / 4,
                Atom::Zero,
                Atom::one(),
                Atom::num(11),
                Atom::num(7) / 3,
            ),
        ];
        for index in lower.terms.keys() {
            let powered_term = powered.terms.get(index).expect("powered residue exists");
            let powered_sum = powered_term
                .orientations
                .iter()
                .map(|orientation| {
                    Ok(orientation.expression.clone()
                        * powered_term.map_exact_source_numerator(&orientation.orientation)?)
                })
                .collect::<Result<Vec<_>>>()?
                .into_iter()
                .fold(Atom::Zero, |sum, term| sum + term)
                * Atom::num(powered.production_prefactor_factor());
            let powered_sum = powered_sum
                .replace(GS.den(W_.a_, W_.b_, W_.c_, W_.d_))
                .with(W_.d_);
            let lower_term = lower.terms.get(index).expect("lower residue exists");
            let lower_sum = lower_term
                .orientations
                .iter()
                .map(|orientation| {
                    Ok(orientation.expression.clone()
                        * lower_term.map_exact_source_numerator(&orientation.orientation)?)
                })
                .collect::<Result<Vec<_>>>()?
                .into_iter()
                .fold(Atom::Zero, |sum, term| sum + term)
                * Atom::num(lower.production_prefactor_factor());
            let ordinary_sum = ordinary
                .terms
                .get(index)
                .expect("ordinary residue exists")
                .orientations
                .iter()
                .map(|orientation| {
                    orientation.expression.clone()
                        * retained_factor.replace_multiple(
                            orientation.orientation.energy_replacements_gs(&graph),
                        )
                })
                .fold(Atom::Zero, |sum, term| sum + term)
                * Atom::num(ordinary.production_prefactor_factor());

            for (q, eq, q3, e3, external_energy, constant) in &points {
                let fixed_point = |mut expression: Atom| {
                    for (source_edge, spatial, energy) in [(edge, q, eq), (spectator_edge, q3, e3)]
                    {
                        expression = expression
                            .replace(GS.emr_mom(source_edge, GS.cind(1)))
                            .with(spatial.clone())
                            .replace(GS.ose(source_edge))
                            .with(energy.clone());
                        for spatial_index in 2..=3 {
                            expression = expression
                                .replace(GS.emr_mom(source_edge, GS.cind(spatial_index)))
                                .with(Atom::Zero);
                        }
                    }
                    for external_edge in [EdgeIndex(0), EdgeIndex(4)] {
                        expression = expression
                            .replace(GS.emr_mom(external_edge, GS.cind(0)))
                            .with(external_energy.clone());
                    }
                    expression
                        .replace(retained_constant.clone())
                        .with(constant.clone())
                };
                let powered_value = fixed_point(powered_sum.clone());
                let lower_value = fixed_point(lower_sum.clone());
                let ordinary_value = fixed_point(ordinary_sum.clone());
                let powered_difference = (powered_value - &lower_value).together();
                assert!(
                    powered_difference.is_zero(),
                    "uncancelled powered and lower exact LU residues differ for index {index}: {powered_difference}"
                );
                let ordinary_difference = (ordinary_value - &lower_value).together();
                assert!(
                    ordinary_difference.is_zero(),
                    "exact lower and ordinary LU residues differ for index {index}: {ordinary_difference}"
                );
                if index.lu_cut_order == Some(2) {
                    assert!(
                        lower_value.together().is_zero(),
                        "the artificially raised second-order LU residue must vanish for the lower source"
                    );
                }
            }
        }
        Ok(())
    }

    #[test]
    fn exact_cff_alias_edges_share_physical_numerator_energy() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph exact_alias_energy_identity {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]

            incoming -> v1 [id=0]
            v1 -> v2 [id=1 lmb_id=0]
            v2 -> v3 [id=2]
            v1 -> v3 [id=3]
            v3 -> outgoing [id=4]
        })?;
        let alias_edges = [EdgeIndex(1), EdgeIndex(2)];
        let spectator_edge = EdgeIndex(3);
        let denominators = [spectator_edge, alias_edges[0], alias_edges[1]].map(|edge| {
            let momentum = FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(edge))
                .finish();
            let mass_squared = graph.underlying[edge].particle.mass_atom().pow(2);
            FourDDenominator {
                source_edge: edge,
                momentum,
                mass_squared: mass_squared.clone(),
                full_expr: GS.emr_mom(edge, GS.cind(0)).pow(2)
                    - (1..=3).fold(mass_squared, |norm_squared, spatial_index| {
                        norm_squared + GS.emr_mom(edge, GS.cind(spatial_index)).pow(2)
                    }),
            }
        });
        let retained_constant =
            Atom::var(symbolica::symbol!("exact_cff_test::alias_retained_factor"));
        let numerator = (GS.emr_mom(alias_edges[0], GS.cind(0)) + &retained_constant)
            * (GS.emr_mom(alias_edges[1], GS.cind(0)) + &retained_constant + 1);
        let options = graph.denominator_only_cff_3d_expression_options();
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let production = graph.generate_3d_expression_for_integrand(
            &[],
            &canonization,
            &options,
            Some(&numerator),
        )?;
        let mut lu_cut = graph
            .determine_raised_esurfaces_from_expression(&production.expression)
            .raised_groups
            .into_iter()
            .find(|group| {
                group.esurface_ids.iter().any(|esurface_id| {
                    !production.expression.surfaces.esurface_cache[*esurface_id]
                        .external_shift
                        .is_empty()
                })
            })
            .expect("the production CFF contains a physical LU surface");
        lu_cut.max_occurence = 2;
        let mut cutset = CutSet::empty(graph.n_hedges());
        cutset.residue_selector.lu = Some(crate::graph::cuts::LuCutSelection {
            raised_group: lu_cut.clone(),
            cut_edge_alternatives: lu_cut
                .esurface_ids
                .iter()
                .map(|esurface_id| {
                    production.expression.surfaces.esurface_cache[*esurface_id]
                        .energies
                        .clone()
                })
                .collect(),
        });

        let (exact, _) =
            graph.cff_from_4d_denominators(&denominators, &cutset, &options, &numerator)?;
        let contract: SuBitGraph = graph.empty_subgraph();
        let ordinary = graph.cff(
            &contract,
            &cutset,
            &OrientationPattern::default(),
            &options,
            Some(&numerator),
        )?;
        assert_eq!(
            exact.terms.keys().collect::<Vec<_>>(),
            ordinary.terms.keys().collect::<Vec<_>>(),
            "exact and ordinary alias sources must expose the same LU residues"
        );

        let points = [
            (
                Atom::Zero,
                Atom::one(),
                Atom::num(3) / 4,
                Atom::num(5) / 4,
                Atom::num(7),
                Atom::num(2),
            ),
            (
                Atom::num(3) / 4,
                Atom::num(5) / 4,
                Atom::Zero,
                Atom::one(),
                Atom::num(11),
                Atom::num(7) / 3,
            ),
        ];
        for index in ordinary.terms.keys() {
            let exact_term = exact.terms.get(index).expect("exact residue exists");
            let exact_sum = exact_term
                .orientations
                .iter()
                .map(|orientation| {
                    Ok(orientation.expression.clone()
                        * exact_term.map_exact_source_numerator(&orientation.orientation)?)
                })
                .collect::<Result<Vec<_>>>()?
                .into_iter()
                .fold(Atom::Zero, |sum, term| sum + term)
                * Atom::num(exact.production_prefactor_factor());
            let ordinary_sum = ordinary
                .terms
                .get(index)
                .expect("ordinary residue exists")
                .orientations
                .iter()
                .map(|orientation| {
                    orientation.expression.clone()
                        * numerator.replace_multiple(
                            orientation.orientation.energy_replacements_gs(&graph),
                        )
                })
                .fold(Atom::Zero, |sum, term| sum + term)
                * Atom::num(ordinary.production_prefactor_factor());

            for (alias_spatial, alias_energy, spectator_spatial, spectator_energy, external, c) in
                &points
            {
                let fixed_point = |mut expression: Atom| {
                    for edge in alias_edges {
                        expression = expression
                            .replace(GS.emr_mom(edge, GS.cind(1)))
                            .with(alias_spatial.clone())
                            .replace(GS.ose(edge))
                            .with(alias_energy.clone());
                    }
                    expression = expression
                        .replace(GS.emr_mom(spectator_edge, GS.cind(1)))
                        .with(spectator_spatial.clone())
                        .replace(GS.ose(spectator_edge))
                        .with(spectator_energy.clone());
                    for edge in [alias_edges[0], alias_edges[1], spectator_edge] {
                        for spatial_index in 2..=3 {
                            expression = expression
                                .replace(GS.emr_mom(edge, GS.cind(spatial_index)))
                                .with(Atom::Zero);
                        }
                    }
                    for external_edge in [EdgeIndex(0), EdgeIndex(4)] {
                        expression = expression
                            .replace(GS.emr_mom(external_edge, GS.cind(0)))
                            .with(external.clone());
                    }
                    expression
                        .replace(retained_constant.clone())
                        .with(c.clone())
                };
                let difference =
                    (fixed_point(exact_sum.clone()) - fixed_point(ordinary_sum.clone())).together();
                assert!(
                    difference.is_zero(),
                    "exact alias-edge and ordinary LU residues differ for index {index}: {difference}"
                );
            }
        }
        Ok(())
    }

    #[test]
    fn exact_cff_cubic_uv_rewrite_matches_production_convention() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph exact_cubic_uv_rewrite {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]

            incoming -> v5 [id=0]
            v0 -> v2 [id=1 lmb_id=0]
            v3 -> v0 [id=2]
            v0 -> v5 [id=3]
            v2 -> v1 [id=4]
            v1 -> v3 [id=5 lmb_id=1]
            v1 -> v4 [id=6]
            v2 -> v3 [id=7]
            v4 -> outgoing [id=8]
        })?;
        let owner_edges = [EdgeIndex(4), EdgeIndex(5), EdgeIndex(7)];
        assert_eq!(
            owner_edges
                .iter()
                .map(|edge| &graph.loop_momentum_basis.edge_signatures[*edge])
                .collect::<BTreeSet<_>>()
                .len(),
            owner_edges.len(),
            "the exact UV owners must have genuinely different physical routings"
        );

        let carrier = EdgeIndex(5);
        let momentum = FunctionBuilder::new(GS.emr_mom)
            .add_arg(usize::from(carrier))
            .finish();
        let mass_squared = Atom::var(GS.m_uv_expansion).pow(2);
        let full_expr = GS.emr_mom(carrier, GS.cind(0)).pow(2)
            - (1..=3).fold(mass_squared.clone(), |norm_squared, spatial_index| {
                norm_squared + GS.emr_mom(carrier, GS.cind(spatial_index)).pow(2)
            });
        let denominator_records = owner_edges
            .into_iter()
            .zip([1, 1, -1])
            .map(|(source_edge, routing_sign)| FourDDenominator {
                source_edge,
                momentum: Atom::num(routing_sign) * &momentum,
                mass_squared: mass_squared.clone(),
                full_expr: full_expr.clone(),
            })
            .collect::<Vec<_>>();
        let mut reference_parsed: Option<three_dimensional_reps::ParsedGraph> = None;
        let mut reference_orientation_sum: Option<Atom> = None;
        let mut reference_quadratic_sum: Option<Atom> = None;
        for denominators in denominator_records.iter().cloned().permutations(3) {
            let source = GraphThreeDSource::from_exact_denominators_in_uv_edges(
                &graph,
                &denominators,
                owner_edges,
            )?;
            let parsed = source.to_three_d_parsed_graph()?;
            if let Some(reference) = &reference_parsed {
                assert_eq!(
                    &parsed, reference,
                    "canonical exact topology must be invariant under denominator-factor order"
                );
            } else {
                reference_parsed = Some(parsed.clone());
            }

            let validation = three_dimensional_reps::validate_parsed_graph(&parsed);
            assert!(
                validation.ok,
                "the cubic UV source must remain momentum-balanced: {validation:?}"
            );
            let nodes = parsed
                .internal_edges
                .iter()
                .flat_map(|edge| [edge.tail, edge.head])
                .collect::<BTreeSet<_>>();
            assert_eq!(
                parsed.internal_edges.len() + 1 - nodes.len(),
                1,
                "the cubic UV source must retain one topological loop"
            );
            assert_eq!(parsed.loop_names.len(), 1);
            assert!(
                parsed.external_edges.is_empty(),
                "zero-denominator cograph components must not create causal surfaces"
            );

            let options = graph.denominator_only_cff_3d_expression_options();
            let cutset = CutSet::empty(graph.n_hedges());
            let (cff, _) = graph.cff_from_4d_denominators_in_uv_edges(
                &denominators,
                owner_edges,
                &cutset,
                &options,
                &Atom::one(),
                None,
            )?;
            assert_eq!(cff.terms.len(), 1);
            let term = cff
                .terms
                .values()
                .next()
                .expect("the empty cutset has one exact CFF term");
            let orientation_sum = term
                .orientations
                .iter()
                .fold(Atom::Zero, |sum, orientation| sum + &orientation.expression)
                * Atom::num(cff.production_prefactor_factor());
            let energy = (1..=3)
                .fold(mass_squared.clone(), |norm_squared, spatial_index| {
                    norm_squared + GS.emr_mom(carrier, GS.cind(spatial_index)).pow(2)
                })
                .sqrt();
            // The raw mathematical q0 contour of (q0²-E²)^-3 is
            // -3i/(128π³E⁵), including the spatial (2π)^-3 normalization.
            // GammaLoop production shares the ordinary CFF convention
            // 1/prod(-2E), whose odd three-denominator product supplies the
            // additional minus sign. This test locks that production convention;
            // it does not choose a physical phase convention.
            let production_contour = Atom::num(3) * Atom::i()
                / (Atom::num(128) * Atom::var(GS.pi).pow(3) * energy.pow(5));
            let difference = (&orientation_sum - &production_contour).together();
            assert!(
                difference.is_zero(),
                "the cubic exact UV CFF differs from the production convention: {difference}"
            );
            if let Some(reference) = &reference_orientation_sum {
                let provenance_difference = (&orientation_sum - reference).together();
                assert!(
                    provenance_difference.is_zero(),
                    "the bare cubic exact UV CFF depends on denominator-factor order: {provenance_difference}"
                );
            } else {
                reference_orientation_sum = Some(orientation_sum.clone());
            }

            let quadratic_numerator = GS.emr_mom(carrier, GS.cind(0)).pow(2);
            let (quadratic_cff, _) = graph.cff_from_4d_denominators_in_uv_edges(
                &denominators,
                owner_edges,
                &cutset,
                &options,
                &quadratic_numerator,
                None,
            )?;
            let quadratic_plan = quadratic_cff
                .terms
                .values()
                .next()
                .and_then(|term| term.exact_source_numerator.as_ref())
                .expect("the exact quadratic CFF term retains its numerator plan");
            assert_eq!(
                quadratic_plan
                    .assignment
                    .energy_degree_bounds()
                    .iter()
                    .map(|(_, degree)| *degree)
                    .collect::<Vec<_>>(),
                vec![1, 1],
                "the quadratic energy must stay factorized over two equivalent occurrences"
            );
            let quadratic_sum = quadratic_cff
                .terms
                .values()
                .flat_map(|term| {
                    term.orientations.iter().map(|orientation| {
                        Ok(orientation.expression.clone()
                            * term.map_exact_source_numerator(&orientation.orientation)?)
                    })
                })
                .collect::<Result<Vec<_>>>()?
                .into_iter()
                .fold(Atom::Zero, |sum, term| sum + term)
                * Atom::num(quadratic_cff.production_prefactor_factor());
            // Under the parent three-occurrence convention, q0² = D + E²
            // combines the double- and triple-pole sectors with their relative
            // occurrence parity. The result is -i/(128 pi³ E³), not their
            // unsigned sum 7i/(128 pi³ E³).
            let quadratic_contour =
                -Atom::i() / (Atom::num(128) * Atom::var(GS.pi).pow(3) * energy.pow(3));
            let quadratic_difference = (&quadratic_sum - &quadratic_contour).together();
            assert!(
                quadratic_difference.is_zero(),
                "the quadratic cubic-pole exact CFF loses its lower-sector occurrence parity: {quadratic_difference}"
            );
            if let Some(reference) = &reference_quadratic_sum {
                let provenance_difference = (&quadratic_sum - reference).together();
                assert!(
                    provenance_difference.is_zero(),
                    "the q0^2 cubic exact UV CFF depends on denominator-factor order: {provenance_difference}"
                );
            } else {
                reference_quadratic_sum = Some(quadratic_sum.clone());
            }

            let spatial_norm = (1..=3).fold(Atom::Zero, |norm_squared, spatial_index| {
                norm_squared + GS.emr_mom(carrier, GS.cind(spatial_index)).pow(2)
            });
            let momentum_squared_numerator = &quadratic_numerator - &spatial_norm;
            let (momentum_squared_cff, _) = graph.cff_from_4d_denominators_in_uv_edges(
                &denominators,
                owner_edges,
                &cutset,
                &options,
                &momentum_squared_numerator,
                None,
            )?;
            let momentum_squared_sum = momentum_squared_cff
                .terms
                .values()
                .flat_map(|term| {
                    term.orientations.iter().map(|orientation| {
                        Ok(orientation.expression.clone()
                            * term.map_exact_source_numerator(&orientation.orientation)?)
                    })
                })
                .collect::<Result<Vec<_>>>()?
                .into_iter()
                .fold(Atom::Zero, |sum, term| sum + term)
                * Atom::num(momentum_squared_cff.production_prefactor_factor());
            let momentum_squared_contour = &quadratic_contour - &spatial_norm * &production_contour;
            let momentum_squared_difference = (momentum_squared_sum - momentum_squared_contour)
                .together()
                .expand();
            assert!(
                momentum_squared_difference.is_zero(),
                "the full momentum-square numerator must retain its temporal/spatial cancellation: {momentum_squared_difference}"
            );
        }
        Ok(())
    }

    #[test]
    fn exact_lu_cut_matches_ordinary_cff_per_residue() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph exact_lu_cut_identity {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]

            incoming -> v1 [id=0]
            v1 -> v2 [id=1 lmb_id=0]
            v1 -> v2 [id=2]
            v2 -> outgoing [id=3]
        })?;
        let options = graph.denominator_only_cff_3d_expression_options();
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let denominators = [1, 2].map(|edge| {
            let edge = EdgeIndex(edge);
            FourDDenominator {
                source_edge: edge,
                momentum: FunctionBuilder::new(GS.emr_mom)
                    .add_arg(usize::from(edge))
                    .finish(),
                mass_squared: graph.underlying[edge].particle.mass_atom().pow(2),
                full_expr: Atom::one(),
            }
        });
        let contract: SuBitGraph = graph.empty_subgraph();

        for (expected_order, numerator) in [(1, Atom::one())] {
            let production = graph.generate_3d_expression_for_integrand(
                &[],
                &canonization,
                &options,
                Some(&numerator),
            )?;
            let raised_groups = graph
                .determine_raised_esurfaces_from_expression(&production.expression)
                .raised_groups;
            let lu_cut = raised_groups
                .into_iter()
                .find(|group| {
                    group.max_occurence == expected_order
                        && group.esurface_ids.iter().any(|esurface_id| {
                            !production.expression.surfaces.esurface_cache[*esurface_id]
                                .external_shift
                                .is_empty()
                        })
                })
                .expect("the production CFF contains the requested physical LU surface");
            let mut cutset = CutSet::empty(graph.n_hedges());
            cutset.residue_selector.lu = Some(crate::graph::cuts::LuCutSelection {
                raised_group: lu_cut.clone(),
                cut_edge_alternatives: lu_cut
                    .esurface_ids
                    .iter()
                    .map(|esurface_id| {
                        production.expression.surfaces.esurface_cache[*esurface_id]
                            .energies
                            .clone()
                    })
                    .collect(),
            });
            let ordinary = graph.cff(
                &contract,
                &cutset,
                &OrientationPattern::default(),
                &options,
                Some(&numerator),
            )?;
            let (exact, _) =
                graph.cff_from_4d_denominators(&denominators, &cutset, &options, &numerator)?;
            assert_eq!(
                exact.production_prefactor_factor(),
                ordinary.production_prefactor_factor()
            );
            assert_eq!(
                exact.terms.keys().collect::<Vec<_>>(),
                ordinary.terms.keys().collect::<Vec<_>>()
            );
            for (index, ordinary_term) in &ordinary.terms {
                let mut ordinary_sum = ordinary_term
                    .orientations
                    .iter()
                    .map(|orientation| {
                        orientation.expression.clone()
                            * numerator.replace_multiple(
                                orientation.orientation.energy_replacements_gs(&graph),
                            )
                    })
                    .fold(Atom::Zero, |sum, term| sum + term)
                    * Atom::num(ordinary.production_prefactor_factor());
                for edge in [EdgeIndex(1), EdgeIndex(2)] {
                    let on_shell_energy = (1..=3)
                        .fold(
                            graph.underlying[edge].particle.mass_atom().pow(2),
                            |norm_squared, spatial_index| {
                                norm_squared + GS.emr_mom(edge, GS.cind(spatial_index)).pow(2)
                            },
                        )
                        .sqrt();
                    ordinary_sum = ordinary_sum.replace(GS.ose(edge)).with(on_shell_energy);
                }
                let exact_term = exact
                    .terms
                    .get(index)
                    .expect("exact and ordinary residue keys agree");
                let exact_sum = exact_term
                    .orientations
                    .iter()
                    .map(|orientation| {
                        Ok(orientation.expression.clone()
                            * exact_term.map_exact_source_numerator(&orientation.orientation)?)
                    })
                    .collect::<Result<Vec<_>>>()?
                    .into_iter()
                    .fold(Atom::Zero, |sum, term| sum + term)
                    * Atom::num(exact.production_prefactor_factor());
                let difference = (exact_sum - ordinary_sum).together();
                assert!(
                    difference.is_zero(),
                    "exact and ordinary LU residues differ for maximum order {expected_order}, index {index}: {difference}"
                );
            }
        }
        Ok(())
    }

    #[test]
    fn exact_raised_lu_cut_matches_ordinary_cff_per_residue() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph exact_raised_lu_cut_identity {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]

            incoming -> v1 [id=0]
            v1 -> v2 [id=1 lmb_id=0]
            v2 -> v3 [id=2]
            v1 -> v3 [id=3]
            v3 -> outgoing [id=4]
        })?;
        let options = graph.denominator_only_cff_3d_expression_options();
        let numerator = GS.emr_mom(EdgeIndex(1), GS.cind(0)).pow(2);
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let production = graph.generate_3d_expression_for_integrand(
            &[],
            &canonization,
            &options,
            Some(&numerator),
        )?;
        let lu_cut = graph
            .determine_raised_esurfaces_from_expression(&production.expression)
            .raised_groups
            .into_iter()
            .find(|group| {
                group.max_occurence > 1
                    && group.esurface_ids.iter().any(|esurface_id| {
                        !production.expression.surfaces.esurface_cache[*esurface_id]
                            .external_shift
                            .is_empty()
                    })
            })
            .expect("the production CFF contains a physical raised LU surface");
        let mut cutset = CutSet::empty(graph.n_hedges());
        cutset.residue_selector.lu = Some(crate::graph::cuts::LuCutSelection {
            raised_group: lu_cut.clone(),
            cut_edge_alternatives: lu_cut
                .esurface_ids
                .iter()
                .map(|esurface_id| {
                    production.expression.surfaces.esurface_cache[*esurface_id]
                        .energies
                        .clone()
                })
                .collect(),
        });
        let contract: SuBitGraph = graph.empty_subgraph();
        let ordinary = graph.cff(
            &contract,
            &cutset,
            &OrientationPattern::default(),
            &options,
            Some(&numerator),
        )?;
        let denominators = [1, 2, 3].map(|edge| {
            let edge = EdgeIndex(edge);
            FourDDenominator {
                source_edge: edge,
                momentum: FunctionBuilder::new(GS.emr_mom)
                    .add_arg(usize::from(edge))
                    .finish(),
                mass_squared: graph.underlying[edge].particle.mass_atom().pow(2),
                full_expr: Atom::one(),
            }
        });
        let (exact, _) =
            graph.cff_from_4d_denominators(&denominators, &cutset, &options, &numerator)?;
        assert_eq!(
            exact.terms.keys().collect::<Vec<_>>(),
            ordinary.terms.keys().collect::<Vec<_>>()
        );
        for (index, ordinary_term) in &ordinary.terms {
            let mut ordinary_sum = ordinary_term
                .orientations
                .iter()
                .map(|orientation| {
                    orientation.expression.clone()
                        * numerator.replace_multiple(
                            orientation.orientation.energy_replacements_gs(&graph),
                        )
                })
                .fold(Atom::Zero, |sum, term| sum + term)
                * Atom::num(ordinary.production_prefactor_factor());
            for edge in [EdgeIndex(1), EdgeIndex(2), EdgeIndex(3)] {
                let on_shell_energy = (1..=3)
                    .fold(
                        graph.underlying[edge].particle.mass_atom().pow(2),
                        |norm_squared, spatial_index| {
                            norm_squared + GS.emr_mom(edge, GS.cind(spatial_index)).pow(2)
                        },
                    )
                    .sqrt();
                ordinary_sum = ordinary_sum.replace(GS.ose(edge)).with(on_shell_energy);
            }
            let exact_term = exact
                .terms
                .get(index)
                .expect("exact and ordinary residue keys agree");
            let exact_sum = exact_term
                .orientations
                .iter()
                .map(|orientation| {
                    Ok(orientation.expression.clone()
                        * exact_term.map_exact_source_numerator(&orientation.orientation)?)
                })
                .collect::<Result<Vec<_>>>()?
                .into_iter()
                .fold(Atom::Zero, |sum, term| sum + term)
                * Atom::num(exact.production_prefactor_factor());
            let mut difference = exact_sum - ordinary_sum;
            for edge in (0..graph.underlying.n_edges()).map(EdgeIndex) {
                for spatial_index in 1..=3 {
                    difference = difference
                        .replace(GS.emr_mom(edge, GS.cind(spatial_index)))
                        .with(Atom::Zero);
                }
                difference = difference.replace(GS.ose(edge)).with(Atom::one());
            }
            for external_edge in [EdgeIndex(0), EdgeIndex(4)] {
                difference = difference
                    .replace(GS.emr_mom(external_edge, GS.cind(0)))
                    .with(Atom::num(7));
            }
            let difference = difference.together();
            assert!(
                difference.is_zero(),
                "exact and ordinary raised LU residues differ for index {index}: {difference}"
            );
        }
        Ok(())
    }

    #[test]
    fn exact_cff_keeps_opposite_source_routing_without_a_sign_bridge() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph exact_opposite_routing {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0 lmb_id=0]
            a -> b [id=1]
        })?;
        let momentum = FunctionBuilder::new(GS.emr_mom).add_arg(0).finish();
        let mass_squared = Atom::var(GS.m_uv_expansion).pow(2);
        let denominators = [
            FourDDenominator {
                source_edge: EdgeIndex(0),
                momentum: momentum.clone(),
                mass_squared: mass_squared.clone(),
                full_expr: Atom::var(symbolica::symbol!("exact_cff_test::positive_routing")),
            },
            FourDDenominator {
                source_edge: EdgeIndex(1),
                momentum: -momentum,
                mass_squared: mass_squared.clone(),
                full_expr: Atom::var(symbolica::symbol!("exact_cff_test::negative_routing")),
            },
        ];
        let options = graph.denominator_only_cff_3d_expression_options();
        let cutset = CutSet::empty(graph.n_hedges());

        let (cff, _) =
            graph.cff_from_4d_denominators(&denominators, &cutset, &options, &Atom::one())?;
        assert_eq!(cff.terms.len(), 1);
        let on_shell_energy = (1..=3)
            .fold(mass_squared, |norm_squared, spatial_index| {
                norm_squared + GS.emr_mom(EdgeIndex(0), GS.cind(spatial_index)).pow(2)
            })
            .sqrt();
        let expected_orientation_sum =
            Atom::i() / (Atom::num(32) * Atom::var(GS.pi).pow(3) * on_shell_energy.pow(3));
        let terms = &cff
            .terms
            .values()
            .next()
            .expect("the empty cutset has one CFF term")
            .orientations;

        assert_eq!(terms.len(), 2);
        let occurrence_offset = graph.underlying.n_edges();
        let first_occurrence = EdgeIndex(occurrence_offset);
        let second_occurrence = EdgeIndex(occurrence_offset + 1);
        for term in terms {
            let orientation = &term.orientation.data.orientation;
            assert_ne!(orientation[first_occurrence], Orientation::Undirected);
            assert_ne!(orientation[second_occurrence], Orientation::Undirected);
            assert_ne!(
                orientation[second_occurrence],
                orientation[first_occurrence]
            );
            assert_eq!(orientation[EdgeIndex(0)], Orientation::Undirected);
            assert_eq!(orientation[EdgeIndex(1)], Orientation::Undirected);
        }
        assert_ne!(
            terms[0].orientation.data.orientation[first_occurrence],
            terms[1].orientation.data.orientation[first_occurrence]
        );
        let orientation_sum = terms
            .iter()
            .fold(Atom::zero(), |sum, term| sum + &term.expression)
            .replace(GS.ose(EdgeIndex(0)))
            .with(on_shell_energy);
        assert!(
            (orientation_sum - expected_orientation_sum)
                .together()
                .is_zero()
        );
        Ok(())
    }

    #[test]
    fn exact_cff_separates_uv_topology_from_the_cograph() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph exact_uv_cograph {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0 lmb_id=0]
            a -> b [id=1]
            b -> c [id=4]
            c -> d [id=2 lmb_id=1]
            c -> d [id=3]
        })?;
        let mut denominators = Vec::new();
        for edge in [EdgeIndex::from(0), EdgeIndex::from(1)] {
            denominators.push(FourDDenominator {
                source_edge: edge,
                momentum: FunctionBuilder::new(GS.emr_mom)
                    .add_arg(usize::from(edge))
                    .finish(),
                mass_squared: Atom::var(GS.m_uv_expansion).pow(2),
                full_expr: Atom::var(symbolica::symbol!("exact_cff_test::uv_full")),
            });
        }
        for edge in [EdgeIndex::from(2), EdgeIndex::from(3)] {
            denominators.push(FourDDenominator {
                source_edge: edge,
                momentum: FunctionBuilder::new(GS.emr_mom)
                    .add_arg(usize::from(edge))
                    .finish(),
                mass_squared: graph.underlying[edge].particle.mass_atom().pow(2),
                full_expr: Atom::var(symbolica::symbol!("exact_cff_test::cograph_full")),
            });
        }
        let options = graph.denominator_only_cff_3d_expression_options();
        let cutset = CutSet::empty(graph.n_hedges());

        let (cff, _) =
            graph.cff_from_4d_denominators(&denominators, &cutset, &options, &Atom::one())?;
        assert!(cff.terms.values().any(|term| !term.orientations.is_empty()));
        Ok(())
    }

    #[test]
    fn exact_cff_lu_residue_factorizes_from_quadratic_cubic_spectator() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph exact_lu_cubic_product {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0 lmb_id=0]
            a -> b [id=1]
            b -> c [id=4]
            c -> d [id=2 lmb_id=1]
            c -> d [id=3]
        })?;
        let options = graph.denominator_only_cff_3d_expression_options();
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let production = graph.generate_3d_expression_for_integrand(
            &[EdgeIndex(4)],
            &canonization,
            &options,
            Some(&Atom::one()),
        )?;
        let lu_surface = graph
            .determine_raised_esurfaces_from_expression(&production.expression)
            .raised_groups
            .into_iter()
            .find(|group| {
                group.esurface_ids.iter().any(|surface_id| {
                    production.expression.surfaces.esurface_cache[*surface_id]
                        .energies
                        .iter()
                        .copied()
                        .collect::<BTreeSet<_>>()
                        == BTreeSet::from([EdgeIndex(0), EdgeIndex(1)])
                })
            })
            .expect("the cograph bubble supplies an LU residue surface");
        let mut lu_cutset = CutSet::empty(graph.n_hedges());
        lu_cutset.residue_selector.lu = Some(crate::graph::cuts::LuCutSelection {
            raised_group: lu_surface,
            cut_edge_alternatives: vec![vec![EdgeIndex(0), EdgeIndex(1)]],
        });

        let cograph_denominators = [EdgeIndex(0), EdgeIndex(1)].map(|edge| FourDDenominator {
            source_edge: edge,
            momentum: FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(edge))
                .finish(),
            mass_squared: graph.underlying[edge].particle.mass_atom().pow(2),
            full_expr: Atom::one(),
        });
        let uv_edge = EdgeIndex(2);
        let uv_momentum = FunctionBuilder::new(GS.emr_mom)
            .add_arg(usize::from(uv_edge))
            .finish();
        let uv_mass_squared = Atom::var(GS.m_uv_expansion).pow(2);
        let uv_full_expr = GS.emr_mom(uv_edge, GS.cind(0)).pow(2)
            - (1..=3).fold(uv_mass_squared.clone(), |norm_squared, spatial_index| {
                norm_squared + GS.emr_mom(uv_edge, GS.cind(spatial_index)).pow(2)
            });
        let uv_denominators = [1, 1, -1].map(|sign| FourDDenominator {
            source_edge: uv_edge,
            momentum: Atom::num(sign) * &uv_momentum,
            mass_squared: uv_mass_squared.clone(),
            full_expr: uv_full_expr.clone(),
        });
        let quadratic_numerator = GS.emr_mom(uv_edge, GS.cind(0)).pow(2);
        let combined_denominators = cograph_denominators
            .iter()
            .cloned()
            .chain(uv_denominators.iter().cloned())
            .collect::<Vec<_>>();

        let (combined, _) = graph.cff_from_4d_denominators(
            &combined_denominators,
            &lu_cutset,
            &options,
            &quadratic_numerator,
        )?;
        let (cograph, _) = graph.cff_from_4d_denominators(
            &cograph_denominators,
            &lu_cutset,
            &options,
            &Atom::one(),
        )?;
        let (uv, _) = graph.cff_from_4d_denominators(
            &uv_denominators,
            &CutSet::empty(graph.n_hedges()),
            &options,
            &quadratic_numerator,
        )?;
        let combined_term = combined
            .terms
            .values()
            .next()
            .expect("the combined source retains the selected LU residue");
        let combined_sum = combined_term
            .orientations
            .iter()
            .map(|orientation| {
                Ok(orientation.expression.clone()
                    * combined_term.map_exact_source_numerator(&orientation.orientation)?)
            })
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .fold(Atom::Zero, |sum, term| sum + term)
            * Atom::num(combined.production_prefactor_factor());
        let cograph_term = cograph
            .terms
            .values()
            .next()
            .expect("the cograph source retains the selected LU residue");
        let cograph_sum = cograph_term
            .orientations
            .iter()
            .map(|orientation| {
                Ok(orientation.expression.clone()
                    * cograph_term.map_exact_source_numerator(&orientation.orientation)?)
            })
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .fold(Atom::Zero, |sum, term| sum + term)
            * Atom::num(cograph.production_prefactor_factor());
        let uv_term = uv
            .terms
            .values()
            .next()
            .expect("the cubic source has its uncut residue sector");
        let uv_sum = uv_term
            .orientations
            .iter()
            .map(|orientation| {
                Ok(orientation.expression.clone()
                    * uv_term.map_exact_source_numerator(&orientation.orientation)?)
            })
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .fold(Atom::Zero, |sum, term| sum + term)
            * Atom::num(uv.production_prefactor_factor());
        let difference = (combined_sum - cograph_sum * uv_sum).together();
        assert!(
            difference.is_zero(),
            "an LU residue in one exact component must factorize from a quadratic cubic spectator: {difference}"
        );
        Ok(())
    }

    #[test]
    fn exact_cff_quartic_lu_component_factorizes_from_muv_bubble() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(
            digraph exact_quartic_lu_muv_product {
                node [num=1]

                a -> b [id=0 lmb_id=0 particle=scalar_1]
                a -> b [id=1 particle=scalar_2]
                b -> c [id=4 particle=scalar_0]
                c -> d [id=2 lmb_id=1 particle=scalar_0]
                c -> d [id=3 particle=scalar_0]
            },
            "scalars"
        )?;
        let options = graph.denominator_only_cff_3d_expression_options();
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let production = graph.generate_3d_expression_for_integrand(
            &[EdgeIndex(4)],
            &canonization,
            &options,
            Some(&Atom::one()),
        )?;
        let lu_surface = graph
            .determine_raised_esurfaces_from_expression(&production.expression)
            .raised_groups
            .into_iter()
            .find(|group| {
                group.esurface_ids.iter().any(|surface_id| {
                    production.expression.surfaces.esurface_cache[*surface_id]
                        .energies
                        .iter()
                        .copied()
                        .collect::<BTreeSet<_>>()
                        == BTreeSet::from([EdgeIndex(0), EdgeIndex(1)])
                })
            })
            .expect("the nonrepeated cograph bubble supplies an LU residue surface");
        let mut lu_cutset = CutSet::empty(graph.n_hedges());
        lu_cutset.residue_selector.lu = Some(crate::graph::cuts::LuCutSelection {
            raised_group: lu_surface,
            cut_edge_alternatives: vec![vec![EdgeIndex(0), EdgeIndex(1)]],
        });

        let cograph_denominators = [EdgeIndex(0), EdgeIndex(1)].map(|edge| FourDDenominator {
            source_edge: edge,
            momentum: FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(edge))
                .finish(),
            mass_squared: graph.underlying[edge].particle.mass_atom().pow(2),
            full_expr: Atom::one(),
        });
        let muv_mass_squared = Atom::var(GS.m_uv_expansion).pow(2);
        let muv_denominators = [EdgeIndex(2), EdgeIndex(3)].map(|edge| FourDDenominator {
            source_edge: edge,
            momentum: FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(edge))
                .finish(),
            mass_squared: muv_mass_squared.clone(),
            full_expr: Atom::one(),
        });
        let quartic_numerator = GS.emr_mom(EdgeIndex(0), GS.cind(0)).pow(4);
        let combined_denominators = cograph_denominators
            .iter()
            .cloned()
            .chain(muv_denominators.iter().cloned())
            .collect::<Vec<_>>();
        let uv_edges = [EdgeIndex(2), EdgeIndex(3)];

        let cograph_source =
            GraphThreeDSource::from_exact_denominators(&graph, &cograph_denominators)?;
        assert!(
            three_dimensional_reps::repeated_groups(&cograph_source.to_three_d_parsed_graph()?)
                .is_empty(),
            "the quartic LU component must exercise a nonrepeated energy"
        );
        let combined_source = GraphThreeDSource::from_exact_denominators_in_uv_edges(
            &graph,
            &combined_denominators,
            uv_edges,
        )?;
        let (generated, _, _, _) = graph.generate_3d_expression_for_4d_term(
            &combined_source,
            &options,
            &quartic_numerator,
            None,
        )?;
        assert_eq!(
            generated.energy_factor_ownership,
            CffEnergyFactorOwnership::VariantLocal
        );
        assert_eq!(
            generated
                .energy_factor_components
                .iter()
                .filter(|component| {
                    component.ownership == CffEnergyFactorOwnership::VariantLocal
                })
                .count(),
            1,
            "only the quartic LU component owns variant-local energy factors"
        );
        assert_eq!(
            generated
                .energy_factor_components
                .iter()
                .filter(|component| {
                    component.ownership == CffEnergyFactorOwnership::GlobalSourceProduct
                })
                .count(),
            1,
            "the uncut MUV bubble remains a pure/global CFF component"
        );
        assert_eq!(
            generated
                .source_energy_degree_bounds
                .iter()
                .map(|(_, degree)| *degree)
                .collect::<Vec<_>>(),
            vec![4],
            "the nonrepeated cograph energy keeps the physical quartic bound"
        );

        let (combined, _) = graph.cff_from_4d_denominators_in_uv_edges(
            &combined_denominators,
            uv_edges,
            &lu_cutset,
            &options,
            &quartic_numerator,
            None,
        )?;
        let (cograph, _) = graph.cff_from_4d_denominators(
            &cograph_denominators,
            &lu_cutset,
            &options,
            &quartic_numerator,
        )?;
        let (muv, _) = graph.cff_from_4d_denominators_in_uv_edges(
            &muv_denominators,
            uv_edges,
            &CutSet::empty(graph.n_hedges()),
            &options,
            &Atom::one(),
            None,
        )?;
        let exact_sum = |cff: &CutCFF| -> Result<Atom> {
            assert_eq!(cff.terms.len(), 1);
            let term = cff
                .terms
                .values()
                .next()
                .expect("the selected exact source has one residue sector");
            Ok(term
                .orientations
                .iter()
                .map(|orientation| {
                    Ok(orientation.expression.clone()
                        * term.map_exact_source_numerator(&orientation.orientation)?)
                })
                .collect::<Result<Vec<_>>>()?
                .into_iter()
                .fold(Atom::Zero, |sum, term| sum + term)
                * Atom::num(cff.production_prefactor_factor()))
        };
        let difference =
            (exact_sum(&combined)? - exact_sum(&cograph)? * exact_sum(&muv)?).together();
        assert!(
            difference.is_zero(),
            "a nonrepeated quartic LU component must factorize from its pure MUV bubble: {difference}"
        );
        Ok(())
    }

    #[test]
    fn exact_original_source_matches_production_cff_after_projection() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph exact_root {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]

            incoming -> a [id=2]
            a -> b [id=0 lmb_id=0]
            a -> b [id=1]
            b -> outgoing [id=3]
        })?;
        let options = graph.denominator_only_cff_3d_expression_options();
        let cutset = CutSet::empty(graph.n_hedges());
        let contract: SuBitGraph = graph.empty_subgraph();
        let ordinary = graph.cff(
            &contract,
            &cutset,
            &OrientationPattern::default(),
            &options,
            None,
        )?;
        let denominators = [EdgeIndex::from(0), EdgeIndex::from(1)].map(|edge| FourDDenominator {
            source_edge: edge,
            momentum: FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(edge))
                .finish(),
            mass_squared: graph.underlying[edge].particle.mass_atom().pow(2),
            full_expr: Atom::one(),
        });
        let (exact, _) =
            graph.cff_from_4d_denominators(&denominators, &cutset, &options, &Atom::one())?;
        assert_eq!(
            exact.production_prefactor_factor(),
            ordinary.production_prefactor_factor(),
            "direct and exact CFF routes must carry the same production convention bridge"
        );
        let ordinary_sum = ordinary
            .terms
            .values()
            .flat_map(|term| &term.orientations)
            .fold(Atom::Zero, |sum, orientation| sum + &orientation.expression)
            * Atom::num(ordinary.production_prefactor_factor());
        let exact_sum = exact
            .terms
            .values()
            .flat_map(|term| &term.orientations)
            .fold(Atom::Zero, |sum, orientation| sum + &orientation.expression)
            * Atom::num(exact.production_prefactor_factor());
        // The ordinary route keeps its inverse-energy product in OSE notation,
        // while an exact 4D source derives square roots from each literal
        // denominator. Compare them at two exact rational on-shell points.
        // This avoids Symbolica's polynomial divider, which is not robust for
        // function-valued square-root indeterminates, without introducing any
        // floating-point tolerance.
        let rational =
            |numerator: i64, denominator: i64| Atom::num(numerator) / Atom::num(denominator);
        for (spatial_0, energy_0, spatial_1, energy_1, external_energy) in [
            (
                Atom::Zero,
                Atom::one(),
                Atom::Zero,
                Atom::one(),
                rational(1, 3),
            ),
            (
                rational(3, 4),
                rational(5, 4),
                rational(4, 3),
                rational(5, 3),
                rational(2, 5),
            ),
        ] {
            let at_point = |mut expression: Atom| {
                for (edge, spatial, energy) in [
                    (EdgeIndex(0), spatial_0.clone(), energy_0.clone()),
                    (EdgeIndex(1), spatial_1.clone(), energy_1.clone()),
                ] {
                    expression = expression
                        .replace(GS.ose(edge))
                        .with(energy)
                        .replace(GS.emr_mom(edge, GS.cind(1)))
                        .with(spatial);
                    for spatial_index in 2..=3 {
                        expression = expression
                            .replace(GS.emr_mom(edge, GS.cind(spatial_index)))
                            .with(Atom::Zero);
                    }
                }
                expression
                    .replace(GS.emr_mom(EdgeIndex(3), GS.cind(0)))
                    .with(external_energy.clone())
                    .expand()
            };
            let difference =
                (at_point(exact_sum.clone()) - at_point(ordinary_sum.clone())).expand();
            assert!(
                difference.is_zero(),
                "source-local canonical orientations must reproduce the ordinary CFF at exact on-shell points: {difference}"
            );
        }
        Ok(())
    }

    #[test]
    fn exact_two_loop_source_matches_production_prefactor_bridge() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph exact_two_loop_prefactor {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0 lmb_id=0]
            a -> b [id=1 lmb_id=1]
            b -> a [id=2]
        })?;
        let options = graph.denominator_only_cff_3d_expression_options();
        let cutset = CutSet::empty(graph.n_hedges());
        let contract: SuBitGraph = graph.empty_subgraph();
        let ordinary = graph.cff(
            &contract,
            &cutset,
            &OrientationPattern::default(),
            &options,
            None,
        )?;
        let denominators = [0, 1, 2].map(|edge| {
            let edge = EdgeIndex::from(edge);
            FourDDenominator {
                source_edge: edge,
                momentum: FunctionBuilder::new(GS.emr_mom)
                    .add_arg(usize::from(edge))
                    .finish(),
                mass_squared: graph.underlying[edge].particle.mass_atom().pow(2),
                full_expr: Atom::one(),
            }
        });
        let (exact, _) =
            graph.cff_from_4d_denominators(&denominators, &cutset, &options, &Atom::one())?;

        assert_eq!(ordinary.production_prefactor_factor(), -1);
        assert_eq!(
            exact.production_prefactor_factor(),
            ordinary.production_prefactor_factor(),
            "direct and exact two-loop CFF routes must carry the same production convention bridge"
        );
        Ok(())
    }

    #[test]
    fn stored_production_cff_uses_typed_energy_factor_ownership() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph typed_energy_factor_ownership {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]

            incoming -> a [id=2]
            a -> b [id=0 lmb_id=0]
            a -> b [id=1]
            b -> outgoing [id=3]
        })?;
        let options = graph.denominator_only_cff_3d_expression_options();
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let contract_subgraph = graph.tree_edges.subtract(&graph.initial_state_cut);
        let contract_edges = graph.paired_edges(&contract_subgraph);
        let generated = graph.generate_3d_expression_for_integrand(
            &contract_edges,
            &canonization,
            &options,
            None,
        )?;
        assert_eq!(
            generated.energy_factor_ownership,
            CffEnergyFactorOwnership::GlobalSourceProduct
        );
        let mut variant_local = generated.clone();
        variant_local.energy_factor_ownership = CffEnergyFactorOwnership::VariantLocal;
        let cutset = CutSet::empty(graph.n_hedges());
        let global = graph.cff_from_production_expression(
            &generated,
            &cutset,
            &OrientationPattern::default(),
        )?;
        let local = graph.cff_from_production_expression(
            &variant_local,
            &cutset,
            &OrientationPattern::default(),
        )?;
        assert_eq!(generated.core_global_prefactor_sign.factor(), -1);
        assert_eq!(
            global.production_prefactor_factor(),
            generated.core_global_prefactor_sign.factor()
        );
        assert_eq!(
            local.production_prefactor_factor(),
            global.production_prefactor_factor()
        );
        let graph_without_is_cut = graph
            .underlying
            .full_filter()
            .subtract(&graph.initial_state_cut.left)
            .subtract(&graph.initial_state_cut.right);
        let inverse_energy_product =
            get_cff_inverse_energy_product_impl(&graph, &graph_without_is_cut, &contract_edges);
        let global_terms = &global
            .terms
            .values()
            .next()
            .expect("the empty cutset has one global CFF term")
            .orientations;
        let local_terms = &local
            .terms
            .values()
            .next()
            .expect("the empty cutset has one local CFF term")
            .orientations;

        assert_eq!(global_terms.len(), local_terms.len());
        for (orientation_id, (global, local)) in global_terms.iter().zip(local_terms).enumerate() {
            assert_eq!(
                global.production_orientation_id,
                Some(OrientationID(orientation_id))
            );
            assert_eq!(
                local.production_orientation_id,
                global.production_orientation_id
            );
            assert!(
                (global.expression.clone() - &local.expression * &inverse_energy_product)
                    .together()
                    .is_zero()
            );
        }
        Ok(())
    }
}
