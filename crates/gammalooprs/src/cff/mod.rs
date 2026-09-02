use std::{collections::BTreeMap, fmt::Display, sync::Arc};

use bincode_trait_derive::{Decode, Encode};
use linnet::half_edge::{
    involution::{EdgeIndex, Hedge},
    subgraph::{SubGraphLike, SubSetLike, SubSetOps},
};
use serde::{Deserialize, Serialize};
use symbolica::atom::{Atom, AtomCore};

use crate::{
    cff::{
        expression::{
            GammaLoopOrientationExpression, OrientationExpression, OrientationID,
            OrientationSelector, ThreeDExpression,
            normalize_cut_edge_support_with_raised_edge_groups,
            normalize_three_d_expression_cut_support_with_raised_edge_groups,
        },
        orientations::GraphOrientation,
        surface::GammaLoopSurfaceCache,
    },
    graph::{
        ExactUvSubLmbFrame, FeynmanGraph, FourDDenominator, Graph, GraphThreeDSource,
        LoopMomentumBasis, cuts::CutSet, get_cff_inverse_energy_product_impl,
        three_d_source::ExactSourceEnergyMapper,
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
    #[cfg(test)]
    pub(crate) fn map_exact_source_atom(
        &self,
        orientation: &OrientationExpression,
        atom: &Atom,
    ) -> Result<Atom> {
        let planned = self
            .exact_source_numerator
            .as_ref()
            .ok_or_else(|| eyre::eyre!("ordinary CFF term has no exact-source numerator map"))?;
        planned.mapper.map_numerator(
            &orientation.loop_energy_map,
            &orientation.edge_energy_map,
            atom,
        )
    }

    #[cfg(test)]
    pub(crate) fn map_exact_source_physical_loop_lift_energies(
        &self,
        orientation: &OrientationExpression,
        physical_edges: impl IntoIterator<Item = EdgeIndex>,
    ) -> Result<Vec<(EdgeIndex, Atom)>> {
        let planned = self
            .exact_source_numerator
            .as_ref()
            .ok_or_else(|| eyre::eyre!("ordinary CFF term has no exact-source numerator map"))?;
        planned.mapper.map_physical_owner_loop_lift_energies(
            &orientation.loop_energy_map,
            &orientation.edge_energy_map,
            physical_edges,
        )
    }

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
                let selector = term.production_orientation_id.map_or_else(
                    || term.orientation.data.orientation.orientation_thetas(),
                    OrientationID::atom,
                );
                term.expression.clone() * selector
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
    residues: Vec<(CutCFFIndex, ThreeDExpression<OrientationID>)>,
    axis: CutCffResidueAxis,
    mut select: F,
) -> Vec<(CutCFFIndex, ThreeDExpression<OrientationID>)>
where
    F: FnMut(ThreeDExpression<OrientationID>) -> Vec<ThreeDExpression<OrientationID>>,
{
    residues
        .into_iter()
        .flat_map(|(index, expression)| {
            select(expression)
                .into_iter()
                .enumerate()
                .map(move |(i, residue)| {
                    let mut new_index = index;
                    axis.set_order(&mut new_index, i + 1);
                    // One causal coordinate has been consumed regardless of
                    // which raised-order branch `i` labels. Its energy-factor,
                    // contact, and residue signs remain owned by the selected
                    // CFF variant; the index records only the residue order.
                    (new_index, residue)
                })
        })
        .collect()
}

fn select_indexed_cff_residues(
    cff: ThreeDExpression<OrientationID>,
    cutset: &CutSet,
) -> Result<Vec<(CutCFFIndex, ThreeDExpression<OrientationID>)>> {
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
    // `CutCFFIndex`. Energy-factor ownership stays with the generated CFF
    // variants and is deliberately independent of this indexing operation.
    let mut residues = vec![(CutCFFIndex::new_all_none(), cff)];

    if let Some(right_threshold) = cutset.residue_selector.right_th_cut.as_ref() {
        residues = apply_indexed_residue_selection(
            residues,
            CutCffResidueAxis::RightThreshold,
            |expression| expression.select_esurface_residue(right_threshold),
        );
    }

    if let Some(left_threshold) = cutset.residue_selector.left_th_cut.as_ref() {
        residues = apply_indexed_residue_selection(
            residues,
            CutCffResidueAxis::LeftThreshold,
            // Powered poles are completed in the canonical causal basis, so
            // left- and right-threshold selections use the same positive-
            // energy Cutkosky convention as ordinary CFF terms.
            |expression| expression.select_esurface_residue(left_threshold),
        );
    }

    if let Some(lu) = cutset.residue_selector.lu.as_ref() {
        residues =
            apply_indexed_residue_selection(residues, CutCffResidueAxis::LuCut, |expression| {
                // Physical support selection is independent of the 3D
                // representation. CFF then consumes the selected surface in
                // GammaLoop's positive-energy Cutkosky convention; future LTD
                // signs and local-series coordinates remain LTD-owned.
                expression
                    .restrict_to_cut_alternatives(&lu.cut_edge_alternatives)
                    .select_esurface_residue(&lu.raised_group)
            });
    }

    Ok(residues)
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
            CffEnergyFactorOwnership::VariantLocal => Atom::num(1),
        };

        let mut terms = BTreeMap::new();
        for (cut_cff_index, expr) in residues {
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

    #[cfg(test)]
    pub(crate) fn cff_from_4d_denominators_in_uv_edges(
        &mut self,
        denominators: &[FourDDenominator],
        uv_edges: impl IntoIterator<Item = EdgeIndex>,
        cutset: &CutSet,
        options: &Generate3DExpressionOptions,
        analysis_numerator: &Atom,
        generation_cache: Option<&mut generation::ExactCffGenerationCache>,
    ) -> Result<(CutCFF, linnet::half_edge::subgraph::SuBitGraph)> {
        self.cff_from_4d_denominators_in_uv_edges_and_boundaries(
            denominators,
            uv_edges,
            [],
            cutset,
            options,
            analysis_numerator,
            generation_cache,
        )
    }

    /// Generate an exact CFF while retaining the crown of a non-vacuum UV
    /// source. Vacuum Taylor terms use the boundary-free wrapper above.
    #[cfg(test)]
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn cff_from_4d_denominators_in_uv_edges_and_boundaries(
        &mut self,
        denominators: &[FourDDenominator],
        uv_edges: impl IntoIterator<Item = EdgeIndex>,
        uv_boundary_hedges: impl IntoIterator<Item = Hedge>,
        cutset: &CutSet,
        options: &Generate3DExpressionOptions,
        analysis_numerator: &Atom,
        generation_cache: Option<&mut generation::ExactCffGenerationCache>,
    ) -> Result<(CutCFF, linnet::half_edge::subgraph::SuBitGraph)> {
        self.cff_from_4d_denominators_in_uv_coordinates(
            denominators,
            uv_edges,
            uv_boundary_hedges,
            None,
            None,
            cutset,
            options,
            analysis_numerator,
            generation_cache,
        )
    }

    /// Generate a proper-subgraph exact CFF in its sub-LMB coordinates. Crown
    /// momenta are external in this source and remain available to the later
    /// outer CFF instead of being treated as inactive child loop energies.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn cff_from_4d_denominators_in_uv_sub_lmb(
        &mut self,
        denominators: &[FourDDenominator],
        uv_edges: impl IntoIterator<Item = EdgeIndex>,
        uv_boundary_hedges: impl IntoIterator<Item = Hedge>,
        sub_lmb: &LoopMomentumBasis,
        frame: ExactUvSubLmbFrame,
        cutset: &CutSet,
        options: &Generate3DExpressionOptions,
        analysis_numerator: &Atom,
        generation_cache: Option<&mut generation::ExactCffGenerationCache>,
    ) -> Result<(CutCFF, linnet::half_edge::subgraph::SuBitGraph)> {
        self.cff_from_4d_denominators_in_uv_coordinates(
            denominators,
            uv_edges,
            uv_boundary_hedges,
            Some(sub_lmb),
            Some(frame),
            cutset,
            options,
            analysis_numerator,
            generation_cache,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn cff_from_4d_denominators_in_uv_coordinates(
        &mut self,
        denominators: &[FourDDenominator],
        uv_edges: impl IntoIterator<Item = EdgeIndex>,
        uv_boundary_hedges: impl IntoIterator<Item = Hedge>,
        sub_lmb: Option<&LoopMomentumBasis>,
        sub_lmb_frame: Option<ExactUvSubLmbFrame>,
        cutset: &CutSet,
        options: &Generate3DExpressionOptions,
        analysis_numerator: &Atom,
        mut generation_cache: Option<&mut generation::ExactCffGenerationCache>,
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
            physical_cut_support_edges,
        ) = {
            let source = if let Some(sub_lmb) = sub_lmb {
                GraphThreeDSource::from_exact_denominators_in_uv_sub_lmb(
                    self,
                    denominators,
                    uv_edges,
                    uv_boundary_hedges,
                    sub_lmb,
                    sub_lmb_frame.expect("a sub-LMB exact CFF has an incidence frame"),
                )?
            } else {
                GraphThreeDSource::from_exact_denominators_in_uv_edges_and_boundaries(
                    self,
                    denominators,
                    uv_edges,
                    uv_boundary_hedges,
                )?
            };
            let (
                generated,
                exact_source_energy_mapper,
                energy_assignment,
                energy_degree_bound_report,
            ) = self.generate_3d_expression_for_4d_term(
                &source,
                options,
                analysis_numerator,
                generation_cache.as_deref_mut(),
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
                source
                    .physical_cut_support_edge_index_map()
                    .expect("exact 4D source has a physical cut-support projection"),
            )
        };
        let (generated, surface_ownership) =
            self.convert_4d_expression_surfaces(generated, &physical_surfaces)?;
        let production_prefactor_bridge = generated.core_global_prefactor_sign;
        let energy_factor_ownership = generated.energy_factor_ownership;
        // The generated expression retains its own shared-core contour sign.
        // Consume that exact typed convention rather than inferring a second
        // sign from a denominator-only or provenance-collapsed source.
        crate::debug_tags!(#generation, #cff, #inspect;
            denominator_count = denominators.len(),
            production_prefactor_bridge = ?production_prefactor_bridge,
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
        let physical_edges = |edge: linnet::half_edge::involution::EdgeIndex| {
            let edge_id = usize::from(edge);
            if edge_id < physical_energy_edges.orientation_edge_count {
                return Ok(vec![edge]);
            }
            physical_cut_support_edges
                .get(&edge_id)
                .cloned()
                .ok_or_else(|| {
                    eyre::eyre!("exact CFF cut support contains unmapped occurrence edge {edge_id}")
                })
        };
        let raised_edge_groups = self.get_raised_edge_groups();
        let retain_physical_support_with_raised_representatives = |support: &mut Vec<
            linnet::half_edge::involution::EdgeIndex,
        >| {
            let representatives =
                normalize_cut_edge_support_with_raised_edge_groups(support, &raised_edge_groups);
            support.extend(representatives);
            support.sort_unstable();
            support.dedup();
        };
        for orientation in cff.orientations.iter_mut() {
            for variant in &mut orientation.variants {
                let mut denominator_edges = variant
                    .denominator_edges
                    .iter()
                    .copied()
                    .map(&physical_edges)
                    .collect::<Result<Vec<_>>>()?
                    .into_iter()
                    .flatten()
                    .collect::<Vec<_>>();
                retain_physical_support_with_raised_representatives(&mut denominator_edges);
                variant.denominator_edges = denominator_edges;
                variant.denominator_edge_support_signs =
                    std::mem::take(&mut variant.denominator_edge_support_signs)
                        .into_iter()
                        .try_fold(BTreeMap::new(), |mut mapped, (support, sign)| {
                            let mut support = support
                                .into_iter()
                                .map(&physical_edges)
                                .collect::<Result<Vec<_>>>()?
                                .into_iter()
                                .flatten()
                                .collect::<Vec<_>>();
                            retain_physical_support_with_raised_representatives(&mut support);
                            *mapped.entry(support).or_insert(1) *= sign;
                            Ok::<_, color_eyre::Report>(mapped)
                        })?;
            }
        }
        let residues = select_indexed_cff_residues(cff, cutset)?;
        let cff_phase = (-Atom::i()).pow(cff_loop_number as i64);
        let cff_normalization = cff_phase / (Atom::var(GS.pi) * 2).pow(3 * cff_loop_number as i64);
        let mut terms = BTreeMap::new();
        for (cut_cff_index, expr) in residues {
            let cff_energy_factor = match energy_factor_ownership {
                CffEnergyFactorOwnership::GlobalSourceProduct => inverse_energy_product.clone(),
                CffEnergyFactorOwnership::VariantLocal => Atom::num(1),
            };
            crate::debug_tags!(#generation, #cff, #inspect;
                cut_index = ?cut_cff_index,
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
            CffEnergyFactorOwnership::VariantLocal => Atom::num(1),
        };
        crate::debug_tags!(#cff, #trace;
            stage = "graph_cff_normalization",
            graph = %self.name,
            cff_loop_number = cff_loop_number,
            production_prefactor_bridge = ?production_prefactor_bridge,
            log.cff_normalization = cff_normalization,
            "Graph CFF normalization"
        );

        let mut terms = BTreeMap::new();

        for (cut_cff_index, expr) in residues {
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
    use crate::{
        cff::surface::GammaLoopLinearEnergyExpr,
        dot,
        graph::{LMBext, parse::IntoGraph},
        initialisation::test_initialise,
        integrands::{
            evaluation::EvaluationMetaData,
            process::{
                GenericEvaluatorFloat,
                evaluators::{EvaluatorStack, SingleOrAllOrientations},
            },
        },
        momentum::{
            ThreeMomentum,
            sample::{BareMomentumSample, ExternalFourMomenta, LoopMomenta, MomentumSample},
        },
        processes::EvaluatorSettings,
        settings::RuntimeSettings,
        utils::{ArbPrec, F, FloatLike, W_, cut_energy},
        uv::{
            UltravioletGraph,
            approx::{OrientationProjection, local_3d::Localizer, local_4d::Full4dCts},
        },
    };
    use itertools::Itertools;
    use linnet::half_edge::involution::{EdgeIndex, EdgeVec, Orientation};
    use linnet::half_edge::subgraph::{
        InternalSubGraph, ModifySubSet, SuBitGraph, SubSetLike, SubSetOps, subset::SubSet,
    };
    use spenso::algebra::{algebraic_traits::IsZero, complex::Complex};
    use symbolica::{
        atom::FunctionBuilder,
        domains::{
            float::{Complex as SymComplex, Real},
            integer::IntegerRing,
            rational::{Fraction, Rational},
        },
        evaluate::ExpressionEvaluator,
        function,
        id::Replacement,
    };
    use three_dimensional_reps::{
        LinearEnergyExpr, OrientationData, OrientationID, ThreeDGraphSource, repeated_groups,
    };
    use typed_index_collections::TiVec;

    #[test]
    fn cff_term_selects_duplicate_physical_orientations_by_production_map_key() {
        let orientation = || OrientationExpression {
            data: OrientationData::new(EdgeVec::from_iter([Orientation::Default])),
            loop_energy_map: Vec::new(),
            edge_energy_map: Vec::new(),
            variants: Vec::new(),
        };
        let keyed = CFFTerm {
            orientations: vec![
                CFFOrientationTerm {
                    expression: Atom::num(2),
                    orientation: orientation(),
                    production_orientation_id: Some(OrientationID(0)),
                },
                CFFOrientationTerm {
                    expression: Atom::num(3),
                    orientation: orientation(),
                    production_orientation_id: Some(OrientationID(1)),
                },
            ],
            exact_source_numerator: None,
        }
        .expression_with_selectors();

        assert_eq!(OrientationID(0).select(&keyed).expand(), Atom::num(2));
        assert_eq!(OrientationID(1).select(&keyed).expand(), Atom::num(3));
        assert_eq!(
            keyed
                .replace(function!(OrientationID::symbol(), W_.a_))
                .with(Atom::one())
                .expand(),
            Atom::num(5),
            "explicit summation must retain every distinct map key even when their physical directions coincide",
        );

        let source_local = CFFTerm {
            orientations: vec![CFFOrientationTerm {
                expression: Atom::num(7),
                orientation: orientation(),
                production_orientation_id: None,
            }],
            exact_source_numerator: None,
        }
        .expression_with_selectors();
        assert_eq!(
            source_local,
            Atom::num(7) * GS.sign_theta(GS.sign(EdgeIndex(0))),
            "a source-local term without a production map key must retain its physical selector",
        );
    }

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
    fn exact_child_sub_lmb_retains_parent_loop_crown_shift() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph exact_child_sub_lmb {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]

            incoming -> v3 [id=0]
            v0 -> v1 [id=1 lmb_id=0]
            v0 -> v1 [id=2]
            v3 -> v0 [id=3 lmb_id=1]
            v1 -> v2 [id=4]
            v2 -> v3 [id=5]
            v2 -> outgoing [id=6]
        })?;
        let uv_edges = [EdgeIndex(1), EdgeIndex(2)];
        let uv_filter = graph
            .get_edge_subgraph(uv_edges[0])
            .union(&graph.get_edge_subgraph(uv_edges[1]));
        let uv_subgraph = InternalSubGraph::cleaned_filter_optimist(uv_filter, graph.as_ref());
        let crown = graph.dummy_stripped_external_flows_of(&uv_subgraph);
        let boundary_hedges = crown.included_iter().collect::<Vec<_>>();
        let sub_lmb =
            graph.try_compatible_sub_lmb(&uv_subgraph, crown, &graph.loop_momentum_basis)?;
        assert_eq!(sub_lmb.loop_edges.len(), 1);
        assert!(sub_lmb.ext_edges.contains(&EdgeIndex(4)));

        let denominators = uv_edges.map(|edge| FourDDenominator {
            source_edge: edge,
            momentum: FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(edge))
                .finish(),
            mass_squared: graph.underlying[edge].particle.mass_atom().pow(2),
            full_expr: Atom::one(),
        });
        assert!(
            GraphThreeDSource::from_exact_denominators_in_uv_sub_lmb(
                &graph,
                &denominators,
                uv_edges,
                boundary_hedges.iter().copied().take(1),
                &sub_lmb,
                ExactUvSubLmbFrame::RetainedPhysicalCrown,
            )
            .is_err(),
            "a proper child source must not omit part of its non-dummy crown",
        );
        {
            let source = GraphThreeDSource::from_exact_denominators_in_uv_sub_lmb(
                &graph,
                &denominators,
                uv_edges,
                boundary_hedges.iter().copied(),
                &sub_lmb,
                ExactUvSubLmbFrame::RetainedPhysicalCrown,
            )?;
            let parsed = source.to_three_d_parsed_graph()?;
            assert_eq!(parsed.loop_names.len(), 1);
            let local_to_occurrence = source
                .energy_edge_index_map(&parsed)
                .expect("the exact child source exposes occurrence energies")
                .internal;
            let occurrence_to_owner = source
                .physical_energy_edge_index_map()
                .expect("the exact child source exposes physical owners")
                .internal;
            let shifted_local_edge = local_to_occurrence
                .iter()
                .find_map(|(local, occurrence)| {
                    (occurrence_to_owner.get(occurrence) == Some(&usize::from(uv_edges[1])))
                        .then_some(*local)
                })
                .expect("the second child denominator retains its source occurrence");
            let shifted_signature = &parsed.internal_edges[shifted_local_edge].signature;
            assert!(
                shifted_signature
                    .external_signature
                    .iter()
                    .any(|coefficient| *coefficient != 0),
                "D2 must retain the parent-loop crown momentum as a child external shift",
            );
        }

        let options = graph.denominator_only_cff_3d_expression_options();
        let cutset = CutSet::empty(graph.n_hedges());
        let (exact, _) = graph.cff_from_4d_denominators_in_uv_sub_lmb(
            &denominators,
            uv_edges,
            boundary_hedges,
            &sub_lmb,
            ExactUvSubLmbFrame::RetainedPhysicalCrown,
            &cutset,
            &options,
            &Atom::one(),
            None,
        )?;
        let exact_term = exact
            .terms
            .get(&CutCFFIndex::new_all_none())
            .expect("the exact child has one uncut sector");
        let exact_sum = exact_term.orientations.iter().try_fold(
            Atom::Zero,
            |sum, orientation| -> Result<Atom> {
                Ok(sum
                    + &orientation.expression
                        * exact_term.map_exact_source_numerator(&orientation.orientation)?)
            },
        )? * Atom::num(exact.production_prefactor_factor());

        let mut child: Graph = dot!(digraph ordinary_child {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]

            incoming -> a [id=2]
            a -> b [id=0 lmb_id=0]
            a -> b [id=1]
            b -> outgoing [id=3]
        })?;
        let child_options = child.denominator_only_cff_3d_expression_options();
        let child_cutset = CutSet::empty(child.n_hedges());
        let child_contract = child.empty_subgraph::<SuBitGraph>();
        let ordinary = child.cff(
            &child_contract,
            &child_cutset,
            &OrientationPattern::default(),
            &child_options,
            None,
        )?;
        let ordinary_term = ordinary
            .terms
            .get(&CutCFFIndex::new_all_none())
            .expect("the standalone child has one uncut sector");
        let mut ordinary_sum = ordinary_term
            .orientations
            .iter()
            .fold(Atom::Zero, |sum, orientation| sum + &orientation.expression)
            * Atom::num(ordinary.production_prefactor_factor());
        ordinary_sum = ordinary_sum.replace_multiple(
            [(0, 1), (1, 2), (2, 3), (3, 4)]
                .into_iter()
                .flat_map(|(child_edge, parent_edge)| {
                    [
                        Replacement::new(
                            function!(GS.ose, child_edge, W_.x___).to_pattern(),
                            function!(GS.ose, parent_edge, W_.x___).to_pattern(),
                        ),
                        Replacement::new(
                            function!(GS.emr_mom, child_edge, W_.x___).to_pattern(),
                            function!(GS.emr_mom, parent_edge, W_.x___).to_pattern(),
                        ),
                    ]
                })
                .collect::<Vec<_>>(),
        );
        assert_eq!(
            exact_sum, ordinary_sum,
            "the parent-source child CFF must equal its standalone ordinary representation"
        );
        Ok(())
    }

    #[test]
    fn exact_child_sub_lmb_keeps_remote_cut_and_crown_emrs_factorized() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(
            digraph exact_child_lu_aliases {
                num = 1
                edge [particle="scalar_1" num=1]
                node [num=1]

                ext0 [style=invis is_cut=0]
                v2 -> ext0 [id=0]
                ext0 -> v3
                v0 -> v1 [id=1 lmb_id=0]
                v0 -> v1 [id=2]
                v3 -> v0 [id=3 lmb_id=1]
                v1 -> v2 [id=4]
                v2 -> v3 [id=5]
            },
            "scalars"
        )?;
        let uv_edges = [EdgeIndex(1), EdgeIndex(2)];
        let uv_filter = graph
            .get_edge_subgraph(uv_edges[0])
            .union(&graph.get_edge_subgraph(uv_edges[1]));
        let uv_subgraph = InternalSubGraph::cleaned_filter_optimist(uv_filter, graph.as_ref());
        let crown = graph.dummy_stripped_external_flows_of(&uv_subgraph);
        let boundary_hedges = crown.included_iter().collect::<Vec<_>>();
        assert!(
            boundary_hedges
                .iter()
                .all(|hedge| graph.underlying[hedge] != EdgeIndex(0)),
            "the remote Cutkosky carrier is not part of the physical child crown",
        );
        let sub_lmb =
            graph.try_compatible_sub_lmb(&uv_subgraph, crown, &graph.loop_momentum_basis)?;
        assert!(
            sub_lmb.edge_signatures[EdgeIndex(0)]
                .external
                .iter()
                .all(|coefficient| *coefficient == crate::momentum::SignOrZero::Zero),
            "the remote cut must have no child external coordinate",
        );
        let denominators = uv_edges.map(|edge| FourDDenominator {
            source_edge: edge,
            momentum: FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(edge))
                .finish(),
            mass_squared: graph.underlying[edge].particle.mass_atom().pow(2),
            full_expr: Atom::one(),
        });
        {
            let source = GraphThreeDSource::from_exact_denominators_in_uv_sub_lmb(
                &graph,
                &denominators,
                uv_edges,
                boundary_hedges.iter().copied(),
                &sub_lmb,
                ExactUvSubLmbFrame::RetainedPhysicalCrown,
            )?;
            let parsed = source.to_three_d_parsed_graph()?;
            assert_eq!(parsed.loop_names.len(), 1);
            assert!(
                parsed.initial_state_cut_edges.is_empty(),
                "a factorized UV child must not import its parent's Cutkosky carrier",
            );

            let options = graph.denominator_only_cff_3d_expression_options();
            let shell_energy = GS.emr_mom(EdgeIndex(1), GS.cind(0));
            let crown_energy = GS.emr_mom(EdgeIndex(3), GS.cind(0));
            let analysis_numerator = &shell_energy * &crown_energy;
            let (generated, mapper, plan, report) = graph.generate_3d_expression_for_4d_term(
                &source,
                &options,
                &analysis_numerator,
                None,
            )?;
            assert_eq!(report.physical_parent_bounds, vec![(1, 1)]);
            assert_eq!(plan.energy_degree_bounds().len(), 1);
            assert_eq!(plan.energy_degree_bounds()[0].1, 1);
            assert_eq!(
                report.assigned_cff_source_bounds,
                plan.energy_degree_bounds()
            );
            let orientation = generated
                .expression
                .orientations
                .first()
                .expect("the child bubble has a causal orientation");
            let mapped_shell = mapper.map_numerator(
                &orientation.loop_energy_map,
                &orientation.edge_energy_map,
                &shell_energy,
            )?;
            assert_ne!(mapped_shell, shell_energy);
            assert_eq!(
                mapper.map_numerator(
                    &orientation.loop_energy_map,
                    &orientation.edge_energy_map,
                    &analysis_numerator,
                )?,
                mapped_shell * crown_energy,
                "the shell energy is CFF-mapped while the crown energy stays factorized",
            );
            let factorized_parent =
                GS.emr_mom(EdgeIndex(0), GS.cind(0)) + GS.emr_mom(EdgeIndex(3), GS.cind(0));
            assert_eq!(
                mapper.map_numerator(
                    &orientation.loop_energy_map,
                    &orientation.edge_energy_map,
                    &factorized_parent,
                )?,
                factorized_parent,
                "remote cut and crown EMRs must remain T-inert for the outer graph",
            );
        }

        let options = graph.denominator_only_cff_3d_expression_options();
        let cutset = CutSet::empty(graph.n_hedges());
        graph.cff_from_4d_denominators_in_uv_sub_lmb(
            &denominators,
            uv_edges,
            boundary_hedges,
            &sub_lmb,
            ExactUvSubLmbFrame::RetainedPhysicalCrown,
            &cutset,
            &options,
            &Atom::one(),
            None,
        )?;
        Ok(())
    }

    #[test]
    fn exact_child_signed_maps_keep_owner_invariant_production_hosts() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph exact_child_orientation_projection {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]

            incoming -> v3 [id=0]
            v0 -> v1 [id=1 lmb_id=0]
            v0 -> v1 [id=2]
            v3 -> v0 [id=3 lmb_id=1]
            v1 -> v2 [id=4]
            v2 -> v3 [id=5]
            v2 -> outgoing [id=6]
        })?;
        let uv_edges = [EdgeIndex(1), EdgeIndex(2)];
        let uv_filter = graph
            .get_edge_subgraph(uv_edges[0])
            .union(&graph.get_edge_subgraph(uv_edges[1]));
        let uv_subgraph = InternalSubGraph::cleaned_filter_optimist(uv_filter, graph.as_ref());
        let crown = graph.dummy_stripped_external_flows_of(&uv_subgraph);
        let boundary_hedges = crown.included_iter().collect::<Vec<_>>();
        let sub_lmb =
            graph.try_compatible_sub_lmb(&uv_subgraph, crown, &graph.loop_momentum_basis)?;
        let options = graph.denominator_only_cff_3d_expression_options();
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let production =
            graph.generate_3d_expression_for_integrand(&[], &canonization, &options, None)?;
        let pattern = OrientationPattern::default();
        let cutset = CutSet::empty(graph.n_hedges());
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact_expression(&production, &options, &pattern, true),
        );
        let source_denominator = |edge| FourDDenominator {
            source_edge: edge,
            momentum: FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(edge))
                .finish(),
            mass_squared: graph.underlying[edge].particle.mass_atom().pow(2),
            full_expr: Atom::one(),
        };
        let denominators = uv_edges.map(source_denominator);
        let mut opposite_denominators = denominators.clone();
        opposite_denominators[0].momentum = -opposite_denominators[0].momentum.clone();
        let physical_host_classes =
            |denominators: &[FourDDenominator]| -> Result<BTreeSet<Vec<OrientationID>>> {
                let source = GraphThreeDSource::from_exact_denominators_in_uv_sub_lmb(
                    &graph,
                    denominators,
                    uv_edges,
                    boundary_hedges.iter().copied(),
                    &sub_lmb,
                    ExactUvSubLmbFrame::RetainedPhysicalCrown,
                )?;
                let (child, _, _, _) = graph.generate_3d_expression_for_4d_term(
                    &source,
                    &options,
                    &Atom::one(),
                    None,
                )?;
                child
                    .expression
                    .orientations
                    .iter()
                    .map(|child_orientation| {
                        let candidates = localizer.source_selector_representatives(
                            &graph,
                            &sub_lmb,
                            child_orientation,
                            &source.contract_subgraph(),
                            None,
                        )?;
                        assert!(
                            !candidates.is_empty(),
                            "every exact child residue must extend to a production map using only its physical graph-edge prefix; denominators={denominators:?}, child={:?}",
                            child_orientation.data,
                        );
                        Ok(candidates)
                    })
                    .collect()
            };

        let ordinary_hosts = physical_host_classes(&denominators)?;
        assert_eq!(
            physical_host_classes(&opposite_denominators)?,
            ordinary_hosts,
            "D(Q) and D(-Q) must induce identical physical orientation constraints",
        );
        Ok(())
    }

    #[test]
    fn exact_unexpanded_cff_matches_ordinary_for_cross_loop_factorized_energy() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(
            digraph exact_unexpanded_cross_loop_factorized_energy {
                num = 1
                edge [particle="scalar_1" num=1]
                node [num=1]

                ext0 [style=invis is_cut=0]
                v2 -> ext0 [id=0]
                ext0 -> v3
                v0 -> v1 [id=1 lmb_id=0]
                v0 -> v1 [id=2]
                v3 -> v0 [id=3 lmb_id=1]
                v1 -> v2 [id=4]
                v2 -> v3 [id=5]
            },
            "scalars"
        )?;
        let numerator = GS.emr_mom(EdgeIndex(1), GS.cind(0))
            * GS.emr_mom(EdgeIndex(3), GS.cind(0)).pow(2)
            * GS.emr_mom(EdgeIndex(4), GS.cind(0))
            * GS.emr_mom(EdgeIndex(5), GS.cind(0));
        let options = graph.denominator_only_cff_3d_expression_options();
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let production = graph.generate_3d_expression_for_integrand(
            &[],
            &canonization,
            &options,
            Some(&numerator),
        )?;
        let cutset = CutSet::empty(graph.n_hedges());
        let ordinary = graph.cff_from_production_expression(
            &production,
            &cutset,
            &OrientationPattern::default(),
        )?;
        let cograph = graph.full_filter().subtract(&graph.initial_state_cut);
        let denominators = Full4dCts::from_coefficient(&Atom::one(), &graph, &cograph)
            .terms()?
            .into_iter()
            .next()
            .expect("the unexpanded graph has one exact denominator term")
            .denominators;
        let (exact, _) =
            graph.cff_from_4d_denominators(&denominators, &cutset, &options, &numerator)?;

        assert_eq!(
            ordinary.energy_degree_bound_report.physical_parent_bounds,
            exact.energy_degree_bound_report.physical_parent_bounds
        );
        assert_eq!(
            ordinary
                .energy_degree_bound_report
                .assigned_cff_source_bounds,
            exact.energy_degree_bound_report.assigned_cff_source_bounds,
            "a one-occurrence exact source must preserve physical energy IDs"
        );
        let index = CutCFFIndex::new_all_none();
        let ordinary_term = ordinary
            .terms
            .get(&index)
            .expect("the ordinary uncut sector exists");
        let exact_term = exact
            .terms
            .get(&index)
            .expect("the exact uncut sector exists");
        assert_eq!(ordinary_term.orientations.len(), 20);
        assert_eq!(
            exact_term.orientations.len(),
            ordinary_term.orientations.len(),
            "source-preserving exact and ordinary CFF catalogues must agree"
        );
        let ordinary_sum = ordinary_term
            .orientations
            .iter()
            .map(|orientation| {
                orientation.expression.clone()
                    * numerator
                        .replace_multiple(orientation.orientation.energy_replacements_gs(&graph))
            })
            .fold(Atom::Zero, |sum, term| sum + term)
            * Atom::num(ordinary.production_prefactor_factor());
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

        let mass_squared = Atom::num(Rational::from((4, 9)));
        let external_energy = Atom::num(Rational::from((5, 2)));
        let external_edges = graph
            .external_momentum_edge_order()
            .into_iter()
            .collect::<BTreeSet<_>>();
        let evaluate_arb = |mut expression: Atom| -> Result<Complex<F<ArbPrec>>> {
            expression = expression
                .replace(Atom::var(GS.numerator_sampling_scale))
                .with(Atom::num(Rational::from((13, 10))))
                .replace(function!(GS.tree_denom_wrapper, W_.x_))
                .with(W_.x_)
                .replace(GS.den(W_.a_, W_.b_, W_.c_, W_.d_))
                .with(W_.d_);
            for edge in 0..graph.underlying.n_edges() {
                let edge = EdgeIndex(edge);
                expression = expression
                    .replace(graph.underlying[edge].particle.mass_atom())
                    .with(mass_squared.clone().sqrt());
                if external_edges.contains(&edge) {
                    expression = expression
                        .replace(GS.emr_mom(edge, GS.cind(0)))
                        .with(external_energy.clone());
                    for spatial_index in 1..=3 {
                        expression = expression
                            .replace(GS.emr_mom(edge, GS.cind(spatial_index)))
                            .with(Atom::Zero)
                            .replace(GS.emr_vec_index(edge, GS.cind(spatial_index)))
                            .with(Atom::Zero);
                    }
                } else {
                    let spatial = match usize::from(edge) {
                        1 => Atom::num(Rational::from((4, 3))),
                        2 => Atom::num(Rational::from((-7, 12))),
                        3..=5 => Atom::num(Rational::from((3, 4))),
                        _ => Atom::Zero,
                    };
                    let on_shell_energy = (mass_squared.clone() + spatial.clone().pow(2)).sqrt();
                    expression = expression
                        .replace(GS.emr_mom(edge, GS.cind(1)))
                        .with(spatial.clone())
                        .replace(GS.emr_vec_index(edge, GS.cind(1)))
                        .with(spatial)
                        .replace(GS.ose(edge))
                        .with(on_shell_energy.clone())
                        .replace(cut_energy(edge))
                        .with(on_shell_energy);
                    for spatial_index in 2..=3 {
                        expression = expression
                            .replace(GS.emr_mom(edge, GS.cind(spatial_index)))
                            .with(Atom::Zero)
                            .replace(GS.emr_vec_index(edge, GS.cind(spatial_index)))
                            .with(Atom::Zero);
                    }
                }
            }
            let parameters = [Atom::var(GS.pi)];
            let rational: ExpressionEvaluator<SymComplex<Fraction<IntegerRing>>> =
                expression.evaluator(&parameters).build().map_err(|error| {
                    eyre::eyre!("failed to build source-identity evaluator: {error}")
                })?;
            let mut arb: ExpressionEvaluator<Complex<F<ArbPrec>>> =
                rational.map_coeff(&|coefficient| {
                    Complex::new(F::from(&coefficient.re), F::from(&coefficient.im))
                });
            let zero = F(ArbPrec::default());
            Ok(arb.evaluate_single(&[Complex::new(zero.clone().pi(), zero)]))
        };
        let ordinary_value = evaluate_arb(ordinary_sum)?;
        let exact_value = evaluate_arb(exact_sum)?;
        let distance = (ordinary_value.clone() - exact_value.clone()).norm().re;
        let ordinary_norm = ordinary_value.clone().norm().re;
        let exact_norm = exact_value.clone().norm().re;
        let scale = if ordinary_norm > exact_norm {
            ordinary_norm
        } else {
            exact_norm
        };
        let relative_distance = if scale.is_zero() {
            distance
        } else {
            distance / scale
        };
        // The two independently assembled contour sums can lose precision through
        // different, highly asymmetric cancellations, so keep a generous precision-scaled gate.
        let tolerance = F(ArbPrec::default().epsilon()).sqrt().sqrt().sqrt();
        assert!(
            relative_distance <= tolerance,
            "source-preserving exact CFF differs from ordinary CFF: ordinary={ordinary_value:e}, exact={exact_value:e}, relative delta={relative_distance:e}, tolerance={tolerance:e}"
        );
        Ok(())
    }

    #[test]
    fn triangle_affine_energy_identity_survives_production_mapping() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(
            digraph diagnostic_triangle_affine_energy_identity {
                num = 1
                edge [particle="scalar_1" num=1]
                node [num=1]

                ext [style=invis]
                v0;
                v1;
                v2;
                ext -> v0 [id=0]
                ext -> v1 [id=1]
                v2 -> ext [id=2]
                v0 -> v2 [id=3 lmb_id=0]
                v1 -> v0 [id=4]
                v2 -> v1 [id=5]
            },
            "scalars"
        )?;
        let q = GS.emr_mom(EdgeIndex(3), GS.cind(0));
        let shifted = GS.emr_mom(EdgeIndex(4), GS.cind(0));
        let external = GS.emr_mom(EdgeIndex(0), GS.cind(0));
        let squared = q.clone().pow(2);
        let affine = q * (shifted + external);
        let build_sum = |numerator: &Atom| -> Result<Atom> {
            let mut graph = graph.clone();
            let options = graph.denominator_only_cff_3d_expression_options();
            let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
            let production = graph.generate_3d_expression_for_integrand(
                &[],
                &canonization,
                &options,
                Some(numerator),
            )?;
            let cff = graph.cff_from_production_expression(
                &production,
                &CutSet::empty(graph.n_hedges()),
                &OrientationPattern::default(),
            )?;
            let term = cff
                .terms
                .get(&CutCFFIndex::new_all_none())
                .expect("the triangle has one uncut CFF sector");
            Ok(term
                .orientations
                .iter()
                .map(|orientation| {
                    orientation.expression.clone()
                        * numerator.replace_multiple(
                            orientation.orientation.energy_replacements_gs(&graph),
                        )
                })
                .fold(Atom::Zero, |sum, term| sum + term)
                * Atom::num(cff.production_prefactor_factor()))
        };
        let evaluate_arb = |mut expression: Atom| -> Result<Complex<F<ArbPrec>>> {
            expression = expression
                .replace(function!(GS.tree_denom_wrapper, W_.x_))
                .with(W_.x_)
                .replace(GS.den(W_.a_, W_.b_, W_.c_, W_.d_))
                .with(W_.d_);
            for (edge, energy) in [
                (EdgeIndex(0), Atom::num(Rational::from((7, 10)))),
                (EdgeIndex(1), Atom::num(Rational::from((-23, 100)))),
                (EdgeIndex(2), Atom::num(Rational::from((47, 100)))),
            ] {
                expression = expression
                    .replace(GS.emr_mom(edge, GS.cind(0)))
                    .with(energy);
            }
            let mass_squared = Atom::num(1);
            for (edge, spatial) in [
                (EdgeIndex(3), [(31, 100), (-47, 100), (83, 100)]),
                (EdgeIndex(4), [(20, 100), (-39, 100), (78, 100)]),
                (EdgeIndex(5), [(24, 100), (-46, 100), (87, 100)]),
            ] {
                let energy = spatial
                    .into_iter()
                    .fold(mass_squared.clone(), |sum, value| {
                        sum + Atom::num(Rational::from(value)).pow(2)
                    })
                    .sqrt();
                expression = expression
                    .replace(GS.ose(edge))
                    .with(energy.clone())
                    .replace(cut_energy(edge))
                    .with(energy);
            }
            let parameters = [Atom::var(GS.pi)];
            let rational: ExpressionEvaluator<SymComplex<Fraction<IntegerRing>>> = expression
                .evaluator(&parameters)
                .build()
                .map_err(|error| eyre::eyre!("failed to build triangle evaluator: {error}"))?;
            let mut arb: ExpressionEvaluator<Complex<F<ArbPrec>>> =
                rational.map_coeff(&|coefficient| {
                    Complex::new(F::from(&coefficient.re), F::from(&coefficient.im))
                });
            let zero = F(ArbPrec::default());
            let pi = zero.clone().pi();
            Ok(arb.evaluate_single(&[Complex::new(pi, zero)]))
        };
        let squared_sum = build_sum(&squared)?;
        let affine_sum = build_sum(&affine)?;
        let squared_value = evaluate_arb(squared_sum.clone())?;
        let affine_value = evaluate_arb(affine_sum.clone())?;
        let distance = (squared_value.clone() - affine_value.clone()).norm().re;
        let scale = {
            let squared_norm = squared_value.clone().norm().re;
            let affine_norm = affine_value.clone().norm().re;
            if squared_norm > affine_norm {
                squared_norm
            } else {
                affine_norm
            }
        };
        let relative_distance = if scale.is_zero() {
            distance
        } else {
            distance / scale
        };
        let tolerance = F(ArbPrec::default().epsilon()).sqrt().sqrt().sqrt();
        assert!(
            relative_distance <= tolerance,
            "GammaLoop's mapped triangle CFF violates Q3^0=Q4^0+Q0^0: squared={squared_value:e}, affine={affine_value:e}, relative delta={relative_distance:e}, tolerance={tolerance:e}"
        );

        let rational = |numerator: i64, denominator: i64| {
            F::<ArbPrec>::from(&Rational::from((numerator, denominator)))
        };
        let loop_moms: LoopMomenta<F<ArbPrec>> = [ThreeMomentum::new(
            rational(31, 100),
            rational(-47, 100),
            rational(83, 100),
        )]
        .into_iter()
        .collect();
        let external_moms: ExternalFourMomenta<F<ArbPrec>> = [
            [
                rational(7, 10),
                rational(11, 100),
                rational(-8, 100),
                rational(5, 100),
            ]
            .into(),
            [
                rational(-23, 100),
                rational(-4, 100),
                rational(7, 100),
                rational(-9, 100),
            ]
            .into(),
            [
                rational(47, 100),
                rational(7, 100),
                rational(-1, 100),
                rational(-4, 100),
            ]
            .into(),
        ]
        .into_iter()
        .collect();
        let sample = MomentumSample {
            sample: BareMomentumSample {
                loop_moms,
                dual_loop_moms: None,
                loop_mom_cache_id: 0,
                loop_mom_base_cache_id: 0,
                external_moms,
                external_mom_cache_id: 0,
                external_mom_base_cache_id: 0,
                jacobian: rational(1, 1),
                orientation: None,
                parameterization_branch: None,
            },
        };
        let orientations = TiVec::<OrientationID, EdgeVec<Orientation>>::new();
        let orientation_filter = SubSet::full(orientations.len());
        let runtime_settings = RuntimeSettings::default();
        let evaluate_stack = |expression: &Atom| -> Result<Complex<F<ArbPrec>>> {
            let mut param_builder = graph.param_builder.clone();
            let (mut evaluator, _) = EvaluatorStack::new_explicit_sum_with_timings(
                std::slice::from_ref(expression),
                &param_builder,
                None,
                &EvaluatorSettings::default(),
            )?;
            let input = <ArbPrec as GenericEvaluatorFloat>::get_parameters(
                &mut param_builder,
                (false, false),
                &graph,
                &sample,
                &[],
                &[],
                None,
                None,
                None,
            );
            Ok(evaluator
                .evaluate(
                    input,
                    SingleOrAllOrientations::All {
                        all: &orientations,
                        filter: &orientation_filter,
                    },
                    &runtime_settings,
                    &mut EvaluationMetaData::new_empty(),
                    false,
                )?
                .pop()
                .expect("the triangle evaluator should return one value")
                .unwrap_real())
        };
        let squared_stack_value = evaluate_stack(&squared_sum)?;
        let affine_stack_value = evaluate_stack(&affine_sum)?;

        for (comparison, left, right) in [
            (
                "squared direct vs production stack",
                &squared_value,
                &squared_stack_value,
            ),
            (
                "affine direct vs production stack",
                &affine_value,
                &affine_stack_value,
            ),
            (
                "squared vs affine production stack",
                &squared_stack_value,
                &affine_stack_value,
            ),
        ] {
            let distance = (left.clone() - right.clone()).norm().re;
            let left_norm = left.clone().norm().re;
            let right_norm = right.clone().norm().re;
            let scale = if left_norm > right_norm {
                left_norm
            } else {
                right_norm
            };
            let relative_distance = if scale.is_zero() {
                distance
            } else {
                distance / scale
            };
            assert!(
                relative_distance <= tolerance,
                "{comparison} failed: left={left:e}, right={right:e}, relative delta={relative_distance:e}, tolerance={tolerance:e}"
            );
        }
        Ok(())
    }

    #[test]
    fn indexed_residue_selection_tracks_axes_independent_of_raised_order() -> Result<()> {
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
            assert!(residues.iter().all(|(index, _)| {
                let present_axes = [
                    index.right_threshold_order,
                    index.left_threshold_order,
                    index.lu_cut_order,
                ]
                .into_iter()
                .flatten()
                .count();
                present_axes == selected_axis_count
            }));
            if selected_axis_count != 0 {
                assert!(
                    residues.iter().any(|(index, _)| {
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
        let exact_terms = |cff: &CutCFF| -> Result<Vec<Atom>> {
            let prefactor = Atom::num(cff.production_prefactor_factor());
            cff.terms
                .values()
                .flat_map(|term| {
                    term.orientations.iter().map(|orientation| {
                        Ok(orientation.expression.clone()
                            * term.map_exact_source_numerator(&orientation.orientation)?
                            * &prefactor)
                    })
                })
                .collect::<Result<Vec<_>>>()
        };

        let (powered, _) =
            graph.cff_from_4d_denominators(&denominators, &cutset, &options, &numerator)?;
        let powered_terms = exact_terms(&powered)?;

        // Re-spell the complete powered rational component as D(-Q). Its
        // canonical exact graph must be identical, while the physical odd
        // numerator remains Q^0+c rather than acquiring the denominator sign.
        let mut reversed_denominators = denominators.clone();
        for denominator in &mut reversed_denominators[1..] {
            denominator.momentum = -denominator.momentum.clone();
        }
        assert_eq!(
            GraphThreeDSource::from_exact_denominators(&graph, &denominators)?
                .to_three_d_parsed_graph()?,
            GraphThreeDSource::from_exact_denominators(&graph, &reversed_denominators)?
                .to_three_d_parsed_graph()?,
            "globally reversing one rational-routing component must not change its canonical exact graph"
        );
        let (reversed, _) = graph.cff_from_4d_denominators(
            &reversed_denominators,
            &cutset,
            &options,
            &numerator,
        )?;
        let reversed_terms = exact_terms(&reversed)?;
        let (lower, _) = graph.cff_from_4d_denominators(
            &denominators[..2],
            &cutset,
            &options,
            &retained_factor,
        )?;
        let lower_terms = exact_terms(&lower)?;
        let denominator_only = GS.den(
            usize::from(edge),
            denominators[1].momentum.clone(),
            denominators[1].mass_squared.clone(),
            denominators[1].full_expr.clone(),
        );
        let (powered_scalar, _) =
            graph.cff_from_4d_denominators(&denominators, &cutset, &options, &denominator_only)?;
        let powered_scalar_terms = exact_terms(&powered_scalar)?;
        let (lower_scalar, _) =
            graph.cff_from_4d_denominators(&denominators[..2], &cutset, &options, &Atom::one())?;
        let lower_scalar_terms = exact_terms(&lower_scalar)?;
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
        let ordinary_terms = ordinary
            .terms
            .values()
            .flat_map(|term| {
                term.orientations.iter().map(|orientation| {
                    orientation.expression.clone()
                        * retained_factor.replace_multiple(
                            orientation.orientation.energy_replacements_gs(&graph),
                        )
                        * &ordinary_prefactor
                })
            })
            .collect::<Vec<_>>();

        // The same identity at the first nontrivial generalized-contact rank
        // exercises the exact numerator mapper rather than only the linear
        // residue handled by the ordinary repeated-pole path.
        let mut quintic_denominators = vec![denominators[0].clone()];
        quintic_denominators.extend((0..5).map(|_| denominators[1].clone()));
        let retained_quartic = retained_factor.pow(4);
        let quintic_numerator = GS.den(
            usize::from(edge),
            denominators[1].momentum.clone(),
            denominators[1].mass_squared.clone(),
            denominators[1].full_expr.clone(),
        ) * &retained_quartic;
        let (quintic, _) = graph.cff_from_4d_denominators(
            &quintic_denominators,
            &cutset,
            &options,
            &quintic_numerator,
        )?;
        let quintic_terms = exact_terms(&quintic)?;
        let mut reversed_quintic_denominators = quintic_denominators.clone();
        for denominator in &mut reversed_quintic_denominators[1..] {
            denominator.momentum = -denominator.momentum.clone();
        }
        let (reversed_quintic, _) = graph.cff_from_4d_denominators(
            &reversed_quintic_denominators,
            &cutset,
            &options,
            &quintic_numerator,
        )?;
        let reversed_quintic_terms = exact_terms(&reversed_quintic)?;
        let (quartic, _) = graph.cff_from_4d_denominators(
            &quintic_denominators[..5],
            &cutset,
            &options,
            &retained_quartic,
        )?;
        let quartic_terms = exact_terms(&quartic)?;

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
        let evaluate_arb = |expression: Atom| -> Result<Complex<F<ArbPrec>>> {
            let parameters = [Atom::var(GS.pi)];
            let rational: ExpressionEvaluator<SymComplex<Fraction<IntegerRing>>> = expression
                .evaluator(&parameters)
                .build()
                .map_err(|error| eyre::eyre!("failed to build odd-routing oracle: {error}"))?;
            let mut arb: ExpressionEvaluator<Complex<F<ArbPrec>>> =
                rational.map_coeff(&|coefficient| {
                    Complex::new(F::from(&coefficient.re), F::from(&coefficient.im))
                });
            Ok(arb.evaluate_single(&[Complex::new(
                F::<ArbPrec>::from_f64(0.0).pi(),
                F::<ArbPrec>::from_f64(0.0),
            )]))
        };
        // Production retains this component structure and evaluates every
        // orientation separately.  The regression does the same in Arb so it
        // also detects precision loss without first expanding the factorized
        // numerator into one prohibitively large symbolic sum.
        let evaluate_terms = |terms: &[Atom],
                              unwrap_denominator: bool,
                              spatial: &Atom,
                              energy: &Atom,
                              constant: &Atom|
         -> Result<Complex<F<ArbPrec>>> {
            let zero = F::<ArbPrec>::from_f64(0.0);
            terms
                .iter()
                .try_fold(Complex::new(zero.clone(), zero), |sum, expression| {
                    let expression = if unwrap_denominator {
                        expression
                            .replace(GS.den(W_.a_, W_.b_, W_.c_, W_.d_))
                            .with(W_.d_)
                    } else {
                        expression.clone()
                    };
                    Ok(sum + evaluate_arb(fixed_point(expression, spatial, energy, constant))?)
                })
        };
        let tolerance = F(ArbPrec::default().epsilon()).sqrt().sqrt().sqrt();
        let points = [
            (Atom::Zero, Atom::one(), Atom::num(2)),
            (
                Atom::num(symbolica::domains::rational::Rational::from((3, 4))),
                Atom::num(symbolica::domains::rational::Rational::from((5, 4))),
                Atom::num(symbolica::domains::rational::Rational::from((7, 3))),
            ),
        ];
        for (spatial, energy, constant) in points {
            let powered_value = evaluate_terms(&powered_terms, true, &spatial, &energy, &constant)?;
            let reversed_value =
                evaluate_terms(&reversed_terms, true, &spatial, &energy, &constant)?;
            let lower_value = evaluate_terms(&lower_terms, false, &spatial, &energy, &constant)?;
            let powered_scalar_value =
                evaluate_terms(&powered_scalar_terms, true, &spatial, &energy, &constant)?;
            let lower_scalar_value =
                evaluate_terms(&lower_scalar_terms, false, &spatial, &energy, &constant)?;
            for (label, candidate, reference) in [
                (
                    "uncancelled D*(Q0+c)/D^3",
                    powered_value.clone(),
                    lower_value,
                ),
                (
                    "uncancelled scalar D/D^3",
                    powered_scalar_value,
                    lower_scalar_value,
                ),
            ] {
                let distance = (candidate.clone() - reference.clone()).norm().re;
                let candidate_norm = candidate.clone().norm().re;
                let reference_norm = reference.clone().norm().re;
                let scale = if candidate_norm > reference_norm {
                    candidate_norm
                } else {
                    reference_norm
                };
                let relative_distance = if scale.is_zero() {
                    distance
                } else {
                    distance / scale
                };
                assert!(
                    relative_distance <= tolerance,
                    "{label} differs from its independently generated lower source at spatial={spatial}, c={constant}: candidate={candidate:e}, lower={reference:e}, relative delta={relative_distance:e}, tolerance={tolerance:e}"
                );
            }
            let distance = (powered_value.clone() - reversed_value.clone()).norm().re;
            let powered_norm = powered_value.clone().norm().re;
            let reversed_norm = reversed_value.clone().norm().re;
            let scale = if powered_norm > reversed_norm {
                powered_norm
            } else {
                reversed_norm
            };
            let relative_distance = if scale.is_zero() {
                distance
            } else {
                distance / scale
            };
            assert!(
                relative_distance <= tolerance,
                "globally reversed D(-Q) component changes odd physical Q0+c numerator at spatial={spatial}, c={constant}: powered={powered_value:e}, reversed={reversed_value:e}, relative delta={relative_distance:e}, tolerance={tolerance:e}"
            );
            let ordinary_value =
                evaluate_terms(&ordinary_terms, false, &spatial, &energy, &constant)?;
            let lower_value = evaluate_terms(&lower_terms, false, &spatial, &energy, &constant)?;
            let normalization_distance = (lower_value.clone() - ordinary_value.clone()).norm().re;
            let lower_norm = lower_value.clone().norm().re;
            let ordinary_norm = ordinary_value.clone().norm().re;
            let normalization_scale = if lower_norm > ordinary_norm {
                lower_norm
            } else {
                ordinary_norm
            };
            let normalization_relative_distance = if normalization_scale.is_zero() {
                normalization_distance
            } else {
                normalization_distance / normalization_scale
            };
            assert!(
                normalization_relative_distance <= tolerance,
                "exact and ordinary (Q0+c)/(D7*D8) CFF sums differ at spatial={spatial}, c={constant}: exact={lower_value:e}, ordinary={ordinary_value:e}, relative delta={normalization_relative_distance:e}, tolerance={tolerance:e}"
            );

            let quintic_value = evaluate_terms(&quintic_terms, true, &spatial, &energy, &constant)?;
            let reversed_quintic_value =
                evaluate_terms(&reversed_quintic_terms, true, &spatial, &energy, &constant)?;
            let quartic_value =
                evaluate_terms(&quartic_terms, false, &spatial, &energy, &constant)?;
            for (label, candidate) in [
                ("uncancelled", quintic_value),
                ("globally reversed", reversed_quintic_value),
            ] {
                let distance = (candidate.clone() - quartic_value.clone()).norm().re;
                let candidate_norm = candidate.clone().norm().re;
                let quartic_norm = quartic_value.clone().norm().re;
                let scale = if candidate_norm > quartic_norm {
                    candidate_norm
                } else {
                    quartic_norm
                };
                let relative_distance = if scale.is_zero() {
                    distance
                } else {
                    distance / scale
                };
                assert!(
                    relative_distance <= tolerance,
                    "{label} D*(Q0+c)^4/D^5 differs from (Q0+c)^4/D^4 at spatial={spatial}, c={constant}: candidate={candidate:e}, lower={quartic_value:e}, relative delta={relative_distance:e}, tolerance={tolerance:e}"
                );
            }
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
                vec![2],
                "an original quadratic energy must stay on its canonical source occurrence"
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
            let fixed_hard_energy = |owner: EdgeIndex| {
                FunctionBuilder::new(GS.emr_mom)
                    .add_arg(GS.uv_momentum_provenance_tag(
                        Atom::num(usize::from(owner) as i64).as_view(),
                        false,
                        momentum.as_view(),
                    ))
                    .add_arg(GS.cind(0))
                    .finish()
            };
            let tagged_quadratic_numerator =
                fixed_hard_energy(owner_edges[0]) * fixed_hard_energy(owner_edges[1]);
            let (tagged_quadratic_cff, _) = graph.cff_from_4d_denominators_in_uv_edges(
                &denominators,
                owner_edges,
                &cutset,
                &options,
                &tagged_quadratic_numerator,
                None,
            )?;
            let tagged_quadratic_sum = tagged_quadratic_cff
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
                * Atom::num(tagged_quadratic_cff.production_prefactor_factor());
            let tagged_difference = (&tagged_quadratic_sum - &quadratic_sum).together();
            assert!(
                tagged_difference.is_zero(),
                "two fixed provenance factors carrying the same hard Q must reproduce Q0^2: {tagged_difference}",
            );
            // The generalized CFF owns the lower-sector/contact realization of
            // the powered pole. Keep this test on production invariants:
            // tagged and untagged forms agree above, and permuting provenance
            // owners cannot change the resulting Laurent functional below.
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
            let momentum_squared_reference = &quadratic_sum - &spatial_norm * &orientation_sum;
            let momentum_squared_difference = (momentum_squared_sum - momentum_squared_reference)
                .together()
                .expand();
            assert!(
                momentum_squared_difference.is_zero(),
                "the full momentum-square numerator must preserve CFF linearity between its temporal and spatial factors: {momentum_squared_difference}"
            );
        }
        Ok(())
    }

    #[test]
    fn exact_theta_equal_channel_temporal_numerators_are_owner_invariant() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph exact_theta_equal_channel_numerator {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0 lmb_id=0]
            a -> c [id=1 lmb_id=1]
            a -> d [id=2]
            b -> c [id=3]
            b -> d [id=4]
        })?;
        let q0 = FunctionBuilder::new(GS.emr_mom).add_arg(0).finish();
        let q1 = FunctionBuilder::new(GS.emr_mom).add_arg(1).finish();
        let momenta = [q0.clone(), q1.clone(), -&q0 - &q1, -&q1, &q0 + &q1];
        let denominators = momenta
            .into_iter()
            .enumerate()
            .map(|(owner, momentum)| FourDDenominator {
                source_edge: EdgeIndex(owner),
                momentum,
                mass_squared: Atom::one(),
                full_expr: Atom::one(),
            })
            .collect::<Vec<_>>();
        let tagged_temporal = |owner: usize, hard: &Atom| {
            FunctionBuilder::new(GS.emr_mom)
                .add_arg(GS.uv_momentum_provenance_tag(
                    Atom::num(owner as i64).as_view(),
                    false,
                    hard.as_view(),
                ))
                .add_arg(GS.cind(0))
                .finish()
        };
        let options = graph.denominator_only_cff_3d_expression_options();
        let cutset = CutSet::empty(graph.n_hedges());
        let exact_cff = |graph: &mut Graph, numerator: &Atom| -> Result<CutCFF> {
            Ok(graph
                .cff_from_4d_denominators(&denominators, &cutset, &options, numerator)?
                .0)
        };
        let exact_sum = |cff: &CutCFF| -> Result<Atom> {
            Ok(cff
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
                * Atom::num(cff.production_prefactor_factor()))
        };
        let exact_unit_cff = exact_cff(&mut graph, &Atom::one())?;
        let exact_unit = exact_sum(&exact_unit_cff)?;
        let q1_squared = tagged_temporal(1, &q1).pow(2);
        let q3_squared = tagged_temporal(3, &(-q1)).pow(2);
        let q1_cff = exact_cff(&mut graph, &q1_squared)?;
        let q3_cff = exact_cff(&mut graph, &q3_squared)?;
        let empty: linnet::half_edge::subgraph::SuBitGraph = graph.empty_subgraph();
        let ordinary = graph.cff(
            &empty,
            &cutset,
            &OrientationPattern::default(),
            &options,
            Some(&Atom::one()),
        )?;
        let ordinary_unit = ordinary
            .terms
            .values()
            .flat_map(|term| &term.orientations)
            .fold(Atom::Zero, |sum, orientation| sum + &orientation.expression)
            * Atom::num(ordinary.production_prefactor_factor());
        let mut fixed_unit_difference = &exact_unit - &ordinary_unit;
        for (edge, x, energy) in [
            (
                EdgeIndex(0),
                Atom::num(2) / Atom::num(3),
                (Atom::num(13) / Atom::num(9)).sqrt(),
            ),
            (
                EdgeIndex(1),
                Atom::num(3) / Atom::num(4),
                Atom::num(5) / Atom::num(4),
            ),
        ] {
            fixed_unit_difference = fixed_unit_difference
                .replace(GS.emr_mom(edge, GS.cind(1)))
                .with(x)
                .replace(GS.emr_mom(edge, GS.cind(2)))
                .with(Atom::Zero)
                .replace(GS.emr_mom(edge, GS.cind(3)))
                .with(Atom::Zero)
                .replace(GS.ose(edge))
                .with(energy);
        }
        let combined_energy = (Atom::num(433) / Atom::num(144)).sqrt();
        fixed_unit_difference = fixed_unit_difference
            .replace(GS.ose(EdgeIndex(2)))
            .with(combined_energy.clone())
            .replace(GS.ose(EdgeIndex(3)))
            .with(Atom::num(5) / Atom::num(4))
            .replace(GS.ose(EdgeIndex(4)))
            .with(combined_energy);
        assert!(
            fixed_unit_difference.together().is_zero(),
            "the reconstructed and ordinary unit-numerator theta CFFs must agree: {fixed_unit_difference}",
        );
        for [q0x, q1x] in [
            [Atom::num(2) / Atom::num(3), Atom::num(3) / Atom::num(4)],
            [Atom::num(4) / Atom::num(5), Atom::num(5) / Atom::num(6)],
        ] {
            let q0_energy = (Atom::one() + q0x.clone().pow(2)).sqrt();
            let q1_energy = (Atom::one() + q1x.clone().pow(2)).sqrt();
            let fix = |mut expression: Atom| {
                for (edge, x, energy) in [
                    (EdgeIndex(0), q0x.clone(), q0_energy.clone()),
                    (EdgeIndex(1), q1x.clone(), q1_energy.clone()),
                ] {
                    expression = expression
                        .replace(GS.emr_mom(edge, GS.cind(1)))
                        .with(x)
                        .replace(GS.emr_mom(edge, GS.cind(2)))
                        .with(Atom::Zero)
                        .replace(GS.emr_mom(edge, GS.cind(3)))
                        .with(Atom::Zero)
                        .replace(GS.ose(edge))
                        .with(energy);
                }
                expression
            };
            let q1_orientations = q1_cff
                .terms
                .values()
                .flat_map(|term| {
                    term.orientations
                        .iter()
                        .map(move |orientation| (term, orientation))
                })
                .collect::<Vec<_>>();
            let q3_orientations = q3_cff
                .terms
                .values()
                .flat_map(|term| {
                    term.orientations
                        .iter()
                        .map(move |orientation| (term, orientation))
                })
                .collect::<Vec<_>>();
            assert_eq!(q1_orientations.len(), q3_orientations.len());
            for (index, ((q1_term, q1_orientation), (q3_term, q3_orientation))) in
                q1_orientations.iter().zip(&q3_orientations).enumerate()
            {
                assert_eq!(
                    q1_orientation.orientation.loop_energy_map,
                    q3_orientation.orientation.loop_energy_map,
                    "owner-dependent loop-energy map at orientation {index}",
                );
                assert_eq!(
                    q1_orientation.orientation.edge_energy_map,
                    q3_orientation.orientation.edge_energy_map,
                    "owner-dependent edge-energy map at orientation {index}",
                );
                let q1_carrier = fix(q1_orientation.expression.clone()).together();
                let q3_carrier = fix(q3_orientation.expression.clone()).together();
                let carrier_difference = &q1_carrier - &q3_carrier;
                assert!(
                    carrier_difference.together().is_zero(),
                    "owner-dependent GL carrier at orientation {index}: {carrier_difference}",
                );
                let q1_numerator =
                    fix(q1_term.map_exact_source_numerator(&q1_orientation.orientation)?)
                        .together();
                let q3_numerator =
                    fix(q3_term.map_exact_source_numerator(&q3_orientation.orientation)?)
                        .together();
                let numerator_difference = &q1_numerator - &q3_numerator;
                assert!(
                    numerator_difference.together().is_zero(),
                    "owner-dependent mapped numerator at orientation {index}: {numerator_difference}",
                );
                // Preserve the production factorization at the symbolic oracle boundary.
                // Multiplying before applying the shell point can leave equivalent positive
                // radicals in the distinct forms `sqrt(x) / x` and `1 / sqrt(x)`, which
                // `together()` deliberately does not identify through a branch assumption.
                let product_difference = &q1_carrier * &q1_numerator - &q3_carrier * &q3_numerator;
                assert!(
                    product_difference.together().is_zero(),
                    "owner-dependent carrier product at orientation {index}: {product_difference}",
                );
            }
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

        for (expected_order, numerator) in [
            (1, Atom::one()),
            (1, GS.emr_mom(EdgeIndex(1), GS.cind(0)).pow(4)),
        ] {
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
                let mut exact_sum = exact_term
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
                for edge in [EdgeIndex(1), EdgeIndex(2)] {
                    let on_shell_energy = (1..=3)
                        .fold(
                            graph.underlying[edge].particle.mass_atom().pow(2),
                            |norm_squared, spatial_index| {
                                norm_squared + GS.emr_mom(edge, GS.cind(spatial_index)).pow(2)
                            },
                        )
                        .sqrt();
                    exact_sum = exact_sum.replace(GS.ose(edge)).with(on_shell_energy);
                }
                let difference = (&exact_sum - &ordinary_sum).together();
                assert!(
                    difference.is_zero(),
                    "exact and ordinary LU residues differ for maximum order {expected_order}, index {index}: {difference}"
                );
            }
        }
        Ok(())
    }

    #[test]
    fn exact_uncut_cubic_one_loop_source_matches_ordinary_cff() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph exact_uncut_cubic_one_loop {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0 lmb_id=0]
            b -> c [id=1]
            c -> a [id=2]
        })?;
        let options = graph.denominator_only_cff_3d_expression_options();
        let numerator = GS.emr_mom(EdgeIndex(0), GS.cind(0)).pow(2);
        let cutset = CutSet::empty(graph.n_hedges());
        let contract: SuBitGraph = graph.empty_subgraph();
        let ordinary = graph.cff(
            &contract,
            &cutset,
            &OrientationPattern::default(),
            &options,
            Some(&numerator),
        )?;
        let source_momentum = FunctionBuilder::new(GS.emr_mom)
            .add_arg(usize::from(EdgeIndex(0)))
            .finish();
        let denominators = [(EdgeIndex(0), 1), (EdgeIndex(1), -1), (EdgeIndex(2), -1)].map(
            |(source_edge, sign)| FourDDenominator {
                source_edge,
                momentum: Atom::num(sign) * &source_momentum,
                mass_squared: graph.underlying[source_edge].particle.mass_atom().pow(2),
                full_expr: Atom::one(),
            },
        );
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

        let index = CutCFFIndex::new_all_none();
        let ordinary_sum = ordinary
            .terms
            .get(&index)
            .expect("the ordinary uncut source exists")
            .orientations
            .iter()
            .map(|orientation| {
                orientation.expression.clone()
                    * numerator
                        .replace_multiple(orientation.orientation.energy_replacements_gs(&graph))
            })
            .fold(Atom::Zero, |sum, term| sum + term)
            * Atom::num(ordinary.production_prefactor_factor());
        let exact_term = exact
            .terms
            .get(&index)
            .expect("the exact uncut source exists");
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

        let evaluate_arb = |mut expression: Atom| -> Result<Complex<F<ArbPrec>>> {
            for edge in 0..3 {
                let edge = EdgeIndex(edge);
                expression = expression
                    .replace(GS.emr_mom(edge, GS.cind(1)))
                    .with(Atom::num(Rational::from((3, 4))))
                    .replace(GS.emr_mom(edge, GS.cind(2)))
                    .with(Atom::Zero)
                    .replace(GS.emr_mom(edge, GS.cind(3)))
                    .with(Atom::Zero)
                    .replace(GS.ose(edge))
                    .with(Atom::num(Rational::from((5, 4))));
            }
            let parameters = [Atom::var(GS.pi)];
            let rational: ExpressionEvaluator<SymComplex<Fraction<IntegerRing>>> =
                expression.evaluator(&parameters).build().map_err(|error| {
                    eyre::eyre!("failed to build cubic one-loop evaluator: {error}")
                })?;
            let mut arb: ExpressionEvaluator<Complex<F<ArbPrec>>> =
                rational.map_coeff(&|coefficient| {
                    Complex::new(F::from(&coefficient.re), F::from(&coefficient.im))
                });
            Ok(arb.evaluate_single(&[Complex::new(
                F::<ArbPrec>::from_f64(0.0).pi(),
                F::<ArbPrec>::from_f64(0.0),
            )]))
        };
        let ordinary_value = evaluate_arb(ordinary_sum)?;
        let exact_value = evaluate_arb(exact_sum)?;
        let distance = (exact_value.clone() - ordinary_value.clone()).norm().re;
        let exact_norm = exact_value.clone().norm().re;
        let ordinary_norm = ordinary_value.clone().norm().re;
        let scale = if exact_norm > ordinary_norm {
            exact_norm
        } else {
            ordinary_norm
        };
        let relative_distance = if scale.is_zero() {
            distance
        } else {
            distance / scale
        };
        let tolerance = F(ArbPrec::default().epsilon()).sqrt().sqrt().sqrt();
        assert!(
            relative_distance <= tolerance,
            "energy-convergent D(Q)D(-Q)^2 source differs from ordinary CFF: exact={exact_value:e}, ordinary={ordinary_value:e}, relative delta={relative_distance:e}, tolerance={tolerance:e}"
        );
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

        let combined_source =
            GraphThreeDSource::from_exact_denominators(&graph, &combined_denominators)?;
        let cograph_source =
            GraphThreeDSource::from_exact_denominators(&graph, &cograph_denominators)?;
        let uv_source = GraphThreeDSource::from_exact_denominators(&graph, &uv_denominators)?;
        assert_eq!(
            combined_source.active_loop_count(),
            cograph_source.active_loop_count() + uv_source.active_loop_count(),
            "the source phase must factor between the independent loop components",
        );
        let (combined_generated, _, _, _) = graph.generate_3d_expression_for_4d_term(
            &combined_source,
            &options,
            &quadratic_numerator,
            None,
        )?;
        let (cograph_generated, _, _, _) = graph.generate_3d_expression_for_4d_term(
            &cograph_source,
            &options,
            &Atom::one(),
            None,
        )?;
        let (uv_generated, _, _, _) = graph.generate_3d_expression_for_4d_term(
            &uv_source,
            &options,
            &quadratic_numerator,
            None,
        )?;
        let mut component_ownership = combined_generated
            .energy_factor_components
            .iter()
            .map(|component| (component.internal_edge_ids.len(), component.ownership))
            .collect::<Vec<_>>();
        component_ownership.sort_by_key(|(edge_count, _)| *edge_count);
        assert_eq!(
            component_ownership,
            vec![
                (2, CffEnergyFactorOwnership::GlobalSourceProduct),
                (3, CffEnergyFactorOwnership::VariantLocal),
            ],
            "the lifted source must retain its pure cograph and generalized UV ownership independently of canonical component order",
        );
        assert_eq!(cograph_generated.core_global_prefactor_sign.factor(), -1);
        assert_eq!(uv_generated.core_global_prefactor_sign.factor(), 1);
        assert_eq!(
            combined_generated.core_global_prefactor_sign.factor(),
            1,
            "the connected mixed-component source must retain the pure cograph's duplicate-line parity in its shared CFF frame",
        );

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
        let cff_sum = |cff: &CutCFF| -> Result<Atom> {
            let term = cff
                .terms
                .values()
                .next()
                .expect("the exact source retains its requested residue sector");
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
        let combined_sum = cff_sum(&combined)?;
        let cograph_sum = cff_sum(&cograph)?;
        let uv_sum = cff_sum(&uv)?;

        let empty_cutset = CutSet::empty(graph.n_hedges());
        let (combined_uncut, _) = graph.cff_from_4d_denominators(
            &combined_denominators,
            &empty_cutset,
            &options,
            &quadratic_numerator,
        )?;
        let (cograph_uncut, _) = graph.cff_from_4d_denominators(
            &cograph_denominators,
            &empty_cutset,
            &options,
            &Atom::one(),
        )?;
        let uncut_difference =
            (cff_sum(&combined_uncut)? - cff_sum(&cograph_uncut)? * &uv_sum).together();
        assert!(
            uncut_difference.is_zero(),
            "an uncut exact source must factorize between its independent rational components: {uncut_difference}",
        );
        let difference = (combined_sum - cograph_sum * uv_sum).together();
        assert!(
            difference.is_zero(),
            "an LU residue in one exact component must factorize from a quadratic cubic spectator: {difference}"
        );
        Ok(())
    }

    #[test]
    fn selected_highest_pole_commutes_with_factorized_cubic_numerator() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph selected_cubic_repeated_channel {
            num = 1
            edge [particle="scalar_1" num=1]
            node [num=1]

            ext0 [style=invis is_cut=0]
            v2 -> ext0 [id=0]
            ext0 -> v3
            v0 -> v1 [id=1 lmb_id=0]
            v0 -> v1 [id=2]
            v3 -> v0 [id=3 lmb_id=1]
            v1 -> v2 [id=4]
            v2 -> v3 [id=5]
        }, "scalars")?;
        let options = graph.denominator_only_cff_3d_expression_options();
        // Keep the three temporal factors separate. Edges 3 and 4 are the
        // repeated outer channel, so this supplies aggregate degree three
        // without expanding the numerator.
        let numerator =
            GS.emr_mom(EdgeIndex(3), GS.cind(0)).pow(2) * GS.emr_mom(EdgeIndex(4), GS.cind(0));
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let production = graph.generate_3d_expression_for_integrand(
            &[],
            &canonization,
            &options,
            Some(&numerator),
        )?;
        let raised_group = graph
            .determine_raised_esurfaces_from_expression(&production.expression)
            .raised_groups
            .into_iter()
            .find(|group| {
                group.max_occurence == 2
                    && group.esurface_ids.iter().any(|surface_id| {
                        !production.expression.surfaces.esurface_cache[*surface_id]
                            .external_shift
                            .is_empty()
                    })
            })
            .expect("the repeated outer channel supplies a second-order LU surface");
        let mut cutset = CutSet::empty(graph.n_hedges());
        cutset.residue_selector.lu = Some(crate::graph::cuts::LuCutSelection {
            cut_edge_alternatives: raised_group
                .esurface_ids
                .iter()
                .map(|surface_id| {
                    production.expression.surfaces.esurface_cache[*surface_id]
                        .energies
                        .clone()
                })
                .collect(),
            raised_group,
        });
        let contracted = graph
            .get_edge_subgraph(EdgeIndex(1))
            .union(&graph.get_edge_subgraph(EdgeIndex(2)));
        let contract_subgraph = contracted
            .union(&graph.tree_edges)
            .subtract(&graph.initial_state_cut);
        let pattern = OrientationPattern::default();
        let generalized = graph.cff(
            &contract_subgraph,
            &cutset,
            &pattern,
            &options,
            Some(&numerator),
        )?;
        let ordinary = graph.cff(
            &contract_subgraph,
            &cutset,
            &pattern,
            &options,
            Some(&Atom::one()),
        )?;
        let highest_pole = |cff: &CutCFF| {
            let term = cff
                .terms
                .iter()
                .find(|(index, _)| index.lu_cut_order == Some(2))
                .expect("the selected CFF retains its maximum LU order")
                .1;
            term.orientations
                .iter()
                .map(|orientation| {
                    &orientation.expression
                        * numerator.replace_multiple(
                            orientation.orientation.energy_replacements_gs(&graph),
                        )
                })
                .fold(Atom::Zero, |sum, term| sum + term)
                * Atom::num(cff.production_prefactor_factor())
        };
        let difference = (highest_pole(&generalized) - highest_pole(&ordinary)).together();
        assert!(
            difference.is_zero(),
            "the maximum-order LU residue must commute with its factorized cubic numerator: {difference}",
        );
        Ok(())
    }

    #[test]
    fn selected_lower_pole_obeys_repeated_channel_polynomial_division() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph selected_cubic_repeated_channel_contact {
            num = 1
            edge [particle="scalar_1" num=1]
            node [num=1]

            ext0 [style=invis is_cut=0]
            v2 -> ext0 [id=0]
            ext0 -> v3
            v0 -> v1 [id=1 lmb_id=0]
            v0 -> v1 [id=2]
            v3 -> v0 [id=3 lmb_id=1]
            v1 -> v2 [id=4]
            v2 -> v3 [id=5]
        }, "scalars")?;
        let options = graph.denominator_only_cff_3d_expression_options();
        let repeated_energy = GS.emr_mom(EdgeIndex(3), GS.cind(0));
        let alias_energy = GS.emr_mom(EdgeIndex(4), GS.cind(0));
        let numerator = repeated_energy.pow(2) * &alias_energy;
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let production = graph.generate_3d_expression_for_integrand(
            &[],
            &canonization,
            &options,
            Some(&numerator),
        )?;
        let raised_group = graph
            .determine_raised_esurfaces_from_expression(&production.expression)
            .raised_groups
            .into_iter()
            .find(|group| {
                group.max_occurence == 2
                    && group.esurface_ids.iter().any(|surface_id| {
                        !production.expression.surfaces.esurface_cache[*surface_id]
                            .external_shift
                            .is_empty()
                    })
            })
            .expect("the repeated outer channel supplies a second-order LU surface");
        let mut cutset = CutSet::empty(graph.n_hedges());
        cutset.residue_selector.lu = Some(crate::graph::cuts::LuCutSelection {
            cut_edge_alternatives: raised_group
                .esurface_ids
                .iter()
                .map(|surface_id| {
                    production.expression.surfaces.esurface_cache[*surface_id]
                        .energies
                        .clone()
                })
                .collect(),
            raised_group,
        });
        let contracted = graph
            .get_edge_subgraph(EdgeIndex(1))
            .union(&graph.get_edge_subgraph(EdgeIndex(2)));
        let contract_subgraph = contracted
            .union(&graph.tree_edges)
            .subtract(&graph.initial_state_cut);
        let pattern = OrientationPattern::default();
        let generalized = graph.cff(
            &contract_subgraph,
            &cutset,
            &pattern,
            &options,
            Some(&numerator),
        )?;

        // For the repeated channel D_3=D_4=D(q), polynomial division gives
        //
        //   q0^2 q0 / (D_3 D_4) = q0 / D_4 + E_q^2 q0 / (D_3 D_4).
        //
        // This test-only spelling is an exact lower-pole oracle. Production
        // keeps the original three numerator factors untouched and relies on
        // the generalized CFF to construct the same quotient and remainder.
        let remainder_numerator = GS.ose(EdgeIndex(3)).pow(2) * &alias_energy;
        let remainder = graph.cff(
            &contract_subgraph,
            &cutset,
            &pattern,
            &options,
            Some(&remainder_numerator),
        )?;
        let mut contact_contract = contract_subgraph.clone();
        contact_contract.add(graph[&EdgeIndex(3)].1);
        let contact = graph.cff(
            &contact_contract,
            &cutset,
            &pattern,
            &options,
            Some(&alias_energy),
        )?;
        let selected_lu1 = |cff: &CutCFF, mapped_numerator: &Atom| {
            cff.terms
                .iter()
                .filter(|(index, _)| index.lu_cut_order == Some(1))
                .flat_map(|(_, term)| {
                    term.orientations.iter().map(|orientation| {
                        &orientation.expression
                            * mapped_numerator.replace_multiple(
                                orientation.orientation.energy_replacements_gs(&graph),
                            )
                    })
                })
                .fold(Atom::Zero, |sum, term| sum + term)
                * Atom::num(cff.production_prefactor_factor())
        };
        let difference = (selected_lu1(&generalized, &numerator)
            - selected_lu1(&remainder, &remainder_numerator)
            - selected_lu1(&contact, &alias_energy))
        .replace(GS.ose(EdgeIndex(4)))
        .with(GS.ose(EdgeIndex(3)))
        .together();
        assert!(
            difference.is_zero(),
            "the first-order LU residue must preserve exact polynomial division on a repeated channel: {difference}",
        );
        Ok(())
    }

    #[test]
    fn exact_raised_lu_residue_factorizes_from_quadratic_cubic_spectator() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph exact_raised_lu_cubic_product {
            edge [num=1 mass=1]
            node [num=1]

            incoming [style=invis]
            outgoing [style=invis]

            incoming -> v1 [id=0]
            v1 -> v2 [id=1 lmb_id=0]
            v2 -> v3 [id=2]
            v1 -> v3 [id=3]
            v1 -> u [id=5 lmb_id=1]
            u -> v1 [id=6]
            v3 -> outgoing [id=4]
        })?;
        let options = graph.denominator_only_cff_3d_expression_options();
        let cograph_edge = EdgeIndex(1);
        let cograph_numerator = GS.emr_mom(cograph_edge, GS.cind(0)).pow(2);
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let production = graph.generate_3d_expression_for_integrand(
            &[EdgeIndex(5), EdgeIndex(6)],
            &canonization,
            &options,
            Some(&cograph_numerator),
        )?;
        let lu_surface = graph
            .determine_raised_esurfaces_from_expression(&production.expression)
            .raised_groups
            .into_iter()
            .find(|group| {
                group.max_occurence == 2
                    && group.esurface_ids.iter().any(|surface_id| {
                        production.expression.surfaces.esurface_cache[*surface_id]
                            .energies
                            .iter()
                            .copied()
                            .collect::<BTreeSet<_>>()
                            .is_subset(&BTreeSet::from([EdgeIndex(1), EdgeIndex(2), EdgeIndex(3)]))
                    })
            })
            .expect("the quadratic cograph supplies a second-order LU surface");
        let mut lu_cutset = CutSet::empty(graph.n_hedges());
        let physical_cut_support = production.expression.surfaces.esurface_cache
            [lu_surface.esurface_ids[0]]
            .energies
            .clone();
        lu_cutset.residue_selector.lu = Some(crate::graph::cuts::LuCutSelection {
            // A raised physical Cutkosky cut has one cut support even when its
            // two residue orders originate from distinct exact E-surfaces.
            cut_edge_alternatives: vec![physical_cut_support; lu_surface.esurface_ids.len()],
            raised_group: lu_surface,
        });

        let cograph_denominators =
            [EdgeIndex(1), EdgeIndex(2), EdgeIndex(3)].map(|edge| FourDDenominator {
                source_edge: edge,
                momentum: FunctionBuilder::new(GS.emr_mom)
                    .add_arg(usize::from(edge))
                    .finish(),
                mass_squared: graph.underlying[edge].particle.mass_atom().pow(2),
                full_expr: Atom::one(),
            });
        let uv_edge = EdgeIndex(5);
        let uv_momentum = FunctionBuilder::new(GS.emr_mom)
            .add_arg(usize::from(uv_edge))
            .finish();
        let uv_mass_squared = Atom::var(GS.m_uv_expansion).pow(2);
        let uv_full_expr = GS.emr_mom(uv_edge, GS.cind(0)).pow(2)
            - (1..=3).fold(uv_mass_squared.clone(), |norm_squared, spatial_index| {
                norm_squared + GS.emr_mom(uv_edge, GS.cind(spatial_index)).pow(2)
            });
        // Match the UV-expanded self-energy provenance pattern D_2 D_3^2.
        // The owner IDs retain cut support only; all three occurrences carry
        // the same even rational denominator D(Q)=D(-Q).
        let uv_denominators = [(EdgeIndex(5), 1), (EdgeIndex(6), -1), (EdgeIndex(6), -1)].map(
            |(source_edge, sign)| FourDDenominator {
                source_edge,
                momentum: Atom::num(sign) * &uv_momentum,
                mass_squared: uv_mass_squared.clone(),
                full_expr: uv_full_expr.clone(),
            },
        );
        let uv_numerator = GS.emr_mom(uv_edge, GS.cind(0)).pow(2);
        let combined_numerator = &cograph_numerator * &uv_numerator;
        let combined_denominators = cograph_denominators
            .iter()
            .cloned()
            .chain(uv_denominators.iter().cloned())
            .collect::<Vec<_>>();
        let uv_edges = [EdgeIndex(5), EdgeIndex(6)];

        let (combined, _) = graph.cff_from_4d_denominators_in_uv_edges(
            &combined_denominators,
            uv_edges,
            &lu_cutset,
            &options,
            &combined_numerator,
            None,
        )?;
        let (cograph, _) = graph.cff_from_4d_denominators(
            &cograph_denominators,
            &lu_cutset,
            &options,
            &cograph_numerator,
        )?;
        let (uv, _) = graph.cff_from_4d_denominators_in_uv_edges(
            &uv_denominators,
            uv_edges,
            &CutSet::empty(graph.n_hedges()),
            &options,
            &uv_numerator,
            None,
        )?;
        assert_eq!(
            combined.terms.keys().collect::<Vec<_>>(),
            cograph.terms.keys().collect::<Vec<_>>(),
            "tensoring an uncut UV spectator must preserve every raised cograph residue order"
        );
        assert_eq!(
            combined
                .terms
                .keys()
                .map(|index| index.lu_cut_order)
                .collect::<Vec<_>>(),
            vec![Some(1), Some(2)]
        );

        let exact_sum = |cff: &CutCFF, index: &CutCFFIndex| -> Result<Atom> {
            let term = cff
                .terms
                .get(index)
                .expect("the requested exact-CFF residue sector exists");
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
        let ordinary_sum = |graph: &Graph, cff: &CutCFF, index: &CutCFFIndex, numerator: &Atom| {
            cff.terms
                .get(index)
                .expect("the requested ordinary-CFF residue sector exists")
                .orientations
                .iter()
                .map(|orientation| {
                    orientation.expression.clone()
                        * numerator
                            .replace_multiple(orientation.orientation.energy_replacements_gs(graph))
                })
                .fold(Atom::Zero, |sum, term| sum + term)
                * Atom::num(cff.production_prefactor_factor())
        };
        let uv_index = CutCFFIndex::new_all_none();
        let uv_sum = exact_sum(&uv, &uv_index)?;
        let graph_mass_squared = (0..graph.underlying.n_edges())
            .map(|edge| {
                graph.underlying[EdgeIndex(edge)]
                    .particle
                    .mass_atom()
                    .pow(2)
            })
            .collect::<Vec<_>>();
        let external_edges = graph
            .external_momentum_edge_order()
            .into_iter()
            .collect::<BTreeSet<_>>();
        let representative_esurface_id = lu_cutset
            .residue_selector
            .lu
            .as_ref()
            .expect("the raised LU selector exists")
            .raised_group
            .esurface_ids[0];
        let representative_esurface =
            &production.expression.surfaces.esurface_cache[representative_esurface_id];
        let rescale_star = Atom::one();
        let cograph_spatial_momentum = Atom::num(Rational::from((3, 4)));
        let uv_spatial_momentum = Atom::num(Rational::from((4, 3)));
        let base_spatial_momentum = |edge: EdgeIndex| {
            if uv_edges.contains(&edge) {
                uv_spatial_momentum.clone()
            } else {
                cograph_spatial_momentum.clone()
            }
        };
        let scaled_on_shell_energy = |edge: EdgeIndex, rescale: Atom| {
            (graph_mass_squared[usize::from(edge)].clone()
                + (base_spatial_momentum(edge) * rescale).pow(2))
            .sqrt()
        };
        let root_energy_sum = representative_esurface
            .energies
            .iter()
            .map(|edge| scaled_on_shell_energy(*edge, rescale_star.clone()))
            .fold(Atom::Zero, |sum, energy| sum + energy);
        let external_shift_coefficient = representative_esurface
            .external_shift
            .iter()
            .map(|(_, coefficient)| *coefficient)
            .sum::<i64>();
        assert_ne!(
            external_shift_coefficient, 0,
            "the raised LU surface must have a nonzero external shift"
        );
        let external_energy = -root_energy_sum / Atom::num(external_shift_coefficient);
        let rescale_expression = |mut expression: Atom| {
            let rescale = Atom::var(GS.rescale);
            expression = expression
                .replace(Atom::var(GS.m_uv_expansion))
                .with(Atom::one());
            for edge in 0..graph_mass_squared.len() {
                let edge = EdgeIndex(edge);
                if external_edges.contains(&edge) {
                    expression = expression
                        .replace(GS.emr_mom(edge, GS.cind(0)))
                        .with(external_energy.clone());
                    for spatial_index in 1..=3 {
                        expression = expression
                            .replace(GS.emr_mom(edge, GS.cind(spatial_index)))
                            .with(Atom::Zero);
                    }
                } else {
                    let on_shell_energy = scaled_on_shell_energy(edge, rescale.clone());
                    expression = expression
                        .replace(GS.emr_mom(edge, GS.cind(1)))
                        .with(base_spatial_momentum(edge) * &rescale)
                        .replace(GS.ose(edge))
                        .with(on_shell_energy.clone())
                        .replace(cut_energy(edge))
                        .with(on_shell_energy);
                    for spatial_index in 2..=3 {
                        expression = expression
                            .replace(GS.emr_mom(edge, GS.cind(spatial_index)))
                            .with(Atom::Zero);
                    }
                }
            }
            expression
        };
        let eta = rescale_expression(representative_esurface.to_atom(&[]));
        assert!(
            eta.replace(GS.rescale)
                .with(rescale_star.clone())
                .expand()
                .is_zero(),
            "the raised-LU spectator oracle must be evaluated at its selected radial root"
        );
        let value_and_t_derivative = |expression: Atom| -> Result<[Atom; 2]> {
            let series = rescale_expression(expression)
                .series(GS.rescale, rescale_star.clone(), 1)
                .map_err(|error| eyre::eyre!("failed to build raised-LU t jet: {error}"))?;
            Ok([
                series.coefficient(Rational::from(0)),
                series.coefficient(Rational::from(1)),
            ])
        };
        let evaluate_arb = |expression: Atom| -> Result<Complex<F<ArbPrec>>> {
            let parameters = [Atom::var(GS.pi)];
            let rational: ExpressionEvaluator<SymComplex<Fraction<IntegerRing>>> =
                expression.evaluator(&parameters).build().map_err(|error| {
                    eyre::eyre!("failed to build raised-LU spectator evaluator: {error}")
                })?;
            let mut arb: ExpressionEvaluator<Complex<F<ArbPrec>>> =
                rational.map_coeff(&|coefficient| {
                    Complex::new(F::from(&coefficient.re), F::from(&coefficient.im))
                });
            Ok(arb.evaluate_single(&[Complex::new(
                F::<ArbPrec>::from_f64(0.0).pi(),
                F::<ArbPrec>::from_f64(0.0),
            )]))
        };
        // One eighth of ArbPrec's requested precision allows substantial and
        // construction-dependent bit loss while still scaling with precision.
        let tolerance = F(ArbPrec::default().epsilon()).sqrt().sqrt().sqrt();
        let mut failures = Vec::new();
        for index in cograph.terms.keys() {
            let combined_jet = value_and_t_derivative(exact_sum(&combined, index)?)?;
            let factorized_jet = value_and_t_derivative(exact_sum(&cograph, index)? * &uv_sum)?;
            for (jet_component, combined_expression, factorized_expression) in
                ["value", "first t derivative"]
                    .into_iter()
                    .zip(combined_jet)
                    .zip(factorized_jet)
                    .map(|((component, combined), factorized)| (component, combined, factorized))
            {
                let combined_value = evaluate_arb(combined_expression)?;
                let factorized_value = evaluate_arb(factorized_expression)?;
                let distance = (combined_value.clone() - factorized_value.clone())
                    .norm()
                    .re;
                let combined_norm = combined_value.clone().norm().re;
                let factorized_norm = factorized_value.clone().norm().re;
                let scale = if combined_norm > factorized_norm {
                    combined_norm
                } else {
                    factorized_norm
                };
                let relative_distance = if scale.is_zero() {
                    distance
                } else {
                    distance / scale
                };
                if relative_distance > tolerance {
                    failures.push(format!(
                        "raised LU residue {index} {jet_component} does not factorize from its quadratic cubic UV spectator: combined={combined_value:e}, factorized={factorized_value:e}, relative delta={relative_distance:e}, tolerance={tolerance:e}"
                    ));
                }
            }
        }
        assert!(
            failures.is_empty(),
            "raised-LU cubic spectator factorization failures:\n{}",
            failures.join("\n")
        );

        // The self-energy Taylor source contains a cubic dotted term. Isolate
        // its numerator-denominator pinch without cancelling D(Q) before CFF
        // generation: the retained factor raises one aliased outer energy by
        // two and leaves the other two outer energies at rank one, matching
        // the GL0 assignment induced by Q_UV slash (Q_UV dot Q_outer).
        assert_eq!(
            graph.loop_momentum_basis.edge_signatures[EdgeIndex(1)],
            graph.loop_momentum_basis.edge_signatures[EdgeIndex(2)],
            "the outer cograph must contain two owner-distinct edges with the same k routing"
        );
        let outer_numerator = [EdgeIndex(1), EdgeIndex(2), EdgeIndex(3)]
            .into_iter()
            .map(|edge| GS.emr_mom(edge, GS.cind(0)) + Atom::one())
            .fold(Atom::one(), |product, factor| product * factor);
        let retained_dotted_factor =
            Atom::num(2) * GS.emr_mom(cograph_edge, GS.cind(0)) * &outer_numerator;
        let (exact_outer, _) = graph.cff_from_4d_denominators(
            &cograph_denominators,
            &lu_cutset,
            &options,
            &retained_dotted_factor,
        )?;
        let mut uv_contract: SuBitGraph = graph.empty_subgraph();
        for edge in uv_edges {
            uv_contract.add(graph[&edge].1);
        }
        let ordinary_outer = graph.cff(
            &uv_contract,
            &lu_cutset,
            &OrientationPattern::default(),
            &options,
            Some(&retained_dotted_factor),
        )?;
        assert_eq!(
            exact_outer.terms.keys().collect::<Vec<_>>(),
            ordinary_outer.terms.keys().collect::<Vec<_>>(),
            "exact and ordinary alias cographs must expose both raised LU orders"
        );
        failures.clear();
        for index in ordinary_outer.terms.keys() {
            let exact_jet = value_and_t_derivative(exact_sum(&exact_outer, index)?)?;
            let ordinary_term = ordinary_outer
                .terms
                .get(index)
                .expect("the ordinary alias-cograph residue exists");
            let ordinary_sum = ordinary_term
                .orientations
                .iter()
                .map(|orientation| {
                    orientation.expression.clone()
                        * retained_dotted_factor.replace_multiple(
                            orientation.orientation.energy_replacements_gs(&graph),
                        )
                })
                .fold(Atom::Zero, |sum, term| sum + term)
                * Atom::num(ordinary_outer.production_prefactor_factor());
            let ordinary_jet = value_and_t_derivative(ordinary_sum)?;
            for (jet_component, exact_expression, ordinary_expression) in
                ["value", "first t derivative"]
                    .into_iter()
                    .zip(exact_jet)
                    .zip(ordinary_jet)
                    .map(|((component, exact), ordinary)| (component, exact, ordinary))
            {
                let exact_value = evaluate_arb(exact_expression)?;
                let ordinary_value = evaluate_arb(ordinary_expression)?;
                let distance = (exact_value.clone() - ordinary_value.clone()).norm().re;
                let exact_norm = exact_value.clone().norm().re;
                let ordinary_norm = ordinary_value.clone().norm().re;
                let scale = if exact_norm > ordinary_norm {
                    exact_norm
                } else {
                    ordinary_norm
                };
                let relative_distance = if scale.is_zero() {
                    distance
                } else {
                    distance / scale
                };
                if relative_distance > tolerance {
                    failures.push(format!(
                        "raised alias residue {index} {jet_component} differs between exact and ordinary cograph CFF: exact={exact_value:e}, ordinary={ordinary_value:e}, relative delta={relative_distance:e}, tolerance={tolerance:e}"
                    ));
                }
            }
        }
        assert!(
            failures.is_empty(),
            "raised-alias exact-vs-ordinary failures:\n{}",
            failures.join("\n")
        );
        assert_eq!(
            exact_outer.production_prefactor_factor(),
            ordinary_outer.production_prefactor_factor(),
            "exact and ordinary alias cographs must use the same contour convention"
        );

        let denominator_atom = GS.den(
            usize::from(uv_edge),
            uv_momentum.clone(),
            uv_mass_squared.clone(),
            uv_full_expr.clone(),
        );
        let dotted_numerator = &denominator_atom * &retained_dotted_factor;
        // Keep the same component order as the UV-expanded GL0 sources: the
        // UV component precedes the three-denominator outer cograph.
        let lower_denominators = uv_denominators[..2]
            .iter()
            .cloned()
            .chain(cograph_denominators.iter().cloned())
            .collect::<Vec<_>>();
        let dotted_denominators = uv_denominators
            .iter()
            .cloned()
            .chain(cograph_denominators.iter().cloned())
            .collect::<Vec<_>>();

        // The linear self-energy source is odd in the UV loop momentum. Its
        // temporal residue cancels between the two UV energy poles, whereas a
        // spatial component survives the energy residue pointwise. Test both
        // against the factorized exact source and the ordinary production CFF
        // so that neither a shared exact-source continuation defect nor a
        // spuriously vanishing temporal reference can hide a parity error.
        let mut cograph_contract: SuBitGraph = graph.empty_subgraph();
        for edge in [EdgeIndex(1), EdgeIndex(2), EdgeIndex(3)] {
            cograph_contract.add(graph[&edge].1);
        }
        failures.clear();
        for (uv_component, odd_uv_numerator) in [
            ("temporal", GS.emr_mom(uv_edge, GS.cind(0))),
            ("spatial", GS.emr_mom(uv_edge, GS.cind(1))),
        ] {
            let odd_combined_numerator = &retained_dotted_factor * &odd_uv_numerator;
            let (odd_uv_source, _) = graph.cff_from_4d_denominators_in_uv_edges(
                &uv_denominators[..2],
                uv_edges,
                &CutSet::empty(graph.n_hedges()),
                &options,
                &odd_uv_numerator,
                None,
            )?;
            let (odd_combined_source, _) = graph.cff_from_4d_denominators_in_uv_edges(
                &lower_denominators,
                uv_edges,
                &lu_cutset,
                &options,
                &odd_combined_numerator,
                None,
            )?;
            let ordinary_uv = graph.cff(
                &cograph_contract,
                &CutSet::empty(graph.n_hedges()),
                &OrientationPattern::default(),
                &options,
                Some(&odd_uv_numerator),
            )?;
            let exact_uv_sum = exact_sum(&odd_uv_source, &uv_index)?;
            let ordinary_uv_sum = ordinary_sum(&graph, &ordinary_uv, &uv_index, &odd_uv_numerator);
            for index in exact_outer.terms.keys() {
                let combined_jet = value_and_t_derivative(exact_sum(&odd_combined_source, index)?)?;
                let exact_factorized_jet =
                    value_and_t_derivative(exact_sum(&exact_outer, index)? * &exact_uv_sum)?;
                let ordinary_factorized_jet = value_and_t_derivative(
                    ordinary_sum(&graph, &ordinary_outer, index, &retained_dotted_factor)
                        * &ordinary_uv_sum,
                )?;
                for (jet_component, combined_expression, exact_expression, ordinary_expression) in
                    ["value", "first t derivative"]
                        .into_iter()
                        .zip(combined_jet)
                        .zip(exact_factorized_jet)
                        .zip(ordinary_factorized_jet)
                        .map(|(((component, combined), exact), ordinary)| {
                            (component, combined, exact, ordinary)
                        })
                {
                    let combined_value = evaluate_arb(combined_expression)?;
                    for (reference, reference_expression) in [
                        ("factorized exact", exact_expression),
                        ("factorized ordinary", ordinary_expression),
                    ] {
                        let reference_value = evaluate_arb(reference_expression)?;
                        let distance = (combined_value.clone() - reference_value.clone()).norm().re;
                        let combined_norm = combined_value.clone().norm().re;
                        let reference_norm = reference_value.clone().norm().re;
                        let scale = if combined_norm > reference_norm {
                            combined_norm
                        } else {
                            reference_norm
                        };
                        let relative_distance = if scale.is_zero() {
                            distance
                        } else {
                            distance / scale
                        };
                        if relative_distance > tolerance {
                            failures.push(format!(
                                "raised LU residue {index} {jet_component} differs from its {reference} odd-UV-{uv_component} source: combined={combined_value:e}, reference={reference_value:e}, relative delta={relative_distance:e}, tolerance={tolerance:e}"
                            ));
                        }
                    }
                }
            }
        }
        assert!(
            failures.is_empty(),
            "raised-LU odd-UV-momentum factorization failures:\n{}",
            failures.join("\n")
        );

        let (retained_source, _) = graph.cff_from_4d_denominators_in_uv_edges(
            &lower_denominators,
            uv_edges,
            &lu_cutset,
            &options,
            &retained_dotted_factor,
            None,
        )?;
        let (dotted_source, _) = graph.cff_from_4d_denominators_in_uv_edges(
            &dotted_denominators,
            uv_edges,
            &lu_cutset,
            &options,
            &dotted_numerator,
            None,
        )?;
        assert_eq!(
            retained_source.terms.keys().collect::<Vec<_>>(),
            dotted_source.terms.keys().collect::<Vec<_>>(),
            "powered and lower exact sources must expose the same raised LU orders"
        );
        failures.clear();
        for index in retained_source.terms.keys() {
            let dotted_sum = exact_sum(&dotted_source, index)?
                .replace(GS.den(W_.a_, W_.b_, W_.c_, W_.d_))
                .with(W_.d_);
            let retained_sum = exact_sum(&retained_source, index)?;
            let dotted_jet = value_and_t_derivative(dotted_sum)?;
            let retained_jet = value_and_t_derivative(retained_sum)?;
            for (jet_component, dotted_expression, retained_expression) in
                ["value", "first t derivative"]
                    .into_iter()
                    .zip(dotted_jet)
                    .zip(retained_jet)
                    .map(|((component, dotted), retained)| (component, dotted, retained))
            {
                let dotted_value = evaluate_arb(dotted_expression)?;
                let retained_value = evaluate_arb(retained_expression)?;
                let distance = (dotted_value.clone() - retained_value.clone()).norm().re;
                let dotted_norm = dotted_value.clone().norm().re;
                let retained_norm = retained_value.clone().norm().re;
                let scale = if dotted_norm > retained_norm {
                    dotted_norm
                } else {
                    retained_norm
                };
                let relative_distance = if scale.is_zero() {
                    distance
                } else {
                    distance / scale
                };
                if relative_distance > tolerance {
                    failures.push(format!(
                        "raised LU residue {index} {jet_component} differs between D_UV*F/D_UV^3 and F/D_UV^2: dotted={dotted_value:e}, lower={retained_value:e}, relative delta={relative_distance:e}, tolerance={tolerance:e}"
                    ));
                }
            }
        }
        assert!(
            failures.is_empty(),
            "raised-LU powered-alias identity failures:\n{}",
            failures.join("\n")
        );

        // Finally retain the genuinely cross-component numerator of the GL0
        // cubic source. Its temporal and one spatial Minkowski component must
        // equal the corresponding sum of factorized cograph and UV residues;
        // testing the t jet makes the order-two Cutkosky derivative part of
        // the oracle rather than checking only its value at the cut.
        let cograph_temporal_numerator = GS.emr_mom(cograph_edge, GS.cind(0)) * &outer_numerator;
        let cograph_spatial_numerator = GS.emr_mom(cograph_edge, GS.cind(1)) * &outer_numerator;
        let uv_temporal_numerator = GS.emr_mom(uv_edge, GS.cind(0)).pow(2);
        let uv_spatial_numerator =
            GS.emr_mom(uv_edge, GS.cind(0)) * GS.emr_mom(uv_edge, GS.cind(1));
        let cross_component_numerator = &cograph_temporal_numerator * &uv_temporal_numerator
            - &cograph_spatial_numerator * &uv_spatial_numerator;
        let (cross_component_source, _) = graph.cff_from_4d_denominators_in_uv_edges(
            &dotted_denominators,
            uv_edges,
            &lu_cutset,
            &options,
            &cross_component_numerator,
            None,
        )?;
        let (cograph_temporal, _) = graph.cff_from_4d_denominators(
            &cograph_denominators,
            &lu_cutset,
            &options,
            &cograph_temporal_numerator,
        )?;
        let (cograph_spatial, _) = graph.cff_from_4d_denominators(
            &cograph_denominators,
            &lu_cutset,
            &options,
            &cograph_spatial_numerator,
        )?;
        let (uv_temporal, _) = graph.cff_from_4d_denominators_in_uv_edges(
            &uv_denominators,
            uv_edges,
            &CutSet::empty(graph.n_hedges()),
            &options,
            &uv_temporal_numerator,
            None,
        )?;
        let (uv_spatial, _) = graph.cff_from_4d_denominators_in_uv_edges(
            &uv_denominators,
            uv_edges,
            &CutSet::empty(graph.n_hedges()),
            &options,
            &uv_spatial_numerator,
            None,
        )?;
        let uv_temporal_sum = exact_sum(&uv_temporal, &uv_index)?;
        let uv_spatial_sum = exact_sum(&uv_spatial, &uv_index)?;
        assert_eq!(
            cross_component_source.terms.keys().collect::<Vec<_>>(),
            cograph_temporal.terms.keys().collect::<Vec<_>>(),
            "the cross-component source must preserve both raised cograph orders"
        );
        failures.clear();
        for index in cograph_temporal.terms.keys() {
            let combined_jet = value_and_t_derivative(exact_sum(&cross_component_source, index)?)?;
            let factorized_jet = value_and_t_derivative(
                exact_sum(&cograph_temporal, index)? * &uv_temporal_sum
                    - exact_sum(&cograph_spatial, index)? * &uv_spatial_sum,
            )?;
            for (jet_component, combined_expression, factorized_expression) in
                ["value", "first t derivative"]
                    .into_iter()
                    .zip(combined_jet)
                    .zip(factorized_jet)
                    .map(|((component, combined), factorized)| (component, combined, factorized))
            {
                let combined_value = evaluate_arb(combined_expression)?;
                let factorized_value = evaluate_arb(factorized_expression)?;
                let distance = (combined_value.clone() - factorized_value.clone())
                    .norm()
                    .re;
                let combined_norm = combined_value.clone().norm().re;
                let factorized_norm = factorized_value.clone().norm().re;
                let scale = if combined_norm > factorized_norm {
                    combined_norm
                } else {
                    factorized_norm
                };
                let relative_distance = if scale.is_zero() {
                    distance
                } else {
                    distance / scale
                };
                if relative_distance > tolerance {
                    failures.push(format!(
                        "raised LU residue {index} {jet_component} does not factorize for Q_UV^0 (Q_UV^0 Q_outer^0-Q_UV^1 Q_outer^1): combined={combined_value:e}, factorized={factorized_value:e}, relative delta={relative_distance:e}, tolerance={tolerance:e}"
                    ));
                }
            }
        }
        assert!(
            failures.is_empty(),
            "raised-LU cross-component factorization failures:\n{}",
            failures.join("\n")
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
    fn connected_parent_uv_limit_matches_uncorrected_child_times_outer() -> Result<()> {
        test_initialise()?;
        let check = |label: &str,
                     mut graph: Graph,
                     expected_bridges: (i64, i64, i64),
                     expected_loop_count: usize,
                     selected_lu_support: Option<&[EdgeIndex]>|
         -> Result<()> {
            let uv_edges = [EdgeIndex(1), EdgeIndex(2)];
            let mut uv_contract: SuBitGraph = graph.empty_subgraph();
            for edge in uv_edges {
                uv_contract.add(graph[&edge].1);
            }
            let uv_subgraph =
                InternalSubGraph::cleaned_filter_optimist(uv_contract.clone(), graph.as_ref());
            assert_eq!(graph.get_loop_number(), expected_loop_count);
            assert_eq!(graph.n_loops(&uv_contract), 1);

            let options = graph.denominator_only_cff_3d_expression_options();
            let cutset = CutSet::empty(graph.n_hedges());
            let pattern = OrientationPattern::default();
            let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
            let root_generated =
                graph.generate_3d_expression_for_integrand(&[], &canonization, &options, None)?;
            let outer_generated = graph.generate_3d_expression_for_integrand(
                &uv_edges,
                &canonization,
                &options,
                None,
            )?;
            assert_eq!(
                (
                    root_generated.energy_factor_ownership,
                    outer_generated.energy_factor_ownership,
                ),
                (
                    CffEnergyFactorOwnership::GlobalSourceProduct,
                    CffEnergyFactorOwnership::GlobalSourceProduct,
                ),
                "{label}: scalar root and reduced outer must retain global energy-factor ownership",
            );
            let root_contract = graph.empty_subgraph::<SuBitGraph>();
            let root = graph.cff(&root_contract, &cutset, &pattern, &options, None)?;
            let outer = graph.cff(&uv_contract, &cutset, &pattern, &options, None)?;

            let crown = graph.dummy_stripped_external_flows_of(&uv_subgraph);
            let child_lmb = graph.try_compatible_sub_lmb(
                &uv_subgraph,
                crown.clone(),
                &graph.loop_momentum_basis,
            )?;
            assert_eq!(child_lmb.loop_edges.len(), 1);
            let muv_mass_squared = Atom::var(GS.m_uv_expansion).pow(2);
            let child_denominators = uv_edges.map(|edge| FourDDenominator {
                source_edge: edge,
                momentum: FunctionBuilder::new(GS.emr_mom)
                    .add_arg(usize::from(edge))
                    .finish(),
                mass_squared: muv_mass_squared.clone(),
                full_expr: Atom::one(),
            });
            {
                let source = GraphThreeDSource::from_exact_denominators_in_uv_sub_lmb(
                    &graph,
                    &child_denominators,
                    uv_edges,
                    crown.included_iter(),
                    &child_lmb,
                    ExactUvSubLmbFrame::RetainedPhysicalCrown,
                )?;
                let (child_generated, _, _, _) = graph.generate_3d_expression_for_4d_term(
                    &source,
                    &options,
                    &Atom::one(),
                    None,
                )?;
                assert_eq!(
                    child_generated.energy_factor_ownership,
                    CffEnergyFactorOwnership::GlobalSourceProduct,
                    "{label}: the scalar retained-crown child must retain global energy-factor ownership",
                );
            }
            let (child, _) = graph.cff_from_4d_denominators_in_uv_sub_lmb(
                &child_denominators,
                uv_edges,
                crown.included_iter(),
                &child_lmb,
                ExactUvSubLmbFrame::RetainedPhysicalCrown,
                &cutset,
                &options,
                &Atom::one(),
                None,
            )?;

            let full_subgraph =
                InternalSubGraph::cleaned_filter_optimist(graph.full_filter(), graph.as_ref());
            let outer_lmb = graph.shrunken_sub_lmb(
                &full_subgraph.filter,
                &uv_subgraph,
                graph.dummy_stripped_external_flows_of(&full_subgraph),
                Some(&graph.loop_momentum_basis),
            )?;
            let parent_carriers = graph
                .loop_momentum_basis
                .loop_edges
                .iter()
                .copied()
                .collect::<Vec<_>>();
            let appended_carriers = child_lmb
                .loop_edges
                .iter()
                .chain(&outer_lmb.loop_edges)
                .copied()
                .collect::<Vec<_>>();
            let carrier_permutation = appended_carriers
                .iter()
                .map(|edge| {
                    parent_carriers
                        .iter()
                        .position(|parent_edge| parent_edge == edge)
                        .expect("parent-compatible child and outer carriers belong to the root LMB")
                })
                .collect::<Vec<_>>();
            let inversion_count = carrier_permutation
                .iter()
                .enumerate()
                .map(|(left, value)| {
                    carrier_permutation[left + 1..]
                        .iter()
                        .filter(|right| value > *right)
                        .count()
                })
                .sum::<usize>();
            assert_eq!(
                appended_carriers.len(),
                parent_carriers.len(),
                "{label}: child and reduced-outer carriers must span the root coordinates",
            );
            let append_orientation = if inversion_count.is_multiple_of(2) {
                1
            } else {
                -1
            };

            let ordinary_raw_sum = |cff: &CutCFF| {
                cff.terms
                    .values()
                    .flat_map(|term| &term.orientations)
                    .fold(Atom::Zero, |sum, orientation| sum + &orientation.expression)
            };
            let exact_raw_sum = |cff: &CutCFF| -> Result<Atom> {
                Ok(cff
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
                    .fold(Atom::Zero, |sum, term| sum + term))
            };
            let root_bridge = root.production_prefactor_factor();
            let outer_bridge = outer.production_prefactor_factor();
            let child_bridge = child.production_prefactor_factor();
            let parent_raw = ordinary_raw_sum(&root);
            let components_raw = ordinary_raw_sum(&outer) * exact_raw_sum(&child)?;

            assert_eq!(
                (root_bridge, outer_bridge, child_bridge),
                expected_bridges,
                "{label}: the connected parent, reduced outer, and retained-crown child must expose their own standalone core bridges",
            );
            // The one-loop outer case has the parity opposite to GL24, while the
            // two-loop outer case matches it. In both, the relative bridge converts
            // separately generated raw child-times-outer data into the raw
            // connected-parent contour frame. It is not also the production host
            // factor: production consumes the root bridge on that connected frame,
            // while the separate components consume their own bridges.
            let relative_bridge = root_bridge * outer_bridge * child_bridge;
            assert_eq!(relative_bridge, -1);
            if outer_bridge * child_bridge == 1 {
                assert_eq!(
                    relative_bridge, root_bridge,
                    "{label}: with unit child and outer bridges, the connected-frame conversion applies the root bridge exactly once"
                );
            }
            let connected_frame = components_raw.clone() * Atom::num(relative_bridge);
            let parent = parent_raw.clone() * Atom::num(root_bridge);
            let factorized = components_raw.clone() * Atom::num(outer_bridge * child_bridge);

            let zero = F::<ArbPrec>::from_f64(0.0);
            let one = F::<ArbPrec>::from_f64(1.0);
            let muv_mass = F::<ArbPrec>::from(&Rational::from((7, 5)));
            let outer_spatial = F::<ArbPrec>::from(&Rational::from((3, 7)));
            let extra_outer_spatial = F::<ArbPrec>::from(&Rational::from((2, 5)));
            let mut parameters = vec![Atom::var(GS.pi), Atom::var(GS.m_uv_expansion)];
            let mut edge_parameters = Vec::new();
            for edge in 0..graph.underlying.n_edges() {
                let edge = EdgeIndex(edge);
                parameters.extend([GS.ose(edge), cut_energy(edge), GS.emr_mom(edge, GS.cind(0))]);
                for component in 1..=3 {
                    parameters.push(GS.emr_mom(edge, GS.cind(component)));
                }
                edge_parameters.push(edge);
            }
            let values_at = |scale: &F<ArbPrec>| {
                let child_spatial = [scale.clone(), outer_spatial.clone() - scale];
                let mut values = vec![
                    Complex::new(zero.clone().pi(), zero.clone()),
                    Complex::new_re(muv_mass.clone()),
                ];
                for edge in &edge_parameters {
                    let child_index = uv_edges.iter().position(|candidate| candidate == edge);
                    let (energy, spatial) = if let Some(child_index) = child_index {
                        let spatial = child_spatial[child_index].clone();
                        (
                            (&muv_mass * &muv_mass + &spatial * &spatial).sqrt(),
                            spatial,
                        )
                    } else if let Some(spatial) = match (expected_loop_count, *edge) {
                        (_, EdgeIndex(3)) => Some(outer_spatial.clone()),
                        (2, EdgeIndex(4)) => Some(-&outer_spatial),
                        (3, EdgeIndex(4)) => Some(&extra_outer_spatial - &outer_spatial),
                        (3, EdgeIndex(5)) => Some(extra_outer_spatial.clone()),
                        (3, EdgeIndex(6)) => Some(-&extra_outer_spatial),
                        _ => None,
                    } {
                        ((&one + &spatial * &spatial).sqrt(), spatial)
                    } else {
                        (F::<ArbPrec>::from_f64(2.0), zero.clone())
                    };
                    values.extend([
                        Complex::new_re(energy.clone()),
                        Complex::new_re(energy.clone()),
                        Complex::new_re(energy),
                    ]);
                    for component in 1..=3 {
                        values.push(Complex::new_re(if component == 1 {
                            spatial.clone()
                        } else {
                            zero.clone()
                        }));
                    }
                }
                values
            };
            let evaluate_arb =
                |expression: Atom, values: &[Complex<F<ArbPrec>>]| -> Result<Complex<F<ArbPrec>>> {
                    let rational: ExpressionEvaluator<SymComplex<Fraction<IntegerRing>>> =
                    expression.evaluator(&parameters).build().map_err(|error| {
                        eyre::eyre!(
                            "failed to build {label} connected-parent UV-bridge evaluator: {error}"
                        )
                    })?;
                    let mut arb: ExpressionEvaluator<Complex<F<ArbPrec>>> =
                        rational.map_coeff(&|coefficient| {
                            Complex::new(F::from(&coefficient.re), F::from(&coefficient.im))
                        });
                    Ok(arb.evaluate_single(values))
                };
            let scale = F::<ArbPrec>::from_f64(1.0e45);
            let uv_values = values_at(&scale);
            let measure = &scale * &scale * &scale;
            let parent_raw = evaluate_arb(parent_raw, &uv_values)? * &measure;
            let connected_frame = evaluate_arb(connected_frame, &uv_values)? * &measure;
            let parent = evaluate_arb(parent, &uv_values)? * &measure;
            let factorized = evaluate_arb(factorized, &uv_values)? * &measure;
            let raw_distance = (parent_raw.clone() - connected_frame.clone()).norm().re;
            let parent_raw_norm = parent_raw.clone().norm().re;
            let connected_frame_norm = connected_frame.clone().norm().re;
            let raw_norm = if parent_raw_norm > connected_frame_norm {
                parent_raw_norm
            } else {
                connected_frame_norm
            };
            let raw_relative_distance = if raw_norm.is_zero() {
                raw_distance
            } else {
                raw_distance / raw_norm
            };
            let distance = (parent.clone() - factorized.clone()).norm().re;
            let parent_norm = parent.clone().norm().re;
            let factorized_norm = factorized.clone().norm().re;
            let norm = if parent_norm > factorized_norm {
                parent_norm
            } else {
                factorized_norm
            };
            let relative_distance = if norm.is_zero() {
                distance
            } else {
                distance / norm.clone()
            };
            let tolerance = F(ArbPrec::default().epsilon()).sqrt().sqrt().sqrt();
            assert!(
                raw_relative_distance <= tolerance,
                "{label}: raw connected-parent UV frame differs from raw child times outer after one relative/root bridge: parent={parent_raw:e}, connected frame={connected_frame:e}, bridge={relative_bridge}, append orientation={append_orientation}, relative delta={raw_relative_distance:e}, tolerance={tolerance:e}",
            );
            assert!(
                relative_distance <= tolerance,
                "{label}: connected parent UV limit differs from its retained-crown exact child times reduced outer: parent={parent:e}, factorized={factorized:e}, standalone bridge ratio={relative_bridge}, relative delta={relative_distance:e}, tolerance={tolerance:e}",
            );
            let naively_bridged =
                factorized.clone() * F::<ArbPrec>::from_f64(relative_bridge as f64);
            let naive_relative_distance = (parent.clone() - naively_bridged).norm().re / norm;
            assert!(
                naive_relative_distance > F::<ArbPrec>::from_f64(1.0),
                "{label}: the standalone core-sign ratio unexpectedly behaved like an extra connected-parent gluing factor",
            );

            if let Some(lu_support) = selected_lu_support {
                let quartic_owner = EdgeIndex(3);
                let outer_momentum_squared = GS.emr_mom(quartic_owner, GS.cind(0)).pow(2)
                    - (1..=3)
                        .map(|component| GS.emr_mom(quartic_owner, GS.cind(component)).pow(2))
                        .fold(Atom::Zero, |sum, component| sum + component);
                let quartic_numerator = outer_momentum_squared.pow(2);
                // Keep the temporal quartic above only to exercise the same VariantLocal
                // construction as GL24.  At a fixed selected residue, GL24's source map has
                // already reduced that numerator to this outer spatial quartic; evaluating it
                // directly here avoids turning this factorization check into a second residue
                // map whose temporal sample would contain the child UV pole.
                let fixed_selected_spatial_quartic = (1..=3)
                    .map(|component| GS.emr_mom(quartic_owner, GS.cind(component)).pow(2))
                    .fold(Atom::Zero, |sum, component| sum - component)
                    .pow(2);
                let root_generated = graph.generate_3d_expression_for_integrand(
                    &[],
                    &canonization,
                    &options,
                    Some(&quartic_numerator),
                )?;
                let outer_generated = graph.generate_3d_expression_for_integrand(
                    &uv_edges,
                    &canonization,
                    &options,
                    Some(&quartic_numerator),
                )?;
                assert_eq!(
                    (
                        root_generated.energy_factor_ownership,
                        outer_generated.energy_factor_ownership,
                    ),
                    (
                        CffEnergyFactorOwnership::VariantLocal,
                        CffEnergyFactorOwnership::VariantLocal,
                    ),
                    "{label}: quartic connected root and selected reduced outer must use variant-local energy factors",
                );
                let target_support = lu_support.iter().copied().collect::<BTreeSet<_>>();
                let lu_surface = graph
                    .determine_raised_esurfaces_from_expression(&root_generated.expression)
                    .raised_groups
                    .into_iter()
                    .find(|group| {
                        group.esurface_ids.iter().any(|surface_id| {
                            root_generated.expression.surfaces.esurface_cache[*surface_id]
                                .energies
                                .iter()
                                .copied()
                                .collect::<BTreeSet<_>>()
                                .eq(&target_support)
                        })
                    })
                    .ok_or_else(|| {
                        eyre::eyre!(
                            "{label}: no selected LU surface has physical support {target_support:?}"
                        )
                    })?;
                let mut selected_cutset = CutSet::empty(graph.n_hedges());
                selected_cutset.residue_selector.lu = Some(crate::graph::cuts::LuCutSelection {
                    raised_group: lu_surface,
                    cut_edge_alternatives: vec![lu_support.to_vec()],
                });
                let variant_root = graph.cff(
                    &root_contract,
                    &selected_cutset,
                    &pattern,
                    &options,
                    Some(&quartic_numerator),
                )?;
                let variant_outer = graph.cff(
                    &uv_contract,
                    &selected_cutset,
                    &pattern,
                    &options,
                    Some(&quartic_numerator),
                )?;
                let fixed_order_one_raw_sum = |cff: &CutCFF| {
                    cff.terms
                        .iter()
                        .filter(|(index, _)| index.lu_cut_order == Some(1))
                        .flat_map(|(_, term)| {
                            term.orientations.iter().map(|orientation| {
                                orientation.expression.clone()
                                    * fixed_selected_spatial_quartic.clone()
                            })
                        })
                        .fold(Atom::Zero, |sum, term| sum + term)
                };
                let variant_parent_raw = fixed_order_one_raw_sum(&variant_root);
                let variant_components_raw =
                    fixed_order_one_raw_sum(&variant_outer) * exact_raw_sum(&child)?;
                let variant_connected_frame =
                    variant_components_raw.clone() * Atom::num(relative_bridge);
                let variant_parent = variant_parent_raw.clone() * Atom::num(root_bridge);
                let variant_factorized =
                    variant_components_raw * Atom::num(outer_bridge * child_bridge);
                let variant_parent_raw = evaluate_arb(variant_parent_raw, &uv_values)? * &measure;
                let variant_connected_frame =
                    evaluate_arb(variant_connected_frame, &uv_values)? * &measure;
                let variant_parent = evaluate_arb(variant_parent, &uv_values)? * &measure;
                let variant_factorized = evaluate_arb(variant_factorized, &uv_values)? * &measure;
                let relative_delta = |left: &Complex<F<ArbPrec>>, right: &Complex<F<ArbPrec>>| {
                    let distance = (left.clone() - right.clone()).norm().re;
                    let left_norm = left.clone().norm().re;
                    let right_norm = right.clone().norm().re;
                    let scale = if left_norm > right_norm {
                        left_norm
                    } else {
                        right_norm
                    };
                    if scale.is_zero() {
                        distance
                    } else {
                        distance / scale
                    }
                };
                let raw_delta = relative_delta(&variant_parent_raw, &variant_connected_frame);
                let production_delta = relative_delta(&variant_parent, &variant_factorized);
                assert!(
                    raw_delta <= tolerance,
                    "{label}: selected-LU quartic raw connected frame differs from child times reduced outer: parent={variant_parent_raw:e}, connected frame={variant_connected_frame:e}, append orientation={append_orientation}, relative delta={raw_delta:e}, tolerance={tolerance:e}",
                );
                assert!(
                    production_delta <= tolerance,
                    "{label}: selected-LU quartic connected-parent UV limit differs from independently bridged child times reduced outer: parent={variant_parent:e}, factorized={variant_factorized:e}, append orientation={append_orientation}, relative delta={production_delta:e}, tolerance={tolerance:e}",
                );
            }
            Ok(())
        };

        check(
            "opposite outer parity",
            dot!(digraph connected_parent_uv_bridge {
                edge [num=1 mass=1]
                node [num=1]
                incoming [style=invis]
                outgoing [style=invis]

                incoming -> a [id=0]
                a -> u [id=1 lmb_id=0]
                a -> u [id=2]
                u -> b [id=3 lmb_id=1]
                a -> b [id=4]
                b -> outgoing [id=5]
            })?,
            (-1, 1, 1),
            2,
            None,
        )?;
        check(
            "GL24-matching outer parity",
            dot!(digraph connected_parent_uv_bridge_three_loop {
                edge [num=1 mass=1]
                node [num=1]
                incoming [style=invis]
                outgoing [style=invis]

                incoming -> a [id=0]
                a -> u [id=1 lmb_id=0]
                a -> u [id=2]
                u -> b [id=3 lmb_id=1]
                a -> b [id=4]
                b -> c [id=5 lmb_id=2]
                a -> c [id=6]
                c -> outgoing [id=7]
            })?,
            (1, -1, 1),
            3,
            Some(&[EdgeIndex(5), EdgeIndex(6)]),
        )?;
        Ok(())
    }

    #[test]
    fn gl24_selected_zero_sample_uses_the_reduced_outer_frame() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph gl24_selected_zero_sample {
            edge [num=1 mass=1]
            node [num=1]

            v0 -> v1 [id=0 is_cut=0]
            v0 -> v3 [id=1 lmb_id=0]
            v0 -> v5 [id=2]
            v1 -> v2 [id=3 lmb_id=1]
            v1 -> v3 [id=4]
            v2 -> v4 [id=5 lmb_id=2]
            v2 -> v4 [id=6]
            v3 -> v5 [id=7]
            v4 -> v5 [id=8 mass=2]
        })?;
        assert_eq!(graph.get_loop_number(), 3);
        assert_eq!(graph.get_edges_in_initial_state_cut(), [EdgeIndex(0)]);

        let uv_edges = [EdgeIndex(5), EdgeIndex(6)];
        let mut uv_contract: SuBitGraph = graph.empty_subgraph();
        for edge in uv_edges {
            uv_contract.add(graph[&edge].1);
        }
        let production_contract = graph.tree_edges.subtract(&graph.initial_state_cut);
        let production_contract_edges = graph.paired_edges(&production_contract);
        let outer_contract = uv_contract.union(&production_contract);
        let outer_contract_edges = graph.paired_edges(&outer_contract);
        let options = graph.denominator_only_cff_3d_expression_options();
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let quartic_owner = EdgeIndex(1);
        let momentum_squared = GS.emr_mom(quartic_owner, GS.cind(0)).pow(2)
            - (1..=3)
                .map(|component| GS.emr_mom(quartic_owner, GS.cind(component)).pow(2))
                .fold(Atom::Zero, |sum, component| sum + component);
        let quartic_numerator = momentum_squared.pow(2);
        let root_generated = graph.generate_3d_expression_for_integrand(
            &production_contract_edges,
            &canonization,
            &options,
            Some(&quartic_numerator),
        )?;
        let outer_generated = graph.generate_3d_expression_for_integrand(
            &outer_contract_edges,
            &canonization,
            &options,
            Some(&quartic_numerator),
        )?;
        assert_eq!(
            (
                root_generated.energy_factor_ownership,
                outer_generated.energy_factor_ownership,
            ),
            (
                CffEnergyFactorOwnership::VariantLocal,
                CffEnergyFactorOwnership::VariantLocal,
            ),
        );

        let target_support = [EdgeIndex(3), EdgeIndex(4)]
            .into_iter()
            .collect::<BTreeSet<_>>();
        let lu_surface = graph
            .determine_raised_esurfaces_from_expression(&root_generated.expression)
            .raised_groups
            .into_iter()
            .find(|group| {
                group.esurface_ids.iter().any(|surface_id| {
                    root_generated.expression.surfaces.esurface_cache[*surface_id]
                        .energies
                        .iter()
                        .copied()
                        .collect::<BTreeSet<_>>()
                        .eq(&target_support)
                })
            })
            .ok_or_else(|| eyre::eyre!("GL24 q3/q4 LU surface was not generated"))?;
        let mut cutset = CutSet::empty(graph.n_hedges());
        cutset.residue_selector.lu = Some(crate::graph::cuts::LuCutSelection {
            raised_group: lu_surface,
            cut_edge_alternatives: vec![vec![EdgeIndex(3), EdgeIndex(4)]],
        });
        let pattern = OrientationPattern::default();
        let root = graph.cff(
            &production_contract,
            &cutset,
            &pattern,
            &options,
            Some(&quartic_numerator),
        )?;
        let outer = graph.cff(
            &outer_contract,
            &cutset,
            &pattern,
            &options,
            Some(&quartic_numerator),
        )?;
        assert_eq!(
            (
                root.production_prefactor_factor(),
                outer.production_prefactor_factor(),
            ),
            (1, -1),
        );
        let zero_sample_terms = |cff: &CutCFF| {
            cff.terms
                .iter()
                .filter(|(index, _)| index.lu_cut_order == Some(1))
                .flat_map(|(_, term)| &term.orientations)
                .filter(|orientation| {
                    orientation.orientation.edge_energy_map[usize::from(quartic_owner)].is_zero()
                })
                .map(|orientation| {
                    (
                        orientation.expression.clone(),
                        orientation.orientation.clone(),
                    )
                })
                .collect::<Vec<_>>()
        };
        let root_zero = zero_sample_terms(&root);
        let outer_zero = zero_sample_terms(&outer);
        assert!(!root_zero.is_empty() && !outer_zero.is_empty());
        let root_target = EdgeVec::from_iter([
            Orientation::Undirected,
            Orientation::Undirected,
            Orientation::Default,
            Orientation::Default,
            Orientation::Default,
            Orientation::Default,
            Orientation::Default,
            Orientation::Default,
            Orientation::Reversed,
        ]);
        let root_mate = EdgeVec::from_iter([
            Orientation::Undirected,
            Orientation::Undirected,
            Orientation::Default,
            Orientation::Default,
            Orientation::Default,
            Orientation::Reversed,
            Orientation::Reversed,
            Orientation::Default,
            Orientation::Reversed,
        ]);
        let outer_target = EdgeVec::from_iter([
            Orientation::Undirected,
            Orientation::Undirected,
            Orientation::Default,
            Orientation::Default,
            Orientation::Default,
            Orientation::Undirected,
            Orientation::Undirected,
            Orientation::Default,
            Orientation::Reversed,
        ]);
        let find_expression = |terms: &[(Atom, OrientationExpression)], target| {
            terms
                .iter()
                .find(|(_, orientation)| orientation.data.orientation == target)
                .map(|(expression, _)| expression.clone())
                .ok_or_else(|| {
                    eyre::eyre!("GL24 selected zero-sample orientation {target:?} is absent")
                })
        };
        let root_raw =
            find_expression(&root_zero, root_target)? + find_expression(&root_zero, root_mate)?;
        let outer_raw = find_expression(&outer_zero, outer_target)?;

        let uv_subgraph =
            InternalSubGraph::cleaned_filter_optimist(uv_contract.clone(), graph.as_ref());
        let crown = graph.dummy_stripped_external_flows_of(&uv_subgraph);
        let child_lmb = graph.try_compatible_sub_lmb(
            &uv_subgraph,
            crown.clone(),
            &graph.loop_momentum_basis,
        )?;
        assert_eq!(
            child_lmb.loop_edges.iter().copied().collect::<Vec<_>>(),
            [EdgeIndex(5)],
        );
        let muv_mass_squared = Atom::var(GS.m_uv_expansion).pow(2);
        let child_denominators = uv_edges.map(|edge| FourDDenominator {
            source_edge: edge,
            momentum: FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(edge))
                .finish(),
            mass_squared: muv_mass_squared.clone(),
            full_expr: Atom::one(),
        });
        let child_cutset = CutSet::empty(graph.n_hedges());
        let (child, _) = graph.cff_from_4d_denominators_in_uv_sub_lmb(
            &child_denominators,
            uv_edges,
            crown.included_iter(),
            &child_lmb,
            ExactUvSubLmbFrame::RetainedPhysicalCrown,
            &child_cutset,
            &options,
            &Atom::one(),
            None,
        )?;
        assert_eq!(child.production_prefactor_factor(), 1);
        let exact_raw = |cff: &CutCFF| -> Result<Atom> {
            Ok(cff
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
                .fold(Atom::Zero, |sum, term| sum + term))
        };
        let child_raw = exact_raw(&child)?;
        let child_raw_for_samples = child_raw.clone();
        let fixed_spatial_quartic = (1..=3)
            .map(|component| GS.emr_mom(quartic_owner, GS.cind(component)).pow(2))
            .fold(Atom::Zero, |sum, component| sum - component)
            .pow(2);
        let root_raw = root_raw * &fixed_spatial_quartic;
        let components_raw = outer_raw * fixed_spatial_quartic * child_raw;
        let root_production = root_raw.clone() * Atom::num(root.production_prefactor_factor());
        let components_production = components_raw.clone()
            * Atom::num(outer.production_prefactor_factor() * child.production_prefactor_factor());

        let zero = F::<ArbPrec>::from_f64(0.0);
        let one = F::<ArbPrec>::from_f64(1.0);
        let muv_mass = F::<ArbPrec>::from(&Rational::from((7, 5)));
        let k0 = F::<ArbPrec>::from(&Rational::from((3, 7)));
        let k1 = F::<ArbPrec>::from(&Rational::from((2, 5)));
        let mut parameters = vec![Atom::var(GS.pi), Atom::var(GS.m_uv_expansion)];
        for edge in 0..graph.underlying.n_edges() {
            let edge = EdgeIndex(edge);
            parameters.extend([GS.ose(edge), cut_energy(edge), GS.emr_mom(edge, GS.cind(0))]);
            for component in 1..=3 {
                parameters.push(GS.emr_mom(edge, GS.cind(component)));
            }
        }
        let values_at = |scale: &F<ArbPrec>| {
            let mut values = vec![
                Complex::new(zero.clone().pi(), zero.clone()),
                Complex::new_re(muv_mass.clone()),
            ];
            for edge in 0..graph.underlying.n_edges() {
                let edge = EdgeIndex(edge);
                let spatial = match edge {
                    EdgeIndex(1) => k0.clone(),
                    EdgeIndex(2) => -&k0,
                    EdgeIndex(3) | EdgeIndex(8) => k1.clone(),
                    EdgeIndex(4) => -&k1,
                    EdgeIndex(5) => scale.clone(),
                    EdgeIndex(6) => &k1 - scale,
                    EdgeIndex(7) => &k0 - &k1,
                    _ => zero.clone(),
                };
                let mass = match edge {
                    EdgeIndex(5) | EdgeIndex(6) => muv_mass.clone(),
                    EdgeIndex(8) => F::<ArbPrec>::from_f64(2.0),
                    _ => one.clone(),
                };
                let energy = if edge == EdgeIndex(0) {
                    F::<ArbPrec>::from_f64(2.0)
                } else {
                    (&mass * &mass + &spatial * &spatial).sqrt()
                };
                values.extend([
                    Complex::new_re(energy.clone()),
                    Complex::new_re(energy.clone()),
                    Complex::new_re(energy),
                ]);
                for component in 1..=3 {
                    values.push(Complex::new_re(if component == 1 {
                        spatial.clone()
                    } else {
                        zero.clone()
                    }));
                }
            }
            values
        };
        let evaluate_arb =
            |expression: Atom, values: &[Complex<F<ArbPrec>>]| -> Result<Complex<F<ArbPrec>>> {
                let rational: ExpressionEvaluator<SymComplex<Fraction<IntegerRing>>> =
                    expression.evaluator(&parameters).build()?;
                let mut arb: ExpressionEvaluator<Complex<F<ArbPrec>>> =
                    rational.map_coeff(&|coefficient| {
                        Complex::new(F::from(&coefficient.re), F::from(&coefficient.im))
                    });
                Ok(arb.evaluate_single(values))
            };
        let scale = F::<ArbPrec>::from_f64(1.0e45);
        let values = values_at(&scale);
        let measure = &scale * &scale * &scale;
        let root_raw = evaluate_arb(root_raw, &values)? * &measure;
        let components_raw = evaluate_arb(components_raw, &values)? * &measure;
        let root_production = evaluate_arb(root_production, &values)? * &measure;
        let components_production = evaluate_arb(components_production, &values)? * &measure;
        let relative_delta = |left: &Complex<F<ArbPrec>>, right: &Complex<F<ArbPrec>>| {
            let distance = (left.clone() - right.clone()).norm().re;
            let left_norm = left.clone().norm().re;
            let right_norm = right.clone().norm().re;
            let norm = if left_norm > right_norm {
                left_norm
            } else {
                right_norm
            };
            if norm.is_zero() {
                distance
            } else {
                distance / norm
            }
        };
        let tolerance = F(ArbPrec::default().epsilon()).sqrt().sqrt().sqrt();
        let raw_delta = relative_delta(&root_raw, &(-&components_raw));
        let production_delta = relative_delta(&root_production, &components_production);
        assert!(
            raw_delta <= tolerance,
            "GL24 selected zero-sample root does not carry the relative parent/outer/child frame: root={root_raw:e}, child*outer={components_raw:e}, relative delta={raw_delta:e}, tolerance={tolerance:e}",
        );
        assert!(
            production_delta <= tolerance,
            "GL24 selected zero-sample production root differs from child times reduced outer: root={root_production:e}, child*outer={components_production:e}, relative delta={production_delta:e}, tolerance={tolerance:e}",
        );
        let spatial_norm = (1..=3)
            .map(|component| GS.emr_mom(quartic_owner, GS.cind(component)).pow(2))
            .fold(Atom::Zero, |sum, component| sum + component);
        let mapped_quartic = |orientation: &OrientationExpression| {
            (orientation.edge_energy_map[usize::from(quartic_owner)]
                .to_atom_gs(&[])
                .pow(2)
                - &spatial_norm)
                .pow(2)
        };
        let sample_expression = |cff: &CutCFF, target: &LinearEnergyExpr| {
            cff.terms
                .iter()
                .filter(|(index, _)| index.lu_cut_order == Some(1))
                .flat_map(|(_, term)| &term.orientations)
                .filter(|orientation| {
                    orientation.orientation.edge_energy_map[usize::from(quartic_owner)]
                        .clone()
                        .canonical()
                        == target.clone().canonical()
                })
                .fold(Atom::Zero, |sum, orientation| {
                    sum + &orientation.expression * mapped_quartic(&orientation.orientation)
                })
        };
        for sample in -2_i64..=2 {
            let target = if sample == 0 {
                LinearEnergyExpr::zero()
            } else {
                LinearEnergyExpr::ose(quartic_owner, sample)
            };
            let root_sample =
                sample_expression(&root, &target) * Atom::num(root.production_prefactor_factor());
            let factorized_sample = sample_expression(&outer, &target)
                * &child_raw_for_samples
                * Atom::num(
                    outer.production_prefactor_factor() * child.production_prefactor_factor(),
                );
            let root_sample = evaluate_arb(root_sample, &values)? * &measure;
            let factorized_sample = evaluate_arb(factorized_sample, &values)? * &measure;
            let delta = relative_delta(&root_sample, &factorized_sample);
            assert!(
                delta <= tolerance,
                "GL24 production root and factorized child times outer differ at q1 sample {sample:+}: root={root_sample:e}, factorized={factorized_sample:e}, relative delta={delta:e}, tolerance={tolerance:e}",
            );
        }
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
    fn exact_nonliteral_physical_mass_alias_keeps_repeated_frame_and_arb_value() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(
            digraph exact_nonliteral_physical_mass_alias {
                node [num=1]
                a -> b [id=0 lmb_id=0 particle=scalar_1]
                a -> b [id=1 particle=scalar_1]
            },
            "scalars"
        )?;
        let edge = EdgeIndex(0);
        let momentum = FunctionBuilder::new(GS.emr_mom)
            .add_arg(usize::from(edge))
            .finish();
        let mass_squared = graph.underlying[edge].particle.mass_atom().pow(2);
        let numerator = GS.emr_mom(edge, GS.cind(0)).pow(2)
            - (1..=3).fold(mass_squared.clone(), |norm_squared, spatial_index| {
                norm_squared + GS.emr_mom(edge, GS.cind(spatial_index)).pow(2)
            });
        let literal = [EdgeIndex(0), EdgeIndex(1)].map(|source_edge| FourDDenominator {
            source_edge,
            momentum: FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(source_edge))
                .finish(),
            mass_squared: graph.underlying[source_edge].particle.mass_atom().pow(2),
            full_expr: Atom::one(),
        });
        let mut nonliteral = literal.clone();
        nonliteral[1].momentum = -momentum;

        let literal_source = GraphThreeDSource::from_exact_denominators(&graph, &literal)?;
        let nonliteral_source = GraphThreeDSource::from_exact_denominators(&graph, &nonliteral)?;
        let literal_parsed = literal_source.to_three_d_parsed_graph()?;
        let nonliteral_parsed = nonliteral_source.to_three_d_parsed_graph()?;
        let literal_groups = repeated_groups(&literal_parsed);
        let nonliteral_groups = repeated_groups(&nonliteral_parsed);
        assert_eq!(literal_groups.len(), 1);
        assert_eq!(nonliteral_groups.len(), 1);
        assert_eq!(literal_groups[0].key, nonliteral_groups[0].key);
        assert_eq!(literal_groups[0].edge_ids.len(), 2);
        assert_eq!(nonliteral_groups[0].edge_ids.len(), 2);

        let options = graph.denominator_only_cff_3d_expression_options();
        let (literal_generated, _, _, _) = graph.generate_3d_expression_for_4d_term(
            &literal_source,
            &options,
            &numerator,
            None,
        )?;
        let (nonliteral_generated, _, _, _) = graph.generate_3d_expression_for_4d_term(
            &nonliteral_source,
            &options,
            &numerator,
            None,
        )?;
        assert_eq!(
            literal_generated.energy_factor_ownership,
            nonliteral_generated.energy_factor_ownership
        );
        assert_eq!(
            literal_generated.core_global_prefactor_sign,
            nonliteral_generated.core_global_prefactor_sign
        );

        let cutset = CutSet::empty(graph.n_hedges());
        let (literal_cff, _) =
            graph.cff_from_4d_denominators(&literal, &cutset, &options, &numerator)?;
        let (nonliteral_cff, _) =
            graph.cff_from_4d_denominators(&nonliteral, &cutset, &options, &numerator)?;
        let exact_sum = |cff: &CutCFF| -> Result<Atom> {
            Ok(cff
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
                * Atom::num(cff.production_prefactor_factor()))
        };
        let fixed_point = |mut expression: Atom| {
            expression = expression
                .replace(graph.underlying[edge].particle.mass_atom())
                .with(Atom::one());
            for (source_edge, spatial) in [
                (EdgeIndex(0), Atom::num(Rational::from((3, 4)))),
                (EdgeIndex(1), Atom::num(Rational::from((-3, 4)))),
            ] {
                expression = expression
                    .replace(GS.emr_mom(source_edge, GS.cind(1)))
                    .with(spatial)
                    .replace(GS.ose(source_edge))
                    .with(Atom::num(Rational::from((5, 4))));
                for spatial_index in 2..=3 {
                    expression = expression
                        .replace(GS.emr_mom(source_edge, GS.cind(spatial_index)))
                        .with(Atom::Zero);
                }
            }
            expression
        };
        let evaluate_arb = |expression: Atom| -> Result<Complex<F<ArbPrec>>> {
            let parameters = [Atom::var(GS.pi)];
            let rational: ExpressionEvaluator<SymComplex<Fraction<IntegerRing>>> = expression
                .evaluator(&parameters)
                .build()
                .map_err(|error| eyre::eyre!("failed to build Q-spelling evaluator: {error}"))?;
            let mut arb: ExpressionEvaluator<Complex<F<ArbPrec>>> =
                rational.map_coeff(&|coefficient| {
                    Complex::new(F::from(&coefficient.re), F::from(&coefficient.im))
                });
            let zero = F(ArbPrec::default());
            Ok(arb.evaluate_single(&[Complex::new(zero.clone().pi(), zero)]))
        };
        let literal_value = evaluate_arb(fixed_point(exact_sum(&literal_cff)?))?;
        let nonliteral_value = evaluate_arb(fixed_point(exact_sum(&nonliteral_cff)?))?;
        assert!(!literal_value.clone().norm().re.is_zero());
        let distance = (literal_value.clone() - nonliteral_value.clone()).norm().re;
        let literal_norm = literal_value.clone().norm().re;
        let nonliteral_norm = nonliteral_value.clone().norm().re;
        let scale = if literal_norm > nonliteral_norm {
            literal_norm
        } else {
            nonliteral_norm
        };
        let relative_distance = distance / scale;
        let tolerance = F(ArbPrec::default().epsilon()).sqrt().sqrt().sqrt();
        assert!(
            relative_distance <= tolerance,
            "Q(1) and algebraically equivalent -Q(0) sources differ: literal={literal_value:e}, nonliteral={nonliteral_value:e}, relative delta={relative_distance:e}, tolerance={tolerance:e}"
        );
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
    fn exact_cff_powered_rational_identities_match_at_arb_precision() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph exact_powered_rational_identities {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]

            incoming -> a [id=0]
            a -> x [id=1 lmb_id=0]
            x -> y [id=2]
            b -> y [id=3]
            a -> b [id=4]
            b -> outgoing [id=5]
        })?;
        let repeated_edges = [EdgeIndex(1), EdgeIndex(2), EdgeIndex(3)];
        let spectator_edge = EdgeIndex(4);
        let carrier = repeated_edges[0];
        let carrier_momentum = FunctionBuilder::new(GS.emr_mom)
            .add_arg(usize::from(carrier))
            .finish();
        let mass_squared = graph.underlying[carrier].particle.mass_atom().pow(2);
        let full_expr = GS.emr_mom(carrier, GS.cind(0)).pow(2)
            - (1..=3).fold(mass_squared.clone(), |norm_squared, spatial_index| {
                norm_squared + GS.emr_mom(carrier, GS.cind(spatial_index)).pow(2)
            });
        let repeated_denominators = repeated_edges
            .into_iter()
            .zip([1, 1, -1])
            .map(|(source_edge, routing_sign)| FourDDenominator {
                source_edge,
                momentum: Atom::num(routing_sign) * &carrier_momentum,
                mass_squared: mass_squared.clone(),
                full_expr: full_expr.clone(),
            })
            .collect::<Vec<_>>();
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
        let powered_denominators = std::iter::once(spectator.clone())
            .chain(repeated_denominators.iter().cloned())
            .collect::<Vec<_>>();
        let lower_denominators = [
            std::iter::once(spectator.clone())
                .chain([
                    repeated_denominators[0].clone(),
                    repeated_denominators[2].clone(),
                ])
                .collect::<Vec<_>>(),
            std::iter::once(spectator)
                .chain([repeated_denominators[2].clone()])
                .collect::<Vec<_>>(),
        ];
        assert_eq!(
            powered_denominators
                .iter()
                .skip(1)
                .map(|denominator| denominator.source_edge)
                .collect::<BTreeSet<_>>()
                .len(),
            3,
            "the cubic denominator must retain three distinct provenance owners"
        );
        assert_eq!(
            repeated_denominators[2].momentum,
            -carrier_momentum.clone(),
            "one occurrence must exercise the even D(-Q) denominator routing"
        );

        let options = graph.denominator_only_cff_3d_expression_options();
        let mut production_options = options.clone();
        production_options.energy_degree_bounds = None;
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let production_contract_edges =
            graph.paired_edges(&graph.tree_edges.subtract(&graph.initial_state_cut));
        let production = graph.generate_3d_expression_for_integrand(
            &production_contract_edges,
            &canonization,
            &production_options,
            None,
        )?;
        let lu_surface = graph
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
            .expect("the dotted bubble contains an external-shift LU surface");
        assert!(lu_surface.esurface_ids.iter().all(|esurface_id| {
            let alternative = &production.expression.surfaces.esurface_cache[*esurface_id].energies;
            alternative.contains(&spectator_edge)
                && alternative.iter().any(|edge| repeated_edges.contains(edge))
        }));

        let uncut = CutSet::empty(graph.n_hedges());
        let exact_sum = |graph: &mut Graph,
                         denominators: &[FourDDenominator],
                         numerator: &Atom|
         -> Result<Atom> {
            let (cff, _) =
                graph.cff_from_4d_denominators(denominators, &uncut, &options, numerator)?;
            assert_eq!(cff.terms.len(), 1, "an uncut exact CFF has one sector");
            let prefactor = Atom::num(cff.production_prefactor_factor());
            let term = cff
                .terms
                .values()
                .next()
                .expect("the uncut exact-CFF sector exists");
            let sum = term
                .orientations
                .iter()
                .try_fold(Atom::Zero, |sum, orientation| {
                    Ok::<_, color_eyre::Report>(
                        sum + orientation.expression.clone()
                            * term.map_exact_source_numerator(&orientation.orientation)?,
                    )
                })?;
            Ok(sum.replace(GS.den(W_.a_, W_.b_, W_.c_, W_.d_)).with(W_.d_) * prefactor)
        };
        let graph_mass_squared = (0..graph.underlying.n_edges())
            .map(|edge| {
                graph.underlying[EdgeIndex(edge)]
                    .particle
                    .mass_atom()
                    .pow(2)
            })
            .collect::<Vec<_>>();
        let external_edges = graph
            .external_momentum_edge_order()
            .into_iter()
            .collect::<BTreeSet<_>>();
        let representative_esurface =
            &production.expression.surfaces.esurface_cache[lu_surface.esurface_ids[0]];
        let rescale_star = Atom::num(Rational::from((3, 4)));
        let external_spatial_momentum = Atom::num(Rational::from((25, 12)));
        let scaled_spatial_momentum = |edge: EdgeIndex, rescale: Atom| {
            if edge == spectator_edge {
                &external_spatial_momentum - rescale
            } else {
                rescale
            }
        };
        let scaled_on_shell_energy = |edge: EdgeIndex, rescale: Atom| {
            (graph_mass_squared[usize::from(edge)].clone()
                + scaled_spatial_momentum(edge, rescale).pow(2))
            .sqrt()
        };
        let root_energy_sum = representative_esurface
            .energies
            .iter()
            .map(|edge| scaled_on_shell_energy(*edge, rescale_star.clone()))
            .fold(Atom::Zero, |sum, energy| sum + energy);
        let external_shift_coefficient = representative_esurface
            .external_shift
            .iter()
            .map(|(_, coefficient)| *coefficient)
            .sum::<i64>();
        assert_ne!(
            external_shift_coefficient, 0,
            "the selected LU surface must have a nonzero net external shift"
        );
        let external_energy = -root_energy_sum / Atom::num(external_shift_coefficient);
        let rescale_expression = |mut expression: Atom| {
            let rescale = Atom::var(GS.rescale);
            for edge in 0..graph_mass_squared.len() {
                let edge = EdgeIndex(edge);
                if external_edges.contains(&edge) {
                    expression = expression
                        .replace(GS.emr_mom(edge, GS.cind(0)))
                        .with(external_energy.clone())
                        .replace(GS.emr_mom(edge, GS.cind(1)))
                        .with(external_spatial_momentum.clone());
                } else {
                    expression = expression
                        .replace(GS.emr_mom(edge, GS.cind(1)))
                        .with(scaled_spatial_momentum(edge, rescale.clone()))
                        .replace(GS.ose(edge))
                        .with(scaled_on_shell_energy(edge, rescale.clone()))
                        .replace(cut_energy(edge))
                        .with(scaled_on_shell_energy(edge, rescale.clone()));
                }
                for spatial_index in 2..=3 {
                    expression = expression
                        .replace(GS.emr_mom(edge, GS.cind(spatial_index)))
                        .with(Atom::Zero);
                }
            }
            expression
        };
        let eta = rescale_expression(representative_esurface.to_atom(&[]));
        assert!(
            eta.replace(GS.rescale)
                .with(rescale_star.clone())
                .expand()
                .is_zero(),
            "the exact test point must lie on the selected LU surface"
        );
        let raised_laurent_residue = |expression: Atom| -> Result<Atom> {
            // Expanding the complete mapped-numerator times CFF rational
            // function supplies the same t derivatives as raised-cut pass two,
            // without assigning representation-dependent raw channels.
            Ok(rescale_expression(expression)
                .series(GS.rescale, rescale_star.clone(), 0)
                .map_err(|error| eyre::eyre!("failed to expand selected LU residue: {error}"))?
                .coefficient(Rational::from(-1)))
        };
        let evaluate_arb = |mut expression: Atom| -> Result<Complex<F<ArbPrec>>> {
            for (edge, mass_squared) in graph_mass_squared.iter().enumerate() {
                let edge = EdgeIndex(edge);
                let on_shell_energy = (1..=3)
                    .fold(mass_squared.clone(), |norm_squared, spatial_index| {
                        norm_squared + GS.emr_mom(edge, GS.cind(spatial_index)).pow(2)
                    })
                    .sqrt();
                expression = expression.replace(GS.ose(edge)).with(on_shell_energy);
            }
            let mut parameters = vec![Atom::var(GS.pi)];
            parameters.extend(
                (0..graph_mass_squared.len())
                    .flat_map(|edge| (0..=3).map(move |component| (edge, component)))
                    .map(|(edge, component)| GS.emr_mom(EdgeIndex(edge), GS.cind(component))),
            );
            let rational: ExpressionEvaluator<SymComplex<Fraction<IntegerRing>>> = expression
                .evaluator(&parameters)
                .build()
                .map_err(|error| eyre::eyre!("failed to build exact-CFF evaluator: {error}"))?;
            let mut arb: ExpressionEvaluator<Complex<F<ArbPrec>>> =
                rational.map_coeff(&|coefficient| {
                    Complex::new(F::from(&coefficient.re), F::from(&coefficient.im))
                });
            let zero = F(ArbPrec::default());
            let values = parameters
                .iter()
                .enumerate()
                .map(|(index, _)| {
                    if index == 0 {
                        Complex::new(zero.clone().pi(), zero.clone())
                    } else {
                        Complex::new(
                            F::<ArbPrec>::from_f64(0.17 + 0.037 * index as f64),
                            zero.clone(),
                        )
                    }
                })
                .collect::<Vec<_>>();
            Ok(arb.evaluate_single(&values))
        };
        // One eighth of ArbPrec's requested precision allows substantial and
        // construction-dependent bit loss while still scaling with precision.
        let tolerance = F(ArbPrec::default().epsilon()).sqrt().sqrt().sqrt();
        let mut failures = Vec::new();
        for (power, lower) in [1_i64, 2].into_iter().zip(&lower_denominators) {
            let denominator_numerator = GS
                .den(
                    usize::from(carrier),
                    &carrier_momentum,
                    &mass_squared,
                    &full_expr,
                )
                .pow(power);
            let powered_uncut =
                exact_sum(&mut graph, &powered_denominators, &denominator_numerator)?;
            let lower_uncut = exact_sum(&mut graph, lower, &Atom::one())?;
            let comparisons = [
                (
                    "uncut rational function",
                    powered_uncut.clone(),
                    lower_uncut.clone(),
                ),
                (
                    "selected LU Laurent residue",
                    raised_laurent_residue(powered_uncut)?,
                    raised_laurent_residue(lower_uncut)?,
                ),
            ];
            for (comparison, powered_expression, lower_expression) in comparisons {
                let powered_value = evaluate_arb(powered_expression)?;
                let lower_value = evaluate_arb(lower_expression)?;
                assert!(
                    !(powered_value.re.is_nan()
                        || powered_value.im.is_nan()
                        || powered_value.re.is_infinite()
                        || powered_value.im.is_infinite()
                        || lower_value.re.is_nan()
                        || lower_value.im.is_nan()
                        || lower_value.re.is_infinite()
                        || lower_value.im.is_infinite()),
                    "{comparison} must evaluate to finite Arb values: powered={powered_value:e}, lower={lower_value:e}"
                );
                assert!(
                    !lower_value.clone().norm().re.is_zero(),
                    "{comparison} must provide a nonzero identity oracle"
                );
                let distance = (powered_value.clone() - lower_value.clone()).norm().re;
                let powered_norm = powered_value.norm().re;
                let lower_norm = lower_value.norm().re;
                let scale = if powered_norm > lower_norm {
                    powered_norm
                } else {
                    lower_norm
                };
                let relative_distance = if scale.is_zero() {
                    distance.clone()
                } else {
                    distance.clone() / scale
                };
                if relative_distance > tolerance {
                    failures.push(format!(
                        "D(Q)^{power}/D(Q)^3 differs from 1/D(Q)^{} for {comparison}: powered={powered_value:e}, lower={lower_value:e}, relative delta={relative_distance:e}, tolerance={tolerance:e}",
                        3 - power
                    ));
                }
            }
        }
        assert!(
            failures.is_empty(),
            "exact-CFF powered rational identity failures:\n{}",
            failures.join("\n")
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
