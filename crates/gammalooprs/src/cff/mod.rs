use std::{collections::BTreeMap, fmt::Display};

use bincode_trait_derive::{Decode, Encode};
use linnet::half_edge::subgraph::{SubGraphLike, SubSetLike, SubSetOps};
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
    settings::global::OrientationPattern,
    utils::GS,
    uv::Integrands,
};
use color_eyre::Result;
use three_dimensional_reps::{CffEnergyFactorOwnership, Generate3DExpressionOptions};

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
}

pub struct CFFTerm {
    // Ordinary CFF maps retain production identity for direct-3D localization;
    // exact 4D maps instead remain source-local and share their parent mapper.
    pub(crate) orientations: Vec<CFFOrientationTerm>,
    exact_source_energy_mapper: Option<ExactSourceEnergyMapper>,
}

impl CFFTerm {
    pub(crate) fn map_exact_source_numerator(
        &self,
        orientation: &OrientationExpression,
        numerator: &Atom,
    ) -> Result<Atom> {
        self.exact_source_energy_mapper
            .as_ref()
            .ok_or_else(|| eyre::eyre!("ordinary CFF term has no exact-source energy mapper"))?
            .map_numerator(
                &orientation.loop_energy_map,
                &orientation.edge_energy_map,
                numerator,
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

pub struct CutCFF {
    pub terms: BTreeMap<CutCFFIndex, CFFTerm>,
}

impl CutCFF {
    pub fn expression_with_selectors(&self) -> Integrands {
        self.terms
            .iter()
            .map(|(index, term)| (*index, term.expression_with_selectors()))
            .collect()
    }
}

#[derive(Clone, Copy)]
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
    pub(crate) fn cff_from_4d_denominators(
        &mut self,
        denominators: &[FourDDenominator],
        cutset: &CutSet,
        options: &Generate3DExpressionOptions,
        analysis_numerator: &Atom,
    ) -> Result<(CutCFF, linnet::half_edge::subgraph::SuBitGraph)> {
        let (
            generated,
            physical_surfaces,
            physical_energy_edges,
            exact_source_energy_mapper,
            inverse_energy_product,
            cff_loop_number,
            contract_subgraph,
        ) = {
            let source = GraphThreeDSource::from_exact_denominators(self, denominators)?;
            let generated =
                self.generate_3d_expression_for_4d_term(&source, options, analysis_numerator)?;
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
                source
                    .exact_source_energy_mapper()
                    .expect("exact 4D source has an owned parent-energy mapper"),
                source
                    .exact_inverse_energy_product()
                    .expect("exact 4D source has an occurrence-local energy product"),
                source.active_loop_count(),
                source.contract_subgraph(),
            )
        };
        let generated = self.convert_4d_expression_surfaces(generated, &physical_surfaces)?;
        let energy_factor_ownership = generated.energy_factor_ownership;
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
        let residues = select_indexed_cff_residues(cff, cutset)?;
        let cff_phase = (-Atom::i()).pow(cff_loop_number as i64);
        let cff_normalization = cff_phase / (Atom::var(GS.pi) * 2).pow(3 * cff_loop_number as i64);
        let cff_energy_factor = match energy_factor_ownership {
            CffEnergyFactorOwnership::GlobalSourceProduct => inverse_energy_product,
            CffEnergyFactorOwnership::VariantLocal => Atom::one(),
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
                exact_source_energy_mapper: Some(exact_source_energy_mapper.clone()),
            };
            for orientation in expr.orientations {
                let expression = orientation
                    .to_atom_gs()
                    .replace_multiple(&replacement_rules)
                    .replace_multiple(exact_source_energy_mapper.exact_ose_replacements())
                    * &cff_energy_factor
                    * &cff_normalization;
                // Distinct exact factors remain in `expression`. The
                // orientation is deliberately source-local: its affine map
                // evaluates the complete parent numerator directly, rather
                // than being remapped to a production OrientationID.
                cff_term.orientations.push(CFFOrientationTerm {
                    expression,
                    orientation,
                });
            }
            terms.insert(cut_cff_index, cff_term);
        }
        Ok((CutCFF { terms }, contract_subgraph))
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

        for (cut_cff_index, expr) in residues {
            let replacement_rules = if cutset.canonicalize_external_shifts {
                expr.surfaces
                    .get_all_replacements_gs_in_lmb(&[], &self.loop_momentum_basis)
            } else {
                expr.surfaces.get_all_replacements_gs(&[])
            };
            let mut cff_term = CFFTerm {
                orientations: vec![],
                exact_source_energy_mapper: None,
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
                });
            }
            terms.insert(cut_cff_index, cff_term);
        }

        let cut_cff = CutCFF { terms };
        Ok(cut_cff)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{dot, graph::parse::IntoGraph, initialisation::test_initialise};
    use linnet::half_edge::involution::{EdgeIndex, Orientation};
    use linnet::half_edge::subgraph::SuBitGraph;
    use symbolica::atom::FunctionBuilder;

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
    fn exact_cff_handles_opposite_repeated_routing_without_a_sign_bridge() -> Result<()> {
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
            assert_eq!(
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
    fn exact_original_source_matches_production_affine_maps_after_projection() -> Result<()> {
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
        let physical_energy_edges =
            GraphThreeDSource::from_exact_denominators(&graph, &denominators)?
                .physical_energy_edge_index_map()
                .expect("exact source has a physical energy-edge projection");
        let (exact, _) =
            graph.cff_from_4d_denominators(&denominators, &cutset, &options, &Atom::one())?;
        let ordinary = ordinary
            .terms
            .values()
            .flat_map(|term| &term.orientations)
            .collect::<Vec<_>>();

        for exact in exact.terms.values().flat_map(|term| &term.orientations) {
            let mut exact_orientation = exact.orientation.clone();
            exact_orientation.remap_energy_edge_indices(&physical_energy_edges);
            assert!(
                ordinary.iter().any(|ordinary| {
                    ordinary.orientation.data.orientation == exact_orientation.data.orientation
                        && ordinary.orientation.edge_energy_map == exact_orientation.edge_energy_map
                }),
                "exact map {:?} with energies {:?} is absent from ordinary maps {:?}",
                exact_orientation.data.orientation,
                exact_orientation.edge_energy_map,
                ordinary
                    .iter()
                    .map(|term| (
                        &term.orientation.data.orientation,
                        &term.orientation.edge_energy_map,
                    ))
                    .collect::<Vec<_>>()
            );
        }
        Ok(())
    }
}
