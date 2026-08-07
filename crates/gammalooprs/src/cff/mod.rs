use std::{collections::BTreeMap, fmt::Display};

use bincode_trait_derive::{Decode, Encode};
use linnet::half_edge::subgraph::{Inclusion, SubGraphLike, SubSetLike, SubSetOps};
use serde::{Deserialize, Serialize};
use symbolica::atom::{Atom, AtomCore};

use crate::{
    cff::{
        expression::{
            GammaLoopOrientationExpression, OrientationExpression, OrientationID,
            OrientationSelector, ThreeDExpression,
        },
        orientations::GraphOrientation,
        surface::GammaLoopSurfaceCache,
    },
    graph::{
        FeynmanGraph, FourDDenominator, Graph, GraphThreeDSource, cuts::CutSet,
        get_cff_inverse_energy_product_impl,
    },
    settings::global::OrientationPattern,
    utils::GS,
    uv::Integrands,
};
use color_eyre::Result;
use three_dimensional_reps::{Generate3DExpressionOptions, ThreeDGraphSource};

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
    // One denominator expression per reduced-graph energy map. Keep the map
    // until the UV localizer can assign its full-graph OrientationID.
    pub(crate) orientations: Vec<CFFOrientationTerm>,
}

impl CFFTerm {
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
) -> Vec<(CutCFFIndex, ThreeDExpression<OrientationID>)> {
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
            |expression| expression.select_esurface_residue(left_threshold),
        );
    }

    if let Some(lu_cut) = cutset.residue_selector.lu_cut.as_ref() {
        residues =
            apply_indexed_residue_selection(residues, CutCffResidueAxis::LuCut, |expression| {
                // Threshold and LU residues stay in the canonical E-surface family. The removed
                // generated-basis branches belonged to the distinct confluent/LTD
                // representation: re-reading those generated repeated-channel coordinates would
                // have put Laurent coefficients back on unresolved threshold surfaces and
                // changed the selected-denominator sign convention. They also depended on
                // lower-sector edge-set alternatives that main's selector does not expose.
                // Residues here therefore use GammaLoop's positive-energy Cutkosky convention.
                expression.select_esurface_residue(lu_cut)
            });
    }

    residues
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
            use_generated_half_edges,
            physical_surfaces,
            exact_ose_replacements,
            inverse_energy_product,
            physical_edge_map,
            cff_loop_number,
            contract_subgraph,
        ) = {
            let source = GraphThreeDSource::from_exact_denominators(self, denominators)?;
            let generated =
                self.generate_3d_expression_for_4d_term(&source, options, analysis_numerator)?;
            let use_generated_half_edges =
                generation::generated_cff_expression_uses_variant_half_edges(&generated);
            let physical_surfaces = generated
                .surfaces
                .linear_surface_cache
                .iter()
                .map(|surface| source.physical_linear_surface(surface))
                .collect::<Vec<_>>();
            (
                generated,
                use_generated_half_edges,
                physical_surfaces,
                source
                    .exact_ose_replacements()
                    .expect("exact 4D source has occurrence-local energy replacements"),
                source
                    .exact_inverse_energy_product()
                    .expect("exact 4D source has an occurrence-local energy product"),
                source
                    .physical_energy_edge_index_map()
                    .expect("exact 4D source has a physical edge map"),
                source.active_loop_count(),
                source.contract_subgraph(),
            )
        };
        let cff = self.convert_4d_expression_surfaces(
            generated,
            use_generated_half_edges,
            &physical_surfaces,
        )?;
        let residues = select_indexed_cff_residues(cff, cutset);
        let cff_phase = (-Atom::i()).pow(cff_loop_number as i64);
        let cff_normalization = cff_phase / (Atom::var(GS.pi) * 2).pow(3 * cff_loop_number as i64);

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
            };
            for orientation in expr.orientations {
                let expression = orientation
                    .to_atom_gs()
                    .replace_multiple(&replacement_rules)
                    .replace_multiple(&exact_ose_replacements)
                    * if use_generated_half_edges {
                        Atom::one()
                    } else {
                        inverse_energy_product.clone()
                    }
                    * &cff_normalization;
                let mut physical_energy_maps = BTreeMap::new();
                for (local_edge, physical_edge) in &physical_edge_map.internal {
                    let physical_edge_id = linnet::half_edge::involution::EdgeIndex(*physical_edge);
                    if contract_subgraph.includes(&self[&physical_edge_id].1) {
                        continue;
                    }
                    let Some(energy_map) = orientation.edge_energy_map.get(*local_edge) else {
                        continue;
                    };
                    let energy_map = energy_map.clone().remap_energy_edges(
                        &physical_edge_map.internal,
                        &physical_edge_map.external,
                    );
                    if let Some(existing) =
                        physical_energy_maps.insert(*physical_edge, energy_map.clone())
                        && existing != energy_map
                    {
                        return Err(eyre::eyre!(
                            "exact 4D denominator occurrences for physical edge {physical_edge} induce incompatible affine energy maps"
                        ));
                    }
                }
                // Distinct exact factors remain in `expression`; this remap is
                // only the compatible physical representative used to attach
                // their term to a production numerator energy map.
                let mut physical_orientation = orientation;
                physical_orientation.remap_energy_edge_indices(&physical_edge_map);
                cff_term.orientations.push(CFFOrientationTerm {
                    expression,
                    orientation: physical_orientation,
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
        // production graph. The old repeated-channel branch switched only this side to a
        // generated residue coordinate, so
        // its surviving edge-energy maps could not be exact restrictions of production maps.
        // Repeated thresholds are selected below in the canonical E-surface family; a distinct
        // confluent/LTD representation remains deferred rather than being selected implicitly.
        let cff = self.generate_3d_expression_for_integrand(
            &contract_edges,
            &canonize_esurface,
            options,
            analysis_numerator,
            false,
        )?;

        let residues = select_indexed_cff_residues(cff, cutset);

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
        crate::debug_tags!(#cff, #trace;
            stage = "graph_cff_normalization",
            graph = %self.name,
            cff_loop_number = cff_loop_number,
            log.cff_normalization = cff_normalization,
            "Graph CFF normalization"
        );

        let mut terms = BTreeMap::new();

        for (cut_cff_index, expr) in residues {
            let use_generated_half_edges = expr.orientations.iter().any(|orientation| {
                orientation
                    .variants
                    .iter()
                    .any(|variant| !variant.half_edges.is_empty())
            });
            let replacement_rules = if cutset.canonicalize_external_shifts {
                expr.surfaces
                    .get_all_replacements_gs_in_lmb(&[], &self.loop_momentum_basis)
            } else {
                expr.surfaces.get_all_replacements_gs(&[])
            };
            let mut cff_term = CFFTerm {
                orientations: vec![],
            };
            for orientation in expr.orientations.iter().filter(|orientation| {
                orientation_pattern.filter_orientation(&orientation.data.orientation)
            }) {
                let eta_expr = orientation.to_atom_gs();
                let mut ose_expr = eta_expr.replace_multiple(&replacement_rules);
                if !use_generated_half_edges {
                    ose_expr *= get_cff_inverse_energy_product_impl(
                        self,
                        &graph_without_is_cut,
                        &contract_edges,
                    );
                }

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
        assert!(cff.terms.values().any(|term| !term.orientations.is_empty()));
        let on_shell_energy = (1..=3)
            .fold(mass_squared, |norm_squared, spatial_index| {
                norm_squared + GS.emr_mom(edge, GS.cind(spatial_index)).pow(2)
            })
            .sqrt();
        // The confluent residue owns its three half-edge factors, so the exact
        // bridge must not append the pure-CFF global energy product.
        let expected_expression =
            -Atom::i() / (Atom::num(32) * Atom::var(GS.pi).pow(3) * on_shell_energy.pow(3));
        for term in cff.terms.values().flat_map(|term| &term.orientations) {
            let orientation = &term.orientation.data.orientation;
            assert_ne!(orientation[edge], Orientation::Undirected);
            assert_eq!(orientation[EdgeIndex::from(1)], Orientation::Undirected);
            assert_eq!(
                orientation[EdgeIndex::from(2)],
                orientation[EdgeIndex::from(3)]
            );
            assert_eq!(term.expression, expected_expression);
        }
        Ok(())
    }

    #[test]
    fn exact_cff_bridges_opposite_repeated_channel_routing() -> Result<()> {
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
        let on_shell_energy = (1..=3)
            .fold(mass_squared, |norm_squared, spatial_index| {
                norm_squared + GS.emr_mom(EdgeIndex(0), GS.cind(spatial_index)).pow(2)
            })
            .sqrt();
        let expected_expression =
            -Atom::i() / (Atom::num(32) * Atom::var(GS.pi).pow(3) * on_shell_energy.pow(3));
        let terms = cff
            .terms
            .values()
            .flat_map(|term| &term.orientations)
            .collect::<Vec<_>>();

        assert!(!terms.is_empty());
        assert!(
            terms
                .into_iter()
                .all(|term| term.expression == expected_expression)
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
    fn exact_original_source_preserves_production_affine_maps() -> Result<()> {
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
        let ordinary = ordinary
            .terms
            .values()
            .flat_map(|term| &term.orientations)
            .collect::<Vec<_>>();

        for exact in exact.terms.values().flat_map(|term| &term.orientations) {
            assert!(
                ordinary.iter().any(|ordinary| {
                    ordinary.orientation.data.orientation == exact.orientation.data.orientation
                        && ordinary.orientation.edge_energy_map == exact.orientation.edge_energy_map
                }),
                "exact map {:?} with energies {:?} is absent from ordinary maps {:?}",
                exact.orientation.data.orientation,
                exact.orientation.edge_energy_map,
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
