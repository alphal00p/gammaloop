use std::{collections::BTreeMap, fmt::Display};

use bincode_trait_derive::{Decode, Encode};
use linnet::half_edge::{
    involution::{EdgeVec, Orientation},
    subgraph::{SubGraphLike, SubSetLike, SubSetOps},
};
use symbolica::atom::{Atom, AtomCore};

use self::generation::{GammaLoopGraphCffExt, GammaLoopSurfaceCacheExt};
use crate::{
    cff::orientations::GraphOrientation,
    graph::{FeynmanGraph, Graph, cuts::CutSet, get_cff_inverse_energy_product_impl},
    settings::global::OrientationPattern,
    utils::GS,
    uv::Integrands,
};
use color_eyre::Result;

pub mod orientations;
//pub mod cut_expression;
pub mod esurface;
pub mod expression;
mod feynkit;
pub mod generation;
pub mod hsurface;
pub mod surface;

pub struct CFFTerm {
    // One per orientation
    pub expression: Vec<Atom>,
    pub orientations: Vec<EdgeVec<Orientation>>,
}

impl CFFTerm {
    pub fn expression_with_selectors(&self) -> Atom {
        let mut result = Atom::Zero;
        for (expr, orient) in self.expression.iter().zip(self.orientations.iter()) {
            result += expr.clone() * orient.orientation_thetas();
        }
        result
    }
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Hash, Encode, Decode)]
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

impl Graph {
    pub fn cff<S: SubGraphLike + SubSetLike>(
        &mut self,
        contract_subgraph: &S,
        cutset: &CutSet,
        orientation_pattern: &OrientationPattern,
    ) -> Result<CutCFF> {
        let canonize_esurface = self.get_esurface_canonization(&self.loop_momentum_basis);
        let mut contract_edges = vec![];

        for (p, eid, _) in self.iter_edges_of(contract_subgraph) {
            if p.is_paired() {
                contract_edges.push(eid);
            }
        }

        let cff = [(
            CutCFFIndex::new_all_none(),
            self.generate_cff(&contract_edges, &canonize_esurface, orientation_pattern)?,
        )];

        let mut residues = BTreeMap::new();

        cff.into_iter()
            .flat_map(|(index, cff_expression)| {
                if let Some(right_threshold) = cutset.residue_selector.right_th_cut.as_ref() {
                    cff_expression
                        .select_energy_surface_residue(right_threshold)
                        .into_iter()
                        .enumerate()
                        .map(|(i, residue)| {
                            let mut new_index = index;
                            new_index.right_threshold_order = Some(i + 1);
                            (new_index, residue)
                        })
                        .collect()
                } else {
                    vec![(index, cff_expression)]
                }
            })
            .flat_map(|(index, cff_expression)| {
                if let Some(left_threshold) = cutset.residue_selector.left_th_cut.as_ref() {
                    cff_expression
                        .select_energy_surface_residue(left_threshold)
                        .into_iter()
                        .enumerate()
                        .map(|(i, residue)| {
                            let mut new_index = index;
                            new_index.left_threshold_order = Some(i + 1);
                            (new_index, residue)
                        })
                        .collect()
                } else {
                    vec![(index, cff_expression)]
                }
            })
            .flat_map(|(index, cff_expression)| {
                if let Some(lu_cut) = cutset.residue_selector.lu_cut.as_ref() {
                    cff_expression
                        .select_energy_surface_residue(lu_cut)
                        .into_iter()
                        .enumerate()
                        .map(|(i, residue)| {
                            let mut new_index = index;
                            new_index.lu_cut_order = Some(i + 1);
                            (new_index, residue)
                        })
                        .collect()
                } else {
                    vec![(index, cff_expression)]
                }
            })
            .for_each(|(index, residue)| {
                residues.insert(index, residue);
            });

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

        let replacement_rules = if cutset.canonicalize_external_shifts {
            self.surface_cache
                .get_all_replacements_in_lmb(&[], &self.loop_momentum_basis)
        } else {
            self.surface_cache.get_all_replacements(&[])
        };

        for (cut_cff_index, expr) in residues.into_iter() {
            let mut cff_term = CFFTerm {
                expression: vec![],
                orientations: vec![],
            };
            for orientation in expr.orientations.iter() {
                let eta_expr = orientation.to_atom();
                let mut ose_expr = eta_expr.replace_multiple(&replacement_rules);

                let inverse_energies = get_cff_inverse_energy_product_impl(
                    self,
                    &graph_without_is_cut,
                    &contract_edges,
                );

                ose_expr *= inverse_energies;
                ose_expr *= cff_normalization.clone();

                crate::debug_tags!(#cff, #trace;
                    stage = "graph_cff_term_expr",
                    graph = %self.name,
                    cut_index = ?cut_cff_index,
                    log.expr = ose_expr,
                    "Graph CFF term expression"
                );
                // println!("ose expr :{}", ose_expr);
                cff_term.expression.push(ose_expr);
                cff_term
                    .orientations
                    .push(orientation.data.orientation.clone());
            }
            terms.insert(cut_cff_index, cff_term);
        }

        let cut_cff = CutCFF { terms };
        Ok(cut_cff)
    }
}

#[cfg(test)]
mod tests {
    use feynkit_cff::RaisedEnergySurfaceGroup;
    use linnet::half_edge::subgraph::SuBitGraph;

    use super::*;
    use crate::{
        finalized_runtime_dot, graph::parse::IntoFinalizedRuntimeGraph,
        initialisation::test_initialise,
    };

    #[test]
    fn residue_terms_preserve_attached_tadpole_weight_and_normalization() {
        test_initialise().unwrap();
        assert!(symbolica::license::LicenseManager::is_licensed());
        let mut graph: Graph = finalized_runtime_dot!(
            digraph bubble_with_tadpole {
                projector=1
                edge [particle=scalar_1 num=1]
                node [num=1]
                ext [style=invis]
                ext -> A:0 [id=3 sink="{ufo_order:0}"]
                B:1 -> ext [id=2 source="{ufo_order:2}"]
                A:2 -> B:3 [id=1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:0}"]
                A:4 -> B:5 [id=0 source="{ufo_order:2}" sink="{ufo_order:1}"]
                A:6 -> A:7 [id=4 lmb_id=1 source="{ufo_order:3}" sink="{ufo_order:4}"]
            },
            "scalars"
        )
        .unwrap();
        assert_eq!(graph.get_loop_number(), 2);
        let pattern = OrientationPattern::default();
        let rewrite = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let generated = graph.generate_cff(&[], &rewrite, &pattern).unwrap();
        let surface = *generated.expression.energy_surface_ids().first().unwrap();
        let group = RaisedEnergySurfaceGroup {
            surface_ids: vec![surface],
            max_occurrence: 1,
        };
        let residue = generated
            .select_energy_surface_residue(&group)
            .into_iter()
            .next()
            .unwrap();

        let mut cuts = CutSet::empty(graph.n_hedges());
        cuts.residue_selector.lu_cut = Some(group);
        let contract: SuBitGraph = graph.empty_subgraph();
        let actual = graph.cff(&contract, &cuts, &pattern).unwrap();
        let replacements = graph.surface_cache.get_all_replacements(&[]);
        let inverse_energies =
            get_cff_inverse_energy_product_impl(&graph, &graph.underlying.full_filter(), &[]);
        // Keep the original two-loop phase (-i)^2 and spatial measure.
        let normalization = -Atom::num(1) / (Atom::var(GS.pi) * 2).pow(6);
        let expected = residue
            .orientations
            .iter()
            .map(|orientation| {
                // The attached free tadpole contributes exactly one half.
                let weighted_tree = Atom::num((1, 2)) * orientation.expression.to_atom_inverse();
                weighted_tree.replace_multiple(&replacements)
                    * inverse_energies.clone()
                    * normalization.clone()
            })
            .collect::<Vec<_>>();
        assert!(expected.iter().any(|atom| !atom.is_zero()));
        assert_eq!(actual.terms.len(), 1);
        let term = &actual.terms[&CutCFFIndex {
            left_threshold_order: None,
            right_threshold_order: None,
            lu_cut_order: Some(1),
        }];
        assert_eq!(term.expression, expected);
        assert_eq!(
            term.orientations,
            residue
                .orientations
                .iter()
                .map(|orientation| orientation.data.orientation.clone())
                .collect::<Vec<_>>()
        );
    }
}
