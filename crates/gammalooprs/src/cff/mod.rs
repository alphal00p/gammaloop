use linnet::half_edge::{
    involution::{EdgeVec, Orientation},
    subgraph::SubSetLike,
};
use symbolica::atom::{Atom, AtomCore};

use crate::{
    cff::expression::GraphOrientation,
    graph::{FeynmanGraph, Graph, cuts::CutSet},
};
use color_eyre::Result;
use eyre::eyre;

pub mod cff_graph;
pub mod cut_expression;
pub mod esurface;
pub mod expression;
pub mod generation;
pub mod hsurface;
pub mod surface;
pub mod tree;

pub struct CFFTerm {
    // One per orientation
    pub expression: Vec<Atom>,
    pub orientations: Vec<EdgeVec<Orientation>>,
}

impl CFFTerm {
    pub fn expression_with_selectors(&self) -> Atom {
        let mut result = Atom::Zero;
        for (expr, orient) in self.expression.iter().zip(self.orientations.iter()) {
            result *= expr.clone() * orient.orientation_thetas();
        }
        result
    }
}

pub struct CutCFF {
    pub terms: Vec<CFFTerm>, //index is the power of the esurface +1
}

impl CutCFF {
    pub fn expression_with_selectors(&self) -> Vec<Atom> {
        self.terms
            .iter()
            .map(|t| t.expression_with_selectors())
            .collect()
    }
}

impl Graph {
    pub fn cff<S: SubSetLike>(&mut self, contract_subgraph: &S, cutset: &CutSet) -> Result<CutCFF> {
        let canonize_esurface = self.get_esurface_canonization(&self.loop_momentum_basis);
        let mut contract_edges = vec![];

        for (p, eid, _) in self.iter_edges_of(contract_subgraph) {
            if p.is_paired() {
                contract_edges.push(eid);
            }
        }

        let cff = self.generate_cff(&contract_edges, &canonize_esurface)?;
        let residue = cff.select_esurface_residue(&cutset.esurfaces);

        let mut terms = vec![];

        let replacement_rules = self.surface_cache.get_all_replacements(&[]);

        for expr in residue.into_iter() {
            let mut cff_term = CFFTerm {
                expression: vec![],
                orientations: vec![],
            };
            for orientation in expr.orientations.iter() {
                let eta_expr = orientation.expression.to_atom_inv();
                let ose_expr = eta_expr.replace_multiple(&replacement_rules);
                cff_term.expression.push(ose_expr);
                cff_term
                    .orientations
                    .push(orientation.data.orientation.clone());
            }
            terms.push(cff_term);
        }

        let cut_cff = CutCFF { terms };
        Ok(cut_cff)
    }
}
