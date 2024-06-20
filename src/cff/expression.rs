use crate::utils::FloatLike;
use derive_more::{From, Into};
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::ops::Index;
use symbolica::representations::Atom;
use typed_index_collections::TiVec;

#[derive(Debug, From, Into, Copy, Clone, Serialize, Deserialize)]
pub struct TermId(usize);

use super::{
    cff_graph::CFFGenerationGraph,
    esurface::{compute_esurface_cache, Esurface, EsurfaceCache, EsurfaceCollection, EsurfaceID},
    generation::generate_cff_limit,
    hsurface::{compute_hsurface_cache, HsurfaceCache, HsurfaceCollection},
    surface::HybridSurfaceID,
    tree::{NodeCache, NodeId, Tree},
};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OrientationExpression {
    pub orientation: Vec<bool>,
    pub dag: CFFGenerationGraph,
    pub expression: Tree<CFFExpressionNode>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CFFExpression {
    pub orientations: TiVec<TermId, OrientationExpression>,
    pub esurfaces: EsurfaceCollection,
    pub hsurfaces: HsurfaceCollection,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum CFFExpressionNode {
    Data(HybridSurfaceID),
    Pointer { term_id: TermId, node_id: NodeId },
}

impl CFFExpression {
    #[inline]
    pub fn evaluate_orientations<T: FloatLike>(&self, energy_cache: &[T]) -> Vec<T> {
        let esurface_cache = self.compute_esurface_cache(energy_cache);
        let hsurface_cache = self.compute_hsurface_cache(energy_cache);

        let mut term_cache = self.build_term_cache();

        self.iter_term_ids()
            .map(|t| {
                evaluate_tree(
                    &self.orientations[t].expression,
                    t,
                    &esurface_cache,
                    &hsurface_cache,
                    &mut term_cache,
                )
            })
            .collect()
    }

    #[inline]
    pub fn evaluate<T: FloatLike>(&self, energy_cache: &[T]) -> T {
        self.evaluate_orientations(energy_cache)
            .into_iter()
            .sum::<T>()
    }

    #[inline]
    pub fn evaluate_orientations_from_caches<T: FloatLike>(
        &self,
        esurface_cache: &EsurfaceCache<T>,
        hsurface_cache: &HsurfaceCache<T>,
    ) -> Vec<T> {
        let mut term_cache = self.build_term_cache();

        self.iter_term_ids()
            .map(|t| {
                evaluate_tree(
                    &self.orientations[t].expression,
                    t,
                    esurface_cache,
                    hsurface_cache,
                    &mut term_cache,
                )
            })
            .collect()
    }

    #[inline]
    pub fn evaluate_orientations_from_esurface_cache<T: FloatLike>(
        &self,
        esurface_cache: &EsurfaceCache<T>,
        energy_cache: &[T],
    ) -> Vec<T> {
        let hsurface_cache = compute_hsurface_cache(&self.hsurfaces, energy_cache);
        self.evaluate_orientations_from_caches(esurface_cache, &hsurface_cache)
    }

    #[inline]
    pub fn evaluate_from_caches<T: FloatLike>(
        &self,
        esurface_cache: &EsurfaceCache<T>,
        hsurface_cache: &HsurfaceCache<T>,
    ) -> T {
        self.evaluate_orientations_from_caches(esurface_cache, hsurface_cache)
            .into_iter()
            .sum::<T>()
    }

    #[inline]
    pub fn evalauate_from_esurface_cache<T: FloatLike>(
        &self,
        esurface_cache: &EsurfaceCache<T>,
        energy_cache: &[T],
    ) -> T {
        self.evaluate_orientations_from_esurface_cache(esurface_cache, energy_cache)
            .into_iter()
            .sum::<T>()
    }

    #[inline]
    pub fn compute_esurface_cache<T: FloatLike>(&self, energy_cache: &[T]) -> EsurfaceCache<T> {
        compute_esurface_cache(&self.esurfaces, energy_cache)
    }

    #[inline]
    pub fn compute_hsurface_cache<T: FloatLike>(&self, energy_cache: &[T]) -> HsurfaceCache<T> {
        compute_hsurface_cache(&self.hsurfaces, energy_cache)
    }

    /// Helper function to undo the factorisation
    fn recursive_term_builder(
        &self,
        res: &mut Vec<Vec<HybridSurfaceID>>,
        current_path: &mut Vec<HybridSurfaceID>,
        term_id: TermId,
        node_id: NodeId,
    ) {
        let node = &self.orientations[term_id].expression.get_node(node_id);

        match node.data {
            CFFExpressionNode::Data(surface_id) => {
                current_path.push(surface_id);

                if node.children.is_empty() {
                    res.push(current_path.clone());
                } else {
                    for child in node.children.iter() {
                        self.recursive_term_builder(
                            res,
                            &mut current_path.clone(),
                            term_id,
                            *child,
                        );
                    }
                }
            }
            CFFExpressionNode::Pointer { term_id, node_id } => {
                self.recursive_term_builder(res, current_path, term_id, node_id);
            }
        }
    }

    /// Undo the factorisationo of a given term
    fn expand_term(&self, term_id: TermId) -> Vec<Vec<HybridSurfaceID>> {
        let mut res = vec![];
        let mut current_path = vec![];

        self.recursive_term_builder(&mut res, &mut current_path, term_id, NodeId::root());
        res
    }

    /// Undo the factorisation of the whole expression of the whole expression
    pub fn expand_terms(&self) -> Vec<Vec<HybridSurfaceID>> {
        self.iter_term_ids()
            .flat_map(|t| self.expand_term(t))
            .collect()
    }

    /// Take the residue of the expression at the surfaec
    pub fn expand_limit(&self, surface_id: HybridSurfaceID) -> Vec<Vec<HybridSurfaceID>> {
        let mut terms = self.expand_terms();
        terms.retain(|product| product.contains(&surface_id));

        for term in terms.iter_mut() {
            term.retain(|&term_surface_id| term_surface_id != surface_id)
        }

        terms
    }

    /// Take the residue of the expression at the surface and convert it to a Symbolica expression
    pub fn expand_limit_to_atom(&self, surface_id: HybridSurfaceID) -> Atom {
        let limit = self.expand_limit(surface_id);

        limit
            .iter()
            .map(|product| {
                let product_atoms = product.iter().map(|surface| match surface {
                    HybridSurfaceID::Esurface(esurface_id) => {
                        self.esurfaces[*esurface_id].to_atom()
                    }
                    HybridSurfaceID::Hsurface(hsurface_id) => {
                        self.hsurfaces[*hsurface_id].to_atom()
                    }
                    HybridSurfaceID::Unit => Atom::new_num(1),
                });

                let surface_product =
                    product_atoms.fold(Atom::new_num(1), |atom_prod, e| atom_prod * &e);

                Atom::new_num(1) / &surface_product
            })
            .fold(Atom::new(), |sum, e| sum + &e)
    }

    pub fn iter_term_ids(&self) -> impl Iterator<Item = TermId> {
        (0..self.orientations.len()).map(TermId)
    }

    #[inline]
    fn build_term_cache<T: FloatLike>(&self) -> TermCache<Option<T>> {
        self.orientations
            .iter()
            .map(|orientation| vec![None; orientation.expression.get_num_nodes()].into())
            .collect_vec()
            .into()
    }

    #[inline]
    pub fn get_num_trees(&self) -> usize {
        self.orientations.len()
    }

    fn recursive_esurface_search(
        &self,
        term_id: TermId,
        node_id: NodeId,
        esurface_id: EsurfaceID,
    ) -> bool {
        let node = self.orientations[term_id].expression.get_node(node_id);

        match node.data {
            CFFExpressionNode::Data(data_surface_id) => {
                match data_surface_id {
                    HybridSurfaceID::Esurface(data_esurface_id) => {
                        if data_esurface_id == esurface_id {
                            return true;
                        }
                    }
                    HybridSurfaceID::Hsurface(_) => {}
                    HybridSurfaceID::Unit => {}
                }

                if node.children.is_empty() {
                    return false;
                }

                node.children
                    .iter()
                    .any(|child| self.recursive_esurface_search(term_id, *child, esurface_id))
            }
            CFFExpressionNode::Pointer { term_id, node_id } => {
                self.recursive_esurface_search(term_id, node_id, esurface_id)
            }
        }
    }

    pub fn term_has_esurface(&self, term_id: TermId, esurface_id: EsurfaceID) -> bool {
        self.recursive_esurface_search(term_id, NodeId::root(), esurface_id)
    }

    fn recursive_construct_atom(
        &self,
        node_id: NodeId,
        term_id: TermId,
        esurface_collection: &EsurfaceCollection,
        hsurface_collection: &HsurfaceCollection,
        rewriter_esurface: Option<&Esurface>,
    ) -> Atom {
        let tree = &self[term_id].expression;
        let node = tree.get_node(node_id);
        match node.data {
            CFFExpressionNode::Pointer { term_id, node_id } => self.recursive_construct_atom(
                node_id,
                term_id,
                esurface_collection,
                hsurface_collection,
                rewriter_esurface,
            ),
            CFFExpressionNode::Data(surface_id) => {
                let surface_atom = match surface_id {
                    HybridSurfaceID::Esurface(esurface_id) => {
                        esurface_collection[esurface_id].to_atom()
                    }
                    HybridSurfaceID::Hsurface(hsurface_id) => match rewriter_esurface {
                        Some(esurface) => {
                            if let Some(esurface_atom) =
                                hsurface_collection[hsurface_id].to_atom_with_rewrite(esurface)
                            {
                                esurface_atom
                            } else {
                                hsurface_collection[hsurface_id].to_atom()
                            }
                        }
                        None => hsurface_collection[hsurface_id].to_atom(),
                    },
                    HybridSurfaceID::Unit => Atom::new_num(1),
                };

                let inv_surface_atom = Atom::new_num(1) / &surface_atom;

                let res = if !node.children.is_empty() {
                    inv_surface_atom
                        * &node
                            .children
                            .iter()
                            .map(|child_index| {
                                self.recursive_construct_atom(
                                    *child_index,
                                    term_id,
                                    esurface_collection,
                                    hsurface_collection,
                                    rewriter_esurface,
                                )
                            })
                            .fold(Atom::new(), |sum, e| sum + &e)
                } else {
                    inv_surface_atom
                };
                res
            }
        }
    }

    pub fn construct_atom_for_term(
        &self,
        term_id: TermId,
        rewriter_esurface: Option<&Esurface>,
    ) -> Atom {
        self.recursive_construct_atom(
            NodeId::root(),
            term_id,
            &self.esurfaces,
            &self.hsurfaces,
            rewriter_esurface,
        )
    }

    pub fn limit_for_esurface(&self, esurface_id: EsurfaceID) -> Result<CFFLimit, String> {
        let circling = self.esurfaces[esurface_id].circled_vertices;

        let terms_with_esurface = self
            .iter_term_ids()
            .filter(|&term_id| self.term_has_esurface(term_id, esurface_id));

        let (dag_left, dag_right) = terms_with_esurface
            .map(|term_id| {
                let term_dag = &self[term_id].dag;
                term_dag.generate_cut(circling)
            })
            .unzip();

        let ref_to_esurface = &self.esurfaces[esurface_id];

        generate_cff_limit(dag_left, dag_right, &self.esurfaces, ref_to_esurface)
    }
}

fn recursive_eval_from_node<T: FloatLike>(
    tree: &Tree<CFFExpressionNode>,
    node_id: NodeId,
    term_id: TermId,
    esurface_cache: &EsurfaceCache<T>,
    hsurface_cache: &HsurfaceCache<T>,
    term_cache: &mut TermCache<Option<T>>,
) -> T {
    let node = tree.get_node(node_id);
    match node.data {
        CFFExpressionNode::Pointer { term_id, node_id } => term_cache[term_id][node_id].unwrap(),
        CFFExpressionNode::Data(surface_id) => {
            let surface_val = match surface_id {
                HybridSurfaceID::Esurface(esurface_id) => esurface_cache[esurface_id].inv(),
                HybridSurfaceID::Hsurface(hsurface_id) => hsurface_cache[hsurface_id].inv(),
                HybridSurfaceID::Unit => T::one(),
            };

            let res = if !node.children.is_empty() {
                surface_val
                    * node
                        .children
                        .iter()
                        .map(|child_index| {
                            recursive_eval_from_node(
                                tree,
                                *child_index,
                                term_id,
                                esurface_cache,
                                hsurface_cache,
                                term_cache,
                            )
                        })
                        .sum::<T>()
            } else {
                surface_val
            };
            term_cache[term_id][node_id] = Some(res);
            res
        }
    }
}

fn evaluate_tree<T: FloatLike>(
    tree: &Tree<CFFExpressionNode>,
    term_id: TermId,
    esurface_cache: &EsurfaceCache<T>,
    hsurface_cache: &HsurfaceCache<T>,
    term_cache: &mut TermCache<Option<T>>,
) -> T {
    recursive_eval_from_node(
        tree,
        NodeId::root(),
        term_id,
        esurface_cache,
        hsurface_cache,
        term_cache,
    )
}

pub type TermCache<T> = TiVec<TermId, NodeCache<T>>;

impl Index<TermId> for CFFExpression {
    type Output = OrientationExpression;

    fn index(&self, term_id: TermId) -> &Self::Output {
        &self.orientations[term_id]
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CFFLimit {
    pub left: CFFExpression,
    pub right: CFFExpression,
}

impl CFFLimit {
    pub fn evaluate_from_esurface_cache<T: FloatLike>(
        &self,
        esurface_cache: &EsurfaceCache<T>,
        energy_cache: &[T],
    ) -> T {
        let left_orientations = self
            .left
            .evaluate_orientations_from_esurface_cache(esurface_cache, energy_cache);
        let right_orientations = self
            .right
            .evaluate_orientations_from_esurface_cache(esurface_cache, energy_cache);

        left_orientations
            .into_iter()
            .zip(right_orientations)
            .map(|(l, r)| l * r)
            .sum()
    }

    pub fn limit_to_atom_with_rewrite(&self, rewriter_esurface: Option<&Esurface>) -> Atom {
        let left_expressions = self
            .left
            .orientations
            .iter_enumerated()
            .map(|(term_id, _)| {
                self.left
                    .construct_atom_for_term(term_id, rewriter_esurface)
            });

        let right_expressions = self
            .right
            .orientations
            .iter_enumerated()
            .map(|(term_id, _)| {
                self.right
                    .construct_atom_for_term(term_id, rewriter_esurface)
            });

        left_expressions
            .zip(right_expressions)
            .fold(Atom::new(), |res, (left, right)| res + &(left * &right))
    }
}

pub enum HybridNode {
    Data(HybridSurfaceID),
    Pointer { term_id: TermId, node_id: NodeId },
}
