use std::ops::{Index, IndexMut};

use serde::{Deserialize, Serialize};

use crate::utils::FloatLike;

#[derive(Debug, Copy, Clone, Serialize, Deserialize)]
pub struct TermId(usize);

impl From<usize> for TermId {
    fn from(value: usize) -> Self {
        TermId(value)
    }
}

use super::{
    esurface::{EsurfaceCache, EsurfaceCollection, EsurfaceId},
    tree::{NodeCache, NodeId, Tree},
};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CFFExpression {
    pub expression: Vec<Tree<CFFExpressionNode>>,
    pub orientations: Vec<Vec<bool>>,
    pub esurfaces: EsurfaceCollection,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum CFFExpressionNode {
    Data(EsurfaceId),
    Pointer { term_id: TermId, node_id: NodeId },
}

impl CFFExpression {
    #[inline]
    pub fn evaluate_orientations<T: FloatLike>(&self, energy_cache: &[T]) -> Vec<T> {
        let esurface_cache = self.compute_esurface_cache(energy_cache);
        let mut term_cache = self.build_term_cache();

        self.iter_term_ids()
            .map(|t| evaluate_tree(&self.expression[t.0], t, &esurface_cache, &mut term_cache))
            .collect()
    }

    #[inline]
    pub fn evaluate<T: FloatLike>(&self, energy_cache: &[T]) -> T {
        self.evaluate_orientations(energy_cache)
            .into_iter()
            .sum::<T>()
    }

    #[inline]
    pub fn compute_esurface_cache<T: FloatLike>(&self, energy_cache: &[T]) -> EsurfaceCache<T> {
        self.esurfaces.compute_esurface_cache(energy_cache)
    }

    fn recursive_term_builder(
        &self,
        res: &mut Vec<Vec<EsurfaceId>>,
        current_path: &mut Vec<EsurfaceId>,
        term_id: TermId,
        node_id: NodeId,
    ) {
        let node = &self.expression[term_id.0].get_node(node_id);

        match node.data {
            CFFExpressionNode::Data(esurface_id) => {
                current_path.push(esurface_id);

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

    fn expand_term(&self, term_id: TermId) -> Vec<Vec<EsurfaceId>> {
        let mut res = vec![];
        let mut current_path = vec![];

        self.recursive_term_builder(&mut res, &mut current_path, term_id, NodeId::root());
        res
    }

    pub fn expand_terms(&self) -> Vec<Vec<EsurfaceId>> {
        self.iter_term_ids()
            .flat_map(|t| self.expand_term(t))
            .collect()
    }

    fn iter_term_ids(&self) -> impl Iterator<Item = TermId> {
        (0..self.expression.len()).map(TermId)
    }

    #[inline]
    fn build_term_cache<T: FloatLike>(&self) -> TermCache<Option<T>> {
        let cache = self
            .expression
            .iter()
            .map(|tree| NodeCache::<Option<T>>::new_default(tree.get_num_nodes(), None))
            .collect();

        TermCache { cache }
    }

    #[inline]
    pub fn get_num_trees(&self) -> usize {
        self.expression.len()
    }
}

fn recursive_eval_from_node<T: FloatLike>(
    tree: &Tree<CFFExpressionNode>,
    node_id: NodeId,
    term_id: TermId,
    esurface_cache: &EsurfaceCache<T>,
    term_cache: &mut TermCache<Option<T>>,
) -> T {
    let node = tree.get_node(node_id);
    match node.data {
        CFFExpressionNode::Pointer { term_id, node_id } => term_cache[(term_id, node_id)].unwrap(),
        CFFExpressionNode::Data(esurface_id) => {
            let res = if !node.children.is_empty() {
                esurface_cache[esurface_id].inv()
                    * node
                        .children
                        .iter()
                        .map(|child_index| {
                            recursive_eval_from_node(
                                tree,
                                *child_index,
                                term_id,
                                esurface_cache,
                                term_cache,
                            )
                        })
                        .sum::<T>()
            } else {
                esurface_cache[esurface_id].inv()
            };
            term_cache[(term_id, node_id)] = Some(res);
            res
        }
    }
}

fn evaluate_tree<T: FloatLike>(
    tree: &Tree<CFFExpressionNode>,
    term_id: TermId,
    esurface_cache: &EsurfaceCache<T>,
    term_cache: &mut TermCache<Option<T>>,
) -> T {
    recursive_eval_from_node(tree, NodeId::root(), term_id, esurface_cache, term_cache)
}

struct TermCache<T> {
    cache: Vec<NodeCache<T>>,
}

impl<T> Index<(TermId, NodeId)> for TermCache<T> {
    type Output = T;

    fn index(&self, (term_id, node_id): (TermId, NodeId)) -> &Self::Output {
        &self.cache[term_id.0][node_id]
    }
}

impl<T> IndexMut<(TermId, NodeId)> for TermCache<T> {
    fn index_mut(&mut self, (term_id, node_id): (TermId, NodeId)) -> &mut Self::Output {
        &mut self.cache[term_id.0][node_id]
    }
}
