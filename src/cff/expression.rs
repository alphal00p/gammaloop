use crate::{
    graph::{EdgeType, Graph},
    utils::{FloatLike, F},
};
use derive_more::{From, Into};
use itertools::Itertools;
use log::debug;
use serde::{Deserialize, Serialize};
use std::{fmt::Debug, ops::Index};
use symbolica::{
    atom::{Atom, AtomView},
    domains::{float::NumericalFloatLike, rational::Rational},
    evaluate::{
        CompileOptions, CompiledCode, CompiledEvaluator, EvalTree, ExpressionEvaluator, FunctionMap,
    },
};
use typed_index_collections::TiVec;

#[derive(Debug, From, Into, Copy, Clone, Serialize, Deserialize)]
pub struct TermId(usize);

use super::{
    cff_graph::CFFGenerationGraph,
    esurface::{
        compute_esurface_cache, Esurface, EsurfaceCache, EsurfaceCollection, EsurfaceID,
        ExternalShift,
    },
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
    pub compiled: Option<CompiledCFFExpression>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum CFFExpressionNode {
    Data(HybridSurfaceID),
    Pointer { term_id: TermId, node_id: NodeId },
}

impl CFFExpression {
    #[inline]
    pub fn evaluate_orientations<T: FloatLike>(
        &self,
        energy_cache: &[F<T>],
        debug: usize,
    ) -> Vec<F<T>> {
        if debug > 3 {
            println!("evaluating cff orientations in eager mode");
        }

        let esurface_cache = self.compute_esurface_cache(energy_cache);
        let hsurface_cache = self.compute_hsurface_cache(energy_cache);

        let mut term_cache = self.build_term_cache();

        self.iter_term_ids()
            .map(|t| {
                let orientation_result = evaluate_tree(
                    &self.orientations[t].expression,
                    t,
                    &esurface_cache,
                    &hsurface_cache,
                    &mut term_cache,
                );

                if debug > 3 {
                    println!(
                        "result of orientation {:?}, {}",
                        self.orientations[t].orientation, orientation_result,
                    )
                }

                orientation_result
            })
            .collect()
    }

    #[inline]
    pub fn evaluate<T: FloatLike>(&self, energy_cache: &[F<T>], debug: usize) -> F<T> {
        self.evaluate_orientations(energy_cache, debug)
            .into_iter()
            .reduce(|acc, x| &acc + &x)
            .unwrap_or_else(|| energy_cache[0].zero())
    }

    #[inline]
    pub fn evaluate_orientations_from_caches<T: FloatLike>(
        &self,
        esurface_cache: &EsurfaceCache<F<T>>,
        hsurface_cache: &HsurfaceCache<F<T>>,
    ) -> Vec<F<T>> {
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
        esurface_cache: &EsurfaceCache<F<T>>,
        energy_cache: &[F<T>],
    ) -> Vec<F<T>> {
        let hsurface_cache = compute_hsurface_cache(&self.hsurfaces, energy_cache);
        self.evaluate_orientations_from_caches(esurface_cache, &hsurface_cache)
    }

    #[inline]
    pub fn evaluate_from_caches<T: FloatLike>(
        &self,
        esurface_cache: &EsurfaceCache<F<T>>,
        hsurface_cache: &HsurfaceCache<F<T>>,
    ) -> F<T> {
        self.evaluate_orientations_from_caches(esurface_cache, hsurface_cache)
            .into_iter()
            .reduce(|acc, x| &acc + &x)
            .unwrap_or_else(|| esurface_cache[EsurfaceID::from(0usize)].zero())
    }

    #[inline]
    pub fn evalauate_from_esurface_cache<T: FloatLike>(
        &self,
        esurface_cache: &EsurfaceCache<F<T>>,
        energy_cache: &[F<T>],
    ) -> F<T> {
        self.evaluate_orientations_from_esurface_cache(esurface_cache, energy_cache)
            .into_iter()
            .reduce(|acc, x| &acc + &x)
            .unwrap_or_else(|| energy_cache[0].zero())
    }

    #[inline]
    pub fn compute_esurface_cache<T: FloatLike>(
        &self,
        energy_cache: &[F<T>],
    ) -> EsurfaceCache<F<T>> {
        compute_esurface_cache(&self.esurfaces, energy_cache)
    }

    #[inline]
    pub fn compute_hsurface_cache<T: FloatLike>(
        &self,
        energy_cache: &[F<T>],
    ) -> HsurfaceCache<F<T>> {
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
    fn build_term_cache<T: FloatLike>(&self) -> TermCache<Option<F<T>>> {
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

    pub fn limit_for_esurface(
        &self,
        esurface_id: EsurfaceID,
        temp_dep_mom: usize,
        temp_dep_mom_expr: &ExternalShift,
    ) -> Result<CFFLimit, String> {
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

        // let (dep_mom, dep_mom_expr) = todo!("construct these objects from somewhere");

        generate_cff_limit(
            dag_left,
            dag_right,
            &self.esurfaces,
            ref_to_esurface,
            temp_dep_mom,
            temp_dep_mom_expr,
        )
    }

    pub fn build_symbolica_evaluators<T: FloatLike + Default>(
        &self,
        graph: &Graph,
    ) -> Vec<ExpressionEvaluator<F<T>>> {
        let function_map = FunctionMap::new();

        let params = graph
            .edges
            .iter()
            .enumerate()
            .map(|(id, edge)| match edge.edge_type {
                EdgeType::Virtual => Atom::parse(&format!("E{}", id)).unwrap(),
                _ => Atom::parse(&format!("p{}", id)).unwrap(),
            })
            .collect_vec();

        self.orientations
            .iter_enumerated()
            .map(|(term_id, _)| {
                let atom = self.construct_atom_for_term(term_id, None);
                let atom_view = atom.as_view();

                let mut tree = atom_view
                    .to_eval_tree(|r| r.clone(), &function_map, &params)
                    .unwrap();

                tree.horner_scheme();
                tree.common_subexpression_elimination();
                tree.common_pair_elimination();

                let tree_ft = tree.map_coeff::<F<T>, _>(&|r| r.into());
                tree_ft.linearize(params.len())
            })
            .collect()
    }

    /// graph currently needed to determine which edges are external
    pub fn build_joint_symbolica_evaluator<T: FloatLike + Default>(
        &self,
        graph: &Graph,
    ) -> ExpressionEvaluator<F<T>> {
        let orientation_atoms = self
            .orientations
            .iter_enumerated()
            .map(|(term_id, _)| self.construct_atom_for_term(term_id, None))
            .collect_vec();

        let orientation_atom_views = orientation_atoms.iter().map(Atom::as_view).collect_vec();
        let function_map = FunctionMap::new();

        let params = graph
            .edges
            .iter()
            .enumerate()
            .map(|(id, edge)| match edge.edge_type {
                EdgeType::Virtual => Atom::parse(&format!("E{}", id)).unwrap(),
                _ => Atom::parse(&format!("p{}", id)).unwrap(),
            })
            .collect_vec();

        let mut tree: EvalTree<Rational> = AtomView::to_eval_tree_multiple(
            &orientation_atom_views,
            |r| r.clone(),
            &function_map,
            &params,
        )
        .unwrap();

        debug!("optimizing cff");
        tree.horner_scheme();
        debug!("horner scheme completed");
        tree.common_subexpression_elimination();
        debug!("common subexpression elimination completed");
        tree.common_pair_elimination();
        debug!("common pair elimination completed");

        let tree_ft = tree.map_coeff::<F<T>, _>(&|r| r.into());
        tree_ft.linearize(graph.edges.len())
    }

    pub fn build_compiled_experssion(&self, graph: &Graph) -> CompiledCFFExpression {
        //let joint = self.build_joint_symbolica_evaluator(graph);
        //let orientations = self.build_symbolica_evaluators(graph);
        todo!();
    }
}

fn recursive_eval_from_node<T: FloatLike>(
    tree: &Tree<CFFExpressionNode>,
    node_id: NodeId,
    term_id: TermId,
    esurface_cache: &EsurfaceCache<F<T>>,
    hsurface_cache: &HsurfaceCache<F<T>>,
    term_cache: &mut TermCache<Option<F<T>>>,
) -> F<T> {
    let node = tree.get_node(node_id);
    match node.data {
        CFFExpressionNode::Pointer { term_id, node_id } => {
            term_cache[term_id][node_id].clone().unwrap()
        }
        CFFExpressionNode::Data(surface_id) => {
            let surface_val = match surface_id {
                HybridSurfaceID::Esurface(esurface_id) => esurface_cache[esurface_id].inv(),
                HybridSurfaceID::Hsurface(hsurface_id) => hsurface_cache[hsurface_id].inv(),
                HybridSurfaceID::Unit => esurface_cache[EsurfaceID::from(0usize)].one(),
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
                        .reduce(|acc, x| acc + x)
                        .unwrap_or_else(|| unreachable!())
            } else {
                surface_val
            };
            term_cache[term_id][node_id] = Some(res.clone());
            res
        }
    }
}

fn evaluate_tree<T: FloatLike>(
    tree: &Tree<CFFExpressionNode>,
    term_id: TermId,
    esurface_cache: &EsurfaceCache<F<T>>,
    hsurface_cache: &HsurfaceCache<F<T>>,
    term_cache: &mut TermCache<Option<F<T>>>,
) -> F<T> {
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
        esurface_cache: &EsurfaceCache<F<T>>,
        energy_cache: &[F<T>],
    ) -> F<T> {
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
            .reduce(|acc, x| &acc + &x)
            .unwrap_or_else(|| esurface_cache[EsurfaceID::from(0usize)].zero())
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

pub struct CompiledCFFExpression {
    metadata: CompiledCFFExpressionMetaData,
    joint: CompiledEvaluator,
    orientations: TiVec<TermId, CompiledEvaluator>,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
struct CompiledCFFExpressionMetaData {
    name: String,
    num_orientations: usize,
}

impl CompiledCFFExpression {
    pub fn evaluate_orientations(&self, energy_cache: &[F<f64>]) -> Vec<F<f64>> {
        let mut out = vec![F(0.0); self.orientations.len()];
        self.joint.evaluate(energy_cache, &mut out);
        out
    }

    pub fn evaluate_one_orientation(&self, orientation: TermId, energy_cache: &[F<f64>]) -> F<f64> {
        let mut out = [F(0.0)];
        self.orientations[orientation].evaluate(energy_cache, &mut out);
        out[0]
    }

    fn from_metedata(metadata: CompiledCFFExpressionMetaData) -> Self {
        let joint = CompiledEvaluator::load(&format!("{}_joint", metadata.name)).unwrap();
        let orientations = (0..metadata.num_orientations)
            .map(|orientation| {
                CompiledEvaluator::load(&format!("{}_orientation_{}", metadata.name, orientation))
                    .unwrap()
            })
            .collect();

        Self {
            metadata,
            joint,
            orientations,
        }
    }
}

impl Clone for CompiledCFFExpression {
    fn clone(&self) -> Self {
        Self::from_metedata(self.metadata.clone())
    }
}

impl Serialize for CompiledCFFExpression {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        CompiledCFFExpressionMetaData::serialize(&self.metadata, serializer)
    }
}

impl<'de> Deserialize<'de> for CompiledCFFExpression {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let metadata = CompiledCFFExpressionMetaData::deserialize(deserializer)?;
        Ok(Self::from_metedata(metadata))
    }
}

impl Debug for CompiledCFFExpression {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        CompiledCFFExpressionMetaData::fmt(&self.metadata, f)
    }
}
