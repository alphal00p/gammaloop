use crate::{
    debug_info::DEBUG_LOGGER,
    gammaloop_integrand::DefaultSample,
    graph::BareGraph,
    momentum::FourMomentum,
    numerator::{Evaluate, Evaluators, Numerator, RepeatingIteratorTensorOrScalar},
    utils::{FloatLike, VarFloat, F},
    ExportSettings, Settings,
};
use bincode::{Decode, Encode};
use color_eyre::Report;
use derive_more::{From, Into};
use eyre::eyre;
use gat_lending_iterator::LendingIterator;
use itertools::Itertools;
use log::info;
use serde::{Deserialize, Serialize};
use smartstring::{LazyCompact, SmartString};
use spenso::{complex::Complex, parametric::SerializableCompiledEvaluator};
use std::{cell::RefCell, fmt::Debug, ops::Index, path::PathBuf};
use symbolica::{
    atom::{Atom, AtomView},
    domains::{float::NumericalFloatLike, rational::Rational},
    evaluate::{EvalTree, ExportedCode, ExpressionEvaluator, FunctionMap, InlineASM},
};
use typed_index_collections::TiVec;

#[derive(Debug, From, Into, Copy, Clone, Serialize, Deserialize, Eq, PartialEq)]
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

pub trait CFFFloat<T: FloatLike> {
    fn get_evaluator(cff: &CFFExpression) -> impl Fn(&[F<T>], &Settings) -> Vec<F<T>>;
}

impl CFFFloat<f64> for f64 {
    fn get_evaluator(cff: &CFFExpression) -> impl Fn(&[F<f64>], &Settings) -> Vec<F<f64>> {
        |energy_cache, settings| {
            if cff.compiled.is_some() {
                cff.compiled_evaluate_orientations(energy_cache, settings)
            } else {
                cff.eager_evaluate_orientations(energy_cache, settings)
            }
        }
    }
}

impl CFFFloat<VarFloat<113>> for VarFloat<113> {
    fn get_evaluator(
        cff: &CFFExpression,
    ) -> impl Fn(&[F<VarFloat<113>>], &Settings) -> Vec<F<VarFloat<113>>> {
        |energy_cache, settings| cff.eager_evaluate_orientations(energy_cache, settings)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OrientationExpression {
    pub orientation: Vec<bool>,
    pub dag: CFFGenerationGraph,
    pub expression: Tree<CFFExpressionNode>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
pub struct CFFExpression {
    #[bincode(with_serde)]
    pub orientations: TiVec<TermId, OrientationExpression>,
    #[bincode(with_serde)]
    pub esurfaces: EsurfaceCollection,
    #[bincode(with_serde)]
    pub hsurfaces: HsurfaceCollection,
    #[bincode(with_serde)]
    pub compiled: CompiledCFFExpression,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum CFFExpressionNode {
    Data(HybridSurfaceID),
    Pointer { term_id: TermId, node_id: NodeId },
}

impl CFFExpression {
    #[inline]
    pub fn evaluate_orientations<T: FloatLike + CFFFloat<T>>(
        &self,
        energy_cache: &[F<T>],
        settings: &Settings,
    ) -> Vec<F<T>> {
        if settings.general.debug > 0 {
            DEBUG_LOGGER.write("esurface_equation", &self.esurfaces);
            DEBUG_LOGGER.write("onshell_energies", &energy_cache);
        }

        T::get_evaluator(self)(energy_cache, settings)
    }

    #[inline]
    fn eager_evaluate_orientations<T: FloatLike>(
        &self,
        energy_cache: &[F<T>],
        settings: &Settings,
    ) -> Vec<F<T>> {
        let esurface_cache = self.compute_esurface_cache(energy_cache);
        let hsurface_cache = self.compute_hsurface_cache(energy_cache);

        if settings.general.debug > 0 {
            DEBUG_LOGGER.write("esurface_values", &esurface_cache);
        }

        let mut term_cache = self.build_term_cache();

        let res = self
            .iter_term_ids()
            .map(|t| {
                evaluate_tree(
                    &self.orientations[t].expression,
                    t,
                    &esurface_cache,
                    &hsurface_cache,
                    &mut term_cache,
                )
            })
            .collect();

        if settings.general.debug > 0 {
            DEBUG_LOGGER.write("orientations", &res);
        }

        res
    }

    #[inline]
    fn compiled_evaluate_orientations(
        &self,
        energy_cache: &[F<f64>],
        settings: &Settings,
    ) -> Vec<F<f64>> {
        let res = self.compiled.evaluate_orientations(energy_cache, settings);

        if settings.general.debug > 2 {
            DEBUG_LOGGER.write("orientations", &res);

            let esurfaces_per_orientation = (0..self.get_num_trees())
                .map(TermId)
                .map(|orientation_id: TermId| {
                    self.esurfaces
                        .iter_enumerated()
                        .map(|(esurface_id, _)| esurface_id)
                        .filter(|esurface_id| self.term_has_esurface(orientation_id, *esurface_id))
                        .collect_vec()
                })
                .collect_vec();

            println!();
            for (orientation_id, esurfaces) in esurfaces_per_orientation.iter().enumerate() {
                println!("orientation: {}", orientation_id);
                let orientaation = &self.orientations[TermId(orientation_id)];
                let or = &orientaation.orientation;
                println!("esurfaces: {:?}", esurfaces);
                println!("or: {:?}", or);
                println!();
            }
        }

        res
    }

    #[inline]
    pub fn evaluate<T: FloatLike>(&self, energy_cache: &[F<T>], settings: &Settings) -> F<T> {
        self.evaluate_orientations(energy_cache, settings)
            .into_iter()
            .reduce(|acc, x| &acc + &x)
            .unwrap_or_else(|| energy_cache[0].zero())
    }

    #[inline]
    pub fn eager_evaluate<T: FloatLike>(&self, energy_cache: &[F<T>], settings: &Settings) -> F<T> {
        self.eager_evaluate_orientations(energy_cache, settings)
            .into_iter()
            .reduce(|acc, x| &acc + &x)
            .unwrap_or_else(|| energy_cache[0].zero())
    }

    #[inline]
    fn evaluate_orientations_from_caches<T: FloatLike>(
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

        let ((dag_left, dag_right), orientations_in_limit) = terms_with_esurface
            .map(|term_id| {
                let term_dag = &self[term_id].dag;
                (
                    term_dag.generate_cut(circling),
                    (self[term_id].orientation.clone(), term_id),
                )
            })
            .unzip();

        let ref_to_esurface = &self.esurfaces[esurface_id];

        generate_cff_limit(
            dag_left,
            dag_right,
            &self.esurfaces,
            ref_to_esurface,
            temp_dep_mom,
            temp_dep_mom_expr,
            orientations_in_limit,
        )
    }

    pub fn build_symbolica_evaluators<T: FloatLike + Default>(
        &self,
        params: &[Atom],
    ) -> Vec<ExpressionEvaluator<F<T>>> {
        let function_map = FunctionMap::new();

        self.orientations
            .iter_enumerated()
            .map(|(term_id, _)| {
                let atom = self.construct_atom_for_term(term_id, None);
                let atom_view = atom.as_view();

                let mut tree = atom_view.to_evaluation_tree(&function_map, params).unwrap();

                tree.horner_scheme();
                tree.common_subexpression_elimination();

                let tree_ft = tree.map_coeff::<F<T>, _>(&|r| r.into());
                tree_ft.linearize(Some(1))
            })
            .collect()
    }

    pub fn build_joint_symbolica_evaluator<T: FloatLike + Default>(
        &self,
        params: &[Atom],
        cpe_rounds: Option<usize>,
    ) -> ExpressionEvaluator<F<T>> {
        let orientation_atoms = self
            .orientations
            .iter_enumerated()
            .map(|(term_id, _)| self.construct_atom_for_term(term_id, None))
            .collect_vec();

        let orientation_atom_views = orientation_atoms.iter().map(Atom::as_view).collect_vec();
        let function_map = FunctionMap::new();

        let mut tree: EvalTree<Rational> =
            AtomView::to_eval_tree_multiple(&orientation_atom_views, &function_map, params)
                .unwrap();

        tree.horner_scheme();
        tree.common_subexpression_elimination();

        let tree_ft = tree.map_coeff::<F<T>, _>(&|r| r.into());
        tree_ft.linearize(cpe_rounds)
    }

    /// does nothing if compile_cff and compile_separate_orientations are both set to false
    pub fn build_compiled_expression<T: FloatLike + Default>(
        &mut self,
        params: &[Atom],
        path: PathBuf,
        graph_name: SmartString<LazyCompact>,
        export_settings: &ExportSettings,
    ) -> Result<(), Report> {
        if !export_settings.compile_cff && !export_settings.compile_separate_orientations {
            return Ok(());
        }

        let expr_str = format!("expression_{}", graph_name);
        let mut cpp_str = String::new();

        let path_to_compiled = path.join("compiled");
        std::fs::create_dir_all(&path_to_compiled)?;

        let path_to_code = path_to_compiled.join(format!("{}.cpp", expr_str));

        info!(
            "Compiling cff source_code {}",
            path_to_code
                .to_str()
                .ok_or(eyre!("could not convert path to string"))?
        );

        let path_to_so = path_to_compiled.join(format!("{}.so", expr_str));
        let path_to_so_str = path_to_so
            .to_str()
            .ok_or(eyre!("could not convert path to string"))?;

        let path_to_code_clone = path_to_code.clone();
        let path_to_code_str = path_to_code_clone
            .to_str()
            .ok_or(eyre!("could not convert path to string"))?;

        if export_settings.compile_cff {
            let joint =
                self.build_joint_symbolica_evaluator::<T>(params, export_settings.cpe_rounds_cff);

            let source_string = if export_settings.gammaloop_compile_options.inline_asm {
                joint.export_asm_str("joint", true, InlineASM::default())
            } else {
                joint.export_cpp_str("joint", true)
            };

            cpp_str.push_str(&source_string);
        }

        if export_settings.compile_separate_orientations {
            let orientations = self.build_symbolica_evaluators::<T>(params);
            for (orientation_id, orientation_evaluator) in orientations.into_iter().enumerate() {
                let orientation_source_str = if export_settings.gammaloop_compile_options.inline_asm
                {
                    orientation_evaluator.export_asm_str(
                        &format!("orientation_{}", orientation_id),
                        !export_settings.compile_cff && orientation_id == 0,
                        InlineASM::X64,
                    )
                } else {
                    orientation_evaluator.export_cpp_str(
                        &format!("orientation_{}", orientation_id),
                        !export_settings.compile_cff && orientation_id == 0,
                    )
                };

                cpp_str.push_str(&orientation_source_str);
            }
        }

        std::fs::write(path_to_code, cpp_str)?;

        let exported_code = ExportedCode::new(
            path_to_code_str.to_string(),
            "joint".to_string(),
            // export_settings.gammaloop_compile_options.inline_asm(),
        );
        exported_code.compile(
            path_to_so_str,
            export_settings
                .gammaloop_compile_options
                .to_symbolica_compile_options(),
        )?;

        let metadata = CompiledCFFExpressionMetaData {
            name: path_to_compiled,
            graph_name,
            num_orientations: self.get_num_trees(),
            compile_cff_present: export_settings.compile_cff,
            compile_separate_orientations_present: export_settings.compile_separate_orientations,
        };

        self.compiled = CompiledCFFExpression::from_metedata(metadata)?;

        info!("Compilation successful");

        Ok(())
    }

    pub fn load_compiled(
        &mut self,
        path: PathBuf,
        graph_name: SmartString<LazyCompact>,
        settings: &Settings,
    ) -> Result<(), Report> {
        let metadata = CompiledCFFExpressionMetaData {
            name: path.join("compiled"),
            graph_name,
            num_orientations: self.get_num_trees(),
            compile_cff_present: settings.general.load_compiled_cff,
            compile_separate_orientations_present: settings
                .general
                .load_compiled_separate_orientations,
        };

        if !settings.general.load_compiled_cff
            && !settings.general.load_compiled_separate_orientations
        {
            self.compiled = CompiledCFFExpression::None;
            return Ok(());
        }

        self.compiled = CompiledCFFExpression::from_metedata(metadata)?;
        Ok(())
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
    pub orientations_in_limit: (Vec<Vec<bool>>, Vec<TermId>),
}

impl CFFLimit {
    pub fn evaluate_from_esurface_cache<T: FloatLike>(
        &self,
        graph: &BareGraph,
        numerator: &mut Numerator<Evaluators>,
        numerator_sample: &DefaultSample<T>,
        esurface_cache: &EsurfaceCache<F<T>>,
        energy_cache: &[F<T>],
        settings: &Settings,
    ) -> Complex<F<T>> {
        let (numerator_sample, tag) = numerator_sample.numerator_sample(settings);

        let emr_energies = graph
            .compute_onshell_energies(&numerator_sample.loop_moms, &numerator_sample.external_moms);

        let emr_3d =
            graph.compute_emr(&numerator_sample.loop_moms, &numerator_sample.external_moms);

        let emr_4d = emr_energies
            .into_iter()
            .zip(emr_3d)
            .map(|(e, p)| FourMomentum::from_args(e, p.px, p.py, p.pz))
            .collect_vec();

        let left_orientations = self
            .left
            .evaluate_orientations_from_esurface_cache(esurface_cache, energy_cache);
        let right_orientations = self
            .right
            .evaluate_orientations_from_esurface_cache(esurface_cache, energy_cache);

        let num_iter = numerator
            .evaluate_all_orientations(&emr_4d, &numerator_sample.polarizations, tag, settings)
            .unwrap();

        let mut cff = left_orientations
            .into_iter()
            .zip(right_orientations)
            .map(|(l, r)| l * r);

        match num_iter {
            RepeatingIteratorTensorOrScalar::Scalars(mut num) => {
                let mut term = 0;
                let mut _terms_evaluated = 0;
                let mut sum = Complex::new_re(energy_cache[0].zero());

                while let Some(num) = num.next() {
                    let term_in_residue = self.orientations_in_limit.1.contains(&TermId(term));
                    let term_in_force_orientation = settings
                        .general
                        .force_orientations
                        .as_ref()
                        .map(|forced_orientations| forced_orientations.contains(&term))
                        .unwrap_or(true);

                    if term_in_residue && term_in_force_orientation {
                        let cff_term = Complex::new_re(cff.next().unwrap());
                        sum += num * cff_term;
                        _terms_evaluated += 1;
                    }
                    term += 1;
                }
                sum
            }
            RepeatingIteratorTensorOrScalar::Tensors(mut _num_iter) => {
                todo!()
            }
        }
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

// custom option so we have control over serialize/deserialize
#[derive(Clone, Serialize, Deserialize)]
#[allow(clippy::large_enum_variant)]
pub enum CompiledCFFExpression {
    // #[serde(skip)]
    Some(InnerCompiledCFF),
    None,
}

impl Default for CompiledCFFExpression {
    fn default() -> Self {
        Self::None
    }
}

#[derive(Clone, Serialize, Deserialize)]
pub struct InnerCompiledCFF {
    metadata: CompiledCFFExpressionMetaData,
    joint: Option<RefCell<SerializableCompiledEvaluator>>,
    orientations: TiVec<TermId, RefCell<SerializableCompiledEvaluator>>,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
struct CompiledCFFExpressionMetaData {
    name: PathBuf,
    graph_name: SmartString<LazyCompact>,
    num_orientations: usize,
    compile_cff_present: bool,
    compile_separate_orientations_present: bool,
}

impl CompiledCFFExpression {
    pub fn evaluate_orientations(
        &self,
        energy_cache: &[F<f64>],
        settings: &Settings,
    ) -> Vec<F<f64>> {
        let expr = self.unwrap();
        let mut out = vec![F(0.0); expr.metadata.num_orientations];
        match &expr.joint {
            Some(evaluator) => evaluator.borrow_mut().evaluate(energy_cache, &mut out),
            None => match &settings.general.force_orientations {
                None => {
                    for (id, out_elem) in out.iter_mut().enumerate() {
                        *out_elem = self.evaluate_one_orientation(id.into(), energy_cache);
                    }
                }
                Some(orientations) => {
                    for id in orientations.iter() {
                        out[*id] = self.evaluate_one_orientation((*id).into(), energy_cache);
                    }
                }
            },
        }

        out
    }

    pub fn evaluate_one_orientation(&self, orientation: TermId, energy_cache: &[F<f64>]) -> F<f64> {
        let expr = self.unwrap();
        if expr.orientations.is_empty() {
            panic!("no orientations to evaluate, most likely they are not generated, set compile_separate_orientations in config to generate them")
        }

        let mut out = [F(0.0)];
        expr.orientations[orientation]
            .borrow_mut()
            .evaluate(energy_cache, &mut out);
        out[0]
    }

    fn from_metedata(metadata: CompiledCFFExpressionMetaData) -> Result<Self, Report> {
        let path_to_joint = metadata
            .name
            .join(format!("expression_{}.so", metadata.graph_name));
        let path_to_joint_str = path_to_joint
            .to_str()
            .ok_or(eyre!("could not convert path to string"))?;

        if metadata.compile_cff_present {
            let joint = SerializableCompiledEvaluator::load(path_to_joint_str, "joint")
                .map(RefCell::new)
                .map_err(|e| eyre!(e))?;

            let orientations = if metadata.compile_separate_orientations_present {
                (0..metadata.num_orientations)
                    .map(|orientation| {
                        joint
                            .borrow()
                            .load_new_function(&format!("orientation_{}", orientation))
                            .map_err(|e| eyre!(e))
                            .map(RefCell::new)
                    })
                    .collect::<Result<_, Report>>()?
            } else {
                vec![].into()
            };

            let inner = InnerCompiledCFF {
                metadata,
                joint: Some(joint),
                orientations,
            };

            Ok(Self::Some(inner))
        } else if metadata.compile_separate_orientations_present {
            let orientation_zero =
                SerializableCompiledEvaluator::load(path_to_joint_str, "orientation_0")
                    .map(RefCell::new)
                    .map_err(|e| eyre!(e))?;

            let mut orientations = vec![orientation_zero];

            for orienatation_id in 1..metadata.num_orientations {
                let orientation_evaluator = orientations[0]
                    .borrow()
                    .load_new_function(&format!("orientation_{}", orienatation_id))
                    .map(RefCell::new)
                    .map_err(|e| eyre!(e))?;

                orientations.push(orientation_evaluator);
            }
            let inner = InnerCompiledCFF {
                metadata,
                joint: None,
                orientations: orientations.into(),
            };

            return Ok(Self::Some(inner));
        } else {
            return Ok(Self::None);
        }
    }

    fn unwrap(&self) -> &InnerCompiledCFF {
        match self {
            CompiledCFFExpression::Some(inner) => inner,
            CompiledCFFExpression::None => panic!("compiled cff not present"),
        }
    }

    fn is_some(&self) -> bool {
        match self {
            CompiledCFFExpression::Some(_) => true,
            CompiledCFFExpression::None => false,
        }
    }
}

impl Debug for CompiledCFFExpression {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Some(inner) => CompiledCFFExpressionMetaData::fmt(&inner.metadata, f),
            Self::None => str::fmt("None", f),
        }
    }
}
