use crate::{
    cff::{
        CutCFFIndex,
        expression::{
            OrientationID, ThreeDExpression, localize_three_d_expression_on_esurface,
            ltd_lu_local_series_coefficients_from_parametric_atom,
            prepare_ltd_lu_local_series_expression,
            remove_ltd_global_contact_completions_from_local_residue,
            select_lu_cut_residue_for_basis,
        },
    },
    debug_tags,
    graph::{
        FeynmanGraph, Graph, LoopMomentumBasis,
        cuts::{CutSet, LuResidueSelectionBasis},
    },
    momentum::Sign,
    numerator::symbolica_ext::AtomCoreExt,
    settings::global::{GenerationSettings, OrientationPattern, ThreeDRepresentation},
    utils::{GS, W_, symbolica_ext::LogPrint},
    uv::{
        ApproximationType, Spinney, UVgenerationSettings,
        approx::{
            expanded_4d::{
                Expanded4DApprox, Expanded4DSourceContext,
                expanded_4d_terms_to_3d_parametric_integrands, expanded_4d_uv_kernel,
                expanded_4d_uv_start,
            },
            final_integrand::FinalIntegrand,
            integrated::Integrated,
            local_3d::Local3DApproximation,
        },
    },
};
use color_eyre::Result;
use eyre::eyre;
use gammaloop_tracing_filter::{LogMessage, debug_instrument};
use idenso::{color::ColorSimplifier, shorthands::metric::MetricSimplifier};
use std::{collections::BTreeMap, hash::Hash};
use tracing::{debug, instrument};

use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, Symbol},
    function, symbol,
};
use three_dimensional_reps::RepresentationMode;

use linnet::half_edge::involution::{EdgeIndex, EdgeVec, Orientation};
use linnet::half_edge::subgraph::{InternalSubGraph, SuBitGraph, SubSetLike, SubSetOps};

use vakint::Vakint;
// use vakint::{EvaluationOrder, LoopNormalizationFactor, Vakint, VakintSettings};
use super::{IntegrandExpr, UltravioletGraph};

pub mod expanded_4d;
pub mod final_integrand;
pub mod integrated;
pub mod local_3d;

pub trait ForestNodeLike: LogMessage {
    fn subgraph(&self) -> &SuBitGraph;
    fn lmb(&self) -> &LoopMomentumBasis;
    fn lmb_id(&self) -> EdgeIndex {
        *self.lmb().loop_edges.first().unwrap()
    }
    // fn lmb_given(&self, subgraph: &SuBitGraph) -> &LoopMomentumBasis;
    fn dod(&self) -> i32;
    fn renormalization_scheme(&self) -> ApproximationType;
    fn topo_order(&self) -> usize;
    fn reduced_subgraph(&self, given: &Self) -> SuBitGraph;
}

pub trait ApproximationKernel<C> {
    fn kernel<S: ForestNodeLike>(
        &self,
        ctx: &C,
        current: &S,
        given: &S,
        atom: &Atom,
    ) -> Result<Atom>;
}

pub struct UVCtx<'a> {
    pub graph: &'a Graph,
    pub settings: &'a UVgenerationSettings,
}

pub trait ApproxKernel {
    fn apply<'a, A: Into<AtomOrView<'a>>>(&self, atom: A) -> Atom;
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum ApproxOp {
    NotComputed,
    Union {
        t_args: Vec<IntegrandExpr>,
        subgraphs: Vec<InternalSubGraph>,
    },
    Dependent {
        t_arg: IntegrandExpr,
        subgraph: InternalSubGraph,
    },
    Root,
}

#[derive(Clone)]
pub struct SimpleApprox {
    t_args: Vec<Atom>,
    pub sign: Sign,
    graph: InternalSubGraph,
}

impl SimpleApprox {
    fn subgraph_shadow(graph: &SuBitGraph, subgraph: &InternalSubGraph) -> Symbol {
        symbol!(&format!(
            "S_{}⊛{}",
            graph.string_label(),
            subgraph.string_label()
        ))
    }

    pub(crate) fn expr(&self, bigger_graph: &SuBitGraph) -> Atom {
        let reduced = Atom::var(Self::subgraph_shadow(bigger_graph, &self.graph));
        let mut mul = Atom::num(1);
        for i in &self.t_args {
            mul *= i;
        }
        reduced * mul
    }

    pub(crate) fn t_op(&self, bigger_graph: &SuBitGraph) -> Atom {
        function!(GS.top, self.expr(bigger_graph))
    }

    pub(crate) fn root(subgraph: InternalSubGraph) -> Self {
        if !subgraph.is_empty() {
            panic!(
                "Root approximation must be empty {} {:?}",
                subgraph.string_label(),
                subgraph
            )
        }
        SimpleApprox {
            t_args: vec![],
            sign: Sign::Positive,
            graph: subgraph,
        }
    }

    pub(crate) fn dependent(&self, bigger_graph: InternalSubGraph) -> Self {
        Self {
            t_args: vec![self.t_op(&bigger_graph.filter)],
            sign: -self.sign,
            graph: bigger_graph,
        }
    }
}

#[derive(Clone)]
pub struct Approximation {
    pub spinney: Spinney,
    pub local_3d: Local3DApprox,
    pub local_4d_expanded: Option<Expanded4DApprox>,
    pub final_integrand: Option<BTreeMap<CutCFFIndex, Atom>>,
    pub integrated_4d: ApproxOp,
    pub topo_order: usize,
    pub simple_approx: Option<SimpleApprox>,
    pub excluded_by_filters: bool,
}

impl LogMessage for Approximation {
    fn log_display(&self) -> String {
        format!(
            "subgraph={}, topo_order={}, dod={}",
            self.spinney.filter().string_label(),
            self.topo_order,
            self.spinney.dod
        )
    }
}

#[derive(Clone, Copy)]
pub(crate) struct RootResidueContext<'a> {
    pub(crate) expression: Option<&'a ThreeDExpression<OrientationID>>,
    pub(crate) lu_residue_selection_basis: LuResidueSelectionBasis,
    pub(crate) lu_residue_reference_basis: LuResidueSelectionBasis,
}

impl ForestNodeLike for Approximation {
    fn dod(&self) -> i32 {
        self.spinney.dod
    }

    fn renormalization_scheme(&self) -> ApproximationType {
        self.spinney.renormalization_scheme
    }

    fn lmb(&self) -> &LoopMomentumBasis {
        &self.spinney.lmb
    }

    // fn lmb_given(&self, subgraph: &SuBitGraph) -> &LoopMomentumBasis {

    //     self.
    // }

    fn reduced_subgraph(&self, given: &Self) -> SuBitGraph {
        self.spinney
            .subgraph
            .subtract(&given.spinney.subgraph)
            .filter
    }

    fn subgraph(&self) -> &SuBitGraph {
        self.spinney.filter()
    }

    fn topo_order(&self) -> usize {
        self.topo_order
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum Local3DApprox {
    NotComputed,
    Dependent { sign: Sign, t_arg: IntegrandExpr },
}

pub struct CutStructure {
    pub cuts: Vec<CutSet>,
}

impl CutStructure {
    pub(crate) fn empty(graph: &Graph) -> Self {
        Self {
            cuts: vec![CutSet::empty(graph.n_hedges())],
        }
    }
}

impl Local3DApprox {
    pub(crate) fn expr(&self) -> Option<(BTreeMap<CutCFFIndex, Atom>, Sign)> {
        match self {
            Local3DApprox::NotComputed => None,
            Local3DApprox::Dependent { sign, t_arg } => Some((t_arg.integrands.clone(), *sign)),
        }
    }

    pub(crate) fn dependent(
        graph: &mut Graph,
        to_contract: &SuBitGraph,
        cuts: &CutSet,
        _settings: &UVgenerationSettings,
        orientation_pattern: &OrientationPattern,
    ) -> Result<Local3DApprox> {
        let local_3d = graph
            .cff(
                &to_contract
                    .union(&graph.tree_edges)
                    .subtract(&graph.initial_state_cut),
                cuts,
                orientation_pattern,
            )?
            .expression_with_selectors();

        let fourddenoms = GS.wrap_tree_denoms(
            graph.denominator(&graph.tree_edges.subtract(&graph.initial_state_cut), |_| -1),
        );

        Ok(Local3DApprox::Dependent {
            sign: Sign::Positive,
            t_arg: IntegrandExpr {
                integrands: local_3d
                    .iter()
                    .map(|(index, atom)| (*index, atom * &fourddenoms))
                    .collect(),
            },
        })
    }

    pub(crate) fn root(
        graph: &mut Graph,
        cuts: &CutSet,
        settings: &UVgenerationSettings,
        orientation_pattern: &OrientationPattern,
    ) -> Result<Local3DApprox> {
        Self::dependent(
            graph,
            &graph.empty_subgraph::<SuBitGraph>(),
            cuts,
            settings,
            orientation_pattern,
        )
    }
}

#[derive(Clone, Copy)]
enum ResidueIndexAxis {
    RightThreshold,
    LeftThreshold,
    LuCut,
}

impl ResidueIndexAxis {
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
    axis: ResidueIndexAxis,
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

fn select_threshold_residue_for_representation(
    expression: ThreeDExpression<OrientationID>,
    threshold: &crate::cff::esurface::RaisedEsurfaceGroup,
    representation: ThreeDRepresentation,
) -> Vec<ThreeDExpression<OrientationID>> {
    match representation {
        ThreeDRepresentation::Cff => expression.select_esurface_residue(threshold),
        ThreeDRepresentation::Ltd => {
            expression.select_esurface_residue_in_generated_basis(threshold)
        }
    }
}

impl Approximation {
    fn zero_terms(allowed_keys: &[CutCFFIndex]) -> BTreeMap<CutCFFIndex, Atom> {
        allowed_keys.iter().map(|&key| (key, Atom::Zero)).collect()
    }

    fn zero_vec(n_terms: usize) -> Vec<Atom> {
        vec![Atom::Zero; n_terms]
    }

    fn indexed_terms_from_cutset(
        cutset: &CutSet,
        terms: Vec<Atom>,
        context: &str,
    ) -> Result<BTreeMap<CutCFFIndex, Atom>> {
        let allowed_keys = cutset.residue_selector.generate_allowed_keys();
        if terms.len() == allowed_keys.len() {
            return Ok(allowed_keys.into_iter().zip(terms).collect());
        }
        if terms.len() == 1 && terms[0].is_zero() {
            return Ok(Self::zero_terms(&allowed_keys));
        }
        Err(eyre!(
            "{context} generated {} local residue integrands, but the residue selector allows {} indices",
            terms.len(),
            allowed_keys.len()
        ))
    }

    pub(crate) fn set_filter_state(
        &mut self,
        graph: &Graph,
        dependent: &Self,
        settings: &UVgenerationSettings,
    ) {
        self.excluded_by_filters =
            settings
                .add_integrated_uv_with_n_loops
                .is_some_and(|target_loop_count| {
                    dependent.excluded_by_filters
                        || graph.n_loops(self.spinney.filter()) != target_loop_count
                });
    }

    pub(crate) fn root(
        &mut self,
        graph: &mut Graph,
        cuts: &CutSet,
        valid_orientations: &[EdgeVec<Orientation>],
        settings: &GenerationSettings,
        root_residue_context: RootResidueContext<'_>,
    ) -> Result<()> {
        self.excluded_by_filters = false;
        self.simple_approx = Some(SimpleApprox::root(self.spinney.subgraph.clone()));
        if settings.uv.only_integrated {
            self.integrated_4d = ApproxOp::Root;
        } else {
            self.integrated_4d = ApproxOp::Root;
            self.local_3d =
                Local3DApprox::root(graph, cuts, &settings.uv, &settings.orientation_pattern)?;
            if settings.explicit_orientation_sum_only {
                let root_expression = root_residue_context.expression.ok_or_else(|| {
                    eyre!(
                        "explicit orientation-summed local-3D UV generation requires the production root 3D expression"
                    )
                })?;
                self.final_integrand = Some(self.final_integrand_from_root_expression(
                    graph,
                    cuts,
                    root_expression,
                    valid_orientations,
                    settings,
                    ThreeDRepresentation::Cff,
                    root_residue_context.lu_residue_selection_basis,
                    root_residue_context.lu_residue_reference_basis,
                    true,
                )?);
            } else {
                let (local_terms, local_sign) = self
                    .local_3d
                    .expr()
                    .expect("root local CFF should have been computed");
                let uv_marker = (settings.uv.add_sigma && settings.uv.keep_sigma).then(|| {
                    self.simple_approx
                        .as_ref()
                        .unwrap()
                        .expr(&graph.full_filter())
                });
                self.final_integrand = Some(
                    FinalIntegrand::new(
                        valid_orientations,
                        &settings.orientation_pattern,
                        uv_marker.as_ref(),
                        true,
                    )
                    .build(
                        graph,
                        self,
                        &local_terms,
                        local_sign,
                        &self.integrated_4d,
                        cuts,
                    )?,
                );
            }
        }

        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn final_integrand_from_root_expression(
        &self,
        graph: &Graph,
        cutset: &CutSet,
        expression: &ThreeDExpression<OrientationID>,
        valid_orientations: &[EdgeVec<Orientation>],
        settings: &GenerationSettings,
        representation: ThreeDRepresentation,
        lu_residue_selection_basis: LuResidueSelectionBasis,
        lu_residue_reference_basis: LuResidueSelectionBasis,
        average_for_outer_orientation_projection: bool,
    ) -> Result<BTreeMap<CutCFFIndex, Atom>> {
        let numerator = graph.production_numerator_atom_for_full_3d_expression();
        // This root is projected directly from the generated 3D expression and
        // therefore bypasses `Graph::cff`, where the loop-measure factor is
        // normally attached.
        let loop_number = graph.get_loop_number();
        let normalization = (-Atom::i()).pow(loop_number as i64)
            / (Atom::var(GS.pi) * 2).pow(3 * loop_number as i64);

        let mut residues = vec![(CutCFFIndex::new_all_none(), expression.clone())];
        if let Some(right_threshold) = cutset.residue_selector.right_th_cut.as_ref() {
            residues = apply_indexed_residue_selection(
                residues,
                ResidueIndexAxis::RightThreshold,
                |expression| {
                    select_threshold_residue_for_representation(
                        expression,
                        right_threshold,
                        representation,
                    )
                },
            );
        }
        if let Some(left_threshold) = cutset.residue_selector.left_th_cut.as_ref() {
            residues = apply_indexed_residue_selection(
                residues,
                ResidueIndexAxis::LeftThreshold,
                |expression| {
                    select_threshold_residue_for_representation(
                        expression,
                        left_threshold,
                        representation,
                    )
                },
            );
        }
        if representation == ThreeDRepresentation::Ltd
            && cutset.residue_selector.has_lu_cut_residue()
            && let Some(lu_cut) = cutset.residue_selector.lu_cut()
        {
            let residue_prefactor_sign = cutset
                .residue_selector
                .ltd_direct_original_residue_prefactor_sign();
            let mut atoms = BTreeMap::new();
            for (index, mut residue) in residues {
                prepare_ltd_lu_local_series_expression(&mut residue, cutset);
                self.localize_ltd_threshold_residue_if_needed(&mut residue, cutset)?;
                let atom = graph.three_d_expression_parametric_atom_with_numerator_gs(
                    &residue,
                    &numerator,
                    RepresentationMode::Ltd,
                    true,
                    &settings.orientation_pattern,
                ) * &normalization;
                for (i, coefficient) in ltd_lu_local_series_coefficients_from_parametric_atom(
                    &graph.name,
                    &residue,
                    atom,
                    lu_cut,
                    cutset
                        .residue_selector
                        .ltd_lu_cut_local_series_esurface_signs(),
                    cutset
                        .residue_selector
                        .ltd_lu_cut_local_series_coordinates(),
                    false,
                )?
                .into_iter()
                .enumerate()
                {
                    let mut coefficient_index = index;
                    coefficient_index.lu_cut_order = Some(i + 1);
                    atoms.insert(
                        coefficient_index,
                        coefficient * Atom::num(residue_prefactor_sign),
                    );
                }
            }
            return atoms
                .into_iter()
                .map(|(index, mut atom)| {
                    if average_for_outer_orientation_projection {
                        // The root term comes from the explicitly summed
                        // production expression, while CFF local-3D forests
                        // apply the orientation projector once to the whole
                        // forest after all local terms are assembled. Embed
                        // the root as a uniform orientation average so that
                        // this later projector acts as the identity on it.
                        if valid_orientations.is_empty() {
                            return Err(eyre!(
                                "cannot embed explicitly summed root 3D expression into an empty orientation projector"
                            ));
                        }
                        atom /= Atom::num(valid_orientations.len() as i64);
                    }
                    Ok((
                        index,
                        atom.replace(GS.dim)
                            .with(4)
                            .simplify_color()
                            .replace(GS.den(W_.a_, W_.b_, W_.c_, W_.d_))
                            .with(W_.d_)
                            .expand_dots()?
                            .collect_factors(),
                    ))
                })
                .collect();
        }
        if let Some(lu_cut) = cutset.residue_selector.lu_cut() {
            residues =
                apply_indexed_residue_selection(residues, ResidueIndexAxis::LuCut, |expression| {
                    select_lu_cut_residue_for_basis(
                        expression,
                        lu_cut,
                        cutset.residue_selector.lu_cut_edge_sets(),
                        cutset.residue_selector.ltd_lu_cut_esurface_signs(),
                        lu_residue_selection_basis,
                    )
                });
        }
        if representation == ThreeDRepresentation::Ltd {
            for (_, residue) in &mut residues {
                if cutset.residue_selector.lu_cut().is_some() {
                    remove_ltd_global_contact_completions_from_local_residue(residue);
                }
                self.localize_ltd_threshold_residue_if_needed(residue, cutset)?;
            }
        }
        residues
            .into_iter()
            .map(|(index, residue)| {
                let mut atom = graph.three_d_expression_parametric_atom_with_numerator_gs(
                    &residue,
                    &numerator,
                    match representation {
                        ThreeDRepresentation::Cff => RepresentationMode::Cff,
                        ThreeDRepresentation::Ltd => RepresentationMode::Ltd,
                    },
                    true,
                    &settings.orientation_pattern,
                ) * &normalization;
                if representation == ThreeDRepresentation::Ltd {
                    let residue_prefactor_sign = match lu_residue_reference_basis {
                        LuResidueSelectionBasis::GeneratedEsurface => cutset
                            .residue_selector
                            .ltd_local_series_residue_prefactor_sign(),
                        LuResidueSelectionBasis::PositiveEnergyCutkosky => {
                            cutset.residue_selector.ltd_residue_prefactor_sign()
                        }
                    };
                    atom *= Atom::num(residue_prefactor_sign);
                }
                if average_for_outer_orientation_projection {
                    // The root term comes from the explicitly summed production
                    // expression, while CFF local-3D forests apply the orientation
                    // projector once to the whole forest after all local terms are
                    // assembled. Embed the root as a uniform orientation average
                    // so that this later projector acts as the identity on it.
                    if valid_orientations.is_empty() {
                        return Err(eyre!(
                            "cannot embed explicitly summed root 3D expression into an empty orientation projector"
                        ));
                    }
                    atom /= Atom::num(valid_orientations.len() as i64);
                }
                let atom = atom
                    .replace(GS.dim)
                    .with(4)
                    .simplify_color()
                    .replace(GS.den(W_.a_, W_.b_, W_.c_, W_.d_))
                    .with(W_.d_)
                    .expand_dots()?
                    .collect_factors();
                Ok((index, atom))
            })
            .collect()
    }

    fn localize_ltd_threshold_residue_if_needed(
        &self,
        residue: &mut ThreeDExpression<OrientationID>,
        cutset: &CutSet,
    ) -> Result<()> {
        if let Some(right_threshold) = cutset.residue_selector.right_th_cut.as_ref() {
            for esurface_id in &right_threshold.esurface_ids {
                localize_three_d_expression_on_esurface(residue, *esurface_id)?;
            }
        }
        if let Some(left_threshold) = cutset.residue_selector.left_th_cut.as_ref() {
            for esurface_id in &left_threshold.esurface_ids {
                localize_three_d_expression_on_esurface(residue, *esurface_id)?;
            }
        }
        Ok(())
    }

    pub(crate) fn root_expanded_4d(
        &mut self,
        graph: &mut Graph,
        cuts: &CutSet,
        valid_orientations: &[EdgeVec<Orientation>],
        settings: &GenerationSettings,
        representation: ThreeDRepresentation,
        root_residue_context: RootResidueContext<'_>,
    ) -> Result<()> {
        self.excluded_by_filters = false;
        self.simple_approx = Some(SimpleApprox::root(self.spinney.subgraph.clone()));
        if settings.uv.only_integrated {
            self.integrated_4d = ApproxOp::Root;
        } else {
            self.integrated_4d = ApproxOp::Root;
            self.local_4d_expanded = Some(Expanded4DApprox::root(graph));
            let root_expression = root_residue_context.expression.ok_or_else(|| {
                eyre!(
                    "expanded-4D local UV generation requires the production root 3D expression for the original integrand"
                )
            })?;
            self.final_integrand = Some(self.final_integrand_from_root_expression(
                graph,
                cuts,
                root_expression,
                valid_orientations,
                settings,
                representation,
                root_residue_context.lu_residue_selection_basis,
                root_residue_context.lu_residue_reference_basis,
                representation == ThreeDRepresentation::Cff
                    && !settings.explicit_orientation_sum_only,
            )?);
        }

        Ok(())
    }

    pub(crate) fn new(spinney: Spinney) -> Approximation {
        Approximation {
            spinney,
            topo_order: 0,
            final_integrand: None,
            simple_approx: None,
            local_3d: Local3DApprox::NotComputed,
            local_4d_expanded: None,
            integrated_4d: ApproxOp::NotComputed,
            excluded_by_filters: false,
        }
    }

    #[debug_instrument(
        graph = %graph.log_display(),
        current = %self.log_display(),
        given = %dependent.log_display(),
        reduced = ?self.reduced_subgraph(dependent),
    )]
    pub(crate) fn compute_integrated(
        &mut self,
        graph: &Graph,
        vakint: (&Vakint, &vakint::VakintSettings),
        dependent: &Self,
        settings: &UVgenerationSettings,
    ) -> Result<()> {
        if self.excluded_by_filters {
            debug_tags!(#generation, #profile, #uv, #graph, #summary;
                stage = "compute_integrated_filtered_zero",

                "Skipped integrated UV computation for filtered mode"
            );
            self.integrated_4d = ApproxOp::Dependent {
                t_arg: IntegrandExpr {
                    integrands: BTreeMap::from([(CutCFFIndex::new_all_none(), Atom::Zero)]),
                },
                subgraph: unsafe {
                    InternalSubGraph::new_unchecked(self.reduced_subgraph(dependent))
                },
            };
            return Ok(());
        }

        let ctx = UVCtx { graph, settings };

        let Some(current) = &dependent.integrated_4d.expr() else {
            return Err(eyre!("integrated_4d not computed"));
        };

        debug!(pole_part = %settings.pole_part,
            simple = % self.simple_approx
                .as_ref()
                .unwrap()
                .expr(&graph.full_filter()).log_print(None),
            "Computing Integrated",
        );

        let started = std::time::Instant::now();
        debug_tags!(#generation, #profile, #uv, #graph, #summary;
            stage = "compute_integrated_start",
            integrand_count = current.len(),
            "Computing integrated UV CT"
        );
        let integrands = current
            .iter()
            .map(|(index, a)| {
                let term_started = std::time::Instant::now();
                let integrated =
                    Integrated::new(vakint.0, vakint.1).kernel(&ctx, self, dependent, a)?;
                let integrated = if dependent.subgraph().is_empty() {
                    integrated
                } else {
                    -integrated
                };
                debug_tags!(#generation, #profile, #uv, #graph, #term, #summary;
                    stage = "compute_integrated_kernel_done",
                    cut_index = ?index,
                    elapsed_ms = term_started.elapsed().as_secs_f64() * 1000.0,
                    "Computed integrated UV CT term"
                );
                Ok((*index, integrated))
            })
            .collect::<Result<BTreeMap<_, _>>>()?;

        self.integrated_4d = ApproxOp::Dependent {
            t_arg: IntegrandExpr { integrands },
            subgraph: unsafe { InternalSubGraph::new_unchecked(self.reduced_subgraph(dependent)) },
        };

        debug_tags!(#generation, #profile, #uv, #graph, #summary;
            stage = "compute_integrated_done",
            integrand_count = current.len(),
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Computed integrated UV CT"
        );
        Ok(())
    }

    /// Computes the 3d approximation of the UV
    #[debug_instrument(
        graph = %graph.log_display(),
        current = %self.log_display(),
        given = %dependent.log_display(),
    )]
    pub(crate) fn compute(
        &mut self,
        graph: &mut Graph,
        cuts: &CutSet,
        dependent: &Self,
        valid_orientations: &[EdgeVec<Orientation>],
        settings: &GenerationSettings,
    ) -> Result<()> {
        let started = std::time::Instant::now();
        debug_tags!(#generation, #profile, #uv, #graph, #summary;
            stage = "compute_local_3d_start",
            "Computing local 3D UV CT"
        );
        let Some((cff, sign)) = dependent.local_3d.expr() else {
            return Err(eyre!("dependent local 3D approximation was not computed"));
        };
        debug_tags!(#generation, #profile, #uv, #graph, #summary;
            stage = "compute_local_3d_parent_expr_done",
            cff_count = cff.len(),
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Loaded parent local 3D UV CT"
        );
        let final_integrand = FinalIntegrand::new(
            valid_orientations,
            &settings.orientation_pattern,
            None,
            false,
        );

        let localize_started = std::time::Instant::now();
        let integrated_t = final_integrand.localized_integrated_ct(
            graph,
            dependent,
            &dependent.integrated_4d,
            cuts,
        )?;
        debug_tags!(#generation, #profile, #uv, #graph, #summary;
            stage = "compute_local_3d_localize_integrated_done",
            localized_count = integrated_t.active.integrands.len(),
            elapsed_ms = localize_started.elapsed().as_secs_f64() * 1000.0,
            total_elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Localized integrated UV CT for local 3D subtraction"
        );

        let ctx = UVCtx {
            graph,
            settings: &settings.uv,
        };
        let uv_marker = (settings.uv.add_sigma && settings.uv.keep_sigma).then(|| {
            self.simple_approx
                .as_ref()
                .unwrap()
                .expr(&ctx.graph.full_filter())
        });

        if cff.len() != integrated_t.active.integrands.len() {
            return Err(eyre!(
                "local 3D UV approximation produced {} residue integrands, but localized finite integrated UV produced {}",
                cff.len(),
                integrated_t.active.integrands.len()
            ));
        }

        let mut integrands = BTreeMap::new();
        let mut tagged_integrands = uv_marker.as_ref().map(|_| BTreeMap::new());
        for (((index_local, local), (index_integ, integ)), (index_frozen, frozen)) in cff
            .into_iter()
            .zip(integrated_t.active.integrands)
            .zip(integrated_t.frozen_integrands)
        {
            let term_started = std::time::Instant::now();
            if index_local != index_integ {
                return Err(eyre!(
                    "Mismatched indices for local and integrated approximations: {:?} vs {:?}",
                    index_local,
                    index_integ
                ));
            }
            if index_local != index_frozen {
                return Err(eyre!(
                    "Mismatched indices for local integrated active and frozen factors: {:?} vs {:?}",
                    index_local,
                    index_frozen
                ));
            }
            let mut sum_3d = Atom::Zero;

            debug_tags!(#generation, #profile, #uv, #integrated, #local, #graph, #term, #trace;
                stage = "compute_local_3d_term_inputs",
                cut_index = ?index_local,
                log.local = local,
                log.localized_integrated = integ,
                log.frozen = frozen,
                "Local 3D UV CT term inputs"
            );

            debug_tags!(#generation, #profile, #uv, #graph, #term, #summary;
                "adding localT(expr)"
            );
            let local_ct = Local3DApproximation::full().kernel(&ctx, &*self, dependent, &local)?;
            debug_tags!(#generation, #profile, #uv, #local, #graph, #term, #trace;
                stage = "compute_local_3d_term_local_done",
                cut_index = ?index_local,
                log.local_ct = local_ct,
                "Computed localT(expr)"
            );
            debug_tags!(#generation, #profile, #uv, #graph, #term, #summary;
                "subtracting localT(integrated(expr))"
            );
            let localized_integrated_ct =
                Local3DApproximation::reduced().kernel(&ctx, &*self, dependent, &integ)? * frozen;
            debug_tags!(#generation, #profile, #uv, #integrated, #local, #graph, #term, #trace;
                stage = "compute_local_3d_term_integrated_done",
                cut_index = ?index_local,
                log.localized_integrated_ct = localized_integrated_ct,
                "Computed localT(integrated(expr))"
            );
            if let (Some(marker), Some(tagged_integrands)) =
                (uv_marker.as_ref(), tagged_integrands.as_mut())
            {
                let tagged_sum_3d = local_ct.clone() * function!(GS.uv_local, marker.clone())
                    - localized_integrated_ct.clone() * function!(GS.uv_integrated, marker.clone());
                tagged_integrands.insert(index_local, tagged_sum_3d.clone());
                debug_tags!(#generation, #profile, #uv, #local, #graph, #term, #trace;
                    stage = "compute_local_3d_term_tagged_output",
                    cut_index = ?index_local,
                    log.tagged_local_3d_ct = tagged_sum_3d,
                    "Computed tagged local 3D UV CT expression"
                );
            }
            sum_3d += local_ct;
            sum_3d -= localized_integrated_ct;
            debug_tags!(#generation, #profile, #uv, #local, #graph, #term, #trace;
                stage = "compute_local_3d_term_output",
                cut_index = ?index_local,
                log.local_3d_ct = sum_3d,
                "Computed local 3D UV CT expression"
            );
            integrands.insert(index_local, sum_3d);
            debug_tags!(#generation, #profile, #uv, #graph, #term, #summary;
                stage = "compute_local_3d_term_done",
                cut_index = ?index_local,
                elapsed_ms = term_started.elapsed().as_secs_f64() * 1000.0,
                total_elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
                "Computed local 3D UV CT term"
            );
        }

        let output_sign = -sign;
        self.local_3d = Local3DApprox::Dependent {
            sign: output_sign,
            t_arg: IntegrandExpr { integrands },
        };

        let final_started = std::time::Instant::now();
        self.final_integrand = Some({
            let (local_terms, local_sign) = self
                .local_3d
                .expr()
                .ok_or(eyre!("Local3d not yet computed"))?;
            let final_local_terms = match tagged_integrands.as_ref() {
                Some(tagged_integrands) => tagged_integrands,
                None => &local_terms,
            };
            FinalIntegrand::new(
                valid_orientations,
                &settings.orientation_pattern,
                uv_marker.as_ref(),
                false,
            )
            .build(
                graph,
                self,
                final_local_terms,
                local_sign,
                &self.integrated_4d,
                cuts,
            )?
        });
        debug_tags!(#generation, #profile, #uv, #graph, #summary;
            stage = "compute_local_3d_done",
            parent_sign = ?sign,
            output_sign = ?output_sign,
            final_integrand_elapsed_ms = final_started.elapsed().as_secs_f64() * 1000.0,
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Computed local 3D UV CT"
        );
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn compute_expanded_4d(
        &mut self,
        graph: &mut Graph,
        cuts: &CutSet,
        dependent: &Self,
        valid_orientations: &[EdgeVec<Orientation>],
        settings: &GenerationSettings,
        representation: ThreeDRepresentation,
        root_expression: Option<&ThreeDExpression<OrientationID>>,
    ) -> Result<()> {
        let Some(parent) = dependent.local_4d_expanded.as_ref() else {
            return Err(eyre!("expanded 4D local UV parent was not computed"));
        };
        let (parent_terms, parent_sign) = parent.expr();

        let terms = parent_terms
            .iter()
            .map(|parent_term| {
                let start = expanded_4d_uv_start(graph, self, dependent, parent_term)?;
                expanded_4d_uv_kernel(graph, self, dependent, &start)
            })
            .collect::<Result<Vec<_>>>()?;

        self.local_4d_expanded = Some(Expanded4DApprox {
            sign: -parent_sign,
            terms,
        });

        self.final_integrand = Some(self.final_integrand_expanded_4d(
            graph,
            cuts,
            valid_orientations,
            settings,
            representation,
            root_expression,
        )?);
        Ok(())
    }

    #[instrument(skip_all)]
    pub(crate) fn final_integrand_expanded_4d(
        &self,
        graph: &mut Graph,
        cutset: &CutSet,
        valid_orientations: &[EdgeVec<Orientation>],
        settings: &crate::settings::global::GenerationSettings,
        representation: crate::settings::global::ThreeDRepresentation,
        root_expression: Option<&ThreeDExpression<OrientationID>>,
    ) -> Result<BTreeMap<CutCFFIndex, Atom>> {
        let Some(local_4d) = self.local_4d_expanded.as_ref() else {
            return Err(eyre!("expanded 4D local UV term not yet computed"));
        };
        let (local_terms, sign) = local_4d.expr();
        let finite = finite_integrated_part(&self.integrated_4d)?;

        let reduced = graph
            .full_filter()
            .subtract(self.spinney.subgraph.included())
            .subtract(&graph.initial_state_cut);
        let outside_numerator = graph
            .numerator(&reduced, self.spinney.subgraph.included())
            .get_single_atom()
            .map_err(|error| eyre!("graph numerator is not a single symbolic atom: {error}"))?
            * graph.global_atom();

        let mut local_atom = Atom::Zero;
        for term in local_terms {
            local_atom += sign * term.clone() * &outside_numerator;
        }

        let mut integrands = expanded_4d_terms_to_3d_parametric_integrands(
            graph,
            &local_atom,
            cutset,
            representation,
            settings,
            valid_orientations,
            Expanded4DSourceContext::local_uv(self.spinney.filter(), self.spinney.filter()),
            root_expression,
        )?;

        if finite.is_zero() {
            return Self::indexed_terms_from_cutset(cutset, integrands, "expanded 4D local UV");
        }

        let finite_atom = {
            let reduced_denominator = Atom::num(1) / graph.denominator(&reduced, |_| 1);
            -sign * finite * reduced_denominator * &outside_numerator
        };

        let mut finite_integrands = expanded_4d_terms_to_3d_parametric_integrands(
            graph,
            &finite_atom,
            cutset,
            representation,
            settings,
            valid_orientations,
            Expanded4DSourceContext::cograph_only(self.spinney.filter()),
            root_expression,
        )?;
        if integrands.len() != finite_integrands.len() {
            if integrands.len() == 1 && integrands[0].is_zero() {
                integrands = Self::zero_vec(finite_integrands.len());
            } else if finite_integrands.len() == 1 && finite_integrands[0].is_zero() {
                finite_integrands = Self::zero_vec(integrands.len());
            } else {
                return Err(eyre!(
                    "expanded 4D local UV generated {} local residue integrands, but finite integrated UV generated {}",
                    integrands.len(),
                    finite_integrands.len()
                ));
            }
        }
        for (sum, finite) in integrands.iter_mut().zip(finite_integrands) {
            *sum += finite;
        }
        Self::indexed_terms_from_cutset(cutset, integrands, "expanded 4D local UV")
    }

    // pub(crate) fn simple_expr(
    //     &self,
    //     graph: &UVGraph,
    //     amplitude: &InternalSubGraph,
    // ) -> Option<Atom> {
    //     let simple_approx = self.simple_approx.as_ref()?;

    //     Some((simple_approx.sign * simple_approx.expr(&amplitude.filter)).into())
    // }
}

fn finite_integrated_part(integrated_4d: &ApproxOp) -> Result<Atom> {
    let mut t4 = match integrated_4d {
        ApproxOp::Root | ApproxOp::NotComputed => return Ok(Atom::Zero),
        ApproxOp::Union { .. } => {
            return Err(eyre!(
                "cannot extract finite integrated UV part from an unresolved union approximation"
            ));
        }
        ApproxOp::Dependent { t_arg, .. } => t_arg.integrands.clone(),
    };

    if t4.len() != 1 {
        return Err(eyre!("expected one integrated 4D approximation term"));
    }
    let finite_term = t4
        .remove(&CutCFFIndex::new_all_none())
        .ok_or_else(|| eyre!("expected one integrated 4D approximation term"))?;
    Ok(finite_term
        .series(GS.dim_epsilon, Atom::Zero, 0)
        .map_err(|error| eyre!("finite integrated UV epsilon expansion failed: {error}"))?
        .coefficient((0, 1).into())
        .replace(GS.m_uv_expansion)
        .with(GS.m_uv_vacuum)
        .map_mink_dim(4)
        .replace(function!(symbol!("vakint::g"), W_.a__))
        .with(function!(symbol!("spenso::g"), W_.a__)))
}

impl ApproxOp {
    pub(crate) fn expr(&self) -> Option<BTreeMap<CutCFFIndex, Atom>> {
        match self {
            ApproxOp::NotComputed => None,
            ApproxOp::Union { .. } => {
                //Never gets hit now
                panic!(
                    "Should not call expr on a union approximation, need to compute the union first"
                );
            }
            ApproxOp::Dependent { t_arg, .. } => Some(t_arg.integrands.clone()),
            ApproxOp::Root => Some(BTreeMap::from([(
                CutCFFIndex::new_all_none(),
                Atom::num(1),
            )])),
        }
    }
}
