use crate::{
    cff::expression::{
        CFFExpression, OrientationID, localize_numerator_energy_maps_on_esurface,
        remove_ltd_global_contact_completions_from_local_residue,
    },
    graph::{Graph, LMBext, LoopMomentumBasis, cuts::CutSet},
    momentum::Sign,
    numerator::symbolica_ext::AtomCoreExt,
    settings::global::{GenerationSettings, ThreeDRepresentation},
    utils::{
        GS, W_,
        symbolica_ext::{LOGPRINTOPTS, LogPrint},
    },
    uv::{
        ApproximationType, Spinney, UVgenerationSettings,
        approx::{
            expanded_4d::{
                Expanded4DApprox, expanded_4d_terms_to_3d_parametric_integrands,
                expanded_4d_uv_kernel, expanded_4d_uv_start,
            },
            integrated::Integrated,
            local_3d::Local3DApproximation,
        },
    },
};
use color_eyre::Result;
use eyre::eyre;
use idenso::{color::ColorSimplifier, metric::MetricSimplifier};
use std::hash::Hash;
use tracing::debug;

use spenso::network::library::TensorLibraryData;
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, Symbol},
    function, parse_lit, symbol,
};
use three_dimensional_reps::RepresentationMode;

use linnet::half_edge::involution::{EdgeIndex, EdgeVec, HedgePair, Orientation};
use linnet::half_edge::subgraph::{InternalSubGraph, SuBitGraph, SubSetLike, SubSetOps};

use tracing::instrument;
use vakint::{Vakint, vakint_symbol};
// use vakint::{EvaluationOrder, LoopNormalizationFactor, Vakint, VakintSettings};
use super::{IntegrandExpr, UltravioletGraph};

pub mod expanded_4d;
pub mod integrated;
pub mod local_3d;
pub mod orientation_localization;

use orientation_localization::localize_reduced_orientation_term;

pub trait ForestNodeLike {
    fn subgraph(&self) -> &SuBitGraph;
    fn lmb(&self) -> &LoopMomentumBasis;
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
        sign: Sign,
        t_args: Vec<IntegrandExpr>,
        subgraphs: Vec<InternalSubGraph>,
    },
    Dependent {
        sign: Sign,
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
    pub local_3d: CFFapprox, //3d denoms
    pub local_4d_expanded: Option<Expanded4DApprox>,
    pub final_integrand: Option<Vec<Atom>>,
    pub integrated_4d: ApproxOp,
    pub topo_order: usize,
    pub simple_approx: Option<SimpleApprox>,
    pub generate_only_integrated_uv_chain_matches: bool,
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
pub enum CFFapprox {
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

impl CFFapprox {
    pub(crate) fn expr(&self) -> Option<(Vec<Atom>, Sign)> {
        match self {
            CFFapprox::NotComputed => None,
            CFFapprox::Dependent { sign, t_arg } => Some((t_arg.integrands.clone(), *sign)),
        }
    }

    pub(crate) fn dependent(
        graph: &mut Graph,
        to_contract: &SuBitGraph,
        cuts: &CutSet,
        _settings: &UVgenerationSettings,
    ) -> Result<CFFapprox> {
        let cff = graph.cff(to_contract, cuts)?.expression_with_selectors();

        Ok(CFFapprox::Dependent {
            sign: Sign::Positive,
            t_arg: IntegrandExpr { integrands: cff },
        })
    }

    pub(crate) fn root(
        graph: &mut Graph,
        cuts: &CutSet,
        settings: &UVgenerationSettings,
    ) -> Result<CFFapprox> {
        Self::dependent(graph, &graph.empty_subgraph::<SuBitGraph>(), cuts, settings)
    }
}

fn internal_paired_edges_of_subgraph(graph: &Graph, subgraph: &SuBitGraph) -> Vec<EdgeIndex> {
    let mut edge_ids: Vec<_> = graph
        .as_ref()
        .iter_edges_of(subgraph)
        .filter_map(|(pair, edge_id, _)| pair.is_paired().then_some(edge_id))
        .collect();
    edge_ids.sort_by_key(|edge| usize::from(*edge));
    edge_ids
}

fn localized_integrated_reduced_factor(
    graph: &mut Graph,
    to_contract: &SuBitGraph,
    cuts: &CutSet,
    valid_orientations: &[EdgeVec<Orientation>],
    explicit_orientation_sum_only: bool,
) -> Result<IntegrandExpr> {
    let cff = graph.cff(to_contract, cuts)?;

    if explicit_orientation_sum_only {
        if valid_orientations.is_empty() {
            return Err(eyre!(
                "cannot embed explicitly summed reduced integrated-UV factor into an empty orientation projector"
            ));
        }
        let normalization = Atom::num(valid_orientations.len() as i64);
        let integrands = cff
            .terms
            .into_iter()
            .map(|term| {
                term.expression
                    .into_iter()
                    .fold(Atom::Zero, |acc, expr| acc + expr)
                    / normalization.clone()
            })
            .collect();

        return Ok(IntegrandExpr { integrands });
    }

    let internal_edges = internal_paired_edges_of_subgraph(graph, to_contract);

    let integrands = cff
        .terms
        .into_iter()
        .map(|term| {
            let mut localized = Atom::Zero;
            for (expr, orientation) in term.expression.into_iter().zip(term.orientations) {
                localized += localize_reduced_orientation_term(
                    &expr,
                    &orientation,
                    valid_orientations,
                    &internal_edges,
                )?;
            }
            Ok(localized)
        })
        .collect::<Result<Vec<_>>>()?;

    Ok(IntegrandExpr { integrands })
}

impl Approximation {
    fn filtered_integrated_uv_mode_is_active(settings: &UVgenerationSettings) -> bool {
        settings.generate_integrated && settings.filtered_integrated_uv_loop_count().is_some()
    }

    fn initialize_filtered_integrated_uv_root(&mut self, settings: &UVgenerationSettings) {
        self.generate_only_integrated_uv_chain_matches =
            Self::filtered_integrated_uv_mode_is_active(settings);
    }

    pub(crate) fn update_filtered_integrated_uv_chain_state(
        &mut self,
        graph: &Graph,
        dependent: &Self,
        settings: &UVgenerationSettings,
    ) {
        self.generate_only_integrated_uv_chain_matches = settings
            .filtered_integrated_uv_loop_count()
            .is_some_and(|target_loop_count| {
                dependent.generate_only_integrated_uv_chain_matches
                    && graph.n_loops(self.spinney.filter()) == target_loop_count
            });
    }

    fn zero_terms(n_terms: usize) -> Vec<Atom> {
        vec![Atom::Zero; n_terms]
    }

    fn set_zero_local_3d(&mut self, sign: Sign, n_terms: usize) {
        self.local_3d = CFFapprox::Dependent {
            sign,
            t_arg: IntegrandExpr {
                integrands: Self::zero_terms(n_terms),
            },
        };
    }

    fn set_zero_local_4d_expanded(&mut self, sign: Sign, n_terms: usize) {
        self.local_4d_expanded = Some(Expanded4DApprox {
            sign,
            terms: Self::zero_terms(n_terms),
        });
    }

    fn should_keep_only_integrated_uv_terms(&self, settings: &UVgenerationSettings) -> bool {
        Self::filtered_integrated_uv_mode_is_active(settings)
            && self.generate_only_integrated_uv_chain_matches
    }

    pub(crate) fn root(
        &mut self,
        graph: &mut Graph,
        cuts: &CutSet,
        valid_orientations: &[EdgeVec<Orientation>],
        settings: &GenerationSettings,
        root_expression: Option<&CFFExpression<OrientationID>>,
    ) -> Result<()> {
        self.initialize_filtered_integrated_uv_root(&settings.uv);
        self.simple_approx = Some(SimpleApprox::root(self.spinney.subgraph.clone()));
        if settings.uv.only_integrated {
            self.integrated_4d = ApproxOp::Root;
        } else {
            self.integrated_4d = ApproxOp::Root;
            self.local_3d = CFFapprox::root(graph, cuts, &settings.uv)?;
            if Self::filtered_integrated_uv_mode_is_active(&settings.uv) {
                let (integrands, _) = self
                    .local_3d
                    .expr()
                    .expect("root local CFF should have been computed");
                self.final_integrand = Some(Self::zero_terms(integrands.len()));
            } else if settings.explicit_orientation_sum_only {
                let root_expression = root_expression.ok_or_else(|| {
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
                    true,
                )?);
            } else {
                self.final_integrand = Some(self.final_integrand(
                    graph,
                    cuts,
                    valid_orientations,
                    &settings.uv,
                    false,
                )?);
            }
        }

        Ok(())
    }

    fn final_integrand_from_root_expression(
        &self,
        graph: &Graph,
        cutset: &CutSet,
        expression: &CFFExpression<OrientationID>,
        valid_orientations: &[EdgeVec<Orientation>],
        settings: &GenerationSettings,
        representation: ThreeDRepresentation,
        average_for_outer_orientation_projection: bool,
    ) -> Result<Vec<Atom>> {
        let reduced = graph.full_filter().subtract(&graph.initial_state_cut);
        let numerator = graph
            .numerator(&reduced, &graph.empty_subgraph())
            .get_single_atom()
            .map_err(|error| eyre!("graph numerator is not a single symbolic atom: {error}"))?
            * graph.global_atom();
        let momentum_replacements = graph.normal_emr_replacement(
            &graph.full_filter(),
            &graph.loop_momentum_basis,
            &[W_.x___],
            HedgePair::is_paired,
        );
        let numerator = numerator
            .replace_multiple(&momentum_replacements)
            .expand()
            .collect_factors();

        let mut residues = vec![expression.clone()];
        if let Some(right_threshold) = cutset.residue_selector.right_th_cut.as_ref() {
            residues = residues
                .into_iter()
                .flat_map(|expression| expression.select_esurface_residue(right_threshold))
                .collect();
        }
        if let Some(left_threshold) = cutset.residue_selector.left_th_cut.as_ref() {
            residues = residues
                .into_iter()
                .flat_map(|expression| expression.select_esurface_residue(left_threshold))
                .collect();
        }
        let lu_cut_representative_esurface = cutset
            .residue_selector
            .lu_cut
            .as_ref()
            .map(|lu_cut| lu_cut.esurface_ids[0]);
        if let Some(lu_cut) = cutset.residue_selector.lu_cut.as_ref() {
            residues = residues
                .into_iter()
                .flat_map(|expression| {
                    expression.select_esurface_residue_with_cut_edges(
                        lu_cut,
                        &cutset.residue_selector.lu_cut_edge_sets,
                    )
                })
                .collect();
        }
        if representation == ThreeDRepresentation::Ltd {
            for residue in &mut residues {
                if let Some(esurface_id) = lu_cut_representative_esurface {
                    localize_numerator_energy_maps_on_esurface(residue, esurface_id)?;
                }
                remove_ltd_global_contact_completions_from_local_residue(residue);
            }
        }

        residues
            .into_iter()
            .map(|residue| {
                let mut atom = graph.three_d_expression_parametric_atom_with_numerator_gs(
                    &residue,
                    &numerator,
                    match representation {
                        ThreeDRepresentation::Cff => RepresentationMode::Cff,
                        ThreeDRepresentation::Ltd => RepresentationMode::Ltd,
                    },
                    true,
                    &settings.orientation_pattern,
                );
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
                Ok(atom
                    .replace(GS.dim)
                    .with(4)
                    .simplify_color()
                    .replace(GS.den(W_.a_, W_.b_, W_.c_, W_.d_))
                    .with(W_.d_)
                    .expand_dots()?
                    .collect_factors())
            })
            .collect()
    }

    pub(crate) fn root_expanded_4d(
        &mut self,
        graph: &mut Graph,
        cuts: &CutSet,
        valid_orientations: &[EdgeVec<Orientation>],
        settings: &crate::settings::global::GenerationSettings,
        representation: ThreeDRepresentation,
        root_expression: Option<&CFFExpression<OrientationID>>,
    ) -> Result<()> {
        self.initialize_filtered_integrated_uv_root(&settings.uv);
        self.simple_approx = Some(SimpleApprox::root(self.spinney.subgraph.clone()));
        if settings.uv.only_integrated {
            self.integrated_4d = ApproxOp::Root;
        } else {
            self.integrated_4d = ApproxOp::Root;
            self.local_4d_expanded = Some(Expanded4DApprox::root(graph));
            if Self::filtered_integrated_uv_mode_is_active(&settings.uv) {
                self.final_integrand = Some(Self::zero_terms(1));
            } else {
                let root_expression = root_expression.ok_or_else(|| {
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
                    representation == ThreeDRepresentation::Cff
                        && !settings.explicit_orientation_sum_only,
                )?);
            }
        }

        Ok(())
    }

    pub(crate) fn new(spinney: Spinney) -> Approximation {
        Approximation {
            spinney,
            topo_order: 0,
            final_integrand: None,
            simple_approx: None,
            local_3d: CFFapprox::NotComputed,
            local_4d_expanded: None,
            integrated_4d: ApproxOp::NotComputed,
            generate_only_integrated_uv_chain_matches: false,
        }
    }

    #[instrument(skip_all)]
    pub(crate) fn compute_integrated(
        &mut self,
        graph: &Graph,
        vakint: (&Vakint, &vakint::VakintSettings),
        dependent: &Self,
        settings: &UVgenerationSettings,
    ) -> Result<()> {
        if Self::filtered_integrated_uv_mode_is_active(settings)
            && !self.generate_only_integrated_uv_chain_matches
        {
            self.integrated_4d = ApproxOp::Dependent {
                t_arg: IntegrandExpr {
                    integrands: vec![Atom::Zero],
                },
                sign: Sign::Positive,
                subgraph: unsafe {
                    InternalSubGraph::new_unchecked(self.reduced_subgraph(dependent))
                },
            };
            return Ok(());
        }

        let ctx = UVCtx { graph, settings };

        let Some((current, sign)) = &dependent.integrated_4d.expr() else {
            return Err(eyre!("integrated_4d not computed"));
        };

        debug!(pole_part = %settings.pole_part,
            simple = % self.simple_approx
                .as_ref()
                .unwrap()
                .expr(&graph.full_filter()).log_print(None),
            "Computing Integrated",
        );

        let integrands = current
            .iter()
            .map(|a| Integrated::new(vakint.0, vakint.1).kernel(&ctx, self, dependent, a))
            .collect::<Result<Vec<_>>>()?;

        self.integrated_4d = ApproxOp::Dependent {
            t_arg: IntegrandExpr { integrands },
            sign: -*sign,
            subgraph: unsafe { InternalSubGraph::new_unchecked(self.reduced_subgraph(dependent)) },
        };

        Ok(())
    }

    /// Computes the 3d approximation of the UV
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn compute(
        &mut self,
        graph: &mut Graph,
        cuts: &CutSet,
        dependent: &Self,
        valid_orientations: &[EdgeVec<Orientation>],
        settings: &UVgenerationSettings,
        explicit_orientation_sum_only: bool,
    ) -> Result<()> {
        let Some((cff, sign)) = dependent.local_3d.expr() else {
            return Err(eyre!("dependent local 3D approximation was not computed"));
        };
        let finite = finite_integrated_part(&dependent.integrated_4d)?;
        debug!(
            "Integrated 4d finite part: {:#}",
            finite.printer(LOGPRINTOPTS)
        );

        let t_arg = localized_integrated_reduced_factor(
            graph,
            dependent.spinney.filter(),
            cuts,
            valid_orientations,
            explicit_orientation_sum_only,
        )?;

        let ctx = UVCtx { graph, settings };

        if Self::filtered_integrated_uv_mode_is_active(settings) {
            self.set_zero_local_3d(-sign, cff.len());
            self.final_integrand = Some(if self.should_keep_only_integrated_uv_terms(settings) {
                self.final_integrand(
                    graph,
                    cuts,
                    valid_orientations,
                    settings,
                    explicit_orientation_sum_only,
                )?
            } else {
                Self::zero_terms(cff.len())
            });
            return Ok(());
        }

        if cff.len() != t_arg.integrands.len() {
            return Err(eyre!(
                "local 3D UV approximation produced {} residue integrands, but localized finite integrated UV produced {}",
                cff.len(),
                t_arg.integrands.len()
            ));
        }

        let mut integrands = vec![];
        for (local, t_arg) in cff.into_iter().zip(t_arg.integrands) {
            let mut sum_3d = Atom::Zero;
            sum_3d += Local3DApproximation {}.kernel(&ctx, &*self, dependent, &local)?;
            sum_3d -=
                Local3DApproximation {}.kernel(&ctx, &*self, dependent, &(&finite * t_arg))?;
            integrands.push(sum_3d);
        }

        self.local_3d = CFFapprox::Dependent {
            sign: -sign,
            t_arg: IntegrandExpr { integrands },
        };

        self.final_integrand = Some(self.final_integrand(
            graph,
            cuts,
            valid_orientations,
            settings,
            explicit_orientation_sum_only,
        )?);
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn compute_expanded_4d(
        &mut self,
        graph: &mut Graph,
        cuts: &CutSet,
        dependent: &Self,
        valid_orientations: &[EdgeVec<Orientation>],
        settings: &crate::settings::global::GenerationSettings,
        representation: crate::settings::global::ThreeDRepresentation,
    ) -> Result<()> {
        let Some(parent) = dependent.local_4d_expanded.as_ref() else {
            return Err(eyre!("expanded 4D local UV parent was not computed"));
        };
        let (parent_terms, parent_sign) = parent.expr();

        if Self::filtered_integrated_uv_mode_is_active(&settings.uv) {
            self.set_zero_local_4d_expanded(-parent_sign, parent_terms.len());
            self.final_integrand =
                Some(if self.should_keep_only_integrated_uv_terms(&settings.uv) {
                    self.final_integrand_expanded_4d(
                        graph,
                        cuts,
                        valid_orientations,
                        settings,
                        representation,
                    )?
                } else {
                    Self::zero_terms(parent_terms.len())
                });
            return Ok(());
        }

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
        )?);
        Ok(())
    }

    #[instrument(skip_all)]
    pub(crate) fn final_integrand(
        &self,
        graph: &mut Graph,
        cutset: &CutSet,
        valid_orientations: &[EdgeVec<Orientation>],
        _settings: &UVgenerationSettings,
        explicit_orientation_sum_only: bool,
    ) -> Result<Vec<Atom>> {
        let global_num = graph.global_atom();
        let (t, s) = self
            .local_3d
            .expr()
            .ok_or(eyre!("Local3d not yet computed"))?;
        let finite = finite_integrated_part(&self.integrated_4d)?;

        debug!(
            "Integrated 4d finite part: {:#}",
            finite.printer(LOGPRINTOPTS)
        );

        let t_arg = localized_integrated_reduced_factor(
            graph,
            self.spinney.filter(),
            cutset,
            valid_orientations,
            explicit_orientation_sum_only,
        )?;
        if t.len() != t_arg.integrands.len() {
            return Err(eyre!(
                "local 3D UV approximation produced {} residue integrands, but localized finite integrated UV produced {}",
                t.len(),
                t_arg.integrands.len()
            ));
        }

        let reduced = graph
            .full_filter()
            .subtract(self.spinney.subgraph.included())
            .subtract(&graph.initial_state_cut);

        let mut integrands = vec![];

        for (t, t_arg) in t.iter().zip(t_arg.integrands.iter()) {
            let mut cff = s * t.clone() - s * (&finite * t_arg);

            for (p, eid, _) in graph.as_ref().iter_edges_of(&reduced) {
                let eid = usize::from(eid) as i64;
                if p.is_paired() {
                    cff = cff
                        .replace(function!(GS.energy, eid))
                        .with(function!(GS.ose, eid));
                }
            }

            cff = cff.replace(function!(GS.ose, W_.a__, W_.e_)).with(W_.e_);

            let mut resnum = graph
                .numerator(&reduced, self.spinney.subgraph.included())
                .get_single_atom()
                .map_err(|error| eyre!("graph numerator is not a single symbolic atom: {error}"))?;

            let bridgeless_reduced = reduced.subtract(&graph.tree_edges);

            let mut reps = Vec::new();
            // only put edges onshell if they are part of a loop
            for (p, eid, _) in graph.as_ref().iter_edges_of(&bridgeless_reduced) {
                if p.is_paired() {
                    reps.push(GS.add_parametric_sign(eid));
                }
            }

            resnum = resnum.replace_multiple(&reps);
            resnum *= cff * &global_num;

            resnum = resnum.replace(GS.dim).with(4).simplify_color(); //.to_dots();
            // println!(
            //     "Resnum {}",
            //     resnum.printer(
            //         SpensoPrintSettings {
            //             with_dim: true,
            //             ..SpensoPrintSettings::compact()
            //         }
            //         .nice_symbolica()
            //     )
            // );

            resnum = resnum.expand_dots()?;

            let debug_preview = resnum
                .replace(function!(GS.theta, W_.a_))
                .with(Atom::one())
                .replace(GS.m_uv)
                .with(Atom::var(symbol!("gammalooprs::m_uv_preview")))
                .replace(parse_lit!(UFO::mass_scalar_1))
                .with(Atom::var(symbol!("gammalooprs::mass_scalar_1_preview")))
                .collect_factors()
                .collect_num();

            debug!(
                "Integrand before parsing for {} for dod{}:{}",
                self.simple_approx
                    .as_ref()
                    .map(|simple| simple.expr(&graph.full_filter()))
                    .unwrap_or_else(|| Atom::var(symbol!("missing_simple_uv_approximation"))),
                self.spinney.dod,
                // orientations
                //     .first()
                //     .unwrap()
                //     .select(&
                debug_preview.log_print(None) // printer(LOGPRINTOPTS)
            );

            integrands.push(resnum.replace_multiple(&reps))
        }

        // debug!("final_cff {res:>}");
        Ok(integrands)
    }

    #[instrument(skip_all)]
    pub(crate) fn final_integrand_expanded_4d(
        &self,
        graph: &mut Graph,
        cutset: &CutSet,
        valid_orientations: &[EdgeVec<Orientation>],
        settings: &crate::settings::global::GenerationSettings,
        representation: crate::settings::global::ThreeDRepresentation,
    ) -> Result<Vec<Atom>> {
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
            &graph.empty_subgraph::<SuBitGraph>(),
        )?;

        if finite.is_zero() {
            return Ok(integrands);
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
            self.spinney.filter(),
        )?;
        if integrands.len() != finite_integrands.len() {
            if integrands.len() == 1 && integrands[0].is_zero() {
                integrands = Self::zero_terms(finite_integrands.len());
            } else if finite_integrands.len() == 1 && finite_integrands[0].is_zero() {
                finite_integrands = Self::zero_terms(integrands.len());
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
        Ok(integrands)
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
    let (mut t4, t4_sign) = match integrated_4d {
        ApproxOp::Root | ApproxOp::NotComputed => return Ok(Atom::Zero),
        ApproxOp::Union { .. } => {
            return Err(eyre!(
                "cannot extract finite integrated UV part from an unresolved union approximation"
            ));
        }
        ApproxOp::Dependent { t_arg, sign, .. } => (t_arg.integrands.clone(), *sign),
    };

    if t4.len() != 1 {
        return Err(eyre!("expected one integrated 4D approximation term"));
    }
    let finite_term = t4
        .pop()
        .map(|t4| t4_sign * t4)
        .ok_or_else(|| eyre!("expected one integrated 4D approximation term"))?;
    Ok(finite_term
        .series(vakint_symbol!("ε"), Atom::Zero, 0.into(), true)
        .map_err(|error| eyre!("finite integrated UV epsilon expansion failed: {error}"))?
        .coefficient((0, 1).into())
        .replace(GS.m_uv_int)
        .with(GS.m_uv)
        .map_mink_dim(4)
        .replace(function!(symbol!("vakint::g"), W_.a__))
        .with(function!(symbol!("spenso::g"), W_.a__)))
}

impl ApproxOp {
    pub(crate) fn expr(&self) -> Option<(Vec<Atom>, Sign)> {
        match self {
            ApproxOp::NotComputed => None,
            ApproxOp::Union { .. } => {
                //Never gets hit now
                panic!(
                    "Should not call expr on a union approximation, need to compute the union first"
                );
            }
            ApproxOp::Dependent { t_arg, sign, .. } => Some((t_arg.integrands.clone(), *sign)),
            ApproxOp::Root => Some((vec![Atom::num(1)], Sign::Positive)),
        }
    }
}
