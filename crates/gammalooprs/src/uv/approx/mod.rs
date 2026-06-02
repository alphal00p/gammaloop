use crate::{
    cff::{CutCFFIndex, ResidueSelectedTerms},
    debug_tags,
    graph::{Graph, LoopMomentumBasis, cuts::CutSet},
    momentum::Sign,
    settings::global::OrientationPattern,
    utils::{GS, symbolica_ext::LogPrint},
    uv::{
        ApproximationType, Spinney, UVgenerationSettings,
        approx::{
            final_integrand::FinalIntegrand, integrated::Integrated, local_3d::Local3DApproximation,
        },
    },
};
use color_eyre::Result;
use eyre::{WrapErr, eyre};
use gammaloop_tracing_filter::{LogMessage, debug_instrument};
use std::hash::Hash;
use tracing::debug;

use symbolica::{
    atom::{Atom, AtomOrView, Symbol},
    function, symbol,
};

use linnet::half_edge::involution::{EdgeIndex, EdgeVec, Orientation};
use linnet::half_edge::subgraph::{InternalSubGraph, SuBitGraph, SubSetLike, SubSetOps};

use vakint::Vakint;
// use vakint::{EvaluationOrder, LoopNormalizationFactor, Vakint, VakintSettings};
use super::UltravioletGraph;

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
        sign: Sign,
        t_args: Vec<ResidueSelectedTerms>,
        subgraphs: Vec<InternalSubGraph>,
    },
    Dependent {
        sign: Sign,
        integrands: ResidueSelectedTerms,
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
        function!(GS.t_op, self.expr(bigger_graph))
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
    pub final_integrand: Option<ResidueSelectedTerms>,
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
pub enum CFFapprox {
    NotComputed,
    Dependent {
        sign: Sign,
        integrands: ResidueSelectedTerms,
    },
}

#[derive(Clone)]
pub struct CutStructure {
    pub cuts: Vec<CutSet>,
}

#[derive(Clone, Copy)]
pub(crate) struct OrientationProjection<'a> {
    pub(crate) valid_orientations: &'a [EdgeVec<Orientation>],
    pub(crate) orientation_pattern: &'a OrientationPattern,
}

#[derive(Clone, Copy)]
pub(crate) struct ResidueProjection<'a> {
    pub(crate) cutset: &'a CutSet,
    pub(crate) orientation: OrientationProjection<'a>,
}

impl<'a> OrientationProjection<'a> {
    pub(crate) fn new(
        valid_orientations: &'a [EdgeVec<Orientation>],
        orientation_pattern: &'a OrientationPattern,
    ) -> Self {
        Self {
            valid_orientations,
            orientation_pattern,
        }
    }

    pub(crate) fn for_cut(self, cutset: &'a CutSet) -> ResidueProjection<'a> {
        ResidueProjection {
            cutset,
            orientation: self,
        }
    }
}

impl CutStructure {
    pub(crate) fn empty(graph: &Graph) -> Self {
        Self {
            cuts: vec![CutSet::empty(graph.n_hedges())],
        }
    }
}

impl CFFapprox {
    pub(crate) fn expr(&self) -> Option<(ResidueSelectedTerms, Sign)> {
        match self {
            CFFapprox::NotComputed => None,
            CFFapprox::Dependent { sign, integrands } => Some((integrands.clone(), *sign)),
        }
    }

    pub(crate) fn dependent(
        graph: &mut Graph,
        to_contract: &SuBitGraph,
        projection: ResidueProjection<'_>,
    ) -> Result<CFFapprox> {
        let cff = graph
            .cff(
                &to_contract
                    .union(&graph.tree_edges)
                    .subtract(&graph.initial_state_cut),
                projection.cutset,
                projection.orientation.orientation_pattern,
            )?
            .expression_with_selectors();

        let fourddenoms = GS.wrap_tree_denoms(
            graph.denominator(&graph.tree_edges.subtract(&graph.initial_state_cut), |_| -1),
        );

        Ok(CFFapprox::Dependent {
            sign: Sign::Positive,
            integrands: cff.map(|a| a * &fourddenoms),
        })
    }

    pub(crate) fn root(graph: &mut Graph, projection: ResidueProjection<'_>) -> Result<CFFapprox> {
        let empty = graph.empty_subgraph::<SuBitGraph>();
        Self::dependent(graph, &empty, projection)
    }
}

impl Approximation {
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

    fn filtered_integrated_uv_mode_is_active(settings: &UVgenerationSettings) -> bool {
        settings.generate_integrated && settings.add_integrated_uv_with_n_loops.is_some()
    }

    fn zero_terms(allowed_keys: &[CutCFFIndex]) -> ResidueSelectedTerms {
        ResidueSelectedTerms::zero_for_keys(allowed_keys)
    }

    fn set_zero_local_3d(&mut self, sign: Sign, allowed_keys: &[CutCFFIndex]) {
        self.local_3d = CFFapprox::Dependent {
            sign,
            integrands: Self::zero_terms(allowed_keys),
        };
    }

    fn should_keep_only_integrated_uv_terms(&self, settings: &UVgenerationSettings) -> bool {
        Self::filtered_integrated_uv_mode_is_active(settings) && !self.excluded_by_filters
    }

    pub(crate) fn root(
        &mut self,
        graph: &mut Graph,
        cuts: &CutSet,
        orientation: OrientationProjection<'_>,
        settings: &UVgenerationSettings,
    ) -> Result<()> {
        let projection = orientation.for_cut(cuts);
        self.excluded_by_filters = false;
        self.simple_approx = Some(SimpleApprox::root(self.spinney.subgraph.clone()));
        if settings.only_integrated {
            self.integrated_4d = ApproxOp::Root;
        } else {
            self.integrated_4d = ApproxOp::Root;
            self.local_3d = CFFapprox::root(graph, projection)?;
            if Self::filtered_integrated_uv_mode_is_active(settings) {
                self.final_integrand = Some(Self::zero_terms(
                    &cuts.residue_selector.generate_allowed_keys(),
                ));
            } else {
                let (local_terms, local_sign) = self
                    .local_3d
                    .expr()
                    .expect("root local CFF should have been computed");
                let uv_marker = (settings.add_sigma && settings.keep_sigma).then(|| {
                    self.simple_approx
                        .as_ref()
                        .unwrap()
                        .expr(&graph.full_filter())
                });
                self.final_integrand = Some(
                    FinalIntegrand::new(projection, uv_marker.as_ref(), true).build(
                        graph,
                        self,
                        &local_terms,
                        local_sign,
                        &self.integrated_4d,
                    )?,
                );
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
                integrands: ResidueSelectedTerms::all_none(Atom::Zero),
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
                debug_tags!(#generation, #profile, #uv, #graph, #term, #summary;
                    stage = "compute_integrated_kernel_done",
                    cut_index = ?index,
                    elapsed_ms = term_started.elapsed().as_secs_f64() * 1000.0,
                    "Computed integrated UV CT term"
                );
                Ok((*index, integrated))
            })
            .collect::<Result<ResidueSelectedTerms>>()?;

        self.integrated_4d = ApproxOp::Dependent {
            integrands,
            sign: -*sign,
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
    #[allow(clippy::too_many_arguments)]
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
        orientation: OrientationProjection<'_>,
        settings: &UVgenerationSettings,
    ) -> Result<()> {
        let projection = orientation.for_cut(cuts);
        let started = std::time::Instant::now();
        debug_tags!(#generation, #profile, #uv, #graph, #summary;
            stage = "compute_local_3d_start",
            "Computing local 3D UV CT"
        );
        let Some((cff, sign)) = dependent.local_3d.expr() else {
            panic!("Should have computed the dependent cff");
        };
        debug_tags!(#generation, #profile, #uv, #graph, #summary;
            stage = "compute_local_3d_parent_expr_done",
            cff_count = cff.len(),
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Loaded parent local 3D UV CT"
        );
        let final_integrand = FinalIntegrand::new(projection, None, false);

        let localize_started = std::time::Instant::now();
        let integrated_t =
            final_integrand.localized_integrated_ct(graph, dependent, &dependent.integrated_4d)?;
        debug_tags!(#generation, #profile, #uv, #graph, #summary;
            stage = "compute_local_3d_localize_integrated_done",
            localized_count = integrated_t.active.len(),
            elapsed_ms = localize_started.elapsed().as_secs_f64() * 1000.0,
            total_elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Localized integrated UV CT for local 3D subtraction"
        );

        let ctx = UVCtx { graph, settings };
        let uv_marker = (settings.add_sigma && settings.keep_sigma).then(|| {
            self.simple_approx
                .as_ref()
                .unwrap()
                .expr(&ctx.graph.full_filter())
        });

        if Self::filtered_integrated_uv_mode_is_active(settings) {
            self.set_zero_local_3d(
                -sign,
                cuts.residue_selector.generate_allowed_keys().as_slice(),
            );
            self.final_integrand = Some(if self.should_keep_only_integrated_uv_terms(settings) {
                self.final_integrand(graph, projection)?
            } else {
                Self::zero_terms(&cuts.residue_selector.generate_allowed_keys())
            });
            return Ok(());
        }

        let mut tagged_integrands = uv_marker.as_ref().map(|_| ResidueSelectedTerms::empty());
        let output_sign = -sign;
        let integrands = cff
            .try_zip_with(&integrated_t.active, |index, local, integ| {
                let term_started = std::time::Instant::now();
                let frozen = integrated_t.frozen_integrands.get(index)?;
                let mut sum_3d = Atom::Zero;

                debug_tags!(#generation, #profile, #uv, #integrated, #local, #graph, #term, #trace;
                    stage = "compute_local_3d_term_inputs",
                    cut_index = ?index,
                    log.local = local,
                    log.localized_integrated = integ,
                    log.frozen = frozen,
                    "Local 3D UV CT term inputs"
                );

                debug_tags!(#generation, #profile, #uv, #graph, #term, #summary;
                    "adding localT(expr)"
                );
                let local_ct =
                    Local3DApproximation::full().kernel(&ctx, &*self, dependent, local)?;
                debug_tags!(#generation, #profile, #uv, #local, #graph, #term, #trace;
                    stage = "compute_local_3d_term_local_done",
                    cut_index = ?index,
                    log.local_ct = local_ct,
                    "Computed localT(expr)"
                );
                debug_tags!(#generation, #profile, #uv, #graph, #term, #summary;
                    "subtracting localT(integrated(expr))"
                );
                let localized_integrated_ct = Local3DApproximation::reduced()
                    .kernel(&ctx, &*self, dependent, integ)?
                    * frozen;
                debug_tags!(#generation, #profile, #uv, #integrated, #local, #graph, #term, #trace;
                    stage = "compute_local_3d_term_integrated_done",
                    cut_index = ?index,
                    log.localized_integrated_ct = localized_integrated_ct,
                    "Computed localT(integrated(expr))"
                );
                if let (Some(marker), Some(tagged_integrands)) =
                    (uv_marker.as_ref(), tagged_integrands.as_mut())
                {
                    let tagged_sum_3d = local_ct.clone() * function!(GS.uv_local, marker.clone())
                        - localized_integrated_ct.clone()
                            * function!(GS.uv_integrated, marker.clone());
                    tagged_integrands.add_assign_term(index, tagged_sum_3d.clone());
                    debug_tags!(#generation, #profile, #uv, #local, #graph, #term, #trace;
                        stage = "compute_local_3d_term_tagged_output",
                        cut_index = ?index,
                        log.tagged_local_3d_ct = tagged_sum_3d,
                        "Computed tagged local 3D UV CT expression"
                    );
                }
                sum_3d += local_ct;
                sum_3d -= localized_integrated_ct;
                debug_tags!(#generation, #profile, #uv, #local, #graph, #term, #trace;
                    stage = "compute_local_3d_term_output",
                    cut_index = ?index,
                    log.local_3d_ct = sum_3d,
                    "Computed local 3D UV CT expression"
                );
                debug_tags!(#generation, #profile, #uv, #graph, #term, #summary;
                    stage = "compute_local_3d_term_done",
                    cut_index = ?index,
                    elapsed_ms = term_started.elapsed().as_secs_f64() * 1000.0,
                    total_elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
                    "Computed local 3D UV CT term"
                );
                Ok(sum_3d)
            })
            .wrap_err("while combining legacy local and integrated approximations")?;
        self.local_3d = CFFapprox::Dependent {
            sign: output_sign,
            integrands,
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
            FinalIntegrand::new(projection, uv_marker.as_ref(), false).build(
                graph,
                self,
                final_local_terms,
                local_sign,
                &self.integrated_4d,
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

    pub(crate) fn final_integrand(
        &self,
        graph: &mut Graph,
        projection: ResidueProjection<'_>,
    ) -> Result<ResidueSelectedTerms> {
        let (local_terms, local_sign) = self
            .local_3d
            .expr()
            .ok_or(eyre!("Local3d not yet computed"))?;

        FinalIntegrand::new(projection, None, false).build(
            graph,
            self,
            &local_terms,
            local_sign,
            &self.integrated_4d,
        )
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

impl ApproxOp {
    pub(crate) fn expr(&self) -> Option<(ResidueSelectedTerms, Sign)> {
        match self {
            ApproxOp::NotComputed => None,
            ApproxOp::Union { .. } => {
                //Never gets hit now
                panic!(
                    "Should not call expr on a union approximation, need to compute the union first"
                );
            }
            ApproxOp::Dependent {
                integrands, sign, ..
            } => Some((integrands.clone(), *sign)),
            ApproxOp::Root => Some((ResidueSelectedTerms::all_none(Atom::num(1)), Sign::Positive)),
        }
    }
}
