use crate::{
    cff::CutCFFIndex,
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
use eyre::eyre;
use gammaloop_tracing_filter::{LogMessage, debug_instrument};
use std::{collections::BTreeMap, hash::Hash};
use tracing::debug;

use symbolica::{
    atom::{Atom, AtomOrView, Symbol},
    function, symbol,
};

use linnet::half_edge::involution::{EdgeIndex, EdgeVec, Orientation};
use linnet::half_edge::subgraph::{InternalSubGraph, SuBitGraph, SubSetLike, SubSetOps};

use vakint::Vakint;
// use vakint::{EvaluationOrder, LoopNormalizationFactor, Vakint, VakintSettings};
use super::{IntegrandExpr, UltravioletGraph};

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
    pub(crate) fn expr(&self) -> Option<(BTreeMap<CutCFFIndex, Atom>, Sign)> {
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
        orientation_pattern: &OrientationPattern,
    ) -> Result<CFFapprox> {
        let cff = graph
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

        Ok(CFFapprox::Dependent {
            sign: Sign::Positive,
            t_arg: IntegrandExpr {
                integrands: cff
                    .iter()
                    .map(|(index, a)| (*index, a * &fourddenoms))
                    .collect(),
            },
        })
    }

    pub(crate) fn root(
        graph: &mut Graph,
        cuts: &CutSet,
        settings: &UVgenerationSettings,
        orientation_pattern: &OrientationPattern,
    ) -> Result<CFFapprox> {
        Self::dependent(
            graph,
            &graph.empty_subgraph::<SuBitGraph>(),
            cuts,
            settings,
            orientation_pattern,
        )
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
            !settings
                .add_integrated_uv_with_n_loops
                .is_some_and(|target_loop_count| {
                    !dependent.excluded_by_filters
                        && graph.n_loops(self.spinney.filter()) == target_loop_count
                });
    }

    fn zero_terms(allowed_keys: &[CutCFFIndex]) -> BTreeMap<CutCFFIndex, Atom> {
        allowed_keys.iter().map(|&key| (key, Atom::Zero)).collect()
    }

    pub(crate) fn root(
        &mut self,
        graph: &mut Graph,
        cuts: &CutSet,
        valid_orientations: &[EdgeVec<Orientation>],
        settings: &UVgenerationSettings,
        orientation_pattern: &OrientationPattern,
    ) -> Result<()> {
        self.excluded_by_filters = false;
        self.simple_approx = Some(SimpleApprox::root(self.spinney.subgraph.clone()));
        if settings.only_integrated {
            self.integrated_4d = ApproxOp::Root;
        } else {
            self.integrated_4d = ApproxOp::Root;
            self.local_3d = CFFapprox::root(graph, cuts, settings, orientation_pattern)?;
            if settings.generate_integrated {
                let (_, _) = self
                    .local_3d
                    .expr()
                    .expect("root local CFF should have been computed");
                self.final_integrand = Some(Self::zero_terms(
                    &cuts.residue_selector.generate_allowed_keys(),
                ));
            } else {
                let (local_terms, local_sign) = self
                    .local_3d
                    .expr()
                    .expect("root local CFF should have been computed");
                self.final_integrand = Some(
                    FinalIntegrand::new(valid_orientations, orientation_pattern).build(
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
                t_arg: IntegrandExpr {
                    integrands: BTreeMap::from([(CutCFFIndex::new_all_none(), Atom::Zero)]),
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
            .collect::<Result<BTreeMap<_, _>>>()?;

        self.integrated_4d = ApproxOp::Dependent {
            t_arg: IntegrandExpr { integrands },
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
        valid_orientations: &[EdgeVec<Orientation>],
        settings: &UVgenerationSettings,
        orientation_pattern: &OrientationPattern,
    ) -> Result<()> {
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
        let final_integrand = FinalIntegrand::new(valid_orientations, orientation_pattern);

        let localize_started = std::time::Instant::now();
        let integrated_t = final_integrand.localized_integrated_ct(
            graph,
            dependent,
            &dependent.integrated_4d,
            cuts,
        )?;
        debug_tags!(#generation, #profile, #uv, #graph, #summary;
            stage = "compute_local_3d_localize_integrated_done",
            localized_count = integrated_t.integrands.len(),
            elapsed_ms = localize_started.elapsed().as_secs_f64() * 1000.0,
            total_elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Localized integrated UV CT for local 3D subtraction"
        );

        let ctx = UVCtx { graph, settings };

        let mut integrands = BTreeMap::new();
        for ((index_local, local), (index_integ, integ)) in
            cff.into_iter().zip(integrated_t.integrands)
        {
            let term_started = std::time::Instant::now();
            if index_local != index_integ {
                return Err(eyre!(
                    "Mismatched indices for local and integrated approximations: {:?} vs {:?}",
                    index_local,
                    index_integ
                ));
            }
            let mut sum_3d = Atom::Zero;

            debug_tags!(#generation, #profile, #uv, #graph, #term, #summary;
                "adding localT(expr)"
            );
            sum_3d += Local3DApproximation {}.kernel(&ctx, &*self, dependent, &local)?;
            debug_tags!(#generation, #profile, #uv, #graph, #term, #summary;
                "subtracting localT(integrated(expr))"
            );
            sum_3d -= Local3DApproximation {}.kernel(&ctx, &*self, dependent, &integ)?;
            integrands.insert(index_local, sum_3d);
            debug_tags!(#generation, #profile, #uv, #graph, #term, #summary;
                stage = "compute_local_3d_term_done",
                cut_index = ?index_local,
                elapsed_ms = term_started.elapsed().as_secs_f64() * 1000.0,
                total_elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
                "Computed local 3D UV CT term"
            );
        }

        self.local_3d = CFFapprox::Dependent {
            sign: -sign,
            t_arg: IntegrandExpr { integrands },
        };

        let final_started = std::time::Instant::now();
        self.final_integrand = Some({
            let (local_terms, local_sign) = self
                .local_3d
                .expr()
                .ok_or(eyre!("Local3d not yet computed"))?;
            final_integrand.build(
                graph,
                self,
                &local_terms,
                local_sign,
                &self.integrated_4d,
                cuts,
            )?
        });
        debug_tags!(#generation, #profile, #uv, #graph, #summary;
            stage = "compute_local_3d_done",
            final_integrand_elapsed_ms = final_started.elapsed().as_secs_f64() * 1000.0,
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Computed local 3D UV CT"
        );
        Ok(())
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
    pub(crate) fn expr(&self) -> Option<(BTreeMap<CutCFFIndex, Atom>, Sign)> {
        match self {
            ApproxOp::NotComputed => None,
            ApproxOp::Union { .. } => {
                //Never gets hit now
                panic!(
                    "Should not call expr on a union approximation, need to compute the union first"
                );
            }
            ApproxOp::Dependent { t_arg, sign, .. } => Some((t_arg.integrands.clone(), *sign)),
            ApproxOp::Root => Some((
                BTreeMap::from([(CutCFFIndex::new_all_none(), Atom::num(1))]),
                Sign::Positive,
            )),
        }
    }
}
