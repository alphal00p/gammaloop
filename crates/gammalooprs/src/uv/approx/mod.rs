use crate::{
    debug_tags,
    graph::{Graph, LoopMomentumBasis, cuts::CutSet},
    momentum::Sign,
    numerator::symbolica_ext::AtomCoreExt,
    utils::{GS, W_, symbolica_ext::LogPrint},
    uv::{
        ApproximationType, Spinney, UVgenerationSettings,
        approx::{integrated::Integrated, local_3d::Local3DApproximation},
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

use linnet::half_edge::involution::{EdgeIndex, EdgeVec, Orientation};
use linnet::half_edge::subgraph::{InternalSubGraph, SuBitGraph, SubSetLike, SubSetOps};

use tracing::instrument;
use vakint::{Vakint, vakint_symbol};
// use vakint::{EvaluationOrder, LoopNormalizationFactor, Vakint, VakintSettings};
use super::{IntegrandExpr, UltravioletGraph};

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
        let cff = graph
            .cff(
                &to_contract
                    .union(&graph.tree_edges)
                    .subtract(&graph.initial_state_cut),
                cuts,
            )?
            .expression_with_selectors();

        let fourddenoms = GS.wrap_tree_denoms(
            graph.denominator(&graph.tree_edges.subtract(&graph.initial_state_cut), |_| -1),
        );

        Ok(CFFapprox::Dependent {
            sign: Sign::Positive,
            t_arg: IntegrandExpr {
                integrands: cff.iter().map(|a| a * &fourddenoms).collect(),
            },
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
    expr: &Atom,
    cuts: &CutSet,
    valid_orientations: &[EdgeVec<Orientation>],
) -> Result<IntegrandExpr> {
    debug_tags!(#uv;graph.name=%graph.name , expr=%expr.log_print(Some(80)),"Localizing integrated UV CT");
    let cff = graph.cff(
        &to_contract
            .union(&graph.tree_edges)
            .subtract(&graph.initial_state_cut),
        cuts,
    )?;

    let fourddenoms = GS.wrap_tree_denoms(
        graph.denominator(&graph.tree_edges.subtract(&graph.initial_state_cut), |_| -1),
    );

    let internal_edges = internal_paired_edges_of_subgraph(graph, to_contract);

    let integrands = cff
        .terms
        .into_iter()
        .map(|term| {
            let mut localized = Atom::Zero;
            for (expr, orientation) in term.expression.into_iter().zip(term.orientations) {
                localized += localize_reduced_orientation_term(
                    &(expr * &fourddenoms),
                    &orientation,
                    valid_orientations,
                    &internal_edges,
                )?;
            }
            Ok(localized * expr)
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

    fn should_keep_only_integrated_uv_terms(&self, settings: &UVgenerationSettings) -> bool {
        Self::filtered_integrated_uv_mode_is_active(settings)
            && self.generate_only_integrated_uv_chain_matches
    }

    pub(crate) fn root(
        &mut self,
        graph: &mut Graph,
        cuts: &CutSet,
        valid_orientations: &[EdgeVec<Orientation>],
        settings: &UVgenerationSettings,
    ) -> Result<()> {
        self.initialize_filtered_integrated_uv_root(settings);
        self.simple_approx = Some(SimpleApprox::root(self.spinney.subgraph.clone()));
        if settings.only_integrated {
            self.integrated_4d = ApproxOp::Root;
        } else {
            self.integrated_4d = ApproxOp::Root;
            self.local_3d = CFFapprox::root(graph, cuts, settings)?;
            if Self::filtered_integrated_uv_mode_is_active(settings) {
                let (integrands, _) = self
                    .local_3d
                    .expr()
                    .expect("root local CFF should have been computed");
                self.final_integrand = Some(Self::zero_terms(integrands.len()));
            } else {
                self.final_integrand =
                    Some(self.final_integrand(graph, cuts, valid_orientations, settings)?);
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
    ) -> Result<()> {
        let Some((cff, sign)) = dependent.local_3d.expr() else {
            panic!("Should have computed the dependent cff");
        };
        let (mut t4, t4_sign) = if let ApproxOp::Root = dependent.integrated_4d {
            (vec![Atom::Zero], Sign::Positive)
        } else {
            dependent
                .integrated_4d
                .expr()
                .unwrap_or((vec![Atom::Zero], Sign::Positive))
        };
        if t4.len() != 1 {
            panic!("Should only have one t_arg for the 4d approximation");
        }
        let finite = t4
            .pop()
            .map(|t4| t4_sign * t4)
            .unwrap()
            .series(vakint_symbol!("ε"), Atom::Zero, 0.into(), true)
            .unwrap()
            .coefficient((0, 1).into())
            .replace(GS.m_uv_int)
            .with(GS.m_uv)
            .map_mink_dim(4);
        debug_tags!(
            #uv;
            finite = %finite.log_print(Some(80)),
            t4 = %t4.iter().enumerate().map(|(i, t4)| format!("t4_{} = {}", i, t4.log_print(Some(80)))).collect::<Vec<_>>().join("\n"),
            "Computing UV subtraction",
        );

        let integrated_t = localized_integrated_reduced_factor(
            graph,
            dependent.spinney.filter(),
            &finite,
            cuts,
            valid_orientations,
        )?;

        let ctx = UVCtx { graph, settings };

        if Self::filtered_integrated_uv_mode_is_active(settings) {
            self.set_zero_local_3d(-sign, cff.len());
            self.final_integrand = Some(if self.should_keep_only_integrated_uv_terms(settings) {
                self.final_integrand(graph, cuts, valid_orientations, settings)?
            } else {
                Self::zero_terms(cff.len())
            });
            return Ok(());
        }

        let mut integrands = vec![];
        for (local, integ) in cff.into_iter().zip(integrated_t.integrands) {
            let mut sum_3d = Atom::Zero;
            sum_3d += Local3DApproximation {}.kernel(&ctx, &*self, dependent, &local)?;
            sum_3d -= Local3DApproximation {}.kernel(&ctx, &*self, dependent, &integ)?;
            integrands.push(sum_3d);
        }

        self.local_3d = CFFapprox::Dependent {
            sign: -sign,
            t_arg: IntegrandExpr { integrands },
        };

        self.final_integrand =
            Some(self.final_integrand(graph, cuts, valid_orientations, settings)?);
        Ok(())
    }

    #[instrument(skip_all)]
    pub(crate) fn final_integrand(
        &self,
        graph: &mut Graph,
        cutset: &CutSet,
        valid_orientations: &[EdgeVec<Orientation>],
        _settings: &UVgenerationSettings,
    ) -> Result<Vec<Atom>> {
        let global_num = graph.global_atom();
        let (t, s) = self
            .local_3d
            .expr()
            .ok_or(eyre!("Local3d not yet computed"))?;
        let (mut t4, t4_sign) = if let ApproxOp::Root = self.integrated_4d {
            (vec![Atom::Zero], Sign::Positive)
        } else {
            self.integrated_4d
                .expr()
                .unwrap_or((vec![Atom::Zero], Sign::Positive))
        };
        if t4.len() != 1 {
            panic!("Should only have one t_arg for the 4d approximation");
        }
        let finite = t4
            .pop()
            .map(|t4| t4_sign * t4)
            .unwrap()
            .series(vakint_symbol!("ε"), Atom::Zero, 0.into(), true)
            .unwrap()
            .coefficient((0, 1).into())
            .replace(GS.m_uv_int)
            .with(GS.m_uv)
            .map_mink_dim(4)
            .replace(function!(symbol!("vakint::g"), W_.a__))
            .with(function!(symbol!("spenso::g"), W_.a__));

        debug_tags!(
            #uv,#final;
            finite = %finite.log_print(Some(80)),
            t4 = %t4.iter().enumerate().map(|(i, t4)| format!("t4_{} = {}", i, t4.log_print(Some(80)))).collect::<Vec<_>>().join("\n"),
            "Computing Final integrand after uv subtraction",
        );

        let integrated_t = localized_integrated_reduced_factor(
            graph,
            self.spinney.filter(),
            &finite,
            cutset,
            valid_orientations,
        )?;

        let reduced = graph
            .full_filter()
            .subtract(self.spinney.subgraph.included())
            .subtract(&graph.initial_state_cut);

        let mut integrands = vec![];

        for (local, integ) in t.iter().zip(integrated_t.integrands.iter()) {
            let mut cff = s * (local - integ);

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
                .unwrap();

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
                    .unwrap()
                    .expr(&graph.full_filter()),
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
