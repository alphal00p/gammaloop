#![allow(dead_code)]

use crate::{
    cff::{CutCFF, expression::GraphOrientation},
    graph::{Graph, LMBext, LoopMomentumBasis, cuts::CutSet},
    momentum::Sign,
    numerator::{aind::Aind, symbolica_ext::AtomCoreExt},
    utils::{
        GS, W_,
        symbolica_ext::{CallSymbol, LOGPRINTOPTS, LogPrint, TypstFormat},
    },
    uv::{UVgenerationSettings, settings::VakintSettings},
};
use ahash::AHashSet;
use color_eyre::Result;
use eyre::eyre;
use idenso::{color::ColorSimplifier, gamma::GammaSimplifier, metric::MetricSimplifier};
use std::{collections::HashSet, hash::Hash};
use tracing::debug;

use spenso::{
    network::{
        library::{TensorLibraryData, symbolic::ETS},
        parsing::SPENSO_TAG,
    },
    structure::{
        representation::{Minkowski, RepName},
        slot::{DummyAind, IsAbstractSlot, Slot},
    },
};
use symbolica::{
    atom::{Atom, AtomCore, FunctionBuilder, Symbol},
    function,
    id::{MatchSettings, Replacement},
    parse, parse_lit,
    solve::SolveError,
    symbol,
};

use linnet::half_edge::{
    HedgeGraph, NodeIndex,
    builder::HedgeGraphBuilder,
    involution::HedgePair,
    subgraph::{
        Inclusion, InternalSubGraph, ModifySubSet, SuBitGraph, SubGraphLike, SubSetLike, SubSetOps,
    },
};

use tracing::{info, instrument};
use vakint::{Vakint, VakintExpression, vakint_symbol};
// use vakint::{EvaluationOrder, LoopNormalizationFactor, Vakint, VakintSettings};

use super::{IntegrandExpr, UltravioletGraph, uv_graph::UVE};

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

    pub(crate) fn union<'a>(
        subgraph: InternalSubGraph,
        union: impl IntoIterator<Item = &'a Self>,
    ) -> Self {
        let mut t_args = vec![];
        let mut sign = Sign::Positive;
        for u in union {
            if u.t_args.len() != 1 {
                panic!("Union can only be applied to dependent approximations");
            }
            t_args.push(u.t_args[0].clone());
            sign = sign * u.sign;
        }
        SimpleApprox {
            t_args,
            sign,
            graph: subgraph,
        }
    }
}

#[instrument(skip_all)]
pub(crate) fn to_vakint_integrand<
    E: UVE,
    V,
    H,
    S: SubGraphLike + SubSetLike<Base = SuBitGraph>,
    SS: SubGraphLike,
>(
    integrand: &Atom,
    graph: &HedgeGraph<E, V, H>,
    reduced: &S,
    dependent_subgraph: &SS,
    settings: &VakintSettings,
    substitute_masses_to_m_uv: bool,
) -> VakintExpression {
    let mut integrand_vakint = integrand.clone();

    debug!(
        integrand = %integrand.log_print(),
        "Integrand vakint init"
    );
    //Atom::Zero

    // strip the momentum wrapper from the denominator
    integrand_vakint = integrand_vakint
        .replace(function!(
            GS.den,
            W_.prop_,
            function!(GS.emr_mom, W_.prop_, W_.mom_),
            W_.x__
        ))
        .with(function!(GS.den, W_.prop_, W_.mom_, W_.x__));

    // println!("Expanded: {:>}", integrand_vakint.expand());

    integrand_vakint = integrand_vakint.expand();

    let mut propagator_id = 1;

    let vk_prop = vakint::symbols::S.prop;
    let vk_edge = vakint_symbol!("edge");
    let vk_topo = vakint_symbol!("topo");

    // let contracted_nodes: BTreeSet<NodeIndex> = graph
    //     .iter_nodes_of(dependent_subgraph)
    //     .map(|(nid, c, v)| nid)
    //     .collect();

    debug!(reduced = %graph.dot(reduced), "Den to prop for");
    // let first = contracted_nodes.first();

    for (pair, index, _data) in graph.iter_edges_of(reduced) {
        if let HedgePair::Paired { source, sink } = pair {
            // let source = if contracted_nodes.contains(&graph.node_id(source)) {
            //     first.unwrap().clone()
            // } else {
            //     graph.node_id(source)
            // };
            // let sink = if contracted_nodes.contains(&graph.node_id(sink)) {
            //     first.unwrap().clone()
            // } else {
            //     graph.node_id(sink)
            // };
            integrand_vakint = integrand_vakint
                .replace(function!(
                    GS.den,
                    usize::from(index) as i64,
                    W_.mom_,
                    W_.mass_,
                    W_.x___
                ))
                .with(function!(
                    vk_prop,
                    propagator_id,
                    function!(
                        vk_edge,
                        usize::from(graph.node_id(source)),
                        usize::from(graph.node_id(sink))
                    ),
                    W_.mom_,
                    if substitute_masses_to_m_uv {
                        GS.m_uv
                    } else {
                        W_.mass_
                    },
                    1
                ))
                .replace(function!(vk_prop, W_.x___, 1).pow(Atom::var(W_.e_)))
                .with(function!(vk_prop, W_.x___, -Atom::var(W_.e_)));
            propagator_id += 1;
        }
    }

    let mut first: Option<NodeIndex> = None;

    debug!(
        reduced = %graph.dot(dependent_subgraph),
        "Shrinking subgraph for vakint"
    );
    // shrink vertices of the subgraph
    for (id, _crown, _data) in graph.iter_nodes_of(dependent_subgraph) {
        debug!(id=%id,"Shrinking Node");

        if let Some(first) = first {
            integrand_vakint = integrand_vakint
                .replace(function!(
                    vk_prop,
                    W_.x_,
                    function!(vk_edge, id.0, W_.y_),
                    W_.x___
                ))
                .with(function!(
                    vk_prop,
                    W_.x_,
                    function!(vk_edge, first.0, W_.y_),
                    W_.x___
                ))
                .replace(function!(
                    vk_prop,
                    W_.x_,
                    function!(vk_edge, W_.y_, id.0),
                    W_.x___
                ))
                .with(function!(
                    vk_prop,
                    W_.x_,
                    function!(vk_edge, W_.y_, first.0),
                    W_.x___
                ))
        } else {
            first = Some(id);
        }
    }

    // flip edges to positive momentum
    // FIXME: how will this work for sums of momenta?
    integrand_vakint = integrand_vakint
        .replace(function!(
            vk_prop,
            W_.x_,
            function!(vk_edge, W_.a_, W_.b_),
            -Atom::var(W_.y_),
            W_.e___
        ))
        .repeat()
        .with(function!(
            vk_prop,
            W_.x_,
            function!(vk_edge, W_.b_, W_.a_),
            W_.y_,
            W_.e___
        ));

    // fuse raised edges
    integrand_vakint = integrand_vakint
        .replace(
            function!(
                vk_prop,
                W_.x_,
                function!(vk_edge, W_.a_, W_.b_),
                W_.x___,
                W_.e_
            ) * function!(
                vk_prop,
                W_.y_,
                function!(vk_edge, W_.b_, W_.c_),
                W_.x___,
                W_.f_
            ),
        )
        .repeat()
        .with(function!(
            vk_prop,
            W_.y_,
            function!(vk_edge, W_.a_, W_.c_),
            W_.x___,
            W_.e_ + W_.f_
        ));

    // println!(
    //     "Integrand pre vakint: {:}",
    //     VakintExpression::try_from(
    //         integrand_vakint
    //             .replace(function!(vk_prop, W_.x__))
    //             .with(function!(vk_topo, function!(vk_prop, W_.x__)))
    //             .replace(function!(vk_topo, W_.x_) * function!(vk_topo, W_.y_))
    //             .repeat()
    //             .with(function!(vk_topo, W_.x_ * W_.y_))
    //     )
    //     .unwrap()
    // );

    let mut a = VakintExpression::try_from(
        integrand_vakint
            .replace(function!(vk_prop, W_.x__))
            .with(function!(vk_topo, function!(vk_prop, W_.x__)))
            .replace(function!(vk_topo, W_.x_) * function!(vk_topo, W_.y_))
            .repeat()
            .with(function!(vk_topo, W_.x_ * W_.y_)),
    )
    .unwrap();

    for (_i, t) in a.0.iter_mut().enumerate() {
        debug!(integral=%t.integral,"Starting integral");

        let mut graph = HedgeGraphBuilder::new();
        //prop(<id>, edge(<node_source>,<node_sink>), <mom>, <mass>, <power>)
        let pat = function!(
            vk_prop,
            W_.a_,
            function!(vakint::symbols::S.edge, W_.i_, W_.j_),
            W_.c_,
            W_.d_,
            W_.e_
        )
        .to_pattern();
        let mut nodemap = std::collections::HashMap::new();

        struct ContractibleEdge {
            mom: Atom,
            mass: Atom,
            power: i32,
        }

        for m in t.integral.pattern_match(&pat, None, None) {
            let i: usize = m[&W_.i_].as_view().try_into().unwrap();
            let j: usize = m[&W_.j_].as_view().try_into().unwrap();
            nodemap.entry(i).or_insert_with(|| graph.add_node(()));
            nodemap.entry(j).or_insert_with(|| graph.add_node(()));

            graph.add_edge(
                nodemap[&i],
                nodemap[&j],
                ContractibleEdge {
                    mom: m[&W_.c_].clone(),
                    mass: m[&W_.d_].clone(),
                    power: m[&W_.e_].as_view().try_into().unwrap(),
                },
                false,
            );
        }

        let mut system = vec![];
        let mut vars = HashSet::new();

        let mut graph: HedgeGraph<ContractibleEdge, ()> = graph.build();

        while let Some(same_mass_two_bond) = graph.a_bond(&|c| {
            let mut count = 0;
            let mut mass = None;
            for (_, _, d) in graph.iter_edges_of(c) {
                count += 1;
                if mass.is_none() {
                    mass = Some(d.data.mass.clone());
                } else if mass.as_ref().unwrap() != &d.data.mass {
                    return false;
                }
                if count > 2 {
                    return false;
                }
            }
            count == 2
        }) {
            let mut iter = same_mass_two_bond.included_iter();
            let first = graph[&iter.next().unwrap()];
            let second = graph[&iter.next().unwrap()];
            graph[first].power += graph[second].power;
            let mut to_contract: SuBitGraph = graph.empty_subgraph();
            to_contract.add(graph[&second].1);
            graph.contract_subgraph(&to_contract, ());
        }
        let mut nodes_to_merge = vec![];

        for c in graph.connected_components(&graph.full_filter()) {
            let Some((nid, _, _)) = graph.iter_nodes_of(&c).next() else {
                continue;
            };
            nodes_to_merge.push(nid);
        }

        if nodes_to_merge.len() > 0 {
            graph.identify_nodes(&nodes_to_merge, ());
        }

        graph.forget_identification_history();
        debug!(graph = %graph.base_dot(), "Graph");

        let mut new_integral: Atom = 1.into();
        for (p, eid, e) in graph.iter_edges() {
            let HedgePair::Paired { source, sink } = p else {
                continue;
            };
            new_integral *= function!(
                vk_prop,
                eid.0 + 1, //vakint propagator ids are 1-indexed
                function!(
                    vakint::symbols::S.edge,
                    graph.node_id(source).0,
                    graph.node_id(sink).0
                ),
                &e.data.mom,
                &e.data.mass.pow(2).replace(GS.m_uv).with(GS.m_uv_int),
                e.data.power
            )
        }

        // println!("{}->{}", t.integral, new_integral);
        t.integral = function!(vakint::symbols::S.topo, new_integral);

        let nloops = graph.cyclotomatic_number(&graph.full_filter());

        let lmb = graph.lmb();
        let mom_pat = function!(GS.emr_mom, W_.a_).to_pattern();
        for (p, e, ed) in graph.iter_edges() {
            if p.is_paired() {
                // println!("{e}");
                let loop_expr = lmb.loop_atom::<Atom>(e, vakint::symbols::S.k, &[], false);

                ed.data
                    .mom
                    .pattern_match(&mom_pat, None, None)
                    .for_each(|m| {
                        let var = mom_pat.replace_wildcards(&m);
                        vars.insert(var);
                    });

                // println!("{loop_expr}");

                // println!("{external_expr}");
                let is_zero = &ed.data.mom - loop_expr;
                // println!("Momentum check for edge {}: {}", e, is_zero);
                system.push(is_zero);
            }
        }

        let vars = vars.into_iter().collect::<Vec<_>>();
        let add_additional_args = [
            Replacement::new(
                function!(GS.emr_mom, W_.i_).to_pattern(),
                function!(GS.emr_mom, W_.i_, W_.a___),
            )
            .with_settings(MatchSettings {
                allow_new_wildcards_on_rhs: true,
                ..Default::default()
            }),
            Replacement::new(
                function!(vakint::symbols::S.k, W_.i_).to_pattern(),
                function!(vakint::symbols::S.k, W_.i_, W_.a___),
            )
            .with_settings(MatchSettings {
                allow_new_wildcards_on_rhs: true,
                ..Default::default()
            }),
        ];
        let a = Atom::solve_linear_system::<u8, _, _>(&system, &vars);
        match a {
            Ok(a) => {
                for (v, k) in a.iter().zip(vars.iter()) {
                    let lhs = k.replace_multiple(&add_additional_args);
                    let rhs = v.replace_multiple(&add_additional_args);
                    // println!("Momentum solution: {} -> {}", lhs, rhs);
                    t.integral = t.integral.replace(lhs.to_pattern()).with(rhs.to_pattern());
                    t.numerator = t.numerator.replace(lhs.to_pattern()).with(rhs.to_pattern());
                }
            }
            Err(SolveError::Underdetermined {
                partial_solution, ..
            }) => {
                let mut reps = vec![];
                for (p, v) in partial_solution.iter().zip(vars.iter()).rev() {
                    let p = p.replace_multiple(&reps);

                    if &p == v {
                        reps.push(Replacement::new(
                            p.replace_multiple(&add_additional_args).to_pattern(),
                            Atom::Zero,
                        ))
                    } else {
                        reps.push(Replacement::new(
                            v.replace_multiple(&add_additional_args).to_pattern(),
                            p.replace_multiple(&add_additional_args).to_pattern(),
                        ))
                    }

                    // println!("Partial solution: {}->{}", v, p);
                }
                // for r in &reps {
                //     println!("Rep: {:#}", r);
                // }

                t.integral = t.integral.replace_multiple(&reps);
                t.numerator = t.numerator.replace_multiple(&reps);
            }
            Err(a) => {
                panic!(
                    "Could not solve momentum system for vakint integrand: {}",
                    a
                );
            }
        }

        // debug!(
        //     "Graph from vakint expression:\n{}\n{}",
        //     graph.dot_lmb_of(&graph.full_filter(), &graph.lmb()),
        //     term
        // );

        t.numerator *= parse!(&settings.additional_normalization).pow(nloops);

        t.numerator = t.numerator.simplify_metrics().to_dots();
        t.integral = t
            .integral
            .replace(function!(GS.emr_mom, W_.x___))
            .with(function!(vakint::symbols::S.p, W_.x___));
        t.numerator = t
            .numerator
            .replace(function!(GS.emr_mom, W_.x___))
            .with(function!(vakint::symbols::S.p, W_.x___))
            .replace(function!(
                SPENSO_TAG.dot,
                Minkowski {}.new_rep(GS.dim).to_symbolic([]),
                W_.a_,
                W_.b_
            ))
            .with(vakint::symbols::S.dot(W_.a_, W_.b_))
            .replace(function!(
                ETS.metric,
                Minkowski {}.to_symbolic([W_.a__]),
                Minkowski {}.to_symbolic([W_.b__])
            ))
            .with(function!(
                vakint::symbols::S.metric,
                Minkowski {}.to_symbolic([W_.a__]),
                Minkowski {}.to_symbolic([W_.b__])
            ));
    }

    a
}

#[derive(Clone)]
pub struct Approximation {
    // The union of all spinneys, remaining graph is full graph minus subgraph
    pub subgraph: InternalSubGraph,
    pub dod: i32,
    pub lmb: LoopMomentumBasis,
    pub local_3d: CFFapprox, //3d denoms
    pub final_integrand: Option<Atom>,
    pub integrated_4d: ApproxOp, //4d
    pub integrated_pole_part: ApproxOp,
    pub simple_approx: Option<SimpleApprox>,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum CFFapprox {
    NotComputed,
    Dependent { sign: Sign, t_arg: IntegrandExpr },
}

pub struct CutStructure {
    pub cuts: Vec<CutSet>,
}

impl CFFapprox {
    pub(crate) fn expr(&self) -> Option<(Atom, Sign)> {
        match self {
            CFFapprox::NotComputed => None,
            CFFapprox::Dependent { sign, t_arg } => Some((t_arg.integrand.clone(), *sign)),
        }
    }

    pub(crate) fn dependent(
        graph: &Graph,
        to_contract: &SuBitGraph,
        cuts: &CutSet,
        _settings: &UVgenerationSettings,
    ) -> Result<CFFapprox> {
        let cff = graph
            .cff(&to_contract.union(&graph.tree_edges), cuts)?
            .expression_with_selectors();

        let fourddenoms = GS.wrap_tree_denoms(graph.denominator(&graph.tree_edges, |_| -1));

        Ok(CFFapprox::Dependent {
            sign: Sign::Positive,
            t_arg: IntegrandExpr {
                integrand: cff_sum * fourddenoms,
            },
        })
    }

    pub(crate) fn root(
        graph: &Graph,
        cuts: &CutSet,
        settings: &UVgenerationSettings,
    ) -> Result<CFFapprox> {
        Self::dependent(graph, &graph.empty_subgraph::<SuBitGraph>(), cuts, settings)
    }
}
impl Approximation {
    pub(crate) fn dot(&self, graph: &Graph) -> String {
        graph.dot_lmb_of(&self.subgraph, &self.lmb)
    }

    pub(crate) fn root(
        &mut self,
        graph: &Graph,
        cuts: &CutSet,
        settings: &UVgenerationSettings,
    ) -> Result<()> {
        self.local_3d = CFFapprox::root(graph, cuts, settings)?;
        self.integrated_4d = ApproxOp::Root;
        self.integrated_pole_part = ApproxOp::Root;
        self.simple_approx = Some(SimpleApprox::root(graph.as_ref().empty_subgraph()));
        self.final_integrand = Some(self.final_integrand(graph, cuts, settings)?);
        Ok(())
    }

    pub(crate) fn new<G, E, V, H>(
        spinney: InternalSubGraph,
        graph: &G,
        lmb: &LoopMomentumBasis,
    ) -> Approximation
    where
        G: UltravioletGraph + AsRef<HedgeGraph<E, V, H>>,
    {
        let lmb = graph.compatible_sub_lmb(&spinney, graph.dummy_less_full_crown(&spinney), lmb);
        // println!("//lmb for spinney \n{}", graph.dot_lmb(&spinney, &lmb));
        Approximation {
            dod: graph.dod(&spinney),
            subgraph: spinney,
            lmb,
            final_integrand: None,
            simple_approx: None,
            local_3d: CFFapprox::NotComputed,
            integrated_4d: ApproxOp::NotComputed,
            integrated_pole_part: ApproxOp::NotComputed,
        }
    }

    pub(crate) fn integrated_4d(
        &self,
        dependent: &Self,
        vakint: (&Vakint, &vakint::VakintSettings),
        uv_graph: &Graph,
        topo_order: usize,
        settings: &UVgenerationSettings,
        pole_part: bool,
    ) -> Result<ApproxOp> {
        let graph = uv_graph.as_ref();
        let reduced = self.subgraph.subtract(&dependent.subgraph);

        let dep = if pole_part {
            debug!(
                "Computing Integrated pole part of {}",
                self.simple_approx
                    .as_ref()
                    .unwrap()
                    .expr(&graph.full_filter())
            );
            dependent.integrated_pole_part.expr()
        } else {
            debug!(
                "Computing Integrated part of {}",
                self.simple_approx
                    .as_ref()
                    .unwrap()
                    .expr(&graph.full_filter())
            );
            dependent.integrated_4d.expr()
        };
        let Some((inner_t, sign)) = dep else {
            return Ok(ApproxOp::NotComputed);
        };

        //(int + inner)*red
        // K = T - Tb
        // K(A*K(B))= T(A*T(B)) + Tb(A*Tb(B))
        // Tb(A*Tb(B))
        //
        // T(A*T4d(B))
        //Tb(A*(T(B)+Tb(B))=Tb(A*Tb(B))
        // (1 + G_inner)*G_red
        let mut t_arg = uv_graph
            .numerator(&reduced, &dependent.subgraph)
            .to_d_dim(GS.dim)
            .get_single_atom()
            .unwrap();

        if pole_part {
            debug!(t_arg = %t_arg.log_print(),"T arg for pole part 4d CT");
        } else {
            debug!(t_arg = %t_arg.log_print(),"T arg for integrated 4d CT");
        }
        t_arg = t_arg.simplify_metrics().simplify_gamma() / uv_graph.denominator(&reduced, |_| 1);
        if pole_part {
            debug!(t_arg = %t_arg.log_print(),"T arg  gamma simplified for pole part 4d CT");
        } else {
            debug!(t_arg = %t_arg.log_print(),"T arg gamma simplified for integrated 4d CT");
        }

        t_arg = t_arg
            .replace(GS.dim)
            .max_level(0)
            .with(Atom::var(GS.dim_epsilon) * (-2) + 4);

        let n_loops = uv_graph.n_loops(&graph.full_filter());

        let mut atomarg = t_arg * inner_t;

        debug!(atomarg = %atomarg.log_print(),"t_arg * inner_t for 4d CT");

        // only apply replacements for edges in the reduced graph
        let mom_reps = graph.uv_wrapped_replacement(&reduced, &self.lmb, &[W_.x___]);

        // println!("Reps:");
        // for r in &mom_reps {
        //     println!("{r}");
        // }

        // println!(
        //     "Expand-prerep {} with dod={} in {:?}",
        //     atomarg, self.dod, self.lmb.ext_edges
        // );

        // rewrite the inner_t as well
        atomarg = atomarg.replace_multiple(&mom_reps);

        // println!(
        //     "Expand {} with dod={} in {:?}",
        //     atomarg, self.dod, self.lmb.ext_edges
        // );
        for e in &self.lmb.ext_edges {
            atomarg = atomarg
                .replace(function!(GS.emr_mom, usize::from(*e) as i64, W_.x___))
                .with(function!(GS.emr_mom, usize::from(*e) as i64, W_.x___) * GS.rescale);
        }
        // TODO: only enable soft CT if doing OS renormalization
        let soft_ct =
            graph.full_crown(&self.subgraph).n_included() == 2 && self.dod > 0 && settings.softct;

        let mut masses = AHashSet::new();
        masses.insert(Atom::var(GS.m_uv));
        // scale all masses, including UV masses from subgraphs

        for (p, _, e) in graph.iter_edges_of(&self.subgraph) {
            if p.is_paired() {
                let e_mass = e.data.mass_atom();
                masses.insert(e_mass);
            }
        }

        if !soft_ct {
            for m in &masses {
                let rescaled = m.clone() * GS.rescale;
                atomarg = atomarg.replace(m.clone()).with(rescaled);
            }

            // expand the propagator around a propagator with a UV mass
            atomarg = atomarg
                .replace(parse!("den(n_,q_,mass_,prop_)"))
                .with(parse!(
                    "den(n_,q_,mass_ + mUV^2 - t^2*mUV^2, prop_- mUV^2 + t^2*mUV^2)"
                ));
        }

        debug!(atomarg = %atomarg.log_print(),"t_arg * inner_t after rescaling masses for 4d CT");
        // den(..) tags a propagator, its first derivative is 1 and the rest is 0
        let mut a = atomarg
            .series(GS.rescale, Atom::Zero, self.dod.into(), true)
            .unwrap()
            .to_atom()
            .replace(parse!("der(0,0,0,1, den(y__))"))
            .with(Atom::num(1))
            .replace(parse!("der(x__, den(y__))"))
            .with(Atom::num(0));

        debug!(a = %a.log_print(),"Series expanded for 4d CT");

        if soft_ct {
            let coeffs = a.coefficient_list::<u8>(&[Atom::var(GS.rescale)]);
            let mut b = Atom::Zero;
            let dod_pow = Atom::var(GS.rescale).pow(self.dod);
            for (pow, mut i) in coeffs {
                if pow == dod_pow {
                    // set the masses in the t=dod term to 0
                    // UV rearrange the denominators
                    for m in &masses {
                        i = i.replace(m.clone()).with(Atom::Zero);
                    }

                    i = i
                        .replace(parse!("den(n_,q_,mass_,prop_)"))
                        .with(parse!("den(n_,q_,mUV^2,prop_-mUV^2)"));
                }

                b += i;
            }

            a = b;
        } else {
            a = a.replace(GS.rescale).with(Atom::num(1));
        }

        let mut integrand_vakint = to_vakint_integrand(
            &a,
            graph,
            &reduced,
            &dependent.subgraph,
            &settings.vakint,
            true,
        );

        for t in &integrand_vakint.0 {
            debug!(integral = %t.integral.log_print(),numerator = %t.numerator.log_print(),"Vakint term as input");
        }
        debug!(settings = ?&vakint.1,"Vakint args");

        // let mut res = vakint
        //     .0
        //     .evaluate(&vakint.1, integrand_vakint.as_view())
        //     .unwrap();

        integrand_vakint.canonicalize(&vakint.1, &vakint.0.topologies, false)?;
        for t in &integrand_vakint.0 {
            debug!(integral = %t.integral.log_print(),numerator = %t.numerator.log_print(),"Vakint term after canonicalization");
        }
        integrand_vakint.tensor_reduce(&vakint.0, &vakint.1)?;
        for t in &integrand_vakint.0 {
            debug!(integral = %t.integral.log_print(),numerator = %t.numerator.log_print(),"Vakint term after tensor reduction");
        }
        integrand_vakint.evaluate_integral(&vakint.0, &vakint.1)?;
        for t in &integrand_vakint.0 {
            debug!(integral = %t.integral.log_print(),numerator = %t.numerator.log_print(),"Vakint term after evaluation");
        }

        let mut res: Atom = integrand_vakint.into();

        debug!(res = %res.expand().log_print(),"Raw post vakint ");

        res = res
            .replace(parse_lit!(vakint::cl2))
            .with(parse_lit!(cl2))
            .replace(parse_lit!(vakint::sqrt3))
            .with(parse_lit!(sqrt(3)));

        let vk_metric = vakint_symbol!("g");
        let mink = Minkowski {}.new_rep(GS.dim);
        // apply metric
        res = res
            .replace(vakint::symbols::S.p.f(&[W_.i_, W_.j_]))
            .when(W_.j_.filter_single(|r| r.is_integer()))
            .with(vakint::symbols::S.p.f(&[
                Atom::var(W_.i_),
                mink.to_symbolic([GS.uvaind.f(&[Atom::num(topo_order), Atom::var(W_.j_)])]),
            ]))
            .replace(
                vakint::symbols::S
                    .p
                    .f(&[Atom::var(W_.i_), vakint::symbols::S.dot_dummy_ind(W_.j_)]),
            )
            .when(W_.j_.filter_single(|r| r.is_integer()))
            .with(vakint::symbols::S.p.f(&[
                Atom::var(W_.i_),
                mink.to_symbolic([GS.uvaind.f(&[Atom::num(topo_order), Atom::var(W_.j_)])]),
            ]))
            .replace(vakint::symbols::S.p.f(&[W_.x__]))
            .with(GS.emr_mom.f(&[W_.x__]));
        res = res
            .replace(function!(vk_metric, W_.x_, W_.y_) * function!(GS.emr_mom, W_.x___, W_.x_))
            .with(function!(GS.emr_mom, W_.x___, W_.y_))
            .replace(function!(
                vk_metric,
                vakint::symbols::S.dot_dummy_ind(W_.x_),
                W_.y_
            ))
            .when(W_.x_.filter_single(|r| r.is_integer()))
            .with(function!(
                vk_metric,
                mink.to_symbolic([GS.uvaind.f(&[Atom::num(topo_order), Atom::var(W_.x_)])]),
                W_.y_
            ))
            .replace(function!(
                vk_metric,
                W_.x_,
                vakint::symbols::S.dot_dummy_ind(W_.y_)
            ))
            .when(W_.y_.filter_single(|r| r.is_integer()))
            .with(function!(
                vk_metric,
                mink.to_symbolic([GS.uvaind.f(&[Atom::num(topo_order), Atom::var(W_.y_)])]),
                W_.x_
            ))
            .replace(function!(vk_metric, W_.x_, W_.y_))
            .with(function!(ETS.metric, W_.x_, W_.y_));

        res = res.replace(vakint::symbols::S.cmplx_i).with(Atom::i());

        res = res
            .simplify_metrics()
            .replace(GS.dim)
            .max_level(0)
            .with(Atom::var(GS.dim_epsilon) * (-2) + 4);

        debug!(res = %res.expand().log_print(),"Replaced post vakint ");
        let series = res
            .series(
                GS.dim_epsilon,
                Atom::Zero,
                (n_loops as i64 + 1).into(),
                true,
            )
            .unwrap();

        debug!(series = %series.to_atom().log_print(),"Series ");

        let mut pole_stripped = Atom::Zero;

        for (power, p) in series.terms() {
            // println!("Power: {}", power);
            // println!("Coeff: {}", p.printer(LOGPRINTOPTS));
            if pole_part {
                if power < 0 {
                    pole_stripped += p * Atom::var(GS.dim_epsilon).pow(power);
                }
            } else {
                if power >= 0 {
                    pole_stripped += p * Atom::var(GS.dim_epsilon).pow(power);
                }
            }
        }

        res = pole_stripped;

        if !pole_part {
            // multiply the results with a vacuum triangle that integrates to 1
            // 1/(k^2 - m_UV^2)^3 = -i / (4 pi)^2 * 1/2 * 1/mUV^2
            // name the mUV mass mUVi as this one should not be expanded
            for l in &self.lmb.loop_edges {
                if !reduced.includes(&graph[l].1) {
                    continue;
                }

                //TODO: Add orientation localisation prefactor (Sum of valid orientation thetas)/(number of valid orientations)

                res /= parse!("(-1i / (4 𝜋)^2 * 1/2 * 1/mUVI^2)");

                let mink: Slot<Minkowski, Aind> = Minkowski {}.new_rep(4).slot(Aind::new_dummy());

                // multiply CFF triangle
                res *= Atom::num((3, 16))
                    / (GS.emr_vec_index(*l, mink.to_atom()) * GS.emr_vec_index(*l, mink.to_atom())
                        + GS.m_uv_int * GS.m_uv_int)
                        .pow((5, 2));
            }
        }

        debug!(pole_part = %pole_part,res = %res.log_print(),"Final integrated 4d CT");

        if res
            .replace(GS.dim)
            .max_level(0)
            .match_iter()
            .next()
            .is_some()
        {
            panic!(
                "The t_arg should not contain dim after expansion, found {}",
                res
            );
        }

        // println!("\nIntegrated CT:\n{}\n", res);

        Ok(ApproxOp::Dependent {
            t_arg: IntegrandExpr { integrand: res },
            sign: -sign,
            subgraph: reduced,
        })
    }

    pub(crate) fn compute_integrated(
        &mut self,
        graph: &Graph,
        vakint: (&Vakint, &vakint::VakintSettings),
        dependent: &Self,
        topo_order: usize,
        settings: &UVgenerationSettings,
    ) -> Result<()> {
        self.integrated_4d =
            self.integrated_4d(dependent, vakint, graph, topo_order, settings, false)?;
        self.integrated_pole_part =
            self.integrated_4d(dependent, vakint, graph, topo_order, settings, true)?;
        Ok(())
    }

    /// Computes the 3d approximation of the UV
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn compute(
        &mut self,
        graph: &Graph,
        cuts: &CutSet,
        dependent: &Self,
        settings: &UVgenerationSettings,
    ) -> Result<()> {
        let Some((cff, sign)) = dependent.local_3d.expr() else {
            panic!("Should have computed the dependent cff");
        };
        let (t4, _) = if let ApproxOp::Root = dependent.integrated_4d {
            (Atom::num(0), Sign::Positive)
        } else {
            dependent
                .integrated_4d
                .expr()
                .unwrap_or((Atom::num(0), Sign::Positive))
        };

        let CFFapprox::Dependent { t_arg, .. } =
            CFFapprox::dependent(graph, &dependent.subgraph.filter, cuts, settings)?
        else {
            unreachable!()
        };

        let mut sum_3d = Atom::Zero;

        sum_3d += self.local_3d(dependent, graph, cff, settings);

        let finite = t4
            .series(vakint_symbol!("ε"), Atom::Zero, 0.into(), true)
            .unwrap()
            .coefficient((0, 1).into())
            .replace(GS.m_uv_int)
            .with(GS.m_uv)
            .map_mink_dim(4);

        debug!(
            "Integrated 4d finite part: {:#}",
            finite.printer(LOGPRINTOPTS)
        );

        // TODO: multiply by the number of orientations or only apply the counterterm to
        // one orientation
        // subtract integrated CT
        sum_3d -= self.local_3d(dependent, graph, finite * t_arg.integrand, settings);

        self.local_3d = CFFapprox::Dependent {
            sign: -sign,
            t_arg: IntegrandExpr {
                integrand: sum_3d, //* 2, //* Atom::i(),
            },
        };

        self.final_integrand = Some(self.final_integrand(graph, cuts, settings)?);
        Ok(())
    }

    #[instrument(skip_all)]
    pub(crate) fn local_3d<E: UVE, V, H, G: UltravioletGraph + AsRef<HedgeGraph<E, V, H>>>(
        &self,
        dependent: &Self,
        uv_graph: &G,
        mut cff: Atom,
        settings: &UVgenerationSettings,
    ) -> Atom {
        let graph = uv_graph.as_ref();
        let reduced = self.subgraph.subtract(&dependent.subgraph);

        // println!("CFF: {}", cff);

        // add data for OSE computation and add an explicit sqrt
        for (p, ei, e) in graph.iter_edges_of(&self.subgraph) {
            let eid = usize::from(ei) as i64;
            if p.is_paired() {
                // set energies from inner_t on-shell
                cff = cff.replace(function!(GS.energy, eid)).with(GS.ose(ei));

                let e_mass = e.data.mass_atom();
                cff = cff.replace(GS.ose(ei)).with(GS.ose_full(
                    ei,
                    e_mass,
                    None,
                    settings.inner_products,
                ));
            }
        }

        let mut atomarg = cff
            * uv_graph
                .numerator(&reduced, &dependent.subgraph)
                .get_single_atom()
                .unwrap();

        // println!(
        //     "Expand-prerep {} with dod={} in {:?}",
        //     atomarg, self.dod, self.lmb.ext_edges
        // );

        // split numerator momenta into OSEs and spatial parts
        let mut reps = Vec::new();
        for (p, eid, e) in graph.iter_edges_of(&self.subgraph) {
            if p.is_paired() {
                let e_mass = e.data.mass_atom();
                reps.push(GS.split_mom_pattern(eid, e_mass, settings.inner_products));
            }
        }

        // println!("Split reps:");
        // for r in &reps {
        //     println!("{r}");
        // }

        atomarg = atomarg.replace_multiple(&reps);

        // only apply replacements for edges in the reduced graph
        let mom_reps = graph.uv_spatial_wrapped_replacement(&reduced, &self.lmb, &[W_.x___]);

        // println!(
        //     "Mom Reps : for {}",
        //     self.simple_approx
        //         .as_ref()
        //         .unwrap()
        //         .expr(&graph.full_filter())
        // );
        // for r in &mom_reps {
        //     println!("{r}");
        // }

        atomarg = atomarg.replace_multiple(&mom_reps);

        // debug!("Before rescaling loop momenta {}",);

        // rescale the loop momenta in the whole subgraph, including previously expanded cycles
        for e in &self.lmb.loop_edges {
            // println!("Rescale {}", e);
            atomarg = atomarg
                .replace(GS.emr_vec_index(*e, W_.x___))
                .with(GS.emr_vec_index(*e, W_.x___) * GS.rescale);
        }

        // println!(
        //     "Expand {} with dod={} in {:?}",
        //     atomarg, self.dod, self.lmb.ext_edges
        // );

        let soft_ct =
            graph.full_crown(&self.subgraph).n_included() == 2 && self.dod > 0 && settings.softct;

        // (re-)expand OSEs from the subgraph only
        for (_, eid, _) in graph.iter_edges_of(&self.subgraph) {
            let eid = usize::from(eid) as i64;
            if soft_ct {
                info!("DOing soft ct{}", graph.dot(&self.subgraph));
                // TODO: rescale the masses in OSEs
                // TODO: also scale masses in the numerator _only_ for the subgraph
                // expand the OSEs around an OSE with a UV mass
                todo!()
            } else {
                // rescale the whole OSE so that the function itself has no poles during the expansion
                atomarg = atomarg
                    .replace(function!(GS.ose, eid, W_.mom_, W_.mass_, W_.prop_))
                    .with(
                        function!(
                            GS.ose,
                            eid,
                            W_.mom_,
                            GS.m_uv * GS.m_uv,
                            (GS.m_uv * GS.m_uv * GS.rescale * GS.rescale + W_.prop_
                                - GS.m_uv * GS.m_uv)
                                / GS.rescale
                                / GS.rescale
                        ) * GS.rescale
                            * GS.rescale,
                    )
                    .replace(function!(GS.ose, eid, W_.mom_, W_.a___)) //rescale the momenta for the same reason
                    .with_map(move |m| {
                        let mut f = FunctionBuilder::new(GS.ose);
                        f = f.add_arg(eid);
                        f = f.add_arg(
                            (m.get(W_.mom_)
                                .unwrap()
                                .to_atom()
                                .replace(GS.rescale)
                                .with(Atom::num(1) / GS.rescale)
                                * GS.rescale)
                                .expand()
                                .replace(GS.rescale)
                                .with(Atom::Zero),
                        );
                        f = f.add_arg(m.get(W_.a___).unwrap().to_atom());

                        f.finish()
                    });
            }
        }

        // atomarg = atomarg
        //     .replace(function!(MS.dot, GS.rescale * W_.x_, W_.y_))
        //     .repeat()
        //     .with(function!(MS.dot, W_.x_, W_.y_) * GS.rescale);

        atomarg = (atomarg
            * Atom::var(GS.rescale).pow(3 * uv_graph.n_loops(&self.subgraph) as i64))
        .replace(GS.rescale)
        .with(Atom::num(1) / GS.rescale);
        // .replace(Atom::var(GS.rescale).pow(2).sqrt()) //.pow((1, 2)))
        // .with(GS.rescale)
        // .replace((Atom::var(GS.rescale).pow(-2) * W_.a___).pow((1, 2))) //.pow((1, 2)))
        // .repeat()
        // .with((Atom::var(W_.a___)).sqrt() / GS.rescale)
        // .replace((Atom::var(GS.rescale).pow(2) * W_.a___).pow((-1, 2))) //.pow((1, 2)))
        // .repeat()
        // .with((Atom::var(W_.a___)).pow((-1, 2)) / GS.rescale)
        // .replace((Atom::var(GS.rescale).pow(-2) * W_.a___).pow((-1, 2))) //.pow((1, 2)))
        // .repeat()
        // .with((Atom::var(W_.a___)).pow((-1, 2)) * GS.rescale)
        // .replace((Atom::var(GS.rescale).pow(2) * W_.a___).pow((1, 2))) //.pow((1, 2)))
        // .repeat()
        // .with((Atom::var(W_.a___)).pow((1, 2)) * GS.rescale);

        // println!("atomarg:{:>}", atomarg.expand());

        // println!("Series expanding {atomarg}");
        //
        let t = Atom::var(GS.rescale);
        debug!(
            atom = %atomarg,
            "Series expanding {} up to dod {}:{}",
            self.simple_approx
                .as_ref()
                .unwrap()
                .expr(&graph.full_filter()),
            self.dod,
            // orientations
            //     .first()
            //     .unwrap()
            //     .select(&
            atomarg.collect_multiple::<i8>(&[t], None, None)
                .replace(function!(GS.theta, W_.a_))
                .with(Atom::one())
                // .replace(GS.m_uv)
                // .with(Atom::Zero)
                // .replace(parse_lit!(UFO::mass_scalar_1))
                // .with(Atom::Zero)
                .replace((Atom::var(GS.rescale).pow(-2) * W_.a_).pow((1, 2))) //.pow((1, 2)))
                .repeat()
                .with((Atom::var(W_.a_)).sqrt() / GS.rescale)
                .collect_factors()
                .collect_num()
                .typst_string() // printer(LOGPRINTOPTS)
        );

        let a = atomarg
            .series(GS.rescale, Atom::Zero, 0.into(), true)
            .unwrap();

        // debug!("Series: {}", a);

        let mut a = a
            .to_atom()
            .replace(parse!("der(0,0,0,1, OSE(y__))"))
            .with(Atom::num(1))
            .replace(parse!("der(x__, OSE(y__))"))
            .with(Atom::num(0));

        if soft_ct {
            let coeffs = a.coefficient_list::<u8>(&[Atom::var(GS.rescale)]);
            let mut b = Atom::Zero;
            let dod_pow = Atom::var(GS.rescale).pow(self.dod);
            for (pow, mut i) in coeffs {
                if pow == dod_pow {
                    // FIXME: how to do this _only_ for the subgraph masses? the numerator is still
                    // only that of the subgraph
                    // set the masses in the t=dod term to 0
                    // UV rearrange the denominators
                    /*for m in &masses {
                        i = i.replace(m.clone()).with(Atom::Zero);
                    }*/

                    i = i
                        .replace(parse!("OSE(n_,q_,mass_,prop_, x___)"))
                        .with(parse!("OSE(n_,q_,mUV^2,prop_+mUV^2,x___)"));
                }

                b += i;
            }

            a = b;
        } else {
            a = a.replace(GS.rescale).with(Atom::num(1));
        }

        // println!("Expanded: {:>}", a.expand());

        a
    }

    #[instrument(skip_all)]
    pub(crate) fn final_integrand(
        &self,
        graph: &Graph,
        cutset: &CutSet,
        settings: &UVgenerationSettings,
    ) -> Result<Atom> {
        let (t, s) = self
            .local_3d
            .expr()
            .ok_or(eyre!("Local3d not yet computed"))?;
        let (t_int, _) = if let ApproxOp::Root = self.integrated_4d {
            (Atom::num(0), Sign::Positive)
        } else {
            self.integrated_4d
                .expr()
                .unwrap_or((Atom::num(0), Sign::Positive))
        };

        let CFFapprox::Dependent { t_arg, .. } =
            CFFapprox::dependent(graph, &self.subgraph.filter, cutset, settings)?
        else {
            unreachable!()
        };

        let finite = t_int
            .series(vakint_symbol!("ε"), Atom::Zero, 0.into(), true)
            .unwrap()
            .coefficient((0, 1).into())
            .replace(GS.m_uv_int)
            .with(GS.m_uv)
            .map_mink_dim(4)
            .replace(function!(symbol!("vakint::g"), W_.a__))
            .with(function!(symbol!("spenso::g"), W_.a__));

        debug!(
            "Integrated 4d finite part: {:#}",
            finite.printer(LOGPRINTOPTS)
        );
        let reduced = graph.full_filter().subtract(self.subgraph.included());

        // let concrete_red = graph.as_ref().concretize(&reduced).map(
        //     |_, _, v| v.clone(),
        //     |_, _, _, e, d| d.map(Clone::clone),
        //     |h, s| NumHedgeData::de,
        // );

        // graph.as_ref()..contract_subgraph(&reduced, node_data_merge);

        let mut cff = s * t - s * finite * t_arg.integrand;

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
            .numerator(&reduced, &self.subgraph.included())
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

        resnum = resnum.replace_multiple(&reps).replace(GS.dim).with(4);
        resnum *= cff;
        resnum = resnum.wrap_color(GS.color_wrap);

        debug!(
            "Integrand before parsing for {} for dod{}:{}",
            self.simple_approx
                .as_ref()
                .unwrap()
                .expr(&graph.full_filter()),
            self.dod,
            // orientations
            //     .first()
            //     .unwrap()
            //     .select(&
            resnum
                .replace(function!(GS.theta, W_.a_))
                .with(Atom::one())
                .replace(GS.m_uv)
                .with(Atom::Zero)
                .replace(parse_lit!(UFO::mass_scalar_1))
                .with(Atom::Zero)
                .collect_factors()
                .collect_num()
                .log_print() // printer(LOGPRINTOPTS)
        );

        // debug!("final_cff {res:>}");
        Ok(resnum.replace_multiple(&reps))
    }

    // pub(crate) fn simple_expr(
    //     &self,
    //     graph: &UVGraph,
    //     amplitude: &InternalSubGraph,
    // ) -> Option<SerializableAtom> {
    //     let simple_approx = self.simple_approx.as_ref()?;

    //     Some((simple_approx.sign * simple_approx.expr(&amplitude.filter)).into())
    // }
}

impl ApproxOp {
    pub(crate) fn sign(&self) -> Option<Sign> {
        match self {
            ApproxOp::NotComputed => None,
            ApproxOp::Union { sign, .. } => Some(*sign),
            ApproxOp::Dependent { sign, .. } => Some(*sign),
            ApproxOp::Root => Some(Sign::Positive),
        }
    }

    pub(crate) fn expr(&self) -> Option<(Atom, Sign)> {
        match self {
            ApproxOp::NotComputed => None,
            ApproxOp::Union { t_args, sign, .. } => {
                //Never gets hit now
                let mut mul = Atom::num(1);
                for t in t_args {
                    mul *= &t.integrand;
                }
                Some((mul, *sign))
            }
            ApproxOp::Dependent { t_arg, sign, .. } => Some((t_arg.integrand.clone(), *sign)),
            ApproxOp::Root => Some((Atom::num(1), Sign::Positive)),
        }
    }

    pub(crate) fn union(dependent: &[&Approximation]) -> Option<Self> {
        let mut t_args = vec![];
        let mut subgraphs = vec![];

        let mut final_sign = Sign::Positive;
        for d in dependent {
            match &d.integrated_4d {
                ApproxOp::Dependent {
                    t_arg,
                    sign,
                    subgraph,
                } => {
                    t_args.push(t_arg.clone());
                    final_sign = final_sign * *sign;
                    subgraphs.push(subgraph.clone())
                }
                _ => return None,
            }
        }

        Some(Self::Union {
            t_args,
            sign: final_sign,
            subgraphs,
        })
    }

    pub(crate) fn is_computed(&self) -> bool {
        !matches!(self, ApproxOp::NotComputed)
    }
}
