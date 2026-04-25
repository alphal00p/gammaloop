use ahash::{HashSet, HashSetExt};
use color_eyre::Result;
use eyre::eyre;
use idenso::{gamma::GammaSimplifier, metric::MetricSimplifier};
use linnet::half_edge::{
    HedgeGraph, NodeIndex,
    builder::HedgeGraphBuilder,
    involution::HedgePair,
    subgraph::{Inclusion, ModifySubSet, SuBitGraph, SubGraphLike, SubSetLike},
};
use spenso::{
    network::{library::symbolic::ETS, parsing::SPENSO_TAG},
    structure::representation::{Minkowski, RepName},
};
use symbolica::{
    atom::{Atom, AtomCore, Symbol},
    function,
    id::{MatchSettings, Replacement},
    parse, parse_lit,
    solve::SolveError,
};
use tracing::{debug, instrument};
use vakint::{Vakint, VakintExpression, vakint_symbol};

use crate::{
    graph::{Graph, LMBext, LoopMomentumBasis},
    utils::{
        GS, W_,
        symbolica_ext::{CallSymbol, LogPrint},
    },
    uv::{
        ApproximationType, UltravioletGraph,
        approx::{ApproximationKernel, ForestNodeLike, UVCtx},
        settings::VakintSettings,
        uv_graph::UVE,
    },
};

pub struct Integrated<'a> {
    pub vakint: &'a Vakint,
    pub vakint_settings: &'a vakint::VakintSettings,
}

impl Graph {
    pub(crate) fn uv_rescaled(
        &self,
        replacement_subgraph: &SuBitGraph,
        n_loops: usize,
        lmb: &LoopMomentumBasis,
        atom: &Atom,
    ) -> Atom {
        // only apply replacements for edges in the reduced graph
        let mom_reps = self.uv_wrapped_replacement(replacement_subgraph, lmb, &[W_.x___]);
        let mut atomarg = atom.replace_multiple(&mom_reps);

        // rescale the loop momenta in the whole subgraph, including previously expanded cycles
        for edge in &lmb.loop_edges {
            atomarg = atomarg
                .replace(GS.emr_mom(*edge, W_.x___))
                .with(GS.emr_mom(*edge, W_.x___) / GS.rescale);
        }

        let tsquare = Atom::var(GS.rescale).pow(2);

        debug!(res = %atomarg.log_print(None),"Rescaled momenta expanded");
        atomarg = atomarg
            .replace(GS.den(W_.a_, W_.mom_, W_.mass_, W_.prop_))
            .with(
                GS.den(
                    W_.a_,
                    W_.mom_,
                    Atom::var(W_.mass_) + Atom::var(GS.m_uv).pow(2),
                    Atom::var(W_.prop_) * &tsquare + Atom::var(GS.m_uv).pow(2) * &tsquare
                        - (Atom::var(GS.m_uv)).pow(2),
                ) / &tsquare,
            )
            .replace(function!(GS.den, W_.a_, W_.mom_, W_.a___))
            .with_map(move |m| {
                let mut f = symbolica::atom::FunctionBuilder::new(GS.den);
                f = f.add_arg(m.get(W_.a_).unwrap().to_atom());
                f = f.add_arg(
                    (m.get(W_.mom_).unwrap().to_atom() * GS.rescale)
                        .expand()
                        .replace(GS.rescale)
                        .with(Atom::Zero),
                );
                f = f.add_arg(m.get(W_.a___).unwrap().to_atom());
                f.finish()
            });

        atomarg *= Atom::var(GS.rescale).pow(-4 * n_loops as i64);
        atomarg
    }
}

impl Integrated<'_> {
    pub fn new<'a>(
        vakint: &'a Vakint,
        vakint_settings: &'a vakint::VakintSettings,
    ) -> Integrated<'a> {
        Integrated {
            vakint,
            vakint_settings,
        }
    }

    pub(crate) fn start<S: super::ForestNodeLike>(
        &self,
        ctx: &UVCtx<'_>,
        current: &S,
        given: &S,
        integrand: &Atom,
    ) -> Result<Atom> {
        let reduced = current.reduced_subgraph(given);
        let graph = ctx.graph;
        let settings = &ctx.settings.uv;
        let mut t_arg = ctx
            .graph
            .numerator(&reduced, given.subgraph())
            .to_d_dim(GS.dim)
            .get_single_atom()
            .unwrap();

        debug!(t_arg = %t_arg.log_print(None),pole_part=%settings.pole_part,"T arg without denoms");
        t_arg = t_arg.simplify_metrics().simplify_gamma() / graph.denominator(&reduced, |_| 1);
        debug!(t_arg = %t_arg.log_print(None),pole_part=%settings.pole_part,"T arg gamma simplified for integrated 4d CT");

        t_arg = t_arg
            .replace(GS.dim)
            .max_level(0)
            .with(Atom::var(GS.dim_epsilon) * (-2) + 4);

        Ok(t_arg * integrand)
    }

    pub(crate) fn t<S: super::ForestNodeLike>(
        &self,
        ctx: &UVCtx<'_>,
        current: &S,
        given: &S,
        integrand: &Atom,
    ) -> Result<Atom> {
        let graph = ctx.graph;
        let reduced = current.reduced_subgraph(given);
        let n_loops = graph.n_loops(current.subgraph()) - graph.n_loops(given.subgraph());

        let atomarg = graph.uv_rescaled(&reduced, n_loops, current.lmb(), integrand);

        debug!(res = %atomarg.log_print(None),n_loops=%n_loops,"Rescaled expanded");
        let a = atomarg
            .series(GS.rescale, Atom::Zero, 1.into(), true)
            .unwrap();

        let mut a = a.to_atom();

        debug!(res = %a.log_print(None),res_raw = %a.to_plain_string(),"Series expanded");
        a = a.replace(GS.rescale).with(Atom::num(1));
        debug!(res = %a.log_print(None),"Series expanded");

        Ok(a)
    }

    pub(crate) fn integrate_and_truncate<S: super::ForestNodeLike>(
        &self,
        ctx: &UVCtx<'_>,
        current: &S,
        given: &S,
        integrand: &Atom,
    ) -> Result<Atom> {
        let graph = ctx.graph;
        let reduced = current.reduced_subgraph(given);
        let settings = ctx.settings;
        let mut integrand_vakint = to_vakint_integrand(
            integrand,
            graph,
            &reduced,
            given.subgraph(),
            &settings.uv.vakint,
            true,
        );

        for t in &integrand_vakint.0 {
            debug!(integral = %t.integral.log_print(None),numerator = %t.numerator.log_print(None),"Vakint term as input");
        }
        debug!(settings = ?&self.vakint_settings,"Vakint args");

        // let mut res = vakint
        //     .0
        //     .evaluate(&vakint.1, integrand_vakint.as_view())
        //     .unwrap();

        integrand_vakint.canonicalize(self.vakint_settings, &self.vakint.topologies, false)?;
        for t in &integrand_vakint.0 {
            debug!(integral = %t.integral.log_print(None),integral_raw = %t.integral.to_plain_string(), numerator_raw = %t.numerator.to_plain_string(), numerator = %t.numerator.log_print(None),"Vakint term after canonicalization");
        }
        integrand_vakint.tensor_reduce(self.vakint, self.vakint_settings)?;
        for t in &integrand_vakint.0 {
            debug!(integral = %t.integral.log_print(None),numerator = %t.numerator.log_print(None),"Vakint term after tensor reduction");
        }
        integrand_vakint.evaluate_integral(self.vakint, self.vakint_settings)?;
        for t in &integrand_vakint.0 {
            debug!(integral = %t.integral.log_print(None),numerator = %t.numerator.log_print(None),"Vakint term after evaluation");
        }

        let mut res: Atom = integrand_vakint.into();

        debug!(res = %res.expand().log_print(None),"Raw post vakint ");

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
            .with(
                vakint::symbols::S.p.f(&[
                    Atom::var(W_.i_),
                    mink.to_symbolic([GS
                        .uvaind
                        .f(&[Atom::num(current.topo_order()), Atom::var(W_.j_)])]),
                ]),
            )
            .replace(
                vakint::symbols::S
                    .p
                    .f(&[Atom::var(W_.i_), vakint::symbols::S.dot_dummy_ind(W_.j_)]),
            )
            .when(W_.j_.filter_single(|r| r.is_integer()))
            .with(
                vakint::symbols::S.p.f(&[
                    Atom::var(W_.i_),
                    mink.to_symbolic([GS
                        .uvaind
                        .f(&[Atom::num(current.topo_order()), Atom::var(W_.j_)])]),
                ]),
            )
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
                mink.to_symbolic([GS
                    .uvaind
                    .f(&[Atom::num(current.topo_order()), Atom::var(W_.x_)])]),
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
                mink.to_symbolic([GS
                    .uvaind
                    .f(&[Atom::num(current.topo_order()), Atom::var(W_.y_)])]),
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

        debug!(res = %res.expand().log_print(None),"Replaced post vakint ");

        let n_loops = graph.n_loops(&graph.full_filter());
        let series = res
            .series(
                GS.dim_epsilon,
                Atom::Zero,
                (n_loops as i64 + 1).into(),
                true,
            )
            .unwrap();

        debug!(series = %series.to_atom().log_print(None),"Series ");

        let mut pole_stripped = Atom::Zero;

        for (power, p) in series.terms() {
            // println!("Power: {}", power);
            // println!("Coeff: {}", p.printer(LOGPRINTOPTS));
            if settings.uv.pole_part {
                if power < 0 {
                    pole_stripped += p * Atom::var(GS.dim_epsilon).pow(power);
                }
            } else if power >= 0 {
                pole_stripped += p * Atom::var(GS.dim_epsilon).pow(power);
            }
        }

        res = pole_stripped;

        if !settings.uv.pole_part {
            // Multiply by the localized normalized integral \int \vec{k} 1 / (|\vec{k}|^2 + mUV^2)^2, which integrates to \pi^2/ mUV
            let pi_atom = (Symbol::PI).to_atom();
            let mut normalization_term_integral = (pi_atom.pow(2)) / GS.m_uv_int;
            // However, gammaloop adds a factro 1/(2*pi)^3 per loop, and this integrated CT will be subject to it, so we must undo it.
            normalization_term_integral /= (Atom::from(2) * pi_atom).pow(3);

            // We need to correct the Wick rotation `i` per loop
            // TODO: Understand this better: this is *not* part of the normalization really, but probably related to the fact that our
            // UV CT is using minkowski denominators and not euclidean ones.
            normalization_term_integral /= Atom::i();

            for l in &current.lmb().loop_edges {
                if !reduced.includes(&graph[l].1) {
                    continue;
                }

                //TODO: Add orientation localisation prefactor (Sum of valid orientation thetas)/(number of valid orientations)
                res /= normalization_term_integral.as_view();

                let spatial_norm_sq = integrated_triangle_spatial_norm_sq(*l);

                // Per-orientation CFF localizer of the normalized cubic tadpole.
                let denominator = spatial_norm_sq + GS.m_uv_int * GS.m_uv_int;
                res /= denominator.as_view() * denominator.as_view();
            }
        }

        debug!(pole_part = %settings.uv.pole_part,res = %res.log_print(None),"Final integrated 4d CT");

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
        Ok(res)
    }
}

fn integrated_triangle_spatial_norm_sq(
    loop_edge: linnet::half_edge::involution::EdgeIndex,
) -> Atom {
    GS.emr_mom(loop_edge, GS.cind(1)).pow(2)
        + GS.emr_mom(loop_edge, GS.cind(2)).pow(2)
        + GS.emr_mom(loop_edge, GS.cind(3)).pow(2)
}

impl ApproximationKernel<UVCtx<'_>> for Integrated<'_> {
    #[instrument(skip_all)]
    fn kernel<S: ForestNodeLike>(
        &self,
        ctx: &UVCtx<'_>,
        current: &S,
        given: &S,
        integrand: &Atom,
    ) -> Result<Atom> {
        match current.renormalization_scheme() {
            ApproximationType::MUV => self.integrate_and_truncate(
                ctx,
                current,
                given,
                &self.t(
                    ctx,
                    current,
                    given,
                    &self.start(ctx, current, given, integrand)?,
                )?,
            ),
            ApproximationType::IR => Err(eyre!("Not yet implemented IR")),
            ApproximationType::VacuumLimit => Ok(Atom::Zero),
            ApproximationType::OS => Err(eyre!("Not yet implemented OS")),
            ApproximationType::Unsubtracted => {
                panic!("should have been kept out of the wood");
            }
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
        integrand = %integrand.log_print(None),
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

    for t in a.0.iter_mut() {
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

                if let Some(m) = &mass
                    && m != &d.data.mass
                {
                    return false;
                } else {
                    mass = Some(d.data.mass.clone());
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

        if !nodes_to_merge.is_empty() {
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

#[cfg(test)]
mod tests {
    use linnet::half_edge::involution::EdgeIndex;

    use super::*;

    #[test]
    fn integrated_triangle_norm_is_euclidean() {
        let edge = EdgeIndex(7);
        let euclidean_norm = integrated_triangle_spatial_norm_sq(edge);
        let minkowski_norm = Minkowski {}
            .new_rep(4)
            .inner_product(GS.emr_vec(edge), GS.emr_vec(edge));

        assert_eq!(
            euclidean_norm,
            GS.emr_mom(edge, GS.cind(1)).pow(2)
                + GS.emr_mom(edge, GS.cind(2)).pow(2)
                + GS.emr_mom(edge, GS.cind(3)).pow(2)
        );
        assert_ne!(euclidean_norm, minkowski_norm);
    }
}
