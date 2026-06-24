use ahash::{HashSet, HashSetExt};
use color_eyre::Result;
use eyre::eyre;
use gammaloop_tracing_filter::{LogMessage, debug_instrument};
use idenso::{
    color::ColorSimplifier,
    dirac::GammaSimplifier,
    representations::Bispinor,
    shorthands::{
        UndoShorthands,
        chain::Chain,
        metric::MetricSimplifier,
        schoonschip::{Schoonschip, SchoonschipSettings},
    },
};

use linnet::half_edge::{
    HedgeGraph, NodeIndex,
    builder::HedgeGraphBuilder,
    involution::HedgePair,
    subgraph::{ModifySubSet, SuBitGraph, SubGraphLike, SubSetLike},
};
use spenso::{
    network::{library::symbolic::ETS, tags::SPENSO_TAG},
    shadowing::TensorCollectExt,
    structure::representation::{Minkowski, RepName},
};
use symbolica::{atom::AtomCore, prelude::*};
use vakint::{Vakint, VakintExpression, vakint_symbol};

use crate::{
    debug_tags,
    graph::{Graph, LMBext, LoopMomentumBasis},
    numerator::aind::Aind,
    utils::{GS, W_, symbolica_ext::CallSymbol},
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

        debug_tags!(#uv, #integrated, #inspect;
            log.res = atomarg,
            "Rescaled momenta expanded"
        );
        atomarg = atomarg
            .replace(GS.den(W_.a_, W_.mom_, W_.mass_, W_.prop_))
            .with(
                GS.den(
                    W_.a_,
                    W_.mom_,
                    &tsquare * Atom::var(W_.mass_) + Atom::var(GS.m_uv).pow(2),
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

    /// Add the numerator of the reduced subgraph, (without given), to the integrand.
    /// Then, 4d -> d-dim on minkowski indices
    #[debug_instrument(
        current = %current.log_display(),
        given = %given.log_display(),
        integrand = %integrand.log_display(),
    )]
    pub(crate) fn start<S: super::ForestNodeLike>(
        &self,
        ctx: &UVCtx<'_>,
        current: &S,
        given: &S,
        integrand: &Atom,
    ) -> Result<Atom> {
        let reduced = current.reduced_subgraph(given);
        let graph = ctx.graph;

        let mut t_arg = ctx
            .graph
            .numerator(&reduced, given.subgraph())
            .to_d_dim(GS.dim)
            .get_single_atom()
            .unwrap();

        t_arg /= graph.denominator(&reduced, |_| 1);

        t_arg = t_arg
            .replace(GS.dim)
            .max_level(0)
            .with(Atom::var(GS.dim_epsilon) * (-2) + 4);

        debug_tags!(#uv, #integrated, #algebra, #start; log.integrand = integrand, reduced = %reduced.string_label(), "Start");

        Ok((t_arg * integrand).simplify_metrics())
    }

    #[debug_instrument(
        current = %current.log_display(),
        given = %given.log_display(),
        integrand = %integrand.log_display(),
    )]
    pub(crate) fn series_and_truncate<S: super::ForestNodeLike>(
        &self,
        ctx: &UVCtx<'_>,
        current: &S,
        given: &S,
        integrand: &Atom,
    ) -> Result<Atom> {
        let graph = ctx.graph;

        let n_loops = graph.n_loops(current.subgraph()) - graph.n_loops(given.subgraph());
        let series = integrand
            .series(GS.dim_epsilon, Atom::Zero, n_loops as i64 + 1)
            .unwrap();
        let series_atom = series.to_atom();

        debug_tags!(#uv, #integrated, #inspect, #series;
            log.series = series_atom,
            "dim epsilon Series "
        );

        let mut pole_stripped = Atom::Zero;

        for (power, p) in series.terms() {
            if power < 0 {
                pole_stripped += p * Atom::var(GS.dim_epsilon).pow(power);
            }
        }

        Ok(pole_stripped)
    }

    #[debug_instrument(
        current = %current.log_display(),
        given = %given.log_display(),
    )]
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

        let rescaled = graph.uv_rescaled(&reduced, n_loops, current.lmb(), integrand);
        debug_tags!(#uv,#integrated,#rescaled;log.res = rescaled, n_loops=%n_loops,"Rescaled expanded");

        let series = rescaled
            .series(GS.rescale, Atom::Zero, 0)
            .unwrap()
            .to_atom();
        debug_tags!(#uv,#integrated, #series;log.res = series, "Series expanded");

        let evalutated = series.replace(GS.rescale).with(Atom::num(1));
        debug_tags!(#uv,#integrated,#series;log.res = evalutated, "Evaluated at t = 1");

        let collected = evalutated
            .simplify_metrics()
            .collect_rep((Bispinor {}).into())
            .collect_gamma_chains();
        debug_tags!(#uv,#integrated,#collect;log.expr = collected, "After gamma chain collection");

        let schoonschip = collected
            .schoonschip_with_settings(&SchoonschipSettings {
                simplify_chain_like_functions: true,
                schoonschip_rank1_tensors: true,
                ..Default::default()
            })
            .normalize_chains();
        debug_tags!(#uv, #integrated, #profile, #trace, #start, #collect;
            log.expr = schoonschip,
            "After gamma schoonschip"
        );
        let collected = schoonschip
            .collect_chains_and_traces()
            .simplify_metrics()
            .collect_gamma_chains()
            .collect_color()
            .collect_factors();
        debug_tags!(#uv, #integrated, #profile, #trace, #start, #collect;
            log.expr = collected,
            "After gamma collection"
        );

        let simplified = collected.simplify_gamma();
        debug_tags!(#uv, #integrated, #vakint, #profile, #trace, #start, #gamma;
            log.expr = simplified,
            "After gamma simplification"
        );
        let schoonschipped = simplified.schoonschip_net::<Aind>();
        debug_tags!(#uv, #integrated, #vakint, #profile, #trace,#schoonschip, #start;
            log.expr = schoonschipped,
            "After Schoonschip net"
        );
        let dotted = schoonschipped.to_dots().normalize_dots();
        debug_tags!(#uv, #integrated, #vakint, #profile, #trace, #dots;
            log.expr = dotted,
            "After dots"
        );

        Ok(dotted)
    }

    #[debug_instrument(
        current = %current.log_display(),
        given = %given.log_display(),
        reduced,
    )]
    pub(crate) fn integrate<S: super::ForestNodeLike>(
        &self,
        ctx: &UVCtx<'_>,
        current: &S,
        given: &S,
        integrand: &Atom,
    ) -> Result<Atom> {
        let graph = ctx.graph;
        let reduced = current.reduced_subgraph(given);
        let settings = ctx.settings;
        let reduced_label = reduced.string_label();
        tracing::Span::current().record("reduced", reduced_label.as_str());
        debug_tags!(#uv, #integrated, #vakint, #trace, #input;
            log.integrand = integrand,
            "Integrating and truncating"
        );
        // Match the vacuum result to the unintegrated CFF loop measure through
        // Vakint's configured per-loop normalization hook. Forest-subtraction
        // signs are folded into the expression at the integrated CT composition
        // sites.
        let mut integrand_vakint = to_vakint_integrand(
            integrand,
            graph,
            &reduced,
            given.subgraph(),
            &settings.vakint,
            true,
        );

        for (term_index, t) in integrand_vakint.0.iter().enumerate() {
            debug_tags!(#uv,#integrated,#vakint,#trace,#to_vakint;
                term_index = %term_index,
                log.integral = t.integral,
                log.numerator = t.numerator,
                "Vakint term as input"
            );
        }
        debug_tags!(#uv,#integrated,#vakint;settings = ?&self.vakint_settings,"Vakint args");

        // let mut res = vakint
        //     .0
        //     .evaluate(&vakint.1, integrand_vakint.as_view())
        //     .unwrap();

        integrand_vakint.canonicalize(self.vakint_settings, &self.vakint.topologies, false)?;
        for (term_index, t) in integrand_vakint.0.iter().enumerate() {
            debug_tags!(#uv,#integrated,#vakint,#trace,#canonicalize;
                term_index = %term_index,
                log.integral = t.integral,
                log.numerator = t.numerator,
                "Vakint term after canonicalization"
            );
        }
        integrand_vakint.tensor_reduce(self.vakint, self.vakint_settings)?;
        for (term_index, t) in integrand_vakint.0.iter().enumerate() {
            debug_tags!(#uv,#integrated,#vakint,#trace,#tensor_reduce;
                term_index = %term_index,
                log.integral = t.integral,
                log.numerator = t.numerator,
                "Vakint term after tensor reduction"
            );
        }
        integrand_vakint.evaluate_integral(self.vakint, self.vakint_settings)?;
        for (term_index, t) in integrand_vakint.0.iter().enumerate() {
            debug_tags!(#uv,#integrated,#vakint,#trace,#evaluate;
                term_index = %term_index,
                log.integral = t.integral,
                log.numerator = t.numerator,
                "Vakint term after evaluation"
            );
        }

        let mut res: Atom = integrand_vakint.into();

        debug_tags!(#uv,#integrated,#vakint,#trace,#raw;
            log.res = res,
            "Raw post vakint "
        );

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
            .when(W_.j_.filter(|r| r.is_integer()))
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
            .when(W_.j_.filter(|r| r.is_integer()))
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
            .when(W_.x_.filter(|r| r.is_integer()))
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
            .when(W_.y_.filter(|r| r.is_integer()))
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

        debug_tags!(#uv, #integrated, #vakint, #inspect, #trace, #replace;
            log.res = res,
            "Replaced post vakint "
        );

        // This strips as many dummies as possible after undoing chains and traces,
        // so that terms can merge later on.
        let bispinor_rep = Bispinor {}.into();
        let after_chainify = res.chainify(bispinor_rep);
        debug_tags!(#uv, #integrated, #vakint, #profile, #trace, #chainify;
            log.expr = after_chainify,
            "Integrated UV chain cleanup after chainify"
        );

        let after_collect_chains = after_chainify.collect_chains(bispinor_rep);
        debug_tags!(#uv, #integrated, #vakint, #profile, #trace, #collect;
            log.expr = after_collect_chains,
            "Integrated UV chain cleanup after collect_chains"
        );

        res = after_collect_chains.undo_single_length();
        debug_tags!(#uv, #integrated, #vakint, #profile, #trace, #undo_single_length;
            log.expr = res,
            "Integrated UV chain cleanup after undo_single_length"
        );

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

impl ApproximationKernel<UVCtx<'_>> for Integrated<'_> {
    #[debug_instrument(
        current = %current.log_display(),
        given = %given.log_display(),
        integrand = %integrand.log_display(),
    )]
    fn kernel<S: ForestNodeLike>(
        &self,
        ctx: &UVCtx<'_>,
        current: &S,
        given: &S,
        integrand: &Atom,
    ) -> Result<Atom> {
        match current.renormalization_scheme() {
            ApproximationType::MUV => {
                let pole_part = if given.subgraph().is_empty() {
                    integrand.clone()
                } else {
                    // In ordinary CT generation nested integrated insertions enter as
                    // finite - full, i.e. minus their pole part. In pole-part mode the
                    // terminal forest sum performs the projection, so keep the pole part
                    // itself here.
                    let pole_part = self.series_and_truncate(ctx, current, given, integrand)?;
                    if ctx.settings.pole_part {
                        pole_part
                    } else {
                        -pole_part
                    }
                };
                let with_added_expr = self.start(ctx, current, given, &pole_part)?;
                let top = self.t(ctx, current, given, &with_added_expr)?;
                let result = self.integrate(ctx, current, given, &top)?;

                debug_tags!(#uv, #integrated, #vakint, #profile, #trace, #result;
                    log.result = result,
                    "Integrated UV after integrate_and_truncate"
                );
                Ok(result)
            }
            ApproximationType::IR => Err(eyre!("Not yet implemented IR")),
            ApproximationType::VaccuumLimit => Err(eyre!("Not yet implemented VaccuumLimit")),
            ApproximationType::OS => Err(eyre!("Not yet implemented OS")),
            ApproximationType::Unsubtracted => {
                panic!("should have been kept out of the wood");
            }
        }
    }
}

#[debug_instrument]
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
    let reduced_label = reduced.string_label();
    let dependent_subgraph_label = dependent_subgraph.string_label();
    let mut integrand_vakint = integrand
        .undo_schoonschip::<Aind>()
        .undo_chain::<Aind>()
        .undo_trace::<Aind>();
    debug_tags!(#uv, #integrated, #vakint, #trace;
        stage = "to_vakint_integrand_after_undo_shorthands",
        reduced = %reduced_label,
        dependent_subgraph = %dependent_subgraph_label,
        substitute_masses_to_m_uv = substitute_masses_to_m_uv,
        log.integrand = integrand_vakint,
        "Vakint trace after undo shorthands"
    );
    //Atom::Zero

    // strip the momentum wrapper from the denominator
    integrand_vakint = integrand_vakint
        // .replace(function!(
        //     GS.den,
        //     W_.prop_,
        //     function!(GS.emr_mom, W_.prop_, W_.mom_),
        //     W_.x__
        // ))
        // .with(function!(GS.den, W_.prop_, W_.mom_, W_.x__))
        .expand();
    debug_tags!(#uv, #integrated, #vakint, #trace;
        stage = "to_vakint_integrand_after_den_strip_expand",
        reduced = %reduced_label,
        dependent_subgraph = %dependent_subgraph_label,
        log.integrand = integrand_vakint,
        "Vakint trace after denominator strip and expand"
    );

    // Nested counterterms can expose a boundary metric next to the propagator
    // metric of the reduced graph. Contract those metric-only structures before
    // the expression is split into Vakint terms.
    integrand_vakint = integrand_vakint.simplify_metrics();
    debug_tags!(#uv, #integrated, #vakint, #trace;
        stage = "to_vakint_integrand_after_simplify_metrics",
        reduced = %reduced_label,
        dependent_subgraph = %dependent_subgraph_label,
        log.integrand = integrand_vakint,
        "Vakint trace after metric simplification"
    );

    let mut propagator_id = 1;

    let vk_prop = vakint::symbols::S.prop;
    let vk_edge = vakint_symbol!("edge");
    let vk_topo = vakint_symbol!("topo");

    // let contracted_nodes: BTreeSet<NodeIndex> = graph
    //     .iter_nodes_of(dependent_subgraph)
    //     .map(|(nid, c, v)| nid)
    //     .collect();

    debug_tags!(#uv, #integrated, #vakint, #graph, #dump;
        reduced = %graph.dot(reduced),
        "Den to prop for"
    );
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
    debug_tags!(#uv, #integrated, #vakint, #trace;
        stage = "to_vakint_integrand_after_den_to_prop",
        reduced = %reduced_label,
        dependent_subgraph = %dependent_subgraph_label,
        log.integrand = integrand_vakint,
        "Vakint trace after denominator-to-propagator conversion"
    );

    let mut first: Option<NodeIndex> = None;

    debug_tags!(#uv, #integrated, #vakint, #graph, #dump;
        reduced = %graph.dot(dependent_subgraph),
        "Shrinking subgraph for vakint"
    );
    // shrink vertices of the subgraph
    for (id, _crown, _data) in graph.iter_nodes_of(dependent_subgraph) {
        debug_tags!(#uv, #integrated, #vakint, #graph, #inspect;
            id = %id,
            "Shrinking Node"
        );

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
    debug_tags!(#uv, #integrated, #vakint, #trace;
        stage = "to_vakint_integrand_after_shrink_subgraph",
        reduced = %reduced_label,
        dependent_subgraph = %dependent_subgraph_label,
        log.integrand = integrand_vakint,
        "Vakint trace after shrinking subgraph"
    );

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
    debug_tags!(#uv, #integrated, #vakint, #trace;
        stage = "to_vakint_integrand_after_flip_fuse",
        reduced = %reduced_label,
        dependent_subgraph = %dependent_subgraph_label,
        log.integrand = integrand_vakint,
        "Vakint trace after edge flip and fuse"
    );

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

    let vakint_input_atom = integrand_vakint
        .replace(function!(vk_prop, W_.x__))
        .with(function!(vk_topo, function!(vk_prop, W_.x__)))
        .replace(function!(vk_topo, W_.x_) * function!(vk_topo, W_.y_))
        .repeat()
        .with(function!(vk_topo, W_.x_ * W_.y_));
    debug_tags!(#uv, #integrated, #vakint, #trace;
        stage = "to_vakint_integrand_before_split_terms",
        reduced = %reduced_label,
        dependent_subgraph = %dependent_subgraph_label,
        log.integrand = vakint_input_atom,
        "Vakint trace before split terms"
    );

    let mut a = VakintExpression::try_from(vakint_input_atom).unwrap();

    for (term_index, t) in a.0.iter_mut().enumerate() {
        debug_tags!(#uv, #integrated, #vakint, #inspect, #trace;
            stage = "to_vakint_integrand_term_initial",
            term_index = %term_index,
            reduced = %reduced_label,
            dependent_subgraph = %dependent_subgraph_label,
            log.integral = t.integral,
            log.numerator = t.numerator,
            "Starting integral"
        );

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
        let uncontracted_propagator_count =
            graph.iter_edges().filter(|(p, _, _)| p.is_paired()).count();
        let uncontracted_propagator_power_sum = graph
            .iter_edges()
            .filter_map(|(p, _, e)| p.is_paired().then_some(e.data.power))
            .sum::<i32>();

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
        debug_tags!(#uv, #integrated, #vakint, #graph, #dump;
            log.graph = %graph.base_dot(),
            "Graph"
        );

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
        let contracted_propagator_count =
            graph.iter_edges().filter(|(p, _, _)| p.is_paired()).count();
        let contracted_propagator_power_sum = graph
            .iter_edges()
            .filter_map(|(p, _, e)| p.is_paired().then_some(e.data.power))
            .sum::<i32>();
        debug_tags!(#uv, #integrated, #vakint, #trace;
            stage = "to_vakint_integrand_term_after_graph_rebuild",
            term_index = %term_index,
            reduced = %reduced_label,
            dependent_subgraph = %dependent_subgraph_label,
            nloops = nloops,
            uncontracted_propagator_count = uncontracted_propagator_count,
            uncontracted_propagator_power_sum = uncontracted_propagator_power_sum,
            contracted_propagator_count = contracted_propagator_count,
            contracted_propagator_power_sum = contracted_propagator_power_sum,
            log.integral = t.integral,
            log.numerator = t.numerator,
            "Vakint trace"
        );

        let lmb = graph.lmb();
        let mom_pat = function!(GS.emr_mom, W_.a_).to_pattern();
        for (p, e, ed) in graph.iter_edges() {
            if p.is_paired() {
                // println!("{e}");
                let loop_expr = lmb.loop_atom::<Atom>(e, GS.loop_mom, &[], false);

                ed.data
                    .mom
                    .pattern_match(&mom_pat, None, None)
                    .for_each(|m| {
                        let var = mom_pat.replace_wildcards(&m).unwrap();
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
            .allow_new_wildcards_on_rhs(true),
            Replacement::new(
                function!(GS.loop_mom, W_.i_).to_pattern(),
                function!(GS.loop_mom, W_.i_, W_.a___),
            )
            .allow_new_wildcards_on_rhs(true),
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
        debug_tags!(#uv, #integrated, #vakint, #trace;
            stage = "to_vakint_integrand_term_after_momentum_solve",
            term_index = %term_index,
            reduced = %reduced_label,
            dependent_subgraph = %dependent_subgraph_label,
            log.integral = t.integral,
            log.numerator = t.numerator,
            "Vakint trace"
        );

        // debug!(
        //     "Graph from vakint expression:\n{}\n{}",
        //     graph.dot_lmb_of(&graph.full_filter(), &graph.lmb()),
        //     term
        // );

        let additional_normalization = parse!(&settings.additional_normalization);
        t.numerator *= additional_normalization.clone().pow(nloops);
        debug_tags!(#uv, #integrated, #vakint, #trace;
            stage = "to_vakint_integrand_term_after_loop_normalization",
            term_index = %term_index,
            reduced = %reduced_label,
            dependent_subgraph = %dependent_subgraph_label,
            nloops = nloops,
            log.additional_normalization = additional_normalization,
            log.numerator = t.numerator,
            "Vakint trace after loop normalization"
        );

        // Vakint needs explicit tensor indices; only translate metric shorthands
        // to dot notation here, without reintroducing Schoonschip rank-1 factors.
        t.numerator = t.numerator.metric_shorthand_to_dot();
        debug_tags!(#uv, #integrated, #vakint, #trace;
            stage = "to_vakint_integrand_term_after_metric_shorthand_to_dot",
            term_index = %term_index,
            reduced = %reduced_label,
            dependent_subgraph = %dependent_subgraph_label,
            log.integral = t.integral,
            log.numerator = t.numerator,
            "Vakint trace"
        );
        t.integral = t
            .integral
            .replace(function!(GS.loop_mom, W_.x___))
            .with(function!(vakint::symbols::S.k, W_.x___))
            .replace(function!(GS.emr_mom, W_.x___))
            .with(function!(vakint::symbols::S.p, W_.x___));
        t.numerator = t
            .numerator
            .replace(function!(GS.loop_mom, W_.x___))
            .with(function!(vakint::symbols::S.k, W_.x___))
            .replace(function!(GS.emr_mom, W_.x___))
            .with(function!(vakint::symbols::S.p, W_.x___))
            .replace(function!(
                SPENSO_TAG.dot,
                function!(W_.a_, W_.a___, Minkowski {}.new_rep(GS.dim).to_symbolic([])),
                function!(W_.b_, W_.b___, Minkowski {}.new_rep(GS.dim).to_symbolic([]))
            ))
            .with(vakint::symbols::S.dot(function!(W_.a_, W_.a___), function!(W_.b_, W_.b___)))
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
        debug_tags!(#uv, #integrated, #vakint, #trace;
            stage = "to_vakint_integrand_term_after_vakint_symbols",
            term_index = %term_index,
            reduced = %reduced_label,
            dependent_subgraph = %dependent_subgraph_label,
            log.integral = t.integral,
            log.numerator = t.numerator,
            "Vakint trace"
        );
    }

    a
}

#[cfg(test)]
mod tests {
    use crate::initialisation::test_initialise;

    use super::*;

    // #[test]
    // fn integrated_triangle_norm_is_euclidean() {
    //     test_initialise().unwrap();

    //     let edge = EdgeIndex(7);
    //     let euclidean_norm = integrated_triangle_spatial_norm_sq(edge);
    //     let minkowski_norm = Minkowski {}
    //         .new_rep(4)
    //         .inner_product(GS.emr_vec(edge), GS.emr_vec(edge));

    //     assert_eq!(
    //         euclidean_norm,
    //         GS.emr_mom(edge, GS.cind(1)).pow(2)
    //             + GS.emr_mom(edge, GS.cind(2)).pow(2)
    //             + GS.emr_mom(edge, GS.cind(3)).pow(2)
    //     );
    //     assert_ne!(euclidean_norm, minkowski_norm);
    // }

    #[test]
    fn vakint_dot_conversion_keeps_loop_momentum_tagged_until_to_dots() {
        test_initialise().unwrap();

        let mink = Minkowski {}.new_rep(GS.dim).to_symbolic([]);
        let numerator = function!(
            ETS.metric,
            function!(GS.emr_mom, 0, mink.clone()),
            function!(GS.loop_mom, 1, mink.clone())
        );

        let converted = numerator
            .simplify_metrics()
            .to_dots()
            .replace(function!(GS.loop_mom, W_.x___))
            .with(function!(vakint::symbols::S.k, W_.x___))
            .replace(function!(GS.emr_mom, W_.x___))
            .with(function!(vakint::symbols::S.p, W_.x___))
            .replace(function!(
                SPENSO_TAG.dot,
                function!(W_.a_, W_.a___, Minkowski {}.new_rep(GS.dim).to_symbolic([])),
                function!(W_.b_, W_.b___, Minkowski {}.new_rep(GS.dim).to_symbolic([]))
            ))
            .with(vakint::symbols::S.dot(function!(W_.a_, W_.a___), function!(W_.b_, W_.b___)));

        assert_eq!(
            converted,
            vakint::symbols::S.dot(
                function!(vakint::symbols::S.p, 0),
                function!(vakint::symbols::S.k, 1)
            )
        );
    }
}
