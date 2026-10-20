use color_eyre::Result;
use eyre::{WrapErr, eyre};
use gammaloop_tracing_filter::debug_instrument;
use idenso::{
    color::ColorSimplifier,
    dirac::{GammaSimplifier, GammaSimplifySettings},
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
use symbolica::{
    atom::{Atom, AtomCore},
    domains::atom::AtomField,
    function,
    id::Replacement,
    parse, parse_lit,
    poly::series::Series,
    solve::SolveCoverage,
};
use symbolica_utils::ReplaceBuilderExt;
use vakint::{Vakint, VakintExpression, vakint_symbol};

use crate::{
    debug_tags,
    graph::LMBext,
    numerator::aind::Aind,
    utils::{GS, W_},
    uv::{
        ApproximationType, UltravioletGraph,
        approx::{ForestNodeLike, Rooted, UVCtx, local_4d::Local4dCts},
        marker::{UvMarker, UvOperation},
        settings::VakintSettings,
        uv_graph::UVE,
    },
};

/// Laurent projections of an integrated counterterm.
///
/// Connected values store their canonical expansion. Factorized values encode
/// componentwise pole and signed nonnegative-power products under the same accessors.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct IntegratedCts {
    expansion: Series<AtomField>,
    scale_power: i64,
}

impl IntegratedCts {
    pub(crate) fn factorized_product<'a>(
        factors: impl IntoIterator<Item = &'a Self>,
        depth: usize,
    ) -> Result<Self> {
        let mut factors = factors.into_iter();
        let first = factors
            .next()
            .ok_or_else(|| eyre!("a factorized integrated counterterm cannot be empty"))?;
        let mut pole = truncate(&first.expansion, false);
        let mut finite_counterterm = -truncate(&first.expansion, true);
        let mut scale_power = first.scale_power;

        for factor in factors {
            pole *= truncate(&factor.expansion, false);
            finite_counterterm *= -truncate(&factor.expansion, true);
            scale_power += factor.scale_power;
        }

        Ok(Self {
            // The public projections add the finite-counterterm sign, so store
            // its negative while keeping the pole product unchanged.
            expansion: series(&(pole - finite_counterterm), depth)?,
            scale_power,
        })
    }

    fn projected_atom(&self, finite: bool) -> Atom {
        truncate(&self.expansion, finite)
            * Atom::var(GS.integrated_loop_scale).pow(self.scale_power)
    }

    pub(crate) fn pole_atom(&self) -> Atom {
        self.projected_atom(false)
    }

    pub(crate) fn finite_counterterm_atom(&self) -> Atom {
        -self.projected_atom(true)
    }

    pub(crate) fn physical_pole_atom(&self) -> Atom {
        self.pole_atom()
            .replace(GS.integrated_loop_scale)
            .with(Atom::one())
    }

    pub(crate) fn physical_finite_counterterm_atom(&self) -> Atom {
        self.finite_counterterm_atom()
            .replace(GS.integrated_loop_scale)
            .with(Atom::one())
    }
}

fn series(expr: &Atom, depth: usize) -> Result<Series<AtomField>> {
    Ok(expr.series(GS.dim_epsilon, 0, depth)?)
}

fn truncate(series: &Series<AtomField>, finite: bool) -> Atom {
    let mut truncated = Atom::Zero;

    for (power, p) in series.terms() {
        if (power >= 0) == finite {
            truncated += p * Atom::var(GS.dim_epsilon).pow(power);
        }
    }

    truncated
}

impl Rooted for IntegratedCts {
    fn root() -> Self {
        Self {
            expansion: series(&Atom::Zero, 1).expect("zero has a Laurent expansion"),
            scale_power: 0,
        }
    }
}

fn simplify(integrand: &Atom) -> Result<Atom> {
    let collected = integrand
        .collect_rep(Minkowski {}.into())
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

    let simplified = collected
        .simplify_gamma_with(GammaSimplifySettings::canonical())
        .collect_rep(Minkowski {}.into())
        .expand_num();
    debug_tags!(#uv, #integrated, #vakint, #profile, #trace, #start, #gamma;
        log.expr = simplified,
        "After gamma simplification"
    );
    let schoonschipped = simplified.schoonschip_net::<Aind>()?;
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

pub(crate) struct Integrated<'a> {
    pub vakint: &'a Vakint,
    pub vakint_settings: &'a vakint::VakintSettings,
}

impl Integrated<'_> {
    pub(crate) fn new<'a>(
        vakint: &'a Vakint,
        vakint_settings: &'a vakint::VakintSettings,
    ) -> Integrated<'a> {
        Integrated {
            vakint,
            vakint_settings,
        }
    }

    pub(crate) fn run<S: super::ForestNodeLike, M: super::ForestNodeLike>(
        &self,
        integrand: &Local4dCts,
        ctx: &UVCtx<'_>,
        current: &S,
        given: &S,
        marker_current: &M,
        marker_given: &M,
    ) -> Result<IntegratedCts> {
        let graph = ctx.graph;

        let n_loops = graph.n_loops(current.subgraph());

        let scheme = current.renormalization_scheme();
        match scheme {
            ApproximationType::MUV | ApproximationType::PolePart => {
                let integrand = integrand
                    .atom()
                    .replace(GS.integrated_loop_scale)
                    .with(Atom::one());
                let simplified = simplify(&integrand)?;
                let marker = UvMarker::new(ctx.settings);
                let integrated = marker.apply(
                    UvOperation::Integrate,
                    marker_current.subgraph(),
                    marker_given.subgraph(),
                    &self.integrate(&simplified, ctx, current, given)?,
                );
                let expansion_depth =
                    usize::try_from(self.vakint_settings.number_of_terms_in_epsilon_expansion)
                        .wrap_err("Vakint epsilon expansion depth must be nonnegative")?;
                let expanded = series(&integrated, expansion_depth.max(n_loops + 1))?.map_coeff(
                    |coefficient| {
                        marker.apply(
                            UvOperation::Series,
                            marker_current.subgraph(),
                            marker_given.subgraph(),
                            coefficient,
                        )
                    },
                );
                let expansion = expanded.map_coeff(|coefficient| {
                    marker.apply(
                        UvOperation::Truncate,
                        marker_current.subgraph(),
                        marker_given.subgraph(),
                        coefficient,
                    )
                });

                // Retain the consumed loop measures for subsequent UV rescalings. Keep this
                // marker independent of mUV so enclosing limits still rescale the vacuum mass.
                Ok(IntegratedCts {
                    expansion,
                    scale_power: 4 * n_loops as i64,
                })
            }
            ApproximationType::IR => Err(eyre!("Not yet implemented IR")),
            ApproximationType::VaccuumLimit => Err(eyre!("Not yet implemented VaccuumLimit")),
            ApproximationType::OS => Err(eyre!("Not yet implemented OS")),
            ApproximationType::Unsubtracted => {
                panic!("should have been kept out of the wood");
            }
        }
    }

    #[debug_instrument(
        current = %current.log_display(),
        given = %given.log_display(),
        reduced,
    )]
    fn integrate<S: ForestNodeLike>(
        &self,
        integrand: &Atom,
        ctx: &UVCtx<'_>,
        current: &S,
        given: &S,
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
            current.subgraph(),
            given.subgraph(),
            &settings.vakint,
            true,
        )?;

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
        for (term_index, t) in integrand_vakint.0.iter_mut().enumerate() {
            // Projector slots and retained source slots share one native
            // Lorentz representation. Restore that metadata before the
            // backend's Laurent expansion, keeping scalar epsilon factors.
            t.numerator = Self::restore_lorentz_dimensions(&t.numerator);
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

        res = Self::restore_numerator(res, current.topo_order());

        res = Self::dimensionally_regularized(&res.simplify_metrics().metric_shorthand_to_dot());

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

    fn dimensionally_regularized(numerator: &Atom) -> Atom {
        let minkowski = Minkowski {}.to_symbolic([GS.dim]).get_symbol().unwrap();
        let dimension = Atom::var(GS.dim);
        let regulated = Atom::num(4) - Atom::num(2) * Atom::var(GS.dim_epsilon);
        numerator.replace_map(|view, _, out| {
            if view.get_symbol() == Some(minkowski) {
                // A retained open Lorentz basis keeps its representation. All
                // scalar occurrences, including inside function coefficients,
                // must instead participate in the Laurent expansion.
                **out = view.to_owned();
            } else if view == dimension.as_view() {
                **out = regulated.clone();
            }
        })
    }

    fn restore_lorentz_dimensions(numerator: &Atom) -> Atom {
        let mink = Minkowski {}.new_rep(GS.dim);
        let minkowski = mink.to_symbolic([]).get_symbol().unwrap();
        let regulated_dimension = Atom::num(4) - Atom::num(2) * Atom::var(GS.dim_epsilon);
        // The projector works in 4-2*epsilon while opaque source slots retain
        // d. These are the same Lorentz representation; unify only metadata,
        // leaving every scalar evanescent coefficient for the later expansion.
        numerator
            .replace(function!(minkowski, &regulated_dimension, W_.x___))
            .with(function!(minkowski, GS.dim, W_.x___))
    }

    fn restore_numerator(mut res: Atom, topology_order: usize) -> Atom {
        let mink = Minkowski {}.new_rep(GS.dim);
        res = Self::restore_lorentz_dimensions(&res);
        res = res
            .replace(parse_lit!(vakint::cl2))
            .with(parse_lit!(cl2))
            .replace(parse_lit!(vakint::sqrt3))
            .with(parse_lit!(sqrt(3)));

        for (source, target) in [
            (vakint::symbols::S.k, GS.loop_mom),
            (vakint::symbols::S.p, GS.emr_mom),
        ] {
            res = res
                .replace(function!(source, W_.i_, W_.j_))
                .when(W_.j_.filter(|index| index.is_integer().is_true()))
                .with(function!(
                    source,
                    W_.i_,
                    mink.to_symbolic([function!(GS.uvaind, topology_order, W_.j_)])
                ))
                .replace(function!(
                    source,
                    W_.i_,
                    vakint::symbols::S.dot_dummy_ind(W_.j_)
                ))
                .when(W_.j_.filter(|index| index.is_integer().is_true()))
                .with(function!(
                    source,
                    W_.i_,
                    mink.to_symbolic([function!(GS.uvaind, topology_order, W_.j_)])
                ))
                // Native projection retains a source index's identity while
                // omitting its representation. Only numeric backend dummies
                // belong to the topology-local chart above.
                .replace(function!(source, W_.i_, W_.j_))
                .when(W_.j_.filter(|index| Aind::try_from(index).is_ok()))
                .with(function!(
                    source,
                    W_.i_,
                    mink.to_symbolic([Atom::var(W_.j_)])
                ))
                .replace(function!(source, W_.x___))
                .with(function!(target, W_.x___));
        }

        let vk_metric = vakint_symbol!("g");
        res = res
            .replace(function!(vk_metric, W_.x_, W_.y_) * function!(GS.emr_mom, W_.x___, W_.x_))
            .with(function!(GS.emr_mom, W_.x___, W_.y_));
        res = res.replace_map(|view, _, output| {
            let symbolica::atom::AtomView::Fun(metric) = view else {
                return;
            };
            if metric.get_symbol() != vk_metric || metric.get_nargs() != 2 {
                return;
            }
            let slots = [metric.get(0), metric.get(1)].map(|slot| {
                let raw_index = if let symbolica::atom::AtomView::Fun(index) = slot
                    && index.get_symbol() == vakint::symbols::S.dot_dummy_ind
                    && index.get_nargs() == 1
                {
                    index.get(0)
                } else {
                    slot
                };
                if raw_index.is_integer().is_true() {
                    mink.to_symbolic([function!(GS.uvaind, topology_order, raw_index)])
                } else if Aind::try_from(raw_index).is_ok() {
                    mink.to_symbolic([raw_index.to_owned()])
                } else {
                    slot.to_owned()
                }
            });
            **output = function!(ETS.metric, &slots[0], &slots[1]);
        });

        res.replace(vakint::symbols::S.cmplx_i).with(Atom::i())
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
) -> Result<VakintExpression> {
    let reduced_label = reduced.string_label();
    let dependent_subgraph_label = dependent_subgraph.string_label();
    let mut integrand_vakint = integrand
        .undo_schoonschip::<Aind>()?
        .undo_chain::<Aind>()?
        .undo_trace::<Aind>()?;
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
                        GS.m_uv_vacuum
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

    let mut a = VakintExpression::try_from(vakint_input_atom)
        .wrap_err("could not split integrand into Vakint terms")?;

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
        let mut momentum_variables = vec![];

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
                &e.data
                    .mass
                    .pow(2)
                    .replace(GS.m_uv_expansion)
                    .with(GS.m_uv_vacuum),
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
                        if !momentum_variables.iter().any(|existing| existing == &var) {
                            momentum_variables.push(var);
                        }
                    });

                // println!("{loop_expr}");

                // println!("{external_expr}");
                let is_zero = &ed.data.mom - loop_expr;
                // println!("Momentum check for edge {}: {}", e, is_zero);
                system.push(is_zero);
            }
        }

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
        let momentum_solution =
            VakintMomentumSolution::solve(&system, &momentum_variables, &add_additional_args)
                .wrap_err("could not solve momentum system for Vakint integrand")?;
        t.numerator = momentum_solution.rewrite_numerator(&t.numerator);
        t.integral = momentum_solution.rewrite_integral(&t.integral, &add_additional_args);
        momentum_solution.ensure_free_variables_eliminated(
            term_index,
            "integral",
            &t.integral,
            &add_additional_args,
        )?;
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

    Ok(a)
}

struct VakintMomentumSolution {
    replacements: Vec<Replacement>,
    free_variables: Vec<Atom>,
}

impl VakintMomentumSolution {
    fn solve(
        system: &[Atom],
        variables: &[Atom],
        add_additional_args: &[Replacement],
    ) -> Result<Self> {
        if variables.is_empty() {
            return Ok(Self {
                replacements: vec![],
                free_variables: vec![],
            });
        }
        if system.is_empty() {
            return Ok(Self {
                replacements: vec![],
                free_variables: variables.to_vec(),
            });
        }

        // Momentum routing is linear, but a subgraph can leave coordinates
        // free. Preserve those explicitly without discarding native guards.
        Atom::system_to_matrix::<u8, _, _>(system, variables)?;
        let solutions = Atom::solve(system).wrt_with_exponent::<u8, _>(variables)?;
        if solutions.coverage() != SolveCoverage::Complete
            || !solutions.coverage_guard().is_empty()
            || solutions.len() != 1
            || solutions[0].is_conditional()
        {
            return Err(eyre!(
                "momentum routing requires one complete unconditional solution family: {solutions:?}"
            ));
        }
        let branch = &solutions[0];
        let solution = branch
            .coordinate_order()
            .map(|variable| {
                branch
                    .get(variable)
                    .cloned()
                    .unwrap_or_else(|| variable.to_atom())
            })
            .collect::<Vec<_>>();
        Ok(Self::from_solution(
            &solution,
            variables,
            add_additional_args,
        ))
    }

    fn from_solution(
        solution: &[Atom],
        variables: &[Atom],
        add_additional_args: &[Replacement],
    ) -> Self {
        debug_assert_eq!(solution.len(), variables.len());
        let mut replacements = vec![];
        let mut free_variables = vec![];

        for index in (0..variables.len()).rev() {
            let replacement = &solution[index];
            let variable = &variables[index];
            let replacement = replacement.replace_multiple(&replacements);
            if replacement == variable {
                free_variables.push(variable.clone());
            } else {
                replacements.push(Replacement::new(
                    variable.replace_multiple(add_additional_args).to_pattern(),
                    replacement
                        .replace_multiple(add_additional_args)
                        .to_pattern(),
                ));
            }
        }

        Self {
            replacements,
            free_variables,
        }
    }

    fn ensure_free_variables_eliminated(
        &self,
        term_index: usize,
        expression_kind: &str,
        expression: &Atom,
        add_additional_args: &[Replacement],
    ) -> Result<()> {
        self.ensure_variables_eliminated(
            term_index,
            expression_kind,
            expression,
            self.free_variables.iter(),
            add_additional_args,
        )
    }

    fn ensure_variables_eliminated<'a>(
        &self,
        term_index: usize,
        expression_kind: &str,
        expression: &Atom,
        variables: impl IntoIterator<Item = &'a Atom>,
        add_additional_args: &[Replacement],
    ) -> Result<()> {
        for variable in variables {
            let bare_pattern = variable.to_pattern();
            let indexed_pattern = variable.replace_multiple(add_additional_args).to_pattern();
            if expression
                .pattern_match(&bare_pattern, None, None)
                .next()
                .is_some()
                || expression
                    .pattern_match(&indexed_pattern, None, None)
                    .next()
                    .is_some()
            {
                return Err(eyre!(
                    "Underdetermined Vakint momentum solve left free variable {} in {} of term {}",
                    variable,
                    expression_kind,
                    term_index
                ));
            }
        }
        Ok(())
    }

    fn rewrite_numerator(&self, numerator: &Atom) -> Atom {
        numerator
            .replace_multiple(&self.replacements)
            .normalize_dots()
    }

    fn rewrite_integral(&self, integral: &Atom, add_additional_args: &[Replacement]) -> Atom {
        let free_variable_zero_replacements =
            self.free_variable_zero_replacements(add_additional_args);
        Self::normalize_topology_momenta(
            integral
                .replace_multiple(&self.replacements)
                .replace_multiple(&free_variable_zero_replacements),
        )
    }

    fn free_variable_zero_replacements(
        &self,
        add_additional_args: &[Replacement],
    ) -> Vec<Replacement> {
        let mut replacements = Vec::with_capacity(2 * self.free_variables.len());
        for variable in &self.free_variables {
            replacements.push(Replacement::new(variable.to_pattern(), Atom::Zero));
            let indexed_variable = variable.replace_multiple(add_additional_args);
            if indexed_variable != *variable {
                replacements.push(Replacement::new(indexed_variable.to_pattern(), Atom::Zero));
            }
        }
        replacements
    }

    fn normalize_topology_momenta(expression: Atom) -> Atom {
        if !expression
            .replace(function!(vakint::symbols::S.topo, W_.x_))
            .matches()
        {
            return expression;
        }

        expression
            .replace(function!(
                vakint::symbols::S.prop,
                W_.edgeid_,
                W_.x_,
                W_.mom_,
                W_.mass_,
                W_.prop_
            ))
            .with_map(|matches| {
                function!(
                    vakint::symbols::S.prop,
                    matches.get(W_.edgeid_).unwrap().to_atom(),
                    matches.get(W_.x_).unwrap().to_atom(),
                    matches.get(W_.mom_).unwrap().to_atom().expand(),
                    matches.get(W_.mass_).unwrap().to_atom(),
                    matches.get(W_.prop_).unwrap().to_atom()
                )
            })
    }
}

#[cfg(test)]
mod tests {
    use crate::initialisation::test_initialise;
    use spenso::structure::slot::{DummyAind, ParseableAind};

    use super::*;

    #[test]
    fn projector_lorentz_metadata_is_constant_during_laurent_expansion() {
        test_initialise().unwrap();
        let epsilon = Atom::var(GS.dim_epsilon);
        let mink = Minkowski {}.new_rep(GS.dim);
        let spin = Bispinor {}.new_rep(4);
        let a = spin.to_symbolic([Aind::Hedge(3, 1).to_atom()]);
        let b = spin.to_symbolic([Aind::Hedge(4, 1).to_atom()]);
        let mu = mink.to_symbolic([Aind::Hedge(12, 1).to_atom()]);
        let dummy = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let regulated_mu = mu
            .replace(GS.dim)
            .with(Atom::num(4) - Atom::num(2) * &epsilon);
        let gamma = function!(idenso::dirac::AGS.gamma, &a, &b, &dummy);
        let scalar = (Atom::var(GS.dim) - 4) * (&epsilon + 1) / &epsilon;
        let source = &scalar * &gamma * function!(vakint::symbols::S.g, &dummy, regulated_mu);
        assert!(
            series(&source, 1)
                .unwrap()
                .to_atom()
                .contains_symbol(symbolica::atom::Symbol::DERIVATIVE),
            "a regulated tensor slot must be restored before scalar Laurent algebra"
        );
        let normalized = Integrated::restore_lorentz_dimensions(&source);
        let expected = &scalar * gamma * function!(vakint::symbols::S.g, &dummy, &mu);
        assert_eq!(normalized, expected);
        assert_eq!(
            Integrated::restore_lorentz_dimensions(&normalized),
            normalized,
            "metadata restoration is idempotent"
        );
        let expanded = series(&normalized, 1).unwrap().to_atom();
        assert!(!expanded.contains_symbol(symbolica::atom::Symbol::DERIVATIVE));
        let restored = Integrated::restore_numerator(expanded, 0).simplify_metrics();
        let regularized = Integrated::dimensionally_regularized(&restored);
        let expected = Atom::num(-2) * function!(idenso::dirac::AGS.gamma, &a, &b, &mu);
        assert!(
            (regularized.clone() - (&epsilon + 1) * &expected)
                .expand()
                .is_zero()
        );
        let finite = truncate(&series(&regularized, 0).unwrap(), true);
        assert!((finite - expected).expand().is_zero());
        for head in [vakint::symbols::S.k, vakint::symbols::S.p] {
            assert_eq!(
                Integrated::restore_lorentz_dimensions(&function!(head, 1, &mu)),
                function!(head, 1, &mu),
                "dimension normalization must retain the Vakint momentum ABI"
            );
        }
    }

    #[test]
    fn scalar_dimensions_are_substituted_inside_functions_but_not_lorentz_slots() {
        test_initialise().unwrap();
        let slot = Minkowski {}
            .new_rep(GS.dim)
            .to_symbolic([Aind::new_dummy().to_atom()]);
        let logarithm = symbolica::symbol!("log");
        let tensor = symbolica::symbol!("dimension_test_tensor");
        let input = function!(logarithm, GS.dim) * function!(tensor, &slot);
        let expected = function!(
            logarithm,
            Atom::num(4) - Atom::num(2) * Atom::var(GS.dim_epsilon)
        ) * function!(tensor, &slot);
        assert_eq!(Integrated::dimensionally_regularized(&input), expected);
        let dimension = Atom::num(4) - Atom::num(2) * Atom::var(GS.dim_epsilon);
        let projected_slot = slot.replace(GS.dim).with(dimension);
        let evanescent = Atom::var(GS.dim) - Atom::num(4);
        assert_eq!(
            Integrated::restore_numerator(&evanescent * function!(tensor, projected_slot), 0),
            evanescent * function!(tensor, &slot),
            "restoration unifies Lorentz metadata without changing scalar d dependence"
        );
        let other_slot = Minkowski {}
            .new_rep(6)
            .to_symbolic([Aind::new_dummy().to_atom()]);
        let other = function!(tensor, &other_slot);
        assert_eq!(Integrated::restore_numerator(other.clone(), 0), other);
        assert_eq!(Integrated::dimensionally_regularized(&other), other);
    }

    #[test]
    fn restored_projector_metadata_closes_mixed_native_gamma_sum() {
        test_initialise().unwrap();
        let mink = Minkowski {}.new_rep(GS.dim);
        let spin = Bispinor {}.new_rep(4);
        let a = spin.to_symbolic([Aind::Hedge(3, 1).to_atom()]);
        let b = spin.to_symbolic([Aind::Hedge(4, 1).to_atom()]);
        let mu = mink.to_symbolic([Aind::Hedge(12, 1).to_atom()]);
        let dummy = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let regulated_mu = mu
            .replace(GS.dim)
            .with(Atom::num(4) - Atom::num(2) * Atom::var(GS.dim_epsilon));
        let gamma = function!(idenso::dirac::AGS.gamma, &a, &b, &mu);
        let epsilon = Atom::var(GS.dim_epsilon);
        let scalar = (Atom::var(GS.dim) - 4) / &epsilon;
        let source = &scalar
            * function!(idenso::dirac::AGS.gamma, &a, &b, &dummy)
            * function!(vakint::symbols::S.g, &dummy, regulated_mu)
            + &gamma;
        let restored = Integrated::restore_numerator(source, 0).simplify_metrics();
        assert!(
            (restored.clone() - (&scalar + 1) * &gamma)
                .expand()
                .is_zero()
        );
        let regularized = Integrated::dimensionally_regularized(&restored);
        assert!((regularized + &gamma).expand().is_zero());
        assert_eq!(
            series(
                &(Integrated::dimensionally_regularized(&scalar) * &epsilon),
                1
            )
            .unwrap()
            .to_atom()
            .expand(),
            -2 * epsilon,
            "scalar evanescent terms must survive the metadata restoration"
        );
    }

    #[test]
    fn restored_vakint_source_indices_keep_identity_and_representation() {
        test_initialise().unwrap();

        let mink = Minkowski {}.new_rep(GS.dim);
        let other_dimension = Minkowski {}.new_rep(6);
        for (source, target) in [
            (vakint::symbols::S.k, GS.loop_mom),
            (vakint::symbols::S.p, GS.emr_mom),
        ] {
            for index in [
                Aind::Hedge(5, 0),
                Aind::Vertex(2, 1),
                Aind::UVTerm(3, 4),
                Aind::Dummy(7),
                Aind::Symbol(symbolica::symbol!("restored_vakint_source_label")),
            ] {
                let index = Atom::from(index);
                assert_eq!(
                    Integrated::restore_numerator(function!(source, 2, &index), 19),
                    function!(target, 2, mink.to_symbolic([index]))
                );
            }

            let index = Atom::from(Aind::Hedge(5, 0));
            for slot in [
                mink.to_symbolic([index.clone()]),
                other_dimension.to_symbolic([index]),
            ] {
                assert_eq!(
                    Integrated::restore_numerator(function!(source, 2, &slot), 19),
                    function!(target, 2, slot)
                );
            }
            for dummy in [Atom::num(5), vakint::symbols::S.dot_dummy_ind(5)] {
                assert_eq!(
                    Integrated::restore_numerator(function!(source, 2, dummy), 19),
                    function!(target, 2, mink.to_symbolic([function!(GS.uvaind, 19, 5)]))
                );
            }
        }

        let source = Atom::from(Aind::Hedge(5, 0));
        let represented = other_dimension.to_symbolic([Atom::from(Aind::Hedge(7, 0))]);
        assert_eq!(
            Integrated::restore_numerator(
                function!(vakint::symbols::S.g, &source, &represented),
                19
            ),
            function!(ETS.metric, mink.to_symbolic([source]), represented)
        );
        assert_eq!(
            Integrated::restore_numerator(
                function!(vakint::symbols::S.g, 5, vakint::symbols::S.dot_dummy_ind(7)),
                19
            ),
            function!(
                ETS.metric,
                mink.to_symbolic([function!(GS.uvaind, 19, 5)]),
                mink.to_symbolic([function!(GS.uvaind, 19, 7)])
            )
        );
    }

    #[test]
    fn restored_vakint_source_indices_close_mixed_tensor_and_dot_sums() {
        test_initialise().unwrap();

        let mink = Minkowski {}.new_rep(GS.dim);
        let h5 = Atom::from(Aind::Hedge(5, 0));
        let h7 = Atom::from(Aind::Hedge(7, 0));
        let q1 = function!(GS.emr_mom, 1, mink.to_symbolic([]));
        let q4 = function!(GS.emr_mom, 4, mink.to_symbolic([]));
        let source = function!(vakint::symbols::S.p, 1, &h5)
            * function!(vakint::symbols::S.p, 1, &h7)
            * function!(GS.emr_mom, 1, mink.to_symbolic([h5]))
            * function!(GS.emr_mom, 4, mink.to_symbolic([h7]))
            + function!(SPENSO_TAG.dot, &q1, &q4);

        // A nested counterterm has restored source labels beside still-native
        // propagator slots. Both summands must have the same closed interface
        // before the next subgraph is materialized as a tensor network.
        let restored = Integrated::restore_numerator(source, 19)
            .undo_schoonschip::<Aind>()
            .unwrap()
            .to_dots();
        let expected =
            (function!(SPENSO_TAG.dot, &q1, &q1) + 1) * function!(SPENSO_TAG.dot, &q1, &q4);
        assert!((restored - expected).expand().is_zero());
    }

    #[test]
    fn integrated_counterterm_projects_one_laurent_expansion() {
        test_initialise().unwrap();

        let epsilon = Atom::var(GS.dim_epsilon);
        let expansion = series(
            &(Atom::num(2) * epsilon.pow(-2)
                + Atom::num(3) * epsilon.pow(-1)
                + Atom::num(5)
                + Atom::num(7) * &epsilon),
            2,
        )
        .unwrap();
        let integrated = IntegratedCts {
            expansion,
            scale_power: 4,
        };
        let scale = Atom::var(GS.integrated_loop_scale).pow(4);

        assert_eq!(
            integrated.pole_atom().expand(),
            ((Atom::num(2) * epsilon.pow(-2) + Atom::num(3) * epsilon.pow(-1)) * &scale).expand()
        );
        assert_eq!(
            integrated.finite_counterterm_atom().expand(),
            (-(Atom::num(5) + Atom::num(7) * epsilon) * scale).expand()
        );
        assert_eq!(
            integrated.physical_finite_counterterm_atom(),
            -(Atom::num(5) + Atom::num(7) * Atom::var(GS.dim_epsilon))
        );
    }

    #[test]
    fn factorized_product_projects_each_component() {
        test_initialise().unwrap();

        let epsilon = Atom::var(GS.dim_epsilon);
        let factors = [2, 3, 5].map(|finite| IntegratedCts {
            expansion: series(&(epsilon.pow(-1) + Atom::num(finite)), 1).unwrap(),
            scale_power: 0,
        });
        let product = IntegratedCts::factorized_product(&factors[..2], 1).unwrap();

        assert_eq!(product.physical_pole_atom(), epsilon.pow(-2));
        assert_eq!(product.physical_finite_counterterm_atom(), Atom::num(6));

        let product = IntegratedCts::factorized_product(&factors, 1).unwrap();
        assert_eq!(product.physical_pole_atom(), epsilon.pow(-3));
        assert_eq!(product.physical_finite_counterterm_atom(), Atom::num(-30));
    }

    // #[test]
    // fn integrated_triangle_norm_is_euclidean() {
    //     test_initialise().unwrap();

    //     let edge = EdgeIndex(7);
    //     let euclidean_norm = GS.emr_mom(edge, GS.cind(1)).pow(2)
    //         + GS.emr_mom(edge, GS.cind(2)).pow(2)
    //         + GS.emr_mom(edge, GS.cind(3)).pow(2);
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

    #[test]
    fn underdetermined_vakint_momentum_solve_tracks_free_variables() {
        test_initialise().unwrap();

        let q0 = function!(GS.emr_mom, 0);
        let q1 = function!(GS.emr_mom, 1);
        let q2 = function!(GS.emr_mom, 2);
        let k0 = function!(GS.loop_mom, 0);
        let k1 = function!(GS.loop_mom, 1);
        let edge_3 = -&q0 - &q1 - &q2;
        let system = vec![&q0 - &k0, &edge_3 - &k1];
        let variables = vec![q0.clone(), q1.clone(), q2.clone()];

        let solution = VakintMomentumSolution::solve(&system, &variables, &[]).unwrap();

        assert_eq!(solution.free_variables, vec![q2]);
        assert!(
            solution
                .ensure_free_variables_eliminated(0, "numerator", &solution.free_variables[0], &[])
                .is_err()
        );
    }

    #[test]
    fn underdetermined_vakint_momentum_solve_projects_topology() {
        test_initialise().unwrap();

        let q1 = function!(GS.emr_mom, 1);
        let q2 = function!(GS.emr_mom, 2);
        let k0 = function!(GS.loop_mom, 0);
        let system = vec![-&q1 - &q2 - &k0];
        let variables = vec![q1.clone(), q2.clone()];
        let topology = function!(
            vakint::symbols::S.topo,
            function!(
                vakint::symbols::S.prop,
                1,
                function!(vakint::symbols::S.edge, 0, 0),
                -&q1 - &q2,
                GS.m_uv_vacuum,
                1
            )
        );

        let solution = VakintMomentumSolution::solve(&system, &variables, &[]).unwrap();

        solution
            .ensure_free_variables_eliminated(0, "numerator", &Atom::from(1), &[])
            .unwrap();
        assert_eq!(
            solution.rewrite_integral(&topology, &[]),
            function!(
                vakint::symbols::S.topo,
                function!(
                    vakint::symbols::S.prop,
                    1,
                    function!(vakint::symbols::S.edge, 0, 0),
                    k0,
                    GS.m_uv_vacuum,
                    1
                )
            )
        );
    }

    #[test]
    fn empty_vakint_momentum_solve_treats_variables_as_free() {
        test_initialise().unwrap();

        let q0 = function!(GS.emr_mom, 0);
        let no_variables = VakintMomentumSolution::solve(&[], &[], &[]).unwrap();
        let free_variable =
            VakintMomentumSolution::solve(&[], std::slice::from_ref(&q0), &[]).unwrap();

        assert!(no_variables.free_variables.is_empty());
        assert_eq!(free_variable.free_variables, vec![q0]);
    }

    #[test]
    fn vakint_momentum_solve_rejects_conditional_and_nonlinear_routings() {
        test_initialise().unwrap();

        let q = function!(GS.emr_mom, 0);
        let k = function!(GS.loop_mom, 0);
        let parameter = parse!("routing_parameter");
        // q/parameter = 0 retains an input-domain guard even when its solution
        // simplifies to q = 0; neither that guard nor a pivot guard may vanish.
        for equation in [&parameter * &q - &k, &q / &parameter, q.pow(2) - &k] {
            assert!(
                VakintMomentumSolution::solve(&[equation], std::slice::from_ref(&q), &[]).is_err()
            );
        }
    }
}
