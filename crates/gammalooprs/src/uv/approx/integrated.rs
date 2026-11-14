use std::{borrow::Cow, collections::HashMap, time::Instant};

use color_eyre::Result;
use eyre::{WrapErr, eyre};
use gammaloop_tracing_filter::debug_instrument;
use idenso::{
    CookMode, CookSettings, IndexTooling,
    dirac::GammaSimplifySettings,
    shorthands::UndoShorthands,
    tensor::{AlgebraSettings, ContractSettings, SymbolicNetExt, SymbolicNetParse, SymbolicTensor},
};

use linnet::half_edge::{
    HedgeGraph, NodeIndex,
    builder::HedgeGraphBuilder,
    involution::HedgePair,
    subgraph::{ModifySubSet, SuBitGraph, SubGraphLike, SubSetLike},
};
use spenso::{
    network::{
        library::symbolic::ETS,
        parsing::{AtomStructureExt, ParseSettings, SchoonschipExpansionMode, ShorthandParsing},
        tags::SPENSO_TAG,
    },
    shadowing::symbolica_utils::LogPrint,
    structure::representation::{BaseRepName, Minkowski, RepName},
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder},
    domains::atom::AtomField,
    function,
    id::Replacement,
    parse, parse_lit,
    poly::series::Series,
    solve::SolveCoverage,
};
use symbolica_utils::ReplaceBuilderExt;
use vakint::{Vakint, VakintExpression, stage_timing::StageSpan, vakint_symbol};

use crate::{
    debug_tags,
    graph::LMBext,
    numerator::aind::Aind,
    utils::{GS, W_},
    uv::{
        ApproximationType, UltravioletGraph,
        approx::{ForestNodeLike, Rooted, UVCtx, local_4d::Local4dCts, streamed},
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
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let simplified = SymbolicTensor::infer_indexed(cooking.try_cook(integrand.as_view())?)?
        .simplify_algebra(&AlgebraSettings {
            gamma: Some(GammaSimplifySettings::canonical()),
            epsilon: true,
            ..Default::default()
        })?
        .to_dots()?;
    let dotted = cooking
        .uncook(simplified.expression().as_view())
        .expand_num();
    debug_tags!(#uv, #integrated, #vakint, #profile, #trace, #gamma;
        log.expr = dotted,
        "Shared tensor and Dirac simplification of integrated UV numerator"
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
        let _stage = StageSpan::new("gamma.integrated.run");
        let graph = ctx.graph;

        let n_loops = graph.n_loops(current.subgraph());

        let scheme = current.renormalization_scheme();
        match scheme {
            ApproximationType::MUV | ApproximationType::PolePart => {
                let simplify_stage = StageSpan::new("gamma.integrated.pre_simplify");
                let integrand = integrand
                    .atom()
                    .replace(GS.integrated_loop_scale)
                    .with(Atom::one());
                let simplified = simplify(&integrand)?;
                drop(simplify_stage);
                let marker = UvMarker::new(ctx.settings);
                let integrated = marker.apply(
                    UvOperation::Integrate,
                    marker_current.subgraph(),
                    marker_given.subgraph(),
                    &self.integrate(&simplified, ctx, current, given)?,
                );
                let expansion_stage = StageSpan::new("gamma.integrated.final_epsilon_series");
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
                drop(expansion_stage);

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
        let _stage = StageSpan::new("gamma.integrated.integrate");
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

        let topology_stage = StageSpan::new("gamma.integrated.topology_canonicalization");
        integrand_vakint.canonicalize(self.vakint_settings, &self.vakint.topologies, false)?;
        drop(topology_stage);
        for (term_index, t) in integrand_vakint.0.iter().enumerate() {
            debug_tags!(#uv,#integrated,#vakint,#trace,#canonicalize;
                term_index = %term_index,
                log.integral = t.integral,
                log.numerator = t.numerator,
                "Vakint term after canonicalization"
            );
        }
        // GammaLoop's source dimension denotes the regulated Lorentz space.
        // Declare that identification before tensor reduction: opaque source
        // tensors keep their input slots, and RustRed correctly refuses to
        // infer this relation during scalar admission. Scalar d dependence
        // must also be explicit before master-order planning: e.g. 1/(d-4)
        // promotes the finite master coefficient into the requested pole.
        let dimension_stage = StageSpan::new("gamma.integrated.dimensional_preparation");
        let minkowski = Minkowski::selfless_symbol();
        let dimension = Atom::num(4) - Atom::num(2) * Atom::var(GS.dim_epsilon);
        for term in &mut integrand_vakint.0 {
            let numerator = Self::dimensionally_regularized(&term.numerator);
            term.numerator = numerator
                .replace(function!(minkowski, GS.dim, W_.x___))
                .with(function!(minkowski, &dimension, W_.x___));
        }
        drop(dimension_stage);
        let tensor_stage = StageSpan::new("gamma.integrated.tensor_reduce");
        integrand_vakint.tensor_reduce(self.vakint, self.vakint_settings)?;
        drop(tensor_stage);
        for (term_index, t) in integrand_vakint.0.iter().enumerate() {
            debug_tags!(#uv,#integrated,#vakint,#trace,#tensor_reduce;
                term_index = %term_index,
                log.integral = t.integral,
                log.numerator = t.numerator,
                "Vakint term after tensor reduction"
            );
        }
        let evaluation_stage = StageSpan::new("gamma.integrated.evaluate_integral");
        integrand_vakint.evaluate_integral(self.vakint, self.vakint_settings)?;
        drop(evaluation_stage);
        for (term_index, t) in integrand_vakint.0.iter().enumerate() {
            debug_tags!(#uv,#integrated,#vakint,#trace,#evaluate;
                term_index = %term_index,
                log.integral = t.integral,
                log.numerator = t.numerator,
                "Vakint term after evaluation"
            );
        }

        let assembly_stage = StageSpan::new("gamma.integrated.assemble_result");
        let mut res: Atom = integrand_vakint.into();
        drop(assembly_stage);

        debug_tags!(#uv,#integrated,#vakint,#trace,#raw;
            log.res = res,
            "Raw post vakint "
        );

        let restore_stage = StageSpan::new("gamma.integrated.restore_numerator");
        res = Self::restore_numerator(res, current.topo_order());
        drop(restore_stage);

        let cleanup_stage = StageSpan::new("gamma.integrated.final_gamma_metric_algebra");
        res = Self::reduce_integrated_numerator(&res)?;
        drop(cleanup_stage);
        debug_tags!(#uv, #integrated, #vakint, #profile, #trace, #collect;
            log.expr = res,
            "Integrated UV chain cleanup through native tensor contraction"
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

    fn reduce_integrated_numerator(numerator: &Atom) -> Result<Atom> {
        // Vacuum tensor reduction can create contracted Clifford words even
        // when the incoming numerator was already simplified. Reduce those
        // words before taking Laurent coefficients: gamma_mu gamma^mu = d
        // contributes evanescent terms to poles and finite counterterms.
        // Open Lorentz-slot dimensions remain metadata; scalar dimensions
        // produced by the Clifford reduction must become 4 - 2 epsilon.
        let reduced = simplify(numerator)?;
        Self::warn_repeated_lorentz_indices(&reduced, "R", None);
        Ok(Self::dimensionally_regularized(&reduced))
    }

    pub(crate) fn warn_repeated_lorentz_indices(
        numerator: &Atom,
        operation: &'static str,
        reduction_status: Option<idenso::tensor::ReductionStatus>,
    ) {
        if numerator.has_repeated_lorentz_indices() {
            // One event per completed operation, not one per index or summand.
            tracing::warn!(
                uv = true,
                integrated = true,
                operation,
                reduction_status = ?reduction_status,
                stage = "unresolved_lorentz_contraction_before_poles",
                display.preview = %numerator.log_print(Some(160)).chars().take(160).collect::<String>(),
                file.expr = %numerator.to_plain_string(),
                "Repeated abstract Lorentz index remains after tensor simplification; unresolved contractions can generate D = 4 - 2 epsilon factors and change Laurent poles. Review this result before using its counterterm."
            );
        }
    }

    pub(crate) fn dimensionally_regularized(numerator: &Atom) -> Atom {
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

    pub(crate) fn restore_lorentz_dimensions(numerator: &Atom) -> Atom {
        let mink = Minkowski {}.new_rep(GS.dim);
        let minkowski = mink.to_symbolic([]).get_symbol().unwrap();
        let regulated_dimension = Atom::num(4) - Atom::num(2) * Atom::var(GS.dim_epsilon);
        // The projector works in 4-2*epsilon while opaque source slots retain
        // d. These are the same Lorentz representation; unify only metadata,
        // leaving every scalar evanescent coefficient for the later expansion.
        // Backend d/epsilon substitutions can return an unexpanded but equal
        // dimension such as 4-(4-(4-2*epsilon)).
        numerator.replace_map(|view, _, out| {
            let AtomView::Fun(representation) = view else {
                return;
            };
            if representation.get_symbol() == minkowski
                && representation.get_nargs() > 0
                && (representation.get(0).to_owned() - &regulated_dimension)
                    .expand()
                    .is_zero()
            {
                **out = FunctionBuilder::new(minkowski)
                    .add_arg(GS.dim)
                    .add_args(representation.iter().skip(1))
                    .finish();
            }
        })
    }

    /// Expose tensor-valued vectors inside chains before covariant reduction.
    ///
    /// Scalar products remain compact, including inverse external invariants.
    /// Chain endpoints survive and newly opened contractions receive fresh
    /// dummy indices. Materialize before tensor-interface validation because
    /// an identically zero antisymmetric chain may have apparent open slots.
    ///
    /// Materialization is additive over summands and multiplicative over
    /// factors that share no index, so a large expression is materialized in
    /// pieces (see [`streamed::map_multilinear`]) instead of as one tensor
    /// network, which holds about fifty times the bytes of its input. Pieces
    /// that are pure scalars without shorthands are kept as they are, which
    /// is what the parser does with them. If a piece introduces a dummy index
    /// (numbered over the whole expression), the whole expression is
    /// materialized at once.
    pub(crate) fn materialize_tensor(expression: &Atom) -> Result<Atom> {
        Self::materialize_tensor_view(expression.as_view())
    }

    pub(crate) fn materialize_tensor_view(expression: AtomView<'_>) -> Result<Atom> {
        let limit = streamed::piece_limit();
        let materialized = if streamed::nothing_to_materialize(expression) {
            expression.to_owned()
        } else if expression.get_byte_size() >= limit
            && let Some(materialized) = streamed::map_multilinear(
                expression,
                limit,
                &streamed::nothing_to_materialize,
                &Self::materialize_whole,
            )?
        {
            materialized
        } else {
            return Self::materialize_whole(expression);
        };
        // The streamed result has no new dummy index, so the whole one is
        // repeatable and must be identical.
        if streamed::check_streamed(expression) {
            let whole = Self::materialize_whole(expression)?;
            streamed::report_check("materialize", expression, &materialized, &whole);
        }
        Ok(materialized)
    }

    /// One tensor network of the whole expression: the reference for
    /// [`Self::materialize_tensor`].
    fn materialize_whole(expression: AtomView<'_>) -> Result<Atom> {
        let explicit = expression.undo_trace::<Aind>()?.undo_chain::<Aind>()?;
        Ok(explicit
            .parse_to_symbolic_net::<Aind>(&ParseSettings {
                shorthand_parsing: ShorthandParsing::Expand {
                    schoonschip: SchoonschipExpansionMode {
                        inner_products: false,
                        expand_inside_chains: true,
                        expand_schoonship: true,
                    },
                    trace: true,
                    chain: true,
                },
                ..Default::default()
            })?
            .simple_execute::<()>()?)
    }

    pub(crate) fn restore_numerator(mut res: Atom, topology_order: usize) -> Atom {
        let mink = Minkowski {}.new_rep(GS.dim);
        res = Self::restore_lorentz_dimensions(&res);
        // Native tensor reduction also puts these backend dummies inside
        // represented gamma/metric slots. Restore their identity before
        // handling momentum and metric heads, so every occurrence shares the
        // same topology-local index without adding a second representation.
        res = res
            .replace(vakint::symbols::S.dot_dummy_ind(W_.j_))
            .when(W_.j_.filter(|index| index.is_integer().is_true()))
            .with(function!(GS.uvaind, topology_order, W_.j_));
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
                let raw_index = slot;
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

/// An integrand passed to [`to_vakint_integrand`] by reference, or by value
/// so that it can be dropped as soon as it has been materialized.
pub(crate) trait IntoIntegrand<'a> {
    fn into_integrand(self) -> Cow<'a, Atom>;
}

impl<'a> IntoIntegrand<'a> for &'a Atom {
    fn into_integrand(self) -> Cow<'a, Atom> {
        Cow::Borrowed(self)
    }
}

impl IntoIntegrand<'static> for Atom {
    fn into_integrand(self) -> Cow<'static, Atom> {
        Cow::Owned(self)
    }
}

/// The summand-local rewrites of [`to_vakint_integrand`] after expansion.
struct TermRewrites {
    propagators: HashMap<usize, (usize, Atom)>,
    contracted_nodes: HashMap<usize, usize>,
    substitute_masses_to_m_uv: bool,
    /// [`Self::after_contraction_whole`] of each distinct product of
    /// denominator factors met so far.
    denominator_rewrites: std::cell::RefCell<HashMap<Atom, Atom>>,
    /// [`Self::contract_metrics`] of each distinct product of indexed factors.
    contraction_images: std::cell::RefCell<HashMap<Atom, Atom>>,
    /// Census: summands and bytes entering the term rewrites.
    rewritten_summands: std::cell::Cell<(usize, usize)>,
}

impl TermRewrites {
    // Nested counterterms can expose a boundary metric next to the propagator
    // metric of the reduced graph. Contract those metric-only structures before
    // the expression is split into Vakint terms.
    fn contract_metrics(integrand: AtomView<'_>) -> Result<Atom> {
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let contracted = SymbolicTensor::infer_indexed(cooking.try_cook(integrand)?)?.contract(
            ContractSettings {
                rank_one: false,
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            },
        )?;
        Ok(cooking.uncook(contracted.expression().as_view()))
    }

    /// [`Self::contract_metrics`], skipped when the expression has neither a
    /// metric nor a repeated explicit index, so that nothing can contract.
    fn contract_metrics_if_any(&self, integrand: AtomView<'_>) -> Result<Atom> {
        if !integrand.contains_symbol(ETS.metric) && !integrand.has_repeated_explicit_indices() {
            return Ok(integrand.to_owned());
        }
        self.contract_metrics_by_indexed_factors(integrand)
    }

    /// [`Self::contract_metrics`] summand by summand on the indexed factors
    /// only, once per distinct product of them (see
    /// [`streamed::map_indexed_factors`]). The contraction of a five-loop
    /// integrand otherwise builds one tensor network of up to millions of
    /// summands that differ only in index-free factors.
    /// `GAMMALOOP_ALGEBRA_LEGACY=metric_contraction` restores the whole
    /// contraction, `GAMMALOOP_ALGEBRA_CHECK=metric_contraction` compares.
    fn contract_metrics_by_indexed_factors(&self, integrand: AtomView<'_>) -> Result<Atom> {
        if streamed::legacy("metric_contraction") {
            return Self::contract_metrics(integrand);
        }
        let mapped = streamed::map_indexed_factors(
            integrand,
            &mut self.contraction_images.borrow_mut(),
            &|piece| Self::contract_metrics(piece),
        )?;
        let contracted = match mapped {
            Some(contracted) => contracted,
            None => Self::contract_metrics(integrand)?,
        };
        if streamed::check("metric_contraction") {
            let whole = Self::contract_metrics(integrand)?;
            streamed::report_check("metric_contraction", integrand, &contracted, &whole);
        }
        Ok(contracted)
    }

    // Denominator identities do not overlap. Visit the expression once rather
    // than rebuilding its complete numerator for every graph edge. Normalize
    // powers after all leaves have been replaced, using the same rule as before.
    fn denominators_to_propagators(&self, integrand: Atom) -> Atom {
        if self.propagators.is_empty() {
            return integrand;
        }
        let vk_prop = vakint::symbols::S.prop;
        let integrand = integrand.replace_map(|part, _, out| {
            let AtomView::Fun(denominator) = part else {
                return;
            };
            if denominator.get_symbol() != GS.den || denominator.get_nargs() < 3 {
                return;
            }
            let Ok(edge) = usize::try_from(denominator.get(0)) else {
                return;
            };
            let Some((id, endpoints)) = self.propagators.get(&edge) else {
                return;
            };
            let mass = if self.substitute_masses_to_m_uv {
                Atom::var(GS.m_uv_vacuum)
            } else {
                denominator.get(2).to_owned()
            };
            **out = function!(vk_prop, *id, endpoints, denominator.get(1), mass, 1);
        });
        integrand
            .replace(function!(vk_prop, W_.x___, 1).pow(Atom::var(W_.e_)))
            .with(function!(vk_prop, W_.x___, -Atom::var(W_.e_)))
    }

    fn shrink_subgraph(&self, integrand: Atom) -> Atom {
        if self.contracted_nodes.is_empty() {
            return integrand;
        }
        let vk_prop = vakint::symbols::S.prop;
        let vk_edge = vakint_symbol!("edge");
        integrand.replace_map(|part, _, out| {
            let AtomView::Fun(propagator) = part else {
                return;
            };
            if propagator.get_symbol() != vk_prop || propagator.get_nargs() < 2 {
                return;
            }
            let AtomView::Fun(endpoints) = propagator.get(1) else {
                return;
            };
            if endpoints.get_symbol() != vk_edge || endpoints.get_nargs() != 2 {
                return;
            }
            let node = |original: AtomView<'_>| {
                usize::try_from(original)
                    .ok()
                    .and_then(|id| self.contracted_nodes.get(&id).copied())
                    .map_or_else(|| original.to_owned(), Atom::num)
            };
            let edge = function!(vk_edge, node(endpoints.get(0)), node(endpoints.get(1)));
            **out = FunctionBuilder::new(vk_prop)
                .add_arg(propagator.get(0))
                .add_arg(edge)
                .add_args(propagator.iter().skip(2))
                .finish();
        })
    }

    fn flip_and_fuse(integrand: Atom) -> Atom {
        let vk_prop = vakint::symbols::S.prop;
        let vk_edge = vakint_symbol!("edge");
        // flip edges to positive momentum
        // FIXME: how will this work for sums of momenta?
        let integrand = integrand
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
        integrand
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
            ))
    }

    fn wrap_topologies(integrand: &Atom) -> Atom {
        let vk_prop = vakint::symbols::S.prop;
        let vk_topo = vakint_symbol!("topo");
        integrand
            .replace(function!(vk_prop, W_.x__))
            .with(function!(vk_topo, function!(vk_prop, W_.x__)))
            .replace(function!(vk_topo, W_.x_) * function!(vk_topo, W_.y_))
            .repeat()
            .with(function!(vk_topo, W_.x_ * W_.y_))
    }

    /// Every rewrite after the metric contraction.
    ///
    /// The rewrites only touch the `den` factors of a summand (they become
    /// `prop`s, are shrunk, flipped, fused and wrapped into one `topo`); all
    /// other factors pass through. A five-loop integrand has a few hundred
    /// distinct denominator products but up to millions of summands, so each
    /// distinct product is rewritten once (by the whole map) and every summand
    /// is rebuilt as its other factors times that result. Pattern matching no
    /// longer runs over the numerators. `GAMMALOOP_ALGEBRA_LEGACY=term_rewrites`
    /// restores the whole map, `GAMMALOOP_ALGEBRA_CHECK=term_rewrites`
    /// compares both.
    fn after_contraction(&self, integrand: Atom) -> Atom {
        if streamed::legacy("term_rewrites") {
            return self.after_contraction_whole(integrand);
        }
        let rewritten = self.after_contraction_by_denominators(integrand.as_view());
        if streamed::check("term_rewrites") {
            let input = integrand.clone();
            let whole = self.after_contraction_whole(integrand);
            streamed::report_check("term_rewrites", input.as_view(), &rewritten, &whole);
        }
        rewritten
    }

    /// The reference for [`Self::after_contraction`]: every rewrite on the
    /// whole expression.
    fn after_contraction_whole(&self, integrand: Atom) -> Atom {
        let integrand = self.denominators_to_propagators(integrand);
        let integrand = self.shrink_subgraph(integrand);
        let integrand = Self::flip_and_fuse(integrand);
        Self::wrap_topologies(&integrand)
    }

    /// [`Self::after_contraction_whole`] applied summand by summand to the
    /// denominator factors only. A summand is a product `N * D` of non-`den`
    /// factors `N` and `den` factors (or their powers) `D`; every pattern of
    /// the rewrites matches inside `D` only (a `prop`, a product of `prop`s
    /// or `topo`s), so the image of the summand is `N * whole(D)`. Normalized
    /// sums and products are canonical, so rebuilding the sum of these images
    /// gives the atom of the whole map. Summands whose other factors mention
    /// a `den`, `prop`, `edge` or `topo` (inside a function or a sum), or
    /// that have no denominator, take the whole map.
    fn after_contraction_by_denominators<'a>(&self, integrand: AtomView<'a>) -> Atom {
        let vk_prop = vakint::symbols::S.prop;
        let vk_edge = vakint_symbol!("edge");
        let vk_topo = vakint_symbol!("topo");
        let is_denominator = |factor: AtomView<'_>| match factor {
            AtomView::Fun(function) => function.get_symbol() == GS.den,
            AtomView::Pow(power) => matches!(
                power.get_base_exp().0,
                AtomView::Fun(function) if function.get_symbol() == GS.den
            ),
            _ => false,
        };
        if vakint::census::enabled() {
            let (summands, bytes) = self.rewritten_summands.get();
            self.rewritten_summands.set((
                summands + vakint::census::summands(integrand),
                bytes + integrand.get_byte_size(),
            ));
        }
        let mut cache = self.denominator_rewrites.borrow_mut();
        let mut denominators = Vec::new();
        let mut others = Vec::new();
        let mut rewrite = |term: AtomView<'a>| -> Atom {
            denominators.clear();
            others.clear();
            match term {
                AtomView::Mul(product) => {
                    for factor in product.iter() {
                        if is_denominator(factor) {
                            denominators.push(factor);
                        } else {
                            others.push(factor);
                        }
                    }
                }
                factor if is_denominator(factor) => denominators.push(factor),
                factor => others.push(factor),
            }
            if denominators.is_empty()
                || others.iter().any(|factor| {
                    factor.contains_symbol(GS.den)
                        || factor.contains_symbol(vk_prop)
                        || factor.contains_symbol(vk_edge)
                        || factor.contains_symbol(vk_topo)
                })
            {
                return self.after_contraction_whole(term.to_owned());
            }
            let key = match denominators[..] {
                [single] => single.to_owned(),
                _ => Atom::mul_many(denominators.iter().copied()),
            };
            let image = cache
                .entry(key)
                .or_insert_with_key(|key| self.after_contraction_whole(key.clone()));
            if others.is_empty() {
                return image.clone();
            }
            Atom::mul_many(
                others
                    .iter()
                    .copied()
                    .chain(std::iter::once(image.as_view())),
            )
        };
        match integrand {
            AtomView::Add(sum) => {
                let terms = sum.iter().map(&mut rewrite).collect::<Vec<_>>();
                Atom::add_many(terms.iter())
            }
            term => rewrite(term),
        }
    }

    /// The whole map from the materialized integrand to the Vakint input, as
    /// one expression: the reference for [`Self::streamed`].
    fn whole(&self, materialized: AtomView<'_>) -> Result<Atom> {
        let expanded = materialized.expand();
        let contracted = Self::contract_metrics(expanded.as_view())?;
        drop(expanded);
        Ok(self.after_contraction_whole(contracted))
    }

    /// [`Self::whole`] on pieces of the expansion of `materialized`, in
    /// chunks of summands, added up once at the end. Expansion, contraction
    /// and the rewrites act on every summand alone and are linear in its
    /// coefficient, and normalized sums are canonical, so this is the whole
    /// map. Returns `None` if a contraction introduces a dummy index.
    fn streamed(&self, materialized: AtomView<'_>, limit: usize) -> Result<Option<Atom>> {
        let mut accumulated = streamed::SumAccumulator::default();
        let mut introduced = false;
        let outcome = streamed::expand_in_pieces(materialized, limit, &mut |piece| {
            streamed::for_each_summand_chunk(piece.as_view(), limit, &mut |chunk| {
                let contracted = self.contract_metrics_if_any(chunk)?;
                if streamed::introduces_dummies(chunk, contracted.as_view()) {
                    introduced = true;
                    return Err(eyre!("a metric contraction introduced a dummy index"));
                }
                accumulated.push(self.after_contraction(contracted));
                Ok(())
            })
        });
        match outcome {
            Ok(()) => Ok(Some(accumulated.finish())),
            Err(_) if introduced => Ok(None),
            Err(error) => Err(error),
        }
    }
}

#[debug_instrument]
pub(crate) fn to_vakint_integrand<
    'a,
    E: UVE,
    V,
    H,
    S: SubGraphLike + SubSetLike<Base = SuBitGraph>,
    SS: SubGraphLike,
>(
    integrand: impl IntoIntegrand<'a>,
    graph: &HedgeGraph<E, V, H>,
    reduced: &S,
    dependent_subgraph: &SS,
    settings: &VakintSettings,
    substitute_masses_to_m_uv: bool,
) -> Result<VakintExpression> {
    let _stage = StageSpan::new("gamma.adaptor.to_vakint");
    let started = Instant::now();
    let reduced_label = reduced.string_label();
    let dependent_subgraph_label = dependent_subgraph.string_label();
    let materialization_stage = StageSpan::new("gamma.adaptor.materialize_tensor");
    let integrand = integrand.into_integrand();
    let materialized = Integrated::materialize_tensor(&integrand)?;
    // An integrand passed by value is not needed any more.
    drop(integrand);
    let census_materialized_bytes = materialized.as_view().get_byte_size();
    let materialization_seconds = started.elapsed().as_secs_f64();
    drop(materialization_stage);
    debug_tags!(#uv, #integrated, #vakint, #trace;
        stage = "to_vakint_integrand_after_undo_shorthands",
        reduced = %reduced_label,
        dependent_subgraph = %dependent_subgraph_label,
        substitute_masses_to_m_uv = substitute_masses_to_m_uv,
        log.integrand = materialized,
        "Vakint trace after undo shorthands"
    );

    let vk_prop = vakint::symbols::S.prop;
    let vk_edge = vakint_symbol!("edge");
    debug_tags!(#uv, #integrated, #vakint, #graph, #dump;
        reduced = %graph.dot(reduced),
        "Den to prop for"
    );
    let mut propagators = HashMap::new();
    for (pair, index, _data) in graph.iter_edges_of(reduced) {
        if let HedgePair::Paired { source, sink } = pair {
            propagators.insert(
                usize::from(index),
                (
                    propagators.len() + 1,
                    function!(
                        vk_edge,
                        usize::from(graph.node_id(source)),
                        usize::from(graph.node_id(sink))
                    ),
                ),
            );
        }
    }
    debug_tags!(#uv, #integrated, #vakint, #graph, #dump;
        reduced = %graph.dot(dependent_subgraph),
        "Shrinking subgraph for vakint"
    );
    // Each disconnected counterterm has its own contracted vertex. Identifying
    // distinct children with each other would create extra quotient cycles
    // without corresponding independent momentum coordinates.
    let mut contracted_nodes = HashMap::new();
    for component in graph.connected_components(dependent_subgraph) {
        let mut first: Option<NodeIndex> = None;
        for (id, _crown, _data) in graph.iter_nodes_of(&component) {
            debug_tags!(#uv, #integrated, #vakint, #graph, #inspect;
                id = %id,
                "Shrinking Node"
            );

            if let Some(first) = first {
                contracted_nodes.insert(id.0, first.0);
            } else {
                first = Some(id);
            }
        }
    }
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let rewrites = TermRewrites {
        propagators,
        contracted_nodes,
        substitute_masses_to_m_uv,
        denominator_rewrites: Default::default(),
        contraction_images: Default::default(),
        rewritten_summands: Default::default(),
    };

    // A large integrand is expanded, contracted and rewritten in pieces, so
    // that neither its unmerged expansion nor a cooked copy, observation map
    // or contraction of the whole expansion is alive at once.
    let limit = streamed::piece_limit();
    let streamed_input = if materialized.as_view().get_byte_size() >= limit {
        let expansion_stage = StageSpan::new("gamma.adaptor.expand");
        let streamed_input = rewrites.streamed(materialized.as_view(), limit)?;
        drop(expansion_stage);
        if let Some(streamed_input) = &streamed_input
            && streamed::check_streamed(materialized.as_view())
        {
            let whole = rewrites.whole(materialized.as_view())?;
            streamed::report_check("to_vakint", materialized.as_view(), streamed_input, &whole);
        }
        streamed_input
    } else {
        None
    };
    let (vakint_input_atom, expansion_seconds, metric_seconds, topology_seconds) =
        if let Some(streamed_input) = streamed_input {
            drop(materialized);
            let seconds = started.elapsed().as_secs_f64() - materialization_seconds;
            (streamed_input, seconds, 0.0, 0.0)
        } else {
            let expansion_started = Instant::now();
            let expansion_stage = StageSpan::new("gamma.adaptor.expand");
            let integrand_vakint = materialized.expand();
            drop(materialized);
            let expansion_seconds = expansion_started.elapsed().as_secs_f64();
            drop(expansion_stage);
            debug_tags!(#uv, #integrated, #vakint, #trace;
                stage = "to_vakint_integrand_after_den_strip_expand",
                reduced = %reduced_label,
                dependent_subgraph = %dependent_subgraph_label,
                log.integrand = integrand_vakint,
                "Vakint trace after denominator strip and expand"
            );

            let metrics_started = Instant::now();
            let metric_stage = StageSpan::new("gamma.adaptor.metric_contraction");
            streamed::dump("contract_metrics", integrand_vakint.as_view());
            let contracted =
                rewrites.contract_metrics_by_indexed_factors(integrand_vakint.as_view())?;
            drop(integrand_vakint);
            let integrand_vakint = contracted;
            let metric_seconds = metrics_started.elapsed().as_secs_f64();
            drop(metric_stage);
            debug_tags!(#uv, #integrated, #vakint, #trace;
                stage = "to_vakint_integrand_after_simplify_metrics",
                reduced = %reduced_label,
                dependent_subgraph = %dependent_subgraph_label,
                log.integrand = integrand_vakint,
                "Vakint trace after metric simplification"
            );

            let topology_started = Instant::now();
            let topology_stage = StageSpan::new("gamma.adaptor.propagator_topology");
            let vakint_input_atom = if !streamed::legacy("term_rewrites") {
                rewrites.after_contraction(integrand_vakint)
            } else {
                let integrand_vakint = rewrites.denominators_to_propagators(integrand_vakint);
                debug_tags!(#uv, #integrated, #vakint, #trace;
                    stage = "to_vakint_integrand_after_den_to_prop",
                    reduced = %reduced_label,
                    dependent_subgraph = %dependent_subgraph_label,
                    log.integrand = integrand_vakint,
                    "Vakint trace after denominator-to-propagator conversion"
                );
                let integrand_vakint = rewrites.shrink_subgraph(integrand_vakint);
                debug_tags!(#uv, #integrated, #vakint, #trace;
                    stage = "to_vakint_integrand_after_shrink_subgraph",
                    reduced = %reduced_label,
                    dependent_subgraph = %dependent_subgraph_label,
                    log.integrand = integrand_vakint,
                    "Vakint trace after shrinking subgraph"
                );
                let integrand_vakint = TermRewrites::flip_and_fuse(integrand_vakint);
                debug_tags!(#uv, #integrated, #vakint, #trace;
                    stage = "to_vakint_integrand_after_flip_fuse",
                    reduced = %reduced_label,
                    dependent_subgraph = %dependent_subgraph_label,
                    log.integrand = integrand_vakint,
                    "Vakint trace after edge flip and fuse"
                );
                TermRewrites::wrap_topologies(&integrand_vakint)
            };
            let topology_seconds = topology_started.elapsed().as_secs_f64();
            drop(topology_stage);
            (
                vakint_input_atom,
                expansion_seconds,
                metric_seconds,
                topology_seconds,
            )
        };
    debug_tags!(#uv, #integrated, #vakint, #trace;
        stage = "to_vakint_integrand_before_split_terms",
        reduced = %reduced_label,
        dependent_subgraph = %dependent_subgraph_label,
        log.integrand = vakint_input_atom,
        "Vakint trace before split terms"
    );

    let census_input = vakint::census::enabled().then(|| {
        (
            vakint::census::summands(vakint_input_atom.as_view()),
            vakint_input_atom.as_view().get_byte_size(),
        )
    });
    let split_started = Instant::now();
    let split_stage = StageSpan::new("gamma.adaptor.split_terms");
    let mut a = VakintExpression::try_from(vakint_input_atom)
        .wrap_err("could not split integrand into Vakint terms")?;
    let split_seconds = split_started.elapsed().as_secs_f64();
    drop(split_stage);
    let additional_normalization = parse!(&settings.additional_normalization);
    let propagator_pattern = function!(
        vk_prop,
        W_.a_,
        function!(vakint::symbols::S.edge, W_.i_, W_.j_),
        W_.c_,
        W_.d_,
        W_.e_
    )
    .to_pattern();
    let momentum_pattern = function!(GS.emr_mom, W_.a_).to_pattern();
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
    let mut graph_seconds = 0.0;
    let mut notation_images = HashMap::new();
    let mut routing_seconds = 0.0;
    let mut notation_seconds = 0.0;

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

        let graph_started = Instant::now();
        let graph_stage = StageSpan::new("gamma.adaptor.rebuild_graph");
        let mut graph = HedgeGraphBuilder::new();
        //prop(<id>, edge(<node_source>,<node_sink>), <mom>, <mass>, <power>)
        let mut nodemap = HashMap::new();

        struct ContractibleEdge {
            mom: Atom,
            mass: Atom,
            power: i32,
        }

        for m in t.integral.pattern_match(&propagator_pattern, None, None) {
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
        graph_seconds += graph_started.elapsed().as_secs_f64();
        drop(graph_stage);
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

        let routing_started = Instant::now();
        let routing_stage = StageSpan::new("gamma.adaptor.momentum_routing");
        let lmb = graph.lmb();
        for (p, e, ed) in graph.iter_edges() {
            if p.is_paired() {
                // println!("{e}");
                let loop_expr = lmb.loop_atom::<Atom>(e, GS.loop_mom, &[], false);

                ed.data
                    .mom
                    .pattern_match(&momentum_pattern, None, None)
                    .for_each(|m| {
                        let var = momentum_pattern.replace_wildcards(&m).unwrap();
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

        let momentum_solution =
            VakintMomentumSolution::solve(&system, &momentum_variables, &add_additional_args)
                .wrap_err("could not solve momentum system for Vakint integrand")?;
        t.numerator = t
            .numerator
            .replace_multiple(&momentum_solution.replacements);
        t.integral = momentum_solution.rewrite_integral(&t.integral, &add_additional_args);
        momentum_solution.ensure_free_variables_eliminated(
            term_index,
            "integral",
            &t.integral,
            &add_additional_args,
        )?;
        routing_seconds += routing_started.elapsed().as_secs_f64();
        drop(routing_stage);
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

        let notation_started = Instant::now();
        let notation_stage = StageSpan::new("gamma.adaptor.tensor_notation");
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
        // Routing and loop normalization precede this one authoritative pass;
        // neither needs an intermediate cook/infer/to_dots/uncook cycle.
        streamed::dump("notation", t.numerator.as_view());
        let notation_whole = |numerator: AtomView<'_>| -> Result<Atom> {
            let notation =
                SymbolicTensor::infer_indexed(cooking.try_cook(numerator)?)?.to_dots()?;
            Ok(cooking.uncook(notation.expression().as_view()))
        };
        // Once per distinct product of indexed factors rather than one tensor
        // network of the whole numerator (`streamed::map_indexed_factors`).
        // `GAMMALOOP_ALGEBRA_LEGACY=dot_notation` restores the whole pass.
        let notated = if streamed::legacy("dot_notation") {
            notation_whole(t.numerator.as_view())?
        } else {
            let mapped = streamed::map_indexed_factors(
                t.numerator.as_view(),
                &mut notation_images,
                &notation_whole,
            )?;
            let notated = match mapped {
                Some(notated) => notated,
                None => notation_whole(t.numerator.as_view())?,
            };
            if streamed::check("dot_notation") {
                let whole = notation_whole(t.numerator.as_view())?;
                streamed::report_check("dot_notation", t.numerator.as_view(), &notated, &whole);
            }
            notated
        };
        t.numerator = notated;
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
        notation_seconds += notation_started.elapsed().as_secs_f64();
        drop(notation_stage);
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

    if let Some((output_summands, output_bytes)) = census_input {
        let (rewritten_summands, rewritten_bytes) = rewrites.rewritten_summands.get();
        let numerator_summands = a
            .0
            .iter()
            .map(|t| vakint::census::summands(t.numerator.as_view()))
            .sum::<usize>();
        let numerator_bytes = a
            .0
            .iter()
            .map(|t| t.numerator.as_view().get_byte_size())
            .sum::<usize>();
        let max_numerator_summands = a
            .0
            .iter()
            .map(|t| vakint::census::summands(t.numerator.as_view()))
            .max()
            .unwrap_or(0);
        vakint::census::record(
            "to_vakint",
            serde_json::json!({
                "materialized_bytes": census_materialized_bytes,
                "rewritten_summands": rewritten_summands,
                "rewritten_bytes": rewritten_bytes,
                "distinct_denominator_products": rewrites.denominator_rewrites.borrow().len(),
                "distinct_indexed_products": rewrites.contraction_images.borrow().len(),
                "distinct_notation_products": notation_images.len(),
                "output_summands": output_summands,
                "output_bytes": output_bytes,
                "terms": a.0.len(),
                "numerator_summands": numerator_summands,
                "max_numerator_summands": max_numerator_summands,
                "numerator_bytes": numerator_bytes,
                "seconds": started.elapsed().as_secs_f64(),
                "materialization_seconds": materialization_seconds,
                "expansion_seconds": expansion_seconds,
                "metric_seconds": metric_seconds,
                "topology_seconds": topology_seconds,
                "split_seconds": split_seconds,
                "routing_seconds": routing_seconds,
                "notation_seconds": notation_seconds,
            }),
        );
    }
    tracing::info!(
        uv = true,
        integrated = true,
        rstar = true,
        vakint = true,
        timing = true,
        stage = "to_vakint_integrand_detail",
        terms = a.0.len(),
        materialization_seconds,
        expansion_seconds,
        metric_seconds,
        topology_seconds,
        split_seconds,
        graph_seconds,
        routing_seconds,
        notation_seconds,
        seconds = started.elapsed().as_secs_f64(),
        "Converted integrand to Vakint"
    );

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
    use crate::{
        finalized_runtime_dot, graph::parse::IntoFinalizedRuntimeGraph,
        initialisation::test_initialise,
    };
    use idenso::representations::Bispinor;
    use linnet::half_edge::involution::EdgeIndex;
    use spenso::structure::slot::{DummyAind, ParseableAind};

    use super::*;

    #[test]
    fn vakint_bridge_contracts_disconnected_children_separately() -> Result<()> {
        test_initialise()?;
        let scalar: crate::graph::Graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=1];
                node [num=1];
                A -> A [particle="scalar_1" num=1 id=0 lmb_id=0
                    source="{ufo_order:0}" sink="{ufo_order:1}"];
            }, "scalars"
        )?;
        let edge_data = scalar.underlying[EdgeIndex(0)].clone();
        // A three-loop parent with two disjoint bubble children. Contracting
        // both bubbles leaves one cycle; identifying their distinct vertices
        // with each other incorrectly turns the connector into another loop.
        let endpoints = [
            (0, 5),
            (0, 6),
            (1, 3),
            (1, 5),
            (1, 5),
            (2, 3),
            (2, 3),
            (2, 4),
            (4, 6),
        ];
        for reverse_order in [false, true] {
            for reverse_edges in [false, true] {
                let mut builder = HedgeGraphBuilder::new();
                let mut nodes = (0..7).map(|_| builder.add_node(())).collect::<Vec<_>>();
                let mut ordered_edges = endpoints.into_iter().enumerate().collect::<Vec<_>>();
                if reverse_order {
                    nodes.reverse();
                    ordered_edges.rotate_left(3);
                }
                let mut edge_ids = [0; 9];
                for (edge, (original, (mut source, mut sink))) in
                    ordered_edges.into_iter().enumerate()
                {
                    edge_ids[original] = edge;
                    if reverse_edges {
                        std::mem::swap(&mut source, &mut sink);
                    }
                    builder.add_edge(nodes[source], nodes[sink], edge_data.clone(), false);
                }
                let graph: HedgeGraph<crate::graph::Edge, ()> = builder.build();
                let mut children: SuBitGraph = graph.empty_subgraph();
                for edge in [3, 4, 5, 6] {
                    children.add(graph[&EdgeIndex(edge_ids[edge])].1);
                }
                assert_eq!(graph.connected_components(&children).len(), 2);
                let momentum = function!(GS.emr_mom, 3);
                for substitute_masses in [true, false] {
                    // Distinct masses prevent equal-denominator fusion from hiding
                    // the spurious loop created by a global child contraction.
                    let mut input = Atom::one();
                    let mut expected = Atom::one();
                    let loop_momentum = function!(vakint::symbols::S.k, 0);
                    for (edge, sign, power) in
                        [(0, 1, 1), (1, -1, 2), (2, 1, 3), (7, 1, 1), (8, 1, 2)]
                    {
                        let sign = if reverse_edges { -sign } else { sign };
                        let mass = Atom::num(edge + 1);
                        input /=
                            function!(GS.den, edge_ids[edge], Atom::num(sign) * &momentum, &mass)
                                .pow(power);
                        let mass_squared = if substitute_masses {
                            Atom::var(GS.m_uv_vacuum).pow(2)
                        } else {
                            mass.pow(2)
                        };
                        expected /=
                            (function!(vakint::symbols::S.dot, &loop_momentum, &loop_momentum)
                                - mass_squared)
                                .pow(power);
                    }
                    let converted = to_vakint_integrand(
                        &input,
                        &graph,
                        &graph.full_filter(),
                        &children,
                        &VakintSettings {
                            additional_normalization: "1".into(),
                            ..Default::default()
                        },
                        substitute_masses,
                    )
                    .wrap_err_with(|| {
                        format!(
                            "reversed order: {reverse_order}, reversed edges: {reverse_edges}, \
                         mass substitution: {substitute_masses}"
                        )
                    })?;
                    let actual = Atom::from(converted)
                        .replace(function!(vakint::symbols::S.topo, W_.x_))
                        .with(W_.x_)
                        .replace(function!(
                            vakint::symbols::S.prop,
                            W_.a_,
                            W_.b_,
                            W_.mom_,
                            W_.mass_,
                            W_.e_
                        ))
                        .with(
                            (function!(vakint::symbols::S.dot, W_.mom_, W_.mom_)
                                - Atom::var(W_.mass_))
                            .pow(-Atom::var(W_.e_)),
                        );
                    assert!(
                        (actual - expected).expand().is_zero(),
                        "the quotient must retain one loop and all five propagator powers"
                    );
                }
            }
        }
        Ok(())
    }

    #[test]
    fn vakint_bridge_routes_tensor_numerators_once_across_propagator_powers() -> Result<()> {
        test_initialise()?;
        let graph: crate::graph::Graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=1];
                node [num=1];
                A -> A [particle="scalar_1" num=1 id=0 lmb_id=0
                    source="{ufo_order:0}" sink="{ufo_order:1}"];
            }, "scalars"
        )?;
        let mink = Minkowski {}.new_rep(GS.dim).to_symbolic([]);
        let mixed = function!(
            ETS.metric,
            function!(GS.emr_mom, 0, &mink),
            function!(GS.emr_mom, 7, &mink)
        );
        let denominator = function!(GS.den, 0, function!(GS.emr_mom, 0), 2);
        let input = &mixed / &denominator + mixed.pow(2) / denominator.pow(3);
        let converted = to_vakint_integrand(
            &input,
            &graph.underlying,
            &graph.full_filter(),
            &graph.empty_subgraph::<SuBitGraph>(),
            &VakintSettings {
                additional_normalization: "z".into(),
                ..Default::default()
            },
            false,
        )?;
        let actual = Atom::from(converted)
            .replace(function!(vakint::symbols::S.topo, W_.x_))
            .with(W_.x_)
            .replace(function!(
                vakint::symbols::S.prop,
                W_.a_,
                W_.b_,
                W_.mom_,
                W_.mass_,
                W_.e_
            ))
            .with(
                (function!(vakint::symbols::S.dot, W_.mom_, W_.mom_) - Atom::var(W_.mass_))
                    .pow(-Atom::var(W_.e_)),
            );
        let loop_momentum = function!(vakint::symbols::S.k, 0);
        let mixed = function!(
            vakint::symbols::S.dot,
            &loop_momentum,
            function!(vakint::symbols::S.p, 7)
        );
        let denominator =
            function!(vakint::symbols::S.dot, &loop_momentum, &loop_momentum) - Atom::num(4);
        let expected = parse!("z") * (&mixed / &denominator + mixed.pow(2) / denominator.pow(3));
        assert!((actual - expected).expand().is_zero());
        Ok(())
    }

    #[test]
    fn vakint_bridge_materializes_zero_color_chains_before_parsing_sums() -> Result<()> {
        test_initialise()?;
        let graph: crate::graph::Graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=1];
                node [num=1];
                A -> A [particle="scalar_1" num=1 id=0 lmb_id=0
                    source="{ufo_order:0}" sink="{ufo_order:1}"];
            }, "scalars"
        )?;
        // The last generator is f(a,b,b), so both apparent tensor summands
        // vanish. Legacy parsing used to compare their exposed interfaces
        // before opening the chain and discovering that zero.
        let chain = parse!(
            "spenso::chain(spenso::coad(8,hedge(0,1)),spenso::coad(8,hedge(11,1)),spenso::f(spenso::in,spenso::out,spenso::coad(8,hedge(1,1))),spenso::f(spenso::in,spenso::out,spenso::coad(8,hedge(11,1))))"
        );
        let mink = Minkowski {}.new_rep(GS.dim);
        let metric = function!(
            ETS.metric,
            mink.to_symbolic([Aind::Hedge(0, 1).to_atom()]),
            mink.to_symbolic([Aind::Hedge(1, 1).to_atom()])
        );
        let momentum = function!(GS.emr_mom, 0);
        let vector = function!(GS.emr_mom, 0, mink.to_symbolic([]));
        let denominator = function!(
            GS.den,
            0,
            momentum,
            Atom::var(GS.m_uv_expansion).pow(2),
            mink.inner_product(&vector, &vector) - Atom::var(GS.m_uv_vacuum).pow(2)
        );
        let integrand = (Atom::i()
            * chain
            * metric
            * (denominator.pow(-1) - Atom::var(GS.m_uv_expansion).pow(2) / denominator.pow(2)))
        .expand();
        assert!(!integrand.is_zero());
        let result = to_vakint_integrand(
            &integrand,
            &graph.underlying,
            &graph.full_filter(),
            &graph.empty_subgraph::<SuBitGraph>(),
            &VakintSettings::default(),
            false,
        )?;
        assert!(result.0.iter().all(|term| term.numerator.is_zero()));
        Ok(())
    }

    #[test]
    fn integrated_two_traces_contract_shared_lorentz_indices_before_poles() -> Result<()> {
        test_initialise()?;
        let epsilon = Atom::var(GS.dim_epsilon);
        let mink = Minkowski {}.new_rep(GS.dim);
        let spin = Bispinor {}.new_rep(4);
        let ports: [_; 4] =
            std::array::from_fn(|i| spin.to_symbolic([Aind::Hedge((i + 3) as u16, 1).to_atom()]));
        let mu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let nu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let gamma =
            |a: &Atom, b: &Atom, index: &Atom| function!(idenso::dirac::AGS.gamma, a, b, index);
        let source = gamma(&ports[0], &ports[1], &mu)
            * gamma(&ports[1], &ports[0], &nu)
            * gamma(&ports[2], &ports[3], &mu)
            * gamma(&ports[3], &ports[2], &nu)
            / epsilon.pow(2);
        let reduced = Integrated::reduce_integrated_numerator(&source)?;
        assert!(!reduced.has_repeated_lorentz_indices());
        let poles = truncate(&series(&reduced, 1)?, false);
        // Tr(gamma_mu gamma_nu) Tr(gamma^mu gamma^nu) = 16 D.
        let expected = Atom::num(64) / epsilon.pow(2) - Atom::num(32) / epsilon;
        assert!((poles - expected).expand().is_zero());
        Ok(())
    }

    #[test]
    fn integrated_unresolved_lorentz_contractions_emit_one_warning() -> Result<()> {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };
        use tracing_subscriber::{Layer, layer::SubscriberExt};

        struct Warnings(Arc<AtomicUsize>);
        impl<S: tracing::Subscriber> Layer<S> for Warnings {
            fn on_event(
                &self,
                event: &tracing::Event<'_>,
                _context: tracing_subscriber::layer::Context<'_, S>,
            ) {
                if *event.metadata().level() == tracing::Level::WARN
                    && event.metadata().fields().field("operation").is_some()
                    && event.metadata().fields().field("stage").is_some()
                {
                    self.0.fetch_add(1, Ordering::Relaxed);
                }
            }
        }

        test_initialise()?;
        let mink = Minkowski {}.new_rep(GS.dim);
        let spin = Bispinor {}.new_rep(4);
        let ports: [_; 4] =
            std::array::from_fn(|i| spin.to_symbolic([Aind::Hedge((i + 3) as u16, 1).to_atom()]));
        let mu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let source = function!(idenso::dirac::AGS.gamma, &ports[0], &ports[1], &mu)
            * function!(idenso::dirac::AGS.gamma, &ports[2], &ports[3], &mu)
            / Atom::var(GS.dim_epsilon);
        let count = Arc::new(AtomicUsize::new(0));
        let subscriber = tracing_subscriber::registry().with(Warnings(Arc::clone(&count)));
        let reduced = tracing::subscriber::with_default(subscriber, || {
            Integrated::reduce_integrated_numerator(&source)
        })?;
        assert!(reduced.has_repeated_lorentz_indices());
        assert_eq!(count.load(Ordering::Relaxed), 1);
        Ok(())
    }

    #[test]
    fn integrated_clifford_contractions_precede_laurent_projection() -> Result<()> {
        test_initialise()?;
        let epsilon = Atom::var(GS.dim_epsilon);
        let mink = Minkowski {}.new_rep(GS.dim);
        let spin = Bispinor {}.new_rep(4);
        let a = spin.to_symbolic([Aind::Hedge(3, 1).to_atom()]);
        let b = spin.to_symbolic([Aind::Hedge(4, 1).to_atom()]);
        let mu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let momentum = function!(GS.emr_mom, 7, mink.to_symbolic([]));
        let gamma = |argument: &Atom| {
            function!(
                idenso::dirac::AGS.gamma,
                SPENSO_TAG.chain_in,
                SPENSO_TAG.chain_out,
                argument
            )
        };
        let slash = spenso::chain!(&a, &b, gamma(&momentum));
        for (word, coefficient) in [
            (
                spenso::chain!(&a, &b, gamma(&mu), gamma(&mu), gamma(&momentum)),
                Atom::var(GS.dim),
            ),
            (
                spenso::chain!(&a, &b, gamma(&mu), gamma(&momentum), gamma(&mu)),
                Atom::num(2) - Atom::var(GS.dim),
            ),
        ] {
            for pole_order in [1, 2] {
                let reduced = Integrated::reduce_integrated_numerator(
                    &(word.clone() / epsilon.pow(pole_order)),
                )?;
                let expected = Integrated::dimensionally_regularized(
                    &(&coefficient * &slash / epsilon.pow(pole_order)),
                );
                let expansion = series(&reduced, 1)?;
                let expected_expansion = series(&expected, 1)?;
                for finite in [false, true] {
                    assert!(
                        (truncate(&expansion, finite) - truncate(&expected_expansion, finite))
                            .expand()
                            .is_zero(),
                        "finite={finite}, pole_order={pole_order}, reduced={reduced}, expected={expected}"
                    );
                }
            }
        }
        // Projecting the unreduced Clifford word first, as the old integration
        // path did, leaves a finite -2 slash(p) in a purported MS pole. Once
        // inserted in a second logarithmically divergent loop, it changes the
        // overall simple pole by -2 slash(p)/epsilon.
        let unreduced = spenso::chain!(&a, &b, gamma(&mu), gamma(&mu), gamma(&momentum)) / &epsilon;
        let projected_first = truncate(&series(&unreduced, 1)?, false);
        let reduced_late = Integrated::reduce_integrated_numerator(&projected_first)?;
        let correct_pole = Atom::num(4) * &slash / &epsilon;
        let extra_parent_pole = truncate(
            &series(&((reduced_late - correct_pole) / &epsilon), 1)?,
            false,
        );
        assert!(
            (extra_parent_pole + Atom::num(2) * slash / epsilon)
                .expand()
                .is_zero()
        );
        Ok(())
    }

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
        let regulated_dimension = Atom::num(4) - Atom::num(2) * &epsilon;
        let unexpanded_dimension = Atom::num(4) - (Atom::num(4) - &regulated_dimension);
        let unexpanded = source
            .replace(regulated_dimension.to_pattern())
            .with(unexpanded_dimension);
        assert_eq!(
            Integrated::restore_lorentz_dimensions(&unexpanded),
            normalized,
            "backend dimension substitutions need algebraic comparison"
        );
        let expanded = series(&normalized, 1).unwrap().to_atom();
        assert!(!expanded.contains_symbol(symbolica::atom::Symbol::DERIVATIVE));
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let restored = Integrated::restore_numerator(expanded, 0);
        let contracted =
            SymbolicTensor::infer_indexed(cooking.try_cook(restored.as_view()).unwrap())
                .unwrap()
                .contract(ContractSettings {
                    rank_one: false,
                    collect_chains: false,
                    collect_traces: false,
                    ..Default::default()
                })
                .unwrap();
        let restored = cooking.uncook(contracted.expression().as_view());
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
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let restored = Integrated::restore_numerator(source, 0);
        let contracted =
            SymbolicTensor::infer_indexed(cooking.try_cook(restored.as_view()).unwrap())
                .unwrap()
                .contract(ContractSettings {
                    rank_one: false,
                    collect_chains: false,
                    collect_traces: false,
                    ..Default::default()
                })
                .unwrap();
        let restored = cooking.uncook(contracted.expression().as_view());
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
            let wrapped_dummy = mink.to_symbolic([vakint::symbols::S.dot_dummy_ind(5)]);
            assert_eq!(
                Integrated::restore_numerator(function!(source, 2, wrapped_dummy), 19),
                function!(target, 2, mink.to_symbolic([function!(GS.uvaind, 19, 5)]))
            );
        }

        let spinor = Bispinor {}.new_rep(4);
        let a = spinor.to_symbolic([Atom::from(Aind::Hedge(1, 0))]);
        let b = spinor.to_symbolic([Atom::from(Aind::Hedge(2, 0))]);
        let backend_slot = mink.to_symbolic([vakint::symbols::S.dot_dummy_ind(3)]);
        let restored_slot = mink.to_symbolic([function!(GS.uvaind, 19, 3)]);
        assert_eq!(
            Integrated::restore_numerator(
                function!(idenso::dirac::AGS.gamma, &a, &b, backend_slot),
                19
            ),
            function!(idenso::dirac::AGS.gamma, a, b, restored_slot)
        );

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
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let restored = Integrated::restore_numerator(source, 19);
        let tensor = SymbolicTensor::infer_indexed(cooking.try_cook(restored.as_view()).unwrap())
            .unwrap()
            .contract(ContractSettings::default())
            .unwrap()
            .to_dots()
            .unwrap();
        let restored = cooking.uncook(tensor.expression().as_view());
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

        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let tensor = SymbolicTensor::infer_indexed(cooking.try_cook(numerator.as_view()).unwrap())
            .unwrap()
            .contract(ContractSettings::default())
            .unwrap()
            .to_dots()
            .unwrap();
        let converted = cooking
            .uncook(tensor.expression().as_view())
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
